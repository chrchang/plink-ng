// This file is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink2_set.h"

#ifdef __cplusplus
namespace plink2 {
#endif

typedef struct MakeSetRangeStruct {
  NONCOPYABLE(MakeSetRangeStruct);
  struct MakeSetRangeStruct* next;
  uint32_t uidx_start;
  uint32_t uidx_end;
} MakeSetRange;

static_assert(kMaxChrCodeDigits == 5, "LoadIntervalBed() must be updated.");
PglErr LoadIntervalBed(const ChrInfo* cip, const uint32_t* variant_bps, const char* sorted_subset_ids, const char* file_descrip, uint32_t zero_based, uint32_t track_set_names, uint32_t border_extend, uint32_t fail_on_no_sets, uint32_t c_prefix, uintptr_t subset_ct, uintptr_t max_subset_id_blen, TextStream* txsp, uintptr_t* set_ct_ptr, char** set_names_ptr, uintptr_t* max_set_id_blen_ptr, uint64_t** range_sort_buf_ptr, MakeSetRange*** make_set_range_arr_ptr) {
  // In plink 1.9, this was named load_range_list() and called directly by
  // ExtractExcludeRange(), define_sets(), and indirectly by annotate(),
  // gene_report(), and clump_reports().  That function required set IDs in
  // column 4, and interpreted column 5 as a set "group label".
  // However, column 5 wasn't used very often, and without it, we're free to
  // generalize this function to point at UCSC interval-BED files.
  //
  // Assumes caller will reset g_bigstack_end later.
  PglErr reterr = kPglRetSuccess;
  {
    LlStr* make_set_ll = nullptr;
    char* set_names = nullptr;
    uintptr_t set_ct = 0;
    uintptr_t max_set_id_blen = 0;
    // if we need to track set names, put together a sorted list
    if (track_set_names) {
      uintptr_t line_idx = 1;
      for (char* line_iter = TextLineEnd(txsp); TextGetUnsafe(txsp, &line_iter); ++line_idx) {
        char* line_start = line_iter;
        char* first_token_end = CurTokenEnd(line_start);
        char* cur_set_id = NextTokenMult(first_token_end, 3);
        char* last_token = cur_set_id;
        if (unlikely(NoMoreTokensKns(last_token))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, file_descrip);
          goto LoadIntervalBed_ret_MALFORMED_INPUT_WW;
        }
        const uint32_t chr_name_slen = first_token_end - line_start;
        *first_token_end = '\0';
        const uint32_t cur_chr_code = GetChrCode(line_start, cip, chr_name_slen);
        if (IsI32Neg(cur_chr_code)) {
          // kludge (21 Jan 2020): --extract/exclude range should not error out
          // if a line mentions a chromosome code not in the dataset.
          // TODO: better condition for skipping this line.  (not totally
          // trivial since some future callers may need to track empty sets)
          if (likely((!set_ct_ptr) && (cur_chr_code == UINT32_MAX))) {
            continue;
          }
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid chromosome code on line %" PRIuPTR " of %s.\n", line_idx, file_descrip);
          goto LoadIntervalBed_ret_MALFORMED_INPUT_WW;
        }
        // chr_mask check removed, we want to track empty sets
        uint32_t set_id_slen = strlen_se(cur_set_id);
        // we're about to possibly clobber \n, so advance line_iter now
        line_iter = AdvPastDelim(&(cur_set_id[set_id_slen]), '\n');
        cur_set_id[set_id_slen] = '\0';
        if (subset_ct) {
          if (bsearch_str(cur_set_id, sorted_subset_ids, set_id_slen, max_subset_id_blen, subset_ct) == -1) {
            continue;
          }
        }
        // when there are repeats, they are likely to be next to each other
        if (make_set_ll && strequal_overread(make_set_ll->str, last_token)) {
          continue;
        }
        uint32_t set_id_blen = set_id_slen + 1;
        // argh, --clump counts positional overlaps which don't include any
        // variants in the dataset.  So we prefix set IDs with a chromosome
        // index in that case (with leading zeroes) and treat cross-chromosome
        // sets as distinct.
        if (!variant_bps) {
          set_id_blen += kMaxChrCodeDigits;
        }
        if (set_id_blen > max_set_id_blen) {
          max_set_id_blen = set_id_blen;
        }
        LlStr* ll_tmp;
        if (unlikely(bigstack_end_alloc_llstr(set_id_blen, &ll_tmp))) {
          goto LoadIntervalBed_ret_NOMEM;
        }
        ll_tmp->next = make_set_ll;
        if (variant_bps) {
          memcpy(ll_tmp->str, last_token, set_id_blen);
        } else {
          u32toa_z5(cur_chr_code, ll_tmp->str);
          // if first character of gene name is a digit, natural sort has
          // strange effects unless we force [3] to be nonnumeric...
          ll_tmp->str[kMaxChrCodeDigits - 1] -= 15;
          memcpy(&(ll_tmp->str[kMaxChrCodeDigits]), last_token, set_id_blen - kMaxChrCodeDigits);
        }
        make_set_ll = ll_tmp;
        ++set_ct;
      }
      if (unlikely(TextStreamErrcode2(txsp, &reterr))) {
        goto LoadIntervalBed_ret_TSTREAM_FAIL;
      }
      if (!set_ct) {
        if (unlikely(fail_on_no_sets)) {
          if (variant_bps) {
            logerrputs("Error: All variants excluded by --gene[-all], since no sets were defined from\n--make-set file.\n");
            reterr = kPglRetMalformedInput;
            goto LoadIntervalBed_ret_1;
          } else {
            if (subset_ct) {
              logerrputs("Error: No --gene-subset genes present in --gene-report file.\n");
              reterr = kPglRetInconsistentInput;
            } else {
              logerrputs("Error: Empty --gene-report file.\n");
              reterr = kPglRetMalformedInput;
            }
            goto LoadIntervalBed_ret_1;
          }
        }
        logerrprintfww("Warning: No valid ranges in %s.\n", file_descrip);
        goto LoadIntervalBed_ret_1;
      }
      // c_prefix is 0 or 2
      max_set_id_blen += c_prefix;
      if (unlikely(max_set_id_blen > kMaxIdBlen)) {
        logerrputs("Error: Set IDs are limited to " MAX_ID_SLEN_STR " characters.\n");
        goto LoadIntervalBed_ret_MALFORMED_INPUT;
      }
      const char** strptr_arr;
      if (unlikely(
              bigstack_alloc_c(set_ct * max_set_id_blen, set_names_ptr) ||
              bigstack_alloc_kcp(set_ct, &strptr_arr))) {
        goto LoadIntervalBed_ret_NOMEM;
      }
      set_names = *set_names_ptr;
      for (uintptr_t set_idx = 0; set_idx != set_ct; ++set_idx) {
        strptr_arr[set_idx] = make_set_ll->str;
        make_set_ll = make_set_ll->next;
      }
      StrptrArrNsort(set_ct, strptr_arr);
      set_ct = CopyAndDedupSortedStrptrsToStrbox(strptr_arr, set_ct, max_set_id_blen, &(set_names[c_prefix]));
      if (c_prefix) {
        for (uintptr_t set_idx = 0; set_idx != set_ct; ++set_idx) {
          memcpy_k(&(set_names[set_idx * max_set_id_blen]), "C_", 2);
        }
      }
      BigstackShrinkTop(set_names, set_ct * max_set_id_blen);
      reterr = TextRewind(txsp);
      if (unlikely(reterr)) {
        goto LoadIntervalBed_ret_TSTREAM_FAIL;
      }
    } else {
      set_ct = 1;
    }
    MakeSetRange** make_set_range_arr = S_CAST(MakeSetRange**, bigstack_end_alloc(set_ct * sizeof(intptr_t)));
    if (unlikely(!make_set_range_arr)) {
      goto LoadIntervalBed_ret_NOMEM;
    }
    ZeroPtrArr(set_ct, make_set_range_arr);
    uintptr_t line_idx = 0;
    uint32_t chr_start = 0;
    uint32_t chr_end = 0;
    for (char* line_iter = &(TextLineEnd(txsp)[-1]); ; line_iter = AdvToDelim(line_iter, '\n')) {
    LoadIntervalBed_LINE_ITER_ALREADY_ADVANCED:
      ++line_iter;
      ++line_idx;
      reterr = TextGetUnsafe(txsp, &line_iter);
      if (reterr) {
        if (likely(reterr == kPglRetEof)) {
          if (unlikely(track_set_names && (line_idx == 1))) {
            goto LoadIntervalBed_ret_REWIND_FAIL;
          }
          reterr = kPglRetSuccess;
          break;
        }
        goto LoadIntervalBed_ret_TSTREAM_FAIL;
      }
      char* line_start = line_iter;
      char* first_token_end = CurTokenEnd(line_start);
      char* last_token = NextTokenMult(first_token_end, 2 + track_set_names);
      if (unlikely(NoMoreTokensKns(last_token))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, file_descrip);
        goto LoadIntervalBed_ret_MALFORMED_INPUT_WW;
      }
      const uint32_t chr_name_slen = first_token_end - line_start;
      *first_token_end = '\0';
      const uint32_t cur_chr_code = GetChrCode(line_start, cip, chr_name_slen);
      if (IsI32Neg(cur_chr_code)) {
        if (likely((!set_ct_ptr) && (cur_chr_code == UINT32_MAX))) {
          continue;
        }
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid chromosome code on line %" PRIuPTR " of %s.\n", line_idx, file_descrip);
        goto LoadIntervalBed_ret_MALFORMED_INPUT_WW;
      }
      line_iter = CurTokenEnd(last_token);
      if (!IsSet(cip->chr_mask, cur_chr_code)) {
        continue;
      }
      if (variant_bps) {
        const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[cur_chr_code];
        chr_start = cip->chr_fo_vidx_start[chr_fo_idx];
        chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        if (chr_end == chr_start) {
          continue;
        }
        // (track_set_names must be true if subset_ct is nonzero)
        // might need to move this outside the if-statement later
        if (subset_ct && (bsearch_str(last_token, sorted_subset_ids, strlen_se(last_token), max_subset_id_blen, subset_ct) == -1)) {
          continue;
        }
      }
      const char* linebuf_iter = FirstNonTspace(&(first_token_end[1]));
      uint32_t range_first;
      if (unlikely(ScanmovUintDefcap(&linebuf_iter, &range_first))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid range start position on line %" PRIuPTR " of %s.\n", line_idx, file_descrip);
        goto LoadIntervalBed_ret_MALFORMED_INPUT_WW;
      }
      range_first += zero_based;
      linebuf_iter = NextToken(linebuf_iter);
      uint32_t range_last;
      if (unlikely(ScanmovUintDefcap(&linebuf_iter, &range_last))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid range end position on line %" PRIuPTR " of %s.\n", line_idx, file_descrip);
        goto LoadIntervalBed_ret_MALFORMED_INPUT_WW;
      }
      if (unlikely(range_last < range_first)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Range end position smaller than range start on line %" PRIuPTR " of %s.\n", line_idx, file_descrip);
        goto LoadIntervalBed_ret_MALFORMED_INPUT_WW;
      }
      if (border_extend > range_first) {
        range_first = 0;
      } else {
        range_first -= border_extend;
      }
      range_last += border_extend;
      const uint32_t last_token_slen = line_iter - last_token;
      // about to potentially clobber \n, advance now
      line_iter = AdvToDelim(line_iter, '\n');
      uint32_t cur_set_idx = 0;
      if (set_ct > 1) {
        // bugfix: bsearch_str_natural requires null-terminated string
        last_token[last_token_slen] = '\0';
        if (c_prefix) {
          last_token = &(last_token[-2]);
          memcpy_k(last_token, "C_", 2);
        } else if (!variant_bps) {
          last_token = &(last_token[-S_CAST(int32_t, kMaxChrCodeDigits)]);
          u32toa_z5(cur_chr_code, last_token);
          last_token[kMaxChrCodeDigits - 1] -= 15;
        }
        // this should never fail
        cur_set_idx = bsearch_str_natural(last_token, set_names, max_set_id_blen, set_ct);
      }
      if (variant_bps) {
        // translate to within-chromosome uidx
        range_first = CountSortedSmallerU32(&(variant_bps[chr_start]), chr_end - chr_start, range_first);
        range_last = CountSortedSmallerU32(&(variant_bps[chr_start]), chr_end - chr_start, range_last + 1);
        if (range_last > range_first) {
          MakeSetRange* msr_tmp = S_CAST(MakeSetRange*, bigstack_end_alloc(sizeof(MakeSetRange)));
          if (unlikely(!msr_tmp)) {
            goto LoadIntervalBed_ret_NOMEM;
          }
          msr_tmp->next = make_set_range_arr[cur_set_idx];
          // normally, I'd keep chr_idx here since that enables by-chromosome
          // sorting, but that's probably not worth bloating MakeSetRange
          // from 16 to 32 bytes
          msr_tmp->uidx_start = chr_start + range_first;
          msr_tmp->uidx_end = chr_start + range_last;
          make_set_range_arr[cur_set_idx] = msr_tmp;
        }
      } else {
        MakeSetRange* msr_tmp = S_CAST(MakeSetRange*, bigstack_end_alloc(sizeof(MakeSetRange)));
        if (unlikely(!msr_tmp)) {
          goto LoadIntervalBed_ret_NOMEM;
        }
        msr_tmp->next = make_set_range_arr[cur_set_idx];
        msr_tmp->uidx_start = range_first;
        msr_tmp->uidx_end = range_last + 1;
        make_set_range_arr[cur_set_idx] = msr_tmp;
      }
      goto LoadIntervalBed_LINE_ITER_ALREADY_ADVANCED;
    }
    // allocate buffer for sorting ranges later
    uint32_t max_set_range_ct = 0;
    for (uint32_t set_idx = 0; set_idx != set_ct; ++set_idx) {
      uint32_t cur_set_range_ct = 0;
      MakeSetRange* msr_tmp = make_set_range_arr[set_idx];
      while (msr_tmp) {
        ++cur_set_range_ct;
        msr_tmp = msr_tmp->next;
      }
      if (cur_set_range_ct > max_set_range_ct) {
        max_set_range_ct = cur_set_range_ct;
      }
    }
    if (range_sort_buf_ptr) {
      if (unlikely(bigstack_end_alloc_u64(max_set_range_ct, range_sort_buf_ptr))) {
        goto LoadIntervalBed_ret_NOMEM;
      }
    }
    if (set_ct_ptr) {
      *set_ct_ptr = set_ct;
    }
    if (max_set_id_blen_ptr) {
      *max_set_id_blen_ptr = max_set_id_blen;
    }
    *make_set_range_arr_ptr = make_set_range_arr;
  }
  while (0) {
  LoadIntervalBed_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadIntervalBed_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, file_descrip);
    reterr = kPglRetRewindFail;
    break;
  LoadIntervalBed_ret_TSTREAM_FAIL:
    TextStreamErrPrint(file_descrip, txsp);
    break;
  LoadIntervalBed_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  LoadIntervalBed_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 LoadIntervalBed_ret_1:
  return reterr;
}

PglErr ExtractExcludeRange(const char* fnames, const ChrInfo* cip, const uint32_t* variant_bps, uint32_t raw_variant_ct, VfilterType vft, uint32_t zero_based, uint32_t bed_border_bp, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t orig_variant_ct = *variant_ct_ptr;
  if (!orig_variant_ct) {
    return kPglRetSuccess;
  }
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  const char* fname_txs = nullptr;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    const uintptr_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* variant_include_mask = nullptr;
    if (vft != kVfilterExclude) {
      if (unlikely(bigstack_calloc_w(raw_variant_ctl, &variant_include_mask))) {
        goto ExtractExcludeRange_ret_NOMEM;
      }
    }
    const char* fnames_iter = fnames;
    do {
      if (fnames_iter == fnames) {
        fname_txs = fnames_iter;
        reterr = InitTextStream(fnames_iter, kTextStreamBlenFast, MAXV(max_thread_ct - 1, 1), &txs);
        if (unlikely(reterr)) {
          goto ExtractExcludeRange_ret_TSTREAM_FAIL;
        }
      } else {
        reterr = TextRetarget(fnames_iter, &txs);
        if (unlikely(reterr)) {
          goto ExtractExcludeRange_ret_TSTREAM_FAIL;
        }
        fname_txs = fnames_iter;
      }
      MakeSetRange** range_arr = nullptr;
      reterr = LoadIntervalBed(cip, variant_bps, nullptr, fname_txs, zero_based, 0, bed_border_bp, 0, 0, 0, 0, &txs, nullptr, nullptr, nullptr, nullptr, &range_arr);
      if (unlikely(reterr)) {
        goto ExtractExcludeRange_ret_1;
      }
      MakeSetRange* msr_tmp = range_arr[0];
      if (vft == kVfilterExclude) {
        while (msr_tmp) {
          ClearBitsNz(msr_tmp->uidx_start, msr_tmp->uidx_end, variant_include);
          msr_tmp = msr_tmp->next;
        }
      } else {
        while (msr_tmp) {
          FillBitsNz(msr_tmp->uidx_start, msr_tmp->uidx_end, variant_include_mask);
          msr_tmp = msr_tmp->next;
        }
        if (vft == kVfilterExtractIntersect) {
          BitvecAnd(variant_include_mask, raw_variant_ctl, variant_include);
          ZeroWArr(raw_variant_ctl, variant_include_mask);
        }
      }
      fnames_iter = strnul(fnames_iter);
      ++fnames_iter;
    } while (*fnames_iter);
    if (vft == kVfilterExtract) {
      BitvecAnd(variant_include_mask, raw_variant_ctl, variant_include);
    }
    *variant_ct_ptr = PopcountWords(variant_include, raw_variant_ctl);
    const char* vft_name = g_vft_names[vft];
    if (*variant_ct_ptr == orig_variant_ct) {
      logerrprintf("Warning: No variants excluded by '--%s bed%c'.\n", vft_name, '1' - zero_based);
    } else {
      const uint32_t excluded_ct = orig_variant_ct - (*variant_ct_ptr);
      logprintf("--%s bed%c: %u variant%s excluded.\n", vft_name, '1' - zero_based, excluded_ct, (excluded_ct == 1)? "" : "s");
    }
  }
  while (0) {
  ExtractExcludeRange_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExtractExcludeRange_ret_TSTREAM_FAIL:
    TextStreamErrPrint(fname_txs, &txs);
    break;
  }
 ExtractExcludeRange_ret_1:
  if (fname_txs) {
    CleanupTextStream2(fname_txs, &txs, &reterr);
  }
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}
#endif
