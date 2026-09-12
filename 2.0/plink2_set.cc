// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang, Benjamin Demaille.
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

#include "plink2_compress_stream.h"

#include <string.h>

#include "include/plink2_bits.h"
#include "include/plink2_simd.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "plink2_cmdline.h"
#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

typedef struct MakeSetRangeStruct {
  NONCOPYABLE(MakeSetRangeStruct);
  struct MakeSetRangeStruct* next;
  uint32_t uidx_start;
  uint32_t uidx_end;
} MakeSetRange;

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
      for (char* line_iter = TextLineEnd(txsp); TextGetUnsafe2(txsp, &line_iter); ++line_idx) {
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
          if (bsearch_strbox(cur_set_id, sorted_subset_ids, set_id_slen, max_subset_id_blen, subset_ct) == -1) {
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
          u32toa_zchr(cur_chr_code, ll_tmp->str);
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
      if (unlikely(bigstack_alloc_c(set_ct * max_set_id_blen, set_names_ptr) ||
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
      }
      if (subset_ct && (bsearch_strbox(last_token, sorted_subset_ids, strlen_se(last_token), max_subset_id_blen, subset_ct) == -1)) {
        continue;
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
        // bugfix: bsearch_strbox_natural requires null-terminated string
        last_token[last_token_slen] = '\0';
        if (c_prefix) {
          last_token = &(last_token[-2]);
          memcpy_k(last_token, "C_", 2);
        } else if (!variant_bps) {
          last_token = &(last_token[-S_CAST(int32_t, kMaxChrCodeDigits)]);
          u32toa_zchr(cur_chr_code, last_token);
          last_token[kMaxChrCodeDigits - 1] -= 15;
        }
        // this should never fail
        cur_set_idx = bsearch_strbox_natural(last_token, set_names, max_set_id_blen, set_ct);
      }
      if (variant_bps) {
        // translate to within-chromosome uidx
        range_first = LowerBoundNonemptyU32(&(variant_bps[chr_start]), chr_end - chr_start, range_first);
        range_last = LowerBoundNonemptyU32(&(variant_bps[chr_start]), chr_end - chr_start, range_last + 1);
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

uint32_t IntervalInSetdef(const uint32_t* setdef, uint32_t variant_uidx_start, uint32_t variant_uidx_end) {
  // - expects half-open interval coordinates as input
  // - assumes variant_uidx_end > variant_uidx_start
  // - returns 0 if intersection empty, nonzero value (not necessarily 1) if
  //   nonempty
  const uint32_t range_ct = setdef[0];
  if (range_ct != UINT32_MAX) {
    // This overlap query may belong in plink2_common.
    if (!range_ct) {
      return 0;
    }
    const uint32_t* set_ranges_start = &(setdef[1]);
    // Check whether variant_uidx_start is contained within an interval.
    // - If an interval is of the form [x, variant_uidx_start),
    //   variant_uidx_start is not within that interval.  We want start_pos to
    //   be the index past the interval.
    // - If an interval is of the form [x, variant_uidx_start + 1),
    //   variant_uidx_start must be within the interval.
    // std::lower_bound(set_ranges_start, &(set_ranges_start[range_ct * 2]),
    //                  variant_uidx_start + 1)
    // has the correct behavior.
    const uint32_t start_pos = LowerBoundNonemptyU32(set_ranges_start, range_ct * 2, variant_uidx_start + 1);
    if (start_pos & 1) {
      return 1;
    }
    if (start_pos == range_ct * 2) {
      // Last interval ends before variant_uidx_start, intersection must be
      // empty.
      return 0;
    }
    return LowerBoundNonemptyU32(&(set_ranges_start[start_pos]), range_ct * 2 - start_pos, variant_uidx_end);
  }

  const uint32_t set_uidx_base = setdef[1];
  if ((variant_uidx_end <= set_uidx_base) || (variant_uidx_start >= set_uidx_base + setdef[2])) {
    // Interval does not intersect set bitvector.
    return setdef[3];
  }
  uint32_t idx_start;
  if (variant_uidx_start < set_uidx_base) {
    if (setdef[3]) {
      return 1;
    }
    idx_start = 0;
  } else {
    idx_start = variant_uidx_start - set_uidx_base;
  }
  uint32_t idx_end;
  if (variant_uidx_end > set_uidx_base + setdef[2]) {
    if (setdef[3]) {
      return 1;
    }
    idx_end = setdef[2];
  } else {
    idx_end = variant_uidx_end - set_uidx_base;
  }
  const uint32_t first_hit = AdvBoundedTo1Bit(R_CAST(const uintptr_t*, &(setdef[4])), idx_start, idx_end);
  return (first_hit < idx_end);
}

void InitSet(SetInfo* sip) {
  sip->fname = nullptr;
  sip->setnames_flattened = nullptr;
  sip->merged_set_name = nullptr;
  sip->flags = kfSet0;
}

void CleanupSet(SetInfo* sip) {
  free_cond(sip->fname);
  free_cond(sip->setnames_flattened);
  free_cond(sip->merged_set_name);
}

uint32_t InSetdef(const uint32_t* setdef, uint32_t variant_idx) {
  const uint32_t range_ct = setdef[0];
  if (range_ct != UINT32_MAX) {
    for (uint32_t range_idx = 0; range_idx != range_ct; ++range_idx) {
      const uint32_t range_start = setdef[range_idx * 2 + 1];
      if (variant_idx < range_start) {
        return 0;
      }
      if (variant_idx < setdef[range_idx * 2 + 2]) {
        return 1;
      }
    }
    return 0;
  }
  const uint32_t offset = setdef[1];
  const uint32_t bit_ct = setdef[2];
  if ((variant_idx < offset) || (variant_idx >= offset + bit_ct)) {
    return setdef[3];
  }
  return IsSet(R_CAST(const uintptr_t*, &(setdef[4])), variant_idx - offset);
}

uint32_t SetdefSize(const uint32_t* setdef, uint32_t variant_ct) {
  const uint32_t range_ct = setdef[0];
  if (range_ct != UINT32_MAX) {
    uint32_t total = 0;
    for (uint32_t range_idx = 0; range_idx != range_ct; ++range_idx) {
      total += setdef[range_idx * 2 + 2] - setdef[range_idx * 2 + 1];
    }
    return total;
  }
  const uint32_t offset = setdef[1];
  const uint32_t bit_ct = setdef[2];
  uint32_t total = PopcountWords(R_CAST(const uintptr_t*, &(setdef[4])), BitCtToWordCt(bit_ct));
  if (setdef[3]) {
    total += offset + (variant_ct - MINV(variant_ct, offset + bit_ct));
  }
  return total;
}

void UnpackSetdef(const uint32_t* setdef, uint32_t variant_ct, uintptr_t* include_bitvec) {
  const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
  const uint32_t range_ct = setdef[0];
  if (range_ct != UINT32_MAX) {
    ZeroWArr(variant_ctl, include_bitvec);
    for (uint32_t range_idx = 0; range_idx != range_ct; ++range_idx) {
      const uint32_t range_start = setdef[range_idx * 2 + 1];
      const uint32_t range_end = MINV(setdef[range_idx * 2 + 2], variant_ct);
      if (range_end > range_start) {
        FillBitsNz(range_start, range_end, include_bitvec);
      }
    }
    return;
  }
  const uint32_t offset = setdef[1];
  const uint32_t bit_ct = setdef[2];
  if (setdef[3]) {
    SetAllBits(variant_ct, include_bitvec);
    if (offset < variant_ct) {
      ClearBitsNz(offset, MINV(offset + bit_ct, variant_ct), include_bitvec);
    }
  } else {
    ZeroWArr(variant_ctl, include_bitvec);
  }
  const uintptr_t* src = R_CAST(const uintptr_t*, &(setdef[4]));
  const uint32_t copy_end = MINV(offset + bit_ct, variant_ct);
  for (uint32_t variant_idx = offset; variant_idx != copy_end; ++variant_idx) {
    if (IsSet(src, variant_idx - offset)) {
      SetBit(variant_idx, include_bitvec);
    }
  }
  ZeroTrailingBits(variant_ct, include_bitvec);
}

// Emits the smaller of the two setdef forms for one set, whose members are
// given as a bitvector over the filtered variant space.  The range form costs
// eight bytes per run of consecutive members, the bitfield form sixteen bytes
// plus one bit per variant between the first and last member, so which one
// wins depends on how localized the set is.
//
// member_bitvec must have a zeroed word past its variant_ct bits, so that a
// run ending on the last variant still terminates the scan for its end.
static uint32_t* SaveSetdef(const uintptr_t* member_bitvec, uint32_t variant_ct, unsigned char** alloc_basep, unsigned char* alloc_end) {
  const uint32_t first_member = AdvBoundedTo1Bit(member_bitvec, 0, variant_ct);
  uint32_t range_ct = 0;
  uint32_t last_end = 0;
  if (first_member != variant_ct) {
    uint32_t range_start = first_member;
    while (1) {
      const uint32_t range_end = AdvTo0Bit(member_bitvec, range_start);
      ++range_ct;
      last_end = range_end;
      if (range_end >= variant_ct) {
        last_end = variant_ct;
        break;
      }
      const uint32_t next_start = AdvBoundedTo1Bit(member_bitvec, range_end, variant_ct);
      if (next_start == variant_ct) {
        break;
      }
      range_start = next_start;
    }
  }
  const uintptr_t range_form_size = (S_CAST(uintptr_t, range_ct) * 2 + 1) * sizeof(int32_t);
  uintptr_t bitfield_form_size = UINTPTR_MAX;
  uint32_t bitfield_offset = 0;
  uint32_t bitfield_bit_ct = 0;
  if (range_ct) {
    // The offset is a multiple of 128, as the encoding requires.
    bitfield_offset = first_member & (~127);
    bitfield_bit_ct = RoundUpPow2(last_end - bitfield_offset, 128);
    bitfield_form_size = 4 * sizeof(int32_t) + bitfield_bit_ct / CHAR_BIT;
  }
  if ((!range_ct) || (range_form_size <= bitfield_form_size)) {
    const uintptr_t alloc_size = RoundUpPow2(range_form_size, 16);
    if (S_CAST(uintptr_t, alloc_end - (*alloc_basep)) < alloc_size) {
      return nullptr;
    }
    uint32_t* setdef = R_CAST(uint32_t*, *alloc_basep);
    *alloc_basep = &((*alloc_basep)[alloc_size]);
    setdef[0] = range_ct;
    uint32_t* write_iter = &(setdef[1]);
    uint32_t range_start = first_member;
    for (uint32_t range_idx = 0; range_idx != range_ct; ++range_idx) {
      const uint32_t range_end = MINV(AdvTo0Bit(member_bitvec, range_start), variant_ct);
      *write_iter++ = range_start;
      *write_iter++ = range_end;
      if (range_end == variant_ct) {
        break;
      }
      range_start = AdvBoundedTo1Bit(member_bitvec, range_end, variant_ct);
    }
    return setdef;
  }
  const uintptr_t alloc_size = RoundUpPow2(bitfield_form_size, 16);
  if (S_CAST(uintptr_t, alloc_end - (*alloc_basep)) < alloc_size) {
    return nullptr;
  }
  uint32_t* setdef = R_CAST(uint32_t*, *alloc_basep);
  *alloc_basep = &((*alloc_basep)[alloc_size]);
  setdef[0] = UINT32_MAX;
  setdef[1] = bitfield_offset;
  setdef[2] = bitfield_bit_ct;
  setdef[3] = 0;
  uintptr_t* dst = R_CAST(uintptr_t*, &(setdef[4]));
  const uint32_t dst_word_ct = bitfield_bit_ct / kBitsPerWord;
  ZeroWArr(dst_word_ct, dst);
  for (uint32_t variant_idx = first_member; variant_idx != last_end; ++variant_idx) {
    if (IsSet(member_bitvec, variant_idx)) {
      SetBit(variant_idx - bitfield_offset, dst);
    }
  }
  return setdef;
}

// True if a set with this name is wanted.
static uint32_t SetNameWanted(const char* name, uint32_t name_slen, const char* setnames_flattened) {
  if (setnames_flattened) {
    const char* iter = setnames_flattened;
    while (*iter) {
      const uint32_t cur_slen = strlen(iter);
      if ((cur_slen == name_slen) && memequal(iter, name, name_slen)) {
        return 1;
      }
      iter = &(iter[cur_slen + 1]);
    }
    return 0;
  }
  return 1;
}

PglErr DefineSets(const SetInfo* sip, const uintptr_t* variant_include, const char* const* variant_ids, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t max_thread_ct, VariantSets* vsp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  const char* fname = sip->fname;
  PglErr reterr = kPglRetSuccess;
  TokenStream tks;
  PreinitTokenStream(&tks);
  {
    const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
    const char* merged_set_name = sip->merged_set_name;

    uint32_t* variant_id_htable;
    uint32_t* htable_dup_base;
    uint32_t variant_id_htable_size;
    uint32_t dup_ct;
    reterr = AllocAndPopulateIdHtableMt(variant_include, variant_ids, variant_ct, bigstack_left() / 8, max_thread_ct, &variant_id_htable, &htable_dup_base, &variant_id_htable_size, &dup_ct);
    if (unlikely(reterr)) {
      goto DefineSets_ret_1;
    }
    unsigned char* htable_mark = g_bigstack_base;

    // Variant uidxs have to be turned into filtered indexes, which is what a
    // setdef is over.
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uint32_t* variant_include_cumulative_popcounts;
    if (unlikely(bigstack_end_alloc_u32(raw_variant_ctl, &variant_include_cumulative_popcounts))) {
      goto DefineSets_ret_NOMEM;
    }
    FillCumulativePopcounts(variant_include, raw_variant_ctl, variant_include_cumulative_popcounts);

    reterr = InitTokenStream(fname, MAXV(max_thread_ct - 1, 1), &tks);
    if (unlikely(reterr)) {
      goto DefineSets_ret_TKSTREAM_FAIL;
    }
    // Pass 1 sizes everything; pass 2 fills it.  A variant ID that is not in
    // the current variant set is ignored, as in PLINK 1.x, and so is a set
    // that ends up empty.
    uintptr_t kept_set_ct = 0;
    uintptr_t max_set_name_blen = 0;
    uintptr_t membership_ct = 0;
    for (uint32_t pass_idx = 0; pass_idx != 2; ++pass_idx) {
      uint32_t* membership_offsets = nullptr;
      uint32_t* memberships = nullptr;
      char* set_names = nullptr;
      if (pass_idx) {
        if (merged_set_name) {
          kept_set_ct = 1;
          max_set_name_blen = strlen(merged_set_name) + 1;
        }
        if (unlikely(!kept_set_ct)) {
          logerrprintfww("Error: --set file '%s' contains no set that --subset/--set-names keeps.\n", fname);
          goto DefineSets_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(bigstack_end_alloc_u32(kept_set_ct + 1, &membership_offsets) ||
                     bigstack_end_alloc_u32(membership_ct, &memberships) ||
                     bigstack_end_alloc_c(kept_set_ct * max_set_name_blen, &set_names))) {
          goto DefineSets_ret_NOMEM;
        }
        vsp->set_names = set_names;
        vsp->max_set_name_blen = max_set_name_blen;
        membership_offsets[0] = 0;
        reterr = TokenRewind(&tks);
        if (unlikely(reterr)) {
          goto DefineSets_ret_TKSTREAM_FAIL;
        }
      }
      // 0: between sets, 1: inside a set being kept, 2: inside a skipped set
      uint32_t in_set = 0;
      uintptr_t set_idx = 0;
      uintptr_t cur_membership_ct = 0;
      while (1) {
        char* shard_boundaries[2];
        reterr = TksNext(&tks, 1, shard_boundaries);
        if (reterr) {
          break;
        }
        char* shard_iter = shard_boundaries[0];
        char* shard_end = shard_boundaries[1];
        while (1) {
          shard_iter = FirstPostspaceBounded(shard_iter, shard_end);
          if (shard_iter == shard_end) {
            break;
          }
          char* token_end = CurTokenEnd(shard_iter);
          const uint32_t token_slen = token_end - shard_iter;
          if ((token_slen == 3) && memequal_sk(shard_iter, "END")) {
            if (unlikely(!in_set)) {
              logerrprintfww("Error: Unmatched END in --set file '%s'.\n", fname);
              goto DefineSets_ret_MALFORMED_INPUT;
            }
            if ((in_set == 1) && (!merged_set_name)) {
              if (!pass_idx) {
                ++kept_set_ct;
              } else {
                membership_offsets[++set_idx] = cur_membership_ct;
              }
            }
            in_set = 0;
          } else if (!in_set) {
            const char token_end_char = *token_end;
            *token_end = '\0';
            in_set = SetNameWanted(shard_iter, token_slen, sip->setnames_flattened)? 1 : 2;
            if ((in_set == 1) && (!merged_set_name)) {
              if (!pass_idx) {
                if (token_slen >= max_set_name_blen) {
                  max_set_name_blen = token_slen + 1;
                }
              } else {
                memcpy(&(set_names[set_idx * max_set_name_blen]), shard_iter, token_slen + 1);
              }
            }
            *token_end = token_end_char;
          } else if (in_set == 1) {
            const char token_end_char = *token_end;
            *token_end = '\0';
            uint32_t cur_llidx;
            const uint32_t variant_uidx = VariantIdDupHtableFind(shard_iter, variant_ids, variant_id_htable, htable_dup_base, token_slen, variant_id_htable_size, max_variant_id_slen, &cur_llidx);
            *token_end = token_end_char;
            if (variant_uidx != UINT32_MAX) {
              if (!pass_idx) {
                ++membership_ct;
              } else {
                memberships[cur_membership_ct] = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx);
              }
              ++cur_membership_ct;
            }
          }
          shard_iter = token_end;
        }
      }
      if (unlikely(reterr != kPglRetEof)) {
        goto DefineSets_ret_TKSTREAM_FAIL;
      }
      reterr = kPglRetSuccess;
      if (unlikely(in_set)) {
        logerrprintfww("Error: --set file '%s' ends in the middle of a set; every set needs a closing END.\n", fname);
        goto DefineSets_ret_MALFORMED_INPUT;
      }
      if (pass_idx && merged_set_name) {
        memcpy(set_names, merged_set_name, max_set_name_blen);
        membership_offsets[1] = cur_membership_ct;
      }
      if (pass_idx) {
        // Now that the file has been read, the hash table can go, and the
        // setdefs can take its place.
        BigstackReset(htable_mark);
        uintptr_t* member_bitvec;
        if (unlikely(bigstack_end_alloc_w(variant_ctl + 1, &member_bitvec) ||
                     bigstack_alloc_u32p(kept_set_ct, &(vsp->setdefs)))) {
          goto DefineSets_ret_NOMEM;
        }
        unsigned char* alloc_base = g_bigstack_base;
        unsigned char* alloc_end = g_bigstack_end;
        for (uintptr_t cur_set_idx = 0; cur_set_idx != kept_set_ct; ++cur_set_idx) {
          ZeroWArr(variant_ctl + 1, member_bitvec);
          const uint32_t offset_start = membership_offsets[cur_set_idx];
          const uint32_t offset_end = membership_offsets[cur_set_idx + 1];
          for (uint32_t uii = offset_start; uii != offset_end; ++uii) {
            SetBit(memberships[uii], member_bitvec);
          }
          vsp->setdefs[cur_set_idx] = SaveSetdef(member_bitvec, variant_ct, &alloc_base, alloc_end);
          if (unlikely(!vsp->setdefs[cur_set_idx])) {
            goto DefineSets_ret_NOMEM;
          }
        }
        BigstackBaseSet(alloc_base);
        vsp->set_ct = kept_set_ct;
      }
    }
    if (CleanupTokenStream3(fname, &tks, &reterr)) {
      goto DefineSets_ret_1;
    }
    // The names have to outlive the temporaries above them.
    char* set_names_final;
    if (unlikely(bigstack_alloc_c(vsp->set_ct * vsp->max_set_name_blen, &set_names_final))) {
      goto DefineSets_ret_NOMEM;
    }
    memcpy(set_names_final, vsp->set_names, vsp->set_ct * vsp->max_set_name_blen);
    vsp->set_names = set_names_final;
    logprintf("--set: %" PRIuPTR " set%s loaded, %" PRIuPTR " membership%s in all.\n", vsp->set_ct, (vsp->set_ct == 1)? "" : "s", membership_ct, (membership_ct == 1)? "" : "s");
    bigstack_mark = g_bigstack_base;
  }
  while (0) {
  DefineSets_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  DefineSets_ret_TKSTREAM_FAIL:
    TokenStreamErrPrint(fname, &tks);
    reterr = TokenStreamErrcode(&tks);
    break;
  DefineSets_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  DefineSets_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 DefineSets_ret_1:
  CleanupTokenStream2(fname, &tks, &reterr);
  BigstackEndReset(bigstack_end_mark);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr WriteSetList(const VariantSets* vsp, const uintptr_t* variant_include, const char* const* variant_ids, uint32_t variant_ct, SetFlags flags, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t output_zst = (flags / kfSetWriteListZs) & 1;
    OutnameZstSet(".set", output_zst, outname_end);
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 64;
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WriteSetList_ret_1;
    }
    const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
    uintptr_t* member_bitvec;
    uint32_t* variant_idx_to_uidx;
    if (unlikely(bigstack_alloc_w(variant_ctl, &member_bitvec) ||
                 bigstack_alloc_u32(variant_ct, &variant_idx_to_uidx))) {
      goto WriteSetList_ret_NOMEM;
    }
    {
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        variant_idx_to_uidx[variant_idx] = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      }
    }
    const uintptr_t set_ct = vsp->set_ct;
    const uintptr_t max_set_name_blen = vsp->max_set_name_blen;
    for (uintptr_t set_idx = 0; set_idx != set_ct; ++set_idx) {
      cswritep = strcpya(cswritep, &(vsp->set_names[set_idx * max_set_name_blen]));
      AppendBinaryEoln(&cswritep);
      UnpackSetdef(vsp->setdefs[set_idx], variant_ct, member_bitvec);
      uintptr_t variant_idx_base = 0;
      uintptr_t member_bits = member_bitvec[0];
      const uint32_t cur_member_ct = PopcountWords(member_bitvec, variant_ctl);
      for (uint32_t uii = 0; uii != cur_member_ct; ++uii) {
        const uintptr_t variant_idx = BitIter1(member_bitvec, &variant_idx_base, &member_bits);
        cswritep = strcpya(cswritep, variant_ids[variant_idx_to_uidx[variant_idx]]);
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto WriteSetList_ret_WRITE_FAIL;
        }
      }
      cswritep = strcpya_k(cswritep, "END");
      AppendBinaryEoln(&cswritep);
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto WriteSetList_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto WriteSetList_ret_WRITE_FAIL;
    }
    logprintfww("--write-set: %" PRIuPTR " set%s written to %s .\n", set_ct, (set_ct == 1)? "" : "s", outname);
  }
  while (0) {
  WriteSetList_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WriteSetList_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 WriteSetList_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr LoadAndSortIntervalBed(const char* fname, const ChrInfo* cip, const char* sorted_subset_ids, uint32_t zero_based, uint32_t border_extend, uintptr_t subset_ct, uintptr_t max_subset_id_blen, uint32_t max_thread_ct, uintptr_t* gene_ct_ptr, char** gene_names_ptr, uintptr_t* max_gene_id_blen_ptr, uintptr_t** chr_bounds_ptr, uint32_t*** genedefs_ptr, uintptr_t* chr_max_gene_ct_ptr) {
  // --clump-range[0]; in plink 1.9, also --annotate and --gene-report
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    reterr = InitTextStreamEx(fname, 1, kMaxLongLine, kTextStreamBlenFast, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto LoadAndSortIntervalBed_ret_TSTREAM_FAIL;
    }
    uintptr_t gene_ct = 0;
    uintptr_t max_gene_id_blen = 0;
    uint64_t* range_sort_buf;
    MakeSetRange** gene_arr;
    reterr = LoadIntervalBed(cip, nullptr, sorted_subset_ids, fname, zero_based, 1, border_extend, 0, 0, subset_ct, max_subset_id_blen, &txs, &gene_ct, gene_names_ptr, &max_gene_id_blen, &range_sort_buf, &gene_arr);
    if (unlikely(reterr)) {
      goto LoadAndSortIntervalBed_ret_1;
    }
    const char* gene_names = *gene_names_ptr;
    const uint32_t chr_idx_end = cip->max_code + 1 + cip->name_ct;
    if (bigstack_alloc_w(chr_idx_end + 1, chr_bounds_ptr) ||
        bigstack_alloc_u32p(gene_ct, genedefs_ptr)) {
      goto LoadAndSortIntervalBed_ret_NOMEM;
    }
    uintptr_t* chr_bounds = *chr_bounds_ptr;
    chr_bounds[0] = 0;
    uint32_t** genedefs = *genedefs_ptr;
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = BigstackEndRoundedDown();
    uintptr_t chr_max_gene_ct = 0;
    uint32_t chr_idx = 0;
    for (uintptr_t gene_idx = 0; gene_idx != gene_ct; ++gene_idx) {
      const char* chrprefixed_gene_name = &(gene_names[gene_idx * max_gene_id_blen]);
      uint32_t new_chr_idx = 0;
      for (uint32_t uii = 0; uii != kMaxChrCodeDigits - 1; ++uii) {
        new_chr_idx += chrprefixed_gene_name[uii] - 48;
        new_chr_idx *= 10;
      }
      // Last prefix character must be nonnumeric to prevent weird natural-sort
      // interaction, so it's offset by 33 instead of 48.
      new_chr_idx += chrprefixed_gene_name[kMaxChrCodeDigits - 1] - 33;
      if (chr_idx < S_CAST(uint32_t, new_chr_idx)) {
        const uintptr_t chr_gene_ct = gene_idx - chr_bounds[chr_idx];
        if (chr_gene_ct > chr_max_gene_ct) {
          chr_max_gene_ct = chr_gene_ct;
        }
        do {
          chr_bounds[++chr_idx] = gene_idx;
        } while (chr_idx < S_CAST(uint32_t, new_chr_idx));
      }
      MakeSetRange* msr_tmp = gene_arr[gene_idx];
      uint32_t range_ct = 0;
      while (msr_tmp) {
        range_sort_buf[range_ct++] = (S_CAST(uint64_t, msr_tmp->uidx_start) << 32) | S_CAST(uint64_t, msr_tmp->uidx_end);
        msr_tmp = msr_tmp->next;
      }
      if (!range_ct) {
        genedefs[gene_idx] = R_CAST(uint32_t*, tmp_alloc_base);
        tmp_alloc_base = &(tmp_alloc_base[16]);
        if (tmp_alloc_end - tmp_alloc_base < 0) {
          goto LoadAndSortIntervalBed_ret_NOMEM;
        }
        genedefs[gene_idx][0] = 0;
        continue;
      }
      // Sort and merge intervals.  (This logic may belong in plink2_common.)
      STD_SORT_PAR_UNSEQ(range_ct, u64cmp, range_sort_buf);
      uint64_t range_write_entry = range_sort_buf[0];
      uint32_t range_read_idx = 1;
      uint64_t range_read_entry;
      for (; range_read_idx != range_ct; ++range_read_idx) {
        range_read_entry = range_sort_buf[range_read_idx];
        // explicit S_CAST to communicate intentional truncation.
        const uint32_t range_read_first_uidx = S_CAST(uint32_t, range_read_entry >> 32);
        if (range_read_first_uidx <= S_CAST(uint32_t, range_write_entry)) {
          break;
        }
        range_write_entry = range_read_entry;
      }
      uint32_t range_write_idx = range_read_idx;
      if (range_read_idx != range_ct) {
        --range_write_idx;
        uint32_t range_read_first_uidx;
        goto LoadAndSortIntervalBed_merge_start;
        for (; range_read_idx != range_ct; ++range_read_idx) {
          range_read_entry = range_sort_buf[range_read_idx];
          range_read_first_uidx = S_CAST(uint32_t, range_read_entry >> 32);
          if (range_read_first_uidx <= S_CAST(uint32_t, range_write_entry)) {
          LoadAndSortIntervalBed_merge_start:
            ;
            const uint32_t range_read_last_uidx = S_CAST(uint32_t, range_read_entry);
            if (range_read_last_uidx > S_CAST(uint32_t, range_write_entry)) {
              range_write_entry = (range_write_entry & 0xffffffff00000000LLU) | S_CAST(uint64_t, range_read_last_uidx);
            }
          } else {
            range_sort_buf[range_write_idx++] = range_write_entry;
            range_write_entry = range_read_entry;
          }
        }
        range_sort_buf[range_write_idx++] = range_write_entry;
      }

      const uintptr_t genedef_alloc_size = RoundUpPow2((range_write_idx * 2 + 1) * sizeof(int32_t), 16);
      uint32_t* genedef_iter = R_CAST(uint32_t*, tmp_alloc_base);
      tmp_alloc_base = &(tmp_alloc_base[genedef_alloc_size]);
      if (tmp_alloc_end - tmp_alloc_base < 0) {
        goto LoadAndSortIntervalBed_ret_NOMEM;
      }
      genedefs[gene_idx] = genedef_iter;
      *genedef_iter++ = range_write_idx;
      for (uint32_t range_idx = 0; range_idx != range_write_idx; ++range_idx) {
        const uint64_t range_entry = range_sort_buf[range_idx];
        *genedef_iter++ = S_CAST(uint32_t, range_entry >> 32);
        *genedef_iter++ = S_CAST(uint32_t, range_entry);
      }
    }

    BigstackBaseSet(tmp_alloc_base);

    const uintptr_t chr_gene_ct = gene_ct - chr_bounds[chr_idx];
    if (chr_gene_ct > chr_max_gene_ct) {
      chr_max_gene_ct = chr_gene_ct;
    }
    while (chr_idx < chr_idx_end) {
      chr_bounds[++chr_idx] = gene_ct;
    }
    *gene_ct_ptr = gene_ct;
    *max_gene_id_blen_ptr = max_gene_id_blen;
    *chr_max_gene_ct_ptr = chr_max_gene_ct;
  }
  while (0) {
  LoadAndSortIntervalBed_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadAndSortIntervalBed_ret_TSTREAM_FAIL:
    TextStreamErrPrint(fname, &txs);
    break;
  }
 LoadAndSortIntervalBed_ret_1:
  CleanupTextStream2(fname, &txs, &reterr);
  BigstackEndReset(bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}
#endif
