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

#include <string.h>

#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_simd.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "plink2_cmdline.h"
#include "plink2_compress_stream.h"
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

void InitGeneReport(GeneReportInfo* grip) {
  grip->report_fname = nullptr;
  grip->glist_fname = nullptr;
  grip->subset_fname = nullptr;
  grip->chr_field = nullptr;
  grip->pos_field = nullptr;
  grip->id_field = nullptr;
  grip->p_field = nullptr;
  // UINT32_MAX = unset; --gene-list-border cannot exceed 0x7ffffffe.
  grip->border = UINT32_MAX;
  grip->flags = kfGeneReport0;
}

void CleanupGeneReport(GeneReportInfo* grip) {
  free_cond(grip->report_fname);
  free_cond(grip->glist_fname);
  free_cond(grip->subset_fname);
  free_cond(grip->chr_field);
  free_cond(grip->pos_field);
  free_cond(grip->id_field);
  free_cond(grip->p_field);
}

// Loads a whitespace-delimited ID file into a strcmp-sorted, deduplicated
// fixed-width box at the end of the bigstack, so that bsearch_strbox() works
// on it.
PglErr LoadSortedIdBox(const char* fname, const char* file_descrip, uint32_t max_thread_ct, char** sorted_ids_ptr, uintptr_t* id_ct_ptr, uintptr_t* max_id_blen_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    reterr = SizeAndInitTextStream(fname, bigstack_left() / 4, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto LoadSortedIdBox_ret_TSTREAM_FAIL;
    }
    uintptr_t id_ct = 0;
    uintptr_t max_id_blen = 0;
    while (1) {
      const char* line_iter = TextGet(&txs);
      if (!line_iter) {
        break;
      }
      while (!IsEolnKns(*line_iter)) {
        const char* token_end = CurTokenEnd(line_iter);
        const uintptr_t slen = token_end - line_iter;
        if (unlikely(slen > kMaxIdSlen)) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s IDs are limited to " MAX_ID_SLEN_STR " characters.\n", file_descrip);
          goto LoadSortedIdBox_ret_MALFORMED_INPUT_WW;
        }
        if (slen >= max_id_blen) {
          max_id_blen = slen + 1;
        }
        ++id_ct;
        line_iter = FirstNonTspace(token_end);
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto LoadSortedIdBox_ret_TSTREAM_FAIL;
    }
    if (unlikely(!id_ct)) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s file is empty.\n", file_descrip);
      goto LoadSortedIdBox_ret_MALFORMED_INPUT_WW;
    }
    char* sorted_ids;
    uint32_t* id_map;
    if (unlikely(bigstack_end_alloc_c(id_ct * max_id_blen, &sorted_ids) ||
                 bigstack_alloc_u32(id_ct, &id_map))) {
      goto LoadSortedIdBox_ret_NOMEM;
    }
    reterr = TextRewind(&txs);
    if (unlikely(reterr)) {
      goto LoadSortedIdBox_ret_TSTREAM_FAIL;
    }
    uintptr_t id_idx = 0;
    while (id_idx != id_ct) {
      const char* line_iter = TextGet(&txs);
      if (unlikely(!line_iter)) {
        goto LoadSortedIdBox_ret_REWIND_FAIL;
      }
      while (!IsEolnKns(*line_iter)) {
        const char* token_end = CurTokenEnd(line_iter);
        const uintptr_t slen = token_end - line_iter;
        if (unlikely(id_idx == id_ct)) {
          goto LoadSortedIdBox_ret_REWIND_FAIL;
        }
        memcpyx(&(sorted_ids[id_idx * max_id_blen]), line_iter, slen, '\0');
        ++id_idx;
        line_iter = FirstNonTspace(token_end);
      }
    }
    if (unlikely(SortStrboxIndexed(id_ct, max_id_blen, 0, sorted_ids, id_map))) {
      goto LoadSortedIdBox_ret_NOMEM;
    }
    // Remove duplicates, so that bsearch_strbox() results are unambiguous.
    uintptr_t write_idx = 1;
    for (uintptr_t read_idx = 1; read_idx != id_ct; ++read_idx) {
      const char* cur_id = &(sorted_ids[read_idx * max_id_blen]);
      if (strequal_overread(&(sorted_ids[(write_idx - 1) * max_id_blen]), cur_id)) {
        continue;
      }
      if (write_idx != read_idx) {
        strcpy(&(sorted_ids[write_idx * max_id_blen]), cur_id);
      }
      ++write_idx;
    }
    *sorted_ids_ptr = sorted_ids;
    *id_ct_ptr = write_idx;
    *max_id_blen_ptr = max_id_blen;
  }
  while (0) {
  LoadSortedIdBox_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadSortedIdBox_ret_TSTREAM_FAIL:
    TextStreamErrPrint(file_descrip, &txs);
    break;
  LoadSortedIdBox_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, file_descrip);
    reterr = kPglRetRewindFail;
    break;
  LoadSortedIdBox_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  }
  CleanupTextStream2(file_descrip, &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

// Recovers the chromosome index from the fixed-width prefix LoadIntervalBed()
// prepends to each set name.  Mirrors the decoder in
// LoadAndSortIntervalBed().
uint32_t GenePrefixToChrIdx(const char* chrprefixed_gene_name) {
  uint32_t chr_idx = 0;
  for (uint32_t uii = 0; uii != kMaxChrCodeDigits - 1; ++uii) {
    chr_idx += ctou32(chrprefixed_gene_name[uii]) - 48;
    chr_idx *= 10;
  }
  // Last prefix character is offset by 33 instead of 48; see
  // LoadAndSortIntervalBed().
  return chr_idx + ctou32(chrprefixed_gene_name[kMaxChrCodeDigits - 1]) - 33;
}

// Each saved report line is stored as an 8-byte-aligned record:
//   [0..7]: ln(p-value), or kLnPvalError when the report has no P column
//   [8..11]: bp coordinate
//   [12..15]: variant ID length
//   [16..]: variant ID, not null-terminated, padded to a multiple of 8 bytes
CONSTI32(kGeneReportRecordHeaderSize, 16);

HEADER_INLINE uintptr_t GeneReportRecordSize(uintptr_t id_slen) {
  return RoundUpPow2(kGeneReportRecordHeaderSize + id_slen, 8);
}

PglErr GeneReport(const GeneReportInfo* grip, const ChrInfo* cip, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  const char* report_fname = grip->report_fname;
  uintptr_t line_idx = 0;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  CompressStreamState css;
  PreinitTextStream(&txs);
  PreinitCstream(&css);
  {
    char* sorted_subset_ids = nullptr;
    uintptr_t subset_ct = 0;
    uintptr_t max_subset_id_blen = 0;
    if (grip->subset_fname) {
      reterr = LoadSortedIdBox(grip->subset_fname, "--gene-subset file", max_thread_ct, &sorted_subset_ids, &subset_ct, &max_subset_id_blen);
      if (unlikely(reterr)) {
        goto GeneReport_ret_1;
      }
    }
    // The border is applied to the variant interval rather than to the gene
    // ranges, so that the reported gene boundaries and DIST values stay
    // relative to the unextended gene.
    uintptr_t gene_ct;
    char* gene_names;
    uintptr_t max_gene_id_blen;
    uintptr_t* chr_bounds;
    uint32_t** genedefs;
    uintptr_t chr_max_gene_ct;
    reterr = LoadAndSortIntervalBed(grip->glist_fname, cip, sorted_subset_ids, (grip->flags / kfGeneReport0based) & 1, 0, subset_ct, max_subset_id_blen, max_thread_ct, &gene_ct, &gene_names, &max_gene_id_blen, &chr_bounds, &genedefs, &chr_max_gene_ct);
    if (unlikely(reterr)) {
      goto GeneReport_ret_1;
    }
    BigstackEndReset(bigstack_end_mark);
    if (unlikely(gene_ct > 0x80000000LLU)) {
      snprintf(g_logbuf, kLogbufSize, "Error: Too many genes in %s (--gene-report can only handle 2147483648).\n", grip->glist_fname);
      goto GeneReport_ret_MALFORMED_INPUT_WW;
    }
    if (!gene_ct) {
      logerrputs("Warning: No genes remain for --gene-report.\n");
    }

    // gene_names is sorted by (chromosome index, gene name); the output is
    // sorted by (gene name, chromosome index), so build the remapping tables.
    uint32_t* gene_nameidx_to_chridx;
    uint32_t* gene_chridx_to_nameidx;
    if (unlikely(bigstack_alloc_u32(gene_ct, &gene_nameidx_to_chridx) ||
                 bigstack_alloc_u32(gene_ct, &gene_chridx_to_nameidx))) {
      goto GeneReport_ret_NOMEM;
    }
    if (gene_ct) {
      // Sort key: gene name, then a space, then the chromosome-index prefix.
      const uintptr_t sort_key_blen = max_gene_id_blen + kMaxChrCodeDigits + 1;
      char* sort_keys;
      if (unlikely(bigstack_alloc_c(gene_ct * sort_key_blen, &sort_keys))) {
        goto GeneReport_ret_NOMEM;
      }
      for (uintptr_t gene_idx = 0; gene_idx != gene_ct; ++gene_idx) {
        const char* chrprefixed_gene_name = &(gene_names[gene_idx * max_gene_id_blen]);
        char* key_iter = strcpyax(&(sort_keys[gene_idx * sort_key_blen]), &(chrprefixed_gene_name[kMaxChrCodeDigits]), ' ');
        memcpyx(key_iter, chrprefixed_gene_name, kMaxChrCodeDigits, '\0');
        gene_nameidx_to_chridx[gene_idx] = gene_idx;
      }
      if (unlikely(SortStrboxIndexed(gene_ct, sort_key_blen, 1, sort_keys, gene_nameidx_to_chridx))) {
        goto GeneReport_ret_NOMEM;
      }
      for (uintptr_t name_idx = 0; name_idx != gene_ct; ++name_idx) {
        gene_chridx_to_nameidx[gene_nameidx_to_chridx[name_idx]] = name_idx;
      }
      BigstackReset(sort_keys);
    }

    reterr = SizeAndInitTextStream(report_fname, bigstack_left() / 8, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto GeneReport_ret_TSTREAM_FAIL;
    }
    char* header_start;
    do {
      ++line_idx;
      header_start = TextGet(&txs);
      if (unlikely(!header_start)) {
        reterr = TextStreamRawErrcode(&txs);
        if (reterr == kPglRetEof) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", report_fname);
          goto GeneReport_ret_MALFORMED_INPUT_WW;
        }
        goto GeneReport_ret_TSTREAM_FAIL;
      }
    } while (strequal_k_unsafe(header_start, "##"));
    if (*header_start == '#') {
      ++header_start;
    }
    const GeneReportFlags flags = grip->flags;
    // [0] = CHROM, [1] = POS, [2] = ID, [3] = P
    const char* col_search_order[4];
    col_search_order[0] = grip->chr_field? grip->chr_field : "CHROM\0CHR\0";
    col_search_order[1] = grip->pos_field? grip->pos_field : "POS\0BP\0";
    col_search_order[2] = grip->id_field? grip->id_field : "ID\0SNP\0";
    col_search_order[3] = grip->p_field? grip->p_field : "P\0UNADJ\0";
    uint32_t col_skips[4];
    uint32_t col_types[4];
    uint32_t relevant_col_ct;
    uint32_t found_type_bitset;
    reterr = SearchHeaderLine(header_start, col_search_order, "gene-report", 4, &relevant_col_ct, &found_type_bitset, col_skips, col_types);
    if (unlikely(reterr)) {
      goto GeneReport_ret_1;
    }
    if (unlikely((found_type_bitset & 7) != 7)) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s must have chromosome, bp coordinate, and variant ID columns.\n", report_fname);
      goto GeneReport_ret_INCONSISTENT_INPUT_WW;
    }
    const uint32_t p_col_present = (found_type_bitset >> 3) & 1;
    if (unlikely((!p_col_present) && (ln_pfilter != kLnPvalError))) {
      snprintf(g_logbuf, kLogbufSize, "Error: --pfilter requires a p-value column in %s.\n", report_fname);
      goto GeneReport_ret_INCONSISTENT_INPUT_WW;
    }

    // Report lines that match at least one gene are stored at the bottom of
    // the arena (growing upwards); (gene, line) matches are stored at the top
    // (8-byte entries, growing downwards).  A single line can match at most
    // chr_max_gene_ct genes, so that much headroom is kept in reserve.
    unsigned char* record_iter = g_bigstack_base;
    uint64_t* match_list_end = R_CAST(uint64_t*, BigstackEndRoundedDown());
    uint64_t* match_list = match_list_end;
    const uint32_t border = (grip->border == UINT32_MAX)? 0 : grip->border;
    uint32_t max_id_slen = 1;
    uintptr_t saved_line_ct = 0;
    uintptr_t skipped_chr_ct = 0;
    while (1) {
      ++line_idx;
      char* line_start = TextGet(&txs);
      if (!line_start) {
        break;
      }
      const char* token_ptrs[4];
      uint32_t token_slens[4];
      if (unlikely(!TokenLexK0(line_start, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, report_fname);
        goto GeneReport_ret_MALFORMED_INPUT_WW;
      }
      const uint32_t chr_idx = GetChrCode(token_ptrs[0], cip, token_slens[0]);
      if (IsI32Neg(chr_idx)) {
        ++skipped_chr_ct;
        continue;
      }
      if (!IsSet(cip->chr_mask, chr_idx)) {
        continue;
      }
      uint32_t variant_bp;
      if (unlikely(ScanUintDefcap(token_ptrs[1], &variant_bp))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, report_fname);
        goto GeneReport_ret_MALFORMED_INPUT_WW;
      }
      double ln_pval = kLnPvalError;
      if (p_col_present) {
        const char* pval_str = token_ptrs[3];
        if (!ScantokLn(pval_str, &ln_pval)) {
          const uint32_t pval_slen = token_slens[3];
          if (IsNanStr(pval_str, pval_slen)) {
            ln_pval = kLnPvalError;
          } else if (likely(strequal_k(pval_str, "INF", pval_slen))) {
            ln_pval = kLnNormalMin;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid p-value on line %" PRIuPTR " of %s.\n", line_idx, report_fname);
            goto GeneReport_ret_MALFORMED_INPUT_WW;
          }
        }
        if (ln_pval > ln_pfilter) {
          continue;
        }
      }
      const uintptr_t id_slen = token_slens[2];
      const uintptr_t record_size = GeneReportRecordSize(id_slen);
      if (unlikely(S_CAST(uintptr_t, R_CAST(unsigned char*, match_list) - record_iter) < record_size + chr_max_gene_ct * sizeof(int64_t))) {
        goto GeneReport_ret_NOMEM;
      }
      const uint32_t bp_start = (variant_bp > border)? (variant_bp - border) : 0;
      const uint32_t bp_end = (variant_bp > UINT32_MAX - border)? UINT32_MAX : (variant_bp + border);
      const uintptr_t gene_idx_end = chr_bounds[chr_idx + 1];
      uint64_t* match_list_stop = match_list;
      for (uintptr_t gene_idx = chr_bounds[chr_idx]; gene_idx != gene_idx_end; ++gene_idx) {
        if (IntervalInSetdef(genedefs[gene_idx], bp_start, bp_end)) {
          *(--match_list) = (S_CAST(uint64_t, gene_chridx_to_nameidx[gene_idx]) << 32) | saved_line_ct;
        }
      }
      if (match_list == match_list_stop) {
        continue;
      }
      if (unlikely(saved_line_ct == 0x100000000LLU)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Too many matching lines in %s (--gene-report can only handle 4294967296).\n", report_fname);
        goto GeneReport_ret_MALFORMED_INPUT_WW;
      }
      memcpy(record_iter, &ln_pval, sizeof(double));
      memcpy(&(record_iter[8]), &variant_bp, sizeof(int32_t));
      const uint32_t id_slen_u32 = id_slen;
      memcpy(&(record_iter[12]), &id_slen_u32, sizeof(int32_t));
      memcpy(&(record_iter[kGeneReportRecordHeaderSize]), token_ptrs[2], id_slen);
      record_iter = &(record_iter[record_size]);
      if (id_slen > max_id_slen) {
        max_id_slen = id_slen;
      }
      ++saved_line_ct;
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto GeneReport_ret_TSTREAM_FAIL;
    }
    if (unlikely(CleanupTextStream2(report_fname, &txs, &reterr))) {
      goto GeneReport_ret_1;
    }
    if (skipped_chr_ct) {
      logerrprintfww("Warning: %" PRIuPTR " line%s in %s skipped due to unrecognized chromosome code%s.\n", skipped_chr_ct, (skipped_chr_ct == 1)? "" : "s", report_fname, (skipped_chr_ct == 1)? "" : "s");
    }

    // Saved-line index -> record lookup table, placed just below the match
    // list since record_iter is not aligned for pointers.
    if (unlikely(S_CAST(uintptr_t, R_CAST(unsigned char*, match_list) - record_iter) < saved_line_ct * sizeof(intptr_t))) {
      goto GeneReport_ret_NOMEM;
    }
    unsigned char** line_lookup = &(R_CAST(unsigned char**, match_list)[-S_CAST(intptr_t, saved_line_ct)]);
    {
      unsigned char* record_scan = g_bigstack_base;
      for (uintptr_t saved_line_idx = 0; saved_line_idx != saved_line_ct; ++saved_line_idx) {
        line_lookup[saved_line_idx] = record_scan;
        uint32_t id_slen;
        memcpy(&id_slen, &(record_scan[12]), sizeof(int32_t));
        record_scan = &(record_scan[GeneReportRecordSize(id_slen)]);
      }
    }
    const uintptr_t match_ct = match_list_end - match_list;
    STD_SORT(match_ct, u64cmp, match_list);
    BigstackBaseSet(record_iter);
    BigstackEndSet(line_lookup);

    const uint32_t output_zst = (flags / kfGeneReportZs) & 1;
    OutnameZstSet(".gene.report", output_zst, outname_end);
    const uintptr_t overflow_buf_size = kCompressStreamBlock + max_gene_id_blen + max_id_slen + 256;
    reterr = InitCstreamAlloc(outname, 0, output_zst, MAXV(max_thread_ct - 1, 1), overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto GeneReport_ret_1;
    }
    const uint32_t chrom_col = (flags / kfGeneReportColChrom) & 1;
    const uint32_t genepos_col = (flags / kfGeneReportColGenepos) & 1;
    const uint32_t genekb_col = (flags / kfGeneReportColGenekb) & 1;
    const uint32_t dist_col = (flags / kfGeneReportColDist) & 1;
    const uint32_t pos_col = (flags / kfGeneReportColPos) & 1;
    const uint32_t p_col = p_col_present && ((flags / kfGeneReportColP) & 1);
    *cswritep++ = '#';
    cswritep = strcpya_k(cswritep, "GENE");
    if (chrom_col) {
      cswritep = strcpya_k(cswritep, "\tCHROM");
    }
    if (genepos_col) {
      cswritep = strcpya_k(cswritep, "\tGENE_START\tGENE_END");
    }
    if (genekb_col) {
      cswritep = strcpya_k(cswritep, "\tGENE_KB");
    }
    if (dist_col) {
      cswritep = strcpya_k(cswritep, "\tDIST");
    }
    cswritep = strcpya_k(cswritep, "\tID");
    if (pos_col) {
      cswritep = strcpya_k(cswritep, "\tPOS");
    }
    if (p_col) {
      cswritep = strcpya_k(cswritep, "\tP");
    }
    AppendBinaryEoln(&cswritep);

    uintptr_t prev_name_idx = ~k0LU;
    const char* cur_gene_name = nullptr;
    uint32_t cur_gene_start = 0;
    uint32_t cur_gene_end = 0;
    uint32_t cur_chr_idx = 0;
    double cur_gene_kb = 0.0;
    for (uintptr_t match_idx = 0; match_idx != match_ct; ++match_idx) {
      const uint64_t cur_match = match_list[match_idx];
      const uintptr_t name_idx = cur_match >> 32;
      if (name_idx != prev_name_idx) {
        prev_name_idx = name_idx;
        const uintptr_t gene_idx = gene_nameidx_to_chridx[name_idx];
        const char* chrprefixed_gene_name = &(gene_names[gene_idx * max_gene_id_blen]);
        cur_gene_name = &(chrprefixed_gene_name[kMaxChrCodeDigits]);
        cur_chr_idx = GenePrefixToChrIdx(chrprefixed_gene_name);
        const uint32_t* genedef = genedefs[gene_idx];
        const uint32_t range_ct = genedef[0];
        cur_gene_start = genedef[1];
        cur_gene_end = genedef[2 * range_ct];
        uint32_t covered_bp_ct = 0;
        for (uint32_t range_idx = 0; range_idx != range_ct; ++range_idx) {
          covered_bp_ct += genedef[2 * range_idx + 2] - genedef[2 * range_idx + 1];
        }
        cur_gene_kb = u31tod(covered_bp_ct) * 0.001;
      }
      const unsigned char* record = line_lookup[S_CAST(uint32_t, cur_match)];
      double ln_pval;
      memcpy(&ln_pval, record, sizeof(double));
      uint32_t variant_bp;
      memcpy(&variant_bp, &(record[8]), sizeof(int32_t));
      uint32_t id_slen;
      memcpy(&id_slen, &(record[12]), sizeof(int32_t));
      cswritep = strcpya(cswritep, cur_gene_name);
      if (chrom_col) {
        *cswritep++ = '\t';
        cswritep = chrtoa(cip, cur_chr_idx, cswritep);
      }
      if (genepos_col) {
        *cswritep++ = '\t';
        cswritep = u32toa_x(cur_gene_start, '\t', cswritep);
        cswritep = u32toa(cur_gene_end - 1, cswritep);
      }
      if (genekb_col) {
        *cswritep++ = '\t';
        cswritep = dtoa_g(cur_gene_kb, cswritep);
      }
      if (dist_col) {
        *cswritep++ = '\t';
        cswritep = dtoa_g(S_CAST(double, S_CAST(int32_t, variant_bp) - S_CAST(int32_t, cur_gene_start)) * 0.001, cswritep);
      }
      *cswritep++ = '\t';
      cswritep = memcpya(cswritep, &(record[kGeneReportRecordHeaderSize]), id_slen);
      if (pos_col) {
        *cswritep++ = '\t';
        cswritep = u32toa(variant_bp, cswritep);
      }
      if (p_col) {
        *cswritep++ = '\t';
        cswritep = lntoa_g(MAXV(ln_pval, output_min_ln), cswritep);
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto GeneReport_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto GeneReport_ret_WRITE_FAIL;
    }
    logprintfww("--gene-report: %" PRIuPTR " gene-variant pair%s written to %s .\n", match_ct, (match_ct == 1)? "" : "s", outname);
  }
  while (0) {
  GeneReport_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GeneReport_ret_TSTREAM_FAIL:
    TextStreamErrPrint(report_fname, &txs);
    break;
  GeneReport_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GeneReport_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  GeneReport_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  }
 GeneReport_ret_1:
  CswriteCloseCond(&css, cswritep);
  CleanupTextStream2(report_fname, &txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

void InitAnnot(AnnotInfo* aip) {
  aip->report_fname = nullptr;
  aip->attrib_fname = nullptr;
  aip->ranges_fname = nullptr;
  aip->filter_fname = nullptr;
  aip->snps_fname = nullptr;
  aip->subset_fname = nullptr;
  aip->chr_field = nullptr;
  aip->pos_field = nullptr;
  aip->id_field = nullptr;
  aip->p_field = nullptr;
  // UINT32_MAX = unset; --annotate-border cannot exceed 0x7ffffffe.
  aip->border = UINT32_MAX;
  aip->flags = kfAnnot0;
}

void CleanupAnnot(AnnotInfo* aip) {
  free_cond(aip->report_fname);
  free_cond(aip->attrib_fname);
  free_cond(aip->ranges_fname);
  free_cond(aip->filter_fname);
  free_cond(aip->snps_fname);
  free_cond(aip->subset_fname);
  free_cond(aip->chr_field);
  free_cond(aip->pos_field);
  free_cond(aip->id_field);
  free_cond(aip->p_field);
}

uint32_t InSetdef(const uint32_t* setdef, uint32_t pos) {
  // Interval-list form only; that is what LoadAndSortIntervalBed() produces
  // when variant_bps is nullptr.
  const uint32_t range_ct = setdef[0];
  if (!range_ct) {
    return 0;
  }
  return LowerBoundNonemptyU32(&(setdef[1]), range_ct * 2, pos + 1) & 1;
}

uint32_t InSetdefDist(const uint32_t* setdef, uint32_t pos, uint32_t border, int32_t* dist_ptr) {
  // Returns 1 and sets *dist_ptr to the signed distance from pos to the
  // nearest interval boundary (0 when pos is inside an interval) if pos is
  // within border bp of the set; returns 0 otherwise.
  // Ties are broken in favor of negative distances, matching plink 1.07's
  // annot.cpp.
  const uint32_t range_ct = setdef[0];
  if (!range_ct) {
    return 0;
  }
  const uint32_t idx = LowerBoundNonemptyU32(&(setdef[1]), range_ct * 2, pos + 1);
  if (idx & 1) {
    *dist_ptr = 0;
    return 1;
  }
  if (!idx) {
    // Before the first interval.
    if (pos + border >= setdef[1]) {
      *dist_ptr = S_CAST(int32_t, pos) - S_CAST(int32_t, setdef[1]);
      return 1;
    }
    return 0;
  }
  // setdef[idx] is the end of the previous interval (exclusive).
  if (idx == range_ct * 2) {
    // After the last interval.
    if (setdef[idx] + border > pos) {
      *dist_ptr = S_CAST(int32_t, pos + 1 - setdef[idx]);
      return 1;
    }
    return 0;
  }
  // Between two intervals; setdef[idx + 1] is the start of the next one.
  if (setdef[idx] + border > pos) {
    int32_t dist = S_CAST(int32_t, pos + 1 - setdef[idx]);
    if (pos + S_CAST(uint32_t, dist) > setdef[idx + 1]) {
      dist = S_CAST(int32_t, pos) - S_CAST(int32_t, setdef[idx + 1]);
    }
    *dist_ptr = dist;
    return 1;
  }
  if (pos + border >= setdef[idx + 1]) {
    *dist_ptr = S_CAST(int32_t, pos) - S_CAST(int32_t, setdef[idx + 1]);
    return 1;
  }
  return 0;
}

// Open-addressed table of offsets into a variable-width, null-terminated blob.
// UINT32_MAX marks an empty slot.  Returns the offset of an existing match, or
// UINT32_MAX after inserting new_offset.
uint32_t AttrHtableAdd(const char* cur_id, uint32_t cur_id_slen, const char* blob, uint32_t htable_size, uint32_t new_offset, uint32_t* htable) {
  for (uint32_t hashval = Hashceil(cur_id, cur_id_slen, htable_size); ; ) {
    const uint32_t cur_entry = htable[hashval];
    if (cur_entry == UINT32_MAX) {
      htable[hashval] = new_offset;
      return UINT32_MAX;
    }
    const char* cur_str = &(blob[cur_entry]);
    if (memequal(cur_id, cur_str, cur_id_slen) && (!cur_str[cur_id_slen])) {
      return cur_entry;
    }
    if (++hashval == htable_size) {
      hashval = 0;
    }
  }
}

// Loads an attribute file: one line per variant, "<variant ID> <attribute>...".
// Lines with no attribute are ignored, as are variants missing from
// sorted_snplist when that is present.
//
// On success:
// * sorted_attr_ids is a natural-sorted, fixed-width box of the distinct
//   attribute names.
// * attr_var_ids is a strcmp-sorted, fixed-width box of the variant IDs, and
//   attr_var_id_map maps a position in it back to the corresponding row of
//   attr_bitfields.
// * attr_bitfields has one attr_ct-bit row per variant, in file order.
PglErr LoadAttribFile(const char* fname, const char* sorted_snplist, uintptr_t snplist_ct, uintptr_t max_snplist_id_blen, uint32_t max_thread_ct, char** sorted_attr_ids_ptr, uintptr_t* attr_ct_ptr, uintptr_t* max_attr_id_blen_ptr, char** attr_var_ids_ptr, uintptr_t* attr_var_ct_ptr, uintptr_t* max_attr_var_id_blen_ptr, uint32_t** attr_var_id_map_ptr, uintptr_t** attr_bitfields_ptr, uint32_t* max_onevar_attr_ct_ptr) {
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    reterr = SizeAndInitTextStream(fname, bigstack_left() / 8, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto LoadAttribFile_ret_TSTREAM_FAIL;
    }
    // Pass 1: size everything.
    uintptr_t var_ct = 0;
    uintptr_t max_var_id_blen = 2;
    uintptr_t max_attr_id_blen = 2;
    uintptr_t attr_token_ct = 0;
    uintptr_t attr_name_byte_ct = 0;
    uint32_t max_onevar_attr_ct = 0;
    while (1) {
      ++line_idx;
      const char* line_iter = TextGet(&txs);
      if (!line_iter) {
        break;
      }
      const char* var_id = FirstNonTspace(line_iter);
      if (IsEolnKns(*var_id)) {
        continue;
      }
      const char* var_id_end = CurTokenEnd(var_id);
      const uintptr_t var_id_slen = var_id_end - var_id;
      const char* attr_iter = FirstNonTspace(var_id_end);
      if (IsEolnKns(*attr_iter)) {
        continue;
      }
      if (snplist_ct && (bsearch_strbox(var_id, sorted_snplist, var_id_slen, max_snplist_id_blen, snplist_ct) == -1)) {
        continue;
      }
      if (unlikely(var_id_slen > kMaxIdSlen)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Variant ID on line %" PRIuPTR " of %s is longer than " MAX_ID_SLEN_STR " characters.\n", line_idx, fname);
        goto LoadAttribFile_ret_MALFORMED_INPUT_WW;
      }
      if (var_id_slen >= max_var_id_blen) {
        max_var_id_blen = var_id_slen + 1;
      }
      uint32_t cur_attr_ct = 0;
      do {
        const char* attr_end = CurTokenEnd(attr_iter);
        const uintptr_t attr_slen = attr_end - attr_iter;
        if (unlikely(attr_slen > kMaxIdSlen)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Attribute name on line %" PRIuPTR " of %s is longer than " MAX_ID_SLEN_STR " characters.\n", line_idx, fname);
          goto LoadAttribFile_ret_MALFORMED_INPUT_WW;
        }
        if (attr_slen >= max_attr_id_blen) {
          max_attr_id_blen = attr_slen + 1;
        }
        attr_name_byte_ct += attr_slen + 1;
        ++attr_token_ct;
        ++cur_attr_ct;
        attr_iter = FirstNonTspace(attr_end);
      } while (!IsEolnKns(*attr_iter));
      if (cur_attr_ct > max_onevar_attr_ct) {
        max_onevar_attr_ct = cur_attr_ct;
      }
      ++var_ct;
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto LoadAttribFile_ret_TSTREAM_FAIL;
    }
    if (unlikely(!var_ct)) {
      snprintf(g_logbuf, kLogbufSize, "Error: No usable lines in %s.\n", fname);
      goto LoadAttribFile_ret_MALFORMED_INPUT_WW;
    }
    if (unlikely(attr_token_ct > 0x7fffffff)) {
      snprintf(g_logbuf, kLogbufSize, "Error: Too many attribute entries in %s.\n", fname);
      goto LoadAttribFile_ret_MALFORMED_INPUT_WW;
    }

    // Pass 2: collect the distinct attribute names in an append-only blob.
    char* name_blob;
    if (unlikely(bigstack_end_alloc_c(attr_name_byte_ct, &name_blob))) {
      goto LoadAttribFile_ret_NOMEM;
    }
    uint32_t* htable;
    uint32_t htable_size;
    if (unlikely(HtableGoodSizeAlloc(attr_token_ct, bigstack_left() / 4, &htable, &htable_size))) {
      goto LoadAttribFile_ret_NOMEM;
    }
    SetAllU32Arr(htable_size, htable);
    uintptr_t blob_byte_ct = 0;
    uintptr_t attr_ct = 0;
    reterr = TextRewind(&txs);
    if (unlikely(reterr)) {
      goto LoadAttribFile_ret_TSTREAM_FAIL;
    }
    line_idx = 0;
    while (1) {
      ++line_idx;
      const char* line_iter = TextGet(&txs);
      if (!line_iter) {
        break;
      }
      const char* var_id = FirstNonTspace(line_iter);
      if (IsEolnKns(*var_id)) {
        continue;
      }
      const char* var_id_end = CurTokenEnd(var_id);
      const char* attr_iter = FirstNonTspace(var_id_end);
      if (IsEolnKns(*attr_iter)) {
        continue;
      }
      if (snplist_ct && (bsearch_strbox(var_id, sorted_snplist, var_id_end - var_id, max_snplist_id_blen, snplist_ct) == -1)) {
        continue;
      }
      do {
        const char* attr_end = CurTokenEnd(attr_iter);
        const uint32_t attr_slen = attr_end - attr_iter;
        if (AttrHtableAdd(attr_iter, attr_slen, name_blob, htable_size, blob_byte_ct, htable) == UINT32_MAX) {
          memcpyx(&(name_blob[blob_byte_ct]), attr_iter, attr_slen, '\0');
          blob_byte_ct += attr_slen + 1;
          ++attr_ct;
        }
        attr_iter = FirstNonTspace(attr_end);
      } while (!IsEolnKns(*attr_iter));
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto LoadAttribFile_ret_TSTREAM_FAIL;
    }
    BigstackReset(htable);

    char* sorted_attr_ids;
    if (unlikely(bigstack_alloc_c(attr_ct * max_attr_id_blen, &sorted_attr_ids))) {
      goto LoadAttribFile_ret_NOMEM;
    }
    {
      const char* blob_iter = name_blob;
      for (uintptr_t attr_idx = 0; attr_idx != attr_ct; ++attr_idx) {
        const uint32_t slen = strlen(blob_iter);
        memcpyx(&(sorted_attr_ids[attr_idx * max_attr_id_blen]), blob_iter, slen, '\0');
        blob_iter = &(blob_iter[slen + 1]);
      }
    }
    BigstackEndReset(bigstack_end_mark);
    {
      uint32_t* ignored_id_map;
      if (unlikely(bigstack_alloc_u32(attr_ct, &ignored_id_map))) {
        goto LoadAttribFile_ret_NOMEM;
      }
      if (unlikely(SortStrboxIndexed(attr_ct, max_attr_id_blen, 1, sorted_attr_ids, ignored_id_map))) {
        goto LoadAttribFile_ret_NOMEM;
      }
      BigstackReset(ignored_id_map);
    }

    // Pass 3: fill the per-variant bitfields.
    const uintptr_t attr_ctl = BitCtToWordCt(attr_ct);
    char* attr_var_ids;
    uint32_t* attr_var_id_map;
    uintptr_t* attr_bitfields;
    if (unlikely(bigstack_alloc_c(var_ct * max_var_id_blen, &attr_var_ids) ||
                 bigstack_alloc_u32(var_ct, &attr_var_id_map) ||
                 bigstack_calloc_w(var_ct * attr_ctl, &attr_bitfields))) {
      goto LoadAttribFile_ret_NOMEM;
    }
    if (unlikely(HtableGoodSizeAlloc(attr_ct, bigstack_left() / 4, &htable, &htable_size))) {
      goto LoadAttribFile_ret_NOMEM;
    }
    PopulateStrboxHtable(sorted_attr_ids, attr_ct, max_attr_id_blen, htable_size, htable);
    reterr = TextRewind(&txs);
    if (unlikely(reterr)) {
      goto LoadAttribFile_ret_TSTREAM_FAIL;
    }
    line_idx = 0;
    uintptr_t var_idx = 0;
    while (var_idx != var_ct) {
      ++line_idx;
      const char* line_iter = TextGet(&txs);
      if (unlikely(!line_iter)) {
        goto LoadAttribFile_ret_REWIND_FAIL;
      }
      const char* var_id = FirstNonTspace(line_iter);
      if (IsEolnKns(*var_id)) {
        continue;
      }
      const char* var_id_end = CurTokenEnd(var_id);
      const uintptr_t var_id_slen = var_id_end - var_id;
      const char* attr_iter = FirstNonTspace(var_id_end);
      if (IsEolnKns(*attr_iter)) {
        continue;
      }
      if (snplist_ct && (bsearch_strbox(var_id, sorted_snplist, var_id_slen, max_snplist_id_blen, snplist_ct) == -1)) {
        continue;
      }
      memcpyx(&(attr_var_ids[var_idx * max_var_id_blen]), var_id, var_id_slen, '\0');
      attr_var_id_map[var_idx] = var_idx;
      uintptr_t* cur_bitfield = &(attr_bitfields[var_idx * attr_ctl]);
      do {
        const char* attr_end = CurTokenEnd(attr_iter);
        const uint32_t attr_slen = attr_end - attr_iter;
        const uint32_t attr_idx = StrboxHtableFindNnt(attr_iter, sorted_attr_ids, htable, max_attr_id_blen, attr_slen, htable_size);
        if (unlikely(attr_idx == UINT32_MAX)) {
          goto LoadAttribFile_ret_REWIND_FAIL;
        }
        SetBit(attr_idx, cur_bitfield);
        attr_iter = FirstNonTspace(attr_end);
      } while (!IsEolnKns(*attr_iter));
      ++var_idx;
    }
    // No TextStreamErrcode2() check here: this loop stops as soon as every
    // counted variant has been read, which is usually before end-of-file, and
    // that helper treats "no error yet" the same as a failure.
    BigstackReset(htable);
    if (unlikely(SortStrboxIndexed(var_ct, max_var_id_blen, 0, attr_var_ids, attr_var_id_map))) {
      goto LoadAttribFile_ret_NOMEM;
    }
    *sorted_attr_ids_ptr = sorted_attr_ids;
    *attr_ct_ptr = attr_ct;
    *max_attr_id_blen_ptr = max_attr_id_blen;
    *attr_var_ids_ptr = attr_var_ids;
    *attr_var_ct_ptr = var_ct;
    *max_attr_var_id_blen_ptr = max_var_id_blen;
    *attr_var_id_map_ptr = attr_var_id_map;
    *attr_bitfields_ptr = attr_bitfields;
    *max_onevar_attr_ct_ptr = max_onevar_attr_ct;
  }
  while (0) {
  LoadAttribFile_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadAttribFile_ret_TSTREAM_FAIL:
    TextStreamErrPrint(fname, &txs);
    break;
  LoadAttribFile_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, fname);
    reterr = kPglRetRewindFail;
    break;
  LoadAttribFile_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  }
  CleanupTextStream2(fname, &txs, &reterr);
  BigstackEndReset(bigstack_end_mark);
  return reterr;
}

PglErr Annotate(const AnnotInfo* aip, const ChrInfo* cip, double ln_pfilter, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  const char* report_fname = aip->report_fname;
  uintptr_t line_idx = 0;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  CompressStreamState css;
  PreinitTextStream(&txs);
  PreinitCstream(&css);
  {
    const AnnotFlags flags = aip->flags;
    const uint32_t zero_based = (flags / kfAnnot0based) & 1;
    const uint32_t block01 = (flags / kfAnnotBlock) & 1;
    const uint32_t prune = (flags / kfAnnotPrune) & 1;
    const uint32_t range_dist = !(flags & kfAnnotMinimal);
    const uint32_t track_distance = (flags / kfAnnotDistance) & 1;
    const uint32_t border = (aip->border == UINT32_MAX)? 0 : aip->border;
    const char* no_annot_str = (flags & kfAnnotNa)? "NA" : ".";

    char* sorted_snplist = nullptr;
    uintptr_t snplist_ct = 0;
    uintptr_t max_snplist_id_blen = 0;
    if (aip->snps_fname) {
      reterr = LoadSortedIdBox(aip->snps_fname, "--annotate snps= file", max_thread_ct, &sorted_snplist, &snplist_ct, &max_snplist_id_blen);
      if (unlikely(reterr)) {
        goto Annotate_ret_1;
      }
    }
    // sorted_snplist has to outlive the interval loads and the attribute load,
    // so the end-of-arena mark used to release the subset IDs is taken after
    // it, not at function entry.
    unsigned char* snplist_end_mark = g_bigstack_end;
    char* sorted_subset_ids = nullptr;
    uintptr_t subset_ct = 0;
    uintptr_t max_subset_id_blen = 0;
    if (aip->subset_fname) {
      reterr = LoadSortedIdBox(aip->subset_fname, "--annotate subset= file", max_thread_ct, &sorted_subset_ids, &subset_ct, &max_subset_id_blen);
      if (unlikely(reterr)) {
        goto Annotate_ret_1;
      }
    }

    uintptr_t range_ct = 0;
    char* range_names = nullptr;
    uintptr_t max_range_name_blen = 0;
    uintptr_t* chr_bounds = nullptr;
    uint32_t** rangedefs = nullptr;
    uintptr_t chr_max_range_ct = 0;
    if (aip->ranges_fname) {
      reterr = LoadAndSortIntervalBed(aip->ranges_fname, cip, sorted_subset_ids, zero_based, 0, subset_ct, max_subset_id_blen, max_thread_ct, &range_ct, &range_names, &max_range_name_blen, &chr_bounds, &rangedefs, &chr_max_range_ct);
      if (unlikely(reterr)) {
        goto Annotate_ret_1;
      }
      if (unlikely(range_ct > 0x7fffffff)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Too many intervals in %s (--annotate can only handle 2147483647).\n", aip->ranges_fname);
        goto Annotate_ret_MALFORMED_INPUT_WW;
      }
    }
    uintptr_t* filter_chr_bounds = nullptr;
    uint32_t** filter_rangedefs = nullptr;
    if (aip->filter_fname) {
      uintptr_t ignored_range_ct;
      char* ignored_names;
      uintptr_t ignored_max_blen;
      uintptr_t ignored_chr_max;
      reterr = LoadAndSortIntervalBed(aip->filter_fname, cip, nullptr, zero_based, 0, 0, 0, max_thread_ct, &ignored_range_ct, &ignored_names, &ignored_max_blen, &filter_chr_bounds, &filter_rangedefs, &ignored_chr_max);
      if (unlikely(reterr)) {
        goto Annotate_ret_1;
      }
    }
    BigstackEndReset(snplist_end_mark);

    char* sorted_attr_ids = nullptr;
    uintptr_t attr_ct = 0;
    uintptr_t max_attr_id_blen = 0;
    char* attr_var_ids = nullptr;
    uintptr_t attr_var_ct = 0;
    uintptr_t max_attr_var_id_blen = 0;
    uint32_t* attr_var_id_map = nullptr;
    uintptr_t* attr_bitfields = nullptr;
    uint32_t max_onevar_attr_ct = 0;
    if (aip->attrib_fname) {
      reterr = LoadAttribFile(aip->attrib_fname, sorted_snplist, snplist_ct, max_snplist_id_blen, max_thread_ct, &sorted_attr_ids, &attr_ct, &max_attr_id_blen, &attr_var_ids, &attr_var_ct, &max_attr_var_id_blen, &attr_var_id_map, &attr_bitfields, &max_onevar_attr_ct);
      if (unlikely(reterr)) {
        goto Annotate_ret_1;
      }
    }
    const uintptr_t attr_ctl = BitCtToWordCt(attr_ct);

    // In 'block' mode, the single ANNOT column is replaced by one 0/1 column
    // per distinct annotation.  The columns are the natural-sorted union of
    // the interval names and the attribute names.
    uint32_t* range_to_col = nullptr;
    uint32_t* attr_to_col = nullptr;
    char* block_col_names = nullptr;
    uintptr_t block_col_ct = 0;
    uintptr_t max_block_col_name_blen = 0;
    if (block01) {
      // merged_ct can be zero, when every annotation source turned out to be
      // empty; plink 1.9 emits the report unchanged in that case.
      const uintptr_t merged_ct = range_ct + attr_ct;
      max_block_col_name_blen = MAXV(max_attr_id_blen, (max_range_name_blen > kMaxChrCodeDigits)? (max_range_name_blen - kMaxChrCodeDigits) : 1);
      // block_col_names, range_to_col and attr_to_col outlive the sort
      // scratch, so they have to be allocated underneath it.
      char* merged_names;
      uint32_t* merged_id_map;
      if (unlikely(bigstack_alloc_c(merged_ct * max_block_col_name_blen, &block_col_names) ||
                   bigstack_alloc_u32(range_ct, &range_to_col) ||
                   bigstack_alloc_u32(attr_ct, &attr_to_col) ||
                   bigstack_alloc_c(merged_ct * max_block_col_name_blen, &merged_names) ||
                   bigstack_alloc_u32(merged_ct, &merged_id_map))) {
        goto Annotate_ret_NOMEM;
      }
      for (uintptr_t range_idx = 0; range_idx != range_ct; ++range_idx) {
        strcpy(&(merged_names[range_idx * max_block_col_name_blen]), &(range_names[range_idx * max_range_name_blen + kMaxChrCodeDigits]));
        merged_id_map[range_idx] = range_idx;
      }
      for (uintptr_t attr_idx = 0; attr_idx != attr_ct; ++attr_idx) {
        strcpy(&(merged_names[(range_ct + attr_idx) * max_block_col_name_blen]), &(sorted_attr_ids[attr_idx * max_attr_id_blen]));
        merged_id_map[range_ct + attr_idx] = range_ct + attr_idx;
      }
      if (unlikely(SortStrboxIndexed(merged_ct, max_block_col_name_blen, 1, merged_names, merged_id_map))) {
        goto Annotate_ret_NOMEM;
      }
      // Deduplicate; an interval and an attribute can share a name, in which
      // case they share a column.
      for (uintptr_t merged_idx = 0; merged_idx != merged_ct; ++merged_idx) {
        const char* cur_name = &(merged_names[merged_idx * max_block_col_name_blen]);
        if ((!block_col_ct) || (!strequal_overread(&(block_col_names[(block_col_ct - 1) * max_block_col_name_blen]), cur_name))) {
          strcpy(&(block_col_names[block_col_ct * max_block_col_name_blen]), cur_name);
          ++block_col_ct;
        }
        const uint32_t orig_idx = merged_id_map[merged_idx];
        if (orig_idx < range_ct) {
          range_to_col[orig_idx] = block_col_ct - 1;
        } else {
          attr_to_col[orig_idx - range_ct] = block_col_ct - 1;
        }
      }
      BigstackReset(merged_names);
    }

    reterr = SizeAndInitTextStream(report_fname, bigstack_left() / 8, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto Annotate_ret_TSTREAM_FAIL;
    }
    char* line_start;
    do {
      ++line_idx;
      line_start = TextGet(&txs);
      if (unlikely(!line_start)) {
        reterr = TextStreamRawErrcode(&txs);
        if (reterr == kPglRetEof) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", report_fname);
          goto Annotate_ret_MALFORMED_INPUT_WW;
        }
        goto Annotate_ret_TSTREAM_FAIL;
      }
    } while (strequal_k_unsafe(line_start, "##"));
    const char* header_start = line_start;
    if (*header_start == '#') {
      ++header_start;
    }
    const uint32_t need_pos = (range_ct != 0) || (filter_chr_bounds != nullptr);
    const uint32_t need_var_id = (attr_ct != 0) || (snplist_ct != 0);
    const char* col_search_order[4];
    col_search_order[0] = need_pos? (aip->chr_field? aip->chr_field : "CHROM\0CHR\0") : "";
    col_search_order[1] = need_pos? (aip->pos_field? aip->pos_field : "POS\0BP\0") : "";
    col_search_order[2] = need_var_id? (aip->id_field? aip->id_field : "ID\0SNP\0") : "";
    col_search_order[3] = aip->p_field? aip->p_field : "P\0UNADJ\0";
    uint32_t col_skips[4];
    uint32_t col_types[4];
    uint32_t relevant_col_ct;
    uint32_t found_type_bitset;
    reterr = SearchHeaderLine(header_start, col_search_order, "annotate", 4, &relevant_col_ct, &found_type_bitset, col_skips, col_types);
    if (unlikely(reterr)) {
      goto Annotate_ret_1;
    }
    if (unlikely(need_pos && ((found_type_bitset & 3) != 3))) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s must have chromosome and bp coordinate columns.\n", report_fname);
      goto Annotate_ret_INCONSISTENT_INPUT_WW;
    }
    if (unlikely(need_var_id && (!(found_type_bitset & 4)))) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s must have a variant ID column.\n", report_fname);
      goto Annotate_ret_INCONSISTENT_INPUT_WW;
    }
    const uint32_t p_col_present = (found_type_bitset >> 3) & 1;
    if (unlikely((!p_col_present) && (ln_pfilter != kLnPvalError))) {
      snprintf(g_logbuf, kLogbufSize, "Error: --pfilter requires a p-value column in %s.\n", report_fname);
      goto Annotate_ret_INCONSISTENT_INPUT_WW;
    }

    const uint32_t output_zst = (flags / kfAnnotZs) & 1;
    OutnameZstSet(".annot", output_zst, outname_end);
    const uintptr_t overflow_buf_size = kCompressStreamBlock + MAXV(max_block_col_name_blen, MAXV(max_range_name_blen, max_attr_id_blen)) + 256;
    reterr = InitCstreamAlloc(outname, 0, output_zst, MAXV(max_thread_ct - 1, 1), overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto Annotate_ret_1;
    }
    // Header: the original columns, whitespace runs collapsed to single tabs,
    // then the new columns.
    *cswritep++ = '#';
    {
      const char* token_iter = header_start;
      uint32_t is_first = 1;
      while (!IsEolnKns(*token_iter)) {
        const char* token_end = CurTokenEnd(token_iter);
        if (!is_first) {
          *cswritep++ = '\t';
        }
        is_first = 0;
        if (unlikely(CsputsStd(token_iter, token_end - token_iter, &css, &cswritep))) {
          goto Annotate_ret_WRITE_FAIL;
        }
        token_iter = FirstNonTspace(token_end);
      }
    }
    if (track_distance) {
      cswritep = strcpya_k(cswritep, "\tDIST\tSGN");
    }
    if (!block01) {
      cswritep = strcpya_k(cswritep, "\tANNOT");
    } else {
      for (uintptr_t col_idx = 0; col_idx != block_col_ct; ++col_idx) {
        *cswritep++ = '\t';
        const char* cur_name = &(block_col_names[col_idx * max_block_col_name_blen]);
        if (unlikely(CsputsStd(cur_name, strlen(cur_name), &css, &cswritep))) {
          goto Annotate_ret_WRITE_FAIL;
        }
      }
    }
    AppendBinaryEoln(&cswritep);
    if (unlikely(Cswrite(&css, &cswritep))) {
      goto Annotate_ret_WRITE_FAIL;
    }

    uintptr_t* block_hits = nullptr;
    if (block01) {
      if (unlikely(bigstack_alloc_w(BitCtToWordCt(block_col_ct), &block_hits))) {
        goto Annotate_ret_NOMEM;
      }
    }
    uintptr_t annot_row_ct = 0;
    uintptr_t total_row_ct = 0;
    while (1) {
      ++line_idx;
      line_start = TextGet(&txs);
      if (!line_start) {
        break;
      }
      const char* token_ptrs[4];
      uint32_t token_slens[4];
      if (unlikely(!TokenLexK0(line_start, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, report_fname);
        goto Annotate_ret_MALFORMED_INPUT_WW;
      }
      uint32_t chr_idx = 0;
      uint32_t cur_bp = 0;
      if (need_pos) {
        chr_idx = GetChrCode(token_ptrs[0], cip, token_slens[0]);
        if (IsI32Neg(chr_idx) || (!IsSet(cip->chr_mask, chr_idx))) {
          continue;
        }
        if (unlikely(ScanUintDefcap(token_ptrs[1], &cur_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, report_fname);
          goto Annotate_ret_MALFORMED_INPUT_WW;
        }
        if (filter_chr_bounds) {
          const uintptr_t filter_idx_end = filter_chr_bounds[chr_idx + 1];
          uintptr_t filter_idx = filter_chr_bounds[chr_idx];
          for (; filter_idx != filter_idx_end; ++filter_idx) {
            if (InSetdef(filter_rangedefs[filter_idx], cur_bp)) {
              break;
            }
          }
          if (filter_idx == filter_idx_end) {
            continue;
          }
        }
      }
      if (snplist_ct && (bsearch_strbox(token_ptrs[2], sorted_snplist, token_slens[2], max_snplist_id_blen, snplist_ct) == -1)) {
        continue;
      }
      if (p_col_present) {
        const char* pval_str = token_ptrs[3];
        double ln_pval;
        if (!ScantokLn(pval_str, &ln_pval)) {
          const uint32_t pval_slen = token_slens[3];
          if (IsNanStr(pval_str, pval_slen)) {
            ln_pval = kLnPvalError;
          } else if (likely(strequal_k(pval_str, "INF", pval_slen))) {
            ln_pval = kLnNormalMin;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid p-value on line %" PRIuPTR " of %s.\n", line_idx, report_fname);
            goto Annotate_ret_MALFORMED_INPUT_WW;
          }
        }
        if (ln_pval > ln_pfilter) {
          continue;
        }
      }

      // Attribute bitfield for this variant, if any.
      const uintptr_t* cur_attr_bits = nullptr;
      if (attr_ct) {
        const int32_t sorted_idx = bsearch_strbox(token_ptrs[2], attr_var_ids, token_slens[2], max_attr_var_id_blen, attr_var_ct);
        if (sorted_idx != -1) {
          cur_attr_bits = &(attr_bitfields[S_CAST(uintptr_t, attr_var_id_map[S_CAST(uint32_t, sorted_idx)]) * attr_ctl]);
        }
      }

      uint32_t abs_min_dist = UINT32_MAX;
      int32_t min_dist = 0;
      uint32_t at_least_one_annot = 0;
      if (block01) {
        ZeroWArr(BitCtToWordCt(block_col_ct), block_hits);
      }
      // Pass over the intervals on this chromosome.
      uintptr_t range_idx_end = 0;
      uintptr_t range_idx = 0;
      if (range_ct) {
        range_idx = chr_bounds[chr_idx];
        range_idx_end = chr_bounds[chr_idx + 1];
      }
      if (!border) {
        for (; range_idx != range_idx_end; ++range_idx) {
          if (InSetdef(rangedefs[range_idx], cur_bp)) {
            at_least_one_annot = 1;
            abs_min_dist = 0;
            if (block01) {
              SetBit(range_to_col[range_idx], block_hits);
            }
          }
        }
      } else {
        for (; range_idx != range_idx_end; ++range_idx) {
          int32_t cur_dist;
          if (InSetdefDist(rangedefs[range_idx], cur_bp, border, &cur_dist)) {
            at_least_one_annot = 1;
            const uint32_t cur_abs_dist = abs_i32(cur_dist);
            if (cur_abs_dist < abs_min_dist) {
              abs_min_dist = cur_abs_dist;
              min_dist = cur_dist;
            }
            if (block01) {
              SetBit(range_to_col[range_idx], block_hits);
            }
          }
        }
      }
      if (cur_attr_bits) {
        for (uintptr_t widx = 0; widx != attr_ctl; ++widx) {
          uintptr_t cur_word = cur_attr_bits[widx];
          if (!cur_word) {
            continue;
          }
          at_least_one_annot = 1;
          if (block01) {
            do {
              const uint32_t attr_idx = widx * kBitsPerWord + ctzw(cur_word);
              SetBit(attr_to_col[attr_idx], block_hits);
              cur_word &= cur_word - 1;
            } while (cur_word);
          }
        }
      }
      if (at_least_one_annot) {
        ++annot_row_ct;
      } else if (prune) {
        continue;
      }
      ++total_row_ct;

      // Pass through the original columns, with whitespace runs collapsed to
      // single tabs.
      {
        const char* token_iter = FirstNonTspace(line_start);
        uint32_t is_first = 1;
        while (!IsEolnKns(*token_iter)) {
          const char* token_end = CurTokenEnd(token_iter);
          if (!is_first) {
            *cswritep++ = '\t';
          }
          is_first = 0;
          if (unlikely(CsputsStd(token_iter, token_end - token_iter, &css, &cswritep))) {
            goto Annotate_ret_WRITE_FAIL;
          }
          token_iter = FirstNonTspace(token_end);
        }
      }
      if (track_distance) {
        *cswritep++ = '\t';
        if (abs_min_dist != UINT32_MAX) {
          cswritep = dtoa_g(u31tod(abs_min_dist) * 0.001, cswritep);
          *cswritep++ = '\t';
          if (!abs_min_dist) {
            cswritep = strcpya(cswritep, no_annot_str);
          } else {
            *cswritep++ = (min_dist > 0)? '+' : '-';
          }
        } else {
          cswritep = strcpya(cswritep, no_annot_str);
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, no_annot_str);
        }
      }
      if (!block01) {
        *cswritep++ = '\t';
        if (!at_least_one_annot) {
          cswritep = strcpya(cswritep, no_annot_str);
        } else {
          uint32_t is_first_annot = 1;
          if (range_ct) {
            for (uintptr_t rid = chr_bounds[chr_idx]; rid != chr_bounds[chr_idx + 1]; ++rid) {
              int32_t cur_dist = 0;
              if (border) {
                if (!InSetdefDist(rangedefs[rid], cur_bp, border, &cur_dist)) {
                  continue;
                }
              } else if (!InSetdef(rangedefs[rid], cur_bp)) {
                continue;
              }
              if (!is_first_annot) {
                *cswritep++ = '|';
              }
              is_first_annot = 0;
              const char* cur_name = &(range_names[rid * max_range_name_blen + kMaxChrCodeDigits]);
              if (unlikely(CsputsStd(cur_name, strlen(cur_name), &css, &cswritep))) {
                goto Annotate_ret_WRITE_FAIL;
              }
              if (range_dist) {
                if (!cur_dist) {
                  cswritep = strcpya_k(cswritep, "(0)");
                } else {
                  *cswritep++ = '(';
                  if (cur_dist > 0) {
                    *cswritep++ = '+';
                  }
                  cswritep = dtoa_g(S_CAST(double, cur_dist) * 0.001, cswritep);
                  cswritep = strcpya_k(cswritep, "kb)");
                }
              }
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto Annotate_ret_WRITE_FAIL;
              }
            }
          }
          if (cur_attr_bits) {
            uintptr_t attr_idx_base = 0;
            uintptr_t cur_bits = cur_attr_bits[0];
            const uint32_t cur_attr_ct = PopcountWords(cur_attr_bits, attr_ctl);
            for (uint32_t attr_ii = 0; attr_ii != cur_attr_ct; ++attr_ii) {
              const uintptr_t attr_idx = BitIter1(cur_attr_bits, &attr_idx_base, &cur_bits);
              if (!is_first_annot) {
                *cswritep++ = '|';
              }
              is_first_annot = 0;
              const char* cur_name = &(sorted_attr_ids[attr_idx * max_attr_id_blen]);
              if (unlikely(CsputsStd(cur_name, strlen(cur_name), &css, &cswritep))) {
                goto Annotate_ret_WRITE_FAIL;
              }
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto Annotate_ret_WRITE_FAIL;
              }
            }
          }
        }
      } else {
        for (uintptr_t col_idx = 0; col_idx != block_col_ct; ++col_idx) {
          *cswritep++ = '\t';
          *cswritep++ = IsSet(block_hits, col_idx)? '1' : '0';
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto Annotate_ret_WRITE_FAIL;
          }
        }
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto Annotate_ret_WRITE_FAIL;
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto Annotate_ret_TSTREAM_FAIL;
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto Annotate_ret_WRITE_FAIL;
    }
    logprintfww("--annotate: %" PRIuPTR " row%s written to %s , %" PRIuPTR " with at least one annotation.\n", total_row_ct, (total_row_ct == 1)? "" : "s", outname, annot_row_ct);
  }
  while (0) {
  Annotate_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Annotate_ret_TSTREAM_FAIL:
    TextStreamErrPrint(report_fname, &txs);
    break;
  Annotate_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  Annotate_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  Annotate_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  }
 Annotate_ret_1:
  CswriteCloseCond(&css, cswritep);
  CleanupTextStream2(report_fname, &txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}
#endif
