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


#include "plink2_fasta.h"

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr RefFromFaContig(const uintptr_t* variant_include, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const ChrInfo* cip, const char* seqbuf, uint32_t force, uint32_t chr_fo_idx, uint32_t variant_uidx_last, uint32_t bp_end, STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uintptr_t* nonref_flags, uint32_t* changed_ct_ptr, uint32_t* validated_ct_ptr, uint32_t* downgraded_ct_ptr) {
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], &variant_uidx_base, &cur_bits);
  uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
  if (variant_bps[variant_uidx_last] >= bp_end) {
    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
    if (!force) {
      char* write_iter = strcpya_k(g_logbuf, "Error: Contig '");
      write_iter = chrtoa(cip, chr_idx, write_iter);
      snprintf(write_iter, kLogbufSize - kMaxIdSlen - 32, "' in --fa file is too short; it is likely to be mismatched with your data. Add the 'force' modifier if this wasn't a mistake, and you just want to mark all reference alleles past the end as provisional.\n");
      WordWrapB(0);
      logerrputsb();
      return kPglRetInconsistentInput;
    } else {
      char* write_iter = strcpya_k(g_logbuf, "Warning: Contig '");
      write_iter = chrtoa(cip, chr_idx, write_iter);
      snprintf(write_iter, kLogbufSize - kMaxIdSlen - 32, "' in --fa file is too short; it is likely to be mismatched with your data.\n");
      WordWrapB(0);
      logerrputsb();
    }
    uint32_t offset = CountSortedSmallerU32(&(variant_bps[variant_uidx]), variant_uidx_last - variant_uidx, bp_end);

    const uint32_t chr_vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    // set all bits in [variant_uidx + offset, chr_vidx_end), and count how
    // many relevant bits were set.
    const uint32_t widx_full_end = chr_vidx_end / kBitsPerWord;
    const uint32_t trailing_bit_ct = chr_vidx_end % kBitsPerWord;
    const uint32_t chr_vidx_truncate = variant_uidx + offset;
    uint32_t widx = chr_vidx_truncate / kBitsPerWord;
    uint32_t downgraded_incr;
    if (widx == widx_full_end) {
      const uintptr_t cur_mask = (k1LU << trailing_bit_ct) - (k1LU << (chr_vidx_truncate % kBitsPerWord));
      downgraded_incr = PopcountWord(variant_include[widx] & (~nonref_flags[widx]) & cur_mask);
      nonref_flags[widx] |= cur_mask;
    } else {
      downgraded_incr = 0;
      if (chr_vidx_truncate % kBitsPerWord) {
        const uintptr_t cur_mask = (~k0LU) << (chr_vidx_truncate % kBitsPerWord);
        downgraded_incr = PopcountWord(variant_include[widx] & (~nonref_flags[widx]) & cur_mask);
        nonref_flags[widx] |= cur_mask;
        ++widx;
      }
      assert(widx > widx_full_end);
      for (; widx != widx_full_end; ++widx) {
        downgraded_incr += PopcountWord(variant_include[widx] & (~nonref_flags[widx]));
        nonref_flags[widx] = ~k0LU;
      }
      if (trailing_bit_ct) {
        const uintptr_t cur_mask = (k1LU << trailing_bit_ct) - k1LU;
        downgraded_incr += PopcountWord(variant_include[widx] & (~nonref_flags[widx]) & cur_mask);
        nonref_flags[widx] |= cur_mask;
      }
    }
    *downgraded_ct_ptr += downgraded_incr;
    if (!offset) {
      return kPglRetSuccess;
    }
    // We know that variant_bps[variant_uidx + offset - 1] < bp_end,
    // variant_bps[variant_uidx + offset] >= bp_end, and that at least one
    // variant_include[] bit is set in [variant_uidx, variant_uidx + offset)
    // (since offset > 0).  Find the last such bit.
    variant_uidx_last = FindLast1BitBefore(variant_include, chr_vidx_truncate);
  }
  uint32_t changed_ct = *changed_ct_ptr;
  uint32_t validated_ct = *validated_ct_ptr;
  uint32_t allele_ct = 2;
  for (; ; variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits)) {
    const uint32_t cur_bp = variant_bps[variant_uidx];
    uintptr_t allele_idx_offset_base = variant_uidx * 2;
    if (allele_idx_offsets) {
      allele_idx_offset_base = allele_idx_offsets[variant_uidx];
      allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
    }
    const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
    const char* cur_ref = &(seqbuf[cur_bp]);
    int32_t consistent_allele_idx = -1;
    for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
      const char* cur_allele = cur_alleles[allele_idx];
      const uint32_t cur_allele_slen = strlen(cur_allele);
      if (strcaseequal(cur_allele, cur_ref, cur_allele_slen)) {
        if (consistent_allele_idx != -1) {
          // Multiple alleles could be ref (this always happens for deletions).
          // Don't try to do anything.
          consistent_allele_idx = -2;
          break;
        }
        consistent_allele_idx = allele_idx;
      }
    }
    if (consistent_allele_idx >= 0) {
      if (consistent_allele_idx) {
        if ((!IsSet(nonref_flags, variant_uidx)) && (!force)) {
          const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
          char* write_iter = strcpya_k(g_logbuf, "Error: --ref-from-fa wants to change reference allele assignment at ");
          write_iter = chrtoa(cip, chr_idx, write_iter);
          *write_iter++ = ':';
          write_iter = u32toa(cur_bp, write_iter);
          snprintf(write_iter, kLogbufSize - kMaxIdSlen - 128, ", but it's marked as 'known'. Add the 'force' modifier to force this change through.\n");
          WordWrapB(0);
          logerrputsb();
          return kPglRetInconsistentInput;
        }
        R_CAST(DoubleAlleleCode*, refalt1_select)[variant_uidx] = consistent_allele_idx;
        ++changed_ct;
      } else {
        ++validated_ct;
      }
      ClearBit(variant_uidx, nonref_flags);
    } else if ((consistent_allele_idx == -1) && (!IsSet(nonref_flags, variant_uidx))) {
      // okay to have multiple matches, but not zero matches
      if (!force) {
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* write_iter = strcpya_k(g_logbuf, "Error: Reference allele at ");
        write_iter = chrtoa(cip, chr_idx, write_iter);
        *write_iter++ = ':';
        write_iter = u32toa(cur_bp, write_iter);
        snprintf(write_iter, kLogbufSize - kMaxIdSlen - 64, " is marked as 'known', but is inconsistent with .fa file. Add the 'force' modifier to downgrade it to provisional.\n");
        WordWrapB(0);
        logerrputsb();
        return kPglRetInconsistentInput;
      }
      SetBit(variant_uidx, nonref_flags);
      *downgraded_ct_ptr += 1;
    }
    if (variant_uidx == variant_uidx_last) {
      *changed_ct_ptr = changed_ct;
      *validated_ct_ptr = validated_ct;
      return kPglRetSuccess;
    }
  }
}

PglErr VNormalizeContig(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const ChrInfo* cip, const char* seqbuf, uint32_t chr_fo_idx, uint32_t variant_uidx_last, uint32_t bp_end, unsigned char** alloc_endp, UnsortedVar* vpos_sortstatusp, uint32_t* __restrict variant_bps, const char** allele_storage, uint32_t* __restrict nchanged_ct_ptr, char* nlist_flush, FILE* nlist_file, char** nlist_write_iterp, uint32_t* __restrict alen_buf) {
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], &variant_uidx_base, &cur_bits);
  uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
  if (variant_bps[variant_uidx_last] >= bp_end) {
    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
    char* write_iter = strcpya_k(g_logbuf, "Error: Contig '");
    write_iter = chrtoa(cip, chr_idx, write_iter);
    snprintf(write_iter, kLogbufSize - kMaxIdSlen - 32, "' in --fa file is too short; it is likely to be mismatched with your data.\n");
    WordWrapB(0);
    logerrputsb();
    return kPglRetInconsistentInput;
  }
  const char* missing_allele_str = &(g_one_char_strs[92]);
  const uint32_t input_missing_geno_code = ctou32(*g_input_missing_geno_ptr);
  unsigned char* alloc_base = g_bigstack_base;
  unsigned char* alloc_end = *alloc_endp;
  char* nlist_write_iter = *nlist_write_iterp;
  uint32_t nchanged_ct = *nchanged_ct_ptr;
  uint32_t allele_ct = 2;
  uint32_t prev_bp = 0;
  uint32_t is_unsorted = ((*vpos_sortstatusp) / kfUnsortedVarBp) & 1;
  uint32_t cur_bp;
  goto VNormalizeContig_loop_start;
  while (1) {
    if (!is_unsorted) {
      if (prev_bp > cur_bp) {
        is_unsorted = 1;
      }
      prev_bp = cur_bp;
    }
    if (variant_uidx == variant_uidx_last) {
      break;
    }
    variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
  VNormalizeContig_loop_start:
    cur_bp = variant_bps[variant_uidx];
    uintptr_t allele_idx_offset_base = variant_uidx * 2;
    if (allele_idx_offsets) {
      allele_idx_offset_base = allele_idx_offsets[variant_uidx];
      allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
    }
    const char** cur_alleles = &(allele_storage[allele_idx_offset_base]);
    {
      uint32_t aidx = 0;
      // Common case: if all alleles have length 1, variant must already be
      // normalized.
      for (; aidx != allele_ct; ++aidx) {
        if (cur_alleles[aidx][1]) {
          break;
        }
      }
      if (aidx == allele_ct) {
        continue;
      }
    }
    // Ignore missing allele if present, get lengths of each allele, check if
    // all right nucleotides match or all left nucleotides match
    uint32_t left_match = 0;
    uint32_t right_match = 0;
    uint32_t min_alen = 0;
    uint32_t missing_aidx = UINT32_MAX;
    for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
      if (cur_alleles[aidx] == missing_allele_str) {
        if (missing_aidx != UINT32_MAX) {
          const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
          char* write_iter = strcpya_k(g_logbuf, "Error: Variant at ");
          write_iter = chrtoa(cip, chr_idx, write_iter);
          *write_iter++ = ':';
          write_iter = u32toa(cur_bp, write_iter);
          snprintf(write_iter, kLogbufSize - kMaxIdSlen - 128, " has multiple missing alleles.\n");
          WordWrapB(0);
          logerrputsb();
          return kPglRetInconsistentInput;
        }
        missing_aidx = aidx;
        alen_buf[aidx] = 1;
        continue;
      }
      const char* cur_allele = cur_alleles[aidx];
      const uint32_t alen = strlen(cur_allele);
      alen_buf[aidx] = alen;
      const uint32_t first_code = ctou32(cur_allele[0]);
      // Special case: if first character of any allele is '<', don't attempt
      // to normalize.
      if (first_code == '<') {
        left_match = UINT32_MAX;
        right_match = UINT32_MAX;
        break;
      }
      const uint32_t last_code = ctou32(cur_allele[alen - 1]);
      if (!min_alen) {
        // must be first nonmissing allele
        left_match = first_code;
        right_match = last_code;
        min_alen = alen;
      } else {
        if (first_code != left_match) {
          left_match = UINT32_MAX;
        }
        if (last_code != right_match) {
          right_match = UINT32_MAX;
        }
        if (alen < min_alen) {
          min_alen = alen;
        }
      }
    }
    if (((left_match == UINT32_MAX) || (min_alen == 1)) && (right_match == UINT32_MAX)) {
      continue;
    }
    // Sanity check: verify alleles aren't all identical.
    const uint32_t first_aidx = (missing_aidx == 0)? 1 : 0;
    const uint32_t first_alen = alen_buf[first_aidx];
    for (uint32_t aidx = first_aidx + 1; ; ) {
      if (aidx == missing_aidx) {
        continue;
      }
      if ((alen_buf[aidx] != first_alen) || (!memequal(cur_alleles[first_aidx], cur_alleles[aidx], min_alen))) {
        break;
      }
      if (++aidx == allele_ct) {
        // probable todo: report ID instead?
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* write_iter = strcpya_k(g_logbuf, "Error: Variant at ");
        write_iter = chrtoa(cip, chr_idx, write_iter);
        *write_iter++ = ':';
        write_iter = u32toa(cur_bp, write_iter);
        snprintf(write_iter, kLogbufSize - kMaxIdSlen - 128, " has duplicate allele codes.\n");
        WordWrapB(0);
        logerrputsb();
        return kPglRetInconsistentInput;
      }
    }
    ++nchanged_ct;
    if (nlist_write_iter) {
      nlist_write_iter = strcpya(nlist_write_iter, variant_ids[variant_uidx]);
      AppendBinaryEoln(&nlist_write_iter);
      if (unlikely(fwrite_ck(nlist_flush, nlist_file, &nlist_write_iter))) {
        return kPglRetWriteFail;
      }
    }
    // Algorithm from Tan A, Abecasis GR, Kang HM (2015) Unified representation
    // of genetic variants:
    //   while (alleles end with the same nucleotide) {
    //     truncate rightmost
    //     if there's now an empty allele, extend alleles 1 base to the left
    //   }
    //   while (all alleles length > 1 and start with same nucleotide) {
    //     truncate leftmost
    //   }
    // https://genome.sph.umich.edu/wiki/Variant_Normalization also covers
    // this.
    const uint32_t rtrim_stop = cur_bp + min_alen - 1;
    if ((right_match == UINT32_MAX) || (!rtrim_stop)) {
      // Easy case: only left-trimming.  We separate this out because bigstack
      // allocation is never required.
      const uint32_t min_alen_m1 = min_alen - 1;
      uint32_t ltrim = 1;
      uint32_t ltrim_p1;
      for (; ltrim != min_alen_m1; ++ltrim) {
        const char cc = cur_alleles[first_aidx][ltrim];
        for (uint32_t aidx = first_aidx + 1; aidx != allele_ct; ++aidx) {
          if (aidx == missing_aidx) {
            continue;
          }
          if (cur_alleles[aidx][ltrim] != cc) {
            goto VNormalizeContig_ltrim1_done;
          }
        }
      }
    VNormalizeContig_ltrim1_done:
      ltrim_p1 = ltrim + 1;
      for (uint32_t aidx = first_aidx; aidx != allele_ct; ++aidx) {
        if (aidx == missing_aidx) {
          continue;
        }
        if (alen_buf[aidx] == ltrim_p1) {
          const uint32_t cur_code = ctou32(cur_alleles[aidx][ltrim]);
          if ((cur_code == 46) || (cur_code == input_missing_geno_code)) {
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            char* write_iter = strcpya_k(g_logbuf, "Error: Variant at ");
            write_iter = chrtoa(cip, chr_idx, write_iter);
            *write_iter++ = ':';
            write_iter = u32toa(cur_bp, write_iter);
            snprintf(write_iter, kLogbufSize - kMaxIdSlen - 128, " has an invalid normalized allele (conflicts with missing code).\n");
            WordWrapB(0);
            logerrputsb();
            return kPglRetInconsistentInput;
          }
          cur_alleles[aidx] = &(g_one_char_strs[2 * cur_code]);
        } else {
          cur_alleles[aidx] = &(cur_alleles[aidx][ltrim]);
        }
      }
      cur_bp += ltrim;
      variant_bps[variant_uidx] = cur_bp;
      continue;
    }
    // cur_bp == 0 is possible for e.g. an indel at the beginning of chr17.
    const char* prev_ref = &(seqbuf[S_CAST(int32_t, cur_bp) - 1]);
    uint32_t rtrim = 0;
    uint32_t lshift = 0;
    uint32_t ltrim = 0;
    const char* shifted_ref = nullptr;
    do {
      ++rtrim;
      uint32_t cur_alen = alen_buf[first_aidx];
      if (rtrim >= cur_alen) {
        right_match = ctou32(prev_ref[S_CAST(int32_t, cur_alen - rtrim)]);
      } else {
        right_match = ctou32(cur_alleles[first_aidx][cur_alen - 1 - rtrim]);
      }
      for (uint32_t aidx = first_aidx + 1; aidx != allele_ct; ++aidx) {
        if (aidx == missing_aidx) {
          continue;
        }
        cur_alen = alen_buf[aidx];
        uint32_t cur_code;
        if (rtrim >= cur_alen) {
          cur_code = ctou32(prev_ref[S_CAST(int32_t, cur_alen - rtrim)]);
        } else {
          cur_code = ctou32(cur_alleles[aidx][cur_alen - 1 - rtrim]);
        }
        if (right_match != cur_code) {
          goto VNormalizeContig_rtrim_done;
        }
      }
    } while (rtrim != rtrim_stop);
  VNormalizeContig_rtrim_done:
    if (rtrim >= min_alen) {
      lshift = rtrim + 1 - min_alen;
      cur_bp -= lshift;
      shifted_ref = &(seqbuf[cur_bp]);
      // min_alen = 1;
      variant_bps[variant_uidx] = cur_bp;
    } else {
      min_alen -= rtrim;
      if ((left_match != UINT32_MAX) && (min_alen > 1)) {
        const uint32_t min_alen_m1 = min_alen - 1;
        ltrim = 1;
        for (; ltrim != min_alen_m1; ++ltrim) {
          left_match = ctou32(cur_alleles[first_aidx][ltrim]);
          for (uint32_t aidx = first_aidx + 1; aidx != allele_ct; ++aidx) {
            if (aidx == missing_aidx) {
              continue;
            }
            const uint32_t cur_code = ctou32(cur_alleles[aidx][ltrim]);
            if (left_match != cur_code) {
              goto VNormalizeContig_ltrim2_done;
            }
          }
        }
      }
    VNormalizeContig_ltrim2_done:
      cur_bp += ltrim;
      variant_bps[variant_uidx] = cur_bp;
    }
    for (uint32_t aidx = first_aidx; aidx != allele_ct; ++aidx) {
      if (aidx == missing_aidx) {
        continue;
      }
      const uint32_t orig_alen = alen_buf[aidx];
      if (orig_alen <= rtrim) {
        // Must be single-character from .fa.
        cur_alleles[aidx] = &(g_one_char_strs[2 * ctou32(prev_ref[S_CAST(int32_t, orig_alen - rtrim)])]);
      } else {
        uint32_t new_slen = orig_alen + lshift - rtrim - ltrim;
        const char* cur_allele = &(cur_alleles[aidx][ltrim]);
        if (new_slen == 1) {
          cur_alleles[aidx] = &(g_one_char_strs[2 * ctou32(cur_allele[0])]);
        } else {
          if (S_CAST(uintptr_t, alloc_end - alloc_base) <= new_slen) {
            return kPglRetNomem;
          }
          alloc_end -= new_slen + 1;
          char* new_allele = R_CAST(char*, alloc_end);
          char* new_allele_iter = new_allele;
          if (lshift) {
            if (lshift < new_slen) {
              new_allele_iter = memcpya(new_allele_iter, shifted_ref, lshift);
              new_slen -= lshift;
            } else {
              new_allele_iter = memcpya(new_allele_iter, shifted_ref, new_slen);
              new_slen = 0;
            }
          }
          new_allele_iter = memcpya(new_allele_iter, cur_allele, new_slen);
          *new_allele_iter = '\0';
          cur_alleles[aidx] = new_allele;
        }
      }
    }
  }
  *alloc_endp = alloc_end;
  if (is_unsorted) {
    *vpos_sortstatusp |= kfUnsortedVarBp;
  }
  *nchanged_ct_ptr = nchanged_ct;
  *nlist_write_iterp = nlist_write_iter;
  return kPglRetSuccess;
}

PglErr ProcessFa(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const ChrInfo* cip, const char* fname, uint32_t max_allele_ct, uint32_t max_allele_slen, FaFlags flags, uint32_t max_thread_ct, UnsortedVar* vpos_sortstatusp, uint32_t* variant_bps, const char** allele_storage, STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uintptr_t* nonref_flags, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  FILE* nlist_file = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream fa_txs;
  PreinitTextStream(&fa_txs);
  {
    const uint32_t chr_ct = cip->chr_ct;
    char* chr_name_buf;
    uintptr_t* chr_already_seen;
    if (unlikely(
            bigstack_calloc_w(BitCtToWordCt(chr_ct), &chr_already_seen) ||
            bigstack_alloc_c(kMaxIdBlen, &chr_name_buf))) {
      goto ProcessFa_ret_NOMEM;
    }
    uint32_t* alen_buf = nullptr;
    char* nlist_write_iter = nullptr;
    char* nlist_flush = nullptr;
    if (flags & kfFaNormalize) {
      if (unlikely(bigstack_alloc_u32(max_allele_ct, &alen_buf))) {
        goto ProcessFa_ret_NOMEM;
      }
      if (flags & kfFaNormalizeList) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".normalized");
        if (unlikely(fopen_checked(outname, FOPEN_WB, &nlist_file))) {
          goto ProcessFa_ret_OPEN_FAIL;
        }
        nlist_write_iter = g_textbuf;
        nlist_flush = &(nlist_write_iter[kMaxMediumLine]);
      }
    }

    // To simplify indel/complex-variant handling, we load an entire contig at
    // a time.  Determine an upper bound for the size of this buffer.
    const uintptr_t* chr_mask = cip->chr_mask;
    uint32_t seqbuf_size = 0;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
      if (!IsSet(chr_mask, cip->chr_file_order[chr_fo_idx])) {
        continue;
      }
      const int32_t chr_vidx_start_m1 = cip->chr_fo_vidx_start[chr_fo_idx] - 1;
      const int32_t chr_vidx_last = FindLast1BitBeforeBounded(variant_include, cip->chr_fo_vidx_start[chr_fo_idx + 1], chr_vidx_start_m1);
      if (chr_vidx_last != chr_vidx_start_m1) {
        const uint32_t cur_bp = variant_bps[S_CAST(uint32_t, chr_vidx_last)];
        if (cur_bp > seqbuf_size) {
          seqbuf_size = cur_bp;
        }
      }
    }
    seqbuf_size += max_allele_slen + 1;
    char* seqbuf;
    if (unlikely(bigstack_alloc_c(seqbuf_size, &seqbuf))) {
      goto ProcessFa_ret_NOMEM;
    }
    // May as well handle insertion before first contig base, and deletion of
    // first base, correctly.
    seqbuf[0] = 'N';

    reterr = SizeAndInitTextStream(fname, bigstack_left(), MAXV(max_thread_ct - 1, 1), &fa_txs);
    if (unlikely(reterr)) {
      goto ProcessFa_ret_TSTREAM_FAIL;
    }

    unsigned char* tmp_alloc_end = g_bigstack_end;
    char* seqbuf_end = nullptr;
    char* seq_iter = nullptr;
    // possible but low-priority todo: exploit .fai when it's present
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t cur_vidx_last = 0;
    uint32_t skip_chr = 1;
    uint32_t is_first_noncomment_line = 1;
    uint32_t ref_changed_ct = 0;
    uint32_t ref_validated_ct = 0;
    uint32_t ref_downgraded_ct = 0;
    uint32_t nchanged_ct = 0;
    while (1) {
      ++line_idx;
      char* line_iter;
      reterr = TextNextLine(&fa_txs, &line_iter);
      if (reterr) {
        if (likely(reterr == kPglRetEof)) {
          reterr = kPglRetSuccess;
          break;
        }
        goto ProcessFa_ret_TSTREAM_FAIL;
      }
      unsigned char ucc = line_iter[0];
      if (ucc < 'A') {
        // > = ascii 62
        // ; = ascii 59
        if (ucc == ';') {
          continue;
        }
        is_first_noncomment_line = 0;
        if (unlikely(ucc != '>')) {
          snprintf(g_logbuf, kLogbufSize, "Error: Unexpected character at beginning of line %" PRIuPTR " of --fa file.\n", line_idx);
          goto ProcessFa_ret_MALFORMED_INPUT_WW;
        }
        if (chr_fo_idx != UINT32_MAX) {
          *seq_iter = '\0';
          const uint32_t bp_end = seq_iter - seqbuf;
          if (flags & kfFaRefFrom) {
            reterr = RefFromFaContig(variant_include, variant_bps, allele_idx_offsets, allele_storage, cip, seqbuf, flags & kfFaRefFromForce, chr_fo_idx, cur_vidx_last, bp_end, refalt1_select, nonref_flags, &ref_changed_ct, &ref_validated_ct, &ref_downgraded_ct);
            if (unlikely(reterr)) {
              goto ProcessFa_ret_1;
            }
          }
          if (flags & kfFaNormalize) {
            reterr = VNormalizeContig(variant_include, variant_ids, allele_idx_offsets, cip, seqbuf, chr_fo_idx, cur_vidx_last, bp_end, &tmp_alloc_end, vpos_sortstatusp, variant_bps, allele_storage, &nchanged_ct, nlist_flush, nlist_file, &nlist_write_iter, alen_buf);
            if (unlikely(reterr)) {
              goto ProcessFa_ret_1;
            }
          }
        }
        char* chr_name_start = &(line_iter[1]);
        if (unlikely(IsSpaceOrEoln(*chr_name_start))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid contig description on line %" PRIuPTR " of --fa file.%s\n", line_idx, (*chr_name_start == ' ')? " (Spaces are not permitted between the leading '>' and the contig name.)" : "");
          goto ProcessFa_ret_MALFORMED_INPUT_WW;
        }
        char* chr_name_end = CurTokenEnd(chr_name_start);
        line_iter = AdvToDelim(chr_name_end, '\n');
        *chr_name_end = '\0';
        const uint32_t chr_name_slen = chr_name_end - chr_name_start;
        chr_fo_idx = UINT32_MAX;
        const uint32_t chr_idx = GetChrCode(chr_name_start, cip, chr_name_slen);
        if ((!IsI32Neg(chr_idx)) && IsSet(cip->chr_mask, chr_idx)) {
          chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
          if (unlikely(IsSet(chr_already_seen, chr_fo_idx))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Duplicate contig name '%s' in --fa file.\n", chr_name_start);
            goto ProcessFa_ret_MALFORMED_INPUT_WW;
          }
          SetBit(chr_fo_idx, chr_already_seen);
          const int32_t chr_vidx_start_m1 = cip->chr_fo_vidx_start[chr_fo_idx] - 1;
          const int32_t chr_vidx_last = FindLast1BitBeforeBounded(variant_include, cip->chr_fo_vidx_start[chr_fo_idx + 1], chr_vidx_start_m1);
          if (chr_vidx_last == chr_vidx_start_m1) {
            chr_fo_idx = UINT32_MAX;
          } else {
            cur_vidx_last = chr_vidx_last;
            seqbuf_end = &(seqbuf[variant_bps[cur_vidx_last] + max_allele_slen]);
            seq_iter = &(seqbuf[1]);
          }
        }
        skip_chr = (chr_fo_idx == UINT32_MAX);
        *chr_name_end = '\n';  // bugfix (16 Apr 2018)
        continue;
      }
      const char* line_start = line_iter;
      line_iter = CurTokenEnd(line_iter);
      const char* seqline_end = line_iter;
      if (skip_chr) {
        if (is_first_noncomment_line) {
          logerrputs("Error: --fa file is not a valid FASTA.\n");
          goto ProcessFa_ret_MALFORMED_INPUT;
        }
        continue;
      }
      ucc = *seqline_end;
      if (unlikely((ucc == ' ') || (ucc == '\t'))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --fa file is malformed.\n", line_idx);
        goto ProcessFa_ret_MALFORMED_INPUT_2;
      }
      uint32_t cur_seq_slen = seqline_end - line_start;
      const uint32_t seq_rem = seqbuf_end - seq_iter;
      if (seq_rem <= cur_seq_slen) {
        cur_seq_slen = seq_rem;
        skip_chr = 1;
      }
      const char* gap_start = S_CAST(const char*, memchr(line_start, '-', cur_seq_slen));
      if (gap_start) {
        logerrprintfww("Warning: Indeterminate-length gap present on line %" PRIuPTR " of --fa file. Ignoring remainder of contig.\n", line_idx);
        cur_seq_slen = gap_start - line_start;
        skip_chr = 1;
      }
      seq_iter = memcpya(seq_iter, line_start, cur_seq_slen);
      // seq_iter = memcpya_toupper(seq_iter, loadbuf, cur_seq_slen);
    }
    if (flags & kfFaRefFrom) {
      const uint32_t ref_from_fa_force = flags & kfFaRefFromForce;
      if (chr_fo_idx != UINT32_MAX) {
        *seq_iter = '\0';
        const uint32_t bp_end = seq_iter - seqbuf;
        reterr = RefFromFaContig(variant_include, variant_bps, allele_idx_offsets, allele_storage, cip, seqbuf, ref_from_fa_force, chr_fo_idx, cur_vidx_last, bp_end, refalt1_select, nonref_flags, &ref_changed_ct, &ref_validated_ct, &ref_downgraded_ct);
        if (unlikely(reterr)) {
          goto ProcessFa_ret_1;
        }
      }
      logprintf("--ref-from-fa%s: %u variant%s changed, %u validated.\n", ref_from_fa_force? " force" : "", ref_changed_ct, (ref_changed_ct == 1)? "" : "s", ref_validated_ct);
      if (ref_downgraded_ct) {
        logerrprintfww("Warning: %u reference allele%s downgraded from 'known' to 'provisional'.\n", ref_downgraded_ct, (ref_downgraded_ct == 1)? "" : "s");
      }
    }
    if (flags & kfFaNormalize) {
      if (chr_fo_idx != UINT32_MAX) {
        // slightly redundant
        *seq_iter = '\0';
        const uint32_t bp_end = seq_iter - seqbuf;
        reterr = VNormalizeContig(variant_include, variant_ids, allele_idx_offsets, cip, seqbuf, chr_fo_idx, cur_vidx_last, bp_end, &tmp_alloc_end, vpos_sortstatusp, variant_bps, allele_storage, &nchanged_ct, nlist_flush, nlist_file, &nlist_write_iter, alen_buf);
        if (unlikely(reterr)) {
          goto ProcessFa_ret_1;
        }
      }
      BigstackEndSet(tmp_alloc_end);
      logprintf("--normalize: %u variant%s changed.\n", nchanged_ct, (nchanged_ct == 1)? "" : "s");
      if (nlist_file) {
        if (unlikely(fclose_flush_null(nlist_flush, nlist_write_iter, &nlist_file))) {
          goto ProcessFa_ret_WRITE_FAIL;
        }
        logprintfww("Affected variant ID%s written to %s .\n", (nchanged_ct == 1)? "" : "s", outname);
      }
      if ((*vpos_sortstatusp) & kfUnsortedVarBp) {
        logerrprintf("Warning: Base-pair positions are now unsorted!\n");
      }
    }
  }
  while (0) {
  ProcessFa_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ProcessFa_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ProcessFa_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--ref-from-fa file", &fa_txs);
    break;
  ProcessFa_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ProcessFa_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
  ProcessFa_ret_MALFORMED_INPUT_2:
    logerrputsb();
  ProcessFa_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 ProcessFa_ret_1:
  fclose_cond(nlist_file);
  CleanupTextStream2("--ref-from-fa file", &fa_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

// Eventually want a standalone function which uses a .fa to fill in ref
// alleles for .ped and similar all-SNP datasets where the reference allele may
// have disappeared, and some of the resulting variants will become
// multiallelic.  But wait till --pmerge is implemented, since this is mostly a
// subset of what --pmerge needs to do.

#ifdef __cplusplus
}  // namespace plink2
#endif
