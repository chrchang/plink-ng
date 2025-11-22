// This file is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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

#include "plink2_ld.h"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <unistd.h>  // unlink()

#include "include/SFMT.h"
#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_stats.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "plink2_compress_stream.h"
#include "plink2_decompress.h"
#include "plink2_filter.h"
#include "plink2_set.h"

#ifdef __cplusplus
#  include <functional>  // std::greater
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

void InitLd(LdInfo* ldip) {
  ldip->prune_flags = kfLdPrune0;
  ldip->prune_window_size = 0;
  ldip->prune_window_incr = 0;
  ldip->prune_last_param = 0.0;
  ldip->ld_console_flags = kfLdConsole0;
  ldip->ld_console_varids[0] = nullptr;
  ldip->ld_console_varids[1] = nullptr;
}

void CleanupLd(LdInfo* ldip) {
  free_cond(ldip->ld_console_varids[0]);
  free_cond(ldip->ld_console_varids[1]);
}

void InitClump(ClumpInfo* clump_ip) {
  clump_ip->fnames_flattened = nullptr;
  clump_ip->range_fname = nullptr;
  clump_ip->test_name = nullptr;
  clump_ip->id_field = nullptr;
  clump_ip->a1_field = nullptr;
  clump_ip->test_field = nullptr;
  clump_ip->p_field = nullptr;
  clump_ip->ln_bin_boundaries = nullptr;
  clump_ip->ln_p1 = kLn10 * -4.0 * (1.0 - kSmallEpsilon);
  clump_ip->ln_p2 = kLn10 * -2.0 * (1.0 - kSmallEpsilon);
  clump_ip->r2 = 0.5 * (1.0 + kSmallEpsilon);
  clump_ip->bin_bound_ct = 0;
  clump_ip->bp_radius = 249999;
  clump_ip->range_border = 0;
  clump_ip->flags = kfClump0;
}

void CleanupClump(ClumpInfo* clump_ip) {
  free_cond(clump_ip->fnames_flattened);
  free_cond(clump_ip->range_fname);
  free_cond(clump_ip->test_name);
  free_cond(clump_ip->id_field);
  free_cond(clump_ip->a1_field);
  free_cond(clump_ip->test_field);
  free_cond(clump_ip->p_field);
  free_cond(clump_ip->ln_bin_boundaries);
}

void InitVcor(VcorInfo* vcip) {
  vcip->ld_snp_list_fname = nullptr;
  InitRangeList(&(vcip->ld_snp_range_list));
  // We want to error out during command-line parsing when the user specifies
  // matrix output with one of these filters, or inter-chr output with one of
  // the radius filters.
  // Unfortunately, --ld-window... flags are parsed before --r[2]-[un]phased,
  // so directly initializing bp_radius and min_r2 to default values doesn't
  // work out cleanly.  Instead, we initialize them to out-of-range values, and
  // fill in defaults during --r[2]-[un]phased parsing.
  vcip->var_ct_radius = 0x7fffffff; // avoids overflows UINT32_MAX would have
  vcip->bp_radius = UINT32_MAX;
  vcip->cm_radius = -1.0;
  vcip->min_r2 = 2.0;
  vcip->flags = kfVcor0;
}

void CleanupVcor(VcorInfo* vcip) {
  free_cond(vcip->ld_snp_list_fname);
  CleanupRangeList(&(vcip->ld_snp_range_list));
}

void StripUnplacedNoCount(const ChrInfo* cip, uintptr_t* variant_include) {
  if (IsSet(cip->chr_mask, 0)) {
    const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[0];
    const uint32_t start_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
    ClearBitsNz(start_uidx, cip->chr_fo_vidx_start[chr_fo_idx + 1], variant_include);
  }
  if (cip->zero_extra_chrs) {
    const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
    for (uint32_t chr_idx = cip->max_code + 1; chr_idx != chr_code_end; ++chr_idx) {
      const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
      const uint32_t start_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
      ClearBitsNz(start_uidx, cip->chr_fo_vidx_start[chr_fo_idx + 1], variant_include);
    }
  }
}

// Move to plink2_common if any users outside plink2_ld.
BoolErr StripUnplaced(const uintptr_t* orig_variant_include, const ChrInfo* cip, uint32_t raw_variant_ct, uint32_t* skipped_variant_ctp, uintptr_t** new_variant_includep) {
  uint32_t skipped_variant_ct = 0;
  if (IsSet(cip->chr_mask, 0)) {
    skipped_variant_ct = CountChrVariantsUnsafe(orig_variant_include, cip, 0);
  }
  const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
  if (cip->zero_extra_chrs) {
    for (uint32_t chr_idx = cip->max_code + 1; chr_idx != chr_code_end; ++chr_idx) {
      if (IsSet(cip->chr_mask, chr_idx)) {
        skipped_variant_ct += CountChrVariantsUnsafe(orig_variant_include, cip, cip->chr_idx_to_foidx[chr_idx]);
      }
    }
  }
  *skipped_variant_ctp = skipped_variant_ct;
  if (!skipped_variant_ct) {
    *new_variant_includep = nullptr;
    return 0;
  }
  const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
  if (unlikely(bigstack_alloc_w(raw_variant_ctl, new_variant_includep))) {
    return 1;
  }
  uintptr_t* new_variant_include = *new_variant_includep;
  memcpy(new_variant_include, orig_variant_include, raw_variant_ctl * sizeof(intptr_t));
  StripUnplacedNoCount(cip, new_variant_include);
  return 0;
}

static inline const uintptr_t* StripUnplacedK(const uintptr_t* orig_variant_include, const ChrInfo* cip, uint32_t raw_variant_ct, uint32_t* skipped_variant_ctp) {
  uintptr_t* new_variant_include;
  if (unlikely(StripUnplaced(orig_variant_include, cip, raw_variant_ct, skipped_variant_ctp, &new_variant_include))) {
    return nullptr;
  }
  return new_variant_include? new_variant_include : orig_variant_include;
}

// Returns number of skipped variants.
uint32_t StripUnplacedMut(const ChrInfo* cip, uintptr_t* variant_include) {
  uint32_t skipped_variant_ct = 0;
  if (IsSet(cip->chr_mask, 0)) {
    skipped_variant_ct = CountChrVariantsUnsafe(variant_include, cip, 0);
  }
  const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
  if (cip->zero_extra_chrs) {
    for (uint32_t chr_idx = cip->max_code + 1; chr_idx != chr_code_end; ++chr_idx) {
      if (IsSet(cip->chr_mask, chr_idx)) {
        // bugfix (15 Jan 2025): passed wrong last parameter
        skipped_variant_ct += CountChrVariantsUnsafe(variant_include, cip, chr_idx);
      }
    }
  }
  if (!skipped_variant_ct) {
    return 0;
  }
  StripUnplacedNoCount(cip, variant_include);
  return skipped_variant_ct;
}


#ifdef USE_AVX2
// todo: see if either approach in avx_jaccard_index.c in
// github.com/CountOnes/hamming_weight helps here.

static inline int32_t DotprodAvx2(const VecW* __restrict hom1_iter, const VecW* __restrict ref2het1_iter, const VecW* __restrict hom2_iter, const VecW* __restrict ref2het2_iter, uintptr_t vec_ct) {
  // popcount(hom1 & hom2) - 2 * popcount(hom1 & hom2 & (ref2het1 ^ ref2het2))
  // ct must be a multiple of 4.
  VecW cnt = vecw_setzero();
  VecW ones_both = vecw_setzero();
  VecW ones_neg = vecw_setzero();
  VecW twos_both = vecw_setzero();
  VecW twos_neg = vecw_setzero();
  for (uintptr_t vec_idx = 0; vec_idx < vec_ct; vec_idx += 4) {
    VecW count1_both = hom1_iter[vec_idx] & hom2_iter[vec_idx];
    VecW cur_xor = ref2het1_iter[vec_idx] ^ ref2het2_iter[vec_idx];
    VecW count1_neg = count1_both & cur_xor;

    VecW count2_both = hom1_iter[vec_idx + 1] & hom2_iter[vec_idx + 1];
    cur_xor = ref2het1_iter[vec_idx + 1] ^ ref2het2_iter[vec_idx + 1];
    VecW count2_neg = count2_both & cur_xor;
    const VecW twos_both_a = Csa256(count1_both, count2_both, &ones_both);
    const VecW twos_neg_a = Csa256(count1_neg, count2_neg, &ones_neg);

    count1_both = hom1_iter[vec_idx + 2] & hom2_iter[vec_idx + 2];
    cur_xor = ref2het1_iter[vec_idx + 2] ^ ref2het2_iter[vec_idx + 2];
    count1_neg = count1_both & cur_xor;

    count2_both = hom1_iter[vec_idx + 3] & hom2_iter[vec_idx + 3];
    cur_xor = ref2het1_iter[vec_idx + 3] ^ ref2het2_iter[vec_idx + 3];
    count2_neg = count2_both & cur_xor;
    const VecW twos_both_b = Csa256(count1_both, count2_both, &ones_both);
    const VecW twos_neg_b = Csa256(count1_neg, count2_neg, &ones_neg);
    const VecW fours_both = Csa256(twos_both_a, twos_both_b, &twos_both);
    const VecW fours_neg = Csa256(twos_neg_a, twos_neg_b, &twos_neg);
    // tried continuing to eights, not worth it
    // deliberate unsigned-int64 overflow here
    cnt = cnt + PopcountVecAvx2(fours_both) - vecw_slli(PopcountVecAvx2(fours_neg), 1);
  }
  cnt = cnt - PopcountVecAvx2(twos_neg);
  cnt = vecw_slli(cnt, 2);
  const VecW twos_sum = PopcountVecAvx2(twos_both) - PopcountVecAvx2(ones_neg);
  cnt = cnt + vecw_slli(twos_sum, 1);
  cnt = cnt + PopcountVecAvx2(ones_both);
  return HsumW(cnt);
}

int32_t DotprodWords(const uintptr_t* __restrict hom1, const uintptr_t* __restrict ref2het1, const uintptr_t* __restrict hom2, const uintptr_t* __restrict ref2het2, uintptr_t word_ct) {
  int32_t tot_both = 0;
  uint32_t widx = 0;
  if (word_ct >= 16) { // this already pays off with a single block
    const uintptr_t block_ct = word_ct / (kWordsPerVec * 4);
    tot_both = DotprodAvx2(R_CAST(const VecW*, hom1), R_CAST(const VecW*, ref2het1), R_CAST(const VecW*, hom2), R_CAST(const VecW*, ref2het2), block_ct * 4);
    widx = block_ct * (4 * kWordsPerVec);
  }
  uint32_t tot_neg = 0;
  for (; widx != word_ct; ++widx) {
    const uintptr_t hom_word = hom1[widx] & hom2[widx];
    const uintptr_t xor_word = ref2het1[widx] ^ ref2het2[widx];
    tot_both += PopcountWord(hom_word);
    tot_neg += PopcountWord(hom_word & xor_word);
  }
  return tot_both - 2 * tot_neg;
}

static inline void SumSsqAvx2(const VecW* __restrict hom1_iter, const VecW* __restrict ref2het1_iter, const VecW* __restrict hom2_iter, const VecW* __restrict ref2het2_iter, uintptr_t vec_ct, uint32_t* __restrict ssq2_ptr, uint32_t* __restrict plus2_ptr) {
  // popcounts (nm1 & hom2) and (nm1 & hom2 & ref2het2).  ct is multiple of 8.
  VecW cnt_ssq = vecw_setzero();
  VecW cnt_plus = vecw_setzero();
  VecW ones_ssq = vecw_setzero();
  VecW ones_plus = vecw_setzero();
  VecW twos_ssq = vecw_setzero();
  VecW twos_plus = vecw_setzero();
  VecW fours_ssq = vecw_setzero();
  VecW fours_plus = vecw_setzero();
  for (uintptr_t vec_idx = 0; vec_idx < vec_ct; vec_idx += 8) {
    VecW count1_ssq = (hom1_iter[vec_idx] | ref2het1_iter[vec_idx]) & hom2_iter[vec_idx];
    VecW count1_plus = count1_ssq & ref2het2_iter[vec_idx];

    VecW count2_ssq = (hom1_iter[vec_idx + 1] | ref2het1_iter[vec_idx + 1]) & hom2_iter[vec_idx + 1];
    VecW count2_plus = count2_ssq & ref2het2_iter[vec_idx + 1];
    VecW twos_ssq_a = Csa256(count1_ssq, count2_ssq, &ones_ssq);
    VecW twos_plus_a = Csa256(count1_plus, count2_plus, &ones_plus);

    count1_ssq = (hom1_iter[vec_idx + 2] | ref2het1_iter[vec_idx + 2]) & hom2_iter[vec_idx + 2];
    count1_plus = count1_ssq & ref2het2_iter[vec_idx + 2];

    count2_ssq = (hom1_iter[vec_idx + 3] | ref2het1_iter[vec_idx + 3]) & hom2_iter[vec_idx + 3];
    count2_plus = count2_ssq & ref2het2_iter[vec_idx + 3];
    VecW twos_ssq_b = Csa256(count1_ssq, count2_ssq, &ones_ssq);
    VecW twos_plus_b = Csa256(count1_plus, count2_plus, &ones_plus);
    const VecW fours_ssq_a = Csa256(twos_ssq_a, twos_ssq_b, &twos_ssq);
    const VecW fours_plus_a = Csa256(twos_plus_a, twos_plus_b, &twos_plus);

    count1_ssq = (hom1_iter[vec_idx + 4] | ref2het1_iter[vec_idx + 4]) & hom2_iter[vec_idx + 4];
    count1_plus = count1_ssq & ref2het2_iter[vec_idx + 4];

    count2_ssq = (hom1_iter[vec_idx + 5] | ref2het1_iter[vec_idx + 5]) & hom2_iter[vec_idx + 5];
    count2_plus = count2_ssq & ref2het2_iter[vec_idx + 5];
    twos_ssq_a = Csa256(count1_ssq, count2_ssq, &ones_ssq);
    twos_plus_a = Csa256(count1_plus, count2_plus, &ones_plus);

    count1_ssq = (hom1_iter[vec_idx + 6] | ref2het1_iter[vec_idx + 6]) & hom2_iter[vec_idx + 6];
    count1_plus = count1_ssq & ref2het2_iter[vec_idx + 6];

    count2_ssq = (hom1_iter[vec_idx + 7] | ref2het1_iter[vec_idx + 7]) & hom2_iter[vec_idx + 7];
    count2_plus = count2_ssq & ref2het2_iter[vec_idx + 7];
    twos_ssq_b = Csa256(count1_ssq, count2_ssq, &ones_ssq);
    twos_plus_b = Csa256(count1_plus, count2_plus, &ones_plus);
    const VecW fours_ssq_b = Csa256(twos_ssq_a, twos_ssq_b, &twos_ssq);
    const VecW fours_plus_b = Csa256(twos_plus_a, twos_plus_b, &twos_plus);
    const VecW eights_ssq = Csa256(fours_ssq_a, fours_ssq_b, &fours_ssq);
    const VecW eights_plus = Csa256(fours_plus_a, fours_plus_b, &fours_plus);
    // negligible benefit from going to sixteens here
    cnt_ssq = cnt_ssq + PopcountVecAvx2(eights_ssq);
    cnt_plus = cnt_plus + PopcountVecAvx2(eights_plus);
  }
  cnt_ssq = vecw_slli(cnt_ssq, 3);
  cnt_plus = vecw_slli(cnt_plus, 3);
  cnt_ssq = cnt_ssq + vecw_slli(PopcountVecAvx2(fours_ssq), 2);
  cnt_plus = cnt_plus + vecw_slli(PopcountVecAvx2(fours_plus), 2);
  cnt_ssq = cnt_ssq + vecw_slli(PopcountVecAvx2(twos_ssq), 1);
  cnt_plus = cnt_plus + vecw_slli(PopcountVecAvx2(twos_plus), 1);
  cnt_ssq = cnt_ssq + PopcountVecAvx2(ones_ssq);
  cnt_plus = cnt_plus + PopcountVecAvx2(ones_plus);
  *ssq2_ptr = HsumW(cnt_ssq);
  *plus2_ptr = HsumW(cnt_plus);
}

void SumSsqWords(const uintptr_t* hom1, const uintptr_t* ref2het1, const uintptr_t* hom2, const uintptr_t* ref2het2, uint32_t word_ct, int32_t* sum2_ptr, uint32_t* ssq2_ptr) {
  uint32_t ssq2 = 0;
  uint32_t plus2 = 0;
  uint32_t widx = 0;
  // this has high constant overhead at the end; may want to require more than
  // one block?
  if (word_ct >= 32) {
    const uintptr_t block_ct = word_ct / (8 * kWordsPerVec);
    SumSsqAvx2(R_CAST(const VecW*, hom1), R_CAST(const VecW*, ref2het1), R_CAST(const VecW*, hom2), R_CAST(const VecW*, ref2het2), block_ct * 8, &ssq2, &plus2);
    widx = block_ct * (8 * kWordsPerVec);
  }
  for (; widx != word_ct; ++widx) {
    const uintptr_t ssq2_word = (hom1[widx] | ref2het1[widx]) & hom2[widx];
    ssq2 += PopcountWord(ssq2_word);
    plus2 += PopcountWord(ssq2_word & ref2het2[widx]);
  }
  *sum2_ptr = S_CAST(int32_t, 2 * plus2 - ssq2);  // deliberate overflow
  *ssq2_ptr = ssq2;
}
#else  // !USE_AVX2
static inline int32_t DotprodVecsNm(const VecW* __restrict hom1_iter, const VecW* __restrict ref2het1_iter, const VecW* __restrict hom2_iter, const VecW* __restrict ref2het2_iter, uintptr_t vec_ct) {
  // popcount(hom1 & hom2) - 2 * popcount(hom1 & hom2 & (ref2het1 ^ ref2het2))
  // ct must be a multiple of 3.
  assert(!(vec_ct % 3));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  VecW acc_both = vecw_setzero();
  VecW acc_neg = vecw_setzero();
  uintptr_t cur_incr = 30;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 30) {
      if (!vec_ct) {
        return HsumW(acc_both) - 2 * HsumW(acc_neg);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_both = vecw_setzero();
    VecW inner_acc_neg = vecw_setzero();
    const VecW* hom1_stop = &(hom1_iter[cur_incr]);
    do {
      VecW count1_both = (*hom1_iter++) & (*hom2_iter++);
      VecW cur_xor = (*ref2het1_iter++) ^ (*ref2het2_iter++);
      VecW count1_neg = count1_both & cur_xor;

      VecW count2_both = (*hom1_iter++) & (*hom2_iter++);
      cur_xor = (*ref2het1_iter++) ^ (*ref2het2_iter++);
      VecW count2_neg = count2_both & cur_xor;

      VecW cur_hom = (*hom1_iter++) & (*hom2_iter++);
      cur_xor = (*ref2het1_iter++) ^ (*ref2het2_iter++);
      VecW half1_neg = cur_hom & cur_xor;
      const VecW half2_both = vecw_srli(cur_hom, 1) & m1;
      const VecW half2_neg = vecw_srli(half1_neg, 1) & m1;
      const VecW half1_both = cur_hom & m1;
      half1_neg = half1_neg & m1;
      count1_both = count1_both - (vecw_srli(count1_both, 1) & m1);
      count1_neg = count1_neg - (vecw_srli(count1_neg, 1) & m1);
      count2_both = count2_both - (vecw_srli(count2_both, 1) & m1);
      count2_neg = count2_neg - (vecw_srli(count2_neg, 1) & m1);
      count1_both = count1_both + half1_both;
      count1_neg = count1_neg + half1_neg;
      count2_both = count2_both + half2_both;
      count2_neg = count2_neg + half2_neg;
      count1_both = (count1_both & m2) + (vecw_srli(count1_both, 2) & m2);
      count1_neg = (count1_neg & m2) + (vecw_srli(count1_neg, 2) & m2);
      count1_both = count1_both + (count2_both & m2) + (vecw_srli(count2_both, 2) & m2);
      count1_neg = count1_neg + (count2_neg & m2) + (vecw_srli(count2_neg, 2) & m2);
      inner_acc_both = inner_acc_both + (count1_both & m4) + (vecw_srli(count1_both, 4) & m4);
      inner_acc_neg = inner_acc_neg + (count1_neg & m4) + (vecw_srli(count1_neg, 4) & m4);
    } while (hom1_iter < hom1_stop);
    const VecW m0 = vecw_setzero();
    acc_both = acc_both + vecw_bytesum(inner_acc_both, m0);
    acc_neg = acc_neg + vecw_bytesum(inner_acc_neg, m0);
  }
}

int32_t DotprodWords(const uintptr_t* __restrict hom1, const uintptr_t* __restrict ref2het1, const uintptr_t* __restrict hom2, const uintptr_t* __restrict ref2het2, uintptr_t word_ct) {
  int32_t tot_both = 0;
  uint32_t widx = 0;
  if (word_ct >= kWordsPerVec * 3) {
    const uintptr_t block_ct = word_ct / (kWordsPerVec * 3);
    tot_both = DotprodVecsNm(R_CAST(const VecW*, hom1), R_CAST(const VecW*, ref2het1), R_CAST(const VecW*, hom2), R_CAST(const VecW*, ref2het2), block_ct * 3);
    widx = block_ct * (3 * kWordsPerVec);
  }
  uint32_t tot_neg = 0;
  for (; widx != word_ct; ++widx) {
    const uintptr_t hom_word = hom1[widx] & hom2[widx];
    const uintptr_t xor_word = ref2het1[widx] ^ ref2het2[widx];
    tot_both += PopcountWord(hom_word);
    tot_neg += PopcountWord(hom_word & xor_word);
  }
  return tot_both - 2 * tot_neg;
}

#ifndef USE_SSE42
static inline void SumSsqVecs(const VecW* __restrict hom1_iter, const VecW* __restrict ref2het1_iter, const VecW* __restrict hom2_iter, const VecW* __restrict ref2het2_iter, uintptr_t vec_ct, uint32_t* __restrict ssq2_ptr, uint32_t* __restrict plus2_ptr) {
  // popcounts (nm1 & hom2) and (nm1 & hom2 & ref2het2).  ct is multiple of 3.
  assert(!(vec_ct % 3));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  VecW acc_ssq2 = vecw_setzero();
  VecW acc_plus2 = vecw_setzero();
  uint32_t cur_incr = 30;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 30) {
      if (!vec_ct) {
        *ssq2_ptr = HsumW(acc_ssq2);
        *plus2_ptr = HsumW(acc_plus2);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_ssq = vecw_setzero();
    VecW inner_acc_plus = vecw_setzero();
    for (uint32_t vec_idx = 0; vec_idx < cur_incr; vec_idx += 3) {
      VecW count1_ssq = (hom1_iter[vec_idx] | ref2het1_iter[vec_idx]) & hom2_iter[vec_idx];
      VecW count1_plus = count1_ssq & ref2het2_iter[vec_idx];

      VecW count2_ssq = (hom1_iter[vec_idx + 1] | ref2het1_iter[vec_idx + 1]) & hom2_iter[vec_idx + 1];
      VecW count2_plus = count2_ssq & ref2het2_iter[vec_idx + 1];

      VecW half1_ssq = (hom1_iter[vec_idx + 2] | ref2het1_iter[vec_idx + 2]) & hom2_iter[vec_idx + 2];
      VecW half1_plus = half1_ssq & ref2het2_iter[vec_idx + 2];
      VecW half2_ssq = vecw_srli(half1_ssq, 1) & m1;
      VecW half2_plus = vecw_srli(half1_plus, 1) & m1;
      half1_ssq = half1_ssq & m1;
      half1_plus = half1_plus & m1;
      count1_ssq = count1_ssq - (vecw_srli(count1_ssq, 1) & m1);
      count1_plus = count1_plus - (vecw_srli(count1_plus, 1) & m1);
      count2_ssq = count2_ssq - (vecw_srli(count2_ssq, 1) & m1);
      count2_plus = count2_plus - (vecw_srli(count2_plus, 1) & m1);
      count1_ssq = count1_ssq + half1_ssq;
      count1_plus = count1_plus + half1_plus;
      count2_ssq = count2_ssq + half2_ssq;
      count2_plus = count2_plus + half2_plus;
      count1_ssq = (count1_ssq & m2) + (vecw_srli(count1_ssq, 2) & m2);
      count1_plus = (count1_plus & m2) + (vecw_srli(count1_plus, 2) & m2);
      count1_ssq = count1_ssq + (count2_ssq & m2) + (vecw_srli(count2_ssq, 2) & m2);
      count1_plus = count1_plus + (count2_plus & m2) + (vecw_srli(count2_plus, 2) & m2);
      inner_acc_ssq = inner_acc_ssq + (count1_ssq & m4) + (vecw_srli(count1_ssq, 4) & m4);
      inner_acc_plus = inner_acc_plus + (count1_plus & m4) + (vecw_srli(count1_plus, 4) & m4);
    }
    hom1_iter = &(hom1_iter[cur_incr]);
    ref2het1_iter = &(ref2het1_iter[cur_incr]);
    hom2_iter = &(hom2_iter[cur_incr]);
    ref2het2_iter = &(ref2het2_iter[cur_incr]);
    const VecW m0 = vecw_setzero();
    acc_ssq2 = acc_ssq2 + vecw_bytesum(inner_acc_ssq, m0);
    acc_plus2 = acc_plus2 + vecw_bytesum(inner_acc_plus, m0);
  }
}
#endif

void SumSsqWords(const uintptr_t* hom1, const uintptr_t* ref2het1, const uintptr_t* hom2, const uintptr_t* ref2het2, uint32_t word_ct, int32_t* sum2_ptr, uint32_t* ssq2_ptr) {
  uint32_t ssq2 = 0;
  uint32_t plus2 = 0;
  uint32_t widx = 0;
#ifndef USE_SSE42
  if (word_ct >= kWordsPerVec * 3) {
    const uintptr_t block_ct = word_ct / (kWordsPerVec * 3);
    SumSsqVecs(R_CAST(const VecW*, hom1), R_CAST(const VecW*, ref2het1), R_CAST(const VecW*, hom2), R_CAST(const VecW*, ref2het2), block_ct * 3, &ssq2, &plus2);
    widx = block_ct * (3 * kWordsPerVec);
  }
#endif
  for (; widx != word_ct; ++widx) {
    const uintptr_t ssq2_word = (hom1[widx] | ref2het1[widx]) & hom2[widx];
    ssq2 += PopcountWord(ssq2_word);
    plus2 += PopcountWord(ssq2_word & ref2het2[widx]);
  }
  *sum2_ptr = S_CAST(int32_t, 2 * plus2 - ssq2);  // deliberate overflow
  *ssq2_ptr = ssq2;
}
#endif  // !USE_AVX2

#if defined(USE_AVX2) || !defined(USE_SSE42)
static inline void SumSsqNmVecs(const VecW* __restrict hom1_iter, const VecW* __restrict ref2het1_iter, const VecW* __restrict hom2_iter, const VecW* __restrict ref2het2_iter, uintptr_t vec_ct, uint32_t* __restrict nm_ptr, uint32_t* __restrict ssq2_ptr, uint32_t* __restrict plus2_ptr) {
  // vec_ct must be a multiple of 3.
  assert(!(vec_ct % 3));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  VecW acc_nm = vecw_setzero();
  VecW acc_ssq2 = vecw_setzero();
  VecW acc_plus2 = vecw_setzero();
  uint32_t cur_incr = 30;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 30) {
      if (!vec_ct) {
        *nm_ptr = HsumW(acc_nm);
        *ssq2_ptr = HsumW(acc_ssq2);
        *plus2_ptr = HsumW(acc_plus2);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_nm = vecw_setzero();
    VecW inner_acc_ssq = vecw_setzero();
    VecW inner_acc_plus = vecw_setzero();
    for (uint32_t vec_idx = 0; vec_idx < cur_incr; vec_idx += 3) {
      VecW nm1 = hom1_iter[vec_idx] | ref2het1_iter[vec_idx];
      VecW hom2 = hom2_iter[vec_idx];
      VecW ref2het2 = ref2het2_iter[vec_idx];
      VecW count1_ssq = nm1 & hom2;
      VecW count1_nm = nm1 & (hom2 | ref2het2);
      VecW count1_plus = count1_ssq & ref2het2;

      nm1 = hom1_iter[vec_idx + 1] | ref2het1_iter[vec_idx + 1];
      hom2 = hom2_iter[vec_idx + 1];
      ref2het2 = ref2het2_iter[vec_idx + 1];
      VecW count2_ssq = nm1 & hom2;
      VecW count2_nm = nm1 & (hom2 | ref2het2);
      VecW count2_plus = count2_ssq & ref2het2;

      nm1 = hom1_iter[vec_idx + 2] | ref2het1_iter[vec_idx + 2];
      hom2 = hom2_iter[vec_idx + 2];
      ref2het2 = ref2het2_iter[vec_idx + 2];
      VecW half_a_ssq = nm1 & hom2;
      VecW half_a_nm = nm1 & (hom2 | ref2het2);
      VecW half_a_plus = half_a_ssq & ref2het2;
      const VecW half_b_ssq = vecw_srli(half_a_ssq, 1) & m1;
      const VecW half_b_nm = vecw_srli(half_a_nm, 1) & m1;
      const VecW half_b_plus = vecw_srli(half_a_plus, 1) & m1;
      half_a_ssq = half_a_ssq & m1;
      half_a_nm = half_a_nm & m1;
      half_a_plus = half_a_plus & m1;
      count1_ssq = count1_ssq - (vecw_srli(count1_ssq, 1) & m1);
      count1_nm = count1_nm - (vecw_srli(count1_nm, 1) & m1);
      count1_plus = count1_plus - (vecw_srli(count1_plus, 1) & m1);
      count2_ssq = count2_ssq - (vecw_srli(count2_ssq, 1) & m1);
      count2_nm = count2_nm - (vecw_srli(count2_nm, 1) & m1);
      count2_plus = count2_plus - (vecw_srli(count2_plus, 1) & m1);
      count1_ssq = count1_ssq + half_a_ssq;
      count1_nm = count1_nm + half_a_nm;
      count1_plus = count1_plus + half_a_plus;
      count2_ssq = count2_ssq + half_b_ssq;
      count2_nm = count2_nm + half_b_nm;
      count2_plus = count2_plus + half_b_plus;
      count1_ssq = (count1_ssq & m2) + (vecw_srli(count1_ssq, 2) & m2);
      count1_nm = (count1_nm & m2) + (vecw_srli(count1_nm, 2) & m2);
      count1_plus = (count1_plus & m2) + (vecw_srli(count1_plus, 2) & m2);
      count1_ssq = count1_ssq + (count2_ssq & m2) + (vecw_srli(count2_ssq, 2) & m2);
      count1_nm = count1_nm + (count2_nm & m2) + (vecw_srli(count2_nm, 2) & m2);
      count1_plus = count1_plus + (count2_plus & m2) + (vecw_srli(count2_plus, 2) & m2);
      inner_acc_nm = inner_acc_nm + (count1_nm & m4) + (vecw_srli(count1_nm, 4) & m4);
      inner_acc_ssq = inner_acc_ssq + (count1_ssq & m4) + (vecw_srli(count1_ssq, 4) & m4);
      inner_acc_plus = inner_acc_plus + (count1_plus & m4) + (vecw_srli(count1_plus, 4) & m4);
    }
    hom1_iter = &(hom1_iter[cur_incr]);
    ref2het1_iter = &(ref2het1_iter[cur_incr]);
    hom2_iter = &(hom2_iter[cur_incr]);
    ref2het2_iter = &(ref2het2_iter[cur_incr]);
    const VecW m0 = vecw_setzero();
    acc_nm = acc_nm + vecw_bytesum(inner_acc_nm, m0);
    acc_ssq2 = acc_ssq2 + vecw_bytesum(inner_acc_ssq, m0);
    acc_plus2 = acc_plus2 + vecw_bytesum(inner_acc_plus, m0);
  }
}
#endif  // __LP64__

void SumSsqNmWords(const uintptr_t* hom1, const uintptr_t* ref2het1, const uintptr_t* hom2, const uintptr_t* ref2het2, uint32_t word_ct, uint32_t* __restrict nm_ptr, int32_t* sum2_ptr, uint32_t* __restrict ssq2_ptr) {
  uint32_t nm = 0;
  uint32_t ssq2 = 0;
  uint32_t plus2 = 0;
  uint32_t widx = 0;
#if defined(USE_AVX2) || !defined(USE_SSE42)
  if (word_ct >= 3 * kWordsPerVec) {
    const uintptr_t block_ct = word_ct / (3 * kWordsPerVec);
    SumSsqNmVecs(R_CAST(const VecW*, hom1), R_CAST(const VecW*, ref2het1), R_CAST(const VecW*, hom2), R_CAST(const VecW*, ref2het2), block_ct * 3, &nm, &ssq2, &plus2);
    widx = block_ct * (3 * kWordsPerVec);
  }
#endif
  for (; widx != word_ct; ++widx) {
    const uintptr_t nm1_word = hom1[widx] | ref2het1[widx];
    const uintptr_t hom2_word = hom2[widx];
    const uintptr_t ref2het2_word = ref2het2[widx];
    nm += PopcountWord(nm1_word & (hom2_word | ref2het2_word));
    const uintptr_t ssq2_word = nm1_word & hom2_word;
    ssq2 += PopcountWord(ssq2_word);
    plus2 += PopcountWord(ssq2_word & ref2het2_word);
  }
  *nm_ptr = nm;
  *sum2_ptr = S_CAST(int32_t, 2 * plus2 - ssq2);  // deliberate overflow
  *ssq2_ptr = ssq2;
}

// er, subcontig_info probably deserves its own type...
void LdPruneNextSubcontig(const uintptr_t* variant_include, const uint32_t* variant_bps, const uint32_t* subcontig_info, const uint32_t* subcontig_thread_assignments, uint32_t prune_window_size, uint32_t thread_idx, uint32_t* subcontig_idx_ptr, uint32_t* subcontig_end_tvidx_ptr, uint32_t* next_window_end_tvidx_ptr, uint32_t* variant_uidx_winstart_ptr, uint32_t* variant_uidx_winend_ptr) {
  uint32_t subcontig_idx = *subcontig_idx_ptr;
  do {
    ++subcontig_idx;
  } while (subcontig_thread_assignments[subcontig_idx] != thread_idx);
  *subcontig_idx_ptr = subcontig_idx;
  const uint32_t subcontig_first_tvidx = *subcontig_end_tvidx_ptr;
  const uint32_t subcontig_len = subcontig_info[3 * subcontig_idx];
  const uint32_t variant_uidx_winstart = subcontig_info[3 * subcontig_idx + 2];
  const uint32_t subcontig_end_tvidx = subcontig_first_tvidx + subcontig_len;
  *subcontig_end_tvidx_ptr = subcontig_end_tvidx;
  if (variant_bps) {
    const uint32_t variant_bp_thresh = variant_bps[variant_uidx_winstart] + prune_window_size;
    uintptr_t variant_uidx_winend_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, variant_uidx_winstart + 1, &variant_uidx_winend_base, &cur_bits);
    uint32_t first_window_len = 1;
    uint32_t variant_uidx_winend;
    do {
      variant_uidx_winend = BitIter1(variant_include, &variant_uidx_winend_base, &cur_bits);
    } while ((variant_bps[variant_uidx_winend] <= variant_bp_thresh) && (++first_window_len < subcontig_len));
    *next_window_end_tvidx_ptr = subcontig_first_tvidx + first_window_len;
    *variant_uidx_winend_ptr = variant_uidx_winend;
  } else {
    *next_window_end_tvidx_ptr = subcontig_first_tvidx + MINV(subcontig_len, prune_window_size);
  }

  *variant_uidx_winstart_ptr = variant_uidx_winstart;
}

void LdPruneNextWindow(const uintptr_t* __restrict variant_include, const uint32_t* __restrict variant_bps, const uint32_t* __restrict tvidxs, const uintptr_t* __restrict cur_window_removed, uint32_t prune_window_size, uint32_t window_incr, uint32_t window_maxl, uint32_t subcontig_end_tvidx, uint32_t* cur_window_size_ptr, uint32_t* __restrict window_start_tvidx_ptr, uint32_t* __restrict variant_uidx_winstart_ptr, uint32_t* __restrict next_window_end_tvidx_ptr, uint32_t* __restrict variant_uidx_winend_ptr, uintptr_t* __restrict occupied_window_slots, uint32_t* winpos_to_slot_idx) {
  uint32_t next_window_end_tvidx = *next_window_end_tvidx_ptr;
  if (next_window_end_tvidx == subcontig_end_tvidx) {
    // just completed last window in subcontig
    *cur_window_size_ptr = 0;
    *window_start_tvidx_ptr = subcontig_end_tvidx;
    ZeroWArr(window_maxl, occupied_window_slots);
    return;
  }
  uint32_t next_window_start_tvidx = *window_start_tvidx_ptr;
  if (variant_bps) {
    // this is guaranteed to be nonnegative
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, *variant_uidx_winstart_ptr + 1, &variant_uidx_base, &cur_bits);
    uint32_t variant_uidx_winend = *variant_uidx_winend_ptr;
    const uint32_t window_start_min_bp = variant_bps[variant_uidx_winend] - prune_window_size;
    uint32_t window_start_bp;
    uint32_t variant_uidx_winstart;
    do {
      // advance window start by as much as necessary to make end advance by at
      // least 1
      ++next_window_start_tvidx;
      variant_uidx_winstart = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      window_start_bp = variant_bps[variant_uidx_winstart];
    } while (window_start_bp < window_start_min_bp);
    // now advance window end as appropriate
    const uint32_t window_end_thresh = window_start_bp + prune_window_size;
    BitIter1Start(variant_include, variant_uidx_winend + 1, &variant_uidx_base, &cur_bits);
    do {
      if (++next_window_end_tvidx == subcontig_end_tvidx) {
        break;
      }
      variant_uidx_winend = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    } while (variant_bps[variant_uidx_winend] <= window_end_thresh);
    *variant_uidx_winstart_ptr = variant_uidx_winstart;
    *variant_uidx_winend_ptr = variant_uidx_winend;
  } else {
    next_window_start_tvidx += window_incr;
    next_window_end_tvidx = MINV(next_window_start_tvidx + prune_window_size, subcontig_end_tvidx);
  }
  const uint32_t cur_window_size = *cur_window_size_ptr;
  uint32_t winpos_write = 0;
  for (uint32_t winpos_read = 0; winpos_read != cur_window_size; ++winpos_read) {
    const uint32_t slot_idx = winpos_to_slot_idx[winpos_read];
    if (IsSet(cur_window_removed, winpos_read) || (tvidxs[slot_idx] < next_window_start_tvidx)) {
      ClearBit(slot_idx, occupied_window_slots);
    } else {
      winpos_to_slot_idx[winpos_write++] = slot_idx;
    }
  }
  *cur_window_size_ptr = winpos_write;
  *window_start_tvidx_ptr = next_window_start_tvidx;
  *next_window_end_tvidx_ptr = next_window_end_tvidx;
}

typedef struct VariantAggsStruct {
  uint32_t nm_ct;
  int32_t sum;
  uint32_t ssq;
} VariantAggs;

// On entry, {cur_nm_ct, cur_first_sum, cur_first_ssq} must be initialized to
// first_vaggs values.
void ComputeIndepPairwiseR2Components(const uintptr_t* __restrict first_genobufs, const uintptr_t* __restrict second_genobufs, const VariantAggs* second_vaggs, uint32_t founder_ct, uint32_t* cur_nm_ct_ptr, int32_t* cur_first_sum_ptr, uint32_t* cur_first_ssq_ptr, int32_t* second_sum_ptr, uint32_t* second_ssq_ptr, int32_t* cur_dotprod_ptr) {
  const uint32_t founder_ctaw = BitCtToAlignedWordCt(founder_ct);
  const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
  // Three cases:
  // 1. Just need dot product.
  // 2. Also need {first_sum, first_ssq} xor {second_sum, second_ssq}.
  // 3. Need all six variables.
  *cur_dotprod_ptr = DotprodWords(first_genobufs, &(first_genobufs[founder_ctaw]), second_genobufs, &(second_genobufs[founder_ctaw]), founder_ctl);
  if (*cur_nm_ct_ptr != founder_ct) {
    SumSsqWords(first_genobufs, &(first_genobufs[founder_ctaw]), second_genobufs, &(second_genobufs[founder_ctaw]), founder_ctl, second_sum_ptr, second_ssq_ptr);
  } else {
    *second_sum_ptr = second_vaggs->sum;
    *second_ssq_ptr = second_vaggs->ssq;
  }
  const uint32_t second_nm_ct = second_vaggs->nm_ct;
  if (second_nm_ct == founder_ct) {
    return;
  }
  if (*cur_nm_ct_ptr != founder_ct) {
    SumSsqNmWords(second_genobufs, &(second_genobufs[founder_ctaw]), first_genobufs, &(first_genobufs[founder_ctaw]), founder_ctl, cur_nm_ct_ptr, cur_first_sum_ptr, cur_first_ssq_ptr);
  } else {
    SumSsqWords(second_genobufs, &(second_genobufs[founder_ctaw]), first_genobufs, &(first_genobufs[founder_ctaw]), founder_ctl, cur_first_sum_ptr, cur_first_ssq_ptr);
    *cur_nm_ct_ptr = second_nm_ct;
  }
}

void FillVaggs(const uintptr_t* hom_vec, const uintptr_t* ref2het_vec, uintptr_t word_ct, VariantAggs* vaggs, uint32_t* nm_ct_ptr, uint32_t* plusone_ct_ptr, uint32_t* minusone_ct_ptr) {
  uint32_t hom_ct;
  uint32_t ref2het_ct;
  uint32_t ref2_ct;
  PopcountWordsIntersect3val(hom_vec, ref2het_vec, word_ct, &hom_ct, &ref2het_ct, &ref2_ct);
  const uint32_t alt2_ct = hom_ct - ref2_ct;
  const uint32_t nm_ct = alt2_ct + ref2het_ct;
  vaggs->nm_ct = nm_ct;
  vaggs->sum = S_CAST(int32_t, ref2_ct - alt2_ct);  // deliberate overflow
  vaggs->ssq = hom_ct;
  *nm_ct_ptr = nm_ct;
  *plusone_ct_ptr = ref2_ct;
  *minusone_ct_ptr = alt2_ct;
}

void IndepPairwiseUpdateSubcontig(uint32_t variant_uidx_winstart, uint32_t x_start, uint32_t x_len, uint32_t y_start, uint32_t y_len, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t founder_nonfemale_ct, uint32_t* is_x_ptr, uint32_t* is_y_ptr, uint32_t* cur_founder_ct_ptr, uint32_t* cur_founder_ctaw_ptr, uint32_t* cur_founder_ctl_ptr, uintptr_t* entire_variant_buf_word_ct_ptr) {
  // _len is better than _end here since we can exploit unsignedness
  const uint32_t is_x = ((variant_uidx_winstart - x_start) < x_len);
  const uint32_t is_y = ((variant_uidx_winstart - y_start) < y_len);
  if ((is_x != (*is_x_ptr)) || (is_y != (*is_y_ptr))) {
    *is_x_ptr = is_x;
    *is_y_ptr = is_y;
    uint32_t cur_founder_ct = founder_ct;
    if (is_x) {
      cur_founder_ct = founder_male_ct;
    } else if (is_y) {
      cur_founder_ct = founder_nonfemale_ct;
    }
    const uint32_t cur_founder_ctaw = BitCtToAlignedWordCt(cur_founder_ct);
    *cur_founder_ct_ptr = cur_founder_ct;
    *cur_founder_ctaw_ptr = cur_founder_ctaw;
    *cur_founder_ctl_ptr = BitCtToWordCt(cur_founder_ct);
    *entire_variant_buf_word_ct_ptr = 2 * cur_founder_ctaw;
    if (is_x) {
      *entire_variant_buf_word_ct_ptr += 2 * BitCtToAlignedWordCt(founder_ct - founder_male_ct);
    }
  }
}

typedef struct IndepPairwiseCtxStruct {
  const uint32_t* subcontig_info;
  const uint32_t* subcontig_thread_assignments;
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const AlleleCode* maj_alleles;
  const double* all_allele_freqs;
  const uint32_t* variant_bps;
  const uintptr_t* preferred_variants;
  uint32_t* tvidx_end;
  uint32_t x_start;
  uint32_t x_len;
  uint32_t y_start;
  uint32_t y_len;
  uint32_t founder_ct;
  uint32_t founder_male_ct;
  uint32_t founder_nonfemale_ct;
  uint32_t prune_window_size;
  uint32_t window_maxl;
  double prune_ld_thresh;
  uint32_t window_incr;
  uint32_t cur_batch_size;

  uintptr_t** genobufs;
  uintptr_t** occupied_window_slots;
  uintptr_t** cur_window_removed;
  double** cur_maj_freqs;
  uintptr_t** removed_variants_write;
  VariantAggs** vaggs;
  VariantAggs** nonmale_vaggs;
  uint32_t** winpos_to_slot_idx;
  uint32_t** tvidxs;
  uint32_t** first_unchecked_tvidx;
  // 't' stands for thread here
  uintptr_t** raw_tgenovecs[2];
} IndepPairwiseCtx;

THREAD_FUNC_DECL IndepPairwiseThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  IndepPairwiseCtx* ctx = S_CAST(IndepPairwiseCtx*, arg->sharedp->context);

  const uint32_t* subcontig_info = ctx->subcontig_info;
  const uint32_t* subcontig_thread_assignments = ctx->subcontig_thread_assignments;
  const uintptr_t* variant_include = ctx->variant_include;
  const uintptr_t* preferred_variants = ctx->preferred_variants;
  const uint32_t x_start = ctx->x_start;
  const uint32_t x_len = ctx->x_len;
  const uint32_t y_start = ctx->y_start;
  const uint32_t y_len = ctx->y_len;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const AlleleCode* maj_alleles = ctx->maj_alleles;
  const double* all_allele_freqs = ctx->all_allele_freqs;
  const uint32_t* variant_bps = ctx->variant_bps;
  const uint32_t founder_ct = ctx->founder_ct;
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t founder_male_ctl2 = NypCtToWordCt(founder_male_ct);
  const uint32_t founder_nonfemale_ct = ctx->founder_nonfemale_ct;
  const uint32_t nonmale_ct = founder_ct - founder_male_ct;
  const uint32_t nonmale_ctaw = BitCtToAlignedWordCt(nonmale_ct);
  const uint32_t nonmale_ctl = BitCtToWordCt(nonmale_ct);
  const uintptr_t raw_tgenovec_single_variant_word_ct = RoundUpPow2(NypCtToWordCt(nonmale_ct) + founder_male_ctl2, kWordsPerVec);
  const uint32_t prune_window_size = ctx->prune_window_size;
  const uint32_t window_maxl = ctx->window_maxl;
  const double prune_ld_thresh = ctx->prune_ld_thresh;
  const uint32_t window_incr = ctx->window_incr;
  const uint32_t tvidx_end = ctx->tvidx_end[tidx];
  uintptr_t* genobufs = ctx->genobufs[tidx];
  uintptr_t* occupied_window_slots = ctx->occupied_window_slots[tidx];
  uintptr_t* cur_window_removed = ctx->cur_window_removed[tidx];
  uintptr_t* removed_variants_write = ctx->removed_variants_write[tidx];
  double* cur_maj_freqs = ctx->cur_maj_freqs[tidx];
  VariantAggs* vaggs = ctx->vaggs[tidx];
  VariantAggs* nonmale_vaggs = ctx->nonmale_vaggs[tidx];
  uint32_t* winpos_to_slot_idx = ctx->winpos_to_slot_idx[tidx];
  uint32_t* tvidxs = ctx->tvidxs[tidx];
  uint32_t* first_unchecked_tvidx = ctx->first_unchecked_tvidx? ctx->first_unchecked_tvidx[tidx] : nullptr;

  uint32_t subcontig_end_tvidx = 0;
  uint32_t subcontig_idx = UINT32_MAX;  // deliberate overflow
  uint32_t window_start_tvidx = 0;
  uint32_t next_window_end_tvidx = 0;
  uint32_t write_slot_idx = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t cur_window_size = 0;
  uint32_t winpos_split = 0;
  uint32_t tvidx_start = 0;
  uint32_t cur_founder_ct = founder_ct;
  uint32_t cur_founder_ctaw = BitCtToAlignedWordCt(founder_ct);
  uint32_t cur_founder_ctl = BitCtToWordCt(founder_ct);
  uintptr_t variant_uidx_base = 0;
  uintptr_t variant_include_bits = variant_include[0];
  uint32_t variant_uidx_winstart = 0;
  uint32_t variant_uidx_winend = 0;
  uintptr_t entire_variant_buf_word_ct = 2 * cur_founder_ctaw;
  uint32_t cur_allele_ct = 2;
  uint32_t parity = 0;
  do {
    const uint32_t cur_batch_size = ctx->cur_batch_size;
    const uint32_t tvidx_stop = MINV(tvidx_start + cur_batch_size, tvidx_end);
    // main loop has to be variant-, not window-, based due to how datasets too
    // large to fit in memory are handled: we may have to halt in the middle of
    // unpacking data for a window, waiting until the current I/O pass is
    // complete before proceeding
    const uintptr_t* raw_tgenovecs = ctx->raw_tgenovecs[parity][tidx];
    for (uint32_t cur_tvidx = tvidx_start; cur_tvidx < tvidx_stop; ) {
      if (cur_tvidx == subcontig_end_tvidx) {
        LdPruneNextSubcontig(variant_include, variant_bps, subcontig_info, subcontig_thread_assignments, prune_window_size, tidx, &subcontig_idx, &subcontig_end_tvidx, &next_window_end_tvidx, &variant_uidx_winstart, &variant_uidx_winend);
        IndepPairwiseUpdateSubcontig(variant_uidx_winstart, x_start, x_len, y_start, y_len, founder_ct, founder_male_ct, founder_nonfemale_ct, &is_x, &is_y, &cur_founder_ct, &cur_founder_ctaw, &cur_founder_ctl, &entire_variant_buf_word_ct);
        BitIter1Start(variant_include, variant_uidx_winstart, &variant_uidx_base, &variant_include_bits);
        winpos_split = 0;
      }
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      write_slot_idx = AdvTo0Bit(occupied_window_slots, write_slot_idx);
      uintptr_t tvidx_offset = cur_tvidx - tvidx_start;
      const uintptr_t* cur_raw_tgenovecs = &(raw_tgenovecs[tvidx_offset * raw_tgenovec_single_variant_word_ct]);
      uintptr_t* cur_genobuf = &(genobufs[write_slot_idx * entire_variant_buf_word_ct]);
      uintptr_t* cur_genobuf_ref2het = &(cur_genobuf[cur_founder_ctaw]);
      // need local genobuf anyway due to halts, so may as well perform split
      // here.
      SplitHomRef2het(cur_raw_tgenovecs, cur_founder_ct, cur_genobuf, cur_genobuf_ref2het);
      uint32_t nm_ct;
      uint32_t plusone_ct;
      uint32_t minusone_ct;
      FillVaggs(cur_genobuf, cur_genobuf_ref2het, cur_founder_ctl, &(vaggs[write_slot_idx]), &nm_ct, &plusone_ct, &minusone_ct);
      if (is_x) {
        cur_genobuf = &(cur_genobuf[2 * cur_founder_ctaw]);
        cur_genobuf_ref2het = &(cur_genobuf[nonmale_ctaw]);
        SplitHomRef2het(&(cur_raw_tgenovecs[founder_male_ctl2]), nonmale_ct, cur_genobuf, cur_genobuf_ref2het);
        uint32_t x_nonmale_nm_ct;
        uint32_t x_nonmale_plusone_ct;
        uint32_t x_nonmale_minusone_ct;
        FillVaggs(cur_genobuf, cur_genobuf_ref2het, nonmale_ctl, &(nonmale_vaggs[write_slot_idx]), &x_nonmale_nm_ct, &x_nonmale_plusone_ct, &x_nonmale_minusone_ct);
        nm_ct += 2 * x_nonmale_nm_ct;
        plusone_ct += 2 * x_nonmale_plusone_ct;
        minusone_ct += 2 * x_nonmale_minusone_ct;
      }
      if (((!plusone_ct) && (!minusone_ct)) || (plusone_ct == nm_ct) || (minusone_ct == nm_ct)) {
        SetBit(cur_window_size, cur_window_removed);
        SetBit(cur_tvidx, removed_variants_write);
      } else {
        tvidxs[write_slot_idx] = cur_tvidx;
        uintptr_t allele_idx_base;
        if (!allele_idx_offsets) {
          allele_idx_base = variant_uidx;
        } else {
          allele_idx_base = allele_idx_offsets[variant_uidx];
          cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_base;
          allele_idx_base -= variant_uidx;
        }
        cur_maj_freqs[write_slot_idx] = GetAlleleFreq(&(all_allele_freqs[allele_idx_base]), maj_alleles[variant_uidx], cur_allele_ct);
        if (preferred_variants && IsSet(preferred_variants, variant_uidx)) {
          cur_maj_freqs[write_slot_idx] -= 1.0;
        }
        if (first_unchecked_tvidx) {
          first_unchecked_tvidx[write_slot_idx] = cur_tvidx + 1;
        }
      }
      SetBit(write_slot_idx, occupied_window_slots);
      winpos_to_slot_idx[cur_window_size++] = write_slot_idx;
      ++cur_tvidx;
      // are we at the end of a window?  if not, load more variant(s) before
      // proceeding.
      if (cur_tvidx != next_window_end_tvidx) {
        continue;
      }
      if (first_unchecked_tvidx) {
        // PLINK 1.x pruning order

        // possible for cur_window_size == 1, if all variants at the end of the
        // previous window were pruned
        uint32_t cur_removed_ct = PopcountWords(cur_window_removed, BitCtToWordCt(cur_window_size));
        uint32_t prev_removed_ct;
        do {
          prev_removed_ct = cur_removed_ct;
          // const uint32_t debug_print = (!IsSet(cur_window_removed, 0)) && (tvidxs[winpos_to_slot_idx[0]] == 0);
          for (uint32_t first_winpos = 0; ; ++first_winpos) {
            // can't use BitIter0 since we care about changes in this loop to
            // cur_window_removed
            first_winpos = AdvTo0Bit(cur_window_removed, first_winpos);
            // can assume empty trailing bit for cur_window_removed
            if (first_winpos == cur_window_size) {
              break;
            }
            const uint32_t first_slot_idx = winpos_to_slot_idx[first_winpos];
            const uint32_t cur_first_unchecked_tvidx = first_unchecked_tvidx[first_slot_idx];
            if (cur_first_unchecked_tvidx == cur_tvidx) {
              continue;
            }
            // safe to use BitIter0 for second_winpos, though
            uintptr_t second_winpos_base;
            uintptr_t cur_window_removed_inv_bits;
            BitIter0Start(cur_window_removed, first_winpos + 1, &second_winpos_base, &cur_window_removed_inv_bits);
            {
              uint32_t second_winpos;
              uint32_t second_slot_idx;
              do {
                second_winpos = BitIter0(cur_window_removed, &second_winpos_base, &cur_window_removed_inv_bits);
                if (second_winpos == cur_window_size) {
                  first_unchecked_tvidx[first_slot_idx] = cur_tvidx;
                  goto IndepPairwiseThread_next_first;
                }
                second_slot_idx = winpos_to_slot_idx[second_winpos];
              } while (tvidxs[second_slot_idx] < cur_first_unchecked_tvidx);
              const uintptr_t* first_genobufs = &(genobufs[first_slot_idx * entire_variant_buf_word_ct]);
              const uint32_t first_nm_ct = vaggs[first_slot_idx].nm_ct;
              const int32_t first_sum = vaggs[first_slot_idx].sum;
              const uint32_t first_ssq = vaggs[first_slot_idx].ssq;
              while (1) {
                const uintptr_t* second_genobufs = &(genobufs[second_slot_idx * entire_variant_buf_word_ct]);
                uint32_t cur_nm_ct = first_nm_ct;
                int32_t cur_first_sum = first_sum;
                uint32_t cur_first_ssq = first_ssq;
                int32_t second_sum;
                uint32_t second_ssq;
                int32_t cur_dotprod;
                ComputeIndepPairwiseR2Components(first_genobufs, second_genobufs, &(vaggs[second_slot_idx]), cur_founder_ct, &cur_nm_ct, &cur_first_sum, &cur_first_ssq, &second_sum, &second_ssq, &cur_dotprod);
                if (is_x) {
                  uint32_t nonmale_nm_ct = nonmale_vaggs[first_slot_idx].nm_ct;
                  int32_t nonmale_first_sum = nonmale_vaggs[first_slot_idx].sum;
                  uint32_t nonmale_first_ssq = nonmale_vaggs[first_slot_idx].ssq;
                  int32_t nonmale_dotprod;
                  int32_t nonmale_second_sum;
                  uint32_t nonmale_second_ssq;
                  ComputeIndepPairwiseR2Components(&(first_genobufs[2 * cur_founder_ctaw]), &(second_genobufs[2 * cur_founder_ctaw]), &(nonmale_vaggs[second_slot_idx]), nonmale_ct, &nonmale_nm_ct, &nonmale_first_sum, &nonmale_first_ssq, &nonmale_second_sum, &nonmale_second_ssq, &nonmale_dotprod);
                  // only --ld-xchr 3 for now
                  // assumes founder_ct < 2^30
                  cur_nm_ct += 2 * nonmale_nm_ct;
                  cur_first_sum += 2 * nonmale_first_sum;
                  cur_first_ssq += 2 * nonmale_first_ssq;
                  second_sum += 2 * nonmale_second_sum;
                  second_ssq += 2 * nonmale_second_ssq;
                  cur_dotprod += 2 * nonmale_dotprod;
                }
                // these three values are actually cur_nm_ct times their
                // true values, but that cancels out
                const double cov12 = S_CAST(double, cur_dotprod * S_CAST(int64_t, cur_nm_ct) - S_CAST(int64_t, cur_first_sum) * second_sum);
                const double variance1 = S_CAST(double, cur_first_ssq * S_CAST(int64_t, cur_nm_ct) - S_CAST(int64_t, cur_first_sum) * cur_first_sum);
                const double variance2 = S_CAST(double, second_ssq * S_CAST(int64_t, cur_nm_ct) - S_CAST(int64_t, second_sum) * second_sum);
                // > instead of >=, so we don't prune from a pair of
                // variants with zero common observations
                if (cov12 * cov12 > prune_ld_thresh * variance1 * variance2) {
                  // this has a surprisingly large ~3% speed penalty on my
                  // main test scenario, but that's an acceptable price to
                  // pay for reproducibility.
                  if (cur_maj_freqs[first_slot_idx] > cur_maj_freqs[second_slot_idx] * (1 + kSmallEpsilon)) {
                    SetBit(first_winpos, cur_window_removed);
                    SetBit(tvidxs[first_slot_idx], removed_variants_write);
                  } else {
                    SetBit(second_winpos, cur_window_removed);
                    SetBit(tvidxs[second_slot_idx], removed_variants_write);
                    const uint32_t next_start_winpos = BitIter0NoAdv(cur_window_removed, &second_winpos_base, &cur_window_removed_inv_bits);
                    if (next_start_winpos < cur_window_size) {
                      first_unchecked_tvidx[first_slot_idx] = tvidxs[winpos_to_slot_idx[next_start_winpos]];
                    } else {
                      first_unchecked_tvidx[first_slot_idx] = cur_tvidx;
                    }
                  }
                  break;
                }
                second_winpos = BitIter0(cur_window_removed, &second_winpos_base, &cur_window_removed_inv_bits);
                if (second_winpos == cur_window_size) {
                  first_unchecked_tvidx[first_slot_idx] = cur_tvidx;
                  break;
                }
                second_slot_idx = winpos_to_slot_idx[second_winpos];
              }  // while (1)
            }
          IndepPairwiseThread_next_first:
            ;
          }
          cur_removed_ct = PopcountWords(cur_window_removed, BitCtToWordCt(cur_window_size));
        } while (cur_removed_ct > prev_removed_ct);
      } else {
        // Within each window, scan in reverse order.  This way, we tend to
        // check the nearest new pairs first, and this should allow us to exit
        // early more often.
        const uint32_t second_winpos_stop = winpos_split? winpos_split : 1;
        for (uint32_t second_winpos = cur_window_size; second_winpos != second_winpos_stop; ) {
          --second_winpos;
          const uint32_t second_slot_idx = winpos_to_slot_idx[second_winpos];
          const uintptr_t* second_genobufs = &(genobufs[second_slot_idx * entire_variant_buf_word_ct]);
          const uint32_t second_nm_ct = vaggs[second_slot_idx].nm_ct;
          const int32_t second_sum = vaggs[second_slot_idx].sum;
          const uint32_t second_ssq = vaggs[second_slot_idx].ssq;
          for (uint32_t first_winpos = second_winpos; first_winpos; ) {
            --first_winpos;
            // possible todo: faster unset-bit reverse-iterator.  but probably
            // doesn't pay off here.
            if (IsSet(cur_window_removed, first_winpos)) {
              continue;
            }
            const uint32_t first_slot_idx = winpos_to_slot_idx[first_winpos];
            const uintptr_t* first_genobufs = &(genobufs[first_slot_idx * entire_variant_buf_word_ct]);
            uint32_t cur_nm_ct = second_nm_ct;
            int32_t cur_second_sum = second_sum;
            uint32_t cur_second_ssq = second_ssq;
            int32_t first_sum;
            uint32_t first_ssq;
            int32_t cur_dotprod;
            ComputeIndepPairwiseR2Components(second_genobufs, first_genobufs, &(vaggs[first_slot_idx]), cur_founder_ct, &cur_nm_ct, &cur_second_sum, &cur_second_ssq, &first_sum, &first_ssq, &cur_dotprod);
            if (is_x) {
              uint32_t nonmale_nm_ct = nonmale_vaggs[second_slot_idx].nm_ct;
              int32_t nonmale_second_sum = nonmale_vaggs[second_slot_idx].sum;
              uint32_t nonmale_second_ssq = nonmale_vaggs[second_slot_idx].ssq;
              int32_t nonmale_dotprod;
              int32_t nonmale_first_sum;
              uint32_t nonmale_first_ssq;
              ComputeIndepPairwiseR2Components(&(second_genobufs[2 * cur_founder_ctaw]), &(first_genobufs[2 * cur_founder_ctaw]), &(nonmale_vaggs[first_slot_idx]), nonmale_ct, &nonmale_nm_ct, &nonmale_second_sum, &nonmale_second_ssq, &nonmale_first_sum, &nonmale_first_ssq, &nonmale_dotprod);
              // only --ld-xchr 3 for now
              // assumes founder_ct < 2^30
              cur_nm_ct += 2 * nonmale_nm_ct;
              first_sum += 2 * nonmale_first_sum;
              first_ssq += 2 * nonmale_first_ssq;
              cur_second_sum += 2 * nonmale_second_sum;
              cur_second_ssq += 2 * nonmale_second_ssq;
              cur_dotprod += 2 * nonmale_dotprod;
            }
            // these three values are actually cur_nm_ct times their
            // true values, but that cancels out
            const double cov12 = S_CAST(double, cur_dotprod * S_CAST(int64_t, cur_nm_ct) - S_CAST(int64_t, first_sum) * cur_second_sum);
            const double variance1 = S_CAST(double, first_ssq * S_CAST(int64_t, cur_nm_ct) - S_CAST(int64_t, first_sum) * first_sum);
            const double variance2 = S_CAST(double, cur_second_ssq * S_CAST(int64_t, cur_nm_ct) - S_CAST(int64_t, cur_second_sum) * cur_second_sum);
            // > instead of >=, so we don't prune from a pair of
            // variants with zero common observations
            if (cov12 * cov12 > prune_ld_thresh * variance1 * variance2) {
              if (cur_maj_freqs[first_slot_idx] <= cur_maj_freqs[second_slot_idx] * (1 + kSmallEpsilon)) {
                SetBit(second_winpos, cur_window_removed);
                SetBit(tvidxs[second_slot_idx], removed_variants_write);
                break;
              }
              SetBit(first_winpos, cur_window_removed);
              SetBit(tvidxs[first_slot_idx], removed_variants_write);
            }
          }  // while (1)
        }
      }
      const uint32_t prev_window_size = cur_window_size;
      LdPruneNextWindow(variant_include, variant_bps, tvidxs, cur_window_removed, prune_window_size, window_incr, window_maxl, subcontig_end_tvidx, &cur_window_size, &window_start_tvidx, &variant_uidx_winstart, &next_window_end_tvidx, &variant_uidx_winend, occupied_window_slots, winpos_to_slot_idx);
      winpos_split = cur_window_size;
      // clear bits here since we set cur_window_removed bits during loading
      // process in monomorphic case
      ZeroWArr(BitCtToWordCt(prev_window_size), cur_window_removed);
      write_slot_idx = 0;
    }
    parity = 1 - parity;
    tvidx_start = tvidx_stop;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr IndepPairwise(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uint32_t* founder_info_cumulative_popcounts, const uintptr_t* founder_nonmale, const uintptr_t* founder_male, const uintptr_t* founder_nonfemale, const LdInfo* ldip, const uintptr_t* preferred_variants, const uint32_t* subcontig_info, const uint32_t* subcontig_thread_assignments, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t founder_nonfemale_ct, uint32_t subcontig_ct, uintptr_t window_max, uint32_t calc_thread_ct, uint32_t max_load, PgenReader* simple_pgrp, uintptr_t* removed_variants_collapsed) {
  ThreadGroup tg;
  PreinitThreads(&tg);
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t founder_nonmale_ct = founder_ct - founder_male_ct;
    if (unlikely(founder_nonmale_ct * 2 + founder_male_ct > 0x7fffffffU)) {
      // may as well document this
      logerrputs("Error: --indep-pairwise does not support >= 2^30 founders.\n");
      goto IndepPairwise_ret_NOT_YET_SUPPORTED;
    }
    const uint32_t founder_nonmale_ctaw = BitCtToAlignedWordCt(founder_nonmale_ct);
    const uint32_t founder_male_ctaw = BitCtToAlignedWordCt(founder_male_ct);
    // Per-thread allocations:
    // - tvidx_batch_size * raw_tgenovec_single_variant_word_ct *
    //     sizeof(intptr_t) for raw genotype data (raw_tgenovecs)
    // - tvidx_batch_size * sizeof(double) for cur_maj_freqs
    // - if pos-based window, tvidx_batch_size * sizeof(int32_t)
    // - All of the above again, to allow loader thread to operate
    //     independently
    // - window_max * 2 * (founder_nonmale_ctaw + founder_male_ctaw) *
    //     kBytesPerVec for split genotype data
    // - max_loadl * sizeof(intptr_t) for removed-variant bitarray
    // - window_max * 3 * sizeof(int32_t) for main missing_ct, sum(x_i),
    //     sum(x_i^2) array
    // - window_max * 3 * sizeof(int32_t) for chrX founder_male missing_ct,
    //     sum(x_i), sum(x_i^2) array
    // - window_max * sizeof(int32_t) for indexes into genotype data bitarrays
    //     (for now, anyway)
    // - window_max * sizeof(int32_t) for live_indices (variant_idxs?)
    // - window_max * sizeof(int32_t) for start_arr (first uncompared
    //     variant_idx)
    IndepPairwiseCtx ctx;
    uintptr_t* tmp_genovec;
    uint32_t* thread_last_subcontig;
    uint32_t* thread_subcontig_start_tvidx;
    uint32_t* thread_last_tvidx;
    uint32_t* thread_last_uidx;
    if (unlikely(bigstack_alloc_w(NypCtToWordCt(raw_sample_ct), &tmp_genovec) ||
                 bigstack_calloc_u32(calc_thread_ct, &ctx.tvidx_end) ||
                 bigstack_calloc_u32(calc_thread_ct, &thread_last_subcontig) ||
                 bigstack_calloc_u32(calc_thread_ct, &thread_subcontig_start_tvidx) ||
                 bigstack_calloc_u32(calc_thread_ct, &thread_last_tvidx) ||
                 bigstack_calloc_u32(calc_thread_ct, &thread_last_uidx) ||
                 bigstack_alloc_wp(calc_thread_ct, &ctx.genobufs) ||
                 bigstack_alloc_wp(calc_thread_ct, &ctx.occupied_window_slots) ||
                 bigstack_alloc_wp(calc_thread_ct, &ctx.cur_window_removed) ||
                 bigstack_alloc_dp(calc_thread_ct, &ctx.cur_maj_freqs) ||
                 bigstack_alloc_wp(calc_thread_ct, &ctx.removed_variants_write) ||
                 BIGSTACK_ALLOC_X(VariantAggs*, calc_thread_ct, &ctx.vaggs) ||
                 BIGSTACK_ALLOC_X(VariantAggs*, calc_thread_ct, &ctx.nonmale_vaggs) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.winpos_to_slot_idx) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.tvidxs) ||
                 bigstack_alloc_wp(calc_thread_ct, &(ctx.raw_tgenovecs[0])) ||
                 bigstack_alloc_wp(calc_thread_ct, &(ctx.raw_tgenovecs[1])))) {
      goto IndepPairwise_ret_NOMEM;
    }
    const uint32_t plink1_order = (ldip->prune_flags / kfLdPrunePlink1Order) & 1;
    if (plink1_order) {
      if (unlikely(bigstack_alloc_u32p(calc_thread_ct, &ctx.first_unchecked_tvidx))) {
        goto IndepPairwise_ret_NOMEM;
      }
    } else {
      ctx.first_unchecked_tvidx = nullptr;
    }
    for (uint32_t subcontig_idx = 0; subcontig_idx != subcontig_ct; ++subcontig_idx) {
      const uint32_t cur_thread_idx = subcontig_thread_assignments[subcontig_idx];
      ctx.tvidx_end[cur_thread_idx] += subcontig_info[3 * subcontig_idx];
    }
    const uintptr_t entire_variant_buf_word_ct = 2 * (founder_nonmale_ctaw + founder_male_ctaw);
    const uint32_t window_maxl = BitCtToWordCt(window_max);
    const uint32_t max_loadl = BitCtToWordCt(max_load);
    const uintptr_t genobuf_alloc = RoundUpPow2(window_max * entire_variant_buf_word_ct * sizeof(intptr_t), kCacheline);
    const uintptr_t occupied_window_slots_alloc = RoundUpPow2(window_maxl * sizeof(intptr_t), kCacheline);
    const uintptr_t cur_window_removed_alloc = RoundUpPow2((1 + window_max / kBitsPerWord) * sizeof(intptr_t), kCacheline);
    const uintptr_t cur_maj_freqs_alloc = RoundUpPow2(window_max * sizeof(double), kCacheline);
    const uintptr_t removed_variants_write_alloc = RoundUpPow2(max_loadl * sizeof(intptr_t), kCacheline);

    // two of these
    const uintptr_t vaggs_alloc = RoundUpPow2(window_max * sizeof(VariantAggs), kCacheline);

    // (2 + plink1_order) of these
    const uintptr_t window_int32_alloc = RoundUpPow2(window_max * sizeof(int32_t), kCacheline);

    const uintptr_t thread_alloc_base = genobuf_alloc + occupied_window_slots_alloc + cur_window_removed_alloc + cur_maj_freqs_alloc + removed_variants_write_alloc + 2 * vaggs_alloc + (2 + plink1_order) * window_int32_alloc;

    const uint32_t founder_ctl2 = NypCtToWordCt(founder_ct);
    const uint32_t founder_male_ctl2 = NypCtToWordCt(founder_male_ct);
    const uint32_t founder_nonmale_ctl2 = NypCtToWordCt(founder_nonmale_ct);
    const uint32_t founder_nonfemale_ctl2 = NypCtToWordCt(founder_nonfemale_ct);
    const uintptr_t raw_tgenovec_single_variant_word_ct = RoundUpPow2(founder_nonmale_ctl2 + founder_male_ctl2, kWordsPerVec);
    // round down
    uintptr_t bigstack_avail_per_thread = RoundDownPow2(bigstack_left() / calc_thread_ct, kCacheline);
    // may as well require capacity for >= 256 variants per thread per pass
    if (unlikely(bigstack_avail_per_thread <= thread_alloc_base + 2 * 256 * raw_tgenovec_single_variant_word_ct * sizeof(intptr_t))) {
      goto IndepPairwise_ret_NOMEM;
    }
    bigstack_avail_per_thread -= thread_alloc_base;
    uint32_t tvidx_batch_size = DivUp(max_load, 2);
    // tried a bunch of powers of two, this seems to be a good value
    if (tvidx_batch_size > 65536) {
      tvidx_batch_size = 65536;
    }
    // tvidx_batch_size = max_load;  // temporary debugging
    if (2 * tvidx_batch_size * raw_tgenovec_single_variant_word_ct * sizeof(intptr_t) > bigstack_avail_per_thread) {
      tvidx_batch_size = bigstack_avail_per_thread / RoundUpPow2(raw_tgenovec_single_variant_word_ct * 2 * sizeof(intptr_t), kCacheline);
    }
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      ctx.genobufs[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(genobuf_alloc));
      ctx.occupied_window_slots[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(occupied_window_slots_alloc));
      ZeroWArr(window_maxl, ctx.occupied_window_slots[tidx]);
      ctx.cur_window_removed[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(cur_window_removed_alloc));
      ZeroWArr(1 + window_max / kBitsPerWord, ctx.cur_window_removed[tidx]);
      ctx.cur_maj_freqs[tidx] = S_CAST(double*, bigstack_alloc_raw(cur_maj_freqs_alloc));
      ctx.removed_variants_write[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(removed_variants_write_alloc));
      ZeroWArr(max_loadl, ctx.removed_variants_write[tidx]);
      ctx.vaggs[tidx] = S_CAST(VariantAggs*, bigstack_alloc_raw(vaggs_alloc));
      ctx.nonmale_vaggs[tidx] = S_CAST(VariantAggs*, bigstack_alloc_raw(vaggs_alloc));
      ctx.winpos_to_slot_idx[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(window_int32_alloc));
      ctx.tvidxs[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(window_int32_alloc));
      if (ctx.first_unchecked_tvidx) {
        ctx.first_unchecked_tvidx[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(window_int32_alloc));
      }
      ctx.raw_tgenovecs[0][tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(tvidx_batch_size * raw_tgenovec_single_variant_word_ct * sizeof(intptr_t)));
      ctx.raw_tgenovecs[1][tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(tvidx_batch_size * raw_tgenovec_single_variant_word_ct * sizeof(intptr_t)));
    }
    ctx.subcontig_info = subcontig_info;
    ctx.subcontig_thread_assignments = subcontig_thread_assignments;
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.maj_alleles = maj_alleles;
    ctx.all_allele_freqs = allele_freqs;
    ctx.variant_bps = variant_bps;
    ctx.preferred_variants = preferred_variants;
    ctx.founder_ct = founder_ct;
    ctx.founder_male_ct = founder_male_ct;
    ctx.founder_nonfemale_ct = founder_nonfemale_ct;
    ctx.prune_window_size = ldip->prune_window_size;
    ctx.window_maxl = window_maxl;
    ctx.prune_ld_thresh = ldip->prune_last_param * (1 + kSmallEpsilon);
    ctx.window_incr = ldip->prune_window_incr;
    ctx.cur_batch_size = tvidx_batch_size;
    const uint32_t all_haploid = IsSet(cip->haploid_mask, 0);
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    GetXymtStartAndEnd(cip, kChrOffsetX, &x_start, &x_end);
    GetXymtStartAndEnd(cip, kChrOffsetY, &y_start, &y_end);
    const uint32_t x_len = x_end - x_start;
    const uint32_t y_len = y_end - y_start;
    ctx.x_start = x_start;
    ctx.x_len = x_len;
    ctx.y_start = y_start;
    ctx.y_len = y_len;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto IndepPairwise_ret_NOMEM;
    }
    SetThreadFuncAndData(IndepPairwiseThread, &ctx, &tg);

    // Main workflow:
    // 1. Set n=0, load batch 0

    // 2. Spawn threads processing batch n
    // 3. Increment n by 1
    // 4. Load batch n unless eof
    // 5. Join threads
    // 6. Goto step 2 unless eof
    //
    // 7. Assemble final results with CopyBitarrRange()
    uint32_t parity = 0;
    uint32_t pct = 0;
    uint32_t next_print_tvidx_start = (max_load + 99) / 100;
    logprintf("--indep-pairwise (%u compute thread%s): ", calc_thread_ct, (calc_thread_ct == 1)? "" : "s");
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t cur_tvidx_start = 0; ; cur_tvidx_start += tvidx_batch_size) {
      if (!IsLastBlock(&tg)) {
        PgrSampleSubsetIndex pssi;
        PgrSetSampleSubsetIndex(founder_info_cumulative_popcounts, simple_pgrp, &pssi);
        uintptr_t** cur_raw_tgenovecs = ctx.raw_tgenovecs[parity];
        const uint32_t cur_tvidx_end = cur_tvidx_start + tvidx_batch_size;
        uint32_t is_x_or_y = 0;
        for (uint32_t subcontig_idx = 0; subcontig_idx != subcontig_ct; ++subcontig_idx) {
          const uint32_t cur_thread_idx = subcontig_thread_assignments[subcontig_idx];
          if (thread_last_subcontig[cur_thread_idx] > subcontig_idx) {
            continue;
          }
          uint32_t cur_tvidx = thread_last_tvidx[cur_thread_idx];
          if (cur_tvidx == cur_tvidx_end) {
            continue;
          }
          uint32_t subcontig_start_tvidx = thread_subcontig_start_tvidx[cur_thread_idx];
          uint32_t tvidx_end = subcontig_start_tvidx + subcontig_info[3 * subcontig_idx];
          if (tvidx_end > cur_tvidx_end) {
            tvidx_end = cur_tvidx_end;
            thread_last_subcontig[cur_thread_idx] = subcontig_idx;
          } else {
            thread_subcontig_start_tvidx[cur_thread_idx] = tvidx_end;
            thread_last_subcontig[cur_thread_idx] = subcontig_idx + 1;
          }
          uintptr_t tvidx_offset_end = tvidx_end - cur_tvidx_start;
          uint32_t variant_uidx;
          if (subcontig_start_tvidx == cur_tvidx) {
            variant_uidx = subcontig_info[3 * subcontig_idx + 2];
          } else {
            variant_uidx = thread_last_uidx[cur_thread_idx];
          }
          const uint32_t is_haploid = IsSet(cip->haploid_mask, GetVariantChr(cip, variant_uidx));
          uint32_t is_x = ((variant_uidx - x_start) < x_len);
          const uint32_t new_is_x_or_y = is_x || ((variant_uidx - y_start) < y_len);

          // due to nonempty subset requirement (removed?)
          is_x = is_x && founder_nonmale_ct;
          if (is_x_or_y != new_is_x_or_y) {
            is_x_or_y = new_is_x_or_y;
            if (is_x_or_y) {
              PgrClearSampleSubsetIndex(simple_pgrp, &pssi);
            } else {
              PgrSetSampleSubsetIndex(founder_info_cumulative_popcounts, simple_pgrp, &pssi);
            }
          }
          uintptr_t* cur_thread_raw_tgenovec = cur_raw_tgenovecs[cur_thread_idx];
          uintptr_t variant_uidx_base;
          uintptr_t cur_bits;
          BitIter1Start(variant_include, variant_uidx, &variant_uidx_base, &cur_bits);
          --variant_uidx;
          // todo: document whether tvidx_offset is guaranteed to be <=
          // tvidx_offset_end if this code is revisited.
          for (uintptr_t tvidx_offset = cur_tvidx - cur_tvidx_start; tvidx_offset < tvidx_offset_end; ++tvidx_offset) {
            variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
            uintptr_t* cur_raw_tgenovec = &(cur_thread_raw_tgenovec[tvidx_offset * raw_tgenovec_single_variant_word_ct]);
            // There is no generalization of Pearson r^2 to multiallelic
            // variants with real traction, and after looking at the existing
            // options I don't see a reason for this to change anytime soon.
            // So we always load major allele counts.
            // probable todo: switch to PgrGetDifflistOrGenovec() and have a
            // fast path for low-MAF variants.  Though this isn't *that*
            // important because knowledgeable users will have already filtered
            // out the lowest-MAF variants before starting the LD-prune job.
            if (!is_x_or_y) {
              reterr = PgrGetInv1(founder_info, pssi, founder_ct, variant_uidx, maj_alleles[variant_uidx], simple_pgrp, cur_raw_tgenovec);
              if (unlikely(reterr)) {
                PgenErrPrintNV(reterr, variant_uidx);
                goto IndepPairwise_ret_1;
              }
              if (is_haploid) {
                SetHetMissing(founder_ctl2, cur_raw_tgenovec);
              }
            } else {
              reterr = PgrGetInv1(nullptr, pssi, raw_sample_ct, variant_uidx, maj_alleles[variant_uidx], simple_pgrp, tmp_genovec);
              if (unlikely(reterr)) {
                PgenErrPrintNV(reterr, variant_uidx);
                goto IndepPairwise_ret_1;
              }
              if (is_x) {
                if (founder_male_ct) {
                  CopyNyparrNonemptySubset(tmp_genovec, founder_male, raw_sample_ct, founder_male_ct, cur_raw_tgenovec);
                  SetHetMissing(founder_male_ctl2, cur_raw_tgenovec);
                }
                CopyNyparrNonemptySubset(tmp_genovec, founder_nonmale, raw_sample_ct, founder_nonmale_ct, &(cur_raw_tgenovec[founder_male_ctl2]));
                if (all_haploid) {
                  // don't just treat chrX identically to autosomes, since for
                  // doubled haploids we still want to give females 2x the
                  // weight of males.  I think.
                  SetHetMissing(founder_nonmale_ctl2, &(cur_raw_tgenovec[founder_male_ctl2]));
                }
              } else {
                if (founder_nonfemale_ct) {
                  CopyNyparrNonemptySubset(tmp_genovec, founder_nonfemale, raw_sample_ct, founder_nonfemale_ct, cur_raw_tgenovec);
                  SetHetMissing(founder_nonfemale_ctl2, cur_raw_tgenovec);
                }
              }
            }
          }
          thread_last_tvidx[cur_thread_idx] = tvidx_end;
          thread_last_uidx[cur_thread_idx] = variant_uidx + 1;
        }
      }
      if (cur_tvidx_start) {
        JoinThreads(&tg);
        if (IsLastBlock(&tg)) {
          break;
        }
        if (cur_tvidx_start >= next_print_tvidx_start) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (cur_tvidx_start * 100LLU) / max_load;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_tvidx_start = (pct * S_CAST(uint64_t, max_load) + 99) / 100;
        }
      }
      if (cur_tvidx_start + tvidx_batch_size >= max_load) {
        DeclareLastThreadBlock(&tg);
      }
      if (unlikely(SpawnThreads(&tg))) {
        goto IndepPairwise_ret_THREAD_CREATE_FAIL;
      }
      parity = 1 - parity;
    }
    ZeroU32Arr(calc_thread_ct, thread_subcontig_start_tvidx);
    for (uint32_t subcontig_idx = 0; subcontig_idx != subcontig_ct; ++subcontig_idx) {
      const uint32_t cur_thread_idx = subcontig_thread_assignments[subcontig_idx];
      const uintptr_t* cur_removed_variants = ctx.removed_variants_write[cur_thread_idx];
      const uint32_t subcontig_len = subcontig_info[3 * subcontig_idx];
      const uint32_t subcontig_idx_start = subcontig_info[3 * subcontig_idx + 1];
      CopyBitarrRange(cur_removed_variants, thread_subcontig_start_tvidx[cur_thread_idx], subcontig_idx_start, subcontig_len, removed_variants_collapsed);
      thread_subcontig_start_tvidx[cur_thread_idx] += subcontig_len;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
  }
  while (0) {
  IndepPairwise_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  IndepPairwise_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  IndepPairwise_ret_NOT_YET_SUPPORTED:
    reterr = kPglRetNotYetSupported;
    break;
  }
 IndepPairwise_ret_1:
  CleanupThreads(&tg);
  // caller will free memory
  return reterr;
}

typedef struct VariantHapAggsStruct {
  uint32_t nm_ct;
  uint32_t sum;
} VariantHapAggs;

// On entry, {cur_nm_ct, cur_first_sum} must be initialized to first_vhaggs
// values.
void ComputeIndepPairphaseR2Components(const uintptr_t* __restrict first_hap_vec, const uintptr_t* __restrict second_hap_vec, const VariantHapAggs* second_vhaggs, uint32_t nm_vec_woffset, uint32_t cur_hap_ct, uint32_t* cur_nm_ct_ptr, uint32_t* cur_first_sum_ptr, uint32_t* second_sum_ptr, uint32_t* cur_dotprod_ptr) {
  const uint32_t cur_hap_ctl = BitCtToWordCt(cur_hap_ct);
  // Three cases:
  // 1. Just need dot product.
  // 2. Also need first_sum xor second_sum.
  // 3. Need all four variables.
  *cur_dotprod_ptr = PopcountWordsIntersect(first_hap_vec, second_hap_vec, cur_hap_ctl);
  const uintptr_t* first_nm_vec = &(first_hap_vec[nm_vec_woffset]);
  if (*cur_nm_ct_ptr != cur_hap_ct) {
    *second_sum_ptr = PopcountWordsIntersect(first_nm_vec, second_hap_vec, cur_hap_ctl);
  } else {
    *second_sum_ptr = second_vhaggs->sum;
  }
  const uint32_t second_nm_ct = second_vhaggs->nm_ct;
  if (second_nm_ct == cur_hap_ct) {
    return;
  }
  const uintptr_t* second_nm_vec = &(second_hap_vec[nm_vec_woffset]);
  *cur_first_sum_ptr = PopcountWordsIntersect(first_hap_vec, second_nm_vec, cur_hap_ctl);
  if (*cur_nm_ct_ptr != cur_hap_ct) {
    *cur_nm_ct_ptr = PopcountWordsIntersect(first_nm_vec, second_nm_vec, cur_hap_ctl);
  } else {
    *cur_nm_ct_ptr = second_nm_ct;
  }
}

// Returns 1 if monomorphic, 0 otherwise.
uint32_t FillVhaggs(const uintptr_t* hap_vec, const uintptr_t* nm_vec, uintptr_t word_ct, VariantHapAggs* vhaggs) {
  const uint32_t nm_ct = PopcountWords(nm_vec, word_ct);
  const uint32_t sum = PopcountWords(hap_vec, word_ct);
  vhaggs->nm_ct = nm_ct;
  vhaggs->sum = sum;
  return (!sum) || (sum == nm_ct);
}

void IndepPairphaseUpdateSubcontig(const ChrInfo* cip, uint32_t variant_uidx_winstart, uint32_t x_start, uint32_t x_len, uint32_t y_start, uint32_t y_len, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t founder_nonfemale_ct, uint32_t* is_x_ptr, uint32_t* is_y_ptr, uint32_t* is_haploid_ptr, uint32_t* cur_hap_ct_ptr, uint32_t* cur_hap_ctl_ptr) {
  // _len is better than _end here since we can exploit unsignedness
  const uint32_t is_x = ((variant_uidx_winstart - x_start) < x_len);
  const uint32_t is_y = ((variant_uidx_winstart - y_start) < y_len);
  const uint32_t is_haploid = IsSet(cip->haploid_mask, GetVariantChr(cip, variant_uidx_winstart));
  if ((is_x != (*is_x_ptr)) || (is_y != (*is_y_ptr)) || (is_haploid != (*is_haploid_ptr))) {
    *is_x_ptr = is_x;
    *is_y_ptr = is_y;
    *is_haploid_ptr = is_haploid;
    uint32_t cur_hap_ct;
    if (is_x) {
      cur_hap_ct = 2 * founder_ct - founder_male_ct;
    } else if (is_y) {
      cur_hap_ct = founder_nonfemale_ct;
    } else {
      cur_hap_ct = founder_ct * (2 - is_haploid);
    }
    *cur_hap_ct_ptr = cur_hap_ct;
    *cur_hap_ctl_ptr = BitCtToWordCt(cur_hap_ct);
  }
}

typedef struct IndepPairphaseCtxStruct {
  const ChrInfo* cip;
  const uint32_t* subcontig_info;
  const uint32_t* subcontig_thread_assignments;
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const AlleleCode* maj_alleles;
  const double* all_allele_freqs;
  const uint32_t* variant_bps;
  const uintptr_t* preferred_variants;
  uint32_t* tvidx_end;
  uint32_t x_start;
  uint32_t x_len;
  uint32_t y_start;
  uint32_t y_len;
  uint32_t founder_ct;
  uint32_t founder_male_ct;
  uint32_t founder_nonfemale_ct;
  uint32_t prune_window_size;
  uint32_t window_maxl;
  double prune_ld_thresh;
  uint32_t window_incr;
  uint32_t cur_batch_size;

  uintptr_t** hap_then_nm_vecs;
  uintptr_t** occupied_window_slots;
  uintptr_t** cur_window_removed;
  double** cur_maj_freqs;
  uintptr_t** removed_variants_write;
  VariantHapAggs** vhaggs;
  uint32_t** winpos_to_slot_idx;
  uint32_t** tvidxs;
  uint32_t** first_unchecked_tvidx;
  uintptr_t** loader_hap_then_nm_vecs[2];
} IndepPairphaseCtx;

THREAD_FUNC_DECL IndepPairphaseThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  IndepPairphaseCtx* ctx = S_CAST(IndepPairphaseCtx*, arg->sharedp->context);

  const ChrInfo* cip = ctx->cip;
  const uint32_t* subcontig_info = ctx->subcontig_info;
  const uint32_t* subcontig_thread_assignments = ctx->subcontig_thread_assignments;
  const uintptr_t* variant_include = ctx->variant_include;
  const uintptr_t* preferred_variants = ctx->preferred_variants;
  const uint32_t x_start = ctx->x_start;
  const uint32_t x_len = ctx->x_len;
  const uint32_t y_start = ctx->y_start;
  const uint32_t y_len = ctx->y_len;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const AlleleCode* maj_alleles = ctx->maj_alleles;
  const double* all_allele_freqs = ctx->all_allele_freqs;
  const uint32_t* variant_bps = ctx->variant_bps;
  const uint32_t founder_ct = ctx->founder_ct;
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t founder_nonfemale_ct = ctx->founder_nonfemale_ct;
  const uint32_t max_hap_ctaw = BitCtToAlignedWordCt(founder_ct * 2);
  const uintptr_t max_hap_ctaw_x2 = max_hap_ctaw * 2;
  const uint32_t prune_window_size = ctx->prune_window_size;
  const uint32_t window_maxl = ctx->window_maxl;
  const double prune_ld_thresh = ctx->prune_ld_thresh;
  const uint32_t window_incr = ctx->window_incr;
  const uint32_t tvidx_end = ctx->tvidx_end[tidx];
  uintptr_t* hap_then_nm_vecs = ctx->hap_then_nm_vecs[tidx];
  uintptr_t* occupied_window_slots = ctx->occupied_window_slots[tidx];
  uintptr_t* cur_window_removed = ctx->cur_window_removed[tidx];
  uintptr_t* removed_variants_write = ctx->removed_variants_write[tidx];
  double* cur_maj_freqs = ctx->cur_maj_freqs[tidx];
  VariantHapAggs* vhaggs = ctx->vhaggs[tidx];
  uint32_t* winpos_to_slot_idx = ctx->winpos_to_slot_idx[tidx];
  uint32_t* tvidxs = ctx->tvidxs[tidx];
  uint32_t* first_unchecked_tvidx = ctx->first_unchecked_tvidx? ctx->first_unchecked_tvidx[tidx] : nullptr;

  uint32_t subcontig_end_tvidx = 0;
  uint32_t subcontig_idx = UINT32_MAX;  // deliberate overflow
  uint32_t window_start_tvidx = 0;
  uint32_t next_window_end_tvidx = 0;
  uint32_t write_slot_idx = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t is_haploid = IsSet(cip->haploid_mask, 0);
  uint32_t cur_window_size = 0;
  uint32_t winpos_split = 0;
  uint32_t tvidx_start = 0;
  uint32_t cur_hap_ct = founder_ct * (2 - is_haploid);
  uint32_t cur_hap_ctl = BitCtToWordCt(cur_hap_ct);
  uintptr_t variant_uidx_base = 0;
  uintptr_t variant_include_bits = variant_include[0];
  uint32_t variant_uidx_winstart = 0;
  uint32_t variant_uidx_winend = 0;
  uint32_t cur_allele_ct = 2;
  uint32_t parity = 0;
  do {
    const uint32_t cur_batch_size = ctx->cur_batch_size;
    const uint32_t tvidx_stop = MINV(tvidx_start + cur_batch_size, tvidx_end);
    // main loop has to be variant-, not window-, based due to how datasets too
    // large to fit in memory are handled: we may have to halt in the middle of
    // unpacking data for a window, waiting until the current I/O pass is
    // complete before proceeding
    const uintptr_t* loader_hap_then_nm_vecs = ctx->loader_hap_then_nm_vecs[parity][tidx];
    for (uint32_t cur_tvidx = tvidx_start; cur_tvidx < tvidx_stop; ) {
      if (cur_tvidx == subcontig_end_tvidx) {
        LdPruneNextSubcontig(variant_include, variant_bps, subcontig_info, subcontig_thread_assignments, prune_window_size, tidx, &subcontig_idx, &subcontig_end_tvidx, &next_window_end_tvidx, &variant_uidx_winstart, &variant_uidx_winend);
        IndepPairphaseUpdateSubcontig(cip, variant_uidx_winstart, x_start, x_len, y_start, y_len, founder_ct, founder_male_ct, founder_nonfemale_ct, &is_x, &is_y, &is_haploid, &cur_hap_ct, &cur_hap_ctl);
        BitIter1Start(variant_include, variant_uidx_winstart, &variant_uidx_base, &variant_include_bits);
        winpos_split = 0;
      }
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      write_slot_idx = AdvTo0Bit(occupied_window_slots, write_slot_idx);
      uintptr_t tvidx_offset = cur_tvidx - tvidx_start;
      {
        const uintptr_t* cur_loader_hap_then_nm_vecs = &(loader_hap_then_nm_vecs[tvidx_offset * max_hap_ctaw_x2]);
        uintptr_t* cur_hap_vec = &(hap_then_nm_vecs[write_slot_idx * max_hap_ctaw_x2]);
        memcpy(cur_hap_vec, cur_loader_hap_then_nm_vecs, max_hap_ctaw_x2 * sizeof(intptr_t));
        uintptr_t* cur_nm_vec = &(cur_hap_vec[max_hap_ctaw]);
        if (FillVhaggs(cur_hap_vec, cur_nm_vec, cur_hap_ctl, &(vhaggs[write_slot_idx]))) {
          SetBit(cur_window_size, cur_window_removed);
          SetBit(cur_tvidx, removed_variants_write);
        } else {
          tvidxs[write_slot_idx] = cur_tvidx;
          uintptr_t allele_idx_base;
          if (!allele_idx_offsets) {
            allele_idx_base = variant_uidx;
          } else {
            allele_idx_base = allele_idx_offsets[variant_uidx];
            cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_base;
            allele_idx_base -= variant_uidx;
          }
          cur_maj_freqs[write_slot_idx] = GetAlleleFreq(&(all_allele_freqs[allele_idx_base]), maj_alleles[variant_uidx], cur_allele_ct);
          if (preferred_variants && IsSet(preferred_variants, variant_uidx)) {
            cur_maj_freqs[write_slot_idx] -= 1.0;
          }
          if (first_unchecked_tvidx) {
            first_unchecked_tvidx[write_slot_idx] = cur_tvidx + 1;
          }
        }
      }
      SetBit(write_slot_idx, occupied_window_slots);
      winpos_to_slot_idx[cur_window_size++] = write_slot_idx;
      ++cur_tvidx;
      // are we at the end of a window?  if not, load more variant(s) before
      // proceeding.
      if (cur_tvidx != next_window_end_tvidx) {
        continue;
      }
      if (first_unchecked_tvidx) {
        // PLINK 1.x pruning order

        // possible for cur_window_size == 1, if all variants at the end of the
        // previous window were pruned
        uint32_t cur_removed_ct = PopcountWords(cur_window_removed, BitCtToWordCt(cur_window_size));
        uint32_t prev_removed_ct;
        do {
          prev_removed_ct = cur_removed_ct;
          for (uint32_t first_winpos = 0; ; ++first_winpos) {
            // can't use BitIter0 since we care about changes in this loop to
            // cur_window_removed
            first_winpos = AdvTo0Bit(cur_window_removed, first_winpos);
            // can assume empty trailing bit for cur_window_removed
            if (first_winpos == cur_window_size) {
              break;
            }
            const uint32_t first_slot_idx = winpos_to_slot_idx[first_winpos];
            const uint32_t cur_first_unchecked_tvidx = first_unchecked_tvidx[first_slot_idx];
            if (cur_first_unchecked_tvidx == cur_tvidx) {
              continue;
            }
            // safe to use BitIter0 for second_winpos, though
            uintptr_t second_winpos_base;
            uintptr_t cur_window_removed_inv_bits;
            BitIter0Start(cur_window_removed, first_winpos + 1, &second_winpos_base, &cur_window_removed_inv_bits);
            {
              uint32_t second_winpos;
              uint32_t second_slot_idx;
              do {
                second_winpos = BitIter0(cur_window_removed, &second_winpos_base, &cur_window_removed_inv_bits);
                if (second_winpos == cur_window_size) {
                  first_unchecked_tvidx[first_slot_idx] = cur_tvidx;
                  goto IndepPairphaseThread_next_first;
                }
                second_slot_idx = winpos_to_slot_idx[second_winpos];
              } while (tvidxs[second_slot_idx] < cur_first_unchecked_tvidx);
              const uintptr_t* first_hap_vec = &(hap_then_nm_vecs[first_slot_idx * max_hap_ctaw_x2]);
              const uint32_t first_nm_ct = vhaggs[first_slot_idx].nm_ct;
              const uint32_t first_sum = vhaggs[first_slot_idx].sum;
              while (1) {
                const uintptr_t* second_hap_vec = &(hap_then_nm_vecs[second_slot_idx * max_hap_ctaw_x2]);
                uint32_t cur_nm_ct = first_nm_ct;
                uint32_t cur_first_sum = first_sum;
                uint32_t second_sum;
                uint32_t cur_dotprod;
                ComputeIndepPairphaseR2Components(first_hap_vec, second_hap_vec, &(vhaggs[second_slot_idx]), max_hap_ctaw, cur_hap_ct, &cur_nm_ct, &cur_first_sum, &second_sum, &cur_dotprod);
                // these three values are actually cur_nm_ct times their
                // true values, but that cancels out
                const double cov12 = S_CAST(double, S_CAST(int64_t, cur_dotprod * S_CAST(uint64_t, cur_nm_ct) - S_CAST(uint64_t, cur_first_sum) * second_sum));
                const double variance1 = S_CAST(double, cur_first_sum * S_CAST(int64_t, cur_nm_ct - cur_first_sum));
                const double variance2 = S_CAST(double, second_sum * S_CAST(int64_t, cur_nm_ct - second_sum));
                // > instead of >=, so we don't prune from a pair of
                // variants with zero common observations
                if (cov12 * cov12 > prune_ld_thresh * variance1 * variance2) {
                  if (cur_maj_freqs[first_slot_idx] > cur_maj_freqs[second_slot_idx] * (1 + kSmallEpsilon)) {
                    SetBit(first_winpos, cur_window_removed);
                    SetBit(tvidxs[first_slot_idx], removed_variants_write);
                  } else {
                    SetBit(second_winpos, cur_window_removed);
                    SetBit(tvidxs[second_slot_idx], removed_variants_write);
                    const uint32_t next_start_winpos = BitIter0NoAdv(cur_window_removed, &second_winpos_base, &cur_window_removed_inv_bits);
                    if (next_start_winpos < cur_window_size) {
                      first_unchecked_tvidx[first_slot_idx] = tvidxs[winpos_to_slot_idx[next_start_winpos]];
                    } else {
                      first_unchecked_tvidx[first_slot_idx] = cur_tvidx;
                    }
                  }
                  break;
                }
                second_winpos = BitIter0(cur_window_removed, &second_winpos_base, &cur_window_removed_inv_bits);
                if (second_winpos == cur_window_size) {
                  first_unchecked_tvidx[first_slot_idx] = cur_tvidx;
                  break;
                }
                second_slot_idx = winpos_to_slot_idx[second_winpos];
              }  // while (1)
            }
          IndepPairphaseThread_next_first:
            ;
          }
          cur_removed_ct = PopcountWords(cur_window_removed, BitCtToWordCt(cur_window_size));
        } while (cur_removed_ct > prev_removed_ct);
      } else {
        // Within each window, scan in reverse order.  This way, we tend to
        // check the nearest new pairs first, and this should allow us to exit
        // early more often.
        const uint32_t second_winpos_stop = winpos_split? winpos_split : 1;
        for (uint32_t second_winpos = cur_window_size; second_winpos != second_winpos_stop; ) {
          --second_winpos;
          const uint32_t second_slot_idx = winpos_to_slot_idx[second_winpos];
          const uintptr_t* second_hap_vec = &(hap_then_nm_vecs[second_slot_idx * max_hap_ctaw_x2]);
          const uint32_t second_nm_ct = vhaggs[second_slot_idx].nm_ct;
          const int32_t second_sum = vhaggs[second_slot_idx].sum;
          for (uint32_t first_winpos = second_winpos; first_winpos; ) {
            --first_winpos;
            // possible todo: faster unset-bit reverse-iterator.  but probably
            // doesn't pay off here.
            if (IsSet(cur_window_removed, first_winpos)) {
              continue;
            }
            const uint32_t first_slot_idx = winpos_to_slot_idx[first_winpos];
            const uintptr_t* first_hap_vec = &(hap_then_nm_vecs[first_slot_idx * max_hap_ctaw_x2]);
            uint32_t cur_nm_ct = second_nm_ct;
            uint32_t cur_second_sum = second_sum;
            uint32_t first_sum;
            uint32_t cur_dotprod;
            ComputeIndepPairphaseR2Components(second_hap_vec, first_hap_vec, &(vhaggs[first_slot_idx]), max_hap_ctaw, cur_hap_ct, &cur_nm_ct, &cur_second_sum, &first_sum, &cur_dotprod);
            // these three values are actually cur_nm_ct times their
            // true values, but that cancels out
            const double cov12 = S_CAST(double, S_CAST(int64_t, cur_dotprod * S_CAST(uint64_t, cur_nm_ct) - S_CAST(uint64_t, first_sum) * cur_second_sum));
            const double variance1 = S_CAST(double, first_sum * S_CAST(int64_t, cur_nm_ct - first_sum));
            const double variance2 = S_CAST(double, cur_second_sum * S_CAST(int64_t, cur_nm_ct - cur_second_sum));
            // > instead of >=, so we don't prune from a pair of
            // variants with zero common observations
            // printf("cur_dotprod: %u  cur_nm_ct: %u  first_sum: %u  cur_second_sum: %u\n", cur_dotprod, cur_nm_ct, first_sum, cur_second_sum);
            // printf("vidx1: %u  vidx2: %u  r^2: %g\n", tvidxs[first_slot_idx], tvidxs[second_slot_idx], cov12 * cov12 / (variance1 * variance2));
            if (cov12 * cov12 > prune_ld_thresh * variance1 * variance2) {
              if (cur_maj_freqs[first_slot_idx] <= cur_maj_freqs[second_slot_idx] * (1 + kSmallEpsilon)) {
                SetBit(second_winpos, cur_window_removed);
                SetBit(tvidxs[second_slot_idx], removed_variants_write);
                break;
              }
              SetBit(first_winpos, cur_window_removed);
              SetBit(tvidxs[first_slot_idx], removed_variants_write);
            }
          }  // while (1)
        }
      }
      const uint32_t prev_window_size = cur_window_size;
      LdPruneNextWindow(variant_include, variant_bps, tvidxs, cur_window_removed, prune_window_size, window_incr, window_maxl, subcontig_end_tvidx, &cur_window_size, &window_start_tvidx, &variant_uidx_winstart, &next_window_end_tvidx, &variant_uidx_winend, occupied_window_slots, winpos_to_slot_idx);
      winpos_split = cur_window_size;
      // clear bits here since we set cur_window_removed bits during loading
      // process in monomorphic case
      ZeroWArr(BitCtToWordCt(prev_window_size), cur_window_removed);
      write_slot_idx = 0;
    }
    parity = 1 - parity;
    tvidx_start = tvidx_stop;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr IndepPairphase(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uint32_t* founder_info_cumulative_popcounts, const uintptr_t* founder_nonmale, const uintptr_t* founder_male, const uintptr_t* founder_nonfemale, const LdInfo* ldip, const uintptr_t* preferred_variants, const uint32_t* subcontig_info, const uint32_t* subcontig_thread_assignments, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t founder_nonfemale_ct, uint32_t subcontig_ct, uintptr_t window_max, uint32_t calc_thread_ct, uint32_t max_load, PgenReader* simple_pgrp, uintptr_t* removed_variants_collapsed) {
  ThreadGroup tg;
  PreinitThreads(&tg);
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t max_hap_ct = founder_ct * 2;
    const uint32_t hap_ctaw = BitCtToAlignedWordCt(max_hap_ct);
    const uintptr_t hap_ctaw_x2 = hap_ctaw * 2;

    // Per-thread allocations:
    // - tvidx_batch_size * hap_ctaw_x2 * sizeof(intptr_t) for loaded haplotype
    //     data (hap_then_nm_vecs)
    // - tvidx_batch_size * sizeof(double) for cur_maj_freqs
    // - if pos-based window, tvidx_batch_size * sizeof(int32_t)
    // - All of the above again, to allow loader thread to operate
    //     independently
    // - window_max * hap_ctaw_x2 * kBytesPerVec for current-window haplotype
    //     data
    // - max_loadl * sizeof(intptr_t) for removed-variant bitarray
    // - window_max * 2 * sizeof(int32_t) for main missing_ct, sum(x_i) array
    // - window_max * sizeof(int32_t) for indexes into genotype data bitarrays
    //     (for now, anyway)
    // - window_max * sizeof(int32_t) for live_indices (variant_idxs?)
    // - window_max * sizeof(int32_t) for start_arr (first uncompared
    //     variant_idx)
    const uint32_t founder_ctl2 = NypCtToWordCt(founder_ct);
    const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
    const uint32_t raw_sample_ctl2 = NypCtToWordCt(raw_sample_ct);
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    IndepPairphaseCtx ctx;
    uintptr_t* genovec;
    uintptr_t* phasepresent = nullptr; // spurious g++ 4.8 warning
    uintptr_t* phaseinfo = nullptr;
    uintptr_t* raw_genovec = nullptr;
    uintptr_t* raw_phasepresent = nullptr;
    uintptr_t* raw_phaseinfo = nullptr;
    uint32_t* thread_last_subcontig;
    uint32_t* thread_subcontig_start_tvidx;
    uint32_t* thread_last_tvidx;
    uint32_t* thread_last_uidx;
    if (unlikely(bigstack_alloc_w(founder_ctl2, &genovec) ||
                 bigstack_alloc_w(founder_ctl, &phasepresent) ||
                 bigstack_alloc_w(founder_ctl, &phaseinfo) ||
                 bigstack_alloc_w(raw_sample_ctl2, &raw_genovec) ||
                 bigstack_alloc_w(raw_sample_ctl, &raw_phasepresent) ||
                 bigstack_alloc_w(raw_sample_ctl, &raw_phaseinfo) ||
                 bigstack_calloc_u32(calc_thread_ct, &ctx.tvidx_end) ||
                 bigstack_calloc_u32(calc_thread_ct, &thread_last_subcontig) ||
                 bigstack_calloc_u32(calc_thread_ct, &thread_subcontig_start_tvidx) ||
                 bigstack_calloc_u32(calc_thread_ct, &thread_last_tvidx) ||
                 bigstack_calloc_u32(calc_thread_ct, &thread_last_uidx) ||
                 bigstack_alloc_wp(calc_thread_ct, &ctx.hap_then_nm_vecs) ||
                 bigstack_alloc_wp(calc_thread_ct, &ctx.occupied_window_slots) ||
                 bigstack_alloc_wp(calc_thread_ct, &ctx.cur_window_removed) ||
                 bigstack_alloc_dp(calc_thread_ct, &ctx.cur_maj_freqs) ||
                 bigstack_alloc_wp(calc_thread_ct, &ctx.removed_variants_write) ||
                 BIGSTACK_ALLOC_X(VariantHapAggs*, calc_thread_ct, &ctx.vhaggs) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.winpos_to_slot_idx) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.tvidxs) ||
                 bigstack_alloc_wp(calc_thread_ct, &(ctx.loader_hap_then_nm_vecs[0])) ||
                 bigstack_alloc_wp(calc_thread_ct, &(ctx.loader_hap_then_nm_vecs[1])))) {
      goto IndepPairphase_ret_NOMEM;
    }
    const uint32_t plink1_order = (ldip->prune_flags / kfLdPrunePlink1Order) & 1;
    if (plink1_order) {
      if (unlikely(bigstack_alloc_u32p(calc_thread_ct, &ctx.first_unchecked_tvidx))) {
        goto IndepPairphase_ret_NOMEM;
      }
    } else {
      ctx.first_unchecked_tvidx = nullptr;
    }
    for (uint32_t subcontig_idx = 0; subcontig_idx != subcontig_ct; ++subcontig_idx) {
      const uint32_t cur_thread_idx = subcontig_thread_assignments[subcontig_idx];
      ctx.tvidx_end[cur_thread_idx] += subcontig_info[3 * subcontig_idx];
    }
    const uint32_t window_maxl = BitCtToWordCt(window_max);
    const uint32_t max_loadl = BitCtToWordCt(max_load);
    const uintptr_t variant_vec_alloc = RoundUpPow2(window_max * hap_ctaw_x2 * sizeof(intptr_t), kCacheline);

    const uintptr_t occupied_window_slots_alloc = RoundUpPow2(window_maxl * sizeof(intptr_t), kCacheline);
    const uintptr_t cur_window_removed_alloc = RoundUpPow2((1 + window_max / kBitsPerWord) * sizeof(intptr_t), kCacheline);
    const uintptr_t cur_maj_freqs_alloc = RoundUpPow2(window_max * sizeof(double), kCacheline);
    const uintptr_t removed_variants_write_alloc = RoundUpPow2(max_loadl * sizeof(intptr_t), kCacheline);
    const uintptr_t vhaggs_alloc = RoundUpPow2(window_max * sizeof(VariantHapAggs), kCacheline);

    // (2 + plink1_order) of these
    const uintptr_t window_int32_alloc = RoundUpPow2(window_max * sizeof(int32_t), kCacheline);

    const uintptr_t thread_alloc_base = variant_vec_alloc + occupied_window_slots_alloc + cur_window_removed_alloc + cur_maj_freqs_alloc + removed_variants_write_alloc + vhaggs_alloc + (2 + plink1_order) * window_int32_alloc;

    // round down
    uintptr_t bigstack_avail_per_thread = RoundDownPow2(bigstack_left() / calc_thread_ct, kCacheline);
    const uintptr_t loader_single_variant_byte_ct = 2 * (hap_ctaw_x2 * sizeof(intptr_t));
    // may as well require capacity for >= 256 variants per thread per pass
    if (unlikely(bigstack_avail_per_thread <= thread_alloc_base + 256 * loader_single_variant_byte_ct)) {
      goto IndepPairphase_ret_NOMEM;
    }
    bigstack_avail_per_thread -= thread_alloc_base;
    uint32_t tvidx_batch_size = DivUp(max_load, 2);
    // tried a bunch of powers of two, this seems to be a good value
    if (tvidx_batch_size > 65536) {
      tvidx_batch_size = 65536;
    }
    // tvidx_batch_size = max_load;  // temporary debugging
    if (tvidx_batch_size * loader_single_variant_byte_ct > bigstack_avail_per_thread) {
      tvidx_batch_size = bigstack_avail_per_thread / loader_single_variant_byte_ct;
    }
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      ctx.hap_then_nm_vecs[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(variant_vec_alloc));
      ctx.occupied_window_slots[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(occupied_window_slots_alloc));
      ZeroWArr(window_maxl, ctx.occupied_window_slots[tidx]);
      ctx.cur_window_removed[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(cur_window_removed_alloc));
      ZeroWArr(1 + window_max / kBitsPerWord, ctx.cur_window_removed[tidx]);
      ctx.cur_maj_freqs[tidx] = S_CAST(double*, bigstack_alloc_raw(cur_maj_freqs_alloc));
      ctx.removed_variants_write[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(removed_variants_write_alloc));
      ZeroWArr(max_loadl, ctx.removed_variants_write[tidx]);
      ctx.vhaggs[tidx] = S_CAST(VariantHapAggs*, bigstack_alloc_raw(vhaggs_alloc));
      ctx.winpos_to_slot_idx[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(window_int32_alloc));
      ctx.tvidxs[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(window_int32_alloc));
      if (ctx.first_unchecked_tvidx) {
        ctx.first_unchecked_tvidx[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(window_int32_alloc));
      }
      ctx.loader_hap_then_nm_vecs[0][tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(tvidx_batch_size * hap_ctaw_x2 * sizeof(intptr_t)));
      ctx.loader_hap_then_nm_vecs[1][tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(tvidx_batch_size * hap_ctaw_x2 * sizeof(intptr_t)));
    }
    ctx.cip = cip;
    ctx.subcontig_info = subcontig_info;
    ctx.subcontig_thread_assignments = subcontig_thread_assignments;
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.maj_alleles = maj_alleles;
    ctx.all_allele_freqs = allele_freqs;
    ctx.variant_bps = variant_bps;
    ctx.preferred_variants = preferred_variants;
    ctx.founder_ct = founder_ct;
    ctx.founder_male_ct = founder_male_ct;
    ctx.founder_nonfemale_ct = founder_nonfemale_ct;
    ctx.prune_window_size = ldip->prune_window_size;
    ctx.window_maxl = window_maxl;
    ctx.prune_ld_thresh = ldip->prune_last_param * (1 + kSmallEpsilon);
    ctx.window_incr = ldip->prune_window_incr;
    ctx.cur_batch_size = tvidx_batch_size;
    const uint32_t all_haploid = IsSet(cip->haploid_mask, 0);
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    GetXymtStartAndEnd(cip, kChrOffsetX, &x_start, &x_end);
    GetXymtStartAndEnd(cip, kChrOffsetY, &y_start, &y_end);
    const uint32_t x_len = x_end - x_start;
    const uint32_t y_len = y_end - y_start;
    ctx.x_start = x_start;
    ctx.x_len = x_len;
    ctx.y_start = y_start;
    ctx.y_len = y_len;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto IndepPairphase_ret_NOMEM;
    }
    SetThreadFuncAndData(IndepPairphaseThread, &ctx, &tg);

    const uint32_t founder_nonmale_ct = founder_ct - founder_male_ct;
    const uint32_t founder_nonmale_ctl = BitCtToWordCt(founder_nonmale_ct);
    uintptr_t hap_vec_overflow = 0;
    uintptr_t nm_vec_overflow = 0;
    // Main workflow:
    // 1. Set n=0, load batch 0

    // 2. Spawn threads processing batch n
    // 3. Increment n by 1
    // 4. Load batch n unless eof
    // 5. Join threads
    // 6. Goto step 2 unless eof
    //
    // 7. Assemble final results with CopyBitarrRange()
    uint32_t parity = 0;
    uint32_t pct = 0;
    uint32_t next_print_tvidx_start = (max_load + 99) / 100;
    logprintf("--indep-pairphase (%u compute thread%s): ", calc_thread_ct, (calc_thread_ct == 1)? "" : "s");
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t cur_tvidx_start = 0; ; cur_tvidx_start += tvidx_batch_size) {
      if (!IsLastBlock(&tg)) {
        PgrSampleSubsetIndex pssi;
        PgrSetSampleSubsetIndex(founder_info_cumulative_popcounts, simple_pgrp, &pssi);
        uintptr_t** cur_loader_hap_then_nm_vecs = ctx.loader_hap_then_nm_vecs[parity];
        const uint32_t cur_tvidx_end = cur_tvidx_start + tvidx_batch_size;
        uint32_t is_x_or_y = 0;
        for (uint32_t subcontig_idx = 0; subcontig_idx != subcontig_ct; ++subcontig_idx) {
          const uint32_t cur_thread_idx = subcontig_thread_assignments[subcontig_idx];
          if (thread_last_subcontig[cur_thread_idx] > subcontig_idx) {
            continue;
          }
          uint32_t cur_tvidx = thread_last_tvidx[cur_thread_idx];
          if (cur_tvidx == cur_tvidx_end) {
            continue;
          }
          uint32_t subcontig_start_tvidx = thread_subcontig_start_tvidx[cur_thread_idx];
          uint32_t tvidx_end = subcontig_start_tvidx + subcontig_info[3 * subcontig_idx];
          if (tvidx_end > cur_tvidx_end) {
            tvidx_end = cur_tvidx_end;
            thread_last_subcontig[cur_thread_idx] = subcontig_idx;
          } else {
            thread_subcontig_start_tvidx[cur_thread_idx] = tvidx_end;
            thread_last_subcontig[cur_thread_idx] = subcontig_idx + 1;
          }
          uintptr_t tvidx_offset_end = tvidx_end - cur_tvidx_start;
          uint32_t variant_uidx;
          if (subcontig_start_tvidx == cur_tvidx) {
            variant_uidx = subcontig_info[3 * subcontig_idx + 2];
          } else {
            variant_uidx = thread_last_uidx[cur_thread_idx];
          }
          const uint32_t is_haploid = IsSet(cip->haploid_mask, GetVariantChr(cip, variant_uidx));
          uint32_t is_x = ((variant_uidx - x_start) < x_len);
          const uint32_t new_is_x_or_y = is_x || ((variant_uidx - y_start) < y_len);

          // due to nonempty subset requirement (removed?)
          is_x = is_x && founder_nonmale_ct;
          if (is_x_or_y != new_is_x_or_y) {
            is_x_or_y = new_is_x_or_y;
            if (is_x_or_y) {
              PgrClearSampleSubsetIndex(simple_pgrp, &pssi);
            } else {
              PgrSetSampleSubsetIndex(founder_info_cumulative_popcounts, simple_pgrp, &pssi);
            }
          }
          uintptr_t* cur_thread_loader_hap_then_nm_vecs = cur_loader_hap_then_nm_vecs[cur_thread_idx];
          uintptr_t variant_uidx_base;
          uintptr_t cur_bits;
          BitIter1Start(variant_include, variant_uidx, &variant_uidx_base, &cur_bits);
          --variant_uidx;
          for (uintptr_t tvidx_offset = cur_tvidx - cur_tvidx_start; tvidx_offset < tvidx_offset_end; ++tvidx_offset) {
            variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
            uintptr_t* cur_loader_hap_vec = &(cur_thread_loader_hap_then_nm_vecs[tvidx_offset * hap_ctaw_x2]);
            uintptr_t* cur_loader_nm_vec = &(cur_loader_hap_vec[hap_ctaw]);
            if (!is_x_or_y) {
              uint32_t phasepresent_ct;
              reterr = PgrGetInv1P(founder_info, pssi, founder_ct, variant_uidx, maj_alleles[variant_uidx], simple_pgrp, genovec, phasepresent, phaseinfo, &phasepresent_ct);
              if (unlikely(reterr)) {
                PgenErrPrintNV(reterr, variant_uidx);
                goto IndepPairphase_ret_1;
              }
              if (!is_haploid) {
                if (unlikely(HapsplitMustPhased(genovec, phasepresent, phaseinfo, founder_ct, phasepresent_ct, cur_loader_hap_vec, cur_loader_nm_vec))) {
                  logputs("\n");
                  logerrprintf("Error: --indep-pairphase: 0-based variant #%u is not fully phased.\n", variant_uidx);
                  goto IndepPairphase_ret_INCONSISTENT_INPUT;
                }
              } else {
                HapsplitHaploid(genovec, founder_ct, cur_loader_hap_vec, cur_loader_nm_vec);
              }
            } else {
              uint32_t phase_exists;
              reterr = PgrGetInv1P(nullptr, pssi, raw_sample_ct, variant_uidx, maj_alleles[variant_uidx], simple_pgrp, raw_genovec, raw_phasepresent, raw_phaseinfo, &phase_exists);
              if (unlikely(reterr)) {
                PgenErrPrintNV(reterr, variant_uidx);
                goto IndepPairphase_ret_1;
              }
              if (is_x) {
                if (founder_male_ct) {
                  CopyNyparrNonemptySubset(raw_genovec, founder_male, raw_sample_ct, founder_male_ct, genovec);
                  // quasi-bugfix (17 Oct 2023): forgot to remove benchmarking
                  // loop
                  HapsplitHaploid(genovec, founder_male_ct, cur_loader_hap_vec, cur_loader_nm_vec);
                }
                CopyNyparrNonemptySubset(raw_genovec, founder_nonmale, raw_sample_ct, founder_nonmale_ct, genovec);
                if (all_haploid) {
                  SetHetMissing(NypCtToWordCt(founder_nonmale_ct), genovec);
                }
                if (phase_exists) {
                  CopyBitarrSubset(raw_phasepresent, founder_nonmale, founder_nonmale_ct, phasepresent);
                  CopyBitarrSubset(raw_phaseinfo, founder_nonmale, founder_nonmale_ct, phaseinfo);
                  phase_exists = !AllWordsAreZero(phasepresent, founder_nonmale_ctl);
                }
                const uint32_t founder_male_fullword_ct = founder_male_ct / kBitsPerWord;
                const uint32_t founder_male_ct_rem = founder_male_ct % kBitsPerWord;
                if (founder_male_ct_rem) {
                  hap_vec_overflow = cur_loader_hap_vec[founder_male_fullword_ct];
                  nm_vec_overflow = cur_loader_nm_vec[founder_male_fullword_ct];
                }
                if (unlikely(HapsplitMustPhased(genovec, phasepresent, phaseinfo, founder_nonmale_ct, phase_exists, &(cur_loader_hap_vec[founder_male_fullword_ct]), &(cur_loader_nm_vec[founder_male_fullword_ct])))) {
                  logputs("\n");
                  logerrprintf("Error: --indep-pairphase: 0-based variant #%u is not fully phased.\n", variant_uidx);
                  goto IndepPairphase_ret_INCONSISTENT_INPUT;
                }
                if (founder_male_ct_rem) {
                  const uint32_t word_idx = founder_male_fullword_ct + ((founder_nonmale_ct * 2) / kBitsPerWord);
                  const uint32_t lshift_ct = (founder_nonmale_ct * 2) % kBitsPerWord;
                  cur_loader_hap_vec[word_idx] |= hap_vec_overflow << lshift_ct;
                  cur_loader_nm_vec[word_idx] |= nm_vec_overflow << lshift_ct;
                  if (lshift_ct + founder_male_ct_rem > kBitsPerWord) {
                    const uint32_t rshift_ct = kBitsPerWord - lshift_ct;
                    cur_loader_hap_vec[word_idx + 1] = hap_vec_overflow >> rshift_ct;
                    cur_loader_nm_vec[word_idx + 1] = nm_vec_overflow >> rshift_ct;
                  }
                }
              } else {
                // chrY
                if (founder_nonfemale_ct) {
                  CopyNyparrNonemptySubset(raw_genovec, founder_nonfemale, raw_sample_ct, founder_nonfemale_ct, genovec);
                  HapsplitHaploid(genovec, founder_nonfemale_ct, cur_loader_hap_vec, cur_loader_nm_vec);
                }
              }
            }
          }
          thread_last_tvidx[cur_thread_idx] = tvidx_end;
          thread_last_uidx[cur_thread_idx] = variant_uidx + 1;
        }
      }
      if (cur_tvidx_start) {
        JoinThreads(&tg);
        if (IsLastBlock(&tg)) {
          break;
        }
        if (cur_tvidx_start >= next_print_tvidx_start) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (cur_tvidx_start * 100LLU) / max_load;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_tvidx_start = (pct * S_CAST(uint64_t, max_load) + 99) / 100;
        }
      }
      if (cur_tvidx_start + tvidx_batch_size >= max_load) {
        DeclareLastThreadBlock(&tg);
      }
      if (unlikely(SpawnThreads(&tg))) {
        goto IndepPairphase_ret_THREAD_CREATE_FAIL;
      }
      parity = 1 - parity;
    }
    ZeroU32Arr(calc_thread_ct, thread_subcontig_start_tvidx);
    for (uint32_t subcontig_idx = 0; subcontig_idx != subcontig_ct; ++subcontig_idx) {
      const uint32_t cur_thread_idx = subcontig_thread_assignments[subcontig_idx];
      const uintptr_t* cur_removed_variants = ctx.removed_variants_write[cur_thread_idx];
      const uint32_t subcontig_len = subcontig_info[3 * subcontig_idx];
      const uint32_t subcontig_idx_start = subcontig_info[3 * subcontig_idx + 1];
      CopyBitarrRange(cur_removed_variants, thread_subcontig_start_tvidx[cur_thread_idx], subcontig_idx_start, subcontig_len, removed_variants_collapsed);
      thread_subcontig_start_tvidx[cur_thread_idx] += subcontig_len;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
  }
  while (0) {
  IndepPairphase_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  IndepPairphase_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  IndepPairphase_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 IndepPairphase_ret_1:
  CleanupThreads(&tg);
  // caller will free memory
  return reterr;
}

PglErr LdPruneSubcontigSplitAll(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, uint32_t prune_window_size, uint32_t* window_max_ptr, uint32_t** subcontig_info_ptr, uint32_t* subcontig_ct_ptr) {
  // variant_bps must be nullptr if window size is not bp-based
  // chr0 assumed to already be removed from variant_include.
  // this will skip over chromosomes/contigs with only 1 variant.
  const uint32_t chr_ct = cip->chr_ct;
  uint32_t* subcontig_info = R_CAST(uint32_t*, g_bigstack_base);
  uint32_t* subcontig_info_iter = subcontig_info;
  uint32_t* subcontig_info_limit = &(R_CAST(uint32_t*, g_bigstack_end)[-3]);
  uint32_t window_max = 0;
  uint32_t variant_idx = 0;
  if (variant_bps) {
    window_max = 1;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
      const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      const uint32_t initial_variant_uidx = AdvBoundedTo1Bit(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], chr_end);
      const uint32_t chr_variant_ct = PopcountBitRange(variant_include, initial_variant_uidx, chr_end);
      const uint32_t variant_idx_end = variant_idx + chr_variant_ct;
      if (chr_variant_ct > 1) {
        uintptr_t variant_uidx_base;
        uintptr_t variant_include_bits;
        BitIter1Start(variant_include, initial_variant_uidx + 1, &variant_uidx_base, &variant_include_bits);
        uint32_t subcontig_uidx_first = initial_variant_uidx;
        uint32_t subcontig_idx_first = variant_idx;
        uint32_t window_idx_first = variant_idx;
        uint32_t window_uidx_first = initial_variant_uidx;
        uint32_t window_pos_first = variant_bps[initial_variant_uidx];
        uint32_t prev_pos = window_pos_first;
        ++variant_idx;
        do {
          const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
          uint32_t variant_bp_thresh = variant_bps[variant_uidx];
          if (variant_bp_thresh < prune_window_size) {
            prev_pos = variant_bp_thresh;
            variant_bp_thresh = 0;
          } else {
            if (variant_bp_thresh - prune_window_size > prev_pos) {
              if (variant_idx > subcontig_idx_first + 1) {
                if (subcontig_info_iter > subcontig_info_limit) {
                  return kPglRetNomem;
                }
                *subcontig_info_iter++ = variant_idx - subcontig_idx_first;
                *subcontig_info_iter++ = subcontig_idx_first;
                *subcontig_info_iter++ = subcontig_uidx_first;
              }
              subcontig_uidx_first = variant_uidx;
              subcontig_idx_first = variant_idx;
            }
            prev_pos = variant_bp_thresh;
            variant_bp_thresh -= prune_window_size;
          }
          if (variant_bp_thresh > window_pos_first) {
            uintptr_t window_uidx_first_base;
            uintptr_t cur_bits;
            BitIter1Start(variant_include, window_uidx_first + 1, &window_uidx_first_base, &cur_bits);
            do {
              window_uidx_first = BitIter1(variant_include, &window_uidx_first_base, &cur_bits);
              window_pos_first = variant_bps[window_uidx_first];
              ++window_idx_first;
            } while (variant_bp_thresh > window_pos_first);
          } else if (variant_idx - window_idx_first == window_max) {
            ++window_max;
          }
        } while (++variant_idx < variant_idx_end);
        if (variant_idx > subcontig_idx_first + 1) {
          if (subcontig_info_iter > subcontig_info_limit) {
            return kPglRetNomem;
          }
          *subcontig_info_iter++ = variant_idx - subcontig_idx_first;
          *subcontig_info_iter++ = subcontig_idx_first;
          *subcontig_info_iter++ = subcontig_uidx_first;
        }
      }
      variant_idx = variant_idx_end;
    }
  } else {
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
      const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      const uint32_t first_variant_uidx = AdvBoundedTo1Bit(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], chr_end);
      const uint32_t chr_variant_ct = PopcountBitRange(variant_include, first_variant_uidx, chr_end);
      if (chr_variant_ct > 1) {
        if (subcontig_info_iter > subcontig_info_limit) {
          return kPglRetNomem;
        }
        *subcontig_info_iter++ = chr_variant_ct;
        *subcontig_info_iter++ = variant_idx;
        *subcontig_info_iter++ = first_variant_uidx;
        if (window_max < prune_window_size) {
          if (chr_variant_ct > window_max) {
            window_max = chr_variant_ct;
          }
        }
      }
      variant_idx += chr_variant_ct;
    }
    if (window_max > prune_window_size) {
      window_max = prune_window_size;
    }
  }
  *subcontig_ct_ptr = S_CAST(uintptr_t, subcontig_info_iter - subcontig_info) / 3;
  *subcontig_info_ptr = subcontig_info;
  BigstackFinalizeU32(subcontig_info, (*subcontig_ct_ptr) * 3);
  *window_max_ptr = window_max;
  return kPglRetSuccess;
}

// next several functions (including load_balance()) will probably move to
// plink2_common
void Minheap64ReplaceRoot(uint32_t heap_size, uint64_t new_root, uint64_t* minheap64_preroot) {
  uint32_t cur_pos = 1;
  while (1) {
    uint32_t child_pos = cur_pos * 2;
    if (child_pos >= heap_size) {
      if (child_pos == heap_size) {
        // special case: one child at end of heap
        const uint64_t child_val = minheap64_preroot[child_pos];
        if (new_root > child_val) {
          minheap64_preroot[cur_pos] = child_val;
          cur_pos = child_pos;
        }
      }
      break;
    }
    uint64_t min_child_val = minheap64_preroot[child_pos];
    const uint64_t child_val2 = minheap64_preroot[child_pos + 1];
    if (child_val2 < min_child_val) {
      min_child_val = child_val2;
      ++child_pos;
    }
    if (new_root <= min_child_val) {
      break;
    }
    minheap64_preroot[cur_pos] = min_child_val;
    cur_pos = child_pos;
  }
  minheap64_preroot[cur_pos] = new_root;
}

/*
void Minheap64DeleteRoot(uint64_t* minheap64_preroot, uint32_t* heap_size_ptr) {
  uint32_t heap_size = *heap_size_ptr;
  const uint64_t new_root = minheap64_preroot[heap_size];
  Minheap64ReplaceRoot(--heap_size, new_root, minheap64_preroot);
  *heap_size_ptr = heap_size;
}
*/

void Minheap64Insert(uint64_t new_entry, uint64_t* minheap64_preroot, uint32_t* heap_size_ptr) {
  // assumes minheap64_preroot[0] == 0
  const uint32_t heap_size = 1 + (*heap_size_ptr);
  *heap_size_ptr = heap_size;
  uint32_t cur_pos = heap_size;
  while (1) {
    const uint32_t parent_pos = cur_pos / 2;
    const uint64_t parent_val = minheap64_preroot[parent_pos];
    if (new_entry >= parent_val) {
      minheap64_preroot[cur_pos] = new_entry;
      return;
    }
    minheap64_preroot[cur_pos] = parent_val;
    cur_pos = parent_pos;
  }
}

// This is intended to split a relatively small number of contig-like regions
// between threads, but it shouldn't totally fall apart if there are millions
// of regions and hundreds of threads.
// Based on the Longest Processing Time algorithm, but with a few adjustments:
// * max(largest_weight, round_up(total_weight / thread_ct)) is noted, and the
//   first 8 * thread_ct thread assignments are based on best-fit to that
//   capacity.  The constant 8 is chosen to be enough to beat basic LPT's
//   4/3 - 1/{3m} approximation factor by a relevant margin, while keeping
//   runtime under control.  (In the event that there is no fit, the capacity
//   is increased.)
// * If any task assignments remain, we use LPT, but attempt to use a lower
//   number of threads; we only add another thread if we would otherwise have
//   to increase max_load.
PglErr LoadBalance(const uint32_t* task_weights, uint32_t task_ct, uint32_t* thread_ct_ptr, uint32_t* thread_assignments, uint32_t* max_load_ptr) {
  // max_load assumed to be initialized to zero
  assert(task_ct);
  const uint32_t orig_thread_ct = *thread_ct_ptr;
  if (orig_thread_ct == 1) {
    ZeroU32Arr(task_ct, thread_assignments);
    // replace this with an acc_uint32 call?
    uint32_t max_load = task_weights[0];
    for (uint32_t task_idx = 1; task_idx != task_ct; ++task_idx) {
      max_load += task_weights[task_idx];
    }
    *max_load_ptr = max_load;
    return kPglRetSuccess;
  }
  assert(task_ct >= orig_thread_ct);
  uint64_t* sorted_tagged_weights;
  uint64_t* minheap64_preroot;
  if (bigstack_alloc_u64(task_ct, &sorted_tagged_weights) ||
      bigstack_alloc_u64(orig_thread_ct + 2, &minheap64_preroot)) {
    return kPglRetNomem;
  }
  minheap64_preroot[0] = 0;
  uint64_t* minheap64 = &(minheap64_preroot[1]);
  uint32_t total_weight = 0;
  for (uintptr_t task_idx = 0; task_idx != task_ct; ++task_idx) {
    const uintptr_t cur_weight = task_weights[task_idx];
    total_weight += cur_weight;
    sorted_tagged_weights[task_idx] = (S_CAST(uint64_t, cur_weight) << 32) + task_idx;
  }
  uint64_t* sorted_tagged_weights_end = &(sorted_tagged_weights[task_ct]);
  // could try std::nth_element if this is ever a bottleneck
#ifdef __cplusplus
  std::sort(sorted_tagged_weights, sorted_tagged_weights_end, std::greater<uint64_t>());
#else
  qsort(sorted_tagged_weights, task_ct, sizeof(int64_t), u64cmp_decr);
#endif
  const uint64_t largest_tagged_weight = sorted_tagged_weights[0];
  uint32_t initial_max_load = largest_tagged_weight >> 32;
  uint32_t thread_ct = 1 + (total_weight - 1) / initial_max_load;
  if (thread_ct > orig_thread_ct) {
    thread_ct = orig_thread_ct;
    initial_max_load = 1 + (total_weight - 1) / orig_thread_ct;
  }

  for (uintptr_t thread_idx = 1; thread_idx != thread_ct; ++thread_idx) {
    minheap64[thread_idx - 1] = thread_ct - thread_idx;
  }
  minheap64[thread_ct - 1] = largest_tagged_weight & 0xffffffff00000000LLU;
  for (uint32_t thread_idx = thread_ct; thread_idx <= orig_thread_ct; ++thread_idx) {
    minheap64[thread_idx] = 0xffffffffffffffffLLU;
  }
  thread_assignments[S_CAST(uint32_t, largest_tagged_weight)] = 0;
  uint64_t max_load_shifted = (S_CAST(uint64_t, initial_max_load) << 32) | 0xffffffffLLU;
  uint64_t* best_fit_end = sorted_tagged_weights_end;
  if (task_ct > 8 * orig_thread_ct) {
    // stop best-fit here
    best_fit_end = &(sorted_tagged_weights[8 * orig_thread_ct]);
  }
  uint64_t* sorted_tagged_weights_iter = &(sorted_tagged_weights[1]);
  while (sorted_tagged_weights_iter != best_fit_end) {
    // maintain minheap64 as fully sorted list
    uint64_t cur_tagged_weight = *sorted_tagged_weights_iter++;
    const uint32_t task_idx = S_CAST(uint32_t, cur_tagged_weight);
    cur_tagged_weight &= 0xffffffff00000000LLU;
    const uintptr_t idxp1 = LowerBoundNonemptyU64(minheap64, thread_ct, max_load_shifted - cur_tagged_weight);
    if (idxp1) {
      uintptr_t idx = idxp1 - 1;
      const uint64_t new_entry = minheap64[idx] + cur_tagged_weight;
      for (; ; ++idx) {
        const uint64_t next_entry = minheap64[idx + 1];
        if (new_entry < next_entry) {
          break;
        }
        minheap64[idx] = next_entry;
      }
      thread_assignments[task_idx] = S_CAST(uint32_t, new_entry);
      minheap64[idx] = new_entry;
    } else if (thread_ct < orig_thread_ct) {
      const uint64_t new_entry = cur_tagged_weight + thread_ct;
      const uintptr_t insert_pt = LowerBoundNonemptyU64(minheap64, thread_ct, new_entry);
      for (uintptr_t thread_idx = thread_ct; thread_idx != insert_pt; --thread_idx) {
        minheap64[thread_idx] = minheap64[thread_idx - 1];
      }
      minheap64[insert_pt] = new_entry;
      thread_assignments[task_idx] = thread_ct++;
    } else {
      // move lowest entry to end of list, shift everything else down
      const uint64_t new_entry = minheap64[0] + cur_tagged_weight;
      for (uint32_t thread_idx = 1; thread_idx != thread_ct; ++thread_idx) {
        minheap64[thread_idx - 1] = minheap64[thread_idx];
      }
      minheap64[thread_ct - 1] = new_entry;
      max_load_shifted = new_entry | 0xffffffffLLU;
      thread_assignments[task_idx] = S_CAST(uint32_t, new_entry);
    }
  }
  if (best_fit_end != sorted_tagged_weights_end) {
    do {
      const uint64_t cur_heaproot = minheap64[0];
      uint64_t cur_tagged_weight = *sorted_tagged_weights_iter++;
      const uint32_t task_idx = S_CAST(uint32_t, cur_tagged_weight);
      uint32_t cur_thread = S_CAST(uint32_t, cur_heaproot);
      cur_tagged_weight &= 0xffffffff00000000LLU;
      uint64_t new_entry = cur_heaproot + cur_tagged_weight;
      if (new_entry > max_load_shifted) {
        if (thread_ct < orig_thread_ct) {
          thread_assignments[task_idx] = thread_ct;
          Minheap64Insert(cur_tagged_weight + thread_ct, minheap64_preroot, &thread_ct);
          continue;
        } else {
          max_load_shifted = new_entry | 0xffffffffLLU;
        }
      }
      thread_assignments[task_idx] = cur_thread;
      Minheap64ReplaceRoot(thread_ct, new_entry, minheap64_preroot);
    } while (sorted_tagged_weights_iter != sorted_tagged_weights_end);
  }
  BigstackReset(sorted_tagged_weights);
  *thread_ct_ptr = thread_ct;
  *max_load_ptr = max_load_shifted >> 32;
  return kPglRetSuccess;
}

PglErr LdPruneWrite(const uintptr_t* variant_include, const uintptr_t* removed_variants_collapsed, const char* const* variant_ids, uint32_t variant_ct, char* outname, char* outname_end) {
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    fputs("Writing...", stdout);
    fflush(stdout);
    snprintf(outname_end, kMaxOutfnameExtBlen, ".prune.in");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto LdPruneWrite_ret_OPEN_FAIL;
    }
    char* write_iter = g_textbuf;
    char* textbuf_flush = &(write_iter[kMaxMediumLine]);
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (IsSet(removed_variants_collapsed, variant_idx)) {
        continue;
      }
      write_iter = strcpya(write_iter, variant_ids[variant_uidx]);
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto LdPruneWrite_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto LdPruneWrite_ret_WRITE_FAIL;
    }

    snprintf(&(outname_end[7]), kMaxOutfnameExtBlen - 7, "out");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto LdPruneWrite_ret_OPEN_FAIL;
    }
    write_iter = g_textbuf;
    variant_uidx_base = 0;
    cur_bits = variant_include[0];
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (!IsSet(removed_variants_collapsed, variant_idx)) {
        continue;
      }
      write_iter = strcpya(write_iter, variant_ids[variant_uidx]);
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto LdPruneWrite_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto LdPruneWrite_ret_WRITE_FAIL;
    }
    *outname_end = '\0';
    putc_unlocked('\r', stdout);
    logprintfww("Variant lists written to %s.prune.in and %s.prune.out .\n", outname, outname);
  }
  while (0) {
  LdPruneWrite_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  LdPruneWrite_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  return reterr;
}

PglErr LdPrune(const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const LdInfo* ldip, const char* indep_preferred_fname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t nosex_ct, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  // common initialization between --indep-pairwise and --indep-pairphase
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t is_pairphase = (ldip->prune_flags / kfLdPrunePairphase) & 1;
    if (founder_ct < 2) {
      logerrprintfww("Error: --indep-pair%s requires at least two founders. (--make-founders may come in handy here.)\n", is_pairphase? "phase" : "wise");
      goto LdPrune_ret_INCONSISTENT_INPUT;
    }
    uint32_t skipped_variant_ct;
    const uintptr_t* variant_include = StripUnplacedK(orig_variant_include, cip, raw_variant_ct, &skipped_variant_ct);
    if (unlikely(variant_include == nullptr)) {
      goto LdPrune_ret_NOMEM;
    }
    if (skipped_variant_ct) {
      logprintf("--indep-pair%s: Ignoring %u chromosome 0 variant%s.\n", is_pairphase? "phase" : "wise", skipped_variant_ct, (skipped_variant_ct == 1)? "" : "s");
      variant_ct -= skipped_variant_ct;
    }

    if (!(ldip->prune_flags & kfLdPruneWindowBp)) {
      variant_bps = nullptr;
    }
    const uint32_t prune_window_size = ldip->prune_window_size;
    uint32_t* subcontig_info;
    uint32_t window_max;
    uint32_t subcontig_ct;
    if (LdPruneSubcontigSplitAll(variant_include, cip, variant_bps, prune_window_size, &window_max, &subcontig_info, &subcontig_ct)) {
      goto LdPrune_ret_NOMEM;
    }
    if (!subcontig_ct) {
      logerrprintf("Warning: Skipping --indep-pair%s since there are no pairs of variants to\nprocess.\n", is_pairphase? "phase" : "wise");
      goto LdPrune_ret_1;
    }

    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* preferred_variants = nullptr;
    uint32_t dup_found;
    if (!indep_preferred_fname) {
      reterr = CheckIdUniqueness(g_bigstack_base, g_bigstack_end, variant_include, variant_ids, variant_ct, max_thread_ct, &dup_found);
      if (unlikely(reterr)) {
        goto LdPrune_ret_1;
      }
    } else {
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, &preferred_variants))) {
        goto LdPrune_ret_NOMEM;
      }
      memcpy(preferred_variants, variant_include, raw_variant_ctl * sizeof(intptr_t));
      reterr = NondupIdLoad(g_bigstack_base, g_bigstack_end, variant_ids, indep_preferred_fname, raw_variant_ct, variant_ct, max_thread_ct, preferred_variants, &dup_found, g_logbuf);
      if (unlikely(reterr)) {
        if (g_logbuf[0]) {
          logerrputsb();
        }
        goto LdPrune_ret_1;
      }
    }
    if (unlikely(dup_found)) {
      logerrprintfww("Error: --indep-pair%s requires unique variant IDs. (--set-all-var-ids and/or --rm-dup may help.)\n", is_pairphase? "phase" : "wise");
      goto LdPrune_ret_INCONSISTENT_INPUT;
    }
    if (preferred_variants) {
      const uint32_t preferred_variant_ct = PopcountWords(preferred_variants, raw_variant_ctl);
      logprintf("--indep-preferred: %u variant%s loaded.\n", preferred_variant_ct, (preferred_variant_ct == 1)? "" : "s");
    }

    if (max_thread_ct > 2) {
      --max_thread_ct;
    }
    if (max_thread_ct > subcontig_ct) {
      max_thread_ct = subcontig_ct;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
    uint32_t* founder_info_cumulative_popcounts;
    // bugfix (25 Mar 2023, 25 Aug 2023): founder_nonmale/founder_male are NOT
    // supposed to be "collapsed".
    uintptr_t* founder_nonmale;
    uintptr_t* founder_male;
    uintptr_t* removed_variants_collapsed;
    uint32_t* subcontig_thread_assignments;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &founder_info_cumulative_popcounts) ||
                 bigstack_alloc_w(raw_sample_ctl, &founder_nonmale) ||
                 bigstack_alloc_w(raw_sample_ctl, &founder_male) ||
                 bigstack_calloc_w(variant_ctl, &removed_variants_collapsed) ||
                 bigstack_alloc_u32(subcontig_ct, &subcontig_thread_assignments))) {
      goto LdPrune_ret_NOMEM;
    }
    FillCumulativePopcounts(founder_info, raw_sample_ctl, founder_info_cumulative_popcounts);
    BitvecAndCopy(founder_info, sex_male, raw_sample_ctl, founder_male);
    BitvecInvmaskCopy(founder_info, sex_male, raw_sample_ctl, founder_nonmale);
    const uint32_t founder_male_ct = PopcountWords(founder_male, raw_sample_ctl);
    uintptr_t* founder_nonfemale = founder_male;
    uint32_t founder_nonfemale_ct = founder_male_ct;
    if (nosex_ct) {
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &founder_nonfemale))) {
        goto LdPrune_ret_NOMEM;
      }
      AlignedBitarrOrnotCopy(sex_male, sex_nm, raw_sample_ct, founder_nonfemale);
      BitvecAnd(founder_info, raw_sample_ctl, founder_nonfemale);
      founder_nonfemale_ct = PopcountWords(founder_nonfemale, raw_sample_ctl);
    }
    uint32_t* subcontig_weights;
    if (unlikely(bigstack_end_alloc_u32(subcontig_ct, &subcontig_weights))) {
      goto LdPrune_ret_NOMEM;
    }

    // initial window_max-based memory requirement estimate
    const uint32_t plink1_order = (ldip->prune_flags / kfLdPrunePlink1Order) & 1;
    if (is_pairphase) {
      const uint32_t max_hap_ct = founder_ct * 2;
      const uintptr_t hap_ctaw_x2 = 2 * BitCtToAlignedWordCt(max_hap_ct);
      // reserve ~1/2 of space for main variant data buffer,
      //   removed_variant_write
      // everything else:
      //   hap_vecs + nm_vecs: thread_ct * window_max * hap_ctaw_x2 * word
      //   occupied_window_slots: thread_ct * window_maxl * word
      //   cur_window_removed: thread_ct * (1 + window_max / kBitsPerWord) *
      //     word
      //   (ignore removed_variant_write)
      //   maj_freqs: thread_ct * window_max * 8
      //   vhaggs: thread_ct * window_max * VariantHapAggs
      //   winpos_to_slot_idx, tvidxs: window_max * 2 * int32
      //   first_unchecked_tvidx: window_max * int32 if plink1_order
      uintptr_t per_thread_alloc = RoundUpPow2(window_max * hap_ctaw_x2 * sizeof(intptr_t), kCacheline) + 2 * RoundUpPow2((1 + window_max / kBitsPerWord) * sizeof(intptr_t), kCacheline) + RoundUpPow2(window_max * sizeof(double), kCacheline) + RoundUpPow2(window_max * sizeof(VariantHapAggs), kCacheline) + (2 + plink1_order) * RoundUpPow2(window_max * sizeof(int32_t), kCacheline);
      uintptr_t bigstack_left2 = bigstack_left();
      if (per_thread_alloc * max_thread_ct > bigstack_left2) {
        if (unlikely(per_thread_alloc > bigstack_left2)) {
          goto LdPrune_ret_NOMEM;
        }
        max_thread_ct = bigstack_left2 / per_thread_alloc;
      }
    } else {
      const uintptr_t entire_variant_buf_word_ct = 2 * (BitCtToAlignedWordCt(founder_ct - founder_male_ct) + BitCtToAlignedWordCt(founder_male_ct));
      // reserve ~1/2 of space for main variant data buffer,
      //   removed_variant_write
      // everything else:
      //   genobufs: thread_ct * window_max * entire_variant_buf_word_ct * word
      //   occupied_window_slots: thread_ct * window_maxl * word
      //   cur_window_removed: thread_ct * (1 + window_max / kBitsPerWord) *
      //     word
      //   (ignore removed_variant_write)
      //   maj_freqs: thread_ct * window_max * 8
      //   vaggs, nonmale_vaggs: thread_ct * window_max * VariantAggs
      //   winpos_to_slot_idx, tvidxs: window_max * 2 * int32
      //   first_unchecked_tvidx: window_max * int32 if plink1_order
      uintptr_t per_thread_alloc = RoundUpPow2(window_max * entire_variant_buf_word_ct * sizeof(intptr_t), kCacheline) + 2 * RoundUpPow2((1 + window_max / kBitsPerWord) * sizeof(intptr_t), kCacheline) + RoundUpPow2(window_max * sizeof(double), kCacheline) + 2 * RoundUpPow2(window_max * sizeof(VariantAggs), kCacheline) + (2 + plink1_order) * RoundUpPow2(window_max * sizeof(int32_t), kCacheline);
      uintptr_t bigstack_left2 = bigstack_left();
      if (per_thread_alloc * max_thread_ct > bigstack_left2) {
        if (unlikely(per_thread_alloc > bigstack_left2)) {
          goto LdPrune_ret_NOMEM;
        }
        max_thread_ct = bigstack_left2 / per_thread_alloc;
      }
    }


    for (uint32_t subcontig_idx = 0; subcontig_idx != subcontig_ct; ++subcontig_idx) {
      // todo: adjust chrX weights upward, and chrY downward
      subcontig_weights[subcontig_idx] = subcontig_info[3 * subcontig_idx];
      // printf("%u %u %u\n", subcontig_info[3 * subcontig_idx], subcontig_info[3 * subcontig_idx + 1], subcontig_info[3 * subcontig_idx + 2]);
    }
    uint32_t max_load = 0;
    if (unlikely(LoadBalance(subcontig_weights, subcontig_ct, &max_thread_ct, subcontig_thread_assignments, &max_load))) {
      goto LdPrune_ret_NOMEM;
    }
    BigstackEndReset(bigstack_end_mark);

    if (is_pairphase) {
      reterr = IndepPairphase(variant_include, cip, variant_bps, allele_idx_offsets, maj_alleles, allele_freqs, founder_info, founder_info_cumulative_popcounts, founder_nonmale, founder_male, founder_nonfemale, ldip, preferred_variants, subcontig_info, subcontig_thread_assignments, raw_sample_ct, founder_ct, founder_male_ct, founder_nonfemale_ct, subcontig_ct, window_max, max_thread_ct, max_load, simple_pgrp, removed_variants_collapsed);
    } else {
      reterr = IndepPairwise(variant_include, cip, variant_bps, allele_idx_offsets, maj_alleles, allele_freqs, founder_info, founder_info_cumulative_popcounts, founder_nonmale, founder_male, founder_nonfemale, ldip, preferred_variants, subcontig_info, subcontig_thread_assignments, raw_sample_ct, founder_ct, founder_male_ct, founder_nonfemale_ct, subcontig_ct, window_max, max_thread_ct, max_load, simple_pgrp, removed_variants_collapsed);
    }
    if (unlikely(reterr)) {
      goto LdPrune_ret_1;
    }
    const uint32_t removed_ct = PopcountWords(removed_variants_collapsed, variant_ctl);
    logprintf("%u/%u variants removed.\n", removed_ct, variant_ct);
    reterr = LdPruneWrite(variant_include, removed_variants_collapsed, variant_ids, variant_ct, outname, outname_end);
    // if (unlikely(reterr)) {
    //   goto LdPrune_ret_1;
    // }
  }
  while (0) {
  LdPrune_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LdPrune_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 LdPrune_ret_1:
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}


// todo: see if this can also be usefully condensed into two bitarrays
#if defined(USE_SSE2) && !defined(USE_AVX2)
void GenoarrSplit12Nm(const uintptr_t* __restrict genoarr, uint32_t sample_ct, uintptr_t* __restrict one_bitarr, uintptr_t* __restrict two_bitarr, uintptr_t* __restrict nm_bitarr) {
  // ok if trailing bits of genoarr are not zeroed out
  // trailing bits of {one,two,nm}_bitarr are zeroed out
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t out_fullvec_ct = sample_ctl2 / (kWordsPerVec * 2);
  const VecW m1 = VCONST_W(kMask5555);
#  ifdef USE_SHUFFLE8
  const VecW swap12 = vecw_setr8(
      0, 1, 4, 5, 2, 3, 6, 7,
      8, 9, 12, 13, 10, 11, 14, 15);
#  else
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  for (uintptr_t vidx = 0; vidx != out_fullvec_ct; ++vidx) {
    const VecW vec_lo = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx]));
    const VecW vec_hi = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx + kWordsPerVec]));
    VecW inv0_lo = vecw_and_notfirst(vec_lo, m1);
    VecW inv0_hi = vecw_and_notfirst(vec_hi, m1);
    VecW inv1_lo = vecw_and_notfirst(vecw_srli(vec_lo, 1), m1);
    VecW inv1_hi = vecw_and_notfirst(vecw_srli(vec_hi, 1), m1);
#  ifdef USE_SHUFFLE8
    inv0_lo = (inv0_lo | vecw_srli(inv0_lo, 3)) & m4;
    inv0_hi = (inv0_hi | vecw_srli(inv0_hi, 3)) & m4;
    inv1_lo = (inv1_lo | vecw_srli(inv1_lo, 3)) & m4;
    inv1_hi = (inv1_hi | vecw_srli(inv1_hi, 3)) & m4;
    inv0_lo = vecw_shuffle8(swap12, inv0_lo);
    inv0_hi = vecw_shuffle8(swap12, inv0_hi);
    inv1_lo = vecw_shuffle8(swap12, inv1_lo);
    inv1_hi = vecw_shuffle8(swap12, inv1_hi);
#  else
    inv0_lo = (inv0_lo | vecw_srli(inv0_lo, 1)) & m2;
    inv0_hi = (inv0_hi | vecw_srli(inv0_hi, 1)) & m2;
    inv1_lo = (inv1_lo | vecw_srli(inv1_lo, 1)) & m2;
    inv1_hi = (inv1_hi | vecw_srli(inv1_hi, 1)) & m2;
    inv0_lo = (inv0_lo | vecw_srli(inv0_lo, 2)) & m4;
    inv0_hi = (inv0_hi | vecw_srli(inv0_hi, 2)) & m4;
    inv1_lo = (inv1_lo | vecw_srli(inv1_lo, 2)) & m4;
    inv1_hi = (inv1_hi | vecw_srli(inv1_hi, 2)) & m4;
#  endif
    inv0_lo = inv0_lo | vecw_srli(inv0_lo, 4);
    inv0_hi = inv0_hi | vecw_srli(inv0_hi, 4);
    inv1_lo = inv1_lo | vecw_srli(inv1_lo, 4);
    inv1_hi = inv1_hi | vecw_srli(inv1_hi, 4);
    const VecW inv0_packed = vecw_gather_even(inv0_lo, inv0_hi, m8);
    const VecW inv1_packed = vecw_gather_even(inv1_lo, inv1_hi, m8);
    const VecW one_packed = vecw_and_notfirst(inv0_packed, inv1_packed);
    const VecW two_packed = vecw_and_notfirst(inv1_packed, inv0_packed);
    const VecW nm_packed = inv0_packed | inv1_packed;
    vecw_storeu(&(one_bitarr[kWordsPerVec * vidx]), one_packed);
    vecw_storeu(&(two_bitarr[kWordsPerVec * vidx]), two_packed);
    vecw_storeu(&(nm_bitarr[kWordsPerVec * vidx]), nm_packed);
  }
  Halfword* one_bitarr_alias = R_CAST(Halfword*, one_bitarr);
  Halfword* two_bitarr_alias = R_CAST(Halfword*, two_bitarr);
  Halfword* nm_bitarr_alias = R_CAST(Halfword*, nm_bitarr);
  for (uint32_t widx = RoundDownPow2(sample_ctl2, kWordsPerVec * 2); widx != sample_ctl2; ++widx) {
    const uintptr_t cur_geno_word = genoarr[widx];
    const uint32_t low_halfword = PackWordToHalfwordMask5555(cur_geno_word);
    const uint32_t high_halfword = PackWordToHalfwordMaskAAAA(cur_geno_word);
    one_bitarr_alias[widx] = low_halfword & (~high_halfword);
    two_bitarr_alias[widx] = high_halfword & (~low_halfword);
    nm_bitarr_alias[widx] = ~(low_halfword & high_halfword);
  }
  const uint32_t sample_ct_rem = sample_ct % kBitsPerWord;
  if (sample_ct_rem) {
    const uint32_t last_widx = sample_ct / kBitsPerWord;
    const uintptr_t trailing_mask = (k1LU << sample_ct_rem) - 1;
    uintptr_t* __attribute__((may_alias)) one_bitarr_last = &(one_bitarr[last_widx]);
    uintptr_t* __attribute__((may_alias)) two_bitarr_last = &(two_bitarr[last_widx]);
    uintptr_t* __attribute__((may_alias)) nm_bitarr_last = &(nm_bitarr[last_widx]);
    *one_bitarr_last &= trailing_mask;
    *two_bitarr_last &= trailing_mask;
    *nm_bitarr_last &= trailing_mask;
  }
}
#else
void GenoarrSplit12Nm(const uintptr_t* __restrict genoarr, uint32_t sample_ct, uintptr_t* __restrict one_bitarr, uintptr_t* __restrict two_bitarr, uintptr_t* __restrict nm_bitarr) {
  // ok if trailing bits of genoarr are not zeroed out
  // trailing bits of {one,two,nm}_bitarr are zeroed out
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  Halfword* one_bitarr_alias = R_CAST(Halfword*, one_bitarr);
  Halfword* two_bitarr_alias = R_CAST(Halfword*, two_bitarr);
  Halfword* nm_bitarr_alias = R_CAST(Halfword*, nm_bitarr);
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uintptr_t cur_geno_word = genoarr[widx];
    const uint32_t low_halfword = PackWordToHalfwordMask5555(cur_geno_word);
    const uint32_t high_halfword = PackWordToHalfwordMaskAAAA(cur_geno_word);
    one_bitarr_alias[widx] = low_halfword & (~high_halfword);
    two_bitarr_alias[widx] = high_halfword & (~low_halfword);
    nm_bitarr_alias[widx] = ~(low_halfword & high_halfword);
  }

  const uint32_t sample_ct_rem = sample_ct % kBitsPerWord;
  if (sample_ct_rem) {
    const uint32_t last_widx = sample_ct / kBitsPerWord;
    const uintptr_t trailing_mask = (k1LU << sample_ct_rem) - 1;
    uintptr_t* __attribute__((may_alias)) one_bitarr_last = &(one_bitarr[last_widx]);
    uintptr_t* __attribute__((may_alias)) two_bitarr_last = &(two_bitarr[last_widx]);
    uintptr_t* __attribute__((may_alias)) nm_bitarr_last = &(nm_bitarr[last_widx]);
    *one_bitarr_last &= trailing_mask;
    *two_bitarr_last &= trailing_mask;
    *nm_bitarr_last &= trailing_mask;
  }
}
#endif

uint32_t GenoBitvecSumMain(const VecW* one_vvec, const VecW* two_vvec, uint32_t vec_ct) {
  // Analog of popcount_vecs.
  const VecW m0 = vecw_setzero();  // bugfix (15 Aug 2018)
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* one_vvec_iter = one_vvec;
  const VecW* two_vvec_iter = two_vvec;
  VecW prev_sad_result = vecw_setzero();
  VecW acc = vecw_setzero();
  uint32_t cur_incr = 15;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 15) {
      if (!vec_ct) {
        acc = acc + prev_sad_result;
        return HsumW(acc);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc = vecw_setzero();
    const VecW* one_vvec_stop = &(one_vvec_iter[cur_incr]);
    do {
      VecW one_count = *one_vvec_iter++;
      VecW two_count = *two_vvec_iter++;
      one_count = one_count - (vecw_srli(one_count, 1) & m1);
      two_count = two_count - (vecw_srli(two_count, 1) & m1);
      one_count = (one_count & m2) + (vecw_srli(one_count, 2) & m2);
      two_count = (two_count & m2) + (vecw_srli(two_count, 2) & m2);
      // one_count and two_count now contain 4-bit partial bitcounts, each in
      // the range 0..4.  finally enough room to compute
      //   2 * two_count + one_count
      // in parallel and add it to the accumulator.
      one_count = vecw_slli(two_count, 1) + one_count;
      inner_acc = inner_acc + (one_count & m4) + (vecw_srli(one_count, 4) & m4);
    } while (one_vvec_iter < one_vvec_stop);
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

uint32_t GenoBitvecSum(const uintptr_t* one_bitvec, const uintptr_t* two_bitvec, uint32_t word_ct) {
  //   popcount(one_bitvec) + 2 * popcount(two_bitvec)
  uint32_t tot = 0;
#ifdef __LP64__
  if (word_ct >= kWordsPerVec) {
#endif
    const uint32_t remainder = word_ct % kWordsPerVec;
    const uint32_t main_block_word_ct = word_ct - remainder;
    word_ct = remainder;
    tot = GenoBitvecSumMain(R_CAST(const VecW*, one_bitvec), R_CAST(const VecW*, two_bitvec), main_block_word_ct / kWordsPerVec);
#ifdef __LP64__
    one_bitvec = &(one_bitvec[main_block_word_ct]);
    two_bitvec = &(two_bitvec[main_block_word_ct]);
  }
  for (uint32_t trailing_word_idx = 0; trailing_word_idx != word_ct; ++trailing_word_idx) {
    tot += PopcountWord(one_bitvec[trailing_word_idx]) + 2 * PopcountWord(two_bitvec[trailing_word_idx]);
  }
#endif
  return tot;
}

uint32_t GenoBitvecSumSubsetMain(const VecW* subset_vvec, const VecW* one_vvec, const VecW* two_vvec, uint32_t vec_ct) {
  // Same as GenoBitvecSumMain(), just with an additional mask.
  const VecW m0 = vecw_setzero();
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* subset_vvec_iter = subset_vvec;
  const VecW* one_vvec_iter = one_vvec;
  const VecW* two_vvec_iter = two_vvec;
  VecW prev_sad_result = vecw_setzero();
  VecW acc = vecw_setzero();
  uint32_t cur_incr = 15;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 15) {
      if (!vec_ct) {
        acc = acc + prev_sad_result;
        return HsumW(acc);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc = vecw_setzero();
    const VecW* subset_vvec_stop = &(subset_vvec_iter[cur_incr]);
    do {
      VecW maskv = *subset_vvec_iter++;
      VecW one_count = (*one_vvec_iter++) & maskv;
      VecW two_count = (*two_vvec_iter++) & maskv;
      one_count = one_count - (vecw_srli(one_count, 1) & m1);
      two_count = two_count - (vecw_srli(two_count, 1) & m1);
      one_count = (one_count & m2) + (vecw_srli(one_count, 2) & m2);
      two_count = (two_count & m2) + (vecw_srli(two_count, 2) & m2);
      one_count = vecw_slli(two_count, 1) + one_count;
      inner_acc = inner_acc + (one_count & m4) + (vecw_srli(one_count, 4) & m4);
    } while (subset_vvec_iter < subset_vvec_stop);
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

uint32_t GenoBitvecSumSubset(const uintptr_t* subset_mask, const uintptr_t* one_bitvec, const uintptr_t* two_bitvec, uint32_t word_ct) {
  //   popcount(subset_mask & one_bitvec)
  // + 2 * popcount(subset_mask & two_bitvec)
  uint32_t tot = 0;
#ifdef __LP64__
  if (word_ct >= kWordsPerVec) {
#endif
    const uint32_t remainder = word_ct % kWordsPerVec;
    const uint32_t main_block_word_ct = word_ct - remainder;
    word_ct = remainder;
    tot = GenoBitvecSumSubsetMain(R_CAST(const VecW*, subset_mask), R_CAST(const VecW*, one_bitvec), R_CAST(const VecW*, two_bitvec), main_block_word_ct / kWordsPerVec);
#ifdef __LP64__
    subset_mask = &(subset_mask[main_block_word_ct]);
    one_bitvec = &(one_bitvec[main_block_word_ct]);
    two_bitvec = &(two_bitvec[main_block_word_ct]);
  }
  for (uint32_t trailing_word_idx = 0; trailing_word_idx != word_ct; ++trailing_word_idx) {
    const uintptr_t subset_word = subset_mask[trailing_word_idx];
    tot += PopcountWord(subset_word & one_bitvec[trailing_word_idx]) + 2 * PopcountWord(subset_word & two_bitvec[trailing_word_idx]);
  }
#endif
  return tot;
}

// phased-hardcall r^2 computation:
//   definitely-known part of dot product is
//     popcount((one_bitvec0 & two_bitvec1) | (two_bitvec0 & one_bitvec1))
//   + popcount(two_bitvec0 & two_bitvec1) * 2
//   + possible phased-het-het term
//   possibly-unknown part is
//     popcount(one_bitvec0 & one_bitvec1) - phased-het-het count
//   when nm_bitvec0 isn't all-ones, also necessary to compute
//     popcount(nm_bitvec0 & one_bitvec1)
//   + popcount(nm_bitvec0 & two_bitvec1) * 2
//   analogous statement is true for nm_bitvec1
//   if both are incomplete, also need popcount of intersection (compute this
//     first and skip rest of computation when zero).
//
// for possibly-unknown part, --ld reports all solutions when multiple
// solutions exist, everything else uses EM solution

void GenoBitvecPhasedDotprodMain(const VecW* one_vvec0, const VecW* two_vvec0, const VecW* one_vvec1, const VecW* two_vvec1, uint32_t vec_ct, uint32_t* __restrict known_dotprod_ptr, uint32_t* __restrict hethet_ct_ptr) {
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* one_vvec0_iter = one_vvec0;
  const VecW* two_vvec0_iter = two_vvec0;
  const VecW* one_vvec1_iter = one_vvec1;
  const VecW* two_vvec1_iter = two_vvec1;
  VecW acc_dotprod = vecw_setzero();
  VecW acc_hethet = vecw_setzero();
  uint32_t cur_incr = 15;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 15) {
      if (!vec_ct) {
        *known_dotprod_ptr = HsumW(acc_dotprod);
        *hethet_ct_ptr = HsumW(acc_hethet);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_dotprod = vecw_setzero();
    VecW inner_acc_hethet = vecw_setzero();
    const VecW* one_vvec0_stop = &(one_vvec0_iter[cur_incr]);
    do {
      VecW one_vword0 = *one_vvec0_iter++;
      VecW two_vword0 = *two_vvec0_iter++;
      VecW one_vword1 = *one_vvec1_iter++;
      VecW two_vword1 = *two_vvec1_iter++;

      VecW dotprod_1x_bits = (one_vword0 & two_vword1) | (one_vword1 & two_vword0);
      VecW dotprod_2x_bits = two_vword0 & two_vword1;
      VecW hethet_bits = one_vword0 & one_vword1;
      dotprod_1x_bits = dotprod_1x_bits - (vecw_srli(dotprod_1x_bits, 1) & m1);
      dotprod_2x_bits = dotprod_2x_bits - (vecw_srli(dotprod_2x_bits, 1) & m1);
      hethet_bits = hethet_bits - (vecw_srli(hethet_bits, 1) & m1);
      dotprod_1x_bits = (dotprod_1x_bits & m2) + (vecw_srli(dotprod_1x_bits, 2) & m2);
      dotprod_2x_bits = (dotprod_2x_bits & m2) + (vecw_srli(dotprod_2x_bits, 2) & m2);
      hethet_bits = (hethet_bits & m2) + (vecw_srli(hethet_bits, 2) & m2);

      // we now have 4-bit partial bitcounts in the range 0..4.  finally have
      // enough room to compute 2 * dotprod_2x_bits + dotprod_1x_bits.
      dotprod_1x_bits = vecw_slli(dotprod_2x_bits, 1) + dotprod_1x_bits;
      inner_acc_hethet = inner_acc_hethet + ((hethet_bits + vecw_srli(hethet_bits, 4)) & m4);
      inner_acc_dotprod = inner_acc_dotprod + (dotprod_1x_bits & m4) + (vecw_srli(dotprod_1x_bits, 4) & m4);
    } while (one_vvec0_iter < one_vvec0_stop);
    const VecW m0 = vecw_setzero();
    acc_hethet = acc_hethet + vecw_bytesum(inner_acc_hethet, m0);
    acc_dotprod = acc_dotprod + vecw_bytesum(inner_acc_dotprod, m0);
  }
}

void GenoBitvecPhasedDotprod(const uintptr_t* one_bitvec0, const uintptr_t* two_bitvec0, const uintptr_t* one_bitvec1, const uintptr_t* two_bitvec1, uint32_t word_ct, uint32_t* __restrict known_dotprod_ptr, uint32_t* __restrict hethet_ct_ptr) {
  // known_dotprod := popcount((one_bitvec0 & two_bitvec1) |
  //                           (two_bitvec0 & one_bitvec1)) +
  //                  2 * popcount(two_bitvec0 & two_bitvec1)
  // hethet_ct := popcount(one_bitvec0 & one_bitvec1)
  uint32_t known_dotprod = 0;
  uint32_t hethet_ct = 0;
#ifdef __LP64__
  if (word_ct >= kWordsPerVec) {
#endif
    const uint32_t remainder = word_ct % kWordsPerVec;
    const uint32_t main_block_word_ct = word_ct - remainder;
    word_ct = remainder;
    GenoBitvecPhasedDotprodMain(R_CAST(const VecW*, one_bitvec0), R_CAST(const VecW*, two_bitvec0), R_CAST(const VecW*, one_bitvec1), R_CAST(const VecW*, two_bitvec1), main_block_word_ct / kWordsPerVec, &known_dotprod, &hethet_ct);
#ifdef __LP64__
    one_bitvec0 = &(one_bitvec0[main_block_word_ct]);
    two_bitvec0 = &(two_bitvec0[main_block_word_ct]);
    one_bitvec1 = &(one_bitvec1[main_block_word_ct]);
    two_bitvec1 = &(two_bitvec1[main_block_word_ct]);
  }
  for (uint32_t trailing_word_idx = 0; trailing_word_idx != word_ct; ++trailing_word_idx) {
    const uintptr_t one_word0 = one_bitvec0[trailing_word_idx];
    const uintptr_t two_word0 = two_bitvec0[trailing_word_idx];
    const uintptr_t one_word1 = one_bitvec1[trailing_word_idx];
    const uintptr_t two_word1 = two_bitvec1[trailing_word_idx];
    known_dotprod += PopcountWord((one_word0 & two_word1) | (one_word1 & two_word0)) + 2 * PopcountWord(two_word0 & two_word1);
    hethet_ct += PopcountWord(one_word0 & one_word1);
  }
#endif
  *known_dotprod_ptr = known_dotprod;
  *hethet_ct_ptr = hethet_ct;
}

void GenoBitvecPhasedDotprodSubsetMain(const VecW* subset_vvec, const VecW* one_vvec0, const VecW* two_vvec0, const VecW* one_vvec1, const VecW* two_vvec1, uint32_t vec_ct, uint32_t* __restrict known_dotprod_ptr, uint32_t* __restrict hethet_ct_ptr) {
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* subset_vvec_iter = subset_vvec;
  const VecW* one_vvec0_iter = one_vvec0;
  const VecW* two_vvec0_iter = two_vvec0;
  const VecW* one_vvec1_iter = one_vvec1;
  const VecW* two_vvec1_iter = two_vvec1;
  VecW acc_dotprod = vecw_setzero();
  VecW acc_hethet = vecw_setzero();
  uint32_t cur_incr = 15;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 15) {
      if (!vec_ct) {
        *known_dotprod_ptr = HsumW(acc_dotprod);
        *hethet_ct_ptr = HsumW(acc_hethet);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_dotprod = vecw_setzero();
    VecW inner_acc_hethet = vecw_setzero();
    const VecW* one_vvec0_stop = &(one_vvec0_iter[cur_incr]);
    do {
      VecW subset_vword = *subset_vvec_iter++;
      VecW one_vword0 = (*one_vvec0_iter++) & subset_vword;
      VecW two_vword0 = (*two_vvec0_iter++) & subset_vword;
      VecW one_vword1 = *one_vvec1_iter++;
      VecW two_vword1 = *two_vvec1_iter++;

      VecW dotprod_1x_bits = (one_vword0 & two_vword1) | (one_vword1 & two_vword0);
      VecW dotprod_2x_bits = two_vword0 & two_vword1;
      VecW hethet_bits = one_vword0 & one_vword1;
      dotprod_1x_bits = dotprod_1x_bits - (vecw_srli(dotprod_1x_bits, 1) & m1);
      dotprod_2x_bits = dotprod_2x_bits - (vecw_srli(dotprod_2x_bits, 1) & m1);
      hethet_bits = hethet_bits - (vecw_srli(hethet_bits, 1) & m1);
      dotprod_1x_bits = (dotprod_1x_bits & m2) + (vecw_srli(dotprod_1x_bits, 2) & m2);
      dotprod_2x_bits = (dotprod_2x_bits & m2) + (vecw_srli(dotprod_2x_bits, 2) & m2);
      hethet_bits = (hethet_bits & m2) + (vecw_srli(hethet_bits, 2) & m2);

      // we now have 4-bit partial bitcounts in the range 0..4.  finally have
      // enough room to compute 2 * dotprod_2x_bits + dotprod_1x_bits.
      dotprod_1x_bits = vecw_slli(dotprod_2x_bits, 1) + dotprod_1x_bits;
      inner_acc_hethet = inner_acc_hethet + ((hethet_bits + vecw_srli(hethet_bits, 4)) & m4);
      inner_acc_dotprod = inner_acc_dotprod + (dotprod_1x_bits & m4) + (vecw_srli(dotprod_1x_bits, 4) & m4);
    } while (one_vvec0_iter < one_vvec0_stop);
    const VecW m0 = vecw_setzero();
    acc_hethet = acc_hethet + vecw_bytesum(inner_acc_hethet, m0);
    acc_dotprod = acc_dotprod + vecw_bytesum(inner_acc_dotprod, m0);
  }
}

void GenoBitvecPhasedDotprodSubset(const uintptr_t* subset_mask, const uintptr_t* one_bitvec0, const uintptr_t* two_bitvec0, const uintptr_t* one_bitvec1, const uintptr_t* two_bitvec1, uint32_t word_ct, uint32_t* __restrict known_dotprod_ptr, uint32_t* __restrict hethet_ct_ptr) {
  uint32_t known_dotprod = 0;
  uint32_t hethet_ct = 0;
#ifdef __LP64__
  if (word_ct >= kWordsPerVec) {
#endif
    const uint32_t remainder = word_ct % kWordsPerVec;
    const uint32_t main_block_word_ct = word_ct - remainder;
    word_ct = remainder;
    GenoBitvecPhasedDotprodSubsetMain(R_CAST(const VecW*, subset_mask), R_CAST(const VecW*, one_bitvec0), R_CAST(const VecW*, two_bitvec0), R_CAST(const VecW*, one_bitvec1), R_CAST(const VecW*, two_bitvec1), main_block_word_ct / kWordsPerVec, &known_dotprod, &hethet_ct);
#ifdef __LP64__
    subset_mask = &(subset_mask[main_block_word_ct]);
    one_bitvec0 = &(one_bitvec0[main_block_word_ct]);
    two_bitvec0 = &(two_bitvec0[main_block_word_ct]);
    one_bitvec1 = &(one_bitvec1[main_block_word_ct]);
    two_bitvec1 = &(two_bitvec1[main_block_word_ct]);
  }
  for (uint32_t trailing_word_idx = 0; trailing_word_idx != word_ct; ++trailing_word_idx) {
    const uintptr_t subset_word = subset_mask[trailing_word_idx];
    const uintptr_t one_word0 = one_bitvec0[trailing_word_idx] & subset_word;
    const uintptr_t two_word0 = two_bitvec0[trailing_word_idx] & subset_word;
    const uintptr_t one_word1 = one_bitvec1[trailing_word_idx];
    const uintptr_t two_word1 = two_bitvec1[trailing_word_idx];
    known_dotprod += PopcountWord((one_word0 & two_word1) | (one_word1 & two_word0)) + 2 * PopcountWord(two_word0 & two_word1);
    hethet_ct += PopcountWord(one_word0 & one_word1);
  }
#endif
  *known_dotprod_ptr = known_dotprod;
  *hethet_ct_ptr = hethet_ct;
}

// nmaj_cts[] must be initialized to correct values for
// no-missing-values-in-other-variant case.
uint32_t HardcallPhasedR2Stats(const uintptr_t* one_bitvec0, const uintptr_t* two_bitvec0, const uintptr_t* nm_bitvec0, const uintptr_t* one_bitvec1, const uintptr_t* two_bitvec1, const uintptr_t* nm_bitvec1, uint32_t sample_ct, uint32_t nm_ct0, uint32_t nm_ct1, uint32_t* __restrict nmaj_cts, uint32_t* __restrict known_dotprod_ptr, uint32_t* __restrict hethet_ct_ptr) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  uint32_t nm_intersection_ct;
  if ((nm_ct0 != sample_ct) && (nm_ct1 != sample_ct)) {
    nm_intersection_ct = PopcountWordsIntersect(nm_bitvec0, nm_bitvec1, sample_ctl);
    if (!nm_intersection_ct) {
      nmaj_cts[0] = 0;
      nmaj_cts[1] = 0;
      *known_dotprod_ptr = 0;
      *hethet_ct_ptr = 0;
      return 0;
    }
  } else {
    nm_intersection_ct = MINV(nm_ct0, nm_ct1);
  }
  if (nm_ct0 != nm_intersection_ct) {
    nmaj_cts[0] = GenoBitvecSumSubset(nm_bitvec1, one_bitvec0, two_bitvec0, sample_ctl);
  }
  if (nm_ct1 != nm_intersection_ct) {
    nmaj_cts[1] = GenoBitvecSumSubset(nm_bitvec0, one_bitvec1, two_bitvec1, sample_ctl);
  }
  GenoBitvecPhasedDotprod(one_bitvec0, two_bitvec0, one_bitvec1, two_bitvec1, sample_ctl, known_dotprod_ptr, hethet_ct_ptr);
  return nm_intersection_ct;
}

void HardcallPhasedR2RefineMain(const VecW* phasepresent0_vvec, const VecW* phaseinfo0_vvec, const VecW* phasepresent1_vvec, const VecW* phaseinfo1_vvec, uint32_t vec_ct, uint32_t* __restrict hethet_decr_ptr, uint32_t* __restrict not_dotprod_ptr) {
  // vec_ct must be a multiple of 3
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* phasepresent0_vvec_iter = phasepresent0_vvec;
  const VecW* phaseinfo0_vvec_iter = phaseinfo0_vvec;
  const VecW* phasepresent1_vvec_iter = phasepresent1_vvec;
  const VecW* phaseinfo1_vvec_iter = phaseinfo1_vvec;
  VecW acc_hethet_decr = vecw_setzero();
  VecW acc_not_dotprod = vecw_setzero();  // like not_hotdog, but more useful
  uint32_t cur_incr = 30;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 30) {
      if (!vec_ct) {
        *hethet_decr_ptr = HsumW(acc_hethet_decr);
        *not_dotprod_ptr = HsumW(acc_not_dotprod);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_hethet_decr = vecw_setzero();
    VecW inner_acc_not_dotprod = vecw_setzero();
    const VecW* phasepresent0_vvec_stop = &(phasepresent0_vvec_iter[cur_incr]);
    do {
      // todo: benchmark against simpler one-vec-at-a-time loop
      VecW mask1 = (*phasepresent0_vvec_iter++) & (*phasepresent1_vvec_iter++);
      VecW mask2 = (*phasepresent0_vvec_iter++) & (*phasepresent1_vvec_iter++);
      VecW mask_half1 = (*phasepresent0_vvec_iter++) & (*phasepresent1_vvec_iter++);
      VecW mask_half2 = vecw_srli(mask_half1, 1) & m1;
      mask_half1 = mask_half1 & m1;

      VecW not_dotprod_count1 = (*phaseinfo0_vvec_iter++) ^ (*phaseinfo1_vvec_iter++);
      VecW not_dotprod_count2 = (*phaseinfo0_vvec_iter++) ^ (*phaseinfo1_vvec_iter++);
      VecW not_dotprod_half1 = (*phaseinfo0_vvec_iter++) ^ (*phaseinfo1_vvec_iter++);
      // bugfix (4 Nov 2017): incorrectly had mask_half1 here
      VecW not_dotprod_half2 = vecw_srli(not_dotprod_half1, 1) & mask_half2;
      not_dotprod_count1 = not_dotprod_count1 & mask1;
      not_dotprod_count2 = not_dotprod_count2 & mask2;
      not_dotprod_half1 = not_dotprod_half1 & mask_half1;

      mask1 = mask1 - (vecw_srli(mask1, 1) & m1);
      mask2 = mask2 - (vecw_srli(mask2, 1) & m1);
      not_dotprod_count1 = not_dotprod_count1 - (vecw_srli(not_dotprod_count1, 1) & m1);
      not_dotprod_count2 = not_dotprod_count2 - (vecw_srli(not_dotprod_count2, 1) & m1);
      mask1 = mask1 + mask_half1;
      mask2 = mask2 + mask_half2;
      not_dotprod_count1 = not_dotprod_count1 + not_dotprod_half1;
      not_dotprod_count2 = not_dotprod_count2 + not_dotprod_half2;

      mask1 = (mask1 & m2) + (vecw_srli(mask1, 2) & m2);
      not_dotprod_count1 = (not_dotprod_count1 & m2) + (vecw_srli(not_dotprod_count1, 2) & m2);
      mask1 = mask1 + (mask2 & m2) + (vecw_srli(mask2, 2) & m2);
      not_dotprod_count1 = not_dotprod_count1 + (not_dotprod_count2 & m2) + (vecw_srli(not_dotprod_count2, 2) & m2);

      inner_acc_hethet_decr = inner_acc_hethet_decr + (mask1 & m4) + (vecw_srli(mask1, 4) & m4);
      inner_acc_not_dotprod = inner_acc_not_dotprod + (not_dotprod_count1 & m4) + (vecw_srli(not_dotprod_count1, 4) & m4);
    } while (phasepresent0_vvec_iter < phasepresent0_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_hethet_decr = acc_hethet_decr + vecw_bytesum(inner_acc_hethet_decr, m0);
    acc_not_dotprod = acc_not_dotprod + vecw_bytesum(inner_acc_not_dotprod, m0);
  }
}

// only needs to be called when hethet_ct > 0, phasepresent0_ct > 0, and
// phasepresent1_ct > 0.
void HardcallPhasedR2Refine(const uintptr_t* phasepresent0, const uintptr_t* phaseinfo0, const uintptr_t* phasepresent1, const uintptr_t* phaseinfo1, uint32_t word_ct, uint32_t* __restrict known_dotprod_ptr, uint32_t* __restrict unknown_hethet_ct_ptr) {
  // unknown_hethet_ct -= popcount(phasepresent0 & phasepresent1)
  // known_dotprod_ptr += popcount(phasepresent0 & phasepresent1 &
  //                               (~(phaseinfo0 ^ phaseinfo1)))
  uint32_t hethet_decr = 0;
  uint32_t not_dotprod = 0;
  if (word_ct >= 3 * kWordsPerVec) {
    const uint32_t remainder = word_ct % (3 * kWordsPerVec);
    const uint32_t main_block_word_ct = word_ct - remainder;
    word_ct = remainder;
    HardcallPhasedR2RefineMain(R_CAST(const VecW*, phasepresent0), R_CAST(const VecW*, phaseinfo0), R_CAST(const VecW*, phasepresent1), R_CAST(const VecW*, phaseinfo1), main_block_word_ct / kWordsPerVec, &hethet_decr, &not_dotprod);
    phasepresent0 = &(phasepresent0[main_block_word_ct]);
    phaseinfo0 = &(phaseinfo0[main_block_word_ct]);
    phasepresent1 = &(phasepresent1[main_block_word_ct]);
    phaseinfo1 = &(phaseinfo1[main_block_word_ct]);
  }
  for (uint32_t trailing_word_idx = 0; trailing_word_idx != word_ct; ++trailing_word_idx) {
    const uintptr_t mask = phasepresent0[trailing_word_idx] & phasepresent1[trailing_word_idx];
    const uintptr_t xor_word = phaseinfo0[trailing_word_idx] ^ phaseinfo1[trailing_word_idx];
    hethet_decr += PopcountWord(mask);
    not_dotprod += PopcountWord(mask & xor_word);
  }
  *known_dotprod_ptr += hethet_decr - not_dotprod;
  *unknown_hethet_ct_ptr -= hethet_decr;
}

void HardcallPhasedR2RefineSubsetMain(const VecW* subset_vvec, const VecW* phasepresent0_vvec, const VecW* phaseinfo0_vvec, const VecW* phasepresent1_vvec, const VecW* phaseinfo1_vvec, uint32_t vec_ct, uint32_t* __restrict hethet_decr_ptr, uint32_t* __restrict not_dotprod_ptr) {
  // vec_ct must be a multiple of 3
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* subset_vvec_iter = subset_vvec;
  const VecW* phasepresent0_vvec_iter = phasepresent0_vvec;
  const VecW* phaseinfo0_vvec_iter = phaseinfo0_vvec;
  const VecW* phasepresent1_vvec_iter = phasepresent1_vvec;
  const VecW* phaseinfo1_vvec_iter = phaseinfo1_vvec;
  VecW acc_hethet_decr = vecw_setzero();
  VecW acc_not_dotprod = vecw_setzero();
  uint32_t cur_incr = 30;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 30) {
      if (!vec_ct) {
        *hethet_decr_ptr = HsumW(acc_hethet_decr);
        *not_dotprod_ptr = HsumW(acc_not_dotprod);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_hethet_decr = vecw_setzero();
    VecW inner_acc_not_dotprod = vecw_setzero();
    const VecW* phasepresent0_vvec_stop = &(phasepresent0_vvec_iter[cur_incr]);
    do {
      VecW mask1 = (*phasepresent0_vvec_iter++) & (*phasepresent1_vvec_iter++) & (*subset_vvec_iter++);
      VecW mask2 = (*phasepresent0_vvec_iter++) & (*phasepresent1_vvec_iter++) & (*subset_vvec_iter++);
      VecW mask_half1 = (*phasepresent0_vvec_iter++) & (*phasepresent1_vvec_iter++) & (*subset_vvec_iter++);
      VecW mask_half2 = vecw_srli(mask_half1, 1) & m1;
      mask_half1 = mask_half1 & m1;

      VecW not_dotprod_count1 = (*phaseinfo0_vvec_iter++) ^ (*phaseinfo1_vvec_iter++);
      VecW not_dotprod_count2 = (*phaseinfo0_vvec_iter++) ^ (*phaseinfo1_vvec_iter++);
      VecW not_dotprod_half1 = (*phaseinfo0_vvec_iter++) ^ (*phaseinfo1_vvec_iter++);
      VecW not_dotprod_half2 = vecw_srli(not_dotprod_half1, 1) & mask_half2;
      not_dotprod_count1 = not_dotprod_count1 & mask1;
      not_dotprod_count2 = not_dotprod_count2 & mask2;
      not_dotprod_half1 = not_dotprod_half1 & mask_half1;

      mask1 = mask1 - (vecw_srli(mask1, 1) & m1);
      mask2 = mask2 - (vecw_srli(mask2, 1) & m1);
      not_dotprod_count1 = not_dotprod_count1 - (vecw_srli(not_dotprod_count1, 1) & m1);
      not_dotprod_count2 = not_dotprod_count2 - (vecw_srli(not_dotprod_count2, 1) & m1);
      mask1 = mask1 + mask_half1;
      mask2 = mask2 + mask_half2;
      not_dotprod_count1 = not_dotprod_count1 + not_dotprod_half1;
      not_dotprod_count2 = not_dotprod_count2 + not_dotprod_half2;

      mask1 = (mask1 & m2) + (vecw_srli(mask1, 2) & m2);
      not_dotprod_count1 = (not_dotprod_count1 & m2) + (vecw_srli(not_dotprod_count1, 2) & m2);
      mask1 = mask1 + (mask2 & m2) + (vecw_srli(mask2, 2) & m2);
      not_dotprod_count1 = not_dotprod_count1 + (not_dotprod_count2 & m2) + (vecw_srli(not_dotprod_count2, 2) & m2);

      inner_acc_hethet_decr = inner_acc_hethet_decr + (mask1 & m4) + (vecw_srli(mask1, 4) & m4);
      inner_acc_not_dotprod = inner_acc_not_dotprod + (not_dotprod_count1 & m4) + (vecw_srli(not_dotprod_count1, 4) & m4);
    } while (phasepresent0_vvec_iter < phasepresent0_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_hethet_decr = acc_hethet_decr + vecw_bytesum(inner_acc_hethet_decr, m0);
    acc_not_dotprod = acc_not_dotprod + vecw_bytesum(inner_acc_not_dotprod, m0);
  }
}

void HardcallPhasedR2RefineSubset(const uintptr_t* subset_mask, const uintptr_t* phasepresent0, const uintptr_t* phaseinfo0, const uintptr_t* phasepresent1, const uintptr_t* phaseinfo1, uint32_t word_ct, uint32_t* __restrict known_dotprod_ptr, uint32_t* __restrict unknown_hethet_ct_ptr) {
  uint32_t hethet_decr = 0;
  uint32_t not_dotprod = 0;
  if (word_ct >= 3 * kWordsPerVec) {
    const uint32_t remainder = word_ct % (3 * kWordsPerVec);
    const uint32_t main_block_word_ct = word_ct - remainder;
    word_ct = remainder;
    HardcallPhasedR2RefineSubsetMain(R_CAST(const VecW*, subset_mask), R_CAST(const VecW*, phasepresent0), R_CAST(const VecW*, phaseinfo0), R_CAST(const VecW*, phasepresent1), R_CAST(const VecW*, phaseinfo1), main_block_word_ct / kWordsPerVec, &hethet_decr, &not_dotprod);
    subset_mask = &(subset_mask[main_block_word_ct]);
    phasepresent0 = &(phasepresent0[main_block_word_ct]);
    phaseinfo0 = &(phaseinfo0[main_block_word_ct]);
    phasepresent1 = &(phasepresent1[main_block_word_ct]);
    phaseinfo1 = &(phaseinfo1[main_block_word_ct]);
  }
  for (uint32_t trailing_word_idx = 0; trailing_word_idx != word_ct; ++trailing_word_idx) {
    const uintptr_t mask = phasepresent0[trailing_word_idx] & phasepresent1[trailing_word_idx] & subset_mask[trailing_word_idx];
    const uintptr_t xor_word = phaseinfo0[trailing_word_idx] ^ phaseinfo1[trailing_word_idx];
    hethet_decr += PopcountWord(mask);
    not_dotprod += PopcountWord(mask & xor_word);
  }
  *known_dotprod_ptr += hethet_decr - not_dotprod;
  *unknown_hethet_ct_ptr -= hethet_decr;
}

// (Phased-)dosage r^2:
// 1. This must be defined such that, as the phased-dosages approach integers,
//    r^2 converges to the corresponding phased-hardcall value.
// 2. We also want r^2 between a variant and itself to be 1 (unless we have no
//    data at all).  Unfortunately, the formula used before 26 Oct 2023 did not
//    have this property.
//
// Suppose one unphased sample has dosage(var0)=0.2 and dosage(var1)=0.2.
// We treat this as P(var0=0/0)=0.8, P(var0=0/1)=0.2,
//                  P(var1=0/0)=0.8, P(var1=0/1)=0.2.
// Previously, we treated the two variants as independent, acting as if
// P(var0=0/0, var1=0/1) = 0.8 * 0.2, etc.; but that is incompatible with (2),
// and generally contrary to the whole point of a LD estimate.
// Revised computation has unknown-hethet frequency equal to the
// min(P(var0=0/1), P(var1=0/1)) upper bound, and known-dotprod equal to the
// corresponding max(0, P(var0=1/1) + P(var1=1/1) - 1) lower bound.

static_assert(sizeof(Dosage) == 2, "plink2_ld dosage-handling routines must be updated.");
#ifdef __LP64__
#  ifdef USE_AVX2
void FillDosageHet(const Dosage* dosage_vec, uint32_t dosagev_ct, Dosage* dosage_het) {
  const __m256i* dosage_vvec_iter = R_CAST(const __m256i*, dosage_vec);
#    if defined(__clang__) && ((!defined(__cplusplus)) || (__cplusplus < 201103L))
  const __m256i all_n32768 = _mm256_set1_epi16(0x8000);
  const __m256i all_n16384 = _mm256_set1_epi16(0xc000);
#    else
  const __m256i all_n32768 = _mm256_set1_epi64x(-0x7fff7fff7fff8000LL);
  const __m256i all_n16384 = _mm256_set1_epi64x(-0x3fff3fff3fff4000LL);
#    endif
  const __m256i all0 = _mm256_setzero_si256();
  const __m256i all1 = _mm256_cmpeq_epi16(all0, all0);
  // 0-16384: leave unchanged
  // 16385-32768: subtract from 32768
  // 65535: set to 0

  // subtract from 0, cmp_epi16 to produce mask, add 32768
  __m256i* dosage_het_iter = R_CAST(__m256i*, dosage_het);
  for (uint32_t vec_idx = 0; vec_idx != dosagev_ct; ++vec_idx) {
    __m256i dosagev = *dosage_vvec_iter++;

    __m256i cur_mask = _mm256_cmpeq_epi16(dosagev, all1);
    dosagev = _mm256_andnot_si256(cur_mask, dosagev);  // 65535 -> 0

    // xor with -32768 is same as subtracting it
    __m256i dosagev_opp = _mm256_xor_si256(dosagev, all_n32768);
    // anything > -16384 after this subtraction was originally >16384.
    // calling the original value x, we want to flip the sign of (x - 32768)
    cur_mask = _mm256_cmpgt_epi16(dosagev_opp, all_n16384);
    dosagev_opp = _mm256_and_si256(cur_mask, dosagev_opp);

    // has the <= 16384 values
    dosagev = _mm256_andnot_si256(cur_mask, dosagev);

    dosagev_opp = _mm256_sub_epi16(all0, dosagev_opp);
    *dosage_het_iter++ = _mm256_add_epi16(dosagev, dosagev_opp);
  }
}

uint64_t DenseDosageSum(const Dosage* dosage_vec, uint32_t vec_ct) {
  // end of dosage_vec assumed to be missing-padded (0-padded also ok)
  const __m256i* dosage_vvec_iter = R_CAST(const __m256i*, dosage_vec);
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
  const __m256i all1 = _mm256_cmpeq_epi16(m16, m16);
  uint64_t sum = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m256i sumv = _mm256_setzero_si256();
    const __m256i* dosage_vvec_stop;
    // individual values in [0..32768]
    // 32768 * 8191 * 16 dosages per __m256i = just under 2^32
    if (vecs_left < 8191) {
      if (!vecs_left) {
        return sum;
      }
      dosage_vvec_stop = &(dosage_vvec_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec_stop = &(dosage_vvec_iter[8191]);
      vecs_left -= 8191;
    }
    do {
      __m256i dosagev = *dosage_vvec_iter++;
      __m256i invmask = _mm256_cmpeq_epi16(dosagev, all1);
      dosagev = _mm256_andnot_si256(invmask, dosagev);

      dosagev = _mm256_add_epi64(_mm256_and_si256(dosagev, m16), _mm256_and_si256(_mm256_srli_epi64(dosagev, 16), m16));
      sumv = _mm256_add_epi64(sumv, dosagev);
    } while (dosage_vvec_iter < dosage_vvec_stop);
    UniVec acc;
    acc.vw = R_CAST(VecW, sumv);
    sum += UniVecHsum32(acc);
  }
}

uint64_t DenseDosageSumSubset(const Dosage* dosage_vec, const Dosage* dosage_mask_vec, uint32_t vec_ct) {
  // end of dosage_vec assumed to be missing-padded (0-padded also ok)
  const __m256i* dosage_vvec_iter = R_CAST(const __m256i*, dosage_vec);
  const __m256i* dosage_mask_vvec_iter = R_CAST(const __m256i*, dosage_mask_vec);
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
  const __m256i all1 = _mm256_cmpeq_epi16(m16, m16);
  uint64_t sum = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m256i sumv = _mm256_setzero_si256();
    const __m256i* dosage_vvec_stop;
    if (vecs_left < 8191) {
      if (!vecs_left) {
        return sum;
      }
      dosage_vvec_stop = &(dosage_vvec_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec_stop = &(dosage_vvec_iter[8191]);
      vecs_left -= 8191;
    }
    do {
      __m256i invmask = *dosage_mask_vvec_iter++;
      __m256i dosagev = *dosage_vvec_iter++;
      invmask = _mm256_cmpeq_epi16(invmask, all1);
      invmask = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosagev, all1));
      dosagev = _mm256_andnot_si256(invmask, dosagev);

      dosagev = _mm256_add_epi64(_mm256_and_si256(dosagev, m16), _mm256_and_si256(_mm256_srli_epi64(dosagev, 16), m16));
      sumv = _mm256_add_epi64(sumv, dosagev);
    } while (dosage_vvec_iter < dosage_vvec_stop);
    UniVec acc;
    acc.vw = R_CAST(VecW, sumv);
    sum += UniVecHsum32(acc);
  }
}

// 65535 treated as missing
uint64_t DosageUnsignedDotprod(const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const __m256i* dosage_vvec0_iter = R_CAST(const __m256i*, dosage_vec0);
  const __m256i* dosage_vvec1_iter = R_CAST(const __m256i*, dosage_vec1);
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
  const __m256i all1 = _mm256_cmpeq_epi16(m16, m16);
  uint64_t dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m256i dotprod_lo = _mm256_setzero_si256();
    __m256i dotprod_hi = _mm256_setzero_si256();
    const __m256i* dosage_vvec0_stop;
    // Products in [0..2^30]; dotprod_lo part is in 0..65535, dotprod_hi part
    // is in 0..16384.
    // 65535 * 4096 * 16 dosages per __m256i = just under 2^32
    if (vecs_left < 4096) {
      if (!vecs_left) {
        return dotprod;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[4096]);
      vecs_left -= 4096;
    }
    do {
      __m256i dosage0 = *dosage_vvec0_iter++;
      __m256i dosage1 = *dosage_vvec1_iter++;
      __m256i invmask = _mm256_cmpeq_epi16(dosage0, all1);
      invmask = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosage1, all1));
      dosage0 = _mm256_andnot_si256(invmask, dosage0);
      dosage1 = _mm256_andnot_si256(invmask, dosage1);

      // todo: try rewriting this loop without 256-bit multiplication, to avoid
      // ~15% universal clock frequency slowdown (128-bit?  sparse?)
      __m256i lo16 = _mm256_mullo_epi16(dosage0, dosage1);
      __m256i hi16 = _mm256_mulhi_epu16(dosage0, dosage1);
      lo16 = _mm256_add_epi64(_mm256_and_si256(lo16, m16), _mm256_and_si256(_mm256_srli_epi64(lo16, 16), m16));
      hi16 = _mm256_and_si256(_mm256_add_epi64(hi16, _mm256_srli_epi64(hi16, 16)), m16);
      dotprod_lo = _mm256_add_epi64(dotprod_lo, lo16);
      dotprod_hi = _mm256_add_epi64(dotprod_hi, hi16);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec acc_lo;
    UniVec acc_hi;
    acc_lo.vw = R_CAST(VecW, dotprod_lo);
    acc_hi.vw = R_CAST(VecW, dotprod_hi);
    dotprod += UniVecHsum32(acc_lo) + 65536 * UniVecHsum32(acc_hi);
  }
}

uint64_t DosageUnsignedDotprodSubset(const Dosage* dosage_mask_vec, const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const __m256i* dosage_mask_iter = R_CAST(const __m256i*, dosage_mask_vec);
  const __m256i* dosage_vvec0_iter = R_CAST(const __m256i*, dosage_vec0);
  const __m256i* dosage_vvec1_iter = R_CAST(const __m256i*, dosage_vec1);
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
  const __m256i all1 = _mm256_cmpeq_epi16(m16, m16);
  uint64_t dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m256i dotprod_lo = _mm256_setzero_si256();
    __m256i dotprod_hi = _mm256_setzero_si256();
    const __m256i* dosage_vvec0_stop;
    if (vecs_left < 4096) {
      if (!vecs_left) {
        return dotprod;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[4096]);
      vecs_left -= 4096;
    }
    do {
      __m256i cur_mask = *dosage_mask_iter++;
      __m256i dosage0 = *dosage_vvec0_iter++;
      __m256i dosage1 = *dosage_vvec1_iter++;
      __m256i invmask = _mm256_cmpeq_epi16(dosage0, all1);
      invmask = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosage1, all1));
      invmask = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(cur_mask, all1));
      dosage0 = _mm256_andnot_si256(invmask, dosage0);
      dosage1 = _mm256_andnot_si256(invmask, dosage1);

      // todo: try rewriting this loop without 256-bit multiplication, to avoid
      // ~15% universal clock frequency slowdown (128-bit?  sparse?)
      __m256i lo16 = _mm256_mullo_epi16(dosage0, dosage1);
      __m256i hi16 = _mm256_mulhi_epu16(dosage0, dosage1);
      lo16 = _mm256_add_epi64(_mm256_and_si256(lo16, m16), _mm256_and_si256(_mm256_srli_epi64(lo16, 16), m16));
      hi16 = _mm256_and_si256(_mm256_add_epi64(hi16, _mm256_srli_epi64(hi16, 16)), m16);
      dotprod_lo = _mm256_add_epi64(dotprod_lo, lo16);
      dotprod_hi = _mm256_add_epi64(dotprod_hi, hi16);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec acc_lo;
    UniVec acc_hi;
    acc_lo.vw = R_CAST(VecW, dotprod_lo);
    acc_hi.vw = R_CAST(VecW, dotprod_hi);
    dotprod += UniVecHsum32(acc_lo) + 65536 * UniVecHsum32(acc_hi);
  }
}

void DosageUnphasedDotprodComponents(const Dosage* dosage_vec0, const Dosage* dosage_vec1, const Dosage* dosage_het0, const Dosage* dosage_het1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const __m256i* dosage_vvec0_iter = R_CAST(const __m256i*, dosage_vec0);
  const __m256i* dosage_vvec1_iter = R_CAST(const __m256i*, dosage_vec1);
  const __m256i* dosage_het0_iter = R_CAST(const __m256i*, dosage_het0);
  const __m256i* dosage_het1_iter = R_CAST(const __m256i*, dosage_het1);
  const __m256i all_32767 = _mm256_set1_epi16(0x7fff);
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
  const __m256i all0 = _mm256_setzero_si256();
  const __m256i all1 = _mm256_cmpeq_epi16(all0, all0);
  uint64_t uhethet_dosage = 0;
  uint64_t known_dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m256i uhethet_sumv = _mm256_setzero_si256();
    __m256i dotprod_sumv = _mm256_setzero_si256();
    const __m256i* dosage_vvec0_stop;
    // individual dotprod_incr values in [0..32768]
    // 32768 * 8191 * 16 dosages per __m256i = just under 2^32
    if (vecs_left < 8191) {
      if (!vecs_left) {
        *known_dotprod_dosagep = known_dotprod * 2;
        *uhethet_dosagep = uhethet_dosage * 2;
        return;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[8191]);
      vecs_left -= 8191;
    }
    do {
      // We wish to compute max(0, dosage0 + dosage1 - 32768), where dosage0
      // and dosage1 are uint16s in {0..32768, 65535}, where 65535 corresponds
      // to a missing value.
      // An annoying property of this computation is that
      // (dosage0 + dosage1 - 32768) ranges from -32768 to 32768, which just
      // barely exceeds the range of a (u)int16.
      // We work around this by observing that dosage0==0 can be treated as if
      // it were a missing value: it is impossible for dosage0 + dosage1 -
      // 32768 to exceed 0 if dosage0 is 0.  With that case made ignorable,
      // (dosage0 + dosage1 - 32769) is then in the -32768..32767 int16 range,
      // so we compute 1 + max(-1, dosage0 + dosage1 - 32769) and then apply
      // our augmented mask.
      const __m256i dosage0_plus_32767 = _mm256_add_epi16(*dosage_vvec0_iter++, all_32767);
      const __m256i dosage1 = *dosage_vvec1_iter++;
      const __m256i het0 = *dosage_het0_iter++;
      const __m256i het1 = *dosage_het1_iter++;
      // Conveniently, dosage0 + 32767 > 0 in int16 space iff dosage0 is
      // missing or zero.
      const __m256i invmask = _mm256_or_si256(_mm256_cmpgt_epi16(dosage0_plus_32767, all0), _mm256_cmpeq_epi16(dosage1, all1));
      // No need to mask this, het0 or het1 is already equal to 0 when there's
      // a missing value.
      __m256i uhethet_incr = _mm256_min_epi16(het0, het1);
      // Adding 32767 and subtracting 32769 are equivalent in vector int16
      // space.  (Though they're not equivalent when working with single
      // int16_ts!)
      const __m256i dosagesum_m32769 = _mm256_add_epi16(dosage0_plus_32767, dosage1);
      const __m256i unmasked_dotprod_incr = _mm256_sub_epi16(_mm256_max_epi16(dosagesum_m32769, all1), all1);
      __m256i dotprod_incr = _mm256_andnot_si256(invmask, unmasked_dotprod_incr);

      uhethet_incr = _mm256_and_si256(_mm256_add_epi64(uhethet_incr, _mm256_srli_epi64(uhethet_incr, 16)), m16);
      dotprod_incr = _mm256_add_epi64(_mm256_and_si256(dotprod_incr, m16), _mm256_and_si256(_mm256_srli_epi64(dotprod_incr, 16), m16));
      uhethet_sumv = _mm256_add_epi64(uhethet_sumv, uhethet_incr);
      dotprod_sumv = _mm256_add_epi64(dotprod_sumv, dotprod_incr);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec uhethet_acc;
    UniVec dotprod_acc;
    uhethet_acc.vw = R_CAST(VecW, uhethet_sumv);
    dotprod_acc.vw = R_CAST(VecW, dotprod_sumv);
    uhethet_dosage += UniVecHsum32(uhethet_acc);
    known_dotprod += UniVecHsum32(dotprod_acc);
  }
}

void DosageUnphasedDotprodComponentsSubset(const Dosage* subset_invmask, const Dosage* dosage_vec0, const Dosage* dosage_vec1, const Dosage* dosage_het0, const Dosage* dosage_het1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const __m256i* invmask_iter = R_CAST(const __m256i*, subset_invmask);
  const __m256i* dosage_vvec0_iter = R_CAST(const __m256i*, dosage_vec0);
  const __m256i* dosage_vvec1_iter = R_CAST(const __m256i*, dosage_vec1);
  const __m256i* dosage_het0_iter = R_CAST(const __m256i*, dosage_het0);
  const __m256i* dosage_het1_iter = R_CAST(const __m256i*, dosage_het1);
  const __m256i all_32767 = _mm256_set1_epi16(0x7fff);
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
  const __m256i all0 = _mm256_setzero_si256();
  const __m256i all1 = _mm256_cmpeq_epi16(all0, all0);
  uint64_t uhethet_dosage = 0;
  uint64_t known_dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m256i uhethet_sumv = _mm256_setzero_si256();
    __m256i dotprod_sumv = _mm256_setzero_si256();
    const __m256i* dosage_vvec0_stop;
    if (vecs_left < 8191) {
      if (!vecs_left) {
        *known_dotprod_dosagep = known_dotprod * 2;
        *uhethet_dosagep = uhethet_dosage * 2;
        return;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[8191]);
      vecs_left -= 8191;
    }
    do {
      const __m256i raw_invmask = *invmask_iter++;
      const __m256i dosage0_plus_32767 = _mm256_add_epi16(*dosage_vvec0_iter++, all_32767);
      const __m256i dosage1 = *dosage_vvec1_iter++;
      const __m256i het0 = *dosage_het0_iter++;
      const __m256i het1 = *dosage_het1_iter++;
      __m256i uhethet_incr = _mm256_andnot_si256(raw_invmask, _mm256_min_epi16(het0, het1));
      __m256i invmask = _mm256_or_si256(raw_invmask, _mm256_cmpgt_epi16(dosage0_plus_32767, all0));
      invmask = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosage1, all1));
      const __m256i dosagesum_m32769 = _mm256_add_epi16(dosage0_plus_32767, dosage1);
      const __m256i unmasked_dotprod_incr = _mm256_sub_epi16(_mm256_max_epi16(dosagesum_m32769, all1), all1);
      __m256i dotprod_incr = _mm256_andnot_si256(invmask, unmasked_dotprod_incr);

      uhethet_incr = _mm256_and_si256(_mm256_add_epi64(uhethet_incr, _mm256_srli_epi64(uhethet_incr, 16)), m16);
      dotprod_incr = _mm256_add_epi64(_mm256_and_si256(dotprod_incr, m16), _mm256_and_si256(_mm256_srli_epi64(dotprod_incr, 16), m16));
      uhethet_sumv = _mm256_add_epi64(uhethet_sumv, uhethet_incr);
      dotprod_sumv = _mm256_add_epi64(dotprod_sumv, dotprod_incr);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec uhethet_acc;
    UniVec dotprod_acc;
    uhethet_acc.vw = R_CAST(VecW, uhethet_sumv);
    dotprod_acc.vw = R_CAST(VecW, dotprod_sumv);
    uhethet_dosage += UniVecHsum32(uhethet_acc);
    known_dotprod += UniVecHsum32(dotprod_acc);
  }
}

void DosagePhasedDotprodComponents(const Dosage* dosage_vec0, const Dosage* dosage_vec1, const SDosage* dphase_delta0, const SDosage* dphase_delta1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const __m256i* dosage_vvec0_iter = R_CAST(const __m256i*, dosage_vec0);
  const __m256i* dosage_vvec1_iter = R_CAST(const __m256i*, dosage_vec1);
  const __m256i* dphase_delta0_iter = R_CAST(const __m256i*, dphase_delta0);
  const __m256i* dphase_delta1_iter = R_CAST(const __m256i*, dphase_delta1);
  const __m256i all_16384 = _mm256_set1_epi16(0x4000);
  const __m256i all_32767 = _mm256_set1_epi16(0x7fff);
#    if defined(__clang__) && ((!defined(__cplusplus)) || (__cplusplus < 201103L))
  const __m256i all_n32768 = _mm256_set1_epi16(0x8000);
#    else
  const __m256i all_n32768 = _mm256_set1_epi64x(-0x7fff7fff7fff8000LL);
#    endif
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
  const __m256i all0 = _mm256_setzero_si256();
  const __m256i all1 = _mm256_cmpeq_epi16(all0, all0);
  uint64_t uhethet_dosage = 0;
  uint64_t known_dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m256i uhethet_sumv = _mm256_setzero_si256();
    __m256i dotprod_sumv = _mm256_setzero_si256();
    const __m256i* dosage_vvec0_stop;
    // dotprod_incrA and dotprod_incrB values in [0..32768]
    // 65536 * 4095 * 16 dosages per __m256i = just under 2^32
    if (vecs_left < 4095) {
      if (!vecs_left) {
        *known_dotprod_dosagep = known_dotprod;
        *uhethet_dosagep = uhethet_dosage;
        return;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[4095]);
      vecs_left -= 4095;
    }
    do {
      const __m256i dosage0 = *dosage_vvec0_iter++;
      const __m256i dosage1 = *dosage_vvec1_iter++;
      const __m256i delta0 = *dphase_delta0_iter++;
      const __m256i delta1 = *dphase_delta1_iter++;
      const __m256i invmask = _mm256_or_si256(_mm256_cmpeq_epi16(dosage0, all1), _mm256_cmpeq_epi16(dosage1, all1));
      const __m256i dosage0A = _mm256_add_epi16(dosage0, delta0);
      const __m256i dosage1A = _mm256_add_epi16(dosage1, delta1);
      const __m256i dosage0B = _mm256_sub_epi16(dosage0, delta0);
      const __m256i dosage1B = _mm256_sub_epi16(dosage1, delta1);
      const __m256i dosageA_m32769 = _mm256_add_epi16(_mm256_add_epi16(dosage0A, dosage1A), all_32767);
      const __m256i dosageB_m32769 = _mm256_add_epi16(_mm256_add_epi16(dosage0B, dosage1B), all_32767);
      const __m256i invmaskA = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosage0A, all0));
      const __m256i invmaskB = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosage0B, all0));
      const __m256i unmasked_dotprod_incrA = _mm256_sub_epi16(_mm256_max_epi16(dosageA_m32769, all1), all1);
      const __m256i unmasked_dotprod_incrB = _mm256_sub_epi16(_mm256_max_epi16(dosageB_m32769, all1), all1);
      __m256i dotprod_incrA = _mm256_andnot_si256(invmaskA, unmasked_dotprod_incrA);
      __m256i dotprod_incrB = _mm256_andnot_si256(invmaskB, unmasked_dotprod_incrB);
      dotprod_incrA = _mm256_add_epi64(_mm256_and_si256(dotprod_incrA, m16), _mm256_and_si256(_mm256_srli_epi64(dotprod_incrA, 16), m16));
      dotprod_incrB = _mm256_add_epi64(_mm256_and_si256(dotprod_incrB, m16), _mm256_and_si256(_mm256_srli_epi64(dotprod_incrB, 16), m16));
      dotprod_sumv = _mm256_add_epi64(dotprod_sumv, dotprod_incrA);
      dotprod_sumv = _mm256_add_epi64(dotprod_sumv, dotprod_incrB);

      const __m256i known0A = _mm256_abs_epi16(_mm256_sub_epi16(all_16384, dosage0A));
      const __m256i known1A = _mm256_abs_epi16(_mm256_sub_epi16(all_16384, dosage1A));
      const __m256i known0B = _mm256_abs_epi16(_mm256_sub_epi16(all_16384, dosage0B));
      const __m256i known1B = _mm256_abs_epi16(_mm256_sub_epi16(all_16384, dosage1B));
      const __m256i maxknownA = _mm256_max_epi16(known0A, known1A);
      const __m256i maxknownB = _mm256_max_epi16(known0B, known1B);
      __m256i uhethet_incr = _mm256_andnot_si256(invmask, _mm256_sub_epi16(all_n32768, _mm256_add_epi16(maxknownA, maxknownB)));
      uhethet_incr = _mm256_add_epi64(_mm256_and_si256(uhethet_incr, m16), _mm256_and_si256(_mm256_srli_epi64(uhethet_incr, 16), m16));
      uhethet_sumv = _mm256_add_epi64(uhethet_sumv, uhethet_incr);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec dotprod_acc;
    UniVec uhethet_acc;
    dotprod_acc.vw = R_CAST(VecW, dotprod_sumv);
    uhethet_acc.vw = R_CAST(VecW, uhethet_sumv);
    known_dotprod += UniVecHsum32(dotprod_acc);
    uhethet_dosage += UniVecHsum32(uhethet_acc);
  }
}

void DosagePhasedDotprodComponentsSubset(const Dosage* subset_invmask, const Dosage* dosage_vec0, const Dosage* dosage_vec1, const SDosage* dphase_delta0, const SDosage* dphase_delta1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const __m256i* invmask_iter = R_CAST(const __m256i*, subset_invmask);
  const __m256i* dosage_vvec0_iter = R_CAST(const __m256i*, dosage_vec0);
  const __m256i* dosage_vvec1_iter = R_CAST(const __m256i*, dosage_vec1);
  const __m256i* dphase_delta0_iter = R_CAST(const __m256i*, dphase_delta0);
  const __m256i* dphase_delta1_iter = R_CAST(const __m256i*, dphase_delta1);
  const __m256i all_16384 = _mm256_set1_epi16(0x4000);
  const __m256i all_32767 = _mm256_set1_epi16(0x7fff);
#    if defined(__clang__) && ((!defined(__cplusplus)) || (__cplusplus < 201103L))
  const __m256i all_n32768 = _mm256_set1_epi16(0x8000);
#    else
  const __m256i all_n32768 = _mm256_set1_epi64x(-0x7fff7fff7fff8000LL);
#    endif
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
  const __m256i all0 = _mm256_setzero_si256();
  const __m256i all1 = _mm256_cmpeq_epi16(all0, all0);
  uint64_t uhethet_dosage = 0;
  uint64_t known_dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m256i uhethet_sumv = _mm256_setzero_si256();
    __m256i dotprod_sumv = _mm256_setzero_si256();
    const __m256i* dosage_vvec0_stop;
    if (vecs_left < 4095) {
      if (!vecs_left) {
        *known_dotprod_dosagep = known_dotprod;
        *uhethet_dosagep = uhethet_dosage;
        return;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[4095]);
      vecs_left -= 4095;
    }
    do {
      __m256i invmask = *invmask_iter++;
      const __m256i dosage0 = *dosage_vvec0_iter++;
      const __m256i dosage1 = *dosage_vvec1_iter++;
      const __m256i delta0 = *dphase_delta0_iter++;
      const __m256i delta1 = *dphase_delta1_iter++;
      invmask = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosage0, all1));
      invmask = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosage1, all1));
      const __m256i dosage0A = _mm256_add_epi16(dosage0, delta0);
      const __m256i dosage1A = _mm256_add_epi16(dosage1, delta1);
      const __m256i dosage0B = _mm256_sub_epi16(dosage0, delta0);
      const __m256i dosage1B = _mm256_sub_epi16(dosage1, delta1);
      const __m256i dosageA_m32769 = _mm256_add_epi16(_mm256_add_epi16(dosage0A, dosage1A), all_32767);
      const __m256i dosageB_m32769 = _mm256_add_epi16(_mm256_add_epi16(dosage0B, dosage1B), all_32767);
      const __m256i invmaskA = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosage0A, all0));
      const __m256i invmaskB = _mm256_or_si256(invmask, _mm256_cmpeq_epi16(dosage0B, all0));
      const __m256i unmasked_dotprod_incrA = _mm256_sub_epi16(_mm256_max_epi16(dosageA_m32769, all1), all1);
      const __m256i unmasked_dotprod_incrB = _mm256_sub_epi16(_mm256_max_epi16(dosageB_m32769, all1), all1);
      __m256i dotprod_incrA = _mm256_andnot_si256(invmaskA, unmasked_dotprod_incrA);
      __m256i dotprod_incrB = _mm256_andnot_si256(invmaskB, unmasked_dotprod_incrB);
      dotprod_incrA = _mm256_add_epi64(_mm256_and_si256(dotprod_incrA, m16), _mm256_and_si256(_mm256_srli_epi64(dotprod_incrA, 16), m16));
      dotprod_incrB = _mm256_add_epi64(_mm256_and_si256(dotprod_incrB, m16), _mm256_and_si256(_mm256_srli_epi64(dotprod_incrB, 16), m16));
      dotprod_sumv = _mm256_add_epi64(dotprod_sumv, dotprod_incrA);
      dotprod_sumv = _mm256_add_epi64(dotprod_sumv, dotprod_incrB);

      const __m256i known0A = _mm256_abs_epi16(_mm256_sub_epi16(all_16384, dosage0A));
      const __m256i known1A = _mm256_abs_epi16(_mm256_sub_epi16(all_16384, dosage1A));
      const __m256i known0B = _mm256_abs_epi16(_mm256_sub_epi16(all_16384, dosage0B));
      const __m256i known1B = _mm256_abs_epi16(_mm256_sub_epi16(all_16384, dosage1B));
      const __m256i maxknownA = _mm256_max_epi16(known0A, known1A);
      const __m256i maxknownB = _mm256_max_epi16(known0B, known1B);
      __m256i uhethet_incr = _mm256_andnot_si256(invmask, _mm256_sub_epi16(all_n32768, _mm256_add_epi16(maxknownA, maxknownB)));
      uhethet_incr = _mm256_add_epi64(_mm256_and_si256(uhethet_incr, m16), _mm256_and_si256(_mm256_srli_epi64(uhethet_incr, 16), m16));
      uhethet_sumv = _mm256_add_epi64(uhethet_sumv, uhethet_incr);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec dotprod_acc;
    UniVec uhethet_acc;
    dotprod_acc.vw = R_CAST(VecW, dotprod_sumv);
    uhethet_acc.vw = R_CAST(VecW, uhethet_sumv);
    known_dotprod += UniVecHsum32(dotprod_acc);
    uhethet_dosage += UniVecHsum32(uhethet_acc);
  }
}
#  else  // !USE_AVX2
void FillDosageHet(const Dosage* dosage_vec, uint32_t dosagev_ct, Dosage* dosage_het) {
  const __m128i* dosage_vvec_iter = R_CAST(const __m128i*, dosage_vec);
#    if defined(__clang__) && ((!defined(__cplusplus)) || (__cplusplus < 201103L))
  const __m128i all_n32768 = _mm_set1_epi16(0x8000);
  const __m128i all_n16384 = _mm_set1_epi16(0xc000);
#    else
  const __m128i all_n32768 = _mm_set1_epi64x(-0x7fff7fff7fff8000LL);
  const __m128i all_n16384 = _mm_set1_epi64x(-0x3fff3fff3fff4000LL);
#    endif
  const __m128i all0 = _mm_setzero_si128();
  const __m128i all1 = _mm_cmpeq_epi16(all0, all0);
  // 0-16384: leave unchanged
  // 16385-32768: subtract from 32768
  // 65535: set to 0

  // subtract from 0, _mm_cmplt_epi16 to produce mask, add 32768
  __m128i* dosage_het_iter = R_CAST(__m128i*, dosage_het);
  for (uint32_t vec_idx = 0; vec_idx != dosagev_ct; ++vec_idx) {
    __m128i dosagev = *dosage_vvec_iter++;

    __m128i cur_mask = _mm_cmpeq_epi16(dosagev, all1);
    dosagev = _mm_andnot_si128(cur_mask, dosagev);  // 65535 -> 0

    // xor with -32768 is same as subtracting it
    __m128i dosagev_opp = _mm_xor_si128(dosagev, all_n32768);
    // anything > -16384 after this subtraction was originally >16384.
    // calling the original value x, we want to flip the sign of (x - 32768)
    cur_mask = _mm_cmpgt_epi16(dosagev_opp, all_n16384);
    dosagev_opp = _mm_and_si128(cur_mask, dosagev_opp);
    dosagev = _mm_andnot_si128(cur_mask, dosagev);  // has the <= 16384 values
    dosagev_opp = _mm_sub_epi16(all0, dosagev_opp);
    *dosage_het_iter++ = _mm_add_epi16(dosagev, dosagev_opp);
  }
}

uint64_t DenseDosageSum(const Dosage* dosage_vec, uint32_t vec_ct) {
  // end of dosage_vec assumed to be missing-padded (0-padded also ok)
  const __m128i* dosage_vvec_iter = R_CAST(const __m128i*, dosage_vec);
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
  const __m128i all1 = _mm_cmpeq_epi16(m16, m16);
  uint64_t sum = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m128i sumv = _mm_setzero_si128();
    const __m128i* dosage_vvec_stop;
    // individual values in [0..32768]
    // 32768 * 16383 * 8 dosages per __m128i = just under 2^32
    if (vecs_left < 16383) {
      if (!vecs_left) {
        return sum;
      }
      dosage_vvec_stop = &(dosage_vvec_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec_stop = &(dosage_vvec_iter[16383]);
      vecs_left -= 16383;
    }
    do {
      __m128i dosagev = *dosage_vvec_iter++;
      __m128i invmask = _mm_cmpeq_epi16(dosagev, all1);
      dosagev = _mm_andnot_si128(invmask, dosagev);

      dosagev = _mm_add_epi64(_mm_and_si128(dosagev, m16), _mm_and_si128(_mm_srli_epi64(dosagev, 16), m16));
      sumv = _mm_add_epi64(sumv, dosagev);
    } while (dosage_vvec_iter < dosage_vvec_stop);
    UniVec acc;
    acc.vw = R_CAST(VecW, sumv);
    sum += UniVecHsum32(acc);
  }
}

uint64_t DenseDosageSumSubset(const Dosage* dosage_vec, const Dosage* dosage_mask_vec, uint32_t vec_ct) {
  // end of dosage_vec assumed to be missing-padded (0-padded also ok)
  const __m128i* dosage_vvec_iter = R_CAST(const __m128i*, dosage_vec);
  const __m128i* dosage_mask_vvec_iter = R_CAST(const __m128i*, dosage_mask_vec);
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
  const __m128i all1 = _mm_cmpeq_epi16(m16, m16);
  uint64_t sum = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m128i sumv = _mm_setzero_si128();
    const __m128i* dosage_vvec_stop;
    if (vecs_left < 16383) {
      if (!vecs_left) {
        return sum;
      }
      dosage_vvec_stop = &(dosage_vvec_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec_stop = &(dosage_vvec_iter[16383]);
      vecs_left -= 16383;
    }
    do {
      __m128i invmask = *dosage_mask_vvec_iter++;
      __m128i dosagev = *dosage_vvec_iter++;
      invmask = _mm_cmpeq_epi16(invmask, all1);
      invmask = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosagev, all1));
      dosagev = _mm_andnot_si128(invmask, dosagev);

      dosagev = _mm_add_epi64(_mm_and_si128(dosagev, m16), _mm_and_si128(_mm_srli_epi64(dosagev, 16), m16));
      sumv = _mm_add_epi64(sumv, dosagev);
    } while (dosage_vvec_iter < dosage_vvec_stop);
    UniVec acc;
    acc.vw = R_CAST(VecW, sumv);
    sum += UniVecHsum32(acc);
  }
}

// 65535 treated as missing
uint64_t DosageUnsignedDotprod(const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const __m128i* dosage_vvec0_iter = R_CAST(const __m128i*, dosage_vec0);
  const __m128i* dosage_vvec1_iter = R_CAST(const __m128i*, dosage_vec1);
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
  const __m128i all1 = _mm_cmpeq_epi16(m16, m16);
  uint64_t dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m128i dotprod_lo = _mm_setzero_si128();
    __m128i dotprod_hi = _mm_setzero_si128();
    const __m128i* dosage_vvec0_stop;
    if (vecs_left < 8192) {
      if (!vecs_left) {
        return dotprod;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[8192]);
      vecs_left -= 8192;
    }
    do {
      __m128i dosage0 = *dosage_vvec0_iter++;
      __m128i dosage1 = *dosage_vvec1_iter++;
      __m128i invmask = _mm_cmpeq_epi16(dosage0, all1);
      invmask = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosage1, all1));
      dosage0 = _mm_andnot_si128(invmask, dosage0);
      dosage1 = _mm_andnot_si128(invmask, dosage1);

      __m128i lo16 = _mm_mullo_epi16(dosage0, dosage1);
      __m128i hi16 = _mm_mulhi_epu16(dosage0, dosage1);
      lo16 = _mm_add_epi64(_mm_and_si128(lo16, m16), _mm_and_si128(_mm_srli_epi64(lo16, 16), m16));
      hi16 = _mm_and_si128(_mm_add_epi64(hi16, _mm_srli_epi64(hi16, 16)), m16);
      dotprod_lo = _mm_add_epi64(dotprod_lo, lo16);
      dotprod_hi = _mm_add_epi64(dotprod_hi, hi16);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec acc_lo;
    UniVec acc_hi;
    acc_lo.vw = R_CAST(VecW, dotprod_lo);
    acc_hi.vw = R_CAST(VecW, dotprod_hi);
    dotprod += UniVecHsum32(acc_lo) + 65536 * UniVecHsum32(acc_hi);
  }
}

uint64_t DosageUnsignedDotprodSubset(const Dosage* dosage_mask_vec, const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const __m128i* dosage_mask_iter = R_CAST(const __m128i*, dosage_mask_vec);
  const __m128i* dosage_vvec0_iter = R_CAST(const __m128i*, dosage_vec0);
  const __m128i* dosage_vvec1_iter = R_CAST(const __m128i*, dosage_vec1);
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
  const __m128i all1 = _mm_cmpeq_epi16(m16, m16);
  uint64_t dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m128i dotprod_lo = _mm_setzero_si128();
    __m128i dotprod_hi = _mm_setzero_si128();
    const __m128i* dosage_vvec0_stop;
    if (vecs_left < 8192) {
      if (!vecs_left) {
        return dotprod;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[8192]);
      vecs_left -= 8192;
    }
    do {
      __m128i cur_mask = *dosage_mask_iter++;
      __m128i dosage0 = *dosage_vvec0_iter++;
      __m128i dosage1 = *dosage_vvec1_iter++;
      __m128i invmask = _mm_cmpeq_epi16(dosage0, all1);
      invmask = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosage1, all1));
      invmask = _mm_or_si128(invmask, _mm_cmpeq_epi16(cur_mask, all1));
      dosage0 = _mm_andnot_si128(invmask, dosage0);
      dosage1 = _mm_andnot_si128(invmask, dosage1);

      __m128i lo16 = _mm_mullo_epi16(dosage0, dosage1);
      __m128i hi16 = _mm_mulhi_epu16(dosage0, dosage1);
      lo16 = _mm_add_epi64(_mm_and_si128(lo16, m16), _mm_and_si128(_mm_srli_epi64(lo16, 16), m16));
      hi16 = _mm_and_si128(_mm_add_epi64(hi16, _mm_srli_epi64(hi16, 16)), m16);
      dotprod_lo = _mm_add_epi64(dotprod_lo, lo16);
      dotprod_hi = _mm_add_epi64(dotprod_hi, hi16);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec acc_lo;
    UniVec acc_hi;
    acc_lo.vw = R_CAST(VecW, dotprod_lo);
    acc_hi.vw = R_CAST(VecW, dotprod_hi);
    dotprod += UniVecHsum32(acc_lo) + 65536 * UniVecHsum32(acc_hi);
  }
}

void DosageUnphasedDotprodComponents(const Dosage* dosage_vec0, const Dosage* dosage_vec1, const Dosage* dosage_het0, const Dosage* dosage_het1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const __m128i* dosage_vvec0_iter = R_CAST(const __m128i*, dosage_vec0);
  const __m128i* dosage_vvec1_iter = R_CAST(const __m128i*, dosage_vec1);
  const __m128i* dosage_het0_iter = R_CAST(const __m128i*, dosage_het0);
  const __m128i* dosage_het1_iter = R_CAST(const __m128i*, dosage_het1);
  const __m128i all_32767 = _mm_set1_epi16(0x7fff);
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
  const __m128i all0 = _mm_setzero_si128();
  const __m128i all1 = _mm_cmpeq_epi16(all0, all0);
  uint64_t uhethet_dosage = 0;
  uint64_t known_dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m128i uhethet_sumv = _mm_setzero_si128();
    __m128i dotprod_sumv = _mm_setzero_si128();
    const __m128i* dosage_vvec0_stop;
    if (vecs_left < 16383) {
      if (!vecs_left) {
        *known_dotprod_dosagep = known_dotprod * 2;
        *uhethet_dosagep = uhethet_dosage * 2;
        return;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[16383]);
      vecs_left -= 16383;
    }
    do {
      const __m128i dosage0_plus_32767 = _mm_add_epi16(*dosage_vvec0_iter++, all_32767);
      const __m128i dosage1 = *dosage_vvec1_iter++;
      const __m128i het0 = *dosage_het0_iter++;
      const __m128i het1 = *dosage_het1_iter++;
      const __m128i invmask = _mm_or_si128(_mm_cmpgt_epi16(dosage0_plus_32767, all0), _mm_cmpeq_epi16(dosage1, all1));
      __m128i uhethet_incr = _mm_min_epi16(het0, het1);
      const __m128i dosagesum_m32769 = _mm_add_epi16(dosage0_plus_32767, dosage1);
      const __m128i unmasked_dotprod_incr = _mm_sub_epi16(_mm_max_epi16(dosagesum_m32769, all1), all1);
      __m128i dotprod_incr = _mm_andnot_si128(invmask, unmasked_dotprod_incr);

      uhethet_incr = _mm_and_si128(_mm_add_epi64(uhethet_incr, _mm_srli_epi64(uhethet_incr, 16)), m16);
      dotprod_incr = _mm_add_epi64(_mm_and_si128(dotprod_incr, m16), _mm_and_si128(_mm_srli_epi64(dotprod_incr, 16), m16));
      uhethet_sumv = _mm_add_epi64(uhethet_sumv, uhethet_incr);
      dotprod_sumv = _mm_add_epi64(dotprod_sumv, dotprod_incr);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec uhethet_acc;
    UniVec dotprod_acc;
    uhethet_acc.vw = R_CAST(VecW, uhethet_sumv);
    dotprod_acc.vw = R_CAST(VecW, dotprod_sumv);
    uhethet_dosage += UniVecHsum32(uhethet_acc);
    known_dotprod += UniVecHsum32(dotprod_acc);
  }
}

void DosageUnphasedDotprodComponentsSubset(const Dosage* subset_invmask, const Dosage* dosage_vec0, const Dosage* dosage_vec1, const Dosage* dosage_het0, const Dosage* dosage_het1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const __m128i* invmask_iter = R_CAST(const __m128i*, subset_invmask);
  const __m128i* dosage_vvec0_iter = R_CAST(const __m128i*, dosage_vec0);
  const __m128i* dosage_vvec1_iter = R_CAST(const __m128i*, dosage_vec1);
  const __m128i* dosage_het0_iter = R_CAST(const __m128i*, dosage_het0);
  const __m128i* dosage_het1_iter = R_CAST(const __m128i*, dosage_het1);
  const __m128i all_32767 = _mm_set1_epi16(0x7fff);
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
  const __m128i all0 = _mm_setzero_si128();
  const __m128i all1 = _mm_cmpeq_epi16(all0, all0);
  uint64_t uhethet_dosage = 0;
  uint64_t known_dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m128i uhethet_sumv = _mm_setzero_si128();
    __m128i dotprod_sumv = _mm_setzero_si128();
    const __m128i* dosage_vvec0_stop;
    if (vecs_left < 16383) {
      if (!vecs_left) {
        *known_dotprod_dosagep = known_dotprod * 2;
        *uhethet_dosagep = uhethet_dosage * 2;
        return;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[16383]);
      vecs_left -= 16383;
    }
    do {
      const __m128i raw_invmask = *invmask_iter++;
      const __m128i dosage0_plus_32767 = _mm_add_epi16(*dosage_vvec0_iter++, all_32767);
      const __m128i dosage1 = *dosage_vvec1_iter++;
      const __m128i het0 = *dosage_het0_iter++;
      const __m128i het1 = *dosage_het1_iter++;
      __m128i uhethet_incr = _mm_andnot_si128(raw_invmask, _mm_min_epi16(het0, het1));
      __m128i invmask = _mm_or_si128(raw_invmask, _mm_cmpgt_epi16(dosage0_plus_32767, all0));
      invmask = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosage1, all1));
      const __m128i dosagesum_m32769 = _mm_add_epi16(dosage0_plus_32767, dosage1);
      const __m128i unmasked_dotprod_incr = _mm_sub_epi16(_mm_max_epi16(dosagesum_m32769, all1), all1);
      __m128i dotprod_incr = _mm_andnot_si128(invmask, unmasked_dotprod_incr);

      uhethet_incr = _mm_and_si128(_mm_add_epi64(uhethet_incr, _mm_srli_epi64(uhethet_incr, 16)), m16);
      dotprod_incr = _mm_add_epi64(_mm_and_si128(dotprod_incr, m16), _mm_and_si128(_mm_srli_epi64(dotprod_incr, 16), m16));
      uhethet_sumv = _mm_add_epi64(uhethet_sumv, uhethet_incr);
      dotprod_sumv = _mm_add_epi64(dotprod_sumv, dotprod_incr);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec uhethet_acc;
    UniVec dotprod_acc;
    uhethet_acc.vw = R_CAST(VecW, uhethet_sumv);
    dotprod_acc.vw = R_CAST(VecW, dotprod_sumv);
    uhethet_dosage += UniVecHsum32(uhethet_acc);
    known_dotprod += UniVecHsum32(dotprod_acc);
  }
}

void DosagePhasedDotprodComponents(const Dosage* dosage_vec0, const Dosage* dosage_vec1, const SDosage* dphase_delta0, const SDosage* dphase_delta1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const __m128i* dosage_vvec0_iter = R_CAST(const __m128i*, dosage_vec0);
  const __m128i* dosage_vvec1_iter = R_CAST(const __m128i*, dosage_vec1);
  const __m128i* dphase_delta0_iter = R_CAST(const __m128i*, dphase_delta0);
  const __m128i* dphase_delta1_iter = R_CAST(const __m128i*, dphase_delta1);
  const __m128i all_16384 = _mm_set1_epi16(0x4000);
  const __m128i all_32767 = _mm_set1_epi16(0x7fff);
#      if defined(__clang__) && ((!defined(__cplusplus)) || (__cplusplus < 201103L))
  const __m128i all_n32768 = _mm_set1_epi16(0x8000);
#      else
  const __m128i all_n32768 = _mm_set1_epi64x(-0x7fff7fff7fff8000LL);
#      endif
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
  const __m128i all0 = _mm_setzero_si128();
  const __m128i all1 = _mm_cmpeq_epi16(all0, all0);
  uint64_t known_dotprod = 0;
  uint64_t uhethet_dosage = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m128i dotprod_sumv = _mm_setzero_si128();
    __m128i uhethet_sumv = _mm_setzero_si128();
    const __m128i* dosage_vvec0_stop;
    if (vecs_left < 8191) {
      if (!vecs_left) {
        *known_dotprod_dosagep = known_dotprod;
        *uhethet_dosagep = uhethet_dosage;
        return;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[8191]);
      vecs_left -= 8191;
    }
    do {
      const __m128i dosage0 = *dosage_vvec0_iter++;
      const __m128i dosage1 = *dosage_vvec1_iter++;
      const __m128i delta0 = *dphase_delta0_iter++;
      const __m128i delta1 = *dphase_delta1_iter++;
      const __m128i invmask = _mm_or_si128(_mm_cmpeq_epi16(dosage0, all1), _mm_cmpeq_epi16(dosage1, all1));
      const __m128i dosage0A = _mm_add_epi16(dosage0, delta0);
      const __m128i dosage1A = _mm_add_epi16(dosage1, delta1);
      const __m128i dosage0B = _mm_sub_epi16(dosage0, delta0);
      const __m128i dosage1B = _mm_sub_epi16(dosage1, delta1);
      const __m128i dosageA_m32769 = _mm_add_epi16(_mm_add_epi16(dosage0A, dosage1A), all_32767);
      const __m128i dosageB_m32769 = _mm_add_epi16(_mm_add_epi16(dosage0B, dosage1B), all_32767);
      const __m128i invmaskA = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosage0A, all0));
      const __m128i invmaskB = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosage0B, all0));
      const __m128i unmasked_dotprod_incrA = _mm_sub_epi16(_mm_max_epi16(dosageA_m32769, all1), all1);
      const __m128i unmasked_dotprod_incrB = _mm_sub_epi16(_mm_max_epi16(dosageB_m32769, all1), all1);
      __m128i dotprod_incrA = _mm_andnot_si128(invmaskA, unmasked_dotprod_incrA);
      __m128i dotprod_incrB = _mm_andnot_si128(invmaskB, unmasked_dotprod_incrB);
      dotprod_incrA = _mm_add_epi64(_mm_and_si128(dotprod_incrA, m16), _mm_and_si128(_mm_srli_epi64(dotprod_incrA, 16), m16));
      dotprod_incrB = _mm_add_epi64(_mm_and_si128(dotprod_incrB, m16), _mm_and_si128(_mm_srli_epi64(dotprod_incrB, 16), m16));
      dotprod_sumv = _mm_add_epi64(dotprod_sumv, dotprod_incrA);
      dotprod_sumv = _mm_add_epi64(dotprod_sumv, dotprod_incrB);

#    ifdef USE_SSE42
      const __m128i known0A = _mm_abs_epi16(_mm_sub_epi16(all_16384, dosage0A));
      const __m128i known1A = _mm_abs_epi16(_mm_sub_epi16(all_16384, dosage1A));
      const __m128i known0B = _mm_abs_epi16(_mm_sub_epi16(all_16384, dosage0B));
      const __m128i known1B = _mm_abs_epi16(_mm_sub_epi16(all_16384, dosage1B));
#    else
      const __m128i pos0A = _mm_sub_epi16(all_16384, dosage0A);
      const __m128i neg0A = _mm_sub_epi16(dosage0A, all_16384);
      const __m128i pos1A = _mm_sub_epi16(all_16384, dosage1A);
      const __m128i neg1A = _mm_sub_epi16(dosage1A, all_16384);
      const __m128i pos0B = _mm_sub_epi16(all_16384, dosage0B);
      const __m128i neg0B = _mm_sub_epi16(dosage0B, all_16384);
      const __m128i pos1B = _mm_sub_epi16(all_16384, dosage1B);
      const __m128i neg1B = _mm_sub_epi16(dosage1B, all_16384);
      const __m128i known0A = _mm_max_epi16(pos0A, neg0A);
      const __m128i known1A = _mm_max_epi16(pos1A, neg1A);
      const __m128i known0B = _mm_max_epi16(pos0B, neg0B);
      const __m128i known1B = _mm_max_epi16(pos1B, neg1B);
#    endif
      const __m128i maxknownA = _mm_max_epi16(known0A, known1A);
      const __m128i maxknownB = _mm_max_epi16(known0B, known1B);
      __m128i uhethet_incr = _mm_andnot_si128(invmask, _mm_sub_epi16(all_n32768, _mm_add_epi16(maxknownA, maxknownB)));
      uhethet_incr = _mm_add_epi64(_mm_and_si128(uhethet_incr, m16), _mm_and_si128(_mm_srli_epi64(uhethet_incr, 16), m16));
      uhethet_sumv = _mm_add_epi64(uhethet_sumv, uhethet_incr);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec dotprod_acc;
    UniVec uhethet_acc;
    dotprod_acc.vw = R_CAST(VecW, dotprod_sumv);
    uhethet_acc.vw = R_CAST(VecW, uhethet_sumv);
    known_dotprod += UniVecHsum32(dotprod_acc);
    uhethet_dosage += UniVecHsum32(uhethet_acc);
  }
}

void DosagePhasedDotprodComponentsSubset(const Dosage* subset_invmask, const Dosage* dosage_vec0, const Dosage* dosage_vec1, const SDosage* dphase_delta0, const SDosage* dphase_delta1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const __m128i* invmask_iter = R_CAST(const __m128i*, subset_invmask);
  const __m128i* dosage_vvec0_iter = R_CAST(const __m128i*, dosage_vec0);
  const __m128i* dosage_vvec1_iter = R_CAST(const __m128i*, dosage_vec1);
  const __m128i* dphase_delta0_iter = R_CAST(const __m128i*, dphase_delta0);
  const __m128i* dphase_delta1_iter = R_CAST(const __m128i*, dphase_delta1);
  const __m128i all_16384 = _mm_set1_epi16(0x4000);
  const __m128i all_32767 = _mm_set1_epi16(0x7fff);
#      if defined(__clang__) && ((!defined(__cplusplus)) || (__cplusplus < 201103L))
  const __m128i all_n32768 = _mm_set1_epi16(0x8000);
#      else
  const __m128i all_n32768 = _mm_set1_epi64x(-0x7fff7fff7fff8000LL);
#      endif
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
  const __m128i all0 = _mm_setzero_si128();
  const __m128i all1 = _mm_cmpeq_epi16(all0, all0);
  uint64_t known_dotprod = 0;
  uint64_t uhethet_dosage = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m128i dotprod_sumv = _mm_setzero_si128();
    __m128i uhethet_sumv = _mm_setzero_si128();
    const __m128i* dosage_vvec0_stop;
    if (vecs_left < 8191) {
      if (!vecs_left) {
        *known_dotprod_dosagep = known_dotprod;
        *uhethet_dosagep = uhethet_dosage;
        return;
      }
      dosage_vvec0_stop = &(dosage_vvec0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dosage_vvec0_stop = &(dosage_vvec0_iter[8191]);
      vecs_left -= 8191;
    }
    do {
      __m128i invmask = *invmask_iter++;
      const __m128i dosage0 = *dosage_vvec0_iter++;
      const __m128i dosage1 = *dosage_vvec1_iter++;
      const __m128i delta0 = *dphase_delta0_iter++;
      const __m128i delta1 = *dphase_delta1_iter++;
      invmask = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosage0, all1));
      invmask = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosage1, all1));
      const __m128i dosage0A = _mm_add_epi16(dosage0, delta0);
      const __m128i dosage1A = _mm_add_epi16(dosage1, delta1);
      const __m128i dosage0B = _mm_sub_epi16(dosage0, delta0);
      const __m128i dosage1B = _mm_sub_epi16(dosage1, delta1);
      const __m128i dosageA_m32769 = _mm_add_epi16(_mm_add_epi16(dosage0A, dosage1A), all_32767);
      const __m128i dosageB_m32769 = _mm_add_epi16(_mm_add_epi16(dosage0B, dosage1B), all_32767);
      const __m128i invmaskA = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosage0A, all0));
      const __m128i invmaskB = _mm_or_si128(invmask, _mm_cmpeq_epi16(dosage0B, all0));
      const __m128i unmasked_dotprod_incrA = _mm_sub_epi16(_mm_max_epi16(dosageA_m32769, all1), all1);
      const __m128i unmasked_dotprod_incrB = _mm_sub_epi16(_mm_max_epi16(dosageB_m32769, all1), all1);
      __m128i dotprod_incrA = _mm_andnot_si128(invmaskA, unmasked_dotprod_incrA);
      __m128i dotprod_incrB = _mm_andnot_si128(invmaskB, unmasked_dotprod_incrB);
      dotprod_incrA = _mm_add_epi64(_mm_and_si128(dotprod_incrA, m16), _mm_and_si128(_mm_srli_epi64(dotprod_incrA, 16), m16));
      dotprod_incrB = _mm_add_epi64(_mm_and_si128(dotprod_incrB, m16), _mm_and_si128(_mm_srli_epi64(dotprod_incrB, 16), m16));
      dotprod_sumv = _mm_add_epi64(dotprod_sumv, dotprod_incrA);
      dotprod_sumv = _mm_add_epi64(dotprod_sumv, dotprod_incrB);

#    ifdef USE_SSE42
      const __m128i known0A = _mm_abs_epi16(_mm_sub_epi16(all_16384, dosage0A));
      const __m128i known1A = _mm_abs_epi16(_mm_sub_epi16(all_16384, dosage1A));
      const __m128i known0B = _mm_abs_epi16(_mm_sub_epi16(all_16384, dosage0B));
      const __m128i known1B = _mm_abs_epi16(_mm_sub_epi16(all_16384, dosage1B));
#    else
      const __m128i pos0A = _mm_sub_epi16(all_16384, dosage0A);
      const __m128i neg0A = _mm_sub_epi16(dosage0A, all_16384);
      const __m128i pos1A = _mm_sub_epi16(all_16384, dosage1A);
      const __m128i neg1A = _mm_sub_epi16(dosage1A, all_16384);
      const __m128i pos0B = _mm_sub_epi16(all_16384, dosage0B);
      const __m128i neg0B = _mm_sub_epi16(dosage0B, all_16384);
      const __m128i pos1B = _mm_sub_epi16(all_16384, dosage1B);
      const __m128i neg1B = _mm_sub_epi16(dosage1B, all_16384);
      const __m128i known0A = _mm_max_epi16(pos0A, neg0A);
      const __m128i known1A = _mm_max_epi16(pos1A, neg1A);
      const __m128i known0B = _mm_max_epi16(pos0B, neg0B);
      const __m128i known1B = _mm_max_epi16(pos1B, neg1B);
#    endif
      const __m128i maxknownA = _mm_max_epi16(known0A, known1A);
      const __m128i maxknownB = _mm_max_epi16(known0B, known1B);
      __m128i uhethet_incr = _mm_andnot_si128(invmask, _mm_sub_epi16(all_n32768, _mm_add_epi16(maxknownA, maxknownB)));
      uhethet_incr = _mm_add_epi64(_mm_and_si128(uhethet_incr, m16), _mm_and_si128(_mm_srli_epi64(uhethet_incr, 16), m16));
      uhethet_sumv = _mm_add_epi64(uhethet_sumv, uhethet_incr);
    } while (dosage_vvec0_iter < dosage_vvec0_stop);
    UniVec dotprod_acc;
    UniVec uhethet_acc;
    dotprod_acc.vw = R_CAST(VecW, dotprod_sumv);
    uhethet_acc.vw = R_CAST(VecW, uhethet_sumv);
    known_dotprod += UniVecHsum32(dotprod_acc);
    uhethet_dosage += UniVecHsum32(uhethet_acc);
  }
}
#  endif  // !USE_AVX2
#else  // !__LP64__
void FillDosageHet(const Dosage* dosage_vec, uint32_t dosagev_ct, Dosage* dosage_het) {
  const uint32_t sample_ctav = dosagev_ct * kDosagePerVec;
  for (uint32_t sample_idx = 0; sample_idx != sample_ctav; ++sample_idx) {
    const uint32_t cur_dosage = dosage_vec[sample_idx];
    uint32_t cur_hetval = cur_dosage;
    if (cur_hetval > 16384) {
      if (cur_hetval == kDosageMissing) {
        cur_hetval = 0;
      } else {
        cur_hetval = 32768 - cur_hetval;
      }
    }
    dosage_het[sample_idx] = cur_hetval;
  }
}

uint64_t DenseDosageSum(const Dosage* dosage_vec, uint32_t vec_ct) {
  const uint32_t sample_ctav = vec_ct * kDosagePerVec;
  uint64_t sum = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ctav; ++sample_idx) {
    const uint32_t cur_dosage = dosage_vec[sample_idx];
    if (cur_dosage != kDosageMissing) {
      sum += cur_dosage;
    }
  }
  return sum;
}

uint64_t DenseDosageSumSubset(const Dosage* dosage_vec, const Dosage* dosage_mask_vec, uint32_t vec_ct) {
  const uint32_t sample_ctav = vec_ct * kDosagePerVec;
  uint64_t sum = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ctav; ++sample_idx) {
    const uint32_t cur_dosage = dosage_vec[sample_idx];
    const uint32_t other_dosage = dosage_mask_vec[sample_idx];
    if ((cur_dosage != kDosageMissing) && (other_dosage != kDosageMissing)) {
      sum += cur_dosage;
    }
  }
  return sum;
}

uint64_t DosageUnsignedDotprod(const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const uint32_t sample_ctav = vec_ct * kDosagePerVec;
  uint64_t dotprod = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ctav; ++sample_idx) {
    const uint32_t cur_dosage0 = dosage_vec0[sample_idx];
    const uint32_t cur_dosage1 = dosage_vec1[sample_idx];
    if ((cur_dosage0 != kDosageMissing) && (cur_dosage1 != kDosageMissing)) {
      dotprod += cur_dosage0 * cur_dosage1;
    }
  }
  return dotprod;
}

uint64_t DosageUnsignedDotprodSubset(const Dosage* dosage_mask_vec, const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const uint32_t sample_ctav = vec_ct * kDosagePerVec;
  uint64_t dotprod = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ctav; ++sample_idx) {
    const uint32_t cur_dosage_maskval = dosage_mask_vec[sample_idx];
    const uint32_t cur_dosage0 = dosage_vec0[sample_idx];
    const uint32_t cur_dosage1 = dosage_vec1[sample_idx];
    if ((cur_dosage_maskval != kDosageMissing) && (cur_dosage0 != kDosageMissing) && (cur_dosage1 != kDosageMissing)) {
      dotprod += cur_dosage0 * cur_dosage1;
    }
  }
  return dotprod;
}

void DosageUnphasedDotprodComponents(const Dosage* dosage_vec0, const Dosage* dosage_vec1, const Dosage* dosage_het0, const Dosage* dosage_het1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const uint32_t sample_ctav = vec_ct * kDosagePerVec;
  // Multiply by 2 later so unphased and phased cases are on the same scale.
  // (Specifically, known_dotprod is in 1/kDosageMax increments and the logical
  // value ranges from 0..2n, where n is founder_ct; and unknown_hethet is in
  // 1/kDosageMax increments and ranges from 0..n.)
  uint64_t half_known_dotprod_dosage = 0;
  uint64_t half_uhethet_dosage = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ctav; ++sample_idx) {
    const int32_t cur_dosage0 = dosage_vec0[sample_idx];
    const int32_t cur_dosage1 = dosage_vec1[sample_idx];
    if ((cur_dosage0 != kDosageMissing) && (cur_dosage1 != kDosageMissing)) {
      half_known_dotprod_dosage += MAXV(0, cur_dosage0 + cur_dosage1 - S_CAST(int32_t, kDosageMax));
      half_uhethet_dosage += MINV(dosage_het0[sample_idx], dosage_het1[sample_idx]);
    }
  }
  *known_dotprod_dosagep = half_known_dotprod_dosage * 2;
  *uhethet_dosagep = half_uhethet_dosage * 2;
}

// subset_invmask values required to be in {0, 65535}, for the sake of the
// vectorized implementations
void DosageUnphasedDotprodComponentsSubset(const Dosage* subset_invmask, const Dosage* dosage_vec0, const Dosage* dosage_vec1, const Dosage* dosage_het0, const Dosage* dosage_het1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const uint32_t sample_ctav = vec_ct * kDosagePerVec;
  // Multiply by 2 later so unphased and phased cases are on the same scale.
  // (Specifically, known_dotprod is in 1/kDosageMax increments and the logical
  // value ranges from 0..2n, where n is founder_ct; and unknown_hethet is in
  // 1/kDosageMax increments and ranges from 0..n.)
  uint64_t half_known_dotprod_dosage = 0;
  uint64_t half_uhethet_dosage = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ctav; ++sample_idx) {
    const int32_t cur_dosage0 = dosage_vec0[sample_idx] | subset_invmask[sample_idx];
    const int32_t cur_dosage1 = dosage_vec1[sample_idx];
    if ((cur_dosage0 != kDosageMissing) && (cur_dosage1 != kDosageMissing)) {
      half_known_dotprod_dosage += MAXV(0, cur_dosage0 + cur_dosage1 - S_CAST(int32_t, kDosageMax));
      half_uhethet_dosage += MINV(dosage_het0[sample_idx], dosage_het1[sample_idx]);
    }
  }
  *known_dotprod_dosagep = half_known_dotprod_dosage * 2;
  *uhethet_dosagep = half_uhethet_dosage * 2;
}

void DosagePhasedDotprodComponents(const Dosage* dosage_vec0, const Dosage* dosage_vec1, const SDosage* dphase_delta0, const SDosage* dphase_delta1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const uint32_t sample_ctav = vec_ct * kDosagePerVec;
  uint64_t known_dotprod_dosage = 0;
  uint64_t uhethet_dosage = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ctav; ++sample_idx) {
    const uint32_t cur_dosage0 = dosage_vec0[sample_idx];
    const uint32_t cur_dosage1 = dosage_vec1[sample_idx];
    if ((cur_dosage0 != kDosageMissing) && (cur_dosage1 != kDosageMissing)) {
      // dosage0 = a + b
      // delta0 = a - b
      // -> dosage + delta = 2a
      //    dosage - delta = 2b
      const int32_t cur_delta0 = dphase_delta0[sample_idx];
      const int32_t cur_delta1 = dphase_delta1[sample_idx];
      const int32_t dosage0A_x2 = cur_dosage0 + cur_delta0;
      const int32_t dosage1A_x2 = cur_dosage1 + cur_delta1;
      const int32_t dosage0B_x2 = cur_dosage0 - cur_delta0;
      const int32_t dosage1B_x2 = cur_dosage1 - cur_delta1;
      known_dotprod_dosage += MAXV(0, dosage0A_x2 + dosage1A_x2 - S_CAST(int32_t, kDosageMax)) + MAXV(0, dosage0B_x2 + dosage1B_x2 - S_CAST(int32_t, kDosageMax));
      const uint32_t known0A = abs_i32(kDosageMid - dosage0A_x2);
      const uint32_t known1A = abs_i32(kDosageMid - dosage1A_x2);
      const uint32_t known0B = abs_i32(kDosageMid - dosage0B_x2);
      const uint32_t known1B = abs_i32(kDosageMid - dosage1B_x2);
      uhethet_dosage += kDosageMax - MAXV(known0A, known1A) - MAXV(known0B, known1B);
    }
  }
  *known_dotprod_dosagep = known_dotprod_dosage;
  *uhethet_dosagep = uhethet_dosage;
}

void DosagePhasedDotprodComponentsSubset(const Dosage* subset_invmask, const Dosage* dosage_vec0, const Dosage* dosage_vec1, const SDosage* dphase_delta0, const SDosage* dphase_delta1, uint32_t vec_ct, uint64_t* known_dotprod_dosagep, uint64_t* uhethet_dosagep) {
  const uint32_t sample_ctav = vec_ct * kDosagePerVec;
  uint64_t known_dotprod_dosage = 0;
  uint64_t uhethet_dosage = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ctav; ++sample_idx) {
    const uint32_t cur_dosage0 = dosage_vec0[sample_idx] | subset_invmask[sample_idx];
    const uint32_t cur_dosage1 = dosage_vec1[sample_idx];
    if ((cur_dosage0 != kDosageMissing) && (cur_dosage1 != kDosageMissing)) {
      // dosage0 = a + b
      // delta0 = a - b
      // -> dosage + delta = 2a
      //    dosage - delta = 2b
      const int32_t cur_delta0 = dphase_delta0[sample_idx];
      const int32_t cur_delta1 = dphase_delta1[sample_idx];
      const int32_t dosage0A_x2 = cur_dosage0 + cur_delta0;
      const int32_t dosage1A_x2 = cur_dosage1 + cur_delta1;
      const int32_t dosage0B_x2 = cur_dosage0 - cur_delta0;
      const int32_t dosage1B_x2 = cur_dosage1 - cur_delta1;
      known_dotprod_dosage += MAXV(0, dosage0A_x2 + dosage1A_x2 - S_CAST(int32_t, kDosageMax)) + MAXV(0, dosage0B_x2 + dosage1B_x2 - S_CAST(int32_t, kDosageMax));
      const uint32_t known0A = abs_i32(kDosageMid - dosage0A_x2);
      const uint32_t known1A = abs_i32(kDosageMid - dosage1A_x2);
      const uint32_t known0B = abs_i32(kDosageMid - dosage0B_x2);
      const uint32_t known1B = abs_i32(kDosageMid - dosage1B_x2);
      uhethet_dosage += kDosageMax - MAXV(known0A, known1A) - MAXV(known0B, known1B);
    }
  }
  *known_dotprod_dosagep = known_dotprod_dosage;
  *uhethet_dosagep = uhethet_dosage;
}
#endif

// nmaj_dosages assumed to be initialized to whole-vector sums
uint32_t DosageR2Freqs(const Dosage* dosage_vec0, const uintptr_t* nm_bitvec0, const Dosage* dosage_vec1, const uintptr_t* nm_bitvec1, uint32_t sample_ct, uint32_t nm_ct0, uint32_t nm_ct1, uint64_t* __restrict nmaj_dosages) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  uint32_t nm_intersection_ct;
  if ((nm_ct0 != sample_ct) && (nm_ct1 != sample_ct)) {
    nm_intersection_ct = PopcountWordsIntersect(nm_bitvec0, nm_bitvec1, sample_ctl);
    if (!nm_intersection_ct) {
      nmaj_dosages[0] = 0;
      nmaj_dosages[1] = 0;
      return 0;
    }
  } else {
    nm_intersection_ct = MINV(nm_ct0, nm_ct1);
  }
  const uint32_t vec_ct = DivUp(sample_ct, kDosagePerVec);
  if (nm_ct0 != nm_intersection_ct) {
    nmaj_dosages[0] = DenseDosageSumSubset(dosage_vec0, dosage_vec1, vec_ct);
  }
  if (nm_ct1 != nm_intersection_ct) {
    nmaj_dosages[1] = DenseDosageSumSubset(dosage_vec1, dosage_vec0, vec_ct);
  }
  return nm_intersection_ct;
}

// nmaj_dosages assumed to be initialized to subsetted-vector sums.
uint32_t DosageR2FreqsSubset(const Dosage* dosage_vec0, const uintptr_t* nm_bitvec0, const Dosage* dosage_vec1, const uintptr_t* nm_bitvec1, const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t subsetted_nm_ct0, uint32_t subsetted_nm_ct1, uint64_t* __restrict nmaj_dosages, uintptr_t* cur_nm_buf, Dosage* invmask_buf) {
  // bugfix (29 Oct 2023): got raw_sample_ct / sample_ct mixed up
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  if (subsetted_nm_ct0 == sample_ct) {
    if (subsetted_nm_ct1 == sample_ct) {
      return sample_ct;
    }
    BitvecAndCopy(nm_bitvec1, sample_include, raw_sample_ctl, cur_nm_buf);
  } else {
    BitvecAndCopy(nm_bitvec0, sample_include, raw_sample_ctl, cur_nm_buf);
    if (subsetted_nm_ct1 != sample_ct) {
      BitvecAnd(nm_bitvec1, raw_sample_ctl, cur_nm_buf);
    }
  }
  const uint32_t nm_intersection_ct = PopcountWords(cur_nm_buf, raw_sample_ctl);
  if (!nm_intersection_ct) {
    nmaj_dosages[0] = 0;
    nmaj_dosages[1] = 0;
    return 0;
  }
  Expand1bitTo16(cur_nm_buf, RoundUpPow2(raw_sample_ct, kDosagePerVec), 0xffff, invmask_buf);
  const uint32_t vec_ct = DivUp(raw_sample_ct, kDosagePerVec);
  nmaj_dosages[0] = DenseDosageSumSubset(dosage_vec0, invmask_buf, vec_ct);
  nmaj_dosages[1] = DenseDosageSumSubset(dosage_vec1, invmask_buf, vec_ct);
  return nm_intersection_ct;
}


// "unscaled" because you need to multiply by allele count to get the proper
// log-likelihood
double EmPhaseUnscaledLnlike(double freq11, double freq12, double freq21, double freq22, double half_hethet_share, double freq11_incr) {
  // bugfix (11 Dec 2023): can't modify freq11, etc. in-place because then
  // they're no longer scaled with old calc_lnlike() known11, etc.
  const double adj_freq11 = freq11 + freq11_incr;
  const double adj_freq22 = freq22 + freq11_incr;
  const double adj_freq12 = freq12 + half_hethet_share - freq11_incr;
  const double adj_freq21 = freq21 + half_hethet_share - freq11_incr;
  const double cross_sum = adj_freq11 * adj_freq22 + adj_freq12 * adj_freq21;
  double lnlike = 0.0;
  if (cross_sum != 0.0) {
    lnlike = half_hethet_share * log(cross_sum);
  }
  if (adj_freq11 != 0.0) {
    lnlike += freq11 * log(adj_freq11);
  }
  if (adj_freq12 != 0.0) {
    lnlike += freq12 * log(adj_freq12);
  }
  if (adj_freq21 != 0.0) {
    lnlike += freq21 * log(adj_freq21);
  }
  if (adj_freq22 != 0.0) {
    lnlike += freq22 * log(adj_freq22);
  }
  return lnlike;
}

ENUM_U31_DEF_START()
  kLDErrNone,
  kLDMonomorphic0,
  kLDMonomorphic1
ENUM_U31_DEF_END(LDErr);

typedef struct PhasedLDExtraRetStruct {
  uint32_t sol_ct;
  uint32_t best_lnlike_mask;
  double freq_majmaj;
  double freq_majmin;
  double freq_minmaj;
  double freq_minmin;
  double half_unphased_hethet_share;
  double freq_majx;
  double freq_minx;
  double freq_xmaj;
  double freq_xmin;
} PhasedLDExtraRet;

// If extra_retp is non-null, results contains relevant cubic_sols.
// If extra_retp is nullptr, results[0] contains r2, and if
//   compute_d_and_dprime is true, results[1] = D and results[2] = D'.
LDErr PhasedLD(const double* nmajsums_d, double known_dotprod_d, double unknown_hethet_d, double twice_tot_recip, uint32_t compute_d_and_dprime, PhasedLDExtraRet* extra_retp, double* results, uint32_t* is_neg_ptr) {
  // known-diplotype dosages (sum is 2 * (valid_obs_d - unknown_hethet_d)):
  //   var0  var1
  //     0  -  0 : 2 * valid_obs_d - majsums[0] - majsums[1] + known_dotprod
  //     1  -  0 : majsums[0] - known_dotprod - unknown_hethet_d
  //     0  -  1 : majsums[1] - known_dotprod - unknown_hethet_d
  //     1  -  1 : known_dotprod

  // bugfix (15 Sep 2023): otherwise possible for freq_majmaj to be slightly
  // less than zero, and this breaks EmPhaseUnscaledLnlike().
  const double freq_majmaj = MAXV(1.0 - (nmajsums_d[0] + nmajsums_d[1] - known_dotprod_d) * twice_tot_recip, 0.0);
  const double freq_majmin = (nmajsums_d[1] - known_dotprod_d - unknown_hethet_d) * twice_tot_recip;
  const double freq_minmaj = (nmajsums_d[0] - known_dotprod_d - unknown_hethet_d) * twice_tot_recip;
  const double freq_minmin = known_dotprod_d * twice_tot_recip;
  const double half_unphased_hethet_share = unknown_hethet_d * twice_tot_recip;
  const double freq_majx = freq_majmaj + freq_majmin + half_unphased_hethet_share;
  const double freq_minx = 1.0 - freq_majx;
  const double freq_xmaj = freq_majmaj + freq_minmaj + half_unphased_hethet_share;
  const double freq_xmin = 1.0 - freq_xmaj;
  // frequency of ~2^{-46} is actually possible with dosages and 2 billion
  // samples, so set this threshold at 2^{-47}
  if ((freq_majx < (kSmallEpsilon * 0.125)) || (freq_minx < (kSmallEpsilon * 0.125))) {
    return kLDMonomorphic0;
  }
  if ((freq_xmaj < (kSmallEpsilon * 0.125)) || (freq_xmin < (kSmallEpsilon * 0.125))) {
    return kLDMonomorphic1;
  }
  uint32_t cubic_sol_ct = 0;
  uint32_t first_relevant_sol_idx = 0;
  uint32_t best_lnlike_mask = 0;
  STD_ARRAY_DECL(double, 3, cubic_sols);
  if (half_unphased_hethet_share != 0.0) {
    // detect degenerate cases to avoid e-17 ugliness
    if ((freq_majmaj * freq_minmin != 0.0) || (freq_majmin * freq_minmaj != 0.0)) {
      // (f11 + x)(f22 + x)(K - x) = x(f12 + K - x)(f21 + K - x)
      // (x - K)(x + f11)(x + f22) + x(x - K - f12)(x - K - f21) = 0
      //   x^3 + (f11 + f22 - K)x^2 + (f11*f22 - K*f11 - K*f22)x - K*f11*f22
      // + x^3 - (2K + f12 + f21)x^2 + (K + f12)(K + f21)x = 0
      cubic_sol_ct = CubicRealRoots(0.5 * (freq_majmaj + freq_minmin - freq_majmin - freq_minmaj - 3 * half_unphased_hethet_share), 0.5 * (freq_majmaj * freq_minmin + freq_majmin * freq_minmaj + half_unphased_hethet_share * (freq_majmin + freq_minmaj - freq_majmaj - freq_minmin + half_unphased_hethet_share)), -0.5 * half_unphased_hethet_share * freq_majmaj * freq_minmin, cubic_sols);
      if (cubic_sol_ct > 1) {
        // have encountered 7.9e-11 difference in testing.
        while (cubic_sols[cubic_sol_ct - 1] > half_unphased_hethet_share + k2m32) {
          --cubic_sol_ct;
          if (cubic_sol_ct == 1) {
            break;
          }
        }
        if (cubic_sols[cubic_sol_ct - 1] > half_unphased_hethet_share - k2m32) {
          cubic_sols[cubic_sol_ct - 1] = half_unphased_hethet_share;
        }
        while ((cubic_sols[first_relevant_sol_idx] < -k2m32) && (first_relevant_sol_idx + 1 < cubic_sol_ct)) {
          ++first_relevant_sol_idx;
        }
      }
      // bugfix (29 Jan 2025): Also need to clip up to 0 when there's only one
      // solution.
      if (cubic_sols[first_relevant_sol_idx] < k2m32) {
        cubic_sols[first_relevant_sol_idx] = 0.0;
      }
    } else {
      // At least one of {f11, f22} is zero, and one of {f12, f21} is zero.
      // Initially suppose that the zero-values are f11 and f12.  Then the
      // equality becomes
      //   x(f22 + x)(K - x) = x(K - x)(f21 + K - x)
      //   x=0 and x=K are always solutions; the rest becomes
      //     f22 + x = f21 + K - x
      //     2x = K + f21 - f22
      //     x = (K + f21 - f22)/2; in-range iff (f21 - f22) in (-K, K).
      // So far so good.  However, plink 1.9 incorrectly *always* checked
      // (f21 - f22) before 6 Oct 2017, when it needed to use all the nonzero
      // values.
      cubic_sols[0] = 0.0;
      const double nonzero_freq_xx = freq_majmaj + freq_minmin;
      const double nonzero_freq_xy = freq_majmin + freq_minmaj;
      // (current code still works if three or all four values are zero)
      if ((nonzero_freq_xx + k2m35 < half_unphased_hethet_share + nonzero_freq_xy) && (nonzero_freq_xy + k2m35 < half_unphased_hethet_share + nonzero_freq_xx)) {
        cubic_sol_ct = 3;
        cubic_sols[1] = (half_unphased_hethet_share + nonzero_freq_xy - nonzero_freq_xx) * 0.5;
        cubic_sols[2] = half_unphased_hethet_share;
      } else {
        cubic_sol_ct = 2;
        cubic_sols[1] = half_unphased_hethet_share;
      }
    }
    // cubic_sol_ct does not contain trailing too-large solutions
    if (cubic_sol_ct > first_relevant_sol_idx + 1) {
      double best_unscaled_lnlike = -DBL_MAX;
      for (uint32_t sol_idx = first_relevant_sol_idx; sol_idx < cubic_sol_ct; ++sol_idx) {
        const double cur_unscaled_lnlike = EmPhaseUnscaledLnlike(freq_majmaj, freq_majmin, freq_minmaj, freq_minmin, half_unphased_hethet_share, cubic_sols[sol_idx]);
        if (cur_unscaled_lnlike > best_unscaled_lnlike) {
          best_unscaled_lnlike = cur_unscaled_lnlike;
          best_lnlike_mask = 1 << sol_idx;
        } else if (cur_unscaled_lnlike == best_unscaled_lnlike) {
          best_lnlike_mask |= 1 << sol_idx;
        }
      }
    }
  } else {
    cubic_sol_ct = 1;
    cubic_sols[0] = 0.0;
  }
  if (!extra_retp) {
    uint32_t sol_idx = first_relevant_sol_idx;
    if (cubic_sol_ct - first_relevant_sol_idx > 1) {
      sol_idx = ctzu32(best_lnlike_mask);
    }
    const double cur_sol_xx = cubic_sols[sol_idx];
    double dd = freq_majmaj + cur_sol_xx - freq_majx * freq_xmaj;
    if (fabs(dd) < kSmallEpsilon) {
      dd = 0.0;
    }
    results[0] = dd * dd / (freq_majx * freq_xmaj * freq_minx * freq_xmin);
    *is_neg_ptr = (dd < 0.0);
    if (compute_d_and_dprime) {
      // maybe this should just always be computed since it's such a small
      // fraction of the total cost?
      results[1] = dd;
      if (dd >= 0.0) {
        results[2] = dd / MINV(freq_xmaj * freq_minx, freq_xmin * freq_majx);
      } else {
        // note that this preserves sign
        results[2] = dd / MINV(freq_xmaj * freq_majx, freq_xmin * freq_minx);
      }
    }
  } else {
    extra_retp->sol_ct = cubic_sol_ct - first_relevant_sol_idx;
    extra_retp->best_lnlike_mask = best_lnlike_mask >> first_relevant_sol_idx;
    extra_retp->freq_majmaj = freq_majmaj;
    extra_retp->freq_majmin = freq_majmin;
    extra_retp->freq_minmaj = freq_minmaj;
    extra_retp->freq_minmin = freq_minmin;
    extra_retp->half_unphased_hethet_share = half_unphased_hethet_share;
    extra_retp->freq_majx = freq_majx;
    extra_retp->freq_minx = freq_minx;
    extra_retp->freq_xmaj = freq_xmaj;
    extra_retp->freq_xmin = freq_xmin;
    for (uint32_t sol_idx = first_relevant_sol_idx; sol_idx != cubic_sol_ct; ++sol_idx) {
      results[sol_idx - first_relevant_sol_idx] = cubic_sols[sol_idx];
    }
  }
  return kLDErrNone;
}

PglErr LdConsole(const uintptr_t* variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const LdInfo* ldip, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, PgenReader* simple_pgrp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely(founder_ct < 2)) {
      logerrputs("Error: --ld requires at least two founders.  (--make-founders may come in handy\nhere.)\n");
      goto LdConsole_ret_INCONSISTENT_INPUT;
    }
    STD_ARRAY_KREF(char*, 2) ld_console_varids = ldip->ld_console_varids;
    // ok to ignore chr_mask here
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    // is_x:
    // * male het calls treated as missing hardcalls
    // * males only have half weight in all computations (or sqrt(0.5) if one
    //   variant on chrX and one variant elsewhere)
    // * SNPHWEX used for HWE stats
    //
    // is_nonx_haploid:
    // * all het calls treated as missing hardcalls
    uint32_t var_uidxs[2];
    uint32_t chr_idxs[2];
    uint32_t is_xs[2];
    uint32_t is_haploids[2];
    uint32_t y_ct = 0;
    for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
      const char* cur_varid = ld_console_varids[var_idx];
      int32_t ii = GetVariantUidxWithoutHtable(cur_varid, variant_ids, variant_include, variant_ct);
      if (unlikely(ii == -1)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --ld variant '%s' does not appear in dataset.\n", cur_varid);
        goto LdConsole_ret_INCONSISTENT_INPUT_WW;
      } else if (unlikely(ii == -2)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --ld variant '%s' appears multiple times in dataset.\n", cur_varid);
        goto LdConsole_ret_INCONSISTENT_INPUT_WW;
      }
      const uint32_t cur_var_uidx = ii;
      var_uidxs[var_idx] = cur_var_uidx;
      const uint32_t chr_idx = GetVariantChr(cip, cur_var_uidx);
      chr_idxs[var_idx] = chr_idx;
      const uint32_t is_x = (chr_idx == x_code);
      is_xs[var_idx] = is_x;
      uint32_t is_haploid = 0;
      if (IsSet(cip->haploid_mask, chr_idx)) {
        is_haploid = 1;
        y_ct += (chr_idx == y_code);
      }
      is_haploids[var_idx] = is_haploid;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    // if both unplaced, don't count as same-chromosome
    const uint32_t is_same_chr = chr_idxs[0] && (chr_idxs[0] == chr_idxs[1]);
    if (y_ct) {
      // only keep non-female founders
      // (relaxed on 12 Oct 2023 to be consistent with plink2 --set-hh-missing
      // behavior)
      uintptr_t* founder_info_tmp;
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &founder_info_tmp))) {
        goto LdConsole_ret_NOMEM;
      }
      // female = sex_nm & (~sex_male)

      // this leaves trailing bits set, but BitvecAnd with founder_info clears
      // them
      BitvecInvertCopy(sex_nm, raw_sample_ctl, founder_info_tmp);
      BitvecOr(sex_male, raw_sample_ctl, founder_info_tmp);
      BitvecAnd(founder_info, raw_sample_ctl, founder_info_tmp);
      founder_info = founder_info_tmp;
      founder_ct = PopcountWords(founder_info, raw_sample_ctl);
      if (unlikely(founder_ct < 2)) {
        logerrputs("Error: --ld requires at least two non-female founders when a chrY variant is\nspecified.  (--make-founders may come in handy here.)\n");
        goto LdConsole_ret_INCONSISTENT_INPUT;
      }
    }
    const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
    uint32_t* founder_info_cumulative_popcounts;
    PgenVariant pgvs[2];
    PreinitPgv(&(pgvs[0]));
    PreinitPgv(&(pgvs[1]));
    if (unlikely(bigstack_alloc_u32(founder_ctl, &founder_info_cumulative_popcounts) ||
                 BigstackAllocPgv(founder_ct, 0, kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent, &(pgvs[0])) ||
                 BigstackAllocPgv(founder_ct, 0, kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent, &(pgvs[1])))) {
      goto LdConsole_ret_NOMEM;
    }

    const uint32_t x_present = (is_xs[0] || is_xs[1]);
    const uint32_t founder_ctv = BitCtToVecCt(founder_ct);
    const uint32_t founder_ctaw = founder_ctv * kWordsPerVec;
    uintptr_t* sex_male_collapsed = nullptr;
    uintptr_t* sex_male_collapsed_interleaved = nullptr;
    uint32_t x_male_ct = 0;
    if (x_present) {
      if (unlikely(bigstack_alloc_w(founder_ctaw, &sex_male_collapsed) ||
                   bigstack_alloc_w(founder_ctaw, &sex_male_collapsed_interleaved))) {
        goto LdConsole_ret_NOMEM;
      }
      CopyBitarrSubset(sex_male, founder_info, founder_ct, sex_male_collapsed);
      ZeroTrailingWords(founder_ctl, sex_male_collapsed);
      FillInterleavedMaskVec(sex_male_collapsed, founder_ctv, sex_male_collapsed_interleaved);
      x_male_ct = PopcountWords(sex_male_collapsed, founder_ctaw);
    }
    FillCumulativePopcounts(founder_info, founder_ctl, founder_info_cumulative_popcounts);

    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(founder_info_cumulative_popcounts, simple_pgrp, &pssi);
    for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
      const uint32_t variant_uidx = var_uidxs[var_idx];
      // (unconditionally allocating phaseinfo/dosage_main and using the most
      // general-purpose loader makes sense when this loop only executes twice,
      // but --r2 will want to use different pgenlib loaders depending on
      // context.)

      PgenVariant* pgvp = &(pgvs[var_idx]);
      reterr = PgrGetInv1Dp(founder_info, pssi, founder_ct, variant_uidx, maj_alleles[variant_uidx], simple_pgrp, pgvp);
      if (unlikely(reterr)) {
        PgenErrPrintNV(reterr, variant_uidx);
        goto LdConsole_ret_1;
      }
      ZeroTrailingNyps(founder_ct, pgvp->genovec);
      if (is_haploids[var_idx]) {
        if (x_male_ct && is_xs[var_idx]) {
          if (pgvp->phasepresent_ct) {
            BitvecInvmask(sex_male_collapsed, founder_ctl, pgvp->phasepresent);
            pgvp->phasepresent_ct = PopcountWords(pgvp->phasepresent, founder_ctl);
          }
          if (pgvp->dphase_ct) {
            EraseMaleDphases(sex_male_collapsed, &(pgvp->dphase_ct), pgvp->dphase_present, pgvp->dphase_delta);
          }
        } else {
          pgvp->phasepresent_ct = 0;
          pgvp->dphase_ct = 0;
        }
      }
    }
    const uint32_t use_phase = is_same_chr && (pgvs[0].phasepresent_ct || pgvs[0].dphase_ct) && (pgvs[1].phasepresent_ct || pgvs[1].dphase_ct);
    // in haploid case, het -> 0.5
    const uint32_t use_dosage = pgvs[0].dosage_ct || pgvs[1].dosage_ct || is_haploids[0] || is_haploids[1];

    // values of interest:
    //   mutually-nonmissing observation count
    //   (all other values computed over mutually-nonmissing set)
    //   4 known-diplotype dosages (0..2 for each sample, in unphased het-het)
    //   (unphased het-het fractional count can be inferred)
    //   dosage sum for each variant
    uint32_t valid_x_male_ct = 0;
    double nmajsums_d[2];
    double x_male_nmajsums_d[2];
    double known_dotprod_d;
    double unknown_hethet_d;
    uint32_t valid_obs_ct;
    uint32_t hethet_hc_found = 0;
    if (!use_dosage) {
      // While we could theoretically optimize around the fact that we only
      // need to make a single phased-r^2 computation, that's silly; it makes a
      // lot more sense to use this as a testing ground for algorithms and data
      // representations suitable for --r/--r2, etc.
      uintptr_t* one_bitvecs[2];
      uintptr_t* two_bitvecs[2];
      uintptr_t* nm_bitvecs[2];
      if (unlikely(bigstack_alloc_w(founder_ctaw, &one_bitvecs[0]) ||
                   bigstack_alloc_w(founder_ctaw, &two_bitvecs[0]) ||
                   bigstack_alloc_w(founder_ctaw, &nm_bitvecs[0]) ||
                   bigstack_alloc_w(founder_ctaw, &one_bitvecs[1]) ||
                   bigstack_alloc_w(founder_ctaw, &two_bitvecs[1]) ||
                   bigstack_alloc_w(founder_ctaw, &nm_bitvecs[1]))) {
        goto LdConsole_ret_NOMEM;
      }
      uint32_t nmaj_cts[2];
      uint32_t nm_cts[2];  // ugh.
      for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
        GenoarrSplit12Nm(pgvs[var_idx].genovec, founder_ct, one_bitvecs[var_idx], two_bitvecs[var_idx], nm_bitvecs[var_idx]);
        nmaj_cts[var_idx] = GenoBitvecSum(one_bitvecs[var_idx], two_bitvecs[var_idx], founder_ctl);
        nm_cts[var_idx] = PopcountWords(nm_bitvecs[var_idx], founder_ctl);
      }
      uint32_t known_dotprod;
      uint32_t unknown_hethet_ct;
      valid_obs_ct = HardcallPhasedR2Stats(one_bitvecs[0], two_bitvecs[0], nm_bitvecs[0], one_bitvecs[1], two_bitvecs[1], nm_bitvecs[1], founder_ct, nm_cts[0], nm_cts[1], nmaj_cts, &known_dotprod, &unknown_hethet_ct);
      if (unlikely(!valid_obs_ct)) {
        goto LdConsole_ret_NO_VALID_OBSERVATIONS;
      }
      hethet_hc_found = (unknown_hethet_ct != 0);
      if (use_phase && hethet_hc_found) {
        // all that's needed for the hardcall-phase correction is:
        //   popcount(phasepresent0 & phasepresent1)
        //   popcount(phasepresent0 & phasepresent1 &
        //            (phaseinfo0 ^ phaseinfo1))
        HardcallPhasedR2Refine(pgvs[0].phasepresent, pgvs[0].phaseinfo, pgvs[1].phasepresent, pgvs[1].phaseinfo, founder_ctl, &known_dotprod, &unknown_hethet_ct);
      }
      nmajsums_d[0] = u31tod(nmaj_cts[0]);
      nmajsums_d[1] = u31tod(nmaj_cts[1]);
      known_dotprod_d = S_CAST(double, known_dotprod);
      unknown_hethet_d = u31tod(unknown_hethet_ct);
    } else {
      const uint32_t founder_dosagev_ct = DivUp(founder_ct, kDosagePerVec);
      Dosage* dosage_vecs[2];
      Dosage* dosage_hets[2];
      uintptr_t* nm_bitvecs[2];
      // founder_ct automatically rounded up as necessary
      if (unlikely(bigstack_alloc_dosage(founder_ct, &dosage_vecs[0]) ||
                   bigstack_alloc_dosage(founder_ct, &dosage_vecs[1]) ||
                   bigstack_alloc_dosage(founder_ct, &dosage_hets[0]) ||
                   bigstack_alloc_dosage(founder_ct, &dosage_hets[1]) ||
                   bigstack_alloc_w(founder_ctl, &nm_bitvecs[0]) ||
                   bigstack_alloc_w(founder_ctl, &nm_bitvecs[1]))) {
        goto LdConsole_ret_NOMEM;
      }
      uint64_t nmaj_dosages[2];
      uint32_t nm_cts[2];
      for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
        const uint32_t dosage_ct = pgvs[var_idx].dosage_ct;
        PopulateDenseDosage(pgvs[var_idx].genovec, pgvs[var_idx].dosage_present, pgvs[var_idx].dosage_main, founder_ct, dosage_ct, dosage_vecs[var_idx]);
        nmaj_dosages[var_idx] = DenseDosageSum(dosage_vecs[var_idx], founder_dosagev_ct);
        FillDosageHet(dosage_vecs[var_idx], founder_dosagev_ct, dosage_hets[var_idx]);
        GenoarrToNonmissing(pgvs[var_idx].genovec, founder_ct, nm_bitvecs[var_idx]);
        // bugfix (10 Dec 2023)
        if (dosage_ct) {
          BitvecOr(pgvs[var_idx].dosage_present, founder_ctl, nm_bitvecs[var_idx]);
        }
        nm_cts[var_idx] = PopcountWords(nm_bitvecs[var_idx], founder_ctl);
      }
      SDosage* main_dphase_deltas[2];
      main_dphase_deltas[0] = nullptr;
      main_dphase_deltas[1] = nullptr;
      if (use_phase) {
        if (unlikely(bigstack_alloc_dphase(founder_ct, &main_dphase_deltas[0]) ||
                     bigstack_alloc_dphase(founder_ct, &main_dphase_deltas[1]))) {
          goto LdConsole_ret_NOMEM;
        }
        for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
          PopulateDenseDphase(pgvs[var_idx].phasepresent, pgvs[var_idx].phaseinfo, pgvs[var_idx].dosage_present, dosage_vecs[var_idx], pgvs[var_idx].dphase_present, pgvs[var_idx].dphase_delta, founder_ct, pgvs[var_idx].phasepresent_ct, pgvs[var_idx].dosage_ct, pgvs[var_idx].dphase_ct, main_dphase_deltas[var_idx]);
        }
      }
      valid_obs_ct = DosageR2Freqs(dosage_vecs[0], nm_bitvecs[0], dosage_vecs[1], nm_bitvecs[1], founder_ct, nm_cts[0], nm_cts[1], nmaj_dosages);
      if (unlikely(!valid_obs_ct)) {
        goto LdConsole_ret_NO_VALID_OBSERVATIONS;
      }
      nmajsums_d[0] = u63tod(nmaj_dosages[0]) * kRecipDosageMid;
      nmajsums_d[1] = u63tod(nmaj_dosages[1]) * kRecipDosageMid;
      uint64_t known_dotprod_dosage;
      uint64_t uhethet_dosage;
      if (!x_male_ct) {
        if (!use_phase) {
          DosageUnphasedDotprodComponents(dosage_vecs[0], dosage_vecs[1], dosage_hets[0], dosage_hets[1], founder_dosagev_ct, &known_dotprod_dosage, &uhethet_dosage);
        } else {
          DosagePhasedDotprodComponents(dosage_vecs[0], dosage_vecs[1], main_dphase_deltas[0], main_dphase_deltas[1], founder_dosagev_ct, &known_dotprod_dosage, &uhethet_dosage);
        }
        known_dotprod_d = u63tod(known_dotprod_dosage) * kRecipDosageMax;
        unknown_hethet_d = u63tod(uhethet_dosage) * kRecipDosageMax;
      } else {
        Dosage* x_nonmale_dosage_invmask;
        Dosage* x_male_dosage_invmask;
        uintptr_t* sex_nonmale_collapsed;
        uintptr_t* nm_buf;
        Dosage* invmask_buf;
        if (unlikely(bigstack_alloc_dosage(founder_ct, &x_nonmale_dosage_invmask) ||
                     bigstack_alloc_dosage(founder_ct, &x_male_dosage_invmask) ||
                     bigstack_alloc_w(founder_ctl, &sex_nonmale_collapsed) ||
                     bigstack_alloc_w(founder_ctl, &nm_buf) ||
                     bigstack_alloc_dosage(founder_ct, &invmask_buf))) {
          goto LdConsole_ret_NOMEM;
        }
        BitvecInvertCopy(sex_male_collapsed, founder_ctl, sex_nonmale_collapsed);
        ZeroTrailingBits(founder_ct, sex_nonmale_collapsed);
        Expand1bitTo16(sex_nonmale_collapsed, RoundUpPow2(founder_ct, kDosagePerVec), 0xffff, x_nonmale_dosage_invmask);
        Expand1bitTo16(sex_male_collapsed, RoundUpPow2(founder_ct, kDosagePerVec), 0xffff, x_male_dosage_invmask);
        uint32_t x_male_nm_cts[2];
        x_male_nm_cts[0] = PopcountWordsIntersect(sex_male_collapsed, nm_bitvecs[0], founder_ctl);
        x_male_nm_cts[1] = PopcountWordsIntersect(sex_male_collapsed, nm_bitvecs[1], founder_ctl);
        uint64_t x_male_nmaj_dosages[2];
        x_male_nmaj_dosages[0] = DenseDosageSumSubset(dosage_vecs[0], x_male_dosage_invmask, founder_dosagev_ct);
        x_male_nmaj_dosages[1] = DenseDosageSumSubset(dosage_vecs[1], x_male_dosage_invmask, founder_dosagev_ct);
        valid_x_male_ct = DosageR2FreqsSubset(dosage_vecs[0], nm_bitvecs[0], dosage_vecs[1], nm_bitvecs[1], sex_male_collapsed, founder_ct, x_male_ct, x_male_nm_cts[0], x_male_nm_cts[1], x_male_nmaj_dosages, nm_buf, invmask_buf);
        if (!use_phase) {
          DosageUnphasedDotprodComponentsSubset(x_nonmale_dosage_invmask, dosage_vecs[0], dosage_vecs[1], dosage_hets[0], dosage_hets[1], founder_dosagev_ct, &known_dotprod_dosage, &uhethet_dosage);
        } else {
          DosagePhasedDotprodComponentsSubset(x_nonmale_dosage_invmask, dosage_vecs[0], dosage_vecs[1], main_dphase_deltas[0], main_dphase_deltas[1], founder_dosagev_ct, &known_dotprod_dosage, &uhethet_dosage);
        }
        uint64_t x_male_known_dotprod_dosage;
        uint64_t x_male_uhethet_dosage;
        DosageUnphasedDotprodComponentsSubset(x_male_dosage_invmask, dosage_vecs[0], dosage_vecs[1], dosage_hets[0], dosage_hets[1], founder_dosagev_ct, &x_male_known_dotprod_dosage, &x_male_uhethet_dosage);

        // males have sqrt(0.5) weight if one variant is chrX, half-weight if
        // both are chrX
        const double male_mult = (is_xs[0] && is_xs[1])? 0.5 : (0.5 * kSqrt2);
        // unknown_hethet_d = u63tod(uhethet_dosage) * kRecipDosageMax;
        known_dotprod_d = (u63tod(known_dotprod_dosage) + male_mult * u63tod(x_male_known_dotprod_dosage)) * kRecipDosageMax;
        unknown_hethet_d = (u63tod(uhethet_dosage) + male_mult * u63tod(x_male_uhethet_dosage)) * kRecipDosageMax;
        x_male_nmajsums_d[0] = u63tod(x_male_nmaj_dosages[0]) * kRecipDosageMid;
        x_male_nmajsums_d[1] = u63tod(x_male_nmaj_dosages[1]) * kRecipDosageMid;
      }
    }
    double valid_obs_d = u31tod(valid_obs_ct);
    if (valid_x_male_ct) {
      const double male_decr = (is_xs[0] && is_xs[1])? 0.5 : (1.0 - 0.5 * kSqrt2);
      nmajsums_d[0] -= male_decr * x_male_nmajsums_d[0];
      nmajsums_d[1] -= male_decr * x_male_nmajsums_d[1];
      valid_obs_d -= male_decr * u31tod(valid_x_male_ct);
    }

    const double twice_tot_recip = 0.5 / valid_obs_d;
    // in plink 1.9, "freq12" refers to first variant=1, second variant=2
    // this most closely corresponds to freq_majmin here
    PhasedLDExtraRet extra_ret;
    double cubic_sols[3];
    const LDErr ld_err = PhasedLD(nmajsums_d, known_dotprod_d, unknown_hethet_d, twice_tot_recip, 1, &extra_ret, cubic_sols, nullptr);
    if (ld_err) {
      logerrprintfww("Warning: Skipping --ld since %s is monomorphic across all valid observations.\n", ld_console_varids[ld_err - 1]);
      goto LdConsole_ret_1;
    }
    logputs("\n");
    logprintfww("--ld %s %s:\n", ld_console_varids[0], ld_console_varids[1]);
    logputs("\n");

    char* write_poststop = &(g_logbuf[80]);
    uint32_t varid_slens[2];
    uint32_t cur_allele_ct = 2;
    for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
      const uint32_t cur_variant_uidx = var_uidxs[var_idx];
      uintptr_t allele_idx_offset_base = cur_variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[cur_variant_uidx];
        cur_allele_ct = allele_idx_offsets[cur_variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);

      const char* cur_varid = ld_console_varids[var_idx];
      const uint32_t cur_varid_slen = strlen(ld_console_varids[var_idx]);
      varid_slens[var_idx] = cur_varid_slen;
      char* write_iter = memcpya(g_logbuf, cur_varid, cur_varid_slen);
      write_iter = strcpya_k(write_iter, " alleles:\n");
      *write_iter = '\0';
      WordWrapB(0);
      logputsb();
      write_iter = strcpya_k(g_logbuf, "  MAJOR = ");
      const uint32_t maj_idx = maj_alleles[cur_variant_uidx];
      const char* maj_allele = cur_alleles[maj_idx];
      const uint32_t maj_slen = strlen(maj_allele);
      uint32_t slen_limit = 70;
      if (!maj_idx) {
        write_iter = strcpya_k(write_iter, "REF = ");
        slen_limit -= 6;
      }
      if (maj_slen < slen_limit) {
        write_iter = memcpyax(write_iter, maj_allele, maj_slen, '\n');
      } else {
        write_iter = memcpya(write_iter, maj_allele, slen_limit - 3);
        write_iter = strcpya_k(write_iter, "...\n");
      }
      *write_iter = '\0';
      logputsb();
      write_iter = strcpya_k(g_logbuf, "  MINOR = ");
      uint32_t allele_idx = (maj_idx == 0)? 1 : 0;
      while (1) {
        const char* cur_allele = cur_alleles[allele_idx];
        const uint32_t cur_slen = strlen(cur_allele);
        if (S_CAST(uintptr_t, write_poststop - write_iter) <= cur_slen) {
          char* write_ellipsis_start = &(g_logbuf[76]);
          if (write_ellipsis_start > write_iter) {
            const uint32_t final_char_ct = S_CAST(uintptr_t, write_ellipsis_start - write_iter);
            memcpy(write_iter, cur_allele, final_char_ct);
          }
          write_iter = strcpya_k(write_ellipsis_start, "...");
          break;
        }
        write_iter = memcpya(write_iter, cur_allele, cur_slen);
        ++allele_idx;
        if (allele_idx == maj_idx) {
          ++allele_idx;
        }
        if (allele_idx == cur_allele_ct) {
          break;
        }
        *write_iter++ = ',';
      }
      *write_iter++ = '\n';
      *write_iter = '\0';
      logputsb();
      if (maj_idx) {
        write_iter = strcpya_k(g_logbuf, "  (REF = ");
        const char* ref_allele = cur_alleles[0];
        const uint32_t ref_slen = strlen(ref_allele);
        if (ref_slen < 70) {
          write_iter = memcpya(write_iter, ref_allele, ref_slen);
        } else {
          write_iter = memcpya(write_iter, ref_allele, 67);
          write_iter = strcpya_k(write_iter, "...");
        }
        memcpy_k(write_iter, ")\n\0", 4);
        logputsb();
      }
    }
    logputs("\n");
    char* write_iter = u32toa(valid_obs_ct, g_logbuf);
    write_iter = strcpya_k(write_iter, " valid");
    if (y_ct) {
      write_iter = strcpya_k(write_iter, " non-female");
    }
    write_iter = strcpya_k(write_iter, " sample");
    if (valid_obs_ct != 1) {
      *write_iter++ = 's';
    }
    if (valid_x_male_ct && ((!y_ct) || (valid_x_male_ct < valid_obs_ct))) {
      write_iter = strcpya_k(write_iter, " (");
      write_iter = u32toa(valid_x_male_ct, write_iter);
      write_iter = strcpya_k(write_iter, " male)");
    }
    if (unknown_hethet_d == 0.0) {
      // not currently checking for (fractional) het pairs in dosage case.
      if (!use_dosage) {
        write_iter = strcpya_k(write_iter, "; ");
        if (hethet_hc_found) {
          write_iter = strcpya_k(write_iter, "all phased");
        } else {
          // not a literal "het pair" in the haploid case, but same concept, so
          // I'll just leave the language unchanged for now
          write_iter = strcpya_k(write_iter, "no het pairs present");
        }
      }
    } else {
      write_iter = strcpya_k(write_iter, "; ");
      // ddosagetoa() assumes kDosageMax rather than kDosageMid multiplier
      const uint64_t unknown_hethet_int_ddosage = S_CAST(int64_t, unknown_hethet_d * kDosageMax);
      char* write_iter2 = ddosagetoa(unknown_hethet_int_ddosage, write_iter);
      const uint32_t wrote_exactly_1 = (write_iter2 == &(write_iter[1])) && (write_iter[0] == '1');
      write_iter = strcpya_k(write_iter2, " het pair");
      if (!wrote_exactly_1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " statistically phased");
    }
    assert(write_iter - g_logbuf < 78);
    memcpy_k(write_iter, ".\n\0", 4);
    logputsb();

    // sol_ct does not contain trailing too-large solutions
    const uint32_t sol_ct = extra_ret.sol_ct;
    if (sol_ct > 1) {
      logputs("Multiple phasing solutions; sample size, HWE, or random mating assumption may\nbe violated.\n\nHWE exact test p-values\n-----------------------\n");
      // (can't actually get here in nonx_haploid_or_mt case, impossible to
      // have a hethet)

      const uint32_t hwe_midp = (ldip->ld_console_flags / kfLdConsoleHweMidp) & 1;
      uint32_t x_nosex_ct = 0;  // usually shouldn't exist, but...
      uintptr_t* nosex_collapsed = nullptr;
      if (x_present) {
        x_nosex_ct = founder_ct - PopcountWordsIntersect(founder_info, sex_nm, raw_sample_ctl);
        if (x_nosex_ct) {
          if (unlikely(bigstack_alloc_w(founder_ctl, &nosex_collapsed))) {
            goto LdConsole_ret_NOMEM;
          }
          CopyBitarrSubset(sex_nm, founder_info, founder_ct, nosex_collapsed);
          // bugfix (17 Oct 2023): need to pass bit count, not word count
          AlignedBitarrInvert(founder_ct, nosex_collapsed);
        }
      }
      // Unlike plink 1.9, we don't restrict these HWE computations to the
      // nonmissing intersection.
      for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
        const uintptr_t* cur_genovec = pgvs[var_idx].genovec;
        STD_ARRAY_DECL(uint32_t, 4, genocounts);
        GenoarrCountFreqsUnsafe(cur_genovec, founder_ct, genocounts);
        write_iter = strcpya_k(g_logbuf, "  ");
        write_iter = strcpya(write_iter, ld_console_varids[var_idx]);
        write_iter = strcpya_k(write_iter, ": ");
        double hwe_ln_pval;
        if (!is_xs[var_idx]) {
          hwe_ln_pval = HweLnP(genocounts[1], genocounts[0], genocounts[2], hwe_midp);
        } else {
          STD_ARRAY_DECL(uint32_t, 4, male_genocounts);
          GenoarrCountSubsetFreqs(cur_genovec, sex_male_collapsed_interleaved, founder_ct, x_male_ct, male_genocounts);
          assert(!male_genocounts[1]);
          if (x_nosex_ct) {
            STD_ARRAY_DECL(uint32_t, 4, nosex_genocounts);
            GenoarrCountSubsetFreqs2(cur_genovec, nosex_collapsed, founder_ct, x_nosex_ct, nosex_genocounts);
            genocounts[0] -= nosex_genocounts[0];
            genocounts[1] -= nosex_genocounts[1];
            genocounts[2] -= nosex_genocounts[2];
          }
          hwe_ln_pval = HweXchrLnP(genocounts[1], genocounts[0] - male_genocounts[0], genocounts[2] - male_genocounts[2], male_genocounts[0], male_genocounts[2], hwe_midp);
        }
        write_iter = lntoa_g(hwe_ln_pval, write_iter);
        memcpy_k(write_iter, "\n", 2);
        logputsb();
      }
    }
    logputs("\n");

    const uint32_t best_lnlike_mask = extra_ret.best_lnlike_mask;
    const double freq_majmaj = extra_ret.freq_majmaj;
    const double freq_majmin = extra_ret.freq_majmin;
    const double freq_minmaj = extra_ret.freq_minmaj;
    const double freq_minmin = extra_ret.freq_minmin;
    const double half_unphased_hethet_share = extra_ret.half_unphased_hethet_share;
    const double freq_majx = extra_ret.freq_majx;
    const double freq_minx = extra_ret.freq_minx;
    const double freq_xmaj = extra_ret.freq_xmaj;
    const double freq_xmin = extra_ret.freq_xmin;
    for (uint32_t sol_idx = 0; sol_idx < sol_ct; ++sol_idx) {
      if (sol_ct > 1) {
        write_iter = strcpya_k(g_logbuf, "Solution #");
        write_iter = u32toa(sol_idx + 1, write_iter);
        if ((best_lnlike_mask >> sol_idx) & 1) {
          write_iter = strcpya_k(write_iter, " (");
          if (best_lnlike_mask & ((1 << sol_idx) - 1)) {
            write_iter = strcpya_k(write_iter, "tied for ");
          }
          write_iter = strcpya_k(write_iter, "best likelihood)");
        }
        assert(write_iter - g_logbuf < 78);
        memcpy_k(write_iter, ":\n\0", 4);
        logputsb();
      }
      const double cur_sol_xx = cubic_sols[sol_idx];
      double dd = freq_majmaj + cur_sol_xx - freq_majx * freq_xmaj;
      if (fabs(dd) < kSmallEpsilon) {
        dd = 0.0;
      }
      write_iter = strcpya_k(g_logbuf, "  r^2 = ");
      write_iter = dtoa_g(dd * dd / (freq_majx * freq_xmaj * freq_minx * freq_xmin), write_iter);
      write_iter = strcpya_k(write_iter, "    |D'| = ");
      double d_prime;
      if (dd >= 0.0) {
        d_prime = dd / MINV(freq_xmaj * freq_minx, freq_xmin * freq_majx);
      } else {
        d_prime = -dd / MINV(freq_xmaj * freq_majx, freq_xmin * freq_minx);
      }
      write_iter = dtoa_g(d_prime, write_iter);
      assert(write_iter - g_logbuf < 79);
      strcpy_k(write_iter, "\n");
      logputsb();

      logputs("\n");

      // Default layout:
      // [8 spaces]Frequencies      :        [centered varID[1]]
      //     (expectations under LE)          MAJOR       MINOR
      //                                    ----------  ----------
      //                              MAJOR  a.bcdefg    a.bcdefg
      //                                    (a.bcdefg)  (a.bcdefg)
      //       [r-justified varID[0]]
      //                              MINOR  a.bcdefg    a.bcdefg
      //                                    (a.bcdefg)  (a.bcdefg)
      //
      // (decimals are fixed-point, and trailing zeroes are erased iff there is
      // an exact match to ~13-digit precision; this is slightly more stringent
      // than plink 1.9's dtoa_f_w9p6_spaced() since there isn't much room here
      // for floating-point error to accumulate)
      // As for long variant IDs:
      // The default layout uses 55 columns, and stops working when
      // strlen(varID[0]) > 26.  So the right half can be shifted up to 24
      // characters before things get ugly in terminal windows.  Thus, once
      // string length > 50, we print only the first 47 characters of varID and
      // follow it with "...".
      // Similarly, when strlen(varID[1]) <= 51, centering is pretty
      // straightforward; beyond that, we print only the first 48 chars.
      uint32_t extra_initial_spaces = 0;
      const uint32_t varid_slen0 = varid_slens[0];
      if (varid_slen0 > 26) {
        // formatting fix (1 Feb 2018): cap this at 24, not 50
        extra_initial_spaces = MINV(varid_slen0 - 26, 24);
      }
      write_iter = strcpya_k(g_logbuf, "        Frequencies      :  ");
      // default center column index is 43 + extra_initial_spaces; we're
      //   currently at column 28
      // for length-1, we want to occupy just the center column index; for
      //   length-2, both center and (center + 1), etc.
      // ((16 + extra_initial_spaces) * 2 - strlen(varID[1])) / 2
      const uint32_t varid_slen1 = varid_slens[1];
      if (varid_slen1 > 51) {
        write_iter = memcpya(write_iter, ld_console_varids[1], 48);
        write_iter = strcpya_k(write_iter, "...");
      } else {
        uint32_t offset_x2 = (16 + extra_initial_spaces) * 2;
        if (offset_x2 > varid_slen1) {
          uint32_t varid1_padding = (offset_x2 - varid_slen1) / 2;
          if (varid1_padding + varid_slen1 > 51) {
            varid1_padding = 51 - varid_slen1;
          }
          write_iter = memseta(write_iter, 32, varid1_padding);
        }
        write_iter = memcpya(write_iter, ld_console_varids[1], varid_slen1);
      }
      strcpy_k(write_iter, "\n");
      logputsb();

      write_iter = strcpya_k(g_logbuf, "  (expectations under LE)");
      write_iter = memseta(write_iter, 32, extra_initial_spaces + 10);
      snprintf(write_iter, 81 - 25 - 24 - 10, "MAJOR       MINOR\n");
      logputsb();

      write_iter = memseta(g_logbuf, 32, extra_initial_spaces + 33);
      snprintf(write_iter, 81 - 24 - 33, "----------  ----------\n");
      logputsb();

      write_iter = strcpya_k(&(g_logbuf[27 + extra_initial_spaces]), "MAJOR  ");
      write_iter = dtoa_f_probp6_spaced(freq_majmaj + cur_sol_xx, write_iter);
      write_iter = strcpya_k(write_iter, "    ");
      const double cur_sol_xy = half_unphased_hethet_share - cur_sol_xx;
      write_iter = dtoa_f_probp6_clipped(freq_majmin + cur_sol_xy, write_iter);
      strcpy_k(write_iter, "\n");
      logputsb();

      write_iter = strcpya_k(&(g_logbuf[27 + extra_initial_spaces]), "      (");
      char* next_paren_open = &(write_iter[11]);
      write_iter = dtoa_f_probp6_clipped(freq_xmaj * freq_majx, write_iter);
      *write_iter++ = ')';
      memset(write_iter, ' ', next_paren_open - write_iter);
      *next_paren_open++ = '(';
      write_iter = dtoa_f_probp6_clipped(freq_xmin * freq_majx, next_paren_open);
      memcpy_k(write_iter, ")\n\0", 4);
      logputsb();

      write_iter = g_logbuf;
      if (varid_slen0 < 26) {
        write_iter = &(write_iter[26 - varid_slen0]);
      }
      write_iter = memcpya(write_iter, ld_console_varids[0], varid_slen0);
      strcpy_k(write_iter, "\n");
      logputsb();

      write_iter = memseta(g_logbuf, 32, 27 + extra_initial_spaces);
      write_iter = strcpya_k(write_iter, "MINOR  ");
      write_iter = dtoa_f_probp6_spaced(freq_minmaj + cur_sol_xy, write_iter);
      write_iter = strcpya_k(write_iter, "    ");
      write_iter = dtoa_f_probp6_clipped(freq_minmin + cur_sol_xx, write_iter);
      strcpy_k(write_iter, "\n");
      logputsb();

      write_iter = strcpya_k(&(g_logbuf[27 + extra_initial_spaces]), "      (");
      next_paren_open = &(write_iter[11]);
      write_iter = dtoa_f_probp6_clipped(freq_xmaj * freq_minx, write_iter);
      *write_iter++ = ')';
      memset(write_iter, ' ', next_paren_open - write_iter);
      *next_paren_open++ = '(';
      write_iter = dtoa_f_probp6_clipped(freq_xmin * freq_minx, next_paren_open);
      memcpy_k(write_iter, ")\n\0", 4);
      logputsb();

      logputs("\n");
      if (dd > 0.0) {
        logputs("  Major alleles are in phase with each other.\n\n");
      } else if (dd < 0.0) {
        logputs("  Major alleles are out of phase with each other.\n\n");
      }
    }
  }
  while (0) {
  LdConsole_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LdConsole_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  LdConsole_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  LdConsole_ret_NO_VALID_OBSERVATIONS:
    logerrputs("Error: No valid observations for --ld.\n");
    reterr = kPglRetDegenerateData;
    break;
  }
 LdConsole_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

// distinguish fully-unphased r^2 computation from no-phased-calls subcase in
// phased-r^2 computation
ENUM_U31_DEF_START()
  kR2PhaseTypeUnphased,
  kR2PhaseTypeOmit,
  kR2PhaseTypePresent
ENUM_U31_DEF_END(R2PhaseType);

static inline R2PhaseType R2PhaseOmit(R2PhaseType phase_type) {
  // convert kR2PhaseTypePresent to kR2PhaseTypeOmit, leave other values
  // unchanged
  return S_CAST(R2PhaseType, phase_type != kR2PhaseTypeUnphased);
}

static inline R2PhaseType GetR2PhaseType(uint32_t phased_r2, uint32_t phase_present) {
  return S_CAST(R2PhaseType, phased_r2 * (1 + phase_present));
}

// ***** --clump implementation starts here *****
// We want to support --glm's handling of multiallelic variants, where we're
// working in terms of (variant, allele) pairs rather than just variants.
// So we maintain variant-based *and* allele-based bitvectors tracking the
// pairs under consideration.
//
// Current multithreading strategies are restricted to a single index variant
// at a time across all threads.  Possible todo: implement a strategy similar
// to --indep-pair{phase,wise}, where each thread works on its own island, and
// benchmark on typical datasets against the existing implementation.

uint32_t GetNextIsland(const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* icandidate_vbitvec, const uintptr_t* observed_variants, uint32_t raw_variant_ct, uint32_t bp_radius, uint32_t* vidx_endp, uint32_t* chr_fo_idxp) {
  uint32_t vidx_start = *vidx_endp;
  uint32_t index_vidx = AdvBoundedTo1Bit(icandidate_vbitvec, vidx_start, raw_variant_ct);
  if (index_vidx == raw_variant_ct) {
    // ensure chr_fo_idx compares unequal
    *chr_fo_idxp = UINT32_MAX;
    return raw_variant_ct;
  }
  uint32_t chr_fo_idx = *chr_fo_idxp;
  uint32_t chr_end_vidx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  if (index_vidx >= chr_end_vidx) {
    chr_fo_idx = GetVariantChrFoIdx(cip, index_vidx);
    *chr_fo_idxp = chr_fo_idx;
    chr_end_vidx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  }
  uint32_t index_variant_bp = variant_bps[index_vidx];
  const uint32_t chr_start_vidx = cip->chr_fo_vidx_start[chr_fo_idx];

  // Locate the left end of the island.  This is straightforward since we
  // already have the leftmost index variant.
  uint32_t search_bp = (index_variant_bp < bp_radius)? 0 : (index_variant_bp - bp_radius);
  const uint32_t first_vidx_in_range = LowerBoundConstrainedNonemptyU32(variant_bps, chr_start_vidx, chr_end_vidx, search_bp);
  vidx_start = AdvTo1Bit(observed_variants, first_vidx_in_range);

  // Now locate the right end of the island.  This requires either
  // - identifying a pair of index variants that are more than bp_radius bp
  //   apart from each other, which also don't share a non-index variant
  //   candidate; or
  // - reaching the end of the chromosome without observing such a gap.
  uint32_t last_known_bp = index_variant_bp;
  for (uint32_t vidx_end = index_vidx + 1; ; ++vidx_end) {
    // Jump (bp_radius + 1) bp to the right, and identify the index variants
    // flanking that spot.  (We know the end of the island cannot come before
    // this boundary.)
    vidx_end = ExpsearchU32(variant_bps, vidx_end, chr_end_vidx, last_known_bp + bp_radius + 1);
    const uint32_t last_index_vidx = FindLast1BitBefore(icandidate_vbitvec, vidx_end);
    const uint32_t last_index_bp = variant_bps[last_index_vidx];
    vidx_end = AdvBoundedTo1Bit(icandidate_vbitvec, vidx_end, chr_end_vidx);
    // Now vidx_end is either the first possibly-out-of-island index variant,
    // or equal to chr_end_vidx (which indicates end-of-chromosome).
    if (vidx_end < chr_end_vidx) {
      const uint32_t later_index_bp = variant_bps[vidx_end];
      const uint32_t delta = later_index_bp - last_index_bp;
      if (delta <= bp_radius) {
        // Obviously no gap here.
        last_known_bp = later_index_bp;
        continue;
      }
    }
    // Either the next index variant is separated by > bp_radius bp from the
    // last confirmed on-island index variant, or no such next index variant
    // exists at all.
    // In both cases, we want to locate the rightmost index-or-non-index
    // variant in range of the last confirmed on-island index variant.
    const uint32_t search_start_vidx = last_index_vidx + 1;
    const uint32_t first_out_of_range_vidx = ExpsearchU32(variant_bps, search_start_vidx, chr_end_vidx, last_index_bp + bp_radius + 1);
    const uint32_t last_in_range_vidx = FindLast1BitBefore(observed_variants, first_out_of_range_vidx);
    if ((vidx_end == chr_end_vidx) || (variant_bps[vidx_end] - variant_bps[last_in_range_vidx] <= bp_radius)) {
      *vidx_endp = last_in_range_vidx + 1;
      return vidx_start;
    }
    last_known_bp = variant_bps[vidx_end];
  }
}

// allele_idx_first and _last, instead of _start and _end, because some
// allele_idxs are skipped.
uint32_t GetNextIslandIdxs(const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const uintptr_t* icandidate_vbitvec, const uintptr_t* observed_variants, const uintptr_t* observed_alleles, const uintptr_t* observed_alleles_cumulative_popcounts_w, uint32_t raw_variant_ct, uint32_t bp_radius, uintptr_t* oaidx_startp, uintptr_t* oaidx_endp, uint32_t* vidx_endp, uint32_t* chr_fo_idxp, uintptr_t* allele_idx_firstp, uintptr_t* allele_idx_lastp) {
  uint32_t vidx_start = GetNextIsland(cip, variant_bps, icandidate_vbitvec, observed_variants, raw_variant_ct, bp_radius, vidx_endp, chr_fo_idxp);
  if (vidx_start == raw_variant_ct) {
    return raw_variant_ct;
  }
  const uint32_t vidx_end = *vidx_endp;
  uintptr_t allele_idx_first;
  uintptr_t allele_idx_last;
  if (!allele_idx_offsets) {
    allele_idx_first = vidx_start * 2;
    allele_idx_last = vidx_end * 2;
  } else {
    allele_idx_first = allele_idx_offsets[vidx_start];
    allele_idx_last = allele_idx_offsets[vidx_end];
  }
  allele_idx_first = AdvTo1Bit(observed_alleles, allele_idx_first);
  allele_idx_last = FindLast1BitBefore(observed_alleles, allele_idx_last);
  *allele_idx_firstp = allele_idx_first;
  *allele_idx_lastp = allele_idx_last;
  if (oaidx_startp) {
    *oaidx_startp = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, allele_idx_first);
  }
  *oaidx_endp = 1 + RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, allele_idx_last);
  return vidx_start;
}

void ScanPhaseDosage(const uintptr_t* observed_variants, PgenFileInfo* pgfip, uint32_t vidx_start, uint32_t vidx_end, uint32_t check_phase, uint32_t check_dosage, uint32_t* load_phasep, uint32_t* load_dosagep) {
  uint32_t vrtype_mask = 0;
  if (*load_phasep != check_phase) {
    vrtype_mask |= 0x90;
  }
  if (*load_dosagep != check_dosage) {
    vrtype_mask |= 0x60;
  }
  if (!vrtype_mask) {
    return;
  }
  // I suspect this is faster than calling AdvBoundedTo1Bit in the inner loop?
  // Todo: benchmark.
  const uint32_t vidx_ct = PopcountBitRange(observed_variants, vidx_start, vidx_end);
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(observed_variants, vidx_start, &variant_uidx_base, &cur_bits);
  for (uint32_t uii = 0; uii != vidx_ct; ++uii) {
    const uintptr_t variant_uidx = BitIter1(observed_variants, &variant_uidx_base, &cur_bits);
    const uint32_t vrtype = GetPgfiVrtype(pgfip, variant_uidx);
    const uint32_t vrtype_masked = vrtype & vrtype_mask;
    if (vrtype_masked) {
      if (vrtype_masked & 0x90) {
        *load_phasep = 1;
        vrtype_mask &= ~0x90;
      }
      if (vrtype_masked & 0x60) {
        *load_dosagep = 1;
        vrtype_mask &= ~0x60;
      }
      if (!vrtype_mask) {
        return;
      }
    }
  }
}

typedef struct ClumpEntryStruct {
  struct ClumpEntryStruct* next;
  uint32_t pval_bin_x2; // low bit indicates >= --clump-p2 threshold
  uint32_t file_idx1_x2; // low bit indicates A1=ALT in force_a1 case
} ClumpEntry;

typedef struct ClumpPvalStruct {
  double ln_pval;
  uintptr_t allele_idx;
#ifdef __cplusplus
  bool operator<(const struct ClumpPvalStruct& rhs) const {
    if (ln_pval != rhs.ln_pval) {
      return (ln_pval < rhs.ln_pval);
    }
    return allele_idx < rhs.allele_idx;
  }
#endif
} ClumpPval;

#ifndef __cplusplus
int32_t ClumpPvalCmp(const void* aa, const void* bb) {
  const ClumpPval* cp1 = S_CAST(const ClumpPval*, aa);
  const ClumpPval* cp2 = S_CAST(const ClumpPval*, bb);
  const double ln_pval1 = cp1->ln_pval;
  const double ln_pval2 = cp2->ln_pval;
  if (ln_pval1 != ln_pval2) {
    return (ln_pval1 < ln_pval2)? -1 : 1;
  }
  return (cp1->allele_idx < cp2->allele_idx)? -1 : 1;
}
#endif

ENUM_U31_DEF_START()
  kClumpJobNone,
  kClumpJobHighmemUnpack,
  kClumpJobHighmemR2,
  kClumpJobLowmemR2
ENUM_U31_DEF_END(ClumpJobType);

// Unpacked representations:
// 1. If not loading dosage:
//    a. If not loading phase:
//         uint32_t is_sparse, nm_ct, nmaj_ct, ssq:
//           RoundUpPow2(16, kBytesPerVec)
//         {one_bitvec, two_bitvec, nm_bitvec}: founder_ctaw words each, total
//           3 * bitvec_byte_ct
//         uint32_t x_male_nm_ct, x_male_nmaj_ct, x_male_ssq: if chrX, another
//           RoundUpPow2(12, kBytesPerVec)
//         (sparse representation usually smaller, but need to compare if not
//         AVX2)
//       Male-specific values placed in the back so that the unpacker doesn't
//       need to distinguish between no-x-male-stats and
//       x-male-stats-not-needed.
//    b. If loading phase, also need phasepresent and phaseinfo, founder_ctaw
//       words each, for a total of 5 * bitvec_byte_ct +
//       RoundUpPow2(16, kBytesPerVec) + RoundUpPow2(0|12, kBytesPerVec)
// 2. If loading dosage:
//    a. If not loading phase:
//         {dosage_vec, dosage_het}: founder_dosagev_ct vectors each; plus
//           nm_bitvec, i.e. 2 * dosagevec_byte_ct + bitvec_byte_ct
//         uint64_t nmaj_dosage, nmaj_dosage_ssq?, uint32_t nm_ct,
//           x_male_nm_ct?, uint64_t x_male_nmaj_dosage?,
//           x_male_nmaj_dosage_ssq?: RoundUpPow2({12|20|24|40}, kBytesPerVec)
//    b. If loading phase, also need main_dphase_deltas, for a total of
//       3 * dosagevec_byte_ct + bitvec_byte_ct + RoundUpPow2({12|20|24|40},
//       kBytesPerVec)
//
// Possible todo: Unpack the index variant in a way that allows r^2 to be
// computed efficiently against other variants *without* unpacking them.

// In the HighmemUnpack step, the main thread wants to prepare the next copy of
// this while the worker threads are processing the current copy.
typedef struct ClumpCtxAlternatingStruct {
  uintptr_t* oaidx_starts;

  unsigned char padding[kCacheline];
} ClumpCtxAlternating;

typedef struct ClumpCtxStruct {
  // Shared constants.
  const uintptr_t* observed_variants;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* variant_last_alidxs;
  const uint32_t* variant_last_alidxs_cumulative_popcounts;
  const uintptr_t* observed_alleles;
  const uintptr_t* observed_alleles_cumulative_popcounts_w;
  const uintptr_t* founder_info;
  const uint32_t* founder_info_cumulative_popcounts;
  const uintptr_t* founder_male_collapsed;
  const Dosage* male_dosage_invmask;
  const uintptr_t* founder_nonmale_collapsed;
  const Dosage* nonmale_dosage_invmask;
  const uintptr_t* founder_female_collapsed;
  const uintptr_t* founder_female_collapsed_interleaved;
  uint32_t founder_ct;
  uint32_t founder_male_ct;
  // Precomputed for convenience.
  uintptr_t pgv_byte_stride;
  uintptr_t bitvec_byte_ct;
  uintptr_t dosagevec_byte_ct;
  double r2_thresh;
  unsigned char allow_overlap;

  // Remaining non-index variants.
  uintptr_t* candidate_oabitvec;

  // In highmem mode, this is the base of (calc_thread_ct + 1) PgenVariant
  // buffer groups, offset by pgv_byte_stride bytes from each other.  They are
  // used by worker-thread PgrGetDp() calls.
  // In lowmem mode, this is the base of at least cur_nonindex_ct PgenVariant
  // buffer groups, offset by pgv_byte_stride bytes from each other.  They are
  // filled by main-thread PgrGetDp() calls, and interpreted by the worker
  // threads.
  PgenVariant pgv_base;

  // In lowmem mode, these store additional PgenVariant fields.
  uint32_t* phasepresent_cts;
  uint32_t* dosage_cts;
  uint32_t* dphase_cts;

  // Other per-thread resources.

  // calc_thread_ct of these.
  PgenReader** pgr_ptrs;
  // Highmem: HighmemUnpack unpacks all current-island(-group) variants to this
  //          memory region, then HighmemR2 reads from here.  Could have
  //          thousands, or even millions of variants here.
  // LowmemR2: each thread just needs workspace for one variant here.  Index
  //           variant is stored first.
  // Unpacked variant representation can vary based on phase_type and
  // load_dosage.
  unsigned char* unpacked_variants;

  // In unphased case, one of these per thread to enable sparse-optimization.
  uintptr_t** raregeno_bufs;
  uint32_t** difflist_sample_id_bufs;
  // chrX workspaces
  //   nm_buf: max(founder_nonmale_ct, founder_male_ct) bits
  //   invmask_buf: max(founder_nonmale_ct, founder_male_ct) uint16s
  uintptr_t** chrx_workspaces;
  // Precomputed for convenience.
  uintptr_t unpacked_byte_stride;

  // Current island-group.
  uintptr_t igroup_oaidx_start;
  uintptr_t allele_widx_start;
  uintptr_t allele_widx_end;
  unsigned char is_x;
  unsigned char is_y;
  unsigned char phase_type;
  unsigned char load_dosage;
  ClumpJobType job_type;

  uintptr_t index_oaidx_offset;
  uintptr_t cur_nonindex_ct;

  ClumpCtxAlternating a[2];

  // Result buffers.
  // (calc_thread_ct + 1) of these, since main thread joins in
  uintptr_t** ld_idx_found;

  uint64_t err_info;
} ClumpCtx;

// phase_type and load_dosage are not read from ctx, since they can change as
// we try to extend an island group along a single chromosome.
uintptr_t UnpackedByteStride(const ClumpCtx* ctx, R2PhaseType phase_type, uint32_t x_exists, uint32_t load_dosage) {
  const uintptr_t bitvec_byte_ct = ctx->bitvec_byte_ct;
  if (!load_dosage) {
    uintptr_t stride = RoundUpPow2(16, kBytesPerVec) + (3 + 2 * (phase_type == kR2PhaseTypePresent)) * bitvec_byte_ct;
#ifndef USE_AVX2
    const uint32_t founder_ct = ctx->founder_ct;
    const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
    const uint32_t max_difflist_len = founder_ct / 64;
    const uintptr_t sparse_req = RoundUpPow2((6 + max_difflist_len) * sizeof(int32_t), kBytesPerVec) + NypCtToVecCt(max_difflist_len) * kBytesPerVec + RoundUpPow2(founder_ctl * (kBytesPerWord + sizeof(int32_t)), kBytesPerVec);
    if (sparse_req > stride) {
      stride = sparse_req;
    }
#endif
    if (x_exists) {
      stride += RoundUpPow2((2 + (phase_type == kR2PhaseTypeUnphased)) * sizeof(int32_t), kBytesPerVec);
    }
    return stride;
  }
  const uintptr_t trail_byte_ct = RoundUpPow2(((1 + (phase_type == kR2PhaseTypeUnphased)) * sizeof(int64_t) + sizeof(int32_t)) << x_exists, kBytesPerVec);
  return (1 + phase_type) * ctx->dosagevec_byte_ct + bitvec_byte_ct + trail_byte_ct;
}

// does not check multiallelic fields
void ClumpPgenVariantIncr(uintptr_t byte_ct, PgenVariant* pgvp) {
  const uintptr_t word_ct = byte_ct / kBytesPerWord;
  pgvp->genovec = &(pgvp->genovec[word_ct]);
  if (pgvp->phasepresent) {
    pgvp->phasepresent = &(pgvp->phasepresent[word_ct]);
    pgvp->phaseinfo = &(pgvp->phaseinfo[word_ct]);
  }
  if (pgvp->dosage_present) {
    const uintptr_t u16_ct = byte_ct / sizeof(int16_t);
    pgvp->dosage_present = &(pgvp->dosage_present[word_ct]);
    pgvp->dosage_main = &(pgvp->dosage_main[u16_ct]);
    if (pgvp->dphase_present) {
      pgvp->dphase_present = &(pgvp->dphase_present[word_ct]);
      pgvp->dphase_delta = &(pgvp->dphase_delta[u16_ct]);
    }
  }
}

// x_male_nmaj_ct, x_male_ssq, and x_male_nm_ct only filled if
// founder_male_collapsed is non-null.
void LdUnpackNondosageDense(const PgenVariant* pgvp, const uintptr_t* founder_male_collapsed, uint32_t sample_ct, R2PhaseType phase_type, unsigned char* dst_iter) {
  const uintptr_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  uint32_t* dst_u32 = R_CAST(uint32_t*, dst_iter);
  dst_u32[0] = 0; // is_sparse
  uintptr_t* dst_witer = R_CAST(uintptr_t*, &(dst_iter[RoundUpPow2(16, kBytesPerVec)]));
  uintptr_t* one_bitvec = dst_witer;
  dst_witer = &(dst_witer[sample_ctaw]);
  uintptr_t* two_bitvec = dst_witer;
  dst_witer = &(dst_witer[sample_ctaw]);
  uintptr_t* nm_bitvec = dst_witer;
  dst_witer = &(dst_witer[sample_ctaw]);
  uintptr_t* phasepresent = nullptr;
  uintptr_t* phaseinfo = nullptr;
  if (phase_type == kR2PhaseTypePresent) {
    phasepresent = dst_witer;
    dst_witer = &(dst_witer[sample_ctaw]);
    phaseinfo = dst_witer;
    dst_witer = &(dst_witer[sample_ctaw]);
  }
  GenoarrSplit12Nm(pgvp->genovec, sample_ct, one_bitvec, two_bitvec, nm_bitvec);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  dst_u32[1] = PopcountWords(nm_bitvec, sample_ctl);
  dst_u32[2] = GenoBitvecSum(one_bitvec, two_bitvec, sample_ctl);
  dst_u32[3] = 0; // defensive
  if (phase_type == kR2PhaseTypePresent) {
    if (!pgvp->phasepresent_ct) {
      ZeroWArr(sample_ctl, phasepresent);
    } else {
      memcpy(phasepresent, pgvp->phasepresent, sample_ctl * sizeof(intptr_t));
      memcpy(phaseinfo, pgvp->phaseinfo, sample_ctl * sizeof(intptr_t));
    }
  } else if (phase_type == kR2PhaseTypeUnphased) {
    dst_u32[3] = dst_u32[2] + 2 * PopcountWords(two_bitvec, sample_ctl);
  }
  if (founder_male_collapsed) {
    uint32_t* dst_x_u32 = R_CAST(uint32_t*, dst_witer);
    // x_male_nm_ct
    dst_x_u32[0] = PopcountWordsIntersect(founder_male_collapsed, nm_bitvec, sample_ctl);
    dst_x_u32[1] = GenoBitvecSumSubset(founder_male_collapsed, one_bitvec, two_bitvec, sample_ctl);
    if (phase_type == kR2PhaseTypeUnphased) {
      // x_male_ssq
      dst_x_u32[2] = dst_x_u32[1] + 2 * PopcountWordsIntersect(founder_male_collapsed, two_bitvec, sample_ctl);
    }
  }
}

// Ok if trailing bits of raregeno aren't clear.
void LdUnpackNondosageSparse(const uintptr_t* raregeno, const uint32_t* difflist_sample_ids, const uintptr_t* founder_male_collapsed, uint32_t sample_ct, uint32_t male_ct, uint32_t difflist_common_geno, uint32_t difflist_len, unsigned char* dst_iter) {
  uint32_t* dst_u32 = R_CAST(uint32_t*, dst_iter);
  dst_u32[0] = 1; // is_sparse
  dst_u32[4] = difflist_common_geno;
  dst_u32[5] = difflist_len;
  uint32_t nm_ct = sample_ct;
  uint32_t nmaj_ct = 0;
  uint32_t ssq = 0;
  if (difflist_common_geno) {
    if (difflist_common_geno == 3) {
      nm_ct = difflist_len;
    } else {
      const uint32_t two_ct = sample_ct - difflist_len;
      nmaj_ct = 2 * two_ct;
      ssq = 4 * two_ct;
    }
  }
  uintptr_t* raregeno_dst = R_CAST(uintptr_t*, &(dst_iter[RoundUpPow2((6 + difflist_len) * sizeof(int32_t), kBytesPerVec)]));
  const uint32_t raregeno_word_ct = NypCtToWordCt(difflist_len);
  uintptr_t* difflist_include_dst = &(raregeno_dst[RoundUpPow2(raregeno_word_ct, kWordsPerVec)]);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  ZeroWArr(sample_ctl, difflist_include_dst);
  uint32_t* difflist_include_cumulative_popcounts_dst = R_CAST(uint32_t*, &(difflist_include_dst[sample_ctl]));
  if (difflist_len) {
    uint32_t* difflist_sample_ids_dst = &(dst_u32[6]);
    for (uint32_t uii = 0; uii != difflist_len; ++uii) {
      const uint32_t sample_idx = difflist_sample_ids[uii];
      difflist_sample_ids_dst[uii] = sample_idx;
      SetBit(sample_idx, difflist_include_dst);
    }
    FillCumulativePopcounts(difflist_include_dst, sample_ctl, difflist_include_cumulative_popcounts_dst);
    memcpy(raregeno_dst, raregeno, raregeno_word_ct * sizeof(intptr_t));
    ZeroTrailingNyps(difflist_len, raregeno_dst);
    STD_ARRAY_DECL(uint32_t, 4, genocounts);
    GenoarrCountFreqsUnsafe(raregeno_dst, difflist_len, genocounts);
    nm_ct -= genocounts[3];
    nmaj_ct += genocounts[1] + 2 * genocounts[2];
    ssq += genocounts[1] + 4 * genocounts[2];
  } else {
    ZeroU32Arr(sample_ctl, difflist_include_cumulative_popcounts_dst);
  }
  dst_u32[1] = nm_ct;
  dst_u32[2] = nmaj_ct;
  dst_u32[3] = ssq;
  if (founder_male_collapsed) {
    uint32_t x_male_nm_ct = male_ct;
    uint32_t x_male_nmaj_ct = 0;
    uint32_t x_male_ssq = 0;
    if (difflist_common_geno) {
      if (difflist_common_geno == 3) {
        x_male_nm_ct = 0;
      } else {
        x_male_nmaj_ct = male_ct * 2;
        x_male_ssq = male_ct * 4;
      }
    }
    if (difflist_len) {
      const uint32_t word_ct_m1 = raregeno_word_ct - 1;
      uint32_t genocounts[4];
      ZeroU32Arr(4, genocounts);
      uint32_t loop_len = kBitsPerWordD2;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(difflist_len, kBitsPerWordD2);
        }
        const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
        uintptr_t raregeno_word = raregeno[widx];
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          const uint32_t sample_idx = cur_difflist_sample_ids[uii];
          if (IsSet(founder_male_collapsed, sample_idx)) {
            const uintptr_t cur_geno = raregeno_word & 3;
            genocounts[cur_geno] += 1;
          }
          raregeno_word = raregeno_word >> 2;
        }
      }
      if (difflist_common_geno != 3) {
        x_male_nm_ct -= genocounts[3];
      } else {
        x_male_nm_ct = genocounts[0] + genocounts[1] + genocounts[2];
      }
      if (difflist_common_geno != 2) {
        x_male_nmaj_ct = genocounts[1] + 2 * genocounts[2];
        x_male_ssq = genocounts[1] + 4 * genocounts[2];
      } else {
        const uint32_t zero_or_missing_ct = genocounts[0] + genocounts[3];
        x_male_nmaj_ct -= genocounts[1] + 2 * zero_or_missing_ct;
        x_male_ssq -= 3 * genocounts[1] + 4 * zero_or_missing_ct;
      }
    }
    uint32_t* dst_x_u32 = &(difflist_include_cumulative_popcounts_dst[sample_ctl]);
    dst_x_u32[0] = x_male_nm_ct;
    dst_x_u32[1] = x_male_nmaj_ct;
    dst_x_u32[2] = x_male_ssq;
  }
}

static inline uint32_t LdNondosageTrailAlignedByteCt(R2PhaseType phase_type, uint32_t x_exists) {
  return x_exists * RoundUpPow2((2 + (phase_type == kR2PhaseTypeUnphased)) * sizeof(int32_t), kBytesPerVec);
}

static inline uint32_t LdDosageTrailAlignedByteCt(R2PhaseType phase_type, uint32_t x_exists) {
  return RoundUpPow2((sizeof(int64_t) + sizeof(int32_t) + sizeof(int64_t) * (phase_type == kR2PhaseTypeUnphased)) << x_exists, kBytesPerVec);
}

void LdUnpackDosage(const PgenVariant* pgvp, const uintptr_t* founder_male_collapsed, const Dosage* male_dosage_invmask, uint32_t sample_ct, R2PhaseType phase_type, unsigned char* dst_iter) {
  const uintptr_t dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  const uintptr_t dosagevec_byte_ct = dosagev_ct * kBytesPerVec;
  Dosage* dosage_vec = R_CAST(Dosage*, dst_iter);
  dst_iter = &(dst_iter[dosagevec_byte_ct]);
  Dosage* dosage_het = nullptr;
  if (phase_type != kR2PhaseTypeUnphased) {
    dosage_het = R_CAST(Dosage*, dst_iter);
    dst_iter = &(dst_iter[dosagevec_byte_ct]);
  }
  uintptr_t* nm_bitvec = R_CAST(uintptr_t*, dst_iter);
  const uintptr_t bitvec_byte_ct = BitCtToVecCt(sample_ct) * kBytesPerVec;
  dst_iter = &(dst_iter[bitvec_byte_ct]);
  SDosage* dense_dphase_delta = nullptr;
  if (phase_type == kR2PhaseTypePresent) {
    dense_dphase_delta = R_CAST(SDosage*, dst_iter);
    dst_iter = &(dst_iter[dosagevec_byte_ct]);
  }
  // In 32-bit build, no alignment guarantee for uint64s.
  unsigned char* nmaj_dosage_uc_ptr = dst_iter;
  dst_iter = &(dst_iter[sizeof(int64_t)]);
  unsigned char* nmaj_dosage_ssq_uc_ptr = nullptr;
  if (phase_type == kR2PhaseTypeUnphased) {
    nmaj_dosage_ssq_uc_ptr = dst_iter;
    dst_iter = &(dst_iter[sizeof(int64_t)]);
  }
  uint32_t* nm_ct_ptr = R_CAST(uint32_t*, dst_iter);

  PopulateDenseDosage(pgvp->genovec, pgvp->dosage_present, pgvp->dosage_main, sample_ct, pgvp->dosage_ct, dosage_vec);
  const uint64_t nmaj_dosage = DenseDosageSum(dosage_vec, dosagev_ct);
  memcpy(nmaj_dosage_uc_ptr, &nmaj_dosage, sizeof(int64_t));
  GenoarrToNonmissing(pgvp->genovec, sample_ct, nm_bitvec);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  BitvecOr(pgvp->dosage_present, sample_ctl, nm_bitvec);
  *nm_ct_ptr = PopcountWords(nm_bitvec, sample_ctl);
  if (phase_type == kR2PhaseTypeUnphased) {
    const uint64_t nmaj_dosage_ssq = DosageUnsignedDotprod(dosage_vec, dosage_vec, dosagev_ct);
    memcpy(nmaj_dosage_ssq_uc_ptr, &nmaj_dosage_ssq, sizeof(int64_t));
  } else {
    FillDosageHet(dosage_vec, dosagev_ct, dosage_het);
    if (phase_type == kR2PhaseTypePresent) {
      PopulateDenseDphase(pgvp->phasepresent, pgvp->phaseinfo, pgvp->dosage_present, dosage_vec, pgvp->dphase_present, pgvp->dphase_delta, sample_ct, pgvp->phasepresent_ct, pgvp->dosage_ct, pgvp->dphase_ct, dense_dphase_delta);
    }
  }
  if (founder_male_collapsed == nullptr) {
    return;
  }
  dst_iter = &(dst_iter[sizeof(int32_t)]);
  uint32_t* x_male_nm_ct_ptr = R_CAST(uint32_t*, dst_iter);
  *x_male_nm_ct_ptr = PopcountWordsIntersect(founder_male_collapsed, nm_bitvec, sample_ctl);
  dst_iter = &(dst_iter[sizeof(int32_t)]);
  const uint64_t x_male_nmaj_dosage = DenseDosageSumSubset(dosage_vec, male_dosage_invmask, dosagev_ct);
  memcpy(dst_iter, &x_male_nmaj_dosage, sizeof(int64_t));
  if (phase_type == kR2PhaseTypeUnphased) {
    dst_iter = &(dst_iter[sizeof(int64_t)]);
    const uint64_t x_male_nmaj_dosage_ssq = DosageUnsignedDotprodSubset(male_dosage_invmask, dosage_vec, dosage_vec, dosagev_ct);
    memcpy(dst_iter, &x_male_nmaj_dosage_ssq, sizeof(int64_t));
  }
}

typedef struct R2NondosageDenseStruct {
  const uintptr_t* one_bitvec;
  const uintptr_t* two_bitvec;
  const uintptr_t* nm_bitvec;
  const uintptr_t* phasepresent; // may be uninitialized
  const uintptr_t* phaseinfo; // may be uninitialized
} R2NondosageDense;

typedef struct R2NondosageSparseStruct {
  uint32_t difflist_common_geno;
  uint32_t difflist_len;
  const uint32_t* difflist_sample_ids;
  const uintptr_t* raregeno;
  const uintptr_t* difflist_include;
  const uint32_t* difflist_include_cumulative_popcounts;
  // probable todo: support phase
} R2NondosageSparse;

typedef union {
  R2NondosageDense d;
  R2NondosageSparse s;
} R2NondosagePayload;

typedef struct R2NondosageVariantStruct {
  uint32_t is_sparse;
  uint32_t nm_ct;
  uint32_t nmaj_ct;
  uint32_t ssq; // set to 0 unless unphased calc
  R2NondosagePayload p;
  uint32_t x_male_nm_ct; // may be uninitialized
  uint32_t x_male_nmaj_ct; // may be uninitialized
  uint32_t x_male_ssq; // may be uninitialized
} R2NondosageVariant;

typedef struct R2DosageVariantStruct {
  const Dosage* dosage_vec;
  const Dosage* dosage_het; // may be uninitialized
  const uintptr_t* nm_bitvec;
  // probable todo: replace dense_dphase_delta with dosage_vec2 + dosage_het2.
  // takes ~33% more space, but I'd be shocked if it didn't result in a
  // noticeable phased-dosage r^2 speedup since computing the two het-sides
  // from scratch is so annoying.  (and if that's how we're handling
  // dosage_het, we may as well handle dosage_vec the same way even though it
  // needs it less.)
  const SDosage* dense_dphase_delta; // may be uninitialized
  uint64_t nmaj_dosage;
  uint64_t nmaj_dosage_ssq; // may be uninitialized
  uint32_t nm_ct;
  uint32_t x_male_nm_ct; // may be uninitialized
  uint64_t x_male_nmaj_dosage; // may be uninitialized
  uint64_t x_male_nmaj_dosage_ssq; // may be uninitialized
} R2DosageVariant;

typedef union {
  R2NondosageVariant nd;
  R2DosageVariant d;
} R2Variant;

void FillR2Nondosage(const unsigned char* src_iter, uint32_t sample_ct, R2PhaseType phase_type, uint32_t is_x, R2NondosageVariant* ndp) {
  // See LdUnpackNondosage{Dense,Sparse}().
  const uint32_t* src_u32 = R_CAST(const uint32_t*, src_iter);
  const uint32_t is_sparse = src_u32[0];
  ndp->is_sparse = is_sparse;
  ndp->nm_ct = src_u32[1];
  ndp->nmaj_ct = src_u32[2];
  ndp->ssq = src_u32[3];
  if (is_sparse) {
    R2NondosageSparse* ndsp = &(ndp->p.s);
    ndsp->difflist_common_geno = src_u32[4];
    const uint32_t difflist_len = src_u32[5];
    ndsp->difflist_len = difflist_len;
    ndsp->difflist_sample_ids = &(src_u32[6]);
    const uintptr_t* raregeno = R_CAST(const uintptr_t*, &(src_iter[RoundUpPow2(sizeof(int32_t) * (6 + difflist_len), kBytesPerVec)]));
    ndsp->raregeno = raregeno;
    const uint32_t raregeno_word_ct = NypCtToWordCt(difflist_len);
    const uintptr_t* difflist_include = &(raregeno[RoundUpPow2(raregeno_word_ct, kWordsPerVec)]);
    ndsp->difflist_include = difflist_include;
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t* difflist_include_cumulative_popcounts = R_CAST(const uint32_t*, &(difflist_include[sample_ctl]));
    ndsp->difflist_include_cumulative_popcounts = difflist_include_cumulative_popcounts;
    if (is_x) {
      const uint32_t* src_x_u32 = &(difflist_include_cumulative_popcounts[sample_ctl]);
      ndp->x_male_nm_ct = src_x_u32[0];
      ndp->x_male_nmaj_ct = src_x_u32[1];
      ndp->x_male_ssq = src_x_u32[2];
    }
    return;
  }
  const uintptr_t* src_witer = R_CAST(const uintptr_t*, &(src_iter[RoundUpPow2(16, kBytesPerVec)]));
  const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  R2NondosageDense* nddp = &(ndp->p.d);
  nddp->one_bitvec = src_witer;
  src_witer = &(src_witer[sample_ctaw]);
  nddp->two_bitvec = src_witer;
  src_witer = &(src_witer[sample_ctaw]);
  nddp->nm_bitvec = src_witer;
  src_witer = &(src_witer[sample_ctaw]);
  if (phase_type == kR2PhaseTypePresent) {
    nddp->phasepresent = src_witer;
    src_witer = &(src_witer[sample_ctaw]);
    nddp->phaseinfo = src_witer;
    src_witer = &(src_witer[sample_ctaw]);
  }
  if (is_x) {
    const uint32_t* src_x_u32 = R_CAST(const uint32_t*, src_witer);
    ndp->x_male_nm_ct = src_x_u32[0];
    ndp->x_male_nmaj_ct = src_x_u32[1];
    if (phase_type == kR2PhaseTypeUnphased) {
      ndp->x_male_ssq = src_x_u32[2];
    }
  }
}

void FillR2Dosage(const unsigned char* src_iter, uint32_t sample_ct, R2PhaseType phase_type, uint32_t is_x, R2DosageVariant* dp) {
  // See LdUnpackDosage().
  const uintptr_t dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  const uintptr_t dosagevec_byte_ct = dosagev_ct * kBytesPerVec;
  dp->dosage_vec = R_CAST(const Dosage*, src_iter);
  src_iter = &(src_iter[dosagevec_byte_ct]);
  if (phase_type != kR2PhaseTypeUnphased) {
    dp->dosage_het = R_CAST(const Dosage*, src_iter);
    src_iter = &(src_iter[dosagevec_byte_ct]);
  }
  dp->nm_bitvec = R_CAST(const uintptr_t*, src_iter);
  const uintptr_t bitvec_byte_ct = BitCtToVecCt(sample_ct) * kBytesPerVec;
  src_iter = &(src_iter[bitvec_byte_ct]);
  if (phase_type == kR2PhaseTypePresent) {
    dp->dense_dphase_delta = R_CAST(const SDosage*, src_iter);
    src_iter = &(src_iter[dosagevec_byte_ct]);
  }
  memcpy(&(dp->nmaj_dosage), src_iter, sizeof(int64_t));
  const unsigned char* trail_iter = &(src_iter[sizeof(int64_t)]);
  if (phase_type == kR2PhaseTypeUnphased) {
    memcpy(&(dp->nmaj_dosage_ssq), trail_iter, sizeof(int64_t));
    trail_iter = &(trail_iter[sizeof(int64_t)]);
  }
  memcpy(&(dp->nm_ct), trail_iter, sizeof(int32_t));
  if (!is_x) {
    return;
  }
  trail_iter = &(trail_iter[sizeof(int32_t)]);
  memcpy(&(dp->x_male_nm_ct), trail_iter, sizeof(int32_t));
  trail_iter = &(trail_iter[sizeof(int32_t)]);
  memcpy(&(dp->x_male_nmaj_dosage), trail_iter, sizeof(int64_t));
  if (phase_type == kR2PhaseTypeUnphased) {
    trail_iter = &(trail_iter[sizeof(int64_t)]);
    memcpy(&(dp->x_male_nmaj_dosage_ssq), trail_iter, sizeof(int64_t));
  }
}

void FillR2V(const unsigned char* src_iter, uint32_t sample_ct, R2PhaseType phase_type, uint32_t is_x, uint32_t load_dosage, R2Variant* r2vp) {
  if (!load_dosage) {
    FillR2Nondosage(src_iter, sample_ct, phase_type, is_x, &(r2vp->nd));
  } else {
    FillR2Dosage(src_iter, sample_ct, phase_type, is_x, &(r2vp->d));
  }
}

void ClumpHighmemUnpack(uintptr_t tidx, uint32_t parity, ClumpCtx* ctx) {
  // Unpack (variant, aidx)s to unpacked_variants.
  const uintptr_t oaidx_end = ctx->a[parity].oaidx_starts[tidx + 1];
  uintptr_t oaidx = ctx->a[parity].oaidx_starts[tidx];
  if (oaidx == oaidx_end) {
    return;
  }
  const uintptr_t* variant_last_alidxs = ctx->variant_last_alidxs;
  const uint32_t* variant_last_alidxs_cumulative_popcounts = ctx->variant_last_alidxs_cumulative_popcounts;
  const uintptr_t* observed_alleles = ctx->observed_alleles;
  const uintptr_t* observed_alleles_cumulative_popcounts_w = ctx->observed_alleles_cumulative_popcounts_w;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uintptr_t pgv_byte_stride = ctx->pgv_byte_stride;
  PgenVariant pgv = ctx->pgv_base;
  ClumpPgenVariantIncr(pgv_byte_stride * tidx, &pgv);
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  const uintptr_t unpacked_byte_stride = ctx->unpacked_byte_stride;
  unsigned char* write_iter;
  {
    const uintptr_t igroup_oaidx_start = ctx->igroup_oaidx_start;
    write_iter = &(ctx->unpacked_variants[(oaidx - igroup_oaidx_start) * unpacked_byte_stride]);
  }
  const uintptr_t* founder_info = ctx->founder_info;
  const uint32_t* founder_info_cumulative_popcounts = ctx->founder_info_cumulative_popcounts;
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t is_x = ctx->is_x;
  const uintptr_t* founder_male_collapsed = is_x? ctx->founder_male_collapsed : nullptr;
  const Dosage* male_dosage_invmask = is_x? ctx->male_dosage_invmask : nullptr;
  const uintptr_t* founder_female_collapsed = ctx->founder_female_collapsed;
  const uintptr_t* founder_female_collapsed_interleaved = ctx->founder_female_collapsed_interleaved;
  const uint32_t founder_ct = ctx->founder_ct;
  const uint32_t founder_ctv2 = NypCtToVecCt(founder_ct);
  const uint32_t is_y = ctx->is_y;
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(founder_info_cumulative_popcounts, pgrp, &pssi);
  const R2PhaseType phase_type = S_CAST(R2PhaseType, ctx->phase_type);
  const uint32_t load_dosage = ctx->load_dosage;
  const uintptr_t allele_idx_start = IdxToUidxW(observed_alleles, observed_alleles_cumulative_popcounts_w, ctx->allele_widx_start, ctx->allele_widx_end, oaidx);
  const uint32_t max_difflist_len = founder_ct / 64;
  uintptr_t* raregeno = nullptr;
  uint32_t* difflist_sample_ids = nullptr;
  if (phase_type == kR2PhaseTypeUnphased) {
    raregeno = ctx->raregeno_bufs[tidx];
    difflist_sample_ids = ctx->difflist_sample_id_bufs[tidx];
  }
  uintptr_t allele_idx_base;
  uintptr_t cur_bits;
  BitIter1Start(observed_alleles, allele_idx_start, &allele_idx_base, &cur_bits);
  uintptr_t variant_uidx;
  PglErr reterr;
  for (; oaidx != oaidx_end; ++oaidx, write_iter = &(write_iter[unpacked_byte_stride])) {
    const uintptr_t allele_idx = BitIter1(observed_alleles, &allele_idx_base, &cur_bits);
    AlleleCode aidx;
    if (!allele_idx_offsets) {
      variant_uidx = allele_idx / 2;
      aidx = allele_idx % 2;
    } else {
      variant_uidx = RawToSubsettedPos(variant_last_alidxs, variant_last_alidxs_cumulative_popcounts, allele_idx);
      aidx = allele_idx - allele_idx_offsets[variant_uidx];
    }
    if (load_dosage) {
      if (phase_type == kR2PhaseTypePresent) {
        reterr = PgrGetInv1Dp(founder_info, pssi, founder_ct, variant_uidx, aidx, pgrp, &pgv);
      } else {
        reterr = PgrGetInv1D(founder_info, pssi, founder_ct, variant_uidx, aidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
      }
      if (unlikely(reterr)) {
        goto ClumpHighmemUnpack_err;
      }
      if (is_y) {
        InterleavedSetMissingCleardosage(founder_female_collapsed, founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec, &pgv.dosage_ct, pgv.dosage_present, pgv.dosage_main);
      }
      LdUnpackDosage(&pgv, founder_male_collapsed, male_dosage_invmask, founder_ct, phase_type, write_iter);
    } else {
      if ((phase_type == kR2PhaseTypeUnphased) && (!is_y)) {
        uint32_t difflist_common_geno;
        uint32_t difflist_len;
        reterr = PgrGetInv1DifflistOrGenovec(founder_info, pssi, founder_ct, max_difflist_len, variant_uidx, aidx, pgrp, pgv.genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
        if (unlikely(reterr)) {
          goto ClumpHighmemUnpack_err;
        }
        if (difflist_common_geno != UINT32_MAX) {
          if (difflist_len <= max_difflist_len) {
            LdUnpackNondosageSparse(raregeno, difflist_sample_ids, founder_male_collapsed, founder_ct, founder_male_ct, difflist_common_geno, difflist_len, write_iter);
            continue;
          }
          PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, founder_ct, difflist_len, pgv.genovec);
        }
      } else {
        if (phase_type == kR2PhaseTypePresent) {
          reterr = PgrGetInv1P(founder_info, pssi, founder_ct, variant_uidx, aidx, pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
        } else {
          reterr = PgrGetInv1(founder_info, pssi, founder_ct, variant_uidx, aidx, pgrp, pgv.genovec);
        }
        if (unlikely(reterr)) {
          goto ClumpHighmemUnpack_err;
        }
        if (is_y) {
          InterleavedSetMissing(founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec);
        }
      }
      LdUnpackNondosageDense(&pgv, founder_male_collapsed, founder_ct, phase_type, write_iter);
    }
  }
  return;
 ClumpHighmemUnpack_err:
  ;
  const uint64_t new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
  UpdateU64IfSmaller(new_err_info, &ctx->err_info);
}

// assumes sample_ct < 2^30.
static inline uint32_t GenoBitvecUnphasedDotprod(const uintptr_t* one_bitvec0, const uintptr_t* two_bitvec0, const uintptr_t* one_bitvec1, const uintptr_t* two_bitvec1, uint32_t word_ct) {
  uint32_t half_hom_part;
  uint32_t hethet_ct;
  GenoBitvecPhasedDotprod(one_bitvec0, two_bitvec0, one_bitvec1, two_bitvec1, word_ct, &half_hom_part, &hethet_ct);
  return 2 * half_hom_part + hethet_ct;
}

static inline uint32_t GenoBitvecUnphasedDotprodSubset(const uintptr_t* subset_mask, const uintptr_t* one_bitvec0, const uintptr_t* two_bitvec0, const uintptr_t* one_bitvec1, const uintptr_t* two_bitvec1, uint32_t word_ct) {
  uint32_t half_hom_part;
  uint32_t hethet_ct;
  GenoBitvecPhasedDotprodSubset(subset_mask, one_bitvec0, two_bitvec0, one_bitvec1, two_bitvec1, word_ct, &half_hom_part, &hethet_ct);
  return 2 * half_hom_part + hethet_ct;
}

uint32_t ComputeR2NondosageUnphased1SparseStats(const R2NondosageVariant* densevp0, const R2NondosageVariant* sparsevp1, uint32_t* nmaj_ct0_ptr, uint32_t* nmaj_ct1_ptr, uint32_t* ssq0_ptr, uint32_t* ssq1_ptr, uint32_t* dotprod_ptr) {
  const uint32_t difflist_common_geno = sparsevp1->p.s.difflist_common_geno;
  const uint32_t difflist_len = sparsevp1->p.s.difflist_len;
  const uint32_t* difflist_sample_ids = sparsevp1->p.s.difflist_sample_ids;
  uint32_t nmaj_ct0 = 0;
  uint32_t nmaj_ct1 = 0;
  uint32_t ssq0 = 0;
  uint32_t ssq1 = 0;
  uint32_t dotprod = 0;
  uint32_t valid_obs_ct = 0;
  if (difflist_common_geno != 3) {
    nmaj_ct0 = densevp0->nmaj_ct;
    ssq0 = densevp0->ssq;
    dotprod = difflist_common_geno * nmaj_ct0;
    valid_obs_ct = densevp0->nm_ct;
    nmaj_ct1 = difflist_common_geno * valid_obs_ct;
    ssq1 = difflist_common_geno * nmaj_ct1;
  }
  if (difflist_len) {
    // genovec would be more convenient than this representation here, but that
    // shouldn't be a big deal.
    const uintptr_t* nm_bitvec0 = densevp0->p.d.nm_bitvec;
    const uintptr_t* one_bitvec0 = densevp0->p.d.one_bitvec;
    const uintptr_t* two_bitvec0 = densevp0->p.d.two_bitvec;
    const uintptr_t* raregeno = sparsevp1->p.s.raregeno;
    const uint32_t word_ct_m1 = (difflist_len - 1) / kBitsPerWordD2;
    uint32_t joint_counts[16]; // low bits = dense geno, high bits = sparse
    ZeroU32Arr(16, joint_counts);
    uint32_t loop_len = kBitsPerWordD2;
    for (uint32_t widx = 0; ; ++widx) {
      if (widx >= word_ct_m1) {
        if (widx > word_ct_m1) {
          break;
        }
        loop_len = ModNz(difflist_len, kBitsPerWordD2);
      }
      const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
      uintptr_t raregeno_word = raregeno[widx];
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uintptr_t cur_sparse_geno = raregeno_word & 3;
        const uint32_t sample_idx = cur_difflist_sample_ids[uii];
        const uint32_t sample_widx = sample_idx / kBitsPerWord;
        const uint32_t sample_idx_lowbits = sample_idx % kBitsPerWord;
        const uintptr_t nm_bit = (nm_bitvec0[sample_widx] >> sample_idx_lowbits) & 1;
        const uintptr_t one_bit = (one_bitvec0[sample_widx] >> sample_idx_lowbits) & 1;
        // "& 3" takes care of this mask
        const uintptr_t two_bit_unmasked = two_bitvec0[sample_widx] >> sample_idx_lowbits;
        const uintptr_t cur_dense_geno = (nm_bit + one_bit + two_bit_unmasked * 2 - 1) & 3;
        joint_counts[cur_dense_geno + cur_sparse_geno * 4] += 1;
        raregeno_word = raregeno_word >> 2;
      }
    }
    if (difflist_common_geno != 3) {
      nmaj_ct0 -= joint_counts[13] + 2 * joint_counts[14];
      ssq0 -= joint_counts[13] + 4 * joint_counts[14];
      const uint32_t sparse_missing_dense_nm_ct = joint_counts[12] + joint_counts[13] + joint_counts[14];
      valid_obs_ct -= sparse_missing_dense_nm_ct;
    } else {
      const uint32_t dense_one_ct = joint_counts[1] + joint_counts[5] + joint_counts[9];
      const uint32_t dense_two_ct = joint_counts[2] + joint_counts[6] + joint_counts[10];
      nmaj_ct0 = dense_one_ct + 2 * dense_two_ct;
      ssq0 = dense_one_ct + 4 * dense_two_ct;
      const uint32_t dense_missing_sparse_nm_ct = joint_counts[3] + joint_counts[7] + joint_counts[11];
      valid_obs_ct = difflist_len - dense_missing_sparse_nm_ct;
    }
    const uint32_t sparse_one_ct = joint_counts[4] + joint_counts[5] + joint_counts[6];
    if (difflist_common_geno != 2) {
      const uint32_t sparse_two_ct = joint_counts[8] + joint_counts[9] + joint_counts[10];
      nmaj_ct1 = sparse_one_ct + 2 * sparse_two_ct;
      ssq1 = sparse_one_ct + 4 * sparse_two_ct;
      dotprod = joint_counts[5] + 2 * (joint_counts[6] + joint_counts[9]) + 4 * joint_counts[10];
    } else {
      const uint32_t sparse_zmiss_ct = joint_counts[0] + joint_counts[1] + joint_counts[2] + joint_counts[12] + joint_counts[13] + joint_counts[14];
      nmaj_ct1 -= sparse_one_ct + 2 * sparse_zmiss_ct;
      ssq1 -= 3 * sparse_one_ct + 4 * sparse_zmiss_ct;
      dotprod -= joint_counts[5] + 2 * (joint_counts[1] + joint_counts[6] + joint_counts[13]) + 4 * (joint_counts[2] + joint_counts[14]);
    }
  }
  *nmaj_ct0_ptr = nmaj_ct0;
  *nmaj_ct1_ptr = nmaj_ct1;
  *ssq0_ptr = ssq0;
  *ssq1_ptr = ssq1;
  *dotprod_ptr = dotprod;
  return valid_obs_ct;
}

uint32_t ComputeR2NondosageUnphased2SparseStats(const R2NondosageVariant* ndp0, const R2NondosageVariant* ndp1, uint32_t* nmaj_ct0_ptr, uint32_t* nmaj_ct1_ptr, uint32_t* ssq0_ptr, uint32_t* ssq1_ptr, uint32_t* dotprod_ptr) {
  const R2NondosageVariant* longvp;
  const R2NondosageVariant* shortvp;
  uint32_t* nmaj_ctlong_ptr;
  uint32_t* nmaj_ctshort_ptr;
  uint32_t* ssqlong_ptr;
  uint32_t* ssqshort_ptr;
  if (ndp0->p.s.difflist_len <= ndp1->p.s.difflist_len) {
    longvp = ndp1;
    shortvp = ndp0;
    nmaj_ctlong_ptr = nmaj_ct1_ptr;
    nmaj_ctshort_ptr = nmaj_ct0_ptr;
    ssqlong_ptr = ssq1_ptr;
    ssqshort_ptr = ssq0_ptr;
  } else {
    longvp = ndp0;
    shortvp = ndp1;
    nmaj_ctlong_ptr = nmaj_ct0_ptr;
    nmaj_ctshort_ptr = nmaj_ct1_ptr;
    ssqlong_ptr = ssq0_ptr;
    ssqshort_ptr = ssq1_ptr;
  }
  const uint32_t difflist_common_geno_short = shortvp->p.s.difflist_common_geno;
  uint32_t nmaj_ctlong = 0;
  uint32_t nmaj_ctshort = 0;
  uint32_t ssqlong = 0;
  uint32_t ssqshort = 0;
  uint32_t dotprod = 0;
  uint32_t valid_obs_ct = 0;
  if (difflist_common_geno_short != 3) {
    nmaj_ctlong = longvp->nmaj_ct;
    ssqlong = longvp->ssq;
    dotprod = difflist_common_geno_short * nmaj_ctlong;
    valid_obs_ct = longvp->nm_ct;
    nmaj_ctshort = difflist_common_geno_short * valid_obs_ct;
    ssqshort = difflist_common_geno_short * nmaj_ctshort;
  }
  const uint32_t difflist_len_short = shortvp->p.s.difflist_len;
  if (difflist_len_short) {
    // tested genovec in place of {difflist_include,
    // difflist_include_cumulative_popcounts}; that benchmarked worse
    const uintptr_t* difflist_include_long = longvp->p.s.difflist_include;
    const uint32_t* difflist_include_long_cumulative_popcounts = longvp->p.s.difflist_include_cumulative_popcounts;

    const uint32_t* difflist_sample_ids_short = shortvp->p.s.difflist_sample_ids;
    const uintptr_t* raregeno_long = longvp->p.s.raregeno;
    const uintptr_t* raregeno_short = shortvp->p.s.raregeno;
    const uint32_t difflist_common_geno_long = longvp->p.s.difflist_common_geno;
    uint32_t joint_counts[16]; // low bits = long, high bits = short
    ZeroU32Arr(16, joint_counts);
    const uint32_t word_ct_m1 = (difflist_len_short - 1) / kBitsPerWordD2;
    uint32_t loop_len = kBitsPerWordD2;
    for (uint32_t widx = 0; ; ++widx) {
      if (widx >= word_ct_m1) {
        if (widx > word_ct_m1) {
          break;
        }
        loop_len = ModNz(difflist_len_short, kBitsPerWordD2);
      }
      const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids_short[widx * kBitsPerWordD2]);
      uintptr_t raregeno_word = raregeno_short[widx];
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uint32_t sample_idx = cur_difflist_sample_ids[uii];
        const uintptr_t cur_geno_short = raregeno_word & 3;
        uintptr_t cur_geno_long = difflist_common_geno_long;
        if (IsSet(difflist_include_long, sample_idx)) {
          const uint32_t difflist_idx_long = RawToSubsettedPos(difflist_include_long, difflist_include_long_cumulative_popcounts, sample_idx);
          cur_geno_long = GetNyparrEntry(raregeno_long, difflist_idx_long);
        }
        joint_counts[cur_geno_long + 4 * cur_geno_short] += 1;
        raregeno_word = raregeno_word >> 2;
      }
    }
    if (difflist_common_geno_short != 3) {
      nmaj_ctlong -= joint_counts[13] + 2 * joint_counts[14];
      ssqlong -= joint_counts[13] + 4 * joint_counts[14];
      const uint32_t short_missing_long_nm_ct = joint_counts[12] + joint_counts[13] + joint_counts[14];
      valid_obs_ct -= short_missing_long_nm_ct;
    } else {
      const uint32_t long_one_ct = joint_counts[1] + joint_counts[5] + joint_counts[9];
      const uint32_t long_two_ct = joint_counts[2] + joint_counts[6] + joint_counts[10];
      nmaj_ctlong = long_one_ct + 2 * long_two_ct;
      ssqlong = long_one_ct + 4 * long_two_ct;
      const uint32_t long_missing_short_nm_ct = joint_counts[3] + joint_counts[7] + joint_counts[11];
      valid_obs_ct = difflist_len_short - long_missing_short_nm_ct;
    }
    const uint32_t short_one_ct = joint_counts[4] + joint_counts[5] + joint_counts[6];
    if (difflist_common_geno_short != 2) {
      const uint32_t short_two_ct = joint_counts[8] + joint_counts[9] + joint_counts[10];
      nmaj_ctshort = short_one_ct + 2 * short_two_ct;
      ssqshort = short_one_ct + 4 * short_two_ct;
      dotprod = joint_counts[5] + 2 * (joint_counts[6] + joint_counts[9]) + 4 * joint_counts[10];
    } else {
      const uint32_t short_zmiss_ct = joint_counts[0] + joint_counts[1] + joint_counts[2] + joint_counts[12] + joint_counts[13] + joint_counts[14];
      nmaj_ctshort -= short_one_ct + 2 * short_zmiss_ct;
      ssqshort -= 3 * short_one_ct + 4 * short_zmiss_ct;
      dotprod -= joint_counts[5] + 2 * (joint_counts[1] + joint_counts[6] + joint_counts[13]) + 4 * (joint_counts[2] + joint_counts[14]);
    }
  }
  *nmaj_ctlong_ptr = nmaj_ctlong;
  *nmaj_ctshort_ptr = nmaj_ctshort;
  *ssqlong_ptr = ssqlong;
  *ssqshort_ptr = ssqshort;
  *dotprod_ptr = dotprod;
  return valid_obs_ct;
}

// Main return value is valid_obs_ct.  On valid_obs_ct=0, other return values
// may not be filled.  (Same is true for the next three functions.)
uint32_t ComputeR2NondosageUnphasedStats(const R2NondosageVariant* ndp0, const R2NondosageVariant* ndp1, uint32_t sample_ct, uint32_t* nmaj_ct0_ptr, uint32_t* nmaj_ct1_ptr, uint32_t* ssq0_ptr, uint32_t* ssq1_ptr, uint32_t* dotprod_ptr) {
  if (ndp0->is_sparse) {
    if (ndp1->is_sparse) {
      return ComputeR2NondosageUnphased2SparseStats(ndp0, ndp1, nmaj_ct0_ptr, nmaj_ct1_ptr, ssq0_ptr, ssq1_ptr, dotprod_ptr);
    } else {
      return ComputeR2NondosageUnphased1SparseStats(ndp1, ndp0, nmaj_ct1_ptr, nmaj_ct0_ptr, ssq1_ptr, ssq0_ptr, dotprod_ptr);
    }
  }
  if (ndp1->is_sparse) {
    return ComputeR2NondosageUnphased1SparseStats(ndp0, ndp1, nmaj_ct0_ptr, nmaj_ct1_ptr, ssq0_ptr, ssq1_ptr, dotprod_ptr);
  }
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uintptr_t* nm_bitvec0 = ndp0->p.d.nm_bitvec;
  const uintptr_t* nm_bitvec1 = ndp1->p.d.nm_bitvec;
  const uint32_t nm_ct0 = ndp0->nm_ct;
  const uint32_t nm_ct1 = ndp1->nm_ct;
  uint32_t valid_obs_ct;
  if ((nm_ct0 != sample_ct) && (nm_ct1 != sample_ct)) {
    valid_obs_ct = PopcountWordsIntersect(nm_bitvec0, nm_bitvec1, sample_ctl);
    if (!valid_obs_ct) {
      return 0;
    }
  } else {
    valid_obs_ct = MINV(nm_ct0, nm_ct1);
  }
  const uintptr_t* one_bitvec0 = ndp0->p.d.one_bitvec;
  const uintptr_t* two_bitvec0 = ndp0->p.d.two_bitvec;
  if (nm_ct0 == valid_obs_ct) {
    *nmaj_ct0_ptr = ndp0->nmaj_ct;
    *ssq0_ptr = ndp0->ssq;
  } else {
    const uint32_t nmaj_ct0 = GenoBitvecSumSubset(nm_bitvec1, one_bitvec0, two_bitvec0, sample_ctl);
    *nmaj_ct0_ptr = nmaj_ct0;
    // 0, 1, 4 instead of 0, 1, 2
    *ssq0_ptr = nmaj_ct0 + 2 * PopcountWordsIntersect(nm_bitvec1, two_bitvec0, sample_ctl);
  }
  const uintptr_t* one_bitvec1 = ndp1->p.d.one_bitvec;
  const uintptr_t* two_bitvec1 = ndp1->p.d.two_bitvec;
  if (nm_ct1 == valid_obs_ct) {
    *nmaj_ct1_ptr = ndp1->nmaj_ct;
    *ssq1_ptr = ndp1->ssq;
  } else {
    const uint32_t nmaj_ct1 = GenoBitvecSumSubset(nm_bitvec0, one_bitvec1, two_bitvec1, sample_ctl);
    *nmaj_ct1_ptr = nmaj_ct1;
    *ssq1_ptr = nmaj_ct1 + 2 * PopcountWordsIntersect(nm_bitvec0, two_bitvec1, sample_ctl);
  }
  *dotprod_ptr = GenoBitvecUnphasedDotprod(one_bitvec0, two_bitvec0, one_bitvec1, two_bitvec1, sample_ctl);
  return valid_obs_ct;
}

uint32_t ComputeR2NondosagePhasedStats(const R2NondosageVariant* ndp0, const R2NondosageVariant* ndp1, uint32_t sample_ct, R2PhaseType phase_type, double* nmajsums_d, double* known_dotprod_d_ptr, double* unknown_hethet_d_ptr) {
  // See HardcallPhasedR2Stats().  Probable todo: make function names more
  // systematic.
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uintptr_t* nm_bitvec0 = ndp0->p.d.nm_bitvec;
  const uintptr_t* nm_bitvec1 = ndp1->p.d.nm_bitvec;
  const uint32_t nm_ct0 = ndp0->nm_ct;
  const uint32_t nm_ct1 = ndp1->nm_ct;
  uint32_t valid_obs_ct;
  if ((nm_ct0 != sample_ct) && (nm_ct1 != sample_ct)) {
    valid_obs_ct = PopcountWordsIntersect(nm_bitvec0, nm_bitvec1, sample_ctl);
    if (!valid_obs_ct) {
      return 0;
    }
  } else {
    valid_obs_ct = MINV(nm_ct0, nm_ct1);
  }
  const uintptr_t* one_bitvec0 = ndp0->p.d.one_bitvec;
  const uintptr_t* two_bitvec0 = ndp0->p.d.two_bitvec;
  uint32_t nmaj_ct0 = ndp0->nmaj_ct;
  if (nm_ct0 != valid_obs_ct) {
    nmaj_ct0 = GenoBitvecSumSubset(nm_bitvec1, one_bitvec0, two_bitvec0, sample_ctl);
  }
  const uintptr_t* one_bitvec1 = ndp1->p.d.one_bitvec;
  const uintptr_t* two_bitvec1 = ndp1->p.d.two_bitvec;
  uint32_t nmaj_ct1 = ndp1->nmaj_ct;
  if (nm_ct1 != valid_obs_ct) {
    nmaj_ct1 = GenoBitvecSumSubset(nm_bitvec0, one_bitvec1, two_bitvec1, sample_ctl);
  }
  uint32_t known_dotprod;
  uint32_t unknown_hethet_ct;
  GenoBitvecPhasedDotprod(one_bitvec0, two_bitvec0, one_bitvec1, two_bitvec1, sample_ctl, &known_dotprod, &unknown_hethet_ct);
  if ((phase_type == kR2PhaseTypePresent) && (unknown_hethet_ct != 0)) {
    // don't bother with no-phase-here optimization for now
    HardcallPhasedR2Refine(ndp0->p.d.phasepresent, ndp0->p.d.phaseinfo, ndp1->p.d.phasepresent, ndp1->p.d.phaseinfo, sample_ctl, &known_dotprod, &unknown_hethet_ct);
  }
  nmajsums_d[0] = u31tod(nmaj_ct0);
  nmajsums_d[1] = u31tod(nmaj_ct1);
  // bugfix (26 Oct 2023): unknown_hethet treatment shouldn't actually change
  // in haploid case
  *known_dotprod_d_ptr = S_CAST(double, known_dotprod);
  *unknown_hethet_d_ptr = u31tod(unknown_hethet_ct);
  return valid_obs_ct;
}

uint32_t ComputeR2DosageUnphasedStats(const R2DosageVariant* dp0, const R2DosageVariant* dp1, uint32_t sample_ct, uint64_t* nmaj_dosages, uint64_t* dosageprod_ptr, uint64_t* ssq0_ptr, uint64_t* ssq1_ptr) {
  const Dosage* dosage_vec0 = dp0->dosage_vec;
  const Dosage* dosage_vec1 = dp1->dosage_vec;
  const uintptr_t* nm_bitvec0 = dp0->nm_bitvec;
  const uintptr_t* nm_bitvec1 = dp1->nm_bitvec;
  const uint32_t nm_ct0 = dp0->nm_ct;
  const uint32_t nm_ct1 = dp1->nm_ct;
  nmaj_dosages[0] = dp0->nmaj_dosage;
  nmaj_dosages[1] = dp1->nmaj_dosage;
  const uint32_t valid_obs_ct = DosageR2Freqs(dosage_vec0, nm_bitvec0, dosage_vec1, nm_bitvec1, sample_ct, nm_ct0, nm_ct1, nmaj_dosages);
  if (!valid_obs_ct) {
    return 0;
  }
  const uint32_t sample_dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  *dosageprod_ptr = DosageUnsignedDotprod(dosage_vec0, dosage_vec1, sample_dosagev_ct);
  if (nm_ct0 == valid_obs_ct) {
    *ssq0_ptr = dp0->nmaj_dosage_ssq;
  } else {
    // bugfix (24 Oct 2023): this needs to mask out opposite missing values
    *ssq0_ptr = DosageUnsignedDotprodSubset(dosage_vec1, dosage_vec0, dosage_vec0, sample_dosagev_ct);
  }
  if (nm_ct1 == valid_obs_ct) {
    *ssq1_ptr = dp1->nmaj_dosage_ssq;
  } else {
    *ssq1_ptr = DosageUnsignedDotprodSubset(dosage_vec0, dosage_vec1, dosage_vec1, sample_dosagev_ct);
  }
  return valid_obs_ct;
}

uint32_t ComputeR2DosagePhasedStats(const R2DosageVariant* dp0, const R2DosageVariant* dp1, uint32_t sample_ct, R2PhaseType phase_type, double* nmajsums_d, double* known_dotprod_d_ptr, double* unknown_hethet_d_ptr) {
  const Dosage* dosage_vec0 = dp0->dosage_vec;
  const Dosage* dosage_vec1 = dp1->dosage_vec;
  const uintptr_t* nm_bitvec0 = dp0->nm_bitvec;
  const uintptr_t* nm_bitvec1 = dp1->nm_bitvec;
  const uint32_t nm_ct0 = dp0->nm_ct;
  const uint32_t nm_ct1 = dp1->nm_ct;
  uint64_t nmaj_dosages[2];
  nmaj_dosages[0] = dp0->nmaj_dosage;
  nmaj_dosages[1] = dp1->nmaj_dosage;
  const uint32_t valid_obs_ct = DosageR2Freqs(dosage_vec0, nm_bitvec0, dosage_vec1, nm_bitvec1, sample_ct, nm_ct0, nm_ct1, nmaj_dosages);
  if (!valid_obs_ct) {
    return 0;
  }
  const uint32_t sample_dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  uint64_t known_dotprod_dosage;
  uint64_t uhethet_dosage;
  if (phase_type != kR2PhaseTypePresent) {
    const Dosage* dosage_het0 = dp0->dosage_het;
    const Dosage* dosage_het1 = dp1->dosage_het;
    DosageUnphasedDotprodComponents(dosage_vec0, dosage_vec1, dosage_het0, dosage_het1, sample_dosagev_ct, &known_dotprod_dosage, &uhethet_dosage);
  } else {
    const SDosage* dphase_delta0 = dp0->dense_dphase_delta;
    const SDosage* dphase_delta1 = dp1->dense_dphase_delta;
    DosagePhasedDotprodComponents(dosage_vec0, dosage_vec1, dphase_delta0, dphase_delta1, sample_dosagev_ct, &known_dotprod_dosage, &uhethet_dosage);
  }
  nmajsums_d[0] = u63tod(nmaj_dosages[0]) * kRecipDosageMid;
  nmajsums_d[1] = u63tod(nmaj_dosages[1]) * kRecipDosageMid;
  *known_dotprod_d_ptr = u63tod(known_dotprod_dosage) * kRecipDosageMax;
  *unknown_hethet_d_ptr = u63tod(uhethet_dosage) * kRecipDosageMax;
  return valid_obs_ct;
}

// d_ptr and/or dprime_ptr can be nullptr.  Neither are filled in unphased
// case.
double ComputeR2(const R2Variant* r2vp0, const R2Variant* r2vp1, uint32_t sample_ct, R2PhaseType phase_type, uint32_t load_dosage, double* d_ptr, double* dprime_ptr, uint32_t* is_neg_ptr) {
  double nmajsums_d[2];
  double known_dotprod_d;
  double unknown_hethet_d;
  uint32_t valid_obs_ct;
  if (!load_dosage) {
    const R2NondosageVariant* ndp0 = &(r2vp0->nd);
    const R2NondosageVariant* ndp1 = &(r2vp1->nd);
    if (phase_type == kR2PhaseTypeUnphased) {
      uint32_t nmaj_ct0;
      uint32_t nmaj_ct1;
      uint32_t ssq0;
      uint32_t ssq1;
      uint32_t dotprod;
      valid_obs_ct = ComputeR2NondosageUnphasedStats(ndp0, ndp1, sample_ct, &nmaj_ct0, &nmaj_ct1, &ssq0, &ssq1, &dotprod);
      if (!valid_obs_ct) {
        return -DBL_MAX;
      }
      // Previously implemented in e.g. IndepPairwiseThread.
      const int64_t variance0_i64 = ssq0 * S_CAST(int64_t, valid_obs_ct) - S_CAST(int64_t, nmaj_ct0) * nmaj_ct0;
      const int64_t variance1_i64 = ssq1 * S_CAST(int64_t, valid_obs_ct) - S_CAST(int64_t, nmaj_ct1) * nmaj_ct1;
      const double variance_prod = S_CAST(double, variance0_i64) * S_CAST(double, variance1_i64);
      if (variance_prod == 0.0) {
        return -DBL_MAX;
      }
      const double cov01 = S_CAST(double, dotprod * S_CAST(int64_t, valid_obs_ct) - S_CAST(int64_t, nmaj_ct0) * nmaj_ct1);
      *is_neg_ptr = (cov01 < 0.0);
      return cov01 * cov01 / variance_prod;
    }
    valid_obs_ct = ComputeR2NondosagePhasedStats(ndp0, ndp1, sample_ct, phase_type, nmajsums_d, &known_dotprod_d, &unknown_hethet_d);
  } else {
    const R2DosageVariant* dp0 = &(r2vp0->d);
    const R2DosageVariant* dp1 = &(r2vp1->d);
    if (phase_type == kR2PhaseTypeUnphased) {
      uint64_t nmaj_dosages[2];
      uint64_t dosageprod;
      uint64_t ssq0;
      uint64_t ssq1;
      valid_obs_ct = ComputeR2DosageUnphasedStats(dp0, dp1, sample_ct, nmaj_dosages, &dosageprod, &ssq0, &ssq1);
      if (!valid_obs_ct) {
        return -DBL_MAX;
      }
      const double variance0 = u127prod_diff_d(ssq0, valid_obs_ct, nmaj_dosages[0], nmaj_dosages[0]);
      const double variance1 = u127prod_diff_d(ssq1, valid_obs_ct, nmaj_dosages[1], nmaj_dosages[1]);
      const double variance_prod = variance0 * variance1;
      if (variance_prod == 0.0) {
        return -DBL_MAX;
      }
      const double cov01 = i127prod_diff_d(dosageprod, valid_obs_ct, nmaj_dosages[0], nmaj_dosages[1]);
      *is_neg_ptr = (cov01 < 0.0);
      return cov01 * cov01 / variance_prod;
    }
    valid_obs_ct = ComputeR2DosagePhasedStats(dp0, dp1, sample_ct, phase_type, nmajsums_d, &known_dotprod_d, &unknown_hethet_d);
  }
  if (!valid_obs_ct) {
    return -DBL_MAX;
  }
  const double twice_tot_recip = 0.5 / u31tod(valid_obs_ct);
  if ((d_ptr == nullptr) && (dprime_ptr == nullptr)) {
    double r2;
    const LDErr ld_err = PhasedLD(nmajsums_d, known_dotprod_d, unknown_hethet_d, twice_tot_recip, 0, nullptr, &r2, is_neg_ptr);
    return (ld_err == kLDErrNone)? r2 : -DBL_MAX;
  }
  double results[3];
  const LDErr ld_err = PhasedLD(nmajsums_d, known_dotprod_d, unknown_hethet_d, twice_tot_recip, 1, nullptr, results, is_neg_ptr);
  if (ld_err != kLDErrNone) {
    if (d_ptr) {
      *d_ptr = 0.0 / 0.0;
    }
    if (dprime_ptr) {
      *dprime_ptr = 0.0 / 0.0;
    }
    return -DBL_MAX;
  }
  if (d_ptr) {
    *d_ptr = results[1];
  }
  if (dprime_ptr) {
    *dprime_ptr = results[2];
  }
  return results[0];
}

uint32_t ComputeR2NondosageUnphased1SparseSubsetStats(const R2NondosageVariant* densevp0, const R2NondosageVariant* sparsevp1, const uintptr_t* sample_include, uint32_t subsetted_nm_ct0, uint32_t* nmaj_ct0_ptr, uint32_t* nmaj_ct1_ptr, uint32_t* ssq0_ptr, uint32_t* ssq1_ptr, uint32_t* dotprod_ptr) {
  const uint32_t difflist_common_geno = sparsevp1->p.s.difflist_common_geno;
  const uint32_t difflist_len = sparsevp1->p.s.difflist_len;
  const uint32_t* difflist_sample_ids = sparsevp1->p.s.difflist_sample_ids;
  uint32_t nmaj_ct0 = 0;
  uint32_t nmaj_ct1 = 0;
  uint32_t ssq0 = 0;
  uint32_t ssq1 = 0;
  uint32_t dotprod = 0;
  uint32_t valid_obs_ct = 0;
  if (difflist_common_geno != 3) {
    nmaj_ct0 = *nmaj_ct0_ptr;
    ssq0 = *ssq0_ptr;
    dotprod = difflist_common_geno * nmaj_ct0;
    valid_obs_ct = subsetted_nm_ct0;
    nmaj_ct1 = difflist_common_geno * valid_obs_ct;
    ssq1 = difflist_common_geno * nmaj_ct1;
  }
  if (difflist_len) {
    const uintptr_t* nm_bitvec0 = densevp0->p.d.nm_bitvec;
    const uintptr_t* one_bitvec0 = densevp0->p.d.one_bitvec;
    const uintptr_t* two_bitvec0 = densevp0->p.d.two_bitvec;
    const uintptr_t* raregeno = sparsevp1->p.s.raregeno;
    const uint32_t word_ct_m1 = (difflist_len - 1) / kBitsPerWordD2;
    uint32_t joint_counts[16]; // low bits = dense geno, high bits = sparse
    ZeroU32Arr(16, joint_counts);
    uint32_t loop_len = kBitsPerWordD2;
    for (uint32_t widx = 0; ; ++widx) {
      if (widx >= word_ct_m1) {
        if (widx > word_ct_m1) {
          break;
        }
        loop_len = ModNz(difflist_len, kBitsPerWordD2);
      }
      const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
      uintptr_t raregeno_word = raregeno[widx];
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uint32_t sample_idx = cur_difflist_sample_ids[uii];
        const uint32_t sample_widx = sample_idx / kBitsPerWord;
        const uint32_t sample_idx_lowbits = sample_idx % kBitsPerWord;
        if ((sample_include[sample_widx] >> sample_idx_lowbits) & 1) {
          const uintptr_t cur_sparse_geno = raregeno_word & 3;
          const uintptr_t nm_bit = (nm_bitvec0[sample_widx] >> sample_idx_lowbits) & 1;
          const uintptr_t one_bit = (one_bitvec0[sample_widx] >> sample_idx_lowbits) & 1;
          // "& 3" takes care of this mask
          const uintptr_t two_bit_unmasked = two_bitvec0[sample_widx] >> sample_idx_lowbits;
          const uintptr_t cur_dense_geno = (nm_bit + one_bit + two_bit_unmasked * 2 - 1) & 3;
          joint_counts[cur_dense_geno + cur_sparse_geno * 4] += 1;
        }
        raregeno_word = raregeno_word >> 2;
      }
    }
    if (difflist_common_geno != 3) {
      nmaj_ct0 -= joint_counts[13] + 2 * joint_counts[14];
      ssq0 -= joint_counts[13] + 4 * joint_counts[14];
      const uint32_t sparse_missing_dense_nm_ct = joint_counts[12] + joint_counts[13] + joint_counts[14];
      valid_obs_ct -= sparse_missing_dense_nm_ct;
    } else {
      const uint32_t dense_one_ct = joint_counts[1] + joint_counts[5] + joint_counts[9];
      const uint32_t dense_two_ct = joint_counts[2] + joint_counts[6] + joint_counts[10];
      nmaj_ct0 = dense_one_ct + 2 * dense_two_ct;
      ssq0 = dense_one_ct + 4 * dense_two_ct;
      const uint32_t dense_zero_ct = joint_counts[0] + joint_counts[4] + joint_counts[8];
      valid_obs_ct = dense_zero_ct + dense_one_ct + dense_two_ct;
    }
    const uint32_t sparse_one_ct = joint_counts[4] + joint_counts[5] + joint_counts[6];
    if (difflist_common_geno != 2) {
      const uint32_t sparse_two_ct = joint_counts[8] + joint_counts[9] + joint_counts[10];
      nmaj_ct1 = sparse_one_ct + 2 * sparse_two_ct;
      ssq1 = sparse_one_ct + 4 * sparse_two_ct;
      dotprod = joint_counts[5] + 2 * (joint_counts[6] + joint_counts[9]) + 4 * joint_counts[10];
    } else {
      const uint32_t sparse_zmiss_ct = joint_counts[0] + joint_counts[1] + joint_counts[2] + joint_counts[12] + joint_counts[13] + joint_counts[14];
      nmaj_ct1 -= sparse_one_ct + 2 * sparse_zmiss_ct;
      ssq1 -= 3 * sparse_one_ct + 4 * sparse_zmiss_ct;
      dotprod -= joint_counts[5] + 2 * (joint_counts[1] + joint_counts[6] + joint_counts[13]) + 4 * (joint_counts[2] + joint_counts[14]);
    }
  }
  *nmaj_ct0_ptr = nmaj_ct0;
  *nmaj_ct1_ptr = nmaj_ct1;
  *ssq0_ptr = ssq0;
  *ssq1_ptr = ssq1;
  *dotprod_ptr = dotprod;
  return valid_obs_ct;
}

uint32_t ComputeR2NondosageUnphased2SparseSubsetStats(const R2NondosageVariant* ndp0, const R2NondosageVariant* ndp1, const uintptr_t* sample_include, uint32_t subsetted_nm_ct0, uint32_t subsetted_nm_ct1, uint32_t* nmaj_ct0_ptr, uint32_t* nmaj_ct1_ptr, uint32_t* ssq0_ptr, uint32_t* ssq1_ptr, uint32_t* dotprod_ptr) {
  const R2NondosageVariant* longvp;
  const R2NondosageVariant* shortvp;
  uint32_t* nmaj_ctlong_ptr;
  uint32_t* nmaj_ctshort_ptr;
  uint32_t* ssqlong_ptr;
  uint32_t* ssqshort_ptr;
  uint32_t long_nm_ct;
  if (ndp0->p.s.difflist_len <= ndp1->p.s.difflist_len) {
    longvp = ndp1;
    shortvp = ndp0;
    nmaj_ctlong_ptr = nmaj_ct1_ptr;
    nmaj_ctshort_ptr = nmaj_ct0_ptr;
    ssqlong_ptr = ssq1_ptr;
    ssqshort_ptr = ssq0_ptr;
    long_nm_ct = subsetted_nm_ct1;
  } else {
    longvp = ndp0;
    shortvp = ndp1;
    nmaj_ctlong_ptr = nmaj_ct0_ptr;
    nmaj_ctshort_ptr = nmaj_ct1_ptr;
    ssqlong_ptr = ssq0_ptr;
    ssqshort_ptr = ssq1_ptr;
    long_nm_ct = subsetted_nm_ct0;
  }
  const uint32_t difflist_common_geno_short = shortvp->p.s.difflist_common_geno;
  uint32_t nmaj_ctlong = 0;
  uint32_t nmaj_ctshort = 0;
  uint32_t ssqlong = 0;
  uint32_t ssqshort = 0;
  uint32_t dotprod = 0;
  uint32_t valid_obs_ct = 0;
  if (difflist_common_geno_short != 3) {
    nmaj_ctlong = *nmaj_ctlong_ptr;
    ssqlong = *ssqlong_ptr;
    dotprod = difflist_common_geno_short * nmaj_ctlong;
    valid_obs_ct = long_nm_ct;
    nmaj_ctshort = difflist_common_geno_short * valid_obs_ct;
    ssqshort = difflist_common_geno_short * nmaj_ctshort;
  }
  const uint32_t difflist_len_short = shortvp->p.s.difflist_len;
  if (difflist_len_short) {
    const uintptr_t* difflist_include_long = longvp->p.s.difflist_include;
    const uint32_t* difflist_include_long_cumulative_popcounts = longvp->p.s.difflist_include_cumulative_popcounts;

    const uint32_t* difflist_sample_ids_short = shortvp->p.s.difflist_sample_ids;
    const uintptr_t* raregeno_long = longvp->p.s.raregeno;
    const uintptr_t* raregeno_short = shortvp->p.s.raregeno;
    const uint32_t difflist_common_geno_long = longvp->p.s.difflist_common_geno;
    uint32_t joint_counts[16]; // low bits = long, high bits = short
    ZeroU32Arr(16, joint_counts);
    const uint32_t word_ct_m1 = (difflist_len_short - 1) / kBitsPerWordD2;
    uint32_t loop_len = kBitsPerWordD2;
    for (uint32_t widx = 0; ; ++widx) {
      if (widx >= word_ct_m1) {
        if (widx > word_ct_m1) {
          break;
        }
        loop_len = ModNz(difflist_len_short, kBitsPerWordD2);
      }
      const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids_short[widx * kBitsPerWordD2]);
      uintptr_t raregeno_word = raregeno_short[widx];
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uint32_t sample_idx = cur_difflist_sample_ids[uii];
        if (IsSet(sample_include, sample_idx)) {
          const uintptr_t cur_geno_short = raregeno_word & 3;
          uintptr_t cur_geno_long = difflist_common_geno_long;
          if (IsSet(difflist_include_long, sample_idx)) {
            const uint32_t difflist_idx_long = RawToSubsettedPos(difflist_include_long, difflist_include_long_cumulative_popcounts, sample_idx);
            cur_geno_long = GetNyparrEntry(raregeno_long, difflist_idx_long);
          }
          joint_counts[cur_geno_long + 4 * cur_geno_short] += 1;
        }
        raregeno_word = raregeno_word >> 2;
      }
    }
    if (difflist_common_geno_short != 3) {
      nmaj_ctlong -= joint_counts[13] + 2 * joint_counts[14];
      ssqlong -= joint_counts[13] + 4 * joint_counts[14];
      const uint32_t short_missing_long_nm_ct = joint_counts[12] + joint_counts[13] + joint_counts[14];
      valid_obs_ct -= short_missing_long_nm_ct;
    } else {
      const uint32_t long_one_ct = joint_counts[1] + joint_counts[5] + joint_counts[9];
      const uint32_t long_two_ct = joint_counts[2] + joint_counts[6] + joint_counts[10];
      nmaj_ctlong = long_one_ct + 2 * long_two_ct;
      ssqlong = long_one_ct + 4 * long_two_ct;
      const uint32_t long_zero_ct = joint_counts[0] + joint_counts[4] + joint_counts[8];
      valid_obs_ct = long_zero_ct + long_one_ct + long_two_ct;
    }
    const uint32_t short_one_ct = joint_counts[4] + joint_counts[5] + joint_counts[6];
    if (difflist_common_geno_short != 2) {
      const uint32_t short_two_ct = joint_counts[8] + joint_counts[9] + joint_counts[10];
      nmaj_ctshort = short_one_ct + 2 * short_two_ct;
      ssqshort = short_one_ct + 4 * short_two_ct;
      dotprod = joint_counts[5] + 2 * (joint_counts[6] + joint_counts[9]) + 4 * joint_counts[10];
    } else {
      const uint32_t short_zmiss_ct = joint_counts[0] + joint_counts[1] + joint_counts[2] + joint_counts[12] + joint_counts[13] + joint_counts[14];
      nmaj_ctshort -= short_one_ct + 2 * short_zmiss_ct;
      ssqshort -= 3 * short_one_ct + 4 * short_zmiss_ct;
      dotprod -= joint_counts[5] + 2 * (joint_counts[1] + joint_counts[6] + joint_counts[13]) + 4 * (joint_counts[2] + joint_counts[14]);
    }
  }
  *nmaj_ctlong_ptr = nmaj_ctlong;
  *nmaj_ctshort_ptr = nmaj_ctshort;
  *ssqlong_ptr = ssqlong;
  *ssqshort_ptr = ssqshort;
  *dotprod_ptr = dotprod;
  return valid_obs_ct;
}

// nmaj_ct0, nmaj_ct1, ssq0, and ssq1 assumed to be initialized to precomputed
// subsetted values.
uint32_t ComputeR2NondosageUnphasedSubsetStats(const R2NondosageVariant* ndp0, const R2NondosageVariant* ndp1, const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t subsetted_nm_ct0, uint32_t subsetted_nm_ct1, uint32_t* nmaj_ct0_ptr, uint32_t* nmaj_ct1_ptr, uint32_t* ssq0_ptr, uint32_t* ssq1_ptr, uint32_t* dotprod_ptr, uintptr_t* cur_nm_buf) {
  if (ndp0->is_sparse) {
    if (ndp1->is_sparse) {
      return ComputeR2NondosageUnphased2SparseSubsetStats(ndp0, ndp1, sample_include, subsetted_nm_ct0, subsetted_nm_ct1, nmaj_ct0_ptr, nmaj_ct1_ptr, ssq0_ptr, ssq1_ptr, dotprod_ptr);
    } else {
      return ComputeR2NondosageUnphased1SparseSubsetStats(ndp1, ndp0, sample_include, subsetted_nm_ct1, nmaj_ct1_ptr, nmaj_ct0_ptr, ssq1_ptr, ssq0_ptr, dotprod_ptr);
    }
  }
  if (ndp1->is_sparse) {
    return ComputeR2NondosageUnphased1SparseSubsetStats(ndp0, ndp1, sample_include, subsetted_nm_ct0, nmaj_ct0_ptr, nmaj_ct1_ptr, ssq0_ptr, ssq1_ptr, dotprod_ptr);
  }
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uintptr_t* nm_bitvec0 = ndp0->p.d.nm_bitvec;
  const uintptr_t* nm_bitvec1 = ndp1->p.d.nm_bitvec;
  const uintptr_t* cur_nm;
  uint32_t valid_obs_ct;
  if (subsetted_nm_ct0 == sample_ct) {
    valid_obs_ct = subsetted_nm_ct1;
    if (subsetted_nm_ct1 == sample_ct) {
      cur_nm = sample_include;
    } else {
      BitvecAndCopy(nm_bitvec1, sample_include, raw_sample_ctl, cur_nm_buf);
      cur_nm = cur_nm_buf;
    }
  } else {
    valid_obs_ct = subsetted_nm_ct0;
    BitvecAndCopy(nm_bitvec0, sample_include, raw_sample_ctl, cur_nm_buf);
    if (subsetted_nm_ct1 != sample_ct) {
      BitvecAnd(nm_bitvec1, raw_sample_ctl, cur_nm_buf);
      valid_obs_ct = PopcountWords(cur_nm_buf, raw_sample_ctl);
      if (!valid_obs_ct) {
        // bugfix (29 Oct 2023)
        *nmaj_ct0_ptr = 0;
        *nmaj_ct1_ptr = 0;
        *ssq0_ptr = 0;
        *ssq1_ptr = 0;
        *dotprod_ptr = 0;
        return 0;
      }
    }
    cur_nm = cur_nm_buf;
  }
  const uintptr_t* one_bitvec0 = ndp0->p.d.one_bitvec;
  const uintptr_t* two_bitvec0 = ndp0->p.d.two_bitvec;
  if (subsetted_nm_ct0 != valid_obs_ct) {
    const uint32_t nmaj_ct0 = GenoBitvecSumSubset(cur_nm, one_bitvec0, two_bitvec0, raw_sample_ctl);
    *nmaj_ct0_ptr = nmaj_ct0;
    // 0, 1, 4 instead of 0, 1, 2
    *ssq0_ptr = nmaj_ct0 + 2 * PopcountWordsIntersect(cur_nm, two_bitvec0, raw_sample_ctl);
  }

  const uintptr_t* one_bitvec1 = ndp1->p.d.one_bitvec;
  const uintptr_t* two_bitvec1 = ndp1->p.d.two_bitvec;
  if (subsetted_nm_ct1 != valid_obs_ct) {
    const uint32_t nmaj_ct1 = GenoBitvecSumSubset(cur_nm, one_bitvec1, two_bitvec1, raw_sample_ctl);
    *nmaj_ct1_ptr = nmaj_ct1;
    *ssq1_ptr = nmaj_ct1 + 2 * PopcountWordsIntersect(cur_nm, two_bitvec1, raw_sample_ctl);
  }
  *dotprod_ptr = GenoBitvecUnphasedDotprodSubset(cur_nm, one_bitvec0, two_bitvec0, one_bitvec1, two_bitvec1, raw_sample_ctl);
  return valid_obs_ct;
}

uint32_t ComputeR2NondosagePhasedSubsetStats(const R2NondosageVariant* ndp0, const R2NondosageVariant* ndp1, const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t subsetted_nm_ct0, uint32_t subsetted_nm_ct1, uint32_t subsetted_nmaj_ct0, uint32_t subsetted_nmaj_ct1, R2PhaseType phase_type, double* nmajsums_d, double* known_dotprod_d_ptr, double* unknown_hethet_d_ptr, uintptr_t* cur_nm_buf) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uintptr_t* nm_bitvec0 = ndp0->p.d.nm_bitvec;
  const uintptr_t* nm_bitvec1 = ndp1->p.d.nm_bitvec;
  const uintptr_t* cur_nm;
  uint32_t valid_obs_ct;
  if (subsetted_nm_ct0 == sample_ct) {
    valid_obs_ct = subsetted_nm_ct1;
    if (subsetted_nm_ct1 == sample_ct) {
      cur_nm = sample_include;
    } else {
      BitvecAndCopy(nm_bitvec1, sample_include, raw_sample_ctl, cur_nm_buf);
      cur_nm = cur_nm_buf;
    }
  } else {
    valid_obs_ct = subsetted_nm_ct0;
    BitvecAndCopy(nm_bitvec0, sample_include, raw_sample_ctl, cur_nm_buf);
    if (subsetted_nm_ct1 != sample_ct) {
      BitvecAnd(nm_bitvec1, raw_sample_ctl, cur_nm_buf);
      valid_obs_ct = PopcountWords(cur_nm_buf, raw_sample_ctl);
      if (!valid_obs_ct) {
        // bugfix (29 Oct 2023)
        nmajsums_d[0] = 0.0;
        nmajsums_d[1] = 0.0;
        *known_dotprod_d_ptr = 0;
        *unknown_hethet_d_ptr = 0;
        return 0;
      }
    }
    cur_nm = cur_nm_buf;
  }
  const uintptr_t* one_bitvec0 = ndp0->p.d.one_bitvec;
  const uintptr_t* two_bitvec0 = ndp0->p.d.two_bitvec;
  uint32_t nmaj_ct0 = subsetted_nmaj_ct0;
  if (subsetted_nm_ct0 != valid_obs_ct) {
    nmaj_ct0 = GenoBitvecSumSubset(cur_nm, one_bitvec0, two_bitvec0, raw_sample_ctl);
  }

  const uintptr_t* one_bitvec1 = ndp1->p.d.one_bitvec;
  const uintptr_t* two_bitvec1 = ndp1->p.d.two_bitvec;
  uint32_t nmaj_ct1 = subsetted_nmaj_ct1;
  if (subsetted_nm_ct1 != valid_obs_ct) {
    nmaj_ct1 = GenoBitvecSumSubset(cur_nm, one_bitvec1, two_bitvec1, raw_sample_ctl);
  }

  uint32_t known_dotprod;
  uint32_t unknown_hethet_ct;
  GenoBitvecPhasedDotprodSubset(cur_nm, one_bitvec0, two_bitvec0, one_bitvec1, two_bitvec1, raw_sample_ctl, &known_dotprod, &unknown_hethet_ct);
  if ((phase_type == kR2PhaseTypePresent) && (unknown_hethet_ct != 0)) {
    // don't bother with no-phase-here optimization for now
    HardcallPhasedR2RefineSubset(cur_nm, ndp0->p.d.phasepresent, ndp0->p.d.phaseinfo, ndp1->p.d.phasepresent, ndp1->p.d.phaseinfo, raw_sample_ctl, &known_dotprod, &unknown_hethet_ct);
  }
  nmajsums_d[0] = u31tod(nmaj_ct0);
  nmajsums_d[1] = u31tod(nmaj_ct1);
  *known_dotprod_d_ptr = S_CAST(double, known_dotprod);
  *unknown_hethet_d_ptr = u31tod(unknown_hethet_ct);
  return valid_obs_ct;
}

// nmaj_dosages, ssq0, and ssq1 assumed to be initialized to subset values
uint32_t ComputeR2DosageUnphasedSubsetStats(const R2DosageVariant* dp0, const R2DosageVariant* dp1, const uintptr_t* sample_include, const Dosage* dosage_subset_invmask, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t subsetted_nm_ct0, uint32_t subsetted_nm_ct1, uint64_t* nmaj_dosages, uint64_t* ssq0_ptr, uint64_t* ssq1_ptr, uint64_t* dosageprod_ptr, uintptr_t* cur_nm_buf, Dosage* invmask_buf) {
  const Dosage* dosage_vec0 = dp0->dosage_vec;
  const Dosage* dosage_vec1 = dp1->dosage_vec;
  const uintptr_t* nm_bitvec0 = dp0->nm_bitvec;
  const uintptr_t* nm_bitvec1 = dp1->nm_bitvec;
  const uint32_t valid_obs_ct = DosageR2FreqsSubset(dosage_vec0, nm_bitvec0, dosage_vec1, nm_bitvec1, sample_include, raw_sample_ct, sample_ct, subsetted_nm_ct0, subsetted_nm_ct1, nmaj_dosages, cur_nm_buf, invmask_buf);
  if (!valid_obs_ct) {
    // bugfix (29 Oct 2023)
    *dosageprod_ptr = 0.0;
    *ssq0_ptr = 0.0;
    *ssq1_ptr = 0.0;
    return 0;
  }
  const Dosage* subset_invmask = (valid_obs_ct == sample_ct)? dosage_subset_invmask : invmask_buf;
  const uint32_t raw_sample_dosagev_ct = DivUp(raw_sample_ct, kDosagePerVec);
  *dosageprod_ptr = DosageUnsignedDotprodSubset(subset_invmask, dosage_vec0, dosage_vec1, raw_sample_dosagev_ct);
  if (subsetted_nm_ct0 != valid_obs_ct) {
    *ssq0_ptr = DosageUnsignedDotprodSubset(subset_invmask, dosage_vec0, dosage_vec0, raw_sample_dosagev_ct);
  }
  if (subsetted_nm_ct1 != valid_obs_ct) {
    *ssq1_ptr = DosageUnsignedDotprodSubset(subset_invmask, dosage_vec1, dosage_vec1, raw_sample_dosagev_ct);
  }
  return valid_obs_ct;
}

// Caller is now responsible for setting nmajsums_d[k] =
// u63tod(nmaj_dosages[k]) * kRecipDosageMid afterward.
// bugfix (29 Oct 2023): caller still shouldn't be responsible for initializing
// known_dotprod_d and unknown_hethet_d, so we need to zero them out on the
// early-return.
uint32_t ComputeR2DosagePhasedSubsetStats(const R2DosageVariant* dp0, const R2DosageVariant* dp1, const uintptr_t* sample_include, const Dosage* dosage_subset_invmask, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t subsetted_nm_ct0, uint32_t subsetted_nm_ct1, R2PhaseType phase_type, uint64_t* nmaj_dosages, double* known_dotprod_d_ptr, double* unknown_hethet_d_ptr, uintptr_t* cur_nm_buf, Dosage* invmask_buf) {
  const Dosage* dosage_vec0 = dp0->dosage_vec;
  const Dosage* dosage_vec1 = dp1->dosage_vec;
  const uintptr_t* nm_bitvec0 = dp0->nm_bitvec;
  const uintptr_t* nm_bitvec1 = dp1->nm_bitvec;
  const uint32_t valid_obs_ct = DosageR2FreqsSubset(dosage_vec0, nm_bitvec0, dosage_vec1, nm_bitvec1, sample_include, raw_sample_ct, sample_ct, subsetted_nm_ct0, subsetted_nm_ct1, nmaj_dosages, cur_nm_buf, invmask_buf);
  if (!valid_obs_ct) {
    *known_dotprod_d_ptr = 0;
    *unknown_hethet_d_ptr = 0;
    return 0;
  }
  const uint32_t raw_sample_dosagev_ct = DivUp(raw_sample_ct, kDosagePerVec);
  const Dosage* subset_invmask = (valid_obs_ct == sample_ct)? dosage_subset_invmask : invmask_buf;
  uint64_t known_dotprod_dosage;
  uint64_t uhethet_dosage;
  if (phase_type != kR2PhaseTypePresent) {
    const Dosage* dosage_het0 = dp0->dosage_het;
    const Dosage* dosage_het1 = dp1->dosage_het;
    DosageUnphasedDotprodComponentsSubset(subset_invmask, dosage_vec0, dosage_vec1, dosage_het0, dosage_het1, raw_sample_dosagev_ct, &known_dotprod_dosage, &uhethet_dosage);
  } else {
    const SDosage* dphase_delta0 = dp0->dense_dphase_delta;
    const SDosage* dphase_delta1 = dp1->dense_dphase_delta;
    DosagePhasedDotprodComponentsSubset(subset_invmask, dosage_vec0, dosage_vec1, dphase_delta0, dphase_delta1, raw_sample_dosagev_ct, &known_dotprod_dosage, &uhethet_dosage);
  }
  *known_dotprod_d_ptr = u63tod(known_dotprod_dosage) * kRecipDosageMax;
  *unknown_hethet_d_ptr = u63tod(uhethet_dosage) * kRecipDosageMax;
  return valid_obs_ct;
}

// There initially was a separate code path for the non-inter-chr case, where
// the male and nonmale genotypes were pre-separated.  It was useful as a
// somewhat-independent implementation to test ComputeXR2() against.  But it's
// now deleted, despite being faster, since chrX-specific code is a defect
// attractor and does not cover a large fraction of typical computational
// loads.
double ComputeXR2(const R2Variant* r2vp0, const R2Variant* r2vp1, const uintptr_t* founder_male_collapsed, const uintptr_t* founder_nonmale_collapsed, const Dosage* male_dosage_invmask, const Dosage* nonmale_dosage_invmask, uint32_t sample_ct, uint32_t male_ct, R2PhaseType phase_type, uint32_t load_dosage, uint32_t both_x, double* d_ptr, double* dprime_ptr, uint32_t* is_neg_ptr, uintptr_t* cur_nm_buf, Dosage* invmask_buf) {
  const double male_downwt = both_x? 0.5 : (1.0 - 0.5 * kSqrt2);

  double male_nmajsums_d[2];
  double male_known_dotprod_d;
  double male_unknown_hethet_d;
  uint32_t male_obs_ct;

  double nonmale_nmajsums_d[2];
  double nonmale_known_dotprod_d;
  double nonmale_unknown_hethet_d;
  uint32_t nonmale_obs_ct;
  if (!load_dosage) {
    const R2NondosageVariant* ndp0 = &(r2vp0->nd);
    const R2NondosageVariant* ndp1 = &(r2vp1->nd);
    if (phase_type == kR2PhaseTypeUnphased) {
      uint32_t nmaj_ct0;
      uint32_t nmaj_ct1;
      uint32_t ssq0;
      uint32_t ssq1;
      uint32_t dotprod;
      const uint32_t valid_obs_ct = ComputeR2NondosageUnphasedStats(ndp0, ndp1, sample_ct, &nmaj_ct0, &nmaj_ct1, &ssq0, &ssq1, &dotprod);
      if (!valid_obs_ct) {
        return -DBL_MAX;
      }
      /*
      printf("nmaj_ct0: %u\n", nmaj_ct0);
      printf("nmaj_ct1: %u\n", nmaj_ct1);
      printf("ssq0: %u\n", ssq0);
      printf("ssq1: %u\n", ssq1);
      printf("dotprod: %u\n", dotprod);
      printf("valid_obs_ct: %u\n", valid_obs_ct);
      */
      uint32_t male_nmaj_ct0 = ndp0->x_male_nmaj_ct;
      uint32_t male_nmaj_ct1 = ndp1->x_male_nmaj_ct;
      uint32_t male_ssq0 = ndp0->x_male_ssq;
      uint32_t male_ssq1 = ndp1->x_male_ssq;
      uint32_t male_dotprod;
      // printf("male_ssq1 before: %u\n", male_ssq1);
      male_obs_ct = ComputeR2NondosageUnphasedSubsetStats(ndp0, ndp1, founder_male_collapsed, sample_ct, male_ct, ndp0->x_male_nm_ct, ndp1->x_male_nm_ct, &male_nmaj_ct0, &male_nmaj_ct1, &male_ssq0, &male_ssq1, &male_dotprod, cur_nm_buf);
      /*
      printf("male_nmaj_ct0: %u\n", male_nmaj_ct0);
      printf("male_nmaj_ct1: %u\n", male_nmaj_ct1);
      printf("male_ssq0: %u\n", male_ssq0);
      printf("male_ssq1: %u\n", male_ssq1);
      printf("male_dotprod: %u\n", male_dotprod);
      printf("male_obs_ct: %u\n", male_obs_ct);
      */

      const double weighted_obs_ct = u31tod(valid_obs_ct) - male_downwt * u31tod(male_obs_ct);
      const double weighted_nmaj_ct0 = u31tod(nmaj_ct0) - male_downwt * u31tod(male_nmaj_ct0);
      const double weighted_nmaj_ct1 = u31tod(nmaj_ct1) - male_downwt * u31tod(male_nmaj_ct1);
      const double weighted_ssq0 = u63tod(ssq0) - male_downwt * u63tod(male_ssq0);
      const double weighted_ssq1 = u63tod(ssq1) - male_downwt * u63tod(male_ssq1);
      const double weighted_dotprod = u63tod(dotprod) - male_downwt * u63tod(male_dotprod);

      const double variance0 = weighted_ssq0 * weighted_obs_ct - weighted_nmaj_ct0 * weighted_nmaj_ct0;
      const double variance1 = weighted_ssq1 * weighted_obs_ct - weighted_nmaj_ct1 * weighted_nmaj_ct1;
      if ((variance0 <= 0.0) || (variance1 <= 0.0)) {
        return -DBL_MAX;
      }
      const double variance_prod = variance0 * variance1;
      const double cov01 = weighted_dotprod * weighted_obs_ct - weighted_nmaj_ct0 * weighted_nmaj_ct1;
      *is_neg_ptr = (cov01 < 0.0);
      return MINV(1.0, cov01 * cov01 / variance_prod);
    }
    const uint32_t x_male_nm_ct0 = ndp0->x_male_nm_ct;
    const uint32_t x_male_nm_ct1 = ndp1->x_male_nm_ct;
    const uint32_t x_male_nmaj_ct0 = ndp0->x_male_nmaj_ct;
    const uint32_t x_male_nmaj_ct1 = ndp1->x_male_nmaj_ct;
    male_obs_ct = ComputeR2NondosagePhasedSubsetStats(ndp0, ndp1, founder_male_collapsed, sample_ct, male_ct, x_male_nm_ct0, x_male_nm_ct1, x_male_nmaj_ct0, x_male_nmaj_ct1, R2PhaseOmit(phase_type), male_nmajsums_d, &male_known_dotprod_d, &male_unknown_hethet_d, cur_nm_buf);
    nonmale_obs_ct = ComputeR2NondosagePhasedSubsetStats(ndp0, ndp1, founder_nonmale_collapsed, sample_ct, sample_ct - male_ct, ndp0->nm_ct - x_male_nm_ct0, ndp1->nm_ct - x_male_nm_ct1, ndp0->nmaj_ct - x_male_nmaj_ct0, ndp1->nmaj_ct - x_male_nmaj_ct1, phase_type, nonmale_nmajsums_d, &nonmale_known_dotprod_d, &nonmale_unknown_hethet_d, cur_nm_buf);
  } else {
    const R2DosageVariant* dp0 = &(r2vp0->d);
    const R2DosageVariant* dp1 = &(r2vp1->d);
    uint64_t sex_nmaj_dosages[2];
    sex_nmaj_dosages[0] = dp0->x_male_nmaj_dosage;
    sex_nmaj_dosages[1] = dp1->x_male_nmaj_dosage;
    if (phase_type == kR2PhaseTypeUnphased) {
      uint64_t nmaj_dosages[2];
      uint64_t dosageprod;
      uint64_t ssq0;
      uint64_t ssq1;
      const uint32_t valid_obs_ct = ComputeR2DosageUnphasedStats(dp0, dp1, sample_ct, nmaj_dosages, &dosageprod, &ssq0, &ssq1);
      if (!valid_obs_ct) {
        return -DBL_MAX;
      }
      uint64_t male_ssq0 = dp0->x_male_nmaj_dosage_ssq;
      uint64_t male_ssq1 = dp1->x_male_nmaj_dosage_ssq;
      uint64_t male_dosageprod;
      male_obs_ct = ComputeR2DosageUnphasedSubsetStats(dp0, dp1, founder_male_collapsed, male_dosage_invmask, sample_ct, male_ct, dp0->x_male_nm_ct, dp1->x_male_nm_ct, sex_nmaj_dosages, &male_ssq0, &male_ssq1, &male_dosageprod, cur_nm_buf, invmask_buf);

      const double weighted_obs_ct = u31tod(valid_obs_ct) - male_downwt * u31tod(male_obs_ct);
      const double weighted_nmaj0 = u63tod(nmaj_dosages[0]) - male_downwt * u63tod(sex_nmaj_dosages[0]);
      const double weighted_nmaj1 = u63tod(nmaj_dosages[1]) - male_downwt * u63tod(sex_nmaj_dosages[1]);
      const double weighted_dosageprod = u63tod(dosageprod) - male_downwt * u63tod(male_dosageprod);
      const double weighted_ssq0 = u63tod(ssq0) - male_downwt * u63tod(male_ssq0);
      const double weighted_ssq1 = u63tod(ssq1) - male_downwt * u63tod(male_ssq1);

      const double variance0 = weighted_ssq0 * weighted_obs_ct - weighted_nmaj0 * weighted_nmaj0;
      const double variance1 = weighted_ssq1 * weighted_obs_ct - weighted_nmaj1 * weighted_nmaj1;
      if ((variance0 <= 0.0) || (variance1 <= 0.0)) {
        return -DBL_MAX;
      }
      const double variance_prod = variance0 * variance1;
      if (variance_prod == 0.0) {
        return -DBL_MAX;
      }
      const double cov01 = weighted_dosageprod * weighted_obs_ct - weighted_nmaj0 * weighted_nmaj1;
      *is_neg_ptr = (cov01 < 0.0);
      return MINV(1.0, cov01 * cov01 / variance_prod);
    }
    const uint32_t x_male_nm_ct0 = dp0->x_male_nm_ct;
    const uint32_t x_male_nm_ct1 = dp1->x_male_nm_ct;
    male_obs_ct = ComputeR2DosagePhasedSubsetStats(dp0, dp1, founder_male_collapsed, male_dosage_invmask, sample_ct, male_ct, x_male_nm_ct0, x_male_nm_ct1, R2PhaseOmit(phase_type), sex_nmaj_dosages, &male_known_dotprod_d, &male_unknown_hethet_d, cur_nm_buf, invmask_buf);
    male_nmajsums_d[0] = u63tod(sex_nmaj_dosages[0]) * kRecipDosageMid;
    male_nmajsums_d[1] = u63tod(sex_nmaj_dosages[1]) * kRecipDosageMid;

    sex_nmaj_dosages[0] = dp0->nmaj_dosage - dp0->x_male_nmaj_dosage;
    sex_nmaj_dosages[1] = dp1->nmaj_dosage - dp1->x_male_nmaj_dosage;
    nonmale_obs_ct = ComputeR2DosagePhasedSubsetStats(dp0, dp1, founder_nonmale_collapsed, nonmale_dosage_invmask, sample_ct, sample_ct - male_ct, dp0->nm_ct - x_male_nm_ct0, dp1->nm_ct - x_male_nm_ct1, phase_type, sex_nmaj_dosages, &nonmale_known_dotprod_d, &nonmale_unknown_hethet_d, cur_nm_buf, invmask_buf);
    nonmale_nmajsums_d[0] = u63tod(sex_nmaj_dosages[0]) * kRecipDosageMid;
    nonmale_nmajsums_d[1] = u63tod(sex_nmaj_dosages[1]) * kRecipDosageMid;
  }
  if (male_obs_ct + nonmale_obs_ct == 0) {
    return -DBL_MAX;
  }
  const double male_wt = 1.0 - male_downwt;
  double nmajsums_d[2];
  nmajsums_d[0] = nonmale_nmajsums_d[0] + male_wt * male_nmajsums_d[0];
  nmajsums_d[1] = nonmale_nmajsums_d[1] + male_wt * male_nmajsums_d[1];
  const double known_dotprod_d = nonmale_known_dotprod_d + male_wt * male_known_dotprod_d;
  const double unknown_hethet_d = nonmale_unknown_hethet_d + male_wt * male_unknown_hethet_d;
  const double valid_obs_d = u31tod(nonmale_obs_ct) + male_wt * u31tod(male_obs_ct);
  const double twice_tot_recip = 0.5 / valid_obs_d;
  if ((d_ptr == nullptr) && (dprime_ptr == nullptr)) {
    double r2;
    const LDErr ld_err = PhasedLD(nmajsums_d, known_dotprod_d, unknown_hethet_d, twice_tot_recip, 0, nullptr, &r2, is_neg_ptr);
    return (ld_err == kLDErrNone)? r2 : -DBL_MAX;
  }
  double results[3];
  const LDErr ld_err = PhasedLD(nmajsums_d, known_dotprod_d, unknown_hethet_d, twice_tot_recip, 1, nullptr, results, is_neg_ptr);
  if (ld_err != kLDErrNone) {
    if (d_ptr) {
      *d_ptr = 0.0 / 0.0;
    }
    if (dprime_ptr) {
      *dprime_ptr = 0.0 / 0.0;
    }
    return -DBL_MAX;
  }
  if (d_ptr) {
    *d_ptr = results[1];
  }
  if (dprime_ptr) {
    *dprime_ptr = results[2];
  }
  return results[0];
}

void ClumpHighmemR2(uintptr_t tidx, uint32_t thread_ct_p1, uint32_t parity, ClumpCtx* ctx) {
  // Compute r^2 of non-index (variant, aidx)s against the current index
  // variant, with everything already unpacked.
  const uint64_t cur_nonindex_ct = ctx->cur_nonindex_ct;
  const uintptr_t nonindex_start = (cur_nonindex_ct * tidx) / thread_ct_p1;
  const uintptr_t nonindex_end = (cur_nonindex_ct * (tidx + 1)) / thread_ct_p1;
  uintptr_t nonindex_rem = nonindex_end - nonindex_start;
  uintptr_t* write_iter = ctx->ld_idx_found[tidx];
  if (!nonindex_rem) {
    *write_iter = ~k0LU;
    return;
  }
  const uintptr_t igroup_oaidx_start = ctx->igroup_oaidx_start;
  const uintptr_t unpacked_byte_stride = ctx->unpacked_byte_stride;
  const uintptr_t* observed_alleles = ctx->observed_alleles;
  const uintptr_t* observed_alleles_cumulative_popcounts_w = ctx->observed_alleles_cumulative_popcounts_w;
  const unsigned char* unpacked_variants = ctx->unpacked_variants;
  const uintptr_t* candidate_oabitvec = ctx->candidate_oabitvec;
  uintptr_t allele_widx_start = ctx->allele_widx_start;
  const uintptr_t allele_widx_end = ctx->allele_widx_end;
  const uintptr_t oaidx_start = ctx->a[parity].oaidx_starts[tidx];
  const uint32_t founder_ct = ctx->founder_ct;
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t allow_overlap = ctx->allow_overlap;
  const uint32_t is_x = ctx->is_x;
  const uintptr_t* founder_male_collapsed = nullptr;
  const uintptr_t* founder_nonmale_collapsed = nullptr;
  const Dosage* male_dosage_invmask = nullptr;
  const Dosage* nonmale_dosage_invmask = nullptr;
  uintptr_t* chrx_nm_buf = nullptr;
  Dosage* chrx_invmask_buf = nullptr;
  if (is_x) {
    founder_male_collapsed = ctx->founder_male_collapsed;
    male_dosage_invmask = ctx->male_dosage_invmask;
    founder_nonmale_collapsed = ctx->founder_nonmale_collapsed;
    nonmale_dosage_invmask = ctx->nonmale_dosage_invmask;
    chrx_nm_buf = ctx->chrx_workspaces[tidx];
    const uint32_t founder_ctaw = BitCtToAlignedWordCt(founder_ct);
    chrx_invmask_buf = R_CAST(Dosage*, &(chrx_nm_buf[founder_ctaw]));
  }
  const R2PhaseType phase_type = S_CAST(R2PhaseType, ctx->phase_type);
  const uint32_t load_dosage = ctx->load_dosage;
  const double r2_thresh = ctx->r2_thresh;
  const unsigned char* unpacked_index_variant = &(unpacked_variants[ctx->index_oaidx_offset * unpacked_byte_stride]);
  R2Variant index_r2v;
  FillR2V(unpacked_index_variant, founder_ct, phase_type, is_x, load_dosage, &index_r2v);
  uintptr_t oaidx_base;
  uintptr_t cur_oaidx_bits;
  BitIter1Start(candidate_oabitvec, oaidx_start, &oaidx_base, &cur_oaidx_bits);
  for (; nonindex_rem; --nonindex_rem) {
    const uintptr_t oaidx = BitIter1(candidate_oabitvec, &oaidx_base, &cur_oaidx_bits);
    const uintptr_t oaidx_offset = oaidx - igroup_oaidx_start;
    const unsigned char* unpacked_cur_variant = &(unpacked_variants[oaidx_offset * unpacked_byte_stride]);
    R2Variant cur_r2v;
    FillR2V(unpacked_cur_variant, founder_ct, phase_type, is_x, load_dosage, &cur_r2v);
    double cur_r2;
    if (!is_x) {
      uint32_t is_neg;
      cur_r2 = ComputeR2(&index_r2v, &cur_r2v, founder_ct, phase_type, load_dosage, nullptr, nullptr, &is_neg);
    } else {
      uint32_t is_neg;
      cur_r2 = ComputeXR2(&index_r2v, &cur_r2v, founder_male_collapsed, founder_nonmale_collapsed, male_dosage_invmask, nonmale_dosage_invmask, founder_ct, founder_male_ct, phase_type, load_dosage, 1, nullptr, nullptr, &is_neg, chrx_nm_buf, chrx_invmask_buf);
    }
    if (cur_r2 > r2_thresh) {
      if (!allow_overlap) {
        *write_iter = oaidx;
      } else {
        const uintptr_t allele_idx = ExpsearchIdxToUidxW(observed_alleles, observed_alleles_cumulative_popcounts_w, allele_widx_end, oaidx, &allele_widx_start);
        *write_iter = allele_idx;
      }
      ++write_iter;
    }
  }
  *write_iter = ~k0LU;
}

void ClumpLowmemR2(uintptr_t tidx, uint32_t thread_ct_p1, uint32_t parity, ClumpCtx* ctx) {
  // Unpack a few non-index PgenVariants that were read by the main thread,
  // then compute r^2 between them and the index-(variant, aidx).
  const uint64_t cur_nonindex_ct = ctx->cur_nonindex_ct;
  // yeah, this variable name sucks
  uintptr_t nonindex_idx = (cur_nonindex_ct * tidx) / thread_ct_p1;
  const uintptr_t nonindex_end = (cur_nonindex_ct * (tidx + 1)) / thread_ct_p1;
  uintptr_t* write_iter = ctx->ld_idx_found[tidx];
  if (nonindex_idx == nonindex_end) {
    *write_iter = ~k0LU;
    return;
  }
  const uintptr_t* observed_alleles = ctx->observed_alleles;
  const uintptr_t* observed_alleles_cumulative_popcounts_w = ctx->observed_alleles_cumulative_popcounts_w;
  const uintptr_t pgv_byte_stride = ctx->pgv_byte_stride;
  PgenVariant pgv = ctx->pgv_base;
  ClumpPgenVariantIncr(pgv_byte_stride * nonindex_idx, &pgv);
  const uintptr_t* candidate_oabitvec = ctx->candidate_oabitvec;
  uintptr_t allele_widx_start = ctx->allele_widx_start;
  const uintptr_t allele_widx_end = ctx->allele_widx_end;
  const uintptr_t oaidx_start = ctx->a[parity].oaidx_starts[tidx];
  uint32_t founder_ct = ctx->founder_ct;
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t allow_overlap = ctx->allow_overlap;
  const uint32_t is_x = ctx->is_x;
  const uintptr_t* founder_male_collapsed = nullptr;
  const uintptr_t* founder_nonmale_collapsed = nullptr;
  const Dosage* male_dosage_invmask = nullptr;
  const Dosage* nonmale_dosage_invmask = nullptr;
  uintptr_t* chrx_nm_buf = nullptr;
  Dosage* chrx_invmask_buf = nullptr;
  if (is_x) {
    founder_male_collapsed = ctx->founder_male_collapsed;
    founder_nonmale_collapsed = ctx->founder_nonmale_collapsed;
    male_dosage_invmask = ctx->male_dosage_invmask;
    nonmale_dosage_invmask = ctx->nonmale_dosage_invmask;
    chrx_nm_buf = ctx->chrx_workspaces[tidx];
    const uint32_t founder_ctaw = BitCtToAlignedWordCt(founder_ct);
    chrx_invmask_buf = R_CAST(Dosage*, &(chrx_nm_buf[founder_ctaw]));
  }
  const R2PhaseType phase_type = S_CAST(R2PhaseType, ctx->phase_type);
  const uint32_t load_dosage = ctx->load_dosage;
  const double r2_thresh = ctx->r2_thresh;
  unsigned char* unpacked_index_variant = ctx->unpacked_variants;
  R2Variant index_r2v;
  FillR2V(unpacked_index_variant, founder_ct, phase_type, is_x, load_dosage, &index_r2v);
  const uint32_t* phasepresent_cts = ctx->phasepresent_cts;
  const uint32_t* dosage_cts = ctx->dosage_cts;
  const uint32_t* dphase_cts = ctx->dphase_cts;
  const uintptr_t unpacked_byte_stride = ctx->unpacked_byte_stride;
  unsigned char* unpacked_cur_variant = &(unpacked_index_variant[(tidx + 1) * unpacked_byte_stride]);
  uintptr_t oaidx_base;
  uintptr_t cur_oaidx_bits;
  BitIter1Start(candidate_oabitvec, oaidx_start, &oaidx_base, &cur_oaidx_bits);
  for (; nonindex_idx != nonindex_end; ++nonindex_idx, ClumpPgenVariantIncr(pgv_byte_stride, &pgv)) {
    if (phasepresent_cts) {
      pgv.phasepresent_ct = phasepresent_cts[nonindex_idx];
    }
    if (dosage_cts) {
      pgv.dosage_ct = dosage_cts[nonindex_idx];
      if (dphase_cts) {
        pgv.dphase_ct = dphase_cts[nonindex_idx];
      }
    }
    const uintptr_t oaidx = BitIter1(candidate_oabitvec, &oaidx_base, &cur_oaidx_bits);
    if (load_dosage) {
      LdUnpackDosage(&pgv, founder_male_collapsed, male_dosage_invmask, founder_ct, phase_type, unpacked_cur_variant);
    } else {
      LdUnpackNondosageDense(&pgv, founder_male_collapsed, founder_ct, phase_type, unpacked_cur_variant);
    }
    R2Variant cur_r2v;
    FillR2V(unpacked_cur_variant, founder_ct, phase_type, is_x, load_dosage, &cur_r2v);
    double cur_r2;
    if (!is_x) {
      uint32_t is_neg;
      cur_r2 = ComputeR2(&index_r2v, &cur_r2v, founder_ct, phase_type, load_dosage, nullptr, nullptr, &is_neg);
    } else {
      uint32_t is_neg;
      cur_r2 = ComputeXR2(&index_r2v, &cur_r2v, founder_male_collapsed, founder_nonmale_collapsed, male_dosage_invmask, nonmale_dosage_invmask, founder_ct, founder_male_ct, phase_type, load_dosage, 1, nullptr, nullptr, &is_neg, chrx_nm_buf, chrx_invmask_buf);
    }
    if (cur_r2 > r2_thresh) {
      if (!allow_overlap) {
        *write_iter = oaidx;
      } else {
        const uintptr_t allele_idx = ExpsearchIdxToUidxW(observed_alleles, observed_alleles_cumulative_popcounts_w, allele_widx_end, oaidx, &allele_widx_start);
        *write_iter = allele_idx;
      }
      ++write_iter;
    }
  }
  *write_iter = ~k0LU;
}

THREAD_FUNC_DECL ClumpThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  const uint32_t calc_thread_ct_p1 = 1 + GetThreadCt(arg->sharedp);
  ClumpCtx* ctx = S_CAST(ClumpCtx*, arg->sharedp->context);
  uint32_t parity = 0;
  do {
    const ClumpJobType job_type = ctx->job_type;
    if (job_type == kClumpJobHighmemUnpack) {
      ClumpHighmemUnpack(tidx, parity, ctx);
    } else if (job_type == kClumpJobHighmemR2) {
      ClumpHighmemR2(tidx, calc_thread_ct_p1, parity, ctx);
    } else if (job_type == kClumpJobLowmemR2) {
      ClumpLowmemR2(tidx, calc_thread_ct_p1, parity, ctx);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

BoolErr ClumpSpillResults(const uintptr_t* observed_alleles, const uintptr_t* observed_alleles_cumulative_popcounts_w, uintptr_t* const* ld_idx_found, uint32_t cur_thread_ct, uintptr_t* prev_save_allele_idxp, uintptr_t* clump_sizep, uintptr_t* icandidate_oabitvec, FILE* clump_overlap_tmp) {
  uintptr_t prev_save_allele_idx = *prev_save_allele_idxp;
  uintptr_t clump_size = *clump_sizep;
  unsigned char buf[16];
  for (uint32_t tidx = 0; tidx != cur_thread_ct; ++tidx) {
    const uintptr_t* read_iter = ld_idx_found[tidx];
    for (; ; ++read_iter) {
      const uintptr_t save_allele_idx = *read_iter;
      if (save_allele_idx == ~k0LU) {
        break;
      }
      const uintptr_t oaidx = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, save_allele_idx);
      ClearBit(oaidx, icandidate_oabitvec);
      unsigned char* write_iter = Vint64Append(save_allele_idx - prev_save_allele_idx, buf);
      if (unlikely(!fwrite_unlocked(buf, 1, write_iter - buf, clump_overlap_tmp))) {
        return 1;
      }
      prev_save_allele_idx = save_allele_idx;
      ++clump_size;
    }
  }
  *prev_save_allele_idxp = prev_save_allele_idx;
  *clump_sizep = clump_size;
  return 0;
}

// 0.0001, 0.001, 0.01, 0.05 with appropriate epsilons
static const double kClumpDefaultLnBinBounds[4] = {
  -9.210340371976706,
  -6.907755278982529,
  -4.605170185988353,
  -2.995732273554161
};

static_assert(kClumpMaxBinBounds * (kMaxLnGSlen + 1) + 256 <= kMaxLongLine, "ClumpReports() needs to be updated.");
PglErr ClumpReports(const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const ClumpInfo* clump_ip, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t nosex_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  char* fname_iter = clump_ip->fnames_flattened;
  uintptr_t line_idx = 0;
  FILE* clump_overlap_tmp = nullptr;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  ClumpCtx ctx;
  TextStream txs;
  CompressStreamState css;
  ThreadGroup tg;
  PreinitTextStream(&txs);
  PreinitCstream(&css);
  PreinitThreads(&tg);
  {
    if (unlikely(founder_ct < 2)) {
      logerrputs("Error: --clump requires at least two founders.  (--make-founders may come in handy\nhere.)\n");
      goto ClumpReports_ret_INCONSISTENT_INPUT;
    } else if (founder_ct > 0x3fffffff) {
      logerrputs("Error: --clump does not support >= 2^30 founders.\n");
      goto ClumpReports_ret_NOT_YET_SUPPORTED;
    }
    uint32_t skipped_variant_ct;
    const uintptr_t* variant_include = StripUnplacedK(orig_variant_include, cip, raw_variant_ct, &skipped_variant_ct);
    if (unlikely(variant_include == nullptr)) {
      goto ClumpReports_ret_NOMEM;
    }
    const uint32_t variant_ct = orig_variant_ct - skipped_variant_ct;
    if (skipped_variant_ct) {
      logprintf("--clump: Ignoring %u chromosome 0 variant%s.\n", skipped_variant_ct, (skipped_variant_ct == 1)? "" : "s");
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uintptr_t raw_allele_ct = allele_idx_offsets? allele_idx_offsets[raw_variant_ct] : (2 * raw_variant_ct);
    const uintptr_t raw_allele_ctl = BitCtToWordCt(raw_allele_ct);
    const uintptr_t* variant_last_alidxs;
    const uint32_t* variant_last_alidxs_cumulative_popcounts;
    reterr = AllocAndFillVariantLastAlidxs(allele_idx_offsets, raw_variant_ct, max_thread_ct, &variant_last_alidxs, &variant_last_alidxs_cumulative_popcounts);
    if (unlikely(reterr)) {
      goto ClumpReports_ret_1;
    }
    uintptr_t* observed_variants;
    uintptr_t* observed_alleles;
    uintptr_t* observed_alleles_cumulative_popcounts_w;
    double* best_ln_pvals;
    // bugfix (29 Oct 2023): there are RawToSubsettedPosW(observed_alleles,
    // observed_alleles_cumulative_popcounts_w, x) calls with x =
    // raw_allele_ct.  In this case, we may need one more entry.
    if (unlikely(bigstack_calloc_w(raw_variant_ctl, &observed_variants) ||
                 bigstack_calloc_w(raw_allele_ctl, &observed_alleles) ||
                 bigstack_alloc_w(1 + (raw_allele_ct / kBitsPerWord), &observed_alleles_cumulative_popcounts_w) ||
                 bigstack_end_calloc_d(raw_allele_ct, &best_ln_pvals))) {
      goto ClumpReports_ret_NOMEM;
    }
    unsigned char* bigstack_mark2 = g_bigstack_base;
    unsigned char* bigstack_end_mark2 = g_bigstack_end;
    uint32_t file_ct = 0;
    do {
      char* fname_end = strnul(fname_iter);
      ++file_ct;
      fname_iter = &(fname_end[1]);
    } while (*fname_iter);
    assert(file_ct < 0x40000000);
    const ClumpFlags flags = clump_ip->flags;
    const uint32_t force_a1 = (flags / kfClumpForceA1) & 1;
    uint32_t* best_fidx_x2s = nullptr;
    if (((file_ct > 1) && (flags & (kfClumpColMaybeF | kfClumpColF | kfClumpColSp2))) || force_a1) {
      if (unlikely(bigstack_alloc_u32(raw_allele_ct, &best_fidx_x2s))) {
        goto ClumpReports_ret_NOMEM;
      }
    }
    const double* ln_bin_boundaries = nullptr;
    uint32_t bin_bound_ct = 0;
    if (flags & kfClumpColBins) {
      bin_bound_ct = clump_ip->bin_bound_ct;
      if (!bin_bound_ct) {
        ln_bin_boundaries = kClumpDefaultLnBinBounds;
        bin_bound_ct = 4;
      } else {
        ln_bin_boundaries = clump_ip->ln_bin_boundaries;
      }
    }
    const char* range_fname = clump_ip->range_fname;
    const uint32_t sp2_col = flags & kfClumpColSp2;
    const uint32_t ranges_col = !!range_fname;
    const uint32_t bounds_col = (flags & kfClumpColBounds) || ((flags & kfClumpColMaybeBounds) && ranges_col);
    const double ln_p1 = clump_ip->ln_p1;
    const double ln_p2 = (sp2_col || bounds_col || ranges_col)? clump_ip->ln_p2 : -DBL_MAX;
    double load_ln_pthresh = MAXV(ln_p1, ln_p2);
    if (bin_bound_ct && (load_ln_pthresh < ln_bin_boundaries[bin_bound_ct - 1])) {
      load_ln_pthresh = ln_bin_boundaries[bin_bound_ct - 1];
    }
    ClumpEntry** clump_entries = nullptr;
    uintptr_t* nonsig_arr = nullptr;
    uint32_t allow_overlap = (flags / kfClumpAllowOverlap) & 1;
    if ((flags & (kfClumpColTotal | kfClumpColBins | kfClumpColSp2)) || bounds_col || range_fname) {
      if (unlikely(BIGSTACK_ALLOC_X(ClumpEntry*, raw_allele_ct + 1, &clump_entries))) {
        goto ClumpReports_ret_NOMEM;
      }
      ZeroPtrArr(raw_allele_ct, clump_entries);
      const uint32_t nonsig_needed = ((flags & (kfClumpColTotal | kfClumpColBins)) != 0) && (load_ln_pthresh < 0.0);
      if (nonsig_needed) {
        if (unlikely(bigstack_calloc_w(raw_allele_ct, &nonsig_arr))) {
          goto ClumpReports_ret_NOMEM;
        }
      }
    } else if (allow_overlap) {
      allow_overlap = 0;
      logputs("Note: --clump-allow-overlap has no effect when --clump 'total', 'bins', and\n'sp2' column-sets are all absent.\n");
    }

    uint32_t* variant_id_htable;
    uint32_t variant_id_htable_size;
    reterr = AllocAndPopulateIdHtableMt(variant_include, variant_ids, variant_ct, bigstack_left() / 2, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
    if (unlikely(reterr)) {
      goto ClumpReports_ret_1;
    }

    const uint32_t search_a1 = force_a1 || (allele_idx_offsets && (!(flags & kfClumpNoA1)));
    const uint32_t search_test = !(flags & kfClumpNoTest);
    const char* col_search_order[4];
    col_search_order[0] = clump_ip->id_field? clump_ip->id_field : "ID\0SNP\0";
    col_search_order[1] = search_a1? (clump_ip->a1_field? clump_ip->a1_field : "A1\0") : "";
    col_search_order[2] = search_test? (clump_ip->test_field? clump_ip->test_field : "TEST\0") : "";
    const uint32_t input_log10 = flags & kfClumpInputLog10;
    col_search_order[3] = clump_ip->p_field? clump_ip->p_field : (input_log10? "LOG10_P\0NEG_LOG10_P\0P\0" : "P\0");
    const char* test_name_flattened = clump_ip->test_name? clump_ip->test_name : "ADD\0";
    const uint32_t save_all_fidxs = ((file_ct > 1) || force_a1) && sp2_col;

    LlStr* missing_variant_ids = nullptr;
    LlStr* missing_variant_allele_pairs = nullptr;
    uintptr_t missing_variant_id_ct = 0;
    uintptr_t missing_variant_id_max_slen = 0;
    uintptr_t missing_variant_allele_pair_ct = 0;
    uintptr_t missing_variant_allele_pair_max_slen = 0;
    // A1 allele normally doesn't matter for biallelic variants, so we default
    // to treating it as always-REF (may as well avoid a few +1s in the code).
    unsigned char* tmp_alloc_base = nullptr;
    unsigned char* tmp_alloc_end = bigstack_end_mark2;
    const uint32_t two_minus_force_a1 = 2 - force_a1;
    uint32_t cur_allele_ct = 2;
    uint32_t biallelic_forced_a1_alt = 0;
    for (uint32_t file_idx1 = file_ct; file_idx1; --file_idx1) {
      if (file_idx1 == 1) {
        fname_iter = clump_ip->fnames_flattened;
      } else {
        fname_iter = &(fname_iter[-3]);
        while (*fname_iter) {
          --fname_iter;
        }
        ++fname_iter;
      }
      if (file_idx1 == file_ct) {
        reterr = SizeAndInitTextStream(fname_iter, bigstack_left() / 8, MAXV(1, max_thread_ct - 1), &txs);
        tmp_alloc_base = g_bigstack_base;
      } else {
        reterr = TextRetarget(fname_iter, &txs);
      }
      if (unlikely(reterr)) {
        goto ClumpReports_ret_TSTREAM_FAIL;
      }

      line_idx = 0;
      const char* header_start;
      do {
        ++line_idx;
        header_start = TextGet(&txs);
        if (unlikely(!header_start)) {
          reterr = TextStreamRawErrcode(&txs);
          if (reterr == kPglRetEof) {
            snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", fname_iter);
            goto ClumpReports_ret_MALFORMED_INPUT_WW;
          }
          goto ClumpReports_ret_TSTREAM_FAIL;
        }
      } while (strequal_k_unsafe(header_start, "##"));
      if (*header_start == '#') {
        ++header_start;
      }

      // [0] = ID
      // [1] = A1
      // [2] = TEST
      // [3] = P
      uint32_t col_skips[4];
      uint32_t col_types[4];
      uint32_t relevant_col_ct;
      uint32_t found_type_bitset;
      reterr = SearchHeaderLine(header_start, col_search_order, "--clump", 4, &relevant_col_ct, &found_type_bitset, col_skips, col_types);
      if (unlikely(reterr)) {
        goto ClumpReports_ret_1;
      }
      if (unlikely((found_type_bitset & 0x9) != 0x9)) {
        logerrputs("Error: --clump requires ID and P columns.\n");
        goto ClumpReports_ret_INCONSISTENT_INPUT;
      }
      if (unlikely(force_a1 && (!(found_type_bitset & 0x2)))) {
        snprintf(g_logbuf, kLogbufSize, "Error: --clump-force-a1 was specified, but there is no A1 column in %s.\n", fname_iter);
        goto ClumpReports_ret_INCONSISTENT_INPUT_WW;
      }

      while (1) {
        ++line_idx;
        const char* line_start = TextGet(&txs);
        if (!line_start) {
          if (likely(!TextStreamErrcode2(&txs, &reterr))) {
            break;
          }
          goto ClumpReports_ret_TSTREAM_FAIL;
        }
        const char* token_ptrs[4];
        uint32_t token_slens[4];
        if (unlikely(!TokenLexK0(line_start, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens))) {
          goto ClumpReports_ret_MISSING_TOKENS;
        }
        if (found_type_bitset & 0x4) {
          if (!InMultistr(test_name_flattened, token_ptrs[2], token_slens[2])) {
            continue;
          }
        }
        const char* pval_str = token_ptrs[3];
        double ln_pval;
        if (!input_log10) {
          if (!ScantokLn(pval_str, &ln_pval)) {
            uint32_t cur_slen;
          ClumpReports_alphabetic_pval:
            cur_slen = token_slens[3];
            if (IsNanStr(pval_str, cur_slen)) {
              continue;
            }
            if (likely(strequal_k(pval_str, "INF", cur_slen) ||
                       (input_log10 && strequal_k(pval_str, "inf", cur_slen)))) {
              // PLINK 1.x underflow
              ln_pval = kLnNormalMin;
            } else {
              goto ClumpReports_ret_INVALID_PVAL;
            }
          }
        } else {
          double neglog10_pval;
          if (!ScantokDouble(pval_str, &neglog10_pval)) {
            goto ClumpReports_alphabetic_pval;
          }
          ln_pval = neglog10_pval * (-kLn10);
          if (unlikely(ln_pval > 0.0)) {
            goto ClumpReports_ret_INVALID_PVAL;
          }
        }
        const char* variant_id = token_ptrs[0];
        const uint32_t variant_id_slen = token_slens[0];
        const uint32_t variant_uidx = VariantIdDupflagHtableFind(token_ptrs[0], variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
        if (variant_uidx & 0x80000000U) {
          if (unlikely(variant_uidx != UINT32_MAX)) {
            snprintf(g_logbuf, kLogbufSize, "Error: --clump variant ID '%s' appears multiple times in main dataset.\n", variant_ids[variant_uidx & 0x7fffffff]);
            goto ClumpReports_ret_INCONSISTENT_INPUT_WW;
          }
          if (ln_pval <= ln_p1) {
            if (unlikely(PtrWSubCk(tmp_alloc_base, sizeof(LlStr) + RoundUpPow2(variant_id_slen + 1, sizeof(intptr_t)), &tmp_alloc_end))) {
              goto ClumpReports_ret_NOMEM;
            }
            LlStr* new_entry = R_CAST(LlStr*, tmp_alloc_end);
            new_entry->next = missing_variant_ids;
            memcpyx(new_entry->str, variant_id, variant_id_slen, '\0');
            missing_variant_ids = new_entry;
            ++missing_variant_id_ct;
            if (variant_id_slen > missing_variant_id_max_slen) {
              missing_variant_id_max_slen = variant_id_slen;
            }
          }
          continue;
        }
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        uint32_t aidx = 0;
        if (cur_allele_ct > two_minus_force_a1) {
          if (unlikely(!(found_type_bitset & 0x2))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Variant ID on line %" PRIuPTR " of %s is multiallelic, but there is no A1 column.\n", line_idx, fname_iter);
            goto ClumpReports_ret_INCONSISTENT_INPUT_WW;
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          const char* allele_code = token_ptrs[1];
          const uint32_t allele_slen = token_slens[1];
          for (; aidx != cur_allele_ct; ++aidx) {
            if (strequal_unsafe(cur_alleles[aidx], allele_code, allele_slen)) {
              break;
            }
          }
          if (aidx == cur_allele_ct) {
            if (ln_pval <= ln_p1) {
              const uint32_t variant_allele_pair_slen = variant_id_slen + allele_slen + 1;
              if (unlikely(PtrWSubCk(tmp_alloc_base, sizeof(LlStr) + RoundUpPow2(variant_allele_pair_slen + 1, sizeof(intptr_t)), &tmp_alloc_end))) {
                goto ClumpReports_ret_NOMEM;
              }
              LlStr* new_entry = R_CAST(LlStr*, tmp_alloc_end);
              new_entry->next = missing_variant_allele_pairs;
              char* write_iter = memcpyax(new_entry->str, variant_id, variant_id_slen, '\t');
              memcpyx(write_iter, allele_code, allele_slen, '\0');
              missing_variant_allele_pairs = new_entry;
              ++missing_variant_allele_pair_ct;
              if (variant_allele_pair_slen > missing_variant_allele_pair_max_slen) {
                missing_variant_allele_pair_max_slen = variant_allele_pair_slen;
              }
            }
            continue;
          }
          if (cur_allele_ct == 2) {
            biallelic_forced_a1_alt = aidx;
            aidx = 0;
          }
        }
        const uintptr_t allele_idx = allele_idx_offset_base + aidx;
        if (ln_pval > load_ln_pthresh) {
          if (unlikely(ln_pval > 0.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: p-value > 1 on line %" PRIuPTR " of %s.\n", line_idx, fname_iter);
            goto ClumpReports_ret_MALFORMED_INPUT_WW;
          }
          if (nonsig_arr && ((!bin_bound_ct) || (ln_pval > ln_bin_boundaries[bin_bound_ct - 1]))) {
            nonsig_arr[allele_idx] += 1;
            // Still need to determine which index-variant this belongs to.
            SetBit(allele_idx, observed_alleles);
            SetBit(variant_uidx, observed_variants);
          }
          continue;
        }
        const uint32_t file_idx1_x2 = (file_idx1 << 1) + biallelic_forced_a1_alt;
        // >= rather than >, to break ties in favor of file_idx1 == 1
        if (best_ln_pvals[allele_idx] >= ln_pval) {
          best_ln_pvals[allele_idx] = ln_pval;
          if (best_fidx_x2s) {
            best_fidx_x2s[allele_idx] = file_idx1_x2;
          }
        }
        SetBit(allele_idx, observed_alleles);
        SetBit(variant_uidx, observed_variants); // could defer this?
        if (clump_entries) {
          if (unlikely(PtrWSubCk(tmp_alloc_base, RoundUpPow2(sizeof(ClumpEntry), sizeof(intptr_t)), &tmp_alloc_end))) {
            goto ClumpReports_ret_NOMEM;
          }
          uint32_t pval_bin_x2 = (ln_pval > ln_p2);
          if (bin_bound_ct) {
            pval_bin_x2 |= LowerBoundNonemptyD(ln_bin_boundaries, bin_bound_ct, ln_pval) << 1;
          }
          ClumpEntry* new_entry = R_CAST(ClumpEntry*, tmp_alloc_end);
          new_entry->next = clump_entries[allele_idx];
          new_entry->pval_bin_x2 = pval_bin_x2;
          new_entry->file_idx1_x2 = file_idx1_x2;
          clump_entries[allele_idx] = new_entry;
        }
      }
    }
    if (unlikely(CleanupTextStream2(fname_iter, &txs, &reterr))) {
      goto ClumpReports_ret_1;
    }

    BigstackEndSet(tmp_alloc_end);

    FillCumulativePopcountsW(observed_alleles, raw_allele_ctl, observed_alleles_cumulative_popcounts_w);
    const uintptr_t observed_allele_ct = observed_alleles_cumulative_popcounts_w[raw_allele_ctl - 1] + PopcountWord(observed_alleles[raw_allele_ctl - 1]);
    if ((raw_allele_ct % kBitsPerWord) == 0) {
      observed_alleles_cumulative_popcounts_w[raw_allele_ctl] = observed_allele_ct;
    }
    const uint32_t output_zst = (flags / kfClumpZs) & 1;
    // Now we have efficient (variant_uidx, aidx) -> oaidx lookup.
    // Free some memory by compacting the information in best_ln_pvals, etc. to
    // exclude unused (variant_uidx, aidx) slots.
    BigstackReset(bigstack_mark2);
    unsigned char** clump_entry_varints = nullptr;
    {
      if (best_fidx_x2s) {
        const uint32_t* best_fidx_x2s_dying = best_fidx_x2s;
        best_fidx_x2s = S_CAST(uint32_t*, bigstack_alloc_raw_rd(observed_allele_ct * sizeof(int32_t)));
        uintptr_t allele_idx_base = 0;
        uintptr_t cur_bits = observed_alleles[0];
        for (uintptr_t oaidx = 0; oaidx != observed_allele_ct; ++oaidx) {
          const uintptr_t allele_idx = BitIter1(observed_alleles, &allele_idx_base, &cur_bits);
          best_fidx_x2s[oaidx] = best_fidx_x2s_dying[allele_idx];
        }
      }

      if (clump_entries) {
        ClumpEntry** clump_entries_dying = clump_entries;
        clump_entries = S_CAST(ClumpEntry**, bigstack_alloc_raw_rd((observed_allele_ct + 1) * sizeof(intptr_t)));
        uintptr_t allele_idx_base = 0;
        uintptr_t cur_bits = observed_alleles[0];
        for (uintptr_t oaidx = 0; oaidx != observed_allele_ct; ++oaidx) {
          const uintptr_t allele_idx = BitIter1(observed_alleles, &allele_idx_base, &cur_bits);
          clump_entries[oaidx] = clump_entries_dying[allele_idx];
        }

        if (nonsig_arr) {
          const uintptr_t* nonsig_arr_dying = nonsig_arr;
          nonsig_arr = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(observed_allele_ct * sizeof(intptr_t)));
          allele_idx_base = 0;
          cur_bits = observed_alleles[0];
          for (uintptr_t oaidx = 0; oaidx != observed_allele_ct; ++oaidx) {
            const uintptr_t allele_idx = BitIter1(observed_alleles, &allele_idx_base, &cur_bits);
            nonsig_arr[oaidx] = nonsig_arr_dying[allele_idx];
          }
        }
      }

      if (missing_variant_ids) {
        // natural-sort, deduplicate, and write
        bigstack_mark2 = g_bigstack_base;

        const char** strptr_arr;
        if (unlikely(bigstack_alloc_kcp(missing_variant_id_ct, &strptr_arr))) {
          goto ClumpReports_ret_NOMEM;
        }
        for (uintptr_t ulii = 0; ulii != missing_variant_id_ct; ++ulii) {
          strptr_arr[ulii] = missing_variant_ids->str;
          missing_variant_ids = missing_variant_ids->next;
        }
        assert(missing_variant_ids == nullptr);
        OutnameZstSet(".clumps.missing_id", output_zst, outname_end);
        uintptr_t nwrite;
        reterr = NsortDedupAndWrite(outname, missing_variant_id_ct, missing_variant_id_max_slen, output_zst, max_thread_ct, strptr_arr, &nwrite);
        if (unlikely(reterr)) {
          goto ClumpReports_ret_1;
        }
        logerrprintfww("Warning: %" PRIuPTR " top variant ID%s in --clump file%s missing from main dataset.  ID%s written to %s .\n", nwrite, (nwrite == 1)? "" : "s", (file_ct == 1)? "" : "s", (nwrite == 1)? "" : "s", outname);
        BigstackReset(bigstack_mark2);
      }

      if (missing_variant_allele_pairs) {
        bigstack_mark2 = g_bigstack_base;

        const char** strptr_arr;
        if (unlikely(bigstack_alloc_kcp(missing_variant_allele_pair_ct, &strptr_arr))) {
          goto ClumpReports_ret_NOMEM;
        }
        for (uintptr_t ulii = 0; ulii != missing_variant_allele_pair_ct; ++ulii) {
          strptr_arr[ulii] = missing_variant_allele_pairs->str;
          missing_variant_allele_pairs = missing_variant_allele_pairs->next;
        }
        assert(missing_variant_allele_pairs == nullptr);
        OutnameZstSet(".clumps.missing_allele", output_zst, outname_end);
        uintptr_t nwrite;
        reterr = NsortDedupAndWrite(outname, missing_variant_allele_pair_ct, missing_variant_allele_pair_max_slen, output_zst, max_thread_ct, strptr_arr, &nwrite);
        if (unlikely(reterr)) {
          goto ClumpReports_ret_1;
        }
        logerrprintfww("Warning: %" PRIuPTR " top (variant ID, A1 allele) pair%s in --clump file%s missing from main dataset due to allele rather than variant ID.  (Variant ID, A1 allele) pair%s written to %s .\n", nwrite, (nwrite == 1)? "" : "s", (file_ct == 1)? "" : "s", (nwrite == 1)? "" : "s", outname);
        BigstackReset(bigstack_mark2);
      }

      if (clump_entries) {
        // Now save pval_bin_x2 values, as well as fidxs if necessary, as
        // varints.
        unsigned char* varint_write_iter = g_bigstack_base;
        for (uintptr_t oaidx = 0; oaidx != observed_allele_ct; ++oaidx) {
          ClumpEntry* ll_iter = clump_entries[oaidx];
          // assign to clump_entries instead of clump_entry_varints, so we
          // don't break strict-aliasing rule
          clump_entries[oaidx] = R_CAST(ClumpEntry*, varint_write_iter);
          while (ll_iter) {
            varint_write_iter = Vint32Append(ll_iter->pval_bin_x2, varint_write_iter);
            if (save_all_fidxs) {
              varint_write_iter = Vint32Append(ll_iter->file_idx1_x2, varint_write_iter);
            }
            if (unlikely(varint_write_iter > tmp_alloc_end)) {
              goto ClumpReports_ret_NOMEM;
            }
            ll_iter = ll_iter->next;
          }
        }
        clump_entry_varints = R_CAST(unsigned char**, clump_entries);
        clump_entry_varints[observed_allele_ct] = varint_write_iter;
        if (unlikely(BigstackBaseSetChecked(varint_write_iter))) {
          goto ClumpReports_ret_NOMEM;
        }
        // defensive
        clump_entries = nullptr;
      }
    }
    // We're done processing all the linked lists allocated next to
    // g_bigstack_end.
    BigstackEndReset(bigstack_end_mark2);
    // Create sorted list of index-variant candidates, then free best_ln_pvals.
    uintptr_t* icandidate_vbitvec_fill;
    uintptr_t* icandidate_abitvec_fill;
    ClumpPval* index_candidates;
    if (unlikely(bigstack_calloc_w(raw_variant_ctl, &icandidate_vbitvec_fill) ||
                 bigstack_calloc_w(raw_allele_ctl, &icandidate_abitvec_fill) ||
                 BIGSTACK_ALLOC_X(ClumpPval, observed_allele_ct, &index_candidates))) {
      goto ClumpReports_ret_NOMEM;
    }
    uint32_t index_candidate_ct;
    {
      uintptr_t index_candidate_ct_w = 0;
      for (uint32_t variant_uidx = 0; ; ++variant_uidx) {
        variant_uidx = AdvBoundedTo1Bit(observed_variants, variant_uidx, raw_variant_ct);
        if (variant_uidx == raw_variant_ct) {
          break;
        }
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        const uintptr_t allele_idx_stop = allele_idx_offset_base + cur_allele_ct;
        for (uintptr_t allele_idx = allele_idx_offset_base; allele_idx != allele_idx_stop; ++allele_idx) {
          if (!IsSet(observed_alleles, allele_idx)) {
            continue;
          }
          const double cur_ln_pval = best_ln_pvals[allele_idx];
          if (cur_ln_pval <= ln_p1) {
            index_candidates[index_candidate_ct_w].ln_pval = cur_ln_pval;
            index_candidates[index_candidate_ct_w].allele_idx = allele_idx;
            ++index_candidate_ct_w;
            SetBit(variant_uidx, icandidate_vbitvec_fill);
            SetBit(allele_idx, icandidate_abitvec_fill);
          }
        }
      }
      if (!index_candidate_ct_w) {
        logerrputs("Warning: No significant --clump results.  Skipping.\n");
        goto ClumpReports_ret_1;
      }
#ifdef __LP64__
      if (unlikely(index_candidate_ct_w >= 0xffffffffU)) {
        logerrputs("Error: --clump does not support >= 2^32 - 1 index-variant candidates.\n");
        goto ClumpReports_ret_NOT_YET_SUPPORTED;
      }
#endif
      index_candidate_ct = index_candidate_ct_w;
      BigstackShrinkTop(index_candidates, sizeof(ClumpPval) * index_candidate_ct);
      // About to make a bunch of cacheline-aligned allocations at end of
      // bigstack.
      // (This is awkward, should be a single function call...)
      BigstackEndReset(bigstack_end_mark);
      bigstack_end_clalign();
      STD_SORT_PAR_UNSEQ(index_candidate_ct, ClumpPvalCmp, index_candidates);
    }
    const uintptr_t* icandidate_vbitvec = icandidate_vbitvec_fill;
    const uintptr_t* icandidate_abitvec = icandidate_abitvec_fill;
    // Main algorithm:
    // - Identify independent "islands", induced by the --clump-kb setting.
    // - Process one island (or island-group, if the island is inefficiently
    //   small) at a time.
    //   - Usually, all we need to determine is: for each remaining (variant
    //     ID, aidx) pair, what clump was it assigned to?  (clump-index is
    //     defined by index_candidates[] position.  Some clump-indexes are
    //     skipped.)
    //   - If --clump-allow-overlap is in effect, this doesn't work since
    //     non-index variants can belong to multiple clumps.  We spill a
    //     temporary file to disk tracking, for each clump, which (variant ID,
    //     aidx) pairs are a member.  This file could be very big, so we use
    //     varints and delta encoding to reduce waste.
    //     (In principle, seekable zstd should be useful here as well, but I
    //     won't investigate that unless/until I decide to use seekable zstd in
    //     other applications.)
    // - After we're done processing all islands, we have the information
    //   needed to write the final p-value sorted report.

    // Create an icandidate_idx -> index_candidates lookup table, to support
    // island-based processing.
    // "_destructive" because the main loop modifies this array in-place.
    uint32_t* icandidate_idx_to_rank0_destructive;
    {
      uint32_t* icandidate_popcounts;
      if (unlikely(bigstack_alloc_u32(index_candidate_ct, &icandidate_idx_to_rank0_destructive) ||
                   bigstack_alloc_u32(raw_allele_ctl, &icandidate_popcounts))) {
        goto ClumpReports_ret_NOMEM;
      }
      FillCumulativePopcounts(icandidate_abitvec, raw_allele_ctl, icandidate_popcounts);
      for (uint32_t rank0 = 0; rank0 != index_candidate_ct; ++rank0) {
        const uintptr_t allele_idx = index_candidates[rank0].allele_idx;
        const uint32_t icandidate_idx = RawToSubsettedPos(icandidate_abitvec, icandidate_popcounts, allele_idx);
        icandidate_idx_to_rank0_destructive[icandidate_idx] = rank0;
      }
      BigstackReset(icandidate_popcounts);
    }

    uint32_t* oallele_idx_to_clump_idx = nullptr;
    uint64_t* clump_idx_to_overlap_fpos_and_len = nullptr;
    if (!allow_overlap) {
      if (unlikely(bigstack_alloc_u32(observed_allele_ct, &oallele_idx_to_clump_idx))) {
        goto ClumpReports_ret_NOMEM;
      }
      SetAllU32Arr(observed_allele_ct, oallele_idx_to_clump_idx);
    } else {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".clumps.tmp");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &clump_overlap_tmp))) {
        goto ClumpReports_ret_OPEN_FAIL;
      }
      // make starting fpos > 0, so we know fpos == 0 marks no clump
      if (unlikely(putc_unlocked(0, clump_overlap_tmp) == EOF)) {
        goto ClumpReports_ret_WRITE_FAIL;
      }
      if (unlikely(bigstack_calloc_u64(index_candidate_ct * (2 * k1LU), &clump_idx_to_overlap_fpos_and_len))) {
        goto ClumpReports_ret_NOMEM;
      }
    }

    // We switch between two parallelization strategies, depending on
    // available memory.
    // 1. In the best case, we can load the entire set of (variant ID, aidx)
    //    pairs on at least one island into memory, and encode them in an
    //    LD-computation-friendly form.  In that case, we parallelize that
    //    operation first, and then iterate through the on-island index
    //    candidates.
    // 2. Otherwise, the main thread reads and unpacks the index (variant,
    //    aidx), then reads an affordable number of other variants at a time
    //    for the worker threads to unpack / r^2-compute with.
    //
    // We choose the worker thread count to ensure at least strategy #2 works.
    const uintptr_t observed_allele_ctl = BitCtToWordCt(observed_allele_ct);
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t founder_ctv = BitCtToVecCt(founder_ct);
    const uint32_t founder_ctv2 = NypCtToVecCt(founder_ct);
    const uint32_t founder_ctaw = founder_ctv * kWordsPerVec;
    const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
    const uint32_t founder_male_ct = PopcountWordsIntersect(founder_info, sex_male, raw_sample_ctl);
    const uint32_t founder_nonmale_ct = founder_ct - founder_male_ct;
    uint32_t founder_nosex_ct = 0;
    if (nosex_ct) {
      founder_nosex_ct = founder_ct - PopcountWordsIntersect(founder_info, sex_nm, raw_sample_ctl);
    }
    uintptr_t* candidate_oabitvec;
    uint32_t* founder_info_cumulative_popcounts;
    if (unlikely(bigstack_alloc_w(observed_allele_ctl, &candidate_oabitvec) ||
                 bigstack_alloc_u32(raw_sample_ctl, &founder_info_cumulative_popcounts))) {
      goto ClumpReports_ret_NOMEM;
    }
    SetAllBits(observed_allele_ct, candidate_oabitvec);
    uintptr_t* icandidate_oabitvec = candidate_oabitvec;
    if (allow_overlap) {
      if (unlikely(bigstack_alloc_w(observed_allele_ctl, &icandidate_oabitvec))) {
        goto ClumpReports_ret_NOMEM;
      }
      SetAllBits(observed_allele_ct, icandidate_oabitvec);
    }
    FillCumulativePopcounts(founder_info, raw_sample_ctl, founder_info_cumulative_popcounts);
    uint32_t x_code;
    if (XymtExists(cip, kChrOffsetX, &x_code)) {
      const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[x_code];
      const uint32_t x_start = cip->chr_fo_vidx_start[chr_fo_idx];
      const uint32_t x_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      if (AllBitsAreZero(icandidate_vbitvec, x_start, x_end)) {
        x_code = UINT32_MAXM1;
      }
    }
    // If founder_nonmale_ct == 0 or founder_male_ct == 0, main loop sets
    // is_x to 0 and is_haploid appropriately.
    const uint32_t x_exists = (x_code < UINT32_MAXM1) && founder_nonmale_ct && founder_male_ct;
    uint32_t y_code;
    if (XymtExists(cip, kChrOffsetY, &y_code)) {
      const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[y_code];
      const uint32_t y_start = cip->chr_fo_vidx_start[chr_fo_idx];
      const uint32_t y_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      if (AllBitsAreZero(icandidate_vbitvec, y_start, y_end)) {
        y_code = UINT32_MAXM1;
      } else if (unlikely(founder_male_ct + founder_nosex_ct == 0)) {
        // Rather not worry about this case.
        logerrputs("Error: --clump: chrY index variant(s) are present, but all founders in the main\ndataset are females.\n");
        goto ClumpReports_ret_INCONSISTENT_INPUT;
      }
    }
    // If founder_nonfemale_ct == founder_ct, we can just treat as is_y=0,
    // is_haploid=1.
    const uint32_t y_exists = (y_code < UINT32_MAXM1) && (founder_male_ct + founder_nosex_ct != founder_ct);
    const uintptr_t bitvec_byte_ct = BitCtToVecCt(founder_ct) * kBytesPerVec;
    uintptr_t* founder_male_collapsed = nullptr;
    Dosage* male_dosage_invmask = nullptr;
    uintptr_t* founder_nonmale_collapsed = nullptr;
    Dosage* nonmale_dosage_invmask = nullptr;
    if (x_exists) {
      const uint32_t founder_ctad = RoundUpPow2(founder_ct, kDosagePerVec);
      if (unlikely(bigstack_alloc_w(founder_ctl, &founder_male_collapsed) ||
                   bigstack_alloc_dosage(founder_ctad, &male_dosage_invmask) ||
                   bigstack_alloc_w(founder_ctaw, &founder_nonmale_collapsed) ||
                   bigstack_alloc_dosage(founder_ctad, &nonmale_dosage_invmask))) {
        goto ClumpReports_ret_NOMEM;
      }
      CopyBitarrSubset(sex_male, founder_info, founder_ct, founder_male_collapsed);
      BitvecInvertCopy(founder_male_collapsed, founder_ctl, founder_nonmale_collapsed);
      ZeroTrailingBits(founder_ct, founder_nonmale_collapsed);
      // potentially needed for correct founder_female_collapsed_interleaved
      // initialization
      ZeroTrailingWords(founder_ctl, founder_nonmale_collapsed);
      Expand1bitTo16(founder_male_collapsed, founder_ctad, 0xffff, male_dosage_invmask);
      Expand1bitTo16(founder_nonmale_collapsed, founder_ctad, 0xffff, nonmale_dosage_invmask);
    }
    uintptr_t* founder_female_collapsed = nullptr;
    uintptr_t* founder_female_collapsed_interleaved = nullptr;
    if (y_exists) {
      if (founder_nonmale_collapsed && (!nosex_ct)) {
        founder_female_collapsed = founder_nonmale_collapsed;
      } else {
        uintptr_t* founder_female_tmp;
        if (unlikely(bigstack_alloc_w(founder_ctaw, &founder_female_collapsed) ||
                     bigstack_alloc_w(raw_sample_ctl, &founder_female_tmp))) {
          goto ClumpReports_ret_NOMEM;
        }
        BitvecInvmaskCopy(sex_nm, sex_male, raw_sample_ctl, founder_female_tmp);
        CopyBitarrSubset(founder_female_tmp, founder_info, founder_ct, founder_female_collapsed);
        ZeroTrailingWords(founder_ctl, founder_female_collapsed);
        BigstackReset(founder_female_tmp);
      }
      if (unlikely(bigstack_alloc_w(founder_ctaw, &founder_female_collapsed_interleaved))) {
        goto ClumpReports_ret_NOMEM;
      }
      FillInterleavedMaskVec(founder_female_collapsed, founder_ctv, founder_female_collapsed_interleaved);
    }
    const uint32_t check_dosage = (pgfip->gflags / kfPgenGlobalDosagePresent) & 1;
    uintptr_t dosagevec_byte_ct = 0;
    if (check_dosage) {
      dosagevec_byte_ct = DivUp(founder_ct, kDosagePerVec) * kBytesPerVec;
    }

    uint32_t calc_thread_ct = MAXV(1, max_thread_ct - 1);
    // no big deal if these are slightly overallocated
    if (unlikely(BIGSTACK_ALLOC_X(PgenReader*, calc_thread_ct, &ctx.pgr_ptrs) ||
                 bigstack_alloc_w(calc_thread_ct + 2, &(ctx.a[0].oaidx_starts)) ||
                 bigstack_alloc_w(calc_thread_ct + 2, &(ctx.a[1].oaidx_starts)) ||
                 bigstack_alloc_wp(calc_thread_ct + 1, &(ctx.ld_idx_found)))) {
      goto ClumpReports_ret_NOMEM;
    }
    ctx.raregeno_bufs = nullptr;
    ctx.difflist_sample_id_bufs = nullptr;
    ctx.chrx_workspaces = nullptr;
    const uint32_t phased_r2 = !(flags & kfClumpUnphased);
    if (!phased_r2) {
      if (unlikely(bigstack_alloc_wp(calc_thread_ct, &ctx.raregeno_bufs) ||
                   bigstack_alloc_u32p(calc_thread_ct, &ctx.difflist_sample_id_bufs))) {
        goto ClumpReports_ret_NOMEM;
      }
    }
    if (x_exists) {
      if (unlikely(bigstack_alloc_wp(calc_thread_ct + 1, &ctx.chrx_workspaces))) {
        goto ClumpReports_ret_NOMEM;
      }
    }
    // Now determine real calc_thread_ct.
    // This function's logic is sufficiently different from the usual
    // PgfiMultiread use case that we fork PgenMtLoadInit()'s computation here
    // instead of modifying that function.

    const uint32_t all_haploid = IsSet(cip->haploid_mask, 0);
    const uint32_t check_phase = phased_r2 && (!all_haploid) && (pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent));
    PgenGlobalFlags effective_gflags = pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
    if (!check_phase) {
      effective_gflags &= kfPgenGlobalDosagePresent;
    }

    if (unlikely(BigstackAllocPgv(founder_ct, 0, effective_gflags, &ctx.pgv_base))) {
      goto ClumpReports_ret_NOMEM;
    }
    const uintptr_t pgv_byte_stride = g_bigstack_base - R_CAST(unsigned char*, ctx.pgv_base.genovec);

    const uint32_t max_difflist_len = founder_ct / 64;
    uintptr_t unpacked_byte_stride;
    if (check_dosage) {
      const uintptr_t dosage_trail_byte_ct = LdDosageTrailAlignedByteCt(S_CAST(R2PhaseType, phased_r2), x_exists);
      unpacked_byte_stride = dosagevec_byte_ct * (1 + phased_r2 + check_phase) + bitvec_byte_ct + dosage_trail_byte_ct;
    } else {
      unpacked_byte_stride = RoundUpPow2(16, kBytesPerVec) + bitvec_byte_ct * (3 + 2 * check_phase);
#ifndef USE_AVX2
      const uintptr_t sparse_req = RoundUpPow2((6 + max_difflist_len) * sizeof(int32_t), kBytesPerVec) + NypCtToVecCt(max_difflist_len) * kBytesPerVec + RoundUpPow2(founder_ctl * (kBytesPerWord + sizeof(int32_t)), kBytesPerVec);
      if (sparse_req > unpacked_byte_stride) {
        unpacked_byte_stride = sparse_req;
      }
#endif
      const uintptr_t nondosage_trail_byte_ct = LdNondosageTrailAlignedByteCt(S_CAST(R2PhaseType, phased_r2), x_exists);
      unpacked_byte_stride += nondosage_trail_byte_ct;
    }
    const uintptr_t unpacked_byte_stride_cachealign = RoundUpPow2(unpacked_byte_stride, kCacheline);
    const uintptr_t pgr_struct_alloc = RoundUpPow2(sizeof(PgenReader), kCacheline);

    // Haven't counted PgenReader instance, two unpacked_variants slots, or
    // ld_idx_found slots required by main thread yet.

    // FillGaussianDArr() uses a minimum per-thread job size of ~4 MiB of
    // memory writes.  I'm guessing that is also a reasonable unpacked_variants
    // shard size to aim for.
    const uint32_t min_pgv_per_thread = 1 + 4194303 / pgv_byte_stride;
    uintptr_t sparse_alloc = 0;
    uint32_t raregeno_word_ct = 0;
    uintptr_t chrx_alloc = 0;
    if (x_exists) {
      const uint32_t larger_half = MAXV(founder_male_ct, founder_nonmale_ct);
      chrx_alloc = RoundUpPow2(BitCtToVecCt(larger_half) + larger_half * sizeof(Dosage), kCacheline);
    }
    {
      const uintptr_t min_ld_idx_found_alloc = WordCtToCachelineCt(min_pgv_per_thread + 1) * kCacheline;
      const uintptr_t min_u32_alloc = Int32CtToCachelineCt(min_pgv_per_thread) * kCacheline;
      const uintptr_t phasepresent_alloc = check_phase * min_u32_alloc;
      const uintptr_t dosage_ct_alloc = check_dosage * min_u32_alloc;
      const uintptr_t dphase_ct_alloc = dosage_ct_alloc * check_phase;
      if (!phased_r2) {
        const uintptr_t raregeno_vec_ct = DivUp(max_difflist_len, kNypsPerVec);
        const uintptr_t difflist_sample_id_vec_ct = DivUp(max_difflist_len, kInt32PerVec);
        sparse_alloc = RoundUpPow2((raregeno_vec_ct + difflist_sample_id_vec_ct) * kBytesPerVec, kCacheline);
        raregeno_word_ct = raregeno_vec_ct * kWordsPerVec;
      }

      uintptr_t bytes_avail = bigstack_left();
      const uintptr_t more_base_alloc = pgr_struct_alloc + pgr_alloc_cacheline_ct * kCacheline + 2 * unpacked_byte_stride_cachealign + min_ld_idx_found_alloc + phasepresent_alloc + dosage_ct_alloc + dphase_ct_alloc + chrx_alloc;
      if (unlikely(bytes_avail < more_base_alloc)) {
        goto ClumpReports_ret_NOMEM;
      }
      bytes_avail -= more_base_alloc;

      const uintptr_t per_thread_target_alloc = 2 * min_pgv_per_thread * pgv_byte_stride + unpacked_byte_stride_cachealign + pgr_struct_alloc + pgr_alloc_cacheline_ct * kCacheline + min_ld_idx_found_alloc + sparse_alloc + chrx_alloc;
      if (bytes_avail < per_thread_target_alloc * calc_thread_ct) {
        calc_thread_ct = bytes_avail / per_thread_target_alloc;
        if (unlikely(!calc_thread_ct)) {
          goto ClumpReports_ret_NOMEM;
        }
      }
    }
    // Ready to initialize the rest of ctx.
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto ClumpReports_ret_NOMEM;
    }

    // Effective end of pgv-buffer allocation, in highmem case.
    g_bigstack_base += calc_thread_ct * pgv_byte_stride;

    // Make this non-null, because PgrInit() branches on that.  However, exact
    // value doesn't matter here.
    pgfip->block_base = g_bigstack_base;
    for (uint32_t tidx = 0; tidx <= calc_thread_ct; ++tidx) {
      ctx.pgr_ptrs[tidx] = S_CAST(PgenReader*, bigstack_end_alloc_raw(pgr_struct_alloc));
      unsigned char* pgr_alloc = S_CAST(unsigned char*, bigstack_end_alloc_raw(pgr_alloc_cacheline_ct * kCacheline));

      // shouldn't be possible for this to fail
      PgrInit(nullptr, 0, pgfip, ctx.pgr_ptrs[tidx], pgr_alloc);

      if (!phased_r2) {
        uintptr_t* raregeno_buf = S_CAST(uintptr_t*, bigstack_end_alloc_raw(sparse_alloc));
        ctx.raregeno_bufs[tidx] = raregeno_buf;
        ctx.difflist_sample_id_bufs[tidx] = R_CAST(uint32_t*, &(raregeno_buf[raregeno_word_ct]));
      }
      if (x_exists) {
        ctx.chrx_workspaces[tidx] = S_CAST(uintptr_t*, bigstack_end_alloc_raw(chrx_alloc));
      }
    }

    unsigned char* multiread_base[2];
    multiread_base[0] = g_bigstack_base;
    multiread_base[1] = nullptr;
    // Decide on a maximum value upfront, rather than reconfiguring this for
    // each island-group.
    uintptr_t multiread_byte_target;
    {
      // Allow up to 1/8 of remaining workspace, going lower if there's always
      // enough memory for a full Vblock (65536 variants).
      const uint32_t observed_variant_ct = PopcountWords(observed_variants, raw_variant_ctl);
      const uint64_t vblock_based_cacheline_ct_limit = PgfiMultireadGetCachelineReq(observed_variants, pgfip, observed_variant_ct, kPglVblockSize);
      const uint64_t proportional_cacheline_ct_limit = bigstack_left() / (8 * kCacheline);
      // Stick to lowmem mode if this rule can't even guarantee 1/500th of a
      // Vblock.
      if (vblock_based_cacheline_ct_limit <= 500 * proportional_cacheline_ct_limit) {
        multiread_byte_target = kCacheline * MINV(vblock_based_cacheline_ct_limit, proportional_cacheline_ct_limit);
      } else {
        // (this value should guarantee lowmem mode is always selected)
        multiread_byte_target = (~k0LU) >> 1;
      }
    }

    ctx.observed_variants = observed_variants;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.variant_last_alidxs = variant_last_alidxs;
    ctx.variant_last_alidxs_cumulative_popcounts = variant_last_alidxs_cumulative_popcounts;
    ctx.observed_alleles = observed_alleles;
    ctx.observed_alleles_cumulative_popcounts_w = observed_alleles_cumulative_popcounts_w;
    ctx.founder_info = founder_info;
    ctx.founder_info_cumulative_popcounts = founder_info_cumulative_popcounts;
    ctx.founder_male_collapsed = founder_male_collapsed;
    ctx.male_dosage_invmask = male_dosage_invmask;
    ctx.founder_nonmale_collapsed = founder_nonmale_collapsed;
    ctx.nonmale_dosage_invmask = nonmale_dosage_invmask;
    ctx.founder_female_collapsed = founder_female_collapsed;
    ctx.founder_female_collapsed_interleaved = founder_female_collapsed_interleaved;
    ctx.founder_ct = founder_ct;
    ctx.founder_male_ct = founder_male_ct;
    ctx.pgv_byte_stride = pgv_byte_stride;
    ctx.bitvec_byte_ct = bitvec_byte_ct;
    ctx.dosagevec_byte_ct = dosagevec_byte_ct;
    ctx.r2_thresh = clump_ip->r2;
    ctx.allow_overlap = allow_overlap;
    ctx.candidate_oabitvec = candidate_oabitvec;
    ctx.err_info = (~0LLU) << 32;

    unsigned char* lowmem_unpacked_variants = &(g_bigstack_end[(2 + calc_thread_ct) * (-S_CAST(intptr_t, unpacked_byte_stride_cachealign))]);
    uint32_t* phasepresent_cts = nullptr;
    uint32_t* dosage_cts = nullptr;
    uint32_t* dphase_cts = nullptr;
    uintptr_t* lowmem_ld_idx_found_base;
    uintptr_t max_lowmem_nonindex_ct;
    {
      // In lowmem mode, split remaining memory between
      // ctx.lowmem_pgv_base, phasepresent_cts, dosage_cts, dphase_cts, and
      // ld_idx_found.
      unsigned char* alloc_iter = R_CAST(unsigned char*, ctx.pgv_base.genovec);
      const uintptr_t bytes_avail = lowmem_unpacked_variants - alloc_iter;
      const uintptr_t u32_field_ct = check_phase + check_dosage + check_dosage * check_phase;
      const uintptr_t per_variant_cost = pgv_byte_stride + sizeof(intptr_t) + sizeof(int32_t) * u32_field_ct;
      max_lowmem_nonindex_ct = bytes_avail / per_variant_cost;
      while (max_lowmem_nonindex_ct * pgv_byte_stride + RoundUpPow2(max_lowmem_nonindex_ct * sizeof(intptr_t), kCacheline) + u32_field_ct * RoundUpPow2(max_lowmem_nonindex_ct * sizeof(int32_t), kCacheline) > bytes_avail) {
        // pgv_byte_stride is a positive multiple of kCacheline, and we lose
        // less than 4 cachelines to adverse rounding, so this will exit
        // quickly enough.
        --max_lowmem_nonindex_ct;
      }
      const uintptr_t pgv_alloc = max_lowmem_nonindex_ct * pgv_byte_stride;
      alloc_iter = &(alloc_iter[pgv_alloc]);
      const uintptr_t u32_alloc = RoundUpPow2(max_lowmem_nonindex_ct * sizeof(int32_t), kCacheline);
      if (check_phase) {
        phasepresent_cts = R_CAST(uint32_t*, alloc_iter);
        alloc_iter = &(alloc_iter[u32_alloc]);
      }
      if (check_dosage) {
        dosage_cts = R_CAST(uint32_t*, alloc_iter);
        alloc_iter = &(alloc_iter[u32_alloc]);
        if (check_phase) {
          dphase_cts = R_CAST(uint32_t*, alloc_iter);
          alloc_iter = &(alloc_iter[u32_alloc]);
        }
      }
      lowmem_ld_idx_found_base = R_CAST(uintptr_t*, alloc_iter);
    }
    ctx.phasepresent_cts = phasepresent_cts;
    ctx.dosage_cts = dosage_cts;
    ctx.dphase_cts = dphase_cts;

    SetThreadFuncAndData(ClumpThread, &ctx, &tg);

    // Main loop:
    // 1. Identify boundaries of next island.  (An island is a set of variants
    //    that is distant enough from the rest of the variants that it can be
    //    processed independently.)
    // 2. Determine highest memory mode that still fits.
    // 3. Holding the memory mode constant, check if additional small islands
    //    can be appended to the group.
    // 4. Process the island-group.
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(founder_info_cumulative_popcounts, simple_pgrp, &pssi);
    uintptr_t max_overlap_clump_size = 0; // does not count self
    uint32_t clump_ct = 0;
    const uintptr_t bytes_avail = bigstack_left();
    const uint32_t bp_radius = clump_ip->bp_radius;
    uint32_t parity = 0;
    uint32_t next_vidx_end = 0;
    uint32_t next_chr_fo_idx = 0;
    uintptr_t next_allele_idx_first;
    uintptr_t next_allele_idx_last;
    uintptr_t next_oaidx_start;
    uintptr_t next_oaidx_end;
    uint32_t next_vidx_start = GetNextIslandIdxs(cip, variant_bps, allele_idx_offsets, icandidate_vbitvec, observed_variants, observed_alleles, observed_alleles_cumulative_popcounts_w, raw_variant_ct, bp_radius, &next_oaidx_start, &next_oaidx_end, &next_vidx_end, &next_chr_fo_idx, &next_allele_idx_first, &next_allele_idx_last);
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t is_x = 0;
    uint32_t is_y = 0;
    uint32_t is_haploid = 0;
    for (uint32_t icandidate_idx_start = 0; icandidate_idx_start < index_candidate_ct; ) {
      uint32_t vidx_start = next_vidx_start;
      uint32_t vidx_end = next_vidx_end;
      if (chr_fo_idx != next_chr_fo_idx) {
        chr_fo_idx = next_chr_fo_idx;
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        is_haploid = IsSet(cip->haploid_mask, chr_idx);
        is_x = (chr_idx == x_code);
        is_y = (chr_idx == y_code);
        if (is_x) {
          if (!founder_nonmale_ct) {
            is_x = 0;
          } else if (!founder_male_ct) {
            is_x = 0;
            is_haploid = 0;
          }
        }
      }

      const uintptr_t allele_idx_first = next_allele_idx_first;
      uintptr_t allele_idx_last = next_allele_idx_last;
      uintptr_t oaidx_start = next_oaidx_start;
      uintptr_t oaidx_end = next_oaidx_end;
      uintptr_t candidate_ct = oaidx_end - oaidx_start;
      uint32_t relevant_phase_exists = 0;
      uint32_t load_dosage = 0;
      ScanPhaseDosage(observed_variants, pgfip, vidx_start, vidx_end, check_phase && ((!is_haploid) || is_x), check_dosage, &relevant_phase_exists, &load_dosage);
      R2PhaseType phase_type = GetR2PhaseType(phased_r2, relevant_phase_exists);
      ctx.is_x = is_x;
      ctx.is_y = is_y;
      ctx.igroup_oaidx_start = oaidx_start;
      ctx.allele_widx_start = allele_idx_first / kBitsPerWord;

      uint64_t unpacked_variant_byte_stride = UnpackedByteStride(&ctx, phase_type, x_exists, load_dosage);

      next_vidx_end = vidx_end;
      next_chr_fo_idx = chr_fo_idx;
      next_vidx_start = GetNextIslandIdxs(cip, variant_bps, allele_idx_offsets, icandidate_vbitvec, observed_variants, observed_alleles, observed_alleles_cumulative_popcounts_w, raw_variant_ct, bp_radius, &next_oaidx_start, &next_oaidx_end, &next_vidx_end, &next_chr_fo_idx, &next_allele_idx_first, &next_allele_idx_last);

      // Highmem requirement:
      //   candidate_ct * unpacked_variant_byte_stride for unpacked variants
      //   max(multiread_byte_target for PgfiMultiread,
      //       (candidate_ct + calc_thread_ct) * sizeof(intptr_t) for
      //       ld_idx_found)
      const uintptr_t high_result_byte_req = RoundUpPow2((candidate_ct + calc_thread_ct) * sizeof(intptr_t), 2 * kCacheline);
      uintptr_t cur_high_multiread_byte_target = MAXV(high_result_byte_req / 2, multiread_byte_target);
      uint32_t icandidate_ct;
      if (bytes_avail >= candidate_ct * unpacked_variant_byte_stride + 2 * cur_high_multiread_byte_target) {
        // highmem loop
        // Don't continue trying to extend if we've reached the end of the
        // chromosome, or the current island-group already appears to be large
        // enough for good worker-thread utilization.
        while ((next_chr_fo_idx == chr_fo_idx) && (candidate_ct * unpacked_variant_byte_stride < (4194304 * k1LU) * calc_thread_ct)) {
          uint32_t ext_relevant_phase_exists = relevant_phase_exists;
          uint32_t ext_load_dosage = load_dosage;
          ScanPhaseDosage(observed_variants, pgfip, next_vidx_start, next_vidx_end, check_phase && ((!is_haploid) || is_x), check_dosage, &ext_relevant_phase_exists, &ext_load_dosage);
          const R2PhaseType ext_phase_type = GetR2PhaseType(phased_r2, ext_relevant_phase_exists);
          const uintptr_t ext_candidate_ct = candidate_ct + next_oaidx_end - next_oaidx_start;
          const uint64_t ext_unpacked_variant_byte_stride = UnpackedByteStride(&ctx, ext_phase_type, x_exists, ext_load_dosage);
          const uintptr_t ext_high_multiread_byte_target = MAXV(RoundUpPow2((ext_candidate_ct + calc_thread_ct) * sizeof(intptr_t), 2 * kCacheline) / 2, multiread_byte_target);
          if (bytes_avail < ext_candidate_ct * S_CAST(uint64_t, ext_unpacked_variant_byte_stride) + 2 * ext_high_multiread_byte_target) {
            // Insufficient memory to extend.
            break;
          }
          vidx_end = next_vidx_end;
          allele_idx_last = next_allele_idx_last;
          candidate_ct = ext_candidate_ct;
          unpacked_variant_byte_stride = ext_unpacked_variant_byte_stride;
          cur_high_multiread_byte_target = ext_high_multiread_byte_target;
          phase_type = ext_phase_type;
          load_dosage = ext_load_dosage;
          next_vidx_start = GetNextIslandIdxs(cip, variant_bps, allele_idx_offsets, icandidate_vbitvec, observed_variants, observed_alleles, observed_alleles_cumulative_popcounts_w, raw_variant_ct, bp_radius, nullptr, &next_oaidx_end, &next_vidx_end, &next_chr_fo_idx, &next_allele_idx_first, &next_allele_idx_last);
        }
        multiread_base[1] = &(g_bigstack_base[cur_high_multiread_byte_target]);
        unsigned char* unpacked_variants = &(multiread_base[1][cur_high_multiread_byte_target]);
        pgfip->block_base = multiread_base[parity];
        icandidate_ct = PopcountBitRange(icandidate_abitvec, allele_idx_first, allele_idx_last + 1);
        ctx.unpacked_variants = unpacked_variants;
        ctx.unpacked_byte_stride = unpacked_variant_byte_stride;
        ctx.phase_type = phase_type;
        ctx.load_dosage = load_dosage;
        ctx.allele_widx_end = 1 + (allele_idx_last / kBitsPerWord);
        ctx.job_type = kClumpJobHighmemUnpack;
        uint32_t multiread_vidx_start = vidx_start;
        const uint64_t fpos_end = GetPgfiFpos(pgfip, vidx_end);
        while (1) {
          // Similar to MultireadNonempty(), but we don't (i) expect to start
          // on a block boundary or (ii) expect to continue all the way to
          // the end of the file.
          const uint64_t fpos_start = GetPgfiLdbaseFpos(pgfip, multiread_vidx_start);
          uint32_t multiread_vidx_end = vidx_end;
          if (fpos_end - fpos_start > cur_high_multiread_byte_target) {
            if (!pgfip->var_fpos) {
              multiread_vidx_end = multiread_vidx_start + cur_high_multiread_byte_target / pgfip->const_vrec_width;
            } else {
              multiread_vidx_end = LastLeqU64(pgfip->var_fpos, multiread_vidx_start, raw_variant_ct, fpos_start + cur_high_multiread_byte_target);
            }
            multiread_vidx_end = 1 + FindLast1BitBefore(observed_variants, multiread_vidx_end);
          }
          const uint32_t load_variant_ct = PopcountBitRange(observed_variants, multiread_vidx_start, multiread_vidx_end);
          reterr = PgfiMultiread(observed_variants, multiread_vidx_start, multiread_vidx_end, load_variant_ct, pgfip);
          if (unlikely(reterr)) {
            goto ClumpReports_ret_PGR_FAIL;
          }
          uintptr_t multiread_allele_idx_start;
          uintptr_t multiread_allele_idx_end;
          if (allele_idx_offsets) {
            multiread_allele_idx_start = allele_idx_offsets[multiread_vidx_start];
            multiread_allele_idx_end = allele_idx_offsets[multiread_vidx_end];
          } else {
            multiread_allele_idx_start = 2 * multiread_vidx_start;
            multiread_allele_idx_end = 2 * multiread_vidx_end;
          }
          const uintptr_t multiread_oaidx_start = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, multiread_allele_idx_start);
          const uintptr_t multiread_oaidx_end = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, multiread_allele_idx_end);
          FillWStarts(calc_thread_ct, multiread_oaidx_start, multiread_oaidx_end - multiread_oaidx_start, ctx.a[parity].oaidx_starts);
          ctx.a[parity].oaidx_starts[calc_thread_ct] = multiread_oaidx_end;
          if (multiread_vidx_start != vidx_start) {
            JoinThreads(&tg);
            reterr = S_CAST(PglErr, ctx.err_info);
            if (unlikely(reterr)) {
              PgenErrPrintNV(reterr, ctx.err_info >> 32);
              goto ClumpReports_ret_1;
            }
          }

          PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
          if (unlikely(SpawnThreads(&tg))) {
            goto ClumpReports_ret_THREAD_CREATE_FAIL;
          }
          parity = 1 - parity;

          if (multiread_vidx_end == vidx_end) {
            break;
          }
          multiread_vidx_start = AdvTo1Bit(observed_variants, multiread_vidx_end);
          pgfip->block_base = multiread_base[parity];
        }
        STD_SORT(icandidate_ct, u32cmp, &(icandidate_idx_to_rank0_destructive[icandidate_idx_start]));
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, ctx.err_info);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, ctx.err_info >> 32);
          goto ClumpReports_ret_1;
        }
        // Now iterate through index variants.
        ctx.job_type = kClumpJobHighmemR2;
        for (uint32_t icandidate_idx = 0; icandidate_idx != icandidate_ct; ++icandidate_idx) {
          if ((icandidate_idx_start + icandidate_idx) % 1000 == 0) {
            printf("\r--clump: %u/%u index candidate%s processed.", icandidate_idx_start + icandidate_idx, index_candidate_ct, (index_candidate_ct == 1)? "" : "s");
            fflush(stdout);
          }
          const uint32_t rank0 = icandidate_idx_to_rank0_destructive[icandidate_idx_start + icandidate_idx];
          const uintptr_t allele_idx = index_candidates[rank0].allele_idx;
          const uintptr_t index_oaidx = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, allele_idx);
          if (!IsSet(icandidate_oabitvec, index_oaidx)) {
            // Already included in another clump.
            continue;
          }
          ++clump_ct;
          uint32_t variant_uidx;
          if (!variant_last_alidxs) {
            variant_uidx = allele_idx / 2;
          } else {
            variant_uidx = RawToSubsettedPos(variant_last_alidxs, variant_last_alidxs_cumulative_popcounts, allele_idx);
          }
          uint32_t window_start_vidx = vidx_start;
          uint32_t window_end_vidx = vidx_end;
          GetBpWindow(observed_variants, variant_bps, variant_uidx, bp_radius, &window_start_vidx, &window_end_vidx);
          uintptr_t window_allele_idx_start;
          uintptr_t window_allele_idx_end;
          if (!allele_idx_offsets) {
            window_allele_idx_start = window_start_vidx * 2;
            window_allele_idx_end = window_end_vidx * 2;
          } else {
            window_allele_idx_start = allele_idx_offsets[window_start_vidx];
            window_allele_idx_end = allele_idx_offsets[window_end_vidx];
          }
          uintptr_t window_oaidx_start = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, window_allele_idx_start);
          window_oaidx_start = AdvTo1Bit(candidate_oabitvec, window_oaidx_start);
          uintptr_t window_oaidx_end = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, window_allele_idx_end);
          window_oaidx_end = 1 + FindLast1BitBefore(candidate_oabitvec, window_oaidx_end);
          // Wait till this point to clear the bit, to simplify
          // window_oaidx_start and window_oaidx_end initialization.
          ClearBit(index_oaidx, candidate_oabitvec);
          if (oallele_idx_to_clump_idx) {
            oallele_idx_to_clump_idx[index_oaidx] = rank0;
          } else {
            const uint64_t fpos = ftello(clump_overlap_tmp);
            clump_idx_to_overlap_fpos_and_len[rank0 * 2] = fpos;
          }
          uintptr_t nonindex_ct = PopcountBitRange(candidate_oabitvec, window_oaidx_start, window_oaidx_end);
          if (!nonindex_ct) {
            // No remaining non-index variants to clump with this.
            continue;
          }
          // If there's exactly one non-index variant to check, there's no
          // point in waking up the worker threads.
          // Possible todo: if allow_overlap is false, each worker thread could
          // operate independently.  Ideally, the number of worker threads
          // assigned to a single index variant is limited by the job size.
          // (However, in the allow_overlap true case, maybe just tune this
          // threshold and leave everything else the same.)
          const uint32_t cur_thread_ct = (nonindex_ct == 1)? 1 : (calc_thread_ct + 1);
          FillWSubsetStarts(candidate_oabitvec, cur_thread_ct, window_oaidx_start, nonindex_ct, ctx.a[parity].oaidx_starts);
          ctx.index_oaidx_offset = index_oaidx - oaidx_start;
          uintptr_t* ld_idx_found_base = R_CAST(uintptr_t*, multiread_base[0]);
          ctx.cur_nonindex_ct = nonindex_ct;
          ctx.ld_idx_found[0] = ld_idx_found_base;
          for (uint32_t tidx = 1; tidx != cur_thread_ct; ++tidx) {
            const uintptr_t offset = (tidx * S_CAST(uint64_t, nonindex_ct)) / (cur_thread_ct);
            // These are (~k0LU)-terminated sequences of oaidxs (usual case) or
            // allele_idxs (--clump-allow-overlap case).  Leave enough room for
            // all possible indexes to be included, plus the terminator.
            ctx.ld_idx_found[tidx] = &(ld_idx_found_base[offset + tidx]);
          }
          if (cur_thread_ct > 1) {
            SpawnThreads(&tg);
          }
          ClumpHighmemR2(cur_thread_ct - 1, cur_thread_ct, parity, &ctx);
          if (cur_thread_ct > 1) {
            JoinThreads(&tg);
          }
          if (oallele_idx_to_clump_idx) {
            for (uint32_t tidx = 0; tidx != cur_thread_ct; ++tidx) {
              const uintptr_t* read_iter = ctx.ld_idx_found[tidx];
              for (; ; ++read_iter) {
                const uintptr_t oaidx = *read_iter;
                if (oaidx == ~k0LU) {
                  break;
                }
                ClearBit(oaidx, candidate_oabitvec);
                oallele_idx_to_clump_idx[oaidx] = rank0;
              }
            }
          } else {
            uintptr_t prev_save_allele_idx = 0;
            uintptr_t clump_size = 0;
            if (unlikely(ClumpSpillResults(observed_alleles, observed_alleles_cumulative_popcounts_w, ctx.ld_idx_found, cur_thread_ct, &prev_save_allele_idx, &clump_size, icandidate_oabitvec, clump_overlap_tmp))) {
              goto ClumpReports_ret_WRITE_FAIL;
            }
            if (clump_size > max_overlap_clump_size) {
              max_overlap_clump_size = clump_size;
            }
            clump_idx_to_overlap_fpos_and_len[rank0 * 2 + 1] = ftello(clump_overlap_tmp) - clump_idx_to_overlap_fpos_and_len[rank0 * 2];
          }
          if (cur_thread_ct > 1) {
            parity = 1 - parity;
          }
        }
      } else {
        icandidate_ct = PopcountBitRange(icandidate_abitvec, allele_idx_first, allele_idx_last + 1);
        ctx.unpacked_variants = lowmem_unpacked_variants;
        ctx.unpacked_byte_stride = unpacked_byte_stride_cachealign;
        ctx.phase_type = phase_type;
        ctx.load_dosage = load_dosage;
        ctx.allele_widx_end = 1 + (allele_idx_last / kBitsPerWord);
        ctx.job_type = kClumpJobLowmemR2;
        STD_SORT(icandidate_ct, u32cmp, &(icandidate_idx_to_rank0_destructive[icandidate_idx_start]));
        const uintptr_t* cur_founder_male_collapsed = is_x? founder_male_collapsed : nullptr;
        const Dosage* cur_male_dosage_invmask = is_x? male_dosage_invmask : nullptr;
        const uint32_t icandidate_idx_end = icandidate_idx_start + icandidate_ct;
        for (uint32_t icandidate_idx = icandidate_idx_start; icandidate_idx != icandidate_idx_end; ++icandidate_idx) {
          if (icandidate_idx % 1000 == 0) {
            printf("\r--clump: %u/%u index candidate%s processed.", icandidate_idx, index_candidate_ct, (index_candidate_ct == 1)? "" : "s");
            fflush(stdout);
          }
          const uint32_t rank0 = icandidate_idx_to_rank0_destructive[icandidate_idx];
          const uintptr_t index_allele_idx = index_candidates[rank0].allele_idx;
          const uintptr_t index_oaidx = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, index_allele_idx);
          if (!IsSet(icandidate_oabitvec, index_oaidx)) {
            // Already included in another clump.
            continue;
          }
          ++clump_ct;
          uint32_t index_variant_uidx;
          AlleleCode index_aidx;
          if (!variant_last_alidxs) {
            index_variant_uidx = index_allele_idx / 2;
            index_aidx = index_allele_idx % 2;
          } else {
            index_variant_uidx = RawToSubsettedPos(variant_last_alidxs, variant_last_alidxs_cumulative_popcounts, index_allele_idx);
            index_aidx = index_allele_idx - allele_idx_offsets[index_variant_uidx];
          }
          uint32_t window_start_vidx = vidx_start;
          uint32_t window_end_vidx = vidx_end;
          GetBpWindow(observed_variants, variant_bps, index_variant_uidx, bp_radius, &window_start_vidx, &window_end_vidx);
          uintptr_t window_allele_idx_start;
          uintptr_t window_allele_idx_end;
          if (!allele_idx_offsets) {
            window_allele_idx_start = window_start_vidx * 2;
            window_allele_idx_end = window_end_vidx * 2;
          } else {
            window_allele_idx_start = allele_idx_offsets[window_start_vidx];
            window_allele_idx_end = allele_idx_offsets[window_end_vidx];
          }
          uintptr_t window_oaidx_cur = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, window_allele_idx_start);
          window_oaidx_cur = AdvTo1Bit(candidate_oabitvec, window_oaidx_cur);
          uintptr_t window_oaidx_end = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, window_allele_idx_end);
          window_oaidx_end = 1 + FindLast1BitBefore(candidate_oabitvec, window_oaidx_end);
          // Wait till this point to clear the bit, to simplify
          // window_oaidx_cur and window_oaidx_end initialization.
          ClearBit(index_oaidx, candidate_oabitvec);
          if (oallele_idx_to_clump_idx) {
            oallele_idx_to_clump_idx[index_oaidx] = rank0;
          } else {
            const uint64_t fpos = ftello(clump_overlap_tmp);
            clump_idx_to_overlap_fpos_and_len[rank0 * 2] = fpos;
          }
          uintptr_t rem_nonindex_ct = PopcountBitRange(candidate_oabitvec, window_oaidx_cur, window_oaidx_end);
          if (!rem_nonindex_ct) {
            // No remaining non-index variants to clump with this.
            continue;
          }
          const uintptr_t allele_widx_end = DivUp(window_allele_idx_end, kBitsPerWord);
          uintptr_t allele_widx_start = window_allele_idx_start / kBitsPerWord;
          uintptr_t prev_save_allele_idx = 0;
          uintptr_t clump_size = 0;
          uint32_t index_variant_needed = 1;
          do {
            const uintptr_t cur_nonindex_ct = MINV(rem_nonindex_ct, max_lowmem_nonindex_ct);
            const uint32_t cur_thread_ct = (cur_nonindex_ct == 1)? 1 : (calc_thread_ct + 1);
            FillWSubsetStarts(candidate_oabitvec, cur_thread_ct, window_oaidx_cur, cur_nonindex_ct, ctx.a[parity].oaidx_starts);

            uintptr_t next_window_oaidx_cur = window_oaidx_end;
            if (cur_nonindex_ct < rem_nonindex_ct) {
              const uintptr_t last_bit_ct = cur_nonindex_ct - (cur_nonindex_ct * S_CAST(uint64_t, cur_thread_ct - 1)) / cur_thread_ct;
              next_window_oaidx_cur = FindNth1BitFrom(candidate_oabitvec, ctx.a[parity].oaidx_starts[cur_thread_ct - 1], last_bit_ct + 1);
            }

            // Load PgenVariants serially.
            PgenVariant pgv = ctx.pgv_base;
            uintptr_t oaidx_base;
            uintptr_t cur_oaidx_bits;
            BitIter1Start(candidate_oabitvec, window_oaidx_cur, &oaidx_base, &cur_oaidx_bits);
            for (uintptr_t nonindex_idx = 0; nonindex_idx != cur_nonindex_ct; ) {
              uintptr_t variant_uidx;
              AlleleCode aidx;
              if (index_variant_needed) {
                variant_uidx = index_variant_uidx;
                aidx = index_aidx;
              } else {
                const uintptr_t oaidx = BitIter1(candidate_oabitvec, &oaidx_base, &cur_oaidx_bits);
                const uintptr_t allele_idx = ExpsearchIdxToUidxW(observed_alleles, observed_alleles_cumulative_popcounts_w, allele_widx_end, oaidx, &allele_widx_start);
                if (!allele_idx_offsets) {
                  variant_uidx = allele_idx / 2;
                  aidx = allele_idx % 2;
                } else {
                  variant_uidx = RawToSubsettedPos(variant_last_alidxs, variant_last_alidxs_cumulative_popcounts, allele_idx);
                  aidx = allele_idx - allele_idx_offsets[variant_uidx];
                }
              }
              if (load_dosage) {
                if (phase_type == kR2PhaseTypePresent) {
                  reterr = PgrGetInv1Dp(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, &pgv);
                } else {
                  reterr = PgrGetInv1D(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
                }
                if (is_y) {
                  InterleavedSetMissingCleardosage(founder_female_collapsed, founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec, &pgv.dosage_ct, pgv.dosage_present, pgv.dosage_main);
                }
              } else {
                if (phase_type == kR2PhaseTypePresent) {
                  reterr = PgrGetInv1P(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
                } else {
                  reterr = PgrGetInv1(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec);
                }
                if (is_y) {
                  InterleavedSetMissing(founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec);
                }
              }
              if (unlikely(reterr)) {
                goto ClumpReports_ret_PGR_FAIL;
              }
              if (!index_variant_needed) {
                ClumpPgenVariantIncr(pgv_byte_stride, &pgv);
                if (phasepresent_cts) {
                  phasepresent_cts[nonindex_idx] = pgv.phasepresent_ct;
                }
                if (dosage_cts) {
                  dosage_cts[nonindex_idx] = pgv.dosage_ct;
                }
                if (dphase_cts) {
                  dphase_cts[nonindex_idx] = pgv.dphase_ct;
                }
                ++nonindex_idx;
              } else {
                if (load_dosage) {
                  LdUnpackDosage(&pgv, cur_founder_male_collapsed, cur_male_dosage_invmask, founder_ct, phase_type, lowmem_unpacked_variants);
                } else {
                  LdUnpackNondosageDense(&pgv, cur_founder_male_collapsed, founder_ct, phase_type, lowmem_unpacked_variants);
                }
                index_variant_needed = 0;
              }
            }

            // Compute r^2s for PgenVariants against index variant in parallel.
            ctx.cur_nonindex_ct = cur_nonindex_ct;
            ctx.ld_idx_found[0] = lowmem_ld_idx_found_base;
            for (uint32_t tidx = 1; tidx != cur_thread_ct; ++tidx) {
              const uintptr_t offset = (tidx * S_CAST(uint64_t, cur_nonindex_ct)) / cur_thread_ct;
              ctx.ld_idx_found[tidx] = &(lowmem_ld_idx_found_base[offset + tidx]);
            }
            if (cur_thread_ct > 1) {
              if (unlikely(SpawnThreads(&tg))) {
                goto ClumpReports_ret_THREAD_CREATE_FAIL;
              }
            }
            ClumpLowmemR2(cur_thread_ct - 1, cur_thread_ct, parity, &ctx);
            if (cur_thread_ct > 1) {
              JoinThreads(&tg);
            }
            if (oallele_idx_to_clump_idx) {
              for (uint32_t tidx = 0; tidx != cur_thread_ct; ++tidx) {
                const uintptr_t* read_iter = ctx.ld_idx_found[tidx];
                for (; ; ++read_iter) {
                  const uintptr_t oaidx = *read_iter;
                  if (oaidx == ~k0LU) {
                    break;
                  }
                  ClearBit(oaidx, candidate_oabitvec);
                  oallele_idx_to_clump_idx[oaidx] = rank0;
                }
              }
            } else {
              if (unlikely(ClumpSpillResults(observed_alleles, observed_alleles_cumulative_popcounts_w, ctx.ld_idx_found, cur_thread_ct, &prev_save_allele_idx, &clump_size, icandidate_oabitvec, clump_overlap_tmp))) {
                goto ClumpReports_ret_WRITE_FAIL;
              }
            }
            window_oaidx_cur = next_window_oaidx_cur;
            rem_nonindex_ct -= cur_nonindex_ct;
            if (cur_thread_ct > 1) {
              parity = 1 - parity;
            }
          } while (rem_nonindex_ct);
          if (clump_idx_to_overlap_fpos_and_len) {
            if (clump_size > max_overlap_clump_size) {
              max_overlap_clump_size = clump_size;
            }
            clump_idx_to_overlap_fpos_and_len[rank0 * 2 + 1] = ftello(clump_overlap_tmp) - clump_idx_to_overlap_fpos_and_len[rank0 * 2];
          }
        }
      }

      icandidate_idx_start += icandidate_ct;
    }
    ctx.job_type = kClumpJobNone;
    DeclareLastThreadBlock(&tg);
    if (unlikely(SpawnThreads(&tg))) {
      goto ClumpReports_ret_THREAD_CREATE_FAIL;
    }
    JoinThreads(&tg);

    fputs("\r", stdout);
    logprintf("--clump: %u clump%s formed from %u index candidate%s.", clump_ct, (clump_ct == 1)? "" : "s", index_candidate_ct, (index_candidate_ct == 1)? "" : "s");
    fputs("  ", stdout);
    logputs("\n");

    // Usual (i.e. not --clump-allow-overlap) postprocessing algorithm:
    // 1. Count the number of (variant, aidx)s associated with each clump.
    // 2. Allocate ordered_members[]; partial count-sums give each clump's
    //    starting position within ordered_members[].
    // 3. Fill ordered_members.
    BigstackDoubleReset(candidate_oabitvec, bigstack_end_mark);
    uintptr_t* clump_ends = nullptr;
    uintptr_t* ordered_members = nullptr;
    unsigned char* overlap_raw_loadbuf = nullptr;
    uintptr_t* overlap_allele_idxs = nullptr;
    if (!allow_overlap) {
      uintptr_t* clump_sizes;
      if (unlikely(bigstack_calloc_w(index_candidate_ct, &clump_sizes))) {
        goto ClumpReports_ret_NOMEM;
      }
      uintptr_t running_total = 0;
      // could parallelize this, but doesn't realistically matter
      for (uintptr_t oaidx = 0; oaidx != observed_allele_ct; ++oaidx) {
        const uint32_t clump_idx = oallele_idx_to_clump_idx[oaidx];
        if (clump_idx != UINT32_MAX) {
          clump_sizes[clump_idx] += 1;
        }
      }
      // convert clump_sizes to clump_write_offsets
      for (uint32_t clump_idx = 0; clump_idx != index_candidate_ct; ++clump_idx) {
        const uintptr_t cur_clump_size = clump_sizes[clump_idx];
        clump_sizes[clump_idx] = running_total;
        running_total += cur_clump_size;
      }
      uintptr_t* clump_write_offsets = clump_sizes;
      if (unlikely(bigstack_alloc_w(running_total, &ordered_members))) {
        goto ClumpReports_ret_NOMEM;
      }
      // Now fill ordered_members, convert oaidxs to allele_idxs, and convert
      // clump_write_offsets to clump_ends as a side-effect.
      uintptr_t allele_idx_base = 0;
      uintptr_t cur_bits = observed_alleles[0];
      for (uintptr_t oaidx = 0; oaidx != observed_allele_ct; ++oaidx) {
        const uintptr_t allele_idx = BitIter1(observed_alleles, &allele_idx_base, &cur_bits);
        const uint32_t clump_idx = oallele_idx_to_clump_idx[oaidx];
        if (clump_idx != UINT32_MAX) {
          const uintptr_t write_offset = clump_write_offsets[clump_idx];
          ordered_members[write_offset] = allele_idx;
          clump_write_offsets[clump_idx] += 1;
        }
      }
      clump_ends = clump_write_offsets;
    } else {
      uint64_t max_vint_byte_ct = 0;
      for (uint32_t clump_idx = 0; clump_idx != index_candidate_ct; ++clump_idx) {
        const uint64_t vint_byte_ct = clump_idx_to_overlap_fpos_and_len[clump_idx * 2 + 1];
        if (vint_byte_ct > max_vint_byte_ct) {
          max_vint_byte_ct = vint_byte_ct;
        }
      }
      if (unlikely(bigstack_alloc64_uc(max_vint_byte_ct, &overlap_raw_loadbuf) ||
                   bigstack_alloc_w(max_overlap_clump_size + 1, &overlap_allele_idxs))) {
        goto ClumpReports_ret_NOMEM;
      }
      if (unlikely(fclose_null(&clump_overlap_tmp))) {
        goto ClumpReports_ret_WRITE_FAIL;
      }
      if (unlikely(fopen_checked(outname, FOPEN_RB, &clump_overlap_tmp))) {
        goto ClumpReports_ret_OPEN_FAIL;
      }
    }

    char* range_group_names = nullptr;
    uintptr_t* rg_chr_bounds = nullptr;
    uint32_t** rg_setdefs = nullptr;
    uintptr_t max_range_group_id_blen = 0;
    if (range_fname) {
      uintptr_t ignored_group_ct;
      uintptr_t ignored_chr_max_group_ct;
      reterr = LoadAndSortIntervalBed(range_fname, cip, nullptr, (flags / kfClumpRange0) & 1, clump_ip->range_border, 0, 0, max_thread_ct, &ignored_group_ct, &range_group_names, &max_range_group_id_blen, &rg_chr_bounds, &rg_setdefs, &ignored_chr_max_group_ct);
      if (unlikely(reterr)) {
        goto ClumpReports_ret_1;
      }
    }

    const uint32_t overflow_buf_size = kCompressStreamBlock + MAXV(kMaxIdSlen + max_variant_id_slen + max_allele_slen, bin_bound_ct * (kMaxLnGSlen + 1)) + 256;
    OutnameZstSet(".clumps", output_zst, outname_end);
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto ClumpReports_ret_1;
    }
    const uint32_t chr_col = flags & kfClumpColChrom;
    const uint32_t pos_col = flags & kfClumpColPos;
    const uint32_t ref_col = flags & kfClumpColRef;
    const uint32_t alt1_col = flags & kfClumpColAlt1;
    const uint32_t alt_col = flags & kfClumpColAlt;
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t all_nonref = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags, flags / kfClumpColMaybeprovref, raw_variant_ct, all_nonref);
    const uint32_t a1_col = (flags & kfClumpColA1) || ((flags & kfClumpColMaybeA1) && MultiallelicVariantPresent(variant_include, allele_idx_offsets, variant_ct));
    const uint32_t f_col = (flags & kfClumpColF) || ((flags & kfClumpColMaybeF) && (file_ct > 1));
    const uint32_t f_in_sp2 = sp2_col && ((flags & kfClumpColF) || (file_ct > 1));
    const uint32_t output_log10 = flags & kfClumpOutputLog10;
    const uint32_t total_col = flags & kfClumpColTotal;
    uintptr_t* cur_bin_counts = nullptr;
    if (bin_bound_ct) {
      if (unlikely(bigstack_alloc_w(bin_bound_ct + 1, &cur_bin_counts))) {
        goto ClumpReports_ret_NOMEM;
      }
    }
    *cswritep++ = '#';
    if (chr_col) {
      cswritep = strcpya_k(cswritep, "CHROM\t");
    }
    if (pos_col) {
      cswritep = strcpya_k(cswritep, "POS\t");
    }
    cswritep = strcpya_k(cswritep, "ID\t");
    if (ref_col) {
      cswritep = strcpya_k(cswritep, "REF\t");
    }
    if (alt1_col) {
      cswritep = strcpya_k(cswritep, "ALT1\t");
    }
    if (alt_col) {
      cswritep = strcpya_k(cswritep, "ALT\t");
    }
    if (provref_col) {
      cswritep = strcpya_k(cswritep, "PROVISIONAL_REF?\t");
    }
    if (a1_col) {
      cswritep = strcpya_k(cswritep, "A1\t");
    }
    if (f_col) {
      cswritep = strcpya_k(cswritep, "F\t");
    }
    if (output_log10) {
      cswritep = strcpya_k(cswritep, "NEG_LOG10_");
    }
    *cswritep++ = 'P';
    if (total_col) {
      cswritep = strcpya_k(cswritep, "\tTOTAL");
    }
    if (bounds_col) {
      cswritep = strcpya_k(cswritep, "\tCLUMP_FIRST_POS\tCLUMP_LAST_POS");
    }
    if (bin_bound_ct) {
      cswritep = strcpya_k(cswritep, "\tNONSIG");
      for (uint32_t bin_idx = bin_bound_ct; bin_idx; ) {
        --bin_idx;
        cswritep = strcpya_k(cswritep, "\tS");
        cswritep = lntoa_g(ln_bin_boundaries[bin_idx], cswritep);
      }
    }
    if (sp2_col) {
      cswritep = strcpya_k(cswritep, "\tSP2");
    }
    if (ranges_col) {
      cswritep = strcpya_k(cswritep, "\tRANGES");
    }
    AppendBinaryEoln(&cswritep);
    if (unlikely(Cswrite(&css, &cswritep))) {
      goto ClumpReports_ret_WRITE_FAIL;
    }

    uintptr_t* clump_allele_idxs = overlap_allele_idxs;
    uintptr_t prev_clump_end = 0;
    uint32_t chr_idx = 0;
    uint32_t index_allele_ct = 2;
    uint32_t index_file_idx1 = 1;
    uint32_t file_idx1 = 1;
    uint32_t first_bp = 0;
    uint32_t last_bp = 0;
    for (uint32_t clump_idx = 0; clump_idx != index_candidate_ct; ++clump_idx) {
      const uintptr_t index_allele_idx = index_candidates[clump_idx].allele_idx;
      uintptr_t clump_size_including_self;
      if (!allow_overlap) {
        const uintptr_t cur_clump_end = clump_ends[clump_idx];
        if (cur_clump_end == prev_clump_end) {
          continue;
        }
        clump_allele_idxs = &(ordered_members[prev_clump_end]);
        clump_size_including_self = cur_clump_end - prev_clump_end;
        prev_clump_end = cur_clump_end;
      } else {
        const uint64_t fpos = clump_idx_to_overlap_fpos_and_len[clump_idx * 2];
        if (fpos == 0) {
          continue;
        }
        const uint64_t vint_byte_ct = clump_idx_to_overlap_fpos_and_len[clump_idx * 2 + 1];
        if (vint_byte_ct) {
          if (unlikely(fseeko(clump_overlap_tmp, fpos, SEEK_SET) ||
                       fread_checked(overlap_raw_loadbuf, vint_byte_ct, clump_overlap_tmp))) {
            goto ClumpReports_ret_READ_FAIL;
          }
        }
        const unsigned char* read_iter = overlap_raw_loadbuf;
        const unsigned char* read_end = &(read_iter[vint_byte_ct]);
        uintptr_t* write_iter = clump_allele_idxs;
        uintptr_t last_allele_idx = 0;
        while (read_iter != read_end) {
          const uintptr_t allele_idx_incr = GetVint64Unsafe(&read_iter);
          last_allele_idx += allele_idx_incr;
          if (last_allele_idx > index_allele_idx) {
            break;
          }
          *write_iter++ = last_allele_idx;
        }
        *write_iter++ = index_allele_idx;
        if (last_allele_idx > index_allele_idx) {
          *write_iter++ = last_allele_idx;
        }
        while (read_iter != read_end) {
          const uintptr_t allele_idx_incr = GetVint64Unsafe(&read_iter);
          last_allele_idx += allele_idx_incr;
          *write_iter++ = last_allele_idx;
        }
        clump_size_including_self = write_iter - clump_allele_idxs;
      }

      uintptr_t index_allele_idx_offset_base;
      uint32_t index_variant_uidx;
      if (!allele_idx_offsets) {
        index_variant_uidx = index_allele_idx / 2;
        index_allele_idx_offset_base = index_allele_idx & (~k1LU);
      } else {
        index_variant_uidx = RawToSubsettedPos(variant_last_alidxs, variant_last_alidxs_cumulative_popcounts, index_allele_idx);
        index_allele_idx_offset_base = allele_idx_offsets[index_variant_uidx];
        index_allele_ct = allele_idx_offsets[index_variant_uidx + 1] - index_allele_idx_offset_base;
      }
      if (chr_col || ranges_col) {
        chr_idx = GetVariantChr(cip, index_variant_uidx);
        if (chr_col) {
          cswritep = chrtoa(cip, chr_idx, cswritep);
          *cswritep++ = '\t';
        }
      }
      if (pos_col) {
        cswritep = u32toa_x(variant_bps[index_variant_uidx], '\t', cswritep);
      }
      cswritep = strcpyax(cswritep, variant_ids[index_variant_uidx], '\t');
      if (ref_col) {
        cswritep = strcpyax(cswritep, allele_storage[index_allele_idx_offset_base], '\t');
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto ClumpReports_ret_WRITE_FAIL;
        }
      }
      if (alt1_col) {
        cswritep = strcpyax(cswritep, allele_storage[index_allele_idx_offset_base + 1], '\t');
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto ClumpReports_ret_WRITE_FAIL;
        }
      }
      if (alt_col) {
        for (uint32_t aidx = 1; aidx != index_allele_ct; ++aidx) {
          cswritep = strcpya(cswritep, allele_storage[index_allele_idx_offset_base + aidx]);
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto ClumpReports_ret_WRITE_FAIL;
          }
          *cswritep++ = ',';
        }
        cswritep[-1] = '\t';
      }
      if (provref_col) {
        *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, index_variant_uidx)))? 'Y' : 'N';
        *cswritep++ = '\t';
      }
      const uintptr_t index_oaidx = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, index_allele_idx);
      if (a1_col) {
        if (index_allele_ct == 2) {
          if (force_a1) {
            const uint32_t aidx = best_fidx_x2s[index_oaidx] & 1;
            cswritep = strcpyax(cswritep, allele_storage[index_allele_idx_offset_base + aidx], '\t');
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto ClumpReports_ret_WRITE_FAIL;
            }
          } else {
            cswritep = strcpya_k(cswritep, ".\t");
          }
        } else {
          cswritep = strcpyax(cswritep, allele_storage[index_allele_idx], '\t');
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto ClumpReports_ret_WRITE_FAIL;
          }
        }
      }
      if (best_fidx_x2s) {
        index_file_idx1 = best_fidx_x2s[index_oaidx] >> 1;
      }
      if (f_col) {
        cswritep = u32toa_x(index_file_idx1, '\t', cswritep);
      }
      const double index_ln_pval = index_candidates[clump_idx].ln_pval;
      if (!output_log10) {
        const double reported_ln = MAXV(index_ln_pval, output_min_ln);
        cswritep = lntoa_g(reported_ln, cswritep);
      } else {
        const double reported_val = (-kRecipLn10) * index_ln_pval;
        cswritep = dtoa_g(reported_val, cswritep);
      }
      if (total_col || bin_bound_ct) {
        uintptr_t total_ct = 0;
        if (bin_bound_ct) {
          ZeroWArr(bin_bound_ct + 1, cur_bin_counts);
          for (uintptr_t member_idx = 0; member_idx != clump_size_including_self; ++member_idx) {
            const uintptr_t cur_allele_idx = clump_allele_idxs[member_idx];
            const uintptr_t cur_oaidx = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, cur_allele_idx);
            if (nonsig_arr) {
              cur_bin_counts[bin_bound_ct] += nonsig_arr[cur_oaidx];
            }
            const unsigned char* varint_read_iter = clump_entry_varints[cur_oaidx];
            const unsigned char* varint_read_end = clump_entry_varints[cur_oaidx + 1];
            while (varint_read_iter != varint_read_end) {
              const uint32_t pval_bin = GetVint32Unsafe(&varint_read_iter) >> 1;
              cur_bin_counts[pval_bin] += 1;
              if (save_all_fidxs) {
                SkipVintUnsafe(&varint_read_iter);
              }
            }
          }
          // Keep appearances of this (variant, uidx) in other files, but don't
          // count the central appearance.
          const uint32_t index_pval_bin = LowerBoundNonemptyD(ln_bin_boundaries, bin_bound_ct, index_ln_pval);
          cur_bin_counts[index_pval_bin] -= 1;

          for (uint32_t bin_idx = 0; bin_idx <= bin_bound_ct; ++bin_idx) {
            total_ct += cur_bin_counts[bin_idx];
          }
        } else {
          for (uintptr_t member_idx = 0; member_idx != clump_size_including_self; ++member_idx) {
            const uintptr_t cur_allele_idx = clump_allele_idxs[member_idx];
            const uintptr_t cur_oaidx = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, cur_allele_idx);
            if (nonsig_arr) {
              total_ct += nonsig_arr[cur_oaidx];
            }
            const unsigned char* varint_read_start = clump_entry_varints[cur_oaidx];
            const unsigned char* varint_read_end = clump_entry_varints[cur_oaidx + 1];
            const uintptr_t varint_ct = CountVints(varint_read_start, varint_read_end);
            total_ct += varint_ct >> save_all_fidxs;
          }
          --total_ct;
        }
        if (total_col) {
          *cswritep++ = '\t';
          cswritep = wtoa(total_ct, cswritep);
        }
      }
      if (bounds_col || ranges_col) {
        for (intptr_t direction = 1; direction != -3; direction -= 2) {
          uint32_t pval_too_high = 1;
          uintptr_t member_idx = (direction == 1)? 0 : (clump_size_including_self - 1);
          for (; member_idx != clump_size_including_self ; member_idx += direction) {
            const uintptr_t cur_allele_idx = clump_allele_idxs[member_idx];
            const uintptr_t cur_oaidx = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, cur_allele_idx);
            const unsigned char* varint_read_iter = clump_entry_varints[cur_oaidx];
            const unsigned char* varint_read_end = clump_entry_varints[cur_oaidx + 1];
            while (varint_read_iter != varint_read_end) {
              pval_too_high = GetVint32Unsafe(&varint_read_iter) & 1;
              if (!pval_too_high) {
                break;
              }
            }
            if (!pval_too_high) {
              uint32_t variant_uidx;
              if (!allele_idx_offsets) {
                variant_uidx = cur_allele_idx / 2;
              } else {
                variant_uidx = RawToSubsettedPos(variant_last_alidxs, variant_last_alidxs_cumulative_popcounts, cur_allele_idx);
              }
              if (direction == 1) {
                first_bp = variant_bps[variant_uidx];
              } else {
                last_bp = variant_bps[variant_uidx];
              }
              break;
            }
          }
          if (member_idx == clump_size_including_self) {
            // special case: no --clump-p2 hits at all, not even index variant
            assert(ln_p1 > ln_p2);
            first_bp = UINT32_MAX;
            break;
          }
        }
        if (bounds_col) {
          *cswritep++ = '\t';
          if (first_bp != UINT32_MAX) {
            cswritep = u32toa_x(first_bp, '\t', cswritep);
            cswritep = u32toa(last_bp, cswritep);
          } else {
            cswritep = strcpya_k(cswritep, ".\t.");
          }
        }
      }
      if (bin_bound_ct) {
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto ClumpReports_ret_WRITE_FAIL;
        }
        for (uint32_t bin_idx = bin_bound_ct + 1; bin_idx; ) {
          --bin_idx;
          *cswritep++ = '\t';
          cswritep = wtoa(cur_bin_counts[bin_idx], cswritep);
        }
      }
      if (sp2_col) {
        *cswritep++ = '\t';
        uint32_t nonempty = 0;
        for (uintptr_t member_idx = 0; member_idx != clump_size_including_self; ++member_idx) {
          const uintptr_t cur_allele_idx = clump_allele_idxs[member_idx];
          const uintptr_t cur_oaidx = RawToSubsettedPosW(observed_alleles, observed_alleles_cumulative_popcounts_w, cur_allele_idx);
          const unsigned char* varint_read_iter = clump_entry_varints[cur_oaidx];
          const unsigned char* varint_read_end = clump_entry_varints[cur_oaidx + 1];
          while (varint_read_iter != varint_read_end) {
            const uint32_t pval_too_high = GetVint32Unsafe(&varint_read_iter) & 1;
            if (!pval_too_high) {
              if (save_all_fidxs) {
                const uint32_t file_idx1_x2 = GetVint32Unsafe(&varint_read_iter);
                file_idx1 = file_idx1_x2 >> 1;
                biallelic_forced_a1_alt = file_idx1_x2 & 1;
              }
              if ((cur_allele_idx != index_allele_idx) || (file_idx1 != index_file_idx1)) {
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto ClumpReports_ret_WRITE_FAIL;
                }
                nonempty = 1;
                uintptr_t cur_allele_idx_offset_base;
                uint32_t cur_variant_uidx;
                if (!allele_idx_offsets) {
                  cur_variant_uidx = cur_allele_idx / 2;
                  cur_allele_idx_offset_base = cur_allele_idx & (~k1LU);
                } else {
                  cur_variant_uidx = RawToSubsettedPos(variant_last_alidxs, variant_last_alidxs_cumulative_popcounts, cur_allele_idx);
                  cur_allele_idx_offset_base = allele_idx_offsets[cur_variant_uidx];
                  cur_allele_ct = allele_idx_offsets[cur_variant_uidx + 1] - cur_allele_idx_offset_base;
                }
                cswritep = strcpya(cswritep, variant_ids[cur_variant_uidx]);
                if (force_a1 || (cur_allele_ct > 2)) {
                  *cswritep++ = '(';
                  cswritep = strcpyax(cswritep, allele_storage[cur_allele_idx + biallelic_forced_a1_alt], ')');
                }
                if (f_in_sp2) {
                  *cswritep++ = '(';
                  cswritep = u32toa_x(file_idx1, ')', cswritep);
                }
                *cswritep++ = ',';
              }
            } else {
              if (save_all_fidxs) {
                SkipVintUnsafe(&varint_read_iter);
              }
            }
          }
        }
        if (nonempty) {
          --cswritep;
        } else {
          *cswritep++ = '.';
        }
      }
      if (ranges_col) {
        *cswritep++ = '\t';
        uint32_t nonempty = 0;
        if (first_bp != UINT32_MAX) {
          const uintptr_t cur_rg_start_idx = rg_chr_bounds[chr_idx];
          uint32_t** cur_rg_setdefs = &(rg_setdefs[cur_rg_start_idx]);
          const char* cur_rg_names = &(range_group_names[cur_rg_start_idx * max_range_group_id_blen + kMaxChrCodeDigits]);
          const uintptr_t cur_rg_ct = rg_chr_bounds[chr_idx + 1] - cur_rg_start_idx;
          const uint32_t end_bp = last_bp + 1;
          for (uintptr_t rg_idx = 0; rg_idx != cur_rg_ct; ++rg_idx) {
            if (IntervalInSetdef(cur_rg_setdefs[rg_idx], first_bp, end_bp)) {
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto ClumpReports_ret_WRITE_FAIL;
              }
              nonempty = 1;
              cswritep = strcpyax(cswritep, &(cur_rg_names[rg_idx * max_range_group_id_blen]), ',');
            }
          }
        }
        if (nonempty) {
          --cswritep;
        } else {
          *cswritep++ = '.';
        }
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto ClumpReports_ret_WRITE_FAIL;
      }
    }

    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto ClumpReports_ret_WRITE_FAIL;
    }
    logprintf("Results written to %s .\n", outname);
    if (clump_overlap_tmp) {
      if (unlikely(fclose_null(&clump_overlap_tmp))) {
        goto ClumpReports_ret_READ_FAIL;
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".clumps.tmp");
      if (unlikely(unlink(outname))) {
        goto ClumpReports_ret_WRITE_FAIL;
      }
    }
  }
  while (0) {
  ClumpReports_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ClumpReports_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ClumpReports_ret_READ_FAIL:
    if (feof_unlocked(clump_overlap_tmp)) {
      errno = 0;
    }
    logerrprintfww(kErrprintfFread, "--clump-allow-overlap temporary file", rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  ClumpReports_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ClumpReports_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ClumpReports_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--clump file", &txs);
    break;
  ClumpReports_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  ClumpReports_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, fname_iter);
  ClumpReports_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  ClumpReports_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ClumpReports_ret_INVALID_PVAL:
    logerrprintfww("Error: Invalid p-value on line %" PRIuPTR " of %s.\n", line_idx, fname_iter);
    reterr = kPglRetInconsistentInput;
    break;
  ClumpReports_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  ClumpReports_ret_NOT_YET_SUPPORTED:
    reterr = kPglRetNotYetSupported;
    break;
  }
 ClumpReports_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&css, cswritep);
  CleanupTextStream2("--clump file", &txs, &reterr);
  fclose_cond(clump_overlap_tmp);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

// indexes here are all subsetted (founder_idx / variant_idx), not sample_uidx
// / variant_uidx.
typedef struct VcorMatrixCtxStruct {
  // Shared constants.
  const uintptr_t* founder_male_collapsed;
  const uintptr_t* founder_nonmale_collapsed;
  // invmask is the most useful representation, since we have a bunch of
  // dosage-processing functions which interpret 65535 as missing.
  const Dosage* male_dosage_invmask;
  const Dosage* nonmale_dosage_invmask;
  const uint32_t* chr_fo_idx_end;
  uint32_t chr_ct;
  uint32_t chrx_idx; // UINT32_MAX if not present
  uint32_t founder_ct;
  uint32_t founder_male_ct;
  unsigned char is_unsquared;
  unsigned char triangle_calc; // true for both square0 and triangle
  unsigned char phase_type;
  unsigned char check_dosage;
  uint32_t variant_ct;
  uintptr_t unpacked_variant_byte_stride;

  // Input data.
  unsigned char* unpacked_row_variants[2];  // read from [row_parity]
  unsigned char* unpacked_col_variants[2];  // read from [col_parity]
  uint32_t cur_row_variant_idx_start;
  uint32_t row_window_size;
  uint32_t cur_col_variant_idx_start;
  uint32_t col_window_size;

  // per-thread chrX workspaces.
  uintptr_t** cur_nm_bufs;
  Dosage** invmask_bufs;

  // Output double-buffer.  write to [row_parity]
  double* results_d[2];
  float* results_f[2];
} VcorMatrixCtx;

THREAD_FUNC_DECL VcorMatrixThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  VcorMatrixCtx* ctx = S_CAST(VcorMatrixCtx*, arg->sharedp->context);
  const uintptr_t* founder_male_collapsed = ctx->founder_male_collapsed;
  const uintptr_t* founder_nonmale_collapsed = ctx->founder_nonmale_collapsed;
  const Dosage* male_dosage_invmask = ctx->male_dosage_invmask;
  const Dosage* nonmale_dosage_invmask = ctx->nonmale_dosage_invmask;
  uintptr_t* cur_nm_buf = ctx->cur_nm_bufs? ctx->cur_nm_bufs[tidx] : nullptr;
  Dosage* invmask_buf = ctx->invmask_bufs? ctx->invmask_bufs[tidx] : nullptr;
  const uint32_t* chr_fo_idx_end = ctx->chr_fo_idx_end;
  const uint32_t chr_ct = ctx->chr_ct;
  const uint32_t chrx_idx = ctx->chrx_idx;
  const uint32_t x_exists = (chrx_idx < UINT32_MAXM1);
  const uint32_t founder_ct = ctx->founder_ct;
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t is_unsquared = ctx->is_unsquared;
  const uint32_t triangle_calc = ctx->triangle_calc;
  const R2PhaseType unpack_phase_type = S_CAST(R2PhaseType, ctx->phase_type);
  const uint32_t check_dosage = ctx->check_dosage;
  const uintptr_t variant_ct = ctx->variant_ct;
  const uintptr_t unpacked_variant_byte_stride = ctx->unpacked_variant_byte_stride;
  // only flips when moving to next row-window.  detect this with
  // (cur_col_variant_idx_start == 0).
  // initialize to 1 instead of 0 so we don't need to special-case first row.
  uint32_t row_parity = 1;

  // always flips
  uint32_t col_parity = 0;
  do {
    const uint32_t cur_row_variant_idx_start = ctx->cur_row_variant_idx_start;
    const uint32_t cur_col_variant_idx_start = ctx->cur_col_variant_idx_start;
    if (cur_col_variant_idx_start == 0) {
      row_parity = 1 - row_parity;
    }
    const uint64_t row_window_size = ctx->row_window_size;
    const uint32_t row_start_offset = (row_window_size * tidx) / calc_thread_ct;
    const uint32_t row_end_offset = (row_window_size * (tidx + 1)) / calc_thread_ct;
    const uintptr_t shard_row_variant_idx_end = cur_row_variant_idx_start + row_end_offset;
    uintptr_t shard_row_variant_idx_start = cur_row_variant_idx_start + row_start_offset;
    uint64_t start_coord;
    if (triangle_calc) {
      start_coord = (cur_row_variant_idx_start * S_CAST(uint64_t, cur_row_variant_idx_start + 1)) / 2;
      if (cur_col_variant_idx_start > shard_row_variant_idx_start) {
        // shards are rectangle-shaped.  we're computing the lower-left
        // half-triangle of the matrix.
        // so shards closer to the upper-right may have truncated or even empty
        // intersection with the region we're computing.
        shard_row_variant_idx_start = cur_col_variant_idx_start;
      }
    } else {
      start_coord = cur_row_variant_idx_start * S_CAST(uint64_t, variant_ct);
    }
    if (shard_row_variant_idx_end > shard_row_variant_idx_start) {
      const unsigned char* unpacked_row_variants_iter = &(ctx->unpacked_row_variants[row_parity][row_start_offset * unpacked_variant_byte_stride]);
      const unsigned char* unpacked_col_variants = ctx->unpacked_col_variants[col_parity];
      double* results_d = ctx->results_d[row_parity];
      float* results_f = ctx->results_f[row_parity];
      double* results_d_row = nullptr;
      float* results_f_row = nullptr;
      const uint32_t col_initial_chr_idx = LowerBoundNonemptyU32(chr_fo_idx_end, chr_ct, cur_col_variant_idx_start + 1);
      uint32_t row_chr_idx = LowerBoundNonemptyU32(chr_fo_idx_end, chr_ct, shard_row_variant_idx_start + 1);
      uint32_t row_chr_end = chr_fo_idx_end[row_chr_idx];
      uint32_t row_is_chrx = (row_chr_idx == chrx_idx);
      const uint32_t shard_col_variant_idx_end = cur_col_variant_idx_start + ctx->col_window_size;
      uint32_t col_variant_idx_stop = shard_col_variant_idx_end;
      for (uint32_t row_variant_idx = shard_row_variant_idx_start; row_variant_idx != shard_row_variant_idx_end; ++row_variant_idx, unpacked_row_variants_iter = &(unpacked_row_variants_iter[unpacked_variant_byte_stride])) {
        if (row_variant_idx == row_chr_end) {
          ++row_chr_idx;
          row_chr_end = chr_fo_idx_end[row_chr_idx];
          row_is_chrx = (row_chr_idx == chrx_idx);
        }
        uint64_t row_start_coord;
        if (triangle_calc) {
          row_start_coord = (row_variant_idx * S_CAST(uint64_t, row_variant_idx + 1)) / 2;
          col_variant_idx_stop = MINV(row_variant_idx + 1, shard_col_variant_idx_end);
        } else {
          row_start_coord = row_variant_idx * S_CAST(uint64_t, variant_ct);
        }
        {
          const uintptr_t row_start_coord_offset = row_start_coord - start_coord;
          if (results_d) {
            results_d_row = &(results_d[row_start_coord_offset]);
          } else {
            results_f_row = &(results_f[row_start_coord_offset]);
          }
        }
        R2Variant row_r2v;
        FillR2V(unpacked_row_variants_iter, founder_ct, unpack_phase_type, x_exists, check_dosage, &row_r2v);
        uint32_t col_chr_idx = col_initial_chr_idx;
        uint32_t col_chr_end = chr_fo_idx_end[col_chr_idx];
        uint32_t col_is_chrx = (col_chr_idx == chrx_idx);
        uint32_t either_is_chrx = row_is_chrx || col_is_chrx;
        uint32_t same_chr = (row_chr_idx == col_chr_idx);
        R2PhaseType compare_phase_type = same_chr? unpack_phase_type : R2PhaseOmit(unpack_phase_type);
        const unsigned char* unpacked_col_variants_iter = unpacked_col_variants;
        for (uint32_t col_variant_idx = cur_col_variant_idx_start; col_variant_idx != col_variant_idx_stop; ++col_variant_idx, unpacked_col_variants_iter = &(unpacked_col_variants_iter[unpacked_variant_byte_stride])) {
          if (col_variant_idx == col_chr_end) {
            ++col_chr_idx;
            col_chr_end = chr_fo_idx_end[col_chr_idx];
            col_is_chrx = (col_chr_idx == chrx_idx);
            either_is_chrx = row_is_chrx || col_is_chrx;
            same_chr = (row_chr_idx == col_chr_idx);
            compare_phase_type = same_chr? unpack_phase_type : R2PhaseOmit(unpack_phase_type);
          }
          R2Variant col_r2v;
          FillR2V(unpacked_col_variants_iter, founder_ct, unpack_phase_type, either_is_chrx, check_dosage, &col_r2v);
          double result;
          uint32_t is_neg;
          if (!either_is_chrx) {
            result = ComputeR2(&row_r2v, &col_r2v, founder_ct, compare_phase_type, check_dosage, nullptr, nullptr, &is_neg);
          } else {
            result = ComputeXR2(&row_r2v, &col_r2v, founder_male_collapsed, founder_nonmale_collapsed, male_dosage_invmask, nonmale_dosage_invmask, founder_ct, founder_male_ct, compare_phase_type, check_dosage, same_chr, nullptr, nullptr, &is_neg, cur_nm_buf, invmask_buf);
          }
          if (result == -DBL_MAX) {
            result = 0.0 / 0.0;
          } else if (is_unsquared) {
            result = sqrt(result);
            if (is_neg) {
              result = -result;
            }
          }
          if (results_d_row) {
            results_d_row[col_variant_idx] = result;
          } else {
            results_f_row[col_variant_idx] = S_CAST(float, result);
          }
        }
      }
    }
    col_parity = 1 - col_parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

typedef struct VcorMatrixWriteCtxStruct {
  unsigned char* zbuf;
  VcorFlags flags;
  uint32_t orig_variant_ct;

  uint32_t cur_row_variant_idx_start;
  uint32_t row_window_size;

  double* results_d[2];
  float* results_f[2];

  FILE* outfile;
  CompressStreamState css;
  char* cswritep;

  PglErr reterr;
} VcorMatrixWriteCtx;

THREAD_FUNC_DECL VcorMatrixWriteThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  VcorMatrixWriteCtx* ctx = S_CAST(VcorMatrixWriteCtx*, arg->sharedp->context);
  const unsigned char* zbuf = ctx->zbuf;
  const VcorFlags flags = ctx->flags;
  const uintptr_t orig_variant_ct = ctx->orig_variant_ct;
  const uint32_t triangle_calc = ((flags & (kfVcorMatrixSq0 | kfVcorMatrixTri)) != 0);
  const uint32_t is_square0 = (flags / kfVcorMatrixSq0) & 1;
  const uint32_t is_bin = ((flags & (kfVcorBin8 | kfVcorBin4)) != 0);
  const uint32_t is_bin4 = (flags / kfVcorBin4) & 1;
  uint32_t row_parity = 0;
  do {
    double* cur_results_d_iter = ctx->results_d[row_parity];
    const uintptr_t row_idx_start = ctx->cur_row_variant_idx_start;
    const uintptr_t row_window_size = ctx->row_window_size;
    if (is_bin) {
      FILE* outfile = ctx->outfile;
      if (is_bin4) {
        float* cur_results_f_iter = ctx->results_f[row_parity];
        if (is_square0) {
          const uintptr_t row_idx_stop = row_idx_start + row_window_size;
          for (uintptr_t row_idx_p1 = row_idx_start + 1; row_idx_p1 <= row_idx_stop; ++row_idx_p1) {
            if (unlikely(fwrite_checked(cur_results_f_iter, row_idx_p1 * sizeof(float), outfile) ||
                         fwrite_checked(zbuf, (orig_variant_ct - row_idx_p1) * sizeof(float), outfile))) {
              goto VcorMatrixWriteThread_ret_WRITE_FAIL;
            }
            cur_results_f_iter = &(cur_results_f_iter[row_idx_p1]);
          }
        } else {
          uintptr_t entry_ct;
          if (triangle_calc) {
            entry_ct = ((row_idx_start * 2 + row_window_size + 1) * row_window_size) / 2;
          } else {
            entry_ct = row_window_size * orig_variant_ct;
          }
          if (unlikely(fwrite_checked(cur_results_f_iter, entry_ct * sizeof(float), outfile))) {
            goto VcorMatrixWriteThread_ret_WRITE_FAIL;
          }
        }
      } else {
        if (is_square0) {
          const uintptr_t row_idx_stop = row_idx_start + row_window_size;
          for (uintptr_t row_idx_p1 = row_idx_start + 1; row_idx_p1 <= row_idx_stop; ++row_idx_p1) {
            if (unlikely(fwrite_checked(cur_results_d_iter, row_idx_p1 * sizeof(double), outfile) ||
                         fwrite_checked(zbuf, (orig_variant_ct - row_idx_p1) * sizeof(double), outfile))) {
              goto VcorMatrixWriteThread_ret_WRITE_FAIL;
            }
            cur_results_d_iter = &(cur_results_d_iter[row_idx_p1]);
          }
        } else {
          uintptr_t entry_ct;
          if (triangle_calc) {
            entry_ct = ((row_idx_start * 2 + row_window_size + 1) * row_window_size) / 2;
          } else {
            entry_ct = row_window_size * orig_variant_ct;
          }
          if (unlikely(fwrite_checked(cur_results_d_iter, entry_ct * sizeof(double), outfile))) {
            goto VcorMatrixWriteThread_ret_WRITE_FAIL;
          }
        }
      }
    } else {
      char* cswritep = ctx->cswritep;
      CompressStreamState* cssp = &(ctx->css);
      const uintptr_t row_idx_stop = row_idx_start + row_window_size;
      uintptr_t col_idx_stop = triangle_calc? (row_idx_start + 1) : orig_variant_ct;
      for (uintptr_t row_idx = row_idx_start; row_idx != row_idx_stop; ++row_idx) {
        for (uintptr_t col_idx = 0; col_idx != col_idx_stop; ++col_idx) {
          cswritep = dtoa_g(*cur_results_d_iter++, cswritep);
          *cswritep++ = '\t';
        }
        if (triangle_calc) {
          if (is_square0) {
            AppendZerotabsUnsafe(orig_variant_ct - col_idx_stop, &cswritep);
          }
          ++col_idx_stop;
        }
        --cswritep;
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(cssp, &cswritep))) {
          goto VcorMatrixWriteThread_ret_WRITE_FAIL;
        }
      }
      ctx->cswritep = cswritep;
    }
    row_parity = 1 - row_parity;
    while (0) {
    VcorMatrixWriteThread_ret_WRITE_FAIL:
      ctx->reterr = kPglRetWriteFail;
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr VcorMatrix(const uintptr_t* orig_variant_include, const ChrInfo* cip, const char* const* variant_ids, const AlleleCode* maj_alleles, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const VcorInfo* vcip, const char* flagname, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  // todo: take a real look at BLIS
  // todo: if we can't just rework the whole function around BLIS, try using
  // dsyrk() for just the unphased, dosage-present, no-missing-genotype,
  // no-chrX-mixed-sex, no-chrY-nonmale case (see CalcGrm[Part]Thread())
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  VcorMatrixCtx ctx;
  VcorMatrixWriteCtx write_ctx;
  ThreadGroup tg;
  ThreadGroup write_tg;
  write_ctx.outfile = nullptr;
  PreinitCstream(&write_ctx.css);
  write_ctx.cswritep = nullptr;
  PreinitThreads(&tg);
  PreinitThreads(&write_tg);
  {
    // Don't want to do much shared initialization with VcorTable() upfront,
    // since chrX/chrY presence may be affected by
    // --ld-snp/--ld-snps/--ld-snp-list
    const VcorFlags flags = vcip->flags;
    if (unlikely((orig_variant_ct > 400000) && (parallel_tot == 1) && (!(flags & kfVcorYesReally)))) {
      logerrprintfww("Error: Gigantic (over 400k variants) %s unfiltered, non-distributed computation. Rerun with the 'yes-really' modifier if you are SURE you have enough hard drive space and want to do this.\n", flagname);
      goto VcorMatrix_ret_INCONSISTENT_INPUT;
    }
    const uint32_t triangle_calc = ((flags & (kfVcorMatrixSq0 | kfVcorMatrixTri)) != 0);
    const uint32_t first_variant_uidx = FindNth1BitFrom(orig_variant_include, 0, 1);
    uint32_t row_variant_uidx_start = first_variant_uidx;
    uint32_t row_variant_idx_start = 0;
    uint32_t row_variant_idx_stop = orig_variant_ct;
    uint32_t variant_uidx_stop = raw_variant_ct;
    uint32_t variant_ct = orig_variant_ct;
    const uintptr_t* variant_include = orig_variant_include;
    if (parallel_tot != 1) {
      if (unlikely(variant_ct < 2 * parallel_tot)) {
        logerrprintf("Error: Too few variants in %s run for --parallel %u %u.\n", flagname, parallel_idx + 1, parallel_tot);
        goto VcorMatrix_ret_INCONSISTENT_INPUT;
      }
      if (triangle_calc) {
        ParallelBounds(variant_ct, 0, parallel_idx, parallel_tot, R_CAST(int32_t*, &row_variant_idx_start), R_CAST(int32_t*, &row_variant_idx_stop));
        if (row_variant_idx_stop != variant_ct) {
          // Outside of square0's zero-padding, we can just act as if there are
          // fewer variants.
          // (though we need to be a bit careful with chrX/chrY check)
          variant_ct = row_variant_idx_stop;
          variant_uidx_stop = 1 + FindNth1BitFrom(variant_include, 0, variant_ct);
          const uint32_t variant_uidx_stopl = BitCtToWordCt(variant_uidx_stop);
          uintptr_t* new_variant_include;
          if (unlikely(bigstack_alloc_w(variant_uidx_stopl, &new_variant_include))) {
            goto VcorMatrix_ret_INCONSISTENT_INPUT;
          }
          memcpy(new_variant_include, variant_include, variant_uidx_stopl * sizeof(intptr_t));
          ZeroTrailingBits(variant_uidx_stop, new_variant_include);
          variant_include = new_variant_include;
        }
      } else {
        row_variant_idx_start = (S_CAST(uint64_t, variant_ct) * parallel_idx) / parallel_tot;
        row_variant_idx_stop = (S_CAST(uint64_t, variant_ct) * (parallel_idx + 1)) / parallel_tot;
      }
      if (row_variant_idx_start) {
        row_variant_uidx_start = FindNth1BitFrom(variant_include, first_variant_uidx, 1 + row_variant_idx_start);
      }
    }

    const uint32_t phased_calc = (flags / kfVcorPhased) & 1;
    const uint32_t is_unsquared = (flags / kfVcorUnsquared) & 1;
    const uint32_t is_bin = ((flags & (kfVcorBin8 | kfVcorBin4)) != 0);
    const uint32_t is_bin4 = (flags / kfVcorBin4) & 1;
    const uint32_t matrix_output_zst = (flags / kfVcorZs) & 1;
    {
      write_ctx.zbuf = nullptr;
      if ((flags & kfVcorMatrixSq0) && is_bin) {
        const uintptr_t max_zeroes_needed = orig_variant_ct + 1 - row_variant_idx_start;
        const uintptr_t max_zerobytes_needed = (2 - is_bin4) * 4 * max_zeroes_needed;
        if (bigstack_alloc_uc(max_zerobytes_needed, &write_ctx.zbuf)) {
          goto VcorMatrix_ret_NOMEM;
        }
        memset(write_ctx.zbuf, 0, max_zerobytes_needed);
      }

      // Write .vars file (unless parallel_idx > 0), and initialize main
      // filename
      char* outname_write_iter = outname_end;
      *outname_write_iter++ = '.';
      if (!phased_calc) {
        outname_write_iter = strcpya_k(outname_write_iter, "un");
      }
      outname_write_iter = strcpya_k(outname_write_iter, "phased.vcor");
      *outname_write_iter++ = '2' - is_unsquared;
      if (is_bin) {
        outname_write_iter = strcpya_k(outname_write_iter, ".bin");
      }

      if (parallel_idx == 0) {
        strcpy_k(outname_write_iter, ".vars");
        // this is never .zst-compressed since, if it's large enough for that
        // to matter, the matrix itself is ridiculously large
        if (fopen_checked(outname, FOPEN_WB, &write_ctx.outfile)) {
          goto VcorMatrix_ret_OPEN_FAIL;
        }
        // variant IDs limited to 16k chars
        char* write_iter = g_textbuf;
        char* write_flush = &(g_textbuf[kMaxMediumLine]);
        uintptr_t variant_uidx_base = 0;
        uintptr_t cur_bits = variant_include[0];
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          write_iter = strcpya(write_iter, variant_ids[variant_uidx]);
          AppendBinaryEoln(&write_iter);
          if (unlikely(fwrite_ck(write_flush, write_ctx.outfile, &write_iter))) {
            goto VcorMatrix_ret_WRITE_FAIL;
          }
        }
        if (fclose_flush_null(write_flush, write_iter, &write_ctx.outfile)) {
          goto VcorMatrix_ret_WRITE_FAIL;
        }
        logprintfww("%s: Variant IDs written to %s .\n", flagname, outname);
      }
      if (parallel_tot != 1) {
        *outname_write_iter++ = '.';
        outname_write_iter = u32toa(parallel_idx + 1, outname_write_iter);
      }
      if (matrix_output_zst) {
        outname_write_iter = strcpya_k(outname_write_iter, ".zst");
      }
      *outname_write_iter = '\0';
    }

    if (is_bin) {
      if (fopen_checked(outname, FOPEN_WB, &write_ctx.outfile)) {
        goto VcorMatrix_ret_OPEN_FAIL;
      }
    } else {
      const uintptr_t overflow_buf_size = kCompressStreamBlock + (kMaxDoubleGSlen + 1) * variant_ct + strlen(EOLN_STR) - 1;
      uint32_t compress_thread_ct = 1;
      if ((!phased_calc) && (max_thread_ct > 4) && (founder_ct <= 8192)) {
        compress_thread_ct = 2 + ((max_thread_ct > 8) && (founder_ct <= 4096));
      }
      reterr = InitCstreamAlloc(outname, 0, matrix_output_zst, compress_thread_ct, overflow_buf_size, &write_ctx.css, &write_ctx.cswritep);
      if (unlikely(reterr)) {
        goto VcorMatrix_ret_1;
      }
    }

    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
    const uint32_t founder_ctv = BitCtToVecCt(founder_ct);
    const uint32_t founder_ctv2 = NypCtToVecCt(founder_ct);
    const uint32_t founder_ctaw = founder_ctv * kWordsPerVec;
    const uint32_t founder_male_ct = PopcountWordsIntersect(founder_info, sex_male, raw_sample_ctl);
    const uint32_t all_haploid = IsSet(cip->haploid_mask, 0);
    PgenGlobalFlags effective_gflags = PgrGetGflags(simple_pgrp) & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
    if (variant_ct < variant_uidx_stop) {
      effective_gflags &= GflagsVfilter(variant_include, PgrGetVrtypes(simple_pgrp), variant_uidx_stop, effective_gflags);
    }
    const uint32_t check_phase = phased_calc && (!all_haploid) && (effective_gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent));
    if (!check_phase) {
      effective_gflags &= kfPgenGlobalDosagePresent;
    }
    const R2PhaseType phase_type = GetR2PhaseType(phased_calc, check_phase);
    const uint32_t check_dosage = (effective_gflags / kfPgenGlobalDosagePresent) & 1;

    uintptr_t* raregeno = nullptr;
    uint32_t* difflist_sample_ids = nullptr;
    const uint32_t max_difflist_len = founder_ct / 64;
    if (!phased_calc) {
      if (unlikely(bigstack_alloc_w(NypCtToWordCt(max_difflist_len), &raregeno) ||
                   bigstack_alloc_u32(max_difflist_len, &difflist_sample_ids))) {
        goto VcorMatrix_ret_NOMEM;
      }
    }

    // create abbreviated chromosome-view for use by worker threads.
    // they may need to know where each chromosome starts/ends (use phase or
    // not?), and which chromosome is X; that's it
    uint32_t x_code = UINT32_MAX;
    uint32_t x_fo_idx = UINT32_MAX;
    // if all-males, we can ignore phase on chrX, as well as skipping
    // male/nonmale-specific stats
    // if all-nonmales, we initialize x_code to prevent phase from being
    // ignored, but can skip the male/nonmale-specific stats
    if (founder_male_ct != founder_ct) {
      if (XymtExists(cip, kChrOffsetX, &x_code) && (founder_male_ct != 0)) {
        x_fo_idx = cip->chr_idx_to_foidx[x_code];
      }
    }
    const uintptr_t* founder_male_collapsed = nullptr;
    const Dosage* male_dosage_invmask = nullptr;
    ctx.founder_nonmale_collapsed = nullptr;
    ctx.nonmale_dosage_invmask = nullptr;
    {
      const uint32_t chr_ct = cip->chr_ct;
      uint32_t* chr_fo_idx_end;
      if (bigstack_alloc_u32(chr_ct, &chr_fo_idx_end)) {
        goto VcorMatrix_ret_NOMEM;
      }
      // not quite the same as FillSubsetChrFoVidxStart since we prune
      // now-empty chromosomes.  may want to have this logic in plink2_common
      // too
      uint32_t new_chr_ct = 0;
      uint32_t chr_variant_uidx_start = 0;
      uint32_t chr_variant_idx_start = 0;
      ctx.chrx_idx = UINT32_MAX;
      for (uint32_t chr_fo_idx = 0; chr_variant_idx_start < variant_ct; ++chr_fo_idx) {
        uint32_t chr_variant_uidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        if (chr_variant_uidx_end > variant_uidx_stop) {
          chr_variant_uidx_end = variant_uidx_stop;
        }
        const uint32_t chr_variant_ct = PopcountBitRange(variant_include, chr_variant_uidx_start, chr_variant_uidx_end);
        if (chr_variant_ct) {
          if (chr_fo_idx == x_fo_idx) {
            ctx.chrx_idx = new_chr_ct;
          }
          const uint32_t chr_variant_idx_end = chr_variant_idx_start + chr_variant_ct;
          chr_fo_idx_end[new_chr_ct++] = chr_variant_idx_end;
          chr_variant_idx_start = chr_variant_idx_end;
        }
        chr_variant_uidx_start = chr_variant_uidx_end;
      }
      assert(chr_variant_idx_start == variant_ct);
      ctx.chr_ct = new_chr_ct;
      BigstackShrinkTop(chr_fo_idx_end, new_chr_ct * sizeof(int32_t));
      ctx.chr_fo_idx_end = chr_fo_idx_end;
      if (ctx.chrx_idx != UINT32_MAX) {
        uintptr_t* founder_male_collapsed_fill;
        uintptr_t* founder_nonmale_collapsed;
        if (unlikely(bigstack_alloc_w(founder_ctl, &founder_male_collapsed_fill) ||
                     bigstack_alloc_w(founder_ctl, &founder_nonmale_collapsed))) {
          goto VcorMatrix_ret_NOMEM;
        }
        CopyBitarrSubset(sex_male, founder_info, founder_ct, founder_male_collapsed_fill);
        founder_male_collapsed = founder_male_collapsed_fill;
        BitvecInvertCopy(founder_male_collapsed, founder_ctl, founder_nonmale_collapsed);
        ZeroTrailingBits(founder_ct, founder_nonmale_collapsed);
        ctx.founder_nonmale_collapsed = founder_nonmale_collapsed;
        if (check_dosage) {
          const uint32_t founder_ctad = RoundUpPow2(founder_ct, kDosagePerVec);
          Dosage* male_dosage_invmask_fill;
          Dosage* nonmale_dosage_invmask;
          if (bigstack_alloc_dosage(founder_ctad, &male_dosage_invmask_fill) ||
              bigstack_alloc_dosage(founder_ctad, &nonmale_dosage_invmask)) {
            goto VcorMatrix_ret_NOMEM;
          }
          Expand1bitTo16(founder_male_collapsed, founder_ctad, 0xffff, male_dosage_invmask_fill);
          Expand1bitTo16(founder_nonmale_collapsed, founder_ctad, 0xffff, nonmale_dosage_invmask);
          male_dosage_invmask = male_dosage_invmask_fill;
          ctx.nonmale_dosage_invmask = nonmale_dosage_invmask;
        }
      }
    }
    const uint32_t x_exists = (ctx.chrx_idx < UINT32_MAXM1);
    ctx.founder_male_collapsed = founder_male_collapsed;
    ctx.male_dosage_invmask = male_dosage_invmask;

    uint32_t y_start;
    uint32_t y_end;
    GetXymtStartAndEnd(cip, kChrOffsetY, &y_start, &y_end);
    uintptr_t* founder_female_collapsed = nullptr;
    uintptr_t* founder_female_collapsed_interleaved = nullptr;
    if (y_end) {
      if (y_end > variant_uidx_stop) {
        y_end = variant_uidx_stop;
      }
      if ((founder_male_ct == founder_ct) || (y_start >= variant_uidx_stop) || AllBitsAreZero(variant_include, y_start, y_end)) {
        y_start = 0;
        y_end = 0;
      }
      if (y_end) {
        uintptr_t* founder_female;
        if (bigstack_end_alloc_w(raw_sample_ctl, &founder_female)) {
          goto VcorMatrix_ret_NOMEM;
        }
        BitvecInvmaskCopy(sex_nm, sex_male, raw_sample_ctl, founder_female);
        BitvecAnd(founder_info, raw_sample_ctl, founder_female);
        if (AllWordsAreZero(founder_female, raw_sample_ctl)) {
          y_start = 0;
          y_end = 0;
        } else {
          if (bigstack_alloc_w(founder_ctaw, &founder_female_collapsed) ||
              bigstack_alloc_w(founder_ctaw, &founder_female_collapsed_interleaved)) {
            goto VcorMatrix_ret_NOMEM;
          }
          CopyBitarrSubset(founder_female, founder_info, founder_ct, founder_female_collapsed);
          ZeroTrailingWords(founder_ctl, founder_female_collapsed);
          FillInterleavedMaskVec(founder_female_collapsed, founder_ctv, founder_female_collapsed_interleaved);
        }
        BigstackEndReset(bigstack_end_mark);
      }
    }
    ctx.founder_ct = founder_ct;
    ctx.founder_male_ct = founder_male_ct;
    ctx.is_unsquared = is_unsquared;
    ctx.triangle_calc = triangle_calc;

    const uintptr_t bitvec_byte_ct = BitCtToVecCt(founder_ct) * kBytesPerVec;
    uintptr_t dosagevec_byte_ct = 0;
    uintptr_t unpacked_variant_byte_stride;
    if (check_dosage) {
      dosagevec_byte_ct = DivUp(founder_ct, kDosagePerVec) * kBytesPerVec;
      const uintptr_t dosage_trail_byte_ct = LdDosageTrailAlignedByteCt(S_CAST(R2PhaseType, phased_calc), x_exists);
      unpacked_variant_byte_stride = dosagevec_byte_ct * (1 + phased_calc + check_phase) + bitvec_byte_ct + dosage_trail_byte_ct;
    } else {
      unpacked_variant_byte_stride = RoundUpPow2(16, kBytesPerVec) + bitvec_byte_ct * (3 + 2 * check_phase);
#ifndef USE_AVX2
      const uintptr_t sparse_req = RoundUpPow2((6 + max_difflist_len) * sizeof(int32_t), kBytesPerVec) + NypCtToVecCt(max_difflist_len) * kBytesPerVec + RoundUpPow2(founder_ctl * (kBytesPerWord + sizeof(int32_t)), kBytesPerVec);
      if (sparse_req > unpacked_variant_byte_stride) {
        unpacked_variant_byte_stride = sparse_req;
      }
#endif
      const uintptr_t nondosage_trail_byte_ct = LdNondosageTrailAlignedByteCt(S_CAST(R2PhaseType, phased_calc), x_exists);
      unpacked_variant_byte_stride += nondosage_trail_byte_ct;
    }
    uint32_t* founder_info_cumulative_popcounts;
    PgenVariant pgv;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &founder_info_cumulative_popcounts) ||
                 BigstackAllocPgv(founder_ct, 0, effective_gflags, &pgv))) {
      goto VcorMatrix_ret_NOMEM;
    }
    FillCumulativePopcounts(founder_info, raw_sample_ctl, founder_info_cumulative_popcounts);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(founder_info_cumulative_popcounts, simple_pgrp, &pssi);

    ctx.phase_type = phase_type;
    ctx.check_dosage = check_dosage;
    ctx.variant_ct = variant_ct;
    ctx.unpacked_variant_byte_stride = unpacked_variant_byte_stride;

    // Determine row-window size.  Byte cost of each row-window variant:
    //   unpacked_variant_byte_stride for preprocessed variant data
    //   variant_ct * 4 * (2 - is_bin4) for results (this is an overestimate
    //     for triangle case)
    // Assign up to ~half of remaining memory to this (ok to overshoot
    // slightly).
    uint32_t usual_row_window_size = row_variant_idx_stop - row_variant_idx_start;
    uint32_t calc_thread_ct = max_thread_ct - (max_thread_ct > 4) - (max_thread_ct > 8);
    if (calc_thread_ct > usual_row_window_size) {
      calc_thread_ct = usual_row_window_size;
    }
    ctx.cur_nm_bufs = nullptr;
    ctx.invmask_bufs = nullptr;
    if (ctx.chrx_idx != UINT32_MAX) {
      if (unlikely(bigstack_alloc_wp(calc_thread_ct, &ctx.cur_nm_bufs))) {
        goto VcorMatrix_ret_NOMEM;
      }
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        if (unlikely(bigstack_alloc_w(founder_ctl, &(ctx.cur_nm_bufs[tidx])))) {
          goto VcorMatrix_ret_NOMEM;
        }
      }
      if (check_dosage) {
        if (unlikely(bigstack_alloc_dosagep(calc_thread_ct, &ctx.invmask_bufs))) {
          goto VcorMatrix_ret_NOMEM;
        }
        for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
          // allocation automatically rounded up to at least end of vector
          if (unlikely(bigstack_alloc_dosage(founder_ct, &(ctx.invmask_bufs[tidx])))) {
            goto VcorMatrix_ret_NOMEM;
          }
        }
      }
    }
    {
      const uintptr_t row_byte_cost = unpacked_variant_byte_stride + variant_ct * (4 * k1LU) * (2 - is_bin4);
      const uintptr_t row_capacity = bigstack_left() / (4 * row_byte_cost);
      if (row_capacity < usual_row_window_size) {
        // If multipass, may as well make usual_row_window_size a multiple of
        // thread_ct.
        if (row_capacity < calc_thread_ct) {
          calc_thread_ct = row_capacity;
        }
        const uint32_t thread_workload = row_capacity / calc_thread_ct;
        if (unlikely(!thread_workload)) {
          goto VcorMatrix_ret_NOMEM;
        }
        usual_row_window_size = thread_workload * calc_thread_ct;
      }
    }
    ctx.cur_row_variant_idx_start = row_variant_idx_start;
    ctx.row_window_size = usual_row_window_size;
    if (unlikely(bigstack_alloc_uc(usual_row_window_size * unpacked_variant_byte_stride, &(ctx.unpacked_row_variants[0])) ||
                 bigstack_alloc_uc(usual_row_window_size * unpacked_variant_byte_stride, &(ctx.unpacked_row_variants[1])))) {
      goto VcorMatrix_ret_NOMEM;
    }
    {
      const uintptr_t slot_ct = usual_row_window_size * variant_ct;
      if (is_bin4) {
        if (unlikely(bigstack_alloc_f(slot_ct, &(ctx.results_f[0])) ||
                     bigstack_alloc_f(slot_ct, &(ctx.results_f[1])))) {
          goto VcorMatrix_ret_NOMEM;
        }
        ctx.results_d[0] = nullptr;
        ctx.results_d[1] = nullptr;
      } else {
        if (unlikely(bigstack_alloc_d(slot_ct, &(ctx.results_d[0])) ||
                     bigstack_alloc_d(slot_ct, &(ctx.results_d[1])))) {
          goto VcorMatrix_ret_NOMEM;
        }
        ctx.results_f[0] = nullptr;
        ctx.results_f[1] = nullptr;
      }
    }
    write_ctx.flags = flags;
    write_ctx.orig_variant_ct = orig_variant_ct;
    write_ctx.row_window_size = usual_row_window_size;
    write_ctx.results_d[0] = ctx.results_d[0];
    write_ctx.results_d[1] = ctx.results_d[1];
    write_ctx.results_f[0] = ctx.results_f[0];
    write_ctx.results_f[1] = ctx.results_f[1];
    write_ctx.reterr = kPglRetSuccess;
    if (unlikely(SetThreadCt(1, &write_tg))) {
      goto VcorMatrix_ret_NOMEM;
    }
    SetThreadFuncAndData(VcorMatrixWriteThread, &write_ctx, &write_tg);

    uint32_t usual_col_window_size = variant_ct;
    {
      const uintptr_t half_bytes_avail = RoundDownPow2(bigstack_left() / 2, kCacheline);
      if (unlikely(half_bytes_avail < unpacked_variant_byte_stride)) {
        goto VcorMatrix_ret_NOMEM;
      }
      if (variant_ct * S_CAST(uint64_t, unpacked_variant_byte_stride) > half_bytes_avail) {
        usual_col_window_size = half_bytes_avail / unpacked_variant_byte_stride;
      }
    }
    if (unlikely(bigstack_alloc_uc(usual_col_window_size * unpacked_variant_byte_stride, &(ctx.unpacked_col_variants[0])) ||
                 bigstack_alloc_uc(usual_col_window_size * unpacked_variant_byte_stride, &(ctx.unpacked_col_variants[1])))) {
      goto VcorMatrix_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto VcorMatrix_ret_NOMEM;
    }
    SetThreadFuncAndData(VcorMatrixThread, &ctx, &tg);
    uint64_t job_size;
    if (triangle_calc) {
      job_size = ((S_CAST(uint64_t, row_variant_idx_stop) * (row_variant_idx_stop + 1)) - (S_CAST(uint64_t, row_variant_idx_start) * (row_variant_idx_start + 1))) / 2;
    } else {
      job_size = (row_variant_idx_stop - row_variant_idx_start) * S_CAST(uint64_t, variant_ct);
    }
    if (flags & kfVcorRefBased) {
      maj_alleles = nullptr;
    }
    AlleleCode aidx = 0;
    uint64_t job_done = 0;
    uint64_t next_job_done = 0;
    uint64_t next_print_job_idx = (job_size + 99) / 100;
    uint32_t pct = 0;
    uint32_t cur_row_variant_idx_start = row_variant_idx_start;
    uint32_t row_window_size = usual_row_window_size;
    uint32_t col_variant_idx_stop = variant_ct;
    uint32_t row_chr_fo_idx = UINT32_MAX;  // deliberate overflow
    uint32_t row_chr_end = 0;
    uint32_t row_read_phase = 0;
    uint32_t row_parity = 0;
    uint32_t col_parity = 0;
    uintptr_t row_variant_uidx_base;
    uintptr_t row_cur_bits;
    BitIter1Start(variant_include, row_variant_uidx_start, &row_variant_uidx_base, &row_cur_bits);
    printf("%s: 0%%", flagname);
    fflush(stdout);
    do {
      // 1. unpack all variants in current row block.
      // 2. iterate through column blocks.
      if (cur_row_variant_idx_start + row_window_size > row_variant_idx_stop) {
        row_window_size = row_variant_idx_stop - cur_row_variant_idx_start;
      }
      unsigned char* row_load_iter = ctx.unpacked_row_variants[row_parity];
      for (uint32_t vidx_offset = 0; vidx_offset != row_window_size; ++vidx_offset, row_load_iter = &(row_load_iter[unpacked_variant_byte_stride])) {
        const uint32_t variant_uidx = BitIter1(variant_include, &row_variant_uidx_base, &row_cur_bits);
        if (variant_uidx >= row_chr_end) {
          do {
            ++row_chr_fo_idx;
            row_chr_end = cip->chr_fo_vidx_start[row_chr_fo_idx + 1];
          } while (variant_uidx >= row_chr_end);
          if (phase_type == kR2PhaseTypePresent) {
            const uint32_t row_chr_idx = cip->chr_file_order[row_chr_fo_idx];
            row_read_phase = (!IsSet(cip->haploid_mask, row_chr_idx)) || (row_chr_idx == x_code);
            if (!row_read_phase) {
              pgv.phasepresent_ct = 0;
              pgv.dphase_ct = 0;
            }
          }
        }
        if (maj_alleles) {
          aidx = maj_alleles[variant_uidx];
        }
        const uint32_t is_y = (variant_uidx < y_end) && (variant_uidx >= y_start);
        if (check_dosage) {
          if (row_read_phase) {
            reterr = PgrGetInv1Dp(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, &pgv);
          } else {
            reterr = PgrGetInv1D(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
          }
          if (unlikely(reterr)) {
            goto VcorMatrix_ret_PGR_FAIL;
          }
          if (is_y) {
            InterleavedSetMissingCleardosage(founder_female_collapsed, founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec, &pgv.dosage_ct, pgv.dosage_present, pgv.dosage_main);
          }
          LdUnpackDosage(&pgv, founder_male_collapsed, male_dosage_invmask, founder_ct, phase_type, row_load_iter);
        } else {
          if ((!phased_calc) && (!is_y)) {
            uint32_t difflist_common_geno;
            uint32_t difflist_len;
            reterr = PgrGetInv1DifflistOrGenovec(founder_info, pssi, founder_ct, max_difflist_len, variant_uidx, aidx, simple_pgrp, pgv.genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
            if (unlikely(reterr)) {
              goto VcorMatrix_ret_PGR_FAIL;
            }
            if (difflist_common_geno != UINT32_MAX) {
              if (difflist_len <= max_difflist_len) {
                LdUnpackNondosageSparse(raregeno, difflist_sample_ids, founder_male_collapsed, founder_ct, founder_male_ct, difflist_common_geno, difflist_len, row_load_iter);
                continue;
              }
              PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, founder_ct, difflist_len, pgv.genovec);
            }
          } else {
            if (row_read_phase) {
              reterr = PgrGetInv1P(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
            } else {
              reterr = PgrGetInv1(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec);
            }
            if (unlikely(reterr)) {
              goto VcorMatrix_ret_PGR_FAIL;
            }
            if (is_y) {
              InterleavedSetMissing(founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec);
            }
          }
          LdUnpackNondosageDense(&pgv, founder_male_collapsed, founder_ct, phase_type, row_load_iter);
        }
      }
      const uint32_t cur_row_variant_idx_stop = cur_row_variant_idx_start + row_window_size;
      if (triangle_calc) {
        col_variant_idx_stop = cur_row_variant_idx_stop;
      }
      uint32_t cur_col_variant_idx_start = 0;
      uint32_t col_window_size = usual_col_window_size;
      uint32_t col_chr_fo_idx = UINT32_MAX;  // deliberate overflow
      uint32_t col_chr_end = 0;
      uint32_t col_read_phase = 0;
      uintptr_t col_variant_uidx_base;
      uintptr_t col_cur_bits;
      BitIter1Start(variant_include, first_variant_uidx, &col_variant_uidx_base, &col_cur_bits);
      do {
        if (cur_col_variant_idx_start + col_window_size > col_variant_idx_stop) {
          col_window_size = col_variant_idx_stop - cur_col_variant_idx_start;
        }
        unsigned char* col_load_iter = ctx.unpacked_col_variants[col_parity];
        // possible todo: don't duplicate already-unpacked row variants
        for (uint32_t vidx_offset = 0; vidx_offset != col_window_size; ++vidx_offset, col_load_iter = &(col_load_iter[unpacked_variant_byte_stride])) {
          const uint32_t variant_uidx = BitIter1(variant_include, &col_variant_uidx_base, &col_cur_bits);
          if (variant_uidx >= col_chr_end) {
            do {
              ++col_chr_fo_idx;
              col_chr_end = cip->chr_fo_vidx_start[col_chr_fo_idx + 1];
            } while (variant_uidx >= col_chr_end);
            if (phase_type == kR2PhaseTypePresent) {
              const uint32_t col_chr_idx = cip->chr_file_order[col_chr_fo_idx];
              col_read_phase = (!IsSet(cip->haploid_mask, col_chr_idx)) || (col_chr_idx == x_code);
              if (!col_read_phase) {
                pgv.phasepresent_ct = 0;
                pgv.dphase_ct = 0;
              }
            }
          }
          if (maj_alleles) {
            aidx = maj_alleles[variant_uidx];
          }
          const uint32_t is_y = (variant_uidx < y_end) && (variant_uidx >= y_start);
          if (check_dosage) {
            if (col_read_phase) {
              reterr = PgrGetInv1Dp(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, &pgv);
            } else {
              reterr = PgrGetInv1D(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
            }
            if (unlikely(reterr)) {
              goto VcorMatrix_ret_PGR_FAIL;
            }
            if (is_y) {
              InterleavedSetMissingCleardosage(founder_female_collapsed, founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec, &pgv.dosage_ct, pgv.dosage_present, pgv.dosage_main);
            }
            LdUnpackDosage(&pgv, founder_male_collapsed, male_dosage_invmask, founder_ct, phase_type, col_load_iter);
          } else {
            if ((!phased_calc) && (!is_y)) {
              uint32_t difflist_common_geno;
              uint32_t difflist_len;
              reterr = PgrGetInv1DifflistOrGenovec(founder_info, pssi, founder_ct, max_difflist_len, variant_uidx, aidx, simple_pgrp, pgv.genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
              if (unlikely(reterr)) {
                goto VcorMatrix_ret_PGR_FAIL;
              }
              if (difflist_common_geno != UINT32_MAX) {
                if (difflist_len <= max_difflist_len) {
                  LdUnpackNondosageSparse(raregeno, difflist_sample_ids, founder_male_collapsed, founder_ct, founder_male_ct, difflist_common_geno, difflist_len, col_load_iter);
                  continue;
                }
                PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, founder_ct, difflist_len, pgv.genovec);
              }
            } else {
              if (col_read_phase) {
                reterr = PgrGetInv1P(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
              } else {
                reterr = PgrGetInv1(founder_info, pssi, founder_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec);
              }
              if (unlikely(reterr)) {
                goto VcorMatrix_ret_PGR_FAIL;
              }
              if (is_y) {
                InterleavedSetMissing(founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec);
              }
            }
            LdUnpackNondosageDense(&pgv, founder_male_collapsed, founder_ct, phase_type, col_load_iter);
          }
        }
        const uint32_t cur_col_variant_idx_stop = cur_col_variant_idx_start + col_window_size;
        job_done = next_job_done;
        if (triangle_calc) {
          // Shard shape can be generalized as a (possibly-zero-height)
          // trapezoid stacked on a (possibly-zero-height) rectangle.
          //
          // Trapezoid starts on the first nonempty row, and ends either after
          // the first full row or, if that's not in the shard, at the shard
          // bottom.  Rectangle starts where the trapezoid ends.

          // both of these values are guaranteed to be less than
          // cur_row_variant_idx_stop == col_variant_idx_stop
          const uintptr_t trapezoid_start_row_idx = MAXV(cur_row_variant_idx_start, cur_col_variant_idx_start);

          const uintptr_t trapezoid_stop_row_idx = MINV(cur_row_variant_idx_stop, MAXV(cur_row_variant_idx_start, cur_col_variant_idx_stop));
          const uintptr_t trapezoid_height = trapezoid_stop_row_idx - trapezoid_start_row_idx;
          if (trapezoid_height) {
            const uintptr_t first_row_len = trapezoid_start_row_idx + 1 - cur_col_variant_idx_start;
            const uintptr_t last_row_len = trapezoid_stop_row_idx - cur_col_variant_idx_start;
            next_job_done += (S_CAST(uint64_t, trapezoid_height) * (first_row_len + last_row_len)) / 2;
          }
          next_job_done += (cur_row_variant_idx_stop - trapezoid_stop_row_idx) * S_CAST(uint64_t, col_window_size);
        } else {
          next_job_done += row_window_size * S_CAST(uint64_t, col_window_size);
        }
        if ((cur_row_variant_idx_start > row_variant_idx_start) || cur_col_variant_idx_start) {
          JoinThreads(&tg);
          if (!cur_col_variant_idx_start) {
            const uint32_t prev_row_variant_idx_start = cur_row_variant_idx_start - usual_row_window_size;
            if (prev_row_variant_idx_start > row_variant_idx_start) {
              JoinThreads(&write_tg);
              if (unlikely(write_ctx.reterr)) {
                reterr = write_ctx.reterr;
                goto VcorMatrix_ret_1;
              }
            }
            write_ctx.cur_row_variant_idx_start = prev_row_variant_idx_start;
            if (unlikely(SpawnThreads(&write_tg))) {
              goto VcorMatrix_ret_THREAD_CREATE_FAIL;
            }
            ctx.cur_row_variant_idx_start = cur_row_variant_idx_start;
            ctx.row_window_size = row_window_size;
          }
        }
        ctx.cur_col_variant_idx_start = cur_col_variant_idx_start;
        ctx.col_window_size = col_window_size;
        if (next_job_done == job_size) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto VcorMatrix_ret_THREAD_CREATE_FAIL;
        }
        col_parity = 1 - col_parity;

        if (job_done >= next_print_job_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (job_done * 100) / job_size;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_job_idx = (pct * job_size + 99) / 100;
        }
        cur_col_variant_idx_start = cur_col_variant_idx_stop;
      } while (cur_col_variant_idx_start != col_variant_idx_stop);
      cur_row_variant_idx_start = cur_row_variant_idx_stop;
      row_parity = 1 - row_parity;
    } while (cur_row_variant_idx_start != row_variant_idx_stop);
    JoinThreads(&tg);
    const uint32_t prev_row_variant_idx_start = cur_row_variant_idx_start - row_window_size;
    if (prev_row_variant_idx_start > row_variant_idx_start) {
      JoinThreads(&write_tg);
      if (unlikely(write_ctx.reterr)) {
        reterr = write_ctx.reterr;
        goto VcorMatrix_ret_1;
      }
    }
    DeclareLastThreadBlock(&write_tg);
    write_ctx.cur_row_variant_idx_start = prev_row_variant_idx_start;
    write_ctx.row_window_size = row_window_size;
    if (unlikely(SpawnThreads(&write_tg))) {
      goto VcorMatrix_ret_THREAD_CREATE_FAIL;
    }
    JoinThreads(&write_tg);
    if (unlikely(write_ctx.reterr)) {
      reterr = write_ctx.reterr;
      goto VcorMatrix_ret_1;
    }
    fputs("\r", stdout);
    logprintfww("%s: Matrix%s written to %s .\n", flagname, (parallel_tot == 1)? "" : " piece", outname);
  }
  while (0) {
  VcorMatrix_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  VcorMatrix_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  VcorMatrix_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  VcorMatrix_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  VcorMatrix_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  VcorMatrix_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  }
 VcorMatrix_ret_1:
  CleanupThreads(&write_tg);
  CleanupThreads(&tg);
  CswriteCloseCond(&write_ctx.css, write_ctx.cswritep);
  fclose_cond(write_ctx.outfile);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

// indexes here are all subsetted (founder_idx / variant_idx), not sample_uidx
// / variant_uidx.
typedef struct VcorTableCtxStruct {
  // Shared constants.
  const uintptr_t* founder_male_collapsed;
  const uintptr_t* founder_nonmale_collapsed;
  const Dosage* male_dosage_invmask;
  const Dosage* nonmale_dosage_invmask;
  uint32_t chrx_idx; // UINT32_MAX if not present
  uint32_t founder_ct;
  uint32_t founder_male_ct;
  unsigned char is_unsquared;
  unsigned char phase_type;
  unsigned char check_dosage;
  unsigned char report_d;
  unsigned char report_dprime;
  uintptr_t unpacked_variant_byte_stride;

  // Input data.
  unsigned char* unpacked_variants[2];  // read from [col_parity]
  // uint32_t* row_uvidxs[2];  // read from [row_parity]
  uint32_t* col_uvidxs[2];  // read from [col_parity]
  ChrIdx* row_chr_idxs[2];  // read from [row_parity]
  ChrIdx* col_chr_idxs[2];  // read from [col_parity]
  uint32_t* col_offset_starts[2];  // read from [row_parity]
  // read from [row_parity].
  // col_offset_end = write_idx_starts[k+1] - write_idx_starts[k].
  // Note that these intervals usually include col_variant_idx ==
  // row_variant_idx.  Those entries are skipped by the writer and aren't
  // filled.
  // In the --ld-snp/--ld-snps/--ld-snp-list case, if col_variant_idx is in the
  // row-variant set, and col_variant_idx < row_variant_idx, the writer skips
  // the entry (since it was previously written with row/column swapped), but
  // we currently don't optimize out the computation.
  uintptr_t* write_idx_starts[2];

  uint32_t* row_window_sizes;  // read from [row_parity]
  uint32_t col_window_starts[2];  // read from [col_parity]
  uint32_t col_window_ends[2];  // read from [col_parity]

  // per-thread chrX workspaces.
  uintptr_t** cur_nm_bufs;
  Dosage** invmask_bufs;

  // Output double-buffer.  Write to [row_parity].  For better locality,
  // r^2 / d / dprime are adjacent (and in that order) if we're writing more
  // than one of them.
  double* results[2];
} VcorTableCtx;

THREAD_FUNC_DECL VcorTableThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  VcorTableCtx* ctx = S_CAST(VcorTableCtx*, arg->sharedp->context);
  const uintptr_t* founder_male_collapsed = ctx->founder_male_collapsed;
  const uintptr_t* founder_nonmale_collapsed = ctx->founder_nonmale_collapsed;
  const Dosage* male_dosage_invmask = ctx->male_dosage_invmask;
  const Dosage* nonmale_dosage_invmask = ctx->nonmale_dosage_invmask;
  uintptr_t* cur_nm_buf = ctx->cur_nm_bufs? ctx->cur_nm_bufs[tidx] : nullptr;
  Dosage* invmask_buf = ctx->invmask_bufs? ctx->invmask_bufs[tidx] : nullptr;
  const uint32_t chrx_idx = ctx->chrx_idx;
  const uint32_t x_exists = (chrx_idx < UINT32_MAXM1);
  const uint32_t founder_ct = ctx->founder_ct;
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t is_unsquared = ctx->is_unsquared;
  const R2PhaseType unpack_phase_type = S_CAST(R2PhaseType, ctx->phase_type);
  const uint32_t check_dosage = ctx->check_dosage;
  const uint32_t report_d = ctx->report_d;
  const uint32_t report_dprime = ctx->report_dprime;
  const uintptr_t unpacked_variant_byte_stride = ctx->unpacked_variant_byte_stride;
  const uintptr_t result_stride = 1 + report_d + report_dprime;
  // only flips when moving to next row-window.  detect this with
  // (col_window_start == 0).
  // initialize to 1 instead of 0 so we don't need to special-case first row.
  uint32_t row_parity = 1;

  // always flips
  uint32_t col_parity = 0;
  double dd = 0.0;
  double dprime = 0.0;
  double* d_ptr = report_d? (&dd) : nullptr;
  double* dprime_ptr = report_dprime? (&dprime) : nullptr;
  do {
    const uint32_t col_window_start = ctx->col_window_starts[col_parity];
    const uint32_t col_window_end = ctx->col_window_ends[col_parity];
    if (col_window_start == 0) {
      row_parity = 1 - row_parity;
    }
    const uint64_t row_window_size = ctx->row_window_sizes[row_parity];
    const uint32_t row_start_offset = (row_window_size * tidx) / calc_thread_ct;
    const uint32_t row_end_offset = (row_window_size * (tidx + 1)) / calc_thread_ct;
    if (row_end_offset > row_start_offset) {
      unsigned char* unpacked_variants = ctx->unpacked_variants[col_parity];
      // uint32_t* row_uvidxs = ctx->row_uvidxs[row_parity];
      uint32_t* col_uvidxs = ctx->col_uvidxs[col_parity];
      const ChrIdx* row_chr_idxs = ctx->row_chr_idxs[row_parity];
      const ChrIdx* col_chr_idxs = ctx->col_chr_idxs[col_parity];
      const uint32_t* col_offset_starts = ctx->col_offset_starts[row_parity];
      const uintptr_t* write_idx_starts = ctx->write_idx_starts[row_parity];
      double* results = ctx->results[row_parity];
      for (uint32_t row_offset_idx = row_start_offset; row_offset_idx != row_end_offset; ++row_offset_idx) {
        const uint32_t col_offset_start = col_offset_starts[row_offset_idx];
        const uintptr_t write_idx_start = write_idx_starts[row_offset_idx];
        uint32_t col_offset_stop = (write_idx_starts[row_offset_idx + 1] - write_idx_start) + col_offset_start;
        if ((col_offset_start >= col_window_end) || (col_offset_stop <= col_window_start)) {
          continue;
        }
        if (col_offset_stop > col_window_end) {
          col_offset_stop = col_window_end;
        }
        // const uint32_t cur_row_uvidx = row_uvidxs[row_offset_idx];
        const uint32_t cur_row_uvidx = row_offset_idx;
        const unsigned char* unpacked_row_ptr = &(unpacked_variants[row_offset_idx * unpacked_variant_byte_stride]);
        R2Variant row_r2v;
        FillR2V(unpacked_row_ptr, founder_ct, unpack_phase_type, x_exists, check_dosage, &row_r2v);
        const uint32_t row_chr_idx = row_chr_idxs[row_offset_idx];
        const uint32_t row_is_chrx = (row_chr_idx == chrx_idx);
        uint32_t col_offset_idx = MAXV(col_offset_start, col_window_start);
        double* write_iter = &(results[(write_idx_start + col_offset_idx - col_offset_start) * result_stride]);
        for (; col_offset_idx != col_offset_stop; ++col_offset_idx) {
          // bugfix (18 Apr 2024): must subtract col_window_start
          const uint32_t cur_col_uvidx = col_uvidxs[col_offset_idx - col_window_start];
          if (cur_row_uvidx == cur_col_uvidx) {
            write_iter = &(write_iter[result_stride]);
            continue;
          }
          const unsigned char* unpacked_col_ptr = &(unpacked_variants[cur_col_uvidx * unpacked_variant_byte_stride]);
          const uint32_t col_chr_idx = col_chr_idxs[col_offset_idx - col_window_start];
          const uint32_t either_is_chrx = row_is_chrx || (col_chr_idx == chrx_idx);
          R2Variant col_r2v;
          FillR2V(unpacked_col_ptr, founder_ct, unpack_phase_type, either_is_chrx, check_dosage, &col_r2v);
          const uint32_t same_chr = (row_chr_idx == col_chr_idx);
          R2PhaseType compare_phase_type = same_chr? unpack_phase_type : R2PhaseOmit(unpack_phase_type);
          double r_or_r2;
          uint32_t is_neg;
          if (!either_is_chrx) {
            r_or_r2 = ComputeR2(&row_r2v, &col_r2v, founder_ct, compare_phase_type, check_dosage, d_ptr, dprime_ptr, &is_neg);
          } else {
            r_or_r2 = ComputeXR2(&row_r2v, &col_r2v, founder_male_collapsed, founder_nonmale_collapsed, male_dosage_invmask, nonmale_dosage_invmask, founder_ct, founder_male_ct, compare_phase_type, check_dosage, same_chr, d_ptr, dprime_ptr, &is_neg, cur_nm_buf, invmask_buf);
          }
          if (r_or_r2 == -DBL_MAX) {
            r_or_r2 = 0.0 / 0.0;
          } else if (is_unsquared) {
            r_or_r2 = sqrt(r_or_r2) ;
            if (is_neg) {
              r_or_r2 = -r_or_r2;
            }
          }
          *write_iter++ = r_or_r2;
          if (report_d) {
            *write_iter++ = dd;
          }
          if (report_dprime) {
            *write_iter++ = dprime;
          }
        }
      }
    }
    col_parity = 1 - col_parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

typedef struct VcorTableWriteCtxStruct {
  const uintptr_t* variant_include;
  const uint32_t* variant_include_cumulative_popcounts;
  const uintptr_t* row_variant_include;
  const uint32_t* row_variant_include_cumulative_popcounts;
  const uintptr_t* row_subset_exclude;
  const ChrInfo* cip;
  const uint32_t* variant_bps;
  const char* const* variant_ids;
  const uintptr_t* allele_idx_offsets;
  const char* const* allele_storage;
  const uintptr_t* nonref_flags;
  const AlleleCode* maj_alleles;
  const double* allele_freqs;
  double r_or_r2_thresh;
  uint32_t raw_variant_ct;
  VcorFlags flags;
  unsigned char all_nonref;
  unsigned char provref_col;

  // these are relative to row_variant_include
  uint32_t variant_ridx_starts[2];
  uint32_t* row_window_sizes;  // [2]

  // these are relative to variant_include
  uint32_t col_variant_idx_starts[2];
  uint32_t* col_offset_starts[2];
  uintptr_t* write_idx_starts[2];

  double* results[2];

  char* row_chr_buf;
  char* col_chr_buf;

  CompressStreamState css;
  char* cswritep;

  PglErr reterr;
} VcorTableWriteCtx;

THREAD_FUNC_DECL VcorTableWriteThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  VcorTableWriteCtx* ctx = S_CAST(VcorTableWriteCtx*, arg->sharedp->context);
  const uintptr_t* variant_include = ctx->variant_include;
  const uint32_t* variant_include_cumulative_popcounts = ctx->variant_include_cumulative_popcounts;
  const uintptr_t* row_variant_include = ctx->row_variant_include;
  const uint32_t* row_variant_include_cumulative_popcounts = ctx->row_variant_include_cumulative_popcounts;
  const uintptr_t* row_subset_exclude = ctx->row_subset_exclude;
  const ChrInfo* cip = ctx->cip;
  const uint32_t* variant_bps = ctx->variant_bps;
  const char* const* variant_ids = ctx->variant_ids;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const char* const* allele_storage = ctx->allele_storage;
  const uintptr_t* nonref_flags = ctx->nonref_flags;
  const AlleleCode* maj_alleles = ctx->maj_alleles;
  const double* allele_freqs = ctx->allele_freqs;
  const VcorFlags flags = ctx->flags;
  const uint32_t inter_chr = (flags / kfVcorInterChr) & 1;
  const uint32_t chr_col = (flags / kfVcorColChrom) & 1;
  const uint32_t pos_col = (flags / kfVcorColPos) & 1;
  const uint32_t id_col = (flags / kfVcorColId) & 1;
  const uint32_t ref_col = (flags / kfVcorColRef) & 1;
  const uint32_t alt1_col = (flags / kfVcorColAlt1) & 1;
  const uint32_t alt_col = (flags / kfVcorColAlt) & 1;
  const uint32_t all_nonref = ctx->all_nonref;
  const uint32_t provref_col = ctx->provref_col;
  const uint32_t maj_col = (flags / kfVcorColMaj) & 1;
  const uint32_t nonmaj_col = (flags / kfVcorColNonmaj) & 1;
  const uint32_t freq_col = (flags / kfVcorColFreq) & 1;
  const uint32_t d_col = (flags / kfVcorColD) & 1;
  const uint32_t dprime_col = ((flags & (kfVcorColDprime | kfVcorColDprimeAbs)) != 0);
  const uint32_t dprime_abs = (flags / kfVcorColDprimeAbs) & 1;
  const uint32_t results_stride = 1 + d_col + dprime_col;
  const double r_or_r2_thresh = ctx->r_or_r2_thresh;
  const uint32_t raw_variant_ctl = BitCtToWordCt(ctx->raw_variant_ct);

  char* row_chr_buf = ctx->row_chr_buf;
  char* col_chr_buf = ctx->col_chr_buf;
  uint32_t row_chr_fo_idx = UINT32_MAX;
  uint32_t row_chr_end = 0;
  uint32_t row_chr_blen = 0;
  uint32_t col_chr_fo_idx = UINT32_MAX;
  uint32_t col_chr_end = 0;
  uint32_t col_chr_blen = 0;
  uint32_t row_allele_ct = 2;
  uint32_t col_allele_ct = 2;
  uint32_t maj_allele_idx = 0;

  uint32_t row_parity = 0;
  do {
    const uint32_t* col_offset_starts = ctx->col_offset_starts[row_parity];
    const uintptr_t* write_idx_starts = ctx->write_idx_starts[row_parity];
    const double* cur_results_iter = ctx->results[row_parity];
    const uintptr_t variant_ridx_start = ctx->variant_ridx_starts[row_parity];
    const uint32_t row_window_size = ctx->row_window_sizes[row_parity];
    const uint32_t col_variant_idx_start = ctx->col_variant_idx_starts[row_parity];
    char* cswritep = ctx->cswritep;
    CompressStreamState* cssp = &(ctx->css);
    const uint32_t row_variant_uidx_start = IdxToUidx(row_variant_include, row_variant_include_cumulative_popcounts, 0, raw_variant_ctl, variant_ridx_start);
    uintptr_t row_variant_uidx_base;
    uintptr_t row_cur_bits;
    BitIter1Start(row_variant_include, row_variant_uidx_start, &row_variant_uidx_base, &row_cur_bits);
    uint32_t col_variant_widx = 0;
    for (uint32_t row_offset = 0; row_offset != row_window_size; ++row_offset) {
      const uint32_t row_variant_uidx = BitIter1(row_variant_include, &row_variant_uidx_base, &row_cur_bits);
      if (row_variant_uidx >= row_chr_end) {
        do {
          ++row_chr_fo_idx;
          row_chr_end = cip->chr_fo_vidx_start[row_chr_fo_idx + 1];
        } while (row_variant_uidx >= row_chr_end);
        if (row_chr_buf) {
          const uint32_t row_chr_idx = cip->chr_file_order[row_chr_fo_idx];
          char* chr_name_end = chrtoa(cip, row_chr_idx, row_chr_buf);
          *chr_name_end = '\t';
          row_chr_blen = 1 + S_CAST(uintptr_t, chr_name_end - row_chr_buf);
          // row_chr_buf == col_chr_buf except in inter_chr case
          col_chr_blen = row_chr_blen;
        }
      }
      uintptr_t row_allele_idx_offset_base = row_variant_uidx * 2;
      if (allele_idx_offsets) {
        row_allele_idx_offset_base = allele_idx_offsets[row_variant_uidx];
        row_allele_ct = allele_idx_offsets[row_variant_uidx + 1] - row_allele_idx_offset_base;
      }
      const char* const* cur_row_alleles = &(allele_storage[row_allele_idx_offset_base]);
      const uint32_t col_offset_start = col_offset_starts[row_offset];
      const uint32_t col_offset_stop = write_idx_starts[row_offset + 1] - write_idx_starts[row_offset] + col_offset_start;
      if (inter_chr) {
        col_chr_fo_idx = UINT32_MAX;
        col_chr_end = 0;
      }
      const uint32_t col_variant_uidx_start = ExpsearchIdxToUidx(variant_include, variant_include_cumulative_popcounts, raw_variant_ctl, col_variant_idx_start + col_offset_start, &col_variant_widx);
      uintptr_t col_variant_uidx_base;
      uintptr_t col_cur_bits;
      BitIter1Start(variant_include, col_variant_uidx_start, &col_variant_uidx_base, &col_cur_bits);
      for (uint32_t col_offset = col_offset_start; col_offset != col_offset_stop; ++col_offset) {
        const uint32_t col_variant_uidx = BitIter1(variant_include, &col_variant_uidx_base, &col_cur_bits);
        if (col_variant_uidx == row_variant_uidx) {
          cur_results_iter = &(cur_results_iter[results_stride]);
          continue;
        }
        // In the usual case, we can report each pair exactly once by requiring
        // col_variant_uidx > row_variant_uidx.
        // However, if --ld-snp/--ld-snps/--ld-snp-list is in effect, that may
        // cause some variant-pairs to not be reported at all.  In that case,
        // we only enforce the col_variant_uidx > row_variant_uidx rule when
        // col_variant_uidx corresponds to a row-variant.  (row_subset_exclude
        // usually points to row_variant_include; the exception is when
        // --parallel is also in effect.)
        if (row_subset_exclude && (col_variant_uidx < row_variant_uidx) && IsSet(row_subset_exclude, col_variant_uidx)) {
          cur_results_iter = &(cur_results_iter[results_stride]);
          continue;
        }
        const double r_or_r2 = *cur_results_iter;
        // !(a >= b) instead of (a < b) so that NaN is handled properly
        if ((r_or_r2_thresh >= 0.0) && (!(fabs(r_or_r2) >= r_or_r2_thresh))) {
          cur_results_iter = &(cur_results_iter[results_stride]);
          continue;
        }
        ++cur_results_iter;

        if (col_variant_uidx >= col_chr_end) {
          do {
            ++col_chr_fo_idx;
            col_chr_end = cip->chr_fo_vidx_start[col_chr_fo_idx + 1];
          } while (col_variant_uidx >= col_chr_end);
          if (col_chr_buf) {
            const uint32_t col_chr_idx = cip->chr_file_order[col_chr_fo_idx];
            char* chr_name_end = chrtoa(cip, col_chr_idx, col_chr_buf);
            *chr_name_end = '\t';
            col_chr_blen = 1 + S_CAST(uintptr_t, chr_name_end - col_chr_buf);
          }
        }
        if (chr_col) {
          cswritep = memcpya(cswritep, row_chr_buf, row_chr_blen);
        }
        if (pos_col) {
          cswritep = u32toa_x(variant_bps[row_variant_uidx], '\t', cswritep);
        }
        if (id_col) {
          cswritep = strcpyax(cswritep, variant_ids[row_variant_uidx], '\t');
        }
        if (ref_col) {
          cswritep = strcpyax(cswritep, cur_row_alleles[0], '\t');
        }
        if (alt1_col) {
          cswritep = strcpyax(cswritep, cur_row_alleles[1], '\t');
        }
        if (alt_col) {
          for (uint32_t allele_idx = 1; allele_idx != row_allele_ct; ++allele_idx) {
            if (unlikely(Cswrite(cssp, &cswritep))) {
              goto VcorTableWriteThread_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, cur_row_alleles[allele_idx], ',');
          }
          cswritep[-1] = '\t';
        }
        if (provref_col) {
          *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, row_variant_uidx)))? 'Y' : 'N';
          *cswritep++ = '\t';
        }
        if (maj_col || nonmaj_col || freq_col) {
          if (maj_alleles) {
            maj_allele_idx = maj_alleles[row_variant_uidx];
          }
          if (maj_col) {
            cswritep = strcpyax(cswritep, cur_row_alleles[maj_allele_idx], '\t');
          }
          if (nonmaj_col) {
            for (uint32_t allele_idx = 0; allele_idx != row_allele_ct; ++allele_idx) {
              if (allele_idx == maj_allele_idx) {
                continue;
              }
              if (unlikely(Cswrite(cssp, &cswritep))) {
                goto VcorTableWriteThread_ret_WRITE_FAIL;
              }
              cswritep = strcpyax(cswritep, cur_row_alleles[allele_idx], ',');
            }
            cswritep[-1] = '\t';
          }
          if (freq_col) {
            const double maj_freq = GetAlleleFreq(&(allele_freqs[row_allele_idx_offset_base - row_variant_uidx]), maj_allele_idx, row_allele_ct);
            cswritep = dtoa_g(1.0 - maj_freq, cswritep);
            *cswritep++ = '\t';
          }
        }

        if (chr_col) {
          cswritep = memcpya(cswritep, col_chr_buf, col_chr_blen);
        }
        if (pos_col) {
          cswritep = u32toa_x(variant_bps[col_variant_uidx], '\t', cswritep);
        }
        if (id_col) {
          cswritep = strcpyax(cswritep, variant_ids[col_variant_uidx], '\t');
        }
        uintptr_t col_allele_idx_offset_base = col_variant_uidx * 2;
        if (allele_idx_offsets) {
          col_allele_idx_offset_base = allele_idx_offsets[col_variant_uidx];
          col_allele_ct = allele_idx_offsets[col_variant_uidx + 1] - col_allele_idx_offset_base;
        }
        const char* const* cur_col_alleles = &(allele_storage[col_allele_idx_offset_base]);
        if (ref_col) {
          cswritep = strcpyax(cswritep, cur_col_alleles[0], '\t');
        }
        if (alt1_col) {
          cswritep = strcpyax(cswritep, cur_col_alleles[1], '\t');
        }
        if (alt_col) {
          for (uint32_t allele_idx = 1; allele_idx != col_allele_ct; ++allele_idx) {
            if (unlikely(Cswrite(cssp, &cswritep))) {
              goto VcorTableWriteThread_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, cur_col_alleles[allele_idx], ',');
          }
          cswritep[-1] = '\t';
        }
        if (provref_col) {
          *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, col_variant_uidx)))? 'Y' : 'N';
          *cswritep++ = '\t';
        }
        if (maj_col || nonmaj_col || freq_col) {
          if (maj_alleles) {
            maj_allele_idx = maj_alleles[col_variant_uidx];
          }
          if (maj_col) {
            cswritep = strcpyax(cswritep, cur_col_alleles[maj_allele_idx], '\t');
          }
          if (nonmaj_col) {
            for (uint32_t allele_idx = 0; allele_idx != col_allele_ct; ++allele_idx) {
              if (allele_idx == maj_allele_idx) {
                continue;
              }
              if (unlikely(Cswrite(cssp, &cswritep))) {
                goto VcorTableWriteThread_ret_WRITE_FAIL;
              }
              cswritep = strcpyax(cswritep, cur_col_alleles[allele_idx], ',');
            }
            cswritep[-1] = '\t';
          }
          if (freq_col) {
            const double maj_freq = GetAlleleFreq(&(allele_freqs[col_allele_idx_offset_base - col_variant_uidx]), maj_allele_idx, col_allele_ct);
            cswritep = dtoa_g(1.0 - maj_freq, cswritep);
            *cswritep++ = '\t';
          }
        }

        // R or R2
        cswritep = dtoa_g(r_or_r2, cswritep);
        if (d_col) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(*cur_results_iter++, cswritep);
        }
        if (dprime_col) {
          *cswritep++ = '\t';
          double dprime = *cur_results_iter++;
          if (dprime_abs) {
            dprime = fabs(dprime);
          }
          cswritep = dtoa_g(dprime, cswritep);
        }
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(cssp, &cswritep))) {
          goto VcorTableWriteThread_ret_WRITE_FAIL;
        }
      }
    }
    ctx->cswritep = cswritep;
    row_parity = 1 - row_parity;
    while (0) {
    VcorTableWriteThread_ret_WRITE_FAIL:
      ctx->reterr = kPglRetWriteFail;
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

// Determine --ld-window-kb, --ld-window-cm, and --ld-window intersection.
// Assumes chr_start_uidx <= uidx_start <= uidx_end on entry, and that
// uidx_start and uidx_end are not greater than their true values for the
// current variant.
void UpdateVcorWindow(const uintptr_t* variant_include, const uint32_t* variant_bps, const double* variant_cms, uint32_t row_snp_subset, uint32_t var_ct_radius, uint32_t bp_radius, double cm_radius, uint32_t chr_end_uidx, uint32_t center_uidx, uint32_t* uidx_startp, uint32_t* uidx_endp) {
  const uint32_t cur_bp = variant_bps[center_uidx];
  if (row_snp_subset) {
    uint32_t uidx_start = *uidx_startp;
    if (cur_bp > bp_radius) {
      uidx_start = ExpsearchU32(variant_bps, uidx_start, center_uidx, cur_bp - bp_radius);
    }
    if (variant_cms) {
      const double cur_cm = variant_cms[center_uidx];
      uidx_start = ExpsearchD(variant_cms, uidx_start, center_uidx, cur_cm - cm_radius);
    }
    // var_ct_radius <= 0x7fffffff, so no overflow risk
    if (uidx_start + var_ct_radius < center_uidx) {
      const uint32_t leading_var_ct = PopcountBitRange(variant_include, uidx_start, center_uidx);
      if (leading_var_ct > var_ct_radius) {
        uidx_start = FindNth1BitFrom(variant_include, uidx_start + 1, leading_var_ct - var_ct_radius);
      }
    }
    *uidx_startp = uidx_start;
  } else {
    *uidx_startp = center_uidx;
  }

  uint32_t uidx_end = MAXV(*uidx_endp, center_uidx + 1);
  if (uidx_end < chr_end_uidx) {
    uint32_t uidx_end_ceil = ExpsearchU32(variant_bps, uidx_end, chr_end_uidx, cur_bp + bp_radius + 1);
    if (variant_cms) {
      const double cur_cm = variant_cms[center_uidx];
      uidx_end_ceil = ExpsearchD(variant_cms, uidx_end, uidx_end_ceil, cur_cm + cm_radius);
    }
    if (center_uidx + 1 + var_ct_radius < uidx_end_ceil) {
      const uint32_t trailing_var_ct = PopcountBitRange(variant_include, center_uidx + 1, uidx_end_ceil);
      if (trailing_var_ct > var_ct_radius) {
        uidx_end_ceil = 1 + FindNth1BitFrom(variant_include, center_uidx + 1, var_ct_radius);
      }
    }
    uidx_end = uidx_end_ceil;
  }
  *uidx_endp = uidx_end;
}

PglErr VcorTable(const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const double* variant_cms, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const VcorInfo* vcip, const char* flagname, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  VcorTableCtx ctx;
  VcorTableWriteCtx write_ctx;
  ThreadGroup tg;
  ThreadGroup write_tg;
  PreinitCstream(&write_ctx.css);
  write_ctx.cswritep = nullptr;
  PreinitThreads(&tg);
  PreinitThreads(&write_tg);
  {
    const VcorFlags flags = vcip->flags;
    const uint32_t is_unsquared = (flags / kfVcorUnsquared) & 1;
    const uint32_t phased_calc = (flags / kfVcorPhased) & 1;
    const uint32_t ref_based = (flags / kfVcorRefBased) & 1;
    if (!(flags & kfVcorAllowAmbiguousAllele)) {
      if (is_unsquared) {
        uint32_t is_ambiguous_biallelic = 0;
        if (!ref_based) {
          is_ambiguous_biallelic = !(flags & (kfVcorColMaj | kfVcorColNonmaj));
        } else {
          const VcorFlags relevant_allele_cols = flags & (kfVcorColRef | kfVcorColAlt1 | kfVcorColAlt);
          if (relevant_allele_cols != kfVcorColAlt1) {
            is_ambiguous_biallelic = (relevant_allele_cols == kfVcor0);
          } else {
            if (unlikely(MultiallelicVariantPresent(orig_variant_include, allele_idx_offsets, orig_variant_ct))) {
              logerrprintfww("Error: The meaning of r's sign cannot be consistently inferred from just the %s 'alt1' column-set at multiallelic variants. Either filter out multiallelic variants, revise the column-set, or use the 'allow-ambiguous-allele' modifier to override this error.\n", flagname);
              return kPglRetInconsistentInput;
            }
          }
        }
        if (unlikely(is_ambiguous_biallelic)) {
          logerrprintfww("Error: %s column-set doesn't include allele columns which clarify the meaning of r's sign. Either switch to --r2-%sphased, add a disambiguating column-set, or use the 'allow-ambiguous-allele' modifier to override this error.\n", flagname, phased_calc? "" : "un");
          return kPglRetInconsistentInput;
        }
      } else {
        uint32_t is_ambiguous_multiallelic;
        if (!ref_based) {
          is_ambiguous_multiallelic = !(flags & (kfVcorColMaj | kfVcorColNonmaj));
        } else {
          is_ambiguous_multiallelic = !(flags & (kfVcorColRef | kfVcorColAlt));
        }
        if (unlikely(is_ambiguous_multiallelic && MultiallelicVariantPresent(orig_variant_include, allele_idx_offsets, orig_variant_ct))) {
          logerrprintfww("Error: %s column-set doesn't include allele columns which clarify which calculation is being performed at multiallelic variants. Either filter out multiallelic variants, revise the column-set (with e.g. \"cols=+%s\"), or use the 'allow-ambiguous-allele' modifier to override this error.\n", flagname, ref_based? "ref" : "maj");
          return kPglRetInconsistentInput;
        }
      }
    }
    // 1. If not inter-chr, remove unplaced variants.
    // 2. If --parallel and/or --ld-snp/--ld-snps/--ld-snp-list, initialize
    //    row_variant_include, then try to shrink variant_include:
    //    - If inter-chr + --parallel, and this isn't the first piece, we can
    //      remove leading variants.
    //    - If not inter-chr, we can remove all variants that are out of range
    //      of row_variant_include elements.
    const uint32_t inter_chr = (flags / kfVcorInterChr) & 1;
    const char* ld_snp_list_fname = vcip->ld_snp_list_fname;
    const RangeList* ld_snp_range_listp = &(vcip->ld_snp_range_list);
    const uint32_t row_snp_subset = ld_snp_list_fname || (ld_snp_range_listp->name_ct != 0);
    const double min_r2 = vcip->min_r2;
    if (unlikely(inter_chr && (!row_snp_subset) && (min_r2 <= 0.0) && (orig_variant_ct > 400000) && (parallel_tot == 1) && (!(flags & kfVcorYesReally)))) {
      logerrprintfww("Error: Gigantic (over 400k variants) %s unfiltered, non-distributed computation. Rerun with the 'yes-really' modifier if you are SURE you have enough hard drive space and want to do this.\n", flagname);
      goto VcorTable_ret_INCONSISTENT_INPUT;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    // There are conditions under which we can optimize out these allocations,
    // but we don't realistically need the memory under those conditions, so
    // just go ahead and make copies 100% of the time.
    uintptr_t* variant_include_buf;
    uintptr_t* row_variant_include_buf;
    if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include_buf) ||
                 bigstack_alloc_w(raw_variant_ctl, &row_variant_include_buf))) {
      goto VcorTable_ret_NOMEM;
    }
    memcpy(variant_include_buf, orig_variant_include, raw_variant_ctl * sizeof(intptr_t));
    memcpy(row_variant_include_buf, orig_variant_include, raw_variant_ctl * sizeof(intptr_t));
    uintptr_t* row_subset_exclude_buf = nullptr;
    if (row_snp_subset && (parallel_idx != 0)) {
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, &row_subset_exclude_buf))) {
        goto VcorTable_ret_NOMEM;
      }
    }
    uint32_t variant_ct = orig_variant_ct;
    if (!inter_chr) {
      uint32_t skipped_variant_ct = StripUnplacedMut(cip, variant_include_buf);
      if (skipped_variant_ct) {
        logprintf("%s: Ignoring %u chromosome 0 variant%s.\n", flagname, skipped_variant_ct, (skipped_variant_ct == 1)? "" : "s");
        variant_ct -= skipped_variant_ct;
      }
    }
    const uint32_t var_ct_radius = vcip->var_ct_radius;
    const uint32_t bp_radius = vcip->bp_radius;
    const double cm_radius = vcip->cm_radius;
    if (cm_radius == -1.0) {
      variant_cms = nullptr;
    }
    uint32_t row_variant_ct = variant_ct;
    if ((parallel_tot != 1) || row_snp_subset) {
      if (row_snp_subset) {
        unsigned char* bigstack_mark2 = g_bigstack_base;
        uint32_t* variant_id_htable;
        uint32_t* htable_dup_base;
        uint32_t variant_id_htable_size;
        AllocAndPopulateIdHtableMt(variant_include_buf, variant_ids, variant_ct, bigstack_left() / 2, max_thread_ct, &variant_id_htable, &htable_dup_base, &variant_id_htable_size, nullptr);
        if (ld_snp_list_fname) {
          const uint32_t fname_slen = strlen(ld_snp_list_fname);
          char* fnames_tmp;
          if (unlikely(bigstack_alloc_c(fname_slen + 2, &fnames_tmp))) {
            goto VcorTable_ret_NOMEM;
          }
          memcpy(fnames_tmp, ld_snp_list_fname, fname_slen + 1);
          fnames_tmp[fname_slen + 1] = '\0';
          reterr = TokenExtractExclude(variant_ids, variant_id_htable, htable_dup_base, fnames_tmp, "ld-snp-list", raw_variant_ct, max_variant_id_slen, variant_id_htable_size, kVfilterExtract, max_thread_ct, row_variant_include_buf, &row_variant_ct);
          if (unlikely(reterr)) {
            goto VcorTable_ret_1;
          }
        } else {
          uintptr_t* seen_uidxs;
          if (unlikely(bigstack_calloc_w(raw_variant_ctl, &seen_uidxs))) {
            goto VcorTable_ret_NOMEM;
          }
          reterr = InterpretVariantRangeList(variant_ids, variant_id_htable, htable_dup_base, ld_snp_range_listp, "--ld-snps", max_variant_id_slen, variant_id_htable_size, seen_uidxs);
          if (unlikely(reterr)) {
            goto VcorTable_ret_1;
          }
          BitvecAnd(seen_uidxs, raw_variant_ctl, row_variant_include_buf);
          row_variant_ct = PopcountWords(row_variant_include_buf, raw_variant_ctl);
        }
        BigstackReset(bigstack_mark2);
      }
      if (parallel_tot > 1) {
        const uint32_t row_shard_start = (S_CAST(uint64_t, row_variant_ct) * parallel_idx) / parallel_tot;
        const uint32_t row_uidx_start = IdxToUidxBasic(row_variant_include_buf, row_shard_start);
        if ((parallel_idx != 0) && row_snp_subset) {
          memcpy(row_subset_exclude_buf, row_variant_include_buf, raw_variant_ctl * sizeof(intptr_t));
        }
        if (row_uidx_start) {
          ClearBitsNz(0, row_uidx_start, row_variant_include_buf);
        }
        uint32_t row_shard_end = row_variant_ct;
        if (parallel_idx + 1 != parallel_tot) {
          row_shard_end = (S_CAST(uint64_t, row_variant_ct) * (parallel_idx + 1)) / parallel_tot;
          const uint32_t row_vecidx_start = row_uidx_start / kBitsPerVec;
          const uint32_t row_uidx_end = row_vecidx_start * kBitsPerVec + IdxToUidxBasic(&(row_variant_include_buf[row_vecidx_start * kWordsPerVec]), row_shard_end - row_shard_start);
          ClearBitsNz(row_uidx_end, raw_variant_ct, row_variant_include_buf);
        }
        row_variant_ct = row_shard_end - row_shard_start;
        if ((row_shard_start != 0) && (!row_snp_subset)) {
          // Safe to remove column-variants before the current row-shard.
          variant_ct -= PopcountBitRange(variant_include_buf, 0, row_uidx_start);
          ClearBitsNz(0, row_uidx_start, variant_include_buf);
        }
      }
    }
    if (!inter_chr) {
      // Safe to remove column-variants out of range of all row variants, as
      // well as "orphan" row variants.
      uint32_t row_chr_fo_idx = UINT32_MAX;
      uint32_t row_chr_end = 0;

      uint32_t uidx_start = 0;
      uint32_t uidx_end = 0;

      uintptr_t row_variant_uidx_base = 0;
      uintptr_t row_cur_bits = row_variant_include_buf[0];
      for (uint32_t row_variant_idx = 0; row_variant_idx != row_variant_ct; ++row_variant_idx) {
        const uint32_t row_variant_uidx = BitIter1(row_variant_include_buf, &row_variant_uidx_base, &row_cur_bits);
        if (row_variant_uidx >= row_chr_end) {
          do {
            ++row_chr_fo_idx;
            row_chr_end = cip->chr_fo_vidx_start[row_chr_fo_idx + 1];
          } while (row_variant_uidx >= row_chr_end);
          const uint32_t row_chr_start = cip->chr_fo_vidx_start[row_chr_fo_idx];
          uidx_start = row_chr_start;
          uidx_end = row_chr_start;
        }
        const uint32_t prev_uidx_end = uidx_end;
        UpdateVcorWindow(variant_include_buf, variant_bps, variant_cms, row_snp_subset, var_ct_radius, bp_radius, cm_radius, row_chr_end, row_variant_uidx, &uidx_start, &uidx_end);
        if (prev_uidx_end < uidx_start) {
          ClearBitsNz(prev_uidx_end, uidx_start, variant_include_buf);
        }
        if (row_snp_subset && (uidx_start + 1 == uidx_end)) {
          // No variants in range.  This row-variant is an "orphan" and can be
          // skipped.
          ClearBit(uidx_start, variant_include_buf);
        }
      }
      BitvecAnd(variant_include_buf, raw_variant_ctl, row_variant_include_buf);
      variant_ct = PopcountWords(variant_include_buf, raw_variant_ctl);
    }
    const uintptr_t* variant_include = variant_include_buf;
    const uintptr_t* row_variant_include = row_variant_include_buf;
    const uint32_t* variant_include_cumulative_popcounts;
    const uint32_t* row_variant_include_cumulative_popcounts;
    {
      // bugfix (16 May 2024): there are RawToSubsettedPos(variant_include,
      // variant_include_cumulative_popcounts, x) calls with x =
      // raw_variant_ct.  In this case, we may need one more entry.
      uint32_t* variant_include_cumulative_popcounts_buf;
      uint32_t* row_variant_include_cumulative_popcounts_buf;
      if (unlikely(bigstack_alloc_u32(1 + (raw_variant_ct / kBitsPerWord), &variant_include_cumulative_popcounts_buf) ||
                   bigstack_alloc_u32(raw_variant_ctl, &row_variant_include_cumulative_popcounts_buf))) {
        goto VcorTable_ret_NOMEM;
      }

      FillCumulativePopcounts(variant_include, raw_variant_ctl, variant_include_cumulative_popcounts_buf);
      if ((raw_variant_ct % kBitsPerWord) == 0) {
        variant_include_cumulative_popcounts_buf[raw_variant_ctl] = variant_ct;
      }
      FillCumulativePopcounts(row_variant_include, raw_variant_ctl, row_variant_include_cumulative_popcounts_buf);
      variant_include_cumulative_popcounts = variant_include_cumulative_popcounts_buf;
      row_variant_include_cumulative_popcounts = row_variant_include_cumulative_popcounts_buf;
    }
    row_variant_ct = row_variant_include_cumulative_popcounts[raw_variant_ctl - 1] + PopcountWord(row_variant_include[raw_variant_ctl - 1]);
    // variant_ct == 0 is possible here.  In that case, we want to write the
    // header line (or nothing at all if parallel_idx > 0) and skip the main
    // loop, not error out.

    const uintptr_t* nonref_flags = PgrGetNonrefFlags(simple_pgrp);
    // "&& (!nonref_flags)" needed since this is after --ref-allele, etc. in
    // the order of operations.
    const uint32_t all_nonref = (PgrGetGflags(simple_pgrp) & kfPgenGlobalAllNonref) && (!nonref_flags);
    write_ctx.all_nonref = all_nonref;
    const uint32_t provref_col = (flags & kfVcorColRef) && ProvrefCol(variant_include, nonref_flags, flags / kfVcorColMaybeprovref, raw_variant_ct, all_nonref);
    const uint32_t d_col = (flags / kfVcorColD) & 1;
    const uint32_t dprime_col = ((flags & (kfVcorColDprime | kfVcorColDprimeAbs)) != 0);
    write_ctx.row_chr_buf = nullptr;
    write_ctx.col_chr_buf = nullptr;
    {
      uintptr_t overflow_buf_size = kCompressStreamBlock + 512;
      const uint32_t chr_col = (flags / kfVcorColChrom) & 1;
      const uint32_t pos_col = (flags / kfVcorColPos) & 1;
      const uint32_t id_col = (flags / kfVcorColId) & 1;
      const uint32_t ref_col = (flags / kfVcorColRef) & 1;
      const uint32_t alt1_col = (flags / kfVcorColAlt1) & 1;
      const uint32_t alt_col = (flags / kfVcorColAlt) & 1;
      // provref_col defined earlier
      const uint32_t maj_col = (flags / kfVcorColMaj) & 1;
      const uint32_t nonmaj_col = (flags / kfVcorColNonmaj) & 1;
      const uint32_t freq_col = (flags / kfVcorColFreq) & 1;
      // d_col, dprime_col defined earlier
      if (chr_col) {
        const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
        if (unlikely(bigstack_alloc_c(max_chr_blen, &write_ctx.row_chr_buf) ||
                     bigstack_alloc_c(max_chr_blen, &write_ctx.col_chr_buf))) {
          goto VcorTable_ret_NOMEM;
        }
        overflow_buf_size += 2 * max_chr_blen;
      }
      if (id_col) {
        overflow_buf_size += 2 * (max_variant_id_slen + 1);
      }
      uintptr_t n_allele = ref_col + alt1_col + maj_col;
      if (!n_allele) {
        if (alt_col || nonmaj_col) {
          n_allele = 1;
        }
      }
      overflow_buf_size += n_allele * (max_allele_slen + 1);
      const uint32_t output_zst = (flags / kfVcorZs) & 1;
      char* outname_write_iter = strcpya_k(outname_end, ".vcor");
      if (parallel_tot != 1) {
        *outname_write_iter++ = '.';
        outname_write_iter = u32toa(parallel_idx + 1, outname_write_iter);
      }
      uint32_t compress_thread_ct = 1;
      if (output_zst) {
        outname_write_iter = strcpya_k(outname_write_iter, ".zst");
        // more room to tune this, but this is an easy win
        if ((!phased_calc) && (min_r2 <= 0.0) && (max_thread_ct > 4) && (founder_ct <= 65536)) {
          compress_thread_ct = 2 + ((max_thread_ct > 8) && (founder_ct <= 32768));
        }
      }
      *outname_write_iter = '\0';
      reterr = InitCstreamAlloc(outname, 0, output_zst, compress_thread_ct, overflow_buf_size, &write_ctx.css, &write_ctx.cswritep);
      if (unlikely(reterr)) {
        goto VcorTable_ret_1;
      }
      if (parallel_idx == 0) {
        char* cswritep = write_ctx.cswritep;
        *cswritep++ = '#';
        if (chr_col) {
          cswritep = strcpya_k(cswritep, "CHROM_A\t");
        }
        if (pos_col) {
          cswritep = strcpya_k(cswritep, "POS_A\t");
        }
        if (id_col) {
          cswritep = strcpya_k(cswritep, "ID_A\t");
        }
        if (ref_col) {
          cswritep = strcpya_k(cswritep, "REF_A\t");
        }
        if (alt1_col) {
          cswritep = strcpya_k(cswritep, "ALT1_A\t");
        }
        if (alt_col) {
          cswritep = strcpya_k(cswritep, "ALT_A\t");
        }
        if (provref_col) {
          cswritep = strcpya_k(cswritep, "PROVISIONAL_REF_A?\t");
        }
        if (maj_col) {
          cswritep = strcpya_k(cswritep, "MAJ_A\t");
        }
        if (nonmaj_col) {
          cswritep = strcpya_k(cswritep, "NONMAJ_A\t");
        }
        if (freq_col) {
          cswritep = strcpya_k(cswritep, "NONMAJ_FREQ_A\t");
        }

        if (chr_col) {
          cswritep = strcpya_k(cswritep, "CHROM_B\t");
        }
        if (pos_col) {
          cswritep = strcpya_k(cswritep, "POS_B\t");
        }
        if (id_col) {
          cswritep = strcpya_k(cswritep, "ID_B\t");
        }
        if (ref_col) {
          cswritep = strcpya_k(cswritep, "REF_B\t");
        }
        if (alt1_col) {
          cswritep = strcpya_k(cswritep, "ALT1_B\t");
        }
        if (alt_col) {
          cswritep = strcpya_k(cswritep, "ALT_B\t");
        }
        if (provref_col) {
          cswritep = strcpya_k(cswritep, "PROVISIONAL_REF_B?\t");
        }
        if (maj_col) {
          cswritep = strcpya_k(cswritep, "MAJ_B\t");
        }
        if (nonmaj_col) {
          cswritep = strcpya_k(cswritep, "NONMAJ_B\t");
        }
        if (freq_col) {
          cswritep = strcpya_k(cswritep, "NONMAJ_FREQ_B\t");
        }

        if (!phased_calc) {
          cswritep = strcpya_k(cswritep, "UN");
        }
        cswritep = strcpya_k(cswritep, "PHASED_R");
        if (!is_unsquared) {
          *cswritep++ = '2';
        }

        if (d_col) {
          cswritep = strcpya_k(cswritep, "\tD");
        }
        if (dprime_col) {
          *cswritep++ = '\t';
          if (flags & kfVcorColDprimeAbs) {
            cswritep = strcpya_k(cswritep, "ABS_");
          }
          cswritep = strcpya_k(cswritep, "DPRIME");
        }
        AppendBinaryEoln(&cswritep);

        write_ctx.cswritep = cswritep;
      }
    }
    if (!row_variant_ct) {
      logprintf("%s: No variant-pairs to process.\n", flagname);
      goto VcorTable_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
    const uint32_t founder_ctv = BitCtToVecCt(founder_ct);
    const uint32_t founder_ctv2 = NypCtToVecCt(founder_ct);
    const uint32_t founder_ctaw = founder_ctv * kWordsPerVec;
    const uint32_t founder_male_ct = PopcountWordsIntersect(founder_info, sex_male, raw_sample_ctl);
    const uint32_t all_haploid = IsSet(cip->haploid_mask, 0);
    PgenGlobalFlags effective_gflags = PgrGetGflags(simple_pgrp) & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
    const uint32_t check_phase = phased_calc && (!all_haploid) && (effective_gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent));
    if (!check_phase) {
      effective_gflags &= kfPgenGlobalDosagePresent;
    }
    const R2PhaseType phase_type = GetR2PhaseType(phased_calc, check_phase);
    const uint32_t check_dosage = (effective_gflags / kfPgenGlobalDosagePresent) & 1;

    uintptr_t* raregeno = nullptr;
    uint32_t* difflist_sample_ids = nullptr;
    const uint32_t max_difflist_len = founder_ct / 64;
    if (!phased_calc) {
      if (unlikely(bigstack_alloc_w(NypCtToWordCt(max_difflist_len), &raregeno) ||
                   bigstack_alloc_u32(max_difflist_len, &difflist_sample_ids))) {
        goto VcorTable_ret_NOMEM;
      }
    }

    const uintptr_t* founder_male_collapsed = nullptr;
    const Dosage* male_dosage_invmask = nullptr;
    ctx.founder_nonmale_collapsed = nullptr;
    ctx.nonmale_dosage_invmask = nullptr;
    ctx.chrx_idx = UINT32_MAX;
    uint32_t x_code = UINT32_MAX;
    // if all-males, we can ignore phase on chrX, as well as skipping
    // male/nonmale-specific stats
    // if all-nonmales, we initialize x_code to prevent phase from being
    // ignored, but can skip the male/nonmale-specific stats
    if (founder_male_ct != founder_ct) {
      if (XymtExists(cip, kChrOffsetX, &x_code) && (founder_male_ct != 0)) {
        const uint32_t x_fo_idx = cip->chr_idx_to_foidx[x_code];
        const uint32_t start_vidx = cip->chr_fo_vidx_start[x_fo_idx];
        const uint32_t end_vidx = cip->chr_fo_vidx_start[x_fo_idx + 1];
        if (!AllBitsAreZero(variant_include, start_vidx, end_vidx)) {
          ctx.chrx_idx = x_code;
          uintptr_t* founder_male_collapsed_fill;
          uintptr_t* founder_nonmale_collapsed;
          if (unlikely(bigstack_alloc_w(founder_ctl, &founder_male_collapsed_fill) ||
                       bigstack_alloc_w(founder_ctl, &founder_nonmale_collapsed))) {
            goto VcorTable_ret_NOMEM;
          }
          CopyBitarrSubset(sex_male, founder_info, founder_ct, founder_male_collapsed_fill);
          founder_male_collapsed = founder_male_collapsed_fill;
          BitvecInvertCopy(founder_male_collapsed, founder_ctl, founder_nonmale_collapsed);
          ZeroTrailingBits(founder_ct, founder_nonmale_collapsed);
          ctx.founder_nonmale_collapsed = founder_nonmale_collapsed;
          if (check_dosage) {
            const uint32_t founder_ctad = RoundUpPow2(founder_ct, kDosagePerVec);
            Dosage* male_dosage_invmask_fill;
            Dosage* nonmale_dosage_invmask;
            if (bigstack_alloc_dosage(founder_ctad, &male_dosage_invmask_fill) ||
                bigstack_alloc_dosage(founder_ctad, &nonmale_dosage_invmask)) {
              goto VcorTable_ret_NOMEM;
            }
            Expand1bitTo16(founder_male_collapsed, founder_ctad, 0xffff, male_dosage_invmask_fill);
            Expand1bitTo16(founder_nonmale_collapsed, founder_ctad, 0xffff, nonmale_dosage_invmask);
            male_dosage_invmask = male_dosage_invmask_fill;
            ctx.nonmale_dosage_invmask = nonmale_dosage_invmask;
          }
        }
      }
    }
    const uint32_t x_exists = (ctx.chrx_idx < UINT32_MAXM1);
    ctx.founder_male_collapsed = founder_male_collapsed;
    ctx.male_dosage_invmask = male_dosage_invmask;

    uint32_t y_start;
    uint32_t y_end;
    GetXymtStartAndEnd(cip, kChrOffsetY, &y_start, &y_end);
    uintptr_t* founder_female_collapsed = nullptr;
    uintptr_t* founder_female_collapsed_interleaved = nullptr;
    if (y_end) {
      if ((founder_male_ct == founder_ct) || AllBitsAreZero(variant_include, y_start, y_end)) {
        y_start = 0;
        y_end = 0;
      }
      if (y_end) {
        uintptr_t* founder_female;
        if (bigstack_end_alloc_w(raw_sample_ctl, &founder_female)) {
          goto VcorTable_ret_NOMEM;
        }
        BitvecInvmaskCopy(sex_nm, sex_male, raw_sample_ctl, founder_female);
        BitvecAnd(founder_info, raw_sample_ctl, founder_female);
        if (AllWordsAreZero(founder_female, raw_sample_ctl)) {
          y_start = 0;
          y_end = 0;
        } else {
          if (bigstack_alloc_w(founder_ctaw, &founder_female_collapsed) ||
              bigstack_alloc_w(founder_ctaw, &founder_female_collapsed_interleaved)) {
            goto VcorTable_ret_NOMEM;
          }
          CopyBitarrSubset(founder_female, founder_info, founder_ct, founder_female_collapsed);
          ZeroTrailingWords(founder_ctl, founder_female_collapsed);
          FillInterleavedMaskVec(founder_female_collapsed, founder_ctv, founder_female_collapsed_interleaved);
        }
        BigstackEndReset(bigstack_end_mark);
      }
    }
    ctx.founder_ct = founder_ct;
    ctx.founder_male_ct = founder_male_ct;
    ctx.is_unsquared = is_unsquared;

    const uintptr_t bitvec_byte_ct = BitCtToVecCt(founder_ct) * kBytesPerVec;
    uintptr_t dosagevec_byte_ct = 0;
    uintptr_t unpacked_variant_byte_stride;
    if (check_dosage) {
      dosagevec_byte_ct = DivUp(founder_ct, kDosagePerVec) * kBytesPerVec;
      const uintptr_t dosage_trail_byte_ct = LdDosageTrailAlignedByteCt(S_CAST(R2PhaseType, phased_calc), x_exists);
      unpacked_variant_byte_stride = dosagevec_byte_ct * (1 + phased_calc + check_phase) + bitvec_byte_ct + dosage_trail_byte_ct;
    } else {
      unpacked_variant_byte_stride = RoundUpPow2(16, kBytesPerVec) + bitvec_byte_ct * (3 + 2 * check_phase);
#ifndef USE_AVX2
      const uintptr_t sparse_req = RoundUpPow2((6 + max_difflist_len) * sizeof(int32_t), kBytesPerVec) + NypCtToVecCt(max_difflist_len) * kBytesPerVec + RoundUpPow2(founder_ctl * (kBytesPerWord + sizeof(int32_t)), kBytesPerVec);
      if (sparse_req > unpacked_variant_byte_stride) {
        unpacked_variant_byte_stride = sparse_req;
      }
#endif
      const uintptr_t nondosage_trail_byte_ct = LdNondosageTrailAlignedByteCt(S_CAST(R2PhaseType, phased_calc), x_exists);
      unpacked_variant_byte_stride += nondosage_trail_byte_ct;
    }
    uint32_t* founder_info_cumulative_popcounts;
    PgenVariant pgv;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &founder_info_cumulative_popcounts) ||
                 BigstackAllocPgv(founder_ct, 0, effective_gflags, &pgv))) {
      goto VcorTable_ret_NOMEM;
    }
    FillCumulativePopcounts(founder_info, raw_sample_ctl, founder_info_cumulative_popcounts);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(founder_info_cumulative_popcounts, simple_pgrp, &pssi);

    ctx.phase_type = phase_type;
    ctx.check_dosage = check_dosage;
    ctx.report_d = d_col;
    ctx.report_dprime = dprime_col;
    ctx.unpacked_variant_byte_stride = unpacked_variant_byte_stride;
    // Some room to tune this, e.g. if .zst compression is on, founder_ct is
    // small, and computation is unphased, we may want more than 1 compression
    // thread.  But this should be good enough for now.
    uint32_t calc_thread_ct = max_thread_ct - (max_thread_ct > 4) - (max_thread_ct > 8);
    ctx.cur_nm_bufs = nullptr;
    ctx.invmask_bufs = nullptr;
    if (x_exists) {
      if (unlikely(bigstack_alloc_wp(calc_thread_ct, &ctx.cur_nm_bufs) ||
                   bigstack_alloc_dosagep(calc_thread_ct, &ctx.invmask_bufs))) {
        goto VcorTable_ret_NOMEM;
      }
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        if (unlikely(bigstack_alloc_w(founder_ctl, &(ctx.cur_nm_bufs[tidx])) ||
                     bigstack_alloc_dosage(founder_ct, &(ctx.invmask_bufs[tidx])))) {
          goto VcorTable_ret_NOMEM;
        }
      }
    }

    // Now initialize the parts of write_ctx that don't depend on uv_capacity.
    write_ctx.variant_include = variant_include;
    write_ctx.variant_include_cumulative_popcounts = variant_include_cumulative_popcounts;
    write_ctx.row_variant_include = row_variant_include;
    write_ctx.row_variant_include_cumulative_popcounts = row_variant_include_cumulative_popcounts;
    write_ctx.row_subset_exclude = nullptr;
    if (row_snp_subset) {
      write_ctx.row_subset_exclude = (parallel_idx == 0)? row_variant_include : row_subset_exclude_buf;
    }
    write_ctx.cip = cip;
    write_ctx.variant_bps = variant_bps;
    write_ctx.variant_ids = variant_ids;
    write_ctx.allele_idx_offsets = allele_idx_offsets;
    write_ctx.allele_storage = allele_storage;
    write_ctx.nonref_flags = nonref_flags;
    write_ctx.maj_alleles = maj_alleles;
    write_ctx.allele_freqs = allele_freqs;
    if (!is_unsquared) {
      write_ctx.r_or_r2_thresh = min_r2;
    } else {
      write_ctx.r_or_r2_thresh = (min_r2 < 0.0)? -1.0 : sqrt(min_r2);
    }
    write_ctx.raw_variant_ct = raw_variant_ct;
    write_ctx.flags = flags;
    write_ctx.provref_col = provref_col;
    write_ctx.reterr = kPglRetSuccess;

    // Assign ~half of remaining workspace to input buffers and ~half to output
    // buffers.
    const uintptr_t results_stride = 1 + d_col + dprime_col;
    const uintptr_t result_capacity = bigstack_left() / (4 * sizeof(double) * results_stride);
    if (unlikely(bigstack_alloc_d(result_capacity * results_stride, &(ctx.results[0])) ||
                 bigstack_alloc_d(result_capacity * results_stride, &(ctx.results[1])))) {
      goto VcorTable_ret_NOMEM;
    }
    write_ctx.results[0] = ctx.results[0];
    write_ctx.results[1] = ctx.results[1];

    uint32_t* col_offset_starts_storage[2];
    uintptr_t* write_idx_starts_storage[2];
    uint32_t row_window_sizes[2];
    ctx.row_window_sizes = row_window_sizes;
    write_ctx.row_window_sizes = row_window_sizes;
    // All remaining allocation sizes are a function of uv_capacity.  Set its
    // value.
    uintptr_t uv_capacity;
    {
      // ctx:
      //   unpacked_variants: 2 * unpacked_variant_byte_stride
      //   (no row_uvidxs)
      //   col_uvidxs: 2 * sizeof(int32_t)
      //   row_chr_idxs: 2 * sizeof(ChrIdx)
      //   col_chr_idxs: 2 * sizeof(ChrIdx)
      //   col_offset_starts: 2 * sizeof(int32_t)
      //   write_idx_starts: 2 * sizeof(intptr_t) (length (uv_capacity + 1))
      // (could optimize row_chr_idxs and col_chr_idxs out)
      uintptr_t bytes_left = bigstack_left();
      // defend against adverse rounding
      if (unlikely(bytes_left < 11 * kCacheline)) {
        goto VcorTable_ret_NOMEM;
      }
      bytes_left -= 11 * kCacheline;
      const uintptr_t bytes_per_unpacked_variant = 2 * (unpacked_variant_byte_stride + 2 * sizeof(int32_t) + 2 * sizeof(ChrIdx) + sizeof(intptr_t));
      uv_capacity = bytes_left / bytes_per_unpacked_variant;
      if (unlikely(uv_capacity < 2)) {
        goto VcorTable_ret_NOMEM;
      }
      if (uv_capacity > variant_ct) {
        uv_capacity = variant_ct;
      }
      // shouldn't be possible for these allocations to fail
      if (unlikely(bigstack_alloc_uc(uv_capacity * unpacked_variant_byte_stride, &(ctx.unpacked_variants[0])) ||
                   bigstack_alloc_uc(uv_capacity * unpacked_variant_byte_stride, &(ctx.unpacked_variants[1])) ||
                   bigstack_alloc_u32(uv_capacity, &(ctx.col_uvidxs[0])) ||
                   bigstack_alloc_u32(uv_capacity, &(ctx.col_uvidxs[1])) ||
                   bigstack_alloc_chridx(uv_capacity, &(ctx.row_chr_idxs[0])) ||
                   bigstack_alloc_chridx(uv_capacity, &(ctx.row_chr_idxs[1])) ||
                   bigstack_alloc_chridx(uv_capacity, &(ctx.col_chr_idxs[0])) ||
                   bigstack_alloc_chridx(uv_capacity, &(ctx.col_chr_idxs[1])) ||
                   bigstack_alloc_u32(uv_capacity, &(col_offset_starts_storage[0])) ||
                   bigstack_alloc_u32(uv_capacity, &(col_offset_starts_storage[1])) ||
                   bigstack_alloc_w(uv_capacity + 1, &(write_idx_starts_storage[0])) ||
                   bigstack_alloc_w(uv_capacity + 1, &(write_idx_starts_storage[1])))) {
        // shouldn't be possible
        goto VcorTable_ret_NOMEM;
      }
      ctx.col_offset_starts[0] = col_offset_starts_storage[0];
      ctx.col_offset_starts[1] = col_offset_starts_storage[1];
      ctx.write_idx_starts[0] = write_idx_starts_storage[0];
      ctx.write_idx_starts[1] = write_idx_starts_storage[1];
      write_ctx.col_offset_starts[0] = col_offset_starts_storage[0];
      write_ctx.col_offset_starts[1] = col_offset_starts_storage[1];
      write_ctx.write_idx_starts[0] = write_idx_starts_storage[0];
      write_ctx.write_idx_starts[1] = write_idx_starts_storage[1];
    }
    if (calc_thread_ct > uv_capacity) {
      calc_thread_ct = uv_capacity;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg) ||
                 SetThreadCt(1, &write_tg))) {
      goto VcorTable_ret_NOMEM;
    }
    SetThreadFuncAndData(VcorTableThread, &ctx, &tg);
    SetThreadFuncAndData(VcorTableWriteThread, &write_ctx, &write_tg);

    if (ref_based) {
      maj_alleles = nullptr;
    }
    AlleleCode aidx = 0;

    logprintf("Running %s with the following filter%s:\n", flagname, inter_chr? "" : "s");
    if (!inter_chr) {
      if (var_ct_radius < 0x7fffffff) {
        logprintf("  --ld-window: %u\n", var_ct_radius + 1);
      }
      logprintf("  --ld-window-kb: %g\n", 0.001 * u31tod(bp_radius));
      if (cm_radius != -1.0) {
        logprintf("  --ld-window-cm: %g\n", cm_radius);
      }
    }
    logprintf("  --ld-window-r2: %g\n", min_r2);
    if (row_snp_subset) {
      if (ld_snp_list_fname) {
        logputs("  --ld-snp-list\n");
      } else {
        logputs("  --ld-snp[s]\n");
      }
    }
    uint32_t next_print_variant_ridx = (row_variant_ct + 99) / 100;
    uint32_t pct = 0;
    uint32_t prev_variant_ridx_start = 0;
    uint32_t cur_variant_ridx_start = 0;
    uint32_t row_chr_fo_idx = UINT32_MAX;  // deliberate overflow
    uint32_t row_chr_idx = 0;
    uint32_t row_chr_start = 0;
    uint32_t row_chr_end = 0;
    uint32_t row_read_phase = 0;
    uint32_t row_parity = 0;
    uint32_t col_parity = 0;
    uint32_t col_variant_widx = 0;
    uintptr_t row_variant_uidx_base;
    uintptr_t row_cur_bits;
    BitIter1Start(row_variant_include, 0, &row_variant_uidx_base, &row_cur_bits);
    printf("%s: 0%%", flagname);
    fflush(stdout);
    do {
      // Since we're using a double- rather than a triple-buffer, we must wait
      // for the writer to finish flushing results for block (n-2) before we
      // overwrite with information about block n.
      if (prev_variant_ridx_start) {
        JoinThreads(&write_tg);
        if (unlikely(write_ctx.reterr)) {
          reterr = write_ctx.reterr;
          goto VcorTable_ret_1;
        }
      }
      // 1. Determine col_offset_start and col_offset_end for each row-variant
      //    in the new row-window, and what the current row_window_size is.
      //    Constraints:
      //    * Can't run out of results-capacity.
      //    * Unless variant_ct == uv_capacity, there's also:
      //      * row_window_size <= row_variant_ct - cur_variant_ridx_start.
      //      * row_window_size <= uv_capacity / 2.  This is unnecessarily
      //        tight in many cases, but never by a factor of more than 2 (4 if
      //        we take results-space into account), and I did not find that to
      //        be a big deal in my testing.
      //    Fill row_parity-indexed buffers.
      // 2. col_capacity := uv_capacity - row_window_size; iterate through
      //    column-shard(s).
      // uint32_t* row_uvidxs = ctx.row_uvidxs[row_parity];
      ChrIdx* row_chr_idxs = ctx.row_chr_idxs[row_parity];
      uint32_t* col_offset_starts = col_offset_starts_storage[row_parity];
      uintptr_t* write_idx_starts = write_idx_starts_storage[row_parity];
      uintptr_t write_idx = 0;
      uint32_t row_variant_uidx_first = 0;
      uint32_t row_variant_uidx_last = 0;
      uint32_t col_variant_idx_start = 0;
      uint32_t row_window_size;
      uint32_t cur_variant_ridx_stop;
      uint32_t col_variant_idx_ct;
      {
        uint32_t row_offset_limit = row_variant_ct - cur_variant_ridx_start;
        if ((uv_capacity < variant_ct) && ((uv_capacity / 2) < row_offset_limit)) {
          row_offset_limit = uv_capacity / 2;
        }

        // ok for these to be less than true value
        uint32_t uidx_start = row_chr_start;
        uint32_t uidx_end = row_chr_start;

        unsigned char* row_load_iter = ctx.unpacked_variants[col_parity];
        uint32_t row_offset = 0;
        // in inter-chr case, should save col-window size
        //   if row-subset, start at variant_idx=0 instead of at diagonal
        for (; row_offset != row_offset_limit; ++row_offset, row_load_iter = &(row_load_iter[unpacked_variant_byte_stride])) {
          const uint32_t row_variant_uidx = BitIter1(row_variant_include, &row_variant_uidx_base, &row_cur_bits);
          if (row_variant_uidx >= row_chr_end) {
            do {
              ++row_chr_fo_idx;
              row_chr_end = cip->chr_fo_vidx_start[row_chr_fo_idx + 1];
            } while (row_variant_uidx >= row_chr_end);
            row_chr_idx = cip->chr_file_order[row_chr_fo_idx];
            row_chr_start = cip->chr_fo_vidx_start[row_chr_fo_idx];
            uidx_start = row_chr_start;
            uidx_end = row_chr_start;
            if (phase_type == kR2PhaseTypePresent) {
              row_read_phase = (!IsSet(cip->haploid_mask, row_chr_idx)) || (row_chr_idx == x_code);
              if (!row_read_phase) {
                pgv.phasepresent_ct = 0;
                pgv.dphase_ct = 0;
              }
            }
          }
          uint32_t cur_col_offset_start;
          uint32_t cur_col_offset_end;
          if (inter_chr) {
            if (row_offset == 0) {
              row_variant_uidx_first = row_variant_uidx;
              if (!row_snp_subset) {
                col_variant_idx_start = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, row_variant_uidx);
              }
            }
            cur_col_offset_start = row_snp_subset? 0 : row_offset;
            cur_col_offset_end = variant_ct - col_variant_idx_start;
          } else {
            UpdateVcorWindow(variant_include, variant_bps, variant_cms, row_snp_subset, var_ct_radius, bp_radius, cm_radius, row_chr_end, row_variant_uidx, &uidx_start, &uidx_end);
            if (row_offset == 0) {
              row_variant_uidx_first = row_variant_uidx;
              col_variant_idx_start = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, uidx_start);
            }
            cur_col_offset_start = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, uidx_start) - col_variant_idx_start;
            cur_col_offset_end = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, uidx_end) - col_variant_idx_start;
          }
          const uintptr_t next_write_idx = write_idx + cur_col_offset_end - cur_col_offset_start;
          if (next_write_idx > result_capacity) {
            // Need to backtrack.
            row_variant_uidx_base = RoundDownPow2(row_variant_uidx, kBitsPerWord);
            row_cur_bits = row_variant_include[row_variant_uidx / kBitsPerWord] & (-(k1LU << (row_variant_uidx % kBitsPerWord)));
            break;
          }
          row_variant_uidx_last = row_variant_uidx;
          row_chr_idxs[row_offset] = row_chr_idx;
          col_offset_starts[row_offset] = cur_col_offset_start;
          write_idx_starts[row_offset] = write_idx;
          write_idx = next_write_idx;

          if (maj_alleles) {
            aidx = maj_alleles[row_variant_uidx];
          }
          const uint32_t is_y = (row_variant_uidx < y_end) && (row_variant_uidx >= y_start);
          if (check_dosage) {
            if (row_read_phase) {
              reterr = PgrGetInv1Dp(founder_info, pssi, founder_ct, row_variant_uidx, aidx, simple_pgrp, &pgv);
            } else {
              reterr = PgrGetInv1D(founder_info, pssi, founder_ct, row_variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
            }
            if (unlikely(reterr)) {
              goto VcorTable_ret_PGR_FAIL;
            }
            if (is_y) {
              InterleavedSetMissingCleardosage(founder_female_collapsed, founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec, &pgv.dosage_ct, pgv.dosage_present, pgv.dosage_main);
            }
            LdUnpackDosage(&pgv, founder_male_collapsed, male_dosage_invmask, founder_ct, phase_type, row_load_iter);
          } else {
            if ((!phased_calc) && (!is_y)) {
              uint32_t difflist_common_geno;
              uint32_t difflist_len;
              reterr = PgrGetInv1DifflistOrGenovec(founder_info, pssi, founder_ct, max_difflist_len, row_variant_uidx, aidx, simple_pgrp, pgv.genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
              if (unlikely(reterr)) {
                goto VcorTable_ret_PGR_FAIL;
              }
              if (difflist_common_geno != UINT32_MAX) {
                if (difflist_len <= max_difflist_len) {
                  LdUnpackNondosageSparse(raregeno, difflist_sample_ids, founder_male_collapsed, founder_ct, founder_male_ct, difflist_common_geno, difflist_len, row_load_iter);
                  continue;
                }
                PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, founder_ct, difflist_len, pgv.genovec);
              }
            } else {
              if (row_read_phase) {
                reterr = PgrGetInv1P(founder_info, pssi, founder_ct, row_variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
              } else {
                reterr = PgrGetInv1(founder_info, pssi, founder_ct, row_variant_uidx, aidx, simple_pgrp, pgv.genovec);
              }
              if (unlikely(reterr)) {
                goto VcorTable_ret_PGR_FAIL;
              }
              if (is_y) {
                InterleavedSetMissing(founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec);
              }
            }
            LdUnpackNondosageDense(&pgv, founder_male_collapsed, founder_ct, phase_type, row_load_iter);
          }
        }
        row_window_size = row_offset;
        write_idx_starts[row_window_size] = write_idx;
        cur_variant_ridx_stop = cur_variant_ridx_start + row_window_size;
        write_ctx.variant_ridx_starts[row_parity] = cur_variant_ridx_start;
        row_window_sizes[row_parity] = row_window_size;
        write_ctx.col_variant_idx_starts[row_parity] = col_variant_idx_start;
        col_variant_idx_ct = col_offset_starts[row_window_size - 1] + write_idx - write_idx_starts[row_window_size - 1];
      }
      uint32_t col_window_start = 0;
      uint32_t col_offset = 0;
      uint32_t col_read_phase = 0;
      uint32_t col_chr_fo_idx;
      uint32_t col_chr_idx;
      uint32_t col_chr_end;
      uintptr_t col_variant_uidx_base;
      uintptr_t col_cur_bits;
      {
        const uint32_t col_variant_uidx_start = ExpsearchIdxToUidx(variant_include, variant_include_cumulative_popcounts, raw_variant_ctl, col_variant_idx_start, &col_variant_widx);
        col_chr_fo_idx = GetVariantChrFoIdx(cip, col_variant_uidx_start);
        col_chr_idx = cip->chr_file_order[col_chr_fo_idx];
        col_chr_end = cip->chr_fo_vidx_start[col_chr_fo_idx + 1];
        if (phase_type == kR2PhaseTypePresent) {
          col_read_phase = (!IsSet(cip->haploid_mask, col_chr_idx)) || (col_chr_idx == x_code);
          if (!col_read_phase) {
            pgv.phasepresent_ct = 0;
            pgv.dphase_ct = 0;
          }
        }
        BitIter1Start(variant_include, col_variant_uidx_start, &col_variant_uidx_base, &col_cur_bits);
      }
      do {
        unsigned char* col_load_iter;
        {
          unsigned char* unpacked_variants = ctx.unpacked_variants[col_parity];
          const uintptr_t row_unpacked_variants_byte_ct = row_window_size * unpacked_variant_byte_stride;
          if (col_offset != 0) {
            memcpy(unpacked_variants, ctx.unpacked_variants[1 - col_parity], row_unpacked_variants_byte_ct);
          }
          col_load_iter = &(unpacked_variants[row_unpacked_variants_byte_ct]);
        }
        const uint32_t col_window_limit = MINV(col_window_start + uv_capacity, col_variant_idx_ct);
        uint32_t* col_uvidxs = ctx.col_uvidxs[col_parity];
        ChrIdx* col_chr_idxs = ctx.col_chr_idxs[col_parity];
        uint32_t col_uvidx = row_window_size;
        for (; col_offset != col_window_limit; ++col_offset) {
          const uint32_t col_variant_uidx = BitIter1(variant_include, &col_variant_uidx_base, &col_cur_bits);
          if (col_variant_uidx >= col_chr_end) {
            do {
              ++col_chr_fo_idx;
              col_chr_end = cip->chr_fo_vidx_start[col_chr_fo_idx + 1];
            } while (col_variant_uidx >= col_chr_end);
            col_chr_idx = cip->chr_file_order[col_chr_fo_idx];
            if (phase_type == kR2PhaseTypePresent) {
              col_read_phase = (!IsSet(cip->haploid_mask, col_chr_idx)) || (col_chr_idx == x_code);
              if (!col_read_phase) {
                pgv.phasepresent_ct = 0;
                pgv.dphase_ct = 0;
              }
            }
          }
          if (IsSet(row_variant_include, col_variant_uidx) && (col_variant_uidx >= row_variant_uidx_first) && (col_variant_uidx <= row_variant_uidx_last)) {
            // bugfix (18 Apr 2024): must subtract col_window_start
            col_uvidxs[col_offset - col_window_start] = RawToSubsettedPos(row_variant_include, row_variant_include_cumulative_popcounts, col_variant_uidx) - cur_variant_ridx_start;
            col_chr_idxs[col_offset - col_window_start] = col_chr_idx;
            continue;
          }
          if (col_uvidx == uv_capacity) {
            col_variant_uidx_base = RoundDownPow2(col_variant_uidx, kBitsPerWord);
            col_cur_bits = variant_include[col_variant_uidx / kBitsPerWord] & (-(k1LU << (col_variant_uidx % kBitsPerWord)));
            break;
          }
          col_uvidxs[col_offset - col_window_start] = col_uvidx++;
          col_chr_idxs[col_offset - col_window_start] = col_chr_idx;
          if (maj_alleles) {
            aidx = maj_alleles[col_variant_uidx];
          }
          const uint32_t is_y = (col_variant_uidx < y_end) && (col_variant_uidx >= y_start);
          if (check_dosage) {
            if (col_read_phase) {
              reterr = PgrGetInv1Dp(founder_info, pssi, founder_ct, col_variant_uidx, aidx, simple_pgrp, &pgv);
            } else {
              reterr = PgrGetInv1D(founder_info, pssi, founder_ct, col_variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
            }
            if (unlikely(reterr)) {
              goto VcorTable_ret_PGR_FAIL;
            }
            if (is_y) {
              InterleavedSetMissingCleardosage(founder_female_collapsed, founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec, &pgv.dosage_ct, pgv.dosage_present, pgv.dosage_main);
            }
            LdUnpackDosage(&pgv, founder_male_collapsed, male_dosage_invmask, founder_ct, phase_type, col_load_iter);
          } else {
            if ((!phased_calc) && (!is_y)) {
              uint32_t difflist_common_geno;
              uint32_t difflist_len;
              reterr = PgrGetInv1DifflistOrGenovec(founder_info, pssi, founder_ct, max_difflist_len, col_variant_uidx, aidx, simple_pgrp, pgv.genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
              if (unlikely(reterr)) {
                goto VcorTable_ret_PGR_FAIL;
              }
              if (difflist_common_geno != UINT32_MAX) {
                if (difflist_len <= max_difflist_len) {
                  LdUnpackNondosageSparse(raregeno, difflist_sample_ids, founder_male_collapsed, founder_ct, founder_male_ct, difflist_common_geno, difflist_len, col_load_iter);
                  goto VcorTable_col_finish;
                }
                PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, founder_ct, difflist_len, pgv.genovec);
              }
            } else {
              if (col_read_phase) {
                reterr = PgrGetInv1P(founder_info, pssi, founder_ct, col_variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
              } else {
                reterr = PgrGetInv1(founder_info, pssi, founder_ct, col_variant_uidx, aidx, simple_pgrp, pgv.genovec);
              }
              if (unlikely(reterr)) {
                goto VcorTable_ret_PGR_FAIL;
              }
              if (is_y) {
                InterleavedSetMissing(founder_female_collapsed_interleaved, founder_ctv2, pgv.genovec);
              }
            }
            LdUnpackNondosageDense(&pgv, founder_male_collapsed, founder_ct, phase_type, col_load_iter);
          }
        VcorTable_col_finish:
          col_load_iter = &(col_load_iter[unpacked_variant_byte_stride]);
        }
        ctx.col_window_starts[col_parity] = col_window_start;
        ctx.col_window_ends[col_parity] = col_offset;

        if (cur_variant_ridx_start || col_window_start) {
          JoinThreads(&tg);
          if (!col_window_start) {
            if (unlikely(SpawnThreads(&write_tg))) {
              goto VcorTable_ret_THREAD_CREATE_FAIL;
            }
          }
        }
        if ((cur_variant_ridx_stop == row_variant_ct) && (col_offset == col_variant_idx_ct)) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto VcorTable_ret_THREAD_CREATE_FAIL;
        }
        col_window_start = col_offset;
        col_parity = 1 - col_parity;
      } while (col_offset != col_variant_idx_ct);
      if (cur_variant_ridx_start >= next_print_variant_ridx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (cur_variant_ridx_start * 100LLU) / row_variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_ridx = (pct * S_CAST(uint64_t, row_variant_ct) + 99) / 100;
      }
      prev_variant_ridx_start = cur_variant_ridx_start;
      cur_variant_ridx_start = cur_variant_ridx_stop;
      row_parity = 1 - row_parity;
    } while (cur_variant_ridx_start != row_variant_ct);
    JoinThreads(&tg);
    if (prev_variant_ridx_start) {
      JoinThreads(&write_tg);
      if (unlikely(write_ctx.reterr)) {
        reterr = write_ctx.reterr;
        goto VcorTable_ret_1;
      }
    }
    DeclareLastThreadBlock(&write_tg);
    if (unlikely(SpawnThreads(&write_tg))) {
      goto VcorTable_ret_THREAD_CREATE_FAIL;
    }
    JoinThreads(&write_tg);
    if (unlikely(write_ctx.reterr)) {
      reterr = write_ctx.reterr;
      goto VcorTable_ret_1;
    }
    fputs("\r", stdout);
    logprintfww("%s: Results written to %s .\n", flagname, outname);
  }
  while (0) {
  VcorTable_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  VcorTable_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  VcorTable_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  VcorTable_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  }
 VcorTable_ret_1:
  CleanupThreads(&write_tg);
  CleanupThreads(&tg);
  CswriteCloseCond(&write_ctx.css, write_ctx.cswritep);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

PglErr Vcor(const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const double* variant_cms, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const VcorInfo* vcip, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  const VcorFlags flags = vcip->flags;
  const uint32_t phased_calc = (flags / kfVcorPhased) & 1;
  const uint32_t is_unsquared = (flags / kfVcorUnsquared) & 1;
  const char* flagname;
  if (phased_calc) {
    flagname = is_unsquared? "--r-phased" : "--r2-phased";
  } else {
    flagname = is_unsquared? "--r-unphased" : "--r2-unphased";
  }
  if (unlikely(founder_ct < 2)) {
    logerrprintfww("Error: %s requires at least two founders. (--make-founders may come in handy here.)\n", flagname);
    return kPglRetInconsistentInput;
  } else if (founder_ct > 0x3fffffff) {
    logerrprintf("Error: %s does not support >= 2^30 founders.\n", flagname);
    return kPglRetNotYetSupported;
  }

  const uint32_t is_matrix = ((flags & (kfVcorBin8 | kfVcorBin4 | kfVcorMatrixShapemask)) != 0);
  PglErr reterr;
  if (is_matrix) {
    reterr = VcorMatrix(orig_variant_include, cip, variant_ids, maj_alleles, founder_info, sex_nm, sex_male, vcip, flagname, raw_variant_ct, orig_variant_ct, raw_sample_ct, founder_ct, parallel_idx, parallel_tot, max_thread_ct, simple_pgrp, outname, outname_end);
  } else {
    reterr = VcorTable(orig_variant_include, cip, variant_bps, variant_ids, variant_cms, allele_idx_offsets, allele_storage, maj_alleles, allele_freqs, founder_info, sex_nm, sex_male, vcip, flagname, raw_variant_ct, orig_variant_ct, raw_sample_ct, founder_ct, max_variant_id_slen, max_allele_slen, parallel_idx, parallel_tot, max_thread_ct, simple_pgrp, outname, outname_end);
  }
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
