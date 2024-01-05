// This file is part of PLINK 2.00, copyright (C) 2005-2024 Shaun Purcell,
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


#include "include/plink2_stats.h"
#include "plink2_filter.h"
#include "plink2_ld.h"
#include "plink2_compress_stream.h"

#include <unistd.h>  // unlink()

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
  clump_ip->flags = kfClump0;
}

void CleanupClump(ClumpInfo* clump_ip) {
  free_cond(clump_ip->fnames_flattened);
  free_cond(clump_ip->test_name);
  free_cond(clump_ip->id_field);
  free_cond(clump_ip->a1_field);
  free_cond(clump_ip->test_field);
  free_cond(clump_ip->p_field);
  free_cond(clump_ip->ln_bin_boundaries);
}

// Move to plink2_common if any users outside plink2_ld.
const uintptr_t* StripUnplaced(const uintptr_t* orig_variant_include, const ChrInfo* cip, uint32_t raw_variant_ct, uint32_t* skipped_variant_ctp) {
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
    return orig_variant_include;
  }
  const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
  uintptr_t* new_variant_include;
  if (unlikely(bigstack_alloc_w(raw_variant_ctl, &new_variant_include))) {
    return nullptr;
  }
  memcpy(new_variant_include, orig_variant_include, raw_variant_ctl * sizeof(intptr_t));
  if (IsSet(cip->chr_mask, 0)) {
    const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[0];
    const uint32_t start_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
    ClearBitsNz(start_uidx, cip->chr_fo_vidx_start[chr_fo_idx + 1], new_variant_include);
  }
  if (cip->zero_extra_chrs) {
    for (uint32_t chr_idx = cip->max_code + 1; chr_idx != chr_code_end; ++chr_idx) {
      const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
      const uint32_t start_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
      ClearBitsNz(start_uidx, cip->chr_fo_vidx_start[chr_fo_idx + 1], new_variant_include);
    }
  }
  return new_variant_include;
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

void IndepPairwiseUpdateSubcontig(uint32_t variant_uidx_winstart, uint32_t x_start, uint32_t x_len, uint32_t y_start, uint32_t y_len, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t* is_x_ptr, uint32_t* is_y_ptr, uint32_t* cur_founder_ct_ptr, uint32_t* cur_founder_ctaw_ptr, uint32_t* cur_founder_ctl_ptr, uintptr_t* entire_variant_buf_word_ct_ptr) {
  // _len is better than _end here since we can exploit unsignedness
  const uint32_t is_x = ((variant_uidx_winstart - x_start) < x_len);
  const uint32_t is_y = ((variant_uidx_winstart - y_start) < y_len);
  if ((is_x != (*is_x_ptr)) || (is_y != (*is_y_ptr))) {
    *is_x_ptr = is_x;
    *is_y_ptr = is_y;
    const uint32_t cur_founder_ct = (is_x || is_y)? founder_male_ct : founder_ct;
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
        IndepPairwiseUpdateSubcontig(variant_uidx_winstart, x_start, x_len, y_start, y_len, founder_ct, founder_male_ct, &is_x, &is_y, &cur_founder_ct, &cur_founder_ctaw, &cur_founder_ctl, &entire_variant_buf_word_ct);
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

PglErr IndepPairwise(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uint32_t* founder_info_cumulative_popcounts, const uintptr_t* founder_nonmale, const uintptr_t* founder_male, const LdInfo* ldip, const uintptr_t* preferred_variants, const uint32_t* subcontig_info, const uint32_t* subcontig_thread_assignments, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t subcontig_ct, uintptr_t window_max, uint32_t calc_thread_ct, uint32_t max_load, PgenReader* simple_pgrp, uintptr_t* removed_variants_collapsed) {
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
    uint32_t next_print_tvidx_start = max_load / 100;
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
              if (is_haploid) {
                SetHetMissing(founder_ctl2, cur_raw_tgenovec);
              }
            } else {
              reterr = PgrGetInv1(nullptr, pssi, raw_sample_ct, variant_uidx, maj_alleles[variant_uidx], simple_pgrp, tmp_genovec);
              if (founder_male_ct) {
                CopyNyparrNonemptySubset(tmp_genovec, founder_male, raw_sample_ct, founder_male_ct, cur_raw_tgenovec);
                SetHetMissing(founder_male_ctl2, cur_raw_tgenovec);
              }
              if (is_x) {
                CopyNyparrNonemptySubset(tmp_genovec, founder_nonmale, raw_sample_ct, founder_nonmale_ct, &(cur_raw_tgenovec[founder_male_ctl2]));
                if (all_haploid) {
                  // don't just treat chrX identically to autosomes, since for
                  // doubled haploids we still want to give females 2x the
                  // weight of males.  I think.
                  SetHetMissing(founder_nonmale_ctl2, &(cur_raw_tgenovec[founder_male_ctl2]));
                }
              }
            }
            if (unlikely(reterr)) {
              PgenErrPrintNV(reterr, variant_uidx);
              goto IndepPairwise_ret_1;
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
          next_print_tvidx_start = (pct * S_CAST(uint64_t, max_load)) / 100;
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

void IndepPairphaseUpdateSubcontig(const ChrInfo* cip, uint32_t variant_uidx_winstart, uint32_t x_start, uint32_t x_len, uint32_t y_start, uint32_t y_len, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t* is_x_ptr, uint32_t* is_y_ptr, uint32_t* is_haploid_ptr, uint32_t* cur_hap_ct_ptr, uint32_t* cur_hap_ctl_ptr) {
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
      cur_hap_ct = founder_male_ct;
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
        IndepPairphaseUpdateSubcontig(cip, variant_uidx_winstart, x_start, x_len, y_start, y_len, founder_ct, founder_male_ct, &is_x, &is_y, &is_haploid, &cur_hap_ct, &cur_hap_ctl);
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

PglErr IndepPairphase(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uint32_t* founder_info_cumulative_popcounts, const uintptr_t* founder_nonmale, const uintptr_t* founder_male, const LdInfo* ldip, const uintptr_t* preferred_variants, const uint32_t* subcontig_info, const uint32_t* subcontig_thread_assignments, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t subcontig_ct, uintptr_t window_max, uint32_t calc_thread_ct, uint32_t max_load, PgenReader* simple_pgrp, uintptr_t* removed_variants_collapsed) {
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
    uintptr_t* phasepresent;
    uintptr_t* phaseinfo;
    uintptr_t* raw_genovec;
    uintptr_t* raw_phasepresent;
    uintptr_t* raw_phaseinfo;
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
    uint32_t next_print_tvidx_start = max_load / 100;
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
              if (founder_male_ct) {
                CopyNyparrNonemptySubset(raw_genovec, founder_male, raw_sample_ct, founder_male_ct, genovec);
                HapsplitHaploid(genovec, founder_male_ct, cur_loader_hap_vec, cur_loader_nm_vec);
              }
              if (is_x) {
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
          next_print_tvidx_start = (pct * S_CAST(uint64_t, max_load)) / 100;
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

PglErr LdPrune(const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uintptr_t* sex_male, const LdInfo* ldip, const char* indep_preferred_fname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
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
    const uintptr_t* variant_include = StripUnplaced(orig_variant_include, cip, raw_variant_ct, &skipped_variant_ct);
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
    } else {
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, &preferred_variants))) {
        goto LdPrune_ret_NOMEM;
      }
      memcpy(preferred_variants, variant_include, raw_variant_ctl * sizeof(intptr_t));
      reterr = NondupIdLoad(variant_ids, indep_preferred_fname, raw_variant_ct, variant_ct, max_thread_ct, preferred_variants, &dup_found);
    }
    if (unlikely(reterr)) {
      goto LdPrune_ret_1;
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
    const uint32_t founder_male_ct = PopcountWordsIntersect(founder_info, sex_male, raw_sample_ctl);
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
      reterr = IndepPairphase(variant_include, cip, variant_bps, allele_idx_offsets, maj_alleles, allele_freqs, founder_info, founder_info_cumulative_popcounts, founder_nonmale, founder_male, ldip, preferred_variants, subcontig_info, subcontig_thread_assignments, raw_sample_ct, founder_ct, founder_male_ct, subcontig_ct, window_max, max_thread_ct, max_load, simple_pgrp, removed_variants_collapsed);
    } else {
      reterr = IndepPairwise(variant_include, cip, variant_bps, allele_idx_offsets, maj_alleles, allele_freqs, founder_info, founder_info_cumulative_popcounts, founder_nonmale, founder_male, ldip, preferred_variants, subcontig_info, subcontig_thread_assignments, raw_sample_ct, founder_ct, founder_male_ct, subcontig_ct, window_max, max_thread_ct, max_load, simple_pgrp, removed_variants_collapsed);
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

  const uint32_t sample_ct_rem = sample_ct % kBitsPerWordD2;
  if (sample_ct_rem) {
    const Halfword trailing_mask = (1U << sample_ct_rem) - 1;
    one_bitarr_alias[sample_ctl2 - 1] &= trailing_mask;
    two_bitarr_alias[sample_ctl2 - 1] &= trailing_mask;
    nm_bitarr_alias[sample_ctl2 - 1] &= trailing_mask;
  }
  if (sample_ctl2 % 2) {
    one_bitarr_alias[sample_ctl2] = 0;
    two_bitarr_alias[sample_ctl2] = 0;
    nm_bitarr_alias[sample_ctl2] = 0;
  }
}

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
  //                  2 * popcount(subset_mask & two_bitvec)
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

// phased-dosage r^2 computation:
//   just do brute force for now
//   use dense dosage_sum/(optional dphase_delta) representation
//     when either dphase_delta is null, can skip that dot product
//   also have a unphased_het_dosage array pointer.  this is null if all
//     dosages are phased.  otherwise...
//     suppose one unphased sample has dosage(var0)=0.2 and dosage(var1)=1.4.
//     this is stored as strandA0[] = strandB0[] = 0.1,
//                       strandA1[] = strandB1[] = 0.7,
//     so the sum of the two products is 0.14.
//     we treat this as P(var0=0/0)=0.8, P(var0=0/1)=0.2,
//                      P(var1=0/1)=0.6, P(var1=1/1)=0.4,
//     so the diplotype dosages are
//       0-0: 2 * 0.9 * 0.3 = 0.54
//       0-1: 2 * 0.9 * 0.7 = 1.26
//       1-0: 2 * 0.1 * 0.3 = 0.06
//       1-1: 2 * 0.1 * 0.7 = 0.14
//     the no-phasing-error components of this are:
//       var0=0/0, var1=0/1:
//         0-0: 0.8 * 0.6 = 0.48
//         0-1:           = 0.48
//       var0=0/0, var1=1/1:
//         0-1: 2 * 0.8 * 0.4 = 0.64
//       var0=0/1, var1=1/1:
//         0-1: 0.2 * 0.4 = 0.08
//         1-1:           = 0.08
//     the uncertain-phasing component of this is 2 * 0.2 * 0.6 = 0.24; a
//       quarter of this contributes to the sum-of-products.
//     if we save P(var=0/1) = (1 - abs(1 - dosagesum)) in unphased_het_dosage
//       for each unphased dosage (and 0 for each phased dosage), subtracting
//       half of the unphased_het_dosage dot product (0.12/2 = 0.06 in this
//       case) from the main dot product yields the definitely-known portion.
//       the unhalved unphased_het_dosage dot product is the maximum possible
//       value of the unknown portion (half_hethet_share).

static_assert(sizeof(Dosage) == 2, "plink2_ld dosage-handling routines must be updated.");
#ifdef __LP64__
#  ifdef USE_AVX2
void FillDosageUhet(const Dosage* dosage_vec, uint32_t dosagev_ct, Dosage* dosage_uhet) {
  const __m256i* dosage_vvec_iter = R_CAST(const __m256i*, dosage_vec);
#    if defined(__APPLE__) && ((!defined(__cplusplus)) || (__cplusplus < 201103L))
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
  __m256i* dosage_uhet_iter = R_CAST(__m256i*, dosage_uhet);
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
    *dosage_uhet_iter++ = _mm256_add_epi16(dosagev, dosagev_opp);
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

uint64_t DosageUnsignedNomissDotprod(const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const __m256i* dosage_vvec0_iter = R_CAST(const __m256i*, dosage_vec0);
  const __m256i* dosage_vvec1_iter = R_CAST(const __m256i*, dosage_vec1);
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
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
      __m256i dosage0 = *dosage_vvec0_iter++;
      __m256i dosage1 = *dosage_vvec1_iter++;

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

int64_t DosageSignedDotprod(const SDosage* dphase_delta0, const SDosage* dphase_delta1, uint32_t vec_ct) {
  const __m256i* dphase_delta0_iter = R_CAST(const __m256i*, dphase_delta0);
  const __m256i* dphase_delta1_iter = R_CAST(const __m256i*, dphase_delta1);
  const __m256i m16 = _mm256_set1_epi64x(kMask0000FFFF);
  const __m256i all_4096 = _mm256_set1_epi16(0x1000);
  uint64_t dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m256i dotprod_lo = _mm256_setzero_si256();
    __m256i dotprod_hi = _mm256_setzero_si256();
    const __m256i* dphase_delta0_stop;
    if (vecs_left < 4096) {
      if (!vecs_left) {
        // this cancels out the shift-hi16-by-4096 below
        return S_CAST(int64_t, dotprod) - (0x10000000LLU * kDosagePerVec) * vec_ct;
      }
      dphase_delta0_stop = &(dphase_delta0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dphase_delta0_stop = &(dphase_delta0_iter[4096]);
      vecs_left -= 4096;
    }
    do {
      __m256i dosage0 = *dphase_delta0_iter++;
      __m256i dosage1 = *dphase_delta1_iter++;

      __m256i hi16 = _mm256_mulhi_epi16(dosage0, dosage1);
      __m256i lo16 = _mm256_mullo_epi16(dosage0, dosage1);
      // original values are in [-16384, 16384]
      // product is in [-2^28, 2^28], so hi16 is in [-4096, 4096]
      // so if we add 4096 to hi16, we can treat it as an unsigned value in the
      //   rest of this loop
      // todo: try rewriting all these dot products to use _mm256_madd_epi16(),
      // pretty sure that's substantially better
      // however, that still triggers universal ~15% clock frequency slowdown,
      // so also try sparse strategy on real imputed data and be prepared to
      // throw out this entire function
      hi16 = _mm256_add_epi16(hi16, all_4096);
      lo16 = _mm256_add_epi64(_mm256_and_si256(lo16, m16), _mm256_and_si256(_mm256_srli_epi64(lo16, 16), m16));
      hi16 = _mm256_and_si256(_mm256_add_epi64(hi16, _mm256_srli_epi64(hi16, 16)), m16);
      dotprod_lo = _mm256_add_epi64(dotprod_lo, lo16);
      dotprod_hi = _mm256_add_epi64(dotprod_hi, hi16);
    } while (dphase_delta0_iter < dphase_delta0_stop);
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
#  else  // !USE_AVX2
void FillDosageUhet(const Dosage* dosage_vec, uint32_t dosagev_ct, Dosage* dosage_uhet) {
  const __m128i* dosage_vvec_iter = R_CAST(const __m128i*, dosage_vec);
#    if defined(__APPLE__) && ((!defined(__cplusplus)) || (__cplusplus < 201103L))
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
  __m128i* dosage_uhet_iter = R_CAST(__m128i*, dosage_uhet);
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
    *dosage_uhet_iter++ = _mm_add_epi16(dosagev, dosagev_opp);
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

uint64_t DosageUnsignedNomissDotprod(const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const __m128i* dosage_vvec0_iter = R_CAST(const __m128i*, dosage_vec0);
  const __m128i* dosage_vvec1_iter = R_CAST(const __m128i*, dosage_vec1);
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
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

int64_t DosageSignedDotprod(const SDosage* dphase_delta0, const SDosage* dphase_delta1, uint32_t vec_ct) {
  const __m128i* dphase_delta0_iter = R_CAST(const __m128i*, dphase_delta0);
  const __m128i* dphase_delta1_iter = R_CAST(const __m128i*, dphase_delta1);
  const __m128i m16 = _mm_set1_epi64x(kMask0000FFFF);
  const __m128i all_4096 = _mm_set1_epi16(0x1000);
  uint64_t dotprod = 0;
  for (uint32_t vecs_left = vec_ct; ; ) {
    __m128i dotprod_lo = _mm_setzero_si128();
    __m128i dotprod_hi = _mm_setzero_si128();
    const __m128i* dphase_delta0_stop;
    if (vecs_left < 8192) {
      if (!vecs_left) {
        // this cancels out the shift-hi16-by-4096 below
        return S_CAST(int64_t, dotprod) - (0x10000000LLU * kDosagePerVec) * vec_ct;
      }
      dphase_delta0_stop = &(dphase_delta0_iter[vecs_left]);
      vecs_left = 0;
    } else {
      dphase_delta0_stop = &(dphase_delta0_iter[8192]);
      vecs_left -= 8192;
    }
    do {
      __m128i dosage0 = *dphase_delta0_iter++;
      __m128i dosage1 = *dphase_delta1_iter++;

      __m128i hi16 = _mm_mulhi_epi16(dosage0, dosage1);
      __m128i lo16 = _mm_mullo_epi16(dosage0, dosage1);
      // original values are in [-16384, 16384]
      // product is in [-2^28, 2^28], so hi16 is in [-4096, 4096]
      // so if we add 4096 to hi16, we can treat it as an unsigned value in the
      //   rest of this loop
      hi16 = _mm_add_epi16(hi16, all_4096);
      lo16 = _mm_add_epi64(_mm_and_si128(lo16, m16), _mm_and_si128(_mm_srli_epi64(lo16, 16), m16));
      hi16 = _mm_and_si128(_mm_add_epi64(hi16, _mm_srli_epi64(hi16, 16)), m16);
      dotprod_lo = _mm_add_epi64(dotprod_lo, lo16);
      dotprod_hi = _mm_add_epi64(dotprod_hi, hi16);
    } while (dphase_delta0_iter < dphase_delta0_stop);
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
#  endif  // !USE_AVX2
#else  // !__LP64__
void FillDosageUhet(const Dosage* dosage_vec, uint32_t dosagev_ct, Dosage* dosage_uhet) {
  const uint32_t sample_cta2 = dosagev_ct * 2;
  for (uint32_t sample_idx = 0; sample_idx != sample_cta2; ++sample_idx) {
    const uint32_t cur_dosage = dosage_vec[sample_idx];
    uint32_t cur_hetval = cur_dosage;
    if (cur_hetval > 16384) {
      if (cur_hetval == kDosageMissing) {
        cur_hetval = 0;
      } else {
        cur_hetval = 32768 - cur_hetval;
      }
    }
    dosage_uhet[sample_idx] = cur_hetval;
  }
}

uint64_t DenseDosageSum(const Dosage* dosage_vec, uint32_t vec_ct) {
  const uint32_t sample_cta2 = vec_ct * 2;
  uint64_t sum = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_cta2; ++sample_idx) {
    const uint32_t cur_dosage = dosage_vec[sample_idx];
    if (cur_dosage != kDosageMissing) {
      sum += cur_dosage;
    }
  }
  return sum;
}

uint64_t DenseDosageSumSubset(const Dosage* dosage_vec, const Dosage* dosage_mask_vec, uint32_t vec_ct) {
  const uint32_t sample_cta2 = vec_ct * 2;
  uint64_t sum = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_cta2; ++sample_idx) {
    const uint32_t cur_dosage = dosage_vec[sample_idx];
    const uint32_t other_dosage = dosage_mask_vec[sample_idx];
    if ((cur_dosage != kDosageMissing) && (other_dosage != kDosageMissing)) {
      sum += cur_dosage;
    }
  }
  return sum;
}

uint64_t DosageUnsignedDotprod(const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const uint32_t sample_cta2 = vec_ct * 2;
  uint64_t dotprod = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_cta2; ++sample_idx) {
    const uint32_t cur_dosage0 = dosage_vec0[sample_idx];
    const uint32_t cur_dosage1 = dosage_vec1[sample_idx];
    if ((cur_dosage0 != kDosageMissing) && (cur_dosage1 != kDosageMissing)) {
      dotprod += cur_dosage0 * cur_dosage1;
    }
  }
  return dotprod;
}

uint64_t DosageUnsignedNomissDotprod(const Dosage* dosage_vec0, const Dosage* dosage_vec1, uint32_t vec_ct) {
  const uint32_t sample_cta2 = vec_ct * 2;
  uint64_t dotprod = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_cta2; ++sample_idx) {
    const uint32_t cur_dosage0 = dosage_vec0[sample_idx];
    const uint32_t cur_dosage1 = dosage_vec1[sample_idx];
    dotprod += cur_dosage0 * cur_dosage1;
  }
  return dotprod;
}

int64_t DosageSignedDotprod(const SDosage* dphase_delta0, const SDosage* dphase_delta1, uint32_t vec_ct) {
  const uint32_t sample_cta2 = vec_ct * 2;
  int64_t dotprod = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_cta2; ++sample_idx) {
    const int32_t cur_diff0 = dphase_delta0[sample_idx];
    const int32_t cur_diff1 = dphase_delta1[sample_idx];
    dotprod += cur_diff0 * cur_diff1;
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
#endif

uint32_t DosageR2Prod(const Dosage* dosage_vec0, const uintptr_t* nm_bitvec0, const Dosage* dosage_vec1, const uintptr_t* nm_bitvec1, uint32_t sample_ct, uint32_t nm_ct0, uint32_t nm_ct1, uint64_t* __restrict nmaj_dosages, uint64_t* __restrict dosageprod_ptr) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  uint32_t nm_intersection_ct;
  if ((nm_ct0 != sample_ct) && (nm_ct1 != sample_ct)) {
    nm_intersection_ct = PopcountWordsIntersect(nm_bitvec0, nm_bitvec1, sample_ctl);
    if (!nm_intersection_ct) {
      nmaj_dosages[0] = 0;
      nmaj_dosages[1] = 0;
      *dosageprod_ptr = 0;
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
  // could conditionally use dosage_unsigned_nomiss here
  *dosageprod_ptr = DosageUnsignedDotprod(dosage_vec0, dosage_vec1, vec_ct);
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
// If extra_retp is nullptr, results[0] contains r2, and if compute_dprime is
//   true, results[1] = dprime.
LDErr PhasedLD(const double* nmajsums_d, double known_dotprod_d, double unknown_hethet_d, double twice_tot_recip, uint32_t compute_dprime, PhasedLDExtraRet* extra_retp, double* results) {
  // known-diplotype dosages (sum is 2 * (valid_obs_d - unknown_hethet_d)):
  //   var0  var1
  //     0  -  0 : 2 * valid_obs_d - majsums[0] - majsums[1] + known_dotprod
  //     1  -  0 : majsums[0] - known_dotprod - unknown_hethet_d
  //     0  -  1 : majsums[1] - known_dotprod - unknown_hethet_d
  //     1  -  1 : known_dotprod

  // bugfix (15 Sep 2023): otherwise possible for freq_majmaj to be slightly
  // less than zero, and this breaks EmPhasedUnscaledLnlike().
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
        // have encountered 7.9e-11 difference in testing, which is more than
        // twice kSmallishEpsilon.
        while (cubic_sols[cubic_sol_ct - 1] > half_unphased_hethet_share + 8 * kSmallishEpsilon) {
          --cubic_sol_ct;
          if (cubic_sol_ct == 1) {
            break;
          }
        }
        if (cubic_sols[cubic_sol_ct - 1] > half_unphased_hethet_share - 8 * kSmallishEpsilon) {
          cubic_sols[cubic_sol_ct - 1] = half_unphased_hethet_share;
        }
        while ((cubic_sols[first_relevant_sol_idx] < 8 * -kSmallishEpsilon) && (first_relevant_sol_idx + 1 < cubic_sol_ct)) {
          ++first_relevant_sol_idx;
        }
        if (cubic_sols[first_relevant_sol_idx] < 8 * kSmallishEpsilon) {
          cubic_sols[first_relevant_sol_idx] = 0.0;
        }
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
      if ((nonzero_freq_xx + kSmallishEpsilon < half_unphased_hethet_share + nonzero_freq_xy) && (nonzero_freq_xy + kSmallishEpsilon < half_unphased_hethet_share + nonzero_freq_xx)) {
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
    if (compute_dprime) {
      // maybe this should just always be computed since it's such a small
      // fraction of the total cost?
      if (dd >= 0.0) {
        results[1] = dd / MINV(freq_xmaj * freq_minx, freq_xmin * freq_majx);
      } else {
        results[1] = -dd / MINV(freq_xmaj * freq_majx, freq_xmin * freq_minx);
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

PglErr LdConsole(const uintptr_t* variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const char* const* allele_storage, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const LdInfo* ldip, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, PgenReader* simple_pgrp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    if (!founder_ct) {
      logerrputs("Error: --ld requires founders.  (--make-founders may come in handy here.)\n");
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
    uint32_t is_nonx_haploids[2];
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
      uint32_t is_nonx_haploid = 0;
      if (IsSet(cip->haploid_mask, chr_idx)) {
        is_nonx_haploid = 1 - is_x;
        y_ct += (chr_idx == y_code);
      }
      is_nonx_haploids[var_idx] = is_nonx_haploid;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    // if both unplaced, don't count as same-chromosome
    const uint32_t is_same_chr = chr_idxs[0] && (chr_idxs[0] == chr_idxs[1]);
    if (y_ct) {
      // only keep male founders
      uintptr_t* founder_info_tmp;
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &founder_info_tmp))) {
        goto LdConsole_ret_NOMEM;
      }
      BitvecAndCopy(founder_info, sex_male, raw_sample_ctl, founder_info_tmp);
      founder_info = founder_info_tmp;
      founder_ct = PopcountWords(founder_info, raw_sample_ctl);
      if (!founder_ct) {
        logerrprintfww("Warning: Skipping --ld since there are no male founders, and %s specified. (--make-founders may come in handy here.)\n", is_same_chr? "chrY variants were" : "a chrY variant was");
        goto LdConsole_ret_1;
      }
    }
    const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
    const uint32_t founder_ctl2 = NypCtToWordCt(founder_ct);
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
    const uint32_t founder_ctv2 = NypCtToVecCt(founder_ct);
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
    uint32_t use_dosage = ldip->ld_console_flags & kfLdConsoleDosage;
    if (use_dosage) {
      logerrputs("Error: Alpha 5 --ld's handling of dosages is incorrect, and has been disabled.\nUse an alpha 6 or later build for this functionality.\n");
      reterr = kPglRetNotYetSupported;
      goto LdConsole_ret_1;
    }

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
      if (is_nonx_haploids[var_idx]) {
        if (!use_dosage) {
          SetHetMissing(founder_ctl2, pgvp->genovec);
        }
        pgvp->phasepresent_ct = 0;
        pgvp->dphase_ct = 0;
      } else if (x_male_ct && is_xs[var_idx]) {
        if (!use_dosage) {
          // bugfix (23 Feb 2022): need genovec vec count, not
          // sex_male_collapsed_interleaved vec count.
          SetMaleHetMissing(sex_male_collapsed_interleaved, founder_ctv2, pgvp->genovec);
        }
        if (pgvp->phasepresent_ct) {
          BitvecInvmask(sex_male_collapsed, founder_ctl, pgvp->phasepresent);
          pgvp->phasepresent_ct = PopcountWords(pgvp->phasepresent, founder_ctl);
        }
        if (pgvp->dphase_ct) {
          EraseMaleDphases(sex_male_collapsed, &(pgvp->dphase_ct), pgvp->dphase_present, pgvp->dphase_delta);
        }
      }
    }
    const uint32_t use_phase = is_same_chr && (pgvs[0].phasepresent_ct || pgvs[0].dphase_ct) && (pgvs[1].phasepresent_ct || pgvs[1].dphase_ct);
    const uint32_t ignore_hethet = is_nonx_haploids[0] || is_nonx_haploids[1];
    if ((!pgvs[0].dosage_ct) && (!pgvs[1].dosage_ct) && (!ignore_hethet)) {
      // If the 'dosage' modifier was present, no literal dosages are present,
      // but one or both variants is purely haploid, may as well use the dosage
      // code path anyway (converting hets to dosage 0.5).
      use_dosage = 0;
    }

    // values of interest:
    //   mutually-nonmissing observation count
    //   (all other values computed over mutually-nonmissing set)
    //   4 known-diplotype dosages (0..2 for each sample, in unphased het-het)
    //   (unphased het-het fractional count can be inferred)
    //   dosage sum for each variant
    double x_male_known_dotprod_d = 0.0;
    uint32_t valid_x_male_ct = 0;
    double nmajsums_d[2];
    double x_male_nmajsums_d[2];
    double known_dotprod_d;
    double unknown_hethet_d;
    uint32_t valid_obs_ct;
    uint32_t hethet_present;
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
      const uint32_t orig_maj_ct1 = nmaj_cts[1];
      uint32_t known_dotprod;
      uint32_t unknown_hethet_ct;
      valid_obs_ct = HardcallPhasedR2Stats(one_bitvecs[0], two_bitvecs[0], nm_bitvecs[0], one_bitvecs[1], two_bitvecs[1], nm_bitvecs[1], founder_ct, nm_cts[0], nm_cts[1], nmaj_cts, &known_dotprod, &unknown_hethet_ct);
      if (unlikely(!valid_obs_ct)) {
        goto LdConsole_ret_NO_VALID_OBSERVATIONS;
      }
      hethet_present = (unknown_hethet_ct != 0);
      if (use_phase && hethet_present) {
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
      if (x_male_ct) {
        // on chrX, store separate full-size copies of one_bitvec, two_bitvec,
        //   and nm_bitvec with nonmales masked out
        // (we can bitvec-and here because we're not doing any further
        // calculations.  it suffices to bitvec-and one side)
        BitvecAnd(sex_male_collapsed, founder_ctl, one_bitvecs[0]);
        BitvecAnd(sex_male_collapsed, founder_ctl, two_bitvecs[0]);
        BitvecAnd(sex_male_collapsed, founder_ctl, nm_bitvecs[0]);
        uint32_t x_male_nmaj_cts[2];
        x_male_nmaj_cts[0] = GenoBitvecSum(one_bitvecs[0], two_bitvecs[0], founder_ctl);

        x_male_nmaj_cts[1] = orig_maj_ct1;

        const uint32_t x_male_nm_ct0 = PopcountWords(nm_bitvecs[0], founder_ctl);
        uint32_t x_male_known_dotprod;
        uint32_t x_male_unknown_hethet_ct;  // ignore
        valid_x_male_ct = HardcallPhasedR2Stats(one_bitvecs[0], two_bitvecs[0], nm_bitvecs[0], one_bitvecs[1], two_bitvecs[1], nm_bitvecs[1], founder_ct, x_male_nm_ct0, nm_cts[1], x_male_nmaj_cts, &x_male_known_dotprod, &x_male_unknown_hethet_ct);
        x_male_nmajsums_d[0] = u31tod(x_male_nmaj_cts[0]);
        x_male_nmajsums_d[1] = u31tod(x_male_nmaj_cts[1]);
        x_male_known_dotprod_d = S_CAST(double, x_male_known_dotprod);
        // hethet impossible for chrX males
        assert(!x_male_unknown_hethet_ct);
      }
    } else {
      // Current brute-force strategy:
      // 1. Expand each variant to all-Dosage format, with an optional
      //    phased dosage signed-difference track.
      // 2. Given (a0+b0), (a0-b0), (a1+b1), and (a1-b1)
      //    We wish to compute a0*a1 + b0*b1
      //      (a0+b0) * (a1+b1) = a0*a1 + b0*b1 + a0*b1 + a1*b0
      //      (a0-b0) * (a1-b1) = a0*a1 + b0*b1 - a0*b1 - a1*b0
      //      so halving the sum of these two dot products works.
      const uint32_t founder_dosagev_ct = DivUp(founder_ct, kDosagePerVec);
      Dosage* dosage_vecs[2];
      Dosage* dosage_uhets[2];
      uintptr_t* nm_bitvecs[2];
      // founder_ct automatically rounded up as necessary
      if (unlikely(bigstack_alloc_dosage(founder_ct, &dosage_vecs[0]) ||
                   bigstack_alloc_dosage(founder_ct, &dosage_vecs[1]) ||
                   bigstack_alloc_dosage(founder_ct, &dosage_uhets[0]) ||
                   bigstack_alloc_dosage(founder_ct, &dosage_uhets[1]) ||
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
        FillDosageUhet(dosage_vecs[var_idx], founder_dosagev_ct, dosage_uhets[var_idx]);
        GenoarrToNonmissingnessUnsafe(pgvs[var_idx].genovec, founder_ct, nm_bitvecs[var_idx]);
        ZeroTrailingBits(founder_ct, nm_bitvecs[var_idx]);
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
        assert(0);
        if (unlikely(bigstack_alloc_dphase(founder_ct, &main_dphase_deltas[0]) ||
                     bigstack_alloc_dphase(founder_ct, &main_dphase_deltas[1]))) {
          goto LdConsole_ret_NOMEM;
        }
        for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
          PopulateDenseDphase(pgvs[var_idx].phasepresent, pgvs[var_idx].phaseinfo, pgvs[var_idx].dosage_present, dosage_vecs[var_idx], pgvs[var_idx].dphase_present, pgvs[var_idx].dphase_delta, founder_ct, pgvs[var_idx].phasepresent_ct, pgvs[var_idx].dosage_ct, pgvs[var_idx].dphase_ct, dosage_uhets[var_idx], main_dphase_deltas[var_idx]);
        }
      }
      const uint64_t orig_nmaj_dosage1 = nmaj_dosages[1];
      uint64_t dosageprod;
      valid_obs_ct = DosageR2Prod(dosage_vecs[0], nm_bitvecs[0], dosage_vecs[1], nm_bitvecs[1], founder_ct, nm_cts[0], nm_cts[1], nmaj_dosages, &dosageprod);
      if (unlikely(!valid_obs_ct)) {
        goto LdConsole_ret_NO_VALID_OBSERVATIONS;
      }
      uint64_t uhethet_dosageprod = 0;
      if (!ignore_hethet) {
        uhethet_dosageprod = DosageUnsignedNomissDotprod(dosage_uhets[0], dosage_uhets[1], founder_dosagev_ct);
      }
      hethet_present = (uhethet_dosageprod != 0);
      if (use_phase && hethet_present) {
        dosageprod = S_CAST(int64_t, dosageprod) + DosageSignedDotprod(main_dphase_deltas[0], main_dphase_deltas[1], founder_dosagev_ct);
      }
      nmajsums_d[0] = u63tod(nmaj_dosages[0]) * kRecipDosageMid;
      nmajsums_d[1] = u63tod(nmaj_dosages[1]) * kRecipDosageMid;
      known_dotprod_d = u63tod(dosageprod - uhethet_dosageprod) * (kRecipDosageMidSq * 0.5);
      unknown_hethet_d = u63tod(uhethet_dosageprod) * kRecipDosageMidSq;
      if (x_male_ct) {
        Dosage* x_male_dosage_invmask;
        if (unlikely(bigstack_alloc_dosage(founder_ct, &x_male_dosage_invmask))) {
          goto LdConsole_ret_NOMEM;
        }
        SetAllDosageArr(founder_dosagev_ct * kDosagePerVec, x_male_dosage_invmask);
        uintptr_t sample_midx_base = 0;
        uintptr_t cur_bits = sex_male_collapsed[0];
        for (uint32_t uii = 0; uii != x_male_ct; ++uii) {
          const uintptr_t sample_midx = BitIter1(sex_male_collapsed, &sample_midx_base, &cur_bits);
          x_male_dosage_invmask[sample_midx] = 0;
        }
        BitvecOr(R_CAST(uintptr_t*, x_male_dosage_invmask), founder_dosagev_ct * kWordsPerVec, R_CAST(uintptr_t*, dosage_vecs[0]));
        BitvecAnd(sex_male_collapsed, founder_ctl, nm_bitvecs[0]);
        uint64_t x_male_nmaj_dosages[2];
        x_male_nmaj_dosages[0] = DenseDosageSum(dosage_vecs[0], founder_dosagev_ct);
        x_male_nmaj_dosages[1] = orig_nmaj_dosage1;
        const uint32_t x_male_nm_ct0 = PopcountWords(nm_bitvecs[0], founder_ctl);
        uint64_t x_male_dosageprod;
        valid_x_male_ct = DosageR2Prod(dosage_vecs[0], nm_bitvecs[0], dosage_vecs[1], nm_bitvecs[1], founder_ct, x_male_nm_ct0, nm_cts[1], x_male_nmaj_dosages, &x_male_dosageprod);
        if (!ignore_hethet) {
          BitvecInvmask(R_CAST(uintptr_t*, x_male_dosage_invmask), founder_dosagev_ct * kWordsPerVec, R_CAST(uintptr_t*, dosage_uhets[0]));
          const uint64_t invalid_uhethet_dosageprod = DosageUnsignedNomissDotprod(dosage_uhets[0], dosage_uhets[1], founder_dosagev_ct);
          unknown_hethet_d -= u63tod(invalid_uhethet_dosageprod) * kRecipDosageMidSq;
          known_dotprod_d += u63tod(invalid_uhethet_dosageprod) * (kRecipDosageMidSq * 0.5);
        }
        x_male_nmajsums_d[0] = u63tod(x_male_nmaj_dosages[0]) * kRecipDosageMid;
        x_male_nmajsums_d[1] = u63tod(x_male_nmaj_dosages[1]) * kRecipDosageMid;
        x_male_known_dotprod_d = u63tod(x_male_dosageprod) * (kRecipDosageMidSq * 0.5);
      }
    }
    double valid_obs_d = u31tod(valid_obs_ct);
    if (valid_x_male_ct) {
      // males have sqrt(0.5) weight if one variant is chrX, half-weight if
      // both are chrX
      const double male_decr = (is_xs[0] && is_xs[1])? 0.5 : (1.0 - 0.5 * kSqrt2);
      nmajsums_d[0] -= male_decr * x_male_nmajsums_d[0];
      nmajsums_d[1] -= male_decr * x_male_nmajsums_d[1];
      known_dotprod_d -= male_decr * x_male_known_dotprod_d;
      valid_obs_d -= male_decr * u31tod(valid_x_male_ct);
    }

    const double twice_tot_recip = 0.5 / valid_obs_d;
    // in plink 1.9, "freq12" refers to first variant=1, second variant=2
    // this most closely corresponds to freq_majmin here
    PhasedLDExtraRet extra_ret;
    double cubic_sols[3];
    const LDErr ld_err = PhasedLD(nmajsums_d, known_dotprod_d, unknown_hethet_d, twice_tot_recip, 1, &extra_ret, cubic_sols);
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
      write_iter = strcpya_k(write_iter, " male");
    }
    write_iter = strcpya_k(write_iter, " sample");
    if (valid_obs_ct != 1) {
      *write_iter++ = 's';
    }
    if (valid_x_male_ct && (!y_ct)) {
      write_iter = strcpya_k(write_iter, " (");
      write_iter = u32toa(valid_x_male_ct, write_iter);
      write_iter = strcpya_k(write_iter, " male)");
    }
    if ((!is_nonx_haploids[0]) && (!is_nonx_haploids[1])) {
      write_iter = strcpya_k(write_iter, "; ");
      if (unknown_hethet_d == 0.0) {
        if (hethet_present) {
          write_iter = strcpya_k(write_iter, "all phased");
        } else {
          write_iter = strcpya_k(write_iter, "no het pairs present");
        }
      } else {
        // ddosagetoa() assumes kDosageMax rather than kDosageMid multiplier
        const uint64_t unknown_hethet_int_ddosage = S_CAST(int64_t, unknown_hethet_d * kDosageMax);
        write_iter = ddosagetoa(unknown_hethet_int_ddosage, write_iter);
        write_iter = strcpya_k(write_iter, " het pair");
        if (unknown_hethet_int_ddosage != kDosageMax) {
          *write_iter++ = 's';
        }
        write_iter = strcpya_k(write_iter, " statistically phased");
      }
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
          // bugfix (18 Oct 2023): need to pass bit count, not word count
          AlignedBitarrInvert(founder_ct, nosex_collapsed);
        }
      }
      // Unlike plink 1.9, we don't restrict these HWE computations to the
      // nonmissing intersection.
      for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
        const uintptr_t* cur_genovec = pgvs[var_idx].genovec;
        STD_ARRAY_DECL(uint32_t, 4, genocounts);
        GenoarrCountFreqsUnsafe(cur_genovec, founder_ct, genocounts);
        double hwe_pval;
        if (!is_xs[var_idx]) {
          hwe_pval = HweP(genocounts[1], genocounts[0], genocounts[2], hwe_midp);
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
          hwe_pval = HweXchrP(genocounts[1], genocounts[0] - male_genocounts[0], genocounts[2] - male_genocounts[2], male_genocounts[0], male_genocounts[2], hwe_midp);
        }
        logprintf("  %s: %g\n", ld_console_varids[var_idx], hwe_pval);
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
      write_iter = strcpya_k(write_iter, "    D' = ");
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
      write_iter = dtoa_f_probp6_spaced(freq_xmaj * freq_majx, write_iter);
      write_iter = strcpya_k(write_iter, ")  (");
      write_iter = dtoa_f_probp6_clipped(freq_xmin * freq_majx, write_iter);
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
      write_iter = dtoa_f_probp6_spaced(freq_xmaj * freq_minx, write_iter);
      write_iter = strcpya_k(write_iter, ")  (");
      write_iter = dtoa_f_probp6_clipped(freq_xmin * freq_minx, write_iter);
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
//    a. If not loading phase: {one_bitvec, two_bitvec, nm_bitvec},
//                             founder_ctaw words each, i.e. 3 * bitvec_byte_ct
//                             uint32_t nmaj_ct, nm_ct for another
//                             RoundUpPow2(8, kBytesPerVec)
//    b. If loading phase, also need phasepresent and phaseinfo, founder_ctaw
//       words each, for a total of 5 * bitvec_byte_ct +
//       RoundUpPow2(8, kBytesPerVec)
// 2. If loading dosage:
//    a. If not loading phase: {dosage_vec, dosage_uhet} founder_dosagev_ct
//                             vectors each; plus nm_bitvec, i.e.
//                             2 * dosagevec_byte_ct + bitvec_byte_ct
//                             uint64_t nmaj_dosage, uint32_t nm_ct for another
//                             RoundUpPow2(12, kBytesPerVec)
//    b. If loading phase, also need main_dphase_deltas, for a total of
//       3 * dosagevec_byte_ct + bitvec_byte_ct
// 3. If chrX, entire nonmale part first, then unphased male part
// 4. If chrY, just unphased male part
//
// Probable todo: Unpack the index variant in a way that allows r^2 to be
// computed efficiently against other variants *without* unpacking them.

// In the HighmemUnpack step, the main thread wants to prepare the next copy of
// this while the worker threads are processing the current copy.
typedef struct ClumpCtxAlternating {
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
  const uintptr_t* founder_male;
  const uint32_t* founder_male_cumulative_popcounts;
  const uintptr_t* founder_nonmale;
  const uint32_t* founder_nonmale_cumulative_popcounts;
  uint32_t raw_sample_ct;
  uint32_t founder_ct;
  uint32_t founder_male_ct;
  // Precomputed for convenience.
  uintptr_t pgv_byte_stride;
  uintptr_t bitvec_byte_ct;
  uintptr_t dosagevec_byte_ct;
  uintptr_t male_bitvec_byte_ct;
  uintptr_t male_dosagevec_byte_ct;
  uintptr_t nonmale_bitvec_byte_ct;
  uintptr_t nonmale_dosagevec_byte_ct;
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
  // Unpacked variant representation can vary based on is_x, is_y, phase_type,
  // and load_dosage.
  unsigned char* unpacked_variants;

  // chrX unpacking workspaces
  // NypCtToCachelineCt(max(founder_nonmale_ct, founder_male_ct)) cachelines
  // each
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
uintptr_t UnpackedByteStride(const ClumpCtx* ctx, R2PhaseType phase_type, uint32_t load_dosage) {
  uintptr_t trail_byte_ct;
  if (load_dosage) {
    trail_byte_ct = (1 + (phase_type == kR2PhaseTypeUnphased)) * sizeof(int64_t) + sizeof(int32_t);
  } else {
    trail_byte_ct = (2 + (phase_type == kR2PhaseTypeUnphased)) * sizeof(int32_t);
  }
  trail_byte_ct = RoundUpPow2(trail_byte_ct, kBytesPerVec);
  if (ctx->is_x) {
    const uintptr_t nonmale_bitvec_byte_ct = ctx->nonmale_bitvec_byte_ct;
    const uintptr_t male_bitvec_byte_ct = ctx->male_bitvec_byte_ct;
    if (load_dosage) {
      return (1 + phase_type) * ctx->nonmale_dosagevec_byte_ct + (1 + (phase_type != kR2PhaseTypeUnphased)) * ctx->male_dosagevec_byte_ct + nonmale_bitvec_byte_ct + male_bitvec_byte_ct + 2 * trail_byte_ct;
    }
    return (3 + 2 * (phase_type == kR2PhaseTypePresent)) * nonmale_bitvec_byte_ct + 3 * male_bitvec_byte_ct + 2 * trail_byte_ct;
  } else if (ctx->is_y) {
    assert(phase_type != kR2PhaseTypePresent);
    const uintptr_t male_bitvec_byte_ct = ctx->male_bitvec_byte_ct;
    if (load_dosage) {
      return (1 + (phase_type != kR2PhaseTypeUnphased)) * ctx->dosagevec_byte_ct + male_bitvec_byte_ct + trail_byte_ct;
    }
    return 3 * male_bitvec_byte_ct + trail_byte_ct;
  }
  const uintptr_t bitvec_byte_ct = ctx->bitvec_byte_ct;
  if (load_dosage) {
    return (1 + phase_type) * ctx->dosagevec_byte_ct + bitvec_byte_ct + trail_byte_ct;
  }
  return (3 + 2 * (phase_type == kR2PhaseTypePresent)) * bitvec_byte_ct + trail_byte_ct;
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

void LdUnpackNondosage(const PgenVariant* pgvp, uint32_t sample_ct, R2PhaseType phase_type, unsigned char* dst_iter) {
  const uintptr_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  uintptr_t* dst_witer = R_CAST(uintptr_t*, dst_iter);
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
  uint32_t* dst_u32iter = R_CAST(uint32_t*, dst_witer);
  uint32_t* nmaj_ct_ptr = dst_u32iter;
  ++dst_u32iter;
  uint32_t* nm_ct_ptr = dst_u32iter;
  GenoarrSplit12Nm(pgvp->genovec, sample_ct, one_bitvec, two_bitvec, nm_bitvec);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  *nmaj_ct_ptr = GenoBitvecSum(one_bitvec, two_bitvec, sample_ctl);
  *nm_ct_ptr = PopcountWords(nm_bitvec, sample_ctl);
  if (phase_type == kR2PhaseTypePresent) {
    if (!pgvp->phasepresent_ct) {
      ZeroWArr(sample_ctl, phasepresent);
    } else {
      memcpy(phasepresent, pgvp->phasepresent, sample_ctl * sizeof(intptr_t));
      memcpy(phaseinfo, pgvp->phaseinfo, sample_ctl * sizeof(intptr_t));
    }
  } else if (phase_type == kR2PhaseTypeUnphased) {
    uint32_t* ssq_ptr = &(dst_u32iter[1]);
    *ssq_ptr = (*nmaj_ct_ptr) + 2 * PopcountWords(two_bitvec, sample_ctl);
  }
}

// compiler should optimize out the conditional when both values are the same?
static inline uint32_t LdNondosageTrailAlignedByteCt(R2PhaseType phase_type) {
  return (phase_type == kR2PhaseTypeUnphased)? RoundUpPow2(3 * sizeof(int32_t), kBytesPerVec) : RoundUpPow2(2 * sizeof(int32_t), kBytesPerVec);
}

static inline uint32_t LdDosageTrailAlignedByteCt(R2PhaseType phase_type) {
  return (phase_type == kR2PhaseTypeUnphased)? RoundUpPow2(sizeof(int64_t) + sizeof(double) + sizeof(int32_t), kBytesPerVec) : RoundUpPow2(sizeof(int64_t) + sizeof(int32_t), kBytesPerVec);
}

unsigned char* LdUnpackNondosageSubset(const PgenVariant* pgvp, const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, R2PhaseType phase_type, unsigned char* dst_iter, uintptr_t* workspace) {
  const uintptr_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  uintptr_t* dst_witer = R_CAST(uintptr_t*, dst_iter);
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
  uint32_t* final_u32s = R_CAST(uint32_t*, dst_witer);
  uint32_t* nmaj_ct_ptr = &(final_u32s[0]);
  uint32_t* nm_ct_ptr = &(final_u32s[1]);
  dst_iter = R_CAST(unsigned char*, dst_witer);
  dst_iter = &(dst_iter[LdNondosageTrailAlignedByteCt(phase_type)]);

  uintptr_t* genovec_collapsed = workspace;
  CopyNyparrNonemptySubset(pgvp->genovec, sample_include, raw_sample_ct, sample_ct, genovec_collapsed);
  GenoarrSplit12Nm(genovec_collapsed, sample_ct, one_bitvec, two_bitvec, nm_bitvec);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  *nmaj_ct_ptr = GenoBitvecSum(one_bitvec, two_bitvec, sample_ctl);
  *nm_ct_ptr = PopcountWords(nm_bitvec, sample_ctl);
  if (phase_type == kR2PhaseTypePresent) {
    if (!pgvp->phasepresent_ct) {
      ZeroWArr(sample_ctl, phasepresent);
    } else {
      CopyBitarrSubset(pgvp->phasepresent, sample_include, sample_ct, phasepresent);
      CopyBitarrSubset(pgvp->phaseinfo, sample_include, sample_ct, phaseinfo);
    }
  } else if (phase_type == kR2PhaseTypeUnphased) {
    uint32_t* ssq_ptr = &(final_u32s[2]);
    *ssq_ptr = (*nm_ct_ptr) + 2 * PopcountWords(two_bitvec, sample_ctl);
  }
  return dst_iter;
}

// Treat this as two mini-records, one with males and one with nonmales.
// (This doesn't work for inter-chr case.)
void LdUnpackChrXNondosage(const PgenVariant* pgvp, const uintptr_t* founder_male, const uintptr_t* founder_nonmale, uint32_t raw_sample_ct, uint32_t founder_male_ct, uint32_t founder_nonmale_ct, R2PhaseType phase_type, unsigned char* dst_iter, uintptr_t* workspace) {
  // Unpack male data, ignoring phase.
  dst_iter = LdUnpackNondosageSubset(pgvp, founder_male, raw_sample_ct, founder_male_ct, R2PhaseOmit(phase_type), dst_iter, workspace);
  // Then unpack nonmale data.
  LdUnpackNondosageSubset(pgvp, founder_nonmale, raw_sample_ct, founder_nonmale_ct, phase_type, dst_iter, workspace);
}

void LdUnpackDosage(const PgenVariant* pgvp, uint32_t sample_ct, R2PhaseType phase_type, unsigned char* dst_iter) {
  const uintptr_t dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  const uintptr_t dosagevec_byte_ct = dosagev_ct * kBytesPerVec;
  Dosage* dosage_vec = R_CAST(Dosage*, dst_iter);
  dst_iter = &(dst_iter[dosagevec_byte_ct]);
  Dosage* dosage_uhet = nullptr;
  if (phase_type != kR2PhaseTypeUnphased) {
    dosage_uhet = R_CAST(Dosage*, dst_iter);
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
  // In 32-bit build, no alignment guarantee for nmaj_dosage or
  // nmaj_dosage_ssq.
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
  GenoarrToNonmissingnessUnsafe(pgvp->genovec, sample_ct, nm_bitvec);
  ZeroTrailingBits(sample_ct, nm_bitvec);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  BitvecOr(pgvp->dosage_present, sample_ctl, nm_bitvec);
  *nm_ct_ptr = PopcountWords(nm_bitvec, sample_ctl);
  if (phase_type == kR2PhaseTypeUnphased) {
    // todo: check if dedicated function is faster
    const uint64_t nmaj_dosage_ssq = DosageUnsignedDotprod(dosage_vec, dosage_vec, dosagev_ct);
    memcpy(nmaj_dosage_ssq_uc_ptr, &nmaj_dosage_ssq, sizeof(int64_t));
  } else {
    FillDosageUhet(dosage_vec, dosagev_ct, dosage_uhet);
    if (phase_type == kR2PhaseTypePresent) {
      PopulateDenseDphase(pgvp->phasepresent, pgvp->phaseinfo, pgvp->dosage_present, dosage_vec, pgvp->dphase_present, pgvp->dphase_delta, sample_ct, pgvp->phasepresent_ct, pgvp->dosage_ct, pgvp->dphase_ct, dosage_uhet, dense_dphase_delta);
    }
  }
}

unsigned char* LdUnpackDosageSubset(const PgenVariant* pgvp, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t raw_sample_ct, uint32_t sample_ct, R2PhaseType phase_type, unsigned char* dst_iter, uintptr_t* workspace) {
  const uintptr_t dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  const uintptr_t dosagevec_byte_ct = dosagev_ct * kBytesPerVec;
  Dosage* dosage_vec = R_CAST(Dosage*, dst_iter);
  dst_iter = &(dst_iter[dosagevec_byte_ct]);
  Dosage* dosage_uhet = nullptr;
  if (phase_type != kR2PhaseTypeUnphased) {
    dosage_uhet = R_CAST(Dosage*, dst_iter);
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
  // In 32-bit build, no alignment guarantee for nmaj_dosage.
  unsigned char* nmaj_dosage_uc_ptr = dst_iter;
  uint32_t dst_offset = sizeof(int64_t);
  unsigned char* nmaj_dosage_ssq_uc_ptr = nullptr;
  if (phase_type == kR2PhaseTypeUnphased) {
    nmaj_dosage_ssq_uc_ptr = &(dst_iter[dst_offset]);
    dst_offset += sizeof(int64_t);
  }
  uint32_t* nm_ct_ptr = R_CAST(uint32_t*, &(dst_iter[dst_offset]));
  dst_offset += sizeof(int32_t);
  dst_iter = &(dst_iter[RoundUpPow2(dst_offset, kBytesPerVec)]);

  const uintptr_t* dosage_present = pgvp->dosage_present;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  {
    uintptr_t* genovec_collapsed = workspace;
    PopulateDenseDosageNonemptySubset(sample_include, sample_include_cumulative_popcounts, pgvp->genovec, dosage_present, pgvp->dosage_main, raw_sample_ct, sample_ct, pgvp->dosage_ct, dosage_vec, genovec_collapsed);
    const uint64_t nmaj_dosage = DenseDosageSum(dosage_vec, dosagev_ct);
    memcpy(nmaj_dosage_uc_ptr, &nmaj_dosage, sizeof(int64_t));
    GenoarrToNonmissingnessUnsafe(genovec_collapsed, sample_ct, nm_bitvec);
    ZeroTrailingBits(sample_ct, nm_bitvec);
  }
  uintptr_t* dosage_present_collapsed = workspace;
  CopyBitarrSubset(dosage_present, sample_include, sample_ct, dosage_present_collapsed);
  BitvecOr(dosage_present_collapsed, sample_ctl, nm_bitvec);
  *nm_ct_ptr = PopcountWords(nm_bitvec, sample_ctl);
  if (phase_type == kR2PhaseTypeUnphased) {
    const uint64_t nmaj_dosage_ssq = DosageUnsignedDotprod(dosage_vec, dosage_vec, dosagev_ct);
    memcpy(nmaj_dosage_ssq_uc_ptr, &nmaj_dosage_ssq, sizeof(uint64_t));
  } else {
    FillDosageUhet(dosage_vec, dosagev_ct, dosage_uhet);
    if (phase_type == kR2PhaseTypePresent) {
      PopulateDenseDphaseSubset(sample_include, sample_include_cumulative_popcounts, pgvp->phasepresent, pgvp->phaseinfo, pgvp->dosage_present, dosage_vec, pgvp->dphase_present, pgvp->dphase_delta, raw_sample_ct, sample_ct, pgvp->phasepresent_ct, pgvp->dosage_ct, pgvp->dphase_ct, dosage_uhet, dense_dphase_delta);
    }
  }
  return dst_iter;
}

void LdUnpackChrXDosage(const PgenVariant* pgvp, const uintptr_t* founder_male, const uint32_t* founder_male_cumulative_popcounts, const uintptr_t* founder_nonmale, const uint32_t* founder_nonmale_cumulative_popcounts, uint32_t raw_sample_ct, uint32_t founder_male_ct, uint32_t founder_nonmale_ct, R2PhaseType phase_type, unsigned char* dst_iter, uintptr_t* workspace) {
  // Unpack male data, ignoring phase.
  dst_iter = LdUnpackDosageSubset(pgvp, founder_male, founder_male_cumulative_popcounts, raw_sample_ct, founder_male_ct, R2PhaseOmit(phase_type), dst_iter, workspace);
  // Then unpack nonmale data.
  LdUnpackDosageSubset(pgvp, founder_nonmale, founder_nonmale_cumulative_popcounts, raw_sample_ct, founder_nonmale_ct, phase_type, dst_iter, workspace);
}

typedef struct R2NondosageVariantStruct {
  const uintptr_t* one_bitvec;
  const uintptr_t* two_bitvec;
  const uintptr_t* nm_bitvec;
  const uintptr_t* phasepresent; // may be uninitialized
  const uintptr_t* phaseinfo; // may be uninitialized
  uint32_t nmaj_ct;
  uint32_t nm_ct;
  uint32_t ssq; // may be uninitialized
} R2NondosageVariant;

typedef struct R2DosageVariantStruct {
  const Dosage* dosage_vec;
  const Dosage* dosage_uhet; // may be uninitialized
  const uintptr_t* nm_bitvec;
  const SDosage* dense_dphase_delta; // may be uninitialized
  uint64_t nmaj_dosage;
  uint64_t nmaj_dosage_ssq; // may be uninitialized
  uint32_t nm_ct;
} R2DosageVariant;

typedef union {
  R2NondosageVariant nd;
  R2DosageVariant d;
  R2NondosageVariant x_nd[2];
  R2DosageVariant x_d[2];
} R2Variant;

const unsigned char* FillR2Nondosage(const unsigned char* src_iter, uint32_t sample_ct, R2PhaseType phase_type, R2NondosageVariant* ndp) {
  // See LdUnpackNondosage().
  const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  const uintptr_t* src_witer = R_CAST(const uintptr_t*, src_iter);
  ndp->one_bitvec = src_witer;
  src_witer = &(src_witer[sample_ctaw]);
  ndp->two_bitvec = src_witer;
  src_witer = &(src_witer[sample_ctaw]);
  ndp->nm_bitvec = src_witer;
  src_witer = &(src_witer[sample_ctaw]);
  if (phase_type == kR2PhaseTypePresent) {
    ndp->phasepresent = src_witer;
    src_witer = &(src_witer[sample_ctaw]);
    ndp->phaseinfo = src_witer;
    src_witer = &(src_witer[sample_ctaw]);
  }
  const uint32_t* final_u32s = R_CAST(const uint32_t*, src_witer);
  ndp->nmaj_ct = final_u32s[0];
  ndp->nm_ct = final_u32s[1];
  if (phase_type == kR2PhaseTypeUnphased) {
    ndp->ssq = final_u32s[2];
  }
  src_iter = R_CAST(const unsigned char*, src_witer);
  return &(src_iter[LdNondosageTrailAlignedByteCt(phase_type)]);
}

const unsigned char* FillR2Dosage(const unsigned char* src_iter, uint32_t sample_ct, R2PhaseType phase_type, R2DosageVariant* dp) {
  // See LdUnpackDosage().
  const uintptr_t dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  const uintptr_t dosagevec_byte_ct = dosagev_ct * kBytesPerVec;
  dp->dosage_vec = R_CAST(const Dosage*, src_iter);
  src_iter = &(src_iter[dosagevec_byte_ct]);
  if (phase_type != kR2PhaseTypeUnphased) {
    dp->dosage_uhet = R_CAST(const Dosage*, src_iter);
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
  return &(src_iter[LdDosageTrailAlignedByteCt(phase_type)]);
}

void FillR2V(const unsigned char* src_iter, uint32_t sample_ct, R2PhaseType phase_type, uint32_t load_dosage, R2Variant* r2vp) {
  if (!load_dosage) {
    FillR2Nondosage(src_iter, sample_ct, phase_type, &(r2vp->nd));
  } else {
    FillR2Dosage(src_iter, sample_ct, phase_type, &(r2vp->d));
  }
}

void FillXR2V(const unsigned char* unpacked_variant, uint32_t male_ct, uint32_t nonmale_ct, R2PhaseType phase_type, uint32_t load_dosage, R2Variant* r2vp) {
  if (!load_dosage) {
    const unsigned char* src_iter = FillR2Nondosage(unpacked_variant, male_ct, R2PhaseOmit(phase_type), &(r2vp->x_nd[0]));
    FillR2Nondosage(src_iter, nonmale_ct, phase_type, &(r2vp->x_nd[1]));
  } else {
    const unsigned char* src_iter = FillR2Dosage(unpacked_variant, male_ct, R2PhaseOmit(phase_type), &(r2vp->x_d[0]));
    FillR2Dosage(src_iter, nonmale_ct, phase_type, &(r2vp->x_d[1]));
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
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t founder_nonmale_ct = ctx->founder_ct - founder_male_ct;
  const uintptr_t* founder_male = nullptr;
  const uintptr_t* founder_nonmale = nullptr;
  const uint32_t* founder_male_cumulative_popcounts = nullptr;
  const uint32_t* founder_nonmale_cumulative_popcounts = nullptr;
  uintptr_t* chrx_workspace = nullptr;
  const uintptr_t* cur_sample_include = nullptr;
  PgrSampleSubsetIndex pssi;
  uint32_t cur_sample_ct;
  {
    const uint32_t is_x = ctx->is_x;
    const uint32_t is_y = ctx->is_y;
    if (is_x) {
      founder_male = ctx->founder_male;
      founder_male_cumulative_popcounts = ctx->founder_male_cumulative_popcounts;
      founder_nonmale = ctx->founder_nonmale;
      founder_nonmale_cumulative_popcounts = ctx->founder_nonmale_cumulative_popcounts;
      cur_sample_ct = ctx->raw_sample_ct;
      PgrClearSampleSubsetIndex(pgrp, &pssi);
      chrx_workspace = ctx->chrx_workspaces[tidx];
    } else {
      const uint32_t* cur_sample_include_cumulative_popcounts;
      if (is_y) {
        cur_sample_include = ctx->founder_male;
        cur_sample_include_cumulative_popcounts = ctx->founder_male_cumulative_popcounts;
        cur_sample_ct = founder_male_ct;
      } else {
        cur_sample_include = ctx->founder_info;
        cur_sample_include_cumulative_popcounts = ctx->founder_info_cumulative_popcounts;
        cur_sample_ct = ctx->founder_ct;
      }
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
    }
  }
  const R2PhaseType phase_type = S_CAST(R2PhaseType, ctx->phase_type);
  const uint32_t load_dosage = ctx->load_dosage;
  const uintptr_t allele_idx_start = IdxToUidxW(observed_alleles, observed_alleles_cumulative_popcounts_w, ctx->allele_widx_start, ctx->allele_widx_end, oaidx);
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
    if (!chrx_workspace) {
      if (load_dosage) {
        if (phase_type == kR2PhaseTypePresent) {
          reterr = PgrGetInv1Dp(cur_sample_include, pssi, cur_sample_ct, variant_uidx, aidx, pgrp, &pgv);
        } else {
          reterr = PgrGetInv1D(cur_sample_include, pssi, cur_sample_ct, variant_uidx, aidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
        }
        if (unlikely(reterr)) {
          goto ClumpHighmemUnpack_err;
        }
        LdUnpackDosage(&pgv, cur_sample_ct, phase_type, write_iter);
      } else {
        if (phase_type == kR2PhaseTypePresent) {
          reterr = PgrGetInv1P(cur_sample_include, pssi, cur_sample_ct, variant_uidx, aidx, pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
        } else {
          reterr = PgrGetInv1(cur_sample_include, pssi, cur_sample_ct, variant_uidx, aidx, pgrp, pgv.genovec);
        }
        if (unlikely(reterr)) {
          goto ClumpHighmemUnpack_err;
        }
        LdUnpackNondosage(&pgv, cur_sample_ct, phase_type, write_iter);
      }
    } else {
      // chrX
      if (load_dosage) {
        if (phase_type == kR2PhaseTypePresent) {
          reterr = PgrGetInv1Dp(nullptr, pssi, cur_sample_ct, variant_uidx, aidx, pgrp, &pgv);
        } else {
          reterr = PgrGetInv1D(nullptr, pssi, cur_sample_ct, variant_uidx, aidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
        }
        if (unlikely(reterr)) {
          goto ClumpHighmemUnpack_err;
        }
        LdUnpackChrXDosage(&pgv, founder_male, founder_male_cumulative_popcounts, founder_nonmale, founder_nonmale_cumulative_popcounts, cur_sample_ct, founder_male_ct, founder_nonmale_ct, phase_type, write_iter, chrx_workspace);
      } else {
        if (phase_type == kR2PhaseTypePresent) {
          reterr = PgrGetInv1P(nullptr, pssi, cur_sample_ct, variant_uidx, aidx, pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
        } else {
          reterr = PgrGetInv1(nullptr, pssi, cur_sample_ct, variant_uidx, aidx, pgrp, pgv.genovec);
        }
        if (unlikely(reterr)) {
          goto ClumpHighmemUnpack_err;
        }
        LdUnpackChrXNondosage(&pgv, founder_male, founder_nonmale, cur_sample_ct, founder_male_ct, founder_nonmale_ct, phase_type, write_iter, chrx_workspace);
      }
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

// Main return value is valid_obs_ct.  On valid_obs_ct=0, other return values
// are not filled.  (Same is true for the next three functions.)
uint32_t ComputeR2NondosageUnphasedStats(const R2NondosageVariant* ndp0, const R2NondosageVariant* ndp1, uint32_t sample_ct, uint32_t* nmaj_ct0_ptr, uint32_t* nmaj_ct1_ptr, uint32_t* ssq0_ptr, uint32_t* ssq1_ptr, uint32_t* dotprod_ptr) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uintptr_t* nm_bitvec0 = ndp0->nm_bitvec;
  const uintptr_t* nm_bitvec1 = ndp1->nm_bitvec;
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
  const uintptr_t* one_bitvec0 = ndp0->one_bitvec;
  const uintptr_t* two_bitvec0 = ndp0->two_bitvec;
  if (nm_ct0 == valid_obs_ct) {
    *nmaj_ct0_ptr = ndp0->nmaj_ct;
    *ssq0_ptr = ndp0->ssq;
  } else {
    const uint32_t nmaj_ct0 = GenoBitvecSumSubset(nm_bitvec1, one_bitvec0, two_bitvec0, sample_ctl);
    *nmaj_ct0_ptr = nmaj_ct0;
    // 0, 1, 4 instead of 0, 1, 2
    *ssq0_ptr = nmaj_ct0 + 2 * PopcountWordsIntersect(nm_bitvec1, two_bitvec0, sample_ctl);
  }
  const uintptr_t* one_bitvec1 = ndp1->one_bitvec;
  const uintptr_t* two_bitvec1 = ndp1->two_bitvec;
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
  const uintptr_t* nm_bitvec0 = ndp0->nm_bitvec;
  const uintptr_t* nm_bitvec1 = ndp1->nm_bitvec;
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
  const uintptr_t* one_bitvec0 = ndp0->one_bitvec;
  const uintptr_t* two_bitvec0 = ndp0->two_bitvec;
  uint32_t nmaj_ct0 = ndp0->nmaj_ct;
  if (nm_ct0 != valid_obs_ct) {
    nmaj_ct0 = GenoBitvecSumSubset(nm_bitvec1, one_bitvec0, two_bitvec0, sample_ctl);
  }
  const uintptr_t* one_bitvec1 = ndp1->one_bitvec;
  const uintptr_t* two_bitvec1 = ndp1->two_bitvec;
  uint32_t nmaj_ct1 = ndp1->nmaj_ct;
  if (nm_ct1 != valid_obs_ct) {
    nmaj_ct1 = GenoBitvecSumSubset(nm_bitvec0, one_bitvec1, two_bitvec1, sample_ctl);
  }
  uint32_t known_dotprod;
  uint32_t unknown_hethet_ct;
  GenoBitvecPhasedDotprod(one_bitvec0, two_bitvec0, one_bitvec1, two_bitvec1, sample_ctl, &known_dotprod, &unknown_hethet_ct);
  if ((phase_type == kR2PhaseTypePresent) && (unknown_hethet_ct != 0)) {
    // don't bother with no-phase-here optimization for now
    HardcallPhasedR2Refine(ndp0->phasepresent, ndp0->phaseinfo, ndp1->phasepresent, ndp1->phaseinfo, sample_ctl, &known_dotprod, &unknown_hethet_ct);
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
  const uint32_t valid_obs_ct = DosageR2Prod(dosage_vec0, nm_bitvec0, dosage_vec1, nm_bitvec1, sample_ct, nm_ct0, nm_ct1, nmaj_dosages, dosageprod_ptr);
  if (!valid_obs_ct) {
    return 0;
  }
  const uint32_t sample_dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  if (nm_ct0 == valid_obs_ct) {
    *ssq0_ptr = dp0->nmaj_dosage_ssq;
  } else {
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
  uint64_t dosageprod;
  const uint32_t valid_obs_ct = DosageR2Prod(dosage_vec0, nm_bitvec0, dosage_vec1, nm_bitvec1, sample_ct, nm_ct0, nm_ct1, nmaj_dosages, &dosageprod);
  if (!valid_obs_ct) {
    return 0;
  }
  const uint32_t sample_dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  const Dosage* dosage_uhet0 = dp0->dosage_uhet;
  const Dosage* dosage_uhet1 = dp1->dosage_uhet;
  uint64_t uhethet_dosageprod = DosageUnsignedNomissDotprod(dosage_uhet0, dosage_uhet1, sample_dosagev_ct);
  if ((phase_type == kR2PhaseTypePresent) && (uhethet_dosageprod != 0)) {
    const SDosage* dphase_delta0 = dp0->dense_dphase_delta;
    const SDosage* dphase_delta1 = dp1->dense_dphase_delta;
    dosageprod = S_CAST(int64_t, dosageprod) + DosageSignedDotprod(dphase_delta0, dphase_delta1, sample_dosagev_ct);
  }
  nmajsums_d[0] = u63tod(nmaj_dosages[0]) * kRecipDosageMid;
  nmajsums_d[1] = u63tod(nmaj_dosages[1]) * kRecipDosageMid;
  *known_dotprod_d_ptr = u63tod(dosageprod - uhethet_dosageprod) * (kRecipDosageMidSq * 0.5);
  *unknown_hethet_d_ptr = u63tod(uhethet_dosageprod) * kRecipDosageMidSq;
  return valid_obs_ct;
}

double ComputeR2(const R2Variant* r2vp0, const R2Variant* r2vp1, uint32_t sample_ct, R2PhaseType phase_type, uint32_t load_dosage) {
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
      return cov01 * cov01 / variance_prod;
    }
    valid_obs_ct = ComputeR2DosagePhasedStats(dp0, dp1, sample_ct, phase_type, nmajsums_d, &known_dotprod_d, &unknown_hethet_d);
  }
  if (!valid_obs_ct) {
    return -DBL_MAX;
  }
  const double twice_tot_recip = 0.5 / u31tod(valid_obs_ct);
  double r2;
  const LDErr ld_err = PhasedLD(nmajsums_d, known_dotprod_d, unknown_hethet_d, twice_tot_recip, 0, nullptr, &r2);
  return (ld_err == kLDErrNone)? r2 : -DBL_MAX;
}

// todo: single-part X, for --r2 inter-chr case
double ComputeTwoPartXR2(const R2Variant* r2vp0, const R2Variant* r2vp1, uint32_t male_ct, uint32_t nonmale_ct, R2PhaseType phase_type, uint32_t load_dosage) {
  // initially fill these with male values
  double nmajsums_d[2];
  double known_dotprod_d;

  double nonmale_nmajsums_d[2];
  double nonmale_known_dotprod_d;
  double nonmale_unknown_hethet_d;
  uint32_t male_obs_ct;
  uint32_t nonmale_obs_ct;
  if (!load_dosage) {
    const R2NondosageVariant* male_vp0 = &(r2vp0->x_nd[0]);
    const R2NondosageVariant* male_vp1 = &(r2vp1->x_nd[0]);
    const R2NondosageVariant* nonmale_vp0 = &(r2vp0->x_nd[1]);
    const R2NondosageVariant* nonmale_vp1 = &(r2vp1->x_nd[1]);
    if (phase_type == kR2PhaseTypeUnphased) {
      uint32_t nmaj_ct0;
      uint32_t nmaj_ct1;
      uint32_t male_ssq0;
      uint32_t male_ssq1;
      uint32_t male_dotprod;
      male_obs_ct = ComputeR2NondosageUnphasedStats(male_vp0, male_vp1, male_ct, &nmaj_ct0, &nmaj_ct1, &male_ssq0, &male_ssq1, &male_dotprod);
      if (!male_obs_ct) {
        nmaj_ct0 = 0;
        nmaj_ct1 = 0;
        male_ssq0 = 0;
        male_ssq1 = 0;
        male_dotprod = 0;
      }
      uint32_t nonmale_nmaj_ct0;
      uint32_t nonmale_nmaj_ct1;
      uint32_t nonmale_ssq0;
      uint32_t nonmale_ssq1;
      uint32_t nonmale_dotprod;
      nonmale_obs_ct = ComputeR2NondosageUnphasedStats(nonmale_vp0, nonmale_vp1, nonmale_ct, &nonmale_nmaj_ct0, &nonmale_nmaj_ct1, &nonmale_ssq0, &nonmale_ssq1, &nonmale_dotprod);
      const uint32_t weighted_obs_ct = male_obs_ct + 2 * nonmale_obs_ct;
      if (!weighted_obs_ct) {
        return -DBL_MAX;
      }
      // these can overflow uint32
      uint64_t ssq0 = male_ssq0;
      uint64_t ssq1 = male_ssq1;
      uint64_t dotprod = male_dotprod;

      if (nonmale_obs_ct) {
        nmaj_ct0 += 2 * nonmale_nmaj_ct0;
        nmaj_ct1 += 2 * nonmale_nmaj_ct1;
        ssq0 += 2LLU * nonmale_ssq0;
        ssq1 += 2LLU * nonmale_ssq1;
        dotprod += 2LLU * nonmale_dotprod;
      }
      // individual terms can exceed 2^63, but difference cannot
      const uint64_t variance0_u64 = ssq0 * S_CAST(uint64_t, weighted_obs_ct) - S_CAST(int64_t, nmaj_ct0) * nmaj_ct0;
      const uint64_t variance1_u64 = ssq1 * S_CAST(uint64_t, weighted_obs_ct) - S_CAST(int64_t, nmaj_ct1) * nmaj_ct1;
      const double variance_prod = u63tod(variance0_u64) * u63tod(variance1_u64);
      if (variance_prod == 0.0) {
        return -DBL_MAX;
      }
      const double cov01 = u63tod(dotprod * S_CAST(uint64_t, weighted_obs_ct) - S_CAST(uint64_t, nmaj_ct0) * nmaj_ct1);
      return cov01 * cov01 / variance_prod;
    }
    double ignore;
    male_obs_ct = ComputeR2NondosagePhasedStats(male_vp0, male_vp1, male_ct, R2PhaseOmit(phase_type), nmajsums_d, &known_dotprod_d, &ignore);
    nonmale_obs_ct = ComputeR2NondosagePhasedStats(nonmale_vp0, nonmale_vp1, nonmale_ct, phase_type, nonmale_nmajsums_d, &nonmale_known_dotprod_d, &nonmale_unknown_hethet_d);
  } else {
    const R2DosageVariant* male_vp0 = &(r2vp0->x_d[0]);
    const R2DosageVariant* male_vp1 = &(r2vp1->x_d[0]);
    const R2DosageVariant* nonmale_vp0 = &(r2vp0->x_d[1]);
    const R2DosageVariant* nonmale_vp1 = &(r2vp1->x_d[1]);
    if (phase_type == kR2PhaseTypeUnphased) {
      uint64_t nmaj_dosages[2];
      uint64_t dosageprod;
      uint64_t ssq0;
      uint64_t ssq1;
      male_obs_ct = ComputeR2DosageUnphasedStats(male_vp0, male_vp1, male_ct, nmaj_dosages, &dosageprod, &ssq0, &ssq1);
      if (!male_obs_ct) {
        nmaj_dosages[0] = 0;
        nmaj_dosages[1] = 0;
        dosageprod = 0;
        ssq0 = 0;
        ssq1 = 0;
      }
      uint64_t nonmale_nmaj_dosages[2];
      uint64_t nonmale_dosageprod;
      uint64_t nonmale_ssq0;
      uint64_t nonmale_ssq1;
      nonmale_obs_ct = ComputeR2DosageUnphasedStats(nonmale_vp0, nonmale_vp1, nonmale_ct, nonmale_nmaj_dosages, &nonmale_dosageprod, &nonmale_ssq0, &nonmale_ssq1);
      const uint32_t weighted_obs_ct = male_ct + 2 * nonmale_obs_ct;
      if (!weighted_obs_ct) {
        return -DBL_MAX;
      }
      if (nonmale_obs_ct) {
        nmaj_dosages[0] += 2 * nonmale_nmaj_dosages[0];
        nmaj_dosages[1] += 2 * nonmale_nmaj_dosages[1];
        // bugfix (27 Oct 2023)
        dosageprod += 2 * nonmale_dosageprod;
        ssq0 += 2 * nonmale_ssq0;
        ssq1 += 2 * nonmale_ssq1;
      }
      const double variance0 = u127prod_diff_d(ssq0, weighted_obs_ct, nmaj_dosages[0], nmaj_dosages[0]);
      const double variance1 = u127prod_diff_d(ssq1, weighted_obs_ct, nmaj_dosages[1], nmaj_dosages[1]);
      const double variance_prod = variance0 * variance1;
      if (variance_prod == 0.0) {
        return -DBL_MAX;
      }
      const double cov01 = i127prod_diff_d(dosageprod, weighted_obs_ct, nmaj_dosages[0], nmaj_dosages[1]);
      return cov01 * cov01 / variance_prod;
    }
    double ignore;
    male_obs_ct = ComputeR2DosagePhasedStats(male_vp0, male_vp1, male_ct, R2PhaseOmit(phase_type), nmajsums_d, &known_dotprod_d, &ignore);
    nonmale_obs_ct = ComputeR2DosagePhasedStats(nonmale_vp0, nonmale_vp1, nonmale_ct, phase_type, nonmale_nmajsums_d, &nonmale_known_dotprod_d, &nonmale_unknown_hethet_d);
  }
  const uint32_t weighted_obs_ct = male_obs_ct + 2 * nonmale_obs_ct;
  if (!weighted_obs_ct) {
    return -DBL_MAX;
  }
  if (!male_obs_ct) {
    nmajsums_d[0] = 0.0;
    nmajsums_d[1] = 0.0;
    known_dotprod_d = 0.0;
  }
  double unknown_hethet_d = 0.0;
  if (nonmale_obs_ct) {
    nmajsums_d[0] += 2 * nonmale_nmajsums_d[0];
    nmajsums_d[1] += 2 * nonmale_nmajsums_d[1];
    known_dotprod_d += 2 * nonmale_known_dotprod_d;
    unknown_hethet_d = 2 * nonmale_unknown_hethet_d;
  }
  const double twice_tot_recip = 0.5 / u31tod(weighted_obs_ct);
  double r2;
  const LDErr ld_err = PhasedLD(nmajsums_d, known_dotprod_d, unknown_hethet_d, twice_tot_recip, 0, nullptr, &r2);
  return (ld_err == kLDErrNone)? r2 : -DBL_MAX;
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
  uint32_t founder_main_ct = ctx->founder_ct;
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t founder_nonmale_ct = founder_main_ct - founder_male_ct;
  const uint32_t allow_overlap = ctx->allow_overlap;
  const uint32_t is_x = ctx->is_x;
  if (ctx->is_y) {
    founder_main_ct = founder_male_ct;
  }
  const R2PhaseType phase_type = S_CAST(R2PhaseType, ctx->phase_type);
  const uint32_t load_dosage = ctx->load_dosage;
  const double r2_thresh = ctx->r2_thresh;
  const unsigned char* unpacked_index_variant = &(unpacked_variants[ctx->index_oaidx_offset * unpacked_byte_stride]);
  R2Variant index_r2v;
  if (!is_x) {
    FillR2V(unpacked_index_variant, founder_main_ct, phase_type, load_dosage, &index_r2v);
  } else {
    FillXR2V(unpacked_index_variant, founder_male_ct, founder_nonmale_ct, phase_type, load_dosage, &index_r2v);
  }
  uintptr_t oaidx_base;
  uintptr_t cur_oaidx_bits;
  BitIter1Start(candidate_oabitvec, oaidx_start, &oaidx_base, &cur_oaidx_bits);
  for (; nonindex_rem; --nonindex_rem) {
    const uintptr_t oaidx = BitIter1(candidate_oabitvec, &oaidx_base, &cur_oaidx_bits);
    const uintptr_t oaidx_offset = oaidx - igroup_oaidx_start;
    const unsigned char* unpacked_cur_variant = &(unpacked_variants[oaidx_offset * unpacked_byte_stride]);
    double cur_r2;
    if (!is_x) {
      R2Variant cur_r2v;
      FillR2V(unpacked_cur_variant, founder_main_ct, phase_type, load_dosage, &cur_r2v);
      cur_r2 = ComputeR2(&index_r2v, &cur_r2v, founder_main_ct, phase_type, load_dosage);
    } else {
      R2Variant cur_r2v;
      FillXR2V(unpacked_cur_variant, founder_male_ct, founder_nonmale_ct, phase_type, load_dosage, &cur_r2v);
      cur_r2 = ComputeTwoPartXR2(&index_r2v, &cur_r2v, founder_male_ct, founder_nonmale_ct, phase_type, load_dosage);
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
  uint32_t founder_main_ct = ctx->founder_ct;
  const uint32_t founder_male_ct = ctx->founder_male_ct;
  const uint32_t founder_nonmale_ct = founder_main_ct - founder_male_ct;
  const uint32_t allow_overlap = ctx->allow_overlap;
  const uint32_t is_x = ctx->is_x;
  const uintptr_t* founder_male = nullptr;
  const uintptr_t* founder_nonmale = nullptr;
  const uint32_t* founder_male_cumulative_popcounts = nullptr;
  const uint32_t* founder_nonmale_cumulative_popcounts = nullptr;
  uintptr_t* chrx_workspace = nullptr;
  if (is_x) {
    founder_male = ctx->founder_male;
    founder_nonmale = ctx->founder_nonmale;
    founder_male_cumulative_popcounts = ctx->founder_male_cumulative_popcounts;
    founder_nonmale_cumulative_popcounts = ctx->founder_nonmale_cumulative_popcounts;
    chrx_workspace = ctx->chrx_workspaces[tidx];
    founder_main_ct = ctx->raw_sample_ct;
  } else if (ctx->is_y) {
    founder_main_ct = founder_male_ct;
  }
  const R2PhaseType phase_type = S_CAST(R2PhaseType, ctx->phase_type);
  const uint32_t load_dosage = ctx->load_dosage;
  const double r2_thresh = ctx->r2_thresh;
  unsigned char* unpacked_index_variant = ctx->unpacked_variants;
  R2Variant index_r2v;
  if (!is_x) {
    FillR2V(unpacked_index_variant, founder_main_ct, phase_type, load_dosage, &index_r2v);
  } else {
    FillXR2V(unpacked_index_variant, founder_male_ct, founder_nonmale_ct, phase_type, load_dosage, &index_r2v);
  }
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
    double cur_r2;
    if (!is_x) {
      if (load_dosage) {
        LdUnpackDosage(&pgv, founder_main_ct, phase_type, unpacked_cur_variant);
      } else {
        LdUnpackNondosage(&pgv, founder_main_ct, phase_type, unpacked_cur_variant);
      }
      R2Variant cur_r2v;
      FillR2V(unpacked_cur_variant, founder_main_ct, phase_type, load_dosage, &cur_r2v);
      cur_r2 = ComputeR2(&index_r2v, &cur_r2v, founder_main_ct, phase_type, load_dosage);
    } else {
      if (load_dosage) {
        LdUnpackChrXDosage(&pgv, founder_male, founder_male_cumulative_popcounts, founder_nonmale, founder_nonmale_cumulative_popcounts, founder_main_ct, founder_male_ct, founder_nonmale_ct, phase_type, unpacked_cur_variant, chrx_workspace);
      } else {
        LdUnpackChrXNondosage(&pgv, founder_male, founder_nonmale, founder_main_ct, founder_male_ct, founder_nonmale_ct, phase_type, unpacked_cur_variant, chrx_workspace);
      }
      R2Variant cur_r2v;
      FillXR2V(unpacked_cur_variant, founder_male_ct, founder_nonmale_ct, phase_type, load_dosage, &cur_r2v);
      cur_r2 = ComputeTwoPartXR2(&index_r2v, &cur_r2v, founder_male_ct, founder_nonmale_ct, phase_type, load_dosage);
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
PglErr ClumpReports(const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* founder_info, const uintptr_t* sex_male, const ClumpInfo* clump_ip, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
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
    if (unlikely(!founder_ct)) {
      logerrputs("Error: --clump requires at least 1 founder.  (--make-founders may come in handy\nhere.)\n");
      goto ClumpReports_ret_INCONSISTENT_INPUT;
    } else if (founder_ct > 0x3fffffff) {
      logerrputs("Error: --clump does not support >= 2^30 founders.\n");
      goto ClumpReports_ret_NOT_YET_SUPPORTED;
    }
    uint32_t skipped_variant_ct;
    const uintptr_t* variant_include = StripUnplaced(orig_variant_include, cip, raw_variant_ct, &skipped_variant_ct);
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
    const uint32_t sp2_col = flags & kfClumpColSp2;
    const double ln_p1 = clump_ip->ln_p1;
    const double ln_p2 = sp2_col? clump_ip->ln_p2 : -DBL_MAX;
    double load_ln_pthresh = MAXV(ln_p1, ln_p2);
    if (bin_bound_ct && (load_ln_pthresh < ln_bin_boundaries[bin_bound_ct - 1])) {
      load_ln_pthresh = ln_bin_boundaries[bin_bound_ct - 1];
    }
    ClumpEntry** clump_entries = nullptr;
    uintptr_t* nonsig_arr = nullptr;
    uint32_t allow_overlap = (flags / kfClumpAllowOverlap) & 1;
    if (flags & (kfClumpColTotal | kfClumpColBins | kfClumpColSp2)) {
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
      BigstackEndReset(BigstackEndRoundedDown());
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
    const uint32_t founder_male_ct = PopcountWordsIntersect(founder_info, sex_male, raw_sample_ctl);
    const uint32_t founder_nonmale_ct = founder_ct - founder_male_ct;
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
      } else if (unlikely(founder_male_ct == 0)) {
        // Rather not worry about this case.
        logerrputs("Error: --clump: chrY index variant(s) are present, but the main dataset\ncontains no male founders.\n");
        goto ClumpReports_ret_INCONSISTENT_INPUT;
      }
    }
    // If founder_nonmale_ct == 0, we can just treat as is_y=0, is_haploid=1.
    const uint32_t y_exists = (y_code < UINT32_MAXM1) && founder_nonmale_ct;
    uintptr_t* founder_male = nullptr;
    uintptr_t* founder_nonmale = nullptr;
    uint32_t* founder_male_cumulative_popcounts = nullptr;
    uint32_t* founder_nonmale_cumulative_popcounts = nullptr;
    uint32_t max_sample_ct = founder_ct;
    if (x_exists || y_exists) {
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &founder_male) ||
                   bigstack_alloc_u32(raw_sample_ctl, &founder_male_cumulative_popcounts))) {
        goto ClumpReports_ret_NOMEM;
      }
      BitvecAndCopy(founder_info, sex_male, raw_sample_ctl, founder_male);
      FillCumulativePopcounts(founder_male, raw_sample_ctl, founder_male_cumulative_popcounts);
      if (x_exists) {
        max_sample_ct = raw_sample_ct;
        if (unlikely(bigstack_alloc_w(raw_sample_ctl, &founder_nonmale) ||
                     bigstack_alloc_u32(raw_sample_ctl, &founder_nonmale_cumulative_popcounts))) {
          goto ClumpReports_ret_NOMEM;
        }
        BitvecInvmaskCopy(founder_info, sex_male, raw_sample_ctl, founder_nonmale);
        FillCumulativePopcounts(founder_nonmale, raw_sample_ctl, founder_nonmale_cumulative_popcounts);
      }
    }

    const uintptr_t bitvec_byte_ct = BitCtToVecCt(founder_ct) * kBytesPerVec;
    uintptr_t male_bitvec_byte_ct = 0;
    uintptr_t nonmale_bitvec_byte_ct = 0;
    if (x_exists || y_exists) {
      male_bitvec_byte_ct = BitCtToVecCt(founder_male_ct) * kBytesPerVec;
      if (x_exists) {
        nonmale_bitvec_byte_ct = BitCtToVecCt(founder_nonmale_ct) * kBytesPerVec;
      }
    }
    const uint32_t check_dosage = (pgfip->gflags / kfPgenGlobalDosagePresent) & 1;
    uintptr_t dosagevec_byte_ct = 0;
    uintptr_t male_dosagevec_byte_ct = 0;
    uintptr_t nonmale_dosagevec_byte_ct = 0;
    if (check_dosage) {
      dosagevec_byte_ct = DivUp(founder_ct, kDosagePerVec) * kBytesPerVec;
      if (x_exists || y_exists) {
        male_dosagevec_byte_ct = DivUp(founder_male_ct, kDosagePerVec) * kBytesPerVec;
        if (x_exists) {
          nonmale_dosagevec_byte_ct = DivUp(founder_nonmale_ct, kDosagePerVec) * kBytesPerVec;
        }
      }
    }

    uint32_t calc_thread_ct = MAXV(1, max_thread_ct - 1);
    // no big deal if these are slightly overallocated
    if (unlikely(BIGSTACK_ALLOC_X(PgenReader*, calc_thread_ct, &ctx.pgr_ptrs) ||
                 bigstack_alloc_w(calc_thread_ct + 2, &(ctx.a[0].oaidx_starts)) ||
                 bigstack_alloc_w(calc_thread_ct + 2, &(ctx.a[1].oaidx_starts)) ||
                 bigstack_alloc_wp(calc_thread_ct + 1, &(ctx.ld_idx_found)))) {
      goto ClumpReports_ret_NOMEM;
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

    const uint32_t phased_r2 = !(flags & kfClumpUnphased);
    const uint32_t all_haploid = IsSet(cip->haploid_mask, 0);
    const uint32_t check_phase = phased_r2 && (!all_haploid) && (pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent));
    PgenGlobalFlags effective_gflags = pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
    if (!check_phase) {
      effective_gflags &= kfPgenGlobalDosagePresent;
    }

    if (unlikely(BigstackAllocPgv(max_sample_ct, 0, effective_gflags, &ctx.pgv_base))) {
      goto ClumpReports_ret_NOMEM;
    }
    const uintptr_t pgv_byte_stride = g_bigstack_base - R_CAST(unsigned char*, ctx.pgv_base.genovec);

    uintptr_t x_unpacked_byte_stride = 0;
    uintptr_t nonxy_unpacked_byte_stride;
    if (check_dosage) {
      if (phased_r2) {
        logerrputs("Error: Alpha 5 --clump's handling of dosages when computing phased-r^2 is\nincorrect, and has been disabled.  Either use an alpha 6 or later build for\nthat functionality, or use --clump-unphased to clump on unphased-r^2.\n");
        reterr = kPglRetNotYetSupported;
        goto ClumpReports_ret_1;
      }
      const uintptr_t dosage_trail_byte_ct = LdDosageTrailAlignedByteCt(S_CAST(R2PhaseType, phased_r2));
      nonxy_unpacked_byte_stride = dosagevec_byte_ct * (1 + phased_r2 + check_phase) + bitvec_byte_ct + dosage_trail_byte_ct;
      if (x_exists) {
        x_unpacked_byte_stride = nonmale_dosagevec_byte_ct * (1 + phased_r2 + check_phase) + nonmale_bitvec_byte_ct + male_dosagevec_byte_ct * (1 + phased_r2) + male_bitvec_byte_ct + 2 * dosage_trail_byte_ct;
      }
    } else {
      const uintptr_t nondosage_trail_byte_ct = LdNondosageTrailAlignedByteCt(S_CAST(R2PhaseType, phased_r2));
      nonxy_unpacked_byte_stride = bitvec_byte_ct * (3 + 2 * check_phase) + nondosage_trail_byte_ct;
      if (x_exists) {
        x_unpacked_byte_stride = nonmale_bitvec_byte_ct * (3 + 2 * check_phase) + male_bitvec_byte_ct * 3 + 2 * nondosage_trail_byte_ct;
      }
    }
    const uintptr_t max_unpacked_byte_stride_cachealign = RoundUpPow2(MAXV(nonxy_unpacked_byte_stride, x_unpacked_byte_stride), kCacheline);
    const uintptr_t pgr_struct_alloc = RoundUpPow2(sizeof(PgenReader), kCacheline);

    // Haven't counted PgenReader instance, two unpacked_variants slots, or
    // ld_idx_found slots required by main thread yet.

    // FillGaussianDArr() uses a minimum per-thread job size of ~4 MiB of
    // memory writes.  I'm guessing that is also a reasonable unpacked_variants
    // shard size to aim for.
    const uint32_t min_pgv_per_thread = 1 + 4194303 / pgv_byte_stride;
    uintptr_t chrx_alloc = 0;
    if (x_exists) {
      chrx_alloc = NypCtToCachelineCt(MAXV(founder_male_ct, founder_nonmale_ct)) * kCacheline;
    }
    {
      const uintptr_t min_ld_idx_found_alloc = WordCtToCachelineCt(min_pgv_per_thread + 1) * kCacheline;
      const uintptr_t min_u32_alloc = Int32CtToCachelineCt(min_pgv_per_thread) * kCacheline;
      const uintptr_t phasepresent_alloc = check_phase * min_u32_alloc;
      const uintptr_t dosage_ct_alloc = check_dosage * min_u32_alloc;
      const uintptr_t dphase_ct_alloc = dosage_ct_alloc * check_phase;

      uintptr_t bytes_avail = bigstack_left();
      const uintptr_t more_base_alloc = pgr_struct_alloc + pgr_alloc_cacheline_ct * kCacheline + 2 * max_unpacked_byte_stride_cachealign + min_ld_idx_found_alloc + phasepresent_alloc + dosage_ct_alloc + dphase_ct_alloc + chrx_alloc;
      if (unlikely(bytes_avail < more_base_alloc)) {
        goto ClumpReports_ret_NOMEM;
      }
      bytes_avail -= more_base_alloc;

      const uintptr_t per_thread_target_alloc = 2 * min_pgv_per_thread * pgv_byte_stride + max_unpacked_byte_stride_cachealign + pgr_struct_alloc + pgr_alloc_cacheline_ct * kCacheline + min_ld_idx_found_alloc + chrx_alloc;
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

      if (x_exists) {
        ctx.chrx_workspaces[tidx] = S_CAST(uintptr_t*, bigstack_end_alloc_raw(chrx_alloc));
      }
    }
    uintptr_t* chrx_workspace = x_exists? ctx.chrx_workspaces[calc_thread_ct] : nullptr;

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
    ctx.founder_male = founder_male;
    ctx.founder_male_cumulative_popcounts = founder_male_cumulative_popcounts;
    ctx.founder_nonmale = founder_nonmale;
    ctx.founder_nonmale_cumulative_popcounts = founder_nonmale_cumulative_popcounts;
    ctx.raw_sample_ct = raw_sample_ct;
    ctx.founder_ct = founder_ct;
    ctx.founder_male_ct = founder_male_ct;
    ctx.pgv_byte_stride = pgv_byte_stride;
    ctx.bitvec_byte_ct = bitvec_byte_ct;
    ctx.dosagevec_byte_ct = dosagevec_byte_ct;
    ctx.male_bitvec_byte_ct = male_bitvec_byte_ct;
    ctx.male_dosagevec_byte_ct = male_dosagevec_byte_ct;
    ctx.nonmale_bitvec_byte_ct = nonmale_bitvec_byte_ct;
    ctx.nonmale_dosagevec_byte_ct = nonmale_dosagevec_byte_ct;
    ctx.r2_thresh = clump_ip->r2;
    ctx.allow_overlap = allow_overlap;
    ctx.candidate_oabitvec = candidate_oabitvec;
    ctx.err_info = (~0LLU) << 32;

    unsigned char* lowmem_unpacked_variants = &(g_bigstack_end[(2 + calc_thread_ct) * (-S_CAST(intptr_t, max_unpacked_byte_stride_cachealign))]);
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

      uint64_t unpacked_variant_byte_stride = UnpackedByteStride(&ctx, phase_type, load_dosage);

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
          const uint64_t ext_unpacked_variant_byte_stride = UnpackedByteStride(&ctx, ext_phase_type, ext_load_dosage);
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
        ctx.unpacked_byte_stride = max_unpacked_byte_stride_cachealign;
        ctx.phase_type = phase_type;
        ctx.load_dosage = load_dosage;
        ctx.allele_widx_end = 1 + (allele_idx_last / kBitsPerWord);
        ctx.job_type = kClumpJobLowmemR2;
        STD_SORT(icandidate_ct, u32cmp, &(icandidate_idx_to_rank0_destructive[icandidate_idx_start]));
        const uintptr_t* cur_sample_include = nullptr;
        PgrSampleSubsetIndex pssi;
        uint32_t cur_sample_ct;
        if (is_x) {
          PgrClearSampleSubsetIndex(simple_pgrp, &pssi);
          cur_sample_ct = raw_sample_ct;
        } else {
          const uint32_t* cur_sample_include_cumulative_popcounts;
          if (is_y) {
            cur_sample_include = founder_male;
            cur_sample_include_cumulative_popcounts = founder_male_cumulative_popcounts;
            cur_sample_ct = founder_male_ct;
          } else {
            cur_sample_include = founder_info;
            cur_sample_include_cumulative_popcounts = founder_info_cumulative_popcounts;
            cur_sample_ct = founder_ct;
          }
          PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, simple_pgrp, &pssi);
        }
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
              if (!is_x) {
                if (load_dosage) {
                  if (phase_type == kR2PhaseTypePresent) {
                    reterr = PgrGetInv1Dp(cur_sample_include, pssi, cur_sample_ct, variant_uidx, aidx, simple_pgrp, &pgv);
                  } else {
                    reterr = PgrGetInv1D(cur_sample_include, pssi, cur_sample_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
                  }
                } else {
                  if (phase_type == kR2PhaseTypePresent) {
                    reterr = PgrGetInv1P(cur_sample_include, pssi, cur_sample_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
                  } else {
                    reterr = PgrGetInv1(cur_sample_include, pssi, cur_sample_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec);
                  }
                }
              } else {
                // chrX
                if (load_dosage) {
                  if (phase_type == kR2PhaseTypePresent) {
                    reterr = PgrGetInv1Dp(nullptr, pssi, cur_sample_ct, variant_uidx, aidx, simple_pgrp, &pgv);
                  } else {
                    reterr = PgrGetInv1D(nullptr, pssi, cur_sample_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
                  }
                } else {
                  if (phase_type == kR2PhaseTypePresent) {
                    reterr = PgrGetInv1P(nullptr, pssi, cur_sample_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &pgv.phasepresent_ct);
                  } else {
                    reterr = PgrGetInv1(nullptr, pssi, cur_sample_ct, variant_uidx, aidx, simple_pgrp, pgv.genovec);
                  }
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
                if (!is_x) {
                  if (load_dosage) {
                    LdUnpackDosage(&pgv, cur_sample_ct, phase_type, lowmem_unpacked_variants);
                  } else {
                    LdUnpackNondosage(&pgv, cur_sample_ct, phase_type, lowmem_unpacked_variants);
                  }
                } else {
                  if (load_dosage) {
                    LdUnpackChrXDosage(&pgv, founder_male, founder_male_cumulative_popcounts, founder_nonmale, founder_nonmale_cumulative_popcounts, cur_sample_ct, founder_male_ct, founder_nonmale_ct, phase_type, lowmem_unpacked_variants, chrx_workspace);
                  } else {
                    LdUnpackChrXNondosage(&pgv, founder_male, founder_nonmale, cur_sample_ct, founder_male_ct, founder_nonmale_ct, phase_type, lowmem_unpacked_variants, chrx_workspace);
                  }
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
#ifndef __LP64__
      if (max_vint_byte_ct > 0x7fffffff) {
        goto ClumpReports_ret_NOMEM;
      }
#endif
      if (unlikely(bigstack_alloc_uc(max_vint_byte_ct, &overlap_raw_loadbuf) ||
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
    const uint32_t provref_col = ref_col && ProvrefCol(orig_variant_include, nonref_flags, flags / kfClumpColMaybeprovref, raw_variant_ct, all_nonref);
    const uint32_t a1_col = (flags & kfClumpColA1) || ((flags & kfClumpColMaybeA1) && MultiallelicVariantPresent(orig_variant_include, allele_idx_offsets, orig_variant_ct));
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
    AppendBinaryEoln(&cswritep);
    if (unlikely(Cswrite(&css, &cswritep))) {
      goto ClumpReports_ret_WRITE_FAIL;
    }

    uintptr_t* clump_allele_idxs = overlap_allele_idxs;
    uintptr_t prev_clump_end = 0;
    uint32_t index_allele_ct = 2;
    uint32_t index_file_idx1 = 1;
    uint32_t file_idx1 = 1;
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
      if (chr_col) {
        uint32_t chr_idx = GetVariantChr(cip, index_variant_uidx);
        cswritep = chrtoa(cip, chr_idx, cswritep);
        *cswritep++ = '\t';
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

#ifdef __cplusplus
}  // namespace plink2
#endif
