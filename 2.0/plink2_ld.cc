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


#include "include/plink2_stats.h"
#include "plink2_ld.h"

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
void LdPruneNextSubcontig(const uintptr_t* variant_include, const uint32_t* variant_bps, const uint32_t* subcontig_info, const uint32_t* subcontig_thread_assignments, uint32_t x_start, uint32_t x_len, uint32_t y_start, uint32_t y_len, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t prune_window_size, uint32_t thread_idx, uint32_t* subcontig_idx_ptr, uint32_t* subcontig_end_tvidx_ptr, uint32_t* next_window_end_tvidx_ptr, uint32_t* is_x_ptr, uint32_t* is_y_ptr, uint32_t* cur_founder_ct_ptr, uint32_t* cur_founder_ctaw_ptr, uint32_t* cur_founder_ctl_ptr, uintptr_t* entire_variant_buf_word_ct_ptr, uint32_t* variant_uidx_winstart_ptr, uint32_t* variant_uidx_winend_ptr) {
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
    // assumed that cur_first_nm initialized to first_vaggs[0], cur_first_sum
    // initialized to first_vaggs[1], cur_first_ssq to first_vaggs[2]
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

typedef struct IndepPairwiseCtxStruct {
  const uint32_t* subcontig_info;
  const uint32_t* subcontig_thread_assignments;
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const AlleleCode* maj_alleles;
  const double* all_allele_freqs;
  const uint32_t* variant_bps;
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
  uintptr_t** raw_tgenovecs[2];
} IndepPairwiseCtx;

THREAD_FUNC_DECL IndepPairwiseThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  IndepPairwiseCtx* ctx = S_CAST(IndepPairwiseCtx*, arg->sharedp->context);

  const uint32_t* subcontig_info = ctx->subcontig_info;
  const uint32_t* subcontig_thread_assignments = ctx->subcontig_thread_assignments;
  const uintptr_t* variant_include = ctx->variant_include;
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
  uint32_t* first_unchecked_tvidx = ctx->first_unchecked_tvidx[tidx];

  uint32_t subcontig_end_tvidx = 0;
  uint32_t subcontig_idx = UINT32_MAX;  // deliberate overflow
  uint32_t window_start_tvidx = 0;
  uint32_t next_window_end_tvidx = 0;
  uint32_t write_slot_idx = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t cur_window_size = 0;
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
        LdPruneNextSubcontig(variant_include, variant_bps, subcontig_info, subcontig_thread_assignments, x_start, x_len, y_start, y_len, founder_ct, founder_male_ct, prune_window_size, tidx, &subcontig_idx, &subcontig_end_tvidx, &next_window_end_tvidx, &is_x, &is_y, &cur_founder_ct, &cur_founder_ctaw, &cur_founder_ctl, &entire_variant_buf_word_ct, &variant_uidx_winstart, &variant_uidx_winend);
        BitIter1Start(variant_include, variant_uidx_winstart, &variant_uidx_base, &variant_include_bits);
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
        first_unchecked_tvidx[write_slot_idx] = cur_tvidx + 1;
      }
      SetBit(write_slot_idx, occupied_window_slots);
      winpos_to_slot_idx[cur_window_size++] = write_slot_idx;
      // are we at the end of a window?
      if (++cur_tvidx == next_window_end_tvidx) {
        // possible for cur_window_size == 1, if all variants at the end of the
        // previous window were pruned
        uint32_t cur_removed_ct = PopcountWords(cur_window_removed, BitCtToWordCt(cur_window_size));
        uint32_t prev_removed_ct;
        do {
          prev_removed_ct = cur_removed_ct;
          uint32_t first_winpos = 0;
          // const uint32_t debug_print = (!IsSet(cur_window_removed, 0)) && (tvidxs[winpos_to_slot_idx[0]] == 0);
          while (1) {
            // can't use BitIter0 since we care about changes in this loop to
            // cur_window_removed
            first_winpos = AdvTo0Bit(cur_window_removed, first_winpos);
            // can assume empty trailing bit for cur_window_removed
            if (first_winpos == cur_window_size) {
              break;
            }
            uint32_t first_slot_idx = winpos_to_slot_idx[first_winpos];
            const uint32_t cur_first_unchecked_tvidx = first_unchecked_tvidx[first_slot_idx];
            // safe to use BitIter0 for second_winpos, though
            uintptr_t second_winpos_base;
            uintptr_t cur_window_removed_inv_bits;
            BitIter0Start(cur_window_removed, first_winpos + 1, &second_winpos_base, &cur_window_removed_inv_bits);
            while (1) {
              uint32_t second_winpos = BitIter0(cur_window_removed, &second_winpos_base, &cur_window_removed_inv_bits);
              if (second_winpos == cur_window_size) {
                break;
              }
              uint32_t second_slot_idx = winpos_to_slot_idx[second_winpos];
              if (tvidxs[second_slot_idx] >= cur_first_unchecked_tvidx) {
                uintptr_t* first_genobufs = &(genobufs[first_slot_idx * entire_variant_buf_word_ct]);
                const uint32_t first_nm_ct = vaggs[first_slot_idx].nm_ct;
                const int32_t first_sum = vaggs[first_slot_idx].sum;
                const uint32_t first_ssq = vaggs[first_slot_idx].ssq;
                while (1) {
                  uintptr_t* second_genobufs = &(genobufs[second_slot_idx * entire_variant_buf_word_ct]);
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
                      /*
                      if (debug_print) {
                        printf("removing %u, keeping %u, freqs %g/%g, r2 = %g\n", tvidxs[first_slot_idx], tvidxs[second_slot_idx], cur_maj_freqs[first_slot_idx], cur_maj_freqs[second_slot_idx], cov12 * cov12 / (variance1 * variance2));
                      }
                      */
                      SetBit(first_winpos, cur_window_removed);
                      SetBit(tvidxs[first_slot_idx], removed_variants_write);
                    } else {
                      /*
                      if (debug_print) {
                        printf("removing %u (second), keeping %u, freqs %g/%g, r2 = %g\n", tvidxs[second_slot_idx], tvidxs[first_slot_idx], cur_maj_freqs[second_slot_idx], cur_maj_freqs[first_slot_idx], cov12 * cov12 / (variance1 * variance2));
                      }
                      */
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
                break;
              }
            }
            ++first_winpos;
          }
          cur_removed_ct = PopcountWords(cur_window_removed, BitCtToWordCt(cur_window_size));
        } while (cur_removed_ct > prev_removed_ct);
        const uint32_t prev_window_size = cur_window_size;
        LdPruneNextWindow(variant_include, variant_bps, tvidxs, cur_window_removed, prune_window_size, window_incr, window_maxl, subcontig_end_tvidx, &cur_window_size, &window_start_tvidx, &variant_uidx_winstart, &next_window_end_tvidx, &variant_uidx_winend, occupied_window_slots, winpos_to_slot_idx);
        // clear bits here since we set cur_window_removed bits during loading
        // process in monomorphic case
        ZeroWArr(BitCtToWordCt(prev_window_size), cur_window_removed);
        write_slot_idx = 0;
      }
    }
    parity = 1 - parity;
    tvidx_start = tvidx_stop;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr IndepPairwise(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uint32_t* founder_info_cumulative_popcounts, const uintptr_t* founder_nonmale, const uintptr_t* founder_male, const LdInfo* ldip, const uint32_t* subcontig_info, const uint32_t* subcontig_thread_assignments, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t subcontig_ct, uintptr_t window_max, uint32_t calc_thread_ct, uint32_t max_load, PgenReader* simple_pgrp, uintptr_t* removed_variants_collapsed) {
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
    //     sizeof(intptr_t) for raw genotype data (g_raw_tgenovecs)
    // - tvidx_batch_size * sizeof(double) for g_maj_freqs
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
    if (unlikely(
            bigstack_alloc_w(NypCtToWordCt(raw_sample_ct), &tmp_genovec) ||
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
            bigstack_alloc_u32p(calc_thread_ct, &ctx.first_unchecked_tvidx) ||
            bigstack_alloc_wp(calc_thread_ct, &(ctx.raw_tgenovecs[0])) ||
            bigstack_alloc_wp(calc_thread_ct, &(ctx.raw_tgenovecs[1])))) {
      goto IndepPairwise_ret_NOMEM;
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

    // three of these
    const uintptr_t window_int32_alloc = RoundUpPow2(window_max * sizeof(int32_t), kCacheline);

    const uintptr_t thread_alloc_base = genobuf_alloc + occupied_window_slots_alloc + cur_window_removed_alloc + cur_maj_freqs_alloc + removed_variants_write_alloc + 2 * vaggs_alloc + 3 * window_int32_alloc;

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
      ctx.first_unchecked_tvidx[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(window_int32_alloc));
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
    ctx.founder_ct = founder_ct;
    ctx.founder_male_ct = founder_male_ct;
    ctx.prune_window_size = ldip->prune_window_size;
    ctx.window_maxl = window_maxl;
    ctx.prune_ld_thresh = ldip->prune_last_param * (1 + kSmallEpsilon);
    ctx.window_incr = ldip->prune_window_incr;
    ctx.cur_batch_size = tvidx_batch_size;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto IndepPairwise_ret_NOMEM;
    }
    SetThreadFuncAndData(IndepPairwiseThread, &ctx, &tg);

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
              goto IndepPairwise_ret_PGR_FAIL;
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
  IndepPairwise_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  IndepPairwise_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  IndepPairwise_ret_NOT_YET_SUPPORTED:
    reterr = kPglRetNotYetSupported;
    break;
  }
  CleanupThreads(&tg);
  // caller will free memory
  return reterr;
}

PglErr IndepPairphase() {
  logerrputs("Error: --indep-pairphase is currently under development.\n");
  return kPglRetNotYetSupported;
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
void minheap64_delete_root(uint64_t* minheap64_preroot, uint32_t* heap_size_ptr) {
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
    const uintptr_t idxp1 = CountSortedSmallerU64(minheap64, thread_ct, max_load_shifted - cur_tagged_weight);
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
      const uintptr_t insert_pt = CountSortedSmallerU64(minheap64, thread_ct, new_entry);
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

PglErr LdPrune(const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uintptr_t* sex_male, const LdInfo* ldip, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  // common initialization between --indep-pairwise and --indep-pairphase
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t is_pairphase = (ldip->prune_flags / kfLdPrunePairphase) & 1;
    if (founder_ct < 2) {
      logerrprintfww("Error: --indep-pair%s requires at least two founders. (PLINK 1.9 --make-founders may come in handy here.)\n", is_pairphase? "phase" : "wise");
      goto LdPrune_ret_INCONSISTENT_INPUT;
    }
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
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uintptr_t* variant_include;
    if (skipped_variant_ct) {
      uintptr_t* new_variant_include;
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, &new_variant_include))) {
        goto LdPrune_ret_NOMEM;
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
      variant_include = new_variant_include;
      variant_ct -= skipped_variant_ct;
      logprintf("--indep-pair%s: Ignoring %u chromosome 0 variant%s.\n", is_pairphase? "phase" : "wise", skipped_variant_ct, (skipped_variant_ct == 1)? "" : "s");
    } else {
      variant_include = orig_variant_include;
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

    {
      uint32_t dup_found;
      reterr = CheckIdUniqueness(g_bigstack_base, g_bigstack_end, variant_include, variant_ids, variant_ct, max_thread_ct, &dup_found);
      if (unlikely(reterr)) {
        goto LdPrune_ret_1;
      }
      if (unlikely(dup_found)) {
        logerrprintfww("Error: --indep-pair%s requires unique variant IDs. (--set-all-var-ids and/or --rm-dup may help.)\n", is_pairphase? "phase" : "wise");
        goto LdPrune_ret_INCONSISTENT_INPUT;
      }
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
    const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
    uint32_t* founder_info_cumulative_popcounts;
    uintptr_t* founder_nonmale_collapsed;
    uintptr_t* founder_male_collapsed;
    uintptr_t* removed_variants_collapsed;
    uint32_t* subcontig_thread_assignments;
    if (unlikely(
            bigstack_alloc_u32(raw_sample_ctl, &founder_info_cumulative_popcounts) ||
            bigstack_alloc_w(founder_ctl, &founder_nonmale_collapsed) ||
            bigstack_alloc_w(founder_ctl, &founder_male_collapsed) ||
            bigstack_calloc_w(variant_ctl, &removed_variants_collapsed) ||
            bigstack_alloc_u32(subcontig_ct, &subcontig_thread_assignments))) {
      goto LdPrune_ret_NOMEM;
    }
    FillCumulativePopcounts(founder_info, raw_sample_ctl, founder_info_cumulative_popcounts);
    CopyBitarrSubset(sex_male, founder_info, founder_ct, founder_male_collapsed);
    AlignedBitarrInvertCopy(founder_male_collapsed, founder_ct, founder_nonmale_collapsed);
    uint32_t* subcontig_weights;
    if (unlikely(bigstack_end_alloc_u32(subcontig_ct, &subcontig_weights))) {
      goto LdPrune_ret_NOMEM;
    }

    // initial window_max-based memory requirement estimate
    if (is_pairphase) {
      // todo
      // this should allow a mix of phased and unphased calls; see --ld
      // implementation below.
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
      //   vstats, nonmale_vstats: thread_ct * window_max * 3 * int32
      //   winpos_to_slot_idx, tvidxs, first_unchecked_vidx: window_max * 3 *
      //     int32
      uintptr_t per_thread_alloc = RoundUpPow2(window_max * entire_variant_buf_word_ct * sizeof(intptr_t), kCacheline) + 2 * RoundUpPow2((1 + window_max / kBitsPerWord) * sizeof(intptr_t), kCacheline) + RoundUpPow2(window_max * sizeof(double), kCacheline) + 2 * RoundUpPow2(window_max * (3 * sizeof(int32_t)), kCacheline) + 3 * RoundUpPow2(window_max * sizeof(int32_t), kCacheline);
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
      reterr = IndepPairphase();
    } else {
      reterr = IndepPairwise(variant_include, cip, variant_bps, allele_idx_offsets, maj_alleles, allele_freqs, founder_info, founder_info_cumulative_popcounts, founder_nonmale_collapsed, founder_male_collapsed, ldip, subcontig_info, subcontig_thread_assignments, raw_sample_ct, founder_ct, founder_male_ct, subcontig_ct, window_max, max_thread_ct, max_load, simple_pgrp, removed_variants_collapsed);
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
  const __m256i all_n32768 = {0x8000800080008000LLU, 0x8000800080008000LLU, 0x8000800080008000LLU, 0x8000800080008000LLU};
  const __m256i all_n16384 = {0xc000c000c000c000LLU, 0xc000c000c000c000LLU, 0xc000c000c000c000LLU, 0xc000c000c000c000LLU};
#    else
  const __m256i all_n32768 = {-0x7fff7fff7fff8000LL, -0x7fff7fff7fff8000LL, -0x7fff7fff7fff8000LL, -0x7fff7fff7fff8000LL};
  const __m256i all_n16384 = {-0x3fff3fff3fff4000LL, -0x3fff3fff3fff4000LL, -0x3fff3fff3fff4000LL, -0x3fff3fff3fff4000LL};
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
  const __m256i m16 = {kMask0000FFFF, kMask0000FFFF, kMask0000FFFF, kMask0000FFFF};
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
  const __m256i m16 = {kMask0000FFFF, kMask0000FFFF, kMask0000FFFF, kMask0000FFFF};
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
  const __m256i m16 = {kMask0000FFFF, kMask0000FFFF, kMask0000FFFF, kMask0000FFFF};
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
  const __m256i m16 = {kMask0000FFFF, kMask0000FFFF, kMask0000FFFF, kMask0000FFFF};
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
  const __m256i m16 = {kMask0000FFFF, kMask0000FFFF, kMask0000FFFF, kMask0000FFFF};
  const __m256i all_4096 = {0x1000100010001000LLU, 0x1000100010001000LLU, 0x1000100010001000LLU, 0x1000100010001000LLU};
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
#  else  // !USE_AVX2
void FillDosageUhet(const Dosage* dosage_vec, uint32_t dosagev_ct, Dosage* dosage_uhet) {
  const __m128i* dosage_vvec_iter = R_CAST(const __m128i*, dosage_vec);
#    if defined(__APPLE__) && ((!defined(__cplusplus)) || (__cplusplus < 201103L))
  const __m128i all_n32768 = {0x8000800080008000LLU, 0x8000800080008000LLU};
  const __m128i all_n16384 = {0xc000c000c000c000LLU, 0xc000c000c000c000LLU};
#    else
  const __m128i all_n32768 = {-0x7fff7fff7fff8000LL, -0x7fff7fff7fff8000LL};
  const __m128i all_n16384 = {-0x3fff3fff3fff4000LL, -0x3fff3fff3fff4000LL};
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
  const __m128i m16 = {kMask0000FFFF, kMask0000FFFF};
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
  const __m128i m16 = {kMask0000FFFF, kMask0000FFFF};
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
  const __m128i m16 = {kMask0000FFFF, kMask0000FFFF};
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
  const __m128i m16 = {kMask0000FFFF, kMask0000FFFF};
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
  const __m128i m16 = {kMask0000FFFF, kMask0000FFFF};
  const __m128i all_4096 = {0x1000100010001000LLU, 0x1000100010001000LLU};
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
#endif

void DosagePhaseinfoPatch(const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const uintptr_t* dosage_present, const Dosage* dosage_vec, const uintptr_t* dphase_present, uint32_t sample_ct, Dosage* dosage_uhet, SDosage* dphase_delta) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    uintptr_t phasepresent_nodphase_word = phasepresent[widx] & (~dphase_present[widx]);
    if (phasepresent_nodphase_word) {
      const uintptr_t phaseinfo_word = phaseinfo[widx];
      const uintptr_t dosage_present_word = dosage_present[widx];
      const uint32_t sample_idx_offset = widx * kBitsPerWord;
      do {
        const uint32_t sample_idx_lowbits = ctzw(phasepresent_nodphase_word);
        const uint32_t sample_idx = sample_idx_offset + sample_idx_lowbits;
        // should compile to blsi
        const uintptr_t shifted_bit = phasepresent_nodphase_word & (-phasepresent_nodphase_word);
        int32_t cur_diff = kDosageMid;
        if (dosage_present_word & shifted_bit) {
          cur_diff = DosageHomdist(dosage_vec[sample_idx]);
        }
        // todo: verify the compiler optimizes this well
        if (!(phaseinfo_word & shifted_bit)) {
          cur_diff = -cur_diff;
        }
        dosage_uhet[sample_idx] = 0;
        dphase_delta[sample_idx] = cur_diff;
        phasepresent_nodphase_word ^= shifted_bit;
      } while (phasepresent_nodphase_word);
    }
  }
}

uint32_t DosagePhasedR2Prod(const Dosage* dosage_vec0, const uintptr_t* nm_bitvec0, const Dosage* dosage_vec1, const uintptr_t* nm_bitvec1, uint32_t sample_ct, uint32_t nm_ct0, uint32_t nm_ct1, uint64_t* __restrict nmaj_dosages, uint64_t* __restrict dosageprod_ptr) {
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
  freq11 += freq11_incr;
  freq22 += freq11_incr;
  freq12 += half_hethet_share - freq11_incr;
  freq21 += half_hethet_share - freq11_incr;
  const double cross_sum = freq11 * freq22 + freq12 * freq21;
  double lnlike = 0.0;
  if (cross_sum != 0.0) {
    lnlike = half_hethet_share * log(cross_sum);
  }
  if (freq11 != 0.0) {
    lnlike += freq11 * log(freq11);
  }
  if (freq12 != 0.0) {
    lnlike += freq12 * log(freq12);
  }
  if (freq21 != 0.0) {
    lnlike += freq21 * log(freq21);
  }
  if (freq22 != 0.0) {
    lnlike += freq22 * log(freq22);
  }
  return lnlike;
}

PglErr LdConsole(const uintptr_t* variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const char* const* allele_storage, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const LdInfo* ldip, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, PgenReader* simple_pgrp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    if (!founder_ct) {
      logerrputs("Error: --ld requires founders.  (PLINK 1.9 --make-founders\nmay come in handy\nhere.)\n");
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
    if (unlikely(
            bigstack_alloc_u32(founder_ctl, &founder_info_cumulative_popcounts) ||
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
      if (unlikely(
              bigstack_alloc_w(founder_ctaw, &sex_male_collapsed) ||
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
        goto LdConsole_ret_PGR_FAIL;
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
          SetMaleHetMissing(sex_male_collapsed_interleaved, founder_ctv, pgvp->genovec);
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
      if (unlikely(
              bigstack_alloc_w(founder_ctaw, &one_bitvecs[0]) ||
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
      if (unlikely(
              bigstack_alloc_dosage(founder_ct, &dosage_vecs[0]) ||
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
        PopulateDenseDosage(pgvs[var_idx].genovec, pgvs[var_idx].dosage_present, pgvs[var_idx].dosage_main, founder_ct, pgvs[var_idx].dosage_ct, dosage_vecs[var_idx]);
        nmaj_dosages[var_idx] = DenseDosageSum(dosage_vecs[var_idx], founder_dosagev_ct);
        FillDosageUhet(dosage_vecs[var_idx], founder_dosagev_ct, dosage_uhets[var_idx]);
        GenoarrToNonmissingnessUnsafe(pgvs[var_idx].genovec, founder_ct, nm_bitvecs[var_idx]);
        ZeroTrailingBits(founder_ct, nm_bitvecs[var_idx]);
        BitvecOr(pgvs[var_idx].dosage_present, founder_ctl, nm_bitvecs[var_idx]);
        nm_cts[var_idx] = PopcountWords(nm_bitvecs[var_idx], founder_ctl);
      }
      SDosage* main_dphase_deltas[2];
      main_dphase_deltas[0] = nullptr;
      main_dphase_deltas[1] = nullptr;
      if (use_phase) {
        if (unlikely(
                bigstack_alloc_dphase(founder_ct, &main_dphase_deltas[0]) ||
                bigstack_alloc_dphase(founder_ct, &main_dphase_deltas[1]))) {
          goto LdConsole_ret_NOMEM;
        }
        for (uint32_t var_idx = 0; var_idx != 2; ++var_idx) {
          // this should probably be a PopulateDenseDphase() function in
          // plink2_common
          ZeroDphaseArr(founder_dosagev_ct * kDosagePerVec, main_dphase_deltas[var_idx]);
          const SDosage* read_dphase_delta = pgvs[var_idx].dphase_delta;
          uintptr_t* cur_dphase_present = pgvs[var_idx].dphase_present;
          const Dosage* cur_dosage_vec = dosage_vecs[var_idx];
          Dosage* cur_dosage_uhet = dosage_uhets[var_idx];
          SDosage* cur_dphase_delta = main_dphase_deltas[var_idx];
          const uint32_t cur_dphase_ct = pgvs[var_idx].dphase_ct;
          if (cur_dphase_ct) {
            uintptr_t sample_uidx_base = 0;
            uintptr_t cur_bits = cur_dphase_present[0];
            for (uint32_t dphase_idx = 0; dphase_idx != cur_dphase_ct; ++dphase_idx) {
              const uintptr_t sample_uidx = BitIter1(cur_dphase_present, &sample_uidx_base, &cur_bits);
              const uint32_t dosage_int = cur_dosage_vec[sample_uidx];
              const int32_t dphase_delta_val = read_dphase_delta[dphase_idx];
              cur_dosage_uhet[sample_uidx] = DosageHomdist(dosage_int) - abs_i32(dphase_delta_val);
              cur_dphase_delta[sample_uidx] = dphase_delta_val;
            }
          } else {
            ZeroWArr(founder_ctl, cur_dphase_present);
          }
          if (pgvs[var_idx].phasepresent_ct) {
            DosagePhaseinfoPatch(pgvs[var_idx].phasepresent, pgvs[var_idx].phaseinfo, pgvs[var_idx].dosage_present, cur_dosage_vec, cur_dphase_present, founder_ct, cur_dosage_uhet, cur_dphase_delta);
          }
        }
      }
      const uint64_t orig_nmaj_dosage1 = nmaj_dosages[1];
      uint64_t dosageprod;
      valid_obs_ct = DosagePhasedR2Prod(dosage_vecs[0], nm_bitvecs[0], dosage_vecs[1], nm_bitvecs[1], founder_ct, nm_cts[0], nm_cts[1], nmaj_dosages, &dosageprod);
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
      nmajsums_d[0] = S_CAST(int64_t, nmaj_dosages[0]) * kRecipDosageMid;
      nmajsums_d[1] = S_CAST(int64_t, nmaj_dosages[1]) * kRecipDosageMid;
      known_dotprod_d = S_CAST(int64_t, dosageprod - uhethet_dosageprod) * (kRecipDosageMidSq * 0.5);
      unknown_hethet_d = S_CAST(int64_t, uhethet_dosageprod) * kRecipDosageMidSq;
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
        valid_x_male_ct = DosagePhasedR2Prod(dosage_vecs[0], nm_bitvecs[0], dosage_vecs[1], nm_bitvecs[1], founder_ct, x_male_nm_ct0, nm_cts[1], x_male_nmaj_dosages, &x_male_dosageprod);
        if (!ignore_hethet) {
          BitvecInvmask(R_CAST(uintptr_t*, x_male_dosage_invmask), founder_dosagev_ct * kWordsPerVec, R_CAST(uintptr_t*, dosage_uhets[0]));
          const uint64_t invalid_uhethet_dosageprod = DosageUnsignedNomissDotprod(dosage_uhets[0], dosage_uhets[1], founder_dosagev_ct);
          unknown_hethet_d -= S_CAST(int64_t, invalid_uhethet_dosageprod) * kRecipDosageMidSq;
        }
        x_male_nmajsums_d[0] = S_CAST(int64_t, x_male_nmaj_dosages[0]) * kRecipDosageMid;
        x_male_nmajsums_d[1] = S_CAST(int64_t, x_male_nmaj_dosages[1]) * kRecipDosageMid;
        x_male_known_dotprod_d = S_CAST(int64_t, x_male_dosageprod) * (kRecipDosageMidSq * 0.5);
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

    // known-diplotype dosages (sum is 2 * (valid_obs_d - unknown_hethet_d)):
    //   var0  var1
    //     0  -  0 : 2 * valid_obs_d - majsums[0] - majsums[1] + known_dotprod
    //     1  -  0 : majsums[0] - known_dotprod - unknown_hethet_d
    //     0  -  1 : majsums[1] - known_dotprod - unknown_hethet_d
    //     1  -  1 : known_dotprod
    double freq_majmaj = 1.0 - (nmajsums_d[0] + nmajsums_d[1] - known_dotprod_d) * twice_tot_recip;
    double freq_majmin = (nmajsums_d[1] - known_dotprod_d - unknown_hethet_d) * twice_tot_recip;
    double freq_minmaj = (nmajsums_d[0] - known_dotprod_d - unknown_hethet_d) * twice_tot_recip;
    double freq_minmin = known_dotprod_d * twice_tot_recip;
    const double half_unphased_hethet_share = unknown_hethet_d * twice_tot_recip;
    const double freq_majx = freq_majmaj + freq_majmin + half_unphased_hethet_share;
    const double freq_minx = 1.0 - freq_majx;
    const double freq_xmaj = freq_majmaj + freq_minmaj + half_unphased_hethet_share;
    const double freq_xmin = 1.0 - freq_xmaj;
    // frequency of ~2^{-46} is actually possible with dosages and 2 billion
    // samples, so set this threshold at 2^{-47}
    if ((freq_majx < (kSmallEpsilon * 0.125)) || (freq_minx < (kSmallEpsilon * 0.125))) {
      logerrprintfww("Warning: Skipping --ld since %s is monomorphic across all valid observations.\n", ld_console_varids[0]);
      goto LdConsole_ret_1;
    }
    if ((freq_xmaj < (kSmallEpsilon * 0.125)) || (freq_xmin < (kSmallEpsilon * 0.125))) {
      logerrprintfww("Warning: Skipping --ld since %s is monomorphic across all valid observations.\n", ld_console_varids[1]);
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

    uint32_t cubic_sol_ct = 0;
    uint32_t first_relevant_sol_idx = 0;
    uint32_t best_lnlike_mask = 0;
    STD_ARRAY_DECL(double, 3, cubic_sols);
    if (half_unphased_hethet_share != 0.0) {
      // detect degenerate cases to avoid e-17 ugliness
      if ((freq_majmaj * freq_minmin != 0.0) || (freq_majmin * freq_minmaj != 0.0)) {
        // (f11 + x)(f22 + x)(K - x) = x(f12 + K - x)(f21 + K - x)
        // (x - K)(x + f11)(x + f22) + x(x - K - f12)(x - K - f21) = 0
        //   x^3 + (f11 + f22 - K)x^2 + (f11*f22 - K*f11 - K*f22)x
        // - K*f11*f22 + x^3 - (2K + f12 + f21)x^2 + (K + f12)(K + f21)x = 0
        cubic_sol_ct = CubicRealRoots(0.5 * (freq_majmaj + freq_minmin - freq_majmin - freq_minmaj - 3 * half_unphased_hethet_share), 0.5 * (freq_majmaj * freq_minmin + freq_majmin * freq_minmaj + half_unphased_hethet_share * (freq_majmin + freq_minmaj - freq_majmaj - freq_minmin + half_unphased_hethet_share)), -0.5 * half_unphased_hethet_share * freq_majmaj * freq_minmin, cubic_sols);
        if (cubic_sol_ct > 1) {
          while (cubic_sols[cubic_sol_ct - 1] > half_unphased_hethet_share + kSmallishEpsilon) {
            --cubic_sol_ct;
          }
          if (cubic_sols[cubic_sol_ct - 1] > half_unphased_hethet_share - kSmallishEpsilon) {
            cubic_sols[cubic_sol_ct - 1] = half_unphased_hethet_share;
          }
          // todo: document why this is safe (or fix if it isn't)
          while (cubic_sols[first_relevant_sol_idx] < -kSmallishEpsilon) {
            ++first_relevant_sol_idx;
          }
          if (cubic_sols[first_relevant_sol_idx] < kSmallishEpsilon) {
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
            AlignedBitarrInvert(founder_ctl, nosex_collapsed);
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
    logputs("\n");

    for (uint32_t sol_idx = first_relevant_sol_idx; sol_idx < cubic_sol_ct; ++sol_idx) {
      if (cubic_sol_ct - first_relevant_sol_idx > 1) {
        write_iter = strcpya_k(g_logbuf, "Solution #");
        write_iter = u32toa(sol_idx + 1 - first_relevant_sol_idx, write_iter);
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
  LdConsole_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
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

#ifdef __cplusplus
}  // namespace plink2
#endif
