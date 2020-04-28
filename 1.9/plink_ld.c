// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink_common.h"

#include <stddef.h>
#include "plink_assoc.h"
#include "plink_glm.h"
#include "plink_ld.h"
#include "plink_stats.h"
#include "pigz.h"

#define MULTIPLEX_LD 1920
#define MULTIPLEX_2LD (MULTIPLEX_LD * 2)

void ld_epi_init(Ld_info* ldip, Epi_info* epi_ip, Clump_info* clump_ip) {
  ldip->modifier = 0;
  ldip->prune_window_size = 0;
  ldip->prune_window_incr = 0;
  ldip->prune_last_param = 0.0;
  ldip->window_size = 0xffffffffU;
  ldip->window_bp = 0xffffffffU;
  ldip->window_cm = -1;
  ldip->window_r2 = 0.2;
  ldip->blocks_max_bp = 0xffffffffU;
  ldip->blocks_min_maf = 0.05;
  ldip->blocks_strong_lowci_outer = 71;
  ldip->blocks_strong_lowci = 72;
  ldip->blocks_strong_highci = 97;
  ldip->blocks_recomb_highci = 89;
  ldip->blocks_inform_frac = 0.95;
  ldip->flipscan_window_size = 10;
  ldip->flipscan_window_bp = 0xffffffffU;
  ldip->flipscan_thresh = 0.5;
  ldip->show_tags_bp = 250000;
  ldip->show_tags_r2 = 0.8;
  ldip->snpstr = nullptr;
  ldip->show_tags_fname = nullptr;
  range_list_init(&(ldip->snps_rl));
  epi_ip->modifier = 0;
  epi_ip->case_only_gap = 1000000;
  epi_ip->epi1 = 0.0;
  epi_ip->epi2 = 0.01;
  epi_ip->je_cellmin = 5;
  epi_ip->ld_mkr1 = nullptr;
  epi_ip->ld_mkr2 = nullptr;
  epi_ip->twolocus_mkr1 = nullptr;
  epi_ip->twolocus_mkr2 = nullptr;
  epi_ip->summary_merge_prefix = nullptr;
  clump_ip->modifier = 0;
  clump_ip->fname_ct = 0;
  clump_ip->bp_radius = 249999;
  clump_ip->range_border = 0;
  clump_ip->fnames_flattened = nullptr;
  clump_ip->annotate_flattened = nullptr;
  clump_ip->snpfield_search_order = nullptr;
  clump_ip->pfield_search_order = nullptr;
  clump_ip->range_fname = nullptr;
  clump_ip->p1 = 1e-4;
  clump_ip->p2 = 1e-2;
  clump_ip->r2 = 0.5;
}

void ld_epi_cleanup(Ld_info* ldip, Epi_info* epi_ip, Clump_info* clump_ip) {
  free_cond(ldip->snpstr);
  free_cond(ldip->show_tags_fname);
  free_range_list(&(ldip->snps_rl));
  free_cond(epi_ip->ld_mkr1);
  free_cond(epi_ip->ld_mkr2);
  free_cond(epi_ip->twolocus_mkr1);
  free_cond(epi_ip->twolocus_mkr2);
  free_cond(epi_ip->summary_merge_prefix);
  free_cond(clump_ip->fnames_flattened);
  free_cond(clump_ip->annotate_flattened);
  free_cond(clump_ip->snpfield_search_order);
  free_cond(clump_ip->pfield_search_order);
  free_cond(clump_ip->range_fname);
}

#ifdef __LP64__
static inline void ld_dot_prod_batch(__m128i* vec1, __m128i* vec2, __m128i* mask1, __m128i* mask2, int32_t* return_vals, uint32_t iters) {
  // Main routine for computation of \sum_i^M (x_i - \mu_x)(y_i - \mu_y), where
  // x_i, y_i \in \{-1, 0, 1\}, but there are missing values.
  //
  //
  // We decompose this sum into
  //   \sum_i x_iy_i - \mu_y\sum_i x_i - \mu_x\sum_i y_i +
  //   (M - # missing)\mu_x\mu_y.
  // *Without* missing values, this can be handled very cleanly.  The last
  // three terms can all be precomputed, and \sum_i x_iy_i can be handled in a
  // manner very similar to bitwise Hamming distance.  This is several times as
  // fast as the lookup tables used for relationship matrices.
  //
  // Unfortunately, when missing values are present,
  // \mu_y\sum_{i: nonmissing from y} x_i and
  // \mu_x\sum_{i: nonmissing from x} y_i must also be evaluated (and, in
  // practice, \mu_y\sum_{i: nonmissing from y} x_i^2 and
  // \mu_x\sum_{i: nonmissing from x} y_i^2 should be determined here as well);
  // this removes much of the speed advantage, and the best applications of the
  // underlying ternary dot product algorithm used here lie elsewhere.
  // Nevertheless, it is still faster, so we use it.
  // (possible todo: accelerated function when there really are no missing
  // values, similar to what is now done for --fast-epistasis)
  //
  //
  // Input:
  // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
  // * mask1 and mask2 mask out missing values (i.e. 00 for missing, 11 for
  //   nonmissing).
  // * return_vals provides space for return values.
  // * iters is the number of 48-byte windows to process, anywhere from 1 to 10
  //   inclusive.
  //
  // This function performs the update
  //   return_vals[0] += (-N) + \sum_i x_iy_i
  //   return_vals[1] += N_y + \sum_{i: nonmissing from y} x_i
  //   return_vals[2] += N_x + \sum_{i: nonmissing from x} y_i
  //   return_vals[3] += N_y - \sum_{i: nonmissing from y} x_i^2
  //   return_vals[4] += N_x - \sum_{i: nonmissing from x} y_i^2
  // where N is the number of samples processed after applying the missingness
  // masks indicated by the subscripts.
  //
  // Computation of terms [1]-[4] is based on the identity
  //   N_y + \sum_{i: nonmissing from y} x_i = popcount2(vec1 & mask2)
  // where "popcount2" refers to starting with two-bit integers instead of
  // one-bit integers in our summing process (this allows us to skip a few
  // operations).  (Once we can assume the presence of hardware popcount, a
  // slightly different implementation may be better.)
  //
  // The trickier [0] computation currently proceeds as follows:
  //
  // 1. zcheck := (vec1 | vec2) & 0x5555...
  // Detects whether at least one member of the pair has a 0/missing value.
  //
  // 2. popcount2(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
  // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i dot
  // product.
  //
  // MULTIPLEX_LD sets of values are usually handled per function call.  If
  // fewer values are present, the ends of all input vectors should be zeroed
  // out.

  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader1;
  __m128i loader2;
  __m128i sum1;
  __m128i sum2;
  __m128i sum11;
  __m128i sum22;
  __m128i sum12;
  __m128i tmp_sum1;
  __m128i tmp_sum2;
  __m128i tmp_sum12;
  __univec acc;
  __univec acc1;
  __univec acc2;
  __univec acc11;
  __univec acc22;
  acc.vi = _mm_setzero_si128();
  acc1.vi = _mm_setzero_si128();
  acc2.vi = _mm_setzero_si128();
  acc11.vi = _mm_setzero_si128();
  acc22.vi = _mm_setzero_si128();
  do {
    loader1 = *vec1++;
    loader2 = *vec2++;
    sum1 = *mask2++;
    sum2 = *mask1++;
    sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    // sum11 = _mm_and_si128(_mm_and_si128(_mm_xor_si128(sum1, m1), m1), loader1);
    // sum22 = _mm_and_si128(_mm_and_si128(_mm_xor_si128(sum2, m1), m1), loader2);
    sum1 = _mm_and_si128(sum1, loader1);
    sum2 = _mm_and_si128(sum2, loader2);
    sum11 = _mm_and_si128(sum1, m1);
    sum22 = _mm_and_si128(sum2, m1);
    // use andnot to eliminate need for 0xaaaa... to occupy an xmm register
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12), _mm_xor_si128(loader1, loader2));
    sum12 = _mm_or_si128(sum12, loader1);

    // sum1, sum2, and sum12 now store the (biased) two-bit sums of
    // interest; merge to 4 bits to prevent overflow.  this merge can be
    // postponed for sum11 and sum22 because the individual terms are 0/1
    // instead of 0/1/2.
    sum1 = _mm_add_epi64(_mm_and_si128(sum1, m2), _mm_and_si128(_mm_srli_epi64(sum1, 2), m2));
    sum2 = _mm_add_epi64(_mm_and_si128(sum2, m2), _mm_and_si128(_mm_srli_epi64(sum2, 2), m2));
    sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2), _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
    tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
    sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
    sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

    sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
    sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
    tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
    sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
    sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

    sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
    sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
    sum11 = _mm_add_epi64(_mm_and_si128(sum11, m2), _mm_and_si128(_mm_srli_epi64(sum11, 2), m2));
    sum22 = _mm_add_epi64(_mm_and_si128(sum22, m2), _mm_and_si128(_mm_srli_epi64(sum22, 2), m2));
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    acc1.vi = _mm_add_epi64(acc1.vi, _mm_add_epi64(_mm_and_si128(sum1, m4), _mm_and_si128(_mm_srli_epi64(sum1, 4), m4)));
    acc2.vi = _mm_add_epi64(acc2.vi, _mm_add_epi64(_mm_and_si128(sum2, m4), _mm_and_si128(_mm_srli_epi64(sum2, 4), m4)));
    acc11.vi = _mm_add_epi64(acc11.vi, _mm_add_epi64(_mm_and_si128(sum11, m4), _mm_and_si128(_mm_srli_epi64(sum11, 4), m4)));
    acc22.vi = _mm_add_epi64(acc22.vi, _mm_add_epi64(_mm_and_si128(sum22, m4), _mm_and_si128(_mm_srli_epi64(sum22, 4), m4)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sum12, m4), _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
  } while (--iters);
  // moved down because we've almost certainly run out of xmm registers
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
#if MULTIPLEX_LD > 960
  acc1.vi = _mm_add_epi64(_mm_and_si128(acc1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
  acc2.vi = _mm_add_epi64(_mm_and_si128(acc2.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc1.vi = _mm_and_si128(_mm_add_epi64(acc1.vi, _mm_srli_epi64(acc1.vi, 8)), m8);
  acc2.vi = _mm_and_si128(_mm_add_epi64(acc2.vi, _mm_srli_epi64(acc2.vi, 8)), m8);
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc11.vi = _mm_and_si128(_mm_add_epi64(acc11.vi, _mm_srli_epi64(acc11.vi, 8)), m8);
  acc22.vi = _mm_and_si128(_mm_add_epi64(acc22.vi, _mm_srli_epi64(acc22.vi, 8)), m8);

  return_vals[0] -= ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[1] += ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[2] += ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[3] += ((acc11.u8[0] + acc11.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[4] += ((acc22.u8[0] + acc22.u8[1]) * 0x1000100010001LLU) >> 48;
}

void ld_dot_prod(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1, uintptr_t* mask2, int32_t* return_vals, uint32_t batch_ct_m1, uint32_t last_batch_size) {
  while (batch_ct_m1--) {
    ld_dot_prod_batch((__m128i*)vec1, (__m128i*)vec2, (__m128i*)mask1, (__m128i*)mask2, return_vals, MULTIPLEX_LD / 192);
    vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
    vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
    mask1 = &(mask1[MULTIPLEX_LD / BITCT2]);
    mask2 = &(mask2[MULTIPLEX_LD / BITCT2]);
  }
  ld_dot_prod_batch((__m128i*)vec1, (__m128i*)vec2, (__m128i*)mask1, (__m128i*)mask2, return_vals, last_batch_size);
}

static inline int32_t ld_dot_prod_nm_batch(__m128i* vec1, __m128i* vec2, uint32_t iters) {
  // faster ld_dot_prod_batch() for no-missing-calls case.
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader1;
  __m128i loader2;
  __m128i sum12;
  __m128i tmp_sum12;
  __univec acc;
  acc.vi = _mm_setzero_si128();
  do {
    loader1 = *vec1++;
    loader2 = *vec2++;
    sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12), _mm_xor_si128(loader1, loader2));
    sum12 = _mm_or_si128(sum12, loader1);
    sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2), _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sum12, m4), _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
  } while (--iters);
#if MULTIPLEX_LD > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  return (uint32_t)(((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48);
}

int32_t ld_dot_prod_nm(uintptr_t* vec1, uintptr_t* vec2, uint32_t founder_ct, uint32_t batch_ct_m1, uint32_t last_batch_size) {
  // accelerated implementation for no-missing-loci case
  int32_t result = (int32_t)founder_ct;
  while (batch_ct_m1--) {
    result -= ld_dot_prod_nm_batch((__m128i*)vec1, (__m128i*)vec2, MULTIPLEX_LD / 192);
    vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
    vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
  }
  result -= ld_dot_prod_nm_batch((__m128i*)vec1, (__m128i*)vec2, last_batch_size);
  return result;
}
#else
static inline void ld_dot_prod_batch(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1, uintptr_t* mask2, int32_t* return_vals, uint32_t iters) {
  uint32_t final_sum1 = 0;
  uint32_t final_sum2 = 0;
  uint32_t final_sum11 = 0;
  uint32_t final_sum22 = 0;
  uint32_t final_sum12 = 0;
  uintptr_t loader1;
  uintptr_t loader2;
  uintptr_t sum1;
  uintptr_t sum2;
  uintptr_t sum11;
  uintptr_t sum22;
  uintptr_t sum12;
  uintptr_t tmp_sum1;
  uintptr_t tmp_sum2;
  uintptr_t tmp_sum12;
  do {
    // (The important part of the header comment on the 64-bit version is
    // copied below.)
    //
    // Input:
    // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
    // * mask1 and mask2 mask out missing values (i.e. 00 for missing, 11 for
    //   nonmissing).
    // * return_vals provides space for return values.
    // * iters is the number of 12-byte windows to process, anywhere from 1 to
    //   40 inclusive.  (No, this is not the interface you'd use for a
    //   general-purpose library.)  [32- and 64-bit differ here.]
    //
    // This function performs the update
    //   return_vals[0] += (-N) + \sum_i x_iy_i
    //   return_vals[1] += N_y + \sum_{i: nonmissing from y} x_i
    //   return_vals[2] += N_x + \sum_{i: nonmissing from x} y_i
    //   return_vals[3] += N_y - \sum_{i: nonmissing from y} x_i^2
    //   return_vals[4] += N_x - \sum_{i: nonmissing from x} y_i^2
    // where N is the number of samples processed after applying the
    // missingness masks indicated by the subscripts.
    //
    // Computation of terms [1]-[4] is based on the identity
    //   N_y + \sum_{i: nonmissing from y} x_i = popcount2(vec1 & mask2)
    // where "popcount2" refers to starting with two-bit integers instead of
    // one-bit integers in our summing process (this allows us to skip a few
    // operations).  (Once we can assume the presence of hardware popcount, a
    // slightly different implementation may be better.)
    //
    // The trickier [0] computation currently proceeds as follows:
    //
    // 1. zcheck := (vec1 | vec2) & 0x5555...
    // Detects whether at least one member of the pair has a 0/missing value.
    //
    // 2. popcount2(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
    // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i
    // dot product.


    loader1 = *vec1++;
    loader2 = *vec2++;
    sum1 = *mask2++;
    sum2 = *mask1++;
    sum12 = (loader1 | loader2) & FIVEMASK;

    sum1 = sum1 & loader1;
    sum2 = sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
    sum12 = sum12 | loader1;
    sum11 = sum1 & FIVEMASK;
    sum22 = sum2 & FIVEMASK;

    sum1 = (sum1 & 0x33333333) + ((sum1 >> 2) & 0x33333333);
    sum2 = (sum2 & 0x33333333) + ((sum2 >> 2) & 0x33333333);
    sum12 = (sum12 & 0x33333333) + ((sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;
    tmp_sum1 = tmp_sum1 & loader1;
    tmp_sum2 = tmp_sum2 & loader2;

    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;
    sum11 += tmp_sum1 & FIVEMASK;
    sum22 += tmp_sum2 & FIVEMASK;

    sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
    sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;

    tmp_sum1 = tmp_sum1 & loader1;
    tmp_sum2 = tmp_sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;
    sum11 += tmp_sum1 & FIVEMASK;
    sum22 += tmp_sum2 & FIVEMASK;

    sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
    sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
    sum11 = (sum11 & 0x33333333) + ((sum11 >> 2) & 0x33333333);
    sum22 = (sum22 & 0x33333333) + ((sum22 >> 2) & 0x33333333);
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    sum1 = (sum1 & 0x0f0f0f0f) + ((sum1 >> 4) & 0x0f0f0f0f);
    sum2 = (sum2 & 0x0f0f0f0f) + ((sum2 >> 4) & 0x0f0f0f0f);
    sum11 = (sum11 & 0x0f0f0f0f) + ((sum11 >> 4) & 0x0f0f0f0f);
    sum22 = (sum22 & 0x0f0f0f0f) + ((sum22 >> 4) & 0x0f0f0f0f);
    sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);

    // technically could do the multiply-and-shift only once every two rounds
    final_sum1 += (sum1 * 0x01010101) >> 24;
    final_sum2 += (sum2 * 0x01010101) >> 24;
    final_sum11 += (sum11 * 0x01010101) >> 24;
    final_sum22 += (sum22 * 0x01010101) >> 24;
    final_sum12 += (sum12 * 0x01010101) >> 24;
  } while (--iters);
  return_vals[0] -= final_sum12;
  return_vals[1] += final_sum1;
  return_vals[2] += final_sum2;
  return_vals[3] += final_sum11;
  return_vals[4] += final_sum22;
}

void ld_dot_prod(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1, uintptr_t* mask2, int32_t* return_vals, uint32_t batch_ct_m1, uint32_t last_batch_size) {
  while (batch_ct_m1--) {
    ld_dot_prod_batch(vec1, vec2, mask1, mask2, return_vals, MULTIPLEX_LD / 48);
    vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
    vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
    mask1 = &(mask1[MULTIPLEX_LD / BITCT2]);
    mask2 = &(mask2[MULTIPLEX_LD / BITCT2]);
  }
  ld_dot_prod_batch(vec1, vec2, mask1, mask2, return_vals, last_batch_size);
}

static inline int32_t ld_dot_prod_nm_batch(uintptr_t* vec1, uintptr_t* vec2, uint32_t iters) {
  uint32_t final_sum12 = 0;
  uintptr_t loader1;
  uintptr_t loader2;
  uintptr_t sum12;
  uintptr_t tmp_sum12;
  do {
    loader1 = *vec1++;
    loader2 = *vec2++;
    sum12 = (loader1 | loader2) & FIVEMASK;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
    sum12 = sum12 | loader1;
    sum12 = (sum12 & 0x33333333) + ((sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);
    sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);

    final_sum12 += (sum12 * 0x01010101) >> 24;
  } while (--iters);
  return final_sum12;
}

int32_t ld_dot_prod_nm(uintptr_t* vec1, uintptr_t* vec2, uint32_t founder_ct, uint32_t batch_ct_m1, uint32_t last_batch_size) {
  int32_t result = (int32_t)founder_ct;
  while (batch_ct_m1--) {
    result -= ld_dot_prod_nm_batch(vec1, vec2, MULTIPLEX_LD / 48);
    vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
    vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
  }
  result -= ld_dot_prod_nm_batch(vec1, vec2, last_batch_size);
  return result;
}
#endif // __LP64__

uint32_t ld_process_load(uintptr_t* geno_buf, uintptr_t* mask_buf, uintptr_t* missing_buf, uint32_t* missing_ct_ptr, double* sum_ptr, double* variance_recip_ptr, uint32_t founder_ct, uint32_t is_x, uint32_t weighted_x, uint32_t nonmale_founder_ct, uintptr_t* founder_male_include2, uintptr_t* nonmale_geno, uintptr_t* nonmale_masks, uintptr_t nonmale_offset) {
  uintptr_t* geno_ptr = geno_buf;
  uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(founder_ct);
  uintptr_t* geno_end = &(geno_buf[founder_ctl2]);
  uintptr_t* mask_buf_ptr = mask_buf;
  uintptr_t* missing_ptr = missing_buf;
  uintptr_t new_missing = 0;
  int64_t llii;
  uint32_t missing_bit_offset = 0;
  uint32_t ssq = 0;
  uint32_t missing_ct = 0;
  int32_t sum = -founder_ct;
  uintptr_t* nm_mask_ptr;
  uintptr_t cur_geno;
  uintptr_t shifted_masked_geno;
  uintptr_t new_geno;
  uintptr_t new_mask;
  while (1) {
    // Desired encodings:
    // new_geno: nonset homozygote -> 00
    //           het/missing       -> 01
    //           set homozygote    -> 10
    // Given PLINK encoding xx, this is (xx - ((xx >> 1) & FIVEMASK)).
    //
    // new_mask: missing   -> 00
    //           otherwise -> 11
    // ...and this is (((xx >> 1) & FIVEMASK) | ((~xx) & FIVEMASK)) * 3.
    //
    // new_missing: missing   -> 1
    //              otherwise -> 0
    // This can be assembled via repeated CTZLU on ~new_mask.
    cur_geno = *geno_ptr;
    shifted_masked_geno = (cur_geno >> 1) & FIVEMASK;
    new_geno = cur_geno - shifted_masked_geno;
    *geno_ptr++ = new_geno;
    new_mask = (((~cur_geno) & FIVEMASK) | shifted_masked_geno) * 3;
    *mask_buf_ptr++ = new_mask;
    new_mask = (~new_mask) & FIVEMASK;
    while (new_mask) {
      new_missing |= ONELU << (missing_bit_offset + (CTZLU(new_mask) / 2));
      missing_ct++;
      new_mask &= new_mask - 1;
    }
    if (geno_ptr == geno_end) {
      break;
    }
    if (missing_bit_offset) {
      missing_bit_offset = 0;
      *missing_ptr++ = new_missing;
      new_missing = 0;
    } else {
      missing_bit_offset = BITCT2;
    }
  }
  *missing_ptr = new_missing;
  if (is_x && (!weighted_x)) {
    // special case #1: recode male clear homozygotes to 01 on X chromosome,
    // for backwards compatibility
    //
    // this is a bit ugly (e.g. results are actually affected by which allele
    // is A1), so may want to switch the default to mode 3
    geno_ptr = geno_buf;
    do {
      new_geno = *geno_ptr;
      *geno_ptr++ = new_geno + ((~(new_geno | (new_geno >> 1))) & (*founder_male_include2++));
    } while (geno_ptr < geno_end);
  }
  geno_ptr = geno_buf;
  while (1) {
    new_geno = *geno_ptr++;
    sum += popcount2_long(new_geno);
    new_geno = (new_geno ^ FIVEMASK) & FIVEMASK;
    if (geno_ptr == geno_end) {
      break;
    }
    ssq += popcount2_long(new_geno);
  }
  // have to be careful with trailing zeroes here
  ssq += popcount2_long(new_geno << (BITCT - 2 * (1 + ((founder_ct - 1) % BITCT2))));
  if (founder_ct % BITCT2) {
    mask_buf[founder_ct / BITCT2] &= (ONELU << (2 * (founder_ct % BITCT2))) - ONELU;
  }
  if (is_x && weighted_x) {
    // special case #2: double-count nonmales
    geno_ptr = geno_buf;
    sum -= founder_ct;
    nonmale_geno = &(nonmale_geno[nonmale_offset]);
    nonmale_masks = &(nonmale_masks[nonmale_offset]);
    mask_buf_ptr = mask_buf;
    nm_mask_ptr = nonmale_masks;
    while (1) {
      new_mask = ~((*founder_male_include2) * 3);
      new_geno = ((*geno_ptr++) & new_mask) | (*founder_male_include2++);
      *nonmale_geno++ = new_geno;
      *nm_mask_ptr++ = new_mask & (*mask_buf_ptr++);
      sum += popcount2_long(new_geno);
      new_geno = (new_geno ^ FIVEMASK) & FIVEMASK;
      if (geno_ptr == geno_end) {
	break;
      }
      ssq += popcount2_long(new_geno);
    }
    ssq += popcount2_long(new_geno << (BITCT - 2 * (1 + ((founder_ct - 1) % BITCT2))));
    missing_ct += founder_ct - (popcount_longs(nonmale_masks, founder_ctl2) / 2);
    founder_ct *= 2;
  } else if (!missing_ct) {
    // save sum and (n^2)/variance, for faster processing of pairwise
    // no-missing-calls case
    llii = (int64_t)((uint64_t)ssq) * founder_ct - ((int64_t)sum) * sum;
    if (!llii) {
      return 0;
    }
    *missing_ct_ptr = 0;
    *sum_ptr = (double)sum;
    *variance_recip_ptr = 1.0 / ((double)llii);
    return 1;
  }
  *missing_ct_ptr = missing_ct;
  return (((int64_t)((uint64_t)ssq)) * (founder_ct - missing_ct) - ((int64_t)sum) * sum)? 1 : 0;
}

uint32_t ld_prune_next_valid_chrom_start(uintptr_t* marker_exclude, uint32_t cur_uidx, Chrom_info* chrom_info_ptr, uint32_t chrom_code_end, uint32_t unfiltered_marker_ct) {
  uint32_t chrom_idx;
  cur_uidx = next_unset(marker_exclude, cur_uidx, unfiltered_marker_ct);
  while (cur_uidx < unfiltered_marker_ct) {
    chrom_idx = get_variant_chrom(chrom_info_ptr, cur_uidx);
    // --aec 0 support
    if (chrom_idx && (chrom_idx < chrom_code_end)) {
      return cur_uidx;
    }
    cur_uidx = next_unset(marker_exclude, get_chrom_end_vidx(chrom_info_ptr, chrom_idx), unfiltered_marker_ct);
  }
  return cur_uidx;
}

void ld_prune_start_chrom(uint32_t ld_window_kb, uint32_t* cur_chrom_ptr, uint32_t* chrom_end_ptr, uint32_t window_unfiltered_start, uint32_t* live_indices, uint32_t* start_arr, uint32_t* window_unfiltered_end_ptr, uint32_t ld_window_size, uint32_t* cur_window_size_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uint32_t* is_haploid_ptr, uint32_t* is_x_ptr, uint32_t* is_y_ptr) {
  uint32_t cur_chrom = get_variant_chrom(chrom_info_ptr, window_unfiltered_start);
  uint32_t window_unfiltered_end = window_unfiltered_start + 1;
  uint32_t chrom_end = get_chrom_end_vidx(chrom_info_ptr, cur_chrom);
  uint32_t uii = 0;
  uint32_t window_size;
  live_indices[0] = window_unfiltered_start;
  next_unset_ck(marker_exclude, chrom_end, &window_unfiltered_end);
  if (ld_window_kb) {
    window_size = 1;
    uii = window_unfiltered_end;
    while ((uii < chrom_end) && (marker_pos[uii] <= marker_pos[window_unfiltered_start] + (1000 * ld_window_size))) {
      window_size++;
      uii++;
      next_unset_ck(marker_exclude, chrom_end, &uii);
    }
    uii = 0;
  } else {
    window_size = ld_window_size;
  }
  for (uii = 1; uii < window_size; uii++) {
    if (window_unfiltered_end == chrom_end) {
      break;
    }
    start_arr[uii - 1] = window_unfiltered_end;
    live_indices[uii] = window_unfiltered_end;
    window_unfiltered_end++;
    next_unset_ck(marker_exclude, chrom_end, &window_unfiltered_end);
  }
  *cur_window_size_ptr = uii;
  start_arr[uii - 1] = window_unfiltered_end;
  *cur_chrom_ptr = cur_chrom;
  *chrom_end_ptr = chrom_end;
  *window_unfiltered_end_ptr = window_unfiltered_end;
  *is_haploid_ptr = IS_SET(chrom_info_ptr->haploid_mask, cur_chrom);
  *is_x_ptr = (((int32_t)cur_chrom) == chrom_info_ptr->xymt_codes[X_OFFSET]);
  *is_y_ptr = (((int32_t)cur_chrom) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
}

int32_t ld_prune_write(char* outname, char* outname_end, uintptr_t* marker_exclude, uintptr_t* pruned_arr, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, uint32_t chrom_code_end) {
  FILE* outfile = nullptr;
  int32_t retval = 0;
  {
    fputs("Writing...", stdout);
    fflush(stdout);
    strcpy(outname_end, ".prune.in");
    if (fopen_checked(outname, "w", &outfile)) {
      goto ld_prune_write_ret_OPEN_FAIL;
    }
    for (uint32_t cur_chrom = 1; cur_chrom < chrom_code_end; cur_chrom++) {
      if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom)) {
	continue;
      }
      const uint32_t chrom_end = get_chrom_end_vidx(chrom_info_ptr, cur_chrom);
      for (uint32_t marker_uidx = get_chrom_start_vidx(chrom_info_ptr, cur_chrom); marker_uidx < chrom_end; marker_uidx++) {
	// pruned_arr initialized to marker_exclude
	if (!IS_SET(pruned_arr, marker_uidx)) {
	  fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
	  putc_unlocked('\n', outfile);
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto ld_prune_write_ret_WRITE_FAIL;
    }
    strcpy(outname_end, ".prune.out");
    if (fopen_checked(outname, "w", &outfile)) {
      goto ld_prune_write_ret_OPEN_FAIL;
    }
    for (uint32_t cur_chrom = 1; cur_chrom < chrom_code_end; cur_chrom++) {
      if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom)) {
	continue;
      }
      const uint32_t chrom_end = get_chrom_end_vidx(chrom_info_ptr, cur_chrom);
      for (uint32_t marker_uidx = get_chrom_start_vidx(chrom_info_ptr, cur_chrom); marker_uidx < chrom_end; marker_uidx++) {
	if ((!IS_SET(marker_exclude, marker_uidx)) && IS_SET(pruned_arr, marker_uidx)) {
	  fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
	  putc_unlocked('\n', outfile);
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto ld_prune_write_ret_WRITE_FAIL;
    }
    *outname_end = '\0';
    putc_unlocked('\r', stdout);
    LOGPRINTFWW("Marker lists written to %s.prune.in and %s.prune.out .\n", outname, outname);
  }
  while (0) {
  ld_prune_write_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ld_prune_write_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t ld_prune(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  // Results are slightly different from PLINK 1.07 when missing data is
  // present, but that's due to a minor bug in 1.07 (sample per-marker
  // variances don't exclude the missing markers).

  // for future consideration: chromosome-based multithread/parallel?
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct);
  uintptr_t founder_ct = popcount_longs(founder_info, unfiltered_sample_ctv2 / 2);
  uintptr_t founder_ctl = BITCT_TO_WORDCT(founder_ct);
#ifdef __LP64__
  uintptr_t founder_ctv = BITCT_TO_ALIGNED_WORDCT(founder_ct);
#else
  uintptr_t founder_ctv = founder_ctl;
#endif
  uintptr_t founder_ct_mld = (founder_ct + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
  uint32_t founder_ct_mld_m1 = ((uint32_t)founder_ct_mld) - 1;
#ifdef __LP64__
  uint32_t founder_ct_mld_rem = (MULTIPLEX_LD / 192) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 192;
#else
  uint32_t founder_ct_mld_rem = (MULTIPLEX_LD / 48) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 48;
#endif
  uintptr_t founder_ct_192_long = founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2) + founder_ct_mld_rem * (192 / BITCT2);
  uintptr_t final_mask = get_final_mask(founder_ct);
  uint32_t weighted_founder_ct = founder_ct;
  uint32_t founder_trail_ct = founder_ct_192_long - founder_ctl * 2;
  uint32_t pairwise = (ldip->modifier / LD_PRUNE_PAIRWISE) & 1;
  uint32_t ignore_x = (ldip->modifier / LD_IGNORE_X) & 1;
  uint32_t weighted_x = (ldip->modifier / LD_WEIGHTED_X) & 1;
  uint32_t window_is_kb = (ldip->modifier / LD_PRUNE_KB_WINDOW) & 1;
  uint32_t ld_window_size = ldip->prune_window_size;
  uint32_t ld_window_incr = ldip->prune_window_incr;
  double ld_last_param = ldip->prune_last_param;
  uint32_t nonmale_founder_ct = 0;
  uintptr_t window_max = 1;
  uintptr_t* geno = nullptr;
  uintptr_t* founder_include2 = nullptr;
  uintptr_t* founder_male_include2 = nullptr;
  uintptr_t* nonmale_geno = nullptr;
  uintptr_t* nonmale_masks = nullptr;
  double* cov_matrix = nullptr;
  double* new_cov_matrix = nullptr;
  MATRIX_INVERT_BUF1_TYPE* irow = nullptr;
  double* work = nullptr;
  uint32_t* idx_remap = nullptr;
  uint32_t tot_exclude_ct = 0;
  uint32_t at_least_one_prune = 0;
  uint32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  int32_t retval = 0;
  uintptr_t* geno_masks;
  uintptr_t* geno_mmasks;
  uintptr_t* pruned_arr;
  uint32_t* live_indices;
  uint32_t* start_arr;
  uint32_t pct;
  uint32_t pct_thresh;
  uint32_t window_unfiltered_start;
  uint32_t window_unfiltered_end;
  uint32_t cur_window_size;
  uint32_t old_window_size;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  int32_t ii;
  uint32_t cur_chrom;
  uint32_t chrom_start;
  uint32_t chrom_end;
  uint32_t is_haploid;
  uint32_t is_x;
  uint32_t is_y;
  uintptr_t* loadbuf;
  double* sums;
  double* variance_recips; // entries are actually n^2 / variance
  uint32_t* missing_cts;
  uint32_t fixed_missing_ct;
  uintptr_t ulii;
  double dxx;
  double dyy;
  double cov12;
  uint32_t fixed_non_missing_ct;
  uint32_t non_missing_ct;
  int32_t dp_result[5];
  double non_missing_ctd;
  uintptr_t* geno_fixed_vec_ptr;
  uintptr_t* geno_var_vec_ptr;
  uintptr_t* mask_fixed_vec_ptr;
  uintptr_t* mask_var_vec_ptr;
  uintptr_t cur_exclude_ct;
  uint32_t prev_end;
  __CLPK_integer window_rem_li;
  uint32_t old_window_rem;
  uint32_t window_rem;
  uint32_t bsearch_min;
  uint32_t bsearch_max;
  uint32_t bsearch_cur;
  double prune_ld_thresh;

  if (founder_ct < 2) {
    LOGERRPRINTF("Warning: Skipping --indep%s since there are less than two founders.\n(--make-founders may come in handy here.)\n", pairwise? "-pairwise" : "");
    goto ld_prune_ret_1;
  }
  if (is_set(chrom_info_ptr->chrom_mask, 0)) {
    ulii = count_chrom_markers(chrom_info_ptr, marker_exclude, 0);
    if (chrom_info_ptr->zero_extra_chroms) {
      for (uii = chrom_info_ptr->max_code + 1; uii < chrom_code_end; uii++) {
	ulii += count_chrom_markers(chrom_info_ptr, marker_exclude, uii);
      }
      chrom_code_end = chrom_info_ptr->max_code + 1;
    }
    marker_ct -= ulii;
    LOGPRINTF("--indep%s: Ignoring %" PRIuPTR " chromosome 0 variant%s.\n", pairwise? "-pairwise" : "", ulii, (ulii == 1)? "" : "s");
  }
  if (marker_ct < 2) {
    LOGERRPRINTF("Error: Too few valid variants for --indep%s.\n", pairwise? "-pairwise" : "");
    goto ld_prune_ret_INVALID_FORMAT;
  }

  // force founder_male_include2 allocation
  if (alloc_collapsed_haploid_filters(founder_info, sex_male, unfiltered_sample_ct, founder_ct, XMHH_EXISTS | hh_exists, 1, &founder_include2, &founder_male_include2)) {
    goto ld_prune_ret_NOMEM;
  }
  if (weighted_x) {
    nonmale_founder_ct = founder_ct - popcount01_longs(founder_male_include2, founder_ctl);
    if (founder_ct + nonmale_founder_ct > 0x7fffffff) {
      // no, this shouldn't ever happen, but may as well document that there
      // theoretically is a 32-bit integer range issue here
      logerrprint("Error: Too many founders for --indep[-pairwise] + --ld-xchr 3.\n");
      goto ld_prune_ret_1;
    }
  }

  if (window_is_kb) {
    // determine maximum number of markers that may need to be loaded at once
    for (cur_chrom = 1; cur_chrom < chrom_code_end; cur_chrom++) {
      if (is_set(chrom_info_ptr->chrom_mask, cur_chrom)) {
	window_max = chrom_window_max(marker_pos, marker_exclude, chrom_info_ptr, cur_chrom, 0x7fffffff, ld_window_size * 1000, window_max);
      }
    }
  }
  if (pairwise) {
    prune_ld_thresh = ld_last_param * (1 + SMALL_EPSILON);
  } else {
#ifdef __LP64__
    if (window_max > 46340) {
      // todo: check what LAPACK's matrix inversion limit actually is.  Guess
      // sqrt(2^31 - 1) for now.
      logerrprint("Error: --indep does not currently support window sizes > 46340.\n");
      goto ld_prune_ret_INVALID_CMDLINE;
    }
#endif
    // r, not r2, in this case
    prune_ld_thresh = 0.999999;
  }

  window_unfiltered_start = ld_prune_next_valid_chrom_start(marker_exclude, 0, chrom_info_ptr, chrom_code_end, unfiltered_marker_ct);

  if (bigstack_alloc_ul(unfiltered_marker_ctl, &pruned_arr)) {
    goto ld_prune_ret_NOMEM;
  }

  memcpy(pruned_arr, marker_exclude, unfiltered_marker_ctl * sizeof(intptr_t));

  if (!window_is_kb) {
    window_max = ld_window_size;
  }
  if (bigstack_alloc_ui(window_max, &live_indices) ||
      bigstack_alloc_ui(window_max, &start_arr) ||
      bigstack_alloc_ul(unfiltered_sample_ctv2, &loadbuf) ||
      bigstack_alloc_ul(window_max * founder_ct_192_long, &geno) ||
      bigstack_alloc_ul(window_max * founder_ct_192_long, &geno_masks) ||
      bigstack_alloc_ul(window_max * founder_ctv, &geno_mmasks) ||
      bigstack_alloc_ui(window_max, &missing_cts) ||
      bigstack_alloc_d(window_max, &sums) ||
      bigstack_alloc_d(window_max, &variance_recips)) {
    goto ld_prune_ret_NOMEM;
  }
  if (weighted_x) {
    if (bigstack_alloc_ul(window_max * founder_ct_192_long, &nonmale_geno) ||
        bigstack_alloc_ul(window_max * founder_ct_192_long, &nonmale_masks)) {
      goto ld_prune_ret_NOMEM;
    }
  }
  for (ulii = 1; ulii <= window_max; ulii++) {
    fill_ulong_zero(founder_trail_ct + 2, &(geno[ulii * founder_ct_192_long - founder_trail_ct - 2]));
    fill_ulong_zero(founder_trail_ct + 2, &(geno_masks[ulii * founder_ct_192_long - founder_trail_ct - 2]));
    if (weighted_x) {
      fill_ulong_zero(founder_trail_ct + 2, &(nonmale_geno[ulii * founder_ct_192_long - founder_trail_ct - 2]));
      fill_ulong_zero(founder_trail_ct + 2, &(nonmale_masks[ulii * founder_ct_192_long - founder_trail_ct - 2]));
    }
  }
  if (!pairwise) {
    if (bigstack_alloc_d(window_max * window_max, &cov_matrix) ||
        bigstack_alloc_d(window_max * window_max, &new_cov_matrix) ||
        bigstack_alloc_ui(window_max, &idx_remap)) {
      goto ld_prune_ret_NOMEM;
    }

    irow = (MATRIX_INVERT_BUF1_TYPE*)bigstack_alloc(window_max * MATRIX_INVERT_BUF1_CHECKED_ALLOC);
    if (!irow) {
      goto ld_prune_ret_NOMEM;
    }

    if (window_max < 4) {
      ulii = 4;
    } else {
      ulii = window_max;
    }
    if (bigstack_alloc_d(ulii * window_max, &work)) {
      goto ld_prune_ret_NOMEM;
    }
  }
  do {
    prev_end = 0;
    ld_prune_start_chrom(window_is_kb, &cur_chrom, &chrom_end, window_unfiltered_start, live_indices, start_arr, &window_unfiltered_end, ld_window_size, &cur_window_size, unfiltered_marker_ct, pruned_arr, chrom_info_ptr, marker_pos, &is_haploid, &is_x, &is_y);
    if (weighted_x) {
      if (is_x) {
	weighted_founder_ct = 2 * founder_ct;
      } else {
	weighted_founder_ct = founder_ct;
      }
    }
    old_window_size = 0;
    cur_exclude_ct = 0;
    if (cur_window_size > 1) {
      for (ulii = 0; ulii < (uintptr_t)cur_window_size; ulii++) {
	uii = live_indices[ulii];
	if (fseeko(bedfile, bed_offset + (uii * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, uii), bedfile, loadbuf, &(geno[ulii * founder_ct_192_long]))) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(geno[ulii * founder_ct_192_long])));
	}
        if (!ld_process_load(&(geno[ulii * founder_ct_192_long]), &(geno_masks[ulii * founder_ct_192_long]), &(geno_mmasks[ulii * founder_ctv]), &(missing_cts[ulii]), &(sums[ulii]), &(variance_recips[ulii]), founder_ct, is_x && (!ignore_x), weighted_x, nonmale_founder_ct, founder_male_include2, nonmale_geno, nonmale_masks, ulii * founder_ct_192_long)) {
	  SET_BIT(uii, pruned_arr);
          cur_exclude_ct++;
	}
      }
    }
    pct = 1;
    chrom_start = get_chrom_start_vidx(chrom_info_ptr, cur_chrom);
    pct_thresh = window_unfiltered_start + ((uint64_t)pct * (chrom_end - chrom_start)) / 100;
    while ((window_unfiltered_start < chrom_end) || (cur_window_size > 1)) {
      if (cur_window_size > 1) {
	do {
	  at_least_one_prune = 0;
	  for (uii = 0; uii < cur_window_size - 1; uii++) {
	    if (IS_SET(pruned_arr, live_indices[uii])) {
	      continue;
	    }
            fixed_missing_ct = missing_cts[uii];
	    fixed_non_missing_ct = weighted_founder_ct - fixed_missing_ct;
	    geno_fixed_vec_ptr = &(geno[uii * founder_ct_192_long]);
	    mask_fixed_vec_ptr = &(geno_masks[uii * founder_ct_192_long]);
	    ujj = uii + 1;
	    while (live_indices[ujj] < start_arr[uii]) {
	      if (++ujj == cur_window_size) {
		break;
	      }
	    }
	    for (; ujj < cur_window_size; ujj++) {
	      if (IS_SET(pruned_arr, live_indices[ujj])) {
		continue;
	      }
	      geno_var_vec_ptr = &(geno[ujj * founder_ct_192_long]);
	      if ((!fixed_missing_ct) && (!missing_cts[ujj]) && ((!is_x) || (!weighted_x))) {
		cov12 = (double)(ld_dot_prod_nm(geno_fixed_vec_ptr, geno_var_vec_ptr, weighted_founder_ct, founder_ct_mld_m1, founder_ct_mld_rem) * ((int64_t)founder_ct)) - sums[uii] * sums[ujj];
		dxx = variance_recips[uii] * variance_recips[ujj];
	      } else {
		mask_var_vec_ptr = &(geno_masks[ujj * founder_ct_192_long]);
		dp_result[0] = weighted_founder_ct;
		// reversed from what I initially thought because I'm passing
		// the ujj-associated buffers before the uii-associated ones.
		dp_result[1] = -((int32_t)fixed_non_missing_ct);
		dp_result[2] = missing_cts[ujj] - weighted_founder_ct;
		dp_result[3] = dp_result[1];
		dp_result[4] = dp_result[2];
		ld_dot_prod(geno_var_vec_ptr, geno_fixed_vec_ptr, mask_var_vec_ptr, mask_fixed_vec_ptr, dp_result, founder_ct_mld_m1, founder_ct_mld_rem);
		if (is_x && weighted_x) {
		  non_missing_ct = (popcount_longs_intersect(&(nonmale_masks[uii * founder_ct_192_long]), &(nonmale_masks[ujj * founder_ct_192_long]), 2 * founder_ctl) + popcount_longs_intersect(mask_fixed_vec_ptr, mask_var_vec_ptr, 2 * founder_ctl)) / 2;
		  ld_dot_prod(&(nonmale_geno[ujj * founder_ct_192_long]), &(nonmale_geno[uii * founder_ct_192_long]), &(nonmale_masks[ujj * founder_ct_192_long]), &(nonmale_masks[uii * founder_ct_192_long]), dp_result, founder_ct_mld_m1, founder_ct_mld_rem);
		} else {
		  non_missing_ct = fixed_non_missing_ct - missing_cts[ujj];
		  if (fixed_missing_ct && missing_cts[ujj]) {
		    non_missing_ct += popcount_longs_intersect(&(geno_mmasks[uii * founder_ctv]), &(geno_mmasks[ujj * founder_ctv]), founder_ctl);
		  }
		}
		non_missing_ctd = (double)((int32_t)non_missing_ct);
		dxx = dp_result[1];
		dyy = dp_result[2];
		cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
		dxx = 1.0 / ((dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy));
	      }
	      if (!pairwise) {
		dxx = cov12 * sqrt(dxx);
		if (dxx != dxx) {
		  // force prune if 0/0 for now
		  dxx = 1.0;
		}
		cov_matrix[uii * window_max + ujj] = dxx;
	      } else {
		dxx = cov12 * cov12 * dxx;
	      }
	      if (dxx > prune_ld_thresh) {
		at_least_one_prune = 1;
		cur_exclude_ct++;
		// remove marker with lower MAF
		// could cache MAFs of all current-window variants, but
		// get_maf() is too cheap for this to make a noticeable
		// difference
		if (get_maf(set_allele_freqs[live_indices[uii]]) < (1 - SMALL_EPSILON) * get_maf(set_allele_freqs[live_indices[ujj]])) {
		  /*
		  if (debug_print) {
		    printf("removed %u, kept %u, MAFs %g/%g, r2 %g\n", live_indices[uii], live_indices[ujj], get_maf(set_allele_freqs[live_indices[uii]]), get_maf(set_allele_freqs[live_indices[ujj]]), dxx);
		  }
		  */
		  SET_BIT(live_indices[uii], pruned_arr);
		} else {
		  /*
		  if (debug_print) {
		    printf("removed %u, kept %u, MAFs %g/%g, r2 %g\n", live_indices[ujj], live_indices[uii], get_maf(set_allele_freqs[live_indices[ujj]]), get_maf(set_allele_freqs[live_indices[uii]]), dxx);
		  }
		  */
		  SET_BIT(live_indices[ujj], pruned_arr);
		  ujj++;
		  while (ujj < cur_window_size) {
		    if (!IS_SET(pruned_arr, live_indices[ujj])) {
		      break;
		    }
		    ujj++;
		  }
		  if (ujj < cur_window_size) {
		    start_arr[uii] = live_indices[ujj];
		  }
		}
		break;
	      }
	    }
	    if (ujj == cur_window_size) {
	      start_arr[uii] = window_unfiltered_end;
	    }
	  }
	} while (at_least_one_prune);
	if (!pairwise) {
	  window_rem = 0;
	  for (uii = 0; uii < old_window_size; uii++) {
	    if (IS_SET(pruned_arr, live_indices[uii])) {
	      continue;
	    }
            idx_remap[window_rem++] = uii;
	  }
	  old_window_rem = window_rem;
	  for (; uii < cur_window_size; uii++) {
	    if (IS_SET(pruned_arr, live_indices[uii])) {
	      continue;
	    }
            idx_remap[window_rem++] = uii;
	  }
	  while (window_rem > 1) {
	    new_cov_matrix[0] = 1.0;
	    for (uii = 1; uii < window_rem; uii++) {
	      ukk = idx_remap[uii];
	      for (ujj = 0; ujj < uii; ujj++) {
		dxx = cov_matrix[idx_remap[ujj] * window_max + ukk];
		new_cov_matrix[ujj * window_rem + uii] = dxx;
		new_cov_matrix[uii * window_rem + ujj] = dxx;
	      }
	      new_cov_matrix[uii * (window_rem + 1)] = 1.0;
	    }
	    window_rem_li = window_rem;
	    ii = invert_matrix_checked(window_rem_li, new_cov_matrix, irow, work);
	    while (ii) {
#ifdef NOLAPACK
	      if (ii == -1) {
		goto ld_prune_ret_NOMEM;
	      }
#endif
	      // 1. binary search for minimum number of bottom right rows/
	      //    columns that must be trimmed to get a nonsingular matrix
	      bsearch_max = window_rem - 1;
	      if (old_window_rem > bsearch_max) {
		// Normally we can assume that only loci not in the previous
		// window need to be considered here.  But, thanks to numeric
		// instability, we might still need to properly handle an
		// apparently-singular old submatrix?
		old_window_size = 0;
		old_window_rem = 0;
	      }
	      bsearch_min = old_window_rem;
	      while (bsearch_min < bsearch_max) {
	        bsearch_cur = (bsearch_min + bsearch_max) / 2;
                new_cov_matrix[0] = 1.0;
		for (uii = 1; uii < bsearch_cur; uii++) {
		  ukk = idx_remap[uii];
		  for (ujj = 0; ujj < uii; ujj++) {
		    dxx = cov_matrix[idx_remap[ujj] * window_max + ukk];
		    new_cov_matrix[ujj * bsearch_cur + uii] = dxx;
		    new_cov_matrix[uii * bsearch_cur + ujj] = dxx;
		  }
		  new_cov_matrix[uii * (bsearch_cur + 1)] = 1.0;
		}
		if (bsearch_cur) {
		  window_rem_li = bsearch_cur;
		  ii = invert_matrix_checked(window_rem_li, new_cov_matrix, irow, work);
		  if (!ii) {
		    bsearch_min = bsearch_cur + 1;
		  } else {
		    bsearch_max = bsearch_cur;
		  }
		} else {
		  bsearch_min = 1;
		}
	      }

	      // 2. the last trimmed row/column must be part of some linear
	      //    combination.  prune *just* that, and retry.
	      ujj = bsearch_min;
	      // bug reported by Kaustubh was a violation of this:
	      // assert(!IS_SET(pruned_arr, live_indices[idx_remap[ujj]]));
              SET_BIT(live_indices[idx_remap[ujj]], pruned_arr);
	      cur_exclude_ct++;
	      window_rem--;
	      for (uii = ujj; uii < window_rem; uii++) {
		idx_remap[uii] = idx_remap[uii + 1];
	      }
	      new_cov_matrix[0] = 1.0;
	      for (uii = 1; uii < window_rem; uii++) {
		ukk = idx_remap[uii];
		for (ujj = 0; ujj < uii; ujj++) {
		  dxx = cov_matrix[idx_remap[ujj] * window_max + ukk];
		  new_cov_matrix[ujj * window_rem + uii] = dxx;
		  new_cov_matrix[uii * window_rem + ujj] = dxx;
		}
		new_cov_matrix[uii * (window_rem + 1)] = 1.0;
	      }
              window_rem_li = window_rem;
	      ii = invert_matrix_checked(window_rem_li, new_cov_matrix, irow, work);
	    }
	    dxx = new_cov_matrix[0];
	    ujj = 0;
	    for (uii = 1; uii < window_rem; uii++) {
              if (new_cov_matrix[uii * (window_rem + 1)] > dxx) {
		dxx = new_cov_matrix[uii * (window_rem + 1)];
		ujj = uii;
	      }
	    }
	    if (dxx > ld_last_param) {
	      SET_BIT(live_indices[idx_remap[ujj]], pruned_arr);
	      cur_exclude_ct++;
	      window_rem--;
	      if (idx_remap[ujj] < (uint32_t)old_window_size) {
                old_window_rem--;
	      }
	      for (uii = ujj; uii < window_rem; uii++) {
                idx_remap[uii] = idx_remap[uii + 1];
	      }
	    } else {
	      // break out
	      window_rem = 1;
	    }
	  }
	}
      }
      for (uii = 0; uii < ld_window_incr; uii++) {
	if (window_unfiltered_start == chrom_end) {
	  break;
	}
	window_unfiltered_start++;
	next_unset_ck(marker_exclude, chrom_end, &window_unfiltered_start);
      }
      if (window_unfiltered_start == chrom_end) {
	break;
      }
      if (window_unfiltered_start >= pct_thresh) {
	pct = ((window_unfiltered_start - chrom_start) * 100LLU) / (chrom_end - chrom_start);
	printf("\r%u%%", pct++);
	fflush(stdout);
	pct_thresh = chrom_start + (((uint64_t)pct * (chrom_end - chrom_start)) / 100);
      }
      ujj = 0;

      if (window_unfiltered_end < window_unfiltered_start) {
	window_unfiltered_end = window_unfiltered_start;
      }

      // copy back previously loaded/computed results
      while (live_indices[ujj] < window_unfiltered_start) {
	ujj++;
	if (ujj == cur_window_size) {
	  break;
	}
      }
      for (uii = 0; ujj < cur_window_size; ujj++) {
	if (IS_SET(pruned_arr, live_indices[ujj])) {
	  continue;
	}
	memcpy(&(geno[uii * founder_ct_192_long]), &(geno[ujj * founder_ct_192_long]), founder_ct_192_long * sizeof(intptr_t));
	memcpy(&(geno_masks[uii * founder_ct_192_long]), &(geno_masks[ujj * founder_ct_192_long]), founder_ct_192_long * sizeof(intptr_t));
	if (is_x && weighted_x) {
	  memcpy(&(nonmale_geno[uii * founder_ct_192_long]), &(nonmale_geno[ujj * founder_ct_192_long]), founder_ct_192_long * sizeof(intptr_t));
	  memcpy(&(nonmale_masks[uii * founder_ct_192_long]), &(nonmale_masks[ujj * founder_ct_192_long]), founder_ct_192_long * sizeof(intptr_t));
	}
	memcpy(&(geno_mmasks[uii * founder_ctv]), &(geno_mmasks[ujj * founder_ctv]), founder_ctl * sizeof(intptr_t));
	live_indices[uii] = live_indices[ujj];
	start_arr[uii] = start_arr[ujj];
	missing_cts[uii] = missing_cts[ujj];
	sums[uii] = sums[ujj];
        variance_recips[uii] = variance_recips[ujj];
	if (!pairwise) {
	  for (ukk = 0; ukk < uii; ukk++) {
	    cov_matrix[ukk * window_max + uii] = cov_matrix[idx_remap[ukk] * window_max + ujj];
	  }
	  idx_remap[uii] = ujj;
	}
	uii++;
      }

      prev_end = uii;
      cur_window_size = uii;
      if (window_is_kb) {
	ujj = 0;
	ukk = window_unfiltered_end;
	while ((ukk < chrom_end) && (marker_pos[ukk] <= marker_pos[window_unfiltered_start] + (1000 * ld_window_size))) {
	  ujj++;
	  ukk++;
	  next_unset_ck(marker_exclude, chrom_end, &ukk);
	}
      } else {
	ujj = ld_window_incr;
      }
      old_window_size = cur_window_size;
      for (uii = 0; uii < ujj; uii++) {
	if (window_unfiltered_end == chrom_end) {
	  break;
	}
	live_indices[cur_window_size] = window_unfiltered_end;
	if (cur_window_size > prev_end) {
	  start_arr[cur_window_size - 1] = window_unfiltered_end;
	}
	if (fseeko(bedfile, bed_offset + (window_unfiltered_end * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, window_unfiltered_end), bedfile, loadbuf, &(geno[cur_window_size * founder_ct_192_long]))) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(geno[cur_window_size * founder_ct_192_long])));
	}
	if (!ld_process_load(&(geno[cur_window_size * founder_ct_192_long]), &(geno_masks[cur_window_size * founder_ct_192_long]), &(geno_mmasks[cur_window_size * founder_ctv]), &(missing_cts[cur_window_size]), &(sums[cur_window_size]), &(variance_recips[cur_window_size]), founder_ct, is_x && (!ignore_x), weighted_x, nonmale_founder_ct, founder_male_include2, nonmale_geno, nonmale_masks, cur_window_size * founder_ct_192_long)) {
	  SET_BIT(window_unfiltered_end, pruned_arr);
	  cur_exclude_ct++;
	}
	cur_window_size++;
	window_unfiltered_end++;
	next_unset_ck(marker_exclude, chrom_end, &window_unfiltered_end);
      }
      if (cur_window_size > prev_end) {
	start_arr[cur_window_size - 1] = window_unfiltered_end;
      }
    }
    putc_unlocked('\r', stdout);
    LOGPRINTF("Pruned %" PRIuPTR " variant%s from chromosome %u, leaving %" PRIuPTR ".\n", cur_exclude_ct, (cur_exclude_ct == 1)? "" : "s", cur_chrom, chrom_end - chrom_start - popcount_bit_idx(marker_exclude, chrom_start, chrom_end) - cur_exclude_ct);
    tot_exclude_ct += cur_exclude_ct;

    // advance chromosomes as necessary
    window_unfiltered_start = ld_prune_next_valid_chrom_start(pruned_arr, window_unfiltered_start, chrom_info_ptr, chrom_code_end, unfiltered_marker_ct);
  } while (window_unfiltered_start < unfiltered_marker_ct);

  LOGPRINTF("Pruning complete.  %u of %" PRIuPTR " variants removed.\n", tot_exclude_ct, marker_ct);
  retval = ld_prune_write(outname, outname_end, marker_exclude, pruned_arr, marker_ids, max_marker_id_len, chrom_info_ptr, chrom_code_end);
  if (retval) {
    goto ld_prune_ret_1;
  }

  while (0) {
  ld_prune_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_prune_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_prune_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
#ifdef __LP64__
  ld_prune_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
#endif
  }
 ld_prune_ret_1:
  bigstack_reset(bigstack_mark);
  return retval;
}

void ld_process_load2(uintptr_t* geno_buf, uintptr_t* mask_buf, uint32_t* missing_ct_ptr, uint32_t founder_ct, uint32_t is_x, uintptr_t* founder_male_include2) {
  // ld_process_load(), except no missing_buf[] to conserve memory (and no
  // --ld-xchr 3 support yet), and no zero-variance check (we just want to
  // dump nans in that case)
  uintptr_t* geno_ptr = geno_buf;
  uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(founder_ct);
  uintptr_t* geno_end = &(geno_buf[founder_ctl2]);
  uintptr_t* mask_buf_ptr = mask_buf;
  uintptr_t cur_geno;
  uintptr_t shifted_masked_geno;
  uintptr_t new_geno;
  uintptr_t new_mask;
  do {
    cur_geno = *geno_ptr;
    shifted_masked_geno = (cur_geno >> 1) & FIVEMASK;
    new_geno = cur_geno - shifted_masked_geno;
    *geno_ptr++ = new_geno;
    new_mask = (((~cur_geno) & FIVEMASK) | shifted_masked_geno) * 3;
    *mask_buf_ptr++ = new_mask;
  } while (geno_ptr < geno_end);
  if (is_x) {
    geno_ptr = geno_buf;
    do {
      new_geno = *geno_ptr;
      *geno_ptr++ = new_geno + ((~(new_geno | (new_geno >> 1))) & (*founder_male_include2++));
    } while (geno_ptr < geno_end);
  }
  if (founder_ct % BITCT2) {
    mask_buf[founder_ct / BITCT2] &= (ONELU << (2 * (founder_ct % BITCT2))) - ONELU;
  }
  *missing_ct_ptr = founder_ct - (popcount_longs(mask_buf, founder_ctl2) / 2);
}

uint32_t ld_missing_ct_intersect(uintptr_t* lptr1, uintptr_t* lptr2, uintptr_t word12_ct, uintptr_t word12_rem, uintptr_t lshift_last) {
  // variant of popcount_longs_intersect()
  uintptr_t tot = 0;
  uintptr_t* lptr1_end2;
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i* vptr1 = (__m128i*)lptr1;
  __m128i* vptr2 = (__m128i*)lptr2;
  __m128i* vend1;
  __m128i loader1;
  __m128i loader2;
  __univec acc;

  while (word12_ct >= 10) {
    word12_ct -= 10;
    vend1 = &(vptr1[60]);
  ld_missing_ct_intersect_main_loop:
    acc.vi = _mm_setzero_si128();
    do {
      loader1 = _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
      loader2 = _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
      loader1 = _mm_add_epi64(loader1, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
      loader2 = _mm_add_epi64(loader2, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
      loader1 = _mm_add_epi64(loader1, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
      loader2 = _mm_add_epi64(loader2, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
      loader1 = _mm_add_epi64(_mm_and_si128(loader1, m2), _mm_and_si128(_mm_srli_epi64(loader1, 2), m2));
      loader1 = _mm_add_epi64(loader1, _mm_add_epi64(_mm_and_si128(loader2, m2), _mm_and_si128(_mm_srli_epi64(loader2, 2), m2)));
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(loader1, m4), _mm_and_si128(_mm_srli_epi64(loader1, 4), m4)));
    } while (vptr1 < vend1);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (word12_ct) {
    vend1 = &(vptr1[word12_ct * 6]);
    word12_ct = 0;
    goto ld_missing_ct_intersect_main_loop;
  }
  lptr1 = (uintptr_t*)vptr1;
  lptr2 = (uintptr_t*)vptr2;
#else
  uintptr_t* lptr1_end = &(lptr1[word12_ct * 12]);
  uintptr_t tmp_stor;
  uintptr_t loader1;
  uintptr_t loader2;
  while (lptr1 < lptr1_end) {
    loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
    loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
    tmp_stor = (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);

    loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
    loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
    tmp_stor += (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);
    tot += (tmp_stor * 0x01010101) >> 24;
  }
#endif
  lptr1_end2 = &(lptr1[word12_rem]);
  while (lptr1 < lptr1_end2) {
    tot += popcount2_long((~((*lptr1++) | (*lptr2++))) & FIVEMASK);
  }
  if (lshift_last) {
    tot += popcount2_long(((~((*lptr1) | (*lptr2))) & FIVEMASK) << lshift_last);
  }
  return tot;
}

int32_t flipscan(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  FILE* outfile_verbose = nullptr;
  uintptr_t* sample_include2 = nullptr;
  uintptr_t* sample_male_include2 = nullptr;
  double min_corr = ldip->flipscan_thresh * (1 - SMALL_EPSILON);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t marker_idx = 0;
  uintptr_t max_window_size = 1;
  uintptr_t pct = 1;
  uintptr_t pct_thresh = marker_ct / 100;
  uint32_t verbose = (ldip->modifier / LD_FLIPSCAN_VERBOSE) & 1;
  uint32_t ignore_x = (ldip->modifier / LD_IGNORE_X) & 1;
  uint32_t max_window_locus_ct = ldip->flipscan_window_size - 1;
  uint32_t window_bp = ldip->flipscan_window_bp;
  uint32_t problem_ct = 0;
  int32_t retval = 0;
  uintptr_t* founder_phenos[2];
  uintptr_t* pheno_male_include2[2];
  uintptr_t* window_geno[2];
  uintptr_t* window_mask[2];
  uintptr_t pheno_ct[2];
  uintptr_t pheno_ct_192_long[2];
  uint32_t pheno_ctl[2];
  uint32_t pheno_ct_mld_m1[2];
  uint32_t pheno_ct_mld_rem[2];
  int32_t dp_result[5];
  double* r_matrix;
  double* r_matrix_ptr;
  double* r_row_ptr;
  uintptr_t* loadbuf_raw;
  uintptr_t* window_geno_ptr;
  uintptr_t* window_mask_ptr;
  uintptr_t* geno_fixed_vec_ptr;
  uintptr_t* mask_fixed_vec_ptr;
  uintptr_t* geno_var_vec_ptr;
  uintptr_t* mask_var_vec_ptr;
  uint32_t* window_uidxs;
  uint32_t* window_cidx_starts;
  uint32_t* neg_uidx_buf;
  uint32_t* missing_cts;
  uint32_t* missing_cts_ptr;
  char* textbuf;
  char* wptr;
  char* wptr_start;
  char* wptr_start2;
  double pos_r_tot;
  double neg_r_tot;
  double ctrl_pheno;
  double case_pheno;
  double non_missing_ctd;
  double cov12;
  double dxx;
  double dyy;
  uintptr_t marker_uidx;
  uintptr_t cur_pheno_ct;
  uintptr_t window_cidx;
  uintptr_t window_cidx2;
  uintptr_t window_cidx3;
  uintptr_t marker_uidx2;
  uintptr_t marker_uidx3;
  uintptr_t cur_192_long;
  uintptr_t cur_ctwd12;
  uintptr_t cur_ctwd12_rem;
  uintptr_t lshift_last;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t chrom_marker_ct;
  uint32_t chrom_marker_idx;
  uint32_t is_haploid;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_case;
  uint32_t marker_pos_thresh;
  uint32_t pos_r_ct;
  uint32_t neg_r_ct;
  uint32_t fixed_missing_ct;
  uint32_t fixed_non_missing_ct;
  uint32_t non_missing_ct;
  uint32_t cur_mld_m1;
  uint32_t cur_mld_rem;
  uint32_t uii;
  ulii = 2 * (max_marker_allele_len + plink_maxsnp) + 256;
  if (ulii <= MAXLINELEN) {
    textbuf = g_textbuf;
  } else {
    if (bigstack_alloc_c(ulii, &textbuf)) {
      goto flipscan_ret_NOMEM;
    }
  }
  if (bigstack_alloc_ul(unfiltered_sample_ctl, &(founder_phenos[0])) ||
      bigstack_alloc_ul(unfiltered_sample_ctl, &(founder_phenos[1])) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw)) {
    goto flipscan_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  memcpy(founder_phenos[0], founder_info, unfiltered_sample_ctl * sizeof(intptr_t));
  bitvec_and(pheno_nm, unfiltered_sample_ctl, founder_phenos[0]);
  if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_exists, 0, founder_phenos[0], sex_male, &sample_include2, &sample_male_include2)) {
    goto flipscan_ret_NOMEM;
  }
  memcpy(founder_phenos[1], founder_phenos[0], unfiltered_sample_ctl * sizeof(intptr_t));
  bitvec_and(pheno_c, unfiltered_sample_ctl, founder_phenos[1]);
  bitvec_andnot(pheno_c, unfiltered_sample_ctl, founder_phenos[0]);
  pheno_ct[0] = popcount_longs(founder_phenos[0], unfiltered_sample_ctl);
  pheno_ct[1] = popcount_longs(founder_phenos[1], unfiltered_sample_ctl);
  if ((!pheno_ct[0]) || (!pheno_ct[1])) {
    if (popcount_longs(founder_info, unfiltered_sample_ctl)) {
      logerrprint("Error: --flip-scan requires at least one case and one control, and only\nconsiders founders.\n");
    } else {
      logerrprint("Error: --flip-scan requires founders.  (--make-founders may come in handy\nhere.)\n");
    }
    goto flipscan_ret_INVALID_CMDLINE;
  }
  for (is_case = 0; is_case < 2; is_case++) {
    pheno_ctl[is_case] = BITCT_TO_WORDCT(pheno_ct[is_case]);

    // ulii == total number of blocks, all but last is size MULTIPLEX_LD
    ulii = (pheno_ct[is_case] + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
    pheno_ct_mld_m1[is_case] = ulii - 1;

    // number of size-{48,192} sub-blocks in trailing block
#ifdef __LP64__
    pheno_ct_mld_rem[is_case] = (MULTIPLEX_LD / 192) - (ulii * MULTIPLEX_LD - pheno_ct[is_case]) / 192;
#else
    pheno_ct_mld_rem[is_case] = (MULTIPLEX_LD / 48) - (ulii * MULTIPLEX_LD - pheno_ct[is_case]) / 48;
#endif

    // number of genotype words per variant, rounded up to the next 192-sample
    // boundary
    pheno_ct_192_long[is_case] = pheno_ct_mld_m1[is_case] * (MULTIPLEX_LD / BITCT2) + pheno_ct_mld_rem[is_case] * (192 / BITCT2);
  }
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    max_window_size = chrom_window_max(marker_pos, marker_exclude, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], max_window_locus_ct * 2 + 1, window_bp * 2, max_window_size);
  }
  if (bigstack_alloc_ui(max_window_size, &window_uidxs) ||
      bigstack_alloc_ui(max_window_size, &window_cidx_starts) ||
      bigstack_alloc_ui(max_window_size, &neg_uidx_buf) ||
      bigstack_alloc_ul(pheno_ctl[0] * 2, &(pheno_male_include2[0])) ||
      bigstack_alloc_ul(pheno_ctl[1] * 2, &(pheno_male_include2[1])) ||
      bigstack_alloc_ui(max_window_size * 2, &missing_cts) ||
      bigstack_alloc_ul(max_window_size * pheno_ct_192_long[0], &(window_geno[0])) ||
      bigstack_alloc_ul(max_window_size * pheno_ct_192_long[0], &(window_mask[0])) ||
      bigstack_alloc_ul(max_window_size * pheno_ct_192_long[1], &(window_geno[1])) ||
      bigstack_alloc_ul(max_window_size * pheno_ct_192_long[1], &(window_mask[1])) ||
      // not advantageous to choose a very large block size here, so O(n^2)
      // memory is fine (though it can be avoided by calculating each
      // correlation twice).
      bigstack_alloc_d(max_window_size * max_window_size * 2, &r_matrix)) {
    goto flipscan_ret_NOMEM;
  }
  ulii = (max_window_size + 1) * 2;
  for (uljj = 0; uljj < max_window_size; uljj++) {
    neg_uidx_buf[uljj * ulii] = 0.0;
    neg_uidx_buf[uljj * ulii + 1] = 0.0;
    // bugfix: initialize r_matrix diagonal
    r_matrix[uljj * ulii] = 0.0;
    r_matrix[uljj * ulii + 1] = 0.0;
  }
  for (is_case = 0; is_case < 2; is_case++) {
    quaterarr_collapse_init(sex_male, unfiltered_sample_ct, founder_phenos[is_case], pheno_ct[is_case], pheno_male_include2[is_case]);
    window_geno_ptr = window_geno[is_case];
    window_mask_ptr = window_mask[is_case];
    cur_192_long = pheno_ct_192_long[is_case];
    ulii = 2 + pheno_ct_192_long[is_case] - pheno_ctl[is_case] * 2;
    for (uljj = 1; uljj <= max_window_size; uljj++) {
      fill_ulong_zero(ulii, &(window_geno_ptr[uljj * cur_192_long - ulii]));
      fill_ulong_zero(ulii, &(window_mask_ptr[uljj * cur_192_long - ulii]));
    }
  }

  memcpy(outname_end, ".flipscan", 10);
  if (fopen_checked(outname, "w", &outfile)) {
    goto flipscan_ret_OPEN_FAIL;
  }
  wptr = memcpya(textbuf, "   CHR ", 7);
  wptr = fw_strcpyn(plink_maxsnp, 3, "SNP", wptr);
  wptr = strcpya(wptr, "           BP   A1   A2        F    POS    R_POS    NEG    R_NEG NEGSNPS\n");
  if (fwrite_checked(textbuf, wptr - textbuf, outfile)) {
    goto flipscan_ret_WRITE_FAIL;
  }
  if (verbose) {
    memcpy(&(outname_end[9]), ".verbose", 9);
    if (fopen_checked(outname, "w", &outfile_verbose)) {
      goto flipscan_ret_OPEN_FAIL;
    }
    outname_end[9] = '\0';
    // er, this is a misalignment disaster
    wptr = memcpya(textbuf, "CHR_INDX ", 9);
    wptr = fw_strcpyn(plink_maxsnp, 8, "SNP_INDX", wptr);
    wptr = memcpya(wptr, "      BP_INDX A1_INDX ", 22);
    wptr = fw_strcpyn(plink_maxsnp, 8, "SNP_PAIR", wptr);
    wptr = strcpya(wptr, "      BP_PAIR A1_PAIR      R_A      R_U\n");
    if (fwrite_checked(textbuf, wptr - textbuf, outfile_verbose)) {
      goto flipscan_ret_WRITE_FAIL;
    }
  }
  printf("--flip-scan%s: 0%%", verbose? " verbose" : "");
  fflush(stdout);
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
    chrom_marker_ct = chrom_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, chrom_end);
    if (chrom_marker_ct < 2) {
      marker_idx += chrom_marker_ct;
      continue;
    }
    wptr_start = width_force(6, textbuf, chrom_name_write(chrom_info_ptr, chrom_idx, textbuf));
    *wptr_start++ = ' ';
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    is_x = (chrom_idx == ((uint32_t)chrom_info_ptr->xymt_codes[X_OFFSET]));
    is_y = (chrom_idx == ((uint32_t)chrom_info_ptr->xymt_codes[Y_OFFSET]));
    if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto flipscan_ret_READ_FAIL;
    }
    chrom_marker_idx = 0;
    window_cidx = max_window_size - 1;
    window_cidx2 = 0;
    do {
      if (++window_cidx == max_window_size) {
	window_cidx = 0;
      }
      window_uidxs[window_cidx] = marker_uidx;

      // circular index of beginning of window starting at current marker
      window_cidx_starts[window_cidx] = window_cidx2;
      if (load_raw(unfiltered_sample_ct4, bedfile, loadbuf_raw)) {
	goto flipscan_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf_raw);
      }
      if (is_haploid && hh_exists) {
        haploid_fix(hh_exists, sample_include2, sample_male_include2, unfiltered_sample_ct, is_x, is_y, (unsigned char*)loadbuf_raw);
      }
      for (is_case = 0; is_case < 2; is_case++) {
	// similar to ld_block_thread() below
	cur_pheno_ct = pheno_ct[is_case];
	uii = cur_pheno_ct / BITCT2;
        cur_ctwd12 = uii / 12;
	cur_ctwd12_rem = uii - (12 * cur_ctwd12);
	lshift_last = 2 * ((0x7fffffc0 - cur_pheno_ct) % BITCT2);
	cur_mld_m1 = pheno_ct_mld_m1[is_case];
        cur_mld_rem = pheno_ct_mld_rem[is_case];
	cur_192_long = pheno_ct_192_long[is_case];
	window_geno_ptr = window_geno[is_case];
	window_mask_ptr = window_mask[is_case];
	missing_cts_ptr = &(missing_cts[is_case * max_window_size]);
	r_matrix_ptr = &(r_matrix[is_case]);
	geno_fixed_vec_ptr = &(window_geno_ptr[window_cidx * cur_192_long]);
	mask_fixed_vec_ptr = &(window_mask_ptr[window_cidx * cur_192_long]);
        copy_quaterarr_nonempty_subset(loadbuf_raw, founder_phenos[is_case], unfiltered_sample_ct, cur_pheno_ct, geno_fixed_vec_ptr);
        ld_process_load2(geno_fixed_vec_ptr, mask_fixed_vec_ptr, &fixed_missing_ct, cur_pheno_ct, is_x && (!ignore_x), pheno_male_include2[is_case]);
	fixed_non_missing_ct = cur_pheno_ct - fixed_missing_ct;
        missing_cts_ptr[window_cidx] = fixed_missing_ct;
	window_cidx3 = window_cidx2;
	while (window_cidx3 != window_cidx) {
	  geno_var_vec_ptr = &(window_geno_ptr[window_cidx3 * cur_192_long]);
	  mask_var_vec_ptr = &(window_mask_ptr[window_cidx3 * cur_192_long]);
	  non_missing_ct = fixed_non_missing_ct - missing_cts_ptr[window_cidx3];
	  if (fixed_missing_ct && missing_cts_ptr[window_cidx3]) {
            non_missing_ct += ld_missing_ct_intersect(mask_var_vec_ptr, mask_fixed_vec_ptr, cur_ctwd12, cur_ctwd12_rem, lshift_last);
	  }
	  if (non_missing_ct) {
	    dp_result[0] = cur_pheno_ct;
	    dp_result[1] = -((int32_t)fixed_non_missing_ct);
	    dp_result[2] = missing_cts_ptr[window_cidx3] - cur_pheno_ct;
	    dp_result[3] = dp_result[1];
	    dp_result[4] = dp_result[2];
	    ld_dot_prod(geno_var_vec_ptr, geno_fixed_vec_ptr, mask_var_vec_ptr, mask_fixed_vec_ptr, dp_result, cur_mld_m1, cur_mld_rem);
	    non_missing_ctd = (double)((int32_t)non_missing_ct);
            dxx = dp_result[1];
            dyy = dp_result[2];
            cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
            dxx = (dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy);
	    dxx = cov12 / sqrt(dxx);
	  } else {
	    dxx = 0.0;
	  }
	  r_matrix_ptr[2 * (window_cidx3 * max_window_size + window_cidx)] = dxx;
	  r_matrix_ptr[2 * (window_cidx * max_window_size + window_cidx3)] = dxx;
          if (++window_cidx3 == max_window_size) {
            window_cidx3 = 0;
	  }
	}
      }

      if (++chrom_marker_idx < chrom_marker_ct) {
        marker_uidx++;
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	    goto flipscan_ret_READ_FAIL;
	  }
	}
        marker_pos_thresh = marker_pos[marker_uidx];
	if (marker_pos_thresh < window_bp) {
	  marker_pos_thresh = 0;
	} else {
	  marker_pos_thresh -= window_bp;
	}
      } else {
	// close out the chromosome
        marker_pos_thresh = 0x80000000U;
      }
      // only need to enforce window locus count constraint during first loop
      // iteration
      ulii = window_cidx2 + max_window_locus_ct;
      if (ulii >= max_window_size) {
	ulii -= max_window_size;
      }
      marker_uidx2 = window_uidxs[window_cidx2];
      if ((ulii == window_cidx) || (marker_pos[marker_uidx2] < marker_pos_thresh)) {
	do {
	  pos_r_tot = 0.0;
	  neg_r_tot = 0.0;
	  pos_r_ct = 0;
	  neg_r_ct = 0;
          r_row_ptr = &(r_matrix[2 * max_window_size * window_cidx2]);
	  window_cidx3 = window_cidx_starts[window_cidx2];
          while (1) {
	    ctrl_pheno = r_row_ptr[2 * window_cidx3];
	    case_pheno = r_row_ptr[2 * window_cidx3 + 1];
	    if ((fabs(ctrl_pheno) >= min_corr) || (fabs(case_pheno) >= min_corr)) {
	      dxx = fabs(ctrl_pheno) + fabs(case_pheno);
	      if (case_pheno * ctrl_pheno >= 0.0) {
                pos_r_ct++;
		pos_r_tot += dxx;
	      } else {
		neg_uidx_buf[neg_r_ct++] = window_uidxs[window_cidx3];
		neg_r_tot += dxx;
	      }
	    }
	    if (window_cidx3 == window_cidx) {
	      break;
	    }
	    if (++window_cidx3 == max_window_size) {
              window_cidx3 = 0;
	    }
	  }
	  wptr_start2 = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr_start);
	  wptr_start2 = memseta(wptr_start2, 32, 3);
          wptr_start2 = uint32toa_w10x(marker_pos[marker_uidx2], ' ', wptr_start2);
	  wptr_start2 = fw_strcpy(4, marker_allele_ptrs[2 * marker_uidx2], wptr_start2);
	  *wptr_start2++ = ' ';
	  wptr = fw_strcpy(4, marker_allele_ptrs[2 * marker_uidx2 + 1], wptr_start2);
	  *wptr++ = ' ';
	  wptr = dtoa_g_wxp3x(1.0 - set_allele_freqs[marker_uidx2], 8, ' ', wptr);
          wptr = uint32toa_w6x(pos_r_ct, ' ', wptr);
	  if (!pos_r_ct) {
	    wptr = memcpya(wptr, "      NA", 8);
	  } else {
            wptr = dtoa_g_wxp3(pos_r_tot / ((int32_t)(pos_r_ct * 2)), 8, wptr);
	  }
          *wptr++ = ' ';
          wptr = uint32toa_w6x(neg_r_ct, ' ', wptr);
	  if (!neg_r_ct) {
	    wptr = memcpya(wptr, "      NA", 8);
	  } else {
	    wptr = dtoa_g_wxp3(neg_r_tot / ((int32_t)(neg_r_ct * 2)), 8, wptr);
	  }
	  *wptr++ = ' ';
          if (fwrite_checked(textbuf, wptr - textbuf, outfile)) {
	    goto flipscan_ret_WRITE_FAIL;
	  }
	  if (neg_r_ct) {
	    for (ulii = 0; ulii < neg_r_ct; ulii++) {
	      if (ulii) {
		putc_unlocked('|', outfile);
	      }
              fputs(&(marker_ids[neg_uidx_buf[ulii] * max_marker_id_len]), outfile);
	    }
	    problem_ct++;
	    if (verbose) {
	      window_cidx3 = window_cidx_starts[window_cidx2];
	      while (1) {
		ctrl_pheno = r_row_ptr[2 * window_cidx3];
		case_pheno = r_row_ptr[2 * window_cidx3 + 1];
		if ((fabs(ctrl_pheno) >= min_corr) || (fabs(case_pheno) >= min_corr)) {
		  marker_uidx3 = window_uidxs[window_cidx3];
		  wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx3 * max_marker_id_len]), wptr_start2);
		  wptr = memseta(wptr, 32, 3);
		  wptr = uint32toa_w10x(marker_pos[marker_uidx3], ' ', wptr);
		  wptr = fw_strcpy(4, marker_allele_ptrs[2 * marker_uidx3], wptr);
                  *wptr++ = ' ';
		  wptr = dtoa_g_wxp3x(case_pheno, 8, ' ', wptr);
		  wptr = dtoa_g_wxp3x(ctrl_pheno, 8, '\n', wptr);
		  if (fwrite_checked(textbuf, wptr - textbuf, outfile_verbose)) {
		    goto flipscan_ret_WRITE_FAIL;
		  }
		}
		if (window_cidx3 == window_cidx) {
		  break;
		}
		if (++window_cidx3 == max_window_size) {
		  window_cidx3 = 0;
		}
	      }
	    }
	  }
	  putc_unlocked('\n', outfile);
	  if (++marker_idx >= pct_thresh) {
	    if (pct > 10) {
	      putc_unlocked('\b', stdout);
	    }
	    pct = (marker_idx * 100LLU) / marker_ct;
	    if (pct < 100) {
	      printf("\b\b%" PRIuPTR "%%", pct);
	      fflush(stdout);
	      pct_thresh = ((++pct) * ((uint64_t)marker_ct)) / 100;
	    }
	  }
	  // better to perform this comparison first
	  if (window_cidx2 == window_cidx) {
	    if (++window_cidx2 == max_window_size) {
	      window_cidx2 = 0;
	    }
	    break;
	  }
	  if (++window_cidx2 == max_window_size) {
	    window_cidx2 = 0;
	  }
	  marker_uidx2 = window_uidxs[window_cidx2];
	} while (marker_pos[marker_uidx2] < marker_pos_thresh);
      }
    } while (chrom_marker_idx < chrom_marker_ct);
  }
  if (fclose_null(&outfile)) {
    goto flipscan_ret_WRITE_FAIL;
  }
  if (verbose) {
    if (fclose_null(&outfile_verbose)) {
      goto flipscan_ret_WRITE_FAIL;
    }
  }
  putc_unlocked('\r', stdout);
  // not actually possible to have exactly one problem variant, heh
  LOGPRINTF("--flip-scan%s: %u variants with at least one negative LD match.\n", verbose? " verbose" : "", problem_ct);
  if (verbose) {
    LOGPRINTFWW("Report written to %s ; neg-match details written to %s.verbose .\n", outname, outname);
  } else {
    LOGPRINTFWW("Report written to %s .\n", outname);
  }
  while (0) {
  flipscan_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  flipscan_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  flipscan_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  flipscan_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  flipscan_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_verbose);
  return retval;
}

// LD multithread globals
static uintptr_t* g_ld_geno1;
static uintptr_t* g_ld_geno2;
static uintptr_t* g_ld_geno_masks1;
static uintptr_t* g_ld_geno_masks2;
static uint32_t* g_ld_missing_cts1;
static uint32_t* g_ld_missing_cts2;
static uint32_t* g_ld_interval1;
static double* g_ld_results;
static float* g_ld_results_f;
static double* g_ld_set_allele_freqs;
static uintptr_t g_ld_idx1_block_size;
static uintptr_t g_ld_idx2_block_size;
static uintptr_t g_ld_idx2_block_start;
static uintptr_t g_ld_block_idx1;
static uintptr_t g_ld_marker_ct;
static uintptr_t g_ld_marker_ctm8;
static uintptr_t g_ld_founder_ct;
static uintptr_t g_ld_founder_ct_192_long;
static uint32_t g_ld_founder_ct_mld_m1;
static uint32_t g_ld_founder_ct_mld_rem;
static uint32_t g_ld_is_r2;
static uint32_t g_ld_thread_ct;

// with '--r2 dprime', males should be downweighted by a factor of 2 when
// considering two X chromosome variants, and by a factor of sqrt(2) when doing
// an inter-chromosome evaluation involving a single Xchr variant.  (The
// sqrt(2) factor is not implemented by PLINK 1.07, but the math compels its
// use.)
static uintptr_t* g_ld_sex_male;
static uintptr_t* g_ld_thread_wkspace;
static uint32_t g_ld_xstart1;
static uint32_t g_ld_xend1;
static uint32_t g_ld_xstart2;
static uint32_t g_ld_xend2;

static char g_ld_delimiter;
static uint32_t g_ld_plink_maxsnp;
static char* g_ld_marker_ids;
static Chrom_info* g_ld_chrom_info_ptr;
static uint32_t* g_ld_marker_pos;
static double* g_ld_marker_cms;
static uintptr_t* g_ld_marker_exclude_idx1;
static uintptr_t* g_ld_marker_exclude;
static char** g_ld_marker_allele_ptrs;
static uintptr_t g_ld_max_marker_id_len;
static uintptr_t g_ld_marker_uidx1;
static uintptr_t g_ld_uidx2_start;
static uintptr_t g_ld_marker_uidx2;
static uintptr_t g_ld_block_idx2;
static double g_ld_window_r2;
static uint32_t g_ld_is_first_block;
static uint32_t g_ld_is_inter_chr;
static uint32_t g_ld_prefix_len;
static uint32_t g_ld_keep_sign;
static uint32_t g_ld_modifier;

THREAD_RET_TYPE ld_block_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t thread_ct = g_ld_thread_ct;
  uintptr_t block_idx1_start = (tidx * g_ld_idx1_block_size) / thread_ct;
  uintptr_t block_idx1_end = ((tidx + 1) * g_ld_idx1_block_size) / thread_ct;
  uintptr_t marker_idx2_maxw = g_ld_marker_ctm8;
  uintptr_t founder_ct = g_ld_founder_ct;
  uintptr_t founder_ctwd = founder_ct / BITCT2;
  uintptr_t founder_ctwd12 = founder_ctwd / 12;
  uintptr_t founder_ctwd12_rem = founder_ctwd - (12 * founder_ctwd12);
  uintptr_t lshift_last = 2 * ((0x7fffffc0 - founder_ct) % BITCT2);
  uintptr_t founder_ct_192_long = g_ld_founder_ct_192_long;
  uintptr_t* geno1 = g_ld_geno1;
  uintptr_t* geno_masks1 = g_ld_geno_masks1;
  uint32_t* missing_cts1 = g_ld_missing_cts1;
  uint32_t* ld_interval1 = g_ld_interval1;
  uint32_t founder_ct_mld_m1 = g_ld_founder_ct_mld_m1;
  uint32_t founder_ct_mld_rem = g_ld_founder_ct_mld_rem;
  uint32_t is_r2 = g_ld_is_r2;
  uint32_t keep_sign = g_ld_keep_sign;
  double* results = g_ld_results;
  float* results_f = g_ld_results_f;
  double* rptr = nullptr;
  float* rptr_f = nullptr;
  int32_t dp_result[5];
  uintptr_t* geno_fixed_vec_ptr;
  uintptr_t* geno_var_vec_ptr;
  uintptr_t* mask_fixed_vec_ptr;
  uintptr_t* mask_var_vec_ptr;
  uintptr_t* geno2;
  uintptr_t* geno_masks2;
  uint32_t* missing_cts2;
  uintptr_t idx2_block_size;
  uintptr_t idx2_block_start;
  uintptr_t block_idx1;
  uintptr_t block_idx2;
  uintptr_t cur_block_idx2_end;
  double non_missing_ctd;
  double cov12;
  double dxx;
  double dyy;
  float non_missing_ctf;
  float cov12_f;
  float fxx;
  float fyy;
  uint32_t fixed_missing_ct;
  uint32_t fixed_non_missing_ct;
  uint32_t non_missing_ct;
  while (1) {
    idx2_block_size = g_ld_idx2_block_size;
    idx2_block_start = g_ld_idx2_block_start;
    geno2 = g_ld_geno2;
    geno_masks2 = g_ld_geno_masks2;
    missing_cts2 = g_ld_missing_cts2;
    for (block_idx1 = block_idx1_start; block_idx1 < block_idx1_end; block_idx1++) {
      fixed_non_missing_ct = ld_interval1[block_idx1 * 2]; // temporary redefine
      block_idx2 = fixed_non_missing_ct;
      cur_block_idx2_end = ld_interval1[block_idx1 * 2 + 1];
      if (block_idx2 < idx2_block_start) {
	if (cur_block_idx2_end <= idx2_block_start) {
	  continue;
	}
	block_idx2 = 0;
      } else {
	block_idx2 -= idx2_block_start;
	if (block_idx2 >= idx2_block_size) {
	  // nondecreasing, so we can safely exit
	  break;
	}
      }
      cur_block_idx2_end -= idx2_block_start;
      if (cur_block_idx2_end > idx2_block_size) {
	cur_block_idx2_end = idx2_block_size;
      }
      if (results) {
	rptr = &(results[block_idx1 * marker_idx2_maxw + block_idx2 + idx2_block_start - fixed_non_missing_ct]);
      } else {
	rptr_f = &(results_f[block_idx1 * marker_idx2_maxw + block_idx2 + idx2_block_start - fixed_non_missing_ct]);
      }
      fixed_missing_ct = missing_cts1[block_idx1];
      fixed_non_missing_ct = founder_ct - fixed_missing_ct;
      geno_fixed_vec_ptr = &(geno1[block_idx1 * founder_ct_192_long]);
      mask_fixed_vec_ptr = &(geno_masks1[block_idx1 * founder_ct_192_long]);
      for (; block_idx2 < cur_block_idx2_end; block_idx2++) {
	geno_var_vec_ptr = &(geno2[block_idx2 * founder_ct_192_long]);
	mask_var_vec_ptr = &(geno_masks2[block_idx2 * founder_ct_192_long]);
	non_missing_ct = fixed_non_missing_ct - missing_cts2[block_idx2];
	if (fixed_missing_ct && missing_cts2[block_idx2]) {
	  non_missing_ct += ld_missing_ct_intersect(mask_var_vec_ptr, mask_fixed_vec_ptr, founder_ctwd12, founder_ctwd12_rem, lshift_last);
	}
	dp_result[0] = founder_ct;
	dp_result[1] = -fixed_non_missing_ct;
	dp_result[2] = missing_cts2[block_idx2] - founder_ct;
	dp_result[3] = dp_result[1];
	dp_result[4] = dp_result[2];
	ld_dot_prod(geno_var_vec_ptr, geno_fixed_vec_ptr, mask_var_vec_ptr, mask_fixed_vec_ptr, dp_result, founder_ct_mld_m1, founder_ct_mld_rem);
	if (results) {
	  non_missing_ctd = (double)((int32_t)non_missing_ct);
	  dxx = dp_result[1];
	  dyy = dp_result[2];
	  cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
	  dxx = (dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy);
	  if (!is_r2) {
	    dxx = cov12 / sqrt(dxx);
	  } else if (!keep_sign) {
	    dxx = (cov12 * cov12) / dxx;
	  } else {
	    dxx = (fabs(cov12) * cov12) / dxx;
	  }
	  *rptr++ = dxx;
	} else {
	  non_missing_ctf = (float)((int32_t)non_missing_ct);
	  fxx = dp_result[1];
	  fyy = dp_result[2];
	  cov12_f = dp_result[0] * non_missing_ctf - fxx * fyy;
	  fxx = (dp_result[3] * non_missing_ctf + fxx * fxx) * (dp_result[4] * non_missing_ctf + fyy * fyy);
	  if (!is_r2) {
	    fxx = cov12_f / sqrt(fxx);
	  } else if (!keep_sign) {
	    fxx = (cov12_f * cov12_f) / fxx;
	  } else {
	    fxx = (fabs(cov12_f) * cov12_f) / fxx;
	  }
	  *rptr_f++ = fxx;
	}
      }
    }
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

uint32_t ld_matrix_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t block_size1 = g_ld_idx1_block_size;
  uintptr_t marker_ct = g_ld_marker_ct;
  uintptr_t marker_ctm8 = g_ld_marker_ctm8;
  uintptr_t block_idx1 = g_ld_block_idx1;
  uintptr_t marker_idx = g_ld_idx2_block_start;
  uintptr_t marker_idx_end = g_ld_idx2_block_size;
  uint32_t is_square = ((g_ld_modifier & LD_MATRIX_SHAPEMASK) == LD_MATRIX_SQ);
  uint32_t is_square0 = ((g_ld_modifier & LD_MATRIX_SHAPEMASK) == LD_MATRIX_SQ0);
  char delimiter = g_ld_delimiter;
  double* results = g_ld_results;
  double* dptr;
  uintptr_t ulii;
  while (block_idx1 < block_size1) {
    dptr = &(results[block_idx1 * marker_ctm8 + marker_idx]);
    while (marker_idx < marker_idx_end) {
      sptr_cur = dtoa_gx(*dptr++, delimiter, sptr_cur);
      marker_idx++;
      if (sptr_cur > readbuf_end) {
	goto ld_matrix_emitn_ret;
      }
    }
    if (is_square0 && (marker_idx < marker_ct)) {
      ulii = (((uintptr_t)(readbuf_end - sptr_cur)) + 1) / 2;
      // bugfix: can't be <= since tab delimiter wouldn't be handled correctly
      // on subsequent pass
      if (ulii < marker_ct - marker_idx) {
	sptr_cur = memcpya(sptr_cur, g_textbuf, ulii * 2);
	marker_idx += ulii;
	goto ld_matrix_emitn_ret;
      } else {
	sptr_cur = memcpya(sptr_cur, g_textbuf, (marker_ct - marker_idx) * 2);
	marker_idx = marker_ct;
      }
    }
    if (delimiter == '\t') {
      sptr_cur--;
    }
    *sptr_cur++ = '\n';
    marker_idx = 0;
    if (!is_square) {
      marker_idx_end++;
    }
    block_idx1++;
  }
 ld_matrix_emitn_ret:
  g_ld_block_idx1 = block_idx1;
  g_ld_idx2_block_start = marker_idx;
  g_ld_idx2_block_size = marker_idx_end;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t ld_report_matrix(pthread_t* threads, Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, uintptr_t* sex_male, uintptr_t* founder_include2, uintptr_t* founder_male_include2, uintptr_t* loadbuf, char* outname, uint32_t hh_exists) {
  FILE* outfile = nullptr;
  uint32_t ld_modifier = ldip->modifier;
  uint32_t is_binary = ld_modifier & (LD_MATRIX_BIN | LD_MATRIX_BIN4);
  uint32_t is_square = ((ld_modifier & LD_MATRIX_SHAPEMASK) == LD_MATRIX_SQ);
  uint32_t is_square0 = ((ld_modifier & LD_MATRIX_SHAPEMASK) == LD_MATRIX_SQ0);
  uint32_t output_single_prec = (ld_modifier / LD_MATRIX_BIN4) & 1;
  uint32_t output_gz = ld_modifier & LD_REPORT_GZ;
  uint32_t ignore_x = (ld_modifier / LD_IGNORE_X) & 1;
  uintptr_t marker_ct = g_ld_marker_ct;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t marker_ctm8 = round_up_pow2(marker_ct, 8);
  uintptr_t founder_ct = g_ld_founder_ct;
  uintptr_t founder_ctl = BITCT_TO_WORDCT(founder_ct);
  uintptr_t founder_ct_192_long = g_ld_founder_ct_192_long;
  uintptr_t final_mask = get_final_mask(founder_ct);
  uintptr_t marker_uidx_base = next_unset_unsafe(marker_exclude, 0);
  uintptr_t marker_uidx1 = marker_uidx_base;
  uintptr_t marker_idx1_start = (((uint64_t)parallel_idx) * marker_ct) / parallel_tot;
  uintptr_t marker_idx1 = marker_idx1_start;
  uintptr_t marker_idx1_end = (((uint64_t)(parallel_idx + 1)) * marker_ct) / parallel_tot;
  uintptr_t pct = 1;
  uint64_t job_size = marker_idx1_end - marker_idx1_start;
  uint64_t pct_thresh = job_size / 100;
  Chrom_info* chrom_info_ptr = g_ld_chrom_info_ptr;
  uint32_t founder_trail_ct = founder_ct_192_long - founder_ctl * 2;
  uint32_t thread_ct = g_ld_thread_ct;
  uint32_t chrom_fo_idx = 0;
  uint32_t is_haploid = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t not_first_write = 0;
  int32_t retval = 0;
  unsigned char* bigstack_mark2;
  uintptr_t* ulptr;
  unsigned char* overflow_buf;
  uint64_t tests_completed;
  uintptr_t thread_workload;
  uintptr_t cur_idx2_block_size;
  uintptr_t marker_idx2_end;
  uintptr_t block_idx1;
  uintptr_t marker_uidx2;
  uintptr_t marker_idx2;
  uintptr_t block_idx2;
  uintptr_t idx1_block_size;
  uintptr_t idx2_block_size;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t is_last_block;

  if (bigstack_alloc_uc(262144, &overflow_buf)) {
    goto ld_report_matrix_ret_NOMEM;
  }
  if (output_single_prec) {
    // force divisibility by 16 instead (cacheline = 64 bytes, float = 4)
    marker_ctm8 = (marker_ctm8 + 8) & (~15);
  }
  if (is_binary) {
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto ld_report_matrix_ret_OPEN_FAIL;
    }
  }
  // claim up to half of memory with idx1 bufs; each marker costs
  //   founder_ct_192_long * sizeof(intptr_t) for genotype buffer
  // + founder_ct_192_long * sizeof(intptr_t) for missing mask buffer
  // + sizeof(int32_t) for g_ld_missing_cts1 entry
  // + 2 * sizeof(int32_t) for g_ld_interval1
  // + marker_ctm8 * sizeof(double) or marker_ctm16 * sizeof(float) for
  //     g_ld_results buffer
  // round down to multiple of thread_ct for better workload distribution
  ulii = founder_ct_192_long * 2 * sizeof(intptr_t) + 3 * sizeof(int32_t) + marker_ctm8 * (8 - 4 * output_single_prec);
  idx1_block_size = bigstack_left() / (ulii * 2);
  thread_workload = idx1_block_size / thread_ct;
  if (!thread_workload) {
    goto ld_report_matrix_ret_NOMEM;
  }
  idx1_block_size = thread_workload * thread_ct;
  if ((parallel_tot > 1) && (marker_ct < 2 * parallel_tot)) {
    LOGERRPRINTF("Error: Too few variants in --r%s run for --parallel %u %u.\n", g_ld_is_r2? "2" : "", parallel_idx + 1, parallel_tot);
    goto ld_report_matrix_ret_INVALID_CMDLINE;
  }
  if (!is_square) {
    job_size = ((uint64_t)marker_ct) * (marker_ct + 1);
    if (parallel_tot > 1) {
      job_size /= parallel_tot;
      marker_idx1_start = triangle_divide(job_size * parallel_idx, 1);
      if (parallel_idx + 1 < parallel_tot) {
        marker_idx1_end = triangle_divide(job_size * (parallel_idx + 1), 1);
      }
      job_size = ((((uint64_t)marker_idx1_end) * (marker_idx1_end + 1)) - (((uint64_t)marker_idx1_start) * (marker_idx1_start + 1))) / 2;
    } else {
      job_size /= 2;
    }
  }
  pct_thresh = job_size / 100;
  if (idx1_block_size > marker_idx1_end - marker_idx1_start) {
    idx1_block_size = marker_idx1_end - marker_idx1_start;
  }
  bigstack_alloc_ul(founder_ct_192_long * idx1_block_size, &g_ld_geno1);
  bigstack_alloc_ul(founder_ct_192_long * idx1_block_size, &g_ld_geno_masks1);
  bigstack_alloc_ui(idx1_block_size, &g_ld_missing_cts1);
  bigstack_alloc_ui(idx1_block_size * 2, &g_ld_interval1);

  if (!output_single_prec) {
    // may want to set g_ld_results_f to nullptr
    if (bigstack_alloc_d(marker_ctm8 * idx1_block_size, &g_ld_results)) {
      goto ld_report_matrix_ret_NOMEM;
    }
  } else {
    g_ld_results = nullptr;
    if (bigstack_alloc_f(marker_ctm8 * idx1_block_size, &g_ld_results_f)) {
      goto ld_report_matrix_ret_NOMEM;
    }
  }

  // claim the other half with idx2 buffer
  ulii -= marker_ctm8 * (8 - 4 * output_single_prec) + 2 * sizeof(int32_t);
  if (!output_single_prec) {
    idx2_block_size = (bigstack_left() / ulii) & (~(7 * ONELU));
  } else {
    idx2_block_size = (bigstack_left() / ulii) & (~(15 * ONELU));
  }
  if (idx2_block_size > marker_ctm8) {
    idx2_block_size = marker_ctm8;
  }
  bigstack_mark2 = g_bigstack_base;
  while (1) {
    if (!idx2_block_size) {
      goto ld_report_matrix_ret_NOMEM;
    }
    if (!(bigstack_alloc_ul(founder_ct_192_long * idx2_block_size, &g_ld_geno2) ||
          bigstack_alloc_ul(founder_ct_192_long * idx2_block_size, &g_ld_geno_masks2) ||
          bigstack_alloc_ui(idx2_block_size, &g_ld_missing_cts2))) {
      break;
    }
    bigstack_reset(bigstack_mark2);
    if (!output_single_prec) {
      idx2_block_size -= 8;
    } else {
      idx2_block_size -= 16;
    }
  }
  uljj = founder_trail_ct + 2;
  for (ulii = 1; ulii <= idx1_block_size; ulii++) {
    fill_ulong_zero(uljj, &(g_ld_geno1[ulii * founder_ct_192_long - uljj]));
    fill_ulong_zero(uljj, &(g_ld_geno_masks1[ulii * founder_ct_192_long - uljj]));
  }
  for (ulii = 1; ulii <= idx2_block_size; ulii++) {
    fill_ulong_zero(uljj, &(g_ld_geno2[ulii * founder_ct_192_long - uljj]));
    fill_ulong_zero(uljj, &(g_ld_geno_masks2[ulii * founder_ct_192_long - uljj]));
  }
  if (is_square) {
    for (ulii = 0; ulii < idx1_block_size; ulii++) {
      g_ld_interval1[ulii * 2] = 0;
      g_ld_interval1[ulii * 2 + 1] = marker_ct;
    }
    g_ld_marker_ctm8 = marker_ctm8;
  } else {
    for (ulii = 0; ulii < idx1_block_size; ulii++) {
      g_ld_interval1[ulii * 2] = 0;
    }
    if (is_square0) {
      if (is_binary) {
	if (!output_single_prec) {
          fill_double_zero(MAXLINELEN / sizeof(double), (double*)g_textbuf);
	} else {
          fill_float_zero(MAXLINELEN / sizeof(float), (float*)g_textbuf);
	}
      } else {
	ulptr = (uintptr_t*)g_textbuf;
	// assume little-endian
	// 0[delim]0[delim]...
#ifdef __LP64__
	ulii = 0x30003000300030LLU | (0x100010001000100LLU * ((unsigned char)g_ld_delimiter));
#else
	ulii = 0x300030 | (0x1000100 * ((unsigned char)g_ld_delimiter));
#endif
        for (uljj = 0; uljj < MAXLINELEN / sizeof(intptr_t); uljj++) {
	  *ulptr++ = ulii;
	}
      }
    }
  }
  if (marker_idx1) {
    marker_uidx1 = jump_forward_unset_unsafe(marker_exclude, marker_uidx1 + 1, marker_idx1);
  }
  g_ld_keep_sign = 0;
  sprintf(g_logbuf, "--r%s %s%s to %s ... ", g_ld_is_r2? "2" : "", is_square? "square" : (is_square0? "square0" : "triangle"), is_binary? (output_single_prec? " bin4" : " bin") : (output_gz? " gz" : ""), outname);
  wordwrapb(16); // strlen("99% [processing]")
  logprintb();
  fputs("0%", stdout);
  do {
    fputs(" [processing]", stdout);
    fflush(stdout);
    if (idx1_block_size > marker_idx1_end - marker_idx1) {
      idx1_block_size = marker_idx1_end - marker_idx1;
      if (idx1_block_size < thread_ct) {
        thread_ct = idx1_block_size;
        g_ld_thread_ct = thread_ct;
      }
    }
    g_ld_idx1_block_size = idx1_block_size;
    // marker_uidx1_tmp = marker_uidx1;
    if (fseeko(bedfile, bed_offset + (marker_uidx1 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto ld_report_matrix_ret_READ_FAIL;
    }
    chrom_end = 0;
    for (block_idx1 = 0; block_idx1 < idx1_block_size; marker_uidx1++, block_idx1++) {
      if (IS_SET(marker_exclude, marker_uidx1)) {
        marker_uidx1 = next_unset_ul_unsafe(marker_exclude, marker_uidx1);
        if (fseeko(bedfile, bed_offset + (marker_uidx1 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	  goto ld_report_matrix_ret_READ_FAIL;
	}
      }
      if (marker_uidx1 >= chrom_end) {
        chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx1);
        chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
        is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
	is_x = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[X_OFFSET]);
	is_y = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
      }
      if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, marker_uidx1), bedfile, loadbuf, &(g_ld_geno1[block_idx1 * founder_ct_192_long]))) {
	goto ld_report_matrix_ret_READ_FAIL;
      }
      if (is_haploid && hh_exists) {
	haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(g_ld_geno1[block_idx1 * founder_ct_192_long])));
      }
      ld_process_load2(&(g_ld_geno1[block_idx1 * founder_ct_192_long]), &(g_ld_geno_masks1[block_idx1 * founder_ct_192_long]), &(g_ld_missing_cts1[block_idx1]), founder_ct, is_x && (!ignore_x), founder_male_include2);
    }
    marker_uidx2 = marker_uidx_base;
    marker_idx2 = 0;
    if (is_square) {
      marker_idx2_end = marker_ct;
    } else {
      marker_idx2_end = marker_idx1 + idx1_block_size;
      for (ulii = 1; ulii <= idx1_block_size; ulii++) {
	g_ld_interval1[2 * ulii - 1] = ulii + marker_idx1;
      }
      if (!output_single_prec) {
        marker_ctm8 = round_up_pow2(marker_idx2_end, 8);
      } else {
        marker_ctm8 = round_up_pow2(marker_idx2_end, 16);
      }
      g_ld_marker_ctm8 = marker_ctm8;
    }
    chrom_end = 0;
    if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto ld_report_matrix_ret_READ_FAIL;
    }
    cur_idx2_block_size = idx2_block_size;
    do {
      if (cur_idx2_block_size > marker_idx2_end - marker_idx2) {
	cur_idx2_block_size = marker_idx2_end - marker_idx2;
      }
      for (block_idx2 = 0; block_idx2 < cur_idx2_block_size; marker_uidx2++, block_idx2++) {
	if (IS_SET(marker_exclude, marker_uidx2)) {
          marker_uidx2 = next_unset_ul_unsafe(marker_exclude, marker_uidx2);
	  if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	    goto ld_report_matrix_ret_READ_FAIL;
	  }
	}
	if (marker_uidx2 >= chrom_end) {
	  chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx2);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
	  is_x = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[X_OFFSET]);
	  is_y = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
	}
	if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, marker_uidx2), bedfile, loadbuf, &(g_ld_geno2[block_idx2 * founder_ct_192_long]))) {
	  goto ld_report_matrix_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(g_ld_geno2[block_idx2 * founder_ct_192_long])));
	}
	ld_process_load2(&(g_ld_geno2[block_idx2 * founder_ct_192_long]), &(g_ld_geno_masks2[block_idx2 * founder_ct_192_long]), &(g_ld_missing_cts2[block_idx2]), founder_ct, is_x && (!ignore_x), founder_male_include2);
      }
      g_ld_idx2_block_size = cur_idx2_block_size;
      g_ld_idx2_block_start = marker_idx2;
      marker_idx2 += cur_idx2_block_size;
      is_last_block = (marker_idx2 >= marker_idx2_end);
      if (spawn_threads2(threads, &ld_block_thread, thread_ct, is_last_block)) {
	goto ld_report_matrix_ret_THREAD_CREATE_FAIL;
      }
      ld_block_thread((void*)0);
      join_threads2(threads, thread_ct, is_last_block);
    } while (!is_last_block);
    fputs("\b\b\b\b\b\b\b\b\b\b\bwriting]   \b\b\b", stdout);
    fflush(stdout);
    if (is_binary) {
      if (!output_single_prec) {
	if (is_square) {
	  for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
	    if (fwrite_checked(&(g_ld_results[block_idx1 * marker_ctm8]), cur_idx2_block_size * sizeof(double), outfile)) {
	      goto ld_report_matrix_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
	    if (fwrite_checked(&(g_ld_results[block_idx1 * marker_ctm8]), (block_idx1 + marker_idx1 + 1) * sizeof(double), outfile)) {
	      goto ld_report_matrix_ret_WRITE_FAIL;
	    }
	    if (is_square0) {
	      ulii = marker_ct - block_idx1 - marker_idx1 - 1;
	      while (ulii) {
		if (ulii > MAXLINELEN / sizeof(double)) {
		  uljj = MAXLINELEN / sizeof(double);
		  ulii -= MAXLINELEN / sizeof(double);
		} else {
		  uljj = ulii;
		  ulii = 0;
		}
		if (fwrite_checked(g_textbuf, uljj * sizeof(double), outfile)) {
		  goto ld_report_matrix_ret_WRITE_FAIL;
		}
	      }
	    }
	  }
	}
      } else {
	if (is_square) {
	  for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
	    if (fwrite_checked(&(g_ld_results_f[block_idx1 * marker_ctm8]), cur_idx2_block_size * sizeof(float), outfile)) {
	      goto ld_report_matrix_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
	    if (fwrite_checked(&(g_ld_results_f[block_idx1 * marker_ctm8]), (block_idx1 + marker_idx1 + 1) * sizeof(float), outfile)) {
	      goto ld_report_matrix_ret_WRITE_FAIL;
	    }
	    if (is_square0) {
	      ulii = marker_ct - block_idx1 - marker_idx1 - 1;
	      while (ulii) {
		if (ulii > MAXLINELEN / sizeof(float)) {
		  uljj = MAXLINELEN / sizeof(float);
		  ulii -= MAXLINELEN / sizeof(float);
		} else {
		  uljj = ulii;
		  ulii = 0;
		}
		if (fwrite_checked(g_textbuf, uljj * sizeof(float), outfile)) {
		  goto ld_report_matrix_ret_WRITE_FAIL;
		}
	      }
	    }
	  }
	}
      }
    } else {
      g_ld_block_idx1 = 0;
      g_ld_idx2_block_start = 0;
      if (is_square) {
        g_ld_idx2_block_size = marker_ct;
      } else {
	g_ld_idx2_block_size = marker_idx1 + 1;
      }
      if (output_gz) {
        parallel_compress(outname, overflow_buf, not_first_write, ld_matrix_emitn);
      } else {
        write_uncompressed(outname, overflow_buf, not_first_write, ld_matrix_emitn);
      }
      not_first_write = 1;
    }
    marker_idx1 += idx1_block_size;
    fputs("\b\b\b\b\b\b\b\b\b\b          \b\b\b\b\b\b\b\b\b\b", stdout);
    if (is_square) {
      tests_completed = marker_idx1 - marker_idx1_start;
    } else {
      tests_completed = ((((uint64_t)marker_idx1) * (marker_idx1 + 1)) - (((uint64_t)marker_idx1_start) * (marker_idx1_start + 1))) / 2;
    }
    if (tests_completed >= pct_thresh) {
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      pct = (tests_completed * 100LLU) / job_size;
      if (pct < 100) {
	printf("\b\b%" PRIuPTR "%%", pct);
	fflush(stdout);
	pct_thresh = ((++pct) * ((uint64_t)job_size)) / 100;
      }
    }
  } while (marker_idx1 < marker_idx1_end);
  fputs("\b\b", stdout);
  logprint("done.\n");
  if (is_binary) {
    if (fclose_null(&outfile)) {
      goto ld_report_matrix_ret_WRITE_FAIL;
    }
  }
  while (0) {
  ld_report_matrix_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_report_matrix_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ld_report_matrix_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_report_matrix_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  ld_report_matrix_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  ld_report_matrix_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
  fclose_cond(outfile);
  // trust parent to free memory
  return retval;
}

uint32_t ld_regular_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  Chrom_info* chrom_info_ptr = g_ld_chrom_info_ptr;
  uintptr_t* marker_exclude_idx1 = g_ld_marker_exclude_idx1;
  uintptr_t* marker_exclude = g_ld_marker_exclude;
  uint32_t* marker_pos = g_ld_marker_pos;
  char* marker_ids = g_ld_marker_ids;
  char** marker_allele_ptrs = g_ld_marker_allele_ptrs;
  uint32_t* ld_interval1 = g_ld_interval1;
  double* results = g_ld_results;
  double* set_allele_freqs = g_ld_set_allele_freqs;
  char* fixed_a1 = nullptr;
  char* fixed_a2 = nullptr;
  uintptr_t max_marker_id_len = g_ld_max_marker_id_len;
  uintptr_t marker_uidx1 = g_ld_marker_uidx1;
  uintptr_t block_idx1 = g_ld_block_idx1;
  uintptr_t block_size1 = g_ld_idx1_block_size;
  uintptr_t marker_uidx2_start = g_ld_uidx2_start;
  uintptr_t block_idx2_start = g_ld_idx2_block_start;
  uintptr_t block_idx2 = g_ld_block_idx2;
  uintptr_t marker_idx2_maxw = g_ld_marker_ctm8;
  uintptr_t marker_uidx2 = g_ld_marker_uidx2;
  double window_r2 = g_ld_window_r2;
  uint32_t plink_maxsnp = g_ld_plink_maxsnp;
  uint32_t is_inter_chr = g_ld_is_inter_chr;

  // 0 = not d/dprime/dprime-signed
  uint32_t dprime_type = g_ld_modifier & LD_DX;

  uint32_t is_r2 = g_ld_is_r2;
  uint32_t prefix_len = g_ld_prefix_len;
  uint32_t chrom_fo_idx1 = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx1);
  uint32_t chrom_idx1 = chrom_info_ptr->chrom_file_order[chrom_fo_idx1];
  uint32_t chrom_end1 = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx1 + 1];
  uint32_t chrom_fo_idx2 = 0;
  uint32_t chrom_idx2 = 0;
  uint32_t fixed_a1_len = 0;
  uint32_t fixed_a2_len = 0;
  uintptr_t block_end2;
  uint32_t coupling;
  uint32_t chrom_end2;
  char* sptr2;
  double* dptr;
  double dxx;
  if (block_idx1 == block_size1) {
    goto ld_regular_emitn_ret;
  }
  if (block_idx2) {
    goto ld_regular_emitn_start_2;
  }
  // block_idx2 is only zero on initial call, never on reentry
  if (g_ld_is_first_block) {
    sptr_cur = memcpya(sptr_cur, " CHR_A         BP_A ", 20);
    sptr_cur = fw_strcpyn(g_ld_plink_maxsnp, 5, "SNP_A", sptr_cur);
    if (set_allele_freqs) {
      sptr_cur = memcpya(sptr_cur, "      MAF_A", 11);
    }
    sptr_cur = memcpya(sptr_cur, "  CHR_B         BP_B ", 21);
    sptr_cur = fw_strcpyn(g_ld_plink_maxsnp, 5, "SNP_B", sptr_cur);
    if (marker_allele_ptrs) {
      sptr_cur = memcpya(sptr_cur, "      PHASE", 11);
    }
    if (set_allele_freqs) {
      sptr_cur = memcpya(sptr_cur, "      MAF_B", 11);
    }
    sptr_cur = memseta(sptr_cur, 32, 11);
    sptr_cur = memcpyl3a(sptr_cur, is_r2? "R2 " : " R ");
    if (dprime_type) {
      sptr_cur = memcpya(sptr_cur, (dprime_type == LD_D)? "           D " : "          DP ", 13);
    }
    *sptr_cur++ = '\n';
  }
  goto ld_regular_emitn_start;
  do {
    marker_uidx1++;
    next_unset_ul_unsafe_ck(marker_exclude_idx1, &marker_uidx1);
    if (marker_uidx1 >= chrom_end1) {
      chrom_fo_idx1 = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx1);
      chrom_idx1 = chrom_info_ptr->chrom_file_order[chrom_fo_idx1];
      chrom_end1 = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx1 + 1];
    }
    block_idx2 = ld_interval1[2 * block_idx1];
    if (block_idx2_start < block_idx2) {
      marker_uidx2_start = jump_forward_unset_unsafe(marker_exclude, marker_uidx2_start + 1, block_idx2 - block_idx2_start);
      block_idx2_start = block_idx2;
    }
  ld_regular_emitn_start:
    marker_uidx2 = marker_uidx2_start;
    sptr2 = width_force(6, g_textbuf, chrom_name_write(chrom_info_ptr, chrom_idx1, g_textbuf));
    sptr2 = memseta(sptr2, 32, 3);
    sptr2 = uint32toa_w10x(marker_pos[marker_uidx1], ' ', sptr2);
    sptr2 = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx1 * max_marker_id_len]), sptr2);
    *sptr2++ = ' ';
    if (set_allele_freqs) {
      sptr2 = width_force(10, sptr2, dtoa_g(1.0 - set_allele_freqs[marker_uidx1], sptr2));
      *sptr2++ = ' ';
    }
    if (!is_inter_chr) {
      sptr2 = width_force(6, sptr2, chrom_name_write(chrom_info_ptr, chrom_idx1, sptr2));
      sptr2 = memseta(sptr2, 32, 3);
    }
    prefix_len = (uintptr_t)(sptr2 - g_textbuf);
  ld_regular_emitn_start_2:
    if (marker_allele_ptrs) {
      fixed_a1 = marker_allele_ptrs[2 * marker_uidx1];
      fixed_a2 = marker_allele_ptrs[2 * marker_uidx1 + 1];
      fixed_a1_len = strlen(fixed_a1);
      fixed_a2_len = strlen(fixed_a2);
    }
    chrom_end2 = 0;
    block_end2 = ld_interval1[2 * block_idx1 + 1];
    dptr = &(results[(block_idx1 * marker_idx2_maxw + block_idx2 - block_idx2_start) * (1 + (dprime_type != 0))]);
    while (block_idx2 < block_end2) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx2);
      dxx = *dptr++;
      if (fabs(dxx) >= window_r2) {
	sptr_cur = memcpya(sptr_cur, g_textbuf, prefix_len);
	if (is_inter_chr) {
	  if (marker_uidx2 >= chrom_end2) {
	    chrom_fo_idx2 = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx2);
	    chrom_idx2 = chrom_info_ptr->chrom_file_order[chrom_fo_idx2];
	    chrom_end2 = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx2 + 1];
	  }
	  sptr_cur = width_force(6, sptr_cur, chrom_name_write(chrom_info_ptr, chrom_idx2, sptr_cur));
	  sptr_cur = memseta(sptr_cur, 32, 3);
	}
	sptr_cur = uint32toa_w10x(marker_pos[marker_uidx2], ' ', sptr_cur);
	sptr_cur = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), sptr_cur);
	*sptr_cur++ = ' ';
	if (marker_allele_ptrs) {
	  coupling = (dxx > 0);
	  sptr2 = memcpya(sptr_cur, fixed_a1, fixed_a1_len);
	  sptr2 = strcpyax(sptr2, marker_allele_ptrs[2 * marker_uidx2 + (1 - coupling)], '/');
	  sptr2 = memcpya(sptr2, fixed_a2, fixed_a2_len);
          sptr2 = strcpya(sptr2, marker_allele_ptrs[2 * marker_uidx2 + coupling]);
	  sptr_cur = width_force(10, sptr_cur, sptr2);
	  *sptr_cur++ = ' ';
	}
	if (set_allele_freqs) {
	  sptr_cur = width_force(10, sptr_cur, dtoa_g(1.0 - set_allele_freqs[marker_uidx2], sptr_cur));
	  *sptr_cur++ = ' ';
	}
	if (is_r2) {
	  dxx = fabs(dxx);
	}
	sptr_cur = width_force(12, sptr_cur, dtoa_g(dxx, sptr_cur));
	if (dprime_type) {
	  *sptr_cur++ = ' ';
          sptr_cur = width_force(12, sptr_cur, dtoa_g(*dptr++, sptr_cur));
	}
	sptr_cur = memcpya(sptr_cur, " \n", 2);
      } else if (dprime_type) {
	dptr++;
      }
      block_idx2++;
      marker_uidx2++;
      if (sptr_cur > readbuf_end) {
        goto ld_regular_emitn_ret;
      }
    }
  } while (++block_idx1 < block_size1);
 ld_regular_emitn_ret:
  g_ld_marker_uidx1 = marker_uidx1;
  g_ld_block_idx1 = block_idx1;
  g_ld_prefix_len = prefix_len;
  g_ld_uidx2_start = marker_uidx2_start;
  g_ld_idx2_block_start = block_idx2_start;
  g_ld_marker_uidx2 = marker_uidx2;
  g_ld_block_idx2 = block_idx2;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

// The following three functions are built around a data representation
// introduced by Xiang Yan et al.'s BOOST software (the original bitwise
// representation I came up with was less efficient); see
// http://bioinformatics.ust.hk/BOOST.html .
//
// The BOOST implementation just evaluated four contingency table values; when
// there is no missing data, the other five can be determined via subtraction.
// two_locus_3x3_zmiss_tablev() function handles this case.  However, with
// *only* that logic, all sites with missing data must be thrown out.
// two_locus_3x3_tablev() handles the other cases, directly summing 6 or 9
// table values when necessary.
//
// If permutation testing is added later, it should exploit the fact that
// [cell xy value in case 3x3 table] + [cell xy value in ctrl 3x3 table]
// is constant across permutations; i.e. we just need to determine the new case
// contingency table, and then the control table falls out via subtraction.
// Several ideas from PERMORY could also be applied.
uint32_t load_and_split3(FILE* bedfile, uintptr_t* rawbuf, uint32_t unfiltered_sample_ct, uintptr_t* casebuf, uintptr_t* pheno_nm, uintptr_t* pheno_c, uint32_t case_ctv, uint32_t ctrl_ctv, uint32_t do_reverse, uint32_t is_case_only, uintptr_t* nm_info_ptr) {
  uintptr_t* rawbuf_end = &(rawbuf[unfiltered_sample_ct / BITCT2]);
  uintptr_t* ctrlbuf = &(casebuf[3 * case_ctv]);
  uintptr_t case_words[4];
  uintptr_t ctrl_words[4];
  uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uint32_t case_rem = 0;
  uint32_t ctrl_rem = 0;
  uint32_t read_shift_max = BITCT2;
  uint32_t sample_uidx = 0;
  uint32_t offset0_case = do_reverse * 2 * case_ctv;
  uint32_t offset2_case = (1 - do_reverse) * 2 * case_ctv;
  uint32_t offset0_ctrl = do_reverse * 2 * ctrl_ctv;
  uint32_t offset2_ctrl = (1 - do_reverse) * 2 * ctrl_ctv;
  uint32_t read_shift;
  uintptr_t read_word;
  uintptr_t ulii;
  if (bedfile) {
    // ld_report_dprime() preloads this and does het. haploid handling, etc.
    if (load_raw(unfiltered_sample_ct4, bedfile, rawbuf)) {
      return RET_READ_FAIL;
    }
  }
  case_words[0] = 0;
  case_words[1] = 0;
  case_words[2] = 0;
  case_words[3] = 0;
  ctrl_words[0] = 0;
  ctrl_words[1] = 0;
  ctrl_words[2] = 0;
  ctrl_words[3] = 0;
  while (1) {
    while (rawbuf < rawbuf_end) {
      read_word = *rawbuf++;
      for (read_shift = 0; read_shift < read_shift_max; sample_uidx++, read_shift++) {
	if (is_set(pheno_nm, sample_uidx)) {
	  ulii = read_word & 3;
	  if (is_set(pheno_c, sample_uidx)) {
	    case_words[ulii] |= ONELU << case_rem;
	    if (++case_rem == BITCT) {
	      casebuf[offset0_case] = case_words[0];
	      casebuf[case_ctv] = case_words[2];
	      casebuf[offset2_case] = case_words[3];
	      casebuf++;
	      case_words[0] = 0;
	      case_words[2] = 0;
	      case_words[3] = 0;
	      case_rem = 0;
	    }
	  } else if (!is_case_only) {
	    ctrl_words[ulii] |= ONELU << ctrl_rem;
	    if (++ctrl_rem == BITCT) {
	      ctrlbuf[offset0_ctrl] = ctrl_words[0];
	      ctrlbuf[ctrl_ctv] = ctrl_words[2];
	      ctrlbuf[offset2_ctrl] = ctrl_words[3];
	      ctrlbuf++;
	      ctrl_words[0] = 0;
	      ctrl_words[2] = 0;
	      ctrl_words[3] = 0;
	      ctrl_rem = 0;
	    }
	  }
	}
	read_word >>= 2;
      }
    }
    if (sample_uidx == unfiltered_sample_ct) {
      if (case_rem) {
	casebuf[offset0_case] = case_words[0];
	casebuf[case_ctv] = case_words[2];
	casebuf[offset2_case] = case_words[3];
      }
      if (ctrl_rem) {
	ctrlbuf[offset0_ctrl] = ctrl_words[0];
	ctrlbuf[ctrl_ctv] = ctrl_words[2];
	ctrlbuf[offset2_ctrl] = ctrl_words[3];
      }
      ulii = 3;
      if (case_words[1]) {
	ulii -= 1;
      }
      if (ctrl_words[1]) {
	ulii -= 2;
      }
      *nm_info_ptr = ulii;
      return 0;
    }
    rawbuf_end++;
    read_shift_max = unfiltered_sample_ct % BITCT2;
  }
}

#ifdef __LP64__
static void two_locus_3x3_tablev(__m128i* vec1, __m128i* vec2, uint32_t* counts_3x3, uint32_t sample_ctv6, uint32_t iter_ct) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i* vec20;
  __m128i* vec21;
  __m128i* vec22;
  __m128i* vend1;
  __m128i loader1;
  __m128i loader20;
  __m128i loader21;
  __m128i loader22;
  __m128i count10;
  __m128i count11;
  __m128i count12;
  __m128i count20;
  __m128i count21;
  __m128i count22;
  __univec acc0;
  __univec acc1;
  __univec acc2;
  uint32_t ct;
  uint32_t ct2;
  while (iter_ct--) {
    ct = sample_ctv6;
    vec20 = vec2;
    vec21 = &(vec20[sample_ctv6]);
    vec22 = &(vec20[2 * sample_ctv6]);
    while (ct >= 30) {
      ct -= 30;
      vend1 = &(vec1[30]);
      acc0.vi = _mm_setzero_si128();
      acc1.vi = _mm_setzero_si128();
      acc2.vi = _mm_setzero_si128();
      do {
      two_locus_3x3_tablev_outer:
	loader1 = *vec1++;
	loader20 = *vec20++;
	loader21 = *vec21++;
	loader22 = *vec22++;
	count10 = _mm_and_si128(loader1, loader20);
	count11 = _mm_and_si128(loader1, loader21);
	count12 = _mm_and_si128(loader1, loader22);
	count10 = _mm_sub_epi64(count10, _mm_and_si128(_mm_srli_epi64(count10, 1), m1));
	count11 = _mm_sub_epi64(count11, _mm_and_si128(_mm_srli_epi64(count11, 1), m1));
	count12 = _mm_sub_epi64(count12, _mm_and_si128(_mm_srli_epi64(count12, 1), m1));
      two_locus_3x3_tablev_two_left:
        // unlike the zmiss variant, this apparently does not suffer from
	// enough register spill to justify shrinking the inner loop
	loader1 = *vec1++;
	loader20 = *vec20++;
	loader21 = *vec21++;
	loader22 = *vec22++;
	count20 = _mm_and_si128(loader1, loader20);
	count21 = _mm_and_si128(loader1, loader21);
	count22 = _mm_and_si128(loader1, loader22);
	count20 = _mm_sub_epi64(count20, _mm_and_si128(_mm_srli_epi64(count20, 1), m1));
	count21 = _mm_sub_epi64(count21, _mm_and_si128(_mm_srli_epi64(count21, 1), m1));
	count22 = _mm_sub_epi64(count22, _mm_and_si128(_mm_srli_epi64(count22, 1), m1));
      two_locus_3x3_tablev_one_left:
	loader1 = *vec1++;
	loader20 = *vec20++;
	loader21 = _mm_and_si128(loader1, loader20); // half1
	loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1); // half2
	count10 = _mm_add_epi64(count10, _mm_and_si128(loader21, m1));
	count20 = _mm_add_epi64(count20, loader22);
	loader20 = *vec21++;
	loader21 = _mm_and_si128(loader1, loader20);
	loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
	count11 = _mm_add_epi64(count11, _mm_and_si128(loader21, m1));
	count21 = _mm_add_epi64(count21, loader22);
	loader20 = *vec22++;
	loader21 = _mm_and_si128(loader1, loader20);
	loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
	count12 = _mm_add_epi64(count12, _mm_and_si128(loader21, m1));
	count22 = _mm_add_epi64(count22, loader22);

	count10 = _mm_add_epi64(_mm_and_si128(count10, m2), _mm_and_si128(_mm_srli_epi64(count10, 2), m2));
	count11 = _mm_add_epi64(_mm_and_si128(count11, m2), _mm_and_si128(_mm_srli_epi64(count11, 2), m2));
	count12 = _mm_add_epi64(_mm_and_si128(count12, m2), _mm_and_si128(_mm_srli_epi64(count12, 2), m2));
	count10 = _mm_add_epi64(count10, _mm_add_epi64(_mm_and_si128(count20, m2), _mm_and_si128(_mm_srli_epi64(count20, 2), m2)));
	count11 = _mm_add_epi64(count11, _mm_add_epi64(_mm_and_si128(count21, m2), _mm_and_si128(_mm_srli_epi64(count21, 2), m2)));
	count12 = _mm_add_epi64(count12, _mm_add_epi64(_mm_and_si128(count22, m2), _mm_and_si128(_mm_srli_epi64(count22, 2), m2)));
	acc0.vi = _mm_add_epi64(acc0.vi, _mm_add_epi64(_mm_and_si128(count10, m4), _mm_and_si128(_mm_srli_epi64(count10, 4), m4)));
	acc1.vi = _mm_add_epi64(acc1.vi, _mm_add_epi64(_mm_and_si128(count11, m4), _mm_and_si128(_mm_srli_epi64(count11, 4), m4)));
	acc2.vi = _mm_add_epi64(acc2.vi, _mm_add_epi64(_mm_and_si128(count12, m4), _mm_and_si128(_mm_srli_epi64(count12, 4), m4)));
      } while (vec1 < vend1);
      const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
      acc0.vi = _mm_add_epi64(_mm_and_si128(acc0.vi, m8), _mm_and_si128(_mm_srli_epi64(acc0.vi, 8), m8));
      acc1.vi = _mm_add_epi64(_mm_and_si128(acc1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
      acc2.vi = _mm_add_epi64(_mm_and_si128(acc2.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
      counts_3x3[0] += ((acc0.u8[0] + acc0.u8[1]) * 0x1000100010001LLU) >> 48;
      counts_3x3[1] += ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
      counts_3x3[2] += ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
    }
    if (ct) {
      vend1 = &(vec1[ct]);
      ct2 = ct % 3;
      acc0.vi = _mm_setzero_si128();
      acc1.vi = _mm_setzero_si128();
      acc2.vi = _mm_setzero_si128();
      ct = 0;
      if (ct2) {
	count10 = _mm_setzero_si128();
	count11 = _mm_setzero_si128();
	count12 = _mm_setzero_si128();
	if (ct2 == 2) {
	  goto two_locus_3x3_tablev_two_left;
	}
	count20 = _mm_setzero_si128();
	count21 = _mm_setzero_si128();
	count22 = _mm_setzero_si128();
	goto two_locus_3x3_tablev_one_left;
      }
      goto two_locus_3x3_tablev_outer;
    }
    counts_3x3 = &(counts_3x3[3]);
  }
}

static inline void two_locus_3x3_zmiss_tablev(__m128i* veca0, __m128i* vecb0, uint32_t* counts_3x3, uint32_t sample_ctv6) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i* vecb1 = &(vecb0[sample_ctv6]);
  __m128i* veca1 = &(veca0[sample_ctv6]);
  __m128i* vend;
  __m128i loadera0;
  __m128i loaderb0;
  __m128i loaderb1;
  __m128i loadera1;
  __m128i countx00;
  __m128i countx01;
  __m128i countx11;
  __m128i countx10;
  __m128i county00;
  __m128i county01;
  __m128i county11;
  __m128i county10;
  __univec acc00;
  __univec acc01;
  __univec acc11;
  __univec acc10;
  uint32_t ct2;
  while (sample_ctv6 >= 30) {
    sample_ctv6 -= 30;
    vend = &(veca0[30]);
    acc00.vi = _mm_setzero_si128();
    acc01.vi = _mm_setzero_si128();
    acc11.vi = _mm_setzero_si128();
    acc10.vi = _mm_setzero_si128();
    do {
    two_locus_3x3_zmiss_tablev_outer:
      loadera0 = *veca0++;
      loaderb0 = *vecb0++;
      loaderb1 = *vecb1++;
      loadera1 = *veca1++;
      countx00 = _mm_and_si128(loadera0, loaderb0);
      countx01 = _mm_and_si128(loadera0, loaderb1);
      countx11 = _mm_and_si128(loadera1, loaderb1);
      countx10 = _mm_and_si128(loadera1, loaderb0);
      countx00 = _mm_sub_epi64(countx00, _mm_and_si128(_mm_srli_epi64(countx00, 1), m1));
      countx01 = _mm_sub_epi64(countx01, _mm_and_si128(_mm_srli_epi64(countx01, 1), m1));
      countx11 = _mm_sub_epi64(countx11, _mm_and_si128(_mm_srli_epi64(countx11, 1), m1));
      countx10 = _mm_sub_epi64(countx10, _mm_and_si128(_mm_srli_epi64(countx10, 1), m1));
      countx00 = _mm_add_epi64(_mm_and_si128(countx00, m2), _mm_and_si128(_mm_srli_epi64(countx00, 2), m2));
      countx01 = _mm_add_epi64(_mm_and_si128(countx01, m2), _mm_and_si128(_mm_srli_epi64(countx01, 2), m2));
      countx11 = _mm_add_epi64(_mm_and_si128(countx11, m2), _mm_and_si128(_mm_srli_epi64(countx11, 2), m2));
      countx10 = _mm_add_epi64(_mm_and_si128(countx10, m2), _mm_and_si128(_mm_srli_epi64(countx10, 2), m2));
    two_locus_3x3_zmiss_tablev_one_left:
      loadera0 = *veca0++;
      loaderb0 = *vecb0++;
      loaderb1 = *vecb1++;
      loadera1 = *veca1++;
      county00 = _mm_and_si128(loadera0, loaderb0);
      county01 = _mm_and_si128(loadera0, loaderb1);
      county11 = _mm_and_si128(loadera1, loaderb1);
      county10 = _mm_and_si128(loadera1, loaderb0);
      county00 = _mm_sub_epi64(county00, _mm_and_si128(_mm_srli_epi64(county00, 1), m1));
      county01 = _mm_sub_epi64(county01, _mm_and_si128(_mm_srli_epi64(county01, 1), m1));
      county11 = _mm_sub_epi64(county11, _mm_and_si128(_mm_srli_epi64(county11, 1), m1));
      county10 = _mm_sub_epi64(county10, _mm_and_si128(_mm_srli_epi64(county10, 1), m1));
      countx00 = _mm_add_epi64(countx00, _mm_add_epi64(_mm_and_si128(county00, m2), _mm_and_si128(_mm_srli_epi64(county00, 2), m2)));
      countx01 = _mm_add_epi64(countx01, _mm_add_epi64(_mm_and_si128(county01, m2), _mm_and_si128(_mm_srli_epi64(county01, 2), m2)));
      countx11 = _mm_add_epi64(countx11, _mm_add_epi64(_mm_and_si128(county11, m2), _mm_and_si128(_mm_srli_epi64(county11, 2), m2)));
      countx10 = _mm_add_epi64(countx10, _mm_add_epi64(_mm_and_si128(county10, m2), _mm_and_si128(_mm_srli_epi64(county10, 2), m2)));
      acc00.vi = _mm_add_epi64(acc00.vi, _mm_add_epi64(_mm_and_si128(countx00, m4), _mm_and_si128(_mm_srli_epi64(countx00, 4), m4)));
      acc01.vi = _mm_add_epi64(acc01.vi, _mm_add_epi64(_mm_and_si128(countx01, m4), _mm_and_si128(_mm_srli_epi64(countx01, 4), m4)));
      acc11.vi = _mm_add_epi64(acc11.vi, _mm_add_epi64(_mm_and_si128(countx11, m4), _mm_and_si128(_mm_srli_epi64(countx11, 4), m4)));
      acc10.vi = _mm_add_epi64(acc10.vi, _mm_add_epi64(_mm_and_si128(countx10, m4), _mm_and_si128(_mm_srli_epi64(countx10, 4), m4)));
    } while (veca0 < vend);
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    acc00.vi = _mm_add_epi64(_mm_and_si128(acc00.vi, m8), _mm_and_si128(_mm_srli_epi64(acc00.vi, 8), m8));
    acc01.vi = _mm_add_epi64(_mm_and_si128(acc01.vi, m8), _mm_and_si128(_mm_srli_epi64(acc01.vi, 8), m8));
    acc11.vi = _mm_add_epi64(_mm_and_si128(acc11.vi, m8), _mm_and_si128(_mm_srli_epi64(acc11.vi, 8), m8));
    acc10.vi = _mm_add_epi64(_mm_and_si128(acc10.vi, m8), _mm_and_si128(_mm_srli_epi64(acc10.vi, 8), m8));
    counts_3x3[0] += ((acc00.u8[0] + acc00.u8[1]) * 0x1000100010001LLU) >> 48;
    counts_3x3[1] += ((acc01.u8[0] + acc01.u8[1]) * 0x1000100010001LLU) >> 48;
    counts_3x3[4] += ((acc11.u8[0] + acc11.u8[1]) * 0x1000100010001LLU) >> 48;
    counts_3x3[3] += ((acc10.u8[0] + acc10.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (sample_ctv6) {
    vend = &(veca0[sample_ctv6]);
    ct2 = sample_ctv6 % 2;
    sample_ctv6 = 0;
    acc00.vi = _mm_setzero_si128();
    acc01.vi = _mm_setzero_si128();
    acc11.vi = _mm_setzero_si128();
    acc10.vi = _mm_setzero_si128();
    if (ct2) {
      countx00 = _mm_setzero_si128();
      countx01 = _mm_setzero_si128();
      countx11 = _mm_setzero_si128();
      countx10 = _mm_setzero_si128();
      goto two_locus_3x3_zmiss_tablev_one_left;
    }
    goto two_locus_3x3_zmiss_tablev_outer;
  }
}
#endif

static void two_locus_count_table_zmiss1(uintptr_t* lptr1, uintptr_t* lptr2, uint32_t* counts_3x3, uint32_t sample_ctv3, uint32_t is_zmiss2) {
#ifdef __LP64__
  fill_uint_zero(6, counts_3x3);
  if (is_zmiss2) {
    two_locus_3x3_zmiss_tablev((__m128i*)lptr1, (__m128i*)lptr2, counts_3x3, sample_ctv3 / 2);
  } else {
    two_locus_3x3_tablev((__m128i*)lptr1, (__m128i*)lptr2, counts_3x3, sample_ctv3 / 2, 2);
  }
#else
  counts_3x3[0] = popcount_longs_intersect(lptr1, lptr2, sample_ctv3);
  counts_3x3[1] = popcount_longs_intersect(lptr1, &(lptr2[sample_ctv3]), sample_ctv3);
  if (!is_zmiss2) {
    counts_3x3[2] = popcount_longs_intersect(lptr1, &(lptr2[2 * sample_ctv3]), sample_ctv3);
    counts_3x3[5] = popcount_longs_intersect(&(lptr1[sample_ctv3]), &(lptr2[2 * sample_ctv3]), sample_ctv3);
  }
  lptr1 = &(lptr1[sample_ctv3]);
  counts_3x3[3] = popcount_longs_intersect(lptr1, lptr2, sample_ctv3);
  counts_3x3[4] = popcount_longs_intersect(lptr1, &(lptr2[sample_ctv3]), sample_ctv3);
#endif
}

static void two_locus_count_table(uintptr_t* lptr1, uintptr_t* lptr2, uint32_t* counts_3x3, uint32_t sample_ctv3, uint32_t is_zmiss2) {
#ifdef __LP64__
  uint32_t uii;
  fill_uint_zero(9, counts_3x3);
  if (!is_zmiss2) {
    two_locus_3x3_tablev((__m128i*)lptr1, (__m128i*)lptr2, counts_3x3, sample_ctv3 / 2, 3);
  } else {
    two_locus_3x3_tablev((__m128i*)lptr2, (__m128i*)lptr1, counts_3x3, sample_ctv3 / 2, 2);
    uii = counts_3x3[1];
    counts_3x3[1] = counts_3x3[3];
    counts_3x3[3] = uii;
    counts_3x3[6] = counts_3x3[2];
    counts_3x3[7] = counts_3x3[5];
  }
#else
  counts_3x3[0] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
  counts_3x3[3] = popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
  counts_3x3[6] = popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
  lptr2 = &(lptr2[sample_ctv3]);
  counts_3x3[1] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
  counts_3x3[4] = popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
  counts_3x3[7] = popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
  if (!is_zmiss2) {
    lptr2 = &(lptr2[sample_ctv3]);
    counts_3x3[2] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
    counts_3x3[5] = popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
    counts_3x3[8] = popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
  }
#endif
}

void fepi_counts_to_joint_effects_stats(uint32_t group_ct, uint32_t* counts, double* diff_ptr, double* case_var_ptr, double* ctrl_var_ptr) {
  // See JointEffects::evaluateStatistic().  This is slightly reordered to
  // avoid a bit of redundant calculation, but the logic is otherwise
  // identical.
  //
  // Two adjustments to the raw counts are applied:
  // 1. If any cell in either the case or control tables is zero, we add 0.5 to
  //    all cells in both tables.
  // 2. Then, if the [hom A2 x hom B2] cell in either the case or control table
  //    is less than 1% of the total (very unlikely since A2/B2 are normally
  //    major), multiply all other cells by a reduction factor and increase the
  //    [hom A2 x hom B2] cell by the total reduction (choosing the factor such
  //    that the [hom A2 x hom B2] cell ends up at about 1%).
  //
  // Then, we define
  //   i22_case := [hom A1 x hom B1] * [hom A2 x hom B2] /
  //               ([hom A1 x hom B2] * [hom A2 x hom B1])
  //   i21_case := [hom A1 x het] * [hom A2 x hom B2] /
  //               ([hom A1 x hom B2] * [hom A2 x het])
  //   i12_case := [het x hom B1] * [hom A2 x hom B2] /
  //               ([het x hom B2] * [hom A2 x hom B1])
  //   i11_case := [het x het] * [hom A2 x hom B2] /
  //               ([het x hom B2] * [hom A2 x het])
  //   (analogously for controls)
  //
  // At this point, two formulas may be applied to the (adjusted) counts:
  // 1. If i11 is greater than 0.5 for both cases and controls (this is usually
   //    true),
  //      xi0 := 0.5
  //      xi1 := 1.0
  //      xi2 := 1.0
  //      xi3 := 2 * i11_case / (2 * i11_case - 1)
  //      invq00 := 1.0 / [hom A2 x hom B2]
  //      invq01 := 1.0 / [hom A2 x het]
  //      ...
  //      inverse_matrix := [ (invq22+invq02+invq20+invq00)*xi0*xi0   (invq20+invq00)*xi0*xi1   (invq02+invq00)*xi0*xi2   invq00*xi0*xi3 ]^-1
  //                        [ ...   (invq21+invq20+invq01+invq00)*xi1*xi1   invq00*xi1*xi2   (invq01+invq00)*xi1*xi3 ]
  //                        [ ...   ...   (invq12+invq10+invq02+invq00)*xi2*xi2   (invq10+invq00)*xi2*xi3 ]
  //                        [ ...   ...   ... (invq11+invq10+invq01+invq00)*xi3*xi3 ]
  //      (bottom left is symmetric copy of upper right)
  //      row_totals_case[i] := sum(row i of inverse_matrix_case)
  //      total_inv_v_case := 1.0 / (row_totals_case[0] + [1] + [2] + [3])
  //      lambda_case := row_totals_case[0] * log(i22_case) * 0.5 +
  //                     row_totals_case[1] * log(i21_case) +
  //                     row_totals_case[2] * log(i12_case) +
  //                     row_totals_case[3] * log(2 * i11_case - 1)
  //      (analogous formulas for lambda_ctrl)
  //      diff := lambda_case * total_inv_v_case -
  //              lambda_ctrl * total_inv_v_ctrl
  //      chisq := diff * diff / (total_inv_v_case + total_inv_v_ctrl)
  //
  // 2. Otherwise,
  //      xi0 := sqrt(i22) / (2 * sqrt(i22) + 2)
  //      xi1 := i21 / (i21 + 1)
  //      xi2 := i12 / (i12 + 1)
  //      xi3 := 1.0
  //      (inverse_matrix, row_totals, total_inv_v defined as before)
  //      mu_case := row_totals_case[0] * log((sqrt(i22_case) + 1) * 0.5) +
  //                 row_totals_case[1] * log((i21_case + 1) * 0.5) +
  //                 row_totals_case[2] * log((i12_case + 1) * 0.5) +
  //                 row_totals_case[3] * log(i11_case)
  //      (similar for mu_ctrl)
  //      diff := mu_case * total_inv_v_case - mu_ctrl * total_inv_v_ctrl
  double dcounts[18];
  double invcounts[18];
  double ivv[8]; // i22_case in [0], i21_case in [1], ..., i22_ctrl in [4]...
  double xiv[8];
  double row_totals[8];
  double to_invert[16];
  MATRIX_INVERT_BUF1_TYPE int_1d_buf[4];
  double dbl_2d_buf[16];
  double tot_inv_v[2];
  double lambda_or_mu[2];
  double dxx;
  double dyy;
  double* dptr;
  double* dptr2;
  double* dptr3;
  uint32_t use_reg_stat;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  tot_inv_v[0] = 0.0;  // gcc7 maybe-uninitialized warning
  dptr = dcounts;
  if (counts[0] && counts[1] && counts[2] && counts[3] && counts[4] && counts[5] && counts[6] && counts[7] && counts[8] && ((group_ct == 1) || (counts[9] && counts[10] && counts[11] && counts[12] && counts[13] && counts[14] && counts[15] && counts[16] && counts[17]))) {
    for (uii = 0; uii < group_ct; uii++) {
      dxx = 0;
      for (ujj = 0; ujj < 9; ujj++) {
	dyy = (double)((int32_t)(*counts++));
	*dptr++ = dyy;
	dxx += dyy;
      }
      if (dyy * 100 < dxx) {
	// This tends to come up with adjacent pairs of markers where MAF
	// "flips" from one side of 0.5 to the other.  Is this really a good
	// way to handle it?
	dyy = dxx / (1.01 * dxx - dyy);
        dptr = &(dptr[-9]);
	for (ujj = 0; ujj < 8; ujj++) {
	  *dptr *= dyy;
	  dptr++;
	}
	*dptr++ = 0.01 * dyy * dxx;
      }
    }
  } else {
    for (uii = 0; uii < group_ct; uii++) {
      dxx = -4.5;
      for (ujj = 0; ujj < 9; ujj++) {
	dyy = 0.5 + (double)((int32_t)(*counts++));
	*dptr++ = dyy;
	dxx += dyy;
      }
      if (dyy * 100 < dxx) {
	dyy = dxx / (1.01 * dxx - dyy + 4.5);
        dptr = &(dptr[-9]);
	for (ujj = 0; ujj < 8; ujj++) {
	  *dptr *= dyy;
	  dptr++;
	}
	*dptr++ = 0.01 * dyy * dxx;
      }
    }
  }
  dptr = dcounts;
  dptr2 = invcounts;
  for (uii = 0; uii < group_ct; uii++) {
    for (ujj = 0; ujj < 9; ujj++) {
      *dptr2++ = 1.0 / (*dptr++);
    }
  }
  dptr2 = ivv;
  uii = 0;
  do {
    dptr = &(dcounts[uii * 9]);
    dptr3 = &(invcounts[uii * 9]);
    dxx = dptr[8];
    *dptr2++ = dxx * dptr[0] * dptr3[2] * dptr3[6];
    *dptr2++ = dxx * dptr[1] * dptr3[2] * dptr3[7];
    *dptr2++ = dxx * dptr[3] * dptr3[5] * dptr3[6];
    *dptr2++ = dxx * dptr[4] * dptr3[5] * dptr3[7];
  } while (++uii < group_ct);
  use_reg_stat = (ivv[3] > 0.5) && ((group_ct == 1) || (ivv[7] > 0.5));
  if (use_reg_stat) {
    dptr2 = xiv;
    for (uii = 0; uii < group_ct; uii++) {
      dxx = 2 * ivv[3 + 4 * uii];
      *dptr2++ = 0.5;
      *dptr2++ = 1.0;
      *dptr2++ = 1.0;
      *dptr2++ = dxx / (dxx - 1);
    }
  } else {
    for (uii = 0; uii < group_ct; uii++) {
      dptr = &(ivv[uii * 4]);
      dptr2 = &(xiv[uii * 4]);
      dxx = sqrt(dptr[0]);
      dptr2[1] = dptr[1] / (dptr[1] + 1);
      dptr2[2] = dptr[2] / (dptr[2] + 1);
      dptr2[3] = 1.0;
      dptr2[0] = dxx / (2 * dxx + 2);
      dptr[0] = dxx; // original i22 is not used from here on
    }
  }
  for (uii = 0; uii < group_ct; uii++) {
    dptr = &(invcounts[uii * 9]);
    dptr2 = &(xiv[uii * 4]);
    // invq00 = dptr[8]
    // invq01 = dptr[7]
    // ...
    // thank god this code doesn't need to be edited every day
    dxx = dptr[8];
    dyy = dptr2[0];
    to_invert[0] = (dptr[0] + dptr[2] + dptr[6] + dxx) * dyy * dyy;
    to_invert[1] = (dptr[2] + dxx) * dyy * dptr2[1];
    to_invert[2] = (dptr[6] + dxx) * dyy * dptr2[2];
    to_invert[3] = dxx * dyy * dptr2[3];
    dyy = dptr2[1];
    to_invert[4] = to_invert[1];
    to_invert[5] = (dptr[1] + dptr[2] + dptr[7] + dxx) * dyy * dyy;
    to_invert[6] = dxx * dyy * dptr2[2];
    to_invert[7] = (dptr[7] + dxx) * dyy * dptr2[3];
    dyy = dptr2[2];
    to_invert[8] = to_invert[2];
    to_invert[9] = to_invert[6];
    to_invert[10] = (dptr[3] + dptr[5] + dptr[6] + dxx) * dyy * dyy;
    to_invert[11] = (dptr[5] + dxx) * dyy * dptr2[3];
    dyy = dptr2[3];
    to_invert[12] = to_invert[3];
    to_invert[13] = to_invert[7];
    to_invert[14] = to_invert[11];
    to_invert[15] = (dptr[4] + dptr[5] + dptr[7] + dxx) * dyy * dyy;
    invert_matrix(4, to_invert, int_1d_buf, dbl_2d_buf);
    dptr = to_invert;
    dptr2 = &(row_totals[uii * 4]);
    dxx = 0;
    for (ujj = 0; ujj < 4; ujj++) {
      dyy = 0;
      for (ukk = 0; ukk < 4; ukk++) {
	dyy += (*dptr++);
      }
      *dptr2++ = dyy;
      dxx += dyy;
    }
    tot_inv_v[uii] = dxx;
  }
  if (use_reg_stat) {
    for (uii = 0; uii < group_ct; uii++) {
      dptr = &(row_totals[uii * 4]);
      dptr2 = &(ivv[uii * 4]);
      lambda_or_mu[uii] = dptr[0] * log(dptr2[0]) * 0.5 +
	                  dptr[1] * log(dptr2[1]) +
	                  dptr[2] * log(dptr2[2]) +
                          dptr[3] * log(2 * dptr2[3] - 1);
    }
  } else {
    for (uii = 0; uii < group_ct; uii++) {
      dptr = &(row_totals[uii * 4]);
      dptr2 = &(ivv[uii * 4]);
      // note that dptr2[0] has sqrt(i22) instead of i22
      // really minor thing to check: cheaper to subtract log(2) than multiply
      // by 0.5 inside log?  (I wouldn't think so: multiplication-by-0.5 is the
      // sort of thing which looks like it's eligible for automatic
      // optimization.)
      lambda_or_mu[uii] = dptr[0] * log((dptr2[0] + 1) * 0.5) +
	                  dptr[1] * log((dptr2[1] + 1) * 0.5) +
	                  dptr[2] * log((dptr2[2] + 1) * 0.5) +
	                  dptr[3] * log(dptr2[3]);
    }
  }
  dxx = tot_inv_v[0];
  if (group_ct == 1) {
    *case_var_ptr = dxx;
    *diff_ptr = lambda_or_mu[0];
    return;
  }
  dxx = 1.0 / dxx;
  dyy = 1.0 / tot_inv_v[1];
  *diff_ptr = lambda_or_mu[0] * dxx - lambda_or_mu[1] * dyy;
  *case_var_ptr = dxx;
  *ctrl_var_ptr = dyy;
}

// epistasis multithread globals
static uint32_t* g_epi_geno1_offsets;
static double* g_epi_all_chisq;
static uintptr_t* g_epi_geno1;
static uintptr_t* g_epi_zmiss1;
static uint32_t* g_epi_idx1_block_bounds;
static uint32_t* g_epi_idx1_block_bounds16;
static double* g_epi_best_chisq1;
static uint32_t* g_epi_best_id1; // best partner ID
static uint32_t* g_epi_n_sig_ct1;
static uint32_t* g_epi_fail_ct1;
static uintptr_t* g_epi_geno2;
static uintptr_t* g_epi_zmiss2;
static uint32_t* g_epi_tot2;
static double* g_epi_boost_precalc2 = nullptr;
static double* g_epi_best_chisq2;
static uint32_t* g_epi_best_id2;
static uint32_t* g_epi_n_sig_ct2;
static uint32_t* g_epi_fail_ct2;
static double* g_epi_recip_cache;
static uint32_t g_epi_thread_ct;
static uint32_t g_epi_case_ct;
static uint32_t g_epi_ctrl_ct;
static uint32_t g_epi_flag;
static uint32_t g_epi_cellmin;
static uintptr_t g_epi_marker_ct;
static uintptr_t g_epi_marker_idx1;
static uintptr_t g_epi_idx2_block_size;
static uintptr_t g_epi_idx2_block_start;
static double g_epi_alpha1sq[3];
static double g_epi_alpha2sq[3];

// The following two functions are essentially ported from Statistics.cpp in
// Richard Howey's CASSI software
// (http://www.staff.ncl.ac.uk/richard.howey/cassi/index.html).  (CASSI is also
// GPLv3-licensed; just remember to give credit to Howey if you redistribute a
// variant of this code.  This would have been a friggin' nightmare to debug if
// he hadn't already done all the real work.)
static void fepi_counts_to_stats(uint32_t* counts_3x3, uint32_t no_ueki, double* or_ptr, double* var_ptr) {
  double c11 = (double)((int32_t)(4 * counts_3x3[0] + 2 * (counts_3x3[1] + counts_3x3[3]) + counts_3x3[4]));
  double c12 = (double)((int32_t)(4 * counts_3x3[2] + 2 * (counts_3x3[1] + counts_3x3[5]) + counts_3x3[4]));
  double c21 = (double)((int32_t)(4 * counts_3x3[6] + 2 * (counts_3x3[3] + counts_3x3[7]) + counts_3x3[4]));
  double c22 = (double)((int32_t)(4 * counts_3x3[8] + 2 * (counts_3x3[5] + counts_3x3[7]) + counts_3x3[4]));
  double rc11;
  double rc12;
  double rc21;
  double rc22;
  double dxx;
  uint32_t no_adj;
  if (!no_ueki) {
    // See AdjustedFastEpistasis::calculateLogOddsAdjustedVariance().
    no_adj = (counts_3x3[0] && counts_3x3[1] && counts_3x3[2] && counts_3x3[3] && counts_3x3[4] && counts_3x3[5] && counts_3x3[6] && counts_3x3[7] && counts_3x3[8]);
    if (!no_adj) {
      c11 += 4.5;
      c12 += 4.5;
      c21 += 4.5;
      c22 += 4.5;
    }
    rc11 = 1.0 / c11;
    rc12 = 1.0 / c12;
    rc21 = 1.0 / c21;
    rc22 = 1.0 / c22;
    *or_ptr = log(c11 * c22 * rc12 * rc21);

    c11 = rc11 - rc12; // bit2
    c12 = rc11 - rc21; // bit3
    dxx = rc11 - rc12 - rc21 + rc22; // bit5
    c21 = rc22 - rc12; // bit6
    c22 = rc22 - rc21; // bit8

    rc11 *= rc11;
    rc12 *= rc12;
    rc21 *= rc21;
    rc22 *= rc22;
    c11 *= c11;
    c12 *= c12;
    c21 *= c21;
    c22 *= c22;
    dxx *= dxx;

    if (no_adj) {
      *var_ptr = 4 * (4 * (rc11 * (double)((int32_t)counts_3x3[0]) +
			   rc12 * (double)((int32_t)counts_3x3[2]) +
			   rc21 * (double)((int32_t)counts_3x3[6]) +
			   rc22 * (double)((int32_t)counts_3x3[8])) +
		      c11 * (double)((int32_t)counts_3x3[1]) +
		      c12 * (double)((int32_t)counts_3x3[3]) +
		      c21 * (double)((int32_t)counts_3x3[5]) +
		      c22 * (double)((int32_t)counts_3x3[7])) +
                 dxx * (double)((int32_t)counts_3x3[4]);
    } else {
      *var_ptr = 4 * (4 * (rc11 * ((double)((int32_t)counts_3x3[0]) + 0.5) +
			   rc12 * ((double)((int32_t)counts_3x3[2]) + 0.5) +
			   rc21 * ((double)((int32_t)counts_3x3[6]) + 0.5) +
			   rc22 * ((double)((int32_t)counts_3x3[8]) + 0.5)) +
		      c11 * ((double)((int32_t)counts_3x3[1]) + 0.5) +
		      c12 * ((double)((int32_t)counts_3x3[3]) + 0.5) +
		      c21 * ((double)((int32_t)counts_3x3[5]) + 0.5) +
		      c22 * ((double)((int32_t)counts_3x3[7]) + 0.5)) +
                 dxx * ((double)((int32_t)counts_3x3[4]) + 0.5);
    }
  } else {
    rc11 = 1.0 / c11;
    rc12 = 1.0 / c12;
    rc21 = 1.0 / c21;
    rc22 = 1.0 / c22;
    *or_ptr = log(c11 * c22 * rc12 * rc21);
    *var_ptr = rc11 + rc12 + rc21 + rc22;
  }
}

void boost_calc_p_bc(uint32_t case0_ct, uint32_t case1_ct, uint32_t case2_ct, uint32_t ctrl0_ct, uint32_t ctrl1_ct, uint32_t ctrl2_ct, double* p_bc) {
  double* recip_cache = g_epi_recip_cache;
  double tot_recip = recip_cache[case0_ct + case1_ct + case2_ct];
  p_bc[0] = ((int32_t)case0_ct) * tot_recip;
  p_bc[1] = ((int32_t)case1_ct) * tot_recip;
  p_bc[2] = ((int32_t)case2_ct) * tot_recip;
  tot_recip = recip_cache[ctrl0_ct + ctrl1_ct + ctrl2_ct];
  p_bc[3] = ((int32_t)ctrl0_ct) * tot_recip;
  p_bc[4] = ((int32_t)ctrl1_ct) * tot_recip;
  p_bc[5] = ((int32_t)ctrl2_ct) * tot_recip;
}

uint32_t boost_calc_p_ca(uint32_t case0_ct, uint32_t case1_ct, uint32_t case2_ct, uint32_t ctrl0_ct, uint32_t ctrl1_ct, uint32_t ctrl2_ct, double* p_ca, uint32_t* df_adj_ptr) {
  double* recip_cache = g_epi_recip_cache;
  uint32_t uii = case0_ct + ctrl0_ct;
  uint32_t df_adj = 0;
  double tot_recip;
  tot_recip = recip_cache[uii];
  if (!uii) {
    df_adj++;
  }
  p_ca[0] = ((int32_t)case0_ct) * tot_recip;
  p_ca[1] = ((int32_t)ctrl0_ct) * tot_recip;
  uii = case1_ct + ctrl1_ct;
  tot_recip = recip_cache[uii];
  if (!uii) {
    df_adj++;
  }
  p_ca[2] = ((int32_t)case1_ct) * tot_recip;
  p_ca[3] = ((int32_t)ctrl1_ct) * tot_recip;
  uii = case2_ct + ctrl2_ct;
  tot_recip = recip_cache[uii];
  if (!uii) {
    df_adj++;
  }
  p_ca[4] = ((int32_t)case2_ct) * tot_recip;
  p_ca[5] = ((int32_t)ctrl2_ct) * tot_recip;
  *df_adj_ptr = df_adj;
  return (df_adj > 1);
}

double fepi_counts_to_boost_chisq(uint32_t* counts, double* p_bc, double* p_ca, double* alpha1sq_ptr, double* alpha2sq_ptr, uintptr_t df_adj, double* chisq_ptr, uint32_t* sig_ct1_ptr, uint32_t* sig_ct2_ptr) {
  // see BOOSTx64.c lines 625-903.
  double interaction_measure = 0.0;
  double tau = 0.0;
  double* recip_cache = g_epi_recip_cache;
  uint32_t* uiptr = counts;
  uint32_t sum = 0;
  uint32_t uoo = 0;
  double mu_xx[9]; // initially p_ab
  double mu_tmp[18];
  double mu0_tmp[18];
  double* dptr = mu_xx;
  double sum_recip;
  double dxx;
  double dyy;
  double mu_error;

  // dirty hack: encode df adjustment in low bits of *chisq_ptr
  uintptr_t ularr[sizeof(double) / BYTECT];

  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  for (uii = 0; uii < 3; uii++) {
    ujj = counts[uii] + counts[uii + 9];
    ukk = counts[uii + 3] + counts[uii + 12];
    umm = counts[uii + 6] + counts[uii + 15];
    unn = ujj + ukk + umm;
    if (!unn) {
      if (uoo++) {
	return NAN;
      }
      df_adj++;
    }
    sum += unn;
    dxx = recip_cache[unn];
    *dptr++ = ((int32_t)ujj) * dxx;
    *dptr++ = ((int32_t)ukk) * dxx;
    *dptr++ = ((int32_t)umm) * dxx;
  }
  for (ukk = 0; ukk < 2; ukk++) {
    for (uii = 0; uii < 3; uii++) {
      dyy = p_ca[2 * uii + ukk];
      dptr = &(p_bc[3 * ukk]);
      dxx = mu_xx[uii] * (*dptr++) * dyy;
      tau += dxx;
      umm = *uiptr++;
      if (umm) {
	if (dxx != 0.0) {
	  //   Cx * log(Cx / y)
	  // = Cx * (log(C) + log(x / y))
	  // = Cx * log(C) + Cx * log(x / y)

	  // caching entropy as well would merely reduce a multiplication to
	  // an addition, which is almost certainly not worth the cost
	  interaction_measure -= ((int32_t)umm) * log(dxx * recip_cache[umm]);
	} else {
	  dxx = (double)((int32_t)umm);
	  interaction_measure += dxx * log(dxx);
	}
      }
      dxx = mu_xx[uii + 3] * (*dptr++) * dyy;
      tau += dxx;
      umm = *uiptr++;
      if (umm) {
	if (dxx != 0.0) {
	  interaction_measure -= ((int32_t)umm) * log(dxx * recip_cache[umm]);
	} else {
	  dxx = (double)((int32_t)umm);
	  interaction_measure += dxx * log(dxx);
	}
      }
      dxx = mu_xx[uii + 6] * (*dptr++) * dyy;
      tau += dxx;
      umm = *uiptr++;
      if (umm) {
	if (dxx != 0.0) {
	  interaction_measure -= ((int32_t)umm) * log(dxx * recip_cache[umm]);
	} else {
	  dxx = (double)((int32_t)umm);
	  interaction_measure += dxx * log(dxx);
	}
      }
    }
  }
  // interaction_measure = interaction_measure / sum - log(sum);
  // interaction_measure = (interaction_measure + log(tau)) * sum * 2;
  sum_recip = recip_cache[sum];
  interaction_measure = 2 * (interaction_measure + ((int32_t)sum) * log(tau * sum_recip));
  // > instead of >= for maximum compatibility, I guess
  if (interaction_measure > alpha1sq_ptr[df_adj]) {
    for (uii = 0; uii < 18; uii++) {
      mu_tmp[uii] = 1.0;
    }
    do {
      memcpy(mu0_tmp, mu_tmp, 18 * sizeof(double));
      dptr = mu_xx; // mu_ij
      for (uii = 0; uii < 18; uii += 2) {
        *dptr++ = mu_tmp[uii] + mu_tmp[uii + 1];
      }
      dptr = mu_tmp;
      for (uii = 0; uii < 9; uii++) {
	dxx = mu_xx[uii];
	if (dxx != 0.0) {
	  dxx = (double)((int32_t)(counts[uii] + counts[uii + 9])) / dxx;
	}
	*dptr *= dxx;
	dptr++;
	*dptr *= dxx;
	dptr++;
      }
      dptr = mu_xx; // mu_ik
      for (uii = 0; uii < 18; uii += 6) {
	for (ukk = uii; ukk < uii + 2; ukk++) {
          *dptr++ = mu_tmp[ukk] + mu_tmp[ukk + 2] + mu_tmp[ukk + 4];
	}
      }
      for (uii = 0; uii < 3; uii++) {
	for (ukk = 0; ukk < 2; ukk++) {
	  dxx = mu_xx[uii * 2 + ukk];
          if (dxx != 0.0) {
            dxx = ((double)((int32_t)(counts[ukk * 9 + uii * 3] + counts[ukk * 9 + uii * 3 + 1] + counts[ukk * 9 + uii * 3 + 2]))) / dxx;
	  }
	  mu_tmp[uii * 6 + ukk] *= dxx;
	  mu_tmp[uii * 6 + ukk + 2] *= dxx;
	  mu_tmp[uii * 6 + ukk + 4] *= dxx;
	}
      }
      dptr = mu_xx; // mu_jk
      for (uii = 0; uii < 6; uii++) {
        *dptr = mu_tmp[uii] + mu_tmp[uii + 6] + mu_tmp[uii + 12];
	dptr++;
      }
      for (ujj = 0; ujj < 3; ujj++) {
	for (ukk = 0; ukk < 2; ukk++) {
	  dxx = mu_xx[ujj * 2 + ukk];
          if (dxx != 0.0) {
	    dxx = ((double)((int32_t)(counts[ukk * 9 + ujj] + counts[ukk * 9 + ujj + 3] + counts[ukk * 9 + ujj + 6]))) / dxx;
	  }
          mu_tmp[ujj * 2 + ukk] *= dxx;
          mu_tmp[ujj * 2 + ukk + 6] *= dxx;
          mu_tmp[ujj * 2 + ukk + 12] *= dxx;
	}
      }
      mu_error = 0.0;
      for (uii = 0; uii < 18; uii++) {
        mu_error += fabs(mu_tmp[uii] - mu0_tmp[uii]);
      }
    } while (mu_error > 0.001);
    tau = 0.0;
    interaction_measure = 0.0;
    uiptr = counts;
    for (ukk = 0; ukk < 2; ukk++) {
      for (uii = 0; uii < 3; uii++) {
	for (ujj = 0; ujj < 3; ujj++) {
	  dxx = ((double)((int32_t)(*uiptr++))) * sum_recip;
	  dyy = mu_tmp[uii * 6 + ujj * 2 + ukk] * sum_recip;
	  if (dxx != 0.0) {
	    if (dyy != 0.0) {
	      interaction_measure += dxx * log(dxx / dyy);
	    } else {
              interaction_measure += dxx * log(dxx);
	    }
	  }
	  tau += dyy;
	}
      }
    }
    interaction_measure = (interaction_measure + log(tau)) * ((int32_t)(sum * 2));
    memcpy(ularr, &interaction_measure, sizeof(double));
    // save df_adj in low two bits
    ularr[0] &= ~(3 * ONELU);
    ularr[0] |= df_adj;
    memcpy(chisq_ptr, ularr, sizeof(double));
    if (interaction_measure < alpha1sq_ptr[df_adj]) {
      interaction_measure = alpha1sq_ptr[df_adj];
    }
  }
  if (interaction_measure >= alpha2sq_ptr[df_adj]) {
    *sig_ct1_ptr += 1;
    *sig_ct2_ptr += 1;
  }
  return interaction_measure;
}

THREAD_RET_TYPE fast_epi_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t block_idx1_start = g_epi_idx1_block_bounds[tidx];
  uintptr_t block_idx1_end = g_epi_idx1_block_bounds[tidx + 1];
  uintptr_t idx1_block_start16 = g_epi_idx1_block_bounds16[tidx];
  uintptr_t marker_idx1 = g_epi_marker_idx1 + block_idx1_start;
  uintptr_t marker_ct = g_epi_marker_ct;
  uint32_t case_ct = g_epi_case_ct;
  uint32_t ctrl_ct = g_epi_ctrl_ct;
  uint32_t case_ctv3 = BITCT_TO_ALIGNED_WORDCT(case_ct);
  uint32_t ctrl_ctv3 = BITCT_TO_ALIGNED_WORDCT(ctrl_ct);
  uint32_t case_ctsplit = 3 * case_ctv3;
  uint32_t ctrl_ctsplit = 3 * ctrl_ctv3;
  uint32_t tot_ctsplit = case_ctsplit + ctrl_ctsplit;
  uint32_t is_case_only = (g_epi_flag / EPI_FAST_CASE_ONLY) & 1;
  uint32_t group_ct = 2 - is_case_only;
  uint32_t tot_stride = group_ct * 3;
  uint32_t no_ueki = (g_epi_flag / EPI_FAST_NO_UEKI) & 1;
  uint32_t is_boost = (g_epi_flag / EPI_FAST_BOOST) & 1;
  uint32_t do_joint_effects = (g_epi_flag / EPI_FAST_JOINT_EFFECTS) & 1;
  uint32_t cellmin = g_epi_cellmin;
  uint32_t best_id_fixed = 0;
  uint32_t is_first_half = 0;
  uintptr_t* geno1 = g_epi_geno1;
  uintptr_t* zmiss1 = g_epi_zmiss1;
  uintptr_t* cur_geno1 = nullptr;
  uintptr_t* cur_geno1_ctrls = nullptr;
  double* cur_boost_precalc2 = nullptr;
  double* p_bc_ptr = nullptr;
  uint32_t* geno1_offsets = g_epi_geno1_offsets;
  uint32_t* best_id1 = &(g_epi_best_id1[idx1_block_start16]);
  double* alpha1sq_ptr = g_epi_alpha1sq;
  double* alpha2sq_ptr = g_epi_alpha2sq;
  double alpha1sq = alpha1sq_ptr[0];
  double alpha2sq = alpha2sq_ptr[0];
  double ctrl_var = 0;
  uint32_t tot1[6];
  uint32_t counts[18];
  double p_bc_tmp[6];
  double p_ca_fixed[6];
  double p_ca_tmp[6];
  uintptr_t* geno2;
  uintptr_t* zmiss2;
  uintptr_t* cur_geno2;
  double* all_chisq_write;
  double* chisq2_ptr;
  double* boost_precalc2;
  double* all_chisq;
  double* best_chisq1;
  double* best_chisq2;
  double* p_ca_ptr;
  uint32_t* n_sig_ct1;
  uint32_t* fail_ct1;
  uint32_t* best_id2;
  uint32_t* n_sig_ct2;
  uint32_t* fail_ct2;
  uint32_t* tot2;
  uint32_t* cur_tot2;
  uintptr_t idx2_block_size;
  uintptr_t cur_idx2_block_size;
  uintptr_t idx2_block_start;
  uintptr_t idx2_block_end;
  uintptr_t idx2_block_sizea16;
  uintptr_t block_idx1;
  uintptr_t block_delta1;
  uintptr_t block_idx2;
  uintptr_t cur_zmiss2;
  uintptr_t cur_zmiss2_tmp;
  uintptr_t ulii;
  double best_chisq_fixed;
  double case_var;
  double ctrl_or;
  double dxx;
  double zsq;
  uint32_t nm_case_fixed;
  uint32_t nm_ctrl_fixed;
  uint32_t nm_fixed;
  uint32_t n_sig_ct_fixed;
  uint32_t fail_ct_fixed;
  uint32_t df_adj_base;
  uint32_t df_adj;
  tot1[3] = 0; // suppress warning
  tot1[4] = 0;
  tot1[5] = 0;
  while (1) {
    idx2_block_size = g_epi_idx2_block_size;
    cur_idx2_block_size = idx2_block_size;
    idx2_block_start = g_epi_idx2_block_start;
    idx2_block_end = idx2_block_start + idx2_block_size;
    idx2_block_sizea16 = round_up_pow2(idx2_block_size, 16);
    geno2 = g_epi_geno2;
    zmiss2 = g_epi_zmiss2;
    tot2 = g_epi_tot2;
    boost_precalc2 = g_epi_boost_precalc2;
    all_chisq = &(g_epi_all_chisq[idx2_block_start]);
    best_chisq1 = &(g_epi_best_chisq1[idx1_block_start16]);
    best_chisq2 = &(g_epi_best_chisq2[tidx * idx2_block_sizea16]);
    n_sig_ct1 = &(g_epi_n_sig_ct1[idx1_block_start16]);
    fail_ct1 = &(g_epi_fail_ct1[idx1_block_start16]);
    best_id2 = &(g_epi_best_id2[tidx * idx2_block_sizea16]);
    n_sig_ct2 = &(g_epi_n_sig_ct2[tidx * idx2_block_sizea16]);
    fail_ct2 = &(g_epi_fail_ct2[tidx * idx2_block_sizea16]);
    for (block_idx1 = block_idx1_start; block_idx1 < block_idx1_end; block_idx1++, marker_idx1++) {
      ulii = geno1_offsets[2 * block_idx1];
      if (ulii > idx2_block_start) {
	block_idx2 = 0;
	cur_idx2_block_size = ulii - idx2_block_start;
	if (cur_idx2_block_size >= idx2_block_size) {
	  cur_idx2_block_size = idx2_block_size;
	} else {
	  is_first_half = 1;
	}
      } else {
	ulii = geno1_offsets[2 * block_idx1 + 1];
	if (ulii >= idx2_block_end) {
	  // may not be done in set1 x all or set1 x set2 cases
	  continue;
	} else {
	  if (ulii <= idx2_block_start) {
	    block_idx2 = 0;
	  } else {
	    block_idx2 = ulii - idx2_block_start;
	  }
	}
      }
      cur_geno1 = &(geno1[block_idx1 * tot_ctsplit]);
      n_sig_ct_fixed = 0;
      fail_ct_fixed = 0;
      nm_case_fixed = is_set_ul(zmiss1, block_idx1 * 2);
      nm_ctrl_fixed = is_set_ul(zmiss1, block_idx1 * 2 + 1);
      nm_fixed = nm_case_fixed & nm_ctrl_fixed;
      tot1[0] = popcount_longs(cur_geno1, case_ctv3);
      tot1[1] = popcount_longs(&(cur_geno1[case_ctv3]), case_ctv3);
      tot1[2] = popcount_longs(&(cur_geno1[2 * case_ctv3]), case_ctv3);
      if (!is_case_only) {
	cur_geno1_ctrls = &(cur_geno1[case_ctsplit]);
	tot1[3] = popcount_longs(cur_geno1_ctrls, ctrl_ctv3);
	tot1[4] = popcount_longs(&(cur_geno1_ctrls[ctrl_ctv3]), ctrl_ctv3);
	tot1[5] = popcount_longs(&(cur_geno1_ctrls[2 * ctrl_ctv3]), ctrl_ctv3);
	if (is_boost) {
	  if (nm_fixed) {
	    cur_boost_precalc2 = &(boost_precalc2[block_idx2 * 6]);
	  } else {
	    p_bc_ptr = p_bc_tmp;
	  }
	  boost_calc_p_ca(tot1[0], tot1[1], tot1[2], tot1[3], tot1[4], tot1[5], p_ca_fixed, &df_adj_base);
	}
      }
      block_delta1 = block_idx1 - block_idx1_start;
      best_chisq_fixed = best_chisq1[block_delta1];
      all_chisq_write = &(all_chisq[block_idx1 * marker_ct]);
    fast_epi_thread_second_half:
      cur_geno2 = &(geno2[block_idx2 * tot_ctsplit]);
      chisq2_ptr = &(best_chisq2[block_idx2]);
      for (; block_idx2 < cur_idx2_block_size; block_idx2++, chisq2_ptr++, cur_geno2 = &(cur_geno2[tot_ctsplit])) {
	cur_tot2 = &(tot2[block_idx2 * tot_stride]);
	// this operation isn't extracting a 2-bit genotype, so don't use the
	// macro
	cur_zmiss2 = (zmiss2[block_idx2 / BITCT2] >> (2 * (block_idx2 % BITCT2))) & 3;
	cur_zmiss2_tmp = cur_zmiss2 & 1;
	if (nm_case_fixed) {
	  two_locus_count_table_zmiss1(cur_geno1, cur_geno2, counts, case_ctv3, cur_zmiss2_tmp);
	  if (cur_zmiss2_tmp) {
	    counts[2] = tot1[0] - counts[0] - counts[1];
	    counts[5] = tot1[1] - counts[3] - counts[4];
	  }
	  counts[6] = cur_tot2[0] - counts[0] - counts[3];
	  counts[7] = cur_tot2[1] - counts[1] - counts[4];
	  counts[8] = cur_tot2[2] - counts[2] - counts[5];
	} else {
	  two_locus_count_table(cur_geno1, cur_geno2, counts, case_ctv3, cur_zmiss2_tmp);
	  if (cur_zmiss2_tmp) {
	    counts[2] = tot1[0] - counts[0] - counts[1];
	    counts[5] = tot1[1] - counts[3] - counts[4];
	    counts[8] = tot1[2] - counts[6] - counts[7];
	  }
	}
	if (!is_case_only) {
	  cur_zmiss2_tmp = cur_zmiss2 >> 1;
	  if (nm_ctrl_fixed) {
	    two_locus_count_table_zmiss1(cur_geno1_ctrls, &(cur_geno2[case_ctsplit]), &(counts[9]), ctrl_ctv3, cur_zmiss2_tmp);
	    if (cur_zmiss2_tmp) {
	      counts[11] = tot1[3] - counts[9] - counts[10];
	      counts[14] = tot1[4] - counts[12] - counts[13];
	    }
	    counts[15] = cur_tot2[3] - counts[9] - counts[12];
	    counts[16] = cur_tot2[4] - counts[10] - counts[13];
	    counts[17] = cur_tot2[5] - counts[11] - counts[14];
	  } else {
	    two_locus_count_table(cur_geno1_ctrls, &(cur_geno2[case_ctsplit]), &(counts[9]), ctrl_ctv3, cur_zmiss2_tmp);
	    if (cur_zmiss2_tmp) {
	      counts[11] = tot1[3] - counts[9] - counts[10];
	      counts[14] = tot1[4] - counts[12] - counts[13];
	      counts[17] = tot1[5] - counts[15] - counts[16];
	    }
	  }
	}
	if (!is_boost) {
	  if (!do_joint_effects) {
	    fepi_counts_to_stats(counts, no_ueki, &dxx, &case_var);
	    if (!is_case_only) {
	      fepi_counts_to_stats(&(counts[9]), no_ueki, &ctrl_or, &ctrl_var);
	      dxx -= ctrl_or;
	    }
	  } else {
	    if (cellmin) {
	      if ((counts[0] < cellmin) || (counts[1] < cellmin) || (counts[2] < cellmin) || (counts[3] < cellmin) || (counts[4] < cellmin) || (counts[5] < cellmin) || (counts[6] < cellmin) || (counts[7] < cellmin) || (counts[8] < cellmin)) {
		goto fast_epi_thread_fail;
	      }
	      if (!is_case_only) {
		if ((counts[9] < cellmin) || (counts[10] < cellmin) || (counts[11] < cellmin) || (counts[12] < cellmin) || (counts[13] < cellmin) || (counts[14] < cellmin) || (counts[15] < cellmin) || (counts[16] < cellmin) || (counts[17] < cellmin)) {
		  goto fast_epi_thread_fail;
		}
	      }
	    }
	    fepi_counts_to_joint_effects_stats(group_ct, counts, &dxx, &case_var, &ctrl_var);
	  }
	  zsq = dxx * dxx / (case_var + ctrl_var);
	  if (!realnum(zsq)) {
	    goto fast_epi_thread_fail;
	  }
	  if (zsq >= alpha1sq) {
	    all_chisq_write[block_idx2] = zsq;
	  }
	  if (zsq >= alpha2sq) {
	    n_sig_ct_fixed++;
	    n_sig_ct2[block_idx2] += 1;
	  }
	fast_epi_thread_boost_save:
	  if (zsq > best_chisq_fixed) {
	    best_chisq_fixed = zsq;
	    best_id_fixed = block_idx2 + idx2_block_start;
	  }
	  dxx = *chisq2_ptr;
	  if (zsq > dxx) {
	    *chisq2_ptr = zsq;
	    best_id2[block_idx2] = marker_idx1;
	  }
	} else {
	  if (nm_fixed) {
	    p_bc_ptr = cur_boost_precalc2;
	    cur_boost_precalc2 = &(cur_boost_precalc2[6]);
	  } else {
	    boost_calc_p_bc(counts[0] + counts[3] + counts[6], counts[1] + counts[4] + counts[7], counts[2] + counts[5] + counts[8], counts[9] + counts[12] + counts[15], counts[10] + counts[13] + counts[16], counts[11] + counts[14] + counts[17], p_bc_ptr);
	  }
	  if (cur_zmiss2 == 3) {
	    p_ca_ptr = p_ca_fixed;
	    df_adj = df_adj_base;
	  } else {
	    if (boost_calc_p_ca(counts[0] + counts[1] + counts[2], counts[3] + counts[4] + counts[5], counts[6] + counts[7] + counts[8], counts[9] + counts[10] + counts[11], counts[12] + counts[13] + counts[14], counts[15] + counts[16] + counts[17], p_ca_tmp, &df_adj)) {
	      goto fast_epi_thread_fail;
	    }
	    p_ca_ptr = p_ca_tmp;
	  }

	  // if approximate zsq >= epi1 threshold but more accurate value is
	  // not, we still want to save the more accurate value
	  // also, we want epi2 counting to be df-sensitive
	  // (punt on df/best_chisq for now)
	  zsq = fepi_counts_to_boost_chisq(counts, p_bc_ptr, p_ca_ptr, alpha1sq_ptr, alpha2sq_ptr, df_adj, &(all_chisq_write[block_idx2]), &n_sig_ct_fixed, &(n_sig_ct2[block_idx2]));
	  if (realnum(zsq)) {
	    goto fast_epi_thread_boost_save;
	  }
	fast_epi_thread_fail:
	  fail_ct_fixed++;
	  fail_ct2[block_idx2] += 1;
	  if (alpha1sq == 0.0) {
	    // special case: log NA when '--epi1 1' specified
	    all_chisq_write[block_idx2] = NAN;
	  }
	}
      }
      if (is_first_half) {
	is_first_half = 0;
	ulii = geno1_offsets[2 * block_idx1 + 1];
	cur_idx2_block_size = idx2_block_size;
	if (ulii < idx2_block_end) {
	  // guaranteed to be larger than idx2_block_start, otherwise there
	  // would have been no first half
	  block_idx2 = ulii - idx2_block_start;
	  if (is_boost && nm_fixed) {
	    cur_boost_precalc2 = &(boost_precalc2[block_idx2 * 6]);
	  }
	  goto fast_epi_thread_second_half;
	}
      }
      if (best_chisq_fixed > best_chisq1[block_delta1]) {
	best_chisq1[block_delta1] = best_chisq_fixed;
	best_id1[block_delta1] = best_id_fixed;
      }
      n_sig_ct1[block_delta1] = n_sig_ct_fixed;
      if (fail_ct_fixed) {
	fail_ct1[block_delta1] = fail_ct_fixed;
      }
    }
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

// epistasis linear/logistic regression multithread globals

static double* g_epi_pheno_d2;
static double* g_epi_phenogeno1;
static double* g_epi_phenogeno2;
static uint32_t* g_epi_genosums1;
static uint32_t* g_epi_genosums2;
static double g_epi_pheno_sum;
static double g_epi_pheno_ssq;
static double g_epi_vif_thresh;

static uint32_t g_epi_pheno_nm_ct;

typedef struct epi_logistic_multithread_struct {
  float* cur_covars_cov_major;
  float* coef;
  float* pp;
  float* sample_1d_buf;
  float* pheno_buf;
  float* param_1d_buf;
  float* param_1d_buf2;
  float* param_2d_buf;
  float* param_2d_buf2;
} Epi_logistic_multithread;

static Epi_logistic_multithread* g_epi_logistic_mt;
static uintptr_t* g_epi_pheno_c;
static float* g_epi_all_chisq_f;
static float* g_epi_best_chisq_f1;
static float* g_epi_best_chisq_f2;

uint32_t matrix_invert_4x4symm(double* dmatrix) {
  double buf[16];
  double determinant;
  // initially, dww = A_{22}A_{34} - A_{23}A_{24}
  //            dxx = A_{23}A_{34} - A_{24}A_{33}
  //            dyy = A_{23}A_{44} - A_{24}A_{34}
  //            dzz = A_{33}A_{44} - A_{34}A_{34}
  double dww = dmatrix[5] * dmatrix[11] - dmatrix[6] * dmatrix[7];
  double dxx = dmatrix[6] * dmatrix[11] - dmatrix[7] * dmatrix[10];
  double dyy = dmatrix[6] * dmatrix[15] - dmatrix[7] * dmatrix[11];
  double dzz = dmatrix[10] * dmatrix[15] - dmatrix[11] * dmatrix[11];
  double dvv;
  double duu;
  buf[0] = dmatrix[5] * dzz
         - dmatrix[6] * dyy
         + dmatrix[7] * dxx;
  buf[1] = dmatrix[2] * dyy
         - dmatrix[1] * dzz
         - dmatrix[3] * dxx;
  buf[2] = dmatrix[1] * dyy
         + dmatrix[2] * (dmatrix[7] * dmatrix[7] - dmatrix[5] * dmatrix[15])
         + dmatrix[3] * dww;
  duu = dmatrix[5] * dmatrix[10] - dmatrix[6] * dmatrix[6];
  buf[3] = dmatrix[2] * dww
         - dmatrix[1] * dxx
         - dmatrix[3] * duu;
  determinant = dmatrix[0] * buf[0] + dmatrix[1] * buf[1] + dmatrix[2] * buf[2] + dmatrix[3] * buf[3];
  if (fabs(determinant) < EPSILON) {
    return 1;
  }
  buf[5] = dmatrix[0] * dzz
         + dmatrix[2] * (dmatrix[3] * dmatrix[11] - dmatrix[2] * dmatrix[15])
         + dmatrix[3] * (dmatrix[2] * dmatrix[11] - dmatrix[3] * dmatrix[10]);
  dzz = dmatrix[1] * dmatrix[15] - dmatrix[3] * dmatrix[7];
  buf[6] = dmatrix[2] * dzz
         - dmatrix[0] * dyy
         + dmatrix[3] * (dmatrix[3] * dmatrix[6] - dmatrix[1] * dmatrix[11]);
  dyy = dmatrix[1] * dmatrix[11] - dmatrix[2] * dmatrix[7];
  dvv = dmatrix[1] * dmatrix[10] - dmatrix[2] * dmatrix[6];
  buf[7] = dmatrix[0] * dxx
         - dmatrix[2] * dyy
         + dmatrix[3] * dvv;
  buf[10] = dmatrix[0] * (dmatrix[5] * dmatrix[15] - dmatrix[7] * dmatrix[7])
          - dmatrix[1] * dzz
          + dmatrix[3] * (dmatrix[1] * dmatrix[7] - dmatrix[3] * dmatrix[5]);
  dxx = dmatrix[1] * dmatrix[6] - dmatrix[2] * dmatrix[5];
  buf[11] = dmatrix[1] * dyy
          - dmatrix[0] * dww
          - dmatrix[3] * dxx;
  buf[15] = dmatrix[0] * duu
          - dmatrix[1] * dvv
          + dmatrix[2] * dxx;
  determinant = 1.0 / determinant; // now reciprocal
  dmatrix[0] = buf[0] * determinant;
  dmatrix[1] = buf[1] * determinant;
  dmatrix[2] = buf[2] * determinant;
  dmatrix[3] = buf[3] * determinant;
  dmatrix[4] = dmatrix[1];
  dmatrix[5] = buf[5] * determinant;
  dmatrix[6] = buf[6] * determinant;
  dmatrix[7] = buf[7] * determinant;
  dmatrix[8] = dmatrix[2];
  dmatrix[9] = dmatrix[6];
  dmatrix[10] = buf[10] * determinant;
  dmatrix[11] = buf[11] * determinant;
  dmatrix[12] = dmatrix[3];
  dmatrix[13] = dmatrix[7];
  dmatrix[14] = dmatrix[11];
  dmatrix[15] = buf[15] * determinant;
  return 0;
}

THREAD_RET_TYPE epi_linear_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t block_idx1_start = g_epi_idx1_block_bounds[tidx];
  uintptr_t block_idx1_end = g_epi_idx1_block_bounds[tidx + 1];
  uintptr_t idx1_block_start16 = g_epi_idx1_block_bounds16[tidx];
  uintptr_t marker_idx1 = g_epi_marker_idx1 + block_idx1_start;
  uintptr_t marker_ct = g_epi_marker_ct;
  double alpha1sq = g_epi_alpha1sq[0];
  double alpha2sq = g_epi_alpha2sq[0];
  double pheno_sum = g_epi_pheno_sum;
  double pheno_ssq = g_epi_pheno_ssq;
  double vif_thresh = g_epi_vif_thresh;
  uint32_t pheno_nm_ct = g_epi_pheno_nm_ct;
  uint32_t best_id_fixed = 0;
  uint32_t is_first_half = 0;
  uintptr_t pheno_nm_ctl2 = QUATERCT_TO_WORDCT(pheno_nm_ct);
  uintptr_t* geno1 = g_epi_geno1;
  double* pheno_d2 = g_epi_pheno_d2;
  uint32_t* geno1_offsets = g_epi_geno1_offsets;
  uint32_t* best_id1 = &(g_epi_best_id1[idx1_block_start16]);
  const double dconst[] = {1.0, 2.0, 2.0, 4.0};
  double dmatrix_buf[16];
  double dmatrix_buf2[4];

  // sum(aa), sum(ab), sum(bb), sum(aab), sum(abb), and sum(aabb) can all be
  // derived from these four quantities.
  uint32_t cur_minor_cts[4]; // 11, 12, 21, 22

  uintptr_t* cur_geno1;
  uintptr_t* geno2;
  uintptr_t* cur_geno2;
  double* phenogeno1;
  double* phenogeno2;
  double* all_chisq_write;
  double* chisq2_ptr;
  double* all_chisq;
  double* best_chisq1;
  double* best_chisq2;
  double* dptr;
  double* dptr2;
  uint32_t* n_sig_ct1;
  uint32_t* fail_ct1;
  uint32_t* best_id2;
  uint32_t* n_sig_ct2;
  uint32_t* fail_ct2;
  uint32_t* genosums1;
  uint32_t* genosums2;
  uintptr_t idx2_block_size;
  uintptr_t cur_idx2_block_size;
  uintptr_t idx2_block_start;
  uintptr_t idx2_block_end;
  uintptr_t idx2_block_sizea16;
  uintptr_t block_idx1;
  uintptr_t block_delta1;
  uintptr_t block_idx2;
  uintptr_t cur_word1;
  uintptr_t cur_word2;
  uintptr_t active_mask;
  uintptr_t param_idx;
  uintptr_t param_idx2;
  uintptr_t cur_sum_aab;
  uintptr_t cur_sum_abb;
  uintptr_t cur_sum_aabb;
  uintptr_t ulii;
  uintptr_t uljj;
  double best_chisq_fixed;
  double sum_a_pheno_base;
  double cur_pheno_sum;
  double cur_pheno_ssq;
  double cur_sum_a_pheno;
  double cur_sum_b_pheno;
  double cur_sum_ab_pheno;
  double sample_ctd;
  double sample_ct_recip;
  double sample_ct_m1_recip;
  double cur_sum_ad;
  double cur_sum_bd;
  double cur_sum_abd;
  double determinant;
  double min_sigma;
  double sigma;
  double dxx;
  double dyy;
  double dzz;
  double dww;
  double dvv;
  double duu;
  double zsq;
  uint32_t n_sig_ct_fixed;
  uint32_t fail_ct_fixed;

  uint32_t sum_a_base;
  uint32_t sum_aa_base;
  uint32_t cur_sum_a;
  uint32_t cur_sum_aa;
  uint32_t cur_sum_b;
  uint32_t cur_sum_bb;
  uint32_t cur_sum_ab;
  uint32_t widx;
  uint32_t sample_idx;
  uint32_t cur_sample_ct;
  uint32_t woffset;
  while (1) {
    idx2_block_size = g_epi_idx2_block_size;
    cur_idx2_block_size = idx2_block_size;
    idx2_block_start = g_epi_idx2_block_start;
    idx2_block_end = idx2_block_start + idx2_block_size;
    idx2_block_sizea16 = round_up_pow2(idx2_block_size, 16);
    geno2 = g_epi_geno2;
    phenogeno1 = g_epi_phenogeno1;
    phenogeno2 = g_epi_phenogeno2;
    genosums1 = g_epi_genosums1;
    genosums2 = g_epi_genosums2;
    all_chisq = &(g_epi_all_chisq[2 * idx2_block_start]);
    best_chisq1 = &(g_epi_best_chisq1[idx1_block_start16]);
    best_chisq2 = &(g_epi_best_chisq2[tidx * idx2_block_sizea16]);
    n_sig_ct1 = &(g_epi_n_sig_ct1[idx1_block_start16]);
    fail_ct1 = &(g_epi_fail_ct1[idx1_block_start16]);
    best_id2 = &(g_epi_best_id2[tidx * idx2_block_sizea16]);
    n_sig_ct2 = &(g_epi_n_sig_ct2[tidx * idx2_block_sizea16]);
    fail_ct2 = &(g_epi_fail_ct2[tidx * idx2_block_sizea16]);
    for (block_idx1 = block_idx1_start; block_idx1 < block_idx1_end; block_idx1++, marker_idx1++) {
      ulii = geno1_offsets[2 * block_idx1];
      if (ulii > idx2_block_start) {
        block_idx2 = 0;
        cur_idx2_block_size = ulii - idx2_block_start;
	if (cur_idx2_block_size >= idx2_block_size) {
          cur_idx2_block_size = idx2_block_size;
	} else {
	  is_first_half = 1;
        }
      } else {
        ulii = geno1_offsets[2 * block_idx1 + 1];
        if (ulii >= idx2_block_end) {
          // may not be done in set1 x all or set1 x set2 cases
          continue;
	} else {
          if (ulii <= idx2_block_start) {
            block_idx2 = 0;
	  } else {
            block_idx2 = ulii - idx2_block_start;
	  }
	}
      }
      cur_geno1 = &(geno1[block_idx1 * pheno_nm_ctl2]);
      n_sig_ct_fixed = 0;
      fail_ct_fixed = 0;
      block_delta1 = block_idx1 - block_idx1_start;
      best_chisq_fixed = best_chisq1[block_delta1];
      sum_a_pheno_base = phenogeno1[block_idx1];
      sum_a_base = genosums1[2 * block_idx1];
      sum_aa_base = genosums1[2 * block_idx1 + 1];

      // [0] = chisq, [1] = beta
      all_chisq_write = &(all_chisq[block_idx1 * marker_ct * 2]);

    epi_linear_thread_second_half:
      cur_geno2 = &(geno2[block_idx2 * pheno_nm_ctl2]);
      chisq2_ptr = &(best_chisq2[block_idx2]);
      for (; block_idx2 < cur_idx2_block_size; block_idx2++, chisq2_ptr++, cur_geno2 = &(cur_geno2[pheno_nm_ctl2])) {
	// Our covariates are 1, genotype A (in {0, 1, 2}), genotype B, and
	//   [genotype A] * [genotype B].
	// The ordinary least squares solution to this system is
	//   (X^T X)^{-1} X^T Y
	// where X^T X is the following 4x4 matrix (where n = # of samples):
	//   [ n       sum(A)   sum(B)   sum(AB)   ]
	//   [ sum(A)  sum(AA)  sum(AB)  sum(AAB)  ]
	//   [ sum(B)  sum(AB)  sum(BB)  sum(ABB)  ]
	//   [ sum(AB) sum(AAB) sum(ABB) sum(AABB) ]
        // (sum(.) denotes the sum of that (product of) genotypes, across all
	// samples.)
	// Meanwhile, X^T Y is the following 4x1 matrix:
	//   [ sum(pheno)      ]
	//   [ sum(A * pheno)  ]
	//   [ sum(B * pheno)  ]
	//   [ sum(AB * pheno) ]
	// Crucially, the VIF and valid parameters checks can also operate
	// purely on the terms above and sum(pheno * pheno).

	// these nine values can be mostly precomputed; just need to subtract
	// from them sometimes when missing values are present.
	cur_pheno_sum = pheno_sum;
	cur_pheno_ssq = pheno_ssq;
	cur_sum_a_pheno = sum_a_pheno_base;
	cur_sum_b_pheno = phenogeno2[block_idx2];
	cur_sum_a = sum_a_base;
	cur_sum_aa = sum_aa_base;
	cur_sum_b = genosums2[block_idx2 * 2];
	cur_sum_bb = genosums2[block_idx2 * 2 + 1];
	cur_sample_ct = pheno_nm_ct;

	cur_sum_ab_pheno = 0.0;
	fill_uint_zero(4, cur_minor_cts);
	for (widx = 0; widx < pheno_nm_ctl2; widx++) {
	  sample_idx = widx * BITCT2;
          cur_word1 = cur_geno1[widx];
          cur_word2 = cur_geno2[widx];
	  // we can entirely skip 5 common cases: 00/00, 00/01, 00/10, 01/00,
	  // 10/00.
	  active_mask = cur_word1 | cur_word2;
	  active_mask = (active_mask & (active_mask >> 1) & FIVEMASK) | (cur_word1 & cur_word2);
	  dptr = &(pheno_d2[sample_idx]);
	  while (active_mask) {
            woffset = CTZLU(active_mask) / 2;
	    dxx = dptr[woffset];
	    woffset *= 2;
	    ulii = (cur_word1 >> woffset) & (3 * ONELU);
            uljj = (cur_word2 >> woffset) & (3 * ONELU);
	    active_mask &= ~((3 * ONELU) << woffset);
	    if (ulii && uljj) {
	      if (ulii == 3) {
		if (uljj == 1) {
		  cur_sum_b_pheno -= dxx;
		  cur_sum_b--;
		  cur_sum_bb--;
		} else if (uljj == 2) {
		  cur_sum_b_pheno -= 2 * dxx;
		  cur_sum_b -= 2;
		  cur_sum_bb -= 4;
		}
	      } else if (uljj == 3) {
		// ulii must be 1 or 2
		cur_sum_a_pheno -= dxx;
		if (ulii == 2) {
		  cur_sum_a_pheno -= dxx;
		}
		cur_sum_a -= ulii;
		cur_sum_aa -= ulii * ulii;
	      } else {
		ulii = ulii * 2 + uljj - 3;
		cur_sum_ab_pheno += dconst[ulii] * dxx;
		cur_minor_cts[ulii] += 1;
		continue;
	      }
	    }
	    cur_pheno_sum -= dxx;
	    cur_pheno_ssq -= dxx * dxx;
	    cur_sample_ct--;
	  }
	}
	if (cur_sample_ct <= 4) {
          goto epi_linear_thread_regression_fail;
	}
	// VIF check.  Mirrors glm_check_vif(), but param_ct is hardcoded to 4
	// and we avoid additional iteration over the sample_idxs.
	sample_ctd = (double)((int32_t)cur_sample_ct);
	sample_ct_recip = 1.0 / sample_ctd;
	sample_ct_m1_recip = 1.0 / ((double)((int32_t)(cur_sample_ct - 1)));
	cur_sum_ab = cur_minor_cts[0] + 2 * (cur_minor_cts[1] + cur_minor_cts[2]) + 4 * cur_minor_cts[3];
	cur_sum_aab = cur_minor_cts[0] + 2 * cur_minor_cts[1] + 4 * cur_minor_cts[2] + (8 * ONELU) * cur_minor_cts[3];
	cur_sum_abb = cur_minor_cts[0] + 4 * cur_minor_cts[1] + 2 * cur_minor_cts[2] + (8 * ONELU) * cur_minor_cts[3];
	cur_sum_aabb = cur_minor_cts[0] + 4 * (cur_minor_cts[1] + cur_minor_cts[2]) + (16 * ONELU) * cur_minor_cts[3];

	cur_sum_ad = (double)((int32_t)cur_sum_a);
	cur_sum_bd = (double)((int32_t)cur_sum_b);
	cur_sum_abd = (double)((int32_t)cur_sum_ab);

	// some genotype means
	dxx = cur_sum_bd * sample_ct_recip;
	dyy = cur_sum_abd * sample_ct_recip;

	dww = ((double)((int32_t)cur_sum_aa)) - cur_sum_ad * cur_sum_ad * sample_ct_recip;
	dvv = ((double)((int32_t)cur_sum_bb)) - cur_sum_bd * dxx;
	duu = ((double)((intptr_t)cur_sum_aabb)) - cur_sum_abd * dyy;
	if ((dww <= 0) || (dvv <= 0) || (duu <= 0)) {
	  goto epi_linear_thread_regression_fail;
	}
	dww = 1.0 / sqrt(dww * sample_ct_m1_recip);
	dvv = 1.0 / sqrt(dvv * sample_ct_m1_recip);
	duu = 1.0 / sqrt(duu * sample_ct_m1_recip);

	dxx = (cur_sum_abd - cur_sum_ad * dxx) * sample_ct_m1_recip;
	dzz = (((double)((intptr_t)cur_sum_abb)) - cur_sum_bd * dyy) * sample_ct_m1_recip;
	dyy = (((double)((intptr_t)cur_sum_aab)) - cur_sum_ad * dyy) * sample_ct_m1_recip;
	// now dxx = A_{12}, dyy = A_{13}, dzz = A_{23}

	dxx *= dww * dvv;
	dyy *= dww * duu;
	dzz *= dvv * duu;
	if ((dxx > 0.999) || (dyy > 0.999) || (dzz > 0.999)) {
	  goto epi_linear_thread_regression_fail;
	}
	// Use analytic formula for 3x3 symmetric matrix inverse.
	// det A = A_{11}A_{22}A_{33} + 2 * A_{12}A_{13}A_{23}
	//       - A_{11}(A_{23}^2) - A_{22}(A_{13}^2) - A_{33}(A_{12}^2)
	// upper left of inverse = (A_{22}A_{33} - (A_{23}^2))(det A)^{-1}
        //                middle = (A_{11}A_{33} - (A_{13}^2))(det A)^{-1}
	//           lower right = (A_{11}A_{22} - (A_{12}^2))(det A)^{-1}
	dww = dxx * dxx;
	dvv = dyy * dyy;
	duu = dzz * dzz;
	determinant = 1 + 2 * dxx * dyy * dzz - dww - dvv - duu;
	if (fabs(determinant) < EPSILON) {
	  goto epi_linear_thread_regression_fail;
	}
	// (1 - x^2)/det > vif_thresh
	// if det > 0:
	//   1 - x^2 > vif_thresh * det
	//   1 - vif_thresh * det > x^2
	// otherwise:
	//   1 - x^2 < vif_thresh * det
	//   1 - vif_thresh * det < x^2
        dxx = 1 - vif_thresh * determinant; // now a threshold
	if (((determinant > 0) && ((dxx > dww) || (dxx > dvv) || (dxx > duu))) || ((determinant < 0) && ((dxx < dww) || (dxx < dvv) || (dxx < duu)))) {
	  goto epi_linear_thread_regression_fail;
	}

	// VIF check done, now perform linear regression
	dmatrix_buf[0] = sample_ctd;
	dmatrix_buf[1] = cur_sum_ad;
        dmatrix_buf[2] = cur_sum_bd;
	dmatrix_buf[3] = cur_sum_abd;
	dmatrix_buf[5] = (double)((int32_t)cur_sum_aa);
	dmatrix_buf[6] = cur_sum_abd;
	dmatrix_buf[7] = (double)((intptr_t)cur_sum_aab);
	dmatrix_buf[10] = (double)((int32_t)cur_sum_bb);
	dmatrix_buf[11] = (double)((intptr_t)cur_sum_abb);
	dmatrix_buf[15] = (double)((intptr_t)cur_sum_aabb);
	if (matrix_invert_4x4symm(dmatrix_buf)) {
	  goto epi_linear_thread_regression_fail;
	}

	for (param_idx = 0; param_idx < 4; param_idx++) {
	  dmatrix_buf2[param_idx] = sqrt(dmatrix_buf[param_idx * 5]);
	}
        for (param_idx = 1; param_idx < 4; param_idx++) {
          dxx = 0.99999 * dmatrix_buf2[param_idx];
          dptr = &(dmatrix_buf[param_idx * 4]);
          dptr2 = dmatrix_buf2;
          for (param_idx2 = 0; param_idx2 < param_idx; param_idx2++) {
            if ((*dptr++) > dxx * (*dptr2++)) {
              goto epi_linear_thread_regression_fail;
	    }
	  }
	}
        min_sigma = MAXV(dmatrix_buf[5], dmatrix_buf[10]);
	if (dmatrix_buf[15] > min_sigma) {
          min_sigma = dmatrix_buf[15];
	}
	min_sigma = 1e-20 / min_sigma;

	for (param_idx = 0; param_idx < 4; param_idx++) {
	  dptr = &(dmatrix_buf[param_idx * 4]);
	  dmatrix_buf2[param_idx] = cur_pheno_sum * dptr[0] + cur_sum_a_pheno * dptr[1] + cur_sum_b_pheno * dptr[2] + cur_sum_ab_pheno * dptr[3];
	}
	// dmatrix_buf2[0..3] now has linear regression result

	// partial = coef[0] + A * coef[1] + B * coef[2] + AB * coef[3] - pheno
	// sigma = \sum_{all samples} (partial * partial)
	//       = \sum (coef[0]^2
        //               + 2 * A * coef[0] * coef[1]
	//               + 2 * B * coef[0] * coef[2]
	//               + 2 * AB * coef[0] * coef[3]
	//               - 2 * coef[0] * pheno
	//               + AA * coef[1]^2
	//               + 2 * AB * coef[1] * coef[2]
	//               + 2 * AAB * coef[1] * coef[3]
	//               - 2 * A * coef[1] * pheno
	//               + BB * coef[2]^2
	//               + 2 * ABB * coef[2] * coef[3]
	//               - 2 * B * coef[2] * pheno
	//               + AABB * coef[3]^2
	//               - 2 * AB * coef[3] * pheno
	//               + pheno * pheno
	sigma = dmatrix_buf2[0] * dmatrix_buf2[0] * sample_ctd
	      + dmatrix_buf2[1] * dmatrix_buf2[1] * ((double)((int32_t)cur_sum_aa))
	      + dmatrix_buf2[2] * dmatrix_buf2[2] * ((double)((int32_t)cur_sum_bb))
	      + dmatrix_buf2[3] * dmatrix_buf2[3] * ((double)((intptr_t)cur_sum_aabb))
              + cur_pheno_ssq
	      + 2 * (dmatrix_buf2[0] * (dmatrix_buf2[1] * cur_sum_ad
                                      + dmatrix_buf2[2] * cur_sum_bd
                                      + dmatrix_buf2[3] * cur_sum_abd
				      - cur_pheno_sum)
		   + dmatrix_buf2[1] * (dmatrix_buf2[2] * cur_sum_abd
				      + dmatrix_buf2[3] * ((double)((intptr_t)cur_sum_aab))
				      - cur_sum_a_pheno)
		   + dmatrix_buf2[2] * (dmatrix_buf2[3] * ((double)((intptr_t)cur_sum_abb))
				      - cur_sum_b_pheno)
		   - dmatrix_buf2[3] * cur_sum_ab_pheno);
	sigma /= (double)((int32_t)(cur_sample_ct - 4));
        if (sigma < min_sigma) {
          goto epi_linear_thread_regression_fail;
	}

	// dmatrix_buf2[3] = linear regression beta for AB term
	// sqrt(dmatrix_buf[15] * sigma) = standard error for AB term
	dxx = dmatrix_buf2[3];
	zsq = (dxx * dxx) / (dmatrix_buf[15] * sigma);
	if (zsq >= alpha1sq) {
          all_chisq_write[2 * block_idx2] = zsq;
          all_chisq_write[2 * block_idx2 + 1] = dxx;
	}
	if (zsq >= alpha2sq) {
          n_sig_ct_fixed++;
	  n_sig_ct2[block_idx2] += 1;
	}
	if (zsq > best_chisq_fixed) {
          best_chisq_fixed = zsq;
	  best_id_fixed = block_idx2 + idx2_block_start;
	}
        dxx = *chisq2_ptr;
        if (zsq > dxx) {
          *chisq2_ptr = zsq;
	  best_id2[block_idx2] = marker_idx1;
	}
	while (0) {
	epi_linear_thread_regression_fail:
	  zsq = 0;
	  fail_ct_fixed++;
	  fail_ct2[block_idx2] += 1;
	  if (alpha1sq == 0.0) {
	    // special case: log NA when '--epi1 1' specified
	    all_chisq_write[block_idx2 * 2] = NAN;
	    all_chisq_write[block_idx2 * 2 + 1] = NAN;
	  }
	}
      }
      if (is_first_half) {
        is_first_half = 0;
	ulii = geno1_offsets[2 * block_idx1 + 1];
        cur_idx2_block_size = idx2_block_size;
        if (ulii < idx2_block_end) {
	  // guaranteed to be larger than idx2_block_start, otherwise there
	  // would have been no first half
          block_idx2 = ulii - idx2_block_start;
	  goto epi_linear_thread_second_half;
	}
      }
      if (best_chisq_fixed > best_chisq1[block_delta1]) {
        best_chisq1[block_delta1] = best_chisq_fixed;
	best_id1[block_delta1] = best_id_fixed;
      }
      n_sig_ct1[block_delta1] = n_sig_ct_fixed;
      if (fail_ct_fixed) {
        fail_ct1[block_delta1] = fail_ct_fixed;
      }
    }
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE epi_logistic_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t block_idx1_start = g_epi_idx1_block_bounds[tidx];
  uintptr_t block_idx1_end = g_epi_idx1_block_bounds[tidx + 1];
  uintptr_t idx1_block_start16 = g_epi_idx1_block_bounds16[tidx];
  uintptr_t marker_idx1 = g_epi_marker_idx1 + block_idx1_start;
  uintptr_t marker_ct = g_epi_marker_ct;
  float alpha1sq = (float)g_epi_alpha1sq[0];
  float alpha2sq = (float)g_epi_alpha2sq[0];
  uint32_t pheno_nm_ct = g_epi_pheno_nm_ct;
  uint32_t best_id_fixed = 0;
  uint32_t is_first_half = 0;
  uintptr_t pheno_nm_ctl2 = QUATERCT_TO_WORDCT(pheno_nm_ct);
  uintptr_t* geno1 = g_epi_geno1;
  uintptr_t* pheno_c = g_epi_pheno_c;
  float* covars_cov_major = g_epi_logistic_mt[tidx].cur_covars_cov_major;
  float* coef = g_epi_logistic_mt[tidx].coef;
  float* pp = g_epi_logistic_mt[tidx].pp;
  float* sample_1d_buf = g_epi_logistic_mt[tidx].sample_1d_buf;
  float* pheno_buf = g_epi_logistic_mt[tidx].pheno_buf;
  float* param_1d_buf = g_epi_logistic_mt[tidx].param_1d_buf;
  float* param_1d_buf2 = g_epi_logistic_mt[tidx].param_1d_buf2;
  float* param_2d_buf = g_epi_logistic_mt[tidx].param_2d_buf;
  float* param_2d_buf2 = g_epi_logistic_mt[tidx].param_2d_buf2;
  uint32_t* geno1_offsets = g_epi_geno1_offsets;
  uint32_t* best_id1 = &(g_epi_best_id1[idx1_block_start16]);
  uintptr_t* cur_geno1;
  uintptr_t* geno2;
  uintptr_t* cur_geno2;
  float* all_chisq_write;
  float* chisq2_ptr;
  float* all_chisq;
  float* best_chisq1;
  float* best_chisq2;
  float* fptr;
  float* fptr2;
  uint32_t* n_sig_ct1;
  uint32_t* fail_ct1;
  uint32_t* best_id2;
  uint32_t* n_sig_ct2;
  uint32_t* fail_ct2;
  uintptr_t idx2_block_size;
  uintptr_t cur_idx2_block_size;
  uintptr_t idx2_block_start;
  uintptr_t idx2_block_end;
  uintptr_t idx2_block_sizea16;
  uintptr_t block_idx1;
  uintptr_t block_delta1;
  uintptr_t block_idx2;
  uintptr_t cur_word1;
  uintptr_t cur_word2;
  uintptr_t param_idx;
  uintptr_t param_idx2;
  uintptr_t cur_sample_cta4;
  uintptr_t ulii;
  uintptr_t uljj;
  float best_chisq_fixed;
  // todo
  float fxx;
  float fyy;
  float zsq;
  uint32_t n_sig_ct_fixed;
  uint32_t fail_ct_fixed;
  uint32_t widx;
  uint32_t loop_end;
  uint32_t sample_idx;
  uint32_t cur_sample_ct;
  while (1) {
    idx2_block_size = g_epi_idx2_block_size;
    cur_idx2_block_size = idx2_block_size;
    idx2_block_start = g_epi_idx2_block_start;
    idx2_block_end = idx2_block_start + idx2_block_size;
    idx2_block_sizea16 = round_up_pow2(idx2_block_size, 16);
    geno2 = g_epi_geno2;
    all_chisq = &(g_epi_all_chisq_f[2 * idx2_block_start]);
    best_chisq1 = &(g_epi_best_chisq_f1[idx1_block_start16]);
    best_chisq2 = &(g_epi_best_chisq_f2[tidx * idx2_block_sizea16]);
    n_sig_ct1 = &(g_epi_n_sig_ct1[idx1_block_start16]);
    fail_ct1 = &(g_epi_fail_ct1[idx1_block_start16]);
    best_id2 = &(g_epi_best_id2[tidx * idx2_block_sizea16]);
    n_sig_ct2 = &(g_epi_n_sig_ct2[tidx * idx2_block_sizea16]);
    fail_ct2 = &(g_epi_fail_ct2[tidx * idx2_block_sizea16]);
    for (block_idx1 = block_idx1_start; block_idx1 < block_idx1_end; block_idx1++, marker_idx1++) {
      ulii = geno1_offsets[2 * block_idx1];
      if (ulii > idx2_block_start) {
        block_idx2 = 0;
        cur_idx2_block_size = ulii - idx2_block_start;
	if (cur_idx2_block_size >= idx2_block_size) {
          cur_idx2_block_size = idx2_block_size;
	} else {
	  is_first_half = 1;
        }
      } else {
        ulii = geno1_offsets[2 * block_idx1 + 1];
        if (ulii >= idx2_block_end) {
          // may not be done in set1 x all or set1 x set2 cases
          continue;
	} else {
          if (ulii <= idx2_block_start) {
            block_idx2 = 0;
	  } else {
            block_idx2 = ulii - idx2_block_start;
	  }
	}
      }
      cur_geno1 = &(geno1[block_idx1 * pheno_nm_ctl2]);
      n_sig_ct_fixed = 0;
      fail_ct_fixed = 0;
      block_delta1 = block_idx1 - block_idx1_start;
      best_chisq_fixed = best_chisq1[block_delta1];

      // [0] = chisq, [1] = beta
      all_chisq_write = &(all_chisq[block_idx1 * marker_ct * 2]);

    epi_logistic_thread_second_half:
      cur_geno2 = &(geno2[block_idx2 * pheno_nm_ctl2]);
      chisq2_ptr = &(best_chisq2[block_idx2]);
      for (; block_idx2 < cur_idx2_block_size; block_idx2++, chisq2_ptr++, cur_geno2 = &(cur_geno2[pheno_nm_ctl2])) {
        fptr = covars_cov_major;
	fptr2 = pheno_buf;
	cur_sample_ct = pheno_nm_ct;
	// this part is similar to glm_logistic().

	// 1. determine number of samples with at least one missing genotype
	for (widx = 0; widx < pheno_nm_ctl2; widx++) {
	  cur_word1 = cur_geno1[widx];
	  cur_word2 = cur_geno2[widx];
	  cur_word1 = cur_word1 & (cur_word1 >> 1);
	  cur_word2 = cur_word2 & (cur_word2 >> 1);
	  cur_sample_ct -= popcount2_long((cur_word1 | cur_word2) & FIVEMASK);
	}
        unsigned char geno_pair_present[12];
	if (cur_sample_ct <= 4) {
	  goto epi_logistic_thread_regression_fail;
	}
	// 2. now populate covariate-major matrix with 16-byte-aligned,
	//    trailing-entries-zeroed rows
        // quasi-bugfix (13 Sep 2018): reliably detect when this matrix is not
        // of full rank, and skip the regression in that case.
        memset(geno_pair_present, 0, 12);
	cur_sample_cta4 = round_up_pow2(cur_sample_ct, 4);
	for (widx = 0; widx < pheno_nm_ctl2; widx++) {
	  sample_idx = widx * BITCT2;
          cur_word1 = cur_geno1[widx];
          cur_word2 = cur_geno2[widx];
          loop_end = sample_idx + BITCT2;
	  if (loop_end > pheno_nm_ct) {
            loop_end = pheno_nm_ct;
	  }
          for (; sample_idx < loop_end; sample_idx++) {
            ulii = cur_word1 & (3 * ONELU);
            uljj = cur_word2 & (3 * ONELU);
	    if ((ulii != 3) && (uljj != 3)) {
              *fptr = 1.0;
              geno_pair_present[ulii + uljj * 4] = 1;
	      fxx = (float)((intptr_t)ulii);
	      fyy = (float)((intptr_t)uljj);
	      // maybe this is faster with continuous writes instead of
	      // continuous reads?  can experiment later
              fptr[cur_sample_cta4] = fxx;
	      fptr[2 * cur_sample_cta4] = fyy;
              fptr[3 * cur_sample_cta4] = fxx * fyy;
	      fptr++;
	      *fptr2++ = (float)((int32_t)is_set(pheno_c, sample_idx));
	    }
            cur_word1 >>= 2;
            cur_word2 >>= 2;
	  }
	}
        if (!geno_pair_present[5]) {
          // not full rank if any 2x2 square in the 3x3 contingency table is
          // empty.
          if (((!geno_pair_present[0]) && (!geno_pair_present[1]) && (!geno_pair_present[4])) ||
              ((!geno_pair_present[1]) && (!geno_pair_present[2]) && (!geno_pair_present[6])) ||
              ((!geno_pair_present[4]) && (!geno_pair_present[8]) && (!geno_pair_present[9])) ||
              ((!geno_pair_present[6]) && (!geno_pair_present[9]) && (!geno_pair_present[10]))) {
            goto epi_logistic_thread_regression_fail;
          }
        }
	if (cur_sample_ct < cur_sample_cta4) {
	  loop_end = cur_sample_cta4 - cur_sample_ct;
	  fill_float_zero(loop_end, fptr);
	  fill_float_zero(loop_end, &(fptr[cur_sample_cta4]));
	  fill_float_zero(loop_end, &(fptr[2 * cur_sample_cta4]));
	  fill_float_zero(loop_end, &(fptr[3 * cur_sample_cta4]));
	  fill_float_zero(loop_end, fptr2);
	}

	fill_float_zero(4, coef);
	if (logistic_regression(cur_sample_ct, 4, sample_1d_buf, param_2d_buf, param_1d_buf, param_2d_buf2, param_1d_buf2, covars_cov_major, pheno_buf, coef, pp)) {
          goto epi_logistic_thread_regression_fail;
	}

	// compute S
	for (param_idx = 0; param_idx < 4; param_idx++) {
          fill_float_zero(4, param_1d_buf);
          param_1d_buf[param_idx] = 1.0;
	  solve_linear_system(param_2d_buf2, param_1d_buf, param_1d_buf2, 4);
	  memcpy(&(param_2d_buf[param_idx * 4]), param_1d_buf2, 4 * sizeof(float));
	}
	for (param_idx = 1; param_idx < 4; param_idx++) {
          fxx = param_2d_buf[param_idx * 5];
	  if ((fxx < 1e-20) || (!realnum(fxx))) {
	    goto epi_logistic_thread_regression_fail;
	  }
          param_2d_buf2[param_idx] = sqrtf(fxx);
	}
	param_2d_buf2[0] = sqrtf(param_2d_buf[0]);
	for (param_idx = 1; param_idx < 4; param_idx++) {
          fxx = 0.99999 * param_2d_buf2[param_idx];
	  fptr = &(param_2d_buf[param_idx * 4]);
	  fptr2 = param_2d_buf2;
	  for (param_idx2 = 0; param_idx2 < param_idx; param_idx2++) {
            if ((*fptr++) > fxx * (*fptr2++)) {
	      goto epi_logistic_thread_regression_fail;
	    }
	  }
	}

	// coef[3] = logistic regression beta for AB term
	// sqrt(param_2d_buf[15]) = standard error for AB term
	zsq = coef[3] * coef[3] / param_2d_buf[15];
	if (zsq >= alpha1sq) {
          all_chisq_write[2 * block_idx2] = zsq;
          all_chisq_write[2 * block_idx2 + 1] = coef[3];
	}
	if (zsq >= alpha2sq) {
          n_sig_ct_fixed++;
	  n_sig_ct2[block_idx2] += 1;
	}
	if (zsq > best_chisq_fixed) {
          best_chisq_fixed = zsq;
	  best_id_fixed = block_idx2 + idx2_block_start;
	}
        fxx = *chisq2_ptr;
        if (zsq > fxx) {
          *chisq2_ptr = zsq;
	  best_id2[block_idx2] = marker_idx1;
	}
	while (0) {
	epi_logistic_thread_regression_fail:
	  zsq = 0;
	  fail_ct_fixed++;
	  fail_ct2[block_idx2] += 1;
	  if (alpha1sq == 0.0) {
	    // special case: log NA when '--epi1 1' specified
	    all_chisq_write[block_idx2 * 2] = NAN;
	    all_chisq_write[block_idx2 * 2 + 1] = NAN;
	  }
	}
      }
      if (is_first_half) {
        is_first_half = 0;
	ulii = geno1_offsets[2 * block_idx1 + 1];
        cur_idx2_block_size = idx2_block_size;
        if (ulii < idx2_block_end) {
          block_idx2 = ulii - idx2_block_start;
	  goto epi_logistic_thread_second_half;
	}
      }
      if (best_chisq_fixed > best_chisq1[block_delta1]) {
        best_chisq1[block_delta1] = best_chisq_fixed;
	best_id1[block_delta1] = best_id_fixed;
      }
      n_sig_ct1[block_delta1] = n_sig_ct_fixed;
      if (fail_ct_fixed) {
        fail_ct1[block_delta1] = fail_ct_fixed;
      }
    }
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

double calc_lnlike(double known11, double known12, double known21, double known22, double center_ct_d, double freq11, double freq12, double freq21, double freq22, double half_hethet_share, double freq11_incr) {
  double lnlike;
  freq11 += freq11_incr;
  freq22 += freq11_incr;
  freq12 += half_hethet_share - freq11_incr;
  freq21 += half_hethet_share - freq11_incr;
  lnlike = center_ct_d * log(freq11 * freq22 + freq12 * freq21);
  if (known11 != 0.0) {
    lnlike += known11 * log(freq11);
  }
  if (known12 != 0.0) {
    lnlike += known12 * log(freq12);
  }
  if (known21 != 0.0) {
    lnlike += known21 * log(freq21);
  }
  if (known22 != 0.0) {
    lnlike += known22 * log(freq22);
  }
  return lnlike;
}

uint32_t em_phase_hethet(double known11, double known12, double known21, double known22, uint32_t center_ct, double* freq1x_ptr, double* freq2x_ptr, double* freqx1_ptr, double* freqx2_ptr, double* freq11_ptr, uint32_t* onside_sol_ct_ptr) {
  // Returns 1 if at least one SNP is monomorphic over all valid observations;
  // returns 0 otherwise, and fills all frequencies using the maximum
  // likelihood solution to the cubic equation.
  // (We're discontinuing most use of EM phasing since better algorithms have
  // been developed, but the two marker case is mathematically clean and fast
  // enough that it'll probably remain useful as an input for some of those
  // better algorithms...)
  double center_ct_d = (int32_t)center_ct;
  double twice_tot = known11 + known12 + known21 + known22 + 2 * center_ct_d;
  uint32_t sol_start_idx = 0;
  uint32_t sol_end_idx = 1;
  double solutions[3];
  double twice_tot_recip;
  double half_hethet_share;
  double freq11;
  double freq12;
  double freq21;
  double freq22;
  double prod_1122;
  double prod_1221;
  double incr_1122;
  double best_sol;
  double best_lnlike;
  double cur_lnlike;
  double freq1x;
  double freq2x;
  double freqx1;
  double freqx2;
  double lbound;
  double dxx;
  uint32_t cur_sol_idx;
  // shouldn't have to worry about subtractive cancellation problems here
  if (twice_tot == 0.0) {
    return 1;
  }
  twice_tot_recip = 1.0 / twice_tot;
  freq11 = known11 * twice_tot_recip;
  freq12 = known12 * twice_tot_recip;
  freq21 = known21 * twice_tot_recip;
  freq22 = known22 * twice_tot_recip;
  prod_1122 = freq11 * freq22;
  prod_1221 = freq12 * freq21;
  half_hethet_share = center_ct_d * twice_tot_recip;
  // the following four values should all be guaranteed nonzero except in the
  // NAN case
  freq1x = freq11 + freq12 + half_hethet_share;
  freq2x = 1.0 - freq1x;
  freqx1 = freq11 + freq21 + half_hethet_share;
  freqx2 = 1.0 - freqx1;
  if (center_ct) {
    if ((prod_1122 != 0.0) || (prod_1221 != 0.0)) {
      sol_end_idx = cubic_real_roots(0.5 * (freq11 + freq22 - freq12 - freq21 - 3 * half_hethet_share), 0.5 * (prod_1122 + prod_1221 + half_hethet_share * (freq12 + freq21 - freq11 - freq22 + half_hethet_share)), -0.5 * half_hethet_share * prod_1122, solutions);
      while (sol_end_idx && (solutions[sol_end_idx - 1] > half_hethet_share + SMALLISH_EPSILON)) {
	sol_end_idx--;
      }
      while ((sol_start_idx < sol_end_idx) && (solutions[sol_start_idx] < -SMALLISH_EPSILON)) {
	sol_start_idx++;
      }
      if (sol_start_idx == sol_end_idx) {
	// Lost a planet Master Obi-Wan has.  How embarrassing...
	// lost root must be a double root at one of the boundary points, just
	// check their likelihoods
	sol_start_idx = 0;
	sol_end_idx = 2;
	solutions[0] = 0;
	solutions[1] = half_hethet_share;
      } else {
	if (solutions[sol_start_idx] < 0) {
	  solutions[sol_start_idx] = 0;
	}
	if (solutions[sol_end_idx - 1] > half_hethet_share) {
	  solutions[sol_end_idx - 1] = half_hethet_share;
	}
      }
    } else {
      solutions[0] = 0;
      // bugfix (6 Oct 2017): need to use all nonzero values here
      const double nonzero_freq_xx = freq11 + freq22;
      const double nonzero_freq_xy = freq12 + freq21;
      if ((nonzero_freq_xx + SMALLISH_EPSILON < half_hethet_share + nonzero_freq_xy) && (nonzero_freq_xy + SMALLISH_EPSILON < half_hethet_share + nonzero_freq_xx)) {
	sol_end_idx = 3;
	solutions[1] = (half_hethet_share + nonzero_freq_xy - nonzero_freq_xx) * 0.5;
	solutions[2] = half_hethet_share;
      } else {
	sol_end_idx = 2;
	solutions[1] = half_hethet_share;
      }
    }
    best_sol = solutions[sol_start_idx];
    if (sol_end_idx > sol_start_idx + 1) {
      // select largest log likelihood
      best_lnlike = calc_lnlike(known11, known12, known21, known22, center_ct_d, freq11, freq12, freq21, freq22, half_hethet_share, best_sol);
      cur_sol_idx = sol_start_idx + 1;
      do {
	incr_1122 = solutions[cur_sol_idx];
        cur_lnlike = calc_lnlike(known11, known12, known21, known22, center_ct_d, freq11, freq12, freq21, freq22, half_hethet_share, incr_1122);
	if (cur_lnlike > best_lnlike) {
          cur_lnlike = best_lnlike;
          best_sol = incr_1122;
	}
      } while (++cur_sol_idx < sol_end_idx);
    }
    if (onside_sol_ct_ptr && (sol_end_idx > sol_start_idx + 1)) {
      if (freqx1 * freq1x >= freq11) {
	dxx = freq1x * freqx1 - freq11;
	if (dxx > half_hethet_share) {
	  dxx = half_hethet_share;
	}
      } else {
	dxx = 0.0;
      }
      // okay to NOT count suboptimal boundary points because they don't permit
      // direction changes within the main interval
      // this should exactly match haploview_blocks_classify()'s D sign check
      if ((freq11 + best_sol) - freqx1 * freq1x >= 0.0) {
	if (best_sol > dxx + SMALLISH_EPSILON) {
          lbound = dxx + SMALLISH_EPSILON;
	} else {
	  lbound = dxx;
	}
	if (best_sol < half_hethet_share - SMALLISH_EPSILON) {
	  half_hethet_share -= SMALLISH_EPSILON;
	}
      } else {
	if (best_sol > SMALLISH_EPSILON) {
	  lbound = SMALLISH_EPSILON;
	} else {
	  lbound = 0.0;
	}
	if (best_sol < dxx - SMALLISH_EPSILON) {
	  half_hethet_share = dxx - SMALLISH_EPSILON;
	} else {
	  half_hethet_share = dxx;
	}
      }
      for (cur_sol_idx = sol_start_idx; cur_sol_idx < sol_end_idx; cur_sol_idx++) {
	if (solutions[cur_sol_idx] < lbound) {
	  sol_start_idx++;
	}
	if (solutions[cur_sol_idx] > half_hethet_share) {
	  break;
	}
      }
      if (cur_sol_idx >= sol_start_idx + 2) {
	*onside_sol_ct_ptr = cur_sol_idx - sol_start_idx;
      }
    }
    freq11 += best_sol;
  } else if ((prod_1122 == 0.0) && (prod_1221 == 0.0)) {
    return 1;
  }
  *freq1x_ptr = freq1x;
  *freq2x_ptr = freq2x;
  *freqx1_ptr = freqx1;
  *freqx2_ptr = freqx2;
  *freq11_ptr = freq11;
  return 0;
}

uint32_t em_phase_hethet_nobase(uint32_t* counts, uint32_t is_x1, uint32_t is_x2, double* freq1x_ptr, double* freq2x_ptr, double* freqx1_ptr, double* freqx2_ptr, double* freq11_ptr) {
  // if is_x1 and/or is_x2 is set, counts[9]..[17] are male-only counts.
  double known11 = (double)(2 * counts[0] + counts[1] + counts[3]);
  double known12 = (double)(2 * counts[2] + counts[1] + counts[5]);
  double known21 = (double)(2 * counts[6] + counts[3] + counts[7]);
  double known22 = (double)(2 * counts[8] + counts[5] + counts[7]);
  if (is_x1 || is_x2) {
    if (is_x1 && is_x2) {
      known11 -= (double)((int32_t)counts[9]);
      known12 -= (double)((int32_t)counts[11]);
      known21 -= (double)((int32_t)counts[15]);
      known22 -= (double)((int32_t)counts[17]);
    } else if (is_x1) {
      known11 -= ((double)(2 * counts[9] + counts[10])) * (1.0 - SQRT_HALF);
      known12 -= ((double)(2 * counts[11] + counts[10])) * (1.0 - SQRT_HALF);
      known21 -= ((double)(2 * counts[15] + counts[16])) * (1.0 - SQRT_HALF);
      known22 -= ((double)(2 * counts[17] + counts[16])) * (1.0 - SQRT_HALF);
    } else {
      known11 -= ((double)(2 * counts[9] + counts[12])) * (1.0 - SQRT_HALF);
      known12 -= ((double)(2 * counts[11] + counts[12])) * (1.0 - SQRT_HALF);
      known21 -= ((double)(2 * counts[15] + counts[14])) * (1.0 - SQRT_HALF);
      known22 -= ((double)(2 * counts[17] + counts[14])) * (1.0 - SQRT_HALF);
    }
  }
  return em_phase_hethet(known11, known12, known21, known22, counts[4], freq1x_ptr, freq2x_ptr, freqx1_ptr, freqx2_ptr, freq11_ptr, nullptr);
}

THREAD_RET_TYPE ld_dprime_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t block_idx1_start = (tidx * g_ld_idx1_block_size) / g_ld_thread_ct;
  uintptr_t block_idx1_end = ((tidx + 1) * g_ld_idx1_block_size) / g_ld_thread_ct;
  uintptr_t marker_idx2_maxw = g_ld_marker_ctm8;
  uintptr_t founder_ct = g_ld_founder_ct;
  uint32_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(founder_ct);
  uint32_t founder_ctsplit = 3 * founder_ctv3;
  uintptr_t* geno1 = g_ld_geno1;
  uintptr_t* zmiss1 = g_epi_zmiss1;
  uintptr_t* sex_male = g_ld_sex_male;
  uintptr_t* cur_geno1_male = nullptr;
  uint32_t* ld_interval1 = g_ld_interval1;
  uint32_t is_dprime = g_ld_modifier & (LD_DPRIME | LD_DPRIME_SIGNED);
  uint32_t is_dprime_unsigned = g_ld_modifier & LD_DPRIME;
  uint32_t is_r2 = g_ld_is_r2;
  uint32_t xstart1 = g_ld_xstart1;
  uint32_t xend1 = g_ld_xend1;
  double* results = g_ld_results;
  uint32_t tot1[6];
  uint32_t counts[18];
  uintptr_t* cur_geno1;
  uintptr_t* cur_geno2;
  uintptr_t* geno2;
  uintptr_t* zmiss2;
  double* rptr;
  uint32_t* tot2;
  uint32_t* cur_tot2;
  uintptr_t idx2_block_size;
  uintptr_t idx2_block_start;
  uintptr_t block_idx1;
  uintptr_t block_idx2;
  uintptr_t cur_zmiss2;
  uintptr_t cur_block_idx2_end;
  double freq11;
  double freq11_expected;
  double freq1x;
  double freq2x;
  double freqx1;
  double freqx2;
  double dxx;
  uint32_t xstart2;
  uint32_t xend2;
  uint32_t x2_present;
  uint32_t is_x1;
  uint32_t is_x2;
  uint32_t nm_fixed;
  if (g_ld_thread_wkspace) {
    cur_geno1_male = &(g_ld_thread_wkspace[tidx * round_up_pow2(founder_ctsplit, CACHELINE_WORD)]);
  }
  // suppress warning
  fill_uint_zero(3, &(tot1[3]));
  while (1) {
    idx2_block_size = g_ld_idx2_block_size;
    idx2_block_start = g_ld_idx2_block_start;
    geno2 = g_ld_geno2;
    zmiss2 = g_epi_zmiss2;
    tot2 = g_epi_tot2;
    xstart2 = g_ld_xstart2;
    xend2 = g_ld_xend2;
    x2_present = (g_ld_thread_wkspace && (idx2_block_start < xend2) && (idx2_block_start + idx2_block_size > xstart2));
    for (block_idx1 = block_idx1_start; block_idx1 < block_idx1_end; block_idx1++) {
      cur_zmiss2 = ld_interval1[block_idx1 * 2];
      block_idx2 = cur_zmiss2;
      cur_block_idx2_end = ld_interval1[block_idx1 * 2 + 1];
      if (block_idx2 < idx2_block_start) {
	if (cur_block_idx2_end <= idx2_block_start) {
	  continue;
	}
	block_idx2 = 0;
      } else {
	block_idx2 -= idx2_block_start;
	if (block_idx2 >= idx2_block_size) {
	  break;
	}
      }
      cur_block_idx2_end -= idx2_block_start;
      if (cur_block_idx2_end > idx2_block_size) {
	cur_block_idx2_end = idx2_block_size;
      }
      is_x1 = (block_idx1 >= xstart1) && (block_idx1 < xend1);
      nm_fixed = is_set_ul(zmiss1, block_idx1);
      cur_geno1 = &(geno1[block_idx1 * founder_ctsplit]);
      tot1[0] = popcount_longs(cur_geno1, founder_ctv3);
      tot1[1] = popcount_longs(&(cur_geno1[founder_ctv3]), founder_ctv3);
      tot1[2] = popcount_longs(&(cur_geno1[2 * founder_ctv3]), founder_ctv3);
      if (is_x1 || x2_present) {
	memcpy(cur_geno1_male, cur_geno1, founder_ctsplit * sizeof(intptr_t));
        bitvec_and(sex_male, founder_ctv3, cur_geno1_male);
        tot1[3] = popcount_longs(cur_geno1_male, founder_ctv3);
        bitvec_and(sex_male, founder_ctv3, &(cur_geno1_male[founder_ctv3]));
	tot1[4] = popcount_longs(&(cur_geno1_male[founder_ctv3]), founder_ctv3);
        bitvec_and(sex_male, founder_ctv3, &(cur_geno1_male[2 * founder_ctv3]));
	tot1[5] = popcount_longs(&(cur_geno1_male[2 * founder_ctv3]), founder_ctv3);
      }
      cur_geno2 = &(geno2[block_idx2 * founder_ctsplit]);
      rptr = &(results[2 * block_idx1 * marker_idx2_maxw]);
      for (; block_idx2 < cur_block_idx2_end; block_idx2++, cur_geno2 = &(cur_geno2[founder_ctsplit])) {
	cur_tot2 = &(tot2[block_idx2 * 3]);
	cur_zmiss2 = is_set_ul(zmiss2, block_idx2);
	if (nm_fixed) {
	  two_locus_count_table_zmiss1(cur_geno1, cur_geno2, counts, founder_ctv3, cur_zmiss2);
	  if (cur_zmiss2) {
	    counts[2] = tot1[0] - counts[0] - counts[1];
	    counts[5] = tot1[1] - counts[3] - counts[4];
	  }
	  counts[6] = cur_tot2[0] - counts[0] - counts[3];
	  counts[7] = cur_tot2[1] - counts[1] - counts[4];
	  counts[8] = cur_tot2[2] - counts[2] - counts[5];
	} else {
	  two_locus_count_table(cur_geno1, cur_geno2, counts, founder_ctv3, cur_zmiss2);
	  if (cur_zmiss2) {
	    counts[2] = tot1[0] - counts[0] - counts[1];
	    counts[5] = tot1[1] - counts[3] - counts[4];
	    counts[8] = tot1[2] - counts[6] - counts[7];
	  }
	}
	is_x2 = ((block_idx2 < xend2) && (block_idx2 >= xstart2));
	if (is_x1 || is_x2) {
	  two_locus_count_table(cur_geno1_male, cur_geno2, &(counts[9]), founder_ctv3, cur_zmiss2);
	  if (cur_zmiss2) {
	    counts[11] = tot1[3] - counts[9] - counts[10];
	    counts[14] = tot1[4] - counts[12] - counts[13];
	    counts[17] = tot1[5] - counts[15] - counts[16];
	  }
	}
	if (em_phase_hethet_nobase(counts, is_x1, is_x2, &freq1x, &freq2x, &freqx1, &freqx2, &freq11)) {
	  *rptr++ = NAN;
	  *rptr++ = NAN;
	  continue;
	}
	freq11_expected = freqx1 * freq1x; // fA * fB temp var
	// a bit of numeric instability here, but not tragic since this is the
	// end of the calculation
	dxx = freq11 - freq11_expected; // D
	if (fabs(dxx) < SMALL_EPSILON) {
	  *rptr++ = 0;
	  *rptr = 0;
	} else {
	  if (is_r2) {
	    *rptr = fabs(dxx) * dxx / (freq11_expected * freq2x * freqx2);
	  } else {
	    *rptr = dxx / sqrt(freq11_expected * freq2x * freqx2);
	  }
	  rptr++;
	  if (is_dprime) {
	    if (dxx >= 0) {
	      dxx /= MINV(freqx1 * freq2x, freqx2 * freq1x);
	    } else {
	      if (is_dprime_unsigned) {
		dxx = -dxx;
	      }
	      dxx /= MINV(freq11_expected, freqx2 * freq2x);
	    }
	  }
	  *rptr = dxx;
	}
	rptr++;
      }
    }
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

int32_t ld_report_dprime(pthread_t* threads, Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uintptr_t* sex_male, uintptr_t* founder_include2, uintptr_t* founder_male_include2, uintptr_t* loadbuf_raw, char* outname, uint32_t hh_exists, uintptr_t marker_idx1_start, uintptr_t marker_idx1_end) {
  Chrom_info* chrom_info_ptr = g_ld_chrom_info_ptr;
  uintptr_t* marker_exclude_idx1 = g_ld_marker_exclude_idx1;
  uintptr_t* marker_exclude = g_ld_marker_exclude;
  uint32_t* marker_pos = g_ld_marker_pos;
  double* marker_cms = g_ld_marker_cms;
  uintptr_t marker_ct = g_ld_marker_ct;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t founder_ct = g_ld_founder_ct;
  uintptr_t founder_ctl = BITCT_TO_WORDCT(founder_ct);
  uintptr_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(founder_ct);
  uintptr_t founder_ctsplit = 3 * founder_ctv3;
  uintptr_t final_mask = get_final_mask(founder_ct);
  uintptr_t orig_marker_ctm8 = g_ld_marker_ctm8;
  uintptr_t marker_idx2_maxw = orig_marker_ctm8;
  uintptr_t marker_idx1 = marker_idx1_start;
  uintptr_t job_size = marker_idx1_end - marker_idx1_start;
  uintptr_t pct_thresh = job_size / 100;
  uintptr_t pct = 1;
  uintptr_t ulii = founder_ctsplit * sizeof(intptr_t) + 2 * sizeof(int32_t) + marker_idx2_maxw * 2 * sizeof(double);
  uint32_t output_gz = ldip->modifier & LD_REPORT_GZ;
  uint32_t is_inter_chr = g_ld_is_inter_chr;
  uint32_t idx1_subset = (ldip->snpstr || ldip->snps_rl.name_ct);
  uint32_t window_size_m1 = ldip->window_size - 1;
  uint32_t window_bp = ldip->window_bp;
  double window_cm = ldip->window_cm;
  uint32_t thread_ct = g_ld_thread_ct;
  uint32_t chrom_fo_idx = 0;
  uint32_t is_haploid = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t not_first_write = 0;
  uint32_t chrom_last = 0;
  uint32_t marker_uidx2_back = 0;
  uint32_t marker_uidx2_fwd = 0;
  uint32_t marker_uidx2_fwd2 = 0;
  uint32_t window_trail_ct = 0;
  uint32_t window_lead_ct = 0;
  int32_t x_code = chrom_info_ptr->xymt_codes[X_OFFSET];
  uint32_t xstart = 0;
  uint32_t xend = 0;
  int32_t retval = 0;
  uintptr_t* loadbuf;
  uintptr_t* dummy_nm;
  uintptr_t* ulptr;
  uint32_t* uiptr;
  unsigned char* overflow_buf;
  unsigned char* bigstack_mark2;
  uintptr_t cur_bigstack_left;
  uintptr_t thread_workload;
  uintptr_t idx1_block_size;
  uintptr_t idx2_block_size;
  uintptr_t cur_idx2_block_size;
  uintptr_t marker_idx2;
  uintptr_t marker_uidx1;
  uintptr_t marker_uidx1_tmp;
  uintptr_t marker_uidx2_base;
  uintptr_t marker_uidx2;
  uintptr_t marker_idx2_base;
  uintptr_t marker_idx2_end;
  uintptr_t block_idx1;
  uintptr_t block_idx2;
  uintptr_t uljj;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t cur_marker_pos;
  double cur_marker_cm;
  uint32_t is_last_block;
  uint32_t uii;
  if (bigstack_alloc_uc(262144, &overflow_buf) ||
      bigstack_alloc_ul(founder_ctl * 2, &loadbuf) ||
      bigstack_alloc_ul(founder_ctl, &dummy_nm)) {
    goto ld_report_dprime_ret_NOMEM;
  }
  loadbuf[founder_ctl * 2 - 2] = 0;
  loadbuf[founder_ctl * 2 - 1] = 0;
  fill_all_bits(founder_ct, dummy_nm);
  g_ld_thread_wkspace = nullptr;
  if ((x_code != -2) && is_set(chrom_info_ptr->chrom_mask, x_code)) {
    uii = get_chrom_start_vidx(chrom_info_ptr, (uint32_t)x_code);
    chrom_end = get_chrom_end_vidx(chrom_info_ptr, (uint32_t)x_code);
    chrom_end = chrom_end - uii - popcount_bit_idx(marker_exclude, uii, chrom_end);
    if (chrom_end) {
      if (bigstack_alloc_ul(round_up_pow2(founder_ctsplit, CACHELINE_WORD) * thread_ct, &g_ld_thread_wkspace)) {
	goto ld_report_dprime_ret_NOMEM;
      }
      xstart = uii - popcount_bit_idx(marker_exclude, 0, uii);
      xend = xstart + chrom_end;
      g_ld_sex_male = sex_male;
    }
  }
  cur_bigstack_left = bigstack_left();
  if (cur_bigstack_left < 2 * CACHELINE) {
    goto ld_report_dprime_ret_NOMEM;
  }
  idx1_block_size = (cur_bigstack_left - 2 * CACHELINE) / (ulii * 2 + 1);
  thread_workload = idx1_block_size / thread_ct;
  if (!thread_workload) {
    goto ld_report_dprime_ret_NOMEM;
  }
  idx1_block_size = thread_workload * thread_ct;
  if (idx1_block_size > job_size) {
    idx1_block_size = job_size;
  }
  if (bigstack_alloc_ul(founder_ctsplit * idx1_block_size, &g_ld_geno1) ||
      bigstack_alloc_ul(BITCT_TO_WORDCT(idx1_block_size), &g_epi_zmiss1) ||
      bigstack_alloc_ui(idx1_block_size * 2, &g_ld_interval1) ||
      // double size since both r/r^2 and dprime are needed
      // (marker_idx2_maxw only needs to be divisible by 4 as a result)
      bigstack_alloc_d(marker_idx2_maxw * 2 * idx1_block_size, &g_ld_results)) {
    goto ld_report_dprime_ret_NOMEM;
  }
  for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
    g_ld_geno1[block_idx1 * founder_ctsplit + founder_ctv3 - 1] = 0;
    g_ld_geno1[block_idx1 * founder_ctsplit + 2 * founder_ctv3 - 1] = 0;
    g_ld_geno1[block_idx1 * founder_ctsplit + founder_ctsplit - 1] = 0;
  }

  ulii = founder_ctsplit * sizeof(intptr_t) + 1 + 3 * sizeof(int32_t);
  cur_bigstack_left = bigstack_left();
  if (cur_bigstack_left >= CACHELINE) {
    cur_bigstack_left -= CACHELINE;
  }
  idx2_block_size = (cur_bigstack_left / ulii) & (~(7 * ONELU));
  if (idx2_block_size > marker_ct) {
    idx2_block_size = round_up_pow2(marker_ct, 8);
  }
  bigstack_mark2 = g_bigstack_base;
  while (1) {
    if (!idx2_block_size) {
      goto ld_report_dprime_ret_NOMEM;
    }
    if (!(bigstack_alloc_ul(founder_ctsplit * idx2_block_size, &g_ld_geno2) ||
          bigstack_alloc_ul(BITCT_TO_WORDCT(idx2_block_size), &g_epi_zmiss2) ||
          bigstack_alloc_ui(idx2_block_size * 3, &g_epi_tot2))) {
      break;
    }
    bigstack_reset(bigstack_mark2);
    idx2_block_size -= 4;
  }
  for (block_idx2 = 0; block_idx2 < idx2_block_size; block_idx2++) {
    g_ld_geno2[block_idx2 * founder_ctsplit + founder_ctv3 - 1] = 0;
    g_ld_geno2[block_idx2 * founder_ctsplit + 2 * founder_ctv3 - 1] = 0;
    g_ld_geno2[block_idx2 * founder_ctsplit + founder_ctsplit - 1] = 0;
  }
  marker_uidx1 = next_unset_unsafe(marker_exclude_idx1, 0);
  if (marker_idx1) {
    marker_uidx1 = jump_forward_unset_unsafe(marker_exclude_idx1, marker_uidx1 + 1, marker_idx1);
  }
  LOGPRINTF("--r%s%s%s d%s%s...", g_ld_is_r2? "2" : "", is_inter_chr? " inter-chr" : "", g_ld_marker_allele_ptrs? " in-phase" : "", (g_ld_modifier & LD_D)? "" : ((g_ld_modifier & LD_DPRIME)? "prime" : "prime-signed"), g_ld_set_allele_freqs? " with-freqs" : "");
  fputs(" 0%", stdout);
  while (1) {
    fputs(" [processing]", stdout);
    fflush(stdout);
    if (idx1_block_size > marker_idx1_end - marker_idx1) {
      idx1_block_size = marker_idx1_end - marker_idx1;
      if (idx1_block_size < thread_ct) {
        thread_ct = idx1_block_size;
        g_ld_thread_ct = thread_ct;
      }
    }
    g_ld_idx1_block_size = idx1_block_size;
    marker_uidx1_tmp = marker_uidx1;
    if ((marker_idx1 < xend) && (marker_idx1 + idx1_block_size > xstart)) {
      uii = MAXV(marker_idx1, xstart);
      g_ld_xstart1 = uii - marker_idx1;
      g_ld_xend1 = MINV(xend, marker_idx1 + idx1_block_size) - uii;
    }

    if (idx1_subset) {
      if (!is_inter_chr) {
	chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx1_tmp);
	marker_uidx2_base = window_back(marker_pos, marker_cms, marker_exclude, next_unset_unsafe(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx]), marker_uidx1, window_size_m1, window_bp, window_cm, &uii);
	marker_idx2_base = marker_uidx2_base - popcount_bit_idx(marker_exclude, 0, marker_uidx2_base);
	marker_idx2 = marker_idx2_base + uii;
      } else {
	marker_uidx2_base = next_unset_unsafe(marker_exclude, 0);
	marker_idx2_base = 0;
	marker_idx2 = 0;
      }
    } else {
      marker_idx2_base = marker_uidx1 + 1 - popcount_bit_idx(marker_exclude, 0, marker_uidx1);
      if (marker_idx2_base == marker_ct) {
	goto ld_report_dprime_done;
      }
      marker_idx2 = marker_idx2_base - 1;
      marker_uidx2_base = next_unset_unsafe(marker_exclude, marker_uidx1 + 1);
    }
    if (fseeko(bedfile, bed_offset + (marker_uidx1 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto ld_report_dprime_ret_READ_FAIL;
    }
    chrom_end = 0;
    fill_ulong_zero(BITCT_TO_WORDCT(idx1_block_size), g_epi_zmiss1);
    for (block_idx1 = 0; block_idx1 < idx1_block_size; marker_uidx1_tmp++, block_idx1++, marker_idx2++) {
      if (IS_SET(marker_exclude_idx1, marker_uidx1_tmp)) {
        ulii = next_unset_ul_unsafe(marker_exclude_idx1, marker_uidx1_tmp);
	uljj = ulii - marker_uidx1_tmp - popcount_bit_idx(marker_exclude, marker_uidx1_tmp, ulii);
	if (uljj) {
	  uii = 1;
	  marker_idx2 += uljj;
	}
	marker_uidx1_tmp = ulii;
        if (fseeko(bedfile, bed_offset + (marker_uidx1_tmp * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
          goto ld_report_dprime_ret_READ_FAIL;
	}
      }
      if (marker_uidx1_tmp >= chrom_end) {
        chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx1_tmp);
        chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
        chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
	chrom_last = prev_unset_unsafe(marker_exclude, chrom_end);
        is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
	is_x = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[X_OFFSET]);
	is_y = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
	uii = 1;
      }
      if (!is_inter_chr) {
	// uii == 0 if we can perform an incremental update, 1 if we need
	// fully-powered window_back()/window_forward()
	if (uii) {
	  if (idx1_subset) {
	    marker_uidx2_back = window_back(marker_pos, marker_cms, marker_exclude, next_unset_unsafe(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx]), marker_uidx1_tmp, window_size_m1, window_bp, window_cm, &window_trail_ct);
	  }
	  marker_uidx2_fwd = window_forward(marker_pos, marker_cms, marker_exclude, marker_uidx1_tmp, chrom_last, window_size_m1, window_bp, window_cm, &window_lead_ct);
	  marker_uidx2_fwd2 = marker_uidx2_fwd;
	  if (marker_uidx2_fwd < chrom_last) {
	    marker_uidx2_fwd2++;
	    next_unset_unsafe_ck(marker_exclude, &marker_uidx2_fwd2);
	  }
	  uii = 0;
	} else {
	  if (idx1_subset) {
	    if (window_trail_ct == window_size_m1) {
	      marker_uidx2_back++;
	      next_unset_unsafe_ck(marker_exclude, &marker_uidx2_back);
	    } else {
	      window_trail_ct++;
	    }
	    cur_marker_pos = marker_pos[marker_uidx1_tmp];
	    if (cur_marker_pos > window_bp) {
	      cur_marker_pos -= window_bp;
	      while (marker_pos[marker_uidx2_back] < cur_marker_pos) {
		window_trail_ct--;
		marker_uidx2_back++;
		next_unset_unsafe_ck(marker_exclude, &marker_uidx2_back);
	      }
	    }
	    if (marker_cms) {
	      cur_marker_cm = marker_cms[marker_uidx1_tmp] - window_cm;
	      while (marker_cms[marker_uidx2_back] < cur_marker_cm) {
		window_trail_ct--;
		marker_uidx2_back++;
		next_unset_unsafe_ck(marker_exclude, &marker_uidx2_back);
	      }
	    }
	  }
	  if (marker_uidx2_fwd < chrom_last) {
	    cur_marker_pos = marker_pos[marker_uidx1_tmp] + window_bp;
	    if (!marker_cms) {
	      while (marker_pos[marker_uidx2_fwd2] <= cur_marker_pos) {
		marker_uidx2_fwd = marker_uidx2_fwd2;
		window_lead_ct++;
		if (marker_uidx2_fwd == chrom_last) {
		  break;
		}
		marker_uidx2_fwd2++;
		next_unset_unsafe_ck(marker_exclude, &marker_uidx2_fwd2);
		if (window_lead_ct > window_size_m1) {
		  break;
		}
	      }
	    } else {
	      cur_marker_cm = marker_cms[marker_uidx1_tmp] + window_cm;
	      while ((marker_pos[marker_uidx2_fwd2] <= cur_marker_pos) && (marker_cms[marker_uidx2_fwd2] <= window_cm)) {
		marker_uidx2_fwd = marker_uidx2_fwd2;
		window_lead_ct++;
		if (marker_uidx2_fwd == chrom_last) {
		  break;
		}
		marker_uidx2_fwd2++;
		next_unset_unsafe_ck(marker_exclude, &marker_uidx2_fwd2);
		if (window_lead_ct > window_size_m1) {
		  break;
		}
	      }
	    }
	  }
	  window_lead_ct--;
	}
      }
      if (!is_inter_chr) {
	if (idx1_subset) {
	  g_ld_interval1[block_idx1 * 2] = marker_idx2 - window_trail_ct - marker_idx2_base;
	} else {
	  g_ld_interval1[block_idx1 * 2] = marker_idx2 + 1 - marker_idx2_base;
	}
        g_ld_interval1[block_idx1 * 2 + 1] = marker_idx2 + window_lead_ct + 1 - marker_idx2_base;
      } else {
	if (!idx1_subset) {
          g_ld_interval1[block_idx1 * 2] = marker_idx2 + 1 - marker_idx2_base;
	} else {
	  g_ld_interval1[block_idx1 * 2] = 0;
	}
	g_ld_interval1[block_idx1 * 2 + 1] = marker_ct - marker_idx2_base;
      }

      if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, marker_uidx1_tmp), bedfile, loadbuf_raw, loadbuf)) {
	goto ld_report_dprime_ret_READ_FAIL;
      }
      if (is_haploid && hh_exists) {
        haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)loadbuf);
      }
      load_and_split3(nullptr, loadbuf, founder_ct, &(g_ld_geno1[block_idx1 * founder_ctsplit]), dummy_nm, dummy_nm, founder_ctv3, 0, 0, 1, &ulii);
      if (ulii == 3) {
        SET_BIT(block_idx1, g_epi_zmiss1);
      }
    }
    marker_uidx2 = marker_uidx2_base;
    if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto ld_report_dprime_ret_READ_FAIL;
    }

    cur_idx2_block_size = idx2_block_size;
    uljj = g_ld_interval1[2 * idx1_block_size - 1];
    marker_idx2_end = uljj + marker_idx2_base;
    marker_idx2_maxw = round_up_pow2(uljj, 4);
    if (marker_idx2_maxw > orig_marker_ctm8) {
      marker_idx2_maxw = orig_marker_ctm8;
    }
    g_ld_marker_ctm8 = marker_idx2_maxw;
    marker_idx2 = marker_idx2_base;
    do {
      if (cur_idx2_block_size > marker_idx2_end - marker_idx2) {
	cur_idx2_block_size = marker_idx2_end - marker_idx2;
      }
      if ((marker_idx2 < xend) && (marker_idx2 + cur_idx2_block_size > xstart)) {
	uii = MAXV(marker_idx2, xstart);
	g_ld_xstart2 = uii - marker_idx2;
	g_ld_xend2 = MINV(xend, marker_idx2 + cur_idx2_block_size) - uii;
      }
      fill_ulong_zero(BITCT_TO_WORDCT(cur_idx2_block_size), g_epi_zmiss2);
      for (block_idx2 = 0; block_idx2 < cur_idx2_block_size; marker_uidx2++, block_idx2++) {
	if (IS_SET(marker_exclude, marker_uidx2)) {
          marker_uidx2 = next_unset_ul_unsafe(marker_exclude, marker_uidx2);
          if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	    goto ld_report_dprime_ret_READ_FAIL;
	  }
	}
	if (marker_uidx2 >= chrom_end) {
	  chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx2);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
	  is_x = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[X_OFFSET]);
	  is_y = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
	}
	if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, marker_uidx2), bedfile, loadbuf_raw, loadbuf)) {
	  goto ld_report_dprime_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)loadbuf);
	}
	ulptr = &(g_ld_geno2[block_idx2 * founder_ctsplit]);
	load_and_split3(nullptr, loadbuf, founder_ct, ulptr, dummy_nm, dummy_nm, founder_ctv3, 0, 0, 1, &ulii);
	uiptr = &(g_epi_tot2[block_idx2 * 3]);
	uiptr[0] = popcount_longs(ulptr, founder_ctv3);
	uiptr[1] = popcount_longs(&(ulptr[founder_ctv3]), founder_ctv3);
        uiptr[2] = popcount_longs(&(ulptr[2 * founder_ctv3]), founder_ctv3);
	if (ulii == 3) {
	  SET_BIT(block_idx2, g_epi_zmiss2);
	}
      }
      g_ld_idx2_block_size = cur_idx2_block_size;
      g_ld_idx2_block_start = marker_idx2 - marker_idx2_base;
      marker_idx2 += cur_idx2_block_size;
      is_last_block = (marker_idx2 >= marker_idx2_end);
      if (spawn_threads2(threads, &ld_dprime_thread, thread_ct, is_last_block)) {
	goto ld_report_dprime_ret_THREAD_CREATE_FAIL;
      }
      ld_dprime_thread((void*)0);
      join_threads2(threads, thread_ct, is_last_block);
    } while (!is_last_block);

    fputs("\b\b\b\b\b\b\b\b\b\b\bwriting]   \b\b\b", stdout);
    fflush(stdout);
    g_ld_marker_uidx1 = marker_uidx1;
    g_ld_block_idx1 = 0;
    g_ld_uidx2_start = marker_uidx2_base;
    g_ld_idx2_block_start = 0;
    g_ld_block_idx2 = 0;
    if (output_gz) {
      parallel_compress(outname, overflow_buf, not_first_write, ld_regular_emitn);
    } else {
      write_uncompressed(outname, overflow_buf, not_first_write, ld_regular_emitn);
    }
    not_first_write = 1;
    g_ld_is_first_block = 0;
  ld_report_dprime_done:
    marker_idx1 += idx1_block_size;
    fputs("\b\b\b\b\b\b\b\b\b\b          \b\b\b\b\b\b\b\b\b\b", stdout);
    if (marker_idx1 >= pct_thresh) {
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      pct = ((marker_idx1 - marker_idx1_start) * 100LLU) / job_size;
      if (pct < 100) {
	printf("\b\b%" PRIuPTR "%%", pct);
	fflush(stdout);
	pct_thresh = marker_idx1_start + ((++pct) * ((uint64_t)job_size)) / 100;
      }
    }
    if (marker_idx1 == marker_idx1_end) {
      break;
    }
    marker_uidx1 = jump_forward_unset_unsafe(marker_exclude_idx1, marker_uidx1 + 1, idx1_block_size);
  }
  fputs("\b\b\b", stdout);
  logprint(" done.\n");
  LOGPRINTFWW("Results written to %s .\n", outname);

  while (0) {
  ld_report_dprime_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_report_dprime_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_report_dprime_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
  return retval;
}

int32_t ld_report_regular(pthread_t* threads, Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, uintptr_t* sex_male, uintptr_t* founder_include2, uintptr_t* founder_male_include2, uintptr_t* loadbuf, char* outname, uint32_t hh_exists) {
  FILE* infile = nullptr;
  uintptr_t* marker_exclude = g_ld_marker_exclude;
  char* marker_ids = g_ld_marker_ids;
  uintptr_t max_marker_id_len = g_ld_max_marker_id_len;
  uint32_t ld_modifier = ldip->modifier;
  uint32_t output_gz = ld_modifier & LD_REPORT_GZ;
  uint32_t ignore_x = (ld_modifier & LD_IGNORE_X) & 1;
  uint32_t is_inter_chr = ld_modifier & LD_INTER_CHR;
  uint32_t snp_list_file = ld_modifier & LD_SNP_LIST_FILE;
  uintptr_t marker_ct = g_ld_marker_ct;
  uintptr_t marker_ct1 = marker_ct;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t founder_ct = g_ld_founder_ct;
  uintptr_t founder_ctl = BITCT_TO_WORDCT(founder_ct);
  uintptr_t founder_ct_192_long = g_ld_founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2) + g_ld_founder_ct_mld_rem * (192 / BITCT2);
  uintptr_t final_mask = get_final_mask(founder_ct);
  uintptr_t pct = 1;
  uintptr_t marker_idx2_maxw = 1;
  Chrom_info* chrom_info_ptr = g_ld_chrom_info_ptr;
  uintptr_t* marker_exclude_idx1 = marker_exclude;
  uint32_t* marker_pos = g_ld_marker_pos;
  double* marker_cms = g_ld_marker_cms;
  uint32_t founder_trail_ct = founder_ct_192_long - founder_ctl * 2;
  uint32_t idx1_subset = (ldip->snpstr || ldip->snps_rl.name_ct);
  uint32_t window_size_m1 = ldip->window_size - 1;
  uint32_t window_bp = ldip->window_bp;
  double window_cm = ldip->window_cm;
  uint32_t thread_ct = g_ld_thread_ct;
  uint32_t chrom_fo_idx = 0;
  uint32_t chrom_fo_idx2 = 0;
  uint32_t is_haploid = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t not_first_write = 0;
  uint32_t marker_uidx2_back = 0;
  uint32_t marker_uidx2_fwd = 0;
  uint32_t marker_uidx2_fwd2 = 0;
  uint32_t window_trail_ct = 0;
  uint32_t window_lead_ct = 0;
  uint32_t chrom_last = 0;
  int32_t retval = 0;
  unsigned char* bigstack_mark2;
  unsigned char* overflow_buf;
  uint32_t* id_map;
  char* sorted_ids;
  char* bufptr;
  uintptr_t thread_workload;
  uintptr_t idx1_block_size;
  uintptr_t idx2_block_size;
  uintptr_t cur_idx2_block_size;
  uintptr_t orig_marker_ctm8;
  uintptr_t marker_idx1_start;
  uintptr_t marker_idx1;
  uintptr_t marker_idx1_end;
  uintptr_t marker_idx2;
  uintptr_t job_size;
  uintptr_t pct_thresh;
  uintptr_t marker_uidx1;
  uintptr_t marker_uidx1_tmp;
  uintptr_t marker_uidx2_base;
  uintptr_t marker_uidx2;
  uintptr_t marker_idx2_base;
  uintptr_t marker_idx2_end;
  uintptr_t block_idx1;
  uintptr_t block_idx2;
  uintptr_t snplist_ct;
  uintptr_t max_snplist_id_len;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t window_size_ceil;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t chrom_idx2;
  uint32_t chrom_end2;
  uint32_t cur_marker_pos;
  double cur_marker_cm;
  uint32_t is_last_block;
  uint32_t uii;
  int32_t ii;
  if (bigstack_alloc_uc(262144, &overflow_buf)) {
    goto ld_report_regular_ret_NOMEM;
  }
  if (idx1_subset) {
    if (bigstack_alloc_ul(unfiltered_marker_ctl, &marker_exclude_idx1)) {
      goto ld_report_regular_ret_NOMEM;
    }
    fill_all_bits(unfiltered_marker_ct, marker_exclude_idx1);
    marker_uidx1 = next_unset_unsafe(marker_exclude, 0);
    if (ldip->snpstr && (!snp_list_file)) {
      bufptr = ldip->snpstr;
      uii = strlen(bufptr) + 1;
      if (uii > max_marker_id_len) {
	goto ld_report_regular_ret_EMPTY_SET1;
      }
      for (marker_idx1 = 0; marker_idx1 < marker_ct; marker_uidx1++, marker_idx1++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx1);
        if (!memcmp(&(marker_ids[marker_uidx1 * max_marker_id_len]), bufptr, uii)) {
	  break;
	}
      }
      if (marker_idx1 == marker_ct) {
	goto ld_report_regular_ret_EMPTY_SET1;
      }
      clear_bit_ul(marker_uidx1, marker_exclude_idx1);
      marker_ct1 = 1;
    } else {
      marker_ct1 = 0;
      retval = sort_item_ids(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref, &sorted_ids, &id_map);
      if (retval) {
	goto ld_report_regular_ret_1;
      }
      if (snp_list_file) {
        if (fopen_checked(ldip->snpstr, FOPEN_RB, &infile)) {
	  goto ld_report_regular_ret_OPEN_FAIL;
	}
	snplist_ct = 0;
	max_snplist_id_len = 0;
	retval = scan_token_ct_len(MAXLINELEN, infile, g_textbuf, &snplist_ct, &max_snplist_id_len);
	if (retval) {
	  goto ld_report_regular_ret_1;
	}
	if (!snplist_ct) {
	  goto ld_report_regular_ret_EMPTY_SET1;
	}
	if (bigstack_alloc_c(snplist_ct * max_snplist_id_len, &bufptr)) {
	  goto ld_report_regular_ret_NOMEM;
	}
	rewind(infile);
	retval = read_tokens(MAXLINELEN, snplist_ct, max_snplist_id_len, infile, g_textbuf, bufptr);
	if (retval) {
	  goto ld_report_regular_ret_1;
	}
        if (fclose_null(&infile)) {
          goto ld_report_regular_ret_READ_FAIL;
	}
	for (marker_idx1 = 0; marker_idx1 < snplist_ct; marker_idx1++) {
          ii = bsearch_str_nl(&(bufptr[marker_idx1 * max_snplist_id_len]), sorted_ids, max_marker_id_len, marker_ct);
          if (ii != -1) {
            uii = id_map[(uint32_t)ii];
            if (!is_set(marker_exclude_idx1, uii)) {
	      logerrprint("Error: Duplicate variant ID in --ld-snp-list file.\n");
	      goto ld_report_regular_ret_INVALID_FORMAT;
	    }
            clear_bit(uii, marker_exclude_idx1);
            marker_ct1++;
	  }
	}
      } else {
        retval = string_range_list_to_bitarr2(sorted_ids, id_map, marker_ct, max_marker_id_len, &(ldip->snps_rl), "ld-snps", marker_exclude_idx1);
	if (retval) {
	  goto ld_report_regular_ret_1;
	}
        bitvec_or(marker_exclude, unfiltered_marker_ctl, marker_exclude_idx1);
	// bugfix, 13 Jan 2017
	// another bugfix, 28 Mar 2017: popcounted the wrong array...
        marker_ct1 = unfiltered_marker_ct - popcount_longs(marker_exclude_idx1, unfiltered_marker_ctl);
      }
      if (!marker_ct1) {
	goto ld_report_regular_ret_EMPTY_SET1;
      }
      bigstack_reset(id_map);
    }
  }
  if ((parallel_tot > 1) && (marker_ct1 < 2 * parallel_tot)) {
    LOGERRPRINTF("Error: Too few variants in --r%s run for --parallel %u %u.\n", g_ld_is_r2? "2" : "", parallel_idx + 1, parallel_tot);
    goto ld_report_regular_ret_INVALID_CMDLINE;
  }
  // yeah, this is uneven in the inter-chr case
  marker_idx1_start = (((uint64_t)parallel_idx) * marker_ct1) / parallel_tot;
  marker_idx1 = marker_idx1_start;
  marker_idx1_end = (((uint64_t)(parallel_idx + 1)) * marker_ct1) / parallel_tot;
  job_size = marker_idx1_end - marker_idx1_start;
  pct_thresh = job_size / 100;

  if (is_inter_chr) {
    marker_idx2_maxw = marker_ct + idx1_subset - 1;
  } else {
    window_size_ceil = (idx1_subset + 1) * (window_size_m1 + 1) - 1;
    if ((window_size_m1 < 12) || ((!idx1_subset) && (window_size_m1 <= 16))) {
      marker_idx2_maxw = window_size_ceil;
    } else {
      for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
        marker_idx2_maxw = chrom_window_max(marker_pos, marker_exclude, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], window_size_ceil, window_bp * (idx1_subset + 1), marker_idx2_maxw);
      }
    }
  }

  g_ld_marker_exclude_idx1 = marker_exclude_idx1;
  g_ld_marker_exclude = marker_exclude;
  g_ld_is_inter_chr = is_inter_chr;

  g_ld_is_first_block = (!parallel_idx);
  if (g_ld_is_r2) {
    g_ld_window_r2 = ldip->window_r2;
  } else {
    g_ld_window_r2 = sqrt(ldip->window_r2);
  }
  if (ld_modifier & LD_DX) {
    // this is more like --fast-epistasis under the hood, since it requires the
    // entire 3x3 table
    g_ld_marker_ctm8 = round_up_pow2(marker_idx2_maxw, 4);
    retval = ld_report_dprime(threads, ldip, bedfile, bed_offset, marker_reverse, unfiltered_sample_ct, founder_info, sex_male, founder_include2, founder_male_include2, loadbuf, outname, hh_exists, marker_idx1_start, marker_idx1_end);
    goto ld_report_regular_ret_1;
  }
  marker_idx2_maxw = round_up_pow2(marker_idx2_maxw, 8);
  orig_marker_ctm8 = marker_idx2_maxw;
  g_ld_marker_ctm8 = marker_idx2_maxw;
  g_ld_keep_sign = 1;
  // each marker costs
  //   founder_ct_192_long * sizeof(intptr_t) for genotype buffer
  // + founder_ct_192_long * sizeof(intptr_t) for missing mask buffer
  // + sizeof(int32_t) for g_ld_missing_cts1 entry
  // + 2 * sizeof(int32_t) for window offset and size
  // + marker_idx2_maxw * sizeof(double) for g_ld_results buffer
  // round down to multiple of thread_ct for better workload distribution
  ulii = founder_ct_192_long * 2 * sizeof(intptr_t) + 3 * sizeof(int32_t) + marker_idx2_maxw * sizeof(double);
  idx1_block_size = bigstack_left() / (ulii * 2);
  thread_workload = idx1_block_size / thread_ct;
  if (!thread_workload) {
    goto ld_report_regular_ret_NOMEM;
  }
  idx1_block_size = thread_workload * thread_ct;
  if (idx1_block_size > job_size) {
    idx1_block_size = job_size;
  }
  bigstack_alloc_ul(founder_ct_192_long * idx1_block_size, &g_ld_geno1);
  bigstack_alloc_ul(founder_ct_192_long * idx1_block_size, &g_ld_geno_masks1);
  bigstack_alloc_ui(idx1_block_size, &g_ld_missing_cts1);
  bigstack_alloc_ui(idx1_block_size * 2, &g_ld_interval1);
  if (bigstack_alloc_d(marker_idx2_maxw * idx1_block_size, &g_ld_results)) {
    goto ld_report_regular_ret_NOMEM;
  }

  ulii -= 2 * sizeof(int32_t) + marker_idx2_maxw * sizeof(double);
  idx2_block_size = (bigstack_left() / ulii) & (~(7 * ONELU));
  if (idx2_block_size > marker_ct) {
    idx2_block_size = round_up_pow2(marker_ct, 8);
  }
  bigstack_mark2 = g_bigstack_base;
  while (1) {
    if (!idx2_block_size) {
      goto ld_report_regular_ret_NOMEM;
    }
    if (!(bigstack_alloc_ul(founder_ct_192_long * idx2_block_size, &g_ld_geno2) ||
          bigstack_alloc_ul(founder_ct_192_long * idx2_block_size, &g_ld_geno_masks2) ||
          bigstack_alloc_ui(idx2_block_size, &g_ld_missing_cts2))) {
      break;
    }
    bigstack_reset(bigstack_mark2);
    idx2_block_size -= 8;
  }
  uljj = founder_trail_ct + 2;
  for (ulii = 1; ulii <= idx1_block_size; ulii++) {
    fill_ulong_zero(uljj, &(g_ld_geno1[ulii * founder_ct_192_long - uljj]));
    fill_ulong_zero(uljj, &(g_ld_geno_masks1[ulii * founder_ct_192_long - uljj]));
  }
  for (ulii = 1; ulii <= idx2_block_size; ulii++) {
    fill_ulong_zero(uljj, &(g_ld_geno2[ulii * founder_ct_192_long - uljj]));
    fill_ulong_zero(uljj, &(g_ld_geno_masks2[ulii * founder_ct_192_long - uljj]));
  }
  marker_uidx1 = next_unset_unsafe(marker_exclude_idx1, 0);
  if (marker_idx1) {
    marker_uidx1 = jump_forward_unset_unsafe(marker_exclude_idx1, marker_uidx1 + 1, marker_idx1);
  }
  sprintf(g_logbuf, "--r%s%s%s%s to %s ... ", g_ld_is_r2? "2" : "", is_inter_chr? " inter-chr" : "", g_ld_marker_allele_ptrs? " in-phase" : "", g_ld_set_allele_freqs? " with-freqs" : "", outname);
  wordwrapb(16); // strlen("99% [processing]")
  logprintb();
  fputs("0%", stdout);
  while (1) {
    fputs(" [processing]", stdout);
    fflush(stdout);
    if (idx1_block_size > marker_idx1_end - marker_idx1) {
      idx1_block_size = marker_idx1_end - marker_idx1;
      if (idx1_block_size < thread_ct) {
        thread_ct = idx1_block_size;
        g_ld_thread_ct = thread_ct;
      }
    }
    g_ld_idx1_block_size = idx1_block_size;
    marker_uidx1_tmp = marker_uidx1;

    if (idx1_subset) {
      if (!is_inter_chr) {
	chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx1_tmp);
	marker_uidx2_base = window_back(marker_pos, marker_cms, marker_exclude, next_unset_unsafe(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx]), marker_uidx1, window_size_m1, window_bp, window_cm, &uii);
	marker_idx2_base = marker_uidx2_base - popcount_bit_idx(marker_exclude, 0, marker_uidx2_base);
	marker_idx2 = marker_idx2_base + uii;
      } else {
	marker_uidx2_base = next_unset_unsafe(marker_exclude, 0);
	marker_idx2_base = 0;
	marker_idx2 = 0; // ignored
      }
    } else {
      marker_idx2_base = marker_uidx1 + 1 - popcount_bit_idx(marker_exclude, 0, marker_uidx1);
      if (marker_idx2_base == marker_ct) {
	goto ld_report_regular_done;
      }
      marker_idx2 = marker_idx2_base - 1;
      marker_uidx2_base = next_unset_unsafe(marker_exclude, marker_uidx1 + 1);
    }
    if (fseeko(bedfile, bed_offset + (marker_uidx1 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto ld_report_regular_ret_READ_FAIL;
    }
    chrom_end = 0;
    for (block_idx1 = 0; block_idx1 < idx1_block_size; marker_uidx1_tmp++, block_idx1++, marker_idx2++) {
      if (IS_SET(marker_exclude_idx1, marker_uidx1_tmp)) {
        ulii = next_unset_ul_unsafe(marker_exclude_idx1, marker_uidx1_tmp);
	uljj = ulii - marker_uidx1_tmp - popcount_bit_idx(marker_exclude, marker_uidx1_tmp, ulii);
	if (uljj) {
	  uii = 1; // recalculate window beginning/end from scratch
	  marker_idx2 += uljj;
	}
	marker_uidx1_tmp = ulii;
        if (fseeko(bedfile, bed_offset + (marker_uidx1_tmp * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
          goto ld_report_regular_ret_READ_FAIL;
	}
      }
      if (marker_uidx1_tmp >= chrom_end) {
        chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx1_tmp);
        chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
        chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
	chrom_last = prev_unset_unsafe(marker_exclude, chrom_end);
        is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
	is_x = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[X_OFFSET]);
	is_y = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
	uii = 1;
      }
      if (!is_inter_chr) {
	// uii == 0 if we can perform an incremental update, 1 if we need
	// fully-powered window_back()/window_forward()
	if (uii) {
	  if (idx1_subset) {
	    marker_uidx2_back = window_back(marker_pos, marker_cms, marker_exclude, next_unset_unsafe(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx]), marker_uidx1_tmp, window_size_m1, window_bp, window_cm, &window_trail_ct);
	  }
	  marker_uidx2_fwd = window_forward(marker_pos, marker_cms, marker_exclude, marker_uidx1_tmp, chrom_last, window_size_m1, window_bp, window_cm, &window_lead_ct);
	  marker_uidx2_fwd2 = marker_uidx2_fwd;
	  if (marker_uidx2_fwd < chrom_last) {
	    marker_uidx2_fwd2++;
	    next_unset_unsafe_ck(marker_exclude, &marker_uidx2_fwd2);
	  }
	  uii = 0;
	} else {
	  if (idx1_subset) {
	    if (window_trail_ct == window_size_m1) {
	      marker_uidx2_back++;
	      next_unset_unsafe_ck(marker_exclude, &marker_uidx2_back);
	    } else {
	      window_trail_ct++;
	    }
	    cur_marker_pos = marker_pos[marker_uidx1_tmp];
	    if (cur_marker_pos > window_bp) {
	      cur_marker_pos -= window_bp;
	      while (marker_pos[marker_uidx2_back] < cur_marker_pos) {
		window_trail_ct--;
		marker_uidx2_back++;
		next_unset_unsafe_ck(marker_exclude, &marker_uidx2_back);
	      }
	    }
	    if (marker_cms) {
	      cur_marker_cm = marker_cms[marker_uidx1_tmp] - window_cm;
	      while (marker_cms[marker_uidx2_back] < cur_marker_cm) {
		window_trail_ct--;
		marker_uidx2_back++;
		next_unset_unsafe_ck(marker_exclude, &marker_uidx2_back);
	      }
	    }
	  }
	  if (marker_uidx2_fwd < chrom_last) {
	    cur_marker_pos = marker_pos[marker_uidx1_tmp] + window_bp;
	    if (!marker_cms) {
	      while (marker_pos[marker_uidx2_fwd2] <= cur_marker_pos) {
		marker_uidx2_fwd = marker_uidx2_fwd2;
		window_lead_ct++;
		if (marker_uidx2_fwd == chrom_last) {
		  break;
		}
		marker_uidx2_fwd2++;
		next_unset_unsafe_ck(marker_exclude, &marker_uidx2_fwd2);
		if (window_lead_ct > window_size_m1) {
		  break;
		}
	      }
	    } else {
	      cur_marker_cm = marker_cms[marker_uidx1_tmp] + window_cm;
	      while ((marker_pos[marker_uidx2_fwd2] <= cur_marker_pos) && (marker_cms[marker_uidx2_fwd2] <= cur_marker_cm)) {
		marker_uidx2_fwd = marker_uidx2_fwd2;
		window_lead_ct++;
		if (marker_uidx2_fwd == chrom_last) {
		  break;
		}
		marker_uidx2_fwd2++;
		next_unset_unsafe_ck(marker_exclude, &marker_uidx2_fwd2);
		if (window_lead_ct > window_size_m1) {
		  break;
		}
	      }
	    }
	  }
	  window_lead_ct--;
	}
      }
      if (!is_inter_chr) {
	if (idx1_subset) {
	  g_ld_interval1[block_idx1 * 2] = marker_idx2 - window_trail_ct - marker_idx2_base;
	} else {
	  g_ld_interval1[block_idx1 * 2] = marker_idx2 + 1 - marker_idx2_base;
	}
        g_ld_interval1[block_idx1 * 2 + 1] = marker_idx2 + window_lead_ct + 1 - marker_idx2_base;
      } else {
	if (!idx1_subset) {
          g_ld_interval1[block_idx1 * 2] = marker_idx2 + 1 - marker_idx2_base;
	} else {
	  g_ld_interval1[block_idx1 * 2] = 0;
	}
	g_ld_interval1[block_idx1 * 2 + 1] = marker_ct - marker_idx2_base;
      }

      if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, marker_uidx1_tmp), bedfile, loadbuf, &(g_ld_geno1[block_idx1 * founder_ct_192_long]))) {
	goto ld_report_regular_ret_READ_FAIL;
      }
      if (is_haploid && hh_exists) {
	haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(g_ld_geno1[block_idx1 * founder_ct_192_long])));
      }
      ld_process_load2(&(g_ld_geno1[block_idx1 * founder_ct_192_long]), &(g_ld_geno_masks1[block_idx1 * founder_ct_192_long]), &(g_ld_missing_cts1[block_idx1]), founder_ct, is_x && (!ignore_x), founder_male_include2);
    }
    marker_uidx2 = marker_uidx2_base;
    if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto ld_report_regular_ret_READ_FAIL;
    }

    cur_idx2_block_size = idx2_block_size;
    uljj = g_ld_interval1[2 * idx1_block_size - 1];
    marker_idx2_end = uljj + marker_idx2_base;
    marker_idx2_maxw = round_up_pow2(uljj, 8);
    if (marker_idx2_maxw > orig_marker_ctm8) {
      marker_idx2_maxw = orig_marker_ctm8;
    }
    g_ld_marker_ctm8 = marker_idx2_maxw;
    marker_idx2 = marker_idx2_base;
    chrom_end2 = 0;
    do {
      if (cur_idx2_block_size > marker_idx2_end - marker_idx2) {
	cur_idx2_block_size = marker_idx2_end - marker_idx2;
      }

      for (block_idx2 = 0; block_idx2 < cur_idx2_block_size; marker_uidx2++, block_idx2++) {
	// todo: when set has big holes in the middle, do not load everything
	if (IS_SET(marker_exclude, marker_uidx2)) {
          marker_uidx2 = next_unset_ul_unsafe(marker_exclude, marker_uidx2);
          if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	    goto ld_report_regular_ret_READ_FAIL;
	  }
	}
	if (marker_uidx2 >= chrom_end2) {
	  chrom_fo_idx2 = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx2);
	  chrom_idx2 = chrom_info_ptr->chrom_file_order[chrom_fo_idx2];
	  chrom_end2 = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx2 + 1];
	  is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx2);
	  is_x = (((int32_t)chrom_idx2) == chrom_info_ptr->xymt_codes[X_OFFSET]);
	  is_y = (((int32_t)chrom_idx2) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
	}
	if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, marker_uidx2), bedfile, loadbuf, &(g_ld_geno2[block_idx2 * founder_ct_192_long]))) {
	  goto ld_report_regular_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(g_ld_geno2[block_idx2 * founder_ct_192_long])));
	}
	ld_process_load2(&(g_ld_geno2[block_idx2 * founder_ct_192_long]), &(g_ld_geno_masks2[block_idx2 * founder_ct_192_long]), &(g_ld_missing_cts2[block_idx2]), founder_ct, is_x && (!ignore_x), founder_male_include2);
      }

      g_ld_idx2_block_size = cur_idx2_block_size;
      g_ld_idx2_block_start = marker_idx2 - marker_idx2_base;
      marker_idx2 += cur_idx2_block_size;
      is_last_block = (marker_idx2 >= marker_idx2_end);
      if (spawn_threads2(threads, &ld_block_thread, thread_ct, is_last_block)) {
	goto ld_report_regular_ret_THREAD_CREATE_FAIL;
      }
      ld_block_thread((void*)0);
      join_threads2(threads, thread_ct, is_last_block);
    } while (!is_last_block);

    fputs("\b\b\b\b\b\b\b\b\b\b\bwriting]   \b\b\b", stdout);
    fflush(stdout);
    g_ld_marker_uidx1 = marker_uidx1;
    g_ld_block_idx1 = 0;
    g_ld_uidx2_start = marker_uidx2_base;
    g_ld_idx2_block_start = 0;
    g_ld_block_idx2 = 0;
    if (output_gz) {
      parallel_compress(outname, overflow_buf, not_first_write, ld_regular_emitn);
    } else {
      write_uncompressed(outname, overflow_buf, not_first_write, ld_regular_emitn);
    }
    not_first_write = 1;
    g_ld_is_first_block = 0;
  ld_report_regular_done:
    marker_idx1 += idx1_block_size;
    fputs("\b\b\b\b\b\b\b\b\b\b          \b\b\b\b\b\b\b\b\b\b", stdout);
    if (marker_idx1 >= pct_thresh) {
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      pct = ((marker_idx1 - marker_idx1_start) * 100LLU) / job_size;
      if (pct < 100) {
	printf("\b\b%" PRIuPTR "%%", pct);
	fflush(stdout);
	pct_thresh = marker_idx1_start + ((++pct) * ((uint64_t)job_size)) / 100;
      }
    }
    if (marker_idx1 == marker_idx1_end) {
      break;
    }
    marker_uidx1 = jump_forward_unset_unsafe(marker_exclude_idx1, marker_uidx1 + 1, idx1_block_size);
  }
  fputs("\b\b", stdout);
  logprint("done.\n");
  while (0) {
  ld_report_regular_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_report_regular_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ld_report_regular_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_report_regular_ret_EMPTY_SET1:
    logerrprint("Error: No valid variants specified by --ld-snp/--ld-snps/--ld-snp-list.\n");
  ld_report_regular_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  ld_report_regular_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  ld_report_regular_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 ld_report_regular_ret_1:
  fclose_cond(infile);
  // trust parent to free memory
  return retval;
}

int32_t ld_report(pthread_t* threads, Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, double* set_allele_freqs, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, double* marker_cms, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct);
  uintptr_t founder_ct = popcount_longs(founder_info, unfiltered_sample_ctv2 / 2);
  uintptr_t* founder_include2 = nullptr;
  uintptr_t* founder_male_include2 = nullptr;
  uintptr_t founder_ct_mld = (founder_ct + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
  uint32_t founder_ct_mld_m1 = ((uint32_t)founder_ct_mld) - 1;
#ifdef __LP64__
  uint32_t founder_ct_mld_rem = (MULTIPLEX_LD / 192) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 192;
#else
  uint32_t founder_ct_mld_rem = (MULTIPLEX_LD / 48) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 48;
#endif
  uintptr_t founder_ct_192_long = founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2) + founder_ct_mld_rem * (192 / BITCT2);
  uint32_t ld_modifier = ldip->modifier;
  uint32_t is_binary = ld_modifier & (LD_MATRIX_BIN | LD_MATRIX_BIN4);
  uint32_t output_gz = ld_modifier & LD_REPORT_GZ;
  char* bufptr = memcpyl3a(outname_end, ".ld");
  int32_t retval = 0;
  uintptr_t* loadbuf;

  g_ld_modifier = ld_modifier;
  g_ld_founder_ct = founder_ct;
  g_ld_founder_ct_192_long = founder_ct_192_long;
  g_ld_founder_ct_mld_m1 = founder_ct_mld_m1;
  g_ld_founder_ct_mld_rem = founder_ct_mld_rem;
  g_ld_is_r2 = ld_modifier & LD_R2;
  g_ld_marker_ct = marker_ct;
  g_ld_chrom_info_ptr = chrom_info_ptr;
  g_ld_thread_ct = g_thread_ct;
  g_ld_set_allele_freqs = (ld_modifier & LD_WITH_FREQS)? set_allele_freqs : nullptr;
  if (founder_ct < 2) {
    LOGERRPRINTF("Warning: Skipping --r%s since there are less than two founders.\n(--make-founders may come in handy here.)\n", g_ld_is_r2? "2" : "");
    goto ld_report_ret_1;
  } else if (founder_ct >= 0x20000000) {
    logerrprint("Error: --r/--r2 does not support >= 2^29 samples.\n");
    goto ld_report_ret_INVALID_CMDLINE;
  }
  if ((marker_ct > 400000) && (!(ld_modifier & LD_YES_REALLY)) && (parallel_tot == 1) && ((ld_modifier & LD_MATRIX_SHAPEMASK) || ((ld_modifier & LD_INTER_CHR) && (!ldip->snpstr) && (!ldip->snps_rl.name_ct) && ((!g_ld_is_r2) || (ldip->window_r2 == 0.0))))) {
    logerrprint("Error: Gigantic (over 400k loci) --r/--r2 unfiltered, non-distributed\ncomputation.  Rerun with the 'yes-really' modifier if you are SURE you have\nenough hard drive space and want to do this.\n");
    goto ld_report_ret_INVALID_CMDLINE;
  }
  if (alloc_collapsed_haploid_filters(founder_info, sex_male, unfiltered_sample_ct, founder_ct, XMHH_EXISTS | hh_exists, 1, &founder_include2, &founder_male_include2)) {
    goto ld_report_ret_NOMEM;
  }
  if (bigstack_alloc_ul(unfiltered_sample_ctv2, &loadbuf)) {
    goto ld_report_ret_NOMEM;
  }
  loadbuf[unfiltered_sample_ctv2 - 2] = 0;
  loadbuf[unfiltered_sample_ctv2 - 1] = 0;
  // possible todo: throw out all monomorphic sites (and, in at least the
  // matrix case, dump a list of expelled site IDs)
  if (is_binary) {
    bufptr = memcpya(bufptr, ".bin", 4);
  }
  if (parallel_tot > 1) {
    *bufptr++ = '.';
    bufptr = uint32toa(parallel_idx + 1, bufptr);
  }
  if (!is_binary) {
    g_ld_delimiter = (ld_modifier & LD_MATRIX_SPACES)? ' ' : '\t';
    if (output_gz) {
      bufptr = memcpyl3a(bufptr, ".gz");
    }
  }
  *bufptr = '\0';
  if (ld_modifier & LD_INPHASE) {
    if (max_marker_allele_len * 4 + plink_maxsnp * 2 + get_max_chrom_slen(chrom_info_ptr) * 2 + 128 > MAXLINELEN) {
      logerrprint("Error: --r/--r2 in-phase does not support very long allele codes.\n");
      goto ld_report_ret_INVALID_CMDLINE;
    }
    g_ld_marker_allele_ptrs = marker_allele_ptrs;
  } else {
    g_ld_marker_allele_ptrs = nullptr;
  }
  if (ld_modifier & (LD_MATRIX_SQ | LD_MATRIX_SQ0 | LD_MATRIX_TRI)) {
    retval = ld_report_matrix(threads, ldip, bedfile, bed_offset, unfiltered_marker_ct, marker_exclude, marker_reverse, unfiltered_sample_ct, founder_info, parallel_idx, parallel_tot, sex_male, founder_include2, founder_male_include2, loadbuf, outname, hh_exists);
  } else {
    g_ld_plink_maxsnp = plink_maxsnp;
    g_ld_marker_ids = marker_ids;
    g_ld_marker_pos = marker_pos;
    g_ld_marker_cms = (ldip->window_cm == -1)? nullptr : marker_cms;
    g_ld_marker_exclude = marker_exclude;
    g_ld_max_marker_id_len = max_marker_id_len;
    retval = ld_report_regular(threads, ldip, bedfile, bed_offset, unfiltered_marker_ct, marker_reverse, unfiltered_sample_ct, founder_info, parallel_idx, parallel_tot, sex_male, founder_include2, founder_male_include2, loadbuf, outname, hh_exists);
  }
  while (0) {
  ld_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_report_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 ld_report_ret_1:
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t show_tags(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  // Similar to ld_prune() and flipscan().
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t founder_ct = popcount_longs(founder_info, unfiltered_sample_ctl);
  uintptr_t founder_ctl = BITCT_TO_WORDCT(founder_ct);
  uintptr_t final_mask = get_final_mask(founder_ct);
  uintptr_t marker_idx = 0;
  uintptr_t max_window_size = 1;
  uintptr_t pct = 1;
  uintptr_t pct_thresh = marker_ct / 100;
  FILE* infile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t* final_set = nullptr;
  uintptr_t* founder_include2 = nullptr;
  uintptr_t* founder_male_include2 = nullptr;
  char* chrom_name_ptr = nullptr;
  double tag_thresh = ldip->show_tags_r2 * (1 - SMALL_EPSILON);
  uint32_t tags_list = (ldip->modifier & LD_SHOW_TAGS_LIST_ALL) || (!ldip->show_tags_fname);
  uint32_t twocolumn = ldip->modifier & LD_SHOW_TAGS_MODE2;
  uint32_t ignore_x = (ldip->modifier & LD_IGNORE_X) & 1;
  uint32_t window_bp = ldip->show_tags_bp;
  uint32_t target_ct = 0;
  uint32_t chrom_name_len = 0;
  int32_t retval = 0;
  char chrom_name_buf[3 + MAX_CHROM_TEXTNUM_SLEN];
  int32_t dp_result[5];
  uintptr_t founder_ct_192_long;
  uintptr_t founder_ctwd12;
  uintptr_t founder_ctwd12_rem;
  uintptr_t lshift_last;
  uintptr_t line_idx;
  uintptr_t unrecog_ct;
  uintptr_t max_window_ctal;
  uintptr_t max_window_ctl;
  uintptr_t marker_uidx;
  uintptr_t window_cidx;
  uintptr_t window_cidx2;
  uintptr_t window_cidx3;
  uintptr_t marker_uidx2;
  uintptr_t ulii;
  uintptr_t* targets;
  uintptr_t* loadbuf_raw;
  uintptr_t* geno;
  uintptr_t* geno_masks;
  uintptr_t* geno_fixed_vec_ptr;
  uintptr_t* mask_fixed_vec_ptr;
  uintptr_t* geno_var_vec_ptr;
  uintptr_t* mask_var_vec_ptr;
  uintptr_t* cur_targets;
  uintptr_t* tag_matrix;
  uintptr_t* tag_matrix_row_ptr;
  char* sorted_marker_ids;
  char* bufptr;
  char* bufptr2;
  uint32_t* marker_id_map;
  uint32_t* window_uidxs;
  uint32_t* window_cidx_starts;
  uint32_t* missing_cts;
  double non_missing_ctd;
  double cov12;
  double dxx;
  double dyy;
  uint32_t founder_ct_mld_m1;
  uint32_t founder_ct_mld_rem;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t chrom_marker_ct;
  uint32_t chrom_marker_idx;
  uint32_t is_haploid;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_target;
  uint32_t marker_pos_thresh;
  uint32_t fixed_missing_ct;
  uint32_t fixed_non_missing_ct;
  uint32_t non_missing_ct;
  uint32_t slen;
  uint32_t tag_ct;
  uint32_t marker_uidx3;
  uint32_t min_bp;
  uint32_t max_bp;
  uint32_t cur_bp;
  uint32_t uii;
  int32_t ii;
  if (founder_ct < 2) {
    logerrprint("Warning: Skipping --show-tags since there are less than two founders.\n(--make-founders may come in handy here.)\n");
    goto show_tags_ret_1;
  }
  if (bigstack_alloc_ul(unfiltered_marker_ctl, &targets) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw)) {
    goto show_tags_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  if (ldip->show_tags_fname) {
    fill_ulong_zero(unfiltered_marker_ctl, targets);
    retval = sort_item_ids(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref, &sorted_marker_ids, &marker_id_map);
    if (retval) {
      goto show_tags_ret_1;
    }
    if (fopen_checked(ldip->show_tags_fname, "r", &infile)) {
      goto show_tags_ret_OPEN_FAIL;
    }
    g_textbuf[MAXLINELEN - 1] = ' ';
    line_idx = 0;
    unrecog_ct = 0;
    while (fgets(g_textbuf, MAXLINELEN, infile)) {
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 1]) {
	LOGERRPRINTF("Error: Line %" PRIuPTR " of --show-tags file is pathologically long.\n", line_idx);
	goto show_tags_ret_INVALID_FORMAT;
      }
      bufptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      slen = strlen_se(bufptr);
      if (twocolumn) {
	bufptr2 = skip_initial_spaces(&(bufptr[slen]));
        if (!bufptr2) {
	  LOGERRPRINTF("Error: Line %" PRIuPTR " of --show-tags file has fewer tokens than expected.\n", line_idx);
	  goto show_tags_ret_INVALID_FORMAT;
	}
        if ((*bufptr2 != '1') || (!is_space_or_eoln(bufptr2[1]))) {
	  continue;
	}
      }
      ii = bsearch_str(bufptr, slen, sorted_marker_ids, max_marker_id_len, marker_ct);
      if (ii == -1) {
	unrecog_ct++;
	continue;
      }
      marker_uidx = marker_id_map[(uint32_t)ii];
      if (IS_SET(targets, marker_uidx)) {
        bufptr[slen] = '\0';
        LOGERRPRINTF("Error: Duplicate variant ID '%s' in --show-tags file.\n", bufptr);
	goto show_tags_ret_INVALID_FORMAT;
      }
      SET_BIT(marker_uidx, targets);
    }
    if (fclose_null(&infile)) {
      goto show_tags_ret_READ_FAIL;
    }
    bigstack_reset((unsigned char*)marker_id_map);
    target_ct = popcount_longs(targets, unfiltered_marker_ctl);
    if (!target_ct) {
      logerrprint("Error: No recognized variant IDs in --show-tags file.\n");
      goto show_tags_ret_INVALID_FORMAT;
    }
    if (bigstack_alloc_ul(unfiltered_marker_ctl, &final_set)) {
      goto show_tags_ret_NOMEM;
    }
    memcpy(final_set, targets, unfiltered_marker_ctl * sizeof(intptr_t));
    LOGPRINTF("--show-tags: %u target variant%s loaded.\n", target_ct, (target_ct == 1)? "" : "s");
    if (unrecog_ct) {
      LOGERRPRINTF("Warning: %" PRIuPTR " unrecognized variant ID%s in --show-tags file.\n", unrecog_ct, (unrecog_ct == 1)? "" : "s");
    }
  } else {
    bitarr_invert_copy(marker_exclude, unfiltered_marker_ct, targets);
  }
  // force founder_male_include2 allocation
  if (alloc_collapsed_haploid_filters(founder_info, sex_male, unfiltered_sample_ct, founder_ct, XMHH_EXISTS | hh_exists, 1, &founder_include2, &founder_male_include2)) {
    goto show_tags_ret_NOMEM;
  }
  founder_ct_mld_m1 = (founder_ct - 1) / MULTIPLEX_LD;
  ulii = founder_ct_mld_m1 + 1;
#ifdef __LP64__
  founder_ct_mld_rem = (MULTIPLEX_LD / 192) - (ulii * MULTIPLEX_LD - founder_ct) / 192;
#else
  founder_ct_mld_rem = (MULTIPLEX_LD / 48) - (ulii * MULTIPLEX_LD - founder_ct) / 48;
#endif
  founder_ct_192_long = founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2) + founder_ct_mld_rem * (192 / BITCT2);
  uii = founder_ct / BITCT2;
  founder_ctwd12 = uii / 12;
  founder_ctwd12_rem = uii - (12 * founder_ctwd12);
  lshift_last = 2 * ((0x7fffffc0 - founder_ct) % BITCT2);
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    max_window_size = chrom_window_max(marker_pos, marker_exclude, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], 0x7fffffff, window_bp * 2, max_window_size);
  }
  max_window_ctl = BITCT_TO_WORDCT(max_window_size);
  max_window_ctal = max_window_ctl * BITCT;
  if (bigstack_alloc_ui(max_window_size, &window_uidxs) ||
      bigstack_alloc_ui(max_window_size, &window_cidx_starts) ||
      bigstack_alloc_ui(max_window_size, &missing_cts) ||
      bigstack_alloc_ul(max_window_size * founder_ct_192_long, &geno) ||
      bigstack_alloc_ul(max_window_size * founder_ct_192_long, &geno_masks) ||
      bigstack_alloc_ul(max_window_ctl, &cur_targets) ||
      bigstack_alloc_ul(max_window_size * max_window_ctl, &tag_matrix)) {
    goto show_tags_ret_NOMEM;
  }
  uii = 2 + founder_ct_192_long - founder_ctl * 2;
  for (ulii = 1; ulii <= max_window_size; ulii++) {
    fill_ulong_zero(uii, &(geno[ulii * founder_ct_192_long - uii]));
    fill_ulong_zero(uii, &(geno_masks[ulii * founder_ct_192_long - uii]));
  }

  if (tags_list) {
    memcpy(outname_end, ".tags.list", 11);
    if (fopen_checked(outname, "w", &outfile)) {
      goto show_tags_ret_WRITE_FAIL;
    }
    sprintf(g_textbuf, "%%%us  CHR         BP NTAG       LEFT      RIGHT   KBSPAN TAGS\n", plink_maxsnp);
    fprintf(outfile, g_textbuf, "SNP");
  }
  printf("--show-tags%s: 0%%", final_set? "" : " all");
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
    chrom_marker_ct = chrom_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, chrom_end);
    if (chrom_marker_ct < 2) {
      marker_idx += chrom_marker_ct;
      continue;
    }
    chrom_name_ptr = chrom_name_buf5w4write(chrom_info_ptr, chrom_idx, &chrom_name_len, chrom_name_buf);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    is_x = (chrom_idx == ((uint32_t)chrom_info_ptr->xymt_codes[X_OFFSET]));
    is_y = (chrom_idx == ((uint32_t)chrom_info_ptr->xymt_codes[Y_OFFSET]));
    if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto show_tags_ret_READ_FAIL;
    }
    chrom_marker_idx = 0;
    window_cidx = max_window_size - 1;
    window_cidx2 = 0;
    do {
      if (++window_cidx == max_window_size) {
        window_cidx = 0;
      }
      window_uidxs[window_cidx] = marker_uidx;
      is_target = IS_SET(targets, marker_uidx);
      if (is_target) {
	SET_BIT(window_cidx, cur_targets);
      } else {
	CLEAR_BIT(window_cidx, cur_targets);
      }

      // circular index of beginning of window starting at current marker
      window_cidx_starts[window_cidx] = window_cidx2;
      geno_fixed_vec_ptr = &(geno[window_cidx * founder_ct_192_long]);
      mask_fixed_vec_ptr = &(geno_masks[window_cidx * founder_ct_192_long]);
      fill_ulong_zero(max_window_ctl, &(tag_matrix[window_cidx * max_window_ctl]));
      if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, geno_fixed_vec_ptr)) {
        goto show_tags_ret_READ_FAIL;
      }
      if (is_haploid && hh_exists) {
	haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)loadbuf_raw);
      }
      ld_process_load2(geno_fixed_vec_ptr, mask_fixed_vec_ptr, &fixed_missing_ct, founder_ct, is_x && (!ignore_x), founder_male_include2);
      fixed_non_missing_ct = founder_ct - fixed_missing_ct;
      missing_cts[window_cidx] = fixed_missing_ct;
      window_cidx3 = window_cidx2;
      while (window_cidx3 != window_cidx) {
	if (is_target || IS_SET(cur_targets, window_cidx3)) {
	  // don't bother computing r^2 if no target variant involved
	  geno_var_vec_ptr = &(geno[window_cidx3 * founder_ct_192_long]);
	  mask_var_vec_ptr = &(geno_masks[window_cidx3 * founder_ct_192_long]);
	  non_missing_ct = fixed_non_missing_ct - missing_cts[window_cidx3];
	  if (fixed_missing_ct && missing_cts[window_cidx3]) {
	    non_missing_ct += ld_missing_ct_intersect(mask_var_vec_ptr, mask_fixed_vec_ptr, founder_ctwd12, founder_ctwd12_rem, lshift_last);
	  }
	  if (non_missing_ct) {
	    dp_result[0] = founder_ct;
	    dp_result[1] = -((int32_t)fixed_non_missing_ct);
	    dp_result[2] = missing_cts[window_cidx3] - founder_ct;
	    dp_result[3] = dp_result[1];
	    dp_result[4] = dp_result[2];
	    ld_dot_prod(geno_var_vec_ptr, geno_fixed_vec_ptr, mask_var_vec_ptr, mask_fixed_vec_ptr, dp_result, founder_ct_mld_m1, founder_ct_mld_rem);
	    non_missing_ctd = (double)((int32_t)non_missing_ct);
	    dxx = dp_result[1];
	    dyy = dp_result[2];
	    cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
	    dxx = (dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy);
	    if (cov12 * cov12 > dxx * tag_thresh) {
	      set_bit_ul(window_cidx * max_window_ctal + window_cidx3, tag_matrix);
	      set_bit_ul(window_cidx3 * max_window_ctal + window_cidx, tag_matrix);
	    }
	  }
	}
        if (++window_cidx3 == max_window_size) {
	  window_cidx3 = 0;
	}
      }
      if (++chrom_marker_idx < chrom_marker_ct) {
        marker_uidx++;
        if (IS_SET(marker_exclude, marker_uidx)) {
          marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
          if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	    goto show_tags_ret_READ_FAIL;
	  }
	}
        marker_pos_thresh = marker_pos[marker_uidx];
        if (marker_pos_thresh < window_bp) {
	  marker_pos_thresh = 0;
	} else {
	  marker_pos_thresh -= window_bp;
	}
      } else {
	// close out the chromosome
	marker_pos_thresh = 0x80000000U;
      }
      marker_uidx2 = window_uidxs[window_cidx2];
      while (marker_pos[marker_uidx2] < marker_pos_thresh) {
	if (IS_SET(cur_targets, window_cidx2)) {
	  // bugfix: tag_matrix_row_ptr is not always 16-byte aligned.
	  tag_ct = popcount_longs_nzbase(tag_matrix, window_cidx2 * max_window_ctl, (window_cidx2 + 1) * max_window_ctl);
	  tag_matrix_row_ptr = &(tag_matrix[window_cidx2 * max_window_ctl]);
	  min_bp = marker_pos[marker_uidx2];
	  max_bp = marker_pos[marker_uidx2];
	  window_cidx3 = window_cidx_starts[window_cidx2];
	  for (uii = 0; uii < tag_ct; uii++, window_cidx3++) {
	    next_set_ul_ck(tag_matrix_row_ptr, max_window_size, &window_cidx3);
	    if (window_cidx3 == max_window_size) {
	      window_cidx3 = next_set_unsafe(tag_matrix_row_ptr, 0);
	    }
	    marker_uidx3 = window_uidxs[window_cidx3];
	    if (final_set) {
	      SET_BIT(marker_uidx3, final_set);
	    }
	    if (tags_list) {
	      cur_bp = marker_pos[marker_uidx3];
	      if (cur_bp < min_bp) {
		min_bp = cur_bp;
	      } else if (cur_bp > max_bp) {
		max_bp = cur_bp;
	      }
	    }
	  }
	  if (tags_list) {
	    bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), g_textbuf);
	    *bufptr++ = ' ';
	    bufptr = memcpyax(bufptr, chrom_name_ptr, chrom_name_len, ' ');
	    bufptr = uint32toa_w10x(marker_pos[marker_uidx2], ' ', bufptr);
	    bufptr = uint32toa_w4x(tag_ct, ' ', bufptr);
	    bufptr = uint32toa_w10x(min_bp, ' ', bufptr);
	    bufptr = uint32toa_w10x(max_bp, ' ', bufptr);
	    bufptr = width_force(8, bufptr, dtoa_g(((int32_t)(max_bp - min_bp + 1)) * 0.001, bufptr));
	    *bufptr++ = ' ';
	    if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	      goto show_tags_ret_WRITE_FAIL;
	    }
	    window_cidx3 = window_cidx_starts[window_cidx2];
	    for (uii = 0; uii < tag_ct; uii++, window_cidx3++) {
	      next_set_ul_ck(tag_matrix_row_ptr, max_window_size, &window_cidx3);
	      if (window_cidx3 == max_window_size) {
		window_cidx3 = next_set_unsafe(tag_matrix_row_ptr, 0);
	      }
	      if (uii) {
		putc_unlocked('|', outfile);
	      }
	      fputs(&(marker_ids[window_uidxs[window_cidx3] * max_marker_id_len]), outfile);
	    }
	    if (!tag_ct) {
	      fputs("NONE", outfile);
	    }
	    putc_unlocked('\n', outfile);
	  }
	}
	if (++marker_idx >= pct_thresh) {
	  if (pct > 10) {
	    putc_unlocked('\b', stdout);
	  }
          pct = (marker_idx * 100LLU) / marker_ct;
          if (pct < 100) {
            printf("\b\b%" PRIuPTR "%%", pct);
            fflush(stdout);
            pct_thresh = ((++pct) * ((uint64_t)marker_ct)) / 100;
	  }
	}
	if (window_cidx2 == window_cidx) {
	  if (++window_cidx2 == max_window_size) {
	    window_cidx2 = 0;
	  }
	  break;
	}
	if (++window_cidx2 == max_window_size) {
	  window_cidx2 = 0;
	}
	marker_uidx2 = window_uidxs[window_cidx2];
      }
    } while (chrom_marker_idx < chrom_marker_ct);
  }
  putc_unlocked('\r', stdout);
  if (tags_list) {
    if (fclose_null(&outfile)) {
      goto show_tags_ret_WRITE_FAIL;
    }
    if (!final_set) {
      LOGPRINTFWW("--show-tags all: Report written to %s .\n", outname);
    }
  }
  if (final_set) {
    memcpy(outname_end, ".tags", 6);
    if (fopen_checked(outname, "w", &outfile)) {
      goto show_tags_ret_OPEN_FAIL;
    }
    if (!twocolumn) {
      marker_uidx = next_set(final_set, 0, unfiltered_marker_ct);
      while (marker_uidx < unfiltered_marker_ct) {
	fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
	putc_unlocked('\n', outfile);
	marker_uidx++;
	next_set_ul_ck(final_set, unfiltered_marker_ct, &marker_uidx);
      }
    } else {
      for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
        next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
	putc_unlocked('\t', outfile);
        putc_unlocked('0' + IS_SET(final_set, marker_uidx), outfile);
        putc_unlocked('\n', outfile);
      }
    }
    if (fclose_null(&outfile)) {
      goto show_tags_ret_WRITE_FAIL;
    }
    uii = popcount_longs(final_set, unfiltered_marker_ctl) - target_ct;
    if (tags_list) {
      LOGPRINTFWW("--show-tags: Main report written to %s.list , and simple tag ID list (%u tag%s added) written to %s .\n", outname, uii, (uii == 1)? "" : "s", outname);
    } else {
      LOGPRINTFWW("--show-tags: Simple tag ID list (%u tag%s added) written to %s .\n", uii, (uii == 1)? "" : "s", outname);
    }
  }
  while (0) {
  show_tags_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  show_tags_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  show_tags_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  show_tags_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  show_tags_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 show_tags_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(infile);
  fclose_cond(outfile);
  return retval;
}

double calc_lnlike_quantile(double known11, double known12, double known21, double known22, double unknown_dh, double freqx1, double freq1x, double freq2x, double freq11_expected, double denom, int32_t quantile) {
  // almost identical to calc_lnlike, but we can skip the equal-to-zero checks
  // when quantile isn't 100
  double tmp11 = quantile * denom + freq11_expected;
  double tmp12 = freq1x - tmp11;
  double tmp21 = freqx1 - tmp11;
  double tmp22 = freq2x - tmp21;
  if (quantile == 100) {
    // One of these values will be ~zero, and we want to ensure its logarithm
    // is treated as a very negative number instead of nan.  May as well do it
    // the same way as Haploview.
    if (tmp11 < 1e-10) {
      tmp11 = 1e-10;
    }
    if (tmp12 < 1e-10) {
      tmp12 = 1e-10;
    }
    if (tmp21 < 1e-10) {
      tmp21 = 1e-10;
    }
    if (tmp22 < 1e-10) {
      tmp22 = 1e-10;
    }
  }
  return known11 * log(tmp11) + known12 * log(tmp12) + known21 * log(tmp21) + known22 * log(tmp22) + unknown_dh * log(tmp11 * tmp22 + tmp12 * tmp21);
}

uint32_t haploview_blocks_classify(uint32_t* counts, uint32_t lowci_max, uint32_t lowci_min, uint32_t recomb_highci, uint32_t strong_highci, uint32_t strong_lowci, uint32_t strong_lowci_outer, uint32_t is_x, double recomb_fast_ln_thresh) {
  // See comments in the middle of haploview_blocks().  The key insight is that
  // we only need to classify the D' confidence intervals into a few types, and
  // this almost never requires evaluation of all 101 log likelihoods.

  // Note that lowCI and highCI are *one-sided* 95% confidence bounds, i.e.
  // together, they form a 90% confidence interval.
  double known11 = (double)(2 * counts[0] + counts[1] + counts[3]);
  double known12 = (double)(2 * counts[2] + counts[1] + counts[5]);
  double known21 = (double)(2 * counts[6] + counts[3] + counts[7]);
  double known22 = (double)(2 * counts[8] + counts[5] + counts[7]);
  double total_prob = 0.0;
  double lnsurf_highstrong_thresh = 0.0;
  uint32_t onside_sol_ct = 1;
  double right_sum[83];
  double freq1x;
  double freq2x;
  double freqx1;
  double freqx2;
  double freq11_expected;
  double unknown_dh;
  double denom;
  double lnlike1;
  double lnsurf_highindiff_thresh;
  double dxx;
  double dyy;
  double dzz;
  uint32_t quantile;
  uint32_t center;
  if (is_x) {
    known11 -= (double)((int32_t)counts[9]);
    known12 -= (double)((int32_t)counts[11]);
    known21 -= (double)((int32_t)counts[12]);
    known22 -= (double)((int32_t)counts[14]);
  }
  if (em_phase_hethet(known11, known12, known21, known22, counts[4], &freq1x, &freq2x, &freqx1, &freqx2, &dzz, &onside_sol_ct)) {
    return 1;
  }
  freq11_expected = freqx1 * freq1x;
  dxx = dzz - freq11_expected;
  if (dxx < 0.0) {
    // D < 0, flip (1,1)<->(1,2) and (2,1)<->(2,2) to make D positive
    dyy = known11;
    known11 = known12;
    known12 = dyy;
    dyy = known21;
    known21 = known22;
    known22 = dyy;
    freq11_expected = freqx2 * freq1x;
    dyy = freqx1;
    freqx1 = freqx2;
    freqx2 = dyy;
    dxx = -dxx;
  }
  dyy = MINV(freqx1 * freq2x, freqx2 * freq1x);
  // this will always be in a term with a 0.01 multiplier from now on, so may
  // as well premultiply.
  denom = 0.01 * dyy;
  unknown_dh = (double)((int32_t)counts[4]);

  // force this to an actual likelihood array entry, so we know for sure
  // total_prob >= 1.0 and can use that inequality for both early exit and
  // determining the "futility threshold" (terms smaller than 2^{-53} / 19 are
  // too small to matter).
  center = (int32_t)(((dxx / dyy) * 100) + 0.5);

  lnlike1 = calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, center);

  // Previously assumed log likelihood was always concave, and used geometric
  // series bounds... then I realized this was NOT a safe assumption to make.
  // See e.g. rs9435793 and rs7531410 in 1000 Genomes phase 1.
  // So, instead, we only use an aggressive approach when onside_sol_ct == 1
  // (fortunately, that is almost always the case).
  if (onside_sol_ct == 1) {
    // It's not actually necessary to keep the entire likelihood array in
    // memory.  This is similar to the HWE and Fisher's exact test
    // calculations: we can get away with tracking a few partial sums, and
    // exploit unimodality, fixed direction on both sides of the center,
    // knowledge of the center's location, and the fact that we only need to
    // classify the CI rather than fully evaluate it.
    //
    // Specifically, we need to determine the following:
    // 1. Is highCI >= 0.98?  Or < 0.90?
    // 2. If highCI >= 0.98, is lowCI >= 0.81?  In [0.71, 0.81)?  Equal to
    //    0.70?  In [0.51, 0.70)?  In [0.01, 0.51)?  Or < 0.01?
    //    (Crucially, if highCI < 0.98, we don't actually need to determine
    //    lowCI at all.)
    // To make this classification with as few relative likelihood evaluations
    // as possible (5 logs, an exp call, 8 multiplies, 9 adds... that's kinda
    // heavy for an inner loop operation), we distinguish the following cases:
    // a. D' >= 0.41.  We first try to quickly rule out highCI >= 0.98 by
    //    inspection of f(0.97).  Then,
    //    * If it's below the futility threshold, jump to case (b).
    //    * Otherwise, sum f(0.98)..f(1.00), and then sum other likelihoods
    //      from f(0.96) on down.
    // b. D' < 0.41.  highCI >= 0.98 is impossible since f(0.41) >= f(0.42) >=
    //    ...; goal is to quickly establish highCI < 0.90.  A large fraction of
    //    the time, this can be accomplished simply by inspecting f(0.89); if
    //    it's less than 1/220, we're done because we know there's a 1
    //    somewhere in the array, and the sum of the likelihoods between
    //    f(0.89) and whatever that 1 entry is is bounded above by 12 * (1/220)
    //    due to fixed direction.  Otherwise, we sum from the top down.
    // This should be good for a ~10x speedup on the larger datasets where it's
    // most wanted.
    if (100 - center < 20 * (100 - strong_highci)) {
      dxx = calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, strong_highci) - lnlike1;
      // ln(2^{-53} / 19) is just under -39.6812
      if ((center > strong_highci) || (dxx > -39.6812)) {
	total_prob = exp(dxx);
	for (quantile = 100; quantile > strong_highci; quantile--) {
	  total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
	}
	if (total_prob > (1.0 / 19.0)) {
	  // branch 1: highCI might be >= 0.98
	  lnsurf_highstrong_thresh = total_prob * 20;
	  for (quantile = strong_highci - 1; quantile >= recomb_highci; quantile--) {
	    total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
	  }
	  lnsurf_highindiff_thresh = total_prob * 20;
	  while (1) {
	    dxx = exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
	    total_prob += dxx;
	    // see comments on branch 2.  this is more complicated because we
	    // still have work to do after resolving whether highCI >= 0.98,
	    // but the reasoning is similar.
	    if (total_prob >= lnsurf_highstrong_thresh) {
	      if (quantile >= center) {
	        goto haploview_blocks_classify_no_highstrong_1;
	      }
	      goto haploview_blocks_classify_no_highstrong_2;
	    }
	    if ((quantile <= lowci_max) && (quantile >= lowci_min)) {
	      // We actually only need the [52..100], [71..100], [72..100], and
	      // [82..100] right sums, but saving a few extra values is
	      // probably more efficient than making this if-statement more
	      // complicated.  [99 - quantile] rather than e.g. [quantile]
	      // is used so memory writes go to sequentially increasing rather
	      // than decreasing addresses.  (okay, this shouldn't matter since
	      // everything should be in L1 cache, but there's negligible
	      // opportunity cost)
	      right_sum[quantile] = total_prob;
	    }
	    dxx *= ((int32_t)quantile);
	    if (total_prob + dxx < lnsurf_highstrong_thresh) {
	      while (1) {
		// Now we want to bound lowCI, optimizing for being able to
		// quickly establish lowCI >= 0.71.
		if (dxx * 19 < total_prob) {
		  // less than 5% remaining on left tail
		  if (quantile >= lowci_max) {
		    return 6;
		  }
		  while (quantile > lowci_min) {
		    quantile--;
		    total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
		    if (quantile <= lowci_max) {
		      right_sum[quantile] = total_prob;
		    }
		  }
		  dyy = right_sum[lowci_min] * (20.0 / 19.0);
		  while (total_prob < dyy) {
		    if ((!quantile) || (dxx <= RECIP_2_53)) {
		      total_prob *= 0.95;
		      if (total_prob >= right_sum[strong_lowci_outer]) {
			// lowCI < 0.70
			// -> f(0.00) + f(0.01) + ... + f(0.70) > 0.05 * total
			return 3;
		      } else if (total_prob < right_sum[lowci_max]) {
			return 6;
		      } else if ((lowci_max > strong_lowci) && (total_prob < right_sum[strong_lowci])) {
			return 5;
		      }
		      return 4;
		    }
		    quantile--;
		    dxx = exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
		    total_prob += dxx;
		  }
		  return 2;
		}
		quantile--;
		dxx = exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
		total_prob += dxx;
		if ((quantile <= lowci_max) && (quantile >= lowci_min)) {
		  right_sum[quantile] = total_prob;
		}
		dxx *= ((int32_t)quantile);
	      }
	    }
	    quantile--;
	  }
	}
      }
      quantile = strong_highci - 1;
    } else {
      quantile = 100;
    }
    // branch 2: highCI guaranteed less than 0.98.  If D' <= 0.875, try to
    // quickly establish highCI < 0.90.
    dxx = calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, recomb_highci) - lnlike1;
    if ((center < recomb_highci) && (dxx < recomb_fast_ln_thresh)) {
      return 0;
    }
    // okay, we'll sum the whole right tail.  May as well sum from the outside
    // in here for a bit more numerical stability, instead of adding exp(dxx)
    // first.
    do {
      total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
    } while (--quantile > recomb_highci);
    total_prob += exp(dxx);
    lnsurf_highindiff_thresh = total_prob * 20;
  haploview_blocks_classify_no_highstrong_1:
    quantile--;
    if (center < recomb_highci) {
      // if we know there's a 1.0 ahead in the likelihood array, may as well
      // take advantage of that
      lnsurf_highstrong_thresh = lnsurf_highindiff_thresh - 1.0;
      while (quantile > center) {
	total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
	if (total_prob >= lnsurf_highstrong_thresh) {
	  return 0;
	}
	quantile--;
      }
      if (!center) {
	return 1;
      }
      total_prob += 1;
      quantile--;
    }
    // likelihoods are now declining, try to exploit that to exit early
    // (it's okay if the first likelihood does not represent a decline)
    while (1) {
      dxx = exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
      total_prob += dxx;
    haploview_blocks_classify_no_highstrong_2:
      if (total_prob >= lnsurf_highindiff_thresh) {
	return 0;
      }
      if (total_prob + ((int32_t)quantile) * dxx < lnsurf_highindiff_thresh) {
	// guaranteed to catch quantile == 0
	return 1;
      }
      quantile--;
    }
  }
  for (quantile = 100; quantile >= recomb_highci; quantile--) {
    total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
    if (quantile == strong_highci) {
      lnsurf_highstrong_thresh = total_prob * 20;
    }
  }
  if (total_prob < (1.0 / 19.0)) {
    return 0;
  }
  lnsurf_highindiff_thresh = total_prob * 20;
  while (1) {
    total_prob += exp(calc_lnlike_quantile(known11, known12, known21, known22, unknown_dh, freqx1, freq1x, freq2x, freq11_expected, denom, quantile) - lnlike1);
    if (total_prob >= lnsurf_highindiff_thresh) {
      return 0;
    }
    if (quantile <= lowci_max) {
      if (quantile >= lowci_min) {
        right_sum[quantile] = total_prob;
      } else if (!quantile) {
	break;
      }
    }
    quantile--;
  }
  if (total_prob >= lnsurf_highstrong_thresh) {
    return 1;
  }
  total_prob *= 0.95;
  if (total_prob < right_sum[strong_lowci]) {
    if ((lowci_max > strong_lowci) && (total_prob >= right_sum[lowci_max])) {
      return 5;
    }
    return 6;
  }
  if (total_prob >= right_sum[strong_lowci_outer]) {
    if ((lowci_min < strong_lowci_outer) && (total_prob >= right_sum[lowci_min])) {
      return 2;
    }
    return 3;
  }
  return 4;
}

int32_t haploview_blocks(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uintptr_t* pheno_nm, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  // See Plink::mkBlks() in blox.cpp (which is, in turn, a port of doGabriel()
  // in FindBlocks.java and computeDPrime() in HaploData.java from Haploview).
  // No unwindowed/inter-chr mode, so little point in bothering with
  // multithreading.
  //
  // MAF < 0.05 markers have a minor effect on PLINK 1.07 --blocks's behavior
  // when present, while Haploview completely ignores them.  We replicate
  // Haploview's behavior.
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  FILE* outfile = nullptr;
  FILE* outfile_det = nullptr;
  // circular.  [2n] = numStrong, [2n+1] = numRec
  uintptr_t* strong_rec_cts = nullptr;
  uintptr_t* founder_include2 = nullptr;
  uintptr_t* founder_male_include2 = nullptr;
  uintptr_t marker_uidx = 0;
  uintptr_t block_idx_first = 0;
  uintptr_t block_uidx_first = 0;
  uintptr_t block_pos_first = 0;
  uintptr_t prev_strong = 0;
  uintptr_t prev_rec = 0;
  uintptr_t markers_done = 0;
  uint32_t no_pheno_req = ldip->modifier & LD_BLOCKS_NO_PHENO_REQ;
  uint32_t max_window_bp = ldip->blocks_max_bp;
  uint32_t max_window_bp1 = 20000;
  uint32_t max_window_bp2 = 30000;
  uint32_t recomb_highci = ldip->blocks_recomb_highci;
  uint32_t strong_highci = ldip->blocks_strong_highci;
  uint32_t strong_lowci = ldip->blocks_strong_lowci;
  uint32_t strong_lowci_outer = ldip->blocks_strong_lowci_outer;
  uint32_t block_ct = 0;
  uint32_t maxspan = 0;
  uint32_t pct = 0;
  int32_t retval = 0;
  double recomb_fast_ln_thresh = -log((int32_t)((100 - recomb_highci) * 20));
  double inform_frac = ldip->blocks_inform_frac + SMALLISH_EPSILON;
  uint32_t inform_thresh_two = 1 + ((int32_t)(3 * inform_frac));
  uint32_t inform_thresh_three = (int32_t)(6 * inform_frac);
  uint32_t counts[15];
  // [0]: (m, m-1)
  // [1]: (m, m-2)
  // [2]: (m-1, m-2)
  // [3]: (m-1, m-3)
  // [4]: (m-2, m-3)
  uint32_t recent_ci_types[5];
  uint32_t index_tots[5];
  uintptr_t* founder_pnm;
  uintptr_t* marker_exclude;
  uintptr_t* in_haploblock;
  uintptr_t* loadbuf_raw;
  uintptr_t* index_data;
  uintptr_t* window_data;
  uintptr_t* window_data_ptr;
  unsigned char* bigstack_mark2;
  uint32_t* block_uidxs;
  uint32_t* forward_block_sizes;
  uint32_t* candidate_pairs;
  char* wptr_start;
  char* wptr;
  char* sptr;
  uintptr_t cur_marker_ct;
  uintptr_t max_block_size;
  uintptr_t marker_idx;
  uintptr_t cur_block_size;
  uintptr_t last_block_size;
  uintptr_t founder_ct;
  uintptr_t founder_ctl2;
  uintptr_t founder_ctv2;
  uintptr_t final_mask;
  uintptr_t futility_rec;
  uintptr_t max_candidates;
  uintptr_t candidate_ct;
  uintptr_t candidate_idx;
  uintptr_t delta;
  uintptr_t pct_thresh;
  uintptr_t ulii;
  double min_maf;
  double max_maf;
  double dxx;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_start;
  uint32_t chrom_end;
  uint32_t is_haploid;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t marker_pos_thresh;
  uint32_t forward_scan_uidx;
  uint32_t block_cidx;
  uint32_t block_cidx2;
  uint32_t cur_strong;
  uint32_t cur_rec;
  uint32_t lowci_max;
  uint32_t lowci_min;
  uint32_t cur_ci_type;
  uint32_t cur_marker_pos;
  uint32_t uii;
  uint32_t ujj;
  // suppress warning
  index_tots[3] = 0;
  index_tots[4] = 0;
  if (ldip->modifier & LD_BLOCKS_NO_SMALL_MAX_SPAN) {
    max_window_bp1 = 0x7fffffff;
    max_window_bp2 = 0x7fffffff;
  }

  // First enforce MAF 0.05 minimum; then, on each chromosome:
  // 1. Determine maximum number of markers that might need to be loaded at
  //    once on current chromosome, and then (re)allocate memory buffers.
  // 2. Find all pairs of markers satisfying the "strong LD" and informative
  //    fraction criteria.  (The original algorithm deferred the informative
  //    fraction calculation; we don't do that because it forces nonsequential
  //    file access.)
  // 3. Sort the pairs in decreasing order primarily by bp distance, and
  //    secondarily by start uidx.
  // 4. Greedily construct blocks from the sorted list (i.e. form largest
  //    blocks first).
  if (bigstack_alloc_ul(unfiltered_sample_ctl, &founder_pnm) ||
      bigstack_alloc_ul(unfiltered_marker_ctl, &marker_exclude) ||
      bigstack_alloc_ul(unfiltered_marker_ctl, &in_haploblock) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw)) {
    goto haploview_blocks_ret_NOMEM;
  }
  memcpy(founder_pnm, founder_info, unfiltered_sample_ctl * sizeof(intptr_t));
  if (!no_pheno_req) {
    bitvec_and(pheno_nm, unfiltered_sample_ctl, founder_pnm);
  }
  founder_ct = popcount_longs(founder_pnm, unfiltered_sample_ctl);
  if (founder_ct < 2) {
    if ((!no_pheno_req) && (!popcount_longs(pheno_nm, unfiltered_sample_ctl))) {
      logerrprint("Warning: Skipping --blocks, since there are less than two founders with\nnonmissing phenotypes.  (The 'no-pheno-req' modifier removes the phenotype\nrestriction.)\n");
    } else {
      logerrprint("Warning: Skipping --blocks, since there are less than two founders with\nnonmissing phenotypes.  (--make-founders may come in handy here.)\n");
    }
    goto haploview_blocks_ret_1;
  }
  final_mask = get_final_mask(founder_ct);
  memcpy(marker_exclude, marker_exclude_orig, unfiltered_marker_ctl * sizeof(intptr_t));
  if (ldip->blocks_min_maf > 0.0) {
    min_maf = ldip->blocks_min_maf * (1 - SMALL_EPSILON);
    max_maf = 1 - min_maf;
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude_orig, &marker_uidx);
      dxx = set_allele_freqs[marker_uidx];
      if ((dxx < min_maf) || (dxx > max_maf)) {
	set_bit_ul(marker_uidx, marker_exclude);
      }
    }
    marker_ct = unfiltered_marker_ct - popcount_longs(marker_exclude, unfiltered_marker_ctl);
  }
  if (marker_ct < 2) {
    logerrprint("Warning: Skipping --blocks since there are too few variants with MAF >= 0.05.\n");
    goto haploview_blocks_ret_1;
  }
  pct_thresh = marker_ct / 100;
  fill_ulong_zero(unfiltered_marker_ctl, in_haploblock);
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  founder_ctl2 = QUATERCT_TO_WORDCT(founder_ct);
  founder_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(founder_ct);
  if (bigstack_alloc_ul(5 * founder_ctv2, &index_data)) {
    goto haploview_blocks_ret_NOMEM;
  }
  if (alloc_collapsed_haploid_filters(founder_info, sex_male, unfiltered_sample_ct, founder_ct, Y_FIX_NEEDED, 1, &founder_include2, &founder_male_include2)) {
    goto haploview_blocks_ret_NOMEM;
  }
  memcpy(outname_end, ".blocks.det", 12);
  if (fopen_checked(outname, "w", &outfile_det)) {
    goto haploview_blocks_ret_OPEN_FAIL;
  }
  if (fputs_checked(" CHR          BP1          BP2           KB  NSNPS SNPS\n", outfile_det)) {
    goto haploview_blocks_ret_WRITE_FAIL;
  }
  outname_end[7] = '\0';
  if (fopen_checked(outname, "w", &outfile)) {
    goto haploview_blocks_ret_OPEN_FAIL;
  }
  bigstack_mark2 = g_bigstack_base;
  fputs("--blocks: 0%", stdout);
  fflush(stdout);
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++, markers_done += cur_marker_ct) {
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    chrom_start = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
    cur_marker_ct = chrom_end - chrom_start - popcount_bit_idx(marker_exclude, chrom_start, chrom_end);
    if (cur_marker_ct < 2) {
      continue;
    }
    marker_uidx = chrom_start;
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    max_block_size = chrom_window_max(marker_pos, marker_exclude, chrom_info_ptr, chrom_idx, 0x7fffffff, max_window_bp, 1);
    if (max_block_size < 2) {
      continue;
    }
#ifndef __LP64__
    if (max_block_size > 65536) {
      logprint("\n");
      logerrprint("Error: 32-bit --blocks cannot analyze potential blocks with more than 65536\nvariants.  Use a 64-bit PLINK build or a smaller --blocks-window-kb value.\n");
      goto haploview_blocks_ret_INVALID_CMDLINE;
    }
#endif
    is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
    is_x = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[X_OFFSET]);
    is_y = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
    bigstack_reset(bigstack_mark2);
    // Need to compute full 3x3 count tables, but only for a limited window;
    // more similar to --clump than --fast-epistasis, so we don't bother with
    // precomputing 0-only/1-only/2-only bitfields or multithreading for now.

    // For each pair, we just need to know 100x the Haploview lowCI and highCI
    // values; lod and dp are unnecessary since the CI value also tracks bad
    // pairs.  More precisely, there are only seven types of CIs worth
    // distinguishing:
    // 0. non-bad pair, and highCI < recHighCI (0.90)
    // 1. "null" (bad pair, or highCI in [0.90, 0.98))
    // 2. highCI >= 0.98, and lowCI < 0.51
    //    (treated the same as type 1, but it takes no additional effort to
    //    distinguish this case)
    // 3. highCI >= 0.98, and lowCI in [0.51, 0.70)
    // 4. highCI >= 0.98, and lowCI == 0.70
    //    (turns out (double)70 / 100.0 compares exactly equal to 0.70, so
    //    Haploview's use of < cutLowCI in its initial "strong LD" check
    //    actually behaves differently from the later "not > cutLowCI" check)
    // 5. highCI >= 0.98, and lowCI in [0.71, 0.81)
    // 6. highCI >= 0.98, and lowCI in [0.81, 1]
    // And it gets better than that: given block size n, we just need to
    // maintain #(type 0) and #(type 4/5/6) arrays (and a tiny array with more
    // detailed information on the most recent blocks) to find all potentially
    // valid blocks in a single pass.  So we can use practically all our memory
    // to track and sort those blocks by bp length.
    if (bigstack_alloc_ui(max_block_size, &block_uidxs) ||
        bigstack_alloc_ui(max_block_size, &forward_block_sizes) ||
        bigstack_alloc_ul(max_block_size * founder_ctv2, &window_data)) {
      goto haploview_blocks_ret_NOMEM;
    }
    if (max_block_size >= 4) {
      // After marker m is fully processed,
      //   strong_rec_cts[(block_cidx + delta) * 2] = numStrong, and
      //   strong_rec_cts[(block_cidx + delta) * 2 + 1] = numRec
      // for the potential [m - delta, m] block, taking array indices modulo
      // max_block_size * 2.
      if (bigstack_alloc_ul(max_block_size * 2, &strong_rec_cts)) {
	goto haploview_blocks_ret_NOMEM;
      }
    }
    window_data_ptr = &(window_data[founder_ctv2 - 2]);
    for (ulii = 0; ulii < max_block_size; ulii++) {
      window_data_ptr[0] = 0;
      window_data_ptr[1] = 0;
      window_data_ptr = &(window_data_ptr[founder_ctv2]);
    }
    block_idx_first = 0;
    block_uidx_first = chrom_start;
    marker_uidx = chrom_start;
    block_pos_first = marker_pos[chrom_start];
    max_candidates = bigstack_left() / (3 * sizeof(int32_t));
    bigstack_alloc_ui(max_candidates * 3, &candidate_pairs);
    candidate_ct = 0;
    cur_block_size = 0;
    fill_uint_zero(3, recent_ci_types);
    // count down instead of up so more memory accesses are sequential
    block_cidx = max_block_size;
    forward_scan_uidx = marker_uidx;
    if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto haploview_blocks_ret_READ_FAIL;
    }
    for (marker_idx = 0; marker_idx < cur_marker_ct; marker_uidx++, marker_idx++) {
      if (block_cidx) {
        block_cidx--;
      } else {
	block_cidx = max_block_size - 1;
      }
      window_data_ptr = &(window_data[block_cidx * founder_ctv2]);
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
        if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
          goto haploview_blocks_ret_READ_FAIL;
	}
      }
      block_uidxs[block_cidx] = marker_uidx;
      if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_pnm, final_mask, 0, bedfile, loadbuf_raw, window_data_ptr)) {
	goto haploview_blocks_ret_READ_FAIL;
      }
      if (is_haploid) {
	haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)window_data_ptr);
      }
      cur_marker_pos = marker_pos[marker_uidx];
      marker_pos_thresh = cur_marker_pos;
      if (marker_pos_thresh < max_window_bp) {
	marker_pos_thresh = 0;
      } else {
	marker_pos_thresh -= max_window_bp;
      }
      if (marker_pos_thresh > block_pos_first) {
	do {
	  block_uidx_first++;
	  next_unset_ul_unsafe_ck(marker_exclude, &block_uidx_first);
	  block_pos_first = marker_pos[block_uidx_first];
	  block_idx_first++;
	} while (marker_pos_thresh > block_pos_first);
      }
      last_block_size = cur_block_size;
      cur_block_size = marker_idx - block_idx_first;
      recent_ci_types[4] = recent_ci_types[2];
      recent_ci_types[2] = recent_ci_types[0];
      recent_ci_types[3] = recent_ci_types[1];
      if (cur_block_size > last_block_size) {
	cur_block_size = last_block_size + 1;
      }
      // now determine maximum local block size, so we can set futility_rec
      // efficiently.
      marker_pos_thresh = cur_marker_pos + max_window_bp;
      if (forward_scan_uidx < marker_uidx) {
	forward_scan_uidx = marker_uidx;
      }
      while (marker_pos_thresh >= marker_pos[forward_scan_uidx]) {
	uii = forward_scan_uidx + 1;
	next_unset_ck(marker_exclude, chrom_end, &uii);
	if (uii == chrom_end) {
	  break;
	}
        forward_scan_uidx = uii;
      }
      uii = forward_scan_uidx + 1 - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, forward_scan_uidx);
      forward_block_sizes[block_cidx] = uii;
      if (!cur_block_size) {
	continue;
      }
      block_cidx2 = block_cidx + 1;
      for (delta = 1; delta <= cur_block_size; delta++, block_cidx2++) {
	if (block_cidx2 == max_block_size) {
	  block_cidx2 = 0;
	}
        if (forward_block_sizes[block_cidx2] > uii) {
	  uii = forward_block_sizes[block_cidx2];
	}
      }
      ulii = uii;
      // If numRec ever reaches this value, we can just move on to the next
      // marker (even skipping the remaining D' evaluations).
      futility_rec = 1 + (((double)((intptr_t)((ulii * (ulii - 1)) / 2))) * (1.0 - inform_frac));
      block_cidx2 = block_cidx + 1;
      cur_strong = 0;
      cur_rec = 0;
      vec_datamask(founder_ct, 0, window_data_ptr, founder_include2, index_data);
      index_tots[0] = popcount2_longs(index_data, founder_ctl2);
      vec_datamask(founder_ct, 2, window_data_ptr, founder_include2, &(index_data[founder_ctv2]));
      index_tots[1] = popcount2_longs(&(index_data[founder_ctv2]), founder_ctl2);
      vec_datamask(founder_ct, 3, window_data_ptr, founder_include2, &(index_data[2 * founder_ctv2]));
      index_tots[2] = popcount2_longs(&(index_data[2 * founder_ctv2]), founder_ctl2);
      if (is_x) {
	vec_datamask(founder_ct, 0, window_data_ptr, founder_male_include2, &(index_data[3 * founder_ctv2]));
	index_tots[3] = popcount2_longs(&(index_data[3 * founder_ctv2]), founder_ctl2);
	vec_datamask(founder_ct, 3, window_data_ptr, founder_male_include2, &(index_data[4 * founder_ctv2]));
	index_tots[4] = popcount2_longs(&(index_data[4 * founder_ctv2]), founder_ctl2);
      }
      lowci_max = 82;
      lowci_min = 52;
      for (delta = 1; delta <= cur_block_size; delta++, block_cidx2++) {
	if (block_cidx2 == max_block_size) {
	  block_cidx2 = 0;
	}
	if (delta >= 4) {
	  prev_rec = strong_rec_cts[block_cidx2 * 2 + 1];
	  if (cur_rec + prev_rec >= futility_rec) {
	    cur_block_size = delta - 1;
	    break;
	  }
          prev_strong = strong_rec_cts[block_cidx2 * 2];
	}
	window_data_ptr = &(window_data[block_cidx2 * founder_ctv2]);
	genovec_3freq(window_data_ptr, index_data, founder_ctl2, &(counts[0]), &(counts[1]), &(counts[2]));
	counts[0] = index_tots[0] - counts[0] - counts[1] - counts[2];
	genovec_3freq(window_data_ptr, &(index_data[founder_ctv2]), founder_ctl2, &(counts[3]), &(counts[4]), &(counts[5]));
	counts[3] = index_tots[1] - counts[3] - counts[4] - counts[5];
	genovec_3freq(window_data_ptr, &(index_data[2 * founder_ctv2]), founder_ctl2, &(counts[6]), &(counts[7]), &(counts[8]));
	counts[6] = index_tots[2] - counts[6] - counts[7] - counts[8];
	if (is_x) {
	  genovec_3freq(window_data_ptr, &(index_data[3 * founder_ctv2]), founder_ctl2, &(counts[9]), &(counts[10]), &(counts[11]));
	  // counts[10] should always be zero
	  counts[9] = index_tots[3] - counts[9] - counts[11];
	  genovec_3freq(window_data_ptr, &(index_data[4 * founder_ctv2]), founder_ctl2, &(counts[12]), &(counts[13]), &(counts[14]));
	  counts[12] = index_tots[4] - counts[12] - counts[14];
	}
	cur_ci_type = haploview_blocks_classify(counts, lowci_max, lowci_min, recomb_highci, strong_highci, strong_lowci, strong_lowci_outer, is_x, recomb_fast_ln_thresh);
	if (cur_ci_type > 4) {
	  cur_strong++;
	} else if (!cur_ci_type) {
	  cur_rec++;
	}
	if (delta < 4) {
	  if (delta == 1) {
	    lowci_max = strong_lowci;
	    recent_ci_types[0] = cur_ci_type;
	    if ((cur_ci_type == 6) && (cur_marker_pos - marker_pos[block_uidxs[block_cidx2]] <= max_window_bp1)) {
	      goto haploview_blocks_save_candidate;
	    }
	  } else if (delta == 2) {
	    recent_ci_types[1] = cur_ci_type;
	    if ((cur_ci_type >= 4) && (cur_marker_pos - marker_pos[block_uidxs[block_cidx2]] <= max_window_bp2)) {
	      uii = 1;
	      if (recent_ci_types[0] >= 3) {
		uii++;
	      }
	      if (recent_ci_types[2] >= 3) {
		uii++;
	      }
	      if (uii >= inform_thresh_two) {
	        goto haploview_blocks_save_candidate;
	      }
	    }
	  } else {
	    lowci_min = strong_lowci_outer;
	    prev_strong = 0; // 5+
	    uii = 0; // 3+, not counting cur_ci_type
	    prev_rec = 0;
	    if (cur_ci_type > 4) {
	      prev_strong++;
	    } else if (!cur_ci_type) {
	      prev_rec++;
	    }
	    for (ujj = 0; ujj < 5; ujj++) {
	      if (recent_ci_types[ujj] >= 3) {
		uii++;
		if (recent_ci_types[ujj] > 4) {
		  prev_strong++;
		}
	      } else if (!recent_ci_types[ujj]) {
		prev_rec++;
	      }
	    }
	    strong_rec_cts[block_cidx2 * 2] = prev_strong;
	    strong_rec_cts[block_cidx2 * 2 + 1] = prev_rec;
	    if ((cur_ci_type >= 4) && (uii >= inform_thresh_three)) {
	      goto haploview_blocks_save_candidate;
	    }
	  }
	} else {
	  prev_strong += cur_strong;
	  prev_rec += cur_rec;
	  strong_rec_cts[block_cidx2 * 2] = prev_strong;
	  strong_rec_cts[block_cidx2 * 2 + 1] = prev_rec;
	  ulii = prev_strong + prev_rec;
	  if ((cur_ci_type >= 4) && (ulii >= 6) && (((intptr_t)ulii) * inform_frac < ((double)((intptr_t)prev_strong)))) {
	  haploview_blocks_save_candidate:
	    if (candidate_ct == max_candidates) {
	      goto haploview_blocks_ret_NOMEM;
	    }
	    uii = block_uidxs[block_cidx2];
	    candidate_pairs[3 * candidate_ct] = cur_marker_pos - marker_pos[uii];
	    candidate_pairs[3 * candidate_ct + 1] = uii;
	    candidate_pairs[3 * candidate_ct + 2] = marker_uidx;
	    candidate_ct++;
	  }
	}
      }
      if (markers_done + marker_idx >= pct_thresh) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	pct = ((markers_done + marker_idx) * 100LLU) / marker_ct;
	printf("\b\b%u%%", pct++);
	fflush(stdout);
	pct_thresh = (pct * marker_ct) / 100;
      }
    }
    if (!candidate_ct) {
      continue;
    }
    qsort(candidate_pairs, candidate_ct, 12, intcmp3_decr);
    if (candidate_pairs[0] > maxspan) {
      maxspan = candidate_pairs[0];
    }
    ulii = 0; // final haploblock count
    for (candidate_idx = 0; candidate_idx < candidate_ct; candidate_idx++) {
      block_cidx = candidate_pairs[candidate_idx * 3 + 1];
      if (is_set(in_haploblock, block_cidx)) {
	continue;
      }
      block_cidx2 = candidate_pairs[candidate_idx * 3 + 2];
      if (is_set(in_haploblock, block_cidx2)) {
	continue;
      }
      candidate_pairs[2 * ulii] = block_cidx;
      candidate_pairs[2 * ulii + 1] = block_cidx2;
      fill_bits(block_cidx, block_cidx2 + 1 - block_cidx, in_haploblock);
      ulii++;
    }
#ifdef __cplusplus
    std::sort((int64_t*)candidate_pairs, (int64_t*)(&(candidate_pairs[ulii * 2])));
#else
    qsort(candidate_pairs, ulii, sizeof(int64_t), llcmp);
#endif
    wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf));
    wptr_start = memseta(wptr_start, 32, 3);
    for (candidate_idx = 0; candidate_idx < ulii; candidate_idx++) {
      putc_unlocked('*', outfile);
      block_cidx = candidate_pairs[2 * candidate_idx];
      block_cidx2 = candidate_pairs[2 * candidate_idx + 1];
      marker_uidx = block_cidx;
      wptr = uint32toa_w10(marker_pos[block_cidx], wptr_start);
      wptr = memseta(wptr, 32, 3);
      wptr = uint32toa_w10x(marker_pos[block_cidx2], ' ', wptr);
      wptr = width_force(12, wptr, dtoa_g(((int32_t)(marker_pos[block_cidx2] + 1 - marker_pos[block_cidx])) * 0.001, wptr));
      *wptr++ = ' ';
      wptr = uint32toa_w6x(block_cidx2 + 1 - block_cidx - popcount_bit_idx(marker_exclude, block_cidx, block_cidx2), ' ', wptr);
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile_det)) {
	goto haploview_blocks_ret_WRITE_FAIL;
      }
      for (marker_uidx = block_cidx; marker_uidx <= block_cidx2; marker_uidx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
        sptr = &(marker_ids[marker_uidx * max_marker_id_len]);
        putc_unlocked(' ', outfile);
        fputs(sptr, outfile);
        if (marker_uidx != block_cidx) {
	  putc_unlocked('|', outfile_det);
	}
	fputs(sptr, outfile_det);
      }
      putc_unlocked('\n', outfile);
      putc_unlocked('\n', outfile_det);
    }
    block_ct += ulii;
  }
  if (fclose_null(&outfile)) {
    goto haploview_blocks_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile_det)) {
    goto haploview_blocks_ret_WRITE_FAIL;
  }
  putc_unlocked('\r', stdout);
  LOGPRINTFWW("--blocks: %u haploblock%s written to %s .\n", block_ct, (block_ct == 1)? "" : "s", outname);
  LOGPRINTFWW("Extra block details written to %s.det .\n", outname);
  if (block_ct) {
    LOGPRINTF("Longest span: %gkb.\n", ((double)(maxspan + 1)) * 0.001);
  }
  while (0) {
  haploview_blocks_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  haploview_blocks_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  haploview_blocks_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  haploview_blocks_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
#ifndef __LP64__
  haploview_blocks_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
#endif
  }
 haploview_blocks_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_det);
  return retval;
}

void twolocus_write_table(FILE* outfile, uint32_t* counts, uint32_t plink_maxsnp, char* mkr1, char* mkr2, char* allele00, char* allele01, char* allele10, char* allele11, uint32_t alen00, uint32_t alen01, uint32_t alen10, uint32_t alen11) {
  // PLINK 1.07's print settings for this function don't handle large numbers
  // well so we break byte-for-byte compatibility.
  char* bufptr = memseta(g_textbuf, 32, plink_maxsnp + 14);
  uint32_t* uiptr = counts;
  uint32_t total = 0;
  uint32_t marg_a[4];
  uint32_t marg_b[4];
  char spaces[7];
  double tot_recip;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  fill_uint_zero(4, marg_b);
  memset(spaces, 32, 7);
  for (uii = 0; uii < 4; uii++) {
    ukk = 0;
    for (ujj = 0; ujj < 4; ujj++) {
      umm = *uiptr++;
      ukk += umm;
      marg_b[ujj] += umm;
    }
    marg_a[uii] = ukk;
    total += ukk;
  }
  tot_recip = 1.0 / ((double)((int32_t)total));
  bufptr = strcpyax(bufptr, mkr2, '\n');
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
  fwrite(g_textbuf, 1, plink_maxsnp + 7, outfile);
  if (alen10 < 4) {
    fwrite(spaces, 1, 9 - 2 * alen10, outfile);
  }
  fputs(allele10, outfile);
  putc_unlocked('/', outfile);
  fputs(allele10, outfile);
  putc_unlocked(' ', outfile);
  if (alen10 + alen11 < 7) {
    fwrite(spaces, 1, 9 - alen10 - alen11, outfile);
  }
  fputs(allele10, outfile);
  putc_unlocked('/', outfile);
  fputs(allele11, outfile);
  putc_unlocked(' ', outfile);
  if (alen11 < 4) {
    fwrite(spaces, 1, 9 - 2 * alen11, outfile);
  }
  fputs(allele11, outfile);
  putc_unlocked('/', outfile);
  fputs(allele11, outfile);
  fputs("        0/0        */*\n", outfile);

  bufptr = fw_strcpy(plink_maxsnp, mkr1, g_textbuf);
  *bufptr++ = ' ';
  if (alen00 == 1) {
    bufptr = memseta(bufptr, 32, 2);
  }
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
  fputs(allele00, outfile);
  putc_unlocked('/', outfile);
  fputs(allele00, outfile);
  bufptr = g_textbuf;
  *bufptr++ = ' ';
  bufptr = uint32toa_w10x(counts[0], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[2], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[3], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[1], ' ', bufptr);
  bufptr = uint32toa_w10x(marg_a[0], '\n', bufptr);
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);

  bufptr = memseta(g_textbuf, 32, plink_maxsnp + 1);
  if (alen00 + alen01 < 4) {
    bufptr = memseta(bufptr, 32, 4 - alen00 - alen01);
  }
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
  fputs(allele00, outfile);
  putc_unlocked('/', outfile);
  fputs(allele01, outfile);
  bufptr = g_textbuf;
  *bufptr++ = ' ';
  bufptr = uint32toa_w10x(counts[8], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[10], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[11], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[9], ' ', bufptr);
  bufptr = uint32toa_w10x(marg_a[2], '\n', bufptr);
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);

  bufptr = memseta(g_textbuf, 32, plink_maxsnp + 1);
  if (alen01 == 1) {
    bufptr = memseta(bufptr, 32, 2);
  }
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
  fputs(allele01, outfile);
  putc_unlocked('/', outfile);
  fputs(allele01, outfile);
  bufptr = g_textbuf;
  *bufptr++ = ' ';
  bufptr = uint32toa_w10x(counts[12], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[14], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[15], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[13], ' ', bufptr);
  bufptr = uint32toa_w10x(marg_a[3], '\n', bufptr);
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);

  bufptr = memseta(g_textbuf, 32, plink_maxsnp + 3);
  bufptr = memcpya(bufptr, "0/0 ", 4);
  bufptr = uint32toa_w10x(counts[4], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[6], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[7], ' ', bufptr);
  bufptr = uint32toa_w10x(counts[5], ' ', bufptr);
  bufptr = uint32toa_w10x(marg_a[1], '\n', bufptr);
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);

  bufptr = memseta(g_textbuf, 32, plink_maxsnp + 3);
  bufptr = memcpya(bufptr, "*/* ", 4);
  bufptr = uint32toa_w10x(marg_b[0], ' ', bufptr);
  bufptr = uint32toa_w10x(marg_b[2], ' ', bufptr);
  bufptr = uint32toa_w10x(marg_b[3], ' ', bufptr);
  bufptr = uint32toa_w10x(marg_b[1], ' ', bufptr);
  bufptr = uint32toa_w10x(total, '\n', bufptr);
  *bufptr++ = '\n';
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);

  bufptr = memseta(g_textbuf, 32, plink_maxsnp + 14);
  bufptr = strcpyax(bufptr, mkr2, '\n');
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
  fwrite(g_textbuf, 1, plink_maxsnp + 9, outfile);
  fputs(allele10, outfile);
  putc_unlocked('/', outfile);
  fputs(allele10, outfile);
  if (alen10 < 4) {
    fwrite(spaces, 1, 9 - 2 * alen10, outfile);
  }
  putc_unlocked(' ', outfile);
  fputs(allele10, outfile);
  putc_unlocked('/', outfile);
  fputs(allele11, outfile);
  if (alen10 + alen11 < 7) {
    fwrite(spaces, 1, 9 - alen10 - alen11, outfile);
  }
  putc_unlocked(' ', outfile);
  fputs(allele11, outfile);
  putc_unlocked('/', outfile);
  fputs(allele11, outfile);
  if (alen11 < 4) {
    fwrite(spaces, 1, 9 - 2 * alen11, outfile);
  }
  fputs(" 0/0        */*\n", outfile);

  bufptr = fw_strcpy(plink_maxsnp, mkr1, g_textbuf);
  *bufptr++ = ' ';
  if (alen00 == 1) {
    bufptr = memseta(bufptr, 32, 2);
  }
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
  fputs(allele00, outfile);
  putc_unlocked('/', outfile);
  fputs(allele00, outfile);
  bufptr = memseta(g_textbuf, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[0]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[2]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[3]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[1]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_clipped(((int32_t)marg_a[0]) * tot_recip, bufptr);
  *bufptr++ = '\n';
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);

  bufptr = memseta(g_textbuf, 32, plink_maxsnp + 1);
  if (alen00 + alen01 < 4) {
    bufptr = memseta(bufptr, 32, 4 - alen00 - alen01);
  }
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
  fputs(allele00, outfile);
  putc_unlocked('/', outfile);
  fputs(allele01, outfile);
  bufptr = memseta(g_textbuf, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[8]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[10]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[11]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[9]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_clipped(((int32_t)marg_a[2]) * tot_recip, bufptr);
  *bufptr++ = '\n';
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);

  bufptr = memseta(g_textbuf, 32, plink_maxsnp + 1);
  if (alen01 == 1) {
    bufptr = memseta(bufptr, 32, 2);
  }
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
  fputs(allele01, outfile);
  putc_unlocked('/', outfile);
  fputs(allele01, outfile);
  bufptr = memseta(g_textbuf, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[12]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[14]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[15]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[13]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_clipped(((int32_t)marg_a[3]) * tot_recip, bufptr);
  *bufptr++ = '\n';
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);

  bufptr = memseta(g_textbuf, 32, plink_maxsnp + 3);
  bufptr = memcpya(bufptr, "0/0  ", 5);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[4]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[6]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[7]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)counts[5]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_clipped(((int32_t)marg_a[1]) * tot_recip, bufptr);
  *bufptr++ = '\n';
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);

  bufptr = memseta(g_textbuf, 32, plink_maxsnp + 3);
  bufptr = memcpya(bufptr, "*/*  ", 5);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)marg_b[0]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)marg_b[2]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)marg_b[3]) * tot_recip, bufptr);
  bufptr = memseta(bufptr, 32, 2);
  bufptr = dtoa_f_w9p6_spaced(((int32_t)marg_b[1]) * tot_recip, bufptr);
  bufptr = memcpya(bufptr, "   1\n\n", 6);
  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
}

int32_t twolocus(Epi_info* epi_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uint32_t pheno_ctrl_ct, uintptr_t* pheno_c, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* mkr1 = outname? epi_ip->twolocus_mkr1 : epi_ip->ld_mkr1;
  char* mkr2 = outname? epi_ip->twolocus_mkr2 : epi_ip->ld_mkr2;
  uintptr_t* sample_include2 = nullptr;
  uintptr_t* sample_male_include2 = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t ulii = strlen(mkr1) + 1;
  uintptr_t uljj = strlen(mkr2) + 1;
  uint32_t hwe_midp = epi_ip->modifier & EPI_HWE_MIDP;
  int32_t retval = 0;
  uint32_t counts_all[16];
  uint32_t counts_cc[32];
  uintptr_t* loadbufs[2];
  uintptr_t marker_uidxs[2];
  double solutions[3];
  uint32_t is_haploid[2];
  uint32_t is_x[2];
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf0_ptr;
  uintptr_t* loadbuf1_ptr;
  uintptr_t* loadbuf0_end;
  char* bufptr;
  char* bufptr2;
  uintptr_t sample_ctl2;
  uintptr_t final_mask;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t sample_idx_end;
  uintptr_t ulkk;
  double twice_tot_recip;
  double half_hethet_share;
  double freq11;
  double freq12;
  double freq21;
  double freq22;
  double freq1x;
  double freq2x;
  double freqx1;
  double freqx2;
  double dxx;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t is_y;
  uint32_t alen00;
  uint32_t alen01;
  uint32_t alen10;
  uint32_t alen11;
  uint32_t count_total;
  if (!outname) {
    ulkk = BITCT_TO_WORDCT(unfiltered_sample_ct);
    // ulkk = (unfiltered_sample_ctl2 + 1) / 2;
    sample_ct = popcount_longs(sample_exclude, ulkk);
    if (!sample_ct) {
      logerrprint("Warning: Skipping --ld since there are no founders.  (--make-founders may come\nin handy here.)\n");
      goto twolocus_ret_1;
    }
    if (bigstack_alloc_ul(ulkk, &loadbuf_raw)) {
      goto twolocus_ret_NOMEM;
    }
    bitarr_invert_copy(sample_exclude, unfiltered_sample_ct, loadbuf_raw);
    sample_exclude = loadbuf_raw;
  }
  sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  final_mask = get_final_mask(sample_ct);
  if ((ulii > max_marker_id_len) || (uljj > max_marker_id_len)) {
    goto twolocus_ret_MARKER_NOT_FOUND;
  }
  marker_uidxs[0] = 0;
  marker_uidxs[1] = 0;
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
    bufptr = &(marker_ids[marker_uidx * max_marker_id_len]);
    if (ulii && (!memcmp(mkr1, bufptr, ulii))) {
      marker_uidxs[0] = marker_uidx;
      if (!uljj) {
	break;
      }
      ulii = 0;
    } else if (uljj && (!memcmp(mkr2, bufptr, uljj))) {
      marker_uidxs[1] = marker_uidx;
      if (!ulii) {
	break;
      }
      uljj = 0;
    }
  }
  if (marker_idx == marker_ct) {
    goto twolocus_ret_MARKER_NOT_FOUND;
  }
  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(sample_ctl2, &loadbufs[0]) ||
      bigstack_alloc_ul(sample_ctl2, &loadbufs[1])) {
    goto twolocus_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  loadbufs[0][sample_ctl2 - 1] = 0;
  loadbufs[1][sample_ctl2 - 1] = 0;
  if (alloc_collapsed_haploid_filters(sample_exclude, sex_male, unfiltered_sample_ct, sample_ct, hh_exists, 0, &sample_include2, &sample_male_include2)) {
    goto twolocus_ret_NOMEM;
  }
  is_haploid[0] = 0;
  is_haploid[1] = 0;
  is_x[0] = 0;
  is_x[1] = 0;
  for (marker_idx = 0; marker_idx < 2; marker_idx++) {
    marker_uidx = marker_uidxs[marker_idx];
    if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto twolocus_ret_READ_FAIL;
    }
    if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, loadbufs[marker_idx])) {
      goto twolocus_ret_READ_FAIL;
    }
    chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx);
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    is_haploid[marker_idx] = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
    if (is_haploid[marker_idx]) {
      is_x[marker_idx] = (chrom_idx == (uint32_t)chrom_info_ptr->xymt_codes[X_OFFSET]);
      is_y = (chrom_idx == (uint32_t)chrom_info_ptr->xymt_codes[Y_OFFSET]);
      haploid_fix(hh_exists, sample_include2, sample_male_include2, sample_ct, is_x[marker_idx], is_y, (unsigned char*)(loadbufs[marker_idx]));
    }
  }
  if (!outname) {
    // --ld needs X chromosome sex stratification instead of --twolocus's
    // case/control stratification
    if (is_x[0] || is_x[1]) {
      pheno_c = sex_male;
      pheno_nm = sex_male;
    } else {
      pheno_c = nullptr;
    }
  }
  fill_uint_zero(16, counts_all);
  fill_uint_zero(32, counts_cc);
  loadbuf0_ptr = loadbufs[0];
  loadbuf1_ptr = loadbufs[1];
  loadbuf0_end = &(loadbuf0_ptr[sample_ct / BITCT2]);
  sample_uidx = 0;
  sample_idx = 0;
  sample_idx_end = BITCT2;
  while (1) {
    while (loadbuf0_ptr < loadbuf0_end) {
      ulii = *loadbuf0_ptr++;
      uljj = *loadbuf1_ptr++;
      if (pheno_c) {
	for (; sample_idx < sample_idx_end; sample_uidx++, sample_idx++) {
          next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
	  ulkk = ((ulii & 3) << 2) | (uljj & 3);
	  ulii >>= 2;
	  uljj >>= 2;
	  counts_all[ulkk] += 1;
          if (IS_SET(pheno_nm, sample_uidx)) {
            counts_cc[(16 * IS_SET(pheno_c, sample_uidx)) + ulkk] += 1;
	  }
	}
      } else {
	for (; sample_idx < sample_idx_end; sample_idx++) {
	  ulkk = ((ulii & 3) << 2) | (uljj & 3);
	  ulii >>= 2;
	  uljj >>= 2;
	  counts_all[ulkk] += 1;
	}
      }
      sample_idx_end += BITCT2;
    }
    if (sample_idx == sample_ct) {
      break;
    }
    loadbuf0_end++;
    sample_idx_end = sample_ct;
  }

  alen00 = strlen(marker_allele_ptrs[2 * marker_uidxs[0]]);
  alen01 = strlen(marker_allele_ptrs[2 * marker_uidxs[0] + 1]);
  alen10 = strlen(marker_allele_ptrs[2 * marker_uidxs[1]]);
  alen11 = strlen(marker_allele_ptrs[2 * marker_uidxs[1] + 1]);
  if (outname) {
    memcpy(outname_end, ".twolocus", 10);
    if (fopen_checked(outname, "w", &outfile)) {
      goto twolocus_ret_OPEN_FAIL;
    }
    fputs("\nAll individuals\n===============\n", outfile);
    twolocus_write_table(outfile, counts_all, plink_maxsnp, mkr1, mkr2, marker_allele_ptrs[2 * marker_uidxs[0]], marker_allele_ptrs[2 * marker_uidxs[0] + 1], marker_allele_ptrs[2 * marker_uidxs[1]], marker_allele_ptrs[2 * marker_uidxs[1] + 1], alen00, alen01, alen10, alen11);
    if (pheno_c) {
      if (pheno_nm_ct != pheno_ctrl_ct) {
	fputs("\nCases\n=====\n", outfile);
	twolocus_write_table(outfile, &(counts_cc[16]), plink_maxsnp, mkr1, mkr2, marker_allele_ptrs[2 * marker_uidxs[0]], marker_allele_ptrs[2 * marker_uidxs[0] + 1], marker_allele_ptrs[2 * marker_uidxs[1]], marker_allele_ptrs[2 * marker_uidxs[1] + 1], alen00, alen01, alen10, alen11);
      }
      if (pheno_ctrl_ct) {
	fputs("\nControls\n========\n", outfile);
	twolocus_write_table(outfile, counts_cc, plink_maxsnp, mkr1, mkr2, marker_allele_ptrs[2 * marker_uidxs[0]], marker_allele_ptrs[2 * marker_uidxs[0] + 1], marker_allele_ptrs[2 * marker_uidxs[1]], marker_allele_ptrs[2 * marker_uidxs[1] + 1], alen00, alen01, alen10, alen11);
      }
    }
    putc_unlocked('\n', outfile);
    if (fclose_null(&outfile)) {
      goto twolocus_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("--twolocus: Report written to %s .\n", outname);
  } else {
    // low counts_cc[] values aren't used, so may as well store marginal counts
    // there
    counts_cc[0] = counts_all[0] + counts_all[2] + counts_all[3];
    counts_cc[2] = counts_all[8] + counts_all[10] + counts_all[11];
    counts_cc[3] = counts_all[12] + counts_all[14] + counts_all[15];
    counts_cc[4] = counts_all[0] + counts_all[8] + counts_all[12];
    counts_cc[6] = counts_all[2] + counts_all[10] + counts_all[14];
    counts_cc[7] = counts_all[3] + counts_all[11] + counts_all[15];
    count_total = counts_cc[0] + counts_cc[2] + counts_cc[3];
    if (!count_total) {
      logerrprint("Error: No valid observations for --ld.\n");
      goto twolocus_ret_INVALID_CMDLINE;
    }
    if ((!counts_cc[2]) && ((!counts_cc[0]) || (!counts_cc[3]))) {
      LOGPREPRINTFWW("Error: %s is monomorphic across all valid observations.\n", mkr1);
      goto twolocus_ret_INVALID_CMDLINE_2;
    } else if ((!counts_cc[6]) && ((!counts_cc[4]) || (!counts_cc[7]))) {
      LOGPREPRINTFWW("Error: %s is monomorphic across all valid observations.\n", mkr2);
      goto twolocus_ret_INVALID_CMDLINE_2;
    } else if ((alen00 > (MAXLINELEN / 4) - 16) || (alen01 > (MAXLINELEN / 4) - 16)) {
      LOGPREPRINTFWW("Error: %s has a pathologically long allele code.\n", mkr1);
      goto twolocus_ret_INVALID_CMDLINE_2;
    } else if ((alen10 > (MAXLINELEN / 4) - 16) || (alen11 > (MAXLINELEN / 4) - 16)) {
      LOGPREPRINTFWW("Error: %s has a pathologically long allele code.\n", mkr2);
      goto twolocus_ret_INVALID_CMDLINE_2;
    }
    LOGPRINTF("\n--ld %s %s:\n", mkr1, mkr2);
    ulii = 0;
    // possible todo: factor out redundancy with other D-prime calculations
    freq11 = (double)(2 * counts_all[0] + counts_all[2] + counts_all[8]);
    freq12 = (double)(2 * counts_all[3] + counts_all[2] + counts_all[11]);
    freq21 = (double)(2 * counts_all[12] + counts_all[8] + counts_all[14]);
    freq22 = (double)(2 * counts_all[15] + counts_all[11] + counts_all[14]);
    if (is_x[0] || is_x[1]) {
      if (is_x[0] && is_x[1]) {
        freq11 -= (double)((int32_t)counts_cc[16]);
        freq12 -= (double)((int32_t)counts_cc[19]);
        freq21 -= (double)((int32_t)counts_cc[28]);
        freq22 -= (double)((int32_t)counts_cc[31]);
      } else if (is_x[0]) {
        freq11 -= ((double)(2 * counts_cc[16] + counts_cc[18])) * (1.0 - SQRT_HALF);
        freq12 -= ((double)(2 * counts_cc[19] + counts_cc[18])) * (1.0 - SQRT_HALF);
        freq21 -= ((double)(2 * counts_cc[28] + counts_cc[30])) * (1.0 - SQRT_HALF);
        freq22 -= ((double)(2 * counts_cc[31] + counts_cc[30])) * (1.0 - SQRT_HALF);
      } else {
        freq11 -= ((double)(2 * counts_cc[16] + counts_cc[24])) * (1.0 - SQRT_HALF);
        freq12 -= ((double)(2 * counts_cc[19] + counts_cc[27])) * (1.0 - SQRT_HALF);
        freq21 -= ((double)(2 * counts_cc[28] + counts_cc[24])) * (1.0 - SQRT_HALF);
        freq22 -= ((double)(2 * counts_cc[31] + counts_cc[27])) * (1.0 - SQRT_HALF);
      }
    }
    twice_tot_recip = 1.0 / (freq11 + freq12 + freq21 + freq22 + 2 * ((int32_t)counts_all[10]));
    freq11 *= twice_tot_recip;
    freq12 *= twice_tot_recip;
    freq21 *= twice_tot_recip;
    freq22 *= twice_tot_recip;
    half_hethet_share = ((int32_t)counts_all[10]) * twice_tot_recip;
    freq1x = freq11 + freq12 + half_hethet_share;
    freq2x = 1.0 - freq1x;
    freqx1 = freq11 + freq21 + half_hethet_share;
    freqx2 = 1.0 - freqx1;
    if (counts_all[10]) {
      // detect degenerate cases to avoid e-17 ugliness
      // possible todo: when there are multiple solutions, compute log
      // likelihood for each and mark the EM solution in some manner
      if ((freq11 * freq22 != 0.0) || (freq12 * freq21 != 0.0)) {
	// (f11 + x)(f22 + x)(K - x) = x(f12 + K - x)(f21 + K - x)
	// (x - K)(x + f11)(x + f22) + x(x - K - f12)(x - K - f21) = 0
	//   x^3 + (f11 + f22 - K)x^2 + (f11*f22 - K*f11 - K*f22)x
	// - K*f11*f22 + x^3 - (2K + f12 + f21)x^2 + (K + f12)(K + f21)x = 0
	uljj = cubic_real_roots(0.5 * (freq11 + freq22 - freq12 - freq21 - 3 * half_hethet_share), 0.5 * (freq11 * freq22 + freq12 * freq21 + half_hethet_share * (freq12 + freq21 - freq11 - freq22 + half_hethet_share)), -0.5 * half_hethet_share * freq11 * freq22, solutions);
	if (uljj > 1) {
	  while (solutions[uljj - 1] > half_hethet_share + SMALLISH_EPSILON) {
	    uljj--;
	  }
	  if (solutions[uljj - 1] > half_hethet_share - SMALLISH_EPSILON) {
	    solutions[uljj - 1] = half_hethet_share;
	  }
	  while (solutions[ulii] < -SMALLISH_EPSILON) {
	    ulii++;
	  }
	  if (solutions[ulii] < SMALLISH_EPSILON) {
	    solutions[ulii] = 0;
	  }
	}
      } else {
        // bugfix (6 Oct 2017):
        // At least one of {f11, f22} is zero, and one of {f12, f21} is zero.
        // Initially suppose that the zero-values are f11 and f12.  Then the
        // equality becomes
        //   x(f22 + x)(K - x) = x(K - x)(f21 + K - x)
        //   x=0 and x=K are always solutions; the rest becomes
        //     f22 + x = f21 + K - x
        //     2x = K + f21 - f22
        //     x = (K + f21 - f22)/2; in-range iff (f21 - f22) in (-K, K).
        // So far so good.  However, this code used to *always* check
        // (f21 - f22), when it's necessary to use all the nonzero values.
        // (this still works if three or all four values are zero)
	solutions[0] = 0;
        const double nonzero_freq_xx = freq11 + freq22;
        const double nonzero_freq_xy = freq12 + freq21;
	if ((nonzero_freq_xx + SMALLISH_EPSILON < half_hethet_share + nonzero_freq_xy) && (nonzero_freq_xy + SMALLISH_EPSILON < half_hethet_share + nonzero_freq_xx)) {
	  uljj = 3;
	  solutions[1] = (half_hethet_share + nonzero_freq_xy - nonzero_freq_xx) * 0.5;
	  solutions[2] = half_hethet_share;
	} else {
	  uljj = 2;
	  solutions[1] = half_hethet_share;
	}
      }
      if (uljj > ulii + 1) {
	// not Xchr/haploid-sensitive yet
	logprint("Multiple haplotype phasing solutions; sample size, HWE, or random mating\nassumption may be violated.\n\nHWE exact test p-values\n-----------------------\n");
	if (is_haploid[0] && (!is_x[0])) {
          LOGPRINTF("   %s: n/a\n", mkr1);
	} else {
	  LOGPRINTF("   %s: %g\n", mkr1, SNPHWE2(counts_cc[2] + counts_all[9], counts_cc[0] + counts_all[1] - 2 * (counts_cc[16] + counts_cc[19]), counts_cc[3] + counts_all[13] - 2 * (counts_cc[28] + counts_cc[31]), hwe_midp));
	}
	if (is_haploid[1] && (!is_x[1])) {
	  LOGPRINTF("   %s: n/a\n", mkr2);
	} else {
	  LOGPRINTF("   %s: %g\n\n", mkr2, SNPHWE2(counts_cc[6] + counts_all[6], counts_cc[4] + counts_all[4] - 2 * (counts_cc[16] + counts_cc[28]), counts_cc[7] + counts_all[7] - 2 * (counts_cc[19] + counts_cc[31]), hwe_midp));
	}
      }
    } else {
      uljj = 1;
      solutions[0] = 0.0;
    }
    if (uljj == ulii + 1) {
      logprint("\n");
    }
    for (ulkk = ulii; ulkk < uljj; ulkk++) {
      if (uljj - ulii > 1) {
        LOGPRINTF("Solution #%" PRIuPTR ":\n", ulkk + 1 - ulii);
      }
      dxx = freq11 + solutions[ulkk] - freqx1 * freq1x; // D
      if (fabs(dxx) < SMALL_EPSILON) {
	dxx = 0;
      }
      bufptr = memcpya(g_logbuf, "   R-sq = ", 10);
      bufptr2 = dtoa_g(dxx * dxx / (freq1x * freqx1 * freq2x * freqx2), bufptr);
      // assumes bufptr2 - bufptr < 15
      bufptr = memseta(bufptr2, 32, 15 - ((uintptr_t)(bufptr2 - bufptr)));
      bufptr = memcpya(bufptr, "D' = ", 5);
      if (dxx >= 0) {
	bufptr = dtoa_g(dxx / MINV(freqx1 * freq2x, freqx2 * freq1x), bufptr);
      } else {
	bufptr = dtoa_g(-dxx / MINV(freqx1 * freq1x, freqx2 * freq2x), bufptr);
      }
      bufptr = memcpya(bufptr, "\n\n", 3);
      logprintb();
      logprint("   Haplotype     Frequency    Expectation under LE\n");
      logprint("   ---------     ---------    --------------------\n");
      bufptr = memseta(g_logbuf, 32, 3);
      if (alen00 + alen10 < 9) {
	bufptr = memseta(bufptr, 32, 9 - alen00 - alen10);
      }
      bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[0]], alen00);
      bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[1]], alen10);
      bufptr = memseta(bufptr, 32, 5);
      bufptr = dtoa_f_w9p6_spaced(freq11 + solutions[ulkk], bufptr);
      bufptr = memseta(bufptr, 32, 15);
      bufptr = dtoa_f_w9p6_clipped(freqx1 * freq1x, bufptr);
      bufptr = memcpya(bufptr, "\n", 2);
      logprintb();
      bufptr = &(g_logbuf[3]);
      if (alen01 + alen10 < 9) {
	bufptr = memseta(bufptr, 32, 9 - alen01 - alen10);
      }
      bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[0] + 1], alen01);
      bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[1]], alen10);
      bufptr = memseta(bufptr, 32, 5);
      bufptr = dtoa_f_w9p6_spaced(freq21 + half_hethet_share - solutions[ulkk], bufptr);
      bufptr = memseta(bufptr, 32, 15);
      bufptr = dtoa_f_w9p6_clipped(freqx1 * freq2x, bufptr);
      bufptr = memcpya(bufptr, "\n", 2);
      logprintb();
      bufptr = &(g_logbuf[3]);
      if (alen00 + alen11 < 9) {
	bufptr = memseta(bufptr, 32, 9 - alen00 - alen11);
      }
      bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[0]], alen00);
      bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[1] + 1], alen11);
      bufptr = memseta(bufptr, 32, 5);
      bufptr = dtoa_f_w9p6_spaced(freq12 + half_hethet_share - solutions[ulkk], bufptr);
      bufptr = memseta(bufptr, 32, 15);
      bufptr = dtoa_f_w9p6_clipped(freqx2 * freq1x, bufptr);
      bufptr = memcpya(bufptr, "\n", 2);
      logprintb();
      bufptr = &(g_logbuf[3]);
      if (alen01 + alen11 < 9) {
	bufptr = memseta(bufptr, 32, 9 - alen01 - alen11);
      }
      bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[0] + 1], alen01);
      bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[1] + 1], alen11);
      bufptr = memseta(bufptr, 32, 5);
      bufptr = dtoa_f_w9p6_spaced(freq22 + solutions[ulkk], bufptr);
      bufptr = memseta(bufptr, 32, 15);
      bufptr = dtoa_f_w9p6_clipped(freqx2 * freq2x, bufptr);
      bufptr = memcpyl3a(bufptr, "\n\n");
      logprintb();
      bufptr = &(g_logbuf[3]);
      bufptr = memcpya(bufptr, "In phase alleles are ", 21);
      if (dxx > 0) {
	bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[0]], alen00);
	bufptr = memcpyax(bufptr, marker_allele_ptrs[2 * marker_uidxs[1]], alen10, '/');
	bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[0] + 1], alen01);
	bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[1] + 1], alen11);
      } else {
	bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[0]], alen00);
	bufptr = memcpyax(bufptr, marker_allele_ptrs[2 * marker_uidxs[1] + 1], alen11, '/');
	bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[0] + 1], alen01);
	bufptr = memcpya(bufptr, marker_allele_ptrs[2 * marker_uidxs[1]], alen10);
      }
      bufptr = memcpyl3a(bufptr, "\n\n");
      logprintb();
    }
  }
  while (0) {
  twolocus_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  twolocus_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  twolocus_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  twolocus_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  twolocus_ret_MARKER_NOT_FOUND:
    if (outname) {
      logerrprint("Error: --twolocus variant name not found.\n");
    } else {
      logerrprint("Error: --ld variant name not found.\n");
    }
    retval = RET_INVALID_CMDLINE;
    break;
  twolocus_ret_INVALID_CMDLINE_2:
    logerrprintb();
  twolocus_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 twolocus_ret_1:
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

void rotate_loadbuf_and_compute_phenogeno(uintptr_t* loadbuf, double* pheno_d2, uint32_t pheno_nm_ct, uintptr_t* loadbuf_write, double* phenogeno, uint32_t* genosums) {
  double cur_phenogeno = 0;
  uint32_t geno1_ct = 0;
  uint32_t geno2_ct = 0;
  uintptr_t cur_word;
  uintptr_t ulii;
  uintptr_t uljj;
  double dxx;
  uint32_t sample_idx;
  uint32_t sample_idx_base;
  uint32_t uii;
  for (sample_idx = 0; sample_idx < pheno_nm_ct;) {
    // we're interested in hom A1 (non-trailing 00) and het (10) bit
    // values here.
    cur_word = ~(*loadbuf++);
    sample_idx_base = sample_idx;
    sample_idx += BITCT2;
    if (sample_idx > pheno_nm_ct) {
      cur_word &= (ONELU << (2 * (pheno_nm_ct % BITCT2))) - ONELU;
    }
    // now hom A1 = 11 and het = 01.  Temporarily erase the 10s.
    uljj = cur_word & FIVEMASK;
    ulii = uljj | (cur_word & (uljj << 1));
    while (ulii) {
      uii = CTZLU(ulii);
      dxx = pheno_d2[sample_idx_base + uii / 2];
      if (!((ulii >> (uii + 1)) & 1)) {
	// het
	cur_phenogeno += dxx;
	geno1_ct++;
      } else {
	// hom A1
	cur_phenogeno += 2 * dxx;
	geno2_ct++;
      }
      ulii &= ~((3 * ONELU) << uii);
    }
    // currently hom A1 = 11, missing = 10, het = 01, hom A2 = 00
    // rotate to hom A1 = 10, missing = 11, het = 01, hom A2 = 00
    // to allow inner loop to use ordinary multiplication
    *loadbuf_write++ = cur_word ^ ((cur_word >> 1) & FIVEMASK);
  }
  *phenogeno = cur_phenogeno;
  genosums[0] = geno1_ct + 2 * geno2_ct;
  genosums[1] = geno1_ct + 4 * geno2_ct;
}

int32_t epistasis_linear_regression(pthread_t* threads, Epi_info* epi_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t marker_uidx_base, uintptr_t marker_ct1, uintptr_t* marker_exclude1, uintptr_t marker_idx1_start, uintptr_t marker_idx1_end, uintptr_t marker_ct2, uintptr_t* marker_exclude2, uint32_t is_triangular, uintptr_t job_size, uint64_t tests_expected, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, double* pheno_d, uint32_t parallel_idx, uint32_t parallel_tot, char* outname, char* outname_end, double output_min_p, double glm_vif_thresh, uintptr_t* loadbuf_raw, uintptr_t* loadbuf, double* best_chisq, uint32_t* best_ids, uint32_t* n_sig_cts, uint32_t* fail_cts, uint32_t* gap_cts) {
  // We use QT --assoc's strategy for speeding up linear regression, since we
  // do not need to support arbitrary covariates.  It's more complicated here
  // because we have 3 covariates instead of one, but two of them are still
  // restricted to the values {0, 1, 2} and the last is the product of the
  // first two.  So we're able to use variations of the QT --assoc bit hacks.
  FILE* outfile = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t pheno_nm_ctl2 = QUATERCT_TO_WORDCT(pheno_nm_ct);
  uintptr_t final_mask = get_final_mask(pheno_nm_ct);
  uintptr_t marker_uidx = marker_uidx_base;
  uintptr_t pct = 1;
  uintptr_t marker_uidx2 = 0;
  uintptr_t marker_idx1 = marker_idx1_start;
  uintptr_t marker_idx2 = 0;
  uint64_t pct_thresh = tests_expected / 100;
  uint64_t tests_complete = 0;
  uint32_t max_thread_ct = g_epi_thread_ct;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t chrom_end = 0;
  uint32_t chrom_idx = 0;
  uint32_t chrom_idx2 = 0;
  int32_t retval = 0;
  unsigned char* bigstack_mark2;
  char* wptr_start;
  char* wptr_start2;
  char* wptr;
  double* pheno_d2;
  double* dptr;
  double* dptr2;
  uint32_t* uiptr;
  uint32_t* uiptr2;
  uint32_t* uiptr3;
  uint32_t* uiptr4;
  uint32_t* uiptr5;
  uintptr_t cur_bigstack_left;
  uintptr_t cur_workload;
  uintptr_t idx1_block_size;
  uintptr_t idx2_block_size;
  uintptr_t idx2_block_sizea16;
  uintptr_t marker_uidx_tmp;
  uintptr_t block_idx1;
  uintptr_t block_idx2;
  uintptr_t cur_idx2_block_size;
  uintptr_t chrom_end2;
  uintptr_t tidx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  double dxx;
  uint32_t chrom_fo_idx;
  uint32_t chrom_fo_idx2;
  uint32_t is_last_block;
  uint32_t sample_uidx;
  uint32_t sample_idx;
  uint32_t uii;
  uint32_t ujj;
  if (bigstack_alloc_d(pheno_nm_ct, &pheno_d2)) {
    goto epistasis_linear_regression_ret_NOMEM;
  }
  g_epi_pheno_d2 = pheno_d2;
  g_epi_pheno_nm_ct = pheno_nm_ct;
  dptr = pheno_d2;
  g_epi_pheno_sum = 0;
  g_epi_pheno_ssq = 0;
  for (sample_uidx = 0, sample_idx = 0; sample_idx < pheno_nm_ct; sample_uidx++, sample_idx++) {
    next_set_unsafe_ck(pheno_nm, &sample_uidx);
    dxx = pheno_d[sample_uidx];
    g_epi_pheno_sum += dxx;
    g_epi_pheno_ssq += dxx * dxx;
    pheno_d2[sample_idx] = dxx;
  }
  // could add an epsilon here, but this is good enough to catch the most
  // common case (all phenotypes are the same integer near zero).
  if (g_epi_pheno_ssq * ((double)((int32_t)pheno_nm_ct)) == g_epi_pheno_sum * g_epi_pheno_sum) {
    logerrprint("Error: Phenotype is constant.\n");
    goto epistasis_linear_regression_ret_INVALID_CMDLINE;
  }
  g_epi_vif_thresh = glm_vif_thresh;

  // claim up to half of memory with idx1 bufs; each marker currently costs:
  //   pheno_nm_ctl2 * sizeof(intptr_t) for geno buf
  //   sizeof(double) for precomputed sum(phenotype * genotype) values
  //   2 * sizeof(int32_t) for precomputed sum(genotype) and sum(genotype^2)
  //     values
  //   4 * sizeof(int32_t) + sizeof(double) + marker_ct2 * 2 * sizeof(double)
  //     for other stuff (see epistasis_report() comment, starting from
  //     "offset"; main result buffer must be double-size to store both beta
  //     and chi-square stat)
  cur_bigstack_left = bigstack_left();
  ulii = 6 * CACHELINE + max_thread_ct * (5 * (CACHELINE - 4)) - 5 * sizeof(int32_t) - sizeof(double);
  if (cur_bigstack_left >= ulii) {
    cur_bigstack_left -= ulii;
  }
  ulii = pheno_nm_ctl2 * sizeof(intptr_t) + 6 * sizeof(int32_t) + 2 * sizeof(double) + marker_ct2 * 2 * sizeof(double);
  idx1_block_size = cur_bigstack_left / (ulii * 2 + 1);
  if (!idx1_block_size) {
    goto epistasis_linear_regression_ret_NOMEM;
  }
  if (idx1_block_size > job_size) {
    idx1_block_size = job_size;
  }
  // pad to avoid threads writing to same cacheline
  ulii = (max_thread_ct - 1) * 15 + idx1_block_size;
  bigstack_alloc_ui(idx1_block_size * 2, &g_epi_geno1_offsets);
  bigstack_alloc_ul(pheno_nm_ctl2 * idx1_block_size, &g_epi_geno1);
  bigstack_alloc_d(idx1_block_size, &g_epi_phenogeno1);
  // may be better to just recompute genosums values in inner loop?  can test
  // this later
  bigstack_alloc_ui(idx1_block_size * 2, &g_epi_genosums1);
  bigstack_alloc_d(idx1_block_size * marker_ct2 * 2, &g_epi_all_chisq);
  bigstack_alloc_d(ulii, &g_epi_best_chisq1);
  bigstack_alloc_ui(ulii, &g_epi_best_id1);
  bigstack_alloc_ui(ulii, &g_epi_n_sig_ct1);
  bigstack_alloc_ui(ulii, &g_epi_fail_ct1);
  for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
    g_epi_geno1[block_idx1 * pheno_nm_ctl2 + pheno_nm_ctl2 - 1] = 0;
  }
  if (is_triangular) {
    fill_uint_zero(2 * idx1_block_size, g_epi_geno1_offsets);
  }

  ulii = pheno_nm_ctl2 * sizeof(intptr_t) + 2 * sizeof(int32_t) + sizeof(double) + max_thread_ct * (3 * sizeof(int32_t) + sizeof(double));
  idx2_block_size = (bigstack_left() - (3 * CACHELINE - sizeof(intptr_t) - 2 * sizeof(int32_t) - sizeof(double)) - max_thread_ct * (3 * (CACHELINE - sizeof(int32_t)) + (CACHELINE - sizeof(double)))) / ulii;
  if (idx2_block_size > marker_ct2) {
    idx2_block_size = marker_ct2;
  }
  idx2_block_size = round_up_pow2(idx2_block_size, 16);

  memcpy(outname_end, ".epi.qt", 8);
  if (parallel_tot > 1) {
    outname_end[7] = '.';
    uint32toa_x(parallel_idx + 1, '\0', &(outname_end[8]));
  }
  if (fopen_checked(outname, "w", &outfile)) {
    goto epistasis_linear_regression_ret_OPEN_FAIL;
  }
  if (!parallel_idx) {
    wptr = memcpya(g_textbuf, "CHR1 ", 5);
    wptr = fw_strcpyn(plink_maxsnp, 4, "SNP1", wptr);
    wptr = memcpya(wptr, " CHR2 ", 6);
    wptr = fw_strcpyn(plink_maxsnp, 4, "SNP2", wptr);
    wptr = memcpya(wptr, "     BETA_INT         STAT            P \n", 41);
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
      goto epistasis_linear_regression_ret_WRITE_FAIL;
    }
  }

  bigstack_mark2 = g_bigstack_base;
  while (1) {
    if (!idx2_block_size) {
      goto epistasis_linear_regression_ret_NOMEM;
    }
    if (!(bigstack_alloc_ul(pheno_nm_ctl2 * idx2_block_size, &g_epi_geno2) ||
          bigstack_alloc_d(idx2_block_size, &g_epi_phenogeno2) ||
          bigstack_alloc_ui(idx2_block_size * 2, &g_epi_genosums2) ||
          bigstack_alloc_d(max_thread_ct * idx2_block_size, &g_epi_best_chisq2) ||
          bigstack_alloc_ui(max_thread_ct * idx2_block_size, &g_epi_best_id2) ||
          bigstack_alloc_ui(max_thread_ct * idx2_block_size, &g_epi_n_sig_ct2) ||
          bigstack_alloc_ui(max_thread_ct * idx2_block_size, &g_epi_fail_ct2))) {
      break;
    }
    bigstack_reset(bigstack_mark2);
    idx2_block_size -= 16;
  }
  for (block_idx2 = 0; block_idx2 < idx2_block_size; block_idx2++) {
    g_epi_geno2[block_idx2 * pheno_nm_ctl2 + pheno_nm_ctl2 - 1] = 0;
  }
  marker_uidx = next_unset_ul_unsafe(marker_exclude1, marker_uidx_base);
  if (marker_idx1) {
    marker_uidx = jump_forward_unset_unsafe(marker_exclude1, marker_uidx + 1, marker_idx1);
  }
  wptr = memcpya(g_logbuf, "QT --epistasis to ", 18);
  wptr = strcpya(wptr, outname);
  memcpy(wptr, " ... ", 6);
  wordwrapb(16); // strlen("99% [processing]")
  logprintb();
  fputs("0%", stdout);
  do {
    fputs(" [processing]", stdout);
    fflush(stdout);
    if (idx1_block_size > marker_idx1_end - marker_idx1) {
      idx1_block_size = marker_idx1_end - marker_idx1;
      if (idx1_block_size < max_thread_ct) {
        max_thread_ct = idx1_block_size;
        g_epi_thread_ct = max_thread_ct;
      }
    }
    g_epi_marker_idx1 = marker_idx1;
    dptr = g_epi_all_chisq;
    dptr2 = &(g_epi_all_chisq[idx1_block_size * marker_ct2 * 2]);
    do {
      *dptr = -1;
      dptr = &(dptr[2]);
    } while (dptr < dptr2);
    marker_uidx_tmp = marker_uidx;
    if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto epistasis_linear_regression_ret_READ_FAIL;
    }
    cur_workload = idx1_block_size * marker_ct2;
    if (is_triangular) {
      for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
        ulii = block_idx1 + marker_idx1 + 1;
        cur_workload -= ulii;
        g_epi_geno1_offsets[2 * block_idx1 + 1] = ulii;
      }
    } else {
      fill_uint_zero(2 * idx1_block_size, g_epi_geno1_offsets);
      marker_uidx2 = marker_uidx_base;
      marker_idx2 = 0;
    }
    tests_complete += cur_workload;
    ulii = 0; // total number of tests
    g_epi_idx1_block_bounds[0] = 0;
    g_epi_idx1_block_bounds16[0] = 0;
    block_idx1 = 0;
    for (tidx = 1; tidx < max_thread_ct; tidx++) {
      uljj = (((uint64_t)cur_workload) * tidx) / max_thread_ct;
      if (is_triangular) {
        do {
          ulii += marker_ct2 - g_epi_geno1_offsets[2 * block_idx1 + 1];
          block_idx1++;
	} while (ulii < uljj);
      } else {
        do {
	  ulii += marker_ct2;
	  block_idx1++;
	} while (ulii < uljj);
      }
      uii = block_idx1 - g_epi_idx1_block_bounds[tidx - 1];
      g_epi_idx1_block_bounds[tidx] = block_idx1;
      g_epi_idx1_block_bounds16[tidx] = g_epi_idx1_block_bounds16[tidx - 1] + round_up_pow2_ui(uii, 16);
    }
    g_epi_idx1_block_bounds[max_thread_ct] = idx1_block_size;
    chrom_end = 0;
    for (block_idx1 = 0; block_idx1 < idx1_block_size; marker_uidx_tmp++, block_idx1++) {
      if (IS_SET(marker_exclude1, marker_uidx_tmp)) {
	marker_uidx_tmp = next_unset_ul_unsafe(marker_exclude1, marker_uidx_tmp);
        if (fseeko(bedfile, bed_offset + (marker_uidx_tmp * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
          goto epistasis_linear_regression_ret_READ_FAIL;
	}
      }
      if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, IS_SET(marker_reverse, marker_uidx_tmp), bedfile, loadbuf_raw, loadbuf)) {
        goto epistasis_linear_regression_ret_READ_FAIL;
      }
      rotate_loadbuf_and_compute_phenogeno(loadbuf, pheno_d2, pheno_nm_ct, &(g_epi_geno1[block_idx1 * pheno_nm_ctl2]), &(g_epi_phenogeno1[block_idx1]), &(g_epi_genosums1[block_idx1 * 2]));
      if (!is_triangular) {
	if (!IS_SET(marker_exclude2, marker_uidx_tmp)) {
          // do not compare against self
          marker_idx2 += marker_uidx_tmp - marker_uidx2 - popcount_bit_idx(marker_exclude2, marker_uidx2, marker_uidx_tmp);
          marker_uidx2 = marker_uidx_tmp;
          g_epi_geno1_offsets[2 * block_idx1] = marker_idx2;
          g_epi_geno1_offsets[2 * block_idx1 + 1] = marker_idx2 + 1;
          gap_cts[block_idx1 + marker_idx1] = 1;
	}
      }
    }
    marker_uidx2 = next_unset_ul_unsafe(marker_exclude2, marker_uidx_base);
    if (is_triangular) {
      marker_idx2 = marker_idx1 + 1;
      marker_uidx2 = jump_forward_unset_unsafe(marker_exclude2, marker_uidx2 + 1, marker_idx2);
    } else {
      marker_idx2 = 0;
    }
    if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto epistasis_linear_regression_ret_READ_FAIL;
    }
    cur_idx2_block_size = idx2_block_size;
    do {
      if (cur_idx2_block_size > marker_ct2 - marker_idx2) {
        cur_idx2_block_size = marker_ct2 - marker_idx2;
      }
      for (block_idx2 = 0; block_idx2 < cur_idx2_block_size; marker_uidx2++, block_idx2++) {
        if (IS_SET(marker_exclude2, marker_uidx2)) {
          marker_uidx2 = next_unset_ul_unsafe(marker_exclude2, marker_uidx2);
          if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
            goto epistasis_linear_regression_ret_READ_FAIL;
	  }
	}
	if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, IS_SET(marker_reverse, marker_uidx2), bedfile, loadbuf_raw, loadbuf)) {
	  goto epistasis_linear_regression_ret_READ_FAIL;
	}
        rotate_loadbuf_and_compute_phenogeno(loadbuf, pheno_d2, pheno_nm_ct, &(g_epi_geno2[block_idx2 * pheno_nm_ctl2]), &(g_epi_phenogeno2[block_idx2]), &(g_epi_genosums2[block_idx2 * 2]));
      }
      g_epi_idx2_block_size = cur_idx2_block_size;
      g_epi_idx2_block_start = marker_idx2;
      idx2_block_sizea16 = round_up_pow2(cur_idx2_block_size, 16);
      fill_uint_zero(idx1_block_size + 15 * (max_thread_ct - 1), g_epi_n_sig_ct1);
      fill_uint_zero(idx1_block_size + 15 * (max_thread_ct - 1), g_epi_fail_ct1);
      fill_uint_zero(idx2_block_sizea16 * max_thread_ct, g_epi_n_sig_ct2);
      fill_uint_zero(idx2_block_sizea16 * max_thread_ct, g_epi_fail_ct2);
      for (tidx = 0; tidx < max_thread_ct; tidx++) {
        ulii = g_epi_idx1_block_bounds[tidx];
        uljj = g_epi_idx1_block_bounds[tidx + 1] - ulii;
        dptr = &(g_epi_best_chisq1[g_epi_idx1_block_bounds[tidx]]);
	dptr2 = &(g_epi_all_chisq[(marker_idx1 + ulii) * 2]);
	for (ulkk = 0; ulkk < uljj; ulkk++) {
          *dptr++ = dptr2[ulkk * 2];
	}
        ulii = g_epi_geno1_offsets[2 * ulii + 1];
        if (ulii < marker_idx2 + cur_idx2_block_size) {
          if (ulii <= marker_idx2) {
            ulii = 0;
	  } else {
            ulii -= marker_idx2;
	  }
          uljj = cur_idx2_block_size - ulii;
	  dptr = &(g_epi_best_chisq2[tidx * idx2_block_sizea16 + ulii]);
	  dptr2 = &(g_epi_all_chisq[(marker_idx2 + ulii) * 2]);
          for (ulkk = 0; ulkk < uljj; ulkk++) {
            *dptr++ = dptr2[ulkk * 2];
	  }
	}
      }
      is_last_block = (marker_idx2 + cur_idx2_block_size >= marker_ct2);
      if (spawn_threads2(threads, &epi_linear_thread, max_thread_ct, is_last_block)) {
	goto epistasis_linear_regression_ret_THREAD_CREATE_FAIL;
      }
      epi_linear_thread((void*)0);
      join_threads2(threads, max_thread_ct, is_last_block);
      // merge best_chisq, best_ids, fail_cts
      for (tidx = 0; tidx < max_thread_ct; tidx++) {
	ulii = g_epi_idx1_block_bounds[tidx];
	uljj = g_epi_idx1_block_bounds[tidx + 1] - ulii;
	uii = g_epi_idx1_block_bounds16[tidx];
	dptr = &(g_epi_best_chisq1[uii]);
	uiptr = &(g_epi_best_id1[uii]);
	uiptr2 = &(g_epi_n_sig_ct1[uii]);
	uiptr3 = &(g_epi_fail_ct1[uii]);
	ulii += marker_idx1;
	dptr2 = &(best_chisq[ulii]);
	uiptr4 = &(n_sig_cts[ulii]);
	uiptr5 = &(fail_cts[ulii]);
	for (block_idx1 = 0; block_idx1 < uljj; block_idx1++, dptr2++, uiptr4++, uiptr5++) {
	  dxx = *dptr++;
	  if (dxx > (*dptr2)) {
	    *dptr2 = dxx;
	    best_ids[block_idx1 + ulii] = uiptr[block_idx1];
	  }
	  *uiptr4 += *uiptr2++;
	  *uiptr5 += *uiptr3++;
	}
      }
      if (is_triangular) {
	for (tidx = 0; tidx < max_thread_ct; tidx++) {
	  block_idx2 = g_epi_geno1_offsets[2 * g_epi_idx1_block_bounds[tidx] + 1];
	  if (block_idx2 <= marker_idx2) {
	    block_idx2 = 0;
	  } else {
	    block_idx2 -= marker_idx2;
	  }
	  dptr = &(g_epi_best_chisq2[tidx * idx2_block_sizea16 + block_idx2]);
	  uiptr = &(g_epi_best_id2[tidx * idx2_block_sizea16]);
	  uiptr2 = &(g_epi_n_sig_ct2[tidx * idx2_block_sizea16 + block_idx2]);
	  uiptr3 = &(g_epi_fail_ct2[tidx * idx2_block_sizea16 + block_idx2]);
	  dptr2 = &(best_chisq[block_idx2 + marker_idx2]);
	  uiptr4 = &(n_sig_cts[block_idx2 + marker_idx2]);
	  uiptr5 = &(fail_cts[block_idx2 + marker_idx2]);
	  for (; block_idx2 < cur_idx2_block_size; block_idx2++, dptr2++, uiptr4++, uiptr5++) {
	    dxx = *dptr++;
	    if (dxx > (*dptr2)) {
	      *dptr2 = dxx;
	      best_ids[block_idx2 + marker_idx2] = uiptr[block_idx2];
	    }
	    *uiptr4 += *uiptr2++;
	    *uiptr5 += *uiptr3++;
	  }
	}
      }
      marker_idx2 += cur_idx2_block_size;
    } while (marker_idx2 < marker_ct2);
    fputs("\b\b\b\b\b\b\b\b\b\b\bwriting]   \b\b\b", stdout);
    fflush(stdout);
    chrom_end = 0;
    block_idx1 = 0;
    while (1) {
      next_unset_ul_unsafe_ck(marker_exclude1, &marker_uidx);
      ujj = g_epi_geno1_offsets[2 * block_idx1];
      marker_idx2 = 0;
      dptr = &(g_epi_all_chisq[block_idx1 * 2 * marker_ct2]);
      if (marker_uidx >= chrom_end) {
	chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx);
	chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
      }
      wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf));
      *wptr_start++ = ' ';
      wptr_start = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
      *wptr_start++ = ' ';
      marker_uidx2 = next_unset_ul_unsafe(marker_exclude2, marker_uidx_base);
      for (chrom_fo_idx2 = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx2); chrom_fo_idx2 < chrom_ct; chrom_fo_idx2++) {
	chrom_idx2 = chrom_info_ptr->chrom_file_order[chrom_fo_idx2];
	chrom_end2 = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx2 + 1];
	wptr_start2 = width_force(4, wptr_start, chrom_name_write(chrom_info_ptr, chrom_idx2, wptr_start));
	*wptr_start2++ = ' ';
	for (; marker_uidx2 < chrom_end2; ++marker_uidx2, next_unset_ul_ck(marker_exclude2, unfiltered_marker_ct, &marker_uidx2), ++marker_idx2, dptr = &(dptr[2])) {
	  if (marker_idx2 == ujj) {
	    marker_idx2 = g_epi_geno1_offsets[2 * block_idx1 + 1];
	    if (marker_idx2 == marker_ct2) {
	      goto epistasis_linear_regression_write_loop;
	    }
	    if (marker_idx2 > ujj) {
	      marker_uidx2 = jump_forward_unset_unsafe(marker_exclude2, marker_uidx2 + 1, marker_idx2 - ujj);
	      dptr = &(dptr[2 * (marker_idx2 - ujj)]);
	      if (marker_uidx2 >= chrom_end2) {
		break;
	      }
	    }
	  } else if (marker_idx2 == marker_ct2) {
	    goto epistasis_linear_regression_write_loop;
	  }
	  dxx = *dptr;
	  if (dxx != -1) {
	    wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr_start2);
	    *wptr++ = ' ';
	    // beta
	    wptr = width_force(12, wptr, dtoa_g(dptr[1], wptr));
            *wptr++ = ' ';
	    wptr = width_force(12, wptr, dtoa_g(dxx, wptr));
	    *wptr++ = ' ';
	    dxx = normdist(-sqrt(dxx)) * 2;
	    wptr = dtoa_g_wxp4x(MAXV(dxx, output_min_p), 12, ' ', wptr);
	    *wptr++ = '\n';
	    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	      goto epistasis_linear_regression_ret_WRITE_FAIL;
	    }
	    // could remove this writeback in --epi1 1 case
	    *dptr = -1;
	  }
	}
      }
    epistasis_linear_regression_write_loop:
      block_idx1++;
      marker_uidx++;
      if (block_idx1 >= idx1_block_size) {
        break;
      }
    }
    marker_idx1 += idx1_block_size;
    fputs("\b\b\b\b\b\b\b\b\b\b          \b\b\b\b\b\b\b\b\b\b", stdout);
    if (tests_complete >= pct_thresh) {
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      pct = (tests_complete * 100LLU) / tests_expected;
      if (pct < 100) {
        printf("\b\b%" PRIuPTR "%%", pct);
        fflush(stdout);
        pct_thresh = ((++pct) * ((uint64_t)tests_expected)) / 100;
      }
    }
  } while (marker_idx1 < marker_idx1_end);
  if (fclose_null(&outfile)) {
    goto epistasis_linear_regression_ret_WRITE_FAIL;
  }
  while (0) {
  epistasis_linear_regression_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  epistasis_linear_regression_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  epistasis_linear_regression_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  epistasis_linear_regression_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  epistasis_linear_regression_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  epistasis_linear_regression_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
  fclose_cond(outfile);
  // caller will free memory
  return retval;
}

int32_t epistasis_logistic_regression(pthread_t* threads, Epi_info* epi_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t marker_uidx_base, uintptr_t marker_ct1, uintptr_t* marker_exclude1, uintptr_t marker_idx1_start, uintptr_t marker_idx1_end, uintptr_t marker_ct2, uintptr_t* marker_exclude2, uint32_t is_triangular, uintptr_t job_size, uint64_t tests_expected, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uintptr_t* pheno_c, uint32_t parallel_idx, uint32_t parallel_tot, char* outname, char* outname_end, double output_min_p, uintptr_t* loadbuf_raw, uintptr_t* loadbuf, double* best_chisq, uint32_t* best_ids, uint32_t* n_sig_cts, uint32_t* fail_cts, uint32_t* gap_cts) {
  FILE* outfile = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t pheno_nm_cta4 = round_up_pow2(pheno_nm_ct, 4);
  uintptr_t pheno_nm_ctl2 = QUATERCT_TO_WORDCT(pheno_nm_ct);
  uintptr_t final_mask = get_final_mask(pheno_nm_ct);
  uintptr_t marker_uidx = marker_uidx_base;
  uintptr_t pct = 1;
  uintptr_t marker_uidx2 = 0;
  uintptr_t marker_idx1 = marker_idx1_start;
  uintptr_t marker_idx2 = 0;
  uint64_t pct_thresh = tests_expected / 100;
  uint64_t tests_complete = 0;
  uint32_t max_thread_ct = g_epi_thread_ct;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t chrom_end = 0;
  uint32_t chrom_idx = 0;
  uint32_t chrom_idx2 = 0;
  int32_t retval = 0;
  unsigned char* bigstack_mark2;
  uintptr_t* ulptr;
  char* wptr_start;
  char* wptr_start2;
  char* wptr;
  float* fptr;
  float* fptr2;
  double* dptr;
  uint32_t* uiptr;
  uint32_t* uiptr2;
  uint32_t* uiptr3;
  uint32_t* uiptr4;
  uint32_t* uiptr5;
  uintptr_t cur_bigstack_left;
  uintptr_t cur_workload;
  uintptr_t idx1_block_size;
  uintptr_t idx2_block_size;
  uintptr_t idx2_block_sizea16;
  uintptr_t marker_uidx_tmp;
  uintptr_t block_idx1;
  uintptr_t block_idx2;
  uintptr_t cur_idx2_block_size;
  uintptr_t chrom_end2;
  uintptr_t tidx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  double dxx;
  float fxx;
  uint32_t chrom_fo_idx;
  uint32_t chrom_fo_idx2;
  uint32_t is_last_block;
  uint32_t uii;
  uint32_t ujj;
  if (bigstack_alloc_ul(pheno_nm_ctl2, &g_epi_pheno_c)) {
    goto epistasis_logistic_regression_ret_NOMEM;
  }
  copy_bitarr_subset(pheno_c, pheno_nm, unfiltered_sample_ct, pheno_nm_ct, g_epi_pheno_c);
  g_epi_pheno_nm_ct = pheno_nm_ct;
  // per-thread buffers
  g_epi_logistic_mt = (Epi_logistic_multithread*)bigstack_alloc(max_thread_ct * sizeof(Epi_logistic_multithread));
  if (!g_epi_logistic_mt) {
    goto epistasis_logistic_regression_ret_NOMEM;
  }
  // param_ct_max = 4 (intercept, A, B, AB)
  for (tidx = 0; tidx < max_thread_ct; tidx++) {
    if (bigstack_alloc_f(pheno_nm_cta4 * 4, &(g_epi_logistic_mt[tidx].cur_covars_cov_major)) ||
        bigstack_alloc_f(4, &(g_epi_logistic_mt[tidx].coef)) ||
        bigstack_alloc_f(pheno_nm_cta4, &(g_epi_logistic_mt[tidx].pp)) ||
        bigstack_alloc_f(pheno_nm_ct, &(g_epi_logistic_mt[tidx].sample_1d_buf)) ||
        bigstack_alloc_f(pheno_nm_ct, &(g_epi_logistic_mt[tidx].pheno_buf)) ||
        bigstack_alloc_f(pheno_nm_ct * 4, &(g_epi_logistic_mt[tidx].param_1d_buf)) ||
        bigstack_alloc_f(pheno_nm_ct, &(g_epi_logistic_mt[tidx].param_1d_buf2)) ||
	bigstack_alloc_f(4 * 4, &(g_epi_logistic_mt[tidx].param_2d_buf)) ||
        bigstack_alloc_f(4 * 4, &(g_epi_logistic_mt[tidx].param_2d_buf2))) {
      goto epistasis_logistic_regression_ret_NOMEM;
    }
  }

  // claim up to half of memory with idx1 bufs; each marker currently costs:
  //   pheno_nm_ctl2 * sizeof(intptr_t) for geno buf
  //   4 * sizeof(int32_t) + sizeof(float) + marker_ct2 * 2 * sizeof(float)
  //     for other stuff (see epistasis_report() comment, starting from
  //     "offset"; main result buffer must be double-size to store both beta
  //     and chi-square stat)
  cur_bigstack_left = bigstack_left();
  ulii = 4 * CACHELINE - 3 * sizeof(int32_t) + max_thread_ct * (5 * (CACHELINE - 4));
  if (cur_bigstack_left >= ulii) {
    cur_bigstack_left -= ulii;
  }
  ulii = pheno_nm_ctl2 * sizeof(intptr_t) + 4 * sizeof(int32_t) + sizeof(float) + marker_ct2 * 2 * sizeof(float);
  idx1_block_size = cur_bigstack_left / (ulii * 2 + 1);
  if (!idx1_block_size) {
    goto epistasis_logistic_regression_ret_NOMEM;
  }
  if (idx1_block_size > job_size) {
    idx1_block_size = job_size;
  }
  // pad to avoid threads writing to same cacheline
  ulii = (max_thread_ct - 1) * 15 + idx1_block_size;
  bigstack_alloc_ui(idx1_block_size * 2, &g_epi_geno1_offsets);
  bigstack_alloc_ul(pheno_nm_ctl2 * idx1_block_size, &g_epi_geno1);
  bigstack_alloc_f(idx1_block_size * marker_ct2 * 2, &g_epi_all_chisq_f);
  bigstack_alloc_f(ulii, &g_epi_best_chisq_f1);
  bigstack_alloc_ui(ulii, &g_epi_best_id1);
  bigstack_alloc_ui(ulii, &g_epi_n_sig_ct1);
  bigstack_alloc_ui(ulii, &g_epi_fail_ct1);
  for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
    g_epi_geno1[block_idx1 * pheno_nm_ctl2 + pheno_nm_ctl2 - 1] = 0;
  }
  if (is_triangular) {
    fill_uint_zero(2 * idx1_block_size, g_epi_geno1_offsets);
  }

  ulii = pheno_nm_ctl2 * sizeof(intptr_t) + max_thread_ct * (3 * sizeof(int32_t) + sizeof(double));
  idx2_block_size = (bigstack_left() - (CACHELINE - sizeof(intptr_t)) - max_thread_ct * (3 * (CACHELINE - sizeof(int32_t)) + (CACHELINE - sizeof(float)))) / ulii;
  if (idx2_block_size > marker_ct2) {
    idx2_block_size = marker_ct2;
  }
  idx2_block_size = round_up_pow2(idx2_block_size, 16);

  memcpy(outname_end, ".epi.cc", 8);
  if (parallel_tot > 1) {
    outname_end[7] = '.';
    uint32toa_x(parallel_idx + 1, '\0', &(outname_end[8]));
  }
  if (fopen_checked(outname, "w", &outfile)) {
    goto epistasis_logistic_regression_ret_OPEN_FAIL;
  }
  if (!parallel_idx) {
    wptr = memcpya(g_textbuf, "CHR1 ", 5);
    wptr = fw_strcpyn(plink_maxsnp, 4, "SNP1", wptr);
    wptr = memcpya(wptr, " CHR2 ", 6);
    wptr = fw_strcpyn(plink_maxsnp, 4, "SNP2", wptr);
    wptr = memcpya(wptr, "       OR_INT         STAT            P \n", 41);
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
      goto epistasis_logistic_regression_ret_WRITE_FAIL;
    }
  }

  bigstack_mark2 = g_bigstack_base;
  while (1) {
    if (!idx2_block_size) {
      goto epistasis_logistic_regression_ret_NOMEM;
    }
    if (!(bigstack_alloc_ul(pheno_nm_ctl2 * idx2_block_size, &g_epi_geno2) ||
          bigstack_alloc_f(max_thread_ct * idx2_block_size, &g_epi_best_chisq_f2) ||
          bigstack_alloc_ui(max_thread_ct * idx2_block_size, &g_epi_best_id2) ||
          bigstack_alloc_ui(max_thread_ct * idx2_block_size, &g_epi_n_sig_ct2) ||
          bigstack_alloc_ui(max_thread_ct * idx2_block_size, &g_epi_fail_ct2))) {
      break;
    }
    bigstack_reset(bigstack_mark2);
    idx2_block_size -= 16;
  }
  for (block_idx2 = 0; block_idx2 < idx2_block_size; block_idx2++) {
    g_epi_geno2[block_idx2 * pheno_nm_ctl2 + pheno_nm_ctl2 - 1] = 0;
  }
  marker_uidx = next_unset_ul_unsafe(marker_exclude1, marker_uidx_base);
  if (marker_idx1) {
    marker_uidx = jump_forward_unset_unsafe(marker_exclude1, marker_uidx + 1, marker_idx1);
  }
  wptr = memcpya(g_logbuf, "C/C --epistasis to ", 19);
  wptr = strcpya(wptr, outname);
  memcpy(wptr, " ... ", 6);
  wordwrapb(16); // strlen("99% [processing]")
  logprintb();
  fputs("0%", stdout);
  do {
    fputs(" [processing]", stdout);
    fflush(stdout);
    if (idx1_block_size > marker_idx1_end - marker_idx1) {
      idx1_block_size = marker_idx1_end - marker_idx1;
      if (idx1_block_size < max_thread_ct) {
        max_thread_ct = idx1_block_size;
        g_epi_thread_ct = max_thread_ct;
      }
    }
    g_epi_marker_idx1 = marker_idx1;
    fptr = g_epi_all_chisq_f;
    fptr2 = &(g_epi_all_chisq_f[idx1_block_size * marker_ct2 * 2]);
    do {
      *fptr = -1;
      fptr = &(fptr[2]);
    } while (fptr < fptr2);
    marker_uidx_tmp = marker_uidx;
    if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto epistasis_logistic_regression_ret_READ_FAIL;
    }
    cur_workload = idx1_block_size * marker_ct2;
    if (is_triangular) {
      for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
        ulii = block_idx1 + marker_idx1 + 1;
        cur_workload -= ulii;
        g_epi_geno1_offsets[2 * block_idx1 + 1] = ulii;
      }
    } else {
      fill_uint_zero(2 * idx1_block_size, g_epi_geno1_offsets);
      marker_uidx2 = marker_uidx_base;
      marker_idx2 = 0;
    }
    tests_complete += cur_workload;
    ulii = 0; // total number of tests
    g_epi_idx1_block_bounds[0] = 0;
    g_epi_idx1_block_bounds16[0] = 0;
    block_idx1 = 0;
    for (tidx = 1; tidx < max_thread_ct; tidx++) {
      uljj = (((uint64_t)cur_workload) * tidx) / max_thread_ct;
      if (is_triangular) {
        do {
          ulii += marker_ct2 - g_epi_geno1_offsets[2 * block_idx1 + 1];
          block_idx1++;
	} while (ulii < uljj);
      } else {
        do {
	  ulii += marker_ct2;
	  block_idx1++;
	} while (ulii < uljj);
      }
      uii = block_idx1 - g_epi_idx1_block_bounds[tidx - 1];
      g_epi_idx1_block_bounds[tidx] = block_idx1;
      g_epi_idx1_block_bounds16[tidx] = g_epi_idx1_block_bounds16[tidx - 1] + round_up_pow2_ui(uii, 16);
    }
    g_epi_idx1_block_bounds[max_thread_ct] = idx1_block_size;
    chrom_end = 0;
    for (block_idx1 = 0; block_idx1 < idx1_block_size; marker_uidx_tmp++, block_idx1++) {
      if (IS_SET(marker_exclude1, marker_uidx_tmp)) {
	marker_uidx_tmp = next_unset_ul_unsafe(marker_exclude1, marker_uidx_tmp);
        if (fseeko(bedfile, bed_offset + (marker_uidx_tmp * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
          goto epistasis_logistic_regression_ret_READ_FAIL;
	}
      }
      // marker_reverse deliberately flipped
      if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, !IS_SET(marker_reverse, marker_uidx_tmp), bedfile, loadbuf_raw, loadbuf)) {
        goto epistasis_logistic_regression_ret_READ_FAIL;
      }
      // rotate to hom A1 = 10, het = 01, hom A2 = 00, missing = 11, to allow
      // inner loop to use ordinary multiplication
      // this is a bit redundant with the forced reverse, but it's not a
      // bottleneck
      rotate_plink1_to_a2ct_and_copy(loadbuf, &(g_epi_geno1[block_idx1 * pheno_nm_ctl2]), pheno_nm_ctl2);
      if (!is_triangular) {
	if (!IS_SET(marker_exclude2, marker_uidx_tmp)) {
          // do not compare against self
          marker_idx2 += marker_uidx_tmp - marker_uidx2 - popcount_bit_idx(marker_exclude2, marker_uidx2, marker_uidx_tmp);
          marker_uidx2 = marker_uidx_tmp;
          g_epi_geno1_offsets[2 * block_idx1] = marker_idx2;
          g_epi_geno1_offsets[2 * block_idx1 + 1] = marker_idx2 + 1;
          gap_cts[block_idx1 + marker_idx1] = 1;
	}
      }
    }
    marker_uidx2 = next_unset_ul_unsafe(marker_exclude2, marker_uidx_base);
    if (is_triangular) {
      marker_idx2 = marker_idx1 + 1;
      marker_uidx2 = jump_forward_unset_unsafe(marker_exclude2, marker_uidx2 + 1, marker_idx2);
    } else {
      marker_idx2 = 0;
    }
    if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto epistasis_logistic_regression_ret_READ_FAIL;
    }
    cur_idx2_block_size = idx2_block_size;
    do {
      if (cur_idx2_block_size > marker_ct2 - marker_idx2) {
        cur_idx2_block_size = marker_ct2 - marker_idx2;
      }
      for (block_idx2 = 0; block_idx2 < cur_idx2_block_size; marker_uidx2++, block_idx2++) {
        if (IS_SET(marker_exclude2, marker_uidx2)) {
          marker_uidx2 = next_unset_ul_unsafe(marker_exclude2, marker_uidx2);
          if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
            goto epistasis_logistic_regression_ret_READ_FAIL;
	  }
	}
        ulptr = &(g_epi_geno2[block_idx2 * pheno_nm_ctl2]);
	// marker_reverse deliberately flipped
	if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, !IS_SET(marker_reverse, marker_uidx2), bedfile, loadbuf_raw, loadbuf)) {
	  goto epistasis_logistic_regression_ret_READ_FAIL;
	}
	rotate_plink1_to_a2ct_and_copy(loadbuf, ulptr, pheno_nm_ctl2);
      }
      g_epi_idx2_block_size = cur_idx2_block_size;
      g_epi_idx2_block_start = marker_idx2;
      idx2_block_sizea16 = round_up_pow2(cur_idx2_block_size, 16);
      fill_uint_zero(idx1_block_size + 15 * (max_thread_ct - 1), g_epi_n_sig_ct1);
      fill_uint_zero(idx1_block_size + 15 * (max_thread_ct - 1), g_epi_fail_ct1);
      fill_uint_zero(idx2_block_sizea16 * max_thread_ct, g_epi_n_sig_ct2);
      fill_uint_zero(idx2_block_sizea16 * max_thread_ct, g_epi_fail_ct2);
      for (tidx = 0; tidx < max_thread_ct; tidx++) {
        ulii = g_epi_idx1_block_bounds[tidx];
        uljj = g_epi_idx1_block_bounds[tidx + 1] - ulii;
        fptr = &(g_epi_best_chisq_f1[g_epi_idx1_block_bounds[tidx]]);
	fptr2 = &(g_epi_all_chisq_f[(marker_idx1 + ulii) * 2]);
	for (ulkk = 0; ulkk < uljj; ulkk++) {
          *fptr++ = fptr2[ulkk * 2];
	}
        ulii = g_epi_geno1_offsets[2 * ulii + 1];
        if (ulii < marker_idx2 + cur_idx2_block_size) {
          if (ulii <= marker_idx2) {
            ulii = 0;
	  } else {
            ulii -= marker_idx2;
	  }
          uljj = cur_idx2_block_size - ulii;
	  fptr = &(g_epi_best_chisq_f2[tidx * idx2_block_sizea16 + ulii]);
	  fptr2 = &(g_epi_all_chisq_f[(marker_idx2 + ulii) * 2]);
          for (ulkk = 0; ulkk < uljj; ulkk++) {
            *fptr++ = fptr2[ulkk * 2];
	  }
	}
      }
      is_last_block = (marker_idx2 + cur_idx2_block_size >= marker_ct2);
      if (spawn_threads2(threads, &epi_logistic_thread, max_thread_ct, is_last_block)) {
	goto epistasis_logistic_regression_ret_THREAD_CREATE_FAIL;
      }
      epi_logistic_thread((void*)0);
      join_threads2(threads, max_thread_ct, is_last_block);
      // merge best_chisq, best_ids, fail_cts
      for (tidx = 0; tidx < max_thread_ct; tidx++) {
	ulii = g_epi_idx1_block_bounds[tidx];
	uljj = g_epi_idx1_block_bounds[tidx + 1] - ulii;
	uii = g_epi_idx1_block_bounds16[tidx];
	fptr = &(g_epi_best_chisq_f1[uii]);
	uiptr = &(g_epi_best_id1[uii]);
	uiptr2 = &(g_epi_n_sig_ct1[uii]);
	uiptr3 = &(g_epi_fail_ct1[uii]);
	ulii += marker_idx1;
	dptr = &(best_chisq[ulii]);
	uiptr4 = &(n_sig_cts[ulii]);
	uiptr5 = &(fail_cts[ulii]);
	for (block_idx1 = 0; block_idx1 < uljj; block_idx1++, dptr++, uiptr4++, uiptr5++) {
	  dxx = (double)(*fptr++);
	  if (dxx > (*dptr)) {
	    *dptr = dxx;
	    best_ids[block_idx1 + ulii] = uiptr[block_idx1];
	  }
	  *uiptr4 += *uiptr2++;
	  *uiptr5 += *uiptr3++;
	}
      }
      if (is_triangular) {
	for (tidx = 0; tidx < max_thread_ct; tidx++) {
	  block_idx2 = g_epi_geno1_offsets[2 * g_epi_idx1_block_bounds[tidx] + 1];
	  if (block_idx2 <= marker_idx2) {
	    block_idx2 = 0;
	  } else {
	    block_idx2 -= marker_idx2;
	  }
	  fptr = &(g_epi_best_chisq_f2[tidx * idx2_block_sizea16 + block_idx2]);
	  uiptr = &(g_epi_best_id2[tidx * idx2_block_sizea16]);
	  uiptr2 = &(g_epi_n_sig_ct2[tidx * idx2_block_sizea16 + block_idx2]);
	  uiptr3 = &(g_epi_fail_ct2[tidx * idx2_block_sizea16 + block_idx2]);
	  dptr = &(best_chisq[block_idx2 + marker_idx2]);
	  uiptr4 = &(n_sig_cts[block_idx2 + marker_idx2]);
	  uiptr5 = &(fail_cts[block_idx2 + marker_idx2]);
	  for (; block_idx2 < cur_idx2_block_size; block_idx2++, dptr++, uiptr4++, uiptr5++) {
	    dxx = (double)(*fptr++);
	    if (dxx > (*dptr)) {
	      *dptr = dxx;
	      best_ids[block_idx2 + marker_idx2] = uiptr[block_idx2];
	    }
	    *uiptr4 += *uiptr2++;
	    *uiptr5 += *uiptr3++;
	  }
	}
      }
      marker_idx2 += cur_idx2_block_size;
    } while (marker_idx2 < marker_ct2);
    fputs("\b\b\b\b\b\b\b\b\b\b\bwriting]   \b\b\b", stdout);
    fflush(stdout);
    chrom_end = 0;
    block_idx1 = 0;
    while (1) {
      next_unset_ul_unsafe_ck(marker_exclude1, &marker_uidx);
      ujj = g_epi_geno1_offsets[2 * block_idx1];
      marker_idx2 = 0;
      fptr = &(g_epi_all_chisq_f[block_idx1 * 2 * marker_ct2]);
      if (marker_uidx >= chrom_end) {
	chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx);
	chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
      }
      wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf));
      *wptr_start++ = ' ';
      wptr_start = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
      *wptr_start++ = ' ';
      marker_uidx2 = next_unset_ul_unsafe(marker_exclude2, marker_uidx_base);
      for (chrom_fo_idx2 = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx2); chrom_fo_idx2 < chrom_ct; chrom_fo_idx2++) {
	chrom_idx2 = chrom_info_ptr->chrom_file_order[chrom_fo_idx2];
	chrom_end2 = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx2 + 1];
	wptr_start2 = width_force(4, wptr_start, chrom_name_write(chrom_info_ptr, chrom_idx2, wptr_start));
	*wptr_start2++ = ' ';
	for (; marker_uidx2 < chrom_end2; ++marker_uidx2, next_unset_ul_ck(marker_exclude2, unfiltered_marker_ct, &marker_uidx2), ++marker_idx2, fptr = &(fptr[2])) {
	  if (marker_idx2 == ujj) {
	    marker_idx2 = g_epi_geno1_offsets[2 * block_idx1 + 1];
	    if (marker_idx2 == marker_ct2) {
	      goto epistasis_logistic_regression_write_loop;
	    }
	    if (marker_idx2 > ujj) {
	      marker_uidx2 = jump_forward_unset_unsafe(marker_exclude2, marker_uidx2 + 1, marker_idx2 - ujj);
	      fptr = &(fptr[2 * (marker_idx2 - ujj)]);
	      if (marker_uidx2 >= chrom_end2) {
		break;
	      }
	    }
	  } else if (marker_idx2 == marker_ct2) {
	    goto epistasis_logistic_regression_write_loop;
	  }
	  fxx = *fptr;
	  if (fxx != -1) {
	    dxx = (double)fxx;
	    wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr_start2);
	    *wptr++ = ' ';
	    // odds ratio
	    wptr = width_force(12, wptr, dtoa_g(exp((double)fptr[1]), wptr));
            *wptr++ = ' ';
	    wptr = width_force(12, wptr, ftoa_g(fxx, wptr));
	    *wptr++ = ' ';
	    dxx = normdist(-sqrt(dxx)) * 2;
	    wptr = dtoa_g_wxp4x(MAXV(dxx, output_min_p), 12, ' ', wptr);
	    *wptr++ = '\n';
	    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	      goto epistasis_logistic_regression_ret_WRITE_FAIL;
	    }
	    // could remove this writeback in --epi1 1 case
	    *fptr = -1;
	  }
	}
      }
    epistasis_logistic_regression_write_loop:
      block_idx1++;
      marker_uidx++;
      if (block_idx1 >= idx1_block_size) {
        break;
      }
    }
    marker_idx1 += idx1_block_size;
    fputs("\b\b\b\b\b\b\b\b\b\b          \b\b\b\b\b\b\b\b\b\b", stdout);
    if (tests_complete >= pct_thresh) {
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      pct = (tests_complete * 100LLU) / tests_expected;
      if (pct < 100) {
        printf("\b\b%" PRIuPTR "%%", pct);
        fflush(stdout);
        pct_thresh = ((++pct) * ((uint64_t)tests_expected)) / 100;
      }
    }
  } while (marker_idx1 < marker_idx1_end);
  if (fclose_null(&outfile)) {
    goto epistasis_logistic_regression_ret_WRITE_FAIL;
  }
  while (0) {
  epistasis_logistic_regression_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  epistasis_logistic_regression_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  epistasis_logistic_regression_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  epistasis_logistic_regression_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  epistasis_logistic_regression_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
  fclose_cond(outfile);
  // caller will free memory
  return retval;
}

int32_t epistasis_report(pthread_t* threads, Epi_info* epi_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct2, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uint32_t ctrl_ct, uintptr_t* pheno_c, double* pheno_d, uint32_t parallel_idx, uint32_t parallel_tot, char* outname, char* outname_end, double output_min_p, double glm_vif_thresh, Set_info* sip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct);
  uintptr_t final_mask = get_final_mask(pheno_nm_ct);
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t marker_uidx_base = next_unset_unsafe(marker_exclude, 0);
  uintptr_t marker_uidx = marker_uidx_base;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t modifier = epi_ip->modifier;
  uint32_t is_fast = modifier & EPI_FAST;
  uint32_t is_boost = (modifier / EPI_FAST_BOOST) & 1;
  uint32_t do_joint_effects = modifier & EPI_FAST_JOINT_EFFECTS;
  uint32_t no_ueki = modifier & EPI_FAST_NO_UEKI;
  uint32_t is_case_only = (modifier / EPI_FAST_CASE_ONLY) & 1;
  uint32_t is_triangular = 1;
  uint32_t is_custom_set1 = modifier & (EPI_SET_BY_SET | EPI_SET_BY_ALL)? 1 : 0;
  uint32_t is_set_by_set = modifier & EPI_SET_BY_SET;
  uint32_t tot_stride = 6 - 3 * is_case_only;
  uint32_t no_p_value = modifier & EPI_FAST_NO_P_VALUE;
  uint32_t case_only_gap = epi_ip->case_only_gap;
  uint32_t is_case_only_window = (is_case_only && case_only_gap);
  uint32_t case_ct = pheno_nm_ct - ctrl_ct;
  uint32_t cellminx3 = 0;
  uintptr_t case_ctl2 = QUATERCT_TO_WORDCT(case_ct);
  uintptr_t case_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(case_ct);
  uintptr_t ctrl_ctl2 = QUATERCT_TO_WORDCT(ctrl_ct);
  uintptr_t case_ctv3 = BITCT_TO_ALIGNED_WORDCT(case_ct);
  uintptr_t ctrl_ctv3 = BITCT_TO_ALIGNED_WORDCT(ctrl_ct);
  uintptr_t case_ctsplit = 3 * case_ctv3;
  uintptr_t ctrl_ctsplit = 3 * ctrl_ctv3;
  uintptr_t pct = 1;
  uintptr_t marker_uidx2 = 0;
  uintptr_t marker_uidx2_trail = 0;
  uintptr_t marker_idx2 = 0;
  uintptr_t marker_idx2_trail = 0;
  uint64_t tests_thrown_out = 0;
  uint64_t tests_complete = 0;
  uint32_t max_thread_ct = g_thread_ct;
  uint32_t chrom_idx = 0;
  uint32_t chrom_end = 0;
  uint32_t last_pos = 0;
  uint32_t first_pos = 0;
  uint32_t uii = 0;
  int32_t retval = 0;
  uint32_t* gap_cts = nullptr;
  uintptr_t* ctrlbuf = nullptr;
  uintptr_t* marker_exclude1 = nullptr;
  uintptr_t* ulptr = nullptr;
  uintptr_t ularr[sizeof(double) / BYTECT];
  uintptr_t* casebuf;
  uintptr_t* loadbuf;
  uintptr_t* marker_exclude2;
  double* best_chisq;
  uint32_t* best_ids;
  uint32_t* n_sig_cts;
  uint32_t* fail_cts;
  uint32_t* marker_idx_to_uidx;
  unsigned char* bigstack_mark2;
  unsigned char* bigstack_mark3;
  char* wptr_start;
  char* wptr_start2;
  char* wptr;
  double* dptr;
  double* dptr2;
  uint32_t* uiptr;
  uint32_t* uiptr2;
  uint32_t* uiptr3;
  uint32_t* uiptr4;
  uint32_t* uiptr5;
  uint64_t tests_expected;
  uint64_t pct_thresh;
  double dxx;
  uintptr_t marker_ct1;
  uintptr_t tot_ctsplit;
  uintptr_t job_size;
  uintptr_t cur_bigstack_left;
  uintptr_t cur_workload;
  uintptr_t marker_idx1_start;
  uintptr_t marker_idx1;
  uintptr_t marker_idx1_end;
  uintptr_t idx1_block_size;
  uintptr_t idx2_block_size;
  uintptr_t idx2_block_sizea16;
  uintptr_t marker_uidx_tmp;
  uintptr_t block_idx1;
  uintptr_t block_idx2;
  uintptr_t cur_idx2_block_size;
  uintptr_t tidx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t chrom_end2;
  uint32_t chrom_fo_idx;
  uint32_t chrom_fo_idx2;
  uint32_t chrom_idx2;
  uint32_t cur_window_end;
  uint32_t is_last_block;
  uint32_t missing_ct;
  uint32_t ujj;

  // common initialization between --epistasis and --fast-epistasis: remove
  // monomorphic and non-autosomal diploid sites
  if (is_custom_set1) {
    if (!sip->ct) {
      sprintf(g_logbuf, "Error: --%sepistasis set-by-%s requires a variant set to be loaded.\n", is_fast? "fast-" : "", is_set_by_set? "set" : "all");
      goto epistasis_report_ret_INVALID_CMDLINE_2;
    } else if (!is_set_by_set) {
      if (sip->ct > 1) {
	logerrprint("Error: --{fast-}epistasis set-by-all requires exactly one set.  (--set-names or\n--set-collapse-all may be handy here.\n");
	goto epistasis_report_ret_INVALID_CMDLINE;
      }
    } else if (sip->ct > 2) {
      logerrprint("Error: --{fast-}epistasis set-by-set requires exactly one or two sets.\n(--set-names or --set-collapse-all may be handy here.)\n");
      goto epistasis_report_ret_INVALID_CMDLINE;
    }
    if (bigstack_alloc_ul(unfiltered_marker_ctl, &marker_exclude1)) {
      goto epistasis_report_ret_NOMEM;
    }
    unpack_set_unfiltered(marker_ct2, unfiltered_marker_ct, marker_exclude, sip->setdefs[0], marker_exclude1);
    if (is_set_by_set && (sip->ct == 1)) {
      marker_ct2 = unfiltered_marker_ct - popcount_longs(marker_exclude1, unfiltered_marker_ctl);
    } else {
      is_triangular = 0;
    }
    // if set-by-set with two sets, wait till after monomorphic sites are
    // removed to unpack 2nd set
  }
  if (pheno_nm_ct >= 0x20000000) {
    // may as well document the existence of sub-2b overflow conditions even
    // though they'll never come up
    logerrprint("Error: --{fast-}epistasis does not support >= 2^29 samples.\n");
    goto epistasis_report_ret_INVALID_CMDLINE;
  }
  if (!pheno_d) {
    if ((case_ct < 2) || ((!is_case_only) && (ctrl_ct < 2))) {
      sprintf(g_logbuf, "Error: --%sepistasis requires at least two cases%s.\n", is_fast? "fast-" : "", is_case_only? "" : " and two controls");
      goto epistasis_report_ret_INVALID_CMDLINE_2;
    }
    if (bigstack_alloc_ul(case_ctv2 + ctrl_ctl2, &casebuf)) {
      goto epistasis_report_ret_NOMEM;
    }
    ctrlbuf = &(casebuf[case_ctv2]);
    ctrlbuf[ctrl_ctl2 - 1] = 0;
  } else {
    case_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
    if (bigstack_alloc_ul(case_ctv2, &casebuf)) {
      goto epistasis_report_ret_NOMEM;
    }
  }
  casebuf[case_ctv2 - 2] = 0;
  casebuf[case_ctv2 - 1] = 0;
  // marker_exclude2 should be on top since we might free it
  if (bigstack_alloc_ul(unfiltered_sample_ctv2, &loadbuf) ||
      bigstack_alloc_ul(unfiltered_marker_ctl, &marker_exclude2)) {
    goto epistasis_report_ret_NOMEM;
  }
  loadbuf[unfiltered_sample_ctv2 - 2] = 0;
  loadbuf[unfiltered_sample_ctv2 - 1] = 0;
  if ((!is_set_by_set) || (sip->ct == 2)) {
    memcpy(marker_exclude2, marker_exclude, unfiltered_marker_ctl * sizeof(intptr_t));
  } else {
    memcpy(marker_exclude2, marker_exclude1, unfiltered_marker_ctl * sizeof(intptr_t));
  }
  if (do_joint_effects && epi_ip->je_cellmin) {
    cellminx3 = epi_ip->je_cellmin * 3;
    if ((case_ct < cellminx3 * 3) || ((!is_case_only) && (ctrl_ct < cellminx3 * 3))) {
      sprintf(g_logbuf, "Error: Too few cases or controls for --je-cellmin %u.\n", epi_ip->je_cellmin);
      goto epistasis_report_ret_INVALID_CMDLINE_2;
    }
    ulii = case_ctl2;
    if ((!is_case_only) && (ctrl_ctl2 > case_ctl2)) {
      ulii = ctrl_ctl2;
    }
    if (bigstack_alloc_ul(ulii, &ulptr)) {
      goto epistasis_report_ret_NOMEM;
    }
    fill_quatervec_55(ulii * BITCT2, ulptr);
  }
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    if (is_set(chrom_info_ptr->haploid_mask, chrom_idx)) {
      uii = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
      fill_bits(uii, chrom_end - uii, marker_exclude2);
      marker_uidx = chrom_end;
      continue;
    }
    // may want to keep two window sizes' raw data loaded for marker 1, to
    // halve the number of non-sequential seeks?
    if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto epistasis_report_ret_READ_FAIL;
    }
    while (marker_uidx < chrom_end) {
      if (is_set(marker_exclude2, marker_uidx)) {
	marker_uidx = next_unset(marker_exclude2, marker_uidx, chrom_end);
	if (marker_uidx == chrom_end) {
	  break;
	}
	if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	  goto epistasis_report_ret_READ_FAIL;
	}
      }
      if ((!no_ueki) && (!cellminx3)) {
	if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, 0, bedfile, loadbuf, casebuf)) {
	  goto epistasis_report_ret_READ_FAIL;
	}
	if (is_boost) {
	  if (less_than_two_genotypes(casebuf, pheno_nm_ct)) {
	    SET_BIT(marker_uidx, marker_exclude2);
	  }
	} else {
	  if (is_monomorphic(casebuf, pheno_nm_ct)) {
	    SET_BIT(marker_uidx, marker_exclude2);
	  }
	}
      } else {
        if (load_and_split(unfiltered_sample_ct, pheno_nm, pheno_c, bedfile, loadbuf, casebuf, ctrlbuf)) {
          goto epistasis_report_ret_READ_FAIL;
	}
	if (no_ueki) {
	  if (is_monomorphic(casebuf, case_ct) || ((!is_case_only) && is_monomorphic(ctrlbuf, ctrl_ct))) {
	    SET_BIT(marker_uidx, marker_exclude2);
	  }
	} else {
	  genovec_3freq(casebuf, ulptr, case_ctl2, &missing_ct, &uii, &ujj);
	  if ((uii < cellminx3) || (ujj < cellminx3) || (case_ct - uii - ujj - missing_ct < cellminx3)) {
	    SET_BIT(marker_uidx, marker_exclude2);
	  } else if (!is_case_only) {
	    genovec_3freq(ctrlbuf, ulptr, ctrl_ctl2, &missing_ct, &uii, &ujj);
	    if ((uii < cellminx3) || (ujj < cellminx3) || (ctrl_ct - uii - ujj - missing_ct < cellminx3)) {
	      SET_BIT(marker_uidx, marker_exclude2);
	    }
	  }
	}
      }
      marker_uidx++;
    }
  }
  ulii = unfiltered_marker_ct - popcount_longs(marker_exclude2, unfiltered_marker_ctl);
  if ((!ulii) || ((ulii == 1) && is_triangular)) {
    goto epistasis_report_ret_TOO_FEW_MARKERS;
  }
  if (ulii != marker_ct2) {
    if (!cellminx3) {
      LOGPRINTF("--%sepistasis: Skipping %" PRIuPTR " monomorphic/non-autosomal site%s.\n", is_fast? "fast-" : "", marker_ct2 - ulii, (marker_ct2 - ulii == 1)? "" : "s");
    } else {
      LOGPRINTF("--%sepistasis: Skipping %" PRIuPTR " site%s due to --je-cellmin setting.\n", is_fast? "fast-" : "", marker_ct2 - ulii, (marker_ct2 - ulii == 1)? "" : "s");
      bigstack_reset(ulptr);
    }
    marker_uidx_base = next_unset_ul_unsafe(marker_exclude2, marker_uidx_base);
  } else if ((!is_custom_set1) || (!is_set_by_set)) {
    bigstack_reset(marker_exclude2);
    marker_exclude2 = marker_exclude;
  }
  if (is_triangular) {
    if (!marker_exclude1) {
      marker_exclude1 = marker_exclude2;
    }
    marker_ct1 = ulii;
    marker_ct2 = ulii;
    tests_expected = ((((uint64_t)marker_ct1) * (marker_ct1 - 1)) / 2);
  } else {
    bitvec_or(marker_exclude2, unfiltered_marker_ctl, marker_exclude1);
    marker_ct1 = unfiltered_marker_ct - popcount_longs(marker_exclude1, unfiltered_marker_ctl);
    if (sip->ct == 2) {
      if (bigstack_alloc_ul(unfiltered_marker_ctl, &ulptr)) {
	goto epistasis_report_ret_NOMEM;
      }
      memcpy(ulptr, marker_exclude2, unfiltered_marker_ctl * sizeof(intptr_t));
      unpack_set_unfiltered(marker_ct2, unfiltered_marker_ct, marker_exclude, sip->setdefs[1], marker_exclude2);
      bitvec_or(ulptr, unfiltered_marker_ctl, marker_exclude2);
      bigstack_reset(ulptr);
      marker_ct2 = unfiltered_marker_ct - popcount_longs(marker_exclude2, unfiltered_marker_ctl);
    } else {
      marker_ct2 = ulii;
    }
    tests_expected = ((uint64_t)marker_ct1) * marker_ct2;
    if (!tests_expected) {
      goto epistasis_report_ret_TOO_FEW_MARKERS;
    }
  }
  if (parallel_tot > 1) {
    if (marker_ct1 < (1 + is_triangular) * parallel_tot) {
      sprintf(g_logbuf, "Error: Too few loci remaining for --parallel %u %u + --%sepistasis.\n", parallel_idx + 1, parallel_tot, is_fast? "fast-" : "");
      goto epistasis_report_ret_INVALID_CMDLINE_2;
    }
    if (is_triangular) {
      // If there are n markers, and we're computing the usual upper right
      // triangle, first row has n-1 entries, second row has n-2, etc.
      // Total entry count is n(n-1)/2; total entry count starting from row r
      // is (n-r)(n-r-1)/2... upside-down triangle_divide() calls produce a
      // good partition.
      // Divide first to avoid 64-bit integer overflow (!) on really huge jobs.
      // (Multiply-by-2 is there because triangle_divide() takes n(n-1) instead
      // of n(n-1)/2 as first parameter.)
      pct_thresh = (2 * tests_expected) / parallel_tot;
      // If parallel_idx == 0, the marker_ct >= 2 * parallel_tot condition
      // ensures the precision loss from dividing and remultiplying does not
      // cause the first marker to be dropped.
      marker_idx1_start = triangle_divide(pct_thresh * (parallel_tot - parallel_idx), -1);
      marker_idx1_end = triangle_divide(pct_thresh * (parallel_tot - parallel_idx - 1), -1);
      tests_expected = ((((uint64_t)marker_idx1_start) * (marker_idx1_start - 1)) - (((uint64_t)marker_idx1_end) * (marker_idx1_end - 1))) / 2;
      marker_idx1_start = marker_ct1 - marker_idx1_start;
      marker_idx1_end = marker_ct1 - marker_idx1_end;
    } else {
      marker_idx1_start = (parallel_idx * ((uint64_t)marker_ct1)) / parallel_tot;
      marker_idx1_end = ((parallel_idx + 1) * ((uint64_t)marker_ct1)) / parallel_tot;
      tests_expected = (marker_idx1_end - marker_idx1_start) * ((uint64_t)marker_ct2);
    }
  } else {
    marker_idx1_start = 0;
    marker_idx1_end = marker_ct1;
  }
  marker_idx1 = marker_idx1_start;
  job_size = marker_idx1_end - marker_idx1_start;
  if (max_thread_ct > job_size) {
    max_thread_ct = job_size;
  }
  if (bigstack_calloc_d(marker_ct1, &best_chisq) ||
      bigstack_calloc_ui(marker_ct1, &best_ids) ||
      bigstack_calloc_ui(marker_ct1, &n_sig_cts) ||
      bigstack_calloc_ui(marker_ct1, &fail_cts) ||
      bigstack_alloc_ui(max_thread_ct + 1, &g_epi_idx1_block_bounds) ||
      bigstack_alloc_ui(max_thread_ct, &g_epi_idx1_block_bounds16)) {
    goto epistasis_report_ret_NOMEM;
  }
  if (is_case_only_window || (!is_triangular)) {
    if (bigstack_calloc_ui(marker_ct1, &gap_cts)) {
      goto epistasis_report_ret_NOMEM;
    }
  }
  bigstack_mark3 = g_bigstack_base;

  g_epi_thread_ct = max_thread_ct;
  g_epi_case_ct = case_ct;
  g_epi_flag = modifier;
  g_epi_marker_ct = marker_ct2;
  g_epi_cellmin = cellminx3 / 3;
  // might want to provide a Bonferroni correction interface...
  if (is_boost) {
    if (epi_ip->epi1 == 0.0) {
      dxx = 0.000005;
    } else {
      dxx = epi_ip->epi1;
    }
    g_epi_alpha1sq[0] = inverse_chiprob(dxx, 4);
    g_epi_alpha1sq[1] = inverse_chiprob(dxx, 2);
    g_epi_alpha1sq[2] = inverse_chiprob(dxx, 1);
    g_epi_alpha2sq[0] = inverse_chiprob(epi_ip->epi2, 4);
    if (g_epi_alpha1sq[0] == g_epi_alpha2sq[0]) {
      // count final instead of screening p-value hits
      g_epi_alpha2sq[0] *= 1 + SMALL_EPSILON;
      g_epi_alpha2sq[1] = g_epi_alpha1sq[1] * (1 + SMALL_EPSILON);
      g_epi_alpha2sq[2] = g_epi_alpha1sq[2] * (1 + SMALL_EPSILON);
    } else {
      g_epi_alpha2sq[1] = inverse_chiprob(epi_ip->epi2, 2);
      g_epi_alpha2sq[2] = inverse_chiprob(epi_ip->epi2, 1);
    }
    if (bigstack_alloc_d(pheno_nm_ct + 1, &g_epi_recip_cache)) {
      goto epistasis_report_ret_NOMEM;
    }
    g_epi_recip_cache[0] = 0.0;
    for (uii = 1; uii <= pheno_nm_ct; uii++) {
      g_epi_recip_cache[uii] = 1.0 / ((double)((int32_t)uii));
    }
  } else {
    if (epi_ip->epi1 == 0.0) {
      dxx = 0.00005;
    } else {
      dxx = epi_ip->epi1 * 0.5;
    }
    dxx = ltqnorm(dxx);
    g_epi_alpha1sq[0] = dxx * dxx;
    dxx = ltqnorm(epi_ip->epi2 / 2);
    g_epi_alpha2sq[0] = dxx * dxx;
  }
  if (!is_fast) {
    if (pheno_d) {
      retval = epistasis_linear_regression(threads, epi_ip, bedfile, bed_offset, unfiltered_marker_ct, marker_reverse, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, marker_uidx_base, marker_ct1, marker_exclude1, marker_idx1_start, marker_idx1_end, marker_ct2, marker_exclude2, is_triangular, job_size, tests_expected, unfiltered_sample_ct, pheno_nm, pheno_nm_ct, pheno_d, parallel_idx, parallel_tot, outname, outname_end, output_min_p, glm_vif_thresh, loadbuf, casebuf, best_chisq, best_ids, n_sig_cts, fail_cts, gap_cts);
    } else {
      retval = epistasis_logistic_regression(threads, epi_ip, bedfile, bed_offset, unfiltered_marker_ct, marker_reverse, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, marker_uidx_base, marker_ct1, marker_exclude1, marker_idx1_start, marker_idx1_end, marker_ct2, marker_exclude2, is_triangular, job_size, tests_expected, unfiltered_sample_ct, pheno_nm, pheno_nm_ct, pheno_c, parallel_idx, parallel_tot, outname, outname_end, output_min_p, loadbuf, casebuf, best_chisq, best_ids, n_sig_cts, fail_cts, gap_cts);
    }
    if (retval) {
      goto epistasis_report_ret_1;
    }
  } else {
    pct_thresh = tests_expected / 100;
    if (is_case_only) {
      g_epi_ctrl_ct = 0;
      ctrl_ctv3 = 0;
      ctrl_ctsplit = 0;
      memcpy(outname_end, ".epi.co", 8);
    } else {
      g_epi_ctrl_ct = ctrl_ct;
      memcpy(outname_end, ".epi.cc", 8);
    }
    if (parallel_tot > 1) {
      outname_end[7] = '.';
      uint32toa_x(parallel_idx + 1, '\0', &(outname_end[8]));
    }
    tot_ctsplit = case_ctsplit + ctrl_ctsplit;
    if (fopen_checked(outname, "w", &outfile)) {
      goto epistasis_report_ret_OPEN_FAIL;
    }
    if (!parallel_idx) {
      wptr = memcpya(g_textbuf, "CHR1 ", 5);
      wptr = fw_strcpyn(plink_maxsnp, 4, "SNP1", wptr);
      wptr = memcpya(wptr, " CHR2 ", 6);
      wptr = fw_strcpyn(plink_maxsnp, 4, "SNP2", wptr);
      wptr = memcpya(wptr, "         STAT ", 14);
      if (is_boost) {
	wptr = memcpya(wptr, "  DF ", 5);
      }
      if (!no_p_value) {
        wptr = memcpya(wptr, "           P ", 13);
      }
      *wptr++ = '\n';
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	goto epistasis_report_ret_WRITE_FAIL;
      }
    }
    // claim up to half of memory with idx1 bufs; each marker currently costs:
    //   (case_ctsplit + ctrl_ctsplit) * sizeof(intptr_t) for loose geno buf
    //   0.25 for missing tracker
    //   sizeof(int32_t) for offset (to skip bottom left triangle, and/or
    //     too-close pairs for case-only tests; will sometimes need to be
    //     larger when sets come into the picture
    //   sizeof(double) for best chisq,
    //   sizeof(int32_t) for best opposite ID,
    //   sizeof(int32_t) for N_SIG count,
    //   sizeof(int32_t) for per-site fail counts, and (bleah)
    //   marker_ct2 * sizeof(double) for the usually oversized results space
    cur_bigstack_left = bigstack_left();
    ulii = 4 * CACHELINE - 3 * sizeof(int32_t) + max_thread_ct * (5 * (CACHELINE - 4));
    if (cur_bigstack_left >= ulii) {
      cur_bigstack_left -= ulii;
    }
    ulii = tot_ctsplit * sizeof(intptr_t) + 4 * sizeof(int32_t) + sizeof(double) + marker_ct2 * sizeof(double);
    idx1_block_size = cur_bigstack_left / (ulii * 2 + 1);
    if (!idx1_block_size) {
      goto epistasis_report_ret_NOMEM;
    }
    if (idx1_block_size > job_size) {
      idx1_block_size = job_size;
    }
    // pad to avoid threads writing to same cacheline
    ulii = (max_thread_ct - 1) * 15 + idx1_block_size;
    // offsets[] isn't really needed, but barely takes any memory
    // if 'case-only', want two more offsets columns to store where the "too
    // close" variants are
    bigstack_alloc_ui(idx1_block_size * 2, &g_epi_geno1_offsets);
    bigstack_alloc_ul(tot_ctsplit * idx1_block_size, &g_epi_geno1);
    bigstack_alloc_ul(QUATERCT_TO_WORDCT(idx1_block_size), &g_epi_zmiss1);
    bigstack_alloc_d(idx1_block_size * marker_ct2, &g_epi_all_chisq);
    bigstack_alloc_d(ulii, &g_epi_best_chisq1);
    bigstack_alloc_ui(ulii, &g_epi_best_id1);
    bigstack_alloc_ui(ulii, &g_epi_n_sig_ct1);
    bigstack_alloc_ui(ulii, &g_epi_fail_ct1);
    for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
      g_epi_geno1[block_idx1 * tot_ctsplit + case_ctv3 - 1] = 0;
      g_epi_geno1[block_idx1 * tot_ctsplit + 2 * case_ctv3 - 1] = 0;
      g_epi_geno1[block_idx1 * tot_ctsplit + case_ctsplit - 1] = 0;
      g_epi_geno1[block_idx1 * tot_ctsplit + case_ctsplit + ctrl_ctv3 - 1] = 0;
      g_epi_geno1[block_idx1 * tot_ctsplit + case_ctsplit + 2 * ctrl_ctv3 - 1] = 0;
      g_epi_geno1[block_idx1 * tot_ctsplit + tot_ctsplit - 1] = 0;
    }
    if (is_triangular) {
      fill_uint_zero(2 * idx1_block_size, g_epi_geno1_offsets);
    }
    // don't actually need best_chisq2, best_id2, n_sig_ct2, fail_ct2 if not
    // triangular, but rather not complicate/duplicate the common case inner
    // loop for now
    ulii = tot_ctsplit * sizeof(intptr_t) + 1 + is_boost * 6 * sizeof(double) + tot_stride * sizeof(int32_t) + max_thread_ct * (3 * sizeof(int32_t) + sizeof(double));
    idx2_block_size = (bigstack_left() - CACHELINE - is_boost * (CACHELINE - 8) - max_thread_ct * (5 * (CACHELINE - 4))) / ulii;
    if (idx2_block_size > marker_ct2) {
      idx2_block_size = marker_ct2;
    }
    idx2_block_size = round_up_pow2(idx2_block_size, 16);
    bigstack_mark2 = g_bigstack_base;
    while (1) {
      if (!idx2_block_size) {
	goto epistasis_report_ret_NOMEM;
      }
      if (!(bigstack_alloc_ul(tot_ctsplit * idx2_block_size, &g_epi_geno2) ||
            bigstack_alloc_ul(QUATERCT_TO_WORDCT(idx2_block_size), &g_epi_zmiss2) ||
	    bigstack_alloc_ui(idx2_block_size * tot_stride, &g_epi_tot2) ||
	    bigstack_alloc_d(max_thread_ct * idx2_block_size, &g_epi_best_chisq2) ||
	    bigstack_alloc_ui(max_thread_ct * idx2_block_size, &g_epi_best_id2) ||
	    bigstack_alloc_ui(max_thread_ct * idx2_block_size, &g_epi_n_sig_ct2) ||
	    bigstack_alloc_ui(max_thread_ct * idx2_block_size, &g_epi_fail_ct2))) {
	if ((!is_boost) || (!bigstack_alloc_d(6 * idx2_block_size, &g_epi_boost_precalc2))) {
	  break;
	}
      }
      bigstack_reset(bigstack_mark2);
      idx2_block_size -= 16;
    }
    for (block_idx2 = 0; block_idx2 < idx2_block_size; block_idx2++) {
      g_epi_geno2[block_idx2 * tot_ctsplit + case_ctv3 - 1] = 0;
      g_epi_geno2[block_idx2 * tot_ctsplit + 2 * case_ctv3 - 1] = 0;
      g_epi_geno2[block_idx2 * tot_ctsplit + case_ctsplit - 1] = 0;
      g_epi_geno2[block_idx2 * tot_ctsplit + case_ctsplit + ctrl_ctv3 - 1] = 0;
      g_epi_geno2[block_idx2 * tot_ctsplit + case_ctsplit + 2 * ctrl_ctv3 - 1] = 0;
      g_epi_geno2[block_idx2 * tot_ctsplit + tot_ctsplit - 1] = 0;
    }
    marker_uidx = next_unset_ul_unsafe(marker_exclude1, marker_uidx_base);
    if (marker_idx1) {
      marker_uidx = jump_forward_unset_unsafe(marker_exclude1, marker_uidx + 1, marker_idx1);
    }
    wptr = memcpya(g_logbuf, "--fast-epistasis", 16);
    if (is_boost) {
      wptr = memcpya(wptr, " boost", 6);
    } else if (no_ueki) {
      wptr = memcpya(wptr, " no-ueki", 8);
    } else if (do_joint_effects) {
      wptr = memcpya(wptr, " joint-effects", 14);
    }
    if (is_case_only) {
      wptr = memcpya(wptr, " case-only", 10);
    }
    wptr = memcpya(wptr, " to ", 4);
    wptr = strcpya(wptr, outname);
    memcpy(wptr, " ... ", 6);
    wordwrapb(16); // strlen("99% [processing]")
    logprintb();
    fputs("0%", stdout);
    do {
      fputs(" [processing]", stdout);
      fflush(stdout);
      if (idx1_block_size > marker_idx1_end - marker_idx1) {
        idx1_block_size = marker_idx1_end - marker_idx1;
        if (idx1_block_size < max_thread_ct) {
	  max_thread_ct = idx1_block_size;
	  g_epi_thread_ct = max_thread_ct;
	}
      }
      g_epi_marker_idx1 = marker_idx1;
      dptr = g_epi_all_chisq;
      dptr2 = &(g_epi_all_chisq[idx1_block_size * marker_ct2]);
      do {
	*dptr++ = -1;
      } while (dptr < dptr2);
      marker_uidx_tmp = marker_uidx;
      if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	goto epistasis_report_ret_READ_FAIL;
      }
      cur_workload = idx1_block_size * marker_ct2;
      if (is_triangular) {
	for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
	  ulii = block_idx1 + marker_idx1 + 1;
	  cur_workload -= ulii;
	  // edit this during loading, when we have to know marker_uidx anyway,
	  // if case-only
	  g_epi_geno1_offsets[2 * block_idx1 + 1] = ulii;
	}
      } else {
        fill_uint_zero(2 * idx1_block_size, g_epi_geno1_offsets);
	marker_uidx2 = marker_uidx_base;
	marker_idx2 = 0;
      }
      tests_complete += cur_workload;
      ulii = 0; // total number of tests
      g_epi_idx1_block_bounds[0] = 0;
      g_epi_idx1_block_bounds16[0] = 0;
      block_idx1 = 0;
      for (tidx = 1; tidx < max_thread_ct; tidx++) {
	uljj = (((uint64_t)cur_workload) * tidx) / max_thread_ct;
	if (is_triangular) {
	  do {
	    // slightly inaccurate for case-only due to the way --gap is
	    // supported, but this doesn't affect any calculation results, only
	    // the progress display
	    ulii += marker_ct2 - g_epi_geno1_offsets[2 * block_idx1 + 1];
	    block_idx1++;
	  } while (ulii < uljj);
	} else {
	  do {
	    ulii += marker_ct2;
	    block_idx1++;
	  } while (ulii < uljj);
	}
	uii = block_idx1 - g_epi_idx1_block_bounds[tidx - 1];
        g_epi_idx1_block_bounds[tidx] = block_idx1;
        g_epi_idx1_block_bounds16[tidx] = g_epi_idx1_block_bounds16[tidx - 1] + round_up_pow2_ui(uii, 16);
      }
      g_epi_idx1_block_bounds[max_thread_ct] = idx1_block_size;
      fill_ulong_zero(QUATERCT_TO_WORDCT(idx1_block_size), g_epi_zmiss1);
      chrom_end = 0;
      for (block_idx1 = 0; block_idx1 < idx1_block_size; marker_uidx_tmp++, block_idx1++) {
        if (IS_SET(marker_exclude1, marker_uidx_tmp)) {
	  marker_uidx_tmp = next_unset_ul_unsafe(marker_exclude1, marker_uidx_tmp);
          if (fseeko(bedfile, bed_offset + (marker_uidx_tmp * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	    goto epistasis_report_ret_READ_FAIL;
	  }
	}
	if (load_and_split3(bedfile, loadbuf, unfiltered_sample_ct, &(g_epi_geno1[block_idx1 * tot_ctsplit]), pheno_nm, pheno_c, case_ctv3, ctrl_ctv3, IS_SET(marker_reverse, marker_uidx_tmp), is_case_only, &ulii)) {
	  goto epistasis_report_ret_READ_FAIL;
	}
	if (ulii) {
	  g_epi_zmiss1[block_idx1 / BITCT2] |= ulii << (2 * (block_idx1 % BITCT2));
	  // g_epi_tot1 doesn't need to exist, better for each thread to
	  // determine those totals on the fly
	}
	if (is_case_only_window) {
	  cur_window_end = marker_pos[marker_uidx_tmp] + case_only_gap;
	  if (marker_uidx_tmp >= chrom_end) {
	    chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx_tmp);
	    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
	    if (is_triangular) {
	      marker_uidx2 = marker_uidx_tmp;
	      marker_idx2 = block_idx1 + marker_idx1;
	      last_pos = marker_pos[marker_uidx_tmp];
	    } else {
	      uii = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
	      if (marker_pos[marker_uidx_tmp] < case_only_gap) {
		ujj = 0;
	      } else {
		ujj = marker_pos[marker_uidx_tmp] + 1 - case_only_gap;
	      }
	      marker_uidx2_trail = next_unset(marker_exclude2, uii + uint32arr_greater_than(&(marker_pos[uii]), marker_uidx_tmp + 1 - uii, ujj), chrom_end);
	      marker_idx2_trail = marker_uidx2_trail - popcount_bit_idx(marker_exclude2, 0, marker_uidx2_trail);
	      if (marker_uidx2_trail < chrom_end) {
		first_pos = marker_pos[marker_uidx2_trail];
		// this could be more efficient, but not a big deal since
		// there aren't many chromosomes
	        marker_uidx2 = next_unset(marker_exclude2, uii + uint32arr_greater_than(&(marker_pos[marker_uidx_tmp]), chrom_end - marker_uidx_tmp, cur_window_end), chrom_end);
	      } else {
		first_pos = 0x7fffffffU;
		marker_uidx2 = chrom_end;
	      }
	      marker_idx2 = marker_idx2_trail + marker_uidx2 - marker_uidx2_trail - popcount_bit_idx(marker_exclude2, marker_uidx2_trail, marker_uidx2);
	      if (marker_uidx2 < chrom_end) {
		last_pos = marker_pos[marker_uidx2];
	      } else {
		last_pos = 0xffffffffU;
	      }
	    }
	  }
	  while (last_pos < cur_window_end) {
	    marker_idx2++;
	    marker_uidx2++;
	    next_unset_ul_ck(marker_exclude2, chrom_end, &marker_uidx2);
	    if (marker_uidx2 != chrom_end) {
	      last_pos = marker_pos[marker_uidx2];
	    } else {
	      last_pos = 0xffffffffU;
	    }
	  }
	  if (is_triangular) {
	    ulii = block_idx1 + marker_idx1;
            gap_cts[ulii] += marker_idx2 - ulii - 1;
	    while (++ulii < marker_idx2) {
	      gap_cts[ulii] += 1;
	    }
	    g_epi_geno1_offsets[2 * block_idx1 + 1] = marker_idx2;
	  } else {
	    uii = marker_pos[marker_uidx_tmp];
	    while (first_pos + case_only_gap <= uii) {
	      marker_idx2_trail++;
	      marker_uidx2_trail++;
	      next_unset_ul_ck(marker_exclude2, chrom_end, &marker_uidx2_trail);
              if (marker_uidx2_trail != chrom_end) {
		first_pos = marker_pos[marker_uidx2_trail];
	      } else {
		first_pos = 0x7fffffffU;
	      }
	    }
	    if (marker_idx2 > marker_idx2_trail) {
	      g_epi_geno1_offsets[2 * block_idx1] = marker_idx2_trail;
	      g_epi_geno1_offsets[2 * block_idx1 + 1] = marker_idx2;
	      gap_cts[block_idx1 + marker_idx1] = marker_idx2 - marker_idx2_trail;
	    }
	  }
	} else if (!is_triangular) {
          if (!IS_SET(marker_exclude2, marker_uidx_tmp)) {
	    // do not compare against self
	    marker_idx2 += marker_uidx_tmp - marker_uidx2 - popcount_bit_idx(marker_exclude2, marker_uidx2, marker_uidx_tmp);
	    marker_uidx2 = marker_uidx_tmp;
	    g_epi_geno1_offsets[2 * block_idx1] = marker_idx2;
	    g_epi_geno1_offsets[2 * block_idx1 + 1] = marker_idx2 + 1;
	    gap_cts[block_idx1 + marker_idx1] = 1;
	  }
	}
      }
      marker_uidx2 = next_unset_ul_unsafe(marker_exclude2, marker_uidx_base);
      if (is_triangular) {
	marker_idx2 = marker_idx1 + 1;
        marker_uidx2 = jump_forward_unset_unsafe(marker_exclude2, marker_uidx2 + 1, marker_idx2);
      } else {
        marker_idx2 = 0;
      }
      if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	goto epistasis_report_ret_READ_FAIL;
      }
      cur_idx2_block_size = idx2_block_size;
      do {
	if (cur_idx2_block_size > marker_ct2 - marker_idx2) {
	  cur_idx2_block_size = marker_ct2 - marker_idx2;
	}
	fill_ulong_zero(QUATERCT_TO_WORDCT(cur_idx2_block_size), g_epi_zmiss2);
        for (block_idx2 = 0; block_idx2 < cur_idx2_block_size; marker_uidx2++, block_idx2++) {
          if (IS_SET(marker_exclude2, marker_uidx2)) {
	    marker_uidx2 = next_unset_ul_unsafe(marker_exclude2, marker_uidx2);
            if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	      goto epistasis_report_ret_READ_FAIL;
	    }
	  }
	  ulptr = &(g_epi_geno2[block_idx2 * tot_ctsplit]);
	  if (load_and_split3(bedfile, loadbuf, unfiltered_sample_ct, ulptr, pheno_nm, pheno_c, case_ctv3, ctrl_ctv3, IS_SET(marker_reverse, marker_uidx2), is_case_only, &ulii)) {
	    goto epistasis_report_ret_READ_FAIL;
	  }
	  uiptr = &(g_epi_tot2[block_idx2 * tot_stride]);
	  uiptr[0] = popcount_longs(ulptr, case_ctv3);
	  uiptr[1] = popcount_longs(&(ulptr[case_ctv3]), case_ctv3);
	  uiptr[2] = popcount_longs(&(ulptr[2 * case_ctv3]), case_ctv3);
	  if (!is_case_only) {
	    ulptr = &(ulptr[case_ctv3 * 3]);
	    uiptr[3] = popcount_longs(ulptr, ctrl_ctv3);
	    uiptr[4] = popcount_longs(&(ulptr[ctrl_ctv3]), ctrl_ctv3);
	    uiptr[5] = popcount_longs(&(ulptr[2 * ctrl_ctv3]), ctrl_ctv3);
	    if (is_boost) {
	      boost_calc_p_bc(uiptr[0], uiptr[1], uiptr[2], uiptr[3], uiptr[4], uiptr[5], &(g_epi_boost_precalc2[block_idx2 * 6]));
	    }
	  }
	  if (ulii) {
	    g_epi_zmiss2[block_idx2 / BITCT2] |= ulii << (2 * (block_idx2 % BITCT2));
	  }
	}
	g_epi_idx2_block_size = cur_idx2_block_size;
	g_epi_idx2_block_start = marker_idx2;
	idx2_block_sizea16 = round_up_pow2(cur_idx2_block_size, 16);
        fill_uint_zero(idx1_block_size + 15 * (max_thread_ct - 1), g_epi_n_sig_ct1);
	fill_uint_zero(idx1_block_size + 15 * (max_thread_ct - 1), g_epi_fail_ct1);
        fill_uint_zero(idx2_block_sizea16 * max_thread_ct, g_epi_n_sig_ct2);
	fill_uint_zero(idx2_block_sizea16 * max_thread_ct, g_epi_fail_ct2);
	for (tidx = 0; tidx < max_thread_ct; tidx++) {
	  ulii = g_epi_idx1_block_bounds[tidx];
	  uljj = g_epi_idx1_block_bounds[tidx + 1];
	  memcpy(&(g_epi_best_chisq1[g_epi_idx1_block_bounds16[tidx]]), &(g_epi_all_chisq[marker_idx1 + ulii]), (uljj - ulii) * sizeof(double));
	  ulii = g_epi_geno1_offsets[2 * ulii + 1];
	  if (ulii < marker_idx2 + cur_idx2_block_size) {
	    if (ulii <= marker_idx2) {
	      ulii = 0;
	    } else {
	      ulii -= marker_idx2;
	    }
	    memcpy(&(g_epi_best_chisq2[tidx * idx2_block_sizea16 + ulii]), &(g_epi_all_chisq[marker_idx2 + ulii]), (cur_idx2_block_size - ulii) * sizeof(double));
	  }
	  // no need to initialize IDs since they are only referenced when a
	  // higher chisq value is present, and when that happens an ID is
          // always written
	}
	is_last_block = (marker_idx2 + cur_idx2_block_size >= marker_ct2);
	if (spawn_threads2(threads, &fast_epi_thread, max_thread_ct, is_last_block)) {
	  goto epistasis_report_ret_THREAD_CREATE_FAIL;
	}
	fast_epi_thread((void*)0);
	join_threads2(threads, max_thread_ct, is_last_block);
	// merge best_chisq, best_ids, fail_cts
	for (tidx = 0; tidx < max_thread_ct; tidx++) {
	  ulii = g_epi_idx1_block_bounds[tidx];
	  uljj = g_epi_idx1_block_bounds[tidx + 1] - ulii;
	  uii = g_epi_idx1_block_bounds16[tidx];
	  dptr = &(g_epi_best_chisq1[uii]);
	  uiptr = &(g_epi_best_id1[uii]);
	  uiptr2 = &(g_epi_n_sig_ct1[uii]);
	  uiptr3 = &(g_epi_fail_ct1[uii]);
	  ulii += marker_idx1;
          dptr2 = &(best_chisq[ulii]);
          uiptr4 = &(n_sig_cts[ulii]);
          uiptr5 = &(fail_cts[ulii]);
	  for (block_idx1 = 0; block_idx1 < uljj; block_idx1++, dptr2++, uiptr4++, uiptr5++) {
	    dxx = *dptr++;
	    if (dxx > (*dptr2)) {
	      *dptr2 = dxx;
	      best_ids[block_idx1 + ulii] = uiptr[block_idx1];
	    }
            *uiptr4 += *uiptr2++;
            *uiptr5 += *uiptr3++;
	  }
	}
	if (is_triangular) {
	  for (tidx = 0; tidx < max_thread_ct; tidx++) {
	    block_idx2 = g_epi_geno1_offsets[2 * g_epi_idx1_block_bounds[tidx] + 1];
	    if (block_idx2 <= marker_idx2) {
	      block_idx2 = 0;
	    } else {
	      block_idx2 -= marker_idx2;
	    }
	    dptr = &(g_epi_best_chisq2[tidx * idx2_block_sizea16 + block_idx2]);
	    uiptr = &(g_epi_best_id2[tidx * idx2_block_sizea16]);
	    uiptr2 = &(g_epi_n_sig_ct2[tidx * idx2_block_sizea16 + block_idx2]);
	    uiptr3 = &(g_epi_fail_ct2[tidx * idx2_block_sizea16 + block_idx2]);
	    dptr2 = &(best_chisq[block_idx2 + marker_idx2]);
	    uiptr4 = &(n_sig_cts[block_idx2 + marker_idx2]);
	    uiptr5 = &(fail_cts[block_idx2 + marker_idx2]);
	    for (; block_idx2 < cur_idx2_block_size; block_idx2++, dptr2++, uiptr4++, uiptr5++) {
	      dxx = *dptr++;
	      if (dxx > (*dptr2)) {
		*dptr2 = dxx;
		best_ids[block_idx2 + marker_idx2] = uiptr[block_idx2];
	      }
	      *uiptr4 += *uiptr2++;
	      *uiptr5 += *uiptr3++;
	    }
	  }
	}
        marker_idx2 += cur_idx2_block_size;
      } while (marker_idx2 < marker_ct2);
      fputs("\b\b\b\b\b\b\b\b\b\b\bwriting]   \b\b\b", stdout);
      fflush(stdout);
      chrom_end = 0;
      block_idx1 = 0;
      while (1) {
	next_unset_ul_unsafe_ck(marker_exclude1, &marker_uidx);
	ujj = g_epi_geno1_offsets[2 * block_idx1];
	marker_idx2 = 0;
	dptr = &(g_epi_all_chisq[block_idx1 * marker_ct2]);
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
	}
        wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf));
	*wptr_start++ = ' ';
	wptr_start = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
	*wptr_start++ = ' ';
	marker_uidx2 = next_unset_ul_unsafe(marker_exclude2, marker_uidx_base);
	for (chrom_fo_idx2 = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx2); chrom_fo_idx2 < chrom_ct; chrom_fo_idx2++) {
	  chrom_end2 = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx2 + 1];
	  if (marker_uidx2 >= chrom_end2) {
	    continue;
	  }
          chrom_idx2 = chrom_info_ptr->chrom_file_order[chrom_fo_idx2];
          wptr_start2 = width_force(4, wptr_start, chrom_name_write(chrom_info_ptr, chrom_idx2, wptr_start));
	  *wptr_start2++ = ' ';
	  for (; marker_uidx2 < chrom_end2; ++marker_uidx2, next_unset_ul_ck(marker_exclude2, unfiltered_marker_ct, &marker_uidx2), ++marker_idx2, ++dptr) {
	    if (marker_idx2 == ujj) {
	      marker_idx2 = g_epi_geno1_offsets[2 * block_idx1 + 1];
	      if (marker_idx2 == marker_ct2) {
		goto epistasis_report_write_loop;
	      }
	      if (marker_idx2 > ujj) {
	        marker_uidx2 = jump_forward_unset_unsafe(marker_exclude2, marker_uidx2 + 1, marker_idx2 - ujj);
	        dptr = &(dptr[marker_idx2 - ujj]);
	        if (marker_uidx2 >= chrom_end2) {
		  break;
	        }
	      }
	    } else if (marker_idx2 == marker_ct2) {
	      goto epistasis_report_write_loop;
	    }
	    dxx = *dptr;
	    if (dxx != -1) {
	      wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr_start2);
	      *wptr++ = ' ';
	      if (is_boost) {
		if (dxx == dxx) { // not nan
		  memcpy(ularr, &dxx, sizeof(double));
		  uii = 4 >> (ularr[0] & 3);
		  // don't want ugly e-324s when zero belongs
		  ularr[0] &= ~(3 * ONELU);
		  memcpy(&dxx, ularr, sizeof(double));
		  wptr = width_force(12, wptr, dtoa_g(dxx, wptr));
		  wptr = memseta(wptr, 32, 4);
		  *wptr++ = '0' + uii;
		  *wptr++ = ' ';
		} else {
                  wptr = memcpya(wptr, "         nan    0 ", 18);
		  uii = 0;
		}
	      } else if (!no_ueki) {
		wptr = width_force(12, wptr, dtoa_g(dxx, wptr));
		*wptr++ = ' ';
	      } else {
		// lower precision compatibility mode
                wptr = dtoa_g_wxp4x(dxx, 12, ' ', wptr);
	      }
	      if (!no_p_value) {
		if (!is_boost) {
		  dxx = normdist(-sqrt(dxx)) * 2;
		  wptr = dtoa_g_wxp4x(MAXV(dxx, output_min_p), 12, ' ', wptr);
		} else if (uii) {
		  dxx = chiprob_p(dxx, uii);
		  wptr = dtoa_g_wxp4x(MAXV(dxx, output_min_p), 12, ' ', wptr);
		} else {
		  wptr = memcpya(wptr, "          NA ", 13);
		}
	      }
	      *wptr++ = '\n';
	      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
		goto epistasis_report_ret_WRITE_FAIL;
	      }
	      // could remove this writeback in --epi1 1 case
	      *dptr = -1;
	    }
	  }
	}
      epistasis_report_write_loop:
	block_idx1++;
	marker_uidx++;
	if (block_idx1 >= idx1_block_size) {
	  break;
	}
      }
      marker_idx1 += idx1_block_size;
      fputs("\b\b\b\b\b\b\b\b\b\b          \b\b\b\b\b\b\b\b\b\b", stdout);
      if (tests_complete >= pct_thresh) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	pct = (tests_complete * 100LLU) / tests_expected;
	if (pct < 100) {
	  printf("\b\b%" PRIuPTR "%%", pct);
	  fflush(stdout);
	  pct_thresh = ((++pct) * ((uint64_t)tests_expected)) / 100;
	}
      }
    } while (marker_idx1 < marker_idx1_end);
    if (fclose_null(&outfile)) {
      goto epistasis_report_ret_WRITE_FAIL;
    }
  }
  memcpy(&(outname_end[7]), ".summary", 9);
  if (parallel_tot > 1) {
    outname_end[15] = '.';
    uint32toa_x(parallel_idx + 1, '\0', &(outname_end[16]));
  }
  if (fopen_checked(outname, "w", &outfile)) {
    goto epistasis_report_ret_OPEN_FAIL;
  }
  wptr = memcpya(g_textbuf, " CHR ", 5);
  wptr = fw_strcpyn(plink_maxsnp, 3, "SNP", wptr);
  if (parallel_tot == 1) {
    wptr = strcpya(wptr, "        N_SIG        N_TOT         PROP   BEST_CHISQ BEST_CHR ");
  } else {
    wptr = strcpya(wptr, "        N_SIG        N_TOT   BEST_CHISQ BEST_CHR ");
  }
  wptr = fw_strcpyn(plink_maxsnp, 8, "BEST_SNP", wptr);
  wptr = memcpya(wptr, " \n", 2);
  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
    goto epistasis_report_ret_WRITE_FAIL;
  }
  bigstack_reset(bigstack_mark3);
  if (bigstack_alloc_ui(marker_ct1, &marker_idx_to_uidx)) {
    goto epistasis_report_ret_NOMEM;
  }
  fill_idx_to_uidx(marker_exclude2, unfiltered_marker_ct, marker_ct2, marker_idx_to_uidx);
  marker_idx1 = marker_idx1_start;
  marker_uidx = next_unset_ul_unsafe(marker_exclude1, marker_uidx_base);
  if (marker_idx1) {
    marker_uidx = jump_forward_unset_unsafe(marker_exclude1, marker_uidx + 1, marker_idx1);
  }
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    if (marker_uidx >= chrom_end) {
      continue;
    }
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf));
    *wptr_start++ = ' ';
    for (; marker_uidx < chrom_end; marker_uidx++, next_unset_ul_ck(marker_exclude1, unfiltered_marker_ct, &marker_uidx), marker_idx1++) {
      uii = n_sig_cts[marker_idx1];
      ujj = fail_cts[marker_idx1];
      if (gap_cts) {
	ujj += gap_cts[marker_idx1];
      }
      tests_thrown_out += (uint64_t)ujj;
      // number of tests attempted in this run:
      // * if set1 and set2 are identical, there are
      //   marker_ct2 - 1 - marker_idx1_start cells between the row and the
      //   same-index column
      // * otherwise, gap_cts[] counted the number of skipped cells
      if (marker_idx1 < marker_idx1_end) {
	if (is_triangular) {
	  ujj = marker_ct2 - 1 - marker_idx1_start - ujj;
	} else {
	  ujj = marker_ct2 - ujj;
	}
      } else {
	// --parallel bugfix
	ujj = job_size - ujj;
      }
      wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
      wptr = memseta(wptr, 32, 3);
      wptr = uint32toa_w10(uii, wptr);
      wptr = memseta(wptr, 32, 3);
      wptr = uint32toa_w10x(ujj, ' ', wptr);
      if (parallel_tot == 1) {
        wptr = dtoa_g_wxp4x(((double)((int32_t)uii)) / ((double)((int32_t)ujj)), 12, ' ', wptr);
      }
      if (ujj) {
	if (parallel_tot == 1) {
	  // or cat mode
	  wptr = dtoa_g_wxp4x(best_chisq[marker_idx1], 12, ' ', wptr);
	} else {
	  // greater precision for accurate merges
	  wptr = dtoa_g_wxp8x(best_chisq[marker_idx1], 12, ' ', wptr);
	}
	uii = marker_idx_to_uidx[best_ids[marker_idx1]];
	wptr = width_force(4, wptr, chrom_name_write(chrom_info_ptr, get_variant_chrom(chrom_info_ptr, uii), wptr));
	*wptr++ = ' ';
	wptr = fw_strcpy(plink_maxsnp, &(marker_ids[uii * max_marker_id_len]), wptr);
      } else {
	wptr = memcpya(wptr, "          NA   NA", 17);
	wptr = memseta(wptr, 32, plink_maxsnp - 1);
	wptr = memcpya(wptr, "NA", 2);
      }
      wptr = memcpya(wptr, " \n", 2);
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	goto epistasis_report_ret_WRITE_FAIL;
      }
    }
  }
  if (is_triangular) {
    tests_thrown_out /= 2; // all fails double-counted in triangle case
  }
  fputs("\b\b", stdout);
  LOGPRINTF("done.\n");
  LOGPRINTFWW("%" PRIu64 " valid test%s performed, summary written to %s .\n", tests_expected - tests_thrown_out, (tests_expected - tests_thrown_out == 1)? "" : "s", outname);

  while (0) {
  epistasis_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  epistasis_report_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  epistasis_report_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  epistasis_report_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  epistasis_report_ret_TOO_FEW_MARKERS:
    if (pheno_d) {
      if (is_triangular) {
        logerrprint("Error: --epistasis requires 2+ non-monomorphic autosomal diploid loci.\n");
      } else {
        logerrprint("Error: Each --epistasis set must contain at least one non-monomorphic autosomal\ndiploid site.\n");
      }
    } else {
      if (is_triangular) {
        logerrprint("Error: --{fast-}epistasis requires 2+ autosomal diploid loci not monomorphic in\neither cases or controls.\n");
      } else {
        logerrprint("Error: Each --{fast-}epistasis set must contain at least one autosomal diploid\nlocus not monomorphic in either cases or controls.\n");
      }
    }
    retval = RET_INVALID_CMDLINE;
    break;
  epistasis_report_ret_INVALID_CMDLINE_2:
    logerrprintb();
  epistasis_report_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  epistasis_report_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 epistasis_report_ret_1:
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t indep_pairphase(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  // Like ld_prune(), except that it computes the full 3x3 contingency table,
  // and is always in pairwise mode.
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct);
  uintptr_t founder_ct = popcount_longs(founder_info, unfiltered_sample_ctv2 / 2);
  uintptr_t founder_ctl = BITCT_TO_WORDCT(founder_ct);
  uintptr_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(founder_ct);
  // no actual case/control split here, but keep the variables the same to
  // minimize divergence from ld_report_dprime()
  uintptr_t founder_ctsplit = 3 * founder_ctv3;
  uintptr_t final_mask = get_final_mask(founder_ct);
  uintptr_t window_max = 1;
  uintptr_t* founder_include2 = nullptr;
  uintptr_t* founder_male_include2 = nullptr;
  uintptr_t* sex_male_collapsed = nullptr;
  uintptr_t* cur_geno1_male = nullptr;
  double prune_ld_thresh = ldip->prune_last_param * (1 + SMALL_EPSILON);
  uint32_t window_is_kb = (ldip->modifier / LD_PRUNE_KB_WINDOW) & 1;
  uint32_t ld_window_size = ldip->prune_window_size;
  uint32_t ld_window_incr = ldip->prune_window_incr;
  uint32_t tot_exclude_ct = 0;
  uint32_t at_least_one_prune = 0;
  uint32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  int32_t retval = 0;
  uint32_t tot1[6];
  uint32_t counts[18];
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* dummy_nm;
  uintptr_t* pruned_arr;
  uintptr_t* geno;
  uintptr_t* zmiss;
  uintptr_t* cur_geno1;
  uintptr_t* cur_geno2;
  uint32_t* live_indices;
  uint32_t* start_arr;
  uint32_t* cur_tots;
  uint32_t* cur_tot2;
  uintptr_t window_maxl;
  uintptr_t cur_exclude_ct;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  double freq1x;
  double freq2x;
  double freqx1;
  double freqx2;
  double freq11;
  double freq11_expected;
  double rsq;
  uint32_t pct_thresh;
  uint32_t window_unfiltered_start;
  uint32_t window_unfiltered_end;
  uint32_t cur_window_size;
  uint32_t cur_chrom;
  uint32_t chrom_start;
  uint32_t chrom_end;
  uint32_t is_haploid;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t prev_end;
  uint32_t nm_fixed;
  uint32_t cur_zmiss2;
  uint32_t pct;
  uint32_t uii;
  if (founder_ct < 2) {
    logerrprint("Warning: Skipping --indep-pairphase since there are less than two founders.\n(--make-founders may come in handy here.)\n");
    goto indep_pairphase_ret_1;
  }
  if (is_set(chrom_info_ptr->chrom_mask, 0)) {
    ulii = count_chrom_markers(chrom_info_ptr, marker_exclude, 0);
    if (chrom_info_ptr->zero_extra_chroms) {
      for (uii = chrom_info_ptr->max_code + 1; uii < chrom_code_end; uii++) {
	ulii += count_chrom_markers(chrom_info_ptr, marker_exclude, uii);
      }
      chrom_code_end = chrom_info_ptr->max_code + 1;
    }
    marker_ct -= ulii;
    LOGPRINTF("--indep-pairphase: Ignoring %" PRIuPTR " chromosome 0 variant%s.\n", ulii, (ulii == 1)? "" : "s");
  }
  if (marker_ct < 2) {
    logerrprint("Error: Too few variants for --indep-pairphase.\n");
    goto indep_pairphase_ret_INVALID_FORMAT;
  }

  // no need to force founder_male_include2 initialization here
  if (alloc_collapsed_haploid_filters(founder_info, sex_male, unfiltered_sample_ct, founder_ct, hh_exists, 1, &founder_include2, &founder_male_include2)) {
    goto indep_pairphase_ret_NOMEM;
  }

  if (window_is_kb) {
    // determine maximum number of markers that may need to be loaded at once
    for (cur_chrom = 1; cur_chrom < chrom_code_end; cur_chrom++) {
      if (is_set(chrom_info_ptr->chrom_mask, cur_chrom)) {
	window_max = chrom_window_max(marker_pos, marker_exclude, chrom_info_ptr, cur_chrom, 0x7fffffff, ld_window_size * 1000, window_max);
      }
    }
  }

  window_unfiltered_start = ld_prune_next_valid_chrom_start(marker_exclude, 0, chrom_info_ptr, chrom_code_end, unfiltered_marker_ct);

  if (bigstack_alloc_ul(unfiltered_marker_ctl, &pruned_arr)) {
    goto indep_pairphase_ret_NOMEM;
  }

  memcpy(pruned_arr, marker_exclude, unfiltered_marker_ctl * sizeof(intptr_t));

  if (!window_is_kb) {
    window_max = ld_window_size;
  }
  window_maxl = BITCT_TO_WORDCT(window_max);
  if (bigstack_alloc_ui(window_max, &live_indices) ||
      bigstack_alloc_ui(window_max, &start_arr) ||
      bigstack_alloc_ul(unfiltered_sample_ctv2, &loadbuf_raw) ||
      bigstack_alloc_ul(founder_ctl * 2, &loadbuf) ||
      bigstack_alloc_ul(founder_ctl, &dummy_nm) ||
      bigstack_alloc_ul(founder_ctsplit * window_max, &geno) ||
      bigstack_alloc_ul(window_maxl, &zmiss) ||
      bigstack_alloc_ui(window_max * 3, &cur_tots)) {
    goto indep_pairphase_ret_NOMEM;
  }
  loadbuf[founder_ctl * 2 - 2] = 0;
  loadbuf[founder_ctl * 2 - 1] = 0;
  fill_all_bits(founder_ct, dummy_nm);
  // bugfix: this loop must start at 0, not 1
  for (ulii = 0; ulii < window_max; ulii++) {
    geno[ulii * founder_ctsplit + founder_ctv3 - 1] = 0;
    geno[ulii * founder_ctsplit + 2 * founder_ctv3 - 1] = 0;
    geno[ulii * founder_ctsplit + founder_ctsplit - 1] = 0;
  }
  if ((chrom_info_ptr->xymt_codes[X_OFFSET] != -2) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->xymt_codes[X_OFFSET])) {
    if (bigstack_alloc_ul(founder_ctl, &sex_male_collapsed) ||
        bigstack_alloc_ul(founder_ctsplit, &cur_geno1_male)) {
      goto indep_pairphase_ret_NOMEM;
    }
    copy_bitarr_subset(sex_male, founder_info, unfiltered_sample_ct, founder_ct, sex_male_collapsed);
  }
  do {
    prev_end = 0;
    ld_prune_start_chrom(window_is_kb, &cur_chrom, &chrom_end, window_unfiltered_start, live_indices, start_arr, &window_unfiltered_end, ld_window_size, &cur_window_size, unfiltered_marker_ct, pruned_arr, chrom_info_ptr, marker_pos, &is_haploid, &is_x, &is_y);
    cur_exclude_ct = 0;
    fill_ulong_zero(window_maxl, zmiss);
    if (cur_window_size > 1) {
      for (ulii = 0; ulii < (uintptr_t)cur_window_size; ulii++) {
	uljj = live_indices[ulii];
	if (fseeko(bedfile, bed_offset + (uljj * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	  goto indep_pairphase_ret_READ_FAIL;
	}
	if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, uljj), bedfile, loadbuf_raw, loadbuf)) {
	  goto indep_pairphase_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)loadbuf);
	}
	cur_geno1 = &(geno[ulii * founder_ctsplit]);
	load_and_split3(nullptr, loadbuf, founder_ct, cur_geno1, dummy_nm, dummy_nm, founder_ctv3, 0, 0, 1, &ulkk);
	cur_tots[ulii * 3] = popcount_longs(cur_geno1, founder_ctv3);
	cur_tots[ulii * 3 + 1] = popcount_longs(&(cur_geno1[founder_ctv3]), founder_ctv3);
	cur_tots[ulii * 3 + 2] = popcount_longs(&(cur_geno1[2 * founder_ctv3]), founder_ctv3);
	if ((!cur_tots[ulii * 3 + 1]) && ((!cur_tots[ulii * 3]) || (!cur_tots[ulii * 3 + 2]))) {
	  SET_BIT(uljj, pruned_arr);
	  cur_exclude_ct++;
	} else if (ulkk == 3) {
	  SET_BIT(ulii, zmiss);
	}
      }
    }
    pct = 1;
    chrom_start = get_chrom_start_vidx(chrom_info_ptr, cur_chrom);
    pct_thresh = window_unfiltered_start + ((uint64_t)pct * (chrom_end - chrom_start)) / 100;
    while ((window_unfiltered_start < chrom_end) || (cur_window_size > 1)) {
      if (cur_window_size > 1) {
	do {
	  at_least_one_prune = 0;
	  for (ulii = 0; ulii < cur_window_size - 1; ulii++) {
	    if (IS_SET(pruned_arr, live_indices[ulii])) {
	      continue;
	    }
	    uljj = ulii + 1;
	    while (live_indices[uljj] < start_arr[ulii]) {
	      if (++uljj == cur_window_size) {
		goto indep_pairphase_skip_marker;
	      }
	    }
	    cur_geno1 = &(geno[ulii * founder_ctsplit]);
	    memcpy(tot1, &(cur_tots[ulii * 3]), 3 * sizeof(int32_t));
	    nm_fixed = is_set_ul(zmiss, ulii);
	    if (is_x) {
	      memcpy(cur_geno1_male, cur_geno1, founder_ctsplit * sizeof(intptr_t));
              bitvec_and(sex_male_collapsed, founder_ctv3, cur_geno1_male);
	      tot1[3] = popcount_longs(cur_geno1_male, founder_ctv3);
              bitvec_and(sex_male_collapsed, founder_ctv3, &(cur_geno1_male[founder_ctv3]));
	      tot1[4] = popcount_longs(&(cur_geno1_male[founder_ctv3]), founder_ctv3);
              bitvec_and(sex_male_collapsed, founder_ctv3, &(cur_geno1_male[2 * founder_ctv3]));
	      tot1[5] = popcount_longs(&(cur_geno1_male[2 * founder_ctv3]), founder_ctv3);
	    }
	    for (; uljj < cur_window_size; uljj++) {
	      if (IS_SET(pruned_arr, live_indices[uljj])) {
		continue;
	      }
	      cur_geno2 = &(geno[uljj * founder_ctsplit]);
	      cur_tot2 = &(cur_tots[uljj * 3]);
	      cur_zmiss2 = IS_SET(zmiss, uljj);
	      if (nm_fixed) {
		two_locus_count_table_zmiss1(cur_geno1, cur_geno2, counts, founder_ctv3, cur_zmiss2);
		if (cur_zmiss2) {
                  counts[2] = tot1[0] - counts[0] - counts[1];
                  counts[5] = tot1[1] - counts[3] - counts[4];
		}
                counts[6] = cur_tot2[0] - counts[0] - counts[3];
                counts[7] = cur_tot2[1] - counts[1] - counts[4];
                counts[8] = cur_tot2[2] - counts[2] - counts[5];
	      } else {
                two_locus_count_table(cur_geno1, cur_geno2, counts, founder_ctv3, cur_zmiss2);
                if (cur_zmiss2) {
                  counts[2] = tot1[0] - counts[0] - counts[1];
                  counts[5] = tot1[1] - counts[3] - counts[4];
                  counts[8] = tot1[2] - counts[6] - counts[7];
		}
	      }
              if (is_x) {
                two_locus_count_table(cur_geno1_male, cur_geno2, &(counts[9]), founder_ctv3, cur_zmiss2);
                if (cur_zmiss2) {
                  counts[11] = tot1[3] - counts[9] - counts[10];
                  counts[14] = tot1[4] - counts[12] - counts[13];
                  counts[17] = tot1[5] - counts[15] - counts[16];
		}
	      }
	      if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x, &freqx1, &freqx2, &freq11)) {
		freq11_expected = freqx1 * freq1x;
		rsq = freq11 - freq11_expected;
		rsq = rsq * rsq / (freq11_expected * freq2x * freqx2);
		if (rsq > prune_ld_thresh) {
		  at_least_one_prune = 1;
		  cur_exclude_ct++;
		  // remove marker with lower MAF
		  if (get_maf(set_allele_freqs[live_indices[ulii]]) < (1 - SMALL_EPSILON) * get_maf(set_allele_freqs[live_indices[uljj]])) {
		    SET_BIT(live_indices[ulii], pruned_arr);
		  } else {
		    SET_BIT(live_indices[uljj], pruned_arr);
		    uljj++;
		    while (uljj < cur_window_size) {
		      if (!IS_SET(pruned_arr, live_indices[uljj])) {
			break;
		      }
		      uljj++;
		    }
		    if (uljj < cur_window_size) {
		      start_arr[ulii] = live_indices[uljj];
		    }
		  }
		  break;
		}
	      }
	    }
	    if (uljj == cur_window_size) {
	    indep_pairphase_skip_marker:
	      start_arr[ulii] = window_unfiltered_end;
	    }
	  }
	} while (at_least_one_prune);
      }
      for (uii = 0; uii < ld_window_incr; uii++) {
	if (window_unfiltered_start == chrom_end) {
	  break;
	}
	window_unfiltered_start++;
	next_unset_ck(marker_exclude, chrom_end, &window_unfiltered_start);
      }
      if (window_unfiltered_start == chrom_end) {
	break;
      }
      if (window_unfiltered_start >= pct_thresh) {
	pct = ((window_unfiltered_start - chrom_start) * 100LLU) / (chrom_end - chrom_start);
	printf("\r%u%%", pct++);
	fflush(stdout);
	pct_thresh = chrom_start + (((uint64_t)pct * (chrom_end - chrom_start)) / 100);
      }
      uljj = 0;
      if (window_unfiltered_end < window_unfiltered_start) {
	window_unfiltered_end = window_unfiltered_start;
      }
      // copy back previously loaded/computed results
      while (live_indices[uljj] < window_unfiltered_start) {
	uljj++;
	if (uljj == cur_window_size) {
	  break;
	}
      }
      for (ulii = 0; uljj < cur_window_size; uljj++) {
	if (IS_SET(pruned_arr, live_indices[uljj])) {
	  continue;
	}
	memcpy(&(geno[ulii * founder_ctsplit]), &(geno[uljj * founder_ctsplit]), founder_ctsplit * sizeof(intptr_t));
	live_indices[ulii] = live_indices[uljj];
	start_arr[ulii] = start_arr[uljj];
	memcpy(&(cur_tots[ulii * 3]), &(cur_tots[uljj * 3]), 3 * sizeof(int32_t));
	// bugfix: forgot to update zmiss
	if (IS_SET(zmiss, uljj)) {
	  SET_BIT(ulii, zmiss);
	} else {
	  CLEAR_BIT(ulii, zmiss);
	}
	ulii++;
      }
      clear_bits(ulii, window_max, zmiss);

      prev_end = ulii;
      cur_window_size = ulii;
      if (window_is_kb) {
	uljj = 0;
	ulkk = window_unfiltered_end;
	while ((window_unfiltered_end < chrom_end) && (marker_pos[window_unfiltered_end] <= marker_pos[window_unfiltered_start] + (1000 * ld_window_size))) {
	  uljj++;
	  window_unfiltered_end++;
	  next_unset_ck(marker_exclude, chrom_end, &window_unfiltered_end);
	}
	window_unfiltered_end = ulkk;
      } else {
	uljj = ld_window_incr;
      }
      for (ulii = 0; ulii < uljj; ulii++) {
	if (window_unfiltered_end == chrom_end) {
	  break;
	}
	live_indices[cur_window_size] = window_unfiltered_end;
	if (cur_window_size > prev_end) {
	  start_arr[cur_window_size - 1] = window_unfiltered_end;
	}
	if (fseeko(bedfile, bed_offset + (window_unfiltered_end * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	  goto indep_pairphase_ret_READ_FAIL;
	}
	if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, IS_SET(marker_reverse, window_unfiltered_end), bedfile, loadbuf_raw, loadbuf)) {
	  goto indep_pairphase_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)loadbuf);
	}
	cur_geno1 = &(geno[cur_window_size * founder_ctsplit]);
	load_and_split3(nullptr, loadbuf, founder_ct, cur_geno1, dummy_nm, dummy_nm, founder_ctv3, 0, 0, 1, &ulkk);
	cur_tots[((uintptr_t)cur_window_size) * 3] = popcount_longs(cur_geno1, founder_ctv3);
	cur_tots[((uintptr_t)cur_window_size) * 3 + 1] = popcount_longs(&(cur_geno1[founder_ctv3]), founder_ctv3);
	cur_tots[((uintptr_t)cur_window_size) * 3 + 2] = popcount_longs(&(cur_geno1[2 * founder_ctv3]), founder_ctv3);
	if ((!cur_tots[((uintptr_t)cur_window_size) * 3 + 1]) && ((!cur_tots[((uintptr_t)cur_window_size) * 3]) || (!cur_tots[((uintptr_t)cur_window_size) * 3 + 2]))) {
	  SET_BIT(window_unfiltered_end, pruned_arr);
	  cur_exclude_ct++;
	} else if (ulkk == 3) {
	  SET_BIT(cur_window_size, zmiss);
	}
	cur_window_size++;
	window_unfiltered_end++;
	next_unset_ck(marker_exclude, chrom_end, &window_unfiltered_end);
      }
      if (cur_window_size > prev_end) {
	start_arr[cur_window_size - 1] = window_unfiltered_end;
      }
    }
    putc_unlocked('\r', stdout);
    LOGPRINTF("Pruned %" PRIuPTR " variant%s from chromosome %u, leaving %" PRIuPTR ".\n", cur_exclude_ct, (cur_exclude_ct == 1)? "" : "s", cur_chrom, chrom_end - chrom_start - popcount_bit_idx(marker_exclude, chrom_start, chrom_end) - cur_exclude_ct);
    tot_exclude_ct += cur_exclude_ct;

    // advance chromosomes as necessary
    window_unfiltered_start = ld_prune_next_valid_chrom_start(pruned_arr, window_unfiltered_start, chrom_info_ptr, chrom_code_end, unfiltered_marker_ct);
  } while (window_unfiltered_start < unfiltered_marker_ct);

  LOGPRINTF("Pruning complete.  %u of %" PRIuPTR " variants removed.\n", tot_exclude_ct, marker_ct);
  retval = ld_prune_write(outname, outname_end, marker_exclude, pruned_arr, marker_ids, max_marker_id_len, chrom_info_ptr, chrom_code_end);
  if (retval) {
    goto indep_pairphase_ret_1;
  }

  while (0) {
  indep_pairphase_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  indep_pairphase_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  indep_pairphase_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 indep_pairphase_ret_1:
  bigstack_reset(bigstack_mark);
  return retval;
}

typedef struct ll_epi_summary_struct {
  struct ll_epi_summary_struct* next;
  double best_chisq;
  char* best_chr_and_snp; // separate allocation; tab-delimited
  uint32_t n_sig;
  uint32_t n_tot;
  uint32_t id_len; // variant ID NOT null-terminated
  char strbuf[];
} Ll_epi_summary;

// N.B. moves g_bigstack_base in word-size instead of cacheline increments
Ll_epi_summary* lle_alloc(char* chrom_id, uint32_t chrom_len, char* marker_id, uint32_t marker_id_len, uint32_t nsig, uint32_t ntot, double chisq) {
  uintptr_t alloc_size = (sizeof(Ll_epi_summary) + chrom_len + marker_id_len + sizeof(intptr_t)) & (~(sizeof(intptr_t) - ONELU));
  Ll_epi_summary* newptr = (Ll_epi_summary*)g_bigstack_base;
  if (bigstack_left() < alloc_size) {
    return nullptr;
  }
  g_bigstack_base = &(g_bigstack_base[alloc_size]);
  newptr->next = nullptr;
  newptr->best_chisq = chisq;
  newptr->n_sig = nsig;
  newptr->n_tot = ntot;
  newptr->id_len = marker_id_len;
  memcpy(newptr->strbuf, marker_id, marker_id_len);
  memcpyx(&(newptr->strbuf[marker_id_len]), chrom_id, chrom_len, '\0');
  return newptr;
}

int32_t validate_epistasis_summary_header(char* bufptr) {
  uint32_t slen = strlen_se(bufptr);
  int32_t retval = 0;
  if ((slen != 3) || memcmp(bufptr, "CHR", 3)) {
    return RET_INVALID_FORMAT;
  }
  bufptr = skip_initial_spaces(&(bufptr[3]));
  slen = strlen_se(bufptr);
  if ((slen != 3) || memcmp(bufptr, "SNP", 3)) {
    return RET_INVALID_FORMAT;
  }
  bufptr = skip_initial_spaces(&(bufptr[3]));
  slen = strlen_se(bufptr);
  if ((slen != 5) || memcmp(bufptr, "N_SIG", 5)) {
    return RET_INVALID_FORMAT;
  }
  bufptr = skip_initial_spaces(&(bufptr[5]));
  slen = strlen_se(bufptr);
  if ((slen != 5) || memcmp(bufptr, "N_TOT", 5)) {
    return RET_INVALID_FORMAT;
  }
  bufptr = skip_initial_spaces(&(bufptr[5]));
  slen = strlen_se(bufptr);
  if (slen == 4) {
    if (memcmp(bufptr, "PROP", 4)) {
      return RET_INVALID_FORMAT;
    }
    retval = -1;
    bufptr = skip_initial_spaces(&(bufptr[4]));
    slen = strlen_se(bufptr);
  }
  if ((slen != 10) || memcmp(bufptr, "BEST_CHISQ", 10)) {
    return RET_INVALID_FORMAT;
  }
  bufptr = skip_initial_spaces(&(bufptr[10]));
  slen = strlen_se(bufptr);
  if ((slen != 8) || memcmp(bufptr, "BEST_CHR", 8)) {
    return RET_INVALID_FORMAT;
  }
  bufptr = skip_initial_spaces(&(bufptr[8]));
  slen = strlen_se(bufptr);
  if ((slen != 8) || memcmp(bufptr, "BEST_SNP", 8)) {
    return RET_INVALID_FORMAT;
  }
  bufptr = skip_initial_spaces(&(bufptr[8]));
  if (!is_eoln_kns(*bufptr)) {
    return RET_INVALID_FORMAT;
  }
  return retval;
}

int32_t epi_summary_merge(Epi_info* epi_ip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  FILE* outfile = nullptr;
  char* inprefix = epi_ip->summary_merge_prefix;
  char* inprefix_end = (char*)memchr(inprefix, 0, FNAMESIZE);
  Ll_epi_summary* list_start = nullptr;
  // first .3 entry is later than first .2 entry, etc., so we can save
  // ourselves some linked list traversal time by starting the first-entry scan
  // after where the last one left off.
  Ll_epi_summary* last_start = nullptr;
  Ll_epi_summary** lle_pp = &list_start; // end-of-list pointer for first file
  uint32_t file_ct = epi_ip->summary_merge_ct;
  int32_t retval = 0;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* id_ptr;
  char* nsig_ptr;
  char* ntot_ptr;
  char* best_chisq_ptr;
  char* best_chr_ptr;
  char* best_marker_ptr;
  Ll_epi_summary* lle_ptr; // traverser for remaining files
  uintptr_t line_idx;
  uintptr_t ulii;
  double cur_chisq;
  uint32_t plink_maxsnp;
  uint32_t file_idx;
  uint32_t chrom_len;
  uint32_t id_len;
  uint32_t is_first_entry;
  int32_t nsig;
  int32_t ntot;
  if (inprefix_end[-1] == '.') {
    inprefix_end--;
  }
  ulii = (uintptr_t)(inprefix_end - inprefix);
  if ((ulii >= 16) && (!memcmp(".summary", &(inprefix[ulii - 8]), 8))) {
    inprefix_end -= 8;
    ulii -= 8;
  }
  bufptr = &(inprefix[ulii - 2]);
  if (memcmp(".epi.", &(inprefix[ulii - 7]), 5) || (memcmp("cc", bufptr, 2) && memcmp("co", bufptr, 2) && memcmp("qt", bufptr, 2))) {
    LOGERRPRINTFWW("Error: Invalid --epistasis-summary-merge filename prefix '%s'. (*.epi.cc, *.epi.co, or *.epi.qt expected.)\n", inprefix);
    goto epi_summary_merge_ret_INVALID_CMDLINE;
  }
  inprefix_end = memcpya(inprefix_end, ".summary.", 9);
  memcpyx(outname_end, &(inprefix[ulii - 7]), 15, '\0');
  // Started out using a hash table, but on second thought, it's unnecessary
  // given the possibilities for distributed .summary files.
  // 1. ALL x ALL, SET x SET: First file lists all marker IDs in the final
  //    order; first entry in remaining files should match an entry in the
  //    middle and the rest match sequentially from there.
  // 2. SET1 x ALL, SET1 x SET2: No duplication whatsoever.  Output will be
  //    such that cat actually works, but asking users to conditionally use cat
  //    would add confusion for little reason; instead we detect the telltale
  //    "PROP" in the first file's header line and switch to cat.

  g_textbuf[MAXLINELEN - 1] = ' ';
  memcpy(inprefix_end, "1", 2);
  if (fopen_checked(inprefix, "r", &infile)) {
    goto epi_summary_merge_ret_OPEN_FAIL;
  }
  retval = load_to_first_token(infile, MAXLINELEN, '\0', "--epistasis-summary-merge file", g_textbuf, &bufptr, &line_idx);
  if (retval) {
    goto epi_summary_merge_ret_1;
  }
  retval = validate_epistasis_summary_header(bufptr);
  if (retval) {
    if (retval == -1) {
      // switch to cat mode.  meow.
      fclose_null(&infile);
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto epi_summary_merge_ret_OPEN_FAIL;
      }
      for (file_idx = 1; file_idx <= file_ct; file_idx++) {
        uint32toa_x(file_idx, '\0', inprefix_end);
	if (fopen_checked(inprefix, FOPEN_RB, &infile)) {
	  goto epi_summary_merge_ret_OPEN_FAIL;
	}
	while (1) {
	  ulii = fread(g_textbuf, 1, MAXLINELEN, infile);
          if (!ulii) {
	    break;
	  }
	  if (fwrite_checked(g_textbuf, ulii, outfile)) {
	    goto epi_summary_merge_ret_WRITE_FAIL;
	  }
	}
	if (fclose_null(&infile)) {
	  goto epi_summary_merge_ret_READ_FAIL;
	}
      }
      retval = 0;
      goto epi_summary_merge_success;
    }
    goto epi_summary_merge_ret_INVALID_HEADER;
  }
  bufptr2 = token_endnn(bufptr);
  bufptr = skip_initial_spaces(bufptr2);
  plink_maxsnp = ((uintptr_t)(token_endnn(bufptr) - bufptr2)) - 1;
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      goto epi_summary_merge_ret_LONG_LINE;
    }
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    chrom_len = strlen_se(bufptr);
    id_ptr = skip_initial_spaces(&(bufptr[chrom_len]));
    id_len = strlen_se(id_ptr);
    nsig_ptr = skip_initial_spaces(&(id_ptr[id_len]));
    ntot_ptr = next_token(nsig_ptr);
    best_chisq_ptr = next_token(ntot_ptr);
    best_chr_ptr = next_token(best_chisq_ptr);
    if (no_more_tokens_kns(best_chr_ptr)) {
      goto epi_summary_merge_ret_MISSING_TOKENS;
    }
    if (scan_uint_icap(nsig_ptr, (uint32_t*)&nsig)) {
      goto epi_summary_merge_ret_INVALID_NSIG;
    }
    if (scan_uint_icap(ntot_ptr, (uint32_t*)&ntot)) {
      goto epi_summary_merge_ret_INVALID_NTOT;
    }
    if (ntot) {
      if (scan_double(best_chisq_ptr, &cur_chisq)) {
	goto epi_summary_merge_ret_INVALID_CHISQ;
      }
    } else {
      cur_chisq = 0;
    }
    *lle_pp = lle_alloc(bufptr, chrom_len, id_ptr, id_len, nsig, ntot, cur_chisq);
    if (!(*lle_pp)) {
      goto epi_summary_merge_ret_NOMEM;
    }
    chrom_len = strlen_se(best_chr_ptr);
    best_marker_ptr = skip_initial_spaces(&(best_chr_ptr[chrom_len]));
    id_len = strlen_se(best_marker_ptr);
    if (!id_len) {
      goto epi_summary_merge_ret_MISSING_TOKENS;
    }
    // throw in an extra word, to reduce the need for reallocation
    ulii = (chrom_len + id_len + 1 + 2 * sizeof(intptr_t)) & (~(sizeof(intptr_t) - 1));
    if (ulii > bigstack_left()) {
      goto epi_summary_merge_ret_NOMEM;
    }
    bufptr = (char*)g_bigstack_base;
    memcpyx(bufptr, best_chr_ptr, chrom_len, '\t');
    memcpy(&(bufptr[chrom_len + 1]), best_marker_ptr, id_len);
    // pad with nulls then tab-terminate, so we can find the buffer end later
    memset(&(bufptr[chrom_len + id_len + 1]), 0, ulii - chrom_len - id_len - 2);
    bufptr[ulii - 1] = '\t';
    (*lle_pp)->best_chr_and_snp = bufptr;
    lle_pp = &((*lle_pp)->next);
    g_bigstack_base = &(g_bigstack_base[ulii]);
  }
  if (fclose_null(&infile)) {
    goto epi_summary_merge_ret_READ_FAIL;
  }
  if (!list_start) {
    LOGPREPRINTFWW("Error: %s has no entries.\n", inprefix);
    goto epi_summary_merge_ret_INVALID_FORMAT_2;
  }
  last_start = list_start->next;
  for (file_idx = 2; file_idx <= file_ct; file_idx++) {
    uint32toa_x(file_idx, '\0', inprefix_end);
    if (fopen_checked(inprefix, "r", &infile)) {
      goto epi_summary_merge_ret_OPEN_FAIL;
    }
    retval = load_to_first_token(infile, MAXLINELEN, '\0', "--epistasis-summary-merge file", g_textbuf, &bufptr, &line_idx);
    if (retval) {
      goto epi_summary_merge_ret_1;
    }
    retval = validate_epistasis_summary_header(bufptr);
    if (retval) {
      goto epi_summary_merge_ret_INVALID_HEADER;
    }
    lle_ptr = last_start;
    is_first_entry = 1;
    while (fgets(g_textbuf, MAXLINELEN, infile)) {
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 1]) {
	goto epi_summary_merge_ret_LONG_LINE;
      }
      bufptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      if (!lle_ptr) {
        LOGPREPRINTFWW("Error: More lines than expected in %s.\n", inprefix);
	goto epi_summary_merge_ret_INVALID_FORMAT_2;
      }
      chrom_len = strlen_se(bufptr);
      id_ptr = skip_initial_spaces(&(bufptr[chrom_len]));
      id_len = strlen_se(id_ptr);
      nsig_ptr = skip_initial_spaces(&(id_ptr[id_len]));
      ntot_ptr = next_token(nsig_ptr);
      best_chisq_ptr = next_token(ntot_ptr);
      best_chr_ptr = next_token(best_chisq_ptr);
      if (no_more_tokens_kns(best_chr_ptr)) {
	goto epi_summary_merge_ret_MISSING_TOKENS;
      }
      if (scan_uint_icap(nsig_ptr, (uint32_t*)&nsig)) {
	goto epi_summary_merge_ret_INVALID_NSIG;
      }
      if (scan_uint_icap(ntot_ptr, (uint32_t*)&ntot)) {
	goto epi_summary_merge_ret_INVALID_NTOT;
      }
      if (!is_first_entry) {
	if ((lle_ptr->id_len != id_len) || memcmp(lle_ptr->strbuf, id_ptr, id_len) || (strlen(&(lle_ptr->strbuf[id_len])) != chrom_len) || memcmp(&(lle_ptr->strbuf[id_len]), bufptr, chrom_len)) {
	  goto epi_summary_merge_ret_MISMATCH;
	}
      } else {
	while (1) {
	  if (!lle_ptr) {
	    goto epi_summary_merge_ret_MISMATCH;
	  }
	  if ((lle_ptr->id_len == id_len) && (!memcmp(lle_ptr->strbuf, id_ptr, id_len))) {
	    break;
	  }
          lle_ptr = lle_ptr->next;
	}
        if ((strlen(&(lle_ptr->strbuf[id_len])) != chrom_len) || memcmp(&(lle_ptr->strbuf[id_len]), bufptr, chrom_len)) {
	  goto epi_summary_merge_ret_MISMATCH;
	}
	last_start = lle_ptr->next;
	is_first_entry = 0;
      }
      if (ntot) {
	if (scan_double(best_chisq_ptr, &cur_chisq)) {
	  goto epi_summary_merge_ret_INVALID_CHISQ;
	}
	lle_ptr->n_sig += nsig;
	lle_ptr->n_tot += ntot;
        if (cur_chisq > lle_ptr->best_chisq) {
	  chrom_len = strlen_se(best_chr_ptr);
          best_marker_ptr = skip_initial_spaces(&(best_chr_ptr[chrom_len]));
          id_len = strlen_se(best_marker_ptr);
	  if (!id_len) {
	    goto epi_summary_merge_ret_MISSING_TOKENS;
	  }
          lle_ptr->best_chisq = cur_chisq;
	  bufptr = lle_ptr->best_chr_and_snp;
          bufptr2 = (char*)memchr(bufptr, '\t', MAXLINELEN);
	  bufptr3 = (char*)memchr(++bufptr2, 0, MAXLINELEN);
	  bufptr4 = (char*)memchr(bufptr3, '\t', MAXLINELEN);
	  ulii = (uintptr_t)(bufptr4 - bufptr);
	  if (ulii <= chrom_len + id_len + 1) {
	    ulii = (chrom_len + id_len + 1 + sizeof(intptr_t)) & (~(sizeof(intptr_t) - 1));
            if (ulii > bigstack_left()) {
	      goto epi_summary_merge_ret_NOMEM;
	    }
            bufptr = (char*)g_bigstack_base;
	    bufptr3 = &(bufptr[ulii - 1]);
	    *bufptr3 = '\t';
            lle_ptr->best_chr_and_snp = bufptr;
            g_bigstack_base = &(g_bigstack_base[ulii]);
	  }
	  bufptr = memcpyax(bufptr, best_chr_ptr, chrom_len, '\t');
	  bufptr = memcpya(bufptr, best_marker_ptr, id_len);
	  if (bufptr < bufptr3) {
	    memset(bufptr, 0, bufptr3 - bufptr);
	  }
	}
      }
      lle_ptr = lle_ptr->next;
    }
    if (fclose_null(&infile)) {
      goto epi_summary_merge_ret_READ_FAIL;
    }
  }

  if (fopen_checked(outname, "w", &outfile)) {
    goto epi_summary_merge_ret_OPEN_FAIL;
  }
  bufptr = memcpya(g_textbuf, " CHR ", 5);
  bufptr = fw_strcpyn(plink_maxsnp, 3, "SNP", bufptr);
  bufptr = strcpya(bufptr, "        N_SIG        N_TOT         PROP   BEST_CHISQ BEST_CHR ");
  bufptr = fw_strcpyn(plink_maxsnp, 8, "BEST_SNP", bufptr);
  bufptr = memcpya(bufptr, " \n", 2);
  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
    goto epi_summary_merge_ret_WRITE_FAIL;
  }
  lle_ptr = list_start;
  do {
    bufptr2 = lle_ptr->strbuf;
    id_len = lle_ptr->id_len;
    bufptr3 = &(bufptr2[id_len]);
    bufptr = fw_strcpy(4, bufptr3, g_textbuf);
    *bufptr++ = ' ';
    bufptr = fw_strcpyn(plink_maxsnp, id_len, bufptr2, bufptr);
    nsig = lle_ptr->n_sig;
    ntot = lle_ptr->n_tot;
    bufptr = memseta(bufptr, 32, 3);
    bufptr = uint32toa_w10(nsig, bufptr);
    bufptr = memseta(bufptr, 32, 3);
    bufptr = uint32toa_w10x(ntot, ' ', bufptr);
    bufptr = dtoa_g_wxp4x(((double)((int32_t)nsig)) / ((double)((int32_t)ntot)), 12, ' ', bufptr);
    bufptr = dtoa_g_wxp4x(lle_ptr->best_chisq, 12, ' ', bufptr);
    // no need to special-case ntot == 0, this code correctly copies 'NA'
    bufptr2 = lle_ptr->best_chr_and_snp;
    bufptr3 = (char*)memchr(bufptr2, '\t', MAXLINELEN);
    ulii = (uintptr_t)(bufptr3 - bufptr2);
    bufptr = fw_strcpyn(4, ulii, bufptr2, bufptr);
    *bufptr++ = ' ';
    bufptr = fw_strcpy(plink_maxsnp, &(bufptr3[1]), bufptr);
    bufptr = memcpya(bufptr, " \n", 2);
    if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
      goto epi_summary_merge_ret_WRITE_FAIL;
    }
    lle_ptr = lle_ptr->next;
  } while (lle_ptr);

 epi_summary_merge_success:
  if (fclose_null(&outfile)) {
    // just kidding!  no success
    goto epi_summary_merge_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--epistasis-summary-merge: Merged summary written to %s .\n", outname);
  while (0) {
  epi_summary_merge_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  epi_summary_merge_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  epi_summary_merge_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  epi_summary_merge_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  epi_summary_merge_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  epi_summary_merge_ret_MISMATCH:
    logerrprint("Error: --epistasis-summary-merge files were generated from different datasets\nand/or settings.\n");
    retval = RET_INVALID_FORMAT;
    break;
  epi_summary_merge_ret_INVALID_NSIG:
    LOGERRPRINTFWW("Error: Invalid N_SIG value on line %" PRIuPTR " of %s .\n", line_idx, inprefix);
    retval = RET_INVALID_FORMAT;
    break;
  epi_summary_merge_ret_INVALID_NTOT:
    LOGERRPRINTFWW("Error: Invalid N_SIG value on line %" PRIuPTR " of %s .\n", line_idx, inprefix);
    retval = RET_INVALID_FORMAT;
    break;
  epi_summary_merge_ret_INVALID_CHISQ:
    LOGERRPRINTFWW("Error: Invalid BEST_CHISQ value on line %" PRIuPTR " of %s .\n", line_idx, inprefix);
    retval = RET_INVALID_FORMAT;
    break;
  epi_summary_merge_ret_MISSING_TOKENS:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, inprefix);
    retval = RET_INVALID_FORMAT;
    break;
  epi_summary_merge_ret_LONG_LINE:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, inprefix);
    retval = RET_INVALID_FORMAT;
    break;
  epi_summary_merge_ret_INVALID_HEADER:
    LOGPREPRINTFWW(g_logbuf, "Error: Invalid --epistasis-summary-merge header in %s.\n", inprefix);
  epi_summary_merge_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 epi_summary_merge_ret_1:
  fclose_cond(infile);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

void test_mishap_write_line(FILE* outfile, char* wptr, uint32_t prev_alen, uint32_t next_alen, const char* prev_aptr, const char* next_aptr, double* total_cts, double* curhap_cts, double tot_recip, double output_min_p, char* flankstr, uint32_t flanklen) {
  // total_cts[0] = caseN[0] + caseN[1]
  // total_cts[1] = controlN[0] + controlN[1]
  char* tbuf_cur = g_textbuf;
  double casen_1 = total_cts[0] - curhap_cts[0];
  double ctrln_1 = total_cts[1] - curhap_cts[1];
  uint32_t uii = prev_alen + next_alen;
  char* wptr2;
  double row_mult;
  double cur_expected;
  double dxx;
  double chisq;
  if (uii <= 10) {
    wptr = memseta(wptr, 32, 10 - uii);
    if (prev_alen) {
      wptr = memcpya(wptr, prev_aptr, prev_alen);
    }
    if (next_alen) {
      wptr = memcpya(wptr, next_aptr, next_alen);
    }
  } else {
    fwrite(g_textbuf, 1, (uintptr_t)(wptr - g_textbuf), outfile);
    if (prev_alen) {
      fputs(prev_aptr, outfile);
    }
    if (next_alen) {
      fputs(next_aptr, outfile);
    }
    tbuf_cur = wptr;
  }
  *wptr++ = ' ';
  if (total_cts[0] > 0.0) {
    wptr = dtoa_g_wxp3(curhap_cts[0] / total_cts[0], 8, wptr);
  } else {
    wptr = memcpya(wptr, "      NA", 8);
  }
  *wptr++ = ' ';
  if (total_cts[1] > 0.0) {
    wptr = dtoa_g_wxp3(curhap_cts[1] / total_cts[1], 8, wptr);
  } else {
    wptr = memcpya(wptr, "      NA", 8);
  }
  *wptr++ = ' ';
  wptr2 = dtoa_g(curhap_cts[0], wptr);
  *wptr2++ = '/';
  wptr2 = dtoa_g(curhap_cts[1], wptr2);
  wptr = width_force(20, wptr, wptr2);
  *wptr++ = ' ';
  wptr2 = dtoa_g(casen_1, wptr);
  *wptr2++ = '/';
  wptr2 = dtoa_g(ctrln_1, wptr2);
  wptr = width_force(20, wptr, wptr2);
  *wptr++ = ' ';
  if ((curhap_cts[0] > 0.0) && (curhap_cts[1] > 0.0) && (casen_1 > 0.0) && (ctrln_1 > 0.0)) {
    row_mult = (curhap_cts[0] + curhap_cts[1]) * tot_recip;
    cur_expected = row_mult * total_cts[0];
    dxx = curhap_cts[0] - cur_expected;
    chisq = dxx * dxx / cur_expected;
    cur_expected = row_mult * total_cts[1];
    dxx = curhap_cts[1] - cur_expected;
    chisq += dxx * dxx / cur_expected;
    row_mult = (total_cts[0] + total_cts[1]) * tot_recip - row_mult;
    cur_expected = row_mult * total_cts[0];
    dxx = casen_1 - cur_expected;
    chisq += dxx * dxx / cur_expected;
    cur_expected = row_mult * total_cts[1];
    dxx = ctrln_1 - cur_expected;
    chisq += dxx * dxx / cur_expected;
    wptr = dtoa_g_wxp3(chisq, 8, wptr);
    *wptr++ = ' ';
    dxx = chiprob_p(chisq, 1);
    wptr = dtoa_g_wxp3(MAXV(dxx, output_min_p), 8, wptr);
  } else {
    wptr = memcpya(wptr, "      NA       NA", 17);
  }
  wptr = memcpya(wptr, flankstr, flanklen);
  fwrite(tbuf_cur, 1, (uintptr_t)(wptr - tbuf_cur), outfile);
}

int32_t test_mishap(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, double min_maf, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  uintptr_t sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  uintptr_t final_mask = get_final_mask(sample_ct);
  char* tbuf2 = &(g_textbuf[MAXLINELEN]);
  char* wptr2 = nullptr;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t inspected_ct = 0;
  uint32_t missing_ct_next = 0;
  uint32_t prev_a1len = 0;
  uint32_t prev_a2len = 0;
  int32_t retval = 0;
  // need following counts:
  //   all 9 flanking hap combinations | current missing
  //     [0]: prev = 00, next = 00
  //     [1]: prev = 00, next = 10
  //     [2]: prev = 00, next = 11
  //     [3]: prev = 10, next = 00
  //     ...
  //   all 9 flanking hap combinations | current nonmissing [9..17]
  uint32_t counts[27];
  // [0]: central call missing, all haps clearing --maf (caseN[0] + caseN[1])
  // [1]: central call nonmissing, all haps clear --maf (ctrlN[0] + ctrlN[1])
  // [2k], k in 1..4: caseN[0] for current hap
  // [2k+1]: controlN[0] for current hap
  // note that all numbers are actually double raw counts
  double hap_ct_table[10];
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_end;
  uintptr_t* prevsnp_ptr;
  uintptr_t* cursnp_ptr;
  uintptr_t* nextsnp_ptr;
  uintptr_t* maskbuf_mid;
  uintptr_t* maskbuf;
  char* wptr;
  uint32_t* uiptr;
  uintptr_t marker_uidx_prev;
  uintptr_t marker_uidx_cur;
  uintptr_t marker_uidx_next;
  double freq11;
  double tot_recip;
  double dxx;
  double dyy;
  double dzz;
  double dww;
  double orig_cmiss_tot;
  double orig_cnm_tot;
  uint32_t flanklen;
  uint32_t missing_ct_cur;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t next_a1len;
  uint32_t next_a2len;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  if (is_set(chrom_info_ptr->haploid_mask, 1)) {
    logerrprint("Error: --test-mishap can only be used on diploid genomes.\n");
    goto test_mishap_ret_INVALID_CMDLINE;
  }
  if (sample_ct >= 0x40000000) {
    logerrprint("Error: --test-mishap does not support >= 2^30 samples.\n");
    goto test_mishap_ret_INVALID_CMDLINE;
  }
  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(sample_ctv2 * 3, &loadbuf) ||
      bigstack_alloc_ul(sample_ctv2, &maskbuf_mid) ||
      bigstack_alloc_ul(sample_ctv2, &maskbuf)) {
    goto test_mishap_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  loadbuf[sample_ctv2 - 2] = 0;
  loadbuf[sample_ctv2 - 1] = 0;
  loadbuf[2 * sample_ctv2 - 2] = 0;
  loadbuf[2 * sample_ctv2 - 1] = 0;
  loadbuf[3 * sample_ctv2 - 2] = 0;
  loadbuf[3 * sample_ctv2 - 1] = 0;
  loadbuf_end = &(loadbuf[sample_ctv2 * 3]);
  tbuf2[0] = ' ';
  memcpy(outname_end, ".missing.hap", 13);
  if (fopen_checked(outname, "w", &outfile)) {
    goto test_mishap_ret_OPEN_FAIL;
  }
  sprintf(g_textbuf, "%%%us  HAPLOTYPE      F_0      F_1                 M_H1                 M_H2    CHISQ        P FLANKING\n", plink_maxsnp);
  fprintf(outfile, g_textbuf, "SNP");
  min_maf *= 1 - SMALL_EPSILON;
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    if (is_set(chrom_info_ptr->haploid_mask, chrom_idx)) {
      continue;
    }
    marker_uidx_cur = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    marker_uidx_cur = next_unset_ul(marker_exclude, marker_uidx_cur, chrom_end);
    if (marker_uidx_cur == chrom_end) {
      continue;
    }
    marker_uidx_next = next_unset_ul(marker_exclude, marker_uidx_cur + 1, chrom_end);
    if (marker_uidx_next == chrom_end) {
      continue;
    }
    prevsnp_ptr = loadbuf;
    fill_ulong_zero(sample_ctl2, prevsnp_ptr);
    cursnp_ptr = &(loadbuf[sample_ctv2]);
    if (fseeko(bedfile, bed_offset + marker_uidx_cur * ((uint64_t)unfiltered_sample_ct4), SEEK_SET)) {
      goto test_mishap_ret_READ_FAIL;
    }
    if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx_cur), bedfile, loadbuf_raw, cursnp_ptr)) {
      goto test_mishap_ret_READ_FAIL;
    }
    missing_ct_cur = count_01(cursnp_ptr, sample_ctl2);
    marker_uidx_prev = ~ZEROLU;
    for (; marker_uidx_cur < chrom_end; marker_uidx_prev = marker_uidx_cur, marker_uidx_cur = marker_uidx_next, prevsnp_ptr = cursnp_ptr, cursnp_ptr = nextsnp_ptr, missing_ct_cur = missing_ct_next, marker_uidx_next++) {
      nextsnp_ptr = &(cursnp_ptr[sample_ctv2]);
      if (nextsnp_ptr == loadbuf_end) {
	nextsnp_ptr = loadbuf;
      }
      if (marker_uidx_next < chrom_end) {
	if (IS_SET(marker_exclude, marker_uidx_next)) {
	  marker_uidx_next = next_unset_ul(marker_exclude, marker_uidx_next, chrom_end);
	  if (marker_uidx_next == chrom_end) {
	    goto test_mishap_last_chrom_snp;
	  }
	  if (fseeko(bedfile, bed_offset + marker_uidx_next * ((uint64_t)unfiltered_sample_ct4), SEEK_SET)) {
	    goto test_mishap_ret_READ_FAIL;
	  }
	}
        if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx_next), bedfile, loadbuf_raw, nextsnp_ptr)) {
          goto test_mishap_ret_READ_FAIL;
	}
        missing_ct_next = count_01(nextsnp_ptr, sample_ctl2);
      } else {
      test_mishap_last_chrom_snp:
        fill_ulong_zero(sample_ctl2, nextsnp_ptr);
      }
      if (missing_ct_cur < 5) {
	continue;
      }
      quatervec_copy_only_01(cursnp_ptr, unfiltered_sample_ct, maskbuf_mid);
      uiptr = counts;
      for (uii = 0; uii < 2; uii++) {
	if (uii) {
	  quatervec_01_invert(unfiltered_sample_ct, maskbuf_mid);
	}
        for (ujj = 0; ujj < 3; ujj++) {
          vec_datamask(unfiltered_sample_ct, ujj + (ujj + 1) / 2, prevsnp_ptr, maskbuf_mid, maskbuf);
	  ukk = popcount01_longs(maskbuf, sample_ctl2);
	  genovec_3freq(nextsnp_ptr, maskbuf, sample_ctl2, &umm, &(uiptr[1]), &(uiptr[2]));
	  uiptr[0] = ukk - umm - uiptr[1] - uiptr[2];
	  uiptr = &(uiptr[3]);
	}
      }
      wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx_cur * max_marker_id_len]), g_textbuf);
      *wptr++ = ' ';
      if (marker_uidx_prev != ~ZEROLU) {
	prev_a1len = strlen(marker_allele_ptrs[2 * marker_uidx_prev]);
	prev_a2len = strlen(marker_allele_ptrs[2 * marker_uidx_prev + 1]);
	wptr2 = strcpya(&(tbuf2[1]), &(marker_ids[marker_uidx_prev * max_marker_id_len]));
      }
      if (marker_uidx_next < chrom_end) {
	next_a1len = strlen(marker_allele_ptrs[2 * marker_uidx_next]);
	next_a2len = strlen(marker_allele_ptrs[2 * marker_uidx_next + 1]);
	if (marker_uidx_prev != ~ZEROLU) {
	  hap_ct_table[0] = (int32_t)(2 * (counts[0] + counts[1] + counts[2] + counts[3] + counts[4] + counts[5] + counts[6] + counts[7] + counts[8]));
	  hap_ct_table[1] = (int32_t)(2 * (counts[9] + counts[10] + counts[11] + counts[12] + counts[13] + counts[14] + counts[15] + counts[16] + counts[17]));
	  tot_recip = hap_ct_table[0] + hap_ct_table[1];
	  if (tot_recip == 0.0) {
	    // minor change: skip markers with zero observations.  (output
	    // wouldn't match PLINK 1.07 anyway, due to EM phasing differences)
	    continue;
	  }
	  orig_cmiss_tot = hap_ct_table[0];
	  orig_cnm_tot = hap_ct_table[1];
	  *wptr2++ = '|';
	  wptr2 = strcpya(wptr2, &(marker_ids[marker_uidx_next * max_marker_id_len]));
	  *wptr2++ = '\n';
	  flanklen = (uintptr_t)(wptr2 - tbuf2);
	  hap_ct_table[2] = (int32_t)(2 * counts[0] + counts[1] + counts[3]);
	  hap_ct_table[3] = (int32_t)(2 * counts[9] + counts[10] + counts[12]);
	  hap_ct_table[4] = (int32_t)(2 * counts[2] + counts[1] + counts[5]);
	  hap_ct_table[5] = (int32_t)(2 * counts[11] + counts[10] + counts[14]);
	  hap_ct_table[6] = (int32_t)(2 * counts[6] + counts[3] + counts[7]);
	  hap_ct_table[7] = (int32_t)(2 * counts[15] + counts[12] + counts[16]);
	  hap_ct_table[8] = (int32_t)(2 * counts[8] + counts[5] + counts[7]);
	  hap_ct_table[9] = (int32_t)(2 * counts[17] + counts[14] + counts[16]);
	  if (counts[4] + counts[13]) {
	    for (uii = 0; uii < 9; uii++) {
	      counts[18 + uii] = counts[uii] + counts[9 + uii];
	    }
	    // no need to check return value
	    em_phase_hethet_nobase(&(counts[18]), 0, 0, &dxx, &dyy, &dzz, &dww, &freq11);
	    // share of counts[4]/counts[13] which goes to 11 or 22 haplotype
	    // (0.5 - dxx) is share which goes to 12/21 haps
	    // (conveniently, there's a 0.5 and a 2 which cancel out here)
	    dxx = (freq11 * tot_recip - (hap_ct_table[2] + hap_ct_table[3])) / ((double)((int32_t)(counts[4] + counts[13])));
	    dyy = ((int32_t)counts[4]) * dxx;
	    dzz = ((int32_t)counts[13]) * dxx;
	    hap_ct_table[2] += dyy;
	    hap_ct_table[3] += dzz;
	    hap_ct_table[8] += dyy;
	    hap_ct_table[9] += dzz;
	    dxx = 1.0 - dxx;
	    dyy = ((int32_t)counts[4]) * dxx;
	    dzz = ((int32_t)counts[13]) * dxx;
	    hap_ct_table[4] += dyy;
	    hap_ct_table[5] += dzz;
	    hap_ct_table[6] += dyy;
	    hap_ct_table[7] += dzz;
	  }
	  dxx = min_maf * tot_recip;
	  if (hap_ct_table[2] + hap_ct_table[3] < dxx) {
	    hap_ct_table[0] -= hap_ct_table[2];
	    hap_ct_table[1] -= hap_ct_table[3];
	    tot_recip -= hap_ct_table[2] + hap_ct_table[3];
	  }
	  if (hap_ct_table[4] + hap_ct_table[5] < dxx) {
	    hap_ct_table[0] -= hap_ct_table[4];
	    hap_ct_table[1] -= hap_ct_table[5];
	    tot_recip -= hap_ct_table[4] + hap_ct_table[5];
	  }
	  if (hap_ct_table[6] + hap_ct_table[7] < dxx) {
	    hap_ct_table[0] -= hap_ct_table[6];
	    hap_ct_table[1] -= hap_ct_table[7];
	    tot_recip -= hap_ct_table[6] + hap_ct_table[7];
	  }
	  if (hap_ct_table[8] + hap_ct_table[9] < dxx) {
	    hap_ct_table[0] -= hap_ct_table[8];
	    hap_ct_table[1] -= hap_ct_table[9];
	    tot_recip -= hap_ct_table[8] + hap_ct_table[9];
	  }
	  tot_recip = 1.0 / tot_recip;
	  if (hap_ct_table[2] + hap_ct_table[3] >= dxx) {
	    test_mishap_write_line(outfile, wptr, prev_a1len, next_a1len, marker_allele_ptrs[2 * marker_uidx_prev], marker_allele_ptrs[2 * marker_uidx_next], hap_ct_table, &(hap_ct_table[2]), tot_recip, output_min_p, tbuf2, flanklen);
	  }
	  if (hap_ct_table[6] + hap_ct_table[7] >= dxx) {
	    test_mishap_write_line(outfile, wptr, prev_a2len, next_a1len, marker_allele_ptrs[2 * marker_uidx_prev + 1], marker_allele_ptrs[2 * marker_uidx_next], hap_ct_table, &(hap_ct_table[6]), tot_recip, output_min_p, tbuf2, flanklen);
	  }
	  if (hap_ct_table[4] + hap_ct_table[5] >= dxx) {
	    test_mishap_write_line(outfile, wptr, prev_a1len, next_a2len, marker_allele_ptrs[2 * marker_uidx_prev], marker_allele_ptrs[2 * marker_uidx_next + 1], hap_ct_table, &(hap_ct_table[4]), tot_recip, output_min_p, tbuf2, flanklen);
	  }
	  if (hap_ct_table[8] + hap_ct_table[9] >= dxx) {
	    test_mishap_write_line(outfile, wptr, prev_a2len, next_a2len, marker_allele_ptrs[2 * marker_uidx_prev + 1], marker_allele_ptrs[2 * marker_uidx_next + 1], hap_ct_table, &(hap_ct_table[8]), tot_recip, output_min_p, tbuf2, flanklen);
	  }
	} else {
	  hap_ct_table[0] = (int32_t)(2 * (counts[0] + counts[1] + counts[2]));
	  hap_ct_table[1] = (int32_t)(2 * (counts[9] + counts[10] + counts[11]));
	  tot_recip = hap_ct_table[0] + hap_ct_table[1];
	  if (tot_recip == 0.0) {
	    continue;
	  }
	  orig_cmiss_tot = hap_ct_table[0];
	  orig_cnm_tot = hap_ct_table[1];
	  wptr2 = strcpya(&(tbuf2[1]), &(marker_ids[marker_uidx_next * max_marker_id_len]));
	  *wptr2++ = '\n';
	  flanklen = (uintptr_t)(wptr2 - tbuf2);
	  dxx = min_maf * tot_recip;
	  hap_ct_table[2] = (int32_t)(counts[0] * 2 + counts[1]);
	  hap_ct_table[3] = (int32_t)(counts[9] * 2 + counts[10]);
	  hap_ct_table[4] = (int32_t)(counts[2] * 2 + counts[1]);
	  hap_ct_table[5] = (int32_t)(counts[11] * 2 + counts[10]);
	  if (hap_ct_table[4] + hap_ct_table[5] < dxx) {
	    hap_ct_table[0] = hap_ct_table[2];
	    hap_ct_table[1] = hap_ct_table[3];
	    tot_recip = hap_ct_table[2] + hap_ct_table[3];
	  } else if (hap_ct_table[2] + hap_ct_table[3] < dxx) {
	    hap_ct_table[0] = hap_ct_table[4];
	    hap_ct_table[1] = hap_ct_table[5];
	    tot_recip = hap_ct_table[4] + hap_ct_table[5];
	  }
	  tot_recip = 1.0 / tot_recip;
	  if (hap_ct_table[2] + hap_ct_table[3] >= dxx) {
	    test_mishap_write_line(outfile, wptr, 0, next_a1len, nullptr, marker_allele_ptrs[2 * marker_uidx_next], hap_ct_table, &(hap_ct_table[2]), tot_recip, output_min_p, tbuf2, flanklen);
	  }
	  if (hap_ct_table[4] + hap_ct_table[5] >= dxx) {
	    test_mishap_write_line(outfile, wptr, 0, next_a2len, nullptr, marker_allele_ptrs[2 * marker_uidx_next + 1], hap_ct_table, &(hap_ct_table[4]), tot_recip, output_min_p, tbuf2, flanklen);
	  }
	}
      } else {
	hap_ct_table[0] = (int32_t)(2 * (counts[0] + counts[3] + counts[6]));
	hap_ct_table[1] = (int32_t)(2 * (counts[9] + counts[12] + counts[15]));
	tot_recip = hap_ct_table[0] + hap_ct_table[1];
	if (tot_recip == 0.0) {
	  continue;
	}
	orig_cmiss_tot = hap_ct_table[0];
	orig_cnm_tot = hap_ct_table[1];
	*wptr2++ = '\n';
	flanklen = (uintptr_t)(wptr2 - tbuf2);
	dxx = min_maf * tot_recip;
	hap_ct_table[2] = (int32_t)(counts[0] * 2 + counts[3]);
	hap_ct_table[3] = (int32_t)(counts[9] * 2 + counts[12]);
	hap_ct_table[4] = (int32_t)(counts[6] * 2 + counts[3]);
	hap_ct_table[5] = (int32_t)(counts[15] * 2 + counts[12]);
	if (hap_ct_table[4] + hap_ct_table[5] < dxx) {
	  hap_ct_table[0] = hap_ct_table[2];
	  hap_ct_table[1] = hap_ct_table[3];
	  tot_recip = hap_ct_table[2] + hap_ct_table[3];
	} else if (hap_ct_table[2] + hap_ct_table[3] < dxx) {
	  hap_ct_table[0] = hap_ct_table[4];
	  hap_ct_table[1] = hap_ct_table[5];
	  tot_recip = hap_ct_table[4] + hap_ct_table[5];
	}
	tot_recip = 1.0 / tot_recip;
	if (hap_ct_table[2] + hap_ct_table[3] >= dxx) {
	  test_mishap_write_line(outfile, wptr, prev_a1len, 0, marker_allele_ptrs[2 * marker_uidx_prev], nullptr, hap_ct_table, &(hap_ct_table[2]), tot_recip, output_min_p, tbuf2, flanklen);
	}
	if (hap_ct_table[4] + hap_ct_table[5] >= dxx) {
	  test_mishap_write_line(outfile, wptr, prev_a2len, 0, marker_allele_ptrs[2 * marker_uidx_prev + 1], nullptr, hap_ct_table, &(hap_ct_table[4]), tot_recip, output_min_p, tbuf2, flanklen);
	}
      }
      hap_ct_table[0] = orig_cmiss_tot * 0.5;
      hap_ct_table[1] = orig_cnm_tot * 0.5;
      hap_ct_table[2] = (int32_t)(counts[1] + counts[3] + counts[4] + counts[5] + counts[7]);
      hap_ct_table[3] = (int32_t)(counts[10] + counts[12] + counts[13] + counts[14] + counts[16]);
      test_mishap_write_line(outfile, wptr, 6, 0, "HETERO", nullptr, hap_ct_table, &(hap_ct_table[2]), 1.0 / (hap_ct_table[0] + hap_ct_table[1]), output_min_p, tbuf2, flanklen);
      inspected_ct++;
      if (!(inspected_ct % 1000)) {
        printf("\r--test-mishap: %uk loci checked.", inspected_ct / 1000);
        fflush(stdout);
      }
    }
  }

  if (fclose_null(&outfile)) {
    goto test_mishap_ret_WRITE_FAIL;
  }
  putc_unlocked('\r', stdout);
  if (inspected_ct < marker_ct) {
    LOGPRINTF("--test-mishap: %u loc%s checked (%" PRIuPTR " skipped).\n", inspected_ct, (inspected_ct == 1)? "us" : "i", marker_ct - inspected_ct);
    LOGPREPRINTFWW("Report written to %s .\n", outname);
  } else {
    LOGPREPRINTFWW("--test-mishap: %u loc%s checked, report written to %s .\n", inspected_ct, (inspected_ct == 1)? "us" : "i", outname);
  }
  logprintb();

  while (0) {
  test_mishap_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  test_mishap_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  test_mishap_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  test_mishap_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  test_mishap_ret_INVALID_CMDLINE:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

static uintptr_t* g_ld_load2_bitfield;
static uintptr_t* g_ld_result_bitfield;

THREAD_RET_TYPE ld_map_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t thread_ct = g_ld_thread_ct;
  // er, this use of "ctv" is nonstandard, probably want to fix this later
  uintptr_t marker_ctv = ((g_ld_marker_ct + 127) / 128) * (128 / BITCT);
  uintptr_t idx1_offset = g_ld_block_idx1;
  uintptr_t block_idx1_start = (tidx * g_ld_idx1_block_size) / thread_ct;
  uintptr_t block_idx1_end = ((tidx + 1) * g_ld_idx1_block_size) / thread_ct;
  uintptr_t founder_ct = g_ld_founder_ct;
  uintptr_t founder_ctwd = founder_ct / BITCT2;
  uintptr_t founder_ctwd12 = founder_ctwd / 12;
  uintptr_t founder_ctwd12_rem = founder_ctwd - (12 * founder_ctwd12);
  uintptr_t lshift_last = 2 * ((0x7fffffc0 - founder_ct) % BITCT2);
  uintptr_t founder_ct_192_long = g_ld_founder_ct_192_long;
  uintptr_t* geno1 = g_ld_geno1;
  uintptr_t* geno_masks1 = g_ld_geno_masks1;
  uint32_t* missing_cts1 = g_ld_missing_cts1;
  uint32_t founder_ct_mld_m1 = g_ld_founder_ct_mld_m1;
  uint32_t founder_ct_mld_rem = g_ld_founder_ct_mld_rem;
  uintptr_t* load2_bitfield = g_ld_load2_bitfield;
  uintptr_t* result_bitfield = g_ld_result_bitfield;
  double r2_thresh = g_ld_window_r2;
  int32_t dp_result[5];
  uintptr_t* geno_fixed_vec_ptr;
  uintptr_t* geno_var_vec_ptr;
  uintptr_t* mask_fixed_vec_ptr;
  uintptr_t* mask_var_vec_ptr;
  uintptr_t* geno2;
  uintptr_t* geno_masks2;
  uintptr_t* rb_cur;
  uint32_t* missing_cts2;
  uintptr_t block_idx1;
  uintptr_t block_idx2;
  double non_missing_ctd;
  double cov12;
  double dxx;
  double dyy;
  uint32_t marker_idx2_start;
  uint32_t marker_idx2;
  uint32_t marker_idx2_end;
  uint32_t fixed_missing_ct;
  uint32_t fixed_non_missing_ct;
  uint32_t non_missing_ct;
  uint32_t uii;
  while (1) {
    marker_idx2_start = g_ld_idx2_block_start;
    marker_idx2_end = g_ld_marker_ctm8;
    geno2 = g_ld_geno2;
    geno_masks2 = g_ld_geno_masks2;
    missing_cts2 = g_ld_missing_cts2;
    rb_cur = &(result_bitfield[block_idx1_start * marker_ctv]);
    for (block_idx1 = block_idx1_start; block_idx1 < block_idx1_end; block_idx1++, rb_cur = &(rb_cur[marker_ctv])) {
      marker_idx2 = block_idx1 + idx1_offset + 1;
      if (marker_idx2 < marker_idx2_start) {
	marker_idx2 = marker_idx2_start;
      } else if (marker_idx2 >= marker_idx2_end) {
        break;
      }
      marker_idx2 = next_set(rb_cur, marker_idx2, marker_idx2_end);
      if (marker_idx2 == marker_idx2_end) {
	continue;
      }
      fixed_missing_ct = missing_cts1[block_idx1];
      fixed_non_missing_ct = founder_ct - fixed_missing_ct;
      geno_fixed_vec_ptr = &(geno1[block_idx1 * founder_ct_192_long]);
      mask_fixed_vec_ptr = &(geno_masks1[block_idx1 * founder_ct_192_long]);
      block_idx2 = popcount_bit_idx(load2_bitfield, marker_idx2_start, marker_idx2);
      while (1) {
        geno_var_vec_ptr = &(geno2[block_idx2 * founder_ct_192_long]);
        mask_var_vec_ptr = &(geno_masks2[block_idx2 * founder_ct_192_long]);
        non_missing_ct = fixed_non_missing_ct - missing_cts2[block_idx2];
        if (fixed_missing_ct && missing_cts2[block_idx2]) {
          non_missing_ct += ld_missing_ct_intersect(mask_var_vec_ptr, mask_fixed_vec_ptr, founder_ctwd12, founder_ctwd12_rem, lshift_last);
	}
        dp_result[0] = founder_ct;
        dp_result[1] = -fixed_non_missing_ct;
        dp_result[2] = missing_cts2[block_idx2] - founder_ct;
        dp_result[3] = dp_result[1];
        dp_result[4] = dp_result[2];
	ld_dot_prod(geno_var_vec_ptr, geno_fixed_vec_ptr, mask_var_vec_ptr, mask_fixed_vec_ptr, dp_result, founder_ct_mld_m1, founder_ct_mld_rem);
	non_missing_ctd = (double)((int32_t)non_missing_ct);
        dxx = dp_result[1];
        dyy = dp_result[2];
        cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
        if (cov12 * cov12 <= r2_thresh * ((dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy))) {
          clear_bit(marker_idx2, rb_cur);
	}
	uii = marker_idx2++;
	if (is_set(rb_cur, marker_idx2)) {
	  if (marker_idx2 == marker_idx2_end) {
	    break;
	  }
	  block_idx2++;
	} else {
          marker_idx2 = next_set(rb_cur, marker_idx2, marker_idx2_end);
	  if (marker_idx2 == marker_idx2_end) {
	    break;
	  }
          block_idx2 += popcount_bit_idx(load2_bitfield, uii, marker_idx2);
	}
      }
    }
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

int32_t construct_ld_map(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t* marker_reverse, uint32_t* marker_idx_to_uidx, uintptr_t unfiltered_sample_ct, uintptr_t* founder_pnm, Set_info* sip, uintptr_t* set_incl, uintptr_t set_ct, uint32_t** setdefs, char* outname, char* outname_end, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* sex_male, Chrom_info* chrom_info_ptr, uint32_t ignore_x, uint32_t hh_exists, uint32_t*** ld_map_ptr) {
  // Takes a bunch of set definitions, and determines which pairs of same-set
  // markers reach/exceed the --set-r2 threshold, saving them (in setdef
  // format) to a newly stack-allocated ld_map[].
  // If --set-r2 write was specified, the map's contents are written to {output
  // prefix}.ldset.
  // Note that, when very large set(s) are present, and there's a moderate
  // amount of "random" long-range LD, the memory requirement may be huge.
  FILE* outfile = nullptr;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t marker_ctv = ((marker_ct + 127) / 128) * (128 / BITCT);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t max_set_id_len = sip->max_name_len;
  uintptr_t founder_ct = popcount_longs(founder_pnm, unfiltered_sample_ctl);
  uintptr_t founder_ctl = BITCT_TO_WORDCT(founder_ct);
  uintptr_t founder_ctv2 = founder_ctl * 2;
  uintptr_t founder_ct_mld = (founder_ct + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
  uint32_t founder_ct_mld_m1 = ((uint32_t)founder_ct_mld) - 1;
#ifdef __LP64__
  uintptr_t founder_ct_mld_rem = (MULTIPLEX_LD / 192) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 192;
#else
  uintptr_t founder_ct_mld_rem = (MULTIPLEX_LD / 48) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 48;
#endif
  uintptr_t founder_ct_192_long = founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2) + founder_ct_mld_rem * (192 / BITCT2);
  uintptr_t final_mask = get_final_mask(founder_ct);
  uint32_t founder_trail_ct = founder_ct_192_long - founder_ctl * 2;
  uint32_t marker_idx = 0;
  uint32_t chrom_fo_idx = 0;
  uint32_t chrom_idx = 0;
  uint32_t is_haploid = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t range_end = 0;
  int32_t retval = 0;
  char charbuf[8];
  uintptr_t* loadbuf;
  uintptr_t* load2_bitfield;
  uintptr_t* founder_include2;
  uintptr_t* founder_male_include2;
  uintptr_t* tmp_set_bitfield;
  uintptr_t* geno1;
  uintptr_t* geno_masks1;
  uintptr_t* geno2;
  uintptr_t* geno_masks2;
  uintptr_t* result_bitfield;
  uintptr_t* rb_ptr;
  uintptr_t* loadbuf_ptr;
  uint32_t** ld_map;
  uint32_t* cur_setdef;
  uint32_t* cur_setdef2;
  char* sptr;
  char* wptr_start;
  char* wptr;
  uintptr_t memreq1;
  uintptr_t memreq2;
  uintptr_t minmem;
  uintptr_t idx1_block_size;
  uintptr_t idx2_block_size;
  uintptr_t cur_idx2_block_size;
  uintptr_t firstw;
  uintptr_t wlen;
  uintptr_t marker_uidx;
  uintptr_t marker_uidx2;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t thread_ct;
  uint32_t chrom_end;
  uint32_t set_idx;
  uint32_t set_uidx;
  uint32_t idx1_block_end;
  uint32_t marker_idx2;
  uint32_t load_idx2_tot;
  uint32_t marker_load_idx2;
  uint32_t block_idx1;
  uint32_t block_idx2;
  uint32_t setdef_incr_aux;
  uint32_t setdef_incr_aux2;
  uint32_t is_last_block;
  uint32_t range_start;
  uint32_t uii;
  if (!founder_ct) {
    logerrprint("Error: Cannot construct LD map, since there are no founders with nonmissing\nphenotypes.  (--make-founders may come in handy here.)\n");
    goto construct_ld_map_ret_INVALID_CMDLINE;
  }
  ld_map = (uint32_t**)bigstack_alloc(marker_ct * sizeof(intptr_t));
  if (!ld_map) {
    goto construct_ld_map_ret_NOMEM;
  }
  *ld_map_ptr = ld_map;
  // To avoid too much back-and-forth disk seeking for large datasets, we
  // construct the LD map in blocks, using similar logic to the --r/--r2 and
  // --fast-epistasis computations.
  // 1. bigstack_end_alloc space for main window markers' raw data, bitfields
  //    for them listing intersecting markers in front (i.e. we only look at
  //    the upper right triangle of the LD matrix), and another union bitfield.
  //    Break the union into secondary windows, and for each secondary window:
  //    a. bigstack_end_alloc secondary window markers' raw data
  //    b. perform multithreaded LD calculations, saving results via in-place
  //       clearing of the first markers' bitfields
  //    Memory requirement per main window marker is:
  //      96 bytes per 192 founders for raw data (rounded up)
  //      32 bytes per 128 filtered markers (rounded up), for the results (16
  //      working, 16 final)
  //      4 bytes for missing_ct
  //      16 extra bytes to ensure enough setdef compression workspace
  //    Memory req. per secondary window marker is 4 + 96 bytes/192 founders.
  //    To reduce false sharing risk, each thread is assigned at least 4
  //    markers.
  // 2. populate the bottom left triangle of the result matrix by referring to
  //    earlier results
  // 3. save final results for each marker in compressed setdef format at the
  //    current workspace bottom
  // 4. dump .ldset file if necessary
  loadbuf = (uintptr_t*)bigstack_end_alloc(unfiltered_sample_ct4);
  if (!loadbuf) {
    // separate since unfiltered_sample_ct4 is a byte, not word, count
    goto construct_ld_map_ret_NOMEM;
  }
  if (bigstack_end_alloc_ul(marker_ctv, &load2_bitfield) ||
      bigstack_end_alloc_ul(marker_ctv, &tmp_set_bitfield) ||
      bigstack_end_alloc_ul(founder_ctv2, &founder_include2) ||
      bigstack_end_alloc_ul(founder_ctv2, &founder_male_include2)) {
    goto construct_ld_map_ret_NOMEM;
  }
  // bugfix: last word might not be initialized by unpack_set().  Also
  // initialize second-to-last word to defend against an unpack_set()
  // implementation change.
#ifndef __LP64__
  // oh, this also matters in 32-bit case
  tmp_set_bitfield[marker_ctv - 4] = 0;
  tmp_set_bitfield[marker_ctv - 3] = 0;
#endif
  tmp_set_bitfield[marker_ctv - 2] = 0;
  tmp_set_bitfield[marker_ctv - 1] = 0;
  g_ld_load2_bitfield = load2_bitfield;
  alloc_collapsed_haploid_filters(founder_pnm, sex_male, unfiltered_sample_ct, founder_ct, XMHH_EXISTS | hh_exists, 1, &founder_include2, &founder_male_include2);
  memreq2 = founder_ct_192_long * sizeof(intptr_t) * 2 + 4;

  // this guarantees enough room for save_set_bitfield() worst case
  memreq1 = memreq2 + marker_ctv * sizeof(intptr_t) * 2 + 16;

  minmem = memreq2 * BITCT;
  if (minmem < memreq1 * 4) {
    minmem = memreq1 * 4;
  }
  g_ld_marker_ct = marker_ct;
  g_ld_founder_ct = founder_ct;
  g_ld_founder_ct_192_long = founder_ct_192_long;
  g_ld_founder_ct_mld_m1 = founder_ct_mld_m1;
  g_ld_founder_ct_mld_rem = founder_ct_mld_rem;
  g_ld_window_r2 = sip->set_r2 * (1 - SMALL_EPSILON);
  do {
    ulii = bigstack_left() / 2;
    if (ulii < minmem) {
      goto construct_ld_map_ret_NOMEM;
    }
    idx1_block_size = (ulii / memreq1) & (~(3 * ONELU));
    if (idx1_block_size > marker_ct - marker_idx) {
      idx1_block_size = marker_ct - marker_idx;
    }
    thread_ct = g_thread_ct;
    if (thread_ct > idx1_block_size / 4) {
      thread_ct = idx1_block_size / 4;
      if (!thread_ct) {
	thread_ct = 1;
      }
    }
    g_ld_thread_ct = thread_ct;
    idx2_block_size = (ulii / memreq2) & (~(BITCT - ONELU));
    if (idx2_block_size > marker_ct) {
      idx2_block_size = marker_ct;
    }
    g_ld_block_idx1 = marker_idx;
    g_ld_idx1_block_size = idx1_block_size;
    bigstack_end_alloc_ul(idx1_block_size * founder_ct_192_long, &geno1);
    bigstack_end_alloc_ul(idx1_block_size * founder_ct_192_long, &geno_masks1);
    bigstack_end_alloc_ui(idx1_block_size, &g_ld_missing_cts1);
    bigstack_end_alloc_ul(idx2_block_size * founder_ct_192_long, &geno2);
    bigstack_end_alloc_ul(idx2_block_size * founder_ct_192_long, &geno_masks2);
    bigstack_end_alloc_ui(idx2_block_size, &g_ld_missing_cts2);
    bigstack_end_alloc_ul(idx1_block_size * marker_ctv, &result_bitfield);
    uljj = founder_trail_ct + 2;
    for (ulii = 1; ulii <= idx1_block_size; ulii++) {
      fill_ulong_zero(uljj, &(geno1[ulii * founder_ct_192_long - uljj]));
      fill_ulong_zero(uljj, &(geno_masks1[ulii * founder_ct_192_long - uljj]));
    }
    for (ulii = 1; ulii <= idx2_block_size; ulii++) {
      fill_ulong_zero(uljj, &(geno2[ulii * founder_ct_192_long - uljj]));
      fill_ulong_zero(uljj, &(geno_masks2[ulii * founder_ct_192_long - uljj]));
    }
    fill_ulong_zero(idx1_block_size * marker_ctv, result_bitfield);
    g_ld_geno1 = geno1;
    g_ld_geno_masks1 = geno_masks1;
    g_ld_geno2 = geno2;
    g_ld_geno_masks2 = geno_masks2;
    g_ld_result_bitfield = result_bitfield;
    idx1_block_end = marker_idx + idx1_block_size;
    fill_ulong_zero(marker_ctv, load2_bitfield);
    fill_ulong_zero(idx1_block_size * marker_ctv, result_bitfield);
    for (set_idx = 0; set_idx < set_ct; set_idx++) {
      cur_setdef = setdefs[set_idx];
      setdef_iter_init(cur_setdef, marker_ct, marker_idx, &marker_idx2, &setdef_incr_aux);
      if (setdef_iter(cur_setdef, &marker_idx2, &setdef_incr_aux) && (marker_idx2 < idx1_block_end)) {
	unpack_set(marker_ct, cur_setdef, tmp_set_bitfield);
        get_set_wrange_align(tmp_set_bitfield, marker_ctv, &firstw, &wlen);
	if (wlen) {
	  uii = marker_idx2;
	  do {
	    bitvec_or(&(tmp_set_bitfield[firstw]), wlen, &(result_bitfield[((marker_idx2 - marker_idx) * marker_ctv + firstw)]));
	    marker_idx2++;
	    next_set_ck(tmp_set_bitfield, idx1_block_end, &marker_idx2);
	  } while (marker_idx2 < idx1_block_end);
	  // don't need to load the first intersecting member or anything
	  // before it, since we're only traversing the upper right triangle
	  wlen += firstw;
#ifdef __LP64__
	  firstw = 2 * (uii / 128);
#else
	  firstw = uii / 32;
#endif
	  clear_bits(0, uii + 1 - firstw * BITCT, &(tmp_set_bitfield[firstw]));
	  bitvec_or(&(tmp_set_bitfield[firstw]), wlen - firstw, &(load2_bitfield[firstw]));
	}
      }
    }
    load_idx2_tot = popcount_longs(load2_bitfield, marker_ctv);
    if (!load_idx2_tot) {
      // no new r^2 computations to make at all!
      goto construct_ld_map_no_new;
    }
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    if (marker_idx) {
      marker_uidx = jump_forward_unset_unsafe(marker_exclude, marker_uidx + 1, marker_idx);
    }
    marker_uidx2 = marker_uidx;
    if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto construct_ld_map_ret_READ_FAIL;
    }
    chrom_end = 0;
    for (block_idx1 = 0; block_idx1 < idx1_block_size; marker_uidx++, block_idx1++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
        if (fseeko(bedfile, bed_offset + (marker_uidx * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
          goto construct_ld_map_ret_READ_FAIL;
	}
      }
      if (marker_uidx >= chrom_end) {
        chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx);
        chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
        is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
        is_x = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[X_OFFSET]);
        is_y = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
      }
      ulii = block_idx1 * founder_ct_192_long;
      loadbuf_ptr = &(geno1[ulii]);
      if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_pnm, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf, loadbuf_ptr)) {
	goto construct_ld_map_ret_READ_FAIL;
      }
      if (is_haploid && hh_exists) {
        haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)loadbuf_ptr);
      }
      ld_process_load2(loadbuf_ptr, &(geno_masks1[ulii]), &(g_ld_missing_cts1[block_idx1]), founder_ct, is_x && (!ignore_x), founder_male_include2);
    }
    chrom_end = 0;
    cur_idx2_block_size = idx2_block_size;
    marker_idx2 = next_set_unsafe(load2_bitfield, 0);
    marker_uidx2 = jump_forward_unset_unsafe(marker_exclude, marker_uidx2 + 1, marker_idx2 - marker_idx);
    marker_load_idx2 = 0;
    if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
      goto construct_ld_map_ret_READ_FAIL;
    }
    do {
      if (cur_idx2_block_size > load_idx2_tot - marker_load_idx2) {
	cur_idx2_block_size = load_idx2_tot - marker_load_idx2;
      }
      g_ld_idx2_block_start = marker_idx2;
      block_idx2 = 0;
      while (1) {
	if (marker_uidx2 >= chrom_end) {
	  chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx2);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
	  is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
	  is_x = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[X_OFFSET]);
	  is_y = (((int32_t)chrom_idx) == chrom_info_ptr->xymt_codes[Y_OFFSET]);
	}
	ulii = block_idx2 * founder_ct_192_long;
	loadbuf_ptr = &(geno2[ulii]);
	if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_pnm, final_mask, IS_SET(marker_reverse, marker_uidx2), bedfile, loadbuf, loadbuf_ptr)) {
	  goto construct_ld_map_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)loadbuf_ptr);
	}
	ld_process_load2(loadbuf_ptr, &(geno_masks2[ulii]), &(g_ld_missing_cts2[block_idx2]), founder_ct, is_x && (!ignore_x), founder_male_include2);
	if (++block_idx2 == cur_idx2_block_size) {
	  break;
	}
        uii = marker_idx2++;
	ulii = ++marker_uidx2;
        if (is_set(load2_bitfield, marker_idx2)) {
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx2);
	} else {
          marker_idx2 = next_set_unsafe(load2_bitfield, marker_idx2);
          marker_uidx2 = jump_forward_unset_unsafe(marker_exclude, marker_uidx2, marker_idx2 - uii);
	}
	if (ulii < marker_uidx2) {
	  if (fseeko(bedfile, bed_offset + (marker_uidx2 * ((uint64_t)unfiltered_sample_ct4)), SEEK_SET)) {
	    goto construct_ld_map_ret_READ_FAIL;
	  }
	}
      }
      g_ld_marker_ctm8 = marker_idx2 + 1; // repurposed
      marker_load_idx2 += cur_idx2_block_size;
      is_last_block = (marker_load_idx2 == load_idx2_tot);
      if (spawn_threads2(threads, &ld_map_thread, thread_ct, is_last_block)) {
	goto construct_ld_map_ret_THREAD_CREATE_FAIL;
      }
      ld_map_thread((void*)0);
      join_threads2(threads, thread_ct, is_last_block);
    } while (!is_last_block);
  construct_ld_map_no_new:
    for (block_idx1 = marker_idx; block_idx1 < idx1_block_end; block_idx1++) {
      rb_ptr = &(result_bitfield[(block_idx1 - marker_idx) * marker_ctv]);
      marker_idx2 = 0;
      while (1) {
	marker_idx2 = next_set(rb_ptr, marker_idx2, block_idx1);
	if (marker_idx2 == block_idx1) {
	  clear_bit(block_idx1, rb_ptr);
	  break;
	}
	if (!in_setdef(ld_map[marker_idx2], block_idx1)) {
	  clear_bit(marker_idx2, rb_ptr);
	}
	marker_idx2++;
      }
      range_start = next_set(rb_ptr, 0, marker_ct);
      if (range_start != marker_ct) {
	range_end = last_set_bit(rb_ptr, marker_ctv) + 1;
      }
      save_set_bitfield(rb_ptr, marker_ct, range_start, range_end, 0, &(ld_map[block_idx1]));
    }
    // free previous round of allocations
    bigstack_end_reset(founder_male_include2);
    marker_idx = idx1_block_end;
  } while (marker_idx < marker_ct);
  if (sip->modifier & SET_R2_WRITE) {
    memcpy(charbuf, outname_end, 8);
    memcpy(outname_end, ".ldset", 7);
    if (fopen_checked(outname, "w", &outfile)) {
      goto construct_ld_map_ret_OPEN_FAIL;
    }
    set_uidx = 0;
    for (set_idx = 0; set_idx < set_ct; set_uidx++, set_idx++) {
      next_set_unsafe_ck(set_incl, &set_uidx);
      sptr = &(sip->names[set_uidx * max_set_id_len]);
      uii = strlen(sptr);
      wptr_start = memcpyax(g_textbuf, sptr, uii, ' ');
      cur_setdef = setdefs[set_idx];
      setdef_iter_init(cur_setdef, marker_ct, 0, &marker_idx, &setdef_incr_aux);

      while (setdef_iter(cur_setdef, &marker_idx, &setdef_incr_aux)) {
        cur_setdef2 = ld_map[marker_idx];
	// cur_setdef2 can contain variants outside of the current set, so we
	// need to look at the intersection.
        setdef_iter_init(cur_setdef2, marker_ct, 0, &marker_idx2, &setdef_incr_aux2);
        uii = 0; // now this tracks whether a first match has been found
	while (setdef_iter(cur_setdef2, &marker_idx2, &setdef_incr_aux2)) {
	  if (in_setdef(cur_setdef, marker_idx2)) {
	    if (!uii) {
	      uii = 1;
	      wptr = strcpyax(wptr_start, &(marker_ids[marker_idx_to_uidx[marker_idx] * max_marker_id_len]), ' ');
	      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
		goto construct_ld_map_ret_WRITE_FAIL;
	      }
	    }
	    fputs(&(marker_ids[marker_idx_to_uidx[marker_idx2] * max_marker_id_len]), outfile);
	    putc_unlocked(' ', outfile);
	  }
	  marker_idx2++;
	}
	if (uii) {
	  if (putc_checked('\n', outfile)) {
	    goto construct_ld_map_ret_WRITE_FAIL;
	  }
	}
        marker_idx++;
      }
    }
    if (fclose_null(&outfile)) {
      goto construct_ld_map_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("--set-r2 write: LD map written to %s .\n", outname);
    memcpy(outname_end, charbuf, 8);
  } else {
    logprint("LD map constructed.\n");
  }
  while (0) {
  construct_ld_map_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  construct_ld_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  construct_ld_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  construct_ld_map_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  construct_ld_map_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  construct_ld_map_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
  fclose_cond(outfile);
  bigstack_end_reset(bigstack_end_mark);
  return retval;
}

void set_test_score(uintptr_t marker_ct, double chisq_threshold, uint32_t set_max, double* chisq_arr, uint32_t** ld_map, uint32_t* cur_setdef, double* sorted_chisq_buf, uint32_t* sorted_marker_idx_buf, uint32_t* proxy_arr, uint32_t* raw_sig_ct_ptr, uint32_t* final_sig_ct_ptr, double* set_score_ptr) {
  // set score statistic = mean of chi-square statistics of set
  // representatives.  --linear t statistics are converted to same-p-value 1df
  // chi-square stats out of necessity; in theory, this hack could be applied
  // to e.g. Fisher's exact test and the variable-df genotypic test as well,
  // but I'll hold off on that until/unless it's specifically requested.

  // sort variants by p-value, then iterate over setdefs, greedily selecting up
  // to sip->set_max significant independent variants from each.
  double chi_sum = 0.0;
  uint32_t raw_sig_ct = 0;
  uint32_t final_sig_ct = 0;
  uint32_t marker_idx;
  uint32_t setdef_incr_aux;
  uint32_t raw_idx;
  uint32_t ld_conflict;
  uint32_t uii;
  setdef_iter_init(cur_setdef, marker_ct, 0, &marker_idx, &setdef_incr_aux);
  while (setdef_iter(cur_setdef, &marker_idx, &setdef_incr_aux)) {
    if (chisq_arr[marker_idx] >= chisq_threshold) {
      sorted_chisq_buf[raw_sig_ct] = chisq_arr[marker_idx];
      sorted_marker_idx_buf[raw_sig_ct] = marker_idx;
      raw_sig_ct++;
    }
    marker_idx++;
  }
  if (!raw_sig_ct) {
    // not possible for initial pass, so no need to set raw_sig_ct_ptr, etc.
    // bugfix: actually, that comment was incorrect
    *set_score_ptr = 0.0;
    return;
  }
  qsort_ext2((char*)sorted_chisq_buf, raw_sig_ct, sizeof(double), double_cmp_deref, (char*)sorted_marker_idx_buf, sizeof(int32_t), (char*)proxy_arr, sizeof(double) + sizeof(int32_t));
  raw_idx = raw_sig_ct;
  do {
    raw_idx--;
    ld_conflict = 0;
    marker_idx = sorted_marker_idx_buf[raw_idx];
    for (uii = 0; uii < final_sig_ct; uii++) {
      if (in_setdef(ld_map[proxy_arr[uii]], marker_idx)) {
	ld_conflict = 1;
	break;
      }
    }
    if (!ld_conflict) {
      proxy_arr[final_sig_ct] = marker_idx;
      chi_sum += sorted_chisq_buf[raw_idx];
      if (++final_sig_ct == set_max) {
	break;
      }
    }
  } while (raw_idx);
  *set_score_ptr = chi_sum / ((double)((int32_t)final_sig_ct));
  if (raw_sig_ct_ptr) {
    *raw_sig_ct_ptr = raw_sig_ct;
    *final_sig_ct_ptr = final_sig_ct;
  }
}

int32_t set_test_common_init(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_reverse, double* orig_chisq, Set_info* sip, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sex_male, uintptr_t* founder_pnm, uint32_t ld_ignore_x, uint32_t hh_exists, const char* flag_descrip, uintptr_t* marker_ct_ptr, uintptr_t** marker_exclude_ptr, uintptr_t** set_incl_ptr, uint32_t** marker_idx_to_uidx_ptr, uint32_t*** setdefs_ptr, uintptr_t* set_ct_ptr, uint32_t* max_sigset_size_ptr, uint32_t*** ld_map_ptr, double* chisq_threshold_ptr, double** orig_set_scores_ptr, double** sorted_chisq_buf_ptr, uint32_t** sorted_marker_idx_buf_ptr, uint32_t** proxy_arr_ptr, uintptr_t** perm_adapt_set_unstopped_ptr, uint32_t** perm_2success_ct_ptr, uint32_t** perm_attempt_ct_ptr, uintptr_t** unstopped_markers_ptr) {
  // Assumes *marker_ct_ptr has value marker_ct_mid, and marker_exclude_ptr
  // initially points to marker_exclude_mid.
  // Side effect: allocates set_incl, marker_idx_to_uidx, ld_map, and several
  // other arrays on stack
  uintptr_t marker_ct_mid = *marker_ct_ptr;
  uintptr_t marker_ct = marker_ct_mid;
  uintptr_t raw_set_ct = sip->ct;
  uintptr_t raw_set_ctl = BITCT_TO_WORDCT(raw_set_ct);
  uintptr_t set_ct = 0;
  uintptr_t* marker_exclude_mid = *marker_exclude_ptr;
  double chisq_threshold = inverse_chiprob(sip->set_p, 1);
  uint32_t max_sigset_size = 0;
  int32_t retval = 0;
  uintptr_t marker_midx;
  uintptr_t set_uidx;
  uintptr_t* set_incl;
  uintptr_t* cur_bitfield;
  double* chisq_ptr;
  double* chisq_end;
  double* orig_set_scores;
  uint32_t** setdefs;
  uint32_t* marker_midx_to_idx;
  uint32_t* cur_setdef;
  uintptr_t marker_idx;
  uintptr_t set_idx;
  uint32_t range_ct;
  uint32_t range_idx;
  uint32_t range_offset;
  uint32_t range_stop;
  uint32_t include_out_of_bounds;
  uint32_t cur_set_size;
  uint32_t cur_range_size;
  uint32_t uii;
  if (bigstack_calloc_ul(raw_set_ctl, set_incl_ptr) ||
      bigstack_alloc_ui(marker_ct_orig, &marker_midx_to_idx)) {
    goto set_test_common_init_ret_NOMEM;
  }
  set_incl = *set_incl_ptr;
  fill_midx_to_idx(marker_exclude_orig, marker_exclude_mid, marker_ct, marker_midx_to_idx);

  // determine which sets contain at least one significant marker.  do not
  // attempt to calculate the sum statistic yet: we need the LD map for that.
  for (set_uidx = 0; set_uidx < raw_set_ct; set_uidx++) {
    cur_setdef = sip->setdefs[set_uidx];
    range_ct = cur_setdef[0];
    cur_set_size = 0;
    uii = 0; // found a significant marker?
    if (range_ct != 0xffffffffU) {
      for (range_idx = 0; range_idx < range_ct; range_idx++) {
        marker_midx = *(++cur_setdef);
        range_stop = *(++cur_setdef);
	cur_range_size = range_stop - marker_midx;
        cur_set_size += cur_range_size;
        if (!uii) {
          chisq_ptr = &(orig_chisq[marker_midx_to_idx[marker_midx]]);
          chisq_end = &(chisq_ptr[cur_range_size]);
	  for (; chisq_ptr < chisq_end; chisq_ptr++) {
	    if (*chisq_ptr >= chisq_threshold) {
	      uii = 1;
              break;
	    }
	  }
	}
      }
    } else {
      range_offset = cur_setdef[1];
      range_stop = cur_setdef[2];
      include_out_of_bounds = cur_setdef[3];
      cur_bitfield = (uintptr_t*)(&(cur_setdef[4]));
      if (include_out_of_bounds && range_offset) {
        for (marker_midx = 0; marker_midx < range_offset; marker_midx++) {
	  // all initial markers guaranteed to be in union, no
	  // marker_midx_to_idx lookup needed
          if (orig_chisq[marker_midx] >= chisq_threshold) {
            uii = 1;
            break;
	  }
	}
        cur_set_size += range_offset;
      }
      cur_set_size += popcount_longs(cur_bitfield, ((range_stop + 127) / 128) * (128 / BITCT));
      if (!uii) {
        for (marker_midx = 0; marker_midx < range_stop; marker_midx++) {
          if (IS_SET(cur_bitfield, marker_midx)) {
            if (orig_chisq[marker_midx_to_idx[marker_midx + range_offset]] >= chisq_threshold) {
              uii = 1;
              break;
	    }
	  }
	}
      }
      if (include_out_of_bounds && (range_offset + range_stop < marker_ct_orig)) {
        cur_set_size += marker_ct_orig - range_offset - range_stop;
        if (!uii) {
          for (marker_idx = marker_midx_to_idx[range_offset + range_stop]; marker_idx < marker_ct; marker_idx++) {
	    // all trailing markers guaranteed to be in union
            if (orig_chisq[marker_idx] >= chisq_threshold) {
              uii = 1;
              break;
	    }
	  }
	}
      }
    }
    if (uii) {
      SET_BIT(set_uidx, set_incl);
      set_ct++;
      if (cur_set_size > max_sigset_size) {
	max_sigset_size = cur_set_size;
      }
    }
  }
  if (!set_ct) {
    logerrprint("Warning: No significant variants in any set.  Skipping permutation-based set\ntest.\n");
    goto set_test_common_init_ret_1;
  }
  LOGPRINTFWW("%s set test: Testing %" PRIuPTR " set%s with at least one significant variant.\n", flag_descrip, set_ct, (set_ct == 1)? "" : "s");
  bigstack_reset((unsigned char*)marker_midx_to_idx);
  if (set_ct < raw_set_ct) {
    marker_ct = marker_ct_orig;
    if (extract_set_union_unfiltered(sip, set_incl, unfiltered_marker_ct, marker_exclude_orig, marker_exclude_ptr, &marker_ct)) {
      goto set_test_common_init_ret_NOMEM;
    }
  }
  // Okay, we've pruned all we can, now it's time to suck it up and construct
  // the potentially huge LD map
  if (bigstack_alloc_ui(marker_ct, marker_idx_to_uidx_ptr)) {
    goto set_test_common_init_ret_NOMEM;
  }
  fill_idx_to_uidx(*marker_exclude_ptr, unfiltered_marker_ct, marker_ct, *marker_idx_to_uidx_ptr);
  if (marker_ct < marker_ct_orig) {
    if (setdefs_compress(sip, set_incl, set_ct, unfiltered_marker_ct, marker_exclude_orig, marker_ct_orig, *marker_exclude_ptr, marker_ct, setdefs_ptr)) {
      goto set_test_common_init_ret_NOMEM;
    }
  } else {
    *setdefs_ptr = sip->setdefs;
  }
  setdefs = *setdefs_ptr;
  retval = construct_ld_map(threads, bedfile, bed_offset, *marker_exclude_ptr, marker_ct, marker_reverse, *marker_idx_to_uidx_ptr, unfiltered_sample_ct, founder_pnm, sip, set_incl, set_ct, setdefs, outname, outname_end, marker_ids, max_marker_id_len, sex_male, chrom_info_ptr, ld_ignore_x, hh_exists, ld_map_ptr);
  if (retval) {
    goto set_test_common_init_ret_1;
  }
  if (marker_ct_mid != marker_ct) {
    // caller needs to collapse other arrays
    inplace_delta_collapse_arr((char*)orig_chisq, sizeof(double), marker_ct_mid, marker_ct, marker_exclude_mid, *marker_exclude_ptr);
  }
  if (bigstack_alloc_d(set_ct, orig_set_scores_ptr) ||
      bigstack_alloc_d(max_sigset_size, sorted_chisq_buf_ptr) ||
      bigstack_alloc_ui(max_sigset_size, sorted_marker_idx_buf_ptr) ||
      // 3 int32s = max(sizeof(double), sizeof(intptr_t)) + sizeof(int32_t)
      bigstack_alloc_ui(max_sigset_size * 3LU, proxy_arr_ptr)) {
    goto set_test_common_init_ret_NOMEM;
  }
  orig_set_scores = *orig_set_scores_ptr;
  for (set_idx = 0; set_idx < set_ct; set_idx++) {
    // we're calling this again during final write anyway, so don't bother
    // saving raw_sig_ct or final_sig_ct now
    set_test_score(marker_ct, chisq_threshold, sip->set_max, orig_chisq, *ld_map_ptr, setdefs[set_idx], *sorted_chisq_buf_ptr, *sorted_marker_idx_buf_ptr, *proxy_arr_ptr, nullptr, nullptr, &(orig_set_scores[set_idx]));
  }
  // just treat --mperm as --perm with min_perms == max_perms, since this isn't
  // a proper max(T) test
  if (bigstack_alloc_ul(BITCT_TO_WORDCT(set_ct), perm_adapt_set_unstopped_ptr) ||
      bigstack_calloc_ui(set_ct, perm_2success_ct_ptr) ||
      bigstack_alloc_ui(set_ct, perm_attempt_ct_ptr) ||
      bigstack_alloc_ul(BITCT_TO_WORDCT(marker_ct), unstopped_markers_ptr)) {
    goto set_test_common_init_ret_NOMEM;
  }
  fill_all_bits(set_ct, *perm_adapt_set_unstopped_ptr);
  fill_all_bits(marker_ct, *unstopped_markers_ptr);
  while (0) {
  set_test_common_init_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
 set_test_common_init_ret_1:
  *marker_ct_ptr = marker_ct;
  *set_ct_ptr = set_ct;
  *max_sigset_size_ptr = max_sigset_size;
  *chisq_threshold_ptr = chisq_threshold;
  return retval;
}

void compute_set_scores(uintptr_t marker_ct, uintptr_t perm_vec_ct, uintptr_t set_ct, double* chisq_matrix, double* orig_set_scores, double* sorted_chisq_buf, uint32_t* sorted_marker_idx_buf, uint32_t* proxy_arr, uint32_t** setdefs, uint32_t** ld_map, Aperm_info* apip, double chisq_threshold, double adaptive_ci_zt, uint32_t first_adapt_check, uint32_t perms_done, uint32_t set_max, uintptr_t* perm_adapt_set_unstopped, uint32_t* perm_2success_ct, uint32_t* perm_attempt_ct) {
  // compute set stats for the just-completed permutations
  uint32_t pidx_offset = perms_done - perm_vec_ct;
  uintptr_t set_idx;
  double stat_high;
  double stat_low;
  double cur_score;
  double pval;
  double dxx;
  uint32_t next_adapt_check;
  uint32_t pidx;
  uint32_t uii;
  for (set_idx = 0; set_idx < set_ct; set_idx++) {
    if (IS_SET(perm_adapt_set_unstopped, set_idx)) {
      next_adapt_check = first_adapt_check;
      uii = perm_2success_ct[set_idx];
      stat_high = orig_set_scores[set_idx] + EPSILON;
      stat_low = orig_set_scores[set_idx] - EPSILON;
      for (pidx = 0; pidx < perm_vec_ct;) {
	set_test_score(marker_ct, chisq_threshold, set_max, &(chisq_matrix[pidx * marker_ct]), ld_map, setdefs[set_idx], sorted_chisq_buf, sorted_marker_idx_buf, proxy_arr, nullptr, nullptr, &cur_score);
	if (cur_score > stat_high) {
	  uii += 2;
	} else if (cur_score > stat_low) {
	  uii++;
	}
	if (++pidx == next_adapt_check - pidx_offset) {
	  if (uii) {
	    pval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	    if ((pval - dxx > apip->alpha) || (pval + dxx < apip->alpha)) {
	      CLEAR_BIT(set_idx, perm_adapt_set_unstopped);
	      perm_attempt_ct[set_idx] = next_adapt_check;
	      break;
	    }
	  }
	  next_adapt_check += (int32_t)(apip->init_interval + ((int32_t)next_adapt_check) * apip->interval_slope);
	}
      }
      perm_2success_ct[set_idx] = uii;
    }
  }
}

int32_t write_set_test_results(char* outname, char* outname_end2, Set_info* sip, uint32_t** ld_map, uint32_t** setdefs, uintptr_t* set_incl, uintptr_t set_ct, uintptr_t marker_ct_orig, uintptr_t marker_ct, uint32_t* marker_idx_to_uidx, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* perm_2success_ct, uint32_t* perm_attempt_ct, uint32_t mtest_adjust, uint32_t perm_count, double pfilter, double output_min_p, double chisq_threshold, double* orig_stats, double* sorted_chisq_buf, uint32_t* sorted_marker_idx_buf, uint32_t* proxy_arr) {
  // assumes caller will free memory from stack
  FILE* outfile = nullptr;
  uintptr_t* nonempty_set_incl = nullptr;
  double* empirical_pvals = nullptr;
  uintptr_t raw_set_ct = sip->ct;
  uintptr_t max_set_id_len = sip->max_name_len;
  uint32_t nonempty_set_ct = 0;
  int32_t retval = 0;
  uintptr_t set_uidx;
  uintptr_t set_idx;
  char* bufptr;
  uint32_t* nonempty_set_idx_to_uidx;
  double cur_score;
  double pval;
  uint32_t raw_sig_ct;
  uint32_t final_sig_ct;
  uint32_t set_midx;
  uint32_t uii;
  if (set_ct && mtest_adjust) {
    if (alloc_and_populate_nonempty_set_incl(sip, &nonempty_set_ct, &nonempty_set_incl)) {
      goto write_set_test_results_ret_NOMEM;
    }
    if (bigstack_alloc_d(nonempty_set_ct, &empirical_pvals)) {
      goto write_set_test_results_ret_NOMEM;
    }
  }
  if (fopen_checked(outname, "w", &outfile)) {
    goto write_set_test_results_ret_OPEN_FAIL;
  }
  fprintf(outfile, "         SET   NSNP   NSIG   ISIG         EMP1 %sSNPS\n", perm_count? "          NP " : "");
  for (set_uidx = 0, set_midx = 0, set_idx = 0; set_uidx < raw_set_ct; set_uidx++) {
    bufptr = fw_strcpy(12, &(sip->names[set_uidx * max_set_id_len]), g_textbuf);
    *bufptr++ = ' ';
    bufptr = uint32toa_w6x(setdef_size(sip->setdefs[set_uidx], marker_ct_orig), ' ', bufptr);
    if (IS_SET(set_incl, set_uidx)) {
      set_test_score(marker_ct, chisq_threshold, sip->set_max, orig_stats, ld_map, setdefs[set_idx], sorted_chisq_buf, sorted_marker_idx_buf, proxy_arr, &raw_sig_ct, &final_sig_ct, &cur_score);
      bufptr = uint32toa_w6x(raw_sig_ct, ' ', bufptr);
      bufptr = uint32toa_w6x(final_sig_ct, ' ', bufptr);
      pval = ((double)(perm_2success_ct[set_idx] + 2)) / ((double)(2 * (perm_attempt_ct[set_idx] + 1)));
      if (empirical_pvals) {
	empirical_pvals[set_midx] = pval;
      }
      if (pval <= pfilter) {
	if (!perm_count) {
	  bufptr = dtoa_g_wxp4x(MAXV(pval, output_min_p), 12, ' ', bufptr);
	} else {
	  bufptr = dtoa_g_wxp4(((double)perm_2success_ct[set_idx]) * 0.5, 12, bufptr);
	  bufptr = memseta(bufptr, 32, 3);
	  bufptr = uint32toa_w10x(perm_attempt_ct[set_idx], ' ', bufptr);
	}
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	  goto write_set_test_results_ret_WRITE_FAIL;
	}
	fputs(&(marker_ids[marker_idx_to_uidx[proxy_arr[0]] * max_marker_id_len]), outfile);
	for (uii = 1; uii < final_sig_ct; uii++) {
	  putc_unlocked('|', outfile);
	  fputs(&(marker_ids[marker_idx_to_uidx[proxy_arr[uii]] * max_marker_id_len]), outfile);
	}
	if (putc_checked('\n', outfile)) {
	  goto write_set_test_results_ret_WRITE_FAIL;
	}
      }
      set_midx++;
      set_idx++;
    } else {
      if (!perm_count) {
        bufptr = memcpya(bufptr, "     0      0            1 NA\n", 30);
      } else {
        bufptr = memcpya(bufptr, "     0      0            0            0 NA\n", 43);
      }
      if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	goto write_set_test_results_ret_WRITE_FAIL;
      }
      if (nonempty_set_incl && is_set(nonempty_set_incl, set_uidx)) {
	empirical_pvals[set_midx] = 1.0;
	set_midx++;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto write_set_test_results_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("Set test results written to %s .\n", outname);
  if (empirical_pvals) {
    if (bigstack_alloc_ui(nonempty_set_ct, &nonempty_set_idx_to_uidx)) {
      goto write_set_test_results_ret_NOMEM;
    }
    fill_idx_to_uidx_incl(nonempty_set_incl, raw_set_ct, nonempty_set_ct, nonempty_set_idx_to_uidx);
    // .qassoc.set.adjusted instead of .set.mperm.adjusted, etc.
    *outname_end2 = '\0';
    retval = multcomp(outname, outname_end2, nonempty_set_idx_to_uidx, nonempty_set_ct, sip->names, max_set_id_len, 0, nullptr, nullptr, pfilter, output_min_p, mtest_adjust, 1, 0.0, nullptr, empirical_pvals);
  }
  while (0) {
  write_set_test_results_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_set_test_results_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_set_test_results_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

typedef struct clump_entry_struct {
  double pval;
  struct clump_entry_struct* next;
  uint32_t fidx;
  char annot[];
} Clump_entry;

typedef struct cur_clump_info_struct {
  double r2;
  uint32_t marker_idx;
  uint32_t fidx;
} Cur_clump_info;

typedef struct clump_missing_id_struct {
  double pval;
  struct clump_missing_id_struct* next;
  char idstr[];
} Clump_missing_id;

void update_clump_histo(double pval, uintptr_t* histo) {
  if (pval < 0.001) {
    if (pval < 0.0001) {
      histo[4] += 1;
    } else {
      histo[3] += 1;
    }
  } else if (pval < 0.01) {
    histo[2] += 1;
  } else if (pval < 0.05) {
    histo[1] += 1;
  } else {
    histo[0] += 1;
  }
}

int32_t clump_reports(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, Clump_info* clump_ip, uintptr_t* sex_male, uint32_t hh_exists) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  gzFile gz_infile = nullptr;
  FILE* outfile = nullptr;
  FILE* outfile_ranges = nullptr;
  FILE* outfile_best = nullptr;
  uintptr_t marker_ctl = BITCT_TO_WORDCT(marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t founder_ct = popcount_longs(founder_info, unfiltered_sample_ctl);
  uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(founder_ct);
  uintptr_t founder_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(founder_ct);
  uintptr_t final_mask = get_final_mask(founder_ct);
  uintptr_t range_group_ct = 0;
  uintptr_t max_range_group_id_len = 0;
  uintptr_t max_header_len = 2;
  uintptr_t snpfield_search_ct = 1;
  uintptr_t pfield_search_ct = 1;
  uintptr_t annot_ct = 0;
  uintptr_t missing_variant_ct = 0;
  uintptr_t cur_rg_ct = 0;
  uintptr_t range_chrom_max = 0;
  uintptr_t unmatched_group_ct = 0;
  uintptr_t* haploid_mask = chrom_info_ptr->haploid_mask;
  char* range_group_names = nullptr;
  char* fname_ptr = nullptr;
  char* annot_flattened = clump_ip->annotate_flattened;
  char* tbuf2 = &(g_textbuf[MAXLINELEN]);
  char* header2_ptr = nullptr;
  char* annot_ptr = nullptr;
  char* cur_rg_names = nullptr;
  uintptr_t* founder_include2 = nullptr;
  uintptr_t* founder_male_include2 = nullptr;
  uintptr_t* rg_chrom_bounds = nullptr;
  uint32_t** rg_setdefs = nullptr;
  uint32_t** cur_rg_setdefs = nullptr;
  Clump_missing_id* not_found_list = nullptr;
  uintptr_t* rangematch_bitfield = nullptr;
  double p1_thresh = clump_ip->p1;
  double p2_thresh = clump_ip->p2;
  double load_pthresh = 0.05;
  double r2_thresh = clump_ip->r2;
  uint32_t allow_overlap = clump_ip->modifier & CLUMP_ALLOW_OVERLAP;
  uint32_t clump_index_first = clump_ip->modifier & CLUMP_INDEX_FIRST;
  uint32_t clump_best = clump_ip->modifier & CLUMP_BEST;
  uint32_t clump_verbose = clump_ip->modifier & CLUMP_VERBOSE;
  uint32_t bp_radius = clump_ip->bp_radius;
  uint32_t best_fidx_match = 0xffffffffU;
  uint32_t require_multifile = clump_ip->modifier & CLUMP_REPLICATE;
  uint32_t index_eligible = 1;
  uint32_t header1_len = 0;
  uint32_t header2_len = 0;
  uint32_t file_ct = 0;
  uint32_t final_clump_ct = 0;
  uint32_t max_missing_id_len = 0;
  int32_t retval = 0;
  uintptr_t histo[5]; // NSIG, S05, S01, S001, S0001
  uint32_t index_tots[5];
  uint32_t counts[18];
  Clump_entry** clump_entries;
  Clump_entry* clump_entry_ptr;
  Clump_entry* best_entry_ptr;
  Cur_clump_info* cur_clump_base;
  Cur_clump_info* cc_ptr;
  uintptr_t* col_bitfield;
  uintptr_t* cur_bitfield;
  uintptr_t* loadbuf_raw;
  uintptr_t* index_data;
  uintptr_t* window_data;
  uintptr_t* window_data_ptr;
  char* sorted_missing_variant_ids;
  char* sorted_header_dict;
  char* loadbuft; // t is for text
  char* cur_a1;
  char* cur_a2;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  uint32_t* header_id_map;
  uint32_t* marker_id_htable;
  uint32_t* parse_table;
  uint32_t* cur_parse_info;
  uint32_t* nsig_arr;
  uint32_t* pval_map;
  uint32_t* marker_uidx_to_idx;
  uint32_t* marker_idx_to_uidx;
  double* sorted_pvals;
  Clump_missing_id* cm_ptr;
  uintptr_t header_dict_ct;
  uintptr_t extra_annot_space;
  uintptr_t cur_bigstack_left;
  uintptr_t loadbuft_size;
  uintptr_t marker_idx;
  uintptr_t last_marker_idx;
  uintptr_t max_window_size; // universal bound
  uintptr_t cur_window_size;
  uintptr_t line_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  double pval;
  double freq1x;
  double freq2x;
  double freqx1;
  double freqx2;
  double freq11;
  double freq11_expected;
  double cur_r2;
  double max_r2;
  double dxx;
  uint32_t marker_id_htable_size;
  uint32_t annot_ct_p2;
  uint32_t annot_ct_p2_ctl;
  uint32_t cur_read_ct;
  uint32_t index_ct;
  uint32_t sp_idx;
  uint32_t file_idx;
  uint32_t ivar_idx;
  uint32_t ivar_uidx;
  uint32_t cur_bp;
  uint32_t min_bp;
  uint32_t max_bp;
  uint32_t clump_chrom_idx;
  uint32_t clump_uidx_first;
  uint32_t clump_uidx_last;
  uint32_t index_fidx;
  uint32_t marker_uidx;
  uint32_t max_r2_uidx;
  uint32_t is_haploid;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t a1_len;
  uint32_t a2_len;
  uint32_t allele_padding;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  int32_t ii;
  // suppress warning
  index_tots[3] = 0;
  index_tots[4] = 0;

  if (annot_flattened && (!clump_verbose) && (!clump_best)) {
    logerrprint("Error: --clump-annotate must be used with --clump-verbose or --clump-best.\n");
    goto clump_reports_ret_INVALID_CMDLINE;
  }
  if (!founder_ct) {
    logerrprint("Warning: Skipping --clump since there are no founders.  (--make-founders may\ncome in handy here.)\n");
    goto clump_reports_ret_1;
  }
  if (clump_best) {
    load_pthresh = 1.0;
  } else {
    if (p2_thresh > load_pthresh) {
      load_pthresh = p2_thresh;
    }
    if (p1_thresh >= load_pthresh) {
      // may as well maximize backwards compatibility re: which comparisons are
      // > vs. >=
      load_pthresh = p1_thresh * (1 + SMALL_EPSILON);
    }
  }
  if (clump_ip->range_fname) {
    // 1. load range file, sort, etc.
    retval = load_range_list_sortpos(clump_ip->range_fname, clump_ip->range_border, 0, nullptr, 0, chrom_info_ptr, &range_group_ct, &range_group_names, &max_range_group_id_len, &rg_chrom_bounds, &rg_setdefs, &range_chrom_max, "--clump-range");
    if (retval) {
      goto clump_reports_ret_1;
    }
  }
  // 2. create marker ID hash table, allocate index-tracking bitfield
  retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, 0, &marker_id_htable_size, &marker_id_htable);
  if (retval) {
    goto clump_reports_ret_1;
  }
  if (bigstack_calloc_ul(marker_ctl, &cur_bitfield) ||
      bigstack_alloc_ui(unfiltered_marker_ct, &marker_uidx_to_idx)) {
    goto clump_reports_ret_NOMEM;
  }
  fill_uidx_to_idx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_uidx_to_idx);
  if (clump_ip->snpfield_search_order) {
    snpfield_search_ct = count_and_measure_multistr(clump_ip->snpfield_search_order, &max_header_len);
  } else {
    max_header_len = 4; // 'SNP' + null terminator
  }
  if (clump_ip->pfield_search_order) {
    pfield_search_ct = count_and_measure_multistr(clump_ip->pfield_search_order, &max_header_len);
  }
  if (annot_flattened) {
    annot_ct = count_and_measure_multistr(annot_flattened, &max_header_len);
  }
  header_dict_ct = snpfield_search_ct + pfield_search_ct + annot_ct;
  // parse_table[2k + 1] stores the number of additional fields to skip before
  // reading that particular entry.  For example, if variant IDs are in the
  // second column in the current file, while p-values are in the fifth column,
  // parse_table[1] is 1 and parse_table[3] = 2.
  // parse_table[2k] stores the type of field contents (0 = variant ID, 1 =
  // P-value, 2 or more = annotation).
  // In the main loop, cur_parse_info[2k] stores the in-loadbuft offset of the
  // the string with that parse_table[2k] index, and cur_parse_info[2k + 1]
  // stores string length.
  annot_ct_p2 = 2 + annot_ct;
  annot_ct_p2_ctl = (annot_ct + (BITCT + 1)) / BITCT;
  if (bigstack_alloc_c(max_header_len * header_dict_ct, &sorted_header_dict) ||
      bigstack_alloc_ui(header_dict_ct, &header_id_map) ||
      bigstack_alloc_ul(annot_ct_p2_ctl, &col_bitfield) ||
      bigstack_alloc_ui(annot_ct_p2 * 2, &parse_table) ||
      bigstack_alloc_ui(annot_ct_p2 * 2, &cur_parse_info)) {
    goto clump_reports_ret_NOMEM;
  }
  ulii = 0; // write position
  if (clump_ip->snpfield_search_order) {
    bufptr = clump_ip->snpfield_search_order;
    uii = 0x40000000;
    do {
      ujj = strlen(bufptr) + 1;
      memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, ujj);
      header_id_map[ulii++] = uii++;
      bufptr = &(bufptr[ujj]);
    } while (*bufptr);
  } else {
    memcpy(sorted_header_dict, "SNP", 4);
    header_id_map[0] = 0x40000000;
    ulii++;
  }
  if (clump_ip->pfield_search_order) {
    bufptr = clump_ip->pfield_search_order;
    uii = 0x20000000;
    do {
      ujj = strlen(bufptr) + 1;
      memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, ujj);
      header_id_map[ulii++] = uii++;
      bufptr = &(bufptr[ujj]);
    } while (*bufptr);
  } else {
    memcpy(&(sorted_header_dict[ulii * max_header_len]), "P", 2);
    header_id_map[ulii++] = 0x20000000;
  }
  if (annot_flattened) {
    bufptr = annot_flattened;
    uii = 2;
    do {
      ujj = strlen(bufptr) + 1;
      memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, ujj);
      header_id_map[ulii++] = uii++;
      bufptr = &(bufptr[ujj]);
    } while (*bufptr);
  }
  if (qsort_ext(sorted_header_dict, header_dict_ct, max_header_len, strcmp_deref, (char*)header_id_map, sizeof(int32_t))) {
    goto clump_reports_ret_NOMEM;
  }
  if (scan_for_duplicate_ids(sorted_header_dict, header_dict_ct, max_header_len)) {
    logerrprint("Error: Duplicate --clump-snp-field/--clump-field/--clump-annotate field name.\n");
    goto clump_reports_ret_INVALID_CMDLINE;
  }

  if (bigstack_calloc_ui(marker_ct, &nsig_arr)) {
    goto clump_reports_ret_NOMEM;
  }
  clump_entries = (Clump_entry**)bigstack_alloc(marker_ct * sizeof(intptr_t));
  if (!clump_entries) {
    goto clump_reports_ret_NOMEM;
  }
  fill_ulong_zero(marker_ct, (uintptr_t*)clump_entries);
  // 3. load file(s) in sequence.  start with array of null pointers, allocate
  //    from bottom of stack (possibly need to save p-val, file number,
  //    annotations, and/or pointer to next entry) while updating
  //    p-val/reverse-lookup array
  bufptr = clump_ip->fnames_flattened;
  do {
    fname_ptr = bufptr;
    bufptr = strchr(bufptr, '\0');
    bufptr++;
    file_ct++;
  } while (*bufptr);
  loadbuft = (char*)g_bigstack_base;
  if (clump_best) {
    if (file_ct == 2) {
      if (!clump_index_first) {
        logerrprint("Error: --clump-best can no longer be used with two --clump files unless\n--clump-index-first is also specified.  (Contact the developers if this is\nproblematic.)\n");
        goto clump_reports_ret_INVALID_CMDLINE;
      }
    } else if (file_ct > 2) {
      logerrprint("Error: --clump-best can no longer be used with more than two --clump files.\n(Contact the developers if this is problematic.)\n");
      goto clump_reports_ret_INVALID_CMDLINE;
    }
    // only draw proxies from this file
    best_fidx_match = file_ct;
  }
  // Suppose the current line has a super-long allele code which must be saved
  // (since it will go into the ANNOT field).  Then the new allocation may need
  // to be the size of the entire line.  So, to be safe, we require the current
  // line to fit in ~half of available workspace.
  // To reduce the risk of 32-bit integer overflow bugs, we cap line length at
  // a bit under 2^30 instead of 2^31 here.
  extra_annot_space = (48 + 2 * annot_ct) & (~(15 * ONELU));
  cur_bigstack_left = bigstack_left();
  if (cur_bigstack_left <= 2 * MAXLINELEN + extra_annot_space) {
    goto clump_reports_ret_NOMEM;
  } else if (cur_bigstack_left - extra_annot_space >= MAXLINEBUFLEN) {
    loadbuft[(MAXLINEBUFLEN / 2) - 1] = ' ';
  }
  if (clump_index_first && (file_ct > 1)) {
    index_eligible = 0;
  }
  // load in reverse order since we're adding to the front of the linked lists
  for (file_idx = file_ct; file_idx; file_idx--) {
    retval = gzopen_read_checked(fname_ptr, &gz_infile);
    if (retval) {
      goto clump_reports_ret_1;
    }
    loadbuft_size = bigstack_left();
    if (loadbuft_size <= 2 * MAXLINELEN + extra_annot_space) {
      goto clump_reports_ret_NOMEM;
    }
    loadbuft_size = (loadbuft_size - extra_annot_space) / 2;
    if (loadbuft_size >= MAXLINEBUFLEN / 2) {
      loadbuft_size = MAXLINEBUFLEN / 2;
      // no space-termination needed
    } else {
      loadbuft[loadbuft_size - 1] = ' ';
    }
    ukk = 0x7fffffff; // highest-precedence variant ID header seen so far
    umm = 0x7fffffff; // highest-precedence p-value header seen so far
    // load_to_first_token() with potentially gzipped input.  Move this to
    // plink_common if anything else needs it.
    line_idx = 0;
    while (1) {
      line_idx++;
      if (!gzgets(gz_infile, loadbuft, loadbuft_size)) {
	if (!gzeof(gz_infile)) {
	  goto clump_reports_ret_READ_FAIL;
	} else {
          LOGPREPRINTFWW("Error: Empty %s.\n", fname_ptr);
	  goto clump_reports_ret_INVALID_FORMAT_2;
	}
      }
      if (!(loadbuft[loadbuft_size - 1])) {
	if (loadbuft_size == MAXLINEBUFLEN / 2) {
	  LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, fname_ptr);
	  goto clump_reports_ret_INVALID_FORMAT_2;
	} else {
	  goto clump_reports_ret_NOMEM;
	}
      }
      bufptr = skip_initial_spaces(loadbuft);
      if (!is_eoln_kns(*bufptr)) {
	break;
      }
    }
    fill_ulong_zero(annot_ct_p2_ctl, col_bitfield);
    uii = 0; // current 0-based column number
    // We don't know in advance when the highest-precedence SNP/p-val columns
    // will appear, so we initially populate parse_table with
    //   [2k]: header type index (0 = variant ID, 1 = p-val, 2+ = annot)
    //   [2k + 1]: 0-based column number
    // and then sort at the end.
    cur_read_ct = 2;
    parse_table[0] = 0;
    parse_table[2] = 1;
    do {
      bufptr2 = token_endnn(bufptr);
      ii = bsearch_str(bufptr, (uintptr_t)(bufptr2 - bufptr), sorted_header_dict, max_header_len, header_dict_ct);
      if (ii != -1) {
	ujj = header_id_map[(uint32_t)ii];
        if (ujj >= 0x40000000) {
          if (ujj < ukk) {
	    // ignore title if higher-precedence title already seen
	    set_bit(0, col_bitfield);
	    ukk = ujj;
	    parse_table[1] = uii; // temporary storage
	  } else if (ujj == ukk) {
	    goto clump_reports_ret_DUPLICATE_HEADER_COL;
	  }
	} else if (ujj >= 0x20000000) {
	  if (ujj < umm) {
	    set_bit(1, col_bitfield);
            umm = ujj;
	    parse_table[3] = uii;
	  } else if (ujj == umm) {
	    goto clump_reports_ret_DUPLICATE_HEADER_COL;
	  }
	} else {
	  if (is_set(col_bitfield, ujj)) {
	    goto clump_reports_ret_DUPLICATE_HEADER_COL;
	  }
	  set_bit(ujj, col_bitfield);
          parse_table[cur_read_ct * 2 + 1] = uii;
	  parse_table[cur_read_ct * 2] = ujj;
	  cur_read_ct++;
	}
      }
      bufptr = skip_initial_spaces(bufptr2);
      uii++;
    } while (!is_eoln_kns(*bufptr));
    if (!is_set(col_bitfield, 0)) {
      LOGPREPRINTFWW("Error: No variant ID field found in %s.\n", fname_ptr);
      goto clump_reports_ret_INVALID_FORMAT_2;
    } else if (!is_set(col_bitfield, 1)) {
      LOGPREPRINTFWW("Error: No p-value field found in %s.\n", fname_ptr);
      goto clump_reports_ret_INVALID_FORMAT_2;
    }
#ifdef __cplusplus
    std::sort((int64_t*)parse_table, (int64_t*)(&(parse_table[cur_read_ct * 2])));
#else
    qsort((int64_t*)parse_table, cur_read_ct, sizeof(int64_t), llcmp);
#endif
    for (uii = cur_read_ct - 1; uii; uii--) {
      parse_table[uii * 2 + 1] -= parse_table[uii * 2 - 1] + 1;
    }
  clump_reports_load_loop:
    while (1) {
      line_idx++;
      if (!gzgets(gz_infile, loadbuft, loadbuft_size)) {
	if (!gzeof(gz_infile)) {
	  goto clump_reports_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuft[loadbuft_size - 1]) {
	if (loadbuft_size == MAXLINEBUFLEN / 2) {
	  LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, fname_ptr);
	  goto clump_reports_ret_INVALID_FORMAT_2;
	}
	goto clump_reports_ret_NOMEM;
      }
      bufptr = skip_initial_spaces(loadbuft);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      fill_uint_zero(annot_ct_p2 * 2, cur_parse_info);
      uii = 0;
      ukk = annot_ct * 2; // annotation string length
      for (; uii < cur_read_ct; uii++) {
	bufptr = next_token_multz(bufptr, parse_table[uii * 2 + 1]);
        if (no_more_tokens_kns(bufptr)) {
	  // PLINK 1.07 --clump just skips the line in this situation, instead
	  // of erroring out, so we replicate that
	  goto clump_reports_load_loop;
	}
	bufptr2 = token_endnn(bufptr);
	ujj = parse_table[uii * 2] * 2;
	cur_parse_info[ujj] = (uintptr_t)(bufptr - loadbuft);
        cur_parse_info[ujj + 1] = (uintptr_t)(bufptr2 - bufptr);
	if (ujj > 2) {
	  ukk += cur_parse_info[ujj + 1];
	}
	bufptr = skip_initial_spaces(bufptr2);
      }
      if (scan_double(&(loadbuft[cur_parse_info[2]]), &pval)) {
	continue;
      }
      if (pval < 0.0) {
	LOGPREPRINTFWW("Error: Negative p-value on line %" PRIuPTR " of %s.\n", line_idx, fname_ptr);
	goto clump_reports_ret_INVALID_FORMAT_2;
      }
      marker_uidx = id_htable_find(&(loadbuft[cur_parse_info[0]]), cur_parse_info[1], marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
      if (marker_uidx == 0xffffffffU) {
	// variant ID not in current fileset
	if ((pval <= p1_thresh) && index_eligible) {
	  // actually a top variant, track it
	  missing_variant_ct++;
	  // screw it, just allocate these outside the workspace
	  uii = cur_parse_info[1];
	  if (uii >= max_missing_id_len) {
	    max_missing_id_len = uii + 1;
	  }
	  cm_ptr = (Clump_missing_id*)malloc(offsetof(Clump_missing_id, idstr) + uii + 1);
	  cm_ptr->pval = pval;
	  cm_ptr->next = not_found_list;
	  not_found_list = cm_ptr;
	  memcpyx(cm_ptr->idstr, &(loadbuft[cur_parse_info[0]]), uii, '\0');
	}
	continue;
      }
      marker_idx = marker_uidx_to_idx[marker_uidx];
      if (pval > load_pthresh) {
	if (pval >= 0.05) {
	  if (pval > 1) {
	    LOGPREPRINTFWW("Error: p-value > 1 on line %" PRIuPTR " of %s.\n", line_idx, fname_ptr);
	    goto clump_reports_ret_INVALID_FORMAT_2;
	  }
	  nsig_arr[marker_idx] += 1;
	}
	continue;
      }
      clump_entry_ptr = (Clump_entry*)bigstack_end_alloc(offsetof(Clump_entry, annot) + ukk - 1);
      if (!clump_entry_ptr) {
	goto clump_reports_ret_NOMEM;
      }
      clump_entry_ptr->pval = pval;
      clump_entry_ptr->next = clump_entries[marker_idx];
      clump_entry_ptr->fidx = file_idx;
      if (annot_ct) {
	bufptr = clump_entry_ptr->annot;
	uii = 2;
	while (1) {
          bufptr = memcpya(bufptr, &(loadbuft[cur_parse_info[uii * 2]]), cur_parse_info[uii * 2 + 1]);
	  if (++uii == annot_ct_p2) {
	    break;
	  }
	  bufptr = memcpya(bufptr, ", ", 2);
	}
	*bufptr = '\0';
      }
      clump_entries[marker_idx] = clump_entry_ptr;
      if ((pval <= p1_thresh) && index_eligible) {
	set_bit(marker_idx, cur_bitfield);
      }
      loadbuft_size = bigstack_left();
      if (loadbuft_size <= 2 * MAXLINELEN + extra_annot_space) {
	goto clump_reports_ret_NOMEM;
      }
      loadbuft_size = (loadbuft_size - extra_annot_space) / 2;
      if (loadbuft_size >= MAXLINEBUFLEN / 2) {
	loadbuft_size = MAXLINEBUFLEN / 2;
	// no space-termination needed
      } else {
	loadbuft[loadbuft_size - 1] = ' ';
      }
    }
    if (gzclose_null(&gz_infile)) {
      goto clump_reports_ret_READ_FAIL;
    }
    if (file_idx > 1) {
      fname_ptr = &(fname_ptr[-3]);
      while (*fname_ptr) {
	fname_ptr--;
      }
      fname_ptr++;
      if (clump_index_first && (file_idx == 2)) {
	index_eligible = 1;
      }
    }
  }
  // 4. sort p-val array, greedily form clumps
  index_ct = popcount_longs(cur_bitfield, marker_ctl);
  if (!index_ct) {
    logerrprint("Warning: No significant --clump results.  Skipping.\n");
    goto clump_reports_ret_1;
  }
  if (bigstack_alloc_d(index_ct, &sorted_pvals) ||
      bigstack_alloc_ui(index_ct, &pval_map)) {
    goto clump_reports_ret_NOMEM;
  }
  marker_idx = 0;
  for (uii = 0; uii < index_ct; uii++, marker_idx++) {
    marker_idx = next_set_unsafe(cur_bitfield, marker_idx);
    clump_entry_ptr = clump_entries[marker_idx];
    pval = clump_entry_ptr->pval;
    if (!clump_index_first) {
      while (clump_entry_ptr->next) {
	clump_entry_ptr = clump_entry_ptr->next;
	if (clump_entry_ptr->pval < pval) {
	  pval = clump_entry_ptr->pval;
	}
      }
    }
    sorted_pvals[uii] = pval;
    pval_map[uii] = marker_idx;
  }
  if (qsort_ext((char*)sorted_pvals, index_ct, sizeof(double), double_cmp_deref, (char*)pval_map, sizeof(int32_t))) {
    goto clump_reports_ret_NOMEM;
  }
  if (bigstack_alloc_ui(marker_ct, &marker_idx_to_uidx) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(5 * founder_ctv2, &index_data)) {
    goto clump_reports_ret_NOMEM;
  }
  for (uii = 1; uii <= 5; uii++) {
    index_data[uii * founder_ctv2 - 2] = 0;
    index_data[uii * founder_ctv2 - 1] = 0;
  }
  if (alloc_collapsed_haploid_filters(founder_info, sex_male, unfiltered_sample_ct, founder_ct, Y_FIX_NEEDED, 1, &founder_include2, &founder_male_include2)) {
    goto clump_reports_ret_NOMEM;
 }
  if (clump_verbose && rg_setdefs) {
    if (bigstack_alloc_ul(BITCT_TO_WORDCT(range_chrom_max), &rangematch_bitfield)) {
      goto clump_reports_ret_NOMEM;
    }
  }
  window_data = (uintptr_t*)g_bigstack_base;
  max_window_size = bigstack_left() / (founder_ctv2 * sizeof(intptr_t) + sizeof(Cur_clump_info));
  if (!max_window_size) {
    goto clump_reports_ret_NOMEM;
  }
  fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_idx_to_uidx);
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  // now this indicates whether a variant has previously been in a clump
  fill_ulong_zero(marker_ctl, cur_bitfield);
  // 5. iterate through clumps, calculate r^2 and write output
  memcpy(outname_end, ".clumped", 9);
  if (fopen_checked(outname, "w", &outfile)) {
    goto clump_reports_ret_OPEN_FAIL;
  }
  bufptr = tbuf2;
  if (clump_verbose) {
    *bufptr++ = '\n';
  }
  bufptr = memcpya(bufptr, " CHR    F ", 10);
  bufptr = fw_strcpyn(plink_maxsnp, 3, "SNP", bufptr);
  // replicate the misaligned non-verbose header for now
  bufptr = memcpya(bufptr, "         BP          ", clump_verbose? 21 : 19);
  bufptr = strcpya(bufptr, "P    TOTAL   NSIG    S05    S01   S001  S0001");
  if (!clump_verbose) {
    bufptr = memcpya(bufptr, "    SP2\n", 8);
    if (fwrite_checked(tbuf2, bufptr - tbuf2, outfile)) {
      goto clump_reports_ret_WRITE_FAIL;
    }
    if (rg_setdefs) {
      memcpy(&(outname_end[8]), ".ranges", 8);
      if (fopen_checked(outname, "w", &outfile_ranges)) {
	goto clump_reports_ret_OPEN_FAIL;
      }
      bufptr = fw_strcpyn(plink_maxsnp, 3, "SNP", &(tbuf2[5]));
      bufptr = strcpya(bufptr, "          P      N                          POS         KB RANGES\n");
      if (fwrite_checked(tbuf2, bufptr - tbuf2, outfile_ranges)) {
	goto clump_reports_ret_WRITE_FAIL;
      }
    }
  } else {
    *bufptr++ = '\n';
    header2_ptr = bufptr;
    header1_len = (uintptr_t)(header2_ptr - tbuf2);
    *bufptr++ = '\n';
    bufptr = memseta(bufptr, 32, 19 + plink_maxsnp);
    bufptr = strcpya(bufptr, "KB      RSQ  ALLELES    F            P ");
    if (annot_flattened) {
      bufptr = memcpya(bufptr, "       ANNOT", 12);
    }
    bufptr = memcpya(bufptr, "\n  (INDEX) ", 11);
    header2_len = (uintptr_t)(bufptr - header2_ptr);
  }
  if (clump_best) {
    memcpy(&(outname_end[8]), ".best", 6);
    if (fopen_checked(outname, "w", &outfile_best)) {
      goto clump_reports_ret_OPEN_FAIL;
    }
    bufptr = fw_strcpyn(plink_maxsnp, 5, "INDEX", g_textbuf);
    *bufptr++ = ' ';
    bufptr = fw_strcpyn(plink_maxsnp, 4, "PSNP", bufptr);
    bufptr = strcpya(bufptr, "    RSQ       KB        P  ALLELES        F\n");
    if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile_best)) {
      goto clump_reports_ret_WRITE_FAIL;
    }
  }
  for (sp_idx = 0; sp_idx < index_ct; sp_idx++) {
    ivar_idx = pval_map[sp_idx];
    if ((!clump_best) && is_set(cur_bitfield, ivar_idx)) {
      continue;
    }
    ivar_uidx = marker_idx_to_uidx[ivar_idx];
    cur_bp = marker_pos[ivar_uidx];
    uii = get_variant_chrom_fo_idx(chrom_info_ptr, ivar_uidx);
    clump_chrom_idx = chrom_info_ptr->chrom_file_order[uii];
    ujj = chrom_info_ptr->chrom_fo_vidx_start[uii];
    if (cur_bp < bp_radius) {
      clump_uidx_first = ujj;
    } else {
      clump_uidx_first = ujj + uint32arr_greater_than(&(marker_pos[ujj]), ivar_uidx + 1 - ujj, cur_bp - bp_radius);
    }
    next_unset_unsafe_ck(marker_exclude, &clump_uidx_first);
    clump_uidx_last = ivar_uidx + uint32arr_greater_than(&(marker_pos[ivar_uidx]), chrom_info_ptr->chrom_fo_vidx_start[uii + 1] - ivar_uidx, cur_bp + bp_radius + 1);
    prev_unset_unsafe_ck(marker_exclude, &clump_uidx_last);
    marker_uidx = clump_uidx_first;
    marker_idx = ivar_idx + popcount_bit_idx(marker_exclude, clump_uidx_first, ivar_uidx) + clump_uidx_first - ivar_uidx;
    // Don't want to seek backwards in the file any more than necessary, so
    // 1. load all clump-inclusion candidates before index variant
    // 2. load index variant, compute pairwise r^2s
    // 3. load one clump-inclusion at a time after index variant, compute r^2
    // 4. write main result
    cur_window_size = 0;
    is_haploid = is_set(haploid_mask, clump_chrom_idx);
    is_x = (clump_chrom_idx == (uint32_t)chrom_info_ptr->xymt_codes[X_OFFSET]);
    is_y = (clump_chrom_idx == (uint32_t)chrom_info_ptr->xymt_codes[Y_OFFSET]);
    window_data_ptr = window_data;
    for (; marker_idx < ivar_idx; marker_uidx++, marker_idx++) {
      next_unset_unsafe_ck(marker_exclude, &marker_uidx);
      if (((!allow_overlap) && is_set(cur_bitfield, marker_idx)) || ((!clump_entries[marker_idx]) && (!nsig_arr[marker_idx]))) {
	continue;
      }
      if (++cur_window_size == max_window_size) {
	goto clump_reports_ret_NOMEM;
      }
      if (fseeko(bedfile, bed_offset + marker_uidx * ((uint64_t)unfiltered_sample_ct4), SEEK_SET)) {
	goto clump_reports_ret_READ_FAIL;
      }
      window_data_ptr[founder_ctv2 - 2] = 0;
      window_data_ptr[founder_ctv2 - 1] = 0;
      if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, is_set(marker_reverse, marker_uidx), bedfile, loadbuf_raw, window_data_ptr)) {
	goto clump_reports_ret_READ_FAIL;
      }
      if (is_haploid) {
	haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)window_data_ptr);
      }
      window_data_ptr = &(window_data_ptr[founder_ctv2]);
    }
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    if (fseeko(bedfile, bed_offset + marker_uidx * ((uint64_t)unfiltered_sample_ct4), SEEK_SET)) {
      goto clump_reports_ret_READ_FAIL;
    }
    window_data_ptr[founder_ctv2 - 2] = 0;
    window_data_ptr[founder_ctv2 - 1] = 0;
    if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, is_set(marker_reverse, marker_uidx), bedfile, loadbuf_raw, window_data_ptr)) {
      goto clump_reports_ret_READ_FAIL;
    }
    if (is_haploid) {
      haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)window_data_ptr);
    }
    vec_datamask(founder_ct, 0, window_data_ptr, founder_include2, index_data);
    index_tots[0] = popcount2_longs(index_data, founder_ctl2);
    vec_datamask(founder_ct, 2, window_data_ptr, founder_include2, &(index_data[founder_ctv2]));
    index_tots[1] = popcount2_longs(&(index_data[founder_ctv2]), founder_ctl2);
    vec_datamask(founder_ct, 3, window_data_ptr, founder_include2, &(index_data[2 * founder_ctv2]));
    index_tots[2] = popcount2_longs(&(index_data[2 * founder_ctv2]), founder_ctl2);
    if (is_x) {
      vec_datamask(founder_ct, 0, window_data_ptr, founder_male_include2, &(index_data[3 * founder_ctv2]));
      index_tots[3] = popcount2_longs(&(index_data[3 * founder_ctv2]), founder_ctl2);
      vec_datamask(founder_ct, 3, window_data_ptr, founder_male_include2, &(index_data[4 * founder_ctv2]));
      index_tots[4] = popcount2_longs(&(index_data[4 * founder_ctv2]), founder_ctl2);
    }
    if (!cur_window_size) {
      cur_clump_base = (Cur_clump_info*)(&(window_data[founder_ctv2]));
    } else {
      cur_clump_base = (Cur_clump_info*)window_data_ptr;
    }
    cc_ptr = cur_clump_base;
    window_data_ptr = window_data;
    marker_uidx = clump_uidx_first;
    marker_idx = ivar_idx + popcount_bit_idx(marker_exclude, clump_uidx_first, ivar_uidx) + clump_uidx_first - ivar_uidx;
    max_r2 = -1;
    max_r2_uidx = 0xffffffffU;
    fill_ulong_zero(5, histo);
    best_entry_ptr = nullptr;
    for (; marker_idx < ivar_idx; marker_uidx++, marker_idx++) {
      marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
      clump_entry_ptr = clump_entries[marker_idx];
      if (((!allow_overlap) && is_set(cur_bitfield, marker_idx)) || ((!clump_entry_ptr) && (!nsig_arr[marker_idx]))) {
	continue;
      }
      genovec_3freq(window_data_ptr, index_data, founder_ctl2, &(counts[0]), &(counts[1]), &(counts[2]));
      counts[0] = index_tots[0] - counts[0] - counts[1] - counts[2];
      genovec_3freq(window_data_ptr, &(index_data[founder_ctv2]), founder_ctl2, &(counts[3]), &(counts[4]), &(counts[5]));
      counts[3] = index_tots[1] - counts[3] - counts[4] - counts[5];
      genovec_3freq(window_data_ptr, &(index_data[2 * founder_ctv2]), founder_ctl2, &(counts[6]), &(counts[7]), &(counts[8]));
      counts[6] = index_tots[2] - counts[6] - counts[7] - counts[8];
      if (is_x) {
        genovec_3freq(window_data_ptr, &(index_data[3 * founder_ctv2]), founder_ctl2, &(counts[9]), &(counts[10]), &(counts[11]));
        counts[9] = index_tots[3] - counts[9] - counts[11];
        genovec_3freq(window_data_ptr, &(index_data[4 * founder_ctv2]), founder_ctl2, &(counts[15]), &(counts[16]), &(counts[17]));
        counts[15] = index_tots[4] - counts[15] - counts[17];
      }
      if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x, &freqx1, &freqx2, &freq11)) {
	freq11_expected = freqx1 * freq1x;
	dxx = freq11 - freq11_expected;
	cur_r2 = fabs(dxx);
        // if r^2 threshold is 0, let everything else through but exclude the
        // apparent zeroes.  Zeroes *are* included if r2_thresh is negative,
	// though (only nans are rejected then).
        if (cur_r2 >= SMALL_EPSILON) {
	  cur_r2 = cur_r2 * dxx / (freq11_expected * freq2x * freqx2);
	} else {
	  cur_r2 = 0;
	}
	if (fabs(cur_r2) > r2_thresh) {
	  while (clump_entry_ptr) {
	    dxx = clump_entry_ptr->pval;
	    update_clump_histo(dxx, histo);
	    if (dxx < p2_thresh) {
	      if (((unsigned char*)cc_ptr) >= g_bigstack_end) {
		goto clump_reports_ret_NOMEM;
	      }
	      cc_ptr->r2 = cur_r2;
	      cc_ptr->marker_idx = marker_idx;
	      uii = clump_entry_ptr->fidx;
	      cc_ptr->fidx = uii;
	      if ((uii == best_fidx_match) && (fabs(cur_r2) > max_r2)) {
		max_r2 = cur_r2;
		max_r2_uidx = marker_uidx;
		best_entry_ptr = clump_entry_ptr;
	      }
	      cc_ptr++;
	    }
	    clump_entry_ptr = clump_entry_ptr->next;
	  }
	  histo[0] += nsig_arr[marker_idx];
	  set_bit(marker_idx, cur_bitfield);
	}
      }
      window_data_ptr = &(window_data_ptr[founder_ctv2]);
    }
    pval = sorted_pvals[sp_idx];
    clump_entry_ptr = clump_entries[ivar_idx];
    uii = 0;
    if (clump_entry_ptr->pval != pval) {
      uii = 1;
      do {
	dxx = clump_entry_ptr->pval;
	update_clump_histo(dxx, histo);
	if (dxx < p2_thresh) {
	  if (((unsigned char*)cc_ptr) >= g_bigstack_end) {
	    goto clump_reports_ret_NOMEM;
	  }
	  cc_ptr->r2 = 1;
	  cc_ptr->marker_idx = ivar_idx;
	  cc_ptr->fidx = clump_entry_ptr->fidx;
	  // clump_best match should be impossible here
	  cc_ptr++;
	}
	clump_entry_ptr = clump_entry_ptr->next;
      } while (clump_entry_ptr->pval != pval);
    }
    index_fidx = clump_entry_ptr->fidx;
    if (annot_flattened) {
      annot_ptr = clump_entry_ptr->annot;
    }
    if ((!clump_best) || allow_overlap || (!is_set(cur_bitfield, ivar_idx))) {
      if (clump_entry_ptr->next) {
	uii = 1;
	do {
	  clump_entry_ptr = clump_entry_ptr->next;
	  dxx = clump_entry_ptr->pval;
	  update_clump_histo(dxx, histo);
	  if (dxx < p2_thresh) {
	    if (((unsigned char*)cc_ptr) >= g_bigstack_end) {
	      goto clump_reports_ret_NOMEM;
	    }
	    cc_ptr->r2 = 1;
	    cc_ptr->marker_idx = ivar_idx;
	    cc_ptr->fidx = clump_entry_ptr->fidx;
	    if (clump_best) {
	      max_r2 = 1;
	      max_r2_uidx = ivar_uidx;
	      best_entry_ptr = clump_entry_ptr;
	    }
	    cc_ptr++;
	  }
	} while (clump_entry_ptr->next);
      }
    }
    // include co-located entries in the clump and mark the position as clumped
    // iff
    //   i. there were co-located entries in the first place, and either
    //     ii-a. overlaps are permitted or
    //     ii-b. index variant position was not previously clumped
    if ((uii || nsig_arr[ivar_idx]) && (allow_overlap || (!is_set(cur_bitfield, ivar_idx)))) {
      histo[0] += nsig_arr[ivar_idx];
      set_bit(ivar_idx, cur_bitfield);
    }
    marker_uidx = ivar_uidx;
    marker_idx = ivar_idx;
    while (marker_uidx < clump_uidx_last) {
      marker_uidx++;
      next_unset_unsafe_ck(marker_exclude, &marker_uidx);
      marker_idx++;
      clump_entry_ptr = clump_entries[marker_idx];
      if (((!allow_overlap) && is_set(cur_bitfield, marker_idx)) || ((!clump_entry_ptr) && (!nsig_arr[marker_idx]))) {
	continue;
      }
      if (fseeko(bedfile, bed_offset + marker_uidx * ((uint64_t)unfiltered_sample_ct4), SEEK_SET)) {
	goto clump_reports_ret_READ_FAIL;
      }
      window_data[founder_ctv2 - 2] = 0;
      window_data[founder_ctv2 - 1] = 0;
      if (load_and_collapse_incl(unfiltered_sample_ct, founder_ct, founder_info, final_mask, is_set(marker_reverse, marker_uidx), bedfile, loadbuf_raw, window_data)) {
	goto clump_reports_ret_READ_FAIL;
      }
      if (is_haploid) {
        haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)window_data);
      }
      genovec_3freq(window_data, index_data, founder_ctl2, &(counts[0]), &(counts[1]), &(counts[2]));
      counts[0] = index_tots[0] - counts[0] - counts[1] - counts[2];
      genovec_3freq(window_data, &(index_data[founder_ctv2]), founder_ctl2, &(counts[3]), &(counts[4]), &(counts[5]));
      counts[3] = index_tots[1] - counts[3] - counts[4] - counts[5];
      genovec_3freq(window_data, &(index_data[2 * founder_ctv2]), founder_ctl2, &(counts[6]), &(counts[7]), &(counts[8]));
      counts[6] = index_tots[2] - counts[6] - counts[7] - counts[8];
      if (is_x) {
        genovec_3freq(window_data, &(index_data[3 * founder_ctv2]), founder_ctl2, &(counts[9]), &(counts[10]), &(counts[11]));
        counts[9] = index_tots[3] - counts[9] - counts[11];
        genovec_3freq(window_data, &(index_data[4 * founder_ctv2]), founder_ctl2, &(counts[15]), &(counts[16]), &(counts[17]));
        counts[15] = index_tots[4] - counts[15] - counts[17];
      }
      if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x, &freqx1, &freqx2, &freq11)) {
        freq11_expected = freqx1 * freq1x;
	dxx = freq11 - freq11_expected;
	cur_r2 = fabs(dxx);
	if (cur_r2 >= SMALL_EPSILON) {
	  cur_r2 = cur_r2 * dxx / (freq11_expected * freq2x * freqx2);
	} else {
	  cur_r2 = 0;
	}
	if (fabs(cur_r2) > r2_thresh) {
	  while (clump_entry_ptr) {
	    dxx = clump_entry_ptr->pval;
            update_clump_histo(dxx, histo);
	    if (dxx < p2_thresh) {
	      if (((unsigned char*)cc_ptr) >= g_bigstack_end) {
		goto clump_reports_ret_NOMEM;
	      }
	      cc_ptr->r2 = cur_r2;
	      cc_ptr->marker_idx = marker_idx;
	      uii = clump_entry_ptr->fidx;
	      cc_ptr->fidx = uii;
	      if ((uii == best_fidx_match) && (fabs(cur_r2) > max_r2)) {
		max_r2 = cur_r2;
		max_r2_uidx = marker_uidx;
		best_entry_ptr = clump_entry_ptr;
	      }
	      cc_ptr++;
	    }
	    clump_entry_ptr = clump_entry_ptr->next;
	  }
	  histo[0] += nsig_arr[marker_idx];
	  set_bit(marker_idx, cur_bitfield);
	}
      }
    }
    cur_window_size = (uintptr_t)(cc_ptr - cur_clump_base);
    if (require_multifile) {
      if (cur_window_size < 2) {
	continue;
      }
      uii = cur_clump_base[0].fidx;
      for (ulii = 1; ulii < cur_window_size; ulii++) {
        if (uii != cur_clump_base[ulii].fidx) {
	  break;
	}
      }
      if (ulii == cur_window_size) {
	continue;
      }
    }
    if (clump_best) {
      bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[ivar_uidx * max_marker_id_len]), g_textbuf);
      *bufptr++ = ' ';
      if (best_entry_ptr) {
	bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[max_r2_uidx * max_marker_id_len]), bufptr);
	*bufptr++ = ' ';
        if (max_r2_uidx == ivar_uidx) {
	  bufptr = memcpya(bufptr, "     *", 6);
	} else {
	  bufptr = dtoa_g_wxp3(fabs(max_r2), 6, bufptr);
	}
	*bufptr++ = ' ';
	bufptr = dtoa_g_wxp3x(((double)((int32_t)(marker_pos[max_r2_uidx] - cur_bp))) * 0.001, 8, ' ', bufptr);
	bufptr = dtoa_g_wxp3x(best_entry_ptr->pval, 8, ' ', bufptr);
	if (max_r2 > 0) {
	  uii = 0;
	} else {
	  uii = 1;
	}
        cur_a1 = marker_allele_ptrs[2 * ivar_uidx];
        cur_a2 = marker_allele_ptrs[2 * ivar_uidx + 1];
        bufptr2 = marker_allele_ptrs[2 * max_r2_uidx + uii];
        bufptr3 = marker_allele_ptrs[2 * max_r2_uidx + 1 - uii];
	bufptr4 = cur_a1;
	for (uii = 3; uii; uii--) {
	  if (!(*(++bufptr4))) {
	    bufptr4 = cur_a2;
	    for (; uii; uii--) {
	      if (!(*(++bufptr4))) {
		bufptr4 = bufptr2;
		for (; uii; uii--) {
		  if (!(*(++bufptr4))) {
		    bufptr4 = bufptr3;
		    for (; uii; uii--) {
		      if (!(*(++bufptr4))) {
			bufptr = memseta(bufptr, 32, uii);
			break;
		      }
		    }
		    break;
		  }
		}
		break;
	      }
	    }
	    break;
	  }
	}
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile_best)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
        fputs(cur_a1, outfile_best);
	fputs(bufptr2, outfile_best);
	putc_unlocked('/', outfile_best);
        fputs(cur_a2, outfile_best);
        fputs(bufptr3, outfile_best);
        g_textbuf[0] = ' ';
        bufptr = uint32toa_w8x(best_fidx_match, ' ', &(g_textbuf[1]));
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile_best)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
	if (annot_flattened) {
          fputs(best_entry_ptr->annot, outfile_best);
	}
        putc_unlocked('\n', outfile_best);
      } else {
	bufptr = fw_strcpyn(plink_maxsnp, 2, "NA", bufptr);
        bufptr = memcpya(bufptr, "     NA       NA       NA       NA       NA \n", 45);
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile_best)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
      }
    }
    bufptr = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, clump_chrom_idx, g_textbuf));
    *bufptr++ = ' ';
    bufptr = uint32toa_w4(index_fidx, bufptr);
    *bufptr++ = ' ';
    bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[ivar_uidx * max_marker_id_len]), bufptr);
    *bufptr++ = ' ';
    bufptr = uint32toa_w10x(cur_bp, ' ', bufptr);
    bufptr = dtoa_g_wxp3x(pval, 10, ' ', bufptr);
#ifdef __LP64__
    // may as well be paranoid
    bufptr = width_force(8, bufptr, int64toa((int64_t)(histo[0] + histo[1] + histo[2] + histo[3] + histo[4]), bufptr));
    *bufptr++ = ' ';
    for (uii = 0; uii < 5; uii++) {
      bufptr = width_force(6, bufptr, int64toa((int64_t)((uintptr_t)histo[uii]), bufptr));
      *bufptr++ = ' ';
    }
#else
    bufptr = uint32toa_w8x(histo[0] + histo[1] + histo[2] + histo[3] + histo[4], ' ', bufptr);
    for (uii = 0; uii < 5; uii++) {
      bufptr = uint32toa_w6x(histo[uii], ' ', bufptr);
    }
#endif
    final_clump_ct++;
    min_bp = cur_bp;
    max_bp = cur_bp;
    if (cur_window_size) {
      marker_idx = cur_clump_base[0].marker_idx;
      if (marker_idx < ivar_idx) {
	min_bp = marker_pos[marker_idx_to_uidx[marker_idx]];
      }
      marker_idx = cur_clump_base[cur_window_size - 1].marker_idx;
      if (marker_idx > ivar_idx) {
	max_bp = marker_pos[marker_idx_to_uidx[marker_idx]];
      }
    }
    if (rg_setdefs) {
      ulii = rg_chrom_bounds[clump_chrom_idx];
      cur_rg_setdefs = &(rg_setdefs[ulii]);
      cur_rg_names = &(range_group_names[ulii * max_range_group_id_len + 4]);
      cur_rg_ct = rg_chrom_bounds[clump_chrom_idx + 1] - ulii;
    }
    if (!clump_verbose) {
      if (!cur_window_size) {
	bufptr = memcpya(bufptr, "NONE\n", 5);
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
      } else {
	// avoid buffer overflow
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
	g_textbuf[0] = '(';
	for (ulii = 0; ulii < cur_window_size;) {
          fputs(&(marker_ids[marker_idx_to_uidx[cur_clump_base[ulii].marker_idx] * max_marker_id_len]), outfile);
	  bufptr = uint32toa_x(cur_clump_base[ulii].fidx, ')', &(g_textbuf[1]));
	  ulii++;
	  if (ulii != cur_window_size) {
	    *bufptr++ = ',';
	  }
	  fwrite(g_textbuf, 1, (uintptr_t)(bufptr - g_textbuf), outfile);
	}
	if (putc_checked('\n', outfile)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
      }
      if (rg_setdefs) {
        bufptr = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, clump_chrom_idx, g_textbuf));
        *bufptr++ = ' ';
        bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[ivar_uidx * max_marker_id_len]), bufptr);
        *bufptr++ = ' ';
        bufptr = dtoa_g_wxp4x(pval, 10, ' ', bufptr);
        bufptr = uint32toa_w6x(cur_window_size + 1, ' ', bufptr);
	if (clump_chrom_idx <= chrom_info_ptr->max_code) {
	  bufptr2 = memcpyl3a(bufptr, "chr");
	  bufptr2 = uint32toa(clump_chrom_idx, bufptr2);
	} else if (chrom_info_ptr->zero_extra_chroms) {
	  bufptr2 = memcpya(bufptr, "chr0", 4);
	} else {
	  bufptr2 = strcpya(bufptr, chrom_info_ptr->nonstd_names[clump_chrom_idx]);
	}
        *bufptr2++ = ':';
        bufptr2 = uint32toa(min_bp, bufptr2);
        bufptr2 = memcpya(bufptr2, "..", 2);
        bufptr2 = uint32toa(max_bp, bufptr2);
        bufptr = width_force(28, bufptr, bufptr2);
        *bufptr++ = ' ';
        bufptr = width_force(10, bufptr, dtoa_g(((int32_t)(max_bp - min_bp + 1)) * 0.001, bufptr));
	bufptr = memcpya(bufptr, " [", 2);
        if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile_ranges)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
	uljj = 0;
	for (ulii = 0; ulii < cur_rg_ct; ulii++) {
	  if (interval_in_setdef(cur_rg_setdefs[ulii], min_bp, max_bp)) {
            if (uljj) {
	      putc_unlocked(',', outfile_ranges);
	    } else {
	      uljj = 1;
	    }
            fputs(&(cur_rg_names[ulii * max_range_group_id_len]), outfile_ranges);
	  }
	}
	fputs("]\n", outfile_ranges);
      }
    } else {
      if (fwrite_checked(tbuf2, header1_len, outfile)) {
	goto clump_reports_ret_WRITE_FAIL;
      }
      *bufptr++ = '\n';
      if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	goto clump_reports_ret_WRITE_FAIL;
      }
      if (cur_window_size) {
	if (fwrite_checked(header2_ptr, header2_len, outfile)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
        bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[ivar_uidx * max_marker_id_len]), g_textbuf);
	bufptr = memcpya(bufptr, "          0    1.000 ", 21);
	cur_a1 = marker_allele_ptrs[2 * ivar_uidx];
	a1_len = strlen(cur_a1);
	if (a1_len < 8) {
	  bufptr = memseta(bufptr, 32, 8 - a1_len);
	}
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
	fwrite(cur_a1, 1, a1_len, outfile);
	cur_a2 = marker_allele_ptrs[2 * ivar_uidx + 1];
        a2_len = strlen(cur_a2);
	if (a1_len + a2_len < 5) {
	  allele_padding = 5 - a1_len - a2_len;
	} else {
	  allele_padding = 0;
	}
	g_textbuf[0] = ' ';
        bufptr = uint32toa_w4x(index_fidx, ' ', &(g_textbuf[1]));
	bufptr = dtoa_g_wxp3x(pval, 12, ' ', bufptr);
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
	if (annot_flattened) {
	  bufptr2 = annot_ptr;
	  for (uii = 11; uii; uii--) {
            if (!(*(++bufptr2))) {
	      fwrite("           ", 1, uii, outfile);
	      break;
	    }
	  }
	  fputs(annot_ptr, outfile);
	}
	fputs("\n\n", outfile);
	last_marker_idx = ~ZEROLU;
	if (rg_setdefs) {
	  fill_ulong_zero(BITCT_TO_WORDCT(cur_rg_ct), rangematch_bitfield);
	  unmatched_group_ct = cur_rg_ct;
	}
	for (ulii = 0; ulii < cur_window_size; ulii++) {
	  bufptr = memseta(g_textbuf, 32, 10);
	  marker_idx = cur_clump_base[ulii].marker_idx;
	  if (last_marker_idx != marker_idx) {
	    marker_uidx = marker_idx_to_uidx[marker_idx];
	    clump_entry_ptr = clump_entries[marker_idx];
	    if (rg_setdefs) {
	      uii = marker_pos[marker_uidx];
	      uljj = 0; // range group idx
	      ulkk = 0; // number of new matches
	      for (ulmm = 0; ulmm < unmatched_group_ct; uljj++, ulmm++) {
		next_unset_ul_unsafe_ck(rangematch_bitfield, &uljj);
		if (interval_in_setdef(cur_rg_setdefs[uljj], uii, uii + 1)) {
		  set_bit(uljj, rangematch_bitfield);
		  ulkk++;
		}
	      }
	      unmatched_group_ct -= ulkk;
	    }
	  }
	  ukk = cur_clump_base[ulii].fidx;
	  while (clump_entry_ptr->fidx != ukk) {
	    clump_entry_ptr = clump_entry_ptr->next;
	  }
	  bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
	  *bufptr++ = ' ';
	  bufptr = dtoa_g_wxp3x(((double)(((int32_t)marker_pos[marker_uidx]) - ((int32_t)cur_bp))) * 0.001, 10, ' ', bufptr);
	  cur_r2 = cur_clump_base[ulii].r2;
	  if (cur_r2 > 0) {
	    ujj = 0;
	  } else {
	    ujj = 1; // reversed phase
	  }
	  bufptr = dtoa_g_wxp3x(fabs(cur_r2), 8, ' ', bufptr);
	  bufptr2 = marker_allele_ptrs[marker_uidx * 2 + ujj];
	  bufptr3 = marker_allele_ptrs[marker_uidx * 2 + 1 - ujj];
	  if (allele_padding) {
	    bufptr4 = bufptr2;
	    for (uii = allele_padding; uii; uii--) {
	      // fast in common case, don't bother to compute strlen for long
	      // indels
	      if (!(*(++bufptr4))) {
		bufptr4 = bufptr3;
		for (; uii; uii--) {
		  if (!(*(++bufptr4))) {
		    bufptr = memseta(bufptr, 32, uii);
		    break;
		  }
		}
		break;
	      }
	    }
	  }
	  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	    goto clump_reports_ret_WRITE_FAIL;
	  }
	  fwrite(cur_a1, 1, a1_len, outfile);
	  fputs(bufptr2, outfile);
	  putc_unlocked('/', outfile);
	  fwrite(cur_a2, 1, a2_len, outfile);
	  fputs(bufptr3, outfile);
	  g_textbuf[0] = ' ';
	  bufptr = uint32toa_w4x(cur_clump_base[ulii].fidx, ' ', &(g_textbuf[1]));
	  bufptr = dtoa_g_wxp3x(clump_entry_ptr->pval, 12, ' ', bufptr);
	  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	    goto clump_reports_ret_WRITE_FAIL;
	  }
	  if (annot_flattened) {
	    bufptr2 = clump_entry_ptr->annot;
	    bufptr3 = bufptr2;
	    for (uii = 11; uii; uii--) {
	      if (!(*(++bufptr3))) {
		fwrite("           ", 1, uii, outfile);
		break;
	      }
	    }
	    fputs(bufptr2, outfile);
	  }
	  putc_unlocked('\n', outfile);
	  last_marker_idx = marker_idx;
	}
	bufptr = memcpya(g_textbuf, "\n          RANGE: ", 18);
	if (clump_chrom_idx <= chrom_info_ptr->max_code) {
	  bufptr = memcpyl3a(bufptr, "chr");
	  bufptr = uint32toa(clump_chrom_idx, bufptr);
	} else if (chrom_info_ptr->zero_extra_chroms) {
	  bufptr = memcpya(bufptr, "chr0", 4);
	} else {
	  bufptr = strcpya(bufptr, chrom_info_ptr->nonstd_names[clump_chrom_idx]);
	}
	*bufptr++ = ':';
	bufptr = uint32toa(min_bp, bufptr);
	bufptr = memcpya(bufptr, "..", 2);
	bufptr = uint32toa(max_bp, bufptr);
	bufptr = memcpya(bufptr, "\n           SPAN: ", 18);
	bufptr = uint32toa((max_bp - min_bp + 1) / 1000, bufptr);
	bufptr = memcpyl3a(bufptr, "kb\n");
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	  goto clump_reports_ret_WRITE_FAIL;
	}
	if (rg_setdefs) {
	  fputs("     GENES w/SNPs: ", outfile);
	  ulii = 0;
	  uljj = 0;
	  unmatched_group_ct = cur_rg_ct - unmatched_group_ct;
	  if (unmatched_group_ct) {
	    while (1) {
	      uljj = next_set_ul_unsafe(rangematch_bitfield, uljj);
	      fputs(&(cur_rg_names[uljj * max_range_group_id_len]), outfile);
	      if (!(--unmatched_group_ct)) {
		break;
	      }
	      uljj++;
	      putc_unlocked(',', outfile);
	    }
	  }
	  putc_unlocked('\n', outfile);
	}
      }
      if (rg_setdefs) {
	if (!cur_window_size) {
	  putc_unlocked('\n', outfile);
	}
	fputs("            GENES: ", outfile);
	uljj = 0;
	for (ulii = 0; ulii < cur_rg_ct; ulii++) {
	  if (interval_in_setdef(cur_rg_setdefs[ulii], min_bp, max_bp)) {
            if (uljj) {
	      if (uljj & 7) {
		putc_unlocked(',', outfile);
	      } else {
		putc_unlocked('\n', outfile);
	      }
	    }
            fputs(&(cur_rg_names[ulii * max_range_group_id_len]), outfile);
	    uljj++;
	  }
	}
	putc_unlocked('\n', outfile);
      }
      if (fwrite_checked("\n------------------------------------------------------------------\n\n", 69, outfile)) {
	goto clump_reports_ret_WRITE_FAIL;
      }
    }
  }
  putc_unlocked('\n', outfile);
  if (missing_variant_ct) {
    // 1. sort by ID (could switch this to hash table-based too)
    // 2. pick smallest pval when duplicates present
    // 3. sort by pval
    // 4. write results
    bigstack_double_reset(bigstack_mark, bigstack_end_mark);
    if (bigstack_alloc_c(missing_variant_ct * max_missing_id_len, &sorted_missing_variant_ids) ||
	bigstack_alloc_d(missing_variant_ct, &sorted_pvals)) {
      goto clump_reports_ret_NOMEM;
    }
    for (ulii = 0; ulii < missing_variant_ct; ulii++) {
      cm_ptr = not_found_list;
      strcpy(&(sorted_missing_variant_ids[ulii * max_missing_id_len]), cm_ptr->idstr);
      sorted_pvals[ulii] = cm_ptr->pval;
      not_found_list = not_found_list->next;
      free(cm_ptr);
    }
    if (qsort_ext(sorted_missing_variant_ids, missing_variant_ct, max_missing_id_len, strcmp_deref, (char*)sorted_pvals, sizeof(double))) {
      goto clump_reports_ret_NOMEM;
    }
    bufptr = sorted_missing_variant_ids;
    uii = strlen(sorted_missing_variant_ids);
    for (ulii = 1; ulii < missing_variant_ct; ulii++) {
      bufptr2 = &(bufptr[max_missing_id_len]);
      ujj = strlen(bufptr2);
      if ((uii == ujj) && (!memcmp(bufptr, bufptr2, uii))) {
	uljj = ulii - 1; // write index
	pval = sorted_pvals[uljj];
	if (pval > sorted_pvals[ulii]) {
	  pval = sorted_pvals[ulii];
	}
        while (++ulii < missing_variant_ct) {
	  bufptr2 = &(bufptr2[max_missing_id_len]);
	  ujj = strlen(bufptr2);
	  if ((uii == ujj) && (!memcmp(bufptr, bufptr2, uii))) {
	    if (pval > sorted_pvals[ulii]) {
	      pval = sorted_pvals[ulii];
	    }
	  } else {
	    sorted_pvals[uljj++] = pval;
	    bufptr = &(bufptr[max_missing_id_len]);
	    memcpy(bufptr, bufptr2, ujj + 1);
	    pval = sorted_pvals[ulii];
	    uii = ujj;
	  }
	}
	sorted_pvals[uljj] = pval;
	ulii = uljj + 1; // save final array length
	break;
      }
      bufptr = bufptr2;
      uii = ujj;
    }
    missing_variant_ct = ulii;
    if (qsort_ext((char*)sorted_pvals, missing_variant_ct, sizeof(double), double_cmp_deref, sorted_missing_variant_ids, max_missing_id_len)) {
      goto clump_reports_ret_NOMEM;
    }
    if (clump_verbose) {
      for (ulii = 0; ulii < missing_variant_ct; ulii++) {
	fputs(&(sorted_missing_variant_ids[ulii * max_missing_id_len]), outfile);
	fputs(" not found in dataset\n", outfile);
      }
      LOGPRINTF("%" PRIuPTR " top variant ID%s missing; see the end of the .clumped file.\n", missing_variant_ct, (missing_variant_ct == 1)? "" : "s");
    } else {
      uljj = MINV(missing_variant_ct, 3);
      for (ulii = 0; ulii < uljj; ulii++) {
	LOGERRPRINTFWW("Warning: '%s' is missing from the main dataset, and is a top variant.\n", &(sorted_missing_variant_ids[ulii * max_missing_id_len]));
      }
      if (missing_variant_ct > 3) {
        fprintf(stderr, "%" PRIuPTR " more top variant ID%s missing; see log file.\n", missing_variant_ct - 3, (missing_variant_ct == 4)? "" : "s");
	for (ulii = 3; ulii < missing_variant_ct; ulii++) {
	  LOGPREPRINTFWW("Warning: '%s' is missing from the main dataset, and is a top variant.\n", &(sorted_missing_variant_ids[ulii * max_missing_id_len]));
	  logstr(g_logbuf);
	}
      }
    }
  }
  putc_unlocked('\n', outfile);
  if (fclose_null(&outfile)) {
    goto clump_reports_ret_WRITE_FAIL;
  }
  outname_end[8] = '\0';
  LOGPRINTF("--clump: %u clump%s formed from %u top variant%s.\n", final_clump_ct, (final_clump_ct == 1)? "" : "s", index_ct, (index_ct == 1)? "" : "s");
  LOGPRINTFWW("Results written to %s .\n", outname);
  if (rg_setdefs && (!clump_verbose)) {
    memcpy(&(outname_end[8]), ".ranges", 8);
    LOGPRINTFWW("--clump-range: Clump/region overlaps reported in %s .\n", outname);
  }
  if (clump_best) {
    memcpy(&(outname_end[8]), ".best", 6);
    LOGPRINTFWW("--clump-best: Best proxies written to %s .\n", outname);
  }
  while (0) {
  clump_reports_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  clump_reports_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  clump_reports_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  clump_reports_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  clump_reports_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  clump_reports_ret_DUPLICATE_HEADER_COL:
    *bufptr2 = '\0';
    LOGPREPRINTFWW("Error: Duplicate column header '%s' in %s.\n", bufptr, fname_ptr);
  clump_reports_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 clump_reports_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  gzclose_cond(gz_infile);
  fclose_cond(outfile);
  fclose_cond(outfile_ranges);
  fclose_cond(outfile_best);
  while (not_found_list) {
    cm_ptr = not_found_list;
    not_found_list = not_found_list->next;
    free(cm_ptr);
  }
  return retval;
}
