// This file is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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

#ifdef __cplusplus
#include <functional> // std::greater
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

void init_ld(ld_info_t* ldip) {
  ldip->prune_modifier = kfLdPrune0;
  ldip->prune_window_size = 0;
  ldip->prune_window_incr = 0;
  ldip->prune_last_param = 0.0;
}

void cleanup_ld(__attribute__((unused)) ld_info_t* ldip) {
}


static inline void popcount_vecs_2intersect(const vul_t* __restrict vvec1_iter, const vul_t* __restrict vvec2a_iter, const vul_t* __restrict vvec2b_iter, uintptr_t vec_ct, uint32_t* popcount_1_2a_ptr, uint32_t* popcount_1_2b_ptr) {
  // popcounts (vvec1 AND vvec2a[0..(ct-1)]) as well as (vvec1 AND vvec2b).  ct
  // is a multiple of 3.
  assert(!(vec_ct % 3));
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t m2 = VCONST_UL(kMask3333);
  const vul_t m4 = VCONST_UL(kMask0F0F);

  // todo: check if moving this right before usage is better, looks like we
  // barely have enough registers...
  const vul_t m8 = VCONST_UL(kMask00FF);
  uint32_t popcount_1_2a = 0;
  uint32_t popcount_1_2b = 0;

  while (1) {
    univec_t acc_a;
    univec_t acc_b;
    acc_a.vi = vul_setzero();
    acc_b.vi = vul_setzero();

    const vul_t* vvec1_stop;
    if (vec_ct < 30) {
      if (!vec_ct) {
	*popcount_1_2a_ptr = popcount_1_2a;
	*popcount_1_2b_ptr = popcount_1_2b;
	return;
      }
      vvec1_stop = &(vvec1_iter[vec_ct]);
      vec_ct = 0;
    } else {
      vvec1_stop = &(vvec1_iter[30]);
      vec_ct -= 30;
    }
    do {
      vul_t loader = *vvec1_iter++;
      vul_t count1a = loader & (*vvec2a_iter++);
      vul_t count1b = loader & (*vvec2b_iter++);
      loader = *vvec1_iter++;
      vul_t count2a = loader & (*vvec2a_iter++);
      vul_t count2b = loader & (*vvec2b_iter++);
      loader = *vvec1_iter++;
      vul_t half1a = loader & (*vvec2a_iter++);
      vul_t half1b = loader & (*vvec2b_iter++);
      const vul_t half2a = vul_rshift(half1a, 1) & m1;
      const vul_t half2b = vul_rshift(half1b, 1) & m1;
      half1a = half1a & m1;
      half1b = half1b & m1;
      count1a = count1a - (vul_rshift(count1a, 1) & m1);
      count1b = count1b - (vul_rshift(count1b, 1) & m1);
      count2a = count2a - (vul_rshift(count2a, 1) & m1);
      count2b = count2b - (vul_rshift(count2b, 1) & m1);
      count1a = count1a + half1a;
      count1b = count1b + half1b;
      count2a = count2a + half2a;
      count2b = count2b + half2b;
      count1a = (count1a & m2) + (vul_rshift(count1a, 2) & m2);
      count1b = (count1b & m2) + (vul_rshift(count1b, 2) & m2);
      count1a = count1a + (count2a & m2) + (vul_rshift(count2a, 2) & m2);
      count1b = count1b + (count2b & m2) + (vul_rshift(count2b, 2) & m2);
      acc_a.vi = acc_a.vi + (count1a & m4) + (vul_rshift(count1a, 4) & m4);
      acc_b.vi = acc_b.vi + (count1b & m4) + (vul_rshift(count1b, 4) & m4);
    } while (vvec1_iter < vvec1_stop);
    acc_a.vi = (acc_a.vi & m8) + (vul_rshift(acc_a.vi, 8) & m8);
    acc_b.vi = (acc_b.vi & m8) + (vul_rshift(acc_b.vi, 8) & m8);
    popcount_1_2a += univec_hsum_16bit(acc_a);
    popcount_1_2b += univec_hsum_16bit(acc_b);
  }
}

// don't bother with popcount_vecs_3intersect for now, but test later

void popcount_longs_2intersect(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2a_iter, const uintptr_t* __restrict bitvec2b_iter, uintptr_t word_ct, uint32_t* popcount_1_2a_ptr, uint32_t* popcount_1_2b_ptr) {
  const uintptr_t* bitvec1_end = &(bitvec1_iter[word_ct]);
  uintptr_t trivec_ct = word_ct / (3 * kWordsPerVec);
  uint32_t popcount_1_2a;
  uint32_t popcount_1_2b;
  popcount_vecs_2intersect((const vul_t*)bitvec1_iter, (const vul_t*)bitvec2a_iter, (const vul_t*)bitvec2b_iter, trivec_ct * 3, &popcount_1_2a, &popcount_1_2b);
  bitvec1_iter = &(bitvec1_iter[trivec_ct * (3 * kWordsPerVec)]);
  bitvec2a_iter = &(bitvec2a_iter[trivec_ct * (3 * kWordsPerVec)]);
  bitvec2b_iter = &(bitvec2b_iter[trivec_ct * (3 * kWordsPerVec)]);
  while (bitvec1_iter < bitvec1_end) {
    const uintptr_t loader1 = *bitvec1_iter++;
    popcount_1_2a += popcount_long(loader1 & (*bitvec2a_iter++));
    popcount_1_2b += popcount_long(loader1 & (*bitvec2b_iter++));
  }
  *popcount_1_2a_ptr = popcount_1_2a;
  *popcount_1_2b_ptr = popcount_1_2b;
}


static inline int32_t dotprod_vecs(const vul_t* __restrict vvec1a_iter, const vul_t* __restrict vvec1b_iter, const vul_t* __restrict vvec2a_iter, const vul_t* __restrict vvec2b_iter, uintptr_t vec_ct) {
  // assumes vvec1a/vvec2a represesent +1s, vvec1b/vvec2b represent -1s, and
  // everything else is 0.  computes
  //   popcount(vvec1a & vvec2a) + popcount(vvec1b & vvec2b)
  //   - popcount(vvec1a & vvec2b) - popcount(vvec1b & vvec2a).
  // ct must be a multiple of 3.
  assert(!(vec_ct % 3));
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t m2 = VCONST_UL(kMask3333);
  const vul_t m4 = VCONST_UL(kMask0F0F);
  int32_t tot = 0;
  while (1) {
    univec_t acc_plus;
    univec_t acc_minus;
    acc_plus.vi = vul_setzero();
    acc_minus.vi = vul_setzero();

    const vul_t* vvec1a_stop;
    if (vec_ct < 30) {
      if (!vec_ct) {
	return tot;
      }
      vvec1a_stop = &(vvec1a_iter[vec_ct]);
      vec_ct = 0;
    } else {
      vvec1a_stop = &(vvec1a_iter[30]);
      vec_ct -= 30;
    }
    do {
      vul_t loader1a = *vvec1a_iter++;
      vul_t loader1b = *vvec1b_iter++;
      vul_t loader2a = *vvec2a_iter++;
      vul_t loader2b = *vvec2b_iter++;
      // loader1a and loader1b are disjoint, etc.; take advantage of that
      vul_t count1_plus = (loader1a & loader2a) | (loader1b & loader2b);
      vul_t count1_minus = (loader1a & loader2b) | (loader1b & loader2a);

      loader1a = *vvec1a_iter++;
      loader1b = *vvec1b_iter++;
      loader2a = *vvec2a_iter++;
      loader2b = *vvec2b_iter++;
      vul_t count2_plus = (loader1a & loader2a) | (loader1b & loader2b);
      vul_t count2_minus = (loader1a & loader2b) | (loader1b & loader2a);

      loader1a = *vvec1a_iter++;
      loader1b = *vvec1b_iter++;
      loader2a = *vvec2a_iter++;
      loader2b = *vvec2b_iter++;
      vul_t half1_plus = (loader1a & loader2a) | (loader1b & loader2b);
      vul_t half1_minus = (loader1a & loader2b) | (loader1b & loader2a);
      const vul_t half2_plus = vul_rshift(half1_plus, 1) & m1;
      const vul_t half2_minus = vul_rshift(half1_minus, 1) & m1;
      half1_plus = half1_plus & m1;
      half1_minus = half1_minus & m1;
      count1_plus = count1_plus - (vul_rshift(count1_plus, 1) & m1);
      count1_minus = count1_minus - (vul_rshift(count1_minus, 1) & m1);
      count2_plus = count2_plus - (vul_rshift(count2_plus, 1) & m1);
      count2_minus = count2_minus - (vul_rshift(count2_minus, 1) & m1);
      count1_plus = count1_plus + half1_plus;
      count1_minus = count1_minus + half1_minus;
      count2_plus = count2_plus + half2_plus;
      count2_minus = count2_minus + half2_minus;
      count1_plus = (count1_plus & m2) + (vul_rshift(count1_plus, 2) & m2);
      count1_minus = (count1_minus & m2) + (vul_rshift(count1_minus, 2) & m2);
      count1_plus = count1_plus + (count2_plus & m2) + (vul_rshift(count2_plus, 2) & m2);
      count1_minus = count1_minus + (count2_minus & m2) + (vul_rshift(count2_minus, 2) & m2);
      acc_plus.vi = acc_plus.vi + (count1_plus & m4) + (vul_rshift(count1_plus, 4) & m4);
      acc_minus.vi = acc_minus.vi + (count1_minus & m4) + (vul_rshift(count1_minus, 4) & m4);
    } while (vvec1a_iter < vvec1a_stop);
    const vul_t m8 = VCONST_UL(kMask00FF);
    acc_plus.vi = (acc_plus.vi & m8) + (vul_rshift(acc_plus.vi, 8) & m8);
    acc_minus.vi = (acc_minus.vi & m8) + (vul_rshift(acc_minus.vi, 8) & m8);
    tot += (uint32_t)univec_hsum_16bit(acc_plus);
    tot -= (uint32_t)univec_hsum_16bit(acc_minus);
  }
}

int32_t dotprod_longs(const uintptr_t* __restrict bitvec1a_iter, const uintptr_t* __restrict bitvec1b_iter, const uintptr_t* __restrict bitvec2a_iter, const uintptr_t* __restrict bitvec2b_iter, uintptr_t word_ct) {
  const uintptr_t* bitvec1a_end = &(bitvec1a_iter[word_ct]);
  uintptr_t trivec_ct = word_ct / (kWordsPerVec * 3);
  int32_t tot = dotprod_vecs((const vul_t*)bitvec1a_iter, (const vul_t*)bitvec1b_iter, (const vul_t*)bitvec2a_iter, (const vul_t*)bitvec2b_iter, trivec_ct * 3);
  bitvec1a_iter = &(bitvec1a_iter[trivec_ct * (3 * kWordsPerVec)]);
  bitvec1b_iter = &(bitvec1b_iter[trivec_ct * (3 * kWordsPerVec)]);
  bitvec2a_iter = &(bitvec2a_iter[trivec_ct * (3 * kWordsPerVec)]);
  bitvec2b_iter = &(bitvec2b_iter[trivec_ct * (3 * kWordsPerVec)]);
  while (bitvec1a_iter < bitvec1a_end) {
    uintptr_t loader1a = *bitvec1a_iter++;
    uintptr_t loader1b = *bitvec1b_iter++;
    uintptr_t loader2a = *bitvec2a_iter++;
    uintptr_t loader2b = *bitvec2b_iter++;
    tot += popcount_long((loader1a & loader2a) | (loader1b & loader2b));
    tot -= popcount_long((loader1a & loader2b) | (loader1b & loader2a));
  }
  return tot;
}

void ldprune_next_subcontig(const uintptr_t* variant_include, const uint32_t* variant_bps, const uint32_t* subcontig_info, const uint32_t* subcontig_thread_assignments, uint32_t x_start, uint32_t x_len, uint32_t y_start, uint32_t y_len, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t prune_window_size, uint32_t thread_idx, uint32_t* subcontig_idx_ptr, uint32_t* subcontig_end_tvidx_ptr, uint32_t* next_window_end_tvidx_ptr, uint32_t* is_x_ptr, uint32_t* is_y_ptr, uint32_t* cur_founder_ct_ptr, uint32_t* cur_founder_ctaw_ptr, uint32_t* cur_founder_ctl_ptr, uintptr_t* entire_variant_buf_word_ct_ptr, uint32_t* variant_uidx_winstart_ptr, uint32_t* variant_uidx_winend_ptr) {
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
    uint32_t variant_uidx_winend = variant_uidx_winstart;
    uint32_t first_window_len = 1;
    do {
      ++variant_uidx_winend;
      next_set_unsafe_ck(variant_include, &variant_uidx_winend);
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
    const uint32_t cur_founder_ctaw = BITCT_TO_ALIGNED_WORDCT(cur_founder_ct);
    *cur_founder_ct_ptr = cur_founder_ct;
    *cur_founder_ctaw_ptr = cur_founder_ctaw;
    *cur_founder_ctl_ptr = BITCT_TO_WORDCT(cur_founder_ct);
    *entire_variant_buf_word_ct_ptr = 3 * cur_founder_ctaw;
    if (is_x) {
      *entire_variant_buf_word_ct_ptr += 3 * BITCT_TO_ALIGNED_WORDCT(founder_ct - founder_male_ct);
    }
  }
}

void genoarr_split_02nm(const uintptr_t* __restrict genoarr, uint32_t sample_ct, uintptr_t* __restrict zero_bitarr, uintptr_t* __restrict two_bitarr, uintptr_t* __restrict nm_bitarr) {
  // ok if trailing bits of genoarr are not zeroed out
  // trailing bits of {zero,two,nm}_bitarr are zeroed out
  const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  halfword_t* zero_bitarr_alias = (halfword_t*)zero_bitarr;
  halfword_t* two_bitarr_alias = (halfword_t*)two_bitarr;
  halfword_t* nm_bitarr_alias = (halfword_t*)nm_bitarr;
  for (uint32_t widx = 0; widx < sample_ctl2; ++widx) {
    const uintptr_t cur_geno_word = genoarr[widx];
    const uint32_t low_halfword = pack_word_to_halfword(cur_geno_word & kMask5555);
    const uint32_t high_halfword = pack_word_to_halfword((cur_geno_word >> 1) & kMask5555);
    zero_bitarr_alias[widx] = ~(low_halfword | high_halfword);
    two_bitarr_alias[widx] = high_halfword & (~low_halfword);
    nm_bitarr_alias[widx] = ~(low_halfword & high_halfword);
  }
  const uint32_t sample_ct_rem = sample_ct % kBitsPerWord;
  if (sample_ct_rem) {
    const uintptr_t trailing_mask = (~k0LU) >> (kBitsPerWord - sample_ct_rem);
    // note that we don't use the halfword aliases here
    const uint32_t last_write_word_idx = sample_ct / kBitsPerWord;
    zero_bitarr[last_write_word_idx] &= trailing_mask;
    two_bitarr[last_write_word_idx] &= trailing_mask;
    nm_bitarr[last_write_word_idx] &= trailing_mask;
  }
}

void ldprune_next_window(const uintptr_t* __restrict variant_include, const uint32_t* __restrict variant_bps, const uint32_t* __restrict tvidxs, const uintptr_t* __restrict cur_window_removed, uint32_t prune_window_size, uint32_t window_incr, uint32_t window_maxl, uint32_t subcontig_end_tvidx, uint32_t* cur_window_size_ptr, uint32_t* __restrict window_start_tvidx_ptr, uint32_t* __restrict variant_uidx_winstart_ptr, uint32_t* __restrict next_window_end_tvidx_ptr, uint32_t* __restrict variant_uidx_winend_ptr, uintptr_t* __restrict occupied_window_slots, uint32_t* winpos_to_slot_idx) {
  uint32_t next_window_end_tvidx = *next_window_end_tvidx_ptr;
  if (next_window_end_tvidx == subcontig_end_tvidx) {
    // just completed last window in subcontig
    *cur_window_size_ptr = 0;
    *window_start_tvidx_ptr = subcontig_end_tvidx;
    fill_ulong_zero(window_maxl, occupied_window_slots);
    return;
  }
  uint32_t next_window_start_tvidx = *window_start_tvidx_ptr;
  if (variant_bps) {
    // this is guaranteed to be nonnegative
    uint32_t variant_uidx_winstart = *variant_uidx_winstart_ptr;
    uint32_t variant_uidx_winend = *variant_uidx_winend_ptr;
    const uint32_t window_start_min_bp = variant_bps[variant_uidx_winend] - prune_window_size;
    uint32_t window_start_bp;
    do {
      // advance window start by as much as necessary to make end advance by at
      // least 1
      ++next_window_start_tvidx;
      ++variant_uidx_winstart;
      next_set_unsafe_ck(variant_include, &variant_uidx_winstart);
      window_start_bp = variant_bps[variant_uidx_winstart];
    } while (window_start_bp < window_start_min_bp);
    // now advance window end as appropriate
    const uint32_t window_end_thresh = window_start_bp + prune_window_size;
    do {
      if (++next_window_end_tvidx == subcontig_end_tvidx) {
	break;
      }
      ++variant_uidx_winend;
      next_set_unsafe_ck(variant_include, &variant_uidx_winend);
    } while (variant_bps[variant_uidx_winend] <= window_end_thresh);
    *variant_uidx_winstart_ptr = variant_uidx_winstart;
    *variant_uidx_winend_ptr = variant_uidx_winend;
  } else {
    next_window_start_tvidx += window_incr;
    next_window_end_tvidx = MINV(next_window_start_tvidx + prune_window_size, subcontig_end_tvidx);
  }
  const uint32_t cur_window_size = *cur_window_size_ptr;
  uint32_t winpos_write = 0;
  for (uint32_t winpos_read = 0; winpos_read < cur_window_size; ++winpos_read) {
    const uint32_t slot_idx = winpos_to_slot_idx[winpos_read];
    if (IS_SET(cur_window_removed, winpos_read) || (tvidxs[slot_idx] < next_window_start_tvidx)) {
      CLEAR_BIT(slot_idx, occupied_window_slots);
    } else {
      winpos_to_slot_idx[winpos_write++] = slot_idx;
    }
  }
  *cur_window_size_ptr = winpos_write;
  *window_start_tvidx_ptr = next_window_start_tvidx;
  *next_window_end_tvidx_ptr = next_window_end_tvidx;
}

void compute_indep_pairwise_r2_components(const uintptr_t* __restrict first_genobufs, const uintptr_t* __restrict second_genobufs, const int32_t* __restrict second_vstats, uint32_t founder_ct, uint32_t* cur_nm_ct_ptr, int32_t* cur_first_sum_ptr, uint32_t* cur_first_ssq_ptr, int32_t* second_sum_ptr, uint32_t* second_ssq_ptr, int32_t* cur_dotprod_ptr) {
  const uint32_t founder_ctaw = BITCT_TO_ALIGNED_WORDCT(founder_ct);
  const uint32_t founder_ctl = BITCT_TO_WORDCT(founder_ct);
  *cur_dotprod_ptr = dotprod_longs(first_genobufs, &(first_genobufs[founder_ctaw]), second_genobufs, &(second_genobufs[founder_ctaw]), founder_ctl);
  if (*cur_nm_ct_ptr != founder_ct) {
    uint32_t plusone_ct;
    uint32_t minusone_ct;
    popcount_longs_2intersect(&(first_genobufs[2 * founder_ctaw]), second_genobufs, &(second_genobufs[founder_ctaw]), founder_ctl, &plusone_ct, &minusone_ct);
    *second_sum_ptr = ((int32_t)plusone_ct) - ((int32_t)minusone_ct);
    *second_ssq_ptr = plusone_ct + minusone_ct;
  } else {
    *second_sum_ptr = second_vstats[1];
    *second_ssq_ptr = second_vstats[2];
  }
  const uint32_t second_nm_ct = second_vstats[0];
  if (second_nm_ct == founder_ct) {
    // assumed that cur_first_nm initialized to first_vstats[0], cur_first_sum
    // initialized to first_vstats[1], cur_first_ssq to first_vstats[2]
    return;
  }
  uint32_t plusone_ct;
  uint32_t minusone_ct;
  popcount_longs_2intersect(&(second_genobufs[2 * founder_ctaw]), first_genobufs, &(first_genobufs[founder_ctaw]), founder_ctl, &plusone_ct, &minusone_ct);
  *cur_first_sum_ptr = ((int32_t)plusone_ct) - ((int32_t)minusone_ct);
  *cur_first_ssq_ptr = plusone_ct + minusone_ct;
  if (*cur_nm_ct_ptr == founder_ct) {
    *cur_nm_ct_ptr = second_nm_ct;
    return;
  }
  *cur_nm_ct_ptr = popcount_longs_intersect(&(first_genobufs[2 * founder_ctaw]), &(second_genobufs[2 * founder_ctaw]), founder_ctl);
}

// multithread globals
static const uint32_t* g_subcontig_info = nullptr;
static const uint32_t* g_subcontig_thread_assignments = nullptr;
static const uintptr_t* g_variant_include = nullptr;
static const uintptr_t* g_variant_allele_idxs = nullptr;
static const alt_allele_ct_t* g_maj_alleles = nullptr;
static const double* g_all_allele_freqs = nullptr;
static const uint32_t* g_variant_bps = nullptr;
static uint32_t* g_tvidx_end = nullptr;
static uint32_t g_x_start = 0;
static uint32_t g_x_len = 0;
static uint32_t g_y_start = 0;
static uint32_t g_y_len = 0;
static uint32_t g_founder_ct = 0;
static uint32_t g_founder_male_ct = 0;
static uint32_t g_prune_window_size = 0;
static uint32_t g_window_maxl = 0;
static double g_prune_ld_thresh = 0.0;
static uint32_t g_window_incr = 0;
static uint32_t g_cur_batch_size = 0;
static uintptr_t** g_genobufs = nullptr;
static uintptr_t** g_occupied_window_slots = nullptr;
static uintptr_t** g_cur_window_removed = nullptr;
static double** g_cur_maj_freqs = nullptr;
static uintptr_t** g_removed_variants_write = nullptr;
static int32_t** g_vstats = nullptr;
static int32_t** g_nonmale_vstats = nullptr;
static uint32_t** g_winpos_to_slot_idx = nullptr;
static uint32_t** g_tvidxs = nullptr;
static uint32_t** g_first_unchecked_tvidx = nullptr;
static uintptr_t** g_raw_tgenovecs[2] = {nullptr, nullptr};

THREAD_FUNC_DECL indep_pairwise_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t* subcontig_info = g_subcontig_info;
  const uint32_t* subcontig_thread_assignments = g_subcontig_thread_assignments;
  const uintptr_t* variant_include = g_variant_include;
  const uint32_t x_start = g_x_start;
  const uint32_t x_len = g_x_len;
  const uint32_t y_start = g_y_start;
  const uint32_t y_len = g_y_len;
  const uintptr_t* variant_allele_idxs = g_variant_allele_idxs;
  const alt_allele_ct_t* maj_alleles = g_maj_alleles;
  const double* all_allele_freqs = g_all_allele_freqs;
  const uint32_t* variant_bps = g_variant_bps;
  const uint32_t founder_ct = g_founder_ct;
  const uint32_t founder_male_ct = g_founder_male_ct;
  const uint32_t founder_male_ctl2 = QUATERCT_TO_WORDCT(founder_male_ct);
  const uint32_t nonmale_ct = founder_ct - founder_male_ct;
  const uint32_t nonmale_ctaw = BITCT_TO_ALIGNED_WORDCT(nonmale_ct);
  const uint32_t nonmale_ctl = BITCT_TO_WORDCT(nonmale_ct);
  const uintptr_t raw_tgenovec_single_variant_word_ct = round_up_pow2(QUATERCT_TO_WORDCT(nonmale_ct) + founder_male_ctl2, kWordsPerVec);
  const uint32_t prune_window_size = g_prune_window_size;
  const uint32_t window_maxl = g_window_maxl;
  const double prune_ld_thresh = g_prune_ld_thresh;
  const uint32_t window_incr = g_window_incr;
  const uint32_t tvidx_end = g_tvidx_end[tidx];
  uintptr_t* genobufs = g_genobufs[tidx];
  uintptr_t* occupied_window_slots = g_occupied_window_slots[tidx];
  uintptr_t* cur_window_removed = g_cur_window_removed[tidx];
  uintptr_t* removed_variants_write = g_removed_variants_write[tidx];
  double* cur_maj_freqs = g_cur_maj_freqs[tidx];
  int32_t* vstats = g_vstats[tidx];
  int32_t* nonmale_vstats = g_nonmale_vstats[tidx];
  uint32_t* winpos_to_slot_idx = g_winpos_to_slot_idx[tidx];
  uint32_t* tvidxs = g_tvidxs[tidx];
  uint32_t* first_unchecked_tvidx = g_first_unchecked_tvidx[tidx];
  
  uint32_t subcontig_end_tvidx = 0;
  uint32_t subcontig_idx = 0xffffffffU; // deliberate overflow
  uint32_t window_start_tvidx = 0;
  uint32_t next_window_end_tvidx = 0;
  uint32_t write_slot_idx = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t cur_window_size = 0;
  uint32_t tvidx_start = 0;
  uint32_t cur_founder_ct = founder_ct;
  uint32_t cur_founder_ctaw = BITCT_TO_ALIGNED_WORDCT(founder_ct);
  uint32_t cur_founder_ctl = BITCT_TO_WORDCT(founder_ct);
  uint32_t variant_uidx = 0;
  uint32_t variant_uidx_winstart = 0;
  uint32_t variant_uidx_winend = 0;
  uintptr_t entire_variant_buf_word_ct = 3 * cur_founder_ctaw;
  uint32_t cur_allele_ct = 2;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uint32_t cur_batch_size = g_cur_batch_size;
    const uint32_t tvidx_stop = MINV(tvidx_start + cur_batch_size, tvidx_end);
    // main loop has to be variant-, not window-, based due to how datasets too
    // large to fit in memory are handled: we may have to halt in the middle of
    // unpacking data for a window, waiting until the current I/O pass is
    // complete before proceeding
    const uintptr_t* raw_tgenovecs = g_raw_tgenovecs[parity][tidx];
    for (uint32_t cur_tvidx = tvidx_start; cur_tvidx < tvidx_stop; ++variant_uidx) {
      if (cur_tvidx == subcontig_end_tvidx) {
	ldprune_next_subcontig(variant_include, variant_bps, subcontig_info, subcontig_thread_assignments, x_start, x_len, y_start, y_len, founder_ct, founder_male_ct, prune_window_size, tidx, &subcontig_idx, &subcontig_end_tvidx, &next_window_end_tvidx, &is_x, &is_y, &cur_founder_ct, &cur_founder_ctaw, &cur_founder_ctl, &entire_variant_buf_word_ct, &variant_uidx_winstart, &variant_uidx_winend);
	variant_uidx = variant_uidx_winstart;
      }
      next_set_unsafe_ck(variant_include, &variant_uidx);
      write_slot_idx = next_unset_unsafe(occupied_window_slots, write_slot_idx);
      uintptr_t tvidx_offset = cur_tvidx - tvidx_start;
      const uintptr_t* cur_raw_tgenovecs = &(raw_tgenovecs[tvidx_offset * raw_tgenovec_single_variant_word_ct]);
      uintptr_t* cur_genobuf = &(genobufs[write_slot_idx * entire_variant_buf_word_ct]);
      uintptr_t* cur_genobuf_minus = &(cur_genobuf[cur_founder_ctaw]);
      uintptr_t* cur_genobuf_nm = &(cur_genobuf_minus[cur_founder_ctaw]);
      genoarr_split_02nm(cur_raw_tgenovecs, cur_founder_ct, cur_genobuf, cur_genobuf_minus, cur_genobuf_nm);
      uint32_t nm_ct = popcount_longs(cur_genobuf_nm, cur_founder_ctl);
      uint32_t plusone_ct = popcount_longs(cur_genobuf, cur_founder_ctl);
      uint32_t minusone_ct = popcount_longs(cur_genobuf_minus, cur_founder_ctl);
      vstats[3 * write_slot_idx] = nm_ct;
      vstats[3 * write_slot_idx + 1] = ((int32_t)plusone_ct) - ((int32_t)minusone_ct);
      vstats[3 * write_slot_idx + 2] = plusone_ct + minusone_ct;
      if (is_x) {
	cur_genobuf = &(cur_genobuf[3 * cur_founder_ctaw]);
	cur_genobuf_minus = &(cur_genobuf[nonmale_ctaw]);
	cur_genobuf_nm = &(cur_genobuf_minus[nonmale_ctaw]);
	genoarr_split_02nm(&(cur_raw_tgenovecs[founder_male_ctl2]), nonmale_ct, cur_genobuf, cur_genobuf_minus, cur_genobuf_nm);
	const uint32_t x_nonmale_nm_ct = popcount_longs(cur_genobuf_nm, nonmale_ctl);
	const uint32_t x_nonmale_plusone_ct = popcount_longs(cur_genobuf, nonmale_ctl);
	const uint32_t x_nonmale_minusone_ct = popcount_longs(cur_genobuf_minus, nonmale_ctl);
	nonmale_vstats[3 * write_slot_idx] = x_nonmale_nm_ct;
	nonmale_vstats[3 * write_slot_idx + 1] = ((int32_t)x_nonmale_plusone_ct) - ((int32_t)x_nonmale_minusone_ct);
	nonmale_vstats[3 * write_slot_idx + 2] = x_nonmale_plusone_ct + x_nonmale_minusone_ct;
	nm_ct += 2 * x_nonmale_nm_ct;
	plusone_ct += 2 * x_nonmale_plusone_ct;
	minusone_ct += 2 * x_nonmale_minusone_ct;
      }
      if (((!plusone_ct) && (!minusone_ct)) || (plusone_ct == nm_ct) || (minusone_ct == nm_ct)) {
	SET_BIT(cur_window_size, cur_window_removed);
	SET_BIT(cur_tvidx, removed_variants_write);
      } else {
	tvidxs[write_slot_idx] = cur_tvidx;
	uintptr_t allele_idx_base;
	if (!variant_allele_idxs) {
	  allele_idx_base = variant_uidx;
	} else {
	  allele_idx_base = variant_allele_idxs[variant_uidx];
	  cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - allele_idx_base;
	  allele_idx_base -= variant_uidx;
	}
	cur_maj_freqs[write_slot_idx] = get_allele_freq(&(all_allele_freqs[allele_idx_base]), maj_alleles[variant_uidx], cur_allele_ct);
	first_unchecked_tvidx[write_slot_idx] = cur_tvidx + 1;
      }
      SET_BIT(write_slot_idx, occupied_window_slots);
      winpos_to_slot_idx[cur_window_size++] = write_slot_idx;
      // are we at the end of a window?
      if (++cur_tvidx == next_window_end_tvidx) {
	// possible for cur_window_size == 1, if all variants at the end of the
	// previous window were pruned
	uint32_t cur_removed_ct = popcount_longs(cur_window_removed, BITCT_TO_WORDCT(cur_window_size));
	uint32_t prev_removed_ct;
	do {
	  prev_removed_ct = cur_removed_ct;
	  uint32_t first_winpos = 0;
	  // const uint32_t debug_print = (!IS_SET(cur_window_removed, 0)) && (tvidxs[winpos_to_slot_idx[0]] == 0);
	  while (1) {
	    next_unset_unsafe_ck(cur_window_removed, &first_winpos);
	    // can assume empty trailing bit for cur_window_removed
	    if (first_winpos == cur_window_size) {
	      break;
	    }
	    uint32_t first_slot_idx = winpos_to_slot_idx[first_winpos];
	    const uint32_t cur_first_unchecked_tvidx = first_unchecked_tvidx[first_slot_idx];
	    uint32_t second_winpos = first_winpos;
	    while (1) {
	      ++second_winpos;
	      next_unset_unsafe_ck(cur_window_removed, &second_winpos);
	      if (second_winpos == cur_window_size) {
		break;
	      }
	      uint32_t second_slot_idx = winpos_to_slot_idx[second_winpos];
	      if (tvidxs[second_slot_idx] >= cur_first_unchecked_tvidx) {
		uintptr_t* first_genobufs = &(genobufs[first_slot_idx * entire_variant_buf_word_ct]);
		const uint32_t first_nm_ct = vstats[3 * first_slot_idx];
		const int32_t first_sum = vstats[3 * first_slot_idx + 1];
		const uint32_t first_ssq = vstats[3 * first_slot_idx + 2];
		while (1) {
		  uintptr_t* second_genobufs = &(genobufs[second_slot_idx * entire_variant_buf_word_ct]);
		  uint32_t cur_nm_ct = first_nm_ct;
		  int32_t cur_first_sum = first_sum;
		  uint32_t cur_first_ssq = first_ssq;
		  int32_t second_sum;
		  uint32_t second_ssq;
		  int32_t cur_dotprod;
		  compute_indep_pairwise_r2_components(first_genobufs, second_genobufs, &(vstats[3 * second_slot_idx]), cur_founder_ct, &cur_nm_ct, &cur_first_sum, &cur_first_ssq, &second_sum, &second_ssq, &cur_dotprod);
		  if (is_x) {
		    uint32_t nonmale_nm_ct = nonmale_vstats[3 * first_slot_idx];
		    int32_t nonmale_first_sum = nonmale_vstats[3 * first_slot_idx + 1];
		    uint32_t nonmale_first_ssq = nonmale_vstats[3 * first_slot_idx + 2];
		    int32_t nonmale_dotprod;
		    int32_t nonmale_second_sum;
		    uint32_t nonmale_second_ssq;
		    compute_indep_pairwise_r2_components(&(first_genobufs[3 * cur_founder_ctaw]), &(second_genobufs[3 * cur_founder_ctaw]), &(nonmale_vstats[3 * second_slot_idx]), nonmale_ct, &nonmale_nm_ct, &nonmale_first_sum, &nonmale_first_ssq, &nonmale_second_sum, &nonmale_second_ssq, &nonmale_dotprod);
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
		  const double cov12 = (double)(cur_dotprod * ((int64_t)cur_nm_ct) - ((int64_t)cur_first_sum) * second_sum);
		  const double variance1 = (double)(cur_first_ssq * ((int64_t)cur_nm_ct) - ((int64_t)cur_first_sum) * cur_first_sum);
		  const double variance2 = (double)(second_ssq * ((int64_t)cur_nm_ct) - ((int64_t)second_sum) * second_sum);
		  // > instead of >=, so we don't prune from a pair of
		  // variants with zero common observations
		  if (cov12 * cov12 > prune_ld_thresh * variance1 * variance2) {
		    // strictly speaking, the (1 + kSmallEpsilon) tolerance
		    // does not appear to be needed yet, but it will be once
		    // --read-freq is implemented.
		    // this has a surprisingly large ~3% speed penalty on my
		    // main test scenario, but that's an acceptable price to
		    // pay for reproducibility.
		    if (cur_maj_freqs[first_slot_idx] > cur_maj_freqs[second_slot_idx] * (1 + kSmallEpsilon)) {
		      /*
		      if (debug_print) {
			printf("removing %u, keeping %u, freqs %g/%g, r2 = %g\n", tvidxs[first_slot_idx], tvidxs[second_slot_idx], cur_maj_freqs[first_slot_idx], cur_maj_freqs[second_slot_idx], cov12 * cov12 / (variance1 * variance2));
		      }
		      */
		      SET_BIT(first_winpos, cur_window_removed);
		      SET_BIT(tvidxs[first_slot_idx], removed_variants_write);
		    } else {
		      /*
		      if (debug_print) {
		        printf("removing %u (second), keeping %u, freqs %g/%g, r2 = %g\n", tvidxs[second_slot_idx], tvidxs[first_slot_idx], cur_maj_freqs[second_slot_idx], cur_maj_freqs[first_slot_idx], cov12 * cov12 / (variance1 * variance2));
		      }
		      */
		      SET_BIT(second_winpos, cur_window_removed);
		      SET_BIT(tvidxs[second_slot_idx], removed_variants_write);
		      const uint32_t next_start_winpos = next_unset_unsafe(cur_window_removed, second_winpos);
		      if (next_start_winpos < cur_window_size) {
			first_unchecked_tvidx[first_slot_idx] = tvidxs[winpos_to_slot_idx[next_start_winpos]];
		      } else {
			first_unchecked_tvidx[first_slot_idx] = cur_tvidx;
		      }
		    }
		    break;
		  }
		  ++second_winpos;
		  next_unset_unsafe_ck(cur_window_removed, &second_winpos);
		  if (second_winpos == cur_window_size) {
		    first_unchecked_tvidx[first_slot_idx] = cur_tvidx;
		    break;
		  }
		  second_slot_idx = winpos_to_slot_idx[second_winpos];
		} // while (1)
		break;
	      }
	    }
	    ++first_winpos;
	  }
	  cur_removed_ct = popcount_longs(cur_window_removed, BITCT_TO_WORDCT(cur_window_size));
	} while (cur_removed_ct > prev_removed_ct);
	const uint32_t prev_window_size = cur_window_size;
	ldprune_next_window(variant_include, variant_bps, tvidxs, cur_window_removed, prune_window_size, window_incr, window_maxl, subcontig_end_tvidx, &cur_window_size, &window_start_tvidx, &variant_uidx_winstart, &next_window_end_tvidx, &variant_uidx_winend, occupied_window_slots, winpos_to_slot_idx);
	// clear bits here since we set cur_window_removed bits during loading
	// process in monomorphic case
	fill_ulong_zero(BITCT_TO_WORDCT(prev_window_size), cur_window_removed);
	write_slot_idx = 0;
      }
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
    tvidx_start = tvidx_stop;
  }
}

pglerr_t indep_pairwise(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, const uintptr_t* variant_allele_idxs, const alt_allele_ct_t* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uint32_t* founder_info_cumulative_popcounts, const uintptr_t* founder_nonmale, const uintptr_t* founder_male, const ld_info_t* ldip, const uint32_t* subcontig_info, const uint32_t* subcontig_thread_assignments, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t founder_male_ct, uint32_t subcontig_ct, uintptr_t window_max, uint32_t calc_thread_ct, uint32_t max_load, pgen_reader_t* simple_pgrp, uintptr_t* removed_variants_collapsed) {
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t founder_nonmale_ct = founder_ct - founder_male_ct;
    if (founder_nonmale_ct * 2 + founder_male_ct > 0x7fffffffU) {
      // may as well document this
      logerrprint("Error: --indep-pairwise does not support >= 2^30 founders.\n");
      goto indep_pairwise_ret_NOT_YET_SUPPORTED;
    }
    const uint32_t founder_nonmale_ctaw = BITCT_TO_ALIGNED_WORDCT(founder_nonmale_ct);
    const uint32_t founder_male_ctaw = BITCT_TO_ALIGNED_WORDCT(founder_male_ct);
    // Per-thread allocations:
    // - tvidx_batch_size * raw_tgenovec_single_variant_word_ct *
    //     sizeof(intptr_t) for raw genotype data (g_raw_tgenovecs)
    // - tvidx_batch_size * sizeof(double) for g_maj_freqs
    // - if pos-based window, tvidx_batch_size * sizeof(int32_t)
    // - All of the above again, to allow loader thread to operate
    //     independently
    // - window_max * 3 * (founder_nonmale_ctaw + founder_male_ctaw) *
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
    uintptr_t* tmp_genovec;
    uint32_t* thread_last_subcontig;
    uint32_t* thread_subcontig_start_tvidx;
    uint32_t* thread_last_tvidx;
    uint32_t* thread_last_uidx;
    pthread_t* threads = nullptr;
    if (bigstack_alloc_ul(QUATERCT_TO_WORDCT(raw_sample_ct), &tmp_genovec) ||
	bigstack_calloc_ui(calc_thread_ct, &g_tvidx_end) ||
	bigstack_calloc_ui(calc_thread_ct, &thread_last_subcontig) ||
	bigstack_calloc_ui(calc_thread_ct, &thread_subcontig_start_tvidx) ||
	bigstack_calloc_ui(calc_thread_ct, &thread_last_tvidx) ||
	bigstack_calloc_ui(calc_thread_ct, &thread_last_uidx) ||
	bigstack_alloc_ulp(calc_thread_ct, &g_genobufs) ||
	bigstack_alloc_ulp(calc_thread_ct, &g_occupied_window_slots) ||
        bigstack_alloc_ulp(calc_thread_ct, &g_cur_window_removed) ||
	bigstack_alloc_dp(calc_thread_ct, &g_cur_maj_freqs) ||
	bigstack_alloc_ulp(calc_thread_ct, &g_removed_variants_write) ||
	bigstack_alloc_ip(calc_thread_ct, &g_vstats) ||
	bigstack_alloc_ip(calc_thread_ct, &g_nonmale_vstats) ||
	bigstack_alloc_uip(calc_thread_ct, &g_winpos_to_slot_idx) ||
	bigstack_alloc_uip(calc_thread_ct, &g_tvidxs) ||
	bigstack_alloc_uip(calc_thread_ct, &g_first_unchecked_tvidx) ||
	bigstack_alloc_ulp(calc_thread_ct, &(g_raw_tgenovecs[0])) ||
        bigstack_alloc_ulp(calc_thread_ct, &(g_raw_tgenovecs[1])) ||
	bigstack_alloc_thread(calc_thread_ct, &threads)) {
      goto indep_pairwise_ret_NOMEM;
    }
    for (uint32_t subcontig_idx = 0; subcontig_idx < subcontig_ct; ++subcontig_idx) {
      const uint32_t cur_thread_idx = subcontig_thread_assignments[subcontig_idx];
      g_tvidx_end[cur_thread_idx] += subcontig_info[3 * subcontig_idx];
    }
    const uintptr_t entire_variant_buf_word_ct = 3 * (founder_nonmale_ctaw + founder_male_ctaw);
    const uint32_t window_maxl = BITCT_TO_WORDCT(window_max);
    const uint32_t max_loadl = BITCT_TO_WORDCT(max_load);
    const uintptr_t genobuf_alloc = round_up_pow2(window_max * entire_variant_buf_word_ct * sizeof(intptr_t), kCacheline);
    const uintptr_t occupied_window_slots_alloc = round_up_pow2(window_maxl * sizeof(intptr_t), kCacheline);
    const uintptr_t cur_window_removed_alloc = round_up_pow2((1 + window_max / kBitsPerWord) * sizeof(intptr_t), kCacheline);
    const uintptr_t cur_maj_freqs_alloc = round_up_pow2(window_max * sizeof(double), kCacheline);
    const uintptr_t removed_variants_write_alloc = round_up_pow2(max_loadl * sizeof(intptr_t), kCacheline);
    const uintptr_t vstats_alloc = round_up_pow2(3 * window_max * sizeof(int32_t), kCacheline); // two of these
    const uintptr_t window_int32_alloc = round_up_pow2(window_max * sizeof(int32_t), kCacheline); // three of these
    const uintptr_t thread_alloc_base = genobuf_alloc + occupied_window_slots_alloc + cur_window_removed_alloc + cur_maj_freqs_alloc + removed_variants_write_alloc + 2 * vstats_alloc + 3 * window_int32_alloc;

    const uint32_t founder_ctl2 = QUATERCT_TO_WORDCT(founder_ct);
    const uint32_t founder_male_ctl2 = QUATERCT_TO_WORDCT(founder_male_ct);
    const uint32_t founder_nonmale_ctl2 = QUATERCT_TO_WORDCT(founder_nonmale_ct);
    const uintptr_t raw_tgenovec_single_variant_word_ct = round_up_pow2(founder_nonmale_ctl2 + founder_male_ctl2, kWordsPerVec);
    // round down
    uintptr_t bigstack_avail_per_thread = round_down_pow2(bigstack_left() / calc_thread_ct, kCacheline);
    // may as well require capacity for >= 256 variants per thread per pass
    if (bigstack_avail_per_thread <= thread_alloc_base + 2 * 256 * raw_tgenovec_single_variant_word_ct * sizeof(intptr_t)) {
      goto indep_pairwise_ret_NOMEM;
    }
    bigstack_avail_per_thread -= thread_alloc_base;
    uint32_t tvidx_batch_size = DIV_UP(max_load, 2);
    // tried a bunch of powers of two, this seems to be a good value
    if (tvidx_batch_size > 65536) {
      tvidx_batch_size = 65536;
    }
    // tvidx_batch_size = max_load; // temporary debugging
    if (2 * tvidx_batch_size * raw_tgenovec_single_variant_word_ct * sizeof(intptr_t) > bigstack_avail_per_thread) {
      tvidx_batch_size = bigstack_avail_per_thread / round_up_pow2(raw_tgenovec_single_variant_word_ct * 2 * sizeof(intptr_t), kCacheline);
    }
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_genobufs[tidx] = (uintptr_t*)bigstack_alloc_raw(genobuf_alloc);
      g_occupied_window_slots[tidx] = (uintptr_t*)bigstack_alloc_raw(occupied_window_slots_alloc);
      fill_ulong_zero(window_maxl, g_occupied_window_slots[tidx]);
      g_cur_window_removed[tidx] = (uintptr_t*)bigstack_alloc_raw(cur_window_removed_alloc);
      fill_ulong_zero(1 + window_max / kBitsPerWord, g_cur_window_removed[tidx]);
      g_cur_maj_freqs[tidx] = (double*)bigstack_alloc_raw(cur_maj_freqs_alloc);
      g_removed_variants_write[tidx] = (uintptr_t*)bigstack_alloc_raw(removed_variants_write_alloc);
      fill_ulong_zero(max_loadl, g_removed_variants_write[tidx]);
      g_vstats[tidx] = (int32_t*)bigstack_alloc_raw(vstats_alloc);
      g_nonmale_vstats[tidx] = (int32_t*)bigstack_alloc_raw(vstats_alloc);
      g_winpos_to_slot_idx[tidx] = (uint32_t*)bigstack_alloc_raw(window_int32_alloc);
      g_tvidxs[tidx] = (uint32_t*)bigstack_alloc_raw(window_int32_alloc);
      g_first_unchecked_tvidx[tidx] = (uint32_t*)bigstack_alloc_raw(window_int32_alloc);
      g_raw_tgenovecs[0][tidx] = (uintptr_t*)bigstack_alloc_raw_rd(tvidx_batch_size * raw_tgenovec_single_variant_word_ct * sizeof(intptr_t));
      g_raw_tgenovecs[1][tidx] = (uintptr_t*)bigstack_alloc_raw_rd(tvidx_batch_size * raw_tgenovec_single_variant_word_ct * sizeof(intptr_t));
    }
    g_subcontig_info = subcontig_info;
    g_subcontig_thread_assignments = subcontig_thread_assignments;
    g_variant_include = variant_include;
    g_variant_allele_idxs = variant_allele_idxs;
    g_maj_alleles = maj_alleles;
    g_all_allele_freqs = allele_freqs;
    g_variant_bps = variant_bps;
    g_founder_ct = founder_ct;
    g_founder_male_ct = founder_male_ct;
    g_prune_ld_thresh = ldip->prune_last_param * (1 + kSmallEpsilon);
    g_prune_window_size = ldip->prune_window_size;
    g_window_maxl = window_maxl;
    g_window_incr = ldip->prune_window_incr;
    g_cur_batch_size = tvidx_batch_size;

    const uint32_t all_haploid = IS_SET(cip->haploid_mask, 0);
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    get_xymt_start_and_end(cip, kChrOffsetX, &x_start, &x_end);
    get_xymt_start_and_end(cip, kChrOffsetY, &y_start, &y_end);
    const uint32_t x_len = x_end - x_start;
    const uint32_t y_len = y_end - y_start;
    g_x_start = x_start;
    g_x_len = x_len;
    g_y_start = y_start;
    g_y_len = y_len;
    // Main workflow:
    // 1. Set n=0, load batch 0
    
    // 2. Spawn threads processing batch n
    // 3. Increment n by 1
    // 4. Load batch n unless eof
    // 5. Join threads
    // 6. Goto step 2 unless eof
    //
    // 7. Assemble final results with copy_bitarr_range()
    uint32_t cur_tvidx_start = 0;
    uint32_t is_last_batch = 0;
    uint32_t parity = 0;
    uint32_t pct = 0;
    uint32_t next_print_tvidx_start = max_load / 100;
    LOGPRINTF("--indep-pairwise (%u compute thread%s): ", calc_thread_ct, (calc_thread_ct == 1)? "" : "s");
    fputs("0%", stdout);
    fflush(stdout);
    while (1) {
      if (!is_last_batch) {
	pgr_clear_ld_cache(simple_pgrp);
	uintptr_t** cur_raw_tgenovecs = g_raw_tgenovecs[parity];
	const uint32_t cur_tvidx_end = cur_tvidx_start + tvidx_batch_size;
	uint32_t is_x_or_y = 0;
	for (uint32_t subcontig_idx = 0; subcontig_idx < subcontig_ct; ++subcontig_idx) {
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
	  const uint32_t is_haploid = IS_SET(cip->haploid_mask, get_variant_chr(cip, variant_uidx));
	  uint32_t is_x = ((variant_uidx - x_start) < x_len);
	  const uint32_t new_is_x_or_y = is_x || ((variant_uidx - y_start) < y_len);

	  // due to nonempty subset requirement (removed?)
	  is_x = is_x && founder_nonmale_ct;
	  if (is_x_or_y != new_is_x_or_y) {
	    is_x_or_y = new_is_x_or_y;
	    pgr_clear_ld_cache(simple_pgrp);
	  }
	  uintptr_t* cur_thread_raw_tgenovec = cur_raw_tgenovecs[cur_thread_idx];
	  for (uintptr_t tvidx_offset = cur_tvidx - cur_tvidx_start; tvidx_offset < tvidx_offset_end; ++tvidx_offset, ++variant_uidx) {
	    next_set_unsafe_ck(variant_include, &variant_uidx);
	    uintptr_t* cur_raw_tgenovec = &(cur_thread_raw_tgenovec[tvidx_offset * raw_tgenovec_single_variant_word_ct]);
	    if (!is_x_or_y) {
	      reterr = pgr_read_allele_countvec_subset_unsafe(founder_info, founder_info_cumulative_popcounts, founder_ct, variant_uidx, maj_alleles[variant_uidx], simple_pgrp, cur_raw_tgenovec);
	      if (is_haploid) {
		set_het_missing(founder_ctl2, cur_raw_tgenovec);
	      }
	    } else {
	      reterr = pgr_read_allele_countvec_subset_unsafe(nullptr, nullptr, raw_sample_ct, variant_uidx, maj_alleles[variant_uidx], simple_pgrp, tmp_genovec);
	      if (founder_male_ct) {
		copy_quaterarr_nonempty_subset(tmp_genovec, founder_male, raw_sample_ct, founder_male_ct, cur_raw_tgenovec);
		set_het_missing(founder_male_ctl2, cur_raw_tgenovec);
	      }
	      if (is_x) {
	        copy_quaterarr_nonempty_subset(tmp_genovec, founder_nonmale, raw_sample_ct, founder_nonmale_ct, &(cur_raw_tgenovec[founder_male_ctl2]));
		if (all_haploid) {
		  // don't just treat chrX identically to autosomes, since for
		  // doubled haploids we still want to give females 2x the
		  // weight of males.  I think.
		  set_het_missing(founder_nonmale_ctl2, &(cur_raw_tgenovec[founder_male_ctl2]));
		}
	      }
	    }
	    if (reterr) {
	      if (cur_tvidx_start) {
		join_threads2z(calc_thread_ct, 0, threads);
		g_cur_batch_size = 0;
		error_cleanup_threads2z(indep_pairwise_thread, calc_thread_ct, threads);
	      }
	      if (reterr != kPglRetReadFail) {
		logprint("\n");
		logerrprint("Error: Malformed .pgen file.\n");
	      }
	      goto indep_pairwise_ret_1;
	    }
	  }
	  thread_last_tvidx[cur_thread_idx] = tvidx_end;
	  thread_last_uidx[cur_thread_idx] = variant_uidx;
	}
      }
      if (cur_tvidx_start) {
	join_threads2z(calc_thread_ct, is_last_batch, threads);
	if (is_last_batch) {
	  break;
	}
	if (cur_tvidx_start >= next_print_tvidx_start) {
	  if (pct > 10) {
	    putc_unlocked('\b', stdout);
	  }
	  pct = (cur_tvidx_start * 100LLU) / max_load;
	  printf("\b\b%u%%", pct++);
	  fflush(stdout);
	  next_print_tvidx_start = (pct * ((uint64_t)max_load)) / 100;
	}
      }
      is_last_batch = (cur_tvidx_start + tvidx_batch_size >= max_load);
      if (spawn_threads2z(indep_pairwise_thread, calc_thread_ct, is_last_batch, threads)) {
	goto indep_pairwise_ret_THREAD_CREATE_FAIL;
      }
      parity = 1 - parity;
      cur_tvidx_start += tvidx_batch_size;
    }
    fill_uint_zero(calc_thread_ct, thread_subcontig_start_tvidx);
    for (uint32_t subcontig_idx = 0; subcontig_idx < subcontig_ct; ++subcontig_idx) {
      const uint32_t cur_thread_idx = subcontig_thread_assignments[subcontig_idx];
      const uintptr_t* cur_removed_variants = g_removed_variants_write[cur_thread_idx];
      const uint32_t subcontig_len = subcontig_info[3 * subcontig_idx];
      const uint32_t subcontig_idx_start = subcontig_info[3 * subcontig_idx + 1];
      copy_bitarr_range(cur_removed_variants, thread_subcontig_start_tvidx[cur_thread_idx], subcontig_idx_start, subcontig_len, removed_variants_collapsed);
      thread_subcontig_start_tvidx[cur_thread_idx] += subcontig_len;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
  }
  while (0) {
  indep_pairwise_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  indep_pairwise_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  indep_pairwise_ret_NOT_YET_SUPPORTED:
    reterr = kPglRetNotYetSupported;
    break;
  }
 indep_pairwise_ret_1:
  // caller will free memory
  return reterr;
}

pglerr_t indep_pairphase() {
  logerrprint("Error: --indep-pairphase is currently under development.\n");
  return kPglRetNotYetSupported;
}

pglerr_t ld_prune_subcontig_split_all(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, uint32_t prune_window_size, uint32_t* window_max_ptr, uint32_t** subcontig_info_ptr, uint32_t* subcontig_ct_ptr) {
  // variant_bps must be nullptr if window size is not bp-based
  // chr0 assumed to already be removed from variant_include.
  // this will skip over chromosomes/contigs with only 1 variant.
  const uint32_t chr_ct = cip->chr_ct;
  uint32_t* subcontig_info = (uint32_t*)g_bigstack_base;
  uint32_t* subcontig_info_iter = subcontig_info;
  uint32_t* subcontig_info_limit = &(((uint32_t*)g_bigstack_end)[-3]);
  uint32_t window_max = 0;
  uint32_t variant_idx = 0;
  if (variant_bps) {
    window_max = 1;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx < chr_ct; ++chr_fo_idx) {
      const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      uint32_t variant_uidx = next_set(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], chr_end);
      const uint32_t chr_variant_ct = popcount_bit_idx(variant_include, variant_uidx, chr_end);
      const uint32_t variant_idx_end = variant_idx + chr_variant_ct;
      if (chr_variant_ct > 1) {
	uint32_t subcontig_uidx_first = variant_uidx;
	uint32_t subcontig_idx_first = variant_idx;
	uint32_t window_idx_first = variant_idx;
	uint32_t window_uidx_first = variant_uidx;
	uint32_t window_pos_first = variant_bps[variant_uidx];
	uint32_t prev_pos = window_pos_first;
	++variant_idx;
	do {
	  ++variant_uidx;
	  next_set_unsafe_ck(variant_include, &variant_uidx);
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
	    do {
	      ++window_uidx_first;
	      next_set_unsafe_ck(variant_include, &window_uidx_first);
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
    for (uint32_t chr_fo_idx = 0; chr_fo_idx < chr_ct; ++chr_fo_idx) {
      const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      const uint32_t first_variant_uidx = next_set(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], chr_end);
      const uint32_t chr_variant_ct = popcount_bit_idx(variant_include, first_variant_uidx, chr_end);
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
  *subcontig_ct_ptr = ((uintptr_t)(subcontig_info_iter - subcontig_info)) / 3;
  *subcontig_info_ptr = subcontig_info;
  bigstack_finalize_ui(subcontig_info, (*subcontig_ct_ptr) * 3);
  *window_max_ptr = window_max;
  return kPglRetSuccess;
}

// next several functions (including load_balance()) will probably move to
// plink2_common
void minheap64_replace_root(uint32_t heap_size, uint64_t new_root, uint64_t* minheap64_preroot) {
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
  minheap64_replace_root(--heap_size, new_root, minheap64_preroot);
  *heap_size_ptr = heap_size;
}
*/

void minheap64_insert(uint64_t new_entry, uint64_t* minheap64_preroot, uint32_t* heap_size_ptr) {
  // assumes minheap64_preroot[0] == 0
  const uint32_t heap_size = 1 + (*heap_size_ptr);
  *heap_size_ptr = heap_size;
  uint32_t cur_pos = heap_size;
  while (1) {
    uint32_t parent_pos = cur_pos / 2;
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
pglerr_t load_balance(const uint32_t* task_weights, uint32_t task_ct, uint32_t* thread_ct_ptr, uint32_t* thread_assignments, uint32_t* max_load_ptr) {
  // max_load assumed to be initialized to zero
  assert(task_ct);
  const uint32_t orig_thread_ct = *thread_ct_ptr;
  if (orig_thread_ct == 1) {
    fill_uint_zero(task_ct, thread_assignments);
    // replace this with an acc_uint32 call?
    uint32_t max_load = task_weights[0];
    for (uint32_t task_idx = 1; task_idx < task_ct; ++task_idx) {
      max_load += task_weights[task_idx];
    }
    *max_load_ptr = max_load;
    return kPglRetSuccess;
  }
  assert(task_ct >= orig_thread_ct);
  uint64_t* sorted_tagged_weights;
  uint64_t* minheap64_preroot;
  if (bigstack_alloc_ull(task_ct, &sorted_tagged_weights) ||
      bigstack_alloc_ull(orig_thread_ct + 2, &minheap64_preroot)) {
    return kPglRetNomem;
  }
  minheap64_preroot[0] = 0;
  uint64_t* minheap64 = &(minheap64_preroot[1]);
  uint32_t total_weight = 0;
  for (uintptr_t task_idx = 0; task_idx < task_ct; ++task_idx) {
    const uintptr_t cur_weight = task_weights[task_idx];
    total_weight += cur_weight;
    sorted_tagged_weights[task_idx] = (((uint64_t)cur_weight) << 32) + (uint64_t)task_idx;
  }
  uint64_t* sorted_tagged_weights_end = &(sorted_tagged_weights[task_ct]);
#ifdef __cplusplus
  // could try std::nth_element if this is ever a bottleneck
  std::sort(sorted_tagged_weights, sorted_tagged_weights_end, std::greater<uint64_t>());
#else
  qsort(sorted_tagged_weights, task_ct, sizeof(int64_t), uint64cmp_decr);
#endif
  const uint64_t largest_tagged_weight = sorted_tagged_weights[0];
  uint32_t initial_max_load = largest_tagged_weight >> 32;
  uint32_t thread_ct = 1 + (total_weight - 1) / initial_max_load;
  if (thread_ct > orig_thread_ct) {
    thread_ct = orig_thread_ct;
    initial_max_load = 1 + (total_weight - 1) / orig_thread_ct;
  }
  
  for (uintptr_t thread_idx = 1; thread_idx < thread_ct; ++thread_idx) {
    minheap64[thread_idx - 1] = thread_ct - thread_idx;
  }
  minheap64[thread_ct - 1] = largest_tagged_weight & 0xffffffff00000000LLU;
  for (uint32_t thread_idx = thread_ct; thread_idx <= orig_thread_ct; ++thread_idx) {
    minheap64[thread_idx] = 0xffffffffffffffffLLU;
  }
  thread_assignments[(uint32_t)largest_tagged_weight] = 0;
  uint64_t max_load_shifted = (((uint64_t)initial_max_load) << 32) | 0xffffffffLLU;
  uint64_t* best_fit_end = sorted_tagged_weights_end;
  if (task_ct > 8 * orig_thread_ct) {
    // stop best-fit here
    best_fit_end = &(sorted_tagged_weights[8 * orig_thread_ct]);
  }
  uint64_t* sorted_tagged_weights_iter = &(sorted_tagged_weights[1]);
  while (sorted_tagged_weights_iter != best_fit_end) {
    // maintain minheap64 as fully sorted list
    uint64_t cur_tagged_weight = *sorted_tagged_weights_iter++;
    const uint32_t task_idx = (uint32_t)cur_tagged_weight;
    cur_tagged_weight &= 0xffffffff00000000LLU;
    const uintptr_t idxp1 = uint64arr_greater_than(minheap64, thread_ct, max_load_shifted - cur_tagged_weight);
    if (idxp1) {
      uintptr_t idx = idxp1 - 1;
      const uint64_t new_entry = minheap64[idx] + cur_tagged_weight;
      while (1) {
	const uint64_t next_entry = minheap64[idx + 1];
	if (new_entry < next_entry) {
	  break;
	}
	minheap64[idx++] = next_entry;
      }
      thread_assignments[task_idx] = (uint32_t)new_entry;
      minheap64[idx] = new_entry;
    } else if (thread_ct < orig_thread_ct) {
      const uint64_t new_entry = cur_tagged_weight + thread_ct;
      const uintptr_t insert_pt = uint64arr_greater_than(minheap64, thread_ct, new_entry);
      for (uintptr_t thread_idx = thread_ct; thread_idx > insert_pt; --thread_idx) {
	minheap64[thread_idx] = minheap64[thread_idx - 1];
      }
      minheap64[insert_pt] = new_entry;
      thread_assignments[task_idx] = thread_ct++;
    } else {
      // move lowest entry to end of list, shift everything else down
      const uint64_t new_entry = minheap64[0] + cur_tagged_weight;
      for (uint32_t thread_idx = 1; thread_idx < thread_ct; ++thread_idx) {
	minheap64[thread_idx - 1] = minheap64[thread_idx];
      }
      minheap64[thread_ct - 1] = new_entry;
      max_load_shifted = new_entry | 0xffffffffLLU;
      thread_assignments[task_idx] = (uint32_t)new_entry;
    }
  }
  if (best_fit_end != sorted_tagged_weights_end) {
    do {
      const uint64_t cur_heaproot = minheap64[0];
      uint64_t cur_tagged_weight = *sorted_tagged_weights_iter++;
      const uint32_t task_idx = (uint32_t)cur_tagged_weight;
      uint32_t cur_thread = (uint32_t)cur_heaproot;
      cur_tagged_weight &= 0xffffffff00000000LLU;
      uint64_t new_entry = cur_heaproot + cur_tagged_weight;
      if (new_entry > max_load_shifted) {
	if (thread_ct < orig_thread_ct) {
	  thread_assignments[task_idx] = thread_ct;
	  minheap64_insert(cur_tagged_weight + thread_ct, minheap64_preroot, &thread_ct);
	  continue;
	} else {
	  max_load_shifted = new_entry | 0xffffffffLLU;
	}
      }
      thread_assignments[task_idx] = cur_thread;
      minheap64_replace_root(thread_ct, new_entry, minheap64_preroot);
    } while (sorted_tagged_weights_iter != sorted_tagged_weights_end);
  }  
  bigstack_reset(sorted_tagged_weights);
  *thread_ct_ptr = thread_ct;
  *max_load_ptr = max_load_shifted >> 32;
  return kPglRetSuccess;
}

pglerr_t ld_prune_write(const uintptr_t* variant_include, const uintptr_t* removed_variants_collapsed, char** variant_ids, uint32_t variant_ct, char* outname, char* outname_end) {
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    fputs("Writing...", stdout);
    fflush(stdout);
    strcpy(outname_end, ".prune.in");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto ld_prune_write_ret_OPEN_FAIL;
    }
    char* textbuf = g_textbuf;
    char* write_iter = textbuf;
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);
    uint32_t variant_uidx = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (is_set(removed_variants_collapsed, variant_idx)) {
	continue;
      }
      write_iter = strcpya(write_iter, variant_ids[variant_uidx]);
      append_binary_eoln(&write_iter);
      if (write_iter >= textbuf_flush) {
        if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	  goto ld_prune_write_ret_WRITE_FAIL;
	}
	write_iter = textbuf;
      }
    }
    if (write_iter > textbuf) {
      if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	goto ld_prune_write_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto ld_prune_write_ret_WRITE_FAIL;
    }

    strcpy(&(outname_end[7]), "out");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto ld_prune_write_ret_OPEN_FAIL;
    }
    write_iter = textbuf;
    variant_uidx = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (!is_set(removed_variants_collapsed, variant_idx)) {
	continue;
      }
      write_iter = strcpya(write_iter, variant_ids[variant_uidx]);
      append_binary_eoln(&write_iter);
      if (write_iter >= textbuf_flush) {
        if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	  goto ld_prune_write_ret_WRITE_FAIL;
	}
	write_iter = textbuf;
      }
    }
    if (write_iter > textbuf) {
      if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	goto ld_prune_write_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto ld_prune_write_ret_WRITE_FAIL;
    }
    *outname_end = '\0';
    putc_unlocked('\r', stdout);
    LOGPRINTFWW("Variant lists written to %s.prune.in and %s.prune.out .\n", outname, outname);
  }
  while (0) {
  ld_prune_write_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ld_prune_write_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  return reterr;
}

pglerr_t ld_prune(const uintptr_t* orig_variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, const alt_allele_ct_t* maj_alleles, const double* allele_freqs, const uintptr_t* founder_info, const uintptr_t* sex_male, const ld_info_t* ldip, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t founder_ct, uint32_t max_thread_ct, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  // common initialization between --indep-pairwise and --indep-pairphase
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t is_pairphase = (ldip->prune_modifier / kfLdPrunePairphase) & 1;
    if (founder_ct < 2) {
      LOGERRPRINTF("Warning: Skipping --indep-pair%s since there are less than two founders.\n(--make-founders may come in handy here.)\n", is_pairphase? "phase" : "wise");
      goto ld_prune_ret_1;
    }
    uint32_t skipped_variant_ct = 0;
    if (is_set(cip->chr_mask, 0)) {
      skipped_variant_ct = count_chr_variants_unsafe(orig_variant_include, cip, 0);
    }
    const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
    if (cip->zero_extra_chrs) {
      for (uint32_t chr_idx = cip->max_code + 1; chr_idx < chr_code_end; ++chr_idx) {
	if (is_set(cip->chr_mask, chr_idx)) {
	  skipped_variant_ct += count_chr_variants_unsafe(orig_variant_include, cip, cip->chr_idx_to_foidx[chr_idx]);
	}
      }
    }
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    const uintptr_t* variant_include;
    if (skipped_variant_ct) {
      uintptr_t* new_variant_include;
      if (bigstack_alloc_ul(raw_variant_ctl, &new_variant_include)) {
	goto ld_prune_ret_NOMEM;
      }
      memcpy(new_variant_include, orig_variant_include, raw_variant_ctl * sizeof(intptr_t));
      if (is_set(cip->chr_mask, 0)) {
	const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[0];
	const uint32_t start_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
	clear_bits_nz(start_uidx, cip->chr_fo_vidx_start[chr_fo_idx + 1], new_variant_include);
      }
      if (cip->zero_extra_chrs) {
        for (uint32_t chr_idx = cip->max_code + 1; chr_idx < chr_code_end; ++chr_idx) {
	  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
	  const uint32_t start_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
	  clear_bits_nz(start_uidx, cip->chr_fo_vidx_start[chr_fo_idx + 1], new_variant_include);
	}
      }
      variant_include = new_variant_include;
      variant_ct -= skipped_variant_ct;
      LOGPRINTF("--indep-pair%s: Ignoring %u chromosome 0 variant%s.\n", is_pairphase? "phase" : "wise", skipped_variant_ct, (skipped_variant_ct == 1)? "" : "s");
    } else {
      variant_include = orig_variant_include;
    }

    if (!(ldip->prune_modifier & kfLdPruneWindowBp)) {
      variant_bps = nullptr;
    }
    const uint32_t prune_window_size = ldip->prune_window_size;
    uint32_t* subcontig_info;
    uint32_t window_max;
    uint32_t subcontig_ct;
    if (ld_prune_subcontig_split_all(variant_include, cip, variant_bps, prune_window_size, &window_max, &subcontig_info, &subcontig_ct)) {
      return kPglRetNomem;
    }
    if (!subcontig_ct) {
      LOGERRPRINTF("Warning: Skipping --indep-pair%s since there are no pairs of variants to\nprocess.\n", is_pairphase? "phase" : "wise");
      goto ld_prune_ret_1;
    }
    if (max_thread_ct > 2) {
      --max_thread_ct;
    }
    if (max_thread_ct > subcontig_ct) {
      max_thread_ct = subcontig_ct;
    }
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    const uint32_t variant_ctl = BITCT_TO_WORDCT(variant_ct);
    const uint32_t founder_male_ct = popcount_longs_intersect(founder_info, sex_male, raw_sample_ctl);
    const uint32_t founder_ctl = BITCT_TO_WORDCT(founder_ct);
    uint32_t* founder_info_cumulative_popcounts;
    uintptr_t* founder_nonmale_collapsed;
    uintptr_t* founder_male_collapsed;
    uintptr_t* removed_variants_collapsed;
    uint32_t* subcontig_thread_assignments;
    if (bigstack_alloc_ui(raw_sample_ctl, &founder_info_cumulative_popcounts) ||
	bigstack_alloc_ul(founder_ctl, &founder_nonmale_collapsed) ||
	bigstack_alloc_ul(founder_ctl, &founder_male_collapsed) ||
	bigstack_calloc_ul(variant_ctl, &removed_variants_collapsed) ||
	bigstack_alloc_ui(subcontig_ct, &subcontig_thread_assignments)) {
      goto ld_prune_ret_NOMEM;
    }
    fill_cumulative_popcounts(founder_info, raw_sample_ctl, founder_info_cumulative_popcounts);
    copy_bitarr_subset(sex_male, founder_info, founder_ct, founder_male_collapsed);
    bitarr_invert_copy(founder_male_collapsed, founder_ct, founder_nonmale_collapsed);
    uint32_t* subcontig_weights;
    if (bigstack_end_alloc_ui(subcontig_ct, &subcontig_weights)) {
      goto ld_prune_ret_NOMEM;
    }

    // initial window_max-based memory requirement estimate
    if (is_pairphase) {
      // todo
    } else {
      const uintptr_t entire_variant_buf_word_ct = 3 * (BITCT_TO_ALIGNED_WORDCT(founder_ct - founder_male_ct) + BITCT_TO_ALIGNED_WORDCT(founder_male_ct));
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
      uintptr_t per_thread_alloc = round_up_pow2(window_max * entire_variant_buf_word_ct * sizeof(intptr_t), kCacheline) + 2 * round_up_pow2((1 + window_max / kBitsPerWord) * sizeof(intptr_t), kCacheline) + round_up_pow2(window_max * sizeof(double), kCacheline) + 2 * round_up_pow2(window_max * (3 * sizeof(int32_t)), kCacheline) + 3 * round_up_pow2(window_max * sizeof(int32_t), kCacheline);
      uintptr_t bigstack_left2 = bigstack_left();
      if (per_thread_alloc * max_thread_ct > bigstack_left2) {
	if (per_thread_alloc > bigstack_left2) {
	  goto ld_prune_ret_NOMEM;
	}
	max_thread_ct = bigstack_left2 / per_thread_alloc;
      }
    }

    
    for (uint32_t subcontig_idx = 0; subcontig_idx < subcontig_ct; ++subcontig_idx) {
      // todo: adjust chrX weights upward, and chrY downward
      subcontig_weights[subcontig_idx] = subcontig_info[3 * subcontig_idx];
      // printf("%u %u %u\n", subcontig_info[3 * subcontig_idx], subcontig_info[3 * subcontig_idx + 1], subcontig_info[3 * subcontig_idx + 2]);
    }
    uint32_t max_load = 0;
    if (load_balance(subcontig_weights, subcontig_ct, &max_thread_ct, subcontig_thread_assignments, &max_load)) {
      goto ld_prune_ret_NOMEM;
    }
    bigstack_end_reset(bigstack_end_mark);
    
    if (is_pairphase) {
      reterr = indep_pairphase();
    } else {
      reterr = indep_pairwise(variant_include, cip, variant_bps, variant_allele_idxs, maj_alleles, allele_freqs, founder_info, founder_info_cumulative_popcounts, founder_nonmale_collapsed, founder_male_collapsed, ldip, subcontig_info, subcontig_thread_assignments, raw_sample_ct, founder_ct, founder_male_ct, subcontig_ct, window_max, max_thread_ct, max_load, simple_pgrp, removed_variants_collapsed);
    }
    if (reterr) {
      goto ld_prune_ret_1;
    }
    const uint32_t removed_ct = popcount_longs(removed_variants_collapsed, variant_ctl);
    LOGPRINTF("%u/%u variants removed.\n", removed_ct, variant_ct);
    reterr = ld_prune_write(variant_include, removed_variants_collapsed, variant_ids, variant_ct, outname, outname_end);
    if (reterr) {
      goto ld_prune_ret_1;
    }
  }
  while (0) {
  ld_prune_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 ld_prune_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
