#include "plink_ld.h"
#include "plink_matrix.h"
#include "pigz.h"

#define MULTIPLEX_LD 1920
#define MULTIPLEX_2LD (MULTIPLEX_LD * 2)

void ld_init(Ld_info* ldip) {
  ldip->modifier = 0;
  ldip->prune_window_size = 0;
  ldip->prune_window_incr = 0;
  ldip->prune_window_kb = 0;
  ldip->prune_last_param = 0.0;
  ldip->window_size = 10;
  ldip->window_bp = 200000;
  ldip->window_r2 = 0.2;
  ldip->snpstr = NULL;
  range_list_init(&(ldip->snps_rl));
}

void ld_cleanup(Ld_info* ldip) {
  free_cond(ldip->snpstr);
  free_range_list(&(ldip->snps_rl));
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
  // \mu_x\sum_{i: nonmissing from x} y_i must be handled in the main loop, and
  // this removes much of the speed advantage.  So the best applications of the
  // underlying ternary dot product algorithm used here lie elsewhere.
  // Nevertheless, it is still faster, so we use it.
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
  //   return_vals[1] += N_y + \sum_{i: nonmissing from y}x_i
  //   return_vals[2] += N_x + \sum_{i: nonmissing from x}y_i
  // where N is the number of individuals processed after applying the
  // missingness masks indicated by the subscripts.  The calculation currently
  // proceeds as follows:
  //
  // 1. N + \sum_i x_i = popcount2(vec1 & mask2)
  // The "2" suffix refers to starting with two-bit integers instead of one-bit
  // integers in our summing process, so we get to skip a few operations.
  // (Once we can assume the presence of hardware popcount, a slightly
  // different implementation may be better.)
  //
  // 2. zcheck := (vec1 | vec2) & 0x5555...
  // Detects whether at least one member of the pair has a 0/missing value.
  //
  // 3. popcount2(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
  // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i dot
  // product.
  //
  // MULTIPLEX_LD sets of values are handled per function call.  If fewer
  // values are present, it is currently safe to zero out the ends of all
  // input vectors.

  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader1;
  __m128i loader2;
  __m128i sum1;
  __m128i sum2;
  __m128i sum12;
  __m128i tmp_sum1;
  __m128i tmp_sum2;
  __m128i tmp_sum12;
  __uni16 acc;
  __uni16 acc1;
  __uni16 acc2;
  acc.vi = _mm_setzero_si128();
  acc1.vi = _mm_setzero_si128();
  acc2.vi = _mm_setzero_si128();
  do {
    loader1 = *vec1++;
    loader2 = *vec2++;
    sum1 = *mask2++;
    sum2 = *mask1++;
    sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    sum1 = _mm_and_si128(sum1, loader1);
    sum2 = _mm_and_si128(sum2, loader2);
    // use andnot to eliminate need for 0xaaaa... to occupy an xmm register
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12), _mm_xor_si128(loader1, loader2));
    sum12 = _mm_or_si128(sum12, loader1);

    // sum1, sum2, and sum12 now store the (biased) two-bit sums of
    // interest
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
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

    sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
    sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    acc1.vi = _mm_add_epi64(acc1.vi, _mm_add_epi64(_mm_and_si128(sum1, m4), _mm_and_si128(_mm_srli_epi64(sum1, 4), m4)));
    acc2.vi = _mm_add_epi64(acc2.vi, _mm_add_epi64(_mm_and_si128(sum2, m4), _mm_and_si128(_mm_srli_epi64(sum2, 4), m4)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sum12, m4), _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
  } while (--iters);
  // moved down since we're out of xmm registers
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
  return_vals[0] -= ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[1] += ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[2] += ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
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
#else
static inline void ld_dot_prod_batch(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1, uintptr_t* mask2, int32_t* return_vals, uint32_t iters) {
  uint32_t final_sum1 = 0;
  uint32_t final_sum2 = 0;
  uint32_t final_sum12 = 0;
  uintptr_t loader1;
  uintptr_t loader2;
  uintptr_t sum1;
  uintptr_t sum2;
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
    //   return_vals[1] += N_y + \sum_i x_i
    //   return_vals[2] += N_x + \sum_i y_i
    // where N is the number of individuals processed after applying the
    // missingness masks indicated by the subscripts.  The calculation
    // currently proceeds as follows:
    //
    // 1. N + \sum_i x_i = popcount_variant(vec1 & mask2)
    // The "variant" suffix refers to starting with two-bit integers instead of
    // one-bit integers in our summing process, so we get to skip a few
    // operations.  (Once all reserachers are using machines with fast hardware
    // popcount, a slightly different implementation may be better.)
    //
    // 2. zcheck := (vec1 | vec2) & 0x5555...
    // Detects whether at least one member of the pair has a 0/missing value.
    //
    // 3. popcount_variant(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
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

    sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
    sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    sum1 = (sum1 & 0x0f0f0f0f) + ((sum1 >> 4) & 0x0f0f0f0f);
    sum2 = (sum2 & 0x0f0f0f0f) + ((sum2 >> 4) & 0x0f0f0f0f);
    sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);

    // technically could do the multiply-and-shift only once every two rounds
    final_sum1 += (sum1 * 0x01010101) >> 24;
    final_sum2 += (sum2 * 0x01010101) >> 24;
    final_sum12 += (sum12 * 0x01010101) >> 24;
  } while (--iters);
  return_vals[0] -= final_sum12;
  return_vals[1] += final_sum1;
  return_vals[2] += final_sum2;
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
#endif // __LP64__

uint32_t ld_process_load(uintptr_t* geno_buf, uintptr_t* mask_buf, uintptr_t* missing_buf, double* marker_stdev_ptr, uint32_t founder_ct, uint32_t is_x, uint32_t weighted_x, uint32_t nonmale_founder_ct, uintptr_t* founder_male_include2, uintptr_t* nonmale_geno, uintptr_t* nonmale_masks, uintptr_t nonmale_offset) {
  uintptr_t* geno_ptr = geno_buf;
  uintptr_t founder_ctl2 = (founder_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t* geno_end = &(geno_buf[founder_ctl2]);
  uintptr_t* mask_buf_ptr = mask_buf;
  uintptr_t* missing_ptr = missing_buf;
  uintptr_t new_missing = 0;
  uint32_t missing_bit_offset = 0;
  uint32_t ssq = 0;
  int32_t sum = -founder_ct;
  uintptr_t* nm_mask_ptr;
  double non_missing_recip;
  uintptr_t cur_geno;
  uintptr_t shifted_masked_geno;
  uintptr_t new_geno;
  uintptr_t new_mask;
  uint32_t missing_ct = 0;
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
    missing_ct += founder_ct - (popcount_longs(nonmale_masks, 0, founder_ctl2) / 2);
    founder_ct *= 2;
  }
  non_missing_recip = 1.0 / (founder_ct - missing_ct);
  *marker_stdev_ptr = non_missing_recip * sqrt(((int64_t)((uint64_t)ssq)) * (founder_ct - missing_ct) - ((int64_t)sum) * sum);
  return missing_ct;
}

uint32_t ld_prune_next_valid_chrom_start(uintptr_t* marker_exclude, uint32_t cur_uidx, Chrom_info* chrom_info_ptr, uint32_t unfiltered_marker_ct) {
  uint32_t max_code = chrom_info_ptr->max_code;
  uint32_t chrom_idx;
  cur_uidx = next_unset(marker_exclude, cur_uidx, unfiltered_marker_ct);
  while (cur_uidx < unfiltered_marker_ct) {
    chrom_idx = get_marker_chrom(chrom_info_ptr, cur_uidx);
    if (chrom_idx && (chrom_idx <= max_code)) {
      return cur_uidx;
    }
    cur_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_end[chrom_idx], unfiltered_marker_ct);
  }
  return cur_uidx;
}

void ld_prune_start_chrom(uint32_t ld_window_kb, uint32_t* cur_chrom_ptr, uint32_t* chrom_end_ptr, uint32_t window_unfiltered_start, uint32_t* live_indices, uint32_t* start_arr, uint32_t* window_unfiltered_end_ptr, uint32_t ld_window_size, uint32_t* cur_window_size_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uint32_t* is_haploid_ptr, uint32_t* is_x_ptr, uint32_t* is_y_ptr) {
  uint32_t cur_chrom = get_marker_chrom(chrom_info_ptr, window_unfiltered_start);
  uint32_t window_unfiltered_end = window_unfiltered_start + 1;
  uint32_t chrom_end = chrom_info_ptr->chrom_end[cur_chrom];
  uint32_t uii = 0;
  uint32_t window_size;
  live_indices[0] = window_unfiltered_start;
  if (ld_window_kb) {
    window_size = 0;
    while ((window_unfiltered_start + window_size < chrom_end) && (marker_pos[window_unfiltered_start + window_size] <= marker_pos[window_unfiltered_start] + (1000 * ld_window_size))) {
      window_size++;
    }
  } else {
    window_size = ld_window_size;
  }
  for (uii = 1; uii < window_size; window_unfiltered_end++, uii++) {
    next_unset_ck(marker_exclude, &window_unfiltered_end, chrom_end);
    if (window_unfiltered_end == chrom_end) {
      break;
    }
    start_arr[uii - 1] = window_unfiltered_end;
    live_indices[uii] = window_unfiltered_end;
  }
  *cur_window_size_ptr = uii;
  start_arr[uii - 1] = window_unfiltered_end;
  *cur_chrom_ptr = cur_chrom;
  *chrom_end_ptr = chrom_end;
  *window_unfiltered_end_ptr = window_unfiltered_end;
  *is_haploid_ptr = IS_SET(chrom_info_ptr->haploid_mask, cur_chrom);
  *is_x_ptr = (((int32_t)cur_chrom) == chrom_info_ptr->x_code)? 1 : 0;
  *is_y_ptr = (((int32_t)cur_chrom) == chrom_info_ptr->y_code)? 1 : 0;
}

int32_t ld_prune(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_indiv_ct, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  // for future consideration: chromosome-based multithread/parallel?
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile_in = NULL;
  FILE* outfile_out = NULL;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = 2 * ((unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  uintptr_t founder_ct = popcount_longs(founder_info, 0, unfiltered_indiv_ctl2 / 2);
  uint32_t weighted_founder_ct = founder_ct;
  uintptr_t founder_ctl = (founder_ct + BITCT - 1) / BITCT;
#ifdef __LP64__
  uintptr_t founder_ctv = 2 * ((founder_ct + 127) / 128);
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
  uintptr_t founder_ct_mld_long = founder_ct_mld * (MULTIPLEX_LD / BITCT2);
  uint32_t founder_trail_ct = founder_ct_mld_long - founder_ctl * 2;
  uint32_t pairwise = (ldip->modifier / LD_PRUNE_PAIRWISE) & 1;
  uint32_t ignore_x = (ldip->modifier / LD_IGNORE_X) & 1;
  uint32_t weighted_x = (ldip->modifier / LD_WEIGHTED_X) & 1;
  uint32_t ld_window_size = ldip->prune_window_size;
  uint32_t ld_window_incr = ldip->prune_window_incr;
  double ld_last_param = ldip->prune_last_param;
  uint32_t nonmale_founder_ct = 0;
  uintptr_t window_max = 0;
  uintptr_t* geno = NULL;
  uintptr_t* founder_include2 = NULL;
  uintptr_t* founder_male_include2 = NULL;
  uintptr_t* nonmale_geno = NULL;
  uintptr_t* nonmale_masks = NULL;
  double* cov_matrix = NULL;
  double* new_cov_matrix = NULL;
  MATRIX_INVERT_BUF1_TYPE* irow = NULL;
  double* work = NULL;
  uint32_t* idx_remap = NULL;
  uint32_t tot_exclude_ct = 0;
  uint32_t at_least_one_prune = 0;
  uint32_t max_code = chrom_info_ptr->max_code;
  int32_t retval = 0;
  uintptr_t* geno_masks;
  uintptr_t* geno_mmasks;
  uintptr_t* pruned_arr;
  uint32_t* live_indices;
  uint32_t* start_arr;
  uint32_t marker_unfiltered_idx;
  uintptr_t marker_idx;
  int32_t pct;
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
  uint32_t chrom_end;
  uint32_t is_haploid;
  uint32_t is_x;
  uint32_t is_y;
  double* marker_stdevs;
  uintptr_t* loadbuf;
  uint32_t* missing_cts;
  uint32_t fixed_missing_ct;
  uintptr_t ulii;
  double dxx;
  double cov12;
  uint32_t fixed_non_missing_ct;
  uint32_t non_missing_ct;
  int32_t dp_result[3];
  double non_missing_recip;
  uintptr_t* geno_fixed_vec_ptr;
  uintptr_t* geno_var_vec_ptr;
  uintptr_t* mask_fixed_vec_ptr;
  uintptr_t* mask_var_vec_ptr;
  uintptr_t cur_exclude_ct;
  uint32_t prev_end;
  char* sptr;
  FILE* fptr;
  __CLPK_integer window_rem_li;
  __CLPK_integer old_window_rem_li;
  uint32_t window_rem;
  double prune_ld_r1;
  if (!founder_ct) {
    sprintf(logbuf, "Warning: Skipping --indep%s since there are no founders.\n", pairwise? "-pairwise" : "");
    logprintb();
    goto ld_prune_ret_1;
  }

  // force founder_male_include2 allocation
  if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, founder_ct, XMHH_EXISTS | hh_exists, 1, founder_info, sex_male, &founder_include2, &founder_male_include2)) {
    goto ld_prune_ret_NOMEM;
  }
  if (weighted_x) {
    nonmale_founder_ct = founder_ct - popcount_longs(founder_male_include2, 0, founder_ctl);
    if (founder_ct + nonmale_founder_ct > 0x7fffffff) {
      // no, this shouldn't ever happen, but may as well document that there
      // theoretically is a 32-bit integer range issue here
      logprint("Error: Too many founders for --indep[-pairwise] + --ld-xchr 3.\n");
      goto ld_prune_ret_1;
    }
  }

  if (ldip->prune_window_kb) {
    // determine maximum number of markers that may need to be loaded at once
    for (cur_chrom = 0; cur_chrom <= max_code; cur_chrom++) {
      if (chrom_exists(chrom_info_ptr, cur_chrom)) {
        uii = chrom_info_ptr->chrom_start[cur_chrom];
	chrom_end = chrom_info_ptr->chrom_end[cur_chrom];
        do {
	  ujj = uii + 1;
	  while ((ujj < chrom_end) && (marker_pos[ujj] <= marker_pos[uii] + (1000 * ld_window_size))) {
	    ujj++;
	  }
          if (ujj - uii > window_max) {
	    window_max = ujj - uii;
	  }
	  uii++;
	} while (ujj < chrom_end);
      }
    }
  }
  if (pairwise) {
    prune_ld_r1 = sqrt(ld_last_param);
  } else {
    prune_ld_r1 = 0.999999;
  }

  window_unfiltered_start = ld_prune_next_valid_chrom_start(marker_exclude, 0, chrom_info_ptr, unfiltered_marker_ct);
  if (window_unfiltered_start == unfiltered_marker_ct) {
    sprintf(logbuf, "Error: No valid variants for --indep%s.\n", pairwise? "-pairwise" : "");
    logprintb();
    goto ld_prune_ret_INVALID_FORMAT;
  }

  if (wkspace_alloc_ul_checked(&pruned_arr, unfiltered_marker_ctl * sizeof(intptr_t))) {
    goto ld_prune_ret_NOMEM;
  }

  memcpy(pruned_arr, marker_exclude, unfiltered_marker_ctl * sizeof(intptr_t));

  if (!ldip->prune_window_kb) {
    window_max = ld_window_size;
  }
  ulii = window_max;
  if (wkspace_alloc_ui_checked(&live_indices, ulii * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&start_arr, ulii * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&marker_stdevs, ulii * sizeof(double)) ||
      wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&geno, ulii * founder_ct_mld_long * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&geno_masks, ulii * founder_ct_mld_long * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&geno_mmasks, ulii * founder_ctv * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&missing_cts, ulii * sizeof(int32_t))) {
    goto ld_prune_ret_NOMEM;
  }
  if (weighted_x) {
    if (wkspace_alloc_ul_checked(&nonmale_geno, ulii * founder_ct_mld_long * sizeof(intptr_t)) ||
        wkspace_alloc_ul_checked(&nonmale_masks, ulii * founder_ct_mld_long * sizeof(intptr_t))) {
      goto ld_prune_ret_NOMEM;
    }
  }
  if (founder_trail_ct) {
    for (ulii = 1; ulii <= window_max; ulii++) {
      fill_ulong_zero(&(geno[ulii * founder_ct_mld_long - founder_trail_ct - 2]), founder_trail_ct + 2);
      fill_ulong_zero(&(geno_masks[ulii * founder_ct_mld_long - founder_trail_ct - 2]), founder_trail_ct + 2);
      if (weighted_x) {
	fill_ulong_zero(&(nonmale_geno[ulii * founder_ct_mld_long - founder_trail_ct - 2]), founder_trail_ct + 2);
	fill_ulong_zero(&(nonmale_masks[ulii * founder_ct_mld_long - founder_trail_ct - 2]), founder_trail_ct + 2);
      }
    }
  }
  if (!pairwise) {
    if (wkspace_alloc_d_checked(&cov_matrix, window_max * window_max * sizeof(double)) ||
        wkspace_alloc_d_checked(&new_cov_matrix, window_max * window_max * sizeof(double)) ||
        wkspace_alloc_ui_checked(&idx_remap, window_max * sizeof(int32_t))) {
      goto ld_prune_ret_NOMEM;
    }

    irow = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(window_max * 2 * sizeof(MATRIX_INVERT_BUF1_TYPE));
    if (!irow) {
      goto ld_prune_ret_NOMEM;
    }

    if (window_max < 4) {
      ulii = 4;
    } else {
      ulii = window_max;
    }
    if (wkspace_alloc_d_checked(&work, ulii * window_max * sizeof(double))) {
      goto ld_prune_ret_NOMEM;
    }
  }
  do {
    prev_end = 0;
    ld_prune_start_chrom(ldip->prune_window_kb, &cur_chrom, &chrom_end, window_unfiltered_start, live_indices, start_arr, &window_unfiltered_end, ld_window_size, &cur_window_size, unfiltered_marker_ct, pruned_arr, chrom_info_ptr, marker_pos, &is_haploid, &is_x, &is_y);
    if (weighted_x) {
      if (is_x) {
	weighted_founder_ct = 2 * founder_ct;
      } else {
	weighted_founder_ct = founder_ct;
      }
    }
    old_window_size = 1;
    if (cur_window_size > 1) {
      for (ulii = 0; ulii < (uintptr_t)cur_window_size; ulii++) {
	uii = live_indices[ulii];
	if (fseeko(bedfile, bed_offset + (uii * unfiltered_indiv_ct4), SEEK_SET)) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (load_and_collapse_incl(bedfile, loadbuf, unfiltered_indiv_ct, &(geno[ulii * founder_ct_mld_long]), founder_ct, founder_info, IS_SET(marker_reverse, uii))) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(geno[ulii * founder_ct_mld_long])));
	}
        missing_cts[ulii] = ld_process_load(&(geno[ulii * founder_ct_mld_long]), &(geno_masks[ulii * founder_ct_mld_long]), &(geno_mmasks[ulii * founder_ctv]), &(marker_stdevs[ulii]), founder_ct, is_x && (!ignore_x), weighted_x, nonmale_founder_ct, founder_male_include2, nonmale_geno, nonmale_masks, ulii * founder_ct_mld_long);
      }
    }
    pct = 1;
    pct_thresh = window_unfiltered_start + ((int64_t)pct * (chrom_end - chrom_info_ptr->chrom_start[cur_chrom])) / 100;
    cur_exclude_ct = 0;
    while ((window_unfiltered_start < chrom_end) || (cur_window_size > 1)) {
      if (cur_window_size > 1) {
	for (uii = 0; uii < cur_window_size; uii++) {
	  if (marker_stdevs[uii] == 0.0) {
	    SET_BIT(pruned_arr, live_indices[uii]);
	    cur_exclude_ct++;
	  }
	}
	do {
	  at_least_one_prune = 0;
	  for (uii = 0; uii < cur_window_size - 1; uii++) {
	    if (IS_SET(pruned_arr, live_indices[uii])) {
	      continue;
	    }
            fixed_missing_ct = missing_cts[uii];
	    fixed_non_missing_ct = weighted_founder_ct - fixed_missing_ct;
	    geno_fixed_vec_ptr = &(geno[uii * founder_ct_mld_long]);
	    mask_fixed_vec_ptr = &(geno_masks[uii * founder_ct_mld_long]);
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
	      geno_var_vec_ptr = &(geno[ujj * founder_ct_mld_long]);
	      mask_var_vec_ptr = &(geno_masks[ujj * founder_ct_mld_long]);

	      dp_result[0] = weighted_founder_ct;
	      // reversed from what I initially thought because I'm passing the
	      // ujj-associated buffers before the uii-associated ones.
	      dp_result[1] = -fixed_non_missing_ct;
	      dp_result[2] = missing_cts[ujj] - weighted_founder_ct;
	      ld_dot_prod(geno_var_vec_ptr, geno_fixed_vec_ptr, mask_var_vec_ptr, mask_fixed_vec_ptr, dp_result, founder_ct_mld_m1, founder_ct_mld_rem);
	      if (is_x && weighted_x) {
		non_missing_ct = (popcount_longs_intersect(&(nonmale_masks[uii * founder_ct_mld_long]), &(nonmale_masks[ujj * founder_ct_mld_long]), 2 * founder_ctl) + popcount_longs_intersect(mask_fixed_vec_ptr, mask_var_vec_ptr, 2 * founder_ctl)) / 2;
		ld_dot_prod(&(nonmale_geno[ujj * founder_ct_mld_long]), &(nonmale_geno[uii * founder_ct_mld_long]), &(nonmale_masks[ujj * founder_ct_mld_long]), &(nonmale_masks[uii * founder_ct_mld_long]), dp_result, founder_ct_mld_m1, founder_ct_mld_rem);
	      } else {
	        non_missing_ct = fixed_non_missing_ct - missing_cts[ujj];
		if (fixed_missing_ct && missing_cts[ujj]) {
		  non_missing_ct += popcount_longs_intersect(&(geno_mmasks[uii * founder_ctv]), &(geno_mmasks[ujj * founder_ctv]), founder_ctl);
		}
	      }
	      non_missing_recip = 1.0 / ((double)((int32_t)non_missing_ct));
	      cov12 = non_missing_recip * (dp_result[0] - (non_missing_recip * dp_result[1]) * dp_result[2]);
	      // r, not squared
	      dxx = cov12 / (marker_stdevs[uii] * marker_stdevs[ujj]);
	      if (!pairwise) {
		cov_matrix[uii * window_max + ujj] = dxx;
	      }
	      if (fabs(dxx) > prune_ld_r1) {
		at_least_one_prune = 1;
		cur_exclude_ct++;
		// remove marker with lower MAF
		if (get_maf(set_allele_freqs[live_indices[uii]]) < get_maf(set_allele_freqs[live_indices[ujj]])) {
		  SET_BIT(pruned_arr, live_indices[uii]);
		} else {
		  SET_BIT(pruned_arr, live_indices[ujj]);
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
	  old_window_rem_li = 0;
	  for (uii = 0; uii < old_window_size; uii++) {
	    if (IS_SET(pruned_arr, live_indices[uii])) {
	      continue;
	    }
            idx_remap[window_rem++] = uii;
	  }
	  old_window_rem_li = window_rem;
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
	    ii = invert_matrix_trunc_singular(window_rem_li, new_cov_matrix, irow, work, old_window_rem_li);
	    while (ii) {
	      if (ii == -1) {
		goto ld_prune_ret_NOMEM;
	      }
	      ujj = ii;
              SET_BIT(pruned_arr, live_indices[idx_remap[ujj]]);
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
	      ii = invert_matrix_trunc_singular(window_rem_li, new_cov_matrix, irow, work, old_window_rem_li);
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
	      SET_BIT(pruned_arr, live_indices[idx_remap[ujj]]);
	      cur_exclude_ct++;
	      window_rem--;
	      if (idx_remap[ujj] < (uint32_t)old_window_size) {
		old_window_rem_li--;
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
	while (IS_SET(marker_exclude, window_unfiltered_start)) {
	  if (window_unfiltered_start == chrom_end) {
	    break;
	  }
	  window_unfiltered_start++;
	}
	if (window_unfiltered_start == chrom_end) {
	  break;
	}
	window_unfiltered_start++;
      }
      if (window_unfiltered_start == chrom_end) {
	break;
      }
      if (window_unfiltered_start >= pct_thresh) {
	pct = (((int64_t)(window_unfiltered_start - chrom_info_ptr->chrom_start[cur_chrom])) * 100) / (chrom_end - chrom_info_ptr->chrom_start[cur_chrom]);
	printf("\r%d%%", pct++);
	fflush(stdout);
	pct_thresh = chrom_info_ptr->chrom_start[cur_chrom] + (((int64_t)pct * (chrom_end - chrom_info_ptr->chrom_start[cur_chrom])) / 100);
      }
      ujj = 0;
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
	memcpy(&(geno[uii * founder_ct_mld_long]), &(geno[ujj * founder_ct_mld_long]), founder_ct_mld_long * sizeof(intptr_t));
	memcpy(&(geno_masks[uii * founder_ct_mld_long]), &(geno_masks[ujj * founder_ct_mld_long]), founder_ct_mld_long * sizeof(intptr_t));
	if (is_x && weighted_x) {
	  memcpy(&(nonmale_geno[uii * founder_ct_mld_long]), &(nonmale_geno[ujj * founder_ct_mld_long]), founder_ct_mld_long * sizeof(intptr_t));
	  memcpy(&(nonmale_masks[uii * founder_ct_mld_long]), &(nonmale_masks[ujj * founder_ct_mld_long]), founder_ct_mld_long * sizeof(intptr_t));
	}
	memcpy(&(geno_mmasks[uii * founder_ctv]), &(geno_mmasks[ujj * founder_ctv]), founder_ctl * sizeof(intptr_t));
	marker_stdevs[uii] = marker_stdevs[ujj];
	live_indices[uii] = live_indices[ujj];
	start_arr[uii] = start_arr[ujj];
	missing_cts[uii] = missing_cts[ujj];
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
      if (ldip->prune_window_kb) {
	ujj = 0;
	while ((window_unfiltered_end + ujj < chrom_end) && (marker_pos[window_unfiltered_end + ujj] <= marker_pos[window_unfiltered_start] + (1000 * ld_window_size))) {
	  ujj++;
	}
      } else {
	ujj = ld_window_incr;
      }
      old_window_size = cur_window_size;
      for (uii = 0; uii < ujj; window_unfiltered_end++, uii++) {
	next_unset_ck(marker_exclude, &window_unfiltered_end, chrom_end);
	if (window_unfiltered_end == chrom_end) {
	  break;
	}
	live_indices[cur_window_size] = window_unfiltered_end;
	if (cur_window_size > prev_end) {
	  start_arr[cur_window_size - 1] = window_unfiltered_end;
	}
	if (fseeko(bedfile, bed_offset + (window_unfiltered_end * unfiltered_indiv_ct4), SEEK_SET)) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (load_and_collapse_incl(bedfile, loadbuf, unfiltered_indiv_ct, &(geno[cur_window_size * founder_ct_mld_long]), founder_ct, founder_info, IS_SET(marker_reverse, window_unfiltered_end))) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(geno[cur_window_size * founder_ct_mld_long])));
	}
	missing_cts[cur_window_size] = ld_process_load(&(geno[cur_window_size * founder_ct_mld_long]), &(geno_masks[cur_window_size * founder_ct_mld_long]), &(geno_mmasks[cur_window_size * founder_ctv]), &(marker_stdevs[cur_window_size]), founder_ct, is_x && (!ignore_x), weighted_x, nonmale_founder_ct, founder_male_include2, nonmale_geno, nonmale_masks, cur_window_size * founder_ct_mld_long);
	cur_window_size++;
      }
      if (cur_window_size > prev_end) {
	start_arr[cur_window_size] = window_unfiltered_end;
      }
    }
    uii = get_marker_chrom(chrom_info_ptr, window_unfiltered_start - 1);
    putchar('\r');
    sprintf(logbuf, "Pruned %" PRIuPTR " variant%s from chromosome %u, leaving %" PRIuPTR ".\n", cur_exclude_ct, (cur_exclude_ct == 1)? "" : "s", uii, chrom_info_ptr->chrom_end[uii] - chrom_info_ptr->chrom_start[uii] - cur_exclude_ct);
    logprintb();
    tot_exclude_ct += cur_exclude_ct;

    // advance chromosomes as necessary
    window_unfiltered_start = ld_prune_next_valid_chrom_start(pruned_arr, window_unfiltered_start, chrom_info_ptr, unfiltered_marker_ct);
  } while (window_unfiltered_start < unfiltered_marker_ct);

  sprintf(logbuf, "Pruning complete.  %u of %" PRIuPTR " variants removed.\n", tot_exclude_ct, marker_ct);
  logprintb();
  strcpy(outname_end, ".prune.in");
  if (fopen_checked(&outfile_in, outname, "w")) {
    goto ld_prune_ret_OPEN_FAIL;
  }
  strcpy(outname_end, ".prune.out");
  if (fopen_checked(&outfile_out, outname, "w")) {
    goto ld_prune_ret_OPEN_FAIL;
  }
  marker_unfiltered_idx = 0;
  marker_idx = 0;
  pct = 1;
  uii = 0;
  for (cur_chrom = 1; cur_chrom <= chrom_info_ptr->max_code; cur_chrom++) {
    if (!IS_SET(chrom_info_ptr->chrom_mask, cur_chrom)) {
      continue;
    }
    if (chrom_info_ptr->chrom_end[cur_chrom]) {
      uii += chrom_info_ptr->chrom_end[cur_chrom] - chrom_info_ptr->chrom_start[cur_chrom];
    }
  }
  pct_thresh = ((int64_t)pct * uii) / 100;
  for (cur_chrom = 1; cur_chrom <= chrom_info_ptr->max_code; cur_chrom++) {
    chrom_end = chrom_info_ptr->chrom_end[cur_chrom];
    if (!chrom_end) {
      continue;
    }
    marker_unfiltered_idx = chrom_info_ptr->chrom_start[cur_chrom];
    for (; marker_unfiltered_idx < chrom_end; marker_unfiltered_idx++) {
      if (!IS_SET(marker_exclude, marker_unfiltered_idx)) {
	sptr = &(marker_ids[marker_unfiltered_idx * max_marker_id_len]);
	fptr = IS_SET(pruned_arr, marker_unfiltered_idx)? outfile_out : outfile_in;
	fwrite(sptr, 1, strlen(sptr), fptr);
	if (putc_checked('\n', fptr)) {
	  goto ld_prune_ret_WRITE_FAIL;
	}
      }
      marker_idx++;
      if (marker_idx == pct_thresh) {
	printf("\rWriting... %d%%", pct);
	fflush(stdout);
	pct = ((int64_t)marker_idx * 100) / uii + 1;
        pct_thresh = ((int64_t)pct * uii) / 100;
      }
    }
  }
  if (fclose_null(&outfile_in)) {
    goto ld_prune_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile_out)) {
    goto ld_prune_ret_WRITE_FAIL;
  }
  *outname_end = '\0';
  putchar('\r');
  sprintf(logbuf, "Marker lists written to %s.prune.in and %s.prune.out.\n", outname, outname);
  logprintb();

  while (0) {
  ld_prune_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_prune_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ld_prune_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_prune_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  ld_prune_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 ld_prune_ret_1:
  fclose_cond(outfile_in);
  fclose_cond(outfile_out);
  wkspace_reset(wkspace_mark);
  return retval;
}

/*
uint32_t ld_process_load2(uintptr_t* geno_buf, uintptr_t* mask_buf, uintptr_t* missing_buf, double* marker_stdev_ptr, uint32_t founder_ct, uint32_t is_x, uintptr_t* founder_male_include2) {
  // ld_process_load(), except no missing_buf[] to conserve memory (and no
  // --ld-xchr 3 support yet)
  uintptr_t* geno_ptr = geno_buf;
  uintptr_t founder_ctl2 = (founder_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t* geno_end = &(geno_buf[founder_ctl2]);
  uintptr_t* mask_buf_ptr = mask_buf;
  uintptr_t new_missing = 0;
  uint32_t missing_bit_offset = 0;
  uint32_t ssq = 0;
  int32_t sum = -founder_ct;
  uintptr_t* nm_mask_ptr;
  double non_missing_recip;
  uintptr_t cur_geno;
  uintptr_t shifted_masked_geno;
  uintptr_t new_geno;
  uintptr_t new_mask;
  uint32_t missing_ct = 0;
  while (1) {
    cur_geno = *geno_ptr;
    shifted_masked_geno = (cur_geno >> 1) & FIVEMASK;
    new_geno = cur_geno - shifted_masked_geno;
    *geno_ptr++ = new_geno;
    new_mask = (((~cur_geno) & FIVEMASK) | shifted_masked_geno) * 3;
    *mask_buf_ptr++ = new_mask;
    if (geno_ptr == geno_end) {
      break;
    }
  }
  if (is_x) {
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
  ssq += popcount2_long(new_geno << (BITCT - 2 * (1 + ((founder_ct - 1) % BITCT2))));
  if (founder_ct % BITCT2) {
    mask_buf[founder_ct / BITCT2] &= (ONELU << (2 * (founder_ct % BITCT2))) - ONELU;
  }
  missing_ct = founder_ct - (popcount_longs(mask_buf, 0, (founder_ct + (BITCT - 1)) / BITCT2) / 2);
  non_missing_recip = 1.0 / (founder_ct - missing_ct);
  *marker_stdev_ptr = non_missing_recip * sqrt(((int64_t)((uint64_t)ssq)) * (founder_ct - missing_ct) - ((int64_t)sum) * sum);
  return missing_ct;
}

int32_t ld_report_square(pthread_t* threads, Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, uintptr_t* sex_male, uintptr_t founder_ct, uintptr_t* founder_include2, uintptr_t* founder_male_include2, char* outname, uint32_t hh_exists) {
  FILE* outfile = NULL;
  gzFile gz_outfile = NULL;
  uint32_t ld_modifier = ldip->modifier;
  uint32_t is_binary = ld_modifier & LD_MATRIX_BIN;
  uint32_t output_gz = ld_modifier & LD_REPORT_GZ;
  uint32_t is_r2 = ld_modifier & LD_R2;
  uint32_t ignore_x = (ldip->modifier / LD_IGNORE_X) & 1;
  uintptr_t marker_ctm8 = (marker_ct + 7) & (~(7 * ONELU));
  uintptr_t founder_ctl = (founder_ct + BITCT - 1) / BITCT;
#ifdef __LP64__
  uintptr_t founder_ctv = 2 * ((founder_ct + 127) / 128);
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
  uintptr_t founder_ct_mld_long = founder_ct_mld * (MULTIPLEX_LD / BITCT2);
  int32_t retval = 0;
  uintptr_t* geno = NULL;
  uintptr_t* geno_masks1;
  uintptr_t* geno_masks2;
  double* marker_stdevs1;
  double* marker_stdevs2;
  uint32_t* missing_cts1;
  uint32_t* missing_cts2;
  double* results;
  uintptr_t marker_idx1;
  uintptr_t marker_idx1_end;
  uintptr_t idx1_block_size;
  if (!output_gz) {
    if (fopen_checked(&outfile, outname, is_binary? "wb" : "w")) {
      goto ld_report_square_ret_OPEN_FAIL;
    }
  } else {
    if (gzopen_checked(&gz_outfile, outname, "wb")) {
      goto ld_report_square_ret_OPEN_FAIL;
    }
  }
  if (g_thread_ct > 1) {
    logprint("Note: --r/--r2 square multithreading hasn't been implemented yet.\n");
  }
  marker_idx1 = (((uint64_t)parallel_idx) * marker_ct) / parallel_tot;
  marker_idx1_end = (((uint64_t)(parallel_idx + 1)) * marker_ct) / parallel_tot;
  // claim up to half of memory with idx1 bufs; each marker costs
  //   founder_ct_mld_long * sizeof(intptr_t) for genotype buffer
  // + founder_ct_mld_long * sizeof(intptr_t) for missing mask buffer
  // + 
  // + marker_ctm8 * sizeof(double) for results buffer
  // small idx1 buffer, with missing1_cts
  // large idx2 buffer, with missing2_cts
  // results: small * all * 8
  do {
    if (idx1_block_size > marker_idx1_end - marker_idx1) {
      idx1_block_size = marker_idx1_end - marker_idx1;
    }
    // ld_process_load(&(geno[ulii * founder_ct_mld_long]), &(geno_masks[ulii * founder_ct_mld_long]), &(geno_mmasks[ulii * founder_ctv]), &(marker_stdevs[ulii]), founder_ct, is_x && (!ignore_x), weighted_x, nonmale_founder_ct, founder_male_include2, nonmale_geno, nonmale_masks, ulii * founder_ct_mld_long);
    
    marker_idx1 += idx1_block_size;
  } while (marker_idx1 < marker_idx1_end);
  // ...
  if (!output_gz) {
    if (fclose_null(&outfile)) {
      goto ld_report_square_ret_WRITE_FAIL;
    }
  } else {
    gzclose(gz_outfile);
    gz_outfile = NULL;
  }
  while (0) {
  ld_report_square_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ld_report_square_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_report_square_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  gzclose_cond(gz_outfile);
  // parent will free memory
  return retval;
}
*/

int32_t ld_report(pthread_t* threads, Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t unfiltered_indiv_ct, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  logprint("Error: --r/--r2 is currently under development.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
  /*
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl2 = 2 * ((unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  uintptr_t founder_ct = popcount_longs(founder_info, 0, unfiltered_indiv_ctl2 / 2);
  uintptr_t* founder_include2 = NULL;
  uintptr_t* founder_male_include2 = NULL;
  uint32_t ld_modifier = ldip->modifier;
  uint32_t is_binary = ld_modifier & LD_MATRIX_BIN;
  uint32_t output_gz = ld_modifier & LD_REPORT_GZ;
  uint32_t is_r2 = ld_modifier & LD_R2;
  char* bufptr = memcpyl3a(outname_end, ".ld");
  int32_t retval = 0;
  if (!founder_ct) {
    sprintf(logbuf, "Warning: Skiping --r%s since there are no founders.\n", is_r2? "2" : "");
    logprintb();
    goto ld_report_ret_1;
  }
  if ((marker_ct > 400000) && (!(ld_modifier & LD_YES_REALLY)) && (parallel_tot == 1) && ((ld_modifier & LD_MATRIX_SHAPEMASK) || ((ld_modifier & LD_INTER_CHR) && (!ldip->snpstr) && (!ldip->snps_rl.name_ct) && ((!is_r2) || (ldip->window_r2 == 0.0))))) {
    logprint("Error: Gigantic (over 400k sites) --r/--r2 unfiltered, non-parallel\ncomputation.  Rerun with the 'yes-really' modifier if you are SURE you have\nenough hard drive space and want to do this.\n");
    goto ld_report_ret_INVALID_CMDLINE;
  }
  if (ld_modifier & LD_SINGLE_PREC) {
    // this will be little more than a copy-and-search-replace job when the
    // other flags are finished
    logprint("Error: --r/--r2 'single-prec' has not been implemented yet.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto ld_report_ret_1;
  }
  if (output_gz) {
    logprint("Error: --r/--r2 'gz' is not implemented yet.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto ld_report_ret_1;
  }
  if (!(ld_modifier & LD_MATRIX_SQ)) {
    logprint("Error: --r/--r2 only supports 'square' output for now.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto ld_report_ret_1;
  }
  if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, founder_ct, XMHH_EXISTS | hh_exists, 1, founder_info, sex_male, &founder_include2, &founder_male_include2)) {
    goto ld_report_ret_NOMEM;
  }
  if (parallel_tot > 1) {
    *bufptr++ = '.';
    bufptr = uint32_write(bufptr, parallel_idx + 1);
  }
  if (output_gz) {
    bufptr = memcpyl3a(bufptr, ".gz");
  } else if (is_binary) {
    bufptr = memcpya(bufptr, ".bin", 4);
  }
  *bufptr = '\0';
  if (ld_modifier & LD_MATRIX_SQ) {
    retval = ld_report_square(threads, ldip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, chrom_info_ptr, unfiltered_indiv_ct, founder_info, parallel_idx, parallel_tot, sex_male, founder_ct, founder_include2, founder_male_include2, outname, hh_exists);
  } else {
    // ...
  }
  if (retval) {
    goto ld_report_ret_1;
  }
  // great success
  while (0) {
  ld_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_report_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 ld_report_ret_1:
  wkspace_reset(wkspace_mark);
  return retval;
  */
}
