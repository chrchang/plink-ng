#include "plink_ld.h"
#include "plink_matrix.h"
#include "pigz.h"

#define MULTIPLEX_LD 1920
#define MULTIPLEX_2LD (MULTIPLEX_LD * 2)

void ld_epi_init(Ld_info* ldip, Epi_info* epi_ip) {
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
  epi_ip->modifier = 0;
  epi_ip->case_only_gap = 1000000;
  epi_ip->epi1 = 0.0001;
  epi_ip->epi2 = 0.05;
  epi_ip->twolocus_mkr1 = NULL;
  epi_ip->twolocus_mkr2 = NULL;
}

void ld_epi_cleanup(Ld_info* ldip, Epi_info* epi_ip) {
  free_cond(ldip->snpstr);
  free_range_list(&(ldip->snps_rl));
  free_cond(epi_ip->twolocus_mkr1);
  free_cond(epi_ip->twolocus_mkr2);
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
  // where N is the number of individuals processed after applying the
  // missingness masks indicated by the subscripts.
  //
  // Calculation of terms [1]-[4] are  [0] term currently proceeds as follows:
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
  __uni16 acc;
  __uni16 acc1;
  __uni16 acc2;
  __uni16 acc11;
  __uni16 acc22;
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
    // postponed for sum11 and sum12 because the individual terms are 0/1
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
    //   return_vals[1] += N_y + \sum_i x_i
    //   return_vals[2] += N_x + \sum_i y_i
    //   return_vals[3] += N_y - \sum_i x_i^2
    //   return_vals[4] += N_x - \sum_i y_i^2
    // where N is the number of individuals processed after applying the
    // missingness masks indicated by the subscripts.  The [0] calculation
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
#endif // __LP64__

uint32_t ld_process_load(uintptr_t* geno_buf, uintptr_t* mask_buf, uintptr_t* missing_buf, uint32_t* missing_ct_ptr, uint32_t founder_ct, uint32_t is_x, uint32_t weighted_x, uint32_t nonmale_founder_ct, uintptr_t* founder_male_include2, uintptr_t* nonmale_geno, uintptr_t* nonmale_masks, uintptr_t nonmale_offset) {
  uintptr_t* geno_ptr = geno_buf;
  uintptr_t founder_ctl2 = (founder_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t* geno_end = &(geno_buf[founder_ctl2]);
  uintptr_t* mask_buf_ptr = mask_buf;
  uintptr_t* missing_ptr = missing_buf;
  uintptr_t new_missing = 0;
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
    missing_ct += founder_ct - (popcount_longs(nonmale_masks, 0, founder_ctl2) / 2);
    founder_ct *= 2;
  }
  *missing_ct_ptr = missing_ct;
  return (((int64_t)((uint64_t)ssq)) * (founder_ct - missing_ct) - ((int64_t)sum) * sum)? 1 : 0;
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
  uintptr_t founder_ct_192_long = founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2) + founder_ct_mld_rem * (192 / BITCT2);
  uint32_t founder_trail_ct = founder_ct_192_long - founder_ctl * 2;
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
  uintptr_t* loadbuf;
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
  char* sptr;
  FILE* fptr;
  __CLPK_integer window_rem_li;
  __CLPK_integer old_window_rem_li;
  uint32_t window_rem;
  double prune_ld_thresh;
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
    nonmale_founder_ct = founder_ct - popcount01_longs(founder_male_include2, 0, founder_ctl);
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
    prune_ld_thresh = ld_last_param;
  } else {
    // r, not r2, in this case
    prune_ld_thresh = 0.999999;
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
      wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&geno, ulii * founder_ct_192_long * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&geno_masks, ulii * founder_ct_192_long * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&geno_mmasks, ulii * founder_ctv * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&missing_cts, ulii * sizeof(int32_t))) {
    goto ld_prune_ret_NOMEM;
  }
  if (weighted_x) {
    if (wkspace_alloc_ul_checked(&nonmale_geno, ulii * founder_ct_192_long * sizeof(intptr_t)) ||
        wkspace_alloc_ul_checked(&nonmale_masks, ulii * founder_ct_192_long * sizeof(intptr_t))) {
      goto ld_prune_ret_NOMEM;
    }
  }
  for (ulii = 1; ulii <= window_max; ulii++) {
    fill_ulong_zero(&(geno[ulii * founder_ct_192_long - founder_trail_ct - 2]), founder_trail_ct + 2);
    fill_ulong_zero(&(geno_masks[ulii * founder_ct_192_long - founder_trail_ct - 2]), founder_trail_ct + 2);
    if (weighted_x) {
      fill_ulong_zero(&(nonmale_geno[ulii * founder_ct_192_long - founder_trail_ct - 2]), founder_trail_ct + 2);
      fill_ulong_zero(&(nonmale_masks[ulii * founder_ct_192_long - founder_trail_ct - 2]), founder_trail_ct + 2);
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
    cur_exclude_ct = 0;
    if (cur_window_size > 1) {
      for (ulii = 0; ulii < (uintptr_t)cur_window_size; ulii++) {
	uii = live_indices[ulii];
	if (fseeko(bedfile, bed_offset + (uii * unfiltered_indiv_ct4), SEEK_SET)) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (load_and_collapse_incl(bedfile, loadbuf, unfiltered_indiv_ct, &(geno[ulii * founder_ct_192_long]), founder_ct, founder_info, IS_SET(marker_reverse, uii))) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(geno[ulii * founder_ct_192_long])));
	}
        if (!ld_process_load(&(geno[ulii * founder_ct_192_long]), &(geno_masks[ulii * founder_ct_192_long]), &(geno_mmasks[ulii * founder_ctv]), &(missing_cts[ulii]), founder_ct, is_x && (!ignore_x), weighted_x, nonmale_founder_ct, founder_male_include2, nonmale_geno, nonmale_masks, ulii * founder_ct_192_long)) {
	  SET_BIT(pruned_arr, uii);
          cur_exclude_ct++;
	}
      }
    }
    pct = 1;
    pct_thresh = window_unfiltered_start + ((int64_t)pct * (chrom_end - chrom_info_ptr->chrom_start[cur_chrom])) / 100;
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
	      mask_var_vec_ptr = &(geno_masks[ujj * founder_ct_192_long]);

	      dp_result[0] = weighted_founder_ct;
	      // reversed from what I initially thought because I'm passing the
	      // ujj-associated buffers before the uii-associated ones.
	      dp_result[1] = -fixed_non_missing_ct;
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
	      dxx = (dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy);
	      if (!pairwise) {
		dxx = cov12 / sqrt(dxx);
		cov_matrix[uii * window_max + ujj] = dxx;
	      } else {
		dxx = (cov12 * cov12) / dxx;
	      }
	      if (dxx > prune_ld_thresh) {
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
	if (load_and_collapse_incl(bedfile, loadbuf, unfiltered_indiv_ct, &(geno[cur_window_size * founder_ct_192_long]), founder_ct, founder_info, IS_SET(marker_reverse, window_unfiltered_end))) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(geno[cur_window_size * founder_ct_192_long])));
	}
	if (!ld_process_load(&(geno[cur_window_size * founder_ct_192_long]), &(geno_masks[cur_window_size * founder_ct_192_long]), &(geno_mmasks[cur_window_size * founder_ctv]), &(missing_cts[cur_window_size]), founder_ct, is_x && (!ignore_x), weighted_x, nonmale_founder_ct, founder_male_include2, nonmale_geno, nonmale_masks, cur_window_size * founder_ct_192_long)) {
	  SET_BIT(pruned_arr, window_unfiltered_end);
	  cur_exclude_ct++;
	}
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

void ld_process_load2(uintptr_t* geno_buf, uintptr_t* mask_buf, uint32_t* missing_ct_ptr, uint32_t founder_ct, uint32_t is_x, uintptr_t* founder_male_include2) {
  // ld_process_load(), except no missing_buf[] to conserve memory (and no
  // --ld-xchr 3 support yet), and no zero-variance check (we just want to
  // dump nans in that case)
  uintptr_t* geno_ptr = geno_buf;
  uintptr_t founder_ctl2 = (founder_ct + (BITCT2 - 1)) / BITCT2;
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
  *missing_ct_ptr = founder_ct - (popcount_longs(mask_buf, 0, (founder_ct + (BITCT - 1)) / BITCT2) / 2);
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
  __uni16 acc;

  while (word12_ct >= 10) {
    word12_ct -= 10;
    acc.vi = _mm_setzero_si128();
    vend1 = &(vptr1[60]);
  ld_missing_ct_intersect_main_loop:
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
    acc.vi = _mm_setzero_si128();
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

// multithread globals
static uintptr_t* g_ld_geno1;
static uintptr_t* g_ld_geno2;
static uintptr_t* g_ld_geno_masks1;
static uintptr_t* g_ld_geno_masks2;
static uint32_t* g_ld_missing_cts1;
static uint32_t* g_ld_missing_cts2;
static double* g_ld_results;
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
static char g_ld_delimiter;
static uint32_t g_ld_plink_maxsnp;
static char** g_ld_marker_allele_ptrs;
static Chrom_info* g_ld_chrom_info_ptr;
static uint32_t* g_ld_marker_pos;

// [4n]: self uidx
// [4n + 1]: uidx of first variant this is paired against
// [4n + 2]: idx of that variant
// [4n + 3]: total number of variants this paired against
// static uint32_t* g_ld_idx1_info;

THREAD_RET_TYPE ld_block_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t block_idx1_start = (tidx * g_ld_idx1_block_size) / g_ld_thread_ct;
  uintptr_t block_idx1_end = ((tidx + 1) * g_ld_idx1_block_size) / g_ld_thread_ct;
  uintptr_t idx2_block_size = g_ld_idx2_block_size;
  uintptr_t idx2_block_start = g_ld_idx2_block_start;
  uintptr_t marker_ctm8 = g_ld_marker_ctm8;
  uintptr_t founder_ct = g_ld_founder_ct;
  uintptr_t founder_ctwd = founder_ct / BITCT2;
  uintptr_t founder_ctwd12 = founder_ctwd / 12;
  uintptr_t founder_ctwd12_rem = founder_ctwd - (12 * founder_ctwd12);
  uintptr_t lshift_last = 2 * ((0x7fffffc0 - founder_ct) % BITCT2);
  uintptr_t founder_ct_192_long = g_ld_founder_ct_192_long;
  uintptr_t* geno1 = g_ld_geno1;
  uintptr_t* geno2 = g_ld_geno2;
  uintptr_t* geno_masks1 = g_ld_geno_masks1;
  uintptr_t* geno_masks2 = g_ld_geno_masks2;
  uint32_t* missing_cts1 = g_ld_missing_cts1;
  uint32_t* missing_cts2 = g_ld_missing_cts2;
  uint32_t founder_ct_mld_m1 = g_ld_founder_ct_mld_m1;
  uint32_t founder_ct_mld_rem = g_ld_founder_ct_mld_rem;
  uint32_t is_r2 = g_ld_is_r2;
  double* results = g_ld_results;

  // todo: add another slot for missing ct and get rid of missing_ct_intersect?
  int32_t dp_result[5];
  double* rptr;
  uintptr_t* geno_fixed_vec_ptr;
  uintptr_t* geno_var_vec_ptr;
  uintptr_t* mask_fixed_vec_ptr;
  uintptr_t* mask_var_vec_ptr;
  uintptr_t block_idx1;
  uintptr_t block_idx2;
  double non_missing_ctd;
  double cov12;
  double dxx;
  double dyy;
  uint32_t fixed_missing_ct;
  uint32_t fixed_non_missing_ct;
  uint32_t non_missing_ct;
  for (block_idx1 = block_idx1_start; block_idx1 < block_idx1_end; block_idx1++) {
    rptr = &(results[block_idx1 * marker_ctm8 + idx2_block_start]);
    fixed_missing_ct = missing_cts1[block_idx1];
    fixed_non_missing_ct = founder_ct - fixed_missing_ct;
    geno_fixed_vec_ptr = &(geno1[block_idx1 * founder_ct_192_long]);
    mask_fixed_vec_ptr = &(geno_masks1[block_idx1 * founder_ct_192_long]);
    for (block_idx2 = 0; block_idx2 < idx2_block_size; block_idx2++, rptr++) {
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
      dxx = (dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy);
      if (!is_r2) {
	dxx = cov12 / sqrt(dxx);
      } else {
	dxx = (cov12 * cov12) / dxx;
      }
      *rptr = dxx;
    }
  }
  THREAD_RETURN;
}

uint32_t ld_square_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t block_size1 = g_ld_idx1_block_size;
  uintptr_t marker_ct = g_ld_marker_ct;
  uintptr_t marker_ctm8 = g_ld_marker_ctm8;
  uintptr_t block_idx1 = g_ld_block_idx1;
  uintptr_t marker_idx = g_ld_idx2_block_start;
  char delimiter = g_ld_delimiter;
  double* results = g_ld_results;
  double* dptr;
  while (block_idx1 < block_size1) {
    dptr = &(results[block_idx1 * marker_ctm8 + marker_idx]);
    while (marker_idx < marker_ct) {
      sptr_cur = double_g_writex(sptr_cur, *dptr++, delimiter);
      marker_idx++;
      if (sptr_cur > readbuf_end) {
	goto ld_square_emitn_ret;
      }
    }
    if (delimiter == '\t') {
      sptr_cur--;
    }
    *sptr_cur++ = '\n';
    marker_idx = 0;
    block_idx1++;
  }
 ld_square_emitn_ret:
  g_ld_block_idx1 = block_idx1;
  g_ld_idx2_block_start = marker_idx;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t ld_report_square(pthread_t* threads, Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, uintptr_t unfiltered_indiv_ct, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, uintptr_t* sex_male, uintptr_t* founder_include2, uintptr_t* founder_male_include2, uintptr_t* loadbuf, char* outname, uint32_t hh_exists) {
  FILE* outfile = NULL;
  uint32_t ld_modifier = ldip->modifier;
  uint32_t is_binary = ld_modifier & LD_MATRIX_BIN;
  uint32_t output_gz = ld_modifier & LD_REPORT_GZ;
  uint32_t ignore_x = (ld_modifier / LD_IGNORE_X) & 1;
  uintptr_t marker_ct = g_ld_marker_ct;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t marker_ctm8 = (marker_ct + 7) & (~(7 * ONELU));
  uintptr_t idx2_rem = marker_ctm8 - marker_ct;
  uintptr_t founder_ct = g_ld_founder_ct;
  uintptr_t founder_ctl = (founder_ct + BITCT - 1) / BITCT;
  uintptr_t founder_ct_192_long = g_ld_founder_ct_192_long;
  uintptr_t marker_uidx_base = next_unset_unsafe(marker_exclude, 0);
  uintptr_t marker_uidx1 = marker_uidx_base;
  uintptr_t marker_idx1_start = (((uint64_t)parallel_idx) * marker_ct) / parallel_tot;
  uintptr_t marker_idx1 = marker_idx1_start;
  uintptr_t marker_idx1_end = (((uint64_t)(parallel_idx + 1)) * marker_ct) / parallel_tot;
  uintptr_t job_size = marker_idx1_end - marker_idx1;
  uintptr_t pct = 1;
  Chrom_info* chrom_info_ptr = g_ld_chrom_info_ptr;
  uintptr_t pct_thresh = job_size / 100;
  uint32_t founder_trail_ct = founder_ct_192_long - founder_ctl * 2;
  uint32_t thread_ct = g_thread_ct;
  uint32_t chrom_fo_idx = 0;
  uint32_t is_haploid = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t not_first_write = 0;
  int32_t retval = 0;
  unsigned char* wkspace_mark2;
  uintptr_t thread_workload;
  uintptr_t cur_idx2_block_size;
  uintptr_t marker_uidx1_tmp;
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

  if (is_binary) {
    if (fopen_checked(&outfile, outname, "wb")) {
      goto ld_report_square_ret_OPEN_FAIL;
    }
  }
  g_ld_marker_ctm8 = marker_ctm8;
  // claim up to half of memory with idx1 bufs; each marker costs
  //   founder_ct_192_long * sizeof(intptr_t) for genotype buffer
  // + founder_ct_192_long * sizeof(intptr_t) for missing mask buffer
  // + sizeof(int32_t) for g_ld_missing_cts1 entry
  // + marker_ctm8 * sizeof(double) for g_ld_results buffer
  // round down to multiple of thread_ct for better workload distribution
  ulii = founder_ct_192_long * 2 * sizeof(intptr_t) + sizeof(int32_t) + marker_ctm8 * sizeof(double);
  idx1_block_size = wkspace_left / (ulii * 2);
  thread_workload = idx1_block_size / thread_ct;
  if (!thread_workload) {
    goto ld_report_square_ret_NOMEM;
  }
  idx1_block_size = thread_workload * thread_ct;
  if (idx1_block_size > job_size) {
    idx1_block_size = job_size;
  }
  g_ld_geno1 = (uintptr_t*)wkspace_alloc(founder_ct_192_long * idx1_block_size * sizeof(intptr_t));
  g_ld_geno_masks1 = (uintptr_t*)wkspace_alloc(founder_ct_192_long * idx1_block_size * sizeof(intptr_t));
  g_ld_missing_cts1 = (uint32_t*)wkspace_alloc(idx1_block_size * sizeof(int32_t));
  if (wkspace_alloc_d_checked(&g_ld_results, marker_ctm8 * idx1_block_size * sizeof(double))) {
    goto ld_report_square_ret_NOMEM;
  }

  // claim the other half with idx2 buffer; all but the g_ld_results buffer
  // cost apply
  ulii -= (marker_ctm8 * sizeof(double) - 4);
  idx2_block_size = (wkspace_left / ulii) & (~(7 * ONELU));
  if (idx2_block_size > marker_ctm8) {
    idx2_block_size = marker_ctm8;
  }
  wkspace_mark2 = wkspace_base;
  while (1) {
    if (!idx2_block_size) {
      goto ld_report_square_ret_NOMEM;
    }
    if (!(wkspace_alloc_ul_checked(&g_ld_geno2, founder_ct_192_long * idx2_block_size * sizeof(intptr_t)) ||
          wkspace_alloc_ul_checked(&g_ld_geno_masks2, founder_ct_192_long * idx2_block_size * sizeof(intptr_t)) ||
          wkspace_alloc_ui_checked(&g_ld_missing_cts2, idx2_block_size * sizeof(int32_t)))) {
      break;
    }
    wkspace_reset(wkspace_mark2);
    idx2_block_size -= 8;
  }
  uljj = founder_trail_ct + 2;
  for (ulii = 1; ulii <= idx1_block_size; ulii++) {
    fill_ulong_zero(&(g_ld_geno1[ulii * founder_ct_192_long - uljj]), uljj);
    fill_ulong_zero(&(g_ld_geno_masks1[ulii * founder_ct_192_long - uljj]), uljj);
  }
  for (ulii = 1; ulii <= idx2_block_size; ulii++) {
    fill_ulong_zero(&(g_ld_geno2[ulii * founder_ct_192_long - uljj]), uljj);
    fill_ulong_zero(&(g_ld_geno_masks2[ulii * founder_ct_192_long - uljj]), uljj);
  }
  if (marker_idx1) {
    marker_uidx1 = jump_forward_unset_unsafe(marker_exclude, marker_uidx1 + 1, marker_idx1);
  }
  sprintf(logbuf, "--r%s square to %s...", g_ld_is_r2? "2" : "", outname);
  logprintb();
  fputs(" 0%", stdout);
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
    marker_uidx1_tmp = marker_uidx1;
    if (fseeko(bedfile, bed_offset + (marker_uidx1 * unfiltered_indiv_ct4), SEEK_SET)) {
      goto ld_report_square_ret_READ_FAIL;
    }
    chrom_end = 0;
    for (block_idx1 = 0; block_idx1 < idx1_block_size; marker_uidx1_tmp++, block_idx1++) {
      if (IS_SET(marker_exclude, marker_uidx1_tmp)) {
        marker_uidx1_tmp = next_unset_ul_unsafe(marker_exclude, marker_uidx1_tmp);
        if (fseeko(bedfile, bed_offset + (marker_uidx1_tmp * unfiltered_indiv_ct4), SEEK_SET)) {
	  goto ld_report_square_ret_READ_FAIL;
	}
      }
      if (marker_uidx1_tmp >= chrom_end) {
        chrom_fo_idx = get_marker_chrom_fo_idx(chrom_info_ptr, marker_uidx1_tmp);
        chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
        is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
	is_x = (((int32_t)chrom_idx) == chrom_info_ptr->x_code)? 1 : 0;
	is_y = (((int32_t)chrom_idx) == chrom_info_ptr->y_code)? 1 : 0;
      }
      if (load_and_collapse_incl(bedfile, loadbuf, unfiltered_indiv_ct, &(g_ld_geno1[block_idx1 * founder_ct_192_long]), founder_ct, founder_info, IS_SET(marker_reverse, marker_uidx1_tmp))) {
	goto ld_report_square_ret_READ_FAIL;
      }
      if (is_haploid && hh_exists) {
	haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(g_ld_geno1[block_idx1 * founder_ct_192_long])));
      }
      ld_process_load2(&(g_ld_geno1[block_idx1 * founder_ct_192_long]), &(g_ld_geno_masks1[block_idx1 * founder_ct_192_long]), &(g_ld_missing_cts1[block_idx1]), founder_ct, is_x && (!ignore_x), founder_male_include2);
    }
    marker_uidx2 = marker_uidx_base;
    marker_idx2 = 0;
    chrom_end = 0;
    if (fseeko(bedfile, bed_offset + (marker_uidx2 * unfiltered_indiv_ct4), SEEK_SET)) {
      goto ld_report_square_ret_READ_FAIL;
    }
    cur_idx2_block_size = idx2_block_size;
    do {
      if (cur_idx2_block_size > marker_ct - marker_idx2) {
	cur_idx2_block_size = marker_ct - marker_idx2;
      }
      for (block_idx2 = 0; block_idx2 < cur_idx2_block_size; marker_uidx2++, block_idx2++) {
	if (IS_SET(marker_exclude, marker_uidx2)) {
          marker_uidx2 = next_unset_ul_unsafe(marker_exclude, marker_uidx2);
	  if (fseeko(bedfile, bed_offset + (marker_uidx2 * unfiltered_indiv_ct4), SEEK_SET)) {
	    goto ld_report_square_ret_READ_FAIL;
	  }
	}
	if (marker_uidx2 >= chrom_end) {
	  chrom_fo_idx = get_marker_chrom_fo_idx(chrom_info_ptr, marker_uidx2);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
	  is_x = (((int32_t)chrom_idx) == chrom_info_ptr->x_code)? 1 : 0;
	  is_y = (((int32_t)chrom_idx) == chrom_info_ptr->y_code)? 1 : 0;
	}
	if (load_and_collapse_incl(bedfile, loadbuf, unfiltered_indiv_ct, &(g_ld_geno2[block_idx2 * founder_ct_192_long]), founder_ct, founder_info, IS_SET(marker_reverse, marker_uidx2))) {
	  goto ld_report_square_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, founder_include2, founder_male_include2, founder_ct, is_x, is_y, (unsigned char*)(&(g_ld_geno2[block_idx2 * founder_ct_192_long])));
	}
	ld_process_load2(&(g_ld_geno2[block_idx2 * founder_ct_192_long]), &(g_ld_geno_masks2[block_idx2 * founder_ct_192_long]), &(g_ld_missing_cts2[block_idx2]), founder_ct, is_x && (!ignore_x), founder_male_include2);
      }
      g_ld_idx2_block_size = cur_idx2_block_size;
      g_ld_idx2_block_start = marker_idx2;
      if (spawn_threads(threads, &ld_block_thread, g_ld_thread_ct)) {
	goto ld_report_square_ret_THREAD_CREATE_FAIL;
      }
      ld_block_thread((void*)0);
      join_threads(threads, g_ld_thread_ct);
      marker_idx2 += cur_idx2_block_size;
    } while (marker_idx2 < marker_ct);
    fputs("\b\b\b\b\b\b\b\b\b\b\bwriting]   \b\b\b", stdout);
    fflush(stdout);
    if (is_binary) {
      if (!idx2_rem) {
	if (fwrite_checked(g_ld_results, marker_ct * idx1_block_size * sizeof(double), outfile)) {
	  goto ld_report_square_ret_WRITE_FAIL;
	}
      } else {
	for (block_idx1 = 0; block_idx1 < idx1_block_size; block_idx1++) {
	  if (fwrite_checked(&(g_ld_results[block_idx1 * marker_ctm8]), cur_idx2_block_size * sizeof(double), outfile)) {
	    goto ld_report_square_ret_WRITE_FAIL;
	  }
	}
      }
    } else {
      g_ld_block_idx1 = 0;
      g_ld_idx2_block_start = 0;
      if (output_gz) {
        parallel_compress(outname, not_first_write, ld_square_emitn);
      } else {
        write_uncompressed(outname, not_first_write, ld_square_emitn);
      }
      not_first_write = 1;
    }
    marker_idx1 += idx1_block_size;
    fputs("\b\b\b\b\b\b\b\b\b\b          \b\b\b\b\b\b\b\b\b\b", stdout);
    if (marker_idx1 >= pct_thresh) {
      if (pct > 10) {
	putchar('\b');
      }
      pct = ((marker_idx1 - marker_idx1_start) * 100LLU) / job_size;
      if (pct < 100) {
	printf("\b\b%" PRIuPTR "%%", pct);
	fflush(stdout);
	pct_thresh = marker_idx1_start + ((++pct) * ((uint64_t)job_size)) / 100;
      }
    }
  } while (marker_idx1 < marker_idx1_end);
  fputs("\b\b\b", stdout);
  logprint(" done.\n");
  if (is_binary) {
    if (fclose_null(&outfile)) {
      goto ld_report_square_ret_WRITE_FAIL;
    }
  }
  while (0) {
  ld_report_square_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_report_square_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ld_report_square_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_report_square_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  ld_report_square_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
  fclose_cond(outfile);
  // trust parent to free memory
  return retval;
}

/*
int32_t ld_report_regular(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, uintptr_t unfiltered_indiv_ct, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, uintptr_t* sex_male, uintptr_t* founder_include2, uintptr_t* founder_male_include2, uintptr_t* loadbuf, char* outname, uint32_t hh_exists) {
  FILE* outfile = NULL;
  gzFile gz_outfile = NULL;
  uint32_t ld_modifier = ldip->modifier;
  uint32_t output_gz = ld_modifier & LD_REPORT_GZ;
  uint32_t ignore_x = (ld_modifier & LD_IGNORE_X) & 1;
  uint32_t is_inter_chr = ld_modifier & LD_INTER_CHR;
  uintptr_t marker_ct = g_ld_marker_ct;
  uintptr_t marker_ct1 = marker_ct;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t founder_ct = g_ld_founder_ct;
  uintptr_t founder_ctl = (founder_ct + BITCT - 1) / BITCT;
  uintptr_t founder_ct_192_long = founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2) + founder_ct_mld_rem * (192 / BITCT2);
  uintptr_t pct = 1;
  uintptr_t marker_idx2_maxw = 1;
  Chrom_info* chrom_info_ptr = g_ld_chrom_info_ptr;
  uintptr_t* marker_exclude_idx1 = marker_exclude;
  uint32_t* marker_pos = g_ld_marker_pos;
  uint32_t founder_trail_ct = founder_ct_192_long - founder_ctl * 2;
  uint32_t idx1_subset = (ldip->snpstr || ldip->snps_rl)? 1 : 0;
  uint32_t window_size = ldip->window_size;
  uint32_t window_bp = ldip->window_bp;
  uint32_t thread_ct = g_thread_ct;
  uint32_t chrom_fo_idx = 0;
  uint32_t is_haploid = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t not_first_write = 0;
  int32_t retval = 0;
  uintptr_t marker_idx1_start;
  uintptr_t marker_idx1;
  uintptr_t marker_idx1_end;
  uintptr_t job_size;
  uintptr_t pct_thresh;
  uintptr_t marker_uidx1;
  uintptr_t marker_uidx2_base;
  uintptr_t marker_uidx2;
  uintptr_t marker_idx2_maxw_m8;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t chrom_size;
  uint32_t cur_marker_pos;
  if (idx1_subset) {
    // todo: allocate marker_exclude mirror if --ld-snp/--ld-snps/--ld-snp-list
    // present
    logprint("Error: --ld-snp, --ld-snps, and --ld-snp-list are not implemented yet.\n");
    return RET_CALC_NOT_YET_SUPPORTED;
    if (!marker_ct1) {
      logprint("Error: No valid variants specified by --ld-snp/--ld-snps/--ld-snp-list.\n");
      goto ld_report_regular_ret_INVALID_CMDLINE;
    }
  }
  if (marker_ct1 < parallel_tot) {
    sprintf(logbuf, "Error: Too few variants in --r%s run for --parallel %u %u.\n", g_ld_is_r2? "2" : "", parallel_idx + 1, parallel_tot);
    logprintb();
    goto ld_report_regular_ret_INVALID_CMDLINE;
  }
  marker_idx1_start = (((uint64_t)parallel_idx) * marker_ct1) / parallel_tot;
  marker_idx1 = marker_idx1_start;
  marker_idx1_end = (((uint64_t)(parallel_idx + 1)) * marker_ct1) / parallel_tot;
  job_size = marker_idx1_end - marker_idx1;
  pct_thresh = job_size / 100;

  if (is_inter_chr) {
    if (idx1_subset) {
      marker_idx2_maxw = marker_ct;
    } else {
      marker_idx2_maxw = marker_ct - 1;
    }
  } else {
    if (idx1_subset) {
      // look backwards and forwards
      if (window_size < 13) {
	marker_idx2_maxw = 2 * window_size - 1;
      } else {
	// ...
      }
    } else if (window_size < 18) {
      // just don't bother scanning in default case
      marker_idx2_maxw = window_size - 1;
    } else {
      for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx) {
	chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
	marker_uidx1 = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
	chrom_size = chrom_end - marker_uidx1 - popcount_bit_idx(marker_exclude, marker_uidx1, chrom_end);
	if (chrom_size <= marker_idx2_maxw) {
	  continue;
	}
	cur_marker_pos = marker_pos[marker_uidx1];
	marker_uidx2 = jump_forward_unset_unsafe(marker_exclude, marker_uidx1, marker_idx2_maxw);
	chrom_size -= ;
	while (1) {
	  // invariant: marker_uidx2 is marker_idx2_maxw sites after
	  // marker_uidx1.
	  if (marker_pos[marker_uidx2] - cur_marker_pos <= window_bp) {
	    marker_idx2_maxw++;
	    if (marker_idx2_maxw + 1 >= window_size) {
	      marker_idx2_maxw = window_size - 1;
	      goto ld_report_regular_maxw_ceil;
	    }
	    if (!(--chrom_size)) {
	      break;
	    }
	  }
	}
      }
    }
  }
 ld_report_regular_maxw_ceil:
  marker_idx2_maxw = (marker_idx2_maxw + 7) & (~(7 * ONELU));
  sprintf(logbuf, "--r%s%s to %s...", g_ld_is_r2? "2" : "", is_inter_chr? " inter-chr" : "", outname);
  logprintb();
  fputs(" 0%", stdout);
  do {
    fputs(" [processing]", stdout);
    fflush(stdout);
    if (idx1_block_size > marker_idx1_end - marker_idx1) {
      // ...
    }
  }
  while (0) {
  ld_report_regular_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_report_regular_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_report_regular_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  // trust parent to free memory
  return retval;
}
*/

int32_t ld_report(pthread_t* threads, Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t unfiltered_indiv_ct, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl2 = 2 * ((unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  uintptr_t founder_ct = popcount_longs(founder_info, 0, unfiltered_indiv_ctl2 / 2);
  uintptr_t* founder_include2 = NULL;
  uintptr_t* founder_male_include2 = NULL;
  uintptr_t founder_ct_mld = (founder_ct + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
  uint32_t founder_ct_mld_m1 = ((uint32_t)founder_ct_mld) - 1;
#ifdef __LP64__
  uint32_t founder_ct_mld_rem = (MULTIPLEX_LD / 192) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 192;
#else
  uint32_t founder_ct_mld_rem = (MULTIPLEX_LD / 48) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 48;
#endif
  uintptr_t founder_ct_192_long = founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2) + founder_ct_mld_rem * (192 / BITCT2);
  uint32_t ld_modifier = ldip->modifier;
  uint32_t is_binary = ld_modifier & LD_MATRIX_BIN;
  uint32_t output_gz = ld_modifier & LD_REPORT_GZ;
  char* bufptr = memcpyl3a(outname_end, ".ld");
  int32_t retval = 0;
  uintptr_t* loadbuf;

  g_ld_founder_ct = founder_ct;
  g_ld_founder_ct_192_long = founder_ct_192_long;
  g_ld_founder_ct_mld_m1 = founder_ct_mld_m1;
  g_ld_founder_ct_mld_rem = founder_ct_mld_rem;
  g_ld_is_r2 = ld_modifier & LD_R2;
  g_ld_marker_ct = marker_ct;
  g_ld_chrom_info_ptr = chrom_info_ptr;
  g_ld_thread_ct = g_thread_ct;
  if (!founder_ct) {
    sprintf(logbuf, "Warning: Skipping --r%s since there are no founders.\n", g_ld_is_r2? "2" : "");
    logprintb();
    goto ld_report_ret_1;
  }
  if ((marker_ct > 400000) && (!(ld_modifier & LD_YES_REALLY)) && (parallel_tot == 1) && ((ld_modifier & LD_MATRIX_SHAPEMASK) || ((ld_modifier & LD_INTER_CHR) && (!ldip->snpstr) && (!ldip->snps_rl.name_ct) && ((!g_ld_is_r2) || (ldip->window_r2 == 0.0))))) {
    logprint("Error: Gigantic (over 400k sites) --r/--r2 unfiltered, non-distributed\ncomputation.  Rerun with the 'yes-really' modifier if you are SURE you have\nenough hard drive space and want to do this.\n");
    goto ld_report_ret_INVALID_CMDLINE;
  }
  if (ld_modifier & LD_SINGLE_PREC) {
    // this will be little more than a copy-and-search-replace job when the
    // other flags are finished
    logprint("Error: --r/--r2 'single-prec' has not been implemented yet.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto ld_report_ret_1;
  }
  if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, founder_ct, XMHH_EXISTS | hh_exists, 1, founder_info, sex_male, &founder_include2, &founder_male_include2)) {
    goto ld_report_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto ld_report_ret_NOMEM;
  }
  if (parallel_tot > 1) {
    *bufptr++ = '.';
    bufptr = uint32_write(bufptr, parallel_idx + 1);
  }
  if (is_binary) {
    bufptr = memcpya(bufptr, ".bin", 4);
  } else {
    g_ld_delimiter = (ld_modifier & LD_MATRIX_SPACES)? ' ' : '\t';
    if (output_gz) {
      bufptr = memcpyl3a(bufptr, ".gz");
    }
  }
  *bufptr = '\0';
  if (ld_modifier & LD_MATRIX_SQ) {
    retval = ld_report_square(threads, ldip, bedfile, bed_offset, unfiltered_marker_ct, marker_exclude, marker_reverse, unfiltered_indiv_ct, founder_info, parallel_idx, parallel_tot, sex_male, founder_include2, founder_male_include2, loadbuf, outname, hh_exists);
  } else if (ld_modifier & (LD_MATRIX_SQ0 | LD_MATRIX_TRI)) {
    logprint("Error: --r/--r2 does not support square0 or triangle output yet.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto ld_report_ret_1;
  } else {
    logprint("Error: --r/--r2 table reports are currently under development.  (Square matrix\noutput is working.)");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto ld_report_ret_1;
    g_ld_plink_maxsnp = plink_maxsnp;
    g_ld_marker_allele_ptrs = marker_allele_ptrs;
    g_ld_marker_pos = marker_pos;
    // retval = ld_report_regular(ldip, bedfile, bed_offset, unfiltered_marker_ct, marker_exclude, marker_reverse, unfiltered_indiv_ct, founder_info, parallel_idx, parallel_tot, sex_male, founder_include2, founder_male_include2, loadbuf, outname, hh_exists);
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
}

int32_t twolocus(Epi_info* epi_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists) {
  logprint("Error: --twolocus is currently under development.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
}

int32_t epistasis_report(pthread_t* threads, Epi_info* epi_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* outname, char* outname_end) {
  logprint("Error: --epistasis and --fast-epistasis are currently under development.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
}
