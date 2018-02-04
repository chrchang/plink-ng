// This library is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.


#include "plink2_random.h"

#ifdef __cplusplus
namespace plink2 {
#endif

sfmt_t g_sfmt;

double rand_normal(sfmt_t* sfmtp, double* secondval_ptr) {
  // Box-Muller.  try changing this to e.g. ziggurat if it's ever a serious
  // bottleneck.
  const double dxx = sqrt(-2 * log(rand_unif(sfmtp)));
  const double dyy = (2 * kPi) * rand_unif(sfmtp);
  *secondval_ptr = dxx * cos(dyy);
  return dxx * sin(dyy);
}

sfmt_t** g_sfmtp_arr;

boolerr_t bigstack_init_sfmtp(uint32_t thread_ct, uint32_t use_main_sfmt_as_element_zero) {
  g_sfmtp_arr = S_CAST(sfmt_t**, bigstack_alloc(thread_ct * sizeof(intptr_t)));
  if (!g_sfmtp_arr) {
    return 1;
  }
  if (use_main_sfmt_as_element_zero) {
    g_sfmtp_arr[0] = &g_sfmt;
  }
  if (thread_ct > use_main_sfmt_as_element_zero) {
    uint32_t uibuf[4];
    for (uint32_t tidx = use_main_sfmt_as_element_zero; tidx < thread_ct; ++tidx) {
      g_sfmtp_arr[tidx] = S_CAST(sfmt_t*, bigstack_alloc(sizeof(sfmt_t)));
      if (!g_sfmtp_arr[tidx]) {
        return 1;
      }
      for (uint32_t uii = 0; uii < 4; ++uii) {
        uibuf[uii] = sfmt_genrand_uint32(&g_sfmt);
      }
      sfmt_init_by_array(g_sfmtp_arr[tidx], uibuf, 4);
    }
  }
  return 0;
}

// multithread globals
static double* g_darray = nullptr;
static uint32_t g_calc_thread_ct = 0;
static uintptr_t g_entry_pair_ct = 0;

THREAD_FUNC_DECL fill_gaussian_darray_thread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uintptr_t entry_pair_ct = g_entry_pair_ct;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uintptr_t idx_start = (tidx * entry_pair_ct) / calc_thread_ct;
  uintptr_t idx_ct = (((tidx + 1) * entry_pair_ct) / calc_thread_ct) - idx_start;
  double* darray_iter = &(g_darray[idx_start * 2]);
  for (uintptr_t ulii = 0; ulii < idx_ct; ++ulii) {
    double dxx;
    *darray_iter++ = rand_normal(sfmtp, &dxx);
    *darray_iter++ = dxx;
  }
  THREAD_RETURN;
}

pglerr_t fill_gaussian_darray(uintptr_t entry_pair_ct, uint32_t thread_ct, double* darray) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uintptr_t max_useful_thread_ct = DIV_UP(entry_pair_ct, 262144);
    if (thread_ct > max_useful_thread_ct) {
      thread_ct = max_useful_thread_ct;
    }
    pthread_t* threads;
    if (bigstack_init_sfmtp(thread_ct, 1) ||
        bigstack_alloc_thread(thread_ct, &threads)) {
      goto fill_gaussian_darray_ret_NOMEM;
    }
    g_darray = darray;
    g_entry_pair_ct = entry_pair_ct;
    g_calc_thread_ct = thread_ct;
    if (spawn_threads(fill_gaussian_darray_thread, thread_ct, threads)) {
      goto fill_gaussian_darray_ret_THREAD_CREATE_FAIL;
    }
    fill_gaussian_darray_thread(S_CAST(void*, 0));
    join_threads(thread_ct, threads);
  }
  while (0) {
  fill_gaussian_darray_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  fill_gaussian_darray_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}


THREAD_FUNC_DECL randomize_bigstack_thread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint64_t bigstack_int64_ct = S_CAST(uintptr_t, g_bigstack_end - g_bigstack_base) / sizeof(int64_t);
  uint64_t* bigstack_int64 = R_CAST(uint64_t*, g_bigstack_base);
  assert(bigstack_int64_ct >= calc_thread_ct);
  const uint64_t start_idx = round_down_pow2((tidx * bigstack_int64_ct) / calc_thread_ct, kInt64PerCacheline);
  uint64_t end_idx = ((tidx + 1) * bigstack_int64_ct) / calc_thread_ct;
  if (tidx + 1 != calc_thread_ct) {
    end_idx = round_down_pow2(end_idx, kInt64PerCacheline);
  }
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  for (uintptr_t ulii = start_idx; ulii < end_idx; ++ulii) {
    bigstack_int64[ulii] = sfmt_genrand_uint64(sfmtp);
  }
  THREAD_RETURN;
}

pglerr_t randomize_bigstack(uint32_t thread_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (thread_ct > 16) {
      thread_ct = 16;
    }
    if (bigstack_init_sfmtp(thread_ct, 1)) {
      goto randomize_bigstack_ret_NOMEM;
    }
    g_calc_thread_ct = thread_ct;
    pthread_t threads[16];
    if (spawn_threads(randomize_bigstack_thread, thread_ct, threads)) {
      goto randomize_bigstack_ret_THREAD_CREATE_FAIL;
    }
    randomize_bigstack_thread(S_CAST(void*, 0));
    join_threads(thread_ct, threads);
    // now ensure the bytes reserved by bigstack_init_sfmtp() are also properly
    // randomized (some of them already are, but there are gaps)
    uint64_t* initial_segment_end = R_CAST(uint64_t*, g_bigstack_base);
    for (uint64_t* initial_segment_iter = R_CAST(uint64_t*, bigstack_mark); initial_segment_iter != initial_segment_end; ++initial_segment_iter) {
      *initial_segment_iter = sfmt_genrand_uint64(&g_sfmt);
    }
  }
  while (0) {
  randomize_bigstack_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  randomize_bigstack_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

void generate_perm1_interleaved(uint32_t tot_bit_ct, uint32_t set_bit_ct, uintptr_t perm_start_idx, uintptr_t perm_end_idx, uintptr_t* perm_buf, sfmt_t* sfmtp) {
  assert(tot_bit_ct > 1);
  const uintptr_t tot_bit_ctl = BITCT_TO_WORDCT(tot_bit_ct);
  const uintptr_t perm_ct = perm_end_idx - perm_start_idx;
  const uint32_t tot_quotient = 0x100000000LLU / tot_bit_ct;
  const uint32_t upper_bound = tot_bit_ct * tot_quotient - 1;
  uint32_t totq_preshift;
  uint64_t totq_magic;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  // seeing as how we're gonna divide by the same number a billion times or so,
  // it just might be worth optimizing that division...
  magic_num(tot_quotient, &totq_magic, &totq_preshift, &totq_postshift, &totq_incr);
  if (set_bit_ct * 2 < tot_bit_ct) {
    for (uintptr_t widx = 0; widx < tot_bit_ctl; ++widx) {
      fill_ulong_zero(perm_ct, &(perm_buf[perm_start_idx + (widx * perm_end_idx)]));
    }
    for (uintptr_t perm_idx = perm_start_idx; perm_idx < perm_end_idx; ++perm_idx) {
      uintptr_t* pbptr = &(perm_buf[perm_idx]);
      for (uint32_t num_set = 0; num_set < set_bit_ct; ++num_set) {
        uintptr_t widx;
        uintptr_t lowbits;
	do {
          uint32_t urand;
	  do {
	    urand = sfmt_genrand_uint32(sfmtp);
	  } while (urand > upper_bound);
	  // this is identical to lowbits = urand / tot_quotient
	  lowbits = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	  widx = lowbits / kBitsPerWord;
	  lowbits &= (kBitsPerWord - 1);
	} while ((pbptr[widx * perm_end_idx] >> lowbits) & 1);
	pbptr[widx * perm_end_idx] |= (k1LU << lowbits);
      }
    }
  } else {
    for (uintptr_t widx = 0; widx < tot_bit_ctl; ++widx) {
      fill_ulong_one(perm_ct, &(perm_buf[perm_start_idx + (widx * perm_end_idx)]));
    }
    // "set" has reversed meaning here
    set_bit_ct = tot_bit_ct - set_bit_ct;
    for (uintptr_t perm_idx = perm_start_idx; perm_idx < perm_end_idx; ++perm_idx) {
      uintptr_t* pbptr = &(perm_buf[perm_idx]);
      for (uint32_t num_set = 0; num_set < set_bit_ct; num_set++) {
        uintptr_t widx;
        uintptr_t lowbits;
	do {
          uint32_t urand;
	  do {
	    urand = sfmt_genrand_uint32(sfmtp);
	  } while (urand > upper_bound);
	  lowbits = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	  widx = lowbits / kBitsPerWord;
	  lowbits &= (kBitsPerWord - 1);
	} while (!((pbptr[widx * perm_end_idx] >> lowbits) & 1));
	pbptr[widx * perm_end_idx] &= ~(k1LU << lowbits);
      }
    }
    const uint32_t remaining_bit_ct = tot_bit_ct % kBitsPerWord;
    if (remaining_bit_ct) {
      const uintptr_t final_mask = (~k0LU) >> (kBitsPerWord - remaining_bit_ct);
      uintptr_t* pbptr = &(perm_buf[(tot_bit_ctl - 1) * perm_end_idx + perm_start_idx]);
      for (uintptr_t perm_idx = perm_start_idx; perm_idx < perm_end_idx; ++perm_idx) {
	*pbptr &= final_mask;
	pbptr++;
      }
    }
  }
}

#ifdef __cplusplus
}
#endif
