// This library is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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
  g_sfmtp_arr = (sfmt_t**)bigstack_alloc(thread_ct * sizeof(intptr_t));
  if (!g_sfmtp_arr) {
    return 1;
  }
  if (use_main_sfmt_as_element_zero) {
    g_sfmtp_arr[0] = &g_sfmt;
  }
  if (thread_ct > use_main_sfmt_as_element_zero) {
    uint32_t uibuf[4];
    for (uint32_t tidx = use_main_sfmt_as_element_zero; tidx < thread_ct; ++tidx) {
      g_sfmtp_arr[tidx] = (sfmt_t*)bigstack_alloc(sizeof(sfmt_t));
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
  const uintptr_t tidx = (uintptr_t)arg;
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
    fill_gaussian_darray_thread((void*)0);
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

#ifdef __cplusplus
}
#endif
