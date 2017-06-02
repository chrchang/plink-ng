#ifndef __PLINK2_RANDOM_H__
#define __PLINK2_RANDOM_H__

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


#include "plink2_common.h"
#include "SFMT.h"

#ifdef __cplusplus
namespace plink2 {
#endif

extern sfmt_t g_sfmt;

HEADER_INLINE double rand_unif(sfmt_t* sfmtp) {
  return (sfmt_genrand_uint32(sfmtp) + 0.5) * kRecip2m32;
}

double rand_normal(sfmt_t* sfmtp, double* secondval_ptr);

extern sfmt_t** g_sfmtp_arr;

boolerr_t bigstack_init_sfmtp(uint32_t thread_ct, uint32_t use_main_sfmt_as_element_zero);

pglerr_t fill_gaussian_darray(uintptr_t entry_pair_ct, uint32_t thread_ct, double* darray);

#ifdef __cplusplus
} // namespace plink2
#endif
 
#endif // __PLINK2_RANDOM_H__
