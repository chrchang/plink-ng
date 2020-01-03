#ifndef __PLINK2_RANDOM_H__
#define __PLINK2_RANDOM_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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


#include "plink2_cmdline.h"
#include "include/SFMT.h"

#ifdef __cplusplus
namespace plink2 {
#endif

HEADER_INLINE double RandUnif(sfmt_t* sfmtp) {
  return (sfmt_genrand_uint32(sfmtp) + 0.5) * kRecip2m32;
}

double RandNormal(sfmt_t* sfmtp, double* secondval_ptr);

BoolErr InitAllocSfmtpArr(uint32_t thread_ct, uint32_t use_main_sfmt_as_element_zero, sfmt_t* sfmtp, sfmt_t*** sfmtp_arrp);

PglErr FillGaussianDArr(uintptr_t entry_pair_ct, uint32_t thread_ct, sfmt_t* sfmtp, double* darray);

PglErr RandomizeBigstack(uint32_t thread_ct, sfmt_t* sfmtp);

// might not need plink 1.9-style interleaving, but I'll postpone ripping that
// out.
// currently requires tot_bit_ct > 1, due to quotient operation.
void GeneratePerm1Interleaved(uint32_t tot_bit_ct, uint32_t set_bit_ct, uintptr_t perm_start_idx, uintptr_t perm_end_idx, uintptr_t* perm_buf, sfmt_t* sfmtp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_RANDOM_H__
