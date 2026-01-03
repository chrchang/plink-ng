#ifndef __PLINK2_RANDOM_H__
#define __PLINK2_RANDOM_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
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

#include "include/SFMT.h"
#include "include/plink2_base.h"
#include "include/plink2_string.h"

#ifdef __cplusplus
namespace plink2 {
#endif

HEADER_INLINE double RandUnif(sfmt_t* sfmtp) {
  return (sfmt_genrand_uint32(sfmtp) + 0.5) * k2m32;
}

double RandNormal(sfmt_t* sfmtp, double* secondval_ptr);

BoolErr InitAllocSfmtpArr(uint32_t thread_ct, uint32_t use_main_sfmt_as_element_zero, sfmt_t* sfmtp, sfmt_t*** sfmtp_arrp);

PglErr FillGaussianDArr(uintptr_t entry_pair_ct, uint32_t thread_ct, sfmt_t* sfmtp, double* darray);

PglErr RandomizeBigstack(uint32_t thread_ct, sfmt_t* sfmtp);

HEADER_INLINE uint32_t RandU32(uint32_t range, sfmt_t* sfmtp) {
  // [0, range)
  // https://lemire.me/blog/2016/06/30/fast-random-shuffling/
  uint64_t random32bit = sfmt_genrand_uint32(sfmtp);
  uint64_t multiresult = random32bit * range;
  uint32_t leftover = S_CAST(uint32_t, multiresult);
  if (leftover < range) {
    const uint32_t threshold = -range % range;
    while (leftover < threshold) {
      random32bit = sfmt_genrand_uint32(sfmtp);
      multiresult = random32bit * range;
      leftover = S_CAST(uint32_t, multiresult);
    }
  }
  return multiresult >> 32;
}

void PermuteU32(uint32_t entry_ct, sfmt_t* sfmtp, uint32_t* u32arr);

// might not need plink 1.9-style interleaving, but I'll postpone ripping that
// out.
// currently requires tot_bit_ct > 1, due to quotient operation.
void GeneratePerm1Interleaved(uint32_t tot_bit_ct, uint32_t set_bit_ct, uintptr_t perm_start_idx, uintptr_t perm_end_idx, uintptr_t* perm_buf, sfmt_t* sfmtp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_RANDOM_H__
