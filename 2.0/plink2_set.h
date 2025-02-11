#ifndef __PLINK2_SET_H__
#define __PLINK2_SET_H__

// This file is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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

#include "include/plink2_base.h"
#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr ExtractExcludeRange(const char* fnames, const ChrInfo* cip, const uint32_t* variant_bps, uint32_t raw_variant_ct, VfilterType vft, uint32_t zero_based, uint32_t bed_border_bp, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

// setdef is guaranteed to be 16-byte aligned.
// Encoding (same as plink 1.9):
// [0]: either number of ranges in set, or UINT32_MAX
// If [0] is not UINT32_MAX:
//   [2k+1], [2k+2]: start and end of 0-based range #k (half-open); sorted
// If [0] is UINT32_MAX:
//   [1]: offset of first bit (always divisible by 128)
//   [2]: number of bits (divisible by 128 unless a variant in the last block
//        is included)
//   [3]: 1 if all out-of-bounds bits are set, 0 otherwise (other flags may be
//        added later)
// Empty sets are always stored in the first format, with [0] == 0.

uint32_t IntervalInSetdef(const uint32_t* setdef, uint32_t variant_uidx_start, uint32_t variant_uidx_end);

PglErr LoadAndSortIntervalBed(const char* fname, const ChrInfo* cip, const char* sorted_subset_ids, uint32_t zero_based, uint32_t border_extend, uintptr_t subset_ct, uintptr_t max_subset_id_blen, uint32_t max_thread_ct, uintptr_t* gene_ct_ptr, char** gene_names_ptr, uintptr_t* max_gene_id_blen_ptr, uintptr_t** chr_bounds_ptr, uint32_t*** genedefs_ptr, uintptr_t* chr_max_gene_ct_ptr);

#ifdef __cplusplus
}
#endif

#endif  // __PLINK2_SET_H__
