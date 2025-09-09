#ifndef __PLINK2_FASTA_H__
#define __PLINK2_FASTA_H__

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

#include "include/pgenlib_misc.h"
#include "include/plink2_base.h"
#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfFa0,
  kfFaRefFrom = (1 << 0),
  kfFaRefFromForce = (1 << 1),
  kfFaNormalize = (1 << 2),
  kfFaNormalizeList = (1 << 3),
  kfFaNormalizeAdjustOverlappingDeletions = (1 << 4),
  kfFaScrapeLengths = (1 << 5)
FLAGSET_DEF_END(FaFlags);

PglErr ProcessFa(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const ChrInfo* cip, const char* fname, uint32_t max_allele_ct, uint32_t max_allele_slen, FaFlags flags, uint32_t output_missing_geno_code, uint32_t max_thread_ct, UnsortedVar* vpos_sortstatusp, uint32_t* variant_bps, const char** allele_storage, AlleleCode* allele_permute, uintptr_t* nonref_flags, uint32_t* contig_lens, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_FASTA_H__
