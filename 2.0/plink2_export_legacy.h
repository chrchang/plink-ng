#ifndef __PLINK2_EXPORT_LEGACY_H__
#define __PLINK2_EXPORT_LEGACY_H__

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
#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "include/plink2_thread.h"
#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// This may belong in a more central location.
typedef struct MTPgenReadCtxStruct {
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const AlleleCode* allele_permute;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  uint32_t sample_ct;

  PgenReader** pgr_ptrs;

  uint32_t* variant_uidx_starts;
  uint32_t cur_block_write_ct;

  uintptr_t* vmaj_readbuf;

  uint64_t err_info;
} MTPgenReadCtx;

THREAD_FUNC_DECL MTPgenReadThread(void* raw_arg);

PglErr ExportIndMajorBed(const uintptr_t* orig_sample_include, const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const AlleleCode* allele_permute, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end);

PglErr ExportTped(const char* outname, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const double* variant_cms, uint32_t sample_ct, uint32_t variant_ct, uint32_t max_allele_slen, char exportf_delim, char lomg_char, PgenReader* simple_pgrp);

PglErr ExportPed(const char* outname, const uintptr_t* orig_sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const char* legacy_output_missing_pheno, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t compound_genotypes, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, char exportf_delim, char lomg_char, PgenFileInfo* pgfip);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_EXPORT_LEGACY_H__
