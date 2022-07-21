#ifndef __PLINK2_GLM_H__
#define __PLINK2_GLM_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2022 Shaun Purcell,
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


#include "plink2_adjust.h"
#include "plink2_glm_shared.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitGlm(GlmInfo* glm_info_ptr);

void CleanupGlm(GlmInfo* glm_info_ptr);

PglErr GlmMain(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const AdjustInfo* adjust_info_ptr, const APerm* aperm_ptr, const char* local_covar_fname, const char* local_pvar_fname, const char* local_psam_fname, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t orig_covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t xchr_model, double ci_size, double vif_thresh, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end);

void LogisticTest();

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_GLM_H__
