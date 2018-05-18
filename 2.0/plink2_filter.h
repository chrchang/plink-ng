#ifndef __PLINK2_FILTER_H__
#define __PLINK2_FILTER_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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


#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr FromToFlag(const char* const* variant_ids, const uint32_t* variant_id_htable, const char* varid_from, const char* varid_to, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uintptr_t* variant_include, ChrInfo* cip, uint32_t* variant_ct_ptr);

PglErr SnpFlag(const uint32_t* variant_bps, const char* const* variant_ids, const uint32_t* variant_id_htable, const char* varid_snp, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, int32_t window_bp, uintptr_t* variant_include, ChrInfo* cip, uint32_t* variant_ct_ptr);

PglErr SnpsFlag(const char* const* variant_ids, const uint32_t* variant_id_htable, const RangeList* snps_range_list_ptr, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

PglErr ExtractExcludeFlagNorange(const char* const* variant_ids, const uint32_t* variant_id_htable, const char* fname, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

void RandomThinProb(const char* flagname_p, const char* unitname, double thin_keep_prob, uint32_t raw_item_ct, uintptr_t* item_include, uint32_t* item_ct_ptr);

PglErr RandomThinCt(const char* flagname_p, const char* unitname, uint32_t thin_keep_ct, uint32_t raw_item_ct, uintptr_t* item_include, uint32_t* item_ct_ptr);

FLAGSET_DEF_START()
  kfKeep0,
  kfKeepRemove = (1 << 0),
  kfKeepFam = (1 << 1)
FLAGSET_DEF_END(KeepFlags);

PglErr KeepOrRemove(const char* fnames, const SampleIdInfo* siip, uint32_t raw_sample_ct, KeepFlags flags, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

PglErr KeepFcol(const char* fname, const SampleIdInfo* siip, const char* strs_flattened, const char* col_name, uint32_t raw_sample_ct, uint32_t col_num, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

PglErr RequirePheno(const PhenoCol* pheno_cols, const char* pheno_names, const char* require_pheno_flattened, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

PglErr KeepRemoveIf(const CmpExpr* cmp_expr, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t affection_01, uint32_t is_remove, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

PglErr KeepRemoveCats(const char* cats_fname, const char* cat_names_flattened, const char* cat_phenoname, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t is_remove, uint32_t max_thread_ct, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

void ComputeAlleleFreqs(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const uint64_t* founder_allele_dosages, uint32_t variant_ct, uint32_t maf_succ, double* allele_freqs);

PglErr ReadAlleleFreqs(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const char* read_freq_fname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t maf_succ, uint32_t max_thread_ct, double* allele_freqs);

void ComputeMajAlleles(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const double* allele_freqs, uint32_t variant_ct, AlleleCode* maj_alleles);

PglErr LoadSampleMissingCts(const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uint32_t* sample_missing_hc_cts, uint32_t* sample_missing_dosage_cts, uint32_t* sample_hethap_cts);

PglErr MindFilter(const uint32_t* sample_missing_cts, const uint32_t* sample_hethap_cts, const SampleIdInfo* siip, uint32_t raw_sample_ct, uint32_t variant_ct, uint32_t variant_ct_y, double mind_thresh, uintptr_t* sample_include, uintptr_t* sex_male, uint32_t* sample_ct_ptr, char* outname, char* outname_end);

void EnforceGenoThresh(const ChrInfo* cip, const uint32_t* variant_missing_cts, const uint32_t* variant_hethap_cts, uint32_t sample_ct, uint32_t male_ct, uint32_t first_hap_uidx, double geno_thresh, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

void EnforceHweThresh(const ChrInfo* cip, const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), const double* hwe_x_pvals, MiscFlags misc_flags, double hwe_thresh, uint32_t nonfounders, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

void EnforceMinorFreqConstraints(const uintptr_t* allele_idx_offsets, const uint64_t* founder_allele_dosages, const double* allele_freqs, double min_maf, double max_maf, uint64_t min_allele_dosage, uint64_t max_allele_dosage, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

void EnforceMachR2Thresh(const ChrInfo* cip, const double* mach_r2_vals, double mach_r2_min, double mach_r2_max, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

void EnforceMinBpSpace(const ChrInfo* cip, const uint32_t* variant_bps, uint32_t min_bp_space, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

PglErr SetRefalt1FromFile(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const TwoColParams* allele_flag_info, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t is_alt1, uint32_t force, uint32_t max_thread_ct, const char** allele_storage, uint32_t* max_allele_slen_ptr, STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uintptr_t* nonref_flags, uintptr_t* previously_seen);

PglErr RefFromFa(const uintptr_t* variant_include, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const ChrInfo* cip, const char* fname, uint32_t max_allele_slen, uint32_t force, STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uintptr_t* nonref_flags);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_FILTER_H__
