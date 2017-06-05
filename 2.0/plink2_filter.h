#ifndef __PLINK2_FILTER_H__
#define __PLINK2_FILTER_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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

// order is such that (kCmpOperatorEq - x) is the inverse of x
ENUM_U31_DEF_START()
  kCmpOperatorNoteq,
  kCmpOperatorLe,
  kCmpOperatorLeq,

  kCmpOperatorGe,
  kCmpOperatorGeq,
  kCmpOperatorEq
ENUM_U31_DEF_END(cmp_binary_op_t);

typedef struct {
  // Restrict to [pheno/covar name] [operator] [pheno val] for now.  could
  // support or/and, parentheses, etc. later.

  // Currently stores null-terminated pheno/covar name, followed by
  // null-terminated value string.  Storage format needs to be synced with
  // plink2.cpp validate_and_alloc_cmp_expr().
  char* pheno_name;
  cmp_binary_op_t binary_op;
} cmp_expr_t;

void init_cmp_expr(cmp_expr_t* cmp_expr_ptr);

void cleanup_cmp_expr(cmp_expr_t* cmp_expr_ptr);

pglerr_t from_to_flag(char** variant_ids, const uint32_t* variant_id_htable, const char* varid_from, const char* varid_to, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uintptr_t* variant_include, chr_info_t* cip, uint32_t* variant_ct_ptr);

pglerr_t snp_flag(const uint32_t* variant_bps, char** variant_ids, const uint32_t* variant_id_htable, const char* varid_snp, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, int32_t window_bp, uintptr_t* variant_include, chr_info_t* cip, uint32_t* variant_ct_ptr);

pglerr_t snps_flag(char** variant_ids, const uint32_t* variant_id_htable, const range_list_t* snps_range_list_ptr, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

pglerr_t extract_exclude_flag_norange(char** variant_ids, const uint32_t* variant_id_htable, const char* fname, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

FLAGSET_DEF_START()
  kfKeep0,
  kfKeepRemove = (1 << 0),
  kfKeepFam = (1 << 1),
  kfKeepForceSid = (1 << 2)
FLAGSET_DEF_END(keep_flags_t);

pglerr_t keep_or_remove(const char* fname, const char* sample_ids, const char* sids, uint32_t raw_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, keep_flags_t flags, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

pglerr_t require_pheno(const pheno_col_t* pheno_cols, const char* pheno_names, char* require_pheno_flattened, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

pglerr_t keep_remove_if(const cmp_expr_t* cmp_expr, const pheno_col_t* pheno_cols, const char* pheno_names, const pheno_col_t* covar_cols, const char* covar_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t affection_01, uint32_t is_remove, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

pglerr_t keep_remove_cats(const char* cats_fname, const char* cat_names_flattened, const char* cat_phenoname, const pheno_col_t* pheno_cols, const char* pheno_names, const pheno_col_t* covar_cols, const char* covar_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t is_remove, uint32_t max_thread_ct, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

void compute_allele_freqs(const uintptr_t* variant_include, const uintptr_t* variant_allele_idxs, const uint64_t* founder_allele_dosages, uint32_t variant_ct, uint32_t maf_succ, double* allele_freqs);

pglerr_t read_allele_freqs(const uintptr_t* variant_include, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const char* read_freq_fname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_alt_allele_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t maf_succ, uint32_t max_thread_ct, double* allele_freqs);

void compute_maj_alleles(const uintptr_t* variant_include, const uintptr_t* variant_allele_idxs, const double* allele_freqs, uint32_t variant_ct, alt_allele_ct_t* maj_alleles);

pglerr_t load_sample_missing_cts(const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, uint32_t* sample_missing_hc_cts, uint32_t* sample_missing_dosage_cts, uint32_t* sample_hethap_cts);

pglerr_t mind_filter(const uint32_t* sample_missing_cts, const uint32_t* sample_hethap_cts, const char* sample_ids, const char* sids, uint32_t raw_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t variant_ct, uint32_t variant_ct_y, double mind_thresh, uintptr_t* sample_include, uintptr_t* sex_male, uint32_t* sample_ct_ptr, char* outname, char* outname_end);

void enforce_geno_thresh(const chr_info_t* cip, const uint32_t* variant_missing_cts, const uint32_t* variant_hethap_cts, uint32_t sample_ct, uint32_t male_ct, uint32_t first_hap_uidx, double geno_thresh, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

void enforce_hwe_thresh(const chr_info_t* cip, const uint32_t* founder_raw_geno_cts, const uint32_t* founder_x_male_geno_cts, const uint32_t* founder_x_nosex_geno_cts, const double* hwe_x_pvals, misc_flags_t misc_flags, double hwe_thresh, uint32_t nonfounders, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

void enforce_minor_freq_constraints(const uintptr_t* variant_allele_idxs, const uint64_t* founder_allele_dosages, const double* allele_freqs, double min_maf, double max_maf, uint64_t min_allele_dosage, uint64_t max_allele_dosage, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

void enforce_mach_r2_thresh(const chr_info_t* cip, const double* mach_r2_vals, double mach_r2_min, double mach_r2_max, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_FILTER_H__
