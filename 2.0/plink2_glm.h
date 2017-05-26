#ifndef __PLINK2_GLM_H__
#define __PLINK2_GLM_H__

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


#include "plink2_adjust.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfGlm0,
  kfGlmZs = (1 << 0),

  // mutually exclusive
  kfGlmSex = (1 << 1),
  kfGlmNoXSex = (1 << 2),

  // mutually exclusive
  kfGlmGenotypic = (1 << 3),
  kfGlmHethom = (1 << 4),
  kfGlmDominant = (1 << 5),
  kfGlmRecessive = (1 << 6),

  kfGlmInteraction = (1 << 7),
  kfGlmHideCovar = (1 << 8),
  kfGlmIntercept = (1 << 9),
  kfGlmFirthFallback = (1 << 10),
  kfGlmFirth = (1 << 11),
  kfGlmPerm = (1 << 12),
  kfGlmPermCount = (1 << 13),
  kfGlmConditionDominant = (1 << 14),
  kfGlmConditionRecessive = (1 << 15),
  kfGlmLocalOmitLast = (1 << 16),
  kfGlmTestsAll = (1 << 17)
FLAGSET_DEF_END(glm_flags_t);

FLAGSET_DEF_START()
  kfGlmCol0,
  kfGlmColChrom = (1 << 0),
  kfGlmColPos = (1 << 1),
  kfGlmColRef = (1 << 2),
  kfGlmColAlt1 = (1 << 3),
  kfGlmColAlt = (1 << 4),
  kfGlmColAltcount = (1 << 5),
  kfGlmColTotallele = (1 << 6),
  kfGlmColAltcountcc = (1 << 7),
  kfGlmColTotallelecc = (1 << 8),
  kfGlmColAltfreq = (1 << 9),
  kfGlmColAltfreqcc = (1 << 10),
  kfGlmColMachR2 = (1 << 11),
  kfGlmColFirthYn = (1 << 12),
  kfGlmColTest = (1 << 13),
  kfGlmColNobs = (1 << 14),

  // if beta specified, ignore orbeta
  kfGlmColBeta = (1 << 15),
  kfGlmColOrbeta = (1 << 16),
  
  kfGlmColSe = (1 << 17),
  kfGlmColCi = (1 << 18),
  kfGlmColT = (1 << 19),
  kfGlmColP = (1 << 20),
  kfGlmColDefault = (kfGlmColChrom | kfGlmColPos | kfGlmColRef | kfGlmColAlt | kfGlmColFirthYn | kfGlmColTest | kfGlmColNobs | kfGlmColOrbeta | kfGlmColSe | kfGlmColCi | kfGlmColT | kfGlmColP),
  kfGlmColAll = ((kfGlmColCi * 2) - kfGlmColChrom)
FLAGSET_DEF_END(glm_cols_t);

typedef struct glm_info_struct {
  glm_flags_t flags;
  glm_cols_t cols;
  uint32_t mperm_ct;
  uint32_t local_cat_ct;
  double max_corr;
  char* condition_varname;
  char* condition_list_fname;
  range_list_t parameters_range_list;
  range_list_t tests_range_list;
} glm_info_t;

void init_glm(glm_info_t* glm_info_ptr);

void cleanup_glm(glm_info_t* glm_info_ptr);

// for testing purposes
// plink2_matrix.h must be included in this file
// boolerr_t logistic_regression(const float* yy, const float* xx, uint32_t sample_ct, uint32_t predictor_ct, float* coef, float* ll, float* pp, float* vv, float* hh, float* grad, float* dcoef);

// boolerr_t firth_regression(const float* yy, const float* xx, uint32_t sample_ct, uint32_t predictor_ct, float* coef, float* hh, matrix_finvert_buf1_t* inv_1d_buf, float* flt_2d_buf, float* pp, float* vv, float* grad, float* dcoef, float* ww, float* tmpnxk_buf);

pglerr_t glm_main(const uintptr_t* orig_sample_include, const char* sample_ids, const char* sids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const pheno_col_t* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const glm_info_t* glm_info_ptr, const adjust_info_t* adjust_info_ptr, const aperm_t* aperm_ptr, const char* local_covar_fname, const char* local_pvar_fname, const char* local_psam_fname, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t orig_covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t xchr_model, double ci_size, double vif_thresh, double pfilter, double output_min_p, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end);

#ifdef __cplusplus
} // namespace plink2
#endif
 
#endif // __PLINK2_GLM_H__
