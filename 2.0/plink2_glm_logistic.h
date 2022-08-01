#ifndef __PLINK2_GLM_LOGISTIC_H__
#define __PLINK2_GLM_LOGISTIC_H__

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


#include "plink2_glm_shared.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// could split this into per-variant and per-tested-allele parts, but only 12
// bytes (sample_obs_ct, multiallelic mach_r2) can go into the former; probably
// unimportant
typedef struct {
  // double beta;
  //   odds ratio = exp(beta)
  // double se;
  //   zval = beta / se
  //   width of asymptotic CI (beta units) = ci_zt * se
  //   T-statistic = zval
  //   pval = ZscoreToP(zval)

  uint32_t sample_obs_ct;

  uint32_t allele_obs_ct;
  double a1_dosage;

  uint16_t firth_fallback;
  uint16_t is_unfinished;
  uint32_t case_allele_obs_ct;
  double a1_case_dosage;

  double mach_r2;

  // case hom-ref, case ref-alt, case alt-alt, ctrl hom-ref, ...
  STD_ARRAY_DECL(uint32_t, 6, geno_hardcall_cts);
} LogisticAuxResult;

typedef struct CcResidualizeCtxStruct {
  float* logistic_nm_sample_offsets;
  float* firth_nm_sample_offsets;
  uint32_t prefitted_pred_ct;
  uint32_t domdev_present_p1;
  uint32_t sample_ct;
} CcResidualizeCtx;

typedef struct GlmLogisticCtxStruct {
  GlmCtx *common;

  uintptr_t* pheno_cc;
  uintptr_t* pheno_x_cc;
  uintptr_t* pheno_y_cc;
  uintptr_t* gcount_case_interleaved_vec;
  uintptr_t* gcount_case_interleaved_vec_x;
  uintptr_t* gcount_case_interleaved_vec_y;
  const float* pheno_f;
  float* pheno_x_f;
  float* pheno_y_f;
  const float* covars_cmaj_f;
  float* covars_cmaj_x_f;
  float* covars_cmaj_y_f;
  CcResidualizeCtx* cc_residualize;
  CcResidualizeCtx* cc_residualize_x;
  CcResidualizeCtx* cc_residualize_y;
  uint16_t separation_found;
  uint16_t separation_found_x;
  uint16_t separation_found_y;
  float* local_covars_vcmaj_f[2];
  LogisticAuxResult* block_aux;
} GlmLogisticCtx;

BoolErr GlmAllocFillAndTestPhenoCovarsCc(const uintptr_t* sample_include, const uintptr_t* pheno_cc, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uint32_t domdev_present_p1, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, GlmFlags glm_flags, uintptr_t** pheno_cc_collapsed_ptr, uintptr_t** gcount_case_interleaved_vec_ptr, float** pheno_f_ptr, RegressionNmPrecomp** nm_precomp_ptr, float** covars_cmaj_f_ptr, CcResidualizeCtx** cc_residualize_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr);

PglErr GlmLogistic(const char* cur_pheno_name, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const char* outname, uint32_t raw_variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLogisticCtx* ctx, TextStream* local_covar_txsp, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, double* orig_permstat, uintptr_t* valid_allele_ct_ptr);

// void LogisticTestInternal();

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_GLM_LOGISTIC_H__
