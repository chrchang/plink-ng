#ifndef __PLINK2_GLM_LINEAR_H__
#define __PLINK2_GLM_LINEAR_H__

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

#include "include/SFMT.h"
#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "include/plink2_text.h"
#include "plink2_cmdline.h"
#include "plink2_common.h"
#include "plink2_glm_shared.h"

#ifdef __cplusplus
namespace plink2 {
#endif

typedef struct {
  // double beta;
  // double se;
  //   zval = beta / se
  //   width of asymptotic CI = ci_zt * se
  //   T-statistic = zval
  //   pval = TstatToP(zval, sample_obs_ct - predictor_ct)

  uint32_t sample_obs_ct;

  uint32_t allele_obs_ct;
  double a1_dosage;

  double mach_r2;
} LinearAuxResult;

typedef struct GlmLinearCtxStruct {
  GlmCtx* common;

  const double* pheno_d;
  const double* pheno_x_d;
  const double* pheno_y_d;
  const double* covars_cmaj_d;
  const double* covars_cmaj_x_d;
  const double* covars_cmaj_y_d;
  double* local_covars_vcmaj_d[2];
  LinearAuxResult* block_aux;

  uint32_t max_difflist_len;
  uint32_t subbatch_size;
} GlmLinearCtx;

BoolErr GlmAllocFillAndTestCovarsQt(const uintptr_t* sample_include, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, RegressionNmPrecomp** nm_precomp_ptr, double** covars_cmaj_d_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr);

BoolErr GlmAllocFillAndTestPhenoCovarsQt(const uintptr_t* sample_include, const double* pheno_qt, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uint32_t is_qt_residualize, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, double** pheno_d_ptr, RegressionNmPrecomp** nm_precomp_ptr, double** covars_cmaj_d_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr);

PglErr GlmLinear(const char* cur_pheno_name, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const char* outname, uint32_t raw_variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLinearCtx* ctx, TextStream* local_covar_txsp, LlStr** outfnames_ll_ptr, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, uintptr_t* valid_allele_ct_ptr);

PglErr GlmLinearBatch(const uintptr_t* pheno_batch, const PhenoCol* pheno_cols, const char* pheno_names, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, uint32_t raw_variant_ct, uint32_t completed_pheno_ct, uint32_t batch_size, uintptr_t max_pheno_name_blen, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLinearCtx* ctx, TextStream* local_covar_txsp, LlStr** outfnames_ll_ptr, char* outname, char* outname_end);

PglErr GlmLinearPerm(const char* cur_pheno_name, const PhenoCol* orig_pheno_col, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* valid_alleles, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const APerm* aperm_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const double* orig_ln_pvals, const PhenoCol* permute_within_phenocol, uint32_t raw_sample_ct, uint32_t raw_variant_ct, uintptr_t valid_allele_ct, uint32_t max_chr_blen, double ln_pfilter, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, MiscFlags misc_flags, PgenFileInfo* pgfip, GlmLinearCtx* ctx, sfmt_t* sfmtp, TextStream* local_covar_txsp, uintptr_t* remaining_variants, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_GLM_LINEAR_H__
