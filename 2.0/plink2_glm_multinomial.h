#ifndef __PLINK2_GLM_MULTINOMIAL_H__
#define __PLINK2_GLM_MULTINOMIAL_H__

// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang, Benjamin Demaille.
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

// Multinomial logistic regression on a categorical phenotype (--glm with
// --mnl-ref): per variant, the model is fitted with and without the genotype
// term, and a likelihood-ratio test with (category count - 1) degrees of
// freedom is reported along with per-category Wald statistics from the full
// fit.

#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "plink2_cmdline.h"
#include "plink2_common.h"
#include "plink2_glm_shared.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Phenotype and covariate data for one sample set (the main one, or the
// chrX/chrY one when those have different samples or covariates).
typedef struct GlmMnlSetStruct {
  // Category of each sample, 0 = the --mnl-ref reference category, 1..(cat_ct
  // - 1) = the other categories present, in their natural-sorted order.
  uint32_t* pheno_cats;
  // Covariate-major, sample_ct stride (not vector-aligned).
  double* covars_cmaj;
  // Covariate-only fit on the whole sample set, used unchanged for every
  // variant without missing calls and as the starting point otherwise.
  // (cat_ct - 1) rows of (1 + covar_ct) coefficients, intercept first.
  double* null_coefs;
  // Per-sample log-likelihood contributions of that fit.
  double* null_lli;
  // Names of the nonreference categories, in output order.
  const char** cat_names;
  uint32_t cat_ct;
} GlmMnlSet;

typedef struct {
  uint32_t sample_obs_ct;
  uint32_t allele_obs_ct;
  double a1_dosage;
  double mach_r2;
} MnlAuxResult;

typedef struct GlmMultinomialCtxStruct {
  GlmCtx* common;

  GlmMnlSet mnl_set;
  GlmMnlSet mnl_set_x;
  GlmMnlSet mnl_set_y;
  MnlAuxResult* block_aux;
} GlmMultinomialCtx;

// Looks up the --mnl-ref category name for the given phenotype; returns
// nullptr if the phenotype has no entry.
const char* MnlRefCatname(const char* mnl_ref_flattened, const char* pheno_name);

// Sets *cat_ct_ptr to the number of categories present among the given
// samples, or to 0 if ref_cat_idx is not among them.  Returns 1 on
// out-of-memory.
BoolErr MnlCountCats(const uintptr_t* sample_include, const PhenoCol* pheno_col, uint32_t sample_ct, uint32_t ref_cat_idx, uint32_t* cat_ct_ptr);

// Fills *mnl_set_ptr, checks the covariates the same way as the other --glm
// regressions do, and fits the covariate-only model.  On failure of either,
// *glm_err_ptr is set and the return value is still 0.
BoolErr GlmAllocFillAndTestPhenoCovarsMnl(const uintptr_t* sample_include, const PhenoCol* pheno_col, uint32_t ref_cat_idx, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, GlmMnlSet* mnl_set_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr);

PglErr GlmMultinomial(const char* cur_pheno_name, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const char* outname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, PgenFileInfo* pgfip, GlmMultinomialCtx* ctx, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, uintptr_t* valid_allele_ct_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_GLM_MULTINOMIAL_H__
