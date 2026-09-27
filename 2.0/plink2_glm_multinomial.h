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

#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "plink2_cmdline.h"
#include "plink2_common.h"
#include "plink2_glm_shared.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Multinomial logistic regression of a categorical phenotype with K levels on
// the covariates and m genotype columns G, one omnibus test per variant (or
// per allele) on (K-1) m degrees of freedom:
//   log(P(y = k) / P(y = ref)) = Z beta_k + G gamma_k,  k != ref
//   H0: gamma_k = 0 for every k.
// G is the additive dosage of each non-omitted allele (all tested jointly),
// or one allele's dominant/recessive/hetonly column or genotypic/hethom pair,
// with the other non-omitted alleles as additive nuisance columns in Z.
//
// "Level" indices are global to the phenotype, with 0 always the reference
// level.  A given sample set may not contain every level (e.g. chrY), and a
// variant's missing calls can remove a level from its regression; "class"
// indices are the compacted level indices of a single regression, with class
// 0 the lowest-indexed level present, which serves as that regression's
// reference.

typedef struct MultinomialSetStruct {
  // Per-sample level index.
  uint32_t* levels;
  // Intercept row, then the covariates (each centered and scaled to unit
  // variance, which leaves the likelihood unchanged); stride is sample_ctav.
  double* covars_pmaj;
  // Covariate-only model fitted to the whole sample set: class-major
  // (null_class_ct - 1) x (covar_ct + 1) coefficients.
  double* null_coefs;
  // level -> class index in the covariate-only fit, UINT32_MAX if absent.
  uint32_t* null_level_to_class;
  double null_ln_lik;
  uint32_t sample_ct;
  uint32_t covar_ct;  // excluding the intercept
  uint32_t null_class_ct;
  // Set when the covariate-only model has no maximum-likelihood estimate (a
  // covariate separates the levels), so that null_coefs holds the
  // Firth-penalized fit instead and every variant is fitted with Firth
  // regression.  null_ln_lik is then unused.
  uint32_t firth_null;
} MultinomialSet;

// One per output row.  In the additive model, a variant has one row, which
// tests all its non-omitted alleles jointly; in the other models, it has one
// row per non-omitted allele (A1), with the other non-omitted alleles as
// additive nuisance predictors.
typedef struct {
  uint32_t sample_obs_ct;
  uint32_t allele_obs_ct;
  double mach_r2;
  double min_expected;
  double chisq;
  // (non-reference classes) x (tested columns)
  uint32_t df;
  uint32_t is_unfinished;
  // Firth regression was used
  uint32_t is_firth;
  // nonzero if chisq is invalid
  uint64_t glm_err;
} MultinomialAuxResult;

typedef struct GlmMultinomialCtxStruct {
  GlmCtx* common;

  // main, chrX, chrY; the latter two are only used when common->sample_ct_x
  // (resp. sample_ct_y) is nonzero
  MultinomialSet sets[3];
  uint32_t level_ct;
  GlmMultinomialTest test_type;
  // 0: 'no-firth' (and always with the score test), 1: 'firth-fallback', 2:
  // 'firth'
  uint32_t firth_mode;
  // fit the full model and save coefficients + standard errors?
  uint32_t save_coefs;
  uint32_t is_additive;
  // genotype columns per allele: 1, or 2 for 'genotypic' (ADD, DOMDEV) and
  // 'hethom' (HOM, HET)
  uint32_t model_col_ct;
  // per-variant slot sizes: rows, A1 alleles per row, coefficient slots per
  // row and non-reference level
  uint32_t max_row_ct;
  uint32_t max_a1_ct;
  uint32_t max_tested_ct;

  // max_row_ct entries per variant
  MultinomialAuxResult* block_aux;
  // max_row_ct * max_a1_ct entries per variant
  double* block_a1_dosage;
  // max_row_ct * level_ct * max_a1_ct entries per variant
  double* block_level_a1;
  // level_ct entries per variant
  uint32_t* block_level_allele_obs;
  // max_row_ct * (level_ct - 1) * max_tested_ct (beta, se) pairs per variant,
  // se -9.0 if unavailable
  double* block_beta_se;
} GlmMultinomialCtx;

// Determines the levels of a categorical phenotype among the given samples.
// cat_to_level[] is indexed by category index (0 = missing), and is
// UINT32_MAX for categories absent from the samples; level_cat_idxs[] is the
// inverse map.  The reference level is ref_name if non-null (and present),
// the first category in natural-sort order otherwise.
// *ref_found_ptr is zero iff ref_name was specified but not found.
BoolErr GlmMultinomialInitLevels(const uintptr_t* sample_include, const PhenoCol* pheno_col, uint32_t sample_ct, const char* ref_name, uint32_t* level_ct_ptr, uint32_t** cat_to_level_ptr, uint32_t** level_cat_idxs_ptr, uint32_t* ref_found_ptr);

// Number of distinct non-missing categories among the given samples.
uint32_t CountPhenoCats(const uintptr_t* sample_include, const PhenoCol* pheno_col, uint32_t sample_ct);

// Fills *setp for the given samples, checks the covariates, and fits the
// covariate-only model.  If that fit fails, it is refitted with Firth
// regression when firth_mode is nonzero (setting setp->firth_null), and
// *glm_err_ptr is set to kGlmErrcodeLogisticConvergeFail (firth_mode zero) or
// kGlmErrcodeFirthConvergeFail if no fit succeeds.
BoolErr GlmAllocFillAndTestPhenoCovarsMultinomial(const uintptr_t* sample_include, const PhenoCol* pheno_col, const uint32_t* cat_to_level, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uint32_t level_ct, uintptr_t covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uint32_t firth_mode, MultinomialSet* setp, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr);

PglErr GlmMultinomial(const char* cur_pheno_name, const char* const* level_names, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const char* outname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, PgenFileInfo* pgfip, GlmMultinomialCtx* ctx, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, uintptr_t* valid_allele_ct_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_GLM_MULTINOMIAL_H__
