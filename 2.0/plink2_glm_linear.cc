// This file is part of PLINK 2.00, copyright (C) 2005-2023 Shaun Purcell,
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

#include "include/plink2_stats.h"
#include "plink2_compress_stream.h"
#include "plink2_glm_linear.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// xtx_state 0: either interactions or local covariates present, no xtx_image
// xtx_state 1: only additive effect
// xtx_state 2: additive and domdev effects
// Note that the immediate return value only corresponds to out-of-memory.
// Other errors are indicated by glm_err_ptr on a *false* return
// value, due to how the caller is expected to handle them.
BoolErr GlmAllocFillAndTestCovarsQt(const uintptr_t* sample_include, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, RegressionNmPrecomp** nm_precomp_ptr, double** covars_cmaj_d_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  if (unlikely(bigstack_alloc_kcp(new_covar_ct, cur_covar_names_ptr) ||
               bigstack_alloc_d(new_nonlocal_covar_ct * sample_ct, covars_cmaj_d_ptr))) {
    return 1;
  }
  double* corr_buf = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  *nm_precomp_ptr = nullptr;
  if (xtx_state) {
    assert(!local_covar_ct);
    if (unlikely(BIGSTACK_ALLOC_X(RegressionNmPrecomp, 1, nm_precomp_ptr))) {
      return 1;
    }
    // x^2 + (2*xtx_state + 2)x + (5*xtx_state - 1)
    // 2x^2 + (2*xtx_state + 4)x + (5*xtx_state)
    // 3x^2 + (2*xtx_state + 4)x + (5*xtx_state)
    // 4x^2 + (4*xtx_state + 4)x + (8*xtx_state - 2)
    // 4x^2 + (4*xtx_state + 5)x + (8*xtx_state - 2)
    //
    // xt_y_image: (1 + xtx_state + x) elements per phenotype in subbatch;
    //   now allocated later
    if (unlikely(bigstack_alloc_d(8 * xtx_state - 2 +
                                  new_covar_ct * (5 + 4 * xtx_state + 4 * new_covar_ct),
                                  &((*nm_precomp_ptr)->xtx_image)) ||
                 bigstack_alloc_d(new_covar_ct * (new_covar_ct + 1), &corr_buf))) {
      return 1;
    }
    (*nm_precomp_ptr)->covarx_dotprod_inv = &((*nm_precomp_ptr)->xtx_image[(new_covar_ct + xtx_state + 1) * (new_covar_ct + xtx_state + 1)]);
    (*nm_precomp_ptr)->corr_inv = &((*nm_precomp_ptr)->covarx_dotprod_inv[(new_covar_ct + 1) * (new_covar_ct + 1)]);
    (*nm_precomp_ptr)->corr_image = &((*nm_precomp_ptr)->corr_inv[new_covar_ct * new_covar_ct]);
    (*nm_precomp_ptr)->corr_inv_sqrts = &((*nm_precomp_ptr)->corr_image[(new_covar_ct + xtx_state) * (new_covar_ct + xtx_state)]);
    (*nm_precomp_ptr)->xt_y_image = nullptr;  // defensive
    bigstack_mark = R_CAST(unsigned char*, corr_buf);
  }
  double* covar_dotprod;
  double* inverse_corr_buf;
  if (unlikely(bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &covar_dotprod) ||
               bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &inverse_corr_buf))) {
    return 1;
  }
  PglErr reterr = GlmFillAndTestCovars(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, max_corr, vif_thresh, covar_dotprod, corr_buf, inverse_corr_buf, *covars_cmaj_d_ptr, *cur_covar_names_ptr, glm_err_ptr);
  if (unlikely(reterr)) {
    return (reterr == kPglRetNomem);
  }
  if (xtx_state) {
    if (InitNmPrecomp(*covars_cmaj_d_ptr, covar_dotprod, corr_buf, inverse_corr_buf, sample_ct, 1, new_covar_ct, xtx_state, *nm_precomp_ptr)) {
      *glm_err_ptr = SetGlmErr0(kGlmErrcodeVifInfinite);
      return 0;  // not out-of-memory error
    }
  }
  BigstackReset(bigstack_mark);
  return 0;
}

void FillPhenoAndXtY(const uintptr_t* sample_include, const double* __restrict pheno_qt, const double* __restrict covars_cmaj_d, uintptr_t sample_ct, uintptr_t domdev_present_p1, uintptr_t covar_ct, double* xt_y_image, double* __restrict pheno_d) {
  double* pheno_d_iter = pheno_d;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  double pheno_sum = 0.0;
  for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    const double cur_pheno = pheno_qt[sample_uidx];
    *pheno_d_iter++ = cur_pheno;
    pheno_sum += cur_pheno;
  }
  if (!xt_y_image) {
    return;
  }
  xt_y_image[0] = pheno_sum;
  ZeroDArr(domdev_present_p1, &(xt_y_image[1]));
  ColMajorVectorMatrixMultiplyStrided(pheno_d, covars_cmaj_d, sample_ct, sample_ct, covar_ct, &(xt_y_image[1 + domdev_present_p1]));
}

BoolErr GlmAllocFillAndTestPhenoCovarsQt(const uintptr_t* sample_include, const double* pheno_qt, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, double** pheno_d_ptr, RegressionNmPrecomp** nm_precomp_ptr, double** covars_cmaj_d_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  if (unlikely(GlmAllocFillAndTestCovarsQt(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, max_corr, vif_thresh, xtx_state, nm_precomp_ptr, covars_cmaj_d_ptr, cur_covar_names_ptr, glm_err_ptr))) {
    return 1;
  }
  if (*glm_err_ptr) {
    // this is a bit messy
    return 0;
  }
  if (unlikely(bigstack_alloc_d(sample_ct, pheno_d_ptr))) {
    return 1;
  }
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  double* xt_y_image = nullptr;
  if (xtx_state) {
    if (unlikely(bigstack_alloc_d(1 + xtx_state + new_covar_ct, &xt_y_image))) {
      return 1;
    }
    (*nm_precomp_ptr)->xt_y_image = xt_y_image;
  }
  FillPhenoAndXtY(sample_include, pheno_qt, *covars_cmaj_d_ptr, sample_ct, xtx_state, new_covar_ct, xt_y_image, *pheno_d_ptr);
  return 0;
}

uintptr_t GetLinearWorkspaceSize(uint32_t sample_ct, uint32_t biallelic_predictor_ct, uint32_t max_extra_allele_ct, uint32_t constraint_ct, uint32_t xmain_ct) {
  // sample_ct * max_predictor_ct < 2^31, and max_predictor_ct < sqrt(2^31), so
  // no overflows

  // sample_nm, tmp_nm = sample_ctl words
  uintptr_t workspace_size = 2 * RoundUpPow2(BitCtToWordCt(sample_ct) * sizeof(intptr_t), kCacheline);

  // nm_pheno_buf = sample_ct doubles
  workspace_size += RoundUpPow2(sample_ct * sizeof(double), kCacheline);

  const uint32_t max_predictor_ct = biallelic_predictor_ct + max_extra_allele_ct;
  // predictors_pmaj = (max_predictor_ct + main_mutated + main_omitted) * sample_ct
  // doubles
  workspace_size += RoundUpPow2((max_predictor_ct + xmain_ct) * sample_ct * sizeof(double), kCacheline);

  // xtx_inv = max_predictor_ct * max_predictor_ct doubles
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ct * sizeof(double), kCacheline);

  // dbl_2d_buf = max_predictor_ct * max(max_predictor_ct, 7) doubles
  workspace_size += RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);

  // inverse_corr_buf = (max_predictor_ct - 1) * max(max_predictor_ct - 1, 4) doubles
  workspace_size += RoundUpPow2((max_predictor_ct - 1) * MAXV((max_predictor_ct - 1), 4) * sizeof(double), kCacheline);

  // semicomputed_biallelic_corr_matrix = (max_predictor_ct - 1)^2 doubles
  workspace_size += RoundUpPow2((biallelic_predictor_ct - 1) * (biallelic_predictor_ct - 1) * sizeof(double), kCacheline);

  // semicomputed_biallelic_inv_corr_sqrts = biallelic_predictor_ct doubles
  workspace_size += RoundUpPow2(biallelic_predictor_ct * sizeof(double), kCacheline);

  // fitted_coefs, xt_y = max_predictor_ct doubles
  workspace_size += 2 * RoundUpPow2(max_predictor_ct * sizeof(double), kCacheline);

  // inv_1d_buf
  workspace_size += RoundUpPow2(MAXV(max_predictor_ct, constraint_ct) * kMatrixInvertBuf1CheckedAlloc, kCacheline);

  // machr2_dosage_sums, machr2_dosage_ssqs
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(uint64_t) * 2, kCacheline);

  if (constraint_ct) {
    // tmphxs_buf, h_transpose_buf = constraint_ct * max_predictor_ct doubles
    workspace_size += 2 * RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(double), kCacheline);

    // inner_buf = constraint_ct * constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * constraint_ct * sizeof(double), kCacheline);

    // cur_constraints_con_major = constraint_ct * max_predictor_ct doubles
    workspace_size += RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(double), kCacheline);
  }
  return workspace_size;
}

// possible todo: delete this, and GlmLinear(), if GlmLinearBatchThread is good
// enough at the same job.
THREAD_FUNC_DECL GlmLinearThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GlmLinearCtx* ctx = S_CAST(GlmLinearCtx*, arg->sharedp->context);
  GlmCtx* common = ctx->common;

  PgenReader* pgrp = common->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = common->genovecs[tidx];
  pgv.dosage_present = nullptr;
  pgv.dosage_main = nullptr;
  if (common->dosage_presents) {
    pgv.dosage_present = common->dosage_presents[tidx];
    pgv.dosage_main = common->dosage_mains[tidx];
  }
  unsigned char* workspace_buf = common->workspace_bufs[tidx];
  const uintptr_t* variant_include = common->variant_include;
  const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
  const AlleleCode* omitted_alleles = common->omitted_alleles;
  const uintptr_t* sex_male_collapsed = common->sex_male_collapsed;
  const ChrInfo* cip = common->cip;
  const uint32_t* subset_chr_fo_vidx_start = common->subset_chr_fo_vidx_start;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const GlmFlags glm_flags = common->glm_flags;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t model_hetonly = (glm_flags / kfGlmHetonly) & 1;
  const uint32_t joint_genotypic = (glm_flags / kfGlmGenotypic) & 1;
  const uint32_t joint_hethom = (glm_flags / kfGlmHethom) & 1;
  const uint32_t domdev_present = joint_genotypic || joint_hethom;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t reported_pred_uidx_start = 1 - include_intercept;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uintptr_t max_reported_test_ct = common->max_reported_test_ct;
  const uintptr_t local_covar_ct = common->local_covar_ct;
  const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
  const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!model_dominant) && (!model_recessive) && (!model_hetonly) && (!common->tests_flag) && (!add_interactions);
  uintptr_t max_sample_ct = MAXV(common->sample_ct, common->sample_ct_x);
  if (max_sample_ct < common->sample_ct_y) {
    max_sample_ct = common->sample_ct_y;
  }
  pgv.patch_01_set = nullptr;
  pgv.patch_01_vals = nullptr;
  pgv.patch_10_set = nullptr;
  pgv.patch_10_vals = nullptr;
  if (common->thread_mhc) {
    const uint32_t max_sample_ctl = BitCtToWordCt(max_sample_ct);
    pgv.patch_01_set = common->thread_mhc[tidx];
    pgv.patch_01_vals = R_CAST(AlleleCode*, &(pgv.patch_01_set[max_sample_ctl]));
    AlleleCode* patch_01_vals_end = &(pgv.patch_01_vals[max_sample_ct]);
    AlignACToVec(&patch_01_vals_end);
    pgv.patch_10_set = R_CAST(uintptr_t*, patch_01_vals_end);
    pgv.patch_10_vals = R_CAST(AlleleCode*, &(pgv.patch_10_set[max_sample_ctl]));
  }
  pgv.patch_01_ct = 0;
  pgv.patch_10_ct = 0;
  pgv.multidosage_sample_ct = 0;
  uint32_t variant_idx_offset = 0;
  uint32_t allele_ct = 2;
  uint32_t omitted_allele_idx = 0;
  uint32_t extra_regression_ct = 0;
  double main_dosage_sum = 0.0;
  double main_dosage_ssq = 0.0;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);

    uintptr_t allele_bidx = variant_bidx;
    if (max_extra_allele_ct) {
      allele_bidx = variant_bidx + CountExtraAlleles(variant_include, allele_idx_offsets, common->read_variant_uidx_starts[0], common->read_variant_uidx_starts[tidx], 0);
    }
    double* beta_se_iter = common->block_beta_se;
    if (beta_se_multiallelic_fused) {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * variant_bidx]);
    } else {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * allele_bidx]);
    }

    LinearAuxResult* block_aux_iter = &(ctx->block_aux[allele_bidx]);
    const double* local_covars_iter = nullptr;
    if (local_covar_ct) {
      // &(nullptr[0]) is okay in C++, but undefined in C
      local_covars_iter = &(ctx->local_covars_vcmaj_d[parity][variant_bidx * max_sample_ct * local_covar_ct]);
    }
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = CountSortedSmallerU32(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      uint32_t cur_variant_bidx_end = subset_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
      if (cur_variant_bidx_end > variant_bidx_end) {
        cur_variant_bidx_end = variant_bidx_end;
      }
      const uint32_t is_haploid = IsSet(cip->haploid_mask, chr_idx);
      const uint32_t is_regular_x = is_haploid && (chr_idx == x_code);
      const uint32_t is_y = (chr_idx == y_code);
      const uint32_t is_nonx_haploid = is_haploid && (!is_regular_x);
      const uintptr_t* cur_sample_include;
      const uint32_t* cur_sample_include_cumulative_popcounts;
      const double* cur_pheno;
      const RegressionNmPrecomp* nm_precomp;
      const double* cur_covars_cmaj;
      const uintptr_t* cur_parameter_subset;
      const uintptr_t* cur_joint_test_params;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      uint32_t cur_constraint_ct;
      if (is_y && common->sample_include_y) {
        cur_sample_include = common->sample_include_y;
        cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
        cur_pheno = ctx->pheno_y_d;
        nm_precomp = common->nm_precomp_y;
        cur_covars_cmaj = ctx->covars_cmaj_y_d;
        cur_parameter_subset = common->parameter_subset_y;
        cur_joint_test_params = common->joint_test_params_y;
        cur_sample_ct = common->sample_ct_y;
        cur_covar_ct = common->covar_ct_y;
        cur_constraint_ct = common->constraint_ct_y;
      } else if (is_regular_x && common->sample_include_x) {
        cur_sample_include = common->sample_include_x;
        cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
        cur_pheno = ctx->pheno_x_d;
        nm_precomp = common->nm_precomp_x;
        cur_covars_cmaj = ctx->covars_cmaj_x_d;
        cur_parameter_subset = common->parameter_subset_x;
        cur_joint_test_params = common->joint_test_params_x;
        cur_sample_ct = common->sample_ct_x;
        cur_covar_ct = common->covar_ct_x;
        cur_constraint_ct = common->constraint_ct_x;
      } else {
        cur_sample_include = common->sample_include;
        cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
        cur_pheno = ctx->pheno_d;
        nm_precomp = common->nm_precomp;
        cur_covars_cmaj = ctx->covars_cmaj_d;
        cur_parameter_subset = common->parameter_subset;
        cur_joint_test_params = common->joint_test_params;
        cur_sample_ct = common->sample_ct;
        cur_covar_ct = common->covar_ct;
        cur_constraint_ct = common->constraint_ct;
      }
      const uint32_t sample_ctl = BitCtToWordCt(cur_sample_ct);
      const uint32_t sample_ctl2 = NypCtToWordCt(cur_sample_ct);
      const uint32_t cur_biallelic_predictor_ct_base = 2 + domdev_present + cur_covar_ct * (1 + add_interactions * domdev_present_p1);
      uint32_t cur_biallelic_predictor_ct = cur_biallelic_predictor_ct_base;
      uint32_t literal_covar_ct = cur_covar_ct;
      if (cur_parameter_subset) {
        cur_biallelic_predictor_ct = PopcountWords(cur_parameter_subset, BitCtToWordCt(cur_biallelic_predictor_ct_base));
        literal_covar_ct = PopcountBitRange(cur_parameter_subset, 2 + domdev_present, 2 + domdev_present + cur_covar_ct);
      }
      const uint32_t max_predictor_ct = cur_biallelic_predictor_ct + max_extra_allele_ct;
      uint32_t reported_pred_uidx_biallelic_end;
      if (hide_covar) {
        if (!cur_parameter_subset) {
          reported_pred_uidx_biallelic_end = 2 + domdev_present;
        } else {
          reported_pred_uidx_biallelic_end = IsSet(cur_parameter_subset, 1) + domdev_present_p1;
        }
      } else {
        reported_pred_uidx_biallelic_end = cur_biallelic_predictor_ct;
      }
      // nm_predictors_pmaj_buf may require up to two extra columns omitted
      // from the main regression.
      // 1. In the multiallelic dominant/recessive/hetonly/hethom cases, the
      //    original genotype column does not appear in the regression, and
      //    we'd rather not reconstruct it from genovec, etc. when we need to
      //    swap it out for another allele, so we keep the original genotype in
      //    an extra column.
      //    To reduce code bloat, we now handle the biallelic cases in the same
      //    way; this is one of the more peripheral code paths so adding more
      //    complexity to speed it up is less justifiable.
      // 2. If --parameters excludes the main (possibly
      //    dominant/recessive/hetonly) genotype column but does care about an
      //    interaction, we want a copy of what the main genotype column's
      //    contents would have been to refer to.
      const uint32_t main_omitted = cur_parameter_subset && (!IsSet(cur_parameter_subset, 1));
      const uint32_t main_mutated = model_dominant || model_recessive || model_hetonly || joint_hethom;
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* tmp_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      double* nm_pheno_buf = S_CAST(double*, arena_alloc_raw_rd(cur_sample_ct * sizeof(double), &workspace_iter));
      double* nm_predictors_pmaj_buf = S_CAST(double*, arena_alloc_raw_rd((max_predictor_ct + main_mutated + main_omitted) * cur_sample_ct * sizeof(double), &workspace_iter));
      double* xtx_inv = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ct * sizeof(double), &workspace_iter));
      double* fitted_coefs = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * sizeof(double), &workspace_iter));
      double* xt_y = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_corr_matrix = S_CAST(double*, arena_alloc_raw_rd((cur_biallelic_predictor_ct - 1) * (cur_biallelic_predictor_ct - 1) * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_inv_corr_sqrts = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));
      MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, arena_alloc_raw_rd(MAXV(max_predictor_ct, cur_constraint_ct) * kMatrixInvertBuf1CheckedAlloc, &workspace_iter));
      double* dbl_2d_buf = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), &workspace_iter));
      uint64_t* machr2_dosage_sums = S_CAST(uint64_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(uint64_t) * 2, &workspace_iter));
      uint64_t* machr2_dosage_ssqs = &(machr2_dosage_sums[max_extra_allele_ct + 2]);

      // could technically have this overlap fitted_coefs/xt_y, but that sets
      // the stage for future bugs
      double* inverse_corr_buf = S_CAST(double*, arena_alloc_raw_rd((max_predictor_ct - 1) * MAXV((max_predictor_ct - 1), 4) * sizeof(double), &workspace_iter));

      // joint test only
      double* tmphxs_buf = nullptr;
      double* h_transpose_buf = nullptr;
      double* inner_buf = nullptr;
      double* cur_constraints_con_major = nullptr;
      if (cur_constraint_ct) {
        tmphxs_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        h_transpose_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        inner_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * cur_constraint_ct * sizeof(double), &workspace_iter));
        cur_constraints_con_major = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        ZeroDArr(cur_constraint_ct * max_predictor_ct, cur_constraints_con_major);
        const uint32_t first_joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
        cur_constraints_con_major[first_joint_test_idx] = 1.0;
        // Rest of this matrix must be updated later, since cur_predictor_ct
        // changes at multiallelic variants.
      }
      assert(S_CAST(uintptr_t, workspace_iter - workspace_buf) == GetLinearWorkspaceSize(cur_sample_ct, cur_biallelic_predictor_ct, max_extra_allele_ct, cur_constraint_ct, main_mutated + main_omitted));
      const double pheno_ssq_base = DotprodD(cur_pheno, cur_pheno, cur_sample_ct);
      const double cur_sample_ct_recip = 1.0 / u31tod(cur_sample_ct);
      const double cur_sample_ct_m1_recip = 1.0 / u31tod(cur_sample_ct - 1);
      const uint32_t sparse_optimization_eligible = (!is_regular_x) && nm_precomp;
      double geno_d_lookup[2];
      if (sparse_optimization_eligible) {
        if (model_hetonly) {
          geno_d_lookup[0] = 1.0;
          geno_d_lookup[1] = 0.0;
        } else {
          geno_d_lookup[1] = 1.0;
          if (is_nonx_haploid) {
            geno_d_lookup[0] = 0.5;
          } else if (model_recessive || joint_hethom) {
            geno_d_lookup[0] = 0.0;
          } else {
            geno_d_lookup[0] = 1.0;
            if (!model_dominant) {
              geno_d_lookup[1] = 2.0;
            }
          }
        }
      }
      const double* xtx_image = nullptr;
      const double* covarx_dotprod_inv = nullptr;
      const double* corr_inv = nullptr;
      const double* xt_y_image = nullptr;
      if (nm_precomp) {
        xtx_image = nm_precomp->xtx_image;
        covarx_dotprod_inv = nm_precomp->covarx_dotprod_inv;
        corr_inv = nm_precomp->corr_inv;
        const uintptr_t nongeno_pred_ct = cur_biallelic_predictor_ct - domdev_present - 2;
        const uintptr_t nonintercept_biallelic_pred_ct = cur_biallelic_predictor_ct - 1;
        memcpy(semicomputed_biallelic_corr_matrix, nm_precomp->corr_image, nonintercept_biallelic_pred_ct * nonintercept_biallelic_pred_ct * sizeof(double));
        memcpy(&(semicomputed_biallelic_inv_corr_sqrts[domdev_present_p1]), nm_precomp->corr_inv_sqrts, nongeno_pred_ct * sizeof(double));
        xt_y_image = nm_precomp->xt_y_image;
      }
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      // when this is set, the last fully-processed variant had no missing
      // genotypes, and if the current variant also has no missing genotypes we
      // may be able to skip reinitialization of most of
      // nm_predictors_pmaj_buf.
      uint32_t prev_nm = 0;

      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        if (allele_idx_offsets) {
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
          if (!beta_se_multiallelic_fused) {
            extra_regression_ct = allele_ct - 2;
          }
        }
        const uint32_t allele_ct_m2 = allele_ct - 2;
        const uint32_t expected_predictor_ct = cur_biallelic_predictor_ct + allele_ct_m2;
        PglErr reterr;
        if (!allele_ct_m2) {
          reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
        } else {
          reterr = PgrGetMD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, &pgv);
          // todo: proper multiallelic dosage support
        }
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto GlmLinearThread_err;
        }
        ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
        GenoarrCountFreqsUnsafe(pgv.genovec, cur_sample_ct, genocounts);
        uint32_t missing_ct = genocounts[3];
        if (!missing_ct) {
          SetAllBits(cur_sample_ct, sample_nm);
        } else {
          GenoarrToNonmissing(pgv.genovec, cur_sample_ct, sample_nm);
          if (pgv.dosage_ct) {
            BitvecOr(pgv.dosage_present, sample_ctl, sample_nm);
            missing_ct = cur_sample_ct - PopcountWords(sample_nm, sample_ctl);
          }
        }
        if (omitted_alleles) {
          omitted_allele_idx = omitted_alleles[variant_uidx];
        }
        uintptr_t const_alleles[DivUp(kPglMaxAlleleCt, kBitsPerWord)];
        const uint32_t allele_ctl = DivUp(allele_ct, kBitsPerWord);
        ZeroWArr(allele_ctl, const_alleles);
        const uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
        const uint32_t nm_sample_ctl = BitCtToWordCt(nm_sample_ct);
        // first predictor column: intercept
        if (!prev_nm) {
          for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
            nm_predictors_pmaj_buf[sample_idx] = 1.0;
          }
        }
        // second predictor column: main genotype
        double* genotype_vals = &(nm_predictors_pmaj_buf[nm_sample_ct]);
        if (main_mutated || main_omitted) {
          genotype_vals = &(nm_predictors_pmaj_buf[expected_predictor_ct * nm_sample_ct]);
        }
        double cur_pheno_ssq = pheno_ssq_base;
        uintptr_t sample_midx_base = 0;
        uintptr_t sample_nm_inv_bits = ~sample_nm[0];
        for (uint32_t missing_idx = 0; missing_idx != missing_ct; ++missing_idx) {
          const uintptr_t sample_midx = BitIter0(sample_nm, &sample_midx_base, &sample_nm_inv_bits);
          cur_pheno_ssq -= cur_pheno[sample_midx] * cur_pheno[sample_midx];
        }
        uint32_t sparse_optimization = 0;
        double* multi_start = nullptr;
        if (!allele_ct_m2) {
          // When prev_nm is set and missing_ct is zero, we don't need to call
          // MultiplySelfTranspose() on nm_predictors_pmaj_buf to get all the
          // predictor x predictor dot products; instead we patch in the
          // genotype x predictor dot products that may change, copying the
          // rest from xtx_image.
          // As a side effect, it is no longer strictly necessary to fill the
          // genotype row of nm_predictors_pmaj_buf.  sparse_optimization
          // indicates that plink 1.9's QT --assoc sparse dot product algorithm
          // will be used instead.
          // probable todos: allow a few dosages to be present, cover chrX
          // case.
          if (omitted_allele_idx) {
            GenovecInvertUnsafe(cur_sample_ct, pgv.genovec);
            ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
            }
            const uint32_t uii = genocounts[0];
            genocounts[0] = genocounts[2];
            genocounts[2] = uii;
          }
          uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
          uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
          if (!missing_ct) {
            // originally had genocounts[0] > 0.875 * nm_sample_ct threshold,
            // but then tried this on high-MAF data and it was still
            // substantially faster
            sparse_optimization = sparse_optimization_eligible && (!pgv.dosage_ct) && prev_nm;
            if (!sparse_optimization) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, genotype_vals);
              if (pgv.dosage_ct) {
                uintptr_t sample_idx_base = 0;
                uintptr_t dosage_present_bits = pgv.dosage_present[0];
                for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                  const uintptr_t sample_idx = BitIter1(pgv.dosage_present, &sample_idx_base, &dosage_present_bits);
                  const uintptr_t dosage_val = pgv.dosage_main[dosage_idx];
                  // 32768 -> 2, 16384 -> 1, 0 -> 0
                  genotype_vals[sample_idx] = kRecipDosageMid * swtod(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_idx);
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                }
              }
            }
          } else {
            if (!pgv.dosage_ct) {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, genotype_vals);
            } else {
              sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              uint32_t dosage_idx = 0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_midx);
                double cur_val;
                if (IsSet(pgv.dosage_present, sample_midx)) {
                  const uintptr_t dosage_val = pgv.dosage_main[dosage_idx++];
                  cur_val = kRecipDosageMid * swtod(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                } else {
                  cur_val = kSmallDoubles[cur_geno];
                }
                genotype_vals[sample_idx] = cur_val;
              }
            }
          }
          // Check for constant genotype column.
          // (Technically, we should recheck later in the chrX no-sex-covariate
          // --xchr-model 1 corner case.)
          if (!pgv.dosage_ct) {
            if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          } else if (pgv.dosage_ct == nm_sample_ct) {
            if (DosageIsConstant(dosage_sum, dosage_ssq, nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          }
          machr2_dosage_sums[1 - omitted_allele_idx] = dosage_sum;
          machr2_dosage_ssqs[1 - omitted_allele_idx] = dosage_ssq;
          machr2_dosage_sums[omitted_allele_idx] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[omitted_allele_idx] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
        } else {
          // multiallelic.
          // If some but not all alleles have constant dosages, we remove just
          // those alleles from the regressions; trim-alts is not necessary to
          // see what's going on with the other alleles.  To reduce parsing
          // complexity, the number of output lines is not affected by this;
          // the ones corresponding to the constant alleles have NA values.
          // Punt on sparse_optimization for now; may be worth revisiting after
          // multiallelic dosage implemented.
          // dosage_ct == 0 temporarily guaranteed if we reach here.
          assert(!pgv.dosage_ct);
          multi_start = &(nm_predictors_pmaj_buf[(expected_predictor_ct - allele_ct_m2) * nm_sample_ct]);
          ZeroU64Arr(allele_ct, machr2_dosage_sums);
          ZeroU64Arr(allele_ct, machr2_dosage_ssqs);
          // postpone multiply for now, since no multiallelic dosages
          machr2_dosage_sums[0] = genocounts[1] + 2 * genocounts[0];
          machr2_dosage_ssqs[0] = genocounts[1] + 4LLU * genocounts[0];
          if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
            SetBit(0, const_alleles);
          }
          if (omitted_allele_idx) {
            // Main genotype column starts as REF.
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallInvDoublePairs, nm_sample_ct, genotype_vals);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallInvDoubles, cur_sample_ct, genotype_vals);
            }
          }
          uint32_t rare_allele_ct = allele_ct_m2;
          double* alt1_start = nullptr;
          double* rarealt_start = multi_start;
          if (omitted_allele_idx != 1) {
            if (omitted_allele_idx) {
              alt1_start = multi_start;
              rarealt_start = &(rarealt_start[nm_sample_ct]);
              --rare_allele_ct;
            } else {
              alt1_start = genotype_vals;
            }
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, alt1_start);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, alt1_start);
            }
          }
          ZeroDArr(rare_allele_ct * nm_sample_ct, rarealt_start);
          // Use sums as ones[] and ssqs as twos[] for rarealts; transform to
          // actual sums/ssqs later.
          if (pgv.patch_01_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_01_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_01_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ct + sample_idx] = 1.0;
                alt1_start[sample_idx] = 0.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ct + sample_idx] = 1.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                alt1_start[sample_idx] = 0.0;
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                machr2_dosage_sums[allele_code] += 1;
                if (allele_code == omitted_allele_idx) {
                  continue;
                }
                const uint32_t cur_col = allele_code - 2 - (allele_code > omitted_allele_idx);
                rarealt_start[cur_col * nm_sample_ct + sample_idx] = 1.0;
              }
            }
          }
          uintptr_t alt1_het_ct = genocounts[1] - pgv.patch_01_ct;
          if (pgv.patch_10_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_10_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_10_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] = 2.0;
                  alt1_start[sample_idx] = 0.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ct + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] += 1.0;
                    alt1_start[sample_idx] = 0.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] = 2.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ct + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] += 1.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  machr2_dosage_ssqs[ac0] += 1;
                  alt1_start[sample_idx] = 0.0;
                  if (ac0 != omitted_allele_idx) {
                    const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                    rarealt_start[ac0_col * nm_sample_ct + sample_idx] = 2.0;
                  }
                } else {
                  machr2_dosage_sums[ac1] += 1;
                  if (ac1 != omitted_allele_idx) {
                    const uint32_t ac1_col = ac1 - 2 - (ac1 > omitted_allele_idx);
                    rarealt_start[ac1_col * nm_sample_ct + sample_idx] = 1.0;
                  }
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    machr2_dosage_sums[ac0] += 1;
                    alt1_start[sample_idx] = 0.0;
                    if (ac0 != omitted_allele_idx) {
                      const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                      rarealt_start[ac0_col * nm_sample_ct + sample_idx] += 1.0;
                    }
                  }
                }
              }
            }
          }
          for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
            const uintptr_t one_ct = machr2_dosage_sums[allele_idx];
            const uintptr_t two_ct = machr2_dosage_ssqs[allele_idx];
            machr2_dosage_sums[allele_idx] = one_ct + 2 * two_ct;
            machr2_dosage_ssqs[allele_idx] = one_ct + 4LLU * two_ct;
            if ((one_ct == nm_sample_ct) || (two_ct == nm_sample_ct) || ((!one_ct) && (!two_ct))) {
              SetBit(allele_idx, const_alleles);
            }
          }
          const uintptr_t alt1_hom_ct = genocounts[2] - pgv.patch_10_ct;
          machr2_dosage_sums[1] = alt1_het_ct + 2 * alt1_hom_ct;
          machr2_dosage_ssqs[1] = alt1_het_ct + 4LLU * alt1_hom_ct;
          if ((alt1_het_ct == nm_sample_ct) || (alt1_hom_ct == nm_sample_ct) || ((!alt1_het_ct) && (!alt1_hom_ct))) {
            SetBit(1, const_alleles);
          }
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            machr2_dosage_sums[allele_idx] *= 0x4000LLU;
            machr2_dosage_ssqs[allele_idx] *= 0x10000000LLU;
          }
        }
        // usually need to save some of {sample_obs_ct, allele_obs_ct,
        // a1_dosage, mach_r2 even for skipped variants
        // compute them all for now, could conditionally skip later
        uint32_t allele_obs_ct = nm_sample_ct * 2;
        if (!is_regular_x) {
          if (is_nonx_haploid) {
            allele_obs_ct = nm_sample_ct;
            // everything is on 0..1 scale, not 0..2
            if (!sparse_optimization) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                genotype_vals[sample_idx] *= 0.5;
              }
              const uint32_t high_ct = nm_sample_ct * allele_ct_m2;
              for (uint32_t uii = 0; uii != high_ct; ++uii) {
                multi_start[uii] *= 0.5;
              }
            }
          }
        } else {
          CopyBitarrSubset(sex_male_collapsed, sample_nm, nm_sample_ct, tmp_nm);
          const uintptr_t* male_nm = tmp_nm;
          const uint32_t nm_male_ct = PopcountWords(male_nm, nm_sample_ctl);
          if (is_xchr_model_1) {
            // special case: multiply male values by 0.5
            uintptr_t sample_idx_base = 0;
            uintptr_t male_nm_bits = male_nm[0];
            for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
              const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
              genotype_vals[sample_idx] *= 0.5;
              // could insert multiallelic loop here instead, but I'm guessing
              // that's worse due to locality of writes?
            }
            for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
              double* cur_start = &(multi_start[extra_allele_idx * nm_sample_ct]);
              sample_idx_base = 0;
              male_nm_bits = male_nm[0];
              for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
                const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
                cur_start[sample_idx] *= 0.5;
              }
            }
            allele_obs_ct -= nm_male_ct;
          }
        }
        const double mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, allele_ct);
        for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
          if (allele_idx == omitted_allele_idx) {
            continue;
          }
          block_aux_iter->sample_obs_ct = nm_sample_ct;
          block_aux_iter->allele_obs_ct = allele_obs_ct;
          double a1_dosage = u63tod(machr2_dosage_sums[allele_idx]) * kRecipDosageMid;
          if (is_xchr_model_1) {
            // ugh.
            double* geno_col = genotype_vals;
            if (allele_idx > (!omitted_allele_idx)) {
              geno_col = &(nm_predictors_pmaj_buf[(expected_predictor_ct - (allele_ct - allele_idx) + (allele_idx < omitted_allele_idx)) * nm_sample_ct]);
            }
            a1_dosage = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              a1_dosage += geno_col[sample_idx];
            }
            if (!allele_ct_m2) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = 0.0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const double cur_dosage = geno_col[sample_idx];
                main_dosage_ssq += cur_dosage * cur_dosage;
              }
            }
          } else {
            if (is_nonx_haploid) {
              a1_dosage *= 0.5;
            }
            if (!allele_ct_m2) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = u63tod(machr2_dosage_ssqs[allele_idx]) * kRecipDosageMidSq;
              if (is_nonx_haploid) {
                main_dosage_ssq *= 0.25;
              }
            }
          }
          block_aux_iter->a1_dosage = a1_dosage;
          block_aux_iter->mach_r2 = mach_r2;
          ++block_aux_iter;
        }
        // now free to skip the actual regression if there are too few samples,
        // or there's a zero-variance genotype column
        GlmErr glm_err = 0;
        if (nm_sample_ct <= expected_predictor_ct) {
          glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
        } else if (IsSet(const_alleles, omitted_allele_idx)) {
          glm_err = SetGlmErr0(kGlmErrcodeConstOmittedAllele);
        }
        if (glm_err) {
          if (missing_ct) {
            // covariates have not been copied yet, so we can't usually change
            // prev_nm from 0 to 1 when missing_ct == 0 (and there's little
            // reason to optimize the zero-covariate case).
            prev_nm = 0;
          }
          uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
          if (beta_se_multiallelic_fused || (!hide_covar)) {
            reported_ct += allele_ct_m2;
          }
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            for (uint32_t uii = 0; uii != reported_ct; ++uii) {
              memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
              beta_se_iter[uii * 2 + 1] = -9.0;
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        } else {
          uint32_t parameter_uidx = 2 + domdev_present;
          double* nm_predictors_pmaj_istart = nullptr;
          if (!sparse_optimization) {
            // only need to do this part once per variant in multiallelic case
            double* nm_predictors_pmaj_iter = &(nm_predictors_pmaj_buf[nm_sample_ct * (parameter_uidx - main_omitted)]);
            if (missing_ct || (!prev_nm)) {
              // fill phenotype
              sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                nm_pheno_buf[sample_idx] = cur_pheno[sample_midx];
              }

              // fill covariates
              for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx, ++parameter_uidx) {
                // strictly speaking, we don't need cur_covars_cmaj to be
                // vector-aligned
                if (cur_parameter_subset && (!IsSet(cur_parameter_subset, parameter_uidx))) {
                  continue;
                }
                const double* cur_covar_col;
                if (covar_idx < local_covar_ct) {
                  cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                } else {
                  cur_covar_col = &(cur_covars_cmaj[(covar_idx - local_covar_ct) * cur_sample_ct]);
                }
                sample_midx_base = 0;
                sample_nm_bits = sample_nm[0];
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                  *nm_predictors_pmaj_iter++ = cur_covar_col[sample_midx];
                }
              }
              nm_predictors_pmaj_istart = nm_predictors_pmaj_iter;
              prev_nm = !(missing_ct || local_covar_ct);
            } else {
              // bugfix (15 Aug 2018): this was not handling --parameters
              // correctly when a covariate was only needed as part of an
              // interaction
              parameter_uidx += cur_covar_ct;
              nm_predictors_pmaj_istart = &(nm_predictors_pmaj_iter[literal_covar_ct * nm_sample_ct]);
            }
          }
          const uint32_t const_allele_ct = PopcountWords(const_alleles, allele_ctl);
          if (const_allele_ct) {
            // Must delete constant-allele columns from nm_predictors_pmaj, and
            // shift later columns back.
            double* read_iter = genotype_vals;
            double* write_iter = genotype_vals;
            for (uint32_t read_allele_idx = 0; read_allele_idx != allele_ct; ++read_allele_idx) {
              if (read_allele_idx == omitted_allele_idx) {
                continue;
              }
              if (!IsSet(const_alleles, read_allele_idx)) {
                if (write_iter != read_iter) {
                  memcpy(write_iter, read_iter, nm_sample_ct * sizeof(double));
                }
                if (write_iter == genotype_vals) {
                  write_iter = multi_start;
                } else {
                  write_iter = &(write_iter[nm_sample_ct]);
                }
              }
              if (read_iter == genotype_vals) {
                read_iter = multi_start;
              } else {
                read_iter = &(read_iter[nm_sample_ct]);
              }
            }
          }
          const uint32_t cur_predictor_ct = expected_predictor_ct - const_allele_ct;
          uint32_t nonconst_extra_regression_idx = UINT32_MAX;  // deliberate overflow
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            if (extra_regression_ct) {
              if (IsSet(const_alleles, extra_regression_idx + (extra_regression_idx >= omitted_allele_idx))) {
                glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
                goto GlmLinearThread_skip_regression;
              }
              ++nonconst_extra_regression_idx;
              if (nonconst_extra_regression_idx) {
                double* swap_target = &(multi_start[(nonconst_extra_regression_idx - 1) * nm_sample_ct]);
                for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
                  double dxx = genotype_vals[uii];
                  genotype_vals[uii] = swap_target[uii];
                  swap_target[uii] = dxx;
                }
              }
            }
            if (!sparse_optimization) {
              double* main_vals = &(nm_predictors_pmaj_buf[nm_sample_ct]);
              double* domdev_vals = nullptr;
              if (main_omitted) {
                // if main_mutated, this will be filled below
                // if not, this aliases genotype_vals
                main_vals = &(nm_predictors_pmaj_buf[(cur_predictor_ct + main_mutated) * nm_sample_ct]);
              } else if (joint_genotypic || joint_hethom) {
                domdev_vals = &(main_vals[nm_sample_ct]);
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 2.0 - cur_genotype_val;
                  }
                  domdev_vals[sample_idx] = cur_genotype_val;
                }
              }
              if (model_dominant) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..1..1
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 1.0;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              } else if (model_recessive || joint_hethom) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..0..1
                  if (cur_genotype_val < 1.0) {
                    cur_genotype_val = 0.0;
                  } else {
                    cur_genotype_val -= 1.0;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              } else if (model_hetonly) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..1..0
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 2.0 - cur_genotype_val;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              }

              // fill interaction terms
              if (add_interactions) {
                double* nm_predictors_pmaj_iter = nm_predictors_pmaj_istart;
                for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
                  const double* cur_covar_col;
                  if (covar_idx < local_covar_ct) {
                    cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                  } else {
                    cur_covar_col = &(cur_covars_cmaj[covar_idx * cur_sample_ct]);
                  }
                  if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                    sample_midx_base = 0;
                    uintptr_t sample_nm_bits = sample_nm[0];
                    for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                      const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                      *nm_predictors_pmaj_iter++ = main_vals[sample_idx] * cur_covar_col[sample_midx];
                    }
                  }
                  ++parameter_uidx;
                  if (domdev_present) {
                    if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                      sample_midx_base = 0;
                      uintptr_t sample_nm_bits = sample_nm[0];
                      for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                        const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                        *nm_predictors_pmaj_iter++ = domdev_vals[sample_idx] * cur_covar_col[sample_midx];
                      }
                    }
                    ++parameter_uidx;
                  }
                }
              }
            }

            // bugfix (12 Sep 2017): forgot to implement per-variant VIF and
            // max-corr checks
            if (xtx_image && prev_nm && (!allele_ct_m2)) {
              // only need to fill in additive and possibly domdev dot
              // products
              memcpy(xtx_inv, xtx_image, cur_predictor_ct * cur_predictor_ct * sizeof(double));
              memcpy(xt_y, xt_y_image, cur_predictor_ct * sizeof(double));
              if (sparse_optimization) {
                // currently does not handle chrX
                double geno_pheno_prod = 0.0;
                double domdev_pheno_prod = 0.0;
                double domdev_geno_prod = 0.0;
                double* geno_dotprod_row = &(xtx_inv[cur_predictor_ct]);
                double* domdev_dotprod_row = &(xtx_inv[2 * cur_predictor_ct]);
                for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
                  uintptr_t geno_word = pgv.genovec[widx];
                  if (geno_word) {
                    const uint32_t sample_idx_base = widx * kBitsPerWordD2;
                    do {
                      const uint32_t lowest_set_bit = ctzw(geno_word);
                      // since there are no missing values, we have a het if
                      // (lowest_set_bit & 1) is zero, and a hom-alt when
                      // it's one.
                      const uint32_t sample_idx = sample_idx_base + (lowest_set_bit / 2);
                      const double geno_d = geno_d_lookup[lowest_set_bit & 1];
                      const double cur_pheno_val = nm_pheno_buf[sample_idx];
                      geno_pheno_prod += geno_d * cur_pheno_val;
                      for (uintptr_t pred_idx = domdev_present + 2; pred_idx != cur_predictor_ct; ++pred_idx) {
                        geno_dotprod_row[pred_idx] += geno_d * nm_predictors_pmaj_buf[pred_idx * nm_sample_ct + sample_idx];
                      }
                      // can have a separate categorical loop here

                      if (domdev_present && (!(lowest_set_bit & 1))) {
                        // domdev = 1
                        domdev_pheno_prod += cur_pheno_val;
                        domdev_geno_prod += geno_d;
                        for (uintptr_t pred_idx = 3; pred_idx != cur_predictor_ct; ++pred_idx) {
                          domdev_dotprod_row[pred_idx] += nm_predictors_pmaj_buf[pred_idx * nm_sample_ct + sample_idx];
                        }
                        // categorical optimization possible here
                      }
                      geno_word &= geno_word - 1;
                    } while (geno_word);
                  }
                }
                xt_y[1] = geno_pheno_prod;
                const double het_ctd = u31tod(genocounts[1]);
                const double homalt_ctd = u31tod(genocounts[2]);
                xtx_inv[cur_predictor_ct] = het_ctd * geno_d_lookup[0] + homalt_ctd * geno_d_lookup[1];
                xtx_inv[cur_predictor_ct + 1] = het_ctd * geno_d_lookup[0] * geno_d_lookup[0] + homalt_ctd * geno_d_lookup[1] * geno_d_lookup[1];
                if (domdev_present) {
                  xt_y[2] = domdev_pheno_prod;
                  xtx_inv[cur_predictor_ct + 2] = domdev_geno_prod;
                  xtx_inv[2 * cur_predictor_ct] = het_ctd;
                  xtx_inv[2 * cur_predictor_ct + 2] = het_ctd;
                }
              } else {
                // !sparse_optimization
                xt_y[1] = DotprodD(&(nm_predictors_pmaj_buf[nm_sample_ct]), nm_pheno_buf, nm_sample_ct);
                uintptr_t start_pred_idx = 0;
                if (!(model_dominant || model_recessive || model_hetonly || joint_hethom)) {
                  start_pred_idx = domdev_present + 2;
                  xtx_inv[cur_predictor_ct] = main_dosage_sum;
                  xtx_inv[cur_predictor_ct + 1] = main_dosage_ssq;
                }
                if (cur_predictor_ct > start_pred_idx) {
                  // categorical optimization possible here
                  ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[nm_sample_ct]), &(nm_predictors_pmaj_buf[start_pred_idx * nm_sample_ct]), nm_sample_ct, nm_sample_ct, cur_predictor_ct - start_pred_idx, &(xtx_inv[cur_predictor_ct + start_pred_idx]));
                }
                if (domdev_present) {
                  xt_y[2] = DotprodD(&(nm_predictors_pmaj_buf[2 * nm_sample_ct]), nm_pheno_buf, nm_sample_ct);
                  // categorical optimization possible here
                  ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[2 * nm_sample_ct]), nm_predictors_pmaj_buf, nm_sample_ct, nm_sample_ct, cur_predictor_ct, &(xtx_inv[2 * cur_predictor_ct]));
                  xtx_inv[cur_predictor_ct + 2] = xtx_inv[2 * cur_predictor_ct + 1];
                }
              }
              glm_err = CheckMaxCorrAndVifNm(xtx_inv, corr_inv, cur_predictor_ct, domdev_present_p1, cur_sample_ct_recip, cur_sample_ct_m1_recip, max_corr, vif_thresh, semicomputed_biallelic_corr_matrix, semicomputed_biallelic_inv_corr_sqrts, dbl_2d_buf, &(dbl_2d_buf[2 * cur_predictor_ct]), &(dbl_2d_buf[3 * cur_predictor_ct]));
              if (glm_err) {
                goto GlmLinearThread_skip_regression;
              }
              const double geno_ssq = xtx_inv[1 + cur_predictor_ct];
              if (!domdev_present) {
                xtx_inv[1 + cur_predictor_ct] = xtx_inv[cur_predictor_ct];
                if (InvertRank1Symm(covarx_dotprod_inv, &(xtx_inv[1 + cur_predictor_ct]), cur_predictor_ct - 1, 1, geno_ssq, dbl_2d_buf, inverse_corr_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                  goto GlmLinearThread_skip_regression;
                }
              } else {
                const double domdev_geno_prod = xtx_inv[2 + cur_predictor_ct];
                const double domdev_ssq = xtx_inv[2 + 2 * cur_predictor_ct];
                xtx_inv[2 + cur_predictor_ct] = xtx_inv[cur_predictor_ct];
                xtx_inv[2 + 2 * cur_predictor_ct] = xtx_inv[2 * cur_predictor_ct];
                if (InvertRank2Symm(covarx_dotprod_inv, &(xtx_inv[2 + cur_predictor_ct]), cur_predictor_ct - 2, cur_predictor_ct, 1, geno_ssq, domdev_geno_prod, domdev_ssq, dbl_2d_buf, inverse_corr_buf, &(inverse_corr_buf[2 * (cur_predictor_ct - 2)]))) {
                  glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                  goto GlmLinearThread_skip_regression;
                }
              }
              // need to make sure xtx_inv remains reflected in NOLAPACK case
              memcpy(xtx_inv, dbl_2d_buf, cur_predictor_ct * cur_predictor_ct * sizeof(double));
              ReflectMatrix(cur_predictor_ct, xtx_inv);
              ColMajorVectorMatrixMultiplyStrided(xt_y, xtx_inv, cur_predictor_ct, cur_predictor_ct, cur_predictor_ct, fitted_coefs);
            } else {
              // generic case
              // major categorical optimization possible here, some
              // multiallelic optimizations possible
              MultiplySelfTranspose(nm_predictors_pmaj_buf, cur_predictor_ct, nm_sample_ct, xtx_inv);

              for (uint32_t pred_idx = 1; pred_idx != cur_predictor_ct; ++pred_idx) {
                dbl_2d_buf[pred_idx] = xtx_inv[pred_idx * cur_predictor_ct];
              }
              glm_err = CheckMaxCorrAndVif(xtx_inv, 1, cur_predictor_ct, nm_sample_ct, max_corr, vif_thresh, dbl_2d_buf, nullptr, inverse_corr_buf, inv_1d_buf);
              if (glm_err) {
                goto GlmLinearThread_skip_regression;
              }
              if (LinearRegressionInv(nm_pheno_buf, nm_predictors_pmaj_buf, cur_predictor_ct, nm_sample_ct, 1, xtx_inv, fitted_coefs, xt_y, inv_1d_buf, dbl_2d_buf)) {
                glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                goto GlmLinearThread_skip_regression;
              }
            }
            {
              // RSS = y^T y - y^T X (X^T X)^{-1} X^T y
              //     = cur_pheno_ssq - xt_y * fitted_coefs
              // s^2 = RSS / df
              // possible todo: improve numerical stability of this computation
              // in non-mean-centered phenotype case
              const double sigma = (cur_pheno_ssq - DotprodxD(xt_y, fitted_coefs, cur_predictor_ct)) / u31tod(nm_sample_ct - cur_predictor_ct);
              for (uint32_t uii = 0; uii != cur_predictor_ct; ++uii) {
                double* s_iter = &(xtx_inv[uii * cur_predictor_ct]);
#ifdef NOLAPACK
                for (uint32_t ujj = 0; ujj != cur_predictor_ct; ++ujj) {
                  s_iter[ujj] *= sigma;
                }
#else
                for (uint32_t ujj = 0; ujj <= uii; ++ujj) {
                  s_iter[ujj] *= sigma;
                }
#endif
              }
              // validParameters() check
              for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
                const double xtx_inv_diag_element = xtx_inv[pred_uidx * (cur_predictor_ct + 1)];
                if (xtx_inv_diag_element < 1e-20) {
                  glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                  goto GlmLinearThread_skip_regression;
                }
                // use dbl_2d_buf[] to store diagonal square roots
                dbl_2d_buf[pred_uidx] = sqrt(xtx_inv_diag_element);
              }
              dbl_2d_buf[0] = sqrt(xtx_inv[0]);
              for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
                const double cur_xtx_inv_diag_sqrt = 0.99999 * dbl_2d_buf[pred_uidx];
                const double* xtx_inv_row = &(xtx_inv[pred_uidx * cur_predictor_ct]);
                for (uint32_t pred_uidx2 = 0; pred_uidx2 != pred_uidx; ++pred_uidx2) {
                  if (xtx_inv_row[pred_uidx2] > cur_xtx_inv_diag_sqrt * dbl_2d_buf[pred_uidx2]) {
                    glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                    goto GlmLinearThread_skip_regression;
                  }
                }
              }
              double* beta_se_iter2 = beta_se_iter;
              for (uint32_t pred_uidx = reported_pred_uidx_start; pred_uidx != reported_pred_uidx_biallelic_end; ++pred_uidx) {
                // In the multiallelic-fused case, if the first allele is
                // constant, this writes the beta/se values for the first
                // nonconstant, non-omitted allele where the results for the
                // first allele belong.  We correct that below.
                *beta_se_iter2++ = fitted_coefs[pred_uidx];
                *beta_se_iter2++ = dbl_2d_buf[pred_uidx];
              }
              // move this up since dbl_2d_buf may be clobbered
              if (cur_constraint_ct) {
                const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeRankDeficient);
                memcpy(beta_se_iter2, &glm_err2, 8);
                ++beta_se_iter2;
                *beta_se_iter2++ = -9.0;
              }
              if (!const_allele_ct) {
                if (beta_se_multiallelic_fused || (!hide_covar)) {
                  for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
                    beta_se_iter2[2 * extra_allele_idx] = fitted_coefs[cur_biallelic_predictor_ct + extra_allele_idx];
                    beta_se_iter2[2 * extra_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + extra_allele_idx];
                  }
                }
              } else if (!beta_se_multiallelic_fused) {
                if (!hide_covar) {
                  // Need to insert some {CONST_ALLELE, -9} entries.
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                  const uint32_t cur_raw_allele_idx = extra_regression_idx + (extra_regression_idx >= omitted_allele_idx);
                  uint32_t extra_read_allele_idx = 0;
                  uint32_t extra_write_allele_idx = 0;
                  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    if ((allele_idx == omitted_allele_idx) || (allele_idx == cur_raw_allele_idx)) {
                      continue;
                    }
                    if (IsSet(const_alleles, allele_idx)) {
                      memcpy(&(beta_se_iter2[2 * extra_write_allele_idx]), &glm_err2, 8);
                      beta_se_iter2[2 * extra_write_allele_idx + 1] = -9.0;
                    } else {
                      beta_se_iter2[2 * extra_write_allele_idx] = fitted_coefs[cur_biallelic_predictor_ct + extra_read_allele_idx];
                      beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + extra_read_allele_idx];
                      ++extra_read_allele_idx;
                    }
                    ++extra_write_allele_idx;
                  }
                }
              } else {
                const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                // Special-case first nonconst allele since it's positioned
                // discontinuously, and its BETA/SE may already be correctly
                // filled.
                uint32_t allele_idx = omitted_allele_idx? 0 : 1;
                uint32_t extra_write_allele_idx = 0;
                if (IsSet(const_alleles, allele_idx)) {
                  memcpy(&(beta_se_iter[2 * include_intercept]), &glm_err2, 8);
                  beta_se_iter[2 * include_intercept + 1] = -9.0;
                  allele_idx = AdvTo0Bit(const_alleles, 1);
                  if (allele_idx == omitted_allele_idx) {
                    allele_idx = AdvTo0Bit(const_alleles, omitted_allele_idx + 1);
                  }
                  extra_write_allele_idx = allele_idx - 1 - (allele_idx > omitted_allele_idx);
                  for (uint32_t uii = 0; uii != extra_write_allele_idx; ++uii) {
                    memcpy(&(beta_se_iter2[2 * uii]), &glm_err2, 8);
                    beta_se_iter2[2 * uii + 1] = -9.0;
                  }
                  beta_se_iter2[2 * extra_write_allele_idx] = fitted_coefs[1];
                  beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[1];
                  ++extra_write_allele_idx;
                }
                ++allele_idx;
                uint32_t nonconst_allele_idx_m1 = 0;
                for (; allele_idx != allele_ct; ++allele_idx) {
                  if (allele_idx == omitted_allele_idx) {
                    continue;
                  }
                  if (!IsSet(const_alleles, allele_idx)) {
                    beta_se_iter2[2 * extra_write_allele_idx] = fitted_coefs[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                    beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                    ++nonconst_allele_idx_m1;
                  } else {
                    memcpy(&(beta_se_iter2[2 * extra_write_allele_idx]), &glm_err2, 8);
                    beta_se_iter2[2 * extra_write_allele_idx + 1] = -9.0;
                  }
                  ++extra_write_allele_idx;
                }
              }
              if (cur_constraint_ct) {
                uint32_t joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 1.0;
                }
#ifndef NOLAPACK
                // xtx_inv upper triangle was not filled
                ReflectMatrix(cur_predictor_ct, xtx_inv);
#endif
                double chisq;
                if (!LinearHypothesisChisq(fitted_coefs, cur_constraints_con_major, xtx_inv, cur_constraint_ct, cur_predictor_ct, cur_predictor_ct, &chisq, tmphxs_buf, h_transpose_buf, inner_buf, inv_1d_buf, dbl_2d_buf)) {
                  beta_se_iter2[-1] = chisq;
                }
                // next test may have different alt allele count
                joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 0.0;
                }
              }
            }
            while (0) {
            GlmLinearThread_skip_regression:
              {
                uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
                if (beta_se_multiallelic_fused || (!hide_covar)) {
                  reported_ct += allele_ct_m2;
                }
                for (uint32_t uii = 0; uii != reported_ct; ++uii) {
                  memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
                  beta_se_iter[uii * 2 + 1] = -9.0;
                }
              }
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        }
        // bugfix (1 Apr 2019): this needs to execute when regression is
        // skipped
        if (local_covars_iter) {
          local_covars_iter = &(local_covars_iter[local_covar_ct * max_sample_ct]);
        }
      }
    }
    parity = 1 - parity;
    variant_idx_offset += cur_block_variant_ct;
    while (0) {
    GlmLinearThread_err:
      UpdateU64IfSmaller(new_err_info, &common->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr GlmLinear(const char* cur_pheno_name, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const char* outname, uint32_t raw_variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLinearCtx* ctx, TextStream* local_covar_txsp, LlStr** outfnames_ll_ptr, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, uintptr_t* valid_allele_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  ThreadGroup tg;
  PreinitCstream(&css);
  PreinitThreads(&tg);
  {
    GlmCtx* common = ctx->common;
    const uintptr_t* variant_include = common->variant_include;
    const ChrInfo* cip = common->cip;
    const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
    const AlleleCode* omitted_alleles = common->omitted_alleles;

    const uint32_t sample_ct = common->sample_ct;
    const uint32_t sample_ct_x = common->sample_ct_x;
    const uint32_t sample_ct_y = common->sample_ct_y;
    const uint32_t covar_ct = common->covar_ct;
    const uintptr_t local_covar_ct = common->local_covar_ct;
    const uint32_t covar_ct_x = common->covar_ct_x;
    const uint32_t covar_ct_y = common->covar_ct_y;

    uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
    if (max_sample_ct < sample_ct_y) {
      max_sample_ct = sample_ct_y;
    }
    uint32_t* local_sample_idx_order = nullptr;
    uint32_t local_line_idx = 0;
    uint32_t local_xy = 0;  // 1 = chrX, 2 = chrY

    const char* local_line_iter = nullptr;
    uint32_t local_prev_chr_code = UINT32_MAX;
    uint32_t local_chr_code = UINT32_MAX;
    uint32_t local_bp = UINT32_MAX;
    uint32_t local_skip_chr = 1;
    if (local_covar_ct) {
      reterr = TextRewind(local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLinear_ret_TSTREAM_FAIL;
      }
      local_line_idx = glm_info_ptr->local_header_line_ct;
      reterr = TextSkip(local_line_idx, local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLinear_ret_TSTREAM_FAIL;
      }
      if (unlikely(bigstack_alloc_u32(local_sample_ct, &local_sample_idx_order))) {
        goto GlmLinear_ret_NOMEM;
      }
      for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
        const uint32_t cur_uidx = local_sample_uidx_order[uii];
        uint32_t cur_idx = UINT32_MAX;
        if ((cur_uidx != UINT32_MAX) && IsSet(common->sample_include, cur_uidx)) {
          cur_idx = RawToSubsettedPos(common->sample_include, common->sample_include_cumulative_popcounts, cur_uidx);
        }
        local_sample_idx_order[uii] = cur_idx;
      }
    }

    const uint32_t variant_ct = common->variant_ct;

    const GlmFlags glm_flags = glm_info_ptr->flags;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    // forced-singlethreaded
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto GlmLinear_ret_1;
    }
    if (outfnames_ll_ptr) {
      if (unlikely(PushLlStr(outname, outfnames_ll_ptr))) {
        goto GlmLinear_ret_NOMEM;
      }
    }
    const uint32_t report_neglog10p = (glm_flags / kfGlmLog10) & 1;
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t domdev_present_p1 = domdev_present + 1;

    const uint32_t constraint_ct = common->constraint_ct;
    const uint32_t constraint_ct_x = common->constraint_ct_x;
    const uint32_t constraint_ct_y = common->constraint_ct_y;

    const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
    uint32_t biallelic_predictor_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
    uint32_t biallelic_predictor_ct_x = 2 + covar_ct_x * (1 + add_interactions);
    uint32_t biallelic_predictor_ct_y = 2 + covar_ct_y * (1 + add_interactions);
    const uintptr_t* parameter_subset = common->parameter_subset;
    const uintptr_t* parameter_subset_x = common->parameter_subset_x;
    const uintptr_t* parameter_subset_y = common->parameter_subset_y;
    if (parameter_subset) {
      biallelic_predictor_ct = PopcountWords(parameter_subset, BitCtToWordCt(biallelic_predictor_ct));
      if (sample_ct_x) {
        biallelic_predictor_ct_x = PopcountWords(parameter_subset_x, BitCtToWordCt(biallelic_predictor_ct_x));
      } else {
        biallelic_predictor_ct_x = 0;
      }
      if (sample_ct_y) {
        // bugfix (7 Feb 2018): had biallelic_predictor_ct_x on right side
        // here, oops
        biallelic_predictor_ct_y = PopcountWords(parameter_subset_y, BitCtToWordCt(biallelic_predictor_ct_y));
      } else {
        biallelic_predictor_ct_y = 0;
      }
    }
    uint32_t biallelic_reported_test_ct = GetBiallelicReportedTestCt(parameter_subset, glm_flags, covar_ct, common->tests_flag);
    uintptr_t max_reported_test_ct = biallelic_reported_test_ct;
    uint32_t biallelic_reported_test_ct_x = 0;
    if (sample_ct_x) {
      biallelic_reported_test_ct_x = GetBiallelicReportedTestCt(parameter_subset_x, glm_flags, covar_ct_x, common->tests_flag);
      if (biallelic_reported_test_ct_x > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_x;
      }
    }
    uint32_t biallelic_reported_test_ct_y = 0;
    if (sample_ct_y) {
      biallelic_reported_test_ct_y = GetBiallelicReportedTestCt(parameter_subset_y, glm_flags, covar_ct_y, common->tests_flag);
      if (biallelic_reported_test_ct_y > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_y;
      }
    }
    const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
    const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
    const GlmColFlags glm_cols = glm_info_ptr->cols;
    const uint32_t test_col = glm_cols & kfGlmColTest;
    if (unlikely((!test_col) && (max_reported_test_ct > 1))) {
      // this is okay in plain multiallelic case due to A1 column
      logerrputs("Error: --glm's 'test' column cannot be omitted when results for multiple\npredictors are reported.  (Did you forget 'hide-covar'?)\n");
      goto GlmLinear_ret_INCONSISTENT_INPUT;
    }
    const uint32_t main_mutated = ((glm_flags & (kfGlmDominant | kfGlmRecessive | kfGlmHetonly | kfGlmHethom)) != kfGlm0);
    // bugfix (4 Mar 2019): forgot to update this for --tests
    const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!main_mutated) && (!common->tests_flag) && (!add_interactions);
    if (beta_se_multiallelic_fused || (!hide_covar)) {
      max_reported_test_ct += max_extra_allele_ct;
    }
    common->max_reported_test_ct = max_reported_test_ct;

    uint32_t x_code = UINT32_MAXM1;
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    if (sample_ct_x) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetX, &x_code, &x_start, &x_end);
    }
    uint32_t y_code = UINT32_MAXM1;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    if (sample_ct_y) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetY, &y_code, &y_start, &y_end);
    }
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t chr_col = glm_cols & kfGlmColChrom;

    // includes trailing tab
    char* chr_buf = nullptr;
    if (chr_col) {
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto GlmLinear_ret_NOMEM;
      }
    }

    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > variant_ct) {
      calc_thread_ct = variant_ct;
    }

    const uint32_t main_omitted = (parameter_subset && (!IsSet(parameter_subset, 1)));
    const uint32_t xmain_ct = main_mutated + main_omitted;
    uintptr_t workspace_alloc = GetLinearWorkspaceSize(sample_ct, biallelic_predictor_ct, max_extra_allele_ct, constraint_ct, xmain_ct);
    if (sample_ct_x) {
      const uintptr_t workspace_alloc_x = GetLinearWorkspaceSize(sample_ct_x, biallelic_predictor_ct_x, max_extra_allele_ct, constraint_ct_x, xmain_ct);
      if (workspace_alloc_x > workspace_alloc) {
        workspace_alloc = workspace_alloc_x;
      }
    }
    if (sample_ct_y) {
      const uintptr_t workspace_alloc_y = GetLinearWorkspaceSize(sample_ct_y, biallelic_predictor_ct_y, max_extra_allele_ct, constraint_ct_y, xmain_ct);
      if (workspace_alloc_y > workspace_alloc) {
        workspace_alloc = workspace_alloc_y;
      }
    }
    // +1 is for top-level common->workspace_bufs
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;

    uintptr_t per_variant_xalloc_byte_ct = max_sample_ct * local_covar_ct * sizeof(double);
    uintptr_t per_alt_allele_xalloc_byte_ct = sizeof(LinearAuxResult);
    if (beta_se_multiallelic_fused) {
      per_variant_xalloc_byte_ct += 2 * max_reported_test_ct * sizeof(double);
    } else {
      per_alt_allele_xalloc_byte_ct += 2 * max_reported_test_ct * sizeof(double);
    }
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    common->thread_mhc = nullptr;
    common->dosage_presents = nullptr;
    common->dosage_mains = nullptr;
    uint32_t read_block_size;
    uintptr_t max_alt_allele_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, max_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, per_alt_allele_xalloc_byte_ct, pgfip, &calc_thread_ct, &common->genovecs, max_extra_allele_ct? (&common->thread_mhc) : nullptr, nullptr, nullptr, dosage_is_present? (&common->dosage_presents) : nullptr, dosage_is_present? (&common->dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, &max_alt_allele_block_size, main_loadbufs, &common->pgr_ptrs, &common->read_variant_uidx_starts))) {
      goto GlmLinear_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto GlmLinear_ret_NOMEM;
    }
    LinearAuxResult* linear_block_aux_bufs[2];
    double* block_beta_se_bufs[2];

    for (uint32_t uii = 0; uii != 2; ++uii) {
      if (unlikely(BIGSTACK_ALLOC_X(LinearAuxResult, max_alt_allele_block_size, &(linear_block_aux_bufs[uii])))) {
        // shouldn't be possible for these to fail?
        goto GlmLinear_ret_NOMEM;
      }
      if (beta_se_multiallelic_fused) {
        if (unlikely(bigstack_alloc_d(read_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii])))) {
          goto GlmLinear_ret_NOMEM;
        }
      } else {
        if (unlikely(bigstack_alloc_d(max_alt_allele_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii])))) {
          goto GlmLinear_ret_NOMEM;
        }
      }
      if (local_covar_ct) {
        // bugfix (5 Mar 2018): don't want sizeof(double) here
        if (unlikely(bigstack_alloc_d(read_block_size * max_sample_ct * local_covar_ct, &(ctx->local_covars_vcmaj_d[uii])))) {
          goto GlmLinear_ret_NOMEM;
        }
      } else {
        ctx->local_covars_vcmaj_d[uii] = nullptr;
      }
    }

    common->workspace_bufs = S_CAST(unsigned char**, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t)));
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      common->workspace_bufs[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(workspace_alloc));
    }
    common->err_info = (~0LLU) << 32;
    SetThreadFuncAndData(GlmLinearThread, ctx, &tg);

    const uint32_t ref_col = glm_cols & kfGlmColRef;
    const uint32_t alt1_col = glm_cols & kfGlmColAlt1;
    const uint32_t alt_col = glm_cols & kfGlmColAlt;
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t all_nonref = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags, glm_cols / kfGlmColMaybeprovref, raw_variant_ct, all_nonref);
    const uint32_t omitted_col = glm_cols & kfGlmColOmitted;
    const uint32_t ax_col = glm_cols & kfGlmColAx;
    const uint32_t a1_ct_col = glm_cols & kfGlmColA1count;
    const uint32_t tot_allele_col = glm_cols & kfGlmColTotallele;
    const uint32_t a1_freq_col = glm_cols & kfGlmColA1freq;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t beta_col = glm_cols & (kfGlmColBeta | kfGlmColOrbeta);
    const uint32_t se_col = glm_cols & kfGlmColSe;
    const uint32_t ci_col = (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t t_col = glm_cols & kfGlmColTz;
    const uint32_t p_col = glm_cols & kfGlmColP;
    const uint32_t err_col = glm_cols & kfGlmColErr;
    *cswritep++ = '#';
    if (chr_col) {
      cswritep = strcpya_k(cswritep, "CHROM\t");
    }
    if (variant_bps) {
      cswritep = strcpya_k(cswritep, "POS\t");
    }
    cswritep = strcpya_k(cswritep, "ID");
    if (ref_col) {
      cswritep = strcpya_k(cswritep, "\tREF");
    }
    if (alt1_col) {
      cswritep = strcpya_k(cswritep, "\tALT1");
    }
    if (alt_col) {
      cswritep = strcpya_k(cswritep, "\tALT");
    }
    if (provref_col) {
      cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
    }
    cswritep = strcpya_k(cswritep, "\tA1");
    if (omitted_col) {
      cswritep = strcpya_k(cswritep, "\tOMITTED");
    }
    if (ax_col) {
      cswritep = strcpya_k(cswritep, "\tAX");
    }
    if (a1_ct_col) {
      cswritep = strcpya_k(cswritep, "\tA1_CT");
    }
    if (tot_allele_col) {
      cswritep = strcpya_k(cswritep, "\tALLELE_CT");
    }
    if (a1_freq_col) {
      cswritep = strcpya_k(cswritep, "\tA1_FREQ");
    }
    if (mach_r2_col) {
      cswritep = strcpya_k(cswritep, "\tMACH_R2");
    }
    if (test_col) {
      cswritep = strcpya_k(cswritep, "\tTEST");
    }
    if (nobs_col) {
      cswritep = strcpya_k(cswritep, "\tOBS_CT");
    }
    if (beta_col) {
      cswritep = strcpya_k(cswritep, "\tBETA");
    }
    if (se_col) {
      cswritep = strcpya_k(cswritep, "\tSE");
    }
    double ci_zt = 0.0;
    if (ci_col) {
      cswritep = strcpya_k(cswritep, "\tL");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      cswritep = strcpya_k(cswritep, "\tU");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      ci_zt = QuantileToZscore((ci_size + 1.0) * 0.5);
    }
    if (t_col) {
      if (!constraint_ct) {
        cswritep = strcpya_k(cswritep, "\tT_STAT");
      } else {
        // F-statistic for joint tests.
        cswritep = strcpya_k(cswritep, "\tT_OR_F_STAT");
      }
    }
    if (p_col) {
      if (report_neglog10p) {
        // TODO: change to NEG_LOG10_P for a5
        cswritep = strcpya_k(cswritep, "\tLOG10_P");
      } else {
        cswritep = strcpya_k(cswritep, "\tP");
      }
    }
    if (err_col) {
      cswritep = strcpya_k(cswritep, "\tERRCODE");
    }
    AppendBinaryEoln(&cswritep);

    // Main workflow:
    // 1. Set n=0, load/skip block 0
    //
    // 2. Spawn threads processing block n
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Load/skip block n unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8, Write results for last block
    uintptr_t write_variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t suppress_mach_r2 = 0;

    uint32_t cur_biallelic_reported_test_ct = 0;
    uint32_t primary_reported_test_idx = include_intercept;
    uint32_t cur_biallelic_predictor_ct = 0;
    uint32_t cur_constraint_ct = 0;

    const char* const* cur_test_names = nullptr;
    uint32_t prev_block_variant_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t allele_ct = 2;
    uint32_t omitted_allele_idx = 0;
    uintptr_t valid_allele_ct = 0;
    logprintfww5("--glm linear regression on phenotype '%s': ", cur_pheno_name);
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_variant_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto GlmLinear_ret_PGR_FAIL;
      }
      if (local_covar_ct && cur_block_variant_ct) {
        const uint32_t uidx_start = read_block_idx * read_block_size;
        const uint32_t uidx_end = MINV(raw_variant_ct, uidx_start + read_block_size);
        if (local_variant_include) {
          reterr = ReadLocalCovarBlock(common, local_sample_uidx_order, local_variant_include, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, local_covar_txsp, &local_line_idx, &local_xy, nullptr, ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
        } else {
          double* prev_local_covar_row_d = nullptr;
          if (variant_idx) {
            prev_local_covar_row_d = &(ctx->local_covars_vcmaj_d[1 - parity][S_CAST(uintptr_t, read_block_size - 1) * max_sample_ct * local_covar_ct]);
          }
          reterr = ReadRfmix2Block(common, variant_bps, local_sample_uidx_order, nullptr, prev_local_covar_row_d, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, glm_info_ptr->local_chrom_col, glm_info_ptr->local_bp_col, glm_info_ptr->local_first_covar_col, local_covar_txsp, &local_line_iter, &local_line_idx, &local_prev_chr_code, &local_chr_code, &local_bp, &local_skip_chr, nullptr, ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
          /*
          for (uint32_t uii = 0; uii < max_sample_ct; ++uii) {
            printf("%g ", ctx->local_covars_vcmaj_d[parity][uii]);
          }
          printf("\n");
          exit(1);
          */
        }
        if (unlikely(reterr)) {
          goto GlmLinear_ret_1;
        }
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, common->err_info);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, common->err_info >> 32);
          goto GlmLinear_ret_1;
        }
      }
      if (!IsLastBlock(&tg)) {
        common->cur_block_variant_ct = cur_block_variant_ct;
        const uint32_t uidx_start = read_block_idx * read_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_variant_ct, calc_thread_ct, uidx_start, common->read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, common->pgr_ptrs);
        ctx->block_aux = linear_block_aux_bufs[parity];
        common->block_beta_se = block_beta_se_bufs[parity];
        if (variant_idx + cur_block_variant_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto GlmLinear_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const double* beta_se_iter = block_beta_se_bufs[parity];
        const LinearAuxResult* cur_block_aux = linear_block_aux_bufs[parity];
        uintptr_t allele_bidx = 0;
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            if ((chr_idx == x_code) && sample_ct_x) {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct_x;
              cur_biallelic_predictor_ct = biallelic_predictor_ct_x;
              cur_constraint_ct = constraint_ct_x;
              cur_test_names = test_names_x;
            } else if ((chr_idx == y_code) && sample_ct_y) {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct_y;
              cur_biallelic_predictor_ct = biallelic_predictor_ct_y;
              cur_constraint_ct = constraint_ct_y;
              cur_test_names = test_names_y;
            } else {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct;
              cur_biallelic_predictor_ct = biallelic_predictor_ct;
              cur_constraint_ct = constraint_ct;
              cur_test_names = test_names;
            }
            suppress_mach_r2 = (chr_idx == x_code) || (chr_idx == mt_code);
            if (cur_constraint_ct) {
              primary_reported_test_idx = cur_biallelic_reported_test_ct - 1;
            }
            if (chr_col) {
              char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
              *chr_name_end = '\t';
              chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
            }
          }
          uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offsets[write_variant_uidx];
          }
          const uint32_t allele_ct_m1 = allele_ct - 1;
          const uint32_t extra_allele_ct = allele_ct - 2;
          if (omitted_alleles) {
            omitted_allele_idx = omitted_alleles[write_variant_uidx];
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          uint32_t variant_is_valid = 0;
          uint32_t a1_allele_idx = 0;
          for (uint32_t nonomitted_allele_idx = 0; nonomitted_allele_idx != allele_ct_m1; ++nonomitted_allele_idx, ++a1_allele_idx) {
            if (beta_se_multiallelic_fused) {
              if (!nonomitted_allele_idx) {
                primary_reported_test_idx = include_intercept;
              } else {
                primary_reported_test_idx = cur_biallelic_reported_test_ct + nonomitted_allele_idx - 1;
              }
            }
            if (nonomitted_allele_idx == omitted_allele_idx) {
              ++a1_allele_idx;
            }
            const double primary_beta = beta_se_iter[primary_reported_test_idx * 2];
            const double primary_se = beta_se_iter[primary_reported_test_idx * 2 + 1];
            const uint32_t allele_is_valid = (primary_se != -9.0);
            variant_is_valid |= allele_is_valid;
            {
              const LinearAuxResult* auxp = &(cur_block_aux[allele_bidx]);
              if (ln_pfilter <= 0.0) {
                if (!allele_is_valid) {
                  goto GlmLinear_allele_iterate;
                }
                double primary_ln_pval;
                if (!cur_constraint_ct) {
                  if (primary_beta == 0.0) {
                    primary_ln_pval = 0.0;
                  } else if (primary_se == 0.0) {
                    primary_ln_pval = -DBL_MAX;
                  } else {
                    const double primary_tstat = primary_beta / primary_se;
                    primary_ln_pval = TstatToLnP(primary_tstat, auxp->sample_obs_ct - cur_biallelic_predictor_ct - extra_allele_ct);
                  }
                } else {
                  primary_ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                }
                if (primary_ln_pval > ln_pfilter) {
                  if (orig_ln_pvals) {
                    orig_ln_pvals[valid_allele_ct] = primary_ln_pval;
                  }
                  goto GlmLinear_allele_iterate;
                }
              }
              uint32_t inner_reported_test_ct = cur_biallelic_reported_test_ct;
              if (extra_allele_ct) {
                if (beta_se_multiallelic_fused) {
                  if (!nonomitted_allele_idx) {
                    inner_reported_test_ct = 1 + include_intercept;
                  } else if (nonomitted_allele_idx == extra_allele_ct) {
                    inner_reported_test_ct -= include_intercept;
                  } else {
                    inner_reported_test_ct = 1;
                  }
                } else if (!hide_covar) {
                  inner_reported_test_ct += extra_allele_ct;
                }
              }
              // possible todo: make number-to-string operations, strlen(),
              // etc. happen only once per variant.
              for (uint32_t allele_test_idx = 0; allele_test_idx != inner_reported_test_ct; ++allele_test_idx) {
                uint32_t test_idx = allele_test_idx;
                if (beta_se_multiallelic_fused && nonomitted_allele_idx) {
                  if (!allele_test_idx) {
                    test_idx = primary_reported_test_idx;
                  } else {
                    // bugfix (26 Jun 2019): only correct to add 1 here in
                    // include_intercept case
                    test_idx += include_intercept;
                  }
                }
                if (chr_col) {
                  cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
                }
                if (variant_bps) {
                  cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
                }
                cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
                if (ref_col) {
                  *cswritep++ = '\t';
                  cswritep = strcpya(cswritep, cur_alleles[0]);
                }
                if (alt1_col) {
                  *cswritep++ = '\t';
                  cswritep = strcpya(cswritep, cur_alleles[1]);
                }
                if (alt_col) {
                  *cswritep++ = '\t';
                  for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLinear_ret_WRITE_FAIL;
                    }
                    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                  }
                  --cswritep;
                }
                *cswritep++ = '\t';
                if (provref_col) {
                  *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, write_variant_uidx)))? 'Y' : 'N';
                  *cswritep++ = '\t';
                }
                const uint32_t multi_a1 = extra_allele_ct && beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx);
                if (multi_a1) {
                  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    if (allele_idx == omitted_allele_idx) {
                      continue;
                    }
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLinear_ret_WRITE_FAIL;
                    }
                    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                  }
                  --cswritep;
                } else {
                  cswritep = strcpya(cswritep, cur_alleles[a1_allele_idx]);
                }
                if (omitted_col) {
                  *cswritep++ = '\t';
                  cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
                }
                if (ax_col) {
                  *cswritep++ = '\t';
                  if (beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx)) {
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLinear_ret_WRITE_FAIL;
                    }
                    cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
                  } else {
                    for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                      if (allele_idx == a1_allele_idx) {
                        continue;
                      }
                      if (unlikely(Cswrite(&css, &cswritep))) {
                        goto GlmLinear_ret_WRITE_FAIL;
                      }
                      cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                    }
                    --cswritep;
                  }
                }
                if (a1_ct_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_dosage, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (tot_allele_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa(auxp->allele_obs_ct, cswritep);
                }
                if (a1_freq_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_dosage / S_CAST(double, auxp->allele_obs_ct), cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (mach_r2_col) {
                  *cswritep++ = '\t';
                  if (!suppress_mach_r2) {
                    cswritep = dtoa_g(auxp->mach_r2, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (test_col) {
                  *cswritep++ = '\t';
                  if (test_idx < cur_biallelic_reported_test_ct) {
                    cswritep = strcpya(cswritep, cur_test_names[test_idx]);
                  } else {
                    // always use basic dosage for untested alleles
                    cswritep = strcpya_k(cswritep, "ADD");
                    if (!beta_se_multiallelic_fused) {
                      // extra alt allele covariate.
                      uint32_t test_xallele_idx = test_idx - cur_biallelic_reported_test_ct;
                      if (omitted_allele_idx < a1_allele_idx) {
                        test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                      }
                      test_xallele_idx = test_xallele_idx + (test_xallele_idx >= a1_allele_idx);
                      if (a1_allele_idx < omitted_allele_idx) {
                        test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                      }
                      if (!test_xallele_idx) {
                        cswritep = strcpya_k(cswritep, "_REF");
                      } else {
                        cswritep = strcpya_k(cswritep, "_ALT");
                        cswritep = u32toa(test_xallele_idx, cswritep);
                      }
                    }
                  }
                }
                if (nobs_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa(auxp->sample_obs_ct, cswritep);
                }
                double ln_pval = kLnPvalError;
                double tstat = 0.0;
                uint32_t test_is_valid;
                if ((!cur_constraint_ct) || (test_idx != primary_reported_test_idx)) {
                  double beta = beta_se_iter[2 * test_idx];
                  double se = beta_se_iter[2 * test_idx + 1];
                  test_is_valid = (se != -9.0);
                  if (test_is_valid) {
                    tstat = beta / se;
                    ln_pval = TstatToLnP(tstat, auxp->sample_obs_ct - cur_biallelic_predictor_ct - extra_allele_ct);
                  }
                  if (beta_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(beta, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  if (se_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(se, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  if (ci_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      const double ci_halfwidth = ci_zt * se;
                      cswritep = dtoa_g(beta - ci_halfwidth, cswritep);
                      *cswritep++ = '\t';
                      cswritep = dtoa_g(beta + ci_halfwidth, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA\tNA");
                    }
                  }
                  if (t_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(tstat, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                } else {
                  // joint test
                  test_is_valid = allele_is_valid;
                  if (beta_col) {
                    cswritep = strcpya_k(cswritep, "\tNA");
                  }
                  if (se_col) {
                    cswritep = strcpya_k(cswritep, "\tNA");
                  }
                  if (ci_col) {
                    cswritep = strcpya_k(cswritep, "\tNA\tNA");
                  }
                  if (t_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(primary_se / u31tod(cur_constraint_ct), cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  // could avoid recomputing
                  if (test_is_valid) {
                    ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                  }
                }
                if (p_col) {
                  *cswritep++ = '\t';
                  if (test_is_valid) {
                    if (report_neglog10p) {
                      const double reported_val = (-kRecipLn10) * ln_pval;
                      cswritep = dtoa_g(reported_val, cswritep);
                    } else {
                      const double reported_ln = MAXV(ln_pval, output_min_ln);
                      cswritep = lntoa_g(reported_ln, cswritep);
                    }
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (err_col) {
                  *cswritep++ = '\t';
                  if (test_is_valid) {
                    *cswritep++ = '.';
                  } else {
                    uint64_t glm_errcode;
                    memcpy(&glm_errcode, &(beta_se_iter[2 * test_idx]), 8);
                    cswritep = AppendGlmErrstr(glm_errcode, cswritep);
                  }
                }
                AppendBinaryEoln(&cswritep);
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto GlmLinear_ret_WRITE_FAIL;
                }
                if ((test_idx == primary_reported_test_idx) && allele_is_valid) {
                  if (orig_ln_pvals) {
                    orig_ln_pvals[valid_allele_ct] = ln_pval;
                  }
                }
              }
            }
          GlmLinear_allele_iterate:
            ++allele_bidx;
            valid_allele_ct += allele_is_valid;
            if (valid_alleles && allele_is_valid) {
              SetBit(allele_idx_offset_base + a1_allele_idx, valid_alleles);
            }
            if (!beta_se_multiallelic_fused) {
              beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
            }
          }
          if (beta_se_multiallelic_fused) {
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
          if ((!variant_is_valid) && valid_alleles) {
            ClearBit(write_variant_uidx, valid_variants);
          }
        }
      }
      if (variant_idx == variant_ct) {
        break;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }
      ++read_block_idx;
      prev_block_variant_ct = cur_block_variant_ct;
      variant_idx += cur_block_variant_ct;
      // crucially, this is independent of the PgenReader block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto GlmLinear_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    logprintf("Results written to %s .\n", outname);
    *valid_allele_ct_ptr = valid_allele_ct;
  }
  while (0) {
  GlmLinear_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GlmLinear_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
    break;
  GlmLinear_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  GlmLinear_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GlmLinear_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  GlmLinear_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 GlmLinear_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

uintptr_t GetLinearSubbatchWorkspaceSize(uint32_t sample_ct, uint32_t subbatch_size, uint32_t biallelic_predictor_ct, uint32_t max_extra_allele_ct, uint32_t constraint_ct, uint32_t xmain_ct) {
  // sample_ct * max_predictor_ct < 2^31, sample_ct * subbatch_size < 2^31,
  // subbatch_size <= 240, and max_predictor_ct < sqrt(2^31), so no overflows

  // sample_nm, tmp_nm = sample_ctl words
  uintptr_t workspace_size = 2 * RoundUpPow2(BitCtToWordCt(sample_ct) * sizeof(intptr_t), kCacheline);

  // nm_pheno_buf = sample_ct * subbatch_size doubles
  workspace_size += RoundUpPow2(sample_ct * subbatch_size * sizeof(double), kCacheline);

  const uint32_t max_predictor_ct = biallelic_predictor_ct + max_extra_allele_ct;
  // predictors_pmaj = (max_predictor_ct + main_mutated + main_omitted) *
  //                   sample_ct doubles
  workspace_size += RoundUpPow2((max_predictor_ct + xmain_ct) * sample_ct * sizeof(double), kCacheline);

  // xtx_inv, xtx_inv2 = max_predictor_ct * max_predictor_ct doubles
  workspace_size += 2 * RoundUpPow2(max_predictor_ct * max_predictor_ct * sizeof(double), kCacheline);

  // dbl_2d_buf = max_predictor_ct * max(max_predictor_ct, 7) doubles
  workspace_size += RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);

  // inverse_corr_buf = (max_predictor_ct - 1) * max(max_predictor_ct - 1, 4) doubles
  workspace_size += RoundUpPow2((max_predictor_ct - 1) * MAXV((max_predictor_ct - 1), 4) * sizeof(double), kCacheline);

  // semicomputed_biallelic_corr_matrix = (max_predictor_ct - 1)^2 doubles
  workspace_size += RoundUpPow2((biallelic_predictor_ct - 1) * (biallelic_predictor_ct - 1) * sizeof(double), kCacheline);

  // semicomputed_biallelic_inv_corr_sqrts = biallelic_predictor_ct doubles
  workspace_size += RoundUpPow2(biallelic_predictor_ct * sizeof(double), kCacheline);

  // fitted_coefs, xt_y = max_predictor_ct * subbatch_size doubles
  workspace_size += 2 * RoundUpPow2(max_predictor_ct * subbatch_size * sizeof(double), kCacheline);

  // inv_1d_buf
  workspace_size += RoundUpPow2(MAXV(max_predictor_ct, constraint_ct) * kMatrixInvertBuf1CheckedAlloc, kCacheline);

  // machr2_dosage_sums, machr2_dosage_ssqs
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(uint64_t) * 2, kCacheline);

  // pheno_ssq_bases, geno_pheno_prods, domdev_pheno_prods: subbatch_size
  workspace_size += 3 * RoundUpPow2(subbatch_size * sizeof(double), kCacheline);

  if (constraint_ct) {
    // tmphxs_buf, h_transpose_buf = constraint_ct * max_predictor_ct doubles
    workspace_size += 2 * RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(double), kCacheline);

    // inner_buf = constraint_ct * constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * constraint_ct * sizeof(double), kCacheline);

    // cur_constraints_con_major = constraint_ct * max_predictor_ct doubles
    workspace_size += RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(double), kCacheline);
  }
  return workspace_size;
}

THREAD_FUNC_DECL GlmLinearSubbatchThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GlmLinearCtx* ctx = S_CAST(GlmLinearCtx*, arg->sharedp->context);
  GlmCtx* common = ctx->common;

  PgenReader* pgrp = common->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = common->genovecs[tidx];
  pgv.dosage_present = nullptr;
  pgv.dosage_main = nullptr;
  if (common->dosage_presents) {
    pgv.dosage_present = common->dosage_presents[tidx];
    pgv.dosage_main = common->dosage_mains[tidx];
  }
  unsigned char* workspace_buf = common->workspace_bufs[tidx];
  const uintptr_t* variant_include = common->variant_include;
  const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
  const AlleleCode* omitted_alleles = common->omitted_alleles;
  const uintptr_t* sex_male_collapsed = common->sex_male_collapsed;
  const ChrInfo* cip = common->cip;
  const uint32_t* subset_chr_fo_vidx_start = common->subset_chr_fo_vidx_start;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const GlmFlags glm_flags = common->glm_flags;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t model_hetonly = (glm_flags / kfGlmHetonly) & 1;
  const uint32_t joint_genotypic = (glm_flags / kfGlmGenotypic) & 1;
  const uint32_t joint_hethom = (glm_flags / kfGlmHethom) & 1;
  const uint32_t domdev_present = joint_genotypic || joint_hethom;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t reported_pred_uidx_start = 1 - include_intercept;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uintptr_t max_reported_test_ct = common->max_reported_test_ct;
  const uintptr_t local_covar_ct = common->local_covar_ct;
  const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
  const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!model_dominant) && (!model_recessive) && (!model_hetonly) && (!common->tests_flag) && (!add_interactions);
  const uint32_t subbatch_size = ctx->subbatch_size;
  uintptr_t max_sample_ct = MAXV(common->sample_ct, common->sample_ct_x);
  if (max_sample_ct < common->sample_ct_y) {
    max_sample_ct = common->sample_ct_y;
  }
  pgv.patch_01_set = nullptr;
  pgv.patch_01_vals = nullptr;
  pgv.patch_10_set = nullptr;
  pgv.patch_10_vals = nullptr;
  if (common->thread_mhc) {
    const uint32_t max_sample_ctl = BitCtToWordCt(max_sample_ct);
    pgv.patch_01_set = common->thread_mhc[tidx];
    pgv.patch_01_vals = R_CAST(AlleleCode*, &(pgv.patch_01_set[max_sample_ctl]));
    AlleleCode* patch_01_vals_end = &(pgv.patch_01_vals[max_sample_ct]);
    AlignACToVec(&patch_01_vals_end);
    pgv.patch_10_set = R_CAST(uintptr_t*, patch_01_vals_end);
    pgv.patch_10_vals = R_CAST(AlleleCode*, &(pgv.patch_10_set[max_sample_ctl]));
  }
  pgv.patch_01_ct = 0;
  pgv.patch_10_ct = 0;
  pgv.multidosage_sample_ct = 0;
  uint32_t variant_idx_offset = 0;
  uint32_t allele_ct = 2;
  uint32_t omitted_allele_idx = 0;
  uint32_t extra_regression_ct = 0;
  double main_dosage_sum = 0.0;
  double main_dosage_ssq = 0.0;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);

    uintptr_t allele_bidx = variant_bidx;
    if (max_extra_allele_ct) {
      allele_bidx = variant_bidx + CountExtraAlleles(variant_include, allele_idx_offsets, common->read_variant_uidx_starts[0], common->read_variant_uidx_starts[tidx], 0);
    }
    double* beta_se_iter = common->block_beta_se;
    if (beta_se_multiallelic_fused) {
      beta_se_iter = &(beta_se_iter[(2 * k1LU * subbatch_size) * max_reported_test_ct * variant_bidx]);
    } else {
      beta_se_iter = &(beta_se_iter[(2 * k1LU * subbatch_size) * max_reported_test_ct * allele_bidx]);
    }

    LinearAuxResult* block_aux_iter = &(ctx->block_aux[allele_bidx]);
    const double* local_covars_iter = nullptr;
    if (local_covar_ct) {
      // &(nullptr[0]) is okay in C++, but undefined in C
      local_covars_iter = &(ctx->local_covars_vcmaj_d[parity][variant_bidx * max_sample_ct * local_covar_ct]);
    }
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = CountSortedSmallerU32(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      uint32_t cur_variant_bidx_end = subset_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
      if (cur_variant_bidx_end > variant_bidx_end) {
        cur_variant_bidx_end = variant_bidx_end;
      }
      const uint32_t is_haploid = IsSet(cip->haploid_mask, chr_idx);
      const uint32_t is_regular_x = is_haploid && (chr_idx == x_code);
      const uint32_t is_y = (chr_idx == y_code);
      const uint32_t is_nonx_haploid = is_haploid && (!is_regular_x);
      const uintptr_t* cur_sample_include;
      const uint32_t* cur_sample_include_cumulative_popcounts;
      const double* cur_pheno_pmaj;
      const RegressionNmPrecomp* nm_precomp;
      const double* cur_covars_cmaj;
      const uintptr_t* cur_parameter_subset;
      const uintptr_t* cur_joint_test_params;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      uint32_t cur_constraint_ct;
      if (is_y && common->sample_include_y) {
        cur_sample_include = common->sample_include_y;
        cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
        cur_pheno_pmaj = ctx->pheno_y_d;
        nm_precomp = common->nm_precomp_y;
        cur_covars_cmaj = ctx->covars_cmaj_y_d;
        cur_parameter_subset = common->parameter_subset_y;
        cur_joint_test_params = common->joint_test_params_y;
        cur_sample_ct = common->sample_ct_y;
        cur_covar_ct = common->covar_ct_y;
        cur_constraint_ct = common->constraint_ct_y;
      } else if (is_regular_x && common->sample_include_x) {
        cur_sample_include = common->sample_include_x;
        cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
        cur_pheno_pmaj = ctx->pheno_x_d;
        nm_precomp = common->nm_precomp_x;
        cur_covars_cmaj = ctx->covars_cmaj_x_d;
        cur_parameter_subset = common->parameter_subset_x;
        cur_joint_test_params = common->joint_test_params_x;
        cur_sample_ct = common->sample_ct_x;
        cur_covar_ct = common->covar_ct_x;
        cur_constraint_ct = common->constraint_ct_x;
      } else {
        cur_sample_include = common->sample_include;
        cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
        cur_pheno_pmaj = ctx->pheno_d;
        nm_precomp = common->nm_precomp;
        cur_covars_cmaj = ctx->covars_cmaj_d;
        cur_parameter_subset = common->parameter_subset;
        cur_joint_test_params = common->joint_test_params;
        cur_sample_ct = common->sample_ct;
        cur_covar_ct = common->covar_ct;
        cur_constraint_ct = common->constraint_ct;
      }
      const uint32_t sample_ctl = BitCtToWordCt(cur_sample_ct);
      const uint32_t sample_ctl2 = NypCtToWordCt(cur_sample_ct);
      const uint32_t cur_biallelic_predictor_ct_base = 2 + domdev_present + cur_covar_ct * (1 + add_interactions * domdev_present_p1);
      uint32_t cur_biallelic_predictor_ct = cur_biallelic_predictor_ct_base;
      uint32_t literal_covar_ct = cur_covar_ct;
      if (cur_parameter_subset) {
        cur_biallelic_predictor_ct = PopcountWords(cur_parameter_subset, BitCtToWordCt(cur_biallelic_predictor_ct_base));
        literal_covar_ct = PopcountBitRange(cur_parameter_subset, 2 + domdev_present, 2 + domdev_present + cur_covar_ct);
      }
      const uint32_t max_predictor_ct = cur_biallelic_predictor_ct + max_extra_allele_ct;
      uint32_t reported_pred_uidx_biallelic_end;
      if (hide_covar) {
        if (!cur_parameter_subset) {
          reported_pred_uidx_biallelic_end = 2 + domdev_present;
        } else {
          reported_pred_uidx_biallelic_end = IsSet(cur_parameter_subset, 1) + domdev_present_p1;
        }
      } else {
        reported_pred_uidx_biallelic_end = cur_biallelic_predictor_ct;
      }
      // nm_predictors_pmaj_buf may require up to two extra columns omitted
      // from the main regression.
      // 1. In the multiallelic dominant/recessive/hetonly/hethom cases, the
      //    original genotype column does not appear in the regression, and
      //    we'd rather not reconstruct it from genovec, etc. when we need to
      //    swap it out for another allele, so we keep the original genotype in
      //    an extra column.
      //    To reduce code bloat, we now handle the biallelic cases in the same
      //    way; this is one of the more peripheral code paths so adding more
      //    complexity to speed it up is less justifiable.
      // 2. If --parameters excludes the main (possibly
      //    dominant/recessive/hetonly) genotype column but does care about an
      //    interaction, we want a copy of what the main genotype column's
      //    contents would have been to refer to.
      const uint32_t main_omitted = cur_parameter_subset && (!IsSet(cur_parameter_subset, 1));
      const uint32_t main_mutated = model_dominant || model_recessive || model_hetonly || joint_hethom;
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* tmp_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      double* nm_pheno_buf = S_CAST(double*, arena_alloc_raw_rd(cur_sample_ct * subbatch_size * sizeof(double), &workspace_iter));
      double* nm_predictors_pmaj_buf = S_CAST(double*, arena_alloc_raw_rd((max_predictor_ct + main_mutated + main_omitted) * cur_sample_ct * sizeof(double), &workspace_iter));
      double* xtx_inv = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ct * sizeof(double), &workspace_iter));
      double* fitted_coefs = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * subbatch_size * sizeof(double), &workspace_iter));
      double* xt_y = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * subbatch_size * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_corr_matrix = S_CAST(double*, arena_alloc_raw_rd((cur_biallelic_predictor_ct - 1) * (cur_biallelic_predictor_ct - 1) * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_inv_corr_sqrts = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));
      MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, arena_alloc_raw_rd(MAXV(max_predictor_ct, cur_constraint_ct) * kMatrixInvertBuf1CheckedAlloc, &workspace_iter));
      double* dbl_2d_buf = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), &workspace_iter));
      uint64_t* machr2_dosage_sums = S_CAST(uint64_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(uint64_t) * 2, &workspace_iter));
      uint64_t* machr2_dosage_ssqs = &(machr2_dosage_sums[max_extra_allele_ct + 2]);
      double* pheno_ssq_bases = S_CAST(double*, arena_alloc_raw_rd(subbatch_size * sizeof(double), &workspace_iter));
      double* geno_pheno_prods = S_CAST(double*, arena_alloc_raw_rd(subbatch_size * sizeof(double), &workspace_iter));
      double* domdev_pheno_prods = S_CAST(double*, arena_alloc_raw_rd(subbatch_size * sizeof(double), &workspace_iter));
      double* xtx_inv2 = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ct * sizeof(double), &workspace_iter));

      // could technically have this overlap fitted_coefs/xt_y, but that sets
      // the stage for future bugs
      double* inverse_corr_buf = S_CAST(double*, arena_alloc_raw_rd((max_predictor_ct - 1) * MAXV((max_predictor_ct - 1), 4) * sizeof(double), &workspace_iter));

      // joint test only
      double* tmphxs_buf = nullptr;
      double* h_transpose_buf = nullptr;
      double* inner_buf = nullptr;
      double* cur_constraints_con_major = nullptr;
      if (cur_constraint_ct) {
        tmphxs_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        h_transpose_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        inner_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * cur_constraint_ct * sizeof(double), &workspace_iter));
        cur_constraints_con_major = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        ZeroDArr(cur_constraint_ct * max_predictor_ct, cur_constraints_con_major);
        const uint32_t first_joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
        cur_constraints_con_major[first_joint_test_idx] = 1.0;
        // Rest of this matrix must be updated later, since cur_predictor_ct
        // changes at multiallelic variants.
      }
      assert(S_CAST(uintptr_t, workspace_iter - workspace_buf) == GetLinearSubbatchWorkspaceSize(cur_sample_ct, subbatch_size, cur_biallelic_predictor_ct, max_extra_allele_ct, cur_constraint_ct, main_mutated + main_omitted));
      for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
        const double* cur_pheno = &(cur_pheno_pmaj[pheno_idx * cur_sample_ct]);
        pheno_ssq_bases[pheno_idx] = DotprodD(cur_pheno, cur_pheno, cur_sample_ct);
      }
      const double cur_sample_ct_recip = 1.0 / u31tod(cur_sample_ct);
      const double cur_sample_ct_m1_recip = 1.0 / u31tod(cur_sample_ct - 1);
      const uint32_t sparse_optimization_eligible = (!is_regular_x) && nm_precomp;
      double geno_d_lookup[2];
      if (sparse_optimization_eligible) {
        if (model_hetonly) {
          geno_d_lookup[0] = 1.0;
          geno_d_lookup[1] = 0.0;
        } else {
          geno_d_lookup[1] = 1.0;
          if (is_nonx_haploid) {
            geno_d_lookup[0] = 0.5;
          } else if (model_recessive || joint_hethom) {
            geno_d_lookup[0] = 0.0;
          } else {
            geno_d_lookup[0] = 1.0;
            if (!model_dominant) {
              geno_d_lookup[1] = 2.0;
            }
          }
        }
      }
      const double* xtx_image = nullptr;
      const double* covarx_dotprod_inv = nullptr;
      const double* corr_inv = nullptr;
      const double* xt_y_image = nullptr;
      if (nm_precomp) {
        xtx_image = nm_precomp->xtx_image;
        covarx_dotprod_inv = nm_precomp->covarx_dotprod_inv;
        corr_inv = nm_precomp->corr_inv;
        const uintptr_t nongeno_pred_ct = cur_biallelic_predictor_ct - domdev_present - 2;
        const uintptr_t nonintercept_biallelic_pred_ct = cur_biallelic_predictor_ct - 1;
        memcpy(semicomputed_biallelic_corr_matrix, nm_precomp->corr_image, nonintercept_biallelic_pred_ct * nonintercept_biallelic_pred_ct * sizeof(double));
        memcpy(&(semicomputed_biallelic_inv_corr_sqrts[domdev_present_p1]), nm_precomp->corr_inv_sqrts, nongeno_pred_ct * sizeof(double));
        xt_y_image = nm_precomp->xt_y_image;
      }
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      // when this is set, the last fully-processed variant had no missing
      // genotypes, and if the current variant also has no missing genotypes we
      // may be able to skip reinitialization of most of
      // nm_predictors_pmaj_buf.
      uint32_t prev_nm = 0;

      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        DPrintf("\nWorker thread %" PRIuPTR " starting unfiltered variant %u", tidx, variant_uidx);
        if (allele_idx_offsets) {
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
          if (!beta_se_multiallelic_fused) {
            extra_regression_ct = allele_ct - 2;
          }
        }
        const uint32_t allele_ct_m2 = allele_ct - 2;
        const uint32_t expected_predictor_ct = cur_biallelic_predictor_ct + allele_ct_m2;
        PglErr reterr;
        if (!allele_ct_m2) {
          reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
        } else {
          reterr = PgrGetMD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, &pgv);
          // todo: proper multiallelic dosage support
        }
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto GlmLinearSubbatchThread_err;
        }
        ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
        GenoarrCountFreqsUnsafe(pgv.genovec, cur_sample_ct, genocounts);
        uint32_t missing_ct = genocounts[3];
        if (!missing_ct) {
          SetAllBits(cur_sample_ct, sample_nm);
        } else {
          GenoarrToNonmissing(pgv.genovec, cur_sample_ct, sample_nm);
          if (pgv.dosage_ct) {
            BitvecOr(pgv.dosage_present, sample_ctl, sample_nm);
            missing_ct = cur_sample_ct - PopcountWords(sample_nm, sample_ctl);
          }
        }
        if (omitted_alleles) {
          omitted_allele_idx = omitted_alleles[variant_uidx];
        }
        uintptr_t const_alleles[DivUp(kPglMaxAlleleCt, kBitsPerWord)];
        const uint32_t allele_ctl = DivUp(allele_ct, kBitsPerWord);
        ZeroWArr(allele_ctl, const_alleles);
        const uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
        const uint32_t nm_sample_ctl = BitCtToWordCt(nm_sample_ct);
        // first predictor column: intercept
        if (!prev_nm) {
          for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
            nm_predictors_pmaj_buf[sample_idx] = 1.0;
          }
        }
        // second predictor column: main genotype
        double* genotype_vals = &(nm_predictors_pmaj_buf[nm_sample_ct]);
        if (main_mutated || main_omitted) {
          genotype_vals = &(nm_predictors_pmaj_buf[expected_predictor_ct * nm_sample_ct]);
        }
        uint32_t sparse_optimization = 0;
        double* multi_start = nullptr;
        if (!allele_ct_m2) {
          // When prev_nm is set and missing_ct is zero, we don't need to call
          // MultiplySelfTranspose() on nm_predictors_pmaj_buf to get all the
          // predictor x predictor dot products; instead we patch in the
          // genotype x predictor dot products that may change, copying the
          // rest from xtx_image.
          // As a side effect, it is no longer strictly necessary to fill the
          // genotype row of nm_predictors_pmaj_buf.  sparse_optimization
          // indicates that plink 1.9's QT --assoc sparse dot product algorithm
          // will be used instead.
          // probable todos: allow a few dosages to be present, cover chrX
          // case.
          if (omitted_allele_idx) {
            GenovecInvertUnsafe(cur_sample_ct, pgv.genovec);
            ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
            }
            const uint32_t uii = genocounts[0];
            genocounts[0] = genocounts[2];
            genocounts[2] = uii;
          }
          uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
          uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
          if (!missing_ct) {
            // originally had genocounts[0] > 0.875 * nm_sample_ct threshold,
            // but then tried this on high-MAF data and it was still
            // substantially faster
            sparse_optimization = sparse_optimization_eligible && (!pgv.dosage_ct) && prev_nm;
            if (!sparse_optimization) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, genotype_vals);
              if (pgv.dosage_ct) {
                uintptr_t sample_idx_base = 0;
                uintptr_t dosage_present_bits = pgv.dosage_present[0];
                for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                  const uintptr_t sample_idx = BitIter1(pgv.dosage_present, &sample_idx_base, &dosage_present_bits);
                  const uintptr_t dosage_val = pgv.dosage_main[dosage_idx];
                  // 32768 -> 2, 16384 -> 1, 0 -> 0
                  genotype_vals[sample_idx] = kRecipDosageMid * swtod(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_idx);
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                }
              }
            }
          } else {
            if (!pgv.dosage_ct) {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, genotype_vals);
            } else {
              uintptr_t sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              uint32_t dosage_idx = 0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_midx);
                double cur_val;
                if (IsSet(pgv.dosage_present, sample_midx)) {
                  const uintptr_t dosage_val = pgv.dosage_main[dosage_idx++];
                  cur_val = kRecipDosageMid * swtod(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                } else {
                  cur_val = kSmallDoubles[cur_geno];
                }
                genotype_vals[sample_idx] = cur_val;
              }
            }
          }
          // Check for constant genotype column.
          // (Technically, we should recheck later in the chrX no-sex-covariate
          // --xchr-model 1 corner case.)
          if (!pgv.dosage_ct) {
            if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          } else if (pgv.dosage_ct == nm_sample_ct) {
            if (DosageIsConstant(dosage_sum, dosage_ssq, nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          }
          machr2_dosage_sums[1 - omitted_allele_idx] = dosage_sum;
          machr2_dosage_ssqs[1 - omitted_allele_idx] = dosage_ssq;
          machr2_dosage_sums[omitted_allele_idx] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[omitted_allele_idx] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
        } else {
          // multiallelic.
          // If some but not all alleles have constant dosages, we remove just
          // those alleles from the regressions; trim-alts is not necessary to
          // see what's going on with the other alleles.  To reduce parsing
          // complexity, the number of output lines is not affected by this;
          // the ones corresponding to the constant alleles have NA values.
          // Punt on sparse_optimization for now; may be worth revisiting after
          // multiallelic dosage implemented.
          // dosage_ct == 0 temporarily guaranteed if we reach here.
          assert(!pgv.dosage_ct);
          multi_start = &(nm_predictors_pmaj_buf[(expected_predictor_ct - allele_ct_m2) * nm_sample_ct]);
          ZeroU64Arr(allele_ct, machr2_dosage_sums);
          ZeroU64Arr(allele_ct, machr2_dosage_ssqs);
          // postpone multiply for now, since no multiallelic dosages
          machr2_dosage_sums[0] = genocounts[1] + 2 * genocounts[0];
          machr2_dosage_ssqs[0] = genocounts[1] + 4LLU * genocounts[0];
          if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
            SetBit(0, const_alleles);
          }
          if (omitted_allele_idx) {
            // Main genotype column starts as REF.
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallInvDoublePairs, nm_sample_ct, genotype_vals);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallInvDoubles, cur_sample_ct, genotype_vals);
            }
          }
          uint32_t rare_allele_ct = allele_ct_m2;
          double* alt1_start = nullptr;
          double* rarealt_start = multi_start;
          if (omitted_allele_idx != 1) {
            if (omitted_allele_idx) {
              alt1_start = multi_start;
              rarealt_start = &(rarealt_start[nm_sample_ct]);
              --rare_allele_ct;
            } else {
              alt1_start = genotype_vals;
            }
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, alt1_start);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, alt1_start);
            }
          }
          ZeroDArr(rare_allele_ct * nm_sample_ct, rarealt_start);
          // Use sums as ones[] and ssqs as twos[] for rarealts; transform to
          // actual sums/ssqs later.
          if (pgv.patch_01_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_01_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_01_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ct + sample_idx] = 1.0;
                alt1_start[sample_idx] = 0.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ct + sample_idx] = 1.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                alt1_start[sample_idx] = 0.0;
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                machr2_dosage_sums[allele_code] += 1;
                if (allele_code == omitted_allele_idx) {
                  continue;
                }
                const uint32_t cur_col = allele_code - 2 - (allele_code > omitted_allele_idx);
                rarealt_start[cur_col * nm_sample_ct + sample_idx] = 1.0;
              }
            }
          }
          uintptr_t alt1_het_ct = genocounts[1] - pgv.patch_01_ct;
          if (pgv.patch_10_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_10_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_10_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] = 2.0;
                  alt1_start[sample_idx] = 0.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ct + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] += 1.0;
                    alt1_start[sample_idx] = 0.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] = 2.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ct + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] += 1.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  machr2_dosage_ssqs[ac0] += 1;
                  alt1_start[sample_idx] = 0.0;
                  if (ac0 != omitted_allele_idx) {
                    const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                    rarealt_start[ac0_col * nm_sample_ct + sample_idx] = 2.0;
                  }
                } else {
                  machr2_dosage_sums[ac1] += 1;
                  if (ac1 != omitted_allele_idx) {
                    const uint32_t ac1_col = ac1 - 2 - (ac1 > omitted_allele_idx);
                    rarealt_start[ac1_col * nm_sample_ct + sample_idx] = 1.0;
                  }
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    machr2_dosage_sums[ac0] += 1;
                    alt1_start[sample_idx] = 0.0;
                    if (ac0 != omitted_allele_idx) {
                      const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                      rarealt_start[ac0_col * nm_sample_ct + sample_idx] += 1.0;
                    }
                  }
                }
              }
            }
          }
          for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
            const uintptr_t one_ct = machr2_dosage_sums[allele_idx];
            const uintptr_t two_ct = machr2_dosage_ssqs[allele_idx];
            machr2_dosage_sums[allele_idx] = one_ct + 2 * two_ct;
            machr2_dosage_ssqs[allele_idx] = one_ct + 4LLU * two_ct;
            if ((one_ct == nm_sample_ct) || (two_ct == nm_sample_ct) || ((!one_ct) && (!two_ct))) {
              SetBit(allele_idx, const_alleles);
            }
          }
          const uintptr_t alt1_hom_ct = genocounts[2] - pgv.patch_10_ct;
          machr2_dosage_sums[1] = alt1_het_ct + 2 * alt1_hom_ct;
          machr2_dosage_ssqs[1] = alt1_het_ct + 4LLU * alt1_hom_ct;
          if ((alt1_het_ct == nm_sample_ct) || (alt1_hom_ct == nm_sample_ct) || ((!alt1_het_ct) && (!alt1_hom_ct))) {
            SetBit(1, const_alleles);
          }
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            machr2_dosage_sums[allele_idx] *= 0x4000LLU;
            machr2_dosage_ssqs[allele_idx] *= 0x10000000LLU;
          }
        }
        // usually need to save some of {sample_obs_ct, allele_obs_ct,
        // a1_dosage, mach_r2 even for skipped variants
        // compute them all for now, could conditionally skip later
        uint32_t allele_obs_ct = nm_sample_ct * 2;
        if (!is_regular_x) {
          if (is_nonx_haploid) {
            allele_obs_ct = nm_sample_ct;
            // everything is on 0..1 scale, not 0..2
            if (!sparse_optimization) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                genotype_vals[sample_idx] *= 0.5;
              }
              const uint32_t high_ct = nm_sample_ct * allele_ct_m2;
              for (uint32_t uii = 0; uii != high_ct; ++uii) {
                multi_start[uii] *= 0.5;
              }
            }
          }
        } else {
          CopyBitarrSubset(sex_male_collapsed, sample_nm, nm_sample_ct, tmp_nm);
          const uintptr_t* male_nm = tmp_nm;
          const uint32_t nm_male_ct = PopcountWords(male_nm, nm_sample_ctl);
          if (is_xchr_model_1) {
            // special case: multiply male values by 0.5
            uintptr_t sample_idx_base = 0;
            uintptr_t male_nm_bits = male_nm[0];
            for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
              const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
              genotype_vals[sample_idx] *= 0.5;
              // could insert multiallelic loop here instead, but I'm guessing
              // that's worse due to locality of writes?
            }
            for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
              double* cur_start = &(multi_start[extra_allele_idx * nm_sample_ct]);
              sample_idx_base = 0;
              male_nm_bits = male_nm[0];
              for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
                const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
                cur_start[sample_idx] *= 0.5;
              }
            }
            allele_obs_ct -= nm_male_ct;
          }
        }
        const double mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, allele_ct);
        for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
          if (allele_idx == omitted_allele_idx) {
            continue;
          }
          block_aux_iter->sample_obs_ct = nm_sample_ct;
          block_aux_iter->allele_obs_ct = allele_obs_ct;
          double a1_dosage = u63tod(machr2_dosage_sums[allele_idx]) * kRecipDosageMid;
          if (is_xchr_model_1) {
            // ugh.
            double* geno_col = genotype_vals;
            if (allele_idx > (!omitted_allele_idx)) {
              geno_col = &(nm_predictors_pmaj_buf[(expected_predictor_ct - (allele_ct - allele_idx) + (allele_idx < omitted_allele_idx)) * nm_sample_ct]);
            }
            a1_dosage = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              a1_dosage += geno_col[sample_idx];
            }
            if (!allele_ct_m2) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = 0.0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const double cur_dosage = geno_col[sample_idx];
                main_dosage_ssq += cur_dosage * cur_dosage;
              }
            }
          } else {
            if (is_nonx_haploid) {
              a1_dosage *= 0.5;
            }
            if (!allele_ct_m2) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = u63tod(machr2_dosage_ssqs[allele_idx]) * kRecipDosageMidSq;
              if (is_nonx_haploid) {
                main_dosage_ssq *= 0.25;
              }
            }
          }
          block_aux_iter->a1_dosage = a1_dosage;
          block_aux_iter->mach_r2 = mach_r2;
          ++block_aux_iter;
        }
        // now free to skip the actual regression if there are too few samples,
        // or there's a zero-variance genotype column
        GlmErr glm_err = 0;
        if (nm_sample_ct <= expected_predictor_ct) {
          glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
        } else if (IsSet(const_alleles, omitted_allele_idx)) {
          glm_err = SetGlmErr0(kGlmErrcodeConstOmittedAllele);
        }
        if (glm_err) {
          if (missing_ct) {
            // covariates have not been copied yet, so we can't usually change
            // prev_nm from 0 to 1 when missing_ct == 0 (and there's little
            // reason to optimize the zero-covariate case).
            prev_nm = 0;
          }
          uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
          if (beta_se_multiallelic_fused || (!hide_covar)) {
            reported_ct += allele_ct_m2;
          }
          const uintptr_t result_row_ct = (extra_regression_ct + 1) * subbatch_size;
          for (uintptr_t ulii = 0; ulii != result_row_ct; ++ulii) {
            for (uint32_t uii = 0; uii != reported_ct; ++uii) {
              memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
              beta_se_iter[uii * 2 + 1] = -9.0;
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        } else {
          uint32_t parameter_uidx = 2 + domdev_present;
          double* nm_predictors_pmaj_istart = nullptr;
          if (!sparse_optimization) {
            // only need to do this part once per variant in multiallelic case
            double* nm_predictors_pmaj_iter = &(nm_predictors_pmaj_buf[nm_sample_ct * (parameter_uidx - main_omitted)]);
            if (missing_ct || (!prev_nm)) {
              // fill phenotypes
              for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                const double* cur_pheno = &(cur_pheno_pmaj[pheno_idx * cur_sample_ct]);
                double* nm_pheno_col = &(nm_pheno_buf[pheno_idx * nm_sample_ct]);
                uintptr_t sample_midx_base = 0;
                uintptr_t sample_nm_bits = sample_nm[0];
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                  nm_pheno_col[sample_idx] = cur_pheno[sample_midx];
                }
              }

              // fill covariates
              for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx, ++parameter_uidx) {
                // strictly speaking, we don't need cur_covars_cmaj to be
                // vector-aligned
                if (cur_parameter_subset && (!IsSet(cur_parameter_subset, parameter_uidx))) {
                  continue;
                }
                const double* cur_covar_col;
                if (covar_idx < local_covar_ct) {
                  cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                } else {
                  cur_covar_col = &(cur_covars_cmaj[(covar_idx - local_covar_ct) * cur_sample_ct]);
                }
                uintptr_t sample_midx_base = 0;
                uintptr_t sample_nm_bits = sample_nm[0];
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                  *nm_predictors_pmaj_iter++ = cur_covar_col[sample_midx];
                }
              }
              nm_predictors_pmaj_istart = nm_predictors_pmaj_iter;
              prev_nm = !(missing_ct || local_covar_ct);
            } else {
              // bugfix (15 Aug 2018): this was not handling --parameters
              // correctly when a covariate was only needed as part of an
              // interaction
              parameter_uidx += cur_covar_ct;
              nm_predictors_pmaj_istart = &(nm_predictors_pmaj_iter[literal_covar_ct * nm_sample_ct]);
            }
          }
          const uint32_t const_allele_ct = PopcountWords(const_alleles, allele_ctl);
          if (const_allele_ct) {
            // Must delete constant-allele columns from nm_predictors_pmaj, and
            // shift later columns back.
            double* read_iter = genotype_vals;
            double* write_iter = genotype_vals;
            for (uint32_t read_allele_idx = 0; read_allele_idx != allele_ct; ++read_allele_idx) {
              if (read_allele_idx == omitted_allele_idx) {
                continue;
              }
              if (!IsSet(const_alleles, read_allele_idx)) {
                if (write_iter != read_iter) {
                  memcpy(write_iter, read_iter, nm_sample_ct * sizeof(double));
                }
                if (write_iter == genotype_vals) {
                  write_iter = multi_start;
                } else {
                  write_iter = &(write_iter[nm_sample_ct]);
                }
              }
              if (read_iter == genotype_vals) {
                read_iter = multi_start;
              } else {
                read_iter = &(read_iter[nm_sample_ct]);
              }
            }
          }
          const uint32_t cur_predictor_ct = expected_predictor_ct - const_allele_ct;
          uint32_t nonconst_extra_regression_idx = UINT32_MAX;  // deliberate overflow
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            if (extra_regression_ct) {
              if (IsSet(const_alleles, extra_regression_idx + (extra_regression_idx >= omitted_allele_idx))) {
                glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
                goto GlmLinearSubbatchThread_skip_regression;
              }
              ++nonconst_extra_regression_idx;
              if (nonconst_extra_regression_idx) {
                double* swap_target = &(multi_start[(nonconst_extra_regression_idx - 1) * nm_sample_ct]);
                for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
                  double dxx = genotype_vals[uii];
                  genotype_vals[uii] = swap_target[uii];
                  swap_target[uii] = dxx;
                }
              }
            }
            if (!sparse_optimization) {
              double* main_vals = &(nm_predictors_pmaj_buf[nm_sample_ct]);
              double* domdev_vals = nullptr;
              if (main_omitted) {
                // if main_mutated, this will be filled below
                // if not, this aliases genotype_vals
                main_vals = &(nm_predictors_pmaj_buf[(cur_predictor_ct + main_mutated) * nm_sample_ct]);
              } else if (joint_genotypic || joint_hethom) {
                domdev_vals = &(main_vals[nm_sample_ct]);
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 2.0 - cur_genotype_val;
                  }
                  domdev_vals[sample_idx] = cur_genotype_val;
                }
              }
              if (model_dominant) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..1..1
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 1.0;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              } else if (model_recessive || joint_hethom) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..0..1
                  if (cur_genotype_val < 1.0) {
                    cur_genotype_val = 0.0;
                  } else {
                    cur_genotype_val -= 1.0;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              } else if (model_hetonly) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..1..0
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 2.0 - cur_genotype_val;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              }

              // fill interaction terms
              if (add_interactions) {
                double* nm_predictors_pmaj_iter = nm_predictors_pmaj_istart;
                for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
                  const double* cur_covar_col;
                  if (covar_idx < local_covar_ct) {
                    cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                  } else {
                    cur_covar_col = &(cur_covars_cmaj[covar_idx * cur_sample_ct]);
                  }
                  if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                    uintptr_t sample_midx_base = 0;
                    uintptr_t sample_nm_bits = sample_nm[0];
                    for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                      const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                      *nm_predictors_pmaj_iter++ = main_vals[sample_idx] * cur_covar_col[sample_midx];
                    }
                  }
                  ++parameter_uidx;
                  if (domdev_present) {
                    if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                      uintptr_t sample_midx_base = 0;
                      uintptr_t sample_nm_bits = sample_nm[0];
                      for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                        const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                        *nm_predictors_pmaj_iter++ = domdev_vals[sample_idx] * cur_covar_col[sample_midx];
                      }
                    }
                    ++parameter_uidx;
                  }
                }
              }
            }

            // bugfix (12 Sep 2017): forgot to implement per-variant VIF and
            // max-corr checks
            DPrintf("\nWorker thread %" PRIuPTR " starting main regression for %u", tidx, variant_uidx);
            if (xtx_image && prev_nm && (!allele_ct_m2)) {
              // only need to fill in additive and possibly domdev dot
              // products
              memcpy(xtx_inv, xtx_image, cur_predictor_ct * cur_predictor_ct * sizeof(double));
              memcpy(xt_y, xt_y_image, cur_predictor_ct * subbatch_size * sizeof(double));
              if (sparse_optimization) {
                // currently does not handle chrX
                ZeroDArr(subbatch_size, geno_pheno_prods);
                if (domdev_pheno_prods) {
                  ZeroDArr(subbatch_size, domdev_pheno_prods);
                }
                double domdev_geno_prod = 0.0;
                double* geno_dotprod_row = &(xtx_inv[cur_predictor_ct]);
                double* domdev_dotprod_row = &(xtx_inv[2 * cur_predictor_ct]);
                for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
                  uintptr_t geno_word = pgv.genovec[widx];
                  if (geno_word) {
                    const uint32_t sample_idx_base = widx * kBitsPerWordD2;
                    do {
                      const uint32_t lowest_set_bit = ctzw(geno_word);
                      // since there are no missing values, we have a het if
                      // (lowest_set_bit & 1) is zero, and a hom-alt when
                      // it's one.
                      const uint32_t sample_idx = sample_idx_base + (lowest_set_bit / 2);
                      const double geno_d = geno_d_lookup[lowest_set_bit & 1];
                      for (uintptr_t pred_idx = domdev_present + 2; pred_idx != cur_predictor_ct; ++pred_idx) {
                        geno_dotprod_row[pred_idx] += geno_d * nm_predictors_pmaj_buf[pred_idx * nm_sample_ct + sample_idx];
                      }
                      // can have a separate categorical loop here

                      if (domdev_present && (!(lowest_set_bit & 1))) {
                        // domdev = 1
                        for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                          const double cur_pheno_val = nm_pheno_buf[pheno_idx * nm_sample_ct + sample_idx];
                          geno_pheno_prods[pheno_idx] += geno_d * cur_pheno_val;
                          domdev_pheno_prods[pheno_idx] += cur_pheno_val;
                        }
                        domdev_geno_prod += geno_d;
                        for (uintptr_t pred_idx = 3; pred_idx != cur_predictor_ct; ++pred_idx) {
                          domdev_dotprod_row[pred_idx] += nm_predictors_pmaj_buf[pred_idx * nm_sample_ct + sample_idx];
                        }
                        // categorical optimization possible here
                      } else {
                        for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                          geno_pheno_prods[pheno_idx] += geno_d * nm_pheno_buf[pheno_idx * nm_sample_ct + sample_idx];
                        }
                      }
                      geno_word &= geno_word - 1;
                    } while (geno_word);
                  }
                }
                const double het_ctd = u31tod(genocounts[1]);
                const double homalt_ctd = u31tod(genocounts[2]);
                xtx_inv[cur_predictor_ct] = het_ctd * geno_d_lookup[0] + homalt_ctd * geno_d_lookup[1];
                xtx_inv[cur_predictor_ct + 1] = het_ctd * geno_d_lookup[0] * geno_d_lookup[0] + homalt_ctd * geno_d_lookup[1] * geno_d_lookup[1];
                if (domdev_present) {
                  xtx_inv[cur_predictor_ct + 2] = domdev_geno_prod;
                  xtx_inv[2 * cur_predictor_ct] = het_ctd;
                  xtx_inv[2 * cur_predictor_ct + 2] = het_ctd;
                }
                for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                  double* cur_xt_y = &(xt_y[cur_predictor_ct * pheno_idx]);
                  cur_xt_y[1] = geno_pheno_prods[pheno_idx];
                  if (domdev_present) {
                    cur_xt_y[2] = domdev_pheno_prods[pheno_idx];
                  }
                }
              } else {
                // !sparse_optimization
                uintptr_t start_pred_idx = 0;
                if (!(model_dominant || model_recessive || model_hetonly || joint_hethom)) {
                  start_pred_idx = domdev_present + 2;
                  xtx_inv[cur_predictor_ct] = main_dosage_sum;
                  xtx_inv[cur_predictor_ct + 1] = main_dosage_ssq;
                }
                if (cur_predictor_ct > start_pred_idx) {
                  // categorical optimization possible here
                  ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[nm_sample_ct]), &(nm_predictors_pmaj_buf[start_pred_idx * nm_sample_ct]), nm_sample_ct, nm_sample_ct, cur_predictor_ct - start_pred_idx, &(xtx_inv[cur_predictor_ct + start_pred_idx]));
                }
                RowMajorMatrixMultiplyStrided(nm_pheno_buf, &(nm_predictors_pmaj_buf[nm_sample_ct]), subbatch_size, nm_sample_ct, 1, 1, nm_sample_ct, cur_predictor_ct, &(xt_y[1]));
                if (domdev_present) {
                  RowMajorMatrixMultiplyStrided(nm_pheno_buf, &(nm_predictors_pmaj_buf[2 * nm_sample_ct]), subbatch_size, nm_sample_ct, 1, 1, nm_sample_ct, cur_predictor_ct, &(xt_y[2]));
                  // categorical optimization possible here
                  ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[2 * nm_sample_ct]), nm_predictors_pmaj_buf, nm_sample_ct, nm_sample_ct, cur_predictor_ct, &(xtx_inv[2 * cur_predictor_ct]));
                  xtx_inv[cur_predictor_ct + 2] = xtx_inv[2 * cur_predictor_ct + 1];
                }
              }
              glm_err = CheckMaxCorrAndVifNm(xtx_inv, corr_inv, cur_predictor_ct, domdev_present_p1, cur_sample_ct_recip, cur_sample_ct_m1_recip, max_corr, vif_thresh, semicomputed_biallelic_corr_matrix, semicomputed_biallelic_inv_corr_sqrts, dbl_2d_buf, &(dbl_2d_buf[2 * cur_predictor_ct]), &(dbl_2d_buf[3 * cur_predictor_ct]));
              if (glm_err) {
                goto GlmLinearSubbatchThread_skip_regression;
              }
              const double geno_ssq = xtx_inv[1 + cur_predictor_ct];
              if (!domdev_present) {
                xtx_inv[1 + cur_predictor_ct] = xtx_inv[cur_predictor_ct];
                if (InvertRank1Symm(covarx_dotprod_inv, &(xtx_inv[1 + cur_predictor_ct]), cur_predictor_ct - 1, 1, geno_ssq, dbl_2d_buf, inverse_corr_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                  goto GlmLinearSubbatchThread_skip_regression;
                }
              } else {
                const double domdev_geno_prod = xtx_inv[2 + cur_predictor_ct];
                const double domdev_ssq = xtx_inv[2 + 2 * cur_predictor_ct];
                xtx_inv[2 + cur_predictor_ct] = xtx_inv[cur_predictor_ct];
                xtx_inv[2 + 2 * cur_predictor_ct] = xtx_inv[2 * cur_predictor_ct];
                if (InvertRank2Symm(covarx_dotprod_inv, &(xtx_inv[2 + cur_predictor_ct]), cur_predictor_ct - 2, cur_predictor_ct, 1, geno_ssq, domdev_geno_prod, domdev_ssq, dbl_2d_buf, inverse_corr_buf, &(inverse_corr_buf[2 * (cur_predictor_ct - 2)]))) {
                  glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                  goto GlmLinearSubbatchThread_skip_regression;
                }
              }
              // need to make sure xtx_inv remains reflected in NOLAPACK case
              memcpy(xtx_inv, dbl_2d_buf, cur_predictor_ct * cur_predictor_ct * sizeof(double));
              ReflectMatrix(cur_predictor_ct, xtx_inv);
              RowMajorMatrixMultiply(xt_y, xtx_inv, subbatch_size, cur_predictor_ct, cur_predictor_ct, fitted_coefs);
            } else {
              // generic case
              // major categorical optimization possible here, some
              // multiallelic optimizations possible
              MultiplySelfTranspose(nm_predictors_pmaj_buf, cur_predictor_ct, nm_sample_ct, xtx_inv);

              for (uint32_t pred_idx = 1; pred_idx != cur_predictor_ct; ++pred_idx) {
                dbl_2d_buf[pred_idx] = xtx_inv[pred_idx * cur_predictor_ct];
              }
              glm_err = CheckMaxCorrAndVif(xtx_inv, 1, cur_predictor_ct, nm_sample_ct, max_corr, vif_thresh, dbl_2d_buf, nullptr, inverse_corr_buf, inv_1d_buf);
              if (glm_err) {
                goto GlmLinearSubbatchThread_skip_regression;
              }
              if (LinearRegressionInv(nm_pheno_buf, nm_predictors_pmaj_buf, cur_predictor_ct, nm_sample_ct, subbatch_size, xtx_inv, fitted_coefs, xt_y, inv_1d_buf, dbl_2d_buf)) {
                glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                goto GlmLinearSubbatchThread_skip_regression;
              }
            }
            DPrintf("\nWorker thread %" PRIuPTR " finished main regression for %u", tidx, variant_uidx);
            // RSS = y^T y - y^T X (X^T X)^{-1} X^T y
            //     = cur_pheno_ssq - xt_y * fitted_coefs
            // s^2 = RSS / df
            // possible todo: improve numerical stability of this computation
            // in non-mean-centered phenotype case

            for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
              {
                double tmp_pheno_ssq = pheno_ssq_bases[pheno_idx];
                if (missing_ct) {
                  const double* cur_pheno = &(cur_pheno_pmaj[pheno_idx * cur_sample_ct]);
                  uintptr_t sample_midx_base = 0;
                  uintptr_t sample_nm_inv_bits = ~sample_nm[0];
                  for (uint32_t missing_idx = 0; missing_idx != missing_ct; ++missing_idx) {
                    const uintptr_t sample_midx = BitIter0(sample_nm, &sample_midx_base, &sample_nm_inv_bits);
                    tmp_pheno_ssq -= cur_pheno[sample_midx] * cur_pheno[sample_midx];
                  }
                }
                const double* tmp_fitted_coefs = &(fitted_coefs[pheno_idx * cur_predictor_ct]);
                const double sigma = (tmp_pheno_ssq - DotprodxD(&(xt_y[pheno_idx * cur_predictor_ct]), tmp_fitted_coefs, cur_predictor_ct)) / u31tod(nm_sample_ct - cur_predictor_ct);
                for (uint32_t uii = 0; uii != cur_predictor_ct; ++uii) {
                  const double* s_iter = &(xtx_inv[uii * cur_predictor_ct]);
                  double* s_iter2 = &(xtx_inv2[uii * cur_predictor_ct]);
#ifdef NOLAPACK
                  for (uint32_t ujj = 0; ujj != cur_predictor_ct; ++ujj) {
                    s_iter2[ujj] = s_iter[ujj] * sigma;
                  }
#else
                  for (uint32_t ujj = 0; ujj <= uii; ++ujj) {
                    s_iter2[ujj] = s_iter[ujj] * sigma;
                  }
#endif
                }
                // validParameters() check
                for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
                  const double xtx_inv2_diag_element = xtx_inv2[pred_uidx * (cur_predictor_ct + 1)];
                  if (xtx_inv2_diag_element < 1e-20) {
                    glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                    goto GlmLinearSubbatchThread_skip_regression_for_one_pheno;
                  }
                  // use dbl_2d_buf[] to store diagonal square roots
                  dbl_2d_buf[pred_uidx] = sqrt(xtx_inv2_diag_element);
                }
                dbl_2d_buf[0] = sqrt(xtx_inv2[0]);
                for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
                  const double cur_xtx_inv2_diag_sqrt = 0.99999 * dbl_2d_buf[pred_uidx];
                  const double* xtx_inv2_row = &(xtx_inv2[pred_uidx * cur_predictor_ct]);
                  for (uint32_t pred_uidx2 = 0; pred_uidx2 != pred_uidx; ++pred_uidx2) {
                    if (xtx_inv2_row[pred_uidx2] > cur_xtx_inv2_diag_sqrt * dbl_2d_buf[pred_uidx2]) {
                      glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                      goto GlmLinearSubbatchThread_skip_regression_for_one_pheno;
                    }
                  }
                }
                double* beta_se_iter2 = beta_se_iter;
                for (uint32_t pred_uidx = reported_pred_uidx_start; pred_uidx != reported_pred_uidx_biallelic_end; ++pred_uidx) {
                  // In the multiallelic-fused case, if the first allele is
                  // constant, this writes the beta/se values for the first
                  // nonconstant, non-omitted allele where the results for the
                  // first allele belong.  We correct that below.
                  *beta_se_iter2++ = tmp_fitted_coefs[pred_uidx];
                  *beta_se_iter2++ = dbl_2d_buf[pred_uidx];
                }
                // move this up since dbl_2d_buf may be clobbered
                if (cur_constraint_ct) {
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeRankDeficient);
                  memcpy(beta_se_iter2, &glm_err2, 8);
                  ++beta_se_iter2;
                  *beta_se_iter2++ = -9.0;
                }
                if (!const_allele_ct) {
                  if (beta_se_multiallelic_fused || (!hide_covar)) {
                    for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
                      beta_se_iter2[2 * extra_allele_idx] = tmp_fitted_coefs[cur_biallelic_predictor_ct + extra_allele_idx];
                      beta_se_iter2[2 * extra_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + extra_allele_idx];
                    }
                  }
                } else if (!beta_se_multiallelic_fused) {
                  if (!hide_covar) {
                    // Need to insert some {CONST_ALLELE, -9} entries.
                    const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                    const uint32_t cur_raw_allele_idx = extra_regression_idx + (extra_regression_idx >= omitted_allele_idx);
                    uint32_t extra_read_allele_idx = 0;
                    uint32_t extra_write_allele_idx = 0;
                    for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                      if ((allele_idx == omitted_allele_idx) || (allele_idx == cur_raw_allele_idx)) {
                        continue;
                      }
                      if (IsSet(const_alleles, allele_idx)) {
                        memcpy(&(beta_se_iter2[2 * extra_write_allele_idx]), &glm_err2, 8);
                        beta_se_iter2[2 * extra_write_allele_idx + 1] = -9.0;
                      } else {
                        beta_se_iter2[2 * extra_write_allele_idx] = tmp_fitted_coefs[cur_biallelic_predictor_ct + extra_read_allele_idx];
                        beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + extra_read_allele_idx];
                        ++extra_read_allele_idx;
                      }
                      ++extra_write_allele_idx;
                    }
                  }
                } else {
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                  // Special-case first nonconst allele since it's positioned
                  // discontinuously, and its BETA/SE may already be correctly
                  // filled.
                  uint32_t allele_idx = omitted_allele_idx? 0 : 1;
                  uint32_t extra_write_allele_idx = 0;
                  if (IsSet(const_alleles, allele_idx)) {
                    memcpy(&(beta_se_iter[2 * include_intercept]), &glm_err2, 8);
                    beta_se_iter[2 * include_intercept + 1] = -9.0;
                    allele_idx = AdvTo0Bit(const_alleles, 1);
                    if (allele_idx == omitted_allele_idx) {
                      allele_idx = AdvTo0Bit(const_alleles, omitted_allele_idx + 1);
                    }
                    extra_write_allele_idx = allele_idx - 1 - (allele_idx > omitted_allele_idx);
                    for (uint32_t uii = 0; uii != extra_write_allele_idx; ++uii) {
                      memcpy(&(beta_se_iter2[2 * uii]), &glm_err2, 8);
                      beta_se_iter2[2 * uii + 1] = -9.0;
                    }
                    beta_se_iter2[2 * extra_write_allele_idx] = tmp_fitted_coefs[1];
                    beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[1];
                    ++extra_write_allele_idx;
                  }
                  ++allele_idx;
                  uint32_t nonconst_allele_idx_m1 = 0;
                  for (; allele_idx != allele_ct; ++allele_idx) {
                    if (allele_idx == omitted_allele_idx) {
                      continue;
                    }
                    if (!IsSet(const_alleles, allele_idx)) {
                      beta_se_iter2[2 * extra_write_allele_idx] = tmp_fitted_coefs[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                      beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                      ++nonconst_allele_idx_m1;
                    } else {
                      memcpy(&(beta_se_iter2[2 * extra_write_allele_idx]), &glm_err2, 8);
                      beta_se_iter2[2 * extra_write_allele_idx + 1] = -9.0;
                    }
                    ++extra_write_allele_idx;
                  }
                }
                if (cur_constraint_ct) {
                  // probable todo: cur_constraints_con_major manipulation can
                  // go outside the pheno_idx loop
                  uint32_t joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                  for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                    joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                    cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 1.0;
                  }
#ifndef NOLAPACK
                  // xtx_inv2 upper triangle was not filled
                  ReflectMatrix(cur_predictor_ct, xtx_inv2);
#endif
                  double chisq;
                  if (!LinearHypothesisChisq(tmp_fitted_coefs, cur_constraints_con_major, xtx_inv2, cur_constraint_ct, cur_predictor_ct, cur_predictor_ct, &chisq, tmphxs_buf, h_transpose_buf, inner_buf, inv_1d_buf, dbl_2d_buf)) {
                    beta_se_iter2[-1] = chisq;
                  }
                  // next test may have different alt allele count
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                  for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                    joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                    cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 0.0;
                  }
                }
              }
              while (0) {
              GlmLinearSubbatchThread_skip_regression_for_one_pheno:
                {
                  uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
                  if (beta_se_multiallelic_fused || (!hide_covar)) {
                    reported_ct += allele_ct_m2;
                  }
                  for (uint32_t uii = 0; uii != reported_ct; ++uii) {
                    memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
                    beta_se_iter[uii * 2 + 1] = -9.0;
                  }
                }
              }
              beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
            }
            while (0) {
            GlmLinearSubbatchThread_skip_regression:
              {
                uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
                if (beta_se_multiallelic_fused || (!hide_covar)) {
                  reported_ct += allele_ct_m2;
                }
                for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                  for (uint32_t uii = 0; uii != reported_ct; ++uii) {
                    memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
                    beta_se_iter[uii * 2 + 1] = -9.0;
                  }
                  beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
                }
              }
            }
          }
        }
        if (local_covars_iter) {
          local_covars_iter = &(local_covars_iter[local_covar_ct * max_sample_ct]);
        }
      }
    }
    parity = 1 - parity;
    variant_idx_offset += cur_block_variant_ct;
    while (0) {
    GlmLinearSubbatchThread_err:
      UpdateU64IfSmaller(new_err_info, &common->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

CONSTI32(kMaxLinearSubbatchSize, 240);
static_assert(kMaxLinearSubbatchSize + 12 <= kMaxOpenFiles, "kMaxLinearSubbatchSize can't be too close to or larger than kMaxOpenFiles.");

PglErr GlmLinearBatch(const uintptr_t* pheno_batch, const PhenoCol* pheno_cols, const char* pheno_names, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, uint32_t raw_variant_ct, uint32_t completed_pheno_ct, uint32_t batch_size, uintptr_t max_pheno_name_blen, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLinearCtx* ctx, TextStream* local_covar_txsp, LlStr** outfnames_ll_ptr, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char** cswritep_arr = nullptr;
  CompressStreamState* css_arr = nullptr;
  uint32_t subbatch_size = 0;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  {
    GlmCtx* common = ctx->common;
    const uintptr_t* variant_include = common->variant_include;
    const ChrInfo* cip = common->cip;
    const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
    const AlleleCode* omitted_alleles = common->omitted_alleles;

    const uint32_t sample_ct = common->sample_ct;
    const uint32_t sample_ct_x = common->sample_ct_x;
    const uint32_t sample_ct_y = common->sample_ct_y;
    const uint32_t covar_ct = common->covar_ct;
    const uintptr_t local_covar_ct = common->local_covar_ct;
    const uint32_t covar_ct_x = common->covar_ct_x;
    const uint32_t covar_ct_y = common->covar_ct_y;

    uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
    if (max_sample_ct < sample_ct_y) {
      max_sample_ct = sample_ct_y;
    }
    const uintptr_t* sample_include = common->sample_include;
    uint32_t* local_sample_idx_order = nullptr;
    uint32_t local_line_idx = 0;
    uint32_t local_xy = 0;  // 1 = chrX, 2 = chrY

    const char* local_line_iter = nullptr;
    uint32_t local_prev_chr_code = UINT32_MAX;
    uint32_t local_chr_code = UINT32_MAX;
    uint32_t local_bp = UINT32_MAX;
    uint32_t local_skip_chr = 1;
    if (local_covar_ct) {
      reterr = TextRewind(local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLinearBatch_ret_TSTREAM_FAIL;
      }
      local_line_idx = glm_info_ptr->local_header_line_ct;
      reterr = TextSkip(local_line_idx, local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLinearBatch_ret_TSTREAM_FAIL;
      }
      if (unlikely(bigstack_alloc_u32(local_sample_ct, &local_sample_idx_order))) {
        goto GlmLinearBatch_ret_NOMEM;
      }
      for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
        const uint32_t cur_uidx = local_sample_uidx_order[uii];
        uint32_t cur_idx = UINT32_MAX;
        if ((cur_uidx != UINT32_MAX) && IsSet(sample_include, cur_uidx)) {
          cur_idx = RawToSubsettedPos(sample_include, common->sample_include_cumulative_popcounts, cur_uidx);
        }
        local_sample_idx_order[uii] = cur_idx;
      }
    }

    const uint32_t variant_ct = common->variant_ct;

    const GlmFlags glm_flags = glm_info_ptr->flags;
    const uint32_t report_neglog10p = (glm_flags / kfGlmLog10) & 1;
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t domdev_present_p1 = domdev_present + 1;

    const uint32_t constraint_ct = common->constraint_ct;
    const uint32_t constraint_ct_x = common->constraint_ct_x;
    const uint32_t constraint_ct_y = common->constraint_ct_y;

    const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
    uint32_t biallelic_predictor_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
    uint32_t biallelic_predictor_ct_x = 2 + covar_ct_x * (1 + add_interactions);
    uint32_t biallelic_predictor_ct_y = 2 + covar_ct_y * (1 + add_interactions);
    const uintptr_t* parameter_subset = common->parameter_subset;
    const uintptr_t* parameter_subset_x = common->parameter_subset_x;
    const uintptr_t* parameter_subset_y = common->parameter_subset_y;
    if (parameter_subset) {
      biallelic_predictor_ct = PopcountWords(parameter_subset, BitCtToWordCt(biallelic_predictor_ct));
      if (sample_ct_x) {
        biallelic_predictor_ct_x = PopcountWords(parameter_subset_x, BitCtToWordCt(biallelic_predictor_ct_x));
      } else {
        biallelic_predictor_ct_x = 0;
      }
      if (sample_ct_y) {
        biallelic_predictor_ct_y = PopcountWords(parameter_subset_y, BitCtToWordCt(biallelic_predictor_ct_y));
      } else {
        biallelic_predictor_ct_y = 0;
      }
    }
    uint32_t biallelic_reported_test_ct = GetBiallelicReportedTestCt(parameter_subset, glm_flags, covar_ct, common->tests_flag);
    uintptr_t max_reported_test_ct = biallelic_reported_test_ct;
    uint32_t biallelic_reported_test_ct_x = 0;
    if (sample_ct_x) {
      biallelic_reported_test_ct_x = GetBiallelicReportedTestCt(parameter_subset_x, glm_flags, covar_ct_x, common->tests_flag);
      if (biallelic_reported_test_ct_x > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_x;
      }
    }
    uint32_t biallelic_reported_test_ct_y = 0;
    if (sample_ct_y) {
      biallelic_reported_test_ct_y = GetBiallelicReportedTestCt(parameter_subset_y, glm_flags, covar_ct_y, common->tests_flag);
      if (biallelic_reported_test_ct_y > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_y;
      }
    }
    const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
    const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
    const GlmColFlags glm_cols = glm_info_ptr->cols;
    const uint32_t test_col = glm_cols & kfGlmColTest;
    if (unlikely((!test_col) && (max_reported_test_ct > 1))) {
      logerrputs("Error: --glm's 'test' column cannot be omitted when results for multiple\npredictors are reported.  (Did you forget 'hide-covar'?)\n");
      goto GlmLinearBatch_ret_INCONSISTENT_INPUT;
    }
    const uint32_t main_mutated = ((glm_flags & (kfGlmDominant | kfGlmRecessive | kfGlmHetonly | kfGlmHethom)) != kfGlm0);
    const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!main_mutated) && (!common->tests_flag) && (!add_interactions);
    if (beta_se_multiallelic_fused || (!hide_covar)) {
      max_reported_test_ct += max_extra_allele_ct;
    }
    common->max_reported_test_ct = max_reported_test_ct;

    uint32_t x_code = UINT32_MAXM1;
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    if (sample_ct_x) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetX, &x_code, &x_start, &x_end);
    }
    uint32_t y_code = UINT32_MAXM1;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    if (sample_ct_y) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetY, &y_code, &y_start, &y_end);
    }
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t chr_col = glm_cols & kfGlmColChrom;

    char* chr_buf = nullptr;
    if (chr_col) {
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto GlmLinearBatch_ret_NOMEM;
      }
    }

    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > variant_ct) {
      calc_thread_ct = variant_ct;
    }

    subbatch_size = MINV(kMaxLinearSubbatchSize, batch_size);
    if (sample_ct * S_CAST(uint64_t, subbatch_size) > 0x7fffffff) {
      subbatch_size = 0x7fffffff / sample_ct;
    }
    if (unlikely(bigstack_calloc_cp(subbatch_size, &cswritep_arr) ||
                 BIGSTACK_ALLOC_X(CompressStreamState, subbatch_size, &css_arr))) {
      goto GlmLinearBatch_ret_NOMEM;
    }
    for (uint32_t fidx = 0; fidx != subbatch_size; ++fidx) {
      PreinitCstream(&(css_arr[fidx]));
    }
    // Subbatch-size-dependent memory allocations:
    //   nm_pheno_buf
    //   fitted_coefs
    //   xt_y
    //   block_beta_se_bufs
    //   css_arr

    const uint32_t main_omitted = (parameter_subset && (!IsSet(parameter_subset, 1)));
    const uint32_t xmain_ct = main_mutated + main_omitted;
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    common->thread_mhc = nullptr;
    common->dosage_presents = nullptr;
    common->dosage_mains = nullptr;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    // This cannot be less than what InitCstreamAlloc() actually allocates.
    uintptr_t cstream_alloc_size = RoundUpPow2(overflow_buf_size, kCacheline);
    if (output_zst) {
      cstream_alloc_size += RoundUpPow2(CstreamWkspaceReq(overflow_buf_size), kCacheline);
    }
    unsigned char* bigstack_mark2 = g_bigstack_base;
    double* pheno_d = nullptr;
    double* pheno_x_d = nullptr;
    double* pheno_y_d = nullptr;
    uintptr_t workspace_alloc;
    uint32_t read_block_size;
    uintptr_t max_alt_allele_block_size;
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    for (; ; subbatch_size = (subbatch_size + 1) / 2) {
      // may permit size 1 here and get rid of GlmLinear() later
      if (subbatch_size == 1) {
        goto GlmLinearBatch_ret_NOMEM;
      }
      BigstackReset(bigstack_mark2);
      if (bigstack_alloc_d(sample_ct * subbatch_size, &pheno_d)) {
        continue;
      }
      if (common->nm_precomp) {
        if (bigstack_alloc_d((2 + domdev_present + covar_ct) * subbatch_size, &(common->nm_precomp->xt_y_image))) {
          continue;
        }
      }
      if (sample_ct_x) {
        if (bigstack_alloc_d(sample_ct_x * subbatch_size, &pheno_x_d)) {
          continue;
        }
        if (common->nm_precomp_x) {
          // domdev_present can't be set here.
          if (bigstack_alloc_d((2 + covar_ct_x) * subbatch_size, &(common->nm_precomp_x->xt_y_image))) {
            continue;
          }
        }
      }
      if (sample_ct_y) {
        if (bigstack_alloc_d(sample_ct_y * subbatch_size, &pheno_y_d)) {
          continue;
        }
        if (common->nm_precomp_y) {
          if (bigstack_alloc_d((2 + covar_ct_y) * subbatch_size, &(common->nm_precomp_y->xt_y_image))) {
            continue;
          }
        }
      }
      workspace_alloc = GetLinearSubbatchWorkspaceSize(sample_ct, subbatch_size, biallelic_predictor_ct, max_extra_allele_ct, constraint_ct, xmain_ct);
      if (sample_ct_x) {
        const uintptr_t workspace_alloc_x = GetLinearSubbatchWorkspaceSize(sample_ct_x, subbatch_size, biallelic_predictor_ct_x, max_extra_allele_ct, constraint_ct_x, xmain_ct);
        if (workspace_alloc_x > workspace_alloc) {
          workspace_alloc = workspace_alloc_x;
        }
      }
      if (sample_ct_y) {
        const uintptr_t workspace_alloc_y = GetLinearSubbatchWorkspaceSize(sample_ct_y, subbatch_size, biallelic_predictor_ct_y, max_extra_allele_ct, constraint_ct_y, xmain_ct);
        if (workspace_alloc_y > workspace_alloc) {
          workspace_alloc = workspace_alloc_y;
        }
      }
      uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;
      uintptr_t per_variant_xalloc_byte_ct = max_sample_ct * local_covar_ct * sizeof(double);
      uintptr_t per_alt_allele_xalloc_byte_ct = sizeof(LinearAuxResult);
      if (beta_se_multiallelic_fused) {
        per_variant_xalloc_byte_ct += 2 * max_reported_test_ct * subbatch_size * sizeof(double);
      } else {
        per_alt_allele_xalloc_byte_ct += 2 * max_reported_test_ct * subbatch_size * sizeof(double);
      }

      uintptr_t bytes_avail = bigstack_left();
      if (bytes_avail < cstream_alloc_size * subbatch_size) {
        continue;
      }
      bytes_avail -= cstream_alloc_size * subbatch_size;
      if (!PgenMtLoadInit(variant_include, max_sample_ct, variant_ct, bytes_avail, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, per_alt_allele_xalloc_byte_ct, pgfip, &calc_thread_ct, &common->genovecs, max_extra_allele_ct? (&common->thread_mhc) : nullptr, nullptr, nullptr, dosage_is_present? (&common->dosage_presents) : nullptr, dosage_is_present? (&common->dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, &max_alt_allele_block_size, main_loadbufs, &common->pgr_ptrs, &common->read_variant_uidx_starts)) {
        break;
      }
    }
    g_failed_alloc_attempt_size = 0;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto GlmLinearBatch_ret_NOMEM;
    }
    LinearAuxResult* linear_block_aux_bufs[2];
    double* block_beta_se_bufs[2];

    for (uint32_t uii = 0; uii != 2; ++uii) {
      if (unlikely(BIGSTACK_ALLOC_X(LinearAuxResult, max_alt_allele_block_size, &(linear_block_aux_bufs[uii])))) {
        // shouldn't be possible for these to fail?
        goto GlmLinearBatch_ret_NOMEM;
      }
      if (beta_se_multiallelic_fused) {
        if (unlikely(bigstack_alloc_d(read_block_size * (2 * k1LU) * max_reported_test_ct * subbatch_size, &(block_beta_se_bufs[uii])))) {
          goto GlmLinearBatch_ret_NOMEM;
        }
      } else {
        if (unlikely(bigstack_alloc_d(max_alt_allele_block_size * 2 * max_reported_test_ct * subbatch_size, &(block_beta_se_bufs[uii])))) {
          goto GlmLinearBatch_ret_NOMEM;
        }
      }
      if (local_covar_ct) {
        if (unlikely(bigstack_alloc_d(read_block_size * max_sample_ct * local_covar_ct, &(ctx->local_covars_vcmaj_d[uii])))) {
          goto GlmLinearBatch_ret_NOMEM;
        }
      } else {
        ctx->local_covars_vcmaj_d[uii] = nullptr;
      }
    }

    common->workspace_bufs = S_CAST(unsigned char**, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t)));
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      common->workspace_bufs[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(workspace_alloc));
    }
    common->err_info = (~0LLU) << 32;
    SetThreadFuncAndData(GlmLinearSubbatchThread, ctx, &tg);

    const uint32_t subbatch_ct = 1 + (batch_size - 1) / subbatch_size;

    const uint32_t ref_col = glm_cols & kfGlmColRef;
    const uint32_t alt1_col = glm_cols & kfGlmColAlt1;
    const uint32_t alt_col = glm_cols & kfGlmColAlt;
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t all_nonref = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags, glm_cols / kfGlmColMaybeprovref, raw_variant_ct, all_nonref);
    const uint32_t omitted_col = glm_cols & kfGlmColOmitted;
    const uint32_t ax_col = glm_cols & kfGlmColAx;
    const uint32_t a1_ct_col = glm_cols & kfGlmColA1count;
    const uint32_t tot_allele_col = glm_cols & kfGlmColTotallele;
    const uint32_t a1_freq_col = glm_cols & kfGlmColA1freq;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t beta_col = glm_cols & (kfGlmColBeta | kfGlmColOrbeta);
    const uint32_t se_col = glm_cols & kfGlmColSe;
    const uint32_t ci_col = (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t t_col = glm_cols & kfGlmColTz;
    const uint32_t p_col = glm_cols & kfGlmColP;
    const uint32_t err_col = glm_cols & kfGlmColErr;
    double ci_zt = 0.0;
    if (ci_col) {
      ci_zt = QuantileToZscore((ci_size + 1.0) * 0.5);
    }
    for (uint32_t subbatch_idx = 0; subbatch_idx != subbatch_ct; ++subbatch_idx) {
      const uint32_t pheno_uidx_start = IdxToUidxBasic(pheno_batch, subbatch_idx * subbatch_size);
      if (subbatch_idx == subbatch_ct - 1) {
        subbatch_size = batch_size - subbatch_idx * subbatch_size;
      }
      ctx->subbatch_size = subbatch_size;
      uint32_t pheno_uidx = pheno_uidx_start;
      for (uint32_t fidx = 0; fidx != subbatch_size; ++fidx, ++pheno_uidx) {
        pheno_uidx = AdvTo1Bit(pheno_batch, pheno_uidx);
        const char* cur_pheno_name = &(pheno_names[pheno_uidx * max_pheno_name_blen]);
        char* outname_end2 = strcpya(&(outname_end[1]), cur_pheno_name);
        outname_end2 = strcpya_k(outname_end2, ".glm.linear");
        if (output_zst) {
          snprintf(outname_end2, 22, ".zst");
        } else {
          *outname_end2 = '\0';
        }

        // forced-singlethreaded
        char* cswritep;
        reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &(css_arr[fidx]), &cswritep);
        if (unlikely(reterr)) {
          goto GlmLinearBatch_ret_1;
        }
        if (outfnames_ll_ptr) {
          if (unlikely(PushLlStr(outname, outfnames_ll_ptr))) {
            goto GlmLinearBatch_ret_NOMEM;
          }
        }
        *cswritep++ = '#';
        if (chr_col) {
          cswritep = strcpya_k(cswritep, "CHROM\t");
        }
        if (variant_bps) {
          cswritep = strcpya_k(cswritep, "POS\t");
        }
        cswritep = strcpya_k(cswritep, "ID");
        if (ref_col) {
          cswritep = strcpya_k(cswritep, "\tREF");
        }
        if (alt1_col) {
          cswritep = strcpya_k(cswritep, "\tALT1");
        }
        if (alt_col) {
          cswritep = strcpya_k(cswritep, "\tALT");
        }
        if (provref_col) {
          cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
        }
        cswritep = strcpya_k(cswritep, "\tA1");
        if (omitted_col) {
          cswritep = strcpya_k(cswritep, "\tOMITTED");
        }
        if (ax_col) {
          cswritep = strcpya_k(cswritep, "\tAX");
        }
        if (a1_ct_col) {
          cswritep = strcpya_k(cswritep, "\tA1_CT");
        }
        if (tot_allele_col) {
          cswritep = strcpya_k(cswritep, "\tALLELE_CT");
        }
        if (a1_freq_col) {
          cswritep = strcpya_k(cswritep, "\tA1_FREQ");
        }
        if (mach_r2_col) {
          cswritep = strcpya_k(cswritep, "\tMACH_R2");
        }
        if (test_col) {
          cswritep = strcpya_k(cswritep, "\tTEST");
        }
        if (nobs_col) {
          cswritep = strcpya_k(cswritep, "\tOBS_CT");
        }
        if (beta_col) {
          cswritep = strcpya_k(cswritep, "\tBETA");
        }
        if (se_col) {
          cswritep = strcpya_k(cswritep, "\tSE");
        }
        if (ci_col) {
          cswritep = strcpya_k(cswritep, "\tL");
          cswritep = dtoa_g(ci_size * 100, cswritep);
          cswritep = strcpya_k(cswritep, "\tU");
          cswritep = dtoa_g(ci_size * 100, cswritep);
        }
        if (t_col) {
          if (!constraint_ct) {
            cswritep = strcpya_k(cswritep, "\tT_STAT");
          } else {
            // F-statistic for joint tests.
            cswritep = strcpya_k(cswritep, "\tT_OR_F_STAT");
          }
        }
        if (p_col) {
          if (report_neglog10p) {
            // TODO: change to NEG_LOG10_P for a5
            cswritep = strcpya_k(cswritep, "\tLOG10_P");
          } else {
            cswritep = strcpya_k(cswritep, "\tP");
          }
        }
        if (err_col) {
          cswritep = strcpya_k(cswritep, "\tERRCODE");
        }
        AppendBinaryEoln(&cswritep);
        cswritep_arr[fidx] = cswritep;

        FillPhenoAndXtY(sample_include, pheno_cols[pheno_uidx].data.qt, ctx->covars_cmaj_d, sample_ct, domdev_present_p1, covar_ct, common->nm_precomp? (&(common->nm_precomp->xt_y_image[fidx * (1 + domdev_present_p1 + covar_ct)])) : nullptr, &(pheno_d[fidx * sample_ct]));
        if (sample_ct_x) {
          FillPhenoAndXtY(common->sample_include_x, pheno_cols[pheno_uidx].data.qt, ctx->covars_cmaj_x_d, sample_ct_x, domdev_present_p1, covar_ct_x, common->nm_precomp_x? (&(common->nm_precomp_x->xt_y_image[fidx * (1 + domdev_present_p1 + covar_ct_x)])) : nullptr, &(pheno_x_d[fidx * sample_ct_x]));
        }
        if (sample_ct_y) {
          FillPhenoAndXtY(common->sample_include_y, pheno_cols[pheno_uidx].data.qt, ctx->covars_cmaj_y_d, sample_ct_y, domdev_present_p1, covar_ct_y, common->nm_precomp_y? (&(common->nm_precomp_y->xt_y_image[fidx * (1 + domdev_present_p1 + covar_ct_y)])) : nullptr, &(pheno_y_d[fidx * sample_ct_y]));
        }
        ctx->pheno_d = pheno_d;
        ctx->pheno_x_d = pheno_x_d;
        ctx->pheno_y_d = pheno_y_d;
      }

      // Main workflow:
      // 1. Set n=0, load/skip block 0
      //
      // 2. Spawn threads processing block n
      // 3. If n>0, write results for block (n-1)
      // 4. Increment n by 1
      // 5. Load/skip block n unless eof
      // 6. Join threads
      // 7. Goto step 2 unless eof
      //
      // 8, Write results for last block
      uintptr_t write_variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      uint32_t parity = 0;
      uint32_t read_block_idx = 0;
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t suppress_mach_r2 = 0;

      uint32_t cur_biallelic_reported_test_ct = 0;
      uint32_t primary_reported_test_idx = include_intercept;
      uint32_t cur_biallelic_predictor_ct = 0;
      uint32_t cur_constraint_ct = 0;

      const char* const* cur_test_names = nullptr;
      uint32_t prev_block_variant_ct = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uint32_t allele_ct = 2;
      uint32_t omitted_allele_idx = 0;
      if (subbatch_size > 1) {
        logprintfww5("--glm linear regression on quantitative phenotypes #%u-%u: ", completed_pheno_ct + 1, completed_pheno_ct + subbatch_size);
      } else {
        logprintfww5("--glm linear regression on phenotype '%s': ", &(pheno_names[pheno_uidx_start * max_pheno_name_blen]));
      }
      fputs("0%", stdout);
      fflush(stdout);
      // bugfix (12 May 2019): need to reinitialize this
      pgfip->block_base = main_loadbufs[0];
      ReinitThreads(&tg);
      for (uint32_t variant_idx = 0; ; ) {
        const uint32_t cur_block_variant_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
        if (unlikely(reterr)) {
          goto GlmLinearBatch_ret_PGR_FAIL;
        }
        if (local_covar_ct && cur_block_variant_ct) {
          const uint32_t uidx_start = read_block_idx * read_block_size;
          const uint32_t uidx_end = MINV(raw_variant_ct, uidx_start + read_block_size);
          if (local_variant_include) {
            reterr = ReadLocalCovarBlock(common, local_sample_uidx_order, local_variant_include, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, local_covar_txsp, &local_line_idx, &local_xy, nullptr, ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
          } else {
            double* prev_local_covar_row_d = nullptr;
            if (variant_idx) {
              prev_local_covar_row_d = &(ctx->local_covars_vcmaj_d[1 - parity][S_CAST(uintptr_t, read_block_size - 1) * max_sample_ct * local_covar_ct]);
            }
            reterr = ReadRfmix2Block(common, variant_bps, local_sample_uidx_order, nullptr, prev_local_covar_row_d, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, glm_info_ptr->local_chrom_col, glm_info_ptr->local_bp_col, glm_info_ptr->local_first_covar_col, local_covar_txsp, &local_line_iter, &local_line_idx, &local_prev_chr_code, &local_chr_code, &local_bp, &local_skip_chr, nullptr, ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
          }
          if (unlikely(reterr)) {
            goto GlmLinearBatch_ret_1;
          }
        }
        if (variant_idx) {
          DPrintf("\nWaiting for worker threads up to filtered variant %u", variant_idx);
          JoinThreads(&tg);
          reterr = S_CAST(PglErr, common->err_info);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, common->err_info >> 32);
            goto GlmLinearBatch_ret_1;
          }
        }
        if (!IsLastBlock(&tg)) {
          common->cur_block_variant_ct = cur_block_variant_ct;
          const uint32_t uidx_start = read_block_idx * read_block_size;
          ComputeUidxStartPartition(variant_include, cur_block_variant_ct, calc_thread_ct, uidx_start, common->read_variant_uidx_starts);
          PgrCopyBaseAndOffset(pgfip, calc_thread_ct, common->pgr_ptrs);
          ctx->block_aux = linear_block_aux_bufs[parity];
          common->block_beta_se = block_beta_se_bufs[parity];
          if (variant_idx + cur_block_variant_ct == variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          DPrintf("\nLaunching worker threads for filtered variants [%u,%u)", variant_idx, variant_idx + cur_block_variant_ct);
          if (unlikely(SpawnThreads(&tg))) {
            goto GlmLinearBatch_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (variant_idx) {
          // write *previous* block results
          const double* beta_se_iter = block_beta_se_bufs[parity];
          const LinearAuxResult* cur_block_aux = linear_block_aux_bufs[parity];
          uintptr_t allele_bidx = 0;
          for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx) {
            const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
            if (write_variant_uidx >= chr_end) {
              do {
                ++chr_fo_idx;
                chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
              } while (write_variant_uidx >= chr_end);
              const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
              if ((chr_idx == x_code) && sample_ct_x) {
                cur_biallelic_reported_test_ct = biallelic_reported_test_ct_x;
                cur_biallelic_predictor_ct = biallelic_predictor_ct_x;
                cur_constraint_ct = constraint_ct_x;
                cur_test_names = test_names_x;
              } else if ((chr_idx == y_code) && sample_ct_y) {
                cur_biallelic_reported_test_ct = biallelic_reported_test_ct_y;
                cur_biallelic_predictor_ct = biallelic_predictor_ct_y;
                cur_constraint_ct = constraint_ct_y;
                cur_test_names = test_names_y;
              } else {
                cur_biallelic_reported_test_ct = biallelic_reported_test_ct;
                cur_biallelic_predictor_ct = biallelic_predictor_ct;
                cur_constraint_ct = constraint_ct;
                cur_test_names = test_names;
              }
              suppress_mach_r2 = (chr_idx == x_code) || (chr_idx == mt_code);
              if (cur_constraint_ct) {
                primary_reported_test_idx = cur_biallelic_reported_test_ct - 1;
              }
              if (chr_col) {
                char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
                *chr_name_end = '\t';
                chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
              }
            }
            uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
            if (allele_idx_offsets) {
              allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
              allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offsets[write_variant_uidx];
            }
            const uint32_t allele_ct_m1 = allele_ct - 1;
            const uint32_t extra_allele_ct = allele_ct - 2;
            if (omitted_alleles) {
              omitted_allele_idx = omitted_alleles[write_variant_uidx];
            }
            const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
            for (uint32_t fidx = 0; fidx != subbatch_size; ++fidx) {
              char* cswritep = cswritep_arr[fidx];
              uint32_t a1_allele_idx = 0;
              uintptr_t allele_bidx_tmp = allele_bidx;
              // bugfix (20 Mar 2020): In multiallelic unfused case, fidx is
              // the *inner* index, allele index is on the outside.
              const double* beta_se_iter2 = beta_se_iter;
              for (uint32_t nonomitted_allele_idx = 0; nonomitted_allele_idx != allele_ct_m1; ++nonomitted_allele_idx, ++a1_allele_idx) {
                if (beta_se_multiallelic_fused) {
                  if (!nonomitted_allele_idx) {
                    primary_reported_test_idx = include_intercept;
                  } else {
                    primary_reported_test_idx = cur_biallelic_reported_test_ct + nonomitted_allele_idx - 1;
                  }
                }
                if (nonomitted_allele_idx == omitted_allele_idx) {
                  ++a1_allele_idx;
                }
                const double primary_beta = beta_se_iter2[primary_reported_test_idx * 2];
                const double primary_se = beta_se_iter2[primary_reported_test_idx * 2 + 1];
                const uint32_t allele_is_valid = (primary_se != -9.0);
                {
                  const LinearAuxResult* auxp = &(cur_block_aux[allele_bidx_tmp]);
                  if (ln_pfilter <= 0.0) {
                    if (!allele_is_valid) {
                      goto GlmLinearBatch_allele_iterate;
                    }
                    double primary_ln_pval;
                    if (!cur_constraint_ct) {
                      if (primary_beta == 0.0) {
                        primary_ln_pval = 0.0;
                      } else if (primary_se == 0.0) {
                        primary_ln_pval = -DBL_MAX;
                      } else {
                        const double primary_tstat = primary_beta / primary_se;
                        primary_ln_pval = TstatToLnP(primary_tstat, auxp->sample_obs_ct - cur_biallelic_predictor_ct - extra_allele_ct);
                      }
                    } else {
                      primary_ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                    }
                    if (primary_ln_pval > ln_pfilter) {
                      goto GlmLinearBatch_allele_iterate;
                    }
                  }
                  uint32_t inner_reported_test_ct = cur_biallelic_reported_test_ct;
                  if (extra_allele_ct) {
                    if (beta_se_multiallelic_fused) {
                      if (!nonomitted_allele_idx) {
                        inner_reported_test_ct = 1 + include_intercept;
                      } else if (nonomitted_allele_idx == extra_allele_ct) {
                        inner_reported_test_ct -= include_intercept;
                      } else {
                        inner_reported_test_ct = 1;
                      }
                    } else if (!hide_covar) {
                      inner_reported_test_ct += extra_allele_ct;
                    }
                  }
                  // possible todo: make number-to-string operations, strlen(),
                  // etc. happen only once per variant.
                  for (uint32_t allele_test_idx = 0; allele_test_idx != inner_reported_test_ct; ++allele_test_idx) {
                    uint32_t test_idx = allele_test_idx;
                    if (beta_se_multiallelic_fused && nonomitted_allele_idx) {
                      if (!allele_test_idx) {
                        test_idx = primary_reported_test_idx;
                      } else {
                        ++test_idx;
                      }
                    }
                    if (chr_col) {
                      cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
                    }
                    if (variant_bps) {
                      cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
                    }
                    cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
                    if (ref_col) {
                      *cswritep++ = '\t';
                      cswritep = strcpya(cswritep, cur_alleles[0]);
                    }
                    if (alt1_col) {
                      *cswritep++ = '\t';
                      cswritep = strcpya(cswritep, cur_alleles[1]);
                    }
                    if (alt_col) {
                      *cswritep++ = '\t';
                      for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
                        if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                          // might not need this assignment, but play it safe
                          // for now
                          cswritep_arr[fidx] = cswritep;
                          goto GlmLinearBatch_ret_WRITE_FAIL;
                        }
                        cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                      }
                      --cswritep;
                    }
                    *cswritep++ = '\t';
                    if (provref_col) {
                      *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, write_variant_uidx)))? 'Y' : 'N';
                      *cswritep++ = '\t';
                    }
                    const uint32_t multi_a1 = extra_allele_ct && beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx);
                    if (multi_a1) {
                      for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                        if (allele_idx == omitted_allele_idx) {
                          continue;
                        }
                        if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                          cswritep_arr[fidx] = cswritep;
                          goto GlmLinearBatch_ret_WRITE_FAIL;
                        }
                        cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                      }
                      --cswritep;
                    } else {
                      cswritep = strcpya(cswritep, cur_alleles[a1_allele_idx]);
                    }
                    if (omitted_col) {
                      *cswritep++ = '\t';
                      cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
                    }
                    if (ax_col) {
                      *cswritep++ = '\t';
                      if (beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx)) {
                        if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                          cswritep_arr[fidx] = cswritep;
                          goto GlmLinearBatch_ret_WRITE_FAIL;
                        }
                        cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
                      } else {
                        for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                          if (allele_idx == a1_allele_idx) {
                            continue;
                          }
                          if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                            cswritep_arr[fidx] = cswritep;
                            goto GlmLinearBatch_ret_WRITE_FAIL;
                          }
                          cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                        }
                        --cswritep;
                      }
                    }
                    if (a1_ct_col) {
                      *cswritep++ = '\t';
                      if (!multi_a1) {
                        cswritep = dtoa_g(auxp->a1_dosage, cswritep);
                      } else {
                        cswritep = strcpya_k(cswritep, "NA");
                      }
                    }
                    if (tot_allele_col) {
                      *cswritep++ = '\t';
                      cswritep = u32toa(auxp->allele_obs_ct, cswritep);
                    }
                    if (a1_freq_col) {
                      *cswritep++ = '\t';
                      if (!multi_a1) {
                        cswritep = dtoa_g(auxp->a1_dosage / S_CAST(double, auxp->allele_obs_ct), cswritep);
                      } else {
                        cswritep = strcpya_k(cswritep, "NA");
                      }
                    }
                    if (mach_r2_col) {
                      *cswritep++ = '\t';
                      if (!suppress_mach_r2) {
                        cswritep = dtoa_g(auxp->mach_r2, cswritep);
                      } else {
                        cswritep = strcpya_k(cswritep, "NA");
                      }
                    }
                    if (test_col) {
                      *cswritep++ = '\t';
                      if (test_idx < cur_biallelic_reported_test_ct) {
                        cswritep = strcpya(cswritep, cur_test_names[test_idx]);
                      } else {
                        // always use basic dosage for untested alleles
                        cswritep = strcpya_k(cswritep, "ADD");
                        if (!beta_se_multiallelic_fused) {
                          // extra alt allele covariate.
                          uint32_t test_xallele_idx = test_idx - cur_biallelic_reported_test_ct;
                          if (omitted_allele_idx < a1_allele_idx) {
                            test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                          }
                          test_xallele_idx = test_xallele_idx + (test_xallele_idx >= a1_allele_idx);
                          if (a1_allele_idx < omitted_allele_idx) {
                            test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                          }
                          if (!test_xallele_idx) {
                            cswritep = strcpya_k(cswritep, "_REF");
                          } else {
                            cswritep = strcpya_k(cswritep, "_ALT");
                            cswritep = u32toa(test_xallele_idx, cswritep);
                          }
                        }
                      }
                    }
                    if (nobs_col) {
                      *cswritep++ = '\t';
                      cswritep = u32toa(auxp->sample_obs_ct, cswritep);
                    }
                    double ln_pval = kLnPvalError;
                    double tstat = 0.0;
                    uint32_t test_is_valid;
                    if ((!cur_constraint_ct) || (test_idx != primary_reported_test_idx)) {
                      double beta = beta_se_iter2[2 * test_idx];
                      double se = beta_se_iter2[2 * test_idx + 1];
                      test_is_valid = (se != -9.0);
                      if (test_is_valid) {
                        tstat = beta / se;
                        ln_pval = TstatToLnP(tstat, auxp->sample_obs_ct - cur_biallelic_predictor_ct - extra_allele_ct);
                      }
                      if (beta_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          cswritep = dtoa_g(beta, cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA");
                        }
                      }
                      if (se_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          cswritep = dtoa_g(se, cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA");
                        }
                      }
                      if (ci_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          const double ci_halfwidth = ci_zt * se;
                          cswritep = dtoa_g(beta - ci_halfwidth, cswritep);
                          *cswritep++ = '\t';
                          cswritep = dtoa_g(beta + ci_halfwidth, cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA\tNA");
                        }
                      }
                      if (t_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          cswritep = dtoa_g(tstat, cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA");
                        }
                      }
                    } else {
                      // joint test
                      test_is_valid = allele_is_valid;
                      if (beta_col) {
                        cswritep = strcpya_k(cswritep, "\tNA");
                      }
                      if (se_col) {
                        cswritep = strcpya_k(cswritep, "\tNA");
                      }
                      if (ci_col) {
                        cswritep = strcpya_k(cswritep, "\tNA\tNA");
                      }
                      if (t_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          cswritep = dtoa_g(primary_se / u31tod(cur_constraint_ct), cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA");
                        }
                      }
                      // could avoid recomputing
                      if (test_is_valid) {
                        ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                      }
                    }
                    if (p_col) {
                      *cswritep++ = '\t';
                      if (test_is_valid) {
                        if (report_neglog10p) {
                          const double reported_val = (-kRecipLn10) * ln_pval;
                          cswritep = dtoa_g(reported_val, cswritep);
                        } else {
                          const double reported_ln = MAXV(ln_pval, output_min_ln);
                          cswritep = lntoa_g(reported_ln, cswritep);
                        }
                      } else {
                        cswritep = strcpya_k(cswritep, "NA");
                      }
                    }
                    if (err_col) {
                      *cswritep++ = '\t';
                      if (test_is_valid) {
                        *cswritep++ = '.';
                      } else {
                        uint64_t glm_errcode;
                        memcpy(&glm_errcode, &(beta_se_iter2[2 * test_idx]), 8);
                        cswritep = AppendGlmErrstr(glm_errcode, cswritep);
                      }
                    }
                    AppendBinaryEoln(&cswritep);
                    if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                      cswritep_arr[fidx] = cswritep;
                      goto GlmLinearBatch_ret_WRITE_FAIL;
                    }
                  }
                }
              GlmLinearBatch_allele_iterate:
                ++allele_bidx_tmp;
                if (!beta_se_multiallelic_fused) {
                  beta_se_iter2 = &(beta_se_iter2[subbatch_size * (2 * k1LU) * max_reported_test_ct]);
                }
              }  // for nonomitted_allele_idx
              beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
              cswritep_arr[fidx] = cswritep;
            }  // for fidx
            if (!beta_se_multiallelic_fused) {
              beta_se_iter = &(beta_se_iter[extra_allele_ct * subbatch_size * (2 * k1LU) * max_reported_test_ct]);
            }
            allele_bidx += allele_ct_m1;
          }  // for variant_bidx
        }
        if (variant_idx == variant_ct) {
          break;
        }
        if (variant_idx >= next_print_variant_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (variant_idx * 100LLU) / variant_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
        }
        ++read_block_idx;
        prev_block_variant_ct = cur_block_variant_ct;
        variant_idx += cur_block_variant_ct;
        // crucially, this is independent of the PgenReader block_base
        // pointers
        pgfip->block_base = main_loadbufs[parity];
      }
      for (uintptr_t fidx = 0; fidx != subbatch_size; ++fidx) {
        if (unlikely(CswriteCloseNull(&(css_arr[fidx]), cswritep_arr[fidx]))) {
          goto GlmLinearBatch_ret_WRITE_FAIL;
        }
      }
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      logputs("done.\n");
      // bugfix (12 May 2019): added batch_size instead of subbatch_size here
      completed_pheno_ct += subbatch_size;
    }
    outname_end[1] = '\0';
    logprintfww("Results written to %s<phenotype name>.glm.linear%s .\n", outname, output_zst? ".zst" : "");
  }
  while (0) {
  GlmLinearBatch_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GlmLinearBatch_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
    break;
  GlmLinearBatch_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  GlmLinearBatch_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GlmLinearBatch_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  GlmLinearBatch_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 GlmLinearBatch_ret_1:
  CleanupThreads(&tg);
  if (css_arr) {
    for (uintptr_t fidx = 0; fidx != subbatch_size; ++fidx) {
      CswriteCloseCond(&(css_arr[fidx]), cswritep_arr[fidx]);
    }
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
