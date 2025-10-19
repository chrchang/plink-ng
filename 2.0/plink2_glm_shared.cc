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

#include "plink2_glm_shared.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "include/plink2_bits.h"
#include "include/plink2_string.h"
#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static const char kGlmErrcodeStrs[][24] = {"", "SAMPLE_CT<=PREDICTOR_CT", "CONST_OMITTED_ALLELE", "CONST_ALLELE", "CORR_TOO_HIGH", "VIF_INFINITE", "VIF_TOO_HIGH", "SEPARATION", "RANK_DEFICIENT", "LOGISTIC_CONVERGE_FAIL", "FIRTH_CONVERGE_FAIL", "INVALID_RESULT"};

char* AppendGlmErrstr(GlmErr glm_err, char* write_iter) {
  // todo: support predictor args
  GlmErrcode errcode = GetGlmErrCode(glm_err);
  write_iter = strcpya(write_iter, kGlmErrcodeStrs[errcode]);
  if (errcode == kGlmErrcodeSeparation) {
    *write_iter++ = ',';
    const uint32_t allele_idx = GetGlmErrArg1(glm_err);
    if (!allele_idx) {
      write_iter = strcpya_k(write_iter, "REF");
    } else {
      write_iter = strcpya_k(write_iter, "ALT");
      write_iter = u32toa(allele_idx, write_iter);
    }
  }
  return write_iter;
}

// * first_predictor_idx should be 1 if first term is constant-1 intercept, 0
//   otherwise
// * dbl_2d_buf[n] expected to be sum of predictor row n on input, is destroyed
// * if corr_buf is not nullptr, lower left of (relevant_predictor_ct x
//   relevant_predictor_ct) body is filled with uninverted correlation matrix
//   on exit, then the next relevant_predictor_ct elements are inverse-stdevs
//   of the relevant predictors
// * lower left of inverse_corr_buf is filled on return
GlmErr CheckMaxCorrAndVif(const double* predictor_dotprods, uint32_t first_predictor_idx, uint32_t predictor_ct, uintptr_t sample_ct, double max_corr, double vif_thresh, double* dbl_2d_buf, double* corr_buf, double* inverse_corr_buf, MatrixInvertBuf1* matrix_invert_buf1) {
  // we have dot products, now determine
  //   (dotprod - sum(a)mean(b)) / (N-1)
  // to get small-sample covariance
  const uintptr_t relevant_predictor_ct = predictor_ct - first_predictor_idx;
  // bugfix (19 Oct 2025): relevant_predictor_ct == 1 fast-path incorrectly
  // initialized corr_buf[1] to 1 instead of the predictor's inverse-stdev.
  // Fast-path has been deleted since exactly-1-covariate is very rare (which
  // is why the bug wasn't caught for >6 years...) and the regular-branch logic
  // doesn't break.
  const uintptr_t relevant_predictor_ct_p1 = relevant_predictor_ct + 1;
  const double sample_ct_recip = 1.0 / u31tod(sample_ct);
  const double sample_ct_m1_d = u31tod(sample_ct - 1);
  const double sample_ct_m1_recip = 1.0 / sample_ct_m1_d;
  for (uintptr_t pred_idx1 = 0; pred_idx1 != relevant_predictor_ct; ++pred_idx1) {
    double* sample_cov_row = &(inverse_corr_buf[pred_idx1 * relevant_predictor_ct]);
    const uintptr_t input_pred_idx1 = pred_idx1 + first_predictor_idx;
    const double* predictor_dotprods_row = &(predictor_dotprods[input_pred_idx1 * predictor_ct]);
    const double pred1_mean_adj = dbl_2d_buf[input_pred_idx1] * sample_ct_recip;
    for (uintptr_t pred_idx2 = 0; pred_idx2 <= pred_idx1; ++pred_idx2) {
      const uintptr_t input_pred_idx2 = pred_idx2 + first_predictor_idx;
      sample_cov_row[pred_idx2] = (predictor_dotprods_row[input_pred_idx2] - pred1_mean_adj * dbl_2d_buf[input_pred_idx2]) * sample_ct_m1_recip;
    }
  }
  // now use dbl_2d_buf to store inverse-sqrts, to get to correlation matrix
  for (uintptr_t pred_idx = 0; pred_idx != relevant_predictor_ct; ++pred_idx) {
    dbl_2d_buf[pred_idx] = 1.0 / sqrt(inverse_corr_buf[pred_idx * relevant_predictor_ct_p1]);
  }
  // invert_symmdef_matrix only cares about bottom left of inverse_corr_buf[]
  for (uintptr_t pred_idx1 = 1; pred_idx1 != relevant_predictor_ct; ++pred_idx1) {
    const double inverse_stdev1 = dbl_2d_buf[pred_idx1];
    // convert from covariances to correlations here
    double* corr_row_iter = &(inverse_corr_buf[pred_idx1 * relevant_predictor_ct]);
    const double* inverse_stdev2_iter = dbl_2d_buf;
    for (uintptr_t pred_idx2 = 0; pred_idx2 != pred_idx1; ++pred_idx2) {
      const double cur_corr = (*corr_row_iter) * inverse_stdev1 * (*inverse_stdev2_iter++);
      // bugfix (14 Sep 2017): need to take absolute value here
      if (fabs(cur_corr) > max_corr) {
        return SetGlmErr2(kGlmErrcodeCorrTooHigh, pred_idx2, pred_idx1);
      }
      *corr_row_iter++ = cur_corr;
    }
  }
  for (uintptr_t pred_idx = 0; pred_idx != relevant_predictor_ct; ++pred_idx) {
    inverse_corr_buf[pred_idx * relevant_predictor_ct_p1] = 1.0;
  }
  if (corr_buf) {
    memcpy(corr_buf, inverse_corr_buf, relevant_predictor_ct * relevant_predictor_ct * sizeof(double));
    memcpy(&(corr_buf[relevant_predictor_ct * relevant_predictor_ct]), dbl_2d_buf, relevant_predictor_ct * sizeof(double));
  }
  if (InvertSymmdefMatrixChecked(relevant_predictor_ct, inverse_corr_buf, matrix_invert_buf1, dbl_2d_buf)) {
    return SetGlmErr0(kGlmErrcodeVifInfinite);
  }
  // VIFs = diagonal elements of inverse correlation matrix
  for (uintptr_t pred_idx = 0; pred_idx != relevant_predictor_ct; ++pred_idx) {
    if (inverse_corr_buf[pred_idx * relevant_predictor_ct_p1] > vif_thresh) {
      return SetGlmErr1(kGlmErrcodeVifTooHigh, pred_idx);
    }
  }
  return 0;
}

// no-missing-genotype optimizations:
// * most of the inter-predictor correlation matrix can be initialized once
//   from an image, doesn't even need to be refreshed; same goes for
//   inv_corr_sqrts
// * inverse of that matrix can be precomputed, and used in rank 1 inverse
//   update (rank 2 for domdev case)
//
// Other notes:
// * dbl_2d_buf does not need to store row sums, since we can get that from
//   predictor_dotprods (first actual, though "irrelevant", predictor is
//   all-1).
// * predictor_ct includes intercept.
// * geno_pred_ct must be 1 or 2.
// * ainv_b_buf[] must have size at least 4 * nongeno_pred_ct.
GlmErr CheckMaxCorrAndVifNm(const double* predictor_dotprods, const double* corr_inv, uint32_t predictor_ct, uint32_t geno_pred_ct, double sample_ct_recip, double sample_ct_m1_recip, double max_corr, double vif_thresh, double* __restrict semicomputed_corr_matrix, double* __restrict semicomputed_inv_corr_sqrts, double* __restrict corr_row_buf, double* __restrict inverse_corr_diag, double* __restrict ainv_b_buf) {
  // we have dot products, now determine
  //   (dotprod - sum(a)mean(b)) / (N-1)
  // to get small-sample covariance
  if (predictor_ct == 2) {
    return 0;
  }
  const uintptr_t relevant_predictor_ct = predictor_ct - 1;
  const uintptr_t relevant_predictor_ct_p1 = predictor_ct;
  // predictor_dotprods[] *rows* 1 (and 2, if geno_pred_ct == 2) are filled,
  // rather than columns
  {
    const double pred1_mean_adj = predictor_dotprods[predictor_ct] * sample_ct_recip;
    semicomputed_corr_matrix[0] = (predictor_dotprods[predictor_ct + 1] - pred1_mean_adj * predictor_dotprods[predictor_ct]) * sample_ct_m1_recip;
  }
  for (uintptr_t pred_idx1 = 1; pred_idx1 != relevant_predictor_ct; ++pred_idx1) {
    double* sample_cov_row = &(semicomputed_corr_matrix[pred_idx1 * relevant_predictor_ct]);
    const uintptr_t input_pred_idx1 = pred_idx1 + 1;
    const double* predictor_dotprods_row = &(predictor_dotprods[input_pred_idx1 * predictor_ct]);
    const double pred1_mean_adj = predictor_dotprods_row[0] * sample_ct_recip;
    uintptr_t pred_idx2 = 0;
    for (; pred_idx2 != geno_pred_ct; ++pred_idx2) {
      const uintptr_t input_pred_idx2 = pred_idx2 + 1;
      sample_cov_row[pred_idx2] = (predictor_dotprods[input_pred_idx2 * predictor_ct + input_pred_idx1] - pred1_mean_adj * predictor_dotprods[input_pred_idx2 * predictor_ct]) * sample_ct_m1_recip;
    }
    // document whether pred_idx2 guaranteed to be <= input_pred_idx1 if this
    // code is revisited
    for (; pred_idx2 < input_pred_idx1; ++pred_idx2) {
      const uintptr_t input_pred_idx2 = pred_idx2 + 1;
      sample_cov_row[pred_idx2] = (predictor_dotprods_row[input_pred_idx2] - pred1_mean_adj * predictor_dotprods[input_pred_idx2 * predictor_ct]) * sample_ct_m1_recip;
    }
  }
  // assumes semicomputed_inv_corr_sqrts[geno_pred_ct..] is pre-initialized
  semicomputed_inv_corr_sqrts[0] = 1.0 / sqrt(semicomputed_corr_matrix[0]);
  const uintptr_t nongeno_pred_ct = relevant_predictor_ct - geno_pred_ct;
  double inverse_stdev1 = semicomputed_inv_corr_sqrts[0];
  const double* nongeno_inverse_stdevs = &(semicomputed_inv_corr_sqrts[geno_pred_ct]);
  const double* corr_col = &(semicomputed_corr_matrix[geno_pred_ct * relevant_predictor_ct]);
  for (uintptr_t nongeno_pred_idx = 0; nongeno_pred_idx != nongeno_pred_ct; ++nongeno_pred_idx) {
    const double inverse_stdev2 = nongeno_inverse_stdevs[nongeno_pred_idx];
    const double cur_corr = inverse_stdev1 * inverse_stdev2 * corr_col[nongeno_pred_idx * relevant_predictor_ct];
    if (fabs(cur_corr) > max_corr) {
      return SetGlmErr2(kGlmErrcodeCorrTooHigh, 0, geno_pred_ct + nongeno_pred_idx);
    }
    corr_row_buf[nongeno_pred_idx] = cur_corr;
  }
  if (geno_pred_ct == 1) {
    if (InvertRank1SymmDiag(corr_inv, corr_row_buf, nongeno_pred_ct, 1.0, inverse_corr_diag, ainv_b_buf)) {
      return SetGlmErr0(kGlmErrcodeVifInfinite);
    }
  } else {
    inverse_stdev1 = 1.0 / sqrt(semicomputed_corr_matrix[relevant_predictor_ct_p1]);
    corr_col = &(semicomputed_corr_matrix[geno_pred_ct * relevant_predictor_ct + 1]);
    double* corr_row2 = &(corr_row_buf[nongeno_pred_ct]);
    for (uintptr_t nongeno_pred_idx = 0; nongeno_pred_idx != nongeno_pred_ct; ++nongeno_pred_idx) {
      const double inverse_stdev2 = nongeno_inverse_stdevs[nongeno_pred_idx];
      const double cur_corr = inverse_stdev1 * inverse_stdev2 * corr_col[nongeno_pred_idx * relevant_predictor_ct];
      if (fabs(cur_corr) > max_corr) {
        return SetGlmErr2(kGlmErrcodeCorrTooHigh, 1, 2 + nongeno_pred_idx);
      }
      corr_row2[nongeno_pred_idx] = cur_corr;
    }
    const double inverse_stdev2 = semicomputed_inv_corr_sqrts[0];
    const double cur_corr = inverse_stdev1 * inverse_stdev2 * semicomputed_corr_matrix[relevant_predictor_ct];
    if (fabs(cur_corr) > max_corr) {
      return SetGlmErr2(kGlmErrcodeCorrTooHigh, 0, 1);
    }
    // do we want special handling of nongeno_pred_ct == 0?
    if (InvertRank2SymmDiag(corr_inv, corr_row_buf, nongeno_pred_ct, 1.0, cur_corr, 1.0, inverse_corr_diag, ainv_b_buf, &(ainv_b_buf[2 * nongeno_pred_ct]))) {
      return SetGlmErr0(kGlmErrcodeVifInfinite);
    }
  }
  // VIFs = diagonal elements of inverse correlation matrix
  for (uintptr_t pred_idx = 0; pred_idx != relevant_predictor_ct; ++pred_idx) {
    if (inverse_corr_diag[pred_idx] > vif_thresh) {
      return SetGlmErr1(kGlmErrcodeVifTooHigh, pred_idx);
    }
  }
  return 0;
}

PglErr GlmFillAndTestCovars(const uintptr_t* sample_include, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, double* covar_dotprod, double* corr_buf, double* inverse_corr_buf, double* covars_cmaj, const char** cur_covar_names, GlmErr* glm_err_ptr) {
  *glm_err_ptr = 0;
  if (covar_ct == local_covar_ct) {
    // bugfix (5 Mar 2018): need to copy local-covar names
    for (uintptr_t local_covar_read_idx = 0; local_covar_read_idx != covar_ct; ++local_covar_read_idx) {
      cur_covar_names[local_covar_read_idx] = &(covar_names[local_covar_read_idx * max_covar_name_blen]);
    }
    return kPglRetSuccess;
  }
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  const uint32_t covar_max_cat_ctl = 1 + (covar_max_nonnull_cat_ct / kBitsPerWord);
  MatrixInvertBuf1* matrix_invert_buf1;
  uintptr_t* cat_covar_wkspace;
  double* dbl_2d_buf;
  if (unlikely(BIGSTACK_ALLOC_X(MatrixInvertBuf1, kMatrixInvertBuf1CheckedAlloc * new_nonlocal_covar_ct, &matrix_invert_buf1) ||
               bigstack_alloc_w(covar_max_cat_ctl, &cat_covar_wkspace) ||
               bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &dbl_2d_buf))) {
    return kPglRetNomem;
  }
  uint32_t* cat_obs_buf;
  // bugfix (11 May 2020): we were previously only initializing cat_obs_buf
  // when extra_cat_ct > 0; this resulted in segfaults when every categorical
  // covariate had only two categories.
  if (unlikely(bigstack_alloc_u32(covar_max_nonnull_cat_ct + 1, &cat_obs_buf))) {
    return kPglRetNomem;
  }
  unsigned char* alloc_base = g_bigstack_base;
  unsigned char* new_covar_name_alloc = g_bigstack_end;
  const uint32_t first_sample_uidx = AdvTo1Bit(sample_include, 0);
  uintptr_t covar_read_uidx_base = 0;
  uintptr_t covar_include_bits = covar_include[0];
  const char** cur_covar_names_iter = cur_covar_names;
  double* covar_write_iter = covars_cmaj;
  double* sum_iter = dbl_2d_buf;
  // Quasi-bugfix (22 Jan 2019): Main loop is too numerically unstable,
  // especially in the single-precision logistic/Firth regression case, when a
  // covariate is on a very different scale from the main genotype column (or
  // other covariates).  We now check sum-of-squares and (rescaled) sample
  // variance of each fixed covariate column, and error out when the range of
  // either value exceeds ~6 orders of magnitude.  (Conveniently, this catches
  // the fairly common "year of birth" problem.)
  // Initial ssq and variance values correspond to genotype columns with a MAF
  // range of [0.01, 0.5].
  double max_ssq = 1.5 * u31tod(sample_ct);
  double min_ssq_minus_sqmean = 0.0198 * u31tod(sample_ct - 1);
  for (uintptr_t covar_read_idx = 0; covar_read_idx != covar_ct; ++covar_read_idx) {
    const uintptr_t covar_read_uidx = BitIter1(covar_include, &covar_read_uidx_base, &covar_include_bits);
    const PhenoCol* cur_covar_col = &(covar_cols[covar_read_uidx]);
    const char* covar_name_base = &(covar_names[covar_read_uidx * max_covar_name_blen]);
    if (cur_covar_col->type_code == kPhenoDtypeOther) {
      // local covariate
      *cur_covar_names_iter++ = covar_name_base;
    } else if (cur_covar_col->type_code == kPhenoDtypeQt) {
      *cur_covar_names_iter++ = covar_name_base;
      const double* covar_vals = cur_covar_col->data.qt;
      uintptr_t sample_uidx_base;
      uintptr_t sample_include_bits;
      BitIter1Start(sample_include, first_sample_uidx, &sample_uidx_base, &sample_include_bits);
      double covar_sum = 0.0;
      double covar_ssq = 0.0;
      for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
        const double cur_covar_val = covar_vals[sample_uidx];
        covar_sum += cur_covar_val;
        covar_ssq += cur_covar_val * cur_covar_val;
        *covar_write_iter++ = cur_covar_val;
      }
      *sum_iter++ = covar_sum;
      if (covar_ssq > max_ssq) {
        max_ssq = covar_ssq;
      }
      covar_ssq -= covar_sum * covar_sum / u31tod(sample_ct);
      if (covar_ssq < min_ssq_minus_sqmean) {
        min_ssq_minus_sqmean = covar_ssq;
      }
    } else {
      // this is equivalent to "--split-cat-pheno omit-most covar-01"
      const uint32_t cur_cat_ct = cur_covar_col->nonnull_category_ct + 1;
      const uint32_t cur_cat_ctl = BitCtToWordCt(cur_cat_ct);
      const uint32_t largest_cat_uidx = IdentifyRemainingCatsAndMostCommon(sample_include, cur_covar_col, sample_ct, cat_covar_wkspace, cat_obs_buf);
      const uint32_t remaining_cat_ct = PopcountWords(cat_covar_wkspace, cur_cat_ctl);
      assert(remaining_cat_ct >= 2);
      ClearBit(largest_cat_uidx, cat_covar_wkspace);
      const uint32_t* covar_vals = cur_covar_col->data.cat;
      const char* const* cur_category_names = cur_covar_col->category_names;
      const uint32_t covar_name_base_slen = strlen(covar_name_base);
      uintptr_t cat_uidx_base;
      uintptr_t cat_covar_wkspace_bits;
      BitIter1Start(cat_covar_wkspace, 1, &cat_uidx_base, &cat_covar_wkspace_bits);
      for (uint32_t cat_idx = 1; cat_idx != remaining_cat_ct; ++cat_idx) {
        const uintptr_t cat_uidx = BitIter1(cat_covar_wkspace, &cat_uidx_base, &cat_covar_wkspace_bits);
        const char* catname = cur_category_names[cat_uidx];
        const uint32_t catname_slen = strlen(catname);
        new_covar_name_alloc -= covar_name_base_slen + catname_slen + 2;
        if (unlikely(new_covar_name_alloc < alloc_base)) {
          return kPglRetNomem;
        }
        char* new_covar_name_write = memcpyax(new_covar_name_alloc, covar_name_base, covar_name_base_slen, '=');
        memcpy(new_covar_name_write, catname, catname_slen + 1);
        *cur_covar_names_iter++ = R_CAST(const char*, new_covar_name_alloc);

        uintptr_t sample_uidx_base;
        uintptr_t sample_include_bits;
        BitIter1Start(sample_include, first_sample_uidx, &sample_uidx_base, &sample_include_bits);
        uint32_t cur_cat_obs_ct = 0;
        for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
          const uint32_t cur_sample_is_in_cat = (covar_vals[sample_uidx] == cat_uidx);
          cur_cat_obs_ct += cur_sample_is_in_cat;
          *covar_write_iter++ = u31tod(cur_sample_is_in_cat);
        }
        double covar_ssq = u31tod(cur_cat_obs_ct);
        *sum_iter++ = covar_ssq;
        if (covar_ssq > max_ssq) {
          max_ssq = covar_ssq;
        }
        covar_ssq -= covar_ssq * covar_ssq / u31tod(sample_ct);
        if (covar_ssq < min_ssq_minus_sqmean) {
          min_ssq_minus_sqmean = covar_ssq;
        }
      }
    }
  }
  BigstackEndSet(new_covar_name_alloc);
  if (min_ssq_minus_sqmean * 1048576.0 < max_ssq) {
    // Could automatically variance-standardize, while keeping track of the
    // linear transformations so we can translate results back to original
    // units in the final output.
    // However, this can have a surprising interaction with the VIF check, so
    // it should come with an option to disable the auto-standardization.
    *glm_err_ptr = SetGlmErr0(kGlmErrcodeUnstableScale);
    return kPglRetSkipped;
  }
  assert(covar_write_iter == &(covars_cmaj[new_nonlocal_covar_ct * sample_ct]));
  MultiplySelfTranspose(covars_cmaj, new_nonlocal_covar_ct, sample_ct, covar_dotprod);
  // intentionally ignore error code, since all callers check glm_err
  *glm_err_ptr = CheckMaxCorrAndVif(covar_dotprod, 0, new_nonlocal_covar_ct, sample_ct, max_corr, vif_thresh, dbl_2d_buf, corr_buf, inverse_corr_buf, matrix_invert_buf1);
  BigstackReset(matrix_invert_buf1);
  return kPglRetSuccess;
}

BoolErr InitNmPrecomp(const double* covars_cmaj, const double* covar_dotprod, const double* corr_buf, const double* corr_inv_tri, uint32_t sample_ct, uint32_t is_qt, uintptr_t new_covar_ct, uintptr_t xtx_state, RegressionNmPrecomp* nm_precomp) {
  uintptr_t stride = new_covar_ct + 1 + xtx_state;
  double* xtx_image = nm_precomp->xtx_image;
  ZeroDArr(stride * stride, xtx_image);
  xtx_image[0] = u31tod(sample_ct);
  for (uintptr_t covar_idx = 0; covar_idx != new_covar_ct; ++covar_idx) {
    const double* cur_covar = &(covars_cmaj[covar_idx * sample_ct]);
    double dxx = 0.0;
    for (uint32_t uii = 0; uii != sample_ct; ++uii) {
      dxx += cur_covar[uii];
    }
    double* xtx_image_row = &(xtx_image[(covar_idx + 1 + xtx_state) * stride]);
    xtx_image_row[0] = dxx;
    memcpy(&(xtx_image_row[1 + xtx_state]), &(covar_dotprod[covar_idx * new_covar_ct]), (covar_idx + 1) * sizeof(double));
  }
  if (is_qt) {
    // also save precomputed inverse of covar dotprods, with intercept included
    // (not currently needed in logistic case)
    double* covarx_dotprod_inv = nm_precomp->covarx_dotprod_inv;
    covarx_dotprod_inv[0] = xtx_image[0];
    const uintptr_t new_covar_ct_p1 = new_covar_ct + 1;
    for (uintptr_t row_idx = 1; row_idx != new_covar_ct_p1; ++row_idx) {
      covarx_dotprod_inv[row_idx * new_covar_ct_p1] = xtx_image[(row_idx + xtx_state) * stride];
      memcpy(&(covarx_dotprod_inv[row_idx * new_covar_ct_p1 + 1]), &(covar_dotprod[(row_idx - 1) * new_covar_ct]), row_idx * sizeof(double));
    }
    // this makes assumptions about amount of space past corr_inv...
    if (InvertSymmdefMatrixChecked(new_covar_ct_p1, covarx_dotprod_inv, R_CAST(MatrixInvertBuf1*, &(nm_precomp->corr_inv[new_covar_ct_p1 * MAXV(new_covar_ct_p1, 3)])), nm_precomp->corr_inv)) {
      return 1;
    }
    ReflectMatrix(new_covar_ct_p1, covarx_dotprod_inv);
    // xt_y_image fill no longer happens in this function (since that's
    // inappropriate in the QT-batch case)
  }

  memcpy(nm_precomp->corr_inv, corr_inv_tri, new_covar_ct * new_covar_ct * sizeof(double));
  ReflectMatrix(new_covar_ct, nm_precomp->corr_inv);
  --stride;
  double* corr_image = nm_precomp->corr_image;
  // also store uninverted correlation matrix
  ZeroDArr(stride * stride, corr_image);
  for (uintptr_t orig_row_idx = 0; orig_row_idx != new_covar_ct; ++orig_row_idx) {
    memcpy(&(corr_image[(orig_row_idx + xtx_state) * stride + xtx_state]), &(corr_buf[orig_row_idx * new_covar_ct]), (orig_row_idx + 1) * sizeof(double));
  }
  // and store inverse-sqrts after the correlation matrix
  memcpy(nm_precomp->corr_inv_sqrts, &(corr_buf[new_covar_ct * new_covar_ct]), new_covar_ct * sizeof(double));
  return 0;
}

uint32_t DosageIsConstant(uint64_t dosage_sum, uint64_t dosage_ssq, uint32_t nm_sample_ct) {
  // Dosages are all identical iff dosage_sum * dosage_sum ==
  // dosage_ssq * nm_sample_ct.
  // Unfortunately, uint64 isn't enough for this comparison.
  const uint64_t dosage_sum_hi = dosage_sum >> 32;
  const uint64_t dosage_sum_lo = S_CAST(uint32_t, dosage_sum);
  const uint64_t dosage_ssq_hi = dosage_ssq >> 32;
  const uint64_t dosage_ssq_lo = S_CAST(uint32_t, dosage_ssq);

  // (a * 2^32 + b)^2 = a^2 * 2^64 + 2ab * 2^32 + b^2
  const uint64_t lhs_ab = dosage_sum_hi * dosage_sum_lo;
  const uint64_t lhs_b2 = dosage_sum_lo * dosage_sum_lo;
  const uint64_t lhs_lo = lhs_b2 + (lhs_ab << 33);
  const uint64_t lhs_hi = (lhs_lo < lhs_b2) + (lhs_ab >> 31) + dosage_sum_hi * dosage_sum_hi;

  // (a * 2^32 + b) * c = ac * 2^32 + bc
  const uint64_t rhs_ac = dosage_ssq_hi * nm_sample_ct;
  const uint64_t rhs_bc = dosage_ssq_lo * nm_sample_ct;
  const uint64_t rhs_lo = rhs_bc + (rhs_ac << 32);
  const uint64_t rhs_hi = (rhs_lo < rhs_bc) + (rhs_ac >> 32);
  return (lhs_hi == rhs_hi) && (lhs_lo == rhs_lo);
}

uint32_t GetBiallelicReportedTestCt(const uintptr_t* parameter_subset, GlmFlags glm_flags, uint32_t covar_ct, uint32_t tests_flag) {
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
  const uint32_t joint_test = domdev_present || tests_flag;

  if (hide_covar) {
    if (!parameter_subset) {
      return 1 + include_intercept + domdev_present + joint_test;
    }
    return include_intercept + domdev_present + joint_test + IsSet(parameter_subset, 1);
  }

  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t biallelic_predictor_ct_base = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
  uint32_t biallelic_predictor_ct = biallelic_predictor_ct_base;
  if (parameter_subset) {
    biallelic_predictor_ct = PopcountWords(parameter_subset, BitCtToWordCt(biallelic_predictor_ct_base));
  }
  return biallelic_predictor_ct + joint_test + include_intercept - 1;
}

typedef struct LocalCovarCoeffparseCtxStruct {
  uint32_t* sample_idx_order;
  uint32_t max_sample_ct;
  uint32_t cur_sample_ct;
  uint32_t tokens_per_sample;
  uint32_t local_covar_ct;
  uint32_t omit_last;
  uint32_t local_haps;
  uint32_t local_cat_ct;
  uint32_t local_cats_1based;
} LocalCovarCoeffparseCtx;

// Processes a single line's worth of payload.
PglErr LoadLocalCovarCoeffs(const LocalCovarCoeffparseCtx* ctx, const char* local_line_iter, uint32_t local_line_idx, float* local_covars_vcmaj_f_iter, double* local_covars_vcmaj_d_iter) {
  const uint32_t* sample_idx_order = ctx->sample_idx_order;
  const uint32_t max_sample_ct = ctx->max_sample_ct;
  const uint32_t cur_sample_ct = ctx->cur_sample_ct;
  const uint32_t tokens_per_sample = ctx->tokens_per_sample;
  const uint32_t local_covar_ct = ctx->local_covar_ct;
  const uint32_t omit_last = ctx->omit_last;
  const uint32_t local_haps = ctx->local_haps;
  const uint32_t local_cat_ct = ctx->local_cat_ct;
  const uint32_t local_cats_1based = ctx->local_cats_1based;
  const uint32_t max_cat_idx = local_cat_ct + local_cats_1based - 1;
  uint32_t sample_idx = 0;
  for (uint32_t local_sample_idx = 0; sample_idx != cur_sample_ct; ++local_sample_idx) {
    const uint32_t cur_sample_idx = sample_idx_order[local_sample_idx];
    if (cur_sample_idx == UINT32_MAX) {
      local_line_iter = NextTokenMult(local_line_iter, tokens_per_sample);
      if (unlikely(!local_line_iter)) {
        logputs("\n");
        logerrprintfww("Error: Fewer tokens than expected on line %u of --glm local-covar= file.\n", local_line_idx);
        return kPglRetMalformedInput;
      }
      continue;
    }
    if (local_cat_ct) {
      uint32_t cat_idx;
      if (unlikely(ScanmovUintCapped(max_cat_idx, &local_line_iter, &cat_idx) || (cat_idx < local_cats_1based))) {
        logputs("\n");
        logerrprintf("Error: Invalid category index on line %u of --glm local-covar= file.\n", local_line_idx);
        return kPglRetMalformedInput;
      }
      cat_idx -= local_cats_1based;
      local_line_iter = FirstNonTspace(FirstSpaceOrEoln(local_line_iter));
      if (!local_haps) {
        if (cat_idx != max_cat_idx) {
          const uint32_t offset = cat_idx * max_sample_ct + cur_sample_idx;
          if (local_covars_vcmaj_f_iter) {
            local_covars_vcmaj_f_iter[offset] = 1.0;
          } else {
            local_covars_vcmaj_d_iter[offset] = 1.0;
          }
        }
      } else {
        if (cat_idx != max_cat_idx) {
          const uint32_t offset = cat_idx * max_sample_ct + cur_sample_idx;
          if (local_covars_vcmaj_f_iter) {
            local_covars_vcmaj_f_iter[offset] = 0.5;
          } else {
            local_covars_vcmaj_d_iter[offset] = 0.5;
          }
        }
        if (unlikely(ScanmovUintCapped(max_cat_idx, &local_line_iter, &cat_idx) || (cat_idx < local_cats_1based))) {
          logputs("\n");
          logerrprintf("Error: Invalid category index on line %u of --glm local-covar= file.\n", local_line_idx);
          return kPglRetMalformedInput;
        }
        cat_idx -= local_cats_1based;
        local_line_iter = FirstNonTspace(FirstSpaceOrEoln(local_line_iter));
        if (cat_idx != max_cat_idx) {
          const uint32_t offset = cat_idx * max_sample_ct + cur_sample_idx;
          if (local_covars_vcmaj_f_iter) {
            local_covars_vcmaj_f_iter[offset] += S_CAST(float, 0.5);
          } else {
            local_covars_vcmaj_d_iter[offset] += 0.5;
          }
        }
      }
    } else {
      if (local_covars_vcmaj_f_iter) {
        float* local_covars_f_iter2 = &(local_covars_vcmaj_f_iter[cur_sample_idx]);
        for (uint32_t covar_idx = 0; covar_idx != local_covar_ct; ++covar_idx) {
          double dxx;
          local_line_iter = ScantokDouble(local_line_iter, &dxx);
          if (unlikely((!local_line_iter) || (fabs(dxx) > 3.4028235677973362e38))) {
            logputs("\n");
            logerrprintf("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
            return kPglRetMalformedInput;
          }
          *local_covars_f_iter2 = S_CAST(float, dxx);
          local_covars_f_iter2 = &(local_covars_f_iter2[max_sample_ct]);
          local_line_iter = FirstNonTspace(local_line_iter);
        }
      } else {
        double* local_covars_d_iter2 = &(local_covars_vcmaj_d_iter[cur_sample_idx]);
        for (uint32_t covar_idx = 0; covar_idx != local_covar_ct; ++covar_idx) {
          double dxx;
          local_line_iter = ScantokDouble(local_line_iter, &dxx);
          if (unlikely(!local_line_iter)) {
            logputs("\n");
            logerrprintf("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
            return kPglRetMalformedInput;
          }
          *local_covars_d_iter2 = dxx;
          local_covars_d_iter2 = &(local_covars_d_iter2[max_sample_ct]);
          local_line_iter = FirstNonTspace(local_line_iter);
        }
      }
      if (omit_last) {
        local_line_iter = FirstNonTspace(FirstSpaceOrEoln(local_line_iter));
      }
      if (local_haps) {
        if (local_covars_vcmaj_f_iter) {
          float* local_covars_f_iter2 = &(local_covars_vcmaj_f_iter[cur_sample_idx]);
          for (uint32_t covar_idx = 0; covar_idx != local_covar_ct; ++covar_idx) {
            double dxx;
            local_line_iter = ScantokDouble(local_line_iter, &dxx);
            if (unlikely((!local_line_iter) || (fabs(dxx) > 3.4028235677973362e38))) {
              logputs("\n");
              logerrprintf("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
              return kPglRetMalformedInput;
            }
            *local_covars_f_iter2 = S_CAST(float, (S_CAST(double, *local_covars_f_iter2) + dxx) * 0.5);
            local_covars_f_iter2 = &(local_covars_f_iter2[max_sample_ct]);
            local_line_iter = FirstNonTspace(local_line_iter);
          }
        } else {
          double* local_covars_d_iter2 = &(local_covars_vcmaj_d_iter[cur_sample_idx]);
          for (uint32_t covar_idx = 0; covar_idx != local_covar_ct; ++covar_idx) {
            double dxx;
            local_line_iter = ScantokDouble(local_line_iter, &dxx);
            if (unlikely(!local_line_iter)) {
              logputs("\n");
              logerrprintf("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
              return kPglRetMalformedInput;
            }
            // may as well defend against overflow
            *local_covars_d_iter2 = (*local_covars_d_iter2) * 0.5 + dxx * 0.5;
            local_covars_d_iter2 = &(local_covars_d_iter2[max_sample_ct]);
            local_line_iter = FirstNonTspace(local_line_iter);
          }
        }
        if (omit_last) {
          local_line_iter = FirstNonTspace(FirstSpaceOrEoln(local_line_iter));
        }
      }
    }
    ++sample_idx;
  }
  return kPglRetSuccess;
}

PglErr ReadLocalCovarBlock(const GlmCtx* common, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t cur_block_variant_ct, uint32_t local_sample_ct, uint32_t local_cat_ct, TextStream* local_covar_txsp, uint32_t* local_line_idx_ptr, uint32_t* local_xy_ptr, float* local_covars_vcmaj_f_iter, double* local_covars_vcmaj_d_iter, uint32_t* local_sample_idx_order) {
  const ChrInfo* cip = common->cip;
  const uintptr_t* variant_include = common->variant_include;
  const uint32_t sample_ct = common->sample_ct;
  const uint32_t sample_ct_x = common->sample_ct_x;
  const uint32_t sample_ct_y = common->sample_ct_y;
  const uint32_t local_covar_ct = common->local_covar_ct;
  const GlmFlags flags = common->glm_flags;
  const uint32_t omit_last = (flags / kfGlmLocalOmitLast) & 1;
  const uint32_t local_haps = (flags / kfGlmLocalHaps) & 1;

  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
  if (max_sample_ct < sample_ct_y) {
    max_sample_ct = sample_ct_y;
  }
  LocalCovarCoeffparseCtx coeffparse_ctx;
  coeffparse_ctx.sample_idx_order = local_sample_idx_order;
  coeffparse_ctx.max_sample_ct = max_sample_ct;
  // cur_sample_ct filled a bit later
  coeffparse_ctx.tokens_per_sample = (local_cat_ct? 1 : (local_covar_ct + omit_last)) << local_haps;
  coeffparse_ctx.local_covar_ct = local_covar_ct;
  coeffparse_ctx.omit_last = omit_last;
  coeffparse_ctx.local_haps = local_haps;
  coeffparse_ctx.local_cat_ct = local_cat_ct;
  coeffparse_ctx.local_cats_1based = (flags / kfGlmLocalCats1based) & 1;
  uint32_t variant_bidx = 0;
  if (local_cat_ct) {
    // assert(local_covar_ct == local_cat_ct - 1);
    if (local_covars_vcmaj_f_iter) {
      ZeroFArr(local_covar_ct * max_sample_ct * S_CAST(uintptr_t, cur_block_variant_ct), local_covars_vcmaj_f_iter);
    } else {
      ZeroDArr(local_covar_ct * max_sample_ct * S_CAST(uintptr_t, cur_block_variant_ct), local_covars_vcmaj_d_iter);
    }
  }
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, variant_uidx_start, &variant_uidx_base, &cur_bits);
  uint32_t local_line_idx = *local_line_idx_ptr;
  while (variant_bidx < cur_block_variant_ct) {
    const uint32_t variant_uidx1 = BitIter1NoAdv(variant_include, &variant_uidx_base, &cur_bits);
    const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx1);
    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
    const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    uint32_t cur_variant_bidx_end = cur_block_variant_ct;
    if (chr_end < variant_uidx_end) {
      cur_variant_bidx_end = variant_bidx + PopcountBitRange(variant_include, variant_uidx1, chr_end);
      assert(cur_variant_bidx_end <= cur_block_variant_ct);
    }
    const uint32_t is_x = (chr_idx == x_code);
    const uint32_t is_y = (chr_idx == y_code);
    const uintptr_t* cur_sample_include;
    const uint32_t* cur_sample_include_cumulative_popcounts;
    if (is_y && common->sample_include_y) {
      cur_sample_include = common->sample_include_y;
      cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
      coeffparse_ctx.cur_sample_ct = sample_ct_y;
    } else if (is_x && common->sample_include_x) {
      cur_sample_include = common->sample_include_x;
      cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
      coeffparse_ctx.cur_sample_ct = sample_ct_x;
    } else {
      cur_sample_include = common->sample_include;
      cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
      coeffparse_ctx.cur_sample_ct = sample_ct;
    }
    const uint32_t new_local_xy = is_x + 2 * is_y;
    if (new_local_xy != *local_xy_ptr) {
      for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
        const uint32_t cur_uidx = local_sample_uidx_order[uii];
        uint32_t cur_idx = UINT32_MAX;
        if ((cur_uidx != UINT32_MAX) && IsSet(cur_sample_include, cur_uidx)) {
          cur_idx = RawToSubsettedPos(cur_sample_include, cur_sample_include_cumulative_popcounts, cur_uidx);
        }
        local_sample_idx_order[uii] = cur_idx;
      }
      *local_xy_ptr = new_local_xy;
    }
    for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
      BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (!IsSet(local_variant_include, local_line_idx)) {
        uint32_t local_line_idx_target_m1 = AdvTo1Bit(local_variant_include, local_line_idx);
        PglErr reterr = TextSkipNz(local_line_idx_target_m1 - local_line_idx, local_covar_txsp);
        if (unlikely(reterr)) {
          if (reterr == kPglRetEof) {
            logputs("\n");
            logerrputs("Error: --glm local-covar= file has fewer lines than local-pvar= file.\n");
            return kPglRetInconsistentInput;
          }
          TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
          return reterr;
        }
        local_line_idx = local_line_idx_target_m1;
      }
      ++local_line_idx;
      PglErr reterr = kPglRetSuccess;
      char* local_covar_line_start = TextGet(local_covar_txsp);
      if (unlikely(!local_covar_line_start)) {
        if (!TextStreamErrcode2(local_covar_txsp, &reterr)) {
          logputs("\n");
          logerrputs("Error: --glm local-covar= file has fewer lines than local-pvar= file.\n");
          return kPglRetInconsistentInput;
        }
        TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
        return reterr;
      }
      reterr = LoadLocalCovarCoeffs(&coeffparse_ctx, local_covar_line_start, local_line_idx, local_covars_vcmaj_f_iter, local_covars_vcmaj_d_iter);
      if (unlikely(reterr)) {
        return reterr;
      }
      if (local_covars_vcmaj_f_iter) {
        local_covars_vcmaj_f_iter += max_sample_ct * local_covar_ct;
      } else {
        local_covars_vcmaj_d_iter += max_sample_ct * local_covar_ct;
      }
    }
  }
  *local_line_idx_ptr = local_line_idx;
  return kPglRetSuccess;
}

static inline void ZeroLocalCovarRows(uint32_t row_ct, uintptr_t row_width, uint32_t* variant_bidx_ptr, float** local_covars_vcmaj_f_iterp, double** local_covars_vcmaj_d_iterp) {
  *variant_bidx_ptr += row_ct;
  const uintptr_t elem_ct = row_ct * row_width;
  if (*local_covars_vcmaj_f_iterp) {
    ZeroFArr(elem_ct, *local_covars_vcmaj_f_iterp);
    *local_covars_vcmaj_f_iterp += elem_ct;
  } else {
    ZeroDArr(elem_ct, *local_covars_vcmaj_d_iterp);
    *local_covars_vcmaj_d_iterp += elem_ct;
  }
}

void DuplicateLocalCovarRow(uint32_t row_ct, uintptr_t row_width, uint32_t* variant_bidx_ptr, float** local_covars_vcmaj_f_iterp, double** local_covars_vcmaj_d_iterp) {
  if (*local_covars_vcmaj_f_iterp) {
    float* local_covars_vcmaj_f_iter = *local_covars_vcmaj_f_iterp;
    // Could also copy one row at a time, but I'd expect doubling to be
    // slightly more efficient.
    for (uint32_t row_idx = 1; row_idx < row_ct; row_idx *= 2) {
      uint32_t row_copy_ct = row_ct - row_idx;
      if (row_copy_ct > row_idx) {
        row_copy_ct = row_idx;
      }
      memcpy(&(local_covars_vcmaj_f_iter[row_idx * row_width]), local_covars_vcmaj_f_iter, row_copy_ct * row_width * sizeof(float));
    }
    local_covars_vcmaj_f_iter = &(local_covars_vcmaj_f_iter[row_ct * row_width]);
  } else {
    double* local_covars_vcmaj_d_iter = *local_covars_vcmaj_d_iterp;
    for (uint32_t row_idx = 1; row_idx < row_ct; row_idx *= 2) {
      uint32_t row_copy_ct = row_ct - row_idx;
      if (row_copy_ct > row_idx) {
        row_copy_ct = row_idx;
      }
      memcpy(&(local_covars_vcmaj_d_iter[row_idx * row_width]), local_covars_vcmaj_d_iter, row_copy_ct * row_width * sizeof(double));
    }
    *local_covars_vcmaj_d_iterp = &(local_covars_vcmaj_d_iter[row_ct * row_width]);
  }
  *variant_bidx_ptr += row_ct;
}


PglErr ReadRfmix2Block(const GlmCtx* common, const uint32_t* variant_bps, const uint32_t* local_sample_uidx_order, const float* prev_local_covar_row_f, const double* prev_local_covar_row_d, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t cur_block_variant_ct, uint32_t local_sample_ct, uint32_t local_cat_ct, uint32_t local_chrom_col, uint32_t local_bp_col, uint32_t local_first_covar_col, TextStream* local_covar_txsp, const char** local_line_iterp, uint32_t* local_line_idx_ptr, uint32_t* local_prev_chr_code_ptr, uint32_t* local_chr_code_ptr, uint32_t* local_bp_ptr, uint32_t* local_skip_chr_ptr, float* local_covars_vcmaj_f_iter, double* local_covars_vcmaj_d_iter, uint32_t* local_sample_idx_order) {
  const ChrInfo* cip = common->cip;
  const uintptr_t* variant_include = common->variant_include;
  const uint32_t sample_ct = common->sample_ct;
  const uint32_t sample_ct_x = common->sample_ct_x;
  const uint32_t sample_ct_y = common->sample_ct_y;
  const uint32_t local_covar_ct = common->local_covar_ct;
  const GlmFlags flags = common->glm_flags;
  const uint32_t omit_last = (flags / kfGlmLocalOmitLast) & 1;
  const uint32_t local_haps = (flags / kfGlmLocalHaps) & 1;
  // There are several complications here:
  // 1. Until we've seen the beginning of the next line, we don't know the end
  //    of the (possibly empty) variant range the current set of local
  //    covariates applies to.
  // 2. The aforementioned variant range may go past the end of the current
  //    variant block.
  // 3. There may be variants which aren't contained in any local-covar
  //    interval (either before the first local-covar position on that
  //    chromosome, or on a chromosome absent from the local-covar file).
  // To address (1), we structure the main loop such that the first part
  // duplicates row contents and advances the matrix-filling iterator when
  // necessary, the second part speculatively parses local covariates to the
  // next matrix row, and the last part reads the position fields on the next
  // line.
  // To manage (2), we track local_prev_chr_code, local_chr_code, local_bp, and
  // prev_local_covar_row across function calls.  When local_prev_chr_code !=
  // local_chr_code on function restart, we know we need to replicate
  // prev_local_covar_row to the end of that chromosome, exiting early if the
  // end of the chromosome is past variant_uidx_end; otherwise, we need to
  // replicate it up to the first variant with position >= local_bp, possibly
  // exiting early.  This is no different from within-variant-block iteration,
  // so the main loop is written such that function reentry isn't
  // special-cased outside of needing to know prev_local_covar_row.
  // For (3), we zero-fill the affected local covariate rows.
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
  if (max_sample_ct < sample_ct_y) {
    max_sample_ct = sample_ct_y;
  }

  LocalCovarCoeffparseCtx coeffparse_ctx;
  coeffparse_ctx.sample_idx_order = local_sample_idx_order;
  coeffparse_ctx.max_sample_ct = max_sample_ct;
  // cur_sample_ct filled a bit later
  coeffparse_ctx.tokens_per_sample = (local_cat_ct? 1 : (local_covar_ct + omit_last)) << local_haps;
  coeffparse_ctx.local_covar_ct = local_covar_ct;
  coeffparse_ctx.omit_last = omit_last;
  coeffparse_ctx.local_haps = local_haps;
  coeffparse_ctx.local_cat_ct = local_cat_ct;
  coeffparse_ctx.local_cats_1based = (flags / kfGlmLocalCats1based) & 1;

  const char* local_line_iter = *local_line_iterp;
  uint32_t local_line_idx = *local_line_idx_ptr;
  uint32_t first_skip;
  uint32_t second_skip;
  uint32_t last_skip;
  if (local_chrom_col < local_bp_col) {
    first_skip = local_chrom_col - 1;
    second_skip = local_bp_col - local_chrom_col;
    last_skip = local_first_covar_col - local_bp_col;
  } else {
    first_skip = local_bp_col - 1;
    second_skip = local_chrom_col - local_bp_col;
    last_skip = local_first_covar_col - local_chrom_col;
  }
  const uintptr_t row_width = local_covar_ct * max_sample_ct;
  uint32_t local_prev_chr_code = *local_prev_chr_code_ptr;
  uint32_t local_chr_code = *local_chr_code_ptr;
  uint32_t local_bp = *local_bp_ptr;
  uint32_t local_skip_chr = *local_skip_chr_ptr;
  // might be inaccurate at file-start and EOF, but that's okay since these
  // variables are reinitialized before use in the first case, and irrelevant
  // in the second
  uint32_t is_x = (local_chr_code == x_code);
  uint32_t is_y = (local_chr_code == y_code);
  uint32_t local_xy = is_x + 2 * is_y;
  const uintptr_t* cur_sample_include;
  const uint32_t* cur_sample_include_cumulative_popcounts;
  if (is_y && common->sample_include_y) {
    cur_sample_include = common->sample_include_y;
    cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
    coeffparse_ctx.cur_sample_ct = sample_ct_y;
  } else if (is_x && common->sample_include_x) {
    cur_sample_include = common->sample_include_x;
    cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
    coeffparse_ctx.cur_sample_ct = sample_ct_x;
  } else {
    cur_sample_include = common->sample_include;
    cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
    coeffparse_ctx.cur_sample_ct = sample_ct;
  }
  // Necessary for cross-block duplication to work.
  if (prev_local_covar_row_f) {
    memcpy(local_covars_vcmaj_f_iter, prev_local_covar_row_f, row_width * sizeof(float));
  } else if (prev_local_covar_row_d) {
    memcpy(local_covars_vcmaj_d_iter, prev_local_covar_row_d, row_width * sizeof(double));
  }

  uint32_t variant_bidx = 0;
  uint32_t variant_uidx = AdvTo1Bit(variant_include, variant_uidx_start);
  uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
  uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
  uint32_t variant_uidx_chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  if (variant_uidx_chr_end > variant_uidx_end) {
    variant_uidx_chr_end = variant_uidx_end;
  }
  uint32_t variant_bp = (chr_idx == local_prev_chr_code)? variant_bps[variant_uidx] : UINT32_MAX;
  while (1) {
    // Loop invariants:
    //   local_chr_code and local_bp are from the current line, if we're not at
    //     the beginning or end of the file.  In the latter two cases,
    //     local_chr_code is UINT32_MAX.
    //   local_prev_chr_code is from the previous relevant line, or UINT32_MAX
    //     at the beginning of the file.
    //   variant_uidx = next variant we need to fill local covariates for.
    //   chr_fo_idx and chr_idx correspond to variant_uidx.
    //   variant_uidx_chr_end is min(end of chr_fo_idx, variant_uidx_end).
    //   variant_bp is UINT32_MAX if we've already iterated through all
    //     variants on local_prev_chr_code (or we're at the beginning of the
    //     file); otherwise it also corresponds to variant_uidx.
    //   local_skip_chr is set iff either all variants on local_chr_code were
    //     filtered out, or local_chr_code == UINT32_MAX.
    //   If local_skip_chr isn't true, variant_uidx is not before the start of
    //     local_prev_chr_code, or after the end of local_chr_code.
    if (!local_skip_chr) {
      // Part 1.
      const uint32_t backfill_stop_bp = (local_chr_code == chr_idx)? local_bp : UINT32_MAX;
      if (backfill_stop_bp > variant_bp) {
        // At least one variant in [prev pos, current pos) (or [prev pos, end
        // of chromosome) if we're moving on to a new chromosome).  Count how
        // many there are, and fill that many lines of the matrix.
        const uint32_t next_variant_uidx = ExpsearchU32(variant_bps, variant_uidx, variant_uidx_chr_end, backfill_stop_bp);
        const uint32_t row_ct = PopcountBitRange(variant_include, variant_uidx, next_variant_uidx);
        DuplicateLocalCovarRow(row_ct, row_width, &variant_bidx, &local_covars_vcmaj_f_iter, &local_covars_vcmaj_d_iter);
        if (variant_bidx == cur_block_variant_ct) {
          break;
        }
        variant_uidx = AdvTo1Bit(variant_include, next_variant_uidx);
        if (variant_uidx >= variant_uidx_chr_end) {
          // probable todo: helper function for this
          chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
          chr_idx = cip->chr_file_order[chr_fo_idx];
          variant_uidx_chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          if (variant_uidx_chr_end > variant_uidx_end) {
            variant_uidx_chr_end = variant_uidx_end;
          }
          variant_bp = UINT32_MAX;
        }
      }
      if (local_chr_code != local_prev_chr_code) {
        const uint32_t local_chr_fo_idx = cip->chr_idx_to_foidx[local_chr_code];
        if (local_chr_code != chr_idx) {
          // Some variants are on chromosomes entirely absent from the
          // local-covar file.  Zero-fill the covariate values for these
          // variants; NA results will be reported for them.
          const uint32_t local_chr_start_vidx = cip->chr_fo_vidx_start[local_chr_fo_idx];
          uint32_t row_ct;
          if (local_chr_start_vidx >= variant_uidx_chr_end) {
            row_ct = cur_block_variant_ct - variant_bidx;
          } else {
            // Verify that we aren't going backwards.  (This condition works
            // since we don't enter the block in the local_skip_chr case.)
            if (unlikely(local_chr_start_vidx < variant_uidx)) {
              logputs("\n");
              logerrputs("Error: --glm local-covar= file has a different chromosome order than the main\ndataset.\n");
              return kPglRetInconsistentInput;
            }
            row_ct = PopcountBitRange(variant_include, variant_uidx, local_chr_start_vidx);
          }
          ZeroLocalCovarRows(row_ct, row_width, &variant_bidx, &local_covars_vcmaj_f_iter, &local_covars_vcmaj_d_iter);
          if (variant_bidx == cur_block_variant_ct) {
            break;
          }
          variant_uidx = AdvTo1Bit(variant_include, local_chr_start_vidx);
          chr_fo_idx = local_chr_fo_idx;
          chr_idx = local_chr_code;
          variant_uidx_chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          if (variant_uidx_chr_end > variant_uidx_end) {
            variant_uidx_chr_end = variant_uidx_end;
          }
        }
        variant_bp = variant_bps[variant_uidx];
      }
      if (local_cat_ct) {
        if (local_covars_vcmaj_f_iter) {
          ZeroFArr(local_covar_ct * max_sample_ct, local_covars_vcmaj_f_iter);
        } else {
          ZeroDArr(local_covar_ct * max_sample_ct, local_covars_vcmaj_d_iter);
        }
      }
      // Part 2.
      PglErr reterr = LoadLocalCovarCoeffs(&coeffparse_ctx, local_line_iter, local_line_idx, local_covars_vcmaj_f_iter, local_covars_vcmaj_d_iter);
      if (unlikely(reterr)) {
        return reterr;
      }
    }
    // Part 3.

    // Not const for now, due to limitation of GetChrCodeCounted().
    char* line_start = TextGet(local_covar_txsp);
    if (!line_start) {
      PglErr reterr = TextStreamErrcode(local_covar_txsp);
      if (unlikely(reterr)) {
        return reterr;
      }
      // EOF.
      local_line_iter = nullptr;
      local_chr_code = UINT32_MAX;
      local_skip_chr = 1;
      if (local_prev_chr_code == chr_idx) {
        uint32_t row_ct = cur_block_variant_ct - variant_bidx;
        if (variant_uidx_chr_end < variant_uidx_end) {
          row_ct = PopcountBitRange(variant_include, variant_uidx, variant_uidx_chr_end);
        }
        if (row_ct) {
          DuplicateLocalCovarRow(row_ct, row_width, &variant_bidx, &local_covars_vcmaj_f_iter, &local_covars_vcmaj_d_iter);
          if (variant_bidx == cur_block_variant_ct) {
            break;
          }
        }
        variant_uidx = AdvTo1Bit(variant_include, variant_uidx_chr_end);
        chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
        chr_idx = cip->chr_file_order[chr_fo_idx];
        // no need to update variant_uidx_chr_end or variant_bp
      }
      const uintptr_t remaining_variant_ct = cur_block_variant_ct - variant_bidx;
      ZeroLocalCovarRows(remaining_variant_ct, row_width, &variant_bidx, &local_covars_vcmaj_f_iter, &local_covars_vcmaj_d_iter);
      break;
    }
    ++local_line_idx;
    char* tok1_start = NextTokenMult0(line_start, first_skip);
    if (unlikely(!tok1_start)) {
      logputs("\n");
      logerrprintf("Error: Line %u of --glm local-covar= file has fewer tokens than expected.\n", local_line_idx);
      return kPglRetMalformedInput;
    }
    char* tok1_end = CurTokenEnd(tok1_start);
    char* tok2_start = NextTokenMult(tok1_end, second_skip);
    if (unlikely(!tok2_start)) {
      logputs("\n");
      logerrprintf("Error: Line %u of --glm local-covar= file has fewer tokens than expected.\n", local_line_idx);
      return kPglRetMalformedInput;
    }
    char* tok2_end = CurTokenEnd(tok2_start);
    local_line_iter = NextTokenMult(tok2_end, last_skip);
    if (unlikely(!local_line_iter)) {
      logputs("\n");
      logerrprintf("Error: Line %u of --glm local-covar= file has fewer tokens than expected.\n", local_line_idx);
      return kPglRetMalformedInput;
    }
    char* chr_code_start = tok1_start;
    char* chr_code_end = tok1_end;
    char* bp_col_start = tok2_start;
    if (local_chrom_col > local_bp_col) {
      chr_code_start = tok2_start;
      chr_code_end = tok2_end;
      bp_col_start = tok1_start;
    }
    const uint32_t next_chr_code = GetChrCodeCounted(cip, chr_code_end - chr_code_start, chr_code_start);
    if (next_chr_code != local_chr_code) {
      local_skip_chr = IsI32Neg(next_chr_code) || (!IsSet(cip->chr_mask, next_chr_code));
      if (local_skip_chr) {
        continue;
      }
      // Still possible for all variants in this chromosome to have been
      // filtered out.
      const uint32_t new_chr_fo_idx = cip->chr_idx_to_foidx[next_chr_code];
      if (local_chr_code != UINT32_MAX) {
        // May as well sanity-check this.
        const uint32_t old_chr_fo_idx = cip->chr_idx_to_foidx[local_chr_code];
        if (unlikely(old_chr_fo_idx > new_chr_fo_idx)) {
          logputs("\n");
          logerrputs("Error: --glm local-covar= file has a different chromosome order than the main\ndataset.\n");
          return kPglRetInconsistentInput;
        }
      }
      const uint32_t chr_start_vidx = cip->chr_fo_vidx_start[new_chr_fo_idx];
      const uint32_t chr_end_vidx = cip->chr_fo_vidx_start[new_chr_fo_idx + 1];
      local_skip_chr = !PopcountBitRange(variant_include, chr_start_vidx, chr_end_vidx);
      if (local_skip_chr) {
        continue;
      }
      // May as well ensure that local_prev_chr_code is either UINT32_MAX,
      // or corresponds to a not-totally-excluded chromosome, so failures
      // of the sanity check above are more likely to be meaningful.
      local_prev_chr_code = local_chr_code;
      local_chr_code = next_chr_code;
      local_bp = UINT32_MAX;
      is_x = (local_chr_code == x_code);
      is_y = (local_chr_code == y_code);
      const uint32_t new_local_xy = is_x + 2 * is_y;
      if (new_local_xy != local_xy) {
        local_xy = new_local_xy;
        if (is_y && common->sample_include_y) {
          cur_sample_include = common->sample_include_y;
          cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
          coeffparse_ctx.cur_sample_ct = sample_ct_y;
        } else if (is_x && common->sample_include_x) {
          cur_sample_include = common->sample_include_x;
          cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
          coeffparse_ctx.cur_sample_ct = sample_ct_x;
        } else {
          cur_sample_include = common->sample_include;
          cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
          coeffparse_ctx.cur_sample_ct = sample_ct;
        }
        for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
          const uint32_t cur_uidx = local_sample_uidx_order[uii];
          uint32_t cur_idx = UINT32_MAX;
          if ((cur_uidx != UINT32_MAX) && IsSet(cur_sample_include, cur_uidx)) {
            cur_idx = RawToSubsettedPos(cur_sample_include, cur_sample_include_cumulative_popcounts, cur_uidx);
          }
          local_sample_idx_order[uii] = cur_idx;
        }
      }
    } else if (local_skip_chr) {
      continue;
    }
    uint32_t next_bp;
    if (unlikely(ScanPosintDefcap(bp_col_start, &next_bp))) {
      logputs("\n");
      logerrprintf("Error: Line %u of --glm local-covar= file has fewer tokens than expected.\n", local_line_idx);
      return kPglRetMalformedInput;
    }
    // could conditionally prohibit duplicate positions
    if (unlikely(S_CAST(int32_t, next_bp) < S_CAST(int32_t, local_bp))) {
      logputs("\n");
      logerrputs("Error: Positions in --glm local-covar= file are not sorted.\n");
      return kPglRetMalformedInput;
    }
    local_bp = next_bp;
  }
  *local_line_iterp = local_line_iter;
  *local_line_idx_ptr = local_line_idx;
  *local_prev_chr_code_ptr = local_prev_chr_code;
  *local_chr_code_ptr = local_chr_code;
  *local_bp_ptr = local_bp;
  *local_skip_chr_ptr = local_skip_chr;
  return kPglRetSuccess;
}

const double kSmallDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, 3.0);

const double kSmallInvDoublePairs[32] ALIGNV16 = PAIR_TABLE16(2.0, 1.0, 0.0, 3.0);

const double kSmallInvDoubles[4] = {2.0, 1.0, 0.0, 3.0};

uint32_t GenoarrToDoublesRemoveMissing(const uintptr_t* genoarr, const double* __restrict table, uint32_t sample_ct, double* __restrict dst) {
  assert(sample_ct);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t subgroup_len = kBitsPerWordD2;
  double* dst_iter = dst;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        return dst_iter - dst;
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      const uintptr_t cur_geno = geno_word & 3;
      if (cur_geno < 3) {
        // *dst_iter++ = u31tod(cur_geno);
        *dst_iter++ = table[cur_geno];
      }
      geno_word >>= 2;
    }
  }
}

BoolErr LinearHypothesisChisq(const double* coef, const double* constraints_con_major, const double* cov_matrix, uintptr_t constraint_ct, uintptr_t predictor_ct, uintptr_t cov_stride, double* chisq_ptr, double* tmphxs_buf, double* h_transpose_buf, double* inner_buf, MatrixInvertBuf1* mi_buf, double* outer_buf) {
  // See PLINK model.cpp Model::linearHypothesis().
  //
  // outer_buf = constraint_ct
  // inner_buf = constraint_ct x constraint_ct
  // tmphxs_buf and h_transpose_buf are constraint_ct x predictor_ct
  // mi_buf only needs to be of length 2 * constraint_ct
  //
  // Since no PLINK function ever calls this with nonzero h[] values, this just
  // takes a df (constraint_ct) parameter for now; it's trivial to switch to
  // the more general interface later.
  ColMajorVectorMatrixMultiplyStrided(coef, constraints_con_major, predictor_ct, predictor_ct, constraint_ct, outer_buf);
  MatrixTransposeCopy(constraints_con_major, constraint_ct, predictor_ct, h_transpose_buf);
  ColMajorMatrixMultiplyStrided(h_transpose_buf, cov_matrix, constraint_ct, constraint_ct, predictor_ct, cov_stride, predictor_ct, constraint_ct, tmphxs_buf);
  // tmp[][] is now predictor-major
  ColMajorMatrixMultiply(tmphxs_buf, constraints_con_major, constraint_ct, constraint_ct, predictor_ct, inner_buf);

  // don't need H-transpose any more, so we can use h_transpose_buf for matrix
  // inversion
  if (InvertMatrix(constraint_ct, inner_buf, mi_buf, h_transpose_buf)) {
    return 1;
  }
  double result = 0.0;
  const double* inner_iter = inner_buf;
  if (constraint_ct > kDotprodDThresh) {
    for (uintptr_t constraint_idx = 0; constraint_idx != constraint_ct; ++constraint_idx) {
      result += DotprodD(inner_iter, outer_buf, constraint_ct) * outer_buf[constraint_idx];
      inner_iter = &(inner_iter[constraint_ct]);
    }
  } else {
    for (uintptr_t constraint_idx = 0; constraint_idx != constraint_ct; ++constraint_idx) {
      result += DotprodDShort(inner_iter, outer_buf, constraint_ct) * outer_buf[constraint_idx];
      inner_iter = &(inner_iter[constraint_ct]);
    }
  }
  if (result < 0.0) {
    // guard against floating point error
    result = 0.0;
  }
  *chisq_ptr = result;
  return 0;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
