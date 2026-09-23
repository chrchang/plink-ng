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

#include "plink2_glm_multinomial.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include "include/pgenlib_misc.h"
#include "include/plink2_bits.h"
#include "include/plink2_stats.h"
#include "include/plink2_string.h"
#include "include/plink2_thread.h"
#include "plink2_compress_stream.h"
#include "plink2_matrix.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Newton-Raphson settings.  A fit has converged once a step changes no
// coefficient by more than kMultinomialStepTol and raises the log-likelihood
// by less than kMultinomialLnLikTol (plus a rounding allowance proportional
// to |log-likelihood|, which matters for very large sample sizes).
static const uint32_t kMultinomialMaxIter = 30;
// The covariate-only fit starts from the intercepts alone, so it gets more
// room.
static const uint32_t kMultinomialNullMaxIter = 100;
static const uint32_t kMultinomialMaxHalvings = 30;
static const double kMultinomialStepTol = 1e-6;
static const double kMultinomialLnLikTol = 1e-8;
static const double kMultinomialLnLikRelRounding = 1e-13;

// The Fisher information is accumulated over blocks of this many samples, so
// its per-sample products take O(param_ct) memory per block instead of
// O(param_ct * sample_ct).
CONSTI32(kMultinomialChunkSize, 256);

uint32_t CountPhenoCats(const uintptr_t* sample_include, const PhenoCol* pheno_col, uint32_t sample_ct) {
  const uint32_t nonnull_cat_ct = pheno_col->nonnull_category_ct;
  const uint32_t* cats = pheno_col->data.cat;
  // A bitset on the stack covers the usual case; category counts are
  // unbounded, so fall back to one scan per category when there are many.
  uint32_t seen_ct = 0;
  if (nonnull_cat_ct < 4096) {
    uintptr_t seen[4096 / kBitsPerWord];
    ZeroWArr(BitCtToWordCt(nonnull_cat_ct + 1), seen);
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      const uint32_t cat_idx = cats[sample_uidx];
      if (!IsSet(seen, cat_idx)) {
        SetBit(cat_idx, seen);
        ++seen_ct;
      }
    }
    return seen_ct - IsSet(seen, 0);
  }
  for (uint32_t cat_idx = 1; cat_idx <= nonnull_cat_ct; ++cat_idx) {
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      if (cats[sample_uidx] == cat_idx) {
        ++seen_ct;
        break;
      }
    }
  }
  return seen_ct;
}

BoolErr GlmMultinomialInitLevels(const uintptr_t* sample_include, const PhenoCol* pheno_col, uint32_t sample_ct, const char* ref_name, uint32_t* level_ct_ptr, uint32_t** cat_to_level_ptr, uint32_t** level_cat_idxs_ptr, uint32_t* ref_found_ptr) {
  const uint32_t nonnull_cat_ct = pheno_col->nonnull_category_ct;
  uint32_t* cat_to_level;
  uint32_t* level_cat_idxs;
  if (unlikely(bigstack_calloc_u32(nonnull_cat_ct + 1, &cat_to_level) ||
               bigstack_alloc_u32(nonnull_cat_ct + 1, &level_cat_idxs))) {
    return 1;
  }
  // count samples per category first
  const uint32_t* cats = pheno_col->data.cat;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    cat_to_level[cats[sample_uidx]] += 1;
  }
  uint32_t ref_cat_idx = 0;
  *ref_found_ptr = 1;
  if (ref_name) {
    *ref_found_ptr = 0;
    for (uint32_t cat_idx = 1; cat_idx <= nonnull_cat_ct; ++cat_idx) {
      if (cat_to_level[cat_idx] && (!strcmp(pheno_col->category_names[cat_idx], ref_name))) {
        ref_cat_idx = cat_idx;
        *ref_found_ptr = 1;
        break;
      }
    }
  }
  if (!ref_cat_idx) {
    for (uint32_t cat_idx = 1; cat_idx <= nonnull_cat_ct; ++cat_idx) {
      if (cat_to_level[cat_idx]) {
        ref_cat_idx = cat_idx;
        break;
      }
    }
  }
  uint32_t level_ct = 0;
  if (ref_cat_idx) {
    level_cat_idxs[0] = ref_cat_idx;
    level_ct = 1;
    for (uint32_t cat_idx = 1; cat_idx <= nonnull_cat_ct; ++cat_idx) {
      if (cat_idx == ref_cat_idx) {
        cat_to_level[cat_idx] = 0;
      } else if (cat_to_level[cat_idx]) {
        cat_to_level[cat_idx] = level_ct;
        level_cat_idxs[level_ct++] = cat_idx;
      } else {
        cat_to_level[cat_idx] = UINT32_MAX;
      }
    }
  }
  cat_to_level[0] = UINT32_MAX;
  *level_ct_ptr = level_ct;
  *cat_to_level_ptr = cat_to_level;
  *level_cat_idxs_ptr = level_cat_idxs;
  return 0;
}

// pred_rows: pred_ct pointers to rows of sample_ct predictor values.
// classes: per-sample class index in [0, class_ct], 0 = reference class.
// coefs: class-major, class_ct x pred_ct.
// Sets *ln_lik_ptr.  If grad is non-null, also fills grad[] with the gradient
// of the log-likelihood, and the lower triangle of info[] (row-major, stride
// param_ct) with the Fisher information
//   I[(k,j),(l,m)] = sum_i x_ij x_im p_ik (delta_kl - p_il),
// computed as the block-diagonal part minus (P o X)(P o X)^T.
// chunk_buf must have room for MultinomialChunkBufDoubleCt() doubles, and
// info_tmp for param_ct^2 doubles when grad is non-null.
static inline uintptr_t MultinomialChunkBufDoubleCt(uint32_t pred_ct, uint32_t class_ct) {
  return (class_ct * (pred_ct + 1) + 1) * S_CAST(uintptr_t, kMultinomialChunkSize);
}

static void MultinomialEval(const double* const* pred_rows, const uint32_t* classes, const double* coefs, uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, double* ln_lik_ptr, double* grad, double* info, double* chunk_buf, double* info_tmp) {
  const uintptr_t param_ct = class_ct * pred_ct;
  // probs: class_ct x kMultinomialChunkSize, holds linear predictors first
  double* probs = chunk_buf;
  double* resid = &(probs[class_ct * kMultinomialChunkSize]);
  // wx[(k,j), i] = p_ik x_ij
  double* wx = &(resid[kMultinomialChunkSize]);
  if (grad) {
    ZeroDArr(param_ct, grad);
    ZeroDArr(param_ct * param_ct, info);
  }
  double ln_lik = 0.0;
  for (uint32_t chunk_start = 0; chunk_start < sample_ct; chunk_start += kMultinomialChunkSize) {
    const uint32_t cur_len = MINV(kMultinomialChunkSize, sample_ct - chunk_start);
    const uint32_t* cur_classes = &(classes[chunk_start]);
    for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
      double* eta = &(probs[class_idx * kMultinomialChunkSize]);
      const double* cur_coefs = &(coefs[class_idx * pred_ct]);
      const double* pred_row = &(pred_rows[0][chunk_start]);
      const double coef0 = cur_coefs[0];
      for (uint32_t uii = 0; uii != cur_len; ++uii) {
        eta[uii] = coef0 * pred_row[uii];
      }
      for (uint32_t pred_idx = 1; pred_idx != pred_ct; ++pred_idx) {
        pred_row = &(pred_rows[pred_idx][chunk_start]);
        const double cur_coef = cur_coefs[pred_idx];
        for (uint32_t uii = 0; uii != cur_len; ++uii) {
          eta[uii] += cur_coef * pred_row[uii];
        }
      }
    }
    for (uint32_t uii = 0; uii != cur_len; ++uii) {
      // the reference class has linear predictor 0
      double max_eta = 0.0;
      for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
        const double cur_eta = probs[class_idx * kMultinomialChunkSize + uii];
        if (cur_eta > max_eta) {
          max_eta = cur_eta;
        }
      }
      const uint32_t cur_class = cur_classes[uii];
      const double obs_eta = cur_class? probs[(cur_class - 1) * kMultinomialChunkSize + uii] : 0.0;
      double denom = exp(-max_eta);
      for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
        double* cur_prob_ptr = &(probs[class_idx * kMultinomialChunkSize + uii]);
        const double cur_exp = exp((*cur_prob_ptr) - max_eta);
        *cur_prob_ptr = cur_exp;
        denom += cur_exp;
      }
      ln_lik += obs_eta - max_eta - log(denom);
      const double denom_recip = 1.0 / denom;
      for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
        probs[class_idx * kMultinomialChunkSize + uii] *= denom_recip;
      }
    }
    if (!grad) {
      continue;
    }
    for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
      const double* cur_probs = &(probs[class_idx * kMultinomialChunkSize]);
      const uint32_t class_code = class_idx + 1;
      for (uint32_t uii = 0; uii != cur_len; ++uii) {
        resid[uii] = S_CAST(double, cur_classes[uii] == class_code) - cur_probs[uii];
      }
      for (uint32_t pred_idx = 0; pred_idx != pred_ct; ++pred_idx) {
        const uintptr_t param_idx = class_idx * pred_ct + pred_idx;
        const double* pred_row = &(pred_rows[pred_idx][chunk_start]);
        grad[param_idx] += DotprodD(resid, pred_row, cur_len);
        double* wx_row = &(wx[param_idx * kMultinomialChunkSize]);
        for (uint32_t uii = 0; uii != cur_len; ++uii) {
          wx_row[uii] = cur_probs[uii] * pred_row[uii];
        }
      }
    }
    MultiplySelfTransposeStrided(wx, param_ct, cur_len, kMultinomialChunkSize, info_tmp);
    for (uintptr_t row_idx = 0; row_idx != param_ct; ++row_idx) {
      double* info_row = &(info[row_idx * param_ct]);
      const double* tmp_row = &(info_tmp[row_idx * param_ct]);
      for (uintptr_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
        info_row[col_idx] -= tmp_row[col_idx];
      }
    }
    for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
      const uintptr_t block_start = class_idx * pred_ct;
      for (uint32_t pred_idx1 = 0; pred_idx1 != pred_ct; ++pred_idx1) {
        const double* wx_row = &(wx[(block_start + pred_idx1) * kMultinomialChunkSize]);
        double* info_row = &(info[(block_start + pred_idx1) * param_ct + block_start]);
        for (uint32_t pred_idx2 = 0; pred_idx2 <= pred_idx1; ++pred_idx2) {
          info_row[pred_idx2] += DotprodD(wx_row, &(pred_rows[pred_idx2][chunk_start]), cur_len);
        }
      }
    }
  }
  *ln_lik_ptr = ln_lik;
}

// Workspace needed by MultinomialFit(), in doubles, beyond chunk_buf and the
// matrix-inversion buffers.
static inline uintptr_t MultinomialFitWkspaceDoubleCt(uintptr_t param_ct) {
  // grad, step, trial, info, info_tmp
  return 3 * param_ct + 2 * param_ct * param_ct;
}

// Newton-Raphson with step halving, starting from coefs[].
// * max_iter == 0: only computes *score_ptr, grad' I^{-1} grad at the
//   starting point.  That is the score statistic when the starting point is
//   the maximum-likelihood estimate under the null hypothesis.
// * Otherwise, coefs[] and *ln_lik_ptr hold the final estimate on return,
//   and *is_unfinished_ptr is set if the iteration limit was reached.  If cov
//   is non-null and the fit converged, cov[] is filled with the inverse of the
//   Fisher information at the estimate (full param_ct x param_ct matrix); if
//   that matrix is singular, cov[0] is set to -1 instead.
// Returns nonzero if the information matrix was singular at an iterate.
// wkspace needs MultinomialFitWkspaceDoubleCt(param_ct) doubles; mi_buf and
// dbl_2d_buf are the usual InvertSymmdefMatrixChecked() buffers.
static BoolErr MultinomialFit(const double* const* pred_rows, const uint32_t* classes, uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, uint32_t max_iter, double* coefs, double* ln_lik_ptr, uint32_t* is_unfinished_ptr, double* score_ptr, double* cov, double* wkspace, double* chunk_buf, MatrixInvertBuf1* mi_buf, double* dbl_2d_buf) {
  const uintptr_t param_ct = class_ct * pred_ct;
  double* grad = wkspace;
  double* step = &(grad[param_ct]);
  double* trial = &(step[param_ct]);
  double* info = &(trial[param_ct]);
  double* info_tmp = &(info[param_ct * param_ct]);
  *is_unfinished_ptr = 0;
  double ln_lik;
  MultinomialEval(pred_rows, classes, coefs, sample_ct, pred_ct, class_ct, &ln_lik, grad, info, chunk_buf, info_tmp);
  for (uint32_t iter_idx = 0; ; ++iter_idx) {
    if (unlikely(!isfinite(ln_lik))) {
      return 1;
    }
    if (InvertSymmdefMatrixChecked(param_ct, info, mi_buf, dbl_2d_buf)) {
      return 1;
    }
    ReflectMatrix(param_ct, info);
    for (uintptr_t param_idx = 0; param_idx != param_ct; ++param_idx) {
      step[param_idx] = DotprodD(&(info[param_idx * param_ct]), grad, param_ct);
    }
    if (!max_iter) {
      *score_ptr = DotprodD(step, grad, param_ct);
      *ln_lik_ptr = ln_lik;
      return 0;
    }
    const double rounding_allowance = fabs(ln_lik) * kMultinomialLnLikRelRounding;
    double step_mult = 1.0;
    double new_ln_lik;
    uint32_t halving_ct = 0;
    while (1) {
      for (uintptr_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        trial[param_idx] = coefs[param_idx] + step_mult * step[param_idx];
      }
      MultinomialEval(pred_rows, classes, trial, sample_ct, pred_ct, class_ct, &new_ln_lik, nullptr, nullptr, chunk_buf, nullptr);
      // (NaN fails this comparison too)
      if (new_ln_lik >= ln_lik - rounding_allowance) {
        break;
      }
      if (++halving_ct == kMultinomialMaxHalvings) {
        // No ascent direction left at double precision: we are at the
        // maximum.
        break;
      }
      step_mult *= 0.5;
    }
    uint32_t converged = (halving_ct == kMultinomialMaxHalvings);
    if (!converged) {
      double max_abs_step = 0.0;
      for (uintptr_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        const double cur_abs_step = fabs(step[param_idx]);
        if (cur_abs_step > max_abs_step) {
          max_abs_step = cur_abs_step;
        }
      }
      max_abs_step *= step_mult;
      converged = (max_abs_step < kMultinomialStepTol) && (new_ln_lik - ln_lik < kMultinomialLnLikTol + rounding_allowance);
      memcpy(coefs, trial, param_ct * sizeof(double));
      ln_lik = new_ln_lik;
    }
    if (converged) {
      break;
    }
    if (iter_idx + 1 == max_iter) {
      *is_unfinished_ptr = 1;
      break;
    }
    MultinomialEval(pred_rows, classes, coefs, sample_ct, pred_ct, class_ct, &ln_lik, grad, info, chunk_buf, info_tmp);
  }
  *ln_lik_ptr = ln_lik;
  if (cov && (!(*is_unfinished_ptr))) {
    double dummy_ln_lik;
    MultinomialEval(pred_rows, classes, coefs, sample_ct, pred_ct, class_ct, &dummy_ln_lik, grad, cov, chunk_buf, info_tmp);
    if (InvertSymmdefMatrixChecked(param_ct, cov, mi_buf, dbl_2d_buf)) {
      cov[0] = -1.0;
    } else {
      ReflectMatrix(param_ct, cov);
    }
  }
  return 0;
}

// Maps the covariate-only fit of a sample set onto a subset of its samples
// in which some levels may be absent.  Both class maps are level ->
// class-index maps; every level present in the subset must be present in the
// full set.  The subset's class 0 is its lowest-indexed level, so when the
// full set's reference class is absent from the subset, the coefficients are
// re-expressed relative to the new reference.
static void RemapNullCoefs(const double* src_coefs, const uint32_t* src_level_to_class, const uint32_t* dst_class_levels, uint32_t dst_class_ct, uint32_t pred_ct, double* dst_coefs) {
  const uint32_t src_ref_class = src_level_to_class[dst_class_levels[0]];
  const double* src_ref_coefs = src_ref_class? (&(src_coefs[(src_ref_class - 1) * pred_ct])) : nullptr;
  for (uint32_t dst_class_idx = 1; dst_class_idx != dst_class_ct; ++dst_class_idx) {
    const uint32_t src_class = src_level_to_class[dst_class_levels[dst_class_idx]];
    double* dst_row = &(dst_coefs[(dst_class_idx - 1) * pred_ct]);
    if (src_class) {
      memcpy(dst_row, &(src_coefs[(src_class - 1) * pred_ct]), pred_ct * sizeof(double));
    } else {
      ZeroDArr(pred_ct, dst_row);
    }
    if (src_ref_coefs) {
      for (uint32_t pred_idx = 0; pred_idx != pred_ct; ++pred_idx) {
        dst_row[pred_idx] -= src_ref_coefs[pred_idx];
      }
    }
  }
}

static inline uintptr_t MultinomialInvertDbl2dCt(uintptr_t dim) {
  return dim * MAXV(dim, 7);
}

BoolErr GlmAllocFillAndTestPhenoCovarsMultinomial(const uintptr_t* sample_include, const PhenoCol* pheno_col, const uint32_t* cat_to_level, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uint32_t level_ct, uintptr_t covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, MultinomialSet* setp, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  *glm_err_ptr = 0;
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uint32_t pred_ct = new_covar_ct + 1;
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  if (unlikely(bigstack_alloc_kcp(new_covar_ct, cur_covar_names_ptr) ||
               bigstack_alloc_u32(sample_ct, &setp->levels) ||
               bigstack_alloc_d(pred_ct * sample_ctav, &setp->covars_pmaj) ||
               bigstack_alloc_u32(level_ct, &setp->null_level_to_class) ||
               bigstack_alloc_d((level_ct - 1) * S_CAST(uintptr_t, pred_ct), &setp->null_coefs))) {
    return 1;
  }
  setp->sample_ct = sample_ct;
  setp->covar_ct = new_covar_ct;
  uint32_t* levels = setp->levels;
  {
    const uint32_t* cats = pheno_col->data.cat;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      levels[sample_idx] = cat_to_level[cats[sample_uidx]];
    }
  }
  double* covars_pmaj = setp->covars_pmaj;
  FillDVec(sample_ct, 1.0, covars_pmaj);
  ZeroDArr(sample_ctav - sample_ct, &(covars_pmaj[sample_ct]));
  unsigned char* bigstack_mark = g_bigstack_base;
  if (new_covar_ct) {
    double* covars_cmaj;
    double* covar_dotprod;
    double* inverse_corr_buf;
    if (unlikely(bigstack_alloc_d(new_covar_ct * sample_ct, &covars_cmaj) ||
                 bigstack_alloc_d(new_covar_ct * new_covar_ct, &covar_dotprod) ||
                 bigstack_alloc_d(new_covar_ct * new_covar_ct, &inverse_corr_buf))) {
      return 1;
    }
    PglErr reterr = GlmFillAndTestCovars(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, 0, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, max_corr, vif_thresh, covar_dotprod, nullptr, inverse_corr_buf, covars_cmaj, *cur_covar_names_ptr, glm_err_ptr);
    if (unlikely(reterr)) {
      return (reterr == kPglRetNomem);
    }
    if (*glm_err_ptr) {
      return 0;
    }
    // Center and scale each covariate.  The likelihood is invariant to this
    // (the model has an intercept, and covariate coefficients are not
    // reported), and it keeps the Newton iterations well-conditioned.
    const double sample_ct_recip = 1.0 / u31tod(sample_ct);
    for (uintptr_t covar_idx = 0; covar_idx != new_covar_ct; ++covar_idx) {
      const double* src = &(covars_cmaj[covar_idx * sample_ct]);
      double* dst = &(covars_pmaj[(covar_idx + 1) * sample_ctav]);
      double sum = 0.0;
      for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        sum += src[sample_idx];
      }
      const double mean = sum * sample_ct_recip;
      double ssq = 0.0;
      for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const double centered = src[sample_idx] - mean;
        dst[sample_idx] = centered;
        ssq += centered * centered;
      }
      // nonconstant covariates guaranteed by GlmDetermineCovars()
      const double scale = 1.0 / sqrt(ssq * sample_ct_recip);
      for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        dst[sample_idx] *= scale;
      }
      ZeroDArr(sample_ctav - sample_ct, &(dst[sample_ct]));
    }
    BigstackReset(bigstack_mark);
  }

  // Covariate-only fit.
  uint32_t* level_to_class = setp->null_level_to_class;
  uint32_t* classes;
  uint32_t* class_sample_cts;
  const double** pred_rows;
  if (unlikely(bigstack_calloc_u32(level_ct, &class_sample_cts) ||
               bigstack_alloc_u32(sample_ct, &classes) ||
               BIGSTACK_ALLOC_X(const double*, pred_ct, &pred_rows))) {
    return 1;
  }
  for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    class_sample_cts[levels[sample_idx]] += 1;
  }
  uint32_t class_ct = 0;
  for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
    if (class_sample_cts[level_idx]) {
      class_sample_cts[class_ct] = class_sample_cts[level_idx];
      level_to_class[level_idx] = class_ct++;
    } else {
      level_to_class[level_idx] = UINT32_MAX;
    }
  }
  // caller has verified class_ct >= 2
  setp->null_class_ct = class_ct;
  for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    classes[sample_idx] = level_to_class[levels[sample_idx]];
  }
  for (uint32_t pred_idx = 0; pred_idx != pred_ct; ++pred_idx) {
    pred_rows[pred_idx] = &(covars_pmaj[pred_idx * sample_ctav]);
  }
  const uint32_t nonref_class_ct = class_ct - 1;
  double* null_coefs = setp->null_coefs;
  // Start from the intercept-only estimate.
  ZeroDArr(nonref_class_ct * S_CAST(uintptr_t, pred_ct), null_coefs);
  const double ref_ct_d = u31tod(class_sample_cts[0]);
  for (uint32_t class_idx = 1; class_idx != class_ct; ++class_idx) {
    null_coefs[(class_idx - 1) * pred_ct] = log(u31tod(class_sample_cts[class_idx]) / ref_ct_d);
  }
  const uintptr_t param_ct = nonref_class_ct * S_CAST(uintptr_t, pred_ct);
  double* fit_wkspace;
  double* chunk_buf;
  MatrixInvertBuf1* mi_buf;
  double* dbl_2d_buf;
  if (unlikely(bigstack_alloc_d(MultinomialFitWkspaceDoubleCt(param_ct), &fit_wkspace) ||
               bigstack_alloc_d(MultinomialChunkBufDoubleCt(pred_ct, nonref_class_ct), &chunk_buf) ||
               BIGSTACK_ALLOC_X(MatrixInvertBuf1, param_ct * kMatrixInvertBuf1CheckedAlloc, &mi_buf) ||
               bigstack_alloc_d(MultinomialInvertDbl2dCt(param_ct), &dbl_2d_buf))) {
    return 1;
  }
  uint32_t is_unfinished;
  if (MultinomialFit(pred_rows, classes, sample_ct, pred_ct, nonref_class_ct, kMultinomialNullMaxIter, null_coefs, &setp->null_ln_lik, &is_unfinished, nullptr, nullptr, fit_wkspace, chunk_buf, mi_buf, dbl_2d_buf) || is_unfinished) {
    *glm_err_ptr = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
  }
  BigstackReset(class_sample_cts);
  return 0;
}

// Per-thread workspace layout, shared by GetMultinomialWorkspaceSize() and
// the compute thread so the two can't drift apart.
typedef struct {
  uintptr_t* sample_nm;
  uint32_t* nm_classes;
  uint32_t* level_sample_cts;
  uint32_t* level_to_class;
  uint32_t* class_levels;
  // per-level dosage range
  double* level_min_geno;
  double* level_max_geno;
  const double** pred_rows;
  // (covar_ct + 1) rows, stride nm_sample_ctav; only used when calls are
  // missing
  double* nm_covars;
  double* geno_vals;
  double* null_coefs;
  double* alt_coefs;
  double* cov;
  double* fit_wkspace;
  double* chunk_buf;
  MatrixInvertBuf1* mi_buf;
  double* dbl_2d_buf;
  // VIF check
  double* dotprods;
  double* inverse_corr_buf;
} MultinomialWorkspace;

// With workspace_buf == nullptr, only returns the size.
static uintptr_t MultinomialWorkspaceLayout(uint32_t sample_ct, uint32_t covar_ct, uint32_t level_ct, unsigned char* workspace_buf, MultinomialWorkspace* wsp) {
  MultinomialWorkspace size_only_ws;
  if (!workspace_buf) {
    wsp = &size_only_ws;
  }
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uint32_t alt_pred_ct = covar_ct + 2;
  const uint32_t nonref_level_ct = level_ct - 1;
  const uintptr_t max_param_ct = nonref_level_ct * S_CAST(uintptr_t, alt_pred_ct);
  const uint32_t nongeno_nonintercept_ct = covar_ct + 1;
  const uintptr_t byte_cts[] = {
    BitCtToWordCt(sample_ct) * sizeof(intptr_t),
    sample_ct * sizeof(int32_t),
    level_ct * sizeof(int32_t),
    level_ct * sizeof(int32_t),
    level_ct * sizeof(int32_t),
    level_ct * sizeof(double),
    level_ct * sizeof(double),
    alt_pred_ct * sizeof(intptr_t),
    (covar_ct + 1) * sample_ctav * sizeof(double),
    sample_ctav * sizeof(double),
    nonref_level_ct * (covar_ct + 1) * sizeof(double),
    max_param_ct * sizeof(double),
    max_param_ct * max_param_ct * sizeof(double),
    MultinomialFitWkspaceDoubleCt(max_param_ct) * sizeof(double),
    MultinomialChunkBufDoubleCt(alt_pred_ct, nonref_level_ct) * sizeof(double),
    max_param_ct * kMatrixInvertBuf1CheckedAlloc,
    MultinomialInvertDbl2dCt(max_param_ct) * sizeof(double),
    nongeno_nonintercept_ct * nongeno_nonintercept_ct * sizeof(double),
    nongeno_nonintercept_ct * nongeno_nonintercept_ct * sizeof(double)
  };
  void** dsts[] = {
    R_CAST(void**, &wsp->sample_nm),
    R_CAST(void**, &wsp->nm_classes),
    R_CAST(void**, &wsp->level_sample_cts),
    R_CAST(void**, &wsp->level_to_class),
    R_CAST(void**, &wsp->class_levels),
    R_CAST(void**, &wsp->level_min_geno),
    R_CAST(void**, &wsp->level_max_geno),
    R_CAST(void**, &wsp->pred_rows),
    R_CAST(void**, &wsp->nm_covars),
    R_CAST(void**, &wsp->geno_vals),
    R_CAST(void**, &wsp->null_coefs),
    R_CAST(void**, &wsp->alt_coefs),
    R_CAST(void**, &wsp->cov),
    R_CAST(void**, &wsp->fit_wkspace),
    R_CAST(void**, &wsp->chunk_buf),
    R_CAST(void**, &wsp->mi_buf),
    R_CAST(void**, &wsp->dbl_2d_buf),
    R_CAST(void**, &wsp->dotprods),
    R_CAST(void**, &wsp->inverse_corr_buf)
  };
  static_assert(sizeof(byte_cts) / sizeof(byte_cts[0]) == sizeof(dsts) / sizeof(dsts[0]), "MultinomialWorkspaceLayout() arrays out of sync.");
  uintptr_t tot_byte_ct = 0;
  for (uint32_t uii = 0; uii != sizeof(byte_cts) / sizeof(byte_cts[0]); ++uii) {
    if (workspace_buf) {
      *(dsts[uii]) = &(workspace_buf[tot_byte_ct]);
    }
    tot_byte_ct += RoundUpPow2(byte_cts[uii], kCacheline);
  }
  return tot_byte_ct;
}

THREAD_FUNC_DECL GlmMultinomialThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GlmMultinomialCtx* ctx = S_CAST(GlmMultinomialCtx*, arg->sharedp->context);
  GlmCtx* common = ctx->common;

  PgenReader* pgrp = common->pgr_ptrs[tidx];
  uintptr_t* genovec = common->genovecs[tidx];
  uintptr_t* dosage_present = nullptr;
  Dosage* dosage_main = nullptr;
  if (common->dosage_presents) {
    dosage_present = common->dosage_presents[tidx];
    dosage_main = common->dosage_mains[tidx];
  }
  unsigned char* workspace_buf = common->workspace_bufs[tidx];
  const uintptr_t* variant_include = common->variant_include;
  const AlleleCode* omitted_alleles = common->omitted_alleles;
  const uintptr_t* sex_male_collapsed = common->sex_male_collapsed;
  const ChrInfo* cip = common->cip;
  const uint32_t* subset_chr_fo_vidx_start = common->subset_chr_fo_vidx_start;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const uint32_t level_ct = ctx->level_ct;
  const uint32_t nonref_level_ct = level_ct - 1;
  const GlmMultinomialTest test_type = ctx->test_type;
  const uint32_t is_score_test = (test_type == kGlmMultinomialTestScore);
  const uint32_t save_coefs = ctx->save_coefs;
  // the Wald test and the coefficient columns need the full fit's inverse
  // information matrix
  const uint32_t need_cov = save_coefs || (test_type == kGlmMultinomialTestWald);
  uint32_t variant_idx_offset = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);

    MultinomialAuxResult* block_aux_iter = &(ctx->block_aux[variant_bidx]);
    double* level_a1_iter = &(ctx->block_level_a1[variant_bidx * S_CAST(uintptr_t, level_ct)]);
    uint32_t* level_allele_obs_iter = &(ctx->block_level_allele_obs[variant_bidx * S_CAST(uintptr_t, level_ct)]);
    double* beta_se_iter = &(ctx->block_beta_se[variant_bidx * 2 * S_CAST(uintptr_t, nonref_level_ct)]);
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = LastLeqU32(subset_chr_fo_vidx_start, 0, cip->chr_ct, variant_idx);
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
      const MultinomialSet* setp;
      if (is_y && common->sample_include_y) {
        cur_sample_include = common->sample_include_y;
        cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
        setp = &(ctx->sets[2]);
      } else if (is_regular_x && common->sample_include_x) {
        cur_sample_include = common->sample_include_x;
        cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
        setp = &(ctx->sets[1]);
      } else {
        cur_sample_include = common->sample_include;
        cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
        setp = &(ctx->sets[0]);
      }
      const uint32_t cur_sample_ct = setp->sample_ct;
      const uint32_t sample_ctl = BitCtToWordCt(cur_sample_ct);
      const uintptr_t sample_ctav = RoundUpPow2(cur_sample_ct, kDoublePerDVec);
      const uint32_t covar_ct = setp->covar_ct;
      const uint32_t null_pred_ct = covar_ct + 1;
      const uint32_t alt_pred_ct = covar_ct + 2;
      const uint32_t* set_levels = setp->levels;
      MultinomialWorkspace ws;
      MultinomialWorkspaceLayout(cur_sample_ct, covar_ct, level_ct, workspace_buf, &ws);
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        uint32_t dosage_ct;
        PglErr reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto GlmMultinomialThread_err;
        }
        ZeroTrailingNyps(cur_sample_ct, genovec);
        GenoarrCountFreqsUnsafe(genovec, cur_sample_ct, genocounts);
        uintptr_t* sample_nm = ws.sample_nm;
        uint32_t missing_ct = genocounts[3];
        if (!missing_ct) {
          SetAllBits(cur_sample_ct, sample_nm);
        } else {
          GenoarrToNonmissing(genovec, cur_sample_ct, sample_nm);
          if (dosage_ct) {
            BitvecOr(dosage_present, sample_ctl, sample_nm);
            missing_ct = cur_sample_ct - PopcountWords(sample_nm, sample_ctl);
          }
        }
        const uint32_t omitted_allele_idx = omitted_alleles? omitted_alleles[variant_uidx] : 0;
        if (omitted_allele_idx) {
          GenovecInvertUnsafe(cur_sample_ct, genovec);
          if (dosage_ct) {
            BiallelicDosage16Invert(dosage_ct, dosage_main);
          }
          const uint32_t uii = genocounts[0];
          genocounts[0] = genocounts[2];
          genocounts[2] = uii;
        }
        const uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
        const uintptr_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kDoublePerDVec);
        // A1 dosages of the samples with a call, in 0..2 units.
        double* geno_vals = ws.geno_vals;
        uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
        uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
        if (!missing_ct) {
          GenoarrLookup16x8bx2(genovec, kSmallDoublePairs, nm_sample_ct, geno_vals);
          if (dosage_ct) {
            uintptr_t sample_idx_base = 0;
            uintptr_t dosage_present_bits = dosage_present[0];
            for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
              const uintptr_t sample_idx = BitIter1(dosage_present, &sample_idx_base, &dosage_present_bits);
              const uint32_t dosage_val = dosage_main[dosage_idx];
              geno_vals[sample_idx] = kRecipDosageMid * u31tod(dosage_val);
              dosage_sum += dosage_val;
              dosage_ssq += dosage_val * dosage_val;
              const uintptr_t cur_geno = GetNyparrEntry(genovec, sample_idx);
              if (cur_geno && (cur_geno != 3)) {
                const uintptr_t prev_val = cur_geno * kDosageMid;
                dosage_sum -= prev_val;
                dosage_ssq -= prev_val * prev_val;
              }
            }
          }
        } else if (!dosage_ct) {
          GenoarrToDoublesRemoveMissing(genovec, kSmallDoubles, cur_sample_ct, geno_vals);
        } else {
          uintptr_t sample_midx_base = 0;
          uintptr_t sample_nm_bits = sample_nm[0];
          uint32_t dosage_idx = 0;
          for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
            const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
            const uintptr_t cur_geno = GetNyparrEntry(genovec, sample_midx);
            double cur_val;
            if (IsSet(dosage_present, sample_midx)) {
              const uint32_t dosage_val = dosage_main[dosage_idx++];
              cur_val = kRecipDosageMid * u31tod(dosage_val);
              dosage_sum += dosage_val;
              dosage_ssq += dosage_val * dosage_val;
              if (cur_geno && (cur_geno != 3)) {
                const uintptr_t prev_val = cur_geno * kDosageMid;
                dosage_sum -= prev_val;
                dosage_ssq -= prev_val * prev_val;
              }
            } else {
              // cur_geno != 3 guaranteed
              cur_val = kSmallDoubles[cur_geno];
            }
            geno_vals[sample_idx] = cur_val;
          }
        }
        {
          uint64_t machr2_dosage_sums[2];
          uint64_t machr2_dosage_ssqs[2];
          machr2_dosage_sums[1 - omitted_allele_idx] = dosage_sum;
          machr2_dosage_ssqs[1 - omitted_allele_idx] = dosage_ssq;
          machr2_dosage_sums[omitted_allele_idx] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[omitted_allele_idx] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
          block_aux_iter->mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, 2);
        }
        // Haploid scaling, per-level counts, and the covariates of the
        // samples with a call.
        uint32_t* level_sample_cts = ws.level_sample_cts;
        ZeroU32Arr(level_ct, level_sample_cts);
        ZeroDArr(level_ct, level_a1_iter);
        ZeroU32Arr(level_ct, level_allele_obs_iter);
        uint32_t* nm_classes = ws.nm_classes;
        double* level_min_geno = ws.level_min_geno;
        double* level_max_geno = ws.level_max_geno;
        {
          uintptr_t sample_idx_base = 0;
          uintptr_t sample_nm_bits = sample_nm[0];
          for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
            const uintptr_t sample_idx = BitIter1(sample_nm, &sample_idx_base, &sample_nm_bits);
            const uint32_t cur_level = set_levels[sample_idx];
            double cur_val = geno_vals[nm_idx];
            uint32_t ploidy = 2;
            if (is_nonx_haploid || (is_regular_x && is_xchr_model_1 && IsSet(sex_male_collapsed, sample_idx))) {
              cur_val *= 0.5;
              ploidy = 1;
              geno_vals[nm_idx] = cur_val;
            }
            if (!level_sample_cts[cur_level]) {
              level_min_geno[cur_level] = cur_val;
              level_max_geno[cur_level] = cur_val;
            } else if (cur_val < level_min_geno[cur_level]) {
              level_min_geno[cur_level] = cur_val;
            } else if (cur_val > level_max_geno[cur_level]) {
              level_max_geno[cur_level] = cur_val;
            }
            // temporarily holds levels
            nm_classes[nm_idx] = cur_level;
            level_sample_cts[cur_level] += 1;
            level_a1_iter[cur_level] += cur_val;
            level_allele_obs_iter[cur_level] += ploidy;
            if (missing_ct) {
              for (uint32_t pred_idx = 0; pred_idx != null_pred_ct; ++pred_idx) {
                ws.nm_covars[pred_idx * nm_sample_ctav + nm_idx] = setp->covars_pmaj[pred_idx * sample_ctav + sample_idx];
              }
            }
          }
        }
        double a1_dosage = 0.0;
        uint32_t allele_obs_ct = 0;
        uint32_t min_level_allele_obs = UINT32_MAX;
        uint32_t class_ct = 0;
        uint32_t* level_to_class = ws.level_to_class;
        uint32_t* class_levels = ws.class_levels;
        for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
          const uint32_t cur_allele_obs = level_allele_obs_iter[level_idx];
          a1_dosage += level_a1_iter[level_idx];
          allele_obs_ct += cur_allele_obs;
          if (level_sample_cts[level_idx]) {
            if (cur_allele_obs < min_level_allele_obs) {
              min_level_allele_obs = cur_allele_obs;
            }
            class_levels[class_ct] = level_idx;
            level_to_class[level_idx] = class_ct++;
          } else {
            level_to_class[level_idx] = UINT32_MAX;
          }
        }
        block_aux_iter->sample_obs_ct = nm_sample_ct;
        block_aux_iter->allele_obs_ct = allele_obs_ct;
        block_aux_iter->a1_dosage = a1_dosage;
        block_aux_iter->min_expected = -9.0;
        if (allele_obs_ct) {
          // smallest expected cell of the (A1, other) x level table under
          // independence
          const double allele_obs_d = u31tod(allele_obs_ct);
          block_aux_iter->min_expected = MINV(a1_dosage, allele_obs_d - a1_dosage) * u31tod(min_level_allele_obs) / allele_obs_d;
        }
        block_aux_iter->chisq = -9.0;
        block_aux_iter->glm_err = 0;
        block_aux_iter->df = class_ct? (class_ct - 1) : 0;
        block_aux_iter->is_unfinished = 0;
        for (uint32_t uii = 0; uii != nonref_level_ct; ++uii) {
          beta_se_iter[2 * uii + 1] = -9.0;
        }
        GlmErr glm_err = 0;
        const uint32_t nonref_class_ct = class_ct - 1;
        const double* const* pred_rows = ws.pred_rows;
        double* null_coefs = ws.null_coefs;
        double* alt_coefs = ws.alt_coefs;
        const uintptr_t null_param_ct = nonref_class_ct * S_CAST(uintptr_t, null_pred_ct);
        const uintptr_t alt_param_ct = nonref_class_ct * S_CAST(uintptr_t, alt_pred_ct);
        double null_ln_lik;
        uint32_t is_unfinished;
        uint32_t coefs_unavailable = 0;
        if (nm_sample_ct <= alt_pred_ct) {
          glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
          goto GlmMultinomialThread_skip_regression;
        }
        if (class_ct < 2) {
          glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
          goto GlmMultinomialThread_skip_regression;
        }
        {
          double min_geno = level_min_geno[class_levels[0]];
          double max_geno = level_max_geno[class_levels[0]];
          for (uint32_t class_idx = 1; class_idx != class_ct; ++class_idx) {
            const uint32_t level_idx = class_levels[class_idx];
            min_geno = MINV(min_geno, level_min_geno[level_idx]);
            max_geno = MAXV(max_geno, level_max_geno[level_idx]);
          }
          if (min_geno == max_geno) {
            glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
            goto GlmMultinomialThread_skip_regression;
          }
        }
        for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
          nm_classes[nm_idx] = level_to_class[nm_classes[nm_idx]];
        }
        {
          const double** pred_rows_w = ws.pred_rows;
          for (uint32_t pred_idx = 0; pred_idx != null_pred_ct; ++pred_idx) {
            pred_rows_w[pred_idx] = missing_ct? (&(ws.nm_covars[pred_idx * nm_sample_ctav])) : (&(setp->covars_pmaj[pred_idx * sample_ctav]));
          }
          pred_rows_w[null_pred_ct] = geno_vals;
        }
        {
          // Center the dosage; this only reparameterizes the intercepts.
          double geno_sum = 0.0;
          for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
            geno_sum += geno_vals[nm_idx];
          }
          const double geno_mean = geno_sum / u31tod(nm_sample_ct);
          for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
            geno_vals[nm_idx] -= geno_mean;
          }
        }
        {
          // Correlation and VIF check on the non-intercept predictors.
          const uint32_t check_pred_ct = alt_pred_ct - 1;
          double* dotprods = ws.dotprods;
          double* row_sums = ws.dbl_2d_buf;
          for (uint32_t pred_idx1 = 0; pred_idx1 != check_pred_ct; ++pred_idx1) {
            const double* row1 = pred_rows[pred_idx1 + 1];
            double* dotprod_row = &(dotprods[pred_idx1 * check_pred_ct]);
            for (uint32_t pred_idx2 = 0; pred_idx2 <= pred_idx1; ++pred_idx2) {
              dotprod_row[pred_idx2] = DotprodD(row1, pred_rows[pred_idx2 + 1], nm_sample_ct);
            }
            double row_sum = 0.0;
            for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
              row_sum += row1[nm_idx];
            }
            row_sums[pred_idx1] = row_sum;
          }
          glm_err = CheckMaxCorrAndVif(dotprods, 0, check_pred_ct, nm_sample_ct, max_corr, vif_thresh, row_sums, nullptr, ws.inverse_corr_buf, ws.mi_buf);
          if (glm_err) {
            goto GlmMultinomialThread_skip_regression;
          }
        }
        {
          // If no dosage in some level exceeds (or falls below) every dosage in
          // the other levels, raising (or lowering) that level's dosage
          // coefficient, with an offsetting intercept, never lowers any
          // sample's likelihood: the maximum-likelihood estimate does not
          // exist.  This includes the common case of a level with no copy of
          // one allele.  The score test does not need the estimate.
          uint32_t separation_allele_idx = UINT32_MAX;
          for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
            double rest_min = DBL_MAX;
            double rest_max = -DBL_MAX;
            for (uint32_t class_idx2 = 0; class_idx2 != class_ct; ++class_idx2) {
              if (class_idx2 != class_idx) {
                const uint32_t level_idx2 = class_levels[class_idx2];
                rest_min = MINV(rest_min, level_min_geno[level_idx2]);
                rest_max = MAXV(rest_max, level_max_geno[level_idx2]);
              }
            }
            const uint32_t level_idx = class_levels[class_idx];
            if (level_max_geno[level_idx] <= rest_min) {
              separation_allele_idx = 1 - omitted_allele_idx;
              break;
            }
            if (level_min_geno[level_idx] >= rest_max) {
              separation_allele_idx = omitted_allele_idx;
              break;
            }
          }
          if (separation_allele_idx != UINT32_MAX) {
            if (!is_score_test) {
              glm_err = SetGlmErr1(kGlmErrcodeSeparation, separation_allele_idx);
              goto GlmMultinomialThread_skip_regression;
            }
            coefs_unavailable = 1;
          }
        }
        // Covariate-only fit to the samples with a call.
        if (!missing_ct) {
          // same samples as the precomputed fit
          memcpy(null_coefs, setp->null_coefs, null_param_ct * sizeof(double));
          null_ln_lik = setp->null_ln_lik;
        } else {
          RemapNullCoefs(setp->null_coefs, setp->null_level_to_class, class_levels, class_ct, null_pred_ct, null_coefs);
          if (MultinomialFit(pred_rows, nm_classes, nm_sample_ct, null_pred_ct, nonref_class_ct, kMultinomialMaxIter, null_coefs, &null_ln_lik, &is_unfinished, nullptr, nullptr, ws.fit_wkspace, ws.chunk_buf, ws.mi_buf, ws.dbl_2d_buf) || is_unfinished) {
            glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
            goto GlmMultinomialThread_skip_regression;
          }
        }
        // Full model, starting from (null estimate, gamma = 0); its first
        // Newton step is the score step.
        for (uint32_t class_idx = 0; class_idx != nonref_class_ct; ++class_idx) {
          double* alt_row = &(alt_coefs[class_idx * alt_pred_ct]);
          memcpy(alt_row, &(null_coefs[class_idx * null_pred_ct]), null_pred_ct * sizeof(double));
          alt_row[null_pred_ct] = 0.0;
        }
        {
          double alt_ln_lik;
          double chisq = 0.0;
          if (is_score_test) {
            double dummy_ln_lik;
            if (MultinomialFit(pred_rows, nm_classes, nm_sample_ct, alt_pred_ct, nonref_class_ct, 0, alt_coefs, &dummy_ln_lik, &is_unfinished, &chisq, nullptr, ws.fit_wkspace, ws.chunk_buf, ws.mi_buf, ws.dbl_2d_buf)) {
              glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
              goto GlmMultinomialThread_skip_regression;
            }
            if (!save_coefs) {
              coefs_unavailable = 1;
            } else if (!coefs_unavailable) {
              if (MultinomialFit(pred_rows, nm_classes, nm_sample_ct, alt_pred_ct, nonref_class_ct, kMultinomialMaxIter, alt_coefs, &alt_ln_lik, &is_unfinished, nullptr, ws.cov, ws.fit_wkspace, ws.chunk_buf, ws.mi_buf, ws.dbl_2d_buf) || is_unfinished) {
                coefs_unavailable = 1;
                block_aux_iter->is_unfinished = 1;
              } else if (ws.cov[0] < 0.0) {
                coefs_unavailable = 1;
              }
            }
          } else {
            if (MultinomialFit(pred_rows, nm_classes, nm_sample_ct, alt_pred_ct, nonref_class_ct, kMultinomialMaxIter, alt_coefs, &alt_ln_lik, &is_unfinished, nullptr, need_cov? ws.cov : nullptr, ws.fit_wkspace, ws.chunk_buf, ws.mi_buf, ws.dbl_2d_buf)) {
              glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
              goto GlmMultinomialThread_skip_regression;
            }
            if (is_unfinished) {
              // Nearly always quasi-separation that the per-level allele
              // counts did not reveal.  The likelihood is still rising, so
              // neither statistic is trustworthy.
              block_aux_iter->is_unfinished = 1;
              coefs_unavailable = 1;
              goto GlmMultinomialThread_skip_regression;
            }
            if (test_type == kGlmMultinomialTestLrt) {
              chisq = 2 * (alt_ln_lik - null_ln_lik);
              if (chisq < 0.0) {
                // rounding
                chisq = 0.0;
              }
            } else {
              // Wald: gamma' Cov(gamma)^{-1} gamma
              if (ws.cov[0] < 0.0) {
                glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                goto GlmMultinomialThread_skip_regression;
              }
              double* wald_cov = ws.fit_wkspace;
              double* wald_coefs = &(wald_cov[nonref_class_ct * nonref_class_ct]);
              const double* cov = ws.cov;
              for (uint32_t class_idx1 = 0; class_idx1 != nonref_class_ct; ++class_idx1) {
                const uintptr_t param_idx1 = class_idx1 * alt_pred_ct + null_pred_ct;
                wald_coefs[class_idx1] = alt_coefs[param_idx1];
                for (uint32_t class_idx2 = 0; class_idx2 != nonref_class_ct; ++class_idx2) {
                  wald_cov[class_idx1 * nonref_class_ct + class_idx2] = cov[param_idx1 * alt_param_ct + class_idx2 * alt_pred_ct + null_pred_ct];
                }
              }
              if (InvertSymmdefMatrixChecked(nonref_class_ct, wald_cov, ws.mi_buf, ws.dbl_2d_buf)) {
                glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                goto GlmMultinomialThread_skip_regression;
              }
              ReflectMatrix(nonref_class_ct, wald_cov);
              for (uint32_t class_idx1 = 0; class_idx1 != nonref_class_ct; ++class_idx1) {
                chisq += wald_coefs[class_idx1] * DotprodD(&(wald_cov[class_idx1 * nonref_class_ct]), wald_coefs, nonref_class_ct);
              }
            }
          }
          if (!isfinite(chisq)) {
            glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
            goto GlmMultinomialThread_skip_regression;
          }
          block_aux_iter->chisq = chisq;
        }
        if (save_coefs && (!coefs_unavailable) && (ws.cov[0] >= 0.0) && (!class_levels[0])) {
          // Coefficients are only comparable across variants when they are
          // relative to the requested reference level.
          const double* cov = ws.cov;
          for (uint32_t class_idx = 1; class_idx != class_ct; ++class_idx) {
            const uintptr_t param_idx = (class_idx - 1) * alt_pred_ct + null_pred_ct;
            double* dst = &(beta_se_iter[2 * (class_levels[class_idx] - 1)]);
            dst[0] = alt_coefs[param_idx];
            dst[1] = sqrt(cov[param_idx * alt_param_ct + param_idx]);
          }
        }
        while (0) {
        GlmMultinomialThread_skip_regression:
          // is_unfinished may be set instead, with glm_err still zero
          memcpy(&(block_aux_iter->glm_err), &glm_err, 8);
          block_aux_iter->chisq = -9.0;
        }
        ++block_aux_iter;
        level_a1_iter = &(level_a1_iter[level_ct]);
        level_allele_obs_iter = &(level_allele_obs_iter[level_ct]);
        beta_se_iter = &(beta_se_iter[2 * nonref_level_ct]);
      }
    }
    variant_idx_offset += cur_block_variant_ct;
    while (0) {
    GlmMultinomialThread_err:
      UpdateU64IfSmaller(new_err_info, &common->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr GlmMultinomial(const char* cur_pheno_name, const char* const* level_names, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const char* outname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, PgenFileInfo* pgfip, GlmMultinomialCtx* ctx, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, uintptr_t* valid_allele_ct_ptr) {
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
    const AlleleCode* omitted_alleles = common->omitted_alleles;
    const uint32_t sample_ct = common->sample_ct;
    const uint32_t sample_ct_x = common->sample_ct_x;
    const uint32_t sample_ct_y = common->sample_ct_y;
    const uint32_t level_ct = ctx->level_ct;
    const uint32_t nonref_level_ct = level_ct - 1;

    uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
    if (max_sample_ct < sample_ct_y) {
      max_sample_ct = sample_ct_y;
    }
    const GlmFlags glm_flags = glm_info_ptr->flags;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    // forced-singlethreaded
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto GlmMultinomial_ret_1;
    }
    const uint32_t report_neglog10p = (glm_flags / kfGlmLog10) & 1;
    const uint32_t true_x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const GlmColFlags glm_cols = glm_info_ptr->cols;
    const uint32_t chr_col = glm_cols & kfGlmColChrom;

    // includes trailing tab
    char* chr_buf = nullptr;
    if (chr_col) {
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto GlmMultinomial_ret_NOMEM;
      }
    }

    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > variant_ct) {
      calc_thread_ct = variant_ct;
    }
    uintptr_t workspace_alloc = MultinomialWorkspaceLayout(sample_ct, ctx->sets[0].covar_ct, level_ct, nullptr, nullptr);
    if (sample_ct_x) {
      const uintptr_t workspace_alloc_x = MultinomialWorkspaceLayout(sample_ct_x, ctx->sets[1].covar_ct, level_ct, nullptr, nullptr);
      if (workspace_alloc_x > workspace_alloc) {
        workspace_alloc = workspace_alloc_x;
      }
    }
    if (sample_ct_y) {
      const uintptr_t workspace_alloc_y = MultinomialWorkspaceLayout(sample_ct_y, ctx->sets[2].covar_ct, level_ct, nullptr, nullptr);
      if (workspace_alloc_y > workspace_alloc) {
        workspace_alloc = workspace_alloc_y;
      }
    }
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    // +1 is for top-level common->workspace_bufs
    const uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;
    const uintptr_t per_variant_xalloc_byte_ct = sizeof(MultinomialAuxResult) + level_ct * (sizeof(double) + sizeof(int32_t)) + 2 * nonref_level_ct * sizeof(double);
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    common->thread_mhc = nullptr;
    common->dosage_presents = nullptr;
    common->dosage_mains = nullptr;
    uint32_t read_block_size;
    uintptr_t max_alt_allele_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, max_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, 0, pgfip, &calc_thread_ct, &common->genovecs, nullptr, nullptr, nullptr, dosage_is_present? (&common->dosage_presents) : nullptr, dosage_is_present? (&common->dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, &max_alt_allele_block_size, main_loadbufs, &common->pgr_ptrs, &common->read_variant_uidx_starts))) {
      goto GlmMultinomial_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto GlmMultinomial_ret_NOMEM;
    }
    MultinomialAuxResult* block_aux_bufs[2];
    double* block_level_a1_bufs[2];
    uint32_t* block_level_allele_obs_bufs[2];
    double* block_beta_se_bufs[2];
    for (uint32_t uii = 0; uii != 2; ++uii) {
      if (unlikely(BIGSTACK_ALLOC_X(MultinomialAuxResult, read_block_size, &(block_aux_bufs[uii])) ||
                   bigstack_alloc_d(read_block_size * S_CAST(uintptr_t, level_ct), &(block_level_a1_bufs[uii])) ||
                   bigstack_alloc_u32(read_block_size * S_CAST(uintptr_t, level_ct), &(block_level_allele_obs_bufs[uii])) ||
                   bigstack_alloc_d(read_block_size * 2 * S_CAST(uintptr_t, nonref_level_ct), &(block_beta_se_bufs[uii])))) {
        goto GlmMultinomial_ret_NOMEM;
      }
    }
    common->workspace_bufs = S_CAST(unsigned char**, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t)));
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      common->workspace_bufs[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(workspace_alloc));
    }
    common->err_info = (~0LLU) << 32;
    SetThreadFuncAndData(GlmMultinomialThread, ctx, &tg);

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
    const uint32_t a1_ct_level_col = glm_cols & kfGlmColA1countcc;
    const uint32_t tot_allele_level_col = glm_cols & kfGlmColTotallelecc;
    const uint32_t a1_freq_col = glm_cols & kfGlmColA1freq;
    const uint32_t a1_freq_level_col = glm_cols & kfGlmColA1freqcc;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t test_col = glm_cols & kfGlmColTest;
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t beta_col = ctx->save_coefs;
    const uint32_t se_col = beta_col && (glm_cols & kfGlmColSe);
    const uint32_t ci_col = beta_col && (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t chisq_col = glm_cols & kfGlmColTz;
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
    if (a1_ct_level_col) {
      for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
        cswritep = strcpya_k(cswritep, "\tA1_CT_");
        cswritep = strcpya(cswritep, level_names[level_idx]);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto GlmMultinomial_ret_WRITE_FAIL;
        }
      }
    }
    if (tot_allele_level_col) {
      for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
        cswritep = strcpya_k(cswritep, "\tALLELE_CT_");
        cswritep = strcpya(cswritep, level_names[level_idx]);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto GlmMultinomial_ret_WRITE_FAIL;
        }
      }
    }
    if (a1_freq_col) {
      cswritep = strcpya_k(cswritep, "\tA1_FREQ");
    }
    if (a1_freq_level_col) {
      for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
        cswritep = strcpya_k(cswritep, "\tA1_FREQ_");
        cswritep = strcpya(cswritep, level_names[level_idx]);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto GlmMultinomial_ret_WRITE_FAIL;
        }
      }
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
    cswritep = strcpya_k(cswritep, "\tMIN_EXPECTED");
    double ci_zt = 0.0;
    if (ci_col) {
      ci_zt = QuantileToZscore((ci_size + 1.0) * 0.5);
    }
    if (beta_col) {
      for (uint32_t level_idx = 1; level_idx != level_ct; ++level_idx) {
        const char* level_name = level_names[level_idx];
        cswritep = strcpya_k(cswritep, "\tBETA_");
        cswritep = strcpya(cswritep, level_name);
        if (se_col) {
          cswritep = strcpya_k(cswritep, "\tSE_");
          cswritep = strcpya(cswritep, level_name);
        }
        if (ci_col) {
          cswritep = strcpya_k(cswritep, "\tL");
          cswritep = dtoa_g(ci_size * 100, cswritep);
          *cswritep++ = '_';
          cswritep = strcpya(cswritep, level_name);
          cswritep = strcpya_k(cswritep, "\tU");
          cswritep = dtoa_g(ci_size * 100, cswritep);
          *cswritep++ = '_';
          cswritep = strcpya(cswritep, level_name);
        }
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto GlmMultinomial_ret_WRITE_FAIL;
        }
      }
    }
    if (chisq_col) {
      cswritep = strcpya_k(cswritep, "\tCHISQ\tDF");
    }
    if (p_col) {
      if (report_neglog10p) {
        cswritep = strcpya_k(cswritep, "\tNEG_LOG10_P");
      } else {
        cswritep = strcpya_k(cswritep, "\tP");
      }
    }
    if (err_col) {
      cswritep = strcpya_k(cswritep, "\tERRCODE");
    }
    AppendBinaryEoln(&cswritep);

    // Same block pipeline as GlmLogistic().
    uintptr_t write_variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t suppress_mach_r2 = 0;
    uint32_t prev_block_variant_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    uintptr_t valid_allele_ct = 0;
    const char* test_type_str = "LRT";
    if (ctx->test_type == kGlmMultinomialTestScore) {
      test_type_str = "score test";
    } else if (ctx->test_type == kGlmMultinomialTestWald) {
      test_type_str = "Wald test";
    }
    logprintfww5("--glm multinomial logistic regression (%s, %u levels) on phenotype '%s': ", test_type_str, level_ct, cur_pheno_name);
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_variant_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto GlmMultinomial_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, common->err_info);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, common->err_info >> 32);
          goto GlmMultinomial_ret_1;
        }
      }
      if (!IsLastBlock(&tg)) {
        common->cur_block_variant_ct = cur_block_variant_ct;
        const uint32_t uidx_start = read_block_idx * read_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_variant_ct, calc_thread_ct, uidx_start, common->read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, common->pgr_ptrs);
        ctx->block_aux = block_aux_bufs[parity];
        ctx->block_level_a1 = block_level_a1_bufs[parity];
        ctx->block_level_allele_obs = block_level_allele_obs_bufs[parity];
        ctx->block_beta_se = block_beta_se_bufs[parity];
        if (variant_idx + cur_block_variant_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto GlmMultinomial_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const MultinomialAuxResult* auxp = block_aux_bufs[parity];
        const double* level_a1_iter = block_level_a1_bufs[parity];
        const uint32_t* level_allele_obs_iter = block_level_allele_obs_bufs[parity];
        const double* beta_se_iter = block_beta_se_bufs[parity];
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            suppress_mach_r2 = (chr_idx == true_x_code) || (chr_idx == mt_code);
            if (chr_col) {
              char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
              *chr_name_end = '\t';
              chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
            }
          }
          const uintptr_t allele_idx_offset_base = common->allele_idx_offsets? common->allele_idx_offsets[write_variant_uidx] : (2 * S_CAST(uintptr_t, write_variant_uidx));
          const uint32_t omitted_allele_idx = omitted_alleles? omitted_alleles[write_variant_uidx] : 0;
          const uint32_t a1_allele_idx = 1 - omitted_allele_idx;
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          const uint32_t is_valid = (auxp->chisq != -9.0);
          double ln_pval = kLnPvalError;
          if (is_valid) {
            ln_pval = ChisqToLnP(auxp->chisq, auxp->df);
            if (orig_ln_pvals) {
              orig_ln_pvals[valid_allele_ct] = ln_pval;
            }
            ++valid_allele_ct;
            if (valid_alleles) {
              SetBit(allele_idx_offset_base + a1_allele_idx, valid_alleles);
            }
          } else if (valid_alleles) {
            ClearBit(write_variant_uidx, valid_variants);
          }
          if ((ln_pfilter <= 0.0) && ((!is_valid) || (ln_pval > ln_pfilter))) {
            goto GlmMultinomial_variant_iterate;
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
            cswritep = strcpya(cswritep, cur_alleles[1]);
          }
          *cswritep++ = '\t';
          if (provref_col) {
            *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, write_variant_uidx)))? 'Y' : 'N';
            *cswritep++ = '\t';
          }
          cswritep = strcpya(cswritep, cur_alleles[a1_allele_idx]);
          if (omitted_col) {
            *cswritep++ = '\t';
            cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
          }
          if (ax_col) {
            *cswritep++ = '\t';
            cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
          }
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto GlmMultinomial_ret_WRITE_FAIL;
          }
          if (a1_ct_col) {
            *cswritep++ = '\t';
            cswritep = dtoa_g(auxp->a1_dosage, cswritep);
          }
          if (tot_allele_col) {
            *cswritep++ = '\t';
            cswritep = u32toa(auxp->allele_obs_ct, cswritep);
          }
          if (a1_ct_level_col) {
            for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
              *cswritep++ = '\t';
              cswritep = dtoa_g(level_a1_iter[level_idx], cswritep);
            }
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto GlmMultinomial_ret_WRITE_FAIL;
            }
          }
          if (tot_allele_level_col) {
            for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
              *cswritep++ = '\t';
              cswritep = u32toa(level_allele_obs_iter[level_idx], cswritep);
            }
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto GlmMultinomial_ret_WRITE_FAIL;
            }
          }
          if (a1_freq_col) {
            *cswritep++ = '\t';
            if (auxp->allele_obs_ct) {
              cswritep = dtoa_g(auxp->a1_dosage / u31tod(auxp->allele_obs_ct), cswritep);
            } else {
              cswritep = strcpya_k(cswritep, "NA");
            }
          }
          if (a1_freq_level_col) {
            for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
              *cswritep++ = '\t';
              const uint32_t cur_allele_obs = level_allele_obs_iter[level_idx];
              if (cur_allele_obs) {
                cswritep = dtoa_g(level_a1_iter[level_idx] / u31tod(cur_allele_obs), cswritep);
              } else {
                cswritep = strcpya_k(cswritep, "NA");
              }
            }
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto GlmMultinomial_ret_WRITE_FAIL;
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
            cswritep = strcpya_k(cswritep, "\tADD");
          }
          if (nobs_col) {
            *cswritep++ = '\t';
            cswritep = u32toa(auxp->sample_obs_ct, cswritep);
          }
          *cswritep++ = '\t';
          if (auxp->min_expected != -9.0) {
            cswritep = dtoa_g(auxp->min_expected, cswritep);
          } else {
            cswritep = strcpya_k(cswritep, "NA");
          }
          if (beta_col) {
            for (uint32_t nonref_level_idx = 0; nonref_level_idx != nonref_level_ct; ++nonref_level_idx) {
              const double beta = beta_se_iter[2 * nonref_level_idx];
              const double se = beta_se_iter[2 * nonref_level_idx + 1];
              const uint32_t coef_is_valid = (se != -9.0);
              *cswritep++ = '\t';
              if (coef_is_valid) {
                cswritep = dtoa_g(beta, cswritep);
              } else {
                cswritep = strcpya_k(cswritep, "NA");
              }
              if (se_col) {
                *cswritep++ = '\t';
                if (coef_is_valid) {
                  cswritep = dtoa_g(se, cswritep);
                } else {
                  cswritep = strcpya_k(cswritep, "NA");
                }
              }
              if (ci_col) {
                *cswritep++ = '\t';
                if (coef_is_valid) {
                  const double ci_radius = ci_zt * se;
                  cswritep = dtoa_g(beta - ci_radius, cswritep);
                  *cswritep++ = '\t';
                  cswritep = dtoa_g(beta + ci_radius, cswritep);
                } else {
                  cswritep = strcpya_k(cswritep, "NA\tNA");
                }
              }
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto GlmMultinomial_ret_WRITE_FAIL;
              }
            }
          }
          if (chisq_col) {
            *cswritep++ = '\t';
            if (is_valid) {
              cswritep = dtoa_g(auxp->chisq, cswritep);
              *cswritep++ = '\t';
              cswritep = u32toa(auxp->df, cswritep);
            } else {
              cswritep = strcpya_k(cswritep, "NA\tNA");
            }
          }
          if (p_col) {
            *cswritep++ = '\t';
            if (is_valid) {
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
            if (auxp->is_unfinished) {
              cswritep = strcpya_k(cswritep, "UNFINISHED");
            } else if (is_valid) {
              *cswritep++ = '.';
            } else {
              uint64_t glm_errcode;
              memcpy(&glm_errcode, &(auxp->glm_err), 8);
              cswritep = AppendGlmErrstr(glm_errcode, cswritep);
            }
          }
          AppendBinaryEoln(&cswritep);
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto GlmMultinomial_ret_WRITE_FAIL;
          }
        GlmMultinomial_variant_iterate:
          ++auxp;
          level_a1_iter = &(level_a1_iter[level_ct]);
          level_allele_obs_iter = &(level_allele_obs_iter[level_ct]);
          beta_se_iter = &(beta_se_iter[2 * nonref_level_ct]);
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
      }
      ++read_block_idx;
      prev_block_variant_ct = cur_block_variant_ct;
      variant_idx += cur_block_variant_ct;
      // crucially, this is independent of the PgenReader block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto GlmMultinomial_ret_WRITE_FAIL;
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
  GlmMultinomial_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GlmMultinomial_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  GlmMultinomial_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GlmMultinomial_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 GlmMultinomial_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
