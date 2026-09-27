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
#include "include/plink2_fmath.h"
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

// Firth-penalized fit (Firth 1993; multinomial form of Bull, Mak and
// Greenwood 2002).  The basic iteration (a quasi-Newton step on the modified
// score, see MultinomialFirthFit()) converges linearly, and can be slow near
// a (quasi-)separated solution, so once the step sizes show slow linear
// convergence it is replaced by BFGS updates.
static const uint32_t kMultinomialFirthMaxIter = 100;
static const double kMultinomialFirthSlowRatio = 0.25;

// Scratch space for MultinomialFirthFit().  With p = pred_ct, J = class_ct
// (non-reference classes), P = J * p, and n rounded up to a vector boundary:
//   prob, ww: J * n doubles each
//   max_eta, denom, sqrt_wts: n doubles each
//   zz: (P + p) * kMultinomialChunkSize doubles
//   ss: (P + p)^2 doubles
//   info, chol, iinv, linv, info_aug, bfgs: P^2 doubles each
//   xty, grad, step, coef_old, ustar, ustar_f, ustar_old, step_f, sdiff,
//     tmpv: P doubles each
//   hab: kMultinomialChunkSize doubles
//   xtw: p doubles
typedef struct {
  double* prob;
  double* max_eta;
  double* denom;
  double* zz;
  double* ss;
  double* info;
  double* chol;
  double* iinv;
  double* linv;
  double* info_aug;
  double* bfgs;
  double* ww;
  double* sqrt_wts;
  double* hab;
  double* xty;
  double* grad;
  double* step;
  double* coef_old;
  double* ustar;
  double* ustar_f;
  double* ustar_old;
  double* step_f;
  double* sdiff;
  double* tmpv;
  double* xtw;
} MultinomialFirthBufs;

// With buf == nullptr, only returns the size.
static uintptr_t MultinomialFirthBufsLayout(uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, unsigned char* buf, MultinomialFirthBufs* fbp) {
  MultinomialFirthBufs size_only_fb;
  if (!buf) {
    fbp = &size_only_fb;
  }
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uintptr_t param_ct = class_ct * S_CAST(uintptr_t, pred_ct);
  const uintptr_t stacked_ct = param_ct + pred_ct;
  const uintptr_t dbl_cts[] = {
    class_ct * sample_ctav,
    sample_ctav,
    sample_ctav,
    stacked_ct * kMultinomialChunkSize,
    stacked_ct * stacked_ct,
    param_ct * param_ct,
    param_ct * param_ct,
    param_ct * param_ct,
    param_ct * param_ct,
    param_ct * param_ct,
    param_ct * param_ct,
    class_ct * sample_ctav,
    sample_ctav,
    kMultinomialChunkSize,
    param_ct,
    param_ct,
    param_ct,
    param_ct,
    param_ct,
    param_ct,
    param_ct,
    param_ct,
    param_ct,
    param_ct,
    pred_ct
  };
  double** dsts[] = {
    &fbp->prob,
    &fbp->max_eta,
    &fbp->denom,
    &fbp->zz,
    &fbp->ss,
    &fbp->info,
    &fbp->chol,
    &fbp->iinv,
    &fbp->linv,
    &fbp->info_aug,
    &fbp->bfgs,
    &fbp->ww,
    &fbp->sqrt_wts,
    &fbp->hab,
    &fbp->xty,
    &fbp->grad,
    &fbp->step,
    &fbp->coef_old,
    &fbp->ustar,
    &fbp->ustar_f,
    &fbp->ustar_old,
    &fbp->step_f,
    &fbp->sdiff,
    &fbp->tmpv,
    &fbp->xtw
  };
  static_assert(sizeof(dbl_cts) / sizeof(dbl_cts[0]) == sizeof(dsts) / sizeof(dsts[0]), "MultinomialFirthBufsLayout() arrays out of sync.");
  uintptr_t tot_byte_ct = 0;
  for (uint32_t uii = 0; uii != sizeof(dbl_cts) / sizeof(dbl_cts[0]); ++uii) {
    if (buf) {
      *(dsts[uii]) = R_CAST(double*, &(buf[tot_byte_ct]));
    }
    tot_byte_ct += RoundUpPow2(dbl_cts[uii] * sizeof(double), kCacheline);
  }
  return tot_byte_ct;
}

// Fills prob[] (class_ct rows of sample_ctav) with the fitted class
// probabilities under coefs[], and lli[] with each sample's log-likelihood
// contribution.  Returns the total log-likelihood, which is not finite when
// the linear predictors overflow.
// xx is predictor-major with vector-aligned rows and zeroed trailing elements;
// coefs is class-major.
static double MultinomialFirthProbs(const double* xx, const double* coefs, const uint32_t* classes, uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, MultinomialFirthBufs* fbp, double* lli) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uintptr_t sample_ct_rem = sample_ctav - sample_ct;
  double* prob = fbp->prob;
  double* max_eta = fbp->max_eta;
  double* denom = fbp->denom;
  // Linear predictors.  The reference class's is zero; each sample's are
  // shifted by their maximum (or zero) before exponentiating.
  ZeroDArr(sample_ctav, max_eta);
  for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
    double* cur_eta = &(prob[class_idx * sample_ctav]);
    ColMajorMatrixVectorMultiplyStrided(xx, &(coefs[class_idx * pred_ct]), sample_ct, sample_ctav, pred_ct, cur_eta);
    ZeroDArr(sample_ct_rem, &(cur_eta[sample_ct]));
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      if (cur_eta[sample_idx] > max_eta[sample_idx]) {
        max_eta[sample_idx] = cur_eta[sample_idx];
      }
    }
  }
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uint32_t cur_class = classes[sample_idx];
    const double cur_eta = cur_class? prob[(cur_class - 1) * sample_ctav + sample_idx] : 0.0;
    lli[sample_idx] = cur_eta - max_eta[sample_idx];
    denom[sample_idx] = -max_eta[sample_idx];
  }
  ZeroDArr(sample_ct_rem, &(denom[sample_ct]));
  expd_v(denom, sample_ctav);
  for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
    double* cur_prob = &(prob[class_idx * sample_ctav]);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      cur_prob[sample_idx] -= max_eta[sample_idx];
    }
    expd_v(cur_prob, sample_ctav);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      denom[sample_idx] += cur_prob[sample_idx];
    }
  }
  double ln_lik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const double cur_lli = lli[sample_idx] - log(denom[sample_idx]);
    lli[sample_idx] = cur_lli;
    ln_lik += cur_lli;
    // denom[] becomes its reciprocal
    denom[sample_idx] = 1.0 / denom[sample_idx];
  }
  for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
    double* cur_prob = &(prob[class_idx * sample_ctav]);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      cur_prob[sample_idx] *= denom[sample_idx];
    }
  }
  return ln_lik;
}

// Lower triangle of gg (dim x dim, row-major) := beta * gg + zz zz^T, where zz
// has dim rows of col_ct elements at the given stride.
static void MultinomialSyrk(const double* zz, uint32_t dim, uint32_t col_ct, uint32_t stride, double beta, double* gg) {
#ifdef NOLAPACK
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    const double* zz_row = &(zz[row_idx * stride]);
    double* gg_row = &(gg[row_idx * dim]);
    for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
      const double dotprod = DotprodD(zz_row, &(zz[col_idx * stride]), col_ct);
      // Like cblas_dsyrk(), don't read gg when beta is zero: it may be
      // uninitialized.
      gg_row[col_idx] = (beta == 0.0)? dotprod : (beta * gg_row[col_idx] + dotprod);
    }
  }
#else
  cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, dim, col_ct, 1.0, zz, stride, beta, gg, dim);
#endif
}

// cc (row_ct x col_ct, column-major with stride ldc) := aa bb, where aa is
// row_ct x common_ct (column-major, stride lda) and bb is common_ct x col_ct
// (column-major, stride ldb).  Never reads cc.
static void MultinomialGemm(const double* aa, const double* bb, uint32_t row_ct, uint32_t lda, uint32_t col_ct, uint32_t ldb, uint32_t common_ct, uint32_t ldc, double* cc) {
#ifdef NOLAPACK
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    double* cc_col = &(cc[col_idx * ldc]);
    ZeroDArr(row_ct, cc_col);
    const double* bb_col = &(bb[col_idx * ldb]);
    for (uint32_t com_idx = 0; com_idx != common_ct; ++com_idx) {
      const double mult = bb_col[com_idx];
      const double* aa_col = &(aa[com_idx * lda]);
      for (uint32_t row_idx = 0; row_idx != row_ct; ++row_idx) {
        cc_col[row_idx] += mult * aa_col[row_idx];
      }
    }
  }
#else
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row_ct, col_ct, common_ct, 1.0, aa, lda, bb, ldb, 0.0, cc, ldc);
#endif
}

// Computes the log-likelihood gradient and the lower triangle of the
// information matrix at the probabilities in fbp->prob.
// If sqrt_wts is not nullptr, sample i's contribution to the information
// matrix is weighted by sqrt_wts[i]^2, and the gradient is not computed (grad
// may then be nullptr).
// The information matrix's (c, d) block, pred_ct square, is
//   X^T diag(P_c (1[c = d] - P_d)) X  =  1[c = d] D_c - Z_c^T Z_d,
// where Z_c = diag(P_c) X and D_c = Z_c^T X.  A chunk of samples at a time, X
// and all the Z_c are stacked into one matrix, whose self-product (one syrk)
// contains every D_c and every Z_c^T Z_d.
// The gradient's block c is
//   X^T (Y_c - P_c)  =  xty_c - (column 0 of D_c),
// since predictor 0 is the intercept.  fbp->xty must hold X^T Y.
static void MultinomialFirthGradAndInfo(const double* xx, const double* sqrt_wts, uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, MultinomialFirthBufs* fbp, double* grad, double* info) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uint32_t param_ct = pred_ct * class_ct;
  const uint32_t stacked_ct = param_ct + pred_ct;
  const double* prob = fbp->prob;
  double* zz = fbp->zz;
  double* ss = fbp->ss;
  for (uint32_t chunk_start = 0; chunk_start < sample_ct; chunk_start += kMultinomialChunkSize) {
    const uint32_t cur_len = MINV(kMultinomialChunkSize, sample_ct - chunk_start);
    double* zz_row = zz;
    for (uint32_t pred_idx = 0; pred_idx != pred_ct; ++pred_idx) {
      const double* xx_row = &(xx[pred_idx * sample_ctav + chunk_start]);
      if (!sqrt_wts) {
        memcpy(zz_row, xx_row, cur_len * sizeof(double));
      } else {
        const double* cur_sqrt_wts = &(sqrt_wts[chunk_start]);
        for (uint32_t uii = 0; uii != cur_len; ++uii) {
          zz_row[uii] = cur_sqrt_wts[uii] * xx_row[uii];
        }
      }
      zz_row = &(zz_row[kMultinomialChunkSize]);
    }
    for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
      const double* cur_prob = &(prob[class_idx * sample_ctav + chunk_start]);
      for (uint32_t pred_idx = 0; pred_idx != pred_ct; ++pred_idx) {
        const double* xx_row = &(zz[pred_idx * kMultinomialChunkSize]);
        for (uint32_t uii = 0; uii != cur_len; ++uii) {
          zz_row[uii] = cur_prob[uii] * xx_row[uii];
        }
        zz_row = &(zz_row[kMultinomialChunkSize]);
      }
    }
    MultinomialSyrk(zz, stacked_ct, cur_len, kMultinomialChunkSize, chunk_start? 1.0 : 0.0, ss);
  }
  // ss (row-major lower triangle): row pred_ct + r holds D in its first
  // pred_ct columns and Z^T Z after that.
  for (uint32_t row_idx = 0; row_idx != param_ct; ++row_idx) {
    const double* ss_row = &(ss[(pred_ct + row_idx) * stacked_ct]);
    const double* gg_row = &(ss_row[pred_ct]);
    if (!sqrt_wts) {
      grad[row_idx] = fbp->xty[row_idx] - ss_row[0];
    }
    double* info_row = &(info[row_idx * param_ct]);
    const uint32_t block_start_col = row_idx - (row_idx % pred_ct);
    for (uint32_t col_idx = 0; col_idx != block_start_col; ++col_idx) {
      info_row[col_idx] = -gg_row[col_idx];
    }
    for (uint32_t col_idx = block_start_col; col_idx <= row_idx; ++col_idx) {
      info_row[col_idx] = ss_row[col_idx - block_start_col] - gg_row[col_idx];
    }
  }
}

// Lower-triangular Cholesky factor of the dim x dim symmetric matrix aa
// (row-major; only the lower triangle is read).  Returns 1 if aa is not
// numerically positive definite.
static BoolErr MultinomialCholesky(const double* aa, uint32_t dim, double* ll) {
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    const double* aa_row = &(aa[row_idx * dim]);
    double* ll_row = &(ll[row_idx * dim]);
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      const double* ll_row2 = &(ll[col_idx * dim]);
      double dxx = aa_row[col_idx];
      for (uint32_t uii = 0; uii != col_idx; ++uii) {
        dxx -= ll_row[uii] * ll_row2[uii];
      }
      ll_row[col_idx] = dxx / ll_row2[col_idx];
    }
    double diag = aa_row[row_idx];
    for (uint32_t uii = 0; uii != row_idx; ++uii) {
      diag -= ll_row[uii] * ll_row[uii];
    }
    // also catches nan
    if (!(diag > kMatrixSingularRcond * aa_row[row_idx])) {
      return 1;
    }
    ll_row[row_idx] = sqrt(diag);
  }
  return 0;
}

// Solves (L L^T) xx = yy.
static void MultinomialCholSolve(const double* ll, const double* yy, uint32_t dim, double* __restrict xx) {
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    const double* ll_row = &(ll[row_idx * dim]);
    double dxx = yy[row_idx];
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      dxx -= ll_row[col_idx] * xx[col_idx];
    }
    xx[row_idx] = dxx / ll_row[row_idx];
  }
  for (uint32_t row_idx = dim; row_idx; ) {
    --row_idx;
    double dxx = xx[row_idx];
    for (uint32_t row_idx2 = row_idx + 1; row_idx2 != dim; ++row_idx2) {
      dxx -= ll[row_idx2 * dim + row_idx] * xx[row_idx2];
    }
    xx[row_idx] = dxx / ll[row_idx * dim + row_idx];
  }
}

// Full inverse (both triangles, row-major) of L L^T, via L^{-1}: the inverse
// is L^{-T} L^{-1}.  linv is scratch space (dim x dim).
static void MultinomialCholInvFull(const double* ll, uint32_t dim, double* __restrict linv, double* __restrict inv) {
  // Row-major lower triangle of L^{-1}, one column at a time.
  for (uint32_t col_idx = 0; col_idx != dim; ++col_idx) {
    for (uint32_t row_idx = col_idx; row_idx != dim; ++row_idx) {
      const double* ll_row = &(ll[row_idx * dim]);
      double dxx = (row_idx == col_idx)? 1.0 : 0.0;
      for (uint32_t uii = col_idx; uii != row_idx; ++uii) {
        dxx -= ll_row[uii] * linv[uii * dim + col_idx];
      }
      linv[row_idx * dim + col_idx] = dxx / ll_row[row_idx];
    }
  }
  // inv[i][j] = sum_{k >= max(i, j)} linv[k][i] * linv[k][j]
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
      double dxx = 0.0;
      for (uint32_t kk = row_idx; kk != dim; ++kk) {
        const double* linv_row = &(linv[kk * dim]);
        dxx += linv_row[row_idx] * linv_row[col_idx];
      }
      inv[row_idx * dim + col_idx] = dxx;
      inv[col_idx * dim + row_idx] = dxx;
    }
  }
}

// Evaluates everything the penalized fit needs at coefs[]: fbp->prob and
// lli[] (see MultinomialFirthProbs()), the gradient and information matrix
// (see MultinomialFirthGradAndInfo()), and the information matrix's Cholesky
// factor.  Sets *ln_lik_ptr and *logdet_ptr (log det of the information
// matrix).  Returns 1 if the log-likelihood is not finite or the information
// matrix is not numerically positive definite.
static BoolErr MultinomialFirthEval(const double* xx, const double* coefs, const uint32_t* classes, uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, MultinomialFirthBufs* fbp, double* lli, double* ln_lik_ptr, double* logdet_ptr) {
  const uint32_t param_ct = pred_ct * class_ct;
  const double ln_lik = MultinomialFirthProbs(xx, coefs, classes, sample_ct, pred_ct, class_ct, fbp, lli);
  if (!isfinite(ln_lik)) {
    return 1;
  }
  MultinomialFirthGradAndInfo(xx, nullptr, sample_ct, pred_ct, class_ct, fbp, fbp->grad, fbp->info);
  if (MultinomialCholesky(fbp->info, param_ct, fbp->chol)) {
    return 1;
  }
  double half_logdet = 0.0;
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
    half_logdet += log(fbp->chol[param_idx * (param_ct + 1)]);
  }
  *ln_lik_ptr = ln_lik;
  *logdet_ptr = 2 * half_logdet;
  return 0;
}

// Firth's modified score at the point last evaluated by
// MultinomialFirthEval():
//   U*_(c,j) = U_(c,j) + 0.5 tr(I^{-1} dI/db_(c,j)).
// With W_i = diag(p_i) - p_i p_i^T and the per-sample "hat" blocks
//   h_i^(ab) = x_i^T (I^{-1})_(ab) x_i,
// the trace is sum_i x_ij sum_(a,b) [dW_i/deta_ic]_(ab) h_i^(ab), which works
// out to sum_i x_ij p_ic (r_ic - sum_a p_ia r_ia), where
//   r_ia = h_i^(aa) - 2 sum_b p_ib h_i^(ab).
// (With two classes, 0.5 p (r - p r) = v h (0.5 - p), logistf's
// correction.)  The h_i^(ab) are computed a chunk of samples and a class a at
// a time: G = X (I^{-1})_(a, b >= a) is one matrix product, and
// h_i^(ab) = sum_l G[i][(b - a) * p + l] x_il (h_i^(ba) = h_i^(ab)).
//
// Also leaves sqrt(1 + tr(W_i h_i)) in fbp->sqrt_wts, for
// MultinomialFirthAugInfo(), and the full inverse of the information matrix
// in fbp->iinv.
static void MultinomialFirthScore(const double* xx, uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, MultinomialFirthBufs* fbp) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uint32_t param_ct = pred_ct * class_ct;
  double* iinv = fbp->iinv;
  MultinomialCholInvFull(fbp->chol, param_ct, fbp->linv, iinv);
  const double* prob = fbp->prob;
  double* ww = fbp->ww;
  double* sqrt_wts = fbp->sqrt_wts;
  double* hab = fbp->hab;
  // G, cur_len x param_ct with stride kMultinomialChunkSize, fits in zz.
  double* gg = fbp->zz;
  for (uint32_t chunk_start = 0; chunk_start < sample_ct; chunk_start += kMultinomialChunkSize) {
    const uint32_t cur_len = MINV(kMultinomialChunkSize, sample_ct - chunk_start);
    const double* xx_chunk = &(xx[chunk_start]);
    // tr(W_i h_i) = sum_(a,b) p_a (1[a = b] - p_b) h^(ab), accumulated here
    // before becoming sqrt(1 + tr(W_i h_i))
    double* cur_sqrt_wts = &(sqrt_wts[chunk_start]);
    ZeroDArr(cur_len, cur_sqrt_wts);
    for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
      ZeroDArr(cur_len, &(ww[class_idx * sample_ctav + chunk_start]));
    }
    for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
      const double* prob_a = &(prob[class_idx * sample_ctav + chunk_start]);
      double* rr_a = &(ww[class_idx * sample_ctav + chunk_start]);
      // Columns b >= a of (I^{-1})_(a, *), i.e. pred_ct rows of iinv from
      // column a * pred_ct on, are (by symmetry) also the transpose of the
      // corresponding columns: a column-major pred_ct x ((J - a) * pred_ct)
      // matrix with stride param_ct.
      const uint32_t col_start = class_idx * pred_ct;
      MultinomialGemm(xx_chunk, &(iinv[col_start * param_ct + col_start]), cur_len, sample_ctav, param_ct - col_start, param_ct, pred_ct, kMultinomialChunkSize, gg);
      for (uint32_t class_idx2 = class_idx; class_idx2 != class_ct; ++class_idx2) {
        const double* gg_iter = &(gg[(class_idx2 - class_idx) * pred_ct * kMultinomialChunkSize]);
        for (uint32_t uii = 0; uii != cur_len; ++uii) {
          hab[uii] = gg_iter[uii] * xx_chunk[uii];
        }
        for (uint32_t pred_idx = 1; pred_idx != pred_ct; ++pred_idx) {
          const double* gg_col = &(gg_iter[pred_idx * kMultinomialChunkSize]);
          const double* xx_col = &(xx_chunk[pred_idx * sample_ctav]);
          for (uint32_t uii = 0; uii != cur_len; ++uii) {
            hab[uii] += gg_col[uii] * xx_col[uii];
          }
        }
        const double* prob_b = &(prob[class_idx2 * sample_ctav + chunk_start]);
        if (class_idx2 == class_idx) {
          for (uint32_t uii = 0; uii != cur_len; ++uii) {
            rr_a[uii] += hab[uii] * (1.0 - 2 * prob_a[uii]);
            cur_sqrt_wts[uii] += hab[uii] * prob_a[uii] * (1.0 - prob_a[uii]);
          }
        } else {
          // h^(ab) and h^(ba) both
          double* rr_b = &(ww[class_idx2 * sample_ctav + chunk_start]);
          for (uint32_t uii = 0; uii != cur_len; ++uii) {
            const double twice_hab = 2 * hab[uii];
            rr_a[uii] -= twice_hab * prob_b[uii];
            rr_b[uii] -= twice_hab * prob_a[uii];
            cur_sqrt_wts[uii] -= twice_hab * prob_a[uii] * prob_b[uii];
          }
        }
      }
    }
    for (uint32_t uii = 0; uii != cur_len; ++uii) {
      // (the trace is nonnegative; this only guards against rounding)
      cur_sqrt_wts[uii] = sqrt(1.0 + MAXV(cur_sqrt_wts[uii], 0.0));
    }
    // ww_c := 0.5 p_c (r_c - sum_a p_a r_a); hab is reused for the sums
    ZeroDArr(cur_len, hab);
    for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
      const double* cur_prob = &(prob[class_idx * sample_ctav + chunk_start]);
      const double* rr = &(ww[class_idx * sample_ctav + chunk_start]);
      for (uint32_t uii = 0; uii != cur_len; ++uii) {
        hab[uii] += cur_prob[uii] * rr[uii];
      }
    }
    for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
      const double* cur_prob = &(prob[class_idx * sample_ctav + chunk_start]);
      double* rr = &(ww[class_idx * sample_ctav + chunk_start]);
      for (uint32_t uii = 0; uii != cur_len; ++uii) {
        rr[uii] = 0.5 * cur_prob[uii] * (rr[uii] - hab[uii]);
      }
    }
  }
  double* ustar = fbp->ustar;
  double* xtw = fbp->xtw;
  for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
    ColMajorVectorMatrixMultiplyStrided(&(ww[class_idx * sample_ctav]), xx, sample_ct, sample_ctav, pred_ct, xtw);
    const double* grad_iter = &(fbp->grad[class_idx * pred_ct]);
    double* ustar_iter = &(ustar[class_idx * pred_ct]);
    for (uint32_t pred_idx = 0; pred_idx != pred_ct; ++pred_idx) {
      ustar_iter[pred_idx] = grad_iter[pred_idx] + xtw[pred_idx];
    }
  }
}

// The information matrix of the pseudo-data whose score is U* (Heinze and
// Schemper 2002), treating the hat blocks as fixed: the information matrix
// with sample i weighted by 1 + tr(W_i h_i) (1 + its leverage), at the point
// last passed to MultinomialFirthScore().  With two classes that is
// X^T diag((1 + h_i) v_i) X, FirthRegressionD()'s step matrix, whose inverse
// logistf reports as the covariance of the estimate.  Its lower triangle is
// left in fbp->info_aug.
static void MultinomialFirthAugInfo(const double* xx, uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, MultinomialFirthBufs* fbp) {
  MultinomialFirthGradAndInfo(xx, fbp->sqrt_wts, sample_ct, pred_ct, class_ct, fbp, nullptr, fbp->info_aug);
}

// Copies the entries of a param_ct-vector that belong to free coefficients
// (those of the first free_pred_ct predictors of each class) to dst, in
// order.
static void MultinomialPackFree(const double* src, uint32_t param_ct, uint32_t pred_ct, uint32_t free_pred_ct, double* dst) {
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
    if ((param_idx % pred_ct) < free_pred_ct) {
      *dst++ = src[param_idx];
    }
  }
}

// Solves I_aug step_f = ustar_f on the free coefficients (computing I_aug
// first), leaving the free block of I_aug's lower triangle, packed, in
// fbp->iinv and its Cholesky factor in fbp->chol.  Returns 1 if that block is
// not numerically positive definite.
static BoolErr MultinomialFirthAugStep(const double* xx, uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, uint32_t free_pred_ct, uint32_t free_param_ct, const double* ustar_f, MultinomialFirthBufs* fbp, double* step_f) {
  MultinomialFirthAugInfo(xx, sample_ct, pred_ct, class_ct, fbp);
  const uint32_t param_ct = pred_ct * class_ct;
  double* info_ff = fbp->iinv;
  uint32_t row_f = 0;
  for (uint32_t row_idx = 0; row_idx != param_ct; ++row_idx) {
    if ((row_idx % pred_ct) >= free_pred_ct) {
      continue;
    }
    const double* info_row = &(fbp->info_aug[row_idx * param_ct]);
    double* info_ff_row = &(info_ff[row_f * free_param_ct]);
    uint32_t col_f = 0;
    for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
      if ((col_idx % pred_ct) < free_pred_ct) {
        info_ff_row[col_f++] = info_row[col_idx];
      }
    }
    ++row_f;
  }
  if (MultinomialCholesky(info_ff, free_param_ct, fbp->chol)) {
    return 1;
  }
  MultinomialCholSolve(fbp->chol, ustar_f, free_param_ct, step_f);
  return 0;
}

// Checks the convergence criterion described above MultinomialFirthFit() for
// a packed step; also sets *max_step_ptr.
static uint32_t MultinomialFirthStepConverged(const double* step_f, const double* coefs, uint32_t param_ct, uint32_t pred_ct, uint32_t free_pred_ct, double* max_step_ptr) {
  double max_step = 0.0;
  uint32_t converged = 1;
  const double* step_iter = step_f;
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
    if ((param_idx % pred_ct) >= free_pred_ct) {
      continue;
    }
    const double abs_step = fabs(*step_iter++);
    if (abs_step > max_step) {
      max_step = abs_step;
    }
    if (abs_step >= 1e-10 * (1.0 + fabs(coefs[param_idx]))) {
      converged = 0;
    }
  }
  *max_step_ptr = max_step;
  return converged;
}

// Maximizes the Firth-penalized log-likelihood
//   l*(b) = l(b) + 0.5 log det I(b),
// I being the full information matrix.  Each iteration computes the modified
// score U* (MultinomialFirthScore()) and a step, capped at 5 per coefficient
// (as in logistf) and halved while l* decreases.  The step is first the
// quasi-Newton step I_aug^{-1} U* (I_aug as described above
// MultinomialFirthAugInfo(); this is FirthRegressionD()'s step), which
// converges linearly: fast when the penalty is small next to the information,
// but at a rate approaching 1 near a (quasi-)separated solution, where the
// penalty dominates the curvature of l* in the separated direction.  So from
// the third iteration on, once a step is more than kMultinomialFirthSlowRatio
// times the previous one, the step matrix becomes a BFGS approximation of the
// negated Hessian of l*, started from the current I_aug and updated with each
// step and the change in U* it produced; that converges superlinearly.  None
// of this changes the fixed point, U* = 0.
// Converged when the step computed at the current coefficients is below
// 1e-10 * (1 + |coefficient|) for every coefficient (with BFGS, the I_aug step
// must be too, which guards against a poor BFGS matrix stopping the iteration
// early); that step is then applied without re-evaluating anything: lli[],
// the log-determinant and I_aug are those at the point the step was taken
// from, which moves l* by about U* . step, far below anything reported.  If
// that has not happened after kMultinomialFirthMaxIter iterations, the current
// point is returned with *is_unfinished_ptr set (like FirthRegressionD()'s
// is_unfinished).
//
// xx: predictor-major, rows sample_ctav long, trailing elements zero; row 0
//   must be the intercept.
// classes: per-sample class index in [0, class_ct], 0 = reference class.
// coefs: starting point on input, estimate on output; class-major, class_ct x
//   pred_ct.
// free_pred_ct: the coefficients of predictors free_pred_ct and up are held
//   at zero, and l* (whose penalty is still that of the full model) is
//   maximized over the others: the restricted fit of a penalized
//   likelihood-ratio test, as in logistf.  The steps are then solved on the
//   free coefficients only.
// lli: set to each sample's (unpenalized) log-likelihood contribution, so
//   that l* = sum(lli) + 0.5 * (*logdet_ptr).
// cov: if not nullptr (unrestricted fits only), filled with I_aug^{-1} (full
//   param_ct x param_ct matrix), or cov[0] set to -1 if I_aug is singular.
//   This is the covariance logistf reports (its 'var', the inverse of
//   X^T diag((1 + h) v) X) and FirthRegressionD() returns, rather than the
//   inverse of I itself; for more than two classes, I_aug (each sample's
//   contribution to I weighted by 1 + its leverage) is the natural extension.
//   The likelihood-ratio statistic does not depend on this choice.
// Returns kGlmErrcodeNone or kGlmErrcodeFirthConvergeFail.
static GlmErrcode MultinomialFirthFit(const double* xx, const uint32_t* classes, uint32_t sample_ct, uint32_t pred_ct, uint32_t class_ct, uint32_t free_pred_ct, double* coefs, double* lli, double* logdet_ptr, double* cov, MultinomialFirthBufs* fbp, uint32_t* is_unfinished_ptr) {
  const uint32_t param_ct = pred_ct * class_ct;
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  *is_unfinished_ptr = 0;
  {
    // X^T Y, constant over the iterations
    double* xty = fbp->xty;
    ZeroDArr(param_ct, xty);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t cur_class = classes[sample_idx];
      if (cur_class) {
        double* xty_iter = &(xty[(cur_class - 1) * pred_ct]);
        for (uint32_t pred_idx = 0; pred_idx != pred_ct; ++pred_idx) {
          xty_iter[pred_idx] += xx[pred_idx * sample_ctav + sample_idx];
        }
      }
    }
  }
  for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
    ZeroDArr(pred_ct - free_pred_ct, &(coefs[class_idx * pred_ct + free_pred_ct]));
  }
  const uint32_t free_param_ct = class_ct * free_pred_ct;
  double ln_lik;
  double logdet;
  if (MultinomialFirthEval(xx, coefs, classes, sample_ct, pred_ct, class_ct, fbp, lli, &ln_lik, &logdet)) {
    return kGlmErrcodeFirthConvergeFail;
  }
  double lstar = ln_lik + 0.5 * logdet;
  double* step = fbp->step;
  double* ustar_f = fbp->ustar_f;
  double* ustar_old = fbp->ustar_old;
  double* step_f = fbp->step_f;
  double* sdiff = fbp->sdiff;
  double* bfgs = fbp->bfgs;
  // use_bfgs: bfgs[] (full symmetric, free_param_ct square) is the step
  //   matrix.
  // bfgs_pending: sdiff[] and ustar_old[] hold the last step taken and U*
  //   before it, not yet folded into bfgs[].
  // aug_is_current: fbp->info_aug is I_aug at the last point scored.
  uint32_t use_bfgs = 0;
  uint32_t bfgs_pending = 0;
  uint32_t aug_is_current = 0;
  double prev_max_step = 0.0;
  for (uint32_t iter_idx = 0; ; ++iter_idx) {
    MultinomialFirthScore(xx, sample_ct, pred_ct, class_ct, fbp);
    MultinomialPackFree(fbp->ustar, param_ct, pred_ct, free_pred_ct, ustar_f);
    aug_is_current = 0;
    if (use_bfgs) {
      if (bfgs_pending) {
        // B += y y^T / (y . s) - (B s)(B s)^T / (s . B s), with the gradient
        // of -l*, so y = U*_old - U*_new; skipped unless y . s > 0 (as it is
        // when l* is locally concave along s), which keeps B positive
        // definite.
        double* bs = fbp->tmpv;
        double sy = 0.0;
        double ss = 0.0;
        double yy = 0.0;
        double sbs = 0.0;
        for (uint32_t row_idx = 0; row_idx != free_param_ct; ++row_idx) {
          const double* bfgs_row = &(bfgs[row_idx * free_param_ct]);
          double dxx = 0.0;
          for (uint32_t col_idx = 0; col_idx != free_param_ct; ++col_idx) {
            dxx += bfgs_row[col_idx] * sdiff[col_idx];
          }
          bs[row_idx] = dxx;
          sbs += sdiff[row_idx] * dxx;
          const double cur_y = ustar_old[row_idx] - ustar_f[row_idx];
          ustar_old[row_idx] = cur_y;
          sy += sdiff[row_idx] * cur_y;
          ss += sdiff[row_idx] * sdiff[row_idx];
          yy += cur_y * cur_y;
        }
        if ((sy > 1e-12 * sqrt(ss * yy)) && (sbs > 0.0)) {
          const double* yvec = ustar_old;
          const double sy_recip = 1.0 / sy;
          const double sbs_recip = 1.0 / sbs;
          for (uint32_t row_idx = 0; row_idx != free_param_ct; ++row_idx) {
            double* bfgs_row = &(bfgs[row_idx * free_param_ct]);
            const double ymult = yvec[row_idx] * sy_recip;
            const double bsmult = bs[row_idx] * sbs_recip;
            for (uint32_t col_idx = 0; col_idx != free_param_ct; ++col_idx) {
              bfgs_row[col_idx] += ymult * yvec[col_idx] - bsmult * bs[col_idx];
            }
          }
        }
        bfgs_pending = 0;
      }
      if (MultinomialCholesky(bfgs, free_param_ct, fbp->chol)) {
        // shouldn't happen; fall back on the I_aug step
        use_bfgs = 0;
      } else {
        MultinomialCholSolve(fbp->chol, ustar_f, free_param_ct, step_f);
      }
    }
    if (!use_bfgs) {
      if (MultinomialFirthAugStep(xx, sample_ct, pred_ct, class_ct, free_pred_ct, free_param_ct, ustar_f, fbp, step_f)) {
        return kGlmErrcodeFirthConvergeFail;
      }
      aug_is_current = 1;
    }
    double max_step;
    uint32_t converged = MultinomialFirthStepConverged(step_f, coefs, param_ct, pred_ct, free_pred_ct, &max_step);
    if (converged && use_bfgs) {
      double* aug_step_f = fbp->tmpv;
      if (MultinomialFirthAugStep(xx, sample_ct, pred_ct, class_ct, free_pred_ct, free_param_ct, ustar_f, fbp, aug_step_f)) {
        return kGlmErrcodeFirthConvergeFail;
      }
      aug_is_current = 1;
      double aug_max_step;
      converged = MultinomialFirthStepConverged(aug_step_f, coefs, param_ct, pred_ct, free_pred_ct, &aug_max_step);
    }
    if (converged) {
      const double* step_f_iter = step_f;
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        if ((param_idx % pred_ct) < free_pred_ct) {
          coefs[param_idx] += *step_f_iter++;
        }
      }
      break;
    }
    if (iter_idx == kMultinomialFirthMaxIter) {
      *is_unfinished_ptr = 1;
      break;
    }
    if ((!use_bfgs) && (iter_idx >= 2) && (max_step > kMultinomialFirthSlowRatio * prev_max_step)) {
      // Switch to BFGS, starting from I_aug here (whose packed lower
      // triangle MultinomialFirthAugStep() left in fbp->iinv); this
      // iteration's step is unchanged.
      use_bfgs = 1;
      const double* info_ff = fbp->iinv;
      for (uint32_t row_idx = 0; row_idx != free_param_ct; ++row_idx) {
        for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
          const double cur_val = info_ff[row_idx * free_param_ct + col_idx];
          bfgs[row_idx * free_param_ct + col_idx] = cur_val;
          bfgs[col_idx * free_param_ct + row_idx] = cur_val;
        }
      }
    }
    prev_max_step = max_step;
    const double scale = (max_step > 5.0)? (5.0 / max_step) : 1.0;
    {
      const double* step_f_iter = step_f;
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        step[param_idx] = ((param_idx % pred_ct) >= free_pred_ct)? 0.0 : (scale * (*step_f_iter++));
      }
    }
    memcpy(fbp->coef_old, coefs, param_ct * sizeof(double));
    const double lstar_tol = 1e-10 * (0.1 + fabs(lstar));
    uint32_t halving_ct = 0;
    while (1) {
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        coefs[param_idx] = fbp->coef_old[param_idx] + step[param_idx];
      }
      if (!MultinomialFirthEval(xx, coefs, classes, sample_ct, pred_ct, class_ct, fbp, lli, &ln_lik, &logdet)) {
        const double new_lstar = ln_lik + 0.5 * logdet;
        if (new_lstar > lstar - lstar_tol) {
          lstar = new_lstar;
          break;
        }
      }
      if (++halving_ct == kMultinomialMaxHalvings) {
        return kGlmErrcodeFirthConvergeFail;
      }
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        step[param_idx] *= 0.5;
      }
    }
    if (use_bfgs) {
      MultinomialPackFree(step, param_ct, pred_ct, free_pred_ct, sdiff);
      memcpy(ustar_old, ustar_f, free_param_ct * sizeof(double));
      bfgs_pending = 1;
    }
  }
  *logdet_ptr = logdet;
  if (!cov) {
    return kGlmErrcodeNone;
  }
  assert(free_pred_ct == pred_ct);
  if (!aug_is_current) {
    MultinomialFirthAugInfo(xx, sample_ct, pred_ct, class_ct, fbp);
  }
  if (MultinomialCholesky(fbp->info_aug, param_ct, fbp->chol)) {
    cov[0] = -1.0;
  } else {
    MultinomialCholInvFull(fbp->chol, param_ct, fbp->linv, cov);
  }
  return kGlmErrcodeNone;
}

static inline uintptr_t MultinomialInvertDbl2dCt(uintptr_t dim) {
  return dim * MAXV(dim, 7);
}

BoolErr GlmAllocFillAndTestPhenoCovarsMultinomial(const uintptr_t* sample_include, const PhenoCol* pheno_col, const uint32_t* cat_to_level, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uint32_t level_ct, uintptr_t covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uint32_t firth_mode, MultinomialSet* setp, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  *glm_err_ptr = 0;
  setp->firth_null = 0;
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
    if (!firth_mode) {
      *glm_err_ptr = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
    } else {
      // A covariate separates the levels, so the covariate-only model has no
      // maximum-likelihood estimate.  As GlmLogistic() does for separated
      // covariates, fit it with Firth regression, and use Firth regression
      // for every variant.  (covars_pmaj is already laid out the way
      // MultinomialFirthFit() wants it.)
      BigstackReset(fit_wkspace);
      unsigned char* firth_buf;
      double* lli;
      if (unlikely(bigstack_alloc_uc(MultinomialFirthBufsLayout(sample_ct, pred_ct, nonref_class_ct, nullptr, nullptr), &firth_buf) ||
                   bigstack_alloc_d(sample_ct, &lli))) {
        return 1;
      }
      MultinomialFirthBufs fb;
      MultinomialFirthBufsLayout(sample_ct, pred_ct, nonref_class_ct, firth_buf, &fb);
      ZeroDArr(nonref_class_ct * S_CAST(uintptr_t, pred_ct), null_coefs);
      for (uint32_t class_idx = 1; class_idx != class_ct; ++class_idx) {
        null_coefs[(class_idx - 1) * pred_ct] = log(u31tod(class_sample_cts[class_idx]) / ref_ct_d);
      }
      double logdet;
      if (MultinomialFirthFit(covars_pmaj, classes, sample_ct, pred_ct, nonref_class_ct, pred_ct, null_coefs, lli, &logdet, nullptr, &fb, &is_unfinished) || is_unfinished) {
        *glm_err_ptr = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
      } else {
        setp->firth_null = 1;
      }
    }
  }
  BigstackReset(class_sample_cts);
  return 0;
}

// Per-thread workspace layout, shared by the size computation and the compute
// thread so the two can't drift apart.
typedef struct {
  uintptr_t* sample_nm;
  uintptr_t* tmp_nm;
  // global level, then class index, of each sample with a call
  uint32_t* nm_levels;
  uint32_t* nm_classes;
  uint32_t* level_sample_cts;
  uint32_t* level_to_class;
  uint32_t* class_levels;
  double* class_min;
  double* class_max;
  // for each genotype predictor: allele it describes, and its output slot
  uint32_t* geno_pred_alleles;
  uint32_t* geno_pred_slots;
  const double** pred_rows;
  // (covar_ct + 1) rows, stride nm_sample_ctav; only used when calls are
  // missing
  double* nm_covars;
  // one row per allele (allele counts or dosages, haploid-scaled), stride
  // sample_ctav
  double* allele_cols;
  // centered genotype predictors: nuisance alleles, then tested columns, then
  // one spare row
  double* geno_preds;
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
  // Firth regression only: the predictors laid out as MultinomialFirthFit()
  // wants them, per-sample log-likelihoods of the restricted and full fits,
  // and the MultinomialFirthBufs area
  double* firth_xx;
  double* firth_lli0;
  double* firth_lli1;
  unsigned char* firth_buf;
} MultinomialWorkspace;

// With workspace_buf == nullptr, only returns the size.
static uintptr_t MultinomialWorkspaceLayout(uint32_t sample_ct, uint32_t covar_ct, uint32_t level_ct, uint32_t max_allele_ct, uint32_t geno_pred_max, uint32_t use_firth, unsigned char* workspace_buf, MultinomialWorkspace* wsp) {
  MultinomialWorkspace size_only_ws;
  if (!workspace_buf) {
    wsp = &size_only_ws;
  }
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uint32_t max_pred_ct = covar_ct + 1 + geno_pred_max;
  const uint32_t nonref_level_ct = level_ct - 1;
  const uintptr_t max_param_ct = nonref_level_ct * S_CAST(uintptr_t, max_pred_ct);
  const uint32_t nonintercept_ct = max_pred_ct - 1;
  const uintptr_t byte_cts[] = {
    BitCtToWordCt(sample_ct) * sizeof(intptr_t),
    BitCtToWordCt(sample_ct) * sizeof(intptr_t),
    sample_ct * sizeof(int32_t),
    sample_ct * sizeof(int32_t),
    level_ct * sizeof(int32_t),
    level_ct * sizeof(int32_t),
    level_ct * sizeof(int32_t),
    level_ct * sizeof(double),
    level_ct * sizeof(double),
    geno_pred_max * sizeof(int32_t),
    geno_pred_max * sizeof(int32_t),
    max_pred_ct * sizeof(intptr_t),
    (covar_ct + 1) * sample_ctav * sizeof(double),
    max_allele_ct * sample_ctav * sizeof(double),
    (geno_pred_max + 1) * sample_ctav * sizeof(double),
    max_param_ct * sizeof(double),
    max_param_ct * sizeof(double),
    max_param_ct * max_param_ct * sizeof(double),
    MultinomialFitWkspaceDoubleCt(max_param_ct) * sizeof(double),
    MultinomialChunkBufDoubleCt(max_pred_ct, nonref_level_ct) * sizeof(double),
    max_param_ct * kMatrixInvertBuf1CheckedAlloc,
    MultinomialInvertDbl2dCt(max_param_ct) * sizeof(double),
    nonintercept_ct * nonintercept_ct * sizeof(double),
    nonintercept_ct * nonintercept_ct * sizeof(double),
    use_firth? (max_pred_ct * sample_ctav * sizeof(double)) : 0,
    use_firth? (sample_ct * sizeof(double)) : 0,
    use_firth? (sample_ct * sizeof(double)) : 0,
    use_firth? MultinomialFirthBufsLayout(sample_ct, max_pred_ct, nonref_level_ct, nullptr, nullptr) : 0
  };
  void** dsts[] = {
    R_CAST(void**, &wsp->sample_nm),
    R_CAST(void**, &wsp->tmp_nm),
    R_CAST(void**, &wsp->nm_levels),
    R_CAST(void**, &wsp->nm_classes),
    R_CAST(void**, &wsp->level_sample_cts),
    R_CAST(void**, &wsp->level_to_class),
    R_CAST(void**, &wsp->class_levels),
    R_CAST(void**, &wsp->class_min),
    R_CAST(void**, &wsp->class_max),
    R_CAST(void**, &wsp->geno_pred_alleles),
    R_CAST(void**, &wsp->geno_pred_slots),
    R_CAST(void**, &wsp->pred_rows),
    R_CAST(void**, &wsp->nm_covars),
    R_CAST(void**, &wsp->allele_cols),
    R_CAST(void**, &wsp->geno_preds),
    R_CAST(void**, &wsp->null_coefs),
    R_CAST(void**, &wsp->alt_coefs),
    R_CAST(void**, &wsp->cov),
    R_CAST(void**, &wsp->fit_wkspace),
    R_CAST(void**, &wsp->chunk_buf),
    R_CAST(void**, &wsp->mi_buf),
    R_CAST(void**, &wsp->dbl_2d_buf),
    R_CAST(void**, &wsp->dotprods),
    R_CAST(void**, &wsp->inverse_corr_buf),
    R_CAST(void**, &wsp->firth_xx),
    R_CAST(void**, &wsp->firth_lli0),
    R_CAST(void**, &wsp->firth_lli1),
    R_CAST(void**, &wsp->firth_buf)
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

// Number of genotype predictors a regression can have: in the additive model,
// every non-omitted allele (tested jointly); otherwise the tested columns plus
// the other non-omitted alleles as additive nuisance columns.
static inline uint32_t MultinomialGenoPredMax(const GlmMultinomialCtx* ctx) {
  const uint32_t max_extra_allele_ct = ctx->common->max_extra_allele_ct;
  return ctx->is_additive? (max_extra_allele_ct + 1) : (max_extra_allele_ct + ctx->model_col_ct);
}

// Classifies a predictor column over the samples of the current regression:
//   0: ordinary
//   1: constant
//   2: some class's values all lie at or below every other class's value
//   3: some class's values all lie at or above every other class's value
// In cases 2 and 3, moving that class's coefficient for the column (with an
// offsetting intercept) never lowers any sample's likelihood, so the
// maximum-likelihood estimate does not exist.
static uint32_t ClassifyPredictor(const double* col, const uint32_t* classes, uint32_t sample_ct, uint32_t class_ct, double* class_min, double* class_max) {
  for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
    class_min[class_idx] = DBL_MAX;
    class_max[class_idx] = -DBL_MAX;
  }
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uint32_t class_idx = classes[sample_idx];
    const double cur_val = col[sample_idx];
    if (cur_val < class_min[class_idx]) {
      class_min[class_idx] = cur_val;
    }
    if (cur_val > class_max[class_idx]) {
      class_max[class_idx] = cur_val;
    }
  }
  double tot_min = class_min[0];
  double tot_max = class_max[0];
  for (uint32_t class_idx = 1; class_idx != class_ct; ++class_idx) {
    tot_min = MINV(tot_min, class_min[class_idx]);
    tot_max = MAXV(tot_max, class_max[class_idx]);
  }
  if (tot_min == tot_max) {
    return 1;
  }
  for (uint32_t class_idx = 0; class_idx != class_ct; ++class_idx) {
    double rest_min = DBL_MAX;
    double rest_max = -DBL_MAX;
    for (uint32_t class_idx2 = 0; class_idx2 != class_ct; ++class_idx2) {
      if (class_idx2 != class_idx) {
        rest_min = MINV(rest_min, class_min[class_idx2]);
        rest_max = MAXV(rest_max, class_max[class_idx2]);
      }
    }
    if (class_max[class_idx] <= rest_min) {
      return 2;
    }
    if (class_min[class_idx] >= rest_max) {
      return 3;
    }
  }
  return 0;
}

// Copies src into dst, centered.  Centering only reparameterizes the
// intercepts.
static void CopyCentered(const double* src, uint32_t ct, double* dst) {
  double sum = 0.0;
  for (uint32_t uii = 0; uii != ct; ++uii) {
    sum += src[uii];
  }
  const double mean = sum / u31tod(ct);
  for (uint32_t uii = 0; uii != ct; ++uii) {
    dst[uii] = src[uii] - mean;
  }
}

// Wald statistic gamma' Cov(gamma)^{-1} gamma over the tested columns (the
// last tested_ct predictors of each class), given the full inverse
// information matrix cov[] of a fit with class_ct non-reference classes and
// pred_ct predictors, the first null_pred_ct untested.  wkspace needs
// (class_ct * tested_ct) * (class_ct * tested_ct + 1) doubles.  Returns 1 if
// the tested block of cov[] is singular.
static BoolErr MultinomialWaldChisq(const double* coefs, const double* cov, uint32_t class_ct, uint32_t pred_ct, uint32_t null_pred_ct, uint32_t tested_ct, double* wkspace, MatrixInvertBuf1* mi_buf, double* dbl_2d_buf, double* chisq_ptr) {
  const uintptr_t param_ct = class_ct * S_CAST(uintptr_t, pred_ct);
  const uint32_t wald_dim = class_ct * tested_ct;
  double* wald_cov = wkspace;
  double* wald_coefs = &(wald_cov[wald_dim * wald_dim]);
  for (uint32_t wald_idx1 = 0; wald_idx1 != wald_dim; ++wald_idx1) {
    const uintptr_t param_idx1 = (wald_idx1 / tested_ct) * pred_ct + null_pred_ct + (wald_idx1 % tested_ct);
    wald_coefs[wald_idx1] = coefs[param_idx1];
    for (uint32_t wald_idx2 = 0; wald_idx2 != wald_dim; ++wald_idx2) {
      const uintptr_t param_idx2 = (wald_idx2 / tested_ct) * pred_ct + null_pred_ct + (wald_idx2 % tested_ct);
      wald_cov[wald_idx1 * wald_dim + wald_idx2] = cov[param_idx1 * param_ct + param_idx2];
    }
  }
  if (InvertSymmdefMatrixChecked(wald_dim, wald_cov, mi_buf, dbl_2d_buf)) {
    return 1;
  }
  ReflectMatrix(wald_dim, wald_cov);
  double chisq = 0.0;
  for (uint32_t wald_idx = 0; wald_idx != wald_dim; ++wald_idx) {
    chisq += wald_coefs[wald_idx] * DotprodD(&(wald_cov[wald_idx * wald_dim]), wald_coefs, wald_dim);
  }
  *chisq_ptr = chisq;
  return 0;
}

THREAD_FUNC_DECL GlmMultinomialThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GlmMultinomialCtx* ctx = S_CAST(GlmMultinomialCtx*, arg->sharedp->context);
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
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const GlmFlags glm_flags = common->glm_flags;
  const uint32_t level_ct = ctx->level_ct;
  const uint32_t nonref_level_ct = level_ct - 1;
  const GlmMultinomialTest test_type = ctx->test_type;
  const uint32_t is_score_test = (test_type == kGlmMultinomialTestScore);
  const uint32_t firth_mode = ctx->firth_mode;
  const uint32_t save_coefs = ctx->save_coefs;
  const uint32_t is_additive = ctx->is_additive;
  const uint32_t model_col_ct = ctx->model_col_ct;
  const uint32_t max_row_ct = ctx->max_row_ct;
  const uint32_t max_a1_ct = ctx->max_a1_ct;
  const uint32_t max_tested_ct = ctx->max_tested_ct;
  const uint32_t max_allele_ct = common->max_extra_allele_ct + 2;
  const uint32_t geno_pred_max = MultinomialGenoPredMax(ctx);
  const uintptr_t row_beta_se_stride = 2 * S_CAST(uintptr_t, nonref_level_ct) * max_tested_ct;
  // the Wald test and the coefficient columns need the full fit's inverse
  // information matrix
  const uint32_t need_cov = save_coefs || (test_type == kGlmMultinomialTestWald);
  uintptr_t max_sample_ct = MAXV(common->sample_ct, common->sample_ct_x);
  if (max_sample_ct < common->sample_ct_y) {
    max_sample_ct = common->sample_ct_y;
  }
  SetPgvThreadMhcNull(max_sample_ct, tidx, common->thread_mhc, &pgv);
  pgv.patch_01_ct = 0;
  pgv.patch_10_ct = 0;
  pgv.multidosage_sample_ct = 0;
  uint32_t variant_idx_offset = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);

    MultinomialAuxResult* block_aux_iter = &(ctx->block_aux[variant_bidx * S_CAST(uintptr_t, max_row_ct)]);
    double* a1_dosage_iter = &(ctx->block_a1_dosage[variant_bidx * S_CAST(uintptr_t, max_row_ct) * max_a1_ct]);
    double* level_a1_iter = &(ctx->block_level_a1[variant_bidx * S_CAST(uintptr_t, max_row_ct) * level_ct * max_a1_ct]);
    uint32_t* level_allele_obs_iter = &(ctx->block_level_allele_obs[variant_bidx * S_CAST(uintptr_t, level_ct)]);
    double* beta_se_iter = &(ctx->block_beta_se[variant_bidx * max_row_ct * row_beta_se_stride]);
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
      const uint32_t base_pred_ct = covar_ct + 1;
      const uint32_t* set_levels = setp->levels;
      const uint32_t always_firth = (firth_mode == 2) || setp->firth_null;
      MultinomialWorkspace ws;
      MultinomialWorkspaceLayout(cur_sample_ct, covar_ct, level_ct, max_allele_ct, geno_pred_max, firth_mode != 0, workspace_buf, &ws);
      uintptr_t* sample_nm = ws.sample_nm;
      uint32_t* nm_levels = ws.nm_levels;
      uint32_t* nm_classes = ws.nm_classes;
      uint32_t* level_sample_cts = ws.level_sample_cts;
      uint32_t* level_to_class = ws.level_to_class;
      uint32_t* class_levels = ws.class_levels;
      double* allele_cols = ws.allele_cols;
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        uint32_t allele_ct = 2;
        if (allele_idx_offsets) {
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
        }
        const uint32_t omitted_allele_idx = omitted_alleles? omitted_alleles[variant_uidx] : 0;
        PglErr reterr;
        if (allele_ct == 2) {
          reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
        } else {
          reterr = PgrGetMD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, &pgv);
          // todo: multiallelic dosages
          assert(!pgv.dosage_ct);
        }
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto GlmMultinomialThread_err;
        }
        uintptr_t* genovec = pgv.genovec;
        const uint32_t dosage_ct = pgv.dosage_ct;
        ZeroTrailingNyps(cur_sample_ct, genovec);
        GenoarrCountFreqsUnsafe(genovec, cur_sample_ct, genocounts);
        uint32_t missing_ct = genocounts[3];
        if (!missing_ct) {
          SetAllBits(cur_sample_ct, sample_nm);
        } else {
          GenoarrToNonmissing(genovec, cur_sample_ct, sample_nm);
          if (dosage_ct) {
            BitvecOr(pgv.dosage_present, sample_ctl, sample_nm);
            missing_ct = cur_sample_ct - PopcountWords(sample_nm, sample_ctl);
          }
        }
        const uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
        const uintptr_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kDoublePerDVec);
        uint64_t machr2_dosage_sums[kPglMaxAlleleCt];
        uint64_t machr2_dosage_ssqs[kPglMaxAlleleCt];
        if (allele_ct == 2) {
          // A1 dosage of the samples with a call, in 0..2 units.
          const uint32_t a1_allele_idx = 1 - omitted_allele_idx;
          double* a1_col = &(allele_cols[a1_allele_idx * sample_ctav]);
          if (omitted_allele_idx) {
            GenovecInvertUnsafe(cur_sample_ct, genovec);
            if (dosage_ct) {
              BiallelicDosage16Invert(dosage_ct, pgv.dosage_main);
            }
            const uint32_t uii = genocounts[0];
            genocounts[0] = genocounts[2];
            genocounts[2] = uii;
          }
          uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
          uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
          if (!missing_ct) {
            GenoarrLookup16x8bx2(genovec, kSmallDoublePairs, nm_sample_ct, a1_col);
            if (dosage_ct) {
              uintptr_t sample_idx_base = 0;
              uintptr_t dosage_present_bits = pgv.dosage_present[0];
              for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
                const uintptr_t sample_idx = BitIter1(pgv.dosage_present, &sample_idx_base, &dosage_present_bits);
                const uint32_t dosage_val = pgv.dosage_main[dosage_idx];
                a1_col[sample_idx] = kRecipDosageMid * u31tod(dosage_val);
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
            GenoarrToDoublesRemoveMissing(genovec, kSmallDoubles, cur_sample_ct, a1_col);
          } else {
            uintptr_t sample_midx_base = 0;
            uintptr_t sample_nm_bits = sample_nm[0];
            uint32_t dosage_idx = 0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
              const uintptr_t cur_geno = GetNyparrEntry(genovec, sample_midx);
              double cur_val;
              if (IsSet(pgv.dosage_present, sample_midx)) {
                const uint32_t dosage_val = pgv.dosage_main[dosage_idx++];
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
              a1_col[sample_idx] = cur_val;
            }
          }
          machr2_dosage_sums[a1_allele_idx] = dosage_sum;
          machr2_dosage_ssqs[a1_allele_idx] = dosage_ssq;
          machr2_dosage_sums[omitted_allele_idx] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[omitted_allele_idx] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
        } else {
          // Multiallelic hardcalls: one count column per allele.
          double* ref_col = allele_cols;
          double* alt1_col = &(allele_cols[sample_ctav]);
          if (!missing_ct) {
            GenoarrLookup16x8bx2(genovec, kSmallInvDoublePairs, nm_sample_ct, ref_col);
            GenoarrLookup16x8bx2(genovec, kSmallDoublePairs, nm_sample_ct, alt1_col);
          } else {
            GenoarrToDoublesRemoveMissing(genovec, kSmallInvDoubles, cur_sample_ct, ref_col);
            GenoarrToDoublesRemoveMissing(genovec, kSmallDoubles, cur_sample_ct, alt1_col);
          }
          for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
            ZeroDArr(nm_sample_ct, &(allele_cols[allele_idx * sample_ctav]));
          }
          if (pgv.patch_01_ct) {
            // ref/altx hets
            const uintptr_t* patch_set_nm = pgv.patch_01_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_01_set, sample_nm, nm_sample_ct, ws.tmp_nm);
              patch_set_nm = ws.tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
              alt1_col[sample_idx] -= 1.0;
              allele_cols[pgv.patch_01_vals[uii] * sample_ctav + sample_idx] += 1.0;
            }
          }
          if (pgv.patch_10_ct) {
            // altx/alty
            const uintptr_t* patch_set_nm = pgv.patch_10_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_10_set, sample_nm, nm_sample_ct, ws.tmp_nm);
              patch_set_nm = ws.tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
              alt1_col[sample_idx] -= 2.0;
              allele_cols[pgv.patch_10_vals[2 * uii] * sample_ctav + sample_idx] += 1.0;
              allele_cols[pgv.patch_10_vals[2 * uii + 1] * sample_ctav + sample_idx] += 1.0;
            }
          }
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            const double* col = &(allele_cols[allele_idx * sample_ctav]);
            double sum = 0.0;
            double ssq = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              sum += col[sample_idx];
              ssq += col[sample_idx] * col[sample_idx];
            }
            machr2_dosage_sums[allele_idx] = S_CAST(uint64_t, sum) * 0x4000LLU;
            machr2_dosage_ssqs[allele_idx] = S_CAST(uint64_t, ssq) * 0x10000000LLU;
          }
        }
        const double mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, allele_ct);
        // Haploid scaling, per-level counts, and the covariates of the
        // samples with a call.
        ZeroU32Arr(level_ct, level_sample_cts);
        ZeroU32Arr(level_ct, level_allele_obs_iter);
        {
          uintptr_t sample_idx_base = 0;
          uintptr_t sample_nm_bits = sample_nm[0];
          for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
            const uintptr_t sample_idx = BitIter1(sample_nm, &sample_idx_base, &sample_nm_bits);
            const uint32_t cur_level = set_levels[sample_idx];
            uint32_t ploidy = 2;
            if (is_nonx_haploid || (is_regular_x && is_xchr_model_1 && IsSet(sex_male_collapsed, sample_idx))) {
              ploidy = 1;
              // (only A1's column is filled in the biallelic case)
              for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                if ((allele_ct > 2) || (allele_idx != omitted_allele_idx)) {
                  allele_cols[allele_idx * sample_ctav + nm_idx] *= 0.5;
                }
              }
            }
            nm_levels[nm_idx] = cur_level;
            level_sample_cts[cur_level] += 1;
            level_allele_obs_iter[cur_level] += ploidy;
            if (missing_ct) {
              for (uint32_t pred_idx = 0; pred_idx != base_pred_ct; ++pred_idx) {
                ws.nm_covars[pred_idx * nm_sample_ctav + nm_idx] = setp->covars_pmaj[pred_idx * sample_ctav + sample_idx];
              }
            }
          }
        }
        uint32_t allele_obs_ct = 0;
        uint32_t min_level_allele_obs = UINT32_MAX;
        uint32_t class_ct = 0;
        for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
          const uint32_t cur_allele_obs = level_allele_obs_iter[level_idx];
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
        for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
          nm_classes[nm_idx] = level_to_class[nm_levels[nm_idx]];
        }
        const uint32_t nonref_class_ct = class_ct - 1;
        const double allele_obs_d = u31tod(allele_obs_ct);
        // Allele totals, for MIN_EXPECTED: the smallest expected cell of the
        // allele x level table under independence.  (In the additive model the
        // table has a row per allele; otherwise its rows are A1 and the other
        // alleles combined.)
        double min_allele_total = DBL_MAX;
        if (is_additive) {
          double non_omitted_total = 0.0;
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            if (allele_idx == omitted_allele_idx) {
              continue;
            }
            const double* col = &(allele_cols[allele_idx * sample_ctav]);
            double total = 0.0;
            for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
              total += col[nm_idx];
            }
            non_omitted_total += total;
            min_allele_total = MINV(min_allele_total, total);
          }
          min_allele_total = MINV(min_allele_total, allele_obs_d - non_omitted_total);
        }
        const double* const* pred_rows = ws.pred_rows;
        {
          const double** pred_rows_w = ws.pred_rows;
          for (uint32_t pred_idx = 0; pred_idx != base_pred_ct; ++pred_idx) {
            pred_rows_w[pred_idx] = missing_ct? (&(ws.nm_covars[pred_idx * nm_sample_ctav])) : (&(setp->covars_pmaj[pred_idx * sample_ctav]));
          }
          for (uint32_t geno_pred_idx = 0; geno_pred_idx != geno_pred_max; ++geno_pred_idx) {
            pred_rows_w[base_pred_ct + geno_pred_idx] = &(ws.geno_preds[geno_pred_idx * sample_ctav]);
          }
        }
        const uint32_t row_ct = is_additive? 1 : (allele_ct - 1);
        uint32_t row_a1_allele_idx = (omitted_allele_idx == 0);
        for (uint32_t row_idx = 0; row_idx != row_ct; ++row_idx) {
          MultinomialAuxResult* auxp = &(block_aux_iter[row_idx]);
          double* row_a1_dosages = &(a1_dosage_iter[row_idx * max_a1_ct]);
          double* row_level_a1 = &(level_a1_iter[row_idx * S_CAST(uintptr_t, level_ct) * max_a1_ct]);
          double* row_beta_se = &(beta_se_iter[row_idx * row_beta_se_stride]);
          if (row_idx) {
            ++row_a1_allele_idx;
            row_a1_allele_idx += (row_a1_allele_idx == omitted_allele_idx);
          }
          auxp->sample_obs_ct = nm_sample_ct;
          auxp->allele_obs_ct = allele_obs_ct;
          auxp->mach_r2 = mach_r2;
          auxp->chisq = -9.0;
          auxp->glm_err = 0;
          auxp->is_unfinished = 0;
          auxp->is_firth = 0;
          auxp->df = 0;
          for (uintptr_t ulii = 0; ulii != row_beta_se_stride; ulii += 2) {
            row_beta_se[ulii + 1] = -9.0;
          }
          // A1 totals, overall and per level.
          {
            uint32_t a1_idx = 0;
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if ((allele_idx == omitted_allele_idx) || ((!is_additive) && (allele_idx != row_a1_allele_idx))) {
                continue;
              }
              const double* col = &(allele_cols[allele_idx * sample_ctav]);
              for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
                row_level_a1[level_idx * max_a1_ct + a1_idx] = 0.0;
              }
              double total = 0.0;
              for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
                row_level_a1[nm_levels[nm_idx] * max_a1_ct + a1_idx] += col[nm_idx];
                total += col[nm_idx];
              }
              row_a1_dosages[a1_idx] = total;
              ++a1_idx;
            }
            if (!is_additive) {
              min_allele_total = MINV(row_a1_dosages[0], allele_obs_d - row_a1_dosages[0]);
            }
            auxp->min_expected = allele_obs_ct? (min_allele_total * u31tod(min_level_allele_obs) / allele_obs_d) : -9.0;
          }
          GlmErr glm_err = 0;
          uint32_t coefs_unavailable = 0;
          uint32_t is_unfinished;
          double null_ln_lik;
          uint32_t nuisance_ct = 0;
          uint32_t geno_pred_ct = 0;
          uint32_t tested_ct;
          uint32_t null_pred_ct;
          uint32_t alt_pred_ct;
          uintptr_t alt_param_ct;
          double* null_coefs = ws.null_coefs;
          double* alt_coefs = ws.alt_coefs;
          uint32_t* geno_pred_alleles = ws.geno_pred_alleles;
          uint32_t* geno_pred_slots = ws.geno_pred_slots;
          {
            const uint32_t max_geno_pred_ct = is_additive? (allele_ct - 1) : (allele_ct - 2 + model_col_ct);
            if (nm_sample_ct <= base_pred_ct + max_geno_pred_ct) {
              glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
              goto GlmMultinomialThread_skip_regression;
            }
          }
          if (class_ct < 2) {
            glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
            goto GlmMultinomialThread_skip_regression;
          }
          // Genotype predictors: the nuisance alleles, then the tested
          // columns.  Constant nuisance alleles, and in the additive model
          // constant tested alleles, are left out.
          if (!is_additive) {
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if ((allele_idx == omitted_allele_idx) || (allele_idx == row_a1_allele_idx)) {
                continue;
              }
              const double* col = &(allele_cols[allele_idx * sample_ctav]);
              if (ClassifyPredictor(col, nm_classes, nm_sample_ct, class_ct, ws.class_min, ws.class_max) == 1) {
                continue;
              }
              CopyCentered(col, nm_sample_ct, &(ws.geno_preds[geno_pred_ct * sample_ctav]));
              geno_pred_alleles[geno_pred_ct] = allele_idx;
              geno_pred_slots[geno_pred_ct] = UINT32_MAX;
              ++geno_pred_ct;
            }
            nuisance_ct = geno_pred_ct;
            // Recode A1 per the model.
            const double* a1_col = &(allele_cols[row_a1_allele_idx * sample_ctav]);
            double* main_col = &(ws.geno_preds[geno_pred_ct * sample_ctav]);
            double* second_col = &(main_col[sample_ctav]);
            for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
              const double cur_val = a1_col[nm_idx];
              // 0..1..0 dominance deviation
              const double domdev_val = (cur_val > 1.0)? (2.0 - cur_val) : cur_val;
              if (glm_flags & kfGlmDominant) {
                // 0..1..1
                main_col[nm_idx] = MINV(cur_val, 1.0);
              } else if (glm_flags & (kfGlmRecessive | kfGlmHethom)) {
                // 0..0..1
                main_col[nm_idx] = (cur_val < 1.0)? 0.0 : (cur_val - 1.0);
              } else if (glm_flags & kfGlmHetonly) {
                main_col[nm_idx] = domdev_val;
              } else {
                // genotypic
                main_col[nm_idx] = cur_val;
              }
              if (model_col_ct == 2) {
                second_col[nm_idx] = domdev_val;
              }
            }
            for (uint32_t col_idx = 0; col_idx != model_col_ct; ++col_idx) {
              double* col = &(main_col[col_idx * sample_ctav]);
              if (ClassifyPredictor(col, nm_classes, nm_sample_ct, class_ct, ws.class_min, ws.class_max) == 1) {
                glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
                goto GlmMultinomialThread_skip_regression;
              }
              CopyCentered(col, nm_sample_ct, col);
              geno_pred_alleles[geno_pred_ct] = row_a1_allele_idx;
              geno_pred_slots[geno_pred_ct] = col_idx;
              ++geno_pred_ct;
            }
          } else {
            uint32_t slot_idx = 0;
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if (allele_idx == omitted_allele_idx) {
                continue;
              }
              const double* col = &(allele_cols[allele_idx * sample_ctav]);
              if (ClassifyPredictor(col, nm_classes, nm_sample_ct, class_ct, ws.class_min, ws.class_max) != 1) {
                CopyCentered(col, nm_sample_ct, &(ws.geno_preds[geno_pred_ct * sample_ctav]));
                geno_pred_alleles[geno_pred_ct] = allele_idx;
                geno_pred_slots[geno_pred_ct] = slot_idx;
                ++geno_pred_ct;
              }
              ++slot_idx;
            }
            if (!geno_pred_ct) {
              glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
              goto GlmMultinomialThread_skip_regression;
            }
          }
          tested_ct = geno_pred_ct - nuisance_ct;
          null_pred_ct = base_pred_ct + nuisance_ct;
          alt_pred_ct = null_pred_ct + tested_ct;
          alt_param_ct = nonref_class_ct * S_CAST(uintptr_t, alt_pred_ct);
          auxp->df = nonref_class_ct * tested_ct;
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
          // Separation.  A separating nuisance allele leaves no
          // maximum-likelihood estimate under the null hypothesis either, so
          // no test survives it; a separating tested column only rules out
          // the likelihood ratio and Wald tests.  Firth regression handles
          // both.
          if (!always_firth) {
            uint32_t sep_allele_idx = UINT32_MAX;
            for (uint32_t geno_pred_idx = 0; geno_pred_idx != geno_pred_ct; ++geno_pred_idx) {
              const uint32_t classification = ClassifyPredictor(pred_rows[base_pred_ct + geno_pred_idx], nm_classes, nm_sample_ct, class_ct, ws.class_min, ws.class_max);
              if (classification < 2) {
                continue;
              }
              const uint32_t is_nuisance = (geno_pred_idx < nuisance_ct);
              if (is_nuisance || (!is_score_test)) {
                sep_allele_idx = (classification == 2)? geno_pred_alleles[geno_pred_idx] : omitted_allele_idx;
                break;
              }
              coefs_unavailable = 1;
            }
            if ((sep_allele_idx == UINT32_MAX) && is_additive && (allele_ct > 2)) {
              // With an intercept, the omitted allele's count is in the span
              // of the tested ones, so it can separate too: e.g. a level whose
              // samples all carry two non-omitted alleles.
              const uint32_t classification = ClassifyPredictor(&(allele_cols[omitted_allele_idx * sample_ctav]), nm_classes, nm_sample_ct, class_ct, ws.class_min, ws.class_max);
              if (classification >= 2) {
                if (!is_score_test) {
                  sep_allele_idx = omitted_allele_idx;
                } else {
                  coefs_unavailable = 1;
                }
              }
            }
            if ((sep_allele_idx == UINT32_MAX) && (model_col_ct == 2)) {
              // (ADD, DOMDEV) and (HOM, HET) span the same space, and a
              // separating direction can hide in either basis; e.g. a level
              // with no A1 homozygote separates on HOM alone.  Check the
              // column the model doesn't include too: HOM for 'genotypic', ADD
              // for 'hethom'.
              const double* a1_col = &(allele_cols[row_a1_allele_idx * sample_ctav]);
              double* spare_col = &(ws.geno_preds[geno_pred_max * sample_ctav]);
              const uint32_t spare_is_hom = !(glm_flags & kfGlmHethom);
              for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
                const double cur_val = a1_col[nm_idx];
                spare_col[nm_idx] = spare_is_hom? ((cur_val < 1.0)? 0.0 : (cur_val - 1.0)) : cur_val;
              }
              const uint32_t classification = ClassifyPredictor(spare_col, nm_classes, nm_sample_ct, class_ct, ws.class_min, ws.class_max);
              if (classification >= 2) {
                if (!is_score_test) {
                  sep_allele_idx = (classification == 2)? row_a1_allele_idx : omitted_allele_idx;
                } else {
                  coefs_unavailable = 1;
                }
              }
            }
            if (sep_allele_idx != UINT32_MAX) {
              if (firth_mode) {
                goto GlmMultinomialThread_firth;
              }
              glm_err = SetGlmErr1(kGlmErrcodeSeparation, sep_allele_idx);
              goto GlmMultinomialThread_skip_regression;
            }
          }
          if (!always_firth) {
            // Null model, fitted to the samples with a call.
            if ((!missing_ct) && (!nuisance_ct)) {
              // same samples and predictors as the precomputed fit
              memcpy(null_coefs, setp->null_coefs, nonref_class_ct * S_CAST(uintptr_t, base_pred_ct) * sizeof(double));
              null_ln_lik = setp->null_ln_lik;
            } else {
              RemapNullCoefs(setp->null_coefs, setp->null_level_to_class, class_levels, class_ct, base_pred_ct, alt_coefs);
              for (uint32_t class_idx = 0; class_idx != nonref_class_ct; ++class_idx) {
                double* null_row = &(null_coefs[class_idx * null_pred_ct]);
                memcpy(null_row, &(alt_coefs[class_idx * base_pred_ct]), base_pred_ct * sizeof(double));
                ZeroDArr(nuisance_ct, &(null_row[base_pred_ct]));
              }
              if (MultinomialFit(pred_rows, nm_classes, nm_sample_ct, null_pred_ct, nonref_class_ct, kMultinomialMaxIter, null_coefs, &null_ln_lik, &is_unfinished, nullptr, nullptr, ws.fit_wkspace, ws.chunk_buf, ws.mi_buf, ws.dbl_2d_buf) || is_unfinished) {
                if (firth_mode) {
                  goto GlmMultinomialThread_firth;
                }
                glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
                goto GlmMultinomialThread_skip_regression;
              }
            }
            // Full model, starting from (null estimate, 0); its first Newton
            // step is the score step.
            for (uint32_t class_idx = 0; class_idx != nonref_class_ct; ++class_idx) {
              double* alt_row = &(alt_coefs[class_idx * alt_pred_ct]);
              memcpy(alt_row, &(null_coefs[class_idx * null_pred_ct]), null_pred_ct * sizeof(double));
              ZeroDArr(tested_ct, &(alt_row[null_pred_ct]));
            }
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
                  auxp->is_unfinished = 1;
                } else if (ws.cov[0] < 0.0) {
                  coefs_unavailable = 1;
                }
              }
            } else {
              if (MultinomialFit(pred_rows, nm_classes, nm_sample_ct, alt_pred_ct, nonref_class_ct, kMultinomialMaxIter, alt_coefs, &alt_ln_lik, &is_unfinished, nullptr, need_cov? ws.cov : nullptr, ws.fit_wkspace, ws.chunk_buf, ws.mi_buf, ws.dbl_2d_buf)) {
                if (firth_mode) {
                  goto GlmMultinomialThread_firth;
                }
                glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
                goto GlmMultinomialThread_skip_regression;
              }
              if (is_unfinished) {
                // Nearly always quasi-separation that the per-column checks
                // did not reveal.  The likelihood is still rising, so neither
                // statistic is trustworthy.
                if (firth_mode) {
                  goto GlmMultinomialThread_firth;
                }
                auxp->is_unfinished = 1;
                goto GlmMultinomialThread_skip_regression;
              }
              if (test_type == kGlmMultinomialTestLrt) {
                chisq = 2 * (alt_ln_lik - null_ln_lik);
                if (chisq < 0.0) {
                  // rounding
                  chisq = 0.0;
                }
              } else if ((ws.cov[0] < 0.0) || MultinomialWaldChisq(alt_coefs, ws.cov, nonref_class_ct, alt_pred_ct, null_pred_ct, tested_ct, ws.fit_wkspace, ws.mi_buf, ws.dbl_2d_buf, &chisq)) {
                if (firth_mode) {
                  goto GlmMultinomialThread_firth;
                }
                glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                goto GlmMultinomialThread_skip_regression;
              }
            }
            if (!isfinite(chisq)) {
              if (firth_mode) {
                goto GlmMultinomialThread_firth;
              }
              glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
              goto GlmMultinomialThread_skip_regression;
            }
            auxp->chisq = chisq;
          } else {
          GlmMultinomialThread_firth:
            // Firth regression: the penalized likelihood-ratio test of
            // logistf, whose restricted fit maximizes the same penalized
            // likelihood (the penalty involving the tested columns) with the
            // tested coefficients held at zero, or the Wald test with the
            // covariance described above MultinomialFirthFit().
            auxp->is_firth = 1;
            auxp->is_unfinished = 0;
            MultinomialFirthBufs fb;
            MultinomialFirthBufsLayout(nm_sample_ct, alt_pred_ct, nonref_class_ct, ws.firth_buf, &fb);
            double* firth_xx = ws.firth_xx;
            for (uint32_t pred_idx = 0; pred_idx != alt_pred_ct; ++pred_idx) {
              double* xx_row = &(firth_xx[pred_idx * nm_sample_ctav]);
              memcpy(xx_row, pred_rows[pred_idx], nm_sample_ct * sizeof(double));
              ZeroDArr(nm_sample_ctav - nm_sample_ct, &(xx_row[nm_sample_ct]));
            }
            // Start from the covariate-only fit, with zero nuisance and tested
            // coefficients.
            RemapNullCoefs(setp->null_coefs, setp->null_level_to_class, class_levels, class_ct, base_pred_ct, null_coefs);
            for (uint32_t class_idx = 0; class_idx != nonref_class_ct; ++class_idx) {
              double* alt_row = &(alt_coefs[class_idx * alt_pred_ct]);
              memcpy(alt_row, &(null_coefs[class_idx * base_pred_ct]), base_pred_ct * sizeof(double));
              ZeroDArr(alt_pred_ct - base_pred_ct, &(alt_row[base_pred_ct]));
            }
            double logdet0 = 0.0;
            uint32_t null_is_unfinished = 0;
            if (test_type == kGlmMultinomialTestLrt) {
              // restricted fit; the full fit then starts from its estimate
              if (MultinomialFirthFit(firth_xx, nm_classes, nm_sample_ct, alt_pred_ct, nonref_class_ct, null_pred_ct, alt_coefs, ws.firth_lli0, &logdet0, nullptr, &fb, &null_is_unfinished)) {
                glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
                goto GlmMultinomialThread_skip_regression;
              }
            }
            double logdet1;
            if (MultinomialFirthFit(firth_xx, nm_classes, nm_sample_ct, alt_pred_ct, nonref_class_ct, alt_pred_ct, alt_coefs, ws.firth_lli1, &logdet1, need_cov? ws.cov : nullptr, &fb, &is_unfinished)) {
              glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
              goto GlmMultinomialThread_skip_regression;
            }
            // As with --glm firth, the statistics of a fit that hit the
            // iteration limit are still reported, with error code UNFINISHED.
            auxp->is_unfinished = null_is_unfinished || is_unfinished;
            double chisq = 0.0;
            if (test_type == kGlmMultinomialTestLrt) {
              // Summed from per-sample differences, so that small statistics
              // keep their relative precision.
              const double* lli0 = ws.firth_lli0;
              const double* lli1 = ws.firth_lli1;
              double lli_diff_sum = 0.0;
              for (uint32_t nm_idx = 0; nm_idx != nm_sample_ct; ++nm_idx) {
                lli_diff_sum += lli1[nm_idx] - lli0[nm_idx];
              }
              chisq = 2 * lli_diff_sum + (logdet1 - logdet0);
              if (chisq < 0.0) {
                chisq = 0.0;
              }
            } else if ((ws.cov[0] < 0.0) || MultinomialWaldChisq(alt_coefs, ws.cov, nonref_class_ct, alt_pred_ct, null_pred_ct, tested_ct, ws.fit_wkspace, ws.mi_buf, ws.dbl_2d_buf, &chisq)) {
              glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
              goto GlmMultinomialThread_skip_regression;
            }
            if (!isfinite(chisq)) {
              glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
              goto GlmMultinomialThread_skip_regression;
            }
            auxp->chisq = chisq;
          }
          if (save_coefs && (!coefs_unavailable) && (ws.cov[0] >= 0.0) && (!class_levels[0])) {
            // Coefficients are only comparable across variants when they are
            // relative to the requested reference level.
            const double* cov = ws.cov;
            for (uint32_t class_idx = 1; class_idx != class_ct; ++class_idx) {
              double* level_beta_se = &(row_beta_se[(class_levels[class_idx] - 1) * 2 * S_CAST(uintptr_t, max_tested_ct)]);
              for (uint32_t tested_idx = 0; tested_idx != tested_ct; ++tested_idx) {
                const uintptr_t param_idx = (class_idx - 1) * alt_pred_ct + null_pred_ct + tested_idx;
                double* dst = &(level_beta_se[2 * geno_pred_slots[nuisance_ct + tested_idx]]);
                dst[0] = alt_coefs[param_idx];
                dst[1] = sqrt(cov[param_idx * alt_param_ct + param_idx]);
              }
            }
          }
          while (0) {
          GlmMultinomialThread_skip_regression:
            // is_unfinished may be set instead, with glm_err still zero
            memcpy(&(auxp->glm_err), &glm_err, 8);
            auxp->chisq = -9.0;
          }
        }
        block_aux_iter = &(block_aux_iter[max_row_ct]);
        a1_dosage_iter = &(a1_dosage_iter[max_row_ct * max_a1_ct]);
        level_a1_iter = &(level_a1_iter[max_row_ct * S_CAST(uintptr_t, level_ct) * max_a1_ct]);
        level_allele_obs_iter = &(level_allele_obs_iter[level_ct]);
        beta_se_iter = &(beta_se_iter[max_row_ct * row_beta_se_stride]);
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
    const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
    const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
    const uint32_t is_additive = !(glm_flags & (kfGlmDominant | kfGlmRecessive | kfGlmHetonly | kfGlmGenotypic | kfGlmHethom));
    const uint32_t model_col_ct = 1 + ((glm_flags & (kfGlmGenotypic | kfGlmHethom)) != 0);
    const uint32_t max_row_ct = is_additive? 1 : (max_extra_allele_ct + 1);
    const uint32_t max_a1_ct = is_additive? (max_extra_allele_ct + 1) : 1;
    const uint32_t max_tested_ct = is_additive? (max_extra_allele_ct + 1) : model_col_ct;
    const uintptr_t row_beta_se_stride = 2 * S_CAST(uintptr_t, nonref_level_ct) * max_tested_ct;
    ctx->is_additive = is_additive;
    ctx->model_col_ct = model_col_ct;
    ctx->max_row_ct = max_row_ct;
    ctx->max_a1_ct = max_a1_ct;
    ctx->max_tested_ct = max_tested_ct;
    const char* test_name = "ADD";
    if (glm_flags & kfGlmDominant) {
      test_name = "DOM";
    } else if (glm_flags & kfGlmRecessive) {
      test_name = "REC";
    } else if (glm_flags & kfGlmHetonly) {
      test_name = "HET";
    } else if (model_col_ct == 2) {
      test_name = "GENO_2DF";
    }
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
    const uint32_t max_allele_ct = max_extra_allele_ct + 2;
    const uint32_t geno_pred_max = MultinomialGenoPredMax(ctx);
    uintptr_t workspace_alloc = MultinomialWorkspaceLayout(sample_ct, ctx->sets[0].covar_ct, level_ct, max_allele_ct, geno_pred_max, ctx->firth_mode != 0, nullptr, nullptr);
    if (sample_ct_x) {
      const uintptr_t workspace_alloc_x = MultinomialWorkspaceLayout(sample_ct_x, ctx->sets[1].covar_ct, level_ct, max_allele_ct, geno_pred_max, ctx->firth_mode != 0, nullptr, nullptr);
      if (workspace_alloc_x > workspace_alloc) {
        workspace_alloc = workspace_alloc_x;
      }
    }
    if (sample_ct_y) {
      const uintptr_t workspace_alloc_y = MultinomialWorkspaceLayout(sample_ct_y, ctx->sets[2].covar_ct, level_ct, max_allele_ct, geno_pred_max, ctx->firth_mode != 0, nullptr, nullptr);
      if (workspace_alloc_y > workspace_alloc) {
        workspace_alloc = workspace_alloc_y;
      }
    }
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    // +1 is for top-level common->workspace_bufs
    const uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;
    const uintptr_t per_variant_xalloc_byte_ct = max_row_ct * (sizeof(MultinomialAuxResult) + max_a1_ct * (1 + level_ct) * sizeof(double) + row_beta_se_stride * sizeof(double)) + level_ct * sizeof(int32_t);
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    common->thread_mhc = nullptr;
    common->dosage_presents = nullptr;
    common->dosage_mains = nullptr;
    uint32_t read_block_size;
    uintptr_t max_alt_allele_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, max_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, 0, pgfip, &calc_thread_ct, &common->genovecs, max_extra_allele_ct? (&common->thread_mhc) : nullptr, nullptr, nullptr, dosage_is_present? (&common->dosage_presents) : nullptr, dosage_is_present? (&common->dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, &max_alt_allele_block_size, main_loadbufs, &common->pgr_ptrs, &common->read_variant_uidx_starts))) {
      goto GlmMultinomial_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto GlmMultinomial_ret_NOMEM;
    }
    MultinomialAuxResult* block_aux_bufs[2];
    double* block_a1_dosage_bufs[2];
    double* block_level_a1_bufs[2];
    uint32_t* block_level_allele_obs_bufs[2];
    double* block_beta_se_bufs[2];
    const uintptr_t block_row_ct = read_block_size * S_CAST(uintptr_t, max_row_ct);
    for (uint32_t uii = 0; uii != 2; ++uii) {
      if (unlikely(BIGSTACK_ALLOC_X(MultinomialAuxResult, block_row_ct, &(block_aux_bufs[uii])) ||
                   bigstack_alloc_d(block_row_ct * max_a1_ct, &(block_a1_dosage_bufs[uii])) ||
                   bigstack_alloc_d(block_row_ct * level_ct * max_a1_ct, &(block_level_a1_bufs[uii])) ||
                   bigstack_alloc_u32(read_block_size * S_CAST(uintptr_t, level_ct), &(block_level_allele_obs_bufs[uii])) ||
                   bigstack_alloc_d(block_row_ct * row_beta_se_stride, &(block_beta_se_bufs[uii])))) {
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
    const uint32_t firth_yn_col = (glm_cols & kfGlmColFirthYn) && (ctx->firth_mode == 1);
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
    if (firth_yn_col) {
      cswritep = strcpya_k(cswritep, "\tFIRTH?");
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
    const char* regression_str = "logistic";
    if (ctx->firth_mode == 1) {
      regression_str = "logistic/Firth";
    } else if (ctx->firth_mode == 2) {
      regression_str = "Firth";
    }
    logprintfww5("--glm multinomial %s regression (%s, %u levels) on phenotype '%s': ", regression_str, test_type_str, level_ct, cur_pheno_name);
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
        ctx->block_a1_dosage = block_a1_dosage_bufs[parity];
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
        const MultinomialAuxResult* block_aux = block_aux_bufs[parity];
        const double* block_a1_dosage = block_a1_dosage_bufs[parity];
        const double* block_level_a1 = block_level_a1_bufs[parity];
        const uint32_t* block_level_allele_obs = block_level_allele_obs_bufs[parity];
        const double* block_beta_se = block_beta_se_bufs[parity];
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
          uintptr_t allele_idx_offset_base = 2 * S_CAST(uintptr_t, write_variant_uidx);
          uint32_t allele_ct = 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offset_base;
          }
          const uint32_t omitted_allele_idx = omitted_alleles? omitted_alleles[write_variant_uidx] : 0;
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          const uint32_t row_ct = is_additive? 1 : (allele_ct - 1);
          const uint32_t row_a1_ct = is_additive? (allele_ct - 1) : 1;
          const uint32_t row_tested_ct = is_additive? (allele_ct - 1) : model_col_ct;
          const uint32_t* level_allele_obs = &(block_level_allele_obs[variant_bidx * S_CAST(uintptr_t, level_ct)]);
          uint32_t variant_is_valid = 0;
          // in the additive model, the first non-omitted allele
          uint32_t row_a1_allele_idx = (omitted_allele_idx == 0);
          for (uint32_t row_idx = 0; row_idx != row_ct; ++row_idx) {
            if (row_idx) {
              ++row_a1_allele_idx;
              row_a1_allele_idx += (row_a1_allele_idx == omitted_allele_idx);
            }
            const uintptr_t row_uidx = variant_bidx * S_CAST(uintptr_t, max_row_ct) + row_idx;
            const MultinomialAuxResult* auxp = &(block_aux[row_uidx]);
            const double* a1_dosages = &(block_a1_dosage[row_uidx * max_a1_ct]);
            const double* level_a1 = &(block_level_a1[row_uidx * level_ct * max_a1_ct]);
            const double* beta_se = &(block_beta_se[row_uidx * row_beta_se_stride]);
            const uint32_t is_valid = (auxp->chisq != -9.0);
            double ln_pval = kLnPvalError;
            if (is_valid) {
              variant_is_valid = 1;
              ln_pval = ChisqToLnP(auxp->chisq, auxp->df);
              if (orig_ln_pvals) {
                orig_ln_pvals[valid_allele_ct] = ln_pval;
              }
              ++valid_allele_ct;
              if (valid_alleles) {
                SetBit(allele_idx_offset_base + row_a1_allele_idx, valid_alleles);
              }
            }
            if ((ln_pfilter <= 0.0) && ((!is_valid) || (ln_pval > ln_pfilter))) {
              continue;
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
                  goto GlmMultinomial_ret_WRITE_FAIL;
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
            // A1: every non-omitted allele in the additive model
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if ((allele_idx == omitted_allele_idx) || ((!is_additive) && (allele_idx != row_a1_allele_idx))) {
                continue;
              }
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto GlmMultinomial_ret_WRITE_FAIL;
              }
              cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
            }
            --cswritep;
            if (omitted_col) {
              *cswritep++ = '\t';
              cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
            }
            if (ax_col) {
              *cswritep++ = '\t';
              if (is_additive) {
                cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
              } else {
                for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                  if (allele_idx == row_a1_allele_idx) {
                    continue;
                  }
                  if (unlikely(Cswrite(&css, &cswritep))) {
                    goto GlmMultinomial_ret_WRITE_FAIL;
                  }
                  cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                }
                --cswritep;
              }
            }
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto GlmMultinomial_ret_WRITE_FAIL;
            }
            // Per-A1 values are comma-separated lists, in A1 order.
            if (a1_ct_col) {
              *cswritep++ = '\t';
              for (uint32_t a1_idx = 0; a1_idx != row_a1_ct; ++a1_idx) {
                cswritep = dtoa_g(a1_dosages[a1_idx], cswritep);
                *cswritep++ = ',';
              }
              --cswritep;
            }
            if (tot_allele_col) {
              *cswritep++ = '\t';
              cswritep = u32toa(auxp->allele_obs_ct, cswritep);
            }
            if (a1_ct_level_col) {
              for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
                *cswritep++ = '\t';
                for (uint32_t a1_idx = 0; a1_idx != row_a1_ct; ++a1_idx) {
                  cswritep = dtoa_g(level_a1[level_idx * max_a1_ct + a1_idx], cswritep);
                  *cswritep++ = ',';
                }
                --cswritep;
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto GlmMultinomial_ret_WRITE_FAIL;
                }
              }
            }
            if (tot_allele_level_col) {
              for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
                *cswritep++ = '\t';
                cswritep = u32toa(level_allele_obs[level_idx], cswritep);
              }
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto GlmMultinomial_ret_WRITE_FAIL;
              }
            }
            if (a1_freq_col) {
              *cswritep++ = '\t';
              if (auxp->allele_obs_ct) {
                for (uint32_t a1_idx = 0; a1_idx != row_a1_ct; ++a1_idx) {
                  cswritep = dtoa_g(a1_dosages[a1_idx] / u31tod(auxp->allele_obs_ct), cswritep);
                  *cswritep++ = ',';
                }
                --cswritep;
              } else {
                cswritep = strcpya_k(cswritep, "NA");
              }
            }
            if (a1_freq_level_col) {
              for (uint32_t level_idx = 0; level_idx != level_ct; ++level_idx) {
                *cswritep++ = '\t';
                const uint32_t cur_allele_obs = level_allele_obs[level_idx];
                if (cur_allele_obs) {
                  for (uint32_t a1_idx = 0; a1_idx != row_a1_ct; ++a1_idx) {
                    cswritep = dtoa_g(level_a1[level_idx * max_a1_ct + a1_idx] / u31tod(cur_allele_obs), cswritep);
                    *cswritep++ = ',';
                  }
                  --cswritep;
                } else {
                  cswritep = strcpya_k(cswritep, "NA");
                }
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto GlmMultinomial_ret_WRITE_FAIL;
                }
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
            if (firth_yn_col) {
              *cswritep++ = '\t';
              // 'Y' - 'N' = 11
              *cswritep++ = 'N' + 11 * auxp->is_firth;
            }
            if (test_col) {
              *cswritep++ = '\t';
              cswritep = strcpya(cswritep, test_name);
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
              // one entry per tested column: A1 alleles in the additive
              // model, (ADD, DOMDEV) or (HOM, HET) for genotypic/hethom
              for (uint32_t nonref_level_idx = 0; nonref_level_idx != nonref_level_ct; ++nonref_level_idx) {
                const double* level_beta_se = &(beta_se[nonref_level_idx * 2 * S_CAST(uintptr_t, max_tested_ct)]);
                for (uint32_t col_type = 0; col_type != 4; ++col_type) {
                  if (((col_type == 1) && (!se_col)) || ((col_type >= 2) && (!ci_col))) {
                    continue;
                  }
                  *cswritep++ = '\t';
                  for (uint32_t tested_idx = 0; tested_idx != row_tested_ct; ++tested_idx) {
                    const double beta = level_beta_se[2 * tested_idx];
                    const double se = level_beta_se[2 * tested_idx + 1];
                    if (se == -9.0) {
                      cswritep = strcpya_k(cswritep, "NA");
                    } else if (col_type == 0) {
                      cswritep = dtoa_g(beta, cswritep);
                    } else if (col_type == 1) {
                      cswritep = dtoa_g(se, cswritep);
                    } else {
                      const double ci_radius = ci_zt * se;
                      cswritep = dtoa_g((col_type == 2)? (beta - ci_radius) : (beta + ci_radius), cswritep);
                    }
                    *cswritep++ = ',';
                  }
                  --cswritep;
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
