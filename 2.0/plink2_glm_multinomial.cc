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

const char* MnlRefCatname(const char* mnl_ref_flattened, const char* pheno_name) {
  if (!mnl_ref_flattened) {
    return nullptr;
  }
  // Each entry is split at its first '=', as in the --mnl-ref parser and
  // GlmMain(); a phenotype whose name contains '=' can't be named.
  const uint32_t pheno_name_slen = strlen(pheno_name);
  for (const char* entry_iter = mnl_ref_flattened; *entry_iter; ) {
    const char* eq_ptr = strchr(entry_iter, '=');
    if ((S_CAST(uintptr_t, eq_ptr - entry_iter) == pheno_name_slen) && memequal(entry_iter, pheno_name, pheno_name_slen)) {
      return &(eq_ptr[1]);
    }
    entry_iter = &(eq_ptr[strlen(eq_ptr) + 1]);
  }
  return nullptr;
}

BoolErr MnlCountCats(const uintptr_t* sample_include, const PhenoCol* pheno_col, uint32_t sample_ct, uint32_t ref_cat_idx, uint32_t* cat_ct_ptr) {
  uintptr_t* observed_cats;
  if (unlikely(bigstack_alloc_w(1 + (pheno_col->nonnull_category_ct / kBitsPerWord), &observed_cats))) {
    return 1;
  }
  const uint32_t cat_ct = IdentifyRemainingCats(sample_include, pheno_col, sample_ct, observed_cats);
  *cat_ct_ptr = IsSet(observed_cats, ref_cat_idx)? cat_ct : 0;
  BigstackReset(observed_cats);
  return 0;
}

// Newton-Raphson iteration limits.  Step-halving is only expected in the
// first iteration or two of a fit started far from the optimum; the
// log-likelihood is concave, so a long run of halvings means the fit is
// numerically stuck.
CONSTI32(kMnlMaxIter, 50);
CONSTI32(kMnlMaxHalvings, 30);

// The information matrix is accumulated over blocks of this many samples, so
// the P_c-scaled copy of the predictors stays small (and cache-resident)
// however many samples there are.
CONSTI32(kMnlSampleBlockSize, 2048);

// Scratch space for MultinomialRegressionD().  With p = predictor_ct, J =
// nonref_cat_ct, and n rounded up to a vector boundary (sample_ctav):
//   prob: J * n doubles
//   max_eta, denom: n doubles each
//   zz: (J + 1) * p * kMnlSampleBlockSize doubles
//   ss: ((J + 1) * p)^2 doubles
//   xty, grad, step, coef_old: J * p doubles each
//   info, chol: (J * p)^2 doubles each
typedef struct MnlSolverBufsStruct {
  double* prob;
  double* max_eta;
  double* denom;
  double* zz;
  double* ss;
  double* xty;
  double* grad;
  double* step;
  double* coef_old;
  double* info;
  double* chol;
} MnlSolverBufs;

uintptr_t GetMnlSolverBufsSize(uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uintptr_t param_ct = S_CAST(uintptr_t, predictor_ct) * nonref_cat_ct;
  uintptr_t size = RoundUpPow2(nonref_cat_ct * sample_ctav * sizeof(double), kCacheline);
  size += 2 * RoundUpPow2(sample_ctav * sizeof(double), kCacheline);
  const uintptr_t stacked_ct = param_ct + predictor_ct;
  size += RoundUpPow2(stacked_ct * kMnlSampleBlockSize * sizeof(double), kCacheline);
  size += RoundUpPow2(stacked_ct * stacked_ct * sizeof(double), kCacheline);
  size += 4 * RoundUpPow2(param_ct * sizeof(double), kCacheline);
  size += 2 * RoundUpPow2(param_ct * param_ct * sizeof(double), kCacheline);
  return size;
}

void CarveMnlSolverBufs(uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct, unsigned char** arena_iterp, MnlSolverBufs* bufs) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uintptr_t param_ct = S_CAST(uintptr_t, predictor_ct) * nonref_cat_ct;
  bufs->prob = S_CAST(double*, arena_alloc_raw_rd(nonref_cat_ct * sample_ctav * sizeof(double), arena_iterp));
  bufs->max_eta = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), arena_iterp));
  bufs->denom = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), arena_iterp));
  const uintptr_t stacked_ct = param_ct + predictor_ct;
  bufs->zz = S_CAST(double*, arena_alloc_raw_rd(stacked_ct * kMnlSampleBlockSize * sizeof(double), arena_iterp));
  bufs->ss = S_CAST(double*, arena_alloc_raw_rd(stacked_ct * stacked_ct * sizeof(double), arena_iterp));
  bufs->xty = S_CAST(double*, arena_alloc_raw_rd(param_ct * sizeof(double), arena_iterp));
  bufs->grad = S_CAST(double*, arena_alloc_raw_rd(param_ct * sizeof(double), arena_iterp));
  bufs->step = S_CAST(double*, arena_alloc_raw_rd(param_ct * sizeof(double), arena_iterp));
  bufs->coef_old = S_CAST(double*, arena_alloc_raw_rd(param_ct * sizeof(double), arena_iterp));
  bufs->info = S_CAST(double*, arena_alloc_raw_rd(param_ct * param_ct * sizeof(double), arena_iterp));
  bufs->chol = S_CAST(double*, arena_alloc_raw_rd(param_ct * param_ct * sizeof(double), arena_iterp));
}

// Fills prob[] (nonref_cat_ct rows of sample_ctav) with the fitted category
// probabilities under coef[], and lli[] with each sample's log-likelihood
// contribution.  Returns the total log-likelihood, which is not finite when
// the linear predictors overflow.
// xx is predictor-major with vector-aligned rows; coef is category-major.
double MnlProbsAndLoglik(const double* xx, const double* coef, const uint32_t* cats, uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct, double* __restrict prob, double* __restrict max_eta, double* __restrict denom, double* __restrict lli) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uintptr_t sample_ct_rem = sample_ctav - sample_ct;
  // Linear predictors.  The reference category's is zero; each sample's are
  // shifted by their maximum (or zero) before exponentiating.
  ZeroDArr(sample_ctav, max_eta);
  for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
    double* cur_eta = &(prob[cat_idx * sample_ctav]);
    ColMajorMatrixVectorMultiplyStrided(xx, &(coef[cat_idx * predictor_ct]), sample_ct, sample_ctav, predictor_ct, cur_eta);
    ZeroDArr(sample_ct_rem, &(cur_eta[sample_ct]));
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      if (cur_eta[sample_idx] > max_eta[sample_idx]) {
        max_eta[sample_idx] = cur_eta[sample_idx];
      }
    }
  }
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uint32_t cur_cat = cats[sample_idx];
    const double cur_eta = cur_cat? prob[(cur_cat - 1) * sample_ctav + sample_idx] : 0.0;
    lli[sample_idx] = cur_eta - max_eta[sample_idx];
    denom[sample_idx] = -max_eta[sample_idx];
  }
  ZeroDArr(sample_ct_rem, &(denom[sample_ct]));
  expd_v(denom, sample_ctav);
  for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
    double* cur_prob = &(prob[cat_idx * sample_ctav]);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      cur_prob[sample_idx] -= max_eta[sample_idx];
    }
    expd_v(cur_prob, sample_ctav);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      denom[sample_idx] += cur_prob[sample_idx];
    }
  }
  double loglik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const double cur_lli = lli[sample_idx] - log(denom[sample_idx]);
    lli[sample_idx] = cur_lli;
    loglik += cur_lli;
    // denom[] becomes its reciprocal
    denom[sample_idx] = 1.0 / denom[sample_idx];
  }
  for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
    double* cur_prob = &(prob[cat_idx * sample_ctav]);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      cur_prob[sample_idx] *= denom[sample_idx];
    }
  }
  return loglik;
}

// Lower triangle of gg (dim x dim, row-major) := beta * gg + zz zz^T, where zz
// has dim rows of col_ct elements at the given stride.
void MnlSyrk(const double* zz, uint32_t dim, uint32_t col_ct, uint32_t stride, double beta, double* gg) {
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

// Computes the log-likelihood gradient and the lower triangle of the
// information matrix (negated Hessian) at the probabilities in prob[].
// If sqrt_wts is not nullptr, sample i's contribution to the information
// matrix is weighted by sqrt_wts[i]^2, and the gradient is not computed (grad
// and xty may then be nullptr).
// Parameters are ordered category-major, so the information matrix's (c, d)
// block, predictor_ct square, is
//   X^T diag(P_c (1[c = d] - P_d)) X  =  1[c = d] D_c - Z_c^T Z_d,
// where Z_c = diag(P_c) X and D_c = Z_c^T X.  A block of samples at a time,
// X and all the Z_c are stacked into one matrix, whose self-product (one
// syrk) contains every D_c and every Z_c^T Z_d.
// The gradient's block c is
//   X^T (Y_c - P_c)  =  xty_c - (column 0 of D_c),
// since predictor 0 is the intercept.  xty must hold X^T Y (category-major).
void MnlGradAndInfo(const double* xx, const double* prob, const double* xty, const double* sqrt_wts, uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct, double* __restrict zz, double* __restrict ss, double* __restrict grad, double* __restrict info) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uint32_t param_ct = predictor_ct * nonref_cat_ct;
  const uint32_t stacked_ct = param_ct + predictor_ct;
  for (uint32_t block_start = 0; block_start < sample_ct; block_start += kMnlSampleBlockSize) {
    uint32_t block_len = sample_ct - block_start;
    if (block_len > kMnlSampleBlockSize) {
      block_len = kMnlSampleBlockSize;
    }
    double* zz_row = zz;
    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      const double* xx_row = &(xx[pred_idx * sample_ctav + block_start]);
      if (!sqrt_wts) {
        memcpy(zz_row, xx_row, block_len * sizeof(double));
      } else {
        const double* cur_sqrt_wts = &(sqrt_wts[block_start]);
        for (uint32_t uii = 0; uii != block_len; ++uii) {
          zz_row[uii] = cur_sqrt_wts[uii] * xx_row[uii];
        }
      }
      zz_row = &(zz_row[kMnlSampleBlockSize]);
    }
    for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
      const double* cur_prob = &(prob[cat_idx * sample_ctav + block_start]);
      for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
        const double* xx_row = &(zz[pred_idx * kMnlSampleBlockSize]);
        for (uint32_t uii = 0; uii != block_len; ++uii) {
          zz_row[uii] = cur_prob[uii] * xx_row[uii];
        }
        zz_row = &(zz_row[kMnlSampleBlockSize]);
      }
    }
    MnlSyrk(zz, stacked_ct, block_len, kMnlSampleBlockSize, block_start? 1.0 : 0.0, ss);
  }
  // ss (row-major lower triangle): row predictor_ct + r holds D in its first
  // predictor_ct columns and Z^T Z after that.
  for (uint32_t row_idx = 0; row_idx != param_ct; ++row_idx) {
    const double* ss_row = &(ss[(predictor_ct + row_idx) * stacked_ct]);
    const double* gg_row = &(ss_row[predictor_ct]);
    if (!sqrt_wts) {
      grad[row_idx] = xty[row_idx] - ss_row[0];
    }
    double* info_row = &(info[row_idx * param_ct]);
    const uint32_t block_start_col = row_idx - (row_idx % predictor_ct);
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
BoolErr MnlCholesky(const double* aa, uint32_t dim, double* ll) {
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
void MnlCholSolve(const double* ll, const double* yy, uint32_t dim, double* __restrict xx) {
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

// Diagonal of (L L^T)^{-1}: element k is the squared norm of column k of
// L^{-1}, i.e. of the solution of L x = e_k.  x_buf must have room for dim
// doubles.
void MnlCholInvDiag(const double* ll, uint32_t dim, double* __restrict x_buf, double* __restrict inv_diag) {
  for (uint32_t col_idx = 0; col_idx != dim; ++col_idx) {
    double ssq = 0.0;
    for (uint32_t row_idx = col_idx; row_idx != dim; ++row_idx) {
      const double* ll_row = &(ll[row_idx * dim]);
      double dxx = (row_idx == col_idx)? 1.0 : 0.0;
      for (uint32_t uii = col_idx; uii != row_idx; ++uii) {
        dxx -= ll_row[uii] * x_buf[uii];
      }
      dxx /= ll_row[row_idx];
      x_buf[row_idx] = dxx;
      ssq += dxx * dxx;
    }
    inv_diag[col_idx] = ssq;
  }
}

// Fits the multinomial logit model
//   P(category c | x) = exp(x . b_c) / (1 + sum_d exp(x . b_d))
// for c = 1..nonref_cat_ct (category 0 is the reference) by Newton-Raphson
// with step-halving.  Convergence: the Newton step computed at the current
// coefficients moves none of them by more than 1e-8, and predicts a
// log-likelihood gain below 1e-10 * (0.1 + |loglik|).  That last step is then
// applied without re-evaluating anything, which leaves an error of order
// (step size)^2 in the coefficients; the log-likelihood contributions are
// those at the point the step was taken from, off by about the predicted gain.
// The standard errors come from the information matrix at that point too when
// the step is below 1e-12 (a relative error of that order), and from a fresh
// evaluation at the returned coefficients otherwise.
//
// xx: predictor-major, rows sample_ctav long, trailing elements zero; row 0
//   must be the intercept (all 1s).
// coef: starting point on input, maximum-likelihood estimate on output;
//   category-major (nonref_cat_ct rows of predictor_ct).
// lli: set to each sample's log-likelihood contribution at the estimate.
// se: if not nullptr, set to the Wald standard errors, from the inverse of the
//   information matrix at the estimate.
//
// Returns kGlmErrcodeNone, kGlmErrcodeMnlConvergeFail, or (when se is
// requested and the information matrix is numerically singular or its
// inverse has a nonpositive diagonal element) kGlmErrcodeInvalidResult.
GlmErrcode MultinomialRegressionD(const double* xx, const uint32_t* cats, uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct, double* __restrict coef, double* __restrict lli, double* __restrict se, MnlSolverBufs* bufs) {
  const uint32_t param_ct = predictor_ct * nonref_cat_ct;
  {
    // X^T Y, constant over the iterations
    const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
    double* xty = bufs->xty;
    ZeroDArr(param_ct, xty);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t cur_cat = cats[sample_idx];
      if (cur_cat) {
        double* xty_iter = &(xty[(cur_cat - 1) * predictor_ct]);
        for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
          xty_iter[pred_idx] += xx[pred_idx * sample_ctav + sample_idx];
        }
      }
    }
  }
  double loglik = MnlProbsAndLoglik(xx, coef, cats, sample_ct, predictor_ct, nonref_cat_ct, bufs->prob, bufs->max_eta, bufs->denom, lli);
  if (!isfinite_d(loglik)) {
    return kGlmErrcodeMnlConvergeFail;
  }
  // final_pass: the coefficients have converged, and the information matrix
  // is being re-evaluated there for the standard errors.
  uint32_t final_pass = 0;
  for (uint32_t iter_idx = 0; ; ++iter_idx) {
    MnlGradAndInfo(xx, bufs->prob, bufs->xty, nullptr, sample_ct, predictor_ct, nonref_cat_ct, bufs->zz, bufs->ss, bufs->grad, bufs->info);
    if (MnlCholesky(bufs->info, param_ct, bufs->chol)) {
      return final_pass? kGlmErrcodeInvalidResult : kGlmErrcodeMnlConvergeFail;
    }
    if (final_pass) {
      break;
    }
    double* step = bufs->step;
    MnlCholSolve(bufs->chol, bufs->grad, param_ct, step);
    double max_step = 0.0;
    // grad . step = grad^T info^{-1} grad, twice the predicted gain
    double twice_ll_gain = 0.0;
    for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
      const double cur_step = step[param_idx];
      twice_ll_gain += bufs->grad[param_idx] * cur_step;
      if (fabs(cur_step) > max_step) {
        max_step = fabs(cur_step);
      }
    }
    if ((max_step < 1e-8) && (twice_ll_gain < 2e-10 * (0.1 + fabs(loglik)))) {
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        coef[param_idx] += step[param_idx];
      }
      if (!se) {
        return kGlmErrcodeNone;
      }
      if (max_step <= 1e-12) {
        break;
      }
      loglik = MnlProbsAndLoglik(xx, coef, cats, sample_ct, predictor_ct, nonref_cat_ct, bufs->prob, bufs->max_eta, bufs->denom, lli);
      if (!isfinite_d(loglik)) {
        return kGlmErrcodeMnlConvergeFail;
      }
      final_pass = 1;
      continue;
    }
    if (iter_idx == kMnlMaxIter) {
      return kGlmErrcodeMnlConvergeFail;
    }
    memcpy(bufs->coef_old, coef, param_ct * sizeof(double));
    const double ll_tol = 1e-10 * (0.1 + fabs(loglik));
    double new_loglik;
    uint32_t halving_ct = 0;
    while (1) {
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        coef[param_idx] = bufs->coef_old[param_idx] + step[param_idx];
      }
      new_loglik = MnlProbsAndLoglik(xx, coef, cats, sample_ct, predictor_ct, nonref_cat_ct, bufs->prob, bufs->max_eta, bufs->denom, lli);
      // A decrease below the convergence tolerance is rounding noise near
      // the optimum, not a reason to backtrack.
      if (new_loglik > loglik - ll_tol) {
        break;
      }
      if (++halving_ct == kMnlMaxHalvings) {
        return kGlmErrcodeMnlConvergeFail;
      }
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        step[param_idx] *= 0.5;
      }
    }
    loglik = new_loglik;
  }
  // bufs->step is free now
  MnlCholInvDiag(bufs->chol, param_ct, bufs->step, se);
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
    const double cur_var = se[param_idx];
    if ((cur_var < 1e-20) || (!isfinite_d(cur_var))) {
      return kGlmErrcodeInvalidResult;
    }
    se[param_idx] = sqrt(cur_var);
  }
  return kGlmErrcodeNone;
}

// Firth-penalized fit limits.  The iteration (a quasi-Newton step on the
// modified score) converges linearly, so it needs more iterations than
// Newton-Raphson on the unpenalized likelihood, and a tighter step tolerance
// to reach the same accuracy.
CONSTI32(kMnlFirthMaxIter, 100);

// Additional scratch space for MultinomialFirthRegressionD(), with the same
// notation as MnlSolverBufs:
//   iinv, linv, info_aug: (J * p)^2 doubles each
//   ww: J * n doubles
//   sqrt_wts: n doubles
//   hab: kMnlSampleBlockSize doubles
//   ustar: J * p doubles
//   xtw: p doubles
typedef struct MnlFirthBufsStruct {
  double* iinv;
  double* linv;
  double* info_aug;
  double* ww;
  double* sqrt_wts;
  double* hab;
  double* ustar;
  double* xtw;
} MnlFirthBufs;

uintptr_t GetMnlFirthBufsSize(uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uintptr_t param_ct = S_CAST(uintptr_t, predictor_ct) * nonref_cat_ct;
  uintptr_t size = 3 * RoundUpPow2(param_ct * param_ct * sizeof(double), kCacheline);
  size += RoundUpPow2(nonref_cat_ct * sample_ctav * sizeof(double), kCacheline);
  size += RoundUpPow2(sample_ctav * sizeof(double), kCacheline);
  size += RoundUpPow2(kMnlSampleBlockSize * sizeof(double), kCacheline);
  size += RoundUpPow2(param_ct * sizeof(double), kCacheline);
  size += RoundUpPow2(predictor_ct * sizeof(double), kCacheline);
  return size;
}

void CarveMnlFirthBufs(uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct, unsigned char** arena_iterp, MnlFirthBufs* fbufs) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uintptr_t param_ct = S_CAST(uintptr_t, predictor_ct) * nonref_cat_ct;
  fbufs->iinv = S_CAST(double*, arena_alloc_raw_rd(param_ct * param_ct * sizeof(double), arena_iterp));
  fbufs->linv = S_CAST(double*, arena_alloc_raw_rd(param_ct * param_ct * sizeof(double), arena_iterp));
  fbufs->info_aug = S_CAST(double*, arena_alloc_raw_rd(param_ct * param_ct * sizeof(double), arena_iterp));
  fbufs->ww = S_CAST(double*, arena_alloc_raw_rd(nonref_cat_ct * sample_ctav * sizeof(double), arena_iterp));
  fbufs->sqrt_wts = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), arena_iterp));
  fbufs->hab = S_CAST(double*, arena_alloc_raw_rd(kMnlSampleBlockSize * sizeof(double), arena_iterp));
  fbufs->ustar = S_CAST(double*, arena_alloc_raw_rd(param_ct * sizeof(double), arena_iterp));
  fbufs->xtw = S_CAST(double*, arena_alloc_raw_rd(predictor_ct * sizeof(double), arena_iterp));
}

// Full inverse (both triangles, row-major) of L L^T, via L^{-1}: the inverse is
// L^{-T} L^{-1}.  linv is scratch space (dim x dim).
void MnlCholInvFull(const double* ll, uint32_t dim, double* __restrict linv, double* __restrict inv) {
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

// cc (row_ct x col_ct, column-major with stride ldc) := aa bb, where aa is
// row_ct x common_ct (column-major, stride lda) and bb is common_ct x col_ct
// (column-major, stride ldb).  Unlike the NOLAPACK version of
// ColMajorMatrixMultiplyStrided(), never reads cc.
void MnlGemm(const double* aa, const double* bb, uint32_t row_ct, uint32_t lda, uint32_t col_ct, uint32_t ldb, uint32_t common_ct, uint32_t ldc, double* cc) {
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

// Evaluates everything the penalized fit needs at coef[]: prob[] and lli[] (see
// MnlProbsAndLoglik()), the gradient and information matrix (see
// MnlGradAndInfo()), and the information matrix's Cholesky factor.  Sets
// *loglik_ptr and *logdet_ptr (log det of the information matrix).  Returns 1
// if the log-likelihood is not finite or the information matrix is not
// numerically positive definite.
BoolErr MnlFirthEval(const double* xx, const double* coef, const uint32_t* cats, uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct, double* __restrict lli, MnlSolverBufs* bufs, double* __restrict loglik_ptr, double* __restrict logdet_ptr) {
  const uint32_t param_ct = predictor_ct * nonref_cat_ct;
  const double loglik = MnlProbsAndLoglik(xx, coef, cats, sample_ct, predictor_ct, nonref_cat_ct, bufs->prob, bufs->max_eta, bufs->denom, lli);
  if (!isfinite_d(loglik)) {
    return 1;
  }
  MnlGradAndInfo(xx, bufs->prob, bufs->xty, nullptr, sample_ct, predictor_ct, nonref_cat_ct, bufs->zz, bufs->ss, bufs->grad, bufs->info);
  if (MnlCholesky(bufs->info, param_ct, bufs->chol)) {
    return 1;
  }
  double half_logdet = 0.0;
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
    half_logdet += log(bufs->chol[param_idx * (param_ct + 1)]);
  }
  *loglik_ptr = loglik;
  *logdet_ptr = 2 * half_logdet;
  return 0;
}

// Firth's modified score (Firth 1993; multinomial form of Bull, Mak and
// Greenwood 2002) at the point last evaluated by MnlFirthEval():
//   U*_(c,j) = U_(c,j) + 0.5 tr(I^{-1} dI/db_(c,j)).
// With W_i = diag(p_i) - p_i p_i^T and the per-sample "hat" blocks
//   h_i^(ab) = x_i^T (I^{-1})_(ab) x_i,
// the trace is sum_i x_ij sum_(a,b) [dW_i/deta_ic]_(ab) h_i^(ab), which works
// out to sum_i x_ij p_ic (r_ic - sum_a p_ia r_ia), where
//   r_ia = h_i^(aa) - 2 sum_b p_ib h_i^(ab).
// (With two categories, 0.5 p (r - p r) = v h (0.5 - p), logistf's
// correction.)  The h_i^(ab) are computed a block of samples and a category a
// at a time: G = X (I^{-1})_(a, b >= a) is one matrix product, and
// h_i^(ab) = sum_l G[i][(b - a) * p + l] x_il (h_i^(ba) = h_i^(ab)).
//
// Also computes the matrix the fitting step is solved with: the information
// matrix of the pseudo-data whose score is U* (Heinze and Schemper 2002),
// treating the hat blocks as fixed, i.e. the information matrix with sample i
// weighted by 1 + tr(W_i h_i) (1 + its leverage).  With two categories that is
// X^T diag((1 + h_i) v_i) X, FirthRegressionD()'s step matrix.  Near a
// separated solution the plain information matrix underestimates the
// curvature of l* badly enough for the step to overshoot and oscillate.
// Its lower triangle is left in fbufs->info_aug, and the full inverse of the
// information matrix in fbufs->iinv.
void MnlFirthScore(const double* xx, uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct, MnlSolverBufs* bufs, MnlFirthBufs* fbufs) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uint32_t param_ct = predictor_ct * nonref_cat_ct;
  double* iinv = fbufs->iinv;
  MnlCholInvFull(bufs->chol, param_ct, fbufs->linv, iinv);
  const double* prob = bufs->prob;
  double* ww = fbufs->ww;
  double* sqrt_wts = fbufs->sqrt_wts;
  double* hab = fbufs->hab;
  // G, block_len x param_ct with stride kMnlSampleBlockSize, fits in zz.
  double* gg = bufs->zz;
  for (uint32_t block_start = 0; block_start < sample_ct; block_start += kMnlSampleBlockSize) {
    uint32_t block_len = sample_ct - block_start;
    if (block_len > kMnlSampleBlockSize) {
      block_len = kMnlSampleBlockSize;
    }
    const double* xx_block = &(xx[block_start]);
    // tr(W_i h_i) = sum_(a,b) p_a (1[a = b] - p_b) h^(ab), accumulated here
    // before becoming sqrt(1 + tr(W_i h_i))
    double* cur_sqrt_wts = &(sqrt_wts[block_start]);
    ZeroDArr(block_len, cur_sqrt_wts);
    for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
      ZeroDArr(block_len, &(ww[cat_idx * sample_ctav + block_start]));
    }
    for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
      const double* prob_a = &(prob[cat_idx * sample_ctav + block_start]);
      double* rr_a = &(ww[cat_idx * sample_ctav + block_start]);
      // Columns b >= a of (I^{-1})_(a, *), i.e. predictor_ct rows of iinv
      // from column a * predictor_ct on, are (by symmetry) also the transpose
      // of the corresponding columns: a column-major predictor_ct x
      // ((J - a) * predictor_ct) matrix with stride param_ct.
      const uint32_t col_start = cat_idx * predictor_ct;
      MnlGemm(xx_block, &(iinv[col_start * param_ct + col_start]), block_len, sample_ctav, param_ct - col_start, param_ct, predictor_ct, kMnlSampleBlockSize, gg);
      for (uint32_t cat_idx2 = cat_idx; cat_idx2 != nonref_cat_ct; ++cat_idx2) {
        const double* gg_iter = &(gg[(cat_idx2 - cat_idx) * predictor_ct * kMnlSampleBlockSize]);
        for (uint32_t uii = 0; uii != block_len; ++uii) {
          hab[uii] = gg_iter[uii] * xx_block[uii];
        }
        for (uint32_t pred_idx = 1; pred_idx != predictor_ct; ++pred_idx) {
          const double* gg_col = &(gg_iter[pred_idx * kMnlSampleBlockSize]);
          const double* xx_col = &(xx_block[pred_idx * sample_ctav]);
          for (uint32_t uii = 0; uii != block_len; ++uii) {
            hab[uii] += gg_col[uii] * xx_col[uii];
          }
        }
        const double* prob_b = &(prob[cat_idx2 * sample_ctav + block_start]);
        if (cat_idx2 == cat_idx) {
          for (uint32_t uii = 0; uii != block_len; ++uii) {
            rr_a[uii] += hab[uii] * (1.0 - 2 * prob_a[uii]);
            cur_sqrt_wts[uii] += hab[uii] * prob_a[uii] * (1.0 - prob_a[uii]);
          }
        } else {
          // h^(ab) and h^(ba) both
          double* rr_b = &(ww[cat_idx2 * sample_ctav + block_start]);
          for (uint32_t uii = 0; uii != block_len; ++uii) {
            const double twice_hab = 2 * hab[uii];
            rr_a[uii] -= twice_hab * prob_b[uii];
            rr_b[uii] -= twice_hab * prob_a[uii];
            cur_sqrt_wts[uii] -= twice_hab * prob_a[uii] * prob_b[uii];
          }
        }
      }
    }
    for (uint32_t uii = 0; uii != block_len; ++uii) {
      // (the trace is nonnegative; this only guards against rounding)
      cur_sqrt_wts[uii] = sqrt(1.0 + MAXV(cur_sqrt_wts[uii], 0.0));
    }
    // ww_c := 0.5 p_c (r_c - sum_a p_a r_a); hab is reused for the sums
    ZeroDArr(block_len, hab);
    for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
      const double* cur_prob = &(prob[cat_idx * sample_ctav + block_start]);
      const double* rr = &(ww[cat_idx * sample_ctav + block_start]);
      for (uint32_t uii = 0; uii != block_len; ++uii) {
        hab[uii] += cur_prob[uii] * rr[uii];
      }
    }
    for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
      const double* cur_prob = &(prob[cat_idx * sample_ctav + block_start]);
      double* rr = &(ww[cat_idx * sample_ctav + block_start]);
      for (uint32_t uii = 0; uii != block_len; ++uii) {
        rr[uii] = 0.5 * cur_prob[uii] * (rr[uii] - hab[uii]);
      }
    }
  }
  double* ustar = fbufs->ustar;
  double* xtw = fbufs->xtw;
  for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
    ColMajorVectorMatrixMultiplyStrided(&(ww[cat_idx * sample_ctav]), xx, sample_ct, sample_ctav, predictor_ct, xtw);
    const double* grad_iter = &(bufs->grad[cat_idx * predictor_ct]);
    double* ustar_iter = &(ustar[cat_idx * predictor_ct]);
    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      ustar_iter[pred_idx] = grad_iter[pred_idx] + xtw[pred_idx];
    }
  }
  MnlGradAndInfo(xx, prob, nullptr, sqrt_wts, sample_ct, predictor_ct, nonref_cat_ct, bufs->zz, bufs->ss, nullptr, fbufs->info_aug);
}

// Maximizes the Firth-penalized log-likelihood
//   l*(b) = l(b) + 0.5 log det I(b),
// I being the full information matrix, by quasi-Newton steps on the modified
// score (b += I_aug^{-1} U*, I_aug as described above MnlFirthScore(); this
// is FirthRegressionD()'s step), capped at 5 per coefficient (as in logistf)
// and halved while l* decreases.  Converged when the step computed at the
// current coefficients is below 1e-10 * (1 + |coefficient|) for every
// coefficient; that step is then applied without re-evaluating anything, as
// in MultinomialRegressionD(): lli[], the log-determinant and I_aug are those
// at the point the step was taken from, which moves l* by about U* . step, far
// below anything reported.
//
// fixed_pred_idx: if not UINT32_MAX, that predictor's coefficients (one per
//   nonreference category) are held at zero, and l* (whose penalty is still
//   that of the full model) is maximized over the others: the restricted fit
//   of a penalized likelihood-ratio test, as in logistf.  The step is then
//   solved with the free coefficients' block of I_aug.
// xx, cats, coef: as in MultinomialRegressionD().
// lli: set to each sample's (unpenalized) log-likelihood contribution at the
//   estimate.
// logdet_ptr: set to log det I at the estimate, so that
//   l* = sum(lli) + 0.5 * (*logdet_ptr).
// se: if not nullptr (the full fit only), set to sqrt(diag(I_aug^{-1})).  This
//   is the covariance logistf reports (its 'var', the inverse of
//   X^T diag((1 + h) v) X) and FirthRegressionD() returns, rather than the
//   inverse of I itself.
// lli and the estimate's l* are those at the point the last step was taken
// from, as described above.
//
// Returns kGlmErrcodeNone, kGlmErrcodeFirthConvergeFail, or (se requested,
// nonpositive or tiny variance) kGlmErrcodeInvalidResult.
GlmErrcode MultinomialFirthRegressionD(const double* xx, const uint32_t* cats, uint32_t sample_ct, uint32_t predictor_ct, uint32_t nonref_cat_ct, uint32_t fixed_pred_idx, double* __restrict coef, double* __restrict lli, double* __restrict logdet_ptr, double* __restrict se, MnlSolverBufs* bufs, MnlFirthBufs* fbufs) {
  const uint32_t param_ct = predictor_ct * nonref_cat_ct;
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  {
    double* xty = bufs->xty;
    ZeroDArr(param_ct, xty);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t cur_cat = cats[sample_idx];
      if (cur_cat) {
        double* xty_iter = &(xty[(cur_cat - 1) * predictor_ct]);
        for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
          xty_iter[pred_idx] += xx[pred_idx * sample_ctav + sample_idx];
        }
      }
    }
  }
  const uint32_t is_restricted = (fixed_pred_idx != UINT32_MAX);
  if (is_restricted) {
    for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
      coef[cat_idx * predictor_ct + fixed_pred_idx] = 0.0;
    }
  }
  const uint32_t free_param_ct = is_restricted? (param_ct - nonref_cat_ct) : param_ct;
  double loglik;
  double logdet;
  if (MnlFirthEval(xx, coef, cats, sample_ct, predictor_ct, nonref_cat_ct, lli, bufs, &loglik, &logdet)) {
    return kGlmErrcodeFirthConvergeFail;
  }
  double lstar = loglik + 0.5 * logdet;
  double* step = bufs->step;
  for (uint32_t iter_idx = 0; ; ++iter_idx) {
    MnlFirthScore(xx, sample_ct, predictor_ct, nonref_cat_ct, bufs, fbufs);
    const double* ustar = fbufs->ustar;
    // bufs->chol isn't needed any more at this point: logdet has been
    // extracted, and iinv computed.
    if (!is_restricted) {
      if (MnlCholesky(fbufs->info_aug, param_ct, bufs->chol)) {
        return kGlmErrcodeFirthConvergeFail;
      }
      MnlCholSolve(bufs->chol, ustar, param_ct, step);
    } else {
      // Free coefficients' block of I_aug (lower triangle), and of U*,
      // packed.  iinv and grad aren't needed any more at this point.
      double* info_ff = fbufs->iinv;
      double* ustar_f = bufs->grad;
      uint32_t row_f = 0;
      for (uint32_t row_idx = 0; row_idx != param_ct; ++row_idx) {
        if ((row_idx % predictor_ct) == fixed_pred_idx) {
          continue;
        }
        ustar_f[row_f] = ustar[row_idx];
        const double* info_row = &(fbufs->info_aug[row_idx * param_ct]);
        double* info_ff_row = &(info_ff[row_f * free_param_ct]);
        uint32_t col_f = 0;
        for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
          if ((col_idx % predictor_ct) == fixed_pred_idx) {
            continue;
          }
          info_ff_row[col_f++] = info_row[col_idx];
        }
        ++row_f;
      }
      if (MnlCholesky(info_ff, free_param_ct, bufs->chol)) {
        return kGlmErrcodeFirthConvergeFail;
      }
      // Solve into coef_old, then scatter into step.
      MnlCholSolve(bufs->chol, ustar_f, free_param_ct, bufs->coef_old);
      row_f = free_param_ct;
      for (uint32_t param_idx = param_ct; param_idx; ) {
        --param_idx;
        step[param_idx] = ((param_idx % predictor_ct) == fixed_pred_idx)? 0.0 : bufs->coef_old[--row_f];
      }
    }
    double max_step = 0.0;
    uint32_t converged = 1;
    for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
      const double abs_step = fabs(step[param_idx]);
      if (abs_step > max_step) {
        max_step = abs_step;
      }
      if (abs_step >= 1e-10 * (1.0 + fabs(coef[param_idx]))) {
        converged = 0;
      }
    }
    if (converged) {
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        coef[param_idx] += step[param_idx];
      }
      break;
    }
    if (iter_idx == kMnlFirthMaxIter) {
      return kGlmErrcodeFirthConvergeFail;
    }
    if (max_step > 5.0) {
      const double scale = 5.0 / max_step;
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        step[param_idx] *= scale;
      }
    }
    memcpy(bufs->coef_old, coef, param_ct * sizeof(double));
    const double lstar_tol = 1e-10 * (0.1 + fabs(lstar));
    uint32_t halving_ct = 0;
    while (1) {
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        coef[param_idx] = bufs->coef_old[param_idx] + step[param_idx];
      }
      if (!MnlFirthEval(xx, coef, cats, sample_ct, predictor_ct, nonref_cat_ct, lli, bufs, &loglik, &logdet)) {
        const double new_lstar = loglik + 0.5 * logdet;
        if (new_lstar > lstar - lstar_tol) {
          lstar = new_lstar;
          break;
        }
      }
      if (++halving_ct == kMnlMaxHalvings) {
        return kGlmErrcodeFirthConvergeFail;
      }
      for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
        step[param_idx] *= 0.5;
      }
    }
  }
  *logdet_ptr = logdet;
  if (!se) {
    return kGlmErrcodeNone;
  }
  // (se is only requested for the full fit, so bufs->chol holds the Cholesky
  // factor of I_aug the last step was solved with)
  assert(!is_restricted);
  MnlCholInvDiag(bufs->chol, param_ct, step, se);
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
    const double cur_var = se[param_idx];
    if ((cur_var < 1e-20) || (!isfinite_d(cur_var))) {
      return kGlmErrcodeInvalidResult;
    }
    se[param_idx] = sqrt(cur_var);
  }
  return kGlmErrcodeNone;
}

BoolErr GlmAllocFillAndTestPhenoCovarsMnl(const uintptr_t* sample_include, const PhenoCol* pheno_col, uint32_t ref_cat_idx, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, GlmMnlSet* mnl_set_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uint32_t nonnull_cat_ct = pheno_col->nonnull_category_ct;
  uint32_t* cat_sizes;
  if (unlikely(bigstack_alloc_kcp(new_covar_ct, cur_covar_names_ptr) ||
               bigstack_alloc_d(new_covar_ct * sample_ct, &(mnl_set_ptr->covars_cmaj)) ||
               bigstack_alloc_u32(sample_ct, &(mnl_set_ptr->pheno_cats)) ||
               bigstack_alloc_u32(nonnull_cat_ct + 1, &cat_sizes))) {
    return 1;
  }
  // Collapse the categories present: the reference goes to 0, the others to
  // 1, 2, ... in (natural-sorted) category index order.  cat_sizes[] is
  // reused as the old-index -> new-index map.
  IdentifyRemainingCatsAndMostCommon(sample_include, pheno_col, sample_ct, nullptr, cat_sizes);
  assert(cat_sizes[ref_cat_idx]);
  uint32_t cat_ct = 1;
  for (uint32_t cat_idx = 1; cat_idx <= nonnull_cat_ct; ++cat_idx) {
    if (cat_sizes[cat_idx] && (cat_idx != ref_cat_idx)) {
      ++cat_ct;
    }
  }
  const uint32_t nonref_cat_ct = cat_ct - 1;
  assert(nonref_cat_ct);
  const char** cat_names;
  double* null_coefs;
  const uintptr_t null_predictor_ct = new_covar_ct + 1;
  if (unlikely(bigstack_alloc_kcp(nonref_cat_ct, &cat_names) ||
               bigstack_alloc_d(nonref_cat_ct * null_predictor_ct, &null_coefs) ||
               bigstack_alloc_d(sample_ct, &(mnl_set_ptr->null_lli)))) {
    return 1;
  }
  unsigned char* bigstack_mark = g_bigstack_base;
  // cat_sizes[] becomes the old-index -> new-index map; each entry is read
  // before it is overwritten.
  uint32_t* cat_remap = cat_sizes;
  const double ln_ref_size = log(u31tod(cat_sizes[ref_cat_idx]));
  uint32_t new_cat_idx = 1;
  for (uint32_t cat_idx = 1; cat_idx <= nonnull_cat_ct; ++cat_idx) {
    if (cat_sizes[cat_idx] && (cat_idx != ref_cat_idx)) {
      cat_names[new_cat_idx - 1] = pheno_col->category_names[cat_idx];
      // Starting point for the covariate-only fit: the intercept-only
      // maximum-likelihood estimate.
      double* cur_null_coefs = &(null_coefs[(new_cat_idx - 1) * null_predictor_ct]);
      cur_null_coefs[0] = log(u31tod(cat_sizes[cat_idx])) - ln_ref_size;
      ZeroDArr(new_covar_ct, &(cur_null_coefs[1]));
      cat_remap[cat_idx] = new_cat_idx;
      ++new_cat_idx;
    }
  }
  cat_remap[ref_cat_idx] = 0;
  {
    const uint32_t* pheno_cats_raw = pheno_col->data.cat;
    uint32_t* pheno_cats = mnl_set_ptr->pheno_cats;
    uintptr_t sample_uidx_base = 0;
    uintptr_t sample_include_bits = sample_include[0];
    for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
      pheno_cats[sample_idx] = cat_remap[pheno_cats_raw[sample_uidx]];
    }
  }
  mnl_set_ptr->cat_names = cat_names;
  mnl_set_ptr->null_coefs = null_coefs;
  mnl_set_ptr->cat_ct = cat_ct;

  double* covar_dotprod;
  double* inverse_corr_buf;
  if (unlikely(bigstack_alloc_d(new_covar_ct * new_covar_ct, &covar_dotprod) ||
               bigstack_alloc_d(new_covar_ct * new_covar_ct, &inverse_corr_buf))) {
    return 1;
  }
  PglErr reterr = GlmFillAndTestCovars(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, 0, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, max_corr, vif_thresh, covar_dotprod, nullptr, inverse_corr_buf, mnl_set_ptr->covars_cmaj, *cur_covar_names_ptr, glm_err_ptr);
  if (unlikely(reterr)) {
    BigstackReset(bigstack_mark);
    return (reterr == kPglRetNomem);
  }
  if (*glm_err_ptr) {
    BigstackReset(bigstack_mark);
    return 0;
  }
  BigstackReset(covar_dotprod);

  // Covariate-only fit on the whole sample set.
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  double* xx;
  if (unlikely(bigstack_alloc_d(null_predictor_ct * sample_ctav, &xx))) {
    return 1;
  }
  const uintptr_t bufs_size = GetMnlSolverBufsSize(sample_ct, null_predictor_ct, nonref_cat_ct);
  if (unlikely(bigstack_left() < bufs_size)) {
    return 1;
  }
  unsigned char* arena_iter = g_bigstack_base;
  MnlSolverBufs bufs;
  CarveMnlSolverBufs(sample_ct, null_predictor_ct, nonref_cat_ct, &arena_iter, &bufs);
  FillDVec(sample_ct, 1.0, xx);
  ZeroDArr(sample_ctav - sample_ct, &(xx[sample_ct]));
  for (uintptr_t covar_idx = 0; covar_idx != new_covar_ct; ++covar_idx) {
    double* xx_row = &(xx[(covar_idx + 1) * sample_ctav]);
    memcpy(xx_row, &(mnl_set_ptr->covars_cmaj[covar_idx * sample_ct]), sample_ct * sizeof(double));
    ZeroDArr(sample_ctav - sample_ct, &(xx_row[sample_ct]));
  }
  const GlmErrcode errcode = MultinomialRegressionD(xx, mnl_set_ptr->pheno_cats, sample_ct, null_predictor_ct, nonref_cat_ct, null_coefs, mnl_set_ptr->null_lli, nullptr, &bufs);
  if (errcode) {
    *glm_err_ptr = SetGlmErr0(kGlmErrcodeMnlConvergeFail);
  }
  BigstackReset(bigstack_mark);
  return 0;
}

// Per-thread workspace.  With p = covar_ct + 2 predictors (intercept,
// covariates, genotype, in that order, so that the covariate-only model is
// the first p - 1 rows), J = cat_ct - 1 and n = sample_ct:
//   sample_nm, tmp_nm, cached_nm: n bits each
//   nm_cats: n uint32s
//   xx: p * n doubles
//   lli, cached_null_lli: n doubles each
//   coef, se: J * p doubles each
//   cached_null_coefs: J * (p - 1) doubles
//   cat_dosage_sums: cat_ct doubles
//   predictor_dotprod_buf: (p - 1)^2 doubles
//   dbl_2d_buf: p * max(p, 7) doubles
//   inverse_corr_buf: p * max(p, 3) doubles
//   inv_1d_buf: p MatrixInvertBuf1 slots
//   solver buffers (full model)
// and, when Firth regression may be used:
//   firth_null_coef: J * p doubles
//   firth_null_lli: n doubles
//   Firth solver buffers
uintptr_t GetMnlWorkspaceSize(uint32_t sample_ct, uint32_t covar_ct, uint32_t cat_ct, uint32_t is_sometimes_firth) {
  const uintptr_t predictor_ct = covar_ct + 2;
  const uintptr_t nonref_cat_ct = cat_ct - 1;
  const uintptr_t param_ct = predictor_ct * nonref_cat_ct;
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  uintptr_t workspace_size = 3 * RoundUpPow2(BitCtToWordCt(sample_ct) * sizeof(intptr_t), kCacheline);
  workspace_size += RoundUpPow2(sample_ct * sizeof(int32_t), kCacheline);
  workspace_size += RoundUpPow2(predictor_ct * sample_ctav * sizeof(double), kCacheline);
  workspace_size += 2 * RoundUpPow2(sample_ctav * sizeof(double), kCacheline);
  workspace_size += 2 * RoundUpPow2(param_ct * sizeof(double), kCacheline);
  workspace_size += RoundUpPow2((predictor_ct - 1) * nonref_cat_ct * sizeof(double), kCacheline);
  workspace_size += RoundUpPow2(cat_ct * sizeof(double), kCacheline);
  workspace_size += RoundUpPow2((predictor_ct - 1) * (predictor_ct - 1) * sizeof(double), kCacheline);
  workspace_size += RoundUpPow2(predictor_ct * MAXV(predictor_ct, 7) * sizeof(double), kCacheline);
  workspace_size += RoundUpPow2(predictor_ct * MAXV(predictor_ct, 3) * sizeof(double), kCacheline);
  workspace_size += RoundUpPow2(predictor_ct * kMatrixInvertBuf1CheckedAlloc, kCacheline);
  workspace_size += GetMnlSolverBufsSize(sample_ct, predictor_ct, nonref_cat_ct);
  if (is_sometimes_firth) {
    workspace_size += RoundUpPow2(param_ct * sizeof(double), kCacheline);
    workspace_size += RoundUpPow2(sample_ctav * sizeof(double), kCacheline);
    workspace_size += GetMnlFirthBufsSize(sample_ct, predictor_ct, nonref_cat_ct);
  }
  return workspace_size;
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
  const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
  const AlleleCode* omitted_alleles = common->omitted_alleles;
  const uintptr_t* sex_male_collapsed = common->sex_male_collapsed;
  const ChrInfo* cip = common->cip;
  const uint32_t* subset_chr_fo_vidx_start = common->subset_chr_fo_vidx_start;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const GlmFlags glm_flags = common->glm_flags;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t is_sometimes_firth = !(glm_flags & kfGlmNoFirth);
  const uint32_t is_always_firth = (glm_flags / kfGlmFirth) & 1;
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const uintptr_t max_reported_test_ct = common->max_reported_test_ct;
  uint32_t variant_idx_offset = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);
    double* beta_se_iter = &(common->block_beta_se[2 * max_reported_test_ct * variant_bidx]);
    MnlAuxResult* block_aux_iter = &(ctx->block_aux[variant_bidx]);
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
      const GlmMnlSet* cur_set;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      if (is_y && common->sample_include_y) {
        cur_sample_include = common->sample_include_y;
        cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
        cur_set = &(ctx->mnl_set_y);
        cur_sample_ct = common->sample_ct_y;
        cur_covar_ct = common->covar_ct_y;
      } else if (is_regular_x && common->sample_include_x) {
        cur_sample_include = common->sample_include_x;
        cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
        cur_set = &(ctx->mnl_set_x);
        cur_sample_ct = common->sample_ct_x;
        cur_covar_ct = common->covar_ct_x;
      } else {
        cur_sample_include = common->sample_include;
        cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
        cur_set = &(ctx->mnl_set);
        cur_sample_ct = common->sample_ct;
        cur_covar_ct = common->covar_ct;
      }
      const uint32_t* cur_pheno_cats = cur_set->pheno_cats;
      const double* cur_covars_cmaj = cur_set->covars_cmaj;
      const uint32_t cat_ct = cur_set->cat_ct;
      const uint32_t nonref_cat_ct = cat_ct - 1;
      const uint32_t sample_ctl = BitCtToWordCt(cur_sample_ct);
      const uint32_t predictor_ct = cur_covar_ct + 2;
      const uint32_t null_predictor_ct = predictor_ct - 1;
      const uint32_t param_ct = predictor_ct * nonref_cat_ct;
      const uint32_t null_param_ct = null_predictor_ct * nonref_cat_ct;
      const uint32_t biallelic_reported_test_ct = include_intercept + 1 + (hide_covar? 0 : cur_covar_ct);
      const uint32_t reported_ct = biallelic_reported_test_ct * nonref_cat_ct + 1;
      const uint32_t sample_ctav = RoundUpPow2(cur_sample_ct, kDoublePerDVec);
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* tmp_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* cached_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uint32_t* nm_cats = S_CAST(uint32_t*, arena_alloc_raw_rd(cur_sample_ct * sizeof(int32_t), &workspace_iter));
      double* xx = S_CAST(double*, arena_alloc_raw_rd(predictor_ct * sample_ctav * sizeof(double), &workspace_iter));
      double* lli = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), &workspace_iter));
      double* cached_null_lli = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), &workspace_iter));
      double* coef = S_CAST(double*, arena_alloc_raw_rd(param_ct * sizeof(double), &workspace_iter));
      double* se = S_CAST(double*, arena_alloc_raw_rd(param_ct * sizeof(double), &workspace_iter));
      double* cached_null_coefs = S_CAST(double*, arena_alloc_raw_rd(null_param_ct * sizeof(double), &workspace_iter));
      double* cat_dosage_sums = S_CAST(double*, arena_alloc_raw_rd(cat_ct * sizeof(double), &workspace_iter));
      double* predictor_dotprod_buf = S_CAST(double*, arena_alloc_raw_rd(null_predictor_ct * null_predictor_ct * sizeof(double), &workspace_iter));
      double* dbl_2d_buf = S_CAST(double*, arena_alloc_raw_rd(predictor_ct * MAXV(predictor_ct, 7) * sizeof(double), &workspace_iter));
      double* inverse_corr_buf = S_CAST(double*, arena_alloc_raw_rd(predictor_ct * MAXV(predictor_ct, 3) * sizeof(double), &workspace_iter));
      MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, arena_alloc_raw_rd(predictor_ct * kMatrixInvertBuf1CheckedAlloc, &workspace_iter));
      MnlSolverBufs bufs;
      CarveMnlSolverBufs(cur_sample_ct, predictor_ct, nonref_cat_ct, &workspace_iter, &bufs);
      double* firth_null_coef = nullptr;
      double* firth_null_lli = nullptr;
      MnlFirthBufs fbufs;
      if (is_sometimes_firth) {
        firth_null_coef = S_CAST(double*, arena_alloc_raw_rd(param_ct * sizeof(double), &workspace_iter));
        firth_null_lli = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), &workspace_iter));
        CarveMnlFirthBufs(cur_sample_ct, predictor_ct, nonref_cat_ct, &workspace_iter, &fbufs);
      }
      assert(S_CAST(uintptr_t, workspace_iter - workspace_buf) == GetMnlWorkspaceSize(cur_sample_ct, cur_covar_ct, cat_ct, is_sometimes_firth));
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      // prev_nm: the previous variant processed here had no missing calls, so
      // the intercept and covariate rows of xx, and nm_cats, are already
      // filled for the full sample set.
      uint32_t prev_nm = 0;
      // cached_null_valid: cached_null_coefs/cached_null_lli hold the
      // covariate-only fit on the samples in cached_nm.
      uint32_t cached_null_valid = 0;

      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        GlmErr glm_err = 0;
        block_aux_iter->firth_fallback = 0;
        if (allele_idx_offsets && (allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx] != 2)) {
          glm_err = SetGlmErr0(kGlmErrcodeMultiallelicUnsupported);
          goto GlmMultinomialThread_skip_regression;
        }
        {
          uint32_t dosage_ct;
          PglErr reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
          if (unlikely(reterr)) {
            new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
            goto GlmMultinomialThread_err;
          }
          ZeroTrailingNyps(cur_sample_ct, genovec);
          GenoarrCountFreqsUnsafe(genovec, cur_sample_ct, genocounts);
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
          const uint32_t nm_sample_ctl = BitCtToWordCt(nm_sample_ct);
          const uint32_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kDoublePerDVec);
          const uint32_t nm_sample_ct_rem = nm_sample_ctav - nm_sample_ct;
          // Intercept and covariate rows, and the categories.  Unlike
          // GlmLogisticThreadD(), the intercept row is refilled whenever the
          // covariates are: its length changes with the missing-call count.
          if (missing_ct || (!prev_nm)) {
            FillDVec(nm_sample_ct, 1.0, xx);
            ZeroDArr(nm_sample_ct_rem, &(xx[nm_sample_ct]));
            if (!missing_ct) {
              memcpy(nm_cats, cur_pheno_cats, cur_sample_ct * sizeof(int32_t));
              for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
                double* xx_row = &(xx[(covar_idx + 1) * nm_sample_ctav]);
                memcpy(xx_row, &(cur_covars_cmaj[covar_idx * cur_sample_ct]), cur_sample_ct * sizeof(double));
                ZeroDArr(nm_sample_ct_rem, &(xx_row[nm_sample_ct]));
              }
            } else {
              uintptr_t sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                nm_cats[sample_idx] = cur_pheno_cats[sample_midx];
              }
              for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
                const double* cur_covar_col = &(cur_covars_cmaj[covar_idx * cur_sample_ct]);
                double* xx_row = &(xx[(covar_idx + 1) * nm_sample_ctav]);
                sample_midx_base = 0;
                sample_nm_bits = sample_nm[0];
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                  xx_row[sample_idx] = cur_covar_col[sample_midx];
                }
                ZeroDArr(nm_sample_ct_rem, &(xx_row[nm_sample_ct]));
              }
            }
            prev_nm = !missing_ct;
          }
          // Genotype row, last.
          double* genotype_vals = &(xx[null_predictor_ct * nm_sample_ctav]);
          uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
          uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
          if (!missing_ct) {
            GenoarrLookup16x8bx2(genovec, kSmallDoublePairs, nm_sample_ct, genotype_vals);
            if (dosage_ct) {
              uintptr_t sample_idx_base = 0;
              uintptr_t dosage_present_bits = dosage_present[0];
              for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
                const uintptr_t sample_idx = BitIter1(dosage_present, &sample_idx_base, &dosage_present_bits);
                const uint32_t dosage_val = dosage_main[dosage_idx];
                // 32768 -> 2, 16384 -> 1, 0 -> 0
                genotype_vals[sample_idx] = kRecipDosageMid * u31tod(dosage_val);
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
          } else {
            if (!dosage_ct) {
              GenoarrToDoublesRemoveMissing(genovec, kSmallDoubles, cur_sample_ct, genotype_vals);
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
                genotype_vals[sample_idx] = cur_val;
              }
            }
          }
          ZeroDArr(nm_sample_ct_rem, &(genotype_vals[nm_sample_ct]));
          uint64_t machr2_dosage_sums[2];
          uint64_t machr2_dosage_ssqs[2];
          machr2_dosage_sums[1] = dosage_sum;
          machr2_dosage_ssqs[1] = dosage_ssq;
          machr2_dosage_sums[0] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[0] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
          const double mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, 2);
          uint32_t allele_obs_ct = nm_sample_ct * 2;
          double a1_dosage = u63tod(dosage_sum) * kRecipDosageMid;
          if (is_nonx_haploid) {
            allele_obs_ct = nm_sample_ct;
            a1_dosage *= 0.5;
            // everything is on 0..1 scale, not 0..2
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              genotype_vals[sample_idx] *= 0.5;
            }
          } else if (is_regular_x && is_xchr_model_1) {
            // special case: multiply male values by 0.5
            CopyBitarrSubset(sex_male_collapsed, sample_nm, nm_sample_ct, tmp_nm);
            const uint32_t nm_male_ct = PopcountWords(tmp_nm, nm_sample_ctl);
            uintptr_t sample_idx_base = 0;
            uintptr_t male_nm_bits = tmp_nm[0];
            for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
              const uintptr_t sample_idx = BitIter1(tmp_nm, &sample_idx_base, &male_nm_bits);
              genotype_vals[sample_idx] *= 0.5;
            }
            allele_obs_ct -= nm_male_ct;
            a1_dosage = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              a1_dosage += genotype_vals[sample_idx];
            }
          }
          block_aux_iter->sample_obs_ct = nm_sample_ct;
          block_aux_iter->allele_obs_ct = allele_obs_ct;
          block_aux_iter->a1_dosage = a1_dosage;
          block_aux_iter->mach_r2 = mach_r2;

          // The model has (cat_ct - 1) * predictor_ct parameters; reported
          // with the existing SAMPLE_CT<=PREDICTOR_CT code.
          if (nm_sample_ct <= param_ct) {
            glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
            goto GlmMultinomialThread_skip_regression;
          }
          {
            const double first_val = genotype_vals[0];
            uint32_t sample_idx = 1;
            for (; sample_idx != nm_sample_ct; ++sample_idx) {
              if (genotype_vals[sample_idx] != first_val) {
                break;
              }
            }
            if (sample_idx == nm_sample_ct) {
              glm_err = SetGlmErr0(kGlmErrcodeConstOmittedAllele);
              goto GlmMultinomialThread_skip_regression;
            }
          }
          // A category whose samples all have zero dosage has an infinite
          // maximum-likelihood genotype coefficient.  (This is the multinomial
          // version of the check GlmLogisticThreadD() makes before fitting.)
          // In Firth-fallback mode, the penalized fit is used instead.
          uint32_t use_firth = is_always_firth;
          if (!is_always_firth) {
            ZeroDArr(cat_ct, cat_dosage_sums);
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              cat_dosage_sums[nm_cats[sample_idx]] += genotype_vals[sample_idx];
            }
            for (uint32_t cat_idx = 0; cat_idx != cat_ct; ++cat_idx) {
              if (cat_dosage_sums[cat_idx] == 0.0) {
                if (!is_sometimes_firth) {
                  glm_err = SetGlmErr1(kGlmErrcodeSeparation, 1 - omitted_allele_idx);
                  goto GlmMultinomialThread_skip_regression;
                }
                use_firth = 1;
                break;
              }
            }
          }
          // Same collinearity checks as the other --glm regressions, on
          // everything but the intercept.
          MultiplySelfTransposeStrided(&(xx[nm_sample_ctav]), null_predictor_ct, nm_sample_ct, nm_sample_ctav, predictor_dotprod_buf);
          for (uint32_t pred_idx = 1; pred_idx != predictor_ct; ++pred_idx) {
            const double* xx_row = &(xx[pred_idx * nm_sample_ctav]);
            double row_sum = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              row_sum += xx_row[sample_idx];
            }
            dbl_2d_buf[pred_idx - 1] = row_sum;
          }
          glm_err = CheckMaxCorrAndVif(predictor_dotprod_buf, 0, null_predictor_ct, nm_sample_ct, max_corr, vif_thresh, dbl_2d_buf, nullptr, inverse_corr_buf, inv_1d_buf);
          if (glm_err) {
            goto GlmMultinomialThread_skip_regression;
          }

          double lrt_chisq = 0.0;
          if (!use_firth) {
            // Covariate-only fit on this variant's samples.
            const double* null_coefs;
            const double* null_lli;
            if (!missing_ct) {
              null_coefs = cur_set->null_coefs;
              null_lli = cur_set->null_lli;
            } else {
              if ((!cached_null_valid) || (!wordsequal(sample_nm, cached_nm, sample_ctl))) {
                // Always start from the full-sample fit, so the result
                // depends only on the sample set, not on which variants this
                // thread saw before.
                memcpy(cached_null_coefs, cur_set->null_coefs, null_param_ct * sizeof(double));
                const GlmErrcode errcode = MultinomialRegressionD(xx, nm_cats, nm_sample_ct, null_predictor_ct, nonref_cat_ct, cached_null_coefs, cached_null_lli, nullptr, &bufs);
                if (errcode) {
                  cached_null_valid = 0;
                  if (!is_sometimes_firth) {
                    glm_err = SetGlmErr0(kGlmErrcodeMnlConvergeFail);
                    goto GlmMultinomialThread_skip_regression;
                  }
                  // e.g. a category has lost all its samples to missing calls
                  use_firth = 1;
                } else {
                  memcpy(cached_nm, sample_nm, sample_ctl * sizeof(intptr_t));
                  cached_null_valid = 1;
                }
              }
              null_coefs = cached_null_coefs;
              null_lli = cached_null_lli;
            }

            if (!use_firth) {
              // Full fit, starting from the covariate-only estimates with zero
              // genotype coefficients.
              for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
                memcpy(&(coef[cat_idx * predictor_ct]), &(null_coefs[cat_idx * null_predictor_ct]), null_predictor_ct * sizeof(double));
                coef[cat_idx * predictor_ct + null_predictor_ct] = 0.0;
              }
              const GlmErrcode errcode = MultinomialRegressionD(xx, nm_cats, nm_sample_ct, predictor_ct, nonref_cat_ct, coef, lli, se, &bufs);
              if (errcode) {
                if (!is_sometimes_firth) {
                  glm_err = SetGlmErr0(errcode);
                  goto GlmMultinomialThread_skip_regression;
                }
                use_firth = 1;
              } else {
                // Likelihood-ratio statistic.  Summing the per-sample
                // differences keeps its accuracy relative to the statistic
                // itself rather than to the (much larger) log-likelihoods.
                double lrt_half = 0.0;
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  lrt_half += lli[sample_idx] - null_lli[sample_idx];
                }
                lrt_chisq = 2 * lrt_half;
              }
            }
          }
          if (use_firth) {
            if (!is_always_firth) {
              block_aux_iter->firth_fallback = 1;
            }
            // Penalized likelihood-ratio test (logistf's PLR): the restricted
            // fit maximizes the same penalized likelihood, whose penalty is
            // that of the full model, with the genotype coefficients held at
            // zero.  It starts from the unpenalized covariate-only fit on the
            // whole sample set, and the full fit starts from its result.
            for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
              memcpy(&(firth_null_coef[cat_idx * predictor_ct]), &(cur_set->null_coefs[cat_idx * null_predictor_ct]), null_predictor_ct * sizeof(double));
              firth_null_coef[cat_idx * predictor_ct + null_predictor_ct] = 0.0;
            }
            double null_logdet;
            if (MultinomialFirthRegressionD(xx, nm_cats, nm_sample_ct, predictor_ct, nonref_cat_ct, null_predictor_ct, firth_null_coef, firth_null_lli, &null_logdet, nullptr, &bufs, &fbufs)) {
              glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
              goto GlmMultinomialThread_skip_regression;
            }
            memcpy(coef, firth_null_coef, param_ct * sizeof(double));
            double full_logdet;
            const GlmErrcode errcode = MultinomialFirthRegressionD(xx, nm_cats, nm_sample_ct, predictor_ct, nonref_cat_ct, UINT32_MAX, coef, lli, &full_logdet, se, &bufs, &fbufs);
            if (errcode) {
              glm_err = SetGlmErr0(errcode);
              goto GlmMultinomialThread_skip_regression;
            }
            double lrt_half = 0.5 * (full_logdet - null_logdet);
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              lrt_half += lli[sample_idx] - firth_null_lli[sample_idx];
            }
            lrt_chisq = 2 * lrt_half;
          }
          if (lrt_chisq < 0.0) {
            lrt_chisq = 0.0;
          }
          // Reported order: [intercept,] genotype, [covariates]; category
          // varies fastest.
          double* beta_se_iter2 = beta_se_iter;
          uint32_t pred_idx = include_intercept? 0 : null_predictor_ct;
          for (uint32_t test_idx = 0; test_idx != biallelic_reported_test_ct; ++test_idx) {
            for (uint32_t cat_idx = 0; cat_idx != nonref_cat_ct; ++cat_idx) {
              *beta_se_iter2++ = coef[cat_idx * predictor_ct + pred_idx];
              *beta_se_iter2++ = se[cat_idx * predictor_ct + pred_idx];
            }
            if (pred_idx == null_predictor_ct) {
              pred_idx = 1;
            } else if (!pred_idx) {
              pred_idx = null_predictor_ct;
            } else {
              ++pred_idx;
            }
          }
          *beta_se_iter2++ = lrt_chisq;
          *beta_se_iter2++ = 0.0;
        }
        while (0) {
        GlmMultinomialThread_skip_regression:
          for (uint32_t uii = 0; uii != reported_ct; ++uii) {
            memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
            beta_se_iter[uii * 2 + 1] = -9.0;
          }
        }
        beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
        ++block_aux_iter;
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

PglErr GlmMultinomial(const char* cur_pheno_name, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const char* outname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, PgenFileInfo* pgfip, GlmMultinomialCtx* ctx, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, uintptr_t* valid_allele_ct_ptr) {
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
    const uint32_t covar_ct_x = common->covar_ct_x;
    const uint32_t covar_ct_y = common->covar_ct_y;
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
    const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
    const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
    const uint32_t is_sometimes_firth = !(glm_flags & kfGlmNoFirth);
    const uint32_t is_always_firth = (glm_flags / kfGlmFirth) & 1;
    const GlmColFlags glm_cols = glm_info_ptr->cols;
    // checked up front by GlmMain()
    assert(glm_cols & kfGlmColTest);

    // Per sample set: nonreference category count, reported tests per
    // category, and the omnibus test's name.
    uint32_t nonref_cat_cts[3];
    uint32_t biallelic_reported_test_cts[3];
    char lrt_names[3][16];
    const GlmMnlSet* sets[3] = {&(ctx->mnl_set), &(ctx->mnl_set_x), &(ctx->mnl_set_y)};
    const uint32_t set_sample_cts[3] = {sample_ct, sample_ct_x, sample_ct_y};
    const uint32_t set_covar_cts[3] = {covar_ct, covar_ct_x, covar_ct_y};
    uintptr_t max_reported_test_ct = 0;
    uintptr_t workspace_alloc = 0;
    for (uint32_t set_idx = 0; set_idx != 3; ++set_idx) {
      nonref_cat_cts[set_idx] = 0;
      biallelic_reported_test_cts[set_idx] = 0;
      lrt_names[set_idx][0] = '\0';
      if (!set_sample_cts[set_idx]) {
        continue;
      }
      const uint32_t cur_nonref_cat_ct = sets[set_idx]->cat_ct - 1;
      nonref_cat_cts[set_idx] = cur_nonref_cat_ct;
      biallelic_reported_test_cts[set_idx] = include_intercept + 1 + (hide_covar? 0 : set_covar_cts[set_idx]);
      const uintptr_t cur_reported_test_ct = biallelic_reported_test_cts[set_idx] * cur_nonref_cat_ct + 1;
      if (cur_reported_test_ct > max_reported_test_ct) {
        max_reported_test_ct = cur_reported_test_ct;
      }
      char* write_iter = strcpya_k(lrt_names[set_idx], "LRT_");
      write_iter = u32toa(cur_nonref_cat_ct, write_iter);
      strcpy_k(write_iter, "DF");
      const uintptr_t cur_workspace_alloc = GetMnlWorkspaceSize(set_sample_cts[set_idx], set_covar_cts[set_idx], sets[set_idx]->cat_ct, is_sometimes_firth);
      if (cur_workspace_alloc > workspace_alloc) {
        workspace_alloc = cur_workspace_alloc;
      }
    }
    common->max_reported_test_ct = max_reported_test_ct;

    uint32_t x_code = UINT32_MAXM1;
    uint32_t y_code = UINT32_MAXM1;
    if (sample_ct_x) {
      x_code = cip->xymt_codes[kChrOffsetX];
    }
    if (sample_ct_y) {
      y_code = cip->xymt_codes[kChrOffsetY];
    }
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
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
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    // +1 is for top-level common->workspace_bufs
    const uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;
    const uintptr_t per_variant_xalloc_byte_ct = sizeof(MnlAuxResult) + 2 * max_reported_test_ct * sizeof(double);
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
    MnlAuxResult* block_aux_bufs[2];
    double* block_beta_se_bufs[2];
    for (uint32_t uii = 0; uii != 2; ++uii) {
      if (unlikely(BIGSTACK_ALLOC_X(MnlAuxResult, read_block_size, &(block_aux_bufs[uii])) ||
                   bigstack_alloc_d(read_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii])))) {
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
    const uint32_t a1_freq_col = glm_cols & kfGlmColA1freq;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t firth_yn_col = (glm_cols & kfGlmColFirthYn) && is_sometimes_firth && (!is_always_firth);
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t orbeta_col = glm_cols & (kfGlmColBeta | kfGlmColOrbeta);
    const uint32_t report_beta_instead_of_odds_ratio = glm_cols & kfGlmColBeta;
    const uint32_t se_col = glm_cols & kfGlmColSe;
    const uint32_t ci_col = (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t z_col = glm_cols & kfGlmColTz;
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
    if (firth_yn_col) {
      cswritep = strcpya_k(cswritep, "\tFIRTH?");
    }
    cswritep = strcpya_k(cswritep, "\tTEST\tCATEGORY");
    if (nobs_col) {
      cswritep = strcpya_k(cswritep, "\tOBS_CT");
    }
    if (orbeta_col) {
      if (report_beta_instead_of_odds_ratio) {
        cswritep = strcpya_k(cswritep, "\tBETA");
      } else {
        cswritep = strcpya_k(cswritep, "\tOR");
      }
    }
    if (se_col) {
      if (report_beta_instead_of_odds_ratio) {
        cswritep = strcpya_k(cswritep, "\tSE");
      } else {
        cswritep = strcpya_k(cswritep, "\tLOG(OR)_SE");
      }
    }
    double ci_zt = 0.0;
    if (ci_col) {
      cswritep = strcpya_k(cswritep, "\tL");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      cswritep = strcpya_k(cswritep, "\tU");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      ci_zt = QuantileToZscore((ci_size + 1.0) * 0.5);
    }
    if (z_col) {
      // Wald Z statistic on the per-category rows, likelihood-ratio
      // chi-square on the omnibus row.
      cswritep = strcpya_k(cswritep, "\tZ_OR_CHISQ_STAT");
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

    // Same block pipeline as GlmLogistic():
    // 1. Set n=0, load/skip block 0
    //
    // 2. Spawn threads processing block n
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Load/skip block n unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8. Write results for last block
    uintptr_t write_variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t suppress_mach_r2 = 0;
    uint32_t cur_set_idx = 0;
    const char* const* cur_test_names = nullptr;
    uint32_t prev_block_variant_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    uint32_t allele_ct = 2;
    uint32_t omitted_allele_idx = 0;
    uintptr_t valid_allele_ct = 0;
    logprintfww5("--glm %s regression on phenotype '%s': ", is_always_firth? "Firth multinomial logistic" : (is_sometimes_firth? "multinomial logistic-Firth hybrid" : "multinomial logistic"), cur_pheno_name);
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
        common->block_beta_se = block_beta_se_bufs[parity];
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
        const double* beta_se_iter = block_beta_se_bufs[parity];
        const MnlAuxResult* auxp = block_aux_bufs[parity];
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx, ++auxp, beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct])) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            if (chr_idx == x_code) {
              cur_set_idx = 1;
              cur_test_names = test_names_x;
            } else if (chr_idx == y_code) {
              cur_set_idx = 2;
              cur_test_names = test_names_y;
            } else {
              cur_set_idx = 0;
              cur_test_names = test_names;
            }
            suppress_mach_r2 = (chr_idx == x_code) || (chr_idx == mt_code);
            if (chr_col) {
              char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
              *chr_name_end = '\t';
              chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
            }
          }
          const uint32_t nonref_cat_ct = nonref_cat_cts[cur_set_idx];
          const uint32_t biallelic_reported_test_ct = biallelic_reported_test_cts[cur_set_idx];
          const char* const* cur_cat_names = sets[cur_set_idx]->cat_names;
          const uint32_t lrt_slot_idx = biallelic_reported_test_ct * nonref_cat_ct;
          uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offsets[write_variant_uidx];
          }
          const uint32_t is_multiallelic = (allele_ct != 2);
          if (omitted_alleles) {
            omitted_allele_idx = omitted_alleles[write_variant_uidx];
          }
          const uint32_t a1_allele_idx = (omitted_allele_idx == 0);
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          const uint32_t variant_is_valid = (beta_se_iter[2 * lrt_slot_idx + 1] != -9.0);
          double lrt_ln_pval = kLnPvalError;
          if (variant_is_valid) {
            lrt_ln_pval = ChisqToLnP(beta_se_iter[2 * lrt_slot_idx], nonref_cat_ct);
          }
          if (ln_pfilter <= 0.0) {
            if ((!variant_is_valid) || (lrt_ln_pval > ln_pfilter)) {
              if (variant_is_valid && orig_ln_pvals) {
                orig_ln_pvals[valid_allele_ct] = lrt_ln_pval;
              }
              goto GlmMultinomial_variant_iterate;
            }
          }
          for (uint32_t slot_idx = 0; slot_idx <= lrt_slot_idx; ++slot_idx) {
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
              for (uint32_t tmp_allele_idx = 1; tmp_allele_idx != allele_ct; ++tmp_allele_idx) {
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto GlmMultinomial_ret_WRITE_FAIL;
                }
                cswritep = strcpyax(cswritep, cur_alleles[tmp_allele_idx], ',');
              }
              --cswritep;
            }
            *cswritep++ = '\t';
            if (provref_col) {
              *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, write_variant_uidx)))? 'Y' : 'N';
              *cswritep++ = '\t';
            }
            if (!is_multiallelic) {
              cswritep = strcpya(cswritep, cur_alleles[a1_allele_idx]);
            } else {
              // No genotype coding has been chosen yet for multiallelic
              // variants; list every non-omitted allele.
              for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                if (allele_idx == omitted_allele_idx) {
                  continue;
                }
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto GlmMultinomial_ret_WRITE_FAIL;
                }
                cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
              }
              --cswritep;
            }
            if (omitted_col) {
              *cswritep++ = '\t';
              cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
            }
            if (ax_col) {
              *cswritep++ = '\t';
              if (!is_multiallelic) {
                cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
              } else {
                cswritep = strcpya_k(cswritep, "NA");
              }
            }
            // Per-variant counts are only computed for biallelic variants.
            if (a1_ct_col) {
              *cswritep++ = '\t';
              if (!is_multiallelic) {
                cswritep = dtoa_g(auxp->a1_dosage, cswritep);
              } else {
                cswritep = strcpya_k(cswritep, "NA");
              }
            }
            if (tot_allele_col) {
              *cswritep++ = '\t';
              if (!is_multiallelic) {
                cswritep = u32toa(auxp->allele_obs_ct, cswritep);
              } else {
                cswritep = strcpya_k(cswritep, "NA");
              }
            }
            if (a1_freq_col) {
              *cswritep++ = '\t';
              if (!is_multiallelic) {
                cswritep = dtoa_g(auxp->a1_dosage / S_CAST(double, auxp->allele_obs_ct), cswritep);
              } else {
                cswritep = strcpya_k(cswritep, "NA");
              }
            }
            if (mach_r2_col) {
              *cswritep++ = '\t';
              if ((!is_multiallelic) && (!suppress_mach_r2)) {
                cswritep = dtoa_g(auxp->mach_r2, cswritep);
              } else {
                cswritep = strcpya_k(cswritep, "NA");
              }
            }
            if (firth_yn_col) {
              *cswritep++ = '\t';
              // 'Y' - 'N' = 11
              *cswritep++ = 'N' + 11 * auxp->firth_fallback;
            }
            *cswritep++ = '\t';
            double ln_pval = kLnPvalError;
            uint32_t test_is_valid;
            if (slot_idx != lrt_slot_idx) {
              const uint32_t test_idx = slot_idx / nonref_cat_ct;
              const uint32_t cat_idx = slot_idx - test_idx * nonref_cat_ct;
              cswritep = strcpyax(cswritep, cur_test_names[test_idx], '\t');
              cswritep = strcpya(cswritep, cur_cat_names[cat_idx]);
              if (nobs_col) {
                *cswritep++ = '\t';
                if (!is_multiallelic) {
                  cswritep = u32toa(auxp->sample_obs_ct, cswritep);
                } else {
                  cswritep = strcpya_k(cswritep, "NA");
                }
              }
              const double beta = beta_se_iter[2 * slot_idx];
              const double se = beta_se_iter[2 * slot_idx + 1];
              test_is_valid = (se != -9.0);
              double zscore = 0.0;
              if (test_is_valid) {
                zscore = beta / se;
                ln_pval = ZscoreToLnP(zscore);
              }
              if (orbeta_col) {
                *cswritep++ = '\t';
                if (test_is_valid) {
                  if (report_beta_instead_of_odds_ratio) {
                    cswritep = dtoa_g(beta, cswritep);
                  } else {
                    cswritep = lntoa_g(beta, cswritep);
                  }
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
                  const double ci_radius = ci_zt * se;
                  if (report_beta_instead_of_odds_ratio) {
                    cswritep = dtoa_g(beta - ci_radius, cswritep);
                    *cswritep++ = '\t';
                    cswritep = dtoa_g(beta + ci_radius, cswritep);
                  } else {
                    cswritep = lntoa_g(beta - ci_radius, cswritep);
                    *cswritep++ = '\t';
                    cswritep = lntoa_g(beta + ci_radius, cswritep);
                  }
                } else {
                  cswritep = strcpya_k(cswritep, "NA\tNA");
                }
              }
              if (z_col) {
                *cswritep++ = '\t';
                if (test_is_valid) {
                  cswritep = dtoa_g(zscore, cswritep);
                } else {
                  cswritep = strcpya_k(cswritep, "NA");
                }
              }
            } else {
              cswritep = strcpya(cswritep, lrt_names[cur_set_idx]);
              cswritep = strcpya_k(cswritep, "\tNA");
              if (nobs_col) {
                *cswritep++ = '\t';
                if (!is_multiallelic) {
                  cswritep = u32toa(auxp->sample_obs_ct, cswritep);
                } else {
                  cswritep = strcpya_k(cswritep, "NA");
                }
              }
              test_is_valid = variant_is_valid;
              ln_pval = lrt_ln_pval;
              if (orbeta_col) {
                cswritep = strcpya_k(cswritep, "\tNA");
              }
              if (se_col) {
                cswritep = strcpya_k(cswritep, "\tNA");
              }
              if (ci_col) {
                cswritep = strcpya_k(cswritep, "\tNA\tNA");
              }
              if (z_col) {
                *cswritep++ = '\t';
                if (test_is_valid) {
                  cswritep = dtoa_g(beta_se_iter[2 * lrt_slot_idx], cswritep);
                } else {
                  cswritep = strcpya_k(cswritep, "NA");
                }
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
                memcpy(&glm_errcode, &(beta_se_iter[2 * slot_idx]), 8);
                cswritep = AppendGlmErrstr(glm_errcode, cswritep);
              }
            }
            AppendBinaryEoln(&cswritep);
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto GlmMultinomial_ret_WRITE_FAIL;
            }
          }
          if (variant_is_valid && orig_ln_pvals) {
            orig_ln_pvals[valid_allele_ct] = lrt_ln_pval;
          }
        GlmMultinomial_variant_iterate:
          if (variant_is_valid) {
            ++valid_allele_ct;
            if (valid_alleles) {
              SetBit(allele_idx_offset_base + a1_allele_idx, valid_alleles);
            }
          } else if (valid_variants) {
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
