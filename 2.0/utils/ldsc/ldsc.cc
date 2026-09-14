// Copyright (C) 2026 Christopher Chang, Benjamin Demaille.
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
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.

// LD Score regression: SNP-heritability (Bulik-Sullivan et al., Nat Genet
// 47:291-295, 2015) and genetic correlation (Bulik-Sullivan et al., Nat Genet
// 47:1236-1241, 2015), from GWAS summary statistics plus LD Scores.
//
// plink2's --ld-score computes the LD Scores; this does the regression on top
// of them.  It follows the reference implementation (bulik/ldsc) closely
// enough to reproduce its numbers rather than merely its method: the same
// regression weights, the same two-iteration IRWLS, the same two-step
// estimator, the same block jackknife, and the same merge order (the LD Score
// file's, which is what the jackknife blocks are cut on).
//
// The reference implementation requires Python 2.7, which reached end of life
// in January 2020, and loads everything through pandas.
//
// Only unpartitioned LD Scores are supported so far; the partitioned
// (stratified) regression, which needs one LD Score column per annotation,
// is not implemented, and neither is --h2-cts.

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "../../include/plink2_base.h"
#include "../../include/plink2_string.h"
#include "../../include/plink2_text.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static const char kLdscVersion[] = "ldsc (plink-ng) v0.1";

// Every regression here has at most two parameters: one LD Score coefficient
// and one intercept.  (Partitioned LD Score regression is what needs more.)
static const uint32_t kLdscMaxP = 2;

static const uint32_t kLdscChrCt = 22;

static const double kLdscSqrt2Pi = 2.5066282746310002;

// ***** small-matrix and distribution helpers *****

// Solves mat * out = vec by Gauss-Jordan with partial pivoting.  dim is 1 or
// 2 here, so a dedicated small-matrix routine beats calling out to LAPACK.
BoolErr LdscSolve(uint32_t dim, const double* mat, const double* vec, double* out) {
  double aug[kLdscMaxP * (kLdscMaxP + 1)];
  const uint32_t width = dim + 1;
  for (uint32_t i = 0; i != dim; ++i) {
    for (uint32_t j = 0; j != dim; ++j) {
      aug[i * width + j] = mat[i * dim + j];
    }
    aug[i * width + dim] = vec[i];
  }
  for (uint32_t col = 0; col != dim; ++col) {
    uint32_t pivot = col;
    double best = fabs(aug[col * width + col]);
    for (uint32_t row = col + 1; row != dim; ++row) {
      const double cur = fabs(aug[row * width + col]);
      if (cur > best) {
        best = cur;
        pivot = row;
      }
    }
    if (best == 0.0) {
      return 1;
    }
    if (pivot != col) {
      for (uint32_t j = 0; j != width; ++j) {
        const double tmp = aug[col * width + j];
        aug[col * width + j] = aug[pivot * width + j];
        aug[pivot * width + j] = tmp;
      }
    }
    const double inv_pivot = 1.0 / aug[col * width + col];
    for (uint32_t j = 0; j != width; ++j) {
      aug[col * width + j] *= inv_pivot;
    }
    for (uint32_t row = 0; row != dim; ++row) {
      if (row == col) {
        continue;
      }
      const double factor = aug[row * width + col];
      if (factor == 0.0) {
        continue;
      }
      for (uint32_t j = 0; j != width; ++j) {
        aug[row * width + j] -= factor * aug[col * width + j];
      }
    }
  }
  for (uint32_t i = 0; i != dim; ++i) {
    out[i] = aug[i * width + dim];
  }
  return 0;
}

double LdscNormalPdf(double x) {
  return exp(-0.5 * x * x) / kLdscSqrt2Pi;
}

// P(|Z| > |z|), i.e. the two-sided normal tail, which is also the chi-square
// (1df) survival function at z^2.
double LdscChisq1Sf(double z) {
  return erfc(fabs(z) * M_SQRT1_2);
}

// Inverse survival function of the standard normal: returns the x with
// P(X > x) = p.  Acklam's rational approximation plus one Halley step, which
// is more than accurate enough for a prevalence conversion.
double LdscNormalIsf(double p) {
  if (p <= 0.0) {
    return INFINITY;
  }
  if (p >= 1.0) {
    return -INFINITY;
  }
  // Solve for the lower tail, then negate.
  const double q = p;
  static const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[4] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
  const double p_low = 0.02425;
  double x;
  if (q < p_low) {
    const double u = sqrt(-2 * log(q));
    x = (((((c[0] * u + c[1]) * u + c[2]) * u + c[3]) * u + c[4]) * u + c[5]) / ((((d[0] * u + d[1]) * u + d[2]) * u + d[3]) * u + 1);
  } else if (q <= 1 - p_low) {
    const double u = q - 0.5;
    const double r = u * u;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * u / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  } else {
    const double u = sqrt(-2 * log(1 - q));
    x = -(((((c[0] * u + c[1]) * u + c[2]) * u + c[3]) * u + c[4]) * u + c[5]) / ((((d[0] * u + d[1]) * u + d[2]) * u + d[3]) * u + 1);
  }
  // Halley refinement.
  const double e = 0.5 * erfc(-x * M_SQRT1_2) - q;
  const double u = e * kLdscSqrt2Pi * exp(x * x / 2);
  x = x - u / (1 + x * u / 2);
  return -x;
}

double LdscMedian(const double* vals, uintptr_t ct) {
  std::vector<double> buf(vals, &(vals[ct]));
  const uintptr_t mid = ct / 2;
  std::nth_element(buf.begin(), buf.begin() + mid, buf.end());
  const double hi = buf[mid];
  if (ct % 2) {
    return hi;
  }
  const double lo = *std::max_element(buf.begin(), buf.begin() + mid);
  return 0.5 * (lo + hi);
}

double LdscClamp(double val, double lo, double hi) {
  if (val < lo) {
    return lo;
  }
  if (val > hi) {
    return hi;
  }
  return val;
}

// Observed-scale to liability-scale conversion factor for a case/control
// study with sample prevalence samp_prev and population prevalence pop_prev.
double LdscLiabilityFactor(double samp_prev, double pop_prev) {
  const double thresh = LdscNormalIsf(pop_prev);
  const double dens = LdscNormalPdf(thresh);
  return pop_prev * pop_prev * (1 - pop_prev) * (1 - pop_prev) / (samp_prev * (1 - samp_prev) * dens * dens);
}

// ***** block jackknife *****

typedef struct LdscJknifeStruct {
  uint32_t n_blocks;
  uint32_t p;
  double est[kLdscMaxP];
  double jknife_est[kLdscMaxP];
  double jknife_se[kLdscMaxP];
  double jknife_cov[kLdscMaxP * kLdscMaxP];
  std::vector<double> delete_values;  // n_blocks x p, row-major
} LdscJknife;

// Turns delete values and the whole-data estimate into the jackknife estimate
// and its covariance, via the pseudovalues.
void LdscFinishJknife(LdscJknife* jkp) {
  const uint32_t n_blocks = jkp->n_blocks;
  const uint32_t p = jkp->p;
  std::vector<double> pseudo(S_CAST(uintptr_t, n_blocks) * p);
  for (uint32_t b = 0; b != n_blocks; ++b) {
    for (uint32_t j = 0; j != p; ++j) {
      pseudo[b * p + j] = n_blocks * jkp->est[j] - (n_blocks - 1) * jkp->delete_values[b * p + j];
    }
  }
  for (uint32_t j = 0; j != p; ++j) {
    double acc = 0.0;
    for (uint32_t b = 0; b != n_blocks; ++b) {
      acc += pseudo[b * p + j];
    }
    jkp->jknife_est[j] = acc / u31tod(n_blocks);
  }
  // np.cov(..., ddof=1) / n_blocks
  for (uint32_t j = 0; j != p; ++j) {
    for (uint32_t k = 0; k != p; ++k) {
      double acc = 0.0;
      for (uint32_t b = 0; b != n_blocks; ++b) {
        acc += (pseudo[b * p + j] - jkp->jknife_est[j]) * (pseudo[b * p + k] - jkp->jknife_est[k]);
      }
      jkp->jknife_cov[j * p + k] = acc / (u31tod(n_blocks - 1) * u31tod(n_blocks));
    }
  }
  for (uint32_t j = 0; j != p; ++j) {
    jkp->jknife_se[j] = sqrt(jkp->jknife_cov[j * p + j]);
  }
}

// Evenly-spaced block boundaries, as np.floor(np.linspace(0, n, n_blocks+1)).
void LdscGetSeparators(uint32_t n, uint32_t n_blocks, std::vector<uint32_t>* sep) {
  sep->resize(n_blocks + 1);
  for (uint32_t i = 0; i <= n_blocks; ++i) {
    (*sep)[i] = S_CAST(uint32_t, floor(u31tod(i) * u31tod(n) / u31tod(n_blocks)));
  }
  (*sep)[0] = 0;
  (*sep)[n_blocks] = n;
}

// Fast block jackknife for linear regression: accumulate X'X and X'y per
// block, then solve once for the whole-data estimate and once per block with
// that block's contribution removed.
BoolErr LdscLstsqJknifeFast(const double* x, const double* y, uint32_t p, const std::vector<uint32_t>& sep, LdscJknife* jkp) {
  const uint32_t n_blocks = sep.size() - 1;
  const uint32_t pp = p * p;
  std::vector<double> xtx_blocks(S_CAST(uintptr_t, n_blocks) * pp, 0.0);
  std::vector<double> xty_blocks(S_CAST(uintptr_t, n_blocks) * p, 0.0);
  for (uint32_t b = 0; b != n_blocks; ++b) {
    double* cur_xtx = &(xtx_blocks[S_CAST(uintptr_t, b) * pp]);
    double* cur_xty = &(xty_blocks[S_CAST(uintptr_t, b) * p]);
    for (uint32_t i = sep[b]; i != sep[b + 1]; ++i) {
      const double* cur_x = &(x[S_CAST(uintptr_t, i) * p]);
      for (uint32_t j = 0; j != p; ++j) {
        cur_xty[j] += cur_x[j] * y[i];
        for (uint32_t k = 0; k != p; ++k) {
          cur_xtx[j * p + k] += cur_x[j] * cur_x[k];
        }
      }
    }
  }
  double xtx_tot[kLdscMaxP * kLdscMaxP];
  double xty_tot[kLdscMaxP];
  for (uint32_t j = 0; j != pp; ++j) {
    double acc = 0.0;
    for (uint32_t b = 0; b != n_blocks; ++b) {
      acc += xtx_blocks[S_CAST(uintptr_t, b) * pp + j];
    }
    xtx_tot[j] = acc;
  }
  for (uint32_t j = 0; j != p; ++j) {
    double acc = 0.0;
    for (uint32_t b = 0; b != n_blocks; ++b) {
      acc += xty_blocks[S_CAST(uintptr_t, b) * p + j];
    }
    xty_tot[j] = acc;
  }
  jkp->n_blocks = n_blocks;
  jkp->p = p;
  if (LdscSolve(p, xtx_tot, xty_tot, jkp->est)) {
    return 1;
  }
  jkp->delete_values.assign(S_CAST(uintptr_t, n_blocks) * p, 0.0);
  double xtx_del[kLdscMaxP * kLdscMaxP];
  double xty_del[kLdscMaxP];
  for (uint32_t b = 0; b != n_blocks; ++b) {
    for (uint32_t j = 0; j != pp; ++j) {
      xtx_del[j] = xtx_tot[j] - xtx_blocks[S_CAST(uintptr_t, b) * pp + j];
    }
    for (uint32_t j = 0; j != p; ++j) {
      xty_del[j] = xty_tot[j] - xty_blocks[S_CAST(uintptr_t, b) * p + j];
    }
    if (LdscSolve(p, xtx_del, xty_del, &(jkp->delete_values[S_CAST(uintptr_t, b) * p]))) {
      return 1;
    }
  }
  LdscFinishJknife(jkp);
  return 0;
}

// Block jackknife for a ratio of two jackknifed quantities.
void LdscRatioJknife(double est, const double* numer_delete, const double* denom_delete, uint32_t n_blocks, double* jknife_est_ptr, double* jknife_se_ptr) {
  LdscJknife jk;
  jk.n_blocks = n_blocks;
  jk.p = 1;
  jk.est[0] = est;
  jk.delete_values.resize(n_blocks);
  for (uint32_t b = 0; b != n_blocks; ++b) {
    jk.delete_values[b] = numer_delete[b] / denom_delete[b];
  }
  LdscFinishJknife(&jk);
  *jknife_est_ptr = jk.jknife_est[0];
  *jknife_se_ptr = jk.jknife_se[0];
}

// ***** regression weights and IRWLS *****

typedef enum { kLdscHsq, kLdscGencov } LdscRegKind;

// Everything the weight update needs.  The per-SNP arrays are full-length;
// row_idxs maps a regression row to its entry in them, so the two-step
// estimator's first step can run on a subset without copying the inputs.
typedef struct LdscRegCtxStruct {
  LdscRegKind kind;
  const double* ld;
  const double* w_ld;
  const double* n_vec;
  const double* n1;
  const double* n2;
  const uint32_t* row_idxs;
  double m_tot;
  double nbar;
  double hsq1;
  double hsq2;
  double intercept_hsq1;
  double intercept_hsq2;
  uint32_t constrain_intercept;
  double fixed_intercept;
} LdscRegCtx;

void LdscUpdateWeights(const LdscRegCtx* ctx, const double* coef, uint32_t p, uint32_t n, double* w_out) {
  const double m_tot = ctx->m_tot;
  const double param = m_tot * coef[0] / ctx->nbar;
  double intercept;
  if (ctx->constrain_intercept) {
    intercept = ctx->fixed_intercept;
  } else {
    assert(p == 2);
    intercept = coef[1];
  }
  if (ctx->kind == kLdscHsq) {
    const double hsq = LdscClamp(param, 0.0, 1.0);
    for (uint32_t i = 0; i != n; ++i) {
      const uint32_t uidx = ctx->row_idxs[i];
      const double ld = MAXV(ctx->ld[uidx], 1.0);
      const double w_ld = MAXV(ctx->w_ld[uidx], 1.0);
      const double c = hsq * ctx->n_vec[uidx] / m_tot;
      const double denom = intercept + c * ld;
      w_out[i] = 1.0 / (2 * denom * denom * w_ld);
    }
    return;
  }
  const double rho_g = LdscClamp(param, -1.0, 1.0);
  const double h1 = LdscClamp(ctx->hsq1, 0.0, 1.0);
  const double h2 = LdscClamp(ctx->hsq2, 0.0, 1.0);
  for (uint32_t i = 0; i != n; ++i) {
    const uint32_t uidx = ctx->row_idxs[i];
    const double ld = MAXV(ctx->ld[uidx], 1.0);
    const double w_ld = MAXV(ctx->w_ld[uidx], 1.0);
    const double cur_n1 = ctx->n1[uidx];
    const double cur_n2 = ctx->n2[uidx];
    const double a = cur_n1 * h1 * ld / m_tot + ctx->intercept_hsq1;
    const double b = cur_n2 * h2 * ld / m_tot + ctx->intercept_hsq2;
    const double c = sqrt(cur_n1 * cur_n2) * rho_g * ld / m_tot + intercept;
    w_out[i] = 1.0 / ((a * b + c * c) * w_ld);
  }
}

// Weighted least squares: minimizes sum_i w_i^2 (y_i - x_i b)^2.  (The
// reference implementation normalizes the weights to sum 1 first; that
// cancels out of the solution.)
BoolErr LdscWls(const double* x, const double* y, const double* w, uint32_t n, uint32_t p, double* coef) {
  double xtx[kLdscMaxP * kLdscMaxP];
  double xty[kLdscMaxP];
  for (uint32_t j = 0; j != p * p; ++j) {
    xtx[j] = 0.0;
  }
  for (uint32_t j = 0; j != p; ++j) {
    xty[j] = 0.0;
  }
  for (uint32_t i = 0; i != n; ++i) {
    const double ww = w[i] * w[i];
    const double* cur_x = &(x[S_CAST(uintptr_t, i) * p]);
    for (uint32_t j = 0; j != p; ++j) {
      xty[j] += ww * cur_x[j] * y[i];
      for (uint32_t k = 0; k != p; ++k) {
        xtx[j * p + k] += ww * cur_x[j] * cur_x[k];
      }
    }
  }
  return LdscSolve(p, xtx, xty, coef);
}

// Iteratively reweighted least squares, two updates, then the block jackknife
// on the finally-weighted design.  This is the reference implementation's
// IRWLS.irwls().
BoolErr LdscIrwls(const LdscRegCtx* ctx, const double* x, const double* y, const double* initial_w, uint32_t n, uint32_t p, const std::vector<uint32_t>& sep, LdscJknife* jkp) {
  std::vector<double> w(n);
  for (uint32_t i = 0; i != n; ++i) {
    if (initial_w[i] <= 0.0) {
      return 1;
    }
    w[i] = sqrt(initial_w[i]);
  }
  double coef[kLdscMaxP];
  std::vector<double> new_w(n);
  for (uint32_t iter = 0; iter != 2; ++iter) {
    if (LdscWls(x, y, &(w[0]), n, p, coef)) {
      return 1;
    }
    LdscUpdateWeights(ctx, coef, p, n, &(new_w[0]));
    for (uint32_t i = 0; i != n; ++i) {
      if (!(new_w[i] > 0.0)) {
        return 1;
      }
      w[i] = sqrt(new_w[i]);
    }
  }
  std::vector<double> xw(S_CAST(uintptr_t, n) * p);
  std::vector<double> yw(n);
  for (uint32_t i = 0; i != n; ++i) {
    for (uint32_t j = 0; j != p; ++j) {
      xw[S_CAST(uintptr_t, i) * p + j] = x[S_CAST(uintptr_t, i) * p + j] * w[i];
    }
    yw[i] = y[i] * w[i];
  }
  return LdscLstsqJknifeFast(&(xw[0]), &(yw[0]), p, sep, jkp);
}

// ***** the regressions *****

typedef struct LdscHsqResultStruct {
  double tot;
  double tot_se;
  double coef;
  double intercept;
  double intercept_se;
  double mean_chisq;
  double lambda_gc;
  double ratio;
  double ratio_se;
  uint32_t constrain_intercept;
  uint32_t ratio_valid;
  uint32_t n_blocks;
  std::vector<double> tot_delete_values;
} LdscHsqResult;

// Combines the free-intercept first step and the constrained-intercept second
// step of the two-step estimator into one jackknife, as
// LD_Score_Regression._combine_twostep_jknives().
void LdscCombineTwostep(const LdscJknife* step1, const LdscJknife* step2, double step1_int, double c, LdscJknife* out) {
  const uint32_t n_blocks = step1->n_blocks;
  out->n_blocks = n_blocks;
  out->p = 2;
  out->est[0] = step2->est[0];
  out->est[1] = step1_int;
  out->delete_values.assign(S_CAST(uintptr_t, n_blocks) * 2, 0.0);
  for (uint32_t b = 0; b != n_blocks; ++b) {
    const double step1_int_delete = step1->delete_values[b * 2 + 1];
    out->delete_values[b * 2] = step2->delete_values[b] - c * (step1_int_delete - step1_int);
    out->delete_values[b * 2 + 1] = step1_int_delete;
  }
  LdscFinishJknife(out);
}

// The shared LD Score regression body: aggregate estimate, initial weights,
// IRWLS (optionally two-step), then the jackknife.  y is chi^2 for h2 and
// z1*z2 for genetic covariance.
BoolErr LdscRegress(const LdscRegCtx* base_ctx, const double* y, uint32_t n_snp, uint32_t n_blocks, const double* fixed_intercept, const double* twostep, const double* step1_filter, double null_intercept, LdscJknife* jkp, double* nbar_ptr) {
  const double m_tot = base_ctx->m_tot;
  const double* ld = base_ctx->ld;
  const double* n_vec = base_ctx->n_vec;
  double nbar = 0.0;
  double y_sum = 0.0;
  double ldn_sum = 0.0;
  for (uint32_t i = 0; i != n_snp; ++i) {
    nbar += n_vec[i];
    y_sum += y[i];
    ldn_sum += ld[i] * n_vec[i];
  }
  nbar /= u31tod(n_snp);
  *nbar_ptr = nbar;
  const double intercept_for_agg = fixed_intercept? (*fixed_intercept) : null_intercept;
  const double tot_agg = m_tot * (y_sum / u31tod(n_snp) - intercept_for_agg) / (ldn_sum / u31tod(n_snp));

  std::vector<uint32_t> all_idxs(n_snp);
  for (uint32_t i = 0; i != n_snp; ++i) {
    all_idxs[i] = i;
  }
  LdscRegCtx agg_ctx = *base_ctx;
  agg_ctx.nbar = nbar;
  agg_ctx.row_idxs = &(all_idxs[0]);
  agg_ctx.constrain_intercept = 1;
  agg_ctx.fixed_intercept = intercept_for_agg;
  // aggregate estimate -> initial weights.  The weight update takes a
  // coefficient on the N-scaled scale, so invert that scaling here.
  double agg_coef[kLdscMaxP];
  agg_coef[0] = tot_agg * nbar / m_tot;
  agg_coef[1] = intercept_for_agg;
  std::vector<double> initial_w(n_snp);
  LdscUpdateWeights(&agg_ctx, agg_coef, 1, n_snp, &(initial_w[0]));

  // x is N-scaled to keep the condition number low.
  std::vector<double> x_scaled(n_snp);
  for (uint32_t i = 0; i != n_snp; ++i) {
    x_scaled[i] = n_vec[i] * ld[i] / nbar;
  }

  if (fixed_intercept) {
    // The two-step estimator exists to keep a free intercept from being
    // dragged around by the high-chi^2 variants; with the intercept
    // constrained there is nothing for it to do.  (The reference
    // implementation errors out on this combination instead.)
    std::vector<double> yp(n_snp);
    for (uint32_t i = 0; i != n_snp; ++i) {
      yp[i] = y[i] - (*fixed_intercept);
    }
    LdscRegCtx ctx = *base_ctx;
    ctx.nbar = nbar;
    ctx.row_idxs = &(all_idxs[0]);
    ctx.constrain_intercept = 1;
    ctx.fixed_intercept = *fixed_intercept;
    std::vector<uint32_t> sep;
    LdscGetSeparators(n_snp, n_blocks, &sep);
    return LdscIrwls(&ctx, &(x_scaled[0]), &(yp[0]), &(initial_w[0]), n_snp, 1, sep, jkp);
  }

  std::vector<double> design(S_CAST(uintptr_t, n_snp) * 2);
  for (uint32_t i = 0; i != n_snp; ++i) {
    design[S_CAST(uintptr_t, i) * 2] = x_scaled[i];
    design[S_CAST(uintptr_t, i) * 2 + 1] = 1.0;
  }

  if (!twostep) {
    LdscRegCtx ctx = *base_ctx;
    ctx.nbar = nbar;
    ctx.row_idxs = &(all_idxs[0]);
    ctx.constrain_intercept = 0;
    std::vector<uint32_t> sep;
    LdscGetSeparators(n_snp, n_blocks, &sep);
    return LdscIrwls(&ctx, &(design[0]), y, &(initial_w[0]), n_snp, 2, sep, jkp);
  }

  // Two-step: a free-intercept regression on the low-signal SNPs, then a
  // constrained-intercept regression on all of them.
  std::vector<uint32_t> step1_idxs;
  for (uint32_t i = 0; i != n_snp; ++i) {
    if (step1_filter[i] < *twostep) {
      step1_idxs.push_back(i);
    }
  }
  const uint32_t n1 = step1_idxs.size();
  if (n1 < n_blocks) {
    fprintf(stderr, "Error: only %u variants pass the two-step cutoff, fewer than the %u\njackknife blocks.  Raise --two-step, or lower --n-blocks.\n", n1, n_blocks);
    return 1;
  }
  std::vector<double> design1(S_CAST(uintptr_t, n1) * 2);
  std::vector<double> y1(n1);
  std::vector<double> initial_w1(n1);
  for (uint32_t i = 0; i != n1; ++i) {
    const uint32_t uidx = step1_idxs[i];
    design1[S_CAST(uintptr_t, i) * 2] = x_scaled[uidx];
    design1[S_CAST(uintptr_t, i) * 2 + 1] = 1.0;
    y1[i] = y[uidx];
    initial_w1[i] = initial_w[uidx];
  }
  LdscRegCtx ctx1 = *base_ctx;
  ctx1.nbar = nbar;
  ctx1.row_idxs = &(step1_idxs[0]);
  ctx1.constrain_intercept = 0;
  std::vector<uint32_t> sep1;
  LdscGetSeparators(n1, n_blocks, &sep1);
  LdscJknife step1_jk;
  if (LdscIrwls(&ctx1, &(design1[0]), &(y1[0]), &(initial_w1[0]), n1, 2, sep1, &step1_jk)) {
    return 1;
  }
  const double step1_int = step1_jk.est[1];

  std::vector<double> yp(n_snp);
  for (uint32_t i = 0; i != n_snp; ++i) {
    yp[i] = y[i] - step1_int;
  }
  // The second step keeps the first step's blocks: its separators are the
  // first step's, mapped back through the filter.
  std::vector<uint32_t> sep2(n_blocks + 1);
  sep2[0] = 0;
  sep2[n_blocks] = n_snp;
  for (uint32_t b = 1; b != n_blocks; ++b) {
    sep2[b] = step1_idxs[sep1[b]];
  }
  LdscRegCtx ctx2 = *base_ctx;
  ctx2.nbar = nbar;
  ctx2.row_idxs = &(all_idxs[0]);
  ctx2.constrain_intercept = 1;
  ctx2.fixed_intercept = step1_int;
  LdscJknife step2_jk;
  if (LdscIrwls(&ctx2, &(x_scaled[0]), &(yp[0]), &(initial_w[0]), n_snp, 1, sep2, &step2_jk)) {
    return 1;
  }
  double c_num = 0.0;
  double c_denom = 0.0;
  for (uint32_t i = 0; i != n_snp; ++i) {
    c_num += initial_w[i] * x_scaled[i];
    c_denom += initial_w[i] * x_scaled[i] * x_scaled[i];
  }
  LdscCombineTwostep(&step1_jk, &step2_jk, step1_int, c_num / c_denom, jkp);
  return 0;
}

BoolErr LdscHsqFit(const double* chisq, const double* ld, const double* w_ld, const double* n_vec, uint32_t n_snp, double m_tot, uint32_t n_blocks, const double* fixed_intercept, const double* twostep, LdscHsqResult* out) {
  LdscRegCtx ctx;
  ctx.kind = kLdscHsq;
  ctx.ld = ld;
  ctx.w_ld = w_ld;
  ctx.n_vec = n_vec;
  ctx.n1 = nullptr;
  ctx.n2 = nullptr;
  ctx.row_idxs = nullptr;
  ctx.m_tot = m_tot;
  ctx.nbar = 0.0;
  ctx.hsq1 = 0.0;
  ctx.hsq2 = 0.0;
  ctx.intercept_hsq1 = 1.0;
  ctx.intercept_hsq2 = 1.0;
  ctx.constrain_intercept = 0;
  ctx.fixed_intercept = 0.0;
  LdscJknife jk;
  double nbar;
  if (LdscRegress(&ctx, chisq, n_snp, n_blocks, fixed_intercept, twostep, chisq, 1.0, &jk, &nbar)) {
    return 1;
  }
  out->constrain_intercept = (fixed_intercept != nullptr);
  out->coef = jk.est[0] / nbar;
  out->tot = m_tot * out->coef;
  const double coef_var = jk.jknife_cov[0] / (nbar * nbar);
  out->tot_se = sqrt(m_tot * m_tot * coef_var);
  if (fixed_intercept) {
    out->intercept = *fixed_intercept;
    out->intercept_se = 0.0 / 0.0;
  } else {
    out->intercept = jk.est[1];
    out->intercept_se = jk.jknife_se[1];
  }
  out->n_blocks = jk.n_blocks;
  out->tot_delete_values.resize(jk.n_blocks);
  for (uint32_t b = 0; b != jk.n_blocks; ++b) {
    out->tot_delete_values[b] = jk.delete_values[S_CAST(uintptr_t, b) * jk.p] * m_tot / nbar;
  }
  double chisq_sum = 0.0;
  for (uint32_t i = 0; i != n_snp; ++i) {
    chisq_sum += chisq[i];
  }
  out->mean_chisq = chisq_sum / u31tod(n_snp);
  out->lambda_gc = LdscMedian(chisq, n_snp) / 0.4549;
  out->ratio_valid = 0;
  out->ratio = 0.0 / 0.0;
  out->ratio_se = 0.0 / 0.0;
  if ((!out->constrain_intercept) && (out->mean_chisq > 1.0)) {
    out->ratio_valid = 1;
    out->ratio = (out->intercept - 1.0) / (out->mean_chisq - 1.0);
    out->ratio_se = out->intercept_se / (out->mean_chisq - 1.0);
  }
  return 0;
}

typedef struct LdscGencovResultStruct {
  double tot;
  double tot_se;
  double intercept;
  double intercept_se;
  double mean_z1z2;
  double p;
  double z;
  uint32_t constrain_intercept;
  uint32_t n_blocks;
  std::vector<double> tot_delete_values;
} LdscGencovResult;

BoolErr LdscGencovFit(const double* z1, const double* z2, const double* ld, const double* w_ld, const double* n1, const double* n2, uint32_t n_snp, double m_tot, uint32_t n_blocks, double hsq1, double hsq2, double intercept_hsq1, double intercept_hsq2, const double* fixed_intercept, const double* twostep, LdscGencovResult* out) {
  std::vector<double> y(n_snp);
  std::vector<double> sqrt_n1n2(n_snp);
  std::vector<double> step1_filter(n_snp);
  for (uint32_t i = 0; i != n_snp; ++i) {
    y[i] = z1[i] * z2[i];
    sqrt_n1n2[i] = sqrt(n1[i] * n2[i]);
    // The reference implementation's step-1 filter is
    // (z1^2 < twostep) & (z2^2 < twostep), which is max(z1^2, z2^2).
    step1_filter[i] = MAXV(z1[i] * z1[i], z2[i] * z2[i]);
  }
  LdscRegCtx ctx;
  ctx.kind = kLdscGencov;
  ctx.ld = ld;
  ctx.w_ld = w_ld;
  ctx.n_vec = &(sqrt_n1n2[0]);
  ctx.n1 = n1;
  ctx.n2 = n2;
  ctx.row_idxs = nullptr;
  ctx.m_tot = m_tot;
  ctx.nbar = 0.0;
  ctx.hsq1 = hsq1;
  ctx.hsq2 = hsq2;
  ctx.intercept_hsq1 = intercept_hsq1;
  ctx.intercept_hsq2 = intercept_hsq2;
  ctx.constrain_intercept = 0;
  ctx.fixed_intercept = 0.0;
  LdscJknife jk;
  double nbar;
  if (LdscRegress(&ctx, &(y[0]), n_snp, n_blocks, fixed_intercept, twostep, &(step1_filter[0]), 0.0, &jk, &nbar)) {
    return 1;
  }
  out->constrain_intercept = (fixed_intercept != nullptr);
  const double coef = jk.est[0] / nbar;
  out->tot = m_tot * coef;
  const double coef_var = jk.jknife_cov[0] / (nbar * nbar);
  out->tot_se = sqrt(m_tot * m_tot * coef_var);
  if (fixed_intercept) {
    out->intercept = *fixed_intercept;
    out->intercept_se = 0.0 / 0.0;
  } else {
    out->intercept = jk.est[1];
    out->intercept_se = jk.jknife_se[1];
  }
  out->n_blocks = jk.n_blocks;
  out->tot_delete_values.resize(jk.n_blocks);
  for (uint32_t b = 0; b != jk.n_blocks; ++b) {
    out->tot_delete_values[b] = jk.delete_values[S_CAST(uintptr_t, b) * jk.p] * m_tot / nbar;
  }
  double acc = 0.0;
  for (uint32_t i = 0; i != n_snp; ++i) {
    acc += y[i];
  }
  out->mean_z1z2 = acc / u31tod(n_snp);
  out->z = out->tot / out->tot_se;
  out->p = LdscChisq1Sf(out->z);
  return 0;
}


// ***** input *****

typedef struct LdscSumstatRowStruct {
  double z;
  double n;
  char a1;
  char a2;
} LdscSumstatRow;

// Column lookup is case-insensitive and accepts the aliases ldsc and PLINK
// emit.
uint32_t LdscMatchCol(const char* tok, uint32_t slen, const char* const* names) {
  char buf[32];
  if (slen >= sizeof(buf)) {
    return 0;
  }
  for (uint32_t i = 0; i != slen; ++i) {
    const char c = tok[i];
    buf[i] = ((c >= 'a') && (c <= 'z'))? (c - 32) : c;
  }
  buf[slen] = '\0';
  for (uint32_t i = 0; names[i]; ++i) {
    if (!strcmp(buf, names[i])) {
      return 1;
    }
  }
  return 0;
}

char LdscAcgtComplement(char c) {
  switch (c) {
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  default: return '\0';
  }
}

char LdscUpcase(char c) {
  return ((c >= 'a') && (c <= 'z'))? (c - 32) : c;
}

// Reads a .sumstats file: SNP, Z, N, and (for genetic correlation) A1 and A2.
// Rows with a missing or unparsable value are dropped, as the reference
// implementation's dropna does.
BoolErr LdscReadSumstats(const char* fname, uint32_t need_alleles, std::unordered_map<std::string, LdscSumstatRow>* dst, uint32_t* dropped_ct_ptr) {
  static const char* kIdNames[] = {"SNP", "ID", "RSID", "MARKERNAME", "VARIANT_ID", nullptr};
  static const char* kA1Names[] = {"A1", "EFFECT_ALLELE", "ALLELE1", nullptr};
  static const char* kA2Names[] = {"A2", "OTHER_ALLELE", "ALLELE0", "ALLELE2", nullptr};
  static const char* kZNames[] = {"Z", "ZSCORE", "Z_SCORE", nullptr};
  static const char* kNNames[] = {"N", "OBS_CT", "NMISS", nullptr};

  TextStream txs;
  PreinitTextStream(&txs);
  PglErr reterr = TextStreamOpen(fname, &txs);
  if (reterr) {
    fprintf(stderr, "Error: Failed to open %s.\n", fname);
    return 1;
  }
  const char* header = TextGet(&txs);
  if (!header) {
    fprintf(stderr, "Error: %s is empty.\n", fname);
    return 1;
  }
  uint32_t col_id = UINT32_MAX;
  uint32_t col_a1 = UINT32_MAX;
  uint32_t col_a2 = UINT32_MAX;
  uint32_t col_z = UINT32_MAX;
  uint32_t col_n = UINT32_MAX;
  uint32_t col_ct = 0;
  {
    const char* iter = FirstNonTspace(header);
    if (*iter == '#') {
      ++iter;
    }
    for (; !IsEolnKns(*iter); ++col_ct) {
      const char* token_end = CurTokenEnd(iter);
      const uint32_t slen = token_end - iter;
      if ((col_id == UINT32_MAX) && LdscMatchCol(iter, slen, kIdNames)) {
        col_id = col_ct;
      } else if ((col_a1 == UINT32_MAX) && LdscMatchCol(iter, slen, kA1Names)) {
        col_a1 = col_ct;
      } else if ((col_a2 == UINT32_MAX) && LdscMatchCol(iter, slen, kA2Names)) {
        col_a2 = col_ct;
      } else if ((col_z == UINT32_MAX) && LdscMatchCol(iter, slen, kZNames)) {
        col_z = col_ct;
      } else if ((col_n == UINT32_MAX) && LdscMatchCol(iter, slen, kNNames)) {
        col_n = col_ct;
      }
      iter = FirstNonTspace(token_end);
    }
  }
  if ((col_id == UINT32_MAX) || (col_z == UINT32_MAX) || (col_n == UINT32_MAX)) {
    fprintf(stderr, "Error: %s must have SNP, Z and N columns.\n", fname);
    return 1;
  }
  if (need_alleles && ((col_a1 == UINT32_MAX) || (col_a2 == UINT32_MAX))) {
    fprintf(stderr, "Error: %s must have A1 and A2 columns for genetic-correlation\nestimation.\n", fname);
    return 1;
  }
  uint32_t max_col = MAXV(col_id, col_z);
  max_col = MAXV(max_col, col_n);
  if (need_alleles) {
    max_col = MAXV(max_col, MAXV(col_a1, col_a2));
  }
  uint32_t dropped_ct = 0;
  BoolErr ret = 0;
  while (1) {
    const char* line_start = TextGet(&txs);
    if (!line_start) {
      break;
    }
    const char* iter = FirstNonTspace(line_start);
    if (IsEolnKns(*iter)) {
      continue;
    }
    const char* id_start = nullptr;
    uint32_t id_slen = 0;
    LdscSumstatRow row;
    row.z = 0.0;
    row.n = 0.0;
    row.a1 = '\0';
    row.a2 = '\0';
    uint32_t ok = 1;
    for (uint32_t col_idx = 0; col_idx <= max_col; ++col_idx) {
      if (IsEolnKns(*iter)) {
        ok = 0;
        break;
      }
      const char* token_end = CurTokenEnd(iter);
      const uint32_t slen = token_end - iter;
      if (col_idx == col_id) {
        id_start = iter;
        id_slen = slen;
      }
      if (col_idx == col_z) {
        if (!ScanadvDouble(iter, &row.z)) {
          ok = 0;
          break;
        }
      }
      if (col_idx == col_n) {
        if ((!ScanadvDouble(iter, &row.n)) || (row.n <= 0.0)) {
          ok = 0;
          break;
        }
      }
      if (need_alleles && (col_idx == col_a1)) {
        if (slen != 1) {
          ok = 0;
          break;
        }
        row.a1 = LdscUpcase(*iter);
      }
      if (need_alleles && (col_idx == col_a2)) {
        if (slen != 1) {
          ok = 0;
          break;
        }
        row.a2 = LdscUpcase(*iter);
      }
      iter = FirstNonTspace(token_end);
    }
    if (!ok) {
      ++dropped_ct;
      continue;
    }
    dst->emplace(std::string(id_start, id_slen), row);
  }
  reterr = kPglRetSuccess;
  CleanupTextStream(&txs, &reterr);
  if (reterr) {
    fprintf(stderr, "Error: Failed to read %s.\n", fname);
    ret = 1;
  }
  *dropped_ct_ptr = dropped_ct;
  return ret;
}

// Reads one LD Score file, appending to ids/l2s in file order: that order is
// what the jackknife blocks are cut on, so it has to be preserved.
BoolErr LdscReadLdscoreFile(const char* fname, std::vector<std::string>* ids, std::vector<double>* l2s) {
  static const char* kIdNames[] = {"SNP", "ID", nullptr};
  static const char* kL2Names[] = {"L2", "LDSCORE", nullptr};
  TextStream txs;
  PreinitTextStream(&txs);
  PglErr reterr = TextStreamOpen(fname, &txs);
  if (reterr) {
    fprintf(stderr, "Error: Failed to open %s.\n", fname);
    return 1;
  }
  const char* header = TextGet(&txs);
  if (!header) {
    fprintf(stderr, "Error: %s is empty.\n", fname);
    return 1;
  }
  uint32_t col_id = UINT32_MAX;
  uint32_t col_l2 = UINT32_MAX;
  uint32_t col_ct = 0;
  {
    const char* iter = FirstNonTspace(header);
    if (*iter == '#') {
      ++iter;
    }
    for (; !IsEolnKns(*iter); ++col_ct) {
      const char* token_end = CurTokenEnd(iter);
      const uint32_t slen = token_end - iter;
      if ((col_id == UINT32_MAX) && LdscMatchCol(iter, slen, kIdNames)) {
        col_id = col_ct;
      } else if ((col_l2 == UINT32_MAX) && LdscMatchCol(iter, slen, kL2Names)) {
        col_l2 = col_ct;
      }
      iter = FirstNonTspace(token_end);
    }
  }
  if ((col_id == UINT32_MAX) || (col_l2 == UINT32_MAX)) {
    fprintf(stderr, "Error: %s must have SNP (or ID) and L2 columns.  Note that only\nunpartitioned LD Scores are supported.\n", fname);
    return 1;
  }
  const uint32_t max_col = MAXV(col_id, col_l2);
  BoolErr ret = 0;
  while (1) {
    const char* line_start = TextGet(&txs);
    if (!line_start) {
      break;
    }
    const char* iter = FirstNonTspace(line_start);
    if (IsEolnKns(*iter)) {
      continue;
    }
    const char* id_start = nullptr;
    uint32_t id_slen = 0;
    double l2 = 0.0;
    uint32_t ok = 1;
    for (uint32_t col_idx = 0; col_idx <= max_col; ++col_idx) {
      if (IsEolnKns(*iter)) {
        ok = 0;
        break;
      }
      const char* token_end = CurTokenEnd(iter);
      if (col_idx == col_id) {
        id_start = iter;
        id_slen = token_end - iter;
      }
      if (col_idx == col_l2) {
        if (!ScanadvDouble(iter, &l2)) {
          ok = 0;
          break;
        }
      }
      iter = FirstNonTspace(token_end);
    }
    if (!ok) {
      // NA LD Scores (monomorphic variants, in plink2's output) are dropped.
      continue;
    }
    ids->push_back(std::string(id_start, id_slen));
    l2s->push_back(l2);
  }
  reterr = kPglRetSuccess;
  CleanupTextStream(&txs, &reterr);
  if (reterr) {
    fprintf(stderr, "Error: Failed to read %s.\n", fname);
    ret = 1;
  }
  return ret;
}

// '@' is replaced by the chromosome number, as in the reference
// implementation; without one, the number is appended.
std::string LdscSubChr(const char* base, uint32_t chr_idx) {
  std::string s(base);
  char buf[16];
  snprintf(buf, sizeof(buf), "%u", chr_idx);
  const size_t at_pos = s.find('@');
  if (at_pos == std::string::npos) {
    return s + buf;
  }
  return s.substr(0, at_pos) + buf + s.substr(at_pos + 1);
}

uint32_t LdscFileExists(const std::string& path) {
  FILE* f = fopen(path.c_str(), FOPEN_RB);
  if (!f) {
    return 0;
  }
  fclose(f);
  return 1;
}

// Accepts both the reference implementation's naming (<base>.l2.ldscore, with
// or without a .gz) and a literal path, so plink2 --ld-score output can be
// handed over directly.
BoolErr LdscResolveLdscorePath(const std::string& base, std::string* out) {
  static const char* kSuffixes[] = {".l2.ldscore", ".l2.ldscore.gz", ".l2.ldscore.zst", "", ".ldscore", ".ldscore.gz", ".ldscore.zst", nullptr};
  for (uint32_t i = 0; kSuffixes[i]; ++i) {
    const std::string cand = base + kSuffixes[i];
    if (LdscFileExists(cand)) {
      *out = cand;
      return 0;
    }
  }
  fprintf(stderr, "Error: Could not find LD Scores at %s[.l2.ldscore/.gz/.zst].\n", base.c_str());
  return 1;
}

// Reads the LD Scores named by --ref-ld/--w-ld (a single fileset) or
// --ref-ld-chr/--w-ld-chr (one per chromosome, concatenated in chromosome
// order).
BoolErr LdscReadLdscores(const char* arg, uint32_t is_chr_split, std::vector<std::string>* ids, std::vector<double>* l2s) {
  if (!is_chr_split) {
    std::string path;
    if (LdscResolveLdscorePath(std::string(arg), &path)) {
      return 1;
    }
    return LdscReadLdscoreFile(path.c_str(), ids, l2s);
  }
  uint32_t found_ct = 0;
  for (uint32_t chr_idx = 1; chr_idx <= kLdscChrCt; ++chr_idx) {
    const std::string base = LdscSubChr(arg, chr_idx);
    std::string path;
    static const char* kSuffixes[] = {".l2.ldscore", ".l2.ldscore.gz", ".l2.ldscore.zst", ".ldscore", ".ldscore.gz", ".ldscore.zst", nullptr};
    uint32_t found = 0;
    for (uint32_t i = 0; kSuffixes[i]; ++i) {
      const std::string cand = base + kSuffixes[i];
      if (LdscFileExists(cand)) {
        path = cand;
        found = 1;
        break;
      }
    }
    if (!found) {
      continue;
    }
    if (LdscReadLdscoreFile(path.c_str(), ids, l2s)) {
      return 1;
    }
    ++found_ct;
  }
  if (!found_ct) {
    fprintf(stderr, "Error: No LD Score files found for %s (expected\n%s<chr>.l2.ldscore[.gz]).\n", arg, arg);
    return 1;
  }
  return 0;
}

// Reads the .l2.M_5_50 (or .l2.M) files holding the number of variants the LD
// Scores were computed from.
BoolErr LdscReadM(const char* arg, uint32_t is_chr_split, uint32_t not_m_5_50, double* m_ptr) {
  const char* suffix = not_m_5_50? ".l2.M" : ".l2.M_5_50";
  double acc = 0.0;
  uint32_t found_ct = 0;
  const uint32_t chr_end = is_chr_split? (kLdscChrCt + 1) : 1;
  for (uint32_t chr_idx = 0; chr_idx != chr_end; ++chr_idx) {
    std::string path;
    if (is_chr_split) {
      path = LdscSubChr(arg, chr_idx + 1) + suffix;
    } else {
      path = std::string(arg) + suffix;
    }
    if (!LdscFileExists(path)) {
      continue;
    }
    FILE* f = fopen(path.c_str(), FOPEN_RB);
    if (!f) {
      continue;
    }
    char buf[256];
    if (!fgets(buf, sizeof(buf), f)) {
      fclose(f);
      continue;
    }
    fclose(f);
    double cur;
    if (!ScanadvDouble(buf, &cur)) {
      fprintf(stderr, "Error: Malformed %s.\n", path.c_str());
      return 1;
    }
    // A second value would mean partitioned LD Scores.
    const char* iter = FirstNonTspace(buf);
    iter = FirstNonTspace(CurTokenEnd(iter));
    if (!IsEolnKns(*iter)) {
      fprintf(stderr, "Error: %s has more than one entry; only unpartitioned LD Scores are\nsupported.\n", path.c_str());
      return 1;
    }
    acc += cur;
    ++found_ct;
  }
  if (!found_ct) {
    return 1;
  }
  *m_ptr = acc;
  return 0;
}

// ***** logging *****

FILE* g_ldsc_logfile = nullptr;

void LdscLog(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  char buf[4096];
  vsnprintf(buf, sizeof(buf), fmt, args);
  va_end(args);
  fputs(buf, stdout);
  if (g_ldsc_logfile) {
    fputs(buf, g_ldsc_logfile);
  }
}

// The reference implementation prints these through numpy, which gives four
// significant digits.
const char* LdscFmt(double val, char* buf, uintptr_t buf_size) {
  if (val != val) {
    snprintf(buf, buf_size, "NA");
  } else {
    snprintf(buf, buf_size, "%.4g", val);
  }
  return buf;
}

// ***** merging *****

// One row per variant that made it through every merge, in LD Score file
// order.
typedef struct LdscDataStruct {
  std::vector<double> ld;
  std::vector<double> w_ld;
  std::vector<double> z1;
  std::vector<double> n1;
  std::vector<double> z2;
  std::vector<double> n2;
  std::vector<char> a1;
  std::vector<char> a2;
  std::vector<std::string> ids;
} LdscData;

// Inner-joins the LD Scores, the regression weight LD Scores and one trait's
// summary statistics, keeping the LD Score file's order.
void LdscMerge(const std::vector<std::string>& ref_ids, const std::vector<double>& ref_l2, const std::unordered_map<std::string, double>& w_ld_map, const std::unordered_map<std::string, LdscSumstatRow>& sumstats, uint32_t keep_alleles, LdscData* dst) {
  const uintptr_t ref_ct = ref_ids.size();
  for (uintptr_t i = 0; i != ref_ct; ++i) {
    const std::unordered_map<std::string, LdscSumstatRow>::const_iterator ss_it = sumstats.find(ref_ids[i]);
    if (ss_it == sumstats.end()) {
      continue;
    }
    const std::unordered_map<std::string, double>::const_iterator w_it = w_ld_map.find(ref_ids[i]);
    if (w_it == w_ld_map.end()) {
      continue;
    }
    dst->ids.push_back(ref_ids[i]);
    dst->ld.push_back(ref_l2[i]);
    dst->w_ld.push_back(w_it->second);
    dst->z1.push_back(ss_it->second.z);
    dst->n1.push_back(ss_it->second.n);
    if (keep_alleles) {
      dst->a1.push_back(ss_it->second.a1);
      dst->a2.push_back(ss_it->second.a2);
    }
  }
}

uint32_t LdscIsValidSnp(char a1, char a2) {
  const char c1 = LdscAcgtComplement(a1);
  const char c2 = LdscAcgtComplement(a2);
  if ((!c1) || (!c2) || (a1 == a2)) {
    return 0;
  }
  // Strand-ambiguous variants are dropped: their orientation cannot be
  // recovered, so a sign error would be silent.
  return (c1 != a2);
}

// ***** output *****

void LdscPrintHsq(const LdscHsqResult* hsq, const char* label, const double* samp_prev, const double* pop_prev) {
  char buf1[64];
  char buf2[64];
  double c = 1.0;
  const char* scale = "Observed";
  if (samp_prev && pop_prev) {
    c = LdscLiabilityFactor(*samp_prev, *pop_prev);
    scale = "Liability";
  }
  if (label) {
    LdscLog("\n%s\n", label);
  }
  LdscLog("Total %s scale h2: %s (%s)\n", scale, LdscFmt(c * hsq->tot, buf1, sizeof(buf1)), LdscFmt(c * hsq->tot_se, buf2, sizeof(buf2)));
  LdscLog("Lambda GC: %s\n", LdscFmt(hsq->lambda_gc, buf1, sizeof(buf1)));
  LdscLog("Mean Chi^2: %s\n", LdscFmt(hsq->mean_chisq, buf1, sizeof(buf1)));
  if (hsq->constrain_intercept) {
    LdscLog("Intercept: constrained to %s\n", LdscFmt(hsq->intercept, buf1, sizeof(buf1)));
    return;
  }
  LdscLog("Intercept: %s (%s)\n", LdscFmt(hsq->intercept, buf1, sizeof(buf1)), LdscFmt(hsq->intercept_se, buf2, sizeof(buf2)));
  if (!hsq->ratio_valid) {
    LdscLog("Ratio: NA (mean chi^2 < 1)\n");
  } else if (hsq->ratio < 0.0) {
    LdscLog("Ratio < 0 (usually indicates GC correction).\n");
  } else {
    LdscLog("Ratio: %s (%s)\n", LdscFmt(hsq->ratio, buf1, sizeof(buf1)), LdscFmt(hsq->ratio_se, buf2, sizeof(buf2)));
  }
}

void LdscPrintGencov(const LdscGencovResult* gencov, const double* samp_prev1, const double* pop_prev1, const double* samp_prev2, const double* pop_prev2) {
  char buf1[64];
  char buf2[64];
  double c = 1.0;
  const char* scale = "Observed";
  if (samp_prev1 && pop_prev1 && samp_prev2 && pop_prev2) {
    c = sqrt(LdscLiabilityFactor(*samp_prev1, *pop_prev1)) * sqrt(LdscLiabilityFactor(*samp_prev2, *pop_prev2));
    scale = "Liability";
  }
  LdscLog("\nGenetic Covariance\n");
  LdscLog("Total %s scale gencov: %s (%s)\n", scale, LdscFmt(c * gencov->tot, buf1, sizeof(buf1)), LdscFmt(c * gencov->tot_se, buf2, sizeof(buf2)));
  LdscLog("Mean z1*z2: %s\n", LdscFmt(gencov->mean_z1z2, buf1, sizeof(buf1)));
  if (gencov->constrain_intercept) {
    LdscLog("Intercept: constrained to %s\n", LdscFmt(gencov->intercept, buf1, sizeof(buf1)));
  } else {
    LdscLog("Intercept: %s (%s)\n", LdscFmt(gencov->intercept, buf1, sizeof(buf1)), LdscFmt(gencov->intercept_se, buf2, sizeof(buf2)));
  }
}

// ***** drivers *****

typedef struct LdscOptsStruct {
  uint32_t n_blocks;
  uint32_t have_chisq_max;
  double chisq_max;
  uint32_t have_two_step;
  double two_step;
  uint32_t no_intercept;
  uint32_t not_m_5_50;
  uint32_t no_check_alleles;
  double m_override;
  std::vector<double> intercept_h2;     // NaN where unspecified
  std::vector<double> intercept_gencov;
  std::vector<double> samp_prev;
  std::vector<double> pop_prev;
} LdscOpts;

const double* LdscOptAt(const std::vector<double>& vals, uint32_t idx) {
  if (idx >= vals.size()) {
    return nullptr;
  }
  if (vals[idx] != vals[idx]) {
    return nullptr;
  }
  return &(vals[idx]);
}

// Applies --chisq-max, then reports how many variants it removed, as the
// reference implementation does.
void LdscFilterChisq(double chisq_max, LdscData* data, uint32_t use_both_z) {
  const uint32_t orig_ct = data->ld.size();
  uint32_t write_idx = 0;
  for (uint32_t i = 0; i != orig_ct; ++i) {
    uint32_t keep;
    if (!use_both_z) {
      keep = (data->z1[i] * data->z1[i] < chisq_max);
    } else {
      const double prod = data->z1[i] * data->z1[i] * data->z2[i] * data->z2[i];
      keep = (prod < chisq_max * chisq_max);
    }
    if (!keep) {
      continue;
    }
    data->ld[write_idx] = data->ld[i];
    data->w_ld[write_idx] = data->w_ld[i];
    data->z1[write_idx] = data->z1[i];
    data->n1[write_idx] = data->n1[i];
    if (use_both_z) {
      data->z2[write_idx] = data->z2[i];
      data->n2[write_idx] = data->n2[i];
    }
    ++write_idx;
  }
  data->ld.resize(write_idx);
  data->w_ld.resize(write_idx);
  data->z1.resize(write_idx);
  data->n1.resize(write_idx);
  if (use_both_z) {
    data->z2.resize(write_idx);
    data->n2.resize(write_idx);
  }
  LdscLog("Removed %u variants with chi^2 > %g (%u remain).\n", orig_ct - write_idx, chisq_max, write_idx);
}

void LdscWarnLength(uint32_t n_snp) {
  if (n_snp < 200000) {
    LdscLog("WARNING: number of variants less than 200k; this is almost always bad.\nSuch a small number of variants is usually a sign of a data-munging problem.\n");
  }
}

BoolErr LdscEstimateH2(const char* sumstats_fname, const std::vector<std::string>& ref_ids, const std::vector<double>& ref_l2, const std::unordered_map<std::string, double>& w_ld_map, double m_tot, const LdscOpts* opts, const char* out_prefix) {
  std::unordered_map<std::string, LdscSumstatRow> sumstats;
  uint32_t dropped_ct;
  if (LdscReadSumstats(sumstats_fname, 0, &sumstats, &dropped_ct)) {
    return 1;
  }
  LdscLog("Read summary statistics for %" PRIuPTR " variants from %s", S_CAST(uintptr_t, sumstats.size()), sumstats_fname);
  if (dropped_ct) {
    LdscLog(" (%u incomplete rows dropped)", dropped_ct);
  }
  LdscLog(".\n");
  LdscData data;
  LdscMerge(ref_ids, ref_l2, w_ld_map, sumstats, 0, &data);
  uint32_t n_snp = data.ld.size();
  if (!n_snp) {
    fprintf(stderr, "Error: No variants remain after merging the summary statistics with the LD\nScores.\n");
    return 1;
  }
  LdscLog("After merging with reference panel LD and regression weight LD, %u variants\nremain.\n", n_snp);
  const double* fixed_intercept = nullptr;
  double fixed_intercept_val = 1.0;
  if (opts->no_intercept) {
    fixed_intercept = &fixed_intercept_val;
  } else {
    const double* cur = LdscOptAt(opts->intercept_h2, 0);
    if (cur) {
      fixed_intercept_val = *cur;
      fixed_intercept = &fixed_intercept_val;
    }
  }
  if (opts->have_chisq_max) {
    LdscFilterChisq(opts->chisq_max, &data, 0);
    n_snp = data.ld.size();
    if (!n_snp) {
      fprintf(stderr, "Error: --chisq-max removed every variant.\n");
      return 1;
    }
  }
  LdscWarnLength(n_snp);
  double two_step_val = 30.0;
  const double* two_step = nullptr;
  if (opts->have_two_step) {
    two_step_val = opts->two_step;
    two_step = &two_step_val;
  } else if (!fixed_intercept) {
    two_step = &two_step_val;
  }
  if (two_step && fixed_intercept) {
    LdscLog("Ignoring --two-step: it only applies when the intercept is free.\n");
    two_step = nullptr;
  } else if (two_step) {
    LdscLog("Using two-step estimator with cutoff at %g.\n", *two_step);
  }
  const uint32_t n_blocks = MINV(opts->n_blocks, n_snp);
  std::vector<double> chisq(n_snp);
  for (uint32_t i = 0; i != n_snp; ++i) {
    chisq[i] = data.z1[i] * data.z1[i];
  }
  LdscHsqResult hsq;
  if (LdscHsqFit(&(chisq[0]), &(data.ld[0]), &(data.w_ld[0]), &(data.n1[0]), n_snp, m_tot, n_blocks, fixed_intercept, two_step, &hsq)) {
    fprintf(stderr, "Error: Heritability regression failed (singular design, or nonpositive\nregression weights).\n");
    return 1;
  }
  LdscPrintHsq(&hsq, "Heritability of phenotype 1", LdscOptAt(opts->samp_prev, 0), LdscOptAt(opts->pop_prev, 0));

  // The same numbers at full precision, for programmatic use.
  const std::string h2_path = std::string(out_prefix) + ".h2";
  FILE* h2_file = fopen(h2_path.c_str(), FOPEN_WB);
  if (!h2_file) {
    fprintf(stderr, "Error: Failed to open %s.\n", h2_path.c_str());
    return 1;
  }
  double liab_factor = 1.0;
  const double* samp_prev = LdscOptAt(opts->samp_prev, 0);
  const double* pop_prev = LdscOptAt(opts->pop_prev, 0);
  if (samp_prev && pop_prev) {
    liab_factor = LdscLiabilityFactor(*samp_prev, *pop_prev);
  }
  fputs("h2\th2_se\tintercept\tintercept_se\tratio\tratio_se\tmean_chisq\tlambda_gc\tn_snp\tn_blocks\tm\n", h2_file);
  fprintf(h2_file, "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%u\t%u\t%.9g\n", liab_factor * hsq.tot, liab_factor * hsq.tot_se, hsq.intercept, hsq.intercept_se, hsq.ratio, hsq.ratio_se, hsq.mean_chisq, hsq.lambda_gc, n_snp, n_blocks, m_tot);
  if (fclose(h2_file)) {
    fprintf(stderr, "Error: Failed to write %s.\n", h2_path.c_str());
    return 1;
  }
  LdscLog("\nResults written to %s .\n", h2_path.c_str());
  return 0;
}

// Merges a second trait into the first trait's rows, flipping its Z where the
// alleles are swapped and dropping variants whose alleles do not match.
BoolErr LdscMergeSecondTrait(const LdscData* base, const std::unordered_map<std::string, LdscSumstatRow>& sumstats2, uint32_t no_check_alleles, LdscData* dst) {
  const uint32_t n_base = base->ld.size();
  uint32_t bad_allele_ct = 0;
  for (uint32_t i = 0; i != n_base; ++i) {
    const std::unordered_map<std::string, LdscSumstatRow>::const_iterator it = sumstats2.find(base->ids[i]);
    if (it == sumstats2.end()) {
      continue;
    }
    const char a1 = base->a1[i];
    const char a2 = base->a2[i];
    const char b1 = it->second.a1;
    const char b2 = it->second.a2;
    double z2 = it->second.z;
    if (!no_check_alleles) {
      if ((!LdscIsValidSnp(a1, a2)) || (!LdscIsValidSnp(b1, b2))) {
        ++bad_allele_ct;
        continue;
      }
      const char c1 = LdscAcgtComplement(b1);
      const char c2 = LdscAcgtComplement(b2);
      if ((a1 == b1) && (a2 == b2)) {
        // same orientation
      } else if ((a1 == c1) && (a2 == c2)) {
        // strand flip, same effect allele
      } else if ((a1 == b2) && (a2 == b1)) {
        z2 = -z2;
      } else if ((a1 == c2) && (a2 == c1)) {
        z2 = -z2;
      } else {
        ++bad_allele_ct;
        continue;
      }
    }
    dst->ids.push_back(base->ids[i]);
    dst->ld.push_back(base->ld[i]);
    dst->w_ld.push_back(base->w_ld[i]);
    dst->z1.push_back(base->z1[i]);
    dst->n1.push_back(base->n1[i]);
    dst->z2.push_back(z2);
    dst->n2.push_back(it->second.n);
  }
  if (bad_allele_ct) {
    LdscLog("Dropped %u variants with mismatched or strand-ambiguous alleles.\n", bad_allele_ct);
  }
  if (dst->ld.empty()) {
    fprintf(stderr, "Error: No variants in common between the two summary statistic files.\n");
    return 1;
  }
  return 0;
}

BoolErr LdscEstimateRg(const std::vector<std::string>& rg_paths, const std::vector<std::string>& ref_ids, const std::vector<double>& ref_l2, const std::unordered_map<std::string, double>& w_ld_map, double m_tot, const LdscOpts* opts, const char* out_prefix) {
  const uint32_t pheno_ct = rg_paths.size();
  std::unordered_map<std::string, LdscSumstatRow> sumstats1;
  uint32_t dropped_ct;
  if (LdscReadSumstats(rg_paths[0].c_str(), 1, &sumstats1, &dropped_ct)) {
    return 1;
  }
  LdscLog("Read summary statistics for %" PRIuPTR " variants from %s.\n", S_CAST(uintptr_t, sumstats1.size()), rg_paths[0].c_str());
  LdscData base;
  LdscMerge(ref_ids, ref_l2, w_ld_map, sumstats1, 1, &base);
  if (base.ld.empty()) {
    fprintf(stderr, "Error: No variants remain after merging %s with the LD Scores.\n", rg_paths[0].c_str());
    return 1;
  }
  LdscLog("After merging with reference panel LD and regression weight LD, %" PRIuPTR "\nvariants remain.\n", S_CAST(uintptr_t, base.ld.size()));

  double two_step_val = 30.0;
  const double* two_step = nullptr;
  const uint32_t intercept_h2_free = (!opts->no_intercept) && (!LdscOptAt(opts->intercept_h2, 0));
  if (opts->have_two_step) {
    two_step_val = opts->two_step;
    two_step = &two_step_val;
  } else if (intercept_h2_free) {
    two_step = &two_step_val;
  }
  if (two_step) {
    LdscLog("Using two-step estimator with cutoff at %g.  It applies to whichever of\nthe three regressions has a free intercept.\n", *two_step);
  }

  std::vector<double> rg_vals(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> rg_ses(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> rg_zs(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> rg_ps(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> h2_2_vals(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> h2_2_ses(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> h2_2_ints(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> h2_2_int_ses(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> gcov_vals(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> gcov_ses(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> gcov_ints(pheno_ct - 1, 0.0 / 0.0);
  std::vector<double> gcov_int_ses(pheno_ct - 1, 0.0 / 0.0);
  double h2_1_val = 0.0 / 0.0;
  double h2_1_se = 0.0 / 0.0;
  for (uint32_t pheno_idx = 1; pheno_idx != pheno_ct; ++pheno_idx) {
    LdscLog("\nComputing rg for phenotype %u/%u\n", pheno_idx + 1, pheno_ct);
    std::unordered_map<std::string, LdscSumstatRow> sumstats2;
    if (LdscReadSumstats(rg_paths[pheno_idx].c_str(), 1, &sumstats2, &dropped_ct)) {
      return 1;
    }
    LdscData data;
    if (LdscMergeSecondTrait(&base, sumstats2, opts->no_check_alleles, &data)) {
      return 1;
    }
    uint32_t n_snp = data.ld.size();
    LdscLog("%u variants with valid alleles in both files.\n", n_snp);
    if (opts->have_chisq_max) {
      LdscFilterChisq(opts->chisq_max, &data, 1);
      n_snp = data.ld.size();
      if (!n_snp) {
        fprintf(stderr, "Error: --chisq-max removed every variant.\n");
        return 1;
      }
    }
    LdscWarnLength(n_snp);
    const uint32_t n_blocks = MINV(opts->n_blocks, n_snp);

    double h2_int1 = 1.0;
    double h2_int2 = 1.0;
    const double* fixed_h2_1 = nullptr;
    const double* fixed_h2_2 = nullptr;
    double gencov_int = 0.0;
    const double* fixed_gencov = nullptr;
    if (opts->no_intercept) {
      fixed_h2_1 = &h2_int1;
      fixed_h2_2 = &h2_int2;
      fixed_gencov = &gencov_int;
    } else {
      const double* cur = LdscOptAt(opts->intercept_h2, 0);
      if (cur) {
        h2_int1 = *cur;
        fixed_h2_1 = &h2_int1;
      }
      cur = LdscOptAt(opts->intercept_h2, pheno_idx);
      if (cur) {
        h2_int2 = *cur;
        fixed_h2_2 = &h2_int2;
      }
      cur = LdscOptAt(opts->intercept_gencov, pheno_idx);
      if (cur) {
        gencov_int = *cur;
        fixed_gencov = &gencov_int;
      }
    }

    std::vector<double> chisq1(n_snp);
    std::vector<double> chisq2(n_snp);
    for (uint32_t i = 0; i != n_snp; ++i) {
      chisq1[i] = data.z1[i] * data.z1[i];
      chisq2[i] = data.z2[i] * data.z2[i];
    }
    LdscHsqResult hsq1;
    LdscHsqResult hsq2;
    if (LdscHsqFit(&(chisq1[0]), &(data.ld[0]), &(data.w_ld[0]), &(data.n1[0]), n_snp, m_tot, n_blocks, fixed_h2_1, two_step, &hsq1) ||
        LdscHsqFit(&(chisq2[0]), &(data.ld[0]), &(data.w_ld[0]), &(data.n2[0]), n_snp, m_tot, n_blocks, fixed_h2_2, two_step, &hsq2)) {
      fprintf(stderr, "Error: Heritability regression failed for phenotype pair 1/%u.\n", pheno_idx + 1);
      return 1;
    }
    LdscGencovResult gencov;
    if (LdscGencovFit(&(data.z1[0]), &(data.z2[0]), &(data.ld[0]), &(data.w_ld[0]), &(data.n1[0]), &(data.n2[0]), n_snp, m_tot, n_blocks, hsq1.tot, hsq2.tot, hsq1.intercept, hsq2.intercept, fixed_gencov, two_step, &gencov)) {
      fprintf(stderr, "Error: Genetic covariance regression failed for phenotype pair 1/%u.\n", pheno_idx + 1);
      return 1;
    }
    if (pheno_idx == 1) {
      LdscPrintHsq(&hsq1, "Heritability of phenotype 1", LdscOptAt(opts->samp_prev, 0), LdscOptAt(opts->pop_prev, 0));
      h2_1_val = hsq1.tot;
      h2_1_se = hsq1.tot_se;
    }
    h2_2_vals[pheno_idx - 1] = hsq2.tot;
    h2_2_ses[pheno_idx - 1] = hsq2.tot_se;
    h2_2_ints[pheno_idx - 1] = hsq2.intercept;
    h2_2_int_ses[pheno_idx - 1] = hsq2.intercept_se;
    gcov_vals[pheno_idx - 1] = gencov.tot;
    gcov_ses[pheno_idx - 1] = gencov.tot_se;
    gcov_ints[pheno_idx - 1] = gencov.intercept;
    gcov_int_ses[pheno_idx - 1] = gencov.intercept_se;
    char label[128];
    snprintf(label, sizeof(label), "Heritability of phenotype %u/%u", pheno_idx + 1, pheno_ct);
    LdscPrintHsq(&hsq2, label, LdscOptAt(opts->samp_prev, pheno_idx), LdscOptAt(opts->pop_prev, pheno_idx));
    LdscPrintGencov(&gencov, LdscOptAt(opts->samp_prev, 0), LdscOptAt(opts->pop_prev, 0), LdscOptAt(opts->samp_prev, pheno_idx), LdscOptAt(opts->pop_prev, pheno_idx));

    LdscLog("\nGenetic Correlation\n");
    char buf1[64];
    char buf2[64];
    if ((hsq1.tot <= 0.0) || (hsq2.tot <= 0.0)) {
      LdscLog("Genetic Correlation: NA (h2 out of bounds)\nWARNING: One of the h2's was out of bounds.  This usually indicates a\ndata-munging error or that h2 or N is low.\n");
      continue;
    }
    const double rg_ratio = gencov.tot / sqrt(hsq1.tot * hsq2.tot);
    std::vector<double> denom_delete(n_blocks);
    uint32_t denom_ok = 1;
    for (uint32_t b = 0; b != n_blocks; ++b) {
      const double prod = hsq1.tot_delete_values[b] * hsq2.tot_delete_values[b];
      if (prod <= 0.0) {
        denom_ok = 0;
        break;
      }
      denom_delete[b] = sqrt(prod);
    }
    if (!denom_ok) {
      LdscLog("Genetic Correlation: %s (NA) (a jackknife block had nonpositive h2)\n", LdscFmt(rg_ratio, buf1, sizeof(buf1)));
      rg_vals[pheno_idx - 1] = rg_ratio;
      continue;
    }
    double rg_jknife;
    double rg_se;
    LdscRatioJknife(rg_ratio, &(gencov.tot_delete_values[0]), &(denom_delete[0]), n_blocks, &rg_jknife, &rg_se);
    const double z = rg_ratio / rg_se;
    const double p = LdscChisq1Sf(z);
    rg_vals[pheno_idx - 1] = rg_ratio;
    rg_ses[pheno_idx - 1] = rg_se;
    rg_zs[pheno_idx - 1] = z;
    rg_ps[pheno_idx - 1] = p;
    if ((rg_ratio > 1.2) || (rg_ratio < -1.2)) {
      LdscLog("Genetic Correlation: %s (%s)\nWARNING: rg was out of bounds.  This often means that h2 is not\nsignificantly different from zero.\n", LdscFmt(rg_ratio, buf1, sizeof(buf1)), LdscFmt(rg_se, buf2, sizeof(buf2)));
    } else {
      LdscLog("Genetic Correlation: %s (%s)\n", LdscFmt(rg_ratio, buf1, sizeof(buf1)), LdscFmt(rg_se, buf2, sizeof(buf2)));
      LdscLog("Z-score: %s\n", LdscFmt(z, buf1, sizeof(buf1)));
      LdscLog("P: %s\n", LdscFmt(p, buf1, sizeof(buf1)));
    }
  }

  // The table the reference implementation prints at the end, also written to
  // <out>.rg for programmatic use.
  LdscLog("\nSummary of Genetic Correlation Results\n");
  LdscLog("p1\tp2\trg\tse\tz\tp\n");
  std::string rg_path = std::string(out_prefix) + ".rg";
  FILE* rg_file = fopen(rg_path.c_str(), FOPEN_WB);
  if (!rg_file) {
    fprintf(stderr, "Error: Failed to open %s.\n", rg_path.c_str());
    return 1;
  }
  fputs("p1\tp2\trg\tse\tz\tp\th2_p1\th2_p1_se\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov\tgcov_se\tgcov_int\tgcov_int_se\n", rg_file);
  for (uint32_t pheno_idx = 1; pheno_idx != pheno_ct; ++pheno_idx) {
    char buf1[64];
    char buf2[64];
    char buf3[64];
    char buf4[64];
    LdscLog("%s\t%s\t%s\t%s\t%s\t%s\n", rg_paths[0].c_str(), rg_paths[pheno_idx].c_str(), LdscFmt(rg_vals[pheno_idx - 1], buf1, sizeof(buf1)), LdscFmt(rg_ses[pheno_idx - 1], buf2, sizeof(buf2)), LdscFmt(rg_zs[pheno_idx - 1], buf3, sizeof(buf3)), LdscFmt(rg_ps[pheno_idx - 1], buf4, sizeof(buf4)));
    fprintf(rg_file, "%s\t%s\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\n", rg_paths[0].c_str(), rg_paths[pheno_idx].c_str(), rg_vals[pheno_idx - 1], rg_ses[pheno_idx - 1], rg_zs[pheno_idx - 1], rg_ps[pheno_idx - 1], h2_1_val, h2_1_se, h2_2_vals[pheno_idx - 1], h2_2_ses[pheno_idx - 1], h2_2_ints[pheno_idx - 1], h2_2_int_ses[pheno_idx - 1], gcov_vals[pheno_idx - 1], gcov_ses[pheno_idx - 1], gcov_ints[pheno_idx - 1], gcov_int_ses[pheno_idx - 1]);
  }
  if (fclose(rg_file)) {
    fprintf(stderr, "Error: Failed to write %s.\n", rg_path.c_str());
    return 1;
  }
  LdscLog("\nResults written to %s .\n", rg_path.c_str());
  return 0;
}

// Splits a comma-separated list.
void LdscSplitComma(const char* arg, std::vector<std::string>* dst) {
  const char* iter = arg;
  while (1) {
    const char* comma = strchr(iter, ',');
    if (!comma) {
      dst->push_back(std::string(iter));
      return;
    }
    dst->push_back(std::string(iter, comma - iter));
    iter = comma + 1;
  }
}

// Parses a comma-separated list of numbers, where 'N' means "not specified"
// (the reference implementation's convention for --samp-prev etc.).
BoolErr LdscParseNumList(const char* arg, const char* flagname, std::vector<double>* dst) {
  std::vector<std::string> parts;
  LdscSplitComma(arg, &parts);
  for (uintptr_t i = 0; i != parts.size(); ++i) {
    if ((parts[i] == "N") || (parts[i] == "n") || parts[i].empty()) {
      dst->push_back(0.0 / 0.0);
      continue;
    }
    double cur;
    if (!ScanadvDouble(parts[i].c_str(), &cur)) {
      fprintf(stderr, "Error: Invalid %s value '%s'.\n", flagname, parts[i].c_str());
      return 1;
    }
    dst->push_back(cur);
  }
  return 0;
}

#ifdef __cplusplus
}  // namespace plink2
#endif

int main(int argc, char** argv) {
  using namespace plink2;
  const char* h2_fname = nullptr;
  const char* rg_arg = nullptr;
  const char* ref_ld_arg = nullptr;
  uint32_t ref_ld_chr_split = 0;
  const char* w_ld_arg = nullptr;
  uint32_t w_ld_chr_split = 0;
  const char* out_prefix = nullptr;
  LdscOpts opts;
  opts.n_blocks = 200;
  opts.have_chisq_max = 0;
  opts.chisq_max = 0.0;
  opts.have_two_step = 0;
  opts.two_step = 30.0;
  opts.no_intercept = 0;
  opts.not_m_5_50 = 0;
  opts.no_check_alleles = 0;
  opts.m_override = 0.0;
  for (int argi = 1; argi < argc; ++argi) {
    const char* cur = argv[argi];
    if ((!strcmp(cur, "--h2")) && (argi + 1 < argc)) {
      h2_fname = argv[++argi];
    } else if ((!strcmp(cur, "--rg")) && (argi + 1 < argc)) {
      rg_arg = argv[++argi];
    } else if ((!strcmp(cur, "--ref-ld")) && (argi + 1 < argc)) {
      ref_ld_arg = argv[++argi];
      ref_ld_chr_split = 0;
    } else if ((!strcmp(cur, "--ref-ld-chr")) && (argi + 1 < argc)) {
      ref_ld_arg = argv[++argi];
      ref_ld_chr_split = 1;
    } else if ((!strcmp(cur, "--w-ld")) && (argi + 1 < argc)) {
      w_ld_arg = argv[++argi];
      w_ld_chr_split = 0;
    } else if ((!strcmp(cur, "--w-ld-chr")) && (argi + 1 < argc)) {
      w_ld_arg = argv[++argi];
      w_ld_chr_split = 1;
    } else if ((!strcmp(cur, "--out")) && (argi + 1 < argc)) {
      out_prefix = argv[++argi];
    } else if ((!strcmp(cur, "--M")) && (argi + 1 < argc)) {
      const char* m_str = argv[++argi];
      if ((!ScanadvDouble(m_str, &opts.m_override)) || (opts.m_override <= 0.0)) {
        fprintf(stderr, "Error: --M argument '%s' is not a positive number.\n", m_str);
        return 1;
      }
    } else if (!strcmp(cur, "--not-M-5-50")) {
      opts.not_m_5_50 = 1;
    } else if ((!strcmp(cur, "--n-blocks")) && (argi + 1 < argc)) {
      const int n_blocks = atoi(argv[++argi]);
      if (n_blocks < 2) {
        fprintf(stderr, "Error: --n-blocks must be at least 2.\n");
        return 1;
      }
      opts.n_blocks = n_blocks;
    } else if ((!strcmp(cur, "--chisq-max")) && (argi + 1 < argc)) {
      if ((!ScanadvDouble(argv[++argi], &opts.chisq_max)) || (opts.chisq_max <= 0.0)) {
        fprintf(stderr, "Error: --chisq-max must be positive.\n");
        return 1;
      }
      opts.have_chisq_max = 1;
    } else if ((!strcmp(cur, "--two-step")) && (argi + 1 < argc)) {
      if ((!ScanadvDouble(argv[++argi], &opts.two_step)) || (opts.two_step <= 0.0)) {
        fprintf(stderr, "Error: --two-step must be positive.\n");
        return 1;
      }
      opts.have_two_step = 1;
    } else if (!strcmp(cur, "--no-intercept")) {
      opts.no_intercept = 1;
    } else if (!strcmp(cur, "--no-check-alleles")) {
      opts.no_check_alleles = 1;
    } else if ((!strcmp(cur, "--intercept-h2")) && (argi + 1 < argc)) {
      if (LdscParseNumList(argv[++argi], "--intercept-h2", &opts.intercept_h2)) {
        return 1;
      }
    } else if ((!strcmp(cur, "--intercept-gencov")) && (argi + 1 < argc)) {
      if (LdscParseNumList(argv[++argi], "--intercept-gencov", &opts.intercept_gencov)) {
        return 1;
      }
    } else if ((!strcmp(cur, "--samp-prev")) && (argi + 1 < argc)) {
      if (LdscParseNumList(argv[++argi], "--samp-prev", &opts.samp_prev)) {
        return 1;
      }
    } else if ((!strcmp(cur, "--pop-prev")) && (argi + 1 < argc)) {
      if (LdscParseNumList(argv[++argi], "--pop-prev", &opts.pop_prev)) {
        return 1;
      }
    } else if (!strcmp(cur, "--version")) {
      printf("%s\n", kLdscVersion);
      return 0;
    } else {
      fprintf(stderr, "Error: unrecognized argument '%s'.\n", cur);
      return 1;
    }
  }
  if ((!(h2_fname || rg_arg)) || (!ref_ld_arg) || (!w_ld_arg) || (!out_prefix)) {
    fprintf(stderr,
            "%s\n"
            "LD Score regression (Bulik-Sullivan et al. 2015).\n\n"
            "Usage: ldsc --h2 <sumstats> --ref-ld[-chr] <LD Scores>\n"
            "            --w-ld[-chr] <LD Scores> --out <prefix>\n"
            "       ldsc --rg <sumstats1,sumstats2,...> --ref-ld[-chr] <LD Scores>\n"
            "            --w-ld[-chr] <LD Scores> --out <prefix>\n\n"
            "  --h2       Estimate SNP-heritability from one summary statistic\n"
            "             file.  Columns are detected case-insensitively: an ID\n"
            "             (SNP/ID), Z, and N.\n"
            "  --rg       Estimate genetic correlation between the first file and\n"
            "             each of the others.  A1 and A2 are then required too;\n"
            "             Z is flipped where the effect alleles are swapped, and\n"
            "             strand-ambiguous variants are dropped.\n"
            "  --ref-ld   LD Scores for the regression, with SNP (or ID) and L2\n"
            "             columns.  plink2 --ld-score output can be used directly,\n"
            "             and so can ldsc's <prefix>.l2.ldscore[.gz].\n"
            "  --ref-ld-chr  One LD Score fileset per chromosome; '@' in the\n"
            "             argument is replaced by the chromosome number, and\n"
            "             otherwise it is appended.\n"
            "  --w-ld[-chr]  LD Scores for the regression weights: the same sum of\n"
            "             r^2, but taken over only the regression variants.\n"
            "  --M        Number of variants the LD Scores were computed from.\n"
            "             Defaults to the .l2.M_5_50 files next to --ref-ld, or,\n"
            "             when those are absent, to the number of LD Scores read.\n"
            "  --not-M-5-50  Use the .l2.M files (all variants) rather than the\n"
            "             .l2.M_5_50 files (common variants only).\n"
            "  --n-blocks Jackknife block count (default 200).\n"
            "  --two-step Cutoff for the two-step estimator; the first step uses\n"
            "             the variants below it.  Defaults to 30 when the\n"
            "             intercept is free, and to off when it is constrained.\n"
            "  --chisq-max  Drop variants with chi^2 above this.  Off by default.\n"
            "  --no-intercept  Constrain the h2 intercept to 1 and the genetic\n"
            "             covariance intercept to 0, i.e. assume no confounding\n"
            "             and no sample overlap.\n"
            "  --intercept-h2  Constrain the h2 intercept(s) instead, one value\n"
            "             per phenotype ('N' to leave one free).\n"
            "  --intercept-gencov  Constrain the genetic covariance intercept(s).\n"
            "  --samp-prev, --pop-prev  Sample and population prevalence per\n"
            "             phenotype; both together convert the observed-scale\n"
            "             estimates to the liability scale.\n"
            "  --no-check-alleles  Skip the allele-matching step of --rg.  Only\n"
            "             safe if both files are known to be on the same strand\n"
            "             with the same effect alleles.\n"
            "  --out      Output prefix; writes <prefix>.log plus <prefix>.h2 or\n"
            "             <prefix>.rg with the estimates at full precision.\n",
            kLdscVersion);
    return 1;
  }
  if (h2_fname && rg_arg) {
    fprintf(stderr, "Error: --h2 and --rg cannot be used together.\n");
    return 1;
  }
  std::string log_path = std::string(out_prefix) + ".log";
  g_ldsc_logfile = fopen(log_path.c_str(), FOPEN_WB);
  if (!g_ldsc_logfile) {
    fprintf(stderr, "Error: Failed to open %s.\n", log_path.c_str());
    return 1;
  }
  LdscLog("%s\n", kLdscVersion);

  std::vector<std::string> ref_ids;
  std::vector<double> ref_l2;
  if (LdscReadLdscores(ref_ld_arg, ref_ld_chr_split, &ref_ids, &ref_l2)) {
    return 1;
  }
  LdscLog("Read reference panel LD Scores for %" PRIuPTR " variants.\n", S_CAST(uintptr_t, ref_ids.size()));
  if (ref_ids.empty()) {
    fprintf(stderr, "Error: No LD Scores read.\n");
    return 1;
  }
  std::vector<std::string> w_ids;
  std::vector<double> w_l2;
  if (LdscReadLdscores(w_ld_arg, w_ld_chr_split, &w_ids, &w_l2)) {
    return 1;
  }
  LdscLog("Read regression weight LD Scores for %" PRIuPTR " variants.\n", S_CAST(uintptr_t, w_ids.size()));
  std::unordered_map<std::string, double> w_ld_map;
  for (uintptr_t i = 0; i != w_ids.size(); ++i) {
    w_ld_map.emplace(w_ids[i], w_l2[i]);
  }

  double m_tot;
  if (opts.m_override > 0.0) {
    m_tot = opts.m_override;
  } else if (!LdscReadM(ref_ld_arg, ref_ld_chr_split, opts.not_m_5_50, &m_tot)) {
    LdscLog("Read M = %g from the %s files.\n", m_tot, opts.not_m_5_50? ".l2.M" : ".l2.M_5_50");
  } else {
    m_tot = u31tod(ref_ids.size());
    LdscLog("No %s file found; taking M to be the %g LD Scores read.  Pass --M if the\nLD Scores were computed from a different variant set.\n", opts.not_m_5_50? ".l2.M" : ".l2.M_5_50", m_tot);
  }

  BoolErr ret;
  if (h2_fname) {
    ret = LdscEstimateH2(h2_fname, ref_ids, ref_l2, w_ld_map, m_tot, &opts, out_prefix);
  } else {
    std::vector<std::string> rg_paths;
    LdscSplitComma(rg_arg, &rg_paths);
    if (rg_paths.size() < 2) {
      fprintf(stderr, "Error: --rg needs at least two summary statistic files.\n");
      return 1;
    }
    ret = LdscEstimateRg(rg_paths, ref_ids, ref_l2, w_ld_map, m_tot, &opts, out_prefix);
  }
  if (fclose(g_ldsc_logfile)) {
    fprintf(stderr, "Error: Failed to write %s.\n", log_path.c_str());
    return 1;
  }
  g_ldsc_logfile = nullptr;
  return ret? 1 : 0;
}
