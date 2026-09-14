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
#include "../../include/plink2_stats.h"
#include "../../include/plink2_string.h"
#include "../../include/plink2_text.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static const char kLdscVersion[] = "ldsc (plink-ng) v0.1";

// A partitioned regression has one coefficient per annotation, plus the
// intercept, so the parameter count is only bounded by the annotation count.
static const uint32_t kLdscMaxAnnot = 256;

static const uint32_t kLdscChrCt = 22;

static const double kLdscSqrt2Pi = 2.5066282746310002;

// ***** small-matrix and distribution helpers *****

// Solves mat * out = vec by Gauss-Jordan with partial pivoting.  dim is the
// parameter count, a small number even for the partitioned regression, so a
// dedicated small-matrix routine beats calling out to LAPACK.
BoolErr LdscSolve(uint32_t dim, const double* mat, const double* vec, double* out) {
  std::vector<double> aug_buf(S_CAST(uintptr_t, dim) * (dim + 1));
  double* aug = &(aug_buf[0]);
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
  std::vector<double> est;
  std::vector<double> jknife_est;
  std::vector<double> jknife_se;
  std::vector<double> jknife_cov;      // p x p, row-major
  std::vector<double> delete_values;   // n_blocks x p, row-major
} LdscJknife;

// Turns delete values and the whole-data estimate into the jackknife estimate
// and its covariance, via the pseudovalues.
void LdscFinishJknife(LdscJknife* jkp) {
  const uint32_t n_blocks = jkp->n_blocks;
  const uint32_t p = jkp->p;
  jkp->jknife_est.resize(p);
  jkp->jknife_se.resize(p);
  jkp->jknife_cov.assign(S_CAST(uintptr_t, p) * p, 0.0);
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
  std::vector<double> xtx_tot_buf(pp);
  std::vector<double> xty_tot_buf(p);
  double* xtx_tot = &(xtx_tot_buf[0]);
  double* xty_tot = &(xty_tot_buf[0]);
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
  jkp->est.resize(p);
  if (LdscSolve(p, xtx_tot, xty_tot, &(jkp->est[0]))) {
    return 1;
  }
  jkp->delete_values.assign(S_CAST(uintptr_t, n_blocks) * p, 0.0);
  std::vector<double> xtx_del_buf(pp);
  std::vector<double> xty_del_buf(p);
  double* xtx_del = &(xtx_del_buf[0]);
  double* xty_del = &(xty_del_buf[0]);
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

// Block jackknife for several ratios at once, giving their covariance.  The
// numerators vary by ratio; the denominator is shared.
void LdscRatioJknifeMulti(const double* est, const double* numer_delete, const double* denom_delete, uint32_t n_blocks, uint32_t p, std::vector<double>* cov, std::vector<double>* se) {
  LdscJknife jk;
  jk.n_blocks = n_blocks;
  jk.p = p;
  jk.est.assign(est, &(est[p]));
  jk.delete_values.resize(S_CAST(uintptr_t, n_blocks) * p);
  for (uint32_t b = 0; b != n_blocks; ++b) {
    for (uint32_t j = 0; j != p; ++j) {
      jk.delete_values[S_CAST(uintptr_t, b) * p + j] = numer_delete[S_CAST(uintptr_t, b) * p + j] / denom_delete[b];
    }
  }
  LdscFinishJknife(&jk);
  cov->swap(jk.jknife_cov);
  se->swap(jk.jknife_se);
}

// Block jackknife for a ratio of two jackknifed quantities.
void LdscRatioJknife(double est, const double* numer_delete, const double* denom_delete, uint32_t n_blocks, double* jknife_est_ptr, double* jknife_se_ptr) {
  LdscJknife jk;
  jk.n_blocks = n_blocks;
  jk.p = 1;
  jk.est.assign(1, est);
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
  const double* ld_mat;  // n_snp x n_annot, row-major
  const double* ld_tot;  // row sums of ld_mat; what the weights use
  const double* w_ld;
  const double* n_vec;
  const double* n1;
  const double* n2;
  const uint32_t* row_idxs;
  uint32_t n_annot;
  double m_tot;
  double nbar;
  double hsq1;
  double hsq2;
  double intercept_hsq1;
  double intercept_hsq2;
  uint32_t constrain_intercept;
  double fixed_intercept;
} LdscRegCtx;

// The weights are a function of the total heritability (or genetic
// covariance) and the total LD Score, not of the per-annotation split, so
// they stay scalar even in the partitioned regression.
void LdscUpdateWeights(const LdscRegCtx* ctx, const double* coef, uint32_t p, uint32_t n, double* w_out) {
  const double m_tot = ctx->m_tot;
  double param = 0.0;
  for (uint32_t j = 0; j != ctx->n_annot; ++j) {
    param += coef[j];
  }
  param = m_tot * param / ctx->nbar;
  double intercept;
  if (ctx->constrain_intercept) {
    intercept = ctx->fixed_intercept;
  } else {
    assert(p == ctx->n_annot + 1);
    intercept = coef[p - 1];
  }
  if (ctx->kind == kLdscHsq) {
    const double hsq = LdscClamp(param, 0.0, 1.0);
    for (uint32_t i = 0; i != n; ++i) {
      const uint32_t uidx = ctx->row_idxs[i];
      const double ld = MAXV(ctx->ld_tot[uidx], 1.0);
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
    const double ld = MAXV(ctx->ld_tot[uidx], 1.0);
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
  std::vector<double> xtx_buf(S_CAST(uintptr_t, p) * p);
  std::vector<double> xty_buf(p);
  double* xtx = &(xtx_buf[0]);
  double* xty = &(xty_buf[0]);
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

// Weights the rows of a design matrix and its response, then runs the block
// jackknife on them.
BoolErr LdscWeightedJknife(const double* x, const double* y, const double* w, uint32_t n, uint32_t p, const std::vector<uint32_t>& sep, LdscJknife* jkp) {
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
  std::vector<double> coef(p);
  std::vector<double> new_w(n);
  for (uint32_t iter = 0; iter != 2; ++iter) {
    if (LdscWls(x, y, &(w[0]), n, p, &(coef[0]))) {
      return 1;
    }
    LdscUpdateWeights(ctx, &(coef[0]), p, n, &(new_w[0]));
    for (uint32_t i = 0; i != n; ++i) {
      if (!(new_w[i] > 0.0)) {
        return 1;
      }
      w[i] = sqrt(new_w[i]);
    }
  }
  return LdscWeightedJknife(x, y, &(w[0]), n, p, sep, jkp);
}

// ***** the regressions *****

typedef struct LdscHsqResultStruct {
  double tot;
  double tot_se;
  double intercept;
  double intercept_se;
  double mean_chisq;
  double lambda_gc;
  double ratio;
  double ratio_se;
  uint32_t constrain_intercept;
  uint32_t ratio_valid;
  uint32_t n_blocks;
  uint32_t n_annot;
  std::vector<double> tot_delete_values;
  // Partitioned output, one entry per annotation.
  std::vector<double> coefs;
  std::vector<double> coef_ses;
  std::vector<double> cat;
  std::vector<double> cat_ses;
  std::vector<double> prop;
  std::vector<double> prop_ses;
  std::vector<double> enrichment;
  std::vector<double> m_prop;
  std::vector<double> coef_cov;  // n_annot x n_annot
  std::vector<double> prop_cov;  // n_annot x n_annot
} LdscHsqResult;

// Combines the free-intercept first step and the constrained-intercept second
// step of the two-step estimator into one jackknife, as
// LD_Score_Regression._combine_twostep_jknives().
void LdscCombineTwostep(const LdscJknife* step1, const LdscJknife* step2, double step1_int, double c, LdscJknife* out) {
  const uint32_t n_blocks = step1->n_blocks;
  out->n_blocks = n_blocks;
  out->p = 2;
  out->est.resize(2);
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
//
// With more than one annotation the weights are computed once from the
// aggregate estimate and left there, rather than being iterated: that is what
// the reference implementation does (its old_weights path), since the
// per-annotation coefficients do not pin down a single conditional variance
// to reweight by.
BoolErr LdscRegress(const LdscRegCtx* base_ctx, const double* y, uint32_t n_snp, uint32_t n_blocks, const double* fixed_intercept, const double* twostep, const double* step1_filter, double null_intercept, LdscJknife* jkp, double* nbar_ptr) {
  const double m_tot = base_ctx->m_tot;
  const uint32_t n_annot = base_ctx->n_annot;
  const double* ld_mat = base_ctx->ld_mat;
  const double* ld_tot = base_ctx->ld_tot;
  const double* n_vec = base_ctx->n_vec;
  double nbar = 0.0;
  double y_sum = 0.0;
  double ldn_sum = 0.0;
  for (uint32_t i = 0; i != n_snp; ++i) {
    nbar += n_vec[i];
    y_sum += y[i];
    ldn_sum += ld_tot[i] * n_vec[i];
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
  agg_ctx.n_annot = 1;
  // aggregate estimate -> initial weights.  The weight update takes a
  // coefficient on the N-scaled scale, so invert that scaling here.
  double agg_coef = tot_agg * nbar / m_tot;
  std::vector<double> initial_w(n_snp);
  LdscUpdateWeights(&agg_ctx, &agg_coef, 1, n_snp, &(initial_w[0]));

  // x is N-scaled to keep the condition number low.
  std::vector<double> x_scaled(S_CAST(uintptr_t, n_snp) * n_annot);
  for (uint32_t i = 0; i != n_snp; ++i) {
    for (uint32_t j = 0; j != n_annot; ++j) {
      x_scaled[S_CAST(uintptr_t, i) * n_annot + j] = n_vec[i] * ld_mat[S_CAST(uintptr_t, i) * n_annot + j] / nbar;
    }
  }

  if (fixed_intercept) {
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
    if (n_annot > 1) {
      std::vector<double> w(n_snp);
      for (uint32_t i = 0; i != n_snp; ++i) {
        if (!(initial_w[i] > 0.0)) {
          return 1;
        }
        w[i] = sqrt(initial_w[i]);
      }
      return LdscWeightedJknife(&(x_scaled[0]), &(yp[0]), &(w[0]), n_snp, n_annot, sep, jkp);
    }
    return LdscIrwls(&ctx, &(x_scaled[0]), &(yp[0]), &(initial_w[0]), n_snp, 1, sep, jkp);
  }

  const uint32_t p = n_annot + 1;
  std::vector<double> design(S_CAST(uintptr_t, n_snp) * p);
  for (uint32_t i = 0; i != n_snp; ++i) {
    for (uint32_t j = 0; j != n_annot; ++j) {
      design[S_CAST(uintptr_t, i) * p + j] = x_scaled[S_CAST(uintptr_t, i) * n_annot + j];
    }
    design[S_CAST(uintptr_t, i) * p + n_annot] = 1.0;
  }

  if ((!twostep) || (n_annot > 1)) {
    LdscRegCtx ctx = *base_ctx;
    ctx.nbar = nbar;
    ctx.row_idxs = &(all_idxs[0]);
    ctx.constrain_intercept = 0;
    std::vector<uint32_t> sep;
    LdscGetSeparators(n_snp, n_blocks, &sep);
    if (n_annot > 1) {
      std::vector<double> w(n_snp);
      for (uint32_t i = 0; i != n_snp; ++i) {
        if (!(initial_w[i] > 0.0)) {
          return 1;
        }
        w[i] = sqrt(initial_w[i]);
      }
      return LdscWeightedJknife(&(design[0]), y, &(w[0]), n_snp, p, sep, jkp);
    }
    return LdscIrwls(&ctx, &(design[0]), y, &(initial_w[0]), n_snp, p, sep, jkp);
  }

  // Two-step: a free-intercept regression on the low-signal variants, then a
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

BoolErr LdscHsqFit(const double* chisq, const double* ld_mat, const double* ld_tot, const double* w_ld, const double* n_vec, uint32_t n_snp, uint32_t n_annot, const double* m_vec, uint32_t n_blocks, const double* fixed_intercept, const double* twostep, LdscHsqResult* out) {
  double m_tot = 0.0;
  for (uint32_t j = 0; j != n_annot; ++j) {
    m_tot += m_vec[j];
  }
  LdscRegCtx ctx;
  ctx.kind = kLdscHsq;
  ctx.ld_mat = ld_mat;
  ctx.ld_tot = ld_tot;
  ctx.w_ld = w_ld;
  ctx.n_vec = n_vec;
  ctx.n1 = nullptr;
  ctx.n2 = nullptr;
  ctx.row_idxs = nullptr;
  ctx.n_annot = n_annot;
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
  out->n_annot = n_annot;
  out->n_blocks = jk.n_blocks;
  // Per-annotation coefficients, and the per-annotation heritability they
  // imply: cat_j = M_j * coef_j.
  out->coefs.resize(n_annot);
  out->coef_ses.resize(n_annot);
  out->cat.resize(n_annot);
  out->cat_ses.resize(n_annot);
  const uint32_t p = jk.p;
  for (uint32_t j = 0; j != n_annot; ++j) {
    out->coefs[j] = jk.est[j] / nbar;
    const double coef_var = jk.jknife_cov[j * p + j] / (nbar * nbar);
    out->coef_ses[j] = sqrt(coef_var);
    out->cat[j] = m_vec[j] * out->coefs[j];
    out->cat_ses[j] = sqrt(m_vec[j] * m_vec[j] * coef_var);
  }
  double tot = 0.0;
  double tot_cov = 0.0;
  for (uint32_t j = 0; j != n_annot; ++j) {
    tot += out->cat[j];
    for (uint32_t k = 0; k != n_annot; ++k) {
      tot_cov += m_vec[j] * m_vec[k] * jk.jknife_cov[j * p + k] / (nbar * nbar);
    }
  }
  out->tot = tot;
  out->tot_se = sqrt(tot_cov);
  if (fixed_intercept) {
    out->intercept = *fixed_intercept;
    out->intercept_se = 0.0 / 0.0;
  } else {
    out->intercept = jk.est[n_annot];
    out->intercept_se = jk.jknife_se[n_annot];
  }
  // Delete values for the total, which is what the genetic-correlation ratio
  // jackknife needs.
  out->tot_delete_values.assign(jk.n_blocks, 0.0);
  for (uint32_t b = 0; b != jk.n_blocks; ++b) {
    double acc = 0.0;
    for (uint32_t j = 0; j != n_annot; ++j) {
      acc += jk.delete_values[S_CAST(uintptr_t, b) * p + j] * m_vec[j];
    }
    out->tot_delete_values[b] = acc / nbar;
  }
  // Proportion of heritability per annotation, jackknifed as a ratio, plus
  // the enrichment that proportion implies.
  out->prop.assign(n_annot, 0.0 / 0.0);
  out->prop_ses.assign(n_annot, 0.0 / 0.0);
  out->enrichment.assign(n_annot, 0.0 / 0.0);
  out->m_prop.resize(n_annot);
  for (uint32_t j = 0; j != n_annot; ++j) {
    out->m_prop[j] = m_vec[j] / m_tot;
  }
  out->coef_cov.assign(S_CAST(uintptr_t, n_annot) * n_annot, 0.0);
  for (uint32_t j = 0; j != n_annot; ++j) {
    for (uint32_t k = 0; k != n_annot; ++k) {
      out->coef_cov[j * n_annot + k] = jk.jknife_cov[j * p + k] / (nbar * nbar);
    }
  }
  if (n_annot > 1) {
    std::vector<double> numer(S_CAST(uintptr_t, jk.n_blocks) * n_annot);
    std::vector<double> denom(jk.n_blocks);
    for (uint32_t b = 0; b != jk.n_blocks; ++b) {
      denom[b] = out->tot_delete_values[b];
      for (uint32_t j = 0; j != n_annot; ++j) {
        numer[S_CAST(uintptr_t, b) * n_annot + j] = m_vec[j] * jk.delete_values[S_CAST(uintptr_t, b) * p + j] / nbar;
      }
    }
    for (uint32_t j = 0; j != n_annot; ++j) {
      out->prop[j] = out->cat[j] / tot;
      out->enrichment[j] = (out->cat[j] / m_vec[j]) / (tot / m_tot);
    }
    std::vector<double> prop_se;
    LdscRatioJknifeMulti(&(out->prop[0]), &(numer[0]), &(denom[0]), jk.n_blocks, n_annot, &(out->prop_cov), &prop_se);
    out->prop_ses = prop_se;
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

BoolErr LdscGencovFit(const double* z1, const double* z2, const double* ld_mat, const double* ld_tot, const double* w_ld, const double* n1, const double* n2, uint32_t n_snp, uint32_t n_annot, const double* m_vec, uint32_t n_blocks, double hsq1, double hsq2, double intercept_hsq1, double intercept_hsq2, const double* fixed_intercept, const double* twostep, LdscGencovResult* out) {
  double m_tot = 0.0;
  for (uint32_t j = 0; j != n_annot; ++j) {
    m_tot += m_vec[j];
  }
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
  ctx.ld_mat = ld_mat;
  ctx.ld_tot = ld_tot;
  ctx.w_ld = w_ld;
  ctx.n_vec = &(sqrt_n1n2[0]);
  ctx.n1 = n1;
  ctx.n2 = n2;
  ctx.row_idxs = nullptr;
  ctx.n_annot = n_annot;
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
  const uint32_t p = jk.p;
  double tot = 0.0;
  double tot_cov = 0.0;
  for (uint32_t j = 0; j != n_annot; ++j) {
    tot += m_vec[j] * jk.est[j] / nbar;
    for (uint32_t k = 0; k != n_annot; ++k) {
      tot_cov += m_vec[j] * m_vec[k] * jk.jknife_cov[j * p + k] / (nbar * nbar);
    }
  }
  out->tot = tot;
  out->tot_se = sqrt(tot_cov);
  if (fixed_intercept) {
    out->intercept = *fixed_intercept;
    out->intercept_se = 0.0 / 0.0;
  } else {
    out->intercept = jk.est[n_annot];
    out->intercept_se = jk.jknife_se[n_annot];
  }
  out->n_blocks = jk.n_blocks;
  out->tot_delete_values.assign(jk.n_blocks, 0.0);
  for (uint32_t b = 0; b != jk.n_blocks; ++b) {
    double acc_del = 0.0;
    for (uint32_t j = 0; j != n_annot; ++j) {
      acc_del += jk.delete_values[S_CAST(uintptr_t, b) * p + j] * m_vec[j];
    }
    out->tot_delete_values[b] = acc_del / nbar;
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
// Defined with the other argument handling, below.
void LdscSplitComma(const char* arg, std::vector<std::string>* dst);

// One or more LD Score columns, plus the variant IDs they belong to.
typedef struct LdscScoresStruct {
  std::vector<std::string> ids;
  std::vector<std::string> annot_names;
  std::vector<uint32_t> per_file_annot_ct;  // one entry per --ref-ld fileset
  std::vector<double> l2;  // ids.size() x annot_names.size(), row-major
} LdscScores;

// Rows of one LD Score fileset, before sorting.
typedef struct LdscRawScoresStruct {
  std::vector<std::string> ids;
  std::vector<int64_t> chrom;
  std::vector<int64_t> bp;
  std::vector<double> l2;
  std::vector<std::string> annot_names;
  uint32_t has_pos;
} LdscRawScores;

// Reads one LD Score file, appending its rows.  Every column other than the
// variant ID, the position columns and MAF/CM is an annotation.
BoolErr LdscReadLdscoreFile(const char* fname, LdscRawScores* dst) {
  static const char* kIdNames[] = {"SNP", "ID", "RSID", nullptr};
  static const char* kChromNames[] = {"CHR", "CHROM", nullptr};
  static const char* kPosNames[] = {"BP", "POS", nullptr};
  static const char* kSkipNames[] = {"MAF", "CM", nullptr};
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
  uint32_t col_chrom = UINT32_MAX;
  uint32_t col_pos = UINT32_MAX;
  std::vector<uint32_t> annot_cols;
  std::vector<std::string> annot_names;
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
      } else if ((col_chrom == UINT32_MAX) && LdscMatchCol(iter, slen, kChromNames)) {
        col_chrom = col_ct;
      } else if ((col_pos == UINT32_MAX) && LdscMatchCol(iter, slen, kPosNames)) {
        col_pos = col_ct;
      } else if (!LdscMatchCol(iter, slen, kSkipNames)) {
        annot_cols.push_back(col_ct);
        annot_names.push_back(std::string(iter, slen));
      }
      iter = FirstNonTspace(token_end);
    }
  }
  if (col_id == UINT32_MAX) {
    fprintf(stderr, "Error: %s must have a SNP (or ID) column.\n", fname);
    return 1;
  }
  if (annot_cols.empty()) {
    fprintf(stderr, "Error: %s has no LD Score column.\n", fname);
    return 1;
  }
  const uint32_t n_annot = annot_cols.size();
  if (dst->annot_names.empty()) {
    dst->annot_names = annot_names;
    dst->has_pos = (col_chrom != UINT32_MAX) && (col_pos != UINT32_MAX);
  } else if (dst->annot_names != annot_names) {
    fprintf(stderr, "Error: %s does not have the same LD Score columns as the other files in\nits fileset.\n", fname);
    return 1;
  }
  uint32_t max_col = col_id;
  if (col_chrom != UINT32_MAX) {
    max_col = MAXV(max_col, col_chrom);
  }
  if (col_pos != UINT32_MAX) {
    max_col = MAXV(max_col, col_pos);
  }
  for (uint32_t i = 0; i != n_annot; ++i) {
    max_col = MAXV(max_col, annot_cols[i]);
  }
  std::vector<double> cur_vals(n_annot);
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
    int64_t chrom = 0;
    int64_t bp = 0;
    uint32_t annot_idx = 0;
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
      } else if (col_idx == col_chrom) {
        double cur;
        // Nonnumeric chromosome codes sort after the numeric ones, in the
        // order they are seen, which is the order the file is in.
        chrom = ScanadvDouble(iter, &cur)? S_CAST(int64_t, cur) : INT64_MAX;
      } else if (col_idx == col_pos) {
        double cur;
        if (!ScanadvDouble(iter, &cur)) {
          ok = 0;
          break;
        }
        bp = S_CAST(int64_t, cur);
      }
      if ((annot_idx != n_annot) && (col_idx == annot_cols[annot_idx])) {
        if (!ScanadvDouble(iter, &(cur_vals[annot_idx]))) {
          ok = 0;
          break;
        }
        ++annot_idx;
      }
      iter = FirstNonTspace(token_end);
    }
    if ((!ok) || (annot_idx != n_annot)) {
      // NA LD Scores (monomorphic variants, in plink2's output) are dropped.
      continue;
    }
    dst->ids.push_back(std::string(id_start, id_slen));
    dst->chrom.push_back(chrom);
    dst->bp.push_back(bp);
    for (uint32_t i = 0; i != n_annot; ++i) {
      dst->l2.push_back(cur_vals[i]);
    }
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
// or without a compression suffix) and a literal path, so plink2 --ld-score
// output can be handed over directly.
BoolErr LdscFindLdscorePath(const std::string& base, uint32_t allow_literal, std::string* out) {
  static const char* kSuffixes[] = {".l2.ldscore", ".l2.ldscore.gz", ".l2.ldscore.zst", ".ldscore", ".ldscore.gz", ".ldscore.zst", nullptr};
  for (uint32_t i = 0; kSuffixes[i]; ++i) {
    const std::string cand = base + kSuffixes[i];
    if (LdscFileExists(cand)) {
      *out = cand;
      return 0;
    }
  }
  if (allow_literal && LdscFileExists(base)) {
    *out = base;
    return 0;
  }
  return 1;
}

// Sorts one fileset's rows by position and drops repeated variant IDs, as the
// reference implementation does ("SEs will be wrong unless sorted": the
// jackknife blocks are cut on this order).
void LdscSortAndDedup(LdscRawScores* raw) {
  const uintptr_t row_ct = raw->ids.size();
  const uint32_t n_annot = raw->annot_names.size();
  std::vector<uint32_t> order(row_ct);
  for (uintptr_t i = 0; i != row_ct; ++i) {
    order[i] = i;
  }
  if (raw->has_pos) {
    const std::vector<int64_t>& chrom = raw->chrom;
    const std::vector<int64_t>& bp = raw->bp;
    std::stable_sort(order.begin(), order.end(), [&chrom, &bp](uint32_t a, uint32_t b) {
      if (chrom[a] != chrom[b]) {
        return chrom[a] < chrom[b];
      }
      return bp[a] < bp[b];
    });
  }
  std::vector<std::string> new_ids;
  std::vector<double> new_l2;
  std::unordered_map<std::string, uint32_t> seen;
  new_ids.reserve(row_ct);
  new_l2.reserve(row_ct * n_annot);
  for (uintptr_t i = 0; i != row_ct; ++i) {
    const uint32_t uidx = order[i];
    if (!seen.emplace(raw->ids[uidx], 1).second) {
      continue;
    }
    new_ids.push_back(raw->ids[uidx]);
    for (uint32_t j = 0; j != n_annot; ++j) {
      new_l2.push_back(raw->l2[S_CAST(uintptr_t, uidx) * n_annot + j]);
    }
  }
  raw->ids.swap(new_ids);
  raw->l2.swap(new_l2);
  raw->chrom.clear();
  raw->bp.clear();
}

// Reads one LD Score fileset: a single file, or one per chromosome
// concatenated.
BoolErr LdscReadOneFileset(const std::string& base, uint32_t is_chr_split, LdscRawScores* raw) {
  if (!is_chr_split) {
    std::string path;
    if (LdscFindLdscorePath(base, 1, &path)) {
      fprintf(stderr, "Error: Could not find LD Scores at %s[.l2.ldscore/.gz/.zst].\n", base.c_str());
      return 1;
    }
    if (LdscReadLdscoreFile(path.c_str(), raw)) {
      return 1;
    }
  } else {
    uint32_t found_ct = 0;
    for (uint32_t chr_idx = 1; chr_idx <= kLdscChrCt; ++chr_idx) {
      std::string path;
      if (LdscFindLdscorePath(LdscSubChr(base.c_str(), chr_idx), 0, &path)) {
        continue;
      }
      if (LdscReadLdscoreFile(path.c_str(), raw)) {
        return 1;
      }
      ++found_ct;
    }
    if (!found_ct) {
      fprintf(stderr, "Error: No LD Score files found for %s (expected\n%s<chr>.l2.ldscore[.gz]).\n", base.c_str(), base.c_str());
      return 1;
    }
  }
  LdscSortAndDedup(raw);
  return 0;
}

// Reads the LD Scores named by --ref-ld/--w-ld (which may be a comma-
// separated list of filesets, concatenated sideways into one annotation set)
// or --ref-ld-chr/--w-ld-chr (one fileset per chromosome).
BoolErr LdscReadLdscores(const char* arg, uint32_t is_chr_split, LdscScores* dst) {
  std::vector<std::string> bases;
  LdscSplitComma(arg, &bases);
  for (uintptr_t file_idx = 0; file_idx != bases.size(); ++file_idx) {
    LdscRawScores raw;
    raw.has_pos = 0;
    if (LdscReadOneFileset(bases[file_idx], is_chr_split, &raw)) {
      return 1;
    }
    const uint32_t cur_annot_ct = raw.annot_names.size();
    dst->per_file_annot_ct.push_back(cur_annot_ct);
    if (!file_idx) {
      dst->ids.swap(raw.ids);
      dst->l2.swap(raw.l2);
      dst->annot_names = raw.annot_names;
      if (bases.size() > 1) {
        // The reference implementation suffixes the column names with the
        // fileset index, since separate filesets can reuse a name.
        for (uint32_t j = 0; j != cur_annot_ct; ++j) {
          dst->annot_names[j] += "_0";
        }
      }
      continue;
    }
    if (raw.ids != dst->ids) {
      fprintf(stderr, "Error: LD Score filesets for concatenation must cover the same variants in\nthe same order; %s does not match %s.\n", bases[file_idx].c_str(), bases[0].c_str());
      return 1;
    }
    const uint32_t prev_annot_ct = dst->annot_names.size();
    const uintptr_t row_ct = dst->ids.size();
    std::vector<double> merged(row_ct * (prev_annot_ct + cur_annot_ct));
    for (uintptr_t i = 0; i != row_ct; ++i) {
      for (uint32_t j = 0; j != prev_annot_ct; ++j) {
        merged[i * (prev_annot_ct + cur_annot_ct) + j] = dst->l2[i * prev_annot_ct + j];
      }
      for (uint32_t j = 0; j != cur_annot_ct; ++j) {
        merged[i * (prev_annot_ct + cur_annot_ct) + prev_annot_ct + j] = raw.l2[i * cur_annot_ct + j];
      }
    }
    dst->l2.swap(merged);
    char suffix[16];
    snprintf(suffix, sizeof(suffix), "_%" PRIuPTR, file_idx);
    for (uint32_t j = 0; j != cur_annot_ct; ++j) {
      dst->annot_names.push_back(raw.annot_names[j] + suffix);
    }
  }
  if (dst->annot_names.size() > kLdscMaxAnnot) {
    fprintf(stderr, "Error: %" PRIuPTR " LD Score columns exceeds the %u-annotation limit.\n", S_CAST(uintptr_t, dst->annot_names.size()), kLdscMaxAnnot);
    return 1;
  }
  return 0;
}

// Reads one .l{2}.M[_5_50] file: one value per annotation.
BoolErr LdscReadMFile(const std::string& path, std::vector<double>* dst) {
  FILE* f = fopen(path.c_str(), FOPEN_RB);
  if (!f) {
    return 1;
  }
  char buf[16384];
  if (!fgets(buf, sizeof(buf), f)) {
    fclose(f);
    return 1;
  }
  fclose(f);
  dst->clear();
  const char* iter = FirstNonTspace(buf);
  while (!IsEolnKns(*iter)) {
    double cur;
    if (!ScanadvDouble(iter, &cur)) {
      fprintf(stderr, "Error: Malformed %s.\n", path.c_str());
      return 1;
    }
    dst->push_back(cur);
    iter = FirstNonTspace(CurTokenEnd(iter));
  }
  return dst->empty();
}

// ***** annotation overlap *****

// Reads one .frq file into an ID-keyed map.  (The reference implementation
// aligns the .annot and .frq files by row position; matching on the variant
// ID instead gives the same answer whenever that assumption holds, and a
// correct one when it does not.)
BoolErr LdscReadFrq(const char* fname, std::unordered_map<std::string, double>* dst) {
  static const char* kIdNames[] = {"SNP", "ID", nullptr};
  static const char* kFrqNames[] = {"FRQ", "MAF", nullptr};
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
  uint32_t col_frq = UINT32_MAX;
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
      } else if ((col_frq == UINT32_MAX) && LdscMatchCol(iter, slen, kFrqNames)) {
        col_frq = col_ct;
      }
      iter = FirstNonTspace(token_end);
    }
  }
  if ((col_id == UINT32_MAX) || (col_frq == UINT32_MAX)) {
    fprintf(stderr, "Error: %s must have SNP and FRQ (or MAF) columns.\n", fname);
    return 1;
  }
  const uint32_t max_col = MAXV(col_id, col_frq);
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
    double frq = 0.0;
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
      } else if (col_idx == col_frq) {
        if (!ScanadvDouble(iter, &frq)) {
          ok = 0;
          break;
        }
      }
      iter = FirstNonTspace(token_end);
    }
    if (!ok) {
      continue;
    }
    dst->emplace(std::string(id_start, id_slen), frq);
  }
  reterr = kPglRetSuccess;
  CleanupTextStream(&txs, &reterr);
  if (reterr) {
    fprintf(stderr, "Error: Failed to read %s.\n", fname);
    return 1;
  }
  return 0;
}

// Reads one .annot file, accumulating A'A over its rows and counting them.
// Only the common variants are kept when a .frq map is supplied, matching the
// 5%-50% minor allele frequency band the .l2.M_5_50 counts use.
BoolErr LdscAccumAnnotFile(const char* fname, const std::unordered_map<std::string, double>* frq_map, uint32_t n_annot, uint32_t annot_offset, uint32_t total_annot, std::vector<double>* overlap, double* row_ct_ptr, std::vector<std::vector<double> >* rows_out) {
  static const char* kIdNames[] = {"SNP", "ID", nullptr};
  static const char* kSkipNames[] = {"CHR", "CHROM", "BP", "POS", "CM", nullptr};
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
  std::vector<uint32_t> annot_cols;
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
      } else if (!LdscMatchCol(iter, slen, kSkipNames)) {
        annot_cols.push_back(col_ct);
      }
      iter = FirstNonTspace(token_end);
    }
  }
  if (annot_cols.size() != n_annot) {
    fprintf(stderr, "Error: %s has %" PRIuPTR " annotation columns, but its LD Score fileset has\n%u.\n", fname, S_CAST(uintptr_t, annot_cols.size()), n_annot);
    return 1;
  }
  if (frq_map && (col_id == UINT32_MAX)) {
    fprintf(stderr, "Error: %s needs a SNP column to be matched against the .frq file.\n", fname);
    return 1;
  }
  uint32_t max_col = col_id;
  if (max_col == UINT32_MAX) {
    max_col = 0;
  }
  for (uint32_t j = 0; j != n_annot; ++j) {
    max_col = MAXV(max_col, annot_cols[j]);
  }
  std::vector<double> cur(n_annot);
  uintptr_t row_idx = 0;
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
    uint32_t annot_idx = 0;
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
      if ((annot_idx != n_annot) && (col_idx == annot_cols[annot_idx])) {
        if (!ScanadvDouble(iter, &(cur[annot_idx]))) {
          ok = 0;
          break;
        }
        ++annot_idx;
      }
      iter = FirstNonTspace(token_end);
    }
    if ((!ok) || (annot_idx != n_annot)) {
      fprintf(stderr, "Error: Malformed line in %s.\n", fname);
      return 1;
    }
    if (frq_map) {
      const std::unordered_map<std::string, double>::const_iterator frq_it = frq_map->find(std::string(id_start, id_slen));
      if (frq_it == frq_map->end()) {
        continue;
      }
      const double frq = frq_it->second;
      if ((frq <= 0.05) || (frq >= 0.95)) {
        continue;
      }
    }
    if (rows_out) {
      // First fileset of several: the rows have to be held so the later
      // filesets' columns can be paired with them.
      if (rows_out->size() <= row_idx) {
        rows_out->resize(row_idx + 1);
      }
      std::vector<double>& dst_row = (*rows_out)[row_idx];
      dst_row.resize(total_annot, 0.0);
      for (uint32_t j = 0; j != n_annot; ++j) {
        dst_row[annot_offset + j] = cur[j];
      }
    } else {
      for (uint32_t j = 0; j != n_annot; ++j) {
        for (uint32_t k = 0; k != n_annot; ++k) {
          (*overlap)[S_CAST(uintptr_t, annot_offset + j) * total_annot + annot_offset + k] += cur[j] * cur[k];
        }
      }
    }
    ++row_idx;
  }
  reterr = kPglRetSuccess;
  CleanupTextStream(&txs, &reterr);
  if (reterr) {
    fprintf(stderr, "Error: Failed to read %s.\n", fname);
    return 1;
  }
  *row_ct_ptr = u31tod(row_idx);
  return 0;
}

// Builds the annotation overlap matrix A'A and the variant count behind it,
// from the .annot files next to --ref-ld.
BoolErr LdscReadAnnot(const char* ref_arg, uint32_t is_chr_split, const char* frq_arg, uint32_t frq_is_chr_split, uint32_t total_annot, const std::vector<uint32_t>& per_file_annot_ct, std::vector<double>* overlap, double* m_tot_ptr) {
  std::vector<std::string> bases;
  LdscSplitComma(ref_arg, &bases);
  if (bases.size() != per_file_annot_ct.size()) {
    fprintf(stderr, "Error: internal error: LD Score fileset count mismatch.\n");
    return 1;
  }
  overlap->assign(S_CAST(uintptr_t, total_annot) * total_annot, 0.0);
  double m_tot = 0.0;
  const uint32_t chr_end = is_chr_split? (kLdscChrCt + 1) : 1;
  for (uint32_t chr_idx = 0; chr_idx != chr_end; ++chr_idx) {
    std::unordered_map<std::string, double> frq_map;
    if (frq_arg) {
      std::string frq_path;
      const std::string frq_base = frq_is_chr_split? LdscSubChr(frq_arg, chr_idx + 1) : std::string(frq_arg);
      static const char* kFrqSuffixes[] = {".frq", ".frq.gz", ".frq.zst", "", nullptr};
      uint32_t found = 0;
      for (uint32_t i = 0; kFrqSuffixes[i]; ++i) {
        const std::string cand = frq_base + kFrqSuffixes[i];
        if (LdscFileExists(cand)) {
          frq_path = cand;
          found = 1;
          break;
        }
      }
      if (!found) {
        if (is_chr_split) {
          continue;
        }
        fprintf(stderr, "Error: Could not find %s[.frq].\n", frq_base.c_str());
        return 1;
      }
      if (LdscReadFrq(frq_path.c_str(), &frq_map)) {
        return 1;
      }
    }
    // With several filesets the rows have to be paired up before A'A can be
    // accumulated, since an overlap entry can span two filesets.
    std::vector<std::vector<double> > rows;
    const uint32_t multi = (bases.size() > 1);
    uint32_t annot_offset = 0;
    uint32_t any_found = 0;
    double cur_row_ct = 0.0;
    for (uintptr_t file_idx = 0; file_idx != bases.size(); ++file_idx) {
      const std::string base = is_chr_split? LdscSubChr(bases[file_idx].c_str(), chr_idx + 1) : bases[file_idx];
      std::string path;
      static const char* kSuffixes[] = {".annot", ".annot.gz", ".annot.zst", nullptr};
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
        if (is_chr_split && (!file_idx)) {
          break;
        }
        fprintf(stderr, "Error: Could not find %s.annot[.gz].  --overlap-annot needs one .annot\nfile per LD Score fileset.\n", base.c_str());
        return 1;
      }
      any_found = 1;
      if (LdscAccumAnnotFile(path.c_str(), frq_arg? &frq_map : nullptr, per_file_annot_ct[file_idx], annot_offset, total_annot, overlap, &cur_row_ct, multi? &rows : nullptr)) {
        return 1;
      }
      annot_offset += per_file_annot_ct[file_idx];
    }
    if (!any_found) {
      continue;
    }
    if (multi) {
      for (uintptr_t i = 0; i != rows.size(); ++i) {
        const std::vector<double>& row = rows[i];
        for (uint32_t j = 0; j != total_annot; ++j) {
          if (row[j] == 0.0) {
            continue;
          }
          for (uint32_t k = 0; k != total_annot; ++k) {
            (*overlap)[S_CAST(uintptr_t, j) * total_annot + k] += row[j] * row[k];
          }
        }
      }
      cur_row_ct = u31tod(rows.size());
    }
    m_tot += cur_row_ct;
  }
  if (m_tot == 0.0) {
    fprintf(stderr, "Error: No .annot rows read.\n");
    return 1;
  }
  *m_tot_ptr = m_tot;
  return 0;
}

// Reads the .l2.M_5_50 (or .l2.M) files holding the number of variants the LD
// Scores were computed from, one value per annotation, summed over
// chromosomes and concatenated over filesets.
BoolErr LdscReadM(const char* arg, uint32_t is_chr_split, uint32_t not_m_5_50, std::vector<double>* m_vec) {
  const char* suffix = not_m_5_50? ".l2.M" : ".l2.M_5_50";
  std::vector<std::string> bases;
  LdscSplitComma(arg, &bases);
  m_vec->clear();
  for (uintptr_t file_idx = 0; file_idx != bases.size(); ++file_idx) {
    std::vector<double> acc;
    uint32_t found_ct = 0;
    const uint32_t chr_end = is_chr_split? (kLdscChrCt + 1) : 1;
    for (uint32_t chr_idx = 0; chr_idx != chr_end; ++chr_idx) {
      std::string path;
      if (is_chr_split) {
        path = LdscSubChr(bases[file_idx].c_str(), chr_idx + 1) + suffix;
      } else {
        path = bases[file_idx] + suffix;
      }
      std::vector<double> cur;
      if (LdscReadMFile(path, &cur)) {
        continue;
      }
      if (!found_ct) {
        acc = cur;
      } else {
        if (cur.size() != acc.size()) {
          fprintf(stderr, "Error: %s has %" PRIuPTR " entries, but the other chromosomes have %" PRIuPTR ".\n", path.c_str(), S_CAST(uintptr_t, cur.size()), S_CAST(uintptr_t, acc.size()));
          return 1;
        }
        for (uintptr_t j = 0; j != acc.size(); ++j) {
          acc[j] += cur[j];
        }
      }
      ++found_ct;
    }
    if (!found_ct) {
      return 1;
    }
    for (uintptr_t j = 0; j != acc.size(); ++j) {
      m_vec->push_back(acc[j]);
    }
  }
  return m_vec->empty();
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
  std::vector<double> ld;      // n_snp x n_annot, row-major
  std::vector<double> ld_tot;  // row sums
  uint32_t n_annot;
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
void LdscMerge(const LdscScores& ref, const std::unordered_map<std::string, double>& w_ld_map, const std::unordered_map<std::string, LdscSumstatRow>& sumstats, uint32_t keep_alleles, LdscData* dst) {
  const uintptr_t ref_ct = ref.ids.size();
  const uint32_t n_annot = ref.annot_names.size();
  dst->n_annot = n_annot;
  for (uintptr_t i = 0; i != ref_ct; ++i) {
    const std::unordered_map<std::string, LdscSumstatRow>::const_iterator ss_it = sumstats.find(ref.ids[i]);
    if (ss_it == sumstats.end()) {
      continue;
    }
    const std::unordered_map<std::string, double>::const_iterator w_it = w_ld_map.find(ref.ids[i]);
    if (w_it == w_ld_map.end()) {
      continue;
    }
    dst->ids.push_back(ref.ids[i]);
    double acc = 0.0;
    for (uint32_t j = 0; j != n_annot; ++j) {
      const double cur = ref.l2[i * n_annot + j];
      dst->ld.push_back(cur);
      acc += cur;
    }
    dst->ld_tot.push_back(acc);
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
  if (hsq->n_annot > 1) {
    LdscLog("Categories: %u (see the .results file for the per-category estimates)\n", hsq->n_annot);
  }
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

// ***** munge: raw GWAS summary statistics to .sumstats *****
//
// The reference implementation's munge_sumstats.py: detect the columns, throw
// out what LD Score regression should not see (low imputation quality, rare
// variants, out-of-range p-values, strand-ambiguous or non-SNP alleles,
// duplicated IDs, low sample size), turn the p-value and the signed statistic
// into a signed Z, and write SNP/A1/A2/Z/N.

typedef enum {
  kMungeFieldNone = 0,
  kMungeFieldSnp,
  kMungeFieldP,
  kMungeFieldA1,
  kMungeFieldA2,
  kMungeFieldN,
  kMungeFieldNCas,
  kMungeFieldNCon,
  kMungeFieldNStudy,
  kMungeFieldInfo,
  kMungeFieldFrq,
  kMungeFieldZ,
  kMungeFieldOr,
  kMungeFieldBeta,
  kMungeFieldLogOdds,
  kMungeFieldSigned
} MungeField;

typedef struct MungeAliasStruct {
  const char* name;
  MungeField field;
} MungeAlias;

// The reference implementation's default_cnames, with its aliases.
static const MungeAlias kMungeAliases[] = {
  {"SNP", kMungeFieldSnp}, {"MARKERNAME", kMungeFieldSnp}, {"SNPID", kMungeFieldSnp},
  {"RS", kMungeFieldSnp}, {"RSID", kMungeFieldSnp}, {"RS_NUMBER", kMungeFieldSnp},
  {"RS_NUMBERS", kMungeFieldSnp},
  {"NSTUDY", kMungeFieldNStudy}, {"N_STUDY", kMungeFieldNStudy},
  {"NSTUDIES", kMungeFieldNStudy}, {"N_STUDIES", kMungeFieldNStudy},
  {"P", kMungeFieldP}, {"PVALUE", kMungeFieldP}, {"P_VALUE", kMungeFieldP},
  {"PVAL", kMungeFieldP}, {"P_VAL", kMungeFieldP}, {"GC_PVALUE", kMungeFieldP},
  {"A1", kMungeFieldA1}, {"ALLELE1", kMungeFieldA1}, {"ALLELE_1", kMungeFieldA1},
  {"EFFECT_ALLELE", kMungeFieldA1}, {"REFERENCE_ALLELE", kMungeFieldA1},
  {"INC_ALLELE", kMungeFieldA1}, {"EA", kMungeFieldA1},
  {"A2", kMungeFieldA2}, {"ALLELE2", kMungeFieldA2}, {"ALLELE_2", kMungeFieldA2},
  {"OTHER_ALLELE", kMungeFieldA2}, {"NON_EFFECT_ALLELE", kMungeFieldA2},
  {"DEC_ALLELE", kMungeFieldA2}, {"NEA", kMungeFieldA2},
  {"N", kMungeFieldN}, {"WEIGHT", kMungeFieldN},
  {"NCASE", kMungeFieldNCas}, {"CASES_N", kMungeFieldNCas}, {"N_CASE", kMungeFieldNCas},
  {"N_CASES", kMungeFieldNCas}, {"N_CAS", kMungeFieldNCas},
  {"N_CONTROLS", kMungeFieldNCon}, {"N_CON", kMungeFieldNCon},
  {"NCONTROL", kMungeFieldNCon}, {"CONTROLS_N", kMungeFieldNCon},
  {"N_CONTROL", kMungeFieldNCon},
  {"ZSCORE", kMungeFieldZ}, {"Z_SCORE", kMungeFieldZ}, {"GC_ZSCORE", kMungeFieldZ},
  {"Z", kMungeFieldZ},
  {"OR", kMungeFieldOr},
  {"B", kMungeFieldBeta}, {"BETA", kMungeFieldBeta}, {"EFFECTS", kMungeFieldBeta},
  {"EFFECT", kMungeFieldBeta},
  {"LOG_ODDS", kMungeFieldLogOdds},
  {"INFO", kMungeFieldInfo},
  {"EAF", kMungeFieldFrq}, {"FRQ", kMungeFieldFrq}, {"MAF", kMungeFieldFrq},
  {"FRQ_U", kMungeFieldFrq}, {"F_U", kMungeFieldFrq},
  {nullptr, kMungeFieldNone}
};

// Uppercases, and maps '-' and '.' to '_', as the reference implementation's
// clean_header() does.
std::string LdscCleanHeader(const char* tok, uint32_t slen) {
  std::string out;
  out.reserve(slen);
  for (uint32_t i = 0; i != slen; ++i) {
    char c = tok[i];
    if ((c >= 'a') && (c <= 'z')) {
      c -= 32;
    } else if ((c == '-') || (c == '.')) {
      c = '_';
    }
    out.push_back(c);
  }
  return out;
}

// The value a signed statistic takes when there is no effect.
double LdscSignedNull(MungeField field) {
  return (field == kMungeFieldOr)? 1.0 : 0.0;
}

const char* LdscSignedName(MungeField field) {
  switch (field) {
  case kMungeFieldZ: return "Z";
  case kMungeFieldOr: return "OR";
  case kMungeFieldBeta: return "BETA";
  case kMungeFieldLogOdds: return "LOG_ODDS";
  default: return "SIGNED_SUMSTAT";
  }
}

// The reference implementation writes an integral sample size without a
// fractional part, and a scaled one to three decimals.
void LdscWriteN(double n, uint32_t integral_column, FILE* out_file) {
  if (integral_column) {
    fprintf(out_file, "%.0f", n);
  } else {
    fprintf(out_file, "%.3f", n);
  }
}

typedef struct MungeRowStruct {
  std::string id;
  char a1;
  char a2;
  double p;
  double signed_stat;
  double n;
  double n_cas;
  double n_con;
  double frq;
  double nstudy;
} MungeRow;

typedef struct MungeOptsStruct {
  const char* fname;
  const char* merge_alleles;
  const char* snp_col;
  const char* a1_col;
  const char* a2_col;
  const char* p_col;
  const char* n_col;
  const char* n_cas_col;
  const char* n_con_col;
  const char* frq_col;
  const char* info_col;
  const char* info_list;
  const char* nstudy_col;
  const char* signed_sumstats;  // "<column>,<null value>"
  const char* ignore;
  double n_override;
  double n_cas_override;
  double n_con_override;
  double info_min;
  double maf_min;
  double n_min;
  uint32_t have_n_min;
  double nstudy_min;
  uint32_t have_nstudy_min;
  uint32_t a1_inc;
  uint32_t no_alleles;
  uint32_t keep_maf;
  uint32_t daner;
} MungeOpts;

// An allele pair the reference implementation would keep: two different ACGT
// alleles that are not each other's complement.
uint32_t LdscMungeValidSnp(char a1, char a2) {
  return LdscIsValidSnp(a1, a2);
}

// The .sumstats-format allele list --merge-alleles restricts to.
typedef struct MungeMergeListStruct {
  std::vector<std::string> ids;
  std::vector<char> a1;
  std::vector<char> a2;
  std::unordered_map<std::string, uint32_t> id_to_idx;
} MungeMergeList;

BoolErr LdscReadMergeAlleles(const char* fname, MungeMergeList* dst) {
  static const char* kIdNames[] = {"SNP", "ID", nullptr};
  static const char* kA1Names[] = {"A1", nullptr};
  static const char* kA2Names[] = {"A2", nullptr};
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
      }
      iter = FirstNonTspace(token_end);
    }
  }
  if ((col_id == UINT32_MAX) || (col_a1 == UINT32_MAX) || (col_a2 == UINT32_MAX)) {
    fprintf(stderr, "Error: --merge-alleles file must have SNP, A1 and A2 columns.\n");
    return 1;
  }
  uint32_t max_col = MAXV(col_id, MAXV(col_a1, col_a2));
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
    char a1 = '\0';
    char a2 = '\0';
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
      } else if (col_idx == col_a1) {
        if (slen != 1) {
          ok = 0;
          break;
        }
        a1 = LdscUpcase(*iter);
      } else if (col_idx == col_a2) {
        if (slen != 1) {
          ok = 0;
          break;
        }
        a2 = LdscUpcase(*iter);
      }
      iter = FirstNonTspace(token_end);
    }
    if (!ok) {
      continue;
    }
    const std::string id(id_start, id_slen);
    if (!dst->id_to_idx.emplace(id, dst->ids.size()).second) {
      continue;
    }
    dst->ids.push_back(id);
    dst->a1.push_back(a1);
    dst->a2.push_back(a2);
  }
  reterr = kPglRetSuccess;
  CleanupTextStream(&txs, &reterr);
  if (reterr) {
    fprintf(stderr, "Error: Failed to read %s.\n", fname);
    return 1;
  }
  return dst->ids.empty();
}

BoolErr LdscMungeSumstats(const MungeOpts* mopts, const char* out_prefix) {
  MungeMergeList merge_list;
  if (mopts->merge_alleles) {
    if (mopts->no_alleles) {
      fprintf(stderr, "Error: --no-alleles and --merge-alleles cannot be used together.\n");
      return 1;
    }
    if (LdscReadMergeAlleles(mopts->merge_alleles, &merge_list)) {
      fprintf(stderr, "Error: No usable rows in %s.\n", mopts->merge_alleles);
      return 1;
    }
    LdscLog("Read %" PRIuPTR " variants for allele merge from %s.\n", S_CAST(uintptr_t, merge_list.ids.size()), mopts->merge_alleles);
  }

  // Column detection: the flag overrides first, then the alias table.
  std::unordered_map<std::string, MungeField> cname_map;
  for (uint32_t i = 0; kMungeAliases[i].name; ++i) {
    cname_map[kMungeAliases[i].name] = kMungeAliases[i].field;
  }
  std::vector<std::string> ignore_names;
  if (mopts->ignore) {
    std::vector<std::string> parts;
    LdscSplitComma(mopts->ignore, &parts);
    for (uintptr_t i = 0; i != parts.size(); ++i) {
      ignore_names.push_back(LdscCleanHeader(parts[i].c_str(), parts[i].size()));
    }
  }
  double signed_null = 0.0;
  uint32_t have_signed_flag = 0;
  std::string signed_flag_name;
  if (mopts->signed_sumstats) {
    const char* comma = strrchr(mopts->signed_sumstats, ',');
    if (!comma) {
      fprintf(stderr, "Error: --signed-sumstats takes <column name>,<null value>.\n");
      return 1;
    }
    if (!ScanadvDouble(&(comma[1]), &signed_null)) {
      fprintf(stderr, "Error: Invalid --signed-sumstats null value '%s'.\n", &(comma[1]));
      return 1;
    }
    signed_flag_name = LdscCleanHeader(mopts->signed_sumstats, comma - mopts->signed_sumstats);
    have_signed_flag = 1;
  }
  struct {
    const char* arg;
    MungeField field;
  } flag_cols[] = {
    {mopts->snp_col, kMungeFieldSnp},
    {mopts->a1_col, kMungeFieldA1},
    {mopts->a2_col, kMungeFieldA2},
    {mopts->p_col, kMungeFieldP},
    {mopts->n_col, kMungeFieldN},
    {mopts->n_cas_col, kMungeFieldNCas},
    {mopts->n_con_col, kMungeFieldNCon},
    {mopts->frq_col, kMungeFieldFrq},
    {mopts->info_col, kMungeFieldInfo},
    {mopts->nstudy_col, kMungeFieldNStudy}
  };
  std::unordered_map<std::string, MungeField> flag_map;
  for (uint32_t i = 0; i != sizeof(flag_cols) / sizeof(flag_cols[0]); ++i) {
    if (flag_cols[i].arg) {
      flag_map[LdscCleanHeader(flag_cols[i].arg, strlen(flag_cols[i].arg))] = flag_cols[i].field;
    }
  }
  if (mopts->info_list) {
    std::vector<std::string> parts;
    LdscSplitComma(mopts->info_list, &parts);
    for (uintptr_t i = 0; i != parts.size(); ++i) {
      flag_map[LdscCleanHeader(parts[i].c_str(), parts[i].size())] = kMungeFieldInfo;
    }
  }
  if (have_signed_flag) {
    flag_map[signed_flag_name] = kMungeFieldSigned;
  }

  TextStream txs;
  PreinitTextStream(&txs);
  PglErr reterr = TextStreamOpen(mopts->fname, &txs);
  if (reterr) {
    fprintf(stderr, "Error: Failed to open %s.\n", mopts->fname);
    return 1;
  }
  const char* header = TextGet(&txs);
  if (!header) {
    fprintf(stderr, "Error: %s is empty.\n", mopts->fname);
    return 1;
  }
  std::vector<MungeField> col_fields;
  std::vector<std::string> col_names;
  double daner_n_cas = 0.0;
  double daner_n_con = 0.0;
  {
    const char* iter = FirstNonTspace(header);
    if (*iter == '#') {
      ++iter;
    }
    for (; !IsEolnKns(*iter); ) {
      const char* token_end = CurTokenEnd(iter);
      const uint32_t slen = token_end - iter;
      const std::string cleaned = LdscCleanHeader(iter, slen);
      col_names.push_back(cleaned);
      MungeField field = kMungeFieldNone;
      const std::unordered_map<std::string, MungeField>::const_iterator flag_it = flag_map.find(cleaned);
      if (flag_it != flag_map.end()) {
        field = flag_it->second;
      } else {
        uint32_t ignored = 0;
        for (uintptr_t i = 0; i != ignore_names.size(); ++i) {
          if (ignore_names[i] == cleaned) {
            ignored = 1;
            break;
          }
        }
        if (!ignored) {
          const std::unordered_map<std::string, MungeField>::const_iterator it = cname_map.find(cleaned);
          if (it != cname_map.end()) {
            field = it->second;
          }
        }
      }
      if (mopts->daner) {
        // PGC daner format: the case and control counts are in the FRQ_A_/
        // FRQ_U_ column names, and FRQ_U_ is the frequency column.
        if (!cleaned.compare(0, 6, "FRQ_A_")) {
          double cur;
          if (ScanadvDouble(&(cleaned.c_str()[6]), &cur)) {
            daner_n_cas = cur;
          }
          field = kMungeFieldNone;
        } else if (!cleaned.compare(0, 6, "FRQ_U_")) {
          double cur;
          if (ScanadvDouble(&(cleaned.c_str()[6]), &cur)) {
            daner_n_con = cur;
          }
          field = kMungeFieldFrq;
        }
      }
      col_fields.push_back(field);
      iter = FirstNonTspace(token_end);
    }
  }
  const uint32_t col_ct = col_fields.size();
  if (mopts->daner) {
    if ((daner_n_cas <= 0.0) || (daner_n_con <= 0.0)) {
      fprintf(stderr, "Error: --daner needs FRQ_A_<case count> and FRQ_U_<control count>\ncolumns.\n");
      return 1;
    }
    LdscLog("Inferred N_cas = %g, N_con = %g from the FRQ_[A/U] column names.\n", daner_n_cas, daner_n_con);
  }

  // One signed statistic, and no field claimed twice.
  MungeField signed_field = kMungeFieldNone;
  std::string signed_col_name;
  {
    std::unordered_map<int, uint32_t> field_counts;
    for (uint32_t i = 0; i != col_ct; ++i) {
      if (col_fields[i] == kMungeFieldNone) {
        continue;
      }
      if (++field_counts[S_CAST(int, col_fields[i])] > 1) {
        fprintf(stderr, "Error: %s has two columns mapping to the same field (%s).  Use --ignore, or\nname the column you want with the matching flag.\n", mopts->fname, col_names[i].c_str());
        return 1;
      }
    }
    if (!mopts->a1_inc) {
      for (uint32_t i = 0; i != col_ct; ++i) {
        const MungeField field = col_fields[i];
        const uint32_t is_signed = (field == kMungeFieldSigned) || (field == kMungeFieldZ) || (field == kMungeFieldOr) || (field == kMungeFieldBeta) || (field == kMungeFieldLogOdds);
        if (!is_signed) {
          continue;
        }
        if (signed_field != kMungeFieldNone) {
          fprintf(stderr, "Error: %s has more than one signed summary statistic column (%s and %s).\nPick one with --signed-sumstats, or drop one with --ignore.\n", mopts->fname, signed_col_name.c_str(), col_names[i].c_str());
          return 1;
        }
        signed_field = field;
        signed_col_name = col_names[i];
        if (field != kMungeFieldSigned) {
          signed_null = LdscSignedNull(field);
        }
      }
      if (signed_field == kMungeFieldNone) {
        fprintf(stderr, "Error: Could not find a signed summary statistic column in %s (Z, OR, BETA\nor LOG_ODDS).  Name one with --signed-sumstats, or pass --a1-inc if A1 is\nalways the trait-increasing allele.\n", mopts->fname);
        return 1;
      }
    }
  }
  // Column -> row field.
  uint32_t col_snp = UINT32_MAX;
  uint32_t col_p = UINT32_MAX;
  uint32_t col_a1 = UINT32_MAX;
  uint32_t col_a2 = UINT32_MAX;
  uint32_t col_n = UINT32_MAX;
  uint32_t col_n_cas = UINT32_MAX;
  uint32_t col_n_con = UINT32_MAX;
  uint32_t col_nstudy = UINT32_MAX;
  uint32_t col_frq = UINT32_MAX;
  uint32_t col_signed = UINT32_MAX;
  std::vector<uint32_t> info_cols;
  for (uint32_t i = 0; i != col_ct; ++i) {
    switch (col_fields[i]) {
    case kMungeFieldSnp: col_snp = i; break;
    case kMungeFieldP: col_p = i; break;
    case kMungeFieldA1: col_a1 = i; break;
    case kMungeFieldA2: col_a2 = i; break;
    case kMungeFieldN: col_n = i; break;
    case kMungeFieldNCas: col_n_cas = i; break;
    case kMungeFieldNCon: col_n_con = i; break;
    case kMungeFieldNStudy: col_nstudy = i; break;
    case kMungeFieldFrq: col_frq = i; break;
    case kMungeFieldInfo: info_cols.push_back(i); break;
    default:
      if (col_fields[i] == signed_field) {
        col_signed = i;
      }
      break;
    }
  }
  if (mopts->daner) {
    col_n = UINT32_MAX;
    col_n_cas = UINT32_MAX;
    col_n_con = UINT32_MAX;
  }
  if (col_snp == UINT32_MAX) {
    fprintf(stderr, "Error: Could not find a variant ID column in %s.\n", mopts->fname);
    return 1;
  }
  if (col_p == UINT32_MAX) {
    fprintf(stderr, "Error: Could not find a p-value column in %s.\n", mopts->fname);
    return 1;
  }
  if ((!mopts->no_alleles) && ((col_a1 == UINT32_MAX) || (col_a2 == UINT32_MAX))) {
    fprintf(stderr, "Error: Could not find A1/A2 columns in %s.  Pass --no-alleles if the file\nreally has none.\n", mopts->fname);
    return 1;
  }
  // The reference implementation drops NSTUDY once N is available.
  const uint32_t have_n_cols = (col_n != UINT32_MAX) || ((col_n_cas != UINT32_MAX) && (col_n_con != UINT32_MAX));
  if (have_n_cols) {
    col_nstudy = UINT32_MAX;
  }
  const uint32_t have_n_flags = (mopts->n_override > 0.0) || ((mopts->n_cas_override > 0.0) && (mopts->n_con_override > 0.0)) || mopts->daner;
  if ((!have_n_cols) && (!have_n_flags) && (col_nstudy == UINT32_MAX)) {
    fprintf(stderr, "Error: Could not determine the sample size for %s.  Pass --N, or --N-cas\nwith --N-con.\n", mopts->fname);
    return 1;
  }

  LdscLog("Interpreting %s columns as follows:\n", mopts->fname);
  for (uint32_t i = 0; i != col_ct; ++i) {
    const char* desc = nullptr;
    switch (col_fields[i]) {
    case kMungeFieldSnp: desc = (col_snp == i)? "variant ID" : nullptr; break;
    case kMungeFieldP: desc = (col_p == i)? "p-value" : nullptr; break;
    case kMungeFieldA1: desc = "A1 (the allele the signed statistic refers to)"; break;
    case kMungeFieldA2: desc = "A2"; break;
    case kMungeFieldN: desc = (col_n == i)? "sample size" : nullptr; break;
    case kMungeFieldNCas: desc = (col_n_cas == i)? "case count" : nullptr; break;
    case kMungeFieldNCon: desc = (col_n_con == i)? "control count" : nullptr; break;
    case kMungeFieldNStudy: desc = (col_nstudy == i)? "number of studies" : nullptr; break;
    case kMungeFieldFrq: desc = "allele frequency"; break;
    case kMungeFieldInfo: desc = "imputation INFO score"; break;
    default:
      if (col_signed == i) {
        desc = "signed summary statistic";
      }
      break;
    }
    if (desc) {
      LdscLog("  %s:\t%s\n", col_names[i].c_str(), desc);
    }
  }
  if (col_signed != UINT32_MAX) {
    LdscLog("  (%s is signed, with %g meaning no effect.)\n", signed_col_name.c_str(), signed_null);
  }

  uint32_t max_col = MAXV(col_snp, col_p);
  const uint32_t opt_cols[] = {col_a1, col_a2, col_n, col_n_cas, col_n_con, col_nstudy, col_frq, col_signed};
  for (uint32_t i = 0; i != sizeof(opt_cols) / sizeof(opt_cols[0]); ++i) {
    if (opt_cols[i] != UINT32_MAX) {
      max_col = MAXV(max_col, opt_cols[i]);
    }
  }
  for (uintptr_t i = 0; i != info_cols.size(); ++i) {
    max_col = MAXV(max_col, info_cols[i]);
  }

  std::vector<MungeRow> rows;
  uintptr_t read_ct = 0;
  uintptr_t drop_na = 0;
  uintptr_t drop_merge = 0;
  uintptr_t drop_info = 0;
  uintptr_t drop_frq = 0;
  uintptr_t drop_p = 0;
  uintptr_t drop_alleles = 0;
  uintptr_t bad_info_ct = 0;
  uintptr_t bad_frq_ct = 0;
  uintptr_t bad_p_ct = 0;
  while (1) {
    const char* line_start = TextGet(&txs);
    if (!line_start) {
      break;
    }
    const char* iter = FirstNonTspace(line_start);
    if (IsEolnKns(*iter)) {
      continue;
    }
    ++read_ct;
    MungeRow row;
    row.a1 = '\0';
    row.a2 = '\0';
    row.p = 0.0;
    row.signed_stat = 0.0;
    row.n = -1.0;
    row.n_cas = -1.0;
    row.n_con = -1.0;
    row.frq = -1.0;
    row.nstudy = -1.0;
    double info_sum = 0.0;
    uint32_t info_ct = 0;
    const char* id_start = nullptr;
    uint32_t id_slen = 0;
    uint32_t ok = 1;
    uint32_t info_idx = 0;
    for (uint32_t col_idx = 0; col_idx <= max_col; ++col_idx) {
      if (IsEolnKns(*iter)) {
        ok = 0;
        break;
      }
      const char* token_end = CurTokenEnd(iter);
      const uint32_t slen = token_end - iter;
      const uint32_t is_missing = ((slen == 1) && ((*iter == '.') || (*iter == '?'))) || ((slen == 2) && (((iter[0] == 'N') && (iter[1] == 'A')) || ((iter[0] == 'n') && (iter[1] == 'a'))));
      if (col_idx == col_snp) {
        id_start = iter;
        id_slen = slen;
      } else if (col_idx == col_p) {
        if (is_missing || (!ScanadvDouble(iter, &row.p))) {
          ok = 0;
        }
      } else if (col_idx == col_a1) {
        if (is_missing || (slen != 1)) {
          ok = 0;
        } else {
          row.a1 = LdscUpcase(*iter);
        }
      } else if (col_idx == col_a2) {
        if (is_missing || (slen != 1)) {
          ok = 0;
        } else {
          row.a2 = LdscUpcase(*iter);
        }
      } else if (col_idx == col_signed) {
        if (is_missing || (!ScanadvDouble(iter, &row.signed_stat))) {
          ok = 0;
        }
      } else if (col_idx == col_n) {
        if (is_missing || (!ScanadvDouble(iter, &row.n))) {
          ok = 0;
        }
      } else if (col_idx == col_n_cas) {
        if (is_missing || (!ScanadvDouble(iter, &row.n_cas))) {
          ok = 0;
        }
      } else if (col_idx == col_n_con) {
        if (is_missing || (!ScanadvDouble(iter, &row.n_con))) {
          ok = 0;
        }
      } else if (col_idx == col_nstudy) {
        if (is_missing || (!ScanadvDouble(iter, &row.nstudy))) {
          ok = 0;
        }
      } else if (col_idx == col_frq) {
        if (is_missing || (!ScanadvDouble(iter, &row.frq))) {
          ok = 0;
        }
      } else if ((info_idx != info_cols.size()) && (col_idx == info_cols[info_idx])) {
        double cur;
        // An unparsable INFO is left out of the average rather than dropping
        // the variant, as in the reference implementation.
        if ((!is_missing) && ScanadvDouble(iter, &cur)) {
          info_sum += cur;
          ++info_ct;
        }
        ++info_idx;
      }
      if (!ok) {
        break;
      }
      iter = FirstNonTspace(token_end);
    }
    if (!ok) {
      ++drop_na;
      continue;
    }
    row.id.assign(id_start, id_slen);
    if (mopts->merge_alleles) {
      if (merge_list.id_to_idx.find(row.id) == merge_list.id_to_idx.end()) {
        ++drop_merge;
        continue;
      }
    }
    if (!info_cols.empty()) {
      if (info_ct != info_cols.size()) {
        // Missing INFO is not a reason to drop a variant.
        info_ct = info_cols.size();
      }
      const double info_mean = info_sum / u31tod(info_ct);
      if ((info_mean > 2.0) || (info_mean < 0.0)) {
        ++bad_info_ct;
      }
      if (!(info_mean >= mopts->info_min)) {
        ++drop_info;
        continue;
      }
    }
    if (col_frq != UINT32_MAX) {
      if ((row.frq < 0.0) || (row.frq > 1.0)) {
        ++bad_frq_ct;
        ++drop_frq;
        continue;
      }
      if (!(MINV(row.frq, 1.0 - row.frq) > mopts->maf_min)) {
        ++drop_frq;
        continue;
      }
    }
    if ((!(row.p > 0.0)) || (row.p > 1.0)) {
      ++bad_p_ct;
      ++drop_p;
      continue;
    }
    if (!mopts->no_alleles) {
      if (!LdscMungeValidSnp(row.a1, row.a2)) {
        ++drop_alleles;
        continue;
      }
    }
    rows.push_back(row);
  }
  reterr = kPglRetSuccess;
  CleanupTextStream(&txs, &reterr);
  if (reterr) {
    fprintf(stderr, "Error: Failed to read %s.\n", mopts->fname);
    return 1;
  }
  LdscLog("Read %" PRIuPTR " variants from %s.\n", read_ct, mopts->fname);
  if (mopts->merge_alleles) {
    LdscLog("Removed %" PRIuPTR " variants not in --merge-alleles.\n", drop_merge);
  }
  LdscLog("Removed %" PRIuPTR " variants with missing values.\n", drop_na);
  if (!info_cols.empty()) {
    if (bad_info_ct) {
      LdscLog("WARNING: %" PRIuPTR " variants had INFO outside [0, 2].  The INFO column may be\nmislabeled.\n", bad_info_ct);
    }
    LdscLog("Removed %" PRIuPTR " variants with INFO < %g.\n", drop_info, mopts->info_min);
  }
  if (col_frq != UINT32_MAX) {
    if (bad_frq_ct) {
      LdscLog("WARNING: %" PRIuPTR " variants had a frequency outside [0, 1].  The frequency\ncolumn may be mislabeled.\n", bad_frq_ct);
    }
    LdscLog("Removed %" PRIuPTR " variants with MAF <= %g.\n", drop_frq, mopts->maf_min);
  }
  if (bad_p_ct) {
    LdscLog("WARNING: %" PRIuPTR " variants had p outside (0, 1].  The p-value column may be\nmislabeled.\n", bad_p_ct);
  }
  LdscLog("Removed %" PRIuPTR " variants with out-of-bounds p-values.\n", drop_p);
  if (!mopts->no_alleles) {
    LdscLog("Removed %" PRIuPTR " variants that were not SNPs or were strand-ambiguous.\n", drop_alleles);
  }
  if (rows.empty()) {
    fprintf(stderr, "Error: No variants remain after filtering.\n");
    return 1;
  }

  // Duplicated IDs: keep the first.
  {
    std::unordered_map<std::string, uint32_t> seen;
    uintptr_t write_idx = 0;
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      if (!seen.emplace(rows[i].id, 1).second) {
        continue;
      }
      if (write_idx != i) {
        rows[write_idx] = rows[i];
      }
      ++write_idx;
    }
    const uintptr_t dup_ct = rows.size() - write_idx;
    rows.resize(write_idx);
    LdscLog("Removed %" PRIuPTR " variants with duplicated IDs (%" PRIuPTR " remain).\n", dup_ct, S_CAST(uintptr_t, rows.size()));
  }

  // Sample size.  With case and control counts, the effective size is scaled
  // by the case fraction relative to the best-powered variants, as in the
  // reference implementation.
  if ((col_n_cas != UINT32_MAX) && (col_n_con != UINT32_MAX)) {
    double max_n = 0.0;
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      max_n = MAXV(max_n, rows[i].n_cas + rows[i].n_con);
    }
    double frac_sum = 0.0;
    uint32_t frac_ct = 0;
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      const double cur_n = rows[i].n_cas + rows[i].n_con;
      if (cur_n == max_n) {
        frac_sum += rows[i].n_cas / cur_n;
        ++frac_ct;
      }
    }
    const double ref_frac = frac_sum / u31tod(frac_ct);
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      const double cur_n = rows[i].n_cas + rows[i].n_con;
      rows[i].n = cur_n * (rows[i].n_cas / cur_n) / ref_frac;
    }
  }
  if ((col_n != UINT32_MAX) || ((col_n_cas != UINT32_MAX) && (col_n_con != UINT32_MAX))) {
    double n_min = mopts->n_min;
    if (!mopts->have_n_min) {
      // The reference implementation's default: the 90th percentile over 1.5.
      std::vector<double> ns(rows.size());
      for (uintptr_t i = 0; i != rows.size(); ++i) {
        ns[i] = rows[i].n;
      }
      std::sort(ns.begin(), ns.end());
      // numpy's default quantile interpolation.
      const double pos = 0.9 * u31tod(ns.size() - 1);
      const uintptr_t lo = S_CAST(uintptr_t, pos);
      const double frac = pos - u31tod(lo);
      const double q90 = (lo + 1 < ns.size())? (ns[lo] * (1 - frac) + ns[lo + 1] * frac) : ns[lo];
      n_min = q90 / 1.5;
    }
    uintptr_t write_idx = 0;
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      if (rows[i].n < n_min) {
        continue;
      }
      if (write_idx != i) {
        rows[write_idx] = rows[i];
      }
      ++write_idx;
    }
    const uintptr_t removed = rows.size() - write_idx;
    rows.resize(write_idx);
    LdscLog("Removed %" PRIuPTR " variants with N < %g (%" PRIuPTR " remain).\n", removed, n_min, S_CAST(uintptr_t, rows.size()));
  } else if (col_nstudy != UINT32_MAX) {
    double nstudy_min = mopts->nstudy_min;
    if (!mopts->have_nstudy_min) {
      nstudy_min = 0.0;
      for (uintptr_t i = 0; i != rows.size(); ++i) {
        nstudy_min = MAXV(nstudy_min, rows[i].nstudy);
      }
    }
    uintptr_t write_idx = 0;
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      if (rows[i].nstudy < nstudy_min) {
        continue;
      }
      if (write_idx != i) {
        rows[write_idx] = rows[i];
      }
      ++write_idx;
    }
    const uintptr_t removed = rows.size() - write_idx;
    rows.resize(write_idx);
    LdscLog("Removed %" PRIuPTR " variants genotyped in fewer than %g studies (%" PRIuPTR "\nremain).\n", removed, nstudy_min, S_CAST(uintptr_t, rows.size()));
  }
  if (rows.empty()) {
    fprintf(stderr, "Error: No variants remain after filtering.\n");
    return 1;
  }
  if ((col_n == UINT32_MAX) && ((col_n_cas == UINT32_MAX) || (col_n_con == UINT32_MAX))) {
    double n_val;
    if (mopts->daner) {
      n_val = daner_n_cas + daner_n_con;
    } else if (mopts->n_override > 0.0) {
      n_val = mopts->n_override;
      LdscLog("Using N = %g.\n", n_val);
    } else {
      n_val = mopts->n_cas_override + mopts->n_con_override;
      LdscLog("Using N_cas = %g, N_con = %g.\n", mopts->n_cas_override, mopts->n_con_override);
    }
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      rows[i].n = n_val;
    }
  } else if (mopts->n_override > 0.0) {
    // An explicit --N takes priority over the column.
    LdscLog("Using N = %g (overriding the file's sample size column).\n", mopts->n_override);
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      rows[i].n = mopts->n_override;
    }
  }

  // The signed statistic's median is a check on the column's meaning: a
  // mislabeled column usually has the wrong median.
  if (col_signed != UINT32_MAX) {
    std::vector<double> stats(rows.size());
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      stats[i] = rows[i].signed_stat;
    }
    const double median = LdscMedian(&(stats[0]), stats.size());
    if (fabs(median - signed_null) > 0.1) {
      fprintf(stderr, "Error: the median of %s is %g, but %g means no effect, so the column looks\nmislabeled.  Use --signed-sumstats to say what it is, or --ignore to drop it.\n", signed_col_name.c_str(), median, signed_null);
      return 1;
    }
    LdscLog("Median %s was %g, which seems sensible.\n", signed_col_name.c_str(), median);
  }

  // p-value and sign to Z.
  std::vector<double> zs(rows.size());
  for (uintptr_t i = 0; i != rows.size(); ++i) {
    double z = LdscNormalIsf(rows[i].p / 2.0);
    if ((col_signed != UINT32_MAX) && (rows[i].signed_stat < signed_null)) {
      z = -z;
    }
    zs[i] = z;
  }

  // Whether every sample size came from a column and is a whole number,
  // which is what decides the reference implementation's number format.
  uint32_t n_is_integral_column = (col_n != UINT32_MAX) && (mopts->n_override <= 0.0);
  if (n_is_integral_column) {
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      if (rows[i].n != floor(rows[i].n)) {
        n_is_integral_column = 0;
        break;
      }
    }
  }

  const std::string out_path = std::string(out_prefix) + ".sumstats";
  FILE* out_file = fopen(out_path.c_str(), FOPEN_WB);
  if (!out_file) {
    fprintf(stderr, "Error: Failed to open %s.\n", out_path.c_str());
    return 1;
  }
  const uint32_t write_frq = mopts->keep_maf && (col_frq != UINT32_MAX);
  fputs("SNP", out_file);
  if (!mopts->no_alleles) {
    fputs("\tA1\tA2", out_file);
  }
  fputs("\tZ\tN", out_file);
  if (write_frq) {
    fputs("\tFRQ", out_file);
  }
  fputc('\n', out_file);
  uintptr_t written_ct = 0;
  uintptr_t nonmissing_ct = 0;
  if (!mopts->merge_alleles) {
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      fputs(rows[i].id.c_str(), out_file);
      if (!mopts->no_alleles) {
        fprintf(out_file, "\t%c\t%c", rows[i].a1, rows[i].a2);
      }
      fprintf(out_file, "\t%.3f\t", zs[i]);
      LdscWriteN(rows[i].n, n_is_integral_column, out_file);
      if (write_frq) {
        fprintf(out_file, "\t%.3f", rows[i].frq);
      }
      fputc('\n', out_file);
      ++written_ct;
      ++nonmissing_ct;
    }
  } else {
    // Every variant in the merge list is written, in its order; the ones the
    // input did not supply, or whose alleles disagree with the list, come out
    // missing.  The alleles and the sign stay as the input file had them, as
    // in the reference implementation: --merge-alleles restricts and
    // validates, and the actual reorientation happens when two files are
    // regressed against each other.
    std::vector<uint32_t> row_of(merge_list.ids.size(), UINT32_MAX);
    uintptr_t mismatch_ct = 0;
    for (uintptr_t i = 0; i != rows.size(); ++i) {
      const std::unordered_map<std::string, uint32_t>::const_iterator it = merge_list.id_to_idx.find(rows[i].id);
      if (it == merge_list.id_to_idx.end()) {
        continue;
      }
      const uint32_t merge_idx = it->second;
      const char m1 = merge_list.a1[merge_idx];
      const char m2 = merge_list.a2[merge_idx];
      const char c1 = LdscAcgtComplement(rows[i].a1);
      const char c2 = LdscAcgtComplement(rows[i].a2);
      const uint32_t same = ((rows[i].a1 == m1) && (rows[i].a2 == m2)) || ((c1 == m1) && (c2 == m2));
      const uint32_t flipped = ((rows[i].a1 == m2) && (rows[i].a2 == m1)) || ((c1 == m2) && (c2 == m1));
      if ((!same) && (!flipped)) {
        ++mismatch_ct;
        continue;
      }
      row_of[merge_idx] = i;
    }
    LdscLog("Removed %" PRIuPTR " variants whose alleles did not match --merge-alleles.\n", mismatch_ct);
    for (uintptr_t merge_idx = 0; merge_idx != merge_list.ids.size(); ++merge_idx) {
      fputs(merge_list.ids[merge_idx].c_str(), out_file);
      const uint32_t row_idx = row_of[merge_idx];
      if (row_idx == UINT32_MAX) {
        fputs("\tNA\tNA\tNA\tNA", out_file);
        if (write_frq) {
          fputs("\tNA", out_file);
        }
      } else {
        fprintf(out_file, "\t%c\t%c", rows[row_idx].a1, rows[row_idx].a2);
        fprintf(out_file, "\t%.3f\t", zs[row_idx]);
        LdscWriteN(rows[row_idx].n, n_is_integral_column, out_file);
        if (write_frq) {
          fprintf(out_file, "\t%.3f", rows[row_idx].frq);
        }
        ++nonmissing_ct;
      }
      fputc('\n', out_file);
      ++written_ct;
    }
  }
  if (fclose(out_file)) {
    fprintf(stderr, "Error: Failed to write %s.\n", out_path.c_str());
    return 1;
  }
  LdscLog("Wrote summary statistics for %" PRIuPTR " variants (%" PRIuPTR " with a nonmissing Z)\nto %s .\n", written_ct, nonmissing_ct, out_path.c_str());

  // Metadata, as the reference implementation reports it.
  std::vector<double> chisq(rows.size());
  double chisq_sum = 0.0;
  double chisq_max = 0.0;
  uintptr_t gws_ct = 0;
  for (uintptr_t i = 0; i != rows.size(); ++i) {
    chisq[i] = zs[i] * zs[i];
    chisq_sum += chisq[i];
    chisq_max = MAXV(chisq_max, chisq[i]);
    if (chisq[i] > 29) {
      ++gws_ct;
    }
  }
  const double mean_chisq = chisq_sum / u31tod(rows.size());
  LdscLog("\nMetadata:\n");
  LdscLog("Mean chi^2 = %.3f\n", mean_chisq);
  if (mean_chisq < 1.02) {
    LdscLog("WARNING: mean chi^2 may be too small.\n");
  }
  LdscLog("Lambda GC = %.3f\n", LdscMedian(&(chisq[0]), chisq.size()) / 0.4549);
  LdscLog("Max chi^2 = %.3f\n", chisq_max);
  LdscLog("%" PRIuPTR " genome-wide significant variants (some may have been filtered out).\n", gws_ct);
  return 0;
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
  uint32_t overlap_annot;
  const char* frq_arg;
  uint32_t frq_is_chr_split;
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
  const uint32_t n_annot = data->n_annot;
  const uint32_t orig_ct = data->ld_tot.size();
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
    for (uint32_t j = 0; j != n_annot; ++j) {
      data->ld[S_CAST(uintptr_t, write_idx) * n_annot + j] = data->ld[S_CAST(uintptr_t, i) * n_annot + j];
    }
    data->ld_tot[write_idx] = data->ld_tot[i];
    data->w_ld[write_idx] = data->w_ld[i];
    data->z1[write_idx] = data->z1[i];
    data->n1[write_idx] = data->n1[i];
    if (use_both_z) {
      data->z2[write_idx] = data->z2[i];
      data->n2[write_idx] = data->n2[i];
    }
    ++write_idx;
  }
  data->ld.resize(S_CAST(uintptr_t, write_idx) * n_annot);
  data->ld_tot.resize(write_idx);
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

BoolErr LdscEstimateH2(const char* sumstats_fname, const LdscScores& ref, const char* ref_ld_arg, uint32_t ref_ld_chr_split, const std::unordered_map<std::string, double>& w_ld_map, const std::vector<double>& m_vec, const LdscOpts* opts, const char* out_prefix) {
  const uint32_t n_annot = ref.annot_names.size();
  double m_tot = 0.0;
  for (uint32_t j = 0; j != n_annot; ++j) {
    m_tot += m_vec[j];
  }
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
  LdscMerge(ref, w_ld_map, sumstats, 0, &data);
  uint32_t n_snp = data.ld_tot.size();
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
  uint32_t have_chisq_max = opts->have_chisq_max;
  double chisq_max = opts->chisq_max;
  if ((!have_chisq_max) && (n_annot > 1)) {
    // With more than one annotation the reference implementation drops the
    // high-chi^2 tail instead of running the two-step estimator, since the
    // latter is not defined for a partitioned regression.
    double max_n = 0.0;
    for (uint32_t i = 0; i != n_snp; ++i) {
      max_n = MAXV(max_n, data.n1[i]);
    }
    chisq_max = MAXV(0.001 * max_n, 80.0);
    have_chisq_max = 1;
  }
  if (have_chisq_max) {
    LdscFilterChisq(chisq_max, &data, 0);
    n_snp = data.ld_tot.size();
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
  } else if ((!fixed_intercept) && (n_annot == 1)) {
    two_step = &two_step_val;
  }
  if (two_step && fixed_intercept) {
    LdscLog("Ignoring --two-step: it only applies when the intercept is free.\n");
    two_step = nullptr;
  } else if (two_step && (n_annot > 1)) {
    LdscLog("Ignoring --two-step: it is not defined for a partitioned regression.\n");
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
  if (LdscHsqFit(&(chisq[0]), &(data.ld[0]), &(data.ld_tot[0]), &(data.w_ld[0]), &(data.n1[0]), n_snp, n_annot, &(m_vec[0]), n_blocks, fixed_intercept, two_step, &hsq)) {
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

  if (n_annot > 1) {
    // Per-annotation output, in the reference implementation's .results
    // columns.  Without --overlap-annot the proportions and enrichments are
    // only right for a partition of the variants; with it they are corrected
    // using the annotation overlap matrix from the .annot files.
    std::vector<double> prop = hsq.prop;
    std::vector<double> prop_ses = hsq.prop_ses;
    std::vector<double> enrichment = hsq.enrichment;
    std::vector<double> enrichment_ses;
    std::vector<double> enrichment_ln_p;
    std::vector<double> m_prop = hsq.m_prop;
    if (opts->overlap_annot) {
      std::vector<double> overlap;
      double annot_m_tot;
      if (LdscReadAnnot(ref_ld_arg, ref_ld_chr_split, opts->frq_arg, opts->frq_is_chr_split, n_annot, ref.per_file_annot_ct, &overlap, &annot_m_tot)) {
        return 1;
      }
      LdscLog("Read annotation overlap for %g variants.\n", annot_m_tot);
      // overlap_prop[i][j] is the fraction of annotation j's variants that
      // are also in annotation i.
      std::vector<double> overlap_prop(S_CAST(uintptr_t, n_annot) * n_annot);
      for (uint32_t i = 0; i != n_annot; ++i) {
        for (uint32_t j = 0; j != n_annot; ++j) {
          overlap_prop[S_CAST(uintptr_t, i) * n_annot + j] = overlap[S_CAST(uintptr_t, i) * n_annot + j] / m_vec[j];
        }
      }
      for (uint32_t i = 0; i != n_annot; ++i) {
        double acc = 0.0;
        for (uint32_t j = 0; j != n_annot; ++j) {
          acc += overlap_prop[S_CAST(uintptr_t, i) * n_annot + j] * hsq.prop[j];
        }
        prop[i] = acc;
        double var = 0.0;
        for (uint32_t j = 0; j != n_annot; ++j) {
          for (uint32_t k = 0; k != n_annot; ++k) {
            var += overlap_prop[S_CAST(uintptr_t, i) * n_annot + j] * hsq.prop_cov[S_CAST(uintptr_t, j) * n_annot + k] * overlap_prop[S_CAST(uintptr_t, i) * n_annot + k];
          }
        }
        prop_ses[i] = sqrt(MAXV(var, 0.0));
        m_prop[i] = m_vec[i] / annot_m_tot;
      }
      enrichment.resize(n_annot);
      enrichment_ses.resize(n_annot);
      enrichment_ln_p.assign(n_annot, 0.0 / 0.0);
      for (uint32_t i = 0; i != n_annot; ++i) {
        enrichment[i] = prop[i] / m_prop[i];
        enrichment_ses[i] = prop_ses[i] / m_prop[i];
      }
      // The enrichment test compares each annotation's coefficient against
      // the rest of the genome, which is a linear combination of the
      // coefficients: inside-annotation share minus outside-annotation share.
      std::vector<double> diff(S_CAST(uintptr_t, n_annot) * n_annot, 0.0);
      for (uint32_t i = 0; i != n_annot; ++i) {
        if (annot_m_tot == m_vec[i]) {
          continue;
        }
        for (uint32_t j = 0; j != n_annot; ++j) {
          const double cur_overlap = overlap[S_CAST(uintptr_t, i) * n_annot + j];
          diff[S_CAST(uintptr_t, i) * n_annot + j] = cur_overlap / m_vec[i] - (m_vec[j] - cur_overlap) / (annot_m_tot - m_vec[i]);
        }
      }
      for (uint32_t i = 0; i != n_annot; ++i) {
        double est = 0.0;
        for (uint32_t j = 0; j != n_annot; ++j) {
          est += diff[S_CAST(uintptr_t, i) * n_annot + j] * hsq.coefs[j];
        }
        double var = 0.0;
        for (uint32_t j = 0; j != n_annot; ++j) {
          for (uint32_t k = 0; k != n_annot; ++k) {
            var += diff[S_CAST(uintptr_t, i) * n_annot + j] * hsq.coef_cov[S_CAST(uintptr_t, j) * n_annot + k] * diff[S_CAST(uintptr_t, i) * n_annot + k];
          }
        }
        const double se = sqrt(MAXV(var, 0.0));
        if (se > 0.0) {
          enrichment_ln_p[i] = TstatToLnP(est / se, hsq.n_blocks);
        }
      }
    }
    const std::string results_path = std::string(out_prefix) + ".results";
    FILE* results_file = fopen(results_path.c_str(), FOPEN_WB);
    if (!results_file) {
      fprintf(stderr, "Error: Failed to open %s.\n", results_path.c_str());
      return 1;
    }
    fputs("Category\tProp._SNPs\tProp._h2\tProp._h2_std_error\tEnrichment", results_file);
    if (opts->overlap_annot) {
      fputs("\tEnrichment_std_error\tEnrichment_p", results_file);
    }
    fputs("\tCoefficient\tCoefficient_std_error\tCoefficient_z-score\n", results_file);
    for (uint32_t j = 0; j != n_annot; ++j) {
      fprintf(results_file, "%s\t%.9g\t%.9g\t%.9g\t%.9g", ref.annot_names[j].c_str(), m_prop[j], prop[j], prop_ses[j], enrichment[j]);
      if (opts->overlap_annot) {
        fprintf(results_file, "\t%.9g", enrichment_ses[j]);
        if (enrichment_ln_p[j] == enrichment_ln_p[j]) {
          fprintf(results_file, "\t%.9g", exp(enrichment_ln_p[j]));
        } else {
          // An annotation that covers every variant has no complement to be
          // compared against.
          fputs("\tNA", results_file);
        }
      }
      fprintf(results_file, "\t%.9g\t%.9g\t%.9g\n", liab_factor * hsq.coefs[j], liab_factor * hsq.coef_ses[j], hsq.coefs[j] / hsq.coef_ses[j]);
    }
    if (fclose(results_file)) {
      fprintf(stderr, "Error: Failed to write %s.\n", results_path.c_str());
      return 1;
    }
    LdscLog("Per-category results written to %s .\n", results_path.c_str());
    if (!opts->overlap_annot) {
      LdscLog("Note: the proportions and enrichments there are only right for annotations\nthat do not overlap.  Add --overlap-annot (with the .annot files, and\n--frqfile[-chr] to restrict to common variants) to correct for overlap.\n");
    }
  }
  return 0;
}

// Merges a second trait into the first trait's rows, flipping its Z where the
// alleles are swapped and dropping variants whose alleles do not match.
BoolErr LdscMergeSecondTrait(const LdscData* base, const std::unordered_map<std::string, LdscSumstatRow>& sumstats2, uint32_t no_check_alleles, LdscData* dst) {
  const uint32_t n_base = base->ld_tot.size();
  const uint32_t n_annot = base->n_annot;
  dst->n_annot = n_annot;
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
    for (uint32_t j = 0; j != n_annot; ++j) {
      dst->ld.push_back(base->ld[S_CAST(uintptr_t, i) * n_annot + j]);
    }
    dst->ld_tot.push_back(base->ld_tot[i]);
    dst->w_ld.push_back(base->w_ld[i]);
    dst->z1.push_back(base->z1[i]);
    dst->n1.push_back(base->n1[i]);
    dst->z2.push_back(z2);
    dst->n2.push_back(it->second.n);
  }
  if (bad_allele_ct) {
    LdscLog("Dropped %u variants with mismatched or strand-ambiguous alleles.\n", bad_allele_ct);
  }
  if (dst->ld_tot.empty()) {
    fprintf(stderr, "Error: No variants in common between the two summary statistic files.\n");
    return 1;
  }
  return 0;
}

BoolErr LdscEstimateRg(const std::vector<std::string>& rg_paths, const LdscScores& ref, const std::unordered_map<std::string, double>& w_ld_map, const std::vector<double>& m_vec, const LdscOpts* opts, const char* out_prefix) {
  const uint32_t n_annot = ref.annot_names.size();
  const uint32_t pheno_ct = rg_paths.size();
  std::unordered_map<std::string, LdscSumstatRow> sumstats1;
  uint32_t dropped_ct;
  if (LdscReadSumstats(rg_paths[0].c_str(), 1, &sumstats1, &dropped_ct)) {
    return 1;
  }
  LdscLog("Read summary statistics for %" PRIuPTR " variants from %s.\n", S_CAST(uintptr_t, sumstats1.size()), rg_paths[0].c_str());
  LdscData base;
  LdscMerge(ref, w_ld_map, sumstats1, 1, &base);
  if (base.ld_tot.empty()) {
    fprintf(stderr, "Error: No variants remain after merging %s with the LD Scores.\n", rg_paths[0].c_str());
    return 1;
  }
  LdscLog("After merging with reference panel LD and regression weight LD, %" PRIuPTR "\nvariants remain.\n", S_CAST(uintptr_t, base.ld_tot.size()));

  double two_step_val = 30.0;
  const double* two_step = nullptr;
  const uint32_t intercept_h2_free = (!opts->no_intercept) && (!LdscOptAt(opts->intercept_h2, 0));
  if (opts->have_two_step) {
    two_step_val = opts->two_step;
    two_step = &two_step_val;
  } else if (intercept_h2_free && (n_annot == 1)) {
    two_step = &two_step_val;
  }
  if (two_step && (n_annot > 1)) {
    LdscLog("Ignoring --two-step: it is not defined for a partitioned regression.\n");
    two_step = nullptr;
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
    uint32_t n_snp = data.ld_tot.size();
    LdscLog("%u variants with valid alleles in both files.\n", n_snp);
    if (opts->have_chisq_max) {
      LdscFilterChisq(opts->chisq_max, &data, 1);
      n_snp = data.ld_tot.size();
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
    if (LdscHsqFit(&(chisq1[0]), &(data.ld[0]), &(data.ld_tot[0]), &(data.w_ld[0]), &(data.n1[0]), n_snp, n_annot, &(m_vec[0]), n_blocks, fixed_h2_1, two_step, &hsq1) ||
        LdscHsqFit(&(chisq2[0]), &(data.ld[0]), &(data.ld_tot[0]), &(data.w_ld[0]), &(data.n2[0]), n_snp, n_annot, &(m_vec[0]), n_blocks, fixed_h2_2, two_step, &hsq2)) {
      fprintf(stderr, "Error: Heritability regression failed for phenotype pair 1/%u.\n", pheno_idx + 1);
      return 1;
    }
    LdscGencovResult gencov;
    if (LdscGencovFit(&(data.z1[0]), &(data.z2[0]), &(data.ld[0]), &(data.ld_tot[0]), &(data.w_ld[0]), &(data.n1[0]), &(data.n2[0]), n_snp, n_annot, &(m_vec[0]), n_blocks, hsq1.tot, hsq2.tot, hsq1.intercept, hsq2.intercept, fixed_gencov, two_step, &gencov)) {
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
  MungeOpts mopts;
  mopts.fname = nullptr;
  mopts.merge_alleles = nullptr;
  mopts.snp_col = nullptr;
  mopts.a1_col = nullptr;
  mopts.a2_col = nullptr;
  mopts.p_col = nullptr;
  mopts.n_col = nullptr;
  mopts.n_cas_col = nullptr;
  mopts.n_con_col = nullptr;
  mopts.frq_col = nullptr;
  mopts.info_col = nullptr;
  mopts.info_list = nullptr;
  mopts.nstudy_col = nullptr;
  mopts.signed_sumstats = nullptr;
  mopts.ignore = nullptr;
  mopts.n_override = 0.0;
  mopts.n_cas_override = 0.0;
  mopts.n_con_override = 0.0;
  mopts.info_min = 0.9;
  mopts.maf_min = 0.01;
  mopts.n_min = 0.0;
  mopts.have_n_min = 0;
  mopts.nstudy_min = 0.0;
  mopts.have_nstudy_min = 0;
  mopts.a1_inc = 0;
  mopts.no_alleles = 0;
  mopts.keep_maf = 0;
  mopts.daner = 0;
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
  opts.overlap_annot = 0;
  opts.frq_arg = nullptr;
  opts.frq_is_chr_split = 0;
  opts.m_override = 0.0;
  for (int argi = 1; argi < argc; ++argi) {
    const char* cur = argv[argi];
    if ((!strcmp(cur, "--h2")) && (argi + 1 < argc)) {
      h2_fname = argv[++argi];
    } else if ((!strcmp(cur, "--rg")) && (argi + 1 < argc)) {
      rg_arg = argv[++argi];
    } else if ((!strcmp(cur, "--munge")) && (argi + 1 < argc)) {
      mopts.fname = argv[++argi];
    } else if ((!strcmp(cur, "--merge-alleles")) && (argi + 1 < argc)) {
      mopts.merge_alleles = argv[++argi];
    } else if ((!strcmp(cur, "--snp")) && (argi + 1 < argc)) {
      mopts.snp_col = argv[++argi];
    } else if ((!strcmp(cur, "--a1")) && (argi + 1 < argc)) {
      mopts.a1_col = argv[++argi];
    } else if ((!strcmp(cur, "--a2")) && (argi + 1 < argc)) {
      mopts.a2_col = argv[++argi];
    } else if ((!strcmp(cur, "--p")) && (argi + 1 < argc)) {
      mopts.p_col = argv[++argi];
    } else if ((!strcmp(cur, "--N-col")) && (argi + 1 < argc)) {
      mopts.n_col = argv[++argi];
    } else if ((!strcmp(cur, "--N-cas-col")) && (argi + 1 < argc)) {
      mopts.n_cas_col = argv[++argi];
    } else if ((!strcmp(cur, "--N-con-col")) && (argi + 1 < argc)) {
      mopts.n_con_col = argv[++argi];
    } else if ((!strcmp(cur, "--frq")) && (argi + 1 < argc)) {
      mopts.frq_col = argv[++argi];
    } else if ((!strcmp(cur, "--info")) && (argi + 1 < argc)) {
      mopts.info_col = argv[++argi];
    } else if ((!strcmp(cur, "--info-list")) && (argi + 1 < argc)) {
      mopts.info_list = argv[++argi];
    } else if ((!strcmp(cur, "--nstudy")) && (argi + 1 < argc)) {
      mopts.nstudy_col = argv[++argi];
    } else if ((!strcmp(cur, "--signed-sumstats")) && (argi + 1 < argc)) {
      mopts.signed_sumstats = argv[++argi];
    } else if ((!strcmp(cur, "--ignore")) && (argi + 1 < argc)) {
      mopts.ignore = argv[++argi];
    } else if ((!strcmp(cur, "--N")) && (argi + 1 < argc)) {
      if ((!ScanadvDouble(argv[++argi], &mopts.n_override)) || (mopts.n_override <= 0.0)) {
        fprintf(stderr, "Error: --N must be positive.\n");
        return 1;
      }
    } else if ((!strcmp(cur, "--N-cas")) && (argi + 1 < argc)) {
      if ((!ScanadvDouble(argv[++argi], &mopts.n_cas_override)) || (mopts.n_cas_override <= 0.0)) {
        fprintf(stderr, "Error: --N-cas must be positive.\n");
        return 1;
      }
    } else if ((!strcmp(cur, "--N-con")) && (argi + 1 < argc)) {
      if ((!ScanadvDouble(argv[++argi], &mopts.n_con_override)) || (mopts.n_con_override <= 0.0)) {
        fprintf(stderr, "Error: --N-con must be positive.\n");
        return 1;
      }
    } else if ((!strcmp(cur, "--info-min")) && (argi + 1 < argc)) {
      if (!ScanadvDouble(argv[++argi], &mopts.info_min)) {
        fprintf(stderr, "Error: Invalid --info-min value.\n");
        return 1;
      }
    } else if ((!strcmp(cur, "--maf-min")) && (argi + 1 < argc)) {
      if (!ScanadvDouble(argv[++argi], &mopts.maf_min)) {
        fprintf(stderr, "Error: Invalid --maf-min value.\n");
        return 1;
      }
    } else if ((!strcmp(cur, "--n-min")) && (argi + 1 < argc)) {
      if (!ScanadvDouble(argv[++argi], &mopts.n_min)) {
        fprintf(stderr, "Error: Invalid --n-min value.\n");
        return 1;
      }
      mopts.have_n_min = 1;
    } else if ((!strcmp(cur, "--nstudy-min")) && (argi + 1 < argc)) {
      if (!ScanadvDouble(argv[++argi], &mopts.nstudy_min)) {
        fprintf(stderr, "Error: Invalid --nstudy-min value.\n");
        return 1;
      }
      mopts.have_nstudy_min = 1;
    } else if (!strcmp(cur, "--a1-inc")) {
      mopts.a1_inc = 1;
    } else if (!strcmp(cur, "--no-alleles")) {
      mopts.no_alleles = 1;
    } else if (!strcmp(cur, "--keep-maf")) {
      mopts.keep_maf = 1;
    } else if (!strcmp(cur, "--daner")) {
      mopts.daner = 1;
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
    } else if (!strcmp(cur, "--overlap-annot")) {
      opts.overlap_annot = 1;
    } else if ((!strcmp(cur, "--frqfile")) && (argi + 1 < argc)) {
      opts.frq_arg = argv[++argi];
      opts.frq_is_chr_split = 0;
    } else if ((!strcmp(cur, "--frqfile-chr")) && (argi + 1 < argc)) {
      opts.frq_arg = argv[++argi];
      opts.frq_is_chr_split = 1;
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
  if (mopts.fname) {
    if (h2_fname || rg_arg) {
      fprintf(stderr, "Error: --munge cannot be combined with --h2 or --rg.  Munge first, then\nrun the regression on the .sumstats file.\n");
      return 1;
    }
    if (!out_prefix) {
      fprintf(stderr, "Error: --munge needs --out.\n");
      return 1;
    }
    const std::string munge_log_path = std::string(out_prefix) + ".log";
    g_ldsc_logfile = fopen(munge_log_path.c_str(), FOPEN_WB);
    if (!g_ldsc_logfile) {
      fprintf(stderr, "Error: Failed to open %s.\n", munge_log_path.c_str());
      return 1;
    }
    LdscLog("%s\n", kLdscVersion);
    const BoolErr munge_ret = LdscMungeSumstats(&mopts, out_prefix);
    if (fclose(g_ldsc_logfile)) {
      fprintf(stderr, "Error: Failed to write %s.\n", munge_log_path.c_str());
      return 1;
    }
    g_ldsc_logfile = nullptr;
    return munge_ret? 1 : 0;
  }
  if ((!(h2_fname || rg_arg)) || (!ref_ld_arg) || (!w_ld_arg) || (!out_prefix)) {
    fprintf(stderr,
            "%s\n"
            "LD Score regression (Bulik-Sullivan et al. 2015).\n\n"
            "Usage: ldsc --munge <raw sumstats> --out <prefix>\n"
            "       ldsc --h2 <sumstats> --ref-ld[-chr] <LD Scores>\n"
            "            --w-ld[-chr] <LD Scores> --out <prefix>\n"
            "       ldsc --rg <sumstats1,sumstats2,...> --ref-ld[-chr] <LD Scores>\n"
            "            --w-ld[-chr] <LD Scores> --out <prefix>\n\n"
            "  --munge    Convert a raw GWAS summary statistic file into the\n"
            "             .sumstats format the regressions read, applying the\n"
            "             reference implementation's quality control: INFO,\n"
            "             frequency, p-value range, allele and sample-size\n"
            "             filters, duplicate IDs dropped, and the p-value plus\n"
            "             the signed statistic turned into a signed Z.  Columns\n"
            "             are detected case-insensitively, with the same aliases;\n"
            "             --snp/--a1/--a2/--p/--N-col/--N-cas-col/--N-con-col/\n"
            "             --frq/--info/--info-list/--nstudy name them explicitly,\n"
            "             --signed-sumstats <column>,<null value> picks the\n"
            "             signed one, and --ignore <columns> drops some.\n"
            "             Thresholds: --info-min (0.9), --maf-min (0.01),\n"
            "             --n-min (the 90th percentile of N over 1.5),\n"
            "             --nstudy-min.  --N/--N-cas/--N-con supply a sample size\n"
            "             the file does not have, --daner reads it from PGC\n"
            "             FRQ_A_/FRQ_U_ column names, --a1-inc says A1 is always\n"
            "             the trait-increasing allele, --no-alleles accepts a\n"
            "             file without them, --keep-maf keeps the frequency\n"
            "             column, and --merge-alleles <file> restricts to a\n"
            "             variant list and puts everything on its alleles.\n"
            "             Writes <prefix>.sumstats, with Z and N to three\n"
            "             decimals as the reference implementation does.\n"
            "  --h2       Estimate SNP-heritability from one summary statistic\n"
            "             file.  Columns are detected case-insensitively: an ID\n"
            "             (SNP/ID), Z, and N.\n"
            "  --rg       Estimate genetic correlation between the first file and\n"
            "             each of the others.  A1 and A2 are then required too;\n"
            "             Z is flipped where the effect alleles are swapped, and\n"
            "             strand-ambiguous variants are dropped.\n"
            "  --ref-ld   LD Scores for the regression, with SNP (or ID) and L2\n"
            "             columns.  plink2 --ld-score output can be used directly,\n"
            "             and so can ldsc's <prefix>.l2.ldscore[.gz].  Every column\n"
            "             other than the ID, the position and MAF/CM is taken as an\n"
            "             annotation, and a comma-separated list of filesets is\n"
            "             concatenated sideways, so a partitioned (stratified)\n"
            "             regression is what you get from partitioned LD Scores.\n"
            "             That needs the .l2.M_5_50 files, and turns the two-step\n"
            "             estimator off in favour of a chi^2 ceiling, as in the\n"
            "             reference implementation.\n"
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
            "  --overlap-annot  Correct the partitioned proportions and\n"
            "             enrichments for annotations that overlap, using the\n"
            "             .annot files next to --ref-ld.  Adds an enrichment\n"
            "             standard error and p-value to the .results file.\n"
            "  --frqfile[-chr]  Allele frequencies, to restrict the overlap\n"
            "             counts to variants with 5%% < MAF < 50%%, which is the\n"
            "             band the .l2.M_5_50 counts use.\n"
            "  --no-check-alleles  Skip the allele-matching step of --rg.  Only\n"
            "             safe if both files are known to be on the same strand\n"
            "             with the same effect alleles.\n"
            "  --out      Output prefix; writes <prefix>.log plus <prefix>.h2 or\n"
            "             <prefix>.rg with the estimates at full precision, and\n"
            "             <prefix>.results for a partitioned regression.\n",
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

  LdscScores ref;
  if (LdscReadLdscores(ref_ld_arg, ref_ld_chr_split, &ref)) {
    return 1;
  }
  if (ref.ids.empty()) {
    fprintf(stderr, "Error: No LD Scores read.\n");
    return 1;
  }
  const uint32_t n_annot = ref.annot_names.size();
  LdscLog("Read reference panel LD Scores for %" PRIuPTR " variants", S_CAST(uintptr_t, ref.ids.size()));
  if (n_annot == 1) {
    LdscLog(".\n");
  } else {
    LdscLog(" in %u annotations.\n", n_annot);
  }
  LdscScores w_scores;
  if (LdscReadLdscores(w_ld_arg, w_ld_chr_split, &w_scores)) {
    return 1;
  }
  if (w_scores.annot_names.size() != 1) {
    fprintf(stderr, "Error: --w-ld/--w-ld-chr must name a single LD Score column; %" PRIuPTR " were\nfound.  The regression weights come from one sum of r^2, taken over the\nregression variants.\n", S_CAST(uintptr_t, w_scores.annot_names.size()));
    return 1;
  }
  LdscLog("Read regression weight LD Scores for %" PRIuPTR " variants.\n", S_CAST(uintptr_t, w_scores.ids.size()));
  std::unordered_map<std::string, double> w_ld_map;
  for (uintptr_t i = 0; i != w_scores.ids.size(); ++i) {
    w_ld_map.emplace(w_scores.ids[i], w_scores.l2[i]);
  }

  std::vector<double> m_vec;
  if (opts.m_override > 0.0) {
    if (n_annot != 1) {
      fprintf(stderr, "Error: --M takes one value per annotation; use the .l2.M%s files for a\npartitioned regression.\n", opts.not_m_5_50? "" : "_5_50");
      return 1;
    }
    m_vec.assign(1, opts.m_override);
  } else if (!LdscReadM(ref_ld_arg, ref_ld_chr_split, opts.not_m_5_50, &m_vec)) {
    if (m_vec.size() != n_annot) {
      fprintf(stderr, "Error: the %s files have %" PRIuPTR " entries, but there are %u LD Score\ncolumns.\n", opts.not_m_5_50? ".l2.M" : ".l2.M_5_50", S_CAST(uintptr_t, m_vec.size()), n_annot);
      return 1;
    }
    double m_sum = 0.0;
    for (uint32_t j = 0; j != n_annot; ++j) {
      m_sum += m_vec[j];
    }
    LdscLog("Read M = %g from the %s files.\n", m_sum, opts.not_m_5_50? ".l2.M" : ".l2.M_5_50");
  } else if (n_annot != 1) {
    fprintf(stderr, "Error: A partitioned regression needs the %s files, which name the\nvariant count per annotation; none were found next to --ref-ld.\n", opts.not_m_5_50? ".l2.M" : ".l2.M_5_50");
    return 1;
  } else {
    m_vec.assign(1, u31tod(ref.ids.size()));
    LdscLog("No %s file found; taking M to be the %g LD Scores read.  Pass --M if the\nLD Scores were computed from a different variant set.\n", opts.not_m_5_50? ".l2.M" : ".l2.M_5_50", m_vec[0]);
  }

  BoolErr ret;
  if (h2_fname) {
    ret = LdscEstimateH2(h2_fname, ref, ref_ld_arg, ref_ld_chr_split, w_ld_map, m_vec, &opts, out_prefix);
  } else {
    std::vector<std::string> rg_paths;
    LdscSplitComma(rg_arg, &rg_paths);
    if (rg_paths.size() < 2) {
      fprintf(stderr, "Error: --rg needs at least two summary statistic files.\n");
      return 1;
    }
    ret = LdscEstimateRg(rg_paths, ref, w_ld_map, m_vec, &opts, out_prefix);
  }
  if (fclose(g_ldsc_logfile)) {
    fprintf(stderr, "Error: Failed to write %s.\n", log_path.c_str());
    return 1;
  }
  g_ldsc_logfile = nullptr;
  return ret? 1 : 0;
}
