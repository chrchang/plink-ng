// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
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
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.

// Multi-Trait Analysis of GWAS (Turley et al., Nat Genet 50:229-237, 2018).
//
// Given GWAS summary statistics for several genetically correlated traits,
// MTAG re-estimates each trait's per-variant effect by borrowing information
// from the others.  It is a generalized inverse-variance-weighted meta-analysis
// across distinct traits rather than across studies of one trait.
//
// The whole pipeline is here: summary statistics are read and put on a common
// allele orientation, Sigma comes from LD Score regression, Omega from the
// method of moments, the per-variant GLS follows, and maxFDR bounds the cost
// of the homogeneity assumption failing.
//
// The reference implementation (JonJala/mtag) requires Python 2.7, which
// reached end of life in January 2020, shells out to a bundled copy of ldsc
// for Sigma, and loads everything through pandas.

#include <string.h>

#include <string>
#include <unordered_map>

#include "include/plink2_base.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static const char kMtagVersion[] = "mtag (plink-ng) v0.1";

// Gauss-Jordan with partial pivoting.  The matrices here are trait x trait, so
// dim is single digits and a dedicated small-matrix routine beats calling out
// to LAPACK.
BoolErr MtagInvertMatrix(uint32_t dim, double* mat, double* workbuf) {
  double* aug = workbuf;  // dim x 2*dim
  const uint32_t width = 2 * dim;
  for (uint32_t i = 0; i != dim; ++i) {
    for (uint32_t j = 0; j != dim; ++j) {
      aug[i * width + j] = mat[i * dim + j];
      aug[i * width + dim + j] = (i == j)? 1.0 : 0.0;
    }
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
    for (uint32_t j = 0; j != dim; ++j) {
      mat[i * dim + j] = aug[i * width + dim + j];
    }
  }
  return 0;
}

// Cyclic Jacobi.  Only the smallest eigenvalue is wanted, for the
// positive-semidefiniteness test, and dim is single digits.
double MtagMinEigenvalue(uint32_t dim, const double* mat, double* workbuf) {
  double* a = workbuf;
  memcpy(a, mat, dim * dim * sizeof(double));
  for (uint32_t sweep = 0; sweep != 100; ++sweep) {
    double off = 0.0;
    for (uint32_t p = 0; p != dim; ++p) {
      for (uint32_t q = p + 1; q != dim; ++q) {
        off += a[p * dim + q] * a[p * dim + q];
      }
    }
    if (off <= 1e-30) {
      break;
    }
    for (uint32_t p = 0; p != dim; ++p) {
      for (uint32_t q = p + 1; q != dim; ++q) {
        const double apq = a[p * dim + q];
        if (fabs(apq) < 1e-300) {
          continue;
        }
        const double theta = (a[q * dim + q] - a[p * dim + p]) / (2 * apq);
        const double t = (theta >= 0.0)? (1.0 / (theta + sqrt(1.0 + theta * theta))) : (-1.0 / (-theta + sqrt(1.0 + theta * theta)));
        const double c = 1.0 / sqrt(1.0 + t * t);
        const double s = t * c;
        for (uint32_t k = 0; k != dim; ++k) {
          const double akp = a[k * dim + p];
          const double akq = a[k * dim + q];
          a[k * dim + p] = c * akp - s * akq;
          a[k * dim + q] = s * akp + c * akq;
        }
        for (uint32_t k = 0; k != dim; ++k) {
          const double apk = a[p * dim + k];
          const double aqk = a[q * dim + k];
          a[p * dim + k] = c * apk - s * aqk;
          a[q * dim + k] = s * apk + c * aqk;
        }
      }
    }
  }
  double min_eig = a[0];
  for (uint32_t i = 1; i != dim; ++i) {
    if (a[i * dim + i] < min_eig) {
      min_eig = a[i * dim + i];
    }
  }
  return min_eig;
}

// Mirrors _posDef_adjustment() in the reference: clamp any off-diagonal that
// exceeds the Cauchy-Schwarz bound, then shrink all off-diagonals by 1% at a
// time until the matrix is positive semidefinite.
void MtagPosDefAdjustment(uint32_t dim, double* mat, double* workbuf, uint32_t* iter_ct_ptr) {
  *iter_ct_ptr = 0;
  if (MtagMinEigenvalue(dim, mat, workbuf) >= 0.0) {
    return;
  }
  for (uint32_t i = 0; i != dim; ++i) {
    for (uint32_t j = i; j != dim; ++j) {
      const double bound = sqrt(mat[i * dim + i] * mat[j * dim + j]);
      if (fabs(mat[i * dim + j]) > bound) {
        const double adjusted = (mat[i * dim + j] >= 0.0)? (0.99 * bound) : (-0.99 * bound);
        mat[i * dim + j] = adjusted;
        mat[j * dim + i] = adjusted;
      }
    }
  }
  uint32_t iter_ct = 0;
  while ((MtagMinEigenvalue(dim, mat, workbuf) < 0.0) && (iter_ct < 1000)) {
    for (uint32_t i = 0; i != dim; ++i) {
      for (uint32_t j = 0; j != dim; ++j) {
        if (i != j) {
          mat[i * dim + j] *= 0.99;
        }
      }
    }
    ++iter_ct;
  }
  *iter_ct_ptr = iter_ct;
}

// Method-of-moments estimate of Omega:
//   mean over variants of (z z' - Sigma) / sqrt(N (x) N)
void MtagGmmOmega(const double* zs, const double* ns, const double* sigma, uintptr_t variant_ct, uint32_t trait_ct, double* omega) {
  for (uint32_t i = 0; i != trait_ct * trait_ct; ++i) {
    omega[i] = 0.0;
  }
  for (uintptr_t vidx = 0; vidx != variant_ct; ++vidx) {
    const double* cur_z = &(zs[vidx * trait_ct]);
    const double* cur_n = &(ns[vidx * trait_ct]);
    for (uint32_t i = 0; i != trait_ct; ++i) {
      for (uint32_t j = 0; j != trait_ct; ++j) {
        const double num = cur_z[i] * cur_z[j] - sigma[i * trait_ct + j];
        omega[i * trait_ct + j] += num / sqrt(cur_n[i] * cur_n[j]);
      }
    }
  }
  const double denom = 1.0 / u31tod(1) / S_CAST(double, variant_ct);
  for (uint32_t i = 0; i != trait_ct * trait_ct; ++i) {
    omega[i] *= denom;
  }
}


// Per-variant GLS.  For trait p:
//   gamma = Omega[:,p], tau2 = Omega[p,p], yy = gamma / tau2
//   xx_m  = Omega - gamma gamma' / tau2 + Sigma_m,  Sigma_m[i][j] = Sigma[i][j] / sqrt(N_mi N_mj)
//   beta  = (yy' xx_m^-1 W_m^-1 z_m) / (yy' xx_m^-1 yy),  se = sqrt(1 / (yy' xx_m^-1 yy))
typedef struct MtagCtxStruct {
  const double* zs;
  const double* ns;
  const double* sigma;
  const double* om_min_gam;  // trait_ct x trait_ct, for the current trait
  const double* yy;
  uintptr_t variant_ct;
  uint32_t trait_ct;
  double* betas;   // stride trait_ct, column cur_trait
  double* ses;
  uint32_t cur_trait;
  uint32_t singular_ct;
} MtagCtx;

void MtagComputeRange(const MtagCtx* ctx, uintptr_t vidx_start, uintptr_t vidx_end, double* workbuf, uint32_t* singular_ct_ptr) {
  const double* zs = ctx->zs;
  const double* ns = ctx->ns;
  const double* sigma = ctx->sigma;
  const double* om_min_gam = ctx->om_min_gam;
  const double* yy = ctx->yy;
  const uint32_t trait_ct = ctx->trait_ct;
  const uint32_t cur_trait = ctx->cur_trait;
  double* betas = ctx->betas;
  double* ses = ctx->ses;
  double* xx = workbuf;
  double* inv_workbuf = &(workbuf[trait_ct * trait_ct]);
  double* row = &(inv_workbuf[2 * trait_ct * trait_ct]);
  uint32_t singular_ct = 0;
  for (uintptr_t vidx = vidx_start; vidx != vidx_end; ++vidx) {
    const double* cur_z = &(zs[vidx * trait_ct]);
    const double* cur_n = &(ns[vidx * trait_ct]);
    for (uint32_t i = 0; i != trait_ct; ++i) {
      for (uint32_t j = 0; j != trait_ct; ++j) {
        xx[i * trait_ct + j] = om_min_gam[i * trait_ct + j] + sigma[i * trait_ct + j] / sqrt(cur_n[i] * cur_n[j]);
      }
    }
    if (unlikely(MtagInvertMatrix(trait_ct, xx, inv_workbuf))) {
      betas[vidx * trait_ct + cur_trait] = 0.0 / 0.0;
      ses[vidx * trait_ct + cur_trait] = 0.0 / 0.0;
      ++singular_ct;
      continue;
    }
    // row = yy' * xx^-1
    for (uint32_t j = 0; j != trait_ct; ++j) {
      double acc = 0.0;
      for (uint32_t q = 0; q != trait_ct; ++q) {
        acc += yy[q] * xx[q * trait_ct + j];
      }
      row[j] = acc;
    }
    double denom = 0.0;
    double numer = 0.0;
    for (uint32_t j = 0; j != trait_ct; ++j) {
      denom += row[j] * yy[j];
      numer += row[j] * (cur_z[j] / sqrt(cur_n[j]));
    }
    betas[vidx * trait_ct + cur_trait] = numer / denom;
    ses[vidx * trait_ct + cur_trait] = sqrt(1.0 / denom);
  }
  *singular_ct_ptr = singular_ct;
}

THREAD_FUNC_DECL MtagThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  const uint32_t thread_ct_p1 = 1 + GetThreadCt(arg->sharedp);
  MtagCtx* ctx = S_CAST(MtagCtx*, arg->sharedp->context);
  const uint32_t trait_ct = ctx->trait_ct;
  double* workbuf = S_CAST(double*, malloc((3 * trait_ct * trait_ct + trait_ct) * sizeof(double)));
  do {
    if (workbuf) {
      const uintptr_t variant_ct = ctx->variant_ct;
      const uintptr_t start = (variant_ct * tidx) / thread_ct_p1;
      const uintptr_t stop = (variant_ct * (tidx + 1)) / thread_ct_p1;
      uint32_t singular_ct;
      MtagComputeRange(ctx, start, stop, workbuf, &singular_ct);
      if (singular_ct) {
        __atomic_fetch_add(&ctx->singular_ct, singular_ct, __ATOMIC_RELAXED);
      }
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  free(workbuf);
  THREAD_RETURN;
}

// ---------------------------------------------------------------------------
// LD Score regression, following Bulik-Sullivan et al. (2015).  Only the
// intercepts are needed here: the univariate intercepts form Sigma's diagonal
// and the bivariate ones its off-diagonals.  Standard errors would need the
// block jackknife, which MTAG does not consume, so it is not implemented.
// ---------------------------------------------------------------------------

// Weighted least squares on [ld, 1].  The reference normalizes the weights
// before calling lstsq, which cancels out of the solution.
BoolErr LdscWls2(const double* xs, const double* ys, const double* ws, uintptr_t ct, double* slope_ptr, double* intercept_ptr) {
  double sw = 0.0;
  double sx = 0.0;
  double sy = 0.0;
  double sxx = 0.0;
  double sxy = 0.0;
  for (uintptr_t idx = 0; idx != ct; ++idx) {
    const double w = ws[idx];
    const double x = xs[idx];
    const double y = ys[idx];
    sw += w;
    sx += w * x;
    sy += w * y;
    sxx += w * x * x;
    sxy += w * x * y;
  }
  const double det = sw * sxx - sx * sx;
  if (det == 0.0) {
    return 1;
  }
  *slope_ptr = (sw * sxy - sx * sy) / det;
  *intercept_ptr = (sxx * sy - sx * sxy) / det;
  return 0;
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

// Hsq.weights(): reciprocal of the conditional variance, times the
// overcounting correction 1 / w_ld.
void LdscHsqWeights(const double* ld, const double* w_ld, const double* ns, uintptr_t ct, double m_tot, double hsq, double intercept, double* dst) {
  hsq = LdscClamp(hsq, 0.0, 1.0);
  for (uintptr_t idx = 0; idx != ct; ++idx) {
    const double cur_ld = MAXV(ld[idx], 1.0);
    const double cur_w_ld = MAXV(w_ld[idx], 1.0);
    const double c = hsq * ns[idx] / m_tot;
    const double denom = intercept + c * cur_ld;
    dst[idx] = 1.0 / (2 * denom * denom * cur_w_ld);
  }
}

// Gencov.weights().
void LdscGencovWeights(const double* ld, const double* w_ld, const double* n1, const double* n2, uintptr_t ct, double m_tot, double h1, double h2, double rho_g, double intercept_gencov, double intercept_hsq1, double intercept_hsq2, double* dst) {
  h1 = LdscClamp(h1, 0.0, 1.0);
  h2 = LdscClamp(h2, 0.0, 1.0);
  rho_g = LdscClamp(rho_g, -1.0, 1.0);
  for (uintptr_t idx = 0; idx != ct; ++idx) {
    const double cur_ld = MAXV(ld[idx], 1.0);
    const double cur_w_ld = MAXV(w_ld[idx], 1.0);
    const double a = n1[idx] * h1 * cur_ld / m_tot + intercept_hsq1;
    const double b = n2[idx] * h2 * cur_ld / m_tot + intercept_hsq2;
    const double c = sqrt(n1[idx] * n2[idx]) * rho_g * cur_ld / m_tot + intercept_gencov;
    dst[idx] = 1.0 / ((a * b + c * c) * cur_w_ld);
  }
}

// The common driver: aggregate estimate, initial weights, then the two IRWLS
// reweighting passes the reference performs.
typedef enum { kLdscHsq, kLdscGencov } LdscKind;

typedef struct LdscArgsStruct {
  const double* ld;
  const double* w_ld;
  const double* ys;      // chi^2, or z1*z2
  const double* ns;      // N, or sqrt(N1 N2)
  const double* n1;      // gencov only
  const double* n2;
  uintptr_t variant_ct;
  double m_tot;
  double hsq1;           // gencov only
  double hsq2;
  double intercept_hsq1;
  double intercept_hsq2;
} LdscArgs;

BoolErr LdscRegress(const LdscArgs* args, LdscKind kind, double* weight_buf, double* xbuf, double* est_ptr, double* intercept_ptr) {
  const uintptr_t ct = args->variant_ct;
  const double* ld = args->ld;
  const double* ns = args->ns;
  const double m_tot = args->m_tot;
  const double null_intercept = (kind == kLdscHsq)? 1.0 : 0.0;

  double mean_y = 0.0;
  double mean_xn = 0.0;
  double mean_n = 0.0;
  for (uintptr_t idx = 0; idx != ct; ++idx) {
    mean_y += args->ys[idx];
    mean_xn += ld[idx] * ns[idx];
    mean_n += ns[idx];
  }
  const double ct_d = S_CAST(double, ct);
  mean_y /= ct_d;
  mean_xn /= ct_d;
  mean_n /= ct_d;
  if (mean_xn == 0.0) {
    return 1;
  }
  // aggregate(): the method-of-moments starting value.
  double tot_agg = m_tot * (mean_y - null_intercept) / mean_xn;

  if (kind == kLdscHsq) {
    LdscHsqWeights(ld, args->w_ld, ns, ct, m_tot, tot_agg, null_intercept, weight_buf);
  } else {
    LdscGencovWeights(ld, args->w_ld, args->n1, args->n2, ct, m_tot, args->hsq1, args->hsq2, tot_agg, null_intercept, args->intercept_hsq1, args->intercept_hsq2, weight_buf);
  }
  // The design matrix is N * ld / Nbar, which keeps the condition number down.
  const double nbar = mean_n;
  for (uintptr_t idx = 0; idx != ct; ++idx) {
    xbuf[idx] = ns[idx] * ld[idx] / nbar;
  }
  double slope = 0.0;
  double intercept = null_intercept;
  // The reference runs exactly two reweighting iterations.
  for (uint32_t iter = 0; iter != 3; ++iter) {
    if (unlikely(LdscWls2(xbuf, args->ys, weight_buf, ct, &slope, &intercept))) {
      return 1;
    }
    if (iter == 2) {
      break;
    }
    const double cur_est = m_tot * slope / nbar;
    if (kind == kLdscHsq) {
      // The reference writes max(x[0][1]) here, which takes the max of a
      // one-element array, i.e. the intercept itself.  No clamping.
      LdscHsqWeights(ld, args->w_ld, ns, ct, m_tot, cur_est, intercept, weight_buf);
    } else {
      LdscGencovWeights(ld, args->w_ld, args->n1, args->n2, ct, m_tot, args->hsq1, args->hsq2, cur_est, intercept, args->intercept_hsq1, args->intercept_hsq2, weight_buf);
    }
  }
  *est_ptr = m_tot * slope / nbar;
  *intercept_ptr = intercept;
  return 0;
}

// ---------------------------------------------------------------------------
// maxFDR (Supplementary Note 1.4).  MTAG assumes the variance-covariance
// matrix of effects is the same for every variant; when that fails, the false
// discovery rate rises.  maxFDR bounds how bad it can get, by maximizing the
// FDR over the space of causal-state probabilities.  The paper ships this
// diagnostic precisely because the assumption is violable, so the estimator
// should not be used without it.
// ---------------------------------------------------------------------------

double MtagNormalSf(double x) {
  return 0.5 * erfc(x * M_SQRT1_2);
}

// Inverse survival function, by bisection.  Called once per run.
double MtagNormalIsf(double p) {
  double lo = 0.0;
  double hi = 40.0;
  for (uint32_t iter = 0; iter != 200; ++iter) {
    const double mid = 0.5 * (lo + hi);
    if (MtagNormalSf(mid) > p) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

// scale_omega(): divide each entry by the prior mass of states causal for both
// traits.
BoolErr MtagScaleOmega(const double* omega, const double* probs, const unsigned char* states, uint32_t state_ct, uint32_t trait_ct, double* dst) {
  for (uint32_t p1 = 0; p1 != trait_ct; ++p1) {
    for (uint32_t p2 = 0; p2 != trait_ct; ++p2) {
      double mass = 0.0;
      for (uint32_t k = 0; k != state_ct; ++k) {
        if (states[k * trait_ct + p1] && states[k * trait_ct + p2]) {
          mass += probs[k];
        }
      }
      if (mass == 0.0) {
        return 1;
      }
      dst[p1 * trait_ct + p2] = omega[p1 * trait_ct + p2] / mass;
    }
  }
  return 0;
}

// MTAG_var_Z_jt_c(): variance of the MTAG Z statistic for trait t under
// causal state c, at a single sample-size row.
double MtagVarZ(uint32_t t, const double* omega, const double* omega_c, const double* sigma_j, uint32_t trait_ct, double* workbuf) {
  double* xx = workbuf;
  double* inv_workbuf = &(workbuf[trait_ct * trait_ct]);
  double* yy = &(inv_workbuf[2 * trait_ct * trait_ct]);
  double* row = &(yy[trait_ct]);
  double* tmp = &(row[trait_ct]);
  const double tau2 = omega[t * trait_ct + t];
  for (uint32_t i = 0; i != trait_ct; ++i) {
    const double gamma_i = omega[i * trait_ct + t];
    yy[i] = gamma_i / tau2;
    for (uint32_t j = 0; j != trait_ct; ++j) {
      xx[i * trait_ct + j] = omega[i * trait_ct + j] - gamma_i * omega[j * trait_ct + t] / tau2 + sigma_j[i * trait_ct + j];
    }
  }
  if (MtagInvertMatrix(trait_ct, xx, inv_workbuf)) {
    return -1.0;
  }
  // row = yy' xx^-1, which is also xx^-1 yy by symmetry.
  for (uint32_t j = 0; j != trait_ct; ++j) {
    double acc = 0.0;
    for (uint32_t q = 0; q != trait_ct; ++q) {
      acc += yy[q] * xx[q * trait_ct + j];
    }
    row[j] = acc;
  }
  // tmp = (Omega_c + Sigma_j) row
  for (uint32_t i = 0; i != trait_ct; ++i) {
    double acc = 0.0;
    for (uint32_t j = 0; j != trait_ct; ++j) {
      acc += (omega_c[i * trait_ct + j] + sigma_j[i * trait_ct + j]) * row[j];
    }
    tmp[i] = acc;
  }
  double numer = 0.0;
  double denom = 0.0;
  for (uint32_t i = 0; i != trait_ct; ++i) {
    numer += row[i] * tmp[i];
    denom += yy[i] * row[i];
  }
  if (denom == 0.0) {
    return -1.0;
  }
  return numer / denom;
}

// compute_fdr() for one grid point and one trait.
double MtagComputeFdr(const double* probs, uint32_t t, const double* omega, const double* sigma_j, const unsigned char* states, uint32_t state_ct, uint32_t trait_ct, double z_threshold, double* omega_tt, double* omega_c, double* workbuf) {
  if (MtagScaleOmega(omega, probs, states, state_ct, trait_ct, omega_tt)) {
    return -1.0;
  }
  if (MtagMinEigenvalue(trait_ct, omega_tt, workbuf) < 0.0) {
    return -1.0;
  }
  double power_tot = 0.0;
  double power_null = 0.0;
  for (uint32_t k = 0; k != state_ct; ++k) {
    if (probs[k] == 0.0) {
      continue;
    }
    for (uint32_t i = 0; i != trait_ct; ++i) {
      for (uint32_t j = 0; j != trait_ct; ++j) {
        omega_c[i * trait_ct + j] = (states[k * trait_ct + i] && states[k * trait_ct + j])? omega_tt[i * trait_ct + j] : 0.0;
      }
    }
    const double var = MtagVarZ(t, omega, omega_c, sigma_j, trait_ct, workbuf);
    if (var <= 0.0) {
      return -1.0;
    }
    const double power = 2 * MtagNormalSf(z_threshold / sqrt(var)) * probs[k];
    power_tot += power;
    if (!states[k * trait_ct + t]) {
      power_null += power;
    }
  }
  if (power_tot == 0.0) {
    return -1.0;
  }
  return power_null / power_tot;
}

// Every state must carry nonzero prior mass for every trait, and every trait
// pair must share some causal mass.
uint32_t MtagGridPointUsable(const double* probs, const unsigned char* states, uint32_t state_ct, uint32_t trait_ct) {
  for (uint32_t p1 = 0; p1 != trait_ct; ++p1) {
    for (uint32_t p2 = 0; p2 != trait_ct; ++p2) {
      double mass = 0.0;
      for (uint32_t k = 0; k != state_ct; ++k) {
        if (states[k * trait_ct + p1] && states[k * trait_ct + p2]) {
          mass += probs[k];
        }
      }
      if (mass == 0.0) {
        return 0;
      }
    }
  }
  return 1;
}

typedef struct MaxFdrCtxStruct {
  const double* omega;
  const double* sigma_j;
  const unsigned char* states;
  uint32_t state_ct;
  uint32_t trait_ct;
  uint32_t intervals;
  double z_threshold;
  uint64_t grid_ct;
  double* thread_maxima;  // thread_ct x trait_ct
  uint32_t thread_ct;
} MaxFdrCtx;

// simplex_walk(): the g-th lattice point, via the g-th combination of
// (state_ct - 1) positions out of (intervals + state_ct - 1).
void MtagGridPoint(uint64_t g, uint32_t state_ct, uint32_t intervals, uint32_t* comb_buf, double* probs) {
  const uint32_t num_dims = state_ct - 1;
  const uint32_t max_ = intervals + num_dims;
  // Unrank the combination in lexicographic order.
  uint64_t remaining = g;
  uint32_t prev = 0;
  for (uint32_t d = 0; d != num_dims; ++d) {
    for (uint32_t cand = (d? (comb_buf[d - 1] + 1) : 0); cand != max_; ++cand) {
      // Number of combinations completing this prefix.
      const uint32_t left = num_dims - d - 1;
      const uint32_t pool = max_ - cand - 1;
      uint64_t cnt = 1;
      for (uint32_t i = 0; i != left; ++i) {
        cnt = cnt * (pool - i) / (i + 1);
      }
      if (remaining < cnt) {
        comb_buf[d] = cand;
        prev = cand;
        break;
      }
      remaining -= cnt;
    }
  }
  (void)prev;
  int32_t last = -1;
  for (uint32_t d = 0; d != num_dims; ++d) {
    probs[d] = (S_CAST(int32_t, comb_buf[d]) - last - 1) / S_CAST(double, intervals);
    last = comb_buf[d];
  }
  probs[num_dims] = (S_CAST(int32_t, max_) - last - 1) / S_CAST(double, intervals);
}

void MtagMaxFdrRange(const MaxFdrCtx* ctx, uint64_t g_start, uint64_t g_end, double* maxima) {
  const uint32_t trait_ct = ctx->trait_ct;
  const uint32_t state_ct = ctx->state_ct;
  double* probs = S_CAST(double*, malloc(state_ct * sizeof(double)));
  uint32_t* comb_buf = S_CAST(uint32_t*, malloc(state_ct * sizeof(uint32_t)));
  double* omega_tt = S_CAST(double*, malloc(trait_ct * trait_ct * sizeof(double)));
  double* omega_c = S_CAST(double*, malloc(trait_ct * trait_ct * sizeof(double)));
  double* workbuf = S_CAST(double*, malloc((3 * trait_ct * trait_ct + 3 * trait_ct) * sizeof(double)));
  if ((!probs) || (!comb_buf) || (!omega_tt) || (!omega_c) || (!workbuf)) {
    return;
  }
  for (uint32_t t = 0; t != trait_ct; ++t) {
    maxima[t] = -1.0;
  }
  for (uint64_t g = g_start; g != g_end; ++g) {
    MtagGridPoint(g, state_ct, ctx->intervals, comb_buf, probs);
    if (!MtagGridPointUsable(probs, ctx->states, state_ct, trait_ct)) {
      continue;
    }
    for (uint32_t t = 0; t != trait_ct; ++t) {
      const double fdr = MtagComputeFdr(probs, t, ctx->omega, ctx->sigma_j, ctx->states, state_ct, trait_ct, ctx->z_threshold, omega_tt, omega_c, workbuf);
      if (fdr > maxima[t]) {
        maxima[t] = fdr;
      }
    }
  }
  free(probs); free(comb_buf); free(omega_tt); free(omega_c); free(workbuf);
}

THREAD_FUNC_DECL MtagMaxFdrThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  const uint32_t thread_ct_p1 = 1 + GetThreadCt(arg->sharedp);
  MaxFdrCtx* ctx = S_CAST(MaxFdrCtx*, arg->sharedp->context);
  do {
    const uint64_t grid_ct = ctx->grid_ct;
    MtagMaxFdrRange(ctx, (grid_ct * tidx) / thread_ct_p1, (grid_ct * (tidx + 1)) / thread_ct_p1, &(ctx->thread_maxima[tidx * ctx->trait_ct]));
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

// ---------------------------------------------------------------------------
// Summary-statistics reading, replacing the reference's separate Python
// munging step.  Each trait's file is read independently, then the traits are
// inner-joined on variant ID and put on a common allele orientation.
// ---------------------------------------------------------------------------

typedef struct SumstatRowStruct {
  double z;
  double n;
  char a1;
  char a2;
} SumstatRow;

char MtagAcgtComplement(char c) {
  switch (c) {
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  default: return '\0';
  }
}

uint32_t MtagIsPalindromic(char a1, char a2) {
  return (MtagAcgtComplement(a1) == a2);
}

// Column lookup is case-insensitive and accepts the aliases ldsc and PLINK
// emit.
uint32_t MtagMatchCol(const char* tok, uint32_t slen, const char* const* names) {
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

// Reads one trait's summary statistics.  Z is taken directly when present,
// otherwise derived from BETA/SE or from OR/SE.
BoolErr MtagReadSumstats(const char* fname, std::unordered_map<std::string, SumstatRow>* dst, double n_override) {
  static const char* kIdNames[] = {"SNP", "ID", "RSID", "MARKERNAME", "VARIANT_ID", nullptr};
  static const char* kA1Names[] = {"A1", "EFFECT_ALLELE", "ALLELE1", "EA", nullptr};
  static const char* kA2Names[] = {"A2", "OTHER_ALLELE", "ALLELE0", "ALLELE2", "NEA", nullptr};
  static const char* kZNames[] = {"Z", "ZSCORE", "Z_SCORE", "T_STAT", "STAT", nullptr};
  static const char* kBetaNames[] = {"BETA", "B", "EFFECT", "LOG_ODDS", nullptr};
  static const char* kOrNames[] = {"OR", "ODDS_RATIO", nullptr};
  static const char* kSeNames[] = {"SE", "STDERR", "STANDARD_ERROR", nullptr};
  static const char* kNNames[] = {"N", "NMISS", "OBS_CT", "N_COMPLETE_SAMPLES", "SAMPLESIZE", nullptr};
  static const char* kNCasNames[] = {"N_CAS", "NCASE", "N_CASES", nullptr};
  static const char* kNConNames[] = {"N_CON", "NCONTROL", "N_CONTROLS", nullptr};

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
  uint32_t col_beta = UINT32_MAX;
  uint32_t col_or = UINT32_MAX;
  uint32_t col_se = UINT32_MAX;
  uint32_t col_n = UINT32_MAX;
  uint32_t col_ncas = UINT32_MAX;
  uint32_t col_ncon = UINT32_MAX;
  uint32_t col_ct = 0;
  {
    const char* iter = FirstNonTspace(header);
    for (; !IsEolnKns(*iter); ++col_ct) {
      const char* token_end = CurTokenEnd(iter);
      const uint32_t slen = token_end - iter;
      if ((col_id == UINT32_MAX) && MtagMatchCol(iter, slen, kIdNames)) {
        col_id = col_ct;
      } else if ((col_a1 == UINT32_MAX) && MtagMatchCol(iter, slen, kA1Names)) {
        col_a1 = col_ct;
      } else if ((col_a2 == UINT32_MAX) && MtagMatchCol(iter, slen, kA2Names)) {
        col_a2 = col_ct;
      } else if ((col_z == UINT32_MAX) && MtagMatchCol(iter, slen, kZNames)) {
        col_z = col_ct;
      } else if ((col_beta == UINT32_MAX) && MtagMatchCol(iter, slen, kBetaNames)) {
        col_beta = col_ct;
      } else if ((col_or == UINT32_MAX) && MtagMatchCol(iter, slen, kOrNames)) {
        col_or = col_ct;
      } else if ((col_se == UINT32_MAX) && MtagMatchCol(iter, slen, kSeNames)) {
        col_se = col_ct;
      } else if ((col_n == UINT32_MAX) && MtagMatchCol(iter, slen, kNNames)) {
        col_n = col_ct;
      } else if ((col_ncas == UINT32_MAX) && MtagMatchCol(iter, slen, kNCasNames)) {
        col_ncas = col_ct;
      } else if ((col_ncon == UINT32_MAX) && MtagMatchCol(iter, slen, kNConNames)) {
        col_ncon = col_ct;
      }
      iter = FirstNonTspace(token_end);
    }
  }
  if (col_id == UINT32_MAX) {
    fprintf(stderr, "Error: %s has no variant ID column.\n", fname);
    return 1;
  }
  if ((col_a1 == UINT32_MAX) || (col_a2 == UINT32_MAX)) {
    fprintf(stderr, "Error: %s needs A1 and A2 columns; MTAG cannot align effect directions without them.\n", fname);
    return 1;
  }
  const uint32_t have_z = (col_z != UINT32_MAX);
  const uint32_t have_beta_se = (col_beta != UINT32_MAX) && (col_se != UINT32_MAX);
  const uint32_t have_or_se = (col_or != UINT32_MAX) && (col_se != UINT32_MAX);
  if (!(have_z || have_beta_se || have_or_se)) {
    fprintf(stderr, "Error: %s needs a Z column, or BETA and SE, or OR and SE.\n", fname);
    return 1;
  }
  const uint32_t have_n = (col_n != UINT32_MAX);
  const uint32_t have_cascon = (col_ncas != UINT32_MAX) && (col_ncon != UINT32_MAX);
  if ((!have_n) && (!have_cascon) && (n_override <= 0.0)) {
    fprintf(stderr, "Error: %s has no sample size column; supply one with --n.\n", fname);
    return 1;
  }
  uintptr_t line_idx = 1;
  uintptr_t kept = 0;
  uintptr_t dropped = 0;
  while (1) {
    const char* line_start = TextGet(&txs);
    if (!line_start) {
      break;
    }
    ++line_idx;
    const char* toks[64];
    uint32_t slens[64];
    const uint32_t max_tok = MINV(col_ct, 64);
    const char* iter = FirstNonTspace(line_start);
    uint32_t tok_ct = 0;
    for (; (tok_ct != max_tok) && (!IsEolnKns(*iter)); ++tok_ct) {
      const char* token_end = CurTokenEnd(iter);
      toks[tok_ct] = iter;
      slens[tok_ct] = token_end - iter;
      iter = FirstNonTspace(token_end);
    }
    // Every needed column must be present on this line.
    uint32_t needed = col_id;
    if (col_a1 > needed) { needed = col_a1; }
    if (col_a2 > needed) { needed = col_a2; }
    if (tok_ct <= needed) {
      ++dropped;
      continue;
    }
    if ((slens[col_a1] != 1) || (slens[col_a2] != 1)) {
      // Indels carry no strand information usable here.
      ++dropped;
      continue;
    }
    char a1 = toks[col_a1][0];
    char a2 = toks[col_a2][0];
    if ((a1 >= 'a') && (a1 <= 'z')) { a1 -= 32; }
    if ((a2 >= 'a') && (a2 <= 'z')) { a2 -= 32; }
    if ((!MtagAcgtComplement(a1)) || (!MtagAcgtComplement(a2)) || (a1 == a2)) {
      ++dropped;
      continue;
    }
    double z;
    if (have_z) {
      if ((col_z >= tok_ct) || (!ScanadvDouble(toks[col_z], &z))) {
        ++dropped;
        continue;
      }
    } else {
      double effect;
      double se;
      const uint32_t eff_col = have_beta_se? col_beta : col_or;
      if ((eff_col >= tok_ct) || (col_se >= tok_ct) ||
          (!ScanadvDouble(toks[eff_col], &effect)) || (!ScanadvDouble(toks[col_se], &se)) || (se <= 0.0)) {
        ++dropped;
        continue;
      }
      if (!have_beta_se) {
        if (effect <= 0.0) {
          ++dropped;
          continue;
        }
        effect = log(effect);
      }
      z = effect / se;
    }
    double n = n_override;
    if (have_n && (col_n < tok_ct)) {
      if (!ScanadvDouble(toks[col_n], &n)) {
        n = n_override;
      }
    } else if (have_cascon && (col_ncas < tok_ct) && (col_ncon < tok_ct)) {
      double ncas;
      double ncon;
      if (ScanadvDouble(toks[col_ncas], &ncas) && ScanadvDouble(toks[col_ncon], &ncon)) {
        n = ncas + ncon;
      }
    }
    if (n < 3.0) {
      ++dropped;
      continue;
    }
    SumstatRow row;
    row.z = z;
    row.n = n;
    row.a1 = a1;
    row.a2 = a2;
    dst->emplace(std::string(toks[col_id], slens[col_id]), row);
    ++kept;
  }
  CleanupTextStream(&txs, &reterr);
  printf("  %s: %" PRIuPTR " variants (%" PRIuPTR " skipped)\n", fname, kept, dropped);
  return 0;
}

// Reads a whitespace-delimited trait_ct x trait_ct matrix.
BoolErr MtagReadMatrix(const char* fname, uint32_t trait_ct, double* dst) {
  TextStream txs;
  PreinitTextStream(&txs);
  PglErr reterr = TextStreamOpen(fname, &txs);
  if (reterr) {
    fprintf(stderr, "Error: Failed to open %s.\n", fname);
    return 1;
  }
  uint32_t row_idx = 0;
  BoolErr ret = 0;
  while (1) {
    const char* line_start = TextGet(&txs);
    if (!line_start) {
      break;
    }
    if (row_idx == trait_ct) {
      fprintf(stderr, "Error: %s has more rows than traits.\n", fname);
      ret = 1;
      break;
    }
    const char* iter = line_start;
    for (uint32_t col_idx = 0; col_idx != trait_ct; ++col_idx) {
      iter = FirstNonTspace(iter);
      double dxx;
      const char* next = ScanadvDouble(iter, &dxx);
      if (!next) {
        fprintf(stderr, "Error: %s row %u has fewer than %u numeric entries.\n", fname, row_idx + 1, trait_ct);
        ret = 1;
        break;
      }
      dst[row_idx * trait_ct + col_idx] = dxx;
      iter = next;
    }
    if (ret) {
      break;
    }
    ++row_idx;
  }
  if ((!ret) && (row_idx != trait_ct)) {
    fprintf(stderr, "Error: %s has %u rows, expected %u.\n", fname, row_idx, trait_ct);
    ret = 1;
  }
  CleanupTextStream(&txs, &reterr);
  return ret;
}

// Reads an LD Score file: a header naming SNP and L2 columns, as written by
// ldsc and shipped with the reference panels.
BoolErr MtagReadLdscores(const char* fname, std::unordered_map<std::string, double>* dst) {
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
  uint32_t snp_col = UINT32_MAX;
  uint32_t l2_col = UINT32_MAX;
  {
    const char* iter = FirstNonTspace(header);
    for (uint32_t col_idx = 0; !IsEolnKns(*iter); ++col_idx) {
      const char* token_end = CurTokenEnd(iter);
      const uint32_t slen = token_end - iter;
      if ((slen == 3) && (!memcmp(iter, "SNP", 3))) {
        snp_col = col_idx;
      } else if ((slen == 2) && (!memcmp(iter, "L2", 2))) {
        l2_col = col_idx;
      } else if ((slen == 2) && (!memcmp(iter, "ID", 2)) && (snp_col == UINT32_MAX)) {
        snp_col = col_idx;
      }
      iter = FirstNonTspace(token_end);
    }
  }
  if ((snp_col == UINT32_MAX) || (l2_col == UINT32_MAX)) {
    fprintf(stderr, "Error: %s must have SNP (or ID) and L2 columns.\n", fname);
    return 1;
  }
  const uint32_t max_col = MAXV(snp_col, l2_col);
  while (1) {
    const char* line_start = TextGet(&txs);
    if (!line_start) {
      break;
    }
    const char* iter = FirstNonTspace(line_start);
    const char* id_start = nullptr;
    uint32_t id_slen = 0;
    double l2 = 0.0;
    uint32_t ok = 0;
    for (uint32_t col_idx = 0; col_idx <= max_col; ++col_idx) {
      if (IsEolnKns(*iter)) {
        break;
      }
      const char* token_end = CurTokenEnd(iter);
      if (col_idx == snp_col) {
        id_start = iter;
        id_slen = token_end - iter;
      }
      if (col_idx == l2_col) {
        if (!ScanadvDouble(iter, &l2)) {
          break;
        }
        ok = 1;
      }
      iter = FirstNonTspace(token_end);
    }
    if (ok && id_start) {
      dst->emplace(std::string(id_start, id_slen), l2);
    }
  }
  CleanupTextStream(&txs, &reterr);
  return 0;
}

// Builds Sigma from LD Score regression: univariate intercepts on the
// diagonal, bivariate intercepts off it.  This is what the reference obtains
// by shelling out to its bundled copy of ldsc.
BoolErr MtagEstimateSigma(const double* zs, const double* ns, const double* ld, uintptr_t variant_ct, uint32_t trait_ct, double m_tot, double* sigma) {
  double* ys = S_CAST(double*, malloc(variant_ct * sizeof(double)));
  double* nbuf = S_CAST(double*, malloc(variant_ct * sizeof(double)));
  double* n1buf = S_CAST(double*, malloc(variant_ct * sizeof(double)));
  double* n2buf = S_CAST(double*, malloc(variant_ct * sizeof(double)));
  double* wbuf = S_CAST(double*, malloc(variant_ct * sizeof(double)));
  double* xbuf = S_CAST(double*, malloc(variant_ct * sizeof(double)));
  double* hsqs = S_CAST(double*, malloc(trait_ct * sizeof(double)));
  if ((!ys) || (!nbuf) || (!n1buf) || (!n2buf) || (!wbuf) || (!xbuf) || (!hsqs)) {
    fprintf(stderr, "Error: out of memory.\n");
    return 1;
  }
  LdscArgs args;
  args.ld = ld;
  args.w_ld = ld;
  args.variant_ct = variant_ct;
  args.m_tot = m_tot;
  printf("LD Score regression on %" PRIuPTR " variants (M = %g):\n", variant_ct, m_tot);
  for (uint32_t trait_idx = 0; trait_idx != trait_ct; ++trait_idx) {
    for (uintptr_t vidx = 0; vidx != variant_ct; ++vidx) {
      const double z = zs[vidx * trait_ct + trait_idx];
      ys[vidx] = z * z;
      nbuf[vidx] = ns[vidx * trait_ct + trait_idx];
    }
    args.ys = ys;
    args.ns = nbuf;
    double hsq;
    double intercept;
    if (unlikely(LdscRegress(&args, kLdscHsq, wbuf, xbuf, &hsq, &intercept))) {
      fprintf(stderr, "Error: univariate LD Score regression failed for trait %u.\n", trait_idx + 1);
      return 1;
    }
    hsqs[trait_idx] = hsq;
    sigma[trait_idx * trait_ct + trait_idx] = intercept;
    printf("  trait %u: h2 = %.10g, intercept = %.10g\n", trait_idx + 1, hsq, intercept);
  }
  for (uint32_t i = 0; i != trait_ct; ++i) {
    for (uint32_t j = i + 1; j != trait_ct; ++j) {
      for (uintptr_t vidx = 0; vidx != variant_ct; ++vidx) {
        const double z1 = zs[vidx * trait_ct + i];
        const double z2 = zs[vidx * trait_ct + j];
        ys[vidx] = z1 * z2;
        n1buf[vidx] = ns[vidx * trait_ct + i];
        n2buf[vidx] = ns[vidx * trait_ct + j];
        nbuf[vidx] = sqrt(n1buf[vidx] * n2buf[vidx]);
      }
      args.ys = ys;
      args.ns = nbuf;
      args.n1 = n1buf;
      args.n2 = n2buf;
      args.hsq1 = hsqs[i];
      args.hsq2 = hsqs[j];
      args.intercept_hsq1 = sigma[i * trait_ct + i];
      args.intercept_hsq2 = sigma[j * trait_ct + j];
      double gencov;
      double intercept;
      if (unlikely(LdscRegress(&args, kLdscGencov, wbuf, xbuf, &gencov, &intercept))) {
        fprintf(stderr, "Error: bivariate LD Score regression failed for traits %u and %u.\n", i + 1, j + 1);
        return 1;
      }
      sigma[i * trait_ct + j] = intercept;
      sigma[j * trait_ct + i] = intercept;
      printf("  traits %u,%u: gencov = %.10g, intercept = %.10g\n", i + 1, j + 1, gencov, intercept);
    }
  }
  free(ys); free(nbuf); free(n1buf); free(n2buf); free(wbuf); free(xbuf); free(hsqs);
  return 0;
}

void MtagPrintMatrix(const char* label, uint32_t trait_ct, const double* mat) {
  printf("%s:\n", label);
  for (uint32_t i = 0; i != trait_ct; ++i) {
    fputs(" ", stdout);
    for (uint32_t j = 0; j != trait_ct; ++j) {
      printf(" %g", mat[i * trait_ct + j]);
    }
    fputs("\n", stdout);
  }
}

#ifdef __cplusplus
}  // namespace plink2
#endif

int main(int argc, char** argv) {
  using namespace plink2;
  const char* z_fname = nullptr;
  const char* sumstats_arg = nullptr;
  double n_override = 0.0;
  uint32_t keep_palindromic = 0;
  const char* sigma_fname = nullptr;
  const char* ldscore_fname = nullptr;
  const char* out_fname = nullptr;
  double m_arg = 0.0;
  double p_sig = 5.0e-8;
  uint32_t intervals = 10;
  uint32_t skip_maxfdr = 0;
  uint32_t thread_ct = 0;
  for (int argi = 1; argi < argc; ++argi) {
    const char* cur = argv[argi];
    if ((!strcmp(cur, "--zn")) && (argi + 1 < argc)) {
      z_fname = argv[++argi];
    } else if ((!strcmp(cur, "--sumstats")) && (argi + 1 < argc)) {
      sumstats_arg = argv[++argi];
    } else if ((!strcmp(cur, "--n")) && (argi + 1 < argc)) {
      if ((!ScanadvDouble(argv[++argi], &n_override)) || (n_override < 3.0)) {
        fprintf(stderr, "Error: --n must be at least 3.\n");
        return 1;
      }
    } else if (!strcmp(cur, "--keep-palindromic")) {
      keep_palindromic = 1;
    } else if ((!strcmp(cur, "--sigma")) && (argi + 1 < argc)) {
      sigma_fname = argv[++argi];
    } else if ((!strcmp(cur, "--ldscore")) && (argi + 1 < argc)) {
      ldscore_fname = argv[++argi];
    } else if ((!strcmp(cur, "--m")) && (argi + 1 < argc)) {
      const char* m_str = argv[++argi];
      if ((!ScanadvDouble(m_str, &m_arg)) || (m_arg <= 0.0)) {
        fprintf(stderr, "Error: --m argument '%s' is not a positive number.\n", m_str);
        return 1;
      }
    } else if ((!strcmp(cur, "--out")) && (argi + 1 < argc)) {
      out_fname = argv[++argi];
    } else if ((!strcmp(cur, "--intervals")) && (argi + 1 < argc)) {
      intervals = atoi(argv[++argi]);
      if (intervals < 2) {
        fprintf(stderr, "Error: --intervals must be at least 2.\n");
        return 1;
      }
    } else if ((!strcmp(cur, "--p-sig")) && (argi + 1 < argc)) {
      if ((!ScanadvDouble(argv[++argi], &p_sig)) || (p_sig <= 0.0) || (p_sig >= 1.0)) {
        fprintf(stderr, "Error: --p-sig must be in (0, 1).\n");
        return 1;
      }
    } else if (!strcmp(cur, "--no-maxfdr")) {
      skip_maxfdr = 1;
    } else if ((!strcmp(cur, "--threads")) && (argi + 1 < argc)) {
      thread_ct = atoi(argv[++argi]);
    } else if (!strcmp(cur, "--version")) {
      printf("%s\n", kMtagVersion);
      return 0;
    } else {
      fprintf(stderr, "Error: unrecognized argument '%s'.\n", cur);
      return 1;
    }
  }
  if ((!(z_fname || sumstats_arg)) || (!out_fname) || (!(sigma_fname || ldscore_fname))) {
    fprintf(stderr,
            "%s\n"
            "Multi-Trait Analysis of GWAS (Turley et al. 2018).\n\n"
            "Usage: mtag --sumstats <file1,file2,...> --ldscore <file>\n"
            "            --out <prefix> [--m <count>] [--threads <n>]\n\n"
            "  --sumstats One GWAS summary statistics file per trait, comma\n"
            "             separated.  Columns are detected case-insensitively:\n"
            "             an ID (SNP/ID/MarkerName), A1 and A2, an effect\n"
            "             (Z, or BETA and SE, or OR and SE), and a sample size\n"
            "             (N, or N_CAS and N_CON).  The traits are inner-joined\n"
            "             on ID and put on a common allele orientation; effect\n"
            "             directions are flipped where A1/A2 are swapped, and\n"
            "             strand-ambiguous A/T and C/G variants are dropped,\n"
            "             since their orientation cannot be recovered.\n"
            "  --n        Sample size to use where a file has no N column.\n"
            "  --keep-palindromic   Keep strand-ambiguous variants.  Only safe\n"
            "             if every input is known to be on the same strand.\n"
            "  --zn       Alternative pre-merged input, with a header line\n"
            "             'ID Z1 N1 Z2 N2 ...' and one row per variant.\n"
            "  --ldscore  LD Scores, with SNP and L2 columns.  Sigma is then\n"
            "             estimated by LD Score regression: univariate intercepts\n"
            "             on the diagonal, bivariate intercepts off it.\n"
            "  --sigma    Supply the residual covariance matrix directly instead,\n"
            "             one row per trait.\n"
            "  --m        Number of variants the LD Scores were estimated from;\n"
            "             defaults to the number read from --ldscore.\n"
            "  --out      Output prefix; writes <prefix>.mtag.\n"
            "  --intervals  maxFDR grid resolution per dimension (default 10).\n"
            "  --p-sig      Significance threshold for maxFDR (default 5e-8).\n"
            "  --no-maxfdr  Skip the maxFDR diagnostic.  Not recommended: MTAG\n"
            "               assumes one effect covariance matrix for every\n"
            "               variant, and maxFDR is what bounds the cost of that\n"
            "               assumption being wrong.\n",
            kMtagVersion);
    return 1;
  }
  if (!thread_ct) {
    thread_ct = 1;
  }

  uint32_t trait_ct = 0;
  uintptr_t variant_ct = 0;
  double* zs = nullptr;
  double* ns = nullptr;
  char** ids = nullptr;

  if (sumstats_arg) {
    // Split the comma-separated file list.
    char* arg_copy = strdup(sumstats_arg);
    const char* fnames[64];
    uint32_t fname_ct = 0;
    for (char* iter = arg_copy; ; ) {
      char* comma = strchr(iter, ',');
      if (comma) {
        *comma = '\0';
      }
      if (fname_ct == 64) {
        fprintf(stderr, "Error: at most 64 traits.\n");
        return 1;
      }
      fnames[fname_ct++] = iter;
      if (!comma) {
        break;
      }
      iter = comma + 1;
    }
    if (fname_ct < 2) {
      fprintf(stderr, "Error: MTAG needs at least two traits.\n");
      return 1;
    }
    trait_ct = fname_ct;
    printf("%s\nReading %u summary statistics files:\n", kMtagVersion, trait_ct);
    std::unordered_map<std::string, SumstatRow>* tables = new std::unordered_map<std::string, SumstatRow>[trait_ct];
    for (uint32_t trait_idx = 0; trait_idx != trait_ct; ++trait_idx) {
      if (MtagReadSumstats(fnames[trait_idx], &(tables[trait_idx]), n_override)) {
        return 1;
      }
    }
    // Inner join on ID, using the first trait's alleles as the orientation.
    const uintptr_t max_ct = tables[0].size();
    zs = S_CAST(double*, malloc(max_ct * trait_ct * sizeof(double)));
    ns = S_CAST(double*, malloc(max_ct * trait_ct * sizeof(double)));
    ids = S_CAST(char**, malloc(max_ct * sizeof(char*)));
    if ((!zs) || (!ns) || (!ids)) {
      fprintf(stderr, "Error: out of memory.\n");
      return 1;
    }
    uintptr_t palindromic_ct = 0;
    uintptr_t mismatch_ct = 0;
    for (const auto& entry : tables[0]) {
      const SumstatRow& ref_row = entry.second;
      if ((!keep_palindromic) && MtagIsPalindromic(ref_row.a1, ref_row.a2)) {
        ++palindromic_ct;
        continue;
      }
      double buf_z[64];
      double buf_n[64];
      uint32_t ok = 1;
      for (uint32_t trait_idx = 0; trait_idx != trait_ct; ++trait_idx) {
        const auto match = tables[trait_idx].find(entry.first);
        if (match == tables[trait_idx].end()) {
          ok = 0;
          break;
        }
        const SumstatRow& row = match->second;
        double sign;
        if ((row.a1 == ref_row.a1) && (row.a2 == ref_row.a2)) {
          sign = 1.0;
        } else if ((row.a1 == ref_row.a2) && (row.a2 == ref_row.a1)) {
          sign = -1.0;
        } else if ((MtagAcgtComplement(row.a1) == ref_row.a1) && (MtagAcgtComplement(row.a2) == ref_row.a2)) {
          sign = 1.0;
        } else if ((MtagAcgtComplement(row.a1) == ref_row.a2) && (MtagAcgtComplement(row.a2) == ref_row.a1)) {
          sign = -1.0;
        } else {
          ok = 0;
          ++mismatch_ct;
          break;
        }
        buf_z[trait_idx] = sign * row.z;
        buf_n[trait_idx] = row.n;
      }
      if (!ok) {
        continue;
      }
      char* id_copy = strdup(entry.first.c_str());
      ids[variant_ct] = id_copy;
      for (uint32_t trait_idx = 0; trait_idx != trait_ct; ++trait_idx) {
        zs[variant_ct * trait_ct + trait_idx] = buf_z[trait_idx];
        ns[variant_ct * trait_ct + trait_idx] = buf_n[trait_idx];
      }
      ++variant_ct;
    }
    delete[] tables;
    free(arg_copy);
    if (palindromic_ct) {
      printf("%" PRIuPTR " strand-ambiguous variants dropped.\n", palindromic_ct);
    }
    if (mismatch_ct) {
      printf("%" PRIuPTR " variants dropped for irreconcilable alleles.\n", mismatch_ct);
    }
    if (variant_ct < 100) {
      fprintf(stderr, "Error: only %" PRIuPTR " variants are shared across all traits.\n", variant_ct);
      return 1;
    }
    printf("%" PRIuPTR " variants x %u traits after merging.\n", variant_ct, trait_ct);
    goto mtag_loaded;
  }

  // Read the pre-merged Z/N table.  Header determines the trait count.
  {
  TextStream txs;
  PreinitTextStream(&txs);
  PglErr reterr = TextStreamOpen(z_fname, &txs);
  if (reterr) {
    fprintf(stderr, "Error: Failed to open %s.\n", z_fname);
    return 1;
  }
  const char* header = TextGet(&txs);
  if (!header) {
    fprintf(stderr, "Error: %s is empty.\n", z_fname);
    return 1;
  }
  uint32_t token_ct = 0;
  {
    const char* iter = FirstNonTspace(header);
    while (!IsEolnKns(*iter)) {
      ++token_ct;
      iter = FirstNonTspace(CurTokenEnd(iter));
    }
  }
  if ((token_ct < 3) || ((token_ct - 1) % 2)) {
    fprintf(stderr, "Error: %s header must be 'ID Z1 N1 Z2 N2 ...'.\n", z_fname);
    return 1;
  }
  trait_ct = (token_ct - 1) / 2;

  uintptr_t capacity = 65536;
  zs = S_CAST(double*, malloc(capacity * trait_ct * sizeof(double)));
  ns = S_CAST(double*, malloc(capacity * trait_ct * sizeof(double)));
  ids = S_CAST(char**, malloc(capacity * sizeof(char*)));
  if ((!zs) || (!ns) || (!ids)) {
    fprintf(stderr, "Error: out of memory.\n");
    return 1;
  }
  while (1) {
    const char* line_start = TextGet(&txs);
    if (!line_start) {
      break;
    }
    if (variant_ct == capacity) {
      capacity *= 2;
      zs = S_CAST(double*, realloc(zs, capacity * trait_ct * sizeof(double)));
      ns = S_CAST(double*, realloc(ns, capacity * trait_ct * sizeof(double)));
      ids = S_CAST(char**, realloc(ids, capacity * sizeof(char*)));
      if ((!zs) || (!ns) || (!ids)) {
        fprintf(stderr, "Error: out of memory.\n");
        return 1;
      }
    }
    const char* iter = FirstNonTspace(line_start);
    const char* id_end = CurTokenEnd(iter);
    const uint32_t id_slen = id_end - iter;
    char* id_copy = S_CAST(char*, malloc(id_slen + 1));
    memcpy(id_copy, iter, id_slen);
    id_copy[id_slen] = '\0';
    ids[variant_ct] = id_copy;
    iter = id_end;
    uint32_t parse_fail = 0;
    for (uint32_t trait_idx = 0; trait_idx != trait_ct; ++trait_idx) {
      double z_val;
      double n_val;
      iter = FirstNonTspace(iter);
      const char* next = ScanadvDouble(iter, &z_val);
      if (!next) {
        parse_fail = 1;
        break;
      }
      iter = FirstNonTspace(next);
      next = ScanadvDouble(iter, &n_val);
      if ((!next) || (n_val < 3.0)) {
        parse_fail = 1;
        break;
      }
      iter = next;
      zs[variant_ct * trait_ct + trait_idx] = z_val;
      ns[variant_ct * trait_ct + trait_idx] = n_val;
    }
    if (parse_fail) {
      free(id_copy);
      continue;
    }
    ++variant_ct;
  }
  CleanupTextStream(&txs, &reterr);
  if (!variant_ct) {
    fprintf(stderr, "Error: no usable rows in %s.\n", z_fname);
    return 1;
  }
  printf("%s\n%" PRIuPTR " variants x %u traits loaded from %s.\n", kMtagVersion, variant_ct, trait_ct, z_fname);
  }
 mtag_loaded:

  double* sigma = S_CAST(double*, malloc(trait_ct * trait_ct * sizeof(double)));
  double* omega = S_CAST(double*, malloc(trait_ct * trait_ct * sizeof(double)));
  double* om_min_gam = S_CAST(double*, malloc(trait_ct * trait_ct * sizeof(double)));
  double* yy = S_CAST(double*, malloc(trait_ct * sizeof(double)));
  double* eig_buf = S_CAST(double*, malloc(trait_ct * trait_ct * sizeof(double)));
  double* betas = S_CAST(double*, malloc(variant_ct * trait_ct * sizeof(double)));
  double* ses = S_CAST(double*, malloc(variant_ct * trait_ct * sizeof(double)));
  if ((!sigma) || (!omega) || (!om_min_gam) || (!yy) || (!eig_buf) || (!betas) || (!ses)) {
    fprintf(stderr, "Error: out of memory.\n");
    return 1;
  }
  if (sigma_fname) {
    if (MtagReadMatrix(sigma_fname, trait_ct, sigma)) {
      return 1;
    }
  } else {
    std::unordered_map<std::string, double> ldmap;
    if (MtagReadLdscores(ldscore_fname, &ldmap)) {
      return 1;
    }
    if (ldmap.empty()) {
      fprintf(stderr, "Error: no LD Scores read from %s.\n", ldscore_fname);
      return 1;
    }
    // Keep only variants with an LD Score, compacting in place.
    double* ld = S_CAST(double*, malloc(variant_ct * sizeof(double)));
    if (!ld) {
      fprintf(stderr, "Error: out of memory.\n");
      return 1;
    }
    uintptr_t kept_ct = 0;
    for (uintptr_t vidx = 0; vidx != variant_ct; ++vidx) {
      const auto match = ldmap.find(std::string(ids[vidx]));
      if (match == ldmap.end()) {
        free(ids[vidx]);
        continue;
      }
      ld[kept_ct] = match->second;
      ids[kept_ct] = ids[vidx];
      memmove(&(zs[kept_ct * trait_ct]), &(zs[vidx * trait_ct]), trait_ct * sizeof(double));
      memmove(&(ns[kept_ct * trait_ct]), &(ns[vidx * trait_ct]), trait_ct * sizeof(double));
      ++kept_ct;
    }
    if (kept_ct < 100) {
      fprintf(stderr, "Error: only %" PRIuPTR " variants have LD Scores; too few to regress on.\n", kept_ct);
      return 1;
    }
    if (kept_ct != variant_ct) {
      printf("%" PRIuPTR " of %" PRIuPTR " variants have an LD Score.\n", kept_ct, variant_ct);
      variant_ct = kept_ct;
    }
    const double m_tot = (m_arg > 0.0)? m_arg : S_CAST(double, ldmap.size());
    if (MtagEstimateSigma(zs, ns, ld, variant_ct, trait_ct, m_tot, sigma)) {
      return 1;
    }
    free(ld);
  }
  MtagPrintMatrix("Sigma", trait_ct, sigma);

  MtagGmmOmega(zs, ns, sigma, variant_ct, trait_ct, omega);
  uint32_t psd_iter_ct;
  MtagPosDefAdjustment(trait_ct, omega, eig_buf, &psd_iter_ct);
  if (psd_iter_ct) {
    printf("Omega required %u positive-semidefiniteness iterations.\n", psd_iter_ct);
  }
  MtagPrintMatrix("Omega", trait_ct, omega);

  MtagCtx ctx;
  ctx.zs = zs;
  ctx.ns = ns;
  ctx.sigma = sigma;
  ctx.om_min_gam = om_min_gam;
  ctx.yy = yy;
  ctx.variant_ct = variant_ct;
  ctx.trait_ct = trait_ct;
  ctx.betas = betas;
  ctx.ses = ses;
  ctx.singular_ct = 0;

  ThreadGroup tg;
  PreinitThreads(&tg);
  uint32_t threads_active = 0;
  if (thread_ct >= 2) {
    if (!SetThreadCt(thread_ct - 1, &tg)) {
      SetThreadFuncAndData(MtagThread, &ctx, &tg);
      threads_active = 1;
    }
  }
  double* main_workbuf = S_CAST(double*, malloc((3 * trait_ct * trait_ct + trait_ct) * sizeof(double)));
  if (!main_workbuf) {
    fprintf(stderr, "Error: out of memory.\n");
    return 1;
  }
  for (uint32_t trait_idx = 0; trait_idx != trait_ct; ++trait_idx) {
    const double tau2 = omega[trait_idx * trait_ct + trait_idx];
    if (tau2 <= 0.0) {
      fprintf(stderr, "Error: Omega[%u,%u] is not positive; trait %u carries no estimable signal.\n", trait_idx, trait_idx, trait_idx + 1);
      return 1;
    }
    for (uint32_t i = 0; i != trait_ct; ++i) {
      const double gamma_i = omega[i * trait_ct + trait_idx];
      yy[i] = gamma_i / tau2;
      for (uint32_t j = 0; j != trait_ct; ++j) {
        om_min_gam[i * trait_ct + j] = omega[i * trait_ct + j] - gamma_i * omega[j * trait_ct + trait_idx] / tau2;
      }
    }
    ctx.cur_trait = trait_idx;
    if (threads_active) {
      if (trait_idx + 1 == trait_ct) {
        DeclareLastThreadBlock(&tg);
      }
      if (SpawnThreads(&tg)) {
        fprintf(stderr, "Error: thread creation failed.\n");
        return 1;
      }
      const uintptr_t start = (variant_ct * (thread_ct - 1)) / thread_ct;
      uint32_t singular_ct;
      MtagComputeRange(&ctx, start, variant_ct, main_workbuf, &singular_ct);
      ctx.singular_ct += singular_ct;
      JoinThreads(&tg);
    } else {
      uint32_t singular_ct;
      MtagComputeRange(&ctx, 0, variant_ct, main_workbuf, &singular_ct);
      ctx.singular_ct += singular_ct;
    }
  }
  CleanupThreads(&tg);
  if (ctx.singular_ct) {
    fprintf(stderr, "Warning: %u (variant, trait) pairs had a singular weight matrix; reported as NA.\n", ctx.singular_ct);
  }

  if (!skip_maxfdr) {
    const uint32_t state_ct = 1U << trait_ct;
    unsigned char* states = S_CAST(unsigned char*, malloc(state_ct * trait_ct));
    if (!states) {
      fprintf(stderr, "Error: out of memory.\n");
      return 1;
    }
    // itertools.product([False, True], repeat=T): the last trait varies fastest.
    for (uint32_t k = 0; k != state_ct; ++k) {
      for (uint32_t t = 0; t != trait_ct; ++t) {
        states[k * trait_ct + t] = (k >> (trait_ct - 1 - t)) & 1;
      }
    }
    // The reference's recommended --n_approx: one sample-size row, the mean.
    double* mean_n = S_CAST(double*, malloc(trait_ct * sizeof(double)));
    double* sigma_j = S_CAST(double*, malloc(trait_ct * trait_ct * sizeof(double)));
    if ((!mean_n) || (!sigma_j)) {
      fprintf(stderr, "Error: out of memory.\n");
      return 1;
    }
    for (uint32_t t = 0; t != trait_ct; ++t) {
      double acc = 0.0;
      for (uintptr_t vidx = 0; vidx != variant_ct; ++vidx) {
        acc += ns[vidx * trait_ct + t];
      }
      mean_n[t] = acc / S_CAST(double, variant_ct);
    }
    for (uint32_t i = 0; i != trait_ct; ++i) {
      for (uint32_t j = 0; j != trait_ct; ++j) {
        sigma_j[i * trait_ct + j] = sigma[i * trait_ct + j] / sqrt(mean_n[i] * mean_n[j]);
      }
    }
    // C(intervals + state_ct - 1, state_ct - 1) lattice points.
    uint64_t grid_ct = 1;
    {
      const uint32_t num_dims = state_ct - 1;
      const uint32_t pool = intervals + num_dims;
      for (uint32_t i = 0; i != num_dims; ++i) {
        grid_ct = grid_ct * (pool - i) / (i + 1);
        if (grid_ct > 200000000LLU) {
          break;
        }
      }
    }
    if (grid_ct > 200000000LLU) {
      fprintf(stderr, "Warning: the maxFDR grid would have more than 2e8 points for %u traits at\n--intervals %u; skipping.  Lower --intervals, or pass --no-maxfdr knowing\nthe estimator is then unbounded.\n", trait_ct, intervals);
    } else {
      MaxFdrCtx fctx;
      fctx.omega = omega;
      fctx.sigma_j = sigma_j;
      fctx.states = states;
      fctx.state_ct = state_ct;
      fctx.trait_ct = trait_ct;
      fctx.intervals = intervals;
      fctx.z_threshold = MtagNormalIsf(p_sig / 2.0);
      fctx.grid_ct = grid_ct;
      fctx.thread_ct = thread_ct;
      double* thread_maxima = S_CAST(double*, malloc(thread_ct * trait_ct * sizeof(double)));
      if (!thread_maxima) {
        fprintf(stderr, "Error: out of memory.\n");
        return 1;
      }
      for (uint32_t i = 0; i != thread_ct * trait_ct; ++i) {
        thread_maxima[i] = -1.0;
      }
      fctx.thread_maxima = thread_maxima;
      printf("maxFDR over %" PRIu64 " grid points (p = %g)...\n", grid_ct, p_sig);
      ThreadGroup ftg;
      PreinitThreads(&ftg);
      if ((thread_ct >= 2) && (!SetThreadCt(thread_ct - 1, &ftg))) {
        SetThreadFuncAndData(MtagMaxFdrThread, &fctx, &ftg);
        DeclareLastThreadBlock(&ftg);
        if (SpawnThreads(&ftg)) {
          fprintf(stderr, "Error: thread creation failed.\n");
          return 1;
        }
        MtagMaxFdrRange(&fctx, (grid_ct * (thread_ct - 1)) / thread_ct, grid_ct, &(thread_maxima[(thread_ct - 1) * trait_ct]));
        JoinThreads(&ftg);
      } else {
        MtagMaxFdrRange(&fctx, 0, grid_ct, thread_maxima);
      }
      CleanupThreads(&ftg);
      for (uint32_t t = 0; t != trait_ct; ++t) {
        double best = -1.0;
        for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
          if (thread_maxima[tidx * trait_ct + t] > best) {
            best = thread_maxima[tidx * trait_ct + t];
          }
        }
        if (best < 0.0) {
          printf("  trait %u: maxFDR undefined (no feasible grid point).\n", t + 1);
        } else {
          printf("  trait %u: maxFDR = %.10g\n", t + 1, best);
        }
      }
      free(thread_maxima);
    }
    free(states); free(mean_n); free(sigma_j);
  }

  char* outname = S_CAST(char*, malloc(strlen(out_fname) + 16));
  sprintf(outname, "%s.mtag", out_fname);
  FILE* outfile = fopen(outname, "w");
  if (!outfile) {
    fprintf(stderr, "Error: Failed to open %s for writing.\n", outname);
    return 1;
  }
  fputs("ID", outfile);
  for (uint32_t trait_idx = 0; trait_idx != trait_ct; ++trait_idx) {
    fprintf(outfile, "\tBETA%u\tSE%u\tZ%u", trait_idx + 1, trait_idx + 1, trait_idx + 1);
  }
  fputs("\n", outfile);
  for (uintptr_t vidx = 0; vidx != variant_ct; ++vidx) {
    fputs(ids[vidx], outfile);
    for (uint32_t trait_idx = 0; trait_idx != trait_ct; ++trait_idx) {
      const double beta = betas[vidx * trait_ct + trait_idx];
      const double se = ses[vidx * trait_ct + trait_idx];
      fprintf(outfile, "\t%.12g\t%.12g\t%.12g", beta, se, beta / se);
    }
    fputs("\n", outfile);
  }
  if (fclose(outfile)) {
    fprintf(stderr, "Error: write failure on %s.\n", outname);
    return 1;
  }
  printf("MTAG results written to %s .\n", outname);
  return 0;
}

