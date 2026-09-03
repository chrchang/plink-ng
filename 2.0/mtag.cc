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
// This implements the estimator: the method-of-moments estimate of Omega, the
// positive-semidefiniteness adjustment, and the per-variant GLS.  Sigma, the
// residual covariance from bivariate LD Score regression, is read from a file
// for now.
//
// The reference implementation (JonJala/mtag) requires Python 2.7, which
// reached end of life in January 2020.

#include <string.h>

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
  const char* sigma_fname = nullptr;
  const char* out_fname = nullptr;
  uint32_t thread_ct = 0;
  for (int argi = 1; argi < argc; ++argi) {
    const char* cur = argv[argi];
    if ((!strcmp(cur, "--zn")) && (argi + 1 < argc)) {
      z_fname = argv[++argi];
    } else if ((!strcmp(cur, "--sigma")) && (argi + 1 < argc)) {
      sigma_fname = argv[++argi];
    } else if ((!strcmp(cur, "--out")) && (argi + 1 < argc)) {
      out_fname = argv[++argi];
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
  if ((!z_fname) || (!sigma_fname) || (!out_fname)) {
    fprintf(stderr,
            "%s\n"
            "Multi-Trait Analysis of GWAS (Turley et al. 2018).\n\n"
            "Usage: mtag --zn <file> --sigma <file> --out <prefix> [--threads <n>]\n\n"
            "  --zn     Header line 'ID Z1 N1 Z2 N2 ...', one row per variant.\n"
            "  --sigma  Residual covariance matrix, one row per trait, whitespace\n"
            "           delimited.  This is the bivariate LD Score regression\n"
            "           intercept matrix.\n"
            "  --out    Output prefix; writes <prefix>.mtag.\n",
            kMtagVersion);
    return 1;
  }
  if (!thread_ct) {
    thread_ct = 1;
  }

  // Read the Z/N table.  Header determines the trait count.
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
  const uint32_t trait_ct = (token_ct - 1) / 2;

  uintptr_t capacity = 65536;
  uintptr_t variant_ct = 0;
  double* zs = S_CAST(double*, malloc(capacity * trait_ct * sizeof(double)));
  double* ns = S_CAST(double*, malloc(capacity * trait_ct * sizeof(double)));
  char** ids = S_CAST(char**, malloc(capacity * sizeof(char*)));
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
  if (MtagReadMatrix(sigma_fname, trait_ct, sigma)) {
    return 1;
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

