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

#include "plink2_epistasis.h"

#include <math.h>
#include <string.h>

#include "include/pgenlib_misc.h"
#include "include/plink2_bits.h"
#include "include/plink2_float.h"
#include "include/plink2_simd.h"
#include "include/plink2_stats.h"
#include "include/plink2_string.h"
#include "include/plink2_thread.h"
#include "plink2_cmdline.h"
#include "plink2_compress_stream.h"
#include "plink2_glm_logistic.h"
#include "plink2_matrix.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitEpi(EpiInfo* epi_ip) {
  epi_ip->flags = kfEpi0;
  epi_ip->ln_epi1 = 1.0;
  epi_ip->ln_epi2 = -2 * kLn10;
}

// --epistasis-boost: the two-stage test of Wan X et al. (2010) BOOST: A fast
// approach to detecting gene-gene interactions in genome-wide case-control
// studies.  Every pair is first scored with the Kirkwood superposition
// approximation, which is closed-form; only the pairs clearing the --epi1
// threshold are then fit properly, by iterative proportional fitting of the
// homogeneous association model (all three two-way interactions, no three-way
// term).  The reported statistic is the fitted one, so --epi1 is a screening
// threshold here rather than a report filter, and a looser --epi1 changes
// which pairs are fit rather than only which are printed.
//
// counts is the 2x3x3 table, laid out as [group * 9 + geno1 * 3 + geno2] with
// group 0 the cases.  An empty genotype row or column costs two degrees of
// freedom instead of dropping the pair, so df is 4, 2 or 1; a second empty row
// (or column) leaves nothing to test.

// p_bc[3 * group + geno2] = P(geno2 | group).
static void FepiBoostPBc(const uint32_t* counts, const double* recip_cache, double* p_bc) {
  for (uint32_t group_idx = 0; group_idx != 2; ++group_idx) {
    const uint32_t* cur_counts = &(counts[group_idx * 9]);
    uint32_t col_cts[3];
    for (uint32_t geno2 = 0; geno2 != 3; ++geno2) {
      col_cts[geno2] = cur_counts[geno2] + cur_counts[geno2 + 3] + cur_counts[geno2 + 6];
    }
    const double tot_recip = recip_cache[col_cts[0] + col_cts[1] + col_cts[2]];
    for (uint32_t geno2 = 0; geno2 != 3; ++geno2) {
      p_bc[group_idx * 3 + geno2] = u31tod(col_cts[geno2]) * tot_recip;
    }
  }
}

// p_ca[2 * geno1 + group] = P(group | geno1).  Returns the number of empty
// genotype rows.
static uint32_t FepiBoostPCa(const uint32_t* counts, const double* recip_cache, double* p_ca) {
  uint32_t empty_row_ct = 0;
  for (uint32_t geno1 = 0; geno1 != 3; ++geno1) {
    const uint32_t* case_row = &(counts[geno1 * 3]);
    const uint32_t* ctrl_row = &(counts[9 + geno1 * 3]);
    const uint32_t case_row_ct = case_row[0] + case_row[1] + case_row[2];
    const uint32_t ctrl_row_ct = ctrl_row[0] + ctrl_row[1] + ctrl_row[2];
    const uint32_t row_ct = case_row_ct + ctrl_row_ct;
    empty_row_ct += (row_ct == 0);
    const double tot_recip = recip_cache[row_ct];
    p_ca[geno1 * 2] = u31tod(case_row_ct) * tot_recip;
    p_ca[geno1 * 2 + 1] = u31tod(ctrl_row_ct) * tot_recip;
  }
  return empty_row_ct;
}

// Returns 1 when the table is too degenerate to test.
// *screen_ptr receives the value used for BEST_CHISQ and for --epi2 counting:
// the screening statistic when the pair does not clear --epi1, and otherwise
// the fitted one, floored at the screening threshold.  (PLINK 1.x mixes the
// two this way, and the summary is only comparable if this does too.)
// *stat_ptr receives the fitted statistic, and *do_report_ptr whether the pair
// cleared the screening threshold and so is reported at all.  (The fitted
// statistic of a perfect fit comes back as rounding noise which can be
// negative, so its sign cannot double as that flag.)
// *df_adj_ptr receives the degrees-of-freedom reduction: df is 4 >> df_adj.
static uint32_t FepiBoost(const uint32_t* counts, const double* recip_cache, const double* alpha1sq, double* screen_ptr, double* stat_ptr, uint32_t* df_adj_ptr, uint32_t* do_report_ptr) {
  double p_ca[6];
  uint32_t df_adj = FepiBoostPCa(counts, recip_cache, p_ca);
  if (df_adj > 1) {
    return 1;
  }
  double p_bc[6];
  FepiBoostPBc(counts, recip_cache, p_bc);

  // p_ab[geno2 * 3 + geno1] = P(geno1 | geno2).
  double p_ab[9];
  uint32_t empty_col_ct = 0;
  uint32_t obs_ct = 0;
  for (uint32_t geno2 = 0; geno2 != 3; ++geno2) {
    uint32_t row_cts[3];
    for (uint32_t geno1 = 0; geno1 != 3; ++geno1) {
      row_cts[geno1] = counts[geno1 * 3 + geno2] + counts[9 + geno1 * 3 + geno2];
    }
    const uint32_t col_ct = row_cts[0] + row_cts[1] + row_cts[2];
    if (!col_ct) {
      if (empty_col_ct) {
        return 1;
      }
      empty_col_ct = 1;
      ++df_adj;
    }
    obs_ct += col_ct;
    const double tot_recip = recip_cache[col_ct];
    for (uint32_t geno1 = 0; geno1 != 3; ++geno1) {
      p_ab[geno2 * 3 + geno1] = u31tod(row_cts[geno1]) * tot_recip;
    }
  }
  *df_adj_ptr = df_adj;

  const double obs_ct_recip = recip_cache[obs_ct];
  double tau = 0.0;
  double screen_stat = 0.0;
  {
    const uint32_t* counts_iter = counts;
    for (uint32_t group_idx = 0; group_idx != 2; ++group_idx) {
      const double* cur_p_bc = &(p_bc[group_idx * 3]);
      for (uint32_t geno1 = 0; geno1 != 3; ++geno1) {
        const double cur_p_ca = p_ca[geno1 * 2 + group_idx];
        for (uint32_t geno2 = 0; geno2 != 3; ++geno2) {
          const double mu = p_ab[geno2 * 3 + geno1] * cur_p_bc[geno2] * cur_p_ca;
          tau += mu;
          const uint32_t obs = *counts_iter++;
          if (obs) {
            const double obs_d = u31tod(obs);
            if (mu != 0.0) {
              screen_stat -= obs_d * log(mu * recip_cache[obs]);
            } else {
              // An observed cell the approximation calls impossible would send
              // the sum to infinity; PLINK 1.x scores it as if mu were 1.
              screen_stat += obs_d * log(obs_d);
            }
          }
        }
      }
    }
  }
  screen_stat = 2 * (screen_stat + u31tod(obs_ct) * log(tau * obs_ct_recip));
  if (screen_stat <= alpha1sq[df_adj]) {
    *screen_ptr = screen_stat;
    *do_report_ptr = 0;
    return 0;
  }

  // Iterative proportional fitting, cycling over the three two-way margins
  // until the fitted table stops moving.  mu is [geno1 * 6 + geno2 * 2 +
  // group].
  // (Only 'no-firth' and 'firth-fallback' modes are supported, so there's no
  // need to extend this to mirror Firth regression without covariates.)
  double mu[18];
  for (uint32_t cell_idx = 0; cell_idx != 18; ++cell_idx) {
    mu[cell_idx] = 1.0;
  }
  double mu_err;
  do {
    double prev_mu[18];
    memcpy(prev_mu, mu, 18 * sizeof(double));
    for (uint32_t cell_idx = 0; cell_idx != 9; ++cell_idx) {
      double* cur_mu = &(mu[cell_idx * 2]);
      const double margin = cur_mu[0] + cur_mu[1];
      double scale = 0.0;
      if (margin != 0.0) {
        scale = u31tod(counts[cell_idx] + counts[cell_idx + 9]) / margin;
      }
      cur_mu[0] *= scale;
      cur_mu[1] *= scale;
    }
    for (uint32_t geno1 = 0; geno1 != 3; ++geno1) {
      for (uint32_t group_idx = 0; group_idx != 2; ++group_idx) {
        double* cur_mu = &(mu[geno1 * 6 + group_idx]);
        const double margin = cur_mu[0] + cur_mu[2] + cur_mu[4];
        double scale = 0.0;
        if (margin != 0.0) {
          const uint32_t* cur_counts = &(counts[group_idx * 9 + geno1 * 3]);
          scale = u31tod(cur_counts[0] + cur_counts[1] + cur_counts[2]) / margin;
        }
        cur_mu[0] *= scale;
        cur_mu[2] *= scale;
        cur_mu[4] *= scale;
      }
    }
    for (uint32_t geno2 = 0; geno2 != 3; ++geno2) {
      for (uint32_t group_idx = 0; group_idx != 2; ++group_idx) {
        double* cur_mu = &(mu[geno2 * 2 + group_idx]);
        const double margin = cur_mu[0] + cur_mu[6] + cur_mu[12];
        double scale = 0.0;
        if (margin != 0.0) {
          const uint32_t* cur_counts = &(counts[group_idx * 9 + geno2]);
          scale = u31tod(cur_counts[0] + cur_counts[3] + cur_counts[6]) / margin;
        }
        cur_mu[0] *= scale;
        cur_mu[6] *= scale;
        cur_mu[12] *= scale;
      }
    }
    mu_err = 0.0;
    for (uint32_t cell_idx = 0; cell_idx != 18; ++cell_idx) {
      mu_err += fabs(mu[cell_idx] - prev_mu[cell_idx]);
    }
  } while (mu_err > 0.001);

  double fit_stat = 0.0;
  tau = 0.0;
  {
    const uint32_t* counts_iter = counts;
    for (uint32_t group_idx = 0; group_idx != 2; ++group_idx) {
      for (uint32_t geno1 = 0; geno1 != 3; ++geno1) {
        for (uint32_t geno2 = 0; geno2 != 3; ++geno2) {
          const double obs_frac = u31tod(*counts_iter++) * obs_ct_recip;
          const double fit_frac = mu[geno1 * 6 + geno2 * 2 + group_idx] * obs_ct_recip;
          if (obs_frac != 0.0) {
            fit_stat += obs_frac * ((fit_frac != 0.0)? log(obs_frac / fit_frac) : log(obs_frac));
          }
          tau += fit_frac;
        }
      }
    }
  }
  fit_stat = (fit_stat + log(tau)) * u31tod(2 * obs_ct);
  *stat_ptr = fit_stat;
  *screen_ptr = MAXV(fit_stat, alpha1sq[df_adj]);
  *do_report_ptr = 1;
  return 0;
}

// Per-variant genotype bitvectors, one triple per group: index 0 is hom-REF,
// 1 is het, 2 is hom-ALT, and a missing call is in none of them.  A 3x3 cell
// count is then one PopcountWordsIntersect.
static void GenovecToGenoBits(const uintptr_t* genovec, uint32_t sample_ct, uintptr_t* hom_buf, uintptr_t* ref2het_buf, uintptr_t* dst) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  SplitHomRef2het(genovec, sample_ct, hom_buf, ref2het_buf);
  uintptr_t* hom_ref = dst;
  uintptr_t* het = &(dst[sample_ctl]);
  uintptr_t* hom_alt = &(dst[2 * sample_ctl]);
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    const uintptr_t hom_word = hom_buf[widx];
    const uintptr_t ref2het_word = ref2het_buf[widx];
    hom_ref[widx] = hom_word & ref2het_word;
    het[widx] = (~hom_word) & ref2het_word;
    hom_alt[widx] = hom_word & (~ref2het_word);
  }
  ZeroTrailingBits(sample_ct, hom_ref);
  ZeroTrailingBits(sample_ct, het);
  ZeroTrailingBits(sample_ct, hom_alt);
}

typedef struct EpiSummaryEntryStruct {
  uint32_t n_sig;
  uint32_t n_tot;
  double best_chisq;
  uint32_t best_vidx;
} EpiSummaryEntry;

// --epistasis-boost with covariates: the postprocessing step of Wan X et al.
// (2010).  The screen and the log-linear fit are left alone; only the pairs
// that clear --epi1 reach this, and those are refit by logistic regression of
// the phenotype on genotype dummies for both variants, their products, and the
// covariates.  The statistic is the likelihood ratio between the model
// carrying the product terms and the one without them.
//
// With no covariate this would only reproduce the log-linear statistic: the
// paper's Methods note that the homogeneous log-linear model is the equivalent
// form of the main-effect logistic model, and the saturated one of the full
// model.  The degrees of freedom therefore follow the same rule as the
// log-linear fit, a genotype value with no samples costing its whole row or
// column and leaving (r - 1) * (c - 1).
typedef struct EpiCovarCtxStruct {
  const double* covar_vals;
  const uintptr_t* case_bitvec;
  uint32_t analysis_ct;
  uint32_t covar_ct;
  uintptr_t max_sample_ctav;
  uintptr_t predictor_ctav;
  double* yy;
  double* xx;
  double* coef;
  double* ll;
  double* pp;
  double* vv;
  double* hh;
  double* grad;
  double* dcoef;
  MatrixInvertBuf1* inv_1d_buf;
  double* dbl_2d_buf;
  uint32_t* nm_sample_idxs;
  unsigned char* nm_row_levels;
  unsigned char* nm_col_levels;
} EpiCovarCtx;

// Log-likelihood of the fitted model.  LogisticRegressionD() leaves p - y in
// its own buffer rather than p, so this recomputes the linear predictor from
// the coefficients.
static double EpiCovarLoglik(const double* xx, const double* yy, const double* coef, uintptr_t sample_ctav, uint32_t sample_ct, uint32_t predictor_ct) {
  double loglik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    double eta = 0.0;
    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      eta += coef[pred_idx] * xx[pred_idx * sample_ctav + sample_idx];
    }
    // log(1 + exp(eta)) overflows for a large positive eta, so the dominant
    // term comes out of the logarithm.
    if (eta > 0.0) {
      loglik += (yy[sample_idx] - 1.0) * eta - log1p(exp(-eta));
    } else {
      loglik += yy[sample_idx] * eta - log1p(exp(eta));
    }
  }
  return loglik;
}

// row_geno_bits and col_geno_bits are the two variants' genotype bitvector
// triples over the analysis samples; counts is the 2x3x3 table the screen was
// computed from, whose two groups partition those same samples.
//
// row_geno_bits and col_geno_bits are the two variants' genotype bitvector
// triples over the analysis samples; counts is the 2x3x3 table the screen was
// computed from, whose two groups partition those same samples.
//
// The two models are built separately rather than as a prefix and an
// extension of one design, because the natural extension is not always full
// rank: a genotype combination with no samples leaves its product term with
// nothing to estimate.  The full model is instead the saturated one over the
// occupied cells, which is what the log-linear fit's saturated model is too.
// The degrees of freedom reported are still the table's, (rows - 1) times
// (columns - 1) over the levels that are present, so a sampling zero leaves
// them alone and only an empty row or column reduces them.
//
// Returns 1 if the pair cannot be fit, in which case it is left out of the
// report.
static uint32_t EpiCovarRefit(const uintptr_t* row_geno_bits, const uintptr_t* col_geno_bits, const uint32_t* counts, uint32_t analysis_ctl, EpiCovarCtx* ctx, double* stat_ptr, uint32_t* df_ptr) {
  uint32_t row_levels[3];
  uint32_t col_levels[3];
  uint32_t row_level_ct = 0;
  uint32_t col_level_ct = 0;
  for (uint32_t geno = 0; geno != 3; ++geno) {
    uint32_t row_tot = 0;
    uint32_t col_tot = 0;
    for (uint32_t other = 0; other != 3; ++other) {
      row_tot += counts[geno * 3 + other] + counts[9 + geno * 3 + other];
      col_tot += counts[other * 3 + geno] + counts[9 + other * 3 + geno];
    }
    if (row_tot) {
      row_levels[row_level_ct++] = geno;
    }
    if (col_tot) {
      col_levels[col_level_ct++] = geno;
    }
  }
  if ((row_level_ct < 2) || (col_level_ct < 2)) {
    return 1;
  }
  // Cell 0 is the reference; the rest get a dummy in the full model.
  uint32_t cell_cols[3][3];
  uint32_t occupied_ct = 0;
  for (uint32_t row_level_idx = 0; row_level_idx != row_level_ct; ++row_level_idx) {
    for (uint32_t col_level_idx = 0; col_level_idx != col_level_ct; ++col_level_idx) {
      const uint32_t cell_idx = row_levels[row_level_idx] * 3 + col_levels[col_level_idx];
      const uint32_t cell_ct = counts[cell_idx] + counts[9 + cell_idx];
      cell_cols[row_level_idx][col_level_idx] = cell_ct? occupied_ct++ : UINT32_MAX;
    }
  }
  const uint32_t covar_ct = ctx->covar_ct;
  const uint32_t main_predictor_ct = row_level_ct + col_level_ct - 1 + covar_ct;
  const uint32_t full_predictor_ct = occupied_ct + covar_ct;
  if (full_predictor_ct <= main_predictor_ct) {
    // Every remaining cell is determined by the margins, so there is no
    // interaction left to test.
    return 1;
  }

  // The samples of each occupied cell, in cell order, so that the design can
  // be filled without decoding a genotype.
  uint32_t* nm_sample_idxs = ctx->nm_sample_idxs;
  unsigned char* nm_row_levels = ctx->nm_row_levels;
  unsigned char* nm_col_levels = ctx->nm_col_levels;
  uint32_t nm_ct = 0;
  for (uint32_t row_level_idx = 0; row_level_idx != row_level_ct; ++row_level_idx) {
    const uintptr_t* row_vec = &(row_geno_bits[row_levels[row_level_idx] * analysis_ctl]);
    for (uint32_t col_level_idx = 0; col_level_idx != col_level_ct; ++col_level_idx) {
      const uintptr_t* col_vec = &(col_geno_bits[col_levels[col_level_idx] * analysis_ctl]);
      for (uint32_t widx = 0; widx != analysis_ctl; ++widx) {
        uintptr_t cur_word = row_vec[widx] & col_vec[widx];
        const uint32_t sample_idx_base = widx * kBitsPerWord;
        while (cur_word) {
          nm_sample_idxs[nm_ct] = sample_idx_base + ctzw(cur_word);
          nm_row_levels[nm_ct] = row_level_idx;
          nm_col_levels[nm_ct] = col_level_idx;
          ++nm_ct;
          cur_word &= cur_word - 1;
        }
      }
    }
  }
  if (nm_ct <= full_predictor_ct) {
    return 1;
  }
  // LogisticRegressionD() derives its own column stride from the sample count,
  // so the design has to be laid out at that stride.
  const uintptr_t sample_ctav = RoundUpPow2(nm_ct, kDoublePerDVec);
  double* xx = ctx->xx;
  double* yy = ctx->yy;
  const double* covar_vals = ctx->covar_vals;
  const uintptr_t covar_stride = ctx->analysis_ct;
  ZeroDArr(sample_ctav, yy);
  for (uint32_t sample_idx = 0; sample_idx != nm_ct; ++sample_idx) {
    yy[sample_idx] = u31tod(IsSet(ctx->case_bitvec, nm_sample_idxs[sample_idx]));
  }
  const uintptr_t predictor_ctav = ctx->predictor_ctav;
  double loglik[2];
  const uint32_t predictor_cts[2] = {main_predictor_ct, full_predictor_ct};
  for (uint32_t model_idx = 0; model_idx != 2; ++model_idx) {
    const uint32_t cur_predictor_ct = predictor_cts[model_idx];
    for (uint32_t pred_idx = 0; pred_idx != cur_predictor_ct; ++pred_idx) {
      ZeroDArr(sample_ctav, &(xx[pred_idx * sample_ctav]));
    }
    // Model 0: intercept, one dummy per non-reference genotype of each
    // variant, covariates.  Model 1: intercept, one dummy per non-reference
    // occupied cell, covariates.
    const uint32_t covar_col_start = cur_predictor_ct - covar_ct;
    for (uint32_t sample_idx = 0; sample_idx != nm_ct; ++sample_idx) {
      const uint32_t analysis_idx = nm_sample_idxs[sample_idx];
      const uint32_t row_level_idx = nm_row_levels[sample_idx];
      const uint32_t col_level_idx = nm_col_levels[sample_idx];
      xx[sample_idx] = 1.0;
      if (!model_idx) {
        if (row_level_idx) {
          xx[row_level_idx * sample_ctav + sample_idx] = 1.0;
        }
        if (col_level_idx) {
          xx[(row_level_ct - 1 + col_level_idx) * sample_ctav + sample_idx] = 1.0;
        }
      } else {
        const uint32_t cell_col = cell_cols[row_level_idx][col_level_idx];
        if (cell_col) {
          xx[cell_col * sample_ctav + sample_idx] = 1.0;
        }
      }
      for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
        xx[(covar_col_start + covar_idx) * sample_ctav + sample_idx] = covar_vals[covar_idx * covar_stride + analysis_idx];
      }
    }
    uint32_t is_unfinished = 0;
    ZeroDArr(predictor_ctav, ctx->coef);
    if (LogisticRegressionD(yy, xx, nullptr, nm_ct, cur_predictor_ct, &is_unfinished, ctx->coef, ctx->ll, ctx->pp, ctx->vv, ctx->hh, ctx->grad, ctx->dcoef, ctx->inv_1d_buf, ctx->dbl_2d_buf) || is_unfinished) {
      return 1;
    }
    loglik[model_idx] = EpiCovarLoglik(xx, yy, ctx->coef, sample_ctav, nm_ct, cur_predictor_ct);
  }
  const double stat = 2 * (loglik[1] - loglik[0]);
  if (!isfinite(stat)) {
    return 1;
  }
  // A perfect fit lands on zero from either side.
  *stat_ptr = MAXV(stat, 0.0);
  // A sampling zero costs the full model a parameter, but the degrees of
  // freedom stay those of the table it is a zero in, as on the unadjusted
  // side: only an empty row or column, which drops a level outright, reduces
  // them.
  *df_ptr = (row_level_ct - 1) * (col_level_ct - 1);
  return 0;
}

PglErr CalcEpiBoost(const uintptr_t* orig_sample_include, const PhenoCol* pheno_cols, const PhenoCol* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const EpiInfo* epi_ip, uint32_t raw_sample_ct, uint32_t pheno_ct, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_allele_slen, double output_min_ln, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  char* cswritetp = nullptr;
  CompressStreamState css;
  CompressStreamState csst;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  PreinitCstream(&csst);
  {
    const EpiFlags flags = epi_ip->flags;

    // PLINK 1.x had a single phenotype.  Rather than guess which of several
    // loaded case/control phenotypes an O(variant_ct^2) scan was meant for,
    // ask.
    uint32_t pheno_idx = UINT32_MAX;
    uint32_t cc_pheno_ct = 0;
    for (uint32_t uii = 0; uii != pheno_ct; ++uii) {
      if (pheno_cols[uii].type_code == kPhenoDtypeCc) {
        ++cc_pheno_ct;
        if (pheno_idx == UINT32_MAX) {
          pheno_idx = uii;
        }
      }
    }
    if (unlikely(!cc_pheno_ct)) {
      logerrputs("Error: --epistasis-boost requires a case/control phenotype.\n");
      goto CalcEpiBoost_ret_INCONSISTENT_INPUT;
    }
    if (unlikely(cc_pheno_ct > 1)) {
      logerrputs("Error: --epistasis-boost needs exactly one case/control phenotype; select one\nwith --pheno-name.\n");
      goto CalcEpiBoost_ret_INCONSISTENT_INPUT;
    }
    const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);

    // A sample missing a covariate is left out of the whole scan, not just of
    // the covariate-adjusted refit, so that the screen and the refit are over
    // the same samples; this is what --glm does with a missing covariate.
    uintptr_t* analysis_include;
    uintptr_t* case_include;
    uintptr_t* ctrl_include;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &analysis_include) ||
                 bigstack_alloc_w(raw_sample_ctl, &case_include) ||
                 bigstack_alloc_w(raw_sample_ctl, &ctrl_include))) {
      goto CalcEpiBoost_ret_NOMEM;
    }
    BitvecAndCopy(orig_sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, analysis_include);
    for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
      const PhenoCol* cur_covar_col = &(covar_cols[covar_idx]);
      if (unlikely(cur_covar_col->type_code == kPhenoDtypeCat)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --epistasis-boost does not support categorical covariates yet ('%s').  Split it into binary covariates with --split-cat-pheno first.\n", &(covar_names[covar_idx * max_covar_name_blen]));
        goto CalcEpiBoost_ret_INCONSISTENT_INPUT_WW;
      }
      BitvecAnd(cur_covar_col->nonmiss, raw_sample_ctl, analysis_include);
    }
    memcpy(ctrl_include, analysis_include, raw_sample_ctl * sizeof(intptr_t));
    BitvecAndCopy(ctrl_include, cur_pheno_col->data.cc, raw_sample_ctl, case_include);
    BitvecInvmask(case_include, raw_sample_ctl, ctrl_include);
    const uint32_t analysis_ct = PopcountWords(analysis_include, raw_sample_ctl);
    const uint32_t case_ct = PopcountWords(case_include, raw_sample_ctl);
    const uint32_t ctrl_ct = PopcountWords(ctrl_include, raw_sample_ctl);
    if (unlikely(case_ct < 2)) {
      logerrputs("Error: --epistasis-boost requires at least 2 cases.\n");
      goto CalcEpiBoost_ret_DEGENERATE_DATA;
    }
    if (unlikely(ctrl_ct < 2)) {
      logerrputs("Error: --epistasis-boost requires at least 2 controls.\n");
      goto CalcEpiBoost_ret_DEGENERATE_DATA;
    }
    // The covariates are stored one column at a time, since the refit walks a
    // single covariate across the samples of one pair.
    double* covar_vals = nullptr;
    if (covar_ct) {
      if (unlikely(bigstack_alloc_d(S_CAST(uintptr_t, analysis_ct) * covar_ct, &covar_vals))) {
        goto CalcEpiBoost_ret_NOMEM;
      }
      uintptr_t sample_uidx_base = 0;
      uintptr_t analysis_include_bits = analysis_include[0];
      for (uint32_t sample_idx = 0; sample_idx != analysis_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(analysis_include, &sample_uidx_base, &analysis_include_bits);
        for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
          const PhenoCol* cur_covar_col = &(covar_cols[covar_idx]);
          const double cur_val = (cur_covar_col->type_code == kPhenoDtypeQt)? cur_covar_col->data.qt[sample_uidx] : u31tod(IsSet(cur_covar_col->data.cc, sample_uidx));
          covar_vals[covar_idx * S_CAST(uintptr_t, analysis_ct) + sample_idx] = cur_val;
        }
      }
    }
    // A covariate that is constant across the analysis samples is collinear
    // with the intercept, so it is dropped rather than left to make every fit
    // fail, as in --glm.
    uint32_t cur_covar_ct = 0;
    for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
      const double* cur_col = &(covar_vals[covar_idx * S_CAST(uintptr_t, analysis_ct)]);
      const double first_val = cur_col[0];
      uint32_t sample_idx = 1;
      for (; sample_idx != analysis_ct; ++sample_idx) {
        if (cur_col[sample_idx] != first_val) {
          break;
        }
      }
      if (sample_idx == analysis_ct) {
        logerrprintf("Warning: Excluding constant covariate '%s' from --epistasis-boost.\n", &(covar_names[covar_idx * max_covar_name_blen]));
        continue;
      }
      if (cur_covar_ct != covar_idx) {
        memcpy(&(covar_vals[cur_covar_ct * S_CAST(uintptr_t, analysis_ct)]), cur_col, analysis_ct * sizeof(double));
      }
      ++cur_covar_ct;
    }
    // With no covariate left there is nothing for the refit to adjust for, and
    // the logistic likelihood ratio would only reproduce the log-linear fit,
    // so the screen-and-fit path is used unchanged.
    if (unlikely(cur_covar_ct && (!(flags & kfEpiNoFirth)))) {
      // temporary
      logerrputs("Error: --epistasis-boost: Firth regression is not implemented yet.  Specify\n'no-firth' for now.\n");
      reterr = kPglRetNotYetSupported;
      goto CalcEpiBoost_ret_1;
    }
    const uint32_t covar_refit = (cur_covar_ct != 0);

    const uint32_t group_ct = 2;
    const uint32_t group_cts[3] = {case_ct, ctrl_ct, analysis_ct};
    const uintptr_t* group_includes[3] = {case_include, ctrl_include, analysis_include};
    const uint32_t case_ctl = BitCtToWordCt(case_ct);
    const uint32_t ctrl_ctl = BitCtToWordCt(ctrl_ct);
    const uint32_t analysis_ctl = BitCtToWordCt(analysis_ct);
    // The refit needs each variant's genotypes over the analysis samples, not
    // split by phenotype, so a third bitvector triple is loaded alongside the
    // two the 2x3x3 table is counted from.  It is only allocated when there is
    // a covariate to adjust for.
    const uint32_t load_group_ct = group_ct + covar_refit;
    const uint32_t group_ctls[3] = {case_ctl, ctrl_ctl, analysis_ctl};
    const uintptr_t words_per_variant = 3 * (S_CAST(uintptr_t, case_ctl + ctrl_ctl) + (covar_refit? analysis_ctl : 0));

    // Non-autosomal variants are left out, as in PLINK 1.x.
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* variant_include;
    if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include))) {
      goto CalcEpiBoost_ret_NOMEM;
    }
    memcpy(variant_include, orig_variant_include, raw_variant_ctl * sizeof(intptr_t));
    const uint32_t chr_ct = cip->chr_ct;
    uint32_t nonautosomal_ct = 0;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if (!IsSet(cip->haploid_mask, chr_idx) && (chr_idx != S_CAST(uint32_t, cip->xymt_codes[kChrOffsetMT]))) {
        continue;
      }
      const uint32_t chr_vidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
      const uint32_t chr_vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      nonautosomal_ct += PopcountBitRange(variant_include, chr_vidx_start, chr_vidx_end);
      ClearBitsNz(chr_vidx_start, chr_vidx_end, variant_include);
    }
    uint32_t variant_ct = orig_variant_ct - nonautosomal_ct;
    if (unlikely(variant_ct < 2)) {
      logerrputs("Error: --epistasis-boost requires at least 2 autosomal variants.\n");
      goto CalcEpiBoost_ret_DEGENERATE_DATA;
    }

    // Monomorphic variants are dropped up front, over the whole phenotyped
    // set rather than group by group, as in PLINK 1.x.  An empty genotype row
    // or column within a group is not degenerate here: it costs two degrees
    // of freedom instead.
    uint32_t* sample_include_cumulative_popcounts[3];
    PgrSampleSubsetIndex pssis[3];
    for (uint32_t group_idx = 0; group_idx != group_ct; ++group_idx) {
      if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &(sample_include_cumulative_popcounts[group_idx])))) {
        goto CalcEpiBoost_ret_NOMEM;
      }
      FillCumulativePopcounts(group_includes[group_idx], raw_sample_ctl, sample_include_cumulative_popcounts[group_idx]);
    }
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &(sample_include_cumulative_popcounts[2])))) {
      goto CalcEpiBoost_ret_NOMEM;
    }
    FillCumulativePopcounts(analysis_include, raw_sample_ctl, sample_include_cumulative_popcounts[2]);
    const uint32_t max_group_ct = MAXV(analysis_ct, MAXV(case_ct, ctrl_ct));
    uintptr_t* genovec;
    uintptr_t* hom_buf;
    uintptr_t* ref2het_buf;
    if (unlikely(bigstack_alloc_w(NypCtToWordCt(max_group_ct), &genovec) ||
                 bigstack_alloc_w(BitCtToWordCt(max_group_ct), &hom_buf) ||
                 bigstack_alloc_w(BitCtToWordCt(max_group_ct), &ref2het_buf))) {
      goto CalcEpiBoost_ret_NOMEM;
    }
    {
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      uint32_t skipped_ct = 0;
      PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts[2], simple_pgrp, &(pssis[2]));
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        reterr = PgrGet(analysis_include, pssis[2], analysis_ct, variant_uidx, simple_pgrp, genovec);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, variant_uidx);
          goto CalcEpiBoost_ret_1;
        }
        // GenoarrCountFreqsUnsafe() counts whole words and derives the hom-REF
        // count by subtraction, so the nyps past sample_ct have to be cleared
        // first; PgrGet() leaves them alone.
        ZeroTrailingNyps(analysis_ct, genovec);
        STD_ARRAY_DECL(uint32_t, 4, genocounts);
        GenoarrCountFreqsUnsafe(genovec, analysis_ct, genocounts);
        const uint32_t nonmiss_ct = genocounts[0] + genocounts[1] + genocounts[2];
        if ((genocounts[0] == nonmiss_ct) || (genocounts[1] == nonmiss_ct) || (genocounts[2] == nonmiss_ct)) {
          ClearBit(variant_uidx, variant_include);
          ++skipped_ct;
        }
      }
      if (skipped_ct) {
        variant_ct -= skipped_ct;
        logprintf("--epistasis-boost: Skipping %u monomorphic variant%s.\n", skipped_ct, (skipped_ct == 1)? "" : "s");
      }
      if (unlikely(variant_ct < 2)) {
        logerrputs("Error: --epistasis-boost has fewer than 2 usable variants left.\n");
        goto CalcEpiBoost_ret_DEGENERATE_DATA;
      }
    }
    if (nonautosomal_ct) {
      logprintf("--epistasis-boost: Skipping %u non-autosomal variant%s.\n", nonautosomal_ct, (nonautosomal_ct == 1)? "" : "s");
    }

    uint32_t* variant_uidxs;
    uint32_t* variant_chr_fo_idxs;
    if (unlikely(bigstack_alloc_u32(variant_ct, &variant_uidxs) ||
                 bigstack_alloc_u32(variant_ct, &variant_chr_fo_idxs))) {
      goto CalcEpiBoost_ret_NOMEM;
    }
    {
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      uint32_t chr_fo_idx = 0;
      uint32_t chr_vidx_end = cip->chr_fo_vidx_start[1];
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        while (variant_uidx >= chr_vidx_end) {
          ++chr_fo_idx;
          chr_vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        }
        variant_uidxs[variant_idx] = variant_uidx;
        variant_chr_fo_idxs[variant_idx] = chr_fo_idx;
      }
    }

    // Rows are split across --parallel jobs; every job still needs to reach
    // every column to its right.
    // ParallelBounds() splits a triangle whose row r carries r entries below
    // it; this scan's row r carries variant_ct - 1 - r entries above it.  The
    // two are mirror images, so the bounds come back in the mirrored index and
    // get reflected here.
    uint32_t mirror_start;
    uint32_t mirror_end;
    // The reflection also reverses the job order, so ask for the mirrored job
    // index as well; job 1 then produces the first rows and concatenating the
    // jobs in order reproduces a single run.
    ParallelBounds(variant_ct, 1, parallel_tot - 1 - parallel_idx, parallel_tot, R_CAST(int32_t*, &mirror_start), R_CAST(int32_t*, &mirror_end));
    const uint32_t row_start_idx = variant_ct - mirror_end;
    const uint32_t row_end_idx = variant_ct - mirror_start;
    if (row_start_idx == row_end_idx) {
      logerrputs("Warning: This --parallel job has no rows to scan.\n");
    }

    // boost divides by cell and margin counts constantly, and there are only
    // analysis_ct + 1 possible denominators.
    double* recip_cache;
    if (unlikely(bigstack_alloc_d(analysis_ct + 1, &recip_cache))) {
      goto CalcEpiBoost_ret_NOMEM;
    }
    recip_cache[0] = 0.0;
    for (uint32_t uii = 1; uii <= analysis_ct; ++uii) {
      recip_cache[uii] = 1.0 / u31tod(uii);
    }

    // The refit's workspace, claimed before the variant blocks take what is
    // left of bigstack.
    EpiCovarCtx covar_ctx;
    {
      // Intercept, up to two dummies per variant, up to four products, and the
      // covariates.
      const uint32_t max_predictor_ct = 9 + cur_covar_ct;
      const uintptr_t max_sample_ctav = RoundUpPow2(analysis_ct, kDoublePerDVec);
      const uintptr_t predictor_ctav = RoundUpPow2(max_predictor_ct, kDoublePerDVec);
      // initialize even in !covar_refit case to address spurious compiler
      // warnings for now
      covar_ctx.covar_vals = covar_vals;
      covar_ctx.case_bitvec = nullptr;
      covar_ctx.analysis_ct = analysis_ct;
      covar_ctx.covar_ct = cur_covar_ct;
      covar_ctx.max_sample_ctav = max_sample_ctav;
      covar_ctx.predictor_ctav = predictor_ctav;
      covar_ctx.yy = nullptr;
      covar_ctx.xx = nullptr;
      covar_ctx.coef = nullptr;
      covar_ctx.ll = nullptr;
      covar_ctx.pp = nullptr;
      covar_ctx.vv = nullptr;
      covar_ctx.hh = nullptr;
      covar_ctx.grad = nullptr;
      covar_ctx.dcoef = nullptr;
      covar_ctx.inv_1d_buf = nullptr;
      covar_ctx.dbl_2d_buf = nullptr;
      covar_ctx.nm_sample_idxs = nullptr;
      covar_ctx.nm_row_levels = nullptr;
      covar_ctx.nm_col_levels = nullptr;
      if (covar_refit) {
        covar_ctx.inv_1d_buf = S_CAST(MatrixInvertBuf1*, bigstack_alloc(max_predictor_ct * kMatrixInvertBuf1CheckedAlloc));
        if (unlikely((!covar_ctx.inv_1d_buf) ||
                     bigstack_alloc_d(max_sample_ctav, &covar_ctx.yy) ||
                     bigstack_alloc_d(max_sample_ctav * max_predictor_ct, &covar_ctx.xx) ||
                     bigstack_alloc_d(predictor_ctav, &covar_ctx.coef) ||
                     bigstack_alloc_d(predictor_ctav * max_predictor_ct, &covar_ctx.ll) ||
                     bigstack_alloc_d(max_sample_ctav, &covar_ctx.pp) ||
                     bigstack_alloc_d(max_sample_ctav, &covar_ctx.vv) ||
                     bigstack_alloc_d(predictor_ctav * max_predictor_ct, &covar_ctx.hh) ||
                     bigstack_alloc_d(predictor_ctav, &covar_ctx.grad) ||
                     bigstack_alloc_d(predictor_ctav, &covar_ctx.dcoef) ||
                     bigstack_alloc_d(max_predictor_ct * MAXV(max_predictor_ct, 7), &covar_ctx.dbl_2d_buf) ||
                     bigstack_alloc_u32(analysis_ct, &covar_ctx.nm_sample_idxs) ||
                     bigstack_alloc_uc(analysis_ct, &covar_ctx.nm_row_levels) ||
                     bigstack_alloc_uc(analysis_ct, &covar_ctx.nm_col_levels))) {
          goto CalcEpiBoost_ret_NOMEM;
        }
        // The refit indexes its phenotype by analysis-sample index, so the
        // case set is needed in those coordinates rather than the raw ones.
        uintptr_t* case_collapsed;
        if (unlikely(bigstack_alloc_w(analysis_ctl, &case_collapsed))) {
          goto CalcEpiBoost_ret_NOMEM;
        }
        CopyBitarrSubset(case_include, analysis_include, analysis_ct, case_collapsed);
        covar_ctx.case_bitvec = case_collapsed;
      }
    }

    // Two variant blocks are held at once: the rows, and the columns they are
    // being tested against.  Everything to the right of a row has to be
    // reachable, so the column block sweeps the whole range each time the row
    // block advances.
    uintptr_t bytes_per_variant = words_per_variant * sizeof(intptr_t);
    const uintptr_t max_slot_ct = (bigstack_left() / 2) / bytes_per_variant;
    if (unlikely(max_slot_ct < 4)) {
      goto CalcEpiBoost_ret_NOMEM;
    }
    // The report has to come out in row-major order, so that concatenating
    // --parallel jobs reproduces a single run, as PLINK 1.x promises.  A row
    // block only preserves that order while the whole column range fits in one
    // block; when it does not, the row block drops to a single row, so its
    // columns are still swept in order.  That costs a reread of the column
    // range per row, but only in the case that was already going to be
    // dominated by rereads.
    uint32_t col_block_size;
    uint32_t row_block_size;
    if (max_slot_ct > variant_ct) {
      col_block_size = variant_ct;
      row_block_size = MINV(max_slot_ct - variant_ct, variant_ct);
    } else {
      col_block_size = max_slot_ct - 1;
      row_block_size = 1;
    }
    uintptr_t* row_bits;
    uintptr_t* col_bits;
    if (unlikely(bigstack_alloc_w(row_block_size * words_per_variant, &row_bits) ||
                 bigstack_alloc_w(col_block_size * words_per_variant, &col_bits))) {
      goto CalcEpiBoost_ret_NOMEM;
    }

    EpiSummaryEntry* summary;
    if (unlikely(BIGSTACK_ALLOC_X(EpiSummaryEntry, variant_ct, &summary))) {
      goto CalcEpiBoost_ret_NOMEM;
    }
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      summary[variant_idx].n_sig = 0;
      summary[variant_idx].n_tot = 0;
      summary[variant_idx].best_chisq = -1.0;
      summary[variant_idx].best_vidx = UINT32_MAX;
    }

    // --epi1 is a screening threshold here, deciding which pairs are fit at
    // all, so PLINK 1.x's default for it is much stricter than a report
    // filter's.
    double ln_alpha1 = epi_ip->ln_epi1;
    if (ln_alpha1 > 0.0) {
      ln_alpha1 = -5 * kLn10 - kLn2;  // log(5e-6)
    }
    // df varies from pair to pair, so the thresholds are chi-square quantiles
    // at df 4, 2 and 1 rather than one p-value cutoff.
    double alpha1sq[3];
    double alpha2sq[3];
    alpha1sq[0] = LnPToChisq(ln_alpha1, 4);
    alpha1sq[1] = LnPToChisq(ln_alpha1, 2);
    alpha1sq[2] = LnPToChisq(ln_alpha1, 1);
    alpha2sq[0] = LnPToChisq(epi_ip->ln_epi2, 4);
    if (alpha1sq[0] == alpha2sq[0]) {
      // --epi1 and --epi2 agree: count the pairs that clear the fit rather
      // than the ones that cleared the screen.
      alpha2sq[0] *= 1 + kSmallEpsilon;
      alpha2sq[1] = alpha1sq[1] * (1 + kSmallEpsilon);
      alpha2sq[2] = alpha1sq[2] * (1 + kSmallEpsilon);
    } else {
      alpha2sq[1] = LnPToChisq(epi_ip->ln_epi2, 2);
      alpha2sq[2] = LnPToChisq(epi_ip->ln_epi2, 1);
    }
    const uint32_t output_zst = (flags / kfEpiZs) & 1;
    char* outname_end2 = strcpya_k(outname_end, ".epi.cc");
    // Main report is <prefix>.epi.cc[.<job>], summary is
    // <prefix>.epi.cc.summary[.<job>], as in PLINK 1.x.
    char* main_end = outname_end2;
    if (parallel_tot > 1) {
      *main_end++ = '.';
      main_end = u32toa(parallel_idx + 1, main_end);
    }
    if (output_zst) {
      snprintf(main_end, kMaxOutfnameExtBlen - S_CAST(uintptr_t, main_end - outname_end), ".zst");
    } else {
      *main_end = '\0';
    }
    const uint32_t chrom_col = flags & kfEpiColChrom;
    const uint32_t pos_col = flags & kfEpiColPos;
    const uint32_t a1_col = (flags & kfEpiColA1) || ((flags & kfEpiColMaybeA1) && MultiallelicVariantPresent(variant_include, allele_idx_offsets, variant_ct));
    const uint32_t stat_col = flags & kfEpiColStat;
    const uint32_t df_col = flags & kfEpiColDf;
    const uint32_t p_col = flags & kfEpiColP;
    uint32_t max_chr_slen = 0;
    if (chrom_col) {
      max_chr_slen = GetMaxChrSlen(cip);
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + 2 * max_chr_slen + 2 * a1_col * max_allele_slen + 2 * kMaxIdSlen + 256;
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto CalcEpiBoost_ret_1;
    }
    if (!parallel_idx) {
      *cswritep++ = '#';
      if (chrom_col) {
        cswritep = strcpya_k(cswritep, "CHROM1\t");
      }
      if (pos_col) {
        cswritep = strcpya_k(cswritep, "POS1\t");
      }
      cswritep = strcpya_k(cswritep, "ID1");
      if (a1_col) {
        cswritep = strcpya_k(cswritep, "\tALLELE1");
      }
      if (chrom_col) {
        cswritep = strcpya_k(cswritep, "\tCHROM2");
      }
      if (pos_col) {
        cswritep = strcpya_k(cswritep, "\tPOS2");
      }
      cswritep = strcpya_k(cswritep, "\tID2");
      if (a1_col) {
        cswritep = strcpya_k(cswritep, "\tALLELE2");
      }
      if (stat_col) {
        cswritep = strcpya_k(cswritep, "\tSTAT");
      }
      if (df_col) {
        cswritep = strcpya_k(cswritep, "\tDF");
      }
      if (p_col) {
        cswritep = strcpya_k(cswritep, "\tP");
      }
      AppendBinaryEoln(&cswritep);
    }

    // Loads variant_idxs [block_start, block_end) into dst.
    // (Declared as a lambda-free helper loop to keep the reader in one place.)
    uint64_t pair_ct_total = 0;
    for (uint32_t row_idx = row_start_idx; row_idx != row_end_idx; ++row_idx) {
      pair_ct_total += variant_ct - row_idx - 1;
    }
    uint64_t pairs_reported = 0;
    uint64_t refit_fail_ct = 0;
    fputs("--epistasis-boost: 0%", stdout);
    fflush(stdout);
    uint64_t next_print_pair = pair_ct_total / 100;
    uint32_t pct = 0;
    uint64_t pairs_seen = 0;

    if (flags & kfEpiRefBased) {
      maj_alleles = nullptr;
    }
    AlleleCode aidx = 0;
    for (uint32_t row_block_start = row_start_idx; row_block_start < row_end_idx; row_block_start += row_block_size) {
      const uint32_t row_block_end = MINV(row_block_start + row_block_size, row_end_idx);
      const uint32_t cur_row_ct = row_block_end - row_block_start;
      for (uint32_t slot_idx = 0; slot_idx != cur_row_ct; ++slot_idx) {
        uintptr_t* dst = &(row_bits[slot_idx * words_per_variant]);
        for (uint32_t group_idx = 0; group_idx != load_group_ct; ++group_idx) {
          PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts[group_idx], simple_pgrp, &(pssis[group_idx]));
          const uint32_t variant_uidx = variant_uidxs[row_block_start + slot_idx];
          if (maj_alleles) {
            aidx = maj_alleles[variant_uidx];
          }
          reterr = PgrGetInv1(group_includes[group_idx], pssis[group_idx], group_cts[group_idx], variant_uidx, aidx, simple_pgrp, genovec);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidxs[row_block_start + slot_idx]);
            goto CalcEpiBoost_ret_1;
          }
          GenovecToGenoBits(genovec, group_cts[group_idx], hom_buf, ref2het_buf, dst);
          dst = &(dst[3 * group_ctls[group_idx]]);
        }
      }
      for (uint32_t col_block_start = row_block_start; col_block_start < variant_ct; col_block_start += col_block_size) {
        const uint32_t col_block_end = MINV(col_block_start + col_block_size, variant_ct);
        const uint32_t cur_col_ct = col_block_end - col_block_start;
        for (uint32_t slot_idx = 0; slot_idx != cur_col_ct; ++slot_idx) {
          uintptr_t* dst = &(col_bits[slot_idx * words_per_variant]);
          for (uint32_t group_idx = 0; group_idx != load_group_ct; ++group_idx) {
            PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts[group_idx], simple_pgrp, &(pssis[group_idx]));
            const uint32_t variant_uidx = variant_uidxs[col_block_start + slot_idx];
            if (maj_alleles) {
              aidx = maj_alleles[variant_uidx];
            }
            reterr = PgrGetInv1(group_includes[group_idx], pssis[group_idx], group_cts[group_idx], variant_uidx, aidx, simple_pgrp, genovec);
            if (unlikely(reterr)) {
              PgenErrPrintNV(reterr, variant_uidxs[col_block_start + slot_idx]);
              goto CalcEpiBoost_ret_1;
            }
            GenovecToGenoBits(genovec, group_cts[group_idx], hom_buf, ref2het_buf, dst);
            dst = &(dst[3 * group_ctls[group_idx]]);
          }
        }
        for (uint32_t row_idx = row_block_start; row_idx != row_block_end; ++row_idx) {
          const uintptr_t* row_slot = &(row_bits[(row_idx - row_block_start) * words_per_variant]);
          const uint32_t col_first = MAXV(col_block_start, row_idx + 1);
          for (uint32_t col_idx = col_first; col_idx != col_block_end; ++col_idx) {
            ++pairs_seen;
            const uintptr_t* col_slot = &(col_bits[(col_idx - col_block_start) * words_per_variant]);
            uint32_t counts[18];
            {
              const uintptr_t* row_iter = row_slot;
              const uintptr_t* col_iter = col_slot;
              for (uint32_t group_idx = 0; group_idx != group_ct; ++group_idx) {
                const uint32_t cur_ctl = group_ctls[group_idx];
                uint32_t* cur_counts = &(counts[group_idx * 9]);
                for (uint32_t geno1 = 0; geno1 != 3; ++geno1) {
                  const uintptr_t* row_vec = &(row_iter[geno1 * cur_ctl]);
                  for (uint32_t geno2 = 0; geno2 != 3; ++geno2) {
                    cur_counts[geno1 * 3 + geno2] = PopcountWordsIntersect(row_vec, &(col_iter[geno2 * cur_ctl]), cur_ctl);
                  }
                }
                row_iter = &(row_iter[3 * cur_ctl]);
                col_iter = &(col_iter[3 * cur_ctl]);
              }
            }
            // chisq is the screening statistic, which drives BEST_CHISQ and
            // the --epi2 count; report_stat is the fitted one the report
            // prints.
            double chisq;
            double report_stat = 0.0;
            double ln_pval = 0.0;
            uint32_t df_adj;
            uint32_t do_report;
            if (FepiBoost(counts, recip_cache, alpha1sq, &chisq, &report_stat, &df_adj, &do_report)) {
              continue;
            }
            const uint32_t is_sig = (chisq >= alpha2sq[df_adj]);
            uint32_t report_df = 4 >> df_adj;
            if (do_report && covar_refit) {
              // The covariate-adjusted statistic replaces the log-linear one.
              // The screen is what the summary is built from, so it is left
              // alone, and a pair the refit cannot fit still counts towards
              // N_TOT: it was tested, it just has no adjusted statistic to
              // report.
              const uintptr_t analysis_offset = 3 * S_CAST(uintptr_t, case_ctl + ctrl_ctl);
              double adj_stat;
              uint32_t adj_df;
              if (EpiCovarRefit(&(row_slot[analysis_offset]), &(col_slot[analysis_offset]), counts, analysis_ctl, &covar_ctx, &adj_stat, &adj_df)) {
                ++refit_fail_ct;
                do_report = 0;
              } else {
                report_stat = adj_stat;
                report_df = adj_df;
              }
            }
            if (do_report) {
              // A perfect fit lands on zero from either side; the p-value of a
              // negative statistic is the p-value of zero.
              ln_pval = ChisqToLnP(MAXV(report_stat, 0.0), report_df);
            }
            summary[row_idx].n_tot += 1;
            summary[col_idx].n_tot += 1;
            if (is_sig) {
              summary[row_idx].n_sig += 1;
              summary[col_idx].n_sig += 1;
            }
            if (chisq > summary[row_idx].best_chisq) {
              summary[row_idx].best_chisq = chisq;
              summary[row_idx].best_vidx = col_idx;
            }
            if (chisq > summary[col_idx].best_chisq) {
              summary[col_idx].best_chisq = chisq;
              summary[col_idx].best_vidx = row_idx;
            }
            if (do_report) {
              ++pairs_reported;
              if (chrom_col) {
                cswritep = chrtoa(cip, cip->chr_file_order[variant_chr_fo_idxs[row_idx]], cswritep);
                *cswritep++ = '\t';
              }
              const uint32_t row_variant_uidx = variant_uidxs[row_idx];
              if (pos_col) {
                cswritep = u32toa_x(variant_bps[row_variant_uidx], '\t', cswritep);
              }
              cswritep = strcpyax(cswritep, variant_ids[row_variant_uidx], '\t');
              if (a1_col) {
                uintptr_t allele_idx_offset_base = row_variant_uidx * 2;
                if (allele_idx_offsets) {
                  allele_idx_offset_base = allele_idx_offsets[row_variant_uidx];
                }
                const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
                if (maj_alleles) {
                  aidx = maj_alleles[row_variant_uidx];
                }
                cswritep = strcpyax(cswritep, cur_alleles[aidx], '\t');
              }
              if (chrom_col) {
                cswritep = chrtoa(cip, cip->chr_file_order[variant_chr_fo_idxs[col_idx]], cswritep);
                *cswritep++ = '\t';
              }
              const uint32_t col_variant_uidx = variant_uidxs[col_idx];
              if (pos_col) {
                cswritep = u32toa_x(variant_bps[col_variant_uidx], '\t', cswritep);
              }
              cswritep = strcpya(cswritep, variant_ids[col_variant_uidx]);
              if (a1_col) {
                *cswritep++ = '\t';
                uintptr_t allele_idx_offset_base = col_variant_uidx * 2;
                if (allele_idx_offsets) {
                  allele_idx_offset_base = allele_idx_offsets[col_variant_uidx];
                }
                const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
                if (maj_alleles) {
                  aidx = maj_alleles[col_variant_uidx];
                }
                cswritep = strcpya(cswritep, cur_alleles[aidx]);
              }
              if (stat_col) {
                *cswritep++ = '\t';
                cswritep = dtoa_g(report_stat, cswritep);
              }
              if (df_col) {
                *cswritep++ = '\t';
                cswritep = u32toa(report_df, cswritep);
              }
              if (p_col) {
                *cswritep++ = '\t';
                cswritep = lntoa_g(MAXV(ln_pval, output_min_ln), cswritep);
              }
              AppendBinaryEoln(&cswritep);
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto CalcEpiBoost_ret_WRITE_FAIL;
              }
            }
          }
          if (pairs_seen >= next_print_pair) {
            if (pct > 9) {
              putc_unlocked('\b', stdout);
            }
            pct = (pairs_seen * 100LLU) / pair_ct_total;
            if (pct > 99) {
              pct = 99;
            }
            printf("\b\b%u%%", pct);
            fflush(stdout);
            next_print_pair = ((pct + 1) * pair_ct_total) / 100;
          }
        }
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto CalcEpiBoost_ret_WRITE_FAIL;
    }
    fputs("\b\b\b", stdout);
    logprintf("--epistasis-boost: %" PRIu64 " pair%s tested, %" PRIu64 " written to %s .\n", pairs_seen, (pairs_seen == 1)? "" : "s", pairs_reported, outname);
    if (covar_refit) {
      logprintf("--epistasis-boost: STAT/DF/P adjusted for %u covariate%s; BEST_CHISQ and N_SIG\nare not.\n", cur_covar_ct, (cur_covar_ct == 1)? "" : "s");
      if (refit_fail_ct) {
        logprintfww("Warning: %" PRIu64 " pair%s left out of the report because the covariate-adjusted fit did not converge.\n", refit_fail_ct, (refit_fail_ct == 1)? "" : "s");
      }
    }

    // Summary report: one row per variant with a tested pair.
    char* summary_end = strcpya_k(outname_end2, ".summary");
    if (parallel_tot > 1) {
      *summary_end++ = '.';
      summary_end = u32toa(parallel_idx + 1, summary_end);
    }
    if (output_zst) {
      snprintf(summary_end, kMaxOutfnameExtBlen - S_CAST(uintptr_t, summary_end - outname_end), ".zst");
    } else {
      *summary_end = '\0';
    }
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &csst, &cswritetp);
    if (unlikely(reterr)) {
      goto CalcEpiBoost_ret_1;
    }
    const uint32_t nsig_col = flags & kfEpiColNsig;
    const uint32_t ntot_col = flags & kfEpiColNtot;
    const uint32_t prop_col = (flags & kfEpiColProp) && (parallel_tot == 1);
    *cswritetp++ = '#';
    if (chrom_col) {
      cswritetp = strcpya_k(cswritetp, "CHROM\t");
    }
    if (pos_col) {
      cswritetp = strcpya_k(cswritetp, "POS\t");
    }
    cswritetp = strcpya_k(cswritetp, "ID");
    if (a1_col) {
      cswritetp = strcpya_k(cswritetp, "\tALLELE");
    }
    if (nsig_col) {
      cswritetp = strcpya_k(cswritetp, "\tN_SIG");
    }
    if (ntot_col) {
      cswritetp = strcpya_k(cswritetp, "\tN_TOT");
    }
    if (prop_col) {
      cswritetp = strcpya_k(cswritetp, "\tPROP");
    }
    cswritetp = strcpya_k(cswritetp, "\tBEST_CHISQ");
    if (chrom_col) {
      cswritetp = strcpya_k(cswritetp, "\tBEST_CHROM");
    }
    cswritetp = strcpya_k(cswritetp, "\tBEST_ID" EOLN_STR);
    uint32_t summary_row_ct = 0;
    // after set-by-set / set-by-all implemented, see if this can still be
    // simplified to usual BitIter1() loop with cached chromosome string, etc.
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const EpiSummaryEntry* cur = &(summary[variant_idx]);
      if (!cur->n_tot) {
        continue;
      }
      ++summary_row_ct;
      if (chrom_col) {
        cswritetp = chrtoa(cip, cip->chr_file_order[variant_chr_fo_idxs[variant_idx]], cswritetp);
        *cswritetp++ = '\t';
      }
      const uint32_t variant_uidx = variant_uidxs[variant_idx];
      if (pos_col) {
        cswritetp = u32toa_x(variant_bps[variant_uidx], '\t', cswritetp);
      }
      cswritetp = strcpyax(cswritetp, variant_ids[variant_uidx], '\t');
      if (a1_col) {
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        if (maj_alleles) {
          aidx = maj_alleles[variant_uidx];
        }
        cswritetp = strcpyax(cswritetp, cur_alleles[aidx], '\t');
      }
      if (nsig_col) {
        cswritetp = u32toa_x(cur->n_sig, '\t', cswritetp);
      }
      if (ntot_col) {
        cswritetp = u32toa_x(cur->n_tot, '\t', cswritetp);
      }
      if (prop_col) {
        cswritetp = dtoa_g(S_CAST(double, cur->n_sig) / S_CAST(double, cur->n_tot), cswritetp);
        *cswritetp++ = '\t';
      }
      cswritetp = dtoa_g(cur->best_chisq, cswritetp);
      *cswritetp++ = '\t';
      if (chrom_col) {
        cswritetp = chrtoa(cip, cip->chr_file_order[variant_chr_fo_idxs[cur->best_vidx]], cswritetp);
        *cswritetp++ = '\t';
      }
      cswritetp = strcpya(cswritetp, variant_ids[variant_uidxs[cur->best_vidx]]);
      AppendBinaryEoln(&cswritetp);
      if (unlikely(Cswrite(&csst, &cswritetp))) {
        goto CalcEpiBoost_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&csst, cswritetp))) {
      goto CalcEpiBoost_ret_WRITE_FAIL;
    }
    logprintfww("--epistasis-boost: Summary for %u variant%s written to %s .\n", summary_row_ct, (summary_row_ct == 1)? "" : "s", outname);
  }
  while (0) {
  CalcEpiBoost_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CalcEpiBoost_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  CalcEpiBoost_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  CalcEpiBoost_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  CalcEpiBoost_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 CalcEpiBoost_ret_1:
  CswriteCloseCond(&css, cswritep);
  CswriteCloseCond(&csst, cswritetp);
  BigstackReset(bigstack_mark);
  return reterr;
}

// --epistasis: the interaction term of a linear regression of the phenotype on
// the two genotypes, their product, and any covariates.  PLINK 1.9 supports no
// covariates here, and solves the 4-parameter case with a hardcoded 4x4
// inverse.
//
// Covariates are added the way --glm does them rather than by solving a
// (4 + covariate)-parameter least squares from scratch for every pair: the
// covariate block of X'X is the same for every pair with no missing calls,
// each variant's own row is the same for every pair it appears in, and only
// the product term's row genuinely varies.  So the per-pair work is one pass
// over the samples where both genotypes are nonzero, a correction pass over
// the samples missing either call, and a solve whose dimension does not depend
// on the sample count.
//
// Each variant is held as three bitvectors over the analysis samples: nonzero,
// hom-ALT, and missing.  The genotype is then 0, 1 or 2 without a lookup, and
// the samples that contribute to the product term are one AND away.

// Per-variant summaries, over the samples where this variant is not missing:
// sum of genotypes, sum of squares, dot product with the phenotype, and dot
// product with each covariate.
CONSTI32(kEpiLinearVariantDoubleCt, 3);

static void EpiLinearFillSlot(const uintptr_t* genovec, const double* pheno_vals, const double* covar_vals, uint32_t sample_ct, uint32_t covar_ct, uintptr_t* bits, double* dbls, uint32_t* miss_ct_ptr) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uintptr_t* nonzero = bits;
  uintptr_t* hom_alt = &(bits[sample_ctl]);
  uintptr_t* missing = &(bits[2 * sample_ctl]);
  Halfword* nonzero_alias = R_CAST(Halfword*, nonzero);
  Halfword* hom_alt_alias = R_CAST(Halfword*, hom_alt);
  Halfword* missing_alias = R_CAST(Halfword*, missing);
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    const uintptr_t lo = geno_word & kMask5555;
    const uintptr_t hi = (geno_word >> 1) & kMask5555;
    // 0 = hom-REF, 1 = het, 2 = hom-ALT, 3 = missing.
    nonzero_alias[widx] = PackWordToHalfword(lo ^ hi);
    hom_alt_alias[widx] = PackWordToHalfword(hi & (~lo));
    missing_alias[widx] = PackWordToHalfword(lo & hi);
  }
  if (sample_ctl2 % 2) {
    nonzero_alias[sample_ctl2] = 0;
    hom_alt_alias[sample_ctl2] = 0;
    missing_alias[sample_ctl2] = 0;
  }
  const uint32_t nonzero_ct = PopcountWords(nonzero, sample_ctl);
  const uint32_t hom_alt_ct = PopcountWords(hom_alt, sample_ctl);
  // Each nonzero genotype contributes 1 and each hom-ALT another 1, so the
  // squares are 1 and 4 respectively.
  dbls[0] = u31tod(nonzero_ct + hom_alt_ct);
  dbls[1] = u31tod(nonzero_ct + 3 * hom_alt_ct);
  double pheno_dot = 0.0;
  double* covar_dots = &(dbls[kEpiLinearVariantDoubleCt]);
  for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
    covar_dots[covar_idx] = 0.0;
  }
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    uintptr_t nonzero_word = nonzero[widx];
    if (!nonzero_word) {
      continue;
    }
    const uintptr_t hom_alt_word = hom_alt[widx];
    const uint32_t sample_idx_base = widx * kBitsPerWord;
    do {
      const uint32_t bit_idx = ctzw(nonzero_word);
      nonzero_word &= nonzero_word - 1;
      const uint32_t sample_idx = sample_idx_base + bit_idx;
      const double cur_geno = u31tod(1 + S_CAST(uint32_t, (hom_alt_word >> bit_idx) & 1));
      pheno_dot += cur_geno * pheno_vals[sample_idx];
      for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
        covar_dots[covar_idx] += cur_geno * covar_vals[covar_idx * S_CAST(uintptr_t, sample_ct) + sample_idx];
      }
    } while (nonzero_word);
  }
  dbls[2] = pheno_dot;
  *miss_ct_ptr = PopcountWords(missing, sample_ctl);
}

// Cholesky inverse for the small symmetric positive-definite matrices this
// scan produces.  Only the lower triangle of the input is read, and the full
// inverse is written back.  LAPACK is the right tool at larger dimensions,
// but there is one of these per variant pair, and at this size its per-call
// overhead is most of the cost.
static BoolErr EpiLinearInvertSymmPd(uint32_t dim, double* matrix, double* chol) {
  // Factor: matrix = L L', L lower triangular, into chol.
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    double* chol_row = &(chol[row_idx * S_CAST(uintptr_t, dim)]);
    const double* mat_row = &(matrix[row_idx * S_CAST(uintptr_t, dim)]);
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      const double* chol_col = &(chol[col_idx * S_CAST(uintptr_t, dim)]);
      double cur_sum = mat_row[col_idx];
      for (uint32_t inner_idx = 0; inner_idx != col_idx; ++inner_idx) {
        cur_sum -= chol_row[inner_idx] * chol_col[inner_idx];
      }
      chol_row[col_idx] = cur_sum / chol_col[col_idx];
    }
    double cur_sum = mat_row[row_idx];
    for (uint32_t inner_idx = 0; inner_idx != row_idx; ++inner_idx) {
      cur_sum -= chol_row[inner_idx] * chol_row[inner_idx];
    }
    // A rank-deficient X'X lands on zero or, through rounding, just past it.
    if (!(cur_sum > 1e-12 * fabs(mat_row[row_idx]))) {
      return 1;
    }
    chol_row[row_idx] = sqrt(cur_sum);
  }
  // Invert L in place.
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    double* chol_row = &(chol[row_idx * S_CAST(uintptr_t, dim)]);
    const double cur_recip = 1.0 / chol_row[row_idx];
    chol_row[row_idx] = cur_recip;
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      double cur_sum = 0.0;
      for (uint32_t inner_idx = col_idx; inner_idx != row_idx; ++inner_idx) {
        cur_sum += chol_row[inner_idx] * chol[inner_idx * S_CAST(uintptr_t, dim) + col_idx];
      }
      chol_row[col_idx] = -cur_sum * cur_recip;
    }
  }
  // matrix^{-1} = (L^{-1})' L^{-1}.
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
      double cur_sum = 0.0;
      for (uint32_t inner_idx = row_idx; inner_idx != dim; ++inner_idx) {
        const double* chol_inner = &(chol[inner_idx * S_CAST(uintptr_t, dim)]);
        cur_sum += chol_inner[row_idx] * chol_inner[col_idx];
      }
      matrix[row_idx * S_CAST(uintptr_t, dim) + col_idx] = cur_sum;
      matrix[col_idx * S_CAST(uintptr_t, dim) + row_idx] = cur_sum;
    }
  }
  return 0;
}

typedef struct EpiLinearCtxStruct {
  // Shared, read-only.
  const double* pheno_vals;
  const double* covar_vals;
  const double* base_xtx;
  const double* base_xty;
  double base_pheno_ssq;
  double vif_thresh;
  double max_corr;
  uint32_t sample_ct;
  uint32_t sample_ctl;
  uint32_t covar_ct;
  uint32_t base_dim;
  uint32_t param_ct;
  uint32_t calc_thread_ct;
  uintptr_t words_per_variant;
  uintptr_t doubles_per_variant;

  // Per-thread scratch.
  double** xtxs;
  double** xtys;
  double** beta_bufs;
  double** dbl_2d_bufs;
  double** inv_stdev_bufs;
  double** chol_bufs;
  double** covar_ab_dots;

  // The current work unit: a run of row variants against the loaded column
  // block.  The threads take rows round-robin, since a row's column count
  // falls off across the unit; the main thread then walks the results in row
  // and column order, so the report stays row-major however many threads are
  // running.  Rows are batched rather than dispatched one at a time because a
  // single row is not enough work to cover a thread handoff.
  const uintptr_t* row_bits;
  const double* row_dbls;
  const uint32_t* row_miss_cts;
  const uintptr_t* col_bits;
  const double* col_dbls;
  const uint32_t* col_miss_cts;
  uint32_t row_block_start;
  uint32_t unit_row_slot_first;
  uint32_t unit_row_slot_end;
  uint32_t col_block_start;
  uint32_t col_block_end;
  uint32_t col_block_size;

  // Results, indexed by (row slot within the unit) * col_block_size + column
  // slot.  A zero df means the pair could not be fit.
  double* betas_out;
  double* ses_out;
  double* tstats_out;
  uint32_t* dfs_out;
} EpiLinearCtx;

static void EpiLinearOnePair(EpiLinearCtx* ctx, uint32_t tidx, uint32_t row_slot, uint32_t col_slot, uintptr_t result_idx) {
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t sample_ctl = ctx->sample_ctl;
  const uint32_t cur_covar_ct = ctx->covar_ct;
  const uint32_t param_ct = ctx->param_ct;
  const uint32_t base_dim = ctx->base_dim;
  const uintptr_t covar_stride = sample_ct;
  const double* pheno_vals = ctx->pheno_vals;
  const double* covar_vals = ctx->covar_vals;
  const double* base_xtx = ctx->base_xtx;
  const double* base_xty = ctx->base_xty;
  const double base_pheno_ssq = ctx->base_pheno_ssq;
  const double vif_thresh = ctx->vif_thresh;
  const double max_corr = ctx->max_corr;
  const uintptr_t* row_nonzero = &(ctx->row_bits[row_slot * ctx->words_per_variant]);
  const uintptr_t* row_hom_alt = &(row_nonzero[sample_ctl]);
  const uintptr_t* row_missing = &(row_nonzero[2 * sample_ctl]);
  const double* row_dbl = &(ctx->row_dbls[row_slot * ctx->doubles_per_variant]);
  const uint32_t row_miss_ct = ctx->row_miss_cts[row_slot];
  const uintptr_t* col_nonzero = &(ctx->col_bits[col_slot * ctx->words_per_variant]);
  const uintptr_t* col_hom_alt = &(col_nonzero[sample_ctl]);
  const uintptr_t* col_missing = &(col_nonzero[2 * sample_ctl]);
  const double* col_dbl = &(ctx->col_dbls[col_slot * ctx->doubles_per_variant]);
  const uint32_t col_miss_ct = ctx->col_miss_cts[col_slot];
  double* xtx = ctx->xtxs[tidx];
  double* xty = ctx->xtys[tidx];
  double* betas = ctx->beta_bufs[tidx];
  double* dbl_2d_buf = ctx->dbl_2d_bufs[tidx];
  double* inv_stdevs = ctx->inv_stdev_bufs[tidx];
  double* chol_buf = ctx->chol_bufs[tidx];
  double* covar_ab_dots = cur_covar_ct? ctx->covar_ab_dots[tidx] : nullptr;
  ctx->dfs_out[result_idx] = 0;
  uint32_t joint_cts[4];
  joint_cts[0] = 0;
  joint_cts[1] = 0;
  joint_cts[2] = 0;
  joint_cts[3] = 0;
  double sum_ab_pheno = 0.0;
  for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
    covar_ab_dots[covar_idx] = 0.0;
  }
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    uintptr_t both_word = row_nonzero[widx] & col_nonzero[widx];
    if (!both_word) {
      continue;
    }
    const uintptr_t row_hom_word = row_hom_alt[widx];
    const uintptr_t col_hom_word = col_hom_alt[widx];
    const uint32_t sample_idx_base = widx * kBitsPerWord;
    do {
      const uint32_t bit_idx = ctzw(both_word);
      both_word &= both_word - 1;
      const uint32_t row_is_hom = S_CAST(uint32_t, (row_hom_word >> bit_idx) & 1);
      const uint32_t col_is_hom = S_CAST(uint32_t, (col_hom_word >> bit_idx) & 1);
      joint_cts[row_is_hom * 2 + col_is_hom] += 1;
      const uint32_t sample_idx = sample_idx_base + bit_idx;
      const double cur_ab = u31tod((1 + row_is_hom) * (1 + col_is_hom));
      sum_ab_pheno += cur_ab * pheno_vals[sample_idx];
      for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
        covar_ab_dots[covar_idx] += cur_ab * covar_vals[covar_idx * covar_stride + sample_idx];
      }
    } while (both_word);
  }
  const double sum_ab = u31tod(joint_cts[0] + 2 * (joint_cts[1] + joint_cts[2]) + 4 * joint_cts[3]);
  const double sum_aab = u31tod(joint_cts[0] + 2 * joint_cts[1] + 4 * joint_cts[2] + 8 * joint_cts[3]);
  const double sum_abb = u31tod(joint_cts[0] + 4 * joint_cts[1] + 2 * joint_cts[2] + 8 * joint_cts[3]);
  const double sum_aabb = u31tod(joint_cts[0] + 4 * (joint_cts[1] + joint_cts[2]) + 16 * joint_cts[3]);

  // X'X, lower triangle only, with the predictors ordered intercept,
  // A, B, AB, covariates.
  xtx[0] = u31tod(sample_ct);
  xtx[param_ct] = row_dbl[0];
  xtx[param_ct + 1] = row_dbl[1];
  xtx[2 * param_ct] = col_dbl[0];
  xtx[2 * param_ct + 1] = sum_ab;
  xtx[2 * param_ct + 2] = col_dbl[1];
  xtx[3 * param_ct] = sum_ab;
  xtx[3 * param_ct + 1] = sum_aab;
  xtx[3 * param_ct + 2] = sum_abb;
  xtx[3 * param_ct + 3] = sum_aabb;
  xty[0] = base_xty[0];
  xty[1] = row_dbl[2];
  xty[2] = col_dbl[2];
  xty[3] = sum_ab_pheno;
  for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
    double* cur_row = &(xtx[(4 + covar_idx) * S_CAST(uintptr_t, param_ct)]);
    cur_row[0] = base_xtx[covar_idx + 1];
    cur_row[1] = row_dbl[kEpiLinearVariantDoubleCt + covar_idx];
    cur_row[2] = col_dbl[kEpiLinearVariantDoubleCt + covar_idx];
    cur_row[3] = covar_ab_dots[covar_idx];
    const double* base_row = &(base_xtx[(covar_idx + 1) * S_CAST(uintptr_t, base_dim)]);
    for (uint32_t covar_idx2 = 0; covar_idx2 <= covar_idx; ++covar_idx2) {
      cur_row[4 + covar_idx2] = base_row[covar_idx2 + 1];
    }
    xty[4 + covar_idx] = base_xty[covar_idx + 1];
  }
  double cur_pheno_ssq = base_pheno_ssq;
  uint32_t cur_sample_ct = sample_ct;

  // Everything above except the product term counts every analysis
  // sample; the ones missing either genotype have to come back out.
  if (row_miss_ct || col_miss_ct) {
    for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
      uintptr_t miss_word = row_missing[widx] | col_missing[widx];
      if (!miss_word) {
        continue;
      }
      const uintptr_t row_miss_word = row_missing[widx];
      const uintptr_t col_miss_word = col_missing[widx];
      const uintptr_t row_nonzero_word = row_nonzero[widx];
      const uintptr_t col_nonzero_word = col_nonzero[widx];
      const uintptr_t row_hom_word = row_hom_alt[widx];
      const uintptr_t col_hom_word = col_hom_alt[widx];
      const uint32_t sample_idx_base = widx * kBitsPerWord;
      do {
        const uint32_t bit_idx = ctzw(miss_word);
        miss_word &= miss_word - 1;
        const uint32_t sample_idx = sample_idx_base + bit_idx;
        const double cur_pheno = pheno_vals[sample_idx];
        --cur_sample_ct;
        xtx[0] -= 1.0;
        xty[0] -= cur_pheno;
        cur_pheno_ssq -= cur_pheno * cur_pheno;
        for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
          const double cur_covar = covar_vals[covar_idx * covar_stride + sample_idx];
          double* cur_row = &(xtx[(4 + covar_idx) * S_CAST(uintptr_t, param_ct)]);
          cur_row[0] -= cur_covar;
          xty[4 + covar_idx] -= cur_covar * cur_pheno;
          for (uint32_t covar_idx2 = 0; covar_idx2 <= covar_idx; ++covar_idx2) {
            cur_row[4 + covar_idx2] -= cur_covar * covar_vals[covar_idx2 * covar_stride + sample_idx];
          }
        }
        if (!((row_miss_word >> bit_idx) & 1)) {
          const double cur_geno = u31tod(S_CAST(uint32_t, ((row_nonzero_word >> bit_idx) & 1) + ((row_hom_word >> bit_idx) & 1)));
          if (cur_geno != 0.0) {
            xtx[param_ct] -= cur_geno;
            xtx[param_ct + 1] -= cur_geno * cur_geno;
            xty[1] -= cur_geno * cur_pheno;
            for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
              xtx[(4 + covar_idx) * S_CAST(uintptr_t, param_ct) + 1] -= cur_geno * covar_vals[covar_idx * covar_stride + sample_idx];
            }
          }
        }
        if (!((col_miss_word >> bit_idx) & 1)) {
          const double cur_geno = u31tod(S_CAST(uint32_t, ((col_nonzero_word >> bit_idx) & 1) + ((col_hom_word >> bit_idx) & 1)));
          if (cur_geno != 0.0) {
            xtx[2 * param_ct] -= cur_geno;
            xtx[2 * param_ct + 2] -= cur_geno * cur_geno;
            xty[2] -= cur_geno * cur_pheno;
            for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
              xtx[(4 + covar_idx) * S_CAST(uintptr_t, param_ct) + 2] -= cur_geno * covar_vals[covar_idx * covar_stride + sample_idx];
            }
          }
        }
      } while (miss_word);
    }
  }
  if (cur_sample_ct <= param_ct + 1) {
    return;
  }
  // The same correlation and variance-inflation checks --glm makes, but taken
  // from the cross-products rather than from a separately inverted correlation
  // matrix: the centered cross-product matrix is the Schur complement of the
  // intercept in X'X, so its inverse is the lower right block of (X'X)^{-1},
  // and the variance inflation factors fall out of the inverse the regression
  // needs anyway.  With one of these per variant pair, inverting a second
  // matrix to get them would be most of the arithmetic.
  const uint32_t pred_ct = param_ct - 1;
  const double cur_sample_ct_recip = 1.0 / u31tod(cur_sample_ct);
  for (uint32_t pred_idx = 1; pred_idx != param_ct; ++pred_idx) {
    const double* xtx_row = &(xtx[pred_idx * S_CAST(uintptr_t, param_ct)]);
    double* centered_row = &(dbl_2d_buf[(pred_idx - 1) * S_CAST(uintptr_t, pred_ct)]);
    const double cur_sum = xtx_row[0];
    for (uint32_t pred_idx2 = 1; pred_idx2 <= pred_idx; ++pred_idx2) {
      centered_row[pred_idx2 - 1] = xtx_row[pred_idx2] - cur_sum * xtx[pred_idx2 * S_CAST(uintptr_t, param_ct)] * cur_sample_ct_recip;
    }
  }
  for (uint32_t pred_idx = 0; pred_idx != pred_ct; ++pred_idx) {
    const double cur_var = dbl_2d_buf[pred_idx * S_CAST(uintptr_t, pred_ct) + pred_idx];
    if (!(cur_var > 0.0)) {
      return;
    }
    inv_stdevs[pred_idx] = 1.0 / sqrt(cur_var);
  }
  for (uint32_t pred_idx = 1; pred_idx != pred_ct; ++pred_idx) {
    const double* centered_row = &(dbl_2d_buf[pred_idx * S_CAST(uintptr_t, pred_ct)]);
    const double cur_inv_stdev = inv_stdevs[pred_idx];
    for (uint32_t pred_idx2 = 0; pred_idx2 != pred_idx; ++pred_idx2) {
      if (fabs(centered_row[pred_idx2] * cur_inv_stdev * inv_stdevs[pred_idx2]) > max_corr) {
        return;
      }
    }
  }
  if (EpiLinearInvertSymmPd(param_ct, xtx, chol_buf)) {
    return;
  }
  for (uint32_t pred_idx = 0; pred_idx != pred_ct; ++pred_idx) {
    const double cur_vif = xtx[(pred_idx + 1) * S_CAST(uintptr_t, param_ct + 1)] * dbl_2d_buf[pred_idx * S_CAST(uintptr_t, pred_ct) + pred_idx];
    if (cur_vif > vif_thresh) {
      return;
    }
  }
  double rss = cur_pheno_ssq;
  for (uint32_t pred_idx = 0; pred_idx != param_ct; ++pred_idx) {
    const double* cur_row = &(xtx[pred_idx * S_CAST(uintptr_t, param_ct)]);
    double cur_beta = 0.0;
    for (uint32_t pred_idx2 = 0; pred_idx2 != param_ct; ++pred_idx2) {
      cur_beta += cur_row[pred_idx2] * xty[pred_idx2];
    }
    betas[pred_idx] = cur_beta;
    rss -= cur_beta * xty[pred_idx];
  }
  if (!(rss > 0.0)) {
    return;
  }
  const uint32_t cur_df = cur_sample_ct - param_ct;
  const double sigma_sq = rss / u31tod(cur_df);
  const double se_sq = xtx[3 * S_CAST(uintptr_t, param_ct) + 3] * sigma_sq;
  if (!(se_sq > 0.0)) {
    return;
  }
  const double beta_int = betas[3];
  const double se = sqrt(se_sq);
  const double tstat = beta_int / se;
  if (!isfinite(tstat)) {
    return;
  }
  ctx->betas_out[result_idx] = beta_int;
  ctx->ses_out[result_idx] = se;
  ctx->tstats_out[result_idx] = tstat;
  ctx->dfs_out[result_idx] = cur_df;
}

THREAD_FUNC_DECL EpiLinearThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  EpiLinearCtx* ctx = S_CAST(EpiLinearCtx*, arg->sharedp->context);
  do {
    const uint32_t unit_row_slot_end = ctx->unit_row_slot_end;
    const uint32_t unit_row_slot_first = ctx->unit_row_slot_first;
    const uint32_t calc_thread_ct = ctx->calc_thread_ct;
    const uint32_t row_block_start = ctx->row_block_start;
    const uint32_t col_block_start = ctx->col_block_start;
    const uint32_t col_block_end = ctx->col_block_end;
    const uintptr_t col_block_size = ctx->col_block_size;
    for (uint32_t row_slot = unit_row_slot_first + tidx; row_slot < unit_row_slot_end; row_slot += calc_thread_ct) {
      const uint32_t row_idx = row_block_start + row_slot;
      const uint32_t col_first = MAXV(col_block_start, row_idx + 1);
      if (col_first >= col_block_end) {
        continue;
      }
      const uintptr_t result_base = (row_slot - unit_row_slot_first) * col_block_size;
      for (uint32_t col_idx = col_first; col_idx != col_block_end; ++col_idx) {
        const uint32_t col_slot = col_idx - col_block_start;
        EpiLinearOnePair(ctx, tidx, row_slot, col_slot, result_base + col_slot);
      }
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr CalcEpiLinear(const uintptr_t* orig_sample_include, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const EpiInfo* epi_ip, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_allele_slen, double vif_thresh, double max_corr, double output_min_ln, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  char* cswritetp = nullptr;
  CompressStreamState css;
  CompressStreamState csst;
  ThreadGroup tg;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  PreinitCstream(&csst);
  PreinitThreads(&tg);
  {
    const EpiFlags flags = epi_ip->flags;
    const uint32_t output_zst = (flags / kfEpiZs) & 1;

    // As with --epistasis-boost, PLINK 1.x had a single phenotype, and an
    // O(variant_ct^2) scan is not something to guess about.
    uint32_t pheno_idx = UINT32_MAX;
    uint32_t qt_pheno_ct = 0;
    for (uint32_t uii = 0; uii != pheno_ct; ++uii) {
      if (pheno_cols[uii].type_code == kPhenoDtypeQt) {
        ++qt_pheno_ct;
        if (pheno_idx == UINT32_MAX) {
          pheno_idx = uii;
        }
      }
    }
    if (unlikely(!qt_pheno_ct)) {
      logerrputs("Error: --epistasis requires a quantitative phenotype.  Its case/control\nbranch is not implemented yet; --epistasis-boost covers case/control data.\n");
      goto CalcEpiLinear_ret_INCONSISTENT_INPUT;
    }
    if (unlikely(qt_pheno_ct > 1)) {
      logerrputs("Error: --epistasis needs exactly one quantitative phenotype; select one with\n--pheno-name.\n");
      goto CalcEpiLinear_ret_INCONSISTENT_INPUT;
    }
    const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* sample_include;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &sample_include))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    BitvecAndCopy(orig_sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, sample_include);
    for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
      const PhenoCol* cur_covar_col = &(covar_cols[covar_idx]);
      if (unlikely(cur_covar_col->type_code == kPhenoDtypeCat)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --epistasis does not support categorical covariates yet ('%s'). Split it into binary covariates with --split-cat-pheno first.\n", &(covar_names[covar_idx * max_covar_name_blen]));
        goto CalcEpiLinear_ret_INCONSISTENT_INPUT_WW;
      }
      BitvecAnd(cur_covar_col->nonmiss, raw_sample_ctl, sample_include);
    }
    const uint32_t sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    if (unlikely(sample_ct < 6)) {
      logerrputs("Error: --epistasis needs more samples than regression parameters.\n");
      goto CalcEpiLinear_ret_DEGENERATE_DATA;
    }
    logprintf("--epistasis: Regressing %s on %u sample%s.\n", &(pheno_names[pheno_idx * max_pheno_name_blen]), sample_ct, (sample_ct == 1)? "" : "s");

    // The covariates are stored one column at a time, since the inner loops
    // walk a single covariate across many samples.
    double* pheno_vals;
    double* covar_vals = nullptr;
    if (unlikely(bigstack_alloc_d(sample_ct, &pheno_vals))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    if (covar_ct) {
      if (unlikely(bigstack_alloc_d(S_CAST(uintptr_t, sample_ct) * covar_ct, &covar_vals))) {
        goto CalcEpiLinear_ret_NOMEM;
      }
    }
    {
      uintptr_t sample_uidx_base = 0;
      uintptr_t sample_include_bits = sample_include[0];
      const double* pheno_qt = cur_pheno_col->data.qt;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
        pheno_vals[sample_idx] = pheno_qt[sample_uidx];
        for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
          const PhenoCol* cur_covar_col = &(covar_cols[covar_idx]);
          const double cur_val = (cur_covar_col->type_code == kPhenoDtypeQt)? cur_covar_col->data.qt[sample_uidx] : u31tod(IsSet(cur_covar_col->data.cc, sample_uidx));
          covar_vals[covar_idx * S_CAST(uintptr_t, sample_ct) + sample_idx] = cur_val;
        }
      }
    }

    // A covariate that is constant across the analysis samples is collinear
    // with the intercept, and would fail every pair rather than none.  --glm
    // drops those with a warning; do the same.
    uint32_t cur_covar_ct = 0;
    for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
      const double* cur_col = &(covar_vals[covar_idx * S_CAST(uintptr_t, sample_ct)]);
      const double first_val = cur_col[0];
      uint32_t sample_idx = 1;
      for (; sample_idx != sample_ct; ++sample_idx) {
        if (cur_col[sample_idx] != first_val) {
          break;
        }
      }
      if (sample_idx == sample_ct) {
        logerrprintf("Warning: Excluding constant covariate '%s' from --epistasis.\n", &(covar_names[covar_idx * max_covar_name_blen]));
        continue;
      }
      if (cur_covar_ct != covar_idx) {
        memcpy(&(covar_vals[cur_covar_ct * S_CAST(uintptr_t, sample_ct)]), cur_col, sample_ct * sizeof(double));
      }
      ++cur_covar_ct;
    }
    const uint32_t param_ct = 4 + cur_covar_ct;
    if (unlikely(sample_ct <= param_ct + 1)) {
      logerrputs("Error: --epistasis needs more samples than regression parameters.\n");
      goto CalcEpiLinear_ret_DEGENERATE_DATA;
    }

    // The intercept-and-covariate block of X'X, and the matching part of X'y,
    // over every analysis sample.  A pair only has to correct these for the
    // samples missing one of its two genotypes.
    const uint32_t base_dim = cur_covar_ct + 1;
    double* base_xtx;
    double* base_xty;
    if (unlikely(bigstack_calloc_d(S_CAST(uintptr_t, base_dim) * base_dim, &base_xtx) ||
                 bigstack_calloc_d(base_dim, &base_xty))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    double base_pheno_ssq = 0.0;
    {
      base_xtx[0] = u31tod(sample_ct);
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const double cur_pheno = pheno_vals[sample_idx];
        base_xty[0] += cur_pheno;
        base_pheno_ssq += cur_pheno * cur_pheno;
        for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
          const double cur_covar = covar_vals[covar_idx * S_CAST(uintptr_t, sample_ct) + sample_idx];
          base_xtx[covar_idx + 1] += cur_covar;
          base_xty[covar_idx + 1] += cur_covar * cur_pheno;
          double* cur_row = &(base_xtx[(covar_idx + 1) * S_CAST(uintptr_t, base_dim)]);
          for (uint32_t covar_idx2 = 0; covar_idx2 <= covar_idx; ++covar_idx2) {
            cur_row[covar_idx2 + 1] += cur_covar * covar_vals[covar_idx2 * S_CAST(uintptr_t, sample_ct) + sample_idx];
          }
        }
      }
      // fill in the upper triangle, and the intercept column
      for (uint32_t row_idx = 1; row_idx != base_dim; ++row_idx) {
        base_xtx[row_idx * S_CAST(uintptr_t, base_dim)] = base_xtx[row_idx];
        for (uint32_t col_idx = row_idx + 1; col_idx != base_dim; ++col_idx) {
          base_xtx[row_idx * S_CAST(uintptr_t, base_dim) + col_idx] = base_xtx[col_idx * S_CAST(uintptr_t, base_dim) + row_idx];
        }
      }
    }

    // Non-autosomal variants are left out, as in PLINK 1.x.
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* variant_include;
    if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    memcpy(variant_include, orig_variant_include, raw_variant_ctl * sizeof(intptr_t));
    const uint32_t chr_ct = cip->chr_ct;
    uint32_t nonautosomal_ct = 0;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if (!IsSet(cip->haploid_mask, chr_idx) && (chr_idx != S_CAST(uint32_t, cip->xymt_codes[kChrOffsetMT]))) {
        continue;
      }
      const uint32_t chr_vidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
      const uint32_t chr_vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      nonautosomal_ct += PopcountBitRange(variant_include, chr_vidx_start, chr_vidx_end);
      ClearBitsNz(chr_vidx_start, chr_vidx_end, variant_include);
    }
    uint32_t variant_ct = orig_variant_ct - nonautosomal_ct;
    if (unlikely(variant_ct < 2)) {
      logerrputs("Error: --epistasis requires at least 2 autosomal variants.\n");
      goto CalcEpiLinear_ret_DEGENERATE_DATA;
    }

    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* genovec;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_w(NypCtToWordCt(sample_ct), &genovec))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);

    // A variant with a single genotype value across the analysis samples makes
    // the regression singular in every pair it appears in, so it is dropped up
    // front rather than failing variant_ct times.
    {
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      uint32_t skipped_ct = 0;
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        reterr = PgrGet(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, genovec);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, variant_uidx);
          goto CalcEpiLinear_ret_1;
        }
        ZeroTrailingNyps(sample_ct, genovec);
        STD_ARRAY_DECL(uint32_t, 4, genocounts);
        GenoarrCountFreqsUnsafe(genovec, sample_ct, genocounts);
        const uint32_t nonmiss_ct = genocounts[0] + genocounts[1] + genocounts[2];
        if ((nonmiss_ct <= param_ct) || (genocounts[0] == nonmiss_ct) || (genocounts[1] == nonmiss_ct) || (genocounts[2] == nonmiss_ct)) {
          ClearBit(variant_uidx, variant_include);
          ++skipped_ct;
        }
      }
      if (skipped_ct) {
        variant_ct -= skipped_ct;
        logprintf("--epistasis: Skipping %u monomorphic variant%s.\n", skipped_ct, (skipped_ct == 1)? "" : "s");
      }
      if (unlikely(variant_ct < 2)) {
        logerrputs("Error: --epistasis has fewer than 2 usable variants left.\n");
        goto CalcEpiLinear_ret_DEGENERATE_DATA;
      }
    }
    if (nonautosomal_ct) {
      logprintf("--epistasis: Skipping %u non-autosomal variant%s.\n", nonautosomal_ct, (nonautosomal_ct == 1)? "" : "s");
    }

    uint32_t* variant_uidxs;
    uint32_t* variant_chr_fo_idxs;
    if (unlikely(bigstack_alloc_u32(variant_ct, &variant_uidxs) ||
                 bigstack_alloc_u32(variant_ct, &variant_chr_fo_idxs))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    {
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      uint32_t chr_fo_idx = 0;
      uint32_t chr_vidx_end = cip->chr_fo_vidx_start[1];
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        while (variant_uidx >= chr_vidx_end) {
          ++chr_fo_idx;
          chr_vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        }
        variant_uidxs[variant_idx] = variant_uidx;
        variant_chr_fo_idxs[variant_idx] = chr_fo_idx;
      }
    }

    // Rows are split across --parallel jobs, mirrored the same way as in
    // CalcEpiBoost(): ParallelBounds() splits a triangle whose row r carries r
    // entries below it, and this scan's row r carries variant_ct - 1 - r
    // entries above it.
    uint32_t mirror_start;
    uint32_t mirror_end;
    ParallelBounds(variant_ct, 1, parallel_tot - 1 - parallel_idx, parallel_tot, R_CAST(int32_t*, &mirror_start), R_CAST(int32_t*, &mirror_end));
    const uint32_t row_start_idx = variant_ct - mirror_end;
    const uint32_t row_end_idx = variant_ct - mirror_start;
    if (row_start_idx == row_end_idx) {
      logerrputs("Warning: This --parallel job has no rows to scan.\n");
    }

    // The pairs are independent, so the threads take rows round-robin, and
    // each needs its own solve scratch.
    uint32_t calc_thread_ct = MINV(max_thread_ct, variant_ct - 1);
    if (!calc_thread_ct) {
      calc_thread_ct = 1;
    }
    const uintptr_t param_ct_x2 = S_CAST(uintptr_t, param_ct) * param_ct;
    EpiLinearCtx ctx;
    if (unlikely(bigstack_alloc_dp(calc_thread_ct, &ctx.xtxs) ||
                 bigstack_alloc_dp(calc_thread_ct, &ctx.xtys) ||
                 bigstack_alloc_dp(calc_thread_ct, &ctx.beta_bufs) ||
                 bigstack_alloc_dp(calc_thread_ct, &ctx.dbl_2d_bufs) ||
                 bigstack_alloc_dp(calc_thread_ct, &ctx.inv_stdev_bufs) ||
                 bigstack_alloc_dp(calc_thread_ct, &ctx.chol_bufs) ||
                 bigstack_alloc_dp(calc_thread_ct, &ctx.covar_ab_dots))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      if (unlikely(bigstack_alloc_d(param_ct_x2, &(ctx.xtxs[tidx])) ||
                   bigstack_alloc_d(param_ct, &(ctx.xtys[tidx])) ||
                   bigstack_alloc_d(param_ct, &(ctx.beta_bufs[tidx])) ||
                   bigstack_alloc_d(param_ct_x2, &(ctx.dbl_2d_bufs[tidx])) ||
                   bigstack_alloc_d(param_ct, &(ctx.inv_stdev_bufs[tidx])) ||
                   bigstack_alloc_d(param_ct_x2, &(ctx.chol_bufs[tidx])))) {
        goto CalcEpiLinear_ret_NOMEM;
      }
      ctx.covar_ab_dots[tidx] = nullptr;
      if (cur_covar_ct) {
        if (unlikely(bigstack_alloc_d(cur_covar_ct, &(ctx.covar_ab_dots[tidx])))) {
          goto CalcEpiLinear_ret_NOMEM;
        }
      }
    }

    // Two variant blocks are held at once, rows and columns, as in
    // CalcEpiBoost().
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uintptr_t words_per_variant = 3 * S_CAST(uintptr_t, sample_ctl);
    const uintptr_t doubles_per_variant = kEpiLinearVariantDoubleCt + cur_covar_ct;
    const uintptr_t bytes_per_variant = words_per_variant * sizeof(intptr_t) + doubles_per_variant * sizeof(double) + sizeof(int32_t);
    uintptr_t max_slot_ct = (bigstack_left() / 2) / bytes_per_variant;
    if (unlikely(max_slot_ct < 4)) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    // The report has to come out in row-major order, so that concatenating
    // --parallel jobs reproduces a single run.  A row block only preserves
    // that while the whole column range fits in one block; when it does not,
    // the row block drops to a single row, so its columns are still swept in
    // order.  That costs a reread of the column range per row, but only in the
    // case that was already going to be dominated by rereads.
    uint32_t col_block_size;
    uint32_t row_block_size;
    if (max_slot_ct > variant_ct) {
      col_block_size = variant_ct;
      row_block_size = MINV(max_slot_ct - variant_ct, variant_ct);
    } else {
      col_block_size = max_slot_ct - 1;
      row_block_size = 1;
    }
    uintptr_t* row_bits;
    uintptr_t* col_bits;
    double* row_dbls;
    double* col_dbls;
    uint32_t* row_miss_cts;
    uint32_t* col_miss_cts;
    if (unlikely(bigstack_alloc_w(row_block_size * words_per_variant, &row_bits) ||
                 bigstack_alloc_w(col_block_size * words_per_variant, &col_bits) ||
                 bigstack_alloc_d(row_block_size * doubles_per_variant, &row_dbls) ||
                 bigstack_alloc_d(col_block_size * doubles_per_variant, &col_dbls) ||
                 bigstack_alloc_u32(row_block_size, &row_miss_cts) ||
                 bigstack_alloc_u32(col_block_size, &col_miss_cts))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    // A batch of rows' results, which the threads fill in any order and the
    // main thread then walks in row and column order.  One row is not enough
    // work to cover a thread handoff, so rows go out in batches sized to a
    // memory budget.
    uint32_t rows_per_unit = (16 * 1024 * 1024) / (col_block_size * (3 * sizeof(double) + sizeof(int32_t)));
    if (!rows_per_unit) {
      rows_per_unit = 1;
    }
    if (rows_per_unit > row_block_size) {
      rows_per_unit = row_block_size;
    }
    const uintptr_t result_ct = S_CAST(uintptr_t, rows_per_unit) * col_block_size;
    if (unlikely(bigstack_alloc_d(result_ct, &ctx.betas_out) ||
                 bigstack_alloc_d(result_ct, &ctx.ses_out) ||
                 bigstack_alloc_d(result_ct, &ctx.tstats_out) ||
                 bigstack_alloc_u32(result_ct, &ctx.dfs_out))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    ctx.pheno_vals = pheno_vals;
    ctx.covar_vals = covar_vals;
    ctx.base_xtx = base_xtx;
    ctx.base_xty = base_xty;
    ctx.base_pheno_ssq = base_pheno_ssq;
    ctx.vif_thresh = vif_thresh;
    ctx.max_corr = max_corr;
    ctx.sample_ct = sample_ct;
    ctx.sample_ctl = sample_ctl;
    ctx.covar_ct = cur_covar_ct;
    ctx.base_dim = base_dim;
    ctx.param_ct = param_ct;
    ctx.calc_thread_ct = calc_thread_ct;
    ctx.words_per_variant = words_per_variant;
    ctx.doubles_per_variant = doubles_per_variant;
    ctx.row_bits = row_bits;
    ctx.row_dbls = row_dbls;
    ctx.row_miss_cts = row_miss_cts;
    ctx.col_bits = col_bits;
    ctx.col_dbls = col_dbls;
    ctx.col_miss_cts = col_miss_cts;
    ctx.col_block_size = col_block_size;
    ctx.row_block_start = 0;
    ctx.unit_row_slot_first = 0;
    ctx.unit_row_slot_end = 0;
    ctx.col_block_start = 0;
    ctx.col_block_end = 0;
    SetThreadFuncAndData(EpiLinearThread, &ctx, &tg);
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto CalcEpiLinear_ret_NOMEM;
    }

    EpiSummaryEntry* summary;
    if (unlikely(BIGSTACK_ALLOC_X(EpiSummaryEntry, variant_ct, &summary))) {
      goto CalcEpiLinear_ret_NOMEM;
    }
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      summary[variant_idx].n_sig = 0;
      summary[variant_idx].n_tot = 0;
      summary[variant_idx].best_chisq = -1.0;
      summary[variant_idx].best_vidx = UINT32_MAX;
    }

    double alpha1_ln = epi_ip->ln_epi1;
    if (alpha1_ln > 0.0) {
      alpha1_ln = -4 * kLn10;
    }
    const double alpha2_ln = epi_ip->ln_epi2;
    char* outname_end2 = strcpya_k(outname_end, ".epi.qt");
    char* main_end = outname_end2;
    if (parallel_tot > 1) {
      *main_end++ = '.';
      main_end = u32toa(parallel_idx + 1, main_end);
    }
    if (output_zst) {
      snprintf(main_end, kMaxOutfnameExtBlen - S_CAST(uintptr_t, main_end - outname_end), ".zst");
    } else {
      *main_end = '\0';
    }
    const uint32_t chrom_col = flags & kfEpiColChrom;
    const uint32_t pos_col = flags & kfEpiColPos;
    const uint32_t a1_col = (flags & kfEpiColA1) || ((flags & kfEpiColMaybeA1) && MultiallelicVariantPresent(variant_include, allele_idx_offsets, variant_ct));
    const uint32_t beta_col = flags & (kfEpiColBeta | kfEpiColOrbeta);
    const uint32_t se_col = flags & kfEpiColSe;
    const uint32_t stat_col = flags & kfEpiColStat;
    const uint32_t p_col = flags & kfEpiColP;
    uint32_t max_chr_slen = 0;
    if (chrom_col) {
      max_chr_slen = GetMaxChrSlen(cip);
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + 2 * max_chr_slen + 2 * a1_col * max_allele_slen + 2 * kMaxIdSlen + 256;
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto CalcEpiLinear_ret_1;
    }
    if (!parallel_idx) {
      *cswritep++ = '#';
      if (chrom_col) {
        cswritep = strcpya_k(cswritep, "CHROM1\t");
      }
      if (pos_col) {
        cswritep = strcpya_k(cswritep, "POS1\t");
      }
      cswritep = strcpya_k(cswritep, "ID1");
      if (a1_col) {
        cswritep = strcpya_k(cswritep, "\tALLELE1");
      }
      if (chrom_col) {
        cswritep = strcpya_k(cswritep, "\tCHROM2");
      }
      if (pos_col) {
        cswritep = strcpya_k(cswritep, "\tPOS2");
      }
      cswritep = strcpya_k(cswritep, "\tID2");
      if (a1_col) {
        cswritep = strcpya_k(cswritep, "\tALLELE2");
      }
      if (beta_col) {
        cswritep = strcpya_k(cswritep, "\tBETA_INT");
      }
      if (se_col) {
        cswritep = strcpya_k(cswritep, "\tSE");
      }
      if (stat_col) {
        cswritep = strcpya_k(cswritep, "\tT_STAT");
      }
      if (p_col) {
        cswritep = strcpya_k(cswritep, "\tP");
      }
      AppendBinaryEoln(&cswritep);
    }

    uint64_t pair_ct_total = 0;
    for (uint32_t row_idx = row_start_idx; row_idx != row_end_idx; ++row_idx) {
      pair_ct_total += variant_ct - row_idx - 1;
    }
    uint64_t pairs_reported = 0;
    uint64_t pairs_seen = 0;
    uint64_t pairs_tested = 0;
    fputs("--epistasis: 0%", stdout);
    fflush(stdout);
    uint64_t next_print_pair = pair_ct_total / 100;
    uint32_t pct = 0;

    if (flags & kfEpiRefBased) {
      maj_alleles = nullptr;
    }
    AlleleCode aidx = 0;
    for (uint32_t row_block_start = row_start_idx; row_block_start < row_end_idx; row_block_start += row_block_size) {
      const uint32_t row_block_end = MINV(row_block_start + row_block_size, row_end_idx);
      const uint32_t cur_row_ct = row_block_end - row_block_start;
      for (uint32_t slot_idx = 0; slot_idx != cur_row_ct; ++slot_idx) {
        const uint32_t variant_uidx = variant_uidxs[row_block_start + slot_idx];
        if (maj_alleles) {
          aidx = maj_alleles[variant_uidx];
        }
        reterr = PgrGetInv1(sample_include, pssi, sample_ct, variant_uidx, aidx, simple_pgrp, genovec);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, variant_uidx);
          goto CalcEpiLinear_ret_1;
        }
        ZeroTrailingNyps(sample_ct, genovec);
        EpiLinearFillSlot(genovec, pheno_vals, covar_vals, sample_ct, cur_covar_ct, &(row_bits[slot_idx * words_per_variant]), &(row_dbls[slot_idx * doubles_per_variant]), &(row_miss_cts[slot_idx]));
      }
      for (uint32_t col_block_start = row_block_start; col_block_start < variant_ct; col_block_start += col_block_size) {
        const uint32_t col_block_end = MINV(col_block_start + col_block_size, variant_ct);
        const uint32_t cur_col_ct = col_block_end - col_block_start;
        for (uint32_t slot_idx = 0; slot_idx != cur_col_ct; ++slot_idx) {
          const uint32_t variant_uidx = variant_uidxs[col_block_start + slot_idx];
          if (maj_alleles) {
            aidx = maj_alleles[variant_uidx];
          }
          reterr = PgrGetInv1(sample_include, pssi, sample_ct, variant_uidx, aidx, simple_pgrp, genovec);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto CalcEpiLinear_ret_1;
          }
          ZeroTrailingNyps(sample_ct, genovec);
          EpiLinearFillSlot(genovec, pheno_vals, covar_vals, sample_ct, cur_covar_ct, &(col_bits[slot_idx * words_per_variant]), &(col_dbls[slot_idx * doubles_per_variant]), &(col_miss_cts[slot_idx]));
        }
        ctx.row_block_start = row_block_start;
        ctx.col_block_start = col_block_start;
        ctx.col_block_end = col_block_end;
        for (uint32_t unit_start = 0; unit_start < cur_row_ct; unit_start += rows_per_unit) {
          const uint32_t unit_end = MINV(unit_start + rows_per_unit, cur_row_ct);
          ctx.unit_row_slot_first = unit_start;
          ctx.unit_row_slot_end = unit_end;
          if (unlikely(SpawnThreads(&tg))) {
            goto CalcEpiLinear_ret_THREAD_CREATE_FAIL;
          }
          JoinThreads(&tg);
          for (uint32_t row_slot = unit_start; row_slot != unit_end; ++row_slot) {
            const uint32_t row_idx = row_block_start + row_slot;
            const uint32_t col_first = MAXV(col_block_start, row_idx + 1);
            if (col_first >= col_block_end) {
              continue;
            }
            const uintptr_t result_base = S_CAST(uintptr_t, row_slot - unit_start) * col_block_size;
            for (uint32_t col_idx = col_first; col_idx != col_block_end; ++col_idx) {
              ++pairs_seen;
              const uintptr_t result_idx = result_base + (col_idx - col_block_start);
              const uint32_t cur_df = ctx.dfs_out[result_idx];
              if (!cur_df) {
                continue;
              }
              ++pairs_tested;
              const double beta_int = ctx.betas_out[result_idx];
              const double se = ctx.ses_out[result_idx];
              const double tstat = ctx.tstats_out[result_idx];
              const double ln_pval = TstatToLnP(tstat, cur_df);
              // The summary's BEST_CHISQ is the squared t-statistic, which is
              // what PLINK 1.9 reports there.
              const double chisq = tstat * tstat;
              summary[row_idx].n_tot += 1;
              summary[col_idx].n_tot += 1;
              if (ln_pval <= alpha2_ln) {
                summary[row_idx].n_sig += 1;
                summary[col_idx].n_sig += 1;
              }
              if (chisq > summary[row_idx].best_chisq) {
                summary[row_idx].best_chisq = chisq;
                summary[row_idx].best_vidx = col_idx;
              }
              if (chisq > summary[col_idx].best_chisq) {
                summary[col_idx].best_chisq = chisq;
                summary[col_idx].best_vidx = row_idx;
              }
              if (ln_pval <= alpha1_ln) {
                ++pairs_reported;
                if (chrom_col) {
                  cswritep = chrtoa(cip, cip->chr_file_order[variant_chr_fo_idxs[row_idx]], cswritep);
                  *cswritep++ = '\t';
                }
                const uint32_t row_variant_uidx = variant_uidxs[row_idx];
                if (pos_col) {
                  cswritep = u32toa_x(variant_bps[row_variant_uidx], '\t', cswritep);
                }
                cswritep = strcpyax(cswritep, variant_ids[row_variant_uidx], '\t');
                if (a1_col) {
                  uintptr_t allele_idx_offset_base = row_variant_uidx * 2;
                  if (allele_idx_offsets) {
                    allele_idx_offset_base = allele_idx_offsets[row_variant_uidx];
                  }
                  const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
                  if (maj_alleles) {
                    aidx = maj_alleles[row_variant_uidx];
                  }
                  cswritep = strcpyax(cswritep, cur_alleles[aidx], '\t');
                }
                if (chrom_col) {
                  cswritep = chrtoa(cip, cip->chr_file_order[variant_chr_fo_idxs[col_idx]], cswritep);
                  *cswritep++ = '\t';
                }
                const uint32_t col_variant_uidx = variant_uidxs[col_idx];
                if (pos_col) {
                  cswritep = u32toa_x(variant_bps[col_variant_uidx], '\t', cswritep);
                }
                cswritep = strcpyax(cswritep, variant_ids[col_variant_uidx], '\t');
                if (a1_col) {
                  uintptr_t allele_idx_offset_base = col_variant_uidx * 2;
                  if (allele_idx_offsets) {
                    allele_idx_offset_base = allele_idx_offsets[col_variant_uidx];
                  }
                  const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
                  if (maj_alleles) {
                    aidx = maj_alleles[col_variant_uidx];
                  }
                  cswritep = strcpya(cswritep, cur_alleles[aidx]);
                }
                if (beta_col) {
                  *cswritep++ = '\t';
                  cswritep = dtoa_g(beta_int, cswritep);
                }
                if (se_col) {
                  *cswritep++ = '\t';
                  cswritep = dtoa_g(se, cswritep);
                }
                if (stat_col) {
                  *cswritep++ = '\t';
                  cswritep = dtoa_g(tstat, cswritep);
                }
                if (p_col) {
                  *cswritep++ = '\t';
                  cswritep = lntoa_g(MAXV(ln_pval, output_min_ln), cswritep);
                }
                AppendBinaryEoln(&cswritep);
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto CalcEpiLinear_ret_WRITE_FAIL;
                }
              }
            }
          }
          if (pairs_seen >= next_print_pair) {
            if (pct > 9) {
              putc_unlocked('\b', stdout);
            }
            pct = (pairs_seen * 100LLU) / pair_ct_total;
            if (pct > 99) {
              pct = 99;
            }
            printf("\b\b%u%%", pct);
            fflush(stdout);
            next_print_pair = ((pct + 1) * pair_ct_total) / 100;
          }
        }
      }
    }
    // The thread group needs a last block to shut down on, and the scan does
    // not know which of its units is last until it is past it.
    ctx.unit_row_slot_first = 0;
    ctx.unit_row_slot_end = 0;
    DeclareLastThreadBlock(&tg);
    if (unlikely(SpawnThreads(&tg))) {
      goto CalcEpiLinear_ret_THREAD_CREATE_FAIL;
    }
    JoinThreads(&tg);
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto CalcEpiLinear_ret_WRITE_FAIL;
    }
    fputs("\b\b\b", stdout);
    logprintf("--epistasis: %" PRIu64 " pair%s tested, %" PRIu64 " written to %s .\n", pairs_tested, (pairs_tested == 1)? "" : "s", pairs_reported, outname);
    if (pairs_tested != pairs_seen) {
      const uint64_t fail_ct = pairs_seen - pairs_tested;
      logprintf("--epistasis: %" PRIu64 " pair%s skipped as singular or rank-deficient.\n", fail_ct, (fail_ct == 1)? "" : "s");
    }

    char* summary_end = strcpya_k(outname_end2, ".summary");
    if (parallel_tot > 1) {
      *summary_end++ = '.';
      summary_end = u32toa(parallel_idx + 1, summary_end);
    }
    if (output_zst) {
      snprintf(summary_end, kMaxOutfnameExtBlen - S_CAST(uintptr_t, summary_end - outname_end), ".zst");
    } else {
      *summary_end = '\0';
    }
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &csst, &cswritetp);
    if (unlikely(reterr)) {
      goto CalcEpiLinear_ret_1;
    }
    const uint32_t nsig_col = flags & kfEpiColNsig;
    const uint32_t ntot_col = flags & kfEpiColNtot;
    const uint32_t prop_col = (flags & kfEpiColProp) && (parallel_tot == 1);
    *cswritetp++ = '#';
    if (chrom_col) {
      cswritetp = strcpya_k(cswritetp, "CHROM\t");
    }
    if (pos_col) {
      cswritetp = strcpya_k(cswritetp, "POS\t");
    }
    cswritetp = strcpya_k(cswritetp, "ID");
    if (a1_col) {
      cswritetp = strcpya_k(cswritetp, "\tALLELE");
    }
    if (nsig_col) {
      cswritetp = strcpya_k(cswritetp, "\tN_SIG");
    }
    if (ntot_col) {
      cswritetp = strcpya_k(cswritetp, "\tN_TOT");
    }
    if (prop_col) {
      cswritetp = strcpya_k(cswritetp, "\tPROP");
    }
    cswritetp = strcpya_k(cswritetp, "\tBEST_CHISQ");
    if (chrom_col) {
      cswritetp = strcpya_k(cswritetp, "\tBEST_CHROM");
    }
    cswritetp = strcpya_k(cswritetp, "\tBEST_ID" EOLN_STR);
    uint32_t summary_row_ct = 0;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const EpiSummaryEntry* cur = &(summary[variant_idx]);
      if (!cur->n_tot) {
        continue;
      }
      ++summary_row_ct;
      if (chrom_col) {
        cswritetp = chrtoa(cip, cip->chr_file_order[variant_chr_fo_idxs[variant_idx]], cswritetp);
        *cswritetp++ = '\t';
      }
      const uint32_t variant_uidx = variant_uidxs[variant_idx];
      if (pos_col) {
        cswritetp = u32toa_x(variant_bps[variant_uidx], '\t', cswritetp);
      }
      cswritetp = strcpyax(cswritetp, variant_ids[variant_uidx], '\t');
      if (a1_col) {
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        if (maj_alleles) {
          aidx = maj_alleles[variant_uidx];
        }
        cswritetp = strcpyax(cswritetp, cur_alleles[aidx], '\t');
      }
      if (nsig_col) {
        cswritetp = u32toa_x(cur->n_sig, '\t', cswritetp);
      }
      if (ntot_col) {
        cswritetp = u32toa_x(cur->n_tot, '\t', cswritetp);
      }
      if (prop_col) {
        cswritetp = dtoa_g(S_CAST(double, cur->n_sig) / S_CAST(double, cur->n_tot), cswritetp);
        *cswritetp++ = '\t';
      }
      cswritetp = dtoa_g(cur->best_chisq, cswritetp);
      *cswritetp++ = '\t';
      if (chrom_col) {
        cswritetp = chrtoa(cip, cip->chr_file_order[variant_chr_fo_idxs[cur->best_vidx]], cswritetp);
        *cswritetp++ = '\t';
      }
      cswritetp = strcpya(cswritetp, variant_ids[variant_uidxs[cur->best_vidx]]);
      AppendBinaryEoln(&cswritetp);
      if (unlikely(Cswrite(&csst, &cswritetp))) {
        goto CalcEpiLinear_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&csst, cswritetp))) {
      goto CalcEpiLinear_ret_WRITE_FAIL;
    }
    logprintfww("--epistasis: Summary for %u variant%s written to %s .\n", summary_row_ct, (summary_row_ct == 1)? "" : "s", outname);
  }
  while (0) {
  CalcEpiLinear_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CalcEpiLinear_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  CalcEpiLinear_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  CalcEpiLinear_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  CalcEpiLinear_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  CalcEpiLinear_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 CalcEpiLinear_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&css, cswritep);
  CswriteCloseCond(&csst, cswritetp);
  BigstackReset(bigstack_mark);
  return reterr;
}


#ifdef __cplusplus
}  // namespace plink2
#endif
