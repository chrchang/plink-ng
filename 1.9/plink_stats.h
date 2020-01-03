#ifndef __PLINK_STATS_H__
#define __PLINK_STATS_H__

// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink_matrix.h"

// result of inverse_chiprob(5e-324, 1)
#define MAX_INVERSE_CHIPROB_1DF 1957.4999902125001

// NOT THREAD-SAFE.
double chiprob_p(double xx, double df);

// NOT THREAD-SAFE.
static inline double chiprob_px(double xx, double df) {
  if (xx != -9) {
    return chiprob_p(xx, df);
  } else {
    return -9;
  }
}

// NOT THREAD-SAFE.
double inverse_chiprob(double qq, double df);

// NOT THREAD-SAFE.
double calc_tprob(double tt, double df);

// NOT THREAD-SAFE.
double inverse_tprob(double dbl_qq, double df);

double ltqnorm(double p);

double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp);

int32_t SNPHWE_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);

int32_t SNPHWE_midp_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t midp);

double fisher22_tail_pval(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, int32_t right_offset, double tot_prob, double right_prob, uint32_t midp, uint32_t new_m11);

void fisher22_precomp_pval_bounds(double pval, uint32_t midp, uint32_t row1_sum, uint32_t col1_sum, uint32_t total, uint32_t* bounds, double* tprobs);

double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23, uint32_t midp);

double chi22_eval(intptr_t m11, intptr_t row1_sum, intptr_t col1_sum, intptr_t total);

double chi22_evalx(intptr_t m11, intptr_t row1_sum, intptr_t col1_sum, intptr_t total);

void chi22_precomp_val_bounds(double chisq, intptr_t row1_sum, intptr_t col1_sum, intptr_t total, uint32_t* bounds, double* coeffs);

double chi23_eval(intptr_t m11, intptr_t m12, intptr_t row1_sum, intptr_t col1_sum, intptr_t col2_sum, intptr_t total);

void chi23_evalx(intptr_t m11, intptr_t m12, intptr_t m13, intptr_t m21, intptr_t m22, intptr_t m23, double* chip, uint32_t* dfp);

double ca_trend_eval(intptr_t case_dom_ct, intptr_t case_ct, intptr_t het_ct, intptr_t homdom_ct, intptr_t total);

double ca_trend_evalx(intptr_t case_dom_ct, intptr_t case_ct, intptr_t het_ct, intptr_t homdom_ct, intptr_t total);

void ca_trend_precomp_val_bounds(double chisq, intptr_t case_ct, intptr_t het_ct, intptr_t homdom_ct, intptr_t total, uint32_t* bounds, double* coeffs);

uint32_t linear_hypothesis_chisq(uintptr_t constraint_ct, uintptr_t param_ct, double* constraints_con_major, double* coefs, double* cov_matrix, double* param_df_buf, double* param_df_buf2, double* df_df_buf, MATRIX_INVERT_BUF1_TYPE* mi_buf, double* df_buf, double* chisq_ptr);

double binom_2sided(uint32_t succ, uint32_t obs, uint32_t midp);

#endif // __PLINK_STATS_H__
