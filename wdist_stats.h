#ifndef __WDIST_STATS_H__
#define __WDIST_STATS_H__

#include "wdist_common.h"

double chiprob_p(double xx, double df);

static inline double chiprob_px(double xx, double df) {
  if (xx != -9) {
    return chiprob_p(xx, df);
  } else {
    return -9;
  }
}

double inverse_chiprob(double qq, double df);

double ltqnorm(double p);

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22);

double fisher22_tail_pval(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t right_offset, double tot_prob, double right_prob, double tail_sum, uint32_t new_m11);

void fisher22_precomp_pval_bounds(double pval, uint32_t row1_sum, uint32_t col1_sum, uint32_t total, uint32_t* bounds, double* tprobs);

double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23);

double chi22_eval(intptr_t m11, intptr_t row1_sum, intptr_t col1_sum, intptr_t total);

double chi22_evalx(intptr_t m11, intptr_t row1_sum, intptr_t col1_sum, intptr_t total);

void chi22_precomp_val_bounds(double chisq, intptr_t row1_sum, intptr_t col1_sum, intptr_t total, uint32_t* bounds, double* coeffs);

double chi23_eval(intptr_t m11, intptr_t m12, intptr_t row1_sum, intptr_t col1_sum, intptr_t col2_sum, intptr_t total);

void chi23_evalx(intptr_t m11, intptr_t m12, intptr_t m13, intptr_t m21, intptr_t m22, intptr_t m23, double* chip, uint32_t* dfp);

double ca_trend_eval(intptr_t case_a2_ct, intptr_t case_ct, intptr_t het_ct, intptr_t homa2_ct, intptr_t total);

double ca_trend_evalx(intptr_t case_a2_ct, intptr_t case_ct, intptr_t het_ct, intptr_t homa2_ct, intptr_t total);

void ca_trend_precomp_val_bounds(double chisq, intptr_t case_ct, intptr_t het_ct, intptr_t homa2_ct, intptr_t total, uint32_t* bounds, double* coeffs);

#endif __WDIST_STATS_H__
