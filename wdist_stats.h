#ifndef __WDIST_STATS_H__
#define __WDIST_STATS_H__

#include "wdist_common.h"

double chiprob_p(double xx, double df);

double inverse_chiprob(double qq, double df);

double ltqnorm(double p);

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22);

double fisher22_tail_pval(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t right_offset, double tot_prob, double right_prob, double tail_sum, uint32_t new_m11);

void fisher22_precomp_pval_bounds(double pval, uint32_t row1_sum, uint32_t col1_sum, uint32_t total, uint32_t* bounds, double* tprobs);

double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23);

void chi22_get_coeffs(double row1_sum, double col1_sum, double total, double* expm11p, double* recip_sump);

void chi22_precomp_val_bounds(double chisq, intptr_t row1_sum, intptr_t col1_sum, intptr_t total, uint32_t* bounds);

#endif
