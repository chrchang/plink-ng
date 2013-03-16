#ifndef __WDIST_STATS_H__
#define __WDIST_STATS_H__

#include "wdist_common.h"

double chiprob_p(double xx, double df);

double inverse_chiprob(double qq, double df);

double ltqnorm(double p);

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22);

void fisher22_precomp_thresh(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t* m11_minp, uint32_t* m11_maxp, uint32_t* tiel, uint32_t* tieh);

double fisher22_tail_pval(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t right_offset, double tot_prob, double right_prob, double tail_sum, uint32_t new_m11);

void fisher22_precomp_pval_bounds(double pval, uint32_t row1_sum, uint32_t col1_sum, uint32_t total, uint32_t* m11_minp, uint32_t* m11_maxp, double* low_pvals, double* high_pvals, double* tot_probp, double* right_probp, double* tail_sump, uint32_t* m11p, uint32_t* right_offsetp) {

double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23);

void chi22_precomp_thresh(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t* m11_minp, uint32_t* m11_maxp, uint32_t* tiel, uint32_t* tieh);

void chi22_precomp_coeffs(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, double* coeff_ap, double* coeff_bp, double* coeff_cp);

#endif
