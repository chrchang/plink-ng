#ifndef __WDIST_STATS_H__
#define __WDIST_STATS_H__

#include "wdist_common.h"

double chiprob_p(double xx, double df);

double inverse_chiprob(double qq, double df);

double ltqnorm(double p);

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22);

void fisher22_precomp_thresh(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t* m11_minp, uint32_t* m11_maxp, uint32_t* tiel, uint32_t* tieh);

double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23);

void chi22_precomp_thresh(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t* m11_minp, uint32_t* m11_maxp, uint32_t* tiel, uint32_t* tieh);

#endif
