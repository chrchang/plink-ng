#ifndef __PLINK2_STATS_H__
#define __PLINK2_STATS_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_cmdline.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// result of current inverse_chiprob(5e-324, 1)
const double kMaxInverseChiprob1df = 1480.8852530551483;

// now thread-safe!
double chiprob_p(double chisq, uint32_t df);

double inverse_chiprob(double pval, uint32_t df);

double calc_tprob(double tt, double df);

double calc_tprob2(double tt, double df, double cached_gamma_mult);

double ltqnorm(double p);

double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp);

// returns 0 if close enough to Hardy-Weinberg equilibrium
uint32_t SNPHWE_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);

uint32_t SNPHWE_midp_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t midp);

double SNPHWEX(int32_t female_hets, int32_t female_hom1, int32_t female_hom2, int32_t male1, int32_t male2, uint32_t midp);

#ifdef __cplusplus
}
#endif

#endif // __PLINK2_STATS_H__
