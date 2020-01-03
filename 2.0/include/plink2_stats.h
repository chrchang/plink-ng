#ifndef __PLINK2_STATS_H__
#define __PLINK2_STATS_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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

#include "plink2_string.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// result of current PToChisq(5e-324, 1)
// static const double kMaxChisq1df = 1480.8852530551483;

// now thread-safe!
double ChisqToP(double chisq, uint32_t df);

double ChisqToLnP(double chisq, uint32_t df);

// only handles df=1 and 2 for now, plan to support 4 later
double PToChisq(double pval, uint32_t df);

// only handles df=1 for now
double LnPToChisq(double ln_pval);

// Better than the original TstatToP() when the same df comes up many times.
// (May want a higher-level interface that allocates and incrementally fills a
// table.)
double TstatToP2(double tt, double df, double cached_gamma_mult);

// No -9 error return since that's a legitimate p-value logarithm.  Caller is
// responsible for validating input.
double TstatToLnP(double tt, uint32_t df);

double FstatToLnP(double ff, uint32_t df1, uint32_t df2);

double QuantileToZscore(double pval);

HEADER_INLINE double ZscoreToP(double zz) {
  return ChisqToP(zz * zz, 1);
}

HEADER_INLINE double ZscoreToLnP(double zz) {
  return ChisqToLnP(zz * zz, 1);
}

double HweP(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp);

// returns 0 if close enough to Hardy-Weinberg equilibrium
uint32_t HweThresh(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);

uint32_t HweThreshMidp(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);

double FisherExact2x2P(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t midp);

double HweXchrP(int32_t female_hets, int32_t female_hom1, int32_t female_hom2, int32_t male1, int32_t male2, uint32_t midp);

#ifdef __cplusplus
}
#endif

#endif  // __PLINK2_STATS_H__
