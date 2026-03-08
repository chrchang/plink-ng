#ifndef __PLINK2_STATS_H__
#define __PLINK2_STATS_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
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

#include "plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// result of current PToChisq(5e-324, 1)
// static const double kMaxChisq1df = 1480.8852530551483;

// now thread-safe!
double ChisqToP(double chisq, uint32_t df);

double ChisqToLnP(double chisq, uint32_t df);

// only handles df=1 and 2 for now, plan to support 4 later
// double PToChisq(double pval, uint32_t df);

// only handles df=1 for now
double LnPToChisq(double ln_pval);

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

// Assumes xx is a nonnegative integer.
double Lfact(double xx);

// HweP() has been replaced by HweLnP().  HweThresh() and HweThreshMidp() have
// been replaced by HweThreshLn().

// probable todo: provide HweLnPEx() where caller can provide
// near-tie-resolution workspace, and a workspace query function
BoolErr HweLnP(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp, double* resultp);

// these return 0 if close enough to Hardy-Weinberg equilibrium
BoolErr HweThresh(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double pval_thresh, uint32_t* out_of_eqp);

BoolErr HweThreshMidp(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double pval_thresh, uint32_t* out_of_eqp);

BoolErr HweThreshLnMain(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp, double ln_thresh, uint32_t* out_of_eqp);

HEADER_INLINE BoolErr HweThreshLn(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp, double thresh, double ln_thresh, uint32_t* out_of_eqp) {
  // kLnNormalMin = -708.3964185...
  if (ln_thresh > -708.396) {
    if (!midp) {
      return HweThresh(obs_hets, obs_hom1, obs_hom2, thresh, out_of_eqp);
    } else {
      return HweThreshMidp(obs_hets, obs_hom1, obs_hom2, thresh, out_of_eqp);
    }
  }
  return HweThreshLnMain(obs_hets, obs_hom1, obs_hom2, midp, ln_thresh, out_of_eqp);
}

BoolErr HweXchrLnP(int32_t female_hets, int32_t female_hom1, int32_t female_hom2, int32_t male1, int32_t male2, uint32_t midp, double* resultp);

#ifdef __cplusplus
}
#endif

#endif  // __PLINK2_STATS_H__
