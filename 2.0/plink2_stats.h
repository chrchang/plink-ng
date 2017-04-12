#ifndef __PLINK2_STATS_H__
#define __PLINK2_STATS_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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


#include "plink2_matrix.h"

#ifdef __cplusplus
namespace plink2 {
#endif

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

// outer_buf = constraint_ct
// inner_buf = constraint_ct x constraint_ct
// tmphxs_buf and h_transpose_buf are constraint_ct x predictor_ct
// mi_buf only needs to be of length 2 * constraint_ct
boolerr_t linear_hypothesis_chisq_f(const float* coef, const float* constraints_con_major, const float* cov_matrix, uint32_t constraint_ct, uint32_t predictor_ct, uint32_t cov_stride, double* chisq_ptr, float* tmphxs_buf, float* h_transpose_buf, float* inner_buf, matrix_finvert_buf1_t* mi_buf, float* outer_buf);

boolerr_t linear_hypothesis_chisq(const double* coef, const double* constraints_con_major, const double* cov_matrix, uintptr_t constraint_ct, uintptr_t predictor_ct, double* chisq_ptr, double* tmphxs_buf, double* h_transpose_buf, double* inner_buf, matrix_invert_buf1_t* mi_buf, double* outer_buf);

#ifdef __cplusplus
}
#endif

#endif // __PLINK2_STATS_H__
