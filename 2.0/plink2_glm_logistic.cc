// This file is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
// Christopher Chang.
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

#include "plink2_glm_logistic.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "include/pgenlib_misc.h"
#include "include/plink2_bits.h"
#include "include/plink2_fmath.h"
#include "include/plink2_stats.h"
#include "include/plink2_string.h"
#include "include/plink2_thread.h"
#include "plink2_compress_stream.h"
#include "plink2_decompress.h"
#include "plink2_matrix.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// refugee from plink2_stats.h; important to remove its plink2_matrix.h
// dependency

// outer_buf = constraint_ct
// inner_buf = constraint_ct x constraint_ct
// tmphxs_buf and h_transpose_buf are constraint_ct x predictor_ct
// mi_buf only needs to be of length 2 * constraint_ct
BoolErr LinearHypothesisChisqF(const float* coef, const float* constraints_con_major, const float* cov_matrix, uint32_t constraint_ct, uint32_t predictor_ct, uint32_t cov_stride, double* chisq_ptr, float* tmphxs_buf, float* h_transpose_buf, float* inner_buf, double* half_inverted_buf, MatrixInvertBuf1* mi_buf, double* dbl_2d_buf, float* outer_buf) {
  ColMajorFvectorMatrixMultiplyStrided(coef, constraints_con_major, predictor_ct, predictor_ct, constraint_ct, outer_buf);
  // h-transpose does not have a special stride
  FmatrixTransposeCopy(constraints_con_major, constraint_ct, predictor_ct, predictor_ct, h_transpose_buf);
  ColMajorFmatrixMultiplyStrided(h_transpose_buf, cov_matrix, constraint_ct, constraint_ct, predictor_ct, cov_stride, predictor_ct, constraint_ct, tmphxs_buf);
  // tmp[][] is now predictor-major
  ColMajorFmatrixMultiplyStrided(tmphxs_buf, constraints_con_major, constraint_ct, constraint_ct, constraint_ct, predictor_ct, predictor_ct, constraint_ct, inner_buf);

  if (InvertFmatrixFirstHalf(constraint_ct, constraint_ct, inner_buf, half_inverted_buf, mi_buf, dbl_2d_buf)) {
    return 1;
  }
  InvertFmatrixSecondHalf(constraint_ct, constraint_ct, half_inverted_buf, inner_buf, mi_buf, dbl_2d_buf);
  double result = 0.0;
  const float* inner_iter = inner_buf;
  if (constraint_ct > kDotprodFThresh) {
    for (uint32_t constraint_idx = 0; constraint_idx != constraint_ct; ++constraint_idx) {
      result += S_CAST(double, DotprodF(inner_iter, outer_buf, constraint_ct) * outer_buf[constraint_idx]);
      inner_iter = &(inner_iter[constraint_ct]);
    }
  } else {
    for (uint32_t constraint_idx = 0; constraint_idx != constraint_ct; ++constraint_idx) {
      result += S_CAST(double, DotprodFShort(inner_iter, outer_buf, constraint_ct) * outer_buf[constraint_idx]);
      inner_iter = &(inner_iter[constraint_ct]);
    }
  }
  if (result < 0.0) {
    // guard against floating point error
    result = 0.0;
  }
  *chisq_ptr = result;
  return 0;
}

// Only called by GlmLogisticThreadF(), so there are the following differences
// from CheckMaxCorrAndVif():
// * predictor_dotprods is not precomputed; we start with predictors_pmaj
//   instead.
// * predictors_pmaj already has the intercept stripped off, so we don't need
//   relevant_predictor_ct := predictor_ct - 1, etc.
// * sample_stride parameter added, since predictors_pmaj has vector-aligned
//   rather than packed rows.
// * dbl_2d_buf not assumed to be filled with row sums, we compute them here.
//   (probably want to modify CheckMaxCorrAndVif() to do the same.)
// This now uses double-precision arithmetic since matrix inversion is too
// inconsistent if we stick to single-precision.
GlmErr CheckMaxCorrAndVifF(const float* predictors_pmaj, uint32_t predictor_ct, uint32_t sample_ct, uint32_t sample_stride, double max_corr, double vif_thresh, float* predictor_dotprod_buf, double* dbl_2d_buf, double* inverse_corr_buf, MatrixInvertBuf1* inv_1d_buf) {
  MultiplySelfTransposeStridedF(predictors_pmaj, predictor_ct, sample_ct, sample_stride, predictor_dotprod_buf);
  for (uintptr_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    const float* predictor_row = &(predictors_pmaj[pred_idx * sample_stride]);
    double row_sum = 0.0;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      row_sum += S_CAST(double, predictor_row[sample_idx]);
    }
    dbl_2d_buf[pred_idx] = row_sum;
  }
  const uint32_t predictor_ct_p1 = predictor_ct + 1;
  const double sample_ct_recip = 1.0 / u31tod(sample_ct);
  const double sample_ct_m1_d = u31tod(sample_ct - 1);
  const double sample_ct_m1_recip = 1.0 / sample_ct_m1_d;
  for (uint32_t pred_idx1 = 0; pred_idx1 != predictor_ct; ++pred_idx1) {
    double* sample_cov_row = &(inverse_corr_buf[pred_idx1 * predictor_ct]);
    const float* predictor_dotprod_row = &(predictor_dotprod_buf[pred_idx1 * predictor_ct]);
    const double pred1_mean_adj = dbl_2d_buf[pred_idx1] * sample_ct_recip;
    for (uint32_t pred_idx2 = 0; pred_idx2 <= pred_idx1; ++pred_idx2) {
      sample_cov_row[pred_idx2] = (S_CAST(double, predictor_dotprod_row[pred_idx2]) - pred1_mean_adj * dbl_2d_buf[pred_idx2]) * sample_ct_m1_recip;
    }
  }
  // now use dbl_2d_buf to store inverse-sqrts, to get to correlation matrix
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    dbl_2d_buf[pred_idx] = 1.0 / sqrt(inverse_corr_buf[pred_idx * predictor_ct_p1]);
  }
  // invert_symmdef_matrix only cares about bottom left of inverse_corr_buf[]
  for (uint32_t pred_idx1 = 1; pred_idx1 != predictor_ct; ++pred_idx1) {
    const double inverse_stdev1 = dbl_2d_buf[pred_idx1];
    // convert from covariances to correlations here
    double* corr_row_iter = &(inverse_corr_buf[pred_idx1 * predictor_ct]);
    const double* inverse_stdev2_iter = dbl_2d_buf;
    for (uintptr_t pred_idx2 = 0; pred_idx2 != pred_idx1; ++pred_idx2) {
      const double cur_corr = (*corr_row_iter) * inverse_stdev1 * (*inverse_stdev2_iter++);
      if (fabs(cur_corr) > max_corr) {
        return SetGlmErr2(kGlmErrcodeCorrTooHigh, pred_idx2, pred_idx1);
      }
      *corr_row_iter++ = cur_corr;
    }
  }
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    inverse_corr_buf[pred_idx * predictor_ct_p1] = 1.0;
  }
  if (InvertSymmdefMatrixChecked(predictor_ct, inverse_corr_buf, inv_1d_buf, dbl_2d_buf)) {
    return SetGlmErr0(kGlmErrcodeVifInfinite);
  }

  // VIFs = diagonal elements of inverse correlation matrix
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    if (inverse_corr_buf[pred_idx * predictor_ct_p1] > vif_thresh) {
      return SetGlmErr1(kGlmErrcodeVifTooHigh, pred_idx);
    }
  }
  return 0;
}

static const float kSmallFloatPairs[32] = PAIR_TABLE16(0.0, 1.0, 2.0, 3.0);

static const float kSmallInvFloatPairs[32] = PAIR_TABLE16(2.0, 1.0, 0.0, 3.0);

static const float kSmallInvFloats[4] = {2.0, 1.0, 0.0, 3.0};

uint32_t GenoarrToFloatsRemoveMissing(const uintptr_t* genoarr, const float* __restrict table, uint32_t sample_ct, float* __restrict dst) {
  assert(sample_ct);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t subgroup_len = kBitsPerWordD2;
  float* dst_iter = dst;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        return S_CAST(uint32_t, dst_iter - dst);
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      const uintptr_t cur_geno = geno_word & 3;
      if (cur_geno < 3) {
        *dst_iter++ = table[cur_geno];
      }
      geno_word >>= 2;
    }
  }
}

// #####
// The following code is based on the winning submission of Pascal Pons in the
// "GWASSpeedup" contest run in April 2013 by Babbage Analytics & Innovation
// and TopCoder, who have donated the results to be used in PLINK.  See:
//   Hill A, Loh PR, Bharadwaj RB, Pons P, Shang J, Guinan E, Lakhani K,
//   Kilty I, Jelinsky SA (2017) Stepwise Distributed Open Innovation Contests
//   for Software Development - Acceleration of Genome-Wide Association
//   Analysis.  Gigascience, 6.
// #####

#ifdef __LP64__
// For equivalent "normal" C/C++ code, see the non-__LP64__ versions of these
// functions.

// N.B. This requires all mm[] rows to be zero-padded at the end, and there
// can't be nan values at the end of vect[].  (The other way around works too.)
//
// This is currently a bit faster than sgemm and sgemv on my Mac, so it isn't
// appropriate to throw out this code yet.
#  ifdef FVEC_32
static inline void MultMatrixDxnVectNF(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* __restrict dest) {
  const uintptr_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  uint32_t row_idx = 0;
  __m256 s1;
  __m256 s2;
  __m256 s3;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    // Handle 4 rows at a time in this loop, regardless of vector size.
    for (; row_idx < row_ctm3; row_idx += 4) {
      s1 = _mm256_setzero_ps();
      s2 = _mm256_setzero_ps();
      s3 = _mm256_setzero_ps();
      __m256 s4 = _mm256_setzero_ps();
      for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
        const float* mm_ptr = &(mm[row_idx * col_ctav + col_idx]);
        const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
        __m256 a1 = _mm256_load_ps(mm_ptr);
        __m256 a2 = _mm256_load_ps(&(mm_ptr[col_ctav]));
        __m256 a3 = _mm256_load_ps(&(mm_ptr[2 * col_ctav]));
        __m256 a4 = _mm256_load_ps(&(mm_ptr[3 * col_ctav]));
        s1 = _mm256_fmadd_ps(a1, vv, s1);
        s2 = _mm256_fmadd_ps(a2, vv, s2);
        s3 = _mm256_fmadd_ps(a3, vv, s3);
        s4 = _mm256_fmadd_ps(a4, vv, s4);
      }
      *dest++ = VecFHsum(R_CAST(VecF, s1));
      *dest++ = VecFHsum(R_CAST(VecF, s2));
      *dest++ = VecFHsum(R_CAST(VecF, s3));
      *dest++ = VecFHsum(R_CAST(VecF, s4));
    }
  }
  s1 = _mm256_setzero_ps();
  s2 = _mm256_setzero_ps();
  s3 = _mm256_setzero_ps();
  switch (row_ct % 4) {
  case 3:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
      const float* mm_ptr = &(mm[row_idx * col_ctav + col_idx]);
      const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
      __m256 a1 = _mm256_load_ps(mm_ptr);
      __m256 a2 = _mm256_load_ps(&(mm_ptr[col_ctav]));
      __m256 a3 = _mm256_load_ps(&(mm_ptr[2 * col_ctav]));
      s1 = _mm256_fmadd_ps(a1, vv, s1);
      s2 = _mm256_fmadd_ps(a2, vv, s2);
      s3 = _mm256_fmadd_ps(a3, vv, s3);
    }
    *dest++ = VecFHsum(R_CAST(VecF, s1));
    *dest++ = VecFHsum(R_CAST(VecF, s2));
    *dest = VecFHsum(R_CAST(VecF, s3));
    break;
  case 2:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
      const float* mm_ptr = &(mm[row_idx * col_ctav + col_idx]);
      const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
      __m256 a1 = _mm256_load_ps(mm_ptr);
      __m256 a2 = _mm256_load_ps(&(mm_ptr[col_ctav]));
      s1 = _mm256_fmadd_ps(a1, vv, s1);
      s2 = _mm256_fmadd_ps(a2, vv, s2);
    }
    *dest++ = VecFHsum(R_CAST(VecF, s1));
    *dest = VecFHsum(R_CAST(VecF, s2));
    break;
  case 1:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
      const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
      __m256 a1 = _mm256_load_ps(&(mm[row_idx * col_ctav + col_idx]));
      s1 = _mm256_fmadd_ps(a1, vv, s1);
    }
    *dest = VecFHsum(R_CAST(VecF, s1));
    break;
  }
}

#  else  // !FVEC_32
static inline void MultMatrixDxnVectNF(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* __restrict dest) {
  const uint32_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  ColMajorFvectorMatrixMultiplyStrided(vect, mm, col_ct, col_ctav, row_ct, dest);
}

#  endif  // !FVEC_32
// !__LP64__

static inline void LogisticSseF(uint32_t nn, float* vect) {
  const VecF zero = vecf_setzero();
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF aa = *R_CAST(VecF*, &(vect[uii]));
    aa = zero - aa;
    // tried substituting in vexpf() here on OS X; it was slower without being
    // significantly more accurate.
    aa = fmath_exp_ps(aa);
    aa = aa + one;
    aa = one / aa;
    *R_CAST(VecF*, &(vect[uii])) = aa;
  }
}

static inline void ComputeVAndPMinusYF(const float* yy, uint32_t nn, float* __restrict pp, float* __restrict vv) {
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF ptmp = *R_CAST(VecF*, &(pp[uii]));
    VecF one_minus_ptmp = one - ptmp;
    *R_CAST(VecF*, &(vv[uii])) = ptmp * one_minus_ptmp;
    VecF ytmp = *R_CAST(const VecF*, &(yy[uii]));
    *R_CAST(VecF*, &(pp[uii])) = ptmp - ytmp;
  }
}

static inline void ComputeVF(const float* pp, uint32_t nn, float* __restrict vv) {
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF ptmp = *R_CAST(const VecF*, &(pp[uii]));
    VecF one_minus_ptmp = one - ptmp;
    *R_CAST(VecF*, &(vv[uii])) = ptmp * one_minus_ptmp;
  }
}

static inline float TripleProductF(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  VecF sum = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF aa = *R_CAST(const VecF*, &(v1[uii]));
    VecF bb = *R_CAST(const VecF*, &(v2[uii]));
    VecF cc = *R_CAST(const VecF*, &(v3[uii]));
    sum = sum + aa * bb * cc;
  }
  return VecFHsum(sum);
}

static inline void ComputeTwoDiagTripleProductF(const float* aa, const float* bb, const float* vv, uint32_t nn, float* __restrict raa_ptr, float* __restrict rab_ptr, float* __restrict rbb_ptr) {
  VecF saa = vecf_setzero();
  VecF sab = vecf_setzero();
  VecF sbb = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    const VecF vtmp = *R_CAST(const VecF*, &(vv[uii]));
    const VecF atmp = *R_CAST(const VecF*, &(aa[uii]));
    const VecF btmp = *R_CAST(const VecF*, &(bb[uii]));
    const VecF av = atmp * vtmp;
    const VecF bv = btmp * vtmp;
    saa = saa + atmp * av;
    sab = sab + atmp * bv;
    sbb = sbb + btmp * bv;
  }
  *raa_ptr = VecFHsum(saa);
  *rab_ptr = VecFHsum(sab);
  *rbb_ptr = VecFHsum(sbb);
}

static inline void ComputeThreeTripleProductF(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  VecF s1 = vecf_setzero();
  VecF s2 = vecf_setzero();
  VecF s3 = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    const VecF a1tmp = *R_CAST(const VecF*, &(a1[uii]));
    const VecF a2tmp = *R_CAST(const VecF*, &(a2[uii]));
    const VecF a3tmp = *R_CAST(const VecF*, &(a3[uii]));
    const VecF vtmp = *R_CAST(const VecF*, &(vv[uii]));
    VecF btmp = *R_CAST(const VecF*, &(bb[uii]));
    btmp = btmp * vtmp;
    s1 = s1 + a1tmp * btmp;
    s2 = s2 + a2tmp * btmp;
    s3 = s3 + a3tmp * btmp;
  }
  *r1_ptr = VecFHsum(s1);
  *r2_ptr = VecFHsum(s2);
  *r3_ptr = VecFHsum(s3);
}

static inline void ComputeTwoPlusOneTripleProductF(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  VecF s1 = vecf_setzero();
  VecF s2 = vecf_setzero();
  VecF s3 = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    const VecF a1tmp = *R_CAST(const VecF*, &(a1[uii]));
    const VecF a2tmp = *R_CAST(const VecF*, &(a2[uii]));
    const VecF btmp = *R_CAST(const VecF*, &(bb[uii]));
    const VecF vtmp = *R_CAST(const VecF*, &(vv[uii]));
    const VecF bv = btmp * vtmp;
    s1 = s1 + btmp * bv;
    s2 = s2 + a1tmp * bv;
    s3 = s3 + a2tmp * bv;
  }
  *r1_ptr = VecFHsum(s1);
  *r2_ptr = VecFHsum(s2);
  *r3_ptr = VecFHsum(s3);
}
#else  // no __LP64__ (and hence, unsafe to assume presence of SSE2)
static inline void LogisticSseF(uint32_t nn, float* vect) {
  // We use explicit static_cast<float> instead of e.g. 1.0f because
  // handling of the latter is actually implementation-specific; see
  //   http://nullprogram.com/blog/2018/05/01/
  // In particular, that blog post claims that
  //   int float_compare() {
  //     float x = 1.3f;
  //     return x == 1.3f;
  //   }
  // returns 0 under gcc and 1 under clang (with -std=c99 -m32, which is one of
  // plink2's compilation settings)????!!!!!!!
  // Unless the author is outright mistaken, this suggests that use of the f
  // suffix should be considered a bug ~100% of the time.
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vect[uii] = S_CAST(float, 1.0) / (1 + expf(-vect[uii]));
  }
}

static inline void ComputeVAndPMinusYF(const float* yy, uint32_t nn, float* __restrict pp, float* __restrict vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vv[uii] = pp[uii] * (S_CAST(float, 1.0) - pp[uii]);
    pp[uii] -= yy[uii];
  }
}

static inline void ComputeVF(const float* pp, uint32_t nn, float* __restrict vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vv[uii] = pp[uii] * (S_CAST(float, 1.0) - pp[uii]);
  }
}

static inline void MultMatrixDxnVectNF(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* __restrict dest) {
  const uint32_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  ColMajorFvectorMatrixMultiplyStrided(vect, mm, col_ct, col_ctav, row_ct, dest);
}

static inline float TripleProductF(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  float fxx = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    fxx += v1[uii] * v2[uii] * v3[uii];
  }
  return fxx;
}

static inline void ComputeTwoDiagTripleProductF(const float* aa, const float* bb, const float* vv, uint32_t nn, float* __restrict raa_ptr, float* __restrict rab_ptr, float* __restrict rbb_ptr) {
  float raa = 0.0;
  float rab = 0.0;
  float rbb = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = aa[uii];
    const float fyy = bb[uii];
    float fzz = vv[uii];
    raa += fxx * fxx * fzz;
    fzz *= fyy;
    rab += fxx * fzz;
    rbb += fyy * fzz;
  }
  *raa_ptr = raa;
  *rab_ptr = rab;
  *rbb_ptr = rbb;
}

static inline void ComputeThreeTripleProductF(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = bb[uii] * vv[uii];
    r1 += a1[uii] * fxx;
    r2 += a2[uii] * fxx;
    r3 += a3[uii] * fxx;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}

static inline void ComputeTwoPlusOneTripleProductF(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = bb[uii];
    const float fyy = fxx * vv[uii];
    r1 += fxx * fyy;
    r2 += a1[uii] * fyy;
    r3 += a2[uii] * fyy;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}
#endif
BoolErr ComputeLoglikCheckedF(const float* yy, const float* pp, uint32_t sample_ct, double* loglik_ptr) {
  // possible todo: look for a high-precision way to accelerate this.
  double loglik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const double new_pi = S_CAST(double, pp[sample_idx]);
    if ((new_pi == 0.0) || (new_pi == 1.0)) {
      return 1;
    }
    loglik += (yy[sample_idx] != S_CAST(float, 0.0))? log(new_pi) : log1p(-new_pi);
  }
  *loglik_ptr = loglik;
  return 0;
}

// M V M^T
// This is the biggest logistic/Firth regression bottleneck.
// Tried to replace this with sqrt(v) followed by ssyrk, but that was slower.
// Also tried to take advantage of the first row of M being constant-1, and
// managed to fail.
void ComputeHessianF(const float* mm, const float* vv, uint32_t col_ct, uint32_t row_ct, float* __restrict dest) {
  const uintptr_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  const uintptr_t row_ctav = RoundUpPow2(row_ct, kFloatPerFVec);
  const uintptr_t row_ctavp1 = row_ctav + 1;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    for (uint32_t row_idx = 0; row_idx < row_ctm3; row_idx += 3) {
      const float* mm_cur = &(mm[row_idx * col_ctav]);
      ComputeTwoDiagTripleProductF(mm_cur, &(mm_cur[col_ctav]), vv, col_ct, &(dest[row_idx * row_ctavp1]), &(dest[(row_idx + 1) * row_ctavp1 - 1]), &(dest[(row_idx + 1) * row_ctavp1]));
      ComputeTwoPlusOneTripleProductF(&(mm_cur[2 * col_ctav]), &(mm_cur[col_ctav]), mm_cur, vv, col_ct, &(dest[(row_idx + 2) * row_ctavp1]), &(dest[(row_idx + 2) * row_ctavp1 - 1]), &(dest[(row_idx + 2) * row_ctavp1 - 2]));
      for (uint32_t row_idx2 = row_idx + 3; row_idx2 != row_ct; ++row_idx2) {
        ComputeThreeTripleProductF(&(mm[row_idx2 * col_ctav]), mm_cur, &(mm_cur[col_ctav]), &(mm_cur[2 * col_ctav]), vv, col_ct, &(dest[row_idx2 * row_ctav + row_idx]), &(dest[row_idx2 * row_ctav + row_idx + 1]), &(dest[row_idx2 * row_ctav + row_idx + 2]));
      }
    }
  }
  switch (row_ct % 3) {
  case 0:
    ComputeTwoPlusOneTripleProductF(&(mm[(row_ct - 3) * col_ctav]), &(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 3) * row_ctavp1]), &(dest[(row_ct - 2) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1 - 2]));
    // fall through
  case 2:
    ComputeTwoDiagTripleProductF(&(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 2) * row_ctavp1]), &(dest[(row_ct - 1) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1]));
    break;
  case 1:
    dest[(row_ct - 1) * row_ctavp1] = TripleProductF(&(mm[(row_ct - 1) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct);
  }
}

void CholeskyDecompositionF(const float* aa, uint32_t predictor_ct, float* __restrict ll) {
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  const uintptr_t predictor_ctavp1 = predictor_ctav + 1;
  const float* aa_diag_elem_ptr = aa;
  float* cur_ll_row = ll;
  for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
    float fxx = *aa_diag_elem_ptr;
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      const float fyy = cur_ll_row[col_idx];
      fxx -= fyy * fyy;
    }
    float fyy;
    if (fxx >= S_CAST(float, 0.0)) {
      fyy = sqrtf(fxx);
    } else {
      fyy = S_CAST(float, 1e-6);
    }
    cur_ll_row[row_idx] = fyy;
    fyy = S_CAST(float, 1.0) / fyy;  // now 1.0 / L[j][j]
    const float* aa_col_iter = aa_diag_elem_ptr;
    float* cur_ll_row2 = cur_ll_row;
    for (uint32_t row_idx2 = row_idx + 1; row_idx2 != predictor_ct; ++row_idx2) {
      aa_col_iter = &(aa_col_iter[predictor_ctav]);
      float fxx2 = *aa_col_iter;
      cur_ll_row2 = &(cur_ll_row2[predictor_ctav]);
      for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
        fxx2 -= cur_ll_row[col_idx] * cur_ll_row2[col_idx];
      }
      cur_ll_row2[row_idx] = fxx2 * fyy;
    }
    aa_diag_elem_ptr = &(aa_diag_elem_ptr[predictor_ctavp1]);
    cur_ll_row = &(cur_ll_row[predictor_ctav]);
  }
}

void SolveLinearSystemF(const float* ll, const float* yy, uint32_t predictor_ct, float* __restrict xx) {
  // Finds x such that y = L(L^T)x, via forward and backward substitution
  //
  // might want to use this in NOLAPACK case only, since we can now produce
  // 32-bit Linux builds with statically linked LAPACK
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  {
    const float* cur_ll_row = ll;
    for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
      float fxx = yy[row_idx];
      for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
        // Note that it doesn't matter what xx starts with.  First iteration of
        // outer loop reads no values of xx and sets xx[0], second iteration
        // reads xx[0] and sets xx[1], etc.
        fxx -= cur_ll_row[col_idx] * xx[col_idx];
      }
      xx[row_idx] = fxx / cur_ll_row[row_idx];
      cur_ll_row = &(cur_ll_row[predictor_ctav]);
    }
  }
  for (uint32_t col_idx = predictor_ct; col_idx; ) {
    float* xx_stop = &(xx[--col_idx]);
    float fxx = *xx_stop;
    const float* ll_col_iter = &(ll[(predictor_ct - 1) * predictor_ctav + col_idx]);
    for (float* xx_iter = &(xx[predictor_ct - 1]); xx_iter != xx_stop; --xx_iter) {
      fxx -= (*ll_col_iter) * (*xx_iter);
      ll_col_iter -= predictor_ctav;
    }
    *xx_stop = fxx / (*ll_col_iter);
  }
}

BoolErr LogisticRegressionF(const float* yy, const float* xx, const float* sample_offsets, uint32_t sample_ct, uint32_t predictor_ct, float* __restrict coef, uint32_t* is_unfinished_ptr, float* __restrict ll, float* __restrict pp, float* __restrict vv, float* __restrict hh, float* __restrict grad, float* __restrict dcoef) {
  // Similar to first part of logistic.cpp fitLM(), but incorporates changes
  // from Pascal Pons et al.'s TopCoder code.
  //
  // Preallocated buffers (initial contents irrelevant):
  // vv    = sample variance buffer
  // hh    = hessian matrix buffer, predictor_ct^2, rows vector-aligned
  // grad  = gradient buffer Y[] (length predictor_ct)
  // dcoef = current coefficient change buffer (length predictor_ct)
  //
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         vector-aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype; trailing elements must be zeroed out
  //
  // Input/output:
  // coef  = starting point, overwritten with logistic regression betas.
  //
  // Outputs:
  // ll    = cholesky decomposition matrix, predictor_ct^2, rows vector-aligned
  // pp    = final likelihoods minus Y[] (not currently used by callers)
  //
  // Returns 1 on convergence failure, 0 otherwise.
  // is_unfinished assumed to be initialized to 0, and is set to 1 if we hit
  // the iteration limit without satisfying other convergence criteria; the
  // main return value is 0 in this case.
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  float min_delta_coef = 1e9;

  ZeroFArr(predictor_ct * predictor_ctav, ll);
  // quasi-bugfix (9 Oct 2024)
  // technically wasn't needed because LogisticSseF() implementation mapped nan
  // to not-nan, but we don't want to depend on that
  ZeroFArr(sample_ctav - sample_ct, &(pp[sample_ct]));
  ZeroFArr(sample_ctav - sample_ct, &(vv[sample_ct]));
  for (uint32_t iteration = 0; ; ++iteration) {
    // P[i] = \sum_j X[i][j] * coef[j];
    ColMajorFmatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);
    if (sample_offsets) {
      AddFVec(sample_offsets, sample_ctav, pp);
    }
    // Possible future optimization: suppose categorical covariates are
    // represented as categorical-covariate-major uint16_t* kk, indicating for
    // each sample which raw covariate index is 1 (with one "fallow" covariate
    // index at the end).
    // Then the above expression becomes
    // P[i] = \sum_j^{regular} X[i][j] * coef[j] +
    //        \sum_j^{cats} coef[K[i][j]]

    // P[i] = 1 / (1 + exp(-P[i]));
    LogisticSseF(sample_ct, pp);

    // V[i] = P[i] * (1 - P[i]);
    // P[i] -= Y[i];
    ComputeVAndPMinusYF(yy, sample_ct, pp, vv);

    // Possible categorical optimizations:
    // 1. skip terms between different categories of the same covariate
    // 2. all same-category terms within the same covariate can be handled with
    //    a single loop over the samples
    // 3. terms involving a regular covariate and a categorical covariate can
    //    be handled with one loop over the samples; covers all categories
    // 4. similarly, one loop over the samples is enough to update all category
    //    pairs between two categorical covariates
    ComputeHessianF(xx, vv, sample_ct, predictor_ct, hh);

    // grad = X^T P
    // Separate categorical loop also possible here
    MultMatrixDxnVectNF(xx, pp, sample_ct, predictor_ct, grad);

    CholeskyDecompositionF(hh, predictor_ct, ll);

    SolveLinearSystemF(ll, grad, predictor_ct, dcoef);

    float delta_coef = 0.0;
    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      const float cur_dcoef = dcoef[pred_idx];
      delta_coef += fabsf(cur_dcoef);
      coef[pred_idx] -= cur_dcoef;
    }
    if (delta_coef < min_delta_coef) {
      min_delta_coef = delta_coef;
    }
    if (delta_coef != delta_coef) {
      return 1;
    }
    if (iteration > 3) {
      if (((delta_coef > S_CAST(float, 20.0)) && (delta_coef > 2 * min_delta_coef)) || ((iteration > 6) && (fabsf(S_CAST(float, 1.0) - delta_coef) < S_CAST(float, 1e-3)))) {
        return 1;
      }
      if (iteration > 13) {
        // If fabsf(any coefficient) > 8e3, this is almost certainly a form of
        // convergence failure that didn't get caught by the
        // (delta_coef > 20.0) check due to a precision quirk.  (8e3 threshold
        // ~= 1e-4 * 2^23, since floats have 23 bits of precision)
        for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
          if (fabsf(coef[pred_idx]) > S_CAST(float, 8e3)) {
            return 1;
          }
        }
        *is_unfinished_ptr = 1;
        return 0;
      }
    }
    // Pons reported that 1.1e-3 was dangerous, so I agree with the decision to
    // tighten this threshold from 1e-3 to 1e-4.
    if (delta_coef < S_CAST(float, 1e-4)) {
      // Be more conservative in throwing out results when we don't hit the
      // iteration limit.
      for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
        if (fabsf(coef[pred_idx]) > S_CAST(float, 6e4)) {
          return 1;
        }
      }
      return 0;
    }
  }
}

#ifdef __LP64__
void CopyAndMeanCenterF(const float* src, uintptr_t ct, float* __restrict dst) {
  const uintptr_t fullvec_ct = ct / kFloatPerFVec;
  const VecF* src_alias = R_CAST(const VecF*, src);
  VecF vsum = vecf_setzero();
  for (uintptr_t vidx = 0; vidx != fullvec_ct; ++vidx) {
    vsum += src_alias[vidx];
  }
  float sum = VecFHsum(vsum);
  const uintptr_t trailing_start_idx = fullvec_ct * kFloatPerFVec;
  for (uintptr_t ulii = trailing_start_idx; ulii != ct; ++ulii) {
    sum += src[ulii];
  }

  const float neg_mean = -sum / S_CAST(float, ct);
  const VecF neg_vmean = VCONST_F(neg_mean);
  VecF* dst_alias = R_CAST(VecF*, dst);
  for (uintptr_t vidx = 0; vidx != fullvec_ct; ++vidx) {
    dst_alias[vidx] = src_alias[vidx] + neg_vmean;
  }
  if (trailing_start_idx != ct) {
    for (uintptr_t ulii = trailing_start_idx; ulii != ct; ++ulii) {
      dst[ulii] = src[ulii] + neg_mean;
    }
    const uintptr_t trailing_stop_idx = trailing_start_idx + kFloatPerFVec;
    for (uintptr_t ulii = ct; ulii != trailing_stop_idx; ++ulii) {
      dst[ulii] = S_CAST(float, 0.0);
    }
  }
}
#else
void CopyAndMeanCenterF(const float* src, uintptr_t ct, float* __restrict dst) {
  float sum = 0.0;
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    sum += src[ulii];
  }
  const float mean = sum / u31tof(ct);
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    dst[ulii] = src[ulii] - mean;
  }
}
#endif

BoolErr LogisticRegressionResidualizedF(const float* yy, const float* xx, const uintptr_t* sample_nm, const CcResidualizeCtx* cc_residualize, uint32_t nm_sample_ct, uint32_t orig_predictor_ct, float* coef, uint32_t* is_unfinished_ptr, float* ll, float* pp, float* vv, float* hh, float* grad, float* dcoef, float* mean_centered_pmaj_buf, float* sample_offsets_buf) {
  if (!cc_residualize->logistic_nm_sample_offsets_f) {
    return 1;
  }
  const uintptr_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kFloatPerFVec);
  const uint32_t domdev_present_p1 = cc_residualize->domdev_present_p1;
  for (uint32_t geno_idx = 0; geno_idx != domdev_present_p1; ++geno_idx) {
    CopyAndMeanCenterF(&(xx[(geno_idx + 1) * nm_sample_ctav]), nm_sample_ct, &(mean_centered_pmaj_buf[geno_idx * nm_sample_ctav]));
  }
  const uint32_t prefitted_pred_ct = cc_residualize->prefitted_pred_ct;
  const uint32_t orig_biallelic_predictor_ct = domdev_present_p1 + prefitted_pred_ct;
  const uint32_t extra_allele_ct = orig_predictor_ct - orig_biallelic_predictor_ct;
  for (uint32_t extra_allele_idx = 0; extra_allele_idx != extra_allele_ct; ++extra_allele_idx) {
    CopyAndMeanCenterF(&(xx[(extra_allele_idx + orig_biallelic_predictor_ct) * nm_sample_ctav]), nm_sample_ct, &(mean_centered_pmaj_buf[(extra_allele_idx + domdev_present_p1) * nm_sample_ctav]));
  }
  const float* sample_offsets = cc_residualize->logistic_nm_sample_offsets_f;
  if (nm_sample_ct != cc_residualize->sample_ct) {
    uintptr_t sample_idx_base = 0;
    uintptr_t sample_nm_bits = sample_nm[0];
    for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(sample_nm, &sample_idx_base, &sample_nm_bits);
      sample_offsets_buf[uii] = sample_offsets[sample_idx];
    }
    // todo: check if this is actually needed
    const uint32_t remainder = (-nm_sample_ct) & (kFloatPerFVec - 1);
    ZeroFArr(remainder, &(sample_offsets_buf[nm_sample_ct]));
    sample_offsets = sample_offsets_buf;
  }
  // genotype, domdev?, other alleles
  const uint32_t regressed_predictor_ct = domdev_present_p1 + extra_allele_ct;
  const uint32_t regressed_predictor_ctav = RoundUpPow2(regressed_predictor_ct, kFloatPerFVec);
  if (LogisticRegressionF(yy, mean_centered_pmaj_buf, sample_offsets, nm_sample_ct, regressed_predictor_ct, &(coef[1]), is_unfinished_ptr, ll, pp, vv, hh, grad, dcoef)) {
    return 1;
  }
  // hh and ll are shifted up and to the left from what the caller expects, due
  // to the missing intercept.  Correct that here.
  // bugfix (4 Sep 2021): Initially thought only bottom-left triangle mattered,
  // but that's not true for genotypic/hethom case.  Also, wider stride
  // expected when regressed_predictor_ct is an exact multiple of
  // kFloatPerFVec.
  const uint32_t expected_predictor_ctav = RoundUpPow2(regressed_predictor_ct + 1, kFloatPerFVec);
  for (uint32_t write_row_idx = regressed_predictor_ct; write_row_idx; --write_row_idx) {
    memcpy(&(hh[write_row_idx * expected_predictor_ctav + 1]), &(hh[(write_row_idx - 1) * regressed_predictor_ctav]), regressed_predictor_ct * sizeof(float));
  }
  for (uint32_t write_row_idx = regressed_predictor_ct; write_row_idx; --write_row_idx) {
    memcpy(&(ll[write_row_idx * expected_predictor_ctav + 1]), &(ll[(write_row_idx - 1) * regressed_predictor_ctav]), regressed_predictor_ct * sizeof(float));
  }
  return 0;
}

// tmpNxK, interpreted as column-major, is sample_ct x predictor_ct
// X, interpreted as column-major, is also sample_ct x predictor_ct
// Hdiag[i] = V[i] (\sum_j tmpNxK[i][j] X[i][j])
void FirthComputeHdiagWeightsF(const float* yy, const float* xx, const float* pp, const float* hh, const float* vv, uint32_t predictor_ct, uint32_t predictor_ctav, uint32_t sample_ct, uint32_t sample_ctav, float* hdiag, float* ww, float* tmpnxk_buf) {
  ColMajorFmatrixMultiplyStrided(xx, hh, sample_ct, sample_ctav, predictor_ct, predictor_ctav, predictor_ct, sample_ctav, tmpnxk_buf);
#ifdef __LP64__
  const VecF half = VCONST_F(0.5);
  for (uint32_t sample_offset = 0; sample_offset < sample_ctav; sample_offset += kFloatPerFVec) {
    VecF dotprods = vecf_setzero();
    const float* xx_row = &(xx[sample_offset]);
    const float* tmpnxk_row = &(tmpnxk_buf[sample_offset]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const VecF cur_xx = *R_CAST(const VecF*, &(xx_row[pred_uidx * sample_ctav]));
      const VecF cur_tmpnxk = *R_CAST(const VecF*, &(tmpnxk_row[pred_uidx * sample_ctav]));
      dotprods = dotprods + cur_xx * cur_tmpnxk;
    }
    // Can handle categorical covariates in a separate loop here, and load the
    // dotprods increment into a union, etc.
    const VecF cur_vv = *R_CAST(const VecF*, &(vv[sample_offset]));
    const VecF cur_pi = *R_CAST(const VecF*, &(pp[sample_offset]));
    const VecF cur_yy = *R_CAST(const VecF*, &(yy[sample_offset]));
    const VecF cur_hdiag = cur_vv * dotprods;
    *R_CAST(VecF*, &(hdiag[sample_offset])) = cur_hdiag;
    const VecF half_minus_cur_pis = half - cur_pi;
    const VecF yy_minus_cur_pis = cur_yy - cur_pi;
    *R_CAST(VecF*, &(ww[sample_offset])) = yy_minus_cur_pis + cur_hdiag * half_minus_cur_pis;
  }
#else
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    float dotprod = 0.0;
    const float* xx_row = &(xx[sample_idx]);
    const float* tmpnxk_row = &(tmpnxk_buf[sample_idx]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      dotprod += xx_row[pred_uidx * sample_ctav] * tmpnxk_row[pred_uidx * sample_ctav];
    }
    const float cur_hdiag = vv[sample_idx] * dotprod;
    hdiag[sample_idx] = cur_hdiag;
    const float cur_pi = pp[sample_idx];
    ww[sample_idx] = (yy[sample_idx] - cur_pi) + cur_hdiag * (S_CAST(float, 0.5) - cur_pi);
  }
#endif
}

void FirthComputeSecondWeightsF(const float* hdiag, const float* vv, __maybe_unused uint32_t sample_ct, __maybe_unused uint32_t sample_ctav, float* ww) {
#ifdef __LP64__
  const VecF one = VCONST_F(1.0);
  for (uint32_t sample_offset = 0; sample_offset < sample_ctav; sample_offset += kFloatPerFVec) {
    const VecF cur_hdiag = *R_CAST(const VecF*, &(hdiag[sample_offset]));
    const VecF cur_vv = *R_CAST(const VecF*, &(vv[sample_offset]));
    *R_CAST(VecF*, &(ww[sample_offset])) = (one + cur_hdiag) * cur_vv;
  }
#else
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    ww[sample_idx] = (S_CAST(float, 1.0) + hdiag[sample_idx]) * vv[sample_idx];
  }
#endif
}

BoolErr FirthRegressionF(const float* yy, const float* xx, const float* sample_offsets, uint32_t sample_ct, uint32_t predictor_ct, float* beta, uint32_t* is_unfinished_ptr, float* hh, double* half_inverted_buf, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, float* pp, float* vv, float* ustar, float* delta, float* hdiag, float* ww, float* hh0, float* tmpnxk_buf) {
  // This is a port of Georg Heinze's logistf 1.24.1 R function (called with
  // pl=FALSE), adapted to use many of plink 1.9's optimizations; see
  //   https://github.com/georgheinze/logistf
  // (Before 30 Dec 2022, this function was based on an earlier version of
  // logistf.)
  //
  // Preallocated buffers (initial contents irrelevant):
  // half_inverted_buf, inv_1d_buf, dbl_2d_buf = for matrix inversion
  // pp    = likelihoods minus Y[] (not currently used by callers)
  // vv    = sample variance buffer
  // ustar = gradient buffer (length predictor_ct)
  // delta = current coefficient change buffer (length predictor_ct)
  // hdiag = intermediate multiplier for weights (length sample_ct)
  // ww    = weight buffer (length sample_ct)
  // hh0   = just like hh below, but not returned
  //
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         vector-aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype
  // sample_offsets = in residualized mode, product of pre-fitted covariate
  //                  beta vector with covariate matrix.  otherwise, nullptr.
  //
  // Input/output:
  // beta  = starting point, overwritten with regression betas.
  //
  // Outputs:
  // hh    = variance-covariance matrix buffer, predictor_ct^2, rows
  //         vector-aligned.  (spends some time as pre-inversion Hessian matrix
  //         too)
  //
  // Returns 1 on convergence failure, 0 otherwise.
  // is_unfinished assumed to be initialized to 0, and is set to 1 if we hit
  // the iteration limit without satisfying other convergence criteria; the
  // main return value is 0 in this case.
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);

  // bugfix (4 Nov 2017): ustar[] trailing elements must be zeroed out
  ZeroFArr(predictor_ctav - predictor_ct, &(ustar[predictor_ct]));

  const uint32_t trailing_sample_ct = sample_ctav - sample_ct;
  if (trailing_sample_ct) {
    // bunch of other trailing elements must also be zeroed out
    ZeroFArr(trailing_sample_ct, &(pp[sample_ct]));
    ZeroFArr(trailing_sample_ct, &(vv[sample_ct]));
    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      ZeroFArr(trailing_sample_ct, &(tmpnxk_buf[sample_ct + pred_idx * sample_ctav]));
    }
  }

  // start with 80% of most logistf convergence defaults (some reduction is
  // appropriate to be consistent with single-precision arithmetic).  (Update,
  // 9 Apr 2020: max_iter now matches logistf default since we've been shown a
  // concrete example where it matters in a way that isn't masked by the
  // single- vs. double-precision difference.)
  const uint32_t max_iter = 25;
  const float gconv = S_CAST(float, 0.0001);
  const float xconv = S_CAST(float, 0.0001);
  const double lconv = 0.0001;
  float delta_max = 0.0;
  double loglik_old = 0.0;
  for (uint32_t iter_idx = 0; ; ++iter_idx) {
    // P[i] = \sum_j beta[j] * X[i][j];
    // categorical optimization possible here
    ColMajorFmatrixVectorMultiplyStrided(xx, beta, sample_ct, sample_ctav, predictor_ct, pp);
    if (sample_offsets) {
      AddFVec(sample_offsets, sample_ctav, pp);
    }
    // P[i] = 1 / (1 + exp(-P[i]));
    LogisticSseF(sample_ct, pp);
    double loglik;
    if (ComputeLoglikCheckedF(yy, pp, sample_ct, &loglik)) {
      return 1;
    }
    // V[i] = P[i] * (1 - P[i]);
    ComputeVF(pp, sample_ct, vv);
    // P[i] -= Y[i] NOT done here

    // hessian = X diag(V) X'
    // note that only lower triangle is filled here
    ComputeHessianF(xx, vv, sample_ct, predictor_ct, hh0);

    // we shouldn't need to compute the log directly, since underflow <->
    // regression failure, right?  check this.
    if (InvertSymmdefFmatrixFirstHalf(predictor_ct, predictor_ctav, hh0, half_inverted_buf, inv_1d_buf, dbl_2d_buf)) {
      return 1;
    }
    double dethh = HalfSymmInvertedDet(half_inverted_buf, inv_1d_buf, predictor_ct, predictor_ct);
    loglik += 0.5 * log(dethh);

    InvertSymmdefFmatrixSecondHalf(predictor_ct, predictor_ctav, half_inverted_buf, hh0, inv_1d_buf, dbl_2d_buf);

    // bugfix (13 Oct 2017): trailing elements of hh0[] rows can't be arbitrary
    // for later MultMatrixDxnVectNF() call
    ReflectFmatrix0(predictor_ct, predictor_ctav, hh0);

    FirthComputeHdiagWeightsF(yy, xx, pp, hh0, vv, predictor_ct, predictor_ctav, sample_ct, sample_ctav, hdiag, ww, tmpnxk_buf);

    // gradient = X' W
    // categorical optimization possible here
    MultMatrixDxnVectNF(xx, ww, sample_ct, predictor_ct, ustar);
    if (iter_idx) {
      float ustar_max = 0.0;
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        const float abs_ustar_cur = fabsf(ustar[pred_uidx]);
        if (abs_ustar_cur > ustar_max) {
          ustar_max = abs_ustar_cur;
        }
      }
      const double loglik_change = loglik - loglik_old;
      if ((delta_max <= xconv) && (ustar_max < gconv) && (loglik_change < lconv)) {
        return 0;
      }
      if (iter_idx > max_iter) {
        *is_unfinished_ptr = 1;
        return 0;
      }
    }
    loglik_old = loglik;

    FirthComputeSecondWeightsF(hdiag, vv, sample_ct, sample_ctav, ww);
    ComputeHessianF(xx, ww, sample_ct, predictor_ct, hh);
    if (InvertSymmdefFmatrixFirstHalf(predictor_ct, predictor_ctav, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf)) {
      return 1;
    }
    InvertSymmdefFmatrixSecondHalf(predictor_ct, predictor_ctav, half_inverted_buf, hh, inv_1d_buf, dbl_2d_buf);
    ReflectFmatrix0(predictor_ct, predictor_ctav, hh);

    // delta := hh * ustar (note that hh is inverted already)
    MultMatrixDxnVectNF(hh, ustar, predictor_ct, predictor_ct, delta);

    delta_max = 0.0;
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const float abs_delta_cur = fabsf(delta[pred_uidx]);
      if (abs_delta_cur > delta_max) {
        delta_max = abs_delta_cur;
      }
    }
    const float maxstep = 5.0;
    if (delta_max > maxstep) {
      const float scaling_factor = maxstep / delta_max;
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        delta[pred_uidx] *= scaling_factor;
      }
      delta_max = maxstep;
    }
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      beta[pred_uidx] += delta[pred_uidx];
    }
  }
}

BoolErr FirthRegressionResidualizedF(const float* yy, const float* xx, const uintptr_t* sample_nm, const CcResidualizeCtx* cc_residualize, uint32_t nm_sample_ct, uint32_t orig_predictor_ct, float* beta, uint32_t* is_unfinished_ptr, float* hh, double* half_inverted_buf, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, float* pp, float* vv, float* ustar, float* delta, float* hdiag, float* ww, float* hh0_buf, float* tmpnxk_buf, float* mean_centered_pmaj_buf, float* sample_offsets_buf) {
  // todo: deduplicate with LogisticRegressionResidualizedF()
  const uintptr_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kFloatPerFVec);
  const uint32_t domdev_present_p1 = cc_residualize->domdev_present_p1;
  for (uint32_t geno_idx = 0; geno_idx != domdev_present_p1; ++geno_idx) {
    CopyAndMeanCenterF(&(xx[(geno_idx + 1) * nm_sample_ctav]), nm_sample_ct, &(mean_centered_pmaj_buf[geno_idx * nm_sample_ctav]));
  }
  const uint32_t prefitted_pred_ct = cc_residualize->prefitted_pred_ct;
  const uint32_t orig_biallelic_predictor_ct = domdev_present_p1 + prefitted_pred_ct;
  const uint32_t extra_allele_ct = orig_predictor_ct - orig_biallelic_predictor_ct;
  for (uint32_t extra_allele_idx = 0; extra_allele_idx != extra_allele_ct; ++extra_allele_idx) {
    CopyAndMeanCenterF(&(xx[(extra_allele_idx + orig_biallelic_predictor_ct) * nm_sample_ctav]), nm_sample_ct, &(mean_centered_pmaj_buf[(extra_allele_idx + domdev_present_p1) * nm_sample_ctav]));
  }
  const float* sample_offsets = cc_residualize->firth_nm_sample_offsets_f;
  if (nm_sample_ct != cc_residualize->sample_ct) {
    uintptr_t sample_idx_base = 0;
    uintptr_t sample_nm_bits = sample_nm[0];
    for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(sample_nm, &sample_idx_base, &sample_nm_bits);
      sample_offsets_buf[uii] = sample_offsets[sample_idx];
    }
    const uint32_t remainder = (-nm_sample_ct) & (kFloatPerFVec - 1);
    ZeroFArr(remainder, &(sample_offsets_buf[nm_sample_ct]));
    sample_offsets = sample_offsets_buf;
  }
  const uint32_t regressed_predictor_ct = domdev_present_p1 + extra_allele_ct;
  const uint32_t regressed_predictor_ctav = RoundUpPow2(regressed_predictor_ct, kFloatPerFVec);
  if (FirthRegressionF(yy, mean_centered_pmaj_buf, sample_offsets, nm_sample_ct, regressed_predictor_ct, &(beta[1]), is_unfinished_ptr, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf, pp, vv, ustar, delta, hdiag, ww, hh0_buf, tmpnxk_buf)) {
    return 1;
  }
  // hh is shifted up and to the left from what the caller expects, due to the
  // missing intercept.  Correct that here.
  // bugfix (4 Sep 2021): Initially thought only bottom-left triangle mattered,
  // but that's not true for genotypic/hethom case.  Also, wider stride
  // expected when regressed_predictor_ct is an exact multiple of
  // kFloatPerFVec.
  const uint32_t expected_predictor_ctav = RoundUpPow2(regressed_predictor_ct + 1, kFloatPerFVec);
  for (uint32_t write_row_idx = regressed_predictor_ct; write_row_idx; --write_row_idx) {
    memcpy(&(hh[write_row_idx * expected_predictor_ctav + 1]), &(hh[(write_row_idx - 1) * regressed_predictor_ctav]), regressed_predictor_ct * sizeof(float));
  }
  return 0;
}

uintptr_t GetLogisticWorkspaceSizeF(uint32_t sample_ct, uint32_t biallelic_predictor_ct, uint32_t domdev_present_p1, uint32_t max_extra_allele_ct, uint32_t constraint_ct, uint32_t xmain_ct, uint32_t gcount_cc, uint32_t is_sometimes_firth, uint32_t is_cc_residualize) {
  // sample_ctav * max_predictor_ct < 2^31, and sample_ct >=
  // biallelic_predictor_ct, so no overflows?
  // could round everything up to multiples of 16 instead of 64
  const uint32_t max_predictor_ct = biallelic_predictor_ct + max_extra_allele_ct;
  const uint32_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  const uint32_t max_predictor_ctav = RoundUpPow2(max_predictor_ct, kFloatPerFVec);
  // sample_nm, pheno_cc_nm, tmp_nm = sample_ctl words
  uintptr_t workspace_size = 3 * RoundUpPow2(BitCtToWordCt(sample_ct) * sizeof(intptr_t), kCacheline);

  // yy = sample_ctav floats
  workspace_size += RoundUpPow2(sample_ctav * sizeof(float), kCacheline);

  // xx = (max_predictor_ct + main_mutated + main_omitted) * sample_ctav floats
  workspace_size += RoundUpPow2((max_predictor_ct + xmain_ct) * sample_ctav * sizeof(float), kCacheline);

  // hh = max_predictor_ct * max_predictor_ctav floats
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ctav * sizeof(float), kCacheline);

  // pp, vv = sample_ctav floats
  workspace_size += 2 * RoundUpPow2(sample_ctav * sizeof(float), kCacheline);

  // coef, grad, dcoef = max_predictor_ctav floats
  workspace_size += 3 * RoundUpPow2(max_predictor_ctav * sizeof(float), kCacheline);

  // ll = max_predictor_ct * max_predictor_ctav floats
  // (technically not needed in pure-Firth case)
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ctav * sizeof(float), kCacheline);

  // semicomputed_biallelic_xtx
  workspace_size += RoundUpPow2(biallelic_predictor_ct * biallelic_predictor_ct * sizeof(double), kCacheline);

  // semicomputed_biallelic_corr_matrix
  workspace_size += RoundUpPow2((biallelic_predictor_ct - 1) * (biallelic_predictor_ct - 1) * sizeof(double), kCacheline);

  // semicomputed_biallelic_inv_corr_sqrts
  workspace_size += RoundUpPow2(biallelic_predictor_ct * sizeof(double), kCacheline);

  // inv_1d_buf
  workspace_size += RoundUpPow2(max_predictor_ct * kMatrixInvertBuf1CheckedAlloc, kCacheline);

  // dbl_2d_buf = max_predictor_ct * max_predictor_ct floats, or VIF/Firth
  //              (which can be max_predictor_ct * 7 floats).
  workspace_size += RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);

  // a1_dosages, a1_case_dosages
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(double) * 2, kCacheline);

  // machr2_dosage_sums, machr2_dosage_ssqs
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(uint64_t) * 2, kCacheline);

  if (gcount_cc && max_extra_allele_ct) {
    // case_one_cts, case_two_cts
    workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(int32_t) * 2, kCacheline);
  }

  // predictor_dotprod_buf
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ct * sizeof(float), kCacheline);

  const uintptr_t other_2d_byte_ct = max_predictor_ct * MAXV(max_predictor_ct, 3) * sizeof(double);
  // inverse_corr_buf/half_inverted_buf
  workspace_size += RoundUpPow2(other_2d_byte_ct, kCacheline);

  if (is_sometimes_firth) {
    // hdiag, ww = sample_ctav floats
    workspace_size += 2 * RoundUpPow2(sample_ctav * sizeof(float), kCacheline);

    // hh0 = max_predictor_ct * max_predictor_ctav floats
    workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ctav * sizeof(float), kCacheline);

    // tmpnxk_buf = max_predictor_ct * sample_ctav floats
    workspace_size += RoundUpPow2(max_predictor_ct * sample_ctav * sizeof(float), kCacheline);
  }
  if (is_cc_residualize) {
    // mean_centered_pmaj_buf = (domdev_present_p1 + max_extra_allele_ct) *
    //   sample_ctav floats
    workspace_size += RoundUpPow2((domdev_present_p1 + max_extra_allele_ct) * sample_ctav * sizeof(float), kCacheline);

    // sample_offsets_buf
    workspace_size += RoundUpPow2(sample_ctav * sizeof(float), kCacheline);
  }
  if (constraint_ct) {
    // tmphxs_buf, h_transpose_buf = constraint_ct * max_predictor_ctav floats
    workspace_size += 2 * RoundUpPow2(constraint_ct * max_predictor_ctav * sizeof(float), kCacheline);

    // inner_buf = constraint_ct * constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * constraint_ct * sizeof(float), kCacheline);

    // outer_buf = constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * sizeof(float), kCacheline);

    // constraints_con_major = constraint_ct * max_predictor_ct
    workspace_size += RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(float), kCacheline);
  }
  return workspace_size;
}

static const float kSmallFloats[4] = {0.0, 1.0, 2.0, 3.0};

THREAD_FUNC_DECL GlmLogisticThreadF(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GlmLogisticCtx* ctx = S_CAST(GlmLogisticCtx*, arg->sharedp->context);
  GlmCtx* common = ctx->common;

  PgenReader* pgrp = common->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = common->genovecs[tidx];
  pgv.dosage_present = nullptr;
  pgv.dosage_main = nullptr;
  if (common->dosage_presents) {
    pgv.dosage_present = common->dosage_presents[tidx];
    pgv.dosage_main = common->dosage_mains[tidx];
  }
  unsigned char* workspace_buf = common->workspace_bufs[tidx];
  const uintptr_t* variant_include = common->variant_include;
  const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
  const AlleleCode* omitted_alleles = common->omitted_alleles;
  const uintptr_t* sex_male_collapsed = common->sex_male_collapsed;
  const ChrInfo* cip = common->cip;
  const uint32_t* subset_chr_fo_vidx_start = common->subset_chr_fo_vidx_start;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const GlmFlags glm_flags = common->glm_flags;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t is_sometimes_firth = !(glm_flags & kfGlmNoFirth);
  const uint32_t is_always_firth = (glm_flags / kfGlmFirth) & 1;
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t model_hetonly = (glm_flags / kfGlmHetonly) & 1;
  const uint32_t joint_genotypic = (glm_flags / kfGlmGenotypic) & 1;
  const uint32_t joint_hethom = (glm_flags / kfGlmHethom) & 1;
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uint32_t domdev_present = joint_genotypic || joint_hethom;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t reported_pred_uidx_start = 1 - include_intercept;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const uintptr_t max_reported_test_ct = common->max_reported_test_ct;
  const uintptr_t local_covar_ct = common->local_covar_ct;
  const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
  // bugfix (20 Mar 2020): Also need to exclude dominant/recessive.
  const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!model_dominant) && (!model_recessive) && (!model_hetonly) && (!common->tests_flag) && (!add_interactions);
  uintptr_t max_sample_ct = MAXV(common->sample_ct, common->sample_ct_x);
  if (max_sample_ct < common->sample_ct_y) {
    max_sample_ct = common->sample_ct_y;
  }
  SetPgvThreadMhcNull(max_sample_ct, tidx, common->thread_mhc, &pgv);
  pgv.patch_01_ct = 0;
  pgv.patch_10_ct = 0;
  pgv.multidosage_sample_ct = 0;
  uint32_t variant_idx_offset = 0;
  uint32_t allele_ct = 2;
  uint32_t omitted_allele_idx = 0;
  uint32_t extra_regression_ct = 0;
  double main_dosage_sum = 0.0;
  double main_dosage_ssq = 0.0;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);

    double* beta_se_iter = common->block_beta_se;
    uintptr_t allele_bidx = variant_bidx;
    if (max_extra_allele_ct) {
      allele_bidx = variant_bidx + CountExtraAlleles(variant_include, allele_idx_offsets, common->read_variant_uidx_starts[0], common->read_variant_uidx_starts[tidx], 0);
    }
    if (beta_se_multiallelic_fused) {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * variant_bidx]);
    } else {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * allele_bidx]);
    }

    LogisticAuxResult* block_aux_iter = &(ctx->block_aux[allele_bidx]);
    const float* local_covars_iter = nullptr;
    if (local_covar_ct) {
      // &(nullptr[0]) is okay in C++, but undefined in C
      local_covars_iter = &(ctx->local_covars_vcmaj_f[parity][variant_bidx * max_sample_ct * local_covar_ct]);
    }
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = LastLeqU32(subset_chr_fo_vidx_start, 0, cip->chr_ct, variant_idx);
      // const uint32_t chr_fo_idx = LowerBoundNonemptyU32(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      uint32_t cur_variant_bidx_end = subset_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
      if (cur_variant_bidx_end > variant_bidx_end) {
        cur_variant_bidx_end = variant_bidx_end;
      }
      // "regular" = not all-female special case.
      const uint32_t is_haploid = IsSet(cip->haploid_mask, chr_idx);
      const uint32_t is_regular_x = is_haploid && (chr_idx == x_code);
      const uint32_t is_y = (chr_idx == y_code);
      const uint32_t is_nonx_haploid = is_haploid && (!is_regular_x);
      const uintptr_t* cur_sample_include;
      const uint32_t* cur_sample_include_cumulative_popcounts;
      const uintptr_t* cur_pheno_cc;
      const uintptr_t* cur_gcount_case_interleaved_vec;
      const float* cur_pheno;
      const RegressionNmPrecomp* nm_precomp;
      const float* cur_covars_cmaj;
      const uintptr_t* cur_parameter_subset;
      const uintptr_t* cur_joint_test_params;
      const CcResidualizeCtx* cur_cc_residualize;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      uint32_t cur_constraint_ct;
      uint32_t cur_is_always_firth;
      if (is_y && common->sample_include_y) {
        cur_sample_include = common->sample_include_y;
        cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
        cur_pheno_cc = ctx->pheno_y_cc;
        cur_gcount_case_interleaved_vec = ctx->gcount_case_interleaved_vec_y;
        cur_pheno = ctx->pheno_y_f;
        nm_precomp = common->nm_precomp_y;
        cur_covars_cmaj = ctx->covars_cmaj_y_f;
        cur_parameter_subset = common->parameter_subset_y;
        cur_joint_test_params = common->joint_test_params_y;
        cur_cc_residualize = ctx->cc_residualize_y;
        cur_sample_ct = common->sample_ct_y;
        cur_covar_ct = common->covar_ct_y;
        cur_constraint_ct = common->constraint_ct_y;
        cur_is_always_firth = is_always_firth || ctx->separation_found_y;
      } else if (is_regular_x && common->sample_include_x) {
        cur_sample_include = common->sample_include_x;
        cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
        cur_pheno_cc = ctx->pheno_x_cc;
        cur_gcount_case_interleaved_vec = ctx->gcount_case_interleaved_vec_x;
        cur_pheno = ctx->pheno_x_f;
        nm_precomp = common->nm_precomp_x;
        cur_covars_cmaj = ctx->covars_cmaj_x_f;
        cur_parameter_subset = common->parameter_subset_x;
        cur_joint_test_params = common->joint_test_params_x;
        cur_cc_residualize = ctx->cc_residualize_x;
        cur_sample_ct = common->sample_ct_x;
        cur_covar_ct = common->covar_ct_x;
        cur_constraint_ct = common->constraint_ct_x;
        cur_is_always_firth = is_always_firth || ctx->separation_found_x;
      } else {
        cur_sample_include = common->sample_include;
        cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
        cur_pheno_cc = ctx->pheno_cc;
        cur_gcount_case_interleaved_vec = ctx->gcount_case_interleaved_vec;
        cur_pheno = ctx->pheno_f;
        nm_precomp = common->nm_precomp;
        cur_covars_cmaj = ctx->covars_cmaj_f;
        cur_parameter_subset = common->parameter_subset;
        cur_joint_test_params = common->joint_test_params;
        cur_cc_residualize = ctx->cc_residualize;
        cur_sample_ct = common->sample_ct;
        cur_covar_ct = common->covar_ct;
        cur_constraint_ct = common->constraint_ct;
        cur_is_always_firth = is_always_firth || ctx->separation_found;
      }
      const uint32_t sample_ctl = BitCtToWordCt(cur_sample_ct);
      const uint32_t sample_ctav = RoundUpPow2(cur_sample_ct, kFloatPerFVec);
      const uint32_t cur_case_ct = PopcountWords(cur_pheno_cc, sample_ctl);
      const uint32_t cur_biallelic_predictor_ct_base = 2 + domdev_present + cur_covar_ct * (1 + add_interactions * domdev_present_p1);
      uint32_t cur_biallelic_predictor_ct = cur_biallelic_predictor_ct_base;
      uint32_t literal_covar_ct = cur_covar_ct;
      if (cur_parameter_subset) {
        cur_biallelic_predictor_ct = PopcountWords(cur_parameter_subset, BitCtToWordCt(cur_biallelic_predictor_ct_base));
        literal_covar_ct = PopcountBitRange(cur_parameter_subset, 2 + domdev_present, 2 + domdev_present + cur_covar_ct);
      }
      const uint32_t max_predictor_ct = cur_biallelic_predictor_ct + max_extra_allele_ct;
      const uint32_t max_predictor_ctav = RoundUpPow2(max_predictor_ct, kFloatPerFVec);
      uint32_t reported_pred_uidx_biallelic_end;
      if (hide_covar) {
        if (!cur_parameter_subset) {
          reported_pred_uidx_biallelic_end = 2 + domdev_present;
        } else {
          reported_pred_uidx_biallelic_end = 1 + IsSet(cur_parameter_subset, 1) + domdev_present;
        }
      } else {
        reported_pred_uidx_biallelic_end = cur_biallelic_predictor_ct;
      }
      // nm_predictors_pmaj_buf may require up to two extra columns omitted
      // from the main regression.
      // 1. In the multiallelic dominant/recessive/hetonly/hethom cases, the
      //    original genotype column does not appear in the regression, and
      //    we'd rather not reconstruct it from genovec, etc. when we need to
      //    swap it out for another allele, so we keep the original genotype in
      //    an extra column.
      //    To reduce code bloat, we now handle the biallelic cases in the same
      //    way; this is one of the more peripheral code paths so adding more
      //    complexity to speed it up is less justifiable.
      // 2. If --parameters excludes the main (possibly
      //    dominant/recessive/hetonly) genotype column but does care about an
      //    interaction, we want a copy of what the main genotype column's
      //    contents would have been to refer to.
      const uint32_t main_omitted = cur_parameter_subset && (!IsSet(cur_parameter_subset, 1));
      const uint32_t main_mutated = model_dominant || model_recessive || model_hetonly || joint_hethom;
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* pheno_cc_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* tmp_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      float* nm_pheno_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
      float* nm_predictors_pmaj_buf = S_CAST(float*, arena_alloc_raw_rd((max_predictor_ct + main_mutated + main_omitted) * sample_ctav * sizeof(float), &workspace_iter));
      float* coef_return = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ctav * sizeof(float), &workspace_iter));
      float* hh_return = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ctav * sizeof(float), &workspace_iter));
      float* pp_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
      float* sample_variance_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
      float* gradient_buf = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ctav * sizeof(float), &workspace_iter));
      float* dcoef_buf = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ctav * sizeof(float), &workspace_iter));
      float* cholesky_decomp_return = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ctav * sizeof(float), &workspace_iter));

      double* semicomputed_biallelic_xtx = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));
      // currently overallocates
      double* semicomputed_biallelic_corr_matrix = S_CAST(double*, arena_alloc_raw_rd((cur_biallelic_predictor_ct - 1) * (cur_biallelic_predictor_ct - 1) * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_inv_corr_sqrts = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));

      MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, arena_alloc_raw_rd(max_predictor_ct * kMatrixInvertBuf1CheckedAlloc, &workspace_iter));
      const uintptr_t dbl_2d_byte_ct = RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);
      double* dbl_2d_buf = S_CAST(double*, arena_alloc_raw(dbl_2d_byte_ct, &workspace_iter));
      double* a1_dosages = S_CAST(double*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(double) * 2, &workspace_iter));
      double* a1_case_dosages = &(a1_dosages[max_extra_allele_ct + 2]);
      uint64_t* machr2_dosage_sums = S_CAST(uint64_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(uint64_t) * 2, &workspace_iter));
      uint64_t* machr2_dosage_ssqs = &(machr2_dosage_sums[max_extra_allele_ct + 2]);
      uint32_t* case_one_cts = nullptr;
      uint32_t* case_two_cts = nullptr;
      if (cur_gcount_case_interleaved_vec && max_extra_allele_ct) {
        case_one_cts = S_CAST(uint32_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(int32_t) * 2, &workspace_iter));
        case_two_cts = &(case_one_cts[max_extra_allele_ct + 2]);
      }
      float* predictor_dotprod_buf = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ct * sizeof(float), &workspace_iter));
      const uintptr_t other_2d_byte_ct = RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 3) * sizeof(double), kCacheline);
      double* inverse_corr_buf = S_CAST(double*, arena_alloc_raw(other_2d_byte_ct, &workspace_iter));

      // these could use the same memory, but not a big deal, use the less
      // bug-prone approach for now
      // Firth-only
      float* hdiag_buf = nullptr;
      float* score_buf = nullptr;
      float* hh0_buf = nullptr;
      float* tmpnxk_buf = nullptr;
      if (is_sometimes_firth) {
        hdiag_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
        score_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
        hh0_buf = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ctav * sizeof(float), &workspace_iter));
        tmpnxk_buf = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ct * sample_ctav * sizeof(float), &workspace_iter));
      }
      float* mean_centered_pmaj_buf = nullptr;
      float* sample_offsets_buf = nullptr;
      if (cur_cc_residualize) {
        mean_centered_pmaj_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float) * (domdev_present_p1 + max_extra_allele_ct), &workspace_iter));
        sample_offsets_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
      }

      // joint test only
      float* tmphxs_buf = nullptr;
      float* h_transpose_buf = nullptr;
      float* inner_buf = nullptr;
      float* outer_buf = nullptr;
      float* cur_constraints_con_major = nullptr;
      if (cur_constraint_ct) {
        tmphxs_buf = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ctav * sizeof(float), &workspace_iter));
        h_transpose_buf = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ctav * sizeof(float), &workspace_iter));
        inner_buf = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * cur_constraint_ct * sizeof(float), &workspace_iter));
        outer_buf = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * sizeof(float), &workspace_iter));
        // bugfix (27 Jan 2019): forgot sizeof(float) here
        cur_constraints_con_major = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(float), &workspace_iter));
        ZeroFArr(cur_constraint_ct * max_predictor_ct, cur_constraints_con_major);
        const uint32_t first_joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
        cur_constraints_con_major[first_joint_test_idx] = 1.0;
        // Rest of this matrix must be updated later, since cur_predictor_ct
        // changes at multiallelic variants.
      }
      assert(S_CAST(uintptr_t, workspace_iter - workspace_buf) == GetLogisticWorkspaceSizeF(cur_sample_ct, cur_biallelic_predictor_ct, domdev_present_p1, max_extra_allele_ct, cur_constraint_ct, main_mutated + main_omitted, cur_gcount_case_interleaved_vec != nullptr, is_sometimes_firth, cur_cc_residualize != nullptr));
      const double cur_sample_ct_recip = 1.0 / u31tod(cur_sample_ct);
      const double cur_sample_ct_m1_recip = 1.0 / u31tod(cur_sample_ct - 1);
      const double* corr_inv = nullptr;
      if (nm_precomp) {
        memcpy(semicomputed_biallelic_xtx, nm_precomp->xtx_image, cur_biallelic_predictor_ct * cur_biallelic_predictor_ct * sizeof(double));
        corr_inv = nm_precomp->corr_inv;
        const uintptr_t nongeno_pred_ct = cur_biallelic_predictor_ct - domdev_present - 2;
        const uintptr_t nonintercept_biallelic_pred_ct = cur_biallelic_predictor_ct - 1;
        memcpy(semicomputed_biallelic_corr_matrix, nm_precomp->corr_image, nonintercept_biallelic_pred_ct * nonintercept_biallelic_pred_ct * sizeof(double));
        memcpy(&(semicomputed_biallelic_inv_corr_sqrts[domdev_present_p1]), nm_precomp->corr_inv_sqrts, nongeno_pred_ct * sizeof(double));
      }
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      // when this is set, the last fully-processed variant had no missing
      // genotypes, and if the current variant also has no missing genotypes we
      // may be able to skip reinitialization of most of
      // nm_predictors_pmaj_buf.
      // (todo: do we want to track prev_biallelic_nm?)
      uint32_t prev_nm = 0;

      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        if (allele_idx_offsets) {
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
          if (!beta_se_multiallelic_fused) {
            extra_regression_ct = allele_ct - 2;
          }
        }
        const uint32_t allele_ct_m2 = allele_ct - 2;
        const uint32_t expected_predictor_ct = cur_biallelic_predictor_ct + allele_ct_m2;
        PglErr reterr;
        if (!allele_ct_m2) {
          reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
        } else {
          reterr = PgrGetMD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, &pgv);
          // todo: proper multiallelic dosage support
        }
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto GlmLogisticThreadF_err;
        }
        ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
        GenoarrCountFreqsUnsafe(pgv.genovec, cur_sample_ct, genocounts);
        uint32_t missing_ct = genocounts[3];
        if (!missing_ct) {
          SetAllBits(cur_sample_ct, sample_nm);
        } else {
          GenoarrToNonmissing(pgv.genovec, cur_sample_ct, sample_nm);
          if (pgv.dosage_ct) {
            BitvecOr(pgv.dosage_present, sample_ctl, sample_nm);
            missing_ct = cur_sample_ct - PopcountWords(sample_nm, sample_ctl);
          }
        }
        if (omitted_alleles) {
          omitted_allele_idx = omitted_alleles[variant_uidx];
        }
        // Once sizeof(AlleleCode) > 1, we probably want to allocate this from
        // g_bigstack instead of the thread stack.
        uintptr_t const_alleles[DivUp(kPglMaxAlleleCt, kBitsPerWord)];
        const uint32_t allele_ctl = DivUp(allele_ct, kBitsPerWord);
        ZeroWArr(allele_ctl, const_alleles);
        const uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
        const uint32_t nm_sample_ctl = BitCtToWordCt(nm_sample_ct);
        const uint32_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kFloatPerFVec);
        const uint32_t nm_sample_ct_rem = nm_sample_ctav - nm_sample_ct;
        // first predictor column: intercept
        if (!prev_nm) {
          FillFVec(nm_sample_ct, S_CAST(float, 1.0), nm_predictors_pmaj_buf);
        }
        // second predictor column: genotype
        float* genotype_vals = &(nm_predictors_pmaj_buf[nm_sample_ctav]);
        if (main_mutated || main_omitted) {
          // bugfix (9 Oct 2024)
          ZeroFArr(nm_sample_ct_rem, &(genotype_vals[nm_sample_ct]));
          genotype_vals = &(nm_predictors_pmaj_buf[expected_predictor_ct * nm_sample_ctav]);
        }
        CopyBitarrSubset(cur_pheno_cc, sample_nm, nm_sample_ct, pheno_cc_nm);
        const uint32_t nm_case_ct = PopcountWords(pheno_cc_nm, nm_sample_ctl);
        float* multi_start = nullptr;
        if (!allele_ct_m2) {
          if (omitted_allele_idx) {
            GenovecInvertUnsafe(cur_sample_ct, pgv.genovec);
            // ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
            }
            const uint32_t uii = genocounts[0];
            genocounts[0] = genocounts[2];
            genocounts[2] = uii;
          }
          uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
          uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
          if (!missing_ct) {
            GenoarrLookup16x4bx2(pgv.genovec, kSmallFloatPairs, nm_sample_ct, genotype_vals);
            if (pgv.dosage_ct) {
              uintptr_t sample_idx_base = 0;
              uintptr_t dosage_present_bits = pgv.dosage_present[0];
              for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                const uintptr_t sample_idx = BitIter1(pgv.dosage_present, &sample_idx_base, &dosage_present_bits);
                const uint32_t dosage_val = pgv.dosage_main[dosage_idx];
                // 32768 -> 2, 16384 -> 1, 0 -> 0
                genotype_vals[sample_idx] = kRecipDosageMidf * u31tof(dosage_val);
                dosage_sum += dosage_val;
                dosage_ssq += dosage_val * dosage_val;
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_idx);
                if (cur_geno && (cur_geno != 3)) {
                  const uintptr_t prev_val = cur_geno * kDosageMid;
                  dosage_sum -= prev_val;
                  dosage_ssq -= prev_val * prev_val;
                }
              }
            }
          } else {
            if (!pgv.dosage_ct) {
              GenoarrToFloatsRemoveMissing(pgv.genovec, kSmallFloats, cur_sample_ct, genotype_vals);
            } else {
              uintptr_t sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              uint32_t dosage_idx = 0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_midx);
                float cur_val;
                if (IsSet(pgv.dosage_present, sample_midx)) {
                  const uint32_t dosage_val = pgv.dosage_main[dosage_idx++];
                  cur_val = kRecipDosageMidf * u31tof(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                } else {
                  // cur_geno != 3 guaranteed
                  cur_val = kSmallFloats[cur_geno];
                }
                genotype_vals[sample_idx] = cur_val;
              }
            }
          }
          // Check for constant genotype column.
          // (Technically, we should recheck later in the chrX no-sex-covariate
          // --xchr-model 1 corner case.)
          if (!pgv.dosage_ct) {
            if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
              // bugfix (28 Mar 2020): didn't set the bit that actually
              // mattered last week...
              const_alleles[0] = 3;
            }
          } else if (pgv.dosage_ct == nm_sample_ct) {
            if (DosageIsConstant(dosage_sum, dosage_ssq, nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          }
          machr2_dosage_sums[1 - omitted_allele_idx] = dosage_sum;
          machr2_dosage_ssqs[1 - omitted_allele_idx] = dosage_ssq;
          machr2_dosage_sums[omitted_allele_idx] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[omitted_allele_idx] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
          if (cur_gcount_case_interleaved_vec) {
            // gcountcc
            STD_ARRAY_REF(uint32_t, 6) cur_geno_hardcall_cts = block_aux_iter->geno_hardcall_cts;
            GenoarrCountSubsetFreqs(pgv.genovec, cur_gcount_case_interleaved_vec, cur_sample_ct, cur_case_ct, R_CAST(STD_ARRAY_REF(uint32_t, 4), cur_geno_hardcall_cts));
            for (uint32_t geno_hardcall_idx = 0; geno_hardcall_idx != 3; ++geno_hardcall_idx) {
              cur_geno_hardcall_cts[3 + geno_hardcall_idx] = genocounts[geno_hardcall_idx] - cur_geno_hardcall_cts[geno_hardcall_idx];
            }
          }
        } else {
          // multiallelic.
          // Update (18 Mar 2020): If some but not all alleles have constant
          // dosages, we remove just those alleles from the regressions;
          // trim-alts is not necessary to see what's going on with the other
          // alleles.  To reduce parsing complexity, the number of output lines
          // is not affected by this; the ones corresponding to the constant
          // alleles have NA values.

          // dosage_ct == 0 temporarily guaranteed if we reach here.
          assert(!pgv.dosage_ct);
          multi_start = &(nm_predictors_pmaj_buf[(expected_predictor_ct - allele_ct_m2) * nm_sample_ctav]);
          ZeroU64Arr(allele_ct, machr2_dosage_sums);
          ZeroU64Arr(allele_ct, machr2_dosage_ssqs);
          // postpone multiply for now, since no multiallelic dosages
          // Use sums as ones[] and ssqs as twos[] for rarealts; transform to
          // actual sums/ssqs later.
          machr2_dosage_sums[0] = genocounts[1];
          machr2_dosage_ssqs[0] = genocounts[0];
          if (omitted_allele_idx) {
            // Main genotype column starts as REF.
            if (!missing_ct) {
              GenoarrLookup16x4bx2(pgv.genovec, kSmallInvFloatPairs, nm_sample_ct, genotype_vals);
            } else {
              GenoarrToFloatsRemoveMissing(pgv.genovec, kSmallInvFloats, cur_sample_ct, genotype_vals);
            }
          }
          uint32_t rare_allele_ct = allele_ct_m2;
          float* alt1_start = nullptr;
          float* rarealt_start = multi_start;
          if (omitted_allele_idx != 1) {
            if (omitted_allele_idx) {
              alt1_start = multi_start;
              ZeroFArr(nm_sample_ct_rem, &(alt1_start[nm_sample_ct]));
              rarealt_start = &(rarealt_start[nm_sample_ctav]);
              --rare_allele_ct;
            } else {
              alt1_start = genotype_vals;
            }
            if (!missing_ct) {
              GenoarrLookup16x4bx2(pgv.genovec, kSmallFloatPairs, nm_sample_ct, alt1_start);
            } else {
              GenoarrToFloatsRemoveMissing(pgv.genovec, kSmallFloats, cur_sample_ct, alt1_start);
            }
          }
          ZeroFArr(rare_allele_ct * nm_sample_ctav, rarealt_start);
          if (pgv.patch_01_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_01_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_01_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ctav + sample_idx] = 1.0;
                alt1_start[sample_idx] = 0.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ctav + sample_idx] = 1.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                alt1_start[sample_idx] = 0.0;
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                machr2_dosage_sums[allele_code] += 1;
                if (allele_code == omitted_allele_idx) {
                  continue;
                }
                const uint32_t cur_col = allele_code - 2 - (allele_code > omitted_allele_idx);
                rarealt_start[cur_col * nm_sample_ctav + sample_idx] = 1.0;
              }
            }
          }
          uintptr_t alt1_het_ct = genocounts[1] - pgv.patch_01_ct;
          if (pgv.patch_10_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_10_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_10_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] = 2.0;
                  alt1_start[sample_idx] = 0.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ctav + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] += S_CAST(float, 1.0);
                    alt1_start[sample_idx] = 0.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] = 2.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ctav + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] += S_CAST(float, 1.0);
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  machr2_dosage_ssqs[ac0] += 1;
                  alt1_start[sample_idx] = 0.0;
                  if (ac0 != omitted_allele_idx) {
                    const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                    rarealt_start[ac0_col * nm_sample_ctav + sample_idx] = 2.0;
                  }
                } else {
                  machr2_dosage_sums[ac1] += 1;
                  if (ac1 != omitted_allele_idx) {
                    const uint32_t ac1_col = ac1 - 2 - (ac1 > omitted_allele_idx);
                    rarealt_start[ac1_col * nm_sample_ctav + sample_idx] = 1.0;
                  }
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    machr2_dosage_sums[ac0] += 1;
                    alt1_start[sample_idx] = 0.0;
                    if (ac0 != omitted_allele_idx) {
                      const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                      rarealt_start[ac0_col * nm_sample_ctav + sample_idx] += S_CAST(float, 1.0);
                    }
                  }
                }
              }
            }
          }
          machr2_dosage_sums[1] = alt1_het_ct;
          machr2_dosage_ssqs[1] = genocounts[2] - pgv.patch_10_ct;
          if (cur_gcount_case_interleaved_vec) {
            // gcountcc.  Need case-specific one_cts and two_cts for each
            // allele.
            STD_ARRAY_DECL(uint32_t, 4, case_hardcall_cts);
            GenoarrCountSubsetFreqs(pgv.genovec, cur_gcount_case_interleaved_vec, cur_sample_ct, cur_case_ct, case_hardcall_cts);
            ZeroU32Arr(allele_ct, case_one_cts);
            ZeroU32Arr(allele_ct, case_two_cts);
            uint32_t case_alt1_het_ct = case_hardcall_cts[1];
            case_one_cts[0] = case_alt1_het_ct;
            case_two_cts[0] = case_hardcall_cts[0];
            if (pgv.patch_01_ct) {
              uintptr_t sample_widx = 0;
              uintptr_t cur_bits = pgv.patch_01_set[0];
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t lowbit = BitIter1y(pgv.patch_01_set, &sample_widx, &cur_bits);
                if (cur_pheno_cc[sample_widx] & lowbit) {
                  const uint32_t allele_code = pgv.patch_01_vals[uii];
                  case_one_cts[allele_code] += 1;
                }
              }
              for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
                case_alt1_het_ct -= case_one_cts[allele_idx];
              }
            }
            uint32_t case_alt1_hom_ct = case_hardcall_cts[2];
            if (pgv.patch_10_ct) {
              uintptr_t sample_widx = 0;
              uintptr_t cur_bits = pgv.patch_10_set[0];
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t lowbit = BitIter1y(pgv.patch_10_set, &sample_widx, &cur_bits);
                if (cur_pheno_cc[sample_widx] & lowbit) {
                  const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                  const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                  --case_alt1_hom_ct;
                  if (ac0 == ac1) {
                    case_two_cts[ac0] += 1;
                  } else {
                    case_one_cts[ac1] += 1;
                    if (ac0 == 1) {
                      ++case_alt1_het_ct;
                    } else {
                      case_one_cts[ac0] += 1;
                    }
                  }
                }
              }
            }
            case_one_cts[1] = case_alt1_het_ct;
            case_two_cts[1] = case_alt1_hom_ct;
            uint32_t nonomitted_allele_idx = 0;
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if (allele_idx == omitted_allele_idx) {
                continue;
              }
              const uint32_t one_ct = machr2_dosage_sums[allele_idx];
              const uint32_t two_ct = machr2_dosage_ssqs[allele_idx];
              const uint32_t case_one_ct = case_one_cts[allele_idx];
              const uint32_t case_two_ct = case_two_cts[allele_idx];
              STD_ARRAY_REF(uint32_t, 6) dst = block_aux_iter[nonomitted_allele_idx].geno_hardcall_cts;
              dst[0] = nm_case_ct - case_one_ct - case_two_ct;
              dst[1] = case_one_ct;
              dst[2] = case_two_ct;
              dst[3] = nm_sample_ct - one_ct - two_ct - dst[0];
              dst[4] = one_ct - case_one_ct;
              dst[5] = two_ct - case_two_ct;
              ++nonomitted_allele_idx;
            }
          }
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            const uintptr_t one_ct = machr2_dosage_sums[allele_idx];
            const uintptr_t two_ct = machr2_dosage_ssqs[allele_idx];
            machr2_dosage_sums[allele_idx] = (one_ct + 2 * two_ct) * 0x4000LLU;
            machr2_dosage_ssqs[allele_idx] = (one_ct + 4LLU * two_ct) * 0x10000000LLU;
            if ((one_ct == nm_sample_ct) || (two_ct == nm_sample_ct) || ((!one_ct) && (!two_ct))) {
              SetBit(allele_idx, const_alleles);
            }
          }
        }
        ZeroFArr(nm_sample_ct_rem, &(genotype_vals[nm_sample_ct]));
        // usually need to save some of {sample_obs_ct, allele_obs_ct,
        // a1_dosage, case_allele_obs_ct, a1_case_dosage, mach_r2 even for
        // skipped variants
        // compute them all for now, could conditionally skip later
        uint32_t allele_obs_ct = nm_sample_ct * 2;
        uint32_t case_allele_obs_ct = nm_case_ct * 2;
        if (!is_regular_x) {
          if (is_nonx_haploid) {
            allele_obs_ct = nm_sample_ct;
            case_allele_obs_ct = nm_case_ct;
            // everything is on 0..1 scale, not 0..2
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              genotype_vals[sample_idx] *= S_CAST(float, 0.5);
            }
            const uint32_t high_ct = nm_sample_ct * allele_ct_m2;
            for (uint32_t uii = 0; uii != high_ct; ++uii) {
              multi_start[uii] *= S_CAST(float, 0.5);
            }
          }
        } else {
          CopyBitarrSubset(sex_male_collapsed, sample_nm, nm_sample_ct, tmp_nm);
          const uintptr_t* male_nm = tmp_nm;
          const uint32_t nm_male_ct = PopcountWords(male_nm, nm_sample_ctl);
          if (is_xchr_model_1) {
            // special case: multiply male values by 0.5
            uintptr_t sample_idx_base = 0;
            uintptr_t male_nm_bits = male_nm[0];
            for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
              const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
              genotype_vals[sample_idx] *= S_CAST(float, 0.5);
              // could insert multiallelic loop here isntead, but I'm guessing
              // that's worse due to locality of writes?
            }
            for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
              float* cur_start = &(multi_start[extra_allele_idx * nm_sample_ctav]);
              sample_idx_base = 0;
              male_nm_bits = male_nm[0];
              for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
                const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
                cur_start[sample_idx] *= S_CAST(float, 0.5);
              }
            }
            allele_obs_ct -= nm_male_ct;
            case_allele_obs_ct -= PopcountWordsIntersect(pheno_cc_nm, male_nm, nm_sample_ctl);
          }
        }
        const double mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, allele_ct);
        uint32_t nonomitted_allele_idx = 0;
        for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
          if (allele_idx == omitted_allele_idx) {
            continue;
          }

          float* geno_col = genotype_vals;
          if (allele_idx > (!omitted_allele_idx)) {
            geno_col = &(nm_predictors_pmaj_buf[(expected_predictor_ct - (allele_ct - allele_idx) + (allele_idx < omitted_allele_idx)) * nm_sample_ctav]);
          }
          double a1_dosage = u63tod(machr2_dosage_sums[allele_idx]) * kRecipDosageMid;
          if (is_xchr_model_1) {
            // ugh.
            a1_dosage = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              a1_dosage += S_CAST(double, geno_col[sample_idx]);
            }
          } else {
            if (is_nonx_haploid) {
              a1_dosage *= 0.5;
            }
          }
          a1_dosages[allele_idx] = a1_dosage;

          // todo: shortcut if gcountcc computed and no dosages
          double a1_case_dosage = 0.0;
          uintptr_t sample_idx_base = 0;
          uintptr_t pheno_cc_nm_bits = pheno_cc_nm[0];
          for (uint32_t uii = 0; uii != nm_case_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(pheno_cc_nm, &sample_idx_base, &pheno_cc_nm_bits);
            a1_case_dosage += S_CAST(double, geno_col[sample_idx]);
          }
          a1_case_dosages[allele_idx] = a1_case_dosage;
          block_aux_iter[nonomitted_allele_idx].sample_obs_ct = nm_sample_ct;
          block_aux_iter[nonomitted_allele_idx].allele_obs_ct = allele_obs_ct;
          if (!allele_ct_m2) {
            // Need main_dosage_sum and main_dosage_ssq for now (probably move
            // this computation in-place later).
            if (is_xchr_model_1) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = 0.0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const double cur_dosage = S_CAST(double, geno_col[sample_idx]);
                main_dosage_ssq += cur_dosage * cur_dosage;
              }
            } else {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = u63tod(machr2_dosage_ssqs[allele_idx]) * kRecipDosageMidSq;
              if (is_nonx_haploid) {
                main_dosage_ssq *= 0.25;
              }
            }
          }
          block_aux_iter[nonomitted_allele_idx].a1_dosage = a1_dosage;

          // bugfix (4 Sep 2018): forgot to save this
          block_aux_iter[nonomitted_allele_idx].case_allele_obs_ct = case_allele_obs_ct;

          block_aux_iter[nonomitted_allele_idx].a1_case_dosage = a1_case_dosage;
          block_aux_iter[nonomitted_allele_idx].firth_fallback = 0;
          block_aux_iter[nonomitted_allele_idx].is_unfinished = 0;
          block_aux_iter[nonomitted_allele_idx].mach_r2 = mach_r2;
          ++nonomitted_allele_idx;
        }
        // Now free to skip the actual regression if there are too few samples,
        // or omitted allele corresponds to a zero-variance genotype column.
        // If another allele has zero variance but the omitted allele does not,
        // we now salvage as many alleles as we can.
        GlmErr glm_err = 0;
        if (nm_sample_ct <= expected_predictor_ct) {
          // reasonable for this to override CONST_ALLELE
          glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
        } else if (IsSet(const_alleles, omitted_allele_idx)) {
          glm_err = SetGlmErr0(kGlmErrcodeConstOmittedAllele);
        }
        if (glm_err) {
          if (missing_ct) {
            // covariates have not been copied yet, so we can't usually change
            // prev_nm from 0 to 1 when missing_ct == 0 (and there's little
            // reason to optimize the zero-covariate case)
            prev_nm = 0;
          }
          uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
          if (allele_ct_m2 && (beta_se_multiallelic_fused || (!hide_covar))) {
            reported_ct += allele_ct_m2;
          }
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            for (uint32_t uii = 0; uii != reported_ct; ++uii) {
              memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
              beta_se_iter[uii * 2 + 1] = -9.0;
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        } else {
          {
            double omitted_dosage = u63tod(allele_obs_ct);
            double omitted_case_dosage = u63tod(case_allele_obs_ct);
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if (allele_idx == omitted_allele_idx) {
                continue;
              }
              omitted_dosage -= a1_dosages[allele_idx];
              omitted_case_dosage -= a1_case_dosages[allele_idx];
            }
            a1_dosages[omitted_allele_idx] = omitted_dosage;
            a1_case_dosages[omitted_allele_idx] = omitted_case_dosage;
          }
          uint32_t parameter_uidx = 2 + domdev_present;
          float* nm_predictors_pmaj_istart = nullptr;
          // only need to do this part once per variant in multiallelic case
          float* nm_predictors_pmaj_iter = &(nm_predictors_pmaj_buf[nm_sample_ctav * (parameter_uidx - main_omitted)]);
          if (missing_ct || (!prev_nm)) {
            // fill phenotype
            uintptr_t sample_midx_base = 0;
            uintptr_t sample_nm_bits = sample_nm[0];
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
              nm_pheno_buf[sample_idx] = cur_pheno[sample_midx];
            }
            // bugfix (13 Oct 2017): must guarantee trailing phenotype values
            // are valid (exact contents don't matter since they are multiplied
            // by zero, but they can't be nan)
            ZeroFArr(nm_sample_ct_rem, &(nm_pheno_buf[nm_sample_ct]));

            // fill covariates
            for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx, ++parameter_uidx) {
              // strictly speaking, we don't need cur_covars_cmaj to be
              // vector-aligned
              if (cur_parameter_subset && (!IsSet(cur_parameter_subset, parameter_uidx))) {
                continue;
              }
              const float* cur_covar_col;
              if (covar_idx < local_covar_ct) {
                cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
              } else {
                cur_covar_col = &(cur_covars_cmaj[(covar_idx - local_covar_ct) * sample_ctav]);
              }
              sample_midx_base = 0;
              sample_nm_bits = sample_nm[0];
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                *nm_predictors_pmaj_iter++ = cur_covar_col[sample_midx];
              }
              ZeromovFArr(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
            }
            nm_predictors_pmaj_istart = nm_predictors_pmaj_iter;
            // bugfix (13 Apr 2021): if local covariates are present, we can't
            // optimize as aggressively
            prev_nm = !(missing_ct || local_covar_ct);
          } else {
            // bugfix (15 Aug 2018): this was not handling --parameters
            // correctly when a covariate was only needed as part of an
            // interaction
            parameter_uidx += cur_covar_ct;
            nm_predictors_pmaj_istart = &(nm_predictors_pmaj_iter[literal_covar_ct * nm_sample_ctav]);
          }
          const uint32_t const_allele_ct = PopcountWords(const_alleles, allele_ctl);
          if (const_allele_ct) {
            // Must delete constant-allele columns from nm_predictors_pmaj, and
            // shift later columns back.
            float* read_iter = genotype_vals;
            float* write_iter = genotype_vals;
            for (uint32_t read_allele_idx = 0; read_allele_idx != allele_ct; ++read_allele_idx) {
              if (read_allele_idx == omitted_allele_idx) {
                continue;
              }
              if (!IsSet(const_alleles, read_allele_idx)) {
                if (write_iter != read_iter) {
                  memcpy(write_iter, read_iter, nm_sample_ctav * sizeof(float));
                }
                if (write_iter == genotype_vals) {
                  write_iter = multi_start;
                } else {
                  write_iter = &(write_iter[nm_sample_ctav]);
                }
              }
              if (read_iter == genotype_vals) {
                read_iter = multi_start;
              } else {
                read_iter = &(read_iter[nm_sample_ctav]);
              }
            }
          }
          const uint32_t cur_predictor_ct = expected_predictor_ct - const_allele_ct;
          const uint32_t cur_predictor_ctav = RoundUpPow2(cur_predictor_ct, kFloatPerFVec);
          const uint32_t cur_predictor_ctavp1 = cur_predictor_ctav + 1;
          uint32_t nonconst_extra_regression_idx = UINT32_MAX;  // deliberate overflow
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            float* main_vals = &(nm_predictors_pmaj_buf[nm_sample_ctav]);
            float* domdev_vals = nullptr;
            uint32_t is_unfinished = 0;
            uint32_t is_residualized = 0;
            // _stop instead of _ct since, in the residualized case, the
            // intercept (predictor index 0) is not included; we iterate over
            // the predictor indices in [1, _stop).
            uint32_t cur_regressed_predictor_stop = cur_predictor_ct;
            uint32_t cur_regressed_predictor_ctav = cur_predictor_ctav;
            uint32_t cur_regressed_predictor_ctavp1 = cur_predictor_ctavp1;
            uint32_t cur_biallelic_regressed_predictor_stop = cur_biallelic_predictor_ct;
            if (extra_regression_ct) {
              if (IsSet(const_alleles, extra_regression_idx + (extra_regression_idx >= omitted_allele_idx))) {
                glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
                goto GlmLogisticThreadF_skip_regression;
              }
              ++nonconst_extra_regression_idx;
              if (nonconst_extra_regression_idx) {
                float* swap_target = &(multi_start[(nonconst_extra_regression_idx - 1) * nm_sample_ctav]);
                for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
                  float fxx = genotype_vals[uii];
                  genotype_vals[uii] = swap_target[uii];
                  swap_target[uii] = fxx;
                }
              }
            }
            if (main_omitted) {
              // if main_mutated, this will be filled below
              // if not, this aliases genotype_vals
              main_vals = &(nm_predictors_pmaj_buf[(cur_predictor_ct + main_mutated) * nm_sample_ctav]);
            } else if (joint_genotypic || joint_hethom) {
              // in hethom case, do this before clobbering genotype data
              domdev_vals = &(main_vals[nm_sample_ctav]);
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                float cur_genotype_val = genotype_vals[sample_idx];
                if (cur_genotype_val > S_CAST(float, 1.0)) {
                  cur_genotype_val = S_CAST(float, 2.0) - cur_genotype_val;
                }
                domdev_vals[sample_idx] = cur_genotype_val;
              }
              ZeroFArr(nm_sample_ct_rem, &(domdev_vals[nm_sample_ct]));
            }
            if (model_dominant) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                float cur_genotype_val = genotype_vals[sample_idx];
                // 0..1..1
                if (cur_genotype_val > S_CAST(float, 1.0)) {
                  cur_genotype_val = 1.0;
                }
                main_vals[sample_idx] = cur_genotype_val;
              }
            } else if (model_recessive || joint_hethom) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                float cur_genotype_val = genotype_vals[sample_idx];
                // 0..0..1
                if (cur_genotype_val < S_CAST(float, 1.0)) {
                  cur_genotype_val = 0.0;
                } else {
                  cur_genotype_val -= S_CAST(float, 1.0);
                }
                main_vals[sample_idx] = cur_genotype_val;
              }
            } else if (model_hetonly) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                float cur_genotype_val = genotype_vals[sample_idx];
                // 0..1..0
                if (cur_genotype_val > S_CAST(float, 1.0)) {
                  cur_genotype_val = S_CAST(float, 2.0) - cur_genotype_val;
                }
                main_vals[sample_idx] = cur_genotype_val;
              }
            }

            // fill interaction terms
            if (add_interactions) {
              nm_predictors_pmaj_iter = nm_predictors_pmaj_istart;
              for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
                const float* cur_covar_col;
                if (covar_idx < local_covar_ct) {
                  cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                } else {
                  cur_covar_col = &(cur_covars_cmaj[covar_idx * sample_ctav]);
                }
                if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                  uintptr_t sample_midx_base = 0;
                  uintptr_t sample_nm_bits = sample_nm[0];
                  for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                    const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                    *nm_predictors_pmaj_iter++ = main_vals[sample_idx] * cur_covar_col[sample_midx];
                  }
                  ZeromovFArr(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
                }
                ++parameter_uidx;
                if (domdev_present) {
                  if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                    uintptr_t sample_midx_base = 0;
                    uintptr_t sample_nm_bits = sample_nm[0];
                    for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                      const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                      *nm_predictors_pmaj_iter++ = domdev_vals[sample_idx] * cur_covar_col[sample_midx];
                    }
                    ZeromovFArr(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
                  }
                  ++parameter_uidx;
                }
              }
            }
            if (corr_inv && prev_nm && (!allele_ct_m2)) {
              uintptr_t start_pred_idx = 0;
              if (!(model_dominant || model_recessive || model_hetonly || joint_hethom)) {
                start_pred_idx = domdev_present + 2;
                semicomputed_biallelic_xtx[cur_predictor_ct] = main_dosage_sum;
                semicomputed_biallelic_xtx[cur_predictor_ct + 1] = main_dosage_ssq;
              }
              if (cur_predictor_ct > start_pred_idx) {
                ColMajorFvectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[nm_sample_ctav]), &(nm_predictors_pmaj_buf[start_pred_idx * nm_sample_ctav]), nm_sample_ct, nm_sample_ctav, cur_predictor_ct - start_pred_idx, &(predictor_dotprod_buf[start_pred_idx]));
                for (uint32_t uii = start_pred_idx; uii != cur_predictor_ct; ++uii) {
                  semicomputed_biallelic_xtx[cur_predictor_ct + uii] = S_CAST(double, predictor_dotprod_buf[uii]);
                }
              }
              if (domdev_present) {
                ColMajorFvectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[2 * nm_sample_ctav]), nm_predictors_pmaj_buf, nm_sample_ct, nm_sample_ctav, cur_predictor_ct, predictor_dotprod_buf);
                for (uint32_t uii = 0; uii != cur_predictor_ct; ++uii) {
                  semicomputed_biallelic_xtx[2 * cur_predictor_ct + uii] = S_CAST(double, predictor_dotprod_buf[uii]);
                }
                semicomputed_biallelic_xtx[cur_predictor_ct + 2] = semicomputed_biallelic_xtx[2 * cur_predictor_ct + 1];
              }
              glm_err = CheckMaxCorrAndVifNm(semicomputed_biallelic_xtx, corr_inv, cur_predictor_ct, domdev_present_p1, cur_sample_ct_recip, cur_sample_ct_m1_recip, max_corr, vif_thresh, semicomputed_biallelic_corr_matrix, semicomputed_biallelic_inv_corr_sqrts, dbl_2d_buf, &(dbl_2d_buf[2 * cur_predictor_ct]), &(dbl_2d_buf[3 * cur_predictor_ct]));
              if (glm_err) {
                goto GlmLogisticThreadF_skip_regression;
              }
            } else {
              glm_err = CheckMaxCorrAndVifF(&(nm_predictors_pmaj_buf[nm_sample_ctav]), cur_predictor_ct - 1, nm_sample_ct, nm_sample_ctav, max_corr, vif_thresh, predictor_dotprod_buf, dbl_2d_buf, inverse_corr_buf, inv_1d_buf);
              if (glm_err) {
                goto GlmLogisticThreadF_skip_regression;
              }
            }
            ZeroFArr(cur_predictor_ctav, coef_return);
            if (!cur_is_always_firth) {
              // Does any genotype column have zero case or zero control
              // dosage?  If yes, faster to skip logistic regression than
              // wait for convergence failure.
              for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                if (IsSet(const_alleles, allele_idx)) {
                  continue;
                }
                const double tot_dosage = a1_dosages[allele_idx];
                const double case_dosage = a1_case_dosages[allele_idx];
                if ((case_dosage == 0.0) || (case_dosage == tot_dosage)) {
                  if (is_sometimes_firth) {
                    goto GlmLogisticThreadF_firth_fallback;
                  }
                  glm_err = SetGlmErr1(kGlmErrcodeSeparation, allele_idx);
                  goto GlmLogisticThreadF_skip_regression;
                }
              }
              if (!cur_cc_residualize) {
                if (LogisticRegressionF(nm_pheno_buf, nm_predictors_pmaj_buf, nullptr, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, cholesky_decomp_return, pp_buf, sample_variance_buf, hh_return, gradient_buf, dcoef_buf)) {
                  if (is_sometimes_firth) {
                    ZeroFArr(cur_predictor_ctav, coef_return);
                    goto GlmLogisticThreadF_firth_fallback;
                  }
                  glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
                  goto GlmLogisticThreadF_skip_regression;
                }
              } else {
                if (LogisticRegressionResidualizedF(nm_pheno_buf, nm_predictors_pmaj_buf, sample_nm, cur_cc_residualize, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, cholesky_decomp_return, pp_buf, sample_variance_buf, hh_return, gradient_buf, dcoef_buf, mean_centered_pmaj_buf, sample_offsets_buf)) {
                  if (is_sometimes_firth) {
                    ZeroFArr(cur_predictor_ctav, coef_return);
                    goto GlmLogisticThreadF_firth_fallback;
                  }
                  glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
                  goto GlmLogisticThreadF_skip_regression;
                }
                is_residualized = 1;
                cur_regressed_predictor_stop = domdev_present + allele_ct;
                cur_regressed_predictor_ctav = RoundUpPow2(cur_regressed_predictor_stop, kFloatPerFVec);
                cur_regressed_predictor_ctavp1 = cur_regressed_predictor_ctav + 1;
                cur_biallelic_regressed_predictor_stop = domdev_present + 2;
              }
              // unlike FirthRegressionF(), hh_return isn't inverted yet, do
              // that here
              for (uint32_t pred_uidx = is_residualized; pred_uidx != cur_regressed_predictor_stop; ++pred_uidx) {
                float* hh_inv_row = &(hh_return[pred_uidx * cur_regressed_predictor_ctav]);
                // ZeroFArr(cur_regressed_predictor_stop, gradient_buf);
                // gradient_buf[pred_uidx] = 1.0;
                // (y is gradient_buf, x is dcoef_buf)
                // SolveLinearSystemF(cholesky_decomp_return, &(gradient_buf[is_residualized]), cur_regressed_predictor_stop - is_residualized, &(hh_inv_row[is_residualized]));
                // that works, but doesn't exploit the sparsity of y

                // hh_return does now have vector-aligned rows
                ZeroFArr(pred_uidx, hh_inv_row);

                float fxx = 1.0;
                for (uint32_t row_idx = pred_uidx; row_idx != cur_regressed_predictor_stop; ++row_idx) {
                  const float* ll_row = &(cholesky_decomp_return[row_idx * cur_regressed_predictor_ctav]);
                  for (uint32_t col_idx = pred_uidx; col_idx != row_idx; ++col_idx) {
                    fxx -= ll_row[col_idx] * hh_inv_row[col_idx];
                  }
                  hh_inv_row[row_idx] = fxx / ll_row[row_idx];
                  fxx = 0.0;
                }
                for (uint32_t col_idx = cur_regressed_predictor_stop; col_idx > is_residualized; ) {
                  fxx = hh_inv_row[--col_idx];
                  float* hh_inv_row_iter = &(hh_inv_row[cur_regressed_predictor_stop - 1]);
                  for (uint32_t row_idx = cur_regressed_predictor_stop - 1; row_idx > col_idx; --row_idx) {
                    fxx -= cholesky_decomp_return[row_idx * cur_regressed_predictor_ctav + col_idx] * (*hh_inv_row_iter--);
                  }
                  *hh_inv_row_iter = fxx / cholesky_decomp_return[col_idx * cur_regressed_predictor_ctavp1];
                }
              }
            } else {
              if (!is_always_firth) {
              GlmLogisticThreadF_firth_fallback:
                block_aux_iter[extra_regression_idx].firth_fallback = 1;
                if (allele_ct_m2 && beta_se_multiallelic_fused) {
                  for (uint32_t uii = 1; uii != allele_ct - 1; ++uii) {
                    block_aux_iter[uii].firth_fallback = 1;
                  }
                }
              }
              if (!cur_cc_residualize) {
                if (FirthRegressionF(nm_pheno_buf, nm_predictors_pmaj_buf, nullptr, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, hh_return, inverse_corr_buf, inv_1d_buf, dbl_2d_buf, pp_buf, sample_variance_buf, gradient_buf, dcoef_buf, hdiag_buf, score_buf, hh0_buf, tmpnxk_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
                  goto GlmLogisticThreadF_skip_regression;
                }
              } else {
                if (FirthRegressionResidualizedF(nm_pheno_buf, nm_predictors_pmaj_buf, sample_nm, cur_cc_residualize, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, hh_return, inverse_corr_buf, inv_1d_buf, dbl_2d_buf, pp_buf, sample_variance_buf, gradient_buf, dcoef_buf, hdiag_buf, score_buf, hh0_buf, tmpnxk_buf, mean_centered_pmaj_buf, sample_offsets_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
                  goto GlmLogisticThreadF_skip_regression;
                }
                is_residualized = 1;
                cur_regressed_predictor_stop = domdev_present + allele_ct;
                cur_regressed_predictor_ctav = RoundUpPow2(cur_regressed_predictor_stop, kFloatPerFVec);
                cur_regressed_predictor_ctavp1 = cur_regressed_predictor_ctav + 1;
                cur_biallelic_regressed_predictor_stop = domdev_present + 2;
              }
            }
            // validParameters() check
            for (uint32_t pred_uidx = 1; pred_uidx != cur_regressed_predictor_stop; ++pred_uidx) {
              const float hh_inv_diag_element = hh_return[pred_uidx * cur_regressed_predictor_ctavp1];
              if ((hh_inv_diag_element < S_CAST(float, 1e-20)) || (!isfinite_f(hh_inv_diag_element))) {
                glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                goto GlmLogisticThreadF_skip_regression;
              }
              // use sample_variance_buf[] to store diagonal square roots
              sample_variance_buf[pred_uidx] = sqrtf(hh_inv_diag_element);
            }
            if (!is_residualized) {
              sample_variance_buf[0] = sqrtf(hh_return[0]);
            }
            for (uint32_t pred_uidx = 1 + is_residualized; pred_uidx != cur_regressed_predictor_stop; ++pred_uidx) {
              const float cur_hh_inv_diag_sqrt = S_CAST(float, 0.99999) * sample_variance_buf[pred_uidx];
              const float* hh_inv_row_iter = &(hh_return[pred_uidx * cur_regressed_predictor_ctav + is_residualized]);
              const float* hh_inv_diag_sqrts_iter = &(sample_variance_buf[is_residualized]);
              for (uint32_t pred_uidx2 = is_residualized; pred_uidx2 != pred_uidx; ++pred_uidx2) {
                if ((*hh_inv_row_iter++) > cur_hh_inv_diag_sqrt * (*hh_inv_diag_sqrts_iter++)) {
                  glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                  goto GlmLogisticThreadF_skip_regression;
                }
              }
            }
            if (is_unfinished) {
              block_aux_iter[extra_regression_idx].is_unfinished = 1;
              if (allele_ct_m2 && beta_se_multiallelic_fused) {
                for (uint32_t uii = 1; uii != allele_ct - 1; ++uii) {
                  block_aux_iter[uii].is_unfinished = 1;
                }
              }
            }
            {
              double* beta_se_iter2 = beta_se_iter;
              for (uint32_t pred_uidx = reported_pred_uidx_start; pred_uidx != reported_pred_uidx_biallelic_end; ++pred_uidx) {
                // In the multiallelic-fused case, if the first allele is
                // constant, this writes the beta/se values for the first
                // nonconstant, non-omitted allele where the results for the
                // first allele belong.  We correct that at the end of this
                // block.
                *beta_se_iter2++ = S_CAST(double, coef_return[pred_uidx]);
                *beta_se_iter2++ = S_CAST(double, sample_variance_buf[pred_uidx]);
              }
              if (cur_constraint_ct) {
                // bugfix (4 Sep 2021): forgot to update this for residualize
                // case
                *beta_se_iter2++ = 0.0;

                uint32_t joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_regressed_predictor_stop + joint_test_idx] = 1.0;
                }
                double chisq;
                if (!LinearHypothesisChisqF(coef_return, cur_constraints_con_major, hh_return, cur_constraint_ct, cur_regressed_predictor_stop, cur_regressed_predictor_ctav, &chisq, tmphxs_buf, h_transpose_buf, inner_buf, inverse_corr_buf, inv_1d_buf, dbl_2d_buf, outer_buf)) {
                  *beta_se_iter2++ = chisq;
                } else {
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeRankDeficient);
                  memcpy(&(beta_se_iter2[-1]), &glm_err2, 8);
                  *beta_se_iter2++ = -9.0;
                }
                // next test may have different alt allele count
                joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_regressed_predictor_stop + joint_test_idx] = 0.0;
                }
              }
              if (!const_allele_ct) {
                if (beta_se_multiallelic_fused || (!hide_covar)) {
                  for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
                    *beta_se_iter2++ = S_CAST(double, coef_return[cur_biallelic_regressed_predictor_stop + extra_allele_idx]);
                    *beta_se_iter2++ = S_CAST(double, sample_variance_buf[cur_biallelic_regressed_predictor_stop + extra_allele_idx]);
                  }
                }
              } else if (!beta_se_multiallelic_fused) {
                if (!hide_covar) {
                  // Need to insert some {CONST_ALLELE, -9} entries.
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                  const uint32_t cur_raw_allele_idx = extra_regression_idx + (extra_regression_idx >= omitted_allele_idx);
                  uint32_t extra_read_allele_idx = 0;
                  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    if ((allele_idx == omitted_allele_idx) || (allele_idx == cur_raw_allele_idx)) {
                      continue;
                    }
                    if (IsSet(const_alleles, allele_idx)) {
                      memcpy(beta_se_iter2, &glm_err2, 8);
                      beta_se_iter2[1] = -9.0;
                      beta_se_iter2 = &(beta_se_iter2[2]);
                    } else {
                      *beta_se_iter2++ = S_CAST(double, coef_return[cur_biallelic_regressed_predictor_stop + extra_read_allele_idx]);
                      *beta_se_iter2++ = S_CAST(double, sample_variance_buf[cur_biallelic_regressed_predictor_stop + extra_read_allele_idx]);
                      ++extra_read_allele_idx;
                    }
                  }
                }
              } else {
                const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                // Special-case first nonconst allele since it's positioned
                // discontinuously, and its BETA/SE may already be correctly
                // filled.
                uint32_t allele_idx = omitted_allele_idx? 0 : 1;
                if (IsSet(const_alleles, allele_idx)) {
                  memcpy(&(beta_se_iter[2 * include_intercept]), &glm_err2, 8);
                  beta_se_iter[2 * include_intercept + 1] = -9.0;
                  allele_idx = AdvTo0Bit(const_alleles, 1);
                  if (allele_idx == omitted_allele_idx) {
                    allele_idx = AdvTo0Bit(const_alleles, omitted_allele_idx + 1);
                  }
                  const uint32_t skip_ct = allele_idx - 1 - (allele_idx > omitted_allele_idx);
                  for (uint32_t uii = 0; uii != skip_ct; ++uii) {
                    memcpy(beta_se_iter2, &glm_err2, 8);
                    beta_se_iter2[1] = -9.0;
                    beta_se_iter2 = &(beta_se_iter2[2]);
                  }
                  *beta_se_iter2++ = S_CAST(double, coef_return[1]);
                  *beta_se_iter2++ = S_CAST(double, sample_variance_buf[1]);
                }
                ++allele_idx;
                uint32_t nonconst_allele_idx_m1 = 0;
                for (; allele_idx != allele_ct; ++allele_idx) {
                  if (allele_idx == omitted_allele_idx) {
                    continue;
                  }
                  if (!IsSet(const_alleles, allele_idx)) {
                    *beta_se_iter2++ = S_CAST(double, coef_return[cur_biallelic_predictor_ct + nonconst_allele_idx_m1]);
                    *beta_se_iter2++ = S_CAST(double, sample_variance_buf[cur_biallelic_predictor_ct + nonconst_allele_idx_m1]);
                    ++nonconst_allele_idx_m1;
                  } else {
                    memcpy(beta_se_iter2, &glm_err2, 8);
                    beta_se_iter2[1] = -9.0;
                    beta_se_iter2 = &(beta_se_iter2[2]);
                  }
                }
              }
            }
            while (0) {
            GlmLogisticThreadF_skip_regression:
              {
                uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
                if (allele_ct_m2 && (beta_se_multiallelic_fused || (!hide_covar))) {
                  reported_ct += allele_ct_m2;
                }
                for (uint32_t uii = 0; uii != reported_ct; ++uii) {
                  memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
                  beta_se_iter[uii * 2 + 1] = -9.0;
                }
              }
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        }
        block_aux_iter = &(block_aux_iter[allele_ct - 1]);
        if (local_covars_iter) {
          local_covars_iter = &(local_covars_iter[local_covar_ct * max_sample_ct]);
        }
      }
    }
    parity = 1 - parity;
    variant_idx_offset += cur_block_variant_ct;
    while (0) {
    GlmLogisticThreadF_err:
      UpdateU64IfSmaller(new_err_info, &common->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

// LogisticSseD doesn't exist since logistic_v_unsafe is used instead

static inline void MultMatrixDxnVectND(const double* mm, const double* vect, uint32_t col_ct, uint32_t row_ct, double* __restrict dest) {
  const uint32_t col_ctav = RoundUpPow2(col_ct, kDoublePerDVec);
  ColMajorVectorMatrixMultiplyStrided(vect, mm, col_ct, col_ctav, row_ct, dest);
}

#ifdef __LP64__
static inline void ComputeVD(const double* pp, uint32_t nn, double* __restrict vv) {
  const VecD one = VCONST_D(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kDoublePerDVec) {
    const VecD cur_pp = *R_CAST(const VecD*, &(pp[uii]));
    *R_CAST(VecD*, &(vv[uii])) = cur_pp * (one - cur_pp);
  }
}

static inline void ComputeVAndPMinusYD(const double* yy, uint32_t nn, double* pp, double* __restrict vv) {
  const VecD one = VCONST_D(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kDoublePerDVec) {
    VecD ptmp = *R_CAST(VecD*, &(pp[uii]));
    VecD one_minus_ptmp = one - ptmp;
    *R_CAST(VecD*, &(vv[uii])) = ptmp * one_minus_ptmp;
    VecD ytmp = *R_CAST(const VecD*, &(yy[uii]));
    *R_CAST(VecD*, &(pp[uii])) = ptmp - ytmp;
  }
}

static inline double TripleProductD(const double* v1, const double* v2, const double* v3, uint32_t nn) {
  VecD sum = vecd_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kDoublePerDVec) {
    VecD aa = *R_CAST(const VecD*, &(v1[uii]));
    VecD bb = *R_CAST(const VecD*, &(v2[uii]));
    VecD cc = *R_CAST(const VecD*, &(v3[uii]));
    sum = sum + aa * bb * cc;
  }
  return VecDHsum(sum);
}

static inline void ComputeTwoDiagTripleProductD(const double* aa, const double* bb, const double* vv, uint32_t nn, double* __restrict raa_ptr, double* __restrict rab_ptr, double* __restrict rbb_ptr) {
  VecD saa = vecd_setzero();
  VecD sab = vecd_setzero();
  VecD sbb = vecd_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kDoublePerDVec) {
    const VecD vtmp = *R_CAST(const VecD*, &(vv[uii]));
    const VecD atmp = *R_CAST(const VecD*, &(aa[uii]));
    const VecD btmp = *R_CAST(const VecD*, &(bb[uii]));
    const VecD av = atmp * vtmp;
    const VecD bv = btmp * vtmp;
    saa = saa + atmp * av;
    sab = sab + atmp * bv;
    sbb = sbb + btmp * bv;
  }
  *raa_ptr = VecDHsum(saa);
  *rab_ptr = VecDHsum(sab);
  *rbb_ptr = VecDHsum(sbb);
}

static inline void ComputeThreeTripleProductD(const double* bb, const double* a1, const double* a2, const double* a3, const double* vv, uint32_t nn, double* __restrict r1_ptr, double* __restrict r2_ptr, double* __restrict r3_ptr) {
  VecD s1 = vecd_setzero();
  VecD s2 = vecd_setzero();
  VecD s3 = vecd_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kDoublePerDVec) {
    const VecD a1tmp = *R_CAST(const VecD*, &(a1[uii]));
    const VecD a2tmp = *R_CAST(const VecD*, &(a2[uii]));
    const VecD a3tmp = *R_CAST(const VecD*, &(a3[uii]));
    const VecD vtmp = *R_CAST(const VecD*, &(vv[uii]));
    VecD btmp = *R_CAST(const VecD*, &(bb[uii]));
    btmp = btmp * vtmp;
    s1 = s1 + a1tmp * btmp;
    s2 = s2 + a2tmp * btmp;
    s3 = s3 + a3tmp * btmp;
  }
  *r1_ptr = VecDHsum(s1);
  *r2_ptr = VecDHsum(s2);
  *r3_ptr = VecDHsum(s3);
}

static inline void ComputeTwoPlusOneTripleProductD(const double* bb, const double* a1, const double* a2, const double* vv, uint32_t nn, double* __restrict r1_ptr, double* __restrict r2_ptr, double* __restrict r3_ptr) {
  VecD s1 = vecd_setzero();
  VecD s2 = vecd_setzero();
  VecD s3 = vecd_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kDoublePerDVec) {
    const VecD a1tmp = *R_CAST(const VecD*, &(a1[uii]));
    const VecD a2tmp = *R_CAST(const VecD*, &(a2[uii]));
    const VecD btmp = *R_CAST(const VecD*, &(bb[uii]));
    const VecD vtmp = *R_CAST(const VecD*, &(vv[uii]));
    const VecD bv = btmp * vtmp;
    s1 = s1 + btmp * bv;
    s2 = s2 + a1tmp * bv;
    s3 = s3 + a2tmp * bv;
  }
  *r1_ptr = VecDHsum(s1);
  *r2_ptr = VecDHsum(s2);
  *r3_ptr = VecDHsum(s3);
}
#else  // !__LP64__
static inline void ComputeVD(const double* pp, uint32_t nn, double* __restrict vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const double cur_pp = pp[uii];
    vv[uii] = cur_pp * (1.0 - cur_pp);
  }
}

static inline void ComputeVAndPMinusYD(const double* yy, uint32_t nn, double* __restrict pp, double* __restrict vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vv[uii] = pp[uii] * (1.0 - pp[uii]);
    pp[uii] -= yy[uii];
  }
}

static inline double TripleProductD(const double* v1, const double* v2, const double* v3, uint32_t nn) {
  double dxx = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    dxx += v1[uii] * v2[uii] * v3[uii];
  }
  return dxx;
}

static inline void ComputeTwoDiagTripleProductD(const double* aa, const double* bb, const double* vv, uint32_t nn, double* __restrict raa_ptr, double* __restrict rab_ptr, double* __restrict rbb_ptr) {
  double raa = 0.0;
  double rab = 0.0;
  double rbb = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const double dxx = aa[uii];
    const double dyy = bb[uii];
    double dzz = vv[uii];
    raa += dxx * dxx * dzz;
    dzz *= dyy;
    rab += dxx * dzz;
    rbb += dyy * dzz;
  }
  *raa_ptr = raa;
  *rab_ptr = rab;
  *rbb_ptr = rbb;
}

static inline void ComputeThreeTripleProductD(const double* bb, const double* a1, const double* a2, const double* a3, const double* vv, uint32_t nn, double* __restrict r1_ptr, double* __restrict r2_ptr, double* __restrict r3_ptr) {
  double r1 = 0.0;
  double r2 = 0.0;
  double r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const double dxx = bb[uii] * vv[uii];
    r1 += a1[uii] * dxx;
    r2 += a2[uii] * dxx;
    r3 += a3[uii] * dxx;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}

static inline void ComputeTwoPlusOneTripleProductD(const double* bb, const double* a1, const double* a2, const double* vv, uint32_t nn, double* __restrict r1_ptr, double* __restrict r2_ptr, double* __restrict r3_ptr) {
  double r1 = 0.0;
  double r2 = 0.0;
  double r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const double dxx = bb[uii];
    const double dyy = dxx * vv[uii];
    r1 += dxx * dyy;
    r2 += a1[uii] * dyy;
    r3 += a2[uii] * dyy;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}
#endif

double ComputeLoglikD(const double* yy, const double* pp, uint32_t sample_ct) {
  // possible todo: look for a high-precision way to accelerate this.
  double loglik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const double new_pi = pp[sample_idx];
    loglik += (yy[sample_idx] != 0.0)? log(new_pi) : log1p(-new_pi);
  }
  return loglik;
}

BoolErr ComputeLoglikCheckedD(const double* yy, const double* pp, uint32_t sample_ct, double* loglik_ptr) {
  // possible todo: look for a high-precision way to accelerate this.
  double loglik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const double new_pi = pp[sample_idx];
    // Imitate logistf stopping condition.
    if ((new_pi == 0.0) || (new_pi == 1.0)) {
      return 1;
    }
    loglik += (yy[sample_idx] != 0.0)? log(new_pi) : log1p(-new_pi);
  }
  *loglik_ptr = loglik;
  return 0;
}

void ComputeHessianD(const double* mm, const double* vv, uint32_t col_ct, uint32_t row_ct, double* __restrict dest) {
  const uintptr_t col_ctav = RoundUpPow2(col_ct, kDoublePerDVec);
  const uintptr_t row_ctav = RoundUpPow2(row_ct, kDoublePerDVec);
  const uintptr_t row_ctavp1 = row_ctav + 1;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    for (uint32_t row_idx = 0; row_idx < row_ctm3; row_idx += 3) {
      const double* mm_cur = &(mm[row_idx * col_ctav]);
      ComputeTwoDiagTripleProductD(mm_cur, &(mm_cur[col_ctav]), vv, col_ct, &(dest[row_idx * row_ctavp1]), &(dest[(row_idx + 1) * row_ctavp1 - 1]), &(dest[(row_idx + 1) * row_ctavp1]));
      ComputeTwoPlusOneTripleProductD(&(mm_cur[2 * col_ctav]), &(mm_cur[col_ctav]), mm_cur, vv, col_ct, &(dest[(row_idx + 2) * row_ctavp1]), &(dest[(row_idx + 2) * row_ctavp1 - 1]), &(dest[(row_idx + 2) * row_ctavp1 - 2]));
      for (uint32_t row_idx2 = row_idx + 3; row_idx2 != row_ct; ++row_idx2) {
        ComputeThreeTripleProductD(&(mm[row_idx2 * col_ctav]), mm_cur, &(mm_cur[col_ctav]), &(mm_cur[2 * col_ctav]), vv, col_ct, &(dest[row_idx2 * row_ctav + row_idx]), &(dest[row_idx2 * row_ctav + row_idx + 1]), &(dest[row_idx2 * row_ctav + row_idx + 2]));
      }
    }
  }
  switch (row_ct % 3) {
  case 0:
    ComputeTwoPlusOneTripleProductD(&(mm[(row_ct - 3) * col_ctav]), &(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 3) * row_ctavp1]), &(dest[(row_ct - 2) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1 - 2]));
    // fall through
  case 2:
    ComputeTwoDiagTripleProductD(&(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 2) * row_ctavp1]), &(dest[(row_ct - 1) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1]));
    break;
  case 1:
    dest[(row_ct - 1) * row_ctavp1] = TripleProductD(&(mm[(row_ct - 1) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct);
  }
}

void CholeskyDecompositionD(const double* aa, uint32_t predictor_ct, double* __restrict ll) {
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kDoublePerDVec);
  const uintptr_t predictor_ctavp1 = predictor_ctav + 1;
  const double* aa_diag_elem_ptr = aa;
  double* cur_ll_row = ll;
  for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
    double dxx = *aa_diag_elem_ptr;
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      const double dyy = cur_ll_row[col_idx];
      dxx -= dyy * dyy;
    }
    double dyy;
    if (dxx >= 0.0) {
      dyy = sqrt(dxx);
    } else {
      dyy = 1e-14;
    }
    cur_ll_row[row_idx] = dyy;
    dyy = 1.0 / dyy;  // now 1.0 / L[j][j]
    const double* aa_col_iter = aa_diag_elem_ptr;
    double* cur_ll_row2 = cur_ll_row;
    for (uint32_t row_idx2 = row_idx + 1; row_idx2 != predictor_ct; ++row_idx2) {
      aa_col_iter = &(aa_col_iter[predictor_ctav]);
      double dxx2 = *aa_col_iter;
      cur_ll_row2 = &(cur_ll_row2[predictor_ctav]);
      for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
        dxx2 -= cur_ll_row[col_idx] * cur_ll_row2[col_idx];
      }
      cur_ll_row2[row_idx] = dxx2 * dyy;
    }
    aa_diag_elem_ptr = &(aa_diag_elem_ptr[predictor_ctavp1]);
    cur_ll_row = &(cur_ll_row[predictor_ctav]);
  }
}

void SolveLinearSystemD(const double* ll, const double* yy, uint32_t predictor_ct, double* __restrict xx) {
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kDoublePerDVec);
  {
    const double* cur_ll_row = ll;
    for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
      double dxx = yy[row_idx];
      for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
        dxx -= cur_ll_row[col_idx] * xx[col_idx];
      }
      xx[row_idx] = dxx / cur_ll_row[row_idx];
      cur_ll_row = &(cur_ll_row[predictor_ctav]);
    }
  }
  for (uint32_t col_idx = predictor_ct; col_idx; ) {
    double* xx_stop = &(xx[--col_idx]);
    double dxx = *xx_stop;
    const double* ll_col_iter = &(ll[(predictor_ct - 1) * predictor_ctav + col_idx]);
    for (double* xx_iter = &(xx[predictor_ct - 1]); xx_iter != xx_stop; --xx_iter) {
      dxx -= (*ll_col_iter) * (*xx_iter);
      ll_col_iter -= predictor_ctav;
    }
    *xx_stop = dxx / (*ll_col_iter);
  }
}

BoolErr LogisticRegressionD(const double* yy, const double* xx, const double* sample_offsets, uint32_t sample_ct, uint32_t predictor_ct, uint32_t* is_unfinished_ptr, double* __restrict coef, double* __restrict ll, double* __restrict pp, double* __restrict vv, double* __restrict hh, double* __restrict grad, double* __restrict dcoef, MatrixInvertBuf1* __restrict mi_buf, double* __restrict dbl_2d_buf) {
  // This imitates R glm.fit().  Main differences from LogisticRegressionF(),
  // beyond precision:
  // - Initialization is somewhat different.
  // - Convergence criteria are now
  //     maxit = 25
  //     |dev - dev_{old}| / (|dev| + 0.1) < 1e-8, where dev := -2 *
  //       log-likelihood.
  // Support for sample_offsets was added in Mar 2024.  I had thought that
  // cc-residualize without single-prec didn't make sense since single-prec
  // should offer more of a speed boost for less accuracy cost on the initial
  // margin, but it turns out that the current single-precision functions have
  // unreliable convergence on biobank-scale datasets, so cc-residualize on its
  // own is often the fastest mode that still works.
  //
  // Preallocated buffers (initial contents irrelevant):
  // vv    = sample variance buffer
  // grad  = gradient buffer Y[] (length predictor_ct)
  // dcoef = current coefficient change buffer (length predictor_ct)
  // mi_buf, dbl_2d_buf: only for NOLAPACK case
  //
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         vector-aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype; trailing elements must be zeroed out
  //
  // Outputs:
  // coef  = main result.
  // ll    = cholesky decomposition matrix, predictor_ct^2, rows vector-aligned
  // hh    = hessian matrix buffer, predictor_ct^2, rows vector-aligned
  // pp    = final likelihoods minus Y[] (not currently used by callers).
  //
  // Returns 1 on convergence failure, 0 otherwise.
  // is_unfinished assumed to be initialized to 0, and is set to 1 if we hit
  // the iteration limit without satisfying other convergence criteria; the
  // main return value is 0 in this case.
  const uint32_t maxit = 25;

  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kDoublePerDVec);
  // some room to optimize this initialization
  double loglik_old;
  {
    // Initial loglik_old computation deliberately omitted, since it doesn't
    // correspond to a coef[] setting.  It's very unlikely for R glm.fit to
    // "converge" on the first iteration, but when it does, it doesn't have the
    // same meaning as satisfying the convergence criterion on later
    // iterations.
    double* zz = vv;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      // 2 * (ln 3 + 4/3)
      zz[sample_idx] = (yy[sample_idx] - 0.5) * 4.863891244002886;
    }

    // For the first iteration, we need to initialize coef in the same manner
    // as glm.fit() if we want to generate identical (up to Cholesky rounding)
    // results.
    BoolErr reterr = LinearRegressionDVec(zz, xx, predictor_ct, sample_ct, coef, hh, grad, mi_buf, dbl_2d_buf);
    if (unlikely(reterr)) {
      return 1;
    }

    // P[i] = \sum_j X[i][j] * coef[j];
    ColMajorMatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);
    if (sample_offsets) {
      AddDVec(sample_offsets, sample_ctav, pp);
    }

    // P[i] = 1 / (1 + exp(-P[i]));
    logistic_v_unsafe(pp, sample_ctav);

    const double loglik = ComputeLoglikD(yy, pp, sample_ct);
    if (loglik != loglik) {
      return 1;
    }
    loglik_old = loglik;
  }

  ZeroDArr(predictor_ct * predictor_ctav, ll);
  // bugfix (9 Oct 2024): although trailing elements of pp and vv are mostly
  // irrelevant, they cannot be nan
  ZeroDArr(sample_ctav - sample_ct, &(pp[sample_ct]));
  ZeroDArr(sample_ctav - sample_ct, &(vv[sample_ct]));

  // This index is 1 less than 'iter' in glm.R.
  for (uint32_t iteration = 1; iteration != maxit; ++iteration) {
    // V[i] = P[i] * (1 - P[i]);
    // P[i] -= Y[i];
    ComputeVAndPMinusYD(yy, sample_ctav, pp, vv);
    // V and P may both contain trailing garbage.

    ComputeHessianD(xx, vv, sample_ct, predictor_ct, hh);

    // grad = X^T P
    // Separate categorical loop also possible here
    ColMajorVectorMatrixMultiplyStrided(pp, xx, sample_ct, sample_ctav, predictor_ct, grad);

    // maybe this should use a QR decomposition instead?
    CholeskyDecompositionD(hh, predictor_ct, ll);

    SolveLinearSystemD(ll, grad, predictor_ct, dcoef);

    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      coef[pred_idx] -= dcoef[pred_idx];
    }
    // P[i] = \sum_j X[i][j] * coef[j];
    ColMajorMatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);
    if (sample_offsets) {
      AddDVec(sample_offsets, sample_ctav, pp);
    }

    // P[i] = 1 / (1 + exp(-P[i]));
    logistic_v_unsafe(pp, sample_ctav);
    const double loglik = ComputeLoglikD(yy, pp, sample_ct);
    if (loglik != loglik) {
      return 1;
    }

    // TODO: determine other non-convergence criteria
    if (fabs(loglik - loglik_old) < 1e-8 * (0.05 + fabs(loglik))) {
      return 0;
    }
    loglik_old = loglik;
  }
  *is_unfinished_ptr = 1;
  return 0;
}

#ifdef __LP64__
void CopyAndMeanCenterD(const double* src, uintptr_t ct, double* __restrict dst) {
  const uintptr_t fullvec_ct = ct / kDoublePerDVec;
  const VecD* src_alias = R_CAST(const VecD*, src);
  VecD vsum = vecd_setzero();
  for (uintptr_t vidx = 0; vidx != fullvec_ct; ++vidx) {
    vsum += src_alias[vidx];
  }
  double sum = VecDHsum(vsum);
  const uintptr_t trailing_start_idx = fullvec_ct * kDoublePerDVec;
  for (uintptr_t ulii = trailing_start_idx; ulii != ct; ++ulii) {
    sum += src[ulii];
  }

  const double neg_mean = -sum / u31tod(ct);
  const VecD neg_vmean = VCONST_D(neg_mean);
  VecD* dst_alias = R_CAST(VecD*, dst);
  for (uintptr_t vidx = 0; vidx != fullvec_ct; ++vidx) {
    dst_alias[vidx] = src_alias[vidx] + neg_vmean;
  }
  if (trailing_start_idx != ct) {
    for (uintptr_t ulii = trailing_start_idx; ulii != ct; ++ulii) {
      dst[ulii] = src[ulii] + neg_mean;
    }
    const uintptr_t trailing_stop_idx = trailing_start_idx + kDoublePerDVec;
    for (uintptr_t ulii = ct; ulii != trailing_stop_idx; ++ulii) {
      dst[ulii] = 0.0;
    }
  }
}
#else
void CopyAndMeanCenterD(const double* src, uintptr_t ct, double* __restrict dst) {
  double sum = 0.0;
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    sum += src[ulii];
  }
  const double mean = sum / u31tod(ct);
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    dst[ulii] = src[ulii] - mean;
  }
}
#endif

BoolErr LogisticRegressionResidualizedD(const double* yy, const double* xx, const uintptr_t* sample_nm, const CcResidualizeCtx* cc_residualize, uint32_t nm_sample_ct, uint32_t orig_predictor_ct, uint32_t* is_unfinished_ptr, double* coef, double* ll, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, double* pp, double* vv, double* hh, double* grad, double* dcoef, double* mean_centered_pmaj_buf, double* sample_offsets_buf) {
  if (!cc_residualize->logistic_nm_sample_offsets_d) {
    return 1;
  }
  const uintptr_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kDoublePerDVec);
  const uint32_t domdev_present_p1 = cc_residualize->domdev_present_p1;
  for (uint32_t geno_idx = 0; geno_idx != domdev_present_p1; ++geno_idx) {
    CopyAndMeanCenterD(&(xx[(geno_idx + 1) * nm_sample_ctav]), nm_sample_ct, &(mean_centered_pmaj_buf[geno_idx * nm_sample_ctav]));
  }
  const uint32_t prefitted_pred_ct = cc_residualize->prefitted_pred_ct;
  const uint32_t orig_biallelic_predictor_ct = domdev_present_p1 + prefitted_pred_ct;
  const uint32_t extra_allele_ct = orig_predictor_ct - orig_biallelic_predictor_ct;
  for (uint32_t extra_allele_idx = 0; extra_allele_idx != extra_allele_ct; ++extra_allele_idx) {
    CopyAndMeanCenterD(&(xx[(extra_allele_idx + orig_biallelic_predictor_ct) * nm_sample_ctav]), nm_sample_ct, &(mean_centered_pmaj_buf[(extra_allele_idx + domdev_present_p1) * nm_sample_ctav]));
  }
  const double* sample_offsets = cc_residualize->logistic_nm_sample_offsets_d;
  if (nm_sample_ct != cc_residualize->sample_ct) {
    uintptr_t sample_idx_base = 0;
    uintptr_t sample_nm_bits = sample_nm[0];
    for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(sample_nm, &sample_idx_base, &sample_nm_bits);
      sample_offsets_buf[uii] = sample_offsets[sample_idx];
    }
    // todo: check if this is actually needed
    const uint32_t remainder = (-nm_sample_ct) & (kDoublePerDVec - 1);
    ZeroDArr(remainder, &(sample_offsets_buf[nm_sample_ct]));
    sample_offsets = sample_offsets_buf;
  }
  // genotype, domdev?, other alleles
  const uint32_t regressed_predictor_ct = domdev_present_p1 + extra_allele_ct;
  const uint32_t regressed_predictor_ctav = RoundUpPow2(regressed_predictor_ct, kDoublePerDVec);
  if (LogisticRegressionD(yy, mean_centered_pmaj_buf, sample_offsets, nm_sample_ct, regressed_predictor_ct, is_unfinished_ptr, &(coef[1]), ll, pp, vv, hh, grad, dcoef, inv_1d_buf, dbl_2d_buf)) {
    return 1;
  }
  // hh and ll are shifted up and to the left from what the caller expects, due
  // to the missing intercept.  Correct that here.
  const uint32_t expected_predictor_ctav = RoundUpPow2(regressed_predictor_ct + 1, kDoublePerDVec);
  for (uint32_t write_row_idx = regressed_predictor_ct; write_row_idx; --write_row_idx) {
    memcpy(&(hh[write_row_idx * expected_predictor_ctav + 1]), &(hh[(write_row_idx - 1) * regressed_predictor_ctav]), regressed_predictor_ct * sizeof(double));
  }
  for (uint32_t write_row_idx = regressed_predictor_ct; write_row_idx; --write_row_idx) {
    memcpy(&(ll[write_row_idx * expected_predictor_ctav + 1]), &(ll[(write_row_idx - 1) * regressed_predictor_ctav]), regressed_predictor_ct * sizeof(double));
  }
  return 0;
}

void FirthComputeHdiagWeightsD(const double* yy, const double* xx, const double* pp, const double* hh, const double* vv, uint32_t predictor_ct, uint32_t predictor_ctav, uint32_t sample_ct, uint32_t sample_ctav, double* hdiag, double* ww, double* tmpnxk_buf) {
  ColMajorMatrixMultiplyStrided(xx, hh, sample_ct, sample_ctav, predictor_ct, predictor_ctav, predictor_ct, sample_ctav, tmpnxk_buf);
  // Assumes tau = 0.5.
#ifdef __LP64__
  // * tmpNxK, interpreted as column-major, is sample_ct x predictor_ct, must
  //   have trailing elements of each column zeroed out (otherwise nan may leak
  //   into trailing elements of hdiag, then ww, ...)
  // * X, interpreted as column-major, is also sample_ct x predictor_ct
  // * Hdiag[i] = V[i] (\sum_j tmpNxK[i][j] X[i][j])
  // * W[i] = (Y[i] - P[i]) + Hdiag[i] * (0.5 - P[i])
  const VecD half = VCONST_D(0.5);
  for (uint32_t sample_offset = 0; sample_offset < sample_ctav; sample_offset += kDoublePerDVec) {
    VecD dotprods = vecd_setzero();
    const double* xx_row = &(xx[sample_offset]);
    const double* tmpnxk_row = &(tmpnxk_buf[sample_offset]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const VecD cur_xx = *R_CAST(const VecD*, &(xx_row[pred_uidx * sample_ctav]));
      const VecD cur_tmpnxk = *R_CAST(const VecD*, &(tmpnxk_row[pred_uidx * sample_ctav]));
      dotprods = dotprods + cur_xx * cur_tmpnxk;
    }
    // Can handle categorical covariates in a separate loop here, and load the
    // dotprods increment into a union, etc.
    const VecD cur_vv = *R_CAST(const VecD*, &(vv[sample_offset]));
    const VecD cur_pi = *R_CAST(const VecD*, &(pp[sample_offset]));
    const VecD cur_yy = *R_CAST(const VecD*, &(yy[sample_offset]));
    const VecD cur_hdiag = cur_vv * dotprods;
    *R_CAST(VecD*, &(hdiag[sample_offset])) = cur_hdiag;
    const VecD half_minus_cur_pis = half - cur_pi;
    const VecD yy_minus_cur_pis = cur_yy - cur_pi;
    *R_CAST(VecD*, &(ww[sample_offset])) = yy_minus_cur_pis + cur_hdiag * half_minus_cur_pis;
  }
#else
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    double dotprod = 0.0;
    const double* xx_row = &(xx[sample_idx]);
    const double* tmpnxk_row = &(tmpnxk_buf[sample_idx]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      dotprod += xx_row[pred_uidx * sample_ctav] * tmpnxk_row[pred_uidx * sample_ctav];
    }
    const double cur_hdiag = vv[sample_idx] * dotprod;
    hdiag[sample_idx] = cur_hdiag;
    const double cur_pi = pp[sample_idx];
    ww[sample_idx] = (yy[sample_idx] - cur_pi) + cur_hdiag * (0.5 - cur_pi);
  }
#endif
}

void FirthComputeSecondWeightsD(const double* hdiag, const double* vv, __maybe_unused uint32_t sample_ct, __maybe_unused uint32_t sample_ctav, double* ww) {
#ifdef __LP64__
  const VecD one = VCONST_D(1.0);
  for (uint32_t sample_offset = 0; sample_offset < sample_ctav; sample_offset += kDoublePerDVec) {
    const VecD cur_hdiag = *R_CAST(const VecD*, &(hdiag[sample_offset]));
    const VecD cur_vv = *R_CAST(const VecD*, &(vv[sample_offset]));
    *R_CAST(VecD*, &(ww[sample_offset])) = (one + cur_hdiag) * cur_vv;
  }
#else
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    ww[sample_idx] = (1.0 + hdiag[sample_idx]) * vv[sample_idx];
  }
#endif
}

BoolErr FirthRegressionD(const double* yy, const double* xx, const double* sample_offsets, uint32_t sample_ct, uint32_t predictor_ct, double* beta, uint32_t* is_unfinished_ptr, double* hh, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, double* pp, double* vv, double* ustar, double* delta, double* hdiag, double* ww, double* hh0, double* tmpnxk_buf) {
  // This is a port of Georg Heinze's logistf 1.24.1 R function (called with
  // pl=FALSE), adapted to use many of plink 1.9's optimizations; see
  //   https://github.com/georgheinze/logistf
  //
  // Preallocated buffers (initial contents irrelevant):
  // inv_1d_buf, dbl_2d_buf = for matrix inversion
  // pp    = likelihoods minus Y[] (not currently used by callers)
  // vv    = sample variance buffer
  // ustar = gradient buffer (length predictor_ct)
  // delta = current coefficient change buffer (length predictor_ct)
  // hdiag = intermediate multiplier for weights (length sample_ct)
  // ww    = weight buffer (length sample_ct)
  // hh0   = just like hh below, but not returned
  //
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         vector-aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype
  // sample_offsets = in residualized mode, product of pre-fitted covariate
  //                  beta vector with covariate matrix.  otherwise, nullptr.
  //
  // Input/output:
  // beta  = starting point, overwritten with logistic regression betas.
  //
  // Outputs:
  // hh    = variance-covariance matrix buffer, predictor_ct^2, rows
  //         vector-aligned.  (spends some time as pre-inversion Hessian matrix
  //         too)
  //
  // Returns 1 on convergence failure, 0 otherwise.
  // is_unfinished assumed to be initialized to 0, and is set to 1 if we hit
  // the iteration limit without satisfying other convergence criteria; the
  // main return value is 0 in this case.
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kDoublePerDVec);
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);

  // ustar[] trailing elements must be zeroed out
  ZeroDArr(predictor_ctav - predictor_ct, &(ustar[predictor_ct]));

  const uint32_t trailing_sample_ct = sample_ctav - sample_ct;
  if (trailing_sample_ct) {
    // bunch of other trailing elements must also be zeroed out
    ZeroDArr(trailing_sample_ct, &(pp[sample_ct]));
    ZeroDArr(trailing_sample_ct, &(vv[sample_ct]));
    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      ZeroDArr(trailing_sample_ct, &(tmpnxk_buf[sample_ct + pred_idx * sample_ctav]));
    }
  }

  const uint32_t max_iter = 25;
  const double gconv = 0.00001;
  const double xconv = 0.00001;
  const double lconv = 0.00001;
  double delta_max = 0.0;
  double loglik_old = 0.0;
  for (uint32_t iter_idx = 0; ; ++iter_idx) {
    // P[i] = \sum_j beta[j] * X[i][j];
    // categorical optimization possible here
    ColMajorMatrixVectorMultiplyStrided(xx, beta, sample_ct, sample_ctav, predictor_ct, pp);
    if (sample_offsets) {
      AddDVec(sample_offsets, sample_ctav, pp);
    }
    // P[i] = 1 / (1 + exp(-P[i]));
    logistic_v_unsafe(pp, sample_ct);
    double loglik;
    if (ComputeLoglikCheckedD(yy, pp, sample_ct, &loglik)) {
      return 1;
    }
    // V[i] = P[i] * (1 - P[i]);
    ComputeVD(pp, sample_ct, vv);
    // P[i] -= Y[i] NOT done here

    // hessian = X diag(V) X'
    // note that only lower triangle is filled here
    ComputeHessianD(xx, vv, sample_ct, predictor_ct, hh0);
    // we shouldn't need to compute the log directly, since underflow <->
    // regression failure, right?  check this.
    if (InvertSymmdefMatrixFirstHalf(predictor_ct, predictor_ctav, hh0, inv_1d_buf, dbl_2d_buf)) {
      return 1;
    }
    const double dethh = HalfSymmInvertedDet(hh0, inv_1d_buf, predictor_ct, predictor_ctav);
    loglik += 0.5 * log(dethh);

    InvertSymmdefMatrixSecondHalf(predictor_ct, predictor_ctav, hh0, inv_1d_buf, dbl_2d_buf);
    // trailing elements of hh0[] rows can't be arbitrary for later
    // MultMatrixDxnVectND() call
    ReflectStridedMatrix0(predictor_ct, predictor_ctav, hh0);

    FirthComputeHdiagWeightsD(yy, xx, pp, hh0, vv, predictor_ct, predictor_ctav, sample_ct, sample_ctav, hdiag, ww, tmpnxk_buf);

    // gradient = X' W
    // categorical optimization possible here
    MultMatrixDxnVectND(xx, ww, sample_ct, predictor_ct, ustar);
    if (iter_idx) {
      double ustar_max = 0.0;
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        const double abs_ustar_cur = fabs(ustar[pred_uidx]);
        if (abs_ustar_cur > ustar_max) {
          ustar_max = abs_ustar_cur;
        }
      }
      const double loglik_change = loglik - loglik_old;
      if ((delta_max <= xconv) && (ustar_max < gconv) && (loglik_change < lconv)) {
        return 0;
      }
      if (iter_idx > max_iter) {
        *is_unfinished_ptr = 1;
        return 0;
      }
    }
    loglik_old = loglik;

    FirthComputeSecondWeightsD(hdiag, vv, sample_ct, sample_ctav, ww);
    ComputeHessianD(xx, ww, sample_ct, predictor_ct, hh);
    if (InvertSymmdefStridedMatrix(predictor_ct, predictor_ctav, hh, inv_1d_buf, dbl_2d_buf)) {
      return 1;
    }
    ReflectStridedMatrix0(predictor_ct, predictor_ctav, hh);

    // delta := hh * ustar (note that hh is inverted already)
    MultMatrixDxnVectND(hh, ustar, predictor_ct, predictor_ct, delta);

    delta_max = 0.0;
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const double abs_delta_cur = fabs(delta[pred_uidx]);
      if (abs_delta_cur > delta_max) {
        delta_max = abs_delta_cur;
      }
    }
    const double maxstep = 5.0;
    if (delta_max > maxstep) {
      const double scaling_factor = maxstep / delta_max;
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        delta[pred_uidx] *= scaling_factor;
      }
      delta_max = maxstep;
    }
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      beta[pred_uidx] += delta[pred_uidx];
    }
  }
}

BoolErr FirthRegressionResidualizedD(const double* yy, const double* xx, const uintptr_t* sample_nm, const CcResidualizeCtx* cc_residualize, uint32_t nm_sample_ct, uint32_t orig_predictor_ct, double* beta, uint32_t* is_unfinished_ptr, double* hh, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, double* pp, double* vv, double* ustar, double* delta, double* hdiag, double* ww, double* hh0_buf, double* tmpnxk_buf, double* mean_centered_pmaj_buf, double* sample_offsets_buf) {
  // todo: deduplicate with LogisticRegressionResidualizedD()
  const uintptr_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kDoublePerDVec);
  const uint32_t domdev_present_p1 = cc_residualize->domdev_present_p1;
  for (uint32_t geno_idx = 0; geno_idx != domdev_present_p1; ++geno_idx) {
    CopyAndMeanCenterD(&(xx[(geno_idx + 1) * nm_sample_ctav]), nm_sample_ct, &(mean_centered_pmaj_buf[geno_idx * nm_sample_ctav]));
  }
  const uint32_t prefitted_pred_ct = cc_residualize->prefitted_pred_ct;
  const uint32_t orig_biallelic_predictor_ct = domdev_present_p1 + prefitted_pred_ct;
  const uint32_t extra_allele_ct = orig_predictor_ct - orig_biallelic_predictor_ct;
  for (uint32_t extra_allele_idx = 0; extra_allele_idx != extra_allele_ct; ++extra_allele_idx) {
    CopyAndMeanCenterD(&(xx[(extra_allele_idx + orig_biallelic_predictor_ct) * nm_sample_ctav]), nm_sample_ct, &(mean_centered_pmaj_buf[(extra_allele_idx + domdev_present_p1) * nm_sample_ctav]));
  }
  const double* sample_offsets = cc_residualize->firth_nm_sample_offsets_d;
  if (nm_sample_ct != cc_residualize->sample_ct) {
    uintptr_t sample_idx_base = 0;
    uintptr_t sample_nm_bits = sample_nm[0];
    for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(sample_nm, &sample_idx_base, &sample_nm_bits);
      sample_offsets_buf[uii] = sample_offsets[sample_idx];
    }
    const uint32_t remainder = (-nm_sample_ct) & (kDoublePerDVec - 1);
    ZeroDArr(remainder, &(sample_offsets_buf[nm_sample_ct]));
    sample_offsets = sample_offsets_buf;
  }
  const uint32_t regressed_predictor_ct = domdev_present_p1 + extra_allele_ct;
  const uint32_t regressed_predictor_ctav = RoundUpPow2(regressed_predictor_ct, kDoublePerDVec);
  if (FirthRegressionD(yy, mean_centered_pmaj_buf, sample_offsets, nm_sample_ct, regressed_predictor_ct, &(beta[1]), is_unfinished_ptr, hh, inv_1d_buf, dbl_2d_buf, pp, vv, ustar, delta, hdiag, ww, hh0_buf, tmpnxk_buf)) {
    return 1;
  }
  // hh is shifted up and to the left from what the caller expects, due to the
  // missing intercept.  Correct that here.
  const uint32_t expected_predictor_ctav = RoundUpPow2(regressed_predictor_ct + 1, kDoublePerDVec);
  for (uint32_t write_row_idx = regressed_predictor_ct; write_row_idx; --write_row_idx) {
    memcpy(&(hh[write_row_idx * expected_predictor_ctav + 1]), &(hh[(write_row_idx - 1) * regressed_predictor_ctav]), regressed_predictor_ct * sizeof(double));
  }
  return 0;
}

BoolErr GlmAllocFillAndTestPhenoCovarsCc(const uintptr_t* sample_include, const uintptr_t* pheno_cc, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uint32_t domdev_present_p1, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, GlmFlags glm_flags, uintptr_t** pheno_cc_collapsed_ptr, uintptr_t** gcount_case_interleaved_vec_ptr, float** pheno_f_ptr, double** pheno_d_ptr, RegressionNmPrecomp** nm_precomp_ptr, float** covars_cmaj_f_ptr, double** covars_cmaj_d_ptr, CcResidualizeCtx** cc_residualize_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  const uint32_t is_single_prec = (glm_flags / kfGlmSinglePrecCc) & 1;
  const uintptr_t sample_ctav = is_single_prec? RoundUpPow2(sample_ct, kFloatPerFVec) : RoundUpPow2(sample_ct, kDoublePerDVec);
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
  const uint32_t is_cc_residualize = !!(glm_flags & (kfGlmFirthResidualize | kfGlmCcResidualize));
  if (unlikely(bigstack_alloc_w(sample_ctv * kWordsPerVec, pheno_cc_collapsed_ptr) ||
               bigstack_alloc_kcp(new_covar_ct, cur_covar_names_ptr))) {
      return 1;
  }
  if (is_single_prec) {
    if (unlikely(bigstack_alloc_f(sample_ctav, pheno_f_ptr) ||
                 bigstack_alloc_f(new_nonlocal_covar_ct * sample_ctav, covars_cmaj_f_ptr))) {
      return 1;
    }
  } else {
    if (unlikely(bigstack_alloc_d(sample_ctav, pheno_d_ptr))) {
      return 1;
    }
  }
  if (is_cc_residualize) {
    assert(!local_covar_ct);
    // Interactions and local covariates prohibited with firth-residualize,
    // so xtx_state guaranteed to be nonzero.
    if (unlikely(BIGSTACK_ALLOC_X(CcResidualizeCtx, 1, cc_residualize_ptr))) {
      return 1;
    }
    (*cc_residualize_ptr)->logistic_nm_sample_offsets_f = nullptr;
    (*cc_residualize_ptr)->firth_nm_sample_offsets_f = nullptr;
    (*cc_residualize_ptr)->logistic_nm_sample_offsets_d = nullptr;
    (*cc_residualize_ptr)->firth_nm_sample_offsets_d = nullptr;
    if ((glm_flags & (kfGlmFirth | kfGlmCcResidualize)) == kfGlmCcResidualize) {
      if (is_single_prec) {
        if (unlikely(bigstack_alloc_f(sample_ctav, &((*cc_residualize_ptr)->logistic_nm_sample_offsets_f)))) {
          return 1;
        }
      } else {
        if (unlikely(bigstack_alloc_d(sample_ctav, &((*cc_residualize_ptr)->logistic_nm_sample_offsets_d)))) {
          return 1;
        }
      }
    }
    if (!(glm_flags & kfGlmNoFirth)) {
      if (is_single_prec) {
        if (unlikely(bigstack_alloc_f(sample_ctav, &((*cc_residualize_ptr)->firth_nm_sample_offsets_f)))) {
          return 1;
        }
      } else {
        if (unlikely(bigstack_alloc_d(sample_ctav, &((*cc_residualize_ptr)->firth_nm_sample_offsets_d)))) {
          return 1;
        }
      }
    }
    (*cc_residualize_ptr)->prefitted_pred_ct = 1 + new_covar_ct;
    (*cc_residualize_ptr)->domdev_present_p1 = domdev_present_p1;
    (*cc_residualize_ptr)->sample_ct = sample_ct;
  }
  double* corr_buf = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  *nm_precomp_ptr = nullptr;
  if (xtx_state) {
    if (unlikely(BIGSTACK_ALLOC_X(RegressionNmPrecomp, 1, nm_precomp_ptr))) {
      return 1;
    }
    // x^2 + (2*xtx_state + 2)x + (5*xtx_state - 1)
    // 2x^2 + (2*xtx_state + 2)x + (5*xtx_state - 1)
    // 3x^2 + (4*xtx_state + 2)x + (8*xtx_state - 3)
    // 3x^2 + (4*xtx_state + 3)x + (8*xtx_state - 3)
    if (unlikely(bigstack_alloc_d(8 * xtx_state - 3 +
                                  new_covar_ct * (3 + 4 * xtx_state + 3 * new_covar_ct),
                                  &((*nm_precomp_ptr)->xtx_image)) ||
                 bigstack_alloc_d(new_covar_ct * (new_covar_ct + 1), &corr_buf))) {
      return 1;
    }
    (*nm_precomp_ptr)->covarx_dotprod_inv = nullptr;
    (*nm_precomp_ptr)->corr_inv = &((*nm_precomp_ptr)->xtx_image[(new_covar_ct + xtx_state + 1) * (new_covar_ct + xtx_state + 1)]);
    (*nm_precomp_ptr)->corr_image = &((*nm_precomp_ptr)->corr_inv[new_covar_ct * new_covar_ct]);
    (*nm_precomp_ptr)->corr_inv_sqrts = &((*nm_precomp_ptr)->corr_image[(new_covar_ct + xtx_state) * (new_covar_ct + xtx_state)]);
    (*nm_precomp_ptr)->xt_y_image = nullptr;
    bigstack_mark = R_CAST(unsigned char*, corr_buf);
  }

  double* covars_cmaj_d;
  double* covar_dotprod;
  double* inverse_corr_buf;
  if (unlikely(bigstack_alloc_d(new_nonlocal_covar_ct * sample_ct, &covars_cmaj_d) ||
               bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &covar_dotprod) ||
               bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &inverse_corr_buf))) {
    return 1;
  }
  if (!is_single_prec) {
    *covars_cmaj_d_ptr = covars_cmaj_d;
    bigstack_mark = R_CAST(unsigned char*, covar_dotprod);
  }
  uintptr_t* pheno_cc_collapsed = *pheno_cc_collapsed_ptr;
  CopyBitarrSubset(pheno_cc, sample_include, sample_ct, pheno_cc_collapsed);
  const uint32_t sample_remv = sample_ctav - sample_ct;
  if (is_single_prec) {
    float* pheno_f_iter = *pheno_f_ptr;
    for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      // can use the bitvector equivalent of GenoarrLookup...(), but this isn't
      // in a critical loop so I'll postpone writing those functions for now
      *pheno_f_iter++ = kSmallFloats[IsSet(pheno_cc_collapsed, sample_idx)];
    }
    ZeroFArr(sample_remv, pheno_f_iter);
  } else {
    double* pheno_d_iter = *pheno_d_ptr;
    for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      *pheno_d_iter++ = kSmallDoubles[IsSet(pheno_cc_collapsed, sample_idx)];
    }
    ZeroDArr(sample_remv, pheno_d_iter);
  }
  PglErr reterr = GlmFillAndTestCovars(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, max_corr, vif_thresh, covar_dotprod, corr_buf, inverse_corr_buf, covars_cmaj_d, *cur_covar_names_ptr, glm_err_ptr);
  if (unlikely(reterr)) {
    return (reterr == kPglRetNomem);
  }
  if (is_single_prec) {
    double* covar_read_iter = covars_cmaj_d;
    float* covar_write_iter = *covars_cmaj_f_ptr;
    for (uintptr_t covar_idx = 0; covar_idx != new_nonlocal_covar_ct; ++covar_idx) {
      for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        *covar_write_iter++ = S_CAST(float, *covar_read_iter++);
      }
      ZeroFArr(sample_remv, covar_write_iter);
      covar_write_iter = &(covar_write_iter[sample_remv]);
    }
  }
  if (xtx_state) {
    // error-out should be impossible
    InitNmPrecomp(covars_cmaj_d, covar_dotprod, corr_buf, inverse_corr_buf, sample_ct, 0, new_covar_ct, xtx_state, *nm_precomp_ptr);
    if (is_cc_residualize) {
      BigstackReset(bigstack_mark);
      const uintptr_t pred_ct = new_covar_ct + 1;
      const uintptr_t pred_ctav = RoundUpPow2(pred_ct, kFloatPerFVec);
      double* dbl_2d_buf;
      if (is_single_prec) {
        float* xx;
        float* coefs;
        float* hh;
        double* half_inverted_buf;
        float* ll;
        float* pp;
        float* vv;
        float* grad;
        float* dcoef;
        float* hdiag;
        float* ww;
        float* hh0_buf;
        float* tmpnxk_buf;
        MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, bigstack_alloc(pred_ct * kMatrixInvertBuf1CheckedAlloc));
        // see GetLogisticWorkspaceSizeF() and corresponding code at top of
        // GlmLogisticThread{F,D}()
        if (unlikely((!inv_1d_buf) ||
                     bigstack_alloc_f(sample_ctav * pred_ct, &xx) ||
                     bigstack_calloc_f(pred_ctav, &coefs) ||
                     bigstack_alloc_f(pred_ct * pred_ctav, &hh) ||
                     bigstack_alloc_d(pred_ct * MAXV(pred_ct, 3), &half_inverted_buf) ||
                     bigstack_alloc_d(pred_ct * MAXV(pred_ct, 7), &dbl_2d_buf) ||
                     bigstack_alloc_f(pred_ct * pred_ctav, &ll) ||
                     bigstack_alloc_f(sample_ctav, &pp) ||
                     bigstack_alloc_f(sample_ctav, &vv) ||
                     bigstack_alloc_f(pred_ctav, &grad) ||
                     bigstack_alloc_f(pred_ctav, &dcoef) ||
                     bigstack_alloc_f(sample_ctav, &hdiag) ||
                     bigstack_alloc_f(sample_ctav, &ww) ||
                     bigstack_alloc_f(pred_ct * pred_ctav, &hh0_buf) ||
                     bigstack_alloc_f(sample_ctav * pred_ct, &tmpnxk_buf))) {
          return 1;
        }
        FillFVec(sample_ct, 1.0, xx);
        memcpy(&(xx[sample_ctav]), *covars_cmaj_f_ptr, sample_ctav * (pred_ct - 1) * sizeof(float));
        uint32_t is_unfinished = 0;
        float* logistic_nm_sample_offsets_f = (*cc_residualize_ptr)->logistic_nm_sample_offsets_f;
        if (logistic_nm_sample_offsets_f) {
          if (unlikely(LogisticRegressionF(*pheno_f_ptr, xx, nullptr, sample_ct, pred_ct, coefs, &is_unfinished, ll, pp, vv, hh, grad, dcoef) || is_unfinished)) {
            if (glm_flags & kfGlmNoFirth) {
              *glm_err_ptr = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
              return 0;
            }
            (*cc_residualize_ptr)->logistic_nm_sample_offsets_f = nullptr;
          }
          ColMajorFmatrixVectorMultiplyStrided(xx, coefs, sample_ct, sample_ctav, pred_ct, logistic_nm_sample_offsets_f);
          ZeroFArr(sample_ctav - sample_ct, &(logistic_nm_sample_offsets_f[sample_ct]));
          ZeroFArr(pred_ctav, coefs);
        }
        float* firth_nm_sample_offsets_f = (*cc_residualize_ptr)->firth_nm_sample_offsets_f;
        if (firth_nm_sample_offsets_f) {
          if (unlikely(FirthRegressionF(*pheno_f_ptr, xx, nullptr, sample_ct, pred_ct, coefs, &is_unfinished, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf, pp, vv, grad, dcoef, hdiag, ww, hh0_buf, tmpnxk_buf) || is_unfinished)) {
            *glm_err_ptr = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
            return 0;
          }
          ColMajorFmatrixVectorMultiplyStrided(xx, coefs, sample_ct, sample_ctav, pred_ct, firth_nm_sample_offsets_f);
          ZeroFArr(sample_ctav - sample_ct, &(firth_nm_sample_offsets_f[sample_ct]));
        }
      } else {
        double* xx;
        double* coefs;
        double* hh;
        double* ll;
        double* pp;
        double* vv;
        double* grad;
        double* dcoef;
        double* hdiag;
        double* ww;
        double* hh0_buf;
        double* tmpnxk_buf;
        MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, bigstack_alloc(pred_ct * kMatrixInvertBuf1CheckedAlloc));
        if (unlikely((!inv_1d_buf) ||
                     bigstack_alloc_d(sample_ctav * pred_ct, &xx) ||
                     bigstack_alloc_d(pred_ctav, &coefs) ||
                     bigstack_alloc_d(pred_ct * pred_ctav, &hh) ||
                     bigstack_alloc_d(pred_ct * MAXV(pred_ct, 7), &dbl_2d_buf) ||
                     bigstack_alloc_d(pred_ct * pred_ctav, &ll) ||
                     bigstack_alloc_d(sample_ctav, &pp) ||
                     bigstack_alloc_d(sample_ctav, &vv) ||
                     bigstack_alloc_d(pred_ctav, &grad) ||
                     bigstack_alloc_d(pred_ctav, &dcoef) ||
                     bigstack_alloc_d(sample_ctav, &hdiag) ||
                     bigstack_alloc_d(sample_ctav, &ww) ||
                     bigstack_alloc_d(pred_ct * pred_ctav, &hh0_buf) ||
                     bigstack_alloc_d(sample_ctav * pred_ct, &tmpnxk_buf))) {
          return 1;
        }
        FillDVec(sample_ct, 1.0, xx);
        for (uintptr_t pred_idx = 1; pred_idx != pred_ct; ++pred_idx) {
          double* write_row_start = &(xx[sample_ctav * pred_idx]);
          memcpy(write_row_start, &(covars_cmaj_d[sample_ct * (pred_idx - 1)]), sample_ct * sizeof(double));
          ZeroDArr(sample_ctav - sample_ct, &(write_row_start[sample_ct]));
        }
        uint32_t is_unfinished = 0;
        double* logistic_nm_sample_offsets_d = (*cc_residualize_ptr)->logistic_nm_sample_offsets_d;
        if (logistic_nm_sample_offsets_d) {
          if (unlikely(LogisticRegressionD(*pheno_d_ptr, xx, nullptr, sample_ct, pred_ct, &is_unfinished, coefs, ll, pp, vv, hh, grad, dcoef, inv_1d_buf, dbl_2d_buf) || is_unfinished)) {
            if (glm_flags & kfGlmNoFirth) {
              *glm_err_ptr = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
              return 0;
            }
            (*cc_residualize_ptr)->logistic_nm_sample_offsets_d = nullptr;
          }
          ColMajorMatrixVectorMultiplyStrided(xx, coefs, sample_ct, sample_ctav, pred_ct, logistic_nm_sample_offsets_d);
          ZeroDArr(sample_ctav - sample_ct, &(logistic_nm_sample_offsets_d[sample_ct]));
        }
        double* firth_nm_sample_offsets_d = (*cc_residualize_ptr)->firth_nm_sample_offsets_d;
        if (firth_nm_sample_offsets_d) {
          ZeroDArr(pred_ctav, coefs);
          if (unlikely(FirthRegressionD(*pheno_d_ptr, xx, nullptr, sample_ct, pred_ct, coefs, &is_unfinished, hh, inv_1d_buf, dbl_2d_buf, pp, vv, grad, dcoef, hdiag, ww, hh0_buf, tmpnxk_buf) || is_unfinished)) {
            *glm_err_ptr = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
            return 0;
          }
          ColMajorMatrixVectorMultiplyStrided(xx, coefs, sample_ct, sample_ctav, pred_ct, firth_nm_sample_offsets_d);
          ZeroDArr(sample_ctav - sample_ct, &(firth_nm_sample_offsets_d[sample_ct]));
        }
      }
      // probable todo: print fitted coefficients and standard errors to log
    }
  }
  BigstackReset(bigstack_mark);
  if (gcount_case_interleaved_vec_ptr) {
    if (unlikely(bigstack_alloc_w(sample_ctv * kWordsPerVec, gcount_case_interleaved_vec_ptr))) {
      return 1;
    }
    ZeroTrailingWords(BitCtToWordCt(sample_ct), pheno_cc_collapsed);
    FillInterleavedMaskVec(pheno_cc_collapsed, sample_ctv, *gcount_case_interleaved_vec_ptr);
  }
  return 0;
}

uintptr_t GetLogisticWorkspaceSizeD(uint32_t sample_ct, uint32_t biallelic_predictor_ct, uint32_t domdev_present_p1, uint32_t max_extra_allele_ct, uint32_t constraint_ct, uint32_t xmain_ct, uint32_t gcount_cc, uint32_t is_sometimes_firth, uint32_t is_cc_residualize) {
  // sample_ctav * max_predictor_ct < 2^31, and sample_ct >=
  // biallelic_predictor_ct, so no overflows?
  // could round everything up to multiples of 16 instead of 64
  const uint32_t max_predictor_ct = biallelic_predictor_ct + max_extra_allele_ct;
  const uint32_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uint32_t max_predictor_ctav = RoundUpPow2(max_predictor_ct, kDoublePerDVec);
  // sample_nm, pheno_cc_nm, tmp_nm = sample_ctl words
  uintptr_t workspace_size = 3 * RoundUpPow2(BitCtToWordCt(sample_ct) * sizeof(intptr_t), kCacheline);

  // yy = sample_ctav doubles
  workspace_size += RoundUpPow2(sample_ctav * sizeof(double), kCacheline);

  // xx = (max_predictor_ct + main_mutated + main_omitted) * sample_ctav
  //      doubles
  workspace_size += RoundUpPow2((max_predictor_ct + xmain_ct) * sample_ctav * sizeof(double), kCacheline);

  // hh = max_predictor_ct * max_predictor_ctav doubles
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ctav * sizeof(double), kCacheline);

  // pp, vv = sample_ctav doubles
  workspace_size += 2 * RoundUpPow2(sample_ctav * sizeof(double), kCacheline);

  // coef, grad, dcoef = max_predictor_ctav doubles
  workspace_size += 3 * RoundUpPow2(max_predictor_ctav * sizeof(double), kCacheline);

  // ll = max_predictor_ct * max_predictor_ctav doubles
  // (technically not needed in pure-Firth case)
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ctav * sizeof(double), kCacheline);

  // semicomputed_biallelic_xtx
  workspace_size += RoundUpPow2(biallelic_predictor_ct * biallelic_predictor_ct * sizeof(double), kCacheline);

  // semicomputed_biallelic_corr_matrix
  workspace_size += RoundUpPow2((biallelic_predictor_ct - 1) * (biallelic_predictor_ct - 1) * sizeof(double), kCacheline);

  // semicomputed_biallelic_inv_corr_sqrts
  workspace_size += RoundUpPow2(biallelic_predictor_ct * sizeof(double), kCacheline);

  // inv_1d_buf
  workspace_size += RoundUpPow2(max_predictor_ct * kMatrixInvertBuf1CheckedAlloc, kCacheline);

  // dbl_2d_buf = max_predictor_ct * max_predictor_ct doubles, or VIF/Firth
  //              (which can be max_predictor_ct * 7 doubles)
  workspace_size += RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);

  // a1_dosages, a1_case_dosages
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(double) * 2, kCacheline);

  // machr2_dosage_sums, machr2_dosage_ssqs
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(uint64_t) * 2, kCacheline);

  if (gcount_cc && max_extra_allele_ct) {
    // case_one_cts, case_two_cts
    workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(int32_t) * 2, kCacheline);
  }

  // predictor_dotprod_buf
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ct * sizeof(double), kCacheline);

  const uintptr_t other_2d_byte_ct = max_predictor_ct * MAXV(max_predictor_ct, 3) * sizeof(double);
  // inverse_corr_buf/half_inverted_buf
  workspace_size += RoundUpPow2(other_2d_byte_ct, kCacheline);

  if (is_sometimes_firth) {
    // hdiag, ww = sample_ctav doubles
    workspace_size += 2 * RoundUpPow2(sample_ctav * sizeof(double), kCacheline);

    // hh0 = max_predictor_ct * max_predictor_ctav doubles
    workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ctav * sizeof(double), kCacheline);

    // tmpnxk_buf = max_predictor_ct * sample_ctav doubles
    workspace_size += RoundUpPow2(max_predictor_ct * sample_ctav * sizeof(double), kCacheline);
  }
  if (is_cc_residualize) {
    // mean_centered_pmaj_buf = (domdev_present_p1 + max_extra_allele_ct) *
    //   sample_ctav doubles
    workspace_size += RoundUpPow2((domdev_present_p1 + max_extra_allele_ct) * sample_ctav * sizeof(double), kCacheline);

    // sample_offsets_buf
    workspace_size += RoundUpPow2(sample_ctav * sizeof(double), kCacheline);
  }
  if (constraint_ct) {
    // tmphxs_buf, h_transpose_buf = constraint_ct * max_predictor_ctav doubles
    workspace_size += 2 * RoundUpPow2(constraint_ct * max_predictor_ctav * sizeof(double), kCacheline);

    // inner_buf = constraint_ct * constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * constraint_ct * sizeof(double), kCacheline);

    // outer_buf = constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * sizeof(double), kCacheline);

    // constraints_con_major = constraint_ct * max_predictor_ct
    workspace_size += RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(double), kCacheline);
  }
  return workspace_size;
}


THREAD_FUNC_DECL GlmLogisticThreadD(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GlmLogisticCtx* ctx = S_CAST(GlmLogisticCtx*, arg->sharedp->context);
  GlmCtx* common = ctx->common;

  PgenReader* pgrp = common->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = common->genovecs[tidx];
  pgv.dosage_present = nullptr;
  pgv.dosage_main = nullptr;
  if (common->dosage_presents) {
    pgv.dosage_present = common->dosage_presents[tidx];
    pgv.dosage_main = common->dosage_mains[tidx];
  }
  unsigned char* workspace_buf = common->workspace_bufs[tidx];
  const uintptr_t* variant_include = common->variant_include;
  const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
  const AlleleCode* omitted_alleles = common->omitted_alleles;
  const uintptr_t* sex_male_collapsed = common->sex_male_collapsed;
  const ChrInfo* cip = common->cip;
  const uint32_t* subset_chr_fo_vidx_start = common->subset_chr_fo_vidx_start;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const GlmFlags glm_flags = common->glm_flags;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t is_sometimes_firth = !(glm_flags & kfGlmNoFirth);
  const uint32_t is_always_firth = (glm_flags / kfGlmFirth) & 1;
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t model_hetonly = (glm_flags / kfGlmHetonly) & 1;
  const uint32_t joint_genotypic = (glm_flags / kfGlmGenotypic) & 1;
  const uint32_t joint_hethom = (glm_flags / kfGlmHethom) & 1;
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uint32_t domdev_present = joint_genotypic || joint_hethom;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t reported_pred_uidx_start = 1 - include_intercept;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const uintptr_t max_reported_test_ct = common->max_reported_test_ct;
  const uintptr_t local_covar_ct = common->local_covar_ct;
  const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
  // bugfix (20 Mar 2020): Also need to exclude dominant/recessive.
  const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!model_dominant) && (!model_recessive) && (!model_hetonly) && (!common->tests_flag) && (!add_interactions);
  uintptr_t max_sample_ct = MAXV(common->sample_ct, common->sample_ct_x);
  if (max_sample_ct < common->sample_ct_y) {
    max_sample_ct = common->sample_ct_y;
  }
  SetPgvThreadMhcNull(max_sample_ct, tidx, common->thread_mhc, &pgv);
  pgv.patch_01_ct = 0;
  pgv.patch_10_ct = 0;
  pgv.multidosage_sample_ct = 0;
  uint32_t variant_idx_offset = 0;
  uint32_t allele_ct = 2;
  uint32_t omitted_allele_idx = 0;
  uint32_t extra_regression_ct = 0;
  double main_dosage_sum = 0.0;
  double main_dosage_ssq = 0.0;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);

    double* beta_se_iter = common->block_beta_se;
    uintptr_t allele_bidx = variant_bidx;
    if (max_extra_allele_ct) {
      allele_bidx = variant_bidx + CountExtraAlleles(variant_include, allele_idx_offsets, common->read_variant_uidx_starts[0], common->read_variant_uidx_starts[tidx], 0);
    }
    if (beta_se_multiallelic_fused) {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * variant_bidx]);
    } else {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * allele_bidx]);
    }

    LogisticAuxResult* block_aux_iter = &(ctx->block_aux[allele_bidx]);
    const double* local_covars_iter = nullptr;
    if (local_covar_ct) {
      local_covars_iter = &(ctx->local_covars_vcmaj_d[parity][variant_bidx * max_sample_ct * local_covar_ct]);
    }
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = LastLeqU32(subset_chr_fo_vidx_start, 0, cip->chr_ct, variant_idx);
      // const uint32_t chr_fo_idx = LowerBoundNonemptyU32(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      uint32_t cur_variant_bidx_end = subset_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
      if (cur_variant_bidx_end > variant_bidx_end) {
        cur_variant_bidx_end = variant_bidx_end;
      }
      // "regular" = not all-female special case.
      const uint32_t is_haploid = IsSet(cip->haploid_mask, chr_idx);
      const uint32_t is_regular_x = is_haploid && (chr_idx == x_code);
      const uint32_t is_y = (chr_idx == y_code);
      const uint32_t is_nonx_haploid = is_haploid && (!is_regular_x);
      const uintptr_t* cur_sample_include;
      const uint32_t* cur_sample_include_cumulative_popcounts;
      const uintptr_t* cur_pheno_cc;
      const uintptr_t* cur_gcount_case_interleaved_vec;
      const double* cur_pheno;
      const RegressionNmPrecomp* nm_precomp;
      const double* cur_covars_cmaj;
      const uintptr_t* cur_parameter_subset;
      const uintptr_t* cur_joint_test_params;
      const CcResidualizeCtx* cur_cc_residualize;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      uint32_t cur_constraint_ct;
      uint32_t cur_is_always_firth;
      if (is_y && common->sample_include_y) {
        cur_sample_include = common->sample_include_y;
        cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
        cur_pheno_cc = ctx->pheno_y_cc;
        cur_gcount_case_interleaved_vec = ctx->gcount_case_interleaved_vec_y;
        cur_pheno = ctx->pheno_y_d;
        nm_precomp = common->nm_precomp_y;
        cur_covars_cmaj = ctx->covars_cmaj_y_d;
        cur_parameter_subset = common->parameter_subset_y;
        cur_joint_test_params = common->joint_test_params_y;
        cur_cc_residualize = ctx->cc_residualize_y;
        cur_sample_ct = common->sample_ct_y;
        cur_covar_ct = common->covar_ct_y;
        cur_constraint_ct = common->constraint_ct_y;
        cur_is_always_firth = is_always_firth || ctx->separation_found_y;
      } else if (is_regular_x && common->sample_include_x) {
        cur_sample_include = common->sample_include_x;
        cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
        cur_pheno_cc = ctx->pheno_x_cc;
        cur_gcount_case_interleaved_vec = ctx->gcount_case_interleaved_vec_x;
        cur_pheno = ctx->pheno_x_d;
        nm_precomp = common->nm_precomp_x;
        cur_covars_cmaj = ctx->covars_cmaj_x_d;
        cur_parameter_subset = common->parameter_subset_x;
        cur_joint_test_params = common->joint_test_params_x;
        cur_cc_residualize = ctx->cc_residualize_x;
        cur_sample_ct = common->sample_ct_x;
        cur_covar_ct = common->covar_ct_x;
        cur_constraint_ct = common->constraint_ct_x;
        cur_is_always_firth = is_always_firth || ctx->separation_found_x;
      } else {
        cur_sample_include = common->sample_include;
        cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
        cur_pheno_cc = ctx->pheno_cc;
        cur_gcount_case_interleaved_vec = ctx->gcount_case_interleaved_vec;
        cur_pheno = ctx->pheno_d;
        nm_precomp = common->nm_precomp;
        cur_covars_cmaj = ctx->covars_cmaj_d;
        cur_parameter_subset = common->parameter_subset;
        cur_joint_test_params = common->joint_test_params;
        cur_cc_residualize = ctx->cc_residualize;
        cur_sample_ct = common->sample_ct;
        cur_covar_ct = common->covar_ct;
        cur_constraint_ct = common->constraint_ct;
        cur_is_always_firth = is_always_firth || ctx->separation_found;
      }
      const uint32_t sample_ctl = BitCtToWordCt(cur_sample_ct);
      const uint32_t sample_ctav = RoundUpPow2(cur_sample_ct, kDoublePerDVec);
      const uint32_t cur_case_ct = PopcountWords(cur_pheno_cc, sample_ctl);
      const uint32_t cur_biallelic_predictor_ct_base = 2 + domdev_present + cur_covar_ct * (1 + add_interactions * domdev_present_p1);
      uint32_t cur_biallelic_predictor_ct = cur_biallelic_predictor_ct_base;
      uint32_t literal_covar_ct = cur_covar_ct;
      if (cur_parameter_subset) {
        cur_biallelic_predictor_ct = PopcountWords(cur_parameter_subset, BitCtToWordCt(cur_biallelic_predictor_ct_base));
        literal_covar_ct = PopcountBitRange(cur_parameter_subset, 2 + domdev_present, 2 + domdev_present + cur_covar_ct);
      }
      const uint32_t max_predictor_ct = cur_biallelic_predictor_ct + max_extra_allele_ct;
      const uint32_t max_predictor_ctav = RoundUpPow2(max_predictor_ct, kDoublePerDVec);
      uint32_t reported_pred_uidx_biallelic_end;
      if (hide_covar) {
        if (!cur_parameter_subset) {
          reported_pred_uidx_biallelic_end = 2 + domdev_present;
        } else {
          reported_pred_uidx_biallelic_end = 1 + IsSet(cur_parameter_subset, 1) + domdev_present;
        }
      } else {
        reported_pred_uidx_biallelic_end = cur_biallelic_predictor_ct;
      }
      // nm_predictors_pmaj_buf may require up to two extra columns omitted
      // from the main regression.
      // 1. In the multiallelic dominant/recessive/hetonly/hethom cases, the
      //    original genotype column does not appear in the regression, and
      //    we'd rather not reconstruct it from genovec, etc. when we need to
      //    swap it out for another allele, so we keep the original genotype in
      //    an extra column.
      //    To reduce code bloat, we now handle the biallelic cases in the same
      //    way; this is one of the more peripheral code paths so adding more
      //    complexity to speed it up is less justifiable.
      // 2. If --parameters excludes the main (possibly
      //    dominant/recessive/hetonly) genotype column but does care about an
      //    interaction, we want a copy of what the main genotype column's
      //    contents would have been to refer to.
      const uint32_t main_omitted = cur_parameter_subset && (!IsSet(cur_parameter_subset, 1));
      const uint32_t main_mutated = model_dominant || model_recessive || model_hetonly || joint_hethom;
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* pheno_cc_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* tmp_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      double* nm_pheno_buf = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), &workspace_iter));
      double* nm_predictors_pmaj_buf = S_CAST(double*, arena_alloc_raw_rd((max_predictor_ct + main_mutated + main_omitted) * sample_ctav * sizeof(double), &workspace_iter));
      double* coef_return = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ctav * sizeof(double), &workspace_iter));
      double* hh_return = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ctav * sizeof(double), &workspace_iter));
      double* pp_buf = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), &workspace_iter));
      double* sample_variance_buf = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), &workspace_iter));
      double* gradient_buf = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ctav * sizeof(double), &workspace_iter));
      double* dcoef_buf = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ctav * sizeof(double), &workspace_iter));
      double* cholesky_decomp_return = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ctav * sizeof(double), &workspace_iter));

      double* semicomputed_biallelic_xtx = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));
      // currently overallocates
      double* semicomputed_biallelic_corr_matrix = S_CAST(double*, arena_alloc_raw_rd((cur_biallelic_predictor_ct - 1) * (cur_biallelic_predictor_ct - 1) * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_inv_corr_sqrts = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));

      MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, arena_alloc_raw_rd(max_predictor_ct * kMatrixInvertBuf1CheckedAlloc, &workspace_iter));
      const uintptr_t dbl_2d_byte_ct = RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);
      double* dbl_2d_buf = S_CAST(double*, arena_alloc_raw(dbl_2d_byte_ct, &workspace_iter));
      double* a1_dosages = S_CAST(double*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(double) * 2, &workspace_iter));
      double* a1_case_dosages = &(a1_dosages[max_extra_allele_ct + 2]);
      uint64_t* machr2_dosage_sums = S_CAST(uint64_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(uint64_t) * 2, &workspace_iter));
      uint64_t* machr2_dosage_ssqs = &(machr2_dosage_sums[max_extra_allele_ct + 2]);
      uint32_t* case_one_cts = nullptr;
      uint32_t* case_two_cts = nullptr;
      if (cur_gcount_case_interleaved_vec && max_extra_allele_ct) {
        case_one_cts = S_CAST(uint32_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(int32_t) * 2, &workspace_iter));
        case_two_cts = &(case_one_cts[max_extra_allele_ct + 2]);
      }
      double* predictor_dotprod_buf = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ct * sizeof(double), &workspace_iter));
      const uintptr_t other_2d_byte_ct = RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 3) * sizeof(double), kCacheline);
      double* inverse_corr_buf = S_CAST(double*, arena_alloc_raw(other_2d_byte_ct, &workspace_iter));

      // these could use the same memory, but not a big deal, use the less
      // bug-prone approach for now
      // Firth-only
      double* hdiag_buf = nullptr;
      double* score_buf = nullptr;
      double* hh0_buf = nullptr;
      double* tmpnxk_buf = nullptr;
      if (is_sometimes_firth) {
        hdiag_buf = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), &workspace_iter));
        score_buf = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), &workspace_iter));
        hh0_buf = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ctav * sizeof(double), &workspace_iter));
        tmpnxk_buf = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * sample_ctav * sizeof(double), &workspace_iter));
      }
      double* mean_centered_pmaj_buf = nullptr;
      double* sample_offsets_buf = nullptr;
      if (cur_cc_residualize) {
        mean_centered_pmaj_buf = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double) * (domdev_present_p1 + max_extra_allele_ct), &workspace_iter));
        sample_offsets_buf = S_CAST(double*, arena_alloc_raw_rd(sample_ctav * sizeof(double), &workspace_iter));
      }

      // joint test only
      double* tmphxs_buf = nullptr;
      double* h_transpose_buf = nullptr;
      double* inner_buf = nullptr;
      double* outer_buf = nullptr;
      double* cur_constraints_con_major = nullptr;
      if (cur_constraint_ct) {
        tmphxs_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ctav * sizeof(double), &workspace_iter));
        h_transpose_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ctav * sizeof(double), &workspace_iter));
        inner_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * cur_constraint_ct * sizeof(double), &workspace_iter));
        outer_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * sizeof(double), &workspace_iter));
        cur_constraints_con_major = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        ZeroDArr(cur_constraint_ct * max_predictor_ct, cur_constraints_con_major);
        const uint32_t first_joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
        cur_constraints_con_major[first_joint_test_idx] = 1.0;
        // Rest of this matrix must be updated later, since cur_predictor_ct
        // changes at multiallelic variants.
      }
      assert(S_CAST(uintptr_t, workspace_iter - workspace_buf) == GetLogisticWorkspaceSizeD(cur_sample_ct, cur_biallelic_predictor_ct, domdev_present_p1, max_extra_allele_ct, cur_constraint_ct, main_mutated + main_omitted, cur_gcount_case_interleaved_vec != nullptr, is_sometimes_firth, cur_cc_residualize != nullptr));
      const double cur_sample_ct_recip = 1.0 / u31tod(cur_sample_ct);
      const double cur_sample_ct_m1_recip = 1.0 / u31tod(cur_sample_ct - 1);
      const double* corr_inv = nullptr;
      if (nm_precomp) {
        memcpy(semicomputed_biallelic_xtx, nm_precomp->xtx_image, cur_biallelic_predictor_ct * cur_biallelic_predictor_ct * sizeof(double));
        corr_inv = nm_precomp->corr_inv;
        const uintptr_t nongeno_pred_ct = cur_biallelic_predictor_ct - domdev_present - 2;
        const uintptr_t nonintercept_biallelic_pred_ct = cur_biallelic_predictor_ct - 1;
        memcpy(semicomputed_biallelic_corr_matrix, nm_precomp->corr_image, nonintercept_biallelic_pred_ct * nonintercept_biallelic_pred_ct * sizeof(double));
        memcpy(&(semicomputed_biallelic_inv_corr_sqrts[domdev_present_p1]), nm_precomp->corr_inv_sqrts, nongeno_pred_ct * sizeof(double));
      }
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      // when this is set, the last fully-processed variant had no missing
      // genotypes, and if the current variant also has no missing genotypes we
      // may be able to skip reinitialization of most of
      // nm_predictors_pmaj_buf.
      // (todo: do we want to track prev_biallelic_nm?)
      uint32_t prev_nm = 0;

      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        if (allele_idx_offsets) {
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
          if (!beta_se_multiallelic_fused) {
            extra_regression_ct = allele_ct - 2;
          }
        }
        const uint32_t allele_ct_m2 = allele_ct - 2;
        const uint32_t expected_predictor_ct = cur_biallelic_predictor_ct + allele_ct_m2;
        PglErr reterr;
        if (!allele_ct_m2) {
          reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
        } else {
          reterr = PgrGetMD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, &pgv);
          // todo: proper multiallelic dosage support
        }
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto GlmLogisticThreadD_err;
        }
        ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
        GenoarrCountFreqsUnsafe(pgv.genovec, cur_sample_ct, genocounts);
        uint32_t missing_ct = genocounts[3];
        if (!missing_ct) {
          SetAllBits(cur_sample_ct, sample_nm);
        } else {
          GenoarrToNonmissing(pgv.genovec, cur_sample_ct, sample_nm);
          if (pgv.dosage_ct) {
            BitvecOr(pgv.dosage_present, sample_ctl, sample_nm);
            missing_ct = cur_sample_ct - PopcountWords(sample_nm, sample_ctl);
          }
        }
        if (omitted_alleles) {
          omitted_allele_idx = omitted_alleles[variant_uidx];
        }
        // Once sizeof(AlleleCode) > 1, we probably want to allocate this from
        // g_bigstack instead of the thread stack.
        uintptr_t const_alleles[DivUp(kPglMaxAlleleCt, kBitsPerWord)];
        const uint32_t allele_ctl = DivUp(allele_ct, kBitsPerWord);
        ZeroWArr(allele_ctl, const_alleles);
        const uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
        const uint32_t nm_sample_ctl = BitCtToWordCt(nm_sample_ct);
        const uint32_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kDoublePerDVec);
        const uint32_t nm_sample_ct_rem = nm_sample_ctav - nm_sample_ct;
        // first predictor column: intercept
        if (!prev_nm) {
          FillDVec(nm_sample_ct, 1.0, nm_predictors_pmaj_buf);
        }
        // second predictor column: genotype
        double* genotype_vals = &(nm_predictors_pmaj_buf[nm_sample_ctav]);
        if (main_mutated || main_omitted) {
          // bugfix (9 Oct 2024)
          ZeroDArr(nm_sample_ct_rem, &(genotype_vals[nm_sample_ct]));
          genotype_vals = &(nm_predictors_pmaj_buf[expected_predictor_ct * nm_sample_ctav]);
        }
        CopyBitarrSubset(cur_pheno_cc, sample_nm, nm_sample_ct, pheno_cc_nm);
        const uint32_t nm_case_ct = PopcountWords(pheno_cc_nm, nm_sample_ctl);
        double* multi_start = nullptr;
        if (!allele_ct_m2) {
          if (omitted_allele_idx) {
            GenovecInvertUnsafe(cur_sample_ct, pgv.genovec);
            // ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
            }
            const uint32_t uii = genocounts[0];
            genocounts[0] = genocounts[2];
            genocounts[2] = uii;
          }
          uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
          uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
          if (!missing_ct) {
            GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, genotype_vals);
            if (pgv.dosage_ct) {
              uintptr_t sample_idx_base = 0;
              uintptr_t dosage_present_bits = pgv.dosage_present[0];
              for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                const uintptr_t sample_idx = BitIter1(pgv.dosage_present, &sample_idx_base, &dosage_present_bits);
                const uint32_t dosage_val = pgv.dosage_main[dosage_idx];
                // 32768 -> 2, 16384 -> 1, 0 -> 0
                genotype_vals[sample_idx] = kRecipDosageMid * u31tod(dosage_val);
                dosage_sum += dosage_val;
                dosage_ssq += dosage_val * dosage_val;
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_idx);
                if (cur_geno && (cur_geno != 3)) {
                  const uintptr_t prev_val = cur_geno * kDosageMid;
                  dosage_sum -= prev_val;
                  dosage_ssq -= prev_val * prev_val;
                }
              }
            }
          } else {
            if (!pgv.dosage_ct) {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, genotype_vals);
            } else {
              uintptr_t sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              uint32_t dosage_idx = 0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_midx);
                double cur_val;
                if (IsSet(pgv.dosage_present, sample_midx)) {
                  const uint32_t dosage_val = pgv.dosage_main[dosage_idx++];
                  cur_val = kRecipDosageMid * u31tod(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                } else {
                  // cur_geno != 3 guaranteed
                  cur_val = kSmallDoubles[cur_geno];
                }
                genotype_vals[sample_idx] = cur_val;
              }
            }
          }
          // Check for constant genotype column.
          // (Technically, we should recheck later in the chrX no-sex-covariate
          // --xchr-model 1 corner case.)
          if (!pgv.dosage_ct) {
            if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          } else if (pgv.dosage_ct == nm_sample_ct) {
            if (DosageIsConstant(dosage_sum, dosage_ssq, nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          }
          machr2_dosage_sums[1 - omitted_allele_idx] = dosage_sum;
          machr2_dosage_ssqs[1 - omitted_allele_idx] = dosage_ssq;
          machr2_dosage_sums[omitted_allele_idx] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[omitted_allele_idx] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
          if (cur_gcount_case_interleaved_vec) {
            // gcountcc
            STD_ARRAY_REF(uint32_t, 6) cur_geno_hardcall_cts = block_aux_iter->geno_hardcall_cts;
            GenoarrCountSubsetFreqs(pgv.genovec, cur_gcount_case_interleaved_vec, cur_sample_ct, cur_case_ct, R_CAST(STD_ARRAY_REF(uint32_t, 4), cur_geno_hardcall_cts));
            for (uint32_t geno_hardcall_idx = 0; geno_hardcall_idx != 3; ++geno_hardcall_idx) {
              cur_geno_hardcall_cts[3 + geno_hardcall_idx] = genocounts[geno_hardcall_idx] - cur_geno_hardcall_cts[geno_hardcall_idx];
            }
          }
        } else {
          // multiallelic.
          // Update (18 Mar 2020): If some but not all alleles have constant
          // dosages, we remove just those alleles from the regressions;
          // trim-alts is not necessary to see what's going on with the other
          // alleles.  To reduce parsing complexity, the number of output lines
          // is not affected by this; the ones corresponding to the constant
          // alleles have NA values.

          // dosage_ct == 0 temporarily guaranteed if we reach here.
          assert(!pgv.dosage_ct);
          multi_start = &(nm_predictors_pmaj_buf[(expected_predictor_ct - allele_ct_m2) * nm_sample_ctav]);
          ZeroU64Arr(allele_ct, machr2_dosage_sums);
          ZeroU64Arr(allele_ct, machr2_dosage_ssqs);
          // postpone multiply for now, since no multiallelic dosages
          // Use sums as ones[] and ssqs as twos[] for rarealts; transform to
          // actual sums/ssqs later.
          machr2_dosage_sums[0] = genocounts[1];
          machr2_dosage_ssqs[0] = genocounts[0];
          if (omitted_allele_idx) {
            // Main genotype column starts as REF.
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallInvDoublePairs, nm_sample_ct, genotype_vals);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallInvDoubles, cur_sample_ct, genotype_vals);
            }
          }
          uint32_t rare_allele_ct = allele_ct_m2;
          double* alt1_start = nullptr;
          double* rarealt_start = multi_start;
          if (omitted_allele_idx != 1) {
            if (omitted_allele_idx) {
              alt1_start = multi_start;
              ZeroDArr(nm_sample_ct_rem, &(alt1_start[nm_sample_ct]));
              rarealt_start = &(rarealt_start[nm_sample_ctav]);
              --rare_allele_ct;
            } else {
              alt1_start = genotype_vals;
            }
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, alt1_start);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, alt1_start);
            }
          }
          ZeroDArr(rare_allele_ct * nm_sample_ctav, rarealt_start);
          if (pgv.patch_01_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_01_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_01_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ctav + sample_idx] = 1.0;
                alt1_start[sample_idx] = 0.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ctav + sample_idx] = 1.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                alt1_start[sample_idx] = 0.0;
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                machr2_dosage_sums[allele_code] += 1;
                if (allele_code == omitted_allele_idx) {
                  continue;
                }
                const uint32_t cur_col = allele_code - 2 - (allele_code > omitted_allele_idx);
                rarealt_start[cur_col * nm_sample_ctav + sample_idx] = 1.0;
              }
            }
          }
          uintptr_t alt1_het_ct = genocounts[1] - pgv.patch_01_ct;
          if (pgv.patch_10_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_10_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_10_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] = 2.0;
                  alt1_start[sample_idx] = 0.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ctav + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] += 1.0;
                    alt1_start[sample_idx] = 0.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] = 2.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ctav + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] += 1.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  machr2_dosage_ssqs[ac0] += 1;
                  alt1_start[sample_idx] = 0.0;
                  if (ac0 != omitted_allele_idx) {
                    const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                    rarealt_start[ac0_col * nm_sample_ctav + sample_idx] = 2.0;
                  }
                } else {
                  machr2_dosage_sums[ac1] += 1;
                  if (ac1 != omitted_allele_idx) {
                    const uint32_t ac1_col = ac1 - 2 - (ac1 > omitted_allele_idx);
                    rarealt_start[ac1_col * nm_sample_ctav + sample_idx] = 1.0;
                  }
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    machr2_dosage_sums[ac0] += 1;
                    alt1_start[sample_idx] = 0.0;
                    if (ac0 != omitted_allele_idx) {
                      const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                      rarealt_start[ac0_col * nm_sample_ctav + sample_idx] += 1.0;
                    }
                  }
                }
              }
            }
          }
          machr2_dosage_sums[1] = alt1_het_ct;
          machr2_dosage_ssqs[1] = genocounts[2] - pgv.patch_10_ct;
          if (cur_gcount_case_interleaved_vec) {
            // gcountcc.  Need case-specific one_cts and two_cts for each
            // allele.
            STD_ARRAY_DECL(uint32_t, 4, case_hardcall_cts);
            GenoarrCountSubsetFreqs(pgv.genovec, cur_gcount_case_interleaved_vec, cur_sample_ct, cur_case_ct, case_hardcall_cts);
            ZeroU32Arr(allele_ct, case_one_cts);
            ZeroU32Arr(allele_ct, case_two_cts);
            uint32_t case_alt1_het_ct = case_hardcall_cts[1];
            case_one_cts[0] = case_alt1_het_ct;
            case_two_cts[0] = case_hardcall_cts[0];
            if (pgv.patch_01_ct) {
              uintptr_t sample_widx = 0;
              uintptr_t cur_bits = pgv.patch_01_set[0];
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t lowbit = BitIter1y(pgv.patch_01_set, &sample_widx, &cur_bits);
                if (cur_pheno_cc[sample_widx] & lowbit) {
                  const uint32_t allele_code = pgv.patch_01_vals[uii];
                  case_one_cts[allele_code] += 1;
                }
              }
              for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
                case_alt1_het_ct -= case_one_cts[allele_idx];
              }
            }
            uint32_t case_alt1_hom_ct = case_hardcall_cts[2];
            if (pgv.patch_10_ct) {
              uintptr_t sample_widx = 0;
              uintptr_t cur_bits = pgv.patch_10_set[0];
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t lowbit = BitIter1y(pgv.patch_10_set, &sample_widx, &cur_bits);
                if (cur_pheno_cc[sample_widx] & lowbit) {
                  const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                  const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                  --case_alt1_hom_ct;
                  if (ac0 == ac1) {
                    case_two_cts[ac0] += 1;
                  } else {
                    case_one_cts[ac1] += 1;
                    if (ac0 == 1) {
                      ++case_alt1_het_ct;
                    } else {
                      case_one_cts[ac0] += 1;
                    }
                  }
                }
              }
            }
            case_one_cts[1] = case_alt1_het_ct;
            case_two_cts[1] = case_alt1_hom_ct;
            uint32_t nonomitted_allele_idx = 0;
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if (allele_idx == omitted_allele_idx) {
                continue;
              }
              const uint32_t one_ct = machr2_dosage_sums[allele_idx];
              const uint32_t two_ct = machr2_dosage_ssqs[allele_idx];
              const uint32_t case_one_ct = case_one_cts[allele_idx];
              const uint32_t case_two_ct = case_two_cts[allele_idx];
              STD_ARRAY_REF(uint32_t, 6) dst = block_aux_iter[nonomitted_allele_idx].geno_hardcall_cts;
              dst[0] = nm_case_ct - case_one_ct - case_two_ct;
              dst[1] = case_one_ct;
              dst[2] = case_two_ct;
              dst[3] = nm_sample_ct - one_ct - two_ct - dst[0];
              dst[4] = one_ct - case_one_ct;
              dst[5] = two_ct - case_two_ct;
              ++nonomitted_allele_idx;
            }
          }
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            const uintptr_t one_ct = machr2_dosage_sums[allele_idx];
            const uintptr_t two_ct = machr2_dosage_ssqs[allele_idx];
            machr2_dosage_sums[allele_idx] = (one_ct + 2 * two_ct) * 0x4000LLU;
            machr2_dosage_ssqs[allele_idx] = (one_ct + 4LLU * two_ct) * 0x10000000LLU;
            if ((one_ct == nm_sample_ct) || (two_ct == nm_sample_ct) || ((!one_ct) && (!two_ct))) {
              SetBit(allele_idx, const_alleles);
            }
          }
        }
        ZeroDArr(nm_sample_ct_rem, &(genotype_vals[nm_sample_ct]));
        // usually need to save some of {sample_obs_ct, allele_obs_ct,
        // a1_dosage, case_allele_obs_ct, a1_case_dosage, mach_r2 even for
        // skipped variants
        // compute them all for now, could conditionally skip later
        uint32_t allele_obs_ct = nm_sample_ct * 2;
        uint32_t case_allele_obs_ct = nm_case_ct * 2;
        if (!is_regular_x) {
          if (is_nonx_haploid) {
            allele_obs_ct = nm_sample_ct;
            case_allele_obs_ct = nm_case_ct;
            // everything is on 0..1 scale, not 0..2
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              genotype_vals[sample_idx] *= 0.5;
            }
            const uint32_t high_ct = nm_sample_ct * allele_ct_m2;
            for (uint32_t uii = 0; uii != high_ct; ++uii) {
              multi_start[uii] *= 0.5;
            }
          }
        } else {
          CopyBitarrSubset(sex_male_collapsed, sample_nm, nm_sample_ct, tmp_nm);
          const uintptr_t* male_nm = tmp_nm;
          const uint32_t nm_male_ct = PopcountWords(male_nm, nm_sample_ctl);
          if (is_xchr_model_1) {
            // special case: multiply male values by 0.5
            uintptr_t sample_idx_base = 0;
            uintptr_t male_nm_bits = male_nm[0];
            for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
              const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
              genotype_vals[sample_idx] *= 0.5;
              // could insert multiallelic loop here isntead, but I'm guessing
              // that's worse due to locality of writes?
            }
            for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
              double* cur_start = &(multi_start[extra_allele_idx * nm_sample_ctav]);
              sample_idx_base = 0;
              male_nm_bits = male_nm[0];
              for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
                const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
                cur_start[sample_idx] *= 0.5;
              }
            }
            allele_obs_ct -= nm_male_ct;
            case_allele_obs_ct -= PopcountWordsIntersect(pheno_cc_nm, male_nm, nm_sample_ctl);
          }
        }
        const double mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, allele_ct);
        uint32_t nonomitted_allele_idx = 0;
        for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
          if (allele_idx == omitted_allele_idx) {
            continue;
          }

          double* geno_col = genotype_vals;
          if (allele_idx > (!omitted_allele_idx)) {
            geno_col = &(nm_predictors_pmaj_buf[(expected_predictor_ct - (allele_ct - allele_idx) + (allele_idx < omitted_allele_idx)) * nm_sample_ctav]);
          }
          double a1_dosage = u63tod(machr2_dosage_sums[allele_idx]) * kRecipDosageMid;
          if (is_xchr_model_1) {
            // ugh.
            a1_dosage = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              a1_dosage += geno_col[sample_idx];
            }
          } else {
            if (is_nonx_haploid) {
              a1_dosage *= 0.5;
            }
          }
          a1_dosages[allele_idx] = a1_dosage;

          // todo: shortcut if gcountcc computed and no dosages
          double a1_case_dosage = 0.0;
          uintptr_t sample_idx_base = 0;
          uintptr_t pheno_cc_nm_bits = pheno_cc_nm[0];
          for (uint32_t uii = 0; uii != nm_case_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(pheno_cc_nm, &sample_idx_base, &pheno_cc_nm_bits);
            a1_case_dosage += geno_col[sample_idx];
          }
          a1_case_dosages[allele_idx] = a1_case_dosage;
          block_aux_iter[nonomitted_allele_idx].sample_obs_ct = nm_sample_ct;
          block_aux_iter[nonomitted_allele_idx].allele_obs_ct = allele_obs_ct;
          if (!allele_ct_m2) {
            // Need main_dosage_sum and main_dosage_ssq for now (probably move
            // this computation in-place later).
            if (is_xchr_model_1) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = 0.0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const double cur_dosage = geno_col[sample_idx];
                main_dosage_ssq += cur_dosage * cur_dosage;
              }
            } else {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = u63tod(machr2_dosage_ssqs[allele_idx]) * kRecipDosageMidSq;
              if (is_nonx_haploid) {
                main_dosage_ssq *= 0.25;
              }
            }
          }
          block_aux_iter[nonomitted_allele_idx].a1_dosage = a1_dosage;

          // bugfix (4 Sep 2018): forgot to save this
          block_aux_iter[nonomitted_allele_idx].case_allele_obs_ct = case_allele_obs_ct;

          block_aux_iter[nonomitted_allele_idx].a1_case_dosage = a1_case_dosage;
          block_aux_iter[nonomitted_allele_idx].firth_fallback = 0;
          block_aux_iter[nonomitted_allele_idx].is_unfinished = 0;
          block_aux_iter[nonomitted_allele_idx].mach_r2 = mach_r2;
          ++nonomitted_allele_idx;
        }
        // Now free to skip the actual regression if there are too few samples,
        // or omitted allele corresponds to a zero-variance genotype column.
        // If another allele has zero variance but the omitted allele does not,
        // we now salvage as many alleles as we can.
        GlmErr glm_err = 0;
        if (nm_sample_ct <= expected_predictor_ct) {
          // reasonable for this to override CONST_ALLELE
          glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
        } else if (IsSet(const_alleles, omitted_allele_idx)) {
          glm_err = SetGlmErr0(kGlmErrcodeConstOmittedAllele);
        }
        if (glm_err) {
          if (missing_ct) {
            // covariates have not been copied yet, so we can't usually change
            // prev_nm from 0 to 1 when missing_ct == 0 (and there's little
            // reason to optimize the zero-covariate case)
            prev_nm = 0;
          }
          uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
          if (allele_ct_m2 && (beta_se_multiallelic_fused || (!hide_covar))) {
            reported_ct += allele_ct_m2;
          }
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            for (uint32_t uii = 0; uii != reported_ct; ++uii) {
              memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
              beta_se_iter[uii * 2 + 1] = -9.0;
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        } else {
          {
            double omitted_dosage = u63tod(allele_obs_ct);
            double omitted_case_dosage = u63tod(case_allele_obs_ct);
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if (allele_idx == omitted_allele_idx) {
                continue;
              }
              omitted_dosage -= a1_dosages[allele_idx];
              omitted_case_dosage -= a1_case_dosages[allele_idx];
            }
            a1_dosages[omitted_allele_idx] = omitted_dosage;
            a1_case_dosages[omitted_allele_idx] = omitted_case_dosage;
          }
          uint32_t parameter_uidx = 2 + domdev_present;
          double* nm_predictors_pmaj_istart = nullptr;
          // only need to do this part once per variant in multiallelic case
          double* nm_predictors_pmaj_iter = &(nm_predictors_pmaj_buf[nm_sample_ctav * (parameter_uidx - main_omitted)]);
          if (missing_ct || (!prev_nm)) {
            // fill phenotype
            uintptr_t sample_midx_base = 0;
            uintptr_t sample_nm_bits = sample_nm[0];
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
              nm_pheno_buf[sample_idx] = cur_pheno[sample_midx];
            }
            // bugfix (13 Oct 2017): must guarantee trailing phenotype values
            // are valid (exact contents don't matter since they are multiplied
            // by zero, but they can't be nan)
            ZeroDArr(nm_sample_ct_rem, &(nm_pheno_buf[nm_sample_ct]));

            // fill covariates
            for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx, ++parameter_uidx) {
              // unlike the float case, cur_covars_cmaj is NOT vector-aligned
              if (cur_parameter_subset && (!IsSet(cur_parameter_subset, parameter_uidx))) {
                continue;
              }
              const double* cur_covar_col;
              if (covar_idx < local_covar_ct) {
                cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
              } else {
                cur_covar_col = &(cur_covars_cmaj[(covar_idx - local_covar_ct) * cur_sample_ct]);
              }
              sample_midx_base = 0;
              sample_nm_bits = sample_nm[0];
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                *nm_predictors_pmaj_iter++ = cur_covar_col[sample_midx];
              }
              ZeromovDArr(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
            }
            nm_predictors_pmaj_istart = nm_predictors_pmaj_iter;
            // bugfix (13 Apr 2021): if local covariates are present, we can't
            // optimize as aggressively
            prev_nm = !(missing_ct || local_covar_ct);
          } else {
            // bugfix (15 Aug 2018): this was not handling --parameters
            // correctly when a covariate was only needed as part of an
            // interaction
            parameter_uidx += cur_covar_ct;
            nm_predictors_pmaj_istart = &(nm_predictors_pmaj_iter[literal_covar_ct * nm_sample_ctav]);
          }
          const uint32_t const_allele_ct = PopcountWords(const_alleles, allele_ctl);
          if (const_allele_ct) {
            // Must delete constant-allele columns from nm_predictors_pmaj, and
            // shift later columns back.
            double* read_iter = genotype_vals;
            double* write_iter = genotype_vals;
            for (uint32_t read_allele_idx = 0; read_allele_idx != allele_ct; ++read_allele_idx) {
              if (read_allele_idx == omitted_allele_idx) {
                continue;
              }
              if (!IsSet(const_alleles, read_allele_idx)) {
                if (write_iter != read_iter) {
                  memcpy(write_iter, read_iter, nm_sample_ctav * sizeof(double));
                }
                if (write_iter == genotype_vals) {
                  write_iter = multi_start;
                } else {
                  write_iter = &(write_iter[nm_sample_ctav]);
                }
              }
              if (read_iter == genotype_vals) {
                read_iter = multi_start;
              } else {
                read_iter = &(read_iter[nm_sample_ctav]);
              }
            }
          }
          const uint32_t cur_predictor_ct = expected_predictor_ct - const_allele_ct;
          const uint32_t cur_predictor_ctav = RoundUpPow2(cur_predictor_ct, kDoublePerDVec);
          const uint32_t cur_predictor_ctavp1 = cur_predictor_ctav + 1;
          uint32_t nonconst_extra_regression_idx = UINT32_MAX;  // deliberate overflow
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            double* main_vals = &(nm_predictors_pmaj_buf[nm_sample_ctav]);
            double* domdev_vals = nullptr;
            uint32_t is_unfinished = 0;
            uint32_t is_residualized = 0;
            // _stop instead of _ct since, in the residualized case, the
            // intercept (predictor index 0) is not included; we iterate over
            // the predictor indices in [1, _stop).
            uint32_t cur_regressed_predictor_stop = cur_predictor_ct;
            uint32_t cur_regressed_predictor_ctav = cur_predictor_ctav;
            uint32_t cur_regressed_predictor_ctavp1 = cur_predictor_ctavp1;
            uint32_t cur_biallelic_regressed_predictor_stop = cur_biallelic_predictor_ct;
            if (extra_regression_ct) {
              if (IsSet(const_alleles, extra_regression_idx + (extra_regression_idx >= omitted_allele_idx))) {
                glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
                goto GlmLogisticThreadD_skip_regression;
              }
              ++nonconst_extra_regression_idx;
              if (nonconst_extra_regression_idx) {
                double* swap_target = &(multi_start[(nonconst_extra_regression_idx - 1) * nm_sample_ctav]);
                for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
                  double dxx = genotype_vals[uii];
                  genotype_vals[uii] = swap_target[uii];
                  swap_target[uii] = dxx;
                }
              }
            }
            if (main_omitted) {
              // if main_mutated, this will be filled below
              // if not, this aliases genotype_vals
              main_vals = &(nm_predictors_pmaj_buf[(cur_predictor_ct + main_mutated) * nm_sample_ctav]);
            } else if (joint_genotypic || joint_hethom) {
              // in hethom case, do this before clobbering genotype data
              domdev_vals = &(main_vals[nm_sample_ctav]);
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                double cur_genotype_val = genotype_vals[sample_idx];
                if (cur_genotype_val > 1.0) {
                  cur_genotype_val = 2.0 - cur_genotype_val;
                }
                domdev_vals[sample_idx] = cur_genotype_val;
              }
              ZeroDArr(nm_sample_ct_rem, &(domdev_vals[nm_sample_ct]));
            }
            if (model_dominant) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                double cur_genotype_val = genotype_vals[sample_idx];
                // 0..1..1
                if (cur_genotype_val > 1.0) {
                  cur_genotype_val = 1.0;
                }
                main_vals[sample_idx] = cur_genotype_val;
              }
            } else if (model_recessive || joint_hethom) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                double cur_genotype_val = genotype_vals[sample_idx];
                // 0..0..1
                if (cur_genotype_val < 1.0) {
                  cur_genotype_val = 0.0;
                } else {
                  cur_genotype_val -= 1.0;
                }
                main_vals[sample_idx] = cur_genotype_val;
              }
            } else if (model_hetonly) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                double cur_genotype_val = genotype_vals[sample_idx];
                // 0..1..0
                if (cur_genotype_val > 1.0) {
                  cur_genotype_val = 2.0 - cur_genotype_val;
                }
                main_vals[sample_idx] = cur_genotype_val;
              }
            }

            // fill interaction terms
            if (add_interactions) {
              nm_predictors_pmaj_iter = nm_predictors_pmaj_istart;
              for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
                const double* cur_covar_col;
                if (covar_idx < local_covar_ct) {
                  cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                } else {
                  cur_covar_col = &(cur_covars_cmaj[covar_idx * cur_sample_ct]);
                }
                if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                  uintptr_t sample_midx_base = 0;
                  uintptr_t sample_nm_bits = sample_nm[0];
                  for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                    const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                    *nm_predictors_pmaj_iter++ = main_vals[sample_idx] * cur_covar_col[sample_midx];
                  }
                  ZeromovDArr(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
                }
                ++parameter_uidx;
                if (domdev_present) {
                  if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                    uintptr_t sample_midx_base = 0;
                    uintptr_t sample_nm_bits = sample_nm[0];
                    for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                      const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                      *nm_predictors_pmaj_iter++ = domdev_vals[sample_idx] * cur_covar_col[sample_midx];
                    }
                    ZeromovDArr(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
                  }
                  ++parameter_uidx;
                }
              }
            }
            if (corr_inv && prev_nm && (!allele_ct_m2)) {
              uintptr_t start_pred_idx = 0;
              if (!(model_dominant || model_recessive || model_hetonly || joint_hethom)) {
                start_pred_idx = domdev_present + 2;
                semicomputed_biallelic_xtx[cur_predictor_ct] = main_dosage_sum;
                semicomputed_biallelic_xtx[cur_predictor_ct + 1] = main_dosage_ssq;
              }
              if (cur_predictor_ct > start_pred_idx) {
                ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[nm_sample_ctav]), &(nm_predictors_pmaj_buf[start_pred_idx * nm_sample_ctav]), nm_sample_ct, nm_sample_ctav, cur_predictor_ct - start_pred_idx, &(predictor_dotprod_buf[start_pred_idx]));
                for (uint32_t uii = start_pred_idx; uii != cur_predictor_ct; ++uii) {
                  semicomputed_biallelic_xtx[cur_predictor_ct + uii] = predictor_dotprod_buf[uii];
                }
              }
              if (domdev_present) {
                ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[2 * nm_sample_ctav]), nm_predictors_pmaj_buf, nm_sample_ct, nm_sample_ctav, cur_predictor_ct, predictor_dotprod_buf);
                for (uint32_t uii = 0; uii != cur_predictor_ct; ++uii) {
                  semicomputed_biallelic_xtx[2 * cur_predictor_ct + uii] = predictor_dotprod_buf[uii];
                }
                semicomputed_biallelic_xtx[cur_predictor_ct + 2] = semicomputed_biallelic_xtx[2 * cur_predictor_ct + 1];
              }
              glm_err = CheckMaxCorrAndVifNm(semicomputed_biallelic_xtx, corr_inv, cur_predictor_ct, domdev_present_p1, cur_sample_ct_recip, cur_sample_ct_m1_recip, max_corr, vif_thresh, semicomputed_biallelic_corr_matrix, semicomputed_biallelic_inv_corr_sqrts, dbl_2d_buf, &(dbl_2d_buf[2 * cur_predictor_ct]), &(dbl_2d_buf[3 * cur_predictor_ct]));
              if (glm_err) {
                goto GlmLogisticThreadD_skip_regression;
              }
            } else {
              MultiplySelfTransposeStrided(&(nm_predictors_pmaj_buf[nm_sample_ctav]), cur_predictor_ct - 1, nm_sample_ct, nm_sample_ctav, predictor_dotprod_buf);
              for (uintptr_t pred_idx = 1; pred_idx != cur_predictor_ct; ++pred_idx) {
                const double* predictor_row = &(nm_predictors_pmaj_buf[pred_idx * nm_sample_ctav]);
                double row_sum = 0.0;
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  row_sum += predictor_row[sample_idx];
                }
                dbl_2d_buf[pred_idx - 1] = row_sum;
              }
              glm_err = CheckMaxCorrAndVif(predictor_dotprod_buf, 0, cur_predictor_ct - 1, nm_sample_ct, max_corr, vif_thresh, dbl_2d_buf, nullptr, inverse_corr_buf, inv_1d_buf);
              if (glm_err) {
                goto GlmLogisticThreadD_skip_regression;
              }
            }
            ZeroDArr(cur_predictor_ctav, coef_return);
            if (!cur_is_always_firth) {
              // Does any genotype column have zero case or zero control
              // dosage?  If yes, faster to skip logistic regression than
              // wait for convergence failure.
              for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                if (IsSet(const_alleles, allele_idx)) {
                  continue;
                }
                const double tot_dosage = a1_dosages[allele_idx];
                const double case_dosage = a1_case_dosages[allele_idx];
                if ((case_dosage == 0.0) || (case_dosage == tot_dosage)) {
                  if (is_sometimes_firth) {
                    goto GlmLogisticThreadD_firth_fallback;
                  }
                  glm_err = SetGlmErr1(kGlmErrcodeSeparation, allele_idx);
                  goto GlmLogisticThreadD_skip_regression;
                }
              }
              if (!cur_cc_residualize) {
                if (LogisticRegressionD(nm_pheno_buf, nm_predictors_pmaj_buf, nullptr, nm_sample_ct, cur_predictor_ct, &is_unfinished, coef_return, cholesky_decomp_return, pp_buf, sample_variance_buf, hh_return, gradient_buf, dcoef_buf, inv_1d_buf, dbl_2d_buf)) {
                  if (is_sometimes_firth) {
                    ZeroDArr(cur_predictor_ctav, coef_return);
                    goto GlmLogisticThreadD_firth_fallback;
                  }
                  glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
                  goto GlmLogisticThreadD_skip_regression;
                }
              } else {
                if (LogisticRegressionResidualizedD(nm_pheno_buf, nm_predictors_pmaj_buf, sample_nm, cur_cc_residualize, nm_sample_ct, cur_predictor_ct, &is_unfinished, coef_return, cholesky_decomp_return, inv_1d_buf, dbl_2d_buf, pp_buf, sample_variance_buf, hh_return, gradient_buf, dcoef_buf, mean_centered_pmaj_buf, sample_offsets_buf)) {
                  if (is_sometimes_firth) {
                    ZeroDArr(cur_predictor_ctav, coef_return);
                    goto GlmLogisticThreadD_firth_fallback;
                  }
                  glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
                  goto GlmLogisticThreadD_skip_regression;
                }
                is_residualized = 1;
                cur_regressed_predictor_stop = domdev_present + allele_ct;
                cur_regressed_predictor_ctav = RoundUpPow2(cur_regressed_predictor_stop, kDoublePerDVec);
                cur_regressed_predictor_ctavp1 = cur_regressed_predictor_ctav + 1;
                cur_biallelic_regressed_predictor_stop = domdev_present + 2;
              }
              // unlike FirthRegressionD(), hh_return isn't inverted yet, do
              // that here
              for (uint32_t pred_uidx = is_residualized; pred_uidx != cur_regressed_predictor_stop; ++pred_uidx) {
                double* hh_inv_row = &(hh_return[pred_uidx * cur_regressed_predictor_ctav]);
                // ZeroDArr(cur_regressed_predictor_stop, gradient_buf);
                // gradient_buf[pred_uidx] = 1.0;
                // (y is gradient_buf, x is dcoef_buf)
                // SolveLinearSystemD(cholesky_decomp_return, &(gradient_buf[is_residualized]), cur_regressed_predictor_stop - is_residualized, &(hh_inv_row[is_residualized]));
                // that works, but doesn't exploit the sparsity of y

                // hh_return does now have vector-aligned rows
                ZeroDArr(pred_uidx, hh_inv_row);

                double dxx = 1.0;
                for (uint32_t row_idx = pred_uidx; row_idx != cur_regressed_predictor_stop; ++row_idx) {
                  const double* ll_row = &(cholesky_decomp_return[row_idx * cur_regressed_predictor_ctav]);
                  for (uint32_t col_idx = pred_uidx; col_idx != row_idx; ++col_idx) {
                    dxx -= ll_row[col_idx] * hh_inv_row[col_idx];
                  }
                  hh_inv_row[row_idx] = dxx / ll_row[row_idx];
                  dxx = 0.0;
                }
                for (uint32_t col_idx = cur_regressed_predictor_stop; col_idx > is_residualized; ) {
                  dxx = hh_inv_row[--col_idx];
                  double* hh_inv_row_iter = &(hh_inv_row[cur_regressed_predictor_stop - 1]);
                  for (uint32_t row_idx = cur_regressed_predictor_stop - 1; row_idx > col_idx; --row_idx) {
                    dxx -= cholesky_decomp_return[row_idx * cur_regressed_predictor_ctav + col_idx] * (*hh_inv_row_iter--);
                  }
                  *hh_inv_row_iter = dxx / cholesky_decomp_return[col_idx * cur_regressed_predictor_ctavp1];
                }
              }
            } else {
              if (!is_always_firth) {
              GlmLogisticThreadD_firth_fallback:
                block_aux_iter[extra_regression_idx].firth_fallback = 1;
                if (allele_ct_m2 && beta_se_multiallelic_fused) {
                  for (uint32_t uii = 1; uii != allele_ct - 1; ++uii) {
                    block_aux_iter[uii].firth_fallback = 1;
                  }
                }
              }
              if (!cur_cc_residualize) {
#ifndef NDEBUG
                if (g_debug_on && (variant_uidx == 0)) {
                  // Dump phenotype values and predictor matrix to .inputs .
                  char fname_buf[kPglFnamesize + 8];
                  const uint32_t outname_slen = strlen(common->outname);
                  char* fname_write_iter = memcpya(fname_buf, common->outname, outname_slen);
                  strcpy_k(fname_write_iter, ".inputs");
                  FILE* outfile;
                  if (unlikely(fopen_checked(fname_buf, FOPEN_WB, &outfile))) {
                    fprintf(stderr, "\nPanic: Failed to open .inputs file for writing.\n");
                    exit(S_CAST(int, kPglRetOpenFail));
                  }
                  for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                    fprintf(outfile, "%g ", nm_pheno_buf[sample_idx]);
                    for (uintptr_t pred_idx = 1; pred_idx != cur_predictor_ct; ++pred_idx) {
                      fprintf(outfile, " %g", nm_predictors_pmaj_buf[nm_sample_ctav * pred_idx + sample_idx]);
                    }
                    fprintf(outfile, "\n");
                  }
                  fclose(outfile);
                }
#endif
                if (FirthRegressionD(nm_pheno_buf, nm_predictors_pmaj_buf, nullptr, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, hh_return, inv_1d_buf, dbl_2d_buf, pp_buf, sample_variance_buf, gradient_buf, dcoef_buf, hdiag_buf, score_buf, hh0_buf, tmpnxk_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
                  goto GlmLogisticThreadD_skip_regression;
                }
              } else {
                if (FirthRegressionResidualizedD(nm_pheno_buf, nm_predictors_pmaj_buf, sample_nm, cur_cc_residualize, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, hh_return, inv_1d_buf, dbl_2d_buf, pp_buf, sample_variance_buf, gradient_buf, dcoef_buf, hdiag_buf, score_buf, hh0_buf, tmpnxk_buf, mean_centered_pmaj_buf, sample_offsets_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
                  goto GlmLogisticThreadD_skip_regression;
                }
                is_residualized = 1;
                cur_regressed_predictor_stop = domdev_present + allele_ct;
                cur_regressed_predictor_ctav = RoundUpPow2(cur_regressed_predictor_stop, kDoublePerDVec);
                cur_regressed_predictor_ctavp1 = cur_regressed_predictor_ctav + 1;
                cur_biallelic_regressed_predictor_stop = domdev_present + 2;
              }
            }
            // validParameters() check
            for (uint32_t pred_uidx = 1; pred_uidx != cur_regressed_predictor_stop; ++pred_uidx) {
              const double hh_inv_diag_element = hh_return[pred_uidx * cur_regressed_predictor_ctavp1];
              if ((hh_inv_diag_element < 1e-20) || (!isfinite_d(hh_inv_diag_element))) {
                glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                goto GlmLogisticThreadD_skip_regression;
              }
              // use sample_variance_buf[] to store diagonal square roots
              sample_variance_buf[pred_uidx] = sqrt(hh_inv_diag_element);
            }
            if (!is_residualized) {
              sample_variance_buf[0] = sqrt(hh_return[0]);
            }
            for (uint32_t pred_uidx = 1 + is_residualized; pred_uidx != cur_regressed_predictor_stop; ++pred_uidx) {
              const double cur_hh_inv_diag_sqrt = 0.99999 * sample_variance_buf[pred_uidx];
              const double* hh_inv_row_iter = &(hh_return[pred_uidx * cur_regressed_predictor_ctav + is_residualized]);
              const double* hh_inv_diag_sqrts_iter = &(sample_variance_buf[is_residualized]);
              for (uint32_t pred_uidx2 = is_residualized; pred_uidx2 != pred_uidx; ++pred_uidx2) {
                if ((*hh_inv_row_iter++) > cur_hh_inv_diag_sqrt * (*hh_inv_diag_sqrts_iter++)) {
                  glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                  goto GlmLogisticThreadD_skip_regression;
                }
              }
            }
            if (is_unfinished) {
              block_aux_iter[extra_regression_idx].is_unfinished = 1;
              if (allele_ct_m2 && beta_se_multiallelic_fused) {
                for (uint32_t uii = 1; uii != allele_ct - 1; ++uii) {
                  block_aux_iter[uii].is_unfinished = 1;
                }
              }
            }
            {
              double* beta_se_iter2 = beta_se_iter;
              for (uint32_t pred_uidx = reported_pred_uidx_start; pred_uidx != reported_pred_uidx_biallelic_end; ++pred_uidx) {
                // In the multiallelic-fused case, if the first allele is
                // constant, this writes the beta/se values for the first
                // nonconstant, non-omitted allele where the results for the
                // first allele belong.  We correct that at the end of this
                // block.
                *beta_se_iter2++ = coef_return[pred_uidx];
                *beta_se_iter2++ = sample_variance_buf[pred_uidx];
              }
              if (cur_constraint_ct) {
                // bugfix (4 Sep 2021): forgot to update this for residualize
                // case
                *beta_se_iter2++ = 0.0;

                uint32_t joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_regressed_predictor_stop + joint_test_idx] = 1.0;
                }
                double chisq;
                if (!LinearHypothesisChisq(coef_return, cur_constraints_con_major, hh_return, cur_constraint_ct, cur_regressed_predictor_stop, cur_regressed_predictor_ctav, &chisq, tmphxs_buf, h_transpose_buf, inner_buf, inv_1d_buf, outer_buf)) {
                  *beta_se_iter2++ = chisq;
                } else {
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeRankDeficient);
                  memcpy(&(beta_se_iter2[-1]), &glm_err2, 8);
                  *beta_se_iter2++ = -9.0;
                }
                // next test may have different alt allele count
                joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_regressed_predictor_stop + joint_test_idx] = 0.0;
                }
              }
              if (!const_allele_ct) {
                if (beta_se_multiallelic_fused || (!hide_covar)) {
                  for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
                    *beta_se_iter2++ = coef_return[cur_biallelic_regressed_predictor_stop + extra_allele_idx];
                    *beta_se_iter2++ = sample_variance_buf[cur_biallelic_regressed_predictor_stop + extra_allele_idx];
                  }
                }
              } else if (!beta_se_multiallelic_fused) {
                if (!hide_covar) {
                  // Need to insert some {CONST_ALLELE, -9} entries.
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                  const uint32_t cur_raw_allele_idx = extra_regression_idx + (extra_regression_idx >= omitted_allele_idx);
                  uint32_t extra_read_allele_idx = 0;
                  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    if ((allele_idx == omitted_allele_idx) || (allele_idx == cur_raw_allele_idx)) {
                      continue;
                    }
                    if (IsSet(const_alleles, allele_idx)) {
                      memcpy(beta_se_iter2, &glm_err2, 8);
                      beta_se_iter2[1] = -9.0;
                      beta_se_iter2 = &(beta_se_iter2[2]);
                    } else {
                      *beta_se_iter2++ = coef_return[cur_biallelic_regressed_predictor_stop + extra_read_allele_idx];
                      *beta_se_iter2++ = sample_variance_buf[cur_biallelic_regressed_predictor_stop + extra_read_allele_idx];
                      ++extra_read_allele_idx;
                    }
                  }
                }
              } else {
                const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                // Special-case first nonconst allele since it's positioned
                // discontinuously, and its BETA/SE may already be correctly
                // filled.
                uint32_t allele_idx = omitted_allele_idx? 0 : 1;
                if (IsSet(const_alleles, allele_idx)) {
                  memcpy(&(beta_se_iter[2 * include_intercept]), &glm_err2, 8);
                  beta_se_iter[2 * include_intercept + 1] = -9.0;
                  allele_idx = AdvTo0Bit(const_alleles, 1);
                  if (allele_idx == omitted_allele_idx) {
                    allele_idx = AdvTo0Bit(const_alleles, omitted_allele_idx + 1);
                  }
                  const uint32_t skip_ct = allele_idx - 1 - (allele_idx > omitted_allele_idx);
                  for (uint32_t uii = 0; uii != skip_ct; ++uii) {
                    memcpy(beta_se_iter2, &glm_err2, 8);
                    beta_se_iter2[1] = -9.0;
                    beta_se_iter2 = &(beta_se_iter2[2]);
                  }
                  *beta_se_iter2++ = coef_return[1];
                  *beta_se_iter2++ = sample_variance_buf[1];
                }
                ++allele_idx;
                uint32_t nonconst_allele_idx_m1 = 0;
                for (; allele_idx != allele_ct; ++allele_idx) {
                  if (allele_idx == omitted_allele_idx) {
                    continue;
                  }
                  if (!IsSet(const_alleles, allele_idx)) {
                    *beta_se_iter2++ = coef_return[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                    *beta_se_iter2++ = sample_variance_buf[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                    ++nonconst_allele_idx_m1;
                  } else {
                    memcpy(beta_se_iter2, &glm_err2, 8);
                    beta_se_iter2[1] = -9.0;
                    beta_se_iter2 = &(beta_se_iter2[2]);
                  }
                }
              }
            }
            while (0) {
            GlmLogisticThreadD_skip_regression:
              {
                uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
                if (allele_ct_m2 && (beta_se_multiallelic_fused || (!hide_covar))) {
                  reported_ct += allele_ct_m2;
                }
                for (uint32_t uii = 0; uii != reported_ct; ++uii) {
                  memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
                  beta_se_iter[uii * 2 + 1] = -9.0;
                }
              }
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        }
        block_aux_iter = &(block_aux_iter[allele_ct - 1]);
        if (local_covars_iter) {
          local_covars_iter = &(local_covars_iter[local_covar_ct * max_sample_ct]);
        }
      }
    }
    parity = 1 - parity;
    variant_idx_offset += cur_block_variant_ct;
    while (0) {
    GlmLogisticThreadD_err:
      UpdateU64IfSmaller(new_err_info, &common->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

// only pass the parameters which aren't also needed by the compute threads,
// for now
// valid_variants and valid_alleles are a bit redundant, may want to remove the
// former later, but let's make that decision during/after permutation test
// implementation
PglErr GlmLogistic(const char* cur_pheno_name, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const char* outname, uint32_t raw_variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLogisticCtx* ctx, TextStream* local_covar_txsp, LlStr** outfnames_ll_ptr, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, double* orig_permstat, uintptr_t* valid_allele_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  ThreadGroup tg;
  PreinitCstream(&css);
  PreinitThreads(&tg);
  {
    GlmCtx* common = ctx->common;
    const uintptr_t* variant_include = common->variant_include;
    const ChrInfo* cip = common->cip;
    const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
    const AlleleCode* omitted_alleles = common->omitted_alleles;

    const uint32_t sample_ct = common->sample_ct;
    const uint32_t sample_ct_x = common->sample_ct_x;
    const uint32_t sample_ct_y = common->sample_ct_y;
    const uint32_t covar_ct = common->covar_ct;
    const uintptr_t local_covar_ct = common->local_covar_ct;
    const uint32_t covar_ct_x = common->covar_ct_x;
    const uint32_t covar_ct_y = common->covar_ct_y;

    uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
    if (max_sample_ct < sample_ct_y) {
      max_sample_ct = sample_ct_y;
    }
    // obvious todo: wrap these in structs
    uint32_t* local_sample_idx_order = nullptr;
    uint32_t local_line_idx = 0;
    uint32_t local_xy = 0;  // 1 = chrX, 2 = chrY

    const char* local_line_iter = nullptr;
    uint32_t local_prev_chr_code = UINT32_MAX;
    uint32_t local_chr_code = UINT32_MAX;
    uint32_t local_bp = UINT32_MAX;
    uint32_t local_skip_chr = 1;
    if (local_covar_ct) {
      reterr = TextRewind(local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLogistic_ret_TSTREAM_FAIL;
      }
      local_line_idx = glm_info_ptr->local_header_line_ct;
      reterr = TextSkip(local_line_idx, local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLogistic_ret_TSTREAM_FAIL;
      }
      if (unlikely(bigstack_alloc_u32(local_sample_ct, &local_sample_idx_order))) {
        goto GlmLogistic_ret_NOMEM;
      }
      for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
        const uint32_t cur_uidx = local_sample_uidx_order[uii];
        uint32_t cur_idx = UINT32_MAX;
        if ((cur_uidx != UINT32_MAX) && IsSet(common->sample_include, cur_uidx)) {
          cur_idx = RawToSubsettedPos(common->sample_include, common->sample_include_cumulative_popcounts, cur_uidx);
        }
        local_sample_idx_order[uii] = cur_idx;
      }
    }

    const uint32_t variant_ct = common->variant_ct;

    const GlmFlags glm_flags = glm_info_ptr->flags;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    // forced-singlethreaded
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto GlmLogistic_ret_1;
    }
    if (outfnames_ll_ptr) {
      if (unlikely(PushLlStr(outname, outfnames_ll_ptr))) {
        goto GlmLogistic_ret_NOMEM;
      }
    }
    const uint32_t report_neglog10p = (glm_flags / kfGlmLog10) & 1;
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t domdev_present_p1 = domdev_present + 1;

    const uint32_t constraint_ct = common->constraint_ct;
    const uint32_t constraint_ct_x = common->constraint_ct_x;
    const uint32_t constraint_ct_y = common->constraint_ct_y;

    const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
    uint32_t biallelic_predictor_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
    uint32_t biallelic_predictor_ct_x = 2 + domdev_present + covar_ct_x * (1 + add_interactions * domdev_present_p1);
    uint32_t biallelic_predictor_ct_y = 2 + domdev_present + covar_ct_y * (1 + add_interactions * domdev_present_p1);
    const uintptr_t* parameter_subset = common->parameter_subset;
    const uintptr_t* parameter_subset_x = common->parameter_subset_x;
    const uintptr_t* parameter_subset_y = common->parameter_subset_y;
    if (parameter_subset) {
      biallelic_predictor_ct = PopcountWords(parameter_subset, BitCtToWordCt(biallelic_predictor_ct));
      if (sample_ct_x) {
        biallelic_predictor_ct_x = PopcountWords(parameter_subset_x, BitCtToWordCt(biallelic_predictor_ct_x));
      } else {
        biallelic_predictor_ct_x = 0;
      }
      if (sample_ct_y) {
        biallelic_predictor_ct_y = PopcountWords(parameter_subset_y, BitCtToWordCt(biallelic_predictor_ct_x));
      } else {
        biallelic_predictor_ct_y = 0;
      }
    }
    uint32_t biallelic_reported_test_ct = GetBiallelicReportedTestCt(parameter_subset, glm_flags, covar_ct, common->tests_flag);
    uintptr_t max_reported_test_ct = biallelic_reported_test_ct;
    uint32_t biallelic_reported_test_ct_x = 0;
    if (sample_ct_x) {
      biallelic_reported_test_ct_x = GetBiallelicReportedTestCt(parameter_subset_x, glm_flags, covar_ct_x, common->tests_flag);
      if (biallelic_reported_test_ct_x > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_x;
      }
    }
    uint32_t biallelic_reported_test_ct_y = 0;
    if (sample_ct_y) {
      biallelic_reported_test_ct_y = GetBiallelicReportedTestCt(parameter_subset_y, glm_flags, covar_ct_y, common->tests_flag);
      if (biallelic_reported_test_ct_y > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_y;
      }
    }
    const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
    const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
    const GlmColFlags glm_cols = glm_info_ptr->cols;
    const uint32_t test_col = glm_cols & kfGlmColTest;
    if (unlikely((!test_col) && (max_reported_test_ct > 1))) {
      // this is okay in plain multiallelic case due to A1 column
      logerrputs("Error: --glm's 'test' column cannot be omitted when results for multiple\npredictors are reported.  (Did you forget 'hide-covar'?)\n");
      goto GlmLogistic_ret_INCONSISTENT_INPUT;
    }
    const uint32_t main_mutated = ((glm_flags & (kfGlmDominant | kfGlmRecessive | kfGlmHetonly | kfGlmHethom)) != kfGlm0);
    // if 'fused', one row per variant
    // otherwise, one row per tested allele
    const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!main_mutated) && (!common->tests_flag) && (!add_interactions);
    if (beta_se_multiallelic_fused || (!hide_covar)) {
      max_reported_test_ct += max_extra_allele_ct;
    }
    common->max_reported_test_ct = max_reported_test_ct;

    const uint32_t is_sometimes_firth = !(glm_flags & kfGlmNoFirth);
    const uint32_t is_always_firth = (glm_flags / kfGlmFirth) & 1;
    const uint32_t is_cc_residualize = !!(glm_flags & (kfGlmFirthResidualize | kfGlmCcResidualize));

    uint32_t x_code = UINT32_MAXM1;
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    if (sample_ct_x) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetX, &x_code, &x_start, &x_end);
    }
    uint32_t y_code = UINT32_MAXM1;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    if (sample_ct_y) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetY, &y_code, &y_start, &y_end);
    }
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t chr_col = glm_cols & kfGlmColChrom;

    // includes trailing tab
    char* chr_buf = nullptr;
    if (chr_col) {
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto GlmLogistic_ret_NOMEM;
      }
    }

    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > variant_ct) {
      calc_thread_ct = variant_ct;
    }

    const uint32_t is_single_prec = (glm_flags / kfGlmSinglePrecCc) & 1;
    const uint32_t main_omitted = parameter_subset && (!IsSet(parameter_subset, 1));
    const uint32_t xmain_ct = main_mutated + main_omitted;
    const uint32_t gcount_cc_col = glm_cols & kfGlmColGcountcc;
    // workflow is similar to --make-bed
    uintptr_t workspace_alloc;
    if (is_single_prec) {
      workspace_alloc = GetLogisticWorkspaceSizeF(sample_ct, biallelic_predictor_ct, domdev_present_p1, max_extra_allele_ct, constraint_ct, xmain_ct, gcount_cc_col, is_sometimes_firth, is_cc_residualize);
      if (sample_ct_x) {
        const uintptr_t workspace_alloc_x = GetLogisticWorkspaceSizeF(sample_ct_x, biallelic_predictor_ct_x, domdev_present_p1, max_extra_allele_ct, constraint_ct_x, xmain_ct, gcount_cc_col, is_sometimes_firth, is_cc_residualize);
        if (workspace_alloc_x > workspace_alloc) {
          workspace_alloc = workspace_alloc_x;
        }
      }
      if (sample_ct_y) {
        const uintptr_t workspace_alloc_y = GetLogisticWorkspaceSizeF(sample_ct_y, biallelic_predictor_ct_y, domdev_present_p1, max_extra_allele_ct, constraint_ct_y, xmain_ct, gcount_cc_col, is_sometimes_firth, is_cc_residualize);
        if (workspace_alloc_y > workspace_alloc) {
          workspace_alloc = workspace_alloc_y;
        }
      }
    } else {
      workspace_alloc = GetLogisticWorkspaceSizeD(sample_ct, biallelic_predictor_ct, domdev_present_p1, max_extra_allele_ct, constraint_ct, xmain_ct, gcount_cc_col, is_sometimes_firth, is_cc_residualize);
      if (sample_ct_x) {
        const uintptr_t workspace_alloc_x = GetLogisticWorkspaceSizeD(sample_ct_x, biallelic_predictor_ct_x, domdev_present_p1, max_extra_allele_ct, constraint_ct_x, xmain_ct, gcount_cc_col, is_sometimes_firth, is_cc_residualize);
        if (workspace_alloc_x > workspace_alloc) {
          workspace_alloc = workspace_alloc_x;
        }
      }
      if (sample_ct_y) {
        const uintptr_t workspace_alloc_y = GetLogisticWorkspaceSizeD(sample_ct_y, biallelic_predictor_ct_y, domdev_present_p1, max_extra_allele_ct, constraint_ct_y, xmain_ct, gcount_cc_col, is_sometimes_firth, is_cc_residualize);
        if (workspace_alloc_y > workspace_alloc) {
          workspace_alloc = workspace_alloc_y;
        }
      }
    }
    // +1 is for top-level common->workspace_bufs
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;

    uintptr_t per_variant_xalloc_byte_ct = max_sample_ct * local_covar_ct * (2 - is_single_prec) * sizeof(float);
    uintptr_t per_alt_allele_xalloc_byte_ct = sizeof(LogisticAuxResult);
    if (beta_se_multiallelic_fused) {
      per_variant_xalloc_byte_ct += 2 * max_reported_test_ct * sizeof(double);
    } else {
      per_alt_allele_xalloc_byte_ct += 2 * max_reported_test_ct * sizeof(double);
    }
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    common->thread_mhc = nullptr;
    common->dosage_presents = nullptr;
    common->dosage_mains = nullptr;
    uint32_t read_block_size;
    uintptr_t max_alt_allele_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, max_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, per_alt_allele_xalloc_byte_ct, pgfip, &calc_thread_ct, &common->genovecs, max_extra_allele_ct? (&common->thread_mhc) : nullptr, nullptr, nullptr, dosage_is_present? (&common->dosage_presents) : nullptr, dosage_is_present? (&common->dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, &max_alt_allele_block_size, main_loadbufs, &common->pgr_ptrs, &common->read_variant_uidx_starts))) {
      goto GlmLogistic_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto GlmLogistic_ret_NOMEM;
    }
    LogisticAuxResult* logistic_block_aux_bufs[2];
    double* block_beta_se_bufs[2];

    for (uint32_t uii = 0; uii != 2; ++uii) {
      if (unlikely(BIGSTACK_ALLOC_X(LogisticAuxResult, max_alt_allele_block_size, &(logistic_block_aux_bufs[uii])))) {
        goto GlmLogistic_ret_NOMEM;
      }
      if (beta_se_multiallelic_fused) {
        if (unlikely(bigstack_alloc_d(read_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii])))) {
          goto GlmLogistic_ret_NOMEM;
        }
      } else {
        if (unlikely(bigstack_alloc_d(max_alt_allele_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii])))) {
          goto GlmLogistic_ret_NOMEM;
        }
      }
      if (local_covar_ct) {
        // bugfix (18 May 2018): don't want sizeof(float) here
        if (is_single_prec) {
          if (unlikely(bigstack_alloc_f(read_block_size * max_sample_ct * local_covar_ct, &(ctx->local_covars_vcmaj_f[uii])))) {
            goto GlmLogistic_ret_NOMEM;
          }
        } else {
          if (unlikely(bigstack_alloc_d(read_block_size * max_sample_ct * local_covar_ct, &(ctx->local_covars_vcmaj_d[uii])))) {
            goto GlmLogistic_ret_NOMEM;
          }
        }
      } else {
        ctx->local_covars_vcmaj_f[uii] = nullptr;
        ctx->local_covars_vcmaj_d[uii] = nullptr;
      }
    }

    if (is_single_prec && (max_sample_ct > 100000)) {
      logerrputs("Warning: --glm's 'single-prec-cc' mode is not recommended on more than ~100000\nsamples.\n");
    }
    common->workspace_bufs = S_CAST(unsigned char**, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t)));
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      common->workspace_bufs[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(workspace_alloc));
    }
    common->err_info = (~0LLU) << 32;
#ifndef NDEBUG
    // temporary debug
    common->outname = outname;
#endif
    SetThreadFuncAndData(is_single_prec? GlmLogisticThreadF : GlmLogisticThreadD, ctx, &tg);

    const uint32_t ref_col = glm_cols & kfGlmColRef;
    const uint32_t alt1_col = glm_cols & kfGlmColAlt1;
    const uint32_t alt_col = glm_cols & kfGlmColAlt;
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t all_nonref = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags, glm_cols / kfGlmColMaybeprovref, raw_variant_ct, all_nonref);
    const uint32_t omitted_col = glm_cols & kfGlmColOmitted;
    const uint32_t ax_col = glm_cols & kfGlmColAx;
    const uint32_t a1_ct_col = glm_cols & kfGlmColA1count;
    const uint32_t tot_allele_col = glm_cols & kfGlmColTotallele;
    const uint32_t a1_ct_cc_col = glm_cols & kfGlmColA1countcc;
    const uint32_t tot_allele_cc_col = glm_cols & kfGlmColTotallelecc;
    const uint32_t a1_freq_col = glm_cols & kfGlmColA1freq;
    const uint32_t a1_freq_cc_col = glm_cols & kfGlmColA1freqcc;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t firth_yn_col = (glm_cols & kfGlmColFirthYn) && is_sometimes_firth && (!is_always_firth);
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t orbeta_col = glm_cols & (kfGlmColBeta | kfGlmColOrbeta);
    const uint32_t report_beta_instead_of_odds_ratio = glm_cols & kfGlmColBeta;
    const uint32_t se_col = glm_cols & kfGlmColSe;
    const uint32_t ci_col = (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t z_col = glm_cols & kfGlmColTz;
    const uint32_t p_col = glm_cols & kfGlmColP;
    const uint32_t err_col = glm_cols & kfGlmColErr;
    *cswritep++ = '#';
    if (chr_col) {
      cswritep = strcpya_k(cswritep, "CHROM\t");
    }
    if (variant_bps) {
      cswritep = strcpya_k(cswritep, "POS\t");
    }
    cswritep = strcpya_k(cswritep, "ID");
    if (ref_col) {
      cswritep = strcpya_k(cswritep, "\tREF");
    }
    if (alt1_col) {
      cswritep = strcpya_k(cswritep, "\tALT1");
    }
    if (alt_col) {
      cswritep = strcpya_k(cswritep, "\tALT");
    }
    if (provref_col) {
      cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
    }
    cswritep = strcpya_k(cswritep, "\tA1");
    if (omitted_col) {
      cswritep = strcpya_k(cswritep, "\tOMITTED");
    }
    if (ax_col) {
      cswritep = strcpya_k(cswritep, "\tAX");
    }
    if (a1_ct_col) {
      cswritep = strcpya_k(cswritep, "\tA1_CT");
    }
    if (tot_allele_col) {
      cswritep = strcpya_k(cswritep, "\tALLELE_CT");
    }
    if (a1_ct_cc_col) {
      cswritep = strcpya_k(cswritep, "\tA1_CASE_CT\tA1_CTRL_CT");
    }
    if (tot_allele_cc_col) {
      cswritep = strcpya_k(cswritep, "\tCASE_ALLELE_CT\tCTRL_ALLELE_CT");
    }
    if (gcount_cc_col) {
      cswritep = strcpya_k(cswritep, "\tCASE_NON_A1_CT\tCASE_HET_A1_CT\tCASE_HOM_A1_CT\tCTRL_NON_A1_CT\tCTRL_HET_A1_CT\tCTRL_HOM_A1_CT");
    }
    if (a1_freq_col) {
      cswritep = strcpya_k(cswritep, "\tA1_FREQ");
    }
    if (a1_freq_cc_col) {
      cswritep = strcpya_k(cswritep, "\tA1_CASE_FREQ\tA1_CTRL_FREQ");
    }
    if (mach_r2_col) {
      cswritep = strcpya_k(cswritep, "\tMACH_R2");
    }
    if (firth_yn_col) {
      cswritep = strcpya_k(cswritep, "\tFIRTH?");
    }
    if (test_col) {
      cswritep = strcpya_k(cswritep, "\tTEST");
    }
    if (nobs_col) {
      cswritep = strcpya_k(cswritep, "\tOBS_CT");
    }
    if (orbeta_col) {
      if (report_beta_instead_of_odds_ratio) {
        cswritep = strcpya_k(cswritep, "\tBETA");
      } else {
        cswritep = strcpya_k(cswritep, "\tOR");
      }
    }
    if (se_col) {
      if (report_beta_instead_of_odds_ratio) {
        cswritep = strcpya_k(cswritep, "\tSE");
      } else {
        cswritep = strcpya_k(cswritep, "\tLOG(OR)_SE");
      }
    }
    double ci_zt = 0.0;
    if (ci_col) {
      cswritep = strcpya_k(cswritep, "\tL");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      cswritep = strcpya_k(cswritep, "\tU");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      ci_zt = QuantileToZscore((ci_size + 1.0) * 0.5);
    }
    if (z_col) {
      if (!constraint_ct) {
        cswritep = strcpya_k(cswritep, "\tZ_STAT");
      } else {
        // F-statistic for joint tests.
        cswritep = strcpya_k(cswritep, "\tZ_OR_F_STAT");
      }
    }
    if (p_col) {
      if (report_neglog10p) {
        cswritep = strcpya_k(cswritep, "\tNEG_LOG10_P");
      } else {
        cswritep = strcpya_k(cswritep, "\tP");
      }
    }
    if (err_col) {
      cswritep = strcpya_k(cswritep, "\tERRCODE");
    }
    AppendBinaryEoln(&cswritep);

    // Main workflow:
    // 1. Set n=0, load/skip block 0
    //
    // 2. Spawn threads processing block n
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Load/skip block n unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8, Write results for last block
    uintptr_t write_variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t suppress_mach_r2 = 0;

    uint32_t cur_biallelic_reported_test_ct = 0;
    uint32_t primary_reported_test_idx = include_intercept;
    uint32_t cur_constraint_ct = 0;

    const char* const* cur_test_names = nullptr;
    uint32_t prev_block_variant_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    uint32_t allele_ct = 2;
    uint32_t omitted_allele_idx = 0;
    uintptr_t valid_allele_ct = 0;
    logprintfww5("--glm %s regression on phenotype '%s': ", is_always_firth? "Firth" : (is_sometimes_firth? "logistic-Firth hybrid" : "logistic"), cur_pheno_name);
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_variant_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto GlmLogistic_ret_PGR_FAIL;
      }
      if (local_covar_ct && cur_block_variant_ct) {
        const uint32_t uidx_start = read_block_idx * read_block_size;
        const uint32_t uidx_end = MINV(raw_variant_ct, uidx_start + read_block_size);
        if (local_variant_include) {
          reterr = ReadLocalCovarBlock(common, local_sample_uidx_order, local_variant_include, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, local_covar_txsp, &local_line_idx, &local_xy, is_single_prec? ctx->local_covars_vcmaj_f[parity] : nullptr, is_single_prec? nullptr : ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
        } else {
          float* prev_local_covar_row_f = nullptr;
          if (variant_idx) {
            prev_local_covar_row_f = &(ctx->local_covars_vcmaj_f[1 - parity][S_CAST(uintptr_t, read_block_size - 1) * max_sample_ct * local_covar_ct]);
          }
          reterr = ReadRfmix2Block(common, variant_bps, local_sample_uidx_order, prev_local_covar_row_f, nullptr, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, glm_info_ptr->local_chrom_col, glm_info_ptr->local_bp_col, glm_info_ptr->local_first_covar_col, local_covar_txsp, &local_line_iter, &local_line_idx, &local_prev_chr_code, &local_chr_code, &local_bp, &local_skip_chr, is_single_prec? ctx->local_covars_vcmaj_f[parity] : nullptr, is_single_prec? nullptr : ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
        }
        if (unlikely(reterr)) {
          goto GlmLogistic_ret_1;
        }
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, common->err_info);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, common->err_info >> 32);
          goto GlmLogistic_ret_1;
        }
      }
      if (!IsLastBlock(&tg)) {
        common->cur_block_variant_ct = cur_block_variant_ct;
        const uint32_t uidx_start = read_block_idx * read_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_variant_ct, calc_thread_ct, uidx_start, common->read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, common->pgr_ptrs);
        ctx->block_aux = logistic_block_aux_bufs[parity];
        common->block_beta_se = block_beta_se_bufs[parity];
        if (variant_idx + cur_block_variant_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto GlmLogistic_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const double* beta_se_iter = block_beta_se_bufs[parity];
        const LogisticAuxResult* cur_block_aux = logistic_block_aux_bufs[parity];
        uintptr_t allele_bidx = 0;
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            if ((chr_idx == x_code) && sample_ct_x) {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct_x;
              cur_constraint_ct = constraint_ct_x;
              cur_test_names = test_names_x;
            } else if ((chr_idx == y_code) && sample_ct_y) {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct_y;
              cur_constraint_ct = constraint_ct_y;
              cur_test_names = test_names_y;
            } else {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct;
              cur_constraint_ct = constraint_ct;
              cur_test_names = test_names;
            }
            suppress_mach_r2 = (chr_idx == x_code) || (chr_idx == mt_code);
            if (cur_constraint_ct) {
              // bugfix (17 May 2018): this was using reported_test_ct instead
              // of cur_reported_test_ct.
              primary_reported_test_idx = cur_biallelic_reported_test_ct - 1;
            }
            if (chr_col) {
              char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
              *chr_name_end = '\t';
              chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
            }
          }
          uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offsets[write_variant_uidx];
          }
          const uint32_t allele_ct_m1 = allele_ct - 1;
          const uint32_t extra_allele_ct = allele_ct - 2;
          if (omitted_alleles) {
            omitted_allele_idx = omitted_alleles[write_variant_uidx];
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          uint32_t variant_is_valid = 0;
          uint32_t a1_allele_idx = 0;
          for (uint32_t nonomitted_allele_idx = 0; nonomitted_allele_idx != allele_ct_m1; ++nonomitted_allele_idx, ++a1_allele_idx) {
            if (beta_se_multiallelic_fused) {
              if (!nonomitted_allele_idx) {
                primary_reported_test_idx = include_intercept;
              } else {
                primary_reported_test_idx = cur_biallelic_reported_test_ct + nonomitted_allele_idx - 1;
              }
            }
            if (nonomitted_allele_idx == omitted_allele_idx) {
              ++a1_allele_idx;
            }
            const double primary_beta = beta_se_iter[primary_reported_test_idx * 2];
            const double primary_se = beta_se_iter[primary_reported_test_idx * 2 + 1];
            const uint32_t allele_is_valid = (primary_se != -9.0);
            variant_is_valid |= allele_is_valid;
            {
              const LogisticAuxResult* auxp = &(cur_block_aux[allele_bidx]);
              if (ln_pfilter <= 0.0) {
                if (!allele_is_valid) {
                  goto GlmLogistic_allele_iterate;
                }
                double permstat;
                double primary_ln_pval;
                if (!cur_constraint_ct) {
                  permstat = fabs(primary_beta / primary_se);
                  // could precompute a tstat threshold instead
                  primary_ln_pval = ZscoreToLnP(permstat);
                } else {
                  // cur_constraint_ct may be different on chrX/chrY than it is
                  // on autosomes, so just have permstat be -log(pval) to be
                  // safe
                  primary_ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                  permstat = -primary_ln_pval;
                }
                if (primary_ln_pval > ln_pfilter) {
                  if (orig_ln_pvals) {
                    orig_ln_pvals[valid_allele_ct] = primary_ln_pval;
                  }
                  if (orig_permstat) {
                    orig_permstat[valid_allele_ct] = permstat;
                  }
                  goto GlmLogistic_allele_iterate;
                }
              }
              uint32_t inner_reported_test_ct = cur_biallelic_reported_test_ct;
              if (extra_allele_ct) {
                if (beta_se_multiallelic_fused) {
                  // in fused case, we're only performing a single multiple
                  // regression, so list all additive results together,
                  // possibly with intercept before.
                  if (!nonomitted_allele_idx) {
                    inner_reported_test_ct = 1 + include_intercept;
                  } else if (nonomitted_allele_idx == extra_allele_ct) {
                    inner_reported_test_ct -= include_intercept;
                  } else {
                    inner_reported_test_ct = 1;
                  }
                } else if (!hide_covar) {
                  inner_reported_test_ct += extra_allele_ct;
                }
              }
              // possible todo: make number-to-string operations, strlen(),
              // etc. happen only once per variant.
              for (uint32_t allele_test_idx = 0; allele_test_idx != inner_reported_test_ct; ++allele_test_idx) {
                uint32_t test_idx = allele_test_idx;
                if (beta_se_multiallelic_fused && nonomitted_allele_idx) {
                  if (!allele_test_idx) {
                    test_idx = primary_reported_test_idx;
                  } else {
                    // bugfix (26 Jun 2019): only correct to add 1 here in
                    // include_intercept case
                    test_idx += include_intercept;
                  }
                }
                if (chr_col) {
                  cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
                }
                if (variant_bps) {
                  cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
                }
                cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
                if (ref_col) {
                  *cswritep++ = '\t';
                  cswritep = strcpya(cswritep, cur_alleles[0]);
                }
                if (alt1_col) {
                  *cswritep++ = '\t';
                  cswritep = strcpya(cswritep, cur_alleles[1]);
                }
                if (alt_col) {
                  *cswritep++ = '\t';
                  for (uint32_t tmp_allele_idx = 1; tmp_allele_idx != allele_ct; ++tmp_allele_idx) {
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLogistic_ret_WRITE_FAIL;
                    }
                    cswritep = strcpyax(cswritep, cur_alleles[tmp_allele_idx], ',');
                  }
                  --cswritep;
                }
                *cswritep++ = '\t';
                if (provref_col) {
                  *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, write_variant_uidx)))? 'Y' : 'N';
                  *cswritep++ = '\t';
                }
                const uint32_t multi_a1 = extra_allele_ct && beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx);
                if (multi_a1) {
                  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    if (allele_idx == omitted_allele_idx) {
                      continue;
                    }
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLogistic_ret_WRITE_FAIL;
                    }
                    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                  }
                  --cswritep;
                } else {
                  cswritep = strcpya(cswritep, cur_alleles[a1_allele_idx]);
                }
                if (omitted_col) {
                  *cswritep++ = '\t';
                  cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
                }
                if (ax_col) {
                  *cswritep++ = '\t';
                  if (beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx)) {
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLogistic_ret_WRITE_FAIL;
                    }
                    cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
                  } else {
                    for (uint32_t tmp_allele_idx = 0; tmp_allele_idx != allele_ct; ++tmp_allele_idx) {
                      if (tmp_allele_idx == a1_allele_idx) {
                        continue;
                      }
                      if (unlikely(Cswrite(&css, &cswritep))) {
                        goto GlmLogistic_ret_WRITE_FAIL;
                      }
                      cswritep = strcpyax(cswritep, cur_alleles[tmp_allele_idx], ',');
                    }
                    --cswritep;
                  }
                }
                if (a1_ct_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_dosage, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (tot_allele_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa(auxp->allele_obs_ct, cswritep);
                }
                if (a1_ct_cc_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_case_dosage, cswritep);
                    *cswritep++ = '\t';
                    cswritep = dtoa_g(auxp->a1_dosage - auxp->a1_case_dosage, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA\tNA");
                  }
                }
                if (tot_allele_cc_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa_x(auxp->case_allele_obs_ct, '\t', cswritep);
                  cswritep = u32toa(auxp->allele_obs_ct - auxp->case_allele_obs_ct, cswritep);
                }
                if (gcount_cc_col) {
                  if (!multi_a1) {
                    STD_ARRAY_KREF(uint32_t, 6) cur_geno_hardcall_cts = auxp->geno_hardcall_cts;
                    for (uint32_t uii = 0; uii != 6; ++uii) {
                      *cswritep++ = '\t';
                      cswritep = u32toa(cur_geno_hardcall_cts[uii], cswritep);
                    }
                  } else {
                    cswritep = strcpya_k(cswritep, "\tNA\tNA\tNA\tNA\tNA\tNA");
                  }
                }
                if (a1_freq_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_dosage / S_CAST(double, auxp->allele_obs_ct), cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (a1_freq_cc_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_case_dosage / S_CAST(double, auxp->case_allele_obs_ct), cswritep);
                    *cswritep++ = '\t';
                    cswritep = dtoa_g((auxp->a1_dosage - auxp->a1_case_dosage) / S_CAST(double, auxp->allele_obs_ct - auxp->case_allele_obs_ct), cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA\tNA");
                  }
                }
                if (mach_r2_col) {
                  *cswritep++ = '\t';
                  if (!suppress_mach_r2) {
                    cswritep = dtoa_g(auxp->mach_r2, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (firth_yn_col) {
                  *cswritep++ = '\t';
                  // 'Y' - 'N' = 11
                  *cswritep++ = 'N' + 11 * auxp->firth_fallback;
                }
                if (test_col) {
                  *cswritep++ = '\t';
                  if (test_idx < cur_biallelic_reported_test_ct) {
                    cswritep = strcpya(cswritep, cur_test_names[test_idx]);
                  } else {
                    // always use basic dosage for untested alleles
                    cswritep = strcpya_k(cswritep, "ADD");
                    if (!beta_se_multiallelic_fused) {
                      // extra alt allele covariate.
                      uint32_t test_xallele_idx = test_idx - cur_biallelic_reported_test_ct;
                      // now we have the 0-based relative position in a list
                      // with the omitted_allele_idx and a1_allele_idx removed.
                      // correct this to the absolute index.  (there may be a
                      // cleaner way to do this with nonomitted_allele_idx?)
                      if (omitted_allele_idx < a1_allele_idx) {
                        test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                      }
                      test_xallele_idx = test_xallele_idx + (test_xallele_idx >= a1_allele_idx);
                      if (a1_allele_idx < omitted_allele_idx) {
                        test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                      }
                      if (!test_xallele_idx) {
                        cswritep = strcpya_k(cswritep, "_REF");
                      } else {
                        cswritep = strcpya_k(cswritep, "_ALT");
                        cswritep = u32toa(test_xallele_idx, cswritep);
                      }
                    }
                  }
                }
                if (nobs_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa(auxp->sample_obs_ct, cswritep);
                }
                double ln_pval = kLnPvalError;
                double permstat = 0.0;
                uint32_t test_is_valid;
                if ((!cur_constraint_ct) || (test_idx != primary_reported_test_idx)) {
                  double beta = beta_se_iter[2 * test_idx];
                  double se = beta_se_iter[2 * test_idx + 1];
                  test_is_valid = (se != -9.0);
                  if (test_is_valid) {
                    permstat = beta / se;
                    ln_pval = ZscoreToLnP(permstat);
                  }
                  if (orbeta_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      if (report_beta_instead_of_odds_ratio) {
                        cswritep = dtoa_g(beta, cswritep);
                      } else {
                        cswritep = lntoa_g(beta, cswritep);
                      }
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  if (se_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(se, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  if (ci_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      const double ci_halfwidth = ci_zt * se;
                      if (report_beta_instead_of_odds_ratio) {
                        cswritep = dtoa_g(beta - ci_halfwidth, cswritep);
                        *cswritep++ = '\t';
                        cswritep = dtoa_g(beta + ci_halfwidth, cswritep);
                      } else {
                        cswritep = lntoa_g(beta - ci_halfwidth, cswritep);
                        *cswritep++ = '\t';
                        cswritep = lntoa_g(beta + ci_halfwidth, cswritep);
                      }
                    } else {
                      cswritep = strcpya_k(cswritep, "NA\tNA");
                    }
                  }
                  if (z_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(permstat, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                } else {
                  // joint test: use F-test instead of Wald test
                  test_is_valid = allele_is_valid;
                  if (orbeta_col) {
                    cswritep = strcpya_k(cswritep, "\tNA");
                  }
                  if (se_col) {
                    cswritep = strcpya_k(cswritep, "\tNA");
                  }
                  if (ci_col) {
                    cswritep = strcpya_k(cswritep, "\tNA\tNA");
                  }
                  if (z_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(primary_se / u31tod(cur_constraint_ct), cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  // could avoid recomputing
                  if (test_is_valid) {
                    ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                    permstat = -ln_pval;
                  }
                }
                if (p_col) {
                  *cswritep++ = '\t';
                  if (test_is_valid) {
                    if (report_neglog10p) {
                      double reported_val = (-kRecipLn10) * ln_pval;
                      cswritep = dtoa_g(reported_val, cswritep);
                    } else {
                      double reported_ln = MAXV(ln_pval, output_min_ln);
                      cswritep = lntoa_g(reported_ln, cswritep);
                    }
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (err_col) {
                  *cswritep++ = '\t';
                  if (test_is_valid) {
                    if (!auxp->is_unfinished) {
                      *cswritep++ = '.';
                    } else {
                      cswritep = strcpya_k(cswritep, "UNFINISHED");
                    }
                  } else {
                    uint64_t glm_errcode;
                    memcpy(&glm_errcode, &(beta_se_iter[2 * test_idx]), 8);
                    cswritep = AppendGlmErrstr(glm_errcode, cswritep);
                  }
                }
                AppendBinaryEoln(&cswritep);
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto GlmLogistic_ret_WRITE_FAIL;
                }
                if ((test_idx == primary_reported_test_idx) && allele_is_valid) {
                  if (orig_ln_pvals) {
                    orig_ln_pvals[valid_allele_ct] = ln_pval;
                  }
                  if (orig_permstat) {
                    orig_permstat[valid_allele_ct] = permstat;
                  }
                }
              }
            }
          GlmLogistic_allele_iterate:
            ++allele_bidx;
            valid_allele_ct += allele_is_valid;
            if (valid_alleles && allele_is_valid) {
              SetBit(allele_idx_offset_base + a1_allele_idx, valid_alleles);
            }
            if (!beta_se_multiallelic_fused) {
              beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
            }
          }
          if (beta_se_multiallelic_fused) {
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
          if ((!variant_is_valid) && valid_alleles) {
            ClearBit(write_variant_uidx, valid_variants);
          }
        }
      }
      if (variant_idx == variant_ct) {
        break;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
      }
      ++read_block_idx;
      prev_block_variant_ct = cur_block_variant_ct;
      variant_idx += cur_block_variant_ct;
      // crucially, this is independent of the PgenReader block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto GlmLogistic_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    logprintf("Results written to %s .\n", outname);
#ifndef NDEBUG
    // temporary debug
    if (g_debug_on && (variant_ct == 1) && ctx->separation_found) {
      char fname_buf[kPglFnamesize + 8];
      const uint32_t outname_slen = strlen(outname);
      char* fname_write_iter = memcpya(fname_buf, outname, outname_slen);
      strcpy_k(fname_write_iter, ".inputs");
      logprintf("Debug matrix written to %s .\n", fname_buf);
    }
#endif
    *valid_allele_ct_ptr = valid_allele_ct;
  }
  while (0) {
  GlmLogistic_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GlmLogistic_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
    break;
  GlmLogistic_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  GlmLogistic_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GlmLogistic_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  GlmLogistic_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 GlmLogistic_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

/*
static inline BoolErr vecaligned_malloc2(uintptr_t size, void* aligned_pp) {
  if (vecaligned_malloc(size, aligned_pp)) {
    return 1;
  }
  memset(*S_CAST(unsigned char**, aligned_pp), 85, size);
  return 0;
}

void LogisticTestFirthF() {
  const uint32_t sample_ct = 2777;
  const uint32_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  const uint32_t predictor_ct = 5;
  const uint32_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  float* yy;
  float* xx;
  float* beta;
  float* pp;
  float* vv;
  float* hh;
  float* ustar;
  float* delta;
  MatrixInvertBuf1* inv_1d_buf;
  double* dbl_2d_buf;
  float* hdiag;
  float* ww;
  float* hh0_buf;
  float* tmpnxk_buf;
  double* half_inverted_buf;
  vecaligned_malloc2(sample_ctav * sizeof(float), &yy);
  vecaligned_malloc2(sample_ctav * predictor_ct * sizeof(float), &xx);
  vecaligned_malloc2(predictor_ctav * sizeof(float), &beta);
  vecaligned_malloc2(sample_ctav * sizeof(float), &pp);
  vecaligned_malloc2(sample_ctav * sizeof(float), &vv);
  vecaligned_malloc2(predictor_ct * predictor_ctav * sizeof(float), &hh);
  vecaligned_malloc2(predictor_ctav * sizeof(float), &ustar);
  vecaligned_malloc2(predictor_ctav * sizeof(float), &delta);
  vecaligned_malloc2(predictor_ct * kMatrixInvertBuf1CheckedAlloc, &inv_1d_buf);
  vecaligned_malloc2(predictor_ct * MAXV(predictor_ct, 7) * sizeof(double), &dbl_2d_buf);
  vecaligned_malloc2(sample_ctav * sizeof(float), &hdiag);
  vecaligned_malloc2(sample_ctav * sizeof(float), &ww);
  vecaligned_malloc2(predictor_ct * predictor_ctav * sizeof(float), &hh0_buf);
  vecaligned_malloc2(predictor_ct * sample_ctav * sizeof(float), &tmpnxk_buf);
  vecaligned_malloc2(predictor_ct * MAXV(predictor_ct, 3) * sizeof(double), &half_inverted_buf);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    if (cos(S_CAST(double, sample_idx + 1) * 0.001) > 0.0) {
      yy[sample_idx] = 1.0;
    } else {
      yy[sample_idx] = 0.0;
    }
  }
  ZeroFArr(sample_ctav - sample_ct, &(yy[sample_ct]));
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    float* xx_col = &(xx[pred_idx * sample_ctav]);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      float dxx = S_CAST(float, 1.0);
      float dyy = S_CAST(float, sample_idx + 1) * S_CAST(float, 0.001);
      for (uint32_t uii = 0; uii < pred_idx; ++uii) {
        dxx *= dyy;
      }
      xx_col[sample_idx] = dxx;
    }
    ZeroFArr(sample_ctav - sample_ct, &(xx_col[sample_ct]));
  }
  ZeroFArr(predictor_ctav, beta);
  uint32_t is_unfinished = 0;
  BoolErr reterr = FirthRegressionF(yy, xx, nullptr, sample_ct, predictor_ct, beta, &is_unfinished, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf, pp, vv, ustar, delta, hdiag, ww, hh0_buf, tmpnxk_buf);
  printf("coef: %g %g %g %g %g\n", S_CAST(double, beta[0]), S_CAST(double, beta[1]), S_CAST(double, beta[2]), S_CAST(double, beta[3]), S_CAST(double, beta[4]));
  printf("sqrt(diag(hh)): %g %g %g %g %g\n", S_CAST(double, sqrtf(hh[0])), S_CAST(double, sqrtf(hh[predictor_ctav + 1])), S_CAST(double, sqrtf(hh[2 * (predictor_ctav + 1)])), S_CAST(double, sqrtf(hh[3 * (predictor_ctav + 1)])), S_CAST(double, sqrtf(hh[4 * (predictor_ctav + 1)])));
  printf("p-values:");
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    const double bt = S_CAST(double, beta[pred_idx]);
    const double se = S_CAST(double, sqrtf(hh[pred_idx * (predictor_ctav + 1)]));
    const double permstat = bt / se;
    const double ln_pval = ZscoreToLnP(permstat);
    char outbuf[24];
    char* outbuf_iter = outbuf;
    *outbuf_iter++ = ' ';
    outbuf_iter = lntoa_g(ln_pval, outbuf_iter);
    *outbuf_iter = '\0';
    fputs(outbuf, stdout);
  }
  printf("\n");
  printf("reterr: %u  is_unfinished: %u\n", S_CAST(uint32_t, reterr), is_unfinished);
}

void LogisticTestD() {
  const uint32_t sample_ct = 2777;
  const uint32_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uint32_t predictor_ct = 5;
  const uint32_t predictor_ctav = RoundUpPow2(predictor_ct, kDoublePerDVec);
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
  vecaligned_malloc2(sample_ctav * sizeof(double), &yy);
  vecaligned_malloc2(sample_ctav * predictor_ct * sizeof(double), &xx);
  vecaligned_malloc2(predictor_ctav * sizeof(double), &coef);
  vecaligned_malloc2(predictor_ct * predictor_ctav * sizeof(double), &ll);
  vecaligned_malloc2(sample_ctav * sizeof(double), &pp);
  vecaligned_malloc2(sample_ctav * sizeof(double), &vv);
  vecaligned_malloc2(predictor_ct * predictor_ctav * sizeof(double), &hh);
  vecaligned_malloc2(predictor_ctav * sizeof(double), &grad);
  vecaligned_malloc2(predictor_ctav * sizeof(double), &dcoef);
  vecaligned_malloc2(predictor_ct * kMatrixInvertBuf1CheckedAlloc, &inv_1d_buf);
  vecaligned_malloc2(predictor_ct * MAXV(predictor_ct, 7) * sizeof(double), &dbl_2d_buf);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    if (cos(S_CAST(double, sample_idx + 1) * 0.001) > 0.0) {
      yy[sample_idx] = 1.0;
    } else {
      yy[sample_idx] = 0.0;
    }
  }
  ZeroDArr(sample_ctav - sample_ct, &(yy[sample_ct]));
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    double* xx_col = &(xx[pred_idx * sample_ctav]);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      double dxx = 1.0;
      double dyy = S_CAST(double, sample_idx + 1) * 0.001;
      for (uint32_t uii = 0; uii < pred_idx; ++uii) {
        dxx *= dyy;
      }
      xx_col[sample_idx] = dxx;
    }
    ZeroDArr(sample_ctav - sample_ct, &(xx_col[sample_ct]));
  }
  ZeroDArr(predictor_ctav, coef);
  uint32_t is_unfinished = 0;
  BoolErr reterr = LogisticRegressionD(yy, xx, nullptr, sample_ct, predictor_ct, &is_unfinished, coef, ll, pp, vv, hh, grad, dcoef, inv_1d_buf, dbl_2d_buf);
  printf("coef: %g %g %g %g %g\n", coef[0], coef[1], coef[2], coef[3], coef[4]);
  for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
    double* hh_inv_row = &(hh[pred_uidx * predictor_ctav]);
    ZeroDArr(predictor_ct, grad);
    grad[pred_uidx] = 1.0;
    SolveLinearSystemD(ll, grad, predictor_ct, hh_inv_row);
  }
  printf("sqrt(diag(hh)): %g %g %g %g %g\n", sqrt(hh[0]), sqrt(hh[predictor_ctav + 1]), sqrt(hh[2 * (predictor_ctav + 1)]), sqrt(hh[3 * (predictor_ctav + 1)]), sqrt(hh[4 * (predictor_ctav + 1)]));
  printf("p-values:");
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    const double bt = coef[pred_idx];
    const double se = sqrt(hh[pred_idx * (predictor_ctav + 1)]);
    const double permstat = bt / se;
    const double ln_pval = ZscoreToLnP(permstat);
    char outbuf[24];
    char* outbuf_iter = outbuf;
    *outbuf_iter++ = ' ';
    outbuf_iter = lntoa_g(ln_pval, outbuf_iter);
    *outbuf_iter = '\0';
    fputs(outbuf, stdout);
  }
  printf("\n");
  printf("reterr: %u  is_unfinished: %u\n", S_CAST(uint32_t, reterr), is_unfinished);
}

void LogisticTestFirthD() {
  const uint32_t sample_ct = 2777;
  const uint32_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  const uint32_t predictor_ct = 5;
  const uint32_t predictor_ctav = RoundUpPow2(predictor_ct, kDoublePerDVec);
  double* yy;
  double* xx;
  double* beta;
  double* pp;
  double* vv;
  double* hh;
  double* ustar;
  double* delta;
  MatrixInvertBuf1* inv_1d_buf;
  double* dbl_2d_buf;
  double* hdiag;
  double* ww;
  double* hh0_buf;
  double* tmpnxk_buf;
  vecaligned_malloc2(sample_ctav * sizeof(double), &yy);
  vecaligned_malloc2(sample_ctav * predictor_ct * sizeof(double), &xx);
  vecaligned_malloc2(predictor_ctav * sizeof(double), &beta);
  vecaligned_malloc2(sample_ctav * sizeof(double), &pp);
  vecaligned_malloc2(sample_ctav * sizeof(double), &vv);
  vecaligned_malloc2(predictor_ct * predictor_ctav * sizeof(double), &hh);
  vecaligned_malloc2(predictor_ctav * sizeof(double), &ustar);
  vecaligned_malloc2(predictor_ctav * sizeof(double), &delta);
  vecaligned_malloc2(predictor_ct * kMatrixInvertBuf1CheckedAlloc, &inv_1d_buf);
  vecaligned_malloc2(predictor_ct * MAXV(predictor_ct, 7) * sizeof(double), &dbl_2d_buf);
  vecaligned_malloc2(sample_ctav * sizeof(double), &hdiag);
  vecaligned_malloc2(sample_ctav * sizeof(double), &ww);
  vecaligned_malloc2(predictor_ct * predictor_ctav * sizeof(double), &hh0_buf);
  vecaligned_malloc2(predictor_ct * sample_ctav * sizeof(double), &tmpnxk_buf);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    if (cos(S_CAST(double, sample_idx + 1) * 0.001) > 0.0) {
      yy[sample_idx] = 1.0;
    } else {
      yy[sample_idx] = 0.0;
    }
  }
  ZeroDArr(sample_ctav - sample_ct, &(yy[sample_ct]));
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    double* xx_col = &(xx[pred_idx * sample_ctav]);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      double dxx = 1.0;
      double dyy = S_CAST(double, sample_idx + 1) * 0.001;
      for (uint32_t uii = 0; uii < pred_idx; ++uii) {
        dxx *= dyy;
      }
      xx_col[sample_idx] = dxx;
    }
    ZeroDArr(sample_ctav - sample_ct, &(xx_col[sample_ct]));
  }
  ZeroDArr(predictor_ctav, beta);
  uint32_t is_unfinished = 0;
  BoolErr reterr = FirthRegressionD(yy, xx, nullptr, sample_ct, predictor_ct, &is_unfinished, beta, hh, inv_1d_buf, dbl_2d_buf, pp, vv, ustar, delta, hdiag, ww, hh0_buf, tmpnxk_buf);
  printf("beta: %g %g %g %g %g\n", beta[0], beta[1], beta[2], beta[3], beta[4]);
  printf("sqrt(diag(hh)): %g %g %g %g %g\n", sqrt(hh[0]), sqrt(hh[predictor_ctav + 1]), sqrt(hh[2 * (predictor_ctav + 1)]), sqrt(hh[3 * (predictor_ctav + 1)]), sqrt(hh[4 * (predictor_ctav + 1)]));
  printf("p-values:");
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    const double bt = beta[pred_idx];
    const double se = sqrt(hh[pred_idx * (predictor_ctav + 1)]);
    const double permstat = bt / se;
    const double ln_pval = ZscoreToLnP(permstat);
    char outbuf[24];
    char* outbuf_iter = outbuf;
    *outbuf_iter++ = ' ';
    outbuf_iter = lntoa_g(ln_pval, outbuf_iter);
    *outbuf_iter = '\0';
    fputs(outbuf, stdout);
  }
  printf("\n");
  printf("reterr: %u  is_unfinished: %u\n", S_CAST(uint32_t, reterr), is_unfinished);
}

void LogisticTestInternal() {
  LogisticTestFirthD();
}
*/

#ifdef __cplusplus
}  // namespace plink2
#endif
