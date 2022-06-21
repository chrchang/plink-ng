// This file is part of PLINK 2.00, copyright (C) 2005-2022 Shaun Purcell,
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

#include "include/plink2_stats.h"
#include "plink2_compress_stream.h"
#include "plink2_glm_logistic.h"

#ifdef __LP64__
#  ifdef __x86_64__
#    include <emmintrin.h>
#  else
#    define SIMDE_ENABLE_NATIVE_ALIASES
#    include "x86/sse2.h"
#  endif
#endif

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

// Only called by GlmLogisticThread(), so there are the following differences
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
    double* sample_corr_row = &(inverse_corr_buf[pred_idx1 * predictor_ct]);
    const float* predictor_dotprod_row = &(predictor_dotprod_buf[pred_idx1 * predictor_ct]);
    const double covar1_mean_adj = dbl_2d_buf[pred_idx1] * sample_ct_recip;
    for (uint32_t pred_idx2 = 0; pred_idx2 <= pred_idx1; ++pred_idx2) {
      sample_corr_row[pred_idx2] = (S_CAST(double, predictor_dotprod_row[pred_idx2]) - covar1_mean_adj * dbl_2d_buf[pred_idx2]) * sample_ct_m1_recip;
    }
  }
  // now use dbl_2d_buf to store inverse-sqrts, to get to correlation matrix
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    dbl_2d_buf[pred_idx] = 1.0 / sqrt(inverse_corr_buf[pred_idx * predictor_ct_p1]);
  }
  // invert_symmdef_matrix only cares about bottom left of inverse_corr_buf[]
  for (uint32_t pred_idx1 = 1; pred_idx1 != predictor_ct; ++pred_idx1) {
    const double inverse_stdev1 = dbl_2d_buf[pred_idx1];
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
// The two instances of fmath_exp_ps() are C ports of Shigeo Mitsunari's fast
// math library functions posted at https://github.com/herumi/fmath .  License
// is http://opensource.org/licenses/BSD-3-Clause .
// (I tried porting fmath_log_ps, but it turns out that Firth regression needs
// double-precision log accuracy; logf() actually interferes with convergence.)

// programmatically generated by:
// typedef union {
//   float f4;
//   uint32_t u4;
// } __uni4;
//
// __uni4 u4;
// int32_t ii;
// for (ii = 0; ii < 1024; ii++) {
//   u4.f4 = pow(2.0f, ((float)ii) / 1024.0);
//   printf("0x%08x", u4.u4 & 0x7fffff);
//   if (ii % 4 != 3) {
//     printf(", ");
//   } else {
//     printf(",\n");
//   }
// }

const uint32_t kFloatExpLookupInt[]
#ifdef FVEC_32
  __attribute__((aligned(32)))
#else
  __attribute__((aligned(16)))
#endif
  = {
0x00000000, 0x00001630, 0x00002c64, 0x0000429c,
0x000058d8, 0x00006f17, 0x0000855b, 0x00009ba2,
0x0000b1ed, 0x0000c83c, 0x0000de8f, 0x0000f4e6,
0x00010b41, 0x0001219f, 0x00013802, 0x00014e68,
0x000164d2, 0x00017b40, 0x000191b2, 0x0001a828,
0x0001bea1, 0x0001d51f, 0x0001eba1, 0x00020226,
0x000218af, 0x00022f3c, 0x000245ce, 0x00025c63,
0x000272fc, 0x00028998, 0x0002a039, 0x0002b6de,
0x0002cd87, 0x0002e433, 0x0002fae4, 0x00031198,
0x00032850, 0x00033f0d, 0x000355cd, 0x00036c91,
0x00038359, 0x00039a25, 0x0003b0f5, 0x0003c7c9,
0x0003dea1, 0x0003f57d, 0x00040c5d, 0x00042341,
0x00043a29, 0x00045115, 0x00046804, 0x00047ef8,
0x000495f0, 0x0004aceb, 0x0004c3eb, 0x0004daef,
0x0004f1f6, 0x00050902, 0x00052012, 0x00053725,
0x00054e3d, 0x00056558, 0x00057c78, 0x0005939c,
0x0005aac3, 0x0005c1ef, 0x0005d91f, 0x0005f052,
0x0006078a, 0x00061ec6, 0x00063606, 0x00064d4a,
0x00066491, 0x00067bdd, 0x0006932d, 0x0006aa81,
0x0006c1d9, 0x0006d935, 0x0006f095, 0x000707f9,
0x00071f62, 0x000736ce, 0x00074e3e, 0x000765b3,
0x00077d2b, 0x000794a8, 0x0007ac28, 0x0007c3ad,
0x0007db35, 0x0007f2c2, 0x00080a53, 0x000821e8,
0x00083981, 0x0008511e, 0x000868c0, 0x00088065,
0x0008980f, 0x0008afbc, 0x0008c76e, 0x0008df23,
0x0008f6dd, 0x00090e9b, 0x0009265d, 0x00093e24,
0x000955ee, 0x00096dbc, 0x0009858f, 0x00099d66,
0x0009b541, 0x0009cd20, 0x0009e503, 0x0009fcea,
0x000a14d5, 0x000a2cc5, 0x000a44b9, 0x000a5cb1,
0x000a74ad, 0x000a8cad, 0x000aa4b1, 0x000abcba,
0x000ad4c6, 0x000aecd7, 0x000b04ec, 0x000b1d05,
0x000b3523, 0x000b4d44, 0x000b656a, 0x000b7d94,
0x000b95c2, 0x000badf4, 0x000bc62b, 0x000bde65,
0x000bf6a4, 0x000c0ee7, 0x000c272f, 0x000c3f7a,
0x000c57ca, 0x000c701e, 0x000c8876, 0x000ca0d2,
0x000cb933, 0x000cd198, 0x000cea01, 0x000d026e,
0x000d1adf, 0x000d3355, 0x000d4bcf, 0x000d644d,
0x000d7cd0, 0x000d9556, 0x000dade1, 0x000dc671,
0x000ddf04, 0x000df79c, 0x000e1038, 0x000e28d8,
0x000e417d, 0x000e5a25, 0x000e72d3, 0x000e8b84,
0x000ea43a, 0x000ebcf3, 0x000ed5b2, 0x000eee74,
0x000f073b, 0x000f2006, 0x000f38d5, 0x000f51a9,
0x000f6a81, 0x000f835d, 0x000f9c3e, 0x000fb523,
0x000fce0c, 0x000fe6fa, 0x000fffec, 0x001018e2,
0x001031dc, 0x00104adb, 0x001063de, 0x00107ce6,
0x001095f2, 0x0010af02, 0x0010c816, 0x0010e12f,
0x0010fa4d, 0x0011136e, 0x00112c94, 0x001145be,
0x00115eed, 0x00117820, 0x00119158, 0x0011aa93,
0x0011c3d3, 0x0011dd18, 0x0011f661, 0x00120fae,
0x00122900, 0x00124256, 0x00125bb0, 0x0012750f,
0x00128e72, 0x0012a7da, 0x0012c146, 0x0012dab7,
0x0012f42c, 0x00130da5, 0x00132723, 0x001340a5,
0x00135a2b, 0x001373b6, 0x00138d46, 0x0013a6d9,
0x0013c072, 0x0013da0e, 0x0013f3af, 0x00140d55,
0x001426ff, 0x001440ae, 0x00145a60, 0x00147418,
0x00148dd4, 0x0014a794, 0x0014c159, 0x0014db22,
0x0014f4f0, 0x00150ec2, 0x00152898, 0x00154274,
0x00155c53, 0x00157637, 0x00159020, 0x0015aa0d,
0x0015c3ff, 0x0015ddf5, 0x0015f7ef, 0x001611ee,
0x00162bf2, 0x001645fa, 0x00166006, 0x00167a18,
0x0016942d, 0x0016ae47, 0x0016c866, 0x0016e289,
0x0016fcb1, 0x001716dd, 0x0017310e, 0x00174b43,
0x0017657d, 0x00177fbc, 0x001799ff, 0x0017b446,
0x0017ce92, 0x0017e8e3, 0x00180338, 0x00181d92,
0x001837f0, 0x00185253, 0x00186cbb, 0x00188727,
0x0018a197, 0x0018bc0d, 0x0018d686, 0x0018f105,
0x00190b88, 0x0019260f, 0x0019409c, 0x00195b2c,
0x001975c2, 0x0019905c, 0x0019aafa, 0x0019c59e,
0x0019e046, 0x0019faf2, 0x001a15a3, 0x001a3059,
0x001a4b13, 0x001a65d2, 0x001a8096, 0x001a9b5e,
0x001ab62b, 0x001ad0fd, 0x001aebd3, 0x001b06ae,
0x001b218d, 0x001b3c71, 0x001b575a, 0x001b7248,
0x001b8d3a, 0x001ba831, 0x001bc32c, 0x001bde2c,
0x001bf931, 0x001c143b, 0x001c2f49, 0x001c4a5c,
0x001c6573, 0x001c8090, 0x001c9bb1, 0x001cb6d6,
0x001cd201, 0x001ced30, 0x001d0864, 0x001d239c,
0x001d3eda, 0x001d5a1c, 0x001d7562, 0x001d90ae,
0x001dabfe, 0x001dc753, 0x001de2ad, 0x001dfe0b,
0x001e196e, 0x001e34d6, 0x001e5043, 0x001e6bb4,
0x001e872a, 0x001ea2a5, 0x001ebe25, 0x001ed9a9,
0x001ef532, 0x001f10c0, 0x001f2c53, 0x001f47eb,
0x001f6387, 0x001f7f28, 0x001f9ace, 0x001fb679,
0x001fd228, 0x001feddc, 0x00200996, 0x00202553,
0x00204116, 0x00205cde, 0x002078aa, 0x0020947b,
0x0020b051, 0x0020cc2c, 0x0020e80b, 0x002103f0,
0x00211fd9, 0x00213bc7, 0x002157ba, 0x002173b2,
0x00218faf, 0x0021abb0, 0x0021c7b7, 0x0021e3c2,
0x0021ffd2, 0x00221be7, 0x00223801, 0x0022541f,
0x00227043, 0x00228c6b, 0x0022a899, 0x0022c4cb,
0x0022e102, 0x0022fd3e, 0x0023197f, 0x002335c5,
0x0023520f, 0x00236e5f, 0x00238ab3, 0x0023a70d,
0x0023c36b, 0x0023dfce, 0x0023fc37, 0x002418a4,
0x00243516, 0x0024518d, 0x00246e08, 0x00248a89,
0x0024a70f, 0x0024c39a, 0x0024e029, 0x0024fcbe,
0x00251958, 0x002535f6, 0x00255299, 0x00256f42,
0x00258bef, 0x0025a8a2, 0x0025c559, 0x0025e215,
0x0025fed7, 0x00261b9d, 0x00263868, 0x00265538,
0x0026720e, 0x00268ee8, 0x0026abc7, 0x0026c8ac,
0x0026e595, 0x00270283, 0x00271f76, 0x00273c6f,
0x0027596c, 0x0027766e, 0x00279376, 0x0027b082,
0x0027cd94, 0x0027eaaa, 0x002807c6, 0x002824e6,
0x0028420c, 0x00285f37, 0x00287c66, 0x0028999b,
0x0028b6d5, 0x0028d414, 0x0028f158, 0x00290ea1,
0x00292bef, 0x00294942, 0x0029669b, 0x002983f8,
0x0029a15b, 0x0029bec2, 0x0029dc2f, 0x0029f9a1,
0x002a1718, 0x002a3494, 0x002a5215, 0x002a6f9b,
0x002a8d26, 0x002aaab7, 0x002ac84c, 0x002ae5e7,
0x002b0387, 0x002b212c, 0x002b3ed6, 0x002b5c85,
0x002b7a3a, 0x002b97f3, 0x002bb5b2, 0x002bd376,
0x002bf13f, 0x002c0f0d, 0x002c2ce0, 0x002c4ab9,
0x002c6897, 0x002c867a, 0x002ca462, 0x002cc24f,
0x002ce041, 0x002cfe39, 0x002d1c36, 0x002d3a38,
0x002d583f, 0x002d764b, 0x002d945d, 0x002db274,
0x002dd090, 0x002deeb1, 0x002e0cd8, 0x002e2b03,
0x002e4934, 0x002e676b, 0x002e85a6, 0x002ea3e7,
0x002ec22d, 0x002ee078, 0x002efec8, 0x002f1d1e,
0x002f3b79, 0x002f59d9, 0x002f783e, 0x002f96a9,
0x002fb519, 0x002fd38e, 0x002ff209, 0x00301089,
0x00302f0e, 0x00304d98, 0x00306c28, 0x00308abd,
0x0030a957, 0x0030c7f7, 0x0030e69c, 0x00310546,
0x003123f6, 0x003142aa, 0x00316165, 0x00318024,
0x00319ee9, 0x0031bdb3, 0x0031dc83, 0x0031fb57,
0x00321a32, 0x00323911, 0x003257f6, 0x003276e0,
0x003295d0, 0x0032b4c5, 0x0032d3bf, 0x0032f2bf,
0x003311c4, 0x003330cf, 0x00334fde, 0x00336ef4,
0x00338e0e, 0x0033ad2e, 0x0033cc54, 0x0033eb7e,
0x00340aaf, 0x003429e4, 0x0034491f, 0x00346860,
0x003487a6, 0x0034a6f1, 0x0034c642, 0x0034e598,
0x003504f3, 0x00352454, 0x003543bb, 0x00356327,
0x00358298, 0x0035a20f, 0x0035c18b, 0x0035e10d,
0x00360094, 0x00362020, 0x00363fb2, 0x00365f4a,
0x00367ee7, 0x00369e89, 0x0036be31, 0x0036dddf,
0x0036fd92, 0x00371d4a, 0x00373d08, 0x00375ccc,
0x00377c95, 0x00379c63, 0x0037bc37, 0x0037dc11,
0x0037fbf0, 0x00381bd4, 0x00383bbe, 0x00385bae,
0x00387ba3, 0x00389b9e, 0x0038bb9e, 0x0038dba4,
0x0038fbaf, 0x00391bc0, 0x00393bd7, 0x00395bf3,
0x00397c14, 0x00399c3b, 0x0039bc68, 0x0039dc9a,
0x0039fcd2, 0x003a1d10, 0x003a3d53, 0x003a5d9b,
0x003a7dea, 0x003a9e3e, 0x003abe97, 0x003adef6,
0x003aff5b, 0x003b1fc5, 0x003b4035, 0x003b60aa,
0x003b8126, 0x003ba1a6, 0x003bc22d, 0x003be2b9,
0x003c034a, 0x003c23e2, 0x003c447f, 0x003c6521,
0x003c85ca, 0x003ca678, 0x003cc72b, 0x003ce7e5,
0x003d08a4, 0x003d2968, 0x003d4a33, 0x003d6b03,
0x003d8bd8, 0x003dacb4, 0x003dcd95, 0x003dee7c,
0x003e0f68, 0x003e305a, 0x003e5152, 0x003e7250,
0x003e9353, 0x003eb45c, 0x003ed56b, 0x003ef67f,
0x003f179a, 0x003f38ba, 0x003f59df, 0x003f7b0b,
0x003f9c3c, 0x003fbd73, 0x003fdeb0, 0x003ffff2,
0x0040213b, 0x00404289, 0x004063dc, 0x00408536,
0x0040a695, 0x0040c7fb, 0x0040e966, 0x00410ad6,
0x00412c4d, 0x00414dc9, 0x00416f4b, 0x004190d3,
0x0041b261, 0x0041d3f5, 0x0041f58e, 0x0042172d,
0x004238d2, 0x00425a7d, 0x00427c2e, 0x00429de4,
0x0042bfa1, 0x0042e163, 0x0043032b, 0x004324f9,
0x004346cd, 0x004368a7, 0x00438a86, 0x0043ac6b,
0x0043ce57, 0x0043f048, 0x0044123f, 0x0044343c,
0x0044563f, 0x00447848, 0x00449a56, 0x0044bc6b,
0x0044de85, 0x004500a5, 0x004522cc, 0x004544f8,
0x0045672a, 0x00458962, 0x0045aba0, 0x0045cde4,
0x0045f02e, 0x0046127e, 0x004634d3, 0x0046572f,
0x00467991, 0x00469bf8, 0x0046be66, 0x0046e0d9,
0x00470353, 0x004725d2, 0x00474858, 0x00476ae3,
0x00478d75, 0x0047b00c, 0x0047d2aa, 0x0047f54d,
0x004817f7, 0x00483aa6, 0x00485d5b, 0x00488017,
0x0048a2d8, 0x0048c5a0, 0x0048e86d, 0x00490b41,
0x00492e1b, 0x004950fa, 0x004973e0, 0x004996cc,
0x0049b9be, 0x0049dcb5, 0x0049ffb3, 0x004a22b7,
0x004a45c1, 0x004a68d1, 0x004a8be8, 0x004aaf04,
0x004ad226, 0x004af54f, 0x004b187d, 0x004b3bb2,
0x004b5eed, 0x004b822e, 0x004ba575, 0x004bc8c2,
0x004bec15, 0x004c0f6e, 0x004c32ce, 0x004c5633,
0x004c799f, 0x004c9d11, 0x004cc089, 0x004ce407,
0x004d078c, 0x004d2b16, 0x004d4ea7, 0x004d723d,
0x004d95da, 0x004db97e, 0x004ddd27, 0x004e00d6,
0x004e248c, 0x004e4848, 0x004e6c0a, 0x004e8fd2,
0x004eb3a1, 0x004ed775, 0x004efb50, 0x004f1f31,
0x004f4319, 0x004f6706, 0x004f8afa, 0x004faef4,
0x004fd2f4, 0x004ff6fb, 0x00501b08, 0x00503f1b,
0x00506334, 0x00508753, 0x0050ab79, 0x0050cfa5,
0x0050f3d7, 0x00511810, 0x00513c4f, 0x00516094,
0x005184df, 0x0051a931, 0x0051cd89, 0x0051f1e7,
0x0052164c, 0x00523ab7, 0x00525f28, 0x005283a0,
0x0052a81e, 0x0052cca2, 0x0052f12c, 0x005315bd,
0x00533a54, 0x00535ef2, 0x00538396, 0x0053a840,
0x0053ccf1, 0x0053f1a8, 0x00541665, 0x00543b29,
0x00545ff3, 0x005484c3, 0x0054a99a, 0x0054ce77,
0x0054f35b, 0x00551845, 0x00553d35, 0x0055622c,
0x00558729, 0x0055ac2d, 0x0055d137, 0x0055f647,
0x00561b5e, 0x0056407b, 0x0056659f, 0x00568ac9,
0x0056affa, 0x0056d531, 0x0056fa6e, 0x00571fb2,
0x005744fd, 0x00576a4e, 0x00578fa5, 0x0057b503,
0x0057da67, 0x0057ffd2, 0x00582543, 0x00584abb,
0x00587039, 0x005895be, 0x0058bb49, 0x0058e0db,
0x00590673, 0x00592c12, 0x005951b8, 0x00597763,
0x00599d16, 0x0059c2cf, 0x0059e88e, 0x005a0e54,
0x005a3421, 0x005a59f4, 0x005a7fcd, 0x005aa5ae,
0x005acb94, 0x005af182, 0x005b1776, 0x005b3d70,
0x005b6371, 0x005b8979, 0x005baf87, 0x005bd59c,
0x005bfbb8, 0x005c21da, 0x005c4802, 0x005c6e32,
0x005c9468, 0x005cbaa4, 0x005ce0e7, 0x005d0731,
0x005d2d82, 0x005d53d9, 0x005d7a36, 0x005da09b,
0x005dc706, 0x005ded77, 0x005e13f0, 0x005e3a6f,
0x005e60f5, 0x005e8781, 0x005eae14, 0x005ed4ae,
0x005efb4e, 0x005f21f5, 0x005f48a3, 0x005f6f58,
0x005f9613, 0x005fbcd5, 0x005fe39e, 0x00600a6d,
0x00603143, 0x00605820, 0x00607f03, 0x0060a5ee,
0x0060ccdf, 0x0060f3d7, 0x00611ad5, 0x006141db,
0x006168e7, 0x00618ffa, 0x0061b713, 0x0061de34,
0x0062055b, 0x00622c89, 0x006253be, 0x00627af9,
0x0062a23c, 0x0062c985, 0x0062f0d5, 0x0063182c,
0x00633f89, 0x006366ee, 0x00638e59, 0x0063b5cb,
0x0063dd44, 0x006404c4, 0x00642c4b, 0x006453d8,
0x00647b6d, 0x0064a308, 0x0064caaa, 0x0064f253,
0x00651a03, 0x006541b9, 0x00656977, 0x0065913c,
0x0065b907, 0x0065e0d9, 0x006608b2, 0x00663092,
0x00665879, 0x00668067, 0x0066a85c, 0x0066d058,
0x0066f85b, 0x00672064, 0x00674875, 0x0067708c,
0x006798ab, 0x0067c0d0, 0x0067e8fd, 0x00681130,
0x0068396a, 0x006861ac, 0x006889f4, 0x0068b243,
0x0068da99, 0x006902f7, 0x00692b5b, 0x006953c6,
0x00697c38, 0x0069a4b1, 0x0069cd32, 0x0069f5b9,
0x006a1e47, 0x006a46dd, 0x006a6f79, 0x006a981c,
0x006ac0c7, 0x006ae978, 0x006b1231, 0x006b3af1,
0x006b63b7, 0x006b8c85, 0x006bb55a, 0x006bde36,
0x006c0719, 0x006c3003, 0x006c58f4, 0x006c81ec,
0x006caaec, 0x006cd3f2, 0x006cfd00, 0x006d2614,
0x006d4f30, 0x006d7853, 0x006da17d, 0x006dcaae,
0x006df3e7, 0x006e1d26, 0x006e466d, 0x006e6fbb,
0x006e9910, 0x006ec26c, 0x006eebcf, 0x006f1539,
0x006f3eab, 0x006f6824, 0x006f91a4, 0x006fbb2b,
0x006fe4ba, 0x00700e4f, 0x007037ec, 0x00706190,
0x00708b3b, 0x0070b4ee, 0x0070dea8, 0x00710868,
0x00713231, 0x00715c00, 0x007185d7, 0x0071afb5,
0x0071d99a, 0x00720386, 0x00722d7a, 0x00725775,
0x00728177, 0x0072ab81, 0x0072d592, 0x0072ffaa,
0x007329c9, 0x007353f0, 0x00737e1e, 0x0073a853,
0x0073d290, 0x0073fcd4, 0x0074271f, 0x00745172,
0x00747bcc, 0x0074a62d, 0x0074d096, 0x0074fb06,
0x0075257d, 0x00754ffc, 0x00757a82, 0x0075a50f,
0x0075cfa4, 0x0075fa40, 0x007624e4, 0x00764f8f,
0x00767a41, 0x0076a4fb, 0x0076cfbc, 0x0076fa85,
0x00772555, 0x0077502d, 0x00777b0b, 0x0077a5f2,
0x0077d0df, 0x0077fbd5, 0x007826d1, 0x007851d5,
0x00787ce1, 0x0078a7f4, 0x0078d30e, 0x0078fe30,
0x0079295a, 0x0079548b, 0x00797fc3, 0x0079ab03,
0x0079d64a, 0x007a0199, 0x007a2cf0, 0x007a584d,
0x007a83b3, 0x007aaf20, 0x007ada94, 0x007b0610,
0x007b3194, 0x007b5d1f, 0x007b88b2, 0x007bb44c,
0x007bdfed, 0x007c0b97, 0x007c3748, 0x007c6300,
0x007c8ec0, 0x007cba88, 0x007ce657, 0x007d122e,
0x007d3e0c, 0x007d69f2, 0x007d95e0, 0x007dc1d5,
0x007dedd2, 0x007e19d6, 0x007e45e2, 0x007e71f6,
0x007e9e11, 0x007eca34, 0x007ef65f, 0x007f2291,
0x007f4ecb, 0x007f7b0d, 0x007fa756, 0x007fd3a7
};

#  ifdef FVEC_32
static inline VecF fmath_exp_ps(VecF xxv) {
  __m256 xx = R_CAST(__m256, xxv);
  const __m256i mask7ff = _mm256_set1_epi32(0x7fffffff);
  // 88
  const __m256i max_x = _mm256_set1_epi32(0x42b00000);
  // -88
  // more sensible 0xc2b00000... not used here due to "narrowing conversion"
  // warning
  const __m256i min_x = _mm256_set1_epi64x(-0x3d4fffff3d500000LL);
  // 2^10 / log(2)
  const __m256i const_aa = _mm256_set1_epi32(0x44b8aa3b);
  // log(2) / 2^10
  const __m256i const_bb = _mm256_set1_epi32(0x3a317218);
  const __m256i f1 = _mm256_set1_epi32(0x3f800000);
  const __m256i mask_s = _mm256_set1_epi32(0x000003ff);
  const __m256i i127s = _mm256_set1_epi32(0x0001fc00);
  const __m256i limit = _mm256_castps_si256(_mm256_and_ps(xx, R_CAST(__m256, mask7ff)));
  const int32_t over = _mm256_movemask_epi8(_mm256_cmpgt_epi32(limit, max_x));
  if (over) {
    xx = _mm256_min_ps(xx, R_CAST(__m256, max_x));
    xx = _mm256_max_ps(xx, R_CAST(__m256, min_x));
  }
  const __m256i rr = _mm256_cvtps_epi32(_mm256_mul_ps(xx, R_CAST(__m256, const_aa)));
  __m256 tt = _mm256_fnmadd_ps(_mm256_cvtepi32_ps(rr), R_CAST(__m256, const_bb), xx);
  tt = _mm256_add_ps(R_CAST(__m256, f1), tt);
  const __m256i v8 = _mm256_and_si256(rr, mask_s);
  __m256i u8 = _mm256_add_epi32(rr, i127s);
  u8 = _mm256_srli_epi32(u8, 10);
  u8 = _mm256_slli_epi32(u8, 23);
  __m256i ti = _mm256_i32gather_epi32(R_CAST(const int*, kFloatExpLookupInt), v8, 4);
  __m256 t0 = _mm256_castsi256_ps(ti);
  t0 = _mm256_or_ps(t0, _mm256_castsi256_ps(u8));
  return R_CAST(VecF, _mm256_mul_ps(tt, t0));
}

// For equivalent "normal" C/C++ code, see the non-__LP64__ versions of these
// functions.

// N.B. This requires all mm[] rows to be zero-padded at the end, and there
// can't be nan values at the end of vect[].  (The other way around works too.)
//
// This is currently a bit faster than sgemm and sgemv on my Mac, so it isn't
// appropriate to throw out this code yet.
static inline void MultMatrixDxnVectN(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
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
const float* const kFloatExpLookup = R_CAST(const float*, kFloatExpLookupInt);

static inline VecF fmath_exp_ps(VecF xxv) {
  __m128 xx = xxv;
  const __m128i mask7ff = _mm_set1_epi32(0x7fffffff);

  // 88
  const __m128i max_x = _mm_set1_epi32(0x42b00000);
  // -88
  // more sensible 0xc2b00000... not used here due to "narrowing conversion"
  // warning
  const __m128i min_x = _mm_set1_epi64x(-0x3d4fffff3d500000LL);
  // 2^10 / log(2)
  const __m128i const_aa = _mm_set1_epi32(0x44b8aa3b);
  // log(2) / 2^10
  const __m128i const_bb = _mm_set1_epi32(0x3a317218);

  const __m128i f1 = _mm_set1_epi32(0x3f800000);
  const __m128i mask_s = _mm_set1_epi64x(0x3ff000003ffLLU);
  const __m128i i127s = _mm_set1_epi32(0x0001fc00);
  const __m128i limit = _mm_castps_si128(_mm_and_ps(xx, R_CAST(__m128, mask7ff)));
  const int32_t over = _mm_movemask_epi8(_mm_cmpgt_epi32(limit, max_x));
  if (over) {
    xx = _mm_min_ps(xx, R_CAST(__m128, max_x));
    xx = _mm_max_ps(xx, R_CAST(__m128, min_x));
  }
  const __m128i rr = _mm_cvtps_epi32(_mm_mul_ps(xx, R_CAST(__m128, const_aa)));
  __m128 tt = _mm_sub_ps(xx, _mm_mul_ps(_mm_cvtepi32_ps(rr), R_CAST(__m128, const_bb)));
  tt = _mm_add_ps(tt, R_CAST(__m128, f1));
  const __m128i v4 = _mm_and_si128(rr, mask_s);
  __m128i u4 = _mm_add_epi32(rr, i127s);
  u4 = _mm_srli_epi32(u4, 10);
  u4 = _mm_slli_epi32(u4, 23);
  const uint32_t v0 = _mm_cvtsi128_si32(v4);
  // uint32_t v1 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(2)));
  // uint32_t v2 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(4)));
  // uint32_t v3 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(6)));
  // make this work with LLVM
  const uint32_t v1 = _mm_extract_epi16(R_CAST(__m128i, v4), 2);
  const uint32_t v2 = _mm_extract_epi16(R_CAST(__m128i, v4), 4);
  const uint32_t v3 = _mm_extract_epi16(R_CAST(__m128i, v4), 6);

  __m128 t0 = _mm_set_ss(kFloatExpLookup[v0]);
  __m128 t1 = _mm_set_ss(kFloatExpLookup[v1]);
  const __m128 t2 = _mm_set_ss(kFloatExpLookup[v2]);
  const __m128 t3 = _mm_set_ss(kFloatExpLookup[v3]);
  t1 = _mm_movelh_ps(t1, t3);
  t1 = _mm_castsi128_ps(_mm_slli_epi64(_mm_castps_si128(t1), 32));
  t0 = _mm_movelh_ps(t0, t2);
  t0 = _mm_or_ps(t0, t1);
  t0 = _mm_or_ps(t0, _mm_castsi128_ps(u4));
  tt = _mm_mul_ps(tt, t0);
  return R_CAST(VecF, tt);
}

static inline void MultMatrixDxnVectN(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uint32_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  ColMajorFvectorMatrixMultiplyStrided(vect, mm, col_ct, col_ctav, row_ct, dest);
}

#  endif  // !FVEC_32

static inline void LogisticSse(uint32_t nn, float* vect) {
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

static inline void ComputeVAndPMinusY(const float* yy, uint32_t nn, float* pp, float* vv) {
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF ptmp = *R_CAST(VecF*, &(pp[uii]));
    VecF one_minus_ptmp = one - ptmp;
    *R_CAST(VecF*, &(vv[uii])) = ptmp * one_minus_ptmp;
    VecF ytmp = *R_CAST(const VecF*, &(yy[uii]));
    *R_CAST(VecF*, &(pp[uii])) = ptmp - ytmp;
  }
}

static inline void ComputeV(const float* pp, uint32_t nn, float* vv) {
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF ptmp = *R_CAST(const VecF*, &(pp[uii]));
    VecF one_minus_ptmp = one - ptmp;
    *R_CAST(VecF*, &(vv[uii])) = ptmp * one_minus_ptmp;
  }
}

static inline float TripleProduct(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  VecF sum = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF aa = *R_CAST(const VecF*, &(v1[uii]));
    VecF bb = *R_CAST(const VecF*, &(v2[uii]));
    VecF cc = *R_CAST(const VecF*, &(v3[uii]));
    sum = sum + aa * bb * cc;
  }
  return VecFHsum(sum);
}

static inline void ComputeTwoDiagTripleProduct(const float* aa, const float* bb, const float* vv, uint32_t nn, float* __restrict raa_ptr, float* __restrict rab_ptr, float* __restrict rbb_ptr) {
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

static inline void ComputeThreeTripleProduct(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
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

static inline void ComputeTwoPlusOneTripleProduct(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
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
static inline void LogisticSse(uint32_t nn, float* vect) {
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

static inline void ComputeVAndPMinusY(const float* yy, uint32_t nn, float* pp, float* vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vv[uii] = pp[uii] * (S_CAST(float, 1.0) - pp[uii]);
    pp[uii] -= yy[uii];
  }
}

static inline void ComputeV(const float* pp, uint32_t nn, float* vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vv[uii] = pp[uii] * (S_CAST(float, 1.0) - pp[uii]);
  }
}

static inline void MultMatrixDxnVectN(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uint32_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  ColMajorFvectorMatrixMultiplyStrided(vect, mm, col_ct, col_ctav, row_ct, dest);
}

static inline float TripleProduct(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  float fxx = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    fxx += (*v1++) * (*v2++) * (*v3++);
  }
  return fxx;
}

static inline void ComputeTwoDiagTripleProduct(const float* aa, const float* bb, const float* vv, uint32_t nn, float* raa_ptr, float* rab_ptr, float* rbb_ptr) {
  float raa = 0.0;
  float rab = 0.0;
  float rbb = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = (*aa++);
    const float fyy = (*bb++);
    float fzz = (*vv++);
    raa += fxx * fxx * fzz;
    fzz *= fyy;
    rab += fxx * fzz;
    rbb += fyy * fzz;
  }
  *raa_ptr = raa;
  *rab_ptr = rab;
  *rbb_ptr = rbb;
}

static inline void ComputeThreeTripleProduct(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* r1_ptr, float* r2_ptr, float* r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = (*bb++) * (*vv++);
    r1 += (*a1++) * fxx;
    r2 += (*a2++) * fxx;
    r3 += (*a3++) * fxx;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}

static inline void ComputeTwoPlusOneTripleProduct(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* r1_ptr, float* r2_ptr, float* r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = (*bb++);
    const float fyy = fxx * (*vv++);
    r1 += fxx * fyy;
    r2 += (*a1++) * fyy;
    r3 += (*a2++) * fyy;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}
#endif
double ComputeLoglik(const float* yy, const float* pp, uint32_t sample_ct) {
  // possible todo: look for a high-precision way to accelerate this.
  double loglik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const double new_pi = S_CAST(double, pp[sample_idx]);
    loglik += (yy[sample_idx] != S_CAST(float, 0.0))? log(new_pi) : log1p(-new_pi);
  }
  return loglik;
}

// M V M^T
// This is the biggest logistic/Firth regression bottleneck.
// Tried to replace this with sqrt(v) followed by ssyrk, but that was slower.
// Also tried to take advantage of the first row of M being constant-1, and
// managed to fail.
void ComputeHessian(const float* mm, const float* vv, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uintptr_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  const uintptr_t row_ctav = RoundUpPow2(row_ct, kFloatPerFVec);
  const uintptr_t row_ctavp1 = row_ctav + 1;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    for (uint32_t row_idx = 0; row_idx < row_ctm3; row_idx += 3) {
      const float* mm_cur = &(mm[row_idx * col_ctav]);
      ComputeTwoDiagTripleProduct(mm_cur, &(mm_cur[col_ctav]), vv, col_ct, &(dest[row_idx * row_ctavp1]), &(dest[(row_idx + 1) * row_ctavp1 - 1]), &(dest[(row_idx + 1) * row_ctavp1]));
      ComputeTwoPlusOneTripleProduct(&(mm_cur[2 * col_ctav]), &(mm_cur[col_ctav]), mm_cur, vv, col_ct, &(dest[(row_idx + 2) * row_ctavp1]), &(dest[(row_idx + 2) * row_ctavp1 - 1]), &(dest[(row_idx + 2) * row_ctavp1 - 2]));
      for (uint32_t row_idx2 = row_idx + 3; row_idx2 != row_ct; ++row_idx2) {
        ComputeThreeTripleProduct(&(mm[row_idx2 * col_ctav]), mm_cur, &(mm_cur[col_ctav]), &(mm_cur[2 * col_ctav]), vv, col_ct, &(dest[row_idx2 * row_ctav + row_idx]), &(dest[row_idx2 * row_ctav + row_idx + 1]), &(dest[row_idx2 * row_ctav + row_idx + 2]));
      }
    }
  }
  switch (row_ct % 3) {
  case 0:
    ComputeTwoPlusOneTripleProduct(&(mm[(row_ct - 3) * col_ctav]), &(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 3) * row_ctavp1]), &(dest[(row_ct - 2) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1 - 2]));
    // fall through
  case 2:
    ComputeTwoDiagTripleProduct(&(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 2) * row_ctavp1]), &(dest[(row_ct - 1) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1]));
    break;
  case 1:
    dest[(row_ct - 1) * row_ctavp1] = TripleProduct(&(mm[(row_ct - 1) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct);
  }
}

void CholeskyDecomposition(const float* aa, uint32_t predictor_ct, float* ll) {
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  const uintptr_t predictor_ctavp1 = predictor_ctav + 1;
  for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
    float fxx = aa[row_idx * predictor_ctavp1];
    float* ll_row_iter = &(ll[row_idx * predictor_ctav]);
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      const float fyy = (*ll_row_iter++);
      fxx -= fyy * fyy;
    }
    float fyy;
    if (fxx >= S_CAST(float, 0.0)) {
      fyy = sqrtf(fxx);
    } else {
      fyy = S_CAST(float, 1e-6);
    }
    ll[row_idx * predictor_ctavp1] = fyy;
    fyy = S_CAST(float, 1.0) / fyy;  // now 1.0 / L[j][j]
    for (uint32_t row_idx2 = row_idx + 1; row_idx2 != predictor_ct; ++row_idx2) {
      float fxx2 = aa[row_idx2 * predictor_ctav + row_idx];
      float* ll_row_iter2 = &(ll[row_idx * predictor_ctav]);
      float* ll_row_iter3 = &(ll[row_idx2 * predictor_ctav]);
      for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
        fxx2 -= (*ll_row_iter2++) * (*ll_row_iter3++);
      }
      ll[row_idx2 * predictor_ctav + row_idx] = fxx2 * fyy;
    }
  }
}

void SolveLinearSystem(const float* ll, const float* yy, uint32_t predictor_ct, float* xx) {
  // Finds x such that y = L(L^T)x, via forward and backward substitution
  //
  // might want to use this in NOLAPACK case only, since we can now produce
  // 32-bit Linux builds with statically linked LAPACK
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
    float fxx = yy[row_idx];
    const float* ll_row_iter = &(ll[row_idx * predictor_ctav]);
    float* xx_iter = xx;
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      fxx -= (*ll_row_iter++) * (*xx_iter++);
    }
    *xx_iter = fxx / (*ll_row_iter);
  }
  for (uint32_t col_idx = predictor_ct; col_idx; ) {
    float fxx = xx[--col_idx];
    float* xx_iter = &(xx[predictor_ct - 1]);
    for (uint32_t row_idx = predictor_ct - 1; row_idx > col_idx; --row_idx) {
      fxx -= ll[row_idx * predictor_ctav + col_idx] * (*xx_iter--);
    }
    *xx_iter = fxx / ll[col_idx * (predictor_ctav + 1)];
  }
}

BoolErr LogisticRegression(const float* yy, const float* xx, const float* sample_offsets, uint32_t sample_ct, uint32_t predictor_ct, float* coef, uint32_t* is_unfinished_ptr, float* ll, float* pp, float* vv, float* hh, float* grad, float* dcoef) {
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
  // coef  = starting point, overwritten with logistic regression betas.  Must
  //         be vector-16-byte waligned.
  //
  // Outputs:
  // ll    = cholesky decomposition matrix, predictor_ct^2, rows vector-aligned
  // pp    = final likelihoods minus Y[] (not currently used by callers)
  //
  // Returns 0 on success, 1 on convergence failure.
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  float min_delta_coef = 1e9;

  ZeroFArr(predictor_ct * predictor_ctav, ll);
  for (uint32_t iteration = 0; ; ++iteration) {
    // P[i] = \sum_j X[i][j] * coef[j];
    ColMajorFmatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);
    if (sample_offsets) {
      AddFVec(sample_offsets, sample_ctav, pp);
    }
    // Suppose categorical covariates are represented as
    // categorical-covariate-major uint16_t* kk, indicating for each sample
    // which raw covariate index is 1 (with one "fallow" covariate index at the
    // end).
    // Then the above expression becomes
    // P[i] = \sum_j^{regular} X[i][j] * coef[j] +
    //        \sum_j^{cats} coef[K[i][j]]

    // P[i] = 1 / (1 + exp(-P[i]));
    LogisticSse(sample_ct, pp);

    // V[i] = P[i] * (1 - P[i]);
    // P[i] -= Y[i];
    ComputeVAndPMinusY(yy, sample_ct, pp, vv);

    // Possible categorical optimizations:
    // 1. skip terms between different categories of the same covariate
    // 2. all same-category terms within the same covariate can be handled with
    //    a single loop over the samples
    // 3. terms involving a regular covariate and a categorical covariate can
    //    be handled with one loop over the samples; covers all categories
    // 4. similarly, one loop over the samples is enough to update all category
    //    pairs between two categorical covariates
    ComputeHessian(xx, vv, sample_ct, predictor_ct, hh);

    // grad = X^T P
    // Separate categorical loop also possible here
    MultMatrixDxnVectN(xx, pp, sample_ct, predictor_ct, grad);

    CholeskyDecomposition(hh, predictor_ct, ll);

    // ZeroFArr(predictor_ct, dcoef);
    SolveLinearSystem(ll, grad, predictor_ct, dcoef);

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
void CopyAndMeanCenterF(const float* src, uintptr_t ct, float* dst) {
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
void CopyAndMeanCenterF(const float* src, uintptr_t ct, float* dst) {
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

BoolErr LogisticRegressionResidualized(const float* yy, const float* xx, const uintptr_t* sample_nm, const CcResidualizeCtx* cc_residualize, uint32_t nm_sample_ct, uint32_t orig_predictor_ct, float* coef, uint32_t* is_unfinished_ptr, float* ll, float* pp, float* vv, float* hh, float* grad, float* dcoef, float* mean_centered_pmaj_buf, float* sample_offsets_buf) {
  if (!cc_residualize->logistic_nm_sample_offsets) {
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
  const float* sample_offsets = cc_residualize->logistic_nm_sample_offsets;
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
  if (LogisticRegression(yy, mean_centered_pmaj_buf, sample_offsets, nm_sample_ct, regressed_predictor_ct, &(coef[1]), is_unfinished_ptr, ll, pp, vv, hh, grad, dcoef)) {
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

#ifdef __LP64__
// tmpNxK, interpreted as column-major, is sample_ct x predictor_ct
// X, interpreted as column-major, is also sample_ct x predictor_ct
// Hdiag[i] = V[i] (\sum_j tmpNxK[i][j] X[i][j])
void FirthComputeWeights(const float* yy, const float* xx, const float* pp, const float* vv, const float* tmpnxk, uint32_t predictor_ct, __maybe_unused uint32_t sample_ct, uint32_t sample_ctav, float* ww) {
  const VecF half = VCONST_F(0.5);
  for (uint32_t sample_offset = 0; sample_offset < sample_ctav; sample_offset += kFloatPerFVec) {
    VecF dotprods = vecf_setzero();
    const float* xx_row = &(xx[sample_offset]);
    const float* tmpnxk_row = &(tmpnxk[sample_offset]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const VecF cur_xx = *R_CAST(const VecF*, &(xx_row[pred_uidx * sample_ctav]));
      const VecF cur_tmpnxk = *R_CAST(const VecF*, &(tmpnxk_row[pred_uidx * sample_ctav]));
      dotprods = dotprods + cur_xx * cur_tmpnxk;
    }
    // Can handle categorical covariates in a separate loop here, and load the
    // dotprods increment into a union, etc.
    const VecF cur_weights = *R_CAST(const VecF*, &(vv[sample_offset]));
    const VecF cur_pis = *R_CAST(const VecF*, &(pp[sample_offset]));
    const VecF cur_yy = *R_CAST(const VecF*, &(yy[sample_offset]));
    const VecF half_minus_cur_pis = half - cur_pis;
    const VecF yy_minus_cur_pis = cur_yy - cur_pis;
    const VecF second_term = half_minus_cur_pis * (cur_weights * dotprods);
    const VecF cur_wws = yy_minus_cur_pis + second_term;
    *R_CAST(VecF*, &(ww[sample_offset])) = cur_wws;
  }
}
#else
void FirthComputeWeights(const float* yy, const float* xx, const float* pp, const float* vv, const float* tmpnxk, uint32_t predictor_ct, uint32_t sample_ct, uint32_t sample_ctav, float* ww) {
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    float dotprod = 0.0;
    const float* xx_row = &(xx[sample_idx]);
    const float* tmpnxk_row = &(tmpnxk[sample_idx]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      dotprod += xx_row[pred_uidx * sample_ctav] * tmpnxk_row[pred_uidx * sample_ctav];
    }
    const float cur_weight = vv[sample_idx];
    const float cur_pi = pp[sample_idx];
    ww[sample_idx] = (yy[sample_idx] - cur_pi) + (S_CAST(float, 0.5) - cur_pi) * cur_weight * dotprod;
  }
}
#endif

BoolErr FirthRegression(const float* yy, const float* xx, const float* sample_offsets, uint32_t sample_ct, uint32_t predictor_ct, float* coef, uint32_t* is_unfinished_ptr, float* hh, double* half_inverted_buf, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, float* pp, float* vv, float* grad, float* dcoef, float* ww, float* tmpnxk_buf) {
  // This is a port of Georg Heinze's logistf R function, adapted to use many
  // of plink 1.9's optimizations; see
  //   http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/fllogistf/
  //
  // Preallocated buffers (initial contents irrelevant):
  // half_inverted_buf, inv_1d_buf, dbl_2d_buf = for matrix inversion
  // pp    = likelihoods minus Y[] (not currently used by callers)
  // vv    = sample variance buffer
  // grad  = gradient buffer (length predictor_ct)
  // dcoef = current coefficient change buffer (length predictor_ct)
  // ww    = Firth-adjusted scores, sample_ct
  //
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         vector-aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype
  // sample_offsets = in residualized mode, product of pre-fitted covariate
  //                  beta vector with covariate matrix.  otherwise, nullptr.
  //
  // Input/output:
  // coef  = starting point, overwritten with logistic regression betas.  Must
  //         be vector-aligned.
  //
  // Outputs:
  // hh    = variance-covariance matrix buffer, predictor_ct^2, rows
  //         vector-aligned.  (spends some time as pre-inversion Hessian matrix
  //         too)
  //
  // Returns 0 on success, 1 on convergence failure.
  // is_unfinished assumed to be initialized to 0, and is set to 1 if we hit
  // the iteration limit without satisfying other convergence criteria.
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  uint32_t is_last_iter = 0;

  // pull these out of the start of the loop, since they happen again in the
  // log-likelihood update
  // P[i] = \sum_j coef[j] * X[i][j];
  // categorical optimization possible here
  ColMajorFmatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);
  if (sample_offsets) {
    AddFVec(sample_offsets, sample_ctav, pp);
  }
  // P[i] = 1 / (1 + exp(-P[i]));
  LogisticSse(sample_ct, pp);
  // V[i] = P[i] * (1 - P[i]);
  ComputeV(pp, sample_ct, vv);
  // P[i] -= Y[i] NOT done here

  // hessian = X diag(V) X'
  // note that only lower triangle is filled here
  ComputeHessian(xx, vv, sample_ct, predictor_ct, hh);

  // we shouldn't need to compute the log directly, since underflow <->
  // regression failure, right?  check this.
  if (InvertSymmdefFmatrixFirstHalf(predictor_ct, predictor_ctav, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf)) {
    return 1;
  }
  double dethh = HalfSymmInvertedDet(half_inverted_buf, inv_1d_buf, predictor_ct);
  double loglik = ComputeLoglik(yy, pp, sample_ct);
  // printf("loglik: %g\n", loglik);
  loglik += 0.5 * log(dethh);

  // bugfix (4 Nov 2017): grad[] trailing elements must be zeroed out
  ZeroFArr(predictor_ctav - predictor_ct, &(grad[predictor_ct]));

  // start with 80% of most logistf convergence defaults (some reduction is
  // appropriate to be consistent with single-precision arithmetic).  (Update,
  // 9 Apr 2020: max_iter now matches logistf default since we've been shown a
  // concrete example where it matters in a way that isn't masked by the
  // single- vs. double-precision difference.)
  //
  // see also the hs_bail condition: if we ever try all five halfsteps, when
  // dcoef_max and grad_max aren't that far from the normal convergence
  // conditions, it's probably pointless to continue with single-precision
  // arithmetic.  (possible todo: use a fully-double-precision routine to
  // finish the job when that happens.)
  const uint32_t max_iter_m1 = 24;
  const float gconv = S_CAST(float, 0.0001);
  const float xconv = S_CAST(float, 0.0001);
  const double lconv = 0.0001;
  uint32_t hs_bail = 0;
  for (uint32_t iter_idx = 0; ; ++iter_idx) {
    InvertSymmdefFmatrixSecondHalf(predictor_ct, predictor_ctav, half_inverted_buf, hh, inv_1d_buf, dbl_2d_buf);
    if (is_last_iter) {
      return 0;
    }
    // bugfix (13 Oct 2017): trailing elements of hh[] rows can't be arbitrary
    // for later MultMatrixDxnVectN() call
    ReflectFmatrix0(predictor_ct, predictor_ctav, hh);

    // categorical optimization possible here
    ColMajorFmatrixMultiplyStrided(xx, hh, sample_ct, sample_ctav, predictor_ct, predictor_ctav, predictor_ct, sample_ctav, tmpnxk_buf);

    FirthComputeWeights(yy, xx, pp, vv, tmpnxk_buf, predictor_ct, sample_ct, sample_ctav, ww);

    // trailing elements of ww can't be nan for MultMatrixDxnVectN()
    ZeroFArr(sample_ctav - sample_ct, &(ww[sample_ct]));

    // gradient (Ustar in logistf) = X' W
    // categorical optimization possible here
    MultMatrixDxnVectN(xx, ww, sample_ct, predictor_ct, grad);
    float grad_max = 0.0;
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const float abs_grad_cur = fabsf(grad[pred_uidx]);
      if (abs_grad_cur > grad_max) {
        grad_max = abs_grad_cur;
      }
    }

    // dcoef := hh * grad (note that hh is inverted already)
    MultMatrixDxnVectN(hh, grad, predictor_ct, predictor_ct, dcoef);

    float dcoef_max = 0.0;
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const float abs_dcoef_cur = fabsf(dcoef[pred_uidx]);
      if (abs_dcoef_cur > dcoef_max) {
        dcoef_max = abs_dcoef_cur;
      }
    }
    const float maxstep = 5.0;
    if (dcoef_max > maxstep) {
      const float scaling_factor = maxstep / dcoef_max;
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        dcoef[pred_uidx] *= scaling_factor;
      }
      dcoef_max = maxstep;
    }
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      coef[pred_uidx] += dcoef[pred_uidx];
    }
    const uint32_t delta_and_grad_converged = (dcoef_max <= xconv) && (grad_max < gconv);
    const double loglik_old = loglik;
    double loglik_thresh = loglik_old;
    if (delta_and_grad_converged) {
      // on the last iteration, we would frequently try all 5 halfsteps when
      // the log-likelihood change was effectively random due to floating point
      // error.  detect this and exit the loop earlier.
      loglik_thresh -= 0.999999 * lconv;
    }

    uint32_t maxhs = 5;
    uint32_t halfstep_idx = 1;
    while (1) {
      // categorical optimization possible here
      ColMajorFmatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);
      if (sample_offsets) {
        AddFVec(sample_offsets, sample_ctav, pp);
      }

      LogisticSse(sample_ct, pp);
      loglik = ComputeLoglik(yy, pp, sample_ct);
      ComputeV(pp, sample_ct, vv);
      ComputeHessian(xx, vv, sample_ct, predictor_ct, hh);
      if (InvertSymmdefFmatrixFirstHalf(predictor_ct, predictor_ctav, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf)) {
        return 1;
      }
      dethh = HalfSymmInvertedDet(half_inverted_buf, inv_1d_buf, predictor_ct);
      loglik += 0.5 * log(dethh);
      if (halfstep_idx > maxhs) {
        break;
      }
      if (loglik >= loglik_thresh) {
        if (loglik >= loglik_old) {
          break;
        }
        maxhs = halfstep_idx;
      } else if (halfstep_idx == maxhs) {
        if ((dcoef_max < S_CAST(float, 0.001)) && (grad_max < S_CAST(float, 0.05)) && (loglik >= loglik_old - lconv)) {
          // we've converged as much as we can with single-precision
          // arithmetic, and now we're flailing around.  don't even take the
          // 2^{-maxhs} step, undo it all and bail.
          // (0.001 and 0.05 constants can obviously be tuned; they were chosen
          // based on a test 500k sample/5 covariate regression.)
          --halfstep_idx;
          --maxhs;
          hs_bail = 1;
        }
      }
      const float multiplier = exp2f(-u31tof(halfstep_idx));
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        coef[pred_uidx] -= dcoef[pred_uidx] * multiplier;
      }
      ++halfstep_idx;
    }
    // printf("%.9g %.9g %g %g\n", loglik, loglik_old, dcoef_max, grad_max);
    const double loglik_change = loglik - loglik_old;
    if ((fabs(loglik_change) <= lconv) && (delta_and_grad_converged || hs_bail)) {
      is_last_iter = 1;
    } else if (iter_idx == max_iter_m1) {
      is_last_iter = 1;
      *is_unfinished_ptr = 1;
    }
  }
}

BoolErr FirthRegressionResidualized(const float* yy, const float* xx, const uintptr_t* sample_nm, const CcResidualizeCtx* cc_residualize, uint32_t nm_sample_ct, uint32_t orig_predictor_ct, float* coef, uint32_t* is_unfinished_ptr, float* hh, double* half_inverted_buf, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, float* pp, float* vv, float* grad, float* dcoef, float* ww, float* tmpnxk_buf, float* mean_centered_pmaj_buf, float* sample_offsets_buf) {
  // todo: deduplicate with LogisticRegressionResidualized()
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
  const float* sample_offsets = cc_residualize->firth_nm_sample_offsets;
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
  if (FirthRegression(yy, mean_centered_pmaj_buf, sample_offsets, nm_sample_ct, regressed_predictor_ct, &(coef[1]), is_unfinished_ptr, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf, pp, vv, grad, dcoef, ww, tmpnxk_buf)) {
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

static const float kSmallFloats[4] = {0.0, 1.0, 2.0, 3.0};

BoolErr GlmAllocFillAndTestPhenoCovarsCc(const uintptr_t* sample_include, const uintptr_t* pheno_cc, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uint32_t domdev_present_p1, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, GlmFlags glm_flags, uintptr_t** pheno_cc_collapsed_ptr, uintptr_t** gcount_case_interleaved_vec_ptr, float** pheno_f_ptr, RegressionNmPrecomp** nm_precomp_ptr, float** covars_cmaj_f_ptr, CcResidualizeCtx** cc_residualize_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
  const uint32_t is_cc_residualize = !!(glm_flags & (kfGlmFirthResidualize | kfGlmCcResidualize));
  if (unlikely(bigstack_alloc_w(sample_ctv * kWordsPerVec, pheno_cc_collapsed_ptr) ||
               bigstack_alloc_f(sample_ctav, pheno_f_ptr) ||
               bigstack_alloc_f(new_nonlocal_covar_ct * sample_ctav, covars_cmaj_f_ptr) ||
               bigstack_alloc_kcp(new_covar_ct, cur_covar_names_ptr))) {
      return 1;
  }
  double* corr_buf = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  *nm_precomp_ptr = nullptr;
  if (xtx_state) {
    assert(!local_covar_ct);
    if (is_cc_residualize) {
      // Interactions and local covariates prohibited with firth-residualize,
      // so xtx_state guaranteed to be nonzero.
      if (unlikely(BIGSTACK_ALLOC_X(CcResidualizeCtx, 1, cc_residualize_ptr))) {
        return 1;
      }
      (*cc_residualize_ptr)->logistic_nm_sample_offsets = nullptr;
      if ((glm_flags & (kfGlmFirth | kfGlmCcResidualize)) == kfGlmCcResidualize) {
        if (unlikely(bigstack_alloc_f(sample_ctav, &((*cc_residualize_ptr)->logistic_nm_sample_offsets)))) {
          return 1;
        }
      }
      (*cc_residualize_ptr)->firth_nm_sample_offsets = nullptr;
      if (!(glm_flags & kfGlmNoFirth)) {
        if (unlikely(bigstack_alloc_f(sample_ctav, &((*cc_residualize_ptr)->firth_nm_sample_offsets)))) {
          return 1;
        }
      }
      (*cc_residualize_ptr)->prefitted_pred_ct = 1 + new_covar_ct;
      (*cc_residualize_ptr)->domdev_present_p1 = domdev_present_p1;
      (*cc_residualize_ptr)->sample_ct = sample_ct;
    }
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
  uintptr_t* pheno_cc_collapsed = *pheno_cc_collapsed_ptr;
  CopyBitarrSubset(pheno_cc, sample_include, sample_ct, pheno_cc_collapsed);
  float* pheno_f_iter = *pheno_f_ptr;
  for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    // can use the bitvector equivalent of GenoarrLookup...(), but this isn't
    // in a critical loop so I'll postpone writing those functions for now
    *pheno_f_iter++ = kSmallFloats[IsSet(pheno_cc_collapsed, sample_idx)];
  }
  const uint32_t sample_remv = sample_ctav - sample_ct;
  ZeroFArr(sample_remv, pheno_f_iter);
  PglErr reterr = GlmFillAndTestCovars(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, max_corr, vif_thresh, covar_dotprod, corr_buf, inverse_corr_buf, covars_cmaj_d, *cur_covar_names_ptr, glm_err_ptr);
  if (unlikely(reterr)) {
    return (reterr == kPglRetNomem);
  }
  double* covar_read_iter = covars_cmaj_d;
  float* covar_write_iter = *covars_cmaj_f_ptr;
  for (uintptr_t covar_idx = 0; covar_idx != new_nonlocal_covar_ct; ++covar_idx) {
    for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      *covar_write_iter++ = S_CAST(float, *covar_read_iter++);
    }
    ZeroFArr(sample_remv, covar_write_iter);
    covar_write_iter = &(covar_write_iter[sample_remv]);
  }
  if (xtx_state) {
    // error-out should be impossible
    InitNmPrecomp(covars_cmaj_d, covar_dotprod, corr_buf, inverse_corr_buf, sample_ct, 0, new_covar_ct, xtx_state, *nm_precomp_ptr);
    if (is_cc_residualize) {
      BigstackReset(bigstack_mark);
      const uintptr_t pred_ct = new_covar_ct + 1;
      const uintptr_t pred_ctav = RoundUpPow2(pred_ct, kFloatPerFVec);
      float* xx;
      float* coefs;
      float* hh;
      double* half_inverted_buf;
      double* dbl_2d_buf;
      float* ll;
      float* pp;
      float* vv;
      float* grad;
      float* dcoef;
      float* ww;
      float* tmpnxk_buf;
      MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, bigstack_alloc(pred_ct * kMatrixInvertBuf1CheckedAlloc));
      // see GetLogisticWorkspaceSize() and corresponding code at top of
      // GlmLogisticThread()
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
                   bigstack_alloc_f(sample_ctav, &ww) ||
                   bigstack_alloc_f(sample_ctav * pred_ct, &tmpnxk_buf))) {
        return 1;
      }
      FillFVec(sample_ct, 1.0, xx);
      memcpy(&(xx[sample_ctav]), *covars_cmaj_f_ptr, sample_ctav * (pred_ct - 1) * sizeof(float));
      uint32_t is_unfinished = 0;
      float* logistic_nm_sample_offsets = (*cc_residualize_ptr)->logistic_nm_sample_offsets;
      if (logistic_nm_sample_offsets) {
        if (unlikely(LogisticRegression(*pheno_f_ptr, xx, nullptr, sample_ct, pred_ct, coefs, &is_unfinished, ll, pp, vv, hh, grad, dcoef) || is_unfinished)) {
          if (glm_flags & kfGlmNoFirth) {
            *glm_err_ptr = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
            return 0;
          }
          (*cc_residualize_ptr)->logistic_nm_sample_offsets = nullptr;
        }
        ColMajorFmatrixVectorMultiplyStrided(xx, coefs, sample_ct, sample_ctav, pred_ct, logistic_nm_sample_offsets);
        ZeroFArr(sample_ctav - sample_ct, &(logistic_nm_sample_offsets[sample_ct]));
      }
      float* firth_nm_sample_offsets = (*cc_residualize_ptr)->firth_nm_sample_offsets;
      if (firth_nm_sample_offsets) {
        if (unlikely(FirthRegression(*pheno_f_ptr, xx, nullptr, sample_ct, pred_ct, coefs, &is_unfinished, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf, pp, vv, grad, dcoef, ww, tmpnxk_buf) || is_unfinished)) {
          *glm_err_ptr = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
          return 0;
        }
        ColMajorFmatrixVectorMultiplyStrided(xx, coefs, sample_ct, sample_ctav, pred_ct, firth_nm_sample_offsets);
        ZeroFArr(sample_ctav - sample_ct, &(firth_nm_sample_offsets[sample_ct]));
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

uintptr_t GetLogisticWorkspaceSize(uint32_t sample_ct, uint32_t biallelic_predictor_ct, uint32_t domdev_present_p1, uint32_t max_extra_allele_ct, uint32_t constraint_ct, uint32_t xmain_ct, uint32_t gcount_cc, uint32_t is_sometimes_firth, uint32_t is_cc_residualize) {
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

  // dbl_2d_buf = max_predictor_ct * max_predictor_ctav floats, or VIF/Firth
  // dbl in practice, the latter value is never smaller due to the max(x, 7)
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
    // ww = sample_ctav floats
    workspace_size += RoundUpPow2(sample_ctav * sizeof(float), kCacheline);

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


THREAD_FUNC_DECL GlmLogisticThread(void* raw_arg) {
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
      const uint32_t chr_fo_idx = CountSortedSmallerU32(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
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
      float* score_buf = nullptr;
      float* tmpnxk_buf = nullptr;
      if (is_sometimes_firth) {
        score_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
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
      assert(S_CAST(uintptr_t, workspace_iter - workspace_buf) == GetLogisticWorkspaceSize(cur_sample_ct, cur_biallelic_predictor_ct, domdev_present_p1, max_extra_allele_ct, cur_constraint_ct, main_mutated + main_omitted, cur_gcount_case_interleaved_vec != nullptr, is_sometimes_firth, cur_cc_residualize != nullptr));
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
          goto GlmLogisticThread_err;
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
                goto GlmLogisticThread_skip_regression;
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
                goto GlmLogisticThread_skip_regression;
              }
            } else {
              glm_err = CheckMaxCorrAndVifF(&(nm_predictors_pmaj_buf[nm_sample_ctav]), cur_predictor_ct - 1, nm_sample_ct, nm_sample_ctav, max_corr, vif_thresh, predictor_dotprod_buf, dbl_2d_buf, inverse_corr_buf, inv_1d_buf);
              if (glm_err) {
                goto GlmLogisticThread_skip_regression;
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
                    goto GlmLogisticThread_firth_fallback;
                  }
                  glm_err = SetGlmErr1(kGlmErrcodeSeparation, allele_idx);
                  goto GlmLogisticThread_skip_regression;
                }
              }
              if (!cur_cc_residualize) {
                if (LogisticRegression(nm_pheno_buf, nm_predictors_pmaj_buf, nullptr, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, cholesky_decomp_return, pp_buf, sample_variance_buf, hh_return, gradient_buf, dcoef_buf)) {
                  if (is_sometimes_firth) {
                    ZeroFArr(cur_predictor_ctav, coef_return);
                    goto GlmLogisticThread_firth_fallback;
                  }
                  glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
                  goto GlmLogisticThread_skip_regression;
                }
              } else {
                if (LogisticRegressionResidualized(nm_pheno_buf, nm_predictors_pmaj_buf, sample_nm, cur_cc_residualize, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, cholesky_decomp_return, pp_buf, sample_variance_buf, hh_return, gradient_buf, dcoef_buf, mean_centered_pmaj_buf, sample_offsets_buf)) {
                  if (is_sometimes_firth) {
                    ZeroFArr(cur_predictor_ctav, coef_return);
                    goto GlmLogisticThread_firth_fallback;
                  }
                  glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
                  goto GlmLogisticThread_skip_regression;
                }
                is_residualized = 1;
                cur_regressed_predictor_stop = domdev_present + allele_ct;
                cur_regressed_predictor_ctav = RoundUpPow2(cur_regressed_predictor_stop, kFloatPerFVec);
                cur_regressed_predictor_ctavp1 = cur_regressed_predictor_ctav + 1;
                cur_biallelic_regressed_predictor_stop = domdev_present + 2;
              }
              // unlike FirthRegression(), hh_return isn't inverted yet, do
              // that here
              for (uint32_t pred_uidx = is_residualized; pred_uidx != cur_regressed_predictor_stop; ++pred_uidx) {
                float* hh_inv_row = &(hh_return[pred_uidx * cur_regressed_predictor_ctav]);
                // ZeroFArr(cur_regressed_predictor_stop, gradient_buf);
                // gradient_buf[pred_uidx] = 1.0;
                // (y is gradient_buf, x is dcoef_buf)
                // SolveLinearSystem(cholesky_decomp_return, &(gradient_buf[is_residualized]), cur_regressed_predictor_stop - is_residualized, &(hh_inv_row[is_residualized]));
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
              GlmLogisticThread_firth_fallback:
                block_aux_iter[extra_regression_idx].firth_fallback = 1;
                if (allele_ct_m2 && beta_se_multiallelic_fused) {
                  for (uint32_t uii = 1; uii != allele_ct - 1; ++uii) {
                    block_aux_iter[uii].firth_fallback = 1;
                  }
                }
              }
              if (!cur_cc_residualize) {
                if (FirthRegression(nm_pheno_buf, nm_predictors_pmaj_buf, nullptr, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, hh_return, inverse_corr_buf, inv_1d_buf, dbl_2d_buf, pp_buf, sample_variance_buf, gradient_buf, dcoef_buf, score_buf, tmpnxk_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
                  goto GlmLogisticThread_skip_regression;
                }
              } else {
                if (FirthRegressionResidualized(nm_pheno_buf, nm_predictors_pmaj_buf, sample_nm, cur_cc_residualize, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, hh_return, inverse_corr_buf, inv_1d_buf, dbl_2d_buf, pp_buf, sample_variance_buf, gradient_buf, dcoef_buf, score_buf, tmpnxk_buf, mean_centered_pmaj_buf, sample_offsets_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
                  goto GlmLogisticThread_skip_regression;
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
                goto GlmLogisticThread_skip_regression;
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
                  goto GlmLogisticThread_skip_regression;
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
            GlmLogisticThread_skip_regression:
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
    GlmLogisticThread_err:
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
PglErr GlmLogistic(const char* cur_pheno_name, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const char* outname, uint32_t raw_variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLogisticCtx* ctx, TextStream* local_covar_txsp, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, double* orig_permstat, uintptr_t* valid_allele_ct_ptr) {
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

    const uint32_t main_omitted = parameter_subset && (!IsSet(parameter_subset, 1));
    const uint32_t xmain_ct = main_mutated + main_omitted;
    const uint32_t gcount_cc_col = glm_cols & kfGlmColGcountcc;
    // workflow is similar to --make-bed
    uintptr_t workspace_alloc = GetLogisticWorkspaceSize(sample_ct, biallelic_predictor_ct, domdev_present_p1, max_extra_allele_ct, constraint_ct, xmain_ct, gcount_cc_col, is_sometimes_firth, is_cc_residualize);
    if (sample_ct_x) {
      const uintptr_t workspace_alloc_x = GetLogisticWorkspaceSize(sample_ct_x, biallelic_predictor_ct_x, domdev_present_p1, max_extra_allele_ct, constraint_ct_x, xmain_ct, gcount_cc_col, is_sometimes_firth, is_cc_residualize);
      if (workspace_alloc_x > workspace_alloc) {
        workspace_alloc = workspace_alloc_x;
      }
    }
    if (sample_ct_y) {
      const uintptr_t workspace_alloc_y = GetLogisticWorkspaceSize(sample_ct_y, biallelic_predictor_ct_y, domdev_present_p1, max_extra_allele_ct, constraint_ct_y, xmain_ct, gcount_cc_col, is_sometimes_firth, is_cc_residualize);
      if (workspace_alloc_y > workspace_alloc) {
        workspace_alloc = workspace_alloc_y;
      }
    }
    // +1 is for top-level common->workspace_bufs
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;

    uintptr_t per_variant_xalloc_byte_ct = max_sample_ct * local_covar_ct * sizeof(float);
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
        if (unlikely(bigstack_alloc_f(read_block_size * max_sample_ct * local_covar_ct, &(ctx->local_covars_vcmaj_f[uii])))) {
          goto GlmLogistic_ret_NOMEM;
        }
      } else {
        ctx->local_covars_vcmaj_f[uii] = nullptr;
      }
    }

    if (max_sample_ct > 2000000) {
      // may want a large-matrix double-precision fallback
      logerrputs("Warning: --glm logistic regression is unreliable on more than ~2 million\nsamples, since it uses single-precision arithmetic.\n");
    }
    common->workspace_bufs = S_CAST(unsigned char**, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t)));
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      common->workspace_bufs[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(workspace_alloc));
    }
    common->err_info = (~0LLU) << 32;
    SetThreadFuncAndData(GlmLogisticThread, ctx, &tg);

    const uint32_t ref_col = glm_cols & kfGlmColRef;
    const uint32_t alt1_col = glm_cols & kfGlmColAlt1;
    const uint32_t alt_col = glm_cols & kfGlmColAlt;
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
    cswritep = strcpya_k(cswritep, "\tA1");
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
        cswritep = strcpya_k(cswritep, "\tLOG10_P");
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
    uint32_t next_print_variant_idx = variant_ct / 100;
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
          reterr = ReadLocalCovarBlock(common, local_sample_uidx_order, local_variant_include, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, local_covar_txsp, &local_line_idx, &local_xy, ctx->local_covars_vcmaj_f[parity], nullptr, local_sample_idx_order);
        } else {
          float* prev_local_covar_row_f = nullptr;
          if (variant_idx) {
            prev_local_covar_row_f = &(ctx->local_covars_vcmaj_f[1 - parity][S_CAST(uintptr_t, read_block_size - 1) * max_sample_ct * local_covar_ct]);
          }
          reterr = ReadRfmix2Block(common, variant_bps, local_sample_uidx_order, prev_local_covar_row_f, nullptr, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, glm_info_ptr->local_chrom_col, glm_info_ptr->local_bp_col, glm_info_ptr->local_first_covar_col, local_covar_txsp, &local_line_iter, &local_line_idx, &local_prev_chr_code, &local_chr_code, &local_bp, &local_skip_chr, ctx->local_covars_vcmaj_f[parity], nullptr, local_sample_idx_order);
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
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
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
