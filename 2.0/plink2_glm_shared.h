#ifndef __PLINK2_GLM_SHARED_H__
#define __PLINK2_GLM_SHARED_H__

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


#include "plink2_common.h"
#include "plink2_matrix.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfGlm0,
  kfGlmZs = (1 << 0),

  kfGlmOmitRef = (1 << 1),

  // mutually exclusive
  kfGlmSex = (1 << 2),
  kfGlmNoXSex = (1 << 3),

  kfGlmLog10 = (1 << 4),

  // mutually exclusive
  kfGlmGenotypic = (1 << 5),
  kfGlmHethom = (1 << 6),
  kfGlmDominant = (1 << 7),
  kfGlmRecessive = (1 << 8),
  kfGlmHetonly = (1 << 9),

  kfGlmInteraction = (1 << 10),
  kfGlmHideCovar = (1 << 11),
  kfGlmIntercept = (1 << 12),
  kfGlmSkipInvalidPheno = (1 << 13),
  kfGlmNoFirth = (1 << 14),
  kfGlmFirth = (1 << 15),
  kfGlmPerm = (1 << 16),
  kfGlmPermCount = (1 << 17),
  kfGlmConditionDominant = (1 << 18),
  kfGlmConditionRecessive = (1 << 19),
  kfGlmConditionMultiallelic = (1 << 20),
  kfGlmLocalOmitLast = (1 << 21),
  kfGlmTestsAll = (1 << 22),
  kfGlmPhenoIds = (1 << 23),
  kfGlmLocalHaps = (1 << 24),
  kfGlmLocalCats1based = (1 << 25),
  kfGlmFirthResidualize = (1 << 26),
  kfGlmCcResidualize = (1 << 27),
  kfGlmQtResidualize = (1 << 28),
  kfGlmResidualizeMask = (kfGlmFirthResidualize | kfGlmCcResidualize | kfGlmQtResidualize),
  kfGlmSinglePrecCc = (1 << 29),
  kfGlmAllowNoCovars = (1 << 30)
FLAGSET_DEF_END(GlmFlags);

FLAGSET_DEF_START()
  kfGlmCol0,
  kfGlmColChrom = (1 << 0),
  kfGlmColPos = (1 << 1),
  kfGlmColRef = (1 << 2),
  kfGlmColAlt1 = (1 << 3),
  kfGlmColAlt = (1 << 4),
  kfGlmColMaybeprovref = (1 << 5),
  kfGlmColProvref = (1 << 6),
  kfGlmColOmitted = (1 << 7),
  kfGlmColA1count = (1 << 8),
  kfGlmColTotallele = (1 << 9),
  kfGlmColA1countcc = (1 << 10),
  kfGlmColTotallelecc = (1 << 11),
  kfGlmColGcountcc = (1 << 12),
  kfGlmColA1freq = (1 << 13),
  kfGlmColA1freqcc = (1 << 14),
  kfGlmColMachR2 = (1 << 15),
  kfGlmColFirthYn = (1 << 16),
  kfGlmColTest = (1 << 17),
  kfGlmColNobs = (1 << 18),

  // if beta specified, ignore orbeta
  kfGlmColBeta = (1 << 19),
  kfGlmColOrbeta = (1 << 20),

  kfGlmColSe = (1 << 21),
  kfGlmColCi = (1 << 22),
  kfGlmColTz = (1 << 23),
  kfGlmColP = (1 << 24),
  kfGlmColErr = (1 << 25),

  kfGlmColAx = (1 << 26),
  kfGlmColDefault = (kfGlmColChrom | kfGlmColPos | kfGlmColRef | kfGlmColAlt | kfGlmColProvref | kfGlmColOmitted | kfGlmColA1freq | kfGlmColFirthYn | kfGlmColTest | kfGlmColNobs | kfGlmColOrbeta | kfGlmColSe | kfGlmColCi | kfGlmColTz | kfGlmColP | kfGlmColErr),
  kfGlmColGwasSsfReq = (kfGlmColChrom | kfGlmColPos | kfGlmColRef | kfGlmColAlt | kfGlmColOmitted | kfGlmColA1freq | kfGlmColTest | kfGlmColNobs | kfGlmColSe | kfGlmColCi | kfGlmColTz | kfGlmColP)
FLAGSET_DEF_END(GlmColFlags);

typedef struct GlmInfoStruct {
  NONCOPYABLE(GlmInfoStruct);
  GlmFlags flags;
  GlmColFlags cols;
  uint32_t mperm_ct;
  uint32_t local_cat_ct;
  uint32_t local_header_line_ct;
  uint32_t local_chrom_col;
  uint32_t local_bp_col;
  uint32_t local_first_covar_col;
  double max_corr;
  char* condition_varname;
  char* condition_list_fname;
  RangeList parameters_range_list;
  RangeList tests_range_list;
} GlmInfo;

// Useful precomputed values for linear and logistic regression, for variants
// with no missing genotypes.
typedef struct {
  double* xtx_image;  // (covar_ct + domdev_present + 2)^2, genotype cols empty
  double* covarx_dotprod_inv;  // (covar_ct + 1) x (covar_ct + 1), reflected
  double* corr_inv;  // covar_ct x covar_ct, reflected
  double* corr_image;  // (covar_ct + domdev_present_p1)^2, genotype cols empty
  double* corr_inv_sqrts;  // covar_ct x 1
  double* xt_y_image;  // (covar_ct + domdev_present + 2) x 1
} RegressionNmPrecomp;

typedef struct GlmCtxStruct {
  const uintptr_t* variant_include;
  const ChrInfo* cip;
  const uintptr_t* allele_idx_offsets;
  const AlleleCode* omitted_alleles;
  const uint32_t* subset_chr_fo_vidx_start;
  const uintptr_t* sample_include;
  const uintptr_t* sample_include_x;
  const uintptr_t* sample_include_y;
  uint32_t* sample_include_cumulative_popcounts;
  const uint32_t* sample_include_x_cumulative_popcounts;
  const uint32_t* sample_include_y_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed;
  uintptr_t* parameter_subset;
  uintptr_t* parameter_subset_x;
  uintptr_t* parameter_subset_y;
  uintptr_t* joint_test_params;
  uintptr_t* joint_test_params_x;
  uintptr_t* joint_test_params_y;
  uint32_t variant_ct;
  uint32_t sample_ct;
  uint32_t sample_ct_x;
  uint32_t sample_ct_y;
  uint32_t max_extra_allele_ct;
  uint32_t covar_ct;
  uint32_t local_covar_ct;
  uint32_t covar_ct_x;
  uint32_t covar_ct_y;
  uint32_t constraint_ct;
  uint32_t constraint_ct_x;
  uint32_t constraint_ct_y;
  uint32_t tests_flag;
  GlmFlags glm_flags;
  uint32_t is_xchr_model_1;
  double max_corr;
  double vif_thresh;
  uintptr_t max_reported_test_ct;

#ifndef NDEBUG
  // temporary debug
  const char* outname;
#endif

  RegressionNmPrecomp* nm_precomp;
  RegressionNmPrecomp* nm_precomp_x;
  RegressionNmPrecomp* nm_precomp_y;

  uint32_t cur_block_variant_ct;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** thread_mhc;
  uintptr_t** dosage_presents;
  Dosage** dosage_mains;
  uint32_t* read_variant_uidx_starts;

  unsigned char** workspace_bufs;

  double* block_beta_se;

  uint64_t err_info;
} GlmCtx;

// may want to try to identify a specific linear dependency in rank-deficient
// case
ENUM_U31_DEF_START()
  kGlmErrcodeNone,
  kGlmErrcodeSampleCtLtePredictorCt,
  kGlmErrcodeConstOmittedAllele,
  kGlmErrcodeConstAllele,
  kGlmErrcodeCorrTooHigh, // 2 predictor args
  kGlmErrcodeVifInfinite,
  kGlmErrcodeVifTooHigh, // 1 predictor arg
  kGlmErrcodeSeparation, // 1 allele arg
  kGlmErrcodeRankDeficient,
  kGlmErrcodeLogisticConvergeFail,
  kGlmErrcodeFirthConvergeFail,
  kGlmErrcodeInvalidResult,
  // no codes for logistic-unfinished and firth-unfinished for now since we
  // still report results there

  // only during initial scan
  kGlmErrcodeUnstableScale
ENUM_U31_DEF_END(GlmErrcode);

#if __cplusplus >= 201103L
// see IntErr in plink2_base.h
struct GlmErr {
  GlmErr() {}

  GlmErr(uint64_t source) : value_(source) {}

  explicit operator uint64_t() const {
    return value_;
  }

  explicit operator uint32_t() const {
    return static_cast<uint32_t>(value_);
  }

  explicit operator bool() const {
    return (value_ != 0);
  }

private:
  uint64_t value_;
};
#else
typedef uint64_t GlmErr;
#endif

// The error code, along with up to two predictor/allele-index arguments, must
// fit in 8 bytes for now.  Since kMaxPhenoCt is larger than 2^16, we need a
// custom encoding.
// Current encoding has GlmErrcode in bits 0..7, the first argument in bits
// 8..31, the second argument in bits 32..55, and the top 8 bits are guaranteed
// to be zero.
HEADER_INLINE GlmErr SetGlmErr0(GlmErrcode errcode) {
  return errcode;
}

HEADER_INLINE GlmErr SetGlmErr1(GlmErrcode errcode, uint32_t arg1) {
  return errcode | (arg1 << 8);
}

HEADER_INLINE GlmErr SetGlmErr2(GlmErrcode errcode, uint32_t arg1, uint32_t arg2) {
  return S_CAST(uint64_t, errcode) | (arg1 << 8) | (S_CAST(uint64_t, arg2) << 32);
}

HEADER_INLINE GlmErrcode GetGlmErrCode(GlmErr glm_err) {
  return S_CAST(GlmErrcode, S_CAST(uint64_t, glm_err) & 255);
}

HEADER_INLINE uint32_t GetGlmErrArg1(GlmErr glm_err) {
  return S_CAST(uint32_t, glm_err) >> 8;
}

HEADER_INLINE uint32_t GetGlmErrArg2(GlmErr glm_err) {
  return S_CAST(uint64_t, glm_err) >> 32;
}

char* AppendGlmErrstr(GlmErr glm_err, char* write_iter);

GlmErr CheckMaxCorrAndVif(const double* predictor_dotprods, uint32_t first_predictor_idx, uint32_t predictor_ct, uintptr_t sample_ct, double max_corr, double vif_thresh, double* dbl_2d_buf, double* corr_buf, double* inverse_corr_buf, MatrixInvertBuf1* matrix_invert_buf1);

GlmErr CheckMaxCorrAndVifNm(const double* predictor_dotprods, const double* corr_inv, uint32_t predictor_ct, uint32_t geno_pred_ct, double sample_ct_recip, double sample_ct_m1_recip, double max_corr, double vif_thresh, double* __restrict semicomputed_corr_matrix, double* __restrict semicomputed_inv_corr_sqrts, double* __restrict corr_row_buf, double* __restrict inverse_corr_diag, double* __restrict ainv_b_buf);

PglErr GlmFillAndTestCovars(const uintptr_t* sample_include, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, double* covar_dotprod, double* corr_buf, double* inverse_corr_buf, double* covars_cmaj, const char** cur_covar_names, GlmErr* glm_err_ptr);

BoolErr InitNmPrecomp(const double* covars_cmaj, const double* covar_dotprod, const double* corr_buf, const double* corr_inv_tri, uint32_t sample_ct, uint32_t is_qt, uintptr_t new_covar_ct, uintptr_t xtx_state, RegressionNmPrecomp* nm_precomp);

uint32_t DosageIsConstant(uint64_t dosage_sum, uint64_t dosage_ssq, uint32_t nm_sample_ct);

uint32_t GetBiallelicReportedTestCt(const uintptr_t* parameter_subset, GlmFlags glm_flags, uint32_t covar_ct, uint32_t tests_flag);

PglErr ReadLocalCovarBlock(const GlmCtx* common, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t cur_block_variant_ct, uint32_t local_sample_ct, uint32_t local_cat_ct, TextStream* local_covar_txsp, uint32_t* local_line_idx_ptr, uint32_t* local_xy_ptr, float* local_covars_vcmaj_f_iter, double* local_covars_vcmaj_d_iter, uint32_t* local_sample_idx_order);

PglErr ReadRfmix2Block(const GlmCtx* common, const uint32_t* variant_bps, const uint32_t* local_sample_uidx_order, const float* prev_local_covar_row_f, const double* prev_local_covar_row_d, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t cur_block_variant_ct, uint32_t local_sample_ct, uint32_t local_cat_ct, uint32_t local_chrom_col, uint32_t local_bp_col, uint32_t local_first_covar_col, TextStream* local_covar_txsp, const char** local_line_iterp, uint32_t* local_line_idx_ptr, uint32_t* local_prev_chr_code_ptr, uint32_t* local_chr_code_ptr, uint32_t* local_bp_ptr, uint32_t* local_skip_chr_ptr, float* local_covars_vcmaj_f_iter, double* local_covars_vcmaj_d_iter, uint32_t* local_sample_idx_order);

extern const double kSmallDoublePairs[32];

extern const double kSmallInvDoublePairs[32];

extern const double kSmallInvDoubles[4];

uint32_t GenoarrToDoublesRemoveMissing(const uintptr_t* genoarr, const double* __restrict table, uint32_t sample_ct, double* __restrict dst);

BoolErr LinearHypothesisChisq(const double* coef, const double* constraints_con_major, const double* cov_matrix, uintptr_t constraint_ct, uintptr_t predictor_ct, uintptr_t cov_stride, double* chisq_ptr, double* tmphxs_buf, double* h_transpose_buf, double* inner_buf, MatrixInvertBuf1* mi_buf, double* outer_buf);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_GLM_SHARED_H__
