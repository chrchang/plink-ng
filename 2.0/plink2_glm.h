#ifndef __PLINK2_GLM_H__
#define __PLINK2_GLM_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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


#include "plink2_adjust.h"

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

  kfGlmInteraction = (1 << 9),
  kfGlmHideCovar = (1 << 10),
  kfGlmIntercept = (1 << 11),
  kfGlmSkip = (1 << 12),
  kfGlmNoFirth = (1 << 13),
  kfGlmFirth = (1 << 14),
  kfGlmPerm = (1 << 15),
  kfGlmPermCount = (1 << 16),
  kfGlmConditionDominant = (1 << 17),
  kfGlmConditionRecessive = (1 << 18),
  kfGlmConditionMultiallelic = (1 << 19),
  kfGlmLocalOmitLast = (1 << 20),
  kfGlmTestsAll = (1 << 21),
  kfGlmPhenoIds = (1 << 22),
  kfGlmLocalHaps = (1 << 23),
  kfGlmLocalCats1based = (1 << 24)
FLAGSET_DEF_END(GlmFlags);

FLAGSET_DEF_START()
  kfGlmCol0,
  kfGlmColChrom = (1 << 0),
  kfGlmColPos = (1 << 1),
  kfGlmColRef = (1 << 2),
  kfGlmColAlt1 = (1 << 3),
  kfGlmColAlt = (1 << 4),
  kfGlmColAx = (1 << 5),
  kfGlmColA1count = (1 << 6),
  kfGlmColTotallele = (1 << 7),
  kfGlmColA1countcc = (1 << 8),
  kfGlmColTotallelecc = (1 << 9),
  kfGlmColGcountcc = (1 << 10),
  kfGlmColA1freq = (1 << 11),
  kfGlmColA1freqcc = (1 << 12),
  kfGlmColMachR2 = (1 << 13),
  kfGlmColFirthYn = (1 << 14),
  kfGlmColTest = (1 << 15),
  kfGlmColNobs = (1 << 16),

  // if beta specified, ignore orbeta
  kfGlmColBeta = (1 << 17),
  kfGlmColOrbeta = (1 << 18),

  kfGlmColSe = (1 << 19),
  kfGlmColCi = (1 << 20),
  kfGlmColTz = (1 << 21),
  kfGlmColP = (1 << 22),
  kfGlmColErr = (1 << 23),
  kfGlmColDefault = (kfGlmColChrom | kfGlmColPos | kfGlmColRef | kfGlmColAlt | kfGlmColFirthYn | kfGlmColTest | kfGlmColNobs | kfGlmColOrbeta | kfGlmColSe | kfGlmColCi | kfGlmColTz | kfGlmColP | kfGlmColErr)
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

void InitGlm(GlmInfo* glm_info_ptr);

void CleanupGlm(GlmInfo* glm_info_ptr);

// for testing purposes
// plink2_matrix.h must be included in this file
// BoolErr LogisticRegression(const float* yy, const float* xx, uint32_t sample_ct, uint32_t predictor_ct, float* coef, uint32_t* is_unfinished_ptr, float* ll, float* pp, float* vv, float* hh, float* grad, float* dcoef) {

// BoolErr FirthRegression(const float* yy, const float* xx, uint32_t sample_ct, uint32_t predictor_ct, float* coef, uint32_t* is_unfinished_ptr, float* hh, double* half_inverted_buf, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, float* pp, float* vv, float* grad, float* dcoef, float* ww, float* tmpnxk_buf) {

PglErr GlmMain(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const AdjustInfo* adjust_info_ptr, const APerm* aperm_ptr, const char* local_covar_fname, const char* local_pvar_fname, const char* local_psam_fname, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t orig_covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t xchr_model, double ci_size, double vif_thresh, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_GLM_H__
