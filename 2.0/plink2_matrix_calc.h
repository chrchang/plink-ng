#ifndef __PLINK2_MATRIX_CALC_H__
#define __PLINK2_MATRIX_CALC_H__

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

#include "include/SFMT.h"
#include "include/pgenlib_misc.h"
#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "plink2_cmdline.h"
#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfKing0,
  kfKingMatrixZs = (1 << 0),
  kfKingMatrixBin = (1 << 1),
  kfKingMatrixBin4 = (1 << 2),
  kfKingMatrixEncodemask = (kfKingMatrixZs | kfKingMatrixBin | kfKingMatrixBin4),
  kfKingMatrixSq = (1 << 3),
  kfKingMatrixSq0 = (1 << 4),
  kfKingMatrixTri = (1 << 5),
  kfKingMatrixShapemask = (kfKingMatrixSq0 | kfKingMatrixSq | kfKingMatrixTri),

  kfKingTableZs = (1 << 6),
  kfKingCounts = (1 << 7),
  kfKingRelCheck = (1 << 8),
  kfKingConcordanceCheck = (1 << 9),
  kfKingCutoffTable = (1 << 10),
  kfKingTableRequireXor = (1 << 11),

  kfKingColMaybefid = (1 << 12),
  kfKingColFid = (1 << 13),
  kfKingColId = (1 << 14),
  kfKingColMaybesid = (1 << 15),
  kfKingColSid = (1 << 16),
  kfKingColNsnp = (1 << 17),
  kfKingColHethet = (1 << 18),
  kfKingColIbs0 = (1 << 19),
  kfKingColIbs1 = (1 << 20),
  kfKingColHamming = (1 << 21),
  kfKingColKinship = (1 << 22),
  kfKingColDefault = (kfKingColMaybefid | kfKingColId | kfKingColMaybesid | kfKingColNsnp | kfKingColHethet | kfKingColIbs0 | kfKingColKinship),
  kfKingColAll = ((kfKingColKinship * 2) - kfKingColMaybefid)
FLAGSET_DEF_END(KingFlags);

// --make-rel, --make-grm-list, --make-grm-bin, --make-grm-sparse
//   --make-rel iff kfGrmMatrixShapemask nonzero
//   --make-grm-list iff kfGrmList
//   --make-grm-bin iff kfGrmBin
//   --make-grm-sparse iff kfGrmSparse
//   (no GRM output possible in --pca case)
FLAGSET_DEF_START()
  kfGrm0,
  kfGrmMatrixZs = (1 << 0),
  kfGrmMatrixBin = (1 << 1),
  kfGrmMatrixBin4 = (1 << 2),
  kfGrmMatrixEncodemask = (kfGrmMatrixZs | kfGrmMatrixBin | kfGrmMatrixBin4),
  kfGrmMatrixSq = (1 << 3),
  kfGrmMatrixSq0 = (1 << 4),
  kfGrmMatrixTri = (1 << 5),
  kfGrmMatrixShapemask = (kfGrmMatrixSq0 | kfGrmMatrixSq | kfGrmMatrixTri),
  kfGrmList = (1 << 6),
  kfGrmListZs = (1 << 7),
  kfGrmBin = (1 << 8),
  kfGrmSparse = (1 << 9),
  kfGrmOutputMask = (kfGrmMatrixShapemask | kfGrmList | kfGrmBin | kfGrmSparse),

  kfGrmMeanimpute = (1 << 10),
  kfGrmCov = (1 << 11),
  kfGrmNoIdHeader = (1 << 12),
  kfGrmNoIdHeaderIidOnly = (1 << 13)
FLAGSET_DEF_END(GrmFlags);

FLAGSET_DEF_START()
  kfPca0,
  kfPcaApprox = (1 << 0),
  kfPcaMeanimpute = (1 << 1),
  kfPcaAlleleWts = (1 << 2),
  kfPcaBiallelicVarWts = (1 << 3),
  kfPcaVarZs = (1 << 4),

  kfPcaScolMaybefid = (1 << 5),
  kfPcaScolFid = (1 << 6),
  kfPcaScolMaybesid = (1 << 7),
  kfPcaScolSid = (1 << 8),
  kfPcaScolDefault = (kfPcaScolMaybefid | kfPcaScolMaybesid),
  kfPcaScolAll = ((kfPcaScolSid * 2) - kfPcaScolMaybefid),

  kfPcaVcolChrom = (1 << 9),
  kfPcaVcolPos = (1 << 10),
  kfPcaVcolRef = (1 << 11),
  kfPcaVcolAlt1 = (1 << 12),
  kfPcaVcolAlt = (1 << 13),
  kfPcaVcolMaybeprovref = (1 << 14),
  kfPcaVcolProvref = (1 << 15),
  kfPcaVcolAx = (1 << 16),
  kfPcaVcolMaj = (1 << 17),
  kfPcaVcolNonmaj = (1 << 18),
  kfPcaVcolDefaultA = (kfPcaVcolChrom | kfPcaVcolRef | kfPcaVcolAlt | kfPcaVcolMaybeprovref),
  kfPcaVcolDefaultB = (kfPcaVcolChrom | kfPcaVcolMaybeprovref | kfPcaVcolMaj | kfPcaVcolNonmaj),
  kfPcaVcolAll = ((kfPcaVcolNonmaj * 2) - kfPcaVcolChrom)
FLAGSET_DEF_END(PcaFlags);

FLAGSET_DEF_START()
  kfScore0,
  kfScoreHeaderIgnore = (1 << 0),
  kfScoreHeaderRead = (1 << 1),
  kfScoreNoMeanimpute = (1 << 2),
  kfScoreDominant = (1 << 3),
  kfScoreRecessive = (1 << 4),
  kfScoreCenter = (1 << 5),
  kfScoreVarianceStandardize = (1 << 6),
  kfScoreSe = (1 << 7),
  kfScoreZs = (1 << 8),
  kfScoreIgnoreDupIds = (1 << 9),
  kfScoreListVariants = (1 << 10),
  kfScoreListVariantsZs = (1 << 11),
  kfScoreMultiInput = (1 << 12),

  kfScoreQsrHeader = (1 << 13),
  kfScoreQsrMin = (1 << 14),

  kfScoreColMaybefid = (1 << 15),
  kfScoreColFid = (1 << 16),
  kfScoreColMaybesid = (1 << 17),
  kfScoreColSid = (1 << 18),
  kfScoreColPheno1 = (1 << 19),
  kfScoreColPhenos = (1 << 20),
  kfScoreColNallele = (1 << 21),
  kfScoreColDenom = (1 << 22),
  kfScoreColDosageSum = (1 << 23),
  kfScoreColScoreAvgs = (1 << 24),
  kfScoreColScoreSums = (1 << 25),
  kfScoreListColDefault = (kfScoreColMaybefid | kfScoreColMaybesid | kfScoreColPhenos | kfScoreColScoreAvgs),
  kfScoreColDefault = (kfScoreListColDefault | kfScoreColNallele | kfScoreColDosageSum),
  kfScoreColAll = ((kfScoreColScoreSums * 2) - kfScoreColMaybefid)
FLAGSET_DEF_END(ScoreFlags);

FLAGSET_DEF_START()
  kfVscore0,
  kfVscoreZs = (1 << 0),
  kfVscoreSinglePrec = (1 << 1),
  kfVscoreBin = (1 << 2),
  kfVscoreBin4 = (1 << 3),

  kfVscoreColChrom = (1 << 4),
  kfVscoreColPos = (1 << 5),
  kfVscoreColRef = (1 << 6),
  kfVscoreColAlt1 = (1 << 7),
  kfVscoreColAlt = (1 << 8),
  kfVscoreColMaybeprovref = (1 << 9),
  kfVscoreColProvref = (1 << 10),
  kfVscoreColAltfreq = (1 << 11),
  kfVscoreColNmiss = (1 << 12),
  kfVscoreColNobs = (1 << 13),
  kfVscoreColDefault = (kfVscoreColChrom | kfVscoreColPos | kfVscoreColRef | kfVscoreColAlt | kfVscoreColMaybeprovref),
  kfVscoreColAll = ((kfVscoreColNobs * 2) - kfVscoreColChrom)
FLAGSET_DEF_END(VscoreFlags);

typedef struct ScoreInfoStruct {
  NONCOPYABLE(ScoreInfoStruct);
  ScoreFlags flags;
  uint32_t varid_col_p1;
  uint32_t allele_col_p1;
  char* input_fname;
  RangeList input_col_idx_range_list;

  // --q-score-range
  char* qsr_range_fname;
  char* qsr_data_fname;
  uint32_t qsr_varid_col_p1;
  uint32_t qsr_val_col_p1;
} ScoreInfo;

FLAGSET_DEF_START()
  kfPhenoSvd0,
  kfPhenoSvdForce = (1 << 0),

  kfPhenoSvdScolMaybefid = (1 << 1),
  kfPhenoSvdScolFid = (1 << 2),
  kfPhenoSvdScolMaybesid = (1 << 3),
  kfPhenoSvdScolSid = (1 << 4),
  kfPhenoSvdScolDefault = (kfPhenoSvdScolMaybefid | kfPhenoSvdScolMaybesid),
  kfPhenoSvdScolAll = ((kfPhenoSvdScolSid * 2) - kfPhenoSvdScolMaybefid),

  kfPhenoSvdPcolId = (1 << 5),
  kfPhenoSvdPcolSv = (1 << 6),
  kfPhenoSvdPcolDefault = (kfPhenoSvdPcolId | kfPhenoSvdPcolSv),
  kfPhenoSvdPcolAll = ((kfPhenoSvdPcolSv * 2) - kfPhenoSvdPcolId)
FLAGSET_DEF_END(PhenoSvdFlags);

typedef struct PhenoSvdInfoStruct {
  PhenoSvdFlags flags;
  uint32_t ct;
  double min_variance_explained;
} PhenoSvdInfo;

void InitScore(ScoreInfo* score_info_ptr);

void CleanupScore(ScoreInfo* score_info_ptr);

void InitPhenoSvd(PhenoSvdInfo* pheno_svd_info_ptr);

CONSTI32(kMaxPc, 8000);

PglErr KingCutoffBatchBinary(const SampleIdInfo* siip, uint32_t raw_sample_ct, double king_cutoff, uintptr_t* sample_include, char* king_cutoff_fprefix, uint32_t* sample_ct_ptr);

PglErr KingCutoffBatchTable(const SampleIdInfo* siip, const char* kin0_fname, uint32_t raw_sample_ct, double king_cutoff, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

PglErr CalcKing(const SampleIdInfo* siip, const uintptr_t* variant_include_orig, const ChrInfo* cip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, double king_cutoff, double king_table_filter, KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, uintptr_t* sample_include, uint32_t* sample_ct_ptr, char* outname, char* outname_end);

ENUM_U31_DEF_START()
  kRcCheck0,
  kRcCheckRel,
  kRcCheckConcordance
ENUM_U31_DEF_END(RelConcordanceCheckMode);

PglErr CalcKingTableSubset(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const char* subset_fname, const char* require_fnames, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, double king_table_filter, double king_table_subset_thresh, RelConcordanceCheckMode rel_or_concordance_check, KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr CalcGrm(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const double* allele_freqs, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, GrmFlags grm_flags, double grm_sparse_cutoff, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end, double** grm_ptr);

#ifndef NOLAPACK
PglErr CalcPca(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const double* allele_freqs, uint32_t raw_sample_ct, uintptr_t pca_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t pc_ct, PcaFlags pca_flags, uint32_t max_thread_ct, PgenReader* simple_pgrp, sfmt_t* sfmtp, double* grm, char* outname, char* outname_end);
#endif

PglErr ScoreReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const double* allele_freqs, const ScoreInfo* score_info_ptr, const char* output_missing_pheno, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t nosex_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t xchr_model, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr Vscore(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_male, const double* allele_freqs, const char* in_fname, const RangeList* col_idx_range_listp, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t nosex_ct, uint32_t max_allele_slen, VscoreFlags flags, uint32_t xchr_model, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end);

#ifndef NOLAPACK
PglErr PhenoSvd(const PhenoSvdInfo* psip, const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t max_thread_ct, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr, PhenoCol* pheno_cols, char* outname, char* outname_end);
#endif

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_MATRIX_CALC_H__
