#ifndef __PLINK2_MATRIX_CALC_H__
#define __PLINK2_MATRIX_CALC_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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

  kfKingNoIdheader = (1 << 8),

  kfKingColId = (1 << 9),
  kfKingColMaybesid = (1 << 10),
  kfKingColSid = (1 << 11),
  kfKingColNsnp = (1 << 12),
  kfKingColHethet = (1 << 13),
  kfKingColIbs0 = (1 << 14),
  kfKingColIbs1 = (1 << 15),
  kfKingColKinship = (1 << 16),
  kfKingColDefault = (kfKingColId | kfKingColMaybesid | kfKingColNsnp | kfKingColHethet | kfKingColIbs0 | kfKingColKinship),
  kfKingColAll = ((kfKingColKinship * 2) - kfKingColId)
FLAGSET_DEF_END(KingFlags);

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
  kfGrmListNoGz = (1 << 6),
  kfGrmListZs = (1 << 7),
  kfGrmListmask = (kfGrmListNoGz | kfGrmListZs),
  kfGrmBin = (1 << 8),

  kfGrmMeanimpute = (1 << 9),
  kfGrmCov = (1 << 10),
  kfGrmNoIdheader = (1 << 11)
FLAGSET_DEF_END(GrmFlags);

FLAGSET_DEF_START()
  kfPca0,
  kfPcaApprox = (1 << 0),
  kfPcaMeanimpute = (1 << 1),
  kfPcaSid = (1 << 2),
  kfPcaVarWts = (1 << 3),
  kfPcaVarZs = (1 << 4),

  kfPcaVcolChrom = (1 << 5),
  kfPcaVcolPos = (1 << 6),
  kfPcaVcolRef = (1 << 7),
  kfPcaVcolAlt1 = (1 << 8),
  kfPcaVcolAlt = (1 << 9),
  kfPcaVcolMaj = (1 << 10),
  kfPcaVcolNonmaj = (1 << 11),
  kfPcaVcolDefault = (kfPcaVcolChrom | kfPcaVcolMaj | kfPcaVcolNonmaj),
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
  kfScoreListVariants = (1 << 9),
  kfScoreListVariantsZs = (1 << 10),

  kfScoreColMaybesid = (1 << 11),
  kfScoreColSid = (1 << 12),
  kfScoreColPheno1 = (1 << 13),
  kfScoreColPhenos = (1 << 14),
  kfScoreColNmissAllele = (1 << 15),
  kfScoreColDenom = (1 << 16),
  kfScoreColDosageSum = (1 << 17),
  kfScoreColScoreAvgs = (1 << 18),
  kfScoreColScoreSums = (1 << 19),
  kfScoreColDefault = (kfScoreColMaybesid | kfScoreColPhenos | kfScoreColNmissAllele | kfScoreColDosageSum | kfScoreColScoreAvgs),
  kfScoreColAll = ((kfScoreColScoreSums * 2) - kfScoreColMaybesid)
FLAGSET_DEF_END(ScoreFlags);

typedef struct score_info_struct {
  ScoreFlags flags;
  uint32_t varid_col_p1;
  uint32_t allele_col_p1;
  char* input_fname;
  RangeList input_col_idx_range_list;
} ScoreInfo;

void InitScore(ScoreInfo* score_info_ptr);

void CleanupScore(ScoreInfo* score_info_ptr);

PglErr KingCutoffBatch(const char* sample_ids, const char* sids, uint32_t raw_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, double king_cutoff, uintptr_t* sample_include, char* king_cutoff_fprefix, uint32_t* sample_ct_ptr);

PglErr CalcKing(const char* sample_ids, const char* sids, uintptr_t* variant_include, const ChrInfo* cip, uint32_t raw_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t raw_variant_ct, uint32_t variant_ct, double king_cutoff, double king_table_filter, KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, uintptr_t* sample_include, uint32_t* sample_ct_ptr, char* outname, char* outname_end);

PglErr CalcKingTableSubset(const uintptr_t* orig_sample_include, const char* sample_ids, const char* sids, uintptr_t* variant_include, const ChrInfo* cip, const char* subset_fname, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t raw_variant_ct, uint32_t variant_ct, double king_table_filter, double king_table_subset_thresh, KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr CalcGrm(const uintptr_t* orig_sample_include, const char* sample_ids, const char* sids, uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* variant_allele_idxs, const AltAlleleCt* maj_alleles, const double* allele_freqs, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t raw_variant_ct, uint32_t variant_ct, GrmFlags grm_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end, double** grm_ptr);

#ifndef NOLAPACK
PglErr CalcPca(const uintptr_t* sample_include, const char* sample_ids, const char* sids, uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const AltAlleleCt* maj_alleles, const double* allele_freqs, uint32_t raw_sample_ct, uintptr_t pca_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t pc_ct, PcaFlags pca_flags, uint32_t max_thread_ct, PgenReader* simple_pgrp, double* grm, char* outname, char* outname_end);
#endif

PglErr ScoreReport(const uintptr_t* sample_include, const char* sample_ids, const char* sids, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const double* allele_freqs, const ScoreInfo* score_info_ptr, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t xchr_model, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_MATRIX_CALC_H__
