#ifndef __PLINK2_MISC_H__
#define __PLINK2_MISC_H__

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


#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfRecoverVarIds0,
  kfRecoverVarIdsStrictBimOrder = (1 << 0),
  kfRecoverVarIdsRigid = (1 << 1),
  kfRecoverVarIdsForce = (1 << 2),
  kfRecoverVarIdsPartial = (1 << 3)
FLAGSET_DEF_END(RecoverVarIdsFlags);

FLAGSET_DEF_START()
  kfUpdateSex0,
  kfUpdateSexMale0 = (1 << 0)
FLAGSET_DEF_END(UpdateSexFlags);

FLAGSET_DEF_START()
  kfPhenoTransform0,
  kfPhenoTransformSplitCat = (1 << 0),
  kfPhenoTransformSplitCatOmitMost = (1 << 1),
  kfPhenoTransformSplitCatOmitLast = (1 << 2),
  kfPhenoTransformSplitCatCovar01 = (1 << 3),
  kfPhenoTransformVstdCovar = (1 << 4),
  kfPhenoTransformVstdAll = (1 << 5),
  kfPhenoTransformQuantnormPheno = (1 << 6),
  kfPhenoTransformQuantnormCovar = (1 << 7),
  kfPhenoTransformQuantnormAll = (1 << 8)
FLAGSET_DEF_END(PhenoTransformFlags);

FLAGSET_DEF_START()
  kfWriteCovar0,
  kfWriteCovarColMaybefid = (1 << 0),
  kfWriteCovarColFid = (1 << 1),
  kfWriteCovarColMaybesid = (1 << 2),
  kfWriteCovarColSid = (1 << 3),
  kfWriteCovarColMaybeparents = (1 << 4),
  kfWriteCovarColParents = (1 << 5),
  kfWriteCovarColSex = (1 << 6),
  kfWriteCovarColPheno1 = (1 << 7),
  kfWriteCovarColPhenos = (1 << 8),
  kfWriteCovarColDefault = (kfWriteCovarColMaybefid | kfWriteCovarColMaybesid),
  kfWriteCovarColAll = ((kfWriteCovarColPhenos * 2) - kfWriteCovarColMaybefid)
FLAGSET_DEF_END(WriteCovarFlags);

FLAGSET_DEF_START()
  kfAlleleFreq0,
  kfAlleleFreqZs = (1 << 0),
  kfAlleleFreqCounts = (1 << 1),
  kfAlleleFreqBinsRefFname = (1 << 2),
  kfAlleleFreqBinsAlt1Fname = (1 << 3),
  kfAlleleFreqBinsOnly = (1 << 4),

  kfAlleleFreqColChrom = (1 << 5),
  kfAlleleFreqColPos = (1 << 6),
  kfAlleleFreqColRef = (1 << 7),
  kfAlleleFreqColAlt1 = (1 << 8),
  kfAlleleFreqColAlt = (1 << 9),
  kfAlleleFreqColReffreq = (1 << 10),
  kfAlleleFreqColAlt1freq = (1 << 11),
  kfAlleleFreqColAltfreq = (1 << 12),
  kfAlleleFreqColFreq = (1 << 13),
  kfAlleleFreqColEq = (1 << 14),
  kfAlleleFreqColEqz = (1 << 15),
  kfAlleleFreqColAlteq = (1 << 16),
  kfAlleleFreqColAlteqz = (1 << 17),
  kfAlleleFreqColNumeq = (1 << 18),
  kfAlleleFreqColAltnumeq = (1 << 19),
  kfAlleleFreqColMachR2 = (1 << 20),
  kfAlleleFreqColMinimac3R2 = (1 << 21),
  kfAlleleFreqColNobs = (1 << 22),
  kfAlleleFreqColDefault = (kfAlleleFreqColChrom | kfAlleleFreqColRef | kfAlleleFreqColAlt | kfAlleleFreqColAltfreq | kfAlleleFreqColNobs),
  kfAlleleFreqColAll = ((kfAlleleFreqColNobs * 2) - kfAlleleFreqColChrom),
  // only mutual exclusion is altfreq/freq/eq/eqz/alteq/alteqz/numeq/numeqz
  // don't force alt1freq/altfreq mutual exclusion since the former plays a bit
  // better with shell scripts
  // alt+alteqz is a bit silly, but I won't bother prohibiting it
  kfAlleleFreqColMutex = ((kfAlleleFreqColAltnumeq * 2) - kfAlleleFreqColAltfreq)
FLAGSET_DEF_END(FreqRptFlags);

// Will need to split this into multiple flagsets if we add any more options...
FLAGSET_DEF_START()
  kfMissingRpt0,
  kfMissingRptZs = (1 << 0),
  kfMissingRptSampleOnly = (1 << 1),
  kfMissingRptVariantOnly = (1 << 2),

  kfMissingRptScolMaybefid = (1 << 3),
  kfMissingRptScolFid = (1 << 4),
  kfMissingRptScolMaybesid = (1 << 5),
  kfMissingRptScolSid = (1 << 6),
  kfMissingRptScolMisspheno1 = (1 << 7),
  kfMissingRptScolMissphenos = (1 << 8),
  kfMissingRptScolNmissDosage = (1 << 9),
  kfMissingRptScolNmiss = (1 << 10),
  kfMissingRptScolNmissHh = (1 << 11),
  kfMissingRptScolHethap = (1 << 12),
  kfMissingRptScolNobs = (1 << 13),
  kfMissingRptScolFmissDosage = (1 << 14),
  kfMissingRptScolFmiss = (1 << 15),
  kfMissingRptScolFmissHh = (1 << 16),
  kfMissingRptScolDefault = (kfMissingRptScolMaybefid | kfMissingRptScolMaybesid | kfMissingRptScolMissphenos | kfMissingRptScolNmiss | kfMissingRptScolNobs | kfMissingRptScolFmiss),
  kfMissingRptScolAll = ((kfMissingRptScolFmissHh * 2) - kfMissingRptScolMaybesid),

  kfMissingRptVcolChrom = (1 << 17),
  kfMissingRptVcolPos = (1 << 18),
  kfMissingRptVcolRef = (1 << 19),
  kfMissingRptVcolAlt1 = (1 << 20),
  kfMissingRptVcolAlt = (1 << 21),
  kfMissingRptVcolNmissDosage = (1 << 22),
  kfMissingRptVcolNmiss = (1 << 23),
  kfMissingRptVcolNmissHh = (1 << 24),
  kfMissingRptVcolHethap = (1 << 25),
  kfMissingRptVcolNobs = (1 << 26),
  kfMissingRptVcolFmissDosage = (1 << 27),
  kfMissingRptVcolFmiss = (1 << 28),
  kfMissingRptVcolFmissHh = (1 << 29),
  kfMissingRptVcolFhethap = (1 << 30),
  kfMissingRptVcolDefault = (kfMissingRptVcolChrom | kfMissingRptVcolNmiss | kfMissingRptVcolNobs | kfMissingRptVcolFmiss),
  kfMissingRptVcolAll = ((kfMissingRptVcolFhethap * 2U) - kfMissingRptVcolChrom)
FLAGSET_DEF_END(MissingRptFlags);

FLAGSET_DEF_START()
  kfGenoCounts0,
  kfGenoCountsZs = (1 << 0),

  kfGenoCountsColChrom = (1 << 1),
  kfGenoCountsColPos = (1 << 2),
  kfGenoCountsColRef = (1 << 3),
  kfGenoCountsColAlt1 = (1 << 4),
  kfGenoCountsColAlt = (1 << 5),
  kfGenoCountsColHomref = (1 << 6),
  kfGenoCountsColRefalt1 = (1 << 7),
  kfGenoCountsColRefalt = (1 << 8),
  kfGenoCountsColHomalt1 = (1 << 9),
  kfGenoCountsColAltxy = (1 << 10),
  kfGenoCountsColXy = (1 << 11),
  kfGenoCountsColHapref = (1 << 12),
  kfGenoCountsColHapalt1 = (1 << 13),
  kfGenoCountsColHapalt = (1 << 14),
  kfGenoCountsColHap = (1 << 15),
  kfGenoCountsColNumeq = (1 << 16),
  kfGenoCountsColMissing = (1 << 17),
  kfGenoCountsColNobs = (1 << 18),
  kfGenoCountsColDefault = (kfGenoCountsColChrom | kfGenoCountsColRef | kfGenoCountsColAlt | kfGenoCountsColHomref | kfGenoCountsColRefalt | kfGenoCountsColAltxy | kfGenoCountsColHapref | kfGenoCountsColHapalt | kfGenoCountsColMissing),
  kfGenoCountsColAll = ((kfGenoCountsColNobs * 2) - kfGenoCountsColChrom),

  kfGenoCountsColPairex = (kfGenoCountsColHapalt | kfGenoCountsColHap),
  kfGenoCountsColMutex = (kfGenoCountsColAltxy | kfGenoCountsColXy | kfGenoCountsColNumeq)
FLAGSET_DEF_END(GenoCountsFlags);

FLAGSET_DEF_START()
  kfHardy0,
  kfHardyZs = (1 << 0),
  kfHardyMidp = (1 << 1),
  kfHardyRedundant = (1 << 2),

  kfHardyColChrom = (1 << 3),
  kfHardyColPos = (1 << 4),
  kfHardyColRef = (1 << 5),
  kfHardyColAlt1 = (1 << 6),
  kfHardyColAlt = (1 << 7),
  kfHardyColAx = (1 << 8),
  kfHardyColGcounts = (1 << 9),
  kfHardyColGcount1col = (1 << 10),
  kfHardyColHetfreq = (1 << 11),
  kfHardyColSexaf = (1 << 12),
  kfHardyColFemalep = (1 << 13),
  kfHardyColP = (1 << 14),
  kfHardyColDefault = (kfHardyColChrom | kfHardyColAx | kfHardyColGcounts | kfHardyColHetfreq | kfHardyColSexaf | kfHardyColP),
  kfHardyColAll = ((kfHardyColP * 2) - kfHardyColChrom)
FLAGSET_DEF_END(HardyFlags);

FLAGSET_DEF_START()
  kfHet0,
  kfHetZs = (1 << 0),
  kfHetSmallSample = (1 << 1),

  kfHetColMaybefid = (1 << 3),
  kfHetColFid = (1 << 4),
  kfHetColMaybesid = (1 << 5),
  kfHetColSid = (1 << 6),
  kfHetColHom = (1 << 7),
  kfHetColHet = (1 << 8),
  kfHetColNobs = (1 << 9),
  kfHetColF = (1 << 10),
  kfHetColDefault = (kfHetColMaybefid | kfHetColMaybesid | kfHetColHom | kfHetColNobs | kfHetColF),
  kfHetColAll = ((kfHetColF * 2) - kfHetColMaybefid)
FLAGSET_DEF_END(HetFlags);

typedef struct UpdateSexStruct {
  NONCOPYABLE(UpdateSexStruct);
  UpdateSexFlags flags;
  uint32_t col_num;
  char* fname;
} UpdateSexInfo;

FLAGSET_DEF_START()
  kfSampleCounts0,
  kfSampleCountsZs = (1 << 0),

  kfSampleCountsColMaybefid = (1 << 1),
  kfSampleCountsColFid = (1 << 2),
  kfSampleCountsColMaybesid = (1 << 3),
  kfSampleCountsColSid = (1 << 4),
  kfSampleCountsColSex = (1 << 5),
  kfSampleCountsColHom = (1 << 6),
  kfSampleCountsColHomref = (1 << 7),
  kfSampleCountsColHomalt = (1 << 8),
  kfSampleCountsColHomaltSnp = (1 << 9),
  kfSampleCountsColHet = (1 << 10),
  kfSampleCountsColRefalt = (1 << 11),
  kfSampleCountsColHet2alt = (1 << 12),

  kfSampleCountsColHetSnp = (1 << 13),
  kfSampleCountsColDiploidTs = (1 << 14),
  kfSampleCountsColTs = (1 << 15),
  kfSampleCountsColDiploidTv = (1 << 16),
  kfSampleCountsColTv = (1 << 17),
  kfSampleCountsColDiploidNonsnpsymb = (1 << 18),
  kfSampleCountsColNonsnpsymb = (1 << 19),
  kfSampleCountsColSymbolic = (1 << 20),
  kfSampleCountsColNonsnp = (1 << 21),

  kfSampleCountsColDiploidSingle = (1 << 22),
  kfSampleCountsColSingle = (1 << 23),
  kfSampleCountsColHaprefWithFemaleY = (1 << 24),
  kfSampleCountsColHapref = (1 << 25),
  kfSampleCountsColHapaltWithFemaleY = (1 << 26),
  kfSampleCountsColHapalt = (1 << 27),
  kfSampleCountsColMissingWithFemaleY = (1 << 28),
  kfSampleCountsColMissing = (1 << 29),
  kfSampleCountsColDefault = (kfSampleCountsColMaybefid | kfSampleCountsColMaybesid | kfSampleCountsColHomref | kfSampleCountsColHomaltSnp | kfSampleCountsColHetSnp | kfSampleCountsColDiploidTs | kfSampleCountsColDiploidTv | kfSampleCountsColDiploidNonsnpsymb | kfSampleCountsColDiploidSingle | kfSampleCountsColHaprefWithFemaleY | kfSampleCountsColHapaltWithFemaleY | kfSampleCountsColMissingWithFemaleY),
  kfSampleCountsColAll = ((kfSampleCountsColMissing * 2) - kfSampleCountsColMaybefid)
FLAGSET_DEF_END(SampleCountsFlags);

FLAGSET_DEF_START()
  kfSdiff0,
  kfSdiffIncludeMissing = (1 << 0),
  kfSdiffPairwise = (1 << 1),
  kfSdiffCountsOnly = (1 << 2),
  kfSdiffOneBase = (1 << 3),
  // no need for kfSdiffFile since that's synonymous with other_id_ct == 0
  kfSdiffZs = (1 << 4),

  kfSdiffColChrom = (1 << 5),
  kfSdiffColPos = (1 << 6),
  kfSdiffColRef = (1 << 7),
  kfSdiffColAlt = (1 << 8),
  kfSdiffColMaybefid = (1 << 9),
  kfSdiffColFid = (1 << 10),
  kfSdiffColId = (1 << 11),
  kfSdiffColMaybesid = (1 << 12),
  kfSdiffColSid = (1 << 13),
  kfSdiffColGeno = (1 << 14),
  kfSdiffColPairwiseDefault = (kfSdiffColChrom | kfSdiffColPos | kfSdiffColRef | kfSdiffColAlt | kfSdiffColGeno),
  kfSdiffColDefault = (kfSdiffColPairwiseDefault | kfSdiffColMaybefid | kfSdiffColId | kfSdiffColMaybesid),
  kfSdiffColAll = ((kfSdiffColGeno * 2) - kfSdiffColChrom),

  kfSdiffCountsColMaybefid = (1 << 15),
  kfSdiffCountsColFid = (1 << 16),
  kfSdiffCountsColMaybesid = (1 << 17),
  kfSdiffCountsColSid = (1 << 18),
  kfSdiffCountsColNobs = (1 << 19),
  kfSdiffCountsColNobsIbs = (1 << 20),
  kfSdiffCountsColIbs0 = (1 << 21),
  kfSdiffCountsColIbs1 = (1 << 22),
  kfSdiffCountsColIbs2 = (1 << 23),
  kfSdiffCountsIbsNeeded = (kfSdiffCountsColNobsIbs | kfSdiffCountsColIbs0 | kfSdiffCountsColIbs1 | kfSdiffCountsColIbs2),
  kfSdiffCountsColHalfmiss = (1 << 24),
  kfSdiffCountsColDiff = (1 << 25),
  kfSdiffCountsColDefault = (kfSdiffCountsColMaybefid | kfSdiffCountsColMaybesid | kfSdiffCountsColNobs | kfSdiffCountsColHalfmiss | kfSdiffCountsColDiff),
  kfSdiffCountsColAll = ((kfSdiffCountsColDiff * 2) - kfSdiffCountsColMaybefid)
FLAGSET_DEF_END(SdiffFlags);

typedef struct SdiffStruct {
  NONCOPYABLE(SdiffStruct);
  SdiffFlags flags;
  Dosage dosage_hap_tol; // missing value when 'dosage' not specified
  char fname_id_delim;
  uint32_t other_id_ct;
  char* first_id_or_fname;
  char* other_ids_flattened;
} SdiffInfo;

FLAGSET_DEF_START()
  kfFst0,
  kfFstMethodWc = (1 << 0),
  kfFstReportVariants = (1 << 1),
  kfFstZs = (1 << 2),
  kfFstOneBasePop = (1 << 3),
  kfFstExplicitPopIds = (1 << 4),
  kfFstPopPairFile = (1 << 5),

  kfFstColNobs = (1 << 6),
  kfFstColSe = (1 << 7),
  kfFstColDefault = kfFstColSe,
  kfFstColAll = ((kfFstColSe * 2) - kfFstColNobs),

  kfFstVcolChrom = (1 << 8),
  kfFstVcolPos = (1 << 9),
  kfFstVcolRef = (1 << 10),
  kfFstVcolAlt1 = (1 << 11),
  kfFstVcolAlt = (1 << 12),
  kfFstVcolNobs = (1 << 13),
  kfFstVcolFstfrac = (1 << 14),
  kfFstVcolFst = (1 << 15),
  kfFstVcolDefault = (kfFstVcolChrom | kfFstVcolPos | kfFstVcolNobs | kfFstVcolFst),
  kfFstVcolAll = ((kfFstVcolFst * 2) - kfFstVcolChrom),
FLAGSET_DEF_END(FstFlags);

typedef struct FstInfoStruct {
  NONCOPYABLE(FstInfoStruct);
  FstFlags flags;
  uint32_t blocksize;
  uint32_t other_id_ct;
  char* pheno_name;
  char* first_id_or_fname;
  char* other_ids_flattened;
} FstInfo;

void InitUpdateSex(UpdateSexInfo* update_sex_info_ptr);

void CleanupUpdateSex(UpdateSexInfo* update_sex_info_ptr);

void InitSdiff(SdiffInfo* sdiff_info_ptr);

void CleanupSdiff(SdiffInfo* sdiff_info_ptr);

void InitFst(FstInfo* fst_info_ptr);

void CleanupFst(FstInfo* fst_info_ptr);

PglErr UpdateVarBps(const ChrInfo* cip, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const TwoColParams* params, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uint32_t htable_size, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* __restrict variant_bps, uint32_t* __restrict variant_ct_ptr, UnsortedVar* vpos_sortstatusp);

PglErr UpdateVarNames(const uintptr_t* variant_include, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const TwoColParams* params, uint32_t raw_variant_ct, uint32_t htable_size, uint32_t max_thread_ct, char** variant_ids, uint32_t* max_variant_id_slen_ptr);

PglErr UpdateVarAlleles(const char* fname, const uintptr_t* variant_include, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uint32_t htable_size, uint32_t max_thread_ct, char** allele_storage_mutable, uint32_t* max_allele_slen_ptr, char* outname, char* outname_end);

PglErr RecoverVarIds(const char* fname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const char* missing_varid, uint32_t raw_variant_ct, uint32_t variant_ct, RecoverVarIdsFlags flags, uint32_t max_thread_ct, char** variant_ids, uint32_t* max_variant_id_slen_ptr, char* outname, char* outname_end);

PglErr Plink1ClusterImport(const char* within_fname, const char* catpheno_name, const char* family_missing_catname, const uintptr_t* sample_include, const char* sample_ids, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t mwithin_val, uint32_t max_thread_ct, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr);

// These functions return kPglRetEof on empty files.
PglErr PrescanSampleIds(const char* fname, SampleIdInfo* siip);

PglErr PrescanParentalIds(const char* fname, uint32_t max_thread_ct, ParentalIdInfo* parental_id_infop);


PglErr UpdateSampleIds(const char* fname, const uintptr_t* sample_include, uint32_t raw_sample_ct, uintptr_t sample_ct, SampleIdInfo* siip);

PglErr UpdateSampleParents(const char* fname, const SampleIdInfo* siip, const uintptr_t* sample_include, uint32_t raw_sample_ct, uintptr_t sample_ct, uint32_t max_thread_ct, ParentalIdInfo* parental_id_infop, uintptr_t* founder_info);

PglErr UpdateSampleSexes(const uintptr_t* sample_include, const SampleIdInfo* siip, const UpdateSexInfo* update_sex_info_ptr, uint32_t raw_sample_ct, uintptr_t sample_ct, uint32_t max_thread_ct, uintptr_t* sex_nm, uintptr_t* sex_male);

PglErr SplitCatPheno(const char* split_cat_phenonames_flattened, const uintptr_t* sample_include, uint32_t raw_sample_ct, PhenoTransformFlags pheno_transform_flags, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr, PhenoCol** covar_cols_ptr, char** covar_names_ptr, uint32_t* covar_ct_ptr, uintptr_t* max_covar_name_blen_ptr);

PglErr PhenoVarianceStandardize(const char* vstd_flattened, const uintptr_t* sample_include, const char* pheno_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uint32_t is_covar_flag, PhenoCol* pheno_cols);

PglErr PhenoQuantileNormalize(const char* quantnorm_flattened, const uintptr_t* sample_include, const char* pheno_names, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uint32_t is_subset_flag, PhenoCol* pheno_cols);

PglErr WriteAlleleFreqs(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint64_t* founder_allele_ddosages, const double* imp_r2_vals, const char* ref_binstr, const char* alt1_binstr, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, FreqRptFlags freq_rpt_flags, uint32_t max_thread_ct, uint32_t nonfounders, char* outname, char* outname_end);

PglErr WriteGenoCounts(const uintptr_t* sample_include, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts), uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t variant_ct, uint32_t x_start, uint32_t max_allele_slen, GenoCountsFlags geno_counts_flags, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr WriteMissingnessReports(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* sample_missing_hc_cts, const uint32_t* sample_missing_dosage_cts, const uint32_t* sample_hethap_cts, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint32_t* variant_missing_hc_cts, const uint32_t* variant_missing_dosage_cts, const uint32_t* variant_hethap_cts, uint32_t sample_ct, uint32_t male_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t variant_ct, uintptr_t max_allele_slen, uint32_t first_hap_uidx, MissingRptFlags missing_rpt_flags, uint32_t max_thread_ct, char* outname, char* outname_end);

PglErr GetMultiallelicMarginalCounts(const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_geno_cts), uint32_t raw_sample_ct, uint32_t autosomal_variant_ct, uint32_t autosomal_xallele_ct, uint32_t hwe_x_ct, uint32_t x_xallele_ct, PgenReader* simple_pgrp, STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts), STD_ARRAY_PTR_DECL(uint32_t, 2, autosomal_xgeno_cts), STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts));

PglErr ComputeHweXPvals(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts), uint32_t x_start, uint32_t hwe_x_ct, uintptr_t x_xallele_ct, uint32_t hwe_midp, uint32_t calc_thread_ct, double** hwe_x_pvals_ptr);

PglErr HardyReport(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, autosomal_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_nosex_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts), const double* hwe_x_pvals, uint32_t variant_ct, uint32_t hwe_x_ct, uint32_t max_allele_slen, double output_min_ln, HardyFlags hardy_flags, uint32_t max_thread_ct, uint32_t nonfounders, char* outname, char* outname_end);

PglErr SampleCounts(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, SampleCountsFlags flags, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end);

PglErr Sdiff(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const SdiffInfo* sdip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t variant_ct, uint32_t iid_sid, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr WriteSnplist(const uintptr_t* variant_include, const char* const* variant_ids, uint32_t variant_ct, uint32_t output_zst, uint32_t allow_dups, uint32_t max_thread_ct, char* outname, char* outname_end);

PglErr WriteCovar(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, WriteCovarFlags write_covar_flags, char* outname, char* outname_end);

PglErr HetReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const double* allele_freqs, const uintptr_t* founder_info, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_allele_ct, HetFlags flags, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end);

PglErr FstReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const FstInfo* fst_infop, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_MISC_H__
