#ifndef __PLINK2_MERGE_H__
#define __PLINK2_MERGE_H__

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

#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

ENUM_U31_DEF_START()
  kPmergeListModeBfile,
  kPmergeListModeBpfile,
  kPmergeListModePfile,
  kPmergeListModePfileVzs
ENUM_U31_DEF_END(PmergeListMode);

FLAGSET_DEF_START()
  kfPmerge0,
  kfPmergeSampleInnerJoin = (1 << 0),
  kfPmergeVariantInnerJoin = (1 << 1),
  kfPmergePhenoInnerJoin = (1 << 2),
  kfPmergeSids = (1 << 3),
  kfPmergeMultiallelicsAlreadyJoined = (1 << 4),
  kfPmergeOutputVzs = (1 << 5)
FLAGSET_DEF_END(PmergeFlags);

ENUM_U31_DEF_START()
  kMergeModeNmMatch,
  kMergeModeNmFirst,
  kMergeModeFirst
ENUM_U31_DEF_END(MergeMode);

ENUM_U31_DEF_START()
  kMergePhenoModeNmMatch,
  kMergePhenoModeNmFirst,
  kMergePhenoModeFirst
ENUM_U31_DEF_END(MergePhenoMode);

ENUM_U31_DEF_START()
  kMergeXheaderModeErase,
  kMergeXheaderModeMatch,
  kMergeXheaderModeFirst
ENUM_U31_DEF_END(MergeXheaderMode);

ENUM_U31_DEF_START()
  kMergeQualModeErase,
  kMergeQualModeNmMatch,
  kMergeQualModeNmFirst,
  kMergeQualModeFirst,
  kMergeQualModeMin
ENUM_U31_DEF_END(MergeQualMode);

ENUM_U31_DEF_START()
  kMergeFilterModeErase,
  kMergeFilterModeNmMatch,
  kMergeFilterModeNmFirst,
  kMergeFilterModeFirst,
  kMergeFilterModeNonpassUnion
ENUM_U31_DEF_END(MergeFilterMode);

ENUM_U31_DEF_START()
  kMergeInfoCmModeErase,
  kMergeInfoCmModeNmMatch,
  kMergeInfoCmModeNmFirst,
  kMergeInfoCmModeFirst
ENUM_U31_DEF_END(MergeInfoCmMode);

ENUM_U31_DEF_START()
  kSort0,
  kSortNone,
  kSortNatural,
  kSortAscii,
  kSortFile
ENUM_U31_DEF_END(SortMode);

// --pgen-diff based here due to --merge-mode 6/7 history
FLAGSET_DEF_START()
  kfPgenDiff0,
  kfPgenDiffIncludeMissing = (1 << 0),
  kfPgenDiffZs = (1 << 1),
  kfPgenDiffColChrom = (1 << 2),
  kfPgenDiffColPos = (1 << 3),
  kfPgenDiffColId = (1 << 4),
  kfPgenDiffColRef = (1 << 5),
  kfPgenDiffColAlt = (1 << 6),
  kfPgenDiffColMaybeprovref = (1 << 7),
  kfPgenDiffColProvref = (1 << 8),
  kfPgenDiffColMaybefid = (1 << 9),
  kfPgenDiffColFid = (1 << 10),
  kfPgenDiffColMaybesid = (1 << 11),
  kfPgenDiffColSid = (1 << 12),
  kfPgenDiffColGeno = (1 << 13),
  kfPgenDiffColDefault = (kfPgenDiffColId | kfPgenDiffColMaybeprovref | kfPgenDiffColMaybefid | kfPgenDiffColMaybesid | kfPgenDiffColGeno),
  kfPgenDiffColAll = ((kfPgenDiffColGeno * 2) - kfPgenDiffColPos)
FLAGSET_DEF_END(PgenDiffFlags);

typedef struct PmergeStruct {
  NONCOPYABLE(PmergeStruct);
  PmergeFlags flags;
  PmergeListMode list_mode;
  MergeMode merge_mode;
  MergePhenoMode merge_parents_mode;
  MergePhenoMode merge_sex_mode;
  MergePhenoMode merge_pheno_mode;
  MergeXheaderMode merge_xheader_mode;
  MergeQualMode merge_qual_mode;
  MergeFilterMode merge_filter_mode;
  MergeInfoCmMode merge_info_mode;
  MergeInfoCmMode merge_cm_mode;
  SortMode merge_pheno_sort;
  SortMode merge_info_sort;
  uint32_t max_allele_ct;
  char* pgen_fname;
  char* pvar_fname;
  char* psam_fname;
  char* list_fname;
  char* list_base_dir;
} PmergeInfo;

typedef struct PgenDiffStruct {
  NONCOPYABLE(PgenDiffStruct);
  PgenDiffFlags flags;
  Dosage dosage_hap_tol; // missing value when 'dosage' not specified
  char* pgen_fname;
  char* pvar_fname;
  char* psam_fname;
} PgenDiffInfo;

void InitPmerge(PmergeInfo* pmerge_info_ptr);

void CleanupPmerge(PmergeInfo* pmerge_info_ptr);

void InitPgenDiff(PgenDiffInfo* pgen_diff_info_ptr);

void CleanupPgenDiff(PgenDiffInfo* pgen_diff_info_ptr);

PglErr Pmerge(const PmergeInfo* pmip, const char* sample_sort_fname, const char* missing_catname, const char* varid_template_str, const char* varid_multi_template_str, const char* varid_multi_nonsnp_template_str, const char* missing_varid_match, MiscFlags misc_flags, LoadFilterLogFlags load_filter_log_merge_flags, SortMode sample_sort_mode, FamCol fam_cols, int32_t missing_pheno, uint32_t new_variant_id_max_allele_slen, char input_missing_geno_char, uint32_t max_thread_ct, SortMode sort_vars_mode, char* pgenname, char* psamname, char* pvarname, char* outname, char* outname_end, ChrInfo* cip);

PglErr PgenDiff(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const PgenDiffInfo* pdip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t max_allele_ct1, uint32_t max_allele_slen, char input_missing_geno_char, uint32_t max_thread_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_MERGE_H__
