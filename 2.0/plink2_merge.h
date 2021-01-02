#ifndef __PLINK2_MERGE_H__
#define __PLINK2_MERGE_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2021 Shaun Purcell,
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
  kfPmergePhenoInnerJoin = (1 << 2)
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
  kMergeQualInfoModeErase,
  kMergeQualInfoModeNmMatch,
  kMergeQualInfoModeNmFirst,
  kMergeQualInfoModeFirst
ENUM_U31_DEF_END(MergeQualInfoMode);

ENUM_U31_DEF_START()
  kMergeFilterModeErase,
  kMergeFilterModeNmMatch,
  kMergeFilterModeNmFirst,
  kMergeFilterModeFirst,
  kMergeFilterModeNonpassUnion
ENUM_U31_DEF_END(MergeFilterMode);

// this is a hybrid, only kfSortFileSid is actually a flag
FLAGSET_DEF_START()
  kfSort0,
  kfSortNone = (1 << 0),
  kfSortNatural = (1 << 1),
  kfSortAscii = (1 << 2),
  kfSortFile = (1 << 3),
  kfSortFileSid = (1 << 4)
FLAGSET_DEF_END(SortFlags);

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
  kfPgenDiffColMaybefid = (1 << 7),
  kfPgenDiffColFid = (1 << 8),
  kfPgenDiffColMaybesid = (1 << 9),
  kfPgenDiffColSid = (1 << 10),
  kfPgenDiffColGeno = (1 << 11),
  kfPgenDiffColDefault = (kfPgenDiffColId | kfPgenDiffColMaybefid | kfPgenDiffColMaybesid | kfPgenDiffColGeno),
  kfPgenDiffColAll = ((kfPgenDiffColGeno * 2) - kfPgenDiffColPos)
FLAGSET_DEF_END(PgenDiffFlags);

typedef struct PmergeStruct {
  NONCOPYABLE(PmergeStruct);
  PmergeFlags flags;
  PmergeListMode list_mode;
  MergeMode merge_mode;
  MergePhenoMode merge_pheno_mode;
  MergeXheaderMode merge_xheader_mode;
  MergeQualInfoMode merge_qual_mode;
  MergeFilterMode merge_filter_mode;
  MergeQualInfoMode merge_info_mode;
  SortFlags merge_info_sort;
  uint32_t max_allele_ct;
  char* pgen_fname;
  char* pvar_fname;
  char* psam_fname;
  char* list_fname;
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

PglErr Pmerge(const PmergeInfo* pmip, const char* sample_sort_fname, MiscFlags misc_flags, SortFlags sample_sort_flags, FamCol fam_cols, uint32_t max_thread_ct, char* pgenname, char* psamname, char* pvarname, char* outname, char* outname_end, ChrInfo* cip);

PglErr PgenDiff(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const PgenDiffInfo* pdip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t max_allele_ct1, uint32_t max_allele_slen, uint32_t max_thread_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_MERGE_H__
