#ifndef __PLINK2_FAMILY_H__
#define __PLINK2_FAMILY_H__

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

FLAGSET_DEF_START()
  kfMendel0,
  kfMendelDuos = (1 << 0),
  kfMendelMissingInDenom = (1 << 1),
  kfMendelFilterVarFirst = (1 << 2),
  kfMendelRptZs = (1 << 3),
  kfMendelRptSummariesOnly = (1 << 4),
  kfMendelRptColMaybefid = (1 << 5),
  kfMendelRptColFid = (1 << 6),
  kfMendelRptColMaybesid = (1 << 7),
  kfMendelRptColSid = (1 << 8),
  kfMendelRptColChrom = (1 << 9),
  kfMendelRptColPos = (1 << 10),
  kfMendelRptColRef = (1 << 11),
  kfMendelRptColAlt = (1 << 12),
  kfMendelRptColCode = (1 << 13),
  kfMendelRptColError = (1 << 14),
  kfMendelRptIcolTrionum = (1 << 15),
  kfMendelRptIcolN = (1 << 16),
  kfMendelRptIcolNobs = (1 << 17),
  kfMendelRptIcolFrac = (1 << 18),
  kfMendelRptLcolN = (1 << 19),
  kfMendelRptLcolNobs = (1 << 20),
  kfMendelRptLcolFrac = (1 << 21),
  kfMendelRptColDefault = (kfMendelRptColMaybefid | kfMendelRptColMaybesid | kfMendelRptColChrom | kfMendelRptColCode | kfMendelRptColError | kfMendelRptIcolN | kfMendelRptLcolN),
  kfMendelRptColAll = ((kfMendelRptLcolFrac * 2) - kfMendelRptColMaybefid)
FLAGSET_DEF_END(MendelFlags);

typedef struct MendelStruct {
  MendelFlags flags;
  double max_trio_error;
  double max_var_error;
  double exclude_one_ratio;  // 0 = off, -1 = always child
} MendelInfo;

void InitMendel(MendelInfo* mendel_info_ptr);

// Main difference from PLINK 1.9 is that subsetted instead of raw
// sample-indexes are used (don't really have a choice for --set-me-missing,
// and PgenReader takes care of subsetting for other callers).
typedef struct FamilyStruct {
  // Paternal sample_idxs in low 32 bits, maternal in high 32, sorted by
  // lowest-child-sample_idx.
  // In include_duos case, missing parents are encoded as uidx=sample_ct.
  const uint64_t* family_list;

  // Child sample_idx in low 32 bits, family_idx (into family_list) in high 32,
  // sorted
  const uint64_t* trio_list;

  // [3k]: child ID
  // [3k+1]: paternal ID
  // [3k+2]: maternal ID
  const uint32_t* trio_lookup;

  // Uses trio_list indexes.
  // (may not need family_fids)
  const char* trio_fids;

  // Uses sample_idxs.
  const char* iids;
  const char* sids;

  // Sparse optimization support.
  // Elements [parent_to_trio_offsets[sample_idx],
  // parent_to_trio_offsets[sample_idx+1]) of parent_to_trio_idxs[] are the
  // trio_idxs where sample_idx is a parent.
  // child_to_trio_idxs[sample_idx] is either UINT32_MAX if the sample is not a
  // child in a trio, or the trio_idx if it is.
  const uint32_t* parent_to_trio_idxs;
  const uint32_t* parent_to_trio_offsets;
  const uint32_t* child_to_trio_idxs;

  uint32_t family_ct;
  uint32_t trio_ct;

  uintptr_t max_fid_blen;
  uintptr_t max_iid_blen;
  uintptr_t max_sid_blen;
} FamilyInfo;

FLAGSET_DEF_START()
  kfTrio0,
  kfTrioDuos = (1 << 0),
  kfTrioPopulateIds = (1 << 1),
  kfTrioPopulateSids = (1 << 2)
FLAGSET_DEF_END(TrioFlags);

void PreinitFamilyInfo(FamilyInfo* fip);

PglErr GetTriosAndFamilies(const uintptr_t* orig_sample_include, const PedigreeIdInfo* piip, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, uint32_t raw_sample_ct, TrioFlags flags, uint32_t* sample_ct_ptr, uintptr_t* trio_sample_include, FamilyInfo* fip);

PglErr MendelErrorScan(const PedigreeIdInfo* piip, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const MendelInfo* mip, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t generate_reports, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uintptr_t* sample_include, uintptr_t* variant_include, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_FAMILY_H__
