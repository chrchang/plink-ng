#ifndef __PLINK2_ADJUST_H__
#define __PLINK2_ADJUST_H__

// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang, Benjamin Demaille.
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

#include "include/pgenlib_misc.h"
#include "include/plink2_base.h"
#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfAdjust0,
  kfAdjustGc = (1 << 0),
  kfAdjustLog10 = (1 << 1),
  kfAdjustZs = (1 << 2),
  kfAdjustInputLog10 = (1 << 3),

  kfAdjustColChrom = (1 << 4),
  kfAdjustColPos = (1 << 5),
  kfAdjustColRef = (1 << 6),
  kfAdjustColAlt1 = (1 << 7),
  kfAdjustColAlt = (1 << 8),
  kfAdjustColMaybeprovref = (1 << 9),
  kfAdjustColProvref = (1 << 10),
  kfAdjustColA1 = (1 << 11),
  kfAdjustColUnadj = (1 << 12),
  kfAdjustColGc = (1 << 13),
  kfAdjustColQq = (1 << 14),
  kfAdjustColBonf = (1 << 15),
  kfAdjustColHolm = (1 << 16),
  kfAdjustColSidakss = (1 << 17),
  kfAdjustColSidaksd = (1 << 18),
  kfAdjustColFdrbh = (1 << 19),
  kfAdjustColFdrby = (1 << 20),
  kfAdjustColDefault = (kfAdjustColChrom | kfAdjustColMaybeprovref | kfAdjustColA1 | kfAdjustColUnadj | kfAdjustColGc | kfAdjustColBonf | kfAdjustColHolm | kfAdjustColSidakss | kfAdjustColSidaksd | kfAdjustColFdrbh | kfAdjustColFdrby),
  kfAdjustColAll = ((kfAdjustColFdrby * 2) - kfAdjustColChrom)
FLAGSET_DEF_END(AdjustFlags);

typedef struct AdjustInfoStruct {
  AdjustFlags flags;
  double lambda;
} AdjustInfo;

FLAGSET_DEF_START()
  kfMeta0,
  kfMetaLogscale = (1 << 0),
  kfMetaQt = (1 << 1),
  kfMetaNoMap = (1 << 2),
  kfMetaNoAllele = (1 << 3),
  kfMetaStudy = (1 << 4),
  kfMetaReportAll = (1 << 5),
  kfMetaWeightedZ = (1 << 6),
  kfMetaZs = (1 << 7)
FLAGSET_DEF_END(MetaFlags);

typedef struct MetaInfoStruct {
  NONCOPYABLE(MetaInfoStruct);
  char* fnames;  // multistr
  char* chr_field;
  char* snp_field;
  char* bp_field;
  char* a1_field;
  char* a2_field;
  char* p_field;
  char* se_field;
  char* ess_field;
  MetaFlags flags;
} MetaInfo;

void InitMeta(MetaInfo* mip);

void CleanupMeta(MetaInfo* mip);

PglErr MetaAnalysis(const MetaInfo* mip, uint32_t max_thread_ct, char* outname, char* outname_end);

typedef struct AdjustFileInfoStruct {
  NONCOPYABLE(AdjustFileInfoStruct);
  AdjustInfo base;
  char* fname;
  char* test_name;
  char* chr_field;
  char* pos_field;
  char* id_field;
  char* ref_field;
  char* alt_field;
  char* provref_field;
  char* a1_field;
  char* test_field;
  char* p_field;
  char* chisq_field;
} AdjustFileInfo;

void InitAdjust(AdjustInfo* adjust_info_ptr, AdjustFileInfo* adjust_file_info_ptr);

void CleanupAdjust(AdjustFileInfo* adjust_file_info_ptr);

PglErr Multcomp(const uintptr_t* variant_include, const ChrInfo* cip, const char* const* chr_ids, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_include, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* nonref_flags, const char* const* loaded_a1, const AdjustInfo* adjust_info_ptr, const double* ln_pvals, const double* chisqs, uint32_t raw_variant_ct, uintptr_t orig_allele_ct, uint32_t max_allele_slen, PgenGlobalFlags gflags, double ln_pfilter, double output_min_ln, uint32_t skip_gc, uint32_t max_thread_ct, char* outname, char* outname_end);

PglErr AdjustFile(const AdjustFileInfo* afip, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_ADJUST_H__
