#ifndef __PLINK2_PERM_H__
#define __PLINK2_PERM_H__

// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
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
#include "include/plink2_base.h"
#include "plink2_common.h"
#include "plink2_compress_stream.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfPerm0,
  kfPermZs = (1 << 0),
  kfPermWithin = (1 << 1),
  kfPermMpermSaveBest = (1 << 2),
  kfPermMpermSaveAll = (1 << 3),
  kfPermCounts = (1 << 4),

  kfPermColChrom = (1 << 5),
  kfPermColPos = (1 << 6),
  kfPermColRef = (1 << 7),
  kfPermColAlt1 = (1 << 8),
  kfPermColAlt = (1 << 9),
  kfPermColMaybeprovref = (1 << 10),
  kfPermColProvref = (1 << 11),
  kfPermColOmitted = (1 << 12),
  kfPermColDefault = (kfPermColChrom | kfPermColRef | kfPermColAlt | kfPermColMaybeprovref | kfPermColOmitted),
  kfPermColAll = ((kfPermColOmitted * 2) - kfPermColChrom)
FLAGSET_DEF_END(PermFlags);

typedef struct PermConfigStruct {
  char* within_phenoname;
  uint32_t aperm_min;
  uint32_t aperm_max;
  double aperm_alpha;
  double aperm_beta;
  double aperm_init_interval;
  double aperm_interval_slope;
  PermFlags flags;
} PermConfig;

void InitPermConfig(PermConfig* perm_config_ptr);

void CleanupPermConfig(PermConfig* perm_config_ptr);

// (2^31 - 1000001) / 2
// ...reason for specifically -1000001 is no longer present, but this does need
// to be <2^30 so I'll leave this unchanged from plink 1.9
CONSTI32(kApermMax, 1073241823);

// permute_within_phenocol == nullptr ok
// orig_pheno_cc and cat_case_cts_ptr must either both be null or both be
// non-null
// on exit:
//   empty categories collapsed out
//   cat_cumulative_sizes[k] = cat_sizes[0] + ... + cat_sizes[k-1]
//   perm_idx_to_orig_uidx[cat_cumulative_sizes[k]..
//     (cat_cumulative_sizes[k+1]-1)] has category-k sample_uidx values
// error = out of memory
BoolErr PermuteWithinInit(const uintptr_t* sample_include, const PhenoCol* permute_within_phenocol, const uintptr_t* orig_pheno_cc, uint32_t raw_sample_ctl, uint32_t sample_ct, uint32_t* cat_ct_ptr, uint32_t** perm_idx_to_orig_uidx_ptr, uint32_t** cat_cumulative_sizes_ptr, uint32_t** cat_case_cts_ptr);

void PermuteWithinB(const uint32_t* perm_idx_to_orig_uidx, const uint32_t* cat_cumulative_sizes, const uint32_t* cat_case_cts, uint32_t raw_sample_ctl, uint32_t cat_ct, sfmt_t* sfmtp, uintptr_t* permuted_pheno_cc);

void PermuteWithinD(const uint32_t* perm_idx_to_orig_uidx, const uint32_t* cat_cumulative_sizes, const double* orig_pheno_qt, uint32_t cat_ct, sfmt_t* sfmtp, double* permuted_pheno_qt, uint32_t* perm_buf);

PglErr InitPermReportWriter(PermFlags perm_flags, uint32_t perm_adapt, uint32_t provref_col, uintptr_t overflow_buf_size, char* outname, char** outname_end2_ptr, CompressStreamState* css_ptr, char** cswritep_ptr);

uint32_t UpdateAdaptiveRemainingVariants(const uintptr_t* valid_alleles, const uintptr_t* valid_alleles_cumulative_popcounts_w, const uintptr_t* allele_idx_offsets, const AlleleCode* omitted_alleles, uint32_t raw_variant_ctl, uintptr_t valid_allele_ct, uint32_t perm_idx_start, uint32_t subbatch_size, double adaptive_intercept, double adaptive_slope, uint32_t* valid_allele_emp1_denoms, uintptr_t* remaining_variants, uint32_t* first_adapt_check_pidx_ptr, uint32_t* perms_actual_ptr);

BoolErr WritePermReportBody(const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* valid_alleles, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* omitted_alleles, const uintptr_t* nonref_flags, const double* orig_permstats, const uint32_t* valid_allele_emp1_ctx2s, const uint32_t* valid_allele_emp1_denoms, uint32_t orig_variant_ct, uint32_t all_nonref, uint32_t provref_col, double ln_pfilter, PermFlags perm_flags, uint32_t perms_total, uint32_t lower_stat_is_more_extreme, double* mperm_best_stats, CompressStreamState* css_ptr, char** cswritep_ptr, char* chr_buf);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_PERM_H__
