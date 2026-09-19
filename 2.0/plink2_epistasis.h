#ifndef __PLINK2_EPISTASIS_H__
#define __PLINK2_EPISTASIS_H__

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

#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// --epistasis-boost
FLAGSET_DEF_START()
  kfEpi0,
  kfEpiRegress = (1 << 0),
  kfEpiZs = (1 << 1),
  kfEpiRefBased = (1 << 2),
  kfEpiNoFirth = (1 << 3),
  kfEpiLog10 = (1 << 4),
  kfEpiColChrom = (1 << 5),
  kfEpiColPos = (1 << 6),
  kfEpiColMaybeA1 = (1 << 7),
  kfEpiColA1 = (1 << 8),
  kfEpiColBeta = (1 << 9),
  kfEpiColOrbeta = (1 << 10),
  kfEpiColSe = (1 << 11),
  kfEpiColStat = (1 << 12),
  kfEpiColDf = (1 << 13),
  kfEpiColP = (1 << 14),
  kfEpiColNsig = (1 << 15),
  kfEpiColNtot = (1 << 16),
  kfEpiColProp = (1 << 17),
  kfEpiColDefault = (kfEpiColChrom | kfEpiColMaybeA1 | kfEpiColOrbeta | kfEpiColSe | kfEpiColStat | kfEpiColDf | kfEpiColP | kfEpiColNsig | kfEpiColNtot | kfEpiColProp)
FLAGSET_DEF_END(EpiFlags);

typedef struct EpiInfoStruct {
  EpiFlags flags;
  double ln_epi1;
  double ln_epi2;
} EpiInfo;

void InitEpi(EpiInfo* epi_ip);

PglErr CalcEpiBoost(const uintptr_t* orig_sample_include, const PhenoCol* pheno_cols, const PhenoCol* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const EpiInfo* epi_ip, uint32_t raw_sample_ct, uint32_t pheno_ct, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_allele_slen, double output_min_ln, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr CalcEpiLinear(const uintptr_t* orig_sample_include, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const EpiInfo* epi_ip, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_allele_slen, double vif_thresh, double max_corr, double output_min_ln, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_EPISTASIS_H__
