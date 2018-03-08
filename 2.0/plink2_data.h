#ifndef __PLINK2_DATA_H__
#define __PLINK2_DATA_H__

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
  kfMake0,
  kfMakeBed = (1 << 0),
  kfMakeBim = (1 << 1),
  kfMakeFam = (1 << 2),
  kfMakePgen = (1 << 3),
  kfMakePvar = (1 << 4),
  kfMakePsam = (1 << 5),
  kfMakeBimZs = (1 << 6),
  kfMakePlink2MSplitBase = (1 << 7), // three bits for multiallelic mode
  kfMakePlink2MSplitAll = kfMakePlink2MSplitBase,
  kfMakePlink2MSplitSnps = 2 * kfMakePlink2MSplitBase,
  kfMakePlink2MMerge = (1 << 9),
  kfMakePlink2MMergeBoth = kfMakePlink2MMerge,
  kfMakePlink2MMergeSnps = kfMakePlink2MMerge + kfMakePlink2MSplitBase,
  kfMakePlink2MMergeAny = kfMakePlink2MMerge + 2 * kfMakePlink2MSplitBase,
  // don't support e.g. '+indels' for now due to lack of standardization re:
  // handling of MNP/'other' classes
  kfMakePlink2TrimAlts = (1 << 10),
  kfMakePlink2MMask = kfMakePlink2TrimAlts - kfMakePlink2MSplitBase,
  kfMakePlink2SetHhMissing = (1 << 11),
  kfMakePlink2SetHhMissingKeepDosage = (1 << 12),
  kfMakePlink2SetMixedMtMissing = (1 << 13),
  kfMakePlink2SetMixedMtMissingKeepDosage = (1 << 14),
  kfMakePgenFormatBase = (1 << 15), // two bits
  kfMakePgenEraseAlt2Plus = (1 << 16),
  kfMakePgenErasePhase = (1 << 17),
  kfMakePgenEraseDosage = (1 << 18)
FLAGSET_DEF_END(MakePlink2Flags);

FLAGSET_DEF_START()
  kfPvarPsam0,
  kfPvarZs = (1 << 0),

  kfPvarColXheader = (1 << 1),
  kfPvarColMaybequal = (1 << 2),
  kfPvarColQual = (1 << 3),
  kfPvarColMaybefilter = (1 << 4),
  kfPvarColFilter = (1 << 5),
  kfPvarColMaybeinfo = (1 << 6),
  kfPvarColInfo = (1 << 7),
  kfPvarColXinfo = (kfPvarColInfo * 2) - kfPvarColMaybeinfo,
  kfPvarColMaybecm = (1 << 8),
  kfPvarColCm = (1 << 9),
  kfPvarColDefault = (kfPvarColXheader | kfPvarColMaybequal | kfPvarColMaybefilter | kfPvarColMaybeinfo | kfPvarColMaybecm),
  kfPvarColAll = ((kfPvarColCm * 2) - kfPvarColXheader),
  kfPsamColMaybefid = (1 << 10),
  kfPsamColFid = (1 << 11),
  kfPsamColMaybesid = (1 << 12),
  kfPsamColSid = (1 << 13),
  kfPsamColMaybeparents = (1 << 14),
  kfPsamColParents = (1 << 15),
  kfPsamColSex = (1 << 16),
  kfPsamColPheno1 = (1 << 17),
  kfPsamColPhenos = (1 << 18),
  kfPsamColDefault = (kfPsamColMaybefid | kfPsamColMaybesid | kfPsamColMaybeparents | kfPsamColSex | kfPsamColPhenos),
  kfPsamColAll = ((kfPsamColPhenos * 2) - kfPsamColMaybefid)
FLAGSET_DEF_END(PvarPsamFlags);

PglErr WriteMapOrBim(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const double* variant_cms, uint32_t variant_ct, uint32_t max_allele_slen, char delim, uint32_t output_zst, uint32_t thread_ct);

PglErr PvarInfoOpenAndReloadHeader(const char* pvar_info_reload, ReadLineStream* pvar_reload_rlsp, char** line_iterp, uint32_t* info_col_idx_ptr);

PglErr PvarInfoReloadAndWrite(uint32_t info_pr_flag_present, uint32_t info_col_idx, uint32_t variant_uidx, uint32_t is_pr, ReadLineStream* pvar_reload_rls, char** line_iterp, char** write_iter_ptr, uint32_t* rls_variant_uidx_ptr);

void AppendChrsetLine(const ChrInfo* cip, char** write_iter_ptr);

PglErr WriteFam(const char* outname, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uint32_t pheno_ct, char delim);

uint32_t DataFidColIsRequired(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t maybe_modifier);

uint32_t DataSidColIsRequired(const uintptr_t* sample_include, const char* sids, uint32_t sample_ct, uint32_t max_sid_blen, uint32_t maybe_modifier);

uint32_t DataParentalColsAreRequired(const uintptr_t* sample_include, const PedigreeIdInfo* piip, uint32_t sample_ct, uint32_t maybe_modifier);

char* AppendPhenoStr(const PhenoCol* pheno_col, const char* output_missing_pheno, uint32_t omp_slen, uint32_t sample_uidx, char* write_iter);

// dosage_int = 0..2 value in 16384ths
// returns distance from 0.5 or 1.5 in 16384ths, whichever is closer
HEADER_INLINE uint32_t BiallelicDosageHalfdist(uint32_t dosage_int) {
  const uint32_t dosage_int_rem = dosage_int & (kDosageMid - 1);
  return abs_i32(S_CAST(int32_t, dosage_int_rem) - kDosage4th);
}

PglErr LoadAlleleAndGenoCounts(const uintptr_t* sample_include, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* variant_allele_idxs, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t first_hap_uidx, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uint64_t* allele_dosages, uint64_t* founder_allele_dosages, uint32_t* variant_missing_hc_cts, uint32_t* variant_missing_dosage_cts, uint32_t* variant_hethap_cts, uint32_t* raw_geno_cts, uint32_t* founder_raw_geno_cts, uint32_t* x_male_geno_cts, uint32_t* founder_x_male_geno_cts, uint32_t* x_nosex_geno_cts, uint32_t* founder_x_nosex_geno_cts, double* mach_r2_vals);

PglErr MakePlink2NoVsort(const char* xheader, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, MakePlink2Flags make_plink2_flags, PvarPsamFlags pvar_psam_flags, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr MakePlink2Vsort(const char* xheader, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const ChrIdx* chr_idxs, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, MakePlink2Flags make_plink2_flags, uint32_t use_nsort, PvarPsamFlags pvar_psam_flags, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr SampleSortFileMap(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* sample_sort_fname, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t** new_sample_idx_to_old_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_DATA_H__
