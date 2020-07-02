#ifndef __PLINK2_DATA_H__
#define __PLINK2_DATA_H__

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
  kfMakePlink2MJoin = (1 << 9),
  kfMakePlink2MJoinBoth = kfMakePlink2MJoin,
  kfMakePlink2MJoinSnps = kfMakePlink2MJoin + kfMakePlink2MSplitBase,
  kfMakePlink2MJoinAny = kfMakePlink2MJoin + 2 * kfMakePlink2MSplitBase,
  // don't support e.g. '+indels' for now due to lack of standardization re:
  // handling of MNP/'other' classes
  kfMakePlink2TrimAlts = (1 << 10),
  kfMakePlink2MMask = kfMakePlink2TrimAlts - kfMakePlink2MSplitBase,
  kfMakePlink2EraseAlt2Plus = (1 << 11),
  kfMakePlink2VaridSemicolon = (1 << 12),
  kfMakePlink2VaridDup = (1 << 13),

  kfMakePlink2SetHhMissing = (1 << 14),
  kfMakePlink2SetHhMissingKeepDosage = (1 << 15),
  kfMakePlink2SetMixedMtMissing = (1 << 16),
  kfMakePlink2SetMixedMtMissingKeepDosage = (1 << 17),
  kfMakePgenFormatBase = (1 << 18), // two bits
  kfMakePgenErasePhase = (1 << 20),
  kfMakePgenEraseDosage = (1 << 21),
  kfMakePgenFillMissingFromDosage = (1 << 22)
FLAGSET_DEF_END(MakePlink2Flags);

CONSTI32(kMaxInfoKeySlen, kMaxIdSlen);
#define MAX_INFO_KEY_SLEN_STR MAX_ID_SLEN_STR

PglErr WriteMapOrBim(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const double* variant_cms, uint32_t variant_ct, uint32_t max_allele_slen, char delim, uint32_t output_zst, uint32_t thread_ct);

PglErr PvarInfoOpenAndReloadHeader(const char* pvar_info_reload, uint32_t max_thread_ct, TextStream* pvar_reload_txsp, char** line_iterp, uint32_t* info_col_idx_ptr);

PglErr PvarInfoReload(uint32_t info_col_idx, uint32_t variant_uidx, TextStream* pvar_reload_txsp, char** line_iterp, uint32_t* trs_variant_uidx_ptr);

PglErr PvarInfoReloadAndWrite(uint32_t info_pr_flag_present, uint32_t info_col_idx, uint32_t variant_uidx, uint32_t is_pr, TextStream* pvar_reload_txsp, char** line_iterp, char** write_iter_ptr, uint32_t* rls_variant_uidx_ptr);

void AppendChrsetLine(const ChrInfo* cip, char** write_iter_ptr);

uint32_t ChrLenLbound(const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint32_t* new_variant_idx_to_old, uint32_t chr_fo_idx, uint32_t max_allele_slen, UnsortedVar vpos_sortstatus);

// fileformat, fileDate, source
// assumes adequate buffer space
void AppendVcfHeaderStart(uint32_t v43, char** cswritepp);

PglErr WriteFam(const char* outname, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uint32_t pheno_ct, char delim);

uint32_t DataFidColIsRequired(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t maybe_modifier);

uint32_t DataSidColIsRequired(const uintptr_t* sample_include, const char* sids, uint32_t sample_ct, uint32_t max_sid_blen, uint32_t maybe_modifier);

uint32_t DataParentalColsAreRequired(const uintptr_t* sample_include, const PedigreeIdInfo* piip, uint32_t sample_ct, uint32_t maybe_modifier);

char* AppendPhenoStr(const PhenoCol* pheno_col, const char* output_missing_pheno, uint32_t omp_slen, uint32_t sample_uidx, char* write_iter);

PglErr LoadAlleleAndGenoCounts(const uintptr_t* sample_include, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t first_hap_uidx, uint32_t is_minimac3_r2, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uintptr_t* allele_presents, uint64_t* allele_ddosages, uint64_t* founder_allele_ddosages, uint32_t* variant_missing_hc_cts, uint32_t* variant_missing_dosage_cts, uint32_t* variant_hethap_cts, STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), double* imp_r2_vals);

void ApplyHardCallThresh(const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t dosage_ct, uint32_t hard_call_halfdist, uintptr_t* genovec);

uint32_t ApplyHardCallThreshPhased(const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t dosage_ct, uint32_t hard_call_halfdist, uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo, uintptr_t* dphase_present, SDosage* dphase_delta, SDosage* tmp_dphase_delta);

PglErr MakePlink2NoVsort(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const char* varid_template_str, __maybe_unused const char* varid_multi_template_str, __maybe_unused const char* varid_multi_nonsnp_template_str, const char* missing_varid_match, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, UnsortedVar vpos_sortstatus, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, uint32_t new_variant_id_max_allele_slen, MiscFlags misc_flags, MakePlink2Flags make_plink2_flags, PvarPsamFlags pvar_psam_flags, uintptr_t pgr_alloc_cacheline_ct, char* xheader, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr MakePlink2Vsort(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const ChrIdx* chr_idxs, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, MakePlink2Flags make_plink2_flags, uint32_t use_nsort, PvarPsamFlags pvar_psam_flags, char* xheader, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr SampleSortFileMap(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* sample_sort_fname, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t** new_sample_idx_to_old_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_DATA_H__
