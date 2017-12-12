#ifndef __PLINK2_DATA_H__
#define __PLINK2_DATA_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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
  kfMakePlink2SetMixedMtMissing = (1 << 12),
  kfMakePgenFormatBase = (1 << 13), // two bits
  kfMakePgenEraseAlt2Plus = (1 << 15),
  kfMakePgenErasePhase = (1 << 16),
  kfMakePgenEraseDosage = (1 << 17)
FLAGSET_DEF_END(make_plink2_t);

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
  kfPsamColMaybesid = (1 << 10),
  kfPsamColSid = (1 << 11),
  kfPsamColMaybeparents = (1 << 12),
  kfPsamColParents = (1 << 13),
  kfPsamColSex = (1 << 14),
  kfPsamColPheno1 = (1 << 15),
  kfPsamColPhenos = (1 << 16),
  kfPsamColDefault = (kfPsamColMaybesid | kfPsamColMaybeparents | kfPsamColSex | kfPsamColPhenos),
  kfPsamColAll = ((kfPsamColPhenos * 2) - kfPsamColMaybesid)
FLAGSET_DEF_END(pvar_psam_t);

pglerr_t write_map_or_bim(const char* outname, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const double* variant_cms, uint32_t variant_ct, uint32_t max_allele_slen, char delim, uint32_t output_zst);

pglerr_t pvar_info_open_and_reload_header(const char* pvar_info_reload, gzFile* gz_pvar_reload_ptr, char** loadbuf_ptr, uintptr_t* loadbuf_size_ptr, uint32_t* info_col_idx_ptr);

pglerr_t pvar_info_reload_and_write(uintptr_t loadbuf_size, uint32_t xheader_info_pr, uint32_t info_col_idx, uint32_t variant_uidx, uint32_t is_pr, gzFile gz_pvar_reload, char** write_iter_ptr, uint32_t* gz_variant_uidx_ptr, char* loadbuf);

void append_chrset_line(const chr_info_t* cip, char** write_iter_ptr);

pglerr_t write_fam(const char* outname, const uintptr_t* sample_include, const char* sample_ids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, char delim);

uint32_t is_parental_info_present(const uintptr_t* sample_include, const char* paternal_ids, const char* maternal_ids, uint32_t sample_ct, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen);

char* append_pheno_str(const pheno_col_t* pheno_col, const char* output_missing_pheno, uint32_t omp_slen, uint32_t sample_uidx, char* write_iter);

// dosage_int = 0..2 value in 16384ths
// returns distance from 0.5 or 1.5 in 16384ths, whichever is closer
HEADER_CINLINE2 uint32_t biallelic_dosage_halfdist(uint32_t dosage_int) {
  const uint32_t dosage_int_rem = dosage_int & (kDosageMid - 1);
  return abs_int32(((int32_t)dosage_int_rem) - kDosage4th);
}

pglerr_t load_allele_and_geno_counts(const uintptr_t* sample_include, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uintptr_t* variant_allele_idxs, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t first_hap_uidx, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, uint64_t* allele_dosages, uint64_t* founder_allele_dosages, uint32_t* variant_missing_hc_cts, uint32_t* variant_missing_dosage_cts, uint32_t* variant_hethap_cts, uint32_t* raw_geno_cts, uint32_t* founder_raw_geno_cts, uint32_t* x_male_geno_cts, uint32_t* founder_x_male_geno_cts, uint32_t* x_nosex_geno_cts, uint32_t* founder_x_nosex_geno_cts, double* mach_r2_vals);

pglerr_t make_plink2_no_vsort(const char* xheader, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, make_plink2_t make_plink2_modifier, pvar_psam_t pvar_psam_modifier, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end);

pglerr_t make_plink2_vsort(const char* xheader, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const chr_idx_t* chr_idxs, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, make_plink2_t make_plink2_modifier, uint32_t use_nsort, pvar_psam_t pvar_psam_modifier, pgen_reader_t* simple_pgrp, char* outname, char* outname_end);

pglerr_t sample_sort_file_map(const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* sample_sort_fname, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t sid_col_present, uint32_t** new_sample_idx_to_old_ptr);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_DATA_H__
