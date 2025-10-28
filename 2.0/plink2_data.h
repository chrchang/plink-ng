#ifndef __PLINK2_DATA_H__
#define __PLINK2_DATA_H__

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

#include "include/pgenlib_misc.h"
#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "include/plink2_text.h"
#include "plink2_cmdline.h"
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
  kfMakePlink2VaridSemicolon = (1 << 11),
  kfMakePlink2VaridDup = (1 << 12),

  kfMakePlink2SetInvalidHaploidMissing = (1 << 13),
  kfMakePlink2SetInvalidHaploidMissingKeepDosage = (1 << 14),
  kfMakePlink2SetMixedMtMissing = (1 << 15),
  kfMakePlink2SetMixedMtMissingKeepDosage = (1 << 16),
  kfMakePlink2SetMeMissing = (1 << 17),
  kfMakePlink2FillMissingWithRef = (1 << 18),
  kfMakePgenFormatBase = (1 << 19), // two bits
  kfMakePgenErasePhase = (1 << 21),
  kfMakePgenEraseDosage = (1 << 22),
  kfMakePgenFillMissingFromDosage = (1 << 23),
  kfMakePgenWriterVer = (1 << 24)
FLAGSET_DEF_END(MakePlink2Flags);

FLAGSET_DEF_START()
  kfImport0,
  kfImportKeepAutoconv = (1 << 0),
  kfImportKeepAutoconvVzs = (1 << 1),
  kfImportDoubleId = (1 << 2),
  kfImportVcfRefNMissing = (1 << 3),
  kfImportLaxChrX = (1 << 4),
  kfImportPolyploidMissing = (1 << 5),
  kfImportPolyploidExplicitError = (1 << 6),
  kfImportVcfAllowNoNonvar = (1 << 7),
  kfImportLaxBgen = (1 << 8),
  kfImportEigNohash = (1 << 9)
FLAGSET_DEF_END(ImportFlags);

CONSTI32(kMaxInfoKeySlen, kMaxIdSlen);
#define MAX_INFO_KEY_SLEN_STR MAX_ID_SLEN_STR

// We only need to distinguish between the following INFO-value-type cases:
// Number=0 (flag), Number=<positive integer>, Number=., Number=A, and
// Number=R.  (Number=G is being treated as Number=., since otherwise chrX is a
// nightmare; see https://github.com/samtools/hts-specs/issues/272 for some
// discussion.)  We use negative numbers to represent the last 3 cases in
// InfoVtype.
CONSTI32(kInfoVtypeUnknown, -1);
CONSTI32(kInfoVtypeA, -2);
CONSTI32(kInfoVtypeR, -3);
// Used in plink2_data (though not plink2_merge) for INFO/PR flag.
CONSTI32(kInfoVtypeSkip, -4);

HEADER_INLINE uint32_t IsInfoVtypeARSkip(int32_t info_vtype) {
  return (info_vtype <= kInfoVtypeA);
}

// Main fixed data structure when splitting/joining INFO is a hashmap of keys.
// Behavior when splitting:
// - Field order in the original variant is retained.
// - Number >= 0 and Number=. don't require any special handling, just copy the
//   entire key=value pair (or lone key, in the Flag case).
// - Number=A and Number=R require splitting the value on ',' and verifying the
//   comma count is correct, but is otherwise straightforward since alleles
//   can't be permuted.
// - Number=G is being treated as Number=.
// When joining:
// - Field order is determined by header line order.
// - Number=./G and Number>0 just require a buffer of size ~info_reload_slen,
//   and a boolean indicating whether no mismatch has been found.
// - Number=0 (Flag) requires a single boolean, we perform an or operation.
// - Number=A/R are the messy ones: we need to have enough space for
//   max_write_allele_ct (or that minus 1) comma-separated values in the =A and
//   =R cases.
typedef struct InfoVtypeStruct {
  NONCOPYABLE(InfoVtypeStruct);
  int32_t num;
  char key[];
} InfoVtype;

PglErr WriteMapOrBim(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, const double* variant_cms, uint32_t variant_ct, uint32_t max_allele_slen, char delim, char output_missing_geno_char, uint32_t output_zst, uint32_t thread_ct);

PglErr WritePsam(const char* outname, const uintptr_t* sample_include, const SampleIdInfo* siip, const ParentalIdInfo* parental_id_infop, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const char* output_missing_pheno, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, PvarPsamFlags pvar_psam_flags, uint32_t psam_01);

PglErr PvarInfoOpenAndReloadHeader(const char* pvar_info_reload, uint32_t max_thread_ct, TextStream* pvar_reload_txsp, char** line_iterp, uint32_t* info_col_idx_ptr);

PglErr PvarInfoReload(uint32_t info_col_idx, uint32_t variant_uidx, TextStream* pvar_reload_txsp, char** line_iterp, uint32_t* trs_variant_uidx_ptr);

PglErr PvarInfoReloadAndWrite(uint32_t info_pr_flag_present, uint32_t info_col_idx, uint32_t variant_uidx, uint32_t is_pr, TextStream* pvar_reload_txsp, char** line_iterp, char** write_iter_ptr, uint32_t* rls_variant_uidx_ptr);

void AppendChrsetLine(const ChrInfo* cip, char** write_iter_ptr);

// uint32_t ChrLenLbound(const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint32_t* new_variant_idx_to_old, uint32_t chr_fo_idx, uint32_t max_allele_slen, UnsortedVar vpos_sortstatus);

uint64_t FindContigHeaderLineLengthStr(const char* contig_name_end);

char* WriteRestOfHeaderLine(const char* hkvline_iter, const char* idval, const char* closing_gt, uint32_t id_slen, char* write_iter);

BoolErr FixAndWriteContigHeaderLine(const char* header_line_iter, const char* id_ptr, const char* closing_gt, uint32_t id_slen, uint32_t contig_len, char** writep_ptr);

// fileformat, fileDate, source
// assumes adequate buffer space
void AppendVcfHeaderStart(uint32_t v43, char** cswritepp);

PglErr WriteFam(const char* outname, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uint32_t* new_sample_idx_to_old, const char* legacy_output_missing_pheno, uint32_t sample_ct, uint32_t pheno_ct, char delim);

uint32_t DataFidColIsRequired(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t maybe_modifier);

uint32_t DataSidColIsRequired(const uintptr_t* sample_include, const char* sids, uint32_t sample_ct, uint32_t max_sid_blen, uint32_t maybe_modifier);

uint32_t DataParentalColsAreRequired(const uintptr_t* sample_include, const SampleIdInfo* siip, const ParentalIdInfo* parental_id_infop, uint32_t sample_ct, uint32_t maybe_modifier);

char* AppendPhenoStr(const PhenoCol* pheno_col, const char* output_missing_pheno, uint32_t omp_slen, uint32_t sample_uidx, char* write_iter);

char* AppendPhenoStrEx(const PhenoCol* pheno_col, const char* output_missing_pheno, uint32_t omp_slen, uint32_t sample_uidx, char ctrl_char, char* write_iter);

uint32_t CopyAndPermute8bit(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uint32_t* __restrict old_sample_idx_to_new, uint32_t sample_ct, uint32_t val_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals);

uint32_t CopyAndPermute16bit(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uint32_t* __restrict old_sample_idx_to_new, uint32_t sample_ct, uint32_t val_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals);

uint32_t Dense8bitToSparse(const uintptr_t* __restrict set, uint32_t sample_ctl, void* __restrict vals);

uint32_t Dense16bitToSparse(const uintptr_t* __restrict set, uint32_t sample_ctl, void* __restrict vals);

PglErr LoadAlleleAndGenoCounts(const uintptr_t* sample_include, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t first_hap_uidx, uint32_t is_minimac3_r2, uint32_t y_nosex_missing_stats, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uintptr_t* allele_presents, uint64_t* allele_ddosages, uint64_t* founder_allele_ddosages, uint32_t* variant_missing_hc_cts, uint32_t* variant_missing_dosage_cts, uint32_t* variant_hethap_cts, STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), double* imp_r2_vals);

void ApplyHardCallThresh(const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t dosage_ct, uint32_t hard_call_halfdist, uintptr_t* genovec);

uint32_t ApplyHardCallThreshPhased(const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t dosage_ct, uint32_t hard_call_halfdist, uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo, uintptr_t* dphase_present, SDosage* dphase_delta, SDosage* tmp_dphase_delta);

BoolErr FillInfoVtypeNum(const char* numstr, uint32_t num_slen, uint32_t* info_has_g_keyp, int32_t* info_vtype_num_ptr);

PglErr MakePlink2NoVsort(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const char* varid_template_str, __maybe_unused const char* varid_multi_template_str, __maybe_unused const char* varid_multi_nonsnp_template_str, const char* missing_varid_match, const char* output_missing_pheno, const char* legacy_output_missing_pheno, const uint32_t* contig_lens, const char* writer_ver, const char* zero_cluster_fname, const char* zero_cluster_phenoname, const FlipInfo* flip_info_ptr, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, char output_missing_geno_char, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, uint32_t new_variant_id_max_allele_slen, MiscFlags misc_flags, MakePlink2Flags make_plink2_flags, PvarPsamFlags pvar_psam_flags, uint32_t mendel_duos, uintptr_t pgr_alloc_cacheline_ct, char* xheader, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr MakePlink2Vsort(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const AlleleCode* allele_permute, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const char* output_missing_pheno, const char* legacy_output_missing_pheno, const uint32_t* contig_lens, const char* rename_chrs_fname, const TwoColParams* update_chr_flag, const char* writer_ver, const char* zero_cluster_fname, const char* zero_cluster_phenoname, const FlipInfo* flip_info_ptr, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, char output_missing_geno_char, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, MiscFlags misc_flags, MakePlink2Flags make_plink2_flags, uint32_t use_nsort, PvarPsamFlags pvar_psam_flags, uint32_t mendel_duos, ChrInfo* cip, char* xheader, ChrIdx* chr_idxs, PgenReader* simple_pgrp, char* outname, char* outname_end);

PglErr SampleSortFileMap(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* sample_sort_fname, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t** new_sample_idx_to_old_ptr);

PglErr LoadSampleMissingCts(const uintptr_t* sample_include, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const char* calc_descrip_str, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t y_nosex_missing_stats, uint32_t raw_sample_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uint32_t* sample_missing_hc_cts, uint32_t* sample_missing_dosage_cts, uint32_t* sample_hethap_cts);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_DATA_H__
