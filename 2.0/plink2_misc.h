#ifndef __PLINK2_MISC_H__
#define __PLINK2_MISC_H__

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
  kfPhenoTransform0,
  kfPhenoTransformSplitCat = (1 << 0),
  kfPhenoTransformSplitCatOmitLast = (1 << 1),
  kfPhenoTransformSplitCatCovar01 = (1 << 2),
  kfPhenoTransformVstdCovar = (1 << 3),
  kfPhenoTransformVstdAll = (1 << 4),
  kfPhenoTransformQuantnormPheno = (1 << 5),
  kfPhenoTransformQuantnormCovar = (1 << 6),
  kfPhenoTransformQuantnormAll = (1 << 7),
FLAGSET_DEF_END(pheno_transform_flags_t);

FLAGSET_DEF_START()
  kfWriteCovar0,
  kfWriteCovarColMaybesid = (1 << 0),
  kfWriteCovarColSid = (1 << 1),
  kfWriteCovarColMaybeparents = (1 << 2),
  kfWriteCovarColParents = (1 << 3),
  kfWriteCovarColSex = (1 << 4),
  kfWriteCovarColPheno1 = (1 << 5),
  kfWriteCovarColPhenos = (1 << 6),
  kfWriteCovarColDefault = kfWriteCovarColMaybesid,
  kfWriteCovarColAll = ((kfWriteCovarColPhenos * 2) - kfWriteCovarColMaybesid)
FLAGSET_DEF_END(write_covar_flags_t);

FLAGSET_DEF_START()
  kfAlleleFreq0,
  kfAlleleFreqZs = (1 << 0),
  kfAlleleFreqCounts = (1 << 1),
  kfAlleleFreqBinsRefFname = (1 << 2),
  kfAlleleFreqBinsAlt1Fname = (1 << 3),
  kfAlleleFreqBinsOnly = (1 << 4),

  kfAlleleFreqColChrom = (1 << 5),
  kfAlleleFreqColPos = (1 << 6),
  kfAlleleFreqColRef = (1 << 7),
  kfAlleleFreqColAlt1 = (1 << 8),
  kfAlleleFreqColAlt = (1 << 9),
  kfAlleleFreqColReffreq = (1 << 10),
  kfAlleleFreqColAlt1freq = (1 << 11),
  kfAlleleFreqColAltfreq = (1 << 12),
  kfAlleleFreqColFreq = (1 << 13),
  kfAlleleFreqColEq = (1 << 14),
  kfAlleleFreqColEqz = (1 << 15),
  kfAlleleFreqColAlteq = (1 << 16),
  kfAlleleFreqColAlteqz = (1 << 17),
  kfAlleleFreqColNumeq = (1 << 18),
  kfAlleleFreqColAltnumeq = (1 << 19),
  kfAlleleFreqColMachR2 = (1 << 20),
  kfAlleleFreqColNobs = (1 << 21),
  kfAlleleFreqColDefault = (kfAlleleFreqColChrom | kfAlleleFreqColRef | kfAlleleFreqColAlt | kfAlleleFreqColAltfreq | kfAlleleFreqColNobs),
  kfAlleleFreqColAll = ((kfAlleleFreqColNobs * 2) - kfAlleleFreqColChrom),
  // only mutual exclusion is altfreq/freq/eq/eqz/alteq/alteqz/numeq/numeqz
  // don't force alt1freq/altfreq mutual exclusion since the former plays a bit
  // better with shell scripts
  // alt+alteqz is a bit silly, but I won't bother prohibiting it
  kfAlleleFreqColMutex = ((kfAlleleFreqColAltnumeq * 2) - kfAlleleFreqColAltfreq)
FLAGSET_DEF_END(allele_freq_t);

FLAGSET_DEF_START()
  kfMissingRpt0,
  kfMissingRptZs = (1 << 0),
  kfMissingRptSampleOnly = (1 << 1),
  kfMissingRptVariantOnly = (1 << 2),
  
  kfMissingRptScolMaybesid = (1 << 3),
  kfMissingRptScolSid = (1 << 4),
  kfMissingRptScolMisspheno1 = (1 << 5),
  kfMissingRptScolMissphenos = (1 << 6),
  kfMissingRptScolNmissDosage = (1 << 7),
  kfMissingRptScolNmiss = (1 << 8),
  kfMissingRptScolNmissHh = (1 << 9),
  kfMissingRptScolHethap = (1 << 10),
  kfMissingRptScolNobs = (1 << 11),
  kfMissingRptScolFmissDosage = (1 << 12),
  kfMissingRptScolFmiss = (1 << 13),
  kfMissingRptScolFmissHh = (1 << 14),
  kfMissingRptScolDefault = (kfMissingRptScolMaybesid | kfMissingRptScolMissphenos | kfMissingRptScolNmiss | kfMissingRptScolNobs | kfMissingRptScolFmiss),
  kfMissingRptScolAll = ((kfMissingRptScolFmissHh * 2) - kfMissingRptScolMaybesid),

  kfMissingRptVcolChrom = (1 << 15),
  kfMissingRptVcolPos = (1 << 16),
  kfMissingRptVcolRef = (1 << 17),
  kfMissingRptVcolAlt1 = (1 << 18),
  kfMissingRptVcolAlt = (1 << 19),
  kfMissingRptVcolNmissDosage = (1 << 20),
  kfMissingRptVcolNmiss = (1 << 21),
  kfMissingRptVcolNmissHh = (1 << 22),
  kfMissingRptVcolHethap = (1 << 23),
  kfMissingRptVcolNobs = (1 << 24),
  kfMissingRptVcolFmissDosage = (1 << 25),
  kfMissingRptVcolFmiss = (1 << 26),
  kfMissingRptVcolFmissHh = (1 << 27),
  kfMissingRptVcolFhethap = (1 << 28),
  kfMissingRptVcolDefault = (kfMissingRptVcolChrom | kfMissingRptVcolNmiss | kfMissingRptVcolNobs | kfMissingRptVcolFmiss),
  kfMissingRptVcolAll = ((kfMissingRptVcolFhethap * 2) - kfMissingRptVcolChrom)
FLAGSET_DEF_END(missing_rpt_t);

FLAGSET_DEF_START()
  kfGenoCounts0,
  kfGenoCountsZs = (1 << 0),
  
  kfGenoCountsColChrom = (1 << 1),
  kfGenoCountsColPos = (1 << 2),
  kfGenoCountsColRef = (1 << 3),
  kfGenoCountsColAlt1 = (1 << 4),
  kfGenoCountsColAlt = (1 << 5),
  kfGenoCountsColHomref = (1 << 6),
  kfGenoCountsColRefalt1 = (1 << 7),
  kfGenoCountsColRefalt = (1 << 8),
  kfGenoCountsColHomalt1 = (1 << 9),
  kfGenoCountsColAltxy = (1 << 10),
  kfGenoCountsColXy = (1 << 11),
  kfGenoCountsColHapref = (1 << 12),
  kfGenoCountsColHapalt1 = (1 << 13),
  kfGenoCountsColHapalt = (1 << 14),
  kfGenoCountsColHap = (1 << 15),
  kfGenoCountsColNumeq = (1 << 16),
  kfGenoCountsColMissing = (1 << 17),
  kfGenoCountsColNobs = (1 << 18),
  kfGenoCountsColDefault = (kfGenoCountsColChrom | kfGenoCountsColRef | kfGenoCountsColAlt | kfGenoCountsColHomref | kfGenoCountsColRefalt | kfGenoCountsColAltxy | kfGenoCountsColHapref | kfGenoCountsColHapalt | kfGenoCountsColMissing),
  kfGenoCountsColAll = ((kfGenoCountsColNobs * 2) - kfGenoCountsColChrom),
  
  kfGenoCountsColPairex = (kfGenoCountsColHapalt | kfGenoCountsColHap),
  kfGenoCountsColMutex = (kfGenoCountsColAltxy | kfGenoCountsColXy | kfGenoCountsColNumeq)
FLAGSET_DEF_END(geno_counts_t);

FLAGSET_DEF_START()
  kfHardy0,
  kfHardyZs = (1 << 0),
  kfHardyMidp = (1 << 1),
  
  kfHardyColChrom = (1 << 2),
  kfHardyColPos = (1 << 3),
  kfHardyColRef = (1 << 4),
  kfHardyColAlt1 = (1 << 5),
  kfHardyColAlt = (1 << 6),
  kfHardyColGcounts = (1 << 7),
  kfHardyColGcount1col = (1 << 8),
  kfHardyColHetfreq = (1 << 9),
  kfHardyColSexaf = (1 << 10),
  kfHardyColFemalep = (1 << 11),
  kfHardyColP = (1 << 12),
  kfHardyColDefault = (kfHardyColChrom | kfHardyColRef | kfHardyColAlt | kfHardyColGcounts | kfHardyColHetfreq | kfHardyColSexaf | kfHardyColP),
  kfHardyColAll = ((kfHardyColP * 2) - kfHardyColChrom)
FLAGSET_DEF_END(hardy_flags_t);

pglerr_t plink1_cluster_import(const char* within_fname, const char* catpheno_name, const char* family_missing_catname, const uintptr_t* sample_include, const char* sample_ids, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t mwithin_val, pheno_col_t** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr);

pglerr_t update_sample_sexes(const char* update_sex_fname, const uintptr_t* sample_include, char* sample_ids, uint32_t raw_sample_ct, uintptr_t sample_ct, uintptr_t max_sample_id_blen, uint32_t update_sex_colm2, uintptr_t* sex_nm, uintptr_t* sex_male);

pglerr_t split_cat_pheno(const char* split_cat_phenonames_flattened, const uintptr_t* sample_include, uint32_t raw_sample_ct, pheno_transform_flags_t pheno_transform_flags, pheno_col_t** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr, pheno_col_t** covar_cols_ptr, char** covar_names_ptr, uint32_t* covar_ct_ptr, uintptr_t* max_covar_name_blen_ptr);

pglerr_t pheno_variance_standardize(const char* vstd_flattened, const uintptr_t* sample_include, const char* pheno_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uint32_t is_covar_flag, pheno_col_t* pheno_cols);

pglerr_t pheno_quantile_normalize(const char* quantnorm_flattened, const uintptr_t* sample_include, const char* pheno_names, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uint32_t is_subset_flag, pheno_col_t* pheno_cols);

pglerr_t write_allele_freqs(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* founder_allele_dosages, const double* mach_r2_vals, const char* ref_binstr, const char* alt1_binstr, uint32_t variant_ct, uint32_t max_alt_allele_ct, uint32_t max_allele_slen, allele_freq_t allele_freq_modifier, uint32_t nonfounders, char* outname, char* outname_end);

pglerr_t write_geno_counts(const uintptr_t* sample_include, const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint32_t* raw_geno_cts, const uint32_t* x_male_geno_cts, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t variant_ct, uint32_t x_start, uint32_t max_allele_slen, geno_counts_t geno_counts_modifier, pgen_reader_t* simple_pgrp, char* outname, char* outname_end);

pglerr_t write_missingness_reports(const uintptr_t* sample_include, const uintptr_t* sex_male, const char* sample_ids, const char* sids, const pheno_col_t* pheno_cols, const char* pheno_names, const uint32_t* sample_missing_hc_cts, const uint32_t* sample_missing_dosage_cts, const uint32_t* sample_hethap_cts, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint32_t* variant_missing_cts, const uint32_t* variant_missing_dosage_cts, const uint32_t* variant_hethap_cts, uint32_t sample_ct, uint32_t male_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t variant_ct, uintptr_t max_allele_slen, uint32_t first_hap_uidx, missing_rpt_t missing_rpt_modifier, char* outname, char* outname_end);

pglerr_t compute_hwe_x_pvals(const uintptr_t* variant_include, const uint32_t* founder_raw_geno_cts, const uint32_t* founder_x_male_geno_cts, const uint32_t* founder_x_nosex_geno_cts, uint32_t x_start, uint32_t hwe_x_ct, uint32_t hwe_midp, uint32_t calc_thread_ct, double** hwe_x_pvals_ptr);

pglerr_t hardy_report(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint32_t* founder_raw_geno_cts, const uint32_t* founder_x_male_geno_cts, const uint32_t* founder_x_nosex_geno_cts, const double* hwe_x_pvals, uint32_t variant_ct, uint32_t hwe_x_ct, uint32_t max_allele_slen, double output_min_p, hardy_flags_t hardy_modifier, uint32_t nonfounders, char* outname, char* outname_end);

pglerr_t write_snplist(const uintptr_t* variant_include, char** variant_ids, uint32_t variant_ct, uint32_t output_zst, char* outname, char* outname_end);

pglerr_t write_covar(const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const pheno_col_t* covar_cols, const char* covar_names, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, write_covar_flags_t write_covar_flags, char* outname, char* outname_end);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_MISC_H__
