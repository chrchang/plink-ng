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


#include "plink2_compress_stream.h"
#include "plink2_export.h"
#include "plink2_filter.h"
#include "plink2_glm.h"
#include "plink2_import.h"
#include "plink2_ld.h"
#include "plink2_matrix_calc.h"
#include "plink2_misc.h"
#include "plink2_psam.h"
#include "plink2_pvar.h"
#include "plink2_random.h"
#include "plink2_set.h"

#include <time.h>  // time()
#include <unistd.h>  // unlink()

#include "plink2_help.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static const char ver_str[] = "PLINK v2.00a2"
#ifdef NOLAPACK
  "NL"
#endif
#ifdef __LP64__
#  ifdef LAPACK_ILP64
    "LM"
#  endif
#  ifdef USE_AVX2
    " AVX2"
#  else
#    ifdef USE_SSE42
      " SSE4.2"
#    else
      " 64-bit"
#    endif
#  endif
#else
  " 32-bit"
#endif

#ifdef USE_MKL
  " Intel"
#endif
  " (23 Mar 2018)";
static const char ver_str2[] =
  // include leading space if day < 10, so character length stays the same
  ""
#ifndef LAPACK_ILP64
  "  "
#endif
#ifndef USE_MKL
  "      "
#endif
#ifdef USE_AVX2
  "  "
#endif
#ifndef NOLAPACK
  "  "
#endif
  "   www.cog-genomics.org/plink/2.0/\n"
  "(C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3\n";
static const char errstr_append[] = "For more info, try '" PROG_NAME_STR " --help [flag name]' or '" PROG_NAME_STR " --help | more'.\n";

#ifndef NOLAPACK
static const char notestr_null_calc2[] = "Commands include --make-bpgen, --export, --freq, --geno-counts, --missing,\n--hardy, --indep-pairwise, --ld, --make-king, --king-cutoff, --write-samples,\n--write-snplist, --make-grm-list, --pca, --glm, --adjust-file, --score,\n--genotyping-rate, --validate, and --zst-decompress.\n\n'" PROG_NAME_STR " --help | more' describes all functions.\n";
#else
static const char notestr_null_calc2[] = "Commands include --make-bpgen, --export, --freq, --geno-counts, --missing,\n--hardy, --indep-pairwise, --ld, --make-king, --king-cutoff, --write-samples,\n--write-snplist, --make-grm-list, --glm, --adjust-file, --score,\n--genotyping-rate, --validate, and --zst-decompress.\n\n'" PROG_NAME_STR " --help | more' describes all functions.\n";
#endif

// covar-variance-standardize + terminating null
CONSTU31(kMaxFlagBlen, 27);

FLAGSET_DEF_START()
  kfLoadParams0,
  kfLoadParamsPgen = (1 << 0),
  kfLoadParamsPsam = (1 << 1),
  kfLoadParamsPvar = (1 << 2),
  kfLoadParamsPfileAll = (kfLoadParamsPgen | kfLoadParamsPsam | kfLoadParamsPvar)
FLAGSET_DEF_END(LoadParams);

FLAGSET_DEF_START()
  kfXload0,
  kfXloadVcf = (1 << 0),
  kfXloadBcf = (1 << 1),
  kfXloadOxSample = (1 << 2),
  kfXloadOxGen = (1 << 3),
  kfXloadOxBgen = (1 << 4),
  kfXloadOxHaps = (1 << 5),
  kfXloadOxLegend = (1 << 6),
  kfXloadPlink1Dosage = (1 << 7),
  kfXloadMap = (1 << 8),
  kfXloadGenDummy = (1 << 9)
FLAGSET_DEF_END(Xload);


// assume for now that .pgen must always be accompanied by both .pvar and .psam
// currently does double-duty in tracking dependencies
FLAGSET64_DEF_START()
  kfFilter0,
  kfFilterAllReq = (1 << 0),
  kfFilterPvarReq = (1 << 1),
  kfFilterPsamReq = (1 << 2),
  kfFilterNonrefFlagsNeeded = (1 << 3),
  kfFilterNoSplitChr = (1 << 4),
  kfFilterExclFemales = (1 << 5),
  kfFilterExclMales = (1 << 6),
  kfFilterExclNosex = (1 << 7),
  kfFilterExclFounders = (1 << 8),
  kfFilterExclNonfounders = (1 << 9),
  kfFilterSnpsOnly = (1 << 10),
  kfFilterSnpsOnlyJustAcgt = (1 << 11),
  kfFilterExtractIbed0 = (1 << 12),
  kfFilterExtractIbed1 = (1 << 13),
  kfFilterExcludeIbed0 = (1 << 14),
  kfFilterExcludeIbed1 = (1 << 15)
FLAGSET64_DEF_END(FilterFlags);

FLAGSET64_DEF_START()
  kfCommand10,
  kfCommand1MakePlink2 = (1 << 0),
  kfCommand1Exportf = (1 << 1),
  kfCommand1MakeKing = (1 << 2),
  kfCommand1KingCutoff = (1 << 3),
  kfCommand1MissingReport = (1 << 4),
  kfCommand1WriteSnplist = (1 << 5),
  kfCommand1AlleleFreq = (1 << 6),
  kfCommand1GenoCounts = (1 << 7),
  kfCommand1Hardy = (1 << 8),
  kfCommand1LdPrune = (1 << 9),
  kfCommand1Pca = (1 << 10),
  kfCommand1Glm = (1 << 11),
  kfCommand1MakeRel = (1 << 12),
  kfCommand1Validate = (1 << 13),
  kfCommand1GenotypingRate = (1 << 14),
  kfCommand1Score = (1 << 15),
  kfCommand1WriteCovar = (1 << 16),
  kfCommand1WriteSamples = (1 << 17),
  kfCommand1Ld = (1 << 18)
FLAGSET64_DEF_END(Command1Flags);

// this is a hybrid, only kfSortFileSid is actually a flag
FLAGSET_DEF_START()
  kfSort0,
  kfSortNone = (1 << 0),
  kfSortNatural = (1 << 1),
  kfSortAscii = (1 << 2),
  kfSortFile = (1 << 3),
  kfSortFileSid = (1 << 4)
FLAGSET_DEF_END(SortFlags);

typedef struct Plink2CmdlineStruct {
  MiscFlags misc_flags;

  // filter_flags tracks some info about flags which may cause the .pgen part
  // of --make-pgen's output to change (affects sample_include/variant_include,
  // etc.)  other .psam/.pvar/.pgen dependency info is now tracked by
  // dependency_flags (which is set to a superset of filter_flags at the end of
  // command-line parsing).
  // (filter_flags & kfFilterPsamReq) and (filter_flags & kfFilterPvarReq) can
  // now be used to detect whether sample_include or variant_include may be
  // modified.
  FilterFlags filter_flags;
  FilterFlags dependency_flags;

  Command1Flags command_flags1;
  PvarPsamFlags pvar_psam_flags;
  ExportfFlags exportf_flags;
  SortFlags sample_sort_flags;
  SortFlags sort_vars_flags;
  GrmFlags grm_flags;
  PcaFlags pca_flags;
  WriteCovarFlags write_covar_flags;
  PhenoTransformFlags pheno_transform_flags;
  RangeList snps_range_list;
  RangeList exclude_snps_range_list;
  RangeList pheno_range_list;
  RangeList covar_range_list;
  FamCol fam_cols;
  IdpasteFlags exportf_id_paste;
  UpdateSexInfo update_sex_info;
  LdInfo ld_info;
  KingFlags king_flags;
  double king_cutoff;
  double king_table_filter;
  double king_table_subset_thresh;
  FreqRptFlags freq_rpt_flags;
  MissingRptFlags missing_rpt_flags;
  GenoCountsFlags geno_counts_flags;
  HardyFlags hardy_flags;
  GlmInfo glm_info;
  AdjustInfo adjust_info;
  ScoreInfo score_info;
  APerm aperm;
  CmpExpr keep_if_expr;
  CmpExpr remove_if_expr;
  CmpExpr extract_if_info_expr;
  CmpExpr exclude_if_info_expr;
  double ci_size;
  float var_min_qual;
  uint32_t splitpar_bound1;
  uint32_t splitpar_bound2;
  uint32_t new_variant_id_max_allele_slen;
  uint32_t update_sex_colm2;

  // maybe support BGEN v1.2-style variable-precision dosages later, at which
  // point these should be floating-point numbers; but let's first see what we
  // gain from v1.1 fixed-point arithmetic
  uint32_t hard_call_thresh;
  uint32_t dosage_erase_thresh;

  double pfilter;
  double output_min_p;
  double vif_thresh;
  double mind_thresh;
  double geno_thresh;
  double hwe_thresh;
  double mach_r2_min;
  double mach_r2_max;
  double min_maf;
  double max_maf;
  double thin_keep_prob;
  double thin_keep_sample_prob;
  uint64_t min_allele_dosage;
  uint64_t max_allele_dosage;
  int32_t missing_pheno;
  int32_t from_bp;
  int32_t to_bp;
  int32_t window_bp;
  uint32_t pca_ct;
  uint32_t xchr_model;
  uint32_t max_thread_ct;
  uint32_t parallel_idx;
  uint32_t parallel_tot;
  uint32_t exportf_bits;
  uint32_t mwithin_val;
  uint32_t min_bp_space;
  uint32_t thin_keep_ct;
  uint32_t thin_keep_sample_ct;
  uint32_t keep_fcol_num;
  char exportf_id_delim;

  char* var_filter_exceptions_flattened;
  char* varid_template;
  char* missing_varid_match;
  char* varid_from;
  char* varid_to;
  char* varid_snp;
  char* varid_exclude_snp;
  char* pheno_fname;
  char* covar_fname;
  char* extract_fnames;
  char* exclude_fnames;
  char* keep_fnames;
  char* keepfam_fnames;
  char* remove_fnames;
  char* removefam_fnames;
  char* sample_sort_fname;
  char* freq_ref_binstr;
  char* freq_alt1_binstr;
  char* glm_local_covar_fname;
  char* glm_local_pvar_fname;
  char* glm_local_psam_fname;
  char* read_freq_fname;
  char* within_fname;
  char* catpheno_name;
  char* family_missing_catname;
  char* keep_cats_fname;
  char* keep_cat_names_flattened;
  char* keep_cat_phenoname;
  char* remove_cats_fname;
  char* remove_cat_names_flattened;
  char* remove_cat_phenoname;
  char* split_cat_phenonames_flattened;
  char* require_pheno_flattened;
  char* require_covar_flattened;
  char* vstd_flattened;
  char* quantnorm_flattened;
  char* covar_quantnorm_flattened;
  char* loop_cats_phenoname;
  char* ref_from_fa_fname;
  char* king_table_subset_fname;
  char* require_info_flattened;
  char* require_no_info_flattened;
  char* keep_fcol_fname;
  char* keep_fcol_flattened;
  char* keep_fcol_name;
  TwoColParams* ref_allele_flag;
  TwoColParams* alt1_allele_flag;
  TwoColParams* update_name_flag;
} Plink2Cmdline;

uint32_t SingleVariantLoaderIsNeeded(const char* king_cutoff_fprefix, Command1Flags command_flags1, MakePlink2Flags make_plink2_flags) {
  return (command_flags1 & (kfCommand1Exportf | kfCommand1MakeKing | kfCommand1GenoCounts | kfCommand1LdPrune | kfCommand1Validate | kfCommand1Pca | kfCommand1MakeRel | kfCommand1Glm | kfCommand1Score | kfCommand1Ld)) || ((command_flags1 & kfCommand1MakePlink2) && (make_plink2_flags & kfMakePgen)) || ((command_flags1 & kfCommand1KingCutoff) && (!king_cutoff_fprefix));
}


// might want to combine these two functions
uint32_t AlleleFreqsAreNeeded(Command1Flags command_flags1, GlmFlags glm_flags, double min_maf, double max_maf) {
  return (command_flags1 & (kfCommand1AlleleFreq | kfCommand1LdPrune | kfCommand1Pca | kfCommand1MakeRel | kfCommand1Score)) || (min_maf != 0.0) || (max_maf != 1.0) || ((command_flags1 & kfCommand1Glm) && (!(glm_flags & kfGlmA0Ref)));
}

uint32_t MajAllelesAreNeeded(Command1Flags command_flags1, GlmFlags glm_flags) {
  return (command_flags1 & (kfCommand1LdPrune | kfCommand1Pca | kfCommand1MakeRel)) || ((command_flags1 & kfCommand1Glm) && (!(glm_flags & kfGlmA0Ref)));
}


uint32_t GetFirstHaploidUidx(const ChrInfo* cip, UnsortedVar vpos_sortstatus) {
  // returns 0x7fffffff if no X/haploid chromosomes present
  if (!(vpos_sortstatus & kfUnsortedVarSplitChr)) {
    const uint32_t chr_ct = cip->chr_ct;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx < chr_ct; ++chr_fo_idx) {
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if (IsSet(cip->haploid_mask, chr_idx)) {
        return cip->chr_fo_vidx_start[chr_fo_idx];
      }
    }
  }
  return 0x7fffffff;
}

uint32_t AlleleDosagesAreNeeded(MiscFlags misc_flags, MakePlink2Flags make_plink2_flags, uint32_t afreq_needed, uint64_t min_allele_dosage, uint64_t max_allele_dosage) {
  return (make_plink2_flags & kfMakePlink2TrimAlts) || ((misc_flags & kfMiscNonfounders) && (afreq_needed || (misc_flags & kfMiscMajRef) || min_allele_dosage || (max_allele_dosage != UINT32_MAX)));
}

uint32_t FounderAlleleDosagesAreNeeded(MiscFlags misc_flags, uint32_t afreq_needed, uint64_t min_allele_dosage, uint64_t max_allele_dosage) {
  return (afreq_needed || (misc_flags & kfMiscMajRef) || min_allele_dosage || (max_allele_dosage != (~0LLU))) && (!(misc_flags & kfMiscNonfounders));
}

uint32_t SampleMissingDosageCtsAreNeeded(MiscFlags misc_flags, uint32_t smaj_missing_geno_report_requested, double mind_thresh, MissingRptFlags missing_rpt_flags) {
  return ((mind_thresh != 1.0) && (misc_flags & kfMiscMindDosage)) || (smaj_missing_geno_report_requested && (missing_rpt_flags & (kfMissingRptScolNmissDosage | kfMissingRptScolFmissDosage)));
}

uint32_t VariantMissingHcCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double geno_thresh, MissingRptFlags missing_rpt_flags) {
  return ((command_flags1 & kfCommand1GenotypingRate) && (!(misc_flags & kfMiscGenotypingRateDosage))) || ((command_flags1 & kfCommand1MissingReport) && (missing_rpt_flags & (kfMissingRptVcolNmiss | kfMissingRptVcolNmissHh | kfMissingRptVcolHethap | kfMissingRptVcolFmiss | kfMissingRptVcolFmissHh | kfMissingRptVcolFhethap))) || ((geno_thresh != 1.0) && (!(misc_flags & kfMiscGenoDosage)));
}

uint32_t VariantHethapCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double geno_thresh, MissingRptFlags missing_rpt_flags, uint32_t first_hap_uidx) {
  return (first_hap_uidx != 0x7fffffff) && (((command_flags1 & kfCommand1MissingReport) && (missing_rpt_flags & (kfMissingRptVcolNmissHh | kfMissingRptVcolHethap | kfMissingRptVcolFmissHh | kfMissingRptVcolFhethap))) || ((geno_thresh != 1.0) && (!(misc_flags & kfMiscGenoHhMissing))));
}

uint32_t VariantMissingDosageCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double geno_thresh, MissingRptFlags missing_rpt_flags) {
  return ((command_flags1 & kfCommand1GenotypingRate) && (misc_flags & kfMiscGenotypingRateDosage)) || ((command_flags1 & kfCommand1MissingReport) && (!(missing_rpt_flags & kfMissingRptSampleOnly)) && (missing_rpt_flags & (kfMissingRptVcolNmissDosage | kfMissingRptVcolFmissDosage))) || ((geno_thresh != 1.0) && (misc_flags & kfMiscGenoDosage));
}

// can simplify --geno-counts all-biallelic case, but let's first make sure the
// general case works for multiallelic variants
uint32_t RawGenoCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double hwe_thresh) {
  return (command_flags1 & kfCommand1GenoCounts) || ((misc_flags & kfMiscNonfounders) && ((command_flags1 & kfCommand1Hardy) || (hwe_thresh != 1.0)));
}

uint32_t FounderRawGenoCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double hwe_thresh) {
  return (!(misc_flags & kfMiscNonfounders)) && ((command_flags1 & kfCommand1Hardy) || (hwe_thresh != 1.0));
}

uint32_t InfoReloadIsNeeded(Command1Flags command_flags1, PvarPsamFlags pvar_psam_flags, ExportfFlags exportf_flags) {
  // add kfExportfBcf later
  return ((command_flags1 & kfCommand1MakePlink2) && (pvar_psam_flags & kfPvarColXinfo)) || ((command_flags1 & kfCommand1Exportf) && (exportf_flags & kfExportfVcf));
}

uint32_t GrmKeepIsNeeded(Command1Flags command_flags1, PcaFlags pca_flags) {
  return ((command_flags1 & kfCommand1Pca) && (!(pca_flags & kfPcaApprox)));
}

void ReportGenotypingRate(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_missing_cts, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t variant_ct, uint32_t is_dosage) {
  // defined the same way as PLINK 1.x, to allow this to serve as a sanity
  // check
  // trivial to multithread this if it ever matters
  uint64_t tot_nony_missing = 0;
  uint64_t tot_y_missing = 0;
  uint64_t cur_tot_missing = 0;
  uint32_t y_start = UINT32_MAX;
  uint32_t y_end = UINT32_MAX;
  uint32_t variant_ct_y = 0;
  uint32_t y_code;
  if (XymtExists(cip, kChrOffsetY, &y_code)) {
    const uint32_t y_chr_fo_idx = cip->chr_idx_to_foidx[y_code];
    y_start = cip->chr_fo_vidx_start[y_chr_fo_idx];
    y_end = cip->chr_fo_vidx_start[y_chr_fo_idx + 1];
    variant_ct_y = PopcountBitRange(variant_include, y_start, y_end);
  }
  uint32_t y_thresh = y_start;
  uint32_t variant_uidx = 0;
  uint32_t is_y = 0;
  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
    MovU32To1Bit(variant_include, &variant_uidx);
    if (variant_uidx >= y_thresh) {
      if (is_y) {
        tot_y_missing = cur_tot_missing;
      } else {
        tot_nony_missing = cur_tot_missing;
      }
      is_y = (variant_uidx < y_end);
      cur_tot_missing = 0;
      if (is_y) {
        y_thresh = y_end;
      } else {
        y_thresh = UINT32_MAX;
      }
    }
    cur_tot_missing += variant_missing_cts[variant_uidx];
  }
  if (is_y) {
    tot_y_missing = cur_tot_missing;
  } else {
    tot_nony_missing += cur_tot_missing;
  }
  if ((!tot_y_missing) && (!tot_nony_missing)) {
    logprintf("Total (%s) genotyping rate %sis exactly 1.\n", is_dosage? "dosage" : "hardcall", (raw_sample_ct != sample_ct)? "in remaining samples " : "");
    return;
  }
  double genotyping_rate;
  if (male_ct && variant_ct_y) {
    const uint64_t nony_possible_obs = (variant_ct - variant_ct_y) * S_CAST(uint64_t, sample_ct);
    const uint64_t y_possible_obs = variant_ct_y * S_CAST(uint64_t, male_ct);
    genotyping_rate = u63tod(nony_possible_obs - tot_nony_missing) / u31tod(sample_ct) + u63tod(y_possible_obs - tot_y_missing) / u31tod(male_ct);
    genotyping_rate /= u31tod(variant_ct);
  } else {
    variant_ct -= variant_ct_y;
    const uint64_t denom = variant_ct * S_CAST(uint64_t, sample_ct);
    genotyping_rate = u63tod(denom - tot_nony_missing) / u63tod(denom);
  }
  if (genotyping_rate >= 0.9999995) {
    logprintf("Total (%s) genotyping rate %sis in [0.9999995, 1).\n", is_dosage? "dosage" : "hardcall", (raw_sample_ct != sample_ct)? "in remaining samples " : "");
  } else {
    logprintf("Total (%s) genotyping rate %sis %g.\n", is_dosage? "dosage" : "hardcall", (raw_sample_ct != sample_ct)? "in remaining samples " : "", genotyping_rate);
  }
}

PglErr ApplyVariantBpFilters(const char* extract_fnames, const char* exclude_fnames, const ChrInfo* cip, const uint32_t* variant_bps, int32_t from_bp, int32_t to_bp, uint32_t raw_variant_ct, FilterFlags filter_flags, UnsortedVar vpos_sortstatus, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  if (!(*variant_ct_ptr)) {
    return kPglRetSuccess;
  }
  if ((from_bp != -1) || (to_bp != -1)) {
    if (vpos_sortstatus & kfUnsortedVarBp) {
      logerrputs("Error: --from-bp and --to-bp require a sorted .pvar/.bim.  Retry this command\nafter using --make-pgen/--make-bed + --sort-vars to sort your data.\n");
      return kPglRetInconsistentInput;
    }
    const uint32_t chr_idx = AdvBoundedTo1Bit(cip->chr_mask, 0, kChrRawEnd);

    // this function shouldn't be called unless variant_ct is nonzero
    assert(chr_idx != kChrRawEnd);

    const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
    uint32_t variant_uidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
    uint32_t variant_uidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    if (from_bp != -1) {
      const uint32_t from_offset = CountSortedSmallerU32(&(variant_bps[variant_uidx_start]), variant_uidx_end - variant_uidx_start, from_bp);
      variant_uidx_start += from_offset;
    }
    if ((to_bp != -1) && (variant_uidx_start < variant_uidx_end)) {
      const uint32_t to_offset = CountSortedSmallerU32(&(variant_bps[variant_uidx_start]), variant_uidx_end - variant_uidx_start, 1 + to_bp);
      variant_uidx_end = variant_uidx_start + to_offset;
    }
    if (variant_uidx_start) {
      ClearBitsNz(0, variant_uidx_start, variant_include);
    }
    if (variant_uidx_end < raw_variant_ct) {
      ClearBitsNz(variant_uidx_end, raw_variant_ct, variant_include);
    }
    *variant_ct_ptr = PopcountBitRange(variant_include, variant_uidx_start, variant_uidx_end);
  }
  if (extract_fnames && (filter_flags & (kfFilterExtractIbed0 | kfFilterExtractIbed1))) {
    if (vpos_sortstatus & kfUnsortedVarBp) {
      logerrputs("Error: '--extract ibed0'/'--extract ibed1' requires a sorted .pvar/.bim.  Retry\nthis command after using --make-pgen/--make-bed + --sort-vars to sort your\ndata.\n");
      return kPglRetInconsistentInput;
    }
    PglErr reterr = ExtractExcludeRange(extract_fnames, cip, variant_bps, raw_variant_ct, 0, (filter_flags / kfFilterExtractIbed0) & 1, variant_include, variant_ct_ptr);
    if (reterr) {
      return reterr;
    }
  }
  if (exclude_fnames && (filter_flags & (kfFilterExcludeIbed0 | kfFilterExcludeIbed1))) {
    if (vpos_sortstatus & kfUnsortedVarBp) {
      logerrputs("Error: '--exclude ibed0'/'--exclude ibed1' requires a sorted .pvar/.bim.  Retry\nthis commandafter using --make-pgen/--make-bed + --sort-vars to sort your\ndata.\n");
      return kPglRetInconsistentInput;
    }
    PglErr reterr = ExtractExcludeRange(exclude_fnames, cip, variant_bps, raw_variant_ct, 1, (filter_flags / kfFilterExcludeIbed0) & 1, variant_include, variant_ct_ptr);
    if (reterr) {
      return reterr;
    }
  }
  return kPglRetSuccess;
}

void UpdateSampleSubsets(const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t* founder_info, uint32_t* founder_ct_ptr, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t* male_ct_ptr, uint32_t* nosex_ct_ptr) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  BitvecAnd(sample_include, raw_sample_ctl, founder_info);
  *founder_ct_ptr = PopcountWords(founder_info, raw_sample_ctl);
  BitvecAnd(sample_include, raw_sample_ctl, sex_male);
  *male_ct_ptr = PopcountWords(sex_male, raw_sample_ctl);
  BitvecAnd(sample_include, raw_sample_ctl, sex_nm);
  *nosex_ct_ptr = sample_ct - PopcountWords(sex_nm, raw_sample_ctl);
}

// command_flags2 will probably be needed before we're done
static_assert(kPglMaxAltAlleleCt == 254, "Plink2Core() --maj-ref needs to be updated.");
PglErr Plink2Core(const Plink2Cmdline* pcp, MakePlink2Flags make_plink2_flags, char* pgenname, char* psamname, char* pvarname, char* outname, char* outname_end, char* king_cutoff_fprefix, ChrInfo* cip) {
  PhenoCol* pheno_cols = nullptr;
  PhenoCol* covar_cols = nullptr;
  PhenoCol* loop_cats_pheno_col = nullptr;
  char* pheno_names = nullptr;
  char* covar_names = nullptr;
  uint32_t pheno_ct = 0;
  uint32_t covar_ct = 0;
  PglErr reterr = kPglRetSuccess;
  PgenFileInfo pgfi;
  PgenReader simple_pgr;
  PreinitPgfi(&pgfi);
  PreinitPgr(&simple_pgr);
  {
    // this predicate will need to exclude --merge-list special case later
    uint32_t pvar_renamed = 0;
    if ((make_plink2_flags & (kfMakeBed | kfMakePgen)) || (pcp->exportf_flags & kfExportfIndMajorBed)) {
      uint32_t fname_slen;
#ifdef _WIN32
      fname_slen = GetFullPathName(pgenname, kPglFnamesize, g_textbuf, nullptr);
      if ((!fname_slen) || (fname_slen > kPglFnamesize))
#else
      if (!realpath(pgenname, g_textbuf))
#endif
      {
        logerrprintfww(kErrprintfFopen, pgenname);
        goto Plink2Core_ret_OPEN_FAIL;
      }
      uint32_t pgen_rename = 0;
      if (make_plink2_flags & kfMakePgen) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
        pgen_rename = RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
      }
      if ((!pgen_rename) && ((make_plink2_flags & kfMakeBed) || (pcp->exportf_flags & kfExportfIndMajorBed))) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".bed");
        pgen_rename = RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
      }
      if (pgen_rename) {
        logerrprintf("Warning: --make-%s input and output filenames match.  Appending '~' to input\nfilenames.\n", (make_plink2_flags & kfMakeBed)? "bed" : ((make_plink2_flags & kfMakePvar)? "pgen" : "bpgen"));
        fname_slen = strlen(pgenname);
        memcpy(g_textbuf, pgenname, fname_slen);
        snprintf(&(pgenname[fname_slen]), 2, "~");
        if (rename(g_textbuf, pgenname)) {
          logerrputs("Error: Failed to append '~' to input .bed/.pgen filename.\n");
          goto Plink2Core_ret_OPEN_FAIL;
        }
        fname_slen = strlen(pvarname);
        memcpy(g_textbuf, pvarname, fname_slen);
        snprintf(&(pvarname[fname_slen]), 2, "~");
        if (rename(g_textbuf, pvarname)) {
          logerrputs("Error: Failed to append '~' to input .bim/.pvar filename.\n");
          goto Plink2Core_ret_OPEN_FAIL;
        }
        pvar_renamed = 1;
        fname_slen = strlen(psamname);
        memcpy(g_textbuf, psamname, fname_slen);
        snprintf(&(psamname[fname_slen]), 2, "~");
        if (rename(g_textbuf, psamname)) {
          logerrputs("Error: Failed to append '~' to input .fam/.psam filename.\n");
          goto Plink2Core_ret_OPEN_FAIL;
        }
      }
    }
    uint32_t raw_sample_ct = 0;
    uintptr_t* sample_include = nullptr;
    uintptr_t* sex_nm = nullptr;
    uintptr_t* sex_male = nullptr;
    uintptr_t* founder_info = nullptr;
    uintptr_t max_pheno_name_blen = 0;
    uint32_t raw_sample_ctl = 0;
    uint32_t sample_ct = 0;

    // There's a tradeoff between using structs like this and passing each of
    // the components separately: large structs make for shorter parameter
    // lists and hence slightly prettier code, but going too far in that
    // direction (e.g. having a single huge struct that gets passed to all
    // commands) makes it far harder to see what the actual data dependencies
    // are.  Right now the codebase errs, to an almost comical degree, in the
    // direction of making data dependencies clear and minimizing passing of
    // unneeded information.  I'll scale this back, but I'll try to avoid doing
    // it in ways that obscure data flows.
    // For instance, sample IDs and parental IDs could be stored in the same
    // struct, but if we pass parental IDs separately it's easier to see which
    // functions actually need them.
    // In lower-level functions, a reasonable rule of thumb is that aggregation
    // of parameters that normally appear together is worthwhile when there are
    // more than 6 parameters, since the x86-64 C calling convention uses
    // registers for the first 6 parameters.
    PedigreeIdInfo pii;
    InitPedigreeIdInfo(pcp->misc_flags, &pii);
    if (psamname[0]) {
      // update (26 Nov 2017): change --no-pheno to also apply to .psam file
      const uint32_t ignore_psam_phenos = (!(pcp->fam_cols & kfFamCol6)) || (pcp->pheno_fname && pcp->pheno_range_list.name_ct);
      reterr = LoadPsam(psamname, pcp->pheno_fname? nullptr : &(pcp->pheno_range_list), pcp->fam_cols, ignore_psam_phenos? 0 : 0x7fffffff, pcp->missing_pheno, (pcp->misc_flags / kfMiscAffection01) & 1, &pii, &sample_include, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
      // todo: move this check after --update-ids once that's implemented.
      if ((pii.sii.flags & kfSampleIdFidPresent) && ((pii.sii.flags & kfSampleIdNoIdHeaderIidOnly) || (pcp->grm_flags & kfGrmNoIdHeaderIidOnly))) {
        for (uint32_t sample_idx = 0; sample_idx < raw_sample_ct; ++sample_idx) {
          const char* cur_sample_id = &(pii.sii.sample_ids[sample_idx * pii.sii.max_sample_id_blen]);
          if (memcmp(cur_sample_id, "0\t", 2)) {
            logerrputs("Error: 'iid-only' modifier can only be used when FIDs are missing or all-0.\n");
            goto Plink2Core_ret_INCONSISTENT_INPUT;
          }
        }
      }

      // todo: add option to discard loaded SIDs
      raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
      sample_ct = PopcountWords(sample_include, raw_sample_ctl);
      const uint32_t known_sex_ct = PopcountWords(sex_nm, raw_sample_ctl);
      const uint32_t male_ct = PopcountWords(sex_male, raw_sample_ctl);
      const uint32_t female_ct = known_sex_ct - male_ct;
      const uint32_t founder_ct = PopcountWords(founder_info, raw_sample_ctl);
      if (known_sex_ct == sample_ct) {
        logprintfww("%u sample%s (%u female%s, %u male%s; %u founder%s) loaded from %s.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", founder_ct, (founder_ct == 1)? "" : "s", psamname);
      } else {
        const uint32_t unknown_sex_ct = sample_ct - known_sex_ct;
        logprintfww("%u sample%s (%u female%s, %u male%s, %u ambiguous; %u founder%s) loaded from %s.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", unknown_sex_ct, founder_ct, (founder_ct == 1)? "" : "s", psamname);
      }
    }

    uint32_t max_variant_id_slen = 1;
    uint32_t info_reload_slen = InfoReloadIsNeeded(pcp->command_flags1, pcp->pvar_psam_flags, pcp->exportf_flags);
    uintptr_t* variant_allele_idxs = nullptr;
    uint32_t raw_variant_ct = 0;
    uint32_t variant_ct = 0;
    char* xheader = nullptr;
    uintptr_t xheader_blen = 0;
    uintptr_t* variant_include = nullptr;
    uint32_t* variant_bps = nullptr;
    char** variant_ids_mutable = nullptr;

    // This is actually a hybrid of const char** and char**: length >1
    // pointed-to strings can be modified in-place by --update-alleles, but we
    // must be able to assign const pointers to length-1 strings to it.  So we
    // set the official type to char**, but there's more const_cast-ing than
    // usual.
    char** allele_storage_mutable = nullptr;

    uintptr_t* pvar_qual_present = nullptr;
    float* pvar_quals = nullptr;
    uintptr_t* pvar_filter_present = nullptr;
    uintptr_t* pvar_filter_npass = nullptr;
    const char* const* pvar_filter_storage = nullptr;
    uintptr_t* nonref_flags = nullptr;
    InfoFlags info_flags = kfInfo0;
    uint32_t max_allele_slen = 0;
    uint32_t max_filter_slen = 0;
    UnsortedVar vpos_sortstatus = kfUnsortedVar0;
    double* variant_cms = nullptr;
    ChrIdx* chr_idxs = nullptr;  // split-chromosome case only
    if (pvarname[0]) {
      char** pvar_filter_storage_mutable = nullptr;
      reterr = LoadPvar(pvarname, pcp->var_filter_exceptions_flattened, pcp->varid_template, pcp->missing_varid_match, pcp->require_info_flattened, pcp->require_no_info_flattened, pcp->extract_if_info_expr, pcp->exclude_if_info_expr, pcp->misc_flags, pcp->pvar_psam_flags, pcp->exportf_flags, pcp->var_min_qual, pcp->splitpar_bound1, pcp->splitpar_bound2, pcp->new_variant_id_max_allele_slen, (pcp->filter_flags / kfFilterSnpsOnly) & 3, !(pcp->dependency_flags & kfFilterNoSplitChr), pcp->max_thread_ct, cip, &max_variant_id_slen, &info_reload_slen, &vpos_sortstatus, &xheader, &variant_include, &variant_bps, &variant_ids_mutable, &variant_allele_idxs, K_CAST(const char***, &allele_storage_mutable), &pvar_qual_present, &pvar_quals, &pvar_filter_present, &pvar_filter_npass, &pvar_filter_storage_mutable, &nonref_flags, &variant_cms, &chr_idxs, &raw_variant_ct, &variant_ct, &max_allele_slen, &xheader_blen, &info_flags, &max_filter_slen);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
      pvar_filter_storage = TO_CONSTCPCONSTP(pvar_filter_storage_mutable);
      if (variant_ct == raw_variant_ct) {
        logprintfww("%u variant%s loaded from %s.\n", variant_ct, (variant_ct == 1)? "" : "s", pvarname);
      } else {
        logprintfww("%u out of %u variant%s loaded from %s.\n", variant_ct, raw_variant_ct, (raw_variant_ct == 1)? "" : "s", pvarname);
      }
      if (info_reload_slen && (make_plink2_flags & (kfMakeBim | kfMakePvar)) && (!pvar_renamed)) {
        // need to be careful with .pvar in this case
        uint32_t fname_slen;
#ifdef _WIN32
        fname_slen = GetFullPathName(pvarname, kPglFnamesize, g_textbuf, nullptr);
        if ((!fname_slen) || (fname_slen > kPglFnamesize))
#else
        if (!realpath(pvarname, g_textbuf))
#endif
        {
          logerrprintfww(kErrprintfFopen, pvarname);
          goto Plink2Core_ret_OPEN_FAIL;
        }
        if (make_plink2_flags & kfMakeBim) {
          OutnameZstSet(".bim", make_plink2_flags & kfMakeBimZs, outname_end);
          pvar_renamed = RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
          if (pvar_renamed) {
            logerrputs("Warning: .bim input and output filenames match.  Appending '~' to input\nfilename.\n");
            fname_slen = strlen(pvarname);
            memcpy(g_textbuf, pvarname, fname_slen);
            snprintf(&(pvarname[fname_slen]), 2, "~");
            if (rename(g_textbuf, pvarname)) {
              logerrputs("Error: Failed to append '~' to input .bim filename.\n");
              goto Plink2Core_ret_OPEN_FAIL;
            }
          }
        }
        if ((!pvar_renamed) && (make_plink2_flags & kfMakePvar)) {
          OutnameZstSet(".pvar", pcp->pvar_psam_flags & kfPvarZs, outname_end);
          // pvar_renamed = RealpathIdentical();
          if (RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]))) {
            logerrputs("Warning: .pvar input and output filenames match.  Appending '~' to input\nfilename.\n");
            fname_slen = strlen(pvarname);
            memcpy(g_textbuf, pvarname, fname_slen);
            snprintf(&(pvarname[fname_slen]), 2, "~");
            if (rename(g_textbuf, pvarname)) {
              logerrputs("Error: Failed to append '~' to input .pvar filename.\n");
              goto Plink2Core_ret_OPEN_FAIL;
            }
          }
        }
      }
      if ((pcp->dependency_flags & (kfFilterAllReq | kfFilterNonrefFlagsNeeded)) == kfFilterNonrefFlagsNeeded) {
        // We need INFO:PR, but we don't really care about anything else
        // potentially in the .pgen.
        // So if INFO was present in the variant file, we don't need to look at
        // the .pgen at all (even if there was no PR key; we infer that all
        // reference alleles are accurate in that case).
        // Conversely, if no .pgen filename was provided, and the variant file
        // doesn't contain an INFO field, print a warning.
        if (info_flags & kfInfoAll) {
          pgenname[0] = '\0';
        } else if (!pgenname[0]) {
          logerrputs("Warning: Variant file does not distinguish between provisional and trusted REF\nalleles, and no .pgen input file was provided.\n");
          if (info_flags & kfInfoPrNonrefDefault) {
            logerrputs("Assuming all REF alleles are provisional.\n");
          } else {
            logerrputs("Assuming all REF alleles are trusted.\n");
          }
        }
      }
    }

    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t pgr_alloc_cacheline_ct = 0;
    if (pgenname[0]) {
      PgenHeaderCtrl header_ctrl;
      uintptr_t cur_alloc_cacheline_ct;
      while (1) {
        reterr = PgfiInitPhase1(pgenname, raw_variant_ct, raw_sample_ct, 0, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, g_logbuf);
        if (!reterr) {
          break;
        }
        // detect and autoconvert plink 1 sample-major files, instead of
        // failing (don't bother supporting plink 0.99 files any more)
        if (reterr == kPglRetSampleMajorBed) {
          char* pgenname_end = memcpya(pgenname, outname, outname_end - outname);
          pgenname_end = strcpya(pgenname_end, ".pgen");
          const uint32_t no_vmaj_ext = (pcp->command_flags1 & kfCommand1MakePlink2) && (!pcp->filter_flags) && ((make_plink2_flags & (kfMakePgen | (kfMakePgenFormatBase * 3))) == kfMakePgen);
          if (no_vmaj_ext) {
            *pgenname_end = '\0';
            make_plink2_flags &= ~kfMakePgen;
            // no --make-just-pgen command, so we'll never entirely skip the
            // make_plink2 operation
          } else {
            snprintf(pgenname_end, kMaxOutfnameExtBlen - 5, ".vmaj");
          }
          reterr = Plink1SampleMajorToPgen(pgenname, raw_variant_ct, raw_sample_ct, (pcp->misc_flags / kfMiscRealRefAlleles) & 1, pcp->max_thread_ct, pgfi.shared_ff);
          if (!reterr) {
            fclose(pgfi.shared_ff);
            pgfi.shared_ff = nullptr;
            continue;
          }
        } else {
          if (reterr != kPglRetReadFail) {
            WordWrapB(0);
            logerrputsb();
          }
        }
        goto Plink2Core_ret_1;
      }
      pgfi.allele_idx_offsets = variant_allele_idxs;
      unsigned char* pgfi_alloc;
      if (bigstack_alloc_uc(cur_alloc_cacheline_ct * kCacheline, &pgfi_alloc)) {
        goto Plink2Core_ret_NOMEM;
      }
      const uint32_t nonref_flags_already_loaded = (nonref_flags != nullptr);
      if ((!nonref_flags) && ((header_ctrl & 192) == 192)) {
        if (bigstack_alloc_w(raw_variant_ctl, &nonref_flags)) {
          goto Plink2Core_ret_NOMEM;
        }
      }
      pgfi.nonref_flags = nonref_flags;
      uint32_t max_vrec_width;
      // only practical effect of setting use_blockload to zero here is that
      // pgr_alloc_cacheline_ct is overestimated by
      // DivUp(max_vrec_width, kCacheline).
      reterr = PgfiInitPhase2(header_ctrl, 1, nonref_flags_already_loaded, 1, 0, raw_variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, &pgr_alloc_cacheline_ct, g_logbuf);
      if (reterr) {
        if (reterr != kPglRetReadFail) {
          WordWrapB(0);
          logerrputsb();
        }
        goto Plink2Core_ret_1;
      }
      if (pcp->misc_flags & kfMiscRealRefAlleles) {
        if (nonref_flags && (!AllBitsAreOne(nonref_flags, raw_variant_ct))) {
          // technically a lie, it's okay if a .bed is first converted to .pgen
          // without this flag, and then the user remembers the existence of
          // --real-ref-alleles later.  but to reduce the ease of
          // foot-shooting, we don't allow this to clobber arbitrary
          // nonref_flags arrays.
          logerrputs("Error: --real-ref-alleles must be used on a plink1 fileset.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }

        // wasteful if nonref_flags was allocated, but nonref_flags isn't that
        // large, and --real-ref-alleles + --make-pgen can be run separately
        // from anything truly memory-limited
        nonref_flags = nullptr;
        pgfi.nonref_flags = nullptr;

        pgfi.gflags &= ~kfPgenGlobalAllNonref;
      }
      if (SingleVariantLoaderIsNeeded(king_cutoff_fprefix, pcp->command_flags1, make_plink2_flags)) {
        // ugly kludge, probably want to add pgenlib_internal support for this
        // hybrid use pattern
        FILE* shared_ff_copy = pgfi.shared_ff;
        pgfi.shared_ff = nullptr;
        unsigned char* simple_pgr_alloc;
        if (bigstack_alloc_uc((pgr_alloc_cacheline_ct + DivUp(max_vrec_width, kCacheline)) * kCacheline, &simple_pgr_alloc)) {
          goto Plink2Core_ret_NOMEM;
        }
        reterr = PgrInit(pgenname, max_vrec_width, &pgfi, &simple_pgr, simple_pgr_alloc);
        if (reterr) {
          if (reterr == kPglRetOpenFail) {
            logerrprintf(kErrprintfFopen, pgenname);
          }
          // only other possibility is kPglRetReadFail
          goto Plink2Core_ret_1;
        }
        pgfi.shared_ff = shared_ff_copy;
        if (pcp->command_flags1 & kfCommand1Validate) {
          logprintfww5("Validating %s... ", pgenname);
          fflush(stdout);
          reterr = PgrValidate(&simple_pgr, g_logbuf);
          if (reterr) {
            if (reterr != kPglRetReadFail) {
              logputs("\n");
              WordWrapB(0);
              logerrputsb();
            }
            goto Plink2Core_ret_1;
          }
          logputs("done.\n");
          if (!(pcp->command_flags1 & (~kfCommand1Validate))) {
            goto Plink2Core_ret_1;
          }
        }
      }
      // any functions using blockload must perform its own PgrInit(), etc.
    } else {
      // bugfix (10-11 Feb 2018): these variables may still be accessed
      pgfi.gflags = S_CAST(PgenGlobalFlags, ((info_flags / kfInfoPrNonrefDefault) & 1) * kfPgenGlobalAllNonref);
      pgfi.nonref_flags = nonref_flags;
    }
    if (pcp->pheno_fname) {
      reterr = LoadPhenos(pcp->pheno_fname, &(pcp->pheno_range_list), sample_include, pii.sii.sample_ids, raw_sample_ct, sample_ct, pii.sii.max_sample_id_blen, pcp->missing_pheno, (pcp->misc_flags / kfMiscAffection01) & 1, (pcp->misc_flags / kfMiscPhenoColNums) & 1, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
    }

    // move processing of PLINK 1.x cluster-loading/filtering flags here, since
    // they're now under the categorical-phenotype umbrella
    if ((pcp->misc_flags & kfMiscCatPhenoFamily) || pcp->within_fname) {
      reterr = Plink1ClusterImport(pcp->within_fname, pcp->catpheno_name, pcp->family_missing_catname, sample_include, pii.sii.sample_ids, raw_sample_ct, sample_ct, pii.sii.max_sample_id_blen, pcp->mwithin_val, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
    }

    if (!pheno_ct) {
      logputs("Note: No phenotype data present.\n");
    } else {
      if (pheno_ct == 1) {
        const uint32_t obs_ct = PopcountWords(pheno_cols[0].nonmiss, raw_sample_ctl);
        if (pheno_cols[0].type_code == kPhenoDtypeCc) {
          const uint32_t case_ct = PopcountWords(pheno_cols[0].data.cc, raw_sample_ctl);
          const uint32_t ctrl_ct = obs_ct - case_ct;
          logprintf("1 binary phenotype loaded (%u case%s, %u control%s).\n", case_ct, (case_ct == 1)? "" : "s", ctrl_ct, (ctrl_ct == 1)? "" : "s");
        } else if (pheno_cols[0].type_code == kPhenoDtypeQt) {
          logprintf("1 quantitative phenotype loaded (%u value%s).\n", obs_ct, (obs_ct == 1)? "" : "s");
        } else {
          logprintf("1 categorical phenotype loaded (%u value%s).\n", obs_ct, (obs_ct == 1)? "" : "s");
        }
      } else {
        uint32_t cc_ct = 0;
        uint32_t qt_ct = 0;
        for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
          const PhenoDtype cur_type_code = pheno_cols[pheno_idx].type_code;
          if (pheno_cols[pheno_idx].type_code == kPhenoDtypeCc) {
            ++cc_ct;
          } else if (cur_type_code == kPhenoDtypeQt) {
            ++qt_ct;
          }
        }
        uint32_t cat_ct = pheno_ct - cc_ct - qt_ct;
        // just brute force this for now
        if (!cc_ct) {
          if (!qt_ct) {
            logprintf("%u categorical phenotypes loaded.\n", pheno_ct);
          } else if (!cat_ct) {
            logprintf("%u quantitative phenotypes loaded.\n", pheno_ct);
          } else {
            logprintf("%u phenotypes loaded (%u quantitative, %u categorical).\n", pheno_ct, qt_ct, cat_ct);
          }
        } else if (!qt_ct) {
          if (!cat_ct) {
            logprintf("%u binary phenotypes loaded.\n", pheno_ct);
          } else {
            logprintf("%u phenotypes loaded (%u binary, %u categorical).\n", pheno_ct, cc_ct, cat_ct);
          }
        } else if (!cat_ct) {
          logprintf("%u phenotypes loaded (%u binary, %u quantitative).\n", pheno_ct, cc_ct, qt_ct);
        } else {
          logprintfww("%u phenotypes loaded (%u binary, %u quantitative, %u categorical).\n", pheno_ct, cc_ct, qt_ct, cat_ct);
        }
      }
    }
    // If something like --snps is combined with a position-based filter which
    // may remove some of the named variants, we need to apply --snps first.
    // Otherwise, it may be very advantageous to apply the position-based
    // filters before constructing the variant ID hash table.  So we split this
    // into two cases.
    const uint32_t full_variant_id_htable_needed = variant_ct && (pcp->varid_from || pcp->varid_to || pcp->varid_snp || pcp->varid_exclude_snp || pcp->snps_range_list.name_ct || pcp->exclude_snps_range_list.name_ct || pcp->update_name_flag);
    if (!full_variant_id_htable_needed) {
      reterr = ApplyVariantBpFilters(pcp->extract_fnames, pcp->exclude_fnames, cip, variant_bps, pcp->from_bp, pcp->to_bp, raw_variant_ct, pcp->filter_flags, vpos_sortstatus, variant_include, &variant_ct);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
    }
    if (variant_ct && (full_variant_id_htable_needed || (pcp->extract_fnames && (!(pcp->filter_flags & (kfFilterExtractIbed0 | kfFilterExtractIbed1)))) || (pcp->exclude_fnames && (!(pcp->filter_flags & (kfFilterExcludeIbed0 | kfFilterExcludeIbed1)))))) {
      // don't bother with having different allow_dups vs. no allow_dups hash
      // table structures, just check specific IDs for duplication in the
      // no-duplicates-allowed cases
      unsigned char* bigstack_mark = g_bigstack_base;
      uint32_t* variant_id_htable = nullptr;
      uint32_t* htable_dup_base = nullptr;
      uint32_t variant_id_htable_size;
      reterr = AllocAndPopulateIdHtableMt(variant_include, TO_CONSTCPCONSTP(variant_ids_mutable), variant_ct, pcp->max_thread_ct, &variant_id_htable, &htable_dup_base, &variant_id_htable_size);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
      if (vpos_sortstatus & kfUnsortedVarBp) {
        if (pcp->varid_from || pcp->varid_to) {
          logerrputs("Error: --from/--to require a sorted .pvar/.bim.  Retry this command after using\n--make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        if (pcp->window_bp != -1) {
          logerrputs("Error: --window requires a sorted .pvar/.bim.  Retry this command after using\n--make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
      }
      if (pcp->varid_from || pcp->varid_to) {
        reterr = FromToFlag(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, pcp->varid_from, pcp->varid_to, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, variant_include, cip, &variant_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->varid_snp) {
        reterr = SnpFlag(variant_bps, TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, pcp->varid_snp, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 0, pcp->window_bp, variant_include, cip, &variant_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->snps_range_list.name_ct) {
        reterr = SnpsFlag(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, &(pcp->snps_range_list), raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 0, variant_include, &variant_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->varid_exclude_snp) {
        reterr = SnpFlag(variant_bps, TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, pcp->varid_exclude_snp, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 1, pcp->window_bp, variant_include, cip, &variant_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->exclude_snps_range_list.name_ct) {
        reterr = SnpsFlag(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, &(pcp->exclude_snps_range_list), raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 1, variant_include, &variant_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->update_name_flag && variant_ct) {
        reterr = UpdateVarNames(variant_include, variant_id_htable, pcp->update_name_flag, raw_variant_ct, variant_id_htable_size, variant_ids_mutable, &max_variant_id_slen);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
        if ((pcp->extract_fnames && (!(pcp->filter_flags & (kfFilterExtractIbed0 | kfFilterExtractIbed1)))) || (pcp->exclude_fnames && (!(pcp->filter_flags & (kfFilterExcludeIbed0 | kfFilterExcludeIbed1))))) {
          // Must reconstruct the hash table in this case.
          BigstackReset(bigstack_mark);
          reterr = AllocAndPopulateIdHtableMt(variant_include, TO_CONSTCPCONSTP(variant_ids_mutable), variant_ct, pcp->max_thread_ct, &variant_id_htable, &htable_dup_base, &variant_id_htable_size);
          if (reterr) {
            goto Plink2Core_ret_1;
          }
        }
      }

      if (pcp->extract_fnames && (!(pcp->filter_flags & (kfFilterExtractIbed0 | kfFilterExtractIbed1)))) {
        reterr = ExtractExcludeFlagNorange(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, pcp->extract_fnames, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 0, variant_include, &variant_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->exclude_fnames && (!(pcp->filter_flags & (kfFilterExcludeIbed0 | kfFilterExcludeIbed1)))) {
        reterr = ExtractExcludeFlagNorange(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, pcp->exclude_fnames, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 1, variant_include, &variant_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      BigstackReset(bigstack_mark);
      if (full_variant_id_htable_needed) {
        reterr = ApplyVariantBpFilters(pcp->extract_fnames, pcp->exclude_fnames, cip, variant_bps, pcp->from_bp, pcp->to_bp, raw_variant_ct, pcp->filter_flags, vpos_sortstatus, variant_include, &variant_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      // todo: --update-alleles

      // todo: --attrib; although it isn't a "standard" file format, it can be
      // more convenient than forcing users to generate full-blown sites-only
      // VCF files, etc.
    }
    // variant_ids[] is fixed from this point on.
    const char* const* variant_ids = TO_CONSTCPCONSTP(variant_ids_mutable);
    // SetRefalt1FromFile() can alter pointers-to-missing in
    // allele_storage[] so we can't make this a const char* const*, but nothing
    // alters the pointed-to strings after this point.
    const char** allele_storage = K_CAST(const char**, allele_storage_mutable);

    if (pcp->thin_keep_prob != 1.0) {
      RandomThinProb("thin", "variant", pcp->thin_keep_prob, raw_variant_ct, variant_include, &variant_ct);
    } else if (pcp->thin_keep_ct != UINT32_MAX) {
      reterr = RandomThinCt("thin-count", "variant", pcp->thin_keep_ct, raw_variant_ct, variant_include, &variant_ct);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
    }
    uint32_t* sample_missing_dosage_cts = nullptr;
    uint32_t* sample_missing_hc_cts = nullptr;
    uint32_t* sample_hethap_cts = nullptr;
    uintptr_t max_covar_name_blen = 0;
    if (psamname[0]) {
      // xid_mode may vary between these operations in a single run, and
      // sample-sort is relatively cheap, so we abandon plink 1.9's "construct
      // sample ID map only once" optimization.
      if (pcp->update_sex_info.fname) {
        reterr = UpdateSampleSexes(sample_include, &pii.sii, &(pcp->update_sex_info), raw_sample_ct, sample_ct, sex_nm, sex_male);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->keepfam_fnames) {
        reterr = KeepOrRemove(pcp->keepfam_fnames, &pii.sii, raw_sample_ct, kfKeepFam, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->keep_fnames) {
        reterr = KeepOrRemove(pcp->keep_fnames, &pii.sii, raw_sample_ct, kfKeep0, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->removefam_fnames) {
        reterr = KeepOrRemove(pcp->removefam_fnames, &pii.sii, raw_sample_ct, kfKeepRemove | kfKeepFam, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->remove_fnames) {
        reterr = KeepOrRemove(pcp->remove_fnames, &pii.sii, raw_sample_ct, kfKeepRemove, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }

      // todo: --attrib-indiv

      if (pcp->keep_fcol_fname) {
        reterr = KeepFcol(pcp->keep_fcol_fname, &pii.sii, pcp->keep_fcol_flattened, pcp->keep_fcol_name, raw_sample_ct, pcp->keep_fcol_num, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->misc_flags & kfMiscRequirePheno) {
        reterr = RequirePheno(pheno_cols, pheno_names, pcp->require_pheno_flattened, raw_sample_ct, pheno_ct, max_pheno_name_blen, 0, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->filter_flags & (kfFilterExclFemales | kfFilterExclMales | kfFilterExclNosex)) {
        if (pcp->filter_flags & kfFilterExclFemales) {
          for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
            sample_include[widx] &= (~sex_nm[widx]) | sex_male[widx];
          }
        }
        if (pcp->filter_flags & kfFilterExclMales) {
          BitvecAndNot(sex_male, raw_sample_ctl, sample_include);
        }
        if (pcp->filter_flags & kfFilterExclNosex) {
          BitvecAnd(sex_nm, raw_sample_ctl, sample_include);
        }
        const uint32_t old_sample_ct = sample_ct;
        sample_ct = PopcountWords(sample_include, raw_sample_ctl);
        const uint32_t removed_ct = old_sample_ct - sample_ct;
        logprintf("%u sample%s removed due to sex filter(s).\n", removed_ct, (removed_ct == 1)? "" : "s");
      }
      if (pcp->filter_flags & (kfFilterExclFounders | kfFilterExclNonfounders)) {
        const uint32_t keep_founders = (pcp->filter_flags / kfFilterExclNonfounders) & 1;
        if (keep_founders) {
          BitvecAnd(founder_info, raw_sample_ctl, sample_include);
        } else {
          BitvecAndNot(founder_info, raw_sample_ctl, sample_include);
        }
        const uint32_t old_sample_ct = sample_ct;
        sample_ct = PopcountWords(sample_include, raw_sample_ctl);
        const uint32_t removed_ct = old_sample_ct - sample_ct;
        logprintf("--keep-%sfounders: %u sample%s removed.\n", keep_founders? "" : "non", removed_ct, (removed_ct == 1)? "" : "s");
      }

      if (pcp->thin_keep_sample_prob != 1.0) {
        RandomThinProb("thin-indiv", "sample", pcp->thin_keep_sample_prob, raw_sample_ct, sample_include, &sample_ct);
      } else if (pcp->thin_keep_sample_ct != UINT32_MAX) {
        reterr = RandomThinCt("thin-indiv-count", "sample", pcp->thin_keep_sample_ct, raw_sample_ct, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }

      const uint32_t smaj_missing_geno_report_requested = (pcp->command_flags1 & kfCommand1MissingReport) && (!(pcp->missing_rpt_flags & kfMissingRptVariantOnly));
      if ((pcp->mind_thresh < 1.0) || smaj_missing_geno_report_requested) {
        if (bigstack_alloc_u32(raw_sample_ct, &sample_missing_hc_cts) ||
            bigstack_alloc_u32(raw_sample_ct, &sample_hethap_cts)) {
          goto Plink2Core_ret_NOMEM;
        }
        if (SampleMissingDosageCtsAreNeeded(pcp->misc_flags, smaj_missing_geno_report_requested, pcp->mind_thresh, pcp->missing_rpt_flags)) {
          if (pgfi.gflags & kfPgenGlobalDosagePresent) {
            if (bigstack_alloc_u32(raw_sample_ct, &sample_missing_dosage_cts)) {
              goto Plink2Core_ret_NOMEM;
            }
          } else {
            sample_missing_dosage_cts = sample_missing_hc_cts;
          }
        }
        // could avoid this call and make LoadAlleleAndGenoCounts() do
        // double duty with --missing?
        reterr = LoadSampleMissingCts(sex_male, variant_include, cip, raw_variant_ct, variant_ct, raw_sample_ct, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, sample_missing_hc_cts, (pgfi.gflags & kfPgenGlobalDosagePresent)? sample_missing_dosage_cts : nullptr, sample_hethap_cts);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
        if (pcp->mind_thresh < 1.0) {
          uint32_t variant_ct_y = 0;
          uint32_t y_code;
          if (XymtExists(cip, kChrOffsetY, &y_code)) {
            variant_ct_y = CountChrVariantsUnsafe(variant_include, cip, y_code);
          }
          reterr = MindFilter((pcp->misc_flags & kfMiscMindDosage)? sample_missing_dosage_cts : sample_missing_hc_cts, (pcp->misc_flags & kfMiscMindHhMissing)? sample_hethap_cts : nullptr, &pii.sii, raw_sample_ct, variant_ct, variant_ct_y, pcp->mind_thresh, sample_include, sex_male, &sample_ct, outname, outname_end);
          if (reterr) {
            goto Plink2Core_ret_1;
          }
        }
        if (!smaj_missing_geno_report_requested) {
          BigstackReset(sample_missing_hc_cts);
        }
        // this results in a small "memory leak" when a regular missingness
        // report is requested, not a big deal
      }
      if (pcp->covar_fname || pcp->covar_range_list.name_ct) {
        const char* cur_covar_fname = pcp->covar_fname? pcp->covar_fname : (pcp->pheno_fname? pcp->pheno_fname : psamname);
        reterr = LoadPhenos(cur_covar_fname, &(pcp->covar_range_list), sample_include, pii.sii.sample_ids, raw_sample_ct, sample_ct, pii.sii.max_sample_id_blen, pcp->missing_pheno, 2, (pcp->misc_flags / kfMiscCovarColNums) & 1, &covar_cols, &covar_names, &covar_ct, &max_covar_name_blen);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
        logprintf("%u covariate%s loaded from %s.\n", covar_ct, (covar_ct == 1)? "" : "s", cur_covar_fname);

        // do we still want to clear some main phenotype values here if some
        // covariate values are missing?  (don't think there's a point to
        // preserving that behavior, let the regression functions do it to
        // their local phenotype copies on their own.)
      }

      if (pcp->misc_flags & kfMiscRequireCovar) {
        reterr = RequirePheno(covar_cols, covar_names, pcp->require_covar_flattened, raw_sample_ct, covar_ct, max_covar_name_blen, 1, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->keep_if_expr.pheno_name) {
        reterr = KeepRemoveIf(&(pcp->keep_if_expr), pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, (pcp->misc_flags / kfMiscAffection01) & 1, 0, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->remove_if_expr.pheno_name) {
        reterr = KeepRemoveIf(&(pcp->remove_if_expr), pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, (pcp->misc_flags / kfMiscAffection01) & 1, 1, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      // meow
      if (pcp->keep_cats_fname || pcp->keep_cat_names_flattened) {
        reterr = KeepRemoveCats(pcp->keep_cats_fname, pcp->keep_cat_names_flattened, pcp->keep_cat_phenoname, pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, 0, pcp->max_thread_ct, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->remove_cats_fname || pcp->remove_cat_names_flattened) {
        reterr = KeepRemoveCats(pcp->remove_cats_fname, pcp->remove_cat_names_flattened, pcp->remove_cat_phenoname, pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, 1, pcp->max_thread_ct, sample_include, &sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
    }

    const uint32_t nonfounders = (pcp->misc_flags / kfMiscNonfounders) & 1;
    uint32_t founder_ct = 0;
    uint32_t male_ct = 0;
    uint32_t nosex_ct = 0;
    if (psamname[0]) {
      if (!sample_ct) {
        logerrputs("Error: No samples remaining after main filters.\n");
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      UpdateSampleSubsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
      if (pcp->filter_flags & kfFilterPsamReq) {
        const uint32_t female_ct = sample_ct - male_ct - nosex_ct;
        if (!nosex_ct) {
          logprintfww("%u sample%s (%u female%s, %u male%s; %u founder%s) remaining after main filters.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", founder_ct, (founder_ct == 1)? "" : "s");
        } else {
          logprintfww("%u sample%s (%u female%s, %u male%s, %u ambiguous; %u founder%s) remaining after main filters.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", nosex_ct, founder_ct, (founder_ct == 1)? "" : "s");
        }
        if (pheno_ct == 1) {
          const PhenoDtype pheno_type_code = pheno_cols[0].type_code;
          const uint32_t obs_ct = PopcountWordsIntersect(pheno_cols[0].nonmiss, sample_include, raw_sample_ctl);
          if (pheno_type_code == kPhenoDtypeCc) {
            const uint32_t case_ct = PopcountWordsIntersect(pheno_cols[0].data.cc, sample_include, raw_sample_ctl);
            const uint32_t ctrl_ct = obs_ct - case_ct;
            logprintf("%u case%s and %u control%s remaining after main filters.\n", case_ct, (case_ct == 1)? "" : "s", ctrl_ct, (ctrl_ct == 1)? "" : "s");
          } else {
            logprintf("%u %s phenotype value%s remaining after main filters.\n", obs_ct, (pheno_type_code == kPhenoDtypeQt)? "quantitative" : "categorical", (obs_ct == 1)? "" : "s");
          }
        }
      }
    }
    if (pcp->pheno_transform_flags & kfPhenoTransformSplitCat) {
      reterr = SplitCatPheno(pcp->split_cat_phenonames_flattened, sample_include, raw_sample_ct, pcp->pheno_transform_flags, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen, &covar_cols, &covar_names, &covar_ct, &max_covar_name_blen);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
    }

    // quantnorm before variance-standardize, since at least that has a minor
    // effect, whereas the other order is pointless
    if (pcp->pheno_transform_flags & (kfPhenoTransformQuantnormPheno | kfPhenoTransformQuantnormAll)) {
      reterr = PhenoQuantileNormalize(pcp->quantnorm_flattened, sample_include, pheno_names, raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, 0, (pcp->pheno_transform_flags / kfPhenoTransformQuantnormPheno) & 1, pheno_cols);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
    }
    if (pcp->pheno_transform_flags & (kfPhenoTransformQuantnormCovar | kfPhenoTransformQuantnormAll)) {
      reterr = PhenoQuantileNormalize((pcp->pheno_transform_flags & kfPhenoTransformQuantnormAll)? pcp->quantnorm_flattened : pcp->covar_quantnorm_flattened, sample_include, covar_names, raw_sample_ct, sample_ct, covar_ct, max_covar_name_blen, 1, (pcp->pheno_transform_flags / kfPhenoTransformQuantnormCovar) & 1, covar_cols);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
    }

    if (pcp->pheno_transform_flags & (kfPhenoTransformVstdCovar | kfPhenoTransformVstdAll)) {
      const uint32_t is_covar_flag = (pcp->pheno_transform_flags / kfPhenoTransformVstdCovar) & 1;
      if (!is_covar_flag) {
        reterr = PhenoVarianceStandardize(pcp->vstd_flattened, sample_include, pheno_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, 0, 0, pheno_cols);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      reterr = PhenoVarianceStandardize(pcp->vstd_flattened, sample_include, covar_names, raw_sample_ct, covar_ct, max_covar_name_blen, 1, is_covar_flag, covar_cols);
      if (reterr) {
        goto Plink2Core_ret_1;
      }
    }

    char* loop_cats_outname_endp1_backup = &(outname_end[1]);
    uintptr_t* loop_cats_sample_include_backup = nullptr;
    uintptr_t* loop_cats_founder_info_backup = nullptr;
    uintptr_t* loop_cats_sex_nm_backup = nullptr;
    uintptr_t* loop_cats_sex_male_backup = nullptr;
    uintptr_t* loop_cats_variant_include_backup = nullptr;
    uintptr_t* loop_cats_cat_include = nullptr;
    uint32_t loop_cats_sample_ct = 0;
    uint32_t loop_cats_variant_ct = variant_ct;
    uint32_t loop_cats_uidx = 0;
    uint32_t loop_cats_idx = 0;
    uint32_t loop_cats_ct = 1;
    if (pcp->loop_cats_phenoname) {
      // 1. check phenotype names, verify it's a categorical pheno; fall back
      //    on checking covariates if no phenotype has same name
      // 2. remove phenotype/covariate, back up sample_include
      // 3. iterate through positive integer category numbers, skip empty
      //    categories
      const char* loop_cats_phenoname = pcp->loop_cats_phenoname;
      const uintptr_t name_blen = 1 + strlen(loop_cats_phenoname);
      if (name_blen <= max_pheno_name_blen) {
        // this boilerplate may belong in its own function
        for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
          if (!memcmp(loop_cats_phenoname, &(pheno_names[pheno_idx * max_pheno_name_blen]), name_blen)) {
            PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
            if (cur_pheno_col->type_code != kPhenoDtypeCat) {
              logerrprintfww("Error: '%s' is not a categorical phenotype.\n", loop_cats_phenoname);
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            loop_cats_pheno_col = cur_pheno_col;
            --pheno_ct;
            for (uint32_t pheno_idx2 = pheno_idx; pheno_idx2 < pheno_ct; ++pheno_idx2) {
              memcpy(&(pheno_cols[pheno_idx2]), &(pheno_cols[pheno_idx2 + 1]), sizeof(PhenoCol));
              memcpy(&(pheno_names[pheno_idx2 * max_pheno_name_blen]), &(pheno_names[(pheno_idx2 + 1) * max_pheno_name_blen]), max_pheno_name_blen);
            }
            break;
          }
        }
      }
      if ((!loop_cats_pheno_col) && (name_blen <= max_covar_name_blen)) {
        for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx) {
          if (!memcmp(loop_cats_phenoname, &(covar_names[covar_idx * max_covar_name_blen]), name_blen)) {
            PhenoCol* cur_covar_col = &(covar_cols[covar_idx]);
            if (cur_covar_col->type_code != kPhenoDtypeCat) {
              logerrprintfww("Error: '%s' is not a categorical covariate.\n", loop_cats_phenoname);
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            loop_cats_pheno_col = cur_covar_col;
            --covar_ct;
            for (uint32_t covar_idx2 = covar_idx; covar_idx2 < covar_ct; ++covar_idx2) {
              memcpy(&(covar_cols[covar_idx2]), &(covar_cols[covar_idx2 + 1]), sizeof(PhenoCol));
              memcpy(&(covar_names[covar_idx2 * max_covar_name_blen]), &(covar_names[(covar_idx2 + 1) * max_covar_name_blen]), max_covar_name_blen);
            }
            break;
          }
        }
      }
      if (!loop_cats_pheno_col) {
        logerrprintfww("Error: --loop-cats phenotype '%s' not loaded.\n", loop_cats_phenoname);
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      if (bigstack_alloc_w(raw_sample_ctl, &loop_cats_sample_include_backup) ||
          bigstack_alloc_w(raw_sample_ctl, &loop_cats_founder_info_backup) ||
          bigstack_alloc_w(raw_sample_ctl, &loop_cats_sex_nm_backup) ||
          bigstack_alloc_w(raw_sample_ctl, &loop_cats_sex_male_backup) ||
          bigstack_alloc_w(1 + (loop_cats_pheno_col->nonnull_category_ct / kBitsPerWord), &loop_cats_cat_include)) {
        goto Plink2Core_ret_NOMEM;
      }
      if (variant_ct != raw_variant_ct) {
        if (bigstack_alloc_w(raw_variant_ctl, &loop_cats_variant_include_backup)) {
          goto Plink2Core_ret_NOMEM;
        }
        memcpy(loop_cats_variant_include_backup, variant_include, raw_variant_ctl * sizeof(intptr_t));
      }
      BitvecAndCopy(sample_include, loop_cats_pheno_col->nonmiss, raw_sample_ctl, loop_cats_sample_include_backup);
      loop_cats_sample_ct = PopcountWords(loop_cats_sample_include_backup, raw_sample_ctl);
      loop_cats_ct = IdentifyRemainingCats(sample_include, loop_cats_pheno_col, sample_ct, loop_cats_cat_include);
      if (!loop_cats_ct) {
        logerrputs("Error: All --loop-cats categories are empty.\n");
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      logprintf("--loop-cats: %u categor%s present.\n", loop_cats_ct, (loop_cats_ct == 1)? "y" : "ies");
      BitvecAndCopy(sample_include, founder_info, raw_sample_ctl, loop_cats_founder_info_backup);
      BitvecAndCopy(sample_include, sex_nm, raw_sample_ctl, loop_cats_sex_nm_backup);
      BitvecAndCopy(sample_include, sex_male, raw_sample_ctl, loop_cats_sex_male_backup);
      *outname_end = '.';
    }
    const uintptr_t raw_allele_ct = variant_allele_idxs? variant_allele_idxs[raw_variant_ct] : (2 * raw_variant_ct);
    unsigned char* bigstack_mark_varfilter = g_bigstack_base;
    unsigned char* bigstack_end_mark_varfilter = g_bigstack_end;
    while (1) {
      if (loop_cats_pheno_col) {
        ++loop_cats_uidx;
        MovU32To1Bit(loop_cats_cat_include, &loop_cats_uidx);
        const char* catname = loop_cats_pheno_col->category_names[loop_cats_uidx];
        const uint32_t catname_slen = strlen(catname);
        if (catname_slen + S_CAST(uintptr_t, loop_cats_outname_endp1_backup - outname) > (kPglFnamesize - kMaxOutfnameExtBlen)) {
          logerrputs("Error: --loop-cats category name too long.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        outname_end = memcpya(loop_cats_outname_endp1_backup, catname, catname_slen);
        sample_ct = GetCatSamples(loop_cats_sample_include_backup, loop_cats_pheno_col, raw_sample_ctl, loop_cats_sample_ct, loop_cats_uidx, sample_include);
        memcpy(founder_info, loop_cats_founder_info_backup, raw_sample_ctl * sizeof(intptr_t));
        memcpy(sex_nm, loop_cats_sex_nm_backup, raw_sample_ctl * sizeof(intptr_t));
        memcpy(sex_male, loop_cats_sex_male_backup, raw_sample_ctl * sizeof(intptr_t));
        UpdateSampleSubsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
        variant_ct = loop_cats_variant_ct;
        if (loop_cats_variant_include_backup) {
          memcpy(variant_include, loop_cats_variant_include_backup, raw_variant_ctl * sizeof(intptr_t));
        } else {
          SetAllBits(variant_ct, variant_include);
        }
        logprintf("--loop-cats: Processing category '%s' (%u sample%s).\n", catname, sample_ct, (sample_ct == 1)? "" : "s");
      }

      // dosages are currently in 32768ths
      uint64_t* allele_dosages = nullptr;  // same indexes as allele_storage
      uint64_t* founder_allele_dosages = nullptr;
      AltAlleleCt* maj_alleles = nullptr;
      double* allele_freqs = nullptr;
      uint32_t* raw_geno_cts = nullptr;
      uint32_t* founder_raw_geno_cts = nullptr;
      unsigned char* bigstack_mark_allele_dosages = g_bigstack_base;
      unsigned char* bigstack_mark_founder_allele_dosages = g_bigstack_base;
      const uint32_t keep_grm = GrmKeepIsNeeded(pcp->command_flags1, pcp->pca_flags);
      double* grm = nullptr;

      if (pcp->command_flags1 & kfCommand1WriteSamples) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".id");
        reterr = WriteSampleIds(sample_include, &pii.sii, outname, sample_ct);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
        logprintfww("--write-samples: Sample IDs written to %s .\n", outname);
        if (!(pcp->command_flags1 & (~kfCommand1WriteSamples))) {
          goto Plink2Core_early_complete;
        }
      }

      if (pgenname[0]) {
        if (AlleleFreqsAreNeeded(pcp->command_flags1, pcp->glm_info.flags, pcp->min_maf, pcp->max_maf)) {
          if (MajAllelesAreNeeded(pcp->command_flags1, pcp->glm_info.flags)) {
            maj_alleles = S_CAST(AltAlleleCt*, bigstack_alloc(raw_variant_ct * sizeof(AltAlleleCt)));
            if (!maj_alleles) {
              goto Plink2Core_ret_NOMEM;
            }
          }
          //   allele_freqs[variant_allele_idxs[variant_uidx] - variant_uidx]
          // stores the frequency estimate for the reference allele; if there's
          // more than 1 alt allele, next element stores alt1 freq, etc.  To
          // save memory, we omit the last alt.
          if (bigstack_alloc_d(raw_allele_ct - raw_variant_ct, &allele_freqs)) {
            goto Plink2Core_ret_NOMEM;
          }
        }
        uint32_t x_start = 0;
        uint32_t x_len = 0;
        uint32_t hwe_x_probs_needed = 0;
        uint32_t x_code;
        if ((!(vpos_sortstatus & kfUnsortedVarSplitChr)) && XymtExists(cip, kChrOffsetX, &x_code)) {
          const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
          x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
          const uint32_t x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
          x_len = x_end - x_start;
          if (x_len && ((pcp->command_flags1 & kfCommand1Hardy) || (pcp->hwe_thresh != 1.0)) && (!AllBitsAreZero(variant_include, x_start, x_end))) {
            if (nonfounders) {
              hwe_x_probs_needed = (sample_ct > nosex_ct);
            } else {
              for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
                if (founder_info[widx] & sex_nm[widx]) {
                  // at least one founder with known gender
                  hwe_x_probs_needed = 1;
                  break;
                }
              }
            }
          }
        }
        bigstack_mark_allele_dosages = g_bigstack_base;
        const uint32_t first_hap_uidx = GetFirstHaploidUidx(cip, vpos_sortstatus);
        if (AlleleDosagesAreNeeded(pcp->misc_flags, make_plink2_flags, (allele_freqs != nullptr), pcp->min_allele_dosage, pcp->max_allele_dosage)) {
          if (bigstack_alloc_u64(raw_allele_ct, &allele_dosages)) {
            goto Plink2Core_ret_NOMEM;
          }
        }
        bigstack_mark_founder_allele_dosages = g_bigstack_base;
        if (FounderAlleleDosagesAreNeeded(pcp->misc_flags, (allele_freqs != nullptr), pcp->min_allele_dosage, pcp->max_allele_dosage)) {
          if ((founder_ct == sample_ct) && allele_dosages) {
            founder_allele_dosages = allele_dosages;
          } else {
            if (bigstack_alloc_u64(raw_allele_ct, &founder_allele_dosages)) {
              goto Plink2Core_ret_NOMEM;
            }
          }
        }
        double* mach_r2_vals = nullptr;
        if ((pcp->freq_rpt_flags & kfAlleleFreqColMachR2) || (pcp->mach_r2_max != 0.0)) {
          if (bigstack_alloc_d(raw_variant_ct, &mach_r2_vals)) {
            goto Plink2Core_ret_NOMEM;
          }
        }

        unsigned char* bigstack_mark_geno_cts = g_bigstack_base;

        // no longer includes hethaps by default
        uint32_t* variant_missing_hc_cts = nullptr;
        uint32_t* variant_hethap_cts = nullptr;
        if (VariantMissingHcCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->geno_thresh, pcp->missing_rpt_flags)) {
          if (bigstack_alloc_u32(raw_variant_ct, &variant_missing_hc_cts)) {
            goto Plink2Core_ret_NOMEM;
          }
          if (VariantHethapCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->geno_thresh, pcp->missing_rpt_flags, first_hap_uidx)) {
            // first_hap_uidx offset can save an entire GB...
            if (bigstack_alloc_u32(raw_variant_ct - first_hap_uidx, &variant_hethap_cts)) {
              goto Plink2Core_ret_NOMEM;
            }
          }
        }
        uint32_t* variant_missing_dosage_cts = nullptr;
        if (VariantMissingDosageCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->geno_thresh, pcp->missing_rpt_flags)) {
          if ((!variant_missing_hc_cts) || (pgfi.gflags & kfPgenGlobalDosagePresent)) {
            if (bigstack_alloc_u32(raw_variant_ct, &variant_missing_dosage_cts)) {
              goto Plink2Core_ret_NOMEM;
            }
          } else {
            variant_missing_dosage_cts = variant_missing_hc_cts;
          }
        }
        uint32_t* x_male_geno_cts = nullptr;
        uint32_t* founder_x_male_geno_cts = nullptr;
        uint32_t* x_nosex_geno_cts = nullptr;
        uint32_t* founder_x_nosex_geno_cts = nullptr;
        // [3n] = homref ct, [3n+1] = het ref-altx total, [3n+2] = nonref
        //   diploid total
        // use unfiltered indexes, since we remove more variants later
        if (RawGenoCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->hwe_thresh)) {
          if (bigstack_alloc_u32((3 * k1LU) * raw_variant_ct, &raw_geno_cts)) {
            goto Plink2Core_ret_NOMEM;
          }
          if (x_len) {
            if (male_ct) {
              if (bigstack_alloc_u32((3 * k1LU) * x_len, &x_male_geno_cts)) {
                goto Plink2Core_ret_NOMEM;
              }
            }
            if (nosex_ct && hwe_x_probs_needed && nonfounders) {
              if (bigstack_alloc_u32((3 * k1LU) * x_len, &x_nosex_geno_cts)) {
                goto Plink2Core_ret_NOMEM;
              }
            }
          }
        }
        if (FounderRawGenoCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->hwe_thresh)) {
          if ((founder_ct == sample_ct) && raw_geno_cts) {
            founder_raw_geno_cts = raw_geno_cts;
            founder_x_male_geno_cts = x_male_geno_cts;
          } else {
            if (bigstack_alloc_u32((3 * k1LU) * raw_variant_ct, &founder_raw_geno_cts)) {
              goto Plink2Core_ret_NOMEM;
            }
            if (x_len && male_ct) {
              const uint32_t founder_male_ct = PopcountWordsIntersect(founder_info, sex_male, raw_sample_ctl);
              if (founder_male_ct) {
                if (bigstack_alloc_u32((3 * k1LU) * x_len, &founder_x_male_geno_cts)) {
                  goto Plink2Core_ret_NOMEM;
                }
              }
            }
          }
          if (nosex_ct && hwe_x_probs_needed && (!nonfounders)) {
            const uint32_t founder_knownsex_ct = PopcountWordsIntersect(founder_info, sex_nm, raw_sample_ctl);
            if (founder_knownsex_ct < founder_ct) {
              if ((founder_ct == sample_ct) && x_nosex_geno_cts) {
                // shouldn't be possible for now
                assert(0);
                // founder_x_nosex_geno_cts = x_nosex_geno_cts;
              } else {
                if (bigstack_alloc_u32((3 * k1LU) * x_len, &founder_x_nosex_geno_cts)) {
                  goto Plink2Core_ret_NOMEM;
                }
              }
            }
          }
        }
        if (allele_dosages || founder_allele_dosages || variant_missing_hc_cts || variant_missing_dosage_cts || variant_hethap_cts || raw_geno_cts || founder_raw_geno_cts || mach_r2_vals) {
          // note that --geno depends on different handling of X/Y than --maf.

          // possible todo: "free" these arrays early in some cases
          // todo: oblig-missing

          // bugfix (22 Nov 2017): given "--genotyping-rate dosage" on a file
          // with no dosage data, variant_missing_hc_cts may be nullptr while
          // variant_missing_dosage_cts is zero-initialized.  In this case,
          // it's fine to pass variant_missing_dosage_cts in the
          // hardcall-missing-count slot... and it's NOT fine to pass in
          // nullptrs for both missing-count arrays...
          const uint32_t dosageless_file = !(pgfi.gflags & kfPgenGlobalDosagePresent);
          reterr = LoadAlleleAndGenoCounts(sample_include, founder_info, sex_nm, sex_male, variant_include, cip, variant_allele_idxs, raw_sample_ct, sample_ct, founder_ct, male_ct, nosex_ct, raw_variant_ct, variant_ct, first_hap_uidx, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, allele_dosages, founder_allele_dosages, ((!variant_missing_hc_cts) && dosageless_file)? variant_missing_dosage_cts : variant_missing_hc_cts, dosageless_file? nullptr : variant_missing_dosage_cts, variant_hethap_cts, raw_geno_cts, founder_raw_geno_cts, x_male_geno_cts, founder_x_male_geno_cts, x_nosex_geno_cts, founder_x_nosex_geno_cts, mach_r2_vals);
          if (reterr) {
            goto Plink2Core_ret_1;
          }
          if (pcp->command_flags1 & kfCommand1GenotypingRate) {
            const uint32_t is_dosage = (pcp->misc_flags / kfMiscGenotypingRateDosage) & 1;
            ReportGenotypingRate(variant_include, cip, is_dosage? variant_missing_dosage_cts : variant_missing_hc_cts, raw_sample_ct, sample_ct, male_ct, variant_ct, is_dosage);
            if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1WriteSamples)))) {
              goto Plink2Core_early_complete;
            }
          }
        }
        if (allele_freqs) {
          const uint32_t maf_succ = (pcp->misc_flags / kfMiscMafSucc) & 1;
          ComputeAlleleFreqs(variant_include, variant_allele_idxs, nonfounders? allele_dosages : founder_allele_dosages, variant_ct, maf_succ, allele_freqs);
          if (pcp->read_freq_fname) {
            reterr = ReadAlleleFreqs(variant_include, variant_ids, variant_allele_idxs, allele_storage, pcp->read_freq_fname, raw_variant_ct, variant_ct, pgfi.max_alt_allele_ct, max_variant_id_slen, max_allele_slen, maf_succ, pcp->max_thread_ct, allele_freqs);
            if (reterr) {
              goto Plink2Core_ret_1;
            }
          }
          if (maj_alleles) {
            ComputeMajAlleles(variant_include, variant_allele_idxs, allele_freqs, variant_ct, maj_alleles);
          }
        } else if (pcp->read_freq_fname) {
          logerrprintf("Warning: Ignoring --read-freq since no command would use the frequencies.\n");
        }

        if (pcp->command_flags1 & kfCommand1AlleleFreq) {
          reterr = WriteAlleleFreqs(variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, nonfounders? allele_dosages : founder_allele_dosages, mach_r2_vals, pcp->freq_ref_binstr, pcp->freq_alt1_binstr, variant_ct, pgfi.max_alt_allele_ct, max_allele_slen, pcp->freq_rpt_flags, pcp->max_thread_ct, nonfounders, outname, outname_end);
          if (reterr) {
            goto Plink2Core_ret_1;
          }
          if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1AlleleFreq | kfCommand1WriteSamples)))) {
            goto Plink2Core_early_complete;
          }
        }
        if (pcp->command_flags1 & kfCommand1GenoCounts) {
          reterr = WriteGenoCounts(sample_include, sex_male, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, raw_geno_cts, x_male_geno_cts, raw_sample_ct, sample_ct, male_ct, variant_ct, x_start, max_allele_slen, pcp->geno_counts_flags, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
          if (reterr) {
            goto Plink2Core_ret_1;
          }
          if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1AlleleFreq | kfCommand1GenoCounts | kfCommand1WriteSamples)))) {
            goto Plink2Core_early_complete;
          }
        }

        if (pcp->command_flags1 & kfCommand1MissingReport) {
          reterr = WriteMissingnessReports(sample_include, &pii.sii, sex_male, pheno_cols, pheno_names, sample_missing_hc_cts, sample_missing_dosage_cts, sample_hethap_cts, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, variant_missing_hc_cts, variant_missing_dosage_cts, variant_hethap_cts, sample_ct, male_ct, pheno_ct, max_pheno_name_blen, variant_ct, max_allele_slen, variant_hethap_cts? first_hap_uidx : 0x7fffffff, pcp->missing_rpt_flags, pcp->max_thread_ct, outname, outname_end);
          if (reterr) {
            goto Plink2Core_ret_1;
          }
          if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1AlleleFreq | kfCommand1GenoCounts | kfCommand1MissingReport | kfCommand1WriteSamples)))) {
            goto Plink2Core_early_complete;
          }
        }

        if (pcp->geno_thresh != 1.0) {
          const uint32_t geno_hh_missing = S_CAST(uint32_t, pcp->misc_flags & kfMiscGenoHhMissing);
          EnforceGenoThresh(cip, (pcp->misc_flags & kfMiscGenoDosage)? variant_missing_dosage_cts : variant_missing_hc_cts, geno_hh_missing? variant_hethap_cts : nullptr, sample_ct, male_ct, geno_hh_missing? first_hap_uidx : 0x7fffffff, pcp->geno_thresh, variant_include, &variant_ct);
        }

        double* hwe_x_pvals = nullptr;
        uint32_t hwe_x_ct = 0;
        if (hwe_x_probs_needed) {
          hwe_x_ct = CountChrVariantsUnsafe(variant_include, cip, cip->xymt_codes[kChrOffsetX]);
          // hwe_x_ct == 0 possible, if --geno filters out all remaining chrX
          // variants
          // also support suppression of --hardy p column (with a gigantic
          // dataset, maybe it's reasonable to stick to femalep, etc.)
          if (hwe_x_ct && ((pcp->hwe_thresh != 1.0) || (pcp->hardy_flags & kfHardyColP))) {
            uint32_t hwe_midp;
            if (pcp->command_flags1 & kfCommand1Hardy) {
              hwe_midp = (pcp->hardy_flags / kfHardyMidp) & 1;
              if (pcp->hwe_thresh != 1.0) {
                const uint32_t hwe_midp2 = (pcp->misc_flags / kfMiscHweMidp) & 1;
                if (hwe_midp != hwe_midp2) {
                  // could support this efficiently, but why bother...
                  logerrputs("Error: --hardy and --hwe must have identical midp settings when chrX is\npresent.\n");
                  goto Plink2Core_ret_INVALID_CMDLINE;
                }
              }
            } else {
              hwe_midp = (pcp->misc_flags / kfMiscHweMidp) & 1;
            }
            reterr = ComputeHweXPvals(variant_include, nonfounders? raw_geno_cts : founder_raw_geno_cts, nonfounders? x_male_geno_cts : founder_x_male_geno_cts, nonfounders? x_nosex_geno_cts : founder_x_nosex_geno_cts, x_start, hwe_x_ct, hwe_midp, pcp->max_thread_ct, &hwe_x_pvals);
            if (reterr) {
              goto Plink2Core_ret_1;
            }
          }
        }
        if (pcp->command_flags1 & kfCommand1Hardy) {
          reterr = HardyReport(variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, nonfounders? raw_geno_cts : founder_raw_geno_cts, nonfounders? x_male_geno_cts : founder_x_male_geno_cts, nonfounders? x_nosex_geno_cts : founder_x_nosex_geno_cts, hwe_x_pvals, variant_ct, hwe_x_ct, max_allele_slen, pcp->output_min_p, pcp->hardy_flags, pcp->max_thread_ct, nonfounders, outname, outname_end);
          if (reterr) {
            goto Plink2Core_ret_1;
          }
          if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1AlleleFreq | kfCommand1GenoCounts | kfCommand1MissingReport | kfCommand1Hardy | kfCommand1WriteSamples)))) {
            goto Plink2Core_early_complete;
          }
        }
        if (pcp->hwe_thresh != 1.0) {
          // assumes no filtering between hwe_x_pvals[] computation and here
          EnforceHweThresh(cip, nonfounders? raw_geno_cts : founder_raw_geno_cts, nonfounders? x_male_geno_cts : founder_x_male_geno_cts, nonfounders? x_nosex_geno_cts : founder_x_nosex_geno_cts, hwe_x_pvals, pcp->misc_flags, pcp->hwe_thresh, nonfounders, variant_include, &variant_ct);
        }
        // raw_geno_cts/founder_raw_geno_cts/hwe_x_pvals no longer needed
        BigstackReset(bigstack_mark_geno_cts);

        if ((pcp->min_maf != 0.0) || (pcp->max_maf != 1.0) || pcp->min_allele_dosage || (pcp->max_allele_dosage != (~0LLU))) {
          EnforceMinorFreqConstraints(variant_allele_idxs, nonfounders? allele_dosages : founder_allele_dosages, allele_freqs, pcp->min_maf, pcp->max_maf, pcp->min_allele_dosage, pcp->max_allele_dosage, variant_include, &variant_ct);
        }

        if (mach_r2_vals) {
          if (pcp->mach_r2_max != 0.0) {
            EnforceMachR2Thresh(cip, mach_r2_vals, pcp->mach_r2_min, pcp->mach_r2_max, variant_include, &variant_ct);
          }
          BigstackReset(mach_r2_vals);
        }
      }

      if (pcp->min_bp_space) {
        if (vpos_sortstatus & kfUnsortedVarBp) {
          logerrputs("Error: --bp-space requires a sorted .pvar/.bim.  Retry this command after using\n--make-pgen/--make-bed + --sort-vars to sort your data.\n");
        }
        EnforceMinBpSpace(cip, variant_bps, pcp->min_bp_space, variant_include, &variant_ct);
      }

      if (pcp->filter_flags & kfFilterPvarReq) {
        if (!variant_ct) {
          logerrputs("Error: No variants remaining after main filters.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        logprintf("%u variant%s remaining after main filters.\n", variant_ct, (variant_ct == 1)? "" : "s");
      }

      if (pgenname[0]) {
        if (pcp->command_flags1 & (kfCommand1MakeKing | kfCommand1KingCutoff)) {
          uintptr_t* prev_sample_include = nullptr;
          const uint32_t prev_sample_ct = sample_ct;
          if (pcp->king_cutoff != -1) {
            if (bigstack_alloc_w(raw_sample_ctl, &prev_sample_include)) {
              goto Plink2Core_ret_NOMEM;
            }
            memcpy(prev_sample_include, sample_include, raw_sample_ctl * sizeof(intptr_t));
          }
          if (pcp->king_table_subset_fname) {
            // command-line parser currently guarantees --king-table-subset
            // isn't used with --king-cutoff or --make-king
            // probable todo: --king-cutoff-table which can use .kin0 as input
            reterr = CalcKingTableSubset(sample_include, &pii.sii, variant_include, cip, pcp->king_table_subset_fname, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, pcp->king_table_filter, pcp->king_table_subset_thresh, pcp->king_flags, pcp->parallel_idx, pcp->parallel_tot, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
            if (reterr) {
              goto Plink2Core_ret_1;
            }
          } else {
            if (king_cutoff_fprefix) {
              reterr = KingCutoffBatch(&pii.sii, raw_sample_ct, pcp->king_cutoff, sample_include, king_cutoff_fprefix, &sample_ct);
            } else {
              reterr = CalcKing(&pii.sii, variant_include, cip, raw_sample_ct, raw_variant_ct, variant_ct, pcp->king_cutoff, pcp->king_table_filter, pcp->king_flags, pcp->parallel_idx, pcp->parallel_tot, pcp->max_thread_ct, &simple_pgr, sample_include, &sample_ct, outname, outname_end);
            }
            if (reterr) {
              goto Plink2Core_ret_1;
            }
            if (pcp->king_cutoff != -1) {
              snprintf(outname_end, kMaxOutfnameExtBlen, ".king.cutoff.in.id");
              reterr = WriteSampleIds(sample_include, &pii.sii, outname, sample_ct);
              if (reterr) {
                goto Plink2Core_ret_1;
              }
              snprintf(&(outname_end[13]), kMaxOutfnameExtBlen - 13, "out.id");
              BitvecAndNot(sample_include, raw_sample_ctl, prev_sample_include);
              const uint32_t removed_sample_ct = prev_sample_ct - sample_ct;
              reterr = WriteSampleIds(prev_sample_include, &pii.sii, outname, removed_sample_ct);
              if (reterr) {
                goto Plink2Core_ret_1;
              }
              BigstackReset(prev_sample_include);
              outname_end[13] = '\0';
              logprintfww("--king-cutoff: Excluded sample ID%s written to %sout, and %u remaining sample ID%s written to %sin .\n", (removed_sample_ct == 1)? "" : "s", outname, sample_ct, (sample_ct == 1)? "" : "s", outname);
              UpdateSampleSubsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
            }
          }
        }
      }
      if ((pcp->command_flags1 & kfCommand1MakeRel) || keep_grm) {
        reterr = CalcGrm(sample_include, &pii.sii, variant_include, cip, variant_allele_idxs, maj_alleles, allele_freqs, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, pcp->grm_flags, pcp->parallel_idx, pcp->parallel_tot, pcp->max_thread_ct, &simple_pgr, outname, outname_end, keep_grm? (&grm) : nullptr);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
        // Retire --rel-cutoff, since --king-cutoff is pretty clearly better.
        // KING-robust has significant systematic biases when interracial
        // couples are involved, though.  Still may be okay for first-degree
        // in that context, but something like PC-Relate should be added soon.
        //
        // In addition to even better kinship coefficients, PC-Relate also
        // provides a replacement for --genome's obsolete IBD sharing
        // estimates, and improved inbreeding coefficients.  There may still be
        // work to do on handling of highly-inbred populations, but PC-Relate
        // appears to be enough of an advance to warrant including now without
        // waiting for further refinements.
        //
        // Assuming we follow through with PC-Relate:
        //   1. use basic PC-AiR method to estimate top PCs in a way that's
        //      more resistant to related samples:
        //      a. compute KING-robust matrix.  keep in memory for now; maybe
        //         allow this to spill to temporary file later since this is
        //         a major scaling limitation.  probably store as floats,
        //         because of scaling; even then, 100k samples requires over
        //         280 GB, and it goes up quadratically from there.
        //         (though the problem isn't actually very serious, since your
        //         PCs will practically always be good enough if you select a
        //         size-10k random sample)
        //      b. also compute ancestry divergence matrix, with entries
        //           0.5 * (1 - \frac{\sum (g_1 - g_2)^2}{hets_1 + hets_2})
        //         main CalcKing() loop doesn't need to be changed; numerator
        //         of fraction is 4 * ibs0_ct + het1hom2_ct + het2hom1_ct,
        //         denominator is 2 * hethet_ct + het1hom2_ct + het2hom1_ct
        //      c. use the algorithm in Appendix B of
        //           https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4836868/
        //         to select an unrelated subset.  (unfortunately, step 7 is
        //         random, so we can't always test for perfect concordance with
        //         the R package.)
        //      d. compute PCs on the unrelated subset, project to related
        //         subset in the usual way.  (since overall algorithm has a
        //         random component anyway, there's no longer much of a point
        //         to defaulting to nonrandom PCA for >5k samples; instead,
        //         have 'random-pca'/'nonrandom-pca' modifiers to force, and
        //         default to switching over at 5k)
        //         allow these PCs to be saved to disk, and later steps to
        //         load them as input.
        //         extra MAF cutoff should be applied here.
        //         may want default PC count to be 6-8 rather than 10.
        //   2. for kinship matrix:
        //      a. for each variant,
        //         i. estimate individual-specific allele frequencies, by
        //            performing linear regression with y = genotype vector and
        //            X = constant column + top PCs from step 1; clip to
        //            [epsilon, 1 - epsilon]
        //            (this is also a step in the computations to follow, but
        //            we don't default to saving the intermediate result to
        //            disk because it's huge, and quick to compute from PCs +
        //            genotype matrix)
        //         ii. compute vector of (g-2q)/sqrt(q(1-q)) values
        //      b. same incremental matrix multiplies as GRM computation,
        //         dsyrk() is our friend
        //   3. for inbreeding coefficients, see equation 6 of
        //        https://www.ncbi.nlm.nih.gov/pubmed/26748516
        //   4. for IBD sharing, see equations 9 and 10 (note that inbreeding
        //      coefficient enters in to 9, and kinship enters into 10).
        //      when classifying first-degree relationships (probably want to
        //      add an option to do this; can consider second-degree too if it
        //      passes accuracy tests), parent-offspring corresponds to
        //      P(IBD=0) < 2^{-9/2} (todo: check theoretical justification for
        //      this threshold, if it's arbitrary it should be tunable)
        //   5. we may want to extend other allele-frequency-dependent
        //      commonly-used functions (--check-sex/--impute-sex is the most
        //      important one that comes to mind) to be able to use
        //      individual-specific allele frequencies.

        // possible todo: unrelated heritability?
      }
#ifndef NOLAPACK
      if (pcp->command_flags1 & kfCommand1Pca) {
        // if the GRM is on the stack, this always frees it
        reterr = CalcPca(sample_include, &pii.sii, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, maj_alleles, allele_freqs, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, pcp->pca_ct, pcp->pca_flags, pcp->max_thread_ct, &simple_pgr, grm, outname, outname_end);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
#endif

      if (pcp->command_flags1 & kfCommand1WriteSnplist) {
        reterr = WriteSnplist(variant_include, variant_ids, variant_ct, (pcp->misc_flags / kfMiscWriteSnplistZs) & 1, pcp->max_thread_ct, outname, outname_end);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & (kfCommand1MakePlink2 | kfCommand1Exportf | kfCommand1WriteCovar)) {
        // If non-null, refalt1_select has (2 * raw_variant_ct) entries.  [2n]
        // stores the index of the new ref allele for variant n, and [2n+1]
        // stores the index of the new alt1 allele.  (0 = original ref, 1 =
        // original alt1, etc.)
        // Operations which instantiate this (--maj-ref, --ref-allele,
        // --alt1-allele) are only usable with fileset creation commands.  No
        // more pass-marker_reverse-to-everything nonsense.
        // (Technically, I could also drop support for --export, but that would
        // force too many real-world jobs to require two plink2 runs instead of
        // one.)

        const uint32_t setting_alleles_from_file = pcp->ref_allele_flag || pcp->alt1_allele_flag || pcp->ref_from_fa_fname;
        const char** allele_storage_backup = nullptr;
        uint32_t max_allele_slen_backup = max_allele_slen;
        uintptr_t* nonref_flags_backup = nullptr;
        uint32_t nonref_flags_was_null = (nonref_flags == nullptr);

        AltAlleleCt* refalt1_select = nullptr;
        if ((pcp->misc_flags & kfMiscMajRef) || setting_alleles_from_file) {
          if (loop_cats_idx + 1 < loop_cats_ct) {
            // --ref-allele and --alt1-allele may alter max_allele_slen and
            // allele_storage[]; --loop-cats doesn't like this.  Save a backup
            // copy.
            if (pcp->ref_allele_flag || pcp->alt1_allele_flag) {
              if (bigstack_end_alloc_kcp(raw_allele_ct, &allele_storage_backup)) {
                goto Plink2Core_ret_NOMEM;
              }
              memcpy(allele_storage_backup, allele_storage, raw_allele_ct * sizeof(intptr_t));
            }
            // nonref_flags may be altered by all four flags.
            if (nonref_flags) {
              if (bigstack_end_alloc_w(raw_variant_ctl, &nonref_flags_backup)) {
                goto Plink2Core_ret_NOMEM;
              }
              memcpy(nonref_flags_backup, nonref_flags, raw_variant_ctl * sizeof(intptr_t));
            }
          }
          const uintptr_t refalt1_word_ct = DivUp(2 * raw_variant_ct * sizeof(AltAlleleCt), kBytesPerWord);
          uintptr_t* refalt1_select_ul;
          // no need to track bigstack_end_mark before this is allocated, etc.,
          // due to the restriction to --make-pgen/--export
          if (bigstack_end_alloc_w(refalt1_word_ct, &refalt1_select_ul)) {
            goto Plink2Core_ret_NOMEM;
          }
          const uintptr_t alt_allele_vals = k1LU << (8 * sizeof(AltAlleleCt));
          const uintptr_t fill_word = ((~k0LU) / ((alt_allele_vals - 1) * (alt_allele_vals + 1))) * alt_allele_vals;
          for (uintptr_t widx = 0; widx < refalt1_word_ct; ++widx) {
            refalt1_select_ul[widx] = fill_word;
          }
          refalt1_select = R_CAST(AltAlleleCt*, refalt1_select_ul);
          const uint32_t not_all_nonref = !(pgfi.gflags & kfPgenGlobalAllNonref);
          if ((not_all_nonref || setting_alleles_from_file) && (!nonref_flags)) {
            if (bigstack_end_alloc_w(raw_variant_ctl, &nonref_flags)) {
              goto Plink2Core_ret_NOMEM;
            }
            pgfi.nonref_flags = nonref_flags;
            if (not_all_nonref) {
              ZeroWArr(raw_variant_ctl, nonref_flags);
            } else {
              SetAllBits(raw_variant_ct, nonref_flags);
            }
          }
          uintptr_t* previously_seen = nullptr;
          if (pcp->ref_allele_flag) {
            if (pcp->alt1_allele_flag) {
              if (bigstack_alloc_w(raw_variant_ctl, &previously_seen)) {
                goto Plink2Core_ret_NOMEM;
              }
            }
            reterr = SetRefalt1FromFile(variant_include, variant_ids, variant_allele_idxs, pcp->ref_allele_flag, raw_variant_ct, variant_ct, max_variant_id_slen, 0, (pcp->misc_flags / kfMiscRefAlleleForce) & 1, pcp->max_thread_ct, allele_storage, &max_allele_slen, refalt1_select, nonref_flags, previously_seen);
            if (reterr) {
              goto Plink2Core_ret_1;
            }
          }
          if (pcp->alt1_allele_flag) {
            reterr = SetRefalt1FromFile(variant_include, variant_ids, variant_allele_idxs, pcp->alt1_allele_flag, raw_variant_ct, variant_ct, max_variant_id_slen, 1, (pcp->misc_flags / kfMiscAlt1AlleleForce) & 1, pcp->max_thread_ct, allele_storage, &max_allele_slen, refalt1_select, nonref_flags, previously_seen);
            if (reterr) {
              goto Plink2Core_ret_1;
            }
            if (previously_seen) {
              BigstackReset(previously_seen);
            }
            // for sanity's sake, --maj-ref, --ref-allele/--alt1-allele, and
            // --ref-from-fa are mutually exclusive
            // (though --ref-from-fa + --alt1-allele may be permitted later)
          } else if (pcp->ref_from_fa_fname) {
            if (vpos_sortstatus & kfUnsortedVarBp) {
              logerrputs("Error: --ref-from-fa requires a sorted .pvar/.bim.  Retry this command after\nusing --make-pgen/--make-bed + --sort-vars to sort your data.\n");
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            reterr = RefFromFa(variant_include, variant_bps, variant_allele_idxs, allele_storage, cip, pcp->ref_from_fa_fname, max_allele_slen, (pcp->misc_flags / kfMiscRefFromFaForce) & 1, refalt1_select, nonref_flags);
            if (reterr) {
              goto Plink2Core_ret_1;
            }
          } else if (!setting_alleles_from_file) {
            // --maj-ref; misc_flags & kfMiscMajRef check may be needed later

            // Since this also sets ALT1 to the second-most-common allele, it
            // can't just subscribe to maj_alleles[].
            const uint64_t* main_allele_dosages = nonfounders? allele_dosages : founder_allele_dosages;
            const uint32_t skip_real_ref = not_all_nonref && (!(pcp->misc_flags & kfMiscMajRefForce));
            if (skip_real_ref && (!nonref_flags)) {
              logerrputs("Warning: --maj-ref has no effect, since no provisional reference alleles are\npresent.  (Did you forget to add the 'force' modifier?)\n");
            } else {
              uint32_t variant_uidx = 0;
              for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
                MovU32To1Bit(variant_include, &variant_uidx);
                if (skip_real_ref && IsSet(nonref_flags, variant_uidx)) {
                  continue;
                }
                const uint64_t* cur_allele_dosages = &(main_allele_dosages[variant_allele_idxs? variant_allele_idxs[variant_uidx] : (2 * variant_uidx)]);
                const uint32_t alt_ct_p1 = variant_allele_idxs? (variant_allele_idxs[variant_uidx + 1] - variant_allele_idxs[variant_uidx]) : 2;
                if (alt_ct_p1 == 2) {
                  // optimize common case
                  if (cur_allele_dosages[1] > cur_allele_dosages[0]) {
                    // assumes AltAlleleCt is unsigned char
                    R_CAST(uint16_t*, refalt1_select)[variant_uidx] = 1;
                    if (nonref_flags) {
                      SetBit(variant_uidx, nonref_flags);
                    }
                  }
                } else {
                  uint32_t new_ref_idx = (cur_allele_dosages[1] > cur_allele_dosages[0]);
                  uint32_t new_alt1_idx = 1 - new_ref_idx;
                  uint64_t ref_dosage = cur_allele_dosages[new_ref_idx];
                  uint64_t alt1_dosage = cur_allele_dosages[new_alt1_idx];
                  for (uint32_t alt_idx = 2; alt_idx < alt_ct_p1; ++alt_idx) {
                    const uint64_t cur_alt_dosage = cur_allele_dosages[alt_idx];
                    if (cur_alt_dosage > alt1_dosage) {
                      if (cur_alt_dosage > ref_dosage) {
                        alt1_dosage = ref_dosage;
                        ref_dosage = cur_alt_dosage;
                        new_alt1_idx = new_ref_idx;
                        new_ref_idx = alt_idx;
                      } else {
                        alt1_dosage = cur_alt_dosage;
                        new_alt1_idx = alt_idx;
                      }
                    }
                  }
                  if (new_ref_idx || (new_alt1_idx != 1)) {
                    refalt1_select[2 * variant_uidx] = new_ref_idx;
                    refalt1_select[2 * variant_uidx + 1] = new_alt1_idx;
                    if (nonref_flags) {
                      SetBit(variant_uidx, nonref_flags);
                    }
                  }
                }
              }
            }
          }
        }
        // founder_allele_dosages no longer needed
        // allele_dosages only needed in trim-alts case
        // todo: trim-alts does NOT need to be dosage-sensitive when we're
        //   erasing dosage.  may want a bitarray to handle that case; and once
        //   that's implemented, make dosage-preserving trim-alts also use that
        //   (get rid of its allele_dosages[] dependency).
        if (make_plink2_flags & kfMakePlink2TrimAlts) {
          BigstackReset(bigstack_mark_founder_allele_dosages);
        } else {
          BigstackReset(bigstack_mark_allele_dosages);
        }

        uint32_t* new_sample_idx_to_old = nullptr;
        if (pcp->sample_sort_flags & (kfSortNatural | kfSortAscii | kfSortFile)) {
          if (sample_ct < 2) {
            logerrputs("Warning: Skipping --sample-sort since <2 samples are present.\n");
          } else {
            if (pcp->sample_sort_flags & kfSortFile) {
              reterr = SampleSortFileMap(sample_include, &pii.sii, pcp->sample_sort_fname, raw_sample_ct, sample_ct, &new_sample_idx_to_old);
              if (reterr) {
                goto Plink2Core_ret_1;
              }
            } else {
              // probably more efficient to have --make-{bed,pgen,bpgen}
              // perform an unfiltered load?  but we should have compute power
              // to spare here, so keep the code simpler for now
              char* sorted_xidbox_tmp;
              uintptr_t max_xid_blen;
              reterr = SortedXidboxInitAlloc(sample_include, &pii.sii, sample_ct, 0, pii.sii.sids? kfXidModeFidIidSid : kfXidModeFidIid, (pcp->sample_sort_flags == kfSortNatural), &sorted_xidbox_tmp, &new_sample_idx_to_old, &max_xid_blen);
              if (reterr) {
                goto Plink2Core_ret_1;
              }
              BigstackReset(sorted_xidbox_tmp);
            }
            logprintf("--indiv-sort: %u samples reordered.\n", sample_ct);
          }
        }

        if (covar_ct && ((pcp->command_flags1 & (kfCommand1Exportf | kfCommand1WriteCovar)) || ((pcp->command_flags1 & kfCommand1MakePlink2) && (make_plink2_flags & (kfMakeBed | kfMakeFam | kfMakePgen | kfMakePsam))))) {
          reterr = WriteCovar(sample_include, &pii, sex_nm, sex_male, pheno_cols, pheno_names, covar_cols, covar_names, new_sample_idx_to_old, sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, pcp->write_covar_flags, outname, outname_end);
          if (reterr) {
            goto Plink2Core_ret_1;
          }
        } else if (pcp->command_flags1 & kfCommand1WriteCovar) {
          logerrputs("Warning: Skipping --write-covar, since no covariates are loaded.\n");
        }

        if (pcp->command_flags1 & kfCommand1MakePlink2) {
          // todo: unsorted case (--update-chr, etc.)
          if (pcp->sort_vars_flags != kfSort0) {
            reterr = MakePlink2Vsort(xheader, sample_include, &pii, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, allele_dosages, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, info_reload_slen? pvarname : nullptr, variant_cms, chr_idxs, xheader_blen, info_flags, raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, pcp->max_thread_ct, pcp->hard_call_thresh, pcp->dosage_erase_thresh, make_plink2_flags, (pcp->sort_vars_flags == kfSortNatural), pcp->pvar_psam_flags, &simple_pgr, outname, outname_end);
          } else {
            if (vpos_sortstatus & kfUnsortedVarBp) {
              logerrputs("Warning: Variants are not sorted by position.  Consider rerunning with the\n--sort-vars flag added to remedy this.\n");
            }
            reterr = MakePlink2NoVsort(xheader, sample_include, &pii, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, allele_dosages, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, info_reload_slen? pvarname : nullptr, variant_cms, xheader_blen, info_flags, raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, pcp->max_thread_ct, pcp->hard_call_thresh, pcp->dosage_erase_thresh, make_plink2_flags, pcp->pvar_psam_flags, pgr_alloc_cacheline_ct, &pgfi, &simple_pgr, outname, outname_end);
          }
          if (reterr) {
            goto Plink2Core_ret_1;
          }
          BigstackReset(bigstack_mark_allele_dosages);
        }

        if (pcp->command_flags1 & kfCommand1Exportf) {
          reterr = Exportf(sample_include, &pii, sex_nm, sex_male, pheno_cols, pheno_names, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, info_reload_slen? pvarname : nullptr, variant_cms, xheader_blen, info_flags, raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, pcp->max_thread_ct, make_plink2_flags, pcp->exportf_flags, pcp->exportf_id_paste, pcp->exportf_id_delim, pcp->exportf_bits, pgr_alloc_cacheline_ct, xheader, &pgfi, &simple_pgr, outname, outname_end);
          if (reterr) {
            goto Plink2Core_ret_1;
          }
        }

        if (allele_storage_backup) {
          max_allele_slen = max_allele_slen_backup;
          memcpy(allele_storage, allele_storage_backup, raw_allele_ct * sizeof(intptr_t));
        }
        if (nonref_flags_backup) {
          memcpy(nonref_flags, nonref_flags_backup, raw_variant_ctl * sizeof(intptr_t));
        } else if (nonref_flags_was_null) {
          nonref_flags = nullptr;
          pgfi.nonref_flags = nullptr;
        }
      }
      BigstackReset(bigstack_mark_allele_dosages);

      if (pcp->command_flags1 & kfCommand1LdPrune) {
        if ((pcp->ld_info.prune_flags & kfLdPruneWindowBp) && (vpos_sortstatus & kfUnsortedVarBp)) {
          logerrputs("Error: When the window size is in kb units, LD-based pruning requires a sorted\n.pvar/.bim.  Retry this command after using --make-pgen/--make-bed +\n--sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        reterr = LdPrune(variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, maj_alleles, allele_freqs, founder_info, sex_male, &(pcp->ld_info), raw_variant_ct, variant_ct, raw_sample_ct, founder_ct, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & kfCommand1Ld) {
        reterr = LdConsole(variant_include, cip, variant_ids, variant_allele_idxs, allele_storage, founder_info, sex_nm, sex_male, &(pcp->ld_info), variant_ct, raw_sample_ct, founder_ct, &simple_pgr);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & kfCommand1Score) {
        reterr = ScoreReport(sample_include, &pii.sii, sex_male, pheno_cols, pheno_names, variant_include, cip, variant_ids, variant_allele_idxs, allele_storage, allele_freqs, &(pcp->score_info), raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_variant_id_slen, pcp->xchr_model, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
      // eventually check for nonzero pheno_ct here?

      if (pcp->command_flags1 & kfCommand1Glm) {
        reterr = GlmMain(sample_include, &pii.sii, sex_nm, sex_male, pheno_cols, pheno_names, covar_cols, covar_names, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, maj_alleles, allele_storage, &(pcp->glm_info), &(pcp->adjust_info), &(pcp->aperm), pcp->glm_local_covar_fname, pcp->glm_local_pvar_fname, pcp->glm_local_psam_fname, raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, raw_variant_ct, variant_ct, max_variant_id_slen, max_allele_slen, pcp->xchr_model, pcp->ci_size, pcp->vif_thresh, pcp->pfilter, pcp->output_min_p, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, &simple_pgr, outname, outname_end);
        if (reterr) {
          goto Plink2Core_ret_1;
        }
      }
    Plink2Core_early_complete:
      if (++loop_cats_idx == loop_cats_ct) {
        break;
      }
      BigstackDoubleReset(bigstack_mark_varfilter, bigstack_end_mark_varfilter);
    }
  }
  while (0) {
  Plink2Core_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Plink2Core_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  Plink2Core_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  Plink2Core_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 Plink2Core_ret_1:
  if (loop_cats_pheno_col) {
    vecaligned_free_cond(loop_cats_pheno_col->nonmiss);
  }
  CleanupPhenoCols(covar_ct, covar_cols);
  CleanupPhenoCols(pheno_ct, pheno_cols);
  free_cond(covar_names);
  free_cond(pheno_names);
  if (CleanupPgr(&simple_pgr) && (!reterr)) {
    reterr = kPglRetReadFail;
  }
  if (CleanupPgfi(&pgfi) && (!reterr)) {
    reterr = kPglRetReadFail;
  }
  // no BigstackReset() needed?
  return reterr;
}

PglErr ZstDecompress(const char* in_fname, const char* out_fname) {
  // Since this needs to be able to dump the decompressed data and nothing but
  // the decompressed data to stdout, we have to duplicate a bit of
  // plink2_decompress code and strip out printing/logging.

  // Strictly speaking, this can currently decompress gzipped files too, but
  // that's not its purpose.
  gzFile gz_infile = gzopen(in_fname, FOPEN_RB);
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (!gz_infile) {
      fprintf(stderr, kErrprintfFopen, in_fname);
      goto ZstDecompress_ret_OPEN_FAIL;
    }
    if (gzbuffer(gz_infile, 131072)) {
      goto ZstDecompress_ret_NOMEM;
    }
    if (out_fname) {
      outfile = fopen(out_fname, FOPEN_WB);
      if (!outfile) {
        fprintf(stderr, kErrprintfFopen, out_fname);
        goto ZstDecompress_ret_OPEN_FAIL;
      }
    } else {
      outfile = stdout;
    }
    unsigned char* buf = R_CAST(unsigned char*, g_textbuf);
    while (1) {
      const int32_t bytes_read = gzread(gz_infile, buf, kTextbufMainSize);
      if (bytes_read <= 0) {
        if (!bytes_read) {
          break;
        }
        goto ZstDecompress_ret_READ_FAIL;
      }
      if (!fwrite_unlocked(buf, bytes_read, 1, outfile)) {
        goto ZstDecompress_ret_WRITE_FAIL;
      }
      fflush(outfile);
    }
    if (gzclose_null(&gz_infile)) {
      goto ZstDecompress_ret_READ_FAIL;
    }
    if (out_fname) {
      if (fclose_null(&outfile)) {
        goto ZstDecompress_ret_WRITE_FAIL;
      }
    }
  }
  // we exit from main() immediately, so need to print nomem/read/write error
  // messages here
  while (0) {
  ZstDecompress_ret_NOMEM:
    // in this exceedingly unlikely case, the --memory flag doesn't help, so
    // print a different message
    fputs("Error: Out of memory.\n", stderr);
    reterr = kPglRetNomem;
    break;
  ZstDecompress_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ZstDecompress_ret_READ_FAIL:
    fputs(kErrstrRead, stderr);
    reterr = kPglRetReadFail;
    break;
  ZstDecompress_ret_WRITE_FAIL:
    fputs(kErrstrWrite, stderr);
    reterr = kPglRetWriteFail;
    break;
  }
  if (out_fname) {
    fclose_cond(outfile);
  }
  gzclose_cond(gz_infile);
  return reterr;
}

PglErr Alloc2col(const char* const* sources, const char* flagname_p, uint32_t param_ct, TwoColParams** tcbuf) {
  uint32_t fname_blen = strlen(sources[0]) + 1;
  if (fname_blen > kPglFnamesize) {
    logerrprintf("Error: --%s filename too long.\n", flagname_p);
    return kPglRetOpenFail;
  }
  if (pgl_malloc(offsetof(TwoColParams, fname) + fname_blen, tcbuf)) {
    return kPglRetNomem;
  }
  memcpy((*tcbuf)->fname, sources[0], fname_blen);
  (*tcbuf)->skip_ct = 0;
  (*tcbuf)->skipchar = '\0';
  if (param_ct > 1) {
    if (ScanPosintDefcap(sources[1], &((*tcbuf)->colx))) {
      logerrprintf("Error: Invalid --%s column number.\n", flagname_p);
      return kPglRetInvalidCmdline;
    }
    if (param_ct > 2) {
      if (ScanPosintDefcap(sources[2], &((*tcbuf)->colid))) {
        logerrprintf("Error: Invalid --%s variant ID column number.\n", flagname_p);
        return kPglRetInvalidCmdline;
      }
      if (param_ct == 4) {
        char cc = sources[3][0];
        if ((cc < '0') || (cc > '9')) {
          cc = ExtractCharParam(sources[3]);
          if (!cc) {
            goto Alloc2col_invalid_skip;
          }
          (*tcbuf)->skipchar = cc;
        } else {
          if (ScanUintDefcap(sources[3], &((*tcbuf)->skip_ct))) {
          Alloc2col_invalid_skip:
            logerrprintf("Error: Invalid --%s skip parameter.  This needs to either be a\nsingle character (usually '#') which, when present at the start of a line,\nindicates it should be skipped; or the number of initial lines to skip.  (Note\nthat in shells such as bash, '#' is a special character that must be\nsurrounded by single- or double-quotes to be parsed correctly.)\n", flagname_p);
            return kPglRetInvalidCmdline;
          }
        }
      }
    } else {
      (*tcbuf)->colid = 1;
    }
    if ((*tcbuf)->colx == (*tcbuf)->colid) {
      logerrprintf("Error: Column numbers for --%s cannot be equal.\n%s", flagname_p, errstr_append);
      return kPglRetInvalidCmdline;
    }
  } else {
    (*tcbuf)->colx = 2;
    (*tcbuf)->colid = 1;
  }
  return kPglRetSuccess;
}

PglErr AllocAndFlattenCommaDelim(const char* const* sources, uint32_t param_ct, char** flattened_buf_ptr) {
  uint32_t tot_blen = 1;
  for (uint32_t param_idx = 0; param_idx < param_ct; ++param_idx) {
    const char* cur_param_iter = sources[param_idx];
    while (1) {
      while (*cur_param_iter == ',') {
        ++cur_param_iter;
      }
      const char* cur_token_end = strchr(cur_param_iter, ',');
      if (!cur_token_end) {
        break;
      }
      tot_blen += 1 + S_CAST(uintptr_t, cur_token_end - cur_param_iter);
      cur_param_iter = &(cur_token_end[1]);
    }
    tot_blen += 1 + strlen(cur_param_iter);
  }
  char* write_iter;
  if (pgl_malloc(tot_blen, &write_iter)) {
    return kPglRetNomem;
  }
  *flattened_buf_ptr = write_iter;
  for (uint32_t param_idx = 0; param_idx < param_ct; ++param_idx) {
    const char* cur_param_iter = sources[param_idx];
    while (1) {
      while (*cur_param_iter == ',') {
        ++cur_param_iter;
      }
      const char* cur_token_end = strchr(cur_param_iter, ',');
      if (!cur_token_end) {
        break;
      }
      write_iter = memcpyax(write_iter, cur_param_iter, cur_token_end - cur_param_iter, '\0');
      cur_param_iter = &(cur_token_end[1]);
    }
    write_iter = strcpyax(write_iter, cur_param_iter, '\0');
  }
  *write_iter = '\0';
  return kPglRetSuccess;
}

void PrintVer() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}

uint32_t CmdlineSingleChr(const ChrInfo* cip, MiscFlags misc_flags) {
  if ((misc_flags & (kfMiscAutosomeOnly | kfMiscAutosomePar)) || (!cip->is_include_stack)) {
    return 0;
  }
  const uint32_t main_chr_ct = PopcountWords(cip->chr_mask, kChrExcludeWords) + PopcountWord(cip->chr_mask[kChrMaskWords - 1]);
  if (main_chr_ct > 1) {
    return 0;
  }
  if (main_chr_ct == 1) {
    return (cip->incl_excl_name_stack == nullptr);
  }
  return cip->incl_excl_name_stack && (!(cip->incl_excl_name_stack->next));
}

void GetExportfTargets(const char* const* argvk, uint32_t param_ct, ExportfFlags* exportf_flags_ptr, IdpasteFlags* exportf_id_paste_ptr, uint32_t* format_param_idxs_ptr) {
  // does not error out if no format present, since this is needed for --recode
  // translation
  // supports multiple formats
  uint32_t format_param_idxs = 0;
  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
    const char* cur_modif = argvk[param_idx];
    const char* cur_modif2 = &(cur_modif[1]);
    ExportfFlags cur_format = kfExportf0;
    switch (*cur_modif) {
    case '2':
      if (!strcmp(cur_modif2, "3")) {
        cur_format = kfExportf23;
      }
      break;
    case 'A':
    case 'a':
      if (!cur_modif2[0]) {
        cur_format = kfExportfA;
      } else if (((cur_modif2[0] & 0xdf) == 'D') && (!cur_modif2[1])) {
        cur_format = kfExportfAD;
      } else if (!strcmp(cur_modif2, "-transpose")) {
        cur_format = kfExportfATranspose;
      }
      break;
    case 'b':
      if (!strcmp(cur_modif2, "eagle")) {
        cur_format = kfExportfBeagle;
      } else if (!strcmp(cur_modif2, "eagle-nomap")) {
        cur_format = kfExportfBeagleNomap;
      } else if ((!strcmp(cur_modif2, "gen-1.1")) || (!strcmp(cur_modif2, "gen_1.1"))) {
        cur_format = kfExportfBgen11;
      } else if ((!strcmp(cur_modif2, "gen-1.2")) || (!strcmp(cur_modif2, "gen_1.2"))) {
        cur_format = kfExportfBgen12;
      } else if ((!strcmp(cur_modif2, "gen-1.3")) || (!strcmp(cur_modif2, "gen_1.3"))) {
        cur_format = kfExportfBgen13;
      } else if (!strcmp(cur_modif2, "imbam")) {
        cur_format = kfExportfBimbam;
      } else if (!strcmp(cur_modif2, "imbam-1chr")) {
        cur_format = kfExportfBimbam1chr;
      }
      break;
    case 'c':
      if (!strcmp(cur_modif2, "ompound-genotypes")) {
        cur_format = kfExportfCompound;
      }
      break;
    case 'f':
      if (!strcmp(cur_modif2, "astphase")) {
        cur_format = kfExportfFastphase;
      } else if (!strcmp(cur_modif2, "astphase-1chr")) {
        cur_format = kfExportfFastphase1chr;
      }
      break;
    case 'h':
      if (!strcmp(cur_modif2, "aps")) {
        cur_format = kfExportfHaps;
        break;
      } else if (!strcmp(cur_modif2, "apslegend")) {
        cur_format = kfExportfHapsLegend;
        break;
      }
      // fall through
    case 'H':
      if ((cur_modif2[0] & 0xdf) == 'V') {
        if (!cur_modif2[1]) {
          cur_format = kfExportfHv;
        } else if (!strcmp(&(cur_modif2[1]), "-1chr")) {
          cur_format = kfExportfHv1chr;
        }
      }
      break;
    case 'i':
      if (!strcmp(cur_modif2, "nd-major-bed")) {
        cur_format = kfExportfIndMajorBed;
      }
      break;
    case 'l':
      if (!strcmp(cur_modif2, "gen")) {
        cur_format = kfExportfLgen;
      } else if (!strcmp(cur_modif2, "gen-ref")) {
        cur_format = kfExportfLgenRef;
      } else if (!strcmp(cur_modif2, "ist")) {
        cur_format = kfExportfList;
      }
      break;
    case 'o':
      if (!strcmp(cur_modif2, "xford")) {
        cur_format = kfExportfOxGen;
      }
      break;
    case 'p':
      if (!strcmp(cur_modif2, "ed")) {
        cur_format = kfExportfPed;
      }
      break;
    case 'r':
      if (!strcmp(cur_modif2, "list")) {
        cur_format = kfExportfRlist;
      }
      break;
    case 's':
      if (!strcmp(cur_modif2, "tructure")) {
        cur_format = kfExportfStructure;
      }
      break;
    case 't':
      if (!strcmp(cur_modif2, "ranspose")) {
        cur_format = kfExportfTranspose;
      }
      break;
    case 'v':
      if ((cur_modif2[0] == 'c') && (cur_modif2[1] == 'f')) {
        if (!cur_modif2[2]) {
          cur_format = kfExportfVcf;
        } else if ((!strcmp(&(cur_modif2[2]), "-fid")) || (!strcmp(&(cur_modif2[2]), "-iid"))) {
          snprintf(g_logbuf, kLogbufSize, "Note: --export 'v%s' modifier is deprecated.  Use 'vcf' + 'id-paste=%s'.\n", cur_modif2, &(cur_modif2[3]));
          cur_format = kfExportfVcf;
          *exportf_id_paste_ptr = (cur_modif2[3] == 'f')? kfIdpasteFid : kfIdpasteIid;
        }
      }
      break;
    }
    if (cur_format) {
      format_param_idxs |= 1U << param_idx;
      *exportf_flags_ptr |= cur_format;
    }
  }
  *format_param_idxs_ptr = format_param_idxs;
}

uint32_t VaridTemplateIsValid(const char* varid_str, const char* flagname_p) {
  const char* sptr = strchr(varid_str, '@');
  const char* sptr2 = strchr(varid_str, '#');
  if ((!sptr) || (!sptr2) || strchr(&(sptr[1]), '@') || strchr(&(sptr2[1]), '#')) {
    logerrprintfww("Error: The --%s template string requires exactly one '@' and one '#'.\n", flagname_p);
    return 0;
  }
  // snp/nonsnp is not sufficient for assigning unique IDs to unnamed 1000
  // Genomes phase 3 variants (see e.g. chr22:18078898).  So we now allow the
  // template string to include allele names, where '$r' = reference allele,
  // '$a' = alt1, and '$1'/'$2' refer to ref/alt1 in ASCII-sort order
  // (necessary for interoperation with plink1).
  // For now, either '$' must be entirely absent from the template string, or
  // '$r' and/or '$a' appear exactly once, or '$1' and '$2' both appear exactly
  // once.
  // probable todo: alternate naming scheme for long indels (e.g. first base,
  //   middle length, last base, like "i18n")
  // possible todo: some way to include alt2, etc. in name
  sptr = strchr(varid_str, '$');
  if (sptr) {
    sptr2 = &(sptr[1]);
    uint32_t first_allele_type_code = ctou32(*sptr2);
    if ((first_allele_type_code == 49) || (first_allele_type_code == 50)) {
      sptr2 = strchr(sptr2, '$');
      if ((!sptr2) || strchr(&(sptr2[1]), '$') || ((first_allele_type_code + ctou32(sptr2[1])) != 99)) {
      VaridTemplateIsValid_dollar_error:
        logerrprintfww("Error: The --%s template string requires either no instances of '$', exactly one instance of '$r' and/or '$a', or exactly one '$1' and one '$2'.\n", flagname_p);
        return 0;
      }
    } else {
      first_allele_type_code &= 0xdf;
      if ((first_allele_type_code != 65) && (first_allele_type_code != 82)) {
        goto VaridTemplateIsValid_dollar_error;
      }
      sptr2 = strchr(sptr2, '$');
      if (sptr2) {
        const uint32_t second_allele_type_code = ctou32(*(++sptr2)) & 0xdf;
        if (((first_allele_type_code + second_allele_type_code) != 147) || strchr(sptr2, '$')) {
          goto VaridTemplateIsValid_dollar_error;
        }
      }
    }
  }
  return 1;
}


static_assert(sizeof(int) == sizeof(int32_t), "main() assumes int and int32_t are synonymous.");
static_assert(!kChrOffsetX, "--autosome-num/--chr-set/--cow/etc. assume kChrOffsetX == 0.");
static_assert(kChrOffsetY == 1, "--chr-set/--cow/... assume kChrOffsetY == 1.");
static_assert(kChrOffsetXY == 2, "--chr-set/--cow/... assume kChrOffsetXY == 2.");
static_assert(kChrOffsetMT == 3, "--chr-set/--cow/... assume kChrOffsetMT == 3.");
#ifdef __cplusplus
}  // namespace plink2
#endif

int main(int argc, char** argv) {
#ifdef __cplusplus
  using namespace plink2;
#endif
  // special case, since it may dump to stdout
  if ((argc > 1) && ((!strcmp(argv[1], "--zst-decompress")) || (!strcmp(argv[1], "-zst-decompress")))) {
    if (argc == 2) {
      fprintf(stderr, "Error: Missing %s parameter.\n", argv[1]);
      return S_CAST(uint32_t, kPglRetInvalidCmdline);
    }
    for (int ii = 2; ii < argc; ++ii) {
      if (IsCmdlineFlag(argv[S_CAST(uint32_t, ii)])) {
        fprintf(stderr, "Error: %s cannot be used with other flags.\n", argv[1]);
        return S_CAST(uint32_t, kPglRetInvalidCmdline);
      }
    }
    if (argc > 4) {
      fprintf(stderr, "Error: %s accepts at most 2 parameters.\n", argv[1]);
      return S_CAST(uint32_t, kPglRetInvalidCmdline);
    }
    return S_CAST(uint32_t, ZstDecompress(argv[2], (argc == 4)? argv[3] : nullptr));
  }

  unsigned char* bigstack_ua = nullptr;
  Plink2CmdlineMeta pcm;
  PreinitPlink2CmdlineMeta(&pcm);
  const char* flagname_p = nullptr;
  char* king_cutoff_fprefix = nullptr;
  char* const_fid = nullptr;
  char* import_single_chr_str = nullptr;
  char* ox_missing_code = nullptr;
  char* vcf_dosage_import_field = nullptr;
  uint32_t* rseeds = nullptr;
  LlStr* file_delete_list = nullptr;
  uint32_t arg_idx = 0;
  uint32_t print_end_time = 0;
  uint32_t warning_errcode = 0;
  PglErr reterr = kPglRetSuccess;
  Plink2Cmdline pc;
  pc.filter_flags = kfFilter0;
  pc.dependency_flags = kfFilter0;
  pc.var_filter_exceptions_flattened = nullptr;
  pc.varid_template = nullptr;
  pc.missing_varid_match = nullptr;
  pc.varid_from = nullptr;
  pc.varid_to = nullptr;
  pc.varid_snp = nullptr;
  pc.varid_exclude_snp = nullptr;
  pc.pheno_fname = nullptr;
  pc.covar_fname = nullptr;
  pc.sample_sort_fname = nullptr;
  pc.keep_fnames = nullptr;
  pc.keepfam_fnames = nullptr;
  pc.remove_fnames = nullptr;
  pc.removefam_fnames = nullptr;
  pc.extract_fnames = nullptr;
  pc.exclude_fnames = nullptr;
  pc.freq_ref_binstr = nullptr;
  pc.freq_alt1_binstr = nullptr;
  pc.glm_local_covar_fname = nullptr;
  pc.glm_local_pvar_fname = nullptr;
  pc.glm_local_psam_fname = nullptr;
  pc.read_freq_fname = nullptr;
  pc.within_fname = nullptr;
  pc.catpheno_name = nullptr;
  pc.family_missing_catname = nullptr;
  pc.keep_cats_fname = nullptr;
  pc.keep_cat_names_flattened = nullptr;
  pc.keep_cat_phenoname = nullptr;
  pc.remove_cats_fname = nullptr;
  pc.remove_cat_names_flattened = nullptr;
  pc.remove_cat_phenoname = nullptr;
  pc.split_cat_phenonames_flattened = nullptr;
  pc.require_pheno_flattened = nullptr;
  pc.require_covar_flattened = nullptr;
  pc.vstd_flattened = nullptr;
  pc.quantnorm_flattened = nullptr;
  pc.covar_quantnorm_flattened = nullptr;
  pc.loop_cats_phenoname = nullptr;
  pc.ref_from_fa_fname = nullptr;
  pc.king_table_subset_fname = nullptr;
  pc.require_info_flattened = nullptr;
  pc.require_no_info_flattened = nullptr;
  pc.keep_fcol_fname = nullptr;
  pc.keep_fcol_flattened = nullptr;
  pc.keep_fcol_name = nullptr;
  pc.ref_allele_flag = nullptr;
  pc.alt1_allele_flag = nullptr;
  pc.update_name_flag = nullptr;
  InitRangeList(&pc.snps_range_list);
  InitRangeList(&pc.exclude_snps_range_list);
  InitRangeList(&pc.pheno_range_list);
  InitRangeList(&pc.covar_range_list);
  InitUpdateSex(&pc.update_sex_info);
  InitLd(&pc.ld_info);
  InitGlm(&pc.glm_info);
  InitScore(&pc.score_info);
  InitCmpExpr(&pc.keep_if_expr);
  InitCmpExpr(&pc.remove_if_expr);
  InitCmpExpr(&pc.extract_if_info_expr);
  InitCmpExpr(&pc.exclude_if_info_expr);
  AdjustFileInfo adjust_file_info;
  InitAdjust(&pc.adjust_info, &adjust_file_info);
  ChrInfo chr_info;
  if (InitChrInfo(&chr_info)) {
    goto main_ret_NOMEM_NOLOG;
  }

  {
    // standardize strtod() behavior
    // (not necessary any more since we use dtoa_g() instead)
    // setlocale(LC_NUMERIC, "C");

    uint32_t first_arg_idx;
    uint32_t flag_ct;
    reterr = CmdlineParsePhase1(ver_str, ver_str2, PROG_NAME_STR, notestr_null_calc2, kCmdlineFormatStr, errstr_append, kMaxFlagBlen, DispHelp, &argc, &argv, &pcm, &first_arg_idx, &flag_ct);
    if (reterr) {
      goto main_ret_NOLOG;
    }
    if (!flag_ct) {
      goto main_ret_NULL_CALC_0;
    }
    if (pgl_malloc(flag_ct * kMaxFlagBlen, &pcm.flag_buf) ||
        pgl_malloc(flag_ct * sizeof(int32_t), &pcm.flag_map)) {
      goto main_ret_NOMEM_NOLOG2;
    }

    // No modifications to argv past this point.
    const char* const* argvk = TO_CONSTCPCONSTP(argv);

    char* flagname_write_iter = pcm.flag_buf;
    uint32_t cur_flag_idx = 0;
    for (arg_idx = first_arg_idx; arg_idx < S_CAST(uint32_t, argc); ++arg_idx) {
      flagname_p = IsCmdlineFlagStart(argvk[arg_idx]);
      if (flagname_p) {
        const uint32_t flag_slen = strlen(flagname_p);
        switch (*flagname_p) {
        case '\0':
          // special case, since we reserve empty names for preprocessed flags
          fputs("Error: Unrecognized flag ('--').\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        case 'a':
          if (strequal_k(flagname_p, "aec", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "allow-extra-chr");
          } else if (strequal_k(flagname_p, "autosome-xy", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "autosome-par");
          } else if (strequal_k(flagname_p, "a1-allele", flag_slen)) {
            fputs("Warning: --a1-allele flag deprecated.  Use --alt1-allele instead.\n", stderr);
            snprintf(flagname_write_iter, kMaxFlagBlen, "alt1-allele");
          } else if (strequal_k(flagname_p, "a2-allele", flag_slen)) {
            fputs("Warning: --a2-allele flag deprecated.  Use --ref-allele instead.\n", stderr);
            snprintf(flagname_write_iter, kMaxFlagBlen, "ref-allele");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'b':
          if (strequal_k(flagname_p, "bed", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "pgen");
          } else if (strequal_k(flagname_p, "bim", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "pvar");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'c':
          if (strequal_k(flagname_p, "covar-number", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "covar-col-nums");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'e':
          if (strequal_k(flagname_p, "extract-if", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-if-info");
          } else if (strequal_k(flagname_p, "exclude-if", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "remove-if-info");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'f':
          if (strequal_k(flagname_p, "fam", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "psam");
          } else if (strequal_k(flagname_p, "filter-males", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-males");
          } else if (strequal_k(flagname_p, "filter-females", flag_slen)) {
            fputs("Note: --filter-females flag deprecated.  Use --keep-females or --remove-males\ninstead.\n", stdout);
            snprintf(flagname_write_iter, kMaxFlagBlen, "remove-males");
          } else if (strequal_k(flagname_p, "filter-founders", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-founders");
          } else if (strequal_k(flagname_p, "filter-nonfounders", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-nonfounders");
          } else if (strequal_k(flagname_p, "filter", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-fcol");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'h':
          if (strequal_k(flagname_p, "hound", flag_slen)) {
            // the creature type should be Dog.
            snprintf(flagname_write_iter, kMaxFlagBlen, "dog");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'k':
          if (strequal_k(flagname_p, "keep-clusters", flag_slen)) {
            fputs("Note: --keep-clusters flag deprecated.  Use --keep-cats instead.\n", stdout);
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-cats");
          } else if (strequal_k(flagname_p, "keep-cluster-names", flag_slen)) {
            fputs("Note: --keep-cluster-names flag deprecated.  Use --keep-cat-names instead.\n", stdout);
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-cat-names");
          } else if (strequal_k(flagname_p, "keep-if-info", flag_slen)) {
            fputs("Note: --keep-if-info renamed to --extract-if-info, for consistency with other\nsample/variant filters (keep/remove = sample filter; extract/exclude = variant\nfilter).\n", stdout);
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-if-info");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'l':
          if (strequal_k(flagname_p, "linear", flag_slen) || strequal_k(flagname_p, "logistic", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "glm");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'm':
          if (strequal_k(flagname_p, "min-ac", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "mac");
          } else if (strequal_k(flagname_p, "max-ac", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "max-mac");
          } else if (strequal_k(flagname_p, "make-bfile", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "make-bed");
          } else if (strequal_k(flagname_p, "make-bpfile", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "make-bpgen");
          } else if (strequal_k(flagname_p, "make-pfile", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "make-pgen");
          } else if (strequal_k(flagname_p, "missing_code", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "missing-code");
          } else if (strequal_k(flagname_p, "max-indv", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "thin-indiv-count");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'n':
          if (strequal_k(flagname_p, "num_threads", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "threads");
          } else if (strequal_k(flagname_p, "no-pheno", flag_slen) ||
                     strequal_k(flagname_p, "no-fam-pheno", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "no-psam-pheno");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'p':
          if (strequal_k(flagname_p, "prune", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "require-pheno");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'r':
          if (strequal_k(flagname_p, "recode", flag_slen)) {
            // special case: translate to "export ped" if no format specified
            const uint32_t param_ct = GetParamCt(argvk, argc, arg_idx);
            if (param_ct > 4) {
              fputs("Error: --recode accepts at most 4 parameters.\n", stderr);
              goto main_ret_INVALID_CMDLINE;
            }
            ExportfFlags dummy;
            IdpasteFlags dummy2;
            uint32_t format_param_idxs;
            GetExportfTargets(&(argvk[arg_idx]), param_ct, &dummy, &dummy2, &format_param_idxs);
            if (!format_param_idxs) {
              snprintf(flagname_write_iter, kMaxFlagBlen, "export ped");
            } else {
              snprintf(flagname_write_iter, kMaxFlagBlen, "export");
            }
          } else if (strequal_k(flagname_p, "remove-founders", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-founders");
          } else if (strequal_k(flagname_p, "remove-nonfounders", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-nonfounders");
          } else if (strequal_k(flagname_p, "remove-clusters", flag_slen)) {
            fputs("Note: --remove-clusters flag deprecated.  Use --remove-cats instead.\n", stdout);
            snprintf(flagname_write_iter, kMaxFlagBlen, "remove-cats");
          } else if (strequal_k(flagname_p, "remove-cluster-names", flag_slen)) {
            fputs("Note: --remove-cluster-names flag deprecated.  Use --remove-cat-names instead.\n", stdout);
            snprintf(flagname_write_iter, kMaxFlagBlen, "remove-cat-names");
          } else if (strequal_k(flagname_p, "remove-if-info", flag_slen)) {
            fputs("Note: --remove-if-info renamed to --exclude-if-info, for consistency with other\nsample/variant filters (keep/remove = sample filter; extract/exclude = variant\nfilter).\n", stdout);
            snprintf(flagname_write_iter, kMaxFlagBlen, "exclude-if-info");
          } else {
            goto main_flag_copy;
          }
          break;
        case 't':
          if (strequal_k(flagname_p, "thread-num", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "threads");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'v':
          if (strequal_k(flagname_p, "vcf-filter", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "var-filter");
          } else if (strequal_k(flagname_p, "vcf-min-qual", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "var-min-qual");
          } else {
            goto main_flag_copy;
          }
          break;
        default:
        main_flag_copy:
          memcpy(flagname_write_iter, flagname_p, flag_slen + 1);
        }
        flagname_write_iter = &(flagname_write_iter[kMaxFlagBlen]);
        pcm.flag_map[cur_flag_idx++] = arg_idx;
      }
    }
    char outname[kPglFnamesize];
    memcpy(outname, "plink2", 6);
    char* outname_end = nullptr;
    int32_t known_procs;
    reterr = CmdlineParsePhase2(ver_str, errstr_append, argvk, 6, kMaxFlagBlen, argc, flag_ct, &pcm, outname, &outname_end, &known_procs, &pc.max_thread_ct);
    if (reterr) {
      goto main_ret_NOLOG;
    }

    char pgenname[kPglFnamesize];
    char psamname[kPglFnamesize];
    char pvarname[kPglFnamesize];
    pgenname[0] = '\0';
    psamname[0] = '\0';
    pvarname[0] = '\0';
    InitPheno();
    cur_flag_idx = 0;
    pc.command_flags1 = kfCommand10;
    // uint64_t command_flags2 = 0;
    pc.misc_flags = kfMisc0;
    pc.pvar_psam_flags = kfPvarPsam0;
    pc.exportf_flags = kfExportf0;
    pc.sample_sort_flags = kfSort0;
    pc.sort_vars_flags = kfSort0;
    pc.grm_flags = kfGrm0;
    pc.pca_flags = kfPca0;
    pc.write_covar_flags = kfWriteCovar0;
    pc.pheno_transform_flags = kfPhenoTransform0;
    pc.fam_cols = kfFamCol13456;
    pc.exportf_id_paste = kfIdpaste0;
    pc.king_flags = kfKing0;
    pc.king_cutoff = -1;
    pc.king_table_filter = -DBL_MAX;
    pc.freq_rpt_flags = kfAlleleFreq0;
    pc.missing_rpt_flags = kfMissingRpt0;
    pc.geno_counts_flags = kfGenoCounts0;
    pc.hardy_flags = kfHardy0;
    pc.aperm.min = 6;
    pc.aperm.max = 1000000;
    pc.aperm.alpha = 0.0;
    pc.aperm.beta = 0.0001;
    pc.aperm.init_interval = 1.0;
    pc.aperm.interval_slope = 0.001;
    pc.ci_size = 0.0;

    // Default value is 1638 = 32768 / 20, and that's applied to imported
    // dosages when --hard-call-threshold is not specified.
    // However, when --make-{b}pgen is run on a dosage-containing dataset,
    // explicit --hard-call-threshold will cause the hardcall set to be
    // regenerated, and that won't happen without --hard-call-threshold.  So we
    // need to distinguish between --hard-call-threshold 0.1 and no flag.
    pc.hard_call_thresh = UINT32_MAX;

    pc.dosage_erase_thresh = 0;
    pc.pfilter = 2.0;  // make --pfilter 1 still filter out NAs
    pc.output_min_p = 0.0;
    pc.vif_thresh = 50.0;
    pc.mind_thresh = 1.0;
    pc.geno_thresh = 1.0;
    pc.hwe_thresh = 1.0;
    pc.mach_r2_min = 0.0;
    pc.mach_r2_max = 0.0;
    pc.min_maf = 0.0;
    pc.max_maf = 1.0;
    pc.thin_keep_prob = 1.0;
    pc.thin_keep_sample_prob = 1.0;
    pc.min_allele_dosage = 0;
    pc.max_allele_dosage = (~0LLU);
    pc.var_min_qual = -1;
    pc.update_sex_colm2 = 1;
    pc.new_variant_id_max_allele_slen = 23;
    pc.splitpar_bound1 = 0;
    pc.splitpar_bound2 = 0;
    pc.missing_pheno = -9;
    pc.from_bp = -1;
    pc.to_bp = -1;
    pc.window_bp = -1;
    pc.pca_ct = 0;
    pc.xchr_model = 2;
    pc.parallel_idx = 0;
    pc.parallel_tot = 1;
    pc.exportf_bits = 0;
    pc.mwithin_val = 1;
    pc.min_bp_space = 0;
    pc.thin_keep_ct = UINT32_MAX;
    pc.thin_keep_sample_ct = UINT32_MAX;
    pc.keep_fcol_num = 0;
    pc.exportf_id_delim = '\0';
    double import_dosage_certainty = 0.0;
    int32_t vcf_min_gq = -1;
    int32_t vcf_min_dp = -1;
    intptr_t malloc_size_mb = 0;
    LoadParams load_params = kfLoadParams0;
    Xload xload = kfXload0;
    uint32_t rseed_ct = 0;
    MakePlink2Flags make_plink2_flags = kfMake0;
    OxfordImportFlags oxford_import_flags = kfOxfordImport0;
    VcfHalfCall vcf_half_call = kVcfHalfCallDefault;
    char range_delim = '-';
    char id_delim = '\0';
    char idspace_to = '\0';
    char input_missing_geno_char = '0';
    char output_missing_geno_char = '.';
    ImportFlags import_flags = kfImport0;
    uint32_t aperm_present = 0;
    uint32_t notchr_present = 0;
    uint32_t permit_multiple_inclusion_filters = 0;
    uint32_t memory_require = 0;
    uint32_t randmem = 0;
    GenDummyInfo gendummy_info;
    InitGenDummy(&gendummy_info);
    Plink1DosageInfo plink1_dosage_info;
    InitPlink1Dosage(&plink1_dosage_info);
    do {
      flagname_p = &(pcm.flag_buf[cur_flag_idx * kMaxFlagBlen]);
      if (!(*flagname_p)) {
        // preprocessed; not relevant now, but will need --d later
        continue;
      }
      const char* flagname_p2 = &(flagname_p[1]);
      arg_idx = pcm.flag_map[cur_flag_idx];
      uint32_t param_ct = GetParamCt(argvk, argc, arg_idx);
      switch (*flagname_p) {
      case '1':
        if (*flagname_p2 == '\0') {
          pc.misc_flags |= kfMiscAffection01;
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'a':
        if (strequal_k_unsafe(flagname_p2, "llow-extra-chr")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (!strequal_k_unsafe(cur_modif, "0")) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --allow-extra-chr parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            chr_info.zero_extra_chrs = 1;
          }
          pc.misc_flags |= kfMiscAllowExtraChrs;
        } else if (strequal_k_unsafe(flagname_p2, "utosome")) {
          pc.misc_flags |= kfMiscAutosomeOnly;
          chr_info.is_include_stack = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "utosome-par")) {
          if (pc.misc_flags & kfMiscAutosomeOnly) {
            logerrputs("Error: --autosome-par cannot be used with --autosome.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.misc_flags |= kfMiscAutosomePar;
          chr_info.is_include_stack = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "llow-no-samples")) {
          logerrputs("Error: --allow-no-samples is retired.  (If you are performing a set of\noperations which doesn't require sample information, the sample file won't be\nloaded at all.)\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (strequal_k_unsafe(flagname_p2, "llow-no-vars")) {
          logerrputs("Error: --allow-no-vars is retired.  (If you are performing a set of operations\nwhich doesn't require variant information, the variant file won't be loaded at\nall.)\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (strequal_k_unsafe(flagname_p2, "djust")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "gc", cur_modif_slen)) {
              pc.adjust_info.flags |= kfAdjustGc;
            } else if (strequal_k(cur_modif, "log10", cur_modif_slen)) {
              pc.adjust_info.flags |= kfAdjustLog10;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.adjust_info.flags |= kfAdjustZs;
            } else if (strequal_k(cur_modif, "qq-plot", cur_modif_slen)) {
              logerrputs("Error: 'qq-plot' modifier retired.  Use e.g. \"--adjust cols=+qq\" instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (pc.adjust_info.flags & kfAdjustColAll) {
                logerrputs("Error: Multiple --adjust cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0unadj\0gc\0qq\0bonf\0holm\0sidakss\0sidaksd\0fdrbh\0fdrby\0", "adjust", kfAdjustColChrom, kfAdjustColDefault, 1, &pc.adjust_info.flags);
              if (reterr) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --adjust parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.adjust_info.flags & kfAdjustColAll)) {
            pc.adjust_info.flags |= kfAdjustColDefault;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-file")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 7)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &adjust_file_info.fname);
          if (reterr) {
            goto main_ret_1;
          }
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "gc", cur_modif_slen)) {
              adjust_file_info.base.flags |= kfAdjustGc;
            } else if (strequal_k(cur_modif, "log10", cur_modif_slen)) {
              adjust_file_info.base.flags |= kfAdjustLog10;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              adjust_file_info.base.flags |= kfAdjustZs;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (adjust_file_info.base.flags & kfAdjustColAll) {
                logerrputs("Error: Multiple --adjust-file cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0unadj\0gc\0qq\0bonf\0holm\0sidakss\0sidaksd\0fdrbh\0fdrby\0", "adjust-file", kfAdjustColChrom, kfAdjustColDefault, 1, &adjust_file_info.base.flags);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "test=", cur_modif_slen)) {
              reterr = CmdlineAllocString(&(cur_modif[5]), "--adjust-file test=", kMaxIdSlen, &adjust_file_info.test_name);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (strequal_k(cur_modif, "input-log10", cur_modif_slen)) {
              adjust_file_info.base.flags |= kfAdjustInputLog10;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --adjust-file parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(adjust_file_info.base.flags & kfAdjustColAll)) {
            adjust_file_info.base.flags |= kfAdjustColDefault;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-chr-field")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &adjust_file_info.chr_field);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-alt-field")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &adjust_file_info.alt_field);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-pos-field")) {
          if (!adjust_file_info.fname) {
            logerrputs("Error: --adjust-pos-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &adjust_file_info.pos_field);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-id-field")) {
          if (!adjust_file_info.fname) {
            logerrputs("Error: --adjust-id-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &adjust_file_info.id_field);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-ref-field")) {
          if (!adjust_file_info.fname) {
            logerrputs("Error: --adjust-ref-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &adjust_file_info.ref_field);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-test-field")) {
          if (!adjust_file_info.fname) {
            logerrputs("Error: --adjust-test-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &adjust_file_info.test_field);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-p-field")) {
          if (!adjust_file_info.fname) {
            logerrputs("Error: --adjust-p-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &adjust_file_info.p_field);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "perm")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 6)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (ScanPosintDefcap(cur_modif, &pc.aperm.min)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm min permutation count '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          ++pc.aperm.min;
          if (param_ct > 1) {
            cur_modif = argvk[arg_idx + 2];
            if (ScanPosintCapped(cur_modif, kApermMax, &pc.aperm.max)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm max permutation count '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (pc.aperm.min >= pc.aperm.max) {
            logerrputs("Error: --aperm min permutation count must be smaller than max.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          aperm_present = 1;
          if (param_ct > 2) {
            cur_modif = argvk[arg_idx + 3];
            if (!ScanadvDouble(cur_modif, &pc.aperm.alpha)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm alpha threshold '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (param_ct > 3) {
              cur_modif = argvk[arg_idx + 4];
              if ((!ScanadvDouble(cur_modif, &pc.aperm.beta)) || (pc.aperm.beta <= 0.0)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm beta '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (param_ct > 4) {
                cur_modif = argvk[arg_idx + 5];
                if (!ScanadvDouble(cur_modif, &pc.aperm.init_interval)) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm initial pruning interval '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                if ((pc.aperm.init_interval < 1.0) || (pc.aperm.init_interval > 1000000.0)) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm initial pruning interval '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                if (param_ct == 6) {
                  cur_modif = argvk[arg_idx + 6];
                  if (!ScanadvDouble(cur_modif, &pc.aperm.interval_slope) || (pc.aperm.interval_slope < 0.0) || (pc.aperm.interval_slope > 1.0)) {
                    snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm pruning interval slope '%s'.\n", cur_modif);
                    goto main_ret_INVALID_CMDLINE_WWA;
                  }
                }
              }
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "utosome-num")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t autosome_ct;
          if (ScanPosintCapped(cur_modif, kMaxChrTextnum, &autosome_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --autosome-num parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          // see plink2_common FinalizeChrset()
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.autosome_ct = autosome_ct;
          // assumes first code is X
          chr_info.xymt_codes[0] = autosome_ct + 1;
          for (uint32_t xymt_idx = 1; xymt_idx < kChrOffsetCt; ++xymt_idx) {
            // bugfix: this needs to be UINT32_MAXM1, not UINT32_MAX, for
            // GetChrCode() to work properly
            chr_info.xymt_codes[xymt_idx] = UINT32_MAXM1;
          }
          chr_info.haploid_mask[0] = 0;
          SetBit(autosome_ct + 1, chr_info.haploid_mask);
        } else if (strequal_k_unsafe(flagname_p2, "lt1-allele")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 5)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* const* sources = &(argvk[arg_idx + 1]);
          if (!strcmp(sources[0], "force")) {
            --param_ct;
            if (!param_ct) {
              logputs("Error: Invalid --alt1-allele parameter sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.misc_flags |= kfMiscAlt1AlleleForce;
            ++sources;
          }
          reterr = Alloc2col(sources, flagname_p, param_ct, &pc.alt1_allele_flag);
          if (reterr) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "llow-no-sex")) {
          logputs("Note: --allow-no-sex no longer has any effect.  (Missing-sex samples are\nautomatically excluded from association analysis when sex is a covariate, and\ntreated normally otherwise.)\n");
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'b':
        if (strequal_k_unsafe(flagname_p2, "file")) {
          if (xload) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (CheckExtraParam(&(argvk[arg_idx]), "vzs", &fname_modif_idx)) {
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* fname_prefix = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname_prefix);
          // pathological case bugfix (1 Feb 2018): need to subtract 1 more, to
          // avoid buffer overflow in the case we rename and append '~'.
          if (slen > (kPglFnamesize - 10)) {
            // could use kPglFnamesize - 2 - 3 * param_ct, but that's pointless
            logerrputs("Error: --bfile parameter too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          snprintf(memcpya(pgenname, fname_prefix, slen), 9, ".bed");
          snprintf(memcpya(psamname, fname_prefix, slen), 9, ".fam");
          char* bimname_end = memcpya(pvarname, fname_prefix, slen);
          bimname_end = Stpcpy(bimname_end, ".bim");
          if (param_ct == 2) {
            snprintf(bimname_end, 5, ".zst");
          }
          load_params |= kfLoadParamsPfileAll;
        } else if (strequal_k_unsafe(flagname_p2, "pfile")) {
          if (load_params || xload) {
            // currently only possible with --bcf, --bfile
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (CheckExtraParam(&(argvk[arg_idx]), "vzs", &fname_modif_idx)) {
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* fname_prefix = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname_prefix);
          if (slen > (kPglFnamesize - 10)) {
            logerrputs("Error: --bpfile parameter too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          snprintf(memcpya(pgenname, fname_prefix, slen), 9, ".pgen");
          snprintf(memcpya(psamname, fname_prefix, slen), 9, ".fam");
          char* bimname_end = memcpya(pvarname, fname_prefix, slen);
          bimname_end = Stpcpy(bimname_end, ".bim");
          if (param_ct == 2) {
            snprintf(bimname_end, 5, ".zst");
          }
          load_params |= kfLoadParamsPfileAll;
        } else if (strequal_k_unsafe(flagname_p2, "iallelic-only")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "strict")) {
              pc.misc_flags |= kfMiscBiallelicOnlyStrict;
            } else if (!strcmp(cur_modif, "list")) {
              pc.misc_flags |= kfMiscBiallelicOnlyList;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --biallelic-only parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.misc_flags |= kfMiscBiallelicOnly;
          logerrputs("Error: --biallelic-only is not implemented yet.\n");
          reterr = kPglRetNotYetSupported;
          goto main_ret_1;
        } else if (strequal_k_unsafe(flagname_p2, "cf")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 2) {
            const char* cur_modif = argvk[arg_idx + 2];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (!StrStartsWith(cur_modif, "dosage=", cur_modif_slen)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --bcf parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            reterr = CmdlineAllocString(&(cur_modif[strlen("dosage=")]), argvk[arg_idx], 4095, &vcf_dosage_import_field);
            if (reterr) {
              goto main_ret_1;
            }
            if (!IsAlphanumeric(vcf_dosage_import_field)) {
              logerrputs("Error: --bcf dosage= parameter is not alphanumeric.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            if (!strcmp(vcf_dosage_import_field, "GT")) {
              logerrputs("Error: --bcf dosage= parameter cannot be 'GT'.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_modif);
          if (slen > kPglFnamesize - 1) {
            logerrputs("Error: --bcf filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_modif, slen + 1);
          xload = kfXloadBcf;
        } else if (strequal_k_unsafe(flagname_p2, "gen")) {
          if (load_params || xload) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "snpid-chr")) {
              oxford_import_flags |= kfOxfordImportBgenSnpIdChr;
            } else if (!strcmp(cur_modif, "ref-first")) {
              oxford_import_flags |= kfOxfordImportRefFirst;
            } else if (!strcmp(cur_modif, "ref-second")) {
              oxford_import_flags |= kfOxfordImportRefSecond;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --bgen parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const char* cur_fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_fname);
          if (slen > kPglFnamesize - 1) {
            logerrputs("Error: --bgen filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_fname, slen + 1);
          xload = kfXloadOxBgen;
        } else if (strequal_k_unsafe(flagname_p2, "p-space")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (ScanPosintDefcap(cur_modif, &pc.min_bp_space)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --bp-space minimum bp distance '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'c':
        if (strequal_k_unsafe(flagname_p2, "hr")) {
          if (pc.misc_flags & (kfMiscAutosomePar | kfMiscAutosomeOnly)) {
            logerrputs("Error: --chr cannot be used with --autosome{-par}.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = ParseChrRanges(&(argvk[arg_idx]), flagname_p, errstr_append, param_ct, (pc.misc_flags / kfMiscAllowExtraChrs) & 1, 0, '-', &chr_info, chr_info.chr_mask);
          if (reterr) {
            goto main_ret_1;
          }
          chr_info.is_include_stack = 1;
        } else if (strequal_k_unsafe(flagname_p2, "ovar")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &pc.covar_fname);
          if (reterr) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ovar-col-nums")) {
          // requires --covar or --pheno, but --pheno hasn't been parsed yet so
          // we don't enforce the condition here
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, range_delim, &pc.covar_range_list);
          if (reterr) {
            goto main_ret_1;
          }
          pc.misc_flags |= kfMiscCovarColNums;
        } else if (strequal_k_unsafe(flagname_p2, "ovar-name")) {
          if (pc.covar_range_list.name_ct) {
            logerrputs("Error: --covar-name can't be used with --covar-col-nums.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // can now be used without --covar
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.covar_range_list);
          if (reterr) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "onst-fid")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(param_ct? argvk[arg_idx + 1] : "0", argvk[arg_idx], kMaxIdSlen, &const_fid);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "i")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (!ScanadvDouble(argvk[arg_idx + 1], &pc.ci_size)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --ci parameter '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if ((pc.ci_size < 0.01) || (pc.ci_size >= 1.0)) {
            logerrputs("Error: --ci confidence interval size must be in [0.01, 1).\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ondition")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (!strcmp("dominant", argvk[arg_idx + 2])) {
              pc.glm_info.flags |= kfGlmConditionDominant;
            } else if (!strcmp("recessive", argvk[arg_idx + 2])) {
              pc.glm_info.flags |= kfGlmConditionRecessive;
            } else {
              fname_modif_idx = 2;
              if (!strcmp("dominant", argvk[arg_idx + 1])) {
                pc.glm_info.flags |= kfGlmConditionDominant;
              } else if (!strcmp("recessive", argvk[arg_idx + 1])) {
                pc.glm_info.flags |= kfGlmConditionRecessive;
              } else {
                logerrputs("Error: Invalid --condition parameter sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            }
          }
          reterr = CmdlineAllocString(argvk[arg_idx + fname_modif_idx], argvk[arg_idx], kMaxIdSlen, &pc.glm_info.condition_varname);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ondition-list")) {
          if (pc.glm_info.condition_varname) {
            logerrputs("Error: --condition-list cannot be used with --condition.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (!strcmp("dominant", argvk[arg_idx + 2])) {
              pc.glm_info.flags |= kfGlmConditionDominant;
            } else if (!strcmp("recessive", argvk[arg_idx + 2])) {
              pc.glm_info.flags |= kfGlmConditionRecessive;
            } else {
              fname_modif_idx = 2;
              if (!strcmp("dominant", argvk[arg_idx + 1])) {
                pc.glm_info.flags |= kfGlmConditionDominant;
              } else if (!strcmp("recessive", argvk[arg_idx + 1])) {
                pc.glm_info.flags |= kfGlmConditionRecessive;
              } else {
                logerrputs("Error: Invalid --condition-list parameter sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            }
          }
          reterr = AllocFname(argvk[arg_idx + fname_modif_idx], flagname_p, 0, &pc.glm_info.condition_list_fname);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ow")) {
          if (chr_info.chrset_source) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.chrset_source = kChrsetSourceCmdline;
          // initialize here instead of FinalizeChrset(), to simplify
          // ReadChrsetHeaderLine()
          chr_info.autosome_ct = 29;
          chr_info.xymt_codes[0] = 30;
          chr_info.xymt_codes[1] = 31;
          chr_info.xymt_codes[2] = UINT32_MAXM1;
          chr_info.xymt_codes[3] = 33;
          chr_info.xymt_codes[4] = UINT32_MAXM1;
          chr_info.xymt_codes[5] = UINT32_MAXM1;
#ifdef __LP64__
          chr_info.haploid_mask[0] = 0x2c0000000LLU;
#else
          chr_info.haploid_mask[0] = 0xc0000000U;
          chr_info.haploid_mask[1] = 2;
#endif
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "hr-set")) {
          if (chr_info.chrset_source) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 5)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          int32_t signed_autosome_ct;
          if (ScanIntAbsBounded(cur_modif, kMaxChrTextnum, &signed_autosome_ct) || (!signed_autosome_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --chr-set parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          // see plink2_common FinalizeChrset()
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.haploid_mask[0] = 0;
          if (signed_autosome_ct < 0) {
            // haploid
            if (param_ct > 1) {
              logerrputs("Error: --chr-set does not accept multiple parameters in haploid mode.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            const uint32_t autosome_ct = -signed_autosome_ct;
            chr_info.autosome_ct = autosome_ct;
            for (uint32_t xymt_idx = 0; xymt_idx < kChrOffsetCt; ++xymt_idx) {
              chr_info.xymt_codes[xymt_idx] = UINT32_MAXM1;
            }
            SetAllBits(autosome_ct + 1, chr_info.haploid_mask);
          } else {
            const uint32_t autosome_ct = signed_autosome_ct;
            chr_info.autosome_ct = autosome_ct;
            // assumes first four codes are x, y, xy, mt
            for (uint32_t xymt_idx = 0; xymt_idx < 4; ++xymt_idx) {
              chr_info.xymt_codes[xymt_idx] = autosome_ct + 1 + xymt_idx;
            }
            for (uint32_t xymt_idx = 4; xymt_idx < kChrOffsetCt; ++xymt_idx) {
              chr_info.xymt_codes[xymt_idx] = UINT32_MAXM1;
            }
            SetBit(autosome_ct + 1, chr_info.haploid_mask);
            SetBit(autosome_ct + 2, chr_info.haploid_mask);
            for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
              cur_modif = argvk[arg_idx + param_idx];
              if (!strcmp(cur_modif, "no-x")) {
                chr_info.xymt_codes[0] = UINT32_MAXM1;
                ClearBit(autosome_ct + 1, chr_info.haploid_mask);
              } else if (!strcmp(cur_modif, "no-y")) {
                chr_info.xymt_codes[1] = UINT32_MAXM1;
                ClearBit(autosome_ct + 2, chr_info.haploid_mask);
              } else if (!strcmp(cur_modif, "no-xy")) {
                chr_info.xymt_codes[2] = UINT32_MAXM1;
              } else if (!strcmp(cur_modif, "no-mt")) {
                chr_info.xymt_codes[3] = UINT32_MAXM1;
                ClearBit(autosome_ct + 4, chr_info.haploid_mask);
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --chr-set parameter '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "hr-override")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (!strcmp(cur_modif, "file")) {
              pc.misc_flags |= kfMiscChrOverrideFile;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --chr-override parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.misc_flags |= kfMiscChrOverrideCmdline;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ovar-quantile-normalize")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &pc.covar_quantnorm_flattened);
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.pheno_transform_flags |= kfPhenoTransformQuantnormCovar;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ovar-variance-standardize")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &pc.vstd_flattened);
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.pheno_transform_flags |= kfPhenoTransformVstdCovar;
          pc.dependency_flags |= kfFilterPsamReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'd':
        if (strequal_k_unsafe(flagname_p2, "ouble-id")) {
          if (const_fid) {
            logerrputs("Error: --double-id cannot be used with --const-fid.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          import_flags |= kfImportDoubleId;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ebug")) {
          g_debug_on = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ata")) {
          if (load_params || (xload & (~kfXloadOxBgen))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t is_gzs = 0;
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "ref-first")) {
              oxford_import_flags |= kfOxfordImportRefFirst;
            } else if (!strcmp(cur_modif, "ref-second")) {
              oxford_import_flags |= kfOxfordImportRefSecond;
            } else if (!strcmp(cur_modif, "gzs")) {
              if (xload & kfXloadOxBgen) {
                // may as well permit e.g. --data ref-first + --bgen
                logerrputs("Error: --data 'gzs' modifier cannot be used with .bgen input.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              is_gzs = 1;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --data parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const char* fname_prefix = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname_prefix);
          if (slen > (kPglFnamesize - 9)) {
            logerrputs("Error: --data parameter too long.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (!(xload & kfXloadOxBgen)) {
            // allow --bgen to override this
            char* genname_end = memcpya(pgenname, fname_prefix, slen);
            genname_end = Stpcpy(genname_end, ".gen");
            if (is_gzs) {
              snprintf(genname_end, 5, ".zst");
            }
            xload |= kfXloadOxGen;
          }
          snprintf(memcpya(psamname, fname_prefix, slen), 9, ".sample");
          xload |= kfXloadOxSample;
        } else if (strequal_k_unsafe(flagname_p2, "osage-erase-threshold")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dosage_erase_frac;
          if ((!ScanadvDouble(cur_modif, &dosage_erase_frac)) || (dosage_erase_frac < 0.0) || (dosage_erase_frac >= (0.5 - kSmallEpsilon))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dosage-erase-threshold parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.dosage_erase_thresh = S_CAST(int32_t, dosage_erase_frac * ((1 + kSmallEpsilon) * kDosageMid));
        } else if (strequal_k_unsafe(flagname_p2, "ummy")) {
          if (load_params || xload) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 8)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (ScanPosintDefcap(argvk[arg_idx + 1], &gendummy_info.sample_ct)) {
            logerrputs("Error: Invalid --dummy sample count.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (ScanPosintDefcap(argvk[arg_idx + 2], &gendummy_info.variant_ct)) {
            logerrputs("Error: Invalid --dummy SNP count.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          uint32_t extra_numeric_param_ct = 0;
          for (uint32_t param_idx = 3; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (MatchUpperKLen(cur_modif, "ACGT", cur_modif_slen)) {
              gendummy_info.flags |= kfGenDummyAcgt;
            } else if (strequal_k(cur_modif, "1234", cur_modif_slen)) {
              gendummy_info.flags |= kfGenDummy1234;
            } else if (strequal_k(cur_modif, "12", cur_modif_slen)) {
              gendummy_info.flags |= kfGenDummy12;
            } else if (StrStartsWith(cur_modif, "pheno-ct=", cur_modif_slen)) {
              const char* pheno_ct_start = &(cur_modif[strlen("pheno-ct=")]);
              if (ScanUintCapped(pheno_ct_start, kMaxPhenoCt, &gendummy_info.pheno_ct)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dummy pheno-ct= parameter '%s'.\n", pheno_ct_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "scalar-pheno", cur_modif_slen)) {
              gendummy_info.flags |= kfGenDummyScalarPheno;
            } else if (StrStartsWith(cur_modif, "dosage-freq=", cur_modif_slen)) {
              const char* dosage_freq_start = &(cur_modif[strlen("dosage-freq=")]);
              double dxx;
              if ((!ScanadvDouble(dosage_freq_start, &dxx)) || (dxx < 0.0) || (dxx > 1.0)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dummy dosage-freq= parameter '%s'.\n", dosage_freq_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              gendummy_info.dosage_freq = dxx;
            } else {
              double dxx;
              if ((extra_numeric_param_ct == 2) || (!ScanadvDouble(cur_modif, &dxx)) || (dxx < 0.0) || (dxx > 1.0)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dummy parameter '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (!extra_numeric_param_ct) {
                gendummy_info.geno_mfreq = dxx;
              } else {
                gendummy_info.pheno_mfreq = dxx;
              }
              ++extra_numeric_param_ct;
            }
          }
          const uint32_t mutually_exclusive_flags = gendummy_info.flags & (kfGenDummyAcgt | kfGenDummy1234 | kfGenDummy12);
          if (mutually_exclusive_flags & (mutually_exclusive_flags - 1)) {
            logerrputs("Error: --dummy 'acgt', '1234', and '12' modifiers are mutually exclusive.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          xload |= kfXloadGenDummy;
        } else if (strequal_k_unsafe(flagname_p2, "ummy-coding")) {
          logerrputs("Error: --dummy-coding is retired.  Use --split-cat-pheno instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "og")) {
          if (chr_info.chrset_source) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.autosome_ct = 38;
          chr_info.xymt_codes[0] = 39;
          chr_info.xymt_codes[1] = 40;
          chr_info.xymt_codes[2] = 41;
          chr_info.xymt_codes[3] = 42;
          chr_info.xymt_codes[4] = UINT32_MAXM1;
          chr_info.xymt_codes[5] = UINT32_MAXM1;
#ifdef __LP64__
          chr_info.haploid_mask[0] = 0x58000000000LLU;
#else
          chr_info.haploid_mask[0] = 0;
          chr_info.haploid_mask[1] = 0x580;
#endif
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'e':
        if (strequal_k_unsafe(flagname_p2, "xtract")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          uint32_t is_ibed = 0;
          if (cur_modif_slen == 5) {
            if (!memcmp(cur_modif, "ibed0", 5)) {
              if (param_ct == 1) {
                logerrputs("Error: '--extract ibed0' requires at least one filename.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.filter_flags |= kfFilterExtractIbed0 | kfFilterNoSplitChr;
              is_ibed = 1;
            } else if ((!memcmp(cur_modif, "ibed1", 5)) || (!memcmp(cur_modif, "range", 5))) {
              if (param_ct == 1) {
                logerrputs("Error: '--extract ibed1' requires at least one filename.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.filter_flags |= kfFilterExtractIbed1 | kfFilterNoSplitChr;
              is_ibed = 1;
            }
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1 + is_ibed]), param_ct - is_ibed, kPglFnamesize, &pc.extract_fnames);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xclude")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          uint32_t is_ibed = 0;
          if (cur_modif_slen == 5) {
            if (!memcmp(cur_modif, "ibed0", 5)) {
              if (param_ct == 1) {
                logerrputs("Error: '--exclude ibed0' requires at least one filename.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.filter_flags |= kfFilterExcludeIbed0 | kfFilterNoSplitChr;
              is_ibed = 1;
            } else if ((!memcmp(cur_modif, "ibed1", 5)) || (!memcmp(cur_modif, "range", 5))) {
              if (param_ct == 1) {
                logerrputs("Error: '--exclude ibed1' requires at least one filename.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.filter_flags |= kfFilterExcludeIbed1 | kfFilterNoSplitChr;
              is_ibed = 1;
            }
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1 + is_ibed]), param_ct - is_ibed, kPglFnamesize, &pc.exclude_fnames);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xport") || strequal_k_unsafe(flagname_p2, "xport ped")) {
          // todo: determine actual limit
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 50)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t format_param_idxs = 0;
          if (!flagname_p2[5]) {
            GetExportfTargets(&(argvk[arg_idx]), param_ct, &pc.exportf_flags, &pc.exportf_id_paste, &format_param_idxs);
            if (!format_param_idxs) {
              logerrputs("Error: --export requires at least one output format.  (Did you forget 'ped' or\n'vcf'?)\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            logputsb();
          } else {
            pc.exportf_flags = kfExportfPed;
          }
          // can't have e.g. bgen-1.1 and bgen-1.2 simultaneously, since they
          // have the same extension and different content.
          const uint64_t bgen_flags = S_CAST(uint64_t, pc.exportf_flags & (kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13));
          if (bgen_flags & (bgen_flags - 1)) {
            logerrputs("Error: Multiple --export bgen versions.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if ((pc.exportf_flags & (kfExportfHaps | kfExportfHapsLegend)) == (kfExportfHaps | kfExportfHapsLegend)) {
            logerrputs("Error: 'haps' and 'hapslegend' formats cannot be exported simultaneously.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if ((pc.exportf_flags & (kfExportfA | kfExportfAD)) == (kfExportfA | kfExportfAD)) {
            logerrputs("Error: 'A' and 'AD' formats cannot be exported simultaneously.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            // could use AdvBoundedTo0Bit()...
            if ((format_param_idxs >> param_idx) & 1) {
              continue;
            }
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (StrStartsWith(cur_modif, "id-paste=", cur_modif_slen)) {
              if (!(pc.exportf_flags & (kfExportfVcf | kfExportfBgen12 | kfExportfBgen13))) {
                // todo: bcf
                logerrputs("Error: The 'id-paste' modifier only applies to --export's vcf, bgen-1.2, and\nbgen-1.3 output formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (pc.exportf_id_paste) {
                logerrputs("Error: Multiple --export id-paste= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("id-paste=")]), "maybefid\0fid\0iid\0maybesid\0sid\0", "export", kfIdpasteMaybefid, kfIdpasteDefault, 1, &pc.exportf_id_paste);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "id-delim=", cur_modif_slen)) {
              if (!(pc.exportf_flags & (kfExportfVcf | kfExportfBgen12 | kfExportfBgen13))) {
                logerrputs("Error: The 'id-delim' modifier only applies to --export's vcf, bgen-1.2, and\nbgen-1.3 output formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (pc.exportf_id_delim) {
                logerrputs("Error: Multiple --export id-delim= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              pc.exportf_id_delim = ExtractCharParam(&(cur_modif[strlen("id-delim=")]));
              if (!pc.exportf_id_delim) {
                logerrputs("Error: --export id-delim= value must be a single character.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if ((ctou32(pc.exportf_id_delim) < ' ') || (pc.exportf_id_delim == '0')) {
                logerrputs("Error: --export id-delim= value cannot be tab, newline, '0', or a nonprinting\ncharacter.\n");
                goto main_ret_INVALID_CMDLINE;
              }
            } else if (StrStartsWith(cur_modif, "vcf-dosage=", cur_modif_slen)) {
              if (!(pc.exportf_flags & kfExportfVcf)) {
                logerrputs("Error: The 'vcf-dosage' modifier only applies to --export's vcf output format.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (pc.exportf_flags & (kfExportfVcfDosageGp | kfExportfVcfDosageDs)) {
                logerrputs("Error: Multiple --export vcf-dosage= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* vcf_dosage_start = &(cur_modif[strlen("vcf-dosage=")]);
              if (!strcmp(vcf_dosage_start, "GP")) {
                pc.exportf_flags |= kfExportfVcfDosageGp;
              } else if (!strcmp(vcf_dosage_start, "DS")) {
                pc.exportf_flags |= kfExportfVcfDosageDs;
              } else if (!strcmp(vcf_dosage_start, "DS-force")) {
                pc.exportf_flags |= kfExportfVcfDosageDs | kfExportfVcfDosageForce;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --export vcf-dosage= parameter '%s'.\n", vcf_dosage_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (StrStartsWith(cur_modif, "bits=", cur_modif_slen)) {
              if (!(pc.exportf_flags & (kfExportfBgen12 | kfExportfBgen13))) {
                logerrputs("Error: The 'bits' modifier only applies to --export's bgen-1.2 and bgen-1.3\noutput formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (pc.exportf_bits) {
                logerrputs("Error: Multiple --export bits= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* bits_start = &(cur_modif[strlen("bits=")]);
              if (ScanPosintCapped(bits_start, 32, &pc.exportf_bits)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --export bits= parameter '%s'.\n", bits_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "include-alt", cur_modif_slen)) {
              if (!(pc.exportf_flags & (kfExportfA | kfExportfAD))) {
                logerrputs("Error: The 'include-alt' modifier only applies to --export's A and AD output\nformats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.exportf_flags |= kfExportfIncludeAlt;
            } else if (strequal_k(cur_modif, "omit-nonmale-y", cur_modif_slen)) {
              if (!(pc.exportf_flags & (kfExportfList | kfExportfRlist))) {
                logerrputs("Error: The 'omit-nonmale-y' modifier only applies to --export's list and rlist\noutput formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.exportf_flags |= kfExportfOmitNonmaleY;
            } else if (strequal_k(cur_modif, "01", cur_modif_slen) || strequal_k(cur_modif, "12", cur_modif_slen)) {
              if (pc.exportf_flags & (kfExportfA | kfExportfAD)) {
                snprintf(g_logbuf, kLogbufSize, "Error: The '%s' modifier does not apply to --export's A and AD output formats.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              if (pc.exportf_flags & kfExportfVcf) {
                logerrputs("Error: '01'/'12' cannot be used with --export's vcf output format.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (cur_modif[0] == '0') {
                if (pc.exportf_flags & kfExportf12) {
                  logerrputs("Error: --export '01' and '12' cannot be used together.\n");
                  goto main_ret_INVALID_CMDLINE;
                }
                pc.exportf_flags |= kfExportf01;
              } else {
                if (pc.exportf_flags & kfExportf01) {
                  logerrputs("Error: --export '01' and '12' cannot be used together.\n");
                  goto main_ret_INVALID_CMDLINE;
                }
                pc.exportf_flags |= kfExportf12;
              }
            } else if (strequal_k(cur_modif, "bgz", cur_modif_slen)) {
              if (!(pc.exportf_flags & (kfExportfOxGen | kfExportfVcf))) {
                logerrputs("Error: The 'bgz' modifier only applies to --export's oxford and vcf output\nformats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.exportf_flags |= kfExportfBgz;
            } else if (strequal_k(cur_modif, "spaces", cur_modif_slen)) {
              pc.exportf_flags |= kfExportfSpaces;
            } else if (strequal_k(cur_modif, "ref-first", cur_modif_slen)) {
              pc.exportf_flags |= kfExportfRefFirst;
            } else if (strequal_k(cur_modif, "gen-gz", cur_modif_slen)) {
              logerrputs("Error: 'gen-gz' modifier retired.  Use '--export oxford bgz' instead.\n");
              goto main_ret_INVALID_CMDLINE_WWA;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --export parameter '%s'.%s\n", cur_modif, ((param_idx == param_ct) && (!outname_end))? " (Did you forget '--out'?)" : "");
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (pc.exportf_flags & (kfExportfVcf | kfExportfBgen12 | kfExportfBgen13)) {
            if (!pc.exportf_id_paste) {
              pc.exportf_id_paste = kfIdpasteDefault;
            }
          }
          pc.command_flags1 |= kfCommand1Exportf;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "xclude-snp")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_exclude_snp);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xclude-snps")) {
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.exclude_snps_range_list);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xtract-if-info")) {
          reterr = ValidateAndAllocCmpExpr(&(argvk[arg_idx + 1]), argvk[arg_idx], param_ct, &pc.extract_if_info_expr);
          if (reterr) {
            goto main_ret_1;
          }
          // validator doesn't currently check for ';'.  also theoretically
          // possible for '=' to be in key
          if (strchr(pc.extract_if_info_expr.pheno_name, ';') || strchr(pc.extract_if_info_expr.pheno_name, '=')) {
            logerrputs("Error: Invalid --extract-if-info expression.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          // LoadPvar() currently checks value string if nonnumeric
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xclude-if-info")) {
          reterr = ValidateAndAllocCmpExpr(&(argvk[arg_idx + 1]), argvk[arg_idx], param_ct, &pc.exclude_if_info_expr);
          if (reterr) {
            goto main_ret_1;
          }
          // validator doesn't currently check for ';'.  also theoretically
          // possible for '=' to be in key
          if (strchr(pc.exclude_if_info_expr.pheno_name, ';') || strchr(pc.exclude_if_info_expr.pheno_name, '=')) {
            logerrputs("Error: Invalid --exclude-if-info expression.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          // LoadPvar() currently checks value string if nonnumeric
          pc.filter_flags |= kfFilterPvarReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'f':
        if (strequal_k_unsafe(flagname_p2, "req")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 5)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t bins_only = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.freq_rpt_flags |= kfAlleleFreqZs;
            } else if (strequal_k(cur_modif, "counts", cur_modif_slen)) {
              pc.freq_rpt_flags |= kfAlleleFreqCounts;
            } else if (strequal_k(cur_modif, "case-control", cur_modif_slen)) {
              logerrputs("Error: --freq 'case-control' modifier has been retired.  Use\n--keep-if/--remove-if in conjunction with Unix text-processing utilities\ninstead.\n");
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (pc.freq_rpt_flags & kfAlleleFreqColAll) {
                logerrputs("Error: Multiple --freq cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0reffreq\0alt1freq\0altfreq\0freq\0eq\0eqz\0alteq\0alteqz\0numeq\0altnumeq\0machr2\0nobs\0", "freq", kfAlleleFreqColChrom, kfAlleleFreqColDefault, 1, &pc.freq_rpt_flags);
              if (reterr) {
                goto main_ret_1;
              }
              const uint32_t mutually_exclusive_cols = pc.freq_rpt_flags & kfAlleleFreqColMutex;
              if (mutually_exclusive_cols & (mutually_exclusive_cols - 1)) {
                logerrputs("Error: --freq's altfreq, freq, eq, eqz, alteq, alteqz, numeq, and altnumeq\ncolumns are mutually exclusive.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else if (strequal_k(cur_modif, "bins-only", cur_modif_slen)) {
              bins_only = 1;
            } else if (StrStartsWith(cur_modif, "refbins=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "refbins-file=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "alt1bins=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "alt1bins-file=", cur_modif_slen)) {
              const uint32_t is_alt1 = (cur_modif[0] == 'a');
              char** binstr_ptr = is_alt1? (&pc.freq_alt1_binstr) : (&pc.freq_ref_binstr);
              if (*binstr_ptr) {
                logerrprintf("Error: Multiple --freq %sbins{-file}= modifiers.\n", is_alt1? "alt1" : "ref");
                goto main_ret_INVALID_CMDLINE;
              }
              if (cur_modif[7 + is_alt1] == '=') {
                // "refbins=", "alt1bins="
                reterr = CmdlineAllocString(&(cur_modif[8 + is_alt1]), is_alt1? "--freq alt1bins=" : "--freq refbins=", 0x7fffffff, binstr_ptr);
              } else {
                // "refbins-file=", "alt1bins-file="
                pc.freq_rpt_flags |= is_alt1? kfAlleleFreqBinsAlt1Fname : kfAlleleFreqBinsRefFname;
                reterr = AllocFname(&(cur_modif[13 + is_alt1]), is_alt1? "freq alt1bins-file=" : "freq refbins-file=", 0, binstr_ptr);
              }
              if (reterr) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --freq parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (bins_only) {
            if ((!pc.freq_ref_binstr) && (!pc.freq_alt1_binstr)) {
              logerrputs("Error: --freq 'bins-only' must be used with 'refbins{-file}=' and/or\n'alt1bins{-file}='.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (pc.freq_rpt_flags & (kfAlleleFreqZs | kfAlleleFreqColAll)) {
              logerrputs("Error: --freq 'bins-only' cannot be used with 'zs' or 'cols=' (which only\naffect the main report).\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.freq_rpt_flags |= kfAlleleFreqBinsOnly;
          }
          if (!(pc.freq_rpt_flags & kfAlleleFreqColAll)) {
            pc.freq_rpt_flags |= kfAlleleFreqColDefault;
          }
          pc.command_flags1 |= kfCommand1AlleleFreq;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "rom")) {
          if (chr_info.is_include_stack) {
            logerrputs("Error: --from/--to cannot be used with --autosome{-par} or --{not-}chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_from);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "rom-bp") ||
                   strequal_k_unsafe(flagname_p2, "rom-kb") ||
                   strequal_k_unsafe(flagname_p2, "rom-mb")) {
          if (!CmdlineSingleChr(&chr_info, pc.misc_flags)) {
            logerrputs("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb must be used with --chr, and only\none chromosome.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (pc.from_bp != -1) {
            logerrputs("Error: Multiple --from-bp/-kb/-mb values.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          // permit negative numbers, to simplify shell script windowing logic
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (!ScanadvDouble(cur_modif, &dxx)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --from-bp/-kb/-mb parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          const char unit_char = flagname_p2[4];
          if (unit_char == 'k') {
            dxx *= 1000;
          } else if (unit_char == 'm') {
            dxx *= 1000000;
          }
          if (dxx <= 0.0) {
            pc.from_bp = 0;
          } else {
            // logical to round up rather than down here (this is actually a
            // change from v1.9)
            // don't use ceil() since e.g. ceil(0.001015 * 1000000) is 1016
            if (dxx > 2147483646.0) {
              logerrprintf("Error: --from-bp/-kb/-mb parameter '%s' too large.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.from_bp = 1 + S_CAST(int32_t, dxx * (1 - kSmallEpsilon));
          }
          pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "orce-intersect")) {
          permit_multiple_inclusion_filters = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "amily")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (IsReservedPhenoName(cur_modif, cur_modif_slen)) {
              snprintf(g_logbuf, kLogbufSize, "Error: '%s' cannot be used as a categorical phenotype name.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_2A;
            }
            reterr = CmdlineAllocString(cur_modif, argvk[arg_idx], kMaxIdSlen, &pc.catpheno_name);
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.misc_flags |= kfMiscCatPhenoFamily;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "amily-missing-catname")) {
          if (!(pc.misc_flags & kfMiscCatPhenoFamily)) {
            logerrputs("Error: --family-missing-catname must be used with --family.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.family_missing_catname);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ilter-cases") || strequal_k_unsafe(flagname_p2, "ilter-controls")) {
          logerrputs("Error: --filter-cases and --filter-controls have been retired.  Use\n--keep-if/--remove-if instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "rqx") || strequal_k_unsafe(flagname_p2, "reqx")) {
          logerrputs("Error: --freqx has been retired.  Use --geno-counts instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'g':
        if (strequal_k_unsafe(flagname_p2, "eno")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t geno_thresh_present = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "dosage")) {
              pc.misc_flags |= kfMiscGenoDosage;
            } else if (!strcmp(cur_modif, "hh-missing")) {
              pc.misc_flags |= kfMiscGenoHhMissing;
            } else if (geno_thresh_present) {
              logerrputs("Error: Invalid --geno parameter sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (!ScanadvDouble(cur_modif, &pc.geno_thresh)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --geno parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else if ((pc.geno_thresh < 0.0) || (pc.geno_thresh > 1.0)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --geno parameter '%s' (must be in [0, 1]).\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else {
              geno_thresh_present = 1;
            }
          }
          if (!geno_thresh_present) {
            pc.geno_thresh = 0.1;
          }
          if (pc.geno_thresh < 1.0) {
            pc.filter_flags |= kfFilterPvarReq;
            pc.dependency_flags = kfFilterAllReq | kfFilterNoSplitChr;
          }
        } else if (strequal_k_unsafe(flagname_p2, "eno-counts")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.geno_counts_flags |= kfGenoCountsZs;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (pc.geno_counts_flags & kfGenoCountsColAll) {
                logerrputs("Error: Multiple --geno-counts cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0homref\0refalt1\0refalt\0homalt1\0altxy\0xy\0hapref\0hapalt1\0hapalt\0hap\0numeq\0missing\0nobs\0", "geno-counts", kfGenoCountsColChrom, kfGenoCountsColDefault, 1, &pc.geno_counts_flags);
              if (reterr) {
                goto main_ret_1;
              }
              if ((pc.geno_counts_flags & kfGenoCountsColPairex) == kfGenoCountsColPairex) {
                logerrputs("Error: --geno-counts's hapaltx and hapx columns are mutually exclusive.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              const uint32_t mutually_exclusive_cols = pc.geno_counts_flags & kfGenoCountsColMutex;
              if (mutually_exclusive_cols & (mutually_exclusive_cols - 1)) {
                logerrputs("Error: --geno-counts's altxy, xy, and numeq columns are mutually exclusive.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --geno-counts parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.geno_counts_flags & kfGenoCountsColAll)) {
            pc.geno_counts_flags |= kfGenoCountsColDefault;
          }
          pc.command_flags1 |= kfCommand1GenoCounts;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "lm")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 17)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmZs;
            } else if (strequal_k(cur_modif, "a0-ref", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmA0Ref;
            } else if (strequal_k(cur_modif, "sex", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmSex;
            } else if (strequal_k(cur_modif, "no-x-sex", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmNoXSex;
            } else if (strequal_k(cur_modif, "log10", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmLog10;
            } else if (strequal_k(cur_modif, "genotypic", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmGenotypic;
            } else if (strequal_k(cur_modif, "hethom", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmHethom;
            } else if (strequal_k(cur_modif, "dominant", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmDominant;
            } else if (strequal_k(cur_modif, "recessive", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmRecessive;
            } else if (strequal_k(cur_modif, "interaction", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmInteraction;
            } else if (strequal_k(cur_modif, "hide-covar", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmHideCovar;
            } else if (strequal_k(cur_modif, "intercept", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmIntercept;
            } else if (strequal_k(cur_modif, "firth-fallback", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmFirthFallback;
            } else if (strequal_k(cur_modif, "firth", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmFirth;
            } else if (strequal_k(cur_modif, "standard-beta", cur_modif_slen)) {
              logerrputs("Error: --glm 'standard-beta' modifier has been retired.  Use\n--{covar-}variance-standardize instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (strequal_k(cur_modif, "perm", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmPerm;
            } else if (strequal_k(cur_modif, "perm-count", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmPermCount;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (pc.glm_info.cols) {
                logerrputs("Error: Multiple --glm cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0a0\0a1count\0totallele\0a1countcc\0totallelecc\0gcountcc\0a1freq\0a1freqcc\0machr2\0firth\0test\0nobs\0beta\0orbeta\0se\0ci\0tz\0p\0", flagname_p, kfGlmColChrom, kfGlmColDefault, 1, &pc.glm_info.cols);
              if (reterr) {
                goto main_ret_1;
              }
              if ((!(pc.glm_info.cols & (kfGlmColBeta | kfGlmColOrbeta))) && ((pc.glm_info.cols & kfGlmColSe) || ((pc.glm_info.cols & kfGlmColCi) && (pc.ci_size != 0)))) {
                logerrputs("Error: --glm's 'se' and 'ci' columns require beta/orbeta to be included.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else if (StrStartsWith0(cur_modif, "mperm", cur_modif_slen)) {
              if ((cur_modif_slen < 7) || (cur_modif[5] != '=')) {
                logerrputs("Error: Improper --glm mperm syntax.  (Use --glm mperm=[value]'.)\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (ScanPosintDefcap(cur_modif, &pc.glm_info.mperm_ct)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --glm mperm parameter '%s'.\n", &(cur_modif[6]));
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (StrStartsWith(cur_modif, "local-covar=", cur_modif_slen)) {
              if (pc.glm_local_covar_fname) {
                logerrputs("Error: Multiple --glm local-covar= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = AllocFname(&(cur_modif[strlen("local-covar=")]), "glm local-covar=", 0, &pc.glm_local_covar_fname);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "local-bim=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "local-pvar=", cur_modif_slen)) {
              if (pc.glm_local_pvar_fname) {
                logerrputs("Error: Multiple --glm local-pvar= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const uint32_t is_pvar = (cur_modif[6] == 'p');
              reterr = AllocFname(&(cur_modif[10 + is_pvar]), "glm local-pvar=", 0, &pc.glm_local_pvar_fname);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "local-fam=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "local-psam=", cur_modif_slen)) {
              if (pc.glm_local_psam_fname) {
                logerrputs("Error: Multiple --glm local-psam= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const uint32_t is_psam = (cur_modif[6] == 'p');
              reterr = AllocFname(&(cur_modif[10 + is_psam]), "glm local-psam=", 0, &pc.glm_local_psam_fname);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (strequal_k(cur_modif, "local-omit-last", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmLocalOmitLast;
            } else if (StrStartsWith(cur_modif, "local-cats=", cur_modif_slen)) {
              if (pc.glm_info.local_cat_ct) {
                logerrputs("Error: Multiple --glm local-cats= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              // bugfix (7 Nov 2017): forgot to offset by strlen("local-cats=")
              if (ScanPosintCapped(&(cur_modif[strlen("local-cats=")]), 4095, &pc.glm_info.local_cat_ct) || (pc.glm_info.local_cat_ct == 1)) {
                logerrputs("Error: Invalid --glm local-cats= category count (must be in [2, 4095]).\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --glm parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!pc.glm_info.cols) {
            pc.glm_info.cols = kfGlmColDefault;
          }
          if ((pc.glm_info.flags & (kfGlmSex | kfGlmNoXSex)) == (kfGlmSex | kfGlmNoXSex)) {
            logerrputs("Error: Conflicting --glm parameters.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if ((pc.glm_info.flags & kfGlmPerm) && pc.glm_info.mperm_ct) {
            logerrputs("Error: --glm 'perm' and 'mperm=' cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          uint32_t alternate_genotype_col_flags = S_CAST(uint32_t, pc.glm_info.flags & (kfGlmGenotypic | kfGlmHethom | kfGlmDominant | kfGlmRecessive));
          if (alternate_genotype_col_flags) {
            pc.xchr_model = 0;
            if (alternate_genotype_col_flags & (alternate_genotype_col_flags - 1)) {
              logerrputs("Error: Conflicting --glm parameters.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          if ((pc.glm_info.flags & kfGlmIntercept) && (!(pc.glm_info.cols & kfGlmColTest))) {
            logerrputs("Error: --glm 'intercept' modifier cannot be used with an omitted 'test' column.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (!pc.glm_local_covar_fname) {
            if (pc.glm_local_pvar_fname || pc.glm_local_psam_fname) {
              logerrputs("Error: Either all three --glm local-covar filenames must be specified, or none\nof them.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (pc.glm_info.flags & kfGlmLocalOmitLast) {
              logerrputs("Error: --glm 'local-omit-last' must be used with 'local-covar='.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (pc.glm_info.local_cat_ct) {
              logerrputs("Error: --glm 'local-cats=' must be used with 'local-covar='.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else {
            if ((!pc.glm_local_pvar_fname) || (!pc.glm_local_psam_fname)) {
              logerrputs("Error: Either all three --glm local-covar filenames must be specified, or none\nof them.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          pc.command_flags1 |= kfCommand1Glm;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "en")) {
          if (load_params || (xload & (~kfXloadOxSample))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 2) {
            const char* cur_modif = argvk[arg_idx + 2];
            if (!strcmp(cur_modif, "ref-first")) {
              oxford_import_flags |= kfOxfordImportRefFirst;
            } else if (!strcmp(cur_modif, "ref-second")) {
              oxford_import_flags |= kfOxfordImportRefSecond;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --gen parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const char* cur_fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_fname);
          if (slen > kPglFnamesize - 1) {
            logerrputs("Error: --gen filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_fname, slen + 1);
          xload |= kfXloadOxGen;
        } else if (strequal_k_unsafe(flagname_p2, "enotyping-rate")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (strcmp("dosage", cur_modif)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --genotyping-rate parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            pc.misc_flags |= kfMiscGenotypingRateDosage;
          }
          pc.command_flags1 |= kfCommand1GenotypingRate;
          pc.dependency_flags |= kfFilterAllReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'h':
        if (strequal_k_unsafe(flagname_p2, "ardy")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 3)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.hardy_flags |= kfHardyZs;
            } else if (strequal_k(cur_modif, "midp", cur_modif_slen)) {
              pc.hardy_flags |= kfHardyMidp;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (pc.hardy_flags & kfHardyColAll) {
                logerrputs("Error: Multiple --hardy cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0gcounts\0gcount1col\0hetfreq\0sexaf\0femalep\0p\0", "freq", kfHardyColChrom, kfHardyColDefault, 1, &pc.hardy_flags);
              if (reterr) {
                goto main_ret_1;
              }
              if ((pc.hardy_flags & (kfHardyColGcounts | kfHardyColGcount1col)) == (kfHardyColGcounts | kfHardyColGcount1col)) {
                logerrputs("Error: --hardy's gcounts and gcounts1col column sets are mutually exclusive.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --hardy parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.hardy_flags & kfHardyColAll)) {
            pc.hardy_flags |= kfHardyColDefault;
          }
          pc.command_flags1 |= kfCommand1Hardy;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "we")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "midp")) {
              pc.misc_flags |= kfMiscHweMidp;
            } else if (!strcmp(cur_modif, "keep-fewhet")) {
              pc.misc_flags |= kfMiscHweKeepFewhet;
            } else {
              if ((pc.hwe_thresh != 1.0) || (!ScanadvDouble(cur_modif, &pc.hwe_thresh))) {
                logerrputs("Error: Invalid --hwe parameter sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if ((pc.hwe_thresh < 0.0) || (pc.hwe_thresh >= 1.0)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --hwe threshold '%s' (must be in [0, 1)).\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            }
          }
          if (pc.hwe_thresh == 1.0) {
            logerrputs("Error: --hwe requires a p-value threshold.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if ((pc.misc_flags & kfMiscHweMidp) && (pc.hwe_thresh >= 0.5)) {
            logerrputs("Error: --hwe threshold must be smaller than 0.5 when using mid-p adjustment.\n");
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "ard-call-threshold")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double hard_call_frac;
          if ((!ScanadvDouble(cur_modif, &hard_call_frac)) || (hard_call_frac < 0.0) || (hard_call_frac >= (0.5 - kSmallEpsilon))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --hard-call-threshold parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.hard_call_thresh = S_CAST(int32_t, hard_call_frac * ((1 + kSmallEpsilon) * kDosageMid));
        } else if (strequal_k_unsafe(flagname_p2, "orse")) {
          if (chr_info.chrset_source) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.autosome_ct = 31;
          chr_info.xymt_codes[0] = 32;
          chr_info.xymt_codes[1] = 33;
          chr_info.xymt_codes[2] = UINT32_MAXM1;
          chr_info.xymt_codes[3] = UINT32_MAXM1;
          chr_info.xymt_codes[4] = UINT32_MAXM1;
          chr_info.xymt_codes[5] = UINT32_MAXM1;
#ifdef __LP64__
          chr_info.haploid_mask[0] = 0x300000000LLU;
#else
          chr_info.haploid_mask[0] = 0;
          chr_info.haploid_mask[1] = 3;
#endif
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "uman")) {
          if (chr_info.chrset_source) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.chrset_source = kChrsetSourceCmdline;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "aps")) {
          if (load_params || xload) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 2) {
            const char* cur_modif = argvk[arg_idx + 2];
            if (!strcmp(cur_modif, "ref-first")) {
              oxford_import_flags |= kfOxfordImportRefFirst;
            } else if (!strcmp(cur_modif, "ref-second")) {
              oxford_import_flags |= kfOxfordImportRefSecond;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --haps parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const char* cur_fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_fname);
          if (slen > kPglFnamesize - 1) {
            logerrputs("Error: --haps filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_fname, slen + 1);
          xload |= kfXloadOxHaps;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'i':
        if (strequal_k_unsafe(flagname_p2, "ndiv-sort")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* indiv_sort_mode_str = argvk[arg_idx + 1];
          const char first_char_upcase_match = indiv_sort_mode_str[0] & 0xdf;
          const uint32_t is_short_name = (indiv_sort_mode_str[1] == '\0');
          if ((is_short_name && (indiv_sort_mode_str[0] == '0')) || (!strcmp(indiv_sort_mode_str, "none")))  {
            pc.sample_sort_flags = kfSortNone;
          } else if ((is_short_name && (first_char_upcase_match == 'N')) || (!strcmp(indiv_sort_mode_str, "natural"))) {
            pc.sample_sort_flags = kfSortNatural;
          } else if ((is_short_name && (first_char_upcase_match == 'A')) || (!strcmp(indiv_sort_mode_str, "ascii"))) {
            pc.sample_sort_flags = kfSortAscii;
          } else if ((is_short_name && ((indiv_sort_mode_str[0] & 0xdf) == 'F')) || (!strcmp(indiv_sort_mode_str, "file"))) {
            if (param_ct == 1) {
              snprintf(g_logbuf, kLogbufSize, "Error: Missing '--indiv-sort %s' filename.\n", indiv_sort_mode_str);
              goto main_ret_INVALID_CMDLINE_2A;
            }
            pc.sample_sort_flags = kfSortFile;
            uint32_t fname_modif_idx = 2;
            if (param_ct == 3) {
              if (CheckExtraParam(&(argvk[arg_idx]), "sid", &fname_modif_idx)) {
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.sample_sort_flags |= kfSortFileSid;
            }
            reterr = AllocFname(argvk[arg_idx + fname_modif_idx], flagname_p, 0, &pc.sample_sort_fname);
            if (reterr) {
              goto main_ret_1;
            }
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a valid mode for --indiv-sort.\n", indiv_sort_mode_str);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if ((param_ct > 1) && (!(pc.sample_sort_flags & kfSortFile))) {
            snprintf(g_logbuf, kLogbufSize, "Error: '--indiv-sort %s' does not accept additional parameters.\n", indiv_sort_mode_str);
            goto main_ret_INVALID_CMDLINE_2A;
          }
        } else if (strequal_k_unsafe(flagname_p2, "d-delim")) {
          if (const_fid || (import_flags & kfImportDoubleId)) {
            logerrputs("Error: --id-delim can no longer be used with --const-fid or --double-id.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "sid")) {
              import_flags |= kfImportIdDelimSid;
            } else {
              id_delim = ExtractCharParam(cur_modif);
              if (!id_delim) {
                logerrputs("Error: --id-delim delimiter must be a single character.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (ctou32(id_delim) < ' ') {
                logerrputs("Error: --id-delim parameter cannot be tab, newline, or a nonprinting character.\n");
                goto main_ret_INVALID_CMDLINE;
              }
            }
          }
          if (!id_delim) {
            id_delim = '_';
          }
        } else if (strequal_k_unsafe(flagname_p2, "ndep-pairwise") || strequal_k_unsafe(flagname_p2, "ndep-pairphase")) {
          if (pc.command_flags1 & kfCommand1LdPrune) {
            logerrputs("Error: Multiple LD pruning commands.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 4)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double first_paramd;
          const char* first_param_end = ScanadvDouble(cur_modif, &first_paramd);
          if ((!first_param_end) || (first_paramd < 0.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s window size '%s'.\n", flagname_p, cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          uint32_t is_kb = 0;
          uint32_t next_param_idx = 2;
          if (MatchUpperK(first_param_end, "KB") && (!first_param_end[2])) {
            is_kb = 1;
          } else if (MatchUpperK(argvk[arg_idx + 2], "KB") && (!argvk[arg_idx + 2][2])) {
            is_kb = 1;
            next_param_idx = 3;
          }
          if (is_kb) {
            pc.ld_info.prune_flags |= kfLdPruneWindowBp;
            if (first_paramd > 2147483.646) {
              pc.ld_info.prune_window_size = 2147483646;
            } else {
              pc.ld_info.prune_window_size = S_CAST(int32_t, first_paramd * 1000 * (1 + kSmallEpsilon));
              if (pc.ld_info.prune_window_size < 2) {
                snprintf(g_logbuf, kLogbufSize, "Error: --%s window size cannot be smaller than 2.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
            }
          } else {
            if (first_paramd > 2147483647) {
              pc.ld_info.prune_window_size = 2147483647;
            } else {
              pc.ld_info.prune_window_size = S_CAST(int32_t, first_paramd);
            }
          }
          if (next_param_idx + 2 == param_ct) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s parameter sequence.\n", flagname_p);
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (next_param_idx < param_ct) {
            // explicit step size
            cur_modif = argvk[arg_idx + next_param_idx];
            if (ScanPosintDefcap(cur_modif, &pc.ld_info.prune_window_incr)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s window-increment '%s'.\n", flagname_p, cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (!is_kb) {
              if (pc.ld_info.prune_window_incr > pc.ld_info.prune_window_size) {
                snprintf(g_logbuf, kLogbufSize, "Error: --%s window-increment cannot be larger than window size.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
            } else if (pc.ld_info.prune_window_incr != 1) {
              snprintf(g_logbuf, kLogbufSize, "Error: --%s window-increment must be 1 when window size is in\nkilobase units.\n", flagname_p);
              goto main_ret_INVALID_CMDLINE_2A;
            }
          } else {
            pc.ld_info.prune_window_incr = 1;
          }
          cur_modif = argvk[arg_idx + param_ct];
          if ((!ScanadvDouble(cur_modif, &pc.ld_info.prune_last_param)) || (pc.ld_info.prune_last_param < 0.0) || (pc.ld_info.prune_last_param >= 1.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s r^2 threshold '%s'.\n", flagname_p2, cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.command_flags1 |= kfCommand1LdPrune;
          pc.dependency_flags |= kfFilterAllReq;
          if (flagname_p2[9] == 'p') {
            pc.ld_info.prune_flags |= kfLdPrunePairphase;
          } else {
            pc.ld_info.prune_flags |= kfLdPrunePairwise;
          }
        } else if (strequal_k_unsafe(flagname_p2, "nput-missing-genotype")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          input_missing_geno_char = ExtractCharParam(cur_modif);
          if (ctou32(input_missing_geno_char) <= ' ') {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --input-missing-genotype parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "nput-missing-phenotype")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (ScanInt32(cur_modif, &pc.missing_pheno) || ((pc.missing_pheno >= 0) && (pc.missing_pheno <= 2)) || (!ScanadvDouble(cur_modif, &dxx)) || (dxx != S_CAST(double, pc.missing_pheno))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --input-missing-phenotype parameter '%s' (must be an integer in [-2147483647, -1] or [3, 2147483647]).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "mport-dosage-certainty")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if ((!ScanadvDouble(cur_modif, &import_dosage_certainty)) || (import_dosage_certainty < 0.0) || (import_dosage_certainty > 1.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage-certainty parameter '%s' (must be in [0, 1]).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          import_dosage_certainty *= 1.0 - kSmallEpsilon;
        } else if (strequal_k_unsafe(flagname_p2, "mport-dosage")) {
          if (load_params || xload) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 11)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t format_num_m1 = 4;
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "noheader", cur_modif_slen)) {
              plink1_dosage_info.flags |= kfPlink1DosageNoheader;
            } else if ((cur_modif_slen > 6) && StrStartsWithUnsafe(cur_modif, "skip") && (cur_modif[4] >= '0') && (cur_modif[4] <= '2') && (cur_modif[5] == '=')) {
              const uint32_t skip_idx = ctou32(cur_modif[4]) - 48;
              if (plink1_dosage_info.skips[skip_idx]) {
                logerrprintf("Error: Multiple --import-dosage skip%u= modifiers.\n", skip_idx);
                goto main_ret_INVALID_CMDLINE;
              }
              if (ScanUintCapped(&(cur_modif[6]), kMaxLongLine / 2, &(plink1_dosage_info.skips[skip_idx]))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage skip%u= parameter '%s'.\n", skip_idx, &(cur_modif[6]));
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "dose1", cur_modif_slen)) {
              plink1_dosage_info.flags |= kfPlink1DosageFormatSingle01;
            } else if (StrStartsWithUnsafe(cur_modif, "format=")) {
              // strequal_k() and StrStartsWith() both suboptimal here
              if (format_num_m1 != 4) {
                logerrputs("Error: Multiple --import-dosage format= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              if (cur_modif_slen == 8) {
                format_num_m1 = ctou32(cur_modif[7]) - 49;
                if (format_num_m1 >= 3) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage format= parameter '%c'.\n", cur_modif[7]);
                  goto main_ret_INVALID_CMDLINE_2A;
                }
              } else if (strequal_k(&(cur_modif[7]), "infer", cur_modif_slen - 7)) {
                format_num_m1 = 3;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage format= parameter '%s'.\n", &(cur_modif[7]));
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "ref-first", cur_modif_slen)) {
              plink1_dosage_info.flags |= kfPlink1DosageRefFirst;
            } else if (strequal_k(cur_modif, "ref-second", cur_modif_slen)) {
              plink1_dosage_info.flags |= kfPlink1DosageRefSecond;
            } else if (StrStartsWith(cur_modif, "id-delim=", cur_modif_slen)) {
              if (plink1_dosage_info.id_delim) {
                logerrputs("Error: Multiple --import-dosage id-delim= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* id_delim_str = &(cur_modif[strlen("id-delim=")]);
              char cc = ExtractCharParam(id_delim_str);
              if (!cc) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage id-delim= parameter '%s'.\n", id_delim_str);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              plink1_dosage_info.id_delim = cc;
            } else if (StrStartsWith(cur_modif, "single-chr=", cur_modif_slen)) {
              if (import_single_chr_str) {
                logerrputs("Error: Multiple --import-dosage single-chr= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* chr_code = &(cur_modif[strlen("single-chr=")]);
              if (!(pc.misc_flags & kfMiscAllowExtraChrs)) {
                if (IsI32Neg(GetChrCodeRaw(chr_code))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage single-chr= chromosome code '%s'. (Did you forget --allow-extra-chr?)\n", chr_code);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
              }
              reterr = CmdlineAllocString(chr_code, argvk[arg_idx], kMaxIdSlen, &import_single_chr_str);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "chr-col-num=", cur_modif_slen)) {
              if (plink1_dosage_info.chr_col_idx != UINT32_MAX) {
                logerrputs("Error: Multiple --import-dosage chr-col-num= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* chr_col_num_start = &(cur_modif[strlen("chr-col-num=")]);
              uint32_t uii;
              if (ScanPosintCapped(chr_col_num_start, kMaxLongLine / 2, &uii)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage chr-col-num= parameter '%s'.\n", chr_col_num_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              plink1_dosage_info.chr_col_idx = uii - 1;
            } else if (StrStartsWith(cur_modif, "pos-col-num=", cur_modif_slen)) {
              if (plink1_dosage_info.pos_col_idx != UINT32_MAX) {
                logerrputs("Error: Multiple --import-dosage pos-col-num= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* pos_col_num_start = &(cur_modif[strlen("pos-col-num=")]);
              uint32_t uii;
              if (ScanPosintCapped(pos_col_num_start, kMaxLongLine / 2, &uii)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage pos-col-num= parameter '%s'.\n", pos_col_num_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              plink1_dosage_info.pos_col_idx = uii - 1;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }

          if (!format_num_m1) {
            plink1_dosage_info.flags |= kfPlink1DosageFormatSingle;
          } else {
            if (plink1_dosage_info.flags & kfPlink1DosageFormatSingle01) {
              if (format_num_m1 != 3) {
                logerrputs("Error: --import-dosage 'dose1' modifier must be used with 'format=1'.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              // format_num_m1 = 0;
              plink1_dosage_info.flags |= kfPlink1DosageFormatSingle;
            } else {
              if (format_num_m1 == 1) {
                plink1_dosage_info.flags |= kfPlink1DosageFormatDouble;
              } else if (format_num_m1 == 2) {
                plink1_dosage_info.flags |= kfPlink1DosageFormatTriple;
              }
            }
          }
          if ((plink1_dosage_info.flags & (kfPlink1DosageRefFirst | kfPlink1DosageRefSecond)) == (kfPlink1DosageRefFirst | kfPlink1DosageRefSecond)) {
            logerrputs("Error: --import-dosage 'ref-first' and 'ref-second' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          const uint32_t id_col_idx = plink1_dosage_info.skips[0];
          const uint32_t a1_col_idx = id_col_idx + plink1_dosage_info.skips[1] + 1;
          const uint32_t data_col_idx = a1_col_idx + plink1_dosage_info.skips[2] + 2;
          const uint32_t chr_col_idx = plink1_dosage_info.chr_col_idx;
          if (chr_col_idx != UINT32_MAX) {
            if (import_single_chr_str) {
              logerrputs("Error: --import-dosage 'single-chr=' and 'chr-col-num=' modifiers cannot be\nused together.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if ((chr_col_idx == id_col_idx) || (chr_col_idx == a1_col_idx) || (chr_col_idx == a1_col_idx + 1)) {
              logerrputs("Error: --import-dosage chr-col-num= value collides with another column.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (chr_col_idx >= data_col_idx) {
              logerrputs("Error: --import-dosage chr-col-num= value too large.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const uint32_t pos_col_idx = plink1_dosage_info.pos_col_idx;
          if (pos_col_idx != UINT32_MAX) {
            if ((pos_col_idx == id_col_idx) || (pos_col_idx == a1_col_idx) || (pos_col_idx == a1_col_idx + 1) || (pos_col_idx == chr_col_idx)) {
              logerrputs("Error: --import-dosage pos-col-num= value collides with another column.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (pos_col_idx >= data_col_idx) {
              logerrputs("Error: --import-dosage pos-col-num= value too large.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_modif);
          if (slen > kPglFnamesize - 1) {
            logerrputs("Error: --import-dosage filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_modif, slen + 1);
          xload = kfXloadPlink1Dosage;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'k':
        if (strequal_k_unsafe(flagname_p2, "eep")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, kPglFnamesize, &pc.keep_fnames);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-fam")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, kPglFnamesize, &pc.keepfam_fnames);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-autoconv")) {
          import_flags |= kfImportKeepAutoconv;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "eep-females")) {
          pc.filter_flags |= kfFilterPsamReq | kfFilterExclMales | kfFilterExclNosex;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "eep-males")) {
          pc.filter_flags |= kfFilterPsamReq | kfFilterExclFemales | kfFilterExclNosex;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "eep-nosex")) {
          pc.filter_flags |= kfFilterPsamReq | kfFilterExclFemales | kfFilterExclMales;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "eep-founders")) {
          pc.filter_flags |= kfFilterPsamReq | kfFilterExclNonfounders;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "eep-nonfounders")) {
          if (pc.filter_flags & kfFilterExclNonfounders) {
            logerrputs("Error: --keep-nonfounders cannot be used with --keep-founders.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.filter_flags |= kfFilterPsamReq | kfFilterExclFounders;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ing-cutoff")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 2) {
            // .king.id, .king.bin appended
            reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 9, &king_cutoff_fprefix);
            if (reterr) {
              goto main_ret_1;
            }
          }
          const char* cur_modif = argvk[arg_idx + param_ct];
          if ((!ScanadvDouble(cur_modif, &pc.king_cutoff)) || (pc.king_cutoff < 0.0) || (pc.king_cutoff >= 0.5)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --king-cutoff parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.command_flags1 |= kfCommand1KingCutoff;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ing-table-filter")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if ((!ScanadvDouble(cur_modif, &pc.king_table_filter)) || (pc.king_table_filter > 0.5)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --king-table-filter parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ing-table-subset")) {
          if (pc.king_cutoff != -1) {
            logerrputs("Error: --king-table-subset cannot be used with --king-cutoff.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &pc.king_table_subset_fname);
          if (reterr) {
            goto main_ret_1;
          }
          if (param_ct == 2) {
            const char* cur_modif = argvk[arg_idx + 2];
            if ((!ScanadvDouble(cur_modif, &pc.king_table_subset_thresh)) || (pc.king_table_subset_thresh > 0.5)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --king-table-subset threshold '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.king_table_subset_thresh = -DBL_MAX;
          }
        } else if (strequal_k_unsafe(flagname_p2, "eep-if")) {
          reterr = ValidateAndAllocCmpExpr(&(argvk[arg_idx + 1]), argvk[arg_idx], param_ct, &pc.keep_if_expr);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-cats")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &pc.keep_cats_fname);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-cat-names")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, kMaxIdBlen, &pc.keep_cat_names_flattened);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-cat-pheno")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.keep_cat_phenoname);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "eep-fcol")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &pc.keep_fcol_fname);
          if (reterr) {
            goto main_ret_1;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 2]), param_ct - 1, kMaxIdBlen, &pc.keep_fcol_flattened);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-fcol-name")) {
          if (!pc.keep_fcol_fname) {
            logerrputs("Error: --keep-fcol-name must be used with --keep-fcol.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.keep_fcol_name);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "eep-fcol-num")) {
          if (!pc.keep_fcol_fname) {
            logerrputs("Error: --keep-fcol-num must be used with --keep-fcol.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.keep_fcol_name) {
            logerrputs("Error: --keep-fcol-num can't be used with --keep-fcol-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (ScanPosintDefcap(cur_modif, &pc.keep_fcol_num) || (pc.keep_fcol_num == 1)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --keep-fcol-num parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "eep-allele-order")) {
          logputs("Note: --keep-allele-order no longer has any effect.\n");
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'l':
        if (strequal_k_unsafe(flagname_p2, "ambda")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double lambda;
          if (!ScanadvDouble(argvk[arg_idx + 1], &lambda)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --lambda parameter '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (lambda < 1.0) {
            logputs("Note: --lambda parameter set to 1.\n");
            lambda = 1.0;
          }
          pc.adjust_info.lambda = lambda;
          adjust_file_info.base.lambda = lambda;
        } else if (strequal_k_unsafe(flagname_p2, "egend")) {
          if (load_params || (xload & (~kfXloadOxHaps))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (!xload) {
            logerrputs("Error: --legend must be used with --haps.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_fname = argvk[arg_idx + 1];
          uint32_t slen = strlen(cur_fname);
          if (slen > kPglFnamesize - 1) {
            logerrputs("Error: --legend filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pvarname, cur_fname, slen + 1);
          const char* chr_code = argvk[arg_idx + 2];
          if (!(pc.misc_flags & kfMiscAllowExtraChrs)) {
            if (IsI32Neg(GetChrCodeRaw(chr_code))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --legend chromosome code '%s'. (Did you forget --allow-extra-chr?)\n", chr_code);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          reterr = CmdlineAllocString(chr_code, argvk[arg_idx], kMaxIdSlen, &import_single_chr_str);
          if (reterr) {
            goto main_ret_1;
          }
          xload |= kfXloadOxLegend;
        } else if (strequal_k_unsafe(flagname_p2, "oop-cats")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.loop_cats_phenoname);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "d")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 4)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t uii = 0; uii < 2; ++uii) {
            reterr = CmdlineAllocString(argvk[arg_idx + uii + 1], argvk[arg_idx], kMaxIdSlen, &(pc.ld_info.ld_console_varids[uii]));
            if (reterr) {
              goto main_ret_1;
            }
          }
          for (uint32_t param_idx = 3; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "dosage")) {
              pc.ld_info.ld_console_flags |= kfLdConsoleDosage;
            } else if (!strcmp(cur_modif, "hwe-midp")) {
              pc.ld_info.ld_console_flags |= kfLdConsoleHweMidp;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --ld parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.command_flags1 |= kfCommand1Ld;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "oop-assoc")) {
          logerrputs("Error: --loop-assoc is retired.  Use --within + --split-cat-pheno instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'm':
        if (strequal_k_unsafe(flagname_p2, "emory")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t mb_modif_idx = 1;
          if (param_ct == 2) {
            if (CheckExtraParam(&(argvk[arg_idx]), "require", &mb_modif_idx)) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            memory_require = 1;
          }
          const char* mb_modif = argvk[arg_idx + mb_modif_idx];
          if (ScanPosintptr(mb_modif, R_CAST(uintptr_t*, &malloc_size_mb))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --memory parameter '%s'.\n", mb_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (malloc_size_mb < S_CAST(intptr_t, kBigstackMinMb)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --memory parameter '%s' (minimum %u).\n", mb_modif, kBigstackMinMb);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
#ifndef __LP64__
          if (malloc_size_mb > S_CAST(intptr_t, kMalloc32bitMbMax)) {
            logerrprintf("Error: --memory parameter too large for 32-bit version (max %u).\n", kMalloc32bitMbMax);
            goto main_ret_INVALID_CMDLINE;
          }
#endif
        } else if (strequal_k_unsafe(flagname_p2, "ake-bed")) {
          if (pc.exportf_flags & kfExportfIndMajorBed) {
            logerrputs("Error: --make-bed cannot be used with --export ind-major-bed.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 3)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "vzs", cur_modif_slen)) {
              make_plink2_flags |= kfMakeBimZs;
            } else if (strequal_k(cur_modif, "trim-alts", cur_modif_slen)) {
              make_plink2_flags |= kfMakePlink2TrimAlts;
            } else if (StrStartsWith(cur_modif, "m=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "multiallelics=", cur_modif_slen)) {
              if (make_plink2_flags & kfMakePlink2MMask) {
                logerrputs("Error: Multiple --make-bed multiallelics= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* mode_start = (cur_modif[1] == '=')? (&(cur_modif[strlen("m=")])) : (&(cur_modif[strlen("multiallelics=")]));
              if (!strcmp(mode_start, "-")) {
                make_plink2_flags |= kfMakePlink2MSplitAll;
              } else if (!strcmp(mode_start, "-snps")) {
                make_plink2_flags |= kfMakePlink2MSplitSnps;
              } else if ((!strcmp(mode_start, "+")) || (!strcmp(mode_start, "+both"))) {
                make_plink2_flags |= kfMakePlink2MMergeBoth;
              } else if (!strcmp(mode_start, "+snps")) {
                make_plink2_flags |= kfMakePlink2MMergeSnps;
              } else if (!strcmp(mode_start, "+any")) {
                make_plink2_flags |= kfMakePlink2MMergeAny;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-bed multiallelics= mode '%s'.\n", mode_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else {
              char* write_iter = strcpya(g_logbuf, "Error: Invalid --make-bed parameter '");
              write_iter = memcpya(write_iter, cur_modif, cur_modif_slen);
              write_iter = strcpya(write_iter, "'.");
              if ((param_idx == 1) && (!outname_end)) {
                // the missing --out mistake is so common--I must have made it
                // over a hundred times by now--that a custom error message is
                // worthwhile.
                write_iter = strcpya(write_iter, " (Did you forget '--out'?)");
              }
              write_iter = strcpya(write_iter, "\n");
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          make_plink2_flags |= kfMakeBed | kfMakeBim | kfMakeFam;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-bpgen")) {
          if (make_plink2_flags & kfMakeBed) {
            logerrputs("Error: --make-bpgen cannot be used with --make-bed.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (import_flags & kfImportKeepAutoconv) {
            logerrputs("Error: --make-bpgen cannot be used with --keep-autoconv.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "vzs", cur_modif_slen)) {
              make_plink2_flags |= kfMakeBimZs;
            } else if (StrStartsWith(cur_modif, "format=", cur_modif_slen)) {
              if (make_plink2_flags & (kfMakePgenFormatBase * 3)) {
                logerrputs("Error: Multiple --make-bpgen format= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const uint32_t fcode_minus_2 = ctou32(cur_modif[7]) - 50;
              if ((fcode_minus_2 > 2) || cur_modif[8]) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-bpgen format code '%s'.\n", &(cur_modif[7]));
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (fcode_minus_2) {
                logerrputs("Error: --make-bpgen formats 3 and 4 (unphased/phased dosage) are not\nimplemented yet.\n");
                reterr = kPglRetNotYetSupported;
                goto main_ret_1;
              }
              make_plink2_flags = S_CAST(MakePlink2Flags, make_plink2_flags | (kfMakePgenFormatBase * (1 + fcode_minus_2)));
            } else if (StrStartsWith(cur_modif, "m=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "multiallelics=", cur_modif_slen)) {
              if (make_plink2_flags & kfMakePlink2MMask) {
                logerrputs("Error: Multiple --make-bpgen multiallelics= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* mode_start = (cur_modif[1] == '=')? (&(cur_modif[strlen("m=")])) : (&(cur_modif[strlen("multiallelics=")]));
              if (!strcmp(mode_start, "-")) {
                make_plink2_flags |= kfMakePlink2MSplitAll;
              } else if (!strcmp(mode_start, "-snps")) {
                make_plink2_flags |= kfMakePlink2MSplitSnps;
              } else if ((!strcmp(mode_start, "+")) || (!strcmp(mode_start, "+both"))) {
                make_plink2_flags |= kfMakePlink2MMergeBoth;
              } else if (!strcmp(mode_start, "+snps")) {
                make_plink2_flags |= kfMakePlink2MMergeSnps;
              } else if (!strcmp(mode_start, "+any")) {
                make_plink2_flags |= kfMakePlink2MMergeAny;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-bpgen multiallelics= mode '%s'.\n", mode_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "trim-alts", cur_modif_slen)) {
              make_plink2_flags |= kfMakePlink2TrimAlts;
            } else if (strequal_k(cur_modif, "erase-alt2+", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenEraseAlt2Plus;
            } else if (strequal_k(cur_modif, "erase-phase", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenErasePhase;
            } else if (strequal_k(cur_modif, "erase-dosage", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenEraseDosage;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-bpgen parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          make_plink2_flags |= kfMakePgen | kfMakeBim | kfMakeFam;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-pgen")) {
          if (make_plink2_flags & (kfMakeBed | kfMakePgen)) {
            logerrputs("Error: --make-pgen cannot be used with --make-bed/--make-bpgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (make_plink2_flags & (kfMakeBim | kfMakeFam | kfMakePvar | kfMakePsam)) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (import_flags & kfImportKeepAutoconv) {
            logerrputs("Error: --make-pgen cannot be used with --keep-autoconv.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4)) {
            goto main_ret_INVALID_CMDLINE_A;
          }
          uint32_t explicit_pvar_cols = 0;
          uint32_t explicit_psam_cols = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "vzs", cur_modif_slen)) {
              pc.pvar_psam_flags |= kfPvarZs;
            } else if (StrStartsWith0(cur_modif, "pvar-cols=", cur_modif_slen)) {
              if (explicit_pvar_cols) {
                logerrputs("Error: Multiple --make-pgen pvar-cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_pvar_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[strlen("pvar-cols=")]), "xheader\0maybequal\0qual\0maybefilter\0filter\0maybeinfo\0info\0maybecm\0cm\0", "--make-pgen pvar-cols", kfPvarColXheader, kfPvarColDefault, 0, &pc.pvar_psam_flags);
              if (reterr) {
                goto main_ret_1;
              }
              if ((pc.pvar_psam_flags & kfPvarColXinfo) && (!(pc.pvar_psam_flags & kfPvarColXheader))) {
                logerrputs("Error: --make-pgen pvar-cols= expression cannot exclude xheader when info is\npresent.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else if (StrStartsWith(cur_modif, "format=", cur_modif_slen)) {
              if (make_plink2_flags & (kfMakePgenFormatBase * 3)) {
                logerrputs("Error: Multiple --make-pgen format= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const uint32_t fcode_minus_2 = ctou32(cur_modif[7]) - 50;
              if ((fcode_minus_2 > 2) || cur_modif[8]) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-pgen format code '%s'.\n", &(cur_modif[7]));
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (fcode_minus_2) {
                logerrputs("Error: --make-pgen formats 3 and 4 (unphased/phased dosage) are not implemented\nyet.\n");
                reterr = kPglRetNotYetSupported;
                goto main_ret_1;
              }
              make_plink2_flags = S_CAST(MakePlink2Flags, make_plink2_flags | (kfMakePgenFormatBase * (1 + fcode_minus_2)));
            } else if (StrStartsWith(cur_modif, "m=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "multiallelics=", cur_modif_slen)) {
              if (make_plink2_flags & kfMakePlink2MMask) {
                logerrputs("Error: Multiple --make-pgen multiallelics= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* mode_start = (cur_modif[1] == '=')? (&(cur_modif[2])) : (&(cur_modif[14]));
              if (!strcmp(mode_start, "-")) {
                make_plink2_flags |= kfMakePlink2MSplitAll;
              } else if (!strcmp(mode_start, "-snps")) {
                make_plink2_flags |= kfMakePlink2MSplitSnps;
              } else if ((!strcmp(mode_start, "+")) || (!strcmp(mode_start, "+both"))) {
                make_plink2_flags |= kfMakePlink2MMergeBoth;
              } else if (!strcmp(mode_start, "+snps")) {
                make_plink2_flags |= kfMakePlink2MMergeSnps;
              } else if (!strcmp(mode_start, "+any")) {
                make_plink2_flags |= kfMakePlink2MMergeAny;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-pgen multiallelics= mode '%s'.\n", mode_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "trim-alts", cur_modif_slen)) {
              make_plink2_flags |= kfMakePlink2TrimAlts;
            } else if (strequal_k(cur_modif, "erase-alt2+", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenEraseAlt2Plus;
            } else if (strequal_k(cur_modif, "erase-phase", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenErasePhase;
            } else if (strequal_k(cur_modif, "erase-dosage", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenEraseDosage;
            } else if (StrStartsWith0(cur_modif, "psam-cols=", cur_modif_slen)) {
              if (explicit_psam_cols) {
                logerrputs("Error: Multiple --make-pgen psam-cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_psam_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[strlen("psam-cols=")]), "maybefid\0fid\0maybesid\0sid\0maybeparents\0parents\0sex\0pheno1\0phenos\0", "make-pgen psam-cols", kfPsamColMaybefid, kfPsamColDefault, 0, &pc.pvar_psam_flags);
              if (reterr) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-pgen parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!explicit_pvar_cols) {
            pc.pvar_psam_flags |= kfPvarColDefault;
          }
          if (!explicit_psam_cols) {
            pc.pvar_psam_flags |= kfPsamColDefault;
          }
          make_plink2_flags |= kfMakePgen | kfMakePvar | kfMakePsam;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-just-bim")) {
          if (make_plink2_flags & (kfMakeBed | kfMakePgen)) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (!strcmp(cur_modif, "zs")) {
              make_plink2_flags |= kfMakeBimZs;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-just-bim parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          make_plink2_flags |= kfMakeBim;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-just-fam")) {
          if (make_plink2_flags & (kfMakeBed | kfMakePgen)) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          make_plink2_flags |= kfMakeFam;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterPsamReq;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ake-just-pvar")) {
          if (make_plink2_flags & (kfMakeBed | kfMakePgen)) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t explicit_cols = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.pvar_psam_flags |= kfPvarZs;
            } else if (StrStartsWith0(cur_modif, "cols=", cur_modif_slen)) {
              if (explicit_cols) {
                logerrputs("Error: Multiple --make-just-pvar cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[5]), "xheader\0maybequal\0qual\0maybefilter\0filter\0maybeinfo\0info\0maybecm\0cm\0", "make-just-pvar", kfPvarColXheader, kfPvarColDefault, 0, &pc.pvar_psam_flags);
              if (reterr) {
                goto main_ret_1;
              }
              if ((pc.pvar_psam_flags & kfPvarColXinfo) && (!(pc.pvar_psam_flags & kfPvarColXheader))) {
                logerrputs("Error: --make-just-pvar cols= expression cannot exclude xheader when info is\npresent.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (pc.pvar_psam_flags & kfPvarColInfo) {
                pc.dependency_flags |= kfFilterNonrefFlagsNeeded;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-just-pvar parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!explicit_cols) {
            pc.pvar_psam_flags |= kfPvarColDefault;
          }
          make_plink2_flags |= kfMakePvar;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-just-psam")) {
          if (make_plink2_flags & (kfMakeBed | kfMakePgen)) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (StrStartsWith0(cur_modif, "cols=", cur_modif_slen)) {
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0maybeparents\0parents\0sex\0pheno1\0phenos\0", "make-just-psam", kfPsamColMaybefid, kfPsamColDefault, 0, &pc.pvar_psam_flags);
              if (reterr) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-just-psam parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.pvar_psam_flags |= kfPsamColDefault;
          }
          make_plink2_flags |= kfMakePsam;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-king")) {
          // may want to add options for handling X/Y/MT
          if (king_cutoff_fprefix) {
            logerrputs("Error: --make-king cannot be used with a --king-cutoff input fileset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          } else if (pc.king_table_subset_fname) {
            logerrputs("Error: --make-king cannot be used with --king-table-subset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "zs")) {
              if (pc.king_flags & kfKingMatrixEncodemask) {
                logerrputs("Error: Multiple --make-king encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixZs;
            } else if (!strcmp(cur_modif, "bin")) {
              if (pc.king_flags & kfKingMatrixEncodemask) {
                logerrputs("Error: Multiple --make-king encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixBin;
            } else if (!strcmp(cur_modif, "bin4")) {
              if (pc.king_flags & kfKingMatrixEncodemask) {
                logerrputs("Error: Multiple --make-king encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixBin4;
            } else if (!strcmp(cur_modif, "square")) {
              if (pc.king_flags & kfKingMatrixShapemask) {
                logerrputs("Error: Multiple --make-king shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixSq;
            } else if (!strcmp(cur_modif, "square0")) {
              if (pc.king_flags & kfKingMatrixShapemask) {
                logerrputs("Error: Multiple --make-king shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixSq0;
            } else if (!strcmp(cur_modif, "triangle")) {
              if (pc.king_flags & kfKingMatrixShapemask) {
                logerrputs("Error: Multiple --make-king shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixTri;
            } else if (!strcmp(cur_modif, "no-idheader")) {
              logerrputs("Error: --make-king 'no-idheader' modifier retired.  Use --no-id-header instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-king parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.king_flags & kfKingMatrixShapemask)) {
            if (pc.king_flags & (kfKingMatrixBin | kfKingMatrixBin4)) {
              pc.king_flags |= kfKingMatrixSq;
            } else {
              pc.king_flags |= kfKingMatrixTri;
            }
          }
          pc.command_flags1 |= kfCommand1MakeKing;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-king-table")) {
          if (king_cutoff_fprefix) {
            logerrputs("Error: --make-king-table cannot be used with a --king-cutoff input fileset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 3)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.king_flags |= kfKingTableZs;
            } else if (strequal_k(cur_modif, "counts", cur_modif_slen)) {
              pc.king_flags |= kfKingCounts;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (pc.king_flags & kfKingColAll) {
                logerrputs("Error: Multiple --make-king-table cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0id\0maybesid\0sid\0nsnp\0hethet\0ibs0\0ibs1\0kinship\0", "make-king-table", kfKingColMaybefid, kfKingColDefault, 1, &pc.king_flags);
              if (reterr) {
                goto main_ret_1;
              }
              if (!(pc.king_flags & kfKingColId)) {
                if (pc.king_flags & (kfKingColMaybefid | kfKingColFid | kfKingColMaybesid | kfKingColSid)) {
                  logerrputs("Error: Invalid --make-king-table column set descriptor ('maybefid', 'fid',\n'maybesid', and 'sid' require 'id').\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
                if (pc.king_table_filter != -DBL_MAX) {
                  logerrputs("Error: --king-table-filter requires --make-king-table cols= to include the 'id'\ncolumn set.\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
                if (pc.king_table_subset_fname) {
                  logerrputs("Error: --king-table-subset requires --make-king-table cols= to include the 'id'\ncolumn set.\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
              }
            } else if (strequal_k(cur_modif, "no-idheader", cur_modif_slen)) {
              logerrputs("Error: --make-king-table 'no-idheader' modifier retired.  Use --no-id-header\ninstead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-king-table parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.king_flags & kfKingColAll)) {
            pc.king_flags |= kfKingColDefault;
          }
          pc.command_flags1 |= kfCommand1MakeKing;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "issing")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.missing_rpt_flags |= kfMissingRptZs;
            } else if (strequal_k(cur_modif, "sample-only", cur_modif_slen)) {
              if (pc.missing_rpt_flags & kfMissingRptVariantOnly) {
                logerrputs("Error: --missing 'sample-only' and 'variant-only' cannot be used together.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              pc.missing_rpt_flags |= kfMissingRptSampleOnly;
            } else if (strequal_k(cur_modif, "variant-only", cur_modif_slen)) {
              if (pc.missing_rpt_flags & kfMissingRptSampleOnly) {
                logerrputs("Error: --missing 'sample-only' and 'variant-only' cannot be used together.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              pc.missing_rpt_flags |= kfMissingRptVariantOnly;
            } else if (StrStartsWith(cur_modif, "scols=", cur_modif_slen)) {
              if (pc.missing_rpt_flags & kfMissingRptScolAll) {
                logerrputs("Error: Multiple --missing scols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("scols=")]), "maybefid\0fid\0maybesid\0sid\0misspheno1\0missphenos\0nmissdosage\0nmiss\0nmisshh\0hethap\0nobs\0fmissdosage\0fmiss\0fmisshh\0", "missing scols", kfMissingRptScolMaybefid, kfMissingRptScolDefault, 1, &pc.missing_rpt_flags);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "vcols=", cur_modif_slen)) {
              if (pc.missing_rpt_flags & kfMissingRptVcolAll) {
                logerrputs("Error: Multiple --missing vcols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("vcols=")]), "chrom\0pos\0ref\0alt1\0alt\0nmissdosage\0nmiss\0nmisshh\0hethap\0nobs\0fmissdosage\0fmiss\0fmisshh\0fhethap\0", "missing vcols", kfMissingRptVcolChrom, kfMissingRptVcolDefault, 1, &pc.missing_rpt_flags);
              if (reterr) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --missing parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const uint32_t explicit_scols = pc.missing_rpt_flags & kfMissingRptScolAll;
          if (pc.missing_rpt_flags & kfMissingRptVariantOnly) {
            if (explicit_scols) {
              logerrputs("Error: --missing 'variant-only' and 'scols=' modifiers cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else {
            if (!explicit_scols) {
              pc.missing_rpt_flags |= kfMissingRptScolDefault;
            }
          }
          const uint32_t explicit_vcols = pc.missing_rpt_flags & kfMissingRptVcolAll;
          if (pc.missing_rpt_flags & kfMissingRptSampleOnly) {
            if (explicit_vcols) {
              logerrputs("Error: --missing 'sample-only' and 'vcols=' modifiers cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else if (!explicit_vcols) {
            pc.missing_rpt_flags |= kfMissingRptVcolDefault;
          }
          pc.command_flags1 |= kfCommand1MissingReport;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "aj-ref")) {
          if (pc.alt1_allele_flag) {
            logerrputs("Error: --maj-ref cannot be used with --ref-allele/--alt1-allele.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (!strcmp(cur_modif, "force")) {
              pc.misc_flags |= kfMiscMajRefForce;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --maj-ref parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.misc_flags |= kfMiscMajRef;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "af")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (!ScanadvDouble(cur_modif, &pc.min_maf)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --maf parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (pc.min_maf < 0.0) {
              snprintf(g_logbuf, kLogbufSize, "Error: --maf parameter '%s' too small (must be >= 0).\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else if (pc.min_maf >= 1.0) {
              snprintf(g_logbuf, kLogbufSize, "Error: --maf parameter '%s' too large (must be < 1).\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.min_maf = 0.01;
          }
          if (pc.min_maf != 0.0) {
            pc.filter_flags |= kfFilterPvarReq;
            pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ax-maf")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (!ScanadvDouble(cur_modif, &pc.max_maf)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-maf parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (pc.max_maf < pc.min_maf) {
            snprintf(g_logbuf, kLogbufSize, "Error: --max-maf parameter '%s' too small (must be >= %g).\n", cur_modif, pc.min_maf);
            goto main_ret_INVALID_CMDLINE_WWA;
          } else if (pc.max_maf >= 1.0) {
            snprintf(g_logbuf, kLogbufSize, "Error: --max-maf parameter '%s' too large (must be < 1).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "ac")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if ((!ScanadvDouble(cur_modif, &dxx)) || (dxx < 0.0) || (dxx > 2147483646.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mac parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE;
          }
          if (dxx > 0.0) {
            // round up, but keep as much precision as possible
            int32_t int_part = S_CAST(int32_t, dxx);
            dxx -= int_part;
            pc.min_allele_dosage = int_part * S_CAST(uint64_t, kDosageMax);
            if (dxx > 0.0) {
              pc.min_allele_dosage += 1 + (dxx * (kDosageMax * (1 - kSmallEpsilon)));
            }
            pc.filter_flags |= kfFilterPvarReq;
            pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ax-mac")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if ((!ScanadvDouble(cur_modif, &dxx)) || (dxx < 0.0) || (dxx > 2147483646.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-mac parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          // round down
          pc.max_allele_dosage = S_CAST(int64_t, dxx * kDosageMax);
          if (pc.max_allele_dosage < pc.min_allele_dosage) {
            // yeah, --mac 0.1 --max-mac 0.1 also isn't allowed
            logerrputs("Error: --max-mac parameter cannot be smaller than --mac parameter.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "ind")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t mind_thresh_present = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "dosage")) {
              pc.misc_flags |= kfMiscMindDosage;
            } else if (!strcmp(cur_modif, "hh-missing")) {
              pc.misc_flags |= kfMiscMindHhMissing;
            } else if (mind_thresh_present) {
              logerrputs("Error: Invalid --mind parameter sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (!ScanadvDouble(cur_modif, &pc.mind_thresh)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mind parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else if ((pc.mind_thresh < 0.0) || (pc.mind_thresh > 1.0)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mind parameter '%s' (must be in [0, 1]).\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else {
              mind_thresh_present = 1;
            }
          }
          if (!mind_thresh_present) {
            pc.mind_thresh = 0.1;
          }
          if (pc.mind_thresh < 1.0) {
            pc.filter_flags |= kfFilterPsamReq;
            pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
          }
        } else if (strequal_k_unsafe(flagname_p2, "issing-var-code")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.missing_varid_match);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "erge-par")) {
          if (pc.exportf_flags & kfExportfVcf) {
            logerrputs("Warning: --merge-par should not be used with VCF export.  (The VCF export\nroutine automatically converts PAR1/PAR2 chromosome codes to X, while using\nthe PAR boundaries to get male ploidy right; --merge-par causes VCF export to\nget male ploidy wrong.)\n");
          }
          pc.misc_flags |= kfMiscMergePar;
          pc.dependency_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "erge-x")) {
          if (pc.misc_flags & kfMiscMergePar) {
            logerrputs("Error: --merge-par cannot be used with --merge-x.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.exportf_flags & kfExportfVcf) {
            logerrputs("Warning: --merge-x should not be used in the same run as VCF export; this\ncauses some ploidies to be wrong.  Instead, use --merge-x + --sort-vars +\n--make-{b}pgen in one run, and follow up with --split-par + --export vcf.\n");
          }
          pc.misc_flags |= kfMiscMergeX;
          pc.dependency_flags |= kfFilterPvarReq;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "af-succ")) {
          pc.misc_flags |= kfMiscMafSucc;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ax-corr")) {
          if (!(pc.command_flags1 & kfCommand1Glm)) {
            logerrputs("Error: --max-corr must be used with --glm.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (!ScanadvDouble(cur_modif, &pc.glm_info.max_corr)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-corr parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if ((pc.glm_info.max_corr < 0.0) || (pc.glm_info.max_corr > 1.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-corr parameter '%s' (must be in [0, 1]).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ach-r2-filter")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (!ScanadvDouble(cur_modif, &pc.mach_r2_min)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mach-r2-filter min parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (pc.mach_r2_min < 0.0) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mach-r2-filter min parameter '%s' (must be nonnegative).\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (param_ct == 2) {
              cur_modif = argvk[arg_idx + 2];
              if (!ScanadvDouble(cur_modif, &pc.mach_r2_max)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mach-r2-filter max parameter '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else {
              pc.mach_r2_max = 2.0;
            }
            if (pc.mach_r2_max < pc.mach_r2_min) {
              logerrputs("Error: --mach-r2-filter min parameter cannot be larger than max parameter.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else {
            pc.mach_r2_min = 0.1;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "issing-code")) {
          if (!(xload & (kfXloadOxGen | kfXloadOxBgen))) {
            // could technically support pure .sample -> .fam/.psam, but let's
            // keep this simple
            logerrputs("Error: --missing-code must be used with --data/--gen/--bgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(param_ct? argvk[arg_idx + 1] : "", argvk[arg_idx], 0x7fffffff, &ox_missing_code);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "issing-genotype")) {
          logerrputs("Error: --missing-genotype flag retired.  Use --input-missing-genotype and/or\n--output-missing-genotype.\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (strequal_k_unsafe(flagname_p2, "issing-phenotype")) {
          logerrputs("Error: --missing-phenotype flag retired.  Use --input-missing-phenotype and/or\n--output-missing-phenotype.\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (strequal_k_unsafe(flagname_p2, "issing-catname")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          double dxx;
          if (ScanadvDouble(cur_modif, &dxx) || IsNanStr(cur_modif, cur_modif_slen)) {
            logerrputs("Error: --missing-catname string cannot be 'NA' or start with a number.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (cur_modif_slen > 31) {
            logerrputs("Error: --missing-catname string too long (max 31 chars).\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          memcpy(g_missing_catname, cur_modif, cur_modif_slen + 1);
        } else if (strequal_k_unsafe(flagname_p2, "ouse")) {
          if (chr_info.chrset_source) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.autosome_ct = 19;
          chr_info.xymt_codes[0] = 20;
          chr_info.xymt_codes[1] = 21;
          chr_info.xymt_codes[2] = UINT32_MAXM1;
          chr_info.xymt_codes[3] = UINT32_MAXM1;
          chr_info.xymt_codes[4] = UINT32_MAXM1;
          chr_info.xymt_codes[5] = UINT32_MAXM1;
          chr_info.haploid_mask[0] = 0x300000;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ake-grm")) {
          logerrputs("Error: --make-grm has been retired due to inconsistent meaning across GCTA\nversions.  Use --make-grm-list or --make-grm-bin.\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (strequal_k_unsafe(flagname_p2, "ake-grm-bin")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 3)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          pc.grm_flags |= kfGrmNoIdHeader | kfGrmBin;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "cov")) {
              pc.grm_flags |= kfGrmCov;
            } else if (!strcmp(cur_modif, "meanimpute")) {
              pc.grm_flags |= kfGrmMeanimpute;
            } else if ((!strcmp(cur_modif, "id-header")) || (!strcmp(cur_modif, "idheader"))) {
              pc.grm_flags &= ~kfGrmNoIdHeader;
            } else if (!strcmp(cur_modif, "iid-only")) {
              pc.grm_flags |= kfGrmNoIdHeaderIidOnly;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-grm-bin parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if ((pc.grm_flags & (kfGrmNoIdHeader | kfGrmNoIdHeaderIidOnly)) == kfGrmNoIdHeaderIidOnly) {
            logerrputs("Error: --make-grm-bin 'id-header' and 'iid-only' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.command_flags1 |= kfCommand1MakeRel;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-grm-gz") || strequal_k_unsafe(flagname_p2, "ake-grm-list")) {
          if (pc.command_flags1 & kfCommand1MakeRel) {
            if (pc.grm_flags & kfGrmBin) {
              logerrputs("Error: --make-grm-list cannot be used with --make-grm-bin.\n");
            } else {
              logerrputs("Error: --make-grm-list cannot be used with --make-grm-gz.\n");
            }
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t compress_stream_type = 0;  // 1 = no-gz, 2 = zs
          pc.grm_flags |= kfGrmNoIdHeader;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "cov")) {
              pc.grm_flags |= kfGrmCov;
            } else if (!strcmp(cur_modif, "meanimpute")) {
              pc.grm_flags |= kfGrmMeanimpute;
            } else if (!strcmp(cur_modif, "no-gz")) {
              if (compress_stream_type) {
                logerrputs("Error: Multiple --make-grm-list compression type modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              compress_stream_type = 1;
              pc.grm_flags |= kfGrmListNoGz;
            } else if (!strcmp(cur_modif, "zs")) {
              if (compress_stream_type) {
                logerrputs("Error: Multiple --make-grm-list compression type modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              compress_stream_type = 2;
              pc.grm_flags |= kfGrmListZs;
            } else if ((!strcmp(cur_modif, "id-header")) || (!strcmp(cur_modif, "idheader"))) {
              pc.grm_flags &= ~kfGrmNoIdHeader;
            } else if (!strcmp(cur_modif, "iid-only")) {
              pc.grm_flags |= kfGrmNoIdHeaderIidOnly;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-grm-list parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (flagname_p2[8] == 'g') {
            if (!compress_stream_type) {
              // screw it, life is too much better with multithreaded .zst
              logerrputs("Error: --make-grm-list no longer supports gzipped output.  Use 'zs' for\nzstd-compressed output (much faster), or use PLINK 1.9 for this function.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            logerrputs("Warning: --make-grm-gz has been renamed to --make-grm-list.\n");
          } else if (!compress_stream_type) {
            compress_stream_type = 1;
            pc.grm_flags |= kfGrmListNoGz;
          }
          if ((pc.grm_flags & (kfGrmNoIdHeader | kfGrmNoIdHeaderIidOnly)) == kfGrmNoIdHeaderIidOnly) {
            logerrputs("Error: --make-grm-list 'id-header' and 'iid-only' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.command_flags1 |= kfCommand1MakeRel;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-rel")) {
          if (pc.command_flags1 & kfCommand1MakeRel) {
            logerrputs("Error: --make-rel cannot be used with --make-grm-list/--make-grm-bin.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "cov")) {
              pc.grm_flags |= kfGrmCov;
            } else if (!strcmp(cur_modif, "meanimpute")) {
              pc.grm_flags |= kfGrmMeanimpute;
            } else if (!strcmp(cur_modif, "zs")) {
              if (pc.grm_flags & kfGrmMatrixEncodemask) {
                logerrputs("Error: Multiple --make-rel encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixZs;
            } else if (!strcmp(cur_modif, "no-idheader")) {
              logerrputs("Error: --make-rel 'no-idheader' modifier retired.  Use --no-id-header instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (!strcmp(cur_modif, "bin")) {
              if (pc.grm_flags & kfGrmMatrixEncodemask) {
                logerrputs("Error: Multiple --make-rel encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixBin;
            } else if (!strcmp(cur_modif, "bin4")) {
              if (pc.grm_flags & kfGrmMatrixEncodemask) {
                logerrputs("Error: Multiple --make-rel encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixBin4;
            } else if (!strcmp(cur_modif, "square")) {
              if (pc.grm_flags & kfGrmMatrixShapemask) {
                logerrputs("Error: Multiple --make-rel shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixSq;
            } else if (!strcmp(cur_modif, "square0")) {
              if (pc.grm_flags & kfGrmMatrixShapemask) {
                logerrputs("Error: Multiple --make-rel shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixSq0;
            } else if (!strcmp(cur_modif, "triangle")) {
              if (pc.grm_flags & kfGrmMatrixShapemask) {
                logerrputs("Error: Multiple --make-rel shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixTri;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-rel parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.grm_flags & kfGrmMatrixShapemask)) {
            if (pc.grm_flags & (kfGrmMatrixBin | kfGrmMatrixBin4)) {
              pc.grm_flags |= kfGrmMatrixSq;
            } else {
              pc.grm_flags |= kfGrmMatrixTri;
            }
          }
          pc.command_flags1 |= kfCommand1MakeRel;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ap")) {
          if (load_params || (xload & (~kfXloadPlink1Dosage))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_modif);
          if (slen > kPglFnamesize - 1) {
            logerrputs("Error: --map filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pvarname, cur_modif, slen + 1);
          xload |= kfXloadMap;
        } else if (strequal_k_unsafe(flagname_p2, "within")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (ScanPosintCapped(cur_modif, kMaxLongLine / 2, &pc.mwithin_val)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mwithin parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "filter")) {
          if (!pc.keep_fcol_fname) {
            logerrputs("Error: --mfilter must be used with --keep-fcol.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.keep_fcol_name) {
            logerrputs("Error: --mfilter can't be used with --keep-fcol-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.keep_fcol_num) {
            logerrputs("Error: --mfilter can't be used with --keep-fcol-num.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          logerrputs("Warning: --mfilter flag deprecated.  Use --keep-fcol-num or --keep-fcol-name\ninstead.  (Note that --keep-fcol-num does not add 2 to the column number.)\n");
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t mfilter_arg;
          if (ScanPosintDefcap(cur_modif, &mfilter_arg)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mfilter parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.keep_fcol_num = mfilter_arg + 2;
        } else if (strequal_k_unsafe(flagname_p2, "pheno")) {
          logerrputs("Warning: --mpheno flag deprecated.  Use --pheno-col-nums instead.  (Note that\n--pheno-col-nums does not add 2 to the column number(s).)\n");
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t mpheno_arg;
          if (ScanPosintDefcap(cur_modif, &mpheno_arg)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mpheno parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          // add two to the number, convert it back to a string, and act as if
          // --pheno-col-nums was used.  See ParseNameRanges().
          const uint32_t pheno_col_nums_arg = mpheno_arg + 2;
          const uint32_t name_max_blen = UintSlen(pheno_col_nums_arg) + 1;
          if (pgl_malloc(name_max_blen + 1, &pc.pheno_range_list.names)) {
            goto main_ret_NOMEM;
          }
          pc.pheno_range_list.name_max_blen = name_max_blen;
          pc.pheno_range_list.name_ct = 1;
          char* write_iter = u32toa(pheno_col_nums_arg, pc.pheno_range_list.names);
          *write_iter++ = '\0';
          pc.pheno_range_list.starts_range = R_CAST(unsigned char*, write_iter);
          *write_iter = '\0';
          pc.misc_flags |= kfMiscPhenoColNums;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'n':
        if (strequal_k_unsafe(flagname_p2, "o-fid")) {
          pc.fam_cols &= ~kfFamCol1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "o-parents")) {
          pc.fam_cols &= ~kfFamCol34;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "o-sex")) {
          pc.fam_cols &= ~kfFamCol5;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "o-psam-pheno")) {
          // move this out of fam_cols?
          pc.fam_cols &= ~kfFamCol6;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "onfounders")) {
          pc.misc_flags |= kfMiscNonfounders;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ot-chr")) {
          if (pc.varid_from) {
            logerrputs("Error: --from/--to cannot be used with --autosome{-par} or --{not-}chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.from_bp != -1) {
            logerrputs("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb cannot be used with --not-chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }

          // allowed:
          //   --allow-extra-chr --chr 5-22 bobs_chrom --not-chr 17
          // allowed:
          //   --allow-extra-chr --not-chr 12-17 bobs_chrom
          // does not make sense, disallowed:
          //   --allow-extra-chr --chr 5-22 --not-chr bobs_chrom
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }

          // --allow-extra-chr present, --chr/--autosome{-xy} not present
          const uint32_t aec_and_no_chr_include = ((pc.misc_flags / kfMiscAllowExtraChrs) & 1) && (!chr_info.is_include_stack);
          reterr = ParseChrRanges(&(argvk[arg_idx]), flagname_p, errstr_append, param_ct, aec_and_no_chr_include, kChrRawEnd - (kChrExcludeWords * kBitsPerWord), '-', &chr_info, chr_info.chr_exclude);
          if (reterr) {
            goto main_ret_1;
          }
          notchr_present = 1;
          // remaining processing now postponed to FinalizeChrset()
        } else if (strequal_k_unsafe(flagname_p2, "ew-id-max-allele-len")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (ScanPosintCapped(cur_modif, kMaxIdSlen - 2, &pc.new_variant_id_max_allele_slen)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --new-id-max-allele-len length parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (param_ct == 2) {
            cur_modif = argvk[arg_idx + 2];
            if (!strcmp(cur_modif, "missing")) {
              pc.misc_flags |= kfMiscNewVarIdOverflowMissing;
            } else if (!strcmp(cur_modif, "truncate")) {
              pc.misc_flags |= kfMiscNewVarIdOverflowTruncate;
            } else if (strcmp(cur_modif, "error")) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --new-id-max-allele-len parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "o-id-header")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (strcmp(cur_modif, "iid-only")) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --no-id-header parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            pc.misc_flags |= kfMiscNoIdHeaderIidOnly;
          }
          pc.misc_flags |= kfMiscNoIdHeader;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'o':
        if (strequal_k_unsafe(flagname_p2, "utput-chr")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* mt_code = argvk[arg_idx + 1];
          if (!strcmp(mt_code, "M")) {
            chr_info.output_encoding = kfChrOutputM;
          } else if (!strcmp(mt_code, "MT")) {
            chr_info.output_encoding = kfChrOutputMT;
          } else if (!strcmp(mt_code, "0M")) {
            chr_info.output_encoding = kfChrOutput0M;
          } else if (!strcmp(mt_code, "chr26")) {
            chr_info.output_encoding = kfChrOutputPrefix;
          } else if (!strcmp(mt_code, "chrM")) {
            chr_info.output_encoding = kfChrOutputPrefix | kfChrOutputM;
          } else if (!strcmp(mt_code, "chrMT")) {
            chr_info.output_encoding = kfChrOutputPrefix | kfChrOutputMT;
          } else if (!strcmp(mt_code, "26")) {
            chr_info.output_encoding = kfChrOutput0;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --output-chr parameter '%s'.\n", mt_code);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "utput-min-p")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if ((!ScanadvDouble(cur_modif, &pc.output_min_p)) || (!(pc.output_min_p >= 0.0)) || (pc.output_min_p >= 1.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --output-min-p parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "xford-single-chr")) {
          if (!(xload & (kfXloadOxGen | kfXloadOxBgen))) {
            logerrputs("Error: --oxford-single-chr must be used with .gen/.bgen input.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (oxford_import_flags & kfOxfordImportBgenSnpIdChr) {
            logerrputs("Error: --oxford-single-chr cannot be used with --bgen 'snpid-chr'.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (!(pc.misc_flags & kfMiscAllowExtraChrs)) {
            if (IsI32Neg(GetChrCodeRaw(cur_modif))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --oxford-single-chr chromosome code '%s'. (Did you forget --allow-extra-chr?)\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          reterr = CmdlineAllocString(cur_modif, argvk[arg_idx], kMaxIdSlen, &import_single_chr_str);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "utput-missing-genotype")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          output_missing_geno_char = ExtractCharParam(cur_modif);
          if (ctou32(output_missing_geno_char) <= ' ') {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --output-missing-genotype parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "utput-missing-phenotype")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          if (cur_modif_slen > 31) {
            logerrputs("Error: --output-missing-phenotype string too long (max 31 chars).\n");
            goto main_ret_INVALID_CMDLINE;
          }
          memcpy(g_output_missing_pheno, cur_modif, cur_modif_slen + 1);
        } else if (!strequal_k_unsafe(flagname_p2, "ut")) {
          // --out is a special case due to logging
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'p':
        if (strequal_k_unsafe(flagname_p2, "file")) {
          if (load_params || xload) {
            // currently only possible with --bcf, --bfile, --pfile
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (CheckExtraParam(&(argvk[arg_idx]), "vzs", &fname_modif_idx)) {
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* fname_prefix = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname_prefix);
          if (slen > (kPglFnamesize - 11)) {
            logerrputs("Error: --pfile parameter too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          snprintf(memcpya(pgenname, fname_prefix, slen), 10, ".pgen");
          snprintf(memcpya(psamname, fname_prefix, slen), 10, ".psam");
          char* pvarname_end = memcpya(pvarname, fname_prefix, slen);
          pvarname_end = Stpcpy(pvarname_end, ".pvar");
          if (param_ct == 2) {
            snprintf(pvarname_end, 5, ".zst");
          }
          load_params |= kfLoadParamsPfileAll;
        } else if (strequal_k_unsafe(flagname_p2, "gen")) {
          if (xload) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          load_params |= kfLoadParamsPgen;
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (slen > (kPglFnamesize - 2)) {
            logerrputs("Error: --pgen parameter too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, fname, slen + 1);
        } else if (strequal_k_unsafe(flagname_p2, "sam")) {
          if (xload & (~(kfXloadVcf | kfXloadBcf | kfXloadPlink1Dosage | kfXloadMap))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          load_params |= kfLoadParamsPsam;
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (slen > (kPglFnamesize - 2)) {
            logerrputs("Error: --psam parameter too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(psamname, fname, slen + 1);
        } else if (strequal_k_unsafe(flagname_p2, "var")) {
          if (xload) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          load_params |= kfLoadParamsPvar;
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (slen > (kPglFnamesize - 2)) {
            logerrputs("Error: --pvar parameter too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pvarname, fname, slen + 1);
        } else if (strequal_k_unsafe(flagname_p2, "heno")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &pc.pheno_fname);
          if (reterr) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "heno-col-nums")) {
          if (!pc.pheno_fname) {
            logerrputs("Error: --pheno-col-nums must be used with --pheno.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.pheno_range_list.name_ct) {
            logerrputs("Error: --pheno-col-nums can't be used with --mpheno.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, range_delim, &pc.pheno_range_list);
          if (reterr) {
            goto main_ret_1;
          }
          pc.misc_flags |= kfMiscPhenoColNums;
        } else if (strequal_k_unsafe(flagname_p2, "heno-name")) {
          if (pc.pheno_range_list.name_ct) {
            logerrputs("Error: --pheno-name can't be used with --pheno-col-nums.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // can now be used without --pheno
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.pheno_range_list);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "arallel")) {
          if (pc.king_flags & kfKingMatrixSq) {
            logerrputs("Error: --parallel cannot be used with '--make-king square'.  Use '--make-king\nsquare0' or plain --make-king instead.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if ((pc.king_cutoff != -1) && (!king_cutoff_fprefix)) {
            logerrputs("Error: --parallel cannot be used with --king-cutoff.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.grm_flags & kfGrmMatrixSq) {
            logerrputs("Error: --parallel cannot be used with '--make-rel square'.  Use '--make-rel\nsquare0' or plain --make-rel instead.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (ScanPosintCapped(argvk[arg_idx + 1], kParallelMax, &pc.parallel_idx)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --parallel job index '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (ScanPosintCapped(argvk[arg_idx + 2], kParallelMax, &pc.parallel_tot) || (pc.parallel_tot == 1) || (pc.parallel_tot < pc.parallel_idx)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --parallel total job count '%s'.\n", argvk[arg_idx + 2]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          --pc.parallel_idx;  // internal 0..(n-1) indexing
        } else if (strequal_k_unsafe(flagname_p2, "arameters")) {
          if (!(pc.command_flags1 & kfCommand1Glm)) {
            logerrputs("Error: --parameters must be used with --glm.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.glm_info.parameters_range_list);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "filter")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (!ScanadvDouble(cur_modif, &pc.pfilter)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --pfilter parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if ((pc.pfilter <= 0.0) || (pc.pfilter > 1.0)) {
            logerrputs("Error: --pfilter threshold must be in (0, 1].\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ca")) {
#ifdef NOLAPACK
          logerrputs("Error: --pca requires " PROG_NAME_STR " to be built with LAPACK.\n");
          goto main_ret_INVALID_CMDLINE;
#endif
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 6)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t explicit_scols = 0;
          uint32_t is_var_wts = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "approx", cur_modif_slen)) {
              pc.pca_flags |= kfPcaApprox;
            } else if (strequal_k(cur_modif, "meanimpute", cur_modif_slen)) {
              pc.pca_flags |= kfPcaMeanimpute;
            } else if (strequal_k(cur_modif, "var-wts", cur_modif_slen)) {
              is_var_wts = 1;
            } else if (strequal_k(cur_modif, "vzs", cur_modif_slen)) {
              pc.pca_flags |= kfPcaVarZs;
            } else if (StrStartsWith(cur_modif, "scols=", cur_modif_slen)) {
              if (explicit_scols) {
                logerrputs("Error: Multiple --pca scols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("scols=")]), "maybefid\0fid\0maybesid\0sid\0", "pca scols", kfPcaScolMaybefid, kfPcaScolDefault, 0, &pc.pca_flags);
              if (reterr) {
                goto main_ret_1;
              }
              explicit_scols = 1;
            } else if (StrStartsWith(cur_modif, "vcols=", cur_modif_slen)) {
              if (pc.pca_flags & kfPcaVcolAll) {
                logerrputs("Error: Multiple --pca vcols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("vcols=")]), "chrom\0pos\0ref\0alt1\0alt\0maj\0nonmaj\0", "pca vcols", kfPcaVcolChrom, kfPcaVcolDefault, 1, &pc.pca_flags);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (strequal_k(cur_modif, "sid", cur_modif_slen)) {
              logerrputs("Error: --pca 'sid' modifier retired.  Use --pca scols= instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else {
              if (pc.pca_ct || ScanPosintDefcap(cur_modif, &pc.pca_ct)) {
                logerrputs("Error: Invalid --pca parameter sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (pc.pca_ct > 8000) {
                // this slightly simplifies output buffering.
                // lower limit for randomized algorithm?
                // (just let memory allocation fail for now...)
                logerrputs("Error: --pca does not support more than 8000 PCs.\n");
                goto main_ret_INVALID_CMDLINE;
              }
            }
          }
          if (pc.pca_flags & kfPcaApprox) {
            if (pc.pca_ct > 100) {
              // double-precision overflow too likely
              logerrputs("Error: --pca approx does not support more than 100 PCs.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          } else {
            // todo: if --make-rel/--make-grm present, verify consistency
            if (pc.parallel_tot != 1) {
              logerrputs("Error: Non-approximate --pca cannot be used with --parallel.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            const uint32_t pca_meanimpute = (pc.pca_flags / kfPcaMeanimpute) & 1;
            if (pc.command_flags1 & kfCommand1MakeRel) {
              if (((pc.grm_flags / kfGrmMeanimpute) & 1) != pca_meanimpute) {
                logerrputs("Error: --make-rel/--make-grm-list/--make-grm-bin meanimpute setting must match\n--pca meanimpute setting.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              if (pc.grm_flags & kfGrmCov) {
                logerrputs("Error: --make-rel/--make-grm-list/--make-grm-bin cannot be used to compute a\ncovariance matrix in the same run as non-approximate --pca.\n");
                goto main_ret_INVALID_CMDLINE;
              }
            } else {
              if (pca_meanimpute) {
                pc.grm_flags |= kfGrmMeanimpute;
              }
            }
          }
          if (!pc.pca_ct) {
            pc.pca_ct = 10;
          }
          if (!explicit_scols) {
            pc.pca_flags |= kfPcaScolDefault;
          }
          if (!(pc.pca_flags & kfPcaVcolAll)) {
            if (is_var_wts) {
              pc.pca_flags |= kfPcaVcolDefault;
            }
          } else if (!is_var_wts) {
            logerrputs("Error: --pca 'vcols=' has no effect without 'var-wts'.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (is_var_wts) {
            pc.pca_flags |= kfPcaVarWts;
          } else if (pc.pca_flags & kfPcaVarZs) {
            logerrputs("Error: --pca 'vzs' modifier has no effect without 'var-wts'.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.command_flags1 |= kfCommand1Pca;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "heno-quantile-normalize")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &pc.quantnorm_flattened);
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.pheno_transform_flags |= kfPhenoTransformQuantnormPheno;
          pc.dependency_flags |= kfFilterPsamReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'q':
        if (strequal_k_unsafe(flagname_p2, "uantile-normalize")) {
          if (pc.pheno_transform_flags & (kfPhenoTransformQuantnormPheno | kfPhenoTransformQuantnormCovar)) {
            logerrputs("Error: --quantile-normalize cannot be used with --pheno-quantile-normalize or\n--covar-quantile-normalize.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &pc.quantnorm_flattened);
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.pheno_transform_flags |= kfPhenoTransformQuantnormAll;
          pc.dependency_flags |= kfFilterPsamReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'r':
        if (strequal_k_unsafe(flagname_p2, "eal-ref-alleles")) {
          if (pc.misc_flags & kfMiscMajRef) {
            logerrputs("Error: --real-ref-alleles cannot be used with --maj-ref.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.misc_flags |= kfMiscRealRefAlleles;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "emove")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, kPglFnamesize, &pc.remove_fnames);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-fam")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, kPglFnamesize, &pc.removefam_fnames);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-females")) {
          pc.filter_flags |= kfFilterPsamReq | kfFilterExclFemales;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "emove-males")) {
          pc.filter_flags |= kfFilterPsamReq | kfFilterExclMales;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "emove-nosex")) {
          pc.filter_flags |= kfFilterPsamReq | kfFilterExclNosex;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ead-freq")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &pc.read_freq_fname);
          if (reterr) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "equire-pheno")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &pc.require_pheno_flattened);
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.misc_flags |= kfMiscRequirePheno;
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "equire-covar")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &pc.require_covar_flattened);
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.misc_flags |= kfMiscRequireCovar;
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "equire-info")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenCommaDelim(&(argvk[arg_idx + 1]), param_ct, &pc.require_info_flattened);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "equire-no-info")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenCommaDelim(&(argvk[arg_idx + 1]), param_ct, &pc.require_no_info_flattened);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-if")) {
          reterr = ValidateAndAllocCmpExpr(&(argvk[arg_idx + 1]), argvk[arg_idx], param_ct, &pc.remove_if_expr);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-cats")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &pc.remove_cats_fname);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-cat-names")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, kMaxIdBlen, &pc.remove_cat_names_flattened);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-cat-pheno")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.remove_cat_phenoname);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ef-allele")) {
          if (pc.misc_flags & kfMiscMajRef) {
            logerrputs("Error: --maj-ref cannot be used with --ref-allele/--alt1-allele.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 5)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* const* sources = &(argvk[arg_idx + 1]);
          if (!strcmp(sources[0], "force")) {
            --param_ct;
            if (!param_ct) {
              logputs("Error: Invalid --ref-allele parameter sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.misc_flags |= kfMiscRefAlleleForce;
            ++sources;
          }
          reterr = Alloc2col(sources, flagname_p, param_ct, &pc.ref_allele_flag);
          if (reterr) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "ef-from-fa")) {
          if ((pc.misc_flags & kfMiscMajRef) || pc.alt1_allele_flag) {
            // could allow --alt1-allele later, but keep this simpler for now
            logerrputs("Error: --ref-from-fa cannot be used with --maj-ref or\n--ref-allele/--alt1-allele.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (CheckExtraParam(&(argvk[arg_idx]), "force", &fname_modif_idx)) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.misc_flags |= kfMiscRefFromFaForce;
          }
          reterr = AllocFname(argvk[arg_idx + fname_modif_idx], flagname_p, 0, &pc.ref_from_fa_fname);
          if (reterr) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq | kfFilterNonrefFlagsNeeded;
        } else if (strequal_k_unsafe(flagname_p2, "andmem")) {
          randmem = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ice")) {
          if (chr_info.chrset_source) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.autosome_ct = 12;
          chr_info.xymt_codes[0] = UINT32_MAXM1;
          chr_info.xymt_codes[1] = UINT32_MAXM1;
          chr_info.xymt_codes[2] = UINT32_MAXM1;
          chr_info.xymt_codes[3] = UINT32_MAXM1;
          chr_info.xymt_codes[4] = UINT32_MAXM1;
          chr_info.xymt_codes[5] = UINT32_MAXM1;
          chr_info.haploid_mask[0] = 0x1fff;
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 's':
        if (strequal_k_unsafe(flagname_p2, "eed")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          rseed_ct = param_ct;
          if (pgl_malloc(param_ct * sizeof(int32_t), &rseeds)) {
            goto main_ret_NOMEM;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (ScanUintCapped(cur_modif, UINT32_MAX, &(rseeds[param_idx - 1]))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --seed parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "plit-par")) {
          if (pc.misc_flags & (kfMiscMergePar | kfMiscMergeX)) {
            logerrputs("Error: --split-par cannot be used with --merge-par/--merge-x.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 1) {
            const char* build_code = argvk[arg_idx + 1];
            if ((!strcmp(build_code, "b38")) || (!strcmp(build_code, "hg38"))) {
              pc.splitpar_bound1 = 2781479;
              pc.splitpar_bound2 = 155701383;
            } else if ((!strcmp(build_code, "b37")) || (!strcmp(build_code, "hg19"))) {
              pc.splitpar_bound1 = 2699520;
              pc.splitpar_bound2 = 154931044;
            } else if ((!strcmp(build_code, "b36")) || (!strcmp(build_code, "hg18"))) {
              pc.splitpar_bound1 = 2709521;
              pc.splitpar_bound2 = 154584237;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Unrecognized --split-par build code '%s'.\n", build_code);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            if (ScanUintDefcap(argvk[arg_idx + 1], &pc.splitpar_bound1)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --split-par parameter '%s'.\n", argvk[arg_idx + 1]);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (ScanUintDefcap(argvk[arg_idx + 2], &pc.splitpar_bound2) || (pc.splitpar_bound2 <= pc.splitpar_bound1)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --split-par parameter '%s'.\n", argvk[arg_idx + 2]);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.dependency_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "et-all-var-ids") || strequal_k_unsafe(flagname_p2, "et-missing-var-ids")) {
          if (flagname_p2[3] == 'm') {
            if (pc.varid_template) {
              logerrputs("Error: --set-missing-var-ids cannot be used with --set-all-var-ids.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            pc.misc_flags |= kfMiscSetMissingVarIds;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (!VaridTemplateIsValid(argvk[arg_idx + 1], flagname_p)) {
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_template);
          if (reterr) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "et-hh-missing")) {
          if (!(pc.command_flags1 & kfCommand1MakePlink2)) {
            logerrputs("Error: --set-hh-missing must be used with --make-{b}pgen/--make-bed.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (!strcmp(cur_modif, "keep-dosage")) {
              if (make_plink2_flags & kfMakePgenEraseDosage) {
                logerrputs("Error: --set-hh-missing 'keep-dosage' modifier cannot be used with\n--make-{b}pgen erase-dosage.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              make_plink2_flags |= kfMakePlink2SetHhMissingKeepDosage;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --set-hh-missing parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          make_plink2_flags |= kfMakePlink2SetHhMissing;
        } else if (strequal_k_unsafe(flagname_p2, "et-mixed-mt-missing")) {
          if (!(pc.command_flags1 & kfCommand1MakePlink2)) {
            logerrputs("Error: --set-mixed-mt-missing must be used with --make-{b}pgen/--make-bed.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (!strcmp(cur_modif, "keep-dosage")) {
              if (make_plink2_flags & kfMakePgenEraseDosage) {
                logerrputs("Error: --set-mixed-mt-missing 'keep-dosage' modifier cannot be used with\n--make-{b}pgen erase-dosage.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              make_plink2_flags |= kfMakePlink2SetMixedMtMissingKeepDosage;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --set-mixed-mt-missing parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          make_plink2_flags |= kfMakePlink2SetMixedMtMissing;
        } else if (strequal_k_unsafe(flagname_p2, "ample")) {
          if (load_params || (xload & (~(kfXloadOxGen | kfXloadOxBgen | kfXloadOxHaps | kfXloadOxLegend | kfXloadOxSample)))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (!(xload & (kfXloadOxGen | kfXloadOxBgen | kfXloadOxHaps))) {
            logerrputs("Error: --sample must be used with --gen/--bgen/--data/--haps.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_fname);
          if (slen > kPglFnamesize - 1) {
            logerrputs("Error: --sample filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(psamname, cur_fname, slen + 1);
          xload |= kfXloadOxSample;
        } else if (strequal_k_unsafe(flagname_p2, "heep")) {
          if (chr_info.chrset_source) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.autosome_ct = 26;
          chr_info.xymt_codes[0] = 27;
          chr_info.xymt_codes[1] = 28;
          chr_info.xymt_codes[2] = UINT32_MAXM1;
          chr_info.xymt_codes[3] = UINT32_MAXM1;
          chr_info.xymt_codes[4] = UINT32_MAXM1;
          chr_info.xymt_codes[5] = UINT32_MAXM1;
          chr_info.haploid_mask[0] = 0x18000000;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "np")) {
          if (pc.varid_exclude_snp) {
            // problematic due to --window
            logerrputs("Error: --snp cannot be used with --exclude-snp.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_snp);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "nps")) {
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.snps_range_list);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "nps-only")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (!strcmp(cur_modif, "just-acgt")) {
              pc.filter_flags |= kfFilterSnpsOnlyJustAcgt;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --snps-only parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.filter_flags |= kfFilterPvarReq | kfFilterSnpsOnly;
        } else if (strequal_k_unsafe(flagname_p2, "core")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 11)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &pc.score_info.input_fname);
          if (reterr) {
            goto main_ret_1;
          }
          uint32_t numeric_param_ct = 0;
          uint32_t score_cols[3];
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "header", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreHeaderIgnore;
            } else if (strequal_k(cur_modif, "header-read", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreHeaderRead;
            } else if (strequal_k(cur_modif, "no-mean-imputation", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreNoMeanimpute;
            } else if (strequal_k(cur_modif, "dominant", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreDominant;
            } else if (strequal_k(cur_modif, "recessive", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreRecessive;
            } else if (strequal_k(cur_modif, "center", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreCenter;
            } else if (strequal_k(cur_modif, "variance-standardize", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreVarianceStandardize;
            } else if (strequal_k(cur_modif, "variance-normalize", cur_modif_slen)) {
              logerrputs("Note: --score's 'variance-normalize' modifier has been renamed to the more\nprecise 'variance-standardize'.\n");
              pc.score_info.flags |= kfScoreVarianceStandardize;
            } else if (strequal_k(cur_modif, "se", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreSe;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreZs;
            } else if (strequal_k(cur_modif, "list-variants", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreListVariants;
            } else if (strequal_k(cur_modif, "list-variants-zs", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreListVariants | kfScoreListVariantsZs;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (pc.score_info.flags & kfScoreColAll) {
                logerrputs("Error: Multiple --score cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0pheno1\0phenos\0nmissallele\0denom\0dosagesum\0scoreavgs\0scoresums\0", "score", kfScoreColMaybefid, kfScoreColDefault, 1, &pc.score_info.flags);
              if (reterr) {
                goto main_ret_1;
              }
            } else if (numeric_param_ct == 3) {
              logerrputs("Error: --score takes at most three numeric parameters.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else {
              if (ScanPosintCapped(cur_modif, kMaxLongLine / 2, &(score_cols[numeric_param_ct]))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --score parameter '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              for (uint32_t uii = 0; uii < numeric_param_ct; ++uii) {
                if (score_cols[uii] == score_cols[numeric_param_ct]) {
                  logerrputs("Error: Identical --score column indexes.\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
              }
              ++numeric_param_ct;
            }
          }
          if ((pc.score_info.flags & (kfScoreHeaderIgnore | kfScoreHeaderRead)) == (kfScoreHeaderIgnore | kfScoreHeaderRead)) {
            logerrputs("Error: --score 'header' and 'header-read' modifiers cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          uint32_t model_flags_u = S_CAST(uint32_t, pc.score_info.flags & (kfScoreDominant | kfScoreRecessive | kfScoreCenter | kfScoreVarianceStandardize));
          if (model_flags_u & (model_flags_u - 1)) {
            logerrputs("Error: --score 'dominant', 'recessive', 'center', and 'variance-standardize'\nmodifiers are mutually exclusive.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (!(pc.score_info.flags & kfScoreColAll)) {
            pc.score_info.flags |= kfScoreColDefault;
          }
          if (numeric_param_ct) {
            pc.score_info.varid_col_p1 = score_cols[0];
          }
          if (numeric_param_ct > 1) {
            pc.score_info.allele_col_p1 = score_cols[1];
          } else {
            pc.score_info.allele_col_p1 = pc.score_info.varid_col_p1 + 1;
          }
          if (numeric_param_ct == 3) {
            // a bit artificial, but it works
            const uint32_t col_idx = score_cols[2];
            const uint32_t col_idx_blen = 1 + UintSlen(col_idx);
            char* new_buf;
            if (pgl_malloc(col_idx_blen + 1, &new_buf)) {
              goto main_ret_NOMEM;
            }
            pc.score_info.input_col_idx_range_list.names = new_buf;
            pc.score_info.input_col_idx_range_list.name_max_blen = col_idx_blen;
            pc.score_info.input_col_idx_range_list.name_ct = 1;
            u32toa_x(col_idx, '\0', new_buf);
            new_buf[col_idx_blen] = '\0';
            pc.score_info.input_col_idx_range_list.starts_range = R_CAST(unsigned char*, &(new_buf[col_idx_blen]));
          }
          pc.command_flags1 |= kfCommand1Score;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "core-col-nums")) {
          if (!(pc.command_flags1 & kfCommand1Score)) {
            logerrputs("Error: --score-col-nums must be used with --score.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.score_info.input_col_idx_range_list.name_ct) {
            logerrputs("Error: --score-col-nums cannot be used when three numeric parameters are\nprovided to --score.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.score_info.input_col_idx_range_list);
          if (reterr) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "plit-cat-pheno")) {
          uint32_t first_phenoname_idx = 1;
          for (; first_phenoname_idx <= param_ct; ++first_phenoname_idx) {
            const char* cur_modif = argvk[arg_idx + first_phenoname_idx];
            if (!strcmp(cur_modif, "omit-last")) {
              pc.pheno_transform_flags |= kfPhenoTransformSplitCatOmitLast;
            } else if (!strcmp(cur_modif, "covar-01")) {
              pc.pheno_transform_flags |= kfPhenoTransformSplitCatCovar01;
            } else {
              break;
            }
          }
          if (first_phenoname_idx <= param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + first_phenoname_idx]), param_ct + 1 - first_phenoname_idx, kMaxIdSlen - 1, &pc.split_cat_phenonames_flattened);
            if (reterr) {
              goto main_ret_1;
            }
            // may as well verify that no phenotype name has an '=' in it
            const char* phenonames_iter = pc.split_cat_phenonames_flattened;
            do {
              const uint32_t cur_phenoname_slen = strlen(phenonames_iter);
              if (memchr(phenonames_iter, '=', cur_phenoname_slen)) {
                logerrputs("Error: --split-cat-pheno phenotype names may not contain the '=' character.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              phenonames_iter = &(phenonames_iter[cur_phenoname_slen + 1]);
            } while (*phenonames_iter);
          } else if (pc.pheno_transform_flags & kfPhenoTransformSplitCatCovar01) {
            logerrputs("Error: --split-cat-pheno 'covar-01' modifier cannot be used without any\nphenotype names.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.pheno_transform_flags |= kfPhenoTransformSplitCat;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ort-vars")) {
          if (!(pc.command_flags1 & kfCommand1MakePlink2)) {
            // todo: permit merge
            logerrputs("Error: --sort-vars must be used with --make-{b}pgen/--make-bed or dataset\nmerging.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* sort_vars_mode_str = argvk[arg_idx + 1];
            const char first_char_upcase_match = sort_vars_mode_str[0] & 0xdf;
            const uint32_t is_short_name = (sort_vars_mode_str[1] == '\0');
            if ((is_short_name && (first_char_upcase_match == 'N')) || (!strcmp(sort_vars_mode_str, "natural"))) {
              pc.sort_vars_flags = kfSortNatural;
            } else if ((is_short_name && (first_char_upcase_match == 'A')) || (!strcmp(sort_vars_mode_str, "ascii"))) {
              pc.sort_vars_flags = kfSortAscii;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a valid mode for --sort-vars.\n", sort_vars_mode_str);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.sort_vars_flags = kfSortNatural;
          }
        } else if (strequal_k_unsafe(flagname_p2, "trict-sid0")) {
          pc.misc_flags |= kfMiscStrictSid0;
          goto main_param_zero;
        } else if (!strequal_k_unsafe(flagname_p2, "ilent")) {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 't':
        if (strequal_k_unsafe(flagname_p2, "hreads")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (ScanPosintDefcap(argvk[arg_idx + 1], &pc.max_thread_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --threads parameter '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (pc.max_thread_ct > kMaxThreads) {
            logprintf("Note: Reducing --threads parameter to %u.  (If this is not large enough,\nrecompile with a larger kMaxThreads setting.)\n", kMaxThreads);
            pc.max_thread_ct = kMaxThreads;
          } else if (known_procs == -1) {
            // trigger BLAS/LAPACK warning?
            known_procs = 0;
          }
        } else if (strequal_k_unsafe(flagname_p2, "o")) {
          if (chr_info.is_include_stack || notchr_present) {
            logerrputs("Error: --from/--to cannot be used with --autosome{-par} or --{not-}chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_to);
          if (reterr) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "o-bp") || strequal_k_unsafe(flagname_p2, "o-kb") || strequal_k_unsafe(flagname_p2, "o-mb")) {
          if (!CmdlineSingleChr(&chr_info, pc.misc_flags)) {
            logerrputs("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb must be used with --chr, and only\none chromosome.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (notchr_present) {
            logerrputs("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb cannot be used with --not-chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (pc.to_bp != -1) {
            logerrputs("Error: Multiple --to-bp/-kb/-mb values.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (!ScanadvDouble(cur_modif, &dxx)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --to-bp/-kb/-mb parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          const char unit_char = flagname_p2[2];
          if (unit_char == 'k') {
            dxx *= 1000;
          } else if (unit_char == 'm') {
            dxx *= 1000000;
          }
          if (dxx < 0) {
            logerrprintf("Error: --to-bp/-kb/-mb parameter '%s' too small.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_A;
          } else if (dxx >= 2147483646) {
            pc.to_bp = 0x7ffffffe;
          } else {
            // round down
            pc.to_bp = S_CAST(int32_t, dxx * (1 + kSmallEpsilon));
          }
          if (pc.from_bp > pc.to_bp) {
            // (if we do permit this, rounding must be postponed)
            logerrputs("Error: --to-bp/-kb/-mb parameter is smaller than --from-bp/-kb/-mb parameter.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "hin")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (!ScanadvDouble(cur_modif, &pc.thin_keep_prob)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --thin variant retention probability '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (pc.thin_keep_prob < (0.5 / 4294967296.0)) {
            logerrputs("Error: --thin variant retention probability too small.\n");
            goto main_ret_INVALID_CMDLINE_A;
          } else if (pc.thin_keep_prob >= (4294967295.5 / 4294967296.0)) {
            uint32_t uii;
            if (ScanUintDefcap(cur_modif, &uii)) {
              logerrputs("Error: --thin variant retention probability too large.\n");
            } else {
              // VCFtools --thin = --bp-space...
              logerrputs("Error: --thin variant retention probability too large.  (Did you mean\n--bp-space?)\n");
            }
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "hin-count")) {
          if (pc.thin_keep_prob != 1.0) {
            logerrputs("Error: --thin cannot be used with --thin-count.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (ScanUintDefcap(cur_modif, &pc.thin_keep_ct) || (!pc.thin_keep_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --thin-count parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "hin-indiv")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (!ScanadvDouble(cur_modif, &pc.thin_keep_sample_prob)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --thin-indiv sample retention probability '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (pc.thin_keep_sample_prob < (0.5 / 4294967296.0)) {
            logerrputs("Error: --thin-indiv sample retention probability too small.\n");
            goto main_ret_INVALID_CMDLINE_A;
          } else if (pc.thin_keep_sample_prob >= (4294967295.5 / 4294967296.0)) {
            logerrputs("Error: --thin-indiv sample retention probability too large.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "hin-indiv-count")) {
          if (pc.thin_keep_sample_prob != 1.0) {
            logerrputs("Error: --thin-indiv cannot be used with --thin-indiv-count.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (ScanUintDefcap(cur_modif, &pc.thin_keep_sample_ct) || (!pc.thin_keep_sample_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --thin-indiv-count parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ests")) {
          if (!(pc.command_flags1 & kfCommand1Glm)) {
            logerrputs("Error: --tests must be used with --glm.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if ((param_ct == 1) && (!strcmp(argvk[arg_idx + 1], "all"))) {
            pc.glm_info.flags |= kfGlmTestsAll;
          } else {
            reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.glm_info.tests_range_list);
            if (reterr) {
              goto main_ret_1;
            }
          }
          logerrputs("Error: --tests is not implemented yet.\n");
          reterr = kPglRetNotYetSupported;
          goto main_ret_1;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'u':
        if (strequal_k_unsafe(flagname_p2, "pdate-sex")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, 0, &pc.update_sex_info.fname);
          if (reterr) {
            goto main_ret_1;
          }
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "male0", cur_modif_slen)) {
              pc.update_sex_info.flags |= kfUpdateSexMale0;
            } else if (StrStartsWith(cur_modif, "col-num=", cur_modif_slen)) {
              const char* col_num_start = &(cur_modif[strlen("col-num=")]);
              if (ScanPosintDefcap(col_num_start, &pc.update_sex_info.col_num) || (pc.update_sex_info.col_num == 1)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --update-sex col-num= parameter '%s'.\n", col_num_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (param_ct == 2) {
              // only one extra parameter, try to interpret it the plink 1.9
              // way but print a warning
              if (ScanPosintDefcap(cur_modif, &pc.update_sex_info.col_num)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --update-sex parameter '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              logerrputs("Warning: --update-sex unlabeled column parameter is now deprecated.  Use\n'col-num=' instead (and add 2 to the value).\n");
              pc.update_sex_info.col_num += 2;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --update-sex parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "pdate-name")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = Alloc2col(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.update_name_flag);
          if (reterr) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'v':
        if (strequal_k_unsafe(flagname_p2, "ar-min-qual")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (ScanFloat(argvk[arg_idx + 1], &pc.var_min_qual) || (pc.var_min_qual < 0.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --var-min-qual parameter '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.var_min_qual *= 1 - kSmallEpsilon;
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "ar-filter")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &pc.var_filter_exceptions_flattened);
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.misc_flags |= kfMiscExcludePvarFilterFail;
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "cf")) {
          // permit accompanying .fam/.psam
          // IIDs must match VCF sample line order
          if ((load_params & (~kfLoadParamsPsam)) || xload) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 2) {
            const char* cur_modif = argvk[arg_idx + 2];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (!StrStartsWith(cur_modif, "dosage=", cur_modif_slen)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --vcf parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            reterr = CmdlineAllocString(&(cur_modif[strlen("dosage=")]), argvk[arg_idx], 4095, &vcf_dosage_import_field);
            if (reterr) {
              goto main_ret_1;
            }
            if (!IsAlphanumeric(vcf_dosage_import_field)) {
              logerrputs("Error: --vcf dosage= parameter is not alphanumeric.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            if (!strcmp(vcf_dosage_import_field, "GT")) {
              logerrputs("Error: --vcf dosage= parameter cannot be 'GT'.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_modif);
          if (slen > kPglFnamesize - 1) {
            logerrputs("Error: --vcf filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_modif, slen + 1);
          xload = kfXloadVcf;
        } else if (strequal_k_unsafe(flagname_p2, "cf-min-gp")) {
          logerrputs("Error: --vcf-min-gp is no longer supported.  Use --import-dosage-certainty\ninstead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "cf-min-gq") || strequal_k_unsafe(flagname_p2, "cf-min-dp")) {
          if (!(xload & kfXloadVcf)) {
            // todo: support BCF too
            logerrprintf("Error: --%s must be used with --vcf.\n", flagname_p);
            goto main_ret_INVALID_CMDLINE;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t uii;
          if (ScanUintDefcap(cur_modif, &uii)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s parameter '%s'.\n", flagname_p, cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (flagname_p2[7] == 'g') {
            vcf_min_gq = uii;
          } else {
            vcf_min_dp = uii;
          }
        } else if (strequal_k_unsafe(flagname_p2, "cf-idspace-to")) {
          if (!(xload & (kfXloadVcf | kfXloadBcf))) {
            logerrputs("Error: --vcf-idspace-to must be used with --vcf/--bcf.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (id_delim == ' ') {
            logerrputs("Error: --vcf-idspace-to cannot be used when the --id-delim character is space.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          idspace_to = ExtractCharParam(argvk[arg_idx + 1]);
          if (!idspace_to) {
            logerrputs("Error: --vcf-idspace-to parameter must be a single character.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (ctou32(idspace_to) <= ' ') {
            logerrputs("Error: --vcf-idspace-to parameter must be a nonspace character.\n");
            goto main_ret_INVALID_CMDLINE;
          }
        } else if (strequal_k_unsafe(flagname_p2, "cf-half-call")) {
          if (!(xload & kfXloadVcf)) {
            logerrputs("Error: --vcf-half-call must be used with --vcf.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* half_call_mode_str = argvk[arg_idx + 1];
          const char first_char_upcase_match = half_call_mode_str[0] & 0xdf;
          const uint32_t is_short_name = (half_call_mode_str[1] == '\0');
          if ((is_short_name && (first_char_upcase_match == 'H')) || (!strcmp(half_call_mode_str, "haploid"))) {
            vcf_half_call = kVcfHalfCallHaploid;
          } else if ((is_short_name && (first_char_upcase_match == 'M')) || (!strcmp(half_call_mode_str, "missing"))) {
            vcf_half_call = kVcfHalfCallMissing;
          } else if ((is_short_name && (first_char_upcase_match == 'E')) || (!strcmp(half_call_mode_str, "error"))) {
            vcf_half_call = kVcfHalfCallError;
          } else if ((is_short_name && (first_char_upcase_match == 'R')) || (!strcmp(half_call_mode_str, "reference"))) {
            vcf_half_call = kVcfHalfCallError;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a valid mode for --vcf-half-call.\n", half_call_mode_str);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "cf-require-gt")) {
          if (!(xload & (kfXloadVcf | kfXloadBcf))) {
            logerrputs("Error: --vcf-require-gt must be used with --vcf/--bcf.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          import_flags |= kfImportVcfRequireGt;
          pc.dependency_flags |= kfFilterAllReq;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "if")) {
          if (!(pc.command_flags1 & kfCommand1Glm)) {
            logerrputs("Error: --vif must be used with --glm/--epistasis.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (!ScanadvDouble(cur_modif, &pc.vif_thresh)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --glm/--epistasis VIF threshold '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (pc.vif_thresh < 1.0) {
            snprintf(g_logbuf, kLogbufSize, "Error: --glm/--epistasis VIF threshold '%s' too small (must be >= 1).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ariance-standardize")) {
          if (pc.pheno_transform_flags & kfPhenoTransformVstdCovar) {
            logerrputs("Error: --variance-standardize cannot be used with --covar-variance-standardize.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), param_ct, 0x7fffffff, &pc.vstd_flattened);
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.pheno_transform_flags |= kfPhenoTransformVstdAll;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "alidate")) {
          pc.command_flags1 |= kfCommand1Validate;
          pc.dependency_flags |= kfFilterAllReq;
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'w':
        if (strequal_k_unsafe(flagname_p2, "rite-snplist")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (strcmp(cur_modif, "zs")) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --write-snplist parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            pc.misc_flags |= kfMiscWriteSnplistZs;
          }
          pc.command_flags1 |= kfCommand1WriteSnplist;
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "rite-samples")) {
          if ((param_ct == 1) && (!strcmp(argvk[arg_idx + 1], "noheader"))) {
            logerrputs("Error: --write-samples 'noheader' modifier retired.  Use --no-id-header\ninstead.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.command_flags1 |= kfCommand1WriteSamples;
          pc.dependency_flags |= kfFilterPsamReq;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "indow")) {
          if (!(pc.varid_snp || pc.varid_exclude_snp)) {
            logerrputs("Error: --window must be used with --snp or --exclude-snp.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (!ScanadvDouble(cur_modif, &dxx) || (dxx < 0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --window parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          dxx *= 500 * (1 + kSmallEpsilon);
          if (dxx > 2147483646) {
            pc.window_bp = 0x7ffffffe;
          } else {
            pc.window_bp = S_CAST(int32_t, dxx);
          }
          pc.filter_flags |= kfFilterNoSplitChr;
          // no need to set kfFilterPvarReq due to --snp/--exclude-snp req.
        } else if (strequal_k_unsafe(flagname_p2, "ithin")) {
          if (pc.misc_flags & kfMiscCatPhenoFamily) {
            logerrputs("Error: --within cannot be used with --family.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if ((cur_modif_slen == 7) && StrStartsWithUnsafe(cur_modif, "keep-") && MatchUpperK(&(cur_modif[5]), "NA")) {
              logerrputs("Error: --within's keep-NA modifier has been retired.  Rename that category in\nthe input file if you wish to keep it.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (param_idx == 1) {
              reterr = AllocFname(cur_modif, flagname_p, 0, &pc.within_fname);
            } else {
              if (IsReservedPhenoName(cur_modif, cur_modif_slen)) {
                snprintf(g_logbuf, kLogbufSize, "Error: '%s' cannot be used as a categorical phenotype name.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              reterr = CmdlineAllocString(cur_modif, argvk[arg_idx], kMaxIdSlen, &pc.catpheno_name);
            }
            if (reterr) {
              goto main_ret_1;
            }
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "rite-covar")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (StrStartsWith0(cur_modif, "cols=", cur_modif_slen)) {
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0maybeparents\0parents\0sex\0pheno1\0phenos\0", "write-covar", kfWriteCovarColMaybefid, kfWriteCovarColDefault, 0, &pc.write_covar_flags);
              if (reterr) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --write-covar parameter '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.write_covar_flags |= kfWriteCovarColDefault;
          }
          pc.command_flags1 |= kfCommand1WriteCovar;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "arning-errcode")) {
          warning_errcode = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "rite-cluster")) {
          logerrputs("Error: --write-cluster is retired.  Use e.g. --make-just-psam.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "ith-phenotype")) {
          logerrputs("Error: --with-phenotype is retired.  Use --write-covar cols=... instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'x':
        if (strequal_k_unsafe(flagname_p2, "chr-model")) {
          if (!(pc.command_flags1 & (kfCommand1Glm | kfCommand1Score))) {
            logerrputs("Error: --xchr-model must be used with --glm or --score.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.glm_info.flags & (kfGlmGenotypic | kfGlmHethom | kfGlmDominant | kfGlmRecessive)) {
            snprintf(g_logbuf, kLogbufSize, "Error: --xchr-model cannot be used with --glm %s.\n", (pc.glm_info.flags & kfGlmGenotypic)? "genotypic" : ((pc.glm_info.flags & kfGlmHethom)? "hethom" : ((pc.glm_info.flags & kfGlmDominant)? "dominant" : "recessive")));
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          pc.xchr_model = ctou32(ExtractCharParam(cur_modif)) - 48;
          if (pc.xchr_model > 2) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --xchr-model parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'z':
        if (strequal_k_unsafe(flagname_p2, "st-level")) {
          if (EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1)) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (ScanPosintCapped(cur_modif, 22, &g_zst_level)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --zst-level parameter '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      default:
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      main_param_zero:
        if (param_ct) {
          snprintf(g_logbuf, kLogbufSize, "Error: --%s doesn't accept parameters.\n", flagname_p);
          goto main_ret_INVALID_CMDLINE_2A;
        }
      }
    } while ((++cur_flag_idx) < flag_ct);
    if (!outname_end) {
      outname_end = &(outname[6]);
    }

    pc.dependency_flags |= pc.filter_flags;
    const uint32_t skip_main = (!pc.command_flags1) && (!(xload & (kfXloadVcf | kfXloadBcf | kfXloadOxBgen | kfXloadOxHaps | kfXloadOxSample | kfXloadPlink1Dosage | kfXloadGenDummy)));
    const uint32_t batch_job = (adjust_file_info.fname != nullptr);
    if (skip_main && (!batch_job)) {
      // add command_flags2 when needed
      goto main_ret_NULL_CALC;
    }
    if (!(load_params || xload || batch_job)) {
      logerrputs("Error: No input dataset.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((xload & kfXloadOxGen) && (!(xload & kfXloadOxSample))) {
      // could permit .fam/.psam, but unless Oxford software supports that mode
      // it's pointless
      logerrputs("Error: --gen must be used with --sample or --data.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((xload & kfXloadOxSample) && (pc.misc_flags & kfMiscAffection01)) {
      // necessary for --data and --data --make-pgen to yield the same output
      logerrputs("Error: --data/--sample cannot be used with --1.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((pc.sample_sort_flags != kfSort0) && (!(pc.command_flags1 & (kfCommand1MakePlink2 | kfCommand1WriteCovar)))) {
      // todo: permit merge
      logerrputs("Error: --indiv-sort must be used with --make-{b}pgen/--make-bed/--write-covar\nor dataset merging.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((make_plink2_flags & (kfMakePlink2MMask | kfMakePlink2TrimAlts | kfMakePgenEraseAlt2Plus | kfMakePgenErasePhase | kfMakePgenEraseDosage)) && (pc.command_flags1 & (~kfCommand1MakePlink2))) {
      logerrputs("Error: When the 'multiallelics=', 'trim-alts', and/or 'erase-...' modifier is\npresent, --make-bed/--make-{b}pgen cannot be combined with other commands.\n(Other filters are fine.)\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (((pc.misc_flags & kfMiscMajRef) || pc.ref_allele_flag || pc.alt1_allele_flag || pc.ref_from_fa_fname) && (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Exportf)))) {
      logerrputs("Error: Flags which alter REF/ALT1 allele settings (--maj-ref, --ref-allele,\n--alt1-allele, --ref-from-fa) must be used with\n--make-bed/--make-{b}pgen/--export and no other commands.\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (pc.keep_cat_phenoname && (!pc.keep_cat_names_flattened) && (!pc.keep_cats_fname)) {
      logerrputs("Error: --keep-cat-pheno must be used with --keep-cats and/or --keep-cat-names.\n");
    }
    if (pc.remove_cat_phenoname && (!pc.remove_cat_names_flattened) && (!pc.remove_cats_fname)) {
      logerrputs("Error: --remove-cat-pheno must be used with --remove-cats and/or\n--remove-cat-names.\n");
    }
    if (aperm_present && (pc.command_flags1 & kfCommand1Glm) && (!(pc.glm_info.flags & kfGlmPerm))) {
      // If --aperm is present, at least one association analysis command which
      // supports adaptive permutation testing was also specified, but no
      // actual adaptive permutation test is happening, the user is likely to
      // be confused.  Produce a warning.  (Not an error since a sophisticated
      // user may want to use --script with different --aperm defaults.)
      logerrputs("Warning: --aperm only controls the settings for adaptive permutation tests; it\ndoes not cause such a test to be performed.  (Did you forget to add the 'perm'\nmodifier to an association analysis flag?)\n");
    }
    if ((pc.hard_call_thresh == UINT32_MAX) && (xload & (kfXloadVcf | kfXloadBcf | kfXloadOxGen | kfXloadOxBgen))) {
      if (pc.dosage_erase_thresh > (kDosageMid / 10)) {
        logerrputs("Error: --dosage-erase-threshold value cannot be larger than (default)\n--hard-call-threshold value.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    } else {
      if (pc.dosage_erase_thresh > pc.hard_call_thresh) {
        logerrputs("Error: --dosage-erase-threshold value cannot be larger than\n--hard-call-threshold value.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }
    if ((oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond)) == (kfOxfordImportRefFirst | kfOxfordImportRefSecond)) {
      logerrputs("Error: --data/--{b}gen 'ref-first' and 'ref-second' modifiers cannot be used\ntogether.\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (!strcmp(g_missing_catname, g_output_missing_pheno)) {
      logerrputs("Error: --missing-catname and --output-missing-phenotype strings can't match.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((pc.misc_flags & kfMiscChrOverrideCmdline) && (!chr_info.chrset_source)) {
      logerrputs("Error: --chr-override requires an explicit chromosome set.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((xload & kfXloadPlink1Dosage) && (!(load_params & kfLoadParamsPsam))) {
      logerrputs("Error: --import-dosage requires a .fam file.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((pc.misc_flags & kfMiscCovarColNums) && (!pc.covar_fname) && (!pc.pheno_fname)) {
      logerrputs("Error: --covar-col-nums requires --covar or --pheno.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((pc.grm_flags & kfGrmMatrixShapemask) && (pc.misc_flags & kfMiscNoIdHeader)) {
      pc.grm_flags |= kfGrmNoIdHeader;
    }
    if (!permit_multiple_inclusion_filters) {
      // Permit only one position- or ID-based variant inclusion filter, since
      // it's not immediately obvious whether the union or intersection should
      // be taken with multiple inclusion filters.
      // However, multiple exclusion filters are fine.  (Also,
      // --autosome{-par}/--chr is exempted since it's more obvious how they
      // interact with other filters.)
      const uint32_t inclusion_filter_extract = (pc.extract_fnames != nullptr);
      const uint32_t inclusion_filter_fromto_id = pc.varid_from || pc.varid_to;
      const uint32_t inclusion_filter_fromto_bp = (pc.from_bp != -1) || (pc.to_bp != -1);
      const uint32_t inclusion_filter_snpflag = (pc.varid_snp != nullptr);
      const uint32_t inclusion_filter_snpsflag = !!pc.snps_range_list.name_ct;
      if (inclusion_filter_extract + inclusion_filter_fromto_id + inclusion_filter_fromto_bp + inclusion_filter_snpflag + inclusion_filter_snpsflag > 1) {
        logerrputs("Error: Multiple variant inclusion filters specified (--extract, --from/--to,\n--from-bp/--to-bp, --snp, --snps).  Add --force-intersect if you really want\nthe intersection of these sets.  (If your variant IDs are unique, you can\nextract the union by e.g. running --write-snplist for each set, followed by\n--extract on all the .snplist files.)\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }

    if (!rseeds) {
      uint32_t seed = S_CAST(uint32_t, time(nullptr));
      snprintf(g_logbuf, kLogbufSize, "Random number seed: %u\n", seed);
      logputs_silent(g_logbuf);
      sfmt_init_gen_rand(&g_sfmt, seed);
    } else {
      if (rseed_ct == 1) {
        sfmt_init_gen_rand(&g_sfmt, rseeds[0]);
      } else {
        sfmt_init_by_array(&g_sfmt, rseeds, rseed_ct);
      }
      free(rseeds);
      rseeds = nullptr;
    }

    if (CmdlineParsePhase3(0, malloc_size_mb, memory_require, &pcm, &bigstack_ua)) {
      goto main_ret_NOMEM;
    }
    g_input_missing_geno_ptr = &(g_one_char_strs[2 * ctou32(input_missing_geno_char)]);
    g_output_missing_geno_ptr = &(g_one_char_strs[2 * ctou32(output_missing_geno_char)]);
    // pigz_init(pc.max_thread_ct);

    if (pc.max_thread_ct > 8) {
      logprintf("Using up to %u threads (change this with --threads).\n", pc.max_thread_ct);
    } else {
      // "1 compute thread" instead of "1 thread" since, when
      // max_thread_ct == 2, some code will use one I/O thread and one
      // compute thread.  Not worth the trouble of writing special-case code
      // to avoid that.  (also, with 2 cores, the I/O thread isn't
      // sufficiently busy to justify only 1 compute thread.)
      logprintf("Using %s%u compute thread%s.\n", (pc.max_thread_ct > 1)? "up to " : "", pc.max_thread_ct, (pc.max_thread_ct == 1)? "" : "s");
    }
    if (randmem) {
      reterr = RandomizeBigstack(pc.max_thread_ct);
      if (reterr) {
        goto main_ret_1;
      }
    }

    print_end_time = 1;

    if (batch_job) {
      if (adjust_file_info.fname) {
        reterr = AdjustFile(&adjust_file_info, pc.pfilter, pc.output_min_p, pc.max_thread_ct, outname, outname_end);
        if (reterr) {
          goto main_ret_1;
        }
      }
      if (skip_main) {
        goto main_ret_1;
      }
    }

    if (0) {
      // nonstandard cases (CNV, etc.) here
    } else {
      if (pc.dependency_flags && (!pc.command_flags1)) {
        logerrputs("Error: Basic file conversions do not support regular filter or transform\noperations.  Rerun your command with --make-bed/--make-{b}pgen.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if (xload) {
        char* convname_end = outname_end;
        if (pc.command_flags1) {
          if (import_flags & kfImportKeepAutoconv) {
            if (pc.misc_flags & kfMiscAffection01) {
              logerrputs("Error: --1 cannot be used with --keep-autoconv.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if ((output_missing_geno_char != '.') && (output_missing_geno_char != input_missing_geno_char)) {
              logerrputs("Error: --output-missing-genotype and --input-missing-genotype parameters cannot\nbe inconsistent when --keep-autoconv is specified.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            double dxx;
            const char* num_end = ScanadvDouble(g_output_missing_pheno, &dxx);
            if (num_end) {
              if (dxx != S_CAST(double, pc.missing_pheno)) {
                logerrputs("Error: --output-missing-phenotype and --input-missing-phenotype parameters\ncannot be inconsistent when --keep-autoconv is specified.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else if (!IsNanStr(g_output_missing_pheno, strlen(g_output_missing_pheno))) {
              logerrputs("Error: --output-missing-phenotype parameter must be numeric or 'NA' when\n--keep-autoconv is specified.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else {
            convname_end = Stpcpy(convname_end, "-temporary");
          }
        } else {
          import_flags |= kfImportKeepAutoconv;
        }
        const uint32_t convname_slen = convname_end - outname;
        uint32_t pgen_generated = 1;
        uint32_t psam_generated = 1;
        if (xload & kfXloadVcf) {
          const uint32_t no_samples_ok = !(pc.dependency_flags & (kfFilterAllReq | kfFilterPsamReq));
          if (no_samples_ok && (!(import_flags & kfImportKeepAutoconv)) && pc.command_flags1) {
            // special case: just treat the VCF as a .pvar file
            strcpy(pvarname, pgenname);
            pgenname[0] = '\0';
            goto main_reinterpret_vcf_instead_of_converting;
          } else {
            reterr = VcfToPgen(pgenname, (load_params & kfLoadParamsPsam)? psamname : nullptr, const_fid, vcf_dosage_import_field, pc.misc_flags, import_flags, no_samples_ok, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, id_delim, idspace_to, vcf_min_gq, vcf_min_dp, vcf_half_call, pc.fam_cols, pc.max_thread_ct, outname, convname_end, &chr_info, &pgen_generated, &psam_generated);
          }
        } else if (xload & kfXloadVcf) {
          logerrputs("Error: --bcf is not implemented yet.\n");
          reterr = kPglRetNotYetSupported;
        } else if (xload & kfXloadOxGen) {
          reterr = OxGenToPgen(pgenname, psamname, import_single_chr_str, ox_missing_code, pc.misc_flags, import_flags, oxford_import_flags, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, pc.max_thread_ct, outname, convname_end, &chr_info);
        } else if (xload & kfXloadOxBgen) {
          reterr = OxBgenToPgen(pgenname, psamname, const_fid, import_single_chr_str, ox_missing_code, pc.misc_flags, import_flags, oxford_import_flags, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, id_delim, idspace_to, pc.max_thread_ct, outname, convname_end, &chr_info);
        } else if (xload & kfXloadOxHaps) {
          reterr = OxHapslegendToPgen(pgenname, pvarname, psamname, import_single_chr_str, ox_missing_code, pc.misc_flags, import_flags, oxford_import_flags, outname, convname_end, &chr_info);
        } else if (xload & kfXloadPlink1Dosage) {
          reterr = Plink1DosageToPgen(pgenname, psamname, (xload & kfXloadMap)? pvarname : nullptr, import_single_chr_str, &plink1_dosage_info, pc.misc_flags, import_flags, pc.fam_cols, pc.missing_pheno, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, pc.max_thread_ct, outname, convname_end, &chr_info);
        } else if (xload & kfXloadGenDummy) {
          reterr = GenerateDummy(&gendummy_info, pc.misc_flags, import_flags, pc.hard_call_thresh, pc.dosage_erase_thresh, pc.max_thread_ct, outname, convname_end, &chr_info);
        }
        if (reterr || (!pc.command_flags1)) {
          goto main_ret_1;
        }

        // todo: we have to skip this when merging is involved
        pc.hard_call_thresh = UINT32_MAX;

        if (pgen_generated) {
          snprintf(memcpya(pgenname, outname, convname_slen), kMaxOutfnameExtBlen - 10, ".pgen");
        }
        snprintf(memcpya(pvarname, outname, convname_slen), kMaxOutfnameExtBlen - 10, ".pvar");
        if (psam_generated) {
          snprintf(memcpya(psamname, outname, convname_slen), kMaxOutfnameExtBlen - 10, ".psam");
        }
        if (!(import_flags & kfImportKeepAutoconv)) {
          if (pgen_generated) {
            if (PushLlStr(pgenname, &file_delete_list)) {
              goto main_ret_NOMEM;
            }
          }
          if (PushLlStr(pvarname, &file_delete_list)) {
            goto main_ret_NOMEM;
          }
          if (psam_generated) {
            if (PushLlStr(psamname, &file_delete_list)) {
              goto main_ret_NOMEM;
            }
          }
        }
        *outname_end = '\0';
      }
    main_reinterpret_vcf_instead_of_converting:
      if (pc.dependency_flags & kfFilterAllReq) {
        if ((!xload) && (load_params != kfLoadParamsPfileAll)) {
          logerrputs("Error: A full fileset (.pgen/.bed + .pvar/.bim + .psam/.fam) is required for\nthis.\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
      } else {
        // no genotype file required
        if (!(pc.dependency_flags & kfFilterNonrefFlagsNeeded)) {
          pgenname[0] = '\0';
        }

        if (pc.dependency_flags & kfFilterPvarReq) {
          if ((!xload) && (!(load_params & kfLoadParamsPvar))) {
            logerrputs("Error: A .pvar/.bim file is required for this.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else {
          pvarname[0] = '\0';
        }
        if (pc.dependency_flags & kfFilterPsamReq) {
          if ((!xload) && (!(load_params & kfLoadParamsPsam))) {
            logerrputs("Error: A .psam/.fam file is required for this.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else {
          psamname[0] = '\0';
        }
      }
      if ((pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Validate | kfCommand1WriteSnplist | kfCommand1WriteCovar | kfCommand1WriteSamples))) || ((pc.command_flags1 & kfCommand1MakePlink2) && (pc.sort_vars_flags == kfSort0))) {
        pc.dependency_flags |= kfFilterNoSplitChr;
      }

      BLAS_SET_NUM_THREADS(1);
      reterr = Plink2Core(&pc, make_plink2_flags, pgenname, psamname, pvarname, outname, outname_end, king_cutoff_fprefix, &chr_info);
    }
  }
  while (0) {
  main_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  main_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  main_ret_INVALID_CMDLINE_UNRECOGNIZED:
    InvalidArg(argv[arg_idx]);
    logerrputsb();
    logerrputs(errstr_append);
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_INVALID_CMDLINE_INPUT_CONFLICT:
    logerrprintf("Error: --%s conflicts with another input flag.\n%s", flagname_p, errstr_append);
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_INVALID_CMDLINE_WWA:
    WordWrapB(0);
  main_ret_INVALID_CMDLINE_2A:
    logerrputsb();
  main_ret_INVALID_CMDLINE_A:
    logerrputs(errstr_append);
  main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_NULL_CALC:
    if (pc.dependency_flags) {
      logerrputs("Warning: No output requested.  (Did you forget --make-bed/--make-{b}pgen?)\nExiting.\n");
    } else {
      logerrputs("Warning: No output requested.  Exiting.\n");
    }
  main_ret_NULL_CALC_0:
    fputs(kCmdlineFormatStr, stdout);
    fputs(notestr_null_calc2, stdout);
    reterr = kPglRetSkipped;
    break;
  }
 main_ret_1:
  DispExitMsg(reterr);
  while (0) {
  main_ret_NOMEM_NOLOG:
    PrintVer();
  main_ret_NOMEM_NOLOG2:
    fputs(kErrstrNomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
    }
    reterr = kPglRetNomem;
    break;
  }
 main_ret_NOLOG:
  free_cond(vcf_dosage_import_field);
  free_cond(ox_missing_code);
  free_cond(import_single_chr_str);
  free_cond(const_fid);
  free_cond(rseeds);
  CleanupPlink2CmdlineMeta(&pcm);
  CleanupAdjust(&adjust_file_info);
  free_cond(king_cutoff_fprefix);
  free_cond(pc.update_name_flag);
  free_cond(pc.alt1_allele_flag);
  free_cond(pc.ref_allele_flag);
  free_cond(pc.keep_fcol_name);
  free_cond(pc.keep_fcol_flattened);
  free_cond(pc.keep_fcol_fname);
  free_cond(pc.require_no_info_flattened);
  free_cond(pc.require_info_flattened);
  free_cond(pc.king_table_subset_fname);
  free_cond(pc.ref_from_fa_fname);
  free_cond(pc.loop_cats_phenoname);
  free_cond(pc.covar_quantnorm_flattened);
  free_cond(pc.quantnorm_flattened);
  free_cond(pc.vstd_flattened);
  free_cond(pc.require_covar_flattened);
  free_cond(pc.require_pheno_flattened);
  free_cond(pc.split_cat_phenonames_flattened);
  free_cond(pc.remove_cat_phenoname);
  free_cond(pc.remove_cat_names_flattened);
  free_cond(pc.remove_cats_fname);
  free_cond(pc.keep_cat_phenoname);
  free_cond(pc.keep_cat_names_flattened);
  free_cond(pc.keep_cats_fname);
  free_cond(pc.family_missing_catname);
  free_cond(pc.catpheno_name);
  free_cond(pc.within_fname);
  free_cond(pc.read_freq_fname);
  free_cond(pc.glm_local_covar_fname);
  free_cond(pc.glm_local_pvar_fname);
  free_cond(pc.glm_local_psam_fname);
  free_cond(pc.freq_alt1_binstr);
  free_cond(pc.freq_ref_binstr);
  free_cond(pc.removefam_fnames);
  free_cond(pc.remove_fnames);
  free_cond(pc.keepfam_fnames);
  free_cond(pc.keep_fnames);
  free_cond(pc.exclude_fnames);
  free_cond(pc.extract_fnames);
  free_cond(pc.sample_sort_fname);
  free_cond(pc.covar_fname);
  free_cond(pc.pheno_fname);
  free_cond(pc.varid_exclude_snp);
  free_cond(pc.varid_snp);
  free_cond(pc.varid_to);
  free_cond(pc.varid_from);
  free_cond(pc.missing_varid_match);
  free_cond(pc.varid_template);
  free_cond(pc.var_filter_exceptions_flattened);
  if (file_delete_list) {
    do {
      LlStr* llstr_ptr = file_delete_list->next;
      unlink(file_delete_list->str);
      free(file_delete_list);
      file_delete_list = llstr_ptr;
    } while (file_delete_list);
  }
  CleanupCmpExpr(&pc.exclude_if_info_expr);
  CleanupCmpExpr(&pc.extract_if_info_expr);
  CleanupCmpExpr(&pc.remove_if_expr);
  CleanupCmpExpr(&pc.keep_if_expr);
  CleanupScore(&pc.score_info);
  CleanupGlm(&pc.glm_info);
  CleanupChrInfo(&chr_info);
  CleanupLd(&pc.ld_info);
  CleanupUpdateSex(&pc.update_sex_info);
  CleanupRangeList(&pc.covar_range_list);
  CleanupRangeList(&pc.pheno_range_list);
  CleanupRangeList(&pc.exclude_snps_range_list);
  CleanupRangeList(&pc.snps_range_list);
  if (warning_errcode && g_stderr_written_to && (!reterr)) {
    logerrputs("--warning-errcode: One or more warnings in this run; exiting with code 61.\n");
    reterr = kPglRetWarningErrcode;
  }
  if (CleanupLogfile(print_end_time) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  if (bigstack_ua) {
    free(bigstack_ua);
  }
  return S_CAST(int32_t, reterr);
}
