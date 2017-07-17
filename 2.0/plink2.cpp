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


#include "plink2_data.h"
#include "plink2_decompress.h"
#include "plink2_filter.h"
#include "plink2_glm.h"
#include "plink2_ld.h"
#include "plink2_matrix_calc.h"
#include "plink2_misc.h"
#include "plink2_psam.h"
#include "plink2_pvar.h"
#include "plink2_random.h"
#include "plink2_set.h"

// #include <locale.h>
#include <time.h>
#include <unistd.h> // getcwd(), gethostname(), sysconf(), unlink()

#include "plink2_help.h"

#ifdef __cplusplus
namespace plink2 {
#endif
  
static const char ver_str[] = "PLINK v2.00a"
#ifdef NOLAPACK
  "NL"
#endif
#ifdef __LP64__
  #ifdef LAPACK_ILP64
    "LM"
  #endif
  #ifdef USE_SSE42
    #ifdef USE_AVX2
      #error "version string code needs to be updated"
    #endif
    " SSE4.2"
  #else
    " 64-bit"
  #endif
#else
  " 32-bit"
#endif

#ifdef USE_MKL
  " Intel"
#endif
  " (17 Jul 2017)";
static const char ver_str2[] =
  // include leading space if day < 10, so character length stays the same
  ""
#ifndef LAPACK_ILP64
  "  "
#endif
#ifndef USE_MKL
  "      "
#endif
#ifndef NOLAPACK
  "  "
#endif
  "    www.cog-genomics.org/plink/2.0/\n"
  "(C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3\n";
static const char errstr_append[] = "For more info, try '" PROG_NAME_STR " --help [flag name]' or '" PROG_NAME_STR " --help | more'.\n";

#ifndef NOLAPACK
static const char notestr_null_calc2[] = "Commands include --make-bpgen, --export, --freq, --geno-counts, --missing,\n--hardy, --indep-pairwise, --make-king, --king-cutoff, --write-snplist,\n--make-grm-gz, --pca, --glm, --score, --genotyping-rate, --validate, and\n--zst-decompress.\n\n'" PROG_NAME_STR " --help | more' describes all functions.\n";
#else
static const char notestr_null_calc2[] = "Commands include --make-bpgen, --export, --freq, --geno-counts, --missing,\n--hardy, --indep-pairwise, --make-king, --king-cutoff, --write-snplist,\n--make-grm-gz, --glm, --score, --genotyping-rate, --validate, and\n--zst-decompress.\n\n'" PROG_NAME_STR " --help | more' describes all functions.\n";
#endif

static const char errstr_nomem[] = "Error: Out of memory.  The --memory flag may be helpful.\n";
static const char errstr_write[] = "Error: File write failure.\n";
static const char errstr_read[] = "Error: File read failure.\n";
static const char errstr_thread_create[] = "Error: Failed to create thread.\n";

#ifndef __LP64__
  // 2047 seems to consistently fail on both OS X and Windows
  #ifdef _WIN32
CONSTU31(kMalloc32bitMbMax, 1760);
  #else
    #ifdef __APPLE__
CONSTU31(kMalloc32bitMbMax, 1920);
    #else
CONSTU31(kMalloc32bitMbMax, 2047);
    #endif
  #endif
#endif

// assumes logfile is open
void disp_exit_msg(pglerr_t reterr) {
  if (reterr) {
    if (reterr == kPglRetNomem) {
      logprint("\n");
      logerrprint(errstr_nomem);
      if (g_failed_alloc_attempt_size) {
	LOGERRPRINTF("Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
      }
    } else if (reterr == kPglRetReadFail) {
      logprint("\n");
      logerrprint(errstr_read);
    } else if (reterr == kPglRetWriteFail) {
      logprint("\n");
      logerrprint(errstr_write);
    } else if (reterr == kPglRetThreadCreateFail) {
      logprint("\n");
      logerrprint(errstr_thread_create);
    }
  }
}

// covar-variance-standardize + terminating null
CONSTU31(kMaxFlagBlen, 27);

FLAGSET_DEF_START()
  kfLoadParams0,
  kfLoadParamsPgen = (1 << 0),
  kfLoadParamsPsam = (1 << 1),
  kfLoadParamsPvar = (1 << 2),
  kfLoadParamsPfileAll = (kfLoadParamsPgen | kfLoadParamsPsam | kfLoadParamsPvar)
FLAGSET_DEF_END(load_params_t);

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
FLAGSET_DEF_END(xload_t);

// maximum number of usable cluster computers, this is arbitrary though it
// shouldn't be larger than 2^32 - 1
CONSTU31(kParallelMax, 32768);

uint32_t realpath_identical(const char* outname, const char* read_realpath, char* write_realpath_buf) {
#ifdef _WIN32
  const uint32_t fname_slen = GetFullPathName(outname, kPglFnamesize, write_realpath_buf, nullptr);
  return (fname_slen && (fname_slen <= kPglFnamesize) && (!strcmp(read_realpath, write_realpath_buf)));
#else
  return (realpath(outname, write_realpath_buf) && (!strcmp(read_realpath, write_realpath_buf)));
#endif
}


// assume for now that .pgen must always be accompanied by both .pvar and .psam
FLAGSET64_DEF_START()
  kfFilter0,
  kfFilterAllReq = (1 << 0),
  kfFilterPvarReq = (1 << 1),
  kfFilterPsamReq = (1 << 2),
  kfFilterNoSplitChr = (1 << 3),
  kfFilterExclFemales = (1 << 4),
  kfFilterExclMales = (1 << 5),
  kfFilterExclNosex = (1 << 6),
  kfFilterExclFounders = (1 << 7),
  kfFilterExclNonfounders = (1 << 8),
  kfFilterSnpsOnly = (1 << 9),
  kfFilterSnpsOnlyJustAcgt = (1 << 10)
FLAGSET64_DEF_END(filter_flags_t);

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
  kfCommand1WriteCovar = (1 << 16)
FLAGSET64_DEF_END(command1_flags_t);

// this is a hybrid, only kfSortFileSid is actually a flag
FLAGSET_DEF_START()
  kfSort0,
  kfSortNone = (1 << 0),
  kfSortNatural = (1 << 1),
  kfSortAscii = (1 << 2),
  kfSortFile = (1 << 3),
  kfSortFileSid = (1 << 4)
FLAGSET_DEF_END(sort_flags_t);

typedef struct plink2_cmdline_struct {
  misc_flags_t misc_flags;
  filter_flags_t filter_flags;
  command1_flags_t command_flags1;
  pvar_psam_t pvar_psam_modifier;
  exportf_flags_t exportf_modifier;
  sort_flags_t sample_sort_flags;
  grm_flags_t grm_flags;
  pca_flags_t pca_flags;
  write_covar_flags_t write_covar_flags;
  pheno_transform_flags_t pheno_transform_flags;
  range_list_t snps_range_list;
  range_list_t exclude_snps_range_list;
  range_list_t pheno_range_list;
  range_list_t covar_range_list;
  fam_col_t fam_cols;
  idpaste_t exportf_id_paste;
  ld_info_t ld_info;
  king_flags_t king_modifier;
  double king_cutoff;
  double king_table_filter;
  allele_freq_t allele_freq_modifier;
  missing_rpt_t missing_rpt_modifier;
  geno_counts_t geno_counts_modifier;
  hardy_flags_t hardy_modifier;
  glm_info_t glm_info;
  adjust_info_t adjust_info;
  score_info_t score_info;
  aperm_t aperm;
  cmp_expr_t keep_if_expr;
  cmp_expr_t remove_if_expr;
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
  char exportf_id_delim;
  
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
  char* update_sex_fname;
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
  char* vstd_flattened;
  char* quantnorm_flattened;
  char* covar_quantnorm_flattened;
} plink2_cmdline_t;

uint32_t is_single_variant_loader_needed(const char* king_cutoff_fprefix, command1_flags_t command_flags1, make_plink2_t make_plink2_modifier) {
  return (command_flags1 & (kfCommand1Exportf | kfCommand1MakeKing | kfCommand1GenoCounts | kfCommand1LdPrune | kfCommand1Validate | kfCommand1Pca | kfCommand1MakeRel | kfCommand1Glm | kfCommand1Score)) || ((command_flags1 & kfCommand1MakePlink2) && (make_plink2_modifier & kfMakePgen)) || ((command_flags1 & kfCommand1KingCutoff) && (!king_cutoff_fprefix));
}

uint32_t are_allele_freqs_needed(command1_flags_t command_flags1, double min_maf, double max_maf) {
  return (command_flags1 & (kfCommand1AlleleFreq | kfCommand1LdPrune | kfCommand1Pca | kfCommand1MakeRel | kfCommand1Score)) || (min_maf != 0.0) || (max_maf != 1.0);
}

uint32_t are_maj_alleles_needed(command1_flags_t command_flags1) {
  return (command_flags1 & (kfCommand1LdPrune | kfCommand1Pca | kfCommand1MakeRel));
}

uint32_t get_first_haploid_uidx(const chr_info_t* cip, unsorted_var_t vpos_sortstatus) {
  // returns 0x7fffffff if no X/haploid chromosomes present
  if (!(vpos_sortstatus & kfUnsortedVarSplitChr)) {
    const uint32_t chr_ct = cip->chr_ct;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx < chr_ct; ++chr_fo_idx) {
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if (IS_SET(cip->haploid_mask, chr_idx)) {
	return cip->chr_fo_vidx_start[chr_fo_idx];
      }
    }
  }
  return 0x7fffffff;
}

uint32_t are_allele_dosages_needed(misc_flags_t misc_flags, make_plink2_t make_plink2_modifier, uint32_t afreq_needed, uint64_t min_allele_dosage, uint64_t max_allele_dosage) {
  return (make_plink2_modifier & kfMakePlink2TrimAlts) || ((misc_flags & kfMiscNonfounders) && (afreq_needed || (misc_flags & kfMiscMajRef) || min_allele_dosage || (max_allele_dosage != 0xffffffffU)));
}

uint32_t are_founder_allele_dosages_needed(misc_flags_t misc_flags, uint32_t afreq_needed, uint64_t min_allele_dosage, uint64_t max_allele_dosage) {
  return (afreq_needed || (misc_flags & kfMiscMajRef) || min_allele_dosage || (max_allele_dosage != (~0LLU))) && (!(misc_flags & kfMiscNonfounders));
}

uint32_t are_sample_missing_dosage_cts_needed(misc_flags_t misc_flags, uint32_t smaj_missing_geno_report_requested, double mind_thresh, missing_rpt_t missing_rpt_modifier) {
  return ((mind_thresh != 1.0) && (misc_flags & kfMiscMindDosage)) || (smaj_missing_geno_report_requested && (missing_rpt_modifier & (kfMissingRptScolNmissDosage | kfMissingRptScolFmissDosage)));
}

uint32_t are_variant_missing_hc_cts_needed(command1_flags_t command_flags1, misc_flags_t misc_flags, double geno_thresh, missing_rpt_t missing_rpt_modifier) {
  return ((command_flags1 & kfCommand1GenotypingRate) && (!(misc_flags & kfMiscGenotypingRateDosage))) || ((command_flags1 & kfCommand1MissingReport) && (missing_rpt_modifier & (kfMissingRptVcolNmiss | kfMissingRptVcolNmissHh | kfMissingRptVcolHethap | kfMissingRptVcolFmiss | kfMissingRptVcolFmissHh | kfMissingRptVcolFhethap))) || ((geno_thresh != 1.0) && (!(misc_flags & kfMiscGenoDosage)));
}

uint32_t are_variant_hethap_cts_needed(command1_flags_t command_flags1, misc_flags_t misc_flags, double geno_thresh, missing_rpt_t missing_rpt_modifier, uint32_t first_hap_uidx) {
  return (first_hap_uidx != 0x7fffffff) && (((command_flags1 & kfCommand1MissingReport) && (missing_rpt_modifier & (kfMissingRptVcolNmissHh | kfMissingRptVcolHethap | kfMissingRptVcolFmissHh | kfMissingRptVcolFhethap))) || ((geno_thresh != 1.0) && (!(misc_flags & kfMiscGenoHhMissing))));
}

uint32_t are_variant_missing_dosage_cts_needed(command1_flags_t command_flags1, misc_flags_t misc_flags, double geno_thresh, missing_rpt_t missing_rpt_modifier) {
  return ((command_flags1 & kfCommand1GenotypingRate) && (misc_flags & kfMiscGenotypingRateDosage)) || ((command_flags1 & kfCommand1MissingReport) && (!(missing_rpt_modifier & kfMissingRptSampleOnly)) && (missing_rpt_modifier & (kfMissingRptVcolNmissDosage | kfMissingRptVcolFmissDosage))) || ((geno_thresh != 1.0) && (misc_flags & kfMiscGenoDosage));
}

// can simplify --geno-counts all-biallelic case, but let's first make sure the
// general case works for multiallelic variants
uint32_t are_raw_geno_cts_needed(command1_flags_t command_flags1, misc_flags_t misc_flags, double hwe_thresh) {
  return (command_flags1 & kfCommand1GenoCounts) || ((misc_flags & kfMiscNonfounders) && ((command_flags1 & kfCommand1Hardy) || (hwe_thresh != 1.0)));
}

uint32_t are_founder_raw_geno_cts_needed(command1_flags_t command_flags1, misc_flags_t misc_flags, double hwe_thresh) {
  return (!(misc_flags & kfMiscNonfounders)) && ((command_flags1 & kfCommand1Hardy) || (hwe_thresh != 1.0));
}

uint32_t is_info_reload_needed(command1_flags_t command_flags1, pvar_psam_t pvar_psam_modifier, exportf_flags_t exportf_modifier) {
  // add kfExportfBcf later
  return ((command_flags1 & kfCommand1MakePlink2) && (pvar_psam_modifier & kfPvarColXinfo)) || ((command_flags1 & kfCommand1Exportf) && (exportf_modifier & kfExportfVcf));
}

uint32_t grm_keep_needed(command1_flags_t command_flags1, pca_flags_t pca_flags) {
  return ((command_flags1 & kfCommand1Pca) && (!(pca_flags & kfPcaApprox)));
}

void report_genotyping_rate(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_missing_cts, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t variant_ct, uint32_t is_dosage) {
  // defined the same way as PLINK 1.x, to allow this to serve as a sanity
  // check
  // trivial to multithread this if it ever matters
  uint64_t tot_nony_missing = 0;
  uint64_t tot_y_missing = 0;
  uint64_t cur_tot_missing = 0;
  uint32_t y_start = 0xffffffffU;
  uint32_t y_end = 0xffffffffU;
  uint32_t variant_ct_y = 0;
  int32_t y_code;
  if (xymt_exists(cip, kChrOffsetY, &y_code)) {
    const uint32_t y_chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)y_code];
    y_start = cip->chr_fo_vidx_start[y_chr_fo_idx];
    y_end = cip->chr_fo_vidx_start[y_chr_fo_idx + 1];
    variant_ct_y = popcount_bit_idx(variant_include, y_start, y_end);
  }
  uint32_t y_thresh = y_start;
  uint32_t variant_uidx = 0;
  uint32_t is_y = 0;  
  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
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
	y_thresh = 0xffffffffU;
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
    LOGPRINTF("Total (%s) genotyping rate %sis exactly 1.\n", is_dosage? "dosage" : "hardcall", (raw_sample_ct != sample_ct)? "in remaining samples " : "");
    return;
  }
  double genotyping_rate;
  if (male_ct && variant_ct_y) {
    const uint64_t nony_possible_obs = (variant_ct - variant_ct_y) * ((uint64_t)sample_ct);
    const uint64_t y_possible_obs = variant_ct_y * ((uint64_t)male_ct);
    genotyping_rate = ((double)((int64_t)(nony_possible_obs - tot_nony_missing))) / ((double)((int32_t)sample_ct)) + ((double)((int64_t)(y_possible_obs - tot_y_missing))) / ((double)((int32_t)male_ct));
    genotyping_rate /= (int32_t)variant_ct;
  } else {
    variant_ct -= variant_ct_y;
    const uint64_t denom = variant_ct * ((uint64_t)sample_ct);
    genotyping_rate = (double)((int64_t)(denom - tot_nony_missing)) / ((double)((int64_t)denom));
  }
  if (genotyping_rate >= 0.9999995) {
    LOGPRINTF("Total (%s) genotyping rate %sis in [0.9999995, 1).\n", is_dosage? "dosage" : "hardcall", (raw_sample_ct != sample_ct)? "in remaining samples " : "");
  } else {
    LOGPRINTF("Total (%s) genotyping rate %sis %g.\n", is_dosage? "dosage" : "hardcall", (raw_sample_ct != sample_ct)? "in remaining samples " : "", genotyping_rate);
  }
}

pglerr_t apply_variant_bp_filters(const char* extract_fnames, const char* exclude_fnames, const chr_info_t* cip, const uint32_t* variant_bps, int32_t from_bp, int32_t to_bp, uint32_t raw_variant_ct, misc_flags_t misc_flags, unsorted_var_t vpos_sortstatus, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  // todo: add --from-bp/--to-bp
  if ((from_bp != -1) || (to_bp != -1)) {
    if (vpos_sortstatus & kfUnsortedVarBp) {
      logerrprint("Error: --from-bp and --to-bp require a sorted .pvar/.bim.  Retry this command\nafter using e.g. plink 1.9 --make-bed to sort your data.\n");
      return kPglRetInconsistentInput;
    }
    const uint32_t chr_idx = next_set(cip->chr_mask, 0, kChrRawEnd);

    // this function shouldn't be called unless variant_ct is nonzero
    assert(chr_idx != kChrRawEnd);

    const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
    uint32_t variant_uidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
    uint32_t variant_uidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    if (from_bp != -1) {
      const uint32_t from_offset = uint32arr_greater_than(&(variant_bps[variant_uidx_start]), variant_uidx_end - variant_uidx_start, (uint32_t)from_bp);
      variant_uidx_start += from_offset;
    }
    if ((to_bp != -1) && (variant_uidx_start < variant_uidx_end)) {
      const uint32_t to_offset = uint32arr_greater_than(&(variant_bps[variant_uidx_start]), variant_uidx_end - variant_uidx_start, 1 + ((uint32_t)to_bp));
      variant_uidx_end = variant_uidx_start + to_offset;
    }
    if (variant_uidx_start) {
      clear_bits_nz(0, variant_uidx_start, variant_include);
    }
    if (variant_uidx_end < raw_variant_ct) {
      clear_bits_nz(variant_uidx_end, raw_variant_ct, variant_include);
    }
    *variant_ct_ptr = popcount_bit_idx(variant_include, variant_uidx_start, variant_uidx_end);
  }
  if (extract_fnames && (misc_flags & kfMiscExtractRange)) {
    if (vpos_sortstatus & kfUnsortedVarBp) {
      logerrprint("Error: '--extract range' requires a sorted .pvar/.bim.  Retry this command\nafter using e.g. plink 1.9 --make-bed to sort your data.\n");
      return kPglRetInconsistentInput;
    }
    pglerr_t reterr = extract_exclude_range(extract_fnames, cip, variant_bps, raw_variant_ct, 0, variant_include, variant_ct_ptr);
    if (reterr) {
      return reterr;
    }
  }
  if (exclude_fnames && (misc_flags & kfMiscExcludeRange)) {
    if (vpos_sortstatus & kfUnsortedVarBp) {
      logerrprint("Error: '--exclude range' requires a sorted .pvar/.bim.  Retry this command\nafter using e.g. plink 1.9 --make-bed to sort your data.\n");
      return kPglRetInconsistentInput;
    }
    pglerr_t reterr = extract_exclude_range(exclude_fnames, cip, variant_bps, raw_variant_ct, 1, variant_include, variant_ct_ptr);
    if (reterr) {
      return reterr;
    }
  }
  return kPglRetSuccess;
}

void update_sample_subsets(const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t* founder_info, uint32_t* founder_ct_ptr, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t* male_ct_ptr, uint32_t* nosex_ct_ptr) {
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  bitvec_and(sample_include, raw_sample_ctl, founder_info);
  *founder_ct_ptr = popcount_longs(founder_info, raw_sample_ctl);
  bitvec_and(sample_include, raw_sample_ctl, sex_male);
  *male_ct_ptr = popcount_longs(sex_male, raw_sample_ctl);
  bitvec_and(sample_include, raw_sample_ctl, sex_nm);
  *nosex_ct_ptr = sample_ct - popcount_longs(sex_nm, raw_sample_ctl);
}

// command_flags2 will probably be needed before we're done
static_assert(kPglMaxAltAlleleCt == 254, "plink2() --maj-ref needs to be updated.");
pglerr_t plink2_core(char* var_filter_exceptions_flattened, char* require_pheno_flattened, char* require_covar_flattened, const plink2_cmdline_t* pcp, make_plink2_t make_plink2_modifier, char* pgenname, char* psamname, char* pvarname, char* outname, char* outname_end, char* king_cutoff_fprefix, chr_info_t* cip) {
  pheno_col_t* pheno_cols = nullptr;
  pheno_col_t* covar_cols = nullptr;
  char* pheno_names = nullptr;
  char* covar_names = nullptr;
  uint32_t pheno_ct = 0;
  uint32_t covar_ct = 0;
  pglerr_t reterr = kPglRetSuccess;
  pgen_file_info_t pgfi;
  pgen_reader_t simple_pgr;
  pgfi_preinit(&pgfi);
  pgr_preinit(&simple_pgr);
  {
    // this predicate will need to exclude --merge-list special case later
    uint32_t pvar_renamed = 0;
    if ((make_plink2_modifier & (kfMakeBed | kfMakePgen)) || (pcp->exportf_modifier & kfExportfIndMajorBed)) {
      uint32_t fname_slen;
#ifdef _WIN32
      fname_slen = GetFullPathName(pgenname, kPglFnamesize, g_textbuf, nullptr);
      if ((!fname_slen) || (fname_slen > kPglFnamesize))
#else
      if (!realpath(pgenname, g_textbuf))
#endif
      {
	LOGERRPRINTFWW(g_errstr_fopen, pgenname);
	goto plink2_ret_OPEN_FAIL;
      }
      uint32_t pgen_rename = 0;
      if (make_plink2_modifier & kfMakePgen) {
        strcpy(outname_end, ".pgen");
	pgen_rename = realpath_identical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
      }
      if ((!pgen_rename) && ((make_plink2_modifier & kfMakeBed) || (pcp->exportf_modifier & kfExportfIndMajorBed))) {
	strcpy(outname_end, ".bed");
	pgen_rename = realpath_identical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
      }
      if (pgen_rename) {
	LOGPRINTF("Note: --make-%s input and output filenames match.  Appending '~' to input\nfilenames.\n", (make_plink2_modifier & kfMakeBed)? "bed" : ((make_plink2_modifier & kfMakePvar)? "pgen" : "bpgen"));
	fname_slen = strlen(pgenname);
	memcpy(g_textbuf, pgenname, fname_slen);
	strcpy(&(pgenname[fname_slen]), "~");
	if (rename(g_textbuf, pgenname)) {
	  logerrprint("Error: Failed to append '~' to input .bed/.pgen filename.\n");
	  goto plink2_ret_OPEN_FAIL;
	}
	fname_slen = strlen(pvarname);
	memcpy(g_textbuf, pvarname, fname_slen);
	strcpy(&(pvarname[fname_slen]), "~");
	if (rename(g_textbuf, pvarname)) {
	  logerrprint("Error: Failed to append '~' to input .bim/.pvar filename.\n");
	  goto plink2_ret_OPEN_FAIL;
	}
	pvar_renamed = 1;
	fname_slen = strlen(psamname);
	memcpy(g_textbuf, psamname, fname_slen);
	strcpy(&(psamname[fname_slen]), "~");
	if (rename(g_textbuf, psamname)) {
	  logerrprint("Error: Failed to append '~' to input .fam/.psam filename.\n");
	  goto plink2_ret_OPEN_FAIL;
	}
      }
    }
    uintptr_t max_sample_id_blen = 4;
    uintptr_t max_sid_blen = 0;
    uintptr_t max_paternal_id_blen = 2;
    uintptr_t max_maternal_id_blen = 2;
    uint32_t raw_sample_ct = 0;
    uintptr_t* sample_include = nullptr;
    char* sample_ids = nullptr;
    char* sids = nullptr;
    char* paternal_ids = nullptr;
    char* maternal_ids = nullptr;
    uintptr_t* sex_nm = nullptr;
    uintptr_t* sex_male = nullptr;
    uintptr_t* founder_info = nullptr;
    uintptr_t max_pheno_name_blen = 0;
    uint32_t raw_sample_ctl = 0;
    uint32_t sample_ct = 0;
    if (psamname[0]) {
      reterr = load_psam(psamname, pcp->pheno_fname? nullptr : &(pcp->pheno_range_list), pcp->fam_cols, (pcp->pheno_fname && pcp->pheno_range_list.name_ct)? 0 : 0x7fffffff, pcp->missing_pheno, (pcp->misc_flags / kfMiscAffection01) & 1, &max_sample_id_blen, &max_sid_blen, &max_paternal_id_blen, &max_maternal_id_blen, &sample_include, &sample_ids, &sids, &paternal_ids, &maternal_ids, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
      if (reterr) {
	goto plink2_ret_1;
      }
      // todo: add option to discard loaded SIDs
      raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
      sample_ct = popcount_longs(sample_include, raw_sample_ctl);
      const uint32_t known_sex_ct = popcount_longs(sex_nm, raw_sample_ctl);
      const uint32_t male_ct = popcount_longs(sex_male, raw_sample_ctl);
      const uint32_t female_ct = known_sex_ct - male_ct;
      const uint32_t founder_ct = popcount_longs(founder_info, raw_sample_ctl);
      if (known_sex_ct == sample_ct) {
        LOGPRINTFWW("%u sample%s (%u female%s, %u male%s; %u founder%s) loaded from %s.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", founder_ct, (founder_ct == 1)? "" : "s", psamname);
      } else {
	const uint32_t unknown_sex_ct = sample_ct - known_sex_ct;
        LOGPRINTFWW("%u sample%s (%u female%s, %u male%s, %u ambiguous; %u founder%s) loaded from %s.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", unknown_sex_ct, founder_ct, (founder_ct == 1)? "" : "s", psamname);
      }
    }

    uint32_t max_variant_id_slen = 1;
    uint32_t info_reload_slen = is_info_reload_needed(pcp->command_flags1, pcp->pvar_psam_modifier, pcp->exportf_modifier);
    uintptr_t* variant_allele_idxs = nullptr;
    uint32_t raw_variant_ct = 0;
    uint32_t variant_ct = 0;
    char* xheader = nullptr;
    uintptr_t xheader_blen = 0;
    uintptr_t* variant_include = nullptr;
    uint32_t* variant_bps = nullptr;
    char** variant_ids = nullptr;
    char** allele_storage = nullptr;
    uintptr_t* pvar_qual_present = nullptr;
    float* pvar_quals = nullptr;
    uintptr_t* pvar_filter_present = nullptr;
    uintptr_t* pvar_filter_npass = nullptr;
    char** pvar_filter_storage = nullptr;
    uintptr_t* nonref_flags = nullptr;
    uint32_t xheader_info_pr = 0;
    uint32_t max_allele_slen = 0;
    uint32_t max_filter_slen = 0;
    unsorted_var_t vpos_sortstatus = kfUnsortedVar0;
    double* variant_cms = nullptr;
    chr_idx_t* chr_idxs = nullptr; // split-chromosome case only
    if (pvarname[0]) {
      reterr = load_pvar(pvarname, var_filter_exceptions_flattened, pcp->varid_template, pcp->missing_varid_match, pcp->misc_flags, pcp->pvar_psam_modifier, pcp->exportf_modifier, pcp->var_min_qual, pcp->splitpar_bound1, pcp->splitpar_bound2, pcp->new_variant_id_max_allele_slen, (pcp->filter_flags / kfFilterSnpsOnly) & 3, !(pcp->filter_flags & kfFilterNoSplitChr), cip, &max_variant_id_slen, &info_reload_slen, &vpos_sortstatus, &xheader, &variant_include, &variant_bps, &variant_ids, &variant_allele_idxs, &allele_storage, &pvar_qual_present, &pvar_quals, &pvar_filter_present, &pvar_filter_npass, &pvar_filter_storage, &nonref_flags, &variant_cms, &chr_idxs, &raw_variant_ct, &variant_ct, &max_allele_slen, &xheader_blen, &xheader_info_pr, &max_filter_slen);
      if (reterr) {
	goto plink2_ret_1;
      }
      if (variant_ct == raw_variant_ct) {
	LOGPRINTFWW("%u variant%s loaded from %s.\n", variant_ct, (variant_ct == 1)? "" : "s", pvarname);
      } else {
	LOGPRINTFWW("%u out of %u variant%s loaded from %s.\n", variant_ct, raw_variant_ct, (raw_variant_ct == 1)? "" : "s", pvarname);
      }
      if (info_reload_slen && (make_plink2_modifier & (kfMakeBim | kfMakePvar)) && (!pvar_renamed)) {
	// need to be careful with .pvar in this case
	uint32_t fname_slen;
#ifdef _WIN32
	fname_slen = GetFullPathName(pvarname, kPglFnamesize, g_textbuf, nullptr);
	if ((!fname_slen) || (fname_slen > kPglFnamesize))
#else
	if (!realpath(pvarname, g_textbuf))
#endif
	{
	  LOGERRPRINTFWW(g_errstr_fopen, pvarname);
	  goto plink2_ret_OPEN_FAIL;
	}
	if (make_plink2_modifier & kfMakeBim) {
	  char* bimname_end = strcpya0(outname_end, ".bim");
	  if (make_plink2_modifier & kfMakeBimZs) {
	    strcpy(bimname_end, ".zst");
	  }
	  pvar_renamed = realpath_identical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
	  if (pvar_renamed) {
	    logprint("Note: .bim input and output filenames match.  Appending '~' to input filename.\n");
	    fname_slen = strlen(pvarname);
	    memcpy(g_textbuf, pvarname, fname_slen);
	    strcpy(&(pvarname[fname_slen]), "~");
	    if (rename(g_textbuf, pvarname)) {
	      logerrprint("Error: Failed to append '~' to input .bim filename.\n");
	      goto plink2_ret_OPEN_FAIL;
	    }
	  }
	}
	if ((!pvar_renamed) && (make_plink2_modifier & kfMakePvar)) {
	  char* pvarname_end = strcpya0(outname_end, ".pvar");
	  if (pcp->pvar_psam_modifier & kfPvarZs) {
	    strcpy(pvarname_end, ".zst");
	  }
	  // pvar_renamed = realpath_identical();
	  if (realpath_identical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]))) {
	    logprint("Note: .pvar input and output filenames match.  Appending '~' to input filename.\n");
	    fname_slen = strlen(pvarname);
	    memcpy(g_textbuf, pvarname, fname_slen);
	    strcpy(&(pvarname[fname_slen]), "~");
	    if (rename(g_textbuf, pvarname)) {
	      logerrprint("Error: Failed to append '~' to input .pvar filename.\n");
	      goto plink2_ret_OPEN_FAIL;
	    }
	  }
	}
      }
    }

    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uintptr_t pgr_alloc_cacheline_ct = 0;
    if (pgenname[0]) {
      pgen_header_ctrl_t header_ctrl;
      uintptr_t cur_alloc_cacheline_ct;
      while (1) {
	reterr = pgfi_init_phase1(pgenname, raw_variant_ct, raw_sample_ct, 0, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, g_logbuf);
	if (!reterr) {
	  break;
	}
	// detect and autoconvert plink 1 sample-major files, instead of
	// failing (don't bother supporting plink 0.99 files any more)
	if (reterr == kPglRetSampleMajorBed) {
	  char* pgenname_end = memcpya(pgenname, outname, (uintptr_t)(outname_end - outname));
	  pgenname_end = strcpya(pgenname_end, ".pgen");
	  const uint32_t no_vmaj_ext = (pcp->command_flags1 & kfCommand1MakePlink2) && (!pcp->filter_flags) && ((make_plink2_modifier & (kfMakePgen | (kfMakePgenFormatBase * 3))) == kfMakePgen);
	  if (no_vmaj_ext) {
	    *pgenname_end = '\0';
	    make_plink2_modifier &= ~kfMakePgen;
	    // no --make-just-pgen command, so we'll never entirely skip the
	    // make_plink2 operation
	  } else {
	    strcpy(pgenname_end, ".vmaj");
	  }
	  reterr = plink1_sample_major_to_pgen(pgenname, raw_variant_ct, raw_sample_ct, (pcp->misc_flags / kfMiscRealRefAlleles) & 1, pcp->max_thread_ct, pgfi.shared_ff);
	  if (!reterr) {
	    fclose(pgfi.shared_ff);
	    pgfi.shared_ff = nullptr;
	    continue;
	  }
	} else {
	  if (reterr != kPglRetReadFail) {
	    wordwrapb(0);
	    logerrprintb();
	  }
	}
	goto plink2_ret_1;
      }
      pgfi.allele_idx_offsets = variant_allele_idxs;
      unsigned char* pgfi_alloc;
      if (bigstack_alloc_uc(cur_alloc_cacheline_ct * kCacheline, &pgfi_alloc)) {
	goto plink2_ret_NOMEM;
      }
      const uint32_t nonref_flags_already_loaded = (nonref_flags != nullptr);
      if ((!nonref_flags) && ((header_ctrl & 192) == 192)) {
	if (bigstack_alloc_ul(raw_variant_ctl, &nonref_flags)) {
	  goto plink2_ret_NOMEM;
	}
      }
      pgfi.nonref_flags = nonref_flags;
      uint32_t max_vrec_width;
      // only practical effect of setting use_blockload to zero here is that
      // pgr_alloc_cacheline_ct is overestimated by
      // DIV_UP(max_vrec_width, kCacheline).
      reterr = pgfi_init_phase2(header_ctrl, 1, nonref_flags_already_loaded, 1, 0, raw_variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, &pgr_alloc_cacheline_ct, g_logbuf);
      if (reterr) {
	if (reterr != kPglRetReadFail) {
	  wordwrapb(0);
	  logerrprintb();
	}
	goto plink2_ret_1;
      }
      if (pcp->misc_flags & kfMiscRealRefAlleles) {
	if (nonref_flags && (!are_all_bits_one(nonref_flags, raw_variant_ct))) {
	  // technically a lie, it's okay if a .bed is first converted to .pgen
	  // without this flag, and then the user remembers the existence of
	  // --real-ref-alleles later.  but to reduce the ease of
	  // foot-shooting, we don't allow this to clobber arbitrary
	  // nonref_flags arrays.
	  logerrprint("Error: --real-ref-alleles must be used on a plink1 fileset.\n");
	  goto plink2_ret_INCONSISTENT_INPUT;
	}

	// wasteful if nonref_flags was allocated, but nonref_flags isn't that
	// large, and --real-ref-alleles + --make-pgen can be run separately
	// from anything truly memory-limited
	nonref_flags = nullptr;
	pgfi.nonref_flags = nullptr;
	
	pgfi.gflags &= ~kfPgenGlobalAllNonref;
      }
      if (is_single_variant_loader_needed(king_cutoff_fprefix, pcp->command_flags1, make_plink2_modifier)) {
	// ugly kludge, probably want to add pgenlib_internal support for this
	// hybrid use pattern
	FILE* shared_ff_copy = pgfi.shared_ff;
	pgfi.shared_ff = nullptr;
	unsigned char* simple_pgr_alloc;
	if (bigstack_alloc_uc((pgr_alloc_cacheline_ct + DIV_UP(max_vrec_width, kCacheline)) * kCacheline, &simple_pgr_alloc)) {
	  goto plink2_ret_NOMEM;
	}
	reterr = pgr_init(pgenname, max_vrec_width, &pgfi, &simple_pgr, simple_pgr_alloc);
	if (reterr) {
	  if (reterr == kPglRetOpenFail) {
	    LOGERRPRINTF(g_errstr_fopen, pgenname);
	  }
	  // only other possibility is kPglRetReadFail
	  goto plink2_ret_1;
	}
	pgfi.shared_ff = shared_ff_copy;
	if (pcp->command_flags1 & kfCommand1Validate) {
	  LOGPRINTFWW5("Validating %s... ", pgenname);
	  fflush(stdout);
	  reterr = pgr_validate(&simple_pgr, g_logbuf);
	  if (reterr) {
	    if (reterr != kPglRetReadFail) {
	      logprint("\n");
	      wordwrapb(0);
	      logerrprintb();
	    }
	    goto plink2_ret_1;
	  }
	  logprint("done.\n");
	  if (!(pcp->command_flags1 & (~kfCommand1Validate))) {
	    goto plink2_ret_1;
	  }
	}
      }
      // any functions using blockload must perform its own pgr_init(), etc.
    }
    if (pcp->pheno_fname) {
      reterr = load_phenos(pcp->pheno_fname, &(pcp->pheno_range_list), sample_include, sample_ids, raw_sample_ct, sample_ct, max_sample_id_blen, pcp->missing_pheno, (pcp->misc_flags / kfMiscAffection01) & 1, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen);
      if (reterr) {
	goto plink2_ret_1;
      }
    }

    // move processing of PLINK 1.x cluster-loading/filtering flags here, since
    // they're now under the categorical-phenotype umbrella
    if ((pcp->misc_flags & kfMiscCatPhenoFamily) || pcp->within_fname) {
      reterr = plink1_cluster_import(pcp->within_fname, pcp->catpheno_name, pcp->family_missing_catname, sample_include, sample_ids, raw_sample_ct, sample_ct, max_sample_id_blen, pcp->mwithin_val, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen);
      if (reterr) {
	goto plink2_ret_1;
      }
    }

    if (!pheno_ct) {
      logprint("Note: No phenotype data present.\n");      
    } else {
      if (pheno_ct == 1) {
	if (pheno_cols[0].type_code == kPhenoDtypeCc) {
	  const uint32_t obs_ct = popcount_longs(pheno_cols[0].nonmiss, raw_sample_ctl);
	  const uint32_t case_ct = popcount_longs(pheno_cols[0].data.cc, raw_sample_ctl);
	  const uint32_t ctrl_ct = obs_ct - case_ct;
	  LOGPRINTF("1 binary phenotype loaded (%u case%s, %u control%s).\n", case_ct, (case_ct == 1)? "" : "s", ctrl_ct, (ctrl_ct == 1)? "" : "s");
	} else if (pheno_cols[0].type_code == kPhenoDtypeQt) {
	  LOGPRINTF("1 quantitative phenotype loaded.\n");
	} else {
	  LOGPRINTF("1 categorical phenotype loaded.\n");
	}
      } else {
	uint32_t cc_ct = 0;
	uint32_t qt_ct = 0;
	for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	  const pheno_dtype_t cur_type_code = pheno_cols[pheno_idx].type_code;
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
	    LOGPRINTF("%u categorical phenotypes loaded.\n", pheno_ct);
	  } else if (!cat_ct) {
	    LOGPRINTF("%u quantitative phenotypes loaded.\n", pheno_ct);
	  } else {
	    LOGPRINTF("%u phenotypes loaded (%u quantitative, %u categorical).\n", pheno_ct, qt_ct, cat_ct);
	  }
	} else if (!qt_ct) {
	  if (!cat_ct) {
	    LOGPRINTF("%u binary phenotypes loaded.\n", pheno_ct);
	  } else {
	    LOGPRINTF("%u phenotypes loaded (%u binary, %u categorical).\n", pheno_ct, cc_ct, cat_ct);
	  }
	} else if (!cat_ct) {
	  LOGPRINTF("%u phenotypes loaded (%u binary, %u quantitative).\n", pheno_ct, cc_ct, qt_ct);
	} else {
	  LOGPRINTFWW("%u phenotypes loaded (%u binary, %u quantitative, %u categorical).\n", pheno_ct, cc_ct, qt_ct, cat_ct);
	}
      }
    }
    const uint32_t full_variant_id_htable_needed = variant_ct && (pcp->varid_from || pcp->varid_to || pcp->varid_snp || pcp->varid_exclude_snp || pcp->snps_range_list.name_ct || pcp->exclude_snps_range_list.name_ct);
    if (variant_ct && (!full_variant_id_htable_needed)) {
      reterr = apply_variant_bp_filters(pcp->extract_fnames, pcp->exclude_fnames, cip, variant_bps, pcp->from_bp, pcp->to_bp, raw_variant_ct, pcp->misc_flags, vpos_sortstatus, variant_include, &variant_ct);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
    if (variant_ct && (full_variant_id_htable_needed || (pcp->extract_fnames && (!(pcp->misc_flags & kfMiscExtractRange))) || (pcp->exclude_fnames && (!(pcp->misc_flags & kfMiscExcludeRange))))) {
      // don't bother with having different allow_dups vs. no allow_dups hash
      // table structures, just check specific IDs for duplication in the
      // no-duplicates-allowed cases
      unsigned char* bigstack_mark = g_bigstack_base;
      uint32_t* variant_id_htable = nullptr;
      uint32_t* htable_dup_base = nullptr;
      uint32_t variant_id_htable_size;
      reterr = alloc_and_populate_id_htable_mt(variant_include, variant_ids, variant_ct, pcp->max_thread_ct, &variant_id_htable, &htable_dup_base, &variant_id_htable_size);
      if (reterr) {
	goto plink2_ret_1;
      }
      if (vpos_sortstatus & kfUnsortedVarBp) {
	if (pcp->varid_from || pcp->varid_to) {
	  logerrprint("Error: --from/--to require a sorted .pvar/.bim.  Retry this command after using\ne.g. plink 1.9 --make-bed to sort your data.\n");
	  goto plink2_ret_INCONSISTENT_INPUT;
	}
	if (pcp->window_bp != -1) {
	  logerrprint("Error: --window requires a sorted .pvar/.bim.  Retry this command\nafter using e.g. plink 1.9 --make-bed to sort your data.\n");
	  goto plink2_ret_INCONSISTENT_INPUT;
	}
      }
      if (pcp->varid_from || pcp->varid_to) {
	reterr = from_to_flag(variant_ids, variant_id_htable, pcp->varid_from, pcp->varid_to, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, variant_include, cip, &variant_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->varid_snp) {
	reterr = snp_flag(variant_bps, variant_ids, variant_id_htable, pcp->varid_snp, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 0, pcp->window_bp, variant_include, cip, &variant_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->snps_range_list.name_ct) {
	reterr = snps_flag(variant_ids, variant_id_htable, &(pcp->snps_range_list), raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 0, variant_include, &variant_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->varid_exclude_snp) {
	reterr = snp_flag(variant_bps, variant_ids, variant_id_htable, pcp->varid_exclude_snp, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 1, pcp->window_bp, variant_include, cip, &variant_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->exclude_snps_range_list.name_ct) {
	reterr = snps_flag(variant_ids, variant_id_htable, &(pcp->exclude_snps_range_list), raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 1, variant_include, &variant_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }

      if (pcp->extract_fnames && (!(pcp->misc_flags & kfMiscExtractRange))) {
	reterr = extract_exclude_flag_norange(variant_ids, variant_id_htable, pcp->extract_fnames, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 0, variant_include, &variant_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->exclude_fnames && (!(pcp->misc_flags & kfMiscExcludeRange))) {
	reterr = extract_exclude_flag_norange(variant_ids, variant_id_htable, pcp->exclude_fnames, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 1, variant_include, &variant_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      bigstack_reset(bigstack_mark);
      if (variant_ct && full_variant_id_htable_needed) {
	reterr = apply_variant_bp_filters(pcp->extract_fnames, pcp->exclude_fnames, cip, variant_bps, pcp->from_bp, pcp->to_bp, raw_variant_ct, pcp->misc_flags, vpos_sortstatus, variant_include, &variant_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
    }
    // xid_mode may vary between these operations in a single run, and
    // sample-sort is relatively cheap, so we abandon plink 1.9's "construct
    // sample ID map only once" optimization.
    if (pcp->update_sex_fname) {
      reterr = update_sample_sexes(pcp->update_sex_fname, sample_include, sample_ids, raw_sample_ct, sample_ct, max_sample_id_blen, pcp->update_sex_colm2, sex_nm, sex_male);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
    if (pcp->keepfam_fnames) {
      reterr = keep_or_remove(pcp->keepfam_fnames, sample_ids, sids, raw_sample_ct, max_sample_id_blen, max_sid_blen, kfKeepFam, sample_include, &sample_ct);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
    if (pcp->keep_fnames) {
      reterr = keep_or_remove(pcp->keep_fnames, sample_ids, sids, raw_sample_ct, max_sample_id_blen, max_sid_blen, (keep_flags_t)(kfKeepForceSid * ((pcp->misc_flags / kfMiscKeepfileSid) & 1)), sample_include, &sample_ct);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
    if (pcp->removefam_fnames) {
      reterr = keep_or_remove(pcp->removefam_fnames, sample_ids, sids, raw_sample_ct, max_sample_id_blen, max_sid_blen, kfKeepRemove | kfKeepFam, sample_include, &sample_ct);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
    if (pcp->remove_fnames) {
      reterr = keep_or_remove(pcp->remove_fnames, sample_ids, sids, raw_sample_ct, max_sample_id_blen, max_sid_blen, kfKeepRemove | ((keep_flags_t)(kfKeepForceSid * ((pcp->misc_flags / kfMiscRemovefileSid) & 1))), sample_include, &sample_ct);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
    uint32_t* sample_missing_dosage_cts = nullptr;
    uint32_t* sample_missing_hc_cts = nullptr;
    uint32_t* sample_hethap_cts = nullptr;
    uintptr_t max_covar_name_blen = 0;
    if (psamname[0]) {
      if (pcp->misc_flags & kfMiscRequirePheno) {
        reterr = require_pheno(pheno_cols, pheno_names, require_pheno_flattened, raw_sample_ct, pheno_ct, max_pheno_name_blen, 0, sample_include, &sample_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->filter_flags & (kfFilterExclFemales | kfFilterExclMales | kfFilterExclNosex)) {
	if (pcp->filter_flags & kfFilterExclFemales) {
	  for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
	    sample_include[widx] &= (~sex_nm[widx]) | sex_male[widx];
	  }
	}
	if (pcp->filter_flags & kfFilterExclMales) {
	  bitvec_andnot(sex_male, raw_sample_ctl, sample_include);
	}
	if (pcp->filter_flags & kfFilterExclNosex) {
	  bitvec_and(sex_nm, raw_sample_ctl, sample_include);
	}
	const uint32_t old_sample_ct = sample_ct;
	sample_ct = popcount_longs(sample_include, raw_sample_ctl);
	const uint32_t removed_ct = old_sample_ct - sample_ct;
	LOGPRINTF("%u sample%s removed due to sex filter(s).\n", removed_ct, (removed_ct == 1)? "" : "s");
      }
      if (pcp->filter_flags & (kfFilterExclFounders | kfFilterExclNonfounders)) {
	const uint32_t keep_founders = (pcp->filter_flags / kfFilterExclNonfounders) & 1;
	if (keep_founders) {
	  bitvec_and(founder_info, raw_sample_ctl, sample_include);
	} else {
	  bitvec_andnot(founder_info, raw_sample_ctl, sample_include);
	}
	const uint32_t old_sample_ct = sample_ct;
	sample_ct = popcount_longs(sample_include, raw_sample_ctl);
	const uint32_t removed_ct = old_sample_ct - sample_ct;
	LOGPRINTF("--keep-%sfounders: %u sample%s removed.\n", keep_founders? "" : "non", removed_ct, (removed_ct == 1)? "" : "s");
      }
      const uint32_t smaj_missing_geno_report_requested = (pcp->command_flags1 & kfCommand1MissingReport) && (!(pcp->missing_rpt_modifier & kfMissingRptVariantOnly));
      if ((pcp->mind_thresh < 1.0) || smaj_missing_geno_report_requested) {
	if (bigstack_alloc_ui(raw_sample_ct, &sample_missing_hc_cts) ||
	    bigstack_alloc_ui(raw_sample_ct, &sample_hethap_cts)) {
	  goto plink2_ret_NOMEM;
	}
	if (are_sample_missing_dosage_cts_needed(pcp->misc_flags, smaj_missing_geno_report_requested, pcp->mind_thresh, pcp->missing_rpt_modifier)) {
	  if (pgfi.gflags & kfPgenGlobalDosagePresent) {
	    if (bigstack_alloc_ui(raw_sample_ct, &sample_missing_dosage_cts)) {
	      goto plink2_ret_NOMEM;
	    }
	  } else {
	    sample_missing_dosage_cts = sample_missing_hc_cts;
	  }
	}
	// could avoid this call and make load_allele_and_geno_counts() do
	// double duty with --missing?
	reterr = load_sample_missing_cts(sex_male, variant_include, cip, raw_variant_ct, variant_ct, raw_sample_ct, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, sample_missing_hc_cts, (pgfi.gflags & kfPgenGlobalDosagePresent)? sample_missing_dosage_cts : nullptr, sample_hethap_cts);
	if (reterr) {
	  goto plink2_ret_1;
	}
	if (pcp->mind_thresh < 1.0) {
	  uint32_t variant_ct_y = 0;
	  int32_t y_code;
	  if (xymt_exists(cip, kChrOffsetY, &y_code)) {
	    variant_ct_y = count_chr_variants_unsafe(variant_include, cip, y_code);
	  }
	  reterr = mind_filter((pcp->misc_flags & kfMiscMindDosage)? sample_missing_dosage_cts : sample_missing_hc_cts, (pcp->misc_flags & kfMiscMindHhMissing)? sample_hethap_cts : nullptr, sample_ids, sids, raw_sample_ct, max_sample_id_blen, max_sid_blen, variant_ct, variant_ct_y, pcp->mind_thresh, sample_include, sex_male, &sample_ct, outname, outname_end);
	  if (reterr) {
	    goto plink2_ret_1;
	  }
	}
	if (!smaj_missing_geno_report_requested) {
	  bigstack_reset(sample_missing_hc_cts);
	}
	// this results in a small "memory leak" when a regular missingness
	// report is requested, not a big deal
      }
      if (pcp->covar_fname || pcp->covar_range_list.name_ct) {
	const char* cur_covar_fname = pcp->covar_fname? pcp->covar_fname : (pcp->pheno_fname? pcp->pheno_fname : psamname);
	reterr = load_phenos(cur_covar_fname, &(pcp->covar_range_list), sample_include, sample_ids, raw_sample_ct, sample_ct, max_sample_id_blen, pcp->missing_pheno, 2, &covar_cols, &covar_names, &covar_ct, &max_covar_name_blen);
	if (reterr) {
	  goto plink2_ret_1;
	}
	LOGPRINTF("%u covariate%s loaded from %s.\n", covar_ct, (covar_ct == 1)? "" : "s", cur_covar_fname);

	// do we still want to clear some main phenotype values here if some
	// covariate values are missing?  (don't think there's a point to
	// preserving that behavior, let the regression functions do it to
	// their local phenotype copies on their own.)
      }

      if (pcp->misc_flags & kfMiscRequireCovar) {
        reterr = require_pheno(covar_cols, covar_names, require_covar_flattened, raw_sample_ct, covar_ct, max_covar_name_blen, 1, sample_include, &sample_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->keep_if_expr.pheno_name) {
	reterr = keep_remove_if(&(pcp->keep_if_expr), pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, (pcp->misc_flags / kfMiscAffection01) & 1, 0, sample_include, &sample_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->remove_if_expr.pheno_name) {
	reterr = keep_remove_if(&(pcp->remove_if_expr), pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, (pcp->misc_flags / kfMiscAffection01) & 1, 1, sample_include, &sample_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      // meow
      if (pcp->keep_cats_fname || pcp->keep_cat_names_flattened) {
	reterr = keep_remove_cats(pcp->keep_cats_fname, pcp->keep_cat_names_flattened, pcp->keep_cat_phenoname, pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, 0, pcp->max_thread_ct, sample_include, &sample_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->remove_cats_fname || pcp->remove_cat_names_flattened) {
	reterr = keep_remove_cats(pcp->remove_cats_fname, pcp->remove_cat_names_flattened, pcp->remove_cat_phenoname, pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, 1, pcp->max_thread_ct, sample_include, &sample_ct);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
    }

    const uint32_t nonfounders = (pcp->misc_flags / kfMiscNonfounders) & 1;
    uint32_t founder_ct = 0;
    uint32_t male_ct = 0;
    uint32_t nosex_ct = 0;
    if (psamname[0]) {
      if ((!sample_ct) && (!(pcp->misc_flags & kfMiscAllowNoSamples))) {
	logerrprint("Error: No samples remaining after main filters.  (Add --allow-no-samples to\npermit this.)\n");
	goto plink2_ret_INCONSISTENT_INPUT;
      }
      update_sample_subsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
      if (pcp->filter_flags) {
	const uint32_t female_ct = sample_ct - male_ct - nosex_ct;
	if (!nosex_ct) {
	  LOGPRINTFWW("%u sample%s (%u female%s, %u male%s; %u founder%s) remaining after main filters.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", founder_ct, (founder_ct == 1)? "" : "s");
	} else {
	  LOGPRINTFWW("%u sample%s (%u female%s, %u male%s, %u ambiguous; %u founder%s) remaining after main filters.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", nosex_ct, founder_ct, (founder_ct == 1)? "" : "s");
	}
	if ((pheno_ct == 1) && (pheno_cols[0].type_code == kPhenoDtypeCc)) {
	  const uint32_t obs_ct = popcount_longs_intersect(pheno_cols[0].nonmiss, sample_include, raw_sample_ctl);
	  const uint32_t case_ct = popcount_longs_intersect(pheno_cols[0].data.cc, sample_include, raw_sample_ctl);
	  const uint32_t ctrl_ct = obs_ct - case_ct;
	  LOGPRINTF("%u case%s and %u control%s remaining after main filters.\n", case_ct, (case_ct == 1)? "" : "s", ctrl_ct, (ctrl_ct == 1)? "" : "s");
	}
      }
    }
    if (pcp->pheno_transform_flags & kfPhenoTransformSplitCat) {
      reterr = split_cat_pheno(pcp->split_cat_phenonames_flattened, sample_include, raw_sample_ct, pcp->pheno_transform_flags, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen, &covar_cols, &covar_names, &covar_ct, &max_covar_name_blen);
      if (reterr) {
	goto plink2_ret_1;
      }
    }

    // quantnorm before variance-standardize, since at least that has a minor
    // effect, whereas the other order is pointless
    if (pcp->pheno_transform_flags & (kfPhenoTransformQuantnormPheno | kfPhenoTransformQuantnormAll)) {
      reterr = pheno_quantile_normalize(pcp->quantnorm_flattened, sample_include, pheno_names, raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, 0, (pcp->pheno_transform_flags / kfPhenoTransformQuantnormPheno) & 1, pheno_cols);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
    if (pcp->pheno_transform_flags & (kfPhenoTransformQuantnormCovar | kfPhenoTransformQuantnormAll)) {
      reterr = pheno_quantile_normalize((pcp->pheno_transform_flags & kfPhenoTransformQuantnormAll)? pcp->quantnorm_flattened : pcp->covar_quantnorm_flattened, sample_include, covar_names, raw_sample_ct, sample_ct, covar_ct, max_covar_name_blen, 1, (pcp->pheno_transform_flags / kfPhenoTransformQuantnormCovar) & 1, covar_cols);
      if (reterr) {
	goto plink2_ret_1;
      }
    }

    if (pcp->pheno_transform_flags & (kfPhenoTransformVstdCovar | kfPhenoTransformVstdAll)) {
      const uint32_t is_covar_flag = (pcp->pheno_transform_flags / kfPhenoTransformVstdCovar) & 1;
      if (!is_covar_flag) {
	reterr = pheno_variance_standardize(pcp->vstd_flattened, sample_include, pheno_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, 0, 0, pheno_cols);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      reterr = pheno_variance_standardize(pcp->vstd_flattened, sample_include, covar_names, raw_sample_ct, covar_ct, max_covar_name_blen, 1, is_covar_flag, covar_cols);
      if (reterr) {
	goto plink2_ret_1;
      }
    }

    // dosages are currently in 32768ths
    uint64_t* allele_dosages = nullptr; // same indexes as allele_storage
    uint64_t* founder_allele_dosages = nullptr;
    alt_allele_ct_t* maj_alleles = nullptr;
    double* allele_freqs = nullptr;
    uint32_t* raw_geno_cts = nullptr;
    uint32_t* founder_raw_geno_cts = nullptr;
    unsigned char* bigstack_mark_allele_dosages = g_bigstack_base;
    unsigned char* bigstack_mark_founder_allele_dosages = g_bigstack_base;
    if (pgenname[0]) {
      if (are_allele_freqs_needed(pcp->command_flags1, pcp->min_maf, pcp->max_maf)) {
	if (are_maj_alleles_needed(pcp->command_flags1)) {
	  maj_alleles = (alt_allele_ct_t*)bigstack_alloc(raw_variant_ct * sizeof(alt_allele_ct_t));
	  if (!maj_alleles) {
	    goto plink2_ret_NOMEM;
	  }
	}
	//   allele_freqs[variant_allele_idxs[variant_uidx] - variant_uidx]
	// stores the frequency estimate for the reference allele; if there's
	// more than 1 alt allele, next element stores alt1 freq, etc.  To save
	// memory, we omit the last alt.
	uintptr_t total_alt_allele_ct = raw_variant_ct;
	if (variant_allele_idxs) {
	  total_alt_allele_ct = variant_allele_idxs[raw_variant_ct] - raw_variant_ct;
	}
	if (bigstack_alloc_d(total_alt_allele_ct, &allele_freqs)) {
	  goto plink2_ret_NOMEM;
	}
      }
      uint32_t x_start = 0;
      uint32_t x_len = 0;
      uint32_t hwe_x_probs_needed = 0;
      int32_t x_code;
      if ((!(vpos_sortstatus & kfUnsortedVarSplitChr)) && xymt_exists(cip, kChrOffsetX, &x_code)) {
	const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)x_code];
	x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
	const uint32_t x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
	x_len = x_end - x_start;
	if (x_len && ((pcp->command_flags1 & kfCommand1Hardy) || (pcp->hwe_thresh != 1.0)) && (!are_all_bits_zero(variant_include, x_start, x_end))) {
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
      const uint32_t first_hap_uidx = get_first_haploid_uidx(cip, vpos_sortstatus);
      if (are_allele_dosages_needed(pcp->misc_flags, make_plink2_modifier, (allele_freqs != nullptr), pcp->min_allele_dosage, pcp->max_allele_dosage)) {
	if (bigstack_alloc_ull(variant_allele_idxs? variant_allele_idxs[raw_variant_ct] : (2 * raw_variant_ct), &allele_dosages)) {
	  goto plink2_ret_NOMEM;
	}
      }
      bigstack_mark_founder_allele_dosages = g_bigstack_base;
      if (are_founder_allele_dosages_needed(pcp->misc_flags, (allele_freqs != nullptr), pcp->min_allele_dosage, pcp->max_allele_dosage)) {
	if ((founder_ct == sample_ct) && allele_dosages) {
	  founder_allele_dosages = allele_dosages;
	} else {
	  if (bigstack_alloc_ull(variant_allele_idxs? variant_allele_idxs[raw_variant_ct] : (2 * raw_variant_ct), &founder_allele_dosages)) {
	    goto plink2_ret_NOMEM;
	  }
	}
      }
      double* mach_r2_vals = nullptr;
      if ((pcp->allele_freq_modifier & kfAlleleFreqColMachR2) || (pcp->mach_r2_max != 0.0)) {
	if (bigstack_alloc_d(raw_variant_ct, &mach_r2_vals)) {
	  goto plink2_ret_NOMEM;
	}
      }
      
      unsigned char* bigstack_mark_geno_cts = g_bigstack_base;
      
      // no longer includes hethaps by default
      uint32_t* variant_missing_hc_cts = nullptr;
      uint32_t* variant_hethap_cts = nullptr;
      if (are_variant_missing_hc_cts_needed(pcp->command_flags1, pcp->misc_flags, pcp->geno_thresh, pcp->missing_rpt_modifier)) {
	if (bigstack_alloc_ui(raw_variant_ct, &variant_missing_hc_cts)) {
	  goto plink2_ret_NOMEM;
	}
	if (are_variant_hethap_cts_needed(pcp->command_flags1, pcp->misc_flags, pcp->geno_thresh, pcp->missing_rpt_modifier, first_hap_uidx)) {
	  // first_hap_uidx offset can save an entire GB...
	  if (bigstack_alloc_ui(raw_variant_ct - first_hap_uidx, &variant_hethap_cts)) {
	    goto plink2_ret_NOMEM;
	  }
	}
      }
      uint32_t* variant_missing_dosage_cts = nullptr;
      if (are_variant_missing_dosage_cts_needed(pcp->command_flags1, pcp->misc_flags, pcp->geno_thresh, pcp->missing_rpt_modifier)) {
	if ((!variant_missing_hc_cts) || (pgfi.gflags & kfPgenGlobalDosagePresent)) {
	  if (bigstack_alloc_ui(raw_variant_ct, &variant_missing_dosage_cts)) {
	    goto plink2_ret_NOMEM;
	  }
	} else {
	  variant_missing_dosage_cts = variant_missing_hc_cts;
	}
      }
      uint32_t* x_male_geno_cts = nullptr;
      uint32_t* founder_x_male_geno_cts = nullptr;
      uint32_t* x_nosex_geno_cts = nullptr;
      uint32_t* founder_x_nosex_geno_cts = nullptr;
      // [3n] = homref ct, [3n+1] = het ref-altx total, [3n+2] = nonref diploid
      //   total
      // use unfiltered indexes, since we remove more variants later
      if (are_raw_geno_cts_needed(pcp->command_flags1, pcp->misc_flags, pcp->hwe_thresh)) {
	if (bigstack_alloc_ui((3 * k1LU) * raw_variant_ct, &raw_geno_cts)) {
	  goto plink2_ret_NOMEM;
	}
	if (x_len) {
	  if (male_ct) {
	    if (bigstack_alloc_ui((3 * k1LU) * x_len, &x_male_geno_cts)) {
	      goto plink2_ret_NOMEM;
	    }
	  }
	  if (nosex_ct && hwe_x_probs_needed && nonfounders) {
	    if (bigstack_alloc_ui((3 * k1LU) * x_len, &x_nosex_geno_cts)) {
	      goto plink2_ret_NOMEM;
	    }
	  }
	}
      }
      if (are_founder_raw_geno_cts_needed(pcp->command_flags1, pcp->misc_flags, pcp->hwe_thresh)) {
	if ((founder_ct == sample_ct) && raw_geno_cts) {
	  founder_raw_geno_cts = raw_geno_cts;
	  founder_x_male_geno_cts = x_male_geno_cts;
	} else {
	  if (bigstack_alloc_ui((3 * k1LU) * raw_variant_ct, &founder_raw_geno_cts)) {
	    goto plink2_ret_NOMEM;
	  }
	  if (x_len && male_ct) {
	    const uint32_t founder_male_ct = popcount_longs_intersect(founder_info, sex_male, raw_sample_ctl);
	    if (founder_male_ct) {
	      if (bigstack_alloc_ui((3 * k1LU) * x_len, &founder_x_male_geno_cts)) {
		goto plink2_ret_NOMEM;
	      }
	    }
	  }
	}
	if (nosex_ct && hwe_x_probs_needed && (!nonfounders)) {
	  const uint32_t founder_knownsex_ct = popcount_longs_intersect(founder_info, sex_nm, raw_sample_ctl);
	  if (founder_knownsex_ct < founder_ct) {
	    if ((founder_ct == sample_ct) && x_nosex_geno_cts) {
	      // shouldn't be possible for now
	      assert(0);
	      // founder_x_nosex_geno_cts = x_nosex_geno_cts;
	    } else {
	      if (bigstack_alloc_ui((3 * k1LU) * x_len, &founder_x_nosex_geno_cts)) {
		goto plink2_ret_NOMEM;
	      }
	    }
	  }
	}
      }
      if (allele_dosages || founder_allele_dosages || variant_missing_hc_cts || variant_missing_dosage_cts || variant_hethap_cts || raw_geno_cts || founder_raw_geno_cts || mach_r2_vals) {
	// note that --geno depends on different handling of X/Y than --maf.

	// possible todo: "free" these arrays early in some cases
	// todo: oblig-missing
	reterr = load_allele_and_geno_counts(sample_include, founder_info, sex_nm, sex_male, variant_include, cip, variant_allele_idxs, raw_sample_ct, sample_ct, founder_ct, male_ct, nosex_ct, raw_variant_ct, variant_ct, first_hap_uidx, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, allele_dosages, founder_allele_dosages, variant_missing_hc_cts, (pgfi.gflags & kfPgenGlobalDosagePresent)? variant_missing_dosage_cts : nullptr, variant_hethap_cts, raw_geno_cts, founder_raw_geno_cts, x_male_geno_cts, founder_x_male_geno_cts, x_nosex_geno_cts, founder_x_nosex_geno_cts, mach_r2_vals);
	if (reterr) {
	  goto plink2_ret_1;
	}
	if (pcp->command_flags1 & kfCommand1GenotypingRate) {
	  const uint32_t is_dosage = (pcp->misc_flags / kfMiscGenotypingRateDosage) & 1;
	  report_genotyping_rate(variant_include, cip, is_dosage? variant_missing_dosage_cts : variant_missing_hc_cts, raw_sample_ct, sample_ct, male_ct, variant_ct, is_dosage);
	  if (!(pcp->command_flags1 & (~kfCommand1GenotypingRate))) {
	    goto plink2_ret_1;
	  }
	}
      }
      if (allele_freqs) {
	const uint32_t maf_succ = (pcp->misc_flags / kfMiscMafSucc) & 1;
	compute_allele_freqs(variant_include, variant_allele_idxs, nonfounders? allele_dosages : founder_allele_dosages, variant_ct, maf_succ, allele_freqs);
	if (pcp->read_freq_fname) {
	  reterr = read_allele_freqs(variant_include, variant_ids, variant_allele_idxs, allele_storage, pcp->read_freq_fname, raw_variant_ct, variant_ct, pgfi.max_alt_allele_ct, max_variant_id_slen, max_allele_slen, maf_succ, pcp->max_thread_ct, allele_freqs);
	  if (reterr) {
	    goto plink2_ret_1;
	  }
	}
	if (maj_alleles) {
	  compute_maj_alleles(variant_include, variant_allele_idxs, allele_freqs, variant_ct, maj_alleles);
	}
      } else if (pcp->read_freq_fname) {
	LOGERRPRINTF("Warning: Ignoring --read-freq since no command would use the frequencies.\n");
      }

      if (pcp->command_flags1 & kfCommand1AlleleFreq) {
	reterr = write_allele_freqs(variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, nonfounders? allele_dosages : founder_allele_dosages, mach_r2_vals, pcp->freq_ref_binstr, pcp->freq_alt1_binstr, variant_ct, pgfi.max_alt_allele_ct, max_allele_slen, pcp->allele_freq_modifier, nonfounders, outname, outname_end);
	if (reterr || (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1AlleleFreq))))) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->command_flags1 & kfCommand1GenoCounts) {
	reterr = write_geno_counts(sample_include, sex_male, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, raw_geno_cts, x_male_geno_cts, raw_sample_ct, sample_ct, male_ct, variant_ct, x_start, max_allele_slen, pcp->geno_counts_modifier, &simple_pgr, outname, outname_end);
	if (reterr || (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1AlleleFreq | kfCommand1GenoCounts))))) {
	  goto plink2_ret_1;
	}
      }
      
      if (pcp->command_flags1 & kfCommand1MissingReport) {
	reterr = write_missingness_reports(sample_include, sex_male, sample_ids, sids, pheno_cols, pheno_names, sample_missing_hc_cts, sample_missing_dosage_cts, sample_hethap_cts, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, variant_missing_hc_cts, variant_missing_dosage_cts, variant_hethap_cts, sample_ct, male_ct, max_sample_id_blen, max_sid_blen, pheno_ct, max_pheno_name_blen, variant_ct, max_allele_slen, variant_hethap_cts? first_hap_uidx : 0x7fffffff, pcp->missing_rpt_modifier, outname, outname_end);
	if (reterr || (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1AlleleFreq | kfCommand1GenoCounts | kfCommand1MissingReport))))) {
	  goto plink2_ret_1;
	}
      }

      if (pcp->geno_thresh != 1.0) {
	const uint32_t geno_hh_missing = (uint32_t)(pcp->misc_flags & kfMiscGenoHhMissing);
	enforce_geno_thresh(cip, (pcp->misc_flags & kfMiscGenoDosage)? variant_missing_dosage_cts : variant_missing_hc_cts, geno_hh_missing? variant_hethap_cts : nullptr, sample_ct, male_ct, geno_hh_missing? first_hap_uidx : 0x7fffffff, pcp->geno_thresh, variant_include, &variant_ct);
      }

      double* hwe_x_pvals = nullptr;
      uint32_t hwe_x_ct = 0;
      if (hwe_x_probs_needed) {
	hwe_x_ct = count_chr_variants_unsafe(variant_include, cip, cip->xymt_codes[kChrOffsetX]);
	// hwe_x_ct == 0 possible, if --geno filters out all remaining chrX
	// variants
	// also support suppression of --hardy p column (with a gigantic
	// dataset, maybe it's reasonable to stick to femalep, etc.)
	if (hwe_x_ct && ((pcp->hwe_thresh != 1.0) || (pcp->hardy_modifier & kfHardyColP))) {
	  uint32_t hwe_midp;
	  if (pcp->command_flags1 & kfCommand1Hardy) {
	    hwe_midp = (pcp->hardy_modifier / kfHardyMidp) & 1;
	    if (pcp->hwe_thresh != 1.0) {
	      const uint32_t hwe_midp2 = (pcp->misc_flags / kfMiscHweMidp) & 1;
	      if (hwe_midp != hwe_midp2) {
		// could support this efficiently, but why bother...
		logerrprint("Error: --hardy and --hwe must have identical midp settings when chrX is\npresent.\n");
		goto plink2_ret_INVALID_CMDLINE;
	      }
	    }
	  } else {
	    hwe_midp = (pcp->misc_flags / kfMiscHweMidp) & 1;
	  }
	  reterr = compute_hwe_x_pvals(variant_include, nonfounders? raw_geno_cts : founder_raw_geno_cts, nonfounders? x_male_geno_cts : founder_x_male_geno_cts, nonfounders? x_nosex_geno_cts : founder_x_nosex_geno_cts, x_start, hwe_x_ct, hwe_midp, pcp->max_thread_ct, &hwe_x_pvals);
	  if (reterr) {
	    goto plink2_ret_1;
	  }
	}
      }
      if (pcp->command_flags1 & kfCommand1Hardy) {
	reterr = hardy_report(variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, nonfounders? raw_geno_cts : founder_raw_geno_cts, nonfounders? x_male_geno_cts : founder_x_male_geno_cts, nonfounders? x_nosex_geno_cts : founder_x_nosex_geno_cts, hwe_x_pvals, variant_ct, hwe_x_ct, max_allele_slen, pcp->output_min_p, pcp->hardy_modifier, nonfounders, outname, outname_end);
	if (reterr || (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1AlleleFreq | kfCommand1GenoCounts | kfCommand1MissingReport | kfCommand1Hardy))))) {
	  goto plink2_ret_1;
	}
      }
      if (pcp->hwe_thresh != 1.0) {
	// assumes no filtering between hwe_x_pvals[] computation and here
	enforce_hwe_thresh(cip, nonfounders? raw_geno_cts : founder_raw_geno_cts, nonfounders? x_male_geno_cts : founder_x_male_geno_cts, nonfounders? x_nosex_geno_cts : founder_x_nosex_geno_cts, hwe_x_pvals, pcp->misc_flags, pcp->hwe_thresh, nonfounders, variant_include, &variant_ct);
      }
      // raw_geno_cts/founder_raw_geno_cts/hwe_x_pvals no longer needed
      bigstack_reset(bigstack_mark_geno_cts);

      if ((pcp->min_maf != 0.0) || (pcp->max_maf != 1.0) || pcp->min_allele_dosage || (pcp->max_allele_dosage != (~0LLU))) {
	enforce_minor_freq_constraints(variant_allele_idxs, nonfounders? allele_dosages : founder_allele_dosages, allele_freqs, pcp->min_maf, pcp->max_maf, pcp->min_allele_dosage, pcp->max_allele_dosage, variant_include, &variant_ct);
      }

      if (mach_r2_vals) {
	if (pcp->mach_r2_max != 0.0) {
	  enforce_mach_r2_thresh(cip, mach_r2_vals, pcp->mach_r2_min, pcp->mach_r2_max, variant_include, &variant_ct);
	}
        bigstack_reset(mach_r2_vals);
      }

      if ((!variant_ct) && (!(pcp->misc_flags & kfMiscAllowNoVars))) {
	logerrprint("Error: No variants remaining after main filters.  (Add --allow-no-vars to\npermit this.)\n");
	goto plink2_ret_INCONSISTENT_INPUT;
      }
      if (pcp->filter_flags) {
	LOGPRINTF("%u variant%s remaining after main filters.\n", variant_ct, (variant_ct == 1)? "" : "s");
      }

      if (pcp->command_flags1 & (kfCommand1MakeKing | kfCommand1KingCutoff)) {
	uintptr_t* prev_sample_include = nullptr;
	const uint32_t prev_sample_ct = sample_ct;
	if (pcp->king_cutoff != -1) {
	  if (bigstack_alloc_ul(raw_sample_ctl, &prev_sample_include)) {
	    goto plink2_ret_NOMEM;
	  }
	  memcpy(prev_sample_include, sample_include, raw_sample_ctl * sizeof(intptr_t));
	}
	if (king_cutoff_fprefix) {
	  reterr = king_cutoff_batch(sample_ids, sids, raw_sample_ct, max_sample_id_blen, max_sid_blen, pcp->king_cutoff, sample_include, king_cutoff_fprefix, &sample_ct);
	} else {
	  reterr = calc_king(sample_ids, sids, variant_include, cip, raw_sample_ct, max_sample_id_blen, max_sid_blen, raw_variant_ct, variant_ct, pcp->king_cutoff, pcp->king_table_filter, pcp->king_modifier, pcp->parallel_idx, pcp->parallel_tot, pcp->max_thread_ct, &simple_pgr, sample_include, &sample_ct, outname, outname_end);
	}
	if (reterr) {
	  goto plink2_ret_1;
	}
	if (pcp->king_cutoff != -1) {
	  strcpy(outname_end, ".king.cutoff.in");
	  reterr = write_sample_ids(sample_include, sample_ids, sids, outname, sample_ct, max_sample_id_blen, max_sid_blen);
	  if (reterr) {
	    goto plink2_ret_1;
	  }
	  strcpy(&(outname_end[13]), "out");
	  bitvec_andnot(sample_include, raw_sample_ctl, prev_sample_include);
	  const uint32_t removed_sample_ct = prev_sample_ct - sample_ct;
	  reterr = write_sample_ids(prev_sample_include, sample_ids, sids, outname, removed_sample_ct, max_sample_id_blen, max_sid_blen);
	  if (reterr) {
	    goto plink2_ret_1;
	  }
	  bigstack_reset(prev_sample_include);
	  outname_end[13] = '\0';
	  LOGPRINTFWW("--king-cutoff: Excluded sample ID%s written to %sout, and %u remaining sample ID%s written to %sin .\n", (removed_sample_ct == 1)? "" : "s", outname, sample_ct, (sample_ct == 1)? "" : "s", outname);
	  update_sample_subsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
	}
      }
    }
    double* grm = nullptr;
    const uint32_t keep_grm = grm_keep_needed(pcp->command_flags1, pcp->pca_flags);
    if ((pcp->command_flags1 & kfCommand1MakeRel) || keep_grm) {
      reterr = calc_grm(sample_include, sample_ids, sids, variant_include, cip, variant_allele_idxs, maj_alleles, allele_freqs, raw_sample_ct, sample_ct, max_sample_id_blen, max_sid_blen, raw_variant_ct, variant_ct, pcp->grm_flags, pcp->parallel_idx, pcp->parallel_tot, pcp->max_thread_ct, &simple_pgr, outname, outname_end, keep_grm? (&grm) : nullptr);
      if (reterr) {
	goto plink2_ret_1;
      }
      // don't bother with --rel-cutoff for now, since --king-cutoff seems to
      // work better...
      
      // possible todo: unrelated heritability?
    }
#ifndef NOLAPACK
    if (pcp->command_flags1 & kfCommand1Pca) {
      // if the GRM is on the stack, this always frees it
      reterr = calc_pca(sample_include, sample_ids, sids, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, maj_alleles, allele_freqs, raw_sample_ct, sample_ct, max_sample_id_blen, max_sid_blen, raw_variant_ct, variant_ct, max_allele_slen, pcp->pca_ct, pcp->pca_flags, pcp->max_thread_ct, &simple_pgr, grm, outname, outname_end);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
#endif
    
    if (pcp->command_flags1 & kfCommand1WriteSnplist) {
      reterr = write_snplist(variant_include, variant_ids, variant_ct, (pcp->misc_flags / kfMiscWriteSnplistZs) & 1, outname, outname_end);
      if (reterr) {
	goto plink2_ret_1;
      }
    }

    if (pcp->command_flags1 & (kfCommand1MakePlink2 | kfCommand1Exportf | kfCommand1WriteCovar)) {
      // If non-null, this has (2 * raw_variant_ct) entries.  [2n] stores the
      // index of the new ref allele for variant n, and [2n+1] stores the index
      // of the new alt1 allele.  (0 = original ref, 1 = original alt1, etc.)
      // If at all possible, operations which instantiate this
      // (--ref-allele, --alt1-allele, ...) should only be usable with fileset
      // creation commands.  no more pass-marker_reverse-to-everything
      // nonsense.
      // (Technically, I could also drop support for --export, but that would
      // force too many real-world jobs to require two plink2 runs instead of
      // one.)
      unsigned char* bigstack_end_mark = g_bigstack_end;
      alt_allele_ct_t* refalt1_select = nullptr;
      if (pcp->misc_flags & kfMiscMajRef) {
	// todo: also support updated version of --a2-allele, etc.
	const uintptr_t refalt1_word_ct = DIV_UP(2 * raw_variant_ct * sizeof(alt_allele_ct_t), kBytesPerWord);
	uintptr_t* refalt1_select_ul;
	if (bigstack_end_alloc_ul(refalt1_word_ct, &refalt1_select_ul)) {
	  goto plink2_ret_NOMEM;
	}
	const uintptr_t alt_allele_vals = (uintptr_t)(k1LU << (8 * sizeof(alt_allele_ct_t)));
	const uintptr_t fill_word = ((~k0LU) / ((alt_allele_vals - 1) * (alt_allele_vals + 1))) * alt_allele_vals;
	for (uintptr_t widx = 0; widx < refalt1_word_ct; ++widx) {
	  refalt1_select_ul[widx] = fill_word;
	}
	refalt1_select = (alt_allele_ct_t*)refalt1_select_ul;
	if (pcp->misc_flags & kfMiscMajRef) {
	  // possible todo: make this subscribe to maj_alleles[] instead?
	  // might be pointless due to ALT1 computation, though.

	  // todo: warning if this is specified without file write command, if
	  // this is ever moved out of the file-write subblock
	  const uint64_t* main_allele_dosages = nonfounders? allele_dosages : founder_allele_dosages;
	  const uint32_t not_all_nonref = !(pgfi.gflags & kfPgenGlobalAllNonref);
	  const uint32_t skip_real_ref = not_all_nonref && (!(pcp->misc_flags & kfMiscMajRefForce));
	  if (skip_real_ref && (!nonref_flags)) {
	    logerrprint("Warning: --maj-ref has no effect, since no provisional reference alleles are\npresent.  (Did you forget to add the 'force' modifier?)\n");
	  } else {
	    if (not_all_nonref && (!nonref_flags)) {
	      if (bigstack_end_alloc_ul(raw_variant_ctl, &nonref_flags)) {
		goto plink2_ret_NOMEM;
	      }
	      pgfi.nonref_flags = nonref_flags;
	      fill_ulong_zero(raw_variant_ctl, nonref_flags);
	    }
	    uint32_t variant_uidx = 0;
	    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	      next_set_unsafe_ck(variant_include, &variant_uidx);
	      if (skip_real_ref && IS_SET(nonref_flags, variant_uidx)) {
		continue;
	      }
	      const uint64_t* cur_allele_dosages = &(main_allele_dosages[variant_allele_idxs? variant_allele_idxs[variant_uidx] : (2 * variant_uidx)]);
	      const uint32_t alt_ct_p1 = variant_allele_idxs? (variant_allele_idxs[variant_uidx + 1] - variant_allele_idxs[variant_uidx]) : 2;
	      if (alt_ct_p1 == 2) {
		// optimize common case
		if (cur_allele_dosages[1] > cur_allele_dosages[0]) {
		  // assumes alt_allele_ct_t is unsigned char
		  ((uint16_t*)refalt1_select)[variant_uidx] = 1;
		  if (nonref_flags) {
		    SET_BIT(variant_uidx, nonref_flags);
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
		    SET_BIT(variant_uidx, nonref_flags);
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
      if (make_plink2_modifier & kfMakePlink2TrimAlts) {
	bigstack_reset(bigstack_mark_founder_allele_dosages);
      } else {
        bigstack_reset(bigstack_mark_allele_dosages);
      }

      uint32_t* new_sample_idx_to_old = nullptr;
      if (pcp->sample_sort_flags & (kfSortNatural | kfSortAscii | kfSortFile)) {
	if (sample_ct < 2) {
	  logerrprint("Warning: Skipping --sample-sort since <2 samples are present.\n");
	} else {
	  if (pcp->sample_sort_flags & kfSortFile) {
	    reterr = sample_sort_file_map(sample_include, sample_ids, sids, pcp->sample_sort_fname, raw_sample_ct, sample_ct, max_sample_id_blen, max_sid_blen, pcp->sample_sort_flags & kfSortFileSid, &new_sample_idx_to_old);
	    if (reterr) {
	      goto plink2_ret_1;
	    }
	  } else {
	    // probably more efficient to have --make-{bed,pgen,bpgen} perform
	    // an unfiltered load?  but we should have compute power to spare
	    // here, so keep the code simpler for now
	    char* sorted_xidbox_tmp;
	    uintptr_t max_xid_blen;
	    reterr = sorted_xidbox_init_alloc(sample_include, sample_ids, sids, sample_ct, max_sample_id_blen, max_sid_blen, sids? kfXidModeFidiidSid : kfXidModeFidiid, (pcp->sample_sort_flags == kfSortNatural), &sorted_xidbox_tmp, &new_sample_idx_to_old, &max_xid_blen);
	    if (reterr) {
	      goto plink2_ret_1;
	    }
	    bigstack_reset(sorted_xidbox_tmp);
	  }
	  LOGPRINTF("--indiv-sort: %u samples reordered.\n", sample_ct);
	}
      }

      if (covar_ct && ((pcp->command_flags1 & (kfCommand1Exportf | kfCommand1WriteCovar)) || ((pcp->command_flags1 & kfCommand1MakePlink2) && (make_plink2_modifier & (kfMakeBed | kfMakeFam | kfMakePgen | kfMakePsam))))) {
	reterr = write_covar(sample_include, sample_ids, sids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, pheno_names, covar_cols, covar_names, new_sample_idx_to_old, sample_ct, max_sample_id_blen, max_sid_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, pcp->write_covar_flags, outname, outname_end);
	if (reterr) {
	  goto plink2_ret_1;
	}
      } else if (pcp->command_flags1 & kfCommand1WriteCovar) {
	logerrprint("Warning: Skipping --write-covar, since no covariates are loaded.\n");
      }
      
      if (pcp->command_flags1 & kfCommand1MakePlink2) {
	// todo: unsorted case (--update-chr, etc.)
	if (vpos_sortstatus & kfUnsortedVarSplitChr) {
	  logerrprint("Error: --make-bed/--make-{b}pgen variant sorting is under development.\n");
	  reterr = kPglRetNotYetSupported;
	  goto plink2_ret_1;
	}
	if (vpos_sortstatus & kfUnsortedVarBp) {
	  logerrprint("Warning: --make-bed/--make-{b}pgen variant sorting is not implemented yet.\n");
	}
	reterr = make_plink2_no_vsort(xheader, sample_include, sample_ids, sids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, allele_dosages, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, info_reload_slen? pvarname : nullptr, variant_cms, xheader_blen, xheader_info_pr, raw_sample_ct, sample_ct, max_sample_id_blen, max_sid_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, pcp->max_thread_ct, pcp->hard_call_thresh, pcp->dosage_erase_thresh, make_plink2_modifier, pcp->pvar_psam_modifier, pgr_alloc_cacheline_ct, &pgfi, &simple_pgr, outname, outname_end);
	if (reterr) {
	  goto plink2_ret_1;
	}
        bigstack_reset(bigstack_mark_allele_dosages);
      }

      if (pcp->command_flags1 & kfCommand1Exportf) {
	reterr = exportf(xheader, sample_include, sample_ids, sids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, pheno_names, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, info_reload_slen? pvarname : nullptr, variant_cms, xheader_blen, xheader_info_pr, raw_sample_ct, sample_ct, max_sample_id_blen, max_sid_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, pcp->max_thread_ct, make_plink2_modifier, pcp->exportf_modifier, pcp->exportf_id_paste, pcp->exportf_id_delim, pcp->exportf_bits, pgr_alloc_cacheline_ct, &pgfi, &simple_pgr, outname, outname_end);
	if (reterr) {
	  goto plink2_ret_1;
	}
      }
      bigstack_end_reset(bigstack_end_mark);
    }
    bigstack_reset(bigstack_mark_allele_dosages);
    
    if (pcp->command_flags1 & kfCommand1LdPrune) {
      if ((pcp->ld_info.prune_modifier & kfLdPruneWindowBp) && (vpos_sortstatus & kfUnsortedVarBp)) {
	logerrprint("Error: When the window size is in kb units, LD-based pruning requires a sorted\n.pvar/.bim.  Retry this command after using e.g. plink 1.9 --make-bed to sort\nyour data.\n");
	goto plink2_ret_INCONSISTENT_INPUT;
      }
      reterr = ld_prune(variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, maj_alleles, allele_freqs, founder_info, sex_male, &(pcp->ld_info), raw_variant_ct, variant_ct, raw_sample_ct, founder_ct, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
      if (reterr) {
	goto plink2_ret_1;
      }
    }

    if (pcp->command_flags1 & kfCommand1Score) {
      reterr = score_report(sample_include, sample_ids, sids, sex_male, pheno_cols, pheno_names, variant_include, cip, variant_ids, variant_allele_idxs, allele_storage, allele_freqs, &(pcp->score_info), sample_ct, max_sample_id_blen, max_sid_blen, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_variant_id_slen, pcp->xchr_model, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
    // eventually check for nonzero pheno_ct here?
    
    if (pcp->command_flags1 & kfCommand1Glm) {
      reterr = glm_main(sample_include, sample_ids, sids, sex_nm, sex_male, pheno_cols, pheno_names, covar_cols, covar_names, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, &(pcp->glm_info), &(pcp->adjust_info), &(pcp->aperm), pcp->glm_local_covar_fname, pcp->glm_local_pvar_fname, pcp->glm_local_psam_fname, raw_sample_ct, sample_ct, max_sample_id_blen, max_sid_blen, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, raw_variant_ct, variant_ct, max_variant_id_slen, max_allele_slen, pcp->xchr_model, pcp->ci_size, pcp->vif_thresh, pcp->pfilter, pcp->output_min_p, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, &simple_pgr, outname, outname_end);
      if (reterr) {
	goto plink2_ret_1;
      }
    }
  }
  while (0) {
  plink2_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  plink2_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  plink2_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  plink2_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 plink2_ret_1:
  cleanup_pheno_cols(covar_ct, covar_cols);
  cleanup_pheno_cols(pheno_ct, pheno_cols);
  free_cond(covar_names);
  free_cond(pheno_names);
  if (pgr_cleanup(&simple_pgr) && (!reterr)) {
    reterr = kPglRetReadFail;
  }
  if (pgfi_cleanup(&pgfi) && (!reterr)) {
    reterr = kPglRetReadFail;
  }
  // no bigstack_reset() needed?
  return reterr;
}

pglerr_t zst_decompress(const char* in_fname, const char* out_fname) {
  // Since this needs to be able to dump the decompressed data and nothing but
  // the decompressed data to stdout, we have to duplicate a bit of
  // plink2_common code and strip out printing/logging.

  // Strictly speaking, this can decompress gzipped files too, but that's not
  // its purpose.
  gzFile gz_infile = gzopen(in_fname, FOPEN_RB);
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (!gz_infile) {
      fprintf(stderr, g_errstr_fopen, in_fname);
      goto zst_decompress_ret_OPEN_FAIL;
    }
    if (gzbuffer(gz_infile, 131072)) {
      goto zst_decompress_ret_NOMEM;
    }
    if (out_fname) {
      outfile = fopen(out_fname, FOPEN_WB);
      if (!outfile) {
	fprintf(stderr, g_errstr_fopen, out_fname);
	goto zst_decompress_ret_OPEN_FAIL;
      }
    } else {
      outfile = stdout;
    }
    unsigned char* buf = (unsigned char*)g_textbuf;
    while (1) {
      const int32_t bytes_read = gzread(gz_infile, buf, kTextbufMainSize);
      if (bytes_read <= 0) {
	if (!bytes_read) {
	  break;
	}
	goto zst_decompress_ret_READ_FAIL;
      }
      if (!fwrite(buf, bytes_read, 1, outfile)) {
	goto zst_decompress_ret_WRITE_FAIL;
      }
      fflush(outfile);
    }
    if (gzclose_null(&gz_infile)) {
      goto zst_decompress_ret_READ_FAIL;
    }
    if (out_fname) {
      if (fclose_null(&outfile)) {
	goto zst_decompress_ret_WRITE_FAIL;
      }
    }
  }
  // we exit from main() immediately, so need to print nomem/read/write error
  // messages here
  while (0) {
  zst_decompress_ret_NOMEM:
    // in this exceedingly unlikely case, the --memory flag doesn't help, so
    // print a different message
    fputs("Error: Out of memory.\n", stderr);
    reterr = kPglRetNomem;
    break;
  zst_decompress_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  zst_decompress_ret_READ_FAIL:
    fputs(errstr_read, stderr);
    reterr = kPglRetReadFail;
    break;
  zst_decompress_ret_WRITE_FAIL:
    fputs(errstr_write, stderr);
    reterr = kPglRetWriteFail;
    break;
  }
  if (out_fname) {
    fclose_cond(outfile);
  }
  gzclose_cond(gz_infile);
  return reterr;
}

// useful when there's e.g. a filename and an optional modifier, and we want to
// permit either parmeter ordering
boolerr_t check_extra_param(char** argv, const char* permitted_modif, uint32_t* other_idx_ptr) {
  const uint32_t idx_base = *other_idx_ptr;
  if (!strcmp(argv[idx_base], permitted_modif)) {
    *other_idx_ptr = idx_base + 1;
  } else if (strcmp(argv[idx_base + 1], permitted_modif)) {
    LOGERRPRINTF("Error: Invalid %s parameter sequence.\n", argv[0]);
    return 1;
  }
  return 0;
}

char extract_char_param(const char* ss) {
  // maps c, 'c', and "c" to c, and anything else to the null char.  This is
  // intended to support e.g. always using '#' to designate a # parameter
  // without worrying about differences between shells.
  const char cc = ss[0];
  if (((cc == '\'') || (cc == '"')) && (ss[1]) && (ss[2] == cc) && (!ss[3])) {
    return ss[1];
  }
  if (cc && (!ss[1])) {
    return cc;
  }
  return '\0';
}

pglerr_t cmdline_alloc_string(const char* source, const char* flag_name, uint32_t max_slen, char** sbuf_ptr) {
  const uint32_t slen = strlen(source);
  if (slen > max_slen) {
    LOGERRPRINTF("Error: %s parameter too long.\n", flag_name);
    return kPglRetInvalidCmdline;
  }
  const uint32_t blen = slen + 1;
  if (pgl_malloc(blen, sbuf_ptr)) {
    return kPglRetNomem;
  }
  memcpy(*sbuf_ptr, source, blen);
  return kPglRetSuccess;
}

pglerr_t alloc_fname(const char* source, const char* flagname_p, uint32_t extra_size, char** fnbuf_ptr) {
  const uint32_t blen = strlen(source) + 1;
  if (blen > (kPglFnamesize - extra_size)) {
    LOGERRPRINTF("Error: --%s filename too long.\n", flagname_p);
    return kPglRetOpenFail;
  }
  if (pgl_malloc(blen + extra_size, fnbuf_ptr)) {
    return kPglRetNomem;
  }
  memcpy(*fnbuf_ptr, source, blen);
  return kPglRetSuccess;
}

pglerr_t alloc_and_flatten(char** sources, uint32_t param_ct, uint32_t max_blen, char** flattened_buf_ptr) {
  uintptr_t tot_blen = 1;
  for (uint32_t param_idx = 0; param_idx < param_ct; ++param_idx) {
    const uint32_t cur_blen = 1 + strlen(sources[param_idx]);
    if (cur_blen > max_blen) {
      return kPglRetInvalidCmdline;
    }
    tot_blen += cur_blen;
  }
  char* buf_iter;
  if (pgl_malloc(tot_blen, &buf_iter)) {
    return kPglRetNomem;
  }
  *flattened_buf_ptr = buf_iter;
  for (uint32_t param_idx = 0; param_idx < param_ct; ++param_idx) {
    buf_iter = strcpyax(buf_iter, sources[param_idx], '\0');
  }
  *buf_iter = '\0';
  return kPglRetSuccess;
}


// may move these to plink2_common or plink2_filter
char* parse_next_binary_op(char* expr_str, uint32_t expr_slen, char** op_start_ptr, cmp_binary_op_t* binary_op_ptr) {
  // !=, <>: kCmpOperatorNoteq
  // <: kCmpOperatorLe
  // <=: kCmpOperatorLeq
  // =, ==: kCmpOperatorEq
  // >=: kCmpOperatorGeq
  // >: kCmpOperatorGe
  char* next_eq = (char*)memchr(expr_str, '=', expr_slen);
  char* next_lt = (char*)memchr(expr_str, '<', expr_slen);
  char* next_gt = (char*)memchr(expr_str, '>', expr_slen);
  if (!next_eq) {
    if (!next_lt) {
      if (!next_gt) {
	return nullptr;
      }
      *op_start_ptr = next_gt;
      *binary_op_ptr = kCmpOperatorGe;
      return &(next_gt[1]);
    }
    if (next_gt == (&(next_lt[1]))) {
      *op_start_ptr = next_lt;
      *binary_op_ptr = kCmpOperatorNoteq;
      return &(next_lt[2]);
    }
    if ((!next_gt) || (next_gt > next_lt)) {
      *op_start_ptr = next_lt;
      *binary_op_ptr = kCmpOperatorLe;
      return &(next_lt[1]);
    }
    *op_start_ptr = next_gt;
    *binary_op_ptr = kCmpOperatorGe;
    return &(next_gt[1]);
  }
  if ((!next_lt) || (next_lt > next_eq)) {
    if ((!next_gt) || (next_gt > next_eq)) {
      if ((next_eq != expr_str) && (next_eq[-1] == '!')) {
	*op_start_ptr = &(next_eq[-1]);
	*binary_op_ptr = kCmpOperatorNoteq;
	return &(next_eq[1]);
      }
      *op_start_ptr = next_eq;
      *binary_op_ptr = kCmpOperatorEq;
      return (next_eq[1] == '=')? (&(next_eq[2])) : (&(next_eq[1]));
    }
    *op_start_ptr = next_gt;
    if (next_eq == (&(next_gt[1]))) {
      *binary_op_ptr = kCmpOperatorGeq;
      return &(next_gt[2]);
    }
    *binary_op_ptr = kCmpOperatorGe;
    return &(next_gt[1]);
  }
  if (next_gt == (&(next_lt[1]))) {
    *op_start_ptr = next_lt;
    *binary_op_ptr = kCmpOperatorNoteq;
    return &(next_lt[2]);
  }
  if ((!next_gt) || (next_gt > next_lt)) {
    *op_start_ptr = next_lt;
    if (next_eq == (&(next_lt[1]))) {
      *binary_op_ptr = kCmpOperatorLeq;
      return &(next_lt[2]);
    }
    *binary_op_ptr = kCmpOperatorLe;
    return &(next_lt[1]);
  }
  *op_start_ptr = next_gt;
  if (next_eq == (&(next_gt[1]))) {
    *binary_op_ptr = kCmpOperatorGeq;
    return &(next_gt[2]);
  }
  *binary_op_ptr = kCmpOperatorGe;
  return &(next_gt[1]);
}

pglerr_t validate_and_alloc_cmp_expr(char** sources, const char* flag_name, uint32_t param_ct, cmp_expr_t* cmp_expr_ptr) {
  // restrict to [pheno/covar name] [operator] [pheno val] for now.  could
  // support or/and, parentheses, etc. later.
  pglerr_t reterr = kPglRetSuccess;
  {
    if ((param_ct != 1) && (param_ct != 3)) {
      goto validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC;
    }
    char* pheno_name_start = sources[0];
    char* pheno_val_start;
    uint32_t pheno_name_slen;
    uint32_t pheno_val_slen;
    if (param_ct == 3) {
      pheno_name_slen = strlen(pheno_name_start);
      char* op_str = sources[1];
      uint32_t op_slen = strlen(op_str);
      // ok to have single/double quotes around operator
      if (op_slen > 2) {
	const char cc = op_str[0];
	if (((cc == '\'') || (cc == '"')) && (op_str[op_slen - 1] == cc)) {
	  ++op_str;
	  op_slen -= 2;
	}
      }
      char* op_start;
      char* op_end = parse_next_binary_op(op_str, op_slen, &op_start, &cmp_expr_ptr->binary_op);
      if ((!op_end) || (*op_end) || (op_start != op_str)) {
	goto validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC;
      }
      pheno_val_start = sources[2];
      pheno_val_slen = strlen(pheno_val_start);
    } else {
      // permit param_ct == 1 as long as tokens are unambiguous
      uint32_t expr_slen = strlen(pheno_name_start);
      char* op_start;
      pheno_val_start = parse_next_binary_op(pheno_name_start, expr_slen, &op_start, &cmp_expr_ptr->binary_op);
      if ((!pheno_val_start) || (!(*pheno_val_start)) || (op_start == pheno_name_start)) {
        goto validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC;
      }
      pheno_name_slen = (uintptr_t)(op_start - pheno_name_start);
      pheno_val_slen = expr_slen - ((uintptr_t)(pheno_val_start - pheno_name_start));
    }
    if ((pheno_name_slen > kMaxIdSlen) || (pheno_val_slen > kMaxIdSlen)) {
      LOGERRPRINTF("Error: ID too long in %s expression.\n", flag_name);
      goto validate_and_alloc_cmp_expr_ret_INVALID_CMDLINE;
    }
    char* new_pheno_name_buf;
    if (pgl_malloc(2 + pheno_name_slen + pheno_val_slen, &new_pheno_name_buf)) {
      goto validate_and_alloc_cmp_expr_ret_NOMEM;
    }
    memcpyx(new_pheno_name_buf, pheno_name_start, pheno_name_slen, '\0');
    // pheno_val_start guaranteed to be null-terminated for now
    memcpy(&(new_pheno_name_buf[pheno_name_slen + 1]), pheno_val_start, pheno_val_slen + 1);
    cmp_expr_ptr->pheno_name = new_pheno_name_buf;
  }
  while (0) {
  validate_and_alloc_cmp_expr_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC:
    LOGERRPRINTF("Error: Invalid %s expression.\n", flag_name);
  validate_and_alloc_cmp_expr_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
  return reterr;
}

/*
pglerr_t alloc_and_flatten_comma_delim(char** sources, uint32_t param_ct, char** flattened_buf_ptr) {
  uint32_t totlen = 1;
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
      totlen += 1 + (uintptr_t)(cur_token_end - cur_param_iter);
      cur_param_iter = &(cur_token_end[1]);
    }
    totlen += 1 + strlen(cur_param_iter);
  }
  char* write_iter;
  if (pgl_malloc(totlen, &write_iter)) {
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
      write_iter = memcpyax(write_iter, cur_param_iter, (uintptr_t)(cur_token_end - cur_param_iter), '\0');
      cur_param_iter = &(cur_token_end[1]);
    }
    write_iter = strcpyax(write_iter, cur_param_iter, '\0');
  }
  *write_iter = '\0';
  return kPglRetSuccess;
}
*/

void invalid_arg(const char* cur_arg) {
  LOGPREPRINTFWW("Error: Unrecognized flag ('%s').\n", cur_arg);
}

void print_ver() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}

pglerr_t rerun(uint32_t rerun_argv_pos, uint32_t rerun_parameter_present, int32_t* argc_ptr, uint32_t* first_arg_idx_ptr, char*** argv_ptr, char*** subst_argv_ptr, char** rerun_buf_ptr) {
  // caller is responsible for freeing rerun_buf
  char** subst_argv2 = nullptr;
  uintptr_t line_idx = 1;
  pglerr_t reterr = kPglRetSuccess;
  gzFile gz_rerunfile;
  {
    char** argv = *argv_ptr;
    gz_rerunfile = gzopen(rerun_parameter_present? argv[rerun_argv_pos + 1] : (PROG_NAME_STR ".log"), FOPEN_RB);
    if (!gz_rerunfile) {
      goto rerun_ret_OPEN_FAIL;
    }
    char* textbuf = g_textbuf;
    textbuf[kMaxMediumLine - 1] = ' ';
    if (!gzgets(gz_rerunfile, textbuf, kMaxMediumLine)) {
      print_ver();
      fputs("Error: Empty log file for --rerun.\n", stderr);
      goto rerun_ret_MALFORMED_INPUT;
    }
    if (!textbuf[kMaxMediumLine - 1]) {
      goto rerun_ret_LONG_LINE;
    }
    if (!gzgets(gz_rerunfile, textbuf, kMaxMediumLine)) {
      print_ver();
      fputs("Error: Only one line in --rerun log file.\n", stderr);
      goto rerun_ret_MALFORMED_INPUT;
    }
    line_idx++;
    if (!textbuf[kMaxMediumLine - 1]) {
      goto rerun_ret_LONG_LINE;
    }
    // don't bother supporting "xx arguments: --aa bb --cc --dd" format
    while (memcmp(textbuf, "Options in effect:", 18) || (textbuf[18] >= ' ')) {
      line_idx++;
      if (!gzgets(gz_rerunfile, textbuf, kMaxMediumLine)) {
	print_ver();
	fputs("Error: Invalid log file for --rerun.\n", stderr);
	goto rerun_ret_MALFORMED_INPUT;
      }
    }
    char* all_args_write_iter = textbuf;
    char* textbuf_limit = &(textbuf[kMaxMediumLine]);
    uint32_t loaded_arg_ct = 0;
    // We load each of the option lines in sequence into textbuf, always
    // overwriting the previous line's newline.  (Note that textbuf[] has
    // size > 2 * kMaxMediumLine; this lets us avoid additional
    // dynamic memory allocation as long as we impose the constraint that all
    // lines combined add up to less than kMaxMediumLine.)
    while (1) {
      all_args_write_iter[kMaxMediumLine - 1] = ' ';
      if (!gzgets(gz_rerunfile, all_args_write_iter, kMaxMediumLine)) {
	break;
      }
      line_idx++;
      if (!all_args_write_iter[kMaxMediumLine - 1]) {
	goto rerun_ret_LONG_LINE;
      }
      char* arg_iter = skip_initial_spaces(all_args_write_iter);
      if (is_eoln_kns(*arg_iter)) {
	*all_args_write_iter = '\0';
	break;
      }
      char* token_end;
      do {
	token_end = token_endnn(arg_iter);
	loaded_arg_ct++;
	arg_iter = skip_initial_spaces(token_end);
      } while (!is_eoln_kns(*arg_iter));
      all_args_write_iter = token_end;
      if (all_args_write_iter >= textbuf_limit) {
	print_ver();
	fputs("Error: --rerun argument sequence too long.\n", stderr);
	goto rerun_ret_MALFORMED_INPUT;
      }
    }
    gzclose_null(&gz_rerunfile);
    const uint32_t line_byte_ct = 1 + (uintptr_t)(all_args_write_iter - textbuf);
    char* rerun_buf;
    if (pgl_malloc(line_byte_ct, &rerun_buf)) {
      goto rerun_ret_NOMEM;
    }
    *rerun_buf_ptr = rerun_buf;
    memcpy(rerun_buf, textbuf, line_byte_ct);
    const uint32_t argc = (uint32_t)(*argc_ptr);
    const uint32_t first_arg_idx = *first_arg_idx_ptr;
    char* rerun_first_token = skip_initial_spaces(rerun_buf);
    char* arg_iter = rerun_first_token;
    // now use textbuf as a lame bitfield
    memset(textbuf, 1, loaded_arg_ct);
    uint32_t loaded_arg_idx = 0;
    uint32_t duplicate_arg_ct = 0;
    do {
      if (no_more_tokens_kns(arg_iter)) {
	print_ver();
	fputs("Error: Line 2 of --rerun log file has fewer tokens than expected.\n", stderr);
	goto rerun_ret_MALFORMED_INPUT;
      }
      char* flagname_p = is_flag_start(arg_iter);
      if (flagname_p) {
	const uint32_t slen = strlen_se(flagname_p);
	uint32_t cmdline_arg_idx = first_arg_idx;
	for (; cmdline_arg_idx < argc; cmdline_arg_idx++) {
	  char* later_flagname_p = is_flag_start(argv[cmdline_arg_idx]);
	  if (later_flagname_p) {
	    const uint32_t slen2 = strlen(later_flagname_p);
	    if ((slen == slen2) && (!memcmp(flagname_p, later_flagname_p, slen))) {
	      cmdline_arg_idx = 0xffffffffU;
	      break;
	    }
	  }
	}
	if (cmdline_arg_idx == 0xffffffffU) {
	  // matching flag, override --rerun
	  do {
	    duplicate_arg_ct++;
	    textbuf[loaded_arg_idx++] = 0;
	    if (loaded_arg_idx == loaded_arg_ct) {
	      break;
	    }
	    arg_iter = next_token(arg_iter);
	  } while (!is_flag(arg_iter));
	} else {
	  loaded_arg_idx++;
	  arg_iter = next_token(arg_iter);
	}
      } else {
	loaded_arg_idx++;
	arg_iter = next_token(arg_iter);
      }
    } while (loaded_arg_idx < loaded_arg_ct);
    if (pgl_malloc((argc + loaded_arg_ct - duplicate_arg_ct - rerun_parameter_present - 1 - first_arg_idx) * sizeof(intptr_t), &subst_argv2)) {
      goto rerun_ret_NOMEM;
    }
    uint32_t new_arg_idx = rerun_argv_pos - first_arg_idx;
    memcpy(subst_argv2, &(argv[first_arg_idx]), new_arg_idx * sizeof(intptr_t));
    arg_iter = rerun_first_token;
    for (loaded_arg_idx = 0; loaded_arg_idx < loaded_arg_ct; ++loaded_arg_idx) {
      arg_iter = skip_initial_spaces(arg_iter);
      char* token_end = token_endnn(arg_iter);
      if (textbuf[loaded_arg_idx]) {
	subst_argv2[new_arg_idx++] = arg_iter;
	*token_end = '\0';
      }
      arg_iter = &(token_end[1]);
    }
    const uint32_t final_copy_start_idx = rerun_argv_pos + rerun_parameter_present + 1;
    memcpy(&(subst_argv2[new_arg_idx]), &(argv[final_copy_start_idx]), (argc - final_copy_start_idx) * sizeof(intptr_t));
    *first_arg_idx_ptr = 0;
    *argc_ptr = new_arg_idx + argc - final_copy_start_idx;
    if (*subst_argv_ptr) {
      free(*subst_argv_ptr);
    }
    *subst_argv_ptr = subst_argv2;
    *argv_ptr = subst_argv2;
    subst_argv2 = nullptr;
  }
  while (0) {
  rerun_ret_NOMEM:
    print_ver();
    reterr = kPglRetNomem;
    break;
  rerun_ret_OPEN_FAIL:
    print_ver();
    reterr = kPglRetOpenFail;
    break;
  rerun_ret_LONG_LINE:
    print_ver();
    fprintf(stderr, "Error: Line %" PRIuPTR " of --rerun log file is pathologically long.\n", line_idx);
  rerun_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  free_cond(subst_argv2);
  gzclose_cond(gz_rerunfile);
  return reterr;
}

uint32_t cmdline_single_chr(const chr_info_t* cip, misc_flags_t misc_flags) {
  if ((misc_flags & (kfMiscAutosomeOnly | kfMiscAutosomePar)) || (!cip->is_include_stack)) {
    return 0;
  }
  const uint32_t main_chr_ct = popcount_longs(cip->chr_mask, kChrExcludeWords) + popcount_long(cip->chr_mask[kChrMaskWords - 1]);
  if (main_chr_ct > 1) {
    return 0;
  }
  if (main_chr_ct == 1) {
    return (cip->incl_excl_name_stack == nullptr);
  }
  return cip->incl_excl_name_stack && (!(cip->incl_excl_name_stack->next));
}

pglerr_t parse_col_descriptor(const char* col_descriptor_iter, const char* supported_ids, const char* cur_flag_name, uint32_t first_col_shifted, uint32_t default_cols_mask, uint32_t prohibit_empty, void* result_ptr) {
  // col_descriptor is usually a pointer to argv[...][5] (first five characters
  // are "cols=").  supported_ids is a multistr.
  // may need to switch first_col_shifted/default_cols_mask/result to uint64_t
  pglerr_t reterr = kPglRetSuccess;
  uint32_t* id_map = nullptr;
  {
    uint32_t max_id_blen = 0;
    uint32_t id_ct = 0;

    // work around strchr not returning a const char*?
    const char* supported_ids_iter = supported_ids;
    
    // can precompute this sorted index and avoid the dynamic
    // allocations/deallocations, but this is cheap enough that I'd rather make
    // it easier to extend functionality.
    do {
      const char* tok_end = (const char*)rawmemchr(supported_ids_iter, '\0');
      const uint32_t slen = (uintptr_t)(tok_end - supported_ids_iter);
      if (slen >= max_id_blen) {
        max_id_blen = slen + 1;
      }
      ++id_ct;
      supported_ids_iter = &(tok_end[1]);
    } while (*supported_ids_iter);
    // max_id_blen + 4 extra bytes at the end, to support a "maybe" search
    // (yes, this can also be precomputed)
    if (pgl_malloc((max_id_blen + 4) * (id_ct + 1), &id_map)) {
      goto parse_col_descriptor_ret_NOMEM;
    }
    char* sorted_ids = (char*)(&(id_map[id_ct]));
    supported_ids_iter = (const char*)supported_ids;
    for (uint32_t id_idx = 0; id_idx < id_ct; ++id_idx) {
      const uint32_t blen = strlen(supported_ids_iter) + 1;
      memcpy(&(sorted_ids[id_idx * max_id_blen]), supported_ids_iter, blen);
      id_map[id_idx] = id_idx;
      supported_ids_iter = &(supported_ids_iter[blen]);
    }
    if (sort_strbox_indexed_malloc(id_ct, max_id_blen, sorted_ids, id_map)) {
      goto parse_col_descriptor_ret_NOMEM;
    }
    uint32_t result = *((uint32_t*)result_ptr);
    // might not want to bother splitting this into two loops
    if ((col_descriptor_iter[0] == '+') || (col_descriptor_iter[0] == '-')) {
      result |= default_cols_mask;
      char* maybebuf = &(sorted_ids[max_id_blen * id_ct]);
      memcpy(maybebuf, "maybe", 5);
      while (1) {
	const char* id_start = &(col_descriptor_iter[1]);
	const char* tok_end = strchr(id_start, ',');
	uint32_t slen;
	if (!tok_end) {
	  slen = strlen(id_start);
	} else {
	  slen = (uintptr_t)(tok_end - id_start);
	}
	int32_t alpha_idx = bsearch_str(id_start, sorted_ids, slen, max_id_blen, id_ct);
	if (alpha_idx == -1) {
	  char* write_iter = strcpya(g_logbuf, "Error: Unrecognized ID '");
	  write_iter = memcpya(write_iter, id_start, slen);
	  write_iter = strcpya(write_iter, "' in --");
	  write_iter = strcpya(write_iter, cur_flag_name);
	  write_iter = strcpya(write_iter, " column set descriptor.\n");
	  goto parse_col_descriptor_ret_INVALID_CMDLINE_WW;
	}
	uint32_t shift = id_map[(uint32_t)alpha_idx];
	if (col_descriptor_iter[0] == '+') {
	  result |= first_col_shifted << shift;
	} else {
	  if (result & (first_col_shifted << shift)) {
	    result -= first_col_shifted << shift;
	  } else if (slen + 5 < max_id_blen) {
	    // special case: if default column set includes e.g. "maybesid",
	    // and user types "-sid", that should work
	    memcpy(&(maybebuf[5]), id_start, slen);
	    alpha_idx = bsearch_str(maybebuf, sorted_ids, slen + 5, max_id_blen, id_ct);
	    if (alpha_idx != -1) {
	      shift = id_map[(uint32_t)alpha_idx];
	      result &= ~(first_col_shifted << shift);
	    }
	  }
	}
	if (!tok_end) {
	  break;
	}
	col_descriptor_iter = &(tok_end[1]);
	if ((col_descriptor_iter[0] != '+') && (col_descriptor_iter[0] != '-')) {
	  goto parse_col_descriptor_ret_MIXED_SIGN;
	}
      }
    } else if (*col_descriptor_iter) {
      while (1) {
	const char* tok_end = strchr(col_descriptor_iter, ',');
	uint32_t slen;
	if (!tok_end) {
	  slen = strlen(col_descriptor_iter);
	} else {
	  slen = (uintptr_t)(tok_end - col_descriptor_iter);
	}
	int32_t alpha_idx = bsearch_str(col_descriptor_iter, sorted_ids, slen, max_id_blen, id_ct);
	if (alpha_idx == -1) {
	  char* write_iter = strcpya(g_logbuf, "Error: Unrecognized ID '");
	  write_iter = memcpya(write_iter, col_descriptor_iter, slen);
	  write_iter = strcpya(write_iter, "' in --");
	  write_iter = strcpya(write_iter, cur_flag_name);
	  write_iter = strcpya(write_iter, " column set descriptor.\n");
	  goto parse_col_descriptor_ret_INVALID_CMDLINE_WW;
	}
	uint32_t shift = id_map[(uint32_t)alpha_idx];
	result |= first_col_shifted << shift;
	if (!tok_end) {
	  break;
	}
	col_descriptor_iter = &(tok_end[1]);
	if ((col_descriptor_iter[0] == '+') || (col_descriptor_iter[0] == '-')) {
	  goto parse_col_descriptor_ret_MIXED_SIGN;
	}
      }
    }
    if (prohibit_empty && (!(result & (first_col_shifted * (0xffffffffU >> (32 - id_ct)))))) {
      char* write_iter = strcpya(g_logbuf, "Error: All columns excluded by --");
      write_iter = strcpya(write_iter, cur_flag_name);
      write_iter = strcpya(write_iter, " column set descriptor.\n");
      goto parse_col_descriptor_ret_INVALID_CMDLINE_WW;
    }
    *((uint32_t*)result_ptr) = result;
  }
  while (0) {
  parse_col_descriptor_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  parse_col_descriptor_ret_MIXED_SIGN:
    sprintf(g_logbuf, "Error: Invalid --%s column set descriptor (either all column set IDs must be preceded by +/-, or none of them can be).\n", cur_flag_name);
  parse_col_descriptor_ret_INVALID_CMDLINE_WW:
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetInvalidCmdline;
    break;
  }
  free_cond(id_map);
  return reterr;
}

void get_exportf_targets(char** argv, uint32_t param_ct, exportf_flags_t* exportf_modifier_ptr, idpaste_t* exportf_id_paste_ptr, uint32_t* format_param_idxs_ptr) {
  // does not error out if no format present, since this is needed for --recode
  // translation
  // supports multiple formats
  uint32_t format_param_idxs = 0;
  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
    const char* cur_modif = argv[param_idx];
    const char* cur_modif2 = &(cur_modif[1]);
    exportf_flags_t cur_format = kfExportf0;
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
	  sprintf(g_logbuf, "Note: --export 'v%s' modifier is deprecated.  Use 'vcf' + 'id-paste=%s'.\n", cur_modif2, &(cur_modif2[3]));
	  cur_format = kfExportfVcf;
	  *exportf_id_paste_ptr = (cur_modif2[3] == 'f')? kfIdpasteFid : kfIdpasteIid;
	}
      }
      break;
    }
    if (cur_format) {
      format_param_idxs |= 1U << param_idx;
      *exportf_modifier_ptr |= cur_format;
    }
  }
  *format_param_idxs_ptr = format_param_idxs;
}

uint32_t varid_template_is_valid(char* varid_str, const char* flagname_p) {
  char* sptr = strchr(varid_str, '@');
  char* sptr2 = strchr(varid_str, '#');
  if ((!sptr) || (!sptr2) || strchr(&(sptr[1]), '@') || strchr(&(sptr2[1]), '#')) {
    LOGERRPRINTFWW("Error: The --%s template string requires exactly one '@' and one '#'.\n", flagname_p);
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
    uint32_t first_allele_type_code = (unsigned char)(*sptr2);
    if ((first_allele_type_code == 49) || (first_allele_type_code == 50)) {
      sptr2 = strchr(sptr2, '$');
      if ((!sptr2) || strchr(&(sptr2[1]), '$') || ((first_allele_type_code + ((unsigned char)sptr2[1])) != 99)) {
      varid_template_is_valid_dollar_error:
	LOGERRPRINTFWW("Error: The --%s template string requires either no instances of '$', exactly one instance of '$r' and/or '$a', or exactly one '$1' and one '$2'.\n", flagname_p);
	return 0;
      }
    } else {
      first_allele_type_code &= 0xdf;
      if ((first_allele_type_code != 65) && (first_allele_type_code != 82)) {
	goto varid_template_is_valid_dollar_error;
      }
      sptr2 = strchr(sptr2, '$');
      if (sptr2) {
	const uint32_t second_allele_type_code = (uint32_t)((unsigned char)(*(++sptr2))) & 0xdf;
	if (((first_allele_type_code + second_allele_type_code) != 147) || strchr(sptr2, '$')) {
	  goto varid_template_is_valid_dollar_error;
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
} // namespace plink2
#endif

int main(int argc, char** argv) {
#ifdef __cplusplus
  using namespace plink2;
#endif
  // special case, since it may dump to stdout
  if ((argc > 1) && ((!strcmp(argv[1], "--zst-decompress")) || (!strcmp(argv[1], "-zst-decompress")))) {
    if (argc == 2) {
      fprintf(stderr, "Error: Missing %s parameter.\n", argv[1]);
      return (uint32_t)kPglRetInvalidCmdline;
    }
    for (int ii = 2; ii < argc; ++ii) {
      if (is_flag(argv[(uint32_t)ii])) {
	fprintf(stderr, "Error: %s cannot be used with other flags.\n", argv[1]);
	return (uint32_t)kPglRetInvalidCmdline;
      }
    }
    if (argc > 4) {
      fprintf(stderr, "Error: %s accepts at most 2 parameters.\n", argv[1]);
      return (uint32_t)kPglRetInvalidCmdline;
    }
    return (uint32_t)zst_decompress(argv[2], (argc == 4)? argv[3] : nullptr);
  }
  
  unsigned char* bigstack_ua = nullptr;
  char** subst_argv = nullptr;
  char* script_buf = nullptr;
  char* rerun_buf = nullptr;
  char* flag_buf = nullptr;
  char* flagname_p = nullptr;
  uint32_t* flag_map = nullptr;
  char* king_cutoff_fprefix = nullptr;
  char* const_fid = nullptr;
  char* var_filter_exceptions_flattened = nullptr;
  char* require_pheno_flattened = nullptr;
  char* require_covar_flattened = nullptr;
  char* import_single_chr_str = nullptr;
  char* ox_missing_code = nullptr;
  char* vcf_dosage_import_field = nullptr;
  FILE* scriptfile = nullptr;
  uint32_t* rseeds = nullptr;
  ll_str_t* file_delete_list = nullptr;
  uint32_t arg_idx = 0;
  uint32_t print_end_time = 0;
  uint32_t warning_errcode = 0;
  pglerr_t reterr = kPglRetSuccess;
  plink2_cmdline_t pc;
  pc.filter_flags = kfFilter0;
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
  pc.update_sex_fname = nullptr;
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
  pc.vstd_flattened = nullptr;
  pc.quantnorm_flattened = nullptr;
  pc.covar_quantnorm_flattened = nullptr;
  init_range_list(&pc.snps_range_list);
  init_range_list(&pc.exclude_snps_range_list);
  init_range_list(&pc.pheno_range_list);
  init_range_list(&pc.covar_range_list);
  init_ld(&pc.ld_info);
  init_glm(&pc.glm_info);
  init_adjust(&pc.adjust_info);
  init_score(&pc.score_info);
  init_cmp_expr(&pc.keep_if_expr);
  init_cmp_expr(&pc.remove_if_expr);
  chr_info_t chr_info;
  if (init_chr_info(&chr_info)) {
    goto main_ret_NOMEM_NOLOG;
  }
  
  {
    // standardize strtod() behavior
    // setlocale(LC_NUMERIC, "C");

    uint32_t first_arg_idx = 1;
    for (arg_idx = 1; arg_idx < (uint32_t)argc; ++arg_idx) {
      if ((!strcmp("-script", argv[arg_idx])) || (!strcmp("--script", argv[arg_idx]))) {
	const uint32_t param_ct = param_count(argv, argc, arg_idx);
	if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	  print_ver();
	  fputs(g_logbuf, stderr);
	  fputs(errstr_append, stderr);
	  goto main_ret_INVALID_CMDLINE;
	}
	for (uint32_t arg_idx2 = arg_idx + 2; arg_idx2 < (uint32_t)argc; ++arg_idx2) {
	  if ((!strcmp("-script", argv[arg_idx2])) || (!strcmp("--script", argv[arg_idx2]))) {
	    print_ver();
	    fputs("Error: Multiple --script flags.  Merge the files into one.\n", stderr);
	    fputs(errstr_append, stderr);
	    goto main_ret_INVALID_CMDLINE;
	  }
	}
	// logging not yet active, so don't use fopen_checked()
	scriptfile = fopen(argv[arg_idx + 1], FOPEN_RB);
	if (!scriptfile) {
	  print_ver();
	  fprintf(stderr, g_errstr_fopen, argv[arg_idx + 1]);
	  goto main_ret_OPEN_FAIL;
	}
	if (fseeko(scriptfile, 0, SEEK_END)) {
	  goto main_ret_READ_FAIL_NOLOG;
	}
	int64_t fsize = ftello(scriptfile);
	if (fsize < 0) {
	  goto main_ret_READ_FAIL_NOLOG;
	}
	if (fsize > 0x7ffffffe) {
	  // could actually happen if user enters parameters in the wrong
	  // order, so may as well catch it and print a somewhat informative
	  // error message
	  print_ver();
	  fputs("Error: --script file too large.", stderr);
	  goto main_ret_INVALID_CMDLINE;
	}
	rewind(scriptfile);
	const uint32_t fsize_ui = (uint64_t)fsize;
	if (pgl_malloc(fsize_ui + 1, &script_buf)) {
	  goto main_ret_NOMEM_NOLOG;
	}
	if (!fread(script_buf, fsize_ui, 1, scriptfile)) {
	  goto main_ret_READ_FAIL_NOLOG;
	}
	script_buf[fsize_ui] = '\0';
	fclose_null(&scriptfile);
	uint32_t num_script_params = 0;
	char* script_buf_iter = script_buf;
	uint32_t char_code;
	do {
	  uint32_t char_code_m1;
	  do {
	    char_code_m1 = ((uint32_t)((unsigned char)(*script_buf_iter++))) - 1;
	  } while (char_code_m1 < 32);
	  if (char_code_m1 == 0xffffffffU) {
	    break;
	  }
	  ++num_script_params;
	  do {
	    char_code = (uint32_t)((unsigned char)(*script_buf_iter++));
	  } while (char_code > 32);
	} while (char_code);
	if (script_buf_iter != (&(script_buf[fsize_ui + 1]))) {
	  print_ver();
	  fputs("Error: Null byte in --script file.\n", stderr);
	  goto main_ret_INVALID_CMDLINE;
	}
	const uint32_t new_param_ct = num_script_params + argc - 3;
	if (pgl_malloc(new_param_ct * sizeof(intptr_t), &subst_argv)) {
	  goto main_ret_NOMEM_NOLOG;
	}
	memcpy(subst_argv, &(argv[1]), arg_idx * sizeof(intptr_t));
	const uint32_t load_param_idx_end = arg_idx + num_script_params;
	script_buf_iter = &(script_buf[-1]);
	for (uint32_t param_idx = arg_idx; param_idx < load_param_idx_end; ++param_idx) {
	  while (((unsigned char)(*(++script_buf_iter))) <= 32);
	  subst_argv[param_idx] = script_buf_iter;
	  while (((unsigned char)(*(++script_buf_iter))) > 32);
	  // could enforce some sort of length limit here
	  *script_buf_iter = '\0';
	}
	memcpy(&(subst_argv[load_param_idx_end]), &(argv[arg_idx + 2]), (argc - arg_idx - 2) * sizeof(intptr_t));
	argc = new_param_ct;
	first_arg_idx = 0;
	argv = subst_argv;
	break;
      }
    }
    for (arg_idx = first_arg_idx; arg_idx < (uint32_t)argc; ++arg_idx) {
      if ((!strcmp("-rerun", argv[arg_idx])) || (!strcmp("--rerun", argv[arg_idx]))) {
	const uint32_t param_ct = param_count(argv, argc, arg_idx);
	if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	  print_ver();
	  fputs(g_logbuf, stderr);
	  fputs(errstr_append, stderr);
	  goto main_ret_INVALID_CMDLINE;
	}
	for (uint32_t arg_idx2 = arg_idx + param_ct + 1; arg_idx2 < (uint32_t)argc; ++arg_idx2) {
	  if ((!strcmp("-rerun", argv[arg_idx2])) || (!strcmp("--rerun", argv[arg_idx2]))) {
	    print_ver();
	    fputs("Error: Duplicate --rerun flag.\n", stderr);
	    goto main_ret_INVALID_CMDLINE;
	  }
	}
	reterr = rerun(arg_idx, param_ct, &argc, &first_arg_idx, &argv, &subst_argv, &rerun_buf);
	if (reterr) {
	  goto main_ret_NOLOG;
	}
	break;
      }
    }
    if ((first_arg_idx < (uint32_t)argc) && (!is_flag(argv[first_arg_idx]))) {
      fputs("Error: First parameter must be a flag.\n", stderr);
      fputs(errstr_append, stderr);
      goto main_ret_INVALID_CMDLINE;
    }
    uint32_t flag_ct = 0;
    uint32_t version_present = 0;
    uint32_t silent_present = 0;
    for (arg_idx = first_arg_idx; arg_idx < (uint32_t)argc; ++arg_idx) {
      flagname_p = is_flag_start(argv[arg_idx]);
      if (flagname_p) {
	if (!strcmp("help", flagname_p)) {
	  print_ver();
	  if ((!first_arg_idx) || (arg_idx != 1) || subst_argv) {
	    fputs("--help present, ignoring other flags.\n", stdout);
	  }
	  if ((arg_idx == ((uint32_t)argc) - 1) && flag_ct) {
	    // make "plink [valid flags/parameters] --help" work, and skip the
	    // parameters
	    char** help_argv;
	    if (pgl_malloc(flag_ct * sizeof(intptr_t), &help_argv)) {
	      goto main_ret_NOMEM_NOLOG2;
	    }
	    uint32_t arg_idx2 = 0;
	    for (uint32_t flag_idx = 0; flag_idx < flag_ct; ++flag_idx) {
	      while (!is_flag_start(argv[++arg_idx2]));
	      help_argv[flag_idx] = argv[arg_idx2];
	    }
	    reterr = disp_help(flag_ct, help_argv);
	    free(help_argv);
	  } else {
	    reterr = disp_help(argc - arg_idx - 1, &(argv[arg_idx + 1]));
	  }
	  goto main_ret_1;
	}
	if ((!strcmp("h", flagname_p)) || (!strcmp("?", flagname_p))) {
	  // these just act like the no-parameter case
	  print_ver();
	  if ((!first_arg_idx) || (arg_idx != 1) || subst_argv) {
	    printf("-%c present, ignoring other flags.\n", *flagname_p);
	  }
	  fputs(g_cmdline_format_str, stdout);
	  fputs(notestr_null_calc2, stdout);
	  reterr = kPglRetHelp;
	  goto main_ret_1;
	}
	if (!strcmp("version", flagname_p)) {
	  version_present = 1;
	} else if (!strcmp("silent", flagname_p)) {
	  silent_present = 1;
	}
	if (strlen(flagname_p) >= kMaxFlagBlen) {
	  print_ver();
	  // shouldn't be possible for this to overflow the buffer...
	  sprintf(g_logbuf, "Error: Unrecognized flag ('%s').\n", argv[arg_idx]);
	  wordwrapb(0);
	  fputs(g_logbuf, stderr);
	  fputs(errstr_append, stderr);
	  goto main_ret_INVALID_CMDLINE;
	}
	++flag_ct;
      }
    }
    if (version_present) {
      fputs(ver_str, stdout);
      putc_unlocked('\n', stdout);
      goto main_ret_1;
    }
    if (silent_present) {
      if (!freopen("/dev/null", "w", stdout)) {
	fputs("Warning: --silent failed.", stderr);
	g_stderr_written_to = 1;
      }
    }
    print_ver();
    if (!flag_ct) {
      goto main_ret_NULL_CALC_0;
    }
    if (pgl_malloc(flag_ct * kMaxFlagBlen, &flag_buf) ||
	pgl_malloc(flag_ct * sizeof(int32_t), &flag_map)) {
      goto main_ret_NOMEM_NOLOG2;
    }
    char* flagname_write_iter = flag_buf;
    uint32_t cur_flag_idx = 0;
    for (arg_idx = first_arg_idx; arg_idx < (uint32_t)argc; ++arg_idx) {
      flagname_p = is_flag_start(argv[arg_idx]);
      if (flagname_p) {
	const uint32_t flag_slen = strlen(flagname_p);
	switch (*flagname_p) {
	case '\0':
	  // special case, since we reserve empty names for preprocessed flags
	  fputs("Error: Unrecognized flag ('--').\n", stderr);
	  goto main_ret_INVALID_CMDLINE;
	case 'a':
	  if ((flag_slen == 3) && (!memcmp(flagname_p, "aec", 3))) {
	    strcpy(flagname_write_iter, "allow-extra-chr");
	  } else if ((flag_slen == 11) && (!memcmp(flagname_p, "autosome-xy", 11))) {
	    strcpy(flagname_write_iter, "autosome-par");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 'b':
	  if ((flag_slen == 3) && (!memcmp(flagname_p, "bed", 3))) {
	    strcpy(flagname_write_iter, "pgen");
	  } else if ((flag_slen == 3) && (!memcmp(flagname_p, "bim", 3))) {
	    strcpy(flagname_write_iter, "pvar");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 'f':
	  if ((flag_slen == 3) && (!memcmp(flagname_p, "fam", 3))) {
	    strcpy(flagname_write_iter, "psam");
	  } else if ((flag_slen == 12) && (!memcmp(flagname_p, "filter-males", 12))) {
	    strcpy(flagname_write_iter, "keep-males");
	  } else if ((flag_slen == 14) && (!memcmp(flagname_p, "filter-females", 14))) {
	    fputs("Note: --filter-females flag deprecated.  Use --keep-females or --remove-males\ninstead.\n", stdout);
	    strcpy(flagname_write_iter, "remove-males");
	  } else if ((flag_slen == 15) && (!memcmp(flagname_p, "filter-founders", 15))) {
	    strcpy(flagname_write_iter, "keep-founders");
	  } else if ((flag_slen == 17) && (!memcmp(flagname_p, "filter-nonfounders", 17))) {
	    strcpy(flagname_write_iter, "keep-nonfounders");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 'h':
	  if ((flag_slen == 5) && (!memcmp(flagname_p, "hound", 5))) {
	    // the creature type should be Dog.
	    strcpy(flagname_write_iter, "dog");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 'k':
	  if ((flag_slen == 13) && (!memcmp(flagname_p, "keep-clusters", 13))) {
	    fputs("Note: --keep-clusters flag deprecated.  Use --keep-cats instead.\n", stdout);
	    strcpy(flagname_write_iter, "keep-cats");
	  } else if ((flag_slen == 18) && (!memcmp(flagname_p, "keep-cluster-names", 18))) {
	    fputs("Note: --keep-cluster-names flag deprecated.  Use --keep-cat-names instead.\n", stdout);
	    strcpy(flagname_write_iter, "keep-cat-names");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 'l':
	  if (((flag_slen == 6) && (!memcmp(flagname_p, "linear", 6))) || ((flag_slen == 8) && (!memcmp(flagname_p, "logistic", 8)))) {
	    strcpy(flagname_write_iter, "glm");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 'm':
	  if ((flag_slen == 6) && (!memcmp(flagname_p, "min-ac", 6))) {
	    strcpy(flagname_write_iter, "mac");
	  } else if ((flag_slen == 6) && (!memcmp(flagname_p, "max-ac", 6))) {
	    strcpy(flagname_write_iter, "max-mac");
	  } else if ((flag_slen == 10) && (!memcmp(flagname_p, "make-bfile", 10))) {
	    strcpy(flagname_write_iter, "make-bed");
	  } else if ((flag_slen == 11) && (!memcmp(flagname_p, "make-bpfile", 11))) {
	    strcpy(flagname_write_iter, "make-bpgen");
	  } else if ((flag_slen == 10) && (!memcmp(flagname_p, "make-pfile", 10))) {
	    strcpy(flagname_write_iter, "make-pgen");
	  } else if ((flag_slen == 12) && (!memcmp(flagname_p, "missing_code", 12))) {
	    strcpy(flagname_write_iter, "missing-code");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 'n':
	  if ((flag_slen == 11) && (!memcmp(flagname_p, "num_threads", 11))) {
	    strcpy(flagname_write_iter, "threads");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
        case 'p':
	  if ((flag_slen == 5) && (!memcmp(flagname_p, "prune", 5))) {
	    strcpy(flagname_write_iter, "require-pheno");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 'r':
	  if ((flag_slen == 6) && (!memcmp(flagname_p, "recode", 6))) {
	    // special case: translate to "export ped" if no format specified
	    const uint32_t param_ct = param_count(argv, argc, arg_idx);
	    if (param_ct > 4) {
	      fputs("Error: --recode accepts at most 4 parameters.\n", stderr);
	      goto main_ret_INVALID_CMDLINE;
	    }
	    exportf_flags_t dummy;
	    idpaste_t dummy2;
	    uint32_t format_param_idxs;
	    get_exportf_targets(&(argv[arg_idx]), param_ct, &dummy, &dummy2, &format_param_idxs);
	    if (!format_param_idxs) {
	      strcpy(flagname_write_iter, "export ped");
	    } else {
	      strcpy(flagname_write_iter, "export");
	    }
	  } else if ((flag_slen == 15) && (!memcmp(flagname_p, "remove-founders", 15))) {
	    strcpy(flagname_write_iter, "keep-founders");
	  } else if ((flag_slen == 17) && (!memcmp(flagname_p, "remove-nonfounders", 17))) {
	    strcpy(flagname_write_iter, "keep-nonfounders");
	  } else if ((flag_slen == 15) && (!memcmp(flagname_p, "remove-clusters", 15))) {
	    fputs("Note: --remove-clusters flag deprecated.  Use --remove-cats instead.\n", stdout);	    
	    strcpy(flagname_write_iter, "remove-cats");
	  } else if ((flag_slen == 20) && (!memcmp(flagname_p, "remove-cluster-names", 20))) {
	    fputs("Note: --remove-cluster-names flag deprecated.  Use --remove-cat-names instead.\n", stdout);
	    strcpy(flagname_write_iter, "remove-cat-names");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 't':
	  if ((flag_slen == 10) && (!memcmp(flagname_p, "thread-num", 10))) {
	    strcpy(flagname_write_iter, "threads");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	case 'v':
	  if ((flag_slen == 10) && (!memcmp(flagname_p, "vcf-filter", 10))) {
	    strcpy(flagname_write_iter, "var-filter");
	  } else if ((flag_slen == 12) && (!memcmp(flagname_p, "vcf-min-qual", 12))) {
	    strcpy(flagname_write_iter, "var-min-qual");
	  } else {
	    goto main_flag_copy;
	  }
	  break;
	default:
	main_flag_copy:
	  memcpy(flagname_write_iter, flagname_p, flag_slen + 1);
	}
	flagname_write_iter = &(flagname_write_iter[kMaxFlagBlen]);
	flag_map[cur_flag_idx++] = arg_idx;
      }
    }
    reterr = sort_cmdline_flags(kMaxFlagBlen, flag_ct, flag_buf, flag_map);
    if (reterr) {
      if (reterr == kPglRetNomem) {
	goto main_ret_NOMEM_NOLOG2;
      }
      goto main_ret_NOLOG;
    }
    char outname[kPglFnamesize];
    memcpy(outname, "plink2", 6);
    char* outname_end = nullptr;
    for (cur_flag_idx = 0; cur_flag_idx < flag_ct; ++cur_flag_idx) {
      int32_t memcmp_out_result = memcmp("out", &(flag_buf[cur_flag_idx * kMaxFlagBlen]), 4);
      if (!memcmp_out_result) {
	arg_idx = flag_map[cur_flag_idx];
	const uint32_t param_ct = param_count(argv, argc, arg_idx);
	if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	  fputs(g_logbuf, stderr);
	  fputs(errstr_append, stderr);
	  goto main_ret_INVALID_CMDLINE;
	}
	if (strlen(argv[arg_idx + 1]) > (kPglFnamesize - kMaxOutfnameExtBlen)) {
	  fflush(stdout);
	  fputs("Error: --out parameter too long.\n", stderr);
	  goto main_ret_OPEN_FAIL;
	}
	const uint32_t slen = strlen(argv[arg_idx + 1]);
	memcpy(outname, argv[arg_idx + 1], slen + 1);
	outname_end = &(outname[slen]);
      }
      if (memcmp_out_result <= 0) {
	break;
      }
    }
    if (init_logfile(0, outname, outname_end? outname_end : &(outname[6]))) {
      goto main_ret_OPEN_FAIL;
    }
    logstr(ver_str);
    logstr("\n");
    logprint("Options in effect:\n");
    for (cur_flag_idx = 0; cur_flag_idx < flag_ct; ++cur_flag_idx) {
      logprint("  --");
      logprint(&(flag_buf[cur_flag_idx * kMaxFlagBlen]));
      arg_idx = flag_map[cur_flag_idx] + 1;
      while ((arg_idx < (uint32_t)argc) && (!is_flag(argv[arg_idx]))) {
	logprint(" ");
	logprint(argv[arg_idx++]);
      }
      logprint("\n");
    }
    logprint("\n");

#ifdef _WIN32
    DWORD windows_dw = kTextbufSize;
    if (GetComputerName(g_textbuf, &windows_dw))
#else
    if (gethostname(g_textbuf, kTextbufSize) != -1)
#endif
    {
      logstr("Hostname: ");
      logstr(g_textbuf);
    }
    logstr("\nWorking directory: ");
    if (!getcwd(g_textbuf, kPglFnamesize)) {
      goto main_ret_READ_FAIL;
    }
    logstr(g_textbuf);
    logstr("\n");
    logprint("Start time: ");
    time_t rawtime;
    time(&rawtime);
    logprint(ctime(&rawtime));
    // ctime string always has a newline at the end
    logstr("\n");

    int32_t known_procs;
#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    pc.max_thread_ct = sysinfo.dwNumberOfProcessors;
    known_procs = pc.max_thread_ct;
#else
    known_procs = sysconf(_SC_NPROCESSORS_ONLN);
    pc.max_thread_ct = (known_procs == -1)? 1 : ((uint32_t)known_procs);
#endif
    // don't subtract 1 any more since, when max_thread_ct > 2, one of the
    // (virtual) cores will be dedicated to I/O and have lots of downtime.
    if (pc.max_thread_ct > kMaxThreads) {
      pc.max_thread_ct = kMaxThreads;
    }

    char pgenname[kPglFnamesize];
    char psamname[kPglFnamesize];
    char pvarname[kPglFnamesize];
    pgenname[0] = '\0';
    psamname[0] = '\0';
    pvarname[0] = '\0';
    init_pheno();
    cur_flag_idx = 0;
    pc.command_flags1 = kfCommand10;
    // uint64_t command_flags2 = 0;
    pc.misc_flags = kfMisc0;
    pc.pvar_psam_modifier = kfPvarPsam0;
    pc.exportf_modifier = kfExportf0;
    pc.sample_sort_flags = kfSort0;
    pc.grm_flags = kfGrm0;
    pc.pca_flags = kfPca0;
    pc.write_covar_flags = kfWriteCovar0;
    pc.pheno_transform_flags = kfPhenoTransform0;
    pc.fam_cols = kfFamCol13456;
    pc.exportf_id_paste = kfIdpaste0;
    pc.king_modifier = kfKing0;
    pc.king_cutoff = -1;
    pc.king_table_filter = -DBL_MAX;
    pc.allele_freq_modifier = kfAlleleFreq0;
    pc.missing_rpt_modifier = kfMissingRpt0;
    pc.geno_counts_modifier = kfGenoCounts0;
    pc.hardy_modifier = kfHardy0;
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
    pc.hard_call_thresh = 0xffffffffU;
    
    pc.dosage_erase_thresh = 0;
    pc.pfilter = 2.0; // make --pfilter 1 still filter out NAs
    pc.output_min_p = 0.0;
    pc.vif_thresh = 50.0;
    pc.mind_thresh = 1.0;
    pc.geno_thresh = 1.0;
    pc.hwe_thresh = 1.0;
    pc.mach_r2_min = 0.0;
    pc.mach_r2_max = 0.0;
    pc.min_maf = 0.0;
    pc.max_maf = 1.0;
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
    pc.exportf_id_delim = '\0';
    double import_dosage_certainty = 0.0;
    int32_t vcf_min_gq = -1;
    int32_t vcf_min_dp = -1;
    intptr_t malloc_size_mb = 0;
    load_params_t load_params = kfLoadParams0;
    xload_t xload = kfXload0;
    uint32_t rseed_ct = 0;
    make_plink2_t make_plink2_modifier = kfMake0;
    oxford_import_t oxford_import_flags = kfOxfordImport0;
    vcf_half_call_t vcf_half_call = kVcfHalfCallDefault;
    char range_delim = '-';
    char id_delim = '\0';
    char idspace_to = '\0';
    char input_missing_geno_char = '0';
    char output_missing_geno_char = '.';
    uint32_t aperm_present = 0;
    uint32_t notchr_present = 0;
    uint32_t permit_multiple_inclusion_filters = 0;
    uint32_t memory_require = 0;
    gendummy_info_t gendummy_info;
    init_gendummy(&gendummy_info);
    plink1_dosage_info_t plink1_dosage_info;
    init_plink1_dosage(&plink1_dosage_info);
    do {
      flagname_p = &(flag_buf[cur_flag_idx * kMaxFlagBlen]);
      if (!(*flagname_p)) {
	// preprocessed; not relevant now, but will need --d later
	continue;
      }
      char* flagname_p2 = &(flagname_p[1]);
      arg_idx = flag_map[cur_flag_idx];
      uint32_t param_ct = param_count(argv, argc, arg_idx);
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
	if (!memcmp(flagname_p2, "llow-extra-chr", 15)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    if (memcmp("0", cur_modif, 2)) {
	      sprintf(g_logbuf, "Error: Invalid --allow-extra-chr parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    chr_info.zero_extra_chrs = 1;
	  }
	  pc.misc_flags |= kfMiscAllowExtraChrs;
	} else if (!memcmp(flagname_p2, "utosome", 8)) {
	  pc.misc_flags |= kfMiscAutosomeOnly;
	  chr_info.is_include_stack = 1;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "utosome-par", 12)) {
	  if (pc.misc_flags & kfMiscAutosomeOnly) {
	    logerrprint("Error: --autosome-par cannot be used with --autosome.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  pc.misc_flags |= kfMiscAutosomePar;
	  chr_info.is_include_stack = 1;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "llow-no-samples", 16)) {
	  // er, lock these out until they're at least close to fully supported
	  logerrprint("Error: --allow-no-samples is not implemented yet.\n");
	  reterr = kPglRetNotYetSupported;
	  goto main_ret_1;
	  pc.misc_flags |= kfMiscAllowNoSamples;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "llow-no-vars", 13)) {
	  logerrprint("Error: --allow-no-vars is not implemented yet.\n");
	  reterr = kPglRetNotYetSupported;
	  goto main_ret_1;
	  pc.misc_flags |= kfMiscAllowNoVars;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "djust", 6)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 3)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 2) && (!memcmp(cur_modif, "gc", 2))) {
	      pc.adjust_info.flags |= kfAdjustGc;
	    } else if ((cur_modif_slen == 4) && (!memcmp(cur_modif, "log10", 5))) {
	      pc.adjust_info.flags |= kfAdjustLog10;
	    } else if ((cur_modif_slen > 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      if (pc.adjust_info.flags & kfAdjustColAll) {
		logerrprint("Error: Multiple --adjust cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      reterr = parse_col_descriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0unadj\0gc\0qq\0bonf\0holm\0sidakss\0sidaksd\0fdrbh\0fdrby\0", "adjust", kfAdjustColChrom, kfAdjustColDefault, 1, &pc.adjust_info.flags);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --adjust parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (!(pc.adjust_info.flags & kfAdjustColAll)) {
	    pc.adjust_info.flags |= kfAdjustColDefault;
	  }
	} else if (!memcmp(flagname_p2, "perm", 5)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 6)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  if (scan_posint_defcap(cur_modif, &pc.aperm.min)) {
	    sprintf(g_logbuf, "Error: Invalid --aperm min permutation count '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  ++pc.aperm.min;
	  if (param_ct > 1) {
	    cur_modif = argv[arg_idx + 2];
	    if (scan_posint_capped(cur_modif, kApermMax, &pc.aperm.max)) {
	      sprintf(g_logbuf, "Error: Invalid --aperm max permutation count '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (pc.aperm.min >= pc.aperm.max) {
	    logerrprint("Error: --aperm min permutation count must be smaller than max.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  aperm_present = 1;
	  if (param_ct > 2) {
	    cur_modif = argv[arg_idx + 3];
	    if (!scanadv_double(cur_modif, &pc.aperm.alpha)) {
	      sprintf(g_logbuf, "Error: Invalid --aperm alpha threshold '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (param_ct > 3) {
	      cur_modif = argv[arg_idx + 4];
	      if ((!scanadv_double(cur_modif, &pc.aperm.beta)) || (pc.aperm.beta <= 0.0)) {
		sprintf(g_logbuf, "Error: Invalid --aperm beta '%s'.\n", cur_modif);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      if (param_ct > 4) {
		cur_modif = argv[arg_idx + 5];
		if (!scanadv_double(cur_modif, &pc.aperm.init_interval)) {
		  sprintf(g_logbuf, "Error: Invalid --aperm initial pruning interval '%s'.\n", cur_modif);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
		if ((pc.aperm.init_interval < 1.0) || (pc.aperm.init_interval > 1000000.0)) {
		  sprintf(g_logbuf, "Error: Invalid --aperm initial pruning interval '%s'.\n", cur_modif);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
		if (param_ct == 6) {
		  cur_modif = argv[arg_idx + 6];
		  if (!scanadv_double(cur_modif, &pc.aperm.interval_slope) || (pc.aperm.interval_slope < 0.0) || (pc.aperm.interval_slope > 1.0)) {
		    sprintf(g_logbuf, "Error: Invalid --aperm pruning interval slope '%s'.\n", cur_modif);
		    goto main_ret_INVALID_CMDLINE_WWA;
		  }
		}
	      }
	    }
	  }
	} else if (!memcmp(flagname_p2, "utosome-num", 12)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  uint32_t autosome_ct;
	  if (scan_posint_capped(cur_modif, kMaxChrTextnum, &autosome_ct)) {
	    sprintf(g_logbuf, "Error: Invalid --autosome-num parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  // see plink2_common finalize_chrset()
	  chr_info.chrset_source = kChrsetSourceCmdline;
	  chr_info.autosome_ct = autosome_ct;
	  // assumes first code is X
	  chr_info.xymt_codes[0] = autosome_ct + 1;
	  for (uint32_t xymt_idx = 1; xymt_idx < kChrOffsetCt; ++xymt_idx) {
	    // bugfix: this needs to be -2, not -1, for get_chr_code() to work
	    // properly
	    chr_info.xymt_codes[xymt_idx] = -2;
	  }
	  chr_info.haploid_mask[0] = 0;
	  set_bit(autosome_ct + 1, chr_info.haploid_mask);
	} else if (!memcmp(flagname_p2, "llow-no-sex", 12)) {
	  logprint("Note: --allow-no-sex no longer has any effect.  (Missing-sex samples are\nautomatically excluded from association analysis when sex is a covariate, and\ntreated normally otherwise.)\n");
	  goto main_param_zero;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'b':
	if (!memcmp(flagname_p2, "file", 5)) {
	  if (xload) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t fname_modif_idx = 1;
	  if (param_ct == 2) {
	    if (check_extra_param(&(argv[arg_idx]), "vzs", &fname_modif_idx)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	  const char* fname_prefix = argv[arg_idx + fname_modif_idx];
	  const uint32_t slen = strlen(fname_prefix);
	  if (slen > (kPglFnamesize - 9)) {
	    // could use kPglFnamesize - 2 - 3 * param_ct, but that's pointless
	    logerrprint("Error: --bfile parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(memcpya(pgenname, fname_prefix, slen), ".bed");
	  strcpy(memcpya(psamname, fname_prefix, slen), ".fam");
	  char* bimname_end = memcpya(pvarname, fname_prefix, slen);
	  bimname_end = strcpya0(bimname_end, ".bim");
	  if (param_ct == 2) {
	    strcpy(bimname_end, ".zst");
	  }
	  load_params |= kfLoadParamsPfileAll;
	} else if (!memcmp(flagname_p2, "pfile", 6)) {
	  if (load_params || xload) {
	    // currently only possible with --bcf, --bfile
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t fname_modif_idx = 1;
	  if (param_ct == 2) {
	    if (check_extra_param(&(argv[arg_idx]), "vzs", &fname_modif_idx)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	  const char* fname_prefix = argv[arg_idx + fname_modif_idx];
	  const uint32_t slen = strlen(fname_prefix);
	  if (slen > (kPglFnamesize - 9)) {
	    logerrprint("Error: --bpfile parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(memcpya(pgenname, fname_prefix, slen), ".pgen");
	  strcpy(memcpya(psamname, fname_prefix, slen), ".fam");
	  char* bimname_end = memcpya(pvarname, fname_prefix, slen);
	  bimname_end = strcpya0(bimname_end, ".bim");
	  if (param_ct == 2) {
	    strcpy(bimname_end, ".zst");
	  }
	  load_params |= kfLoadParamsPfileAll;
	} else if (!memcmp(flagname_p2, "iallelic-only", 14)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "strict")) {
	      pc.misc_flags |= kfMiscBiallelicOnlyStrict;
	    } else if (!strcmp(cur_modif, "list")) {
	      pc.misc_flags |= kfMiscBiallelicOnlyList;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --biallelic-only parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  pc.misc_flags |= kfMiscBiallelicOnly;
	  logerrprint("Error: --biallelic-only is not implemented yet.\n");
	  reterr = kPglRetNotYetSupported;
	  goto main_ret_1;
	} else if (!memcmp(flagname_p2, "cf", 3)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct == 2) {
	    const char* cur_modif = argv[arg_idx + 2];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen < 8) || memcmp(cur_modif, "dosage=", 7)) {
	      sprintf(g_logbuf, "Error: Invalid --bcf parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    reterr = cmdline_alloc_string(&(cur_modif[7]), argv[arg_idx], 4095, &vcf_dosage_import_field);
	    if (reterr) {
	      goto main_ret_1;
	    }
	    if (!is_alphanumeric(vcf_dosage_import_field)) {
	      logerrprint("Error: --bcf dosage= parameter is not alphanumeric.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (!strcmp(vcf_dosage_import_field, "GT")) {
	      logerrprint("Error: --bcf dosage= parameter cannot be 'GT'.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  const uint32_t slen = strlen(cur_modif);
	  if (slen > kPglFnamesize - 1) {
	    logerrprint("Error: --bcf filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pgenname, cur_modif, slen + 1);
	  xload = kfXloadBcf;
	} else if (!memcmp(flagname_p2, "gen", 4)) {
	  if (load_params || xload) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 3)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "snpid-chr")) {
	      oxford_import_flags |= kfOxfordImportBgenSnpIdChr;
	    } else if (!strcmp(cur_modif, "ref-first")) {
	      oxford_import_flags |= kfOxfordImportRefFirst;
	    } else if (!strcmp(cur_modif, "ref-second")) {
	      oxford_import_flags |= kfOxfordImportRefSecond;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --bgen parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  const char* cur_fname = argv[arg_idx + 1];
	  const uint32_t slen = strlen(cur_fname);
	  if (slen > kPglFnamesize - 1) {
	    logerrprint("Error: --bgen filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pgenname, cur_fname, slen + 1);
	  xload = kfXloadOxBgen;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'c':
	if (!memcmp(flagname_p2, "hr", 3)) {
	  if (pc.misc_flags & (kfMiscAutosomePar | kfMiscAutosomeOnly)) {
	    logerrprint("Error: --chr cannot be used with --autosome{-par}.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = parse_chr_ranges(flagname_p, errstr_append, param_ct, (pc.misc_flags / kfMiscAllowExtraChrs) & 1, 0, '-', &(argv[arg_idx]), &chr_info, chr_info.chr_mask);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  chr_info.is_include_stack = 1;
	} else if (!memcmp(flagname_p2, "ovar", 5)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_fname(argv[arg_idx + 1], flagname_p, 0, &pc.covar_fname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "ovar-name", 10)) {
	  // can now be used without --covar
	  reterr = parse_name_ranges(&(argv[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.covar_range_list);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "onst-fid", 9)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(param_ct? argv[arg_idx + 1] : "0", argv[arg_idx], kMaxIdSlen, &const_fid);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "i", 2)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (!scanadv_double(argv[arg_idx + 1], &pc.ci_size)) {
	    sprintf(g_logbuf, "Error: Invalid --ci parameter '%s'.\n", argv[arg_idx + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if ((pc.ci_size < 0.01) || (pc.ci_size >= 1.0)) {
	    logerrprint("Error: --ci confidence interval size must be in [0.01, 1).\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	} else if (!memcmp(flagname_p2, "ondition", 9)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t fname_modif_idx = 1;
	  if (param_ct == 2) {
	    if (!strcmp("dominant", argv[arg_idx + 2])) {
	      pc.glm_info.flags |= kfGlmConditionDominant;
	    } else if (!strcmp("recessive", argv[arg_idx + 2])) {
	      pc.glm_info.flags |= kfGlmConditionRecessive;
	    } else {
	      fname_modif_idx = 2;
	      if (!strcmp("dominant", argv[arg_idx + 1])) {
		pc.glm_info.flags |= kfGlmConditionDominant;
	      } else if (!strcmp("recessive", argv[arg_idx + 1])) {
		pc.glm_info.flags |= kfGlmConditionRecessive;
	      } else {
		logerrprint("Error: Invalid --condition parameter sequence.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    }
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + fname_modif_idx], argv[arg_idx], kMaxIdSlen, &pc.glm_info.condition_varname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "ondition-list", 14)) {
	  if (pc.glm_info.condition_varname) {
	    logerrprint("Error: --condition-list cannot be used with --condition.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t fname_modif_idx = 1;
	  if (param_ct == 2) {
	    if (!strcmp("dominant", argv[arg_idx + 2])) {
	      pc.glm_info.flags |= kfGlmConditionDominant;
	    } else if (!strcmp("recessive", argv[arg_idx + 2])) {
	      pc.glm_info.flags |= kfGlmConditionRecessive;
	    } else {
	      fname_modif_idx = 2;
	      if (!strcmp("dominant", argv[arg_idx + 1])) {
		pc.glm_info.flags |= kfGlmConditionDominant;
	      } else if (!strcmp("recessive", argv[arg_idx + 1])) {
		pc.glm_info.flags |= kfGlmConditionRecessive;
	      } else {
		logerrprint("Error: Invalid --condition-list parameter sequence.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    }
	  }
	  reterr = alloc_fname(argv[arg_idx + fname_modif_idx], flagname_p, 0, &pc.glm_info.condition_list_fname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "ow", 3)) {
	  if (chr_info.chrset_source) {
	    logerrprint("Error: Conflicting chromosome-set flags.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  chr_info.chrset_source = kChrsetSourceCmdline;
	  // initialize here instead of finalize_chrset(), to simplify
	  // read_chrset_header_line()
	  chr_info.autosome_ct = 29;
	  chr_info.xymt_codes[0] = 30;
	  chr_info.xymt_codes[1] = 31;
	  chr_info.xymt_codes[2] = -2;
	  chr_info.xymt_codes[3] = 33;
	  chr_info.xymt_codes[4] = -2;
	  chr_info.xymt_codes[5] = -2;
	  chr_info.haploid_mask[0] = 0xc0000000U;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "hr-set", 7)) {
	  if (chr_info.chrset_source) {
	    logerrprint("Error: Conflicting chromosome-set flags.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 5)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  int32_t signed_autosome_ct;
	  if (scan_int_abs_bounded(cur_modif, kMaxChrTextnum, &signed_autosome_ct) || (!signed_autosome_ct)) {
	    sprintf(g_logbuf, "Error: Invalid --chr-set parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  // see plink2_common finalize_chrset()
	  chr_info.chrset_source = kChrsetSourceCmdline;
	  chr_info.haploid_mask[0] = 0;
	  if (signed_autosome_ct < 0) {
	    // haploid
	    if (param_ct > 1) {
	      logerrprint("Error: --chr-set does not accept multiple parameters in haploid mode.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    const uint32_t autosome_ct = -signed_autosome_ct;
	    chr_info.autosome_ct = autosome_ct;
	    for (uint32_t xymt_idx = 0; xymt_idx < kChrOffsetCt; ++xymt_idx) {
	      chr_info.xymt_codes[xymt_idx] = -2;
	    }
	    fill_all_bits(autosome_ct + 1, chr_info.haploid_mask);
	  } else {
	    const uint32_t autosome_ct = signed_autosome_ct;
	    chr_info.autosome_ct = autosome_ct;
	    // assumes first four codes are x, y, xy, mt
	    for (uint32_t xymt_idx = 0; xymt_idx < 4; ++xymt_idx) {
	      chr_info.xymt_codes[xymt_idx] = autosome_ct + 1 + xymt_idx;
	    }
	    for (uint32_t xymt_idx = 4; xymt_idx < kChrOffsetCt; ++xymt_idx) {
	      chr_info.xymt_codes[xymt_idx] = -2;
	    }
	    set_bit(autosome_ct + 1, chr_info.haploid_mask);
	    set_bit(autosome_ct + 2, chr_info.haploid_mask);
	    for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
	      cur_modif = argv[arg_idx + param_idx];
	      if (!strcmp(cur_modif, "no-x")) {
		chr_info.xymt_codes[0] = -2;
		clear_bit(autosome_ct + 1, chr_info.haploid_mask);
	      } else if (!strcmp(cur_modif, "no-y")) {
		chr_info.xymt_codes[1] = -2;
		clear_bit(autosome_ct + 2, chr_info.haploid_mask);
	      } else if (!strcmp(cur_modif, "no-xy")) {
		chr_info.xymt_codes[2] = -2;
	      } else if (!strcmp(cur_modif, "no-mt")) {
		chr_info.xymt_codes[3] = -2;
	      } else {
		sprintf(g_logbuf, "Error: Invalid --chr-set parameter '%s'.\n", cur_modif);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	    }
	  }
	} else if (!memcmp(flagname_p2, "hr-override", 12)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    if (!strcmp(cur_modif, "file")) {
	      pc.misc_flags |= kfMiscChrOverrideFile;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --chr-override parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  } else {
	    pc.misc_flags |= kfMiscChrOverrideCmdline;
	  }
	} else if (!memcmp(flagname_p2, "ovar-quantile-normalize", 24)) {
	  if (param_ct) {
	    reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, 0x7fffffff, &pc.covar_quantnorm_flattened);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.pheno_transform_flags |= kfPhenoTransformQuantnormCovar;
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "ovar-variance-standardize", 26)) {
	  if (param_ct) {
	    reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, 0x7fffffff, &pc.vstd_flattened);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.pheno_transform_flags |= kfPhenoTransformVstdCovar;
	  pc.filter_flags |= kfFilterPsamReq;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'd':
	if (!memcmp(flagname_p2, "ouble-id", 9)) {
	  if (const_fid) {
	    logerrprint("Error: --double-id cannot be used with --const-fid.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  pc.misc_flags |= kfMiscDoubleId;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "ebug", 5)) {
	  g_debug_on = 1;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "ata", 4)) {
	  if (load_params || (xload & (~kfXloadOxBgen))) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 3)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t is_gzs = 0;
	  for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "ref-first")) {
	      oxford_import_flags |= kfOxfordImportRefFirst;
	    } else if (!strcmp(cur_modif, "ref-second")) {
	      oxford_import_flags |= kfOxfordImportRefSecond;
	    } else if (!strcmp(cur_modif, "gzs")) {
	      if (xload & kfXloadOxBgen) {
		// may as well permit e.g. --data ref-first + --bgen
		logerrprint("Error: --data 'gzs' modifier cannot be used with .bgen input.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      is_gzs = 1;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --data parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  const char* fname_prefix = argv[arg_idx + 1];
	  const uint32_t slen = strlen(fname_prefix);
	  if (slen > (kPglFnamesize - 9)) {
	    logerrprint("Error: --data parameter too long.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (!(xload & kfXloadOxBgen)) {
	    // allow --bgen to override this
	    char* genname_end = memcpya(pgenname, fname_prefix, slen);
	    genname_end = strcpya0(genname_end, ".gen");
	    if (is_gzs) {
	      strcpy(genname_end, ".zst");
	    }
	    xload |= kfXloadOxGen;
	  }
	  strcpy(memcpya(psamname, fname_prefix, slen), ".sample");
	  xload |= kfXloadOxSample;
	} else if (!memcmp(flagname_p2, "osage-erase-threshold", 22)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  double dosage_erase_frac;
	  if ((!scanadv_double(cur_modif, &dosage_erase_frac)) || (dosage_erase_frac < 0.0) || (dosage_erase_frac >= (0.5 - kSmallEpsilon))) {
	    sprintf(g_logbuf, "Error: Invalid --dosage-erase-threshold parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  pc.dosage_erase_thresh = (int32_t)(dosage_erase_frac * ((1 + kSmallEpsilon) * kDosageMid));
	} else if (!memcmp(flagname_p2, "ummy", 5)) {
	  if (load_params || xload) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 2, 8)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  // todo: support --allow-no-samples/--allow-no-vars
	  if (scan_posint_defcap(argv[arg_idx + 1], &gendummy_info.sample_ct)) {
	    logerrprint("Error: Invalid --dummy sample count.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (scan_posint_defcap(argv[arg_idx + 2], &gendummy_info.variant_ct)) {
	    logerrprint("Error: Invalid --dummy SNP count.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  uint32_t extra_numeric_param_ct = 0;
	  for (uint32_t param_idx = 3; param_idx <= param_ct; ++param_idx) {
	    char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 4) && match_upper_counted(cur_modif, "ACGT", 4)) {
	      gendummy_info.flags |= kfGenDummyAcgt;
	    } else if ((cur_modif_slen == 4) && (!memcmp(cur_modif, "1234", 4))) {
	      gendummy_info.flags |= kfGenDummy1234;
	    } else if ((cur_modif_slen == 2) && (!memcmp(cur_modif, "12", 2))) {
	      gendummy_info.flags |= kfGenDummy12;
	    } else if ((cur_modif_slen > 9) && (!memcmp(cur_modif, "pheno-ct=", 9))) {
	      const char* pheno_ct_start = &(cur_modif[9]);
	      if (scan_uint_capped(pheno_ct_start, kMaxPhenoCt, &gendummy_info.pheno_ct)) {
		sprintf(g_logbuf, "Error: Invalid --dummy pheno-ct= parameter '%s'.\n", pheno_ct_start);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	    } else if ((cur_modif_slen == 12) && (!memcmp(cur_modif, "scalar-pheno", 12))) {
	      gendummy_info.flags |= kfGenDummyScalarPheno;
	    } else if ((cur_modif_slen > 12) && (!memcmp(cur_modif, "dosage-freq=", 12))) {
	      char* dosage_freq_start = &(cur_modif[12]);
	      double dxx;
	      if ((!scanadv_double(dosage_freq_start, &dxx)) || (dxx < 0.0) || (dxx > 1.0)) {
		sprintf(g_logbuf, "Error: Invalid --dummy dosage-freq= parameter '%s'.\n", dosage_freq_start);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      gendummy_info.dosage_freq = dxx;
	    } else {
	      double dxx;
	      if ((extra_numeric_param_ct == 2) || (!scanadv_double(cur_modif, &dxx)) || (dxx < 0.0) || (dxx > 1.0)) {
		sprintf(g_logbuf, "Error: Invalid --dummy parameter '%s'.\n", cur_modif);
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
	    logerrprint("Error: --dummy 'acgt', '1234', and '12' modifiers are mutually exclusive.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  xload |= kfXloadGenDummy;
	} else if (!memcmp(flagname_p2, "ummy-coding", 12)) {
	  logerrprint("Error: --dummy-coding is retired.  Use --split-cat-pheno instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (!memcmp(flagname_p2, "og", 3)) {
	  if (chr_info.chrset_source) {
	    logerrprint("Error: Conflicting chromosome-set flags.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  chr_info.chrset_source = kChrsetSourceCmdline;
	  chr_info.autosome_ct = 38;
	  chr_info.xymt_codes[0] = 39;
	  chr_info.xymt_codes[1] = 40;
	  chr_info.xymt_codes[2] = 41;
	  chr_info.xymt_codes[3] = 42;
	  chr_info.xymt_codes[4] = -2;
	  chr_info.xymt_codes[5] = -2;
#ifdef __LP64__
	  chr_info.haploid_mask[0] = 0x18000000000LLU;
#else
	  chr_info.haploid_mask[0] = 0;
	  chr_info.haploid_mask[1] = 0x180;
#endif
	  goto main_param_zero;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'e':
	if (!memcmp(flagname_p2, "xtract", 7)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const uint32_t is_range = !strcmp(argv[arg_idx + 1], "range");
	  if (is_range) {
	    if (param_ct == 1) {
	      logerrprint("Error: '--extract range' requires at least one filename.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    pc.misc_flags |= kfMiscExtractRange;
	    pc.filter_flags |= kfFilterNoSplitChr;
	  }
	  reterr = alloc_and_flatten(&(argv[arg_idx + 1 + is_range]), param_ct - is_range, kPglFnamesize, &pc.extract_fnames);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPvarReq;
	} else if (!memcmp(flagname_p2, "xclude", 7)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const uint32_t is_range = !strcmp(argv[arg_idx + 1], "range");
	  if (is_range) {
	    if (param_ct == 1) {
	      logerrprint("Error: '--exclude range' requires at least one filename.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    pc.misc_flags |= kfMiscExcludeRange;
	    pc.filter_flags |= kfFilterNoSplitChr;
	  }
	  reterr = alloc_and_flatten(&(argv[arg_idx + 1 + is_range]), param_ct - is_range, kPglFnamesize, &pc.exclude_fnames);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPvarReq;
	} else if ((!memcmp(flagname_p2, "xport", 6)) || (!memcmp(flagname_p2, "xport ped", 10))) {
	  // todo: determine actual limit
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 50)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t format_param_idxs = 0;
	  if (!flagname_p2[5]) {
	    get_exportf_targets(&(argv[arg_idx]), param_ct, &pc.exportf_modifier, &pc.exportf_id_paste, &format_param_idxs);
	    if (!format_param_idxs) {
	      logerrprint("Error: --export requires at least one output format.  (Did you forget 'ped' or\n'vcf'?)\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    logprintb();
	  } else {
	    pc.exportf_modifier = kfExportfPed;
	  }
	  // can't have e.g. bgen-1.1 and bgen-1.2 simultaneously, since they
	  // have the same extension and different content.
	  const uint64_t bgen_flags = (uint64_t)(pc.exportf_modifier & (kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13));
	  if (bgen_flags & (bgen_flags - 1)) {
	    logerrprint("Error: Multiple --export bgen versions.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if ((pc.exportf_modifier & (kfExportfHaps | kfExportfHapsLegend)) == (kfExportfHaps | kfExportfHapsLegend)) {
	    logerrprint("Error: 'haps' and 'hapslegend' formats cannot be exported simultaneously.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    // could use next_unset()...
	    if ((format_param_idxs >> param_idx) & 1) {
	      continue;
	    }
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if (cur_modif_slen > 9) {
	      if (!memcmp(cur_modif, "id-paste=", 9)) {
		if (!(pc.exportf_modifier & (kfExportfVcf | kfExportfBgen12 | kfExportfBgen13))) {
		  // todo: bcf
		  logerrprint("Error: The 'id-paste' modifier only applies to --export's vcf, bgen-1.2, and\nbgen-1.3 output formats.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
		if (pc.exportf_id_paste) {
		  logerrprint("Error: Multiple --export id-paste= modifiers.\n");
		  goto main_ret_INVALID_CMDLINE;
		}
		reterr = parse_col_descriptor(&(cur_modif[9]), "fid\0iid\0maybesid\0sid\0", "export", kfIdpasteFid, kfIdpasteDefault, 1, &pc.exportf_id_paste);
		if (reterr) {
		  goto main_ret_1;
		}
	      } else if (!memcmp(cur_modif, "id-delim=", 9)) {
		if (!(pc.exportf_modifier & (kfExportfVcf | kfExportfBgen12 | kfExportfBgen13))) {
		  logerrprint("Error: The 'id-delim' modifier only applies to --export's vcf, bgen-1.2, and\nbgen-1.3 output formats.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
		if (pc.exportf_id_delim) {
		  logerrprint("Error: Multiple --export id-delim= modifiers.\n");
		  goto main_ret_INVALID_CMDLINE;
		}
		pc.exportf_id_delim = extract_char_param(&(cur_modif[9]));
		if (!pc.exportf_id_delim) {
		  logerrprint("Error: --export id-delim= value must be a single character.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
		if ((((unsigned char)pc.exportf_id_delim) < ' ') || (pc.exportf_id_delim == '0')) {
		  logerrprint("Error: --export id-delim= value cannot be tab, newline, '0', or a nonprinting\ncharacter.\n");
		  goto main_ret_INVALID_CMDLINE;
		}
	      } else if (!memcmp(cur_modif, "vcf-dosage=", 11)) {
		if (!(pc.exportf_modifier & kfExportfVcf)) {
		  logerrprint("Error: The 'vcf-dosage' modifier only applies to --export's vcf output format.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
		if (pc.exportf_modifier & (kfExportfVcfDosageGp | kfExportfVcfDosageDs)) {
		  logerrprint("Error: Multiple --export vcf-dosage= modifiers.\n");
		  goto main_ret_INVALID_CMDLINE;
		}
		const char* vcf_dosage_start = &(cur_modif[11]);
		if (!strcmp(vcf_dosage_start, "GP")) {
		  pc.exportf_modifier |= kfExportfVcfDosageGp;
		} else if (!strcmp(vcf_dosage_start, "DS")) {
		  pc.exportf_modifier |= kfExportfVcfDosageDs;
		} else {
		  sprintf(g_logbuf, "Error: Invalid --export vcf-dosage= parameter '%s'.\n", vcf_dosage_start);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
	      } else if (!memcmp(cur_modif, "bits=", 5)) {
		if (!(pc.exportf_modifier & (kfExportfBgen12 | kfExportfBgen13))) {
		  logerrprint("Error: The 'bits' modifier only applies to --export's bgen-1.2 and bgen-1.3\noutput formats.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
		if (pc.exportf_bits) {
		  logerrprint("Error: Multiple --export bits= modifiers.\n");
		  goto main_ret_INVALID_CMDLINE;
		}
		const char* bits_start = &(cur_modif[5]);
		if (scan_posint_capped(bits_start, 32, &pc.exportf_bits)) {
		  sprintf(g_logbuf, "Error: Invalid --export bits= parameter '%s'.\n", bits_start);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
	      } else if ((cur_modif_slen == 11) && (!memcmp(cur_modif, "include-alt", 11))) {
		if (!(pc.exportf_modifier & (kfExportfA | kfExportfAD))) {
		  logerrprint("Error: The 'include-alt' modifier only applies to --export's A and AD output\nformats.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
		pc.exportf_modifier |= kfExportfIncludeAlt;
	      } else if ((cur_modif_slen == 14) && (!memcmp(cur_modif, "omit-nonmale-y", 14))) {
		if (!(pc.exportf_modifier & (kfExportfList | kfExportfRlist))) {
		  logerrprint("Error: The 'omit-nonmale-y' modifier only applies to --export's list and rlist\noutput formats.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
		pc.exportf_modifier |= kfExportfOmitNonmaleY;
	      }
	    } else if ((cur_modif_slen == 2) && ((!memcmp(cur_modif, "01", 2)) || (!memcmp(cur_modif, "12", 2)))) {
	      if (pc.exportf_modifier & (kfExportfA | kfExportfAD)) {
		sprintf(g_logbuf, "Error: The '%s' modifier does not apply to --export's A and AD output formats.\n", cur_modif);
		goto main_ret_INVALID_CMDLINE_2A;
	      }
	      if (pc.exportf_modifier & kfExportfVcf) {
		logerrprint("Error: '01'/'12' cannot be used with --export's vcf output format.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      if (cur_modif[0] == '0') {
		if (pc.exportf_modifier & kfExportf12) {
		  logerrprint("Error: --export '01' and '12' cannot be used together.\n");
		  goto main_ret_INVALID_CMDLINE;
		}
		pc.exportf_modifier |= kfExportf01;
	      } else {
		if (pc.exportf_modifier & kfExportf01) {
		  logerrprint("Error: --export '01' and '12' cannot be used together.\n");
		  goto main_ret_INVALID_CMDLINE;
		}
		pc.exportf_modifier |= kfExportf12;
	      }
	    } else if ((cur_modif_slen == 3) && (!memcmp(cur_modif, "bgz", 3))) {
	      if (!(pc.exportf_modifier & kfExportfVcf)) {
		logerrprint("Error: The 'bgz' modifier only applies to --export's vcf output format.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.exportf_modifier |= kfExportfBgz;
	    } else if ((cur_modif_slen == 6) && (!memcmp(cur_modif, "spaces", 6))) {
	      pc.exportf_modifier |= kfExportfSpaces;
	    } else if ((cur_modif_slen == 9) && (!memcmp(cur_modif, "ref-first", 9))) {
	      pc.exportf_modifier |= kfExportfRefFirst;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --export parameter '%s'.%s\n", cur_modif, ((param_idx == param_ct) && (!outname_end))? " (Did you forget '--out'?)" : "");
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (pc.exportf_modifier & (kfExportfVcf | kfExportfBgen12 | kfExportfBgen13)) {
	    if (!pc.exportf_id_paste) {
	      pc.exportf_id_paste = kfIdpasteDefault;
	    }
	  }
	  pc.command_flags1 |= kfCommand1Exportf;
	} else if (!memcmp(flagname_p2, "xclude-snp", 11)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + 1], argv[arg_idx], kMaxIdSlen, &pc.varid_exclude_snp);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPvarReq;
	} else if (!memcmp(flagname_p2, "xclude-snps", 12)) {
	  reterr = parse_name_ranges(&(argv[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.exclude_snps_range_list);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPvarReq;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;

      case 'f':
	if (!memcmp(flagname_p2, "req", 4)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 5)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t bins_only = 0;
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 2) && (!memcmp(cur_modif, "zs", 2))) {
	      pc.allele_freq_modifier |= kfAlleleFreqZs;
	    } else if ((cur_modif_slen == 6) && (!memcmp(cur_modif, "counts", 6))) {
	      pc.allele_freq_modifier |= kfAlleleFreqCounts;
	    } else if ((cur_modif_slen == 12) && (!memcmp(cur_modif, "case-control", 12))) {
	      logerrprint("Error: --freq 'case-control' modifier has been retired.  Use\n--keep-if/--remove-if in conjunction with Unix text-processing utilities\ninstead.\n");
	    } else if ((cur_modif_slen > 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      if (pc.allele_freq_modifier & kfAlleleFreqColAll) {
		logerrprint("Error: Multiple --freq cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      reterr = parse_col_descriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0reffreq\0alt1freq\0altfreq\0freq\0eq\0eqz\0alteq\0alteqz\0numeq\0altnumeq\0machr2\0nobs\0", "freq", kfAlleleFreqColChrom, kfAlleleFreqColDefault, 1, &pc.allele_freq_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	      const uint32_t mutually_exclusive_cols = pc.allele_freq_modifier & kfAlleleFreqColMutex;
	      if (mutually_exclusive_cols & (mutually_exclusive_cols - 1)) {
		logerrprint("Error: --freq's altfreq, freq, eq, eqz, alteq, alteqz, numeq, and altnumeq\ncolumns are mutually exclusive.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else if ((cur_modif_slen == 9) && (!memcmp(cur_modif, "bins-only", 9))) {
	      bins_only = 1;
	    } else if (((cur_modif_slen > 8) && (!memcmp(cur_modif, "refbins=", 8))) || ((cur_modif_slen > 13) && (!memcmp(cur_modif, "refbins-file=", 13))) || ((cur_modif_slen > 9) && (!memcmp(cur_modif, "alt1bins=", 9))) || ((cur_modif_slen > 14) && (!memcmp(cur_modif, "alt1bins-file=", 14)))) {
	      const uint32_t is_alt1 = (cur_modif[0] == 'a');
	      char** binstr_ptr = is_alt1? (&pc.freq_alt1_binstr) : (&pc.freq_ref_binstr);
	      if (*binstr_ptr) {
		LOGERRPRINTF("Error: Multiple --freq %sbins{-file}= modifiers.\n", is_alt1? "alt1" : "ref");
		goto main_ret_INVALID_CMDLINE;
	      }
	      if (cur_modif[7 + is_alt1] == '=') {
		reterr = cmdline_alloc_string(&(cur_modif[8 + is_alt1]), is_alt1? "--freq alt1bins=" : "--freq refbins=", 0x7fffffff, binstr_ptr);
	      } else {
		pc.allele_freq_modifier |= is_alt1? kfAlleleFreqBinsAlt1Fname : kfAlleleFreqBinsRefFname;
		reterr = alloc_fname(&(cur_modif[13 + is_alt1]), is_alt1? "freq alt1bins-file=" : "freq refbins-file=", 0, binstr_ptr);
	      }
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --freq parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (bins_only) {
	    if ((!pc.freq_ref_binstr) && (!pc.freq_alt1_binstr)) {
	      logerrprint("Error: --freq 'bins-only' must be used with 'refbins{-file}=' and/or\n'alt1bins{-file}='.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (pc.allele_freq_modifier & (kfAlleleFreqZs | kfAlleleFreqColAll)) {
	      logerrprint("Error: --freq 'bins-only' cannot be used with 'zs' or 'cols=' (which only\naffect the main report).\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    pc.allele_freq_modifier |= kfAlleleFreqBinsOnly;
	  }
	  if (!(pc.allele_freq_modifier & kfAlleleFreqColAll)) {
	    pc.allele_freq_modifier |= kfAlleleFreqColDefault;
	  }
	  pc.command_flags1 |= kfCommand1AlleleFreq;
	} else if (!memcmp(flagname_p2, "rom", 4)) {
	  if (chr_info.is_include_stack) {
	    logerrprint("Error: --from/--to cannot be used with --autosome{-par} or --{not-}chr.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + 1], argv[arg_idx], kMaxIdSlen, &pc.varid_from);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
	} else if ((!memcmp(flagname_p2, "rom-bp", 7)) || (!memcmp(flagname_p2, "rom-kb", 7)) || (!memcmp(flagname_p2, "rom-mb", 7))) {
	  if (!cmdline_single_chr(&chr_info, pc.misc_flags)) {
	    logerrprint("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb must be used with --chr, and only\none chromosome.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (pc.from_bp != -1) {
	    logerrprint("Error: Multiple --from-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  // permit negative numbers, to simplify shell script windowing logic
	  char* cur_modif = argv[arg_idx + 1];
	  double dxx;
	  if (!scanadv_double(cur_modif, &dxx)) {
	    sprintf(g_logbuf, "Error: Invalid --from-bp/-kb/-mb parameter '%s'.\n", cur_modif);
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
	      LOGERRPRINTF("Error: --from-bp/-kb/-mb parameter '%s' too large.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    pc.from_bp = 1 + (int32_t)(dxx * (1 - kSmallEpsilon));
	  }
	  pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
	} else if (!memcmp(flagname_p2, "orce-intersect", 15)) {
	  permit_multiple_inclusion_filters = 1;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "amily", 6)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if (is_reserved_pheno_name(cur_modif, cur_modif_slen)) {
	      sprintf(g_logbuf, "Error: '%s' cannot be used as a categorical phenotype name.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    reterr = cmdline_alloc_string(cur_modif, argv[arg_idx], kMaxIdSlen, &pc.catpheno_name);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.misc_flags |= kfMiscCatPhenoFamily;
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "amily-missing-catname", 22)) {
	  if (!(pc.misc_flags & kfMiscCatPhenoFamily)) {
	    logerrprint("Error: --family-missing-catname must be used with --family.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + 1], argv[arg_idx], kMaxIdSlen, &pc.family_missing_catname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if ((!memcmp(flagname_p2, "ilter-cases", 12)) || (!memcmp(flagname_p2, "ilter-controls", 15))) {
	  logerrprint("Error: --filter-cases and --filter-controls have been retired.  Use\n--keep-if/--remove-if instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if ((!memcmp(flagname_p2, "rqx", 4)) || (!memcmp(flagname_p2, "reqx", 5))) {
	  logerrprint("Error: --freqx has been retired.  Use --geno-counts instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'g':
	if (!memcmp(flagname_p2, "eno", 4)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t geno_thresh_present = 0;
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "dosage")) {
	      pc.misc_flags |= kfMiscGenoDosage;
	    } else if (!strcmp(cur_modif, "hh-missing")) {
	      pc.misc_flags |= kfMiscGenoHhMissing;
	    } else if (geno_thresh_present) {
	      logerrprint("Error: Invalid --geno parameter sequence.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (!scanadv_double(cur_modif, &pc.geno_thresh)) {
	      sprintf(g_logbuf, "Error: Invalid --geno parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    } else if ((pc.geno_thresh < 0.0) || (pc.geno_thresh > 1.0)) {
	      sprintf(g_logbuf, "Error: Invalid --geno parameter '%s' (must be in [0, 1]).\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    } else {
	      geno_thresh_present = 1;
	    }
	  }
	  if (!geno_thresh_present) {
	    pc.geno_thresh = 0.1;
	  }
	  if (pc.geno_thresh < 1.0) {
	    pc.filter_flags |= kfFilterAllReq | kfFilterNoSplitChr;
	  }
	} else if (!memcmp(flagname_p2, "eno-counts", 11)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 2) && (!memcmp(cur_modif, "zs", 2))) {
	      pc.geno_counts_modifier |= kfGenoCountsZs;
	    } else if ((cur_modif_slen > 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      if (pc.geno_counts_modifier & kfGenoCountsColAll) {
		logerrprint("Error: Multiple --geno-counts cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }	      
	      reterr = parse_col_descriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0homref\0refalt1\0refalt\0homalt1\0altxy\0xy\0hapref\0hapalt1\0hapalt\0hap\0numeq\0missing\0nobs\0", "geno-counts", kfGenoCountsColChrom, kfGenoCountsColDefault, 1, &pc.geno_counts_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	      if ((pc.geno_counts_modifier & kfGenoCountsColPairex) == kfGenoCountsColPairex) {
		logerrprint("Error: --geno-counts's hapaltx and hapx columns are mutually exclusive.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      const uint32_t mutually_exclusive_cols = pc.geno_counts_modifier & kfGenoCountsColMutex;
	      if (mutually_exclusive_cols & (mutually_exclusive_cols - 1)) {
		logerrprint("Error: --geno-counts's altxy, xy, and numeq columns are mutually exclusive.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --geno-counts parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (!(pc.geno_counts_modifier & kfGenoCountsColAll)) {
	    pc.geno_counts_modifier |= kfGenoCountsColDefault;
	  }
	  pc.command_flags1 |= kfCommand1GenoCounts;
	} else if (!memcmp(flagname_p2, "lm", 3)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 15)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 2) && (!memcmp(cur_modif, "zs", 2))) {
	      pc.glm_info.flags |= kfGlmZs;
	    } else if ((cur_modif_slen == 3) && (!memcmp(cur_modif, "sex", 3))) {
	      pc.glm_info.flags |= kfGlmSex;
	    } else if ((cur_modif_slen == 8) && (!memcmp(cur_modif, "no-x-sex", 8))) {
	      pc.glm_info.flags |= kfGlmNoXSex;
	    } else if ((cur_modif_slen == 9) && (!memcmp(cur_modif, "genotypic", 9))) {
	      pc.glm_info.flags |= kfGlmGenotypic;
	    } else if ((cur_modif_slen == 6) && (!memcmp(cur_modif, "hethom", 6))) {
	      pc.glm_info.flags |= kfGlmHethom;
	    } else if ((cur_modif_slen == 8) && (!memcmp(cur_modif, "dominant", 8))) {
	      pc.glm_info.flags |= kfGlmDominant;
	    } else if ((cur_modif_slen == 9) && (!memcmp(cur_modif, "recessive", 9))) {
	      pc.glm_info.flags |= kfGlmRecessive;
	    } else if ((cur_modif_slen == 11) && (!memcmp(cur_modif, "interaction", 11))) {
	      pc.glm_info.flags |= kfGlmInteraction;
	    } else if ((cur_modif_slen == 10) && (!memcmp(cur_modif, "hide-covar", 10))) {
	      pc.glm_info.flags |= kfGlmHideCovar;
	    } else if ((cur_modif_slen == 9) && (!memcmp(cur_modif, "intercept", 9))) {
	      pc.glm_info.flags |= kfGlmIntercept;
	    } else if ((cur_modif_slen == 14) && (!memcmp(cur_modif, "firth-fallback", 14))) {
	      pc.glm_info.flags |= kfGlmFirthFallback;
	    } else if ((cur_modif_slen == 5) && (!memcmp(cur_modif, "firth", 5))) {
	      pc.glm_info.flags |= kfGlmFirth;
	    } else if ((cur_modif_slen == 13) && (!memcmp(cur_modif, "standard-beta", 13))) {
	      logerrprint("Error: --glm 'standard-beta' modifier has been retired.  Use\n--{covar-}variance-standardize instead.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if ((cur_modif_slen == 4) && (!memcmp(cur_modif, "perm", 4))) {
	      pc.glm_info.flags |= kfGlmPerm;
	    } else if ((cur_modif_slen == 10) && (!memcmp(cur_modif, "perm-count", 10))) {
	      pc.glm_info.flags |= kfGlmPermCount;
	    } else if ((cur_modif_slen > 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      if (pc.glm_info.cols) {
		logerrprint("Error: Multiple --glm cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      reterr = parse_col_descriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0altcount\0totallele\0altcountcc\0totallelecc\0altfreq\0altfreqcc\0machr2\0firth\0test\0nobs\0beta\0orbeta\0se\0ci\0t\0p\0", flagname_p, kfGlmColChrom, kfGlmColDefault, 1, &pc.glm_info.cols);
	      if (reterr) {
		goto main_ret_1;
	      }
	      if ((!(pc.glm_info.cols & (kfGlmColBeta | kfGlmColOrbeta))) && ((pc.glm_info.cols & kfGlmColSe) || ((pc.glm_info.cols & kfGlmColCi) && (pc.ci_size != 0)))) {
		logerrprint("Error: --glm's 'se' and 'ci' columns require beta/orbeta to be included.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else if ((cur_modif_slen >= 5) && (!memcmp(cur_modif, "mperm", 5))) {
	      if ((cur_modif_slen < 7) || (cur_modif[5] != '=')) {
		logerrprint("Error: Improper --glm mperm syntax.  (Use --glm mperm=[value]'.)\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      if (scan_posint_defcap(cur_modif, &pc.glm_info.mperm_ct)) {
		sprintf(g_logbuf, "Error: Invalid --glm mperm parameter '%s'.\n", &(cur_modif[6]));
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	    } else if ((cur_modif_slen > 12) && (!memcmp(cur_modif, "local-covar=", 12))) {
	      if (pc.glm_local_covar_fname) {
		logerrprint("Error: Multiple --glm local-covar= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      reterr = alloc_fname(&(cur_modif[12]), "glm local-covar=", 0, &pc.glm_local_covar_fname);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else if ((cur_modif_slen > 10) && ((!memcmp(cur_modif, "local-bim=", 10)) || ((cur_modif_slen > 11) && (!memcmp(cur_modif, "local-pvar=", 11))))) {
	      if (pc.glm_local_pvar_fname) {
		logerrprint("Error: Multiple --glm local-pvar= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      const uint32_t is_pvar = (cur_modif[6] == 'p');
	      reterr = alloc_fname(&(cur_modif[10 + is_pvar]), "glm local-pvar=", 0, &pc.glm_local_pvar_fname);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else if ((cur_modif_slen > 10) && ((!memcmp(cur_modif, "local-fam=", 10)) || ((cur_modif_slen > 11) && (!memcmp(cur_modif, "local-psam=", 11))))) {
	      if (pc.glm_local_psam_fname) {
		logerrprint("Error: Multiple --glm local-psam= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      const uint32_t is_psam = (cur_modif[6] == 'p');
	      reterr = alloc_fname(&(cur_modif[10 + is_psam]), "glm local-psam=", 0, &pc.glm_local_psam_fname);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else if ((cur_modif_slen == 15) && (!memcmp(cur_modif, "local-omit-last", 15))) {
	      pc.glm_info.flags |= kfGlmLocalOmitLast;
	    } else if ((cur_modif_slen > 11) && (!memcmp(cur_modif, "local-cats=", 11))) {
	      if (pc.glm_info.local_cat_ct) {
		logerrprint("Error: Multiple --glm local-cats= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      if (scan_posint_capped(cur_modif, 4095, &pc.glm_info.local_cat_ct) || (pc.glm_info.local_cat_ct == 1)) {
		logerrprint("Error: Invalid --glm local-cats= category count (must be in [2, 4095]).\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --glm parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (!pc.glm_info.cols) {
	    pc.glm_info.cols = kfGlmColDefault;
	  }
	  if ((pc.glm_info.flags & (kfGlmSex | kfGlmNoXSex)) == (kfGlmSex | kfGlmNoXSex)) {
	    logerrprint("Error: Conflicting --glm parameters.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if ((pc.glm_info.flags & kfGlmPerm) && pc.glm_info.mperm_ct) {
	    logerrprint("Error: --glm 'perm' and 'mperm=' cannot be used together.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  uint32_t alternate_genotype_col_flags = (uint32_t)(pc.glm_info.flags & (kfGlmGenotypic | kfGlmHethom | kfGlmDominant | kfGlmRecessive));
	  if (alternate_genotype_col_flags) {
	    pc.xchr_model = 0;
	    if (alternate_genotype_col_flags & (alternate_genotype_col_flags - 1)) {
	      logerrprint("Error: Conflicting --glm parameters.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	  if ((pc.glm_info.flags & kfGlmIntercept) && (!(pc.glm_info.cols & kfGlmColTest))) {
	    logerrprint("Error: --glm 'intercept' modifier cannot be used with an omitted 'test' column.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (!pc.glm_local_covar_fname) {
	    if (pc.glm_local_pvar_fname || pc.glm_local_psam_fname) {
	      logerrprint("Error: Either all three --glm local-covar filenames must be specified, or none\nof them.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (pc.glm_info.flags & kfGlmLocalOmitLast) {
	      logerrprint("Error: --glm 'local-omit-last' must be used with 'local-covar='.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (pc.glm_info.local_cat_ct) {
	      logerrprint("Error: --glm 'local-cats=' must be used with 'local-covar='.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else {
	    if ((!pc.glm_local_pvar_fname) || (!pc.glm_local_psam_fname)) {
	      logerrprint("Error: Either all three --glm local-covar filenames must be specified, or none\nof them.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	  pc.command_flags1 |= kfCommand1Glm;
	} else if (!memcmp(flagname_p2, "en", 3)) {
	  if (load_params || (xload & (~kfXloadOxSample))) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct == 2) {
	    const char* cur_modif = argv[arg_idx + 2];
	    if (!strcmp(cur_modif, "ref-first")) {
	      oxford_import_flags |= kfOxfordImportRefFirst;
	    } else if (!strcmp(cur_modif, "ref-second")) {
	      oxford_import_flags |= kfOxfordImportRefSecond;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --gen parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  const char* cur_fname = argv[arg_idx + 1];
	  const uint32_t slen = strlen(cur_fname);
	  if (slen > kPglFnamesize - 1) {
	    logerrprint("Error: --gen filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pgenname, cur_fname, slen + 1);
	  xload |= kfXloadOxGen;
	} else if (!memcmp(flagname_p2, "enotyping-rate", 15)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    if (strcmp("dosage", cur_modif)) {
	      sprintf(g_logbuf, "Error: Invalid --genotyping-rate parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    pc.misc_flags |= kfMiscGenotypingRateDosage;
	  }
	  pc.command_flags1 |= kfCommand1GenotypingRate;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;

      case 'h':
	if (!memcmp(flagname_p2, "ardy", 5)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 3)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "zs")) {
	      pc.hardy_modifier |= kfHardyZs;
	    } else if (!strcmp(cur_modif, "midp")) {
	      pc.hardy_modifier |= kfHardyMidp;
	    } else if ((strlen(cur_modif) > 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      if (pc.hardy_modifier & kfHardyColAll) {
		logerrprint("Error: Multiple --hardy cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }	      
	      reterr = parse_col_descriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0gcounts\0gcount1col\0hetfreq\0sexaf\0femalep\0p\0", "freq", kfHardyColChrom, kfHardyColDefault, 1, &pc.hardy_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	      if ((pc.hardy_modifier & (kfHardyColGcounts | kfHardyColGcount1col)) == (kfHardyColGcounts | kfHardyColGcount1col)) {
		logerrprint("Error: --hardy's gcounts and gcounts1col column sets are mutually exclusive.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --hardy parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (!(pc.hardy_modifier & kfHardyColAll)) {
	    pc.hardy_modifier |= kfHardyColDefault;
	  }
	  pc.command_flags1 |= kfCommand1Hardy;
	} else if (!memcmp(flagname_p2, "we", 3)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 3)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "midp")) {
	      pc.misc_flags |= kfMiscHweMidp;
	    } else if (!strcmp(cur_modif, "keep-fewhet")) {
	      pc.misc_flags |= kfMiscHweKeepFewhet;
	    } else {
	      if ((pc.hwe_thresh != 1.0) || (!scanadv_double(cur_modif, &pc.hwe_thresh))) {
		logerrprint("Error: Invalid --hwe parameter sequence.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      if ((pc.hwe_thresh < 0.0) || (pc.hwe_thresh >= 1.0)) {
		sprintf(g_logbuf, "Error: Invalid --hwe threshold '%s' (must be in [0, 1)).\n", cur_modif);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	    }
	  }
	  if (pc.hwe_thresh == 1.0) {
	    logerrprint("Error: --hwe requires a p-value threshold.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if ((pc.misc_flags & kfMiscHweMidp) && (pc.hwe_thresh >= 0.5)) {
	    logerrprint("Error: --hwe threshold must be smaller than 0.5 when using mid-p adjustment.\n");
	  }
	  pc.filter_flags |= kfFilterAllReq | kfFilterNoSplitChr;
	} else if (!memcmp(flagname_p2, "ard-call-threshold", 19)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  double hard_call_frac;
	  if ((!scanadv_double(cur_modif, &hard_call_frac)) || (hard_call_frac < 0.0) || (hard_call_frac >= (0.5 - kSmallEpsilon))) {
	    sprintf(g_logbuf, "Error: Invalid --hard-call-threshold parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  pc.hard_call_thresh = (int32_t)(hard_call_frac * ((1 + kSmallEpsilon) * kDosageMid));
	} else if (!memcmp(flagname_p2, "orse", 5)) {
	  if (chr_info.chrset_source) {
	    logerrprint("Error: Conflicting chromosome-set flags.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  chr_info.chrset_source = kChrsetSourceCmdline;
	  chr_info.autosome_ct = 31;
	  chr_info.xymt_codes[0] = 32;
	  chr_info.xymt_codes[1] = 33;
	  chr_info.xymt_codes[2] = -2;
	  chr_info.xymt_codes[3] = -2;
	  chr_info.xymt_codes[4] = -2;
	  chr_info.xymt_codes[5] = -2;
#ifdef __LP64__
	  chr_info.haploid_mask[0] = 0x300000000LLU;
#else
	  chr_info.haploid_mask[0] = 0;
	  chr_info.haploid_mask[1] = 3;
#endif
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "uman", 5)) {
	  if (chr_info.chrset_source) {
	    logerrprint("Error: Conflicting chromosome-set flags.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  chr_info.chrset_source = kChrsetSourceCmdline;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "aps", 4)) {
	  if (load_params || xload) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct == 2) {
	    const char* cur_modif = argv[arg_idx + 2];
	    if (!strcmp(cur_modif, "ref-first")) {
	      oxford_import_flags |= kfOxfordImportRefFirst;
	    } else if (!strcmp(cur_modif, "ref-second")) {
	      oxford_import_flags |= kfOxfordImportRefSecond;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --haps parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  const char* cur_fname = argv[arg_idx + 1];
	  const uint32_t slen = strlen(cur_fname);
	  if (slen > kPglFnamesize - 1) {
	    logerrprint("Error: --haps filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pgenname, cur_fname, slen + 1);
	  xload |= kfXloadOxHaps;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'i':
	if (!memcmp(flagname_p2, "ndiv-sort", 10)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 3)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* indiv_sort_mode_str = argv[arg_idx + 1];
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
	      sprintf(g_logbuf, "Error: Missing '--indiv-sort %s' filename.\n", indiv_sort_mode_str);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    pc.sample_sort_flags = kfSortFile;
	    uint32_t fname_modif_idx = 2;
	    if (param_ct == 3) {
	      if (check_extra_param(&(argv[arg_idx]), "sid", &fname_modif_idx)) {
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.sample_sort_flags |= kfSortFileSid;
	    }
	    reterr = alloc_fname(argv[arg_idx + fname_modif_idx], flagname_p, 0, &pc.sample_sort_fname);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  } else {
	    sprintf(g_logbuf, "Error: '%s' is not a valid mode for --indiv-sort.\n", indiv_sort_mode_str);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if ((param_ct > 1) && (!(pc.sample_sort_flags & kfSortFile))) {
	    sprintf(g_logbuf, "Error: '--indiv-sort %s' does not accept additional parameters.\n", indiv_sort_mode_str);
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	} else if (!memcmp(flagname_p2, "d-delim", 8)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    id_delim = extract_char_param(argv[arg_idx + 1]);
	    if (!id_delim) {
	      logerrprint("Error: --id-delim parameter must be a single character.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (((unsigned char)id_delim) < ' ') {
	      logerrprint("Error: --id-delim parameter cannot be tab, newline, or a nonprinting character.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	  } else {
	    id_delim = '_';
	  }
	} else if ((!memcmp(flagname_p2, "ndep-pairwise", 14)) || (!memcmp(flagname_p2, "ndep-pairphase", 15))) {
	  if (pc.command_flags1 & kfCommand1LdPrune) {
	    logerrprint("Error: Multiple LD pruning commands.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 2, 4)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  double first_paramd;
	  char* first_param_end = scanadv_double(cur_modif, &first_paramd);
	  if ((!first_param_end) || (first_paramd < 0.0)) {
	    sprintf(g_logbuf, "Error: Invalid --%s window size '%s'.\n", flagname_p, cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  uint32_t is_kb = 0;
	  uint32_t next_param_idx = 2;
	  if (match_upper_counted(first_param_end, "KB", 2) && (!first_param_end[2])) {
	    is_kb = 1;
	  } else if (match_upper_counted(argv[arg_idx + 2], "KB", 2) && (!argv[arg_idx + 2][2])) {
	    is_kb = 1;
	    next_param_idx = 3;
	  }
	  if (is_kb) {
	    pc.ld_info.prune_modifier |= kfLdPruneWindowBp;
	    if (first_paramd > 2147483.646) {
	      pc.ld_info.prune_window_size = 2147483646;
	    } else {
	      pc.ld_info.prune_window_size = ((int32_t)(first_paramd * 1000 * (1 + kSmallEpsilon)));
	      if (pc.ld_info.prune_window_size < 2) {
		sprintf(g_logbuf, "Error: --%s window size cannot be smaller than 2.\n", flagname_p);
		goto main_ret_INVALID_CMDLINE_2A;
	      }
	    }
	  } else {
	    if (first_paramd > 2147483647) {
	      pc.ld_info.prune_window_size = 2147483647;
	    } else {
	      pc.ld_info.prune_window_size = ((int32_t)first_paramd);
	    }
	  }
	  if (next_param_idx + 2 == param_ct) {
	    sprintf(g_logbuf, "Error: Invalid --%s parameter sequence.\n", flagname_p);
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (next_param_idx < param_ct) {
	    // explicit step size
	    cur_modif = argv[arg_idx + next_param_idx];
	    if (scan_posint_defcap(cur_modif, &pc.ld_info.prune_window_incr)) {
	      sprintf(g_logbuf, "Error: Invalid --%s window-increment '%s'.\n", flagname_p, cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (!is_kb) {
	      if (pc.ld_info.prune_window_incr > pc.ld_info.prune_window_size) {
	        sprintf(g_logbuf, "Error: --%s window-increment cannot be larger than window size.\n", flagname_p);
	        goto main_ret_INVALID_CMDLINE_2A;
	      }
	    } else if (pc.ld_info.prune_window_incr != 1) {
	      sprintf(g_logbuf, "Error: --%s window-increment must be 1 when window size is in\nkilobase units.\n", flagname_p);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	  } else {
	    pc.ld_info.prune_window_incr = 1;
	  }
	  cur_modif = argv[arg_idx + param_ct];
	  if ((!scanadv_double(cur_modif, &pc.ld_info.prune_last_param)) || (pc.ld_info.prune_last_param < 0.0) || (pc.ld_info.prune_last_param >= 1.0)) {
	    sprintf(g_logbuf, "Error: Invalid --%s r^2 threshold '%s'.\n", flagname_p2, cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  pc.command_flags1 |= kfCommand1LdPrune;
	  if (flagname_p2[9] == 'p') {
	    pc.ld_info.prune_modifier |= kfLdPrunePairphase;
	  } else {
	    pc.ld_info.prune_modifier |= kfLdPrunePairwise;
	  }
	} else if (!memcmp(flagname_p2, "nput-missing-genotype", 22)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  input_missing_geno_char = extract_char_param(cur_modif);
	  if (((unsigned char)input_missing_geno_char) <= ' ') {
	    sprintf(g_logbuf, "Error: Invalid --input-missing-genotype parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else if (!memcmp(flagname_p2, "nput-missing-phenotype", 23)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  double dxx;
	  if (scan_int32(cur_modif, &pc.missing_pheno) || ((pc.missing_pheno >= 0) && (pc.missing_pheno <= 2)) || (!scanadv_double(cur_modif, &dxx)) || (dxx != ((double)pc.missing_pheno))) {
	    sprintf(g_logbuf, "Error: Invalid --input-missing-phenotype parameter '%s' (must be an integer in [-2147483647, -1] or [3, 2147483647]).\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else if (!memcmp(flagname_p2, "mport-dosage-certainty", 23)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  if ((!scanadv_double(cur_modif, &import_dosage_certainty)) || (import_dosage_certainty < 0.0) || (import_dosage_certainty > 1.0)) {
	    sprintf(g_logbuf, "Error: Invalid --import-dosage-certainty parameter '%s' (must be in [0, 1]).\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  import_dosage_certainty *= 1.0 - kSmallEpsilon;
	} else if (!memcmp(flagname_p2, "mport-dosage", 13)) {
	  if (load_params || xload) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 10)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t format_num_m1 = 3;
	  for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 8) && (!memcmp(cur_modif, "noheader", 8))) {
	      plink1_dosage_info.flags |= kfPlink1DosageNoheader;
	    } else if ((cur_modif_slen > 6) && (!memcmp(cur_modif, "skip", 4)) && (cur_modif[4] >= '0') && (cur_modif[4] <= '2') && (cur_modif[5] == '=')) {
	      const uint32_t skip_idx = (uint32_t)((unsigned char)cur_modif[4]) - 48;
	      if (plink1_dosage_info.skips[skip_idx]) {
		LOGERRPRINTF("Error: Multiple --import-dosage skip%u= modifiers.\n", skip_idx);
		goto main_ret_INVALID_CMDLINE;
	      }
	      if (scan_uint_capped(&(cur_modif[6]), kMaxLongLine / 2, &(plink1_dosage_info.skips[skip_idx]))) {
		sprintf(g_logbuf, "Error: Invalid --import-dosage skip%u= parameter '%s'.\n", skip_idx, &(cur_modif[6]));
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	    } else if ((cur_modif_slen == 5) && (!memcmp(cur_modif, "dose1", 5))) {
	      plink1_dosage_info.flags |= kfPlink1DosageFormatSingle01;
	    } else if ((cur_modif_slen == 8) && (!memcmp(cur_modif, "format=", 7))) {
	      if (format_num_m1 != 3) {
	        logerrprint("Error: Multiple --import-dosage format= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      format_num_m1 = (uint32_t)((unsigned char)cur_modif[7]) - 49;
	      if (format_num_m1 >= 3) {
		sprintf(g_logbuf, "Error: Invalid --import-dosage format= parameter '%c'.\n", cur_modif[7]);
		goto main_ret_INVALID_CMDLINE_2A;
	      }
	    } else if ((cur_modif_slen == 9) && (!memcmp(cur_modif, "ref-first", 9))) {
	      plink1_dosage_info.flags |= kfPlink1DosageRefFirst;
	    } else if ((cur_modif_slen == 10) && (!memcmp(cur_modif, "ref-second", 10))) {
	      plink1_dosage_info.flags |= kfPlink1DosageRefSecond;
	    } else if ((cur_modif_slen > 11) && (!memcmp(cur_modif, "single-chr=", 11))) {
	      if (import_single_chr_str) {
		logerrprint("Error: Multiple --import-dosage single-chr= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      const char* chr_code = &(cur_modif[11]);
	      if (!(pc.misc_flags & kfMiscAllowExtraChrs)) {
		if (get_chr_code_raw(chr_code) < 0) {
		  sprintf(g_logbuf, "Error: Invalid --import-dosage single-chr= chromosome code '%s'. (Did you forget --allow-extra-chr?)\n", chr_code);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
	      }
	      reterr = cmdline_alloc_string(chr_code, argv[arg_idx], kMaxIdSlen, &import_single_chr_str);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else if ((cur_modif_slen > 12) && (!memcmp(cur_modif, "chr-col-num=", 12))) {
	      if (plink1_dosage_info.chr_col_idx != 0xffffffffU) {
		logerrprint("Error: Multiple --import-dosage chr-col-num= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      uint32_t uii;
	      if (scan_posint_capped(&(cur_modif[12]), kMaxLongLine / 2, &uii)) {
		sprintf(g_logbuf, "Error: Invalid --import-dosage chr-col-num= parameter '%s'.\n", &(cur_modif[12]));
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      plink1_dosage_info.chr_col_idx = uii - 1;
	    } else if ((cur_modif_slen > 12) && (!memcmp(cur_modif, "pos-col-num=", 12))) {
	      if (plink1_dosage_info.pos_col_idx != 0xffffffffU) {
		logerrprint("Error: Multiple --import-dosage pos-col-num= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      uint32_t uii;
	      if (scan_posint_capped(&(cur_modif[12]), kMaxLongLine / 2, &uii)) {
		sprintf(g_logbuf, "Error: Invalid --import-dosage pos-col-num= parameter '%s'.\n", &(cur_modif[12]));
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      plink1_dosage_info.pos_col_idx = uii - 1;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --import-dosage parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }

	  if (!format_num_m1) {
	    plink1_dosage_info.flags |= kfPlink1DosageFormatSingle;
	  } else {
	    if (plink1_dosage_info.flags & kfPlink1DosageFormatSingle01) {
	      logerrprint("Error: --import-dosage 'dose1' modifier must be used with 'format=1'.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (format_num_m1 == 2) {
	      plink1_dosage_info.flags |= kfPlink1DosageFormatTriple;
	    }
	  }
	  if ((plink1_dosage_info.flags & (kfPlink1DosageRefFirst | kfPlink1DosageRefSecond)) == (kfPlink1DosageRefFirst | kfPlink1DosageRefSecond)) {
	    logerrprint("Error: --import-dosage 'ref-first' and 'ref-second' modifiers cannot be used\ntogether.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  const uint32_t id_col_idx = plink1_dosage_info.skips[0];
	  const uint32_t a1_col_idx = id_col_idx + plink1_dosage_info.skips[1] + 1;
	  const uint32_t data_col_idx = a1_col_idx + plink1_dosage_info.skips[2] + 2;
	  const uint32_t chr_col_idx = plink1_dosage_info.chr_col_idx;
	  if (chr_col_idx != 0xffffffffU) {
	    if (import_single_chr_str) {
	      logerrprint("Error: --import-dosage 'single-chr=' and 'chr-col-num=' modifiers cannot be\nused together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if ((chr_col_idx == id_col_idx) || (chr_col_idx == a1_col_idx) || (chr_col_idx == a1_col_idx + 1)) {
	      logerrprint("Error: --import-dosage chr-col-num= value collides with another column.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (chr_col_idx >= data_col_idx) {
	      logerrprint("Error: --import-dosage chr-col-num= value too large.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	  const uint32_t pos_col_idx = plink1_dosage_info.pos_col_idx;
	  if (pos_col_idx != 0xffffffffU) {
	    if ((pos_col_idx == id_col_idx) || (pos_col_idx == a1_col_idx) || (pos_col_idx == a1_col_idx + 1) || (pos_col_idx == chr_col_idx)) {
	      logerrprint("Error: --import-dosage pos-col-num= value collides with another column.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (pos_col_idx >= data_col_idx) {
	      logerrprint("Error: --import-dosage pos-col-num= value too large.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  const uint32_t slen = strlen(cur_modif);
	  if (slen > kPglFnamesize - 1) {
	    logerrprint("Error: --import-dosage filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pgenname, cur_modif, slen + 1);
	  xload = kfXloadPlink1Dosage;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;

      case 'k':
	if (!memcmp(flagname_p2, "eep", 4)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const uint32_t sid_present = !strcmp(argv[arg_idx + 1], "sid");
	  if (sid_present) {
	    if (param_ct == 1) {
	      logerrprint("Error: '--keep sid' requires at least one filename.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    pc.misc_flags |= kfMiscKeepfileSid;
	  }
	  reterr = alloc_and_flatten(&(argv[arg_idx + 1 + sid_present]), param_ct - sid_present, kPglFnamesize, &pc.keep_fnames);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "eep-fam", 8)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, kPglFnamesize, &pc.keepfam_fnames);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "eep-autoconv", 13)) {
	  pc.misc_flags |= kfMiscKeepAutoconv;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "eep-females", 12)) {
	  pc.filter_flags |= kfFilterPsamReq | kfFilterExclMales | kfFilterExclNosex;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "eep-males", 10)) {
	  pc.filter_flags |= kfFilterPsamReq | kfFilterExclFemales | kfFilterExclNosex;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "eep-nosex", 10)) {
	  pc.filter_flags |= kfFilterPsamReq | kfFilterExclFemales | kfFilterExclMales;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "eep-founders", 13)) {
	  pc.filter_flags |= kfFilterPsamReq | kfFilterExclNonfounders;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "eep-nonfounders", 16)) {
	  if (pc.filter_flags & kfFilterExclNonfounders) {
	    logerrprint("Error: --keep-nonfounders cannot be used with --keep-founders.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  pc.filter_flags |= kfFilterPsamReq | kfFilterExclFounders;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "ing-cutoff", 11)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct == 2) {
	    // .king.id, .king.bin appended
	    reterr = alloc_fname(argv[arg_idx + 1], flagname_p, 9, &king_cutoff_fprefix);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  char* cur_modif = argv[arg_idx + param_ct];
	  if ((!scanadv_double(cur_modif, &pc.king_cutoff)) || (pc.king_cutoff < 0.0) || (pc.king_cutoff >= 0.5)) {
	    sprintf(g_logbuf, "Error: Invalid --king-cutoff parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  pc.command_flags1 |= kfCommand1KingCutoff;
	} else if (!memcmp(flagname_p2, "ing-table-filter", 17)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  if ((!scanadv_double(cur_modif, &pc.king_table_filter)) || (pc.king_table_filter > 0.5)) {
	    sprintf(g_logbuf, "Error: Invalid --king-table-filter parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else if (!memcmp(flagname_p2, "eep-if", 7)) {
	  reterr = validate_and_alloc_cmp_expr(&(argv[arg_idx + 1]), argv[arg_idx], param_ct, &pc.keep_if_expr);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "eep-cats", 9)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_fname(argv[arg_idx + 1], flagname_p, 0, &pc.keep_cats_fname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "eep-cat-names", 14)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, kMaxIdBlen, &pc.keep_cat_names_flattened);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "eep-cat-pheno", 14)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + 1], argv[arg_idx], kMaxIdSlen, &pc.keep_cat_phenoname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "eep-allele-order", 17)) {
	  logprint("Note: --keep-allele-order no longer has any effect.\n");
	  goto main_param_zero;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;

      case 'l':
	if (!memcmp(flagname_p2, "ambda", 6)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (!scanadv_double(argv[arg_idx + 1], &pc.adjust_info.lambda)) {
	    sprintf(g_logbuf, "Error: Invalid --lambda parameter '%s'.\n", argv[arg_idx + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (pc.adjust_info.lambda < 1.0) {
	    logprint("Note: --lambda parameter set to 1.\n");
	    pc.adjust_info.lambda = 1.0;
	  }
	} else if (!memcmp(flagname_p2, "egend", 6)) {
	  if (load_params || (xload & (~kfXloadOxHaps))) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (!xload) {
	    logerrprint("Error: --legend must be used with --haps.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 2, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_fname = argv[arg_idx + 1];
	  uint32_t slen = strlen(cur_fname);
	  if (slen > kPglFnamesize - 1) {
	    logerrprint("Error: --legend filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pvarname, cur_fname, slen + 1);
	  const char* chr_code = argv[arg_idx + 2];
	  if (!(pc.misc_flags & kfMiscAllowExtraChrs)) {
	    if (get_chr_code_raw(chr_code) < 0) {
	      sprintf(g_logbuf, "Error: Invalid --legend chromosome code '%s'. (Did you forget --allow-extra-chr?)\n", chr_code);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  reterr = cmdline_alloc_string(chr_code, argv[arg_idx], kMaxIdSlen, &import_single_chr_str);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  xload |= kfXloadOxLegend;
	} else if (!memcmp(flagname_p2, "oop-assoc", 10)) {
	  logerrprint("Error: --loop-assoc is retired.  Use --within + --split-cat-pheno instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;

      case 'm':
	if (!memcmp(flagname_p2, "emory", 6)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t mb_modif_idx = 1;
	  if (param_ct == 2) {
	    if (check_extra_param(&(argv[arg_idx]), "require", &mb_modif_idx)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    memory_require = 1;
	  }
	  const char* mb_modif = argv[arg_idx + mb_modif_idx];
	  if (scan_posintptr(mb_modif, (uintptr_t*)(&malloc_size_mb))) {
	    sprintf(g_logbuf, "Error: Invalid --memory parameter '%s'.\n", mb_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (malloc_size_mb < (intptr_t)kBigstackMinMb) {
	    sprintf(g_logbuf, "Error: Invalid --memory parameter '%s' (minimum %u).\n", mb_modif, kBigstackMinMb);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
#ifndef __LP64__
	  if (malloc_size_mb > (intptr_t)kMalloc32bitMbMax) {
	    LOGERRPRINTF("Error: --memory parameter too large for 32-bit version (max %u).\n", kMalloc32bitMbMax);
	    goto main_ret_INVALID_CMDLINE;
	  }
#endif
	} else if (!memcmp(flagname_p2, "ake-bed", 8)) {
	  if (pc.exportf_modifier & kfExportfIndMajorBed) {
	    logerrprint("Error: --make-bed cannot be used with --export ind-major-bed.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 3)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 3) && (!memcmp(cur_modif, "vzs", 3))) {
	      make_plink2_modifier |= kfMakeBimZs;
	    } else if ((cur_modif_slen == 9) && (!memcmp(cur_modif, "trim-alts", 9))) {
	      make_plink2_modifier |= kfMakePlink2TrimAlts;
	    } else if (((cur_modif_slen > 2) && (!memcmp(cur_modif, "m=", 2))) || ((cur_modif_slen > 14) && (!memcmp(cur_modif, "multiallelics=", 14)))) {
	      if (make_plink2_modifier & kfMakePlink2MMask) {
		logerrprint("Error: Multiple --make-bed multiallelics= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      const char* mode_start = (cur_modif[1] == '=')? (&(cur_modif[2])) : (&(cur_modif[14]));
	      if (!strcmp(mode_start, "-")) {
		make_plink2_modifier |= kfMakePlink2MSplitAll;
	      } else if (!strcmp(mode_start, "-snps")) {
		make_plink2_modifier |= kfMakePlink2MSplitSnps;
	      } else if ((!strcmp(mode_start, "+")) || (!strcmp(mode_start, "+both"))) {
		make_plink2_modifier |= kfMakePlink2MMergeBoth;
	      } else if (!strcmp(mode_start, "+snps")) {
		make_plink2_modifier |= kfMakePlink2MMergeSnps;
	      } else if (!strcmp(mode_start, "+any")) {
		make_plink2_modifier |= kfMakePlink2MMergeAny;
	      } else {
		sprintf(g_logbuf, "Error: Invalid --make-bed multiallelics= mode '%s'.\n", mode_start);
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
	  make_plink2_modifier |= kfMakeBed | kfMakeBim | kfMakeFam;
	  pc.command_flags1 |= kfCommand1MakePlink2;
	} else if (!memcmp(flagname_p2, "ake-bpgen", 10)) {
	  if (make_plink2_modifier & kfMakeBed) {
	    logerrprint("Error: --make-bpgen cannot be used with --make-bed.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (pc.misc_flags & kfMiscKeepAutoconv) {
	    logerrprint("Error: --make-bpgen cannot be used with --keep-autoconv.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 3) && (!memcmp(cur_modif, "vzs", 3))) {
	      make_plink2_modifier |= kfMakeBimZs;
	    } else if ((cur_modif_slen > 7) && (!memcmp(cur_modif, "format=", 7))) {
	      if (make_plink2_modifier & (kfMakePgenFormatBase * 3)) {
		logerrprint("Error: Multiple --make-bpgen format= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      const uint32_t fcode_minus_2 = ((uint32_t)((unsigned char)cur_modif[7])) - 50;
	      if ((fcode_minus_2 > 2) || cur_modif[8]) {
		sprintf(g_logbuf, "Error: Invalid --make-bpgen format code '%s'.\n", &(cur_modif[7]));
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      if (fcode_minus_2) {
		logerrprint("Error: --make-bpgen formats 3 and 4 (unphased/phased dosage) are not\nimplemented yet.\n");
		reterr = kPglRetNotYetSupported;
		goto main_ret_1;
	      }
	      make_plink2_modifier = (make_plink2_t)(make_plink2_modifier | (kfMakePgenFormatBase * (1 + fcode_minus_2)));
	    } else if (((cur_modif_slen > 2) && (!memcmp(cur_modif, "m=", 2))) && ((cur_modif_slen > 14) && (!memcmp(cur_modif, "multiallelics=", 14)))) {
	      if (make_plink2_modifier & kfMakePlink2MMask) {
		logerrprint("Error: Multiple --make-bpgen multiallelics= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      const char* mode_start = (cur_modif[1] == '=')? (&(cur_modif[2])) : (&(cur_modif[14]));
	      if (!strcmp(mode_start, "-")) {
		make_plink2_modifier |= kfMakePlink2MSplitAll;
	      } else if (!strcmp(mode_start, "-snps")) {
		make_plink2_modifier |= kfMakePlink2MSplitSnps;
	      } else if ((!strcmp(mode_start, "+")) || (!strcmp(mode_start, "+both"))) {
		make_plink2_modifier |= kfMakePlink2MMergeBoth;
	      } else if (!strcmp(mode_start, "+snps")) {
		make_plink2_modifier |= kfMakePlink2MMergeSnps;
	      } else if (!strcmp(mode_start, "+any")) {
		make_plink2_modifier |= kfMakePlink2MMergeAny;
	      } else {
		sprintf(g_logbuf, "Error: Invalid --make-bpgen multiallelics= mode '%s'.\n", mode_start);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	    } else if ((cur_modif_slen == 9) && (!memcmp(cur_modif, "trim-alts", 9))) {
	      make_plink2_modifier |= kfMakePlink2TrimAlts;
	    } else if ((cur_modif_slen == 11) && (!memcmp(cur_modif, "erase-alt2+", 11))) {
	      make_plink2_modifier |= kfMakePgenEraseAlt2Plus;
	    } else if ((cur_modif_slen == 11) && (!memcmp(cur_modif, "erase-phase", 11))) {
	      make_plink2_modifier |= kfMakePgenErasePhase;
	    } else if ((cur_modif_slen == 12) && (!memcmp(cur_modif, "erase-dosage", 12))) {
	      make_plink2_modifier |= kfMakePgenEraseDosage;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-bpgen parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  make_plink2_modifier |= kfMakePgen | kfMakeBim | kfMakeFam;
	  pc.command_flags1 |= kfCommand1MakePlink2;
	} else if (!memcmp(flagname_p2, "ake-pgen", 9)) {
	  if (make_plink2_modifier & (kfMakeBed | kfMakePgen)) {
	    logerrprint("Error: --make-pgen cannot be used with --make-bed/--make-bpgen.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (make_plink2_modifier & (kfMakeBim | kfMakeFam | kfMakePvar | kfMakePsam)) {
	    logerrprint("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (pc.misc_flags & kfMiscKeepAutoconv) {
	    logerrprint("Error: --make-pgen cannot be used with --keep-autoconv.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 4)) {
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  uint32_t explicit_pvar_cols = 0;
	  uint32_t explicit_psam_cols = 0;
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 3) && (!memcmp(cur_modif, "vzs", 3))) {
	      pc.pvar_psam_modifier |= kfPvarZs;
	    } else if ((cur_modif_slen >= 10) && (!memcmp(cur_modif, "pvar-cols=", 10))) {
	      if (explicit_pvar_cols) {
		logerrprint("Error: Multiple --make-pgen pvar-cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      explicit_pvar_cols = 1;
	      reterr = parse_col_descriptor(&(cur_modif[10]), "xheader\0maybequal\0qual\0maybefilter\0filter\0maybeinfo\0info\0maybecm\0cm\0", "--make-pgen pvar-cols", kfPvarColXheader, kfPvarColDefault, 0, &pc.pvar_psam_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	      if ((pc.pvar_psam_modifier & kfPvarColXinfo) && (!(pc.pvar_psam_modifier & kfPvarColXheader))) {
		logerrprint("Error: --make-pgen pvar-cols= expression cannot exclude xheader when info is\npresent.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else if ((cur_modif_slen > 7) && (!memcmp(cur_modif, "format=", 7))) {
	      if (make_plink2_modifier & (kfMakePgenFormatBase * 3)) {
		logerrprint("Error: Multiple --make-pgen format= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      const uint32_t fcode_minus_2 = ((uint32_t)((unsigned char)cur_modif[7])) - 50;
	      if ((fcode_minus_2 > 2) || cur_modif[8]) {
		sprintf(g_logbuf, "Error: Invalid --make-pgen format code '%s'.\n", &(cur_modif[7]));
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      if (fcode_minus_2) {
		logerrprint("Error: --make-pgen formats 3 and 4 (unphased/phased dosage) are not implemented\nyet.\n");
		reterr = kPglRetNotYetSupported;
		goto main_ret_1;
	      }
	      make_plink2_modifier = (make_plink2_t)(make_plink2_modifier | (kfMakePgenFormatBase * (1 + fcode_minus_2)));
	    } else if (((cur_modif_slen > 2) && (!memcmp(cur_modif, "m=", 2))) && ((cur_modif_slen > 14) && (!memcmp(cur_modif, "multiallelics=", 14)))) {
	      if (make_plink2_modifier & kfMakePlink2MMask) {
		logerrprint("Error: Multiple --make-pgen multiallelics= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      const char* mode_start = (cur_modif[1] == '=')? (&(cur_modif[2])) : (&(cur_modif[14]));
	      if (!strcmp(mode_start, "-")) {
		make_plink2_modifier |= kfMakePlink2MSplitAll;
	      } else if (!strcmp(mode_start, "-snps")) {
		make_plink2_modifier |= kfMakePlink2MSplitSnps;
	      } else if ((!strcmp(mode_start, "+")) || (!strcmp(mode_start, "+both"))) {
		make_plink2_modifier |= kfMakePlink2MMergeBoth;
	      } else if (!strcmp(mode_start, "+snps")) {
		make_plink2_modifier |= kfMakePlink2MMergeSnps;
	      } else if (!strcmp(mode_start, "+any")) {
		make_plink2_modifier |= kfMakePlink2MMergeAny;
	      } else {
		sprintf(g_logbuf, "Error: Invalid --make-pgen multiallelics= mode '%s'.\n", mode_start);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	    } else if ((cur_modif_slen == 9) && (!memcmp(cur_modif, "trim-alts", 9))) {
	      make_plink2_modifier |= kfMakePlink2TrimAlts;
	    } else if ((cur_modif_slen == 11) && (!memcmp(cur_modif, "erase-alt2+", 11))) {
	      make_plink2_modifier |= kfMakePgenEraseAlt2Plus;
	    } else if ((cur_modif_slen == 11) && (!memcmp(cur_modif, "erase-phase", 11))) {
	      make_plink2_modifier |= kfMakePgenErasePhase;
	    } else if ((cur_modif_slen == 12) && (!memcmp(cur_modif, "erase-dosage", 12))) {
	      make_plink2_modifier |= kfMakePgenEraseDosage;
	    } else if ((cur_modif_slen >= 10) && (!memcmp(cur_modif, "psam-cols=", 10))) {
	      if (explicit_psam_cols) {
		logerrprint("Error: Multiple --make-pgen psam-cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      explicit_psam_cols = 1;
	      reterr = parse_col_descriptor(&(cur_modif[10]), "maybesid\0sid\0maybeparents\0parents\0sex\0pheno1\0phenos\0", "make-pgen psam-cols", kfPsamColMaybesid, kfPsamColDefault, 0, &pc.pvar_psam_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-pgen parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (!explicit_pvar_cols) {
	    pc.pvar_psam_modifier |= kfPvarColDefault;
	  }
	  if (!explicit_psam_cols) {
	    pc.pvar_psam_modifier |= kfPsamColDefault;
	  }
	  make_plink2_modifier |= kfMakePgen | kfMakePvar | kfMakePsam;
	  pc.command_flags1 |= kfCommand1MakePlink2;
	} else if (!memcmp(flagname_p2, "ake-just-bim", 13)) {
	  if (make_plink2_modifier & (kfMakeBed | kfMakePgen)) {
	    logerrprint("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    if (!strcmp(cur_modif, "zs")) {
	      make_plink2_modifier |= kfMakeBimZs;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-just-bim parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  make_plink2_modifier |= kfMakeBim;
	  pc.command_flags1 |= kfCommand1MakePlink2;
	} else if (!memcmp(flagname_p2, "ake-just-fam", 13)) {
	  if (make_plink2_modifier & (kfMakeBed | kfMakePgen)) {
	    logerrprint("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  make_plink2_modifier |= kfMakeFam;
	  pc.command_flags1 |= kfCommand1MakePlink2;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "ake-just-pvar", 14)) {
	  if (make_plink2_modifier & (kfMakeBed | kfMakePgen)) {
	    logerrprint("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t explicit_cols = 0;
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 2) && (!memcmp(cur_modif, "zs", 2))) {
	      pc.pvar_psam_modifier |= kfPvarZs;
	    } else if ((cur_modif_slen >= 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      if (explicit_cols) {
		logerrprint("Error: Multiple --make-just-pvar cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      explicit_cols = 1;
	      reterr = parse_col_descriptor(&(cur_modif[5]), "xheader\0maybequal\0qual\0maybefilter\0filter\0maybeinfo\0info\0maybecm\0cm\0", "make-just-pvar", kfPvarColXheader, kfPvarColDefault, 0, &pc.pvar_psam_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	      if ((pc.pvar_psam_modifier & kfPvarColXinfo) && (!(pc.pvar_psam_modifier & kfPvarColXheader))) {
		logerrprint("Error: --make-just-pvar cols= expression cannot exclude xheader when info is\npresent.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-just-pvar parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (!explicit_cols) {
	    pc.pvar_psam_modifier |= kfPvarColDefault;
	  }
	  make_plink2_modifier |= kfMakePvar;
	  pc.command_flags1 |= kfCommand1MakePlink2;
	} else if (!memcmp(flagname_p2, "ake-just-psam", 14)) {
	  if (make_plink2_modifier & (kfMakeBed | kfMakePgen)) {
	    logerrprint("Error: --make-just-... cannot be used with --make-bed/--make-{b}pgen.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    if ((strlen(cur_modif) >= 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      reterr = parse_col_descriptor(&(cur_modif[5]), "maybesid\0sid\0maybeparents\0parents\0sex\0pheno1\0phenos\0", "make-just-psam", kfPsamColMaybesid, kfPsamColDefault, 0, &pc.pvar_psam_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-just-psam parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  } else {
	    pc.pvar_psam_modifier |= kfPsamColDefault;
	  }
	  make_plink2_modifier |= kfMakePsam;
	  pc.command_flags1 |= kfCommand1MakePlink2;
	} else if (!memcmp(flagname_p2, "ake-king", 9)) {
	  // may want to add options for handling X/Y/MT
	  if (king_cutoff_fprefix) {
	    logerrprint("Error: --make-king cannot be used with a --king-cutoff input fileset.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "zs")) {
	      if (pc.king_modifier & kfKingMatrixEncodemask) {
		logerrprint("Error: Multiple --make-king encoding modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.king_modifier |= kfKingMatrixZs;
	    } else if (!strcmp(cur_modif, "bin")) {
	      if (pc.king_modifier & kfKingMatrixEncodemask) {
		logerrprint("Error: Multiple --make-king encoding modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.king_modifier |= kfKingMatrixBin;
	    } else if (!strcmp(cur_modif, "bin4")) {
	      if (pc.king_modifier & kfKingMatrixEncodemask) {
		logerrprint("Error: Multiple --make-king encoding modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.king_modifier |= kfKingMatrixBin4;
	    } else if (!strcmp(cur_modif, "square")) {
	      if (pc.king_modifier & kfKingMatrixShapemask) {
		logerrprint("Error: Multiple --make-king shape modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.king_modifier |= kfKingMatrixSq;
	    } else if (!strcmp(cur_modif, "square0")) {
	      if (pc.king_modifier & kfKingMatrixShapemask) {
		logerrprint("Error: Multiple --make-king shape modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.king_modifier |= kfKingMatrixSq0;
	    } else if (!strcmp(cur_modif, "triangle")) {
	      if (pc.king_modifier & kfKingMatrixShapemask) {
		logerrprint("Error: Multiple --make-king shape modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.king_modifier |= kfKingMatrixTri;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-king parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (!(pc.king_modifier & kfKingMatrixShapemask)) {
	    if (pc.king_modifier & (kfKingMatrixBin | kfKingMatrixBin4)) {
	      pc.king_modifier |= kfKingMatrixSq;
	    } else {
	      pc.king_modifier |= kfKingMatrixTri;
	    }
	  }
	  pc.command_flags1 |= kfCommand1MakeKing;
	} else if (!memcmp(flagname_p2, "ake-king-table", 15)) {
	  if (king_cutoff_fprefix) {
	    logerrprint("Error: --make-king-table cannot be used with a --king-cutoff input fileset.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 3)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "zs")) {
	      pc.king_modifier |= kfKingTableZs;
	    } else if (!strcmp(cur_modif, "counts")) {
	      pc.king_modifier |= kfKingCounts;
	    } else if ((strlen(cur_modif) > 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      if (pc.king_modifier & kfKingColAll) {
		logerrprint("Error: Multiple --make-king-table cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      reterr = parse_col_descriptor(&(cur_modif[5]), "id\0maybesid\0sid\0nsnp\0hethet\0ibs0\0ibs1\0kinship\0", "make-king-table", kfKingColId, kfKingColDefault, 1, &pc.king_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	      if ((pc.king_modifier & (kfKingColMaybesid | kfKingColSid)) && (!(pc.king_modifier & kfKingColId))) {
		logerrprint("Error: Invalid --make-king-table column set descriptor ('maybesid' and 'sid'\nrequire 'id').\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-king-table parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (!(pc.king_modifier & kfKingColAll)) {
	    pc.king_modifier |= kfKingColDefault;
	  }
	  pc.command_flags1 |= kfCommand1MakeKing;
	} else if (!memcmp(flagname_p2, "issing", 7)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 4)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 2) && (!memcmp(cur_modif, "zs", 2))) {
	      pc.missing_rpt_modifier |= kfMissingRptZs;
	    } else if ((cur_modif_slen == 11) && (!memcmp(cur_modif, "sample-only", 11))) {
	      if (pc.missing_rpt_modifier & kfMissingRptVariantOnly) {
		logerrprint("Error: --missing 'sample-only' and 'variant-only' cannot be used together.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      pc.missing_rpt_modifier |= kfMissingRptSampleOnly;
	    } else if ((cur_modif_slen == 12) && (!memcmp(cur_modif, "variant-only", 12))) {
	      if (pc.missing_rpt_modifier & kfMissingRptSampleOnly) {
		logerrprint("Error: --missing 'sample-only' and 'variant-only' cannot be used together.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      pc.missing_rpt_modifier |= kfMissingRptVariantOnly;
	    } else if ((cur_modif_slen > 6) && (!memcmp(cur_modif, "scols=", 6))) {
	      if (pc.missing_rpt_modifier & kfMissingRptScolAll) {
		logerrprint("Error: Multiple --missing scols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      reterr = parse_col_descriptor(&(cur_modif[6]), "maybesid\0sid\0misspheno1\0missphenos\0nmissdosage\0nmiss\0nmisshh\0hethap\0nobs\0fmissdosage\0fmiss\0fmisshh\0", "missing scols", kfMissingRptScolMaybesid, kfMissingRptScolDefault, 1, &pc.missing_rpt_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else if ((cur_modif_slen > 6) && (!memcmp(cur_modif, "vcols=", 6))) {
	      if (pc.missing_rpt_modifier & kfMissingRptVcolAll) {
		logerrprint("Error: Multiple --missing vcols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      reterr = parse_col_descriptor(&(cur_modif[6]), "chrom\0pos\0ref\0alt1\0alt\0nmissdosage\0nmiss\0nmisshh\0hethap\0nobs\0fmissdosage\0fmiss\0fmisshh\0fhethap\0", "missing vcols", kfMissingRptVcolChrom, kfMissingRptVcolDefault, 1, &pc.missing_rpt_modifier);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --missing parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  const uint32_t explicit_scols = pc.missing_rpt_modifier & kfMissingRptScolAll;
	  if (pc.missing_rpt_modifier & kfMissingRptVariantOnly) {
	    if (explicit_scols) {
	      logerrprint("Error: --missing 'variant-only' and 'scols=' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else {
	    pc.filter_flags |= kfFilterNoSplitChr;
            if (!explicit_scols) {
	      pc.missing_rpt_modifier |= kfMissingRptScolDefault;
	    }
	  }
	  const uint32_t explicit_vcols = pc.missing_rpt_modifier & kfMissingRptVcolAll;
	  if (pc.missing_rpt_modifier & kfMissingRptSampleOnly) {
	    if (explicit_vcols) {
	      logerrprint("Error: --missing 'sample-only' and 'vcols=' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!explicit_vcols) {
	    pc.missing_rpt_modifier |= kfMissingRptVcolDefault;
	  }
	  pc.command_flags1 |= kfCommand1MissingReport;
	} else if (!memcmp(flagname_p2, "aj-ref", 7)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    if (!strcmp(cur_modif, "force")) {
	      pc.misc_flags |= kfMiscMajRefForce;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --maj-ref parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  pc.misc_flags |= kfMiscMajRef;
	  pc.filter_flags |= kfFilterAllReq | kfFilterNoSplitChr;
	} else if (!memcmp(flagname_p2, "af", 3)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    char* cur_modif = argv[arg_idx + 1];
	    if (!scanadv_double(cur_modif, &pc.min_maf)) {
	      sprintf(g_logbuf, "Error: Invalid --maf parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (pc.min_maf < 0.0) {
	      sprintf(g_logbuf, "Error: --maf parameter '%s' too small (must be >= 0).\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    } else if (pc.min_maf >= 1.0) {
	      sprintf(g_logbuf, "Error: --maf parameter '%s' too large (must be < 1).\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  } else {
	    pc.min_maf = 0.01;
	  }
	  if (pc.min_maf != 0.0) {
	    pc.filter_flags |= kfFilterAllReq | kfFilterNoSplitChr;
	  }
	} else if (!memcmp(flagname_p2, "ax-maf", 7)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  if (!scanadv_double(cur_modif, &pc.max_maf)) {
	    sprintf(g_logbuf, "Error: Invalid --max-maf parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (pc.max_maf < pc.min_maf) {
	    sprintf(g_logbuf, "Error: --max-maf parameter '%s' too small (must be >= %g).\n", cur_modif, pc.min_maf);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  } else if (pc.max_maf >= 1.0) {
	    sprintf(g_logbuf, "Error: --max-maf parameter '%s' too large (must be < 1).\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  pc.filter_flags |= kfFilterAllReq | kfFilterNoSplitChr;
	} else if (!memcmp(flagname_p2, "ac", 3)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  double dxx;
	  if ((!scanadv_double(cur_modif, &dxx)) || (dxx < 0.0) || (dxx > 2147483646.0)) {
	    sprintf(g_logbuf, "Error: Invalid --mac parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (dxx > 0.0) {
	    // round up, but keep as much precision as possible
	    int32_t int_part = (int32_t)dxx;
	    dxx -= int_part;
	    pc.min_allele_dosage = int_part * ((uint64_t)kDosageMax);
	    if (dxx > 0.0) {
	      pc.min_allele_dosage += 1 + (dxx * (kDosageMax * (1 - kSmallEpsilon)));
	    }
	    pc.filter_flags |= kfFilterAllReq | kfFilterNoSplitChr;
	  }
	} else if (!memcmp(flagname_p2, "ax-mac", 7)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  double dxx;
	  if ((!scanadv_double(cur_modif, &dxx)) || (dxx < 0.0) || (dxx > 2147483646.0)) {
	    sprintf(g_logbuf, "Error: Invalid --max-mac parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  // round down
	  pc.max_allele_dosage = (int64_t)(dxx * kDosageMax);
	  if (pc.max_allele_dosage < pc.min_allele_dosage) {
	    // yeah, --mac 0.1 --max-mac 0.1 also isn't allowed
	    logerrprint("Error: --max-mac parameter cannot be smaller than --mac parameter.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  pc.filter_flags |= kfFilterAllReq | kfFilterNoSplitChr;
	} else if (!memcmp(flagname_p2, "ind", 4)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t mind_thresh_present = 0;
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "dosage")) {
	      pc.misc_flags |= kfMiscMindDosage;
	    } else if (!strcmp(cur_modif, "hh-missing")) {
	      pc.misc_flags |= kfMiscMindHhMissing;
	    } else if (mind_thresh_present) {
	      logerrprint("Error: Invalid --mind parameter sequence.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (!scanadv_double(cur_modif, &pc.mind_thresh)) {
	      sprintf(g_logbuf, "Error: Invalid --mind parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    } else if ((pc.mind_thresh < 0.0) || (pc.mind_thresh > 1.0)) {
	      sprintf(g_logbuf, "Error: Invalid --mind parameter '%s' (must be in [0, 1]).\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    } else {
	      mind_thresh_present = 1;
	    }
	  }
	  if (!mind_thresh_present) {
	    pc.mind_thresh = 0.1;
	  }
	  if (pc.mind_thresh < 1.0) {
	    pc.filter_flags |= kfFilterAllReq | kfFilterNoSplitChr;
	  }
	} else if (!memcmp(flagname_p2, "issing-var-code", 16)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + 1], argv[arg_idx], kMaxIdSlen, &pc.missing_varid_match);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "erge-par", 9)) {
	  if (pc.exportf_modifier & kfExportfVcf) {
	    logerrprint("Warning: --merge-par should not be used with VCF export.  (The VCF export\nroutine automatically converts PAR1/PAR2 chromosome codes to X, while using\nthe PAR boundaries to get male ploidy right; --merge-par causes VCF export to\nget male ploidy wrong.)\n");
	  }
	  pc.misc_flags |= kfMiscMergePar;
	  pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "af-succ", 8)) {
	  pc.misc_flags |= kfMiscMafSucc;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "ax-corr", 8)) {
	  if (!(pc.command_flags1 & kfCommand1Glm)) {
	    logerrprint("Error: --max-corr must be used with --glm.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  if (!scanadv_double(cur_modif, &pc.glm_info.max_corr)) {
	    sprintf(g_logbuf, "Error: Invalid --max-corr parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if ((pc.glm_info.max_corr < 0.0) || (pc.glm_info.max_corr > 1.0)) {
	    sprintf(g_logbuf, "Error: Invalid --max-corr parameter '%s' (must be in [0, 1]).\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else if (!memcmp(flagname_p2, "ach-r2-filter", 14)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    char* cur_modif = argv[arg_idx + 1];
	    if (!scanadv_double(cur_modif, &pc.mach_r2_min)) {
	      sprintf(g_logbuf, "Error: Invalid --mach-r2-filter min parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (pc.mach_r2_min < 0.0) {
	      sprintf(g_logbuf, "Error: Invalid --mach-r2-filter min parameter '%s' (must be nonnegative).\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (param_ct == 2) {
	      cur_modif = argv[arg_idx + 2];
	      if (!scanadv_double(cur_modif, &pc.mach_r2_max)) {
		sprintf(g_logbuf, "Error: Invalid --mach-r2-filter max parameter '%s'.\n", cur_modif);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	    } else {
	      pc.mach_r2_max = 2.0;
	    }
	    if (pc.mach_r2_max < pc.mach_r2_min) {
	      logerrprint("Error: --mach-r2-filter min parameter cannot be larger than max parameter.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else {
	    pc.mach_r2_min = 0.1;
	  }
	  pc.filter_flags |= kfFilterAllReq | kfFilterNoSplitChr;
	} else if (!memcmp(flagname_p2, "issing-code", 12)) {
	  if (!(xload & (kfXloadOxGen | kfXloadOxBgen))) {
	    // could technically support pure .sample -> .fam/.psam, but let's
	    // keep this simple
	    logerrprint("Error: --missing-code must be used with --data/--gen/--bgen.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(param_ct? argv[arg_idx + 1] : "", argv[arg_idx], 0x7fffffff, &ox_missing_code);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "issing-genotype", 16)) {
	  logerrprint("Error: --missing-genotype flag retired.  Use --input-missing-genotype and/or\n--output-missing-genotype.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!memcmp(flagname_p2, "issing-phenotype", 17)) {
	  logerrprint("Error: --missing-phenotype flag retired.  Use --input-missing-phenotype and/or\n--output-missing-phenotype.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!memcmp(flagname_p2, "issing-catname", 15)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  const uint32_t cur_modif_slen = strlen(cur_modif);
	  double dxx;
	  if (scanadv_double(cur_modif, &dxx) || is_nan_str(cur_modif, cur_modif_slen)) {
	    logerrprint("Error: --missing-catname string cannot be 'NA' or start with a number.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (cur_modif_slen > 31) {
	    logerrprint("Error: --missing-catname string too long (max 31 chars).\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  memcpy(g_missing_catname, cur_modif, cur_modif_slen + 1);
	} else if (!memcmp(flagname_p2, "ouse", 5)) {
	  if (chr_info.chrset_source) {
	    logerrprint("Error: Conflicting chromosome-set flags.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  chr_info.chrset_source = kChrsetSourceCmdline;
	  chr_info.autosome_ct = 19;
	  chr_info.xymt_codes[0] = 20;
	  chr_info.xymt_codes[1] = 21;
	  chr_info.xymt_codes[2] = -2;
	  chr_info.xymt_codes[3] = -2;
	  chr_info.xymt_codes[4] = -2;
	  chr_info.xymt_codes[5] = -2;
	  chr_info.haploid_mask[0] = 0x300000;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "ake-grm", 8)) {
	  logerrprint("Error: --make-grm has been retired due to inconsistent meaning across GCTA\nversions.  Use --make-grm-gz or --make-grm-bin.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!memcmp(flagname_p2, "ake-grm-bin", 12)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "cov")) {
	      pc.grm_flags |= kfGrmCov;
	    } else if (!strcmp(cur_modif, "meanimpute")) {
	      pc.grm_flags |= kfGrmMeanimpute;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-grm-bin parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  pc.grm_flags |= kfGrmBin;
	  pc.command_flags1 |= kfCommand1MakeRel;
	} else if (!memcmp(flagname_p2, "ake-grm-gz", 11)) {
	  if (pc.command_flags1 & kfCommand1MakeRel) {
	    logerrprint("Error: --make-grm-gz cannot be used with --make-grm-bin.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 3)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t compress_stream_type = 0; // 1 = no-gz, 2 = zs
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "cov")) {
	      pc.grm_flags |= kfGrmCov;
	    } else if (!strcmp(cur_modif, "meanimpute")) {
	      pc.grm_flags |= kfGrmMeanimpute;
	    } else if (!strcmp(cur_modif, "no-gz")) {
	      if (compress_stream_type) {
		logerrprint("Error: Multiple --make-grm-gz compression type modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      compress_stream_type = 1;
	      pc.grm_flags |= kfGrmTableNoGz;
	    } else if (!strcmp(cur_modif, "zs")) {
	      if (compress_stream_type) {
		logerrprint("Error: Multiple --make-grm-gz compression type modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      compress_stream_type = 2;
	      pc.grm_flags |= kfGrmTableZs;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-grm-gz parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if (!compress_stream_type) {
	    pc.grm_flags |= kfGrmTableGz;
	  }
	  pc.command_flags1 |= kfCommand1MakeRel;
	} else if (!memcmp(flagname_p2, "ake-rel", 8)) {
	  if (pc.command_flags1 & kfCommand1MakeRel) {
	    logerrprint("Error: --make-rel cannot be used with --make-grm-gz/--make-grm-bin.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 4)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (!strcmp(cur_modif, "cov")) {
	      pc.grm_flags |= kfGrmCov;
	    } else if (!strcmp(cur_modif, "meanimpute")) {
	      pc.grm_flags |= kfGrmMeanimpute;
	    } else if (!strcmp(cur_modif, "zs")) {
	      if (pc.grm_flags & kfGrmMatrixEncodemask) {
		logerrprint("Error: Multiple --make-rel encoding modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.grm_flags |= kfGrmMatrixZs;
	    } else if (!strcmp(cur_modif, "bin")) {
	      if (pc.grm_flags & kfGrmMatrixEncodemask) {
		logerrprint("Error: Multiple --make-rel encoding modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.grm_flags |= kfGrmMatrixBin;
	    } else if (!strcmp(cur_modif, "bin4")) {
	      if (pc.grm_flags & kfGrmMatrixEncodemask) {
		logerrprint("Error: Multiple --make-rel encoding modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.grm_flags |= kfGrmMatrixBin4;
	    } else if (!strcmp(cur_modif, "square")) {
	      if (pc.grm_flags & kfGrmMatrixShapemask) {
		logerrprint("Error: Multiple --make-rel shape modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.grm_flags |= kfGrmMatrixSq;
	    } else if (!strcmp(cur_modif, "square0")) {
	      if (pc.grm_flags & kfGrmMatrixShapemask) {
		logerrprint("Error: Multiple --make-rel shape modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.grm_flags |= kfGrmMatrixSq0;
	    } else if (!strcmp(cur_modif, "triangle")) {
	      if (pc.grm_flags & kfGrmMatrixShapemask) {
		logerrprint("Error: Multiple --make-rel shape modifiers.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      pc.grm_flags |= kfGrmMatrixTri;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --make-rel parameter '%s'.\n", cur_modif);
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
	} else if (!memcmp(flagname_p2, "ap", 3)) {
	  if (load_params || (xload & (~kfXloadPlink1Dosage))) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  const uint32_t slen = strlen(cur_modif);
	  if (slen > kPglFnamesize - 1) {
	    logerrprint("Error: --map filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pvarname, cur_modif, slen + 1);
	  xload |= kfXloadMap;
	} else if (!memcmp(flagname_p2, "within", 7)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  if (scan_posint_capped(cur_modif, kMaxLongLine / 2, &pc.mwithin_val)) {
	    sprintf(g_logbuf, "Error: Invalid --mwithin parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;

      case 'n':
	if (!memcmp(flagname_p2, "o-fid", 6)) {
	  pc.fam_cols &= ~kfFamCol1;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "o-parents", 10)) {
	  pc.fam_cols &= ~kfFamCol34;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "o-sex", 6)) {
	  pc.fam_cols &= ~kfFamCol5;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "o-pheno", 8)) {
	  pc.fam_cols &= ~kfFamCol6;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "onfounders", 11)) {
	  pc.misc_flags |= kfMiscNonfounders;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "ot-chr", 7)) {
	  if (pc.varid_from) {
	    logerrprint("Error: --from/--to cannot be used with --autosome{-par} or --{not-}chr.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (pc.from_bp != -1) {
	    logerrprint("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb cannot be used with --not-chr.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }

	  // allowed:
	  //   --allow-extra-chr --chr 5-22 bobs_chrom --not-chr 17
	  // allowed:
	  //   --allow-extra-chr --not-chr 12-17 bobs_chrom
	  // does not make sense, disallowed:
	  //   --allow-extra-chr --chr 5-22 --not-chr bobs_chrom
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }

	  // --allow-extra-chr present, --chr/--autosome{-xy} not present
	  const uint32_t aec_and_no_chr_include = ((pc.misc_flags / kfMiscAllowExtraChrs) & 1) && (!chr_info.is_include_stack);
	  reterr = parse_chr_ranges(flagname_p, errstr_append, param_ct, aec_and_no_chr_include, kChrRawEnd - (kChrExcludeWords * kBitsPerWord), '-', &(argv[arg_idx]), &chr_info, chr_info.chr_exclude);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  notchr_present = 1;
	  // remaining processing now postponed to finalize_chrset()
	} else if (!memcmp(flagname_p2, "ew-id-max-allele-len", 21)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  if (scan_posint_capped(cur_modif, kMaxIdSlen - 2, &pc.new_variant_id_max_allele_slen)) {
	    sprintf(g_logbuf, "Error: Invalid --new-id-max-allele-len length parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (param_ct == 2) {
	    cur_modif = argv[arg_idx + 2];
	    if (!strcmp(cur_modif, "missing")) {
	      pc.misc_flags |= kfMiscNewVarIdOverflowMissing;
	    } else if (!strcmp(cur_modif, "truncate")) {
	      pc.misc_flags |= kfMiscNewVarIdOverflowTruncate;
	    } else if (strcmp(cur_modif, "error")) {
	      sprintf(g_logbuf, "Error: Invalid --new-id-max-allele-len parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'o':
	if (!memcmp(flagname_p2, "utput-chr", 10)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* mt_code = argv[arg_idx + 1];
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
	    sprintf(g_logbuf, "Error: Invalid --output-chr parameter '%s'.\n", mt_code);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else if (!memcmp(flagname_p2, "utput-min-p", 12)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  if ((!scanadv_double(cur_modif, &pc.output_min_p)) || (!(pc.output_min_p >= 0.0)) || (pc.output_min_p >= 1.0)) {
	    sprintf(g_logbuf, "Error: Invalid --output-min-p parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else if (!memcmp(flagname_p2, "xford-single-chr", 17)) {
	  if (!(xload & kfXloadOxGen)) {
	    if (xload & kfXloadOxBgen) {
	      logerrprint("Error: --oxford-single-chr must be used with .gen input.  (Single-chromosome\n.bgen files do not require this, since they still contain chromosome codes.)\n");
	    } else {
	      logerrprint("Error: --oxford-single-chr must be used with .gen input.\n");
	    }
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  if (!(pc.misc_flags & kfMiscAllowExtraChrs)) {
	    if (get_chr_code_raw(cur_modif) < 0) {
	      sprintf(g_logbuf, "Error: Invalid --oxford-single-chr chromosome code '%s'. (Did you forget --allow-extra-chr?)\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  reterr = cmdline_alloc_string(cur_modif, argv[arg_idx], kMaxIdSlen, &import_single_chr_str);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "utput-missing-genotype", 23)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  output_missing_geno_char = extract_char_param(cur_modif);
	  if (((unsigned char)output_missing_geno_char) <= ' ') {
	    sprintf(g_logbuf, "Error: Invalid --output-missing-genotype parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else if (!memcmp(flagname_p2, "utput-missing-phenotype", 24)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  const uint32_t cur_modif_slen = strlen(cur_modif);
	  if (cur_modif_slen > 31) {
	    logerrprint("Error: --output-missing-phenotype string too long (max 31 chars).\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  memcpy(g_output_missing_pheno, cur_modif, cur_modif_slen + 1);
	} else if (memcmp(flagname_p2, "ut", 3)) {
	  // --out is a special case due to logging
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'p':
	if (!memcmp(flagname_p2, "file", 5)) {
	  if (load_params || xload) {
	    // currently only possible with --bcf, --bfile, --pfile
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t fname_modif_idx = 1;
	  if (param_ct == 2) {
	    if (check_extra_param(&(argv[arg_idx]), "vzs", &fname_modif_idx)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	  const char* fname_prefix = argv[arg_idx + fname_modif_idx];
	  const uint32_t slen = strlen(fname_prefix);
	  if (slen > (kPglFnamesize - 10)) {
	    logerrprint("Error: --pfile parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(memcpya(pgenname, fname_prefix, slen), ".pgen");
	  strcpy(memcpya(psamname, fname_prefix, slen), ".psam");
	  char* pvarname_end = memcpya(pvarname, fname_prefix, slen);
	  pvarname_end = strcpya0(pvarname_end, ".pvar");
	  if (param_ct == 2) {
	    strcpy(pvarname_end, ".zst");
	  }
	  load_params |= kfLoadParamsPfileAll;
	} else if (!memcmp(flagname_p2, "gen", 4)) {
	  if (xload) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  load_params |= kfLoadParamsPgen;
	  char* fname = argv[arg_idx + 1];
	  const uint32_t slen = strlen(fname);
	  if (slen > (kPglFnamesize - 1)) {
	    logerrprint("Error: --pgen parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pgenname, fname, slen + 1);
	} else if (!memcmp(flagname_p2, "sam", 4)) {
	  if (xload & (~(kfXloadVcf | kfXloadBcf | kfXloadPlink1Dosage | kfXloadMap))) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  load_params |= kfLoadParamsPsam;
	  char* fname = argv[arg_idx + 1];
	  const uint32_t slen = strlen(fname);
	  if (slen > (kPglFnamesize - 1)) {
	    logerrprint("Error: --psam parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(psamname, fname, slen + 1);
	} else if (!memcmp(flagname_p2, "var", 4)) {
	  if (xload) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  load_params |= kfLoadParamsPvar;
	  char* fname = argv[arg_idx + 1];
	  const uint32_t slen = strlen(fname);
	  if (slen > (kPglFnamesize - 1)) {
	    logerrprint("Error: --pvar parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pvarname, fname, slen + 1);
	} else if (!memcmp(flagname_p2, "heno", 5)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_fname(argv[arg_idx + 1], flagname_p, 0, &pc.pheno_fname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "heno-name", 10)) {
	  // can now be used without --pheno
	  reterr = parse_name_ranges(&(argv[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.pheno_range_list);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "arallel", 8)) {
	  if (pc.king_modifier & kfKingMatrixSq) {
	    logerrprint("Error: --parallel cannot be used with '--make-king square'.  Use '--make-king\nsquare0' or plain --make-king instead.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if ((pc.king_cutoff != -1) && (!king_cutoff_fprefix)) {
	    logerrprint("Error: --parallel cannot be used with --king-cutoff.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (pc.grm_flags & kfGrmMatrixSq) {
	    logerrprint("Error: --parallel cannot be used with '--make-rel square'.  Use '--make-rel\nsquare0' or plain --make-rel instead.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 2, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (scan_posint_capped(argv[arg_idx + 1], kParallelMax, &pc.parallel_idx)) {
	    sprintf(g_logbuf, "Error: Invalid --parallel job index '%s'.\n", argv[arg_idx + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (scan_posint_capped(argv[arg_idx + 2], kParallelMax, &pc.parallel_tot) || (pc.parallel_tot == 1) || (pc.parallel_tot < pc.parallel_idx)) {
	    sprintf(g_logbuf, "Error: Invalid --parallel total job count '%s'.\n", argv[arg_idx + 2]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  --pc.parallel_idx; // internal 0..(n-1) indexing
	} else if (!memcmp(flagname_p2, "arameters", 10)) {
	  if (!(pc.command_flags1 & kfCommand1Glm)) {
	    logerrprint("Error: --parameters must be used with --glm.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  reterr = parse_name_ranges(&(argv[arg_idx]), errstr_append, param_ct, 1, '-', &pc.glm_info.parameters_range_list);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "filter", 7)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  if (!scanadv_double(cur_modif, &pc.pfilter)) {
	    sprintf(g_logbuf, "Error: Invalid --pfilter parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if ((pc.pfilter <= 0.0) || (pc.pfilter > 1.0)) {
	    logerrprint("Error: --pfilter threshold must be in (0, 1].\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	} else if (!memcmp(flagname_p2, "ca", 3)) {
#ifdef NOLAPACK
	  logerrprint("Error: --pca requires " PROG_NAME_STR " to be built with LAPACK.\n");
	  goto main_ret_INVALID_CMDLINE;
#endif
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 6)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  uint32_t is_var_wts = 0;
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 6) && (!memcmp(cur_modif, "approx", 6))) {
	      pc.pca_flags |= kfPcaApprox;
	    } else if ((cur_modif_slen == 10) && (!memcmp(cur_modif, "meanimpute", 10))) {
	      pc.pca_flags |= kfPcaMeanimpute;
	    } else if ((cur_modif_slen == 3) && (!memcmp(cur_modif, "sid", 3))) {
	      pc.pca_flags |= kfPcaSid;
	    } else if ((cur_modif_slen == 7) && (!memcmp(cur_modif, "var-wts", 7))) {
	      is_var_wts = 1;
	    } else if ((cur_modif_slen == 3) && (!memcmp(cur_modif, "vzs", 3))) {
	      pc.pca_flags |= kfPcaVarZs;
	    } else if ((cur_modif_slen > 6) && (!memcmp(cur_modif, "vcols=", 6))) {
	      if (pc.pca_flags & kfPcaVcolAll) {
		logerrprint("Error: Multiple --pca vcols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      reterr = parse_col_descriptor(&(cur_modif[6]), "chrom\0pos\0ref\0alt1\0alt\0maj\0nonmaj\0", "pca vcols", kfPcaVcolChrom, kfPcaVcolDefault, 1, &pc.pca_flags);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else {
	      if (pc.pca_ct || scan_posint_defcap(cur_modif, &pc.pca_ct)) {
		logerrprint("Error: Invalid --pca parameter sequence.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      if (pc.pca_ct > 8000) {
		// this slightly simplifies output buffering.
		// lower limit for randomized algorithm?
		// (just let memory allocation fail for now...)
		logerrprint("Error: --pca does not support more than 8000 PCs.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	    }
	  }
	  if (pc.pca_flags & kfPcaApprox) {
	    if (pc.pca_ct > 100) {
	      // double-precision overflow too likely
	      logerrprint("Error: --pca approx does not support more than 100 PCs.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	  } else {
	    // todo: if --make-rel/--make-grm present, verify consistency
	    if (pc.parallel_tot != 1) {
	      logerrprint("Error: Non-approximate --pca cannot be used with --parallel.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    const uint32_t pca_meanimpute = (pc.pca_flags / kfPcaMeanimpute) & 1;
	    if (pc.command_flags1 & kfCommand1MakeRel) {
	      if (((pc.grm_flags / kfGrmMeanimpute) & 1) != pca_meanimpute) {
		logerrprint("Error: --make-rel/--make-grm-gz/--make-grm-bin meanimpute setting must match\n--pca meanimpute setting.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      if (pc.grm_flags & kfGrmCov) {
		logerrprint("Error: --make-rel/--make-grm-gz/--make-grm-bin cannot be used to compute a\ncovariance matrix in the same run as non-approximate --pca.\n");
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
	  if (!(pc.pca_flags & kfPcaVcolAll)) {
	    if (is_var_wts) {
	      pc.pca_flags |= kfPcaVcolDefault;
	    }
	  } else if (!is_var_wts) {
	    logerrprint("Error: --pca 'vcols=' has no effect without 'var-wts'.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (is_var_wts) {
	    pc.pca_flags |= kfPcaVarWts;
	  } else if (pc.pca_flags & kfPcaVarZs) {
	    logerrprint("Error: --pca 'vzs' modifier has no effect without 'var-wts'.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  pc.command_flags1 |= kfCommand1Pca;
	} else if (!memcmp(flagname_p2, "heno-quantile-normalize", 24)) {
	  if (param_ct) {
	    reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, 0x7fffffff, &pc.quantnorm_flattened);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.pheno_transform_flags |= kfPhenoTransformQuantnormPheno;
	  pc.filter_flags |= kfFilterPsamReq;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;

      case 'q':
	if (!memcmp(flagname_p2, "uantile-normalize", 18)) {
	  if (pc.pheno_transform_flags & (kfPhenoTransformQuantnormPheno | kfPhenoTransformQuantnormCovar)) {
	    logerrprint("Error: --quantile-normalize cannot be used with --pheno-quantile-normalize or\n--covar-quantile-normalize.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (param_ct) {
	    reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, 0x7fffffff, &pc.quantnorm_flattened);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.pheno_transform_flags |= kfPhenoTransformQuantnormAll;
	  pc.filter_flags |= kfFilterPsamReq;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'r':
	if (!memcmp(flagname_p2, "eal-ref-alleles", 16)) {
	  if (pc.misc_flags & kfMiscMajRef) {
	    logerrprint("Error: --real-ref-alleles cannot be used with --maj-ref.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  pc.misc_flags |= kfMiscRealRefAlleles;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "emove", 6)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const uint32_t sid_present = !strcmp(argv[arg_idx + 1], "sid");
	  if (sid_present) {
	    if (param_ct == 1) {
	      logerrprint("Error: '--remove sid' requires at least one filename.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    pc.misc_flags |= kfMiscRemovefileSid;
	  }
	  reterr = alloc_and_flatten(&(argv[arg_idx + 1 + sid_present]), param_ct - sid_present, kPglFnamesize, &pc.remove_fnames);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "emove-fam", 10)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, kPglFnamesize, &pc.removefam_fnames);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "emove-females", 14)) {
	  pc.filter_flags |= kfFilterPsamReq | kfFilterExclFemales;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "emove-males", 12)) {
	  pc.filter_flags |= kfFilterPsamReq | kfFilterExclMales;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "emove-nosex", 12)) {
	  pc.filter_flags |= kfFilterPsamReq | kfFilterExclNosex;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "ead-freq", 9)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_fname(argv[arg_idx + 1], flagname_p, 0, &pc.read_freq_fname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterAllReq;
	} else if (!memcmp(flagname_p2, "equire-pheno", 13)) {
	  if (param_ct) {
	    reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, 0x7fffffff, &require_pheno_flattened);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.misc_flags |= kfMiscRequirePheno;
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "equire-covar", 13)) {
	  if (param_ct) {
	    reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, 0x7fffffff, &require_covar_flattened);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.misc_flags |= kfMiscRequireCovar;
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "emove-if", 9)) {
	  reterr = validate_and_alloc_cmp_expr(&(argv[arg_idx + 1]), argv[arg_idx], param_ct, &pc.remove_if_expr);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "emove-cats", 11)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_fname(argv[arg_idx + 1], flagname_p, 0, &pc.remove_cats_fname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "emove-cat-names", 16)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, kMaxIdBlen, &pc.remove_cat_names_flattened);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "emove-cat-pheno", 14)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + 1], argv[arg_idx], kMaxIdSlen, &pc.remove_cat_phenoname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "ice", 4)) {
	  if (chr_info.chrset_source) {
	    logerrprint("Error: Conflicting chromosome-set flags.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  chr_info.chrset_source = kChrsetSourceCmdline;
	  chr_info.autosome_ct = 12;
	  chr_info.xymt_codes[0] = -2;
	  chr_info.xymt_codes[1] = -2;
	  chr_info.xymt_codes[2] = -2;
	  chr_info.xymt_codes[3] = -2;
	  chr_info.xymt_codes[4] = -2;
	  chr_info.xymt_codes[5] = -2;
	  chr_info.haploid_mask[0] = 0x1fff;
	  goto main_param_zero;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 's':
	if (!memcmp(flagname_p2, "eed", 4)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 0x7fffffff)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  rseed_ct = param_ct;
	  if (pgl_malloc(param_ct * sizeof(int32_t), &rseeds)) {
	    goto main_ret_NOMEM;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    if (scan_uint_capped(cur_modif, 0xffffffffU, &(rseeds[param_idx - 1]))) {
	      sprintf(g_logbuf, "Error: Invalid --seed parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	} else if (!memcmp(flagname_p2, "plit-par", 9)) {
	  if (pc.misc_flags & kfMiscMergePar) {
	    logerrprint("Error: --split-par cannot be used with --merge-par.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct == 1) {
	    const char* build_code = argv[arg_idx + 1];
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
	      sprintf(g_logbuf, "Error: Unrecognized --split-par build code '%s'.\n", build_code);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  } else {
	    if (scan_uint_defcap(argv[arg_idx + 1], &pc.splitpar_bound1)) {
	      sprintf(g_logbuf, "Error: Invalid --split-par parameter '%s'.\n", argv[arg_idx + 1]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (scan_uint_defcap(argv[arg_idx + 2], &pc.splitpar_bound2) || (pc.splitpar_bound2 <= pc.splitpar_bound1)) {
	      sprintf(g_logbuf, "Error: Invalid --split-par parameter '%s'.\n", argv[arg_idx + 2]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
	} else if ((!memcmp(flagname_p2, "et-all-var-ids", 15)) || (!memcmp(flagname_p2, "et-missing-var-ids", 19))) {
	  if (flagname_p2[3] == 'm') {
	    if (pc.varid_template) {
	      logerrprint("Error: --set-missing-var-ids cannot be used with --set-all-var-ids.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    pc.misc_flags |= kfMiscSetMissingVarIds;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (!varid_template_is_valid(argv[arg_idx + 1], flagname_p)) {
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + 1], argv[arg_idx], kMaxIdSlen, &pc.varid_template);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPvarReq;
	} else if (!memcmp(flagname_p2, "et-hh-missing", 14)) {
	  if (!(pc.command_flags1 & kfCommand1MakePlink2)) {
	    logerrprint("Error: --set-hh-missing must be used with --make-{b}pgen/--make-bed.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  make_plink2_modifier |= kfMakePlink2SetHhMissing;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "et-mixed-mt-missing", 20)) {
	  if (!(pc.command_flags1 & kfCommand1MakePlink2)) {
	    logerrprint("Error: --set-mixed-mt-missing must be used with --make-{b}pgen/--make-bed.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  make_plink2_modifier |= kfMakePlink2SetMixedMtMissing;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "ample", 6)) {
	  if (load_params || (xload & (~(kfXloadOxGen | kfXloadOxBgen | kfXloadOxHaps | kfXloadOxLegend | kfXloadOxSample)))) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (!(xload & (kfXloadOxGen | kfXloadOxBgen | kfXloadOxHaps))) {
	    logerrprint("Error: --sample must be used with --gen/--bgen/--data/--haps.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_fname = argv[arg_idx + 1];
	  const uint32_t slen = strlen(cur_fname);
	  if (slen > kPglFnamesize - 1) {
	    logerrprint("Error: --sample filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(psamname, cur_fname, slen + 1);
	  xload |= kfXloadOxSample;
	} else if (!memcmp(flagname_p2, "heep", 5)) {
	  if (chr_info.chrset_source) {
	    logerrprint("Error: Conflicting chromosome-set flags.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  chr_info.chrset_source = kChrsetSourceCmdline;
	  chr_info.autosome_ct = 26;
	  chr_info.xymt_codes[0] = 27;
	  chr_info.xymt_codes[1] = 28;
	  chr_info.xymt_codes[2] = -2;
	  chr_info.xymt_codes[3] = -2;
	  chr_info.xymt_codes[4] = -2;
	  chr_info.xymt_codes[5] = -2;
	  chr_info.haploid_mask[0] = 0x18000000;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "np", 3)) {
	  if (pc.varid_exclude_snp) {
	    // problematic due to --window
	    logerrprint("Error: --snp cannot be used with --exclude-snp.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + 1], argv[arg_idx], kMaxIdSlen, &pc.varid_snp);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPvarReq;
	} else if (!memcmp(flagname_p2, "nps", 4)) {
	  reterr = parse_name_ranges(&(argv[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.snps_range_list);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPvarReq;
	} else if (!memcmp(flagname_p2, "nps-only", 9)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    if (!strcmp(cur_modif, "just-acgt")) {
	      pc.filter_flags |= kfFilterSnpsOnlyJustAcgt;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --snps-only parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  pc.filter_flags |= kfFilterPvarReq | kfFilterSnpsOnly;
	} else if (!memcmp(flagname_p2, "core", 5)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 11)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_fname(argv[arg_idx + 1], flagname_p, 0, &pc.score_info.input_fname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  uint32_t numeric_param_ct = 0;
	  uint32_t score_cols[3];
	  for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 6) && (!memcmp(cur_modif, "header", 6))) {
	      pc.score_info.flags |= kfScoreHeaderIgnore;
	    } else if ((cur_modif_slen == 11) && (!memcmp(cur_modif, "header-read", 11))) {
	      pc.score_info.flags |= kfScoreHeaderRead;
	    } else if ((cur_modif_slen == 18) && (!memcmp(cur_modif, "no-mean-imputation", 18))) {
	      pc.score_info.flags |= kfScoreNoMeanimpute;
	    } else if ((cur_modif_slen == 6) && (!memcmp(cur_modif, "center", 6))) {
	      pc.score_info.flags |= kfScoreCenter;
	    } else if ((cur_modif_slen == 20) && (!memcmp(cur_modif, "variance-standardize", 20))) {
	      pc.score_info.flags |= kfScoreVarianceStandardize;
	    } else if ((cur_modif_slen == 18) && (!memcmp(cur_modif, "variance-normalize", 18))) {
	      logerrprint("Note: --score's 'variance-normalize' modifier has been renamed to the more\nprecise 'variance-standardize'.\n");
	      pc.score_info.flags |= kfScoreVarianceStandardize;
	    } else if ((cur_modif_slen == 2) && (!memcmp(cur_modif, "se", 2))) {
	      pc.score_info.flags |= kfScoreSe;
	    } else if ((cur_modif_slen == 2) && (!memcmp(cur_modif, "zs", 2))) {
	      pc.score_info.flags |= kfScoreZs;
	    } else if ((cur_modif_slen == 13) && (!memcmp(cur_modif, "list-variants", 13))) {
	      pc.score_info.flags |= kfScoreListVariants;
	    } else if ((cur_modif_slen == 16) && (!memcmp(cur_modif, "list-variants-zs", 16))) {
	      pc.score_info.flags |= kfScoreListVariants | kfScoreListVariantsZs;
	    } else if ((cur_modif_slen > 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      if (pc.score_info.flags & kfScoreColAll) {
		logerrprint("Error: Multiple --score cols= modifiers.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      reterr = parse_col_descriptor(&(cur_modif[5]), "maybesid\0sid\0pheno1\0phenos\0nmissallele\0denom\0dosagesum\0scoreavgs\0scoresums\0", "score", kfScoreColMaybesid, kfScoreColDefault, 1, &pc.score_info.flags);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else if (numeric_param_ct == 3) {
	      logerrprint("Error: --score takes at most three numeric parameters.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else {
	      if (scan_posint_capped(cur_modif, kMaxLongLine / 2, &(score_cols[numeric_param_ct]))) {
		sprintf(g_logbuf, "Error: Invalid --score parameter '%s'.\n", cur_modif);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      for (uint32_t uii = 0; uii < numeric_param_ct; ++uii) {
		if (score_cols[uii] == score_cols[numeric_param_ct]) {
		  logerrprint("Error: Identical --score column indexes.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
	      }
	      ++numeric_param_ct;
	    }
	  }
	  if ((pc.score_info.flags & (kfScoreHeaderIgnore | kfScoreHeaderRead)) == (kfScoreHeaderIgnore | kfScoreHeaderRead)) {
	    logerrprint("Error: --score 'header' and 'header-read' modifiers cannot be used together.\n");
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
	    const uint32_t col_idx_blen = 1 + int_slen(col_idx);
	    char* new_buf;
	    if (pgl_malloc(col_idx_blen + 1, &new_buf)) {
	      goto main_ret_NOMEM;
	    }
	    pc.score_info.input_col_idx_range_list.names = new_buf;
	    pc.score_info.input_col_idx_range_list.name_max_blen = col_idx_blen;
	    pc.score_info.input_col_idx_range_list.name_ct = 1;
	    uint32toa_x(col_idx, '\0', new_buf);
	    new_buf[col_idx_blen] = '\0';
	    pc.score_info.input_col_idx_range_list.starts_range = (unsigned char*)(&(new_buf[col_idx_blen]));
	  }
	  pc.command_flags1 |= kfCommand1Score;
	} else if (!memcmp(flagname_p2, "core-col-nums", 14)) {
	  if (!(pc.command_flags1 & kfCommand1Score)) {
	    logerrprint("Error: --score-col-nums must be used with --score.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (pc.score_info.input_col_idx_range_list.name_ct) {
	    logerrprint("Error: --score-col-nums cannot be used when three numeric parameters are\nprovided to --score.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  reterr = parse_name_ranges(&(argv[arg_idx]), errstr_append, param_ct, 1, '-', &pc.score_info.input_col_idx_range_list);
	  if (reterr) {
	    goto main_ret_1;
	  }
	} else if (!memcmp(flagname_p2, "plit-cat-pheno", 15)) {
	  uint32_t first_phenoname_idx = 1;
	  for (; first_phenoname_idx <= param_ct; ++first_phenoname_idx) {
	    const char* cur_modif = argv[arg_idx + first_phenoname_idx];
	    if (!strcmp(cur_modif, "omit-last")) {
	      pc.pheno_transform_flags |= kfPhenoTransformSplitCatOmitLast;
	    } else if (!strcmp(cur_modif, "covar-01")) {
	      pc.pheno_transform_flags |= kfPhenoTransformSplitCatCovar01;
	    } else {
	      break;
	    }
	  }
	  if (first_phenoname_idx <= param_ct) {
	    reterr = alloc_and_flatten(&(argv[arg_idx + first_phenoname_idx]), param_ct + 1 - first_phenoname_idx, kMaxIdSlen - 1, &pc.split_cat_phenonames_flattened);
	    if (reterr) {
	      goto main_ret_1;
	    }
	    // may as well verify that no phenotype name has an '=' in it
	    char* phenonames_iter = pc.split_cat_phenonames_flattened;
	    do {
	      const uint32_t cur_phenoname_slen = strlen(phenonames_iter);
	      if (memchr(phenonames_iter, '=', cur_phenoname_slen)) {
		logerrprint("Error: --split-cat-pheno phenotype names may not contain the '=' character.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      phenonames_iter = &(phenonames_iter[cur_phenoname_slen + 1]);
	    } while (*phenonames_iter);
	  } else if (pc.pheno_transform_flags & kfPhenoTransformSplitCatCovar01) {
	    logerrprint("Error: --split-cat-pheno 'covar-01' modifier cannot be used without any\nphenotype names.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  pc.pheno_transform_flags |= kfPhenoTransformSplitCat;
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "ort-vars", 9)) {
	  logerrprint("Error: --sort-vars is not implemented yet.\n");
	  reterr = kPglRetNotYetSupported;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 't':
	if (!memcmp(flagname_p2, "hreads", 7)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (scan_posint_defcap(argv[arg_idx + 1], &pc.max_thread_ct)) {
	    sprintf(g_logbuf, "Error: Invalid --threads parameter '%s'.\n", argv[arg_idx + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (pc.max_thread_ct > kMaxThreads) {
	    LOGPRINTF("Note: Reducing --threads parameter to %u.  (If this is not large enough,\nrecompile with a larger kMaxThreads setting.)\n", kMaxThreads);
	    pc.max_thread_ct = kMaxThreads;
	  } else if (known_procs == -1) {
	    // trigger BLAS/LAPACK warning?
	    known_procs = 0;
	  }
	} else if (!memcmp(flagname_p2, "o", 2)) {
	  if (chr_info.is_include_stack || notchr_present) {
	    logerrprint("Error: --from/--to cannot be used with --autosome{-par} or --{not-}chr.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = cmdline_alloc_string(argv[arg_idx + 1], argv[arg_idx], kMaxIdSlen, &pc.varid_to);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
	} else if ((!memcmp(flagname_p2, "o-bp", 5)) || (!memcmp(flagname_p2, "o-kb", 5)) || (!memcmp(flagname_p2, "o-mb", 5))) {
	  if (!cmdline_single_chr(&chr_info, pc.misc_flags)) {
	    logerrprint("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb must be used with --chr, and only\none chromosome.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (notchr_present) {
	    logerrprint("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb cannot be used with --not-chr.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (pc.to_bp != -1) {
	    logerrprint("Error: Multiple --to-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  double dxx;
	  if (!scanadv_double(cur_modif, &dxx)) {
	    sprintf(g_logbuf, "Error: Invalid --to-bp/-kb/-mb parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  const char unit_char = flagname_p2[2];
	  if (unit_char == 'k') {
	    dxx *= 1000;
	  } else if (unit_char == 'm') {
	    dxx *= 1000000;
	  }
	  if (dxx < 0) {
	    LOGERRPRINTF("Error: --to-bp/-kb/-mb parameter '%s' too small.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_A;
	  } else if (dxx >= 2147483646) {
	    pc.to_bp = 0x7ffffffe;
	  } else {
	    // round down
	    pc.to_bp = (int32_t)(dxx * (1 + kSmallEpsilon));
	  }
	  if (pc.from_bp > pc.to_bp) {
	    // (if we do permit this, rounding must be postponed)
	    logerrprint("Error: --to-bp/-kb/-mb parameter is smaller than --from-bp/-kb/-mb parameter.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  pc.filter_flags |= kfFilterPvarReq;
	} else if (!memcmp(flagname_p2, "ests", 5)) {
	  if (!(pc.command_flags1 & kfCommand1Glm)) {
	    logerrprint("Error: --tests must be used with --glm.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if ((param_ct == 1) && (!strcmp(argv[arg_idx + 1], "all"))) {
	    pc.glm_info.flags |= kfGlmTestsAll;
	  } else {
	    reterr = parse_name_ranges(&(argv[arg_idx]), errstr_append, param_ct, 1, '-', &pc.glm_info.tests_range_list);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  logerrprint("Error: --tests is not implemented yet.\n");
	  reterr = kPglRetNotYetSupported;
	  goto main_ret_1;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;

      case 'u':
	if (!memcmp(flagname_p2, "pdate-sex", 10)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  reterr = alloc_fname(argv[arg_idx + 1], flagname_p, 0, &pc.update_sex_fname);
	  if (reterr) {
	    goto main_ret_1;
	  }
	  if (param_ct == 2) {
	    const char* cur_modif = argv[arg_idx + 2];
	    if (scan_posint_defcap(cur_modif, &pc.update_sex_colm2)) {
	      sprintf(g_logbuf, "Error: Invalid --update-sex column parameter '%s'. (This must be a positive integer.)\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'v':
	if (!memcmp(flagname_p2, "ar-min-qual", 12)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (scan_float(argv[arg_idx + 1], &pc.var_min_qual) || (pc.var_min_qual < 0.0)) {
	    sprintf(g_logbuf, "Error: Invalid --var-min-qual parameter '%s'.\n", argv[arg_idx + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  pc.var_min_qual *= 1 - kSmallEpsilon;
	  pc.filter_flags |= kfFilterPvarReq;
	} else if (!memcmp(flagname_p2, "ar-filter", 10)) {
	  if (param_ct) {
	    reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, 0x7fffffff, &var_filter_exceptions_flattened);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.misc_flags |= kfMiscExcludePvarFilterFail;
	  pc.filter_flags |= kfFilterPvarReq;
        } else if (!memcmp(flagname_p2, "cf", 3)) {
	  // permit accompanying .fam/.psam
	  // IIDs must match VCF sample line order
	  if ((load_params & (~kfLoadParamsPsam)) || xload) {
	    goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct == 2) {
	    const char* cur_modif = argv[arg_idx + 2];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen < 8) || memcmp(cur_modif, "dosage=", 7)) {
	      sprintf(g_logbuf, "Error: Invalid --vcf parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    reterr = cmdline_alloc_string(&(cur_modif[7]), argv[arg_idx], 4095, &vcf_dosage_import_field);
	    if (reterr) {
	      goto main_ret_1;
	    }
	    if (!is_alphanumeric(vcf_dosage_import_field)) {
	      logerrprint("Error: --vcf dosage= parameter is not alphanumeric.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (!strcmp(vcf_dosage_import_field, "GT")) {
	      logerrprint("Error: --vcf dosage= parameter cannot be 'GT'.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  const uint32_t slen = strlen(cur_modif);
	  if (slen > kPglFnamesize - 1) {
	    logerrprint("Error: --vcf filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pgenname, cur_modif, slen + 1);
	  xload = kfXloadVcf;
        } else if (!memcmp(flagname_p2, "cf-min-gp", 10)) {
	  logerrprint("Error: --vcf-min-gp is no longer supported.  Use --import-dosage-certainty\ninstead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if ((!memcmp(flagname_p2, "cf-min-gq", 10)) || (!memcmp(flagname_p2, "cf-min-dp", 10))) {
	  if (!(xload & kfXloadVcf)) {
	    // todo: support BCF too
	    LOGERRPRINTF("Error: --%s must be used with --vcf.\n", flagname_p);
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  uint32_t uii;
	  if (scan_uint_defcap(cur_modif, &uii)) {
	    sprintf(g_logbuf, "Error: Invalid --%s parameter '%s'.\n", flagname_p, cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (flagname_p2[7] == 'g') {
	    vcf_min_gq = uii;
	  } else {
	    vcf_min_dp = uii;
	  }
	} else if (!memcmp(flagname_p2, "cf-idspace-to", 14)) {
	  if (!(xload & (kfXloadVcf | kfXloadBcf))) {
	    logerrprint("Error: --vcf-idspace-to must be used with --vcf/--bcf.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (id_delim == ' ') {
	    logerrprint("Error: --vcf-idspace-to cannot be used when the --id-delim character is space.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  idspace_to = extract_char_param(argv[arg_idx + 1]);
	  if (!idspace_to) {
	    logerrprint("Error: --vcf-idspace-to parameter must be a single character.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (((unsigned char)idspace_to) <= ' ') {
	    logerrprint("Error: --vcf-idspace-to parameter must be a nonspace character.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	} else if (!memcmp(flagname_p2, "cf-half-call", 13)) {
	  if (!(xload & kfXloadVcf)) {
	    logerrprint("Error: --vcf-half-call must be used with --vcf.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* half_call_mode_str = argv[arg_idx + 1];
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
	    sprintf(g_logbuf, "Error: '%s' is not a valid mode for --vcf-half-call.\n", half_call_mode_str);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else if (!memcmp(flagname_p2, "cf-require-gt", 14)) {
	  if (!(xload & (kfXloadVcf | kfXloadBcf))) {
	    logerrprint("Error: --vcf-require-gt must be used with --vcf/--bcf.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  pc.misc_flags |= kfMiscVcfRequireGt;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "if", 3)) {
	  if (!(pc.command_flags1 & kfCommand1Glm)) {
	    logerrprint("Error: --vif must be used with --glm/--epistasis.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  if (!scanadv_double(cur_modif, &pc.vif_thresh)) {
	    sprintf(g_logbuf, "Error: Invalid --glm/--epistasis VIF threshold '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (pc.vif_thresh < 1.0) {
	    sprintf(g_logbuf, "Error: --glm/--epistasis VIF threshold '%s' too small (must be >= 1).\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else if (!memcmp(flagname_p2, "ariance-standardize", 20)) {
	  if (pc.pheno_transform_flags & kfPhenoTransformVstdCovar) {
	    logerrprint("Error: --variance-standardize cannot be used with --covar-variance-standardize.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (param_ct) {
	    reterr = alloc_and_flatten(&(argv[arg_idx + 1]), param_ct, 0x7fffffff, &pc.vstd_flattened);
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.pheno_transform_flags |= kfPhenoTransformVstdAll;
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "alidate", 8)) {
	  pc.command_flags1 |= kfCommand1Validate;
	  goto main_param_zero;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;
	
      case 'w':
	if (!memcmp(flagname_p2, "rite-snplist", 13)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    if (strcmp(cur_modif, "zs")) {
	      sprintf(g_logbuf, "Error: Invalid --write-snplist parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    pc.misc_flags |= kfMiscWriteSnplistZs;
	  }
	  pc.command_flags1 |= kfCommand1WriteSnplist;
	} else if (!memcmp(flagname_p2, "indow", 6)) {
	  if (!(pc.varid_snp || pc.varid_exclude_snp)) {
	    logerrprint("Error: --window must be used with --snp or --exclude-snp.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  char* cur_modif = argv[arg_idx + 1];
	  double dxx;
	  if (!scanadv_double(cur_modif, &dxx) || (dxx < 0)) {
	    sprintf(g_logbuf, "Error: Invalid --window parameter '%s'.\n", cur_modif);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  dxx *= 500 * (1 + kSmallEpsilon);
	  if (dxx > 2147483646) {
	    pc.window_bp = 0x7ffffffe;
	  } else {
	    pc.window_bp = (int32_t)dxx;
	  }
	  pc.filter_flags |= kfFilterNoSplitChr;
	  // no need to set kfFilterPvarReq due to --snp/--exclude-snp req.
	} else if (!memcmp(flagname_p2, "ithin", 6)) {
	  if (pc.misc_flags & kfMiscCatPhenoFamily) {
	    logerrprint("Error: --within cannot be used with --family.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
	    const char* cur_modif = argv[arg_idx + param_idx];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen == 7) && (!memcmp(cur_modif, "keep-", 5)) && match_upper_counted(&(cur_modif[5]), "NA", 2)) {
	      logerrprint("Error: --within's keep-NA modifier has been retired.  Rename that category in\nthe input file if you wish to keep it.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (param_idx == 1) {
	      reterr = alloc_fname(cur_modif, flagname_p, 0, &pc.within_fname);
	    } else {
	      if (is_reserved_pheno_name(cur_modif, cur_modif_slen)) {
		sprintf(g_logbuf, "Error: '%s' cannot be used as a categorical phenotype name.\n", cur_modif);
		goto main_ret_INVALID_CMDLINE_2A;
	      }
	      reterr = cmdline_alloc_string(cur_modif, argv[arg_idx], kMaxIdSlen, &pc.catpheno_name);
	    }
	    if (reterr) {
	      goto main_ret_1;
	    }
	  }
	  pc.filter_flags |= kfFilterPsamReq;
	} else if (!memcmp(flagname_p2, "rite-covar", 11)) {
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 0, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (param_ct) {
	    const char* cur_modif = argv[arg_idx + 1];
	    const uint32_t cur_modif_slen = strlen(cur_modif);
	    if ((cur_modif_slen >= 5) && (!memcmp(cur_modif, "cols=", 5))) {
	      reterr = parse_col_descriptor(&(cur_modif[5]), "maybesid\0sid\0maybeparents\0parents\0sex\0pheno1\0phenos\0", "write-covar", kfWriteCovarColMaybesid, kfWriteCovarColDefault, 0, &pc.write_covar_flags);
	      if (reterr) {
		goto main_ret_1;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --write-covar parameter '%s'.\n", cur_modif);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  } else {
	    pc.write_covar_flags |= kfWriteCovarColDefault;
	  }
	  pc.command_flags1 |= kfCommand1WriteCovar;
	} else if (!memcmp(flagname_p2, "arning-errcode", 15)) {
	  warning_errcode = 1;
	  goto main_param_zero;
	} else if (!memcmp(flagname_p2, "rite-cluster", 13)) {
	  logerrprint("Error: --write-cluster is retired.  Use e.g. --make-just-psam.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (!memcmp(flagname_p2, "ith-phenotype", 14)) {
	  logerrprint("Error: --with-phenotype is retired.  Use --write-covar cols=... instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else {
	  goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
	}
	break;

      case 'x':
	if (!memcmp(flagname_p2, "chr-model", 10)) {
	  if (!(pc.command_flags1 & (kfCommand1Glm | kfCommand1Score))) {
	    logerrprint("Error: --xchr-model must be used with --glm or --score.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (pc.glm_info.flags & (kfGlmGenotypic | kfGlmHethom | kfGlmDominant | kfGlmRecessive)) {
	    sprintf(g_logbuf, "Error: --xchr-model cannot be used with --glm %s.\n", (pc.glm_info.flags & kfGlmGenotypic)? "genotypic" : ((pc.glm_info.flags & kfGlmHethom)? "hethom" : ((pc.glm_info.flags & kfGlmDominant)? "dominant" : "recessive")));
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  if (enforce_param_ct_range(argv[arg_idx], param_ct, 1, 1)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  const char* cur_modif = argv[arg_idx + 1];
	  pc.xchr_model = ((uint32_t)extract_char_param(cur_modif)) - 48;
	  if (pc.xchr_model > 2) {
	    sprintf(g_logbuf, "Error: Invalid --xchr-model parameter '%s'.\n", cur_modif);
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
	  sprintf(g_logbuf, "Error: --%s doesn't accept parameters.\n", flagname_p);
	  goto main_ret_INVALID_CMDLINE_2A;
	}
      }
    } while ((++cur_flag_idx) < flag_ct);
    if (!outname_end) {
      outname_end = &(outname[6]);
    }
    
    if ((!pc.command_flags1) && (!(xload & (kfXloadVcf | kfXloadBcf | kfXloadOxBgen | kfXloadOxHaps | kfXloadOxSample | kfXloadPlink1Dosage | kfXloadGenDummy)))) {
      // add command_flags2 when needed
      goto main_ret_NULL_CALC;
    }
    if (!(load_params || xload)) {
      logerrprint("Error: No input dataset.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((xload & kfXloadOxGen) && (!(xload & kfXloadOxSample))) {
      // could permit .fam/.psam, but unless Oxford software supports that mode
      // it's pointless
      logerrprint("Error: --gen must be used with --sample or --data.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((xload & kfXloadOxSample) && (pc.misc_flags & kfMiscAffection01)) {
      // necessary for --data and --data --make-pgen to yield the same output
      logerrprint("Error: --data/--sample cannot be used with --1.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((pc.sample_sort_flags != kfSort0) && (!(pc.command_flags1 & (kfCommand1MakePlink2 | kfCommand1WriteCovar)))) {
      // todo: permit merge
      logerrprint("Error: --indiv-sort must be used with --make-{b}pgen/--make-bed/--write-covar\nor dataset merging.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((make_plink2_modifier & (kfMakePlink2MMask | kfMakePlink2TrimAlts | kfMakePgenEraseAlt2Plus | kfMakePgenErasePhase | kfMakePgenEraseDosage)) && (pc.command_flags1 & (~kfCommand1MakePlink2))) {
      logerrprint("Error: When the 'multiallelics=', 'trim-alts', and/or 'erase-...' modifier is\npresent, --make-bed/--make-{b}pgen cannot be combined with other commands.\n(Other filters are fine.)\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (aperm_present && (pc.command_flags1 & kfCommand1Glm) && (!(pc.glm_info.flags & kfGlmPerm))) {
      // If --aperm is present, at least one association analysis command which
      // supports adaptive permutation testing was also specified, but no
      // actual adaptive permutation test is happening, the user is likely to
      // be confused.  Produce a warning.  (Not an error since a sophisticated
      // user may want to use --script with different --aperm defaults.)
      logerrprint("Warning: --aperm only controls the settings for adaptive permutation tests; it\ndoes not cause such a test to be performed.  (Did you forget to add the 'perm'\nmodifier to an association analysis flag?)\n");
    }
    if ((pc.hard_call_thresh == 0xffffffffU) && (xload & (kfXloadVcf | kfXloadBcf | kfXloadOxGen | kfXloadOxBgen))) {
      if (pc.dosage_erase_thresh > (kDosageMid / 10)) {
	logerrprint("Error: --dosage-erase-threshold value cannot be larger than (default)\n--hard-call-threshold value.\n");
	goto main_ret_INVALID_CMDLINE_A;
      }
    } else {
      if (pc.dosage_erase_thresh > pc.hard_call_thresh) {
	logerrprint("Error: --dosage-erase-threshold value cannot be larger than\n--hard-call-threshold value.\n");
	goto main_ret_INVALID_CMDLINE_A;
      }
    }
    if ((oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond)) == (kfOxfordImportRefFirst | kfOxfordImportRefSecond)) {
      logerrprint("Error: --data/--{b}gen 'ref-first' and 'ref-second' modifiers cannot be used\ntogether.\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (!strcmp(g_missing_catname, g_output_missing_pheno)) {
      logerrprint("Error: --missing-catname and --output-missing-phenotype strings can't match.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((pc.misc_flags & kfMiscChrOverrideCmdline) && (!chr_info.chrset_source)) {
      logerrprint("Error: --chr-override requires an explicit chromosome set.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((xload & kfXloadPlink1Dosage) && (!(load_params & kfLoadParamsPsam))) {
      logerrprint("Error: --import-dosage requires a .fam file.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (!permit_multiple_inclusion_filters) {
      // Permit only one position- or ID-based variant inclusion filter, since
      // it's not immediately obvious whether the union or intersection should be
      // taken with multiple inclusion filters.
      // However, multiple exclusion filters are fine.  (Also,
      // --autosome{-par}/--chr is exempted since it's more obvious how they
      // interact with other filters.)
      const uint32_t inclusion_filter_extract = (pc.extract_fnames != nullptr);
      const uint32_t inclusion_filter_fromto_id = pc.varid_from || pc.varid_to;
      const uint32_t inclusion_filter_fromto_bp = (pc.from_bp != -1) || (pc.to_bp != -1);
      const uint32_t inclusion_filter_snpflag = (pc.varid_snp != nullptr);
      const uint32_t inclusion_filter_snpsflag = !!pc.snps_range_list.name_ct;
      if (inclusion_filter_extract + inclusion_filter_fromto_id + inclusion_filter_fromto_bp + inclusion_filter_snpflag + inclusion_filter_snpsflag > 1) {
	logerrprint("Error: Multiple variant inclusion filters specified (--extract, --from/--to,\n--from-bp/--to-bp, --snp, --snps).  Add --force-intersect if you really want\nthe intersection of these sets.  (If your variant IDs are unique, you can\nextract the union by e.g. running --write-snplist for each set, followed by\n--extract on all the .snplist files.)\n");
	goto main_ret_INVALID_CMDLINE_A;
      }
    }

    free_cond(subst_argv);
    free_cond(script_buf);
    free_cond(rerun_buf);
    free_cond(flag_buf);
    free_cond(flag_map);
    subst_argv = nullptr;
    script_buf = nullptr;
    rerun_buf = nullptr;
    flag_buf = nullptr;
    flag_map = nullptr;
    if (!rseeds) {
      uint32_t seed = (uint32_t)time(nullptr);
      sprintf(g_logbuf, "Random number seed: %u\n", seed);
      logstr(g_logbuf);
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
    
    uint64_t total_mb = detect_mb();
    if (!malloc_size_mb) {
      if (!total_mb) {
	malloc_size_mb = kBigstackDefaultMb;
      } else if (total_mb < (kBigstackMinMb * 2)) {
	malloc_size_mb = kBigstackMinMb;
      } else {
	malloc_size_mb = total_mb / 2;
      }
    }
    assert(malloc_size_mb >= (intptr_t)kBigstackMinMb);
#ifndef __LP64__
    if (malloc_size_mb > (intptr_t)kMalloc32bitMbMax) {
      malloc_size_mb = kMalloc32bitMbMax;
    }
#endif
    if (total_mb) {
      sprintf(g_logbuf, "%" PRIu64 " MB RAM detected; reserving %" PRIdPTR " MB for main workspace.\n", total_mb, malloc_size_mb);
    } else {
      sprintf(g_logbuf, "Failed to determine total system memory.  Attempting to reserve %" PRIuPTR " MB.\n", malloc_size_mb);
    }
    logprintb();
    uintptr_t malloc_mb_final;
    if (init_bigstack(malloc_size_mb, &malloc_mb_final, &bigstack_ua)) {
      goto main_ret_NOMEM;
    }
    g_input_missing_geno_ptr = &(g_one_char_strs[2 * ((unsigned char)input_missing_geno_char)]);
    g_output_missing_geno_ptr = &(g_one_char_strs[2 * ((unsigned char)output_missing_geno_char)]);
    if (((uintptr_t)malloc_size_mb) != malloc_mb_final) {
      if (memory_require) {
	goto main_ret_NOMEM;
      }
      LOGPRINTF("Allocated %" PRIuPTR " MB successfully, after larger attempt(s) failed.\n", malloc_mb_final);
    }

#ifndef _WIN32
    pthread_attr_init(&g_smallstack_thread_attr);
    pthread_attr_setstacksize(&g_smallstack_thread_attr, kDefaultThreadStack);
#endif
    // pigz_init(pc.max_thread_ct);

    print_end_time = 1;
    if (0) {
      // nonstandard cases (CNV, etc.) here
    } else {
      if (pc.filter_flags) {
	if (!pc.command_flags1) {
	  logerrprint("Error: Basic file conversions do not support regular filtering operations.\nRerun your command with --make-bed/--make-{b}pgen.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      }
      // print this here since some import functions are now multithreaded
      if (pc.max_thread_ct > 8) {
	LOGPRINTF("Using up to %u threads (change this with --threads).\n", pc.max_thread_ct);
      } else {
	// "1 compute thread" instead of "1 thread" since, when
	// max_thread_ct == 2, some code will use one I/O thread and one
	// compute thread.  Not worth the trouble of writing special-case code
	// to avoid that.  (also, with 2 cores, the I/O thread isn't
	// sufficiently busy to justify only 1 compute thread.)
	LOGPRINTF("Using %s%u compute thread%s.\n", (pc.max_thread_ct > 1)? "up to " : "", pc.max_thread_ct, (pc.max_thread_ct == 1)? "" : "s");
      }
      if (xload) {
	char* convname_end = outname_end;
	if (pc.command_flags1) {
	  if (pc.misc_flags & kfMiscKeepAutoconv) {
	    if (pc.misc_flags & kfMiscAffection01) {
	      logerrprint("Error: --1 cannot be used with --keep-autoconv.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if ((output_missing_geno_char != '.') && (output_missing_geno_char != input_missing_geno_char)) {
	      logerrprint("Error: --output-missing-genotype and --input-missing-genotype parameters cannot\nbe inconsistent when --keep-autoconv is specified.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    double dxx;
	    const char* num_end = scanadv_double(g_output_missing_pheno, &dxx);
	    if (num_end) {
	      if (dxx != ((double)pc.missing_pheno)) {
		logerrprint("Error: --output-missing-phenotype and --input-missing-phenotype parameters\ncannot be inconsistent when --keep-autoconv is specified.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else if (!is_nan_str(g_output_missing_pheno, strlen(g_output_missing_pheno))) {
	      logerrprint("Error: --output-missing-phenotype parameter must be numeric or 'NA' when\n--keep-autoconv is specified.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else {
	    convname_end = strcpya0(convname_end, "-temporary");
	  }
	} else {
	  pc.misc_flags |= kfMiscKeepAutoconv;
	}
	const uint32_t convname_slen = (uintptr_t)(convname_end - outname);
	const uint32_t psam_specified = (load_params & kfLoadParamsPsam);
	if (xload & kfXloadVcf) {
	  reterr = vcf_to_pgen(pgenname, psam_specified? psamname : nullptr, const_fid, vcf_dosage_import_field, pc.misc_flags, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, id_delim, idspace_to, vcf_min_gq, vcf_min_dp, vcf_half_call, pc.fam_cols, outname, convname_end, &chr_info);
	} else if (xload & kfXloadVcf) {
	  logerrprint("Error: --bcf is not implemented yet.\n");
	  reterr = kPglRetNotYetSupported;
	} else if (xload & kfXloadOxGen) {
	  reterr = ox_gen_to_pgen(pgenname, psamname, import_single_chr_str, ox_missing_code, pc.misc_flags, oxford_import_flags, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, outname, convname_end, &chr_info);
	} else if (xload & kfXloadOxBgen) {
	  reterr = ox_bgen_to_pgen(pgenname, psamname, const_fid, ox_missing_code, pc.misc_flags, oxford_import_flags, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, id_delim, idspace_to, pc.max_thread_ct, outname, convname_end, &chr_info);
	} else if (xload & kfXloadOxHaps) {
	  reterr = ox_hapslegend_to_pgen(pgenname, pvarname, psamname, import_single_chr_str, ox_missing_code, pc.misc_flags, oxford_import_flags, outname, convname_end, &chr_info);
	} else if (xload & kfXloadPlink1Dosage) {
	  reterr = plink1_dosage_to_pgen(pgenname, psamname, (xload & kfXloadMap)? pvarname : nullptr, import_single_chr_str, &plink1_dosage_info, pc.misc_flags, pc.fam_cols, pc.missing_pheno, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, pc.max_thread_ct, outname, convname_end, &chr_info);
	} else if (xload & kfXloadGenDummy) {
	  reterr = generate_dummy(&gendummy_info, pc.misc_flags, pc.hard_call_thresh, pc.dosage_erase_thresh, pc.max_thread_ct, outname, convname_end, &chr_info);
	}
	if (reterr || (!pc.command_flags1)) {
	  goto main_ret_1;
	}

	// todo: we have to skip this when merging is involved
	pc.hard_call_thresh = 0xffffffffU;
	
	strcpy(memcpya(pgenname, outname, convname_slen), ".pgen");
	strcpy(memcpya(pvarname, outname, convname_slen), ".pvar");
	if (!psam_specified) {
	  strcpy(memcpya(psamname, outname, convname_slen), ".psam");
	}
	if (!(pc.misc_flags & kfMiscKeepAutoconv)) {
	  if (push_llstr(pgenname, &file_delete_list) || push_llstr(pvarname, &file_delete_list)) {
	    goto main_ret_NOMEM;
	  }
	  if (!psam_specified) {
	    if (push_llstr(psamname, &file_delete_list)) {
	      goto main_ret_NOMEM;
	    }
	  }
	}
	*outname_end = '\0';
      }
      const uint32_t calc_all_req = (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Validate | kfCommand1WriteSnplist | kfCommand1WriteCovar))) || ((pc.command_flags1 & kfCommand1MakePlink2) && (make_plink2_modifier & (kfMakeBed | kfMakePgen)));
      if (calc_all_req || (pc.filter_flags & kfFilterAllReq)) {
	if ((!xload) && (load_params != kfLoadParamsPfileAll)) {
	  logerrprint("Error: A full fileset (.pgen/.bed + .pvar/.bim + .psam/.fam) is required for\nthis.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
      } else {
	// no genotype file required
	pgenname[0] = '\0';
	
	const uint32_t calc_pvar_req = (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1WriteCovar))) || ((pc.command_flags1 & kfCommand1MakePlink2) && (make_plink2_modifier & (kfMakeBed | kfMakeBim | kfMakePgen | kfMakePvar)));
	if (calc_pvar_req || (pc.filter_flags & kfFilterPvarReq)) {
	  if ((!xload) && (!(load_params & kfLoadParamsPvar))) {
	    logerrprint("Error: A .pvar/.bim file is required for this.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	} else {
	  pvarname[0] = '\0';
	}
	const uint32_t calc_psam_req = (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1WriteSnplist))) || ((pc.command_flags1 & kfCommand1MakePlink2) && (make_plink2_modifier & (kfMakeBed | kfMakeFam | kfMakePgen | kfMakePsam)));
	if (calc_psam_req || (pc.filter_flags & kfFilterPsamReq)) {
	  if ((!xload) && (!(load_params & kfLoadParamsPsam))) {
	    logerrprint("Error: A .psam/.fam file is required for this.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	} else {
	  psamname[0] = '\0';
	}
      }
      if (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Validate | kfCommand1WriteSnplist | kfCommand1WriteCovar))) {
	pc.filter_flags |= kfFilterNoSplitChr;
      }

      BLAS_SET_NUM_THREADS(1);
      reterr = plink2_core(var_filter_exceptions_flattened, require_pheno_flattened, require_covar_flattened, &pc, make_plink2_modifier, pgenname, psamname, pvarname, outname, outname_end, king_cutoff_fprefix, &chr_info);
    }    
  }
  while (0) {
  main_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  main_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  main_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  main_ret_INVALID_CMDLINE_UNRECOGNIZED:
    invalid_arg(argv[arg_idx]);
    logerrprintb();
    logerrprint(errstr_append);
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_INVALID_CMDLINE_INPUT_CONFLICT:
    LOGERRPRINTF("Error: --%s conflicts with another input flag.\n%s", flagname_p, errstr_append);
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_INVALID_CMDLINE_WWA:
    wordwrapb(0);
  main_ret_INVALID_CMDLINE_2A:
    logerrprintb();
  main_ret_INVALID_CMDLINE_A:
    logerrprint(errstr_append);
  main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_NULL_CALC:
    if (pc.filter_flags) {
      logerrprint("Warning: No output requested.  (Did you forget --make-bed/--make-{b}pgen?)\nExiting.\n");
    } else {
      logerrprint("Warning: No output requested.  Exiting.\n");
    }
  main_ret_NULL_CALC_0:
    fputs(g_cmdline_format_str, stdout);
    fputs(notestr_null_calc2, stdout);
    reterr = kPglRetSkipped;
    break;
  }
 main_ret_1:
  disp_exit_msg(reterr);
  while (0) {
  main_ret_NOMEM_NOLOG:
    print_ver();
  main_ret_NOMEM_NOLOG2:
    fputs(errstr_nomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
    }
    reterr = kPglRetNomem;
    break;
  main_ret_READ_FAIL_NOLOG:
    print_ver();
    fputs(errstr_read, stderr);
    reterr = kPglRetReadFail;
    break;
  }
 main_ret_NOLOG:
  fclose_cond(scriptfile);
  free_cond(vcf_dosage_import_field);
  free_cond(ox_missing_code);
  free_cond(import_single_chr_str);
  free_cond(const_fid);
  free_cond(require_covar_flattened);
  free_cond(require_pheno_flattened);
  free_cond(var_filter_exceptions_flattened);
  free_cond(rseeds);
  free_cond(subst_argv);
  free_cond(script_buf);
  free_cond(rerun_buf);
  free_cond(flag_buf);
  free_cond(flag_map);
  free_cond(king_cutoff_fprefix);
  free_cond(pc.covar_quantnorm_flattened);
  free_cond(pc.quantnorm_flattened);
  free_cond(pc.vstd_flattened);
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
  free_cond(pc.update_sex_fname);
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
  if (file_delete_list) {
    do {
      ll_str_t* llstr_ptr = file_delete_list->next;
      unlink(file_delete_list->ss);
      free(file_delete_list);
      file_delete_list = llstr_ptr;
    } while (file_delete_list);
  }
  cleanup_cmp_expr(&pc.remove_if_expr);
  cleanup_cmp_expr(&pc.keep_if_expr);
  cleanup_score(&pc.score_info);
  cleanup_glm(&pc.glm_info);
  cleanup_chr_info(&chr_info);
  cleanup_ld(&pc.ld_info);
  cleanup_range_list(&pc.covar_range_list);
  cleanup_range_list(&pc.pheno_range_list);
  cleanup_range_list(&pc.exclude_snps_range_list);
  cleanup_range_list(&pc.snps_range_list);
  if (warning_errcode && g_stderr_written_to && (!reterr)) {
    logerrprint("--warning-errcode: One or more warnings in this run; exiting with code 61.\n");
    reterr = kPglRetWarningErrcode;
  }
  if (cleanup_logfile(print_end_time) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  if (bigstack_ua) {
    free(bigstack_ua);
  }
  return (uint32_t)reterr;
}
