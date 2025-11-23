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


#include "include/SFMT.h"
#include "include/pgenlib_misc.h"
#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "include/plink2_bits.h"
#include "include/plink2_string.h"
#include "include/plink2_thread.h"
#include "include/plink2_zstfile.h"
#include "plink2_adjust.h"
#include "plink2_cmdline.h"
#include "plink2_common.h"
#include "plink2_compress_stream.h"
#include "plink2_data.h"
#include "plink2_export.h"
#include "plink2_family.h"
#include "plink2_fasta.h"
#include "plink2_filter.h"
#include "plink2_glm.h"
#include "plink2_glm_shared.h"
#include "plink2_help.h"
#include "plink2_import.h"
#include "plink2_import_legacy.h"
#include "plink2_ld.h"
#include "plink2_matrix.h"
#include "plink2_matrix_calc.h"
#include "plink2_merge.h"
#include "plink2_misc.h"
#include "plink2_psam.h"
#include "plink2_pvar.h"
#include "plink2_random.h"
#include "plink2_set.h"

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>  // time()
#include <unistd.h>  // unlink()

#if defined(__APPLE__) && defined(__x86_64__)
#  include <fenv.h>  // fesetenv()
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

static PREFER_CONSTEXPR char ver_str[] = "PLINK v2.0.0-a.7"
#ifdef NOLAPACK
  "NL"
#elif defined(LAPACK_ILP64)
  "LM"
#endif

#ifdef __LP64__
#  ifdef USE_CUDA
    " CUDA"
#  elif defined(USE_AVX2)
    " AVX2"
#  elif defined(__APPLE__) && !defined(__x86_64__)
    " M1"
#  elif defined(USE_SSE42)
    " SSE4.2"
#  else
    " 64-bit"
#  endif
#else
  " 32-bit"
#endif

#ifdef USE_MKL
  " Intel"
#elif defined(USE_AOCL)
  " AMD"
#endif
  " (23 Nov 2025)";
static PREFER_CONSTEXPR char ver_str2[] =
  // include leading space if day < 10, so character length stays the same
  ""

#ifdef NOLAPACK
#elif defined(LAPACK_ILP64)
#else
  "  "
#endif

#ifdef USE_MKL
#elif defined(USE_AOCL)
  "  "
#else
  "      "
#endif

#ifdef __LP64__
#  ifdef USE_CUDA
    "  "
#  elif defined(USE_AVX2)
    "  "
#  elif defined(__APPLE__) && !defined(__x86_64__)
    "    "
#  endif
#endif

  "       cog-genomics.org/plink/2.0/\n"
  "(C) 2005-2025 Shaun Purcell, Christopher Chang    GNU General Public License v3\n";
#ifdef HAS_CONSTEXPR
static_assert(CompileTimeSlen(ver_str) + CompileTimeSlen(ver_str2) == 160, "ver_str/ver_str2 must be updated");
#endif
static const char errstr_append[] = "For more info, try \"" PROG_NAME_STR " --help <flag name>\" or \"" PROG_NAME_STR " --help | more\".\n";

#ifndef NOLAPACK
static const char notestr_null_calc2[] = "Commands include --rm-dup list, --make-bpgen, --export, --freq, --geno-counts,\n--sample-counts, --missing, --hardy, --mendel, --het, --fst, --indep-pairwise,\n--r2-phased, --sample-diff, --make-king, --king-cutoff, --pmerge, --pgen-diff,\n--check-sex, --write-samples, --write-snplist, --make-grm-list, --pca, --glm,\n--adjust-file, --gwas-ssf, --pheno-svd, --clump, --score-list, --variant-score,\n--genotyping-rate, --pgen-info, --validate, and --zst-decompress.\n\n\"" PROG_NAME_STR " --help | more\" describes all functions.\n";
#else
// no --pca
static const char notestr_null_calc2[] = "Commands include --rm-dup list, --make-bpgen, --export, --freq, --geno-counts,\n--sample-counts, --missing, --hardy, --mendel, --het, --fst, --indep-pairwise,\n--r2-phased, --sample-diff, --make-king, --king-cutoff, --pmerge, --pgen-diff,\n--check-sex, --write-samples, --write-snplist, --make-grm-list, --glm,\n--adjust-file, --gwas-ssf, --clump, --score-list, --variant-score,\n--genotyping-rate, --pgen-info, --validate, and --zst-decompress.\n\n\"" PROG_NAME_STR " --help | more\" describes all functions.\n";
#endif

// multiallelics-already-joined + terminating null
CONSTI32(kMaxFlagBlen, 29);

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
  kfXloadGenDummy = (1 << 9),
  kfXloadPed = (1 << 10),
  kfXloadTped = (1 << 11),
  kfXloadEigGeno = (1 << 12),
  kfXloadEigInd = (1 << 13),
  kfXloadEigSnp = (1 << 14)
FLAGSET_DEF_END(Xload);


// assume for now that .pgen must always be accompanied by both .pvar and .psam
// currently does double-duty in tracking dependencies
FLAGSET64_DEF_START()
  kfFilter0,
  kfFilterAllReq = (1 << 0),
  kfFilterPvarReq = (1 << 1),
  kfFilterPsamReq = (1 << 2),
  kfFilterOpportunisticPgen = (1 << 3),
  kfFilterNonrefFlagsNeeded = (1 << 4),
  kfFilterNonrefFlagsNeededSet = kfFilterNonrefFlagsNeeded | kfFilterOpportunisticPgen,
  kfFilterNoSplitChr = (1 << 5),
  kfFilterExclFemales = (1 << 6),
  kfFilterExclMales = (1 << 7),
  kfFilterExclNosex = (1 << 8),
  kfFilterExclFounders = (1 << 9),
  kfFilterExclNonfounders = (1 << 10),
  kfFilterSnpsOnlyJustAcgt = (1 << 11),
  kfFilterExtractBed0 = (1 << 12),
  kfFilterExtractBed1 = (1 << 13),
  kfFilterExtractIntersectBed0 = (1 << 14),
  kfFilterExtractIntersectBed1 = (1 << 15),
  kfFilterExcludeBed0 = (1 << 16),
  kfFilterExcludeBed1 = (1 << 17),
  kfFilterSelectSidRepresentatives = (1 << 18),
  kfFilterMendel = (1 << 19)
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
  kfCommand1Ld = (1 << 18),
  kfCommand1PgenInfo = (1 << 19),
  kfCommand1RmDupList = (1 << 20),
  kfCommand1Sdiff = (1 << 21),
  kfCommand1SampleCounts = (1 << 22),
  kfCommand1Vscore = (1 << 23),
  kfCommand1Het = (1 << 24),
  kfCommand1Fst = (1 << 25),
  kfCommand1Pmerge = (1 << 26),
  kfCommand1PgenDiff = (1 << 27),
  kfCommand1Clump = (1 << 28),
  kfCommand1Vcor = (1 << 29),
  kfCommand1PhenoSvd = (1 << 30),
  kfCommand1CheckOrImputeSex = (1U << 31),
  kfCommand1MendelReport = (1LLU << 32)
FLAGSET64_DEF_END(Command1Flags);

void PgenInfoPrint(const char* pgenname, const PgenFileInfo* pgfip, PgenExtensionLl* header_exts, PgenHeaderCtrl header_ctrl, uint32_t max_allele_ct) {
  logprintfww("--pgen-info on %s:\n", pgenname);
  logprintf("  Variants: %u\n", pgfip->raw_variant_ct);
  logprintf("  Samples: %u\n", pgfip->raw_sample_ct);
  const uint32_t nonref_flags_status = header_ctrl >> 6;
  if (!nonref_flags_status) {
    logputs("  REF allele known/provisional status not stored in .pgen\n");
  } else if (nonref_flags_status == 1) {
    logputs("  REF alleles are all known\n");
  } else if (nonref_flags_status == 2) {
    logputs("  REF alleles are all provisional\n");
  } else {
    // could report exact counts of each
    logputs("  REF alleles are a mix of known and provisional\n");
  }
  if (max_allele_ct >= UINT32_MAXM1) {
    if (max_allele_ct == UINT32_MAX) {
      logputs("  Maximum allele count for a single variant: >2, not explicitly stored\n");
    } else {
      logputs("  Maximum allele count for a single variant: not explicitly stored\n");
    }
  } else {
    logprintf("  Maximum allele count for a single variant: %u\n", max_allele_ct);
  }
  if (pgfip->gflags & kfPgenGlobalHardcallPhasePresent) {
    logputs("  Explicitly phased hardcalls present\n");
  } else {
    logputs("  No hardcalls are explicitly phased\n");
  }
  if (pgfip->gflags & kfPgenGlobalDosagePresent) {
    if (pgfip->gflags & kfPgenGlobalDosagePhasePresent) {
      logputs("  Explicitly phased dosages present\n");
    } else {
      logputs("  Dosage present, none explicitly phased\n");
    }
  } else {
    logputs("  No dosages present\n");
  }
  assert(header_exts && (header_exts->type_idx == 1));
  if (header_exts->size != ~0LLU) {
    char* write_iter = strcpya_k(g_logbuf, "  Writer: ");
    write_iter = memcpya(write_iter, header_exts->contents, header_exts->size);
    strcpy_k(write_iter, "\n");
    WordWrapB(0);
    logputsb();
  }
}

PglErr PgenInfoStandalone(const char* pgenname, const char* pginame) {
  PgenFileInfo pgfi;
  PglErr reterr = kPglRetSuccess;
  PreinitPgfi(&pgfi);
  PgenExtensionLl ext_slot; // shouldn't have shorter lifetime than pgfi
  ext_slot.contents = nullptr;
  {
    PgenHeaderCtrl header_ctrl;
    uintptr_t cur_alloc_cacheline_ct;
    reterr = PgfiInitPhase1(pgenname, pginame, UINT32_MAX, UINT32_MAX, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, g_logbuf);
    if (unlikely(reterr)) {
      if ((reterr == kPglRetSampleMajorBed) || (reterr == kPglRetImproperFunctionCall)) {
        logerrputs("Warning: Skipping --pgen-info since a .bed file was provided.\n");
        reterr = kPglRetSuccess;
      } else {
        WordWrapB(0);
        logerrputsb();
      }
      goto PgenInfoStandalone_ret_1;
    }
    const uint32_t raw_variant_ct = pgfi.raw_variant_ct;
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    unsigned char* pgfi_alloc;
    if (unlikely(bigstack_alloc_uc(cur_alloc_cacheline_ct * kCacheline, &pgfi_alloc) ||
                 bigstack_alloc_w(raw_variant_ct + 1, &pgfi.allele_idx_offsets) ||
                 bigstack_alloc_w(raw_variant_ctl, &pgfi.nonref_flags))) {
      goto PgenInfoStandalone_ret_NOMEM;
    }
    ext_slot.next = nullptr;
    ext_slot.size = ~0LLU;
    ext_slot.type_idx = 1;
    PgenExtensionLl* header_exts = &ext_slot;
    uintptr_t pgr_alloc_cacheline_ct = 0;
    uint32_t max_vrec_width;
    reterr = PgfiInitPhase2Ex(header_ctrl, 0, 0, 1, 0, raw_variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, header_exts, nullptr, &pgr_alloc_cacheline_ct, g_logbuf);
    if (unlikely(reterr)) {
      WordWrapB(0);
      logerrputsb();
      goto PgenInfoStandalone_ret_1;
    }
    if (header_exts->size != ~0LLU) {
#ifndef __LP64__
      if (unlikely(header_exts->size >= (1LLU << 31))) {
        goto PgenInfoStandalone_ret_NOMEM;
      }
#endif
      if (unlikely(pgl_malloc(header_exts->size, &(header_exts->contents)))) {
        goto PgenInfoStandalone_ret_NOMEM;
      }
      reterr = PgfiInitLoadExts(header_ctrl, &pgfi, header_exts, nullptr, g_logbuf);
      if (unlikely(reterr)) {
        WordWrapB(0);
        logerrputsb();
        goto PgenInfoStandalone_ret_1;
      }
    }
    uint32_t max_allele_ct = 2;
    if (pgfi.gflags & kfPgenGlobalMultiallelicHardcallFound) {
      max_allele_ct = UINT32_MAX;
    } else if (pgfi.gflags & kfPgenGlobalDosagePresent) {
      max_allele_ct = UINT32_MAXM1;
    }
    PgenInfoPrint(pgenname, &pgfi, header_exts, header_ctrl, max_allele_ct);
  }
  while (0) {
  PgenInfoStandalone_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 PgenInfoStandalone_ret_1:
  free_cond(ext_slot.contents);
  CleanupPgfi2(pgenname, &pgfi, &reterr);
  // no BigstackReset() needed?
  return reterr;
}

typedef struct Plink2CmdlineStruct {
  NONCOPYABLE(Plink2CmdlineStruct);
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

  // .log should show where --chr, --exclude-if-info, etc. are applied
  LoadFilterLogFlags load_filter_log_flags;

  Command1Flags command_flags1;
  PvarPsamFlags pvar_psam_flags;
  SortMode sample_sort_mode;
  SortMode sort_vars_mode;
  GrmFlags grm_flags;
  double grm_sparse_cutoff;
  PcaFlags pca_flags;
  WriteCovarFlags write_covar_flags;
  PhenoTransformFlags pheno_transform_flags;
  FaFlags fa_flags;
  RangeList snps_range_list;
  RangeList exclude_snps_range_list;
  RangeList pheno_range_list;
  RangeList covar_range_list;
  RangeList vscore_col_idx_range_list;
  FamCol fam_cols;
  UpdateAllelesInfo update_alleles_info;
  UpdateSexInfo update_sex_info;
  LdInfo ld_info;
  SdiffInfo sdiff_info;
  KingFlags king_flags;
  double king_cutoff;
  double king_table_filter;
  double king_table_subset_thresh;
  FreqRptFlags freq_rpt_flags;
  MissingRptFlags missing_rpt_flags;
  GenoCountsFlags geno_counts_flags;
  HardyFlags hardy_flags;
  HetFlags het_flags;
  SampleCountsFlags sample_counts_flags;
  RecoverVarIdsFlags recover_var_ids_flags;
  VscoreFlags vscore_flags;
  AlleleAlphanumFlags allele_alphanum_flags;
  RmDupMode rmdup_mode;
  STD_ARRAY_DECL(FreqFilterMode, 4, filter_modes);
  SelectSidMissingnessMode select_sid_missingness_mode;
  SelectSidTiebreakMode select_sid_tiebreak_mode;
  GlmInfo glm_info;
  AdjustInfo adjust_info;
  ScoreInfo score_info;
  FstInfo fst_info;
  PgenDiffInfo pgen_diff_info;
  APerm aperm;
  CmpExpr keep_if_expr;
  CmpExpr remove_if_expr;
  CmpExpr extract_if_info_expr;
  CmpExpr exclude_if_info_expr;
  ExtractColCondInfo extract_col_cond_info;
  ExportfInfo exportf_info;
  GwasSsfInfo gwas_ssf_info;
  ClumpInfo clump_info;
  VcorInfo vcor_info;
  PhenoSvdInfo pheno_svd_info;
  CheckSexInfo check_sex_info;
  MendelInfo mendel_info;
  FlipInfo flip_info;
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

  double ln_pfilter;
  double output_min_ln;
  double vif_thresh;
  double mind_thresh;
  double geno_thresh;
  double hwe_ln_thresh;
  double hwe_sample_size_term; // -1 if unspecified
  double mach_r2_min;
  double mach_r2_max;
  double minimac3_r2_min;
  double minimac3_r2_max;
  double af_pseudocount;
  double min_maf;
  double max_maf;
  double thin_keep_prob;
  double thin_keep_sample_prob;
  uint64_t min_allele_ddosage;
  uint64_t max_allele_ddosage;
  int32_t missing_pheno;
  int32_t from_bp;
  int32_t to_bp;
  int32_t window_bp;
  uint32_t pca_ct;
  uint32_t xchr_model;
  uint32_t max_thread_ct;
  uint32_t parallel_idx;
  uint32_t parallel_tot;
  uint32_t mwithin_val;
  uint32_t min_bp_space;
  uint32_t thin_keep_ct;
  uint32_t thin_keep_sample_ct;
  uint32_t keep_col_match_num;
  uint32_t filter_min_allele_ct;
  uint32_t filter_max_allele_ct;
  uint32_t bed_border_bp;
  char input_missing_geno_char;
  char output_missing_geno_char;
  char legacy_output_missing_geno_char;
  // might want input_missing_catname and/or output_missing_catname later, but
  // let's start with the simplest implementation
  char missing_catname[kMaxMissingPhenostrBlen];
  char output_missing_pheno[kMaxMissingPhenostrBlen];
  char legacy_output_missing_pheno[kMaxMissingPhenostrBlen];

  char* pginame;
  char* var_filter_exceptions_flattened;
  char* varid_template_str;
  char* varid_multi_template_str;
  char* varid_multi_nonsnp_template_str;
  char* missing_varid_match;
  char* varid_from;
  char* varid_to;
  char* varid_snp;
  char* varid_exclude_snp;
  char* pheno_fname;
  char* covar_fname;
  char* extract_fnames;
  char* extract_intersect_fnames;
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
  char* fa_fname;
  char* king_table_subset_fname;
  char* king_table_require_fnames;
  char* require_info_flattened;
  char* require_no_info_flattened;
  char* keep_col_match_fname;
  char* keep_col_match_flattened;
  char* keep_col_match_name;
  char* update_sample_ids_fname;
  char* update_parental_ids_fname;
  char* recover_var_ids_fname;
  char* vscore_fname;
  char* indep_preferred_fname;
  char* not_pheno_flattened;
  char* not_covar_flattened;
  char* indv_str;
  char* rename_chrs_fname;
  char* zero_cluster_fname;
  char* zero_cluster_phenoname;
  TwoColParams* ref_allele_flag;
  TwoColParams* alt_allele_flag;
  TwoColParams* update_chr_flag;
  TwoColParams* update_map_flag;
  TwoColParams* update_name_flag;
} Plink2Cmdline;

// er, probably time to just always initialize this...
uint32_t SingleVariantLoaderIsNeeded(const char* king_cutoff_fprefix, Command1Flags command_flags1, MakePlink2Flags make_plink2_flags, RmDupMode rmdup_mode, double hwe_ln_thresh) {
  return (command_flags1 & (kfCommand1Exportf | kfCommand1MakeKing | kfCommand1GenoCounts | kfCommand1LdPrune | kfCommand1Validate | kfCommand1Pca | kfCommand1MakeRel | kfCommand1Glm | kfCommand1Score | kfCommand1Ld | kfCommand1Hardy | kfCommand1Sdiff | kfCommand1PgenDiff | kfCommand1Clump | kfCommand1Vcor)) ||
    ((command_flags1 & kfCommand1MakePlink2) && (make_plink2_flags & kfMakePgen)) ||
    ((command_flags1 & kfCommand1KingCutoff) && (!king_cutoff_fprefix)) ||
    (rmdup_mode != kRmDup0) ||
    (hwe_ln_thresh != -DBL_MAX);
}


uint32_t DecentAlleleFreqsAreNeeded(Command1Flags command_flags1, CheckSexFlags check_sex_flags, HetFlags het_flags, ScoreFlags score_flags) {
  // Keep this in sync with --error-on-freq-calc.
  return (command_flags1 & (kfCommand1Pca | kfCommand1MakeRel)) ||
    (check_sex_flags & kfCheckSexUseX) ||
    ((command_flags1 & kfCommand1Score) && ((!(score_flags & kfScoreNoMeanimpute)) || (score_flags & (kfScoreCenter | kfScoreVarianceStandardize)))) ||
    ((command_flags1 & kfCommand1Het) && (!(het_flags & kfHetSmallSample)));
}

// not actually needed for e.g. --hardy, --hwe, etc. if no multiallelic
// variants are retained, but let's keep this simpler for now
uint32_t MajAllelesAreNeeded(Command1Flags command_flags1, PcaFlags pca_flags, GlmFlags glm_flags, VcorFlags vcor_flags) {
  // Keep this in sync with --error-on-freq-calc.
  return (command_flags1 & (kfCommand1LdPrune | kfCommand1Ld)) ||
    ((command_flags1 & kfCommand1Pca) && (pca_flags & kfPcaBiallelicVarWts)) ||
    ((command_flags1 & kfCommand1Glm) && (!(glm_flags & kfGlmOmitRef))) ||
    ((command_flags1 & kfCommand1Vcor) && ((!(vcor_flags & kfVcorRefBased)) || (vcor_flags & (kfVcorColMaj | kfVcorColNonmaj))));
}

// only needs to cover cases not captured by DecentAlleleFreqsAreNeeded() or
// MajAllelesAreNeeded()
uint32_t IndecentAlleleFreqsAreNeeded(Command1Flags command_flags1, VcorFlags vcor_flags, double min_maf, double max_maf) {
  // Keep this in sync with --error-on-freq-calc.
  // Vscore could go either here or in the decent bucket
  return (command_flags1 & kfCommand1Vscore) ||
    ((command_flags1 & kfCommand1Vcor) && (vcor_flags & kfVcorColFreq)) ||
    (min_maf != 0.0) ||
    (max_maf != 1.0);
}


uint32_t GetFirstHaploidUidx(const ChrInfo* cip, UnsortedVar vpos_sortstatus) {
  // returns 0x7fffffff if no X/haploid chromosomes present
  if (!(vpos_sortstatus & kfUnsortedVarSplitChr)) {
    const uint32_t chr_ct = cip->chr_ct;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if (IsSet(cip->haploid_mask, chr_idx)) {
        return cip->chr_fo_vidx_start[chr_fo_idx];
      }
    }
  }
  return 0x7fffffff;
}

uint32_t AlleleDosagesAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, uint32_t afreq_needed, uint64_t min_allele_ddosage, uint64_t max_allele_ddosage, uint32_t* regular_freqcounts_neededp) {
  if (!(misc_flags & kfMiscNonfounders)) {
    return 0;
  }
  if ((command_flags1 & kfCommand1AlleleFreq) ||
      (misc_flags & kfMiscMajRef) ||
      min_allele_ddosage ||
      (max_allele_ddosage != (~0LLU))) {
    // Keep this in sync with --error-on-freq-calc.
    *regular_freqcounts_neededp = 1;
    return 1;
  }
  return afreq_needed;
}

uint32_t FounderAlleleDosagesAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, uint32_t afreq_needed, uint64_t min_allele_ddosage, uint64_t max_allele_ddosage, uint32_t* regular_freqcounts_neededp) {
  if (misc_flags & kfMiscNonfounders) {
    return 0;
  }
  if ((command_flags1 & kfCommand1AlleleFreq) ||
      (misc_flags & kfMiscMajRef) ||
      min_allele_ddosage ||
      (max_allele_ddosage != (~0LLU))) {
    // Keep this in sync with --error-on-freq-calc.
    *regular_freqcounts_neededp = 1;
    return 1;
  }
  return afreq_needed;
}

uint32_t SampleMissingDosageCtsAreNeeded(MiscFlags misc_flags, uint32_t smaj_missing_geno_report_requested, double mind_thresh, SelectSidMissingnessMode select_sid_missingness_mode, MissingRptFlags missing_rpt_flags) {
  return ((mind_thresh != 1.0) && (misc_flags & kfMiscMindDosage)) ||
    (smaj_missing_geno_report_requested && (missing_rpt_flags & (kfMissingRptScolNmissDosage | kfMissingRptScolFmissDosage))) ||
    (select_sid_missingness_mode == kSelectSidMissingnessDosage);
}

uint32_t VariantMissingHcCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double geno_thresh, MissingRptFlags missing_rpt_flags) {
  // Keep this in sync with --error-on-freq-calc.
  return ((command_flags1 & kfCommand1GenotypingRate) && (!(misc_flags & kfMiscGenotypingRateDosage))) ||
    ((command_flags1 & kfCommand1MissingReport) && (missing_rpt_flags & (kfMissingRptVcolNmiss | kfMissingRptVcolNmissHh | kfMissingRptVcolHethap | kfMissingRptVcolFmiss | kfMissingRptVcolFmissHh | kfMissingRptVcolFhethap))) ||
    ((geno_thresh != 1.0) && (!(misc_flags & kfMiscGenoDosage)));
}

uint32_t VariantHethapCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double geno_thresh, MissingRptFlags missing_rpt_flags, uint32_t first_hap_uidx) {
  // bugfix (20 Oct 2024): --geno hh-missing condition was backwards
  return (first_hap_uidx != 0x7fffffff) &&
    (((command_flags1 & kfCommand1MissingReport) && (missing_rpt_flags & (kfMissingRptVcolNmissHh | kfMissingRptVcolHethap | kfMissingRptVcolFmissHh | kfMissingRptVcolFhethap))) ||
     ((geno_thresh != 1.0) && (misc_flags & kfMiscGenoHhMissing)));
}

uint32_t VariantMissingDosageCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double geno_thresh, MissingRptFlags missing_rpt_flags) {
  // Keep this in sync with --error-on-freq-calc.
  return ((command_flags1 & kfCommand1GenotypingRate) && (misc_flags & kfMiscGenotypingRateDosage)) ||
    ((command_flags1 & kfCommand1MissingReport) && (!(missing_rpt_flags & kfMissingRptSampleOnly)) && (missing_rpt_flags & (kfMissingRptVcolNmissDosage | kfMissingRptVcolFmissDosage))) ||
    ((geno_thresh != 1.0) && (misc_flags & kfMiscGenoDosage));
}

// can simplify --geno-counts all-biallelic case, but let's first make sure the
// general case works for multiallelic variants
uint32_t RawGenoCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double hwe_ln_thresh) {
  // Keep this in sync with --error-on-freq-calc.
  return (command_flags1 & kfCommand1GenoCounts) ||
    ((misc_flags & kfMiscNonfounders) && ((command_flags1 & kfCommand1Hardy) || (hwe_ln_thresh != -DBL_MAX)));
}

uint32_t FounderRawGenoCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double hwe_ln_thresh) {
  // Keep this in sync with --error-on-freq-calc.
  return (!(misc_flags & kfMiscNonfounders)) && ((command_flags1 & kfCommand1Hardy) || (hwe_ln_thresh != -DBL_MAX));
}

uint32_t InfoReloadIsNeeded(Command1Flags command_flags1, PvarPsamFlags pvar_psam_flags, ExportfFlags exportf_flags, RmDupMode rmdup_mode) {
  return ((command_flags1 & kfCommand1MakePlink2) && (pvar_psam_flags & kfPvarColXinfo)) ||
    ((command_flags1 & kfCommand1Exportf) && (exportf_flags & (kfExportfVcf | kfExportfBcf))) ||
    (rmdup_mode != kRmDup0);
}

uint32_t GrmKeepIsNeeded(Command1Flags command_flags1, PcaFlags pca_flags) {
  return ((command_flags1 & kfCommand1Pca) && (!(pca_flags & kfPcaApprox)));
}

uint32_t AlleleUniquenessCheckIsNeeded(const char* update_alleles_fname, const char* read_freq_fname, const TwoColParams* ref_allele_flag, const TwoColParams* alt_allele_flag, Command1Flags command_flags1, RmDupMode rmdup_mode, FaFlags fa_flags) {
  return (command_flags1 & (kfCommand1PgenDiff | kfCommand1Score | kfCommand1Clump)) || update_alleles_fname || read_freq_fname || (rmdup_mode != kRmDup0) || ref_allele_flag || alt_allele_flag || (fa_flags & kfFaRefFrom);
}

void ReportGenotypingRate(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_missing_cts, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t y_sample_ct, uint32_t variant_ct, uint32_t is_dosage) {
  // Each variant has equal weight.
  // By default, only males are considered on chrY; with
  // --y-nosex-missing-stats, missing-sex samples are also counted.  If
  // y_sample_ct == 0, chrY is not considered at all.

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
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  uint32_t is_y = 0;
  for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
    const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
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
  if (y_sample_ct && variant_ct_y) {
    const uint64_t nony_possible_obs = (variant_ct - variant_ct_y) * S_CAST(uint64_t, sample_ct);
    const uint64_t y_possible_obs = variant_ct_y * S_CAST(uint64_t, y_sample_ct);
    genotyping_rate = u63tod(nony_possible_obs - tot_nony_missing) / u31tod(sample_ct) + u63tod(y_possible_obs - tot_y_missing) / u31tod(y_sample_ct);
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

PglErr ApplyVariantBpFilters(const char* extract_fnames, const char* extract_intersect_fnames, const char* exclude_fnames, const ChrInfo* cip, const uint32_t* variant_bps, int32_t from_bp, int32_t to_bp, uint32_t bed_border_bp, uint32_t raw_variant_ct, FilterFlags filter_flags, UnsortedVar vpos_sortstatus, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  if (!(*variant_ct_ptr)) {
    return kPglRetSuccess;
  }
  if ((from_bp != -1) || (to_bp != -1)) {
    if (unlikely(vpos_sortstatus & kfUnsortedVarBp)) {
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
      variant_uidx_start = LowerBoundConstrainedNonemptyU32(variant_bps, variant_uidx_start, variant_uidx_end, from_bp);
    }
    if ((to_bp != -1) && (variant_uidx_start < variant_uidx_end)) {
      variant_uidx_end = LowerBoundConstrainedNonemptyU32(variant_bps, variant_uidx_start, variant_uidx_end, 1 + to_bp);
    }
    if (variant_uidx_start) {
      ClearBitsNz(0, variant_uidx_start, variant_include);
    }
    if (variant_uidx_end < raw_variant_ct) {
      ClearBitsNz(variant_uidx_end, raw_variant_ct, variant_include);
    }
    *variant_ct_ptr = PopcountBitRange(variant_include, variant_uidx_start, variant_uidx_end);
  }
  if (extract_fnames && (filter_flags & (kfFilterExtractBed0 | kfFilterExtractBed1))) {
    if (unlikely(vpos_sortstatus & kfUnsortedVarBp)) {
      logerrputs("Error: '--extract bed0'/'--extract bed1' requires a sorted .pvar/.bim.  Retry\nthis command after using --make-pgen/--make-bed + --sort-vars to sort your\ndata.\n");
      return kPglRetInconsistentInput;
    }
    PglErr reterr = ExtractExcludeRange(extract_fnames, cip, variant_bps, raw_variant_ct, kVfilterExtract, (filter_flags / kfFilterExtractBed0) & 1, bed_border_bp, max_thread_ct, variant_include, variant_ct_ptr);
    if (unlikely(reterr)) {
      return reterr;
    }
  }
  if (extract_intersect_fnames && (filter_flags & (kfFilterExtractIntersectBed0 | kfFilterExtractIntersectBed1))) {
    if (unlikely(vpos_sortstatus & kfUnsortedVarBp)) {
      logerrputs("Error: '--extract-intersect bed0'/'--extract-intersect bed1' requires a sorted\n.pvar/.bim.  Retry this command after using --make-pgen/--make-bed +\n--sort-vars to sort your data.\n");
      return kPglRetInconsistentInput;
    }
    PglErr reterr = ExtractExcludeRange(extract_intersect_fnames, cip, variant_bps, raw_variant_ct, kVfilterExtractIntersect, (filter_flags / kfFilterExtractIntersectBed0) & 1, bed_border_bp, max_thread_ct, variant_include, variant_ct_ptr);
    if (unlikely(reterr)) {
      return reterr;
    }
  }
  if (exclude_fnames && (filter_flags & (kfFilterExcludeBed0 | kfFilterExcludeBed1))) {
    if (unlikely(vpos_sortstatus & kfUnsortedVarBp)) {
      logerrputs("Error: '--exclude bed0'/'--exclude bed1' requires a sorted .pvar/.bim.  Retry\nthis command after using --make-pgen/--make-bed + --sort-vars to sort your\ndata.\n");
      return kPglRetInconsistentInput;
    }
    PglErr reterr = ExtractExcludeRange(exclude_fnames, cip, variant_bps, raw_variant_ct, kVfilterExclude, (filter_flags / kfFilterExcludeBed0) & 1, bed_border_bp, max_thread_ct, variant_include, variant_ct_ptr);
    if (unlikely(reterr)) {
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
static_assert(kPglMaxAlleleCt == 255, "Plink2Core() --maj-ref needs to be updated.");
PglErr Plink2Core(const Plink2Cmdline* pcp, MakePlink2Flags make_plink2_flags, char* pgenname, char* psamname, char* pvarname, char* outname, char* outname_end, char* king_cutoff_fprefix, ChrInfo* cip, sfmt_t* sfmtp) {
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
  PgenExtensionLl ext_slot; // shouldn't have shorter lifetime than pgfi
  ext_slot.contents = nullptr;
  {
    uint32_t pvar_renamed = 0;
    if ((make_plink2_flags & (kfMakeBed | kfMakePgen)) ||
        (pcp->exportf_info.flags & kfExportfIndMajorBed)) {
      uint32_t fname_slen;
#ifdef _WIN32
      fname_slen = GetFullPathName(pgenname, kPglFnamesize, g_textbuf, nullptr);
      if (unlikely((!fname_slen) || (fname_slen > kPglFnamesize)))
#else
      if (unlikely(!realpath(pgenname, g_textbuf)))
#endif
      {
        logerrprintfww(kErrprintfFopen, pgenname, strerror(errno));
        goto Plink2Core_ret_OPEN_FAIL;
      }
      uint32_t pgen_rename = 0;
      if (make_plink2_flags & kfMakePgen) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
        pgen_rename = RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
      }
      if ((!pgen_rename) && ((make_plink2_flags & kfMakeBed) || (pcp->exportf_info.flags & kfExportfIndMajorBed))) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".bed");
        pgen_rename = RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
      }
      if (pgen_rename) {
        logerrprintf("Warning: --make-%s input and output filenames match.  Appending '~' to input\nfilenames.\n", (make_plink2_flags & kfMakeBed)? "bed" : ((make_plink2_flags & kfMakePvar)? "pgen" : "bpgen"));
        fname_slen = strlen(pgenname);
        memcpy(g_textbuf, pgenname, fname_slen + 1);
        snprintf(&(pgenname[fname_slen]), 2, "~");
        if (unlikely(rename(g_textbuf, pgenname))) {
          logerrprintf("Error: Failed to append '~' to input .bed/.pgen filename: %s.\n", strerror(errno));
          goto Plink2Core_ret_OPEN_FAIL;
        }
        // if/when --make-pgen can generate .pgen.pgi files, change pginame to
        // a Plink2Core() parameter, and conditionally rename the file here
        fname_slen = strlen(pvarname);
        memcpy(g_textbuf, pvarname, fname_slen + 1);
        snprintf(&(pvarname[fname_slen]), 2, "~");
        if (unlikely(rename(g_textbuf, pvarname))) {
          logerrprintf("Error: Failed to append '~' to input .bim/.pvar filename: %s.\n", strerror(errno));
          goto Plink2Core_ret_OPEN_FAIL;
        }
        pvar_renamed = 1;
        fname_slen = strlen(psamname);
        memcpy(g_textbuf, psamname, fname_slen + 1);
        snprintf(&(psamname[fname_slen]), 2, "~");
        if (unlikely(rename(g_textbuf, psamname))) {
          logerrprintf("Error: Failed to append '~' to input .fam/.psam filename: %s.\n", strerror(errno));
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

    uint32_t update_sample_ids_empty = 0;
    uint32_t update_parental_ids_empty = 0;

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
      if (pcp->update_sample_ids_fname) {
        reterr = PrescanSampleIds(pcp->update_sample_ids_fname, &pii.sii);
        if (reterr) {
          if (unlikely(reterr != kPglRetEof)) {
            goto Plink2Core_ret_1;
          }
          reterr = kPglRetSuccess;
          update_sample_ids_empty = 1;
        }
      }
      if (pcp->update_parental_ids_fname) {
        reterr = PrescanParentalIds(pcp->update_parental_ids_fname, pcp->max_thread_ct, &pii.parental_id_info);
        if (reterr) {
          if (unlikely(reterr != kPglRetEof)) {
            goto Plink2Core_ret_1;
          }
          reterr = kPglRetSuccess;
          update_parental_ids_empty = 1;
        }
      }
      // update (26 Nov 2017): change --no-pheno to also apply to .psam file
      const uint32_t ignore_psam_phenos = (!(pcp->fam_cols & kfFamCol6)) || (pcp->pheno_fname && pcp->pheno_range_list.name_ct);
      reterr = LoadPsam(psamname, pcp->pheno_fname? nullptr : &(pcp->pheno_range_list), pcp->missing_catname, pcp->fam_cols, ignore_psam_phenos? 0 : 0x7fffffff, pcp->missing_pheno, (pcp->misc_flags / kfMiscAffection01) & 1, (pcp->misc_flags / kfMiscNoCategorical) & 1, (pcp->misc_flags / kfMiscNeg9PhenoReallyMissing) & 1, pcp->max_thread_ct, &pii, &sample_include, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
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
    uint32_t info_reload_slen = InfoReloadIsNeeded(pcp->command_flags1, pcp->pvar_psam_flags, pcp->exportf_info.flags, pcp->rmdup_mode);
    uintptr_t* allele_idx_offsets = nullptr;
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
    uint32_t max_allele_ct = 2;
    uint32_t max_allele_slen = 0;
    uint32_t max_filter_slen = 0;
    UnsortedVar vpos_sortstatus = kfUnsortedVar0;
    double* variant_cms = nullptr;
    ChrIdx* chr_idxs = nullptr;  // split-chromosome case only
    if (pvarname[0]) {
      char** pvar_filter_storage_mutable = nullptr;

      // LoadPvar() uses pvar_psam_flags to determine what's needed for .pvar
      // export.  These booleans are just for tracking requirements beyond
      // that.
      const uint32_t xheader_needed = (pcp->exportf_info.flags & (kfExportfVcf | kfExportfBcf))? 1 : 0;
      const uint32_t qualfilter_needed = xheader_needed || ((pcp->rmdup_mode != kRmDup0) && (pcp->rmdup_mode <= kRmDupExcludeMismatch));

      uint32_t neg_bp_seen = 0;
      reterr = LoadPvar(pvarname, pcp->var_filter_exceptions_flattened, pcp->varid_template_str, pcp->varid_multi_template_str, pcp->varid_multi_nonsnp_template_str, pcp->missing_varid_match, pcp->require_info_flattened, pcp->require_no_info_flattened, &(pcp->extract_if_info_expr), &(pcp->exclude_if_info_expr), pcp->misc_flags, pcp->pvar_psam_flags, pcp->load_filter_log_flags, xheader_needed, qualfilter_needed, pcp->var_min_qual, pcp->splitpar_bound1, pcp->splitpar_bound2, pcp->new_variant_id_max_allele_slen, (pcp->filter_flags / kfFilterSnpsOnlyJustAcgt) & 1, !(pcp->dependency_flags & kfFilterNoSplitChr), pcp->filter_min_allele_ct, pcp->filter_max_allele_ct, pcp->input_missing_geno_char, pcp->max_thread_ct, cip, &max_variant_id_slen, &info_reload_slen, &vpos_sortstatus, &xheader, &variant_include, &variant_bps, &variant_ids_mutable, &allele_idx_offsets, K_CAST(const char***, &allele_storage_mutable), &pvar_qual_present, &pvar_quals, &pvar_filter_present, &pvar_filter_npass, &pvar_filter_storage_mutable, &nonref_flags, &variant_cms, &chr_idxs, &raw_variant_ct, &variant_ct, &neg_bp_seen, &max_allele_ct, &max_allele_slen, &xheader_blen, &info_flags, &max_filter_slen);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
      LoadFilterLogFlags load_filter_log_flags = pcp->load_filter_log_flags;
      if (load_filter_log_flags & kfLoadFilterLogImportMergeAlreadyApplied) {
        load_filter_log_flags &= ~kfLoadFilterLogImportMergeMask;
      }
      if (unlikely(!variant_ct)) {
        // conditionally permit this?
        if (!raw_variant_ct) {
          logerrprintfww("Error: No variants in %s.\n", pvarname);
          goto Plink2Core_ret_MALFORMED_INPUT;
        }
        char* write_iter = strcpya_k(g_logbuf, "Error: All ");
        write_iter = u32toa(raw_variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " variant");
        if (raw_variant_ct != 1) {
          *write_iter++ = 's';
        }
        write_iter = strcpya_k(write_iter, " in ");
        write_iter = strcpya(write_iter, pvarname);
        write_iter = strcpya_k(write_iter, " excluded by ");
        if (load_filter_log_flags) {
          AppendLoadFilterFlagnames(load_filter_log_flags, &write_iter);
        }
        if (neg_bp_seen) {
          if (load_filter_log_flags) {
            write_iter = strcpya_k(write_iter, " and/or ");
          }
          write_iter = strcpya_k(write_iter, "negative bp coordinates");
        }
        strcpy_k(write_iter, ".\n");
        goto Plink2Core_ret_INCONSISTENT_INPUT_WW;
      }
      pvar_filter_storage = TO_CONSTCPCONSTP(pvar_filter_storage_mutable);
      {
        char* write_iter = u32toa(raw_variant_ct, g_logbuf);
        write_iter = strcpya_k(write_iter, " variant");
        if (raw_variant_ct != 1) {
          *write_iter++ = 's';
        }
        if (raw_variant_ct != variant_ct) {
          write_iter = strcpya_k(write_iter, " in ");
          write_iter = strcpya(write_iter, pvarname);
          write_iter = strcpya_k(write_iter, "; ");
          write_iter = u32toa(raw_variant_ct - variant_ct, write_iter);
          write_iter = strcpya_k(write_iter, " excluded by ");
          if (load_filter_log_flags) {
            AppendLoadFilterFlagnames(load_filter_log_flags, &write_iter);
          }
          if (neg_bp_seen) {
            if (load_filter_log_flags) {
              write_iter = strcpya_k(write_iter, " and/or ");
            }
            write_iter = strcpya_k(write_iter, "negative bp coordinates");
          }
          write_iter = strcpya_k(write_iter, ", ");
          write_iter = u32toa(variant_ct, write_iter);
          write_iter = strcpya_k(write_iter, " remaining");
        } else {
          write_iter = strcpya_k(write_iter, " loaded from ");
          write_iter = strcpya(write_iter, pvarname);
          if (load_filter_log_flags) {
            write_iter = strcpya_k(write_iter, " (");
            AppendLoadFilterFlagnames(load_filter_log_flags, &write_iter);
            write_iter = strcpya_k(write_iter, " had no effect)");
          }
        }
        strcpy_k(write_iter, ".\n");
        WordWrapB(0);
        logputsb();
      }
      if (info_reload_slen && (make_plink2_flags & (kfMakeBim | kfMakePvar)) && (!pvar_renamed)) {
        // need to be careful with .pvar in this case
        uint32_t fname_slen;
#ifdef _WIN32
        fname_slen = GetFullPathName(pvarname, kPglFnamesize, g_textbuf, nullptr);
        if (unlikely((!fname_slen) || (fname_slen > kPglFnamesize)))
#else
        if (unlikely(!realpath(pvarname, g_textbuf)))
#endif
        {
          logerrprintfww(kErrprintfFopen, pvarname, strerror(errno));
          goto Plink2Core_ret_OPEN_FAIL;
        }
        if (make_plink2_flags & kfMakeBim) {
          OutnameZstSet(".bim", make_plink2_flags & kfMakeBimZs, outname_end);
          pvar_renamed = RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
          if (pvar_renamed) {
            logerrputs("Warning: .bim input and output filenames match.  Appending '~' to input\nfilename.\n");
            fname_slen = strlen(pvarname);
            memcpy(g_textbuf, pvarname, fname_slen + 1);
            snprintf(&(pvarname[fname_slen]), 2, "~");
            if (unlikely(rename(g_textbuf, pvarname))) {
              logerrprintf("Error: Failed to append '~' to input .bim filename: %s.\n", strerror(errno));
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
            memcpy(g_textbuf, pvarname, fname_slen + 1);
            snprintf(&(pvarname[fname_slen]), 2, "~");
            if (unlikely(rename(g_textbuf, pvarname))) {
              logerrprintf("Error: Failed to append '~' to input .pvar filename: %s.\n", strerror(errno));
              goto Plink2Core_ret_OPEN_FAIL;
            }
          }
        }
      }
      if ((pcp->dependency_flags & (kfFilterAllReq | kfFilterNonrefFlagsNeeded)) == kfFilterNonrefFlagsNeeded) {
        // We need INFO/PR, but we don't really care about anything else
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
        reterr = PgfiInitPhase1(pgenname, pcp->pginame, raw_variant_ct, raw_sample_ct, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, g_logbuf);
        if (!reterr) {
          break;
        }
        // detect and autoconvert plink 1 sample-major files, instead of
        // failing (don't bother supporting plink 0.99 files any more)
        if (unlikely(reterr != kPglRetSampleMajorBed)) {
          WordWrapB(0);
          logerrputsb();
          goto Plink2Core_ret_1;
        }
        char* pgenname_end = memcpya(pgenname, outname, outname_end - outname);
        pgenname_end = strcpya_k(pgenname_end, ".pgen");
        const uint32_t no_vmaj_ext = (pcp->command_flags1 & kfCommand1MakePlink2) && (!pcp->filter_flags) && ((make_plink2_flags & (kfMakePgen | kfMakePlink2MJoin | (kfMakePgenFormatBase * 3) | kfMakePgenWriterVer)) == kfMakePgen);
        if (no_vmaj_ext) {
          *pgenname_end = '\0';
          make_plink2_flags &= ~kfMakePgen;
          // no --make-just-pgen command, so we'll never entirely skip the
          // MakePlink2 operation
        } else {
          snprintf(pgenname_end, kMaxOutfnameExtBlen - 5, ".vmaj");
        }
        reterr = Plink1SampleMajorToPgen(pgenname, nullptr, raw_variant_ct, raw_sample_ct, (pcp->misc_flags / kfMiscRealRefAlleles) & 1, pcp->max_thread_ct, pgfi.shared_ff);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        fclose(pgfi.shared_ff);
        pgfi.shared_ff = nullptr;
      }
      pgfi.allele_idx_offsets = allele_idx_offsets;
      pgfi.max_allele_ct = max_allele_ct;
      unsigned char* pgfi_alloc;
      if (unlikely(bigstack_alloc_uc(cur_alloc_cacheline_ct * kCacheline, &pgfi_alloc))) {
        goto Plink2Core_ret_NOMEM;
      }
      const uint32_t nonref_flags_already_loaded = (nonref_flags != nullptr);
      if (!nonref_flags_already_loaded) {
        const uint32_t nonref_flags_status_shifted = header_ctrl & 192;
        if (nonref_flags_status_shifted == 192) {
          if (unlikely(bigstack_alloc_w(raw_variant_ctl, &nonref_flags))) {
            goto Plink2Core_ret_NOMEM;
          }
        } else if ((pgfi.const_vrtype != kPglVrtypePlink1) && (!nonref_flags_status_shifted)) {
          // This should never happen with a chain of plink2-generated .pgen
          // files, but it's permitted by the specification and makes future
          // merge unsafe.
          logerrputs("Warning: .pgen file indicates that provisional-REF information is stored in the\ncompanion .pvar, but that .pvar lacks either an INFO/PR header or an INFO\ncolumn.  Assuming all REF alleles are correct.\n");
        }
      }
      pgfi.nonref_flags = nonref_flags;
      PgenExtensionLl* header_exts = nullptr;
      if (pcp->command_flags1 & kfCommand1PgenInfo) {
        ext_slot.next = nullptr;
        ext_slot.size = ~0LLU;
        ext_slot.type_idx = 1;
        header_exts = &ext_slot;
      }
      uint32_t max_vrec_width;
      // only practical effect of setting use_blockload to zero here is that
      // pgr_alloc_cacheline_ct is overestimated by
      // DivUp(max_vrec_width, kCacheline).
      reterr = PgfiInitPhase2Ex(header_ctrl, 1, nonref_flags_already_loaded, 1, 0, raw_variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, header_exts, nullptr, &pgr_alloc_cacheline_ct, g_logbuf);
      if (unlikely(reterr)) {
        WordWrapB(0);
        logerrputsb();
        goto Plink2Core_ret_1;
      }
      if (unlikely((!allele_idx_offsets) && (pgfi.gflags & kfPgenGlobalMultiallelicHardcallFound))) {
        logerrputs("Error: .pgen file contains multiallelic variants, while .pvar does not.\n");
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      if (pcp->misc_flags & kfMiscRealRefAlleles) {
        if (unlikely(nonref_flags && (!AllBitsAreOne(nonref_flags, raw_variant_ct)))) {
          // To reduce the ease of foot-shooting, we don't allow this to
          // clobber arbitrary nonref_flags arrays.
          logerrputs("Error: --real-ref-alleles must be used on a plink1 or similar fileset.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }

        // wasteful if nonref_flags was allocated, but nonref_flags isn't that
        // large, and --real-ref-alleles + --make-pgen can be run separately
        // from anything truly memory-limited
        nonref_flags = nullptr;
        pgfi.nonref_flags = nullptr;

        pgfi.gflags &= ~kfPgenGlobalAllNonref;
      }
      if (header_exts && (header_exts->size != ~0LLU)) {
#ifndef __LP64__
        if (unlikely(header_exts->size >= (1LLU << 31))) {
          goto Plink2Core_ret_NOMEM;
        }
#endif
        if (unlikely(pgl_malloc(header_exts->size, &(header_exts->contents)))) {
          goto Plink2Core_ret_NOMEM;
        }
        reterr = PgfiInitLoadExts(header_ctrl, &pgfi, header_exts, nullptr, g_logbuf);
        if (unlikely(reterr)) {
          WordWrapB(0);
          logerrputsb();
          goto Plink2Core_ret_1;
        }
      }
      if (SingleVariantLoaderIsNeeded(king_cutoff_fprefix, pcp->command_flags1, make_plink2_flags, pcp->rmdup_mode, pcp->hwe_ln_thresh)) {
        // ugly kludge, probably want to add pgenlib_internal support for this
        // hybrid use pattern
        FILE* shared_ff_copy = pgfi.shared_ff;
        pgfi.shared_ff = nullptr;
        unsigned char* simple_pgr_alloc;
        if (unlikely(bigstack_alloc_uc((pgr_alloc_cacheline_ct + DivUp(max_vrec_width, kCacheline)) * kCacheline, &simple_pgr_alloc))) {
          goto Plink2Core_ret_NOMEM;
        }
        reterr = PgrInit(pgenname, max_vrec_width, &pgfi, &simple_pgr, simple_pgr_alloc);
        if (unlikely(reterr)) {
          if (reterr == kPglRetOpenFail) {
            logerrprintfww(kErrprintfFopen, pgenname, strerror(errno));
          } else {
            assert(reterr == kPglRetReadFail);
            logerrprintfww(kErrprintfFread, pgenname, rstrerror(errno));
          }
          goto Plink2Core_ret_1;
        }
        pgfi.shared_ff = shared_ff_copy;
        if (pcp->command_flags1 & kfCommand1Validate) {
          uintptr_t* genovec_buf;
          if (unlikely(bigstack_alloc_w(NypCtToWordCt(raw_sample_ct), &genovec_buf))) {
            goto Plink2Core_ret_NOMEM;
          }
          logprintfww5("Validating %s... ", pgenname);
          fflush(stdout);
          reterr = PgrValidate(&simple_pgr, genovec_buf, g_logbuf);
          if (unlikely(reterr)) {
            logputs("\n");
            WordWrapB(0);
            logerrputsb();
            goto Plink2Core_ret_1;
          }
          logputs("done.\n");
          if (pcp->command_flags1 == kfCommand1Validate) {
            goto Plink2Core_ret_1;
          }
          BigstackReset(genovec_buf);
        }
      }
      // any functions using blockload must perform its own PgrInit(), etc.
      if (pcp->command_flags1 & kfCommand1PgenInfo) {
        if (pgfi.const_vrtype == kPglVrtypePlink1) {
          logerrputs("Warning: Skipping --pgen-info since a .bed file was provided.\n");
        } else {
          PgenInfoPrint(pgenname, &pgfi, header_exts, header_ctrl, max_allele_ct);
        }
        if (!(pcp->command_flags1 & (~(kfCommand1Validate | kfCommand1PgenInfo)))) {
          goto Plink2Core_ret_1;
        }
      }

    } else {
      // bugfix (10-11 Feb 2018): these variables may still be accessed
      pgfi.gflags = S_CAST(PgenGlobalFlags, ((info_flags / kfInfoPrNonrefDefault) & 1) * kfPgenGlobalAllNonref);
      pgfi.nonref_flags = nonref_flags;
      // dirty hack (3 Oct 2023): this is needed for --sort-vars
      // --make-just-pvar to work
      GET_PRIVATE(simple_pgr, m).fi.gflags = pgfi.gflags;
      GET_PRIVATE(simple_pgr, m).fi.nonref_flags = nonref_flags;
    }
    if (pcp->misc_flags & kfMiscMakeFoundersFirst) {
      reterr = MakeFounders(sample_include, raw_sample_ct, sample_ct, (pcp->misc_flags / kfMiscMakeFoundersRequire2Missing) & 1, &pii, founder_info);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }
    if (pcp->pheno_fname) {
      reterr = LoadPhenos(pcp->pheno_fname, &(pcp->pheno_range_list), sample_include, &pii.sii, pcp->missing_catname, raw_sample_ct, sample_ct, pcp->missing_pheno, (pcp->misc_flags / kfMiscAffection01) & 1, (pcp->misc_flags / kfMiscNoCategorical) & 1, (pcp->misc_flags / kfMiscPhenoIidOnly) & 1, (pcp->misc_flags / kfMiscPhenoColNums) & 1, (pcp->misc_flags / kfMiscNeg9PhenoReallyMissing) & 1, pcp->max_thread_ct, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }
    if (pcp->not_pheno_flattened) {
      reterr = IgnorePhenosOrCovars(pcp->not_pheno_flattened, 0, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }

    // move processing of PLINK 1.x cluster-loading/filtering flags here, since
    // they're now under the categorical-phenotype umbrella
    if ((pcp->misc_flags & kfMiscCatPhenoFamily) || pcp->within_fname) {
      reterr = Plink1ClusterImport(pcp->within_fname, pcp->catpheno_name, pcp->family_missing_catname, sample_include, pii.sii.sample_ids, pcp->missing_catname, raw_sample_ct, sample_ct, pii.sii.max_sample_id_blen, pcp->mwithin_val, pcp->max_thread_ct, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen);
      if (unlikely(reterr)) {
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
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
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
    if (variant_ct && AlleleUniquenessCheckIsNeeded(pcp->update_alleles_info.fname, pcp->read_freq_fname, pcp->ref_allele_flag, pcp->alt_allele_flag, pcp->command_flags1, pcp->rmdup_mode, pcp->fa_flags)) {
      reterr = CheckAlleleUniqueness(variant_include, cip, chr_idxs, variant_bps, TO_CONSTCPCONSTP(variant_ids_mutable), allele_idx_offsets, TO_CONSTCPCONSTP(allele_storage_mutable), variant_ct, max_allele_ct, pcp->max_thread_ct);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }
    // If something like --snps is combined with a position-based filter which
    // may remove some of the named variants, we need to apply --snps first.
    // Otherwise, it may be very advantageous to apply the position-based
    // filters before constructing the variant ID hash table.  So we split this
    // into two cases.
    // (Actually, hash table construction is now usually fast enough, even with
    // >80m variants, that this optimization isn't a big deal.  But it's
    // already implemented, so let's keep it around until/unless it creates
    // some sort of conflict.)

    // Additional minor optimization: --recover-var-ids doesn't actually need
    // the ID hash table; it's just positioned in this code block anyway
    // because it would be too weird for it to be in a different position than
    // --update-name in the order of operations.
    const uint32_t htable_needed_early = variant_ct && (pcp->varid_from || pcp->varid_to || pcp->varid_snp || pcp->varid_exclude_snp || pcp->snps_range_list.name_ct || pcp->exclude_snps_range_list.name_ct);
    const uint32_t full_variant_id_htable_needed = variant_ct && (htable_needed_early || pcp->update_map_flag || pcp->update_name_flag || pcp->update_alleles_info.fname || (pcp->rmdup_mode != kRmDup0) || pcp->extract_col_cond_info.params || (pcp->flip_info.fname && (!pcp->flip_info.subset_fname)));
    if (!full_variant_id_htable_needed) {
      reterr = ApplyVariantBpFilters(pcp->extract_fnames, pcp->extract_intersect_fnames, pcp->exclude_fnames, cip, variant_bps, pcp->from_bp, pcp->to_bp, pcp->bed_border_bp, raw_variant_ct, pcp->filter_flags, vpos_sortstatus, pcp->max_thread_ct, variant_include, &variant_ct);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }

    const uint32_t extract_exclude_by_id = (pcp->extract_fnames && (!(pcp->filter_flags & (kfFilterExtractBed0 | kfFilterExtractBed1)))) ||
      (pcp->extract_intersect_fnames && (!(pcp->filter_flags & (kfFilterExtractIntersectBed0 | kfFilterExtractIntersectBed1)))) ||
      (pcp->exclude_fnames && (!(pcp->filter_flags & (kfFilterExcludeBed0 | kfFilterExcludeBed1))));
    if (variant_ct && (full_variant_id_htable_needed || extract_exclude_by_id || pcp->recover_var_ids_fname)) {
      if (vpos_sortstatus & kfUnsortedVarBp) {
        if (unlikely(pcp->varid_from || pcp->varid_to)) {
          logerrputs("Error: --from/--to require a sorted .pvar/.bim.  Retry this command after using\n--make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(pcp->window_bp != -1)) {
          logerrputs("Error: --window requires a sorted .pvar/.bim.  Retry this command after using\n--make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(pcp->recover_var_ids_fname)) {
          logerrputs("Error: --recover-var-ids requires a sorted .pvar/.bim.  Retry this command\nafter using --make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
      }
      // don't bother with having different allow_dups vs. no allow_dups hash
      // table structures, just check specific IDs for duplication in the
      // no-duplicates-allowed cases
      unsigned char* bigstack_mark = g_bigstack_base;
      uint32_t* variant_id_htable = nullptr;
      uint32_t* htable_dup_base = nullptr;
      uint32_t dup_ct = 0;
      uint32_t variant_id_htable_size = 0;
      if ((!pcp->recover_var_ids_fname) || htable_needed_early) {
        reterr = AllocAndPopulateIdHtableMt(variant_include, TO_CONSTCPCONSTP(variant_ids_mutable), variant_ct, bigstack_left() / 8, pcp->max_thread_ct, &variant_id_htable, &htable_dup_base, &variant_id_htable_size, &dup_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->varid_from || pcp->varid_to) {
        reterr = FromToFlag(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->varid_from, pcp->varid_to, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, variant_include, cip, &variant_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->varid_snp) {
        reterr = SnpFlag(variant_bps, TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->varid_snp, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 0, pcp->window_bp, variant_include, cip, &variant_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->snps_range_list.name_ct) {
        reterr = SnpsFlag(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, &(pcp->snps_range_list), raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 0, variant_include, &variant_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->varid_exclude_snp) {
        reterr = SnpFlag(variant_bps, TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->varid_exclude_snp, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 1, pcp->window_bp, variant_include, cip, &variant_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->exclude_snps_range_list.name_ct) {
        reterr = SnpsFlag(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, &(pcp->exclude_snps_range_list), raw_variant_ct, max_variant_id_slen, variant_id_htable_size, 1, variant_include, &variant_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (variant_ct) {
        if (pcp->update_map_flag) {
          reterr = UpdateVarBps(cip, TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->update_map_flag, (pcp->sort_vars_mode > kSortNone), raw_variant_ct, max_variant_id_slen, variant_id_htable_size, pcp->max_thread_ct, variant_include, variant_bps, &variant_ct, &vpos_sortstatus);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        } else if (pcp->update_name_flag || pcp->recover_var_ids_fname) {
          if (pcp->update_name_flag) {
            reterr = UpdateVarNames(variant_include, variant_id_htable, htable_dup_base, pcp->update_name_flag, raw_variant_ct, variant_id_htable_size, pcp->max_thread_ct, variant_ids_mutable, &max_variant_id_slen);
          } else {
            reterr = RecoverVarIds(pcp->recover_var_ids_fname, variant_include, cip, variant_bps, allele_idx_offsets, TO_CONSTCPCONSTP(allele_storage_mutable), pcp->missing_varid_match, raw_variant_ct, variant_ct, pcp->recover_var_ids_flags, pcp->max_thread_ct, variant_ids_mutable, &max_variant_id_slen, outname, outname_end);
          }
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
          if (extract_exclude_by_id || pcp->update_alleles_info.fname || (pcp->rmdup_mode != kRmDup0) || (pcp->flip_info.fname && (!pcp->flip_info.subset_fname))) {
            // Must (re)construct the hash table in this case.
            BigstackReset(bigstack_mark);
            dup_ct = 0;
            reterr = AllocAndPopulateIdHtableMt(variant_include, TO_CONSTCPCONSTP(variant_ids_mutable), variant_ct, bigstack_left() / 8, pcp->max_thread_ct, &variant_id_htable, &htable_dup_base, &variant_id_htable_size, &dup_ct);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
          }
        }

        if (pcp->update_alleles_info.fname) {
          reterr = UpdateVarAlleles(variant_include, TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, allele_idx_offsets, &(pcp->update_alleles_info), raw_variant_ct, max_variant_id_slen, variant_id_htable_size, max_allele_ct, pcp->input_missing_geno_char, pcp->max_thread_ct, allele_storage_mutable, &max_allele_slen, outname, outname_end);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }

        if (pcp->flip_info.fname && (!pcp->flip_info.subset_fname)) {
          reterr = FlipAlleles(variant_include, TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, allele_idx_offsets, &(pcp->flip_info), raw_variant_ct, variant_ct, max_variant_id_slen, variant_id_htable_size, pcp->max_thread_ct, allele_storage_mutable);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }

        if (pcp->extract_fnames && (!(pcp->filter_flags & (kfFilterExtractBed0 | kfFilterExtractBed1)))) {
          reterr = TokenExtractExclude(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->extract_fnames, "extract", raw_variant_ct, max_variant_id_slen, variant_id_htable_size, kVfilterExtract, pcp->max_thread_ct, variant_include, &variant_ct);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
        if (pcp->extract_intersect_fnames && (!(pcp->filter_flags & (kfFilterExtractIntersectBed0 | kfFilterExtractIntersectBed1)))) {
          reterr = TokenExtractExclude(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->extract_intersect_fnames, "extract-intersect", raw_variant_ct, max_variant_id_slen, variant_id_htable_size, kVfilterExtractIntersect, pcp->max_thread_ct, variant_include, &variant_ct);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
        if (pcp->exclude_fnames && (!(pcp->filter_flags & (kfFilterExcludeBed0 | kfFilterExcludeBed1)))) {
          reterr = TokenExtractExclude(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->exclude_fnames, "exclude", raw_variant_ct, max_variant_id_slen, variant_id_htable_size, kVfilterExclude, pcp->max_thread_ct, variant_include, &variant_ct);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
        if (pcp->extract_col_cond_info.params) {
          reterr = ExtractColCond(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, &pcp->extract_col_cond_info, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, pcp->max_thread_ct, variant_include, &variant_ct);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
        if (pcp->rmdup_mode != kRmDup0) {
          reterr = RmDup(sample_include, cip, variant_bps, TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, allele_idx_offsets, TO_CONSTCPCONSTP(allele_storage_mutable), pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, info_reload_slen? pvarname : nullptr, variant_cms, pcp->missing_varid_match, raw_sample_ct, sample_ct, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, dup_ct, pcp->rmdup_mode, (pcp->command_flags1 / kfCommand1RmDupList) & 1, pcp->max_thread_ct, pgenname[0]? (&simple_pgr) : nullptr, variant_include, &variant_ct, outname, outname_end);
          if (reterr || (!(pcp->command_flags1 & (~(kfCommand1Validate | kfCommand1PgenInfo | kfCommand1RmDupList))))) {
            goto Plink2Core_ret_1;
          }
        }
      }

      BigstackReset(bigstack_mark);
      if (full_variant_id_htable_needed) {
        reterr = ApplyVariantBpFilters(pcp->extract_fnames, pcp->extract_intersect_fnames, pcp->exclude_fnames, cip, variant_bps, pcp->from_bp, pcp->to_bp, pcp->bed_border_bp, raw_variant_ct, pcp->filter_flags, vpos_sortstatus, pcp->max_thread_ct, variant_include, &variant_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      // todo: --attrib; although it isn't a "standard" file format, it can be
      // more convenient than forcing users to generate full-blown sites-only
      // VCF files, etc.
    }
    // variant_ids[] is fixed from this point on.
    const char* const* variant_ids = TO_CONSTCPCONSTP(variant_ids_mutable);

    if (pcp->allele_alphanum_flags != kfAlleleAlphanum0) {
      reterr = AlleleAlphanumUpdate(variant_include, variant_ids, allele_idx_offsets, variant_ct, pcp->allele_alphanum_flags, pcp->max_thread_ct, allele_storage_mutable);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }

    // SetRefalt1FromFile() can alter pointers-to-missing in
    // allele_storage[] so we can't make this a const char* const*, but nothing
    // alters the pointed-to strings after this point.
    const char** allele_storage = K_CAST(const char**, allele_storage_mutable);

    if (pcp->thin_keep_prob != 1.0) {
      RandomThinProb("thin", "variant", pcp->thin_keep_prob, raw_variant_ct, sfmtp, variant_include, &variant_ct);
    } else if (pcp->thin_keep_ct != UINT32_MAX) {
      reterr = RandomThinCt("thin-count", "variant", pcp->thin_keep_ct, raw_variant_ct, sfmtp, variant_include, &variant_ct);
      if (unlikely(reterr)) {
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
      if (pcp->update_sample_ids_fname) {
        if (update_sample_ids_empty) {
          logputs("--update-ids: 0 samples updated.\n");
        } else {
          reterr = UpdateSampleIds(pcp->update_sample_ids_fname, sample_include, raw_sample_ct, sample_ct, &pii.sii);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
      } else {
        if (pcp->update_parental_ids_fname) {
          if (update_parental_ids_empty) {
            logputs("--update-parents: 0 samples updated.\n");
          } else {
            reterr = UpdateSampleParents(pcp->update_parental_ids_fname, &pii.sii, sample_include, raw_sample_ct, sample_ct, pcp->max_thread_ct, &pii.parental_id_info, founder_info);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
            // bugfix (16 Jan 2021)
            pii.sii.flags |= kfSampleIdParentsPresent;
          }
        }
        // --update-parents goes here
        if (pcp->update_sex_info.fname) {
          reterr = UpdateSampleSexes(sample_include, &pii.sii, &(pcp->update_sex_info), raw_sample_ct, sample_ct, pcp->max_thread_ct, sex_nm, sex_male);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
      }
      if ((pii.sii.flags & kfSampleIdFidPresent) && ((pii.sii.flags & kfSampleIdNoIdHeaderIidOnly) || (pcp->grm_flags & kfGrmNoIdHeaderIidOnly))) {
        for (uint32_t sample_idx = 0; sample_idx != raw_sample_ct; ++sample_idx) {
          const char* cur_sample_id = &(pii.sii.sample_ids[sample_idx * pii.sii.max_sample_id_blen]);
          if (unlikely(!memequal_sk(cur_sample_id, "0\t"))) {
            logerrputs("Error: 'iid-only' modifier can only be used when FIDs are missing or all-0.\n");
            goto Plink2Core_ret_INCONSISTENT_INPUT;
          }
        }
      }
      if (pcp->keepfam_fnames) {
        reterr = KeepOrRemove(pcp->keepfam_fnames, &pii.sii, raw_sample_ct, kfKeepFam, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->keep_fnames) {
        reterr = KeepOrRemove(pcp->keep_fnames, &pii.sii, raw_sample_ct, kfKeep0, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->removefam_fnames) {
        reterr = KeepOrRemove(pcp->removefam_fnames, &pii.sii, raw_sample_ct, kfKeepRemove | kfKeepFam, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->remove_fnames) {
        reterr = KeepOrRemove(pcp->remove_fnames, &pii.sii, raw_sample_ct, kfKeepRemove, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->indv_str) {
        KeepOneId(pcp->indv_str, &pii.sii, raw_sample_ct, (pcp->misc_flags / kfMiscIidSid) & 1, sample_include, &sample_ct);
      }

      // todo: --attrib-indiv

      if (pcp->keep_col_match_fname) {
        reterr = KeepColMatch(pcp->keep_col_match_fname, &pii.sii, pcp->keep_col_match_flattened, pcp->keep_col_match_name, raw_sample_ct, pcp->keep_col_match_num, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->misc_flags & kfMiscRequirePheno) {
        reterr = RequirePheno(pheno_cols, pheno_names, pcp->require_pheno_flattened, raw_sample_ct, pheno_ct, max_pheno_name_blen, 0, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->filter_flags & (kfFilterExclFemales | kfFilterExclMales | kfFilterExclNosex)) {
        if (pcp->filter_flags & kfFilterExclFemales) {
          for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
            sample_include[widx] &= (~sex_nm[widx]) | sex_male[widx];
          }
        }
        if (pcp->filter_flags & kfFilterExclMales) {
          BitvecInvmask(sex_male, raw_sample_ctl, sample_include);
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
          BitvecInvmask(founder_info, raw_sample_ctl, sample_include);
        }
        const uint32_t old_sample_ct = sample_ct;
        sample_ct = PopcountWords(sample_include, raw_sample_ctl);
        const uint32_t removed_ct = old_sample_ct - sample_ct;
        logprintf("--keep-%sfounders: %u sample%s removed.\n", keep_founders? "" : "non", removed_ct, (removed_ct == 1)? "" : "s");
      }

      if (pcp->thin_keep_sample_prob != 1.0) {
        RandomThinProb("thin-indiv", "sample", pcp->thin_keep_sample_prob, raw_sample_ct, sfmtp, sample_include, &sample_ct);
      } else if (pcp->thin_keep_sample_ct != UINT32_MAX) {
        reterr = RandomThinCt("thin-indiv-count", "sample", pcp->thin_keep_sample_ct, raw_sample_ct, sfmtp, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      const uint32_t smaj_missing_geno_report_requested = (pcp->command_flags1 & kfCommand1MissingReport) && (!(pcp->missing_rpt_flags & kfMissingRptVariantOnly));
      if ((pcp->mind_thresh < 1.0) || (pcp->select_sid_missingness_mode != kSelectSidMissingness0) || smaj_missing_geno_report_requested) {
        if (unlikely(bigstack_alloc_u32(raw_sample_ct, &sample_missing_hc_cts) ||
                     bigstack_alloc_u32(raw_sample_ct, &sample_hethap_cts))) {
          goto Plink2Core_ret_NOMEM;
        }
        if (SampleMissingDosageCtsAreNeeded(pcp->misc_flags, smaj_missing_geno_report_requested, pcp->mind_thresh, pcp->select_sid_missingness_mode, pcp->missing_rpt_flags)) {
          if (pgfi.gflags & kfPgenGlobalDosagePresent) {
            if (unlikely(bigstack_alloc_u32(raw_sample_ct, &sample_missing_dosage_cts))) {
              goto Plink2Core_ret_NOMEM;
            }
          } else {
            sample_missing_dosage_cts = sample_missing_hc_cts;
          }
        }
        // could avoid this call and make LoadAlleleAndGenoCounts() do
        // double duty with --missing?
        reterr = LoadSampleMissingCts(sample_include, sex_nm, sex_male, variant_include, cip, "sample missingness", raw_variant_ct, variant_ct, raw_sample_ct, (pcp->misc_flags / kfMiscYNosexMissingStats) & 1, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, sample_missing_hc_cts, (pgfi.gflags & kfPgenGlobalDosagePresent)? sample_missing_dosage_cts : nullptr, sample_hethap_cts);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        if (pcp->mind_thresh < 1.0) {
          uint32_t variant_ct_y = 0;
          uint32_t y_code;
          if (XymtExists(cip, kChrOffsetY, &y_code)) {
            variant_ct_y = CountChrVariantsUnsafe(variant_include, cip, y_code);
          }
          reterr = MindFilter((pcp->misc_flags & kfMiscMindDosage)? sample_missing_dosage_cts : sample_missing_hc_cts, (pcp->misc_flags & kfMiscMindHhMissing)? sample_hethap_cts : nullptr, &pii.sii, raw_sample_ct, variant_ct, variant_ct_y, pcp->mind_thresh, sample_include, sex_male, &sample_ct, outname, outname_end);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
      }
      if (pcp->filter_flags & kfFilterSelectSidRepresentatives) {
        reterr = SelectSidRepresentatives(sex_nm, sex_male, (pcp->select_sid_missingness_mode == kSelectSidMissingnessDosage)? sample_missing_dosage_cts : sample_missing_hc_cts, (pcp->select_sid_missingness_mode == kSelectSidMissingnessHhMissing)? sample_hethap_cts : nullptr, &pii, pcp->select_sid_missingness_mode, pcp->select_sid_tiebreak_mode, (pcp->misc_flags / kfMiscSelectSidParentsOnly) & 1, raw_sample_ct, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (sample_missing_hc_cts && (!smaj_missing_geno_report_requested)) {
        BigstackReset(sample_missing_hc_cts);
      }
      // this results in a small "memory leak" when a regular missingness
      // report is requested, not a big deal
      if (pcp->misc_flags & kfMiscMakeFoundersNotfirst) {
        reterr = MakeFounders(sample_include, raw_sample_ct, sample_ct, (pcp->misc_flags / kfMiscMakeFoundersRequire2Missing) & 1, &pii, founder_info);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->covar_fname || pcp->covar_range_list.name_ct) {
        const char* cur_covar_fname = pcp->covar_fname? pcp->covar_fname : (pcp->pheno_fname? pcp->pheno_fname : psamname);
        reterr = LoadPhenos(cur_covar_fname, &(pcp->covar_range_list), sample_include, &pii.sii, pcp->missing_catname, raw_sample_ct, sample_ct, pcp->missing_pheno, 2, (pcp->misc_flags / kfMiscNoCategorical) & 1, (pcp->misc_flags / kfMiscCovarIidOnly) & 1, (pcp->misc_flags / kfMiscCovarColNums) & 1, (pcp->misc_flags / kfMiscNeg9PhenoReallyMissing) & 1, pcp->max_thread_ct, &covar_cols, &covar_names, &covar_ct, &max_covar_name_blen);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        logprintf("%u covariate%s loaded from %s.\n", covar_ct, (covar_ct == 1)? "" : "s", cur_covar_fname);

        // do we still want to clear some main phenotype values here if some
        // covariate values are missing?  (don't think there's a point to
        // preserving that behavior, let the regression functions do it to
        // their local phenotype copies on their own.)
      }
      if (pcp->not_covar_flattened) {
        reterr = IgnorePhenosOrCovars(pcp->not_covar_flattened, 1, &covar_cols, &covar_names, &covar_ct, &max_covar_name_blen);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }


      if (pcp->misc_flags & kfMiscRequireCovar) {
        reterr = RequirePheno(covar_cols, covar_names, pcp->require_covar_flattened, raw_sample_ct, covar_ct, max_covar_name_blen, 1, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->keep_if_expr.etype != kCmpExprTypeNull) {
        reterr = KeepRemoveIf(&(pcp->keep_if_expr), pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, (pcp->misc_flags / kfMiscAffection01) & 1, 0, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->remove_if_expr.etype != kCmpExprTypeNull) {
        reterr = KeepRemoveIf(&(pcp->remove_if_expr), pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, (pcp->misc_flags / kfMiscAffection01) & 1, 1, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      // meow
      if (pcp->keep_cats_fname || pcp->keep_cat_names_flattened) {
        reterr = KeepRemoveCats(pcp->keep_cats_fname, pcp->keep_cat_names_flattened, pcp->keep_cat_phenoname, pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, 0, pcp->max_thread_ct, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->remove_cats_fname || pcp->remove_cat_names_flattened) {
        reterr = KeepRemoveCats(pcp->remove_cats_fname, pcp->remove_cat_names_flattened, pcp->remove_cat_phenoname, pheno_cols, pheno_names, covar_cols, covar_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, 1, pcp->max_thread_ct, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
    }

    const uint32_t nonfounders = (pcp->misc_flags / kfMiscNonfounders) & 1;
    uint32_t founder_ct = 0;
    uint32_t male_ct = 0;
    uint32_t nosex_ct = 0;
    if (psamname[0]) {
      if (unlikely(!sample_ct)) {
        logerrputs("Error: No samples remaining after main filters.\n");
        goto Plink2Core_ret_DEGENERATE_DATA;
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
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }

    // quantnorm before variance-standardize, since at least that has a minor
    // effect, whereas the other order is pointless
    if (pcp->pheno_transform_flags & (kfPhenoTransformQuantnormPheno | kfPhenoTransformQuantnormAll)) {
      reterr = PhenoQuantileNormalize(pcp->quantnorm_flattened, sample_include, pheno_names, raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, 0, (pcp->pheno_transform_flags / kfPhenoTransformQuantnormPheno) & 1, pheno_cols);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }
    if (pcp->pheno_transform_flags & (kfPhenoTransformQuantnormCovar | kfPhenoTransformQuantnormAll)) {
      reterr = PhenoQuantileNormalize((pcp->pheno_transform_flags & kfPhenoTransformQuantnormAll)? pcp->quantnorm_flattened : pcp->covar_quantnorm_flattened, sample_include, covar_names, raw_sample_ct, sample_ct, covar_ct, max_covar_name_blen, 1, (pcp->pheno_transform_flags / kfPhenoTransformQuantnormCovar) & 1, covar_cols);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }

    if (pcp->pheno_transform_flags & (kfPhenoTransformVstdCovar | kfPhenoTransformVstdAll)) {
      const uint32_t is_covar_flag = (pcp->pheno_transform_flags / kfPhenoTransformVstdCovar) & 1;
      if (!is_covar_flag) {
        reterr = PhenoVarianceStandardize(pcp->vstd_flattened, sample_include, pheno_names, raw_sample_ct, pheno_ct, max_pheno_name_blen, 0, 0, pheno_cols);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      reterr = PhenoVarianceStandardize(pcp->vstd_flattened, sample_include, covar_names, raw_sample_ct, covar_ct, max_covar_name_blen, 1, is_covar_flag, covar_cols);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }

#ifndef NOLAPACK
    if (pcp->command_flags1 & kfCommand1PhenoSvd) {
      reterr = PhenoSvd(&pcp->pheno_svd_info, sample_include, &pii.sii, raw_sample_ct, sample_ct, pcp->max_thread_ct, &pheno_names, &pheno_ct, &max_pheno_name_blen, pheno_cols, outname, outname_end);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }
#endif

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
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          if (memequal(loop_cats_phenoname, &(pheno_names[pheno_idx * max_pheno_name_blen]), name_blen)) {
            PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
            if (unlikely(cur_pheno_col->type_code != kPhenoDtypeCat)) {
              logerrprintfww("Error: '%s' is not a categorical phenotype.\n", loop_cats_phenoname);
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            // bugfix (13 Dec 2019): we need to copy the contents of
            // pheno_cols[pheno_idx] when pheno_idx != pheno_ct - 1.
            // Also note that this makes pheno_ct == 0 and pheno_cols !=
            // nullptr possible simultaneously.
            --pheno_ct;
            if (pheno_idx < pheno_ct) {
              PhenoCol pheno_col_copy = *cur_pheno_col;
              for (uint32_t pheno_idx2 = pheno_idx; pheno_idx2 != pheno_ct; ++pheno_idx2) {
                pheno_cols[pheno_idx2] = pheno_cols[pheno_idx2 + 1];
                memcpy(&(pheno_names[pheno_idx2 * max_pheno_name_blen]), &(pheno_names[(pheno_idx2 + 1) * max_pheno_name_blen]), max_pheno_name_blen);
              }
              pheno_cols[pheno_ct] = pheno_col_copy;
            }
            loop_cats_pheno_col = &(pheno_cols[pheno_ct]);
            break;
          }
        }
      }
      if ((!loop_cats_pheno_col) && (name_blen <= max_covar_name_blen)) {
        for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
          if (memequal(loop_cats_phenoname, &(covar_names[covar_idx * max_covar_name_blen]), name_blen)) {
            PhenoCol* cur_covar_col = &(covar_cols[covar_idx]);
            if (unlikely(cur_covar_col->type_code != kPhenoDtypeCat)) {
              logerrprintfww("Error: '%s' is not a categorical covariate.\n", loop_cats_phenoname);
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            --covar_ct;
            if (covar_idx < covar_ct) {
              PhenoCol covar_col_copy = *cur_covar_col;
              for (uint32_t covar_idx2 = covar_idx; covar_idx2 != covar_ct; ++covar_idx2) {
                covar_cols[covar_idx2] = covar_cols[covar_idx2 + 1];
                memcpy(&(covar_names[covar_idx2 * max_covar_name_blen]), &(covar_names[(covar_idx2 + 1) * max_covar_name_blen]), max_covar_name_blen);
              }
              covar_cols[covar_ct] = covar_col_copy;
            }
            loop_cats_pheno_col = &(covar_cols[covar_ct]);
            break;
          }
        }
      }
      if (unlikely(!loop_cats_pheno_col)) {
        logerrprintfww("Error: --loop-cats phenotype '%s' not loaded.\n", loop_cats_phenoname);
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &loop_cats_sample_include_backup) ||
                   bigstack_alloc_w(raw_sample_ctl, &loop_cats_founder_info_backup) ||
                   bigstack_alloc_w(raw_sample_ctl, &loop_cats_sex_nm_backup) ||
                   bigstack_alloc_w(raw_sample_ctl, &loop_cats_sex_male_backup) ||
                   bigstack_alloc_w(1 + (loop_cats_pheno_col->nonnull_category_ct / kBitsPerWord), &loop_cats_cat_include))) {
        goto Plink2Core_ret_NOMEM;
      }
      if (variant_ct != raw_variant_ct) {
        if (unlikely(bigstack_alloc_w(raw_variant_ctl, &loop_cats_variant_include_backup))) {
          goto Plink2Core_ret_NOMEM;
        }
        memcpy(loop_cats_variant_include_backup, variant_include, raw_variant_ctl * sizeof(intptr_t));
      }
      BitvecAndCopy(sample_include, loop_cats_pheno_col->nonmiss, raw_sample_ctl, loop_cats_sample_include_backup);
      loop_cats_sample_ct = PopcountWords(loop_cats_sample_include_backup, raw_sample_ctl);
      loop_cats_ct = IdentifyRemainingCats(sample_include, loop_cats_pheno_col, sample_ct, loop_cats_cat_include);
      // bugfix (13 Dec 2019): --loop-cats does not iterate over the null
      // category.  Don't include it in the count.
      loop_cats_ct -= loop_cats_cat_include[0] & 1;
      if (unlikely(!loop_cats_ct)) {
        logerrputs("Error: All --loop-cats categories are empty.\n");
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      logprintf("--loop-cats: %u categor%s present.\n", loop_cats_ct, (loop_cats_ct == 1)? "y" : "ies");
      BitvecAndCopy(sample_include, founder_info, raw_sample_ctl, loop_cats_founder_info_backup);
      BitvecAndCopy(sample_include, sex_nm, raw_sample_ctl, loop_cats_sex_nm_backup);
      BitvecAndCopy(sample_include, sex_male, raw_sample_ctl, loop_cats_sex_male_backup);
      *outname_end = '.';
    }
    const uintptr_t raw_allele_ct = allele_idx_offsets? allele_idx_offsets[raw_variant_ct] : (2 * raw_variant_ct);
    unsigned char* bigstack_mark_varfilter = g_bigstack_base;
    unsigned char* bigstack_end_mark_varfilter = g_bigstack_end;
    for (uint32_t loop_cats_idx = 0; loop_cats_idx != loop_cats_ct; ++loop_cats_idx) {
      BigstackDoubleReset(bigstack_mark_varfilter, bigstack_end_mark_varfilter);
      if (loop_cats_pheno_col) {
        loop_cats_uidx = AdvTo1Bit(loop_cats_cat_include, loop_cats_uidx + 1);
        const char* catname = loop_cats_pheno_col->category_names[loop_cats_uidx];
        const uint32_t catname_slen = strlen(catname);
        if (unlikely(catname_slen + S_CAST(uintptr_t, loop_cats_outname_endp1_backup - outname) > (kPglFnamesize - kMaxOutfnameExtBlen))) {
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

      uintptr_t* allele_presents = nullptr;

      // dosages are currently in 32768ths, hence the extra 'd'
      // same indexes as allele_storage.  We can't omit the last allele (in the
      // way allele_freqs does) because the sum isn't constant (missing
      // values).
      uint64_t* allele_ddosages = nullptr;
      uint64_t* founder_allele_ddosages = nullptr;

      AlleleCode* maj_alleles = nullptr;
      double* allele_freqs = nullptr;
      STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts) = nullptr;
      STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts) = nullptr;
      unsigned char* bigstack_mark_allele_ddosages = g_bigstack_base;
      const uint32_t keep_grm = GrmKeepIsNeeded(pcp->command_flags1, pcp->pca_flags);
      double* grm = nullptr;

      if (pcp->command_flags1 & kfCommand1WriteSamples) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".id");
        reterr = WriteSampleIds(sample_include, &pii.sii, outname, sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        logprintfww("--write-samples: Sample IDs written to %s .\n", outname);
        if (!(pcp->command_flags1 & (~(kfCommand1WriteSamples | kfCommand1Validate | kfCommand1PgenInfo | kfCommand1RmDupList)))) {
          continue;
        }
      }

      if (pgenname[0]) {
        if (unlikely((pcp->command_flags1 & (kfCommand1LdPrune | kfCommand1Ld)) && (founder_ct < 50) && (!(pcp->misc_flags & kfMiscAllowBadLd)))) {
          if (sample_ct < 50) {
            logerrputs("Error: This run estimates linkage disequilibrium between variants, but there\nare less than 50 samples to estimate from.  You should perform this operation\non a larger dataset.\n(Strictly speaking, you can also override this error with --bad-ld, but this is\nalmost always a bad idea.)\n");
          } else {
            logerrputs("Error: This run estimates linkage disequilibrium between variants, but there\nare less than 50 founders to estimate from.  --make-founders may help.\n(Strictly speaking, you can also override this error with --bad-ld, but this is\nalmost always a bad idea.)\n");
          }
          goto Plink2Core_ret_DEGENERATE_DATA;
        }
        const uint32_t decent_afreqs_needed = DecentAlleleFreqsAreNeeded(pcp->command_flags1, pcp->check_sex_info.flags, pcp->het_flags, pcp->score_info.flags);
        const uint32_t maj_alleles_needed = MajAllelesAreNeeded(pcp->command_flags1, pcp->pca_flags, pcp->glm_info.flags, pcp->vcor_info.flags);
        if (decent_afreqs_needed || maj_alleles_needed || IndecentAlleleFreqsAreNeeded(pcp->command_flags1, pcp->vcor_info.flags, pcp->min_maf, pcp->max_maf)) {
          if (unlikely((!pcp->read_freq_fname) && ((sample_ct < 50) || ((!nonfounders) && (founder_ct < 50))) && decent_afreqs_needed && (!(pcp->misc_flags & kfMiscAllowBadFreqs)))) {
            if ((!nonfounders) && (sample_ct >= 50)) {
              logerrputs("Error: This run requires decent allele frequencies, but they aren't being\nloaded with --read-freq, and less than 50 founders are available to impute them\nfrom.  Possible solutions:\n* You can use --nonfounders to include nonfounders when imputing allele\n  frequencies.\n* You can generate (with --freq) or obtain an allele frequency file based on a\n  larger similar-population reference dataset, and load it with --read-freq.\n* (Not recommended) You can override this error with --bad-freqs.\n");
            } else {
              logerrputs("Error: This run requires decent allele frequencies, but they aren't being\nloaded with --read-freq, and less than 50 samples are available to impute them\nfrom.\nYou should generate (with --freq) or obtain an allele frequency file based on a\nlarger similar-population reference dataset, and load it with --read-freq.\n");
              if (nonfounders) {
                logerrputs("If you're certain you want to proceed without doing that, use --bad-freqs to\noverride this error.\n");
              } else if (sample_ct != founder_ct) {
                logerrputs("If you're certain you want to proceed without doing that, use --bad-freqs to\noverride this error, and consider using --nonfounders as well.\n");
              }
            }
            goto Plink2Core_ret_DEGENERATE_DATA;
          }
          if (maj_alleles_needed) {
            if (unlikely(bigstack_alloc_ac(raw_variant_ct, &maj_alleles))) {
              goto Plink2Core_ret_NOMEM;
            }
          }
          //   allele_freqs[allele_idx_offsets[variant_uidx] - variant_uidx]
          // stores the frequency estimate for the reference allele; if there's
          // more than 1 alt allele, next element stores alt1 freq, etc.  To
          // save memory, we omit the last alt.
          if (unlikely(bigstack_alloc_d(raw_allele_ct - raw_variant_ct, &allele_freqs))) {
            goto Plink2Core_ret_NOMEM;
          }
        }
        if ((sample_ct != founder_ct) && (pcp->min_allele_ddosage || (pcp->max_allele_ddosage != (~0LLU)) || ((!pcp->read_freq_fname) && (pcp->freq_rpt_flags & kfAlleleFreqCounts))) && (!nonfounders) && (!(pcp->misc_flags & kfMiscAcFounders))) {
          logerrputs("Error: --mac/--max-mac/\"--freq counts\" specified, but with neither\n--ac-founders nor --nonfounders; and nonfounders are present.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        uint32_t x_start = 0;
        uint32_t x_len = 0;
        uint32_t hwe_x_probs_needed = 0;
        uint32_t x_code;
        if (XymtExists(cip, kChrOffsetX, &x_code) && (!(vpos_sortstatus & kfUnsortedVarSplitChr))) {
          const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
          x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
          const uint32_t x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
          x_len = x_end - x_start;
          if (x_len && ((pcp->command_flags1 & kfCommand1Hardy) || (pcp->hwe_ln_thresh != -DBL_MAX)) && (!AllBitsAreZero(variant_include, x_start, x_end))) {
            if (nonfounders) {
              hwe_x_probs_needed = (sample_ct > nosex_ct);
            } else {
              // at least one founder with known gender?
              hwe_x_probs_needed = !IntersectionIsEmpty(founder_info, sex_nm, raw_sample_ctl);
            }
          }
        }
        bigstack_mark_allele_ddosages = g_bigstack_base;
        const uint32_t first_hap_uidx = GetFirstHaploidUidx(cip, vpos_sortstatus);
        if (make_plink2_flags & kfMakePlink2TrimAlts) {
          if (unlikely(bigstack_alloc_w(BitCtToWordCt(raw_allele_ct), &allele_presents))) {
            goto Plink2Core_ret_NOMEM;
          }
        }
        uint32_t regular_freqcounts_needed = (allele_presents != nullptr);
        if (AlleleDosagesAreNeeded(pcp->command_flags1, pcp->misc_flags, (allele_freqs != nullptr), pcp->min_allele_ddosage, pcp->max_allele_ddosage, &regular_freqcounts_needed)) {
          if (unlikely(bigstack_alloc_u64(raw_allele_ct, &allele_ddosages))) {
            goto Plink2Core_ret_NOMEM;
          }
        }
        uint32_t overlapping_allele_ddosages = 0;
        if (FounderAlleleDosagesAreNeeded(pcp->command_flags1, pcp->misc_flags, (allele_freqs != nullptr), pcp->min_allele_ddosage, pcp->max_allele_ddosage, &regular_freqcounts_needed)) {
          if (allele_ddosages && (founder_ct == sample_ct)) {
            overlapping_allele_ddosages = 1;
          } else {
            if (unlikely(bigstack_alloc_u64(raw_allele_ct, &founder_allele_ddosages))) {
              goto Plink2Core_ret_NOMEM;
            }
          }
        }
        double* imp_r2_vals = nullptr;
        const uint32_t is_minimac3_r2 = (pcp->freq_rpt_flags & kfAlleleFreqColMinimac3R2) || (pcp->minimac3_r2_max != 0.0);
        if (is_minimac3_r2 || (pcp->freq_rpt_flags & kfAlleleFreqColMachR2) || (pcp->mach_r2_max != 0.0)) {
          if (unlikely(bigstack_alloc_d(raw_variant_ct, &imp_r2_vals))) {
            goto Plink2Core_ret_NOMEM;
          }
        }

        unsigned char* bigstack_mark_geno_cts = g_bigstack_base;

        // no longer includes hethaps by default
        uint32_t* variant_missing_hc_cts = nullptr;
        uint32_t* variant_hethap_cts = nullptr;
        if (VariantMissingHcCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->geno_thresh, pcp->missing_rpt_flags)) {
          if (unlikely(bigstack_alloc_u32(raw_variant_ct, &variant_missing_hc_cts))) {
            goto Plink2Core_ret_NOMEM;
          }
          if (VariantHethapCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->geno_thresh, pcp->missing_rpt_flags, first_hap_uidx)) {
            // first_hap_uidx offset can save an entire GB...
            if (unlikely(bigstack_alloc_u32(raw_variant_ct - first_hap_uidx, &variant_hethap_cts))) {
              goto Plink2Core_ret_NOMEM;
            }
          }
        }
        uint32_t* variant_missing_dosage_cts = nullptr;
        if (VariantMissingDosageCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->geno_thresh, pcp->missing_rpt_flags)) {
          if ((!variant_missing_hc_cts) || (pgfi.gflags & kfPgenGlobalDosagePresent)) {
            if (unlikely(bigstack_alloc_u32(raw_variant_ct, &variant_missing_dosage_cts))) {
              goto Plink2Core_ret_NOMEM;
            }
          } else {
            variant_missing_dosage_cts = variant_missing_hc_cts;
          }
        }
        STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts) = nullptr;
        STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts) = nullptr;
        STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts) = nullptr;
        STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts) = nullptr;
        // [0] = homref ct, [1] = het ref-altx total, [2] = nonref diploid
        //   total
        // use unfiltered indexes, since we remove more variants later
        if (RawGenoCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->hwe_ln_thresh)) {
          if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, raw_variant_ct, &raw_geno_cts))) {
            goto Plink2Core_ret_NOMEM;
          }
          if (x_len) {
            if (male_ct) {
              if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, x_len, &x_male_geno_cts))) {
                goto Plink2Core_ret_NOMEM;
              }
            }
            if (nosex_ct && hwe_x_probs_needed && nonfounders) {
              if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, x_len, &x_nosex_geno_cts))) {
                goto Plink2Core_ret_NOMEM;
              }
            }
          }
        }
        if (FounderRawGenoCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->hwe_ln_thresh)) {
          if ((founder_ct == sample_ct) && raw_geno_cts) {
            founder_raw_geno_cts = raw_geno_cts;
            founder_x_male_geno_cts = x_male_geno_cts;
          } else {
            if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, raw_variant_ct, &founder_raw_geno_cts))) {
              goto Plink2Core_ret_NOMEM;
            }
            if (x_len && male_ct) {
              if (!IntersectionIsEmpty(founder_info, sex_male, raw_sample_ctl)) {
                if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, x_len, &founder_x_male_geno_cts))) {
                  goto Plink2Core_ret_NOMEM;
                }
              }
            }
          }
          if (nosex_ct && hwe_x_probs_needed && (!nonfounders)) {
            const uint32_t founder_knownsex_ct = PopcountWordsIntersect(founder_info, sex_nm, raw_sample_ctl);
            if (founder_knownsex_ct < founder_ct) {
              if ((founder_ct == sample_ct) && x_nosex_geno_cts) {
                // shouldn't be possible for now, since x_nosex_geno_cts can
                // only be non-null when nonfounders is true
                assert(0);
                // founder_x_nosex_geno_cts = x_nosex_geno_cts;
              } else {
                if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, x_len, &founder_x_nosex_geno_cts))) {
                  goto Plink2Core_ret_NOMEM;
                }
              }
            }
          }
        }
        unsigned char* bigstack_mark_read_freqs = g_bigstack_base;
        regular_freqcounts_needed = regular_freqcounts_needed || variant_missing_hc_cts || variant_missing_dosage_cts || variant_hethap_cts || raw_geno_cts || founder_raw_geno_cts || imp_r2_vals;
        uintptr_t* variant_afreqcalc = nullptr;
        uint32_t afreqcalc_variant_ct = 0;
        if (allele_freqs) {
          if (pcp->read_freq_fname) {
            reterr = ReadAlleleFreqs(variant_include, variant_ids, allele_idx_offsets, allele_storage, pcp->read_freq_fname, raw_variant_ct, variant_ct, max_allele_ct, max_variant_id_slen, max_allele_slen, pcp->af_pseudocount, pcp->max_thread_ct, allele_freqs, &variant_afreqcalc);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
            if (variant_ct) {
              afreqcalc_variant_ct = PopcountWords(variant_afreqcalc, raw_variant_ctl);
            }
          } else {
            variant_afreqcalc = variant_include;
            afreqcalc_variant_ct = variant_ct;
          }
        } else if (pcp->read_freq_fname) {
          logerrputs("Warning: Ignoring --read-freq since no command would use the frequencies.\n");
        }
        if (regular_freqcounts_needed || afreqcalc_variant_ct) {
          if (unlikely(pcp->misc_flags & kfMiscErrorOnFreqCalc)) {
            // This needs to be kept in sync with DecentAlleleFreqsAreNeeded(),
            // MajAllelesAreNeeded(), IndecentAlleleFreqsAreNeeded(),
            // [Founder]AlleleDosagesAreNeeded(),
            // VariantMissingHcCtsAreNeeded(),
            // VariantMissingDosageCtsAreNeeded(),
            // [Founder]RawGenoCtsAreNeeded(), TrimAlts, and is_minimac3_r2.
            logerrputs("Error: --error-on-freq-calc specified, but allele frequency calculation is\nneeded.\nFlags which may invoke the allele frequency calculation include --check-sex,\n--freq, --geno, --geno-counts, --genotyping-rate, --glm (unless 'omit-ref' is\nspecified), --hardy, --het (unless 'small-sample' is specified), --hwe,\n--impute-sex, --indep-pairwise, the 'trim-alts' modifier of\n--make-[b]pgen/--make-bed, --make-grm-{bin,list}, --make-rel, --[max-]mac,\n--[max-]maf, --minimac3-r2-filter, --missing, --pca, --score[-list] (unless\n'no-mean-imputation' is specified, and neither 'center' nor\n'variance-standardize' are), and --variant-score.\n");
            goto Plink2Core_ret_INVALID_CMDLINE;
          }
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
          reterr = LoadAlleleAndGenoCounts(sample_include, founder_info, sex_nm, sex_male, regular_freqcounts_needed? variant_include : variant_afreqcalc, cip, allele_idx_offsets, raw_sample_ct, sample_ct, founder_ct, male_ct, nosex_ct, raw_variant_ct, regular_freqcounts_needed? variant_ct : afreqcalc_variant_ct, first_hap_uidx, is_minimac3_r2, (pcp->misc_flags / kfMiscYNosexMissingStats) & 1, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, allele_presents, allele_ddosages, founder_allele_ddosages, ((!variant_missing_hc_cts) && dosageless_file)? variant_missing_dosage_cts : variant_missing_hc_cts, dosageless_file? nullptr : variant_missing_dosage_cts, variant_hethap_cts, raw_geno_cts, founder_raw_geno_cts, x_male_geno_cts, founder_x_male_geno_cts, x_nosex_geno_cts, founder_x_nosex_geno_cts, imp_r2_vals);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
          if (overlapping_allele_ddosages) {
            founder_allele_ddosages = allele_ddosages;
          }
          if (pcp->command_flags1 & kfCommand1GenotypingRate) {
            // possible todo: also report this opportunistically
            // (variant_missing_hc_cts filled for other reasons).  worth
            // multithreading in that case.
            const uint32_t is_dosage = (pcp->misc_flags / kfMiscGenotypingRateDosage) & 1;
            const uint32_t y_sample_ct = male_ct + ((pcp->misc_flags / kfMiscYNosexMissingStats) & 1) * nosex_ct;
            ReportGenotypingRate(variant_include, cip, is_dosage? variant_missing_dosage_cts : variant_missing_hc_cts, raw_sample_ct, sample_ct, y_sample_ct, variant_ct, is_dosage);
            if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1PhenoSvd | kfCommand1WriteSamples | kfCommand1Validate | kfCommand1PgenInfo | kfCommand1RmDupList)))) {
              continue;
            }
          }
        }
        if (allele_freqs) {
          if (afreqcalc_variant_ct) {
            ComputeAlleleFreqs(variant_afreqcalc, allele_idx_offsets, nonfounders? allele_ddosages : founder_allele_ddosages, afreqcalc_variant_ct, pcp->af_pseudocount, allele_freqs);
          }
          if (maj_alleles) {
            ComputeMajAlleles(variant_include, allele_idx_offsets, allele_freqs, variant_ct, maj_alleles);
          }
        }
        BigstackReset(bigstack_mark_read_freqs);

        if (pcp->command_flags1 & kfCommand1AlleleFreq) {
          reterr = WriteAlleleFreqs(variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, nonref_flags, nonfounders? allele_ddosages : founder_allele_ddosages, imp_r2_vals, pcp->freq_ref_binstr, pcp->freq_alt1_binstr, raw_variant_ct, variant_ct, max_allele_ct, max_allele_slen, pgfi.gflags, pcp->freq_rpt_flags, pcp->max_thread_ct, nonfounders, outname, outname_end);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
          if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1PhenoSvd | kfCommand1AlleleFreq | kfCommand1WriteSamples | kfCommand1Validate | kfCommand1PgenInfo | kfCommand1RmDupList)))) {
            continue;
          }
        }
        if (pcp->command_flags1 & kfCommand1GenoCounts) {
          reterr = WriteGenoCounts(sample_include, sex_nm, sex_male, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, nonref_flags, raw_geno_cts, x_male_geno_cts, raw_sample_ct, sample_ct, male_ct, nosex_ct, raw_variant_ct, variant_ct, x_start, max_allele_slen, pgfi.gflags, pcp->geno_counts_flags, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
          if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1PhenoSvd | kfCommand1AlleleFreq | kfCommand1GenoCounts | kfCommand1WriteSamples | kfCommand1Validate | kfCommand1PgenInfo | kfCommand1RmDupList)))) {
            continue;
          }
        }

        if (pcp->command_flags1 & kfCommand1MissingReport) {
          reterr = WriteMissingnessReports(sample_include, &pii.sii, sex_nm, sex_male, pheno_cols, pheno_names, sample_missing_hc_cts, sample_missing_dosage_cts, sample_hethap_cts, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, nonref_flags, variant_missing_hc_cts, variant_missing_dosage_cts, variant_hethap_cts, raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_allele_slen, nosex_ct && (pcp->misc_flags & kfMiscYNosexMissingStats), pgfi.gflags, variant_hethap_cts? first_hap_uidx : 0x7fffffff, pcp->missing_rpt_flags, pcp->max_thread_ct, outname, outname_end);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
          if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1PhenoSvd | kfCommand1AlleleFreq | kfCommand1GenoCounts | kfCommand1MissingReport | kfCommand1WriteSamples | kfCommand1Validate | kfCommand1PgenInfo | kfCommand1RmDupList)))) {
            continue;
          }
        }

        if (pcp->geno_thresh != 1.0) {
          const uint32_t geno_hh_missing = S_CAST(uint32_t, pcp->misc_flags & kfMiscGenoHhMissing);
          EnforceGenoThresh(cip, (pcp->misc_flags & kfMiscGenoDosage)? variant_missing_dosage_cts : variant_missing_hc_cts, geno_hh_missing? variant_hethap_cts : nullptr, sample_ct, male_ct, geno_hh_missing? first_hap_uidx : 0x7fffffff, pcp->geno_thresh, variant_include, &variant_ct);
        }

        if ((pcp->command_flags1 & kfCommand1Hardy) || (pcp->hwe_ln_thresh != -DBL_MAX)) {
          if (cip->haploid_mask[0] & 1) {
            if (unlikely(pcp->command_flags1 & kfCommand1Hardy)) {
              logerrputs("Error: --hardy is pointless on an all-haploid genome.\n");
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            // could make hwe_ln_thresh non-const and set it to -DBL_MAX
            // earlier on in this case
            logerrputs("Warning: --hwe has no effect since entire genome is haploid.\n");
          } else {
            STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_geno_cts) = nonfounders? raw_geno_cts : founder_raw_geno_cts;
            STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_male_geno_cts) = nonfounders? x_male_geno_cts : founder_x_male_geno_cts;
            STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_nosex_geno_cts) = nonfounders? x_nosex_geno_cts : founder_x_nosex_geno_cts;

            // For multiallelic variants, perform one 'biallelic' test per
            // allele.  This is much, much faster than e.g. triallelic exact
            // tests, and the resulting set of p-values are arguably more
            // useful anyway.

            // One entry per autosomal multiallelic alt allele.
            // Only two cells per entry, since third can be inferred by
            // subtracting from the nonmissing genotype count (this is
            // potentially a memory hog).
            STD_ARRAY_PTR_DECL(uint32_t, 2, autosomal_xgeno_cts) = nullptr;
            // one entry per chrX multiallelic alt allele.
            STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts) = nullptr;
            STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts) = nullptr;
            uint32_t hwe_x_ct = 0;
            if (hwe_x_probs_needed) {
              hwe_x_ct = CountChrVariantsUnsafe(variant_include, cip, cip->xymt_codes[kChrOffsetX]);
              // hwe_x_ct == 0 possible when hwe_x_probs_needed set, if --geno
              // filters out all chrX variants
            }
            uintptr_t autosomal_xallele_ct = 0;
            uintptr_t x_xallele_ct = 0;
            if (allele_idx_offsets) {
              // bugfix (3 Mar 2023): if !hwe_x_probs_needed, need to subtract
              // chrX variant count here
              const uint32_t autosomal_variant_ct = variant_ct - hwe_x_ct - CountNonAutosomalVariants(variant_include, cip, 1 - hwe_x_probs_needed, 1);
              const uint32_t chr_ct = cip->chr_ct;
              for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
                const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
                if ((chr_idx != x_code) && (!IsSet(cip->haploid_mask, chr_idx))) {
                  autosomal_xallele_ct += CountExtraAlleles(variant_include, allele_idx_offsets, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], 1);
                }
              }
              if (hwe_x_ct) {
                x_xallele_ct = CountExtraAlleles(variant_include, allele_idx_offsets, x_start, x_start + x_len, 1);
              }
              if (autosomal_xallele_ct || x_xallele_ct) {
                if (autosomal_xallele_ct) {
                  if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 2, autosomal_xallele_ct, &autosomal_xgeno_cts))) {
                    goto Plink2Core_ret_NOMEM;
                  }
                }
                if (x_xallele_ct) {
                  if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 2, x_xallele_ct, &x_knownsex_xgeno_cts))) {
                    goto Plink2Core_ret_NOMEM;
                  }
                  if (hwe_x_male_geno_cts) {
                    if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 2, x_xallele_ct, &x_male_xgeno_cts))) {
                      goto Plink2Core_ret_NOMEM;
                    }
                  }
                }
                reterr = GetMultiallelicMarginalCounts(nonfounders? sample_include : founder_info, sex_nm, sex_male, variant_include, cip, allele_idx_offsets, hwe_geno_cts, raw_sample_ct, autosomal_variant_ct, autosomal_xallele_ct, hwe_x_ct, x_xallele_ct, &simple_pgr, x_male_xgeno_cts, autosomal_xgeno_cts, x_knownsex_xgeno_cts);
                if (unlikely(reterr)) {
                  goto Plink2Core_ret_1;
                }
              }
            }

            double* hwe_x_ln_pvals = nullptr;
            if (hwe_x_ct && ((pcp->hwe_ln_thresh != -DBL_MAX) || (pcp->hardy_flags & kfHardyColP))) {
              // support suppression of --hardy p column (with a gigantic
              // dataset, maybe it's reasonable to stick to femalep, etc.)
              uint32_t hwe_midp;
              if (pcp->command_flags1 & kfCommand1Hardy) {
                hwe_midp = (pcp->hardy_flags / kfHardyMidp) & 1;
                if (pcp->hwe_ln_thresh != -DBL_MAX) {
                  const uint32_t hwe_midp2 = (pcp->misc_flags / kfMiscHweMidp) & 1;
                  if (unlikely(hwe_midp != hwe_midp2)) {
                    // could support this efficiently, but why bother...
                    logerrputs("Error: --hardy and --hwe must have identical midp settings when chrX is\npresent.\n");
                    goto Plink2Core_ret_INVALID_CMDLINE;
                  }
                }
              } else {
                hwe_midp = (pcp->misc_flags / kfMiscHweMidp) & 1;
              }
              reterr = ComputeHweXLnPvals(variant_include, allele_idx_offsets, hwe_geno_cts, hwe_x_male_geno_cts, hwe_x_nosex_geno_cts, x_knownsex_xgeno_cts, x_male_xgeno_cts, x_start, hwe_x_ct, x_xallele_ct, hwe_midp, pcp->max_thread_ct, &hwe_x_ln_pvals);
              if (unlikely(reterr)) {
                goto Plink2Core_ret_1;
              }
            }
            if (pcp->command_flags1 & kfCommand1Hardy) {
              reterr = HardyReport(variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, nonref_flags, hwe_geno_cts, autosomal_xgeno_cts, hwe_x_male_geno_cts, hwe_x_nosex_geno_cts, x_knownsex_xgeno_cts, x_male_xgeno_cts, hwe_x_ln_pvals, variant_ct, hwe_x_ct, max_allele_slen, pgfi.gflags, pcp->output_min_ln, pcp->hardy_flags, pcp->max_thread_ct, nonfounders, outname, outname_end);
              if (unlikely(reterr)) {
                goto Plink2Core_ret_1;
              }
              if (!(pcp->command_flags1 & (~(kfCommand1GenotypingRate | kfCommand1PhenoSvd | kfCommand1AlleleFreq | kfCommand1GenoCounts | kfCommand1MissingReport | kfCommand1Hardy | kfCommand1WriteSamples | kfCommand1Validate | kfCommand1PgenInfo | kfCommand1RmDupList)))) {
                continue;
              }
            }
            if (pcp->hwe_ln_thresh != -DBL_MAX) {
              // assumes no filtering between hwe_x_ln_pvals[] computation and
              // here
              reterr = EnforceHweThresh(cip, allele_idx_offsets, hwe_geno_cts, autosomal_xgeno_cts, hwe_x_male_geno_cts, hwe_x_nosex_geno_cts, x_knownsex_xgeno_cts, x_male_xgeno_cts, hwe_x_ln_pvals, pcp->misc_flags, pcp->hwe_ln_thresh, pcp->hwe_sample_size_term, nonfounders, variant_include, &variant_ct);
              if (unlikely(reterr)) {
                goto Plink2Core_ret_1;
              }
            }
          }
        }
        // raw_geno_cts/founder_raw_geno_cts/hwe_x_ln_pvals no longer needed
        BigstackReset(bigstack_mark_geno_cts);

        if ((pcp->min_maf != 0.0) || (pcp->max_maf != 1.0) || pcp->min_allele_ddosage || (pcp->max_allele_ddosage != (~0LLU))) {
          EnforceFreqConstraints(allele_idx_offsets, nonfounders? allele_ddosages : founder_allele_ddosages, allele_freqs, pcp->filter_modes, pcp->min_maf, pcp->max_maf, pcp->min_allele_ddosage, pcp->max_allele_ddosage, variant_include, &variant_ct);
        }

        if (imp_r2_vals) {
          if (pcp->mach_r2_max != 0.0) {
            EnforceImpR2Thresh(cip, imp_r2_vals, pcp->mach_r2_min, pcp->mach_r2_max, 0, variant_include, &variant_ct);
          } else if (pcp->minimac3_r2_max != 0.0) {
            EnforceImpR2Thresh(cip, imp_r2_vals, pcp->minimac3_r2_min, pcp->minimac3_r2_max, 1, variant_include, &variant_ct);
          }
          BigstackReset(imp_r2_vals);
        }
      }

      if (pcp->min_bp_space) {
        if (vpos_sortstatus & kfUnsortedVarBp) {
          logerrputs("Error: --bp-space requires a sorted .pvar/.bim.  Retry this command after using\n--make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        EnforceMinBpSpace(cip, variant_bps, pcp->min_bp_space, variant_include, &variant_ct);
      }

      if (pcp->filter_flags & kfFilterPvarReq) {
        if (unlikely(!variant_ct)) {
          // do we want this to be conditionally acceptable?
          logerrputs("Error: No variants remaining after main filters.\n");
          goto Plink2Core_ret_DEGENERATE_DATA;
        }
        logprintf("%u variant%s remaining after main filters.\n", variant_ct, (variant_ct == 1)? "" : "s");
      }

      if ((pcp->command_flags1 & kfCommand1MendelReport) || (pcp->filter_flags & kfFilterMendel)) {
        if (vpos_sortstatus & kfUnsortedVarBp) {
          logerrputs("Error: --me and --mendel require a sorted .pvar/.bim.  Retry this command after\nusing --make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        const uint32_t generate_reports = (pcp->command_flags1 / kfCommand1MendelReport) & 1;
        reterr = MendelErrorScan(&pii, founder_info, sex_nm, sex_male, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, &(pcp->mendel_info), raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, max_allele_ct, max_allele_slen, generate_reports, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, sample_include, variant_include, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        if (pcp->filter_flags & kfFilterMendel) {
          sample_ct = PopcountWords(sample_include, raw_sample_ctl);
          variant_ct = PopcountWords(variant_include, raw_variant_ctl);
          logprintf("--me: %u sample%s and %u variant%s remaining.\n", sample_ct, (sample_ct == 1)? "" : "s", variant_ct, (variant_ct == 1)? "" : "s");
          UpdateSampleSubsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
        }
      }

      if (pcp->command_flags1 & kfCommand1SampleCounts) {
        reterr = SampleCounts(sample_include, &pii.sii, sex_nm, sex_male, variant_include, cip, allele_idx_offsets, allele_storage, raw_sample_ct, sample_ct, male_ct, raw_variant_ct, variant_ct, max_allele_ct, pcp->sample_counts_flags, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & kfCommand1Sdiff) {
        reterr = Sdiff(sample_include, &pii.sii, sex_nm, sex_male, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, &(pcp->sdiff_info), raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, (pcp->misc_flags / kfMiscIidSid) & 1, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & (kfCommand1MakeKing | kfCommand1KingCutoff)) {
        uintptr_t* prev_sample_include = nullptr;
        const uint32_t prev_sample_ct = sample_ct;
        if (pcp->command_flags1 & kfCommand1KingCutoff) {
          if (unlikely(bigstack_alloc_w(raw_sample_ctl, &prev_sample_include))) {
            goto Plink2Core_ret_NOMEM;
          }
          memcpy(prev_sample_include, sample_include, raw_sample_ctl * sizeof(intptr_t));
        }
        RelConcordanceCheckMode rel_or_concordance_check = kRcCheck0;
        if (pcp->king_flags & kfKingRelCheck) {
          if (OnlyOneFid(sample_include, &pii.sii, sample_ct)) {
            logerrputs("Warning: --make-king-table 'rel-check' modifier has no effect since only one\nFID is present.\n");
          } else {
            rel_or_concordance_check = kRcCheckRel;
          }
        } else if (pcp->king_flags & kfKingConcordanceCheck) {
          rel_or_concordance_check = kRcCheckConcordance;
        }
        if (pcp->king_table_require_fnames || pcp->king_table_subset_fname || rel_or_concordance_check) {
          // command-line parser currently guarantees
          // --king-table-subset/--king-table-require[-xor] and
          // "--make-king-table rel-check" aren't used with
          // --king-cutoff[-table] or --make-king
          reterr = CalcKingTableSubset(sample_include, &pii.sii, variant_include, cip, pcp->king_table_subset_fname, pcp->king_table_require_fnames, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, pcp->king_table_filter, pcp->king_table_subset_thresh, rel_or_concordance_check, pcp->king_flags, pcp->parallel_idx, pcp->parallel_tot, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        } else {
          if (king_cutoff_fprefix) {
            if (pcp->king_flags & kfKingCutoffTable) {
              reterr = KingCutoffBatchTable(&pii.sii, king_cutoff_fprefix, raw_sample_ct, pcp->king_cutoff, sample_include, &sample_ct);
            } else {
              reterr = KingCutoffBatchBinary(&pii.sii, raw_sample_ct, pcp->king_cutoff, sample_include, king_cutoff_fprefix, &sample_ct);
            }
          } else {
            reterr = CalcKing(&pii.sii, variant_include, cip, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, pcp->king_cutoff, pcp->king_table_filter, pcp->king_flags, pcp->parallel_idx, pcp->parallel_tot, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, &simple_pgr, sample_include, &sample_ct, outname, outname_end);
          }
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
          if (pcp->command_flags1 & kfCommand1KingCutoff) {
            snprintf(outname_end, kMaxOutfnameExtBlen, ".king.cutoff.in.id");
            reterr = WriteSampleIds(sample_include, &pii.sii, outname, sample_ct);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
            snprintf(&(outname_end[13]), kMaxOutfnameExtBlen - 13, "out.id");
            BitvecInvmask(sample_include, raw_sample_ctl, prev_sample_include);
            const uint32_t removed_sample_ct = prev_sample_ct - sample_ct;
            reterr = WriteSampleIds(prev_sample_include, &pii.sii, outname, removed_sample_ct);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
            BigstackReset(prev_sample_include);
            outname_end[13] = '\0';
            logprintfww("--king-cutoff%s: Excluded sample ID%s written to %sout.id , and %u remaining sample ID%s written to %sin.id .\n", (pcp->king_flags & kfKingCutoffTable)? "-table" : "", (removed_sample_ct == 1)? "" : "s", outname, sample_ct, (sample_ct == 1)? "" : "s", outname);
            UpdateSampleSubsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
          }
        }
      }
      if ((pcp->command_flags1 & kfCommand1MakeRel) || keep_grm) {
        reterr = CalcGrm(sample_include, &pii.sii, variant_include, cip, allele_idx_offsets, allele_freqs, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, max_allele_ct, pcp->grm_flags, pcp->grm_sparse_cutoff, pcp->parallel_idx, pcp->parallel_tot, pcp->max_thread_ct, &simple_pgr, outname, outname_end, keep_grm? (&grm) : nullptr);
        if (unlikely(reterr)) {
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
        //         the R package.)  This minimizes the risk of getting PCs
        //         which map to e.g. a sample-duplicate cluster instead of
        //         broad ancestry.
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
        reterr = CalcPca(sample_include, &pii.sii, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, maj_alleles, allele_freqs, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, max_allele_ct, max_allele_slen, pcp->pca_ct, pcp->pca_flags, pcp->max_thread_ct, &simple_pgr, sfmtp, grm, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
#endif

      if (pcp->command_flags1 & kfCommand1CheckOrImputeSex) {
        reterr = CheckOrImputeSex(sample_include, &pii.sii, variant_include, cip, allele_idx_offsets, allele_freqs, &(pcp->check_sex_info), raw_sample_ct, sample_ct, raw_variant_ct, max_allele_ct, pcp->max_thread_ct, pgr_alloc_cacheline_ct, sex_nm, sex_male, &pgfi, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        if (pcp->check_sex_info.flags & kfCheckSexImpute) {
          nosex_ct = sample_ct - PopcountWords(sex_nm, raw_sample_ctl);
          male_ct = PopcountWords(sex_male, raw_sample_ctl);
        }
      }

      if (pcp->command_flags1 & kfCommand1WriteSnplist) {
        reterr = WriteSnplist(variant_include, variant_ids, variant_ct, (pcp->misc_flags / kfMiscWriteSnplistZs) & 1, (pcp->misc_flags / kfMiscWriteSnplistAllowDups) & 1, pcp->max_thread_ct, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & (kfCommand1MakePlink2 | kfCommand1Exportf | kfCommand1WriteCovar)) {
        // If non-null, allele_permute[allele_idx_offsets[j]] stores the old
        // allele-index of the new ref allele,
        // allele_permute[allele_idx_offsets[j] + 1] stores the old
        // allele-index of the new alt1 allele, etc.
        // Operations which instantiate this (--maj-ref, --ref-allele,
        // --alt[1]-allele) are only usable with fileset creation commands.  No
        // more pass-marker_reverse-to-everything nonsense.
        // (Technically, I could also drop support for --export, but that would
        // force too many real-world jobs to require two plink2 runs instead of
        // one.)

        const uint32_t setting_alleles_from_file = pcp->ref_allele_flag || pcp->alt_allele_flag || pcp->fa_fname;
        uint32_t* variant_bps_backup = nullptr;
        const char** allele_storage_backup = nullptr;
        uint32_t max_allele_slen_backup = max_allele_slen;
        uintptr_t* nonref_flags_backup = nullptr;
        uint32_t nonref_flags_was_null = (nonref_flags == nullptr);
        AlleleCode* allele_permute = nullptr;
        uint32_t* contig_lens = nullptr;
        if (pcp->fa_flags & kfFaScrapeLengths) {
          if (unlikely(bigstack_end_calloc_u32(cip->max_code + 1 + cip->name_ct, &contig_lens))) {
            goto Plink2Core_ret_NOMEM;
          }
        }
        if ((pcp->misc_flags & kfMiscMajRef) || setting_alleles_from_file || contig_lens) {
          const uint32_t need_allele_permute = ((pcp->misc_flags & kfMiscMajRef) || pcp->ref_allele_flag || pcp->alt_allele_flag || (pcp->fa_flags & kfFaRefFrom));
          if (loop_cats_idx + 1 < loop_cats_ct) {
            // --ref-allele/--alt[1]-allele/--normalize may alter
            // allele_storage[] and (in the former two cases) max_allele_slen;
            // --loop-cats doesn't like this.  Save a backup copy.
            if (pcp->ref_allele_flag || pcp->alt_allele_flag || (pcp->fa_flags & kfFaNormalize)) {
              if (pcp->fa_flags & kfFaNormalize) {
                if (unlikely(bigstack_end_alloc_u32(raw_variant_ct, &variant_bps_backup))) {
                  goto Plink2Core_ret_NOMEM;
                }
                memcpy(variant_bps_backup, variant_bps, raw_variant_ct * sizeof(int32_t));
              }
              if (unlikely(bigstack_end_alloc_kcp(raw_allele_ct, &allele_storage_backup))) {
                goto Plink2Core_ret_NOMEM;
              }
              memcpy(allele_storage_backup, allele_storage, raw_allele_ct * sizeof(intptr_t));
            }
            // nonref_flags may be altered by everything but --normalize.
            if (nonref_flags && need_allele_permute) {
              if (unlikely(bigstack_end_alloc_w(raw_variant_ctl, &nonref_flags_backup))) {
                goto Plink2Core_ret_NOMEM;
              }
              memcpy(nonref_flags_backup, nonref_flags, raw_variant_ctl * sizeof(intptr_t));
            }
          }
          const uint32_t not_all_nonref = !(pgfi.gflags & kfPgenGlobalAllNonref);
          if (need_allele_permute) {
            // no need to track bigstack_end_mark before this is allocated,
            // etc., due to the restriction to --make-pgen/--export
            if (unlikely(bigstack_end_alloc_ac(raw_allele_ct, &allele_permute))) {
              goto Plink2Core_ret_NOMEM;
            }
            reterr = InitAllelePermuteUnsafe(allele_idx_offsets, raw_variant_ct, pcp->max_thread_ct, allele_permute);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
            if ((not_all_nonref || setting_alleles_from_file) && (!nonref_flags)) {
              if (unlikely(bigstack_end_alloc_w(raw_variant_ctl, &nonref_flags))) {
                goto Plink2Core_ret_NOMEM;
              }
              pgfi.nonref_flags = nonref_flags;
              // make it clear that this is probably worth cleaning up a bit...
              GET_PRIVATE(simple_pgr, m).fi.nonref_flags = nonref_flags;
              if (not_all_nonref) {
                ZeroWArr(raw_variant_ctl, nonref_flags);
              } else {
                SetAllBits(raw_variant_ct, nonref_flags);
              }
            }
          }
          uintptr_t* previously_seen = nullptr;
          if (pcp->ref_allele_flag) {
            if (pcp->alt_allele_flag) {
              if (unlikely(bigstack_alloc_w(raw_variant_ctl, &previously_seen))) {
                goto Plink2Core_ret_NOMEM;
              }
            }
            reterr = SetRefalt1FromFile(variant_include, variant_ids, allele_idx_offsets, pcp->ref_allele_flag, raw_variant_ct, variant_ct, max_variant_id_slen, max_allele_ct, kRefalt1ModeRef, (pcp->misc_flags / kfMiscRefAlleleForce) & 1, pcp->input_missing_geno_char, pcp->max_thread_ct, allele_storage, &max_allele_slen, allele_permute, nonref_flags, previously_seen);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
          }
          if (pcp->alt_allele_flag) {
            reterr = SetRefalt1FromFile(variant_include, variant_ids, allele_idx_offsets, pcp->alt_allele_flag, raw_variant_ct, variant_ct, max_variant_id_slen, max_allele_ct, (pcp->misc_flags & kfMiscAlt1Allele)? kRefalt1ModeAlt1 : kRefalt1ModeAlt, (pcp->misc_flags / kfMiscAltAlleleForce) & 1, pcp->input_missing_geno_char, pcp->max_thread_ct, allele_storage, &max_allele_slen, allele_permute, nonref_flags, previously_seen);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
            if (previously_seen) {
              BigstackReset(previously_seen);
            }
          }
          if (pcp->fa_flags & (kfFaRefFrom | kfFaNormalize | kfFaScrapeLengths)) {
            // for sanity's sake, --maj-ref, --ref-allele/--alt[1]-allele, and
            // --ref-from-fa are mutually exclusive
            // (though --ref-from-fa + --alt[1]-allele may be permitted later)
            if (unlikely((pcp->fa_flags & (kfFaRefFrom | kfFaNormalize)) && (vpos_sortstatus & kfUnsortedVarBp))) {
              logerrputs("Error: --normalize and --ref-from-fa require a sorted .pvar/.bim.  Retry this\ncommand after using --make-pgen/--make-bed + --sort-vars to sort your data.\n");
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            reterr = ProcessFa(variant_include, variant_ids, allele_idx_offsets, cip, pcp->fa_fname, max_allele_ct, max_allele_slen, pcp->fa_flags, ctou32(pcp->output_missing_geno_char), pcp->max_thread_ct, &vpos_sortstatus, variant_bps, allele_storage, allele_permute, nonref_flags, contig_lens, outname, outname_end);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
          }
          if (pcp->misc_flags & kfMiscMajRef) {
            // Since this also sorts ALT alleles by frequency, it can't just
            // subscribe to maj_alleles[].
            const uint64_t* main_allele_ddosages = nonfounders? allele_ddosages : founder_allele_ddosages;
            const uint32_t skip_real_ref = not_all_nonref && (!(pcp->misc_flags & kfMiscMajRefForce));
            if (skip_real_ref && AllWordsAreZero(nonref_flags, raw_variant_ctl)) {
              logerrputs("Warning: --maj-ref has no effect, since no provisional reference alleles are\npresent.  (Did you forget to add the 'force' modifier?)\n");
            } else {
              reterr = MajRef(variant_include, allele_idx_offsets, main_allele_ddosages, variant_ct, max_allele_ct, skip_real_ref, allele_permute, nonref_flags);
              if (unlikely(reterr)) {
                goto Plink2Core_ret_1;
              }
            }
          }
        }
        if (allele_presents) {
          if (!allele_permute) {
            if (unlikely(bigstack_end_alloc_ac(raw_allele_ct, &allele_permute))) {
              goto Plink2Core_ret_NOMEM;
            }
            reterr = InitAllelePermuteUnsafe(allele_idx_offsets, raw_variant_ct, pcp->max_thread_ct, allele_permute);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
          }
          reterr = TrimAllelePermute(variant_include, allele_idx_offsets, allele_presents, variant_ct, pcp->max_thread_ct, allele_permute);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
        BigstackReset(bigstack_mark_allele_ddosages);

        uint32_t* new_sample_idx_to_old = nullptr;
        if (pcp->sample_sort_mode > kSortNone) {
          if (sample_ct < 2) {
            logerrputs("Warning: Skipping --sample-sort since <2 samples are present.\n");
          } else {
            if (pcp->sample_sort_mode == kSortFile) {
              reterr = SampleSortFileMap(sample_include, &pii.sii, pcp->sample_sort_fname, raw_sample_ct, sample_ct, &new_sample_idx_to_old);
              if (unlikely(reterr)) {
                goto Plink2Core_ret_1;
              }
            } else {
              // probably more efficient to have --make-{bed,pgen,bpgen}
              // perform an unfiltered load?  but we should have compute power
              // to spare here, so keep the code simpler for now
              char* sorted_xidbox_tmp;
              uintptr_t max_xid_blen;
              reterr = SortedXidboxInitAlloc(sample_include, &pii.sii, sample_ct, pii.sii.sids? kfXidModeFidIidSid : kfXidModeFidIid, (pcp->sample_sort_mode == kSortNatural), &sorted_xidbox_tmp, &new_sample_idx_to_old, &max_xid_blen);
              if (unlikely(reterr)) {
                goto Plink2Core_ret_1;
              }
              BigstackReset(sorted_xidbox_tmp);
            }
            logprintf("--indiv-sort: %u samples reordered.\n", sample_ct);
          }
        }

        // update (18 Mar 2018): permit no-covariate --write-covar when used to
        // write phenotypes.
        if ((covar_ct || (pheno_ct && (pcp->write_covar_flags & (kfWriteCovarColPheno1 | kfWriteCovarColPhenos)))) && ((pcp->command_flags1 & (kfCommand1Exportf | kfCommand1WriteCovar)) || ((pcp->command_flags1 & kfCommand1MakePlink2) && (make_plink2_flags & (kfMakeBed | kfMakeFam | kfMakePgen | kfMakePsam))))) {
          reterr = WriteCovar(sample_include, &pii, sex_nm, sex_male, pheno_cols, pheno_names, covar_cols, covar_names, new_sample_idx_to_old, pcp->output_missing_pheno, sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, pcp->write_covar_flags, outname, outname_end);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        } else if (pcp->command_flags1 & kfCommand1WriteCovar) {
          if (pcp->write_covar_flags & (kfWriteCovarColPheno1 | kfWriteCovarColPhenos)) {
            logerrputs("Warning: Skipping --write-covar, since no phenotypes or covariates are loaded.\n");
          } else {
            logerrputs("Warning: Skipping --write-covar, since no covariates are loaded.\n");
          }
        }

        if (pcp->command_flags1 & kfCommand1MakePlink2) {
          if (pcp->sort_vars_mode > kSortNone) {
            reterr = MakePlink2Vsort(sample_include, &pii, founder_info, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, variant_include, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_presents, allele_permute, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, info_reload_slen? pvarname : nullptr, variant_cms, pcp->output_missing_pheno, pcp->legacy_output_missing_pheno, contig_lens, pcp->rename_chrs_fname, pcp->update_chr_flag, (make_plink2_flags & kfMakePgenWriterVer)? ver_str : nullptr, pcp->zero_cluster_fname, pcp->zero_cluster_phenoname, &(pcp->flip_info), xheader_blen, info_flags, raw_sample_ct, sample_ct, male_ct, nosex_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_allele_ct, max_variant_id_slen, max_allele_slen, max_filter_slen, info_reload_slen, pcp->output_missing_geno_char, pcp->max_thread_ct, pcp->hard_call_thresh, pcp->dosage_erase_thresh, pcp->misc_flags, make_plink2_flags, (pcp->sort_vars_mode == kSortNatural), pcp->pvar_psam_flags, (pcp->mendel_info.flags / kfMendelDuos) & 1, cip, xheader, chr_idxs, &simple_pgr, outname, outname_end);
          } else {
            if (vpos_sortstatus & kfUnsortedVarBp) {
              logerrputs("Warning: Variants are not sorted by position.  Consider rerunning with the\n--sort-vars flag added to remedy this.\n");
            }
            reterr = MakePlink2NoVsort(sample_include, &pii, founder_info, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, info_reload_slen? pvarname : nullptr, variant_cms, pcp->varid_template_str, pcp->varid_multi_template_str, pcp->varid_multi_nonsnp_template_str, pcp->missing_varid_match, pcp->output_missing_pheno, pcp->legacy_output_missing_pheno, contig_lens, (make_plink2_flags & kfMakePgenWriterVer)? ver_str : nullptr, pcp->zero_cluster_fname, pcp->zero_cluster_phenoname, &(pcp->flip_info), xheader_blen, info_flags, raw_sample_ct, sample_ct, male_ct, nosex_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_allele_ct, max_variant_id_slen, max_allele_slen, max_filter_slen, info_reload_slen, pcp->output_missing_geno_char, pcp->max_thread_ct, pcp->hard_call_thresh, pcp->dosage_erase_thresh, pcp->new_variant_id_max_allele_slen, pcp->misc_flags, make_plink2_flags, pcp->pvar_psam_flags, (pcp->mendel_info.flags / kfMendelDuos) & 1, pgr_alloc_cacheline_ct, xheader, &pgfi, &simple_pgr, outname, outname_end);
          }
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
          // no BigstackReset needed here, since allele_presents only needed
          // if 'trim-alts', and later operations are prohibited in that case
        }

        if (pcp->command_flags1 & kfCommand1Exportf) {
          reterr = Exportf(sample_include, &pii, sex_nm, sex_male, pheno_cols, pheno_names, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, info_reload_slen? pvarname : nullptr, variant_cms, &(pcp->exportf_info), pcp->legacy_output_missing_pheno, contig_lens, xheader_blen, info_flags, raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_variant_id_slen, max_allele_slen, max_filter_slen, info_reload_slen, pcp->input_missing_geno_char, pcp->output_missing_geno_char, pcp->legacy_output_missing_geno_char, pcp->max_thread_ct, make_plink2_flags, pgr_alloc_cacheline_ct, xheader, &pgfi, &simple_pgr, outname, outname_end);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }

        if (variant_bps_backup) {
          memcpy(variant_bps, variant_bps_backup, raw_variant_ct * sizeof(int32_t));
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
          GET_PRIVATE(simple_pgr, m).fi.nonref_flags = nullptr;
        }
      }
      BigstackReset(bigstack_mark_allele_ddosages);

      if (pcp->command_flags1 & kfCommand1PgenDiff) {
        if (unlikely(vpos_sortstatus & kfUnsortedVarBp)) {
          logerrputs("Error: --pgen-diff requires sorted .pvar/.bim files.  Retry this command after\nusing --make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        reterr = PgenDiff(sample_include, &pii.sii, sex_nm, sex_male, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, &(pcp->pgen_diff_info), raw_sample_ct, sample_ct, raw_variant_ct, max_allele_ct, max_allele_slen, pcp->input_missing_geno_char, pcp->max_thread_ct, &pgfi, &simple_pgr, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & kfCommand1LdPrune) {
        if (unlikely((pcp->ld_info.prune_flags & kfLdPruneWindowBp) && (vpos_sortstatus & kfUnsortedVarBp))) {
          logerrputs("Error: When the window size is in kb units, LD-based pruning requires a sorted\n.pvar/.bim.  Retry this command after using --make-pgen/--make-bed +\n--sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        reterr = LdPrune(variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, maj_alleles, allele_freqs, founder_info, sex_nm, sex_male, &(pcp->ld_info), pcp->indep_preferred_fname, raw_variant_ct, variant_ct, raw_sample_ct, founder_ct, nosex_ct, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & kfCommand1Ld) {
        reterr = LdConsole(variant_include, cip, variant_ids, allele_idx_offsets, allele_storage, maj_alleles, founder_info, sex_nm, sex_male, &(pcp->ld_info), variant_ct, raw_sample_ct, founder_ct, &simple_pgr);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & kfCommand1Vcor) {
        if (unlikely(vpos_sortstatus & kfUnsortedVarBp)) {
          logerrputs("Error: --r[2]-[un]phased runs require a sorted .pvar/.bim.  Retry this command\nafter using --make-pgen/--make-bed + --sort-vars to sort your data.\n");
          return kPglRetInconsistentInput;
        }
        if (unlikely((vpos_sortstatus & kfUnsortedVarCm) && (pcp->vcor_info.cm_radius != -1.0))) {
          logerrputs("Error: --ld-window-cm requires nondecreasing CM values on each chromosome.\nRetry this command after regenerating your CM coordinates.\n");
          return kPglRetInconsistentInput;
        }
        reterr = Vcor(variant_include, cip, variant_bps, variant_ids, variant_cms, allele_idx_offsets, allele_storage, maj_alleles, allele_freqs, founder_info, sex_nm, sex_male, &(pcp->vcor_info), raw_variant_ct, variant_ct, raw_sample_ct, founder_ct, max_variant_id_slen, max_allele_slen, pcp->parallel_idx, pcp->parallel_tot, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & kfCommand1Het) {
        reterr = HetReport(sample_include, &pii.sii, variant_include, cip, allele_idx_offsets, allele_freqs, founder_info, raw_sample_ct, sample_ct, founder_ct, raw_variant_ct, variant_ct, max_allele_ct, pcp->het_flags, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & kfCommand1Fst) {
        reterr = FstReport(sample_include, sex_male, pheno_cols, pheno_names, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, &(pcp->fst_info), raw_sample_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_allele_ct, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (pcp->command_flags1 & kfCommand1Score) {
        reterr = ScoreReport(sample_include, &pii.sii, sex_nm, sex_male, pheno_cols, pheno_names, variant_include, cip, variant_ids, allele_idx_offsets, allele_storage, allele_freqs, &(pcp->score_info), pcp->output_missing_pheno, raw_sample_ct, sample_ct, nosex_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_variant_id_slen, pcp->xchr_model, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->command_flags1 & kfCommand1Vscore) {
        reterr = Vscore(variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, sample_include, &pii.sii, sex_male, allele_freqs, pcp->vscore_fname, &(pcp->vscore_col_idx_range_list), raw_variant_ct, variant_ct, raw_sample_ct, sample_ct, nosex_ct, max_allele_slen, pcp->vscore_flags, pcp->xchr_model, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      // eventually check for nonzero pheno_ct here?

      if (pcp->command_flags1 & kfCommand1Glm) {
        if (unlikely(pcp->glm_info.local_first_covar_col && (vpos_sortstatus & kfUnsortedVarBp))) {
          logerrputs("Error: --glm + local-pos-cols= requires a sorted .pvar/.bim.  Retry this\ncommand after using --make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        reterr = GlmMain(sample_include, &pii.sii, sex_nm, sex_male, pheno_cols, pheno_names, covar_cols, covar_names, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, maj_alleles, allele_storage, &(pcp->glm_info), &(pcp->adjust_info), &(pcp->aperm), pcp->glm_local_covar_fname, pcp->glm_local_pvar_fname, pcp->glm_local_psam_fname, &(pcp->gwas_ssf_info), raw_sample_ct, sample_ct, pheno_ct, max_pheno_name_blen, covar_ct, max_covar_name_blen, raw_variant_ct, variant_ct, max_variant_id_slen, max_allele_slen, pcp->xchr_model, pcp->ci_size, pcp->vif_thresh, pcp->ln_pfilter, pcp->output_min_ln, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, &simple_pgr, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
    }
    if (pcp->command_flags1 & kfCommand1Clump) {
      if (unlikely(vpos_sortstatus & kfUnsortedVarBp)) {
        logerrputs("Error: --clump requires a sorted .pvar/.bim.  Retry this command after using\n--make-pgen/--make-bed + --sort-vars to sort your data.\n");
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      reterr = ClumpReports(variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, founder_info, sex_nm, sex_male, &(pcp->clump_info), raw_variant_ct, variant_ct, raw_sample_ct, founder_ct, nosex_ct, max_variant_id_slen, max_allele_slen, pcp->output_min_ln, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, &simple_pgr, outname, outname_end);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
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
  Plink2Core_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  Plink2Core_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  Plink2Core_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  Plink2Core_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 Plink2Core_ret_1:
  if (loop_cats_pheno_col) {
    // Current implementation requires this to happen before CleanupPhenoCols()
    // on pheno_cols/covar_cols, since loop_cats_pheno_col actually points to
    // the last element of one of those arrays.
    vecaligned_free_cond(loop_cats_pheno_col->nonmiss);
  }
  CleanupPhenoCols(covar_ct, covar_cols);
  CleanupPhenoCols(pheno_ct, pheno_cols);
  free_cond(covar_names);
  free_cond(pheno_names);
  CleanupPgr2(".pgen file", &simple_pgr, &reterr);
  free_cond(ext_slot.contents);
  CleanupPgfi2(".pgen file", &pgfi, &reterr);
  assert(pgfi.block_base == nullptr);
  // no BigstackReset() needed?
  return reterr;
}

PglErr ZstDecompress(const char* in_fname, const char* out_fname) {
  zstRFILE zrf;
  PreinitZstRfile(&zrf);
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    reterr = ZstRfileOpen(in_fname, &zrf);
    if (unlikely(reterr)) {
      if (reterr == kPglRetNomem) {
        goto ZstDecompress_ret_NOMEM;
      }
      if (reterr == kPglRetOpenFail) {
        fprintf(stderr, kErrprintfFopen, in_fname, strerror(errno));
        goto ZstDecompress_ret_OPEN_FAIL;
      }
      goto ZstDecompress_ret_RFILE_FAIL;
    }
    if (out_fname) {
      outfile = fopen(out_fname, FOPEN_WB);
      if (unlikely(!outfile)) {
        fprintf(stderr, kErrprintfFopen, out_fname, strerror(errno));
        goto ZstDecompress_ret_OPEN_FAIL;
      }
    } else {
      outfile = stdout;
    }
    unsigned char* buf = R_CAST(unsigned char*, g_textbuf);
    while (1) {
      const int32_t bytes_read = zstread(&zrf, buf, kTextbufMainSize);
      if (bytes_read <= 0) {
        if (likely(!bytes_read)) {
          break;
        }
        reterr = ZstRfileErrcode(&zrf);
        goto ZstDecompress_ret_RFILE_FAIL;
      }
      if (unlikely(!fwrite_unlocked(buf, bytes_read, 1, outfile))) {
        goto ZstDecompress_ret_WRITE_FAIL;
      }
      fflush(outfile);
    }
    if (out_fname) {
      if (unlikely(fclose_null(&outfile))) {
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
  ZstDecompress_ret_RFILE_FAIL:
    if (reterr == kPglRetReadFail) {
      fprintf(stderr, kErrprintfFread, in_fname, zsterror(&zrf));
    } else {
      fprintf(stderr, kErrprintfDecompress, in_fname, zsterror(&zrf));
    }
    break;
  ZstDecompress_ret_WRITE_FAIL:
    fprintf(stderr, kErrstrWrite, strerror(errno));
    reterr = kPglRetWriteFail;
    break;
  }
  if (out_fname) {
    fclose_cond(outfile);
  }
  if (unlikely(CleanupZstRfile(&zrf, &reterr))) {
    fprintf(stderr, kErrprintfFread, in_fname, rstrerror(errno));
  }
  return reterr;
}

PglErr Alloc2col(const char* const* sources, const char* flagname_p, uint32_t param_ct, TwoColParams** tcbuf) {
  uint32_t fname_blen = strlen(sources[0]) + 1;
  if (unlikely(fname_blen > kPglFnamesize)) {
    logerrprintf("Error: --%s filename too long.\n", flagname_p);
    return kPglRetOpenFail;
  }
  if (unlikely(!FileExists(sources[0]))) {
    logerrprintfww("Error: --%s: %s does not exist.\n", flagname_p, sources[0]);
    return kPglRetOpenFail;
  }
  if (unlikely(pgl_malloc(offsetof(TwoColParams, fname) + fname_blen, tcbuf))) {
    return kPglRetNomem;
  }
  memcpy((*tcbuf)->fname, sources[0], fname_blen);
  (*tcbuf)->skip_ct = 0;
  (*tcbuf)->skipchar = '\0';
  if (param_ct > 1) {
    if (unlikely(ScanPosintDefcapx(sources[1], &((*tcbuf)->colx)))) {
      logerrprintf("Error: Invalid --%s column number.\n", flagname_p);
      return kPglRetInvalidCmdline;
    }
    if (param_ct > 2) {
      if (unlikely(ScanPosintDefcapx(sources[2], &((*tcbuf)->colid)))) {
        logerrprintf("Error: Invalid --%s variant ID column number.\n", flagname_p);
        return kPglRetInvalidCmdline;
      }
      if (param_ct == 4) {
        char cc = sources[3][0];
        if ((cc < '0') || (cc > '9')) {
          cc = ExtractCharParam(sources[3]);
          if (unlikely(!cc)) {
            goto Alloc2col_invalid_skip;
          }
          (*tcbuf)->skipchar = cc;
        } else {
          if (unlikely(ScanUintDefcapx(sources[3], &((*tcbuf)->skip_ct)))) {
          Alloc2col_invalid_skip:
            logerrprintfww("Error: Invalid --%s skip argument.  This needs to either be a single character (usually '#') which, when present at the start of a line, indicates it should be skipped; or the number of initial lines to skip. (Note that in shells such as bash, '#' is a special character that must be surrounded by single- or double-quotes to be parsed correctly.)\n", flagname_p);
            return kPglRetInvalidCmdline;
          }
        }
      }
    } else {
      (*tcbuf)->colid = 1;
    }
    if (unlikely((*tcbuf)->colx == (*tcbuf)->colid)) {
      logerrprintf("Error: Column numbers for --%s cannot be equal.\n%s", flagname_p, errstr_append);
      return kPglRetInvalidCmdline;
    }
  } else {
    (*tcbuf)->colx = 2;
    (*tcbuf)->colid = 1;
  }
  return kPglRetSuccess;
}

// flagname_p only needed when check_file_existence true
PglErr AllocAndFlattenCommaDelimEx(const char* const* sources, const char* flagname_p, uint32_t param_ct, uint32_t check_file_existence, char** flattened_buf_ptr) {
  uint32_t tot_blen = 1;
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
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
  if (unlikely(pgl_malloc(tot_blen, &write_iter))) {
    return kPglRetNomem;
  }
  *flattened_buf_ptr = write_iter;
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
    const char* cur_param_iter = sources[param_idx];
    while (1) {
      while (*cur_param_iter == ',') {
        ++cur_param_iter;
      }
      const char* cur_token_end = strchr(cur_param_iter, ',');
      if (!cur_token_end) {
        break;
      }
      char* write_iter_next = memcpyax(write_iter, cur_param_iter, cur_token_end - cur_param_iter, '\0');
      if (check_file_existence) {
        if (unlikely(!FileExists(write_iter))) {
          logerrprintfww("Error: --%s: %s does not exist.\n", flagname_p, write_iter);
          return kPglRetOpenFail;
        }
      }
      write_iter = write_iter_next;
      cur_param_iter = &(cur_token_end[1]);
    }
    write_iter = strcpyax(write_iter, cur_param_iter, '\0');
  }
  *write_iter = '\0';
  return kPglRetSuccess;
}

static inline PglErr AllocAndFlattenCommaDelim(const char* const* sources, uint32_t param_ct, char** flattened_buf_ptr) {
  return AllocAndFlattenCommaDelimEx(sources, nullptr, param_ct, 0, flattened_buf_ptr);
}

void PrintVer() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}

uint32_t CmdlineSingleChr(const ChrInfo* cip, LoadFilterLogFlags load_filter_log_flags) {
  if ((load_filter_log_flags & (kfLoadFilterLogAutosome | kfLoadFilterLogAutosomePar)) || (!cip->is_include_stack)) {
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

void GetExportfTargets(const char* const* argvk, uint32_t param_ct, ExportfFlags* exportf_flags_ptr, IdpasteFlags* exportf_id_paste_ptr, uint64_t* format_param_idxs_ptr) {
  // supports multiple formats
  uint64_t format_param_idxs = 0;
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
      } else if ((!strcmp(cur_modif2, "v")) || (!strcmp(cur_modif2, "-transpose"))) {
        cur_format = kfExportfAv;
      }
      break;
    case 'b':
      {
        const uint32_t cur_modif2_slen = strlen(cur_modif2);
        if (strequal_k(cur_modif2, "cf", cur_modif2_slen)) {
          cur_format = kfExportfBcf43;
        } else if (strequal_k(cur_modif2, "cf-4.2", cur_modif2_slen)) {
          cur_format = kfExportfBcf42;
        } else if (strequal_k(cur_modif2, "eagle", cur_modif2_slen)) {
          cur_format = kfExportfBeagle;
        } else if (strequal_k(cur_modif2, "eagle-nomap", cur_modif2_slen)) {
          cur_format = kfExportfBeagleNomap;
        } else if (strequal_k(cur_modif2, "gen-1.1", cur_modif2_slen) ||
                   strequal_k(cur_modif2, "gen_1.1", cur_modif2_slen)) {
          cur_format = kfExportfBgen11;
        } else if (strequal_k(cur_modif2, "gen-1.2", cur_modif2_slen) ||
                   strequal_k(cur_modif2, "gen_1.2", cur_modif2_slen)) {
          cur_format = kfExportfBgen12;
        } else if (strequal_k(cur_modif2, "gen-1.3", cur_modif2_slen) ||
                   strequal_k(cur_modif2, "gen_1.3", cur_modif2_slen)) {
          cur_format = kfExportfBgen13;
        } else if (strequal_k(cur_modif2, "imbam", cur_modif2_slen)) {
          cur_format = kfExportfBimbam;
        } else if (strequal_k(cur_modif2, "imbam-1chr", cur_modif2_slen)) {
          cur_format = kfExportfBimbam1chr;
        }
        break;
      }
    case 'c':
      if (!strcmp(cur_modif2, "ompound-genotypes")) {
        cur_format = kfExportfCompound;
      }
      break;
    case 'e':
      if (!strcmp(cur_modif2, "ig")) {
        cur_format = kfExportfEig;
      } else if (!strcmp(cur_modif2, "igt")) {
        cur_format = kfExportfEigt;
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
      {
        const uint32_t cur_modif2_slen = strlen(cur_modif2);
        if (strequal_k(cur_modif2, "gen", cur_modif2_slen)) {
          cur_format = kfExportfLgen;
        } else if (strequal_k(cur_modif2, "gen-ref", cur_modif2_slen)) {
          cur_format = kfExportfLgenRef;
        } else if (strequal_k(cur_modif2, "ist", cur_modif2_slen)) {
          cur_format = kfExportfList;
        }
        break;
      }
    case 'o':
      if (!strcmp(cur_modif2, "xford")) {
        cur_format = kfExportfOxGenV1;
      } else if (!strcmp(cur_modif2, "xford-v2")) {
        cur_format = kfExportfOxGenV2;
      }
      break;
    case 'p':
      {
        const uint32_t cur_modif2_slen = strlen(cur_modif2);
        if (strequal_k(cur_modif2, "ed", cur_modif2_slen)) {
          cur_format = kfExportfPed;
        } else if (strequal_k(cur_modif2, "hylip", cur_modif2_slen)) {
          cur_format = kfExportfPhylip;
        } else if (strequal_k(cur_modif2, "hylip-phased", cur_modif2_slen)) {
          cur_format = kfExportfPhylipPhased;
        }
        break;
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
      if ((!strcmp(cur_modif2, "ped")) || (!strcmp(cur_modif2, "ranspose"))) {
        cur_format = kfExportfTped;
      }
      break;
    case 'v':
      if ((cur_modif2[0] == 'c') && (cur_modif2[1] == 'f')) {
        if (!cur_modif2[2]) {
          cur_format = kfExportfVcf43;
        } else if (!strcmp(&(cur_modif2[2]), "-4.2")) {
          cur_format = kfExportfVcf42;
        } else if ((!strcmp(&(cur_modif2[2]), "-fid")) || (!strcmp(&(cur_modif2[2]), "-iid"))) {
          logprintf("Note: --export 'v%s' modifier is deprecated.  Use 'vcf' + 'id-paste=%s'.\n", cur_modif2, &(cur_modif2[3]));
          cur_format = kfExportfVcf43;
          *exportf_id_paste_ptr = (cur_modif2[3] == 'f')? kfIdpasteFid : kfIdpasteIid;
        }
      }
      break;
    }
    if (cur_format) {
      format_param_idxs |= 1LLU << param_idx;
      *exportf_flags_ptr |= cur_format;
    }
  }
  *format_param_idxs_ptr = format_param_idxs;
}

uint32_t VaridTemplateIsValid(const char* varid_str, const char* flagname_p) {
  const char* sptr = strchr(varid_str, '@');
  const char* sptr2 = strchr(varid_str, '#');
  if (unlikely((!sptr) || (!sptr2) || strchr(&(sptr[1]), '@') || strchr(&(sptr2[1]), '#'))) {
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
      if (unlikely((!sptr2) || strchr(&(sptr2[1]), '$') || ((first_allele_type_code + ctou32(sptr2[1])) != 99))) {
      VaridTemplateIsValid_dollar_error:
        logerrprintfww("Error: The --%s template string requires either no instances of '$', exactly one instance of '$r' and/or '$a', or exactly one '$1' and one '$2'.\n", flagname_p);
        return 0;
      }
    } else {
      first_allele_type_code &= 0xdf;
      if (unlikely((first_allele_type_code != 65) && (first_allele_type_code != 82))) {
        goto VaridTemplateIsValid_dollar_error;
      }
      sptr2 = strchr(sptr2, '$');
      if (sptr2) {
        const uint32_t second_allele_type_code = ctou32(*(++sptr2)) & 0xdf;
        if (unlikely(((first_allele_type_code + second_allele_type_code) != 147) || strchr(sptr2, '$'))) {
          goto VaridTemplateIsValid_dollar_error;
        }
      }
    }
  }
  return 1;
}

BoolErr ParseFreqSelector(const char* mode_str, const char* flagname_p, FreqFilterMode* modep) {
  // assume preinitialization to kFreqFilterNonmajor
  if (mode_str[0] == ':') {
    ++mode_str;
  }
  const uint32_t mode_slen = strlen(mode_str);
  if (strequal_k(mode_str, "nref", mode_slen)) {
    *modep = kFreqFilterNref;
  } else if (strequal_k(mode_str, "alt1", mode_slen)) {
    *modep = kFreqFilterAlt1;
  } else if (strequal_k(mode_str, "minor", mode_slen)) {
    *modep = kFreqFilterMinor;
  } else if (unlikely(!strequal_k(mode_str, "nonmajor", mode_slen))) {
    snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s mode '%s'.\n", flagname_p, mode_str);
    return 1;
  }
  return 0;
}


static_assert(sizeof(int) == sizeof(int32_t), "main() assumes int and int32_t are synonymous.");
static_assert(!kChrOffsetX, "--autosome-num/--chr-set/--cow/etc. assume kChrOffsetX == 0.");
static_assert(kChrOffsetY == 1, "--chr-set/--cow/... assume kChrOffsetY == 1.");
static_assert(kChrOffsetXY == 2, "--chr-set/--cow/... assume kChrOffsetXY == 2.");
static_assert(kChrOffsetMT == 3, "--chr-set/--cow/... assume kChrOffsetMT == 3.");
#ifdef __cplusplus
}  // namespace plink2
#endif

#if defined(CPU_CHECK_SSE42) || defined(CPU_CHECK_AVX2)
int RealMain(int argc, char** argv) {
#else
int main(int argc, char** argv) {
#endif
#ifdef __cplusplus
  using namespace plink2;
#endif

#ifdef __APPLE__
  #if __x86_64__
    fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);
  #endif
#else
#  if defined __LP64__ && defined __x86_64__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#  endif
#endif
  // special case, since it may dump to stdout
  if (argc > 1) {
    const char* argv1 = argv[1];
    const uint32_t argv1_slen = strlen(argv1);
    if (strequal_k(argv1, "--zst-decompress", argv1_slen) ||
        strequal_k(argv1, "-zst-decompress", argv1_slen) ||
        strequal_k(argv1, "--zd", argv1_slen) ||
        strequal_k(argv1, "-zd", argv1_slen)) {
      if (unlikely(argc == 2)) {
        fprintf(stderr, "Error: Missing %s argument.\n", argv[1]);
        return S_CAST(uint32_t, kPglRetInvalidCmdline);
      }
      for (int ii = 2; ii != argc; ++ii) {
        if (unlikely(IsCmdlineFlag(argv[S_CAST(uint32_t, ii)]))) {
          fprintf(stderr, "Error: %s cannot be used with other flags.\n", argv[1]);
          return S_CAST(uint32_t, kPglRetInvalidCmdline);
        }
      }
      if (unlikely(argc > 4)) {
        fprintf(stderr, "Error: %s accepts at most 2 arguments.\n", argv[1]);
        return S_CAST(uint32_t, kPglRetInvalidCmdline);
      }
      return S_CAST(uint32_t, ZstDecompress(argv[2], (argc == 4)? argv[3] : nullptr));
    }
  }

  unsigned char* bigstack_ua = nullptr;
  Plink2CmdlineMeta pcm;
  PreinitPlink2CmdlineMeta(&pcm);
  PmergeInfo pmerge_info;
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
  pc.load_filter_log_flags = kfLoadFilterLog0;
  pc.pginame = nullptr;
  pc.var_filter_exceptions_flattened = nullptr;
  pc.varid_template_str = nullptr;
  pc.varid_multi_template_str = nullptr;
  pc.varid_multi_nonsnp_template_str = nullptr;
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
  pc.extract_intersect_fnames = nullptr;
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
  pc.fa_fname = nullptr;
  pc.king_table_subset_fname = nullptr;
  pc.king_table_require_fnames = nullptr;
  pc.require_info_flattened = nullptr;
  pc.require_no_info_flattened = nullptr;
  pc.keep_col_match_fname = nullptr;
  pc.keep_col_match_flattened = nullptr;
  pc.keep_col_match_name = nullptr;
  pc.ref_allele_flag = nullptr;
  pc.alt_allele_flag = nullptr;
  pc.update_chr_flag = nullptr;
  pc.update_map_flag = nullptr;
  pc.update_name_flag = nullptr;
  pc.update_sample_ids_fname = nullptr;
  pc.update_parental_ids_fname = nullptr;
  pc.recover_var_ids_fname = nullptr;
  pc.vscore_fname = nullptr;
  pc.indep_preferred_fname = nullptr;
  pc.not_pheno_flattened = nullptr;
  pc.not_covar_flattened = nullptr;
  pc.indv_str = nullptr;
  pc.rename_chrs_fname = nullptr;
  pc.zero_cluster_fname = nullptr;
  pc.zero_cluster_phenoname = nullptr;
  InitRangeList(&pc.snps_range_list);
  InitRangeList(&pc.exclude_snps_range_list);
  InitRangeList(&pc.pheno_range_list);
  InitRangeList(&pc.covar_range_list);
  InitRangeList(&pc.vscore_col_idx_range_list);
  InitUpdateAlleles(&pc.update_alleles_info);
  InitUpdateSex(&pc.update_sex_info);
  InitLd(&pc.ld_info);
  InitSdiff(&pc.sdiff_info);
  InitGlm(&pc.glm_info);
  InitScore(&pc.score_info);
  InitFst(&pc.fst_info);
  InitPmerge(&pmerge_info);
  InitPgenDiff(&pc.pgen_diff_info);
  InitCmpExpr(&pc.keep_if_expr);
  InitCmpExpr(&pc.remove_if_expr);
  InitCmpExpr(&pc.extract_if_info_expr);
  InitCmpExpr(&pc.exclude_if_info_expr);
  InitExtractColCond(&pc.extract_col_cond_info);
  InitExportf(&pc.exportf_info);
  InitGwasSsf(&pc.gwas_ssf_info);
  InitClump(&pc.clump_info);
  InitVcor(&pc.vcor_info);
  InitPhenoSvd(&pc.pheno_svd_info);
  InitCheckSex(&pc.check_sex_info);
  InitMendel(&pc.mendel_info);
  InitFlip(&pc.flip_info);
  GenDummyInfo gendummy_info;
  InitGenDummy(&gendummy_info);
  AdjustFileInfo adjust_file_info;
  InitAdjust(&pc.adjust_info, &adjust_file_info);
  ChrInfo chr_info;
  if (unlikely(InitChrInfo(&chr_info))) {
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
      if (reterr == kPglRetHelp) {
        reterr = kPglRetSuccess;
      }
      goto main_ret_NOLOG;
    }
    if (!flag_ct) {
      goto main_ret_NULL_CALC_0;
    }
    // + (kBytesPerWord - 1) to support overread during sort
    if (unlikely(pgl_malloc(flag_ct * kMaxFlagBlen + kBytesPerWord - 1, &pcm.flag_buf) ||
                 pgl_malloc(flag_ct * sizeof(int32_t), &pcm.flag_map))) {
      goto main_ret_NOMEM_NOLOG2;
    }

    // No modifications to argv past this point.
    const char* const* argvk = TO_CONSTCPCONSTP(argv);

    char* flagname_write_iter = pcm.flag_buf;
    uint32_t cur_flag_idx = 0;
    for (arg_idx = first_arg_idx; arg_idx != S_CAST(uint32_t, argc); ++arg_idx) {
      flagname_p = IsCmdlineFlagStart(argvk[arg_idx]);
      if (flagname_p) {
        const uint32_t flag_slen = strlen(flagname_p);
        switch (*flagname_p) {
        case 'F':
          if (strequal_k(flagname_p, "Fst", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "fst");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'a':
          if (strequal_k(flagname_p, "aec", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "allow-extra-chr");
          } else if (strequal_k(flagname_p, "autosome-xy", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "autosome-par");
          } else if (strequal_k(flagname_p, "a1-allele", flag_slen)) {
            fputs("Warning: --a1-allele flag deprecated.  Use --alt1-allele instead.\n", stderr);
            g_stderr_written_to = 1;
            snprintf(flagname_write_iter, kMaxFlagBlen, "alt1-allele");
          } else if (strequal_k(flagname_p, "a2-allele", flag_slen)) {
            fputs("Warning: --a2-allele flag deprecated.  Use --ref-allele instead.\n", stderr);
            g_stderr_written_to = 1;
            snprintf(flagname_write_iter, kMaxFlagBlen, "ref-allele");
          } else if (strequal_k(flagname_p, "allele-1234", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "allele1234");
          } else if (strequal_k(flagname_p, "alleleacgt", flag_slen) ||
                     strequal_k(flagname_p, "allele-ACGT", flag_slen) ||
                     strequal_k(flagname_p, "allele-acgt", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "alleleACGT");
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
          if (strequal_k(flagname_p, "covarExcludeList", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "not-covar");
          } else if (strequal_k(flagname_p, "clump-snp-field", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "clump-id-field");
          } else if (strequal_k(flagname_p, "clump-field", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "clump-p-field");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'e':
          if (strequal_k(flagname_p, "extract-if", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-if-info");
          } else if (strequal_k(flagname_p, "exclude-if", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "exclude-if-info");
          } else if (strequal_k(flagname_p, "extract-fcol", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-col-cond");
          } else if (strequal_k(flagname_p, "extract-fcol-match", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-col-cond-match");
          } else if (strequal_k(flagname_p, "extract-fcol-mismatch", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-col-cond-mismatch");
          } else if (strequal_k(flagname_p, "extract-fcol-substr", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-col-cond-substr");
          } else if (strequal_k(flagname_p, "extract-fcol-min", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-col-cond-min");
          } else if (strequal_k(flagname_p, "extract-fcol-max", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-col-cond-max");
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
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-col-match");
          } else if (strequal_k(flagname_p, "fill-missing-a2", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "fill-missing-with-ref");
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
          } else if (strequal_k(flagname_p, "keep-fcol", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-col-match");
          } else if (strequal_k(flagname_p, "keep-fcol-name", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-col-match-name");
          } else if (strequal_k(flagname_p, "keep-fcol-num", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "keep-col-match-num");
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
          } else if (strequal_k(flagname_p, "min-af", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "maf");
          } else if (strequal_k(flagname_p, "max-af", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "max-maf");
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
          } else if (strequal_k(flagname_p, "merge-max-allele-ct", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "merge-max-alleles");
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
          } else if (strequal_k(flagname_p, "norm", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "normalize");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'p':
          if (strequal_k(flagname_p, "prune", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "require-pheno");
          } else if (strequal_k(flagname_p, "phenoExcludeList", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "not-pheno");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'q':
          if (strequal_k(flagname_p, "qual-scores", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-col-cond");
          } else if (strequal_k(flagname_p, "qual-threshold", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-col-cond-min");
          } else if (strequal_k(flagname_p, "qual-max-threshold", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "extract-col-cond-max");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'r':
          if (strequal_k(flagname_p, "recode", flag_slen)) {
            // don't bother translating no-modifier --recode to "--export ped",
            // just let it fail
            snprintf(flagname_write_iter, kMaxFlagBlen, "export");
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
          } else if (strequal_k(flagname_p, "recode-allele", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "export-allele");
          } else {
            goto main_flag_copy;
          }
          break;
        case 's':
          if (strequal_k(flagname_p, "sdiff", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "sample-diff");
          } else if (strequal_k(flagname_p, "set-hh-missing", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "set-invalid-haploid-missing");
          } else {
            goto main_flag_copy;
          }
          break;
        case 't':
          if (strequal_k(flagname_p, "thread-num", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "threads");
          } else if (strequal_k(flagname_p, "tfam", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "psam");
          } else {
            goto main_flag_copy;
          }
          break;
        case 'v':
          if (strequal_k(flagname_p, "vcf-filter", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "var-filter");
          } else if (strequal_k(flagname_p, "vcf-min-qual", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "var-min-qual");
          } else if (strequal_k(flagname_p, "vscore", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "variant-score");
          } else if (strequal_k(flagname_p, "vcf-idspace-to", flag_slen)) {
            snprintf(flagname_write_iter, kMaxFlagBlen, "idspace-to");
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
    memcpy_k(outname, "plink2", 6);
    char* outname_end = nullptr;
    char range_delim;
    uint32_t strict_extra_chr;
    int32_t known_procs;
    reterr = CmdlineParsePhase2(ver_str, errstr_append, argvk, 6, kMaxFlagBlen, argc, flag_ct, &pcm, outname, &outname_end, &range_delim, &strict_extra_chr, &known_procs, &pc.max_thread_ct);
    if (unlikely(reterr)) {
      goto main_ret_NOLOG;
    }

    char pgenname[kPglFnamesize];
    char psamname[kPglFnamesize];
    char pvarname[kPglFnamesize];
    pgenname[0] = '\0';
    psamname[0] = '\0';
    pvarname[0] = '\0';
    snprintf(pc.missing_catname, kMaxMissingPhenostrBlen, "NONE");
    snprintf(pc.output_missing_pheno, kMaxMissingPhenostrBlen, "NA");
    snprintf(pc.legacy_output_missing_pheno, kMaxMissingPhenostrBlen, "-9");
    cur_flag_idx = 0;
    pc.command_flags1 = kfCommand10;
    // uint64_t command_flags2 = 0;
    pc.misc_flags = strict_extra_chr? kfMiscProhibitExtraChr : kfMisc0;
    pc.pvar_psam_flags = kfPvarPsam0;
    pc.sample_sort_mode = kSort0;
    pc.sort_vars_mode = kSort0;
    pc.grm_flags = kfGrm0;
    pc.grm_sparse_cutoff = -DBL_MAX;
    pc.pca_flags = kfPca0;
    pc.write_covar_flags = kfWriteCovar0;
    pc.pheno_transform_flags = kfPhenoTransform0;
    pc.fa_flags = kfFa0;
    pc.fam_cols = kfFamCol13456;
    pc.king_flags = kfKing0;
    pc.king_cutoff = -1;
    pc.king_table_filter = -DBL_MAX;
    pc.freq_rpt_flags = kfAlleleFreq0;
    pc.missing_rpt_flags = kfMissingRpt0;
    pc.geno_counts_flags = kfGenoCounts0;
    pc.hardy_flags = kfHardy0;
    pc.het_flags = kfHet0;
    pc.sample_counts_flags = kfSampleCounts0;
    pc.recover_var_ids_flags = kfRecoverVarIds0;
    pc.vscore_flags = kfVscore0;
    pc.allele_alphanum_flags = kfAlleleAlphanum0;
    pc.select_sid_missingness_mode = kSelectSidMissingness0;
    pc.select_sid_tiebreak_mode = kSelectSidTiebreak0;
    pc.rmdup_mode = kRmDup0;
    for (uint32_t uii = 0; uii != 4; ++uii) {
      pc.filter_modes[uii] = kFreqFilterNonmajor;
    }
    pc.aperm.min = 6;
    pc.aperm.max = 1000000;
    pc.aperm.alpha = 0.0;
    pc.aperm.beta = 0.0001;
    pc.aperm.init_interval = 1.0;
    pc.aperm.interval_slope = 0.001;
    pc.ci_size = 0.0;

    // Default value is 1638 = 32768 / 20, and that's applied to imported
    // dosages when --hard-call-threshold is not specified.
    // However, when --make-[b]pgen is run on a dosage-containing dataset,
    // explicit --hard-call-threshold will cause the hardcall set to be
    // regenerated, and that won't happen without --hard-call-threshold.  So we
    // need to distinguish between --hard-call-threshold 0.1 and no flag.
    pc.hard_call_thresh = UINT32_MAX;

    pc.dosage_erase_thresh = 0;
    pc.ln_pfilter = kLnPvalError;  // make --pfilter 1 still filter out NAs
    pc.output_min_ln = -DBL_MAX;
    pc.vif_thresh = 50.0;
    pc.mind_thresh = 1.0;
    pc.geno_thresh = 1.0;
    pc.hwe_ln_thresh = -DBL_MAX;
    pc.hwe_sample_size_term = -1;
    pc.mach_r2_min = 0.0;
    pc.mach_r2_max = 0.0;
    pc.minimac3_r2_min = 0.0;
    pc.minimac3_r2_max = 0.0;
    pc.af_pseudocount = 0.0;
    pc.min_maf = 0.0;
    pc.max_maf = 1.0;
    pc.thin_keep_prob = 1.0;
    pc.thin_keep_sample_prob = 1.0;
    pc.min_allele_ddosage = 0;
    pc.max_allele_ddosage = (~0LLU);
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
    pc.mwithin_val = 1;
    pc.min_bp_space = 0;
    pc.thin_keep_ct = UINT32_MAX;
    pc.thin_keep_sample_ct = UINT32_MAX;
    pc.keep_col_match_num = 0;
    pc.filter_min_allele_ct = 0;
    pc.filter_max_allele_ct = UINT32_MAX;
    pc.bed_border_bp = 0;
    pc.input_missing_geno_char = '0';
    pc.output_missing_geno_char = '.';
    pc.legacy_output_missing_geno_char = '0';
    double import_dosage_certainty = 0.0;
    int32_t vcf_min_gq = -1;
    int32_t vcf_min_dp = -1;
    int32_t vcf_max_dp = 0x7fffffff;
    uintptr_t malloc_size_mib = 0;
    LoadParams load_params = kfLoadParams0;
    Xload xload = kfXload0;
    uint32_t rseed_ct = 0;
    MakePlink2Flags make_plink2_flags = kfMake0;
    OxfordImportFlags oxford_import_flags = kfOxfordImport0;
    VcfHalfCall vcf_half_call = kVcfHalfCallDefault;
    char id_delim = '\0';
    char idspace_to = '\0';
    ImportFlags import_flags = kfImport0;
    uint32_t delete_pmerge_result = 0;
    uint32_t aperm_present = 0;
    uint32_t clump_log10_p1_present = 0;
    uint32_t clump_log10_p2_present = 0;
    uint32_t score_col_nums_present = 0;
    uint32_t r2_required = 0;
    uint32_t permit_multiple_inclusion_filters = 0;
    uint32_t memory_require = 0;
#ifdef USE_MKL
    uint32_t mkl_native = 0;
#endif
    uint32_t randmem = 0;
    uint32_t pmerge_required = 0;
    uint32_t allow_misleading_out_arg = 0;
    uint32_t allow_normalize_with_split = 0;
    uint32_t import_max_allele_ct = 0x7ffffffe;
    Plink1DosageInfo plink1_dosage_info;
    InitPlink1Dosage(&plink1_dosage_info);
    do {
      flagname_p = &(pcm.flag_buf[cur_flag_idx * kMaxFlagBlen]);
      const char* flagname_p2 = &(flagname_p[1]);
      arg_idx = pcm.flag_map[cur_flag_idx];
      uint32_t param_ct = GetParamCt(argvk, argc, arg_idx);
      switch (*flagname_p) {
      case '1':
        if (likely(*flagname_p2 == '\0')) {
          pc.misc_flags |= kfMiscAffection01;
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'a':
        if (strequal_k_unsafe(flagname_p2, "llow-extra-chr")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (unlikely(!strequal_k_unsafe(cur_modif, "0"))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --allow-extra-chr argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            chr_info.zero_extra_chrs = 1;
          }
          pc.misc_flags &= ~kfMiscProhibitExtraChr;
        } else if (strequal_k_unsafe(flagname_p2, "utosome")) {
          chr_info.is_include_stack = 1;
          pc.load_filter_log_flags |= kfLoadFilterLogAutosome;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "utosome-par")) {
          if (unlikely(pc.load_filter_log_flags & kfLoadFilterLogAutosome)) {
            logerrputs("Error: --autosome-par cannot be used with --autosome.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.is_include_stack = 1;
          pc.load_filter_log_flags |= kfLoadFilterLogAutosomePar;
          goto main_param_zero;
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "llow-no-samples"))) {
          logerrputs("Error: --allow-no-samples is retired.  (If you are performing a set of\noperations which doesn't require sample information, the sample file won't be\nloaded at all.)\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "llow-no-vars"))) {
          logerrputs("Error: --allow-no-vars is retired.  (If you are performing a set of operations\nwhich doesn't require variant information, the variant file won't be loaded at\nall.)\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (strequal_k_unsafe(flagname_p2, "djust")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4))) {
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
            } else if (unlikely(strequal_k(cur_modif, "qq-plot", cur_modif_slen))) {
              logerrputs("Error: 'qq-plot' modifier retired.  Use e.g. \"--adjust cols=+qq\" instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (likely(StrStartsWith(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(pc.adjust_info.flags & kfAdjustColAll)) {
                logerrputs("Error: Multiple --adjust cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0maybeprovref\0provref\0a1\0unadj\0gc\0qq\0bonf\0holm\0sidakss\0sidaksd\0fdrbh\0fdrby\0", "adjust", kfAdjustColChrom, kfAdjustColDefault, 1, &pc.adjust_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --adjust argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.adjust_info.flags & kfAdjustColAll)) {
            pc.adjust_info.flags |= kfAdjustColDefault;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-file")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 7))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &adjust_file_info.fname);
          if (unlikely(reterr)) {
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
              if (unlikely(adjust_file_info.base.flags & kfAdjustColAll)) {
                logerrputs("Error: Multiple --adjust-file cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0maybeprovref\0provref\0a1\0unadj\0gc\0qq\0bonf\0holm\0sidakss\0sidaksd\0fdrbh\0fdrby\0", "adjust-file", kfAdjustColChrom, kfAdjustColDefault, 1, &adjust_file_info.base.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "test=", cur_modif_slen)) {
              reterr = CmdlineAllocString(&(cur_modif[5]), "--adjust-file test=", kMaxIdSlen, &adjust_file_info.test_name);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (likely(strequal_k(cur_modif, "input-log10", cur_modif_slen))) {
              adjust_file_info.base.flags |= kfAdjustInputLog10;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --adjust-file argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(adjust_file_info.base.flags & kfAdjustColAll)) {
            adjust_file_info.base.flags |= kfAdjustColDefault;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-chr-field")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &adjust_file_info.chr_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-alt-field")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &adjust_file_info.alt_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-provref-field")) {
          if (unlikely(!adjust_file_info.fname)) {
            logerrputs("Error: --adjust-provref-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &adjust_file_info.provref_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-a1-field")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &adjust_file_info.a1_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-pos-field")) {
          if (unlikely(!adjust_file_info.fname)) {
            logerrputs("Error: --adjust-pos-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &adjust_file_info.pos_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-id-field")) {
          if (unlikely(!adjust_file_info.fname)) {
            logerrputs("Error: --adjust-id-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &adjust_file_info.id_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-ref-field")) {
          if (unlikely(!adjust_file_info.fname)) {
            logerrputs("Error: --adjust-ref-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &adjust_file_info.ref_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-test-field")) {
          if (unlikely(!adjust_file_info.fname)) {
            logerrputs("Error: --adjust-test-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &adjust_file_info.test_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "djust-p-field")) {
          if (unlikely(!adjust_file_info.fname)) {
            logerrputs("Error: --adjust-p-field must be used with --adjust-file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &adjust_file_info.p_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "perm")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 6))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(ScanPosintDefcapx(cur_modif, &pc.aperm.min))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm min permutation count '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          ++pc.aperm.min;
          if (param_ct > 1) {
            cur_modif = argvk[arg_idx + 2];
            if (unlikely(ScanPosintCappedx(cur_modif, kApermMax, &pc.aperm.max))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm max permutation count '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (unlikely(pc.aperm.min >= pc.aperm.max)) {
            logerrputs("Error: --aperm min permutation count must be smaller than max.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          aperm_present = 1;
          if (param_ct > 2) {
            cur_modif = argvk[arg_idx + 3];
            if (unlikely(!ScantokDouble(cur_modif, &pc.aperm.alpha))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm alpha threshold '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (param_ct > 3) {
              cur_modif = argvk[arg_idx + 4];
              if (unlikely((!ScantokDouble(cur_modif, &pc.aperm.beta)) || (pc.aperm.beta <= 0.0))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm beta '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (param_ct > 4) {
                cur_modif = argvk[arg_idx + 5];
                if (unlikely(!ScantokDouble(cur_modif, &pc.aperm.init_interval))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm initial pruning interval '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                if (unlikely((pc.aperm.init_interval < 1.0) || (pc.aperm.init_interval > 1000000.0))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm initial pruning interval '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                if (param_ct == 6) {
                  cur_modif = argvk[arg_idx + 6];
                  if (unlikely((!ScantokDouble(cur_modif, &pc.aperm.interval_slope)) || (pc.aperm.interval_slope < 0.0) || (pc.aperm.interval_slope > 1.0))) {
                    snprintf(g_logbuf, kLogbufSize, "Error: Invalid --aperm pruning interval slope '%s'.\n", cur_modif);
                    goto main_ret_INVALID_CMDLINE_WWA;
                  }
                }
              }
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "utosome-num")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t autosome_ct;
          if (unlikely(ScanPosintCappedx(cur_modif, kMaxChrTextnum, &autosome_ct))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --autosome-num argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          // see plink2_common FinalizeChrset()
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.autosome_ct = autosome_ct;
          // assumes first code is X
          chr_info.xymt_codes[0] = autosome_ct + 1;
          for (uint32_t xymt_idx = 1; xymt_idx != kChrOffsetCt; ++xymt_idx) {
            // bugfix: this needs to be UINT32_MAXM1, not UINT32_MAX, for
            // GetChrCode() to work properly
            chr_info.xymt_codes[xymt_idx] = UINT32_MAXM1;
          }
          chr_info.haploid_mask[0] = 0;
          SetBit(autosome_ct + 1, chr_info.haploid_mask);
        } else if (strequal_k_unsafe(flagname_p2, "lt-allele") || strequal_k_unsafe(flagname_p2, "lt1-allele")) {
          if (unlikely(pc.alt_allele_flag)) {
            logerrputs("Error: --alt-allele cannot be used with --alt1-allele.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 5))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* const* sources = &(argvk[arg_idx + 1]);
          if (!strcmp(sources[0], "force")) {
            --param_ct;
            if (unlikely(!param_ct)) {
              logerrputs("Error: Invalid --alt[1]-allele argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.misc_flags |= kfMiscAltAlleleForce;
            ++sources;
          }
          reterr = Alloc2col(sources, flagname_p, param_ct, &pc.alt_allele_flag);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          if (flagname_p2[2] == '1') {
            pc.misc_flags |= kfMiscAlt1Allele;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "f-pseudocount")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double dxx;
          if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &dxx)) || (dxx < 0.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --af-pseudocount argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.af_pseudocount = dxx;
        } else if (strequal_k_unsafe(flagname_p2, "llow-misleading-out-arg")) {
          allow_misleading_out_arg = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "llow-normalize-with-split")) {
          allow_normalize_with_split = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "c-founders")) {
          pc.misc_flags |= kfMiscAcFounders;
        } else if (strequal_k_unsafe(flagname_p2, "llele1234")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            if (unlikely(strcmp(argvk[arg_idx + 1], "multichar"))) {
              logerrputs("Error: Invalid --allele1234 argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.allele_alphanum_flags |= kfAlleleAlphanumMultichar;
          }
          pc.allele_alphanum_flags |= kfAlleleAlphanum1234;
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "lleleACGT")) {
          if (unlikely(pc.allele_alphanum_flags != kfAlleleAlphanum0)) {
            logerrputs("Error: --allele1234 cannot be used with --alleleACGT.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            if (unlikely(strcmp(argvk[arg_idx + 1], "multichar"))) {
              logerrputs("Error: Invalid --alleleACGT argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.allele_alphanum_flags |= kfAlleleAlphanumMultichar;
          }
          pc.allele_alphanum_flags |= kfAlleleAlphanumAcgt;
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "ssoc"))) {
          logerrputs("Error: --assoc is retired.  Use --glm instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (likely(strequal_k_unsafe(flagname_p2, "llow-no-sex"))) {
          logputs("Note: --allow-no-sex no longer has any effect.  (Missing-sex samples are\nautomatically excluded from association analysis when sex is a covariate, and\ntreated normally otherwise.)\n");
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'b':
        if (strequal_k_unsafe(flagname_p2, "file")) {
          if (unlikely(xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "vzs", &fname_modif_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* fname_prefix = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname_prefix);
          // pathological case bugfix (1 Feb 2018): need to subtract 1 more, to
          // avoid buffer overflow in the case we rename and append '~'.
          if (unlikely(slen > (kPglFnamesize - 10))) {
            // could use kPglFnamesize - 2 - 3 * param_ct, but that's pointless
            logerrputs("Error: --bfile argument too long.\n");
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
          if (unlikely(load_params || xload)) {
            // currently only possible with --bcf, --bfile
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "vzs", &fname_modif_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* fname_prefix = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname_prefix);
          if (unlikely(slen > (kPglFnamesize - 10))) {
            logerrputs("Error: --bpfile argument too long.\n");
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
          logerrputs("Error: --biallelic-only is retired.  Use e.g. \"--max-alleles 2\" instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "cf")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 2) {
            const char* cur_modif = argvk[arg_idx + 2];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (unlikely(!StrStartsWith(cur_modif, "dosage=", cur_modif_slen))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --bcf argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            reterr = CmdlineAllocString(&(cur_modif[strlen("dosage=")]), argvk[arg_idx], 4095, &vcf_dosage_import_field);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            const uint32_t vcf_dosage_import_field_slen = strlen(vcf_dosage_import_field);
            if (unlikely(!((strequal_k(vcf_dosage_import_field, "GP-force", vcf_dosage_import_field_slen)) || IsAlphanumeric(vcf_dosage_import_field)))) {
              logerrputs("Error: --bcf dosage= argument is not alphanumeric.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            if (unlikely(strequal_k(vcf_dosage_import_field, "GT", vcf_dosage_import_field_slen))) {
              logerrputs("Error: --bcf dosage= argument cannot be 'GT'.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            // rather not worry about string-index 0 here
            if (unlikely(strequal_k(vcf_dosage_import_field, "PASS", vcf_dosage_import_field_slen))) {
              logerrputs("Error: --bcf dosage= argument cannot be 'PASS'.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_modif);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --bcf filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_modif, slen + 1);
          xload = kfXloadBcf;
        } else if (strequal_k_unsafe(flagname_p2, "gen")) {
          if (unlikely(load_params || xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(param_ct == 1)) {
            logerrputs("Error: --bgen now requires a REF/ALT mode ('ref-first', 'ref-last', or\n'ref-unknown').  As of this writing, raw UK Biobank files are ref-first, while\nolder and PLINK-exported BGEN files are more likely to be ref-last.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "snpid-chr", cur_modif_slen)) {
              oxford_import_flags |= kfOxfordImportBgenSnpIdChr;
            } else if (strequal_k(cur_modif, "ref-first", cur_modif_slen)) {
              oxford_import_flags |= kfOxfordImportRefFirst;
            } else if (strequal_k(cur_modif, "ref-unknown", cur_modif_slen)) {
              oxford_import_flags |= kfOxfordImportRefUnknown;
            } else if (strequal_k(cur_modif, "ref-last", cur_modif_slen)) {
              oxford_import_flags |= kfOxfordImportRefLast;
            } else if (likely(strequal_k(cur_modif, "ref-second", cur_modif_slen))) {
              logerrputs("Warning: --bgen 'ref-second' modifier is deprecated.  Use 'ref-last' instead.\n");
              oxford_import_flags |= kfOxfordImportRefLast;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --bgen argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(oxford_import_flags & kfOxfordImportRefAll)) {
            logerrputs("Error: No --bgen REF/ALT mode specified ('ref-first', 'ref-last', or\n'ref-unknown'); this is now required.\n");
          }
          const char* cur_fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --bgen filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_fname, slen + 1);
          xload = kfXloadOxBgen;
        } else if (strequal_k_unsafe(flagname_p2, "p-space")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(ScanPosintDefcapx(cur_modif, &pc.min_bp_space))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --bp-space minimum bp distance '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "ed-border-bp")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double dxx;
          if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &dxx)) || (dxx < 0.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --bed-border-bp argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (dxx > 2147483646) {
            pc.bed_border_bp = 0x7ffffffe;
          } else {
            pc.bed_border_bp = S_CAST(int32_t, dxx * (1 + kSmallEpsilon));
          }
        } else if (strequal_k_unsafe(flagname_p2, "ed-border-kb")) {
          if (unlikely(pc.bed_border_bp)) {
            logerrputs("Error: --bed-border-kb cannot be used with --bed-border-bp.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double dxx;
          if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &dxx)) || (dxx < 0.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --bed-border-kb argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (dxx > 2147483.646) {
            pc.bed_border_bp = 0x7ffffffe;
          } else {
            pc.bed_border_bp = S_CAST(int32_t, dxx * 1000 * (1 + kSmallEpsilon));
          }
        } else if (strequal_k_unsafe(flagname_p2, "ad-ld")) {
          pc.misc_flags |= kfMiscAllowBadLd;
          goto main_param_zero;
        } else if (likely(strequal_k_unsafe(flagname_p2, "ad-freqs"))) {
          pc.misc_flags |= kfMiscAllowBadFreqs;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "merge")) {
          logerrputs("Error: --bmerge is retired.  Use --pmerge instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'c':
        if (strequal_k_unsafe(flagname_p2, "hr")) {
          if (unlikely(pc.load_filter_log_flags & (kfLoadFilterLogAutosome | kfLoadFilterLogAutosomePar))) {
            logerrputs("Error: --chr cannot be used with --autosome[-par].\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = ParseChrRanges(&(argvk[arg_idx]), flagname_p, errstr_append, param_ct, (pc.misc_flags / kfMiscProhibitExtraChr) & 1, 0, '-', &chr_info, chr_info.chr_mask);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          chr_info.is_include_stack = 1;
          pc.load_filter_log_flags |= kfLoadFilterLogChr;
        } else if (strequal_k_unsafe(flagname_p2, "ovar")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "iid-only", &fname_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.misc_flags |= kfMiscCovarIidOnly;
          }
          reterr = AllocFname(argvk[arg_idx + fname_idx], flagname_p, &pc.covar_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ovar-col-nums")) {
          // requires --covar or --pheno, but --pheno hasn't been parsed yet so
          // we don't enforce the condition here
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.covar_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.misc_flags |= kfMiscCovarColNums;
        } else if (strequal_k_unsafe(flagname_p2, "ovar-name")) {
          if (unlikely(pc.covar_range_list.name_ct)) {
            logerrputs("Error: --covar-name can't be used with --covar-col-nums.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // can now be used without --covar
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.covar_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "onst-fid")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(param_ct? argvk[arg_idx + 1] : "0", argvk[arg_idx], kMaxIdSlen, &const_fid);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "i")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(!ScantokDouble(argvk[arg_idx + 1], &pc.ci_size))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --ci argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely((pc.ci_size < 0.01) || (pc.ci_size >= 1.0))) {
            logerrputs("Error: --ci confidence interval size must be in [0.01, 1).\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ondition")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.glm_info.condition_varname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "dominant", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmConditionDominant;
            } else if (strequal_k(cur_modif, "recessive", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmConditionRecessive;
            } else if (likely(strequal_k(cur_modif, "multiallelic", cur_modif_slen) ||
                              strequal_k(cur_modif, "m", cur_modif_slen))) {
              pc.glm_info.flags |= kfGlmConditionMultiallelic;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --condition argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (unlikely((pc.glm_info.flags & (kfGlmConditionDominant | kfGlmConditionRecessive)) == (kfGlmConditionDominant | kfGlmConditionRecessive))) {
            logerrputs("Error: --condition 'dominant' and 'recessive' modifiers can't be used together.\n");
            goto main_ret_INVALID_CMDLINE;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ondition-list")) {
          if (unlikely(pc.glm_info.condition_varname)) {
            logerrputs("Error: --condition-list cannot be used with --condition.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.glm_info.condition_list_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "dominant", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmConditionDominant;
            } else if (strequal_k(cur_modif, "recessive", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmConditionRecessive;
            } else if (likely(strequal_k(cur_modif, "multiallelic", cur_modif_slen) ||
                              strequal_k(cur_modif, "m", cur_modif_slen))) {
              pc.glm_info.flags |= kfGlmConditionMultiallelic;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --condition-list argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (unlikely((pc.glm_info.flags & (kfGlmConditionDominant | kfGlmConditionRecessive)) == (kfGlmConditionDominant | kfGlmConditionRecessive))) {
            logerrputs("Error: --condition-list 'dominant' and 'recessive' modifiers can't be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ow")) {
          if (unlikely(chr_info.chrset_source)) {
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
          if (unlikely(chr_info.chrset_source)) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 5))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          int32_t signed_autosome_ct;
          if (unlikely(ScanInt32x(cur_modif, &signed_autosome_ct) || (!signed_autosome_ct))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --chr-set argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(abs_i32(signed_autosome_ct) > kMaxChrTextnum)) {
            if (abs_i32(signed_autosome_ct) <= kMaxContigs) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --chr-set argument '%s' (max %u). (If you're dealing with a draft assembly with lots of contigs, rather than actual autosomes, PLINK can handle it if you name your contigs e.g. 'contig1', 'contig2', etc.)\n", cur_modif, kMaxChrTextnum);
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --chr-set argument '%s' (max %u).\n", cur_modif, kMaxChrTextnum);
            }
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          // see plink2_common FinalizeChrset()
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.haploid_mask[0] = 0;
          if (signed_autosome_ct < 0) {
            // haploid
            if (unlikely(param_ct > 1)) {
              logerrputs("Error: --chr-set does not accept multiple arguments in haploid mode.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            const uint32_t autosome_ct = -signed_autosome_ct;
            chr_info.autosome_ct = autosome_ct;
            for (uint32_t xymt_idx = 0; xymt_idx != kChrOffsetCt; ++xymt_idx) {
              chr_info.xymt_codes[xymt_idx] = UINT32_MAXM1;
            }
            SetAllBits(autosome_ct + 1, chr_info.haploid_mask);
          } else {
            const uint32_t autosome_ct = signed_autosome_ct;
            chr_info.autosome_ct = autosome_ct;
            // assumes first four codes are x, y, xy, mt
            // update (11 Aug 2023): unless no-xy is specified, PARs are now
            // enabled.
            for (uint32_t xymt_idx = 0; xymt_idx != kChrOffsetCt; ++xymt_idx) {
              chr_info.xymt_codes[xymt_idx] = autosome_ct + 1 + xymt_idx;
            }
            SetBit(autosome_ct + 1, chr_info.haploid_mask);
            SetBit(autosome_ct + 2, chr_info.haploid_mask);
            SetBit(autosome_ct + 4, chr_info.haploid_mask);
            for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
              cur_modif = argvk[arg_idx + param_idx];
              const uint32_t cur_modif_slen = strlen(cur_modif);
              if (strequal_k(cur_modif, "no-x", cur_modif_slen)) {
                chr_info.xymt_codes[kChrOffsetX] = UINT32_MAXM1;
                ClearBit(autosome_ct + 1, chr_info.haploid_mask);
              } else if (strequal_k(cur_modif, "no-y", cur_modif_slen)) {
                chr_info.xymt_codes[kChrOffsetY] = UINT32_MAXM1;
                ClearBit(autosome_ct + 2, chr_info.haploid_mask);
              } else if (strequal_k(cur_modif, "no-xy", cur_modif_slen)) {
                chr_info.xymt_codes[kChrOffsetXY] = UINT32_MAXM1;
                chr_info.xymt_codes[kChrOffsetPAR1] = UINT32_MAXM1;
                chr_info.xymt_codes[kChrOffsetPAR2] = UINT32_MAXM1;
              } else if (likely(strequal_k(cur_modif, "no-mt", cur_modif_slen))) {
                chr_info.xymt_codes[kChrOffsetMT] = UINT32_MAXM1;
                ClearBit(autosome_ct + 4, chr_info.haploid_mask);
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --chr-set argument '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "hr-override")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (likely(!strcmp(cur_modif, "file"))) {
              pc.misc_flags |= kfMiscChrOverrideFile;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --chr-override argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.misc_flags |= kfMiscChrOverrideCmdline;
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t first_fname_param_idx = 0;
          uint32_t explicit_cols = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              if (unlikely(first_fname_param_idx)) {
                logerrputs("Error: Invalid --clump argument sequence ('zs' must come before filename(s)).\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.clump_info.flags |= kfClumpZs;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (unlikely(first_fname_param_idx)) {
                logerrputs("Error: Invalid --clump argument sequence ('cols=' must come before\nfilename(s)).\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely(explicit_cols)) {
                logerrputs("Error: Multiple --clump cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[strlen("cols=")]), "chrom\0pos\0ref\0alt1\0alt\0maybeprovref\0provref\0maybea1\0a1\0maybef\0f\0total\0maybebounds\0bounds\0bins\0sp2\0", "clump cols=", kfClumpColChrom, kfClumpColDefault, 0, &pc.clump_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              if (!first_fname_param_idx) {
                first_fname_param_idx = param_idx;
              }
            }
          }
          if (unlikely(!first_fname_param_idx)) {
            logerrputs("Error: --clump requires at least one filename.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = AllocAndFlattenCommaDelimEx(&(argvk[arg_idx + first_fname_param_idx]), flagname_p, param_ct + 1 - first_fname_param_idx, 1, &pc.clump_info.fnames_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          if (!explicit_cols) {
            pc.clump_info.flags |= kfClumpColDefault;
          }
          pc.command_flags1 |= kfCommand1Clump;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "lump-id-field")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-id-field must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.clump_info.id_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-p-field")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-p-field must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.clump_info.p_field);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-kb")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-kb must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double dxx;
          if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &dxx)) || (dxx < 0.001))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --clump-kb argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          dxx *= 1000;
          if (dxx > 2147483647.0) {
            pc.clump_info.bp_radius = 0x7ffffffe;
          } else {
            pc.clump_info.bp_radius = S_CAST(int32_t, dxx * (1.0 + kSmallEpsilon) - 1);
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-log10-p1")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-log10-p1 must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double neglog10_p1;
          if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &neglog10_p1)) || (neglog10_p1 < 0.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --clump-log10-p1 argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.clump_info.ln_p1 = neglog10_p1 * (-kLn10 * (1.0 - kSmallEpsilon));
          clump_log10_p1_present = 1;
        } else if (strequal_k_unsafe(flagname_p2, "lump-log10-p2")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-log10-p2 must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double neglog10_p2;
          if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &neglog10_p2)) || (neglog10_p2 < 0.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --clump-log10-p2 argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.clump_info.ln_p2 = neglog10_p2 * (-kLn10 * (1.0 - kSmallEpsilon));
          clump_log10_p2_present = 1;
        } else if (strequal_k_unsafe(flagname_p2, "lump-p1")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-p1 must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(clump_log10_p1_present)) {
            logerrputs("Error: --clump-p1 cannot be used with --clump-log10-p1.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double ln_p1;
          if (unlikely((!ScantokLn(argvk[arg_idx + 1], &ln_p1)) || (ln_p1 > 0.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --clump-p1 argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.clump_info.ln_p1 = ln_p1 * (1.0 - kSmallEpsilon);
        } else if (strequal_k_unsafe(flagname_p2, "lump-p2")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-p2 must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(clump_log10_p2_present)) {
            logerrputs("Error: --clump-p2 cannot be used with --clump-log10-p2.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double ln_p2;
          if (unlikely((!ScantokLn(argvk[arg_idx + 1], &ln_p2)) || (ln_p2 > 0.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --clump-p2 argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.clump_info.ln_p2 = ln_p2 * (1.0 - kSmallEpsilon);
        } else if (strequal_k_unsafe(flagname_p2, "lump-r2")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-r2 must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double dxx;
          if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &dxx)) || (dxx >= 1.0 - kSmallEpsilon))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --clump-r2 argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.clump_info.r2 = dxx * (1.0 + kSmallEpsilon);
        } else if (strequal_k_unsafe(flagname_p2, "lump-a1-field")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-a1-field must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (!param_ct) {
            pc.clump_info.flags |= kfClumpNoA1;
          } else {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.clump_info.a1_field);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-test")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-test must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (!param_ct) {
            pc.clump_info.flags |= kfClumpNoTest;
          } else {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.clump_info.test_name);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-test-field")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-test-field must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (!param_ct) {
            pc.clump_info.flags |= kfClumpNoTest;
          } else {
            if (unlikely(pc.clump_info.flags & kfClumpNoTest)) {
              logerrputs("Error: Nonempty --clump-test-field does not make sense with empty --clump-test\nargument.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.clump_info.test_field);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-range") || strequal_k_unsafe(flagname_p2, "lump-range0")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-range[0] must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(pc.clump_info.range_fname)) {
            logerrputs("Error: --clump-range and --clump-range0 cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.clump_info.range_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          if (flagname_p2[10]) {
            pc.clump_info.flags |= kfClumpRange0;
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-range-border")) {
          if (unlikely(!pc.clump_info.range_fname)) {
            logerrputs("Error: --clump-range-border must be used with --clump-range[0].\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double dxx;
          if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &dxx)) || (dxx < 0.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --clump-range-border argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (dxx > 2147483.646) {
            pc.clump_info.range_border = 0x7ffffffe;
          } else {
            pc.clump_info.range_border = S_CAST(int32_t, dxx * 1000 * (1 + kSmallEpsilon));
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-bins")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-bins must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(!(pc.clump_info.flags & kfClumpColBins))) {
            logerrputs("Error: --clump-bins does not make sense when --clump 'bins' column set has been\nexcluded.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // may as well enforce limit of 2^26 bin-bounds.  this should never
          // come up, and it lets us remove some bounds-checks.
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, kClumpMaxBinBounds))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t comma_ct = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            comma_ct += CountByteInStr(argvk[arg_idx + param_idx], ',');
          }
          const uint32_t bin_bound_ct = comma_ct + param_ct;
          if (bin_bound_ct > kClumpMaxBinBounds) {
            logerrputs("Error: --clump-bins is currently limited to 2^26 bin boundaries.\n");
            reterr = kPglRetNotYetSupported;
            goto main_ret_1;
          }
          pc.clump_info.bin_bound_ct = bin_bound_ct;
          if (unlikely(pgl_malloc(sizeof(double) * bin_bound_ct, &pc.clump_info.ln_bin_boundaries))) {
            goto main_ret_NOMEM;
          }
          double* ln_bin_boundaries_iter = pc.clump_info.ln_bin_boundaries;
          double prev_ln = -DBL_MAX;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* arg_iter = &(argvk[arg_idx + param_idx][0]);
            while (1) {
              double cur_ln;
              arg_iter = ScanadvLn(arg_iter, &cur_ln);
              if (unlikely((!arg_iter) || ((*arg_iter != '\0') && (*arg_iter != ',')))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --clump-bins argument '%s'.\n", argvk[arg_idx + param_idx]);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (unlikely(cur_ln <= prev_ln)) {
                logerrputs("Error: --clump-bins values are not in increasing order.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely(cur_ln >= 0.0)) {
                logerrputs("Error: --clump-bins values >= 1 do not make sense.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              *ln_bin_boundaries_iter++ = cur_ln * (1.0 + kSmallEpsilon);
              if (*arg_iter == '\0') {
                break;
              }
              ++arg_iter;
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-allow-overlap")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-allow-overlap must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.clump_info.flags |= kfClumpAllowOverlap;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "lump-force-a1")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-force-a1 must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(pc.clump_info.flags & kfClumpNoA1)) {
            logerrputs("Error: --clump-force-a1 does not make sense with empty --clump-a1-field\nargument.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.clump_info.flags |= kfClumpForceA1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "lump-log10")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-log10 must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "input-only", cur_modif_slen)) {
              pc.clump_info.flags |= kfClumpInputLog10;
            } else if (likely(strequal_k(cur_modif, "output-only", cur_modif_slen))) {
              pc.clump_info.flags |= kfClumpOutputLog10;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --clump-log10 argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.clump_info.flags |= kfClumpInputLog10 | kfClumpOutputLog10;
          }
        } else if (strequal_k_unsafe(flagname_p2, "lump-unphased")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Clump))) {
            logerrputs("Error: --clump-unphased must be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.clump_info.flags |= kfClumpUnphased;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ovar-quantile-normalize")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.covar_quantnorm_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.pheno_transform_flags |= kfPhenoTransformQuantnormCovar;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "heck-sex")) {
          // parsing is almost identical for --impute-sex
        main_parse_check_sex:
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 7))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t explicit_cols = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (StrStartsWith(cur_modif, "max-female-xf=", cur_modif_slen)) {
              if (pc.check_sex_info.max_female_xf != -1.0) {
                logerrprintf("Error: Multiple --%s max-female-xf= modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE;
              }
              const char* arg_start = &(cur_modif[strlen("max-female-xf=")]);
              double dxx;
              if (unlikely((!ScantokDouble(arg_start, &dxx)) || (dxx < 0.0) || (dxx >= 1.0 - kSmallEpsilon))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s max-female-xf= argument '%s'.\n", flagname_p, arg_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              pc.check_sex_info.flags |= kfCheckSexUseX;
              pc.check_sex_info.max_female_xf = dxx * (1.0 + kSmallEpsilon);
            } else if (StrStartsWith(cur_modif, "min-male-xf=", cur_modif_slen)) {
              if (pc.check_sex_info.min_male_xf != -1.0) {
                logerrprintf("Error: Multiple --%s min-male-xf= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* arg_start = &(cur_modif[strlen("min-male-xf=")]);
              double dxx;
              if (unlikely((!ScantokDouble(arg_start, &dxx)) || (dxx <= 0.0) || (dxx > 1.0))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s min-male-xf= argument '%s'.\n", flagname_p, arg_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              pc.check_sex_info.flags |= kfCheckSexUseX;
              pc.check_sex_info.min_male_xf = dxx * (1.0 - kSmallEpsilon);
            } else if (StrStartsWith(cur_modif, "max-female-ycount=", cur_modif_slen)) {
              if (pc.check_sex_info.max_female_ycount != UINT32_MAX) {
                logerrprintf("Error: Multiple --%s max-female-ycount= modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE;
              }
              pc.check_sex_info.flags |= kfCheckSexUseY;
              const char* arg_start = &(cur_modif[strlen("max-female-ycount=")]);
              if (unlikely(ScanUintDefcapx(arg_start, &pc.check_sex_info.max_female_ycount))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s max-female-ycount= argument '%s'.\n", flagname_p, arg_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (StrStartsWith(cur_modif, "min-male-ycount=", cur_modif_slen)) {
              if (pc.check_sex_info.min_male_ycount != UINT32_MAX) {
                logerrprintf("Error: Multiple --%s min-male-ycount= modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE;
              }
              pc.check_sex_info.flags |= kfCheckSexUseY;
              const char* arg_start = &(cur_modif[strlen("min-male-ycount=")]);
              if (unlikely(ScanPosintDefcapx(arg_start, &pc.check_sex_info.min_male_ycount))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s min-male-ycount= argument '%s'.\n", flagname_p, arg_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (StrStartsWith(cur_modif, "max-female-yrate=", cur_modif_slen)) {
              if (pc.check_sex_info.max_female_yrate != -1.0) {
                logerrprintf("Error: Multiple --%s max-female-yrate= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* arg_start = &(cur_modif[strlen("max-female-yrate=")]);
              double dxx;
              if (unlikely((!ScantokDouble(arg_start, &dxx)) || (dxx < 0.0) || (dxx >= 1.0 - kSmallEpsilon))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s max-female-yrate= argument '%s'.\n", flagname_p, arg_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              pc.check_sex_info.max_female_yrate = dxx * (1.0 + kSmallEpsilon);
              pc.check_sex_info.flags |= kfCheckSexUseY;
            } else if (StrStartsWith(cur_modif, "min-male-yrate=", cur_modif_slen)) {
              if (pc.check_sex_info.min_male_yrate != -1.0) {
                logerrprintf("Error: Multiple --%s min-male-yrate= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* arg_start = &(cur_modif[strlen("min-male-yrate=")]);
              double dxx;
              if (unlikely((!ScantokDouble(arg_start, &dxx)) || (dxx <= 0.0) || (dxx > 1.0))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s min-male-yrate= argument '%s'.\n", flagname_p, arg_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              pc.check_sex_info.min_male_yrate = dxx * (1.0 - kSmallEpsilon);
              pc.check_sex_info.flags |= kfCheckSexUseY;
            } else if (likely(StrStartsWith(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(pc.check_sex_info.flags & kfCheckSexColAll)) {
                logerrprintf("Error: Multiple --%s cols= modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0pedsex\0status\0xf\0ycount\0yrate\0yobs\0", flagname_p, kfCheckSexColMaybefid, kfCheckSexColDefault, 0, &pc.check_sex_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              logerrprintf("Error: Invalid --%s argument sequence.\n", flagname_p);
              double dxx;
              if ((param_idx == 1) && (param_ct == 2) && ScantokDouble(cur_modif, &dxx) && ScantokDouble(argvk[arg_idx + 2], &dxx)) {
                logerrprintfww("(Did you mean \"--%s max-female-xf=%s min-male-xf=%s\"?)\n", flagname_p, cur_modif, argvk[arg_idx + 2]);
              }
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          if (!(pc.check_sex_info.flags & (kfCheckSexUseX | kfCheckSexUseY))) {
            logerrprintfww("Warning: No --%s thresholds specified. Setting min-male-xf=1 and max-female-yrate=0; if you aren't just sanity-checking pre-cleaned data, you should look at the distributions of xf and yrate in the .checksex output file, and then rerun --check-sex with data-derived thresholds.\n", flagname_p);
            pc.check_sex_info.min_male_xf = 1.0 - kSmallEpsilon;
            pc.check_sex_info.max_female_yrate = 0.0;
            pc.check_sex_info.flags |= kfCheckSexUseX | kfCheckSexUseY;
          } else {
            if (unlikely((pc.check_sex_info.max_female_xf != -1.0) && (pc.check_sex_info.min_male_xf != -1.0) && (pc.check_sex_info.min_male_xf <= pc.check_sex_info.max_female_xf))) {
              logerrprintfww("Error: --%s max-female-xf= argument must be smaller than min-male-xf= argument.\n", flagname_p);
              goto main_ret_INVALID_CMDLINE;
            }
            if (pc.check_sex_info.flags & kfCheckSexUseY) {
              if (unlikely((pc.check_sex_info.max_female_ycount != UINT32_MAX) && (pc.check_sex_info.min_male_ycount <= pc.check_sex_info.max_female_ycount))) {
                logerrprintfww("Error: --%s max-female-ycount= argument must be smaller than min-male-ycount= argument.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE;
              }
              if (unlikely((pc.check_sex_info.max_female_yrate != -1.0) && (pc.check_sex_info.min_male_yrate != -1.0) && (pc.check_sex_info.min_male_yrate <= pc.check_sex_info.max_female_yrate))) {
                logerrprintfww("Error: --%s max-female-yrate= argument must be smaller than min-male-yrate= argument.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE;
              }
            }
          }
          if (!explicit_cols) {
            pc.check_sex_info.flags |= kfCheckSexColDefault;
          }
          pc.command_flags1 |= kfCommand1CheckOrImputeSex;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (likely(strequal_k_unsafe(flagname_p2, "ovar-variance-standardize"))) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.vstd_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.pheno_transform_flags |= kfPhenoTransformVstdCovar;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ovar-number")) {
          logerrputs("Error: --covar-number is retired.  Use --covar-col-nums instead (and add 2 to\nconvert from PLINK 1.x covariate-indexes to covariate-column-numbers).\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'd':
        if (strequal_k_unsafe(flagname_p2, "ouble-id")) {
          if (unlikely(const_fid)) {
            logerrputs("Error: --double-id cannot be used with --const-fid.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          import_flags |= kfImportDoubleId;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ebug")) {
          g_debug_on = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ata")) {
          if (unlikely(load_params || (xload & (~kfXloadOxBgen)))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(param_ct == 1)) {
            logerrputs("Error: --data now requires a REF/ALT mode ('ref-first', 'ref-last', or\n'ref-unknown').\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t is_gzs = 0;
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "ref-first", cur_modif_slen)) {
              oxford_import_flags |= kfOxfordImportRefFirst;
            } else if (strequal_k(cur_modif, "ref-unknown", cur_modif_slen)) {
              oxford_import_flags |= kfOxfordImportRefUnknown;
            } else if (strequal_k(cur_modif, "ref-last", cur_modif_slen) ||
                       strequal_k(cur_modif, "ref-second", cur_modif_slen)) {
              oxford_import_flags |= kfOxfordImportRefLast;
            } else if (likely(strequal_k(cur_modif, "gzs", cur_modif_slen))) {
              if (unlikely(xload & kfXloadOxBgen)) {
                // may as well permit e.g. --data ref-first + --bgen
                logerrputs("Error: --data 'gzs' modifier cannot be used with .bgen input.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              is_gzs = 1;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --data argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(oxford_import_flags & kfOxfordImportRefAll)) {
            logerrputs("Error: No --data REF/ALT mode specified ('ref-first', 'ref-last', or\n'ref-unknown'); this is now required.\n");
          }
          const char* fname_prefix = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname_prefix);
          if (unlikely(slen > (kPglFnamesize - 9))) {
            logerrputs("Error: --data argument too long.\n");
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
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dosage_erase_frac;
          if (unlikely((!ScantokDouble(cur_modif, &dosage_erase_frac)) || (dosage_erase_frac < 0.0) || (dosage_erase_frac >= (0.5 - kSmallEpsilon)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dosage-erase-threshold argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.dosage_erase_thresh = S_CAST(int32_t, dosage_erase_frac * ((1 + kSmallEpsilon) * kDosageMid));
        } else if (strequal_k_unsafe(flagname_p2, "ummy")) {
          if (unlikely(load_params || xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 9))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(ScanPosintDefcapx(argvk[arg_idx + 1], &gendummy_info.sample_ct))) {
            logerrputs("Error: Invalid --dummy sample count.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(ScanPosintDefcapx(argvk[arg_idx + 2], &gendummy_info.variant_ct))) {
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
              if (unlikely(ScanUintCappedx(pheno_ct_start, kMaxPhenoCt, &gendummy_info.pheno_ct))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dummy pheno-ct= argument '%s'.\n", pheno_ct_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "scalar-pheno", cur_modif_slen)) {
              gendummy_info.flags |= kfGenDummyScalarPheno;
            } else if (StrStartsWith(cur_modif, "phase-freq=", cur_modif_slen)) {
              const char* phase_freq_start = &(cur_modif[strlen("phase-freq=")]);
              double dxx;
              if (unlikely((!ScantokDouble(phase_freq_start, &dxx)) || (dxx < 0.0) || (dxx > 1.0))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dummy phase-freq= argument '%s'.\n", phase_freq_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              gendummy_info.phase_freq = dxx;
            } else if (StrStartsWith(cur_modif, "dosage-freq=", cur_modif_slen)) {
              const char* dosage_freq_start = &(cur_modif[strlen("dosage-freq=")]);
              double dxx;
              if (unlikely((!ScantokDouble(dosage_freq_start, &dxx)) || (dxx < 0.0) || (dxx > 1.0))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dummy dosage-freq= argument '%s'.\n", dosage_freq_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              gendummy_info.dosage_freq = dxx;
            } else if (!extra_numeric_param_ct) {
              const uint32_t mfreq_ct = 1 + CountByte(cur_modif, ',', cur_modif_slen);
              gendummy_info.geno_mfreq_ct = mfreq_ct;
              if (unlikely(pgl_malloc(mfreq_ct * sizeof(double), &gendummy_info.geno_mfreqs))) {
                goto main_ret_NOMEM;
              }
              const char* cur_modif_iter = cur_modif;
              for (uint32_t mfreq_idx = 0; ; ) {
                double dxx;
                cur_modif_iter = ScanadvDouble(cur_modif_iter, &dxx);
                if (unlikely((!cur_modif_iter) || (dxx < 0.0) || (dxx > 1.0) || ((*cur_modif_iter != ',') && (*cur_modif_iter != '\0')))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dummy argument '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                gendummy_info.geno_mfreqs[mfreq_idx] = dxx;
                if (++mfreq_idx == mfreq_ct) {
                  break;
                }
                ++cur_modif_iter;
              }
              ++extra_numeric_param_ct;
            } else if (likely(extra_numeric_param_ct == 1)) {
              double dxx;
              if (unlikely((!ScantokDouble(cur_modif, &dxx)) || (dxx < 0.0) || (dxx > 1.0))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dummy argument '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              gendummy_info.pheno_mfreq = dxx;
              ++extra_numeric_param_ct;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --dummy argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const uint32_t mutually_exclusive_flags = gendummy_info.flags & (kfGenDummyAcgt | kfGenDummy1234 | kfGenDummy12);
          if (unlikely(mutually_exclusive_flags & (mutually_exclusive_flags - 1))) {
            logerrputs("Error: --dummy 'acgt', '1234', and '12' modifiers are mutually exclusive.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          xload |= kfXloadGenDummy;
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "ummy-coding"))) {
          logerrputs("Error: --dummy-coding is retired.  Use --split-cat-pheno instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "osage"))) {
          logerrputs("Error: --dosage has been replaced with --import-dosage, which converts to .pgen\nformat and provides access to the full range of plink2 flags.  (Run --glm on\nthe imported dataset to invoke the original --dosage linear/logistic\nregression.)\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "og")) {
          if (unlikely(chr_info.chrset_source)) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.chrset_source = kChrsetSourceCmdline;
          chr_info.autosome_ct = 38;
          chr_info.xymt_codes[0] = 39;
          chr_info.xymt_codes[1] = 40;
          chr_info.xymt_codes[2] = 41;
          chr_info.xymt_codes[3] = 42;
          chr_info.xymt_codes[4] = 43;
          chr_info.xymt_codes[5] = 44;
#ifdef __LP64__
          chr_info.haploid_mask[0] = 0x58000000000LLU;
#else
          chr_info.haploid_mask[0] = 0;
          chr_info.haploid_mask[1] = 0x580;
#endif
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "elete-pmerge-result")) {
          delete_pmerge_result = 1;
          pmerge_required = 1;
        } else if (unlikely(*flagname_p2)) {
          // --d is preprocessed
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'e':
        if (strequal_k_unsafe(flagname_p2, "xtract")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          uint32_t is_interval_bed = 0;
          if (strequal_k(cur_modif, "bed0", cur_modif_slen) || strequal_k(cur_modif, "ibed0", cur_modif_slen)) {
            if (unlikely(param_ct == 1)) {
              logerrputs("Error: '--extract bed0' requires at least one filename.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.filter_flags |= kfFilterExtractBed0 | kfFilterNoSplitChr;
            is_interval_bed = 1;
          } else if (strequal_k(cur_modif, "bed1", cur_modif_slen) || strequal_k(cur_modif, "ibed1", cur_modif_slen) || strequal_k(cur_modif, "range", cur_modif_slen)) {
            if (unlikely(param_ct == 1)) {
              logerrputs("Error: '--extract bed1' requires at least one filename.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.filter_flags |= kfFilterExtractBed1 | kfFilterNoSplitChr;
            is_interval_bed = 1;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1 + is_interval_bed]), flagname_p, param_ct - is_interval_bed, &pc.extract_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xtract-intersect")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          uint32_t is_interval_bed = 0;
          if (strequal_k(cur_modif, "bed0", cur_modif_slen)) {
            if (unlikely(param_ct == 1)) {
              logerrputs("Error: '--extract-intersect bed0' requires at least one filename.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.filter_flags |= kfFilterExtractIntersectBed0 | kfFilterNoSplitChr;
            is_interval_bed = 1;
          } else if (strequal_k(cur_modif, "bed1", cur_modif_slen)) {
            if (unlikely(param_ct == 1)) {
              logerrputs("Error: '--extract-intersect bed1' requires at least one filename.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.filter_flags |= kfFilterExtractIntersectBed1 | kfFilterNoSplitChr;
            is_interval_bed = 1;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1 + is_interval_bed]), flagname_p, param_ct - is_interval_bed, &pc.extract_intersect_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xclude")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          uint32_t is_interval_bed = 0;
          if (strequal_k(cur_modif, "bed0", cur_modif_slen) || strequal_k(cur_modif, "ibed0", cur_modif_slen)) {
            if (unlikely(param_ct == 1)) {
              logerrputs("Error: '--exclude bed0' requires at least one filename.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.filter_flags |= kfFilterExcludeBed0 | kfFilterNoSplitChr;
            is_interval_bed = 1;
          } else if (strequal_k(cur_modif, "bed1", cur_modif_slen) || strequal_k(cur_modif, "ibed1", cur_modif_slen) || strequal_k(cur_modif, "range", cur_modif_slen)) {
            if (unlikely(param_ct == 1)) {
              logerrputs("Error: '--exclude bed1' requires at least one filename.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.filter_flags |= kfFilterExcludeBed1 | kfFilterNoSplitChr;
            is_interval_bed = 1;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1 + is_interval_bed]), flagname_p, param_ct - is_interval_bed, &pc.exclude_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xport")) {
          // todo: determine actual limit
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 50))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint64_t format_param_idxs = 0;
          GetExportfTargets(&(argvk[arg_idx]), param_ct, &pc.exportf_info.flags, &pc.exportf_info.idpaste_flags, &format_param_idxs);
          if (unlikely(!format_param_idxs)) {
            logerrputs("Error: --export requires at least one output format.  (Did you forget 'ped' or\n'vcf'?)\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // can't have e.g. bgen-1.1 and bgen-1.2 simultaneously, since they
          // have the same extension and different content.
          const uint64_t bgen_flags = S_CAST(uint64_t, pc.exportf_info.flags & (kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13));
          if (unlikely(bgen_flags & (bgen_flags - 1))) {
            logerrputs("Error: Multiple --export bgen versions.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely((pc.exportf_info.flags & kfExportfVcf) == kfExportfVcf)) {
            logerrputs("Error: 'vcf' and 'vcf-4.2' formats cannot be exported simultaneously.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.exportf_info.flags & kfExportfBcf) == kfExportfBcf)) {
            logerrputs("Error: 'bcf' and 'bcf-4.2' formats cannot be exported simultaneously.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.exportf_info.flags & (kfExportfHaps | kfExportfHapsLegend)) == (kfExportfHaps | kfExportfHapsLegend))) {
            logerrputs("Error: 'haps' and 'hapslegend' formats cannot be exported simultaneously.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely((pc.exportf_info.flags & (kfExportfA | kfExportfAD)) == (kfExportfA | kfExportfAD))) {
            logerrputs("Error: 'A' and 'AD' formats cannot be exported simultaneously.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely((pc.exportf_info.flags & (kfExportfOxGenV1 | kfExportfOxGenV2)) == (kfExportfOxGenV1 | kfExportfOxGenV2))) {
            logerrputs("Error: 'oxford' and 'oxford-v2' formats cannot be exported simultaneously.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely((pc.exportf_info.flags & (kfExportfPed | kfExportfCompound)) == (kfExportfPed | kfExportfCompound))) {
            logerrputs("Error: 'ped' and 'compound-genotypes' formats cannot be exported\nsimultaneously.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely((pc.exportf_info.flags & (kfExportfPhylip | kfExportfPhylipPhased)) == (kfExportfPhylip | kfExportfPhylipPhased))) {
            logerrputs("Error: --export 'phylip' and 'phylip-phased' formats cannot be exported\nsimultaneously.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.exportf_info.flags & (kfExportfEig | kfExportfEigt)) == (kfExportfEig | kfExportfEigt))) {
            logerrputs("Error: --export 'eig' and 'eigt' formats cannot be exported simultaneously.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            // could use AdvBoundedTo0Bit()...
            if ((format_param_idxs >> param_idx) & 1) {
              continue;
            }
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (StrStartsWith(cur_modif, "id-paste=", cur_modif_slen)) {
              if (unlikely(pc.exportf_info.idpaste_flags)) {
                logerrputs("Error: Multiple --export id-paste= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("id-paste=")]), "maybefid\0fid\0iid\0maybesid\0sid\0", "export", kfIdpasteMaybefid, kfIdpasteDefault, 1, &pc.exportf_info.idpaste_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "id-delim=", cur_modif_slen)) {
              if (unlikely(pc.exportf_info.id_delim)) {
                logerrputs("Error: Multiple --export id-delim= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              pc.exportf_info.id_delim = ExtractCharParam(&(cur_modif[strlen("id-delim=")]));
              if (unlikely(!pc.exportf_info.id_delim)) {
                logerrputs("Error: --export id-delim= value must be a single character.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely((ctou32(pc.exportf_info.id_delim) < ' ') || (pc.exportf_info.id_delim == '0'))) {
                logerrputs("Error: --export id-delim= value cannot be tab, newline, '0', or a nonprinting\ncharacter.\n");
                goto main_ret_INVALID_CMDLINE;
              }
            } else if (StrStartsWith(cur_modif, "vcf-dosage=", cur_modif_slen)) {
              if (unlikely(!(pc.exportf_info.flags & (kfExportfVcf | kfExportfBcf)))) {
                logerrputs("Error: The 'vcf-dosage' modifier only applies to --export's vcf and bcf output\nformats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely(pc.exportf_info.vcf_mode != kVcfExport0)) {
                logerrputs("Error: Multiple --export vcf-dosage= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* vcf_dosage_start = &(cur_modif[strlen("vcf-dosage=")]);
              const uint32_t vcf_dosage_start_slen = strlen(vcf_dosage_start);
              if (strequal_k(vcf_dosage_start, "GP", vcf_dosage_start_slen)) {
                if (pc.exportf_info.flags & (kfExportfVcf42 | kfExportfBcf42)) {
                  logerrputs("Error: --export vcf-dosage=GP cannot be used in {v,b}cf-4.2 mode.\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
                pc.exportf_info.vcf_mode = kVcfExportGp;
              } else if (strequal_k(vcf_dosage_start, "DS", vcf_dosage_start_slen)) {
                pc.exportf_info.vcf_mode = kVcfExportDs;
              } else if (strequal_k(vcf_dosage_start, "DS-force", vcf_dosage_start_slen)) {
                pc.exportf_info.vcf_mode = kVcfExportDsForce;
              } else if (strequal_k(vcf_dosage_start, "DS-only", vcf_dosage_start_slen)) {
                pc.exportf_info.vcf_mode = kVcfExportDsOnly;
              } else if (strequal_k(vcf_dosage_start, "HDS", vcf_dosage_start_slen)) {
                pc.exportf_info.vcf_mode = kVcfExportHds;
              } else if (likely(strequal_k(vcf_dosage_start, "HDS-force", vcf_dosage_start_slen))) {
                pc.exportf_info.vcf_mode = kVcfExportHdsForce;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --export vcf-dosage= argument '%s'.\n", vcf_dosage_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (StrStartsWith(cur_modif, "bits=", cur_modif_slen)) {
              if (unlikely(!(pc.exportf_info.flags & (kfExportfBgen12 | kfExportfBgen13)))) {
                logerrputs("Error: The 'bits' modifier only applies to --export's bgen-1.2 and bgen-1.3\noutput formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely(pc.exportf_info.bgen_bits)) {
                logerrputs("Error: Multiple --export bits= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* bits_start = &(cur_modif[strlen("bits=")]);
              uint32_t bgen_bits;
              if (unlikely(ScanPosintCappedx(bits_start, 24, &bgen_bits))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --export bits= argument '%s'.\n", bits_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (bgen_bits > 16) {
                // okay, bits=24 should eventually be permitted because of
                // multiallelic-variant dosage precision issues
                logerrputs("Error: Current --export bits= maximum is 16.\n");
                reterr = kPglRetNotYetSupported;
                goto main_ret_1;
              }
              if (bgen_bits & (bgen_bits - 1)) {
                logerrputs("Warning: Support for most non-power-of-2 bits= export values is likely to be\ndiscontinued, since .bgen size tends to be larger than the\nnext-higher-power-of-2 precision level.\n");
              }
              pc.exportf_info.bgen_bits = bgen_bits;
            } else if (strequal_k(cur_modif, "include-alt", cur_modif_slen)) {
              if (unlikely(!(pc.exportf_info.flags & (kfExportfA | kfExportfAD)))) {
                logerrputs("Error: The 'include-alt' modifier only applies to --export's A and AD output\nformats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.exportf_info.flags |= kfExportfIncludeAlt;
            } else if (strequal_k(cur_modif, "omit-nonmale-y", cur_modif_slen)) {
              if (unlikely(!(pc.exportf_info.flags & (kfExportfList | kfExportfRlist)))) {
                logerrputs("Error: The 'omit-nonmale-y' modifier only applies to --export's list and rlist\noutput formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.exportf_info.flags |= kfExportfOmitNonmaleY;
            } else if (strequal_k(cur_modif, "01", cur_modif_slen) || strequal_k(cur_modif, "12", cur_modif_slen)) {
              if (unlikely(pc.exportf_info.flags & (kfExportfA | kfExportfAD))) {
                snprintf(g_logbuf, kLogbufSize, "Error: The '%s' modifier does not apply to --export's A and AD output formats.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              if (unlikely(pc.exportf_info.flags & (kfExportfVcf | kfExportfBcf | kfExportfPhylip | kfExportfPhylipPhased))) {
                logerrputs("Error: '01'/'12' cannot be used with --export's vcf, bcf, phylip, or\nphylip-phased output formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (cur_modif[0] == '0') {
                if (unlikely(pc.exportf_info.flags & kfExportf12)) {
                  logerrputs("Error: --export '01' and '12' cannot be used together.\n");
                  goto main_ret_INVALID_CMDLINE;
                }
                pc.exportf_info.flags |= kfExportf01;
              } else {
                if (unlikely(pc.exportf_info.flags & kfExportf01)) {
                  logerrputs("Error: --export '01' and '12' cannot be used together.\n");
                  goto main_ret_INVALID_CMDLINE;
                }
                pc.exportf_info.flags |= kfExportf12;
              }
            } else if (strequal_k(cur_modif, "bgz", cur_modif_slen)) {
              if (unlikely(!(pc.exportf_info.flags & (kfExportfHaps | kfExportfHapsLegend | kfExportfOxGen | kfExportfVcf)))) {
                logerrputs("Error: The 'bgz' modifier only applies to --export's haps[legend], oxford, and\nvcf output formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.exportf_info.flags |= kfExportfBgz;
            } else if (strequal_k(cur_modif, "spaces", cur_modif_slen)) {
              pc.exportf_info.flags |= kfExportfSpaces;
            } else if (strequal_k(cur_modif, "sample-v2", cur_modif_slen)) {
              if (unlikely(!(pc.exportf_info.flags & (kfExportfHaps | kfExportfHapsLegend | kfExportfOxGen | kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13)))) {
                logerrputs("Error: 'sample-v2' modifier only applies to --export's haps[legend], oxford,\nand bgen output formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.exportf_info.flags |= kfExportfSampleV2;
            } else if (strequal_k(cur_modif, "bgen-omit-sample-id-block", cur_modif_slen)) {
              if (unlikely(!(pc.exportf_info.flags & (kfExportfBgen12 | kfExportfBgen13)))) {
                logerrputs("Error: 'bgen-omit-sample-id-block' modifier only applies to --export's bgen-1.2\nand bgen-1.3 output formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.exportf_info.flags |= kfExportfBgenOmitSampleIdBlock;
            } else if (strequal_k(cur_modif, "used-sites", cur_modif_slen)) {
              if (unlikely(!(pc.exportf_info.flags & (kfExportfPhylip | kfExportfPhylipPhased)))) {
                logerrputs("Error: 'used-sites' modifier only applies to --export's phylip and\nphylip-phased output formats.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.exportf_info.flags |= kfExportfPhylipUsedSites;
            } else if (strequal_k(cur_modif, "ref-last", cur_modif_slen)) {
              // do nothing for now
            } else if (likely(strequal_k(cur_modif, "ref-first", cur_modif_slen))) {
              pc.exportf_info.flags |= kfExportfRefFirst;
            } else if (strequal_k(cur_modif, "gen-gz", cur_modif_slen)) {
              logerrputs("Error: 'gen-gz' modifier retired.  Use '--export oxford bgz' instead.\n");
              goto main_ret_INVALID_CMDLINE_WWA;
            } else if (StrStartsWith(cur_modif, "dosage=", cur_modif_slen)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --export argument '%s'. (Did you mean 'vcf-%s'?)\n", cur_modif, cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --export argument '%s'.%s\n", cur_modif, ((param_idx == param_ct) && (!outname_end))? " (Did you forget '--out'?)" : "");
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (pc.exportf_info.idpaste_flags || pc.exportf_info.id_delim) {
            if (unlikely(!(pc.exportf_info.flags & (kfExportfVcf | kfExportfBcf | kfExportfBgen12 | kfExportfBgen13 | kfExportfSampleV2 | kfExportfPhylip | kfExportfPhylipPhased | kfExportfEig | kfExportfEigt)))) {
              logerrputs("Error: The 'id-delim' and 'id-paste' modifiers only apply to --export's vcf,\nbcf, bgen-1.2, bgen-1.3, sample-v2, phylip, phylip-phased, eig, and eigt output\nformats.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          if (pc.exportf_info.flags & (kfExportfVcf | kfExportfBcf | kfExportfBgen12 | kfExportfBgen13 | kfExportfSampleV2 | kfExportfPhylip | kfExportfPhylipPhased | kfExportfEig | kfExportfEigt)) {
            if (!pc.exportf_info.idpaste_flags) {
              pc.exportf_info.idpaste_flags = kfIdpasteDefault;
            }
          }
          pc.command_flags1 |= kfCommand1Exportf;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "xport-allele")) {
          if (unlikely((!(pc.command_flags1 & kfCommand1Exportf)) || (!(pc.exportf_info.flags & (kfExportfA | kfExportfAv | kfExportfAD))))) {
            logerrputs("Error: --export-allele must be used with --export A/AD/Av.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.exportf_info.export_allele_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "xclude-snp")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_exclude_snp);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xclude-snps")) {
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.exclude_snps_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xtract-col-cond")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = Alloc2col(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.extract_col_cond_info.params);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xtract-col-cond-match")) {
          if (unlikely(!pc.extract_col_cond_info.params)) {
            logerrputs("Error: --extract-col-cond-match must be used with --extract-col-cond.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, kMaxIdBlen, &pc.extract_col_cond_info.match_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "xtract-col-cond-mismatch")) {
          if (unlikely(!pc.extract_col_cond_info.params)) {
            logerrputs("Error: --extract-col-cond-mismatch must be used with --extract-col-cond.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // could make this check airtight?  right now, there's no error if
          // the user also specifies "--extract-col-cond-min 0"
          if (unlikely((pc.extract_col_cond_info.min != 0.0) || (pc.extract_col_cond_info.max != DBL_MAX))) {
            logerrputs("Error: --extract-col-cond-mismatch cannot be used with --extract-col-cond-max\nor --extract-col-cond-min.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, kMaxIdBlen, &pc.extract_col_cond_info.mismatch_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "xtract-col-cond-substr")) {
          if (unlikely((!pc.extract_col_cond_info.match_flattened) && (!pc.extract_col_cond_info.mismatch_flattened))) {
            logerrputs("Error: --extract-col-cond-substr must be used with --extract-col-cond-match\nand/or --extract-col-cond-mismatch.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.extract_col_cond_info.match_substr = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "xtract-col-cond-max")) {
          if (unlikely(!pc.extract_col_cond_info.params)) {
            logerrputs("Error: --extract-col-cond-match must be used with --extract-col-cond.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.extract_col_cond_info.match_flattened)) {
            logerrputs("Error: --extract-col-cond-max cannot be used with --extract-col-cond-match.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(!ScantokDouble(cur_modif, &pc.extract_col_cond_info.max))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --extract-col-cond-max argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "xtract-col-cond-min")) {
          if (unlikely(!pc.extract_col_cond_info.params)) {
            logerrputs("Error: --extract-col-cond-match must be used with --extract-col-cond.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.extract_col_cond_info.match_flattened)) {
            logerrputs("Error: --extract-col-cond-min cannot be used with --extract-col-cond-match.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(!ScantokDouble(cur_modif, &pc.extract_col_cond_info.min))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --extract-col-cond-min argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "xtract-if-info")) {
          reterr = ValidateAndAllocCmpExpr(&(argvk[arg_idx + 1]), argvk[arg_idx], param_ct, &pc.extract_if_info_expr);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          // TODO: prohibit ';' and '=' in keys
          // LoadPvar() currently checks value strings when nonnumeric
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogExtractIfInfo;
        } else if (strequal_k_unsafe(flagname_p2, "xclude-if-info")) {
          reterr = ValidateAndAllocCmpExpr(&(argvk[arg_idx + 1]), argvk[arg_idx], param_ct, &pc.exclude_if_info_expr);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogExcludeIfInfo;
        } else if (strequal_k_unsafe(flagname_p2, "igfile")) {
          if (unlikely(load_params || xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "nohash", &fname_modif_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            import_flags |= kfImportEigNohash;
          }
          const char* fname_prefix = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname_prefix);
          if (unlikely(slen > (kPglFnamesize - 6))) {
            logerrputs("Error: --eigfile argument too long.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          snprintf(memcpya(pgenname, fname_prefix, slen), 6, ".geno");
          snprintf(memcpya(psamname, fname_prefix, slen), 6, ".ind");
          snprintf(memcpya(pvarname, fname_prefix, slen), 6, ".snp");
          xload = kfXloadEigGeno | kfXloadEigInd | kfXloadEigSnp;
        } else if (strequal_k_unsafe(flagname_p2, "iggeno")) {
          if (unlikely(load_params || (xload & (~(kfXloadEigGeno | kfXloadEigInd | kfXloadEigSnp))))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "nohash", &fname_modif_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            import_flags |= kfImportEigNohash;
          }
          const char* fname = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --eiggeno argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, fname, slen + 1);
          xload |= kfXloadEigGeno;
        } else if (strequal_k_unsafe(flagname_p2, "igind")) {
          if (unlikely(load_params || (xload & (~(kfXloadEigGeno | kfXloadEigInd | kfXloadEigSnp))))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(!(xload & kfXloadEigGeno))) {
            logerrputs("Error: --eigind must be used with --eigfile/--eiggeno.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --eigind argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(psamname, fname, slen + 1);
          xload |= kfXloadEigInd;
        } else if (strequal_k_unsafe(flagname_p2, "igind")) {
          if (unlikely(load_params || (xload & (~(kfXloadEigGeno | kfXloadEigInd | kfXloadEigSnp))))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely((xload & (kfXloadEigGeno | kfXloadEigInd)) != (kfXloadEigGeno | kfXloadEigInd))) {
            logerrputs("Error: --eigsnp must be used with either --eigfile, or --eiggeno + --eigind.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --eigsnp argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pvarname, fname, slen + 1);
          xload |= kfXloadEigSnp;
        } else if (strequal_k_unsafe(flagname_p2, "xclude-palindromic-snps")) {
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogExcludePalindromicSnps;
          goto main_param_zero;
        } else if (likely(strequal_k_unsafe(flagname_p2, "rror-on-freq-calc"))) {
          pc.misc_flags |= kfMiscErrorOnFreqCalc;
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'f':
        if (strequal_k_unsafe(flagname_p2, "req")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 5))) {
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
              if (unlikely(pc.freq_rpt_flags & kfAlleleFreqColAll)) {
                logerrputs("Error: Multiple --freq cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0maybeprovref\0provref\0reffreq\0alt1freq\0altfreq\0freq\0eq\0eqz\0alteq\0alteqz\0numeq\0altnumeq\0machr2\0minimac3r2\0nobs\0", "freq", kfAlleleFreqColChrom, kfAlleleFreqColDefault, 1, &pc.freq_rpt_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              const uint32_t mutually_exclusive_cols = pc.freq_rpt_flags & kfAlleleFreqColMutex;
              if (unlikely(mutually_exclusive_cols & (mutually_exclusive_cols - 1))) {
                logerrputs("Error: --freq's altfreq, freq, eq, eqz, alteq, alteqz, numeq, and altnumeq\ncolumns are mutually exclusive.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if ((pc.freq_rpt_flags & (kfAlleleFreqColMachR2 | kfAlleleFreqColMinimac3R2)) == (kfAlleleFreqColMachR2 | kfAlleleFreqColMinimac3R2)) {
                logerrputs("Error: --freq machr2 and minimac3r2 columns are mutually exclusive.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else if (strequal_k(cur_modif, "bins-only", cur_modif_slen)) {
              bins_only = 1;
            } else if (likely(StrStartsWith(cur_modif, "refbins=", cur_modif_slen) ||
                              StrStartsWith(cur_modif, "refbins-file=", cur_modif_slen) ||
                              StrStartsWith(cur_modif, "alt1bins=", cur_modif_slen) ||
                              StrStartsWith(cur_modif, "alt1bins-file=", cur_modif_slen))) {
              const uint32_t is_alt1 = (cur_modif[0] == 'a');
              char** binstr_ptr = is_alt1? (&pc.freq_alt1_binstr) : (&pc.freq_ref_binstr);
              if (unlikely(*binstr_ptr)) {
                logerrprintf("Error: Multiple --freq %sbins[-file]= modifiers.\n", is_alt1? "alt1" : "ref");
                goto main_ret_INVALID_CMDLINE;
              }
              if (cur_modif[7 + is_alt1] == '=') {
                // "refbins=", "alt1bins="
                reterr = CmdlineAllocString(&(cur_modif[8 + is_alt1]), is_alt1? "--freq alt1bins=" : "--freq refbins=", 0x7fffffff, binstr_ptr);
              } else {
                // "refbins-file=", "alt1bins-file="
                pc.freq_rpt_flags |= is_alt1? kfAlleleFreqBinsAlt1Fname : kfAlleleFreqBinsRefFname;
                reterr = AllocFname(&(cur_modif[13 + is_alt1]), is_alt1? "freq alt1bins-file=" : "freq refbins-file=", binstr_ptr);
              }
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --freq argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (bins_only) {
            if (unlikely((!pc.freq_ref_binstr) && (!pc.freq_alt1_binstr))) {
              logerrputs("Error: --freq 'bins-only' must be used with 'refbins[-file]=' and/or\n'alt1bins[-file]='.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(pc.freq_rpt_flags & (kfAlleleFreqZs | kfAlleleFreqColAll))) {
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
          if (unlikely(chr_info.is_include_stack)) {
            logerrputs("Error: --from/--to cannot be used with --autosome[-par] or --[not-]chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_from);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "rom-bp") ||
                   strequal_k_unsafe(flagname_p2, "rom-kb") ||
                   strequal_k_unsafe(flagname_p2, "rom-mb")) {
          if (unlikely(!CmdlineSingleChr(&chr_info, pc.load_filter_log_flags))) {
            logerrputs("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb must be used with --chr, and only\none chromosome.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(pc.from_bp != -1)) {
            logerrputs("Error: Multiple --from-bp/-kb/-mb values.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          // permit negative numbers, to simplify shell script windowing logic
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (unlikely(!ScantokDouble(cur_modif, &dxx))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --from-bp/-kb/-mb argument '%s'.\n", cur_modif);
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
            if (unlikely(dxx > 2147483646.0)) {
              logerrprintf("Error: --from-bp/-kb/-mb argument '%s' too large.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.from_bp = 1 + S_CAST(int32_t, dxx * (1 - kSmallEpsilon));
          }
          pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "orce-intersect")) {
          permit_multiple_inclusion_filters = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "amily")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (unlikely(IsReservedPhenoName(cur_modif, cur_modif_slen))) {
              snprintf(g_logbuf, kLogbufSize, "Error: '%s' cannot be used as a categorical phenotype name.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_2A;
            }
            reterr = CmdlineAllocString(cur_modif, argvk[arg_idx], kMaxIdSlen, &pc.catpheno_name);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.misc_flags |= kfMiscCatPhenoFamily;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "a")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.fa_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "st")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.fst_info.pheno_name);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          // Similar to --sample-diff.
          uint32_t explicit_method = 0;
          uint32_t explicit_cols = 0;
          uint32_t explicit_vcols = 0;
          uint32_t param_idx = 2;
          for (; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (StrStartsWith(cur_modif, "method=", cur_modif_slen)) {
              if (unlikely(explicit_method)) {
                logerrputs("Error: Multiple --fst method= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* method_name = &(cur_modif[7]);
              const uint32_t method_name_slen = cur_modif_slen - 7;
              if (strequal_k(method_name, "wc", method_name_slen)) {
                explicit_method = 2;
              } else if (likely(strequal_k(method_name, "hudson", method_name_slen))) {
                explicit_method = 1;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Unsupported --fst method '%s'.\n", method_name);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (StrStartsWith(cur_modif, "blocksize=", cur_modif_slen)) {
              if (unlikely(pc.fst_info.blocksize)) {
                logerrputs("Error: Multiple --fst blocksize= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              if (unlikely(ScanPosintDefcapx(&(cur_modif[10]), &pc.fst_info.blocksize))) {
                logerrputs("Error: Invalid --fst blocksize.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else if (StrStartsWith0(cur_modif, "cols=", cur_modif_slen)) {
              if (unlikely(explicit_cols)) {
                logerrputs("Error: Multiple --fst cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[strlen("cols=")]), "nobs\0", "fst cols", kfFstColNobs, kfFstColDefault, 0, &pc.fst_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (strequal_k(cur_modif, "report-variants", cur_modif_slen)) {
              pc.fst_info.flags |= kfFstReportVariants;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.fst_info.flags |= kfFstZs;
            } else if (StrStartsWith0(cur_modif, "vcols=", cur_modif_slen)) {
              if (unlikely(explicit_vcols)) {
                logerrputs("Error: Multiple --fst vcols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_vcols = 1;
              reterr = ParseColDescriptor(&(cur_modif[strlen("vcols=")]), "chrom\0pos\0ref\0alt\0maybeprovref\0provref\0nobs\0nallele\0fstfrac\0fst\0", "fst vcols", kfFstVcolChrom, kfFstVcolDefault, 0, &pc.fst_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "base=", cur_modif_slen)) {
              // a.<basepop>.2.fst.var
              reterr = CmdlineAllocString(&(cur_modif[5]), argvk[arg_idx], kPglFnamesize - 13, &(pc.fst_info.first_id_or_fname));
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              pc.fst_info.flags |= kfFstOneBasePop;
              break;
            } else if (StrStartsWith(cur_modif, "ids=", cur_modif_slen)) {
              if (unlikely(param_idx == param_ct)) {
                logerrputs("Error: --fst 'ids=' requires at least two (space-delimited) populations.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              reterr = CmdlineAllocString(&(cur_modif[4]), argvk[arg_idx], kPglFnamesize - 13, &(pc.fst_info.first_id_or_fname));
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              pc.fst_info.flags |= kfFstExplicitPopIds;
              break;
            } else if (likely(StrStartsWith(cur_modif, "file=", cur_modif_slen))) {
              if (param_idx != param_ct) {
                logerrputs("Error: --fst 'file=' must appear after all other modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              reterr = AllocFname(&(cur_modif[5]), argvk[arg_idx], &pc.fst_info.first_id_or_fname);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              pc.fst_info.flags |= kfFstPopPairFile;
              break;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --fst argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (explicit_method == 2) {
            pc.fst_info.flags |= kfFstMethodWc;
          }
          if (!explicit_cols) {
            pc.fst_info.flags |= kfFstColDefault;
          }
          if (!(pc.fst_info.flags & kfFstReportVariants)) {
            if (unlikely(pc.fst_info.flags & kfFstZs)) {
              logerrputs("Error: --fst 'zs' modifier only makes sense with 'report-variants'.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(explicit_vcols)) {
              logerrputs("Error: --fst 'vcols=' modifier only makes sense with 'report-variants'.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else {
            if (!explicit_vcols) {
              pc.fst_info.flags |= kfFstVcolDefault;
            }
          }
          if (param_idx < param_ct) {
            // base= or id=
            const uint32_t other_id_ct = param_ct - param_idx;
            reterr = AllocAndFlatten(&(argvk[arg_idx + param_idx + 1]), flagname_p, other_id_ct, kPglFnamesize - 13, &pc.fst_info.other_ids_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.command_flags1 |= kfCommand1Fst;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ill-missing-with-ref")) {
          make_plink2_flags |= kfMakePlink2FillMissingWithRef;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "lip")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "permissive", &fname_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.flip_info.flags |= kfFlipPermissive;
          }
          reterr = AllocFname(argvk[arg_idx + fname_idx], flagname_p, &pc.flip_info.fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "lip-subset")) {
          if (unlikely(!pc.flip_info.fname)) {
            logerrputs("Error: --flip-subset must be used with --flip.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(make_plink2_flags & kfMakePlink2FillMissingWithRef)) {
            logerrputs("Error: --flip-subset cannot be used with --fill-missing-with-ref.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.allele_alphanum_flags != kfAlleleAlphanum0)) {
            logerrputs("Error: --flip-subset cannot be used with --allele1234/--alleleACGT.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.flip_info.subset_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterAllReq;
        } else if (likely(strequal_k_unsafe(flagname_p2, "amily-missing-catname"))) {
          if (unlikely(!(pc.misc_flags & kfMiscCatPhenoFamily))) {
            logerrputs("Error: --family-missing-catname must be used with --family.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.family_missing_catname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ilter-cases") || strequal_k_unsafe(flagname_p2, "ilter-controls")) {
          logerrputs("Error: --filter-cases and --filter-controls have been retired.  Use\n--keep-if/--remove-if instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "rqx") || strequal_k_unsafe(flagname_p2, "reqx")) {
          logerrputs("Error: --freqx has been retired.  Use --geno-counts instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "ile")) {
          logerrputs("Error: --file has been retired.\nContinued use of the .ped + .map fileset format is strongly discouraged.  Use\n--pedmap to convert to a more appropriate format.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'g':
        if (strequal_k_unsafe(flagname_p2, "eno")) {
          if (unlikely(make_plink2_flags & kfMakePlink2FillMissingWithRef)) {
            // While --geno/--mind technically makes more sense before
            // --fill-missing-with-ref than after it, it doesn't really make
            // sense to use them together at all.
            logerrputs("Error: --geno/--mind cannot be used with --fill-missing-with-ref.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t geno_thresh_present = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "dosage", cur_modif_slen)) {
              pc.misc_flags |= kfMiscGenoDosage;
            } else if (strequal_k(cur_modif, "hh-missing", cur_modif_slen)) {
              pc.misc_flags |= kfMiscGenoHhMissing;
            } else if (unlikely(geno_thresh_present)) {
              logerrputs("Error: Invalid --geno argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (unlikely(!ScantokDouble(cur_modif, &pc.geno_thresh))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --geno argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else if (unlikely((pc.geno_thresh < 0.0) || (pc.geno_thresh > 1.0))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --geno argument '%s' (must be in [0, 1]).\n", cur_modif);
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
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.geno_counts_flags |= kfGenoCountsZs;
            } else if (likely(StrStartsWith(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(pc.geno_counts_flags & kfGenoCountsColAll)) {
                logerrputs("Error: Multiple --geno-counts cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0maybeprovref\0provref\0homref\0refalt1\0refalt\0homalt1\0altxy\0xy\0hapref\0hapalt1\0hapalt\0hap\0numeq\0missing\0nobs\0", "geno-counts", kfGenoCountsColChrom, kfGenoCountsColDefault, 1, &pc.geno_counts_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              if (unlikely((pc.geno_counts_flags & kfGenoCountsColPairex) == kfGenoCountsColPairex)) {
                logerrputs("Error: --geno-counts's hapaltx and hapx columns are mutually exclusive.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              const uint32_t mutually_exclusive_cols = pc.geno_counts_flags & kfGenoCountsColMutex;
              if (unlikely(mutually_exclusive_cols & (mutually_exclusive_cols - 1))) {
                logerrputs("Error: --geno-counts's altxy, xy, and numeq columns are mutually exclusive.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --geno-counts argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.geno_counts_flags & kfGenoCountsColAll)) {
            pc.geno_counts_flags |= kfGenoCountsColDefault;
          }
          pc.command_flags1 |= kfCommand1GenoCounts;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "lm")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 18))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t explicit_firth_fallback = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmZs;
            } else if (strequal_k(cur_modif, "omit-ref", cur_modif_slen) ||
                       strequal_k(cur_modif, "a0-ref", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmOmitRef;
            } else if (strequal_k(cur_modif, "sex", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmSex;
            } else if (strequal_k(cur_modif, "no-x-sex", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmNoXSex;
            } else if (strequal_k(cur_modif, "log10", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmLog10;
            } else if (strequal_k(cur_modif, "pheno-ids", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmPhenoIds;
            } else if (strequal_k(cur_modif, "genotypic", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmGenotypic;
            } else if (strequal_k(cur_modif, "hethom", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmHethom;
            } else if (strequal_k(cur_modif, "dominant", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmDominant;
            } else if (strequal_k(cur_modif, "recessive", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmRecessive;
            } else if (strequal_k(cur_modif, "hetonly", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmHetonly;
            } else if (strequal_k(cur_modif, "interaction", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmInteraction;
            } else if (strequal_k(cur_modif, "hide-covar", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmHideCovar;
            } else if (strequal_k(cur_modif, "intercept", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmIntercept;
            } else if (strequal_k(cur_modif, "skip-invalid-pheno", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmSkipInvalidPheno;
            } else if (strequal_k(cur_modif, "skip", cur_modif_slen)) {
              logputs("Note: --glm 'skip' modifier has been renamed to 'skip-invalid-pheno'.\n");
              pc.glm_info.flags |= kfGlmSkipInvalidPheno;
            } else if (strequal_k(cur_modif, "no-firth", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmNoFirth;
            } else if (strequal_k(cur_modif, "firth-fallback", cur_modif_slen)) {
              explicit_firth_fallback = 1;
            } else if (strequal_k(cur_modif, "firth", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmFirth;
            } else if (strequal_k(cur_modif, "firth-residualize", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmFirthResidualize;
            } else if (strequal_k(cur_modif, "cc-residualize", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmCcResidualize;
            } else if (strequal_k(cur_modif, "qt-residualize", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmQtResidualize;
            } else if (strequal_k(cur_modif, "single-prec-cc", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmSinglePrecCc;
            } else if (unlikely(strequal_k(cur_modif, "standard-beta", cur_modif_slen))) {
              logerrputs("Error: --glm 'standard-beta' modifier has been retired.  Use\n--{covar-}variance-standardize instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (strequal_k(cur_modif, "perm", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmPerm;
            } else if (strequal_k(cur_modif, "perm-count", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmPermCount;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (unlikely(pc.glm_info.cols)) {
                logerrputs("Error: Multiple --glm cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0maybeprovref\0provref\0omitted\0a1count\0totallele\0a1countcc\0totallelecc\0gcountcc\0a1freq\0a1freqcc\0machr2\0firth\0test\0nobs\0beta\0orbeta\0se\0ci\0tz\0p\0err\0ax\0", "glm", kfGlmColChrom, kfGlmColDefault, 1, &pc.glm_info.cols);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              if (unlikely((!(pc.glm_info.cols & (kfGlmColBeta | kfGlmColOrbeta))) && ((pc.glm_info.cols & kfGlmColSe) || ((pc.glm_info.cols & kfGlmColCi) && (pc.ci_size != 0))))) {
                logerrputs("Error: --glm's 'se' and 'ci' columns require beta/orbeta to be included.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (pc.glm_info.cols & kfGlmColAx) {
                logerrputs("Warning: --glm 'ax' column is deprecated; you probably want 'omitted' instead.\nIf you think you need ax's current behavior, contact the developers and\ndescribe your use case.\n");
              }
            } else if (StrStartsWith0(cur_modif, "mperm", cur_modif_slen)) {
              if (unlikely((cur_modif_slen < 7) || (cur_modif[5] != '='))) {
                logerrputs("Error: Improper --glm mperm syntax.  (Use --glm mperm=[value]'.)\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely(ScanPosintDefcapx(&(cur_modif[6]), &pc.glm_info.mperm_ct))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --glm mperm= argument '%s'.\n", &(cur_modif[6]));
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (StrStartsWith(cur_modif, "local-covar=", cur_modif_slen)) {
              if (unlikely(pc.glm_local_covar_fname)) {
                logerrputs("Error: Multiple --glm local-covar= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = AllocFname(&(cur_modif[strlen("local-covar=")]), "glm local-covar=", &pc.glm_local_covar_fname);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "local-bim=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "local-pvar=", cur_modif_slen)) {
              if (unlikely(pc.glm_local_pvar_fname)) {
                logerrputs("Error: Multiple --glm local-pvar= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const uint32_t is_pvar = (cur_modif[6] == 'p');
              reterr = AllocFname(&(cur_modif[10 + is_pvar]), "glm local-pvar=", &pc.glm_local_pvar_fname);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "local-fam=", cur_modif_slen) ||
                       StrStartsWith(cur_modif, "local-psam=", cur_modif_slen)) {
              if (unlikely(pc.glm_local_psam_fname)) {
                logerrputs("Error: Multiple --glm local-psam= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const uint32_t is_psam = (cur_modif[6] == 'p');
              reterr = AllocFname(&(cur_modif[10 + is_psam]), "glm local-psam=", &pc.glm_local_psam_fname);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (strequal_k(cur_modif, "local-omit-last", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmLocalOmitLast;
            } else if (strequal_k(cur_modif, "local-haps", cur_modif_slen)) {
              pc.glm_info.flags |= kfGlmLocalHaps;
            } else if (StrStartsWith(cur_modif, "local-pos-cols=", cur_modif_slen)) {
              const char* cur_modif_iter = &(cur_modif[strlen("local-pos-cols=")]);
              uint32_t header_line_ct;
              if (unlikely(ScanmovPosintCapped(0x7ffffffe, &cur_modif_iter, &header_line_ct) || (*cur_modif_iter != ','))) {
                logerrputs("Error: Invalid local-pos-cols= value sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.glm_info.local_header_line_ct = header_line_ct;
              ++cur_modif_iter;
              uint32_t chrom_col_idx;
              if (unlikely(ScanmovPosintCapped(kMaxLongLine / 2, &cur_modif_iter, &chrom_col_idx) || (*cur_modif_iter != ','))) {
                logerrputs("Error: Invalid local-pos-cols= value sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.glm_info.local_chrom_col = chrom_col_idx;
              ++cur_modif_iter;
              uint32_t bp_col_idx;
              if (unlikely(ScanmovPosintCapped(kMaxLongLine / 2, &cur_modif_iter, &bp_col_idx) || (bp_col_idx == chrom_col_idx) || (*cur_modif_iter != ','))) {
                logerrputs("Error: Invalid local-pos-cols= value sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.glm_info.local_bp_col = bp_col_idx;
              ++cur_modif_iter;
              uint32_t first_covar_col;
              if (unlikely(ScanmovPosintCapped(kMaxLongLine / 2, &cur_modif_iter, &first_covar_col) || (*cur_modif_iter) || (first_covar_col <= chrom_col_idx) || (first_covar_col <= bp_col_idx))) {
                logerrputs("Error: Invalid local-pos-cols= value sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.glm_info.local_first_covar_col = first_covar_col;
            } else if (StrStartsWith(cur_modif, "local-cats=", cur_modif_slen)) {
              if (unlikely(pc.glm_info.local_cat_ct)) {
                logerrputs("Error: Multiple --glm local-cats[0]= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              // bugfix (7 Nov 2017): forgot to offset by strlen("local-cats=")
              if (unlikely(ScanPosintCappedx(&(cur_modif[strlen("local-cats=")]), 4095, &pc.glm_info.local_cat_ct) || (pc.glm_info.local_cat_ct == 1))) {
                logerrputs("Error: Invalid --glm local-cats= category count (must be in [2, 4095]).\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.glm_info.flags |= kfGlmLocalCats1based;
            } else if (StrStartsWith(cur_modif, "local-cats0=", cur_modif_slen)) {
              if (unlikely(pc.glm_info.local_cat_ct)) {
                logerrputs("Error: Multiple --glm local-cats[0]= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              if (unlikely(ScanPosintCappedx(&(cur_modif[strlen("local-cats0=")]), 4095, &pc.glm_info.local_cat_ct) || (pc.glm_info.local_cat_ct == 1))) {
                logerrputs("Error: Invalid --glm local-cats0= category count (must be in [2, 4095]).\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else if (likely(strequal_k(cur_modif, "allow-no-covars", cur_modif_slen))) {
              pc.glm_info.flags |= kfGlmAllowNoCovars;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --glm argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!pc.glm_info.cols) {
            pc.glm_info.cols = kfGlmColDefault;
          }
          if (unlikely((explicit_firth_fallback && (pc.glm_info.flags & (kfGlmNoFirth | kfGlmFirth))) || ((pc.glm_info.flags & (kfGlmNoFirth | kfGlmFirth)) == (kfGlmNoFirth | kfGlmFirth)))) {
            logerrputs("Error: Conflicting --glm arguments.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.glm_info.flags & (kfGlmSex | kfGlmNoXSex)) == (kfGlmSex | kfGlmNoXSex))) {
            logerrputs("Error: Conflicting --glm arguments.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.glm_info.flags & kfGlmPerm) && pc.glm_info.mperm_ct)) {
            logerrputs("Error: --glm 'perm' and 'mperm=' cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.glm_info.flags & kfGlmResidualizeMask) {
            if ((pc.glm_info.flags & (kfGlmFirthResidualize | kfGlmCcResidualize)) == (kfGlmFirthResidualize | kfGlmCcResidualize)) {
              logputs("Note: 'firth-residualize' is redundant when 'cc-residualize' is already\nspecified.\n");
              pc.glm_info.flags ^= kfGlmFirthResidualize;
            }
            if (unlikely(!(pc.glm_info.flags & kfGlmHideCovar))) {
              logerrputs("Error: --glm '{cc,firth,qt}-residualize' requires 'hide-covar' to be specified\nas well.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(pc.glm_info.flags & kfGlmInteraction)) {
              logerrputs("Error: --glm '{cc,firth,qt}-residualize' cannot be used with 'interaction'.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(pc.glm_info.flags & kfGlmIntercept)) {
              logerrputs("Error: --glm '{cc,firth,qt}-residualize' cannot be used with 'intercept'.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(pc.glm_local_covar_fname)) {
              logerrputs("Error: --glm '{cc,firth,qt}-residualize' cannot be used with local covariates.\n(If you want to include a per-chromosome polygenic effect, you need to include\nthat as a regular covariate and perform per-chromosome --glm runs; sorry about\nthe inconvenience.)\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely((pc.glm_info.flags & (kfGlmFirthResidualize | kfGlmNoFirth)) == (kfGlmFirthResidualize | kfGlmNoFirth))) {
              logerrputs("Error: --glm 'firth-residualize' doesn't make sense with 'no-firth'.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          }
          uint32_t alternate_genotype_col_flags = S_CAST(uint32_t, pc.glm_info.flags & (kfGlmGenotypic | kfGlmHethom | kfGlmDominant | kfGlmRecessive | kfGlmHetonly));
          if (alternate_genotype_col_flags) {
            if (unlikely(alternate_genotype_col_flags & (alternate_genotype_col_flags - 1))) {
              logerrputs("Error: Conflicting --glm arguments.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          if (unlikely((pc.glm_info.flags & kfGlmIntercept) && (!(pc.glm_info.cols & kfGlmColTest)))) {
            logerrputs("Error: --glm 'intercept' modifier cannot be used with an omitted 'test' column.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (!pc.glm_local_covar_fname) {
            if (unlikely(pc.glm_local_pvar_fname || pc.glm_local_psam_fname)) {
              logerrputs("Error: Either all three --glm local-covar filenames must be specified, or none\nof them.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(pc.glm_info.flags & kfGlmLocalOmitLast)) {
              logerrputs("Error: --glm 'local-omit-last' must be used with 'local-covar='.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(pc.glm_info.flags & kfGlmLocalHaps)) {
              logerrputs("Error: --glm 'local-haps' must be used with 'local-covar='.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(pc.glm_info.local_cat_ct)) {
              logerrputs("Error: --glm 'local-cats[0]=' must be used with 'local-covar='.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely((!pc.covar_fname) && (!pc.covar_range_list.name_ct) && (!(pc.glm_info.flags & kfGlmAllowNoCovars)))) {
              logerrputs("Error: --glm invoked without --covar/--covar-name/--covar-col-nums; this is\nusually an analytical mistake.  Add the 'allow-no-covars' modifier if you are\nsure you want this.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else {
            if (unlikely(pc.glm_local_pvar_fname && pc.glm_info.local_first_covar_col)) {
              logerrputs("Error: --glm local-pvar= and local-pos-cols= cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely((!pc.glm_local_psam_fname) || ((!pc.glm_local_pvar_fname) && (!pc.glm_info.local_first_covar_col)))) {
              logerrputs("Error: --glm local-covar= must be used with local-psam= and either local-pvar=\nor local-pos-cols=.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          pc.command_flags1 |= kfCommand1Glm;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "en")) {
          if (unlikely(load_params || (xload & (~kfXloadOxSample)))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(param_ct == 1)) {
            logerrputs("Error: --gen now requires a REF/ALT mode ('ref-first', 'ref-last', or\n'ref-unknown').\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 2];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          if (strequal_k(cur_modif, "ref-first", cur_modif_slen)) {
            oxford_import_flags |= kfOxfordImportRefFirst;
          } else if (strequal_k(cur_modif, "ref-unknown", cur_modif_slen)) {
            oxford_import_flags |= kfOxfordImportRefUnknown;
          } else if (likely(strequal_k(cur_modif, "ref-last", cur_modif_slen) ||
                            strequal_k(cur_modif, "ref-second", cur_modif_slen))) {
            oxford_import_flags |= kfOxfordImportRefLast;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --gen argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          const char* cur_fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --gen filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_fname, slen + 1);
          xload |= kfXloadOxGen;
        } else if (strequal_k_unsafe(flagname_p2, "was-ssf")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 7))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.gwas_ssf_info.flags |= kfGwasSsfZs;
            } else if (strequal_k(cur_modif, "delete-orig-glm", cur_modif_slen)) {
              if (unlikely(!(pc.command_flags1 & kfCommand1Glm))) {
                logerrputs("Error: --gwas-ssf 'delete-orig-glm' was specified, but --glm was not.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              pc.gwas_ssf_info.flags |= kfGwasSsfDeleteOrigGlm;
            } else if (strequal_k(cur_modif, "real-ref-alleles", cur_modif_slen)) {
              pc.gwas_ssf_info.flags |= kfGwasSsfRealRefAlleles;
            } else if (strequal_k(cur_modif, "allow-ambiguous-indels", cur_modif_slen)) {
              pc.gwas_ssf_info.flags |= kfGwasSsfAllowAmbiguousIndels;
            } else if (StrStartsWith(cur_modif, "a1freq-lower-limit=", cur_modif_slen)) {
              const char* lower_limit_start = &(cur_modif[strlen("a1freq-lower-limit=")]);
              double dxx;
              if (unlikely((!ScantokDouble(lower_limit_start, &dxx)) || (dxx < 0.0) || (dxx > 1.0))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --gwas-ssf a1freq-lower-bound= argument '%s'.\n", lower_limit_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              pc.gwas_ssf_info.a1freq_lower_limit = dxx * (1.0 - kSmallEpsilon);
            } else if (StrStartsWith(cur_modif, "file=", cur_modif_slen)) {
              reterr = AllocFname(&(cur_modif[strlen("file=")]), argvk[arg_idx], &pc.gwas_ssf_info.fname);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "file-list=", cur_modif_slen)) {
              reterr = AllocFname(&(cur_modif[strlen("file-list=")]), argvk[arg_idx], &pc.gwas_ssf_info.list_fname);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (likely(StrStartsWith(cur_modif, "rsid=", cur_modif_slen))) {
              if (unlikely(pc.gwas_ssf_info.rsid_mode != kGwasSsfRsidMode0)) {
                logerrputs("Error: Multiple --gwas-ssf rsid= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* mode_start = &(cur_modif[strlen("rsid=")]);
              const uint32_t mode_slen = cur_modif_slen - strlen("rsid=");
              if (strequal_k(mode_start, "no", mode_slen)) {
                pc.gwas_ssf_info.rsid_mode = kGwasSsfRsidModeNo;
              } else if (strequal_k(mode_start, "infer", mode_slen)) {
                pc.gwas_ssf_info.rsid_mode = kGwasSsfRsidModeInfer;
              } else if (likely(strequal_k(mode_start, "yes", mode_slen))) {
                pc.gwas_ssf_info.rsid_mode = kGwasSsfRsidModeYes;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --gwas-ssf rsid= mode '%s' ('no'/'infer'/'yes' expected).\n", mode_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --gwas-ssf argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (unlikely((pc.gwas_ssf_info.flags & (kfGwasSsfRealRefAlleles | kfGwasSsfAllowAmbiguousIndels)) == (kfGwasSsfRealRefAlleles | kfGwasSsfAllowAmbiguousIndels))) {
            logerrputs("Error: --gwas-ssf 'real-ref-alleles' and 'allow-ambiguous-indels' modifiers\ndon't make sense together.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (pc.command_flags1 & kfCommand1Glm) {
            if (unlikely(pc.glm_info.flags & (kfGlmHethom | kfGlmDominant | kfGlmRecessive | kfGlmHetonly))) {
              logerrputs("Error: --glm output does not contain an 'ADD' test, which is required by\n--gwas-ssf.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(((pc.glm_info.cols & (kfGlmColGwasSsfReq | kfGlmColProvref)) != (kfGlmColGwasSsfReq | kfGlmColProvref)) || (!(pc.glm_info.cols & (kfGlmColBeta | kfGlmColOrbeta))))) {
              logerrputs("Error: --glm column set is inappropriate for --gwas-ssf.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (pc.gwas_ssf_info.flags & kfGwasSsfDeleteOrigGlm) {
              if ((pc.glm_info.flags & (kfGlmZs | kfGlmHideCovar)) != (kfGlmZs | kfGlmHideCovar)) {
                logputs("Note: --gwas-ssf delete-orig-glm: automatically setting --glm 'zs' and\n'hide-covar' to reduce temporary file size.\n");
                pc.glm_info.flags |= kfGlmZs | kfGlmHideCovar;
              }
            }
            if (unlikely((pc.gwas_ssf_info.flags & kfGwasSsfRealRefAlleles) && (!pc.gwas_ssf_info.fname) && (!pc.gwas_ssf_info.list_fname))) {
              logerrputs("Error: --gwas-ssf 'real-ref-alleles' modifier was specified without file= or\nfile-list= input.  It has no effect in this context, since the concurrent --glm\nrun includes the provref column in its output.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else if (unlikely((!pc.gwas_ssf_info.fname) && (!pc.gwas_ssf_info.list_fname))) {
            logerrputs("Error: --gwas-ssf specified without any input (--glm, file=, and/or file-list=\nrequired).\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.gwas_ssf_info.rsid_mode == kGwasSsfRsidMode0) {
            pc.gwas_ssf_info.rsid_mode = kGwasSsfRsidModeInfer;
          }
        } else if (likely(strequal_k_unsafe(flagname_p2, "enotyping-rate"))) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (unlikely(strcmp("dosage", cur_modif))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --genotyping-rate argument '%s'.\n", cur_modif);
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
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 5))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.hardy_flags |= kfHardyZs;
            } else if (strequal_k(cur_modif, "midp", cur_modif_slen)) {
              pc.hardy_flags |= kfHardyMidp;
            } else if (strequal_k(cur_modif, "log10", cur_modif_slen)) {
              pc.hardy_flags |= kfHardyLog10;
            } else if (strequal_k(cur_modif, "redundant", cur_modif_slen)) {
              pc.hardy_flags |= kfHardyRedundant;
            } else if (likely(StrStartsWith(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(pc.hardy_flags & kfHardyColAll)) {
                logerrputs("Error: Multiple --hardy cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0maybeprovref\0provref\0ax\0gcounts\0gcount1col\0hetfreq\0sexaf\0femalep\0p\0", "hardy", kfHardyColChrom, kfHardyColDefault, 1, &pc.hardy_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              if (unlikely((pc.hardy_flags & (kfHardyColGcounts | kfHardyColGcount1col)) == (kfHardyColGcounts | kfHardyColGcount1col))) {
                logerrputs("Error: --hardy's gcounts and gcounts1col column sets are mutually exclusive.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --hardy argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.hardy_flags & kfHardyColAll)) {
            pc.hardy_flags |= kfHardyColDefault;
          }
          pc.command_flags1 |= kfCommand1Hardy;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "we")) {
          if (unlikely(pc.flip_info.subset_fname || (make_plink2_flags & kfMakePlink2FillMissingWithRef))) {
            logerrputs("Error: --hwe cannot be used with --fill-missing-with-ref/--flip-subset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t ln_thresh_seen = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "midp", cur_modif_slen)) {
              pc.misc_flags |= kfMiscHweMidp;
            } else if (strequal_k(cur_modif, "keep-fewhet", cur_modif_slen)) {
              pc.misc_flags |= kfMiscHweKeepFewhet;
            } else if (!ln_thresh_seen) {
              if (unlikely(!ScantokLn(cur_modif, &pc.hwe_ln_thresh))) {
                logerrputs("Error: Invalid --hwe argument sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely(pc.hwe_ln_thresh > 0.0)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --hwe threshold '%s' (must be in [0, 1]).\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              ln_thresh_seen = 1;
            } else {
              if (unlikely((pc.hwe_sample_size_term != -1) || (!ScantokDouble(cur_modif, &pc.hwe_sample_size_term)))) {
                logerrputs("Error: Invalid --hwe argument sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely(pc.hwe_sample_size_term < 0.0)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --hwe sample size term '%s' (must be nonnegative).\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            }
          }
          if (unlikely(!ln_thresh_seen)) {
            logerrputs("Error: No --hwe p-value threshold specified.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.hwe_ln_thresh == 0.0) && (pc.hwe_sample_size_term <= 0.0))) {
            logerrputs("Error: --hwe threshold cannot be 1 unless a nonzero sample size term is\nspecified.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if ((pc.misc_flags & kfMiscHweMidp) && (pc.hwe_ln_thresh >= -kLn2)) {
            logerrputs("Error: --hwe threshold must be smaller than 0.5 when using mid-p adjustment.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "ard-call-threshold")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double hard_call_frac;
          if (unlikely((!ScantokDouble(cur_modif, &hard_call_frac)) || (hard_call_frac < 0.0) || (hard_call_frac >= (0.5 - kSmallEpsilon)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --hard-call-threshold argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.hard_call_thresh = S_CAST(int32_t, hard_call_frac * ((1 + kSmallEpsilon) * kDosageMid));
        } else if (strequal_k_unsafe(flagname_p2, "orse")) {
          if (unlikely(chr_info.chrset_source)) {
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
          if (unlikely(chr_info.chrset_source)) {
            logerrputs("Error: Conflicting chromosome-set flags.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          chr_info.chrset_source = kChrsetSourceCmdline;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "et")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.het_flags |= kfHetZs;
            } else if (strequal_k(cur_modif, "small-sample", cur_modif_slen)) {
              pc.het_flags |= kfHetSmallSample;
            } else if (likely(StrStartsWith(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(pc.het_flags & kfHetColAll)) {
                logerrputs("Error: Multiple --het cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0hom\0het\0nobs\0f\0", "het", kfHetColMaybefid, kfHetColDefault, 1, &pc.het_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --het argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.het_flags & kfHetColAll)) {
            pc.het_flags |= kfHetColDefault;
          }
          pc.command_flags1 |= kfCommand1Het;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (likely(strequal_k_unsafe(flagname_p2, "aps"))) {
          if (unlikely(load_params || xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 2) {
            const char* cur_modif = argvk[arg_idx + 2];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "ref-first", cur_modif_slen)) {
              oxford_import_flags |= kfOxfordImportRefFirst;
            } else if (likely(strequal_k(cur_modif, "ref-last", cur_modif_slen) ||
                              strequal_k(cur_modif, "ref-second", cur_modif_slen))) {
              oxford_import_flags |= kfOxfordImportRefLast;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --haps argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const char* cur_fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
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
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* mode_str = argvk[arg_idx + 1];
          const char first_char_upcase_match = mode_str[0] & 0xdf;
          const uint32_t mode_slen = strlen(mode_str);
          if (strequal_k(mode_str, "0", mode_slen) ||
              strequal_k(mode_str, "none", mode_slen)) {
            pc.sample_sort_mode = kSortNone;
          } else if (((mode_slen == 1) && (first_char_upcase_match == 'N')) ||
                     strequal_k(mode_str, "natural", mode_slen)) {
            pc.sample_sort_mode = kSortNatural;
          } else if (((mode_slen == 1) && (first_char_upcase_match == 'A')) ||
                     strequal_k(mode_str, "ascii", mode_slen)) {
            pc.sample_sort_mode = kSortAscii;
          } else if (likely(((mode_slen == 1) && (first_char_upcase_match == 'F')) ||
                            strequal_k(mode_str, "file", mode_slen))) {
            if (unlikely(param_ct == 1)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Missing '--indiv-sort %s' filename.\n", mode_str);
              goto main_ret_INVALID_CMDLINE_2A;
            }
            pc.sample_sort_mode = kSortFile;
            reterr = AllocFname(argvk[arg_idx + 2], flagname_p, &pc.sample_sort_fname);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a valid mode for --indiv-sort.\n", mode_str);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely((param_ct > 1) && (pc.sample_sort_mode != kSortFile))) {
            snprintf(g_logbuf, kLogbufSize, "Error: '--indiv-sort %s' does not accept additional arguments.\n", mode_str);
            goto main_ret_INVALID_CMDLINE_2A;
          }
        } else if (strequal_k_unsafe(flagname_p2, "d-delim")) {
          if (unlikely(const_fid || (import_flags & kfImportDoubleId))) {
            logerrputs("Error: --id-delim can no longer be used with --const-fid or --double-id.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            id_delim = ExtractCharParam(cur_modif);
            if (unlikely(!id_delim)) {
              logerrputs("Error: --id-delim delimiter must be a single character.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely(ctou32(id_delim) < ' ')) {
              logerrputs("Error: --id-delim argument cannot be tab, newline, or a nonprinting character.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          }
          if (!id_delim) {
            id_delim = '_';
          }
        } else if (strequal_k_unsafe(flagname_p2, "ndep-pairwise") || strequal_k_unsafe(flagname_p2, "ndep-pairphase")) {
          if (unlikely(pc.command_flags1 & kfCommand1LdPrune)) {
            logerrputs("Error: Multiple LD pruning commands.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double first_paramd;
          const char* first_param_end = ScanadvDouble(cur_modif, &first_paramd);
          if (unlikely((!first_param_end) || (first_paramd < 0.0))) {
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
              if (unlikely(pc.ld_info.prune_window_size < 2)) {
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
          if (unlikely(next_param_idx + 2 == param_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s argument sequence.\n", flagname_p);
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (next_param_idx < param_ct) {
            // explicit step size
            cur_modif = argvk[arg_idx + next_param_idx];
            if (unlikely(ScanPosintDefcapx(cur_modif, &pc.ld_info.prune_window_incr))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s window-increment '%s'.\n", flagname_p, cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (!is_kb) {
              if (unlikely(pc.ld_info.prune_window_incr > pc.ld_info.prune_window_size)) {
                snprintf(g_logbuf, kLogbufSize, "Error: --%s window-increment cannot be larger than window size.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
            } else if (unlikely(pc.ld_info.prune_window_incr != 1)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --%s window-increment must be 1 when window size is in\nkilobase units.\n", flagname_p);
              goto main_ret_INVALID_CMDLINE_2A;
            }
          } else {
            pc.ld_info.prune_window_incr = 1;
          }
          cur_modif = argvk[arg_idx + param_ct];
          if (unlikely((!ScantokDouble(cur_modif, &pc.ld_info.prune_last_param)) || (pc.ld_info.prune_last_param < 0.0) || (pc.ld_info.prune_last_param >= 1.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s r^2 threshold '%s'.\n", flagname_p, cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.command_flags1 |= kfCommand1LdPrune;
          pc.dependency_flags |= kfFilterAllReq;
          if (flagname_p2[9] == 'p') {
            pc.ld_info.prune_flags |= kfLdPrunePairphase;
          } else {
            pc.ld_info.prune_flags |= kfLdPrunePairwise;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ndep-preferred")) {
          if (!(pc.ld_info.prune_flags & (kfLdPrunePairwise | kfLdPrunePairphase))) {
            logerrputs("Error: --indep-preferred must be used with --indep-pairwise or\n--indep-pairphase.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.indep_preferred_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ndep-order")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          if (strequal_k(cur_modif, "1", cur_modif_slen)) {
            pc.ld_info.prune_flags ^= kfLdPrunePlink1Order;
          } else if (unlikely(!strequal_k(cur_modif, "2", cur_modif_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --indep-order mode '%s' ('1' or '2' expected).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "nput-missing-genotype")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          pc.input_missing_geno_char = ExtractCharParam(cur_modif);
          if (unlikely(ctou32(pc.input_missing_geno_char) <= ' ')) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --input-missing-genotype argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "nput-missing-phenotype")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (unlikely(ScanInt32x(cur_modif, &pc.missing_pheno) || ((pc.missing_pheno >= 0) && (pc.missing_pheno <= 2)) || (!ScantokDouble(cur_modif, &dxx)) || (dxx != S_CAST(double, pc.missing_pheno)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --input-missing-phenotype argument '%s' (must be an integer in [-2147483647, -1] or [3, 2147483647]).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (pc.missing_pheno == -9) {
            // No reason to require --neg9-pheno-really-missing when user
            // explicitly sets --input-missing-phenotype to -9.
            pc.misc_flags |= kfMiscNeg9PhenoReallyMissing;
          }
        } else if (strequal_k_unsafe(flagname_p2, "mport-dosage-certainty")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely((!ScantokDouble(cur_modif, &import_dosage_certainty)) || (import_dosage_certainty < 0.0) || (import_dosage_certainty > 1.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage-certainty argument '%s' (must be in [0, 1]).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          import_dosage_certainty *= 1.0 - kSmallEpsilon;
          // We may as well enforce --hard-call-threshold +
          // --import-dosage-certainty <= 1.
          uint32_t hard_call_thresh = pc.hard_call_thresh;
          if (hard_call_thresh == UINT32_MAX) {
            hard_call_thresh = kDosageMid / 10;
          }
          if (u31tod(hard_call_thresh) + import_dosage_certainty * kDosageMid >= u31tod(kDosageMid)) {
            logerrputs("Error: --hard-call-threshold + --import-dosage-certainty settings cannot add up\nto more than 1.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else if (strequal_k_unsafe(flagname_p2, "id-sid")) {
          pc.misc_flags |= kfMiscIidSid;
        } else if (strequal_k_unsafe(flagname_p2, "ndv")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, kMaxIdBlen, &pc.indv_str);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "mport-max-alleles")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(ScanPosintDefcapx(argvk[arg_idx + 1], &import_max_allele_ct) || (import_max_allele_ct < 2))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-max-alleles argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.load_filter_log_flags |= kfLoadFilterLogImportMaxAlleles;
        } else if (strequal_k_unsafe(flagname_p2, "mpute-sex")) {
          if (unlikely(pc.command_flags1 & kfCommand1CheckOrImputeSex)) {
            logerrputs("Error: --check-sex is redundant with --impute-sex.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.check_sex_info.flags |= kfCheckSexImpute;
          goto main_parse_check_sex;
        } else if (strequal_k_unsafe(flagname_p2, "dspace-to")) {
          if (unlikely(id_delim == ' ')) {
            logerrputs("Error: --idspace-to cannot be used when the --id-delim character is space.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          idspace_to = ExtractCharParam(argvk[arg_idx + 1]);
          if (unlikely(!idspace_to)) {
            logerrputs("Error: --idspace-to argument must be a single character.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(ctou32(idspace_to) <= ' ')) {
            logerrputs("Error: --idspace-to argument must be a nonspace character.\n");
            goto main_ret_INVALID_CMDLINE;
          }
        } else if (likely(strequal_k_unsafe(flagname_p2, "mport-dosage"))) {
          if (unlikely(load_params || xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 11))) {
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
              if (unlikely(plink1_dosage_info.skips[skip_idx])) {
                logerrprintf("Error: Multiple --import-dosage skip%u= modifiers.\n", skip_idx);
                goto main_ret_INVALID_CMDLINE;
              }
              if (unlikely(ScanUintCappedx(&(cur_modif[6]), kMaxLongLine / 2, &(plink1_dosage_info.skips[skip_idx])))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage skip%u= argument '%s'.\n", skip_idx, &(cur_modif[6]));
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "dose1", cur_modif_slen)) {
              plink1_dosage_info.flags |= kfPlink1DosageFormatSingle01;
            } else if (StrStartsWithUnsafe(cur_modif, "format=")) {
              // strequal_k() and StrStartsWith() both suboptimal here
              if (unlikely(format_num_m1 != 4)) {
                logerrputs("Error: Multiple --import-dosage format= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              if (cur_modif_slen == 8) {
                format_num_m1 = ctou32(cur_modif[7]) - 49;
                if (unlikely(format_num_m1 >= 3)) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage format= argument '%c'.\n", cur_modif[7]);
                  goto main_ret_INVALID_CMDLINE_2A;
                }
              } else if (likely(strequal_k(&(cur_modif[7]), "infer", cur_modif_slen - 7))) {
                format_num_m1 = 3;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage format= argument '%s'.\n", &(cur_modif[7]));
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "ref-first", cur_modif_slen)) {
              plink1_dosage_info.flags |= kfPlink1DosageRefFirst;
            } else if (strequal_k(cur_modif, "ref-last", cur_modif_slen) ||
                       strequal_k(cur_modif, "ref-second", cur_modif_slen)) {
              plink1_dosage_info.flags |= kfPlink1DosageRefLast;
            } else if (StrStartsWith(cur_modif, "id-delim=", cur_modif_slen)) {
              if (unlikely(plink1_dosage_info.id_delim)) {
                logerrputs("Error: Multiple --import-dosage id-delim= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* id_delim_str = &(cur_modif[strlen("id-delim=")]);
              char cc = ExtractCharParam(id_delim_str);
              if (unlikely(!cc)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage id-delim= argument '%s'.\n", id_delim_str);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              plink1_dosage_info.id_delim = cc;
            } else if (StrStartsWith(cur_modif, "single-chr=", cur_modif_slen)) {
              if (unlikely(import_single_chr_str)) {
                logerrputs("Error: Multiple --import-dosage single-chr= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* chr_code = &(cur_modif[strlen("single-chr=")]);
              if (pc.misc_flags & kfMiscProhibitExtraChr) {
                if (unlikely(IsI32Neg(GetChrCodeRaw(chr_code)))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage single-chr= chromosome code '%s'.\n", chr_code);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
              }
              reterr = CmdlineAllocString(chr_code, argvk[arg_idx], kMaxIdSlen, &import_single_chr_str);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "chr-col-num=", cur_modif_slen)) {
              if (unlikely(plink1_dosage_info.chr_col_idx != UINT32_MAX)) {
                logerrputs("Error: Multiple --import-dosage chr-col-num= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* chr_col_num_start = &(cur_modif[strlen("chr-col-num=")]);
              uint32_t uii;
              if (unlikely(ScanPosintCappedx(chr_col_num_start, kMaxLongLine / 2, &uii))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage chr-col-num= argument '%s'.\n", chr_col_num_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              plink1_dosage_info.chr_col_idx = uii - 1;
            } else if (likely(StrStartsWith(cur_modif, "pos-col-num=", cur_modif_slen))) {
              if (unlikely(plink1_dosage_info.pos_col_idx != UINT32_MAX)) {
                logerrputs("Error: Multiple --import-dosage pos-col-num= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* pos_col_num_start = &(cur_modif[strlen("pos-col-num=")]);
              uint32_t uii;
              if (unlikely(ScanPosintCappedx(pos_col_num_start, kMaxLongLine / 2, &uii))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage pos-col-num= argument '%s'.\n", pos_col_num_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              plink1_dosage_info.pos_col_idx = uii - 1;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --import-dosage argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }

          if (!format_num_m1) {
            plink1_dosage_info.flags |= kfPlink1DosageFormatSingle;
          } else {
            if (plink1_dosage_info.flags & kfPlink1DosageFormatSingle01) {
              if (unlikely(format_num_m1 < 3)) {
                logerrputs("Error: --import-dosage 'dose1' modifier cannot be used with 'format=2' or\n'format=3'.\n");
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
          if (unlikely((plink1_dosage_info.flags & (kfPlink1DosageRefFirst | kfPlink1DosageRefLast)) == (kfPlink1DosageRefFirst | kfPlink1DosageRefLast))) {
            logerrputs("Error: --import-dosage 'ref-first' and 'ref-last' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          const uint32_t id_col_idx = plink1_dosage_info.skips[0];
          const uint32_t a1_col_idx = id_col_idx + plink1_dosage_info.skips[1] + 1;
          const uint32_t data_col_idx = a1_col_idx + plink1_dosage_info.skips[2] + 2;
          const uint32_t chr_col_idx = plink1_dosage_info.chr_col_idx;
          if (chr_col_idx != UINT32_MAX) {
            if (unlikely(import_single_chr_str)) {
              logerrputs("Error: --import-dosage 'single-chr=' and 'chr-col-num=' modifiers cannot be\nused together.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (unlikely((chr_col_idx == id_col_idx) || (chr_col_idx == a1_col_idx) || (chr_col_idx == a1_col_idx + 1))) {
              logerrputs("Error: --import-dosage chr-col-num= value collides with another column.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (unlikely(chr_col_idx >= data_col_idx)) {
              logerrputs("Error: --import-dosage chr-col-num= value too large.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const uint32_t pos_col_idx = plink1_dosage_info.pos_col_idx;
          if (pos_col_idx != UINT32_MAX) {
            if (unlikely((pos_col_idx == id_col_idx) || (pos_col_idx == a1_col_idx) || (pos_col_idx == a1_col_idx + 1) || (pos_col_idx == chr_col_idx))) {
              logerrputs("Error: --import-dosage pos-col-num= value collides with another column.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (unlikely(pos_col_idx >= data_col_idx)) {
              logerrputs("Error: --import-dosage pos-col-num= value too large.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_modif);
          if (unlikely(slen > kPglFnamesize - 1)) {
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
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.keep_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-fam")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.keepfam_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-autoconv")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (likely(!strcmp(cur_modif, "vzs"))) {
              import_flags |= kfImportKeepAutoconvVzs;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --keep-autoconv argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          import_flags |= kfImportKeepAutoconv;
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
          if (unlikely(pc.filter_flags & kfFilterExclNonfounders)) {
            logerrputs("Error: --keep-nonfounders cannot be used with --keep-founders.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.filter_flags |= kfFilterPsamReq | kfFilterExclFounders;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ing-cutoff") || strequal_k_unsafe(flagname_p2, "ing-cutoff-table")) {
          const uint32_t is_table = (flagname_p2[strlen("ing-cutoff")] != '\0');
          if (is_table) {
            if (unlikely(pc.command_flags1 & kfCommand1KingCutoff)) {
              logerrputs("Error: --king-cutoff cannot be used with --king-cutoff-table.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            pc.king_flags |= kfKingCutoffTable;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1 + is_table, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 2) {
            if (is_table) {
              reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &king_cutoff_fprefix);
            } else {
              reterr = AllocFnamePrefix(argvk[arg_idx + 1], ".king.id\0.king.bin\0", flagname_p, &king_cutoff_fprefix);
            }
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            pc.dependency_flags |= kfFilterPsamReq;
          } else {
            pc.dependency_flags |= kfFilterAllReq;
          }
          const char* cur_modif = argvk[arg_idx + param_ct];
          if (unlikely((!ScantokDouble(cur_modif, &pc.king_cutoff)) || (pc.king_cutoff < 0.0) || (pc.king_cutoff >= 0.5))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --king-cutoff[-table] argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.command_flags1 |= kfCommand1KingCutoff;
        } else if (strequal_k_unsafe(flagname_p2, "ing-table-filter")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely((!ScantokDouble(cur_modif, &pc.king_table_filter)) || (pc.king_table_filter > 0.5))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --king-table-filter argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ing-table-subset")) {
          if (unlikely(pc.command_flags1 & kfCommand1KingCutoff)) {
            logerrputs("Error: --king-table-subset cannot be used with --king-cutoff.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.king_table_subset_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          if (param_ct == 2) {
            const char* cur_modif = argvk[arg_idx + 2];
            if (unlikely((!ScantokDouble(cur_modif, &pc.king_table_subset_thresh)) || (pc.king_table_subset_thresh > 0.5))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --king-table-subset threshold '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.king_table_subset_thresh = -DBL_MAX;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ing-table-require") || strequal_k_unsafe(flagname_p2, "ing-table-require-xor")) {
          if (unlikely(pc.command_flags1 & kfCommand1KingCutoff)) {
            logerrputs("Error: --king-table-require[-xor] cannot be used with --king-cutoff.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.king_table_require_fnames) {
            logerrputs("Error: --king-table-require cannot be used with --king-table-require-xor.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.king_table_require_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          if (flagname_p2[strlen("ing-table-require")] != '\0') {
            pc.king_flags |= kfKingTableRequireXor;
          }
        } else if (strequal_k_unsafe(flagname_p2, "eep-if")) {
          reterr = ValidateAndAllocCmpExpr(&(argvk[arg_idx + 1]), argvk[arg_idx], param_ct, &pc.keep_if_expr);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-cats")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.keep_cats_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-cat-names")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, kMaxIdBlen, &pc.keep_cat_names_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-cat-pheno")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.keep_cat_phenoname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "eep-col-match")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.keep_col_match_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 2]), flagname_p, param_ct - 1, kMaxIdBlen, &pc.keep_col_match_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "eep-col-match-name")) {
          if (unlikely(!pc.keep_col_match_fname)) {
            logerrputs("Error: --keep-col-match-name must be used with --keep-col-match.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.keep_col_match_name);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "eep-col-match-num")) {
          if (unlikely(!pc.keep_col_match_fname)) {
            logerrputs("Error: --keep-col-match-num must be used with --keep-col-match.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.keep_col_match_name)) {
            logerrputs("Error: --keep-col-match-num can't be used with --keep-col-match-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(ScanPosintDefcapx(cur_modif, &pc.keep_col_match_num) || (pc.keep_col_match_num == 1))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --keep-col-match-num argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (likely(strequal_k_unsafe(flagname_p2, "eep-allele-order"))) {
          if (unlikely((pc.command_flags1 & kfCommand1Glm) && (!(pc.glm_info.flags & kfGlmOmitRef)))) {
            // update (17 May 2018): Best to error out instead of ignore if
            // this is used with --linear/--logistic/--glm without 'omit-ref',
            // since in that case the user probably wants to add 'omit-ref'.
            logerrputs("Error: To make --glm always test ALT alleles, you must use --glm's 'omit-ref'\nmodifier, not --keep-allele-order.  (--keep-allele-order no longer has any\neffect, since plink2 always keeps track of REF/ALT alleles; but --glm defaults\nto testing minor instead of ALT alleles, since this can be necessary for\navoiding multicollinearity.)\n");
            goto main_ret_INVALID_CMDLINE;
          }
          logputs("Note: --keep-allele-order no longer has any effect.\n");
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'l':
        if (strequal_k_unsafe(flagname_p2, "ambda")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double lambda;
          if (unlikely(!ScantokDouble(argvk[arg_idx + 1], &lambda))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --lambda argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (lambda < 1.0) {
            logputs("Note: --lambda argument set to 1.\n");
            lambda = 1.0;
          }
          pc.adjust_info.lambda = lambda;
          adjust_file_info.base.lambda = lambda;
        } else if (strequal_k_unsafe(flagname_p2, "egend")) {
          if (unlikely(load_params || (xload & (~kfXloadOxHaps)))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(!xload)) {
            logerrputs("Error: --legend must be used with --haps.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_fname = argvk[arg_idx + 1];
          uint32_t slen = strlen(cur_fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --legend filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pvarname, cur_fname, slen + 1);
          const char* chr_code = argvk[arg_idx + 2];
          if (pc.misc_flags & kfMiscProhibitExtraChr) {
            if (unlikely(IsI32Neg(GetChrCodeRaw(chr_code)))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --legend chromosome code '%s'.\n", chr_code);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          reterr = CmdlineAllocString(chr_code, argvk[arg_idx], kMaxIdSlen, &import_single_chr_str);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          xload |= kfXloadOxLegend;
        } else if (strequal_k_unsafe(flagname_p2, "oop-cats")) {
          if (unlikely(pc.command_flags1 & kfCommand1Clump)) {
            logerrputs("Error: --loop-cats cannot currently be used with --clump.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.loop_cats_phenoname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ax-chrx-import")) {
          import_flags |= kfImportLaxChrX;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ax-bgen-import")) {
          if (unlikely(!(xload & kfXloadOxBgen))) {
            logerrputs("Error: --lax-bgen-import must be used with --bgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          import_flags |= kfImportLaxBgen;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "d-window")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t uii;
          if (unlikely(ScanUintDefcapx(cur_modif, &uii) || (uii < 2))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --ld-window argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.vcor_info.var_ct_radius = uii - 1;
          r2_required = 1;
        } else if (strequal_k_unsafe(flagname_p2, "d-window-kb")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (unlikely((!ScantokDouble(cur_modif, &dxx)) || (dxx < 0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --ld-window-kb argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          dxx *= 1000 * (1 + kSmallEpsilon);
          if (dxx > 2147483646) {
            pc.vcor_info.bp_radius = 2147483646;
          } else {
            pc.vcor_info.bp_radius = S_CAST(int32_t, dxx);
          }
          r2_required = 1;
        } else if (strequal_k_unsafe(flagname_p2, "d-window-cm")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (unlikely((!ScantokDouble(cur_modif, &dxx)) || (dxx < 0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --ld-window-cm argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.vcor_info.cm_radius = dxx * (1 + kSmallEpsilon);
          r2_required = 1;
        } else if (strequal_k_unsafe(flagname_p2, "d-window-r2")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (unlikely((!ScantokDouble(cur_modif, &dxx)) || (dxx > 1.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --ld-window-r2 argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (dxx > 0.0) {
            dxx *= 1 - kSmallEpsilon;
          }
          pc.vcor_info.min_r2 = dxx;
          r2_required = 1;
        } else if (strequal_k_unsafe(flagname_p2, "d-snp")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          // probably want a helper function for this
          const char* variant_id = argvk[arg_idx + 1];
          const uint32_t variant_id_blen = 1 + strlen(variant_id);
          char* new_buf;
          if (unlikely(pgl_malloc(variant_id_blen + 1, &new_buf))) {
            goto main_ret_NOMEM;
          }
          pc.vcor_info.ld_snp_range_list.names = new_buf;
          pc.vcor_info.ld_snp_range_list.name_max_blen = variant_id_blen;
          pc.vcor_info.ld_snp_range_list.name_ct = 1;
          memcpy(new_buf, variant_id, variant_id_blen);
          new_buf[variant_id_blen] = '\0';
          pc.vcor_info.ld_snp_range_list.starts_range = R_CAST(unsigned char*, &(new_buf[variant_id_blen]));
          r2_required = 1;
        } else if (strequal_k_unsafe(flagname_p2, "d-snps")) {
          if (unlikely(pc.vcor_info.ld_snp_list_fname || pc.vcor_info.ld_snp_range_list.name_ct)) {
            logerrputs("Error: --ld-snps cannot be used with --ld-snp or --ld-snp-list.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.vcor_info.ld_snp_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          r2_required = 1;
        } else if (strequal_k_unsafe(flagname_p2, "d-snp-list")) {
          if (unlikely(pc.vcor_info.ld_snp_range_list.name_ct)) {
            logerrputs("Error: --ld-snp-list cannot be used with --ld-snp.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.vcor_info.ld_snp_list_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          r2_required = 1;
        } else if (likely(strequal_k_unsafe(flagname_p2, "d"))) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t uii = 0; uii != 2; ++uii) {
            reterr = CmdlineAllocString(argvk[arg_idx + uii + 1], argvk[arg_idx], kMaxIdSlen, &(pc.ld_info.ld_console_varids[uii]));
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          for (uint32_t param_idx = 3; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (likely(!strcmp(cur_modif, "hwe-midp"))) {
              pc.ld_info.ld_console_flags |= kfLdConsoleHweMidp;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --ld argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.command_flags1 |= kfCommand1Ld;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "oop-assoc")) {
          logerrputs("Error: --loop-assoc is retired.  Use --within + --split-cat-pheno instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "ist-duplicate-vars")) {
          logerrputs("Error: --list-duplicate-vars is retired.  We recommend --set-all-var-ids +\n--rm-dup for variant deduplication.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'm':
        if (strequal_k_unsafe(flagname_p2, "emory")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t mb_modif_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "require", &mb_modif_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            memory_require = 1;
          }
          const char* mb_modif = argvk[arg_idx + mb_modif_idx];
          if (unlikely(ScanPosintptrx(mb_modif, &malloc_size_mib))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --memory argument '%s'.\n", mb_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(malloc_size_mib < kBigstackMinMib)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --memory argument '%s' (minimum %u).\n", mb_modif, kBigstackMinMib);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
#ifndef __LP64__
          if (unlikely(malloc_size_mib > kMalloc32bitMibMax)) {
            logerrprintf("Error: --memory argument too large for 32-bit version (max %u).\n", kMalloc32bitMibMax);
            goto main_ret_INVALID_CMDLINE;
          }
#endif
        } else if (strequal_k_unsafe(flagname_p2, "ake-bed")) {
          if (unlikely(pc.exportf_info.flags & kfExportfIndMajorBed)) {
            logerrputs("Error: --make-bed cannot be used with --export ind-major-bed.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "vzs", cur_modif_slen)) {
              make_plink2_flags |= kfMakeBimZs;
            } else if (strequal_k(cur_modif, "trim-alts", cur_modif_slen)) {
              make_plink2_flags |= kfMakePlink2TrimAlts;
            } else if (likely(StrStartsWith(cur_modif, "m=", cur_modif_slen) ||
                              StrStartsWith(cur_modif, "multiallelics=", cur_modif_slen))) {
              if (unlikely(make_plink2_flags & kfMakePlink2MMask)) {
                logerrputs("Error: Multiple --make-bed multiallelics= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* mode_start = (cur_modif[1] == '=')? (&(cur_modif[strlen("m=")])) : (&(cur_modif[strlen("multiallelics=")]));
              const uint32_t mode_slen = cur_modif_slen - S_CAST(uintptr_t, mode_start - cur_modif);
              if (strequal_k(mode_start, "-", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MSplitAll;
              } else if (strequal_k(mode_start, "-snps", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MSplitSnps;
              } else {
                // merge modes not supported here
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-bed multiallelics= split mode '%s'.\n", mode_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else {
              char* write_iter = strcpya_k(g_logbuf, "Error: Invalid --make-bed argument '");
              write_iter = memcpya(write_iter, cur_modif, cur_modif_slen);
              write_iter = strcpya_k(write_iter, "'.");
              if ((param_idx == 1) && (!outname_end)) {
                // the missing --out mistake is so common--I must have made it
                // over a hundred times by now--that a custom error message is
                // worthwhile.
                write_iter = strcpya_k(write_iter, " (Did you forget '--out'?)");
              }
              *write_iter++ = '\n';
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          make_plink2_flags |= kfMakeBed | kfMakeBim | kfMakeFam;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-bpgen")) {
          if (unlikely(make_plink2_flags & kfMakeBed)) {
            logerrputs("Error: --make-bpgen cannot be used with --make-bed.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(import_flags & kfImportKeepAutoconv)) {
            logerrputs("Error: --make-bpgen cannot be used with --keep-autoconv.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 8))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t varid_semicolon = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "vzs", cur_modif_slen)) {
              make_plink2_flags |= kfMakeBimZs;
            } else if (StrStartsWith(cur_modif, "format=", cur_modif_slen)) {
              if (unlikely(make_plink2_flags & (kfMakePgenFormatBase * 3))) {
                logerrputs("Error: Multiple --make-bpgen format= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const uint32_t fcode_minus_2 = ctou32(cur_modif[7]) - 50;
              if (unlikely((fcode_minus_2 > 2) || cur_modif[8])) {
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
              if (unlikely(make_plink2_flags & kfMakePlink2MMask)) {
                logerrputs("Error: Multiple --make-bpgen multiallelics= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              // er, some of this belongs in its own function...
              const char* mode_start = (cur_modif[1] == '=')? (&(cur_modif[strlen("m=")])) : (&(cur_modif[strlen("multiallelics=")]));
              const uint32_t mode_slen = cur_modif_slen - S_CAST(uintptr_t, mode_start - cur_modif);
              if (strequal_k(mode_start, "-", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MSplitAll;
              } else if (strequal_k(mode_start, "-snps", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MSplitSnps;
              } else if (strequal_k(mode_start, "+", mode_slen) ||
                         strequal_k(mode_start, "+both", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MJoinBoth;
              } else if (strequal_k(mode_start, "+snps", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MJoinSnps;
              } else if (likely(strequal_k(mode_start, "+any", mode_slen))) {
                make_plink2_flags |= kfMakePlink2MJoinAny;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-bpgen multiallelics= mode '%s'.\n", mode_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "trim-alts", cur_modif_slen)) {
              make_plink2_flags |= kfMakePlink2TrimAlts;
            } else if (strequal_k(cur_modif, "varid-split", cur_modif_slen)) {
              varid_semicolon |= 1;
            } else if (strequal_k(cur_modif, "varid-split-dup", cur_modif_slen)) {
              varid_semicolon |= 2;
            } else if (strequal_k(cur_modif, "varid-join", cur_modif_slen)) {
              varid_semicolon |= 4;
            } else if (strequal_k(cur_modif, "varid-dup", cur_modif_slen)) {
              make_plink2_flags |= kfMakePlink2VaridDup;
            } else if (strequal_k(cur_modif, "erase-phase", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenErasePhase;
            } else if (strequal_k(cur_modif, "erase-dosage", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenEraseDosage;
            } else if (strequal_k(cur_modif, "fill-missing-from-dosage", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenFillMissingFromDosage;
            } else if (likely(strequal_k(cur_modif, "writer-ver", cur_modif_slen))) {
              make_plink2_flags |= kfMakePgenWriterVer;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-bpgen argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (varid_semicolon) {
            if (unlikely((make_plink2_flags & kfMakePlink2VaridDup) || (varid_semicolon & (varid_semicolon - 1)))) {
              logerrputs("Error: --make-bpgen 'varid-split', 'varid-split-dup', 'varid-dup', and\n'varid-join' modifiers are mutually exclusive.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (varid_semicolon & 3) {
              if (unlikely((make_plink2_flags & kfMakePlink2MJoin) || (!(make_plink2_flags & kfMakePlink2MMask)))) {
                logerrputs("Error: --make-bpgen 'varid-split' must be used with a multiallelics= split\nmode.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else {
              if (unlikely(!(make_plink2_flags & kfMakePlink2MJoin))) {
                logerrputs("Error: --make-bpgen 'varid-join' must be used with a multiallelics= join mode.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            }
            make_plink2_flags |= kfMakePlink2VaridSemicolon;
            if (varid_semicolon == 2) {
              make_plink2_flags |= kfMakePlink2VaridDup;
            }
          } else if (make_plink2_flags & kfMakePlink2VaridDup) {
            if (unlikely((make_plink2_flags & kfMakePlink2MJoin) || (!(make_plink2_flags & kfMakePlink2MMask)))) {
              logerrputs("Error: --make-bpgen 'varid-dup' must be used with a multiallelics= split mode.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          make_plink2_flags |= kfMakePgen | kfMakeBim | kfMakeFam;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-pgen")) {
          if (unlikely(make_plink2_flags & (kfMakeBed | kfMakePgen))) {
            logerrputs("Error: --make-pgen cannot be used with --make-bed/--make-bpgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(make_plink2_flags & (kfMakeBim | kfMakeFam | kfMakePvar | kfMakePsam))) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-[b]pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(import_flags & kfImportKeepAutoconv)) {
            logerrputs("Error: --make-pgen cannot be used with --keep-autoconv.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 10))) {
            goto main_ret_INVALID_CMDLINE_A;
          }
          uint32_t explicit_pvar_cols = 0;
          uint32_t explicit_psam_cols = 0;
          uint32_t varid_semicolon = 0;
          // bugfix (28 Jun 2025): this needs to be kept in sync with
          // auto-transpose no_vmaj_ext
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "vzs", cur_modif_slen)) {
              pc.pvar_psam_flags |= kfPvarZs;
            } else if (StrStartsWith0(cur_modif, "pvar-cols=", cur_modif_slen)) {
              if (unlikely(explicit_pvar_cols)) {
                logerrputs("Error: Multiple --make-pgen pvar-cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_pvar_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[strlen("pvar-cols=")]), "xheader\0vcfheader\0maybequal\0qual\0maybefilter\0filter\0maybeinfo\0info\0maybecm\0cm\0", "make-pgen pvar-cols", kfPvarColXheader, kfPvarColDefault, 0, &pc.pvar_psam_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              if (unlikely((pc.pvar_psam_flags & kfPvarColXinfo) && (!(pc.pvar_psam_flags & (kfPvarColXheader | kfPvarColVcfheader))))) {
                logerrputs("Error: --make-pgen pvar-cols= expression cannot exclude xheader/vcfheader when\ninfo is present.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else if (StrStartsWith(cur_modif, "format=", cur_modif_slen)) {
              if (unlikely(make_plink2_flags & (kfMakePgenFormatBase * 3))) {
                logerrputs("Error: Multiple --make-pgen format= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const uint32_t fcode_minus_2 = ctou32(cur_modif[7]) - 50;
              if (unlikely((fcode_minus_2 > 2) || cur_modif[8])) {
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
              if (unlikely(make_plink2_flags & kfMakePlink2MMask)) {
                logerrputs("Error: Multiple --make-pgen multiallelics= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              const char* mode_start = (cur_modif[1] == '=')? (&(cur_modif[2])) : (&(cur_modif[14]));
              const uint32_t mode_slen = cur_modif_slen - S_CAST(uintptr_t, mode_start - cur_modif);
              if (strequal_k(mode_start, "-", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MSplitAll;
              } else if (strequal_k(mode_start, "-snps", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MSplitSnps;
              } else if (strequal_k(mode_start, "+", mode_slen) ||
                         strequal_k(mode_start, "+both", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MJoinBoth;
              } else if (strequal_k(mode_start, "+snps", mode_slen)) {
                make_plink2_flags |= kfMakePlink2MJoinSnps;
              } else if (likely(strequal_k(mode_start, "+any", mode_slen))) {
                make_plink2_flags |= kfMakePlink2MJoinAny;
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-pgen multiallelics= mode '%s'.\n", mode_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (strequal_k(cur_modif, "trim-alts", cur_modif_slen)) {
              make_plink2_flags |= kfMakePlink2TrimAlts;
            } else if (strequal_k(cur_modif, "varid-split", cur_modif_slen)) {
              varid_semicolon |= 1;
            } else if (strequal_k(cur_modif, "varid-split-dup", cur_modif_slen)) {
              varid_semicolon |= 2;
            } else if (strequal_k(cur_modif, "varid-join", cur_modif_slen)) {
              varid_semicolon |= 4;
            } else if (strequal_k(cur_modif, "varid-dup", cur_modif_slen)) {
              make_plink2_flags |= kfMakePlink2VaridDup;
            } else if (strequal_k(cur_modif, "erase-phase", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenErasePhase;
            } else if (strequal_k(cur_modif, "erase-dosage", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenEraseDosage;
            } else if (strequal_k(cur_modif, "fill-missing-from-dosage", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenFillMissingFromDosage;
            } else if (strequal_k(cur_modif, "writer-ver", cur_modif_slen)) {
              make_plink2_flags |= kfMakePgenWriterVer;
            } else if (likely(StrStartsWith0(cur_modif, "psam-cols=", cur_modif_slen))) {
              if (unlikely(explicit_psam_cols)) {
                logerrputs("Error: Multiple --make-pgen psam-cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_psam_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[strlen("psam-cols=")]), "maybefid\0fid\0maybesid\0sid\0maybeparents\0parents\0sex\0pheno1\0phenos\0", "make-pgen psam-cols", kfPsamColMaybefid, kfPsamColDefault, 0, &pc.pvar_psam_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-pgen argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (varid_semicolon) {
            if (unlikely((make_plink2_flags & kfMakePlink2VaridDup) || (varid_semicolon & (varid_semicolon - 1)))) {
              logerrputs("Error: --make-pgen 'varid-split', 'varid-split-dup', 'varid-dup', and\n'varid-join' modifiers are mutually exclusive.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (varid_semicolon & 3) {
              if (unlikely((make_plink2_flags & kfMakePlink2MJoin) || (!(make_plink2_flags & kfMakePlink2MMask)))) {
                logerrputs("Error: --make-pgen 'varid-split' must be used with a multiallelics= split mode.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else {
              if (unlikely(!(make_plink2_flags & kfMakePlink2MJoin))) {
                logerrputs("Error: --make-pgen 'varid-join' must be used with a multiallelics= join mode.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            }
            make_plink2_flags |= kfMakePlink2VaridSemicolon;
            if (varid_semicolon == 2) {
              make_plink2_flags |= kfMakePlink2VaridDup;
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
          if (unlikely(make_plink2_flags & (kfMakeBed | kfMakePgen))) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-[b]pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (likely(!strcmp(cur_modif, "zs"))) {
              make_plink2_flags |= kfMakeBimZs;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-just-bim argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          make_plink2_flags |= kfMakeBim;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-just-fam")) {
          if (unlikely(make_plink2_flags & (kfMakeBed | kfMakePgen))) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-[b]pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          make_plink2_flags |= kfMakeFam;
          pc.command_flags1 |= kfCommand1MakePlink2;
          pc.dependency_flags |= kfFilterPsamReq;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ake-just-pvar")) {
          if (unlikely(make_plink2_flags & (kfMakeBed | kfMakePgen))) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-[b]pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t explicit_cols = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.pvar_psam_flags |= kfPvarZs;
            } else if (likely(StrStartsWith0(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(explicit_cols)) {
                logerrputs("Error: Multiple --make-just-pvar cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[5]), "xheader\0vcfheader\0maybequal\0qual\0maybefilter\0filter\0maybeinfo\0info\0maybecm\0cm\0", "make-just-pvar", kfPvarColXheader, kfPvarColDefault, 0, &pc.pvar_psam_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              if (unlikely((pc.pvar_psam_flags & kfPvarColXinfo) && (!(pc.pvar_psam_flags & (kfPvarColXheader | kfPvarColVcfheader))))) {
                logerrputs("Error: --make-just-pvar cols= expression cannot exclude xheader/vcfheader when\ninfo is present.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (pc.pvar_psam_flags & kfPvarColInfo) {
                pc.dependency_flags |= kfFilterNonrefFlagsNeededSet;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-just-pvar argument '%s'.\n", cur_modif);
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
          if (unlikely(make_plink2_flags & (kfMakeBed | kfMakePgen))) {
            logerrputs("Error: --make-just-... cannot be used with --make-bed/--make-[b]pgen.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (likely(StrStartsWith0(cur_modif, "cols=", cur_modif_slen))) {
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0maybeparents\0parents\0sex\0pheno1\0phenos\0", "make-just-psam", kfPsamColMaybefid, kfPsamColDefault, 0, &pc.pvar_psam_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-just-psam argument '%s'.\n", cur_modif);
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
          if (unlikely(king_cutoff_fprefix)) {
            logerrputs("Error: --make-king cannot be used with --king-cutoff[-table] input file(s).\n");
            goto main_ret_INVALID_CMDLINE_A;
          } else if (unlikely(pc.king_table_subset_fname || pc.king_table_require_fnames)) {
            logerrputs("Error: --make-king cannot be used with --king-table-subset or\n--king-table-require[-xor].\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              if (unlikely(pc.king_flags & kfKingMatrixEncodemask)) {
                logerrputs("Error: Multiple --make-king encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixZs;
            } else if (strequal_k(cur_modif, "bin", cur_modif_slen)) {
              if (unlikely(pc.king_flags & kfKingMatrixEncodemask)) {
                logerrputs("Error: Multiple --make-king encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixBin;
            } else if (strequal_k(cur_modif, "bin4", cur_modif_slen)) {
              if (unlikely(pc.king_flags & kfKingMatrixEncodemask)) {
                logerrputs("Error: Multiple --make-king encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixBin4;
            } else if (strequal_k(cur_modif, "square", cur_modif_slen)) {
              if (unlikely(pc.king_flags & kfKingMatrixShapemask)) {
                logerrputs("Error: Multiple --make-king shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixSq;
            } else if (strequal_k(cur_modif, "square0", cur_modif_slen)) {
              if (unlikely(pc.king_flags & kfKingMatrixShapemask)) {
                logerrputs("Error: Multiple --make-king shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixSq0;
            } else if (likely(strequal_k(cur_modif, "triangle", cur_modif_slen))) {
              if (unlikely(pc.king_flags & kfKingMatrixShapemask)) {
                logerrputs("Error: Multiple --make-king shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.king_flags |= kfKingMatrixTri;
            } else if (strequal_k(cur_modif, "no-idheader", cur_modif_slen)) {
              logerrputs("Error: --make-king 'no-idheader' modifier retired.  Use --no-id-header instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-king argument '%s'.\n", cur_modif);
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
          if (unlikely(king_cutoff_fprefix)) {
            logerrputs("Error: --make-king-table cannot be used with --king-cutoff[-table] input file(s).\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.king_flags |= kfKingTableZs;
            } else if (strequal_k(cur_modif, "counts", cur_modif_slen)) {
              pc.king_flags |= kfKingCounts;
            } else if (strequal_k(cur_modif, "rel-check", cur_modif_slen)) {
              pc.king_flags |= kfKingRelCheck;
            } else if (strequal_k(cur_modif, "concordance-check", cur_modif_slen)) {
              pc.king_flags |= kfKingConcordanceCheck;
            } else if (likely(StrStartsWith(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(pc.king_flags & kfKingColAll)) {
                logerrputs("Error: Multiple --make-king-table cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0id\0maybesid\0sid\0nsnp\0hethet\0ibs0\0ibs1\0ibs\0kinship\0", "make-king-table", kfKingColMaybefid, kfKingColDefault, 1, &pc.king_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              if (!(pc.king_flags & kfKingColId)) {
                if (unlikely(pc.king_flags & (kfKingColMaybefid | kfKingColFid | kfKingColMaybesid | kfKingColSid))) {
                  logerrputs("Error: Invalid --make-king-table column set descriptor ('maybefid', 'fid',\n'maybesid', and 'sid' require 'id').\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
                if (unlikely(pc.king_table_filter != -DBL_MAX)) {
                  logerrputs("Error: --king-table-filter requires --make-king-table cols= to include the 'id'\ncolumn set.\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
                if (unlikely(pc.king_table_subset_fname)) {
                  logerrputs("Error: --king-table-subset requires --make-king-table cols= to include the 'id'\ncolumn set.\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
                if (unlikely(pc.king_table_require_fnames)) {
                  logerrputs("Error: --king-table-require[-xor] requires --make-king-table cols= to include\nthe 'id' column set.\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
              }
            } else if (strequal_k(cur_modif, "no-idheader", cur_modif_slen)) {
              logerrputs("Error: --make-king-table 'no-idheader' modifier retired.  Use --no-id-header\ninstead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-king-table argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const uint32_t rel_concordance_flags = pc.king_flags & (kfKingRelCheck | kfKingConcordanceCheck);
          if (rel_concordance_flags) {
            if (unlikely(rel_concordance_flags == (kfKingRelCheck | kfKingConcordanceCheck))) {
              logerrputs("Error: --make-king-table 'concordance-check' and 'rel-check' modifiers cannot\nbe used together.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            if (unlikely(pc.command_flags1 & kfCommand1KingCutoff)) {
              logerrputs("Error: --make-king-table 'concordance-check'/'rel-check' cannot be used with\n--king-cutoff.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          }
          if (!(pc.king_flags & kfKingColAll)) {
            pc.king_flags |= kfKingColDefault;
          }
          pc.command_flags1 |= kfCommand1MakeKing;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "issing")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.missing_rpt_flags |= kfMissingRptZs;
            } else if (strequal_k(cur_modif, "sample-only", cur_modif_slen)) {
              if (unlikely(pc.missing_rpt_flags & kfMissingRptVariantOnly)) {
                logerrputs("Error: --missing 'sample-only' and 'variant-only' cannot be used together.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              pc.missing_rpt_flags |= kfMissingRptSampleOnly;
            } else if (strequal_k(cur_modif, "variant-only", cur_modif_slen)) {
              if (unlikely(pc.missing_rpt_flags & kfMissingRptSampleOnly)) {
                logerrputs("Error: --missing 'sample-only' and 'variant-only' cannot be used together.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              pc.missing_rpt_flags |= kfMissingRptVariantOnly;
            } else if (StrStartsWith(cur_modif, "scols=", cur_modif_slen)) {
              if (unlikely(pc.missing_rpt_flags & kfMissingRptScolAll)) {
                logerrputs("Error: Multiple --missing scols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("scols=")]), "maybefid\0fid\0maybesid\0sid\0misspheno1\0missphenos\0nmissdosage\0nmiss\0nmisshh\0hethap\0nobs\0fmissdosage\0fmiss\0fmisshh\0", "missing scols", kfMissingRptScolMaybefid, kfMissingRptScolDefault, 1, &pc.missing_rpt_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (likely(StrStartsWith(cur_modif, "vcols=", cur_modif_slen))) {
              if (unlikely(pc.missing_rpt_flags & kfMissingRptVcolAll)) {
                logerrputs("Error: Multiple --missing vcols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("vcols=")]), "chrom\0pos\0ref\0alt\0maybeprovref\0provref\0nmissdosage\0nmiss\0nmisshh\0hethap\0nobs\0fmissdosage\0fmiss\0fmisshh\0fhethap\0", "missing vcols", kfMissingRptVcolChrom, kfMissingRptVcolDefault, 1, &pc.missing_rpt_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --missing argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const uint32_t explicit_scols = pc.missing_rpt_flags & kfMissingRptScolAll;
          if (pc.missing_rpt_flags & kfMissingRptVariantOnly) {
            if (unlikely(explicit_scols)) {
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
            if (unlikely(explicit_vcols)) {
              logerrputs("Error: --missing 'sample-only' and 'vcols=' modifiers cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else if (!explicit_vcols) {
            pc.missing_rpt_flags |= kfMissingRptVcolDefault;
          }
          pc.command_flags1 |= kfCommand1MissingReport;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "aj-ref")) {
          if (unlikely(pc.alt_allele_flag)) {
            logerrputs("Error: --maj-ref cannot be used with --ref-allele/--alt[1]-allele.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (likely(!strcmp(cur_modif, "force"))) {
              pc.misc_flags |= kfMiscMajRefForce;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --maj-ref argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.misc_flags |= kfMiscMajRef;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "af")) {
          if (unlikely(pc.flip_info.subset_fname || (make_plink2_flags & kfMakePlink2FillMissingWithRef))) {
            logerrputs("Error: --maf cannot be used with --fill-missing-with-ref/--flip-subset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            const char* mode_str = ScantokDouble(cur_modif, &pc.min_maf);
            if (!mode_str) {
              pc.min_maf = 0.01;
              if (unlikely(param_ct == 2)) {
                logerrputs("Error: Invalid --maf argument sequence.\n");
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              mode_str = cur_modif;
            } else {
              if (unlikely(pc.min_maf < 0.0)) {
                snprintf(g_logbuf, kLogbufSize, "Error: --maf argument '%s' too small (must be >= 0).\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              } else if (unlikely(pc.min_maf > 1.0)) {
                snprintf(g_logbuf, kLogbufSize, "Error: --maf argument '%s' too large (must be <= 1).\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (mode_str[0] == ':') {
                if (unlikely(param_ct == 2)) {
                  logerrputs("Error: Invalid --maf argument sequence.\n");
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
              } else {
                if (unlikely(mode_str[0])) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --maf argument '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                if (param_ct == 2) {
                  mode_str = argvk[arg_idx + 2];
                } else {
                  mode_str = nullptr;
                }
              }
            }
            if (mode_str) {
              if (ParseFreqSelector(mode_str, flagname_p, &(pc.filter_modes[0]))) {
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            }
          } else {
            pc.min_maf = 0.01;
          }
          if (pc.min_maf != 0.0) {
            pc.filter_flags |= kfFilterPvarReq;
            pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ax-maf")) {
          if (unlikely(pc.flip_info.subset_fname || (make_plink2_flags & kfMakePlink2FillMissingWithRef))) {
            logerrputs("Error: --max-maf cannot be used with --fill-missing-with-ref/--flip-subset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const char* mode_str = ScantokDouble(cur_modif, &pc.max_maf);
          if (unlikely(!mode_str)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-maf argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(pc.max_maf >= 1.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: --max-maf argument '%s' too large (must be < 1).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (mode_str[0] == ':') {
            if (unlikely(param_ct == 2)) {
              logerrputs("Error: Invalid --max-maf argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            if (unlikely(mode_str[0])) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-maf argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (param_ct == 2) {
              mode_str = argvk[arg_idx + 2];
            } else {
              mode_str = nullptr;
            }
          }
          if (mode_str) {
            if (ParseFreqSelector(mode_str, flagname_p, &(pc.filter_modes[1]))) {
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (unlikely((pc.filter_modes[0] == pc.filter_modes[1]) && (pc.max_maf < pc.min_maf))) {
            snprintf(g_logbuf, kLogbufSize, "Error: --max-maf argument '%s' too small (must be >= %g).\n", cur_modif, pc.min_maf);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "ac")) {
          if (unlikely(pc.flip_info.subset_fname || (make_plink2_flags & kfMakePlink2FillMissingWithRef))) {
            logerrputs("Error: --mac cannot be used with --fill-missing-with-ref/--flip-subset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          const char* mode_str = ScantokDouble(cur_modif, &dxx);
          if (unlikely((!mode_str) || (dxx < 0.0) || (dxx > 2147483646.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mac argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE;
          }
          if (dxx > 0.0) {
            // round up, but keep as much precision as possible
            int32_t int_part = S_CAST(int32_t, dxx);
            dxx -= int_part;
            pc.min_allele_ddosage = int_part * S_CAST(uint64_t, kDosageMax);
            if (dxx > 0.0) {
              pc.min_allele_ddosage += 1 + S_CAST(uint64_t, dxx * (kDosageMax * (1 - kSmallEpsilon)));
            }
            // yeah, this should be its own function...
            if (mode_str[0] == ':') {
              if (unlikely(param_ct == 2)) {
                logerrputs("Error: Invalid --mac argument sequence.\n");
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else {
              if (unlikely(mode_str[0])) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mac argument '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (param_ct == 2) {
                mode_str = argvk[arg_idx + 2];
              } else {
                mode_str = nullptr;
              }
            }
            if (mode_str) {
              if (ParseFreqSelector(mode_str, flagname_p, &(pc.filter_modes[2]))) {
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            }
            pc.filter_flags |= kfFilterPvarReq;
            pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ax-mac")) {
          if (unlikely(pc.flip_info.subset_fname || (make_plink2_flags & kfMakePlink2FillMissingWithRef))) {
            logerrputs("Error: --max-mac cannot be used with --fill-missing-with-ref/--flip-subset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          const char* mode_str = ScantokDouble(cur_modif, &dxx);
          if (unlikely((!mode_str) || (dxx < 0.0) || (dxx > 2147483646.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-mac argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          // round down
          pc.max_allele_ddosage = S_CAST(int64_t, dxx * kDosageMax);
          if (mode_str[0] == ':') {
            if (unlikely(param_ct == 2)) {
              logerrputs("Error: Invalid --max-mac argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            if (unlikely(mode_str[0])) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-mac argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (param_ct == 2) {
              mode_str = argvk[arg_idx + 2];
            } else {
              mode_str = nullptr;
            }
          }
          if (mode_str) {
            if (ParseFreqSelector(mode_str, flagname_p, &(pc.filter_modes[3]))) {
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (unlikely((pc.filter_modes[2] == pc.filter_modes[3]) && (pc.max_allele_ddosage < pc.min_allele_ddosage))) {
            // yeah, --mac 0.1 --max-mac 0.1 also isn't allowed
            logerrputs("Error: --max-mac argument cannot be smaller than --mac argument when modes are\nidentical.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "ind")) {
          if (unlikely(make_plink2_flags & kfMakePlink2FillMissingWithRef)) {
            logerrputs("Error: --geno/--mind cannot be used with --fill-missing-with-ref.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t mind_thresh_present = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "dosage", cur_modif_slen)) {
              pc.misc_flags |= kfMiscMindDosage;
            } else if (strequal_k(cur_modif, "hh-missing", cur_modif_slen)) {
              pc.misc_flags |= kfMiscMindHhMissing;
            } else if (unlikely(mind_thresh_present)) {
              logerrputs("Error: Invalid --mind argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (unlikely(!ScantokDouble(cur_modif, &pc.mind_thresh))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mind argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else if (unlikely((pc.mind_thresh < 0.0) || (pc.mind_thresh > 1.0))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mind argument '%s' (must be in [0, 1]).\n", cur_modif);
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
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.missing_varid_match);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "erge-par")) {
          if (pc.exportf_info.flags & kfExportfVcf) {
            logerrputs("Warning: --merge-par should not be used with VCF export.  (The VCF export\nroutine automatically converts PAR1/PAR2 chromosome codes to X, while using\nthe PAR boundaries to get male ploidy right; --merge-par causes VCF export to\nget male ploidy wrong.)\n");
          }
          pc.misc_flags |= kfMiscMergePar;
          pc.dependency_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "erge-x")) {
          if (unlikely(pc.misc_flags & kfMiscMergePar)) {
            logerrputs("Error: --merge-par cannot be used with --merge-x.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (pc.exportf_info.flags & kfExportfVcf) {
            logerrputs("Warning: --merge-x should not be used in the same run as VCF export; this\ncauses some ploidies to be wrong.  Instead, use --merge-x + --sort-vars +\n--make-[b]pgen in one run, and follow up with --split-par + --export vcf.\n");
          }
          pc.misc_flags |= kfMiscMergeX;
          pc.dependency_flags |= kfFilterPvarReq;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "af-succ")) {
          if (pc.af_pseudocount != 0.0) {
            logerrputs("Error: --maf-succ cannot be used with --af-pseudocount.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          logputs("Note: --maf-succ flag deprecated.  Use \"--af-pseudocount 1\" instead.\n");
          pc.af_pseudocount = 1.0;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ax-corr")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Glm))) {
            logerrputs("Error: --max-corr must be used with --glm.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(!ScantokDouble(cur_modif, &pc.glm_info.max_corr))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-corr argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely((pc.glm_info.max_corr < 0.0) || (pc.glm_info.max_corr > 1.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-corr argument '%s' (must be in [0, 1]).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ach-r2-filter")) {
          if (unlikely(pc.freq_rpt_flags & kfAlleleFreqColMinimac3R2)) {
            logerrputs("Error: --freq minimac3r2 output and --mach-r2-filter cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.flip_info.subset_fname || (make_plink2_flags & kfMakePlink2FillMissingWithRef))) {
            logerrputs("Error: --mach-r2-filter cannot be used with\n--fill-missing-with-ref/--flip-subset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (unlikely(!ScantokDouble(cur_modif, &pc.mach_r2_min))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mach-r2-filter min argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (unlikely(pc.mach_r2_min < 0.0)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mach-r2-filter min argument '%s' (must be nonnegative).\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (param_ct == 2) {
              cur_modif = argvk[arg_idx + 2];
              if (unlikely((!ScantokDouble(cur_modif, &pc.mach_r2_max)) || (pc.mach_r2_max == 0.0))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mach-r2-filter max argument '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else {
              pc.mach_r2_max = 2.0;
            }
            if (unlikely(pc.mach_r2_max < pc.mach_r2_min)) {
              logerrputs("Error: --mach-r2-filter min argument cannot be larger than max argument.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else {
            pc.mach_r2_min = 0.1;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "inimac3-r2-filter")) {
          if (unlikely(pc.mach_r2_max != 0.0)) {
            logerrputs("Error: --mach-r2-filter and --minimac3-r2-filter cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.freq_rpt_flags & kfAlleleFreqColMachR2)) {
            logerrputs("Error: --freq machr2 output and --minimac3-r2-filter cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.flip_info.subset_fname || (make_plink2_flags & kfMakePlink2FillMissingWithRef))) {
            logerrputs("Error: --minimac3-r2-filter cannot be used with\n--fill-missing-with-ref/--flip-subset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(!ScantokDouble(cur_modif, &pc.minimac3_r2_min))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --minimac3-r2-filter min argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(pc.minimac3_r2_min < 0.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --minimac3-r2-filter min argument '%s' (must be nonnegative).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (param_ct == 2) {
            cur_modif = argvk[arg_idx + 2];
            if (unlikely((!ScantokDouble(cur_modif, &pc.minimac3_r2_max)) || (pc.minimac3_r2_max == 0.0))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --minimac3-r2-filter max argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.minimac3_r2_max = 1.0;
          }
          if (unlikely(pc.minimac3_r2_max < pc.minimac3_r2_min)) {
            logerrputs("Error: --minimac3-r2-filter min argument cannot be larger than max argument.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.dependency_flags |= kfFilterAllReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "issing-code")) {
          if (unlikely(!(xload & (kfXloadOxGen | kfXloadOxBgen | kfXloadOxHaps)))) {
            // could technically support pure .sample -> .fam/.psam, but let's
            // keep this simple
            logerrputs("Error: --missing-code must be used with --data/--gen/--bgen/--haps.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(param_ct? argvk[arg_idx + 1] : "", argvk[arg_idx], 0x7fffffff, &ox_missing_code);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "issing-genotype"))) {
          logerrputs("Error: --missing-genotype flag retired.  Use --input-missing-genotype and/or\n--output-missing-genotype.\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "issing-phenotype"))) {
          logerrputs("Error: --missing-phenotype flag retired.  Use --input-missing-phenotype and/or\n--output-missing-phenotype.\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (strequal_k_unsafe(flagname_p2, "issing-catname")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          double dxx;
          if (unlikely(ScanadvDouble(cur_modif, &dxx) || IsNanStr(cur_modif, cur_modif_slen))) {
            logerrputs("Error: --missing-catname string cannot be 'NA' or start with a number.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(cur_modif_slen > 31)) {
            logerrputs("Error: --missing-catname string too long (max 31 chars).\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          memcpy(pc.missing_catname, cur_modif, cur_modif_slen + 1);
        } else if (strequal_k_unsafe(flagname_p2, "ouse")) {
          if (unlikely(chr_info.chrset_source)) {
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
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "ake-grm"))) {
          logerrputs("Error: --make-grm has been retired due to inconsistent meaning across GCTA\nversions.  Use --make-grm-list or --make-grm-bin.\n");
          goto main_ret_INVALID_CMDLINE;
        } else if (strequal_k_unsafe(flagname_p2, "ake-grm-bin")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          pc.grm_flags |= kfGrmNoIdHeader | kfGrmBin;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "cov", cur_modif_slen)) {
              pc.grm_flags |= kfGrmCov;
            } else if (strequal_k(cur_modif, "meanimpute", cur_modif_slen)) {
              pc.grm_flags |= kfGrmMeanimpute;
            } else if (strequal_k(cur_modif, "id-header", cur_modif_slen) ||
                       strequal_k(cur_modif, "idheader", cur_modif_slen)) {
              pc.grm_flags &= ~kfGrmNoIdHeader;
            } else if (likely(strequal_k(cur_modif, "iid-only", cur_modif_slen))) {
              pc.grm_flags |= kfGrmNoIdHeaderIidOnly;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-grm-bin argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (unlikely((pc.grm_flags & (kfGrmNoIdHeader | kfGrmNoIdHeaderIidOnly)) == kfGrmNoIdHeaderIidOnly)) {
            logerrputs("Error: --make-grm-bin 'id-header' and 'iid-only' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.command_flags1 |= kfCommand1MakeRel;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-grm-gz") || strequal_k_unsafe(flagname_p2, "ake-grm-list")) {
          // While --make-grm-gz is otherwise retired, it's painless to accept
          // and translate "--make-grm-gz no-gz".
          if (unlikely(pc.command_flags1 & kfCommand1MakeRel)) {
            if (pc.grm_flags & kfGrmBin) {
              logerrputs("Error: --make-grm-list cannot be used with --make-grm-bin.\n");
            } else {
              logerrputs("Error: --make-grm-list cannot be used with --make-grm-gz.\n");
            }
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t no_gz = 0;
          pc.grm_flags |= kfGrmNoIdHeader | kfGrmList;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "cov", cur_modif_slen)) {
              pc.grm_flags |= kfGrmCov;
            } else if (strequal_k(cur_modif, "meanimpute", cur_modif_slen)) {
              pc.grm_flags |= kfGrmMeanimpute;
            } else if (strequal_k(cur_modif, "no-gz", cur_modif_slen)) {
              no_gz = 1;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.grm_flags |= kfGrmListZs;
            } else if (strequal_k(cur_modif, "id-header", cur_modif_slen) ||
                       strequal_k(cur_modif, "idheader", cur_modif_slen)) {
              pc.grm_flags &= ~kfGrmNoIdHeader;
            } else if (likely(strequal_k(cur_modif, "iid-only", cur_modif_slen))) {
              pc.grm_flags |= kfGrmNoIdHeaderIidOnly;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-grm-list argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (flagname_p2[8] == 'g') {
            if (unlikely(!no_gz)) {
              // screw it, life is too much better with multithreaded .zst
              logerrputs("Error: --make-grm-list no longer supports gzipped output.  Use 'zs' for\nzstd-compressed output (much faster), or use PLINK 1.9 for this function.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            logerrputs("Warning: --make-grm-gz has been renamed to --make-grm-list.\n");
          }
          if (unlikely((pc.grm_flags & (kfGrmNoIdHeader | kfGrmNoIdHeaderIidOnly)) == kfGrmNoIdHeaderIidOnly)) {
            logerrputs("Error: --make-grm-list 'id-header' and 'iid-only' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.command_flags1 |= kfCommand1MakeRel;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-grm-sparse")) {
          if (unlikely(pc.command_flags1 & kfCommand1MakeRel)) {
            // easy to support, but I don't know why anyone would want it since
            // GCTA is the main sparse-GRM consumer and it has --make-bK-sparse
            // to subset.
            if (pc.grm_flags & kfGrmBin) {
              logerrputs("Error: --make-grm-sparse cannot be used with --make-grm-bin.\n");
            } else {
              logerrputs("Error: --make-grm-sparse cannot be used with --make-grm-list.\n");
            }
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 5))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          {
            const char* cur_modif = argvk[arg_idx + 1];
            double dxx;
            if (unlikely(!ScantokDouble(cur_modif, &dxx))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-grm-sparse threshold '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            pc.grm_sparse_cutoff = dxx * (1 - kSmallEpsilon);
          }
          pc.grm_flags |= kfGrmNoIdHeader | kfGrmSparse;
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "cov", cur_modif_slen)) {
              pc.grm_flags |= kfGrmCov;
            } else if (strequal_k(cur_modif, "meanimpute", cur_modif_slen)) {
              pc.grm_flags |= kfGrmMeanimpute;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.grm_flags |= kfGrmListZs;
            } else if (strequal_k(cur_modif, "id-header", cur_modif_slen) ||
                       strequal_k(cur_modif, "idheader", cur_modif_slen)) {
              pc.grm_flags &= ~kfGrmNoIdHeader;
            } else if (likely(strequal_k(cur_modif, "iid-only", cur_modif_slen))) {
              pc.grm_flags |= kfGrmNoIdHeaderIidOnly;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-grm-sparse argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (unlikely((pc.grm_flags & (kfGrmNoIdHeader | kfGrmNoIdHeaderIidOnly)) == kfGrmNoIdHeaderIidOnly)) {
            logerrputs("Error: --make-grm-sparse 'id-header' and 'iid-only' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.command_flags1 |= kfCommand1MakeRel;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ake-rel")) {
          if (unlikely(pc.command_flags1 & kfCommand1MakeRel)) {
            logerrputs("Error: --make-rel cannot be used with --make-grm-{bin,list,sparse}.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "cov", cur_modif_slen)) {
              pc.grm_flags |= kfGrmCov;
            } else if (strequal_k(cur_modif, "meanimpute", cur_modif_slen)) {
              pc.grm_flags |= kfGrmMeanimpute;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              if (unlikely(pc.grm_flags & kfGrmMatrixEncodemask)) {
                logerrputs("Error: Multiple --make-rel encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixZs;
            } else if (unlikely(strequal_k(cur_modif, "no-idheader", cur_modif_slen))) {
              logerrputs("Error: --make-rel 'no-idheader' modifier retired.  Use --no-id-header instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (strequal_k(cur_modif, "bin", cur_modif_slen)) {
              if (unlikely(pc.grm_flags & kfGrmMatrixEncodemask)) {
                logerrputs("Error: Multiple --make-rel encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixBin;
            } else if (strequal_k(cur_modif, "bin4", cur_modif_slen)) {
              if (unlikely(pc.grm_flags & kfGrmMatrixEncodemask)) {
                logerrputs("Error: Multiple --make-rel encoding modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixBin4;
            } else if (strequal_k(cur_modif, "square", cur_modif_slen)) {
              if (unlikely(pc.grm_flags & kfGrmMatrixShapemask)) {
                logerrputs("Error: Multiple --make-rel shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixSq;
            } else if (strequal_k(cur_modif, "square0", cur_modif_slen)) {
              if (unlikely(pc.grm_flags & kfGrmMatrixShapemask)) {
                logerrputs("Error: Multiple --make-rel shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixSq0;
            } else if (likely(strequal_k(cur_modif, "triangle", cur_modif_slen))) {
              if (unlikely(pc.grm_flags & kfGrmMatrixShapemask)) {
                logerrputs("Error: Multiple --make-rel shape modifiers.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              pc.grm_flags |= kfGrmMatrixTri;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-rel argument '%s'.\n", cur_modif);
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
          if (unlikely(load_params || (xload & (~kfXloadPlink1Dosage)))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --map filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pvarname, fname, slen + 1);
          xload |= kfXloadMap;
        } else if (strequal_k_unsafe(flagname_p2, "within")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(ScanPosintCappedx(cur_modif, kMaxLongLine / 2, &pc.mwithin_val))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mwithin argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "filter")) {
          if (unlikely(!pc.keep_col_match_fname)) {
            logerrputs("Error: --mfilter must be used with --keep-col-match.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.keep_col_match_name)) {
            logerrputs("Error: --mfilter can't be used with --keep-col-match-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.keep_col_match_num)) {
            logerrputs("Error: --mfilter can't be used with --keep-col-match-num.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          logerrputs("Warning: --mfilter flag deprecated.  Use --keep-col-match-num or\n--keep-col-match-name instead.  (Note that --keep-col-match-num does not add 2\nto the column number.)\n");
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t mfilter_arg;
          if (unlikely(ScanPosintDefcapx(cur_modif, &mfilter_arg))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mfilter argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.keep_col_match_num = mfilter_arg + 2;
        } else if (strequal_k_unsafe(flagname_p2, "ax-alleles")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(ScanPosintDefcapx(cur_modif, &pc.filter_max_allele_ct))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --max-alleles argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogMaxAlleles;
        } else if (strequal_k_unsafe(flagname_p2, "in-alleles")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(ScanPosintDefcapx(cur_modif, &pc.filter_min_allele_ct))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --min-alleles argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(pc.filter_min_allele_ct > pc.filter_max_allele_ct)) {
            logerrputs("Error: --min-alleles argument can't be larger than --max-alleles argument.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (pc.filter_min_allele_ct == 1) {
            pc.filter_min_allele_ct = 0;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogMinAlleles;
        } else if (strequal_k_unsafe(flagname_p2, "erge-mode")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          if (strequal_k(cur_modif, "nm-first", cur_modif_slen) || strequal_k(cur_modif, "2", cur_modif_slen)) {
            pmerge_info.merge_mode = kMergeModeNmFirst;
          } else if (strequal_k(cur_modif, "first", cur_modif_slen) || strequal_k(cur_modif, "4", cur_modif_slen)) {
            pmerge_info.merge_mode = kMergeModeFirst;
          } else if (unlikely(strequal_k(cur_modif, "3", cur_modif_slen))) {
            logerrputs("Error: --merge-mode 3 discontinued.  (You can get the same results with\n--merge-mode 2 if you reverse your fileset order.)\n");
            goto main_ret_INVALID_CMDLINE_A;
          } else if (unlikely(strequal_k(cur_modif, "5", cur_modif_slen))) {
            logerrputs("Error: --merge-mode 5 discontinued.  (You can get the same results with\n--merge-mode 4 if you reverse your fileset order.)\n");
            goto main_ret_INVALID_CMDLINE_A;
          } else if (unlikely(strequal_k(cur_modif, "6", cur_modif_slen) || strequal_k(cur_modif, "7", cur_modif_slen))) {
            logerrputs("Error: --merge-mode 6/7 has been discontinued.  Use --pgen-diff instead.\n");
            goto main_ret_INVALID_CMDLINE_A;
          } else if (unlikely(!(strequal_k(cur_modif, "nm-match", cur_modif_slen) || strequal_k(cur_modif, "1", cur_modif_slen)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --merge-mode argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pmerge_required = 1;
        } else if (strequal_k_unsafe(flagname_p2, "erge-parents-mode") ||
                   strequal_k_unsafe(flagname_p2, "erge-sex-mode") ||
                   strequal_k_unsafe(flagname_p2, "erge-pheno-mode")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          MergePhenoMode mode;
          if (strequal_k(cur_modif, "nm-match", cur_modif_slen) || strequal_k(cur_modif, "1", cur_modif_slen)) {
            mode = kMergePhenoModeNmMatch;
          } else if (strequal_k(cur_modif, "nm-first", cur_modif_slen) || strequal_k(cur_modif, "2", cur_modif_slen)) {
            mode = kMergePhenoModeNmFirst;
          } else if (likely(strequal_k(cur_modif, "first", cur_modif_slen) || strequal_k(cur_modif, "4", cur_modif_slen))) {
            mode = kMergePhenoModeFirst;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s argument '%s'.\n", flagname_p, cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pmerge_required = 1;
          if (flagname_p2[5] == 'p') {
            if (flagname_p2[6] == 'a') {
              pmerge_info.merge_parents_mode = mode;
            } else {
              pmerge_info.merge_pheno_mode = mode;
            }
          } else {
            pmerge_info.merge_sex_mode = mode;
          }
        } else if (strequal_k_unsafe(flagname_p2, "erge-xheader-mode")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          if (strequal_k(cur_modif, "erase", cur_modif_slen)) {
            if (unlikely(pmerge_info.merge_info_mode != kMergeInfoCmModeErase)) {
              logerrputs("Error: \"--merge-xheader-mode erase\" requires \"--merge-info-mode erase\".\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pmerge_info.merge_xheader_mode = kMergeXheaderModeErase;
          } else if (strequal_k(cur_modif, "match", cur_modif_slen)) {
            pmerge_info.merge_xheader_mode = kMergeXheaderModeMatch;
          } else if (unlikely(!strequal_k(cur_modif, "first", cur_modif_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --merge-xheader-mode argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pmerge_required = 1;
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "erge-equal-pos"))) {
          logerrputs("Error: --merge-equal-pos has been retired.  Use e.g. --set-all-var-ids before\nmerging instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "erge-qual-mode")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          if (strequal_k(cur_modif, "erase", cur_modif_slen)) {
            pmerge_info.merge_qual_mode = kMergeQualModeErase;
          } else if (strequal_k(cur_modif, "nm-match", cur_modif_slen)) {
            pmerge_info.merge_qual_mode = kMergeQualModeNmMatch;
          } else if (strequal_k(cur_modif, "nm-first", cur_modif_slen)) {
            pmerge_info.merge_qual_mode = kMergeQualModeNmFirst;
          } else if (strequal_k(cur_modif, "first", cur_modif_slen)) {
            pmerge_info.merge_qual_mode = kMergeQualModeFirst;
          } else if (unlikely(!strequal_k(cur_modif, "min", cur_modif_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --merge-qual-mode argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pmerge_required = 1;
        } else if (strequal_k_unsafe(flagname_p2, "erge-filter-mode")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          if (strequal_k(cur_modif, "erase", cur_modif_slen)) {
            pmerge_info.merge_filter_mode = kMergeFilterModeErase;
          } else if (strequal_k(cur_modif, "nm-match", cur_modif_slen)) {
            pmerge_info.merge_filter_mode = kMergeFilterModeNmMatch;
          } else if (strequal_k(cur_modif, "nm-first", cur_modif_slen)) {
            pmerge_info.merge_filter_mode = kMergeFilterModeNmFirst;
          } else if (strequal_k(cur_modif, "first", cur_modif_slen)) {
            pmerge_info.merge_filter_mode = kMergeFilterModeFirst;
          } else if (unlikely(!strequal_k(cur_modif, "np-union", cur_modif_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --merge-filter-mode argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pmerge_required = 1;
        } else if (strequal_k_unsafe(flagname_p2, "erge-info-mode") ||
                   strequal_k_unsafe(flagname_p2, "erge-cm-mode")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          MergeInfoCmMode mode;
          if (strequal_k(cur_modif, "erase", cur_modif_slen)) {
            mode = kMergeInfoCmModeErase;
          } else if (strequal_k(cur_modif, "nm-match", cur_modif_slen)) {
            mode = kMergeInfoCmModeNmMatch;
          } else if (strequal_k(cur_modif, "nm-first", cur_modif_slen)) {
            mode = kMergeInfoCmModeNmFirst;
          } else if (likely(strequal_k(cur_modif, "first", cur_modif_slen))) {
            mode = kMergeInfoCmModeFirst;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s argument '%s'.\n", flagname_p, cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pmerge_required = 1;
          if (flagname_p2[5] == 'i') {
            pmerge_info.merge_info_mode = mode;
          } else {
            pmerge_info.merge_cm_mode = mode;
          }
        } else if (strequal_k_unsafe(flagname_p2, "erge-info-sort") ||
                   strequal_k_unsafe(flagname_p2, "erge-pheno-sort")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* mode_str = argvk[arg_idx + 1];
          const char first_char_upcase_match = mode_str[0] & 0xdf;
          const uint32_t mode_slen = strlen(mode_str);
          SortMode sort_mode = kSortNone;
          if (((mode_slen == 1) && (first_char_upcase_match == 'N')) ||
              strequal_k(mode_str, "natural", mode_slen)) {
            sort_mode = kSortNatural;
          } else if (((mode_slen == 1) && (first_char_upcase_match == 'A')) ||
                     strequal_k(mode_str, "ascii", mode_slen)) {
            sort_mode = kSortAscii;
          } else if (unlikely(!(strequal_k(mode_str, "0", mode_slen) ||
                                strequal_k(mode_str, "none", mode_slen)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s argument '%s'.\n", flagname_p, mode_str);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pmerge_required = 1;
          if (flagname_p2[5] == 'i') {
            pmerge_info.merge_info_sort = sort_mode;
          } else {
            pmerge_info.merge_pheno_sort = sort_mode;
          }
        } else if (strequal_k_unsafe(flagname_p2, "erge-max-alleles")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          // Limit to .pgen maximum - 1, to make it straightforward to
          // represent over-the-limit variants in temporary files.
          uint32_t merge_max_allele_ct;
          if (unlikely(ScanPosintCappedx(cur_modif, kPglMaxAlleleCt - 1, &merge_max_allele_ct) || (merge_max_allele_ct == 1))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --merge-max-alleles argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pmerge_info.max_allele_ct = merge_max_allele_ct;
          pc.load_filter_log_flags |= kfLoadFilterLogMergeMaxAlleles;
        } else if (strequal_k_unsafe(flagname_p2, "ultiallelics-already-joined")) {
          pmerge_info.flags |= kfPmergeMultiallelicsAlreadyJoined;
          pmerge_required = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ake-founders")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t is_first = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "require-2-missing", cur_modif_slen)) {
              pc.misc_flags |= kfMiscMakeFoundersRequire2Missing;
            } else if (likely(strequal_k(cur_modif, "first", cur_modif_slen))) {
              is_first = 1;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --make-founders argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.misc_flags |= is_first? kfMiscMakeFoundersFirst : kfMiscMakeFoundersNotfirst;
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "e")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t numeric_param_ct = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (!strcmp(cur_modif, "var-first")) {
              pc.mendel_info.flags |= kfMendelFilterVarFirst;
            } else {
              double dxx;
              if (unlikely((numeric_param_ct == 2) || (!ScantokDouble(cur_modif, &dxx)))) {
                logerrputs("Error: Invalid --me argument sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely((dxx < 0.0) || (dxx > 1.0))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --me max per-%s error '%s'.\n", numeric_param_ct? "variant" : "trio", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (numeric_param_ct == 0) {
                pc.mendel_info.max_trio_error = dxx;
              } else {
                pc.mendel_info.max_var_error = dxx;
              }
              ++numeric_param_ct;
            }
          }
          if (unlikely(numeric_param_ct < 2)) {
            logerrputs("Error: --me requires max per-trio and per-variant error rates.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if ((pc.mendel_info.max_trio_error == 1.0) && (pc.mendel_info.max_var_error == 1.0)) {
            logputs("\"--me 1 1\" has no effect; ignoring.\n");
            pc.mendel_info.flags &= ~kfMendelFilterVarFirst;
          } else {
            pc.filter_flags |= kfFilterAllReq | kfFilterMendel;
          }
        } else if (strequal_k_unsafe(flagname_p2, "e-exclude-one")) {
          if (unlikely(!(pc.filter_flags & kfFilterMendel))) {
            logerrputs("Error: --me-exclude-one must be used with --me.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          pc.mendel_info.exclude_one_ratio = -1;
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            double dxx;
            if (unlikely(!ScantokDouble(cur_modif, &dxx))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --me-exclude-one ratio '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (unlikely(dxx <= 1.0)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --me-exclude-one ratio '%s' (must be larger than 1).\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            pc.mendel_info.exclude_one_ratio = dxx;
          }
        } else if (strequal_k_unsafe(flagname_p2, "endel")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t explicit_cols = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.mendel_info.flags |= kfMendelRptZs;
            } else if (strequal_k(cur_modif, "summaries-only", cur_modif_slen)) {
              pc.mendel_info.flags |= kfMendelRptSummariesOnly;
            } else if (likely(StrStartsWith(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(explicit_cols)) {
                logerrputs("Error: Multiple --mendel cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("cols=")]), "maybefid\0fid\0maybesid\0sid\0chrom\0pos\0ref\0alt\0code\0error\0trionum\0ni\0nobsi\0fraci\0nl\0nobsl\0fracl\0", "mendel", kfMendelRptColMaybefid, kfMendelRptColDefault, 0, &pc.mendel_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              explicit_cols = 1;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mendel argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!explicit_cols) {
            pc.mendel_info.flags |= kfMendelRptColDefault;
          }
          pc.command_flags1 |= kfCommand1MendelReport;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "endel-duos")) {
          pc.mendel_info.flags |= kfMendelDuos;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "endel-missing-in-denom")) {
          if (unlikely(!((pc.filter_flags & kfFilterMendel) || (pc.command_flags1 & kfCommand1MendelReport)))) {
            logerrputs("Error: --mendel-missing-in-denom must be used with --me/--mendel.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.mendel_info.flags |= kfMendelMissingInDenom;
          goto main_param_zero;
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "endel-multigen"))) {
          logerrputs("Error: --mendel-multigen has been retired, since it has very limited practical\nvalue and is incompatible with how missing genotypes are now accounted for in\nMendel error rate calculations.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "erge-sids")) {
          if (unlikely(pc.sample_sort_mode == kSortFile)) {
            logerrputs("Error: --merge-sids cannot be used with an --indiv-sort file.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pmerge_info.flags |= kfPmergeSids;
          pmerge_required = 1;
          goto main_param_zero;
        } else if (likely(strequal_k_unsafe(flagname_p2, "pheno"))) {
          logerrputs("Warning: --mpheno flag deprecated.  Use --pheno-col-nums instead.  (Note that\n--pheno-col-nums does not add 2 to the column number(s).)\n");
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t mpheno_arg;
          if (unlikely(ScanPosintDefcapx(cur_modif, &mpheno_arg))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --mpheno argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          // add two to the number, convert it back to a string, and act as if
          // --pheno-col-nums was used.  See ParseNameRanges().
          const uint32_t pheno_col_nums_arg = mpheno_arg + 2;
          const uint32_t name_max_blen = UintSlen(pheno_col_nums_arg) + 1;
          if (unlikely(pgl_malloc(name_max_blen + 1, &pc.pheno_range_list.names))) {
            goto main_ret_NOMEM;
          }
          pc.pheno_range_list.name_max_blen = name_max_blen;
          pc.pheno_range_list.name_ct = 1;
          char* write_iter = u32toa(pheno_col_nums_arg, pc.pheno_range_list.names);
          *write_iter++ = '\0';
          pc.pheno_range_list.starts_range = R_CAST(unsigned char*, write_iter);
          *write_iter = '\0';
          pc.misc_flags |= kfMiscPhenoColNums;
        } else if (strequal_k_unsafe(flagname_p2, "erge")) {
          logerrputs("Error: --merge is retired.  Use --pmerge instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "erge-list")) {
          logerrputs("Error: --merge-list is retired.  Use --pmerge-list instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
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
          if (unlikely(pc.misc_flags & kfMiscAcFounders)) {
            logerrputs("Error: --ac-founders and --nonfounders cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.misc_flags |= kfMiscNonfounders;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ot-chr")) {
          if (unlikely(pc.varid_from)) {
            logerrputs("Error: --from/--to cannot be used with --autosome[-par] or --[not-]chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.from_bp != -1)) {
            logerrputs("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb cannot be used with --not-chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }

          // allowed:
          //   --allow-extra-chr --chr 5-22 bobs_chrom --not-chr 17
          // allowed:
          //   --allow-extra-chr --not-chr 12-17 bobs_chrom
          // does not make sense, disallowed:
          //   --allow-extra-chr --chr 5-22 --not-chr bobs_chrom
          // i.e. don't see a practical reason to make --not-chr accept
          // arbitrary chromosome codes when --chr or --autosome[-par] is
          // present
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }

          const uint32_t extra_chrs_prohibited_or_chr_include = (pc.misc_flags & kfMiscProhibitExtraChr) || chr_info.is_include_stack;
          reterr = ParseChrRanges(&(argvk[arg_idx]), flagname_p, errstr_append, param_ct, extra_chrs_prohibited_or_chr_include, kChrRawEnd - (kChrExcludeWords * kBitsPerWord), '-', &chr_info, chr_info.chr_exclude);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.load_filter_log_flags |= kfLoadFilterLogNotChr;
          // remaining processing now postponed to FinalizeChrset()
        } else if (strequal_k_unsafe(flagname_p2, "ew-id-max-allele-len")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(ScanPosintCappedx(cur_modif, kMaxIdSlen - 2, &pc.new_variant_id_max_allele_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --new-id-max-allele-len length argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (param_ct == 2) {
            cur_modif = argvk[arg_idx + 2];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "missing", cur_modif_slen)) {
              pc.misc_flags |= kfMiscNewVarIdOverflowMissing;
            } else if (strequal_k(cur_modif, "truncate", cur_modif_slen)) {
              pc.misc_flags |= kfMiscNewVarIdOverflowTruncate;
            } else if (unlikely(!strequal_k(cur_modif, "error", cur_modif_slen))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --new-id-max-allele-len argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "ormalize")) {
          // Prevent a confusing order-of-operations dependency.
          if (unlikely(pc.alt_allele_flag)) {
            logerrputs("Error: --normalize cannot be used with --ref-allele/--alt[1]-allele.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(!pc.fa_fname)) {
            logerrputs("Error: --normalize requires --fa.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "list", cur_modif_slen)) {
              pc.fa_flags |= kfFaNormalizeList;
            } else if (strequal_k(cur_modif, "shrink-overlapping-deletions", cur_modif_slen)) {
              logerrputs("Warning: 'shrink-overlapping-deletions' modifier has been renamed to\n'adjust-overlapping-deletions'.\n");
              pc.fa_flags |= kfFaNormalizeAdjustOverlappingDeletions;
            } else if (likely(strequal_k(cur_modif, "adjust-overlapping-deletions", cur_modif_slen))) {
              pc.fa_flags |= kfFaNormalizeAdjustOverlappingDeletions;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --normalize argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!allow_normalize_with_split) {
            if (unlikely(make_plink2_flags & (kfMakePlink2MSplitAll | kfMakePlink2MSplitSnps))) {
              logerrputs("Error: --normalize specified with a variant-split operation.  This probably\ndoesn't do what you want, since left-normalization occurs before variant-split\nin the " PROG_NAME_STR " order of operations.  Instead, you probably want to split first,\nand then left-normalize in a subsequent " PROG_NAME_STR " run.\nUse --allow-normalize-with-split to suppress this error.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          pc.fa_flags |= kfFaNormalize;
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "o-input-missing-phenotype")) {
          if (unlikely(pc.missing_pheno != -9)) {
            logerrputs("Error: --input-missing-phenotype and --no-input-missing-phenotype don't make\nsense together.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(pc.misc_flags & kfMiscNeg9PhenoReallyMissing)) {
            logerrputs("Error: --neg9-pheno-really-missing and --no-input-missing-phenotype don't make\nsense together.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.missing_pheno = 0;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ative")) {
#ifdef USE_MKL
          mkl_native = 1;
#endif
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ot-pheno")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenCommaDelim(&(argvk[arg_idx + 1]), param_ct, &pc.not_pheno_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ot-covar")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenCommaDelim(&(argvk[arg_idx + 1]), param_ct, &pc.not_covar_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "o-categorical")) {
          pc.misc_flags |= kfMiscNoCategorical;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "eg9-pheno-really-missing")) {
          pc.misc_flags |= kfMiscNeg9PhenoReallyMissing;
          goto main_param_zero;
        } else if (likely(strequal_k_unsafe(flagname_p2, "o-id-header"))) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (unlikely(strcmp(cur_modif, "iid-only"))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --no-id-header argument '%s'.\n", cur_modif);
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
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* mt_code = argvk[arg_idx + 1];
          const uint32_t code_slen = strlen(mt_code);
          if (strequal_k(mt_code, "M", code_slen)) {
            chr_info.output_encoding = kfChrOutputM;
          } else if (strequal_k(mt_code, "MT", code_slen)) {
            chr_info.output_encoding = kfChrOutputMT;
          } else if (strequal_k(mt_code, "0M", code_slen)) {
            chr_info.output_encoding = kfChrOutput0M;
          } else if (strequal_k(mt_code, "chr26", code_slen)) {
            chr_info.output_encoding = kfChrOutputPrefix;
          } else if (strequal_k(mt_code, "chrM", code_slen)) {
            chr_info.output_encoding = kfChrOutputPrefix | kfChrOutputM;
          } else if (strequal_k(mt_code, "chrMT", code_slen)) {
            chr_info.output_encoding = kfChrOutputPrefix | kfChrOutputMT;
          } else if (likely(strequal_k(mt_code, "26", code_slen))) {
            chr_info.output_encoding = kfChrOutput0;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --output-chr argument '%s'.\n", mt_code);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "utput-min-p")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely((!ScantokLn(cur_modif, &pc.output_min_ln)) || (pc.output_min_ln >= 0.0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --output-min-p argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "xford-single-chr")) {
          if (unlikely(!(xload & (kfXloadOxGen | kfXloadOxBgen)))) {
            logerrputs("Error: --oxford-single-chr must be used with .gen/.bgen input.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(oxford_import_flags & kfOxfordImportBgenSnpIdChr)) {
            logerrputs("Error: --oxford-single-chr cannot be used with --bgen 'snpid-chr'.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t chr_code_raw = GetChrCodeRaw(cur_modif);
          if (unlikely((chr_code_raw == UINT32_MAXM1) || ((pc.misc_flags & kfMiscProhibitExtraChr) && (chr_code_raw == UINT32_MAX)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --oxford-single-chr chromosome code '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          reterr = CmdlineAllocString(cur_modif, argvk[arg_idx], kMaxIdSlen, &import_single_chr_str);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "utput-missing-genotype")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          pc.output_missing_geno_char = ExtractCharParam(cur_modif);
          if (unlikely(ctou32(pc.output_missing_geno_char) <= ' ')) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --output-missing-genotype argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.legacy_output_missing_geno_char = pc.output_missing_geno_char;
        } else if (strequal_k_unsafe(flagname_p2, "utput-missing-phenotype")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t cur_modif_slen = strlen(cur_modif);
          if (unlikely(cur_modif_slen > 31)) {
            logerrputs("Error: --output-missing-phenotype string too long (max 31 chars).\n");
            goto main_ret_INVALID_CMDLINE;
          }
          // quasi-bugfix (2 Jan 2021): make this apply to .fam files
          memcpy(pc.legacy_output_missing_pheno, cur_modif, cur_modif_slen + 1);
          memcpy(pc.output_missing_pheno, cur_modif, cur_modif_slen + 1);
        } else if (unlikely(!strequal_k_unsafe(flagname_p2, "ut"))) {
          // --out is a special case due to logging
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'p':
        if (strequal_k_unsafe(flagname_p2, "file")) {
          if (unlikely(load_params || xload)) {
            // currently only possible with --bcf, --bfile, --pfile
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "vzs", &fname_modif_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* fname_prefix = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname_prefix);
          if (unlikely(slen > (kPglFnamesize - 11))) {
            logerrputs("Error: --pfile argument too long.\n");
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
          if (unlikely(xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          load_params |= kfLoadParamsPgen;
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (unlikely(slen > (kPglFnamesize - 2))) {
            logerrputs("Error: --pgen argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, fname, slen + 1);
        } else if (strequal_k_unsafe(flagname_p2, "sam")) {
          if (unlikely(xload & (~(kfXloadVcf | kfXloadBcf | kfXloadPlink1Dosage | kfXloadMap)))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          load_params |= kfLoadParamsPsam;
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (unlikely(slen > (kPglFnamesize - 2))) {
            logerrputs("Error: --psam argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(psamname, fname, slen + 1);
        } else if (strequal_k_unsafe(flagname_p2, "var")) {
          if (unlikely(xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          load_params |= kfLoadParamsPvar;
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (unlikely(slen > (kPglFnamesize - 2))) {
            logerrputs("Error: --pvar argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pvarname, fname, slen + 1);
        } else if (strequal_k_unsafe(flagname_p2, "heno")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "iid-only", &fname_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.misc_flags |= kfMiscPhenoIidOnly;
          }
          reterr = AllocFname(argvk[arg_idx + fname_idx], flagname_p, &pc.pheno_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "heno-col-nums")) {
          if (unlikely(!pc.pheno_fname)) {
            logerrputs("Error: --pheno-col-nums must be used with --pheno.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.pheno_range_list.name_ct)) {
            logerrputs("Error: --pheno-col-nums can't be used with --mpheno.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.pheno_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.misc_flags |= kfMiscPhenoColNums;
        } else if (strequal_k_unsafe(flagname_p2, "heno-name")) {
          if (unlikely(pc.pheno_range_list.name_ct)) {
            logerrputs("Error: --pheno-name can't be used with --pheno-col-nums.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // can now be used without --pheno
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.pheno_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "arallel")) {
          if (unlikely(pc.king_flags & kfKingMatrixSq)) {
            logerrputs("Error: --parallel cannot be used with \"--make-king square\".  Use \"--make-king\nsquare0\" or plain --make-king instead.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.command_flags1 & kfCommand1KingCutoff) && (!king_cutoff_fprefix))) {
            logerrputs("Error: --parallel cannot be used with no-input-fileset --king-cutoff.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.grm_flags & kfGrmMatrixSq)) {
            logerrputs("Error: --parallel cannot be used with \"--make-rel square\".  Use \"--make-rel\nsquare0\" or plain --make-rel instead.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(ScanPosintCappedx(argvk[arg_idx + 1], kParallelMax, &pc.parallel_idx))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --parallel job index '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(ScanPosintCappedx(argvk[arg_idx + 2], kParallelMax, &pc.parallel_tot) || (pc.parallel_tot == 1) || (pc.parallel_tot < pc.parallel_idx))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --parallel total job count '%s'.\n", argvk[arg_idx + 2]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          --pc.parallel_idx;  // internal 0..(n-1) indexing
        } else if (strequal_k_unsafe(flagname_p2, "arameters")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Glm))) {
            logerrputs("Error: --parameters must be used with --glm.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.glm_info.flags & kfGlmResidualizeMask)) {
            logerrputs("Error: --parameters cannot be used with --glm {cc,firth,qt}-residualize.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.glm_info.parameters_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "filter")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(!ScantokLn(cur_modif, &pc.ln_pfilter))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --pfilter argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely((pc.ln_pfilter == -DBL_MAX) || (pc.ln_pfilter > 0.0))) {
            logerrputs("Error: --pfilter threshold must be in (0, 1].\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ca")) {
#ifdef NOLAPACK
          logerrputs("Error: --pca requires " PROG_NAME_STR " to be built with LAPACK.\n");
          goto main_ret_INVALID_CMDLINE;
#endif
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 6))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t explicit_scols = 0;
          uint32_t vcols_idx = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "approx", cur_modif_slen)) {
              pc.pca_flags |= kfPcaApprox;
            } else if (strequal_k(cur_modif, "meanimpute", cur_modif_slen)) {
              pc.pca_flags |= kfPcaMeanimpute;
            } else if (strequal_k(cur_modif, "allele-wts", cur_modif_slen)) {
              pc.pca_flags |= kfPcaAlleleWts;
            } else if (strequal_k(cur_modif, "biallelic-var-wts", cur_modif_slen)) {
              pc.pca_flags |= kfPcaBiallelicVarWts;
            } else if (strequal_k(cur_modif, "vzs", cur_modif_slen)) {
              pc.pca_flags |= kfPcaVarZs;
            } else if (StrStartsWith0(cur_modif, "scols=", cur_modif_slen)) {
              if (unlikely(explicit_scols)) {
                logerrputs("Error: Multiple --pca scols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("scols=")]), "maybefid\0fid\0maybesid\0sid\0", "pca scols", kfPcaScolMaybefid, kfPcaScolDefault, 0, &pc.pca_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              explicit_scols = 1;
            } else if (StrStartsWith0(cur_modif, "vcols=", cur_modif_slen)) {
              if (unlikely(vcols_idx)) {
                logerrputs("Error: Multiple --pca vcols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              vcols_idx = param_idx;
            } else if (unlikely(strequal_k(cur_modif, "var-wts", cur_modif_slen))) {
              logerrputs("Error: --pca 'var-wts' modifier retired.  If your dataset contains no\nmultiallelic variants, the 'biallelic-var-wts' modifier has the same effect.\nOtherwise, migrate to 'allele-wts'.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (unlikely(strequal_k(cur_modif, "sid", cur_modif_slen))) {
              logerrputs("Error: --pca 'sid' modifier retired.  Use --pca scols= instead.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else {
              if (unlikely(pc.pca_ct || ScanPosintDefcapx(cur_modif, &pc.pca_ct))) {
                logerrputs("Error: Invalid --pca argument sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely(pc.pca_ct > kMaxPc)) {
                // this slightly simplifies output buffering.
                logerrputs("Error: --pca does not support more than 8000 PCs.\n");
                goto main_ret_INVALID_CMDLINE;
              }
            }
          }
          if (pc.pca_flags & kfPcaApprox) {
            if (unlikely(pc.pca_ct > 100)) {
              // double-precision overflow too likely
              logerrputs("Error: --pca approx does not support more than 100 PCs.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          } else {
            // todo: if --make-rel/--make-grm present, verify consistency
            if (unlikely(pc.parallel_tot != 1)) {
              logerrputs("Error: Non-approximate --pca cannot be used with --parallel.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            const uint32_t pca_meanimpute = (pc.pca_flags / kfPcaMeanimpute) & 1;
            if (pc.command_flags1 & kfCommand1MakeRel) {
              if (unlikely(((pc.grm_flags / kfGrmMeanimpute) & 1) != pca_meanimpute)) {
                logerrputs("Error: --make-rel/--make-grm-{bin,list,sparse} meanimpute setting must match\n--pca meanimpute setting.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              if (unlikely(pc.grm_flags & kfGrmCov)) {
                logerrputs("Error: --make-rel/--make-grm-{bin,list,sparse} cannot be used to compute a\ncovariance matrix in the same run as non-approximate --pca.\n");
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
          const PcaFlags wts_flags = pc.pca_flags & (kfPcaAlleleWts | kfPcaBiallelicVarWts);
          if (wts_flags) {
            if (wts_flags == (kfPcaAlleleWts | kfPcaBiallelicVarWts)) {
              // could support this, but don't see the point
              logerrputs("Error: --pca 'allele-wts' and 'biallelic-var-wts' modifiers cannot be used\ntogether.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            const uint32_t is_var_wts = (wts_flags / kfPcaBiallelicVarWts) & 1;
            const PcaFlags default_cols = is_var_wts? kfPcaVcolDefaultB : kfPcaVcolDefaultA;
            if (!vcols_idx) {
              pc.pca_flags |= default_cols;
            } else {
              const char* cur_modif = argvk[arg_idx + vcols_idx];
              reterr = ParseColDescriptor(&(cur_modif[strlen("vcols=")]), "chrom\0pos\0ref\0alt1\0alt\0maybeprovref\0provref\0ax\0maj\0nonmaj\0", "pca vcols", kfPcaVcolChrom, default_cols, is_var_wts, &pc.pca_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              if (is_var_wts) {
                if (pc.pca_flags & kfPcaVcolAx) {
                  logerrputs("Error: \"--pca biallelic-allele-wts\" doesn't support the 'ax' column set.\n");
                  goto main_ret_INVALID_CMDLINE;
                }
              } else {
                if (pc.pca_flags & (kfPcaVcolMaj | kfPcaVcolNonmaj)) {
                  logerrputs("Error: \"--pca allele-wts\" doesn't support the 'maj' and 'nonmaj' column sets.\n");
                  goto main_ret_INVALID_CMDLINE;
                }
              }
            }
          } else {
            if (unlikely(pc.pca_flags & kfPcaVarZs)) {
              logerrputs("Error: --pca 'vzs' modifier has no effect when variant/allele weights have not\nbeen requested.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            if (unlikely(vcols_idx)) {
              logerrputs("Error: --pca 'vcols=' has no effect when variant/allele weights have not been\nrequested.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          }
          pc.command_flags1 |= kfCommand1Pca;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "gen-info")) {
          pc.command_flags1 |= kfCommand1PgenInfo;
          pc.dependency_flags |= kfFilterAllReq;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "merge")) {
          if (unlikely(import_flags & kfImportKeepAutoconv)) {
            logerrputs("Error: --pmerge cannot be used with --keep-autoconv.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct < 3) {
            uint32_t fname_idx = 1;
            if (param_ct == 2) {
              if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "vzs", &fname_idx))) {
                goto main_ret_INVALID_CMDLINE_A;
              }
            }
            const uint32_t vzs = (param_ct == 2);
            const char* prefix = argvk[arg_idx + fname_idx];
            const uint32_t prefix_slen = strlen(prefix);
            const uint32_t fname_alloc = prefix_slen + 6;
            if (unlikely(fname_alloc + vzs * 4 > kPglFnamesize)) {
              logerrputs("Error: --pmerge argument too long.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            if (unlikely(pgl_malloc(fname_alloc, &pmerge_info.pgen_fname) ||
                         pgl_malloc(fname_alloc + vzs * 4, &pmerge_info.pvar_fname) ||
                         pgl_malloc(fname_alloc, &pmerge_info.psam_fname))) {
              goto main_ret_NOMEM;
            }
            memcpy(pmerge_info.pgen_fname, prefix, prefix_slen);
            memcpy(&(pmerge_info.pgen_fname[prefix_slen]), ".pgen", 6);
            memcpy(pmerge_info.pvar_fname, prefix, prefix_slen);
            memcpy(&(pmerge_info.pvar_fname[prefix_slen]), ".pvar", 6);
            if (vzs) {
              memcpy(&(pmerge_info.pvar_fname[prefix_slen + 5]), ".zst", 5);
            }
            memcpy(pmerge_info.psam_fname, prefix, prefix_slen);
            memcpy(&(pmerge_info.psam_fname[prefix_slen]), ".psam", 6);
          } else {
            reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pmerge_info.pgen_fname);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            reterr = AllocFname(argvk[arg_idx + 2], flagname_p, &pmerge_info.pvar_fname);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            reterr = AllocFname(argvk[arg_idx + 3], flagname_p, &pmerge_info.psam_fname);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.command_flags1 |= kfCommand1Pmerge;
        } else if (strequal_k_unsafe(flagname_p2, "merge-list")) {
          if (pc.command_flags1 & kfCommand1Pmerge) {
            logerrputs("Error: --pmerge-list cannot be used with --pmerge.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(import_flags & kfImportKeepAutoconv)) {
            logerrputs("Error: --pmerge-list cannot be used with --keep-autoconv.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_idx = 0;
          uint32_t explicit_mode = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "bfile", cur_modif_slen)) {
              ++explicit_mode;
              pmerge_info.list_mode = kPmergeListModeBfile;
            } else if (strequal_k(cur_modif, "bpfile", cur_modif_slen)) {
              ++explicit_mode;
              pmerge_info.list_mode = kPmergeListModeBpfile;
            } else if (strequal_k(cur_modif, "pfile", cur_modif_slen)) {
              ++explicit_mode;
              pmerge_info.list_mode = kPmergeListModePfile;
            } else if (strequal_k(cur_modif, "pfile-vzs", cur_modif_slen)) {
              ++explicit_mode;
              pmerge_info.list_mode = kPmergeListModePfileVzs;
            } else {
              if (unlikely(fname_idx)) {
                logerrputs("Error: Invalid --pmerge-list argument sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              fname_idx = param_idx;
            }
          }
          if (unlikely(explicit_mode > 1)) {
            logerrputs("Error: Multiple --pmerge-list modes specified.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(!fname_idx)) {
            logerrputs("Error: No filename provided for --pmerge-list.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = AllocFname(argvk[arg_idx + fname_idx], flagname_p, &pmerge_info.list_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.command_flags1 |= kfCommand1Pmerge;
        } else if (strequal_k_unsafe(flagname_p2, "merge-list-dir")) {
          if (unlikely((!(pc.command_flags1 & kfCommand1Pmerge)) || (!pmerge_info.list_fname))) {
            logerrputs("Error: --pmerge-list-dir must be used with --pmerge-list.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* dir_str = argvk[arg_idx + 1];
          const uint32_t dir_slen = strlen(dir_str);
          uint32_t append_path_sep = 0;
#ifdef _WIN32
          const char path_sep = '\\';
          if ((dir_str[dir_slen - 1] != '/') && (dir_str[dir_slen - 1] != '\\')) {
            append_path_sep = 1;
          }
#else
          const char path_sep = '/';
          if (dir_str[dir_slen - 1] != '/') {
            append_path_sep = 1;
          }
#endif
          const uint32_t alloc_blen = dir_slen + 1 + append_path_sep;
          if (unlikely(alloc_blen > (kPglFnamesize - 9))) {
            logerrputs("Error: --pmerge-list-dir argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          if (unlikely(pgl_malloc(alloc_blen, &pmerge_info.list_base_dir))) {
            goto main_ret_NOMEM;
          }
          char* write_iter = memcpya(pmerge_info.list_base_dir, dir_str, dir_slen);
          if (append_path_sep) {
            *write_iter++ = path_sep;
          }
          *write_iter = '\0';
        } else if (strequal_k_unsafe(flagname_p2, "merge-output-vzs")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Pmerge))) {
            logerrputs("Error: --pmerge-output-vzs must be used with --pmerge or --pmerge-list.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pmerge_info.flags |= kfPmergeOutputVzs;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "gen-diff")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 7))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t first_optional_param_idx = (param_ct > 6)? 4 : 2;
          uint32_t explicit_cols = 0;
          uint32_t vzs = 0;
          for (uint32_t param_idx = first_optional_param_idx; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "vzs", cur_modif_slen)) {
              vzs = 1;
            } else if (strequal_k(cur_modif, "include-missing", cur_modif_slen)) {
              pc.pgen_diff_info.flags |= kfPgenDiffIncludeMissing;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.pgen_diff_info.flags |= kfPgenDiffZs;
            } else if (StrStartsWith0(cur_modif, "dosage", cur_modif_slen)) {
              if (unlikely(pc.pgen_diff_info.dosage_hap_tol != kDosageMissing)) {
                logerrputs("Error: Multiple --pgen-diff dosage modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              if (cur_modif[6] == '\0') {
                pc.pgen_diff_info.dosage_hap_tol = 0;
              } else {
                if (unlikely(cur_modif[6] != '=')) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --pgen-diff argument '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                double dxx;
                if (unlikely((!ScantokDouble(&(cur_modif[7]), &dxx)) || (dxx < 0.0) || (dxx > (0.5 - kSmallEpsilon)))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --pgen-diff argument '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                pc.pgen_diff_info.dosage_hap_tol = S_CAST(int32_t, dxx * ((1 + kSmallEpsilon) * kDosageMax));
              }
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (unlikely(explicit_cols)) {
                logerrputs("Error: Multiple --pgen-diff cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0id\0ref\0alt\0maybeprovref\0provref\0maybefid\0fid\0maybesid\0sid\0geno\0", "pgen-diff", kfPgenDiffColChrom, kfPgenDiffColDefault, 0, &pc.pgen_diff_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (likely(param_idx == 2)) {
              first_optional_param_idx = 4;
              param_idx = 3;
            } else {
              logerrputs("Error: Invalid --pgen-diff argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          if (first_optional_param_idx == 2) {
            const char* prefix = argvk[arg_idx + 1];
            const uint32_t prefix_slen = strlen(prefix);
            const uint32_t fname_alloc = prefix_slen + 6;
            if (unlikely(fname_alloc + vzs * 4 > kPglFnamesize)) {
              logerrputs("Error: --pgen-diff argument too long.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            if (unlikely(pgl_malloc(fname_alloc, &pc.pgen_diff_info.pgen_fname) ||
                         pgl_malloc(fname_alloc + vzs * 4, &pc.pgen_diff_info.pvar_fname) ||
                         pgl_malloc(fname_alloc, &pc.pgen_diff_info.psam_fname))) {
              goto main_ret_NOMEM;
            }
            memcpy(pc.pgen_diff_info.pgen_fname, prefix, prefix_slen);
            memcpy(&(pc.pgen_diff_info.pgen_fname[prefix_slen]), ".pgen", 6);
            memcpy(pc.pgen_diff_info.pvar_fname, prefix, prefix_slen);
            memcpy(&(pc.pgen_diff_info.pvar_fname[prefix_slen]), ".pvar", 6);
            if (vzs) {
              memcpy(&(pc.pgen_diff_info.pvar_fname[prefix_slen + 5]), ".zst", 5);
            }
            memcpy(pc.pgen_diff_info.psam_fname, prefix, prefix_slen);
            memcpy(&(pc.pgen_diff_info.psam_fname[prefix_slen]), ".psam", 6);
          } else {
            if (unlikely(vzs)) {
              logerrputs("Error: --pgen-diff 'vzs' modifier can only be used when a single fileset prefix\nis specified.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.pgen_diff_info.pgen_fname);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            reterr = AllocFname(argvk[arg_idx + 2], flagname_p, &pc.pgen_diff_info.pvar_fname);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            reterr = AllocFname(argvk[arg_idx + 3], flagname_p, &pc.pgen_diff_info.psam_fname);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          if (!explicit_cols) {
            pc.pgen_diff_info.flags |= kfPgenDiffColDefault;
          }
          pc.command_flags1 |= kfCommand1PgenDiff;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "heno-inner-join")) {
          pmerge_info.flags |= kfPmergePhenoInnerJoin;
          pmerge_required = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ed")) {
          if (unlikely(load_params || (xload & (~kfXloadMap)))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --ped argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, fname, slen + 1);
          xload |= kfXloadPed;
        } else if (strequal_k_unsafe(flagname_p2, "edmap")) {
          if (unlikely(load_params || (xload & (~(kfXloadPed | kfXloadMap))))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* fname_prefix = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname_prefix);
          if (unlikely(slen > kPglFnamesize - 5)) {
            logerrputs("Error: --pedmap argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          if (!(xload & kfXloadMap)) {
            // allow --map to override this
            char* prefix_end = memcpya(pvarname, fname_prefix, slen);
            snprintf(prefix_end, 5, ".map");
            xload |= kfXloadMap;
          }
          if (!(xload & kfXloadPed)) {
            // allow --ped to override this
            char* prefix_end = memcpya(pgenname, fname_prefix, slen);
            snprintf(prefix_end, 5, ".ped");
            xload |= kfXloadPed;
          }
        } else if (strequal_k_unsafe(flagname_p2, "gi")) {
          if (unlikely(xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(!(pgenname[0]))) {
            // this doesn't catch e.g. --pgi + --bfile, but an imperfect sanity
            // check is better than no check
            logerrputs("Error: --pgi must be used with --pfile or --pgen.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          // AllocFnamePrefix() should be used if there's ever a possibility of
          // appending '~'.
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.pginame);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "olyploid-mode")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* mode_str = argvk[arg_idx + 1];
          const char first_char_upcase_match = mode_str[0] & 0xdf;
          const uint32_t mode_slen = strlen(mode_str);
          if (((mode_slen == 1) && (first_char_upcase_match == 'E')) ||
              strequal_k(mode_str, "error", mode_slen)) {
            import_flags |= kfImportPolyploidExplicitError;
          } else if (likely(((mode_slen == 1) && (first_char_upcase_match == 'M')) ||
                            strequal_k(mode_str, "missing", mode_slen))) {
            import_flags |= kfImportPolyploidMissing;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a valid mode for --polyploid-mode.\n", mode_str);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "heno-svd")) {
#ifdef NOLAPACK
          logerrputs("Error: --pheno-svd requires " PROG_NAME_STR " to be built with LAPACK.\n");
          goto main_ret_INVALID_CMDLINE;
#endif
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* first_arg = argvk[arg_idx + 1];
          const uint32_t first_arg_slen = strlen(first_arg);
          if (StrStartsWith(first_arg, "variance=", first_arg_slen)) {
            double dxx;
            if (unlikely((!ScantokDouble(&(first_arg[9]), &dxx)) || (dxx <= 0.0) || (dxx > 1.0))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --pheno-svd variance= argument '%s'.\n", first_arg);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            pc.pheno_svd_info.min_variance_explained = dxx;
          } else {
            uint32_t uii;
            if (unlikely(ScanPosintDefcapx(first_arg, &uii))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --pheno-svd output phenotype count '%s'.\n", first_arg);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (unlikely(uii > kMaxPc)) {
              logerrputs("Error: --pheno-svd does not support more than 8000 output phenotypes.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            pc.pheno_svd_info.ct = uii;
          }
          uint32_t explicit_scols = 0;
          uint32_t explicit_pcols = 0;
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "force", cur_modif_slen)) {
              pc.pheno_svd_info.flags |= kfPhenoSvdForce;
            } else if (StrStartsWith0(cur_modif, "scols=", cur_modif_slen)) {
              if (unlikely(explicit_scols)) {
                logerrputs("Error: Multiple --pheno-svd scols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_scols = 1;
              reterr = ParseColDescriptor(&(cur_modif[strlen("scols=")]), "maybefid\0fid\0maybesid\0sid\0", "pheno-svd scols", kfPhenoSvdScolMaybefid, kfPhenoSvdScolDefault, 0, &pc.pheno_svd_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (likely(StrStartsWith0(cur_modif, "pcols=", cur_modif_slen))) {
              if (unlikely(explicit_scols)) {
                logerrputs("Error: Multiple --pheno-svd pcols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_scols = 1;
              reterr = ParseColDescriptor(&(cur_modif[strlen("pcols=")]), "id\0sv\0", "pheno-svd pcols", kfPhenoSvdPcolId, kfPhenoSvdPcolDefault, 0, &pc.pheno_svd_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              logerrputs("Error: Invalid --pheno-svd argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          if (!explicit_scols) {
            pc.pheno_svd_info.flags |= kfPhenoSvdScolDefault;
          }
          if (!explicit_pcols) {
            pc.pheno_svd_info.flags |= kfPhenoSvdPcolDefault;
          }
          pc.command_flags1 |= kfCommand1PhenoSvd;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (likely(strequal_k_unsafe(flagname_p2, "heno-quantile-normalize"))) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.quantnorm_flattened);
            if (unlikely(reterr)) {
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
          if (unlikely(pc.pheno_transform_flags & (kfPhenoTransformQuantnormPheno | kfPhenoTransformQuantnormCovar))) {
            logerrputs("Error: --quantile-normalize cannot be used with --pheno-quantile-normalize or\n--covar-quantile-normalize.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.quantnorm_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.pheno_transform_flags |= kfPhenoTransformQuantnormAll;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (likely(strequal_k_unsafe(flagname_p2, "-score-range"))) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 2, 6))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.score_info.qsr_range_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          reterr = AllocFname(argvk[arg_idx + 2], flagname_p, &pc.score_info.qsr_data_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          uint32_t numeric_param_ct = 0;
          uint32_t qsr_cols[3];
          for (uint32_t param_idx = 3; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "header", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreQsrHeader;
            } else if (strequal_k(cur_modif, "min", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreQsrMin;
            } else {
              if (unlikely(ScanPosintCappedx(cur_modif, kMaxLongLine / 2, &(qsr_cols[numeric_param_ct])))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --q-score-range argument '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (unlikely(numeric_param_ct == 2)) {
                logerrputs("Error: --q-score-range takes at most two numeric arguments.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely(numeric_param_ct && (qsr_cols[0] == qsr_cols[1]))) {
                logerrputs("Error: Identical --q-score-range column indexes.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              ++numeric_param_ct;
            }
          }
          if (numeric_param_ct) {
            pc.score_info.qsr_varid_col_p1 = qsr_cols[0];
          }
          if (numeric_param_ct == 2) {
            pc.score_info.qsr_val_col_p1 = qsr_cols[1];
          } else {
            pc.score_info.qsr_val_col_p1 = pc.score_info.qsr_varid_col_p1 + 1;
          }
          // no dependencies since we enforce --score requirement later
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'r':
        if (strequal_k_unsafe(flagname_p2, "eal-ref-alleles")) {
          if (unlikely(pc.misc_flags & kfMiscMajRef)) {
            logerrputs("Error: --real-ref-alleles cannot be used with --maj-ref.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.misc_flags |= kfMiscRealRefAlleles;
          // --real-ref-alleles applies to --pmerge[-list] input.
          if (!(pc.command_flags1 & kfCommand1Pmerge)) {
            pc.dependency_flags |= kfFilterAllReq;
          }
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "emove")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.remove_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-fam")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.removefam_fnames);
          if (unlikely(reterr)) {
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
          if (unlikely(pc.command_flags1 & kfCommand1AlleleFreq)) {
            // --read-freq can't promise OBS_CTs, so simplest to just continue
            // disallowing this
            logerrputs("Error: --freq and --read-freq cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.misc_flags & kfMiscMajRef)) {
            logerrputs("Error: --read-freq cannot be used with --maj-ref.  (We recommend preprocessing\nthe --read-freq file to generate a file that can be used with\n--ref-allele/--alt-allele.)\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.read_freq_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "equire-pheno")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.require_pheno_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.misc_flags |= kfMiscRequirePheno;
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "equire-covar")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.require_covar_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.misc_flags |= kfMiscRequireCovar;
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "equire-info")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenCommaDelim(&(argvk[arg_idx + 1]), param_ct, &pc.require_info_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogRequireInfo;
        } else if (strequal_k_unsafe(flagname_p2, "equire-no-info")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenCommaDelim(&(argvk[arg_idx + 1]), param_ct, &pc.require_no_info_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogRequireNoInfo;
        } else if (strequal_k_unsafe(flagname_p2, "emove-if")) {
          reterr = ValidateAndAllocCmpExpr(&(argvk[arg_idx + 1]), argvk[arg_idx], param_ct, &pc.remove_if_expr);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-cats")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.remove_cats_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-cat-names")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, kMaxIdBlen, &pc.remove_cat_names_flattened);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "emove-cat-pheno")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.remove_cat_phenoname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ef-allele")) {
          if (unlikely(pc.misc_flags & kfMiscMajRef)) {
            logerrputs("Error: --maj-ref cannot be used with --ref-allele/--alt[1]-allele.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.fa_flags & kfFaNormalize)) {
            logerrputs("Error: --normalize cannot be used with --ref-allele/--alt[1]-allele.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 5))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* const* sources = &(argvk[arg_idx + 1]);
          if (!strcmp(sources[0], "force")) {
            --param_ct;
            if (unlikely(!param_ct)) {
              logerrputs("Error: Invalid --ref-allele argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.misc_flags |= kfMiscRefAlleleForce;
            ++sources;
          }
          reterr = Alloc2col(sources, flagname_p, param_ct, &pc.ref_allele_flag);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "ef-from-fa")) {
          if (unlikely((pc.misc_flags & kfMiscMajRef) || pc.alt_allele_flag)) {
            // could allow --alt[1]-allele later, but keep this simpler for now
            logerrputs("Error: --ref-from-fa cannot be used with --maj-ref or\n--ref-allele/--alt[1]-allele.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "force", cur_modif_slen)) {
              pc.fa_flags |= kfFaRefFromForce;
            } else if (likely((!fname_modif_idx) && (!pc.fa_fname))) {
              fname_modif_idx = param_idx;
            } else {
              logerrputs("Error: Invalid --ref-from-fa argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          if (fname_modif_idx) {
            reterr = AllocFname(argvk[arg_idx + fname_modif_idx], flagname_p, &pc.fa_fname);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            logerrputs("Warning: Filename-argument form of --ref-from-fa is deprecated.  Use --fa to\nspecify the .fa file instead.\n");
          } else if (unlikely(!pc.fa_fname)) {
            logerrputs("Error: --ref-from-fa requires --fa.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.fa_flags |= kfFaRefFrom;
          pc.dependency_flags |= kfFilterPvarReq | kfFilterNonrefFlagsNeededSet;
        } else if (strequal_k_unsafe(flagname_p2, "m-dup")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          RmDupMode rmdup_mode = kRmDup0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "list", cur_modif_slen)) {
              pc.command_flags1 |= kfCommand1RmDupList;
            } else if (rmdup_mode != kRmDup0) {
              logerrputs("Error: Invalid --rm-dup argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            } else if (strequal_k(cur_modif, "error", cur_modif_slen)) {
              rmdup_mode = kRmDupError;
            } else if (strequal_k(cur_modif, "retain-mismatch", cur_modif_slen)) {
              rmdup_mode = kRmDupRetainMismatch;
            } else if (strequal_k(cur_modif, "exclude-mismatch", cur_modif_slen)) {
              rmdup_mode = kRmDupExcludeMismatch;
            } else if (strequal_k(cur_modif, "exclude-all", cur_modif_slen)) {
              rmdup_mode = kRmDupExcludeAll;
            } else if (likely(strequal_k(cur_modif, "force-first", cur_modif_slen))) {
              rmdup_mode = kRmDupForceFirst;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --rm-dup argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (rmdup_mode < kRmDupExcludeAll) {
            if (rmdup_mode == kRmDup0) {
              rmdup_mode = kRmDupError;
            }
            pc.dependency_flags |= kfFilterOpportunisticPgen;
          }
          pc.rmdup_mode = rmdup_mode;
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "ecover-var-ids")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_idx = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "strict-bim-order", cur_modif_slen)) {
              pc.recover_var_ids_flags |= kfRecoverVarIdsStrictBimOrder;
            } else if (strequal_k(cur_modif, "rigid", cur_modif_slen)) {
              pc.recover_var_ids_flags |= kfRecoverVarIdsRigid;
            } else if (strequal_k(cur_modif, "force", cur_modif_slen)) {
              pc.recover_var_ids_flags |= kfRecoverVarIdsForce;
            } else if (strequal_k(cur_modif, "partial", cur_modif_slen)) {
              pc.recover_var_ids_flags |= kfRecoverVarIdsPartial;
            } else {
              if (unlikely(fname_idx)) {
                logerrputs("Error: Invalid --recover-var-ids argument sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              fname_idx = param_idx;
            }
          }
          if (unlikely(!fname_idx)) {
            logerrputs("Error: No filename provided for --recover-var-ids.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.recover_var_ids_flags & (kfRecoverVarIdsRigid | kfRecoverVarIdsForce)) == (kfRecoverVarIdsRigid | kfRecoverVarIdsForce))) {
            logerrputs("Error: --recover-var-ids 'rigid' and 'force' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = AllocFname(argvk[arg_idx + fname_idx], flagname_p, &pc.recover_var_ids_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "2-phased") || strequal_k_unsafe(flagname_p2, "2-unphased") || strequal_k_unsafe(flagname_p2, "-phased") || strequal_k_unsafe(flagname_p2, "-unphased")) {
          if (unlikely(pc.command_flags1 & kfCommand1Vcor)) {
            logerrputs("Error: --r-phased, --r-unphased, --r2-phased, and --r2-unphased are mutually\nexclusive.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 6))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const uint32_t is_unsquared = (flagname_p2[0] == '-');
          const uint32_t is_phased = (flagname_p2[2 - is_unsquared] == 'p');
          uint32_t cols_idx = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "square", cur_modif_slen)) {
              if (unlikely(pc.vcor_info.flags & kfVcorMatrixShapemask)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Multiple --%s shape modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              pc.vcor_info.flags |= kfVcorMatrixSq;
            } else if (strequal_k(cur_modif, "square0", cur_modif_slen)) {
              if (unlikely(pc.vcor_info.flags & kfVcorMatrixShapemask)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Multiple --%s shape modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              pc.vcor_info.flags |= kfVcorMatrixSq0;
            } else if (strequal_k(cur_modif, "triangle", cur_modif_slen)) {
              if (unlikely(pc.vcor_info.flags & kfVcorMatrixShapemask)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Multiple --%s shape modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              pc.vcor_info.flags |= kfVcorMatrixTri;
            } else if (strequal_k(cur_modif, "inter-chr", cur_modif_slen)) {
              pc.vcor_info.flags |= kfVcorInterChr;
            } else if (strequal_k(cur_modif, "yes-really", cur_modif_slen)) {
              pc.vcor_info.flags |= kfVcorYesReally;
            } else if (strequal_k(cur_modif, "ref-based", cur_modif_slen)) {
              pc.vcor_info.flags |= kfVcorRefBased;
            } else if (strequal_k(cur_modif, "allow-ambiguous-allele", cur_modif_slen)) {
              pc.vcor_info.flags |= kfVcorAllowAmbiguousAllele;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              if (unlikely(pc.vcor_info.flags & kfVcorEncodemask)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Multiple --%s encoding modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              pc.vcor_info.flags |= kfVcorZs;
            } else if (strequal_k(cur_modif, "bin", cur_modif_slen)) {
              if (unlikely(pc.vcor_info.flags & kfVcorEncodemask)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Multiple --%s encoding modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              pc.vcor_info.flags |= kfVcorBin8;
            } else if (strequal_k(cur_modif, "bin4", cur_modif_slen)) {
              if (unlikely(pc.vcor_info.flags & kfVcorEncodemask)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Multiple --%s encoding modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              pc.vcor_info.flags |= kfVcorBin4;
            } else if (likely(StrStartsWith(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(cols_idx)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Multiple --%s cols= modifiers.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              // Postpone parse, since we need to know whether ref-based is
              // present.
              cols_idx = param_idx;
            } else if (strequal_k(cur_modif, "in-phase", cur_modif_slen)) {
              if (is_unsquared) {
                snprintf(g_logbuf, kLogbufSize, "Error: --%s does not have PLINK 1.x's 'in-phase' modifier. (The same information is available from the major-allele columns and the sign of R.)\n", flagname_p);
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: --%s does not have PLINK 1.x's 'in-phase' modifier. Instead, run --r-%sphased, then look at the major-allele columns and the sign of R.\n", flagname_p, is_phased? "" : "un");
              }
              goto main_ret_INVALID_CMDLINE_WWA;
            } else if (strequal_k(cur_modif, "with-freqs", cur_modif_slen)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --%s does not have PLINK 1.x's 'with-freqs' modifier. Use 'cols=+freq' instead.\n", flagname_p);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else if (strequal_k(cur_modif, "d", cur_modif_slen) || strequal_k(cur_modif, "dprime", cur_modif_slen) || strequal_k(cur_modif, "dprime-signed", cur_modif_slen)) {
              const char* new_col_name;
              if (cur_modif_slen == 1) {
                new_col_name = cur_modif;
              } else if (cur_modif_slen == 6) {
                new_col_name = "dprimeabs";
              } else {
                new_col_name = "dprime";
              }
              if (is_phased) {
                snprintf(g_logbuf, kLogbufSize, "Error: --%s does not have PLINK 1.x's '%s' modifier. Use 'cols=+%s' instead.\n", flagname_p, cur_modif, new_col_name);
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: --%s does not support computation of D or D'. Use --r%s-phased with 'cols=+%s' instead.\n", flagname_p, is_unsquared? "" : "2", new_col_name);
              }
              goto main_ret_INVALID_CMDLINE_WWA;
            } else if (strequal_k(cur_modif, "spaces", cur_modif_slen)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --%s does not support space-delimited output. (You can use e.g. the Unix 'tr' command to convert tabs to spaces.)\n", flagname_p);
              goto main_ret_INVALID_CMDLINE_WWA;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s argument '%s'.\n", flagname_p, cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (cols_idx) {
            VcorFlags default_cols = kfVcorColDefault;
            if (is_unsquared) {
              default_cols |= (pc.vcor_info.flags & kfVcorRefBased)? kfVcorColRef : kfVcorColMaj;
            }
            const char* cur_modif = argvk[arg_idx + cols_idx];
            reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0id\0ref\0alt1\0alt\0maybeprovref\0provref\0maj\0nonmaj\0freq\0d\0dprime\0dprimeabs\0", flagname_p, kfVcorColChrom, default_cols, 0, &pc.vcor_info.flags);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            if (unlikely((!is_phased) && (pc.vcor_info.flags & (kfVcorColD | kfVcorColDprime | kfVcorColDprimeAbs)))) {
              snprintf(g_logbuf, kLogbufSize, "Error: --%s does not support computation of D or D'. Use --r%s-phased instead.\n", flagname_p, is_unsquared? "" : "2");
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          const uint32_t must_matrix = (pc.vcor_info.flags & (kfVcorBin8 | kfVcorBin4 | kfVcorMatrixShapemask));
          const uint32_t must_all_pairs = must_matrix || (pc.vcor_info.flags & kfVcorInterChr);
          if (must_all_pairs) {
            if (unlikely((pc.vcor_info.var_ct_radius != 0x7fffffff) || (pc.vcor_info.bp_radius != UINT32_MAX) || (pc.vcor_info.cm_radius != -1.0))) {
              snprintf(g_logbuf, kLogbufSize, "Error: All-pairs --%s settings cannot be used with --ld-window/--ld-window-kb/--ld-window-cm.\n", flagname_p);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            if (pc.vcor_info.bp_radius == UINT32_MAX) {
              pc.vcor_info.bp_radius = 1000000;
            }
          }
          if (must_matrix) {
            const uint32_t must_table = cols_idx || (pc.vcor_info.flags & kfVcorInterChr) || pc.vcor_info.ld_snp_list_fname || pc.vcor_info.ld_snp_range_list.name_ct || (pc.vcor_info.min_r2 != 2.0);
            if (unlikely(must_matrix && must_table)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Matrix-only and table-only --%s settings cannot be used together.\n", flagname_p);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            if (!cols_idx) {
              pc.vcor_info.flags |= kfVcorColDefault;
              if (is_unsquared) {
                if (pc.vcor_info.flags & kfVcorRefBased) {
                  pc.vcor_info.flags |= kfVcorColRef;
                } else {
                  pc.vcor_info.flags |= kfVcorColMaj;
                }
              }
            }
            if (pc.vcor_info.min_r2 == 2.0) {
              pc.vcor_info.min_r2 = 0.2 * (1 - kSmallEpsilon);
            }
          }
          pc.vcor_info.flags |= S_CAST(VcorFlags, (is_phased * kfVcorPhased) | (is_unsquared * kfVcorUnsquared));
          pc.command_flags1 |= kfCommand1Vcor;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "andmem")) {
          randmem = 1;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ename-chrs")) {
          if (unlikely(chr_info.is_include_stack || (pc.load_filter_log_flags & kfLoadFilterLogNotChr))) {
            // too confusing (does the chromosome filter apply before?  after?
            // both?)
            logerrputs("Error: --rename-chrs cannot be used with --autosome[-par] or --[not-]chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.recover_var_ids_fname)) {
            logerrputs("Error: --rename-chrs cannot be used with --recover-var-ids or --update-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.rename_chrs_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (likely(strequal_k_unsafe(flagname_p2, "ice"))) {
          if (unlikely(chr_info.chrset_source)) {
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
        } else if (strequal_k_unsafe(flagname_p2, "2")) {
          logerrputs("Error: --r2 has been replaced by --r2-phased and --r2-unphased.\n--r2-phased computes the textbook haplotype-frequency-based r^2, and\ncorresponds to PLINK 1.9 --r2's behavior when the 'd', 'dprime', or\n'dprime-signed' modifier was present.  --r2-unphased computes the simpler r^2\nbetween (unphased) dosage vectors, and corresponds to how PLINK 1.x --r2\nbehaved without a D'-related modifier.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "")) {
          logerrputs("Error: --r has been replaced by --r-phased and --r-unphased.\n--r-phased computes the textbook haplotype-frequency-based r, and corresponds\nto PLINK 1.9 --r's behavior when the 'd', 'dprime', or 'dprime-signed' modifier\nwas present.  --r-unphased computes the simpler r between (unphased) dosage\nvectors, and corresponds to how PLINK 1.x --r behaved without a D'-related\nmodifier.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 's':
        if (strequal_k_unsafe(flagname_p2, "eed")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          rseed_ct = param_ct;
          if (unlikely(pgl_malloc(param_ct * sizeof(int32_t), &rseeds))) {
            goto main_ret_NOMEM;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (unlikely(ScanUintCappedx(cur_modif, UINT32_MAX, &(rseeds[param_idx - 1])))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --seed argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "plit-par")) {
          if (unlikely(pc.misc_flags & (kfMiscMergePar | kfMiscMergeX))) {
            logerrputs("Error: --split-par cannot be used with --merge-par/--merge-x.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 1) {
            const char* build_code = argvk[arg_idx + 1];
            const uint32_t code_slen = strlen(build_code);
            if (strequal_k(build_code, "chm13", code_slen)) {
              pc.splitpar_bound1 = 2394410;
              pc.splitpar_bound2 = 153925835;
            } else if (strequal_k(build_code, "b38", code_slen) ||
                strequal_k(build_code, "hg38", code_slen)) {
              pc.splitpar_bound1 = 2781479;
              pc.splitpar_bound2 = 155701383;
            } else if (strequal_k(build_code, "b37", code_slen) ||
                       strequal_k(build_code, "hg19", code_slen)) {
              pc.splitpar_bound1 = 2699520;
              pc.splitpar_bound2 = 154931044;
            } else if (likely(strequal_k(build_code, "b36", code_slen) ||
                              strequal_k(build_code, "hg18", code_slen))) {
              pc.splitpar_bound1 = 2709521;
              pc.splitpar_bound2 = 154584237;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Unrecognized --split-par build code '%s'.\n", build_code);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            if (unlikely(ScanUintDefcapx(argvk[arg_idx + 1], &pc.splitpar_bound1))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --split-par argument '%s'.\n", argvk[arg_idx + 1]);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            if (unlikely(ScanUintDefcapx(argvk[arg_idx + 2], &pc.splitpar_bound2) || (pc.splitpar_bound2 <= pc.splitpar_bound1))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --split-par argument '%s'.\n", argvk[arg_idx + 2]);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.dependency_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "et-all-var-ids") || strequal_k_unsafe(flagname_p2, "et-missing-var-ids")) {
          if (flagname_p2[3] == 'm') {
            if (unlikely(pc.varid_template_str)) {
              logerrputs("Error: --set-missing-var-ids cannot be used with --set-all-var-ids.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            pc.misc_flags |= kfMiscSetMissingVarIds;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(!VaridTemplateIsValid(argvk[arg_idx + 1], flagname_p))) {
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_template_str);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          if (!(pc.command_flags1 & kfCommand1Pmerge)) {
            pc.dependency_flags |= kfFilterPvarReq;
          }
        } else if (strequal_k_unsafe(flagname_p2, "et-invalid-haploid-missing")) {
          if (unlikely(!(make_plink2_flags & (kfMakeBed | kfMakePgen)))) {
            logerrputs("Error: --set-invalid-haploid-missing must be used with\n--make-[b]pgen/--make-bed.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (likely(!strcmp(cur_modif, "keep-dosage"))) {
              if (unlikely(make_plink2_flags & kfMakePgenEraseDosage)) {
                logerrputs("Error: --set-invalid-haploid-missing 'keep-dosage' modifier cannot be used with\n--make-[b]pgen erase-dosage.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              make_plink2_flags |= kfMakePlink2SetInvalidHaploidMissingKeepDosage;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --set-invalid-haploid-missing argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          make_plink2_flags |= kfMakePlink2SetInvalidHaploidMissing;
        } else if (strequal_k_unsafe(flagname_p2, "et-mixed-mt-missing")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1MakePlink2))) {
            logerrputs("Error: --set-mixed-mt-missing must be used with --make-[b]pgen/--make-bed.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(make_plink2_flags & kfMakePlink2SetMeMissing)) {
            // rather not specify order of operations here
            logerrputs("Error: --set-mixed-mt-missing cannot be used with --set-me-missing.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (likely(!strcmp(cur_modif, "keep-dosage"))) {
              if (unlikely(make_plink2_flags & kfMakePgenEraseDosage)) {
                logerrputs("Error: --set-mixed-mt-missing 'keep-dosage' modifier cannot be used with\n--make-[b]pgen erase-dosage.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              make_plink2_flags |= kfMakePlink2SetMixedMtMissingKeepDosage;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --set-mixed-mt-missing argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          make_plink2_flags |= kfMakePlink2SetMixedMtMissing;
        } else if (strequal_k_unsafe(flagname_p2, "ample")) {
          if (unlikely(load_params || (xload & (~(kfXloadOxGen | kfXloadOxBgen | kfXloadOxHaps | kfXloadOxLegend | kfXloadOxSample))))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(!(xload & (kfXloadOxGen | kfXloadOxBgen | kfXloadOxHaps)))) {
            logerrputs("Error: --sample must be used with --gen/--bgen/--data/--haps.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --sample filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(psamname, cur_fname, slen + 1);
          xload |= kfXloadOxSample;
        } else if (strequal_k_unsafe(flagname_p2, "heep")) {
          if (unlikely(chr_info.chrset_source)) {
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
          if (unlikely(pc.varid_exclude_snp)) {
            // problematic due to --window
            logerrputs("Error: --snp cannot be used with --exclude-snp.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_snp);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "nps")) {
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.snps_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "nps-only")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            if (likely(!strcmp(cur_modif, "just-acgt"))) {
              pc.filter_flags |= kfFilterSnpsOnlyJustAcgt;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --snps-only argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogSnpsOnly;
        } else if (strequal_k_unsafe(flagname_p2, "core") || strequal_k_unsafe(flagname_p2, "core-list")) {
          const uint32_t multi_input = (flagname_p2[4] == '-');
          if (multi_input) {
            if (unlikely(pc.command_flags1 & kfCommand1Score)) {
              logerrputs("Error: --score-list cannot be used with --score.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            pc.score_info.flags |= kfScoreMultiInput;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 12))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.score_info.input_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          uint32_t numeric_param_ct = 0;
          uint32_t score_cols[4];
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
              logerrputs("Note: --score[-list]'s 'variance-normalize' modifier has been renamed to the\nmore precise 'variance-standardize'.\n");
              pc.score_info.flags |= kfScoreVarianceStandardize;
            } else if (strequal_k(cur_modif, "se", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreSe;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreZs;
            } else if (strequal_k(cur_modif, "ignore-dup-ids", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreIgnoreDupIds;
            } else if (strequal_k(cur_modif, "list-variants", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreListVariants;
            } else if (strequal_k(cur_modif, "list-variants-zs", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreListVariants | kfScoreListVariantsZs;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (unlikely(pc.score_info.flags & kfScoreColAll)) {
                logerrputs("Error: Multiple --score[-list] cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0pheno1\0phenos\0nallele\0denom\0dosagesum\0scoreavgs\0scoresums\0", "score", kfScoreColMaybefid, multi_input? kfScoreListColDefault : kfScoreColDefault, 1, &pc.score_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              if (unlikely(multi_input && (pc.score_info.flags & (kfScoreColNallele | kfScoreColDenom | kfScoreColDosageSum)))) {
                logerrputs("Error: --score-list cannot be used with the nallele, denom, or dosagesum\ncolumn-sets.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            } else {
              if (unlikely(ScanPosintCappedx(cur_modif, kMaxLongLine / 2, &(score_cols[numeric_param_ct])))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --score[-list] argument '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (unlikely(numeric_param_ct == 3)) {
                logerrputs("Error: --score[-list] takes at most three numeric arguments.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              for (uint32_t uii = 0; uii != numeric_param_ct; ++uii) {
                if (unlikely(score_cols[uii] == score_cols[numeric_param_ct])) {
                  logerrputs("Error: Identical --score[-list] column indexes.\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
              }
              ++numeric_param_ct;
            }
          }
          if (unlikely((pc.score_info.flags & (kfScoreHeaderIgnore | kfScoreHeaderRead)) == (kfScoreHeaderIgnore | kfScoreHeaderRead))) {
            logerrputs("Error: --score[-list] 'header' and 'header-read' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          uint32_t model_flags_u = S_CAST(uint32_t, pc.score_info.flags & (kfScoreDominant | kfScoreRecessive | kfScoreCenter | kfScoreVarianceStandardize));
          if (unlikely(model_flags_u & (model_flags_u - 1))) {
            logerrputs("Error: --score[-list] 'dominant', 'recessive', 'center', and\n'variance-standardize' modifiers are mutually exclusive.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.score_info.flags & (kfScoreListVariants | kfScoreMultiInput)) == (kfScoreListVariants | kfScoreMultiInput))) {
            logerrputs("Error: 'list-variants[-zs]' modifier can only be used with --score, not\n--score-list.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (!(pc.score_info.flags & kfScoreColAll)) {
            pc.score_info.flags |= multi_input? kfScoreListColDefault : kfScoreColDefault;
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
            if (unlikely(pgl_malloc(col_idx_blen + 1, &new_buf))) {
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
          if (unlikely(pc.score_info.input_col_idx_range_list.name_ct)) {
            logerrputs("Error: --score-col-nums cannot be used when three numeric arguments are\nprovided to --score.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.score_info.input_col_idx_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          score_col_nums_present = 1;
        } else if (strequal_k_unsafe(flagname_p2, "plit-cat-pheno")) {
          uint32_t first_phenoname_idx = 1;
          for (; first_phenoname_idx <= param_ct; ++first_phenoname_idx) {
            const char* cur_modif = argvk[arg_idx + first_phenoname_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "omit-most", cur_modif_slen)) {
              pc.pheno_transform_flags |= kfPhenoTransformSplitCatOmitMost;
            } else if (strequal_k(cur_modif, "omit-last", cur_modif_slen)) {
              pc.pheno_transform_flags |= kfPhenoTransformSplitCatOmitLast;
            } else if (strequal_k(cur_modif, "covar-01", cur_modif_slen)) {
              pc.pheno_transform_flags |= kfPhenoTransformSplitCatCovar01;
            } else {
              break;
            }
          }
          if ((pc.pheno_transform_flags & (kfPhenoTransformSplitCatOmitMost | kfPhenoTransformSplitCatOmitLast)) == (kfPhenoTransformSplitCatOmitMost | kfPhenoTransformSplitCatOmitLast)) {
            logerrputs("Error: --split-cat-pheno 'omit-most' and 'omit-last' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (first_phenoname_idx <= param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + first_phenoname_idx]), flagname_p, param_ct + 1 - first_phenoname_idx, kMaxIdSlen - 1, &pc.split_cat_phenonames_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            // may as well verify that no phenotype name has an '=' in it
            const char* phenonames_iter = pc.split_cat_phenonames_flattened;
            do {
              const uint32_t cur_phenoname_slen = strlen(phenonames_iter);
              if (unlikely(memchr(phenonames_iter, '=', cur_phenoname_slen))) {
                logerrputs("Error: --split-cat-pheno phenotype names may not contain the '=' character.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              phenonames_iter = &(phenonames_iter[cur_phenoname_slen + 1]);
            } while (*phenonames_iter);
          } else if (unlikely(pc.pheno_transform_flags & kfPhenoTransformSplitCatCovar01)) {
            logerrputs("Error: --split-cat-pheno 'covar-01' modifier cannot be used without any\nphenotype names.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.pheno_transform_flags |= kfPhenoTransformSplitCat;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ort-vars")) {
          if (unlikely(!(((pc.command_flags1 & kfCommand1MakePlink2) && (make_plink2_flags & (kfMakeBed | kfMakeBim | kfMakePgen | kfMakePvar))) || (pc.command_flags1 & kfCommand1Pmerge)))) {
            logerrputs("Error: --sort-vars must be used with --make-[b]pgen/--make-bed,\n--make-just-{bim,pvar}, or dataset\nmerging.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* mode_str = argvk[arg_idx + 1];
            const char first_char_upcase_match = mode_str[0] & 0xdf;
            const uint32_t mode_slen = strlen(mode_str);
            if (((mode_slen == 1) && (first_char_upcase_match == 'N')) ||
                strequal_k(mode_str, "natural", mode_slen)) {
              pc.sort_vars_mode = kSortNatural;
            } else if (likely(((mode_slen == 1) && (first_char_upcase_match == 'A')) ||
                              strequal_k(mode_str, "ascii", mode_slen))) {
              pc.sort_vars_mode = kSortAscii;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a valid mode for --sort-vars.\n", mode_str);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.sort_vars_mode = kSortNatural;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ample-diff")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t param_idx = 1;

          // defer this since default depends on whether we're in 'pairwise'
          // mode
          uint32_t diff_cols_param_idx = 0;
          uint32_t is_file = 0;
          char sdiff_id_delim = '\0';
          for (; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (StrStartsWith(cur_modif, "id-delim=", cur_modif_slen)) {
              if (unlikely(sdiff_id_delim)) {
                logerrputs("Error: Multiple --sample-diff id-delim= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              sdiff_id_delim = ExtractCharParam(&(cur_modif[strlen("id-delim=")]));
              if (unlikely(!sdiff_id_delim)) {
                logerrputs("Error: --sample-diff id-delim= value must be a single character.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely((ctou32(sdiff_id_delim) <= ' ') || (sdiff_id_delim == '.'))) {
                logerrputs("Error: --sample-diff id-delim= value cannot be tab/space/newline, '.', or a\nnonprinting character.\n");
                goto main_ret_INVALID_CMDLINE;
              }
            } else if (StrStartsWith0(cur_modif, "dosage", cur_modif_slen)) {
              if (unlikely(pc.sdiff_info.dosage_hap_tol != kDosageMissing)) {
                logerrputs("Error: Multiple --sample-diff dosage modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              if (cur_modif[6] == '\0') {
                pc.sdiff_info.dosage_hap_tol = 0;
              } else {
                if (unlikely(cur_modif[6] != '=')) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --sample-diff argument '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                double dxx;
                if (unlikely((!ScantokDouble(&(cur_modif[7]), &dxx)) || (dxx < 0.0) || (dxx > (0.5 - kSmallEpsilon)))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Invalid --sample-diff argument '%s'.\n", cur_modif);
                  goto main_ret_INVALID_CMDLINE_WWA;
                }
                pc.sdiff_info.dosage_hap_tol = S_CAST(int32_t, dxx * ((1 + kSmallEpsilon) * kDosageMax));
              }
            } else if (strequal_k(cur_modif, "include-missing", cur_modif_slen)) {
              pc.sdiff_info.flags |= kfSdiffIncludeMissing;
            } else if (strequal_k(cur_modif, "pairwise", cur_modif_slen)) {
              pc.sdiff_info.flags |= kfSdiffPairwise;
            } else if (strequal_k(cur_modif, "counts-only", cur_modif_slen)) {
              pc.sdiff_info.flags |= kfSdiffCountsOnly;
            } else if (StrStartsWith(cur_modif, "fname-id-delim=", cur_modif_slen)) {
              if (unlikely(pc.sdiff_info.fname_id_delim)) {
                logerrputs("Error: Multiple --sample-diff fname-id-delim= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              pc.sdiff_info.fname_id_delim = ExtractCharParam(&(cur_modif[strlen("fname-id-delim=")]));
              if (unlikely(!pc.sdiff_info.fname_id_delim)) {
                logerrputs("Error: --sample-diff fname-id-delim= value must be a single character.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              if (unlikely((ctou32(pc.sdiff_info.fname_id_delim) <= ' ') || (pc.sdiff_info.fname_id_delim == '.'))) {
                logerrputs("Error: --sample-diff fname-id-delim= value cannot be tab/space/newline, '.', or\na nonprinting character.\n");
                goto main_ret_INVALID_CMDLINE;
              }
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.sdiff_info.flags |= kfSdiffZs;
            } else if (StrStartsWith0(cur_modif, "cols=", cur_modif_slen)) {
              if (unlikely(diff_cols_param_idx)) {
                logerrputs("Error: Multiple --sample-diff cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              diff_cols_param_idx = param_idx;
            } else if (StrStartsWith(cur_modif, "counts-cols=", cur_modif_slen)) {
              if (unlikely(pc.sdiff_info.flags & kfSdiffCountsColAll)) {
                logerrputs("Error: Multiple --sample-diff counts-cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[strlen("counts-cols=")]), "maybefid\0fid\0maybesid\0sid\0nobs\0nobsibs\0ibs0\0ibs1\0ibs2\0halfmiss\0diff\0", "sample-diff counts-cols=", kfSdiffCountsColMaybefid, kfSdiffCountsColDefault, 1, &pc.sdiff_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else if (StrStartsWith(cur_modif, "base=", cur_modif_slen)) {
              reterr = CmdlineAllocString(&(cur_modif[5]), argvk[arg_idx], kPglFnamesize - 11, &(pc.sdiff_info.first_id_or_fname));
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              pc.sdiff_info.flags |= kfSdiffOneBase;
              break;
            } else if (StrStartsWith(cur_modif, "ids=", cur_modif_slen)) {
              reterr = CmdlineAllocString(&(cur_modif[4]), argvk[arg_idx], kMaxIdSlen, &(pc.sdiff_info.first_id_or_fname));
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              break;
            } else if (likely(StrStartsWith(cur_modif, "file=", cur_modif_slen))) {
              reterr = AllocFname(&(cur_modif[5]), argvk[arg_idx], &pc.sdiff_info.first_id_or_fname);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
              is_file = 1;
              break;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --sample-diff argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (param_idx > param_ct) {
            logerrputs("Error: Invalid --sample-diff argument sequence (base=/id= must be\nsecond-to-last argument or earlier, or file= must be last argument).\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (diff_cols_param_idx) {
            reterr = ParseColDescriptor(&(argvk[arg_idx + diff_cols_param_idx][5]), "chrom\0pos\0ref\0alt\0maybeprovref\0provref\0maybefid\0fid\0id\0maybesid\0sid\0geno\0", "sample-diff cols=", kfSdiffColChrom, (pc.sdiff_info.flags & kfSdiffPairwise)? kfSdiffColPairwiseDefault : kfSdiffColDefault, 0, &pc.sdiff_info.flags);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            if (!(pc.sdiff_info.flags & kfSdiffColId)) {
              if (unlikely(pc.sdiff_info.flags & (kfSdiffColMaybefid | kfSdiffColFid | kfSdiffColMaybesid | kfSdiffColSid))) {
                logerrputs("Error: Invalid --sample-diff cols= set ('maybefid', 'fid', 'maybesid', and\n'sid' require 'id').\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
            }
          } else {
            pc.sdiff_info.flags |= (pc.sdiff_info.flags & kfSdiffPairwise)? kfSdiffColPairwiseDefault : kfSdiffColDefault;
          }
          if ((pc.sdiff_info.flags & (kfSdiffPairwise | kfSdiffCountsOnly)) == (kfSdiffPairwise | kfSdiffCountsOnly)) {
            logerrputs("Error: --sample-diff 'pairwise' and 'counts-only' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          const uint32_t other_id_ct = param_ct - param_idx;
          if (!other_id_ct) {
            if (unlikely(!is_file)) {
              logerrputs("Error: --sample-diff 'base='/'ids=' require 2 or more space-separated sample\nIDs.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (sdiff_id_delim) {
              logerrputs("Error: --sample-diff id-delim= modifier does not apply to file= mode.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          } else {
            if (unlikely(is_file)) {
              // this constraint is a bit arbitrary, but may as well have it
              // for consistency with base=/ids=
              logerrputs("Error: --sample-diff 'file=' must appear after all other modifiers.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            reterr = AllocAndFlatten(&(argvk[arg_idx + param_idx + 1]), flagname_p, other_id_ct, kPglFnamesize, &pc.sdiff_info.other_ids_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            if (sdiff_id_delim) {
              char* id_iter = pc.sdiff_info.first_id_or_fname;
              if (unlikely(ReplaceCharAdvChecked(sdiff_id_delim, '\t', &id_iter))) {
                logerrputs("Error: --sample-diff sample IDs cannot include tabs.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              id_iter = pc.sdiff_info.other_ids_flattened;
              for (uint32_t uii = 0; uii != other_id_ct; ++uii) {
                if (unlikely(ReplaceCharAdvChecked(sdiff_id_delim, '\t', &id_iter))) {
                  logerrputs("Error: --sample-diff sample IDs cannot include tabs.\n");
                  goto main_ret_INVALID_CMDLINE;
                }
                ++id_iter;
              }
            }
          }
          pc.sdiff_info.other_id_ct = other_id_ct;
          if (!(pc.sdiff_info.flags & kfSdiffCountsColAll)) {
            pc.sdiff_info.flags |= kfSdiffCountsColDefault;
          }
          if ((pc.sdiff_info.flags & (kfSdiffIncludeMissing | kfSdiffCountsColHalfmiss)) == kfSdiffCountsColHalfmiss) {
            pc.sdiff_info.flags &= ~kfSdiffCountsColHalfmiss;
          }
          if (!pc.sdiff_info.fname_id_delim) {
            pc.sdiff_info.fname_id_delim = '_';
          }
          pc.command_flags1 |= kfCommand1Sdiff;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "ample-counts")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.sample_counts_flags |= kfSampleCountsZs;
            } else if (likely(StrStartsWith(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(pc.sample_counts_flags & kfGenoCountsColAll)) {
                logerrputs("Error: Multiple --sample-counts cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0sex\0hom\0homref\0homalt\0homaltsnp\0het\0refalt\0het2alt\0hetsnp\0dipts\0ts\0diptv\0tv\0dipnonsnpsymb\0nonsnpsymb\0symbolic\0nonsnp\0dipsingle\0single\0haprefwfemaley\0hapref\0hapaltwfemaley\0hapalt\0missingwithfemaley\0missing\0", "sample-counts", kfSampleCountsColMaybefid, kfSampleCountsColDefault, 1, &pc.sample_counts_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --sample-counts argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if (!(pc.sample_counts_flags & kfSampleCountsColAll)) {
            pc.sample_counts_flags |= kfSampleCountsColDefault;
          }
          pc.command_flags1 |= kfCommand1SampleCounts;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "trict-sid0")) {
          pc.misc_flags |= kfMiscStrictSid0;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "ample-inner-join")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Pmerge))) {
            logerrputs("Error: --sample-inner-join must be used with --pmerge[-list].\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pmerge_info.flags & kfPmergeSids)) {
            logerrputs("Error: --sample-inner-join cannot be used with --merge-sids.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pmerge_info.flags |= kfPmergeSampleInnerJoin;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "elect-sid-representatives")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          pc.select_sid_missingness_mode = kSelectSidMissingnessHc;
          uint32_t missingness_modes_seen = 0;
          uint32_t tiebreak_modes_seen = 0;
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "hc", cur_modif_slen) ||
                strequal_k(cur_modif, "hardcall", cur_modif_slen)) {
              // default
              ++missingness_modes_seen;
            } else if (strequal_k(cur_modif, "hh-missing", cur_modif_slen)) {
              pc.select_sid_missingness_mode = kSelectSidMissingnessHhMissing;
              ++missingness_modes_seen;
            } else if (strequal_k(cur_modif, "dosage", cur_modif_slen)) {
              pc.select_sid_missingness_mode = kSelectSidMissingnessDosage;
              ++missingness_modes_seen;
            } else if (strequal_k(cur_modif, "sid-only", cur_modif_slen)) {
              pc.select_sid_missingness_mode = kSelectSidMissingness0;
              ++missingness_modes_seen;
            } else if (strequal_k(cur_modif, "first", cur_modif_slen)) {
              pc.select_sid_tiebreak_mode = kSelectSidTiebreakFirst;
              ++tiebreak_modes_seen;
            } else if (strequal_k(cur_modif, "first-ascii", cur_modif_slen)) {
              pc.select_sid_tiebreak_mode = kSelectSidTiebreakFirstAscii;
              ++tiebreak_modes_seen;
            } else if (strequal_k(cur_modif, "last", cur_modif_slen)) {
              pc.select_sid_tiebreak_mode = kSelectSidTiebreakLast;
              ++tiebreak_modes_seen;
            } else if (strequal_k(cur_modif, "last-ascii", cur_modif_slen)) {
              pc.select_sid_tiebreak_mode = kSelectSidTiebreakLastAscii;
              ++tiebreak_modes_seen;
            } else if (likely(strequal_k(cur_modif, "parents-only", cur_modif_slen))) {
              pc.misc_flags |= kfMiscSelectSidParentsOnly;
            } else {
              logerrputs("Error: Invalid --select-sid-representatives argument sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          if (unlikely(missingness_modes_seen > 1)) {
            logerrputs("Error: Multiple --select-sid-representatives missingness-mode modifiers.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(tiebreak_modes_seen != 1)) {
            if (tiebreak_modes_seen) {
              logerrputs("Error: Multiple --select-sid-representatives tiebreak-mode modifiers.\n");
              goto main_ret_INVALID_CMDLINE;
            }
            logerrputs("Error: --select-sid-representatives requires a tiebreak mode.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.filter_flags |= kfFilterPsamReq | kfFilterSelectSidRepresentatives;
          if (pc.select_sid_missingness_mode != kSelectSidMissingness0) {
            pc.filter_flags |= kfFilterAllReq;
          }
        } else if (strequal_k_unsafe(flagname_p2, "et-me-missing")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1MakePlink2))) {
            logerrputs("Error: --set-me-missing must be used with --make-[b]pgen/--make-bed.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(make_plink2_flags & kfMakePlink2SetInvalidHaploidMissing)) {
            logerrputs("Error: --set-me-missing cannot be used with --set-invalid-haploid-missing.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.flip_info.subset_fname)) {
            logerrputs("Error: --set-me-missing cannot be used with --flip-subset.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          make_plink2_flags |= kfMakePlink2SetMeMissing;
          goto main_param_zero;
        } else if (unlikely(!strequal_k_unsafe(flagname_p2, "ilent"))) {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 't':
        if (strequal_k_unsafe(flagname_p2, "hreads")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(ScanPosintDefcapx(argvk[arg_idx + 1], &pc.max_thread_ct))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --threads argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (pc.max_thread_ct > kMaxThreads) {
            logprintfww("Note: Reducing --threads argument to %u.  (If this is not large enough,\nrecompile with a larger kMaxThreads setting.)\n", kMaxThreads);
            pc.max_thread_ct = kMaxThreads;
          } else if (known_procs == -1) {
            // trigger BLAS/LAPACK warning?
            known_procs = 0;
          }
        } else if (strequal_k_unsafe(flagname_p2, "o")) {
          if (unlikely(chr_info.is_include_stack || (pc.load_filter_log_flags & kfLoadFilterLogNotChr))) {
            logerrputs("Error: --from/--to cannot be used with --autosome[-par] or --[not-]chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, &pc.varid_to);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq | kfFilterNoSplitChr;
        } else if (strequal_k_unsafe(flagname_p2, "o-bp") || strequal_k_unsafe(flagname_p2, "o-kb") || strequal_k_unsafe(flagname_p2, "o-mb")) {
          if (unlikely(!CmdlineSingleChr(&chr_info, pc.load_filter_log_flags))) {
            logerrputs("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb must be used with --chr, and only\none chromosome.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.load_filter_log_flags & kfLoadFilterLogNotChr))) {
            logerrputs("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb cannot be used with --not-chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.to_bp != -1)) {
            logerrputs("Error: Multiple --to-bp/-kb/-mb values.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (unlikely(!ScantokDouble(cur_modif, &dxx))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --to-bp/-kb/-mb argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          const char unit_char = flagname_p2[2];
          if (unit_char == 'k') {
            dxx *= 1000;
          } else if (unit_char == 'm') {
            dxx *= 1000000;
          }
          if (unlikely(dxx < 0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: --to-bp/-kb/-mb argument '%s' too small.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_2A;
          } else if (dxx >= 2147483646) {
            pc.to_bp = 0x7ffffffe;
          } else {
            // round down
            pc.to_bp = S_CAST(int32_t, dxx * (1 + kSmallEpsilon));
          }
          if (unlikely(pc.from_bp > pc.to_bp)) {
            // (if we do permit this, rounding must be postponed)
            logerrputs("Error: --to-bp/-kb/-mb argument is smaller than --from-bp/-kb/-mb argument.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "hin")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(!ScantokDouble(cur_modif, &pc.thin_keep_prob))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --thin variant retention probability '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(pc.thin_keep_prob < (0.5 / 4294967296.0))) {
            logerrputs("Error: --thin variant retention probability too small.\n");
            goto main_ret_INVALID_CMDLINE_A;
          } else if (unlikely(pc.thin_keep_prob >= (4294967295.5 / 4294967296.0))) {
            uint32_t uii;
            if (ScanUintDefcapx(cur_modif, &uii)) {
              logerrputs("Error: --thin variant retention probability too large.\n");
            } else {
              // VCFtools --thin = --bp-space...
              logerrputs("Error: --thin variant retention probability too large.  (Did you mean\n--bp-space?)\n");
            }
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "hin-count")) {
          if (unlikely(pc.thin_keep_prob != 1.0)) {
            logerrputs("Error: --thin cannot be used with --thin-count.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(ScanUintDefcapx(cur_modif, &pc.thin_keep_ct) || (!pc.thin_keep_ct))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --thin-count argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "hin-indiv")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(!ScantokDouble(cur_modif, &pc.thin_keep_sample_prob))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --thin-indiv sample retention probability '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(pc.thin_keep_sample_prob < (0.5 / 4294967296.0))) {
            logerrputs("Error: --thin-indiv sample retention probability too small.\n");
            goto main_ret_INVALID_CMDLINE_A;
          } else if (unlikely(pc.thin_keep_sample_prob >= (4294967295.5 / 4294967296.0))) {
            logerrputs("Error: --thin-indiv sample retention probability too large.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "hin-indiv-count")) {
          if (unlikely(pc.thin_keep_sample_prob != 1.0)) {
            logerrputs("Error: --thin-indiv cannot be used with --thin-indiv-count.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(ScanUintDefcapx(cur_modif, &pc.thin_keep_sample_ct) || (!pc.thin_keep_sample_ct))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --thin-indiv-count argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "file")) {
          if (unlikely((load_params & (~kfLoadParamsPsam)) || xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* fname_prefix = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname_prefix);
          if (unlikely(slen > kPglFnamesize - 6)) {
            logerrputs("Error: --tfile argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          if (!load_params) {
            // allow --tfam (aliased to --psam) to override this
            char* prefix_end = memcpya(psamname, fname_prefix, slen);
            snprintf(prefix_end, 6, ".tfam");
            load_params = kfLoadParamsPsam;
          }
          char* prefix_end = memcpya(pgenname, fname_prefix, slen);
          snprintf(prefix_end, 6, ".tped");
          xload = kfXloadTped;
        } else if (strequal_k_unsafe(flagname_p2, "ped")) {
          if (unlikely((load_params & (~kfLoadParamsPsam)) || (xload & (~kfXloadTped)))) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* fname = argvk[arg_idx + 1];
          const uint32_t slen = strlen(fname);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --tped argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, fname, slen + 1);
          xload = kfXloadTped;
        } else if (likely(strequal_k_unsafe(flagname_p2, "ests"))) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Glm))) {
            logerrputs("Error: --tests must be used with --glm.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.glm_info.flags & kfGlmResidualizeMask)) {
            logerrputs("Error: --tests cannot be used with --glm {cc,firth,qt}-residualize.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.gwas_ssf_info.flags & kfGwasSsfDeleteOrigGlm)) {
            logerrputs("Error: --tests doesn't make sense with \"--gwas-ssf delete-orig-glm\" unless\ntest= is specified for --gwas-ssf.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if ((param_ct == 1) && (!strcmp(argvk[arg_idx + 1], "all"))) {
            pc.glm_info.flags |= kfGlmTestsAll;
          } else {
            reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.glm_info.tests_range_list);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'u':
        if (strequal_k_unsafe(flagname_p2, "pdate-sex")) {
          if (unlikely(pc.update_sample_ids_fname)) {
            logerrputs("Error: --update-sex cannot be used with --update-ids.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.update_sex_info.fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "male0", cur_modif_slen)) {
              pc.update_sex_info.flags |= kfUpdateSexMale0;
            } else if (StrStartsWith(cur_modif, "col-num=", cur_modif_slen)) {
              const char* col_num_start = &(cur_modif[strlen("col-num=")]);
              if (unlikely(ScanPosintDefcapx(col_num_start, &pc.update_sex_info.col_num) || (pc.update_sex_info.col_num == 1))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --update-sex col-num= argument '%s'.\n", col_num_start);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            } else if (likely(param_ct == 2)) {
              // only one extra argument, try to interpret it the plink 1.9 way
              // but print a warning
              if (unlikely(ScanPosintDefcapx(cur_modif, &pc.update_sex_info.col_num))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --update-sex argument '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              logerrputs("Warning: --update-sex unlabeled column argument is now deprecated.  Use\n'col-num=' instead (and add 2 to the value).\n");
              pc.update_sex_info.col_num += 2;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --update-sex argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "pdate-chr")) {
          if (unlikely(pc.sort_vars_mode <= kSortNone)) {
            logerrputs("Error: --update-chr must be used with --sort-vars.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(chr_info.is_include_stack || (pc.load_filter_log_flags & kfLoadFilterLogNotChr))) {
            // too confusing (does the chromosome filter apply before?  after?
            // both?)
            logerrputs("Error: --update-chr cannot be used with --autosome[-par] or --[not-]chr.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.recover_var_ids_fname)) {
            logerrputs("Error: --update-chr cannot be used with --recover-var-ids or --update-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.rename_chrs_fname)) {
            logerrputs("Error: --update-chr cannot be used with --rename-chrs.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = Alloc2col(&(argvk[arg_idx + 1]), flagname_p2, param_ct, &pc.update_chr_flag);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "pdate-map")) {
          if (unlikely(pc.recover_var_ids_fname)) {
            logerrputs("Error: --update-map cannot be used with --recover-var-ids or --update-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = Alloc2col(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.update_map_flag);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "pdate-name")) {
          if (unlikely(pc.recover_var_ids_fname)) {
            logerrputs("Error: --update-name cannot be used with --recover-var-ids.\n");
          }
          if (unlikely(pc.update_chr_flag)) {
            logerrputs("Error: --update-chr cannot be used with --recover-var-ids or --update-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.update_map_flag)) {
            logerrputs("Error: --update-map cannot be used with --recover-var-ids or --update-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = Alloc2col(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.update_name_flag);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "pdate-ids")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.update_sample_ids_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "pdate-parents")) {
          if (unlikely(pc.update_sample_ids_fname)) {
            logerrputs("Error: --update-parents cannot be used with --update-ids.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.update_parental_ids_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (likely(strequal_k_unsafe(flagname_p2, "pdate-alleles"))) {
          // Prevent the most confusing order-of-operations dependencies.  (May
          // want to prevent more combinations later.)
          if (unlikely(pc.varid_template_str)) {
            logerrputs("Error: --update-alleles cannot be used with --set-missing-var-ids or\n--set-all-var-ids.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 3))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "allow-mismatch", cur_modif_slen)) {
              pc.update_alleles_info.flags |= kfUpdateAllelesAllowMismatch;
            } else if (strequal_k(cur_modif, "strict-missing", cur_modif_slen)) {
              pc.update_alleles_info.flags |= kfUpdateAllelesStrictMissing;
            } else {
              if (unlikely(pc.update_alleles_info.fname)) {
                logerrputs("Error: Invalid --update-alleles argument sequence.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              reterr = AllocFname(cur_modif, flagname_p, &pc.update_alleles_info.fname);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            }
          }
          if (unlikely(!pc.update_alleles_info.fname)) {
            logerrputs("Error: Invalid --update-alleles argument sequence.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.dependency_flags |= kfFilterPvarReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'v':
        if (strequal_k_unsafe(flagname_p2, "ar-min-qual")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(ScanFloatAllowInf(argvk[arg_idx + 1], &pc.var_min_qual) || (pc.var_min_qual < S_CAST(float, 0.0)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --var-min-qual argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.var_min_qual *= S_CAST(float, 1 - kBigEpsilon);
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogVarMinQual;
        } else if (strequal_k_unsafe(flagname_p2, "ar-filter")) {
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.var_filter_exceptions_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.filter_flags |= kfFilterPvarReq;
          pc.load_filter_log_flags |= kfLoadFilterLogVarFilter;
        } else if (strequal_k_unsafe(flagname_p2, "cf")) {
          // permit accompanying .fam/.psam
          // IIDs must match VCF sample line order
          if (unlikely((load_params & (~kfLoadParamsPsam)) || xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct == 2) {
            const char* cur_modif = argvk[arg_idx + 2];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            // tolerate vcf-dosage= as well, so it's possible to use the same
            // pattern as --export
            if (unlikely((!StrStartsWith(cur_modif, "dosage=", cur_modif_slen)) && (!StrStartsWith(cur_modif, "vcf-dosage=", cur_modif_slen)))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --vcf argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            const char* dosage_field_start = (cur_modif[0] == 'v')? &(cur_modif[11]) : &(cur_modif[7]);
            reterr = CmdlineAllocString(dosage_field_start, argvk[arg_idx], 4095, &vcf_dosage_import_field);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
            const uint32_t vcf_dosage_import_field_slen = strlen(vcf_dosage_import_field);
            if (unlikely(!((strequal_k(vcf_dosage_import_field, "GP-force", vcf_dosage_import_field_slen)) || IsAlphanumeric(vcf_dosage_import_field)))) {
              // special case: it is reasonble for a user who just used
              // "--export vcf vcf-dosage=DS-force" or DS-only to try DS-force
              // or DS-only here.
              if ((strequal_k(vcf_dosage_import_field, "DS-force", vcf_dosage_import_field_slen)) || (strequal_k(vcf_dosage_import_field, "DS-only", vcf_dosage_import_field_slen))) {
                logerrputs("Error: --vcf dosage= argument is not alphanumeric.  (It should contain just the\nname of the dosage field, i.e. 'DS' instead of 'DS-force'/'DS-only'.)\n");
              } else {
                logerrputs("Error: --vcf dosage= argument is not alphanumeric.\n");
              }
              goto main_ret_INVALID_CMDLINE;
            }
            if (unlikely(strequal_k(vcf_dosage_import_field, "GT", vcf_dosage_import_field_slen))) {
              logerrputs("Error: --vcf dosage= argument cannot be 'GT'.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          }
          const char* cur_modif = argvk[arg_idx + 1];
          const uint32_t slen = strlen(cur_modif);
          if (unlikely(slen > kPglFnamesize - 1)) {
            logerrputs("Error: --vcf filename too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          memcpy(pgenname, cur_modif, slen + 1);
          xload = kfXloadVcf;
        } else if (unlikely(strequal_k_unsafe(flagname_p2, "cf-min-gp"))) {
          logerrputs("Error: --vcf-min-gp is no longer supported.  Use --import-dosage-certainty\ninstead.\n");
          goto main_ret_INVALID_CMDLINE_A;
        } else if (strequal_k_unsafe(flagname_p2, "cf-min-gq") || strequal_k_unsafe(flagname_p2, "cf-min-dp")) {
          if (unlikely(!(xload & (kfXloadVcf | kfXloadBcf)))) {
            logerrprintf("Error: --%s must be used with --vcf/--bcf.\n", flagname_p);
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t uii;
          if (unlikely(ScanUintDefcapx(cur_modif, &uii))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s argument '%s'.\n", flagname_p, cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (flagname_p2[7] == 'g') {
            vcf_min_gq = uii;
          } else {
            vcf_min_dp = uii;
            if (vcf_max_dp < vcf_min_dp) {
              logerrputs("Error: --vcf-min-dp value cannot be larger than --vcf-max-dp value.\n");
              goto main_ret_INVALID_CMDLINE;
            }
          }
        } else if (strequal_k_unsafe(flagname_p2, "cf-max-dp")) {
          if (unlikely(!(xload & (kfXloadVcf | kfXloadBcf)))) {
            logerrputs("Error: --vcf-max-dp must be used with --vcf/--bcf.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          uint32_t uii;
          if (unlikely(ScanUintDefcapx(cur_modif, &uii))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --vcf-max-dp argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          vcf_max_dp = uii;
        } else if (strequal_k_unsafe(flagname_p2, "cf-half-call")) {
          if (unlikely(!(xload & (kfXloadVcf | kfXloadBcf)))) {
            logerrputs("Error: --vcf-half-call must be used with --vcf/--bcf.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* mode_str = argvk[arg_idx + 1];
          const char first_char_upcase_match = mode_str[0] & 0xdf;
          const uint32_t mode_slen = strlen(mode_str);
          if (((mode_slen == 1) && (first_char_upcase_match == 'H')) ||
              strequal_k(mode_str, "haploid", mode_slen)) {
            vcf_half_call = kVcfHalfCallHaploid;
          } else if (((mode_slen == 1) && (first_char_upcase_match == 'M')) ||
                     strequal_k(mode_str, "missing", mode_slen)) {
            vcf_half_call = kVcfHalfCallMissing;
          } else if (((mode_slen == 1) && (first_char_upcase_match == 'E')) ||
                     strequal_k(mode_str, "error", mode_slen)) {
            vcf_half_call = kVcfHalfCallError;
          } else if (likely(((mode_slen == 1) && (first_char_upcase_match == 'R')) ||
                            strequal_k(mode_str, "reference", mode_slen))) {
            vcf_half_call = kVcfHalfCallReference;
          } else {
            snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a valid mode for --vcf-half-call.\n", mode_str);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "cf-require-gt")) {
          if (unlikely(!(xload & (kfXloadVcf | kfXloadBcf)))) {
            logerrputs("Error: --vcf-require-gt must be used with --vcf/--bcf.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          pc.load_filter_log_flags |= kfLoadFilterLogVcfRequireGt;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "cf-ref-n-missing")) {
          if (unlikely(!(xload & (kfXloadVcf | kfXloadBcf)))) {
            logerrputs("Error: --vcf-ref-n-missing must be used with --vcf/--bcf.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          import_flags |= kfImportVcfRefNMissing;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "cf-allow-no-nonvar")) {
          // Don't bother with --bcf, since GATK GenotypeGVCFs only generates
          // VCFs.
          if (unlikely(!(xload & kfXloadVcf))) {
            logerrputs("Error: --vcf-allow-no-nonvar must be used with --vcf.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          import_flags |= kfImportVcfAllowNoNonvar;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "if")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Glm))) {
            logerrputs("Error: --vif must be used with --glm/--epistasis.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(!ScantokDouble(cur_modif, &pc.vif_thresh))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --glm/--epistasis VIF threshold '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(pc.vif_thresh < 1.0)) {
            snprintf(g_logbuf, kLogbufSize, "Error: --glm/--epistasis VIF threshold '%s' too small (must be >= 1).\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ariance-standardize")) {
          if (unlikely(pc.pheno_transform_flags & kfPhenoTransformVstdCovar)) {
            logerrputs("Error: --variance-standardize cannot be used with --covar-variance-standardize.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (param_ct) {
            reterr = AllocAndFlatten(&(argvk[arg_idx + 1]), flagname_p, param_ct, 0x7fffffff, &pc.vstd_flattened);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.pheno_transform_flags |= kfPhenoTransformVstdAll;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ar-id-multi") ||
                   strequal_k_unsafe(flagname_p2, "ar-id-multi-nonsnp")) {
          if (unlikely(!pc.varid_template_str)) {
            logerrputs("Error: --var-id-multi[-nonsnp] must be used with --set-missing-var-ids or\n--set-all-var-ids.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(!VaridTemplateIsValid(argvk[arg_idx + 1], flagname_p))) {
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = CmdlineAllocString(argvk[arg_idx + 1], argvk[arg_idx], kMaxIdSlen, flagname_p2[11]? (&pc.varid_multi_nonsnp_template_str) : (&pc.varid_multi_template_str));
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ariant-score")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 4))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.vscore_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          uint32_t explicit_cols = 0;
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.vscore_flags |= kfVscoreZs;
            } else if (strequal_k(cur_modif, "single-prec", cur_modif_slen)) {
              pc.vscore_flags |= kfVscoreSinglePrec;
            } else if (strequal_k(cur_modif, "bin", cur_modif_slen)) {
              // 'zs' and 'bin' *are* allowed together, since the variant-ID
              // (.vscore.vars) file accompanying the .vscore.bin can be large
              // enough to deserve compression.
              pc.vscore_flags |= kfVscoreBin;
            } else if (strequal_k(cur_modif, "bin4", cur_modif_slen)) {
              pc.vscore_flags |= kfVscoreBin4;
            } else if (likely(StrStartsWith0(cur_modif, "cols=", cur_modif_slen))) {
              if (unlikely(explicit_cols)) {
                logerrputs("Error: Multiple --variant-score cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              explicit_cols = 1;
              reterr = ParseColDescriptor(&(cur_modif[5]), "chrom\0pos\0ref\0alt1\0alt\0maybeprovref\0provref\0altfreq\0nmiss\0nobs\0", "variant-score", kfVscoreColChrom, kfVscoreColDefault, 0, &pc.vscore_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --variant-score argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          if ((pc.vscore_flags & kfVscoreBin) && (pc.vscore_flags & (kfVscoreSinglePrec | kfVscoreBin4))) {
            logerrputs("Error: --variant-score 'bin' modifier cannot be used with 'bin4' or\n'single-prec'.\n");
            goto main_ret_INVALID_CMDLINE;
          }
          if (!explicit_cols) {
            pc.vscore_flags |= kfVscoreColDefault;
          } else if (pc.vscore_flags & (kfVscoreBin | kfVscoreBin4)) {
            logerrputs("Error: --variant-score 'bin'/'bin4' and 'cols=' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.command_flags1 |= kfCommand1Vscore;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "score-col-nums")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Vscore))) {
            logerrputs("Error: --vscore-col-nums must be used with --variant-score.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.vscore_col_idx_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "ariant-inner-join")) {
          if (unlikely(!(pc.command_flags1 & kfCommand1Pmerge))) {
            logerrputs("Error: --variant-inner-join must be used with --pmerge[-list].\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pmerge_info.flags |= kfPmergeVariantInnerJoin;
          goto main_param_zero;
        } else if (likely(strequal_k_unsafe(flagname_p2, "alidate"))) {
          pc.command_flags1 |= kfCommand1Validate;
          pc.dependency_flags |= kfFilterAllReq;
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'w':
        if (strequal_k_unsafe(flagname_p2, "rite-snplist")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.misc_flags |= kfMiscWriteSnplistZs;
            } else if (likely(strequal_k(cur_modif, "allow-dups", cur_modif_slen))) {
              pc.misc_flags |= kfMiscWriteSnplistAllowDups;
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --write-snplist argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          pc.command_flags1 |= kfCommand1WriteSnplist;
          pc.dependency_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "rite-samples")) {
          if (unlikely((param_ct == 1) && (!strcmp(argvk[arg_idx + 1], "noheader")))) {
            logerrputs("Error: --write-samples 'noheader' modifier retired.  Use --no-id-header\ninstead.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.command_flags1 |= kfCommand1WriteSamples;
          pc.dependency_flags |= kfFilterPsamReq;
          goto main_param_zero;
        } else if (strequal_k_unsafe(flagname_p2, "indow")) {
          if (unlikely(!(pc.varid_snp || pc.varid_exclude_snp))) {
            logerrputs("Error: --window must be used with --snp or --exclude-snp.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          double dxx;
          if (unlikely((!ScantokDouble(cur_modif, &dxx)) || (dxx < 0))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --window argument '%s'.\n", cur_modif);
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
          if (unlikely(pc.misc_flags & kfMiscCatPhenoFamily)) {
            logerrputs("Error: --within cannot be used with --family.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (unlikely((cur_modif_slen == 7) && StrStartsWithUnsafe(cur_modif, "keep-") && MatchUpperK(&(cur_modif[5]), "NA"))) {
              logerrputs("Error: --within's keep-NA modifier has been retired.  Rename that category in\nthe input file if you wish to keep it.\n");
              goto main_ret_INVALID_CMDLINE_A;
            }
            if (param_idx == 1) {
              reterr = AllocFname(cur_modif, flagname_p, &pc.within_fname);
            } else {
              if (unlikely(IsReservedPhenoName(cur_modif, cur_modif_slen))) {
                snprintf(g_logbuf, kLogbufSize, "Error: '%s' cannot be used as a categorical phenotype name.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_2A;
              }
              reterr = CmdlineAllocString(cur_modif, argvk[arg_idx], kMaxIdSlen, &pc.catpheno_name);
            }
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "rite-covar")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* cur_modif = argvk[arg_idx + 1];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (likely(StrStartsWith0(cur_modif, "cols=", cur_modif_slen))) {
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0maybeparents\0parents\0sex\0pheno1\0phenos\0", "write-covar", kfWriteCovarColMaybefid, kfWriteCovarColDefault, 0, &pc.write_covar_flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --write-covar argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          } else {
            pc.write_covar_flags |= kfWriteCovarColDefault;
          }
          pc.command_flags1 |= kfCommand1WriteCovar;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (likely(strequal_k_unsafe(flagname_p2, "arning-errcode"))) {
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
        if (likely(strequal_k_unsafe(flagname_p2, "chr-model"))) {
          if (unlikely(!(pc.command_flags1 & (kfCommand1Glm | kfCommand1Score | kfCommand1Vscore)))) {
            logerrputs("Error: --xchr-model must be used with --glm, --score[-list], or\n--variant-score.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // quasi-bugfix (18 Sep 2021): nothing wrong with combining
          // --xchr-model with e.g. --glm genotypic, it gets ignored by --glm
          // but can still apply to other commands.
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          pc.xchr_model = ctou32(ExtractCharParam(cur_modif)) - 48;
          if (unlikely(pc.xchr_model > 2)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --xchr-model argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'y':
        if (likely(strequal_k_unsafe(flagname_p2, "-nosex-missing-stats"))) {
          pc.misc_flags |= kfMiscYNosexMissingStats;
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'z':
        if (strequal_k_unsafe(flagname_p2, "st-level")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          // Note that, in zstd 1.3.4, multithreaded compression is
          // nondeterministic unless level 5+ is explicitly requested, see
          //   https://github.com/facebook/zstd/issues/1077 .
          // I've postponed the decision on whether this sort of nondeterminism
          // is acceptable in plink2 for now (by reverting to 1.3.3).
          if (unlikely(ScanPosintCappedx(cur_modif, 22, &g_zst_level))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --zst-level argument '%s'.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else if (likely(strequal_k_unsafe(flagname_p2, "ero-cluster"))) {
          if (unlikely(!(pc.command_flags1 & kfCommand1MakePlink2))) {
            logerrputs("Error: --zero-cluster must be used with --make-[b]pgen/--make-bed.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely((pc.geno_thresh != 1.0) ||
                       (pc.mind_thresh != 1.0) ||
                       (pc.hwe_ln_thresh != -DBL_MAX) ||
                       (pc.min_maf != 0.0) ||
                       (pc.max_maf != 1.0) ||
                       pc.min_allele_ddosage ||
                       (pc.max_allele_ddosage != (~0LLU)) ||
                       (pc.mach_r2_max != 0.0) ||
                       (pc.minimac3_r2_max != 0.0))) {
            logerrputs("Error: --zero-cluster cannot be used with genotype-based filters.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.zero_cluster_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          if (param_ct == 2) {
            reterr = CmdlineAllocString(argvk[arg_idx + 2], argvk[arg_idx], kMaxIdSlen, &pc.zero_cluster_phenoname);
            if (unlikely(reterr)) {
              goto main_ret_1;
            }
          }
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      default:
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      main_param_zero:
        if (unlikely(param_ct)) {
          snprintf(g_logbuf, kLogbufSize, "Error: You provided %u argument%s to --%s, which has no parameters.\n", param_ct, (param_ct == 1)? "" : "s", flagname_p);
          goto main_ret_INVALID_CMDLINE_WWA;
        }
      }
    } while ((++cur_flag_idx) < flag_ct);
    if (!outname_end) {
      outname_end = &(outname[6]);
    } else if (!allow_misleading_out_arg) {
      char* outname_last_dot = Memrchr(outname, '.', outname_end - outname);
      if (outname_last_dot) {
        char* outname_ext = &(outname_last_dot[1]);
        const uint32_t outname_ext_slen = outname_end - outname_ext;
        if (strequal_k(outname_ext, "bcf", outname_ext_slen) ||
            strequal_k(outname_ext, "bed", outname_ext_slen) ||
            strequal_k(outname_ext, "bgen", outname_ext_slen) ||
            strequal_k(outname_ext, "gz", outname_ext_slen) ||
            strequal_k(outname_ext, "ped", outname_ext_slen) ||
            strequal_k(outname_ext, "pgen", outname_ext_slen) ||
            strequal_k(outname_ext, "vcf", outname_ext_slen)) {
          logerrputs("Error: --out argument is supposed to be a filename *prefix*; the value you've\nprovided looks like it is intended to be an entire filename.\nIf the current --out argument is really what you want, use\n--allow-misleading-out-arg to suppress this error.\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
      }
    }

    pc.dependency_flags |= pc.filter_flags;
    const uint32_t skip_main = (!pc.command_flags1) && (!(xload & (kfXloadVcf | kfXloadBcf | kfXloadOxBgen | kfXloadOxHaps | kfXloadOxSample | kfXloadEigGeno | kfXloadPlink1Dosage | kfXloadGenDummy | kfXloadPed | kfXloadTped)));
    const uint32_t batch_job = (adjust_file_info.fname != nullptr) || (pc.gwas_ssf_info.fname != nullptr) || (pc.gwas_ssf_info.list_fname != nullptr);
    if (skip_main && (!batch_job)) {
      // add command_flags2 when needed
      goto main_ret_NULL_CALC;
    }
    if (unlikely(!(load_params || xload || batch_job || pmerge_info.list_fname))) {
      logerrputs("Error: No input dataset.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((xload & kfXloadMap) && (!(xload & (kfXloadPed | kfXloadPlink1Dosage))))) {
      logerrputs("Error: --map must be used with --import-dosage or --ped.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((xload & kfXloadOxGen) && (!(xload & kfXloadOxSample)))) {
      // could permit .fam/.psam, but unless Oxford software supports that mode
      // it's pointless
      logerrputs("Error: --gen must be used with --sample or --data.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((xload & kfXloadEigGeno) && ((xload & (kfXloadEigInd | kfXloadEigSnp)) != (kfXloadEigInd | kfXloadEigSnp)))) {
      logerrputs("Error: --eiggeno must be used with either --eigfile, or --eigind + --eigsnp.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((pc.sample_sort_mode != kSort0) && (!(pc.command_flags1 & (kfCommand1MakePlink2 | kfCommand1WriteCovar | kfCommand1Pmerge))))) {
      logerrputs("Error: --indiv-sort must be used with --make-[b]pgen/--make-bed/--write-covar\nor dataset merging.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    // may as well permit merge here
    if (unlikely((make_plink2_flags & (kfMakePlink2MMask | kfMakePlink2TrimAlts | kfMakePgenErasePhase | kfMakePgenEraseDosage)) && (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Pmerge))))) {
      logerrputs("Error: When the 'multiallelics=', 'trim-alts', and/or 'erase-...' modifier is\npresent, --make-bed/--make-[b]pgen cannot be combined with other commands.\n(Other filters are fine.)\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (make_plink2_flags & kfMakePlink2MMask) {
      if (unlikely(make_plink2_flags & kfMakePlink2TrimAlts)) {
        logerrputs("Error: --make-bed/--make-[b]pgen 'multiallelics=' cannot be used with\n'trim-alts'.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely((pc.misc_flags & kfMiscMajRef) || pc.ref_allele_flag || pc.alt_allele_flag || (pc.fa_flags & kfFaRefFrom))) {
        logerrputs("Error: When the 'multiallelics=' modifier is present, --make-bed/--make-[b]pgen\ncannot be used with a flag which alters REF/ALT allele settings.\n");
        goto main_ret_INVALID_CMDLINE;
      }
    }
    if (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Pmerge))) {
      if (unlikely(pc.sort_vars_mode > kSortNone)) {
        logerrputs("Error: --sort-vars must be used with --make-[b]pgen/--make-bed and no other\nnon-merge commands.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(make_plink2_flags & (kfMakePlink2SetInvalidHaploidMissing | kfMakePlink2SetMixedMtMissing | kfMakePlink2SetMeMissing | kfMakePlink2FillMissingWithRef))) {
        logerrputs("Errpr: --set-invalid-haploid-missing/--set-mixed-mt-missing/--set-me-missing/\n--fill-missing-with-ref must be used with --make-[b]pgen/--make-bed and no\nother non-merge commands.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(pc.flip_info.subset_fname)) {
        logerrputs("Error: --flip-subset must be used with --flip, --make-[b]pgen/--make-bed, and\nno other non-merge commands.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(pc.zero_cluster_fname)) {
        logerrputs("Error: --zero-cluster must be used with --make-[b]pgen/--make-bed and no other\nnon-merge commands.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Exportf | kfCommand1Pmerge))) {
        if (unlikely((pc.misc_flags & kfMiscMajRef) || pc.ref_allele_flag || pc.alt_allele_flag || (pc.fa_flags & kfFaRefFrom))) {
          logerrputs("Error: Flags which alter REF/ALT allele settings (--maj-ref, --ref-allele,\n--alt[1]-allele, --ref-from-fa) must be used with\n--make-bed/--make-[b]pgen/--export and no other non-merge commands.\n");
          goto main_ret_INVALID_CMDLINE;
        }
        if (unlikely(pc.fa_flags & kfFaNormalize)) {
          logerrputs("Error: --normalize must be used with --make-bed/--make-[b]pgen/--export and no\nother non-merge commands.\n");
          goto main_ret_INVALID_CMDLINE;
        }
      }
    }
    if (pc.flip_info.subset_fname || (make_plink2_flags & kfMakePlink2FillMissingWithRef)) {
      if (unlikely(!(pc.command_flags1 & kfCommand1MakePlink2))) {
        logerrprintf("Error: --%s must be used with --make-[b]pgen/--make-bed.\n", pc.flip_info.subset_fname? "flip-subset" : "fill-missing-with-ref");
        goto main_ret_INVALID_CMDLINE_A;
      }
      if (unlikely(pc.flip_info.subset_fname && (make_plink2_flags & kfMakePlink2TrimAlts))) {
        logerrputs("Error: --flip-subset cannot be used with --make-[b]pgen/--make-bed's\n'trim-alts' modifier.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }
    if (pc.fa_fname && (((make_plink2_flags & kfMakePvar) && (pc.pvar_psam_flags & (kfPvarColXheader | kfPvarColVcfheader))) || (pc.exportf_info.flags & (kfExportfVcf | kfExportfBcf)))) {
      pc.fa_flags |= kfFaScrapeLengths;
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
      if (unlikely(pc.dosage_erase_thresh > (kDosageMid / 10))) {
        logerrputs("Error: --dosage-erase-threshold value cannot be larger than (default)\n--hard-call-threshold value.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    } else {
      if (unlikely(pc.dosage_erase_thresh > pc.hard_call_thresh)) {
        logerrputs("Error: --dosage-erase-threshold value cannot be larger than\n--hard-call-threshold value.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }
    {
      const uint32_t mode_flag = oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefLast | kfOxfordImportRefUnknown);
      if (unlikely(mode_flag & (mode_flag - 1))) {
        logerrputs("Error: --data/--[b]gen 'ref-first', 'ref-last', and 'ref-unknown' modifiers\ncannot be used together.\n");
        goto main_ret_INVALID_CMDLINE;
      }
    }
    if (unlikely(!strcmp(pc.missing_catname, pc.output_missing_pheno))) {
      logerrputs("Error: --missing-catname and --output-missing-phenotype strings can't match.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((pc.misc_flags & kfMiscChrOverrideCmdline) && (!chr_info.chrset_source))) {
      logerrputs("Error: --chr-override requires an explicit chromosome set.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((xload & kfXloadPlink1Dosage) && (!(load_params & kfLoadParamsPsam)))) {
      logerrputs("Error: --import-dosage requires a .fam file.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((pc.misc_flags & kfMiscCovarColNums) && (!pc.covar_fname) && (!pc.pheno_fname))) {
      logerrputs("Error: --covar-col-nums requires --covar or --pheno.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((pc.grm_flags & kfGrmMatrixShapemask) && (pc.misc_flags & kfMiscNoIdHeader)) {
      pc.grm_flags |= kfGrmNoIdHeader;
    }
    if (unlikely(pc.score_info.qsr_range_fname && (!pc.score_info.input_fname))) {
      logerrputs("Error: --q-score-range cannot be used without --score[-list].\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely(score_col_nums_present && (!(pc.command_flags1 & kfCommand1Score)))) {
      logerrputs("Error: --score-col-nums must be used with --score[-list].\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely(r2_required && (!(pc.command_flags1 & kfCommand1Vcor)))) {
      logerrputs("Error: --ld-window.../--ld-snp... must be used with --r[2]-[un]phased.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((pc.mendel_info.flags & kfMendelDuos) && (!((pc.filter_flags & kfFilterMendel) || (pc.command_flags1 & kfCommand1MendelReport) || (make_plink2_flags & kfMakePlink2SetMeMissing))))) {
      logerrputs("Error: --mendel-duos must be used with --me, --mendel, or --set-me-missing.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely(pc.rename_chrs_fname && (pc.sort_vars_mode <= kSortNone))) {
      logerrputs("Error: --rename-chrs must be used with --sort-vars.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (pc.extract_col_cond_info.params) {
      if (unlikely((!pc.extract_col_cond_info.match_substr) && pc.extract_col_cond_info.match_flattened && pc.extract_col_cond_info.mismatch_flattened)) {
        logerrputs("Error: --extract-col-cond-match and --extract-col-cond-mismatch can only be\nused together when --extract-col-cond-substr is specified.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      if (unlikely(pc.extract_col_cond_info.max < pc.extract_col_cond_info.min)) {
        logerrputs("Error: --extract-col-cond-max value can't be smaller than\n--extract-col-cond-min value.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }
    if (pc.check_sex_info.flags & kfCheckSexImpute) {
      if (unlikely(!(pc.command_flags1 & (kfCommand1MakePlink2 | kfCommand1Exportf | kfCommand1WriteCovar)))) {
        logerrputs("Error: --impute-sex must be used with\n--make-[b]pgen/--make-bed/--export/--write-covar.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      if (unlikely(pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Exportf | kfCommand1WriteCovar | kfCommand1CheckOrImputeSex)))) {
        // enforce this to reduce order-of-operations confusion
        logerrputs("Error: --impute-sex cannot be used with any commands other than\n--make-[b]pgen/--make-bed/--export/--write-covar.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }
    if (!permit_multiple_inclusion_filters) {
      // Permit only one position- or ID-based variant inclusion filter, since
      // it's not immediately obvious whether the union or intersection should
      // be taken with multiple inclusion filters.
      // However, multiple exclusion filters are fine.  (Also,
      // --autosome[-par]/--chr is exempted since it's more obvious how they
      // interact with other filters.)
      const uint32_t inclusion_filter_extract = (pc.extract_fnames != nullptr);
      const uint32_t inclusion_filter_extract_col_cond = (pc.extract_col_cond_info.params != nullptr);
      const uint32_t inclusion_filter_extract_intersect = (pc.extract_intersect_fnames != nullptr);
      const uint32_t inclusion_filter_fromto_id = pc.varid_from || pc.varid_to;
      const uint32_t inclusion_filter_fromto_bp = (pc.from_bp != -1) || (pc.to_bp != -1);
      const uint32_t inclusion_filter_snpflag = (pc.varid_snp != nullptr);
      const uint32_t inclusion_filter_snpsflag = !!pc.snps_range_list.name_ct;
      if (unlikely(inclusion_filter_extract + inclusion_filter_extract_col_cond + inclusion_filter_extract_intersect + inclusion_filter_fromto_id + inclusion_filter_fromto_bp + inclusion_filter_snpflag + inclusion_filter_snpsflag > 1)) {
        logerrputs("Error: Multiple variant inclusion filters specified (--extract,\n--extract-col-cond, --extract-intersect, --from/--to, --from-bp/--to-bp, --snp,\n--snps).  Add --force-intersect if you really want the intersection of these\nsets.  (If your variant IDs are unique, you can extract the union by e.g.\nrunning --write-snplist for each set, followed by --extract on all the .snplist\nfiles.)\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }
    if (unlikely(pmerge_required && (!(pc.command_flags1 & kfCommand1Pmerge)))) {
      if (delete_pmerge_result) {
        logerrputs("Error: --delete-pmerge-result must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_cm_mode != kMergeInfoCmModeNmFirst) {
        logerrputs("Error: --merge-cm-mode must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_filter_mode != kMergeFilterModeNmFirst) {
        logerrputs("Error: --merge-filter-mode must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_info_mode != kMergeInfoCmModeNmFirst) {
        logerrputs("Error: --merge-info-mode must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_info_sort != kSortNone) {
        logerrputs("Error: --merge-info-sort must be used with --pmerge[-list].\n");
      } else if (pmerge_info.max_allele_ct != 0) {
        logerrputs("Error: --merge-max-allele-ct must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_mode != kMergeModeNmMatch) {
        logerrputs("Error: --merge-mode must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_parents_mode != kMergePhenoModeNmMatch) {
        logerrputs("Error: --merge-parents-mode must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_pheno_mode != kMergePhenoModeNmMatch) {
        logerrputs("Error: --merge-pheno-mode must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_pheno_sort != kSortNone) {
        logerrputs("Error: --merge-pheno-sort must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_qual_mode != kMergeQualModeMin) {
        logerrputs("Error: --merge-qual-mode must be used with --pmerge[-list].\n");
      } else if (pmerge_info.flags & kfPmergeSids) {
        logerrputs("Error: --merge-sids must be used with --pmerge[-list].\n");
      } else if (pmerge_info.merge_xheader_mode != kMergeXheaderModeFirst) {
        logerrputs("Error: --merge-xheader-mode must be used with --pmerge[-list].\n");
      } else if (pmerge_info.flags & kfPmergeMultiallelicsAlreadyJoined) {
        logerrputs("Error: --multiallelics-already-joined must be used with --pmerge[-list].\n");
      } else if (pmerge_info.flags & kfPmergePhenoInnerJoin) {
        logerrputs("Error: --pheno-inner-join must be used with --pmerge[-list].\n");
      } else {
        assert(0);
      }
      goto main_ret_INVALID_CMDLINE;
    }
    if (delete_pmerge_result) {
      if (unlikely(pc.command_flags1 == kfCommand1Pmerge)) {
        logerrputs("Error: --delete-pmerge-result doesn't make sense when --pmerge[-list] is run\nwithout any other commands.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }

#ifdef USE_MKL
    if (!mkl_native) {
#  ifdef USE_SSE42
#    ifdef USE_AVX2
      int status = mkl_cbwr_set(MKL_CBWR_AVX2);
#    else
      int status = mkl_cbwr_set(MKL_CBWR_SSE4_2);
#    endif
      if (status != MKL_CBWR_SUCCESS) {
        // could happen for AMD processors.  This is the least-bad option?
        mkl_cbwr_set(MKL_CBWR_COMPATIBLE);
      }
#  else
      mkl_cbwr_set(MKL_CBWR_COMPATIBLE);
#  endif
    }
#endif
    sfmt_t main_sfmt;
    if (!rseeds) {
      uint32_t seed = S_CAST(uint32_t, time(nullptr));
      snprintf(g_logbuf, kLogbufSize, "Random number seed: %u\n", seed);
      logputs_silent(g_logbuf);
      sfmt_init_gen_rand(&main_sfmt, seed);
    } else {
      if (rseed_ct == 1) {
        sfmt_init_gen_rand(&main_sfmt, rseeds[0]);
      } else {
        sfmt_init_by_array(&main_sfmt, rseeds, rseed_ct);
      }
      free(rseeds);
      rseeds = nullptr;
    }

    if (unlikely(CmdlineParsePhase3(0, malloc_size_mib, memory_require, &pcm, &bigstack_ua))) {
      goto main_ret_NOMEM;
    }

    if (pc.max_thread_ct > 8) {
      logprintf("Using up to %u threads (change this with --threads).\n", pc.max_thread_ct);
    } else {
      // "x compute threads" instead of "x threads" since e.g. when
      // max_thread_ct == 2, some code will use one I/O thread and two
      // compute threads.  Not worth the trouble of writing special-case code
      // to avoid that: with 2 cores, the I/O thread usually isn't busy enough
      // to justify only 1 compute thread.
      logprintf("Using %s%u compute thread%s.\n", (pc.max_thread_ct > 1)? "up to " : "", pc.max_thread_ct, (pc.max_thread_ct == 1)? "" : "s");
    }
    if (randmem) {
      reterr = RandomizeBigstack(pc.max_thread_ct, &main_sfmt);
      if (unlikely(reterr)) {
        goto main_ret_1;
      }
    }

    print_end_time = 1;

    if (batch_job) {
      if (adjust_file_info.fname) {
        reterr = AdjustFile(&adjust_file_info, pc.ln_pfilter, pc.output_min_ln, pc.max_thread_ct, outname, outname_end);
        if (unlikely(reterr)) {
          goto main_ret_1;
        }
      }
      if (pc.gwas_ssf_info.fname || pc.gwas_ssf_info.list_fname) {
        reterr = GwasSsfStandalone(&pc.gwas_ssf_info, pc.max_thread_ct);
        if (unlikely(reterr)) {
          goto main_ret_1;
        }
      }
      if (skip_main) {
        goto main_ret_1;
      }
    }

    // nonstandard cases (CNV, etc.) here
    if ((pc.command_flags1 == kfCommand1PgenInfo) && (load_params != kfLoadParamsPfileAll) && (!xload)) {
      // special case: don't require .psam/.pvar file
      if (unlikely(!(load_params & kfLoadParamsPgen))) {
        logerrputs("Error: --pgen-info requires a .pgen file.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      reterr = PgenInfoStandalone(pgenname, pc.pginame);
    } else {
      // --real-ref-alleles is an exception since it applies to merge.
      if (unlikely(pc.dependency_flags && (!(pc.command_flags1 & (~kfCommand1Pmerge))))) {
        logerrputs("Error: Basic file conversions do not support regular filter or transform\noperations.  Rerun your command with --make-bed/--make-[b]pgen.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if ((xload && pc.command_flags1 && (import_flags & kfImportKeepAutoconv)) || ((pc.command_flags1 & kfCommand1Pmerge) & (pc.command_flags1 & (~kfCommand1Pmerge)))) {
        // pfile-affecting input/output settings must be consistent since we
        // need to reload what we write.
        if (unlikely(pc.misc_flags & kfMiscAffection01)) {
          logerrputs("Error: --1 cannot be used with --keep-autoconv or non-standalone\n--pmerge[-list].\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
      }
      if (xload) {
        char* convname_end = outname_end;
        if (pc.command_flags1) {
          if (!(import_flags & kfImportKeepAutoconv)) {
            convname_end = Stpcpy(convname_end, "-temporary");
          }
        } else {
          import_flags |= kfImportKeepAutoconv;
        }
        const uint32_t convname_slen = convname_end - outname;
        const LoadFilterLogFlags load_filter_log_import_flags = pc.load_filter_log_flags & kfLoadFilterLogImportMask;
        uint32_t pgen_generated = 1;
        uint32_t pgi_generated = 0;
        uint32_t psam_generated = 1;

        const uint32_t is_update_or_impute_sex = (pc.update_sex_info.fname || (pc.check_sex_info.flags & kfCheckSexImpute));
        // Compress by default for VCF/BCF due to potentially large INFO
        // section; otherwise only do it in "--keep-autoconv vzs" case.
        uint32_t pvar_is_compressed;
        if (xload & (kfXloadVcf | kfXloadBcf)) {
          const uint32_t no_samples_ok = !((pc.dependency_flags & (kfFilterAllReq | kfFilterPsamReq)) || (pc.command_flags1 & kfCommand1Pmerge));
          const uint32_t is_vcf = (xload / kfXloadVcf) & 1;
          pvar_is_compressed = ((import_flags & (kfImportKeepAutoconv | kfImportKeepAutoconvVzs)) != kfImportKeepAutoconv);
          if (no_samples_ok && is_vcf && (!(import_flags & (kfImportKeepAutoconv | kfImportVcfRefNMissing))) && pc.command_flags1) {
            // special case: just treat the VCF as a .pvar file
            strcpy(pvarname, pgenname);
            pgenname[0] = '\0';
            goto main_reinterpret_vcf_instead_of_converting;
          }
          // Default to compression level 1 for temporary .pvar files for now.
          //
          // Level 1 may actually be best in a much wider variety of scenarios
          // as of this writing, but I won't try to tune any other compression
          // defaults for now.  Interestingly, the current setup is actually a
          // bit backwards given observed zstd behavior: the long INFO fields
          // motivating automatic compression here are best handled with level
          // 3, while the other simpler text files generated by plink2 are
          // likely to compress *better*, not just faster, with level 1 for
          // some reason.
          const uint32_t zst_level = g_zst_level;
          if (!(import_flags & kfImportKeepAutoconv)) {
            g_zst_level = 1;
          }
          if (is_vcf) {
            reterr = VcfToPgen(pgenname, (load_params & kfLoadParamsPsam)? psamname : nullptr, const_fid, vcf_dosage_import_field, pc.misc_flags, import_flags, load_filter_log_import_flags, no_samples_ok, is_update_or_impute_sex, !!pc.splitpar_bound2, pc.sort_vars_mode > kSortNone, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, id_delim, idspace_to, vcf_min_gq, vcf_min_dp, vcf_max_dp, vcf_half_call, pc.fam_cols, import_max_allele_ct, pc.max_thread_ct, outname, convname_end, &chr_info, &pgen_generated, &psam_generated);
          } else {
            reterr = BcfToPgen(pgenname, (load_params & kfLoadParamsPsam)? psamname : nullptr, const_fid, vcf_dosage_import_field, pc.misc_flags, import_flags, load_filter_log_import_flags, no_samples_ok, is_update_or_impute_sex, !!pc.splitpar_bound2, pc.sort_vars_mode > kSortNone, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, id_delim, idspace_to, vcf_min_gq, vcf_min_dp, vcf_max_dp, vcf_half_call, pc.fam_cols, import_max_allele_ct, pc.max_thread_ct, outname, convname_end, &chr_info, &pgen_generated, &psam_generated);
          }
          g_zst_level = zst_level;
        } else {
          pvar_is_compressed = (import_flags / kfImportKeepAutoconvVzs) & 1;
          // We have previously verified that --1 and either --keep-autoconv
          // or non-standalone --pmerge[-list] haven't been specified
          // simultaneously.  Therefore:
          // - If --1 is specified, and we're executing either --pmerge[-list]
          //   or at least one regular plink2 command, the temporary .psam file
          //   should use 0/1 coding for binary phenotypes.
          // - However, if --1 is specified and we aren't executing any regular
          //   plink2 commands, the output .psam is not temporary and should
          //   use 1/2 coding.
          const uint32_t psam_01 = (pc.misc_flags & kfMiscAffection01) && pc.command_flags1;
          if (xload & kfXloadPed) {
            reterr = PedmapToPgen(pgenname, pvarname, pc.missing_catname, pc.misc_flags, import_flags, load_filter_log_import_flags, psam_01, pc.fam_cols, pc.missing_pheno, pc.input_missing_geno_char, pc.max_thread_ct, outname, convname_end, &chr_info);
          } else if (xload & kfXloadTped) {
            reterr = TpedToPgen(pgenname, psamname, pc.missing_catname, pc.misc_flags, import_flags, load_filter_log_import_flags, pc.fam_cols, pc.missing_pheno, pc.input_missing_geno_char, pc.max_thread_ct, outname, convname_end, &chr_info, &psam_generated);
          } else if (xload & kfXloadOxGen) {
            reterr = OxGenToPgen(pgenname, psamname, const_fid, import_single_chr_str, ox_missing_code, pc.missing_catname, pc.misc_flags, import_flags, load_filter_log_import_flags, oxford_import_flags, psam_01, !!pc.splitpar_bound2, pc.sort_vars_mode > kSortNone, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, id_delim, pc.max_thread_ct, outname, convname_end, &chr_info);
          } else if (xload & kfXloadOxBgen) {
            reterr = OxBgenToPgen(pgenname, psamname, const_fid, import_single_chr_str, ox_missing_code, pc.missing_catname, pc.misc_flags, import_flags, load_filter_log_import_flags, oxford_import_flags, psam_01, is_update_or_impute_sex, !!pc.splitpar_bound2, pc.sort_vars_mode > kSortNone, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, id_delim, idspace_to, import_max_allele_ct, pc.max_thread_ct, outname, convname_end, &chr_info);
          } else if (xload & kfXloadOxHaps) {
            reterr = OxHapslegendToPgen(pgenname, pvarname, psamname, const_fid, import_single_chr_str, ox_missing_code, pc.missing_catname, pc.misc_flags, import_flags, load_filter_log_import_flags, oxford_import_flags, psam_01, is_update_or_impute_sex, !!pc.splitpar_bound2, pc.sort_vars_mode > kSortNone, id_delim, pc.max_thread_ct, outname, convname_end, &chr_info, &pgi_generated);
          } else if (xload & kfXloadEigGeno) {
            reterr = EigfileToPgen(pgenname, psamname, pvarname, const_fid, pc.missing_catname, pc.misc_flags, import_flags, load_filter_log_import_flags, psam_01, id_delim, pc.max_thread_ct, outname, convname_end, &chr_info);
          } else if (xload & kfXloadPlink1Dosage) {
            reterr = Plink1DosageToPgen(pgenname, psamname, (xload & kfXloadMap)? pvarname : nullptr, import_single_chr_str, &plink1_dosage_info, pc.missing_catname, pc.misc_flags, import_flags, load_filter_log_import_flags, psam_01, pc.fam_cols, pc.missing_pheno, pc.hard_call_thresh, pc.dosage_erase_thresh, import_dosage_certainty, pc.max_thread_ct, outname, convname_end, &chr_info);
          } else if (likely(xload & kfXloadGenDummy)) {
            reterr = GenerateDummy(&gendummy_info, import_flags, load_filter_log_import_flags, psam_01, pc.hard_call_thresh, pc.dosage_erase_thresh, pc.max_thread_ct, &main_sfmt, outname, convname_end, &chr_info);
          } else {
            // We should have errored out with a better message before this
            // point.
            assert(0);
            logerrputs("Error: Invalid input flag combination.\n");
            goto main_ret_INVALID_CMDLINE;
          }
        }
        if (reterr || (!pc.command_flags1)) {
          goto main_ret_1;
        }

        if (load_filter_log_import_flags) {
          pc.load_filter_log_flags |= kfLoadFilterLogImportMergeAlreadyApplied;
        }

        pc.hard_call_thresh = UINT32_MAX;
        pc.dosage_erase_thresh = 0;

        if (pgen_generated) {
          snprintf(memcpya(pgenname, outname, convname_slen), kMaxOutfnameExtBlen - 10, ".pgen");
        }
        snprintf(memcpya(pvarname, outname, convname_slen), kMaxOutfnameExtBlen - 10, pvar_is_compressed? ".pvar.zst" : ".pvar");
        if (psam_generated) {
          snprintf(memcpya(psamname, outname, convname_slen), kMaxOutfnameExtBlen - 10, ".psam");
        }
        if (!(import_flags & kfImportKeepAutoconv)) {
          if (pgen_generated) {
            if (unlikely(PushLlStr(pgenname, &file_delete_list))) {
              goto main_ret_NOMEM;
            }
            if (pgi_generated) {
              char* pgenname_end = &(pgenname[convname_slen + strlen(".pgen")]);
              strcpy_k(pgenname_end, ".pgi");
              if (unlikely(PushLlStr(pgenname, &file_delete_list))) {
                goto main_ret_NOMEM;
              }
              pgenname_end[0] = '\0';
            }
          }
          if (unlikely(PushLlStr(pvarname, &file_delete_list))) {
            goto main_ret_NOMEM;
          }
          if (psam_generated) {
            if (unlikely(PushLlStr(psamname, &file_delete_list))) {
              goto main_ret_NOMEM;
            }
          }
        }
        *outname_end = '\0';
      }
    main_reinterpret_vcf_instead_of_converting:
      if ((pc.dependency_flags & kfFilterOpportunisticPgen) && (pgenname[0] != '\0')) {
        pc.dependency_flags |= kfFilterAllReq;
      }
      if (pmerge_info.list_fname) {
        if (unlikely((!xload) && (load_params != kfLoadParamsPfileAll) && (load_params!= kfLoadParams0))) {
          logerrputs("Error: --pmerge-list cannot be used with a proper subset of {--pgen, --pvar,\n--psam}.  You must specify either none or all of those filenames.\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
      } else if ((pc.dependency_flags & kfFilterAllReq) || (pc.command_flags1 & kfCommand1Pmerge)) {
        if (unlikely((!xload) && (load_params != kfLoadParamsPfileAll))) {
          logerrputs("Error: A full fileset (.pgen/.bed + .pvar/.bim + .psam/.fam) is required for\nthis.\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
      } else {
        // no genotype file required
        if (!(pc.dependency_flags & kfFilterOpportunisticPgen)) {
          pgenname[0] = '\0';
        }

        if (pc.dependency_flags & kfFilterPvarReq) {
          if (unlikely((!xload) && (!(load_params & kfLoadParamsPvar)))) {
            logerrputs("Error: A .pvar/.bim file is required for this.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else {
          pvarname[0] = '\0';
        }
        if (pc.dependency_flags & kfFilterPsamReq) {
          if (unlikely((!xload) && (!(load_params & kfLoadParamsPsam)))) {
            logerrputs("Error: A .psam/.fam file is required for this.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else {
          psamname[0] = '\0';
        }
      }

      if (pc.command_flags1 & kfCommand1Pmerge) {
        const LoadFilterLogFlags load_filter_log_merge_flags = pc.load_filter_log_flags & kfLoadFilterLogMergeMask;
        char* merge_outname_end = outname_end;
        if (make_plink2_flags & (kfMakePgen | kfMakePvar | kfMakePsam)) {
          merge_outname_end = strcpya_k(merge_outname_end, "-merge");
        }
        reterr = Pmerge(&pmerge_info, pc.sample_sort_fname, pc.missing_catname, pc.varid_template_str, pc.varid_multi_template_str, pc.varid_multi_nonsnp_template_str, pc.missing_varid_match, pc.misc_flags, load_filter_log_merge_flags, pc.sample_sort_mode, pc.fam_cols, pc.missing_pheno, pc.new_variant_id_max_allele_slen, pc.input_missing_geno_char, pc.max_thread_ct, pc.sort_vars_mode, pgenname, psamname, pvarname, outname, merge_outname_end, &chr_info);
        if (unlikely(reterr)) {
          goto main_ret_1;
        }
        pc.command_flags1 ^= kfCommand1Pmerge;
        if (!pc.command_flags1) {
          goto main_ret_1;
        }
        if (delete_pmerge_result) {
          if (unlikely(PushLlStr(pgenname, &file_delete_list) ||
                       PushLlStr(pvarname, &file_delete_list) ||
                       PushLlStr(psamname, &file_delete_list))) {
            goto main_ret_NOMEM;
          }
        }
        // Don't need to repeat --set-{all,missing}-var-ids.
        if (pc.varid_template_str) {
          free(pc.varid_template_str);
          pc.varid_template_str = nullptr;
        }
        if (load_filter_log_merge_flags) {
          pc.load_filter_log_flags |= kfLoadFilterLogImportMergeAlreadyApplied;
        }
        // --real-ref-alleles communicates to --pmerge[-list] that all input
        // .bed filesets have real REF alleles.
        // However, we then need to clear this flag since Plink2Core() would
        // otherwise error out (since --real-ref-alleles does not make sense on
        // the merged .pgen).
        pc.misc_flags &= ~kfMiscRealRefAlleles;
      }

      if ((pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Validate | kfCommand1WriteSnplist | kfCommand1WriteCovar | kfCommand1WriteSamples))) || ((pc.command_flags1 & kfCommand1MakePlink2) && (pc.sort_vars_mode <= kSortNone))) {
        // split-chromosome prohibited for all commands unless explicitly
        // permitted here
        pc.dependency_flags |= kfFilterNoSplitChr;
      }

      BLAS_SET_NUM_THREADS(1);
      reterr = Plink2Core(&pc, make_plink2_flags, pgenname, psamname, pvarname, outname, outname_end, king_cutoff_fprefix, &chr_info, &main_sfmt);
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
      logerrputs("Warning: No output requested.  (Did you forget --make-bed/--make-[b]pgen?)\nExiting.\n");
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
  if (reterr == kPglRetNomemCustomMsg) {
    if (g_failed_alloc_attempt_size) {
      logerrprintf("Failed allocation size: %" PRIu64 "\n", g_failed_alloc_attempt_size);
    }
    reterr = kPglRetNomem;
  } else {
    DispExitMsg(reterr);
  }
  while (0) {
  main_ret_NOMEM_NOLOG:
    PrintVer();
  main_ret_NOMEM_NOLOG2:
    fputs(kErrstrNomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIu64 "\n", g_failed_alloc_attempt_size);
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
  free_cond(pc.zero_cluster_phenoname);
  free_cond(pc.zero_cluster_fname);
  free_cond(pc.rename_chrs_fname);
  free_cond(pc.indv_str);
  free_cond(pc.not_covar_flattened);
  free_cond(pc.not_pheno_flattened);
  free_cond(pc.indep_preferred_fname);
  free_cond(pc.vscore_fname);
  free_cond(pc.recover_var_ids_fname);
  free_cond(pc.update_parental_ids_fname);
  free_cond(pc.update_sample_ids_fname);
  free_cond(pc.update_name_flag);
  free_cond(pc.update_map_flag);
  free_cond(pc.update_chr_flag);
  free_cond(pc.alt_allele_flag);
  free_cond(pc.ref_allele_flag);
  free_cond(pc.keep_col_match_name);
  free_cond(pc.keep_col_match_flattened);
  free_cond(pc.keep_col_match_fname);
  free_cond(pc.require_no_info_flattened);
  free_cond(pc.require_info_flattened);
  free_cond(pc.king_table_require_fnames);
  free_cond(pc.king_table_subset_fname);
  free_cond(pc.fa_fname);
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
  free_cond(pc.extract_intersect_fnames);
  free_cond(pc.extract_fnames);
  free_cond(pc.sample_sort_fname);
  free_cond(pc.covar_fname);
  free_cond(pc.pheno_fname);
  free_cond(pc.varid_exclude_snp);
  free_cond(pc.varid_snp);
  free_cond(pc.varid_to);
  free_cond(pc.varid_from);
  free_cond(pc.missing_varid_match);
  free_cond(pc.varid_multi_nonsnp_template_str);
  free_cond(pc.varid_multi_template_str);
  free_cond(pc.varid_template_str);
  free_cond(pc.var_filter_exceptions_flattened);
  free_cond(pc.pginame);
  if (file_delete_list) {
    do {
      LlStr* llstr_ptr = file_delete_list->next;
      if (unlikely(unlink(file_delete_list->str))) {
        logerrprintfww("Error: Failed to delete %s : %s.\n", file_delete_list->str, strerror(errno));
        if (reterr == kPglRetSuccess) {
          reterr = kPglRetWriteFail;
        }
      }
      free(file_delete_list);
      file_delete_list = llstr_ptr;
    } while (file_delete_list);
  }
  CleanupChrInfo(&chr_info);
  CleanupGenDummy(&gendummy_info);
  CleanupFlip(&pc.flip_info);
  CleanupVcor(&pc.vcor_info);
  CleanupClump(&pc.clump_info);
  CleanupGwasSsf(&pc.gwas_ssf_info);
  CleanupExportf(&pc.exportf_info);
  CleanupExtractColCond(&pc.extract_col_cond_info);
  CleanupCmpExpr(&pc.exclude_if_info_expr);
  CleanupCmpExpr(&pc.extract_if_info_expr);
  CleanupCmpExpr(&pc.remove_if_expr);
  CleanupCmpExpr(&pc.keep_if_expr);
  CleanupPgenDiff(&pc.pgen_diff_info);
  CleanupPmerge(&pmerge_info);
  CleanupFst(&pc.fst_info);
  CleanupScore(&pc.score_info);
  CleanupGlm(&pc.glm_info);
  CleanupSdiff(&pc.sdiff_info);
  CleanupLd(&pc.ld_info);
  CleanupUpdateSex(&pc.update_sex_info);
  CleanupUpdateAlleles(&pc.update_alleles_info);
  CleanupRangeList(&pc.vscore_col_idx_range_list);
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
  free_cond(bigstack_ua);
  return S_CAST(int32_t, reterr);
}
