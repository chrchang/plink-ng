// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Uncomment "#define NOLAPACK" in plink_common.h to build without LAPACK.

#include "plink_common.h"
#include <ctype.h>
#include <locale.h>
#include <time.h>
#include <unistd.h>

#ifdef __APPLE__
  #include <sys/types.h>
  #include <sys/sysctl.h>
#endif

#include "plink_assoc.h"
#include "plink_calc.h"
#include "plink_cnv.h"
#include "plink_data.h"
#include "plink_dosage.h"
#include "plink_family.h"
#include "plink_filter.h"
#include "plink_glm.h"
#include "plink_help.h"
#include "plink_homozyg.h"
#include "plink_lasso.h"
#include "plink_ld.h"
#include "plink_misc.h"
#include "plink_perm.h"
#ifdef __cplusplus
  #ifndef _WIN32
    #include "plink_rserve.h"
  #endif
#endif
#include "plink_set.h"
#include "plink_stats.h"
#include "pigz.h"

#define DEFAULT_PPC_GAP 500000

#define DEFAULT_IBS_TEST_PERMS 100000

// should switch this to IMPORT_... and make it include .ped/Oxford too; only
// .bed and future .pgen format don't require an import step
#define LOAD_RARE_GRM 1
#define LOAD_RARE_GRM_BIN 2
#define LOAD_RARE_LGEN 4
#define LOAD_RARE_TRANSPOSE 8
#define LOAD_RARE_TPED 0x10
#define LOAD_RARE_TFAM 0x20
#define LOAD_RARE_TRANSPOSE_MASK (LOAD_RARE_TRANSPOSE | LOAD_RARE_TPED | LOAD_RARE_TFAM)
#define LOAD_RARE_DUMMY 0x40
#define LOAD_RARE_SIMULATE 0x80
#define LOAD_RARE_CNV 0x100
#define LOAD_RARE_GVAR 0x200
#define LOAD_RARE_23 0x400
// er, this won't actually be rare...
#define LOAD_RARE_VCF 0x800
#define LOAD_RARE_BCF 0x1000

#define LOAD_RARE_DOSAGE 0x2000

#define LOAD_PARAMS_PED 1
#define LOAD_PARAMS_MAP 2
#define LOAD_PARAMS_TEXT_ALL 3
#define LOAD_PARAMS_BED 4
#define LOAD_PARAMS_BIM 8
#define LOAD_PARAMS_FAM 0x10
#define LOAD_PARAMS_BFILE_ALL 0x1c
#define LOAD_PARAMS_OXGEN 0x20
#define LOAD_PARAMS_OXBGEN 0x40
#define LOAD_PARAMS_OXSAMPLE 0x80
#define LOAD_PARAMS_OX_ALL 0xe0

// maximum number of usable cluster computers, this is arbitrary though it
// shouldn't be larger than 2^31 - 1
#define PARALLEL_MAX 32768

static const char ver_str[] =
#ifdef STABLE_BUILD
  "PLINK v1.90b6.18"
#else
  "PLINK v1.90p"
#endif
#ifdef NOLAPACK
  "NL"
#endif
#ifdef __LP64__
  " 64-bit"
#else
  " 32-bit"
#endif
  " (16 Jun 2020)";
static const char ver_str2[] =
  // include leading space if day < 10, so character length stays the same
  ""
#ifdef STABLE_BUILD
  "" // (don't want this when version number has two trailing digits)
#else
  "    " // (don't want this when version number has e.g. "b3" before "p")
#endif
#ifndef NOLAPACK
  "  "
#endif
  "        www.cog-genomics.org/plink/1.9/\n"
  "(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3"
#if SPECIES_DEFAULT > 0
  "\nRecompiled with default species = "
  #if SPECIES_DEFAULT == SPECIES_COW
    "cow"
  #elif SPECIES_DEFAULT == SPECIES_DOG
    "dog"
  #elif SPECIES_DEFAULT == SPECIES_HORSE
    "horse"
  #elif SPECIES_DEFAULT == SPECIES_MOUSE
    "mouse"
  #elif SPECIES_DEFAULT == SPECIES_RICE
    "rice"
  #elif SPECIES_DEFAULT == SPECIES_SHEEP
    "sheep"
  #endif
#endif
  "\n";
static const char errstr_append[] = "For more information, try \"" PROG_NAME_STR " --help <flag name>\" or \"" PROG_NAME_STR " --help | more\".\n";
#ifdef STABLE_BUILD
  #ifndef NOLAPACK
static const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --flip-scan, --merge-list,\n--write-snplist, --list-duplicate-vars, --freqx, --missing, --test-mishap,\n--hardy, --mendel, --ibc, --impute-sex, --indep-pairphase, --r2, --show-tags,\n--blocks, --distance, --genome, --homozyg, --make-rel, --make-grm-gz,\n--rel-cutoff, --cluster, --pca, --neighbour, --ibs-test, --regress-distance,\n--model, --bd, --gxe, --logistic, --dosage, --lasso, --test-missing,\n--make-perm-pheno, --tdt, --qfam, --annotate, --clump, --gene-report,\n--meta-analysis, --epistasis, --fast-epistasis, and --score.\n\n\"" PROG_NAME_STR " --help | more\" describes all functions (warning: long).\n";
  #else
static const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --flip-scan, --merge-list,\n--write-snplist, --list-duplicate-vars, --freqx, --missing, --test-mishap,\n--hardy, --mendel, --ibc, --impute-sex, --indep-pairphase, --r2, --show-tags,\n--blocks, --distance, --genome, --homozyg, --make-rel, --make-grm-gz,\n--rel-cutoff, --cluster, --neighbour, --ibs-test, --regress-distance, --model,\n--bd, --gxe, --logistic, --dosage, --lasso, --test-missing, --make-perm-pheno,\n--tdt, --qfam, --annotate, --clump, --gene-report, --meta-analysis,\n--epistasis, --fast-epistasis, and --score.\n\n\"" PROG_NAME_STR " --help | more\" describes all functions (warning: long).\n";
  #endif
#else
  #ifndef NOLAPACK
static const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --flip-scan, --merge-list,\n--write-snplist, --list-duplicate-vars, --freqx, --missing, --test-mishap,\n--hardy, --mendel, --ibc, --impute-sex, --indep-pairphase, --r2, --show-tags,\n--blocks, --distance, --genome, --homozyg, --make-rel, --make-grm-gz,\n--rel-cutoff, --cluster, --pca, --neighbour, --ibs-test, --regress-distance,\n--model, --bd, --gxe, --logistic, --dosage, --lasso, --test-missing,\n--make-perm-pheno, --unrelated-heritability, --tdt, --dfam, --qfam, --tucc,\n--annotate, --clump, --gene-report, --meta-analysis, --epistasis,\n--fast-epistasis, and --score.\n\n\"" PROG_NAME_STR " --help | more\" describes all functions (warning: long).\n";
  #else
static const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --flip-scan, --merge-list,\n--write-snplist, --list-duplicate-vars, --freqx, --missing, --test-mishap,\n--hardy, --mendel, --ibc, --impute-sex, --indep-pairphase, --r2, --show-tags,\n--blocks, --distance, --genome, --homozyg, --make-rel, --make-grm-gz,\n--rel-cutoff, --cluster, --neighbour, --ibs-test, --regress-distance, --model,\n--bd, --gxe, --logistic, --dosage, --lasso, --test-missing, --make-perm-pheno,\n--tdt, --dfam, --qfam, --tucc, --annotate, --clump, --gene-report,\n--meta-analysis, --epistasis, --fast-epistasis, and --score.\n\n\"" PROG_NAME_STR " --help | more\" describes all functions (warning: long).\n";
  #endif
#endif

static const char errstr_nomem[] = "Error: Out of memory.  The --memory flag may be helpful.\n";
static const char errstr_write[] = "Error: File write failure.\n";
static const char errstr_read[] = "Error: File read failure.\n";

void disp_exit_msg(int32_t retval) {
  switch (retval) {
  case RET_NOMEM:
    logprint("\n");
    logerrprint(errstr_nomem);
    if (g_failed_alloc_attempt_size) {
      LOGERRPRINTF("Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
    }
    break;
  case RET_WRITE_FAIL:
    logprint("\n");
    logerrprint(errstr_write);
    break;
  case RET_READ_FAIL:
    logprint("\n");
    logerrprint(errstr_read);
    break;
  case RET_THREAD_CREATE_FAIL:
    logprint("\n");
    logerrprint("Error: Failed to create thread.\n");
    break;
  }
}

// back to our regular program

const unsigned char acgt_reverse_arr[] = "1B2DEF3HIJKLMNOPQRS4";
const unsigned char acgt_arr[] = "ACGT";
// g_one_char_strs offsets (double)
const unsigned char acgt_reverse_arr1[] = "b\204d\210\212\214f\220\222\224\226\230\232\234\236\240\242\244\246h";
const unsigned char acgt_arr1[] = "\202\206\216\250";
// const unsigned char acgt_reverse_arr1[] = "\"D$HJL&PRTVXZ\\^`bdf(";
// const unsigned char acgt_arr1[] = "BFNh";

static inline unsigned char conditional_convert(unsigned char diff, unsigned char ucc2_max, const unsigned char* convert_arr, unsigned char ucc) {
  unsigned char ucc2 = ucc - diff;
  return (ucc2 < ucc2_max)? convert_arr[ucc2] : ucc;
}

static inline void conditional_convert1(unsigned char diff, unsigned char ucc2_max, const unsigned char* convert_arr, char** allele_ptr) {
  unsigned char ucc2 = ((unsigned char)(**allele_ptr)) - diff;
  if (ucc2 < ucc2_max) {
    *allele_ptr = (char*)(&(g_one_char_strs[convert_arr[ucc2]]));
  }
}

void allelexxxx_recode(uint32_t allelexxxx, char** marker_allele_ptrs, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct) {
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t recode_multichar = allelexxxx & ALLELE_RECODE_MULTICHAR;
  const unsigned char* convert_arr;
  const unsigned char* convert_arr1;
  char** map_ptr;
  char** map_ptr_stop;
  char* cptr;
  uint32_t marker_uidx_stop;
  unsigned char diff;
  unsigned char ucc2_max;
  unsigned char ucc;
  if (allelexxxx & ALLELE_RECODE_ACGT) {
    diff = 49;
    ucc2_max = 4;
    convert_arr = acgt_arr;
    convert_arr1 = acgt_arr1;
  } else {
    diff = 65;
    ucc2_max = 20;
    convert_arr = acgt_reverse_arr;
    convert_arr1 = acgt_reverse_arr1;
  }
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    map_ptr = &(marker_allele_ptrs[marker_uidx * 2]);
    map_ptr_stop = &(marker_allele_ptrs[marker_uidx_stop * 2]);
    marker_uidx = marker_uidx_stop;
    if (recode_multichar) {
      do {
	cptr = *map_ptr;
	if (!cptr[1]) {
	  conditional_convert1(diff, ucc2_max, convert_arr1, map_ptr);
	} else {
	  ucc = *cptr;
	  do {
	    *cptr = conditional_convert(diff, ucc2_max, convert_arr, ucc);
	    ucc = *(++cptr);
	  } while (ucc);
	}
      } while (++map_ptr < map_ptr_stop);
    } else {
      do {
        if (!(map_ptr[0][1])) {
          conditional_convert1(diff, ucc2_max, convert_arr1, map_ptr);
	}
      } while (++map_ptr < map_ptr_stop);
    }
  }
}

void calc_marker_reverse_bin(uintptr_t* marker_reverse, uintptr_t* marker_exclude, uint32_t unfiltered_marker_ct, uint32_t marker_ct, double* set_allele_freqs) {
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  double dxx;
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    for (; marker_uidx < marker_uidx_stop; marker_uidx++) {
      dxx = set_allele_freqs[marker_uidx];
      if (dxx < 0.5) {
	SET_BIT(marker_uidx, marker_reverse);
	set_allele_freqs[marker_uidx] = 1.0 - dxx;
      }
    }
  }
}

void swap_reversed_marker_alleles(uintptr_t unfiltered_marker_ct, uintptr_t* marker_reverse, char** marker_allele_ptrs) {
  uintptr_t marker_uidx = 0;
  char* swap_ptr;
  while (1) {
    next_set_ul_ck(marker_reverse, unfiltered_marker_ct, &marker_uidx);
    if (marker_uidx == unfiltered_marker_ct) {
      return;
    }
    swap_ptr = marker_allele_ptrs[marker_uidx * 2];
    marker_allele_ptrs[marker_uidx * 2] = marker_allele_ptrs[marker_uidx * 2 + 1];
    marker_allele_ptrs[marker_uidx * 2 + 1] = swap_ptr;
    marker_uidx++;
  }
}

static inline int32_t bed_suffix_conflict(uint64_t calculation_type, uint32_t recode_modifier) {
  return (calculation_type & CALC_MAKE_BED) || ((calculation_type & CALC_RECODE) && (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF | RECODE_RLIST)));
}

static inline uint32_t are_marker_pos_needed(uint64_t calculation_type, uint64_t misc_flags, char* cm_map_fname, char* set_fname, uint32_t min_bp_space, uint32_t genome_skip_write, uint32_t ld_modifier, uint32_t epi_modifier, uint32_t cluster_modifier) {
  return (calculation_type & (CALC_MAKE_BED | CALC_MAKE_BIM | CALC_RECODE | CALC_GENOME | CALC_HOMOZYG | CALC_LD_PRUNE | CALC_MODEL | CALC_GLM | CALC_CLUMP | CALC_BLOCKS | CALC_FLIPSCAN | CALC_TDT | CALC_QFAM | CALC_FST | CALC_SHOW_TAGS | CALC_DUPVAR | CALC_RPLUGIN | CALC_TUCC)) || (misc_flags & (MISC_EXTRACT_RANGE | MISC_EXCLUDE_RANGE)) || cm_map_fname || set_fname || min_bp_space || genome_skip_write || ((calculation_type & CALC_LD) && (!(ld_modifier & LD_MATRIX_SHAPEMASK))) || ((calculation_type & CALC_EPI) && (epi_modifier & EPI_FAST_CASE_ONLY)) || ((calculation_type & CALC_CMH) && (!(cluster_modifier & CLUSTER_CMH2)));
}

static inline uint32_t are_marker_cms_needed(uint64_t calculation_type, char* cm_map_fname, Two_col_params* update_cm, Ld_info* ldip) {
  if ((calculation_type & CALC_LD) && (!(ldip->modifier & LD_MATRIX_SHAPEMASK)) && (ldip->window_cm != -1)) {
    return MARKER_CMS_FORCED;
  }
  if (calculation_type & (CALC_MAKE_BED | CALC_MAKE_BIM | CALC_RECODE | CALC_TUCC)) {
    if (cm_map_fname || update_cm) {
      return MARKER_CMS_FORCED;
    } else {
      return MARKER_CMS_OPTIONAL;
    }
  }
  return 0;
}

static inline uint32_t are_marker_alleles_needed(uint64_t calculation_type, char* freqname, Homozyg_info* homozyg_ptr, Two_col_params* a1alleles, Two_col_params* a2alleles, uint32_t ld_modifier, uint32_t snps_only, uint32_t clump_modifier, uint32_t cluster_modifier, uint32_t rel_modifier) {
  return (freqname || (calculation_type & (CALC_FREQ | CALC_HARDY | CALC_MAKE_BED | CALC_MAKE_BIM | CALC_RECODE | CALC_MODEL | CALC_GLM | CALC_LASSO | CALC_LIST_23_INDELS | CALC_EPI | CALC_TESTMISHAP | CALC_SCORE | CALC_MENDEL | CALC_TDT | CALC_FLIPSCAN | CALC_QFAM | CALC_HOMOG | CALC_DUPVAR | CALC_RPLUGIN | CALC_DFAM | CALC_TUCC)) || ((calculation_type & CALC_HOMOZYG) && (homozyg_ptr->modifier & HOMOZYG_GROUP_VERBOSE)) || ((calculation_type & CALC_LD) && (ld_modifier & LD_INPHASE)) || ((calculation_type & CALC_PCA) && (rel_modifier & REL_PCA_VAR_WTS)) || ((calculation_type & CALC_CMH) && (!(cluster_modifier & CLUSTER_CMH2))) || a1alleles || a2alleles || snps_only || (clump_modifier & (CLUMP_VERBOSE | CLUMP_BEST)));
}

static inline int32_t relationship_or_ibc_req(uint64_t calculation_type) {
  return (relationship_req(calculation_type) || (calculation_type & CALC_IBC));
}

int32_t plink(char* outname, char* outname_end, char* bedname, char* bimname, char* famname, char* cm_map_fname, char* cm_map_chrname, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* keepfamname, char* removefamname, char* filtername, char* freqname, char* distance_wts_fname, char* read_dists_fname, char* read_dists_id_fname, char* evecname, char* mergename1, char* mergename2, char* mergename3, char* missing_mid_template, char* missing_marker_id_match, char* makepheno_str, char* phenoname_str, Two_col_params* a1alleles, Two_col_params* a2alleles, char* recode_allele_name, char* covar_fname, char* update_alleles_fname, char* read_genome_fname, Two_col_params* qual_filter, Two_col_params* update_chr, Two_col_params* update_cm, Two_col_params* update_map, Two_col_params* update_name, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* loop_assoc_fname, char* flip_fname, char* flip_subset_fname, char* sample_sort_fname, char* filtervals_flattened, char* condition_mname, char* condition_fname, char* filter_attrib_fname, char* filter_attrib_liststr, char* filter_attrib_sample_fname, char* filter_attrib_sample_liststr, char* rplugin_fname, char* rplugin_host_or_socket, int32_t rplugin_port, double qual_min_thresh, double qual_max_thresh, double thin_keep_prob, double thin_keep_sample_prob, uint32_t new_id_max_allele_len, uint32_t thin_keep_ct, uint32_t thin_keep_sample_ct, uint32_t min_bp_space, uint32_t mfilter_col, uint32_t fam_cols, int32_t missing_pheno, char* output_missing_pheno, uint32_t mpheno_col, uint32_t pheno_modifier, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, Family_info* fam_ip, double check_sex_fthresh, double check_sex_mthresh, uint32_t check_sex_f_yobs, uint32_t check_sex_m_yobs, double distance_exp, double min_maf, double max_maf, double geno_thresh, double mind_thresh, double hwe_thresh, double tail_bottom, double tail_top, uint64_t misc_flags, uint64_t filter_flags, uint64_t calculation_type, uint32_t dist_calc_type, uintptr_t groupdist_iters, uint32_t groupdist_d, uintptr_t regress_iters, uint32_t regress_d, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t splitx_bound1, uint32_t splitx_bound2, uint32_t ppc_gap, uint32_t sex_missing_pheno, uint32_t update_sex_col, uint32_t hwe_modifier, uint32_t min_ac, uint32_t max_ac, uint32_t genome_modifier, double genome_min_pi_hat, double genome_max_pi_hat, Homozyg_info* homozyg_ptr, Cluster_info* cluster_ptr, uint32_t neighbor_n1, uint32_t neighbor_n2, Set_info* sip, Ld_info* ldip, Epi_info* epi_ip, Clump_info* clump_ip, Rel_info* relip, Score_info* sc_ip, uint32_t recode_modifier, uint32_t allelexxxx, uint32_t merge_type, uint32_t sample_sort, int32_t marker_pos_start, int32_t marker_pos_end, int32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, Range_list* snps_range_list_ptr, uint32_t write_var_range_ct, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t write_covar_modifier, uint32_t write_covar_dummy_max_categories, uint32_t dupvar_modifier, uint32_t mwithin_col, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, uint32_t glm_mperm_val, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, double ci_size, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uint32_t gxe_mcovar, Aperm_info* apip, uint32_t mperm_save, uintptr_t ibs_test_perms, uint32_t perm_batch_size, double lasso_h2, uint32_t lasso_lambda_iters, double lasso_minlambda, Range_list* lasso_select_covars_range_list_ptr, uint32_t testmiss_modifier, uint32_t testmiss_mperm_val, uint32_t permphe_ct, int32_t known_procs, Ll_str** file_delete_list_ptr) {
  FILE* bedfile = nullptr;
  FILE* phenofile = nullptr;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t* marker_exclude = nullptr;
  uintptr_t marker_exclude_ct = 0;
  uintptr_t marker_ct = 0;
  uintptr_t max_marker_id_blen = 0;
  // set_allele_freqs = .bed set bit frequency in middle of loading process, A2
  //   allele frequency later.
  double* set_allele_freqs = nullptr;
  uintptr_t unfiltered_sample_ct = 0;
  uintptr_t unfiltered_sample_ct4 = 0;
  uintptr_t unfiltered_sample_ctl = 0;
  uintptr_t* sample_exclude = nullptr;
  uintptr_t sample_exclude_ct = 0;
  uintptr_t sample_ct = 0;
  uint32_t* sample_sort_map = nullptr;
  uintptr_t* founder_info = nullptr;
  uintptr_t* sex_nm = nullptr;
  uintptr_t* sex_male = nullptr;
  uint32_t genome_skip_write = (cluster_ptr->ppc != 0.0) && (!(calculation_type & CALC_GENOME)) && (!read_genome_fname);
  uint32_t marker_pos_needed = are_marker_pos_needed(calculation_type, misc_flags, cm_map_fname, sip->fname, min_bp_space, genome_skip_write, ldip->modifier, epi_ip->modifier, cluster_ptr->modifier);
  uint32_t marker_cms_needed = are_marker_cms_needed(calculation_type, cm_map_fname, update_cm, ldip);
  uint32_t marker_alleles_needed = are_marker_alleles_needed(calculation_type, freqname, homozyg_ptr, a1alleles, a2alleles, ldip->modifier, (filter_flags / FILTER_SNPS_ONLY) & 1, clump_ip->modifier, cluster_ptr->modifier, relip->modifier);

  // add other nchrobs subscribers later
  uint32_t nchrobs_needed = (calculation_type & CALC_GENOME) && freqname;
  uint32_t uii = 0;
  int64_t llxx = 0;
  uint32_t nonfounders = (misc_flags / MISC_NONFOUNDERS) & 1;
  uint32_t allow_no_samples = (misc_flags / MISC_ALLOW_NO_SAMPLES) & 1;
  uint32_t allow_no_variants = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  uint32_t pheno_all = pheno_modifier & PHENO_ALL;
  char* marker_ids = nullptr;
  uint32_t* marker_id_htable = nullptr;
  uint32_t marker_id_htable_size = 0;
  double* marker_cms = nullptr;
  // marker_allele_ptrs[2 * i] is id of A1 (usually minor) allele at marker i
  // marker_allele_ptrs[2 * i + 1] is id of A2 allele
  // Single-character allele names point to g_one_char_strs[]; otherwise
  // string allocation occurs on the heap.
  char** marker_allele_ptrs = nullptr;
  uintptr_t max_marker_allele_blen = 2; // er, this probably should be slen
  uintptr_t* marker_reverse = nullptr;
  int32_t retval = 0;
  uint32_t map_is_unsorted = 0;
  uint32_t affection = 0;
  uint32_t gender_unk_ct = 0;
  uintptr_t* pheno_nm = nullptr;
  uintptr_t* pheno_nm_datagen = nullptr; // --make-bed/--recode/--write-covar
  uintptr_t* orig_pheno_nm = nullptr; // --all-pheno + --pheno-merge
  uintptr_t* pheno_c = nullptr;
  uintptr_t* orig_pheno_c = nullptr;
  uintptr_t* geno_excl_bitfield = nullptr;
  uintptr_t* ac_excl_bitfield = nullptr;
  double* pheno_d = nullptr;
  double* orig_pheno_d = nullptr;
  char* sample_ids = nullptr;
  uintptr_t max_sample_id_len = 4;
  char* paternal_ids = nullptr;
  uintptr_t max_paternal_id_len = 2;
  char* maternal_ids = nullptr;
  uintptr_t max_maternal_id_len = 2;
  unsigned char* bigstack_mark = nullptr;
  uintptr_t cluster_ct = 0;
  uint32_t* cluster_map = nullptr; // unfiltered sample IDs
  // index for cluster_map, length (cluster_ct + 1)
  // cluster_starts[n+1] - cluster_starts[n] = length of cluster n (0-based)
  uint32_t* cluster_starts = nullptr;
  char* cluster_ids = nullptr;
  uintptr_t max_cluster_id_len = 2;
  double* mds_plot_dmatrix_copy = nullptr;
  uintptr_t* cluster_merge_prevented = nullptr;
  double* cluster_sorted_ibs = nullptr;
  char* cptr = nullptr;
  double missing_phenod = (double)missing_pheno;
  double ci_zt = 0.0;
  uintptr_t bed_offset = 3;
  uint32_t* marker_pos = nullptr;
  uint32_t hh_exists = 0;
  uint32_t pheno_ctrl_ct = 0;
  uintptr_t covar_ct = 0;
  char* covar_names = nullptr;
  uintptr_t max_covar_name_len = 0;
  uintptr_t* covar_nm = nullptr;
  double* covar_d = nullptr;
  uintptr_t* gxe_covar_nm = nullptr;
  uintptr_t* gxe_covar_c = nullptr;
  uintptr_t* pca_sample_exclude = nullptr;
  uintptr_t pca_sample_ct = 0;
  uintptr_t ulii = 0;
  uint32_t pheno_nm_ct = 0;
  uint32_t plink_maxsnp = 0;
  uint32_t plink_maxfid = 0;
  uint32_t plink_maxiid = 0;
  uint32_t max_bim_linelen = 0;
  uint32_t sample_male_ct = 0;
  uint32_t sample_f_ct = 0;
  uint32_t sample_f_male_ct = 0;
  unsigned char* bigstack_mark2 = nullptr;
  unsigned char* bigstack_mark_precluster = nullptr;
  unsigned char* bigstack_mark_postcluster = nullptr;
  uint32_t* nchrobs = nullptr;
  int32_t* hwe_lls = nullptr;
  int32_t* hwe_lhs = nullptr;
  int32_t* hwe_hhs = nullptr;
  int32_t* hwe_ll_cases = nullptr;
  int32_t* hwe_lh_cases = nullptr;
  int32_t* hwe_hh_cases = nullptr;
  int32_t* hwe_ll_allfs = nullptr;
  int32_t* hwe_lh_allfs = nullptr;
  int32_t* hwe_hh_allfs = nullptr;
  int32_t* hwe_hapl_allfs = nullptr;
  int32_t* hwe_haph_allfs = nullptr;
  pthread_t threads[MAX_THREADS];
  uint32_t* uiptr;
  double* rel_ibc;
  uintptr_t uljj;
  uint32_t ujj;
  uint32_t ukk;
  char* outname_end2;
  int32_t ii;
  int64_t llyy;
  int64_t llzz;
  Pedigree_rel_info pri;
  uintptr_t marker_uidx;

  if ((cm_map_fname || update_cm) && (!marker_cms_needed)) {
    LOGERRPRINTF("Error: --%s results would never be used.  (Did you forget --make-bed?)\n", cm_map_fname? "cm-map" : "update-cm");
    goto plink_ret_INVALID_CMDLINE;
  } else if (update_map && (!marker_pos_needed)) {
    logerrprint("Error: --update-map results would never be used.  (Did you forget --make-bed?)\n");
    goto plink_ret_INVALID_CMDLINE;
  }
  if (ci_size != 0.0) {
    ci_zt = ltqnorm(1 - (1 - ci_size) / 2);
  }
  if (relip->modifier & REL_CALC_COV) {
    relip->ibc_type = -1;
  }

  // famname[0] is nonzero iff we're not in the --merge-list special case
  if ((calculation_type & CALC_MAKE_BED) && famname[0]) {
#ifdef _WIN32
    uii = GetFullPathName(bedname, FNAMESIZE, g_textbuf, nullptr);
    if ((!uii) || (uii > FNAMESIZE))
#else
    if (!realpath(bedname, g_textbuf))
#endif
    {
      uii = strlen(bedname);
      if ((uii > 8) && ((!memcmp(&(bedname[uii - 8]), ".bed.bed", 8)) || (!memcmp(&(bedname[uii - 8]), ".bim.bed", 8)) || (!memcmp(&(bedname[uii - 8]), ".fam.bed", 8)))) {
	LOGERRPRINTFWW("Error: Failed to open %s. (--bfile expects a filename *prefix*; '.bed', '.bim', and '.fam' are automatically appended.)\n", bedname);
      } else {
        LOGERRPRINTFWW(g_errstr_fopen, bedname);
      }
      goto plink_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".bed", 5);
    if (realpath_identical(outname, g_textbuf, &(g_textbuf[FNAMESIZE + 64]))) {
      logerrprint("Warning: --make-bed input and output filenames match.  Appending '~' to input\nfilenames.\n");
      uii = strlen(bedname);
      memcpy(g_textbuf, bedname, uii + 1);
      memcpy(&(bedname[uii]), "~", 2);
      if (rename(g_textbuf, bedname)) {
	logerrprint("Error: Failed to append '~' to input .bed filename.\n");
	goto plink_ret_OPEN_FAIL;
      }
      uii = strlen(bimname);
      memcpy(g_textbuf, bimname, uii + 1);
      memcpy(&(bimname[uii]), "~", 2);
      if (rename(g_textbuf, bimname)) {
	logerrprint("Error: Failed to append '~' to input .bim filename.\n");
	goto plink_ret_OPEN_FAIL;
      }
      uii = strlen(famname);
      memcpy(g_textbuf, famname, uii + 1);
      memcpy(&(famname[uii]), "~", 2);
      if (rename(g_textbuf, famname)) {
	logerrprint("Error: Failed to append '~' to input .fam filename.\n");
	goto plink_ret_OPEN_FAIL;
      }
    }
  }

  if (calculation_type & CALC_MERGE) {
    if (((fam_cols & FAM_COL_13456) != FAM_COL_13456) || (misc_flags & MISC_AFFECTION_01) || (missing_pheno != -9)) {
      logerrprint("Error: --merge/--bmerge/--merge-list cannot be used with an irregularly\nformatted reference fileset (--no-fid, --no-parents, --no-sex, --no-pheno,\n--1, --missing-pheno).  Use --make-bed first.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    // Only append -merge to the filename stem if --make-bed or --recode lgen
    // is specified.
    // bugfix: also no need if --merge-mode 6 or 7.
    ulii = ((merge_type & MERGE_MODE_MASK) < 6) && bed_suffix_conflict(calculation_type, recode_modifier);
    if (ulii) {
      memcpy(outname_end, "-merge", 6);
    }
    retval = merge_datasets(bedname, bimname, famname, outname, ulii? &(outname_end[6]) : outname_end, mergename1, mergename2, mergename3, sample_sort_fname, calculation_type, merge_type, sample_sort, misc_flags, chrom_info_ptr);
    if (retval || (!(calculation_type & (~CALC_MERGE)))) {
      goto plink_ret_1;
    }
    if ((merge_type & MERGE_MODE_MASK) < 6) {
      uljj = (uintptr_t)(outname_end - outname) + (ulii? 6 : 0);
      memcpy(memcpya(bedname, outname, uljj), ".bed", 5);
      memcpy(memcpya(famname, bedname, uljj), ".fam", 5);
      memcpy(memcpya(bimname, bedname, uljj), ".bim", 5);
      if ((calculation_type & CALC_MAKE_BED) && ulii) {
	if (push_ll_str(bedname, file_delete_list_ptr) || push_ll_str(famname, file_delete_list_ptr) || push_ll_str(bimname, file_delete_list_ptr)) {
	  goto plink_ret_NOMEM;
	}
      }
    }
    // do we want to permit other operations with --merge-mode 6/7 at all?
    // seems like a recipe for confusion...
  }

  // don't use fopen_checked() here, since we want to customize the error
  // message.
  if (bedname[0]) {
    bedfile = fopen(bedname, FOPEN_RB);
    if (!bedfile) {
      uii = strlen(bedname);
      if ((uii > 8) && ((!memcmp(&(bedname[uii - 8]), ".bed.bed", 8)) || (!memcmp(&(bedname[uii - 8]), ".bim.bed", 8)) || (!memcmp(&(bedname[uii - 8]), ".fam.bed", 8)))) {
	LOGERRPRINTFWW("Error: Failed to open %s. (--bfile expects a filename *prefix*; '.bed', '.bim', and '.fam' are automatically appended.)\n", bedname);
      } else {
	LOGERRPRINTFWW(g_errstr_fopen, bedname);
      }
      goto plink_ret_OPEN_FAIL;
    }
  }
  // load .bim, count markers, filter chromosomes
  if (bimname[0]) {
    if (update_name) {
      ulii = 0;
      retval = scan_max_strlen(update_name->fname, update_name->colid, update_name->colx, update_name->skip, update_name->skipchar, &max_marker_id_blen, &ulii);
      if (retval) {
	goto plink_ret_1;
      }
      if (ulii > 80) {
	// only warn on long new marker ID, since if there's a long old marker
	// ID and no long new one, it's reasonable to infer that the user is
	// fixing the problem, so we shouldn't spam them.
	logerrprint("Warning: Unusually long new variant ID(s) in --update-name file.  Double-check\nyour file and command-line parameters, and consider changing your naming\nscheme if you encounter memory problems.\n");
      }
      if (ulii > max_marker_id_blen) {
	max_marker_id_blen = ulii;
      }
    }
    if (!marker_alleles_needed) {
      allelexxxx = 0;
    }
    retval = load_bim(bimname, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_blen, &marker_exclude, &set_allele_freqs, nchrobs_needed? (&nchrobs) : nullptr, &marker_allele_ptrs, &max_marker_allele_blen, &marker_ids, missing_mid_template, new_id_max_allele_len, missing_marker_id_match, chrom_info_ptr, &marker_cms, &marker_pos, misc_flags, filter_flags, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, &map_is_unsorted, marker_pos_needed, marker_cms_needed, marker_alleles_needed, ((!(calculation_type & (~(CALC_MAKE_BED | CALC_MAKE_BIM | CALC_MAKE_FAM)))) && (mind_thresh == 1.0) && (geno_thresh == 1.0) && (hwe_thresh == 0.0) && (!update_map) && (!freqname))? nullptr : "make-bed", ".bim file", &max_bim_linelen);
    if (retval) {
      goto plink_ret_1;
    }
    if (!update_chr) {
      if (finalize_chrom_info(chrom_info_ptr)) {
	goto plink_ret_NOMEM;
      }
    }
  }

  // load .fam, count samples
  if (famname[0]) {
    uii = fam_cols & FAM_COL_6;
    if (uii && phenoname) {
      uii = (pheno_modifier & PHENO_MERGE) && (!makepheno_str);
    }
    if (!uii) {
      pheno_modifier &= ~PHENO_MERGE; // nothing to merge
    }
    if (update_ids_fname) {
      ulii = 0;
      retval = scan_max_fam_indiv_strlen(update_ids_fname, 3, &max_sample_id_len);
      if (retval) {
	goto plink_ret_1;
      }
    } else if (update_parents_fname) {
      retval = scan_max_strlen(update_parents_fname, 3, 4, 0, '\0', &max_paternal_id_len, &max_maternal_id_len);
      if (retval) {
	goto plink_ret_1;
      }
    }

    retval = load_fam(famname, fam_cols, uii, missing_pheno, (misc_flags / MISC_AFFECTION_01) & 1, &unfiltered_sample_ct, &sample_ids, &max_sample_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_nm, &sex_male, &affection, &pheno_nm, &pheno_c, &pheno_d, &founder_info, &sample_exclude, allow_no_samples);
    if (retval) {
      goto plink_ret_1;
    }

    unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
    unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);

    if (misc_flags & MISC_MAKE_FOUNDERS_FIRST) {
      if (make_founders(unfiltered_sample_ct, unfiltered_sample_ct, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, sample_exclude, founder_info)) {
	goto plink_ret_NOMEM;
      }
    }

    if ((pheno_modifier & PHENO_MERGE) && pheno_all) {
      if (aligned_malloc(unfiltered_sample_ctl * sizeof(intptr_t), &orig_pheno_nm)) {
	goto plink_ret_NOMEM;
      }
      memcpy(orig_pheno_nm, pheno_nm, unfiltered_sample_ctl * sizeof(intptr_t));
      if (pheno_c) {
	if (aligned_malloc(unfiltered_sample_ctl * sizeof(intptr_t), &orig_pheno_c)) {
	  goto plink_ret_NOMEM;
	}
	memcpy(orig_pheno_c, pheno_c, unfiltered_sample_ctl * sizeof(intptr_t));
      } else if (pheno_d) {
	orig_pheno_d = (double*)malloc(unfiltered_sample_ct * sizeof(double));
	if (!orig_pheno_d) {
	  goto plink_ret_NOMEM;
	}
	memcpy(orig_pheno_d, pheno_d, unfiltered_sample_ct * sizeof(double));
      }
    }
    count_genders(sex_nm, sex_male, sample_exclude, unfiltered_sample_ct, &uii, &ujj, &gender_unk_ct);
    if (gender_unk_ct) {
      LOGPRINTF("%" PRIuPTR " %s (%u male%s, %u female%s, %u ambiguous) loaded from .fam.\n", unfiltered_sample_ct, species_str(unfiltered_sample_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s", gender_unk_ct);
      retval = write_nosex(outname, outname_end, unfiltered_sample_ct, sample_exclude, sex_nm, gender_unk_ct, sample_ids, max_sample_id_len);
      if (retval) {
	goto plink_ret_1;
      }
    } else {
      LOGPRINTF("%" PRIuPTR " %s (%d male%s, %d female%s) loaded from .fam.\n", unfiltered_sample_ct, species_str(unfiltered_sample_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s");
    }
    if (unfiltered_sample_ct) {
      uii = popcount_longs(pheno_nm, unfiltered_sample_ctl);
      if (uii) {
	LOGPRINTF("%u phenotype value%s loaded from .fam.\n", uii, (uii == 1)? "" : "s");
      }

      if (phenoname && fopen_checked(phenoname, "r", &phenofile)) {
	goto plink_ret_OPEN_FAIL;
      }

      if (phenofile || update_ids_fname || update_parents_fname || update_sex_fname || (filter_flags & FILTER_TAIL_PHENO)) {
	bigstack_mark = g_bigstack_base;
	retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, 0, sample_ids, max_sample_id_len, 0, 0, strcmp_deref, &cptr, &uiptr);
	if (retval) {
	  goto plink_ret_1;
	}

	if (makepheno_str) {
	  retval = makepheno_load(phenofile, makepheno_str, unfiltered_sample_ct, cptr, max_sample_id_len, uiptr, pheno_nm, &pheno_c);
	  if (retval) {
	    goto plink_ret_1;
	  }
	} else if (phenofile) {
	  retval = load_pheno(phenofile, unfiltered_sample_ct, 0, cptr, max_sample_id_len, uiptr, missing_pheno, (misc_flags / MISC_AFFECTION_01) & 1, mpheno_col, phenoname_str, pheno_nm, &pheno_c, &pheno_d, nullptr, 0);
	  if (retval) {
	    if (retval == LOAD_PHENO_LAST_COL) {
	      logerrprintb();
	      retval = RET_INVALID_FORMAT;
	      bigstack_reset(bigstack_mark);
	    }
	    goto plink_ret_1;
	  }
	}
	if (filter_flags & FILTER_TAIL_PHENO) {
	  retval = convert_tail_pheno(unfiltered_sample_ct, pheno_nm, &pheno_c, &pheno_d, tail_bottom, tail_top, missing_phenod);
	  if (retval) {
	    goto plink_ret_1;
	  }
	}
	bigstack_reset(bigstack_mark);
      }

      if (pheno_c) {
	if (calculation_type & (CALC_REGRESS_REL | CALC_REGRESS_DISTANCE | CALC_UNRELATED_HERITABILITY | CALC_GXE)) {
	  if (calculation_type & CALC_REGRESS_REL) {
	    logerrprint("Error: --regress-rel calculation requires a scalar phenotype.\n");
	  } else if (calculation_type & CALC_REGRESS_DISTANCE) {
	    logerrprint("Error: --regress-distance calculation requires a scalar phenotype.\n");
	  } else if (calculation_type & CALC_UNRELATED_HERITABILITY) {
	    logerrprint("Error: --unrelated-heritability requires a scalar phenotype.\n");
	  } else if (calculation_type & CALC_GXE) {
	    logerrprint("Error: --gxe requires a scalar phenotype.\n");
	  }
	  goto plink_ret_INVALID_CMDLINE;
	}
      } else {
	if (calculation_type & CALC_CLUSTER) {
	  if (cluster_ptr->modifier & CLUSTER_CC) {
	    logerrprint("Error: --cc requires a case/control phenotype.\n");
	    goto plink_ret_INVALID_CMDLINE;
	  } else if ((cluster_ptr->max_cases != 0xffffffffU) || (cluster_ptr->max_ctrls != 0xffffffffU)) {
	    logerrprint("Error: --mcc requires a case/control phenotype.\n");
	    goto plink_ret_INVALID_CMDLINE;
	  }
	} else if ((calculation_type & CALC_EPI) && (epi_ip->modifier & EPI_FAST)) {
	  logerrprint("Error: --fast-epistasis requires a case/control phenotype.\n");
	  goto plink_ret_INVALID_CMDLINE;
	} else if (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_FLIPSCAN)) {
	  if (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST)) {
	    logerrprint("Error: --ibs-test and --groupdist calculations require a case/control\nphenotype.\n");
	  } else if (calculation_type & CALC_FLIPSCAN) {
	    logerrprint("Error: --flip-scan requires a case/control phenotype.\n");
	  }
	  goto plink_ret_INVALID_CMDLINE;
	} else if ((calculation_type & CALC_RECODE) && (recode_modifier & (RECODE_HV | RECODE_HV_1CHR))) {
	  logerrprint("Error: --recode HV{-1chr} requires a case/control phenotype.\n");
	  goto plink_ret_INVALID_CMDLINE;
	} else if ((calculation_type & CALC_FST) && (misc_flags & MISC_FST_CC)) {
	  logerrprint("Error: \"--fst case-control\" requires a case/control phenotype.\n");
	  goto plink_ret_INVALID_CMDLINE;
	} else if ((calculation_type & CALC_FREQ) && (misc_flags & MISC_FREQ_CC)) {
	  logerrprint("Error: \"--freq case-control\" requires a case/control phenotype.\n");
	  goto plink_ret_INVALID_CMDLINE;
	}
      }

      if (!pheno_all) {
	if (loop_assoc_fname || (!pheno_d)) {
	  if ((calculation_type & CALC_GLM) && (!(glm_modifier & GLM_LOGISTIC))) {
	    logerrprint("Error: --linear without --all-pheno requires a scalar phenotype.\n");
	    goto plink_ret_INVALID_CMDLINE;
	  } else if (calculation_type & CALC_QFAM) {
	    logerrprint("Error: QFAM test requires a scalar phenotype.\n");
	    goto plink_ret_INVALID_CMDLINE;
	  }
	} else if (!pheno_c) {
	  if ((calculation_type & CALC_MODEL) && (!(model_modifier & MODEL_ASSOC))) {
	    logerrprint("Error: --model requires a case/control phenotype.\n");
	    goto plink_ret_INVALID_CMDLINE;
	  } else if ((calculation_type & CALC_GLM) && (glm_modifier & GLM_LOGISTIC)) {
	    logerrprint("Error: --logistic without --all-pheno requires a case/control phenotype.\n");
	    goto plink_ret_INVALID_CMDLINE;
	  } else if (calculation_type & (CALC_CMH | CALC_HOMOG | CALC_TESTMISS | CALC_TDT | CALC_DFAM)) {
	    if (calculation_type & CALC_CMH) {
	      logerrprint("Error: --mh and --mh2 require a case/control phenotype.\n");
	    } else if (calculation_type & CALC_HOMOG) {
	      logerrprint("Error: --homog requires a case/control phenotype.\n");
	    } else if (calculation_type & CALC_TESTMISS) {
	      logerrprint("Error: --test-missing requires a case/control phenotype.\n");
	    } else if (calculation_type & CALC_TDT) {
	      logerrprint("Error: --tdt requires a case/control phenotype.\n");
	    } else {
	      logerrprint("Error: --dfam requires a case/control phenotype.\n");
	    }
	    goto plink_ret_INVALID_CMDLINE;
	  }
	}
      }
    }
  }

  if (cm_map_fname && unfiltered_marker_ct) {
    // need sorted bps, but not marker IDs
    if (map_is_unsorted & UNSORTED_BP) {
      logerrprint("Error: --cm-map requires a sorted .bim file.  Retry this command after using\n--make-bed to sort your data.\n");
      goto plink_ret_INVALID_FORMAT;
    }
    retval = apply_cm_map(cm_map_fname, cm_map_chrname, unfiltered_marker_ct, marker_exclude, marker_pos, marker_cms, chrom_info_ptr);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (unfiltered_marker_ct != marker_exclude_ct) {
    // bugfix (12 Oct 2017): Previous plink 1.9 implementation of
    // --{exclude-}snps could affect unwanted variants at the same position as
    // named variant(s).  Backport plink 2.0's safe implementation.  This has
    // the drawback of sometimes increasing the overhead associated with
    // --snps, but we can live with that.
    uii = update_cm || update_map || update_name || (marker_alleles_needed && (update_alleles_fname || (flip_fname && (!flip_subset_fname)))) || filter_attrib_fname || qual_filter || snps_range_list_ptr->names;
    if (uii || extractname || excludename) {
      // only permit duplicate marker IDs for --extract/--exclude
      bigstack_mark = g_bigstack_base;
      retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_blen, !uii, &marker_id_htable_size, &marker_id_htable);
      if (retval) {
	goto plink_ret_1;
      }
      if (snps_range_list_ptr->names) {
        retval = snps_flag(marker_ids, marker_id_htable, snps_range_list_ptr, unfiltered_marker_ct, max_marker_id_blen, marker_id_htable_size, (filter_flags / FILTER_EXCLUDE_MARKERNAME_SNP) & 1, marker_exclude, &marker_exclude_ct);
        if (retval) {
          goto plink_ret_1;
        }
      }
      if (update_cm) {
	retval = update_marker_cms(update_cm, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_cms);
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if (update_map) {
	retval = update_marker_pos(update_map, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, marker_pos, &map_is_unsorted, chrom_info_ptr);
	if (retval) {
	  goto plink_ret_1;
	}
      } else if (update_name) {
	retval = update_marker_names(update_name, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_blen, unfiltered_marker_ct);
	if (retval) {
	  goto plink_ret_1;
	}
	if (update_alleles_fname || (marker_alleles_needed && flip_fname && (!flip_subset_fname)) || extractname || excludename) {
	  bigstack_reset(bigstack_mark);
	  retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_blen, 0, &marker_id_htable_size, &marker_id_htable);
	  if (retval) {
	    goto plink_ret_1;
	  }
	}
      }
      if (marker_alleles_needed) {
	if (update_alleles_fname) {
	  retval = update_marker_alleles(update_alleles_fname, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, marker_allele_ptrs, &max_marker_allele_blen, outname, outname_end);
	  if (retval) {
	    goto plink_ret_1;
	  }
	}
	if (flip_fname && (!flip_subset_fname)) {
	  retval = flip_strand(flip_fname, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, marker_allele_ptrs);
	  if (retval) {
	    goto plink_ret_1;
	  }
	}
      }
      if (extractname) {
	if (!(misc_flags & MISC_EXTRACT_RANGE)) {
	  retval = extract_exclude_flag_norange(extractname, marker_id_htable, marker_id_htable_size, 0, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, allow_no_variants);
	  if (retval) {
	    goto plink_ret_1;
	  }
	} else {
	  if (map_is_unsorted & UNSORTED_BP) {
	    logerrprint("Error: \"--extract range\" requires a sorted .bim.  Retry this command after\nusing --make-bed to sort your data.\n");
	    goto plink_ret_INVALID_CMDLINE;
	  }
	  retval = extract_exclude_range(extractname, marker_pos, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0, allow_no_variants, chrom_info_ptr);
	  if (retval) {
	    goto plink_ret_1;
	  }
	  uljj = unfiltered_marker_ct - marker_exclude_ct;
	  LOGPRINTF("--extract range: %" PRIuPTR " variant%s remaining.\n", uljj, (uljj == 1)? "" : "s");
	}
      }
      if (excludename) {
	if (!(misc_flags & MISC_EXCLUDE_RANGE)) {
	  retval = extract_exclude_flag_norange(excludename, marker_id_htable, marker_id_htable_size, 1, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, allow_no_variants);
	  if (retval) {
	    goto plink_ret_1;
	  }
	} else {
	  if (map_is_unsorted & UNSORTED_BP) {
	    logerrprint("Error: \"--exclude range\" requires a sorted .bim.  Retry this command after\nusing --make-bed to sort your data.\n");
	    goto plink_ret_INVALID_CMDLINE;
	  }
	  retval = extract_exclude_range(excludename, marker_pos, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 1, allow_no_variants, chrom_info_ptr);
	  if (retval) {
	    goto plink_ret_1;
	  }
	  uljj = unfiltered_marker_ct - marker_exclude_ct;
	  LOGPRINTF("--exclude range: %" PRIuPTR " variant%s remaining.\n", uljj, (uljj == 1)? "" : "s");
	}
      }
      if (filter_attrib_fname) {
	retval = filter_attrib(filter_attrib_fname, filter_attrib_liststr, marker_id_htable, marker_id_htable_size, allow_no_variants, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if (qual_filter) {
	retval = filter_qual_scores(qual_filter, qual_min_thresh, qual_max_thresh, marker_id_htable, marker_id_htable_size, allow_no_variants, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
	if (retval) {
	  goto plink_ret_1;
	}
      }
      bigstack_reset(bigstack_mark);
    }

    if (allelexxxx) {
      allelexxxx_recode(allelexxxx, marker_allele_ptrs, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct);
    }

    if (thin_keep_prob != 1.0) {
      if (random_thin_markers(thin_keep_prob, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, allow_no_variants)) {
	goto plink_ret_ALL_MARKERS_EXCLUDED;
      }
    } else if (thin_keep_ct != 0xffffffffU) {
      retval = random_thin_markers_ct(thin_keep_ct, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
      if (retval) {
	goto plink_ret_1;
      }
    }
  }

  if (bedfile) {
    if (fseeko(bedfile, 0, SEEK_END)) {
      goto plink_ret_READ_FAIL;
    }
    llxx = ftello(bedfile);
    if (!llxx) {
      logerrprint("Error: Empty .bed file.\n");
      goto plink_ret_INVALID_FORMAT;
    }
    rewind(bedfile);
    uii = fread(g_textbuf, 1, 3, bedfile);
    llyy = ((uint64_t)unfiltered_sample_ct4) * unfiltered_marker_ct;
    llzz = ((uint64_t)unfiltered_sample_ct) * ((unfiltered_marker_ct + 3) / 4);
    if ((uii == 3) && (!memcmp(g_textbuf, "l\x1b\x01", 3))) {
      llyy += 3;
    } else if ((uii == 3) && (!memcmp(g_textbuf, "l\x1b", 3))) {
      // v1.00 sample-major
      llyy = llzz + 3;
      bed_offset = 2;
    } else if (uii && (*g_textbuf == '\x01')) {
      // v0.99 SNP-major
      llyy += 1;
      bed_offset = 1;
    } else if (uii && (!(*g_textbuf))) {
      // v0.99 sample-major
      llyy = llzz + 1;
      bed_offset = 2;
    } else {
      // pre-v0.99, sample-major, no header bytes
      if (llxx != llzz) {
	// probably not PLINK-format at all, so give this error instead of
	// "invalid file size"
	logerrprint("Error: Invalid header bytes in .bed file.\n");
	goto plink_ret_INVALID_FORMAT;
      }
      llyy = llzz;
      bed_offset = 2;
    }
    if (llxx != llyy) {
      if ((*g_textbuf == '#') || ((uii == 3) && (!memcmp(g_textbuf, "chr", 3)))) {
	logerrprint("Error: Invalid header bytes in PLINK 1 .bed file.  (Is this a UCSC Genome\nBrowser BED file instead?)\n");
	goto plink_ret_INVALID_FORMAT;
      } else {
	sprintf(g_logbuf, "Error: Invalid .bed file size (expected %" PRId64 " bytes).\n", llyy);
	goto plink_ret_INVALID_FORMAT_2;
      }
    }

    if (bed_offset == 2) {
      strcpy(outname_end, ".bed.vmaj"); // not really temporary
      logprint("Sample-major .bed file detected.  Transposing...\n");
      fclose(bedfile);
      retval = sample_major_to_snp_major(bedname, outname, unfiltered_marker_ct, unfiltered_sample_ct, llxx);
      if (retval) {
	goto plink_ret_1;
      }
      LOGPRINTFWW("Variant-major .bed written to %s .\n", outname);
      strcpy(bedname, outname);
      if (fopen_checked(bedname, FOPEN_RB, &bedfile)) {
	goto plink_ret_OPEN_FAIL;
      }
      bed_offset = 3;
    }
  }

  if (unfiltered_sample_ct && (update_ids_fname || update_parents_fname || update_sex_fname || keepname || keepfamname || removename || removefamname || filter_attrib_sample_fname || om_ip->marker_fname || filtername)) {
    bigstack_mark = g_bigstack_base;
    retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, 0, sample_ids, max_sample_id_len, 0, 0, strcmp_deref, &cptr, &uiptr);
    if (retval) {
      goto plink_ret_1;
    }
    if (update_ids_fname) {
      retval = update_sample_ids(update_ids_fname, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, sample_ids);
      if (retval) {
	goto plink_ret_1;
      }
      bigstack_reset(g_bigstack_base);
      retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, 0, sample_ids, max_sample_id_len, 0, 0, strcmp_deref, &cptr, &uiptr);
      if (retval) {
	goto plink_ret_1;
      }
    } else {
      if (update_parents_fname) {
	retval = update_sample_parents(update_parents_fname, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if (update_sex_fname) {
        retval = update_sample_sexes(update_sex_fname, update_sex_col, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, sex_nm, sex_male);
	if (retval) {
	  goto plink_ret_1;
	}
      }
    }
    // sample_exclude_ct assumed to be 0 before this point
    if (keepfamname) {
      retval = keep_or_remove(keepfamname, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, 2, allow_no_samples);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (keepname) {
      retval = keep_or_remove(keepname, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, 0, allow_no_samples);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (removefamname) {
      retval = keep_or_remove(removefamname, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, 3, allow_no_samples);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (removename) {
      retval = keep_or_remove(removename, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, 1, allow_no_samples);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (filter_attrib_sample_fname) {
      retval = filter_attrib_sample(filter_attrib_sample_fname, filter_attrib_sample_liststr, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, unfiltered_sample_ct, allow_no_samples, sample_exclude, &sample_exclude_ct);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (om_ip->marker_fname) {
      // would rather do this with pre-sorted markers, but that might break
      // order-of-operations assumptions in existing pipelines
      if ((unfiltered_sample_ct == sample_exclude_ct) || (unfiltered_marker_ct == marker_exclude_ct)) {
	// don't need this if everything that refers to om_ip is skipped
        oblig_missing_cleanup(om_ip);
	om_ip->cluster_ct = 0;
	om_ip->entry_ct = 0;
      } else {
	retval = load_oblig_missing(bedfile, bed_offset, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_blen, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, unfiltered_sample_ct, sex_male, chrom_info_ptr, om_ip);
	if (retval) {
	  goto plink_ret_1;
	}
      }
    }
    if (filtername) {
      if (!mfilter_col) {
	mfilter_col = 1;
      }
      retval = filter_samples_file(filtername, cptr, unfiltered_sample_ct, max_sample_id_len, uiptr, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, filtervals_flattened, mfilter_col, allow_no_samples);
      if (retval) {
	goto plink_ret_1;
      }
    }
    bigstack_reset(bigstack_mark);
  }

  if (famname[0] && (unfiltered_sample_ct != sample_exclude_ct)) {
    if (gender_unk_ct && (!(sex_missing_pheno & ALLOW_NO_SEX))) {
      uii = popcount_longs_exclude(pheno_nm, sex_nm, unfiltered_sample_ctl);
      if (uii) {
	if ((!sex_missing_pheno) && (calculation_type & (CALC_MAKE_BED | CALC_MAKE_FAM | CALC_RECODE | CALC_WRITE_COVAR))) {
	  if (calculation_type & (~(CALC_MAKE_BED | CALC_MAKE_FAM | CALC_RECODE | CALC_WRITE_COVAR))) {
	    logerrprint("Error: When ambiguous-sex samples with phenotype data are present,\n--make-bed/--make-just-fam/--recode/--write-covar usually cannot be combined\nwith other commands.  Split them across multiple PLINK runs, or use\n--allow-no-sex/--must-have-sex.\n");
	    goto plink_ret_INVALID_CMDLINE;
	  }
	} else {
	  // either --must-have-sex without --allow-no-sex, or no data
	  // generation command
	  bitvec_and(sex_nm, unfiltered_sample_ctl, pheno_nm);
	}
      }
      if (uii || pheno_all || loop_assoc_fname) {
	logerrprint("Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those\nphenotypes to be ignored, use the --allow-no-sex flag.\n");
      }
    }
    if (filter_flags & FILTER_PRUNE) {
      bitvec_ornot(pheno_nm, unfiltered_sample_ctl, sample_exclude);
      zero_trailing_bits(unfiltered_sample_ct, sample_exclude);
      sample_exclude_ct = popcount_longs(sample_exclude, unfiltered_sample_ctl);
      if ((sample_exclude_ct == unfiltered_sample_ct) && (!allow_no_samples)) {
	LOGERRPRINTF("Error: All %s removed by --prune.\n", g_species_plural);
	goto plink_ret_ALL_SAMPLES_EXCLUDED;
      }
      LOGPRINTF("--prune: %" PRIuPTR " %s remaining.\n", unfiltered_sample_ct - sample_exclude_ct, species_str(unfiltered_sample_ct - sample_exclude_ct));
    }

    if (filter_flags & (FILTER_BINARY_CASES | FILTER_BINARY_CONTROLS)) {
      if (!pheno_c) {
	logerrprint("Error: --filter-cases/--filter-controls requires a case/control phenotype.\n");
	goto plink_ret_INVALID_CMDLINE;
      }
      ii = sample_exclude_ct;
      // fcc == 1: exclude all zeroes in pheno_c
      // fcc == 2: exclude all ones in pheno_c
      // -> flip on fcc == 1
      filter_samples_bitfields(unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, pheno_c, (filter_flags / FILTER_BINARY_CASES) & 1, pheno_nm);
      if ((sample_exclude_ct == unfiltered_sample_ct) && (!allow_no_samples)) {
	LOGERRPRINTF("Error: All %s removed due to case/control status (--filter-%s).\n", g_species_plural, (filter_flags & FILTER_BINARY_CASES)? "cases" : "controls");
	goto plink_ret_ALL_SAMPLES_EXCLUDED;
      }
      ii = sample_exclude_ct - ii;
      LOGPRINTF("%d %s removed due to case/control status (--filter-%s).\n", ii, species_str(ii), (filter_flags & FILTER_BINARY_CASES)? "cases" : "controls");
    }
    if (filter_flags & (FILTER_BINARY_FEMALES | FILTER_BINARY_MALES)) {
      ii = sample_exclude_ct;
      filter_samples_bitfields(unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, sex_male, (filter_flags / FILTER_BINARY_MALES) & 1, sex_nm);
      if ((sample_exclude_ct == unfiltered_sample_ct) && (!allow_no_samples)) {
	LOGERRPRINTF("Error: All %s removed due to gender filter (--filter-%s).\n", g_species_plural, (filter_flags & FILTER_BINARY_MALES)? "males" : "females");
	goto plink_ret_ALL_SAMPLES_EXCLUDED;
      }
      ii = sample_exclude_ct - ii;
      LOGPRINTF("%d %s removed due to gender filter (--filter-%s).\n", ii, species_str(ii), (filter_flags & FILTER_BINARY_MALES)? "males" : "females");
    }
    if (filter_flags & (FILTER_BINARY_FOUNDERS | FILTER_BINARY_NONFOUNDERS)) {
      ii = sample_exclude_ct;
      filter_samples_bitfields(unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, founder_info, (filter_flags / FILTER_BINARY_FOUNDERS) & 1, nullptr);
      if ((sample_exclude_ct == unfiltered_sample_ct) && (!allow_no_samples)) {
	LOGERRPRINTF("Error: All %s removed due to founder status (--filter-%s).\n", g_species_plural, (filter_flags & FILTER_BINARY_FOUNDERS)? "founders" : "nonfounders");
	goto plink_ret_ALL_SAMPLES_EXCLUDED;
      }
      ii = sample_exclude_ct - ii;
      LOGPRINTF("%d %s removed due to founder status (--filter-%s).\n", ii, species_str(ii), (filter_flags & FILTER_BINARY_FOUNDERS)? "founders" : "nonfounders");
    }

    if (thin_keep_sample_prob != 1.0) {
      if (random_thin_samples(thin_keep_sample_prob, unfiltered_sample_ct, allow_no_samples, sample_exclude, &sample_exclude_ct)) {
        goto plink_ret_ALL_SAMPLES_EXCLUDED;
      }
    } else if (thin_keep_sample_ct != 0xffffffffU) {
      retval = random_thin_samples_ct(thin_keep_sample_ct, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct);
      if (retval) {
        goto plink_ret_1;
      }
    }

    if (mind_thresh < 1.0) {
      retval = mind_filter(bedfile, bed_offset, outname, outname_end, mind_thresh, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, sample_ids, max_sample_id_len, sex_male, chrom_info_ptr, om_ip, allow_no_samples);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (cluster_ptr->fname || (misc_flags & MISC_FAMILY_CLUSTERS)) {
      // could save off bigstack_mark here and free immediately after
      // load_clusters(), if clusters are *only* used for filtering.  But not a
      // big deal.
      retval = load_clusters(cluster_ptr->fname, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, sample_ids, max_sample_id_len, mwithin_col, (misc_flags / MISC_LOAD_CLUSTER_KEEP_NA) & 1, &cluster_ct, &cluster_map, &cluster_starts, &cluster_ids, &max_cluster_id_len, cluster_ptr->keep_fname, cluster_ptr->keep_flattened, cluster_ptr->remove_fname, cluster_ptr->remove_flattened, allow_no_samples);
      if (retval) {
	goto plink_ret_1;
      }
    }
    sample_ct = unfiltered_sample_ct - sample_exclude_ct;
    if ((!sample_ct) && (!allow_no_samples)) {
      // defensive; currently shouldn't happen since we're checking at every
      // filter
      LOGERRPRINTF("Error: No %s pass QC.\n", g_species_plural);
      goto plink_ret_ALL_SAMPLES_EXCLUDED;
    }

    if ((sample_ct < 2) && (relationship_or_ibc_req(calculation_type) || distance_req(read_dists_fname, calculation_type) || (calculation_type & (CALC_GENOME | CALC_CLUSTER | CALC_NEIGHBOR)))) {
      sprintf(g_logbuf, "Error: At least 2 %s required for pairwise analysis.\n", g_species_plural);
      goto plink_ret_INVALID_CMDLINE_2;
    }

    if ((parallel_tot > 1) && (calculation_type & (CALC_DISTANCE | CALC_GENOME | CALC_RELATIONSHIP)) && (parallel_tot > sample_ct / 2)) {
      sprintf(g_logbuf, "Error: Too many --parallel jobs (maximum %" PRIuPTR "/2 = %" PRIuPTR ").\n", sample_ct, sample_ct / 2);
      goto plink_ret_INVALID_CMDLINE_2;
    }
  }
#ifndef NOLAPACK
  // Multithreaded BLAS/LAPACK call?  If yes, and either user requested <=
  // [0.5 * nprocs] threads or nprocs is unknown, warn that BLAS/LAPACK
  // multithreading is not under PLINK's control.
  uii = ((!known_procs) || (known_procs >= 2 * ((int32_t)g_thread_ct))) && ((calculation_type & (CALC_LASSO | CALC_PCA | CALC_UNRELATED_HERITABILITY)) || ((calculation_type & CALC_GLM) && pheno_d) || cluster_ptr->mds_dim_ct || ((calculation_type & CALC_LD_PRUNE) && (!(ldip->modifier & (LD_PRUNE_PAIRWISE | LD_PRUNE_PAIRPHASE)))));
#endif
  if (g_thread_ct > 1) {
    if ((calculation_type & (CALC_RELATIONSHIP | CALC_REL_CUTOFF | CALC_GDISTANCE_MASK | CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE | CALC_GENOME | CALC_REGRESS_REL | CALC_UNRELATED_HERITABILITY | CALC_LD | CALC_PCA | CALC_MAKE_PERM_PHENO | CALC_QFAM)) || ((calculation_type & CALC_MODEL) && (model_modifier & (MODEL_PERM | MODEL_MPERM))) || ((calculation_type & CALC_GLM) && (glm_modifier & (GLM_PERM | GLM_MPERM))) || ((calculation_type & CALC_TESTMISS) && (testmiss_modifier & (TESTMISS_PERM | TESTMISS_MPERM))) || ((calculation_type & CALC_TDT) && (fam_ip->tdt_modifier & (TDT_PERM | TDT_MPERM))) || ((calculation_type & CALC_DFAM) && (fam_ip->dfam_modifier & (DFAM_PERM | DFAM_MPERM))) || ((calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) && (!read_genome_fname) && ((cluster_ptr->ppc != 0.0) || (!read_dists_fname))) || ((calculation_type & CALC_EPI) && (epi_ip->modifier & (EPI_FAST | EPI_REG)))
#ifndef _WIN32
        || ((calculation_type & CALC_FREQ) && (misc_flags & MISC_FREQ_GZ))
        || ((calculation_type & CALC_MISSING_REPORT) && (misc_flags & MISC_MISSING_GZ))
        || ((calculation_type & CALC_HARDY) && (hwe_modifier & HWE_GZ))
        || ((calculation_type & CALC_HET) && (misc_flags & MISC_HET_GZ))
	|| ((calculation_type & CALC_RECODE) && (((recode_modifier & (RECODE_VCF | RECODE_BGZ)) == (RECODE_VCF | RECODE_BGZ))))
#endif
) {
      LOGPRINTF("Using up to %u threads (change this with --threads).\n", g_thread_ct);
    } else {
#ifndef NOLAPACK
      if (uii) {
	logprint("Using 1 thread.\n");
      } else {
#endif
	logprint("Using 1 thread (no multithreaded calculations invoked).\n");
#ifndef NOLAPACK
      }
#endif
    }
  } else {
    logprint("Using 1 thread.\n");
  }
#ifndef NOLAPACK
  if (uii) {
    logerrprint("Warning: This run includes BLAS/LAPACK linear algebra operations which\ncurrently disregard the --threads limit.  If this is problematic, you may want\nto recompile against single-threaded BLAS/LAPACK.\n");
  }
#endif

  if (famname[0]) {
    if ((calculation_type & (CALC_SEXCHECK | CALC_MISSING_REPORT | CALC_GENOME | CALC_HOMOZYG | CALC_SCORE | CALC_MENDEL | CALC_HET)) || cluster_ptr->mds_dim_ct) {
      calc_plink_maxfid(unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, &plink_maxfid, &plink_maxiid);
    }

    if (sample_sort & (SAMPLE_SORT_NATURAL | SAMPLE_SORT_ASCII)) {
      retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, sample_exclude_ct, sample_ids, max_sample_id_len, 0, 0, (sample_sort & SAMPLE_SORT_NATURAL)? strcmp_natural_deref : strcmp_deref, &cptr, &uiptr);
      if (retval) {
	goto plink_ret_1;
      }
      sample_sort_map = uiptr;
      bigstack_reset(cptr);
    } else if (sample_sort == SAMPLE_SORT_FILE) {
      retval = sample_sort_file_map(sample_sort_fname, unfiltered_sample_ct, sample_exclude, unfiltered_sample_ct - sample_exclude_ct, sample_ids, max_sample_id_len, &sample_sort_map);
      if (retval) {
	goto plink_ret_1;
      }
    }

    if ((filter_flags & FILTER_MAKE_FOUNDERS) && (!(misc_flags & MISC_MAKE_FOUNDERS_FIRST))) {
      if (make_founders(unfiltered_sample_ct, sample_ct, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, sample_exclude, founder_info)) {
	goto plink_ret_NOMEM;
      }
    }

    if (calculation_type & CALC_WRITE_CLUSTER) {
      retval = write_clusters(outname, outname_end, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, (misc_flags / MISC_WRITE_CLUSTER_OMIT_UNASSIGNED) & 1, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len);
      if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER))))) {
	goto plink_ret_1;
      }
    }

    // this currently has to come last since covar data structures refer to
    // filtered sample indices.
    if (covar_fname) {
      // update this as more covariate-referencing commands are added
      if (!(calculation_type & (CALC_MAKE_BED | CALC_MAKE_FAM | CALC_RECODE | CALC_WRITE_COVAR | CALC_GXE | CALC_GLM | CALC_LASSO | CALC_RPLUGIN))) {
	logerrprint("Warning: Ignoring --covar since no commands reference the covariates.\n");
      } else {
	// if only --gxe, ignore --covar-name/--covar-number
	uii = (calculation_type & (CALC_MAKE_BED | CALC_MAKE_FAM | CALC_RECODE | CALC_WRITE_COVAR | CALC_GLM | CALC_LASSO | CALC_RPLUGIN))? 1 : 0;
	retval = load_covars(covar_fname, unfiltered_sample_ct, sample_exclude, sample_ct, nullptr, nullptr, sample_ids, max_sample_id_len, missing_phenod, uii? covar_modifier : (covar_modifier & COVAR_KEEP_PHENO_ON_MISSING_COV), uii? covar_range_list_ptr : nullptr, gxe_mcovar, &covar_ct, &covar_names, &max_covar_name_len, pheno_nm, &covar_nm, &covar_d, &gxe_covar_nm, &gxe_covar_c);
	if (retval) {
	  goto plink_ret_1;
	}
      }
    }

    bitvec_andnot(sample_exclude, unfiltered_sample_ctl, pheno_nm);
    if (pheno_c) {
      bitvec_and(pheno_nm, unfiltered_sample_ctl, pheno_c);
    }
    bitvec_andnot(sample_exclude, unfiltered_sample_ctl, founder_info);
    bitvec_andnot(sample_exclude, unfiltered_sample_ctl, sex_nm);
    if (gender_unk_ct) {
      gender_unk_ct = sample_ct - popcount_longs(sex_nm, unfiltered_sample_ctl);
    }
    bitvec_and(sex_nm, unfiltered_sample_ctl, sex_male);

    pheno_nm_ct = popcount_longs(pheno_nm, unfiltered_sample_ctl);
    if (!pheno_nm_ct) {
      hwe_modifier |= HWE_THRESH_ALL;
    } else if (pheno_c) {
      pheno_ctrl_ct = pheno_nm_ct - popcount_longs(pheno_c, unfiltered_sample_ctl);
      if (!pheno_ctrl_ct) {
	hwe_modifier |= HWE_THRESH_ALL;
      }
    } else {
      hwe_modifier |= HWE_THRESH_ALL;
    }
    ulii = popcount_longs(founder_info, unfiltered_sample_ctl);
    LOGPRINTFWW("Before main variant filters, %" PRIuPTR " founder%s and %" PRIuPTR " nonfounder%s present.\n", ulii, (ulii == 1)? "" : "s", sample_ct - ulii, (sample_ct - ulii == 1)? "" : "s");

  }

  if (bimname[0] && (unfiltered_marker_ct != marker_exclude_ct)) {
    plink_maxsnp = calc_plink_maxsnp(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_blen);
    uii = BITCT_TO_WORDCT(unfiltered_marker_ct);
    if (bigstack_calloc_ul(uii, &marker_reverse)) {
      goto plink_ret_NOMEM;
    }
    if (bedfile && sample_ct) {
      retval = calc_freqs_and_hwe(bedfile, outname, outname_end, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_blen, unfiltered_sample_ct, sample_exclude, sample_exclude_ct, sample_ids, max_sample_id_len, founder_info, nonfounders, (misc_flags / MISC_MAF_SUCC) & 1, set_allele_freqs, bed_offset, (hwe_thresh > 0.0) || (calculation_type & CALC_HARDY), hwe_modifier & HWE_THRESH_ALL, (pheno_nm_ct && pheno_c)? ((calculation_type / CALC_HARDY) & 1) : 0, min_ac, max_ac, geno_thresh, pheno_nm, pheno_nm_ct? pheno_c : nullptr, &hwe_lls, &hwe_lhs, &hwe_hhs, &hwe_ll_cases, &hwe_lh_cases, &hwe_hh_cases, &hwe_ll_allfs, &hwe_lh_allfs, &hwe_hh_allfs, &hwe_hapl_allfs, &hwe_haph_allfs, &geno_excl_bitfield, &ac_excl_bitfield, &sample_male_ct, &sample_f_ct, &sample_f_male_ct, chrom_info_ptr, om_ip, sex_nm, sex_male, map_is_unsorted & UNSORTED_SPLIT_CHROM, &hh_exists);
      if (retval) {
	goto plink_ret_1;
      }
    } else {
      for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
        set_allele_freqs[marker_uidx] = 1.0;
      }
    }

    if (freqname) {
      retval = read_external_freqs(freqname, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_blen, chrom_info_ptr, marker_allele_ptrs, set_allele_freqs, nchrobs, (misc_flags / MISC_MAF_SUCC) & 1);
      if (retval) {
	goto plink_ret_1;
      }
    }

    if (!(misc_flags & MISC_KEEP_ALLELE_ORDER)) {
      // after this, set_allele_freqs[] has A2 freqs
      calc_marker_reverse_bin(marker_reverse, marker_exclude, unfiltered_marker_ct, unfiltered_marker_ct - marker_exclude_ct, set_allele_freqs);
    }

    if (a1alleles || a2alleles) {
      retval = load_ax_alleles(a1alleles? a1alleles : a2alleles, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_allele_ptrs, &max_marker_allele_blen, marker_reverse, marker_ids, max_marker_id_blen, set_allele_freqs, a2alleles? 1 : 0);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (marker_allele_ptrs) {
      swap_reversed_marker_alleles(unfiltered_marker_ct, marker_reverse, marker_allele_ptrs);
    }

    // contrary to the PLINK 1.07 flowchart, --freq effectively resolves before
    // --geno.
    if (sample_ct) {
      if (calculation_type & CALC_FREQ) {
	if (cluster_ct && (!(misc_flags & MISC_FREQX))) {
	  if (misc_flags & MISC_FREQ_COUNTS) {
	    logprint("Note: --freq 'counts' modifier has no effect on cluster-stratified report.\n");
	  }
	  retval = write_stratified_freqs(bedfile, bed_offset, outname, outname_end, (misc_flags / MISC_FREQ_GZ) & 1, plink_maxsnp, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, marker_ids, max_marker_id_blen, marker_allele_ptrs, max_marker_allele_blen, unfiltered_sample_ct, sample_ct, sample_f_ct, founder_info, nonfounders, sex_male, sample_f_male_ct, marker_reverse, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len);
	} else if (misc_flags & MISC_FREQ_CC) {
	  retval = write_cc_freqs(bedfile, bed_offset, outname, outname_end, (misc_flags / MISC_FREQ_GZ) & 1, plink_maxsnp, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, marker_ids, max_marker_id_blen, marker_allele_ptrs, max_marker_allele_blen, unfiltered_sample_ct, founder_info, nonfounders, sex_male, marker_reverse, pheno_nm, pheno_c);
	} else {
	  retval = write_freqs(outname, outname_end, plink_maxsnp, unfiltered_marker_ct, marker_exclude, set_allele_freqs, chrom_info_ptr, marker_ids, max_marker_id_blen, marker_allele_ptrs, max_marker_allele_blen, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, hwe_hapl_allfs, hwe_haph_allfs, sample_f_ct, sample_f_male_ct, nonfounders, misc_flags, marker_reverse);
	}
	if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ))))) {
	  goto plink_ret_1;
	}
      }
      if (calculation_type & CALC_MISSING_REPORT) {
	retval = write_missingness_reports(bedfile, bed_offset, outname, outname_end, (misc_flags / MISC_MISSING_GZ) & 1, plink_maxfid, plink_maxiid, plink_maxsnp, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, chrom_info_ptr, om_ip, marker_ids, max_marker_id_blen, unfiltered_sample_ct, sample_ct, sample_exclude, pheno_nm, sex_male, sample_male_ct, sample_ids, max_sample_id_len, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len, hh_exists);
	if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ | CALC_MISSING_REPORT))))) {
	  goto plink_ret_1;
	}
      }

      if (geno_excl_bitfield) {
	ulii = marker_exclude_ct;
	uljj = BITCT_TO_WORDCT(unfiltered_marker_ct);
	bitvec_or(geno_excl_bitfield, uljj, marker_exclude);
	marker_exclude_ct = popcount_longs(marker_exclude, uljj);
	if ((marker_exclude_ct == unfiltered_marker_ct) && (!allow_no_variants)) {
	  logerrprint("Error: All variants excluded due to missing genotype data (--geno).\n");
	  goto plink_ret_ALL_MARKERS_EXCLUDED;
	}
	ulii = marker_exclude_ct - ulii;
	LOGPRINTF("%" PRIuPTR " variant%s removed due to missing genotype data (--geno).\n", ulii, (ulii == 1)? "" : "s");
      }
    }
    oblig_missing_cleanup(om_ip);
    if (sample_ct) {
      if (calculation_type & CALC_HARDY) {
	retval = hardy_report(outname, outname_end, output_min_p, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_allele_ptrs, max_marker_allele_blen, marker_reverse, hwe_lls, hwe_lhs, hwe_hhs, hwe_modifier, nonfounders, hwe_ll_cases, hwe_lh_cases, hwe_hh_cases, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, pheno_nm_ct, pheno_c, chrom_info_ptr);
	if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ | CALC_HARDY))))) {
	  goto plink_ret_1;
	}
      }
      if (hwe_thresh > 0.0) {
	if (enforce_hwe_threshold(hwe_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, hwe_lls, hwe_lhs, hwe_hhs, hwe_modifier, allow_no_variants, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, chrom_info_ptr)) {
	  goto plink_ret_ALL_MARKERS_EXCLUDED;
	}
      }
    }
    if ((min_maf != 0.0) || (max_maf != 0.5) || ac_excl_bitfield) {
      if (enforce_minor_allele_thresholds(min_maf, max_maf, unfiltered_marker_ct, marker_exclude, ac_excl_bitfield, &marker_exclude_ct, set_allele_freqs, allow_no_variants)) {
	goto plink_ret_ALL_MARKERS_EXCLUDED;
      }
    }
    if (min_bp_space && (unfiltered_marker_ct > marker_exclude_ct)) {
      if (map_is_unsorted & UNSORTED_BP) {
	logerrprint("Error: --bp-space requires a sorted .bim file.  Retry this command after using\n--make-bed to sort your data.\n");
	goto plink_ret_INVALID_FORMAT;
      }
      enforce_min_bp_space(min_bp_space, unfiltered_marker_ct, marker_exclude, marker_pos, &marker_exclude_ct, chrom_info_ptr);
    }

    if (bedfile && sample_ct && (unfiltered_marker_ct > marker_exclude_ct)) {
      if ((calculation_type & CALC_MENDEL) || (fam_ip->mendel_modifier & MENDEL_FILTER)) {
	retval = mendel_error_scan(fam_ip, bedfile, bed_offset, outname, outname_end, plink_maxfid, plink_maxiid, plink_maxsnp, allow_no_variants, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, marker_reverse, marker_ids, max_marker_id_blen, marker_allele_ptrs, max_marker_allele_blen, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, hh_exists, chrom_info_ptr, (calculation_type / CALC_MENDEL) & 1);
	if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ | CALC_MISSING_REPORT | CALC_MENDEL))))) {
	  goto plink_ret_1;
	}
	if (fam_ip->mendel_modifier & MENDEL_FILTER) {
	  // gah
	  sample_ct = unfiltered_sample_ct - sample_exclude_ct;
	  bitvec_andnot(sample_exclude, unfiltered_sample_ctl, founder_info);
	  bitvec_andnot(sample_exclude, unfiltered_sample_ctl, sex_nm);
	  bitvec_and(sex_nm, unfiltered_sample_ctl, sex_male);
	  if (pheno_nm_ct) {
	    bitvec_andnot(sample_exclude, unfiltered_sample_ctl, pheno_nm);
	    pheno_nm_ct = popcount_longs(pheno_nm, unfiltered_sample_ctl);
	    if (pheno_c) {
	      bitvec_and(pheno_nm, unfiltered_sample_ctl, pheno_c);
	      pheno_ctrl_ct = pheno_nm_ct - popcount_longs(pheno_c, unfiltered_sample_ctl);
	    }
	  }
	}
      }
      bigstack_reset(hwe_lls);
    }
    if (sip->fname) {
      if (map_is_unsorted & UNSORTED_BP) {
	logerrprint("Error: --set/--make-set requires a sorted .bim file.  Retry this command after\nusing --make-bed to sort your data.\n");
	goto plink_ret_INVALID_FORMAT;
      }
      retval = define_sets(sip, unfiltered_marker_ct, marker_exclude, marker_pos, &marker_exclude_ct, marker_ids, max_marker_id_blen, chrom_info_ptr, allow_no_variants);
      if (retval) {
	goto plink_ret_1;
      }
    }

    marker_ct = unfiltered_marker_ct - marker_exclude_ct;
    if ((!marker_ct) && (!allow_no_variants)) {
      // defensive
      logerrprint("Error: All variants fail QC.\n");
      goto plink_ret_ALL_MARKERS_EXCLUDED;
    }
    // could add --parallel + --r/--r2/--{fast-}epistasis check here.  (this
    // currently happens downstream)

    if (bedfile) {
      LOGPRINTFWW("%" PRIuPTR " variant%s and %" PRIuPTR " %s pass filters and QC%s.\n", marker_ct, (marker_ct == 1)? "" : "s", sample_ct, species_str(sample_ct), (calculation_type & CALC_REL_CUTOFF)? " (before --rel-cutoff)": "");
    } else {
      LOGPRINTFWW("%" PRIuPTR " variant%s filters and QC.\n", marker_ct, (marker_ct == 1)? " passes" : "s pass");
    }
  } else if ((!allow_no_variants) && ((calculation_type & (~CALC_ONLY_FAM)) || (filter_flags & FILTER_BIM_REQ))) {
    // defensive
    logerrprint("Error: No variants remaining.\n");
    goto plink_ret_ALL_MARKERS_EXCLUDED;
  }

  if (famname[0] && sample_ct) {
    if (!pheno_nm_ct) {
      logprint("Note: No phenotypes present.\n");
    } else if (pheno_c) {
      if (pheno_nm_ct != sample_ct) {
	sprintf(g_logbuf, "Among remaining phenotypes, %u %s and %u %s.  (%" PRIuPTR " phenotype%s missing.)\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "is a case" : "are cases", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "is a control" : "are controls", sample_ct - pheno_nm_ct, (sample_ct - pheno_nm_ct == 1)? " is" : "s are");
      } else {
	sprintf(g_logbuf, "Among remaining phenotypes, %u %s and %u %s.\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "is a case" : "are cases", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "is a control" : "are controls");
      }
      wordwrapb(0);
      logprintb();
    } else {
      logprint("Phenotype data is quantitative.\n");
    }
  }

  if (sample_ct) {
    if (marker_ct) {
      if (relationship_or_ibc_req(calculation_type)) {
	if (relip->pca_cluster_names_flattened || relip->pca_clusters_fname) {
	  retval = extract_clusters(unfiltered_sample_ct, sample_exclude, sample_ct, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len, relip->pca_cluster_names_flattened, relip->pca_clusters_fname, &pca_sample_exclude, &pca_sample_ct);
	  if (retval) {
	    goto plink_ret_1;
	  }
	  if (pca_sample_ct < 2) {
	    logerrprint("Error: Too few samples specified by --pca-cluster-names/--pca-clusters.\n");
	    goto plink_ret_1;
	  }
	  if (pca_sample_ct == sample_ct) {
	    logerrprint("Warning: --pca-cluster-names/--pca-clusters has no effect since all samples are\nin the named clusters.\n");
	    pca_sample_exclude = nullptr;
	  } else {
	    LOGPRINTF("--pca-cluster-names/--pca-clusters: %" PRIuPTR " samples specified.\n", pca_sample_ct);
	    ulii = unfiltered_sample_ct - pca_sample_ct;
	  }
	}
	retval = calc_rel(threads, parallel_idx, parallel_tot, calculation_type, relip, bedfile, bed_offset, outname, outname_end, distance_wts_fname, (dist_calc_type & DISTANCE_WTS_NOHEADER), unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ct, marker_ids, max_marker_id_blen, unfiltered_sample_ct, pca_sample_exclude? pca_sample_exclude : sample_exclude, pca_sample_exclude? (&ulii) : (&sample_exclude_ct), sample_ids, max_sample_id_len, set_allele_freqs, &rel_ibc, chrom_info_ptr);
	if (retval) {
	  goto plink_ret_1;
	}
	if ((!pca_sample_exclude) && (sample_ct != unfiltered_sample_ct + sample_exclude_ct)) {
	  sample_ct = unfiltered_sample_ct - sample_exclude_ct;
	  if ((sample_ct < 2) && (distance_req(read_dists_fname, calculation_type) || (calculation_type & (CALC_REGRESS_REL | CALC_PCA | CALC_GENOME | CALC_CLUSTER | CALC_NEIGHBOR)))) {
	    // pathological case
	    sprintf(g_logbuf, "Error: Too many %s pruned for additional pairwise analysis steps.\n", g_species_plural);
	    goto plink_ret_INVALID_CMDLINE_2;
	  }
	}
	if (calculation_type & CALC_REL_CUTOFF) {
	  // ugh, probably better to just stop supporting this
	  bitvec_andnot(sample_exclude, unfiltered_sample_ctl, founder_info);
	  bitvec_andnot(sample_exclude, unfiltered_sample_ctl, sex_nm);
	  bitvec_and(sex_nm, unfiltered_sample_ctl, sex_male);
	  if (pheno_nm_ct) {
	    bitvec_andnot(sample_exclude, unfiltered_sample_ctl, pheno_nm);
	    pheno_nm_ct = popcount_longs(pheno_nm, unfiltered_sample_ctl);
	    if (pheno_c) {
	      bitvec_and(pheno_nm, unfiltered_sample_ctl, pheno_c);
	      pheno_ctrl_ct = pheno_nm_ct - popcount_longs(pheno_c, unfiltered_sample_ctl);
	    }
	  }
	}

	if (calculation_type & CALC_REGRESS_REL) {
	  retval = regress_rel_main(unfiltered_sample_ct, sample_exclude, sample_ct, relip, threads, pheno_d);
	  if (retval) {
	    goto plink_ret_1;
	  }
	}
#ifndef NOLAPACK
	if (calculation_type & CALC_PCA) {
	  retval = calc_pca(bedfile, bed_offset, outname, outname_end, calculation_type, relip, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_allele_ptrs, marker_reverse, unfiltered_sample_ct, sample_exclude, sample_ct, pca_sample_exclude? pca_sample_exclude : sample_exclude, pca_sample_exclude? pca_sample_ct : sample_ct, sample_ids, max_sample_id_len, set_allele_freqs, chrom_info_ptr, rel_ibc);
	} else if (calculation_type & CALC_UNRELATED_HERITABILITY) {
	  if (sample_ct != pheno_nm_ct) {
	    logerrprint("Error: --unrelated-heritability requires phenotype data for all samples.\n(--prune should help.)\n");
	    goto plink_ret_INVALID_CMDLINE;
	  }
	  retval = calc_unrelated_herit(calculation_type, relip, unfiltered_sample_ct, sample_exclude, sample_ct, pheno_d, rel_ibc);
	}
#endif
	bigstack_reset(g_sample_missing_unwt);
	if (retval) {
	  goto plink_ret_1;
	}
	g_sample_missing_unwt = nullptr;
	g_missing_dbl_excluded = nullptr;
      }

      if (calculation_type & CALC_SEXCHECK) {
	retval = sexcheck(bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, sex_nm, sex_male, misc_flags, check_sex_fthresh, check_sex_mthresh, check_sex_f_yobs, check_sex_m_yobs, chrom_info_ptr, set_allele_freqs, &gender_unk_ct);
	if (retval) {
	  goto plink_ret_1;
	}
      }

      if (calculation_type & CALC_TUCC) {
	retval = make_pseudocontrols(bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_cms, marker_pos, marker_allele_ptrs, max_marker_allele_blen, marker_reverse, unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, chrom_info_ptr, fam_ip);
	if (retval) {
	  goto plink_ret_1;
	}
      }
    }

    if (calculation_type & CALC_MAKE_PERM_PHENO) {
      retval = make_perm_pheno(threads, outname, outname_end, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, cluster_ct, cluster_map, cluster_starts, pheno_nm_ct, pheno_nm, pheno_c, pheno_d, output_missing_pheno, permphe_ct);
      if (retval) {
	goto plink_ret_1;
      }
    }

    if ((calculation_type & CALC_GENOME) || genome_skip_write) {
      // er, this probably should be moved inside calc_genome(), since we're
      // using get_trios_and_families() instead of pri elsewhere
      retval = populate_pedigree_rel_info(&pri, unfiltered_sample_ct, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
      if (retval) {
	goto plink_ret_1;
      }
    }
  }

  if (marker_ct) {
    if (calculation_type & CALC_WRITE_SET) {
      retval = write_set(sip, outname, outname_end, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_blen, marker_pos, chrom_info_ptr);
      if (retval) {
	goto plink_ret_1;
      }
    }

    if (calculation_type & CALC_WRITE_SNPLIST) {
      retval = write_snplist(outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, nullptr, 0);
      if (retval) {
	goto plink_ret_1;
      }
    }

    if (calculation_type & CALC_WRITE_VAR_RANGES) {
      retval = write_var_ranges(outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, write_var_range_ct);
      if (retval) {
	goto plink_ret_1;
      }
    }

    if (calculation_type & CALC_LIST_23_INDELS) {
      retval = write_snplist(outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_allele_ptrs, 1);
      if (retval) {
	goto plink_ret_1;
      }
    }

    if (calculation_type & CALC_DUPVAR) {
      if (map_is_unsorted & UNSORTED_BP) {
	logerrprint("Error: --list-duplicate-vars requires a sorted .bim file.  Retry this command\nafter using --make-bed to sort your data.\n");
	goto plink_ret_INVALID_FORMAT;
      }
      retval = list_duplicate_vars(outname, outname_end, dupvar_modifier, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_pos, chrom_info_ptr, marker_allele_ptrs);
      if (retval) {
	goto plink_ret_1;
      }
    }
  }

  if (calculation_type & (CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_MAKE_BIM | CALC_MAKE_FAM | CALC_RECODE)) {
    if (gender_unk_ct && (sex_missing_pheno & MUST_HAVE_SEX)) {
      if (aligned_malloc(unfiltered_sample_ctl * sizeof(intptr_t), &pheno_nm_datagen)) {
	goto plink_ret_NOMEM;
      }
      memcpy(pheno_nm_datagen, pheno_nm, unfiltered_sample_ctl * sizeof(intptr_t));
      bitvec_and(sex_nm, unfiltered_sample_ctl, pheno_nm_datagen);
    }
    if (covar_ct && (calculation_type & (CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_MAKE_FAM | CALC_RECODE)) && sample_ct) {
      retval = write_covars(outname, outname_end, write_covar_modifier, write_covar_dummy_max_categories, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm_datagen? pheno_nm_datagen : pheno_nm, pheno_c, pheno_d, missing_phenod, output_missing_pheno, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (calculation_type & (CALC_MAKE_BED | CALC_MAKE_BIM | CALC_MAKE_FAM)) {
      retval = make_bed(bedfile, bed_offset, bimname, outname, outname_end, calculation_type, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_cms, marker_pos, marker_allele_ptrs, marker_reverse, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info, sex_nm, sex_male, pheno_nm_datagen? pheno_nm_datagen : pheno_nm, pheno_c, pheno_d, output_missing_pheno, map_is_unsorted & (UNSORTED_CHROM | UNSORTED_BP | UNSORTED_SPLIT_CHROM), sample_sort_map, misc_flags, splitx_bound1, splitx_bound2, update_chr, flip_fname, flip_subset_fname, cluster_ptr->zerofname, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len, hh_exists, chrom_info_ptr, fam_ip->mendel_modifier, max_bim_linelen);
      if (retval) {
        goto plink_ret_1;
      }
    }
    if (calculation_type & CALC_RECODE) {
      retval = recode(recode_modifier, bedfile, bed_offset, outname, outname_end, recode_allele_name, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_sample_ct, sample_exclude, sample_ct, marker_ids, max_marker_id_blen, marker_cms, marker_allele_ptrs, max_marker_allele_blen, marker_pos, marker_reverse, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm_datagen? pheno_nm_datagen : pheno_nm, pheno_c, pheno_d, output_missing_pheno, map_is_unsorted & (UNSORTED_CHROM | UNSORTED_BP | UNSORTED_SPLIT_CHROM), misc_flags, hh_exists, chrom_info_ptr);
      if (retval) {
        goto plink_ret_1;
      }
    }
    aligned_free_cond_null(&pheno_nm_datagen);
  }

  if ((!marker_ct) || (!sample_ct)) {
    goto plink_ret_1;
  }

  if ((calculation_type & CALC_EPI) && epi_ip->twolocus_mkr1) {
    retval = twolocus(epi_ip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, plink_maxsnp, marker_allele_ptrs, chrom_info_ptr, unfiltered_sample_ct, sample_exclude, sample_ct, pheno_nm, pheno_nm_ct, pheno_ctrl_ct, pheno_c, sex_male, outname, outname_end, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_SHOW_TAGS) {
    if (map_is_unsorted & UNSORTED_BP) {
      logerrprint("Error: --show-tags requires a sorted .bim file.  Retry this command after using\n--make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = show_tags(ldip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, plink_maxsnp, marker_pos, chrom_info_ptr, unfiltered_sample_ct, founder_info, sex_male, outname, outname_end, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_BLOCKS) {
    if (map_is_unsorted & UNSORTED_BP) {
      logerrprint("Error: --blocks requires a sorted .bim file.  Retry this command after using\n--make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = haploview_blocks(ldip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_blen, marker_pos, chrom_info_ptr, set_allele_freqs, unfiltered_sample_ct, founder_info, pheno_nm, sex_male, outname, outname_end, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_HOMOZYG) {
    if (map_is_unsorted & UNSORTED_BP) {
      logerrprint("Error: Run-of-homozygosity scanning requires a sorted .bim.  Retry this command\nafter using --make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = calc_homozyg(homozyg_ptr, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_blen, plink_maxsnp, marker_allele_ptrs, max_marker_allele_blen, marker_reverse, chrom_info_ptr, marker_pos, sample_ct, unfiltered_sample_ct, sample_exclude, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, outname, outname_end, pheno_nm, pheno_c, pheno_d, output_missing_pheno, sex_male);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_LD_PRUNE) {
    if (map_is_unsorted & UNSORTED_BP) {
      logerrprint("Error: LD-based pruning requires a sorted .bim.  Retry this command after using\n--make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    if (!(ldip->modifier & LD_PRUNE_PAIRPHASE)) {
      retval = ld_prune(ldip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, chrom_info_ptr, set_allele_freqs, marker_pos, unfiltered_sample_ct, founder_info, sex_male, outname, outname_end, hh_exists);
    } else {
      retval = indep_pairphase(ldip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, chrom_info_ptr, set_allele_freqs, marker_pos, unfiltered_sample_ct, founder_info, sex_male, outname, outname_end, hh_exists);
    }
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_FLIPSCAN) {
    if (map_is_unsorted & UNSORTED_BP) {
      logerrprint("Error: LD-based strand flip scanning requires a sorted .bim.  Retry this\ncommand after using --make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = flipscan(ldip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, plink_maxsnp, marker_allele_ptrs, max_marker_allele_blen, chrom_info_ptr, set_allele_freqs, marker_pos, unfiltered_sample_ct, pheno_nm, pheno_c, founder_info, sex_male, outname, outname_end, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if ((calculation_type & CALC_EPI) && epi_ip->ld_mkr1) {
    retval = twolocus(epi_ip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, plink_maxsnp, marker_allele_ptrs, chrom_info_ptr, unfiltered_sample_ct, founder_info, 0, nullptr, 0, 0, nullptr, sex_male, nullptr, nullptr, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_LD) {
    if (!(ldip->modifier & (LD_MATRIX_SHAPEMASK | LD_INTER_CHR))) {
      if (map_is_unsorted & UNSORTED_BP) {
	logerrprint("Error: Windowed --r/--r2 runs require a sorted .bim.  Retry this command after\nusing --make-bed to sort your data.\n");
	goto plink_ret_INVALID_CMDLINE;
      }
      if ((map_is_unsorted & UNSORTED_CM) && (ldip->window_cm != -1)) {
	logerrprint("Error: --ld-window-cm requires nondecreasing CM values on each chromosome.\nRetry this command after regenerating your CM coordinates.\n");
	goto plink_ret_INVALID_CMDLINE;
      }
    }
    retval = ld_report(threads, ldip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, plink_maxsnp, marker_allele_ptrs, max_marker_allele_blen, set_allele_freqs, chrom_info_ptr, marker_pos, marker_cms, unfiltered_sample_ct, founder_info, parallel_idx, parallel_tot, sex_male, outname, outname_end, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }
  if (calculation_type & CALC_TESTMISHAP) {
    if (map_is_unsorted & UNSORTED_BP) {
      logerrprint("Error: --test-mishap requires a sorted .bim.  Retry this command after using\n--make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = test_mishap(bedfile, bed_offset, outname, outname_end, output_min_p, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_allele_ptrs, min_maf, chrom_info_ptr, unfiltered_sample_ct, sample_exclude, sample_ct);
    if (retval) {
      goto plink_ret_1;
    }
  }

  // sometimes no more need for marker_ids/marker_allele_ptrs, conditional
  // unload to clear space for IBS matrix, etc.?  (probably want to initially
  // load at far end of stack to make this workable...)

  if (calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) {
    bigstack_mark_postcluster = g_bigstack_base;
    ulii = (sample_ct * (sample_ct - 1)) >> 1;
    if (cluster_ptr->mds_dim_ct) {
#ifndef __LP64__
      // catch 32-bit intptr_t overflow
      if (sample_ct > 23169) {
        goto plink_ret_NOMEM;
      }
#endif
      if (((!read_dists_fname) && (!read_genome_fname)) || (cluster_ptr->modifier & CLUSTER_MISSING)) {
	if ((!(cluster_ptr->modifier & CLUSTER_MDS)) || (!cluster_ct)) {
          if (bigstack_alloc_d(ulii, &mds_plot_dmatrix_copy)) {
            goto plink_ret_NOMEM;
          }
	} else {
	  ulii = cluster_ct + sample_ct - cluster_starts[cluster_ct];
          if (bigstack_alloc_d((ulii * (ulii - 1)) / 2, &mds_plot_dmatrix_copy)) {
            goto plink_ret_NOMEM;
          }
	}
      }
    }

    if (cluster_ct) {
      ulii = cluster_ct + sample_ct - cluster_starts[cluster_ct];
#ifndef __LP64__
      if (ulii > 23169) {
	goto plink_ret_NOMEM;
      }
#endif
      ulii = (ulii * (ulii - 1)) >> 1;
#ifndef __LP64__
    } else if (sample_ct > 23169) {
      goto plink_ret_NOMEM;
#endif
    }
    if (bigstack_alloc_ul(BITCT_TO_WORDCT(ulii), &cluster_merge_prevented)) {
      goto plink_ret_NOMEM;
    }
    if (cluster_ct || (calculation_type & CALC_GENOME) || genome_skip_write) {
      if (bigstack_alloc_d(ulii, &cluster_sorted_ibs)) {
	goto plink_ret_NOMEM;
      }
      if (cluster_ptr->modifier & CLUSTER_GROUP_AVG) {
        fill_double_zero(ulii, cluster_sorted_ibs);
      } else {
	for (uljj = 0; uljj < ulii; uljj++) {
	  cluster_sorted_ibs[uljj] = 1.0;
	}
      }
    }
    bigstack_mark_precluster = g_bigstack_base;
  }

  bigstack_mark2 = g_bigstack_base;

  if (distance_req(read_dists_fname, calculation_type)) {
    retval = calc_distance(threads, parallel_idx, parallel_tot, bedfile, bed_offset, outname, outname_end, read_dists_fname, distance_wts_fname, distance_exp, calculation_type, dist_calc_type, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, set_allele_freqs, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, chrom_info_ptr);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (read_dists_fname && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))) {
    // use delayed and specialized load for --cluster/--neighbour, since PPC
    // values may be needed, and user may want to process a distance matrix too
    // large to be loaded in memory by doing some pre-clustering
    if (bigstack_alloc_d((sample_ct * (sample_ct - 1)) / 2, &g_dists)) {
      goto plink_ret_NOMEM;
    }
    retval = read_dists(read_dists_fname, read_dists_id_fname, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, 0, nullptr, nullptr, 0, 0, g_dists, 0, nullptr, nullptr);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_IBS_TEST) {
    if (cluster_starts) {
      logerrprint("Error: --ibs-test does not currently support within-cluster permutation.\nContact the developers if you need this.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto plink_ret_1;
    }
    retval = ibs_test_calc(threads, read_dists_fname, unfiltered_sample_ct, sample_exclude, sample_ct, ibs_test_perms, pheno_nm_ct, pheno_ctrl_ct, pheno_nm, pheno_c);
    if (retval) {
      goto plink_ret_1;
    }
  }
  if (calculation_type & CALC_GROUPDIST) {
    retval = groupdist_calc(threads, unfiltered_sample_ct, sample_exclude, sample_ct, groupdist_iters, groupdist_d, pheno_nm_ct, pheno_ctrl_ct, pheno_nm, pheno_c);
    if (retval) {
      goto plink_ret_1;
    }
  }
  if (calculation_type & CALC_REGRESS_DISTANCE) {
    if (sample_ct != pheno_nm_ct) {
      logerrprint("Error: --regress-distance requires phenotype data for all samples.  (--prune\nshould help.)\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = regress_distance(threads, calculation_type, pheno_d, unfiltered_sample_ct, sample_exclude, sample_ct, g_thread_ct, regress_iters, regress_d);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (read_dists_fname && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))) {
    bigstack_reset(g_dists);
    g_dists = nullptr;
  }

  if ((calculation_type & CALC_GENOME) || genome_skip_write) {
    bigstack_reset(bigstack_mark2);
    g_dists = nullptr;
    retval = calc_genome(threads, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, marker_pos, set_allele_freqs, nchrobs, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info, parallel_idx, parallel_tot, outname, outname_end, nonfounders, calculation_type, genome_modifier, ppc_gap, genome_min_pi_hat, genome_max_pi_hat, pheno_nm, pheno_c, pri, genome_skip_write);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_HET) {
    retval = het_report(bedfile, bed_offset, outname, outname_end, (misc_flags / MISC_HET_GZ) & 1, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, (misc_flags & MISC_HET_SMALL_SAMPLE)? founder_info : nullptr, chrom_info_ptr, set_allele_freqs);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_FST) {
    retval = fst_report(bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_pos, chrom_info_ptr, unfiltered_sample_ct, sample_exclude, pheno_nm, (misc_flags & MISC_FST_CC)? pheno_c : nullptr, cluster_ct, cluster_map, cluster_starts);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) {
    retval = calc_cluster_neighbor(threads, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, set_allele_freqs, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, read_dists_fname, read_dists_id_fname, read_genome_fname, outname, outname_end, calculation_type, cluster_ct, cluster_map, cluster_starts, cluster_ptr, missing_pheno, neighbor_n1, neighbor_n2, ppc_gap, pheno_c, mds_plot_dmatrix_copy, cluster_merge_prevented, cluster_sorted_ibs, bigstack_mark_precluster, bigstack_mark_postcluster);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if ((calculation_type & CALC_EPI) && (epi_ip->modifier & (EPI_FAST | EPI_REG))) {
    if ((map_is_unsorted & UNSORTED_BP) && (epi_ip->modifier & EPI_FAST_CASE_ONLY)) {
      logerrprint("Error: --fast-epistasis case-only requires a sorted .bim.  Retry this command\nafter using --make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = epistasis_report(threads, epi_ip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, marker_pos, plink_maxsnp, chrom_info_ptr, unfiltered_sample_ct, pheno_nm, pheno_nm_ct, pheno_ctrl_ct, pheno_c, pheno_d, parallel_idx, parallel_tot, outname, outname_end, output_min_p, glm_vif_thresh, sip);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_SCORE) {
    retval = score_report(sc_ip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, marker_allele_ptrs, set_allele_freqs, sample_ct, unfiltered_sample_ct, sample_exclude, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_pheno, hh_exists, chrom_info_ptr, outname, outname_end);
    if (retval) {
      goto plink_ret_1;
    }
  }

#if defined __cplusplus && !defined _WIN32
  if (calculation_type & CALC_RPLUGIN) {
    retval = rserve_call(rplugin_fname, rplugin_host_or_socket, rplugin_port, (misc_flags / MISC_RPLUGIN_DEBUG) & 1, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_blen, marker_allele_ptrs, marker_pos, plink_maxsnp, chrom_info_ptr, unfiltered_sample_ct, pheno_nm, pheno_nm_ct, pheno_c, pheno_d, cluster_ct, cluster_map, cluster_starts, covar_ct, covar_d, outname, outname_end);
    if (retval) {
      goto plink_ret_1;
    }
  }
#endif

  if (calculation_type & CALC_LASSO_LAMBDA) {
    retval = lasso_lambda(marker_exclude, marker_reverse, chrom_info_ptr, sample_exclude, sex_male, pheno_nm, covar_nm, bed_offset, unfiltered_marker_ct, marker_ct, unfiltered_sample_ct, sample_ct, pheno_nm_ct, hh_exists, lasso_lambda_iters, lasso_h2, bedfile, outname, outname_end, &lasso_minlambda);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & (CALC_MODEL | CALC_GXE | CALC_GLM | CALC_LASSO | CALC_CMH | CALC_HOMOG | CALC_TESTMISS | CALC_TDT | CALC_DFAM | CALC_QFAM)) {
    // can't use pheno_ctrl_ct in here since new phenotypes may be loaded, and
    // we don't bother updating it...
    if ((!pheno_all) && (!loop_assoc_fname)) {
      outname_end2 = outname_end;
      goto plink_skip_all_pheno;
    }
    uii = 0; // phenotype/cluster number
    *outname_end = '.';
    if (loop_assoc_fname) {
      retval = load_clusters(loop_assoc_fname, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, sample_ids, max_sample_id_len, mwithin_col, (misc_flags / MISC_LOAD_CLUSTER_KEEP_NA) & 1, &cluster_ct, &cluster_map, &cluster_starts, &cluster_ids, &max_cluster_id_len, nullptr, nullptr, nullptr, nullptr, 0);
      if (retval) {
	goto plink_ret_1;
      }
      // bugfix (17 Jun 2019): need to reset pheno_nm, otherwise missing values
      // in original phenotype remain missing in the --loop-assoc analysis
      bitarr_invert_copy(sample_exclude, unfiltered_sample_ct, pheno_nm);
      pheno_nm_ct = unfiltered_sample_ct - sample_exclude_ct;

      if (pheno_d) {
	free(pheno_d);
	pheno_d = nullptr;
      }
      if (!pheno_c) {
	if (aligned_malloc(unfiltered_sample_ctl * sizeof(intptr_t), &pheno_c)) {
	  goto plink_ret_NOMEM;
	}
      }
    } else {
      bigstack_mark = g_bigstack_base;
      retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, sample_exclude_ct, sample_ids, max_sample_id_len, 0, 0, strcmp_deref, &cptr, &uiptr);
      if (retval) {
	goto plink_ret_1;
      }
    }
    do {
      if (loop_assoc_fname) {
	if (uii == cluster_ct) {
	  break;
	}
	outname_end2 = strcpya(&(outname_end[1]), &(cluster_ids[uii * max_cluster_id_len]));
	fill_ulong_zero(unfiltered_sample_ctl, pheno_c);
	ukk = cluster_starts[uii + 1];
	for (ujj = cluster_starts[uii]; ujj < ukk; ujj++) {
	  SET_BIT(cluster_map[ujj], pheno_c);
	}
	uii++;
      } else {
	// --all-pheno
	if (pheno_modifier & PHENO_MERGE) {
	  memcpy(pheno_nm, orig_pheno_nm, unfiltered_sample_ctl * sizeof(intptr_t));
	  if (orig_pheno_c) {
	    if (!pheno_c) {
	      free(pheno_d);
	      pheno_d = nullptr;
	      if (aligned_malloc(unfiltered_sample_ctl * sizeof(intptr_t), &pheno_c)) {
		goto plink_ret_NOMEM;
	      }
	    }
	    memcpy(pheno_c, orig_pheno_c, unfiltered_sample_ctl * sizeof(intptr_t));
	  } else {
	    memcpy(pheno_d, orig_pheno_d, unfiltered_sample_ct * sizeof(double));
	  }
	} else {
	  fill_ulong_zero(unfiltered_sample_ctl, pheno_nm);
	  aligned_free_cond_null(&pheno_c);
	  if (pheno_d) {
	    free(pheno_d);
	    pheno_d = nullptr;
	  }
	}
        // bugfix (9 Jan 2018): forgot to increment phenotype index when
        // skipping empty pheno
      plink_skip_empty_pheno:
	uii++;
	rewind(phenofile);
	outname_end[1] = '\0';
	retval = load_pheno(phenofile, unfiltered_sample_ct, sample_exclude_ct, cptr, max_sample_id_len, uiptr, missing_pheno, (misc_flags / MISC_AFFECTION_01) & 1, uii, nullptr, pheno_nm, &pheno_c, &pheno_d, &(outname_end[1]), (uintptr_t)((&(outname[FNAMESIZE - 32])) - outname_end));
	if (retval == LOAD_PHENO_LAST_COL) {
	  bigstack_reset(bigstack_mark);
	  retval = 0; // exit code bugfix
	  break;
	} else if (retval) {
	  goto plink_ret_1;
	}
	bitvec_andnot(sample_exclude, unfiltered_sample_ctl, pheno_nm);
	if (gender_unk_ct && (!(sex_missing_pheno & ALLOW_NO_SEX))) {
	  bitvec_and(sex_nm, unfiltered_sample_ctl, pheno_nm);
	}
	pheno_nm_ct = popcount_longs(pheno_nm, unfiltered_sample_ctl);
	if (!pheno_nm_ct) {
	  goto plink_skip_empty_pheno;
	}
	if (!outname_end[1]) {
	  outname_end[1] = 'P';
	  outname_end2 = uint32toa(uii, &(outname_end[2]));
	} else {
          outname_end2 = (char*)memchr(&(outname_end[1]), '\0', FNAMESIZE);
	}
      }
      *outname_end2 = '\0';
      if (pheno_c) {
	bitvec_and(pheno_nm, unfiltered_sample_ctl, pheno_c);
        ujj = popcount_longs(pheno_c, unfiltered_sample_ctl);
	ukk = pheno_nm_ct - ujj;
	ulii = unfiltered_sample_ct - sample_exclude_ct - pheno_nm_ct;
        if (ulii) {
          LOGPRINTF("%s has %u case%s, %u control%s, and %" PRIuPTR " missing phenotype%s.\n", &(outname_end[1]), ujj, (ujj == 1)? "" : "s", ukk, (ukk == 1)? "" : "s", ulii, (ulii == 1)? "" : "s");
	} else {
          LOGPRINTF("%s has %u case%s and %u control%s.\n", &(outname_end[1]), ujj, (ujj == 1)? "" : "s", ukk, (ukk == 1)? "" : "s");
	}
      }
    plink_skip_all_pheno:
      if (calculation_type & CALC_MODEL) {
	if (pheno_d) {
	  if (model_modifier & MODEL_ASSOC) {
	    retval = qassoc(threads, bedfile, bed_offset, outname, outname_end2, model_modifier, model_mperm_val, pfilter, output_min_p, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_pos, marker_allele_ptrs, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, cluster_ct, cluster_map, cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_d, founder_info, sex_male, hh_exists, ldip->modifier & LD_IGNORE_X, perm_batch_size, sip);
	  }
	} else {
	  retval = model_assoc(threads, bedfile, bed_offset, outname, outname_end2, model_modifier, model_cell_ct, model_mperm_val, ci_size, ci_zt, pfilter, output_min_p, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_pos, marker_allele_ptrs, max_marker_allele_blen, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, cluster_ct, cluster_map, loop_assoc_fname? nullptr : cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, founder_info, sex_male, hh_exists, ldip->modifier & LD_IGNORE_X, perm_batch_size, sip);
	}
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if (calculation_type & CALC_GLM) {
	if (!(glm_modifier & GLM_NO_SNP)) {
	  if (pheno_d) {
#ifndef NOLAPACK
	    retval = glm_linear_assoc(threads, bedfile, bed_offset, outname, outname_end2, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, parameters_range_list_ptr, tests_range_list_ptr, ci_size, ci_zt, pfilter, output_min_p, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_pos, marker_allele_ptrs, max_marker_allele_blen, marker_reverse, condition_mname, condition_fname, chrom_info_ptr, unfiltered_sample_ct, sample_ct, sample_exclude, cluster_ct, cluster_map, cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_d, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, founder_info, sex_nm, sex_male, ldip->modifier & LD_IGNORE_X, hh_exists, perm_batch_size, sip);
#else
            logerrprint("Warning: Skipping --logistic on --all-pheno QT since this is a no-LAPACK " PROG_NAME_CAPS"\nbuild.\n");
#endif
	  } else {
	    retval = glm_logistic_assoc(threads, bedfile, bed_offset, outname, outname_end2, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, parameters_range_list_ptr, tests_range_list_ptr, ci_size, ci_zt, pfilter, output_min_p, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_pos, marker_allele_ptrs, max_marker_allele_blen, marker_reverse, condition_mname, condition_fname, chrom_info_ptr, unfiltered_sample_ct, sample_ct, sample_exclude, cluster_ct, cluster_map, cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, founder_info, sex_nm, sex_male, ldip->modifier & LD_IGNORE_X, hh_exists, perm_batch_size, sip);
	  }
	} else {
	  if (pheno_d) {
#ifndef NOLAPACK
	    retval = glm_linear_nosnp(threads, bedfile, bed_offset, outname, outname_end2, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, parameters_range_list_ptr, tests_range_list_ptr, ci_size, ci_zt, pfilter, output_min_p, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_reverse, condition_mname, condition_fname, chrom_info_ptr, unfiltered_sample_ct, sample_ct, sample_exclude, cluster_ct, cluster_map, cluster_starts, mperm_save, pheno_nm_ct, pheno_nm, pheno_d, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_nm, sex_male, hh_exists, perm_batch_size, sip);
#else
            logerrprint("Warning: Skipping --logistic on --all-pheno QT since this is a no-LAPACK " PROG_NAME_CAPS"\nbuild.\n");
#endif

	  } else {
	    retval = glm_logistic_nosnp(threads, bedfile, bed_offset, outname, outname_end2, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, parameters_range_list_ptr, tests_range_list_ptr, ci_size, ci_zt, pfilter, output_min_p, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_reverse, condition_mname, condition_fname, chrom_info_ptr, unfiltered_sample_ct, sample_ct, sample_exclude, cluster_ct, cluster_map, cluster_starts, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_nm, sex_male, hh_exists, perm_batch_size, sip);
	  }
	}
	if (retval) {
	  goto plink_ret_1;
	}
      }
      // if case/control phenotype loaded with --all-pheno, skip --gxe
      if ((calculation_type & CALC_GXE) && pheno_d) {
	retval = gxe_assoc(bedfile, bed_offset, outname, outname_end2, output_min_p, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, sample_ct, sample_exclude, pheno_nm, pheno_d, gxe_covar_nm, gxe_covar_c, sex_male, hh_exists);
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if (calculation_type & CALC_LASSO) {
	retval = lasso(threads, bedfile, bed_offset, outname, outname_end2, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_allele_ptrs, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, sample_ct, pheno_nm_ct, lasso_h2, lasso_minlambda, lasso_select_covars_range_list_ptr, misc_flags, sample_exclude, pheno_nm, pheno_c, pheno_d, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_male, hh_exists);
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if ((calculation_type & CALC_CMH) && pheno_c) {
	if (!(cluster_ptr->modifier & CLUSTER_CMH2)) {
          retval = cmh_assoc(threads, bedfile, bed_offset, outname, outname_end2, cluster_ptr->cmh_mperm_val, cluster_ptr->modifier, ci_size, pfilter, output_min_p, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_pos, marker_allele_ptrs, marker_reverse, chrom_info_ptr, set_allele_freqs, unfiltered_sample_ct, cluster_ct, cluster_map, cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, sex_male, hh_exists, sip);
	} else {
          retval = cmh2_assoc(bedfile, bed_offset, outname, outname_end, output_min_p, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, cluster_ct, cluster_map, cluster_starts, pheno_nm_ct, pheno_nm, pheno_c, sex_male, hh_exists);
	}
        if (retval) {
          goto plink_ret_1;
	}
      }
      if ((calculation_type & CALC_HOMOG) && pheno_c) {
	retval = homog_assoc(bedfile, bed_offset, outname, outname_end, output_min_p, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_allele_ptrs, max_marker_allele_blen, marker_reverse, chrom_info_ptr, set_allele_freqs, unfiltered_sample_ct, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len, pheno_nm_ct, pheno_nm, pheno_c, sex_male, hh_exists);
        if (retval) {
          goto plink_ret_1;
	}
      }
      if ((calculation_type & CALC_TESTMISS) && pheno_c) {
        retval = testmiss(threads, bedfile, bed_offset, outname, outname_end2, testmiss_mperm_val, testmiss_modifier, pfilter, output_min_p, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, chrom_info_ptr, unfiltered_sample_ct, cluster_ct, cluster_map, loop_assoc_fname? nullptr : cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, sex_male, hh_exists);
        if (retval) {
	  goto plink_ret_1;
	}
      }
      if ((calculation_type & CALC_TDT) && pheno_c) {
	retval = tdt(threads, bedfile, bed_offset, outname, outname_end2, ci_size, ci_zt, pfilter, output_min_p, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_pos, marker_allele_ptrs, max_marker_allele_blen, marker_reverse, unfiltered_sample_ct, sample_exclude, sample_ct, apip, mperm_save, pheno_nm, pheno_c, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, chrom_info_ptr, hh_exists, fam_ip);
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if ((calculation_type & CALC_DFAM) && pheno_c) {
	retval = dfam(threads, bedfile, bed_offset, outname, outname_end2, pfilter, output_min_p, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_allele_ptrs, max_marker_allele_blen, marker_reverse, unfiltered_sample_ct, sample_exclude, sample_ct, cluster_ct, cluster_map, loop_assoc_fname? nullptr : cluster_starts, apip, mperm_save, pheno_c, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, chrom_info_ptr, hh_exists, (cluster_ptr->fname != nullptr), perm_batch_size, fam_ip, sip);
	if (retval) {
	  goto plink_ret_1;
	}
      }
#ifndef NOLAPACK
      if ((calculation_type & CALC_QFAM) && pheno_d) {
	if (covar_ct) {
	  logerrprint("Warning: The QFAM test ignores covariates.\n");
	}
	if (cluster_ct) {
	  logerrprint("Warning: Clusters have no effect on the QFAM permutation test.\n");
	}
	if (mtest_adjust && (fam_ip->qfam_modifier & QFAM_PERM)) {
	  logerrprint("Warning: The QFAM test does not support --adjust.  Use max(T) permutation to\nobtain multiple-testing corrected p-values.\n");
	}
        retval = qfam(threads, bedfile, bed_offset, outname, outname_end2, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_pos, marker_allele_ptrs, marker_reverse, unfiltered_sample_ct, sample_exclude, sample_ct, apip, pheno_nm, pheno_d, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, chrom_info_ptr, hh_exists, perm_batch_size, fam_ip);
        if (retval) {
	  goto plink_ret_1;
	}
      }
#endif
    } while (pheno_all || loop_assoc_fname);
  }
  if (calculation_type & CALC_CLUMP) {
    if (map_is_unsorted & UNSORTED_BP) {
      logerrprint("Error: --clump requires a sorted .bim.  Retry this command after using\n--make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = clump_reports(bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, plink_maxsnp, marker_pos, marker_allele_ptrs, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, founder_info, clump_ip, sex_male, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }
  while (0) {
  plink_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  plink_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  plink_ret_INVALID_FORMAT_2:
    logerrprintb();
  plink_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  plink_ret_INVALID_CMDLINE_2:
    logerrprintb();
  plink_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  plink_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  plink_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  plink_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 plink_ret_1:
  aligned_free_cond(pheno_nm_datagen);
  free_cond(orig_pheno_d);
  aligned_free_cond(orig_pheno_c);
  aligned_free_cond(orig_pheno_nm);
  free_cond(pheno_d);
  aligned_free_cond(pheno_c);
  fclose_cond(phenofile);
  fclose_cond(bedfile);
  cleanup_allele_storage(max_marker_allele_blen - 1, 2 * unfiltered_marker_ct, marker_allele_ptrs);
  return retval;
}

// meta-analysis-report-dups + terminating null, or 'recode 01 fastphase-1chr'
#define MAX_FLAG_LEN 26

static inline int32_t is_flag(char* param) {
  unsigned char ucc = param[1];
  return ((*param == '-') && ((ucc > '9') || ((ucc < '0') && (ucc != '.') && (ucc != '\0'))));
}

static inline char* is_flag_start(char* param) {
  unsigned char ucc = param[1];
  if ((*param == '-') && ((ucc > '9') || ((ucc < '0') && (ucc != '.') && (ucc != '\0')))) {
    return (ucc == '-')? (&(param[2])) : (&(param[1]));
  }
  return nullptr;
}

uint32_t param_count(int32_t argc, char** argv, int32_t flag_idx) {
  // Counts the number of optional parameters given to the flag at position
  // flag_idx, treating any parameter not beginning with "--" as optional.
  int32_t opt_params = 0;
  int32_t cur_idx = flag_idx + 1;
  while (cur_idx < argc) {
    if (is_flag(argv[cur_idx])) {
      break;
    }
    opt_params++;
    cur_idx++;
  }
  return opt_params;
}

int32_t enforce_param_ct_range(uint32_t param_ct, char* flag_name, uint32_t min_ct, uint32_t max_ct) {
  if (param_ct > max_ct) {
    if (max_ct > min_ct) {
      sprintf(g_logbuf, "Error: %s accepts at most %u parameter%s.\n", flag_name, max_ct, (max_ct == 1)? "" : "s");
    } else {
      sprintf(g_logbuf, "Error: %s only accepts %u parameter%s.\n", flag_name, max_ct, (max_ct == 1)? "" : "s");
    }
    return -1;
  } else if (param_ct < min_ct) {
    if (min_ct == 1) {
      sprintf(g_logbuf, "Error: Missing %s parameter.\n", flag_name);
    } else {
      sprintf(g_logbuf, "Error: %s requires %s%u parameters.\n", flag_name, (min_ct < max_ct)? "at least " : "", min_ct);
    }
    return -1;
  }
  return 0;
}

int32_t parse_next_range(uint32_t param_ct, char range_delim, char** argv, uint32_t* cur_param_idx_ptr, char** cur_arg_pptr, char** range_start_ptr, uint32_t* rs_len_ptr, char** range_end_ptr, uint32_t* re_len_ptr) {
  // Starts reading from argv[cur_param_idx][cur_pos].  If a valid range is
  // next, range_start + rs_len + range_end + re_len are updated.  If only a
  // single item is next, range_end is set to nullptr and range_start + rs_len
  // are updated.  If there are no items left, range_start is set to nullptr.
  // If the input is not well-formed, -1 is returned instead of 0.
  uint32_t cur_param_idx = *cur_param_idx_ptr;
  char* cur_arg_ptr = *cur_arg_pptr;
  char cc;
  if (cur_param_idx > param_ct) {
    *cur_arg_pptr = nullptr;
    return 0;
  }
  while (1) {
    cc = *cur_arg_ptr;
    if (!cc) {
      *cur_param_idx_ptr = ++cur_param_idx;
      if (cur_param_idx > param_ct) {
	*range_start_ptr = nullptr;
	return 0;
      }
      cur_arg_ptr = argv[cur_param_idx];
      cc = *cur_arg_ptr;
    }
    if (cc == range_delim) {
      return -1;
    }
    if (cc != ',') {
      break;
    }
    cur_arg_ptr++;
  }
  *range_start_ptr = cur_arg_ptr;
  do {
    cc = *(++cur_arg_ptr);
    if ((!cc) || (cc == ',')) {
      *rs_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_start_ptr));
      *cur_arg_pptr = cur_arg_ptr;
      *range_end_ptr = nullptr;
      return 0;
    }
  } while (cc != range_delim);
  *rs_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_start_ptr));
  cc = *(++cur_arg_ptr);
  if ((!cc) || (cc == ',') || (cc == range_delim)) {
    return -1;
  }
  *range_end_ptr = cur_arg_ptr;
  do {
    cc = *(++cur_arg_ptr);
    if (cc == range_delim) {
      return -1;
    }
  } while (cc && (cc != ','));
  *re_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_end_ptr));
  *cur_arg_pptr = cur_arg_ptr;
  return 0;
}

int32_t parse_chrom_ranges(uint32_t param_ct, char range_delim, char** argv, uintptr_t* chrom_mask, Chrom_info* chrom_info_ptr, uint32_t allow_extra_chroms, char* cur_flag_str) {
  uint32_t argct = 0;
  uint32_t cur_param_idx = 1;
  int32_t retval = 0;
  char* cur_arg_ptr;
  char* range_start;
  uint32_t rs_len;
  char* range_end;
  uint32_t re_len;
  int32_t chrom_code_start;
  int32_t chrom_code_end;
  if (param_ct) {
    cur_arg_ptr = argv[1];
    while (1) {
      if (parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	sprintf(g_logbuf, "Error: Invalid --%s parameter '%s'.\n", cur_flag_str, argv[cur_param_idx]);
	goto parse_chrom_ranges_ret_INVALID_CMDLINE_WWA;
      }
      if (!range_start) {
	break;
      }
      chrom_code_start = get_chrom_code_counted(chrom_info_ptr, rs_len, range_start);
      if (chrom_code_start < 0) {
	range_start[rs_len] = '\0';
	if (!allow_extra_chroms) {
	  sprintf(g_logbuf, "Error: Invalid --%s chromosome code '%s'.\n", cur_flag_str, range_start);
	  goto parse_chrom_ranges_ret_INVALID_CMDLINE_WWA;
	} else if (range_end) {
	  goto parse_chrom_ranges_ret_INVALID_CMDLINE_NONSTD;
	}
        if (push_ll_str(range_start, &(chrom_info_ptr->incl_excl_name_stack))) {
	  goto parse_chrom_ranges_ret_NOMEM;
	}
      } else if (range_end) {
        chrom_code_end = get_chrom_code_counted(chrom_info_ptr, re_len, range_end);
	if (chrom_code_end < 0) {
	  if (!allow_extra_chroms) {
	    range_end[re_len] = '\0';
	    sprintf(g_logbuf, "Error: Invalid --%s chromosome code '%s'.\n", cur_flag_str, range_end);
	    goto parse_chrom_ranges_ret_INVALID_CMDLINE_WWA;
	  } else {
	    goto parse_chrom_ranges_ret_INVALID_CMDLINE_NONSTD;
	  }
	}
        if (chrom_code_end <= chrom_code_start) {
	  range_start[rs_len] = '\0';
	  range_end[re_len] = '\0';
	  sprintf(g_logbuf, "Error: --%s chromosome code '%s' is not greater than '%s'.\n", cur_flag_str, range_end, range_start);
	  goto parse_chrom_ranges_ret_INVALID_CMDLINE_WWA;
	}
	fill_bits(chrom_code_start, chrom_code_end + 1 - chrom_code_start, chrom_mask);
      } else {
        set_bit(chrom_code_start, chrom_mask);
      }
      argct++;
    }
  }
  if (!argct) {
    LOGERRPRINTF("Error: --%s requires at least one value.\n%s", cur_flag_str, errstr_append);
    return RET_INVALID_CMDLINE;
  }
  while (0) {
  parse_chrom_ranges_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  parse_chrom_ranges_ret_INVALID_CMDLINE_NONSTD:
    logerrprint("Error: Chromosome ranges cannot include nonstandard names.\n");
    retval = RET_INVALID_CMDLINE;
    break;
  parse_chrom_ranges_ret_INVALID_CMDLINE_WWA:
    wordwrapb(0);
    logerrprintb();
    logerrprint(errstr_append);
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

int32_t parse_name_ranges(uint32_t param_ct, char range_delim, char** argv, Range_list* range_list_ptr, uint32_t require_posint) {
  uint32_t name_ct = 0;
  uint32_t cur_param_idx = 1;
  uint32_t name_max_len = 0;
  char* cur_arg_ptr;
  char* range_start;
  uint32_t rs_len;
  char* range_end;
  uint32_t re_len;
  char* cur_name_str;
  char* dup_check;
  unsigned char* cur_name_starts_range;
  uint32_t last_val;
  uint32_t cur_val;
  // two passes.  first pass: count parameters, determine name_max_len;
  // then allocate memory; then fill it.
  if (param_ct) {
    cur_arg_ptr = argv[1];
    while (1) {
      if (parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	LOGERRPRINTFWW("Error: Invalid %s parameter '%s'.\n", argv[0], argv[cur_param_idx]);
	logerrprint(errstr_append);
        return RET_INVALID_CMDLINE;
      }
      if (!range_start) {
	break;
      }
      name_ct++;
      if (rs_len > name_max_len) {
	name_max_len = rs_len; // does NOT include trailing null yet
      }
      if (range_end) {
	name_ct++;
	if (re_len > name_max_len) {
	  name_max_len = re_len;
	}
      }
    }
  }
  if (!name_ct) {
    LOGERRPRINTF("Error: %s requires at least one value.\n%s", argv[0], errstr_append);
    return RET_INVALID_CMDLINE;
  }
  range_list_ptr->name_max_len = ++name_max_len;
  range_list_ptr->name_ct = name_ct;
  range_list_ptr->names = (char*)malloc(name_ct * ((uintptr_t)name_max_len));
  if (!range_list_ptr->names) {
    return RET_NOMEM;
  }
  range_list_ptr->starts_range = (unsigned char*)malloc(name_ct * sizeof(char));
  if (!range_list_ptr->starts_range) {
    return RET_NOMEM;
  }
  cur_name_str = range_list_ptr->names;
  cur_name_starts_range = range_list_ptr->starts_range;
  cur_param_idx = 1;
  cur_arg_ptr = argv[1];
  while (1) {
    parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len);
    if (!range_start) {
      if (require_posint) {
	last_val = 0;
	for (cur_param_idx = 0; cur_param_idx < name_ct; cur_param_idx++) {
	  cur_name_str = &(range_list_ptr->names[cur_param_idx * ((uintptr_t)name_max_len)]);
	  dup_check = cur_name_str; // actually a numeric check
	  do {
	    if (is_not_digit(*dup_check)) {
	      LOGERRPRINTFWW("Error: Invalid %s parameter '%s'.\n", argv[0], cur_name_str);
	      return RET_INVALID_CMDLINE;
	    }
	  } while (*(++dup_check));
	  if (scan_posint_defcapx(cur_name_str, &cur_val)) {
	    LOGERRPRINTFWW("Error: Invalid %s parameter '%s'.\n", argv[0], cur_name_str);
	    return RET_INVALID_CMDLINE;
	  }
	  if (range_list_ptr->starts_range[cur_param_idx]) {
	    last_val = cur_val;
	  } else {
	    if (cur_val <= last_val) {
	      LOGERRPRINTFWW("Error: Invalid %s range '%s-%s'.\n", argv[0], &(range_list_ptr->names[(cur_param_idx - 1) * name_max_len]), cur_name_str);
	      return RET_INVALID_CMDLINE;
	    }
	    last_val = 0;
	  }
	}
      }
      return 0;
    }
    memcpyx(cur_name_str, range_start, rs_len, 0);
    dup_check = range_list_ptr->names;
    while (dup_check < cur_name_str) {
      if (!memcmp(dup_check, cur_name_str, rs_len + 1)) {
	LOGERRPRINTFWW("Error: Duplicate %s parameter '%s'.\n", argv[0], cur_name_str);
	return RET_INVALID_CMDLINE;
      }
      dup_check = &(dup_check[name_max_len]);
    }
    cur_name_str = &(cur_name_str[name_max_len]);
    if (range_end) {
      *cur_name_starts_range++ = 1;
      memcpyx(cur_name_str, range_end, re_len, 0);
      dup_check = range_list_ptr->names;
      while (dup_check < cur_name_str) {
	if (!memcmp(dup_check, cur_name_str, rs_len + 1)) {
	  LOGERRPRINTFWW("Error: Duplicate %s parameter '%s'.\n", argv[0], cur_name_str);
	  return RET_INVALID_CMDLINE;
	}
        dup_check = &(dup_check[name_max_len]);
      }
      cur_name_str = &(cur_name_str[name_max_len]);
      *cur_name_starts_range++ = 0;
    } else {
      *cur_name_starts_range++ = 0;
    }
  }
}

void invalid_arg(char* argv) {
  LOGPREPRINTFWW("Error: Unrecognized flag ('%s').\n", argv);
}

void print_ver() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}

int32_t rerun(uint32_t rerun_argv_pos, uint32_t rerun_parameter_present, int32_t* argc_ptr, uint32_t* cur_arg_ptr, char*** argv_ptr, char*** subst_argv_ptr, char** rerun_buf_ptr) {
  // caller is responsible for freeing rerun_buf
  char** argv = *argv_ptr;
  FILE* rerunfile = fopen(rerun_parameter_present? argv[rerun_argv_pos + 1] : (PROG_NAME_STR ".log"), "r");
  uintptr_t line_idx = 1;
  char** subst_argv2 = nullptr;
  uint32_t argc = (uint32_t)(*argc_ptr);
  uint32_t cur_arg = *cur_arg_ptr;
  int32_t retval = 0;
  char* rerun_buf;
  char* rerun_start_ptr;
  char* sptr;
  char* argptr;
  char* argptr2;
  char* load_ptr;
  uint32_t line_byte_ct;
  uint32_t loaded_arg_ct;
  uint32_t duplicate_ct;
  uint32_t loaded_arg_idx;
  uint32_t cmdline_arg_idx;
  uint32_t new_arg_idx;
  uint32_t slen;
  uint32_t slen2;
  if (!rerunfile) {
    goto rerun_ret_OPEN_FAIL;
  }
  g_textbuf[MAXLINELEN - 1] = ' ';
  if (!fgets(g_textbuf, MAXLINELEN, rerunfile)) {
    print_ver();
    fputs("Error: Empty log file for --rerun.\n", stderr);
    goto rerun_ret_INVALID_FORMAT;
  }
  if (!g_textbuf[MAXLINELEN - 1]) {
    goto rerun_ret_LONG_LINE;
  }
  if (!fgets(g_textbuf, MAXLINELEN, rerunfile)) {
    print_ver();
    fputs("Error: Only one line in --rerun log file.\n", stderr);
    goto rerun_ret_INVALID_FORMAT;
  }
  line_idx++;
  if (!g_textbuf[MAXLINELEN - 1]) {
    goto rerun_ret_LONG_LINE;
  }
  if ((g_textbuf[0] >= '0') && (g_textbuf[0] <= '9')) {
    // Old "xx arguments: --aa bb --cc --dd" format
    fclose_null(&rerunfile);
    if (scan_posint_cappedx(g_textbuf, (MAXLINELEN / 2), &loaded_arg_ct)) {
      print_ver();
      fputs("Error: Invalid argument count on line 2 of --rerun log file.\n", stderr);
      goto rerun_ret_INVALID_FORMAT;
    }
    line_byte_ct = strlen(g_textbuf) + 1;
    rerun_buf = (char*)malloc(line_byte_ct);
    if (!rerun_buf) {
      goto rerun_ret_NOMEM;
    }
    *rerun_buf_ptr = rerun_buf;
    memcpy(rerun_buf, g_textbuf, line_byte_ct);
    // skip "xx arguments: ", to get to the first flag
    rerun_start_ptr = next_token_mult(rerun_buf, 2);
  } else {
    // Current, and also PLINK 1.07, "Options in effect:"
    while (memcmp(g_textbuf, "Options in effect:", 18) || (g_textbuf[18] >= ' ')) {
      line_idx++;
      if (!fgets(g_textbuf, MAXLINELEN, rerunfile)) {
	print_ver();
	fputs("Error: Invalid log file for --rerun.\n", stderr);
	goto rerun_ret_INVALID_FORMAT;
      }
    }
    load_ptr = g_textbuf;
    loaded_arg_ct = 0;
    // We load each of the option lines in sequence into g_textbuf, always
    // overwriting the previous line's newline.  (Note that g_textbuf[] has
    // size > 2 * MAXLINELEN; this lets us avoid additional dynamic memory
    // allocation as long as we impose the constraint that all lines combined
    // add up to less than MAXLINELEN bytes.)
    while (1) {
      load_ptr[MAXLINELEN - 1] = ' ';
      if (!fgets(load_ptr, MAXLINELEN, rerunfile)) {
	break;
      }
      line_idx++;
      if (!load_ptr[MAXLINELEN - 1]) {
	goto rerun_ret_LONG_LINE;
      }
      sptr = skip_initial_spaces(load_ptr);
      if (is_eoln_kns(*sptr)) {
	*load_ptr = '\0';
	break;
      }
      do {
	argptr = token_endnn(sptr);
	loaded_arg_ct++;
        sptr = skip_initial_spaces(argptr);
      } while (!is_eoln_kns(*sptr));
      load_ptr = argptr;
      if (load_ptr >= &(g_textbuf[MAXLINELEN])) {
	print_ver();
	fputs("Error: --rerun argument sequence too long.\n", stderr);
	goto rerun_ret_INVALID_FORMAT;
      }
    }
    fclose_null(&rerunfile);
    line_byte_ct = 1 + (uintptr_t)(load_ptr - g_textbuf);
    rerun_buf = (char*)malloc(line_byte_ct);
    if (!rerun_buf) {
      goto rerun_ret_NOMEM;
    }
    memcpy(rerun_buf, g_textbuf, line_byte_ct);
    rerun_start_ptr = skip_initial_spaces(rerun_buf);
  }
  sptr = rerun_start_ptr;

  // now use g_textbuf as a lame bitfield
  memset(g_textbuf, 1, loaded_arg_ct);
  loaded_arg_idx = 0;
  duplicate_ct = 0;
  do {
    if (no_more_tokens_kns(sptr)) {
      print_ver();
      fputs("Error: Line 2 of --rerun log file has fewer tokens than expected.\n", stderr);
      goto rerun_ret_INVALID_FORMAT;
    }
    argptr = is_flag_start(sptr);
    if (argptr) {
      slen = strlen_se(argptr);
      for (cmdline_arg_idx = cur_arg; cmdline_arg_idx < argc; cmdline_arg_idx++) {
	argptr2 = is_flag_start(argv[cmdline_arg_idx]);
	if (argptr2) {
	  slen2 = strlen(argptr2);
	  if ((slen == slen2) && (!memcmp(argptr, argptr2, slen))) {
	    cmdline_arg_idx = 0xffffffffU;
	    break;
	  }
	}
      }
      if (cmdline_arg_idx == 0xffffffffU) {
	// matching flag, override --rerun
	do {
	  duplicate_ct++;
	  g_textbuf[loaded_arg_idx++] = 0;
	  if (loaded_arg_idx == loaded_arg_ct) {
	    break;
	  }
	  sptr = next_token(sptr);
	} while (!is_flag(sptr));
      } else {
	loaded_arg_idx++;
	sptr = next_token(sptr);
      }
    } else {
      loaded_arg_idx++;
      sptr = next_token(sptr);
    }
  } while (loaded_arg_idx < loaded_arg_ct);
  subst_argv2 = (char**)malloc((argc + loaded_arg_ct - duplicate_ct - rerun_parameter_present - 1 - cur_arg) * sizeof(char*));
  if (!subst_argv2) {
    goto rerun_ret_NOMEM;
  }
  new_arg_idx = 0;
  for (cmdline_arg_idx = cur_arg; cmdline_arg_idx < rerun_argv_pos; cmdline_arg_idx++) {
    subst_argv2[new_arg_idx++] = argv[cmdline_arg_idx];
  }
  sptr = rerun_start_ptr;
  for (loaded_arg_idx = 0; loaded_arg_idx < loaded_arg_ct; loaded_arg_idx++) {
    if (g_textbuf[loaded_arg_idx]) {
      slen = strlen_se(sptr);
      subst_argv2[new_arg_idx++] = sptr;
      sptr[slen] = '\0';
      if (loaded_arg_idx != loaded_arg_ct - 1) {
	sptr = skip_initial_spaces(&(sptr[slen + 1]));
      }
    } else {
      sptr = next_token(sptr);
    }
  }
  for (cmdline_arg_idx = rerun_argv_pos + rerun_parameter_present + 1; cmdline_arg_idx < argc; cmdline_arg_idx++) {
    subst_argv2[new_arg_idx++] = argv[cmdline_arg_idx];
  }
  *cur_arg_ptr = 0;
  *argc_ptr = new_arg_idx;
  if (*subst_argv_ptr) {
    free(*subst_argv_ptr);
  }
  *subst_argv_ptr = subst_argv2;
  *argv_ptr = subst_argv2;
  subst_argv2 = nullptr;
  while (0) {
  rerun_ret_NOMEM:
    print_ver();
    fputs(errstr_nomem, stderr);
    retval = RET_NOMEM;
    break;
  rerun_ret_OPEN_FAIL:
    print_ver();
    fputs("Error: Failed to open --rerun file.\n", stderr);
    retval = RET_OPEN_FAIL;
    break;
  rerun_ret_LONG_LINE:
    print_ver();
    fprintf(stderr, "Error: Line %" PRIuPTR " of --rerun log file is pathologically long.\n", line_idx);
  rerun_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  free_cond(subst_argv2);
  fclose_cond(rerunfile);
  return retval;
}

char extract_char_param(char* ss) {
  // maps c, 'c', and "c" to c, and anything else to the null char.  This is
  // intended to support e.g. always using '#' to designate a # parameter
  // without worrying about differences between shells.
  char cc = ss[0];
  if (((cc == '\'') || (cc == '"')) && (ss[1]) && (ss[2] == cc) && (!ss[3])) {
    return ss[1];
  } else if (cc && (!ss[1])) {
    return cc;
  } else {
    return '\0';
  }
}

int32_t alloc_string(char** sbuf, const char* source) {
  uint32_t slen = strlen(source) + 1;
  *sbuf = (char*)malloc(slen * sizeof(char));
  if (!(*sbuf)) {
    return -1;
  }
  memcpy(*sbuf, source, slen);
  return 0;
}

int32_t alloc_fname(char** fnbuf, char* source, char* argptr, uint32_t extra_size) {
  uint32_t slen = strlen(source) + 1;
  if (slen > (FNAMESIZE - extra_size)) {
    LOGERRPRINTF("Error: --%s filename too long.\n", argptr);
    return RET_OPEN_FAIL;
  }
  *fnbuf = (char*)malloc((slen + extra_size) * sizeof(char));
  if (!(*fnbuf)) {
    return RET_NOMEM;
  }
  memcpy(*fnbuf, source, slen);
  return 0;
}

int32_t alloc_and_flatten(char** flattened_buf_ptr, char** sources, uint32_t param_ct) {
  uint32_t totlen = 1;
  char* bufptr;
  uint32_t param_idx;
  for (param_idx = 0; param_idx < param_ct; param_idx++) {
    totlen += 1 + strlen(sources[param_idx]);
  }
  bufptr = (char*)malloc(totlen);
  if (!bufptr) {
    return RET_NOMEM;
  }
  *flattened_buf_ptr = bufptr;
  for (param_idx = 0; param_idx < param_ct; param_idx++) {
    bufptr = strcpyax(bufptr, sources[param_idx], '\0');
  }
  *bufptr = '\0';
  return 0;
}

int32_t alloc_and_flatten_comma_delim(char** flattened_buf_ptr, char** sources, uint32_t param_ct) {
  uint32_t totlen = 1;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  uint32_t param_idx;
  for (param_idx = 0; param_idx < param_ct; param_idx++) {
    bufptr = sources[param_idx];
    while (1) {
      while (*bufptr == ',') {
	bufptr++;
      }
      bufptr2 = strchr(bufptr, ',');
      if (!bufptr2) {
	break;
      }
      totlen += 1 + (uintptr_t)(bufptr2 - bufptr);
      bufptr = &(bufptr2[1]);
    }
    totlen += 1 + strlen(bufptr);
  }
  bufptr = (char*)malloc(totlen);
  if (!bufptr) {
    return RET_NOMEM;
  }
  *flattened_buf_ptr = bufptr;
  for (param_idx = 0; param_idx < param_ct; param_idx++) {
    bufptr2 = sources[param_idx];
    while (1) {
      while (*bufptr2 == ',') {
	bufptr2++;
      }
      bufptr3 = strchr(bufptr2, ',');
      if (!bufptr3) {
	break;
      }
      bufptr = memcpyax(bufptr, bufptr2, (uintptr_t)(bufptr3 - bufptr2), '\0');
      bufptr2 = &(bufptr3[1]);
    }
    bufptr = strcpyax(bufptr, bufptr2, '\0');
  }
  *bufptr = '\0';
  return 0;
}

int32_t alloc_2col(Two_col_params** tcbuf, char** params_ptr, char* argptr, uint32_t param_ct) {
  uint32_t slen = strlen(*params_ptr) + 1;
  char cc;
  if (slen > FNAMESIZE) {
    LOGERRPRINTF("Error: --%s filename too long.\n", argptr);
    return RET_OPEN_FAIL;
  }
  *tcbuf = (Two_col_params*)malloc(sizeof(Two_col_params) + slen);
  if (!(*tcbuf)) {
    return RET_NOMEM;
  }
  memcpy((*tcbuf)->fname, params_ptr[0], slen);
  (*tcbuf)->skip = 0;
  (*tcbuf)->skipchar = '\0';
  if (param_ct > 1) {
    if (scan_posint_defcapx(params_ptr[1], &((*tcbuf)->colx))) {
      LOGERRPRINTF("Error: Invalid --%s column number.\n", argptr);
      return RET_INVALID_FORMAT;
    }
    if (param_ct > 2) {
      if (scan_posint_defcapx(params_ptr[2], &((*tcbuf)->colid))) {
	LOGERRPRINTF("Error: Invalid --%s variant ID column number.\n", argptr);
	return RET_INVALID_FORMAT;
      }
      if (param_ct == 4) {
	cc = params_ptr[3][0];
	if ((cc < '0') || (cc > '9')) {
	  cc = extract_char_param(params_ptr[3]);
	  if (!cc) {
            goto alloc_2col_invalid_skip;
	  }
	  (*tcbuf)->skipchar = cc;
	} else {
	  if (scan_uint_defcapx(params_ptr[3], &((*tcbuf)->skip))) {
	  alloc_2col_invalid_skip:
	    LOGERRPRINTF("Error: Invalid --%s skip parameter.  This needs to either be a\nsingle character (usually '#') which, when present at the start of a line,\nindicates it should be skipped; or the number of initial lines to skip.  (Note\nthat in shells such as bash, '#' is a special character that must be\nsurrounded by single- or double-quotes to be parsed correctly.)\n", argptr);
	    return RET_INVALID_FORMAT;
	  }
	}
      }
    } else {
      (*tcbuf)->colid = 1;
    }
    if ((*tcbuf)->colx == (*tcbuf)->colid) {
      LOGERRPRINTF("Error: Column numbers for --%s cannot be equal.\n%s", argptr, errstr_append);
      return RET_INVALID_FORMAT;
    }
  } else {
    (*tcbuf)->colx = 2;
    (*tcbuf)->colid = 1;
  }
  return 0;
}

int32_t flag_match(const char* to_match, uint32_t* cur_flag_ptr, uint32_t flag_ct, char* flag_buf) {
  int32_t ii;
  while (*cur_flag_ptr < flag_ct) {
    ii = strcmp(to_match, &(flag_buf[(*cur_flag_ptr) * MAX_FLAG_LEN]));
    if (ii < 0) {
      return 0;
    }
    *cur_flag_ptr += 1;
    if (!ii) {
      flag_buf[((*cur_flag_ptr) - 1) * MAX_FLAG_LEN] = '\0';
      return 1;
    }
  }
  return 0;
}

uint32_t species_flag(uint32_t* species_code_ptr, uint32_t new_code) {
  if (*species_code_ptr) {
    logerrprint("Error: Multiple chromosome set flags.\n");
    return 1;
  }
  *species_code_ptr = new_code;
  return 0;
}

uint32_t valid_varid_template_string(char* varid_str, const char* flag_name) {
  char* sptr = strchr(varid_str, '@');
  char* sptr2 = strchr(varid_str, '#');
  if ((!sptr) || (!sptr2) || strchr(&(sptr[1]), '@') || strchr(&(sptr2[1]), '#')) {
    LOGERRPRINTFWW("Error: The %s template string requires exactly one '@' and one '#'.\n", flag_name);
    return 0;
  }
  // snp/nonsnp is not sufficient for assigning unique IDs to unnamed 1000
  // Genomes phase 3 variants (see e.g. chr22:18078898).  So we now allow the
  // template string to include allele names, where '$1' = first allele in
  // ASCII-sort order, and '$2' = second allele.
  // Use of sorted order, instead of ref/alt or major/minor, is necessary to
  // minimize the chance of assigning different names to the same variant in
  // different datasets.  (Triallelic variants are still a problem, but no
  // solution is possible there until the file format is extended.)
  // For now, either '$' must be entirely absent from the template string, or
  // it appears exactly twice, once in '$1' and once in '$2'.
  sptr = strchr(varid_str, '$');
  if (sptr) {
    sptr2 = strchr(&(sptr[1]), '$');
    if ((!sptr2) || strchr(&(sptr2[1]), '$') || (!(((sptr[1] == '1') && (sptr2[1] == '2')) || ((sptr[1] == '2') && (sptr2[1] == '1'))))) {
      LOGERRPRINTFWW("Error: The %s template string requires either no instances of '$', or exactly one '$1' and one '$2'.\n", flag_name);
      return 0;
    }
  }
  return 1;
}

int32_t init_delim_and_species(uint32_t flag_ct, char* flag_buf, uint32_t* flag_map, int32_t argc, char** argv, char* range_delim_ptr, Chrom_info* chrom_info_ptr) {
  uint32_t species_code = SPECIES_DEFAULT;
  uint32_t flag_idx = 0;
  uint32_t retval = 0;
  int32_t cur_arg;
  uint32_t param_ct;
  int32_t ii;
  uint32_t param_idx;
  if (flag_match("autosome-num", &flag_idx, flag_ct, flag_buf)) {
    species_code = SPECIES_UNKNOWN;
    cur_arg = flag_map[flag_idx - 1];
    param_ct = param_count(argc, argv, cur_arg);
    if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE_2A;
    }
    if (scan_posint_cappedx(argv[cur_arg + 1], MAX_CHROM_TEXTNUM, (uint32_t*)(&ii))) {
      sprintf(g_logbuf, "Error: Invalid --autosome-num parameter '%s'.\n", argv[cur_arg + 1]);
      goto init_delim_and_species_ret_INVALID_CMDLINE_WWA;
    }
    chrom_info_ptr->xymt_codes[X_OFFSET] = ii + 1;
    // bugfix (28 Apr 2020): We were segfaulting instead of printing an
    // appropriate error message when e.g. --autosome-num was specified yet
    // chrY was encountered.
    // plink 2.0 fixed this issue by using -2 instead of -1 to indicate a null
    // xymt_code (-1 signals that there's nothing to worry about when
    // --allow-extra-chr is specified).  It touches more code than I'd like,
    // but I haven't gotten any related plink 2.0 bug reports for a long time,
    // so we backport that fix.
    chrom_info_ptr->xymt_codes[Y_OFFSET] = -2;
    chrom_info_ptr->xymt_codes[XY_OFFSET] = -2;
    chrom_info_ptr->xymt_codes[MT_OFFSET] = -2;
    chrom_info_ptr->max_code = ii + 1;
    chrom_info_ptr->autosome_ct = ii;
    set_bit(ii + 1, chrom_info_ptr->haploid_mask);
  }
  if (flag_match("chr-set", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_UNKNOWN)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    cur_arg = flag_map[flag_idx - 1];
    param_ct = param_count(argc, argv, cur_arg);
    if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 5)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE_2A;
    }
    if (scan_int_abs_boundedx(argv[cur_arg + 1], MAX_CHROM_TEXTNUM, &ii) || (!ii)) {
      sprintf(g_logbuf, "Error: Invalid --chr-set parameter '%s'.\n", argv[cur_arg + 1]);
      goto init_delim_and_species_ret_INVALID_CMDLINE_WWA;
    }
    if (ii < 0) {
      if (param_ct > 1) {
	logerrprint("Error: --chr-set does not accept multiple parameters in haploid mode.\n");
	goto init_delim_and_species_ret_INVALID_CMDLINE_A;
      }
      ii = -ii;
      chrom_info_ptr->autosome_ct = ii;
      chrom_info_ptr->xymt_codes[X_OFFSET] = -2;
      chrom_info_ptr->xymt_codes[Y_OFFSET] = -2;
      chrom_info_ptr->xymt_codes[XY_OFFSET] = -2;
      chrom_info_ptr->xymt_codes[MT_OFFSET] = -2;
      chrom_info_ptr->max_code = ii;
      fill_all_bits(((uint32_t)ii) + 1, chrom_info_ptr->haploid_mask);
    } else {
      chrom_info_ptr->autosome_ct = ii;
      chrom_info_ptr->xymt_codes[X_OFFSET] = ii + 1;
      chrom_info_ptr->xymt_codes[Y_OFFSET] = ii + 2;
      chrom_info_ptr->xymt_codes[XY_OFFSET] = ii + 3;
      chrom_info_ptr->xymt_codes[MT_OFFSET] = ii + 4;
      set_bit(ii + 1, chrom_info_ptr->haploid_mask);
      set_bit(ii + 2, chrom_info_ptr->haploid_mask);
      for (param_idx = 2; param_idx <= param_ct; param_idx++) {
	if (!strcmp(argv[cur_arg + param_idx], "no-x")) {
	  chrom_info_ptr->xymt_codes[X_OFFSET] = -2;
	  clear_bit(ii + 1, chrom_info_ptr->haploid_mask);
	} else if (!strcmp(argv[cur_arg + param_idx], "no-y")) {
	  chrom_info_ptr->xymt_codes[Y_OFFSET] = -2;
	  clear_bit(ii + 2, chrom_info_ptr->haploid_mask);
	} else if (!strcmp(argv[cur_arg + param_idx], "no-xy")) {
	  chrom_info_ptr->xymt_codes[XY_OFFSET] = -2;
	} else if (!strcmp(argv[cur_arg + param_idx], "no-mt")) {
	  chrom_info_ptr->xymt_codes[MT_OFFSET] = -2;
	} else {
	  sprintf(g_logbuf, "Error: Invalid --chr-set parameter '%s'.\n", argv[cur_arg + param_idx]);
	  goto init_delim_and_species_ret_INVALID_CMDLINE_WWA;
	}
      }
      if (chrom_info_ptr->xymt_codes[MT_OFFSET] != -2) {
	chrom_info_ptr->max_code = ii + 4;
      } else if (chrom_info_ptr->xymt_codes[XY_OFFSET] != -2) {
	chrom_info_ptr->max_code = ii + 3;
      } else if (chrom_info_ptr->xymt_codes[Y_OFFSET] != -2) {
	chrom_info_ptr->max_code = ii + 2;
      } else if (chrom_info_ptr->xymt_codes[X_OFFSET] != -2) {
	chrom_info_ptr->max_code = ii + 1;
      } else {
	chrom_info_ptr->max_code = ii;
      }
    }
  }
  if (flag_match("cow", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_COW)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logerrprint("Error: --cow doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("d", &flag_idx, flag_ct, flag_buf)) {
    // moved here to support --covar-name + --d
    cur_arg = flag_map[flag_idx - 1];
    param_ct = param_count(argc, argv, cur_arg);
    if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE_2A;
    }
    *range_delim_ptr = extract_char_param(argv[cur_arg + 1]);
    if (!(*range_delim_ptr)) {
      logerrprint("Error: --d parameter too long (must be a single character).\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE_A;
    } else if ((*range_delim_ptr == '-') || (*range_delim_ptr == ',')) {
      logerrprint("Error: --d parameter cannot be '-' or ','.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE_A;
    }
  }
  if (flag_match("dog", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_DOG)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logerrprint("Error: --dog doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("horse", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_HORSE)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logerrprint("Error: --horse doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("mouse", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_MOUSE)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logerrprint("Error: --mouse doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("rice", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_RICE)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logerrprint("Error: --rice doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("sheep", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_SHEEP)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logerrprint("Error: --sheep doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  init_species(species_code, chrom_info_ptr);
  while (0) {
  init_delim_and_species_ret_INVALID_CMDLINE_WWA:
    wordwrapb(0);
  init_delim_and_species_ret_INVALID_CMDLINE_2A:
    logerrprintb();
  init_delim_and_species_ret_INVALID_CMDLINE_A:
    logerrprint(errstr_append);
    retval = RET_INVALID_CMDLINE;
    break;
  init_delim_and_species_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

int32_t recode_type_set(uint32_t* recode_modifier_ptr, uint32_t cur_code) {
  if (*recode_modifier_ptr & (RECODE_TYPEMASK - cur_code)) {
    logerrprint("Error: Multiple --recode output formats.\n");
    return -1;
  }
  *recode_modifier_ptr |= cur_code;
  return 0;
}

int32_t main(int32_t argc, char** argv) {
  char* outname_end = nullptr;
  char** subst_argv = nullptr;
  char* script_buf = nullptr;
  char* rerun_buf = nullptr;
  char* flag_buf = nullptr;
  uint32_t* flag_map = nullptr;
  char* makepheno_str = nullptr;
  char* phenoname_str = nullptr;
  Two_col_params* a1alleles = nullptr;
  Two_col_params* a2alleles = nullptr;
  char* sample_sort_fname = nullptr;
  char* filtervals_flattened = nullptr;
  char* evecname = nullptr;
  char* filtername = nullptr;
  char* distance_wts_fname = nullptr;
  char* read_dists_fname = nullptr;
  char* read_dists_id_fname = nullptr;
  char* freqname = nullptr;
  char* extractname = nullptr;
  char* excludename = nullptr;
  char* keepname = nullptr;
  char* removename = nullptr;
  char* keepfamname = nullptr;
  char* removefamname = nullptr;
  char* cm_map_fname = nullptr;
  char* cm_map_chrname = nullptr;
  char* phenoname = nullptr;
  char* recode_allele_name = nullptr;
  char* lgen_reference_fname = nullptr;
  char* covar_fname = nullptr;
  char* update_alleles_fname = nullptr;
  Two_col_params* qual_filter = nullptr;
  Two_col_params* update_chr = nullptr;
  Two_col_params* update_cm = nullptr;
  Two_col_params* update_map = nullptr;
  Two_col_params* update_name = nullptr;
  char* oxford_single_chr = nullptr;
  char* oxford_pheno_name = nullptr;
  char* update_ids_fname = nullptr;
  char* update_parents_fname = nullptr;
  char* update_sex_fname = nullptr;
  char* loop_assoc_fname = nullptr;
  char* flip_fname = nullptr;
  char* flip_subset_fname = nullptr;
  char* read_genome_fname = nullptr;
  char* condition_mname = nullptr;
  char* condition_fname = nullptr;
  char* missing_marker_id_match = nullptr;
  char* filter_attrib_fname = nullptr;
  char* filter_attrib_liststr = nullptr;
  char* filter_attrib_sample_fname = nullptr;
  char* filter_attrib_sample_liststr = nullptr;
  char* const_fid = nullptr;
  char* vcf_filter_exceptions_flattened = nullptr;
  char* gene_report_fname = nullptr;
  char* gene_report_glist = nullptr;
  char* gene_report_subset = nullptr;
  char* gene_report_snp_field = nullptr;
  char* metaanal_fnames = nullptr;
  char* metaanal_chrfield_search_order = nullptr;
  char* metaanal_snpfield_search_order = nullptr;
  char* metaanal_bpfield_search_order = nullptr;
  char* metaanal_a1field_search_order = nullptr;
  char* metaanal_a2field_search_order = nullptr;
  char* metaanal_pfield_search_order = nullptr;
  char* metaanal_sefield_search_order = nullptr;
  char* metaanal_essfield_search_order = nullptr;
  char* rplugin_fname = nullptr;
  char* rplugin_host_or_socket = nullptr;
  uint32_t gene_report_border = 0;
  uint32_t metaanal_flags = 0;
  int32_t rplugin_port = -2; // now interpreted by rserve_call() as default
  double vcf_min_qual = -1;
  double vcf_min_gq = -1;
  double vcf_min_gp = -1;
  double qual_min_thresh = 0.0;
  double qual_max_thresh = DBL_MAX;
  char id_delim = '\0';
  char vcf_idspace_to = '\0';
  unsigned char vcf_half_call = 0;
  int32_t retval = 0;
  uint32_t load_params = 0; // describes what file params have been provided
  uint32_t load_rare = 0;
  uint32_t fam_cols = FAM_COL_13456;
  uint32_t mpheno_col = 0;
  uint32_t mwithin_col = 0;
  uint64_t misc_flags = 0;
  uint64_t filter_flags = 0;
  double thin_keep_prob = 1.0;
  double thin_keep_sample_prob = 1.0;
  uint32_t thin_keep_ct = 0xffffffffU;
  uint32_t thin_keep_sample_ct = 0xffffffffU;
  uint32_t min_bp_space = 0;
  uint32_t check_sex_f_yobs = 0;
  uint32_t check_sex_m_yobs = 0;
  double check_sex_fthresh = 0.2;
  double check_sex_mthresh = 0.8;
  double distance_exp = 0.0;
  double min_maf = 0.0;
  double max_maf = 0.5;
  double geno_thresh = 1.0;
  double mind_thresh = 1.0;
  double hwe_thresh = 0.0;
  uint32_t cur_arg = 1;
  uint64_t calculation_type = 0;
  uint32_t dist_calc_type = 0;
  uint32_t mfilter_col = 0;
  uint32_t pheno_modifier = 0;
  int32_t missing_pheno = -9;
  intptr_t malloc_size_mb = 0;
  uintptr_t groupdist_iters = ITERS_DEFAULT;
  uint32_t groupdist_d = 0;
  uintptr_t regress_iters = ITERS_DEFAULT;
  uint32_t regress_d = 0;
  uint32_t parallel_idx = 0;
  uint32_t parallel_tot = 1;
  uint32_t splitx_bound1 = 0;
  uint32_t splitx_bound2 = 0;
  uint32_t sex_missing_pheno = 0;
  uint32_t update_sex_col = 1;
  uint32_t hwe_modifier = 0;
  uint32_t min_ac = 0;
  uint32_t max_ac = 0x7fffffff;
  uint32_t write_covar_modifier = 0;
  uint32_t write_covar_dummy_max_categories = 49;
  uint32_t dupvar_modifier = 0;
  uint32_t model_modifier = 0;
  int32_t model_cell_ct = -1;
  uint32_t gxe_mcovar = 0;
  uint32_t glm_modifier = 0;
  double glm_vif_thresh = 50.0;
  uint32_t glm_xchr_model = 1;
  uint32_t ppc_gap = DEFAULT_PPC_GAP;
  uint32_t* rseeds = nullptr;
  uint32_t rseed_ct = 0;
  uint32_t genome_modifier = 0;
  double genome_min_pi_hat = -1.0;
  double genome_max_pi_hat = 1.0;
  FILE* scriptfile = nullptr;
  uint32_t recode_modifier = 0;
  uint32_t recode_require_format = 0;
  uint32_t allelexxxx = 0;
  uint32_t merge_type = 0;
  uint32_t sample_sort = 0;
  uint32_t cur_flag = 0;
  uint32_t flag_ct = 0;
  uint32_t dummy_marker_ct = 0;
  uint32_t dummy_sample_ct = 0;
  uint32_t dummy_flags = 0;
  double dummy_missing_geno = 0.0;
  double dummy_missing_pheno = 0.0;
  char* simulate_fname = nullptr;
  uint32_t simulate_flags = 0;
  uint32_t simulate_cases = 1000;
  uint32_t simulate_controls = 1000;
  double simulate_prevalence = 0.01;
  char* simulate_label = nullptr;
  double simulate_missing = 0.0;
  uint32_t simulate_qt_samples = 1000;
  char* markername_from = nullptr;
  char* markername_to = nullptr;
  char* markername_snp = nullptr;
  // minor bugfix: --snp and '--window 0 --snp' should actually behave
  // differently when there are other variants at the same bp coordinate
  int32_t snp_window_size = -1;
  int32_t marker_pos_start = -1;
  int32_t marker_pos_end = -1;
  uint32_t write_var_range_ct = 0;
  uint32_t lgen_modifier = 0;
  uint32_t covar_modifier = 0;
  uint32_t update_map_modifier = 0;
  uint32_t model_mperm_val = 0;
  uint32_t glm_mperm_val = 0;
  uint32_t mperm_save = 0;
  uint32_t mperm_val = 0;
  double ci_size = 0.0;
  double pfilter = 2.0; // make --pfilter 1 still filter out NAs
  double output_min_p = 0.0;
  uint32_t perm_batch_size = DEFAULT_PERM_BATCH_SIZE;
  uint32_t mtest_adjust = 0;
  double adjust_lambda = 0.0;
  uintptr_t ibs_test_perms = DEFAULT_IBS_TEST_PERMS;
  uint32_t neighbor_n1 = 0;
  uint32_t neighbor_n2 = 0;
  uint32_t cnv_calc_type = 0;
  uint32_t cnv_sample_mperms = 0;
  uint32_t cnv_test_mperms = 0;
  uint32_t cnv_test_region_mperms = 0;
  uint32_t cnv_enrichment_test_mperms = 0;
  uint32_t cnv_min_seglen = 0;
  uint32_t cnv_max_seglen = 0xffffffffU;
  double cnv_min_score = -DBL_MAX;
  double cnv_max_score = DBL_MAX;
  uint32_t cnv_min_sites = 0;
  uint32_t cnv_max_sites = 0xffffffffU;
  uint32_t cnv_intersect_filter_type = 0;
  char* cnv_intersect_filter_fname = nullptr;
  char* cnv_subset_fname = nullptr;
  uint32_t cnv_overlap_type = 0;
  double cnv_overlap_val = 0.0;
  uint32_t cnv_freq_type = 0;
  uint32_t cnv_freq_val = 0;
  double cnv_freq_val2 = 0.0;
  uint32_t cnv_test_window = 0;
  uint32_t segment_modifier = 0;
  uint32_t permphe_ct = 0;
  uint32_t matrix_flag_state = 0; // 1 = present and unclaimed, 2 = claimed
  double hard_call_threshold = 0.1;
  double tail_bottom = 0.0;
  double tail_top = 0.0;
  double lasso_h2 = 0.0;
  double lasso_minlambda = -1;
  uint32_t lasso_lambda_iters = 0;
  uint32_t testmiss_modifier = 0;
  uint32_t testmiss_mperm_val = 0;

  // this default limit plays well with e.g. fbstring small-string optimization
  uint32_t new_id_max_allele_len = 23;

  uint32_t aperm_present = 0;
  char* segment_spanning_fname = nullptr;
  char* oxford_missing_code = nullptr;
  char range_delim = '-';
  uint32_t modifier_23 = 0;
  double pheno_23 = DBL_MAX;
  char* fid_23 = nullptr;
  char* iid_23 = nullptr;
  char* paternal_id_23 = nullptr;
  char* maternal_id_23 = nullptr;
  Ll_str* file_delete_list = nullptr;
  uint32_t chrom_flag_present = 0;
  uintptr_t chrom_exclude[CHROM_MASK_INITIAL_WORDS];
  // er, except for first four, these should not be preallocated...
  char outname[FNAMESIZE];
  char mapname[FNAMESIZE];
  char pedname[FNAMESIZE];
  char famname[FNAMESIZE];
  char mergename1[FNAMESIZE];
  char mergename2[FNAMESIZE];
  char mergename3[FNAMESIZE];
  char output_missing_pheno[32];
  char* missing_mid_template;
#ifdef __APPLE__
  int32_t mib[2];
  size_t sztmp;
#endif
  unsigned char* bigstack_ua = nullptr; // ua = unaligned
  char* bubble = nullptr;
  unsigned char* bigstack_initial_base;
  uint32_t param_ct;
  time_t rawtime;
  char* argptr;
  char* sptr;
  const char* csptr;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t num_params;
  int32_t in_param;
  Chrom_info chrom_info;
  Family_info family_info;
  Oblig_missing_info oblig_missing_info;
  Aperm_info aperm;
  Cluster_info cluster;
  Set_info set_info;
  Homozyg_info homozyg;
  Ld_info ld_info;
  Epi_info epi_info;
  Annot_info annot_info;
  Clump_info clump_info;
  Rel_info rel_info;
  Score_info score_info;
  Dosage_info dosage_info;
  Range_list snps_range_list;
  Range_list covar_range_list;
  Range_list lasso_select_covars_range_list;
  Range_list parameters_range_list;
  Range_list tests_range_list;
  char* argptr2;
  char* flagptr;
  double dxx;
  char cc;
  int32_t known_procs;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  intptr_t default_alloc_mb;
  int64_t llxx;
  Ll_str* ll_str_ptr;
#ifdef _WIN32
  SYSTEM_INFO sysinfo;
  MEMORYSTATUSEX memstatus;
  DWORD windows_dw; // why the f*** does uint32_t not work?
#endif
  family_init(&family_info);
  oblig_missing_init(&oblig_missing_info);
  aperm_init(&aperm);
  cluster_init(&cluster);
  set_init(&set_info, &annot_info);
  homozyg_init(&homozyg);
  ld_epi_init(&ld_info, &epi_info, &clump_info);
  rel_init(&rel_info);
  misc_init(&score_info);
  dosage_init(&dosage_info);
  range_list_init(&snps_range_list);
  range_list_init(&covar_range_list);
  range_list_init(&lasso_select_covars_range_list);
  range_list_init(&parameters_range_list);
  range_list_init(&tests_range_list);
  missing_mid_template = nullptr;

  // standardize strtod() behavior
  setlocale(LC_NUMERIC, "C");

  if (init_chrom_info(&chrom_info)) {
    goto main_ret_NOMEM_NOLOG;
  }

  for (uii = 1; uii < (uint32_t)argc; uii++) {
    if ((!strcmp("-script", argv[uii])) || (!strcmp("--script", argv[uii]))) {
      ujj = param_count(argc, argv, uii);
      if (enforce_param_ct_range(ujj, argv[uii], 1, 1)) {
	print_ver();
	fputs(g_logbuf, stderr);
	fputs(errstr_append, stderr);
	goto main_ret_INVALID_CMDLINE;
      }
      for (ujj = uii + 2; ujj < (uint32_t)argc; ujj++) {
	if ((!strcmp("-script", argv[ujj])) || (!strcmp("--script", argv[ujj]))) {
	  print_ver();
	  fputs("Error: Multiple --script flags.  Merge the files into one.\n", stderr);
	  fputs(errstr_append, stderr);
	  goto main_ret_INVALID_CMDLINE;
	}
      }
      // logging not yet active, so don't use fopen_checked()
      scriptfile = fopen(argv[uii + 1], FOPEN_RB);
      if (!scriptfile) {
	print_ver();
	fprintf(stderr, g_errstr_fopen, argv[uii + 1]);
	goto main_ret_OPEN_FAIL;
      }
      if (fseeko(scriptfile, 0, SEEK_END)) {
	goto main_ret_READ_FAIL_NOLOG;
      }
      llxx = ftello(scriptfile);
      if (llxx == -1) {
	goto main_ret_READ_FAIL_NOLOG;
      } else if (llxx > 0x7fffffff) {
	// could actually happen if user enters parameters in the wrong order,
	// so may as well catch it and print a somewhat informative error msg
	print_ver();
        fputs("Error: --script file too large.\n", stderr);
        goto main_ret_INVALID_CMDLINE;
      }
      rewind(scriptfile);
      ujj = (uint32_t)((uint64_t)llxx);
      script_buf = (char*)malloc(ujj);
      if (!script_buf) {
	goto main_ret_NOMEM_NOLOG;
      }
      ukk = fread(script_buf, 1, ujj, scriptfile);
      if (ukk < ujj) {
	goto main_ret_READ_FAIL_NOLOG;
      }
      fclose_null(&scriptfile);
      num_params = 0;
      in_param = 0;
      for (ukk = 0; ukk < ujj; ukk++) {
	if (is_space_or_eoln(script_buf[ukk])) {
	  in_param = 0;
	} else if (!in_param) {
	  num_params++;
	  in_param = 1;
	}
      }
      subst_argv = (char**)malloc((num_params + argc - 3) * sizeof(char*));
      if (!subst_argv) {
	goto main_ret_NOMEM_NOLOG;
      }
      num_params = 0;
      in_param = 0;
      for (ukk = 1; ukk < uii; ukk++) {
        subst_argv[num_params++] = argv[ukk];
      }
      for (ukk = 0; ukk < ujj; ukk++) {
	if (is_space_or_eoln(script_buf[ukk])) {
	  if (in_param) {
	    script_buf[ukk] = '\0';
	    in_param = 0;
	  }
	} else if (!in_param) {
	  subst_argv[num_params++] = &(script_buf[ukk]);
	  in_param = 1;
	}
      }
      for (ujj = uii + 2; ujj < (uint32_t)argc; ujj++) {
	subst_argv[num_params++] = argv[ujj];
      }
      argc = num_params;
      cur_arg = 0;
      argv = subst_argv;
      break;
    }
  }
  for (uii = cur_arg; uii < (uint32_t)argc; uii++) {
    if ((!strcmp("-rerun", argv[uii])) || (!strcmp("--rerun", argv[uii]))) {
      ujj = param_count(argc, argv, uii);
      if (enforce_param_ct_range(ujj, argv[uii], 0, 1)) {
	print_ver();
	fputs(g_logbuf, stderr);
	fputs(errstr_append, stderr);
	goto main_ret_INVALID_CMDLINE;
      }
      for (ukk = uii + ujj + 1; ukk < (uint32_t)argc; ukk++) {
	if ((!strcmp("-rerun", argv[ukk])) || (!strcmp("--rerun", argv[ukk]))) {
	  print_ver();
	  fputs("Error: Duplicate --rerun flag.\n", stderr);
	  goto main_ret_INVALID_CMDLINE;
	}
      }
      retval = rerun(uii, ujj, &argc, &cur_arg, &argv, &subst_argv, &rerun_buf);
      if (retval) {
	goto main_ret_NOLOG;
      }
      break;
    }
  }
  if ((cur_arg < (uint32_t)argc) && (!is_flag(argv[cur_arg]))) {
    print_ver();
    fputs("Error: First parameter must be a flag.\n", stderr);
    fputs(errstr_append, stderr);
    goto main_ret_INVALID_CMDLINE;
  }
  flag_ct = 0;
  // (these aren't resolved immediately since we want --help to take
  // precedence)
  ujj = 0; // --version?
  ukk = 0; // silence?
  for (uii = cur_arg; uii < (uint32_t)argc; uii++) {
    argptr = is_flag_start(argv[uii]);
    if (argptr) {
      if (!strcmp("help", argptr)) {
	print_ver();
	if ((cur_arg != 1) || (uii != 1) || subst_argv) {
	  fputs("--help present, ignoring other flags.\n", stdout);
	}
	if ((uii == ((uint32_t)argc) - 1) && flag_ct) {
	  // make "plink <valid flags/parameters> --help" work, and skip the
	  // parameters
	  char** help_argv = (char**)malloc(flag_ct * sizeof(intptr_t));
	  if (!help_argv) {
	    goto main_ret_NOMEM_NOLOG2;
	  }
	  uint32_t arg_idx2 = 0;
	  for (uint32_t flag_idx = 0; flag_idx < flag_ct; ++flag_idx) {
	    while (!is_flag_start(argv[++arg_idx2]));
	    help_argv[flag_idx] = argv[arg_idx2];
	  }
	  retval = disp_help(flag_ct, help_argv);
	  free(help_argv);
	} else {
	  retval = disp_help(argc - uii - 1, &(argv[uii + 1]));
	}
	goto main_ret_1;
      }
      if ((!strcmp("h", argptr)) || (!strcmp("?", argptr))) {
	// these just act like the no-parameter case
	print_ver();
	if ((cur_arg != 1) || (uii != 1) || subst_argv) {
	  printf("-%s present, ignoring other flags.\n", argptr);
	}
	fputs(g_cmdline_format_str, stdout);
	fputs(notestr_null_calc2, stdout);
	goto main_ret_1;
      }
      if (!strcmp("version", argptr)) {
	ujj = 1;
      } else if ((!strcmp("silent", argptr)) || (!strcmp("gplink", argptr))) {
	ukk = 1;
      }
      if (strlen(argptr) >= MAX_FLAG_LEN) {
	print_ver();
	// shouldn't be possible for this to overflow the buffer...
	sprintf(g_logbuf, "Error: Unrecognized flag ('%s').\n", argv[uii]);
	wordwrapb(0);
	fputs(g_logbuf, stderr);
	fputs(errstr_append, stderr);
        goto main_ret_INVALID_CMDLINE;
      }
      flag_ct++;
    }
  }
  if (!flag_ct) {
    print_ver();
    fputs(g_cmdline_format_str, stdout);
    fputs(notestr_null_calc2, stdout);
    retval = RET_NULL_CALC;
    goto main_ret_1;
  }
  if (ujj) {
    fputs(ver_str, stdout);
    putc_unlocked('\n', stdout);
    goto main_ret_1;
  }
  if (ukk) {
    // bugfix: redirection to /dev/null doesn't work on Windows
    if (!freopen(
#ifdef _WIN32
"nul"
#else
"/dev/null"
#endif
, "w", stdout)) {
      fputs("Warning: --silent failed.\n", stderr);
    }
  }
  print_ver();
  flag_buf = (char*)malloc(flag_ct * MAX_FLAG_LEN * sizeof(char));
  flag_map = (uint32_t*)malloc(flag_ct * sizeof(int32_t));
  if ((!flag_buf) || (!flag_map)) {
    goto main_ret_NOMEM_NOLOG2;
  }
  flagptr = flag_buf;
  for (uii = cur_arg; uii < (uint32_t)argc; uii++) {
    argptr = is_flag_start(argv[uii]);
    if (argptr) {
      ukk = strlen(argptr) + 1;
      // handle aliases now, so sorting will have the desired effects
      switch (*argptr) {
      case '\0':
	// special case, since we reserve empty names for preprocessed flags
	fputs("Error: Unrecognized flag ('--').\n", stderr);
	goto main_ret_INVALID_CMDLINE;
      case 'F':
	if (!strcmp(argptr, "Fst")) {
	  memcpy(flagptr, "fst", 4);
	  break;
	}
	goto main_flag_copy;
      case 'a':
	if ((ukk == 11) && (!memcmp(argptr, "allele", 6))) {
	  if (match_upper(&(argptr[6]), "ACGT")) {
	    memcpy(flagptr, "alleleACGT", 11);
	    break;
	  }
	} else if ((ukk == 12) && (!memcmp(argptr, "allele-", 7))) {
          if (!memcmp(&(argptr[7]), "1234", 4)) {
	    memcpy(flagptr, "allele1234", 11);
	    break;
	  } else if (match_upper(&(argptr[7]), "ACGT")) {
	    memcpy(flagptr, "alleleACGT", 11);
	    break;
	  }
	} else if (!strcmp(argptr, "aec")) {
	  memcpy(flagptr, "allow-extra-chr", 16);
	  break;
	}
	goto main_flag_copy;
      case 'c':
        if (!strcmp(argptr, "chr-excl")) {
          fputs("Note: --chr-excl flag has been renamed to --not-chr.\n", stdout);
	  memcpy(flagptr, "not-chr", 8);
	  break;
	} else if (!strcmp(argptr, "cmh")) {
	  memcpy(flagptr, "mh", 3);
	  break;
	} else if (!strcmp(argptr, "chr-output")) {
	  fputs("Note: --chr-output flag has been renamed to --output-chr, for consistency with\n--output-missing-{genotype/phenotype}.\n", stdout);
	  memcpy(flagptr, "output-chr", 11);
	  break;
	}
	goto main_flag_copy;
      case 'e':
	if (!strcmp(argptr, "extract-snp")) {
	  memcpy(flagptr, "snp", 4);
	  break;
	} else if (!strcmp(argptr, "exponent")) {
	  memcpy(flagptr, "distance-exp", 13);
	  break;
	} else if (!strcmp(argptr, "export")) {
	  memcpy(flagptr, "recode", 7);
	  recode_require_format = 1;
	  break;
	}
	goto main_flag_copy;
      case 'f':
	if (!strcmp(argptr, "frqx")) {
	  memcpy(flagptr, "freqx", 6);
	  break;
	} else if (!strcmp(argptr, "flipscan")) {
	  memcpy(flagptr, "flip-scan", 10);
	  break;
	} else if (!strcmp(argptr, "flipscan-window")) {
	  memcpy(flagptr, "flip-scan-window", 17);
	  break;
	} else if (!strcmp(argptr, "flipscan-window-kb")) {
	  memcpy(flagptr, "flip-scan-window-kb", 20);
	  break;
	} else if (!strcmp(argptr, "flipscan-threshold")) {
	  memcpy(flagptr, "flip-scan-threshold", 20);
	  break;
	}
	goto main_flag_copy;
      case 'g':
	if (!strcmp(argptr, "grm-cutoff")) {
          memcpy(flagptr, "rel-cutoff", 11);
	  break;
	}
	goto main_flag_copy;
      case 'h':
        if (!strcmp(argptr, "hwe2")) {
	  fputs("Warning: --hwe2 flag is obsolete, and now treated as an alias for \"--hwe midp\".\n", stderr);
	  memcpy(flagptr, "hwe midp", 9);
	  break;
        } else if (!strcmp(argptr, "hardy2")) {
	  fputs("Warning: --hardy2 flag is obsolete, and now treated as an alias for\n\"--hardy midp\".\n", stderr);
	  memcpy(flagptr, "hardy midp", 11);
	  break;
	}
	goto main_flag_copy;

      case 'k':
	if (!memcmp(argptr, "k", 2)) {
	  memcpy(flagptr, "K", 2);
	  break;
	}
	goto main_flag_copy;
      case 'l':
	if (!strcmp(argptr, "list")) {
	  memcpy(flagptr, "recode list", 12);
	  fputs("Note: --list flag deprecated.  Use \"--recode list\" instead.\n", stdout);
	  recode_modifier |= RECODE_LIST;
	  misc_flags |= MISC_SET_HH_MISSING;
	  break;
	} else if (!strcmp(argptr, "load-dists")) {
          memcpy(flagptr, "read-dists", 11);
          fputs("Note: --load-dists flag has been renamed to --read-dists.\n", stdout);
          break;
	}
	goto main_flag_copy;
      case 'm':
	if (!strcmp(argptr, "missing_code")) {
	  memcpy(flagptr, "missing-code", 13);
	  break;
	} else if (!strcmp(argptr, "mh1")) {
	  memcpy(flagptr, "mh", 3);
	  break;
	} else if (!strcmp(argptr, "make-set-collapse-all")) {
	  memcpy(flagptr, "set-collapse-all", 17);
	  break;
	} else if (!strcmp(argptr, "min-ac")) {
	  memcpy(flagptr, "mac", 4);
	  break;
	} else if (!strcmp(argptr, "max-ac")) {
	  memcpy(flagptr, "max-mac", 8);
	  break;
	} else if (!strcmp(argptr, "max-indv")) {
	  memcpy(flagptr, "thin-indiv-count", 17);
	  break;
	}
	goto main_flag_copy;
      case 'n':
	if (!strcmp(argptr, "neighbor")) {
	  memcpy(flagptr, "neighbour", 10);
	  break;
	} else if (!strcmp(argptr, "num_threads")) {
	  memcpy(flagptr, "threads", 8);
	  break;
	}
	goto main_flag_copy;
      case 'r':
	if ((ukk >= 8) && (!memcmp(argptr, "recode", 6))) {
	  ujj = 0; // alias match?
	  argptr2 = &(argptr[6]);
          switch (ukk) {
	  case 8:
            if (tolower(*argptr2) == 'a') {
	      memcpy(flagptr, "recode A", 9);
	      recode_modifier |= RECODE_A;
	      ujj = 1;
	    }
	    break;
	  case 9:
	    if (!memcmp(argptr2, "12", 2)) {
              memcpy(flagptr, "recode 12", 10);
	      recode_modifier |= RECODE_12;
	      // technically, this should also affect --tucc, but whatever
	      ujj = 1;
            } else if (match_upper(argptr2, "AD")) {
              memcpy(flagptr, "recode AD", 10);
	      recode_modifier |= RECODE_AD;
	      ujj = 1;
	    } else if (match_upper(argptr2, "HV")) {
	      memcpy(flagptr, "recode HV-1chr", 15);
	      recode_modifier |= RECODE_HV_1CHR;
              fputs("Note: --recodeHV flag deprecated.  Use \"--recode HV\" or \"--recode HV-1chr\".\n", stdout);
	      ujj = 2;
	    }
	    break;
	  case 11:
	    if (!memcmp(argptr2, "-vcf", 4)) {
	      memcpy(flagptr, "recode vcf", 11);
	      recode_modifier |= RECODE_VCF;
	      ujj = 1;
	    }
	    break;
          case 12:
            if (!memcmp(argptr2, "-lgen", 5)) {
              memcpy(flagptr, "recode lgen", 12);
	      recode_modifier |= RECODE_LGEN;
	      // backwards compatibility
	      misc_flags |= MISC_SET_HH_MISSING;
              ujj = 1;
	    }
	    break;
	  case 13:
	    if (!memcmp(argptr2, "-rlist", 6)) {
	      memcpy(flagptr, "recode rlist", 13);
	      recode_modifier |= RECODE_RLIST;
	      misc_flags |= MISC_SET_HH_MISSING;
	      ujj = 1;
	    }
	    break;
	  case 14:
	    if (!memcmp(argptr2, "-beagle", 7)) {
	      memcpy(flagptr, "recode beagle", 14);
	      recode_modifier |= RECODE_BEAGLE;
	      ujj = 1;
	    } else if (!memcmp(argptr2, "-bimbam", 7)) {
	      memcpy(flagptr, "recode bimbam-1chr", 19);
	      recode_modifier |= RECODE_BIMBAM_1CHR;
	      misc_flags |= MISC_SET_HH_MISSING;
	      fputs("Note: --recode-bimbam flag deprecated.  Use \"--recode bimbam\" or\n\"--recode bimbam-1chr\".\n", stdout);
	      ujj = 2;
	    }
	    break;
	  case 17:
	    if (!memcmp(argptr2, "-fastphase", 10)) {
	      memcpy(flagptr, "recode 01 fastphase-1chr", 25);
	      recode_modifier |= RECODE_01 | RECODE_FASTPHASE_1CHR;
	      misc_flags |= MISC_SET_HH_MISSING;
	      fputs("Note: --recode-fastphase flag deprecated.  Use e.g. \"--recode 01 fastphase-1chr\".\n", stdout);
	      ujj = 2;
	    } else if (!memcmp(argptr2, "-structure", 10)) {
	      memcpy(flagptr, "recode structure", 17);
	      recode_modifier |= RECODE_STRUCTURE;
	      misc_flags |= MISC_SET_HH_MISSING;
	      ujj = 1;
	    }
	    break;
	  }
	  if (ujj) {
	    if (ujj == 1) {
	      printf("Note: --%s flag deprecated.  Use \"--%s ...\".\n", argptr, flagptr);
	    }
	    break;
	  }
	} else if (!strcmp(argptr, "reference-allele")) {
	  memcpy(flagptr, "a1-allele", 10);
	  break;
	}
	goto main_flag_copy;
      case 't':
        if (!strcmp(argptr, "thread-num")) {
	  memcpy(flagptr, "threads", 8);
	  break;
	}
	goto main_flag_copy;
      case 'u':
	if (!strcmp(argptr, "update-freq")) {
	  memcpy(flagptr, "read-freq", 10);
	  break;
	} else if (!strcmp(argptr, "update-ref-allele")) {
	  // GCTA alias
	  memcpy(flagptr, "a1-allele", 10);
	  break;
	}
	goto main_flag_copy;
      default:
      main_flag_copy:
	memcpy(flagptr, argptr, ukk);
      }
      flagptr = &(flagptr[MAX_FLAG_LEN]);
      flag_map[cur_flag++] = uii;
    }
  }
  // requires MAX_FLAG_LEN to be at least sizeof(void*) + sizeof(int32_t)
  sptr = (char*)malloc(flag_ct * MAX_FLAG_LEN);
  if (!sptr) {
    goto main_ret_NOMEM_NOLOG2;
  }
  qsort_ext2(flag_buf, flag_ct, MAX_FLAG_LEN, strcmp_deref, (char*)flag_map, sizeof(int32_t), sptr, MAX_FLAG_LEN);
  free(sptr);
  ujj = strlen_se(flag_buf);
  for (cur_flag = 1; cur_flag < flag_ct; cur_flag++) {
    ukk = strlen_se(&(flag_buf[cur_flag * MAX_FLAG_LEN]));
    if ((ujj == ukk) && (!memcmp(&(flag_buf[(cur_flag - 1) * MAX_FLAG_LEN]), &(flag_buf[cur_flag * MAX_FLAG_LEN]), ukk))) {
      flag_buf[cur_flag * MAX_FLAG_LEN + ukk] = '\0'; // just in case of aliases
      fprintf(stderr, "Error: Duplicate --%s flag.\n", &(flag_buf[cur_flag * MAX_FLAG_LEN]));
      goto main_ret_INVALID_CMDLINE;
    }
    ujj = ukk;
  }

  uii = 5;
  memcpy(outname, PROG_NAME_STR, 6);
  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    ii = memcmp("out", &(flag_buf[cur_flag * MAX_FLAG_LEN]), 4);
    if (!ii) {
      ujj = flag_map[cur_flag];
      ukk = param_count(argc, argv, ujj);
      if (enforce_param_ct_range(ukk, argv[ujj], 1, 1)) {
	fputs(g_logbuf, stderr);
	fputs(errstr_append, stderr);
	goto main_ret_INVALID_CMDLINE;
      }
      if (strlen(argv[ujj + 1]) > (FNAMESIZE - MAX_POST_EXT)) {
	fputs("Error: --out parameter too long.\n", stderr);
	goto main_ret_OPEN_FAIL;
      }
      uii = strlen(argv[ujj + 1]);
      memcpy(outname, argv[ujj + 1], uii + 1);
      outname_end = &(outname[uii]);
    }
    if (ii <= 0) {
      break;
    }
  }
  memcpy(&(outname[uii]), ".log", 5);
  g_logfile = fopen(outname, "w");
  if (!g_logfile) {
    fprintf(stderr, "Error: Failed to open %s.  Try ", outname);
    if (!memcmp(outname, PROG_NAME_STR, 6)) {
      fputs("using --out.\n", stderr);
    } else {
      fputs("changing the --out parameter.\n", stderr);
    }
    goto main_ret_OPEN_FAIL;
  }
  printf("Logging to %s.\n", outname);
  outname[uii] = '\0';

  logstr(ver_str);
  logstr("\n");
  logprint("Options in effect:\n");
  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    logprint("  --");
    logprint(&(flag_buf[cur_flag * MAX_FLAG_LEN]));
    ii = flag_map[cur_flag] + 1;
    while ((ii < argc) && (!is_flag(argv[ii]))) {
      logprint(" ");
      logprint(argv[ii++]);
    }
    logprint("\n");
  }
  logprint("\n");

#ifdef _WIN32
  windows_dw = TEXTBUF_SIZE;
  if (GetComputerName(g_textbuf, &windows_dw))
#else
  if (gethostname(g_textbuf, TEXTBUF_SIZE) != -1)
#endif
  {
    logstr("Hostname: ");
    logstr(g_textbuf);
  }
  logstr("\nWorking directory: ");
  if (!getcwd(g_textbuf, FNAMESIZE)) {
    goto main_ret_READ_FAIL;
  }
  logstr(g_textbuf);
  logstr("\nStart time: ");
  time(&rawtime);
  logstr(ctime(&rawtime));
  logstr("\n");

#ifdef _WIN32
  GetSystemInfo(&sysinfo);
  g_thread_ct = sysinfo.dwNumberOfProcessors;
  known_procs = g_thread_ct;
#else
  known_procs = sysconf(_SC_NPROCESSORS_ONLN);
  g_thread_ct = (known_procs == -1)? 1 : known_procs;
#endif
  if (g_thread_ct > 8) {
    if (g_thread_ct > MAX_THREADS) {
      g_thread_ct = MAX_THREADS;
    } else {
      g_thread_ct--;
    }
  }
  mapname[0] = '\0';
  pedname[0] = '\0';
  famname[0] = '\0';
  memcpyl3(output_missing_pheno, "-9");
  // stuff that must be processed before regular alphabetical loop
  retval = init_delim_and_species(flag_ct, flag_buf, flag_map, argc, argv, &range_delim, &chrom_info);
  if (retval) {
    goto main_ret_1;
  }
  fill_ulong_zero(CHROM_MASK_INITIAL_WORDS, chrom_exclude);
  cur_flag = 0;
  do {
    argptr = &(flag_buf[cur_flag * MAX_FLAG_LEN]);
    if (!(*argptr)) {
      // preprocessed
      continue;
    }
    argptr2 = &(argptr[1]);
    cur_arg = flag_map[cur_flag];
    param_ct = param_count(argc, argv, cur_arg);
    switch (*argptr) {
    case '1':
      if (*argptr2 == '\0') {
	misc_flags |= MISC_AFFECTION_01;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case '2':
      if (!memcmp(argptr2, "3file", 6)) {
	if (chrom_info.species != SPECIES_HUMAN) {
	  logerrprint("Error: --23file cannot be used with a nonhuman species flag.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 7)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	ii = strlen(argv[cur_arg + 1]);
	if (ii > FNAMESIZE - 1) {
	  logerrprint("Error: --23file filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
        memcpy(pedname, argv[cur_arg + 1], ii + 1);
	if (param_ct > 1) {
	  if (strchr(argv[cur_arg + 2], ' ')) {
	    logerrprint("Error: Space present in --23file family ID.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (alloc_string(&fid_23, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	  if (param_ct > 2) {
	    if (strchr(argv[cur_arg + 3], ' ')) {
	      logerrprint("Error: Space present in --23file sample ID.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (alloc_string(&iid_23, argv[cur_arg + 3])) {
	      goto main_ret_NOMEM;
	    }
	    if (param_ct > 3) {
	      cc = extract_char_param(argv[cur_arg + 4]);
	      if ((cc == 'M') || (cc == 'm') || (cc == '1')) {
		modifier_23 |= M23_MALE;
	      } else if ((cc == 'F') || (cc == 'f') || (cc == '2')) {
		modifier_23 |= M23_FEMALE;
	      } else if (cc == '0') {
		modifier_23 |= M23_FORCE_MISSING_SEX;
	      } else if ((cc != 'I') && (cc != 'i')) {
		logerrprint("Error: Invalid --23file sex parameter (M or 1 = male, F or 2 = female,\nI = infer from data, 0 = force missing).\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      if (param_ct > 4) {
		if (scan_doublex(argv[cur_arg + 5], &pheno_23)) {
		  sprintf(g_logbuf, "Error: Invalid --23file phenotype '%s'.\n", argv[cur_arg + 5]);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
		if (param_ct > 5) {
		  if (strchr(argv[cur_arg + 6], ' ')) {
		    logerrprint("Error: Space present in --23file paternal ID.\n");
		    goto main_ret_INVALID_CMDLINE;
		  }
		  if (alloc_string(&paternal_id_23, argv[cur_arg + 6])) {
		    goto main_ret_NOMEM;
		  }
		  if (param_ct > 6) {
		    if (strchr(argv[cur_arg + 7], ' ')) {
		      logerrprint("Error: Space present in --23file maternal ID.\n");
		      goto main_ret_INVALID_CMDLINE;
		    }
		    if (alloc_string(&maternal_id_23, argv[cur_arg + 7])) {
		      goto main_ret_NOMEM;
		    }
		  }
		}
	      }
	    }
	  }
	}
	load_rare = LOAD_RARE_23;
      } else if ((!memcmp(argptr2, "3file-convert-xy", 17)) || (!memcmp(argptr2, "3file-make-xylist", 18))) {
        sprintf(g_logbuf, "Error: --%s has been retired due to brain-damaged design.  Use\n--split-x instead.\n", argptr);
        goto main_ret_INVALID_CMDLINE_2A;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'D':
      if (*argptr2 == '\0') {
	logprint("Note: --D flag deprecated.  Use e.g. \"--r2 dprime\".\n");
	ld_info.modifier |= LD_DPRIME;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'K':
      if (*argptr2 == '\0') {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &cluster.min_ct)) {
          sprintf(g_logbuf, "Error: Invalid --K cluster count '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'R':
      if (*argptr2 == '\0') {
#if defined __cplusplus && !defined _WIN32
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp(argv[cur_arg + 1], "debug")) {
	    uii = 2;
	  } else if (strcmp(argv[cur_arg + 2], "debug")) {
	    sprintf(g_logbuf, "Error: Invalid --R modifier '%s'.\n", argv[cur_arg + 2]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
          misc_flags |= MISC_RPLUGIN_DEBUG;
	}
	retval = alloc_fname(&rplugin_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	calculation_type |= CALC_RPLUGIN;
#else
  #ifdef _WIN32
        logerrprint("Error: --R does not currently support Windows.\n");
  #else
	logerrprint("Error: --R requires " PROG_NAME_CAPS " to be built with a C++, not a pure C, compiler.\n");
  #endif
	goto main_ret_INVALID_CMDLINE;
#endif

      } else if (!memcmp(argptr2, "-port", 6)) {
	if (!rplugin_fname) {
	  logerrprint("Error: --R-port must be used with --R.\n");
	  goto main_ret_INVALID_CMDLINE;
        }
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_cappedx(argv[cur_arg + 1], 65535, (uint32_t*)(&rplugin_port))) {
	  sprintf(g_logbuf, "Error: Invalid --R-port parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "-debug", 7)) {
	if (!rplugin_fname) {
	  logerrprint("Error: --R-debug must be used with --R.\n");
	  goto main_ret_INVALID_CMDLINE;
        }
	logprint("Note: --R-debug flag deprecated.  Use e.g. \"--R <filename> debug\".\n");
	misc_flags |= MISC_RPLUGIN_DEBUG;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "-host", 6)) {
        if (!rplugin_fname) {
          logerrprint("Error: --R-host must be used with --R.\n");
          goto main_ret_INVALID_CMDLINE;
        }
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_string(&rplugin_host_or_socket, argv[cur_arg + 1]);
        if (retval) {
          goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "-socket", 8)) {
	if (!rplugin_fname) {
	  logerrprint("Error: --R-socket must be used with --R.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (rplugin_port != -2) {
	  logerrprint("Error: --R-socket cannot be used with --R-port.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (rplugin_host_or_socket) {
	  logerrprint("Error: --R-socket cannot be used with --R-host.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_string(&rplugin_host_or_socket, argv[cur_arg + 1]);
        if (retval) {
          goto main_ret_INVALID_CMDLINE_WWA;
	}
	rplugin_port = -1;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'Z':
      if (!memcmp(argptr2, "-genome", 8)) {
	kk = 1;
	genome_modifier |= GENOME_OUTPUT_GZ;
	goto main_genome_flag;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }

    case 'a':
      if (!memcmp(argptr2, "utosome", 8)) {
	fill_bits(1, chrom_info.autosome_ct, chrom_info.chrom_mask);
	chrom_info.is_include_stack = 1;
	chrom_flag_present = 1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "utosome-xy", 11)) {
	if (chrom_flag_present) {
          logerrprint("Error: --autosome-xy cannot be used with --autosome.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (chrom_info.xymt_codes[XY_OFFSET] == -2) {
	  logerrprint("Error: --autosome-xy used with a species lacking an XY region.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	fill_bits(1, chrom_info.autosome_ct, chrom_info.chrom_mask);
	set_bit(chrom_info.xymt_codes[XY_OFFSET], chrom_info.chrom_mask);
	chrom_info.is_include_stack = 1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "llow-extra-chr", 15)) {
	if (load_rare == LOAD_RARE_23) {
	  logerrprint("Error: --allow-extra-chr cannot currently be used with --23file.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (param_ct) {
	  if (memcmp("0", argv[cur_arg + 1], 2)) {
            sprintf(g_logbuf, "Error: Invalid --allow-extra-chr parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  chrom_info.zero_extra_chroms = 1;
	}
	misc_flags |= MISC_ALLOW_EXTRA_CHROMS;
      } else if (!memcmp(argptr2, "llow-no-sex", 12)) {
        sex_missing_pheno |= ALLOW_NO_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "llow-no-samples", 16)) {
	misc_flags |= MISC_ALLOW_NO_SAMPLES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "llow-no-vars", 13)) {
	misc_flags |= MISC_ALLOW_NO_VARS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "llow-no-covars", 15)) {
	covar_modifier |= COVAR_ALLOW_NONE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ll", 3)) {
	logprint("Note: --all flag has no effect.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "llele1234", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct == 1) {
	  if (strcmp("multichar", argv[cur_arg + 1])) {
	    sprintf(g_logbuf, "Error: Invalid --allele1234 parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  allelexxxx = ALLELE_RECODE_MULTICHAR;
	} else {
	  allelexxxx = ALLELE_RECODE;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "lleleACGT", 9)) {
	if (allelexxxx) {
	  logerrprint("Error: --allele1234 and --alleleACGT cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct == 1) {
	  if (strcmp("multichar", argv[cur_arg + 1])) {
	    sprintf(g_logbuf, "Error: Invalid --alleleACGT parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  allelexxxx = ALLELE_RECODE_ACGT | ALLELE_RECODE_MULTICHAR;
	} else {
	  allelexxxx = ALLELE_RECODE_ACGT;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "llele-count", 12)) {
	lgen_modifier |= LGEN_ALLELE_COUNT;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ll-pheno", 9)) {
	pheno_modifier |= PHENO_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ssoc", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 6)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "counts")) {
	    if (model_modifier & MODEL_QMASK) {
	      logerrprint("Error: --assoc 'qt-means' modifier does not make sense with 'counts'.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_ASSOC_COUNTS;
	  } else if (!strcmp(argv[cur_arg + uii], "fisher")) {
	    if (model_modifier & MODEL_QMASK) {
	      logerrprint("Error: --assoc 'qt-means'/'lin' does not make sense with 'fisher'.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_FISHER;
	  } else if (!strcmp(argv[cur_arg + uii], "fisher-midp")) {
            if (model_modifier & MODEL_QMASK) {
              logerrprint("Error: --assoc 'qt-means'/'lin' does not make sense with 'fisher-midp'.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
            model_modifier |= MODEL_FISHER | MODEL_FISHER_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (model_modifier & MODEL_MPERM) {
	      logerrprint("Error: --assoc 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_PERM;
	  } else if (!strcmp(argv[cur_arg + uii], "genedrop")) {
	    if (model_modifier & MODEL_QMASK) {
	      logerrprint("Error: --assoc 'qt-means'/'lin' does not make sense with 'genedrop'.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    model_modifier |= MODEL_PERM_COUNT;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (model_modifier & MODEL_PERM) {
	      logerrprint("Error: --assoc 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (model_modifier & MODEL_MPERM) {
	      logerrprint("Error: Duplicate --assoc 'mperm' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (scan_posint_defcapx(&(argv[cur_arg + uii][6]), &model_mperm_val)) {
	      sprintf(g_logbuf, "Error: Invalid --assoc mperm parameter '%s'.\n", &(argv[cur_arg + uii][6]));
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    model_modifier |= MODEL_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "qt-means")) {
	    if (model_modifier & MODEL_DMASK) {
	      logerrprint("Error: --assoc 'qt-means' does not make sense with a case/control-specific\nmodifier.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_QT_MEANS;
	  } else if (!strcmp(argv[cur_arg + uii], "lin")) {
	    if (model_modifier & MODEL_DMASK) {
	      logerrprint("Error: --assoc 'lin' does not make sense with a case/control-specific modifier.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_LIN;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    logerrprint("Error: Improper --assoc mperm syntax.  (Use \"--assoc mperm=<value>\".)\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
	    model_modifier |= MODEL_SET_TEST;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --assoc parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	model_modifier |= MODEL_ASSOC;
	calculation_type |= CALC_MODEL;
      } else if (!memcmp(argptr2, "djust", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	mtest_adjust = 1;
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "gc")) {
	    mtest_adjust |= ADJUST_GC;
	  } else if (!strcmp(argv[cur_arg + uii], "log10")) {
	    mtest_adjust |= ADJUST_LOG10;
	  } else if (!strcmp(argv[cur_arg + uii], "qq-plot")) {
	    mtest_adjust |= ADJUST_QQ;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --adjust parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
      } else if (!memcmp(argptr2, "perm", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 6)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &aperm.min)) {
	  sprintf(g_logbuf, "Error: Invalid --aperm min permutation count '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	aperm.min++;
	if (param_ct > 1) {
	  if (scan_posint_cappedx(argv[cur_arg + 2], APERM_MAX, &aperm.max)) {
	    sprintf(g_logbuf, "Error: Invalid --aperm max permutation count '%s'.\n", argv[cur_arg + 2]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	// may as well disallow equality since there's no reason not to use
	// max(T) then...
	if (aperm.min >= aperm.max) {
	  logerrprint("Error: --aperm min permutation count must be smaller than max.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (param_ct > 2) {
	  if (scan_doublex(argv[cur_arg + 3], &aperm.alpha)) {
	    sprintf(g_logbuf, "Error: Invalid --aperm alpha threshold '%s'.\n", argv[cur_arg + 3]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (param_ct > 3) {
	    if (scan_doublex(argv[cur_arg + 4], &aperm.beta) || (aperm.beta <= 0)) {
	      sprintf(g_logbuf, "Error: Invalid --aperm beta '%s'.\n", argv[cur_arg + 4]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (param_ct > 4) {
	      if (scan_doublex(argv[cur_arg + 5], &aperm.init_interval)) {
		sprintf(g_logbuf, "Error: Invalid --aperm initial pruning interval '%s'.\n", argv[cur_arg + 5]);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      if ((aperm.init_interval < 1) || (aperm.init_interval > 1000000)) {
		sprintf(g_logbuf, "Error: Invalid --aperm initial pruning interval '%s'.\n", argv[cur_arg + 5]);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      if (param_ct == 6) {
		if (scan_doublex(argv[cur_arg + 6], &aperm.interval_slope)) {
		  sprintf(g_logbuf, "Error: Invalid --aperm pruning interval slope '%s'.\n", argv[cur_arg + 6]);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
		if ((aperm.interval_slope < 0) || (aperm.interval_slope > 1)) {
		  sprintf(g_logbuf, "Error: Invalid --aperm pruning interval slope '%s'.\n", argv[cur_arg + 6]);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
	      }
	    }
	  }
	}
	aperm_present = 1;
      } else if (!memcmp(argptr2, "1-allele", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_2col(&a1alleles, &(argv[cur_arg + 1]), argptr, param_ct);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "2-allele", 9)) {
	if (a1alleles) {
	  logerrprint("Error: --a2-allele cannot be used with --a1-allele.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_2col(&a2alleles, &(argv[cur_arg + 1]), argptr, param_ct);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "nnotate", 8)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 10)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&annot_info.fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  ujj = strlen(argv[cur_arg + uii]);
	  if ((ujj > 7) && (!memcmp(argv[cur_arg + uii], "attrib=", 7))) {
            retval = alloc_fname(&annot_info.attrib_fname, &(argv[cur_arg + uii][7]), argptr, 0);
	    if (retval) {
	      goto main_ret_1;
	    }
	  } else if ((ujj > 7) && (!memcmp(argv[cur_arg + uii], "ranges=", 7))) {
            retval = alloc_fname(&annot_info.ranges_fname, &(argv[cur_arg + uii][7]), argptr, 0);
	    if (retval) {
	      goto main_ret_1;
	    }
	  } else if ((ujj > 7) && (!memcmp(argv[cur_arg + uii], "filter=", 7))) {
            retval = alloc_fname(&annot_info.filter_fname, &(argv[cur_arg + uii][7]), argptr, 0);
	    if (retval) {
	      goto main_ret_1;
	    }
	  } else if ((ujj > 7) && (!memcmp(argv[cur_arg + uii], "subset=", 7))) {
            retval = alloc_fname(&annot_info.subset_fname, &(argv[cur_arg + uii][7]), argptr, 0);
	    if (retval) {
	      goto main_ret_1;
	    }
	  } else if ((ujj > 5) && (!memcmp(argv[cur_arg + uii], "snps=", 5))) {
            retval = alloc_fname(&annot_info.snps_fname, &(argv[cur_arg + uii][5]), argptr, 0);
	    if (retval) {
	      goto main_ret_1;
	    }
	  } else if ((ujj == 2) && (!memcmp(argv[cur_arg + uii], "NA", 2))) {
	    if (annot_info.modifier & ANNOT_PRUNE) {
              logerrprint("Error: --annotate 'NA' and 'prune' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    annot_info.modifier |= ANNOT_NA;
	  } else if ((ujj == 5) && (!memcmp(argv[cur_arg + uii], "prune", 5))) {
	    if (annot_info.modifier & ANNOT_NA) {
              logerrprint("Error: --annotate 'NA' and 'prune' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    annot_info.modifier |= ANNOT_PRUNE;
	  } else if ((ujj == 5) && (!memcmp(argv[cur_arg + uii], "block", 5))) {
            annot_info.modifier |= ANNOT_BLOCK;
	  } else if ((ujj == 7) && (!memcmp(argv[cur_arg + uii], "minimal", 7))) {
	    annot_info.modifier |= ANNOT_MINIMAL;
	  } else if ((ujj == 8) && (!memcmp(argv[cur_arg + uii], "distance", 8))) {
	    annot_info.modifier |= ANNOT_DISTANCE;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --annotate parameter '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	if ((annot_info.modifier & ANNOT_BLOCK) && (annot_info.modifier & (ANNOT_NA | ANNOT_MINIMAL))) {
	  logerrprint("Error: --annotate 'block' cannot be used with 'NA' or 'minimal'.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (!annot_info.attrib_fname) {
	  if (!annot_info.ranges_fname) {
	    logerrprint("Error: --annotate must be used with 'attrib' and/or 'ranges'.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	} else if (!annot_info.ranges_fname) {
	  if (annot_info.subset_fname) {
	    logerrprint("Error: --annotate 'subset' modifier must be used with 'ranges'.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  } else if (annot_info.modifier & (ANNOT_MINIMAL | ANNOT_DISTANCE)) {
	    sprintf(g_logbuf, "Error: --annotate '%s' modifier must be used with 'ranges'.\n", (annot_info.modifier & ANNOT_MINIMAL)? "minimal" : "distance");
            goto main_ret_INVALID_CMDLINE_2A;
	  }
	}
      } else if (!memcmp(argptr2, "nnotate-snp-field", 18)) {
	if (!annot_info.attrib_fname) {
	  logerrprint("Error: --annotate-snp-field must be used with --annotate + 'attrib'.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&annot_info.snpfield, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ttrib", 6)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&filter_attrib_fname, argv[cur_arg + 1], argptr, 0);
        if (retval) {
          goto main_ret_1;
	}
	if (param_ct == 2) {
	  // force comma-terminated string to simplify parsing
	  uii = strlen(argv[cur_arg + 2]);
	  filter_attrib_liststr = (char*)malloc(uii + 2);
	  if (!filter_attrib_liststr) {
	    goto main_ret_NOMEM;
	  }
          memcpy(filter_attrib_liststr, argv[cur_arg + 2], uii);
	  memcpy(&(filter_attrib_liststr[uii]), ",", 2);
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "ttrib-indiv", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&filter_attrib_sample_fname, argv[cur_arg + 1], argptr, 0);
        if (retval) {
          goto main_ret_1;
	}
	if (param_ct == 2) {
	  uii = strlen(argv[cur_arg + 2]);
	  filter_attrib_sample_liststr = (char*)malloc(uii + 2);
	  if (!filter_attrib_sample_liststr) {
	    goto main_ret_NOMEM;
	  }
          memcpy(filter_attrib_sample_liststr, argv[cur_arg + 2], uii);
	  memcpy(&(filter_attrib_sample_liststr[uii]), ",", 2);
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if ((!memcmp(argptr2, "lt-group", 9)) ||
                 (!memcmp(argptr2, "lt-snp", 7))) {
        goto main_hap_disabled_message;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'b':
      if (!memcmp(argptr2, "file", 5)) {
	if (load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  csptr = argv[cur_arg + 1];
	  if (strlen(csptr) > (FNAMESIZE - 5)) {
	    logerrprint("Error: --bfile parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  csptr = PROG_NAME_STR;
	}
	if (!(load_params & LOAD_PARAMS_BED)) {
	  memcpy(strcpya(pedname, csptr), ".bed", 5);
	  load_params |= LOAD_PARAMS_BED;
	}
	memcpy(strcpya(mapname, csptr), ".bim", 5);
	memcpy(strcpya(famname, csptr), ".fam", 5);
	load_params |= LOAD_PARAMS_BIM | LOAD_PARAMS_FAM;
      } else if (!memcmp(argptr2, "ed", 3)) {
	if (load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	load_params |= LOAD_PARAMS_BED;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logerrprint("Error: --bed parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "im", 3)) {
	if (load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	load_params |= LOAD_PARAMS_BIM;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logerrprint("Error: --bim parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "merge", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct == 2) {
	  logerrprint("Error: --bmerge must have exactly 1 or 3 parameters.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (param_ct == 3) {
	  if (++jj > FNAMESIZE) {
	    logerrprint("Error: --bmerge .bed filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  jj = strlen(argv[cur_arg + 2]) + 1;
	  if (jj > FNAMESIZE) {
	    logerrprint("Error: --bmerge .bim filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename2, argv[cur_arg + 2], jj);
	  jj = strlen(argv[cur_arg + 3]) + 1;
	  if (jj > FNAMESIZE) {
	    logerrprint("Error: --bmerge .fam filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename3, argv[cur_arg + 3], jj);
	} else {
	  if (jj > (FNAMESIZE - 5)) {
	    logerrprint("Error: --bmerge filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(memcpya(mergename1, argv[cur_arg + 1], jj), ".bed", 5);
	  memcpy(memcpya(mergename2, argv[cur_arg + 1], jj), ".bim", 5);
	  memcpy(memcpya(mergename3, argv[cur_arg + 1], jj), ".fam", 5);
	}
	calculation_type |= CALC_MERGE;
	merge_type |= MERGE_BINARY;
      } else if (!memcmp(argptr2, "p-space", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &min_bp_space)) {
	  sprintf(g_logbuf, "Error: Invalid --bp-space minimum bp distance '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "eta", 4)) {
	logprint("Note: --beta flag deprecated.  Use e.g. \"--logistic beta\".\n");
	glm_modifier |= GLM_BETA;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "d", 2)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (cluster.modifier & CLUSTER_CMH_MPERM) {
	      logerrprint("Error: --bd 'mperm' and 'perm[-bd]' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (cluster.modifier & CLUSTER_CMH_PERM_BD) {
	      logerrprint("Error: --bd 'perm' and 'perm-bd' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    cluster.modifier |= CLUSTER_CMH_PERM;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-bd")) {
	    if (cluster.modifier & CLUSTER_CMH_MPERM) {
	      logerrprint("Error: --bd 'mperm' and 'perm[-bd]' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if ((cluster.modifier & (CLUSTER_CMH_PERM | CLUSTER_CMH_PERM_BD)) == CLUSTER_CMH_PERM) {
	      logerrprint("Error: --bd 'perm' and 'perm-bd' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (mtest_adjust) {
	      logerrprint("Error: --bd 'perm-bd' mode cannot currently be used with --adjust.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    cluster.modifier |= CLUSTER_CMH_PERM_BD;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (cluster.modifier & CLUSTER_CMH_PERM) {
	      logerrprint("Error: --bd 'mperm' and 'perm[-bd]' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (cluster.modifier & CLUSTER_CMH_MPERM) {
	      logerrprint("Error: Duplicate --bd 'mperm' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (scan_posint_defcapx(&(argv[cur_arg + uii][6]), &(cluster.cmh_mperm_val))) {
	      sprintf(g_logbuf, "Error: Invalid --bd mperm parameter '%s'.\n", &(argv[cur_arg + uii][6]));
              goto main_ret_INVALID_CMDLINE_WWA;
	    }
            cluster.modifier |= CLUSTER_CMH_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
            cluster.modifier |= CLUSTER_CMH_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
	    cluster.modifier |= CLUSTER_CMH_SET_TEST;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
            logerrprint("Error: Improper --bd mperm syntax.  (Use \"--bd mperm=<value>\".)\n");
            goto main_ret_INVALID_CMDLINE_A;
	  } else {
            sprintf(g_logbuf, "Error: Invalid --bd parameter '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_CMH;
	cluster.modifier |= CLUSTER_CMH_BD;
      } else if (!memcmp(argptr2, "iallelic-only", 14)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "strict")) {
	    misc_flags |= MISC_BIALLELIC_ONLY_STRICT;
	  } else if (!strcmp(argv[cur_arg + uii], "list")) {
	    misc_flags |= MISC_BIALLELIC_ONLY_LIST;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --biallelic-only modifier '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        misc_flags |= MISC_BIALLELIC_ONLY;
      } else if (!memcmp(argptr2, "cf", 3)) {
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = strlen(argv[cur_arg + 1]);
	if (uii > FNAMESIZE - 1) {
	  logerrprint("Error: --bcf filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(pedname, argv[cur_arg + 1], uii + 1);
	load_rare = LOAD_RARE_BCF;
      } else if (!memcmp(argptr2, "gen", 4)) {
	if (load_rare || (load_params & LOAD_PARAMS_BFILE_ALL)) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	load_params |= LOAD_PARAMS_OXBGEN;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logerrprint("Error: --bgen parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
	if (param_ct == 2) {
	  if (strcmp(argv[cur_arg + 2], "snpid-chr")) {
	    sprintf(g_logbuf, "Error: Invalid --bgen modifier '%s'.\n", argv[cur_arg + 2]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
          misc_flags |= MISC_OXFORD_SNPID_CHR;
	}
      } else if (!memcmp(argptr2, "locks", 6)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "no-pheno-req")) {
            ld_info.modifier |= LD_BLOCKS_NO_PHENO_REQ;
	  } else if (!strcmp(argv[cur_arg + uii], "no-small-max-span")) {
            ld_info.modifier |= LD_BLOCKS_NO_SMALL_MAX_SPAN;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --blocks parameter '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_BLOCKS;
      } else if (!memcmp(argptr2,"locks-inform-frac", 17)) {
	if (!(calculation_type & CALC_BLOCKS)) {
	  logerrprint("Error: --blocks-inform-frac must be used with --blocks.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0) || (dxx > 1.0 - SMALL_EPSILON)) {
	  sprintf(g_logbuf, "Error: Invalid --blocks-inform-frac parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	ld_info.blocks_inform_frac = dxx;
      } else if (!memcmp(argptr2, "locks-max-kb", 13)) {
	if (!(calculation_type & CALC_BLOCKS)) {
	  logerrprint("Error: --blocks-max-kb must be used with --blocks.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --blocks-max-kb parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (dxx > 2147483.646) {
	  ld_info.blocks_max_bp = 2147483646;
	} else {
	  ld_info.blocks_max_bp = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
	}
      } else if (!memcmp(argptr2, "locks-min-maf", 14)) {
        if (!(calculation_type & CALC_BLOCKS)) {
	  logerrprint("Error: --blocks-min-maf must be used with --blocks.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0) || (dxx > 0.5)) {
	  sprintf(g_logbuf, "Error: Invalid --blocks-min-maf parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        ld_info.blocks_min_maf = dxx;
      } else if (!memcmp(argptr2, "locks-recomb-highci", 20)) {
        if (!(calculation_type & CALC_BLOCKS)) {
	  logerrprint("Error: --blocks-recomb-highci must be used with --blocks.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0) || (dxx > 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --blocks-recomb-highci parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        ld_info.blocks_recomb_highci = ((int32_t)((dxx + SMALL_EPSILON) * 100));
	if (ld_info.blocks_recomb_highci < 2) {
	  logerrprint("Error: --blocks-recomb-highci parameter must be at least 0.02.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	ld_info.blocks_recomb_highci--;
      } else if (!memcmp(argptr2, "locks-strong-highci", 20)) {
        if (!(calculation_type & CALC_BLOCKS)) {
	  logerrprint("Error: --blocks-strong-highci must be used with --blocks.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < SMALL_EPSILON) || (dxx > 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --blocks-strong-highci parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        ld_info.blocks_strong_highci = (int32_t)((dxx - SMALL_EPSILON) * 100);
	// 0.83 lower bound is needed for now for sane handling of size-2
	// special case
	if (ld_info.blocks_strong_highci < 83) {
	  logerrprint("Error: --blocks-strong-highci parameter currently must be larger than 0.83.\nContact the developers if this is problematic.\n");
          goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "locks-strong-lowci", 19)) {
        if (!(calculation_type & CALC_BLOCKS)) {
	  logerrprint("Error: --blocks-strong-lowci must be used with --blocks.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < SMALL_EPSILON) || (dxx >= 1)) {
	  sprintf(g_logbuf, "Error: Invalid --blocks-strong-lowci parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	ld_info.blocks_strong_lowci_outer = 2 + (int32_t)((dxx - SMALL_EPSILON) * 100);
	ld_info.blocks_strong_lowci = 2 + (int32_t)((dxx + SMALL_EPSILON) * 100);
	if ((ld_info.blocks_strong_lowci_outer < 52) || (ld_info.blocks_strong_lowci > 82)) {
	  logerrprint("Error: --blocks-strong-lowci parameter currently must be in (0.5, 0.81).\nContact the developers if this is problematic.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "order", 6)) {
	if (!annot_info.ranges_fname) {
	  logerrprint("Error: --border now must be used with --annotate + 'ranges'.\n");
	  if (!annot_info.fname) {
            logerrprint("Use --make-set-border with --make-set, etc.)\n");
	  }
          goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --border parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WW;
	}
	if (dxx > 2147483.646) {
	  annot_info.border = 0x7ffffffe;
	} else {
	  annot_info.border = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "pfile", 6)) {
	sprintf(g_logbuf, "Error: Unrecognized flag (%s).  (This is PLINK 1.9, not 2.x.)\n", argv[cur_arg]);
	goto main_ret_INVALID_CMDLINE_2;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'c':
      if (!memcmp(argptr2, "hr", 3)) {
	if (chrom_flag_present) {
	  logerrprint("Error: --chr cannot be used with --autosome[-xy].\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        retval = parse_chrom_ranges(param_ct, '-', &(argv[cur_arg]), chrom_info.chrom_mask, &chrom_info, (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1, argptr);
	if (retval) {
	  goto main_ret_1;
	}
	chrom_info.is_include_stack = 1;
	chrom_flag_present = 1;
      } else if (!memcmp(argptr2, "ompound-genotypes", 18)) {
	logprint("Note: --compound-genotypes flag unnecessary (spaces between alleles in .ped\nand .lgen files are optional if all alleles are single-character).\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ompress", 8)) {
	logerrprint("Error: --compress flag retired.  Use e.g. \"gzip <filename>\".\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ounts", 6)) {
	if (model_modifier & MODEL_ASSOC) {
	  if (model_modifier & MODEL_QMASK) {
	    logerrprint("Error: --assoc 'qt-means'/'lin' does not make sense with --counts.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  logprint("Note: --counts flag deprecated.  Use \"--assoc counts\" instead.\n");
          model_modifier |= MODEL_ASSOC_COUNTS;
	} else {
	  logprint("Note: --counts flag deprecated.  Use \"--freq counts\" or --freqx instead.\n");
	}
	misc_flags |= MISC_FREQ_COUNTS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ovar", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp(argv[cur_arg + 1], "keep-pheno-on-missing-cov")) {
	    uii = 2;
	  } else if (strcmp(argv[cur_arg + 2], "keep-pheno-on-missing-cov")) {
	    sprintf(g_logbuf, "Error: Invalid --covar parameter '%s'.\n", argv[cur_arg + 2]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (rplugin_fname) {
	    logerrprint("Error: --covar 'keep-pheno-on-missing-cov' modifier cannot be used with --R.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
          covar_modifier |= COVAR_KEEP_PHENO_ON_MISSING_COV;
	}
	retval = alloc_fname(&covar_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ovar-name", 10)) {
	if (!covar_fname) {
	  logerrprint("Error: --covar-name must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &covar_range_list, 0);
	if (retval) {
	  goto main_ret_1;
	}
	covar_modifier |= COVAR_NAME;
      } else if (!memcmp(argptr2, "ovar-number", 12)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (covar_modifier & COVAR_NAME) {
	  logerrprint("Error: --covar-number cannot be used with --covar-name.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (!covar_fname) {
	  logerrprint("Error: --covar-number must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = parse_name_ranges(param_ct, '-', &(argv[cur_arg]), &covar_range_list, 1);
	if (retval) {
	  goto main_ret_1;
	}
	covar_modifier |= COVAR_NUMBER;
      } else if (!memcmp(argptr2, "ell", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], (uint32_t*)&model_cell_ct)) {
	  sprintf(g_logbuf, "Error: Invalid --cell parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "i", 2)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx)) {
	  sprintf(g_logbuf, "Error: Invalid --ci parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((dxx < 0.01) || (dxx >= 1.0)) {
	  logerrprint("Error: --ci confidence interval size s must satisfy 0.01 <= s < 1.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	ci_size = dxx;
      } else if (!memcmp(argptr2, "luster", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "cc")) {
            cluster.modifier |= CLUSTER_CC;
	  } else if (!strcmp(argv[cur_arg + uii], "group-avg")) {
	    if (cluster.modifier & CLUSTER_OLD_TIEBREAKS) {
              logerrprint("Error: --cluster 'group-avg' and 'old-tiebreaks' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    cluster.modifier |= CLUSTER_GROUP_AVG;
	  } else if (!strcmp(argv[cur_arg + uii], "missing")) {
	    cluster.modifier |= CLUSTER_MISSING;
	  } else if (!strcmp(argv[cur_arg + uii], "only2")) {
	    cluster.modifier |= CLUSTER_ONLY2;
	  } else if (!strcmp(argv[cur_arg + uii], "old-tiebreaks")) {
	    if (cluster.modifier & CLUSTER_GROUP_AVG) {
              logerrprint("Error: --cluster 'group-avg' and 'old-tiebreaks' cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
	    cluster.modifier |= CLUSTER_OLD_TIEBREAKS;
	  } else {
            sprintf(g_logbuf, "Error: Invalid --cluster parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        calculation_type |= CALC_CLUSTER;
      } else if (!memcmp(argptr2, "c", 2)) {
        logprint("Note: --cc flag deprecated.  Use \"--cluster cc\".\n");
        cluster.modifier |= CLUSTER_CC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "luster-missing", 15)) {
	if (calculation_type & CALC_CLUSTER) {
	  logerrprint("Error: --cluster-missing cannot be used with --cluster.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --cluster-missing flag deprecated.  Use \"--cluster missing\".\n");
        calculation_type |= CALC_CLUSTER;
        cluster.modifier |= CLUSTER_MISSING;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "file", 5)) {
        UNSTABLE("cfile");
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (calculation_type & CALC_MERGE) {
	  logerrprint("Error: --cfile cannot be used with --bmerge.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	sptr = argv[cur_arg + 1];
	uii = strlen(sptr);
	if (uii > (FNAMESIZE - 9)) {
	  logerrprint("Error: --cfile parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(memcpya(pedname, sptr, uii), ".cnv", 5);
	memcpy(memcpya(famname, sptr, uii), ".fam", 5);
	memcpy(memcpya(mapname, sptr, uii), ".cnv.map", 9);
	load_rare = LOAD_RARE_CNV;
      } else if (!memcmp(argptr2, "nv-count", 9)) {
	UNSTABLE("cnv-count");
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&cnv_intersect_filter_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	cnv_intersect_filter_type = CNV_COUNT;
      } else if (!memcmp(argptr2, "nv-del", 7)) {
	UNSTABLE("cnv-del");
	cnv_calc_type |= CNV_DEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-disrupt", 11)) {
	UNSTABLE("cnv-disrupt");
	cnv_overlap_type = CNV_DISRUPT;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-dup", 7)) {
	UNSTABLE("cnv-dup");
	if (cnv_calc_type & CNV_DEL) {
	  logerrprint("Error: --cnv-dup cannot be used with --cnv-del.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	cnv_calc_type |= CNV_DUP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-enrichment-test", 19)) {
	UNSTABLE("cnv-enrichment-test");
	if (!cnv_intersect_filter_type) {
	  logerrprint("Error: --cnv-enrichment-test must be used with --cnv-count.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_posint_defcapx(argv[cur_arg + 1], &cnv_enrichment_test_mperms)) {
	    sprintf(g_logbuf, "Error: Invalid --cnv-enrichment-test permutation count '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	cnv_calc_type |= CNV_ENRICHMENT_TEST;
      } else if (!memcmp(argptr2, "nv-exclude", 11)) {
	UNSTABLE("cnv-exclude");
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (cnv_intersect_filter_type) {
	  logerrprint("Error: --cnv-exclude cannot be used with --cnv-count.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	retval = alloc_fname(&cnv_intersect_filter_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	cnv_intersect_filter_type = CNV_EXCLUDE;
      } else if (!memcmp(argptr2, "nv-exclude-off-by-1", 20)) {
	UNSTABLE("cnv-exclude-off-by-1");
        cnv_calc_type |= CNV_EXCLUDE_OFF_BY_1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-freq-exclude-above", 22)) {
	UNSTABLE("cnv-freq-exclude-above");
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &cnv_freq_val)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-freq-exclude-above parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	cnv_freq_type = CNV_FREQ_EXCLUDE_ABOVE;
      } else if (!memcmp(argptr2, "nv-freq-exclude-below", 22)) {
	UNSTABLE("cnv-freq-exclude-below");
	if (cnv_freq_type) {
	  logerrprint("Error: --cnv-freq-exclude-below cannot be used with --cnv-freq-exclude-above.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &cnv_freq_val) || (cnv_freq_val == 1)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-freq-exclude-below parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	cnv_freq_type = CNV_FREQ_EXCLUDE_BELOW;
      } else if (!memcmp(argptr2, "nv-freq-exclude-exact", 22)) {
	UNSTABLE("cnv-freq-exclude-exact");
	if (cnv_freq_type) {
	  logerrprint("Error: --cnv-freq-exclude-exact cannot be used with\n--cnv-freq-exclude-above/-below.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &cnv_freq_val)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-freq-exclude-exact parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	cnv_freq_type = CNV_FREQ_EXCLUDE_EXACT;
      } else if (!memcmp(argptr2, "nv-freq-include-exact", 22)) {
	UNSTABLE("cnv-freq-include-exact");
	if (cnv_freq_type) {
	  logerrprint("Error: --cnv-freq-include-exact cannot be used with\n--cnv-freq-exclude-above/-below/-exact.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &cnv_freq_val)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-freq-include-exact parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	cnv_freq_type = CNV_FREQ_INCLUDE_EXACT;
      } else if (!memcmp(argptr2, "nv-freq-method2", 16)) {
	UNSTABLE("cnv-freq-method2");
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_doublex(argv[cur_arg + 1], &cnv_freq_val2) || (cnv_freq_val2 < 0) || (cnv_freq_val2 > 1)) {
	    sprintf(g_logbuf, "Error: Invalid --cnv-freq-method2 parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	cnv_freq_type |= CNV_FREQ_METHOD2;
	if (cnv_freq_val2 == 0) {
	  // allow >= comparison to be used
	  cnv_freq_val2 = SMALLISH_EPSILON;
	}
      } else if (!memcmp(argptr2, "nv-freq-overlap", 16)) {
	UNSTABLE("cnv-freq-overlap");
	if (!(cnv_freq_type & CNV_FREQ_FILTER)) {
	  logerrprint("Error: --cnv-freq-overlap must be used with --cnv-freq-include-exact or\n--cnv-freq-exclude-above/-below/-exact.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_freq_type & CNV_FREQ_METHOD2) {
	  logerrprint("Error: --cnv-freq-overlap cannot be used with --cnv-freq-method2.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_doublex(argv[cur_arg + 1], &cnv_freq_val2) || (cnv_freq_val2 < 0) || (cnv_freq_val2 > 1)) {
	    sprintf(g_logbuf, "Error: Invalid --cnv-freq-overlap parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	if (cnv_freq_val2 == 0) {
	  cnv_freq_val2 = SMALLISH_EPSILON;
	}
	cnv_freq_type |= CNV_FREQ_OVERLAP;
      } else if (!memcmp(argptr2, "nv-indiv-perm", 14)) {
	UNSTABLE("cnv-indiv-perm");
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_posint_defcapx(argv[cur_arg + 1], &cnv_sample_mperms)) {
	    sprintf(g_logbuf, "Error: Invalid --cnv-indiv-perm permutation count '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	cnv_calc_type |= CNV_SAMPLE_PERM;
      } else if (!memcmp(argptr2, "nv-intersect", 13)) {
	UNSTABLE("cnv-intersect");
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (cnv_intersect_filter_type) {
	  logerrprint("Error: --cnv-intersect cannot be used with --cnv-count/--cnv-exclude.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	retval = alloc_fname(&cnv_intersect_filter_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	cnv_intersect_filter_type = CNV_INTERSECT;
      } else if (!memcmp(argptr2, "nv-kb", 6)) {
	UNSTABLE("cnv-kb");
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0.001) || (dxx > 2147483.646)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-kb size '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	cnv_min_seglen = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
      } else if (!memcmp(argptr2, "nv-list", 8)) {
	UNSTABLE("cnv-list");
	if ((load_rare & (~LOAD_RARE_CNV)) || load_params) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (calculation_type & CALC_MERGE) {
	  logerrprint("Error: --cnv-list cannot be used with --bmerge.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	strcpya(pedname, argv[cur_arg + 1]);
	load_rare = LOAD_RARE_CNV;
      } else if (!memcmp(argptr2, "nv-make-map", 12)) {
	UNSTABLE("cnv-make-map");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-make-map cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "short")) {
            sprintf(g_logbuf, "Error: Invalid --cnv-make-map parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  cnv_calc_type |= CNV_MAKE_MAP;
	} else {
	  cnv_calc_type |= CNV_MAKE_MAP | CNV_MAKE_MAP_LONG;
	}
      } else if (!memcmp(argptr2, "nv-max-kb", 10)) {
	UNSTABLE("cnv-max-kb");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-max-kb cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0.001) || (dxx > 2147483.646)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-max-kb size '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	cnv_max_seglen = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
	if (cnv_min_seglen > cnv_max_seglen) {
	  logerrprint("Error: --cnv-max-kb value cannot be smaller than --cnv-kb value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "nv-max-score", 13)) {
	UNSTABLE("cnv-max-score");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-max-score cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &cnv_max_score)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-max-score value '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "nv-max-sites", 13)) {
	UNSTABLE("cnv-max-sites");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-max-sites cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &cnv_max_sites)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-max-sites parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "nv-overlap", 11)) {
	UNSTABLE("cnv-overlap");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-overlap cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_overlap_type == CNV_DISRUPT) {
	  logerrprint("Error: --cnv-overlap cannot be used with --cnv-disrupt.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &cnv_overlap_val) || (cnv_overlap_val < 0) || (cnv_overlap_val > 1))  {
	  sprintf(g_logbuf, "Error: Invalid --cnv-overlap value '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (cnv_overlap_val != 0) {
	  // ignore --cnv-overlap 0
	  cnv_overlap_type = CNV_OVERLAP;
	}
	if ((cnv_freq_type & CNV_FREQ_FILTER) && (!(cnv_freq_type & (CNV_FREQ_OVERLAP | CNV_FREQ_METHOD2)))) {
	  logprint("Note: --cnv-overlap + --cnv-freq-... deprecated.  Use --cnv-freq-overlap.\n");
	  if (cnv_overlap_val != 0) {
	    cnv_freq_type |= CNV_FREQ_OVERLAP;
	    cnv_freq_val2 = cnv_overlap_val;
	  }
	}
      } else if (!memcmp(argptr2, "nv-region-overlap", 18)) {
	UNSTABLE("cnv-region-overlap");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-region-overlap cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_overlap_type) {
	  logerrprint("Error: --cnv-region-overlap cannot be used with --cnv-overlap/-disrupt.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &cnv_overlap_val) || (cnv_overlap_val <= 0) || (cnv_overlap_val > 1))  {
	  sprintf(g_logbuf, "Error: Invalid --cnv-region-overlap value '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	cnv_overlap_type = CNV_OVERLAP_REGION;
      } else if (!memcmp(argptr2, "nv-score", 9)) {
	UNSTABLE("cnv-score");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-score cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &cnv_min_score)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-score value '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (cnv_min_score > cnv_max_score) {
	  logerrprint("Error: --cnv-score value cannot be greater than --cnv-max-score value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "nv-sites", 9)) {
	UNSTABLE("cnv-sites");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-sites cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &cnv_min_sites)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-sites parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (cnv_min_sites > cnv_max_sites) {
	  logerrprint("Error: --cnv-sites value cannot be greater than --cnv-max-sites value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "nv-subset", 10)) {
	UNSTABLE("cnv-subset");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-subset cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!cnv_intersect_filter_type) {
	  logerrprint("Error: --cnv-subset must be used with --cnv-intersect/-exclude/-count.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&cnv_subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "nv-test", 8)) {
	UNSTABLE("cnv-test");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-test cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct == 2) {
	  if (!strcmp(argv[cur_arg + 1], "1sided")) {
	    uii = 2;
	    cnv_calc_type |= CNV_TEST_FORCE_1SIDED;
	  } else if (!strcmp(argv[cur_arg + 1], "2sided")) {
	    uii = 2;
	    cnv_calc_type |= CNV_TEST_FORCE_2SIDED;
	  } else {
	    uii = 1;
	    if (!strcmp(argv[cur_arg + 2], "1sided")) {
	      cnv_calc_type |= CNV_TEST_FORCE_1SIDED;
	    } else if (!strcmp(argv[cur_arg + 2], "2sided")) {
	      cnv_calc_type |= CNV_TEST_FORCE_2SIDED;
	    } else {
	      logerrprint("Error: Invalid --cnv-test parameter sequence.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	} else {
	  uii = 1;
	}
	if (scan_posint_defcapx(argv[cur_arg + uii], &cnv_test_mperms)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-test permutation count '%s'.\n", argv[cur_arg + uii]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	cnv_calc_type |= CNV_TEST;
      } else if (!memcmp(argptr2, "nv-test-1sided", 15)) {
	UNSTABLE("cnv-test-1sided");
	if (cnv_calc_type & CNV_TEST_FORCE_2SIDED) {
	  logerrprint("Error: --cnv-test cannot be both 1-sided and 2-sided at the same time.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --cnv-test-1sided flag deprecated.  Use \"--cnv-test 1sided\".\n");
	cnv_calc_type |= CNV_TEST_FORCE_1SIDED;
      } else if (!memcmp(argptr2, "nv-test-2sided", 15)) {
	UNSTABLE("cnv-test-2sided");
	if (cnv_calc_type & CNV_TEST_FORCE_1SIDED) {
	  logerrprint("Error: --cnv-test cannot be both 1-sided and 2-sided at the same time.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --cnv-test-2sided flag deprecated.  Use \"--cnv-test 2sided\".\n");
	cnv_calc_type |= CNV_TEST_FORCE_2SIDED;
      } else if (!memcmp(argptr2, "nv-test-region", 15)) {
	UNSTABLE("cnv-test-region");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-test-region cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_posint_defcapx(argv[cur_arg + 1], &cnv_test_region_mperms)) {
	    sprintf(g_logbuf, "Error: Invalid --cnv-test-region permutation count '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	cnv_calc_type |= CNV_TEST_REGION;
      } else if (!memcmp(argptr2, "nv-test-window", 15)) {
	UNSTABLE("cnv-test-window");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-test-window cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0.001)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-test-window size '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	dxx *= 1000;
	if (dxx > 2147483646) {
	  cnv_test_window = 0x7ffffffe;
	} else {
	  cnv_test_window = (int32_t)(dxx * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "nv-union-overlap", 17)) {
	UNSTABLE("cnv-union-overlap");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-union-overlap cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_overlap_type) {
	  logerrprint("Error: --cnv-union-overlap cannot be used with --cnv-{region-}overlap/-disrupt.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &cnv_overlap_val) || (cnv_overlap_val <= 0) || (cnv_overlap_val > 1)) {
	  sprintf(g_logbuf, "Error: Invalid --cnv-union-overlap value '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	cnv_overlap_type = CNV_OVERLAP_UNION;
      } else if (!memcmp(argptr2, "nv-write", 9)) {
	UNSTABLE("cnv-write");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-write cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "freq")) {
            sprintf(g_logbuf, "Error: Invalid --cnv-write parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (!(cnv_freq_val & CNV_FREQ_METHOD2)) {
	    logerrprint("Error: --cnv-write 'freq' modifier must be used with --cnv-freq-method2.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  cnv_calc_type |= CNV_WRITE_FREQ;
	}
	cnv_calc_type |= CNV_WRITE;
      } else if (!memcmp(argptr2, "nv-write-freq", 14)) {
	UNSTABLE("cnv-write-freq");
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logerrprint("Error: --cnv-write freq cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!(cnv_freq_val & CNV_FREQ_METHOD2)) {
	  logerrprint("Error: --cnv-write 'freq' modifier must be used with --cnv-freq-method2.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (!(cnv_calc_type & CNV_WRITE)) {
	  logerrprint("Error: --cnv-write-freq must be used with --cnv-write.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --cnv-write-freq flag deprecated.  Use \"--cnv-write freq\".\n");
	cnv_calc_type |= CNV_WRITE_FREQ;
      } else if (!memcmp(argptr2, "onsensus-match", 15)) {
        logprint("Note: --consensus-match flag deprecated.  Use \"--homozyg consensus-match\".\n");
	homozyg.modifier |= HOMOZYG_CONSENSUS_MATCH;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ondition", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp("dominant", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_DOMINANT;
	  } else if (!strcmp("recessive", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_RECESSIVE;
	  } else {
	    uii = 2;
            if (!strcmp("dominant", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_DOMINANT;
	    } else if (!strcmp("recessive", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_RECESSIVE;
	    } else {
	      logerrprint("Error: Invalid --condition parameter sequence.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	}
	if (alloc_string(&condition_mname, argv[cur_arg + uii])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ondition-list", 14)) {
	if (condition_mname) {
	  logerrprint("Error: --condition-list cannot be used with --condition.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp("dominant", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_DOMINANT;
	  } else if (!strcmp("recessive", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_RECESSIVE;
	  } else {
	    uii = 2;
            if (!strcmp("dominant", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_DOMINANT;
	    } else if (!strcmp("recessive", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_RECESSIVE;
	    } else {
	      logerrprint("Error: Invalid --condition-list parameter sequence.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  }
	}
	retval = alloc_fname(&condition_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "omplement-sets", 15)) {
        set_info.modifier |= SET_COMPLEMENTS | SET_C_PREFIX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "onst-fid", 9)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&const_fid, param_ct? argv[cur_arg + 1] : "0")) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ase-only", 9)) {
	logprint("Note: --case-only flag deprecated.  Use \"--fast-epistasis case-only\".\n");
        epi_info.modifier |= EPI_FAST_CASE_ONLY;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "m-map", 6)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct == 1) {
	  // must contain exactly one '@'
          sptr = strchr(argv[cur_arg + 1], '@');
	  if (!sptr) {
            logerrprint("Error: --cm-map requires either a '@' in the filename pattern, or a chromosome\ncode as the second parameter.\n");
            goto main_ret_INVALID_CMDLINE_A;
	  }
          if (strchr(&(sptr[1]), '@')) {
	    logerrprint("Error: Multiple '@'s in --cm-map filename pattern.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
          if (alloc_string(&cm_map_fname, argv[cur_arg + 1])) {
	    goto main_ret_NOMEM;
	  }
	} else {
          retval = alloc_fname(&cm_map_fname, argv[cur_arg + 1], argptr, 0);
          if (retval) {
            goto main_ret_1;
	  }
          if (alloc_string(&cm_map_chrname, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "heck-sex", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 5)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	// pre-scan for 'y-only', since that changes interpretation of the
	// first two numeric parameters
	ujj = 0;
	if (param_ct) {
	  for (uii = 1; uii <= param_ct; uii++) {
            if (!strcmp(argv[cur_arg + uii], "y-only")) {
	      ujj = uii;
	      break;
	    }
	  }
	}
	if (!ujj) {
	  // now ujj is number of numeric parameters read
	  for (uii = 1; uii <= param_ct; uii++) {
	    if (!strcmp(argv[cur_arg + uii], "ycount")) {
	      misc_flags |= MISC_SEXCHECK_YCOUNT;
	    } else {
	      if (!ujj) {
		if (scan_doublex(argv[cur_arg + uii], &check_sex_fthresh) || (check_sex_fthresh <= 0.0)) {
		  sprintf(g_logbuf, "Error: Invalid --check-sex female F-statistic estimate ceiling '%s'.\n", argv[cur_arg + uii]);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
	      } else if (ujj == 1) {
		if (scan_doublex(argv[cur_arg + uii], &check_sex_mthresh) || (check_sex_mthresh >= 1.0)) {
		  sprintf(g_logbuf, "Error: Invalid --check-sex male F-statistic estimate floor '%s'.\n", argv[cur_arg + uii]);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
	      } else if (ujj == 2) {
		if (scan_uint_defcapx(argv[cur_arg + uii], &check_sex_f_yobs)) {
		  logerrprint("Error: Invalid --check-sex female Ychr maximum nonmissing genotype count.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
	      } else if (ujj == 3) {
		if (scan_uint_defcapx(argv[cur_arg + uii], &check_sex_m_yobs)) {
		  logerrprint("Error: Invalid --check-sex male Ychr minimum nonmissing genotype count.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
	      } else {
		logerrprint("Error: Invalid --check-sex parameter sequence.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      ujj++;
	    }
	  }
	  if ((ujj > 2) && (!(misc_flags & MISC_SEXCHECK_YCOUNT))) {
	    logerrprint("Error: --check-sex only accepts >2 numeric parameters in 'ycount' mode.\n");
            goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (check_sex_fthresh > check_sex_mthresh) {
	    logerrprint("Error: --check-sex female F estimate ceiling cannot be larger than male floor.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  // actually fine if check_sex_f_yobs > check_sex_m_yobs
	} else {
	  check_sex_m_yobs = 1;
	  ukk = 0; // number of numeric parameters
	  // in practice, second numeric parameter is not really optional if
	  // first is provided...
	  for (uii = 1; uii <= param_ct; uii++) {
	    if (uii == ujj) {
	      continue;
	    }
	    // may as well print a more informative error message in this case
	    if (!strcmp(argv[cur_arg + uii], "ycount")) {
	      logerrprint("Error: Conflicting --check-sex modes.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (!ukk) {
	      if (scan_uint_defcapx(argv[cur_arg + uii], &check_sex_f_yobs)) {
		logerrprint("Error: Invalid --check-sex y-only female Ychr maximum nonmissing genotype\ncount.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else if (ukk == 1) {
	      if (scan_uint_defcapx(argv[cur_arg + uii], &check_sex_m_yobs)) {
		logerrprint("Error: Invalid --check-sex y-only male Ychr minimum nonmissing genotype count.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else {
	      logerrprint("Error: Invalid --check-sex y-only parameter sequence.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    ukk++;
	  }
	  if (check_sex_f_yobs >= check_sex_m_yobs) {
	    logerrprint("Error: In y-only mode, --check-sex female Y observation threshold must be\nsmaller than male threshold.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  misc_flags |= MISC_SEXCHECK_YCOUNT | MISC_SEXCHECK_YONLY;
	}
        calculation_type |= CALC_SEXCHECK;
      } else if (!memcmp(argptr2, "lump", 5)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_and_flatten_comma_delim(&clump_info.fnames_flattened, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
	calculation_type |= CALC_CLUMP;
	// only needs to distinguish between one and two
	if ((param_ct > 1) || strchr(argv[cur_arg + 1], ',')) {
	  clump_info.fname_ct = 2;
	} else {
	  clump_info.fname_ct = 1;
	}
      } else if (!memcmp(argptr2, "lump-allow-overlap", 19)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-allow-overlap must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
        clump_info.modifier |= CLUMP_ALLOW_OVERLAP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "lump-annotate", 14)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-annotate must be used with --clump.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten_comma_delim(&clump_info.annotate_flattened, &(argv[cur_arg + 1]), param_ct);
        if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "lump-best", 10)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-best must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
        clump_info.modifier |= CLUMP_BEST;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "lump-field", 11)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-field must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&clump_info.pfield_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "lump-index-first", 17)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-index-first must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE_A;
	} else if (clump_info.fname_ct == 1) {
	  logerrprint("Warning: --clump-index-first has no effect when there is only one --clump file.\n");
	} else {
          clump_info.modifier |= CLUMP_INDEX_FIRST;
	}
        goto main_param_zero;
      } else if (!memcmp(argptr2, "lump-kb", 8)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-kb must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0.001)) {
	  sprintf(g_logbuf, "Error: Invalid --clump-kb parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WW;
	}
	dxx *= 1000;
        if (dxx > 2147483647.0) {
	  clump_info.bp_radius = 0x7ffffffe;
	} else {
	  clump_info.bp_radius = ((int32_t)(dxx * (1 + SMALL_EPSILON) - 1));
	}
      } else if (!memcmp(argptr2, "lump-p1", 8)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-p1 must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx <= 0) || (dxx > 1)) {
	  sprintf(g_logbuf, "Error: Invalid --clump-p1 parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WW;
	}
	clump_info.p1 = dxx;
      } else if (!memcmp(argptr2, "lump-p2", 8)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-p2 must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < clump_info.p1) || (dxx > 1)) {
	  sprintf(g_logbuf, "Error: Invalid --clump-p2 parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WW;
	}
	clump_info.p2 = dxx;
      } else if (!memcmp(argptr2, "lump-r2", 8)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-r2 must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx >= 1)) {
	  sprintf(g_logbuf, "Error: Invalid --clump-r2 parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WW;
	}
	clump_info.r2 = dxx;
      } else if (!memcmp(argptr2, "lump-range", 11)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-range must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&clump_info.range_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "lump-range-border", 18)) {
	if (!clump_info.range_fname) {
	  logerrprint("Error: --clump-range-border must be used with --clump-range.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --clump-range-border parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WW;
	}
	if (dxx > 2147483.646) {
	  clump_info.range_border = 0x7ffffffe;
	} else {
	  clump_info.range_border = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "lump-replicate", 15)) {
        if (clump_info.fname_ct < 2) {
	  logerrprint("Error: --clump-replicate requires multiple --clump files.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	if (clump_info.modifier & CLUMP_BEST) {
	  logprint("Note: --clump-replicate no longer has any effect when used with --clump-best.\n");
	} else {
          clump_info.modifier |= CLUMP_REPLICATE;
	}
        goto main_param_zero;
      } else if (!memcmp(argptr2, "lump-snp-field", 15)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-snp-field must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&clump_info.snpfield_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "lump-verbose", 13)) {
        if (!clump_info.fname_ct) {
	  logerrprint("Error: --clump-verbose must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
        clump_info.modifier |= CLUMP_VERBOSE;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "hap", 4)) {
        goto main_hap_disabled_message;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'd':
      if (!memcmp(argptr2, "ebug", 5)) {
	g_debug_on = 1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ata", 4)) {
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  csptr = argv[cur_arg + 1];
	  if (strlen(csptr) > (FNAMESIZE - 8)) {
	    logerrprint("Error: --data parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  csptr = PROG_NAME_STR;
	}
	if (!(load_params & LOAD_PARAMS_OXBGEN)) {
	  memcpy(strcpya(pedname, csptr), ".gen", 5);
	  load_params |= LOAD_PARAMS_OXGEN;
	}
	// cheating: this is of course more like a .fam file
	memcpy(strcpya(mapname, csptr), ".sample", 8);
	load_params |= LOAD_PARAMS_OXSAMPLE;
      } else if (!memcmp(argptr2, "ecompress", 10)) {
	logerrprint("Error: --decompress flag retired.  Use e.g. \"gunzip <filename>\".\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "istance", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 6)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "square")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      logerrprint("Error: --distance 'square' and 'square0' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      logerrprint("Error: --distance 'square' and 'triangle' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    dist_calc_type |= DISTANCE_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      logerrprint("Error: --distance 'square' and 'square0' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      logerrprint("Error: --distance 'square0' and 'triangle' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    dist_calc_type |= DISTANCE_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      logerrprint("Error: --distance 'square' and 'triangle' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      logerrprint("Error: --distance 'square0' and 'triangle' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    dist_calc_type |= DISTANCE_TRI;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (dist_calc_type & (DISTANCE_BIN | DISTANCE_BIN4)) {
	      logerrprint("Error: Conflicting --distance modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    dist_calc_type |= DISTANCE_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (dist_calc_type & (DISTANCE_GZ | DISTANCE_BIN4)) {
	      logerrprint("Error: Conflicting --distance modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    dist_calc_type |= DISTANCE_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "bin4")) {
	    if (dist_calc_type & (DISTANCE_GZ | DISTANCE_BIN)) {
	      logerrprint("Error: Conflicting --distance modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    dist_calc_type |= DISTANCE_BIN4;
	  } else if (!strcmp(argv[cur_arg + uii], "ibs")) {
	    if (dist_calc_type & DISTANCE_IBS) {
	      logerrprint("Error: Duplicate --distance 'ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "1-ibs")) {
	    if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	      logerrprint("Error: Duplicate --distance '1-ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_1_MINUS_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "allele-ct")) {
	    if (dist_calc_type & DISTANCE_ALCT) {
	      logerrprint("Error: Duplicate --distance 'allele-ct' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_ALCT;
	  } else if (!strcmp(argv[cur_arg + uii], "flat-missing")) {
	    dist_calc_type |= DISTANCE_FLAT_MISSING;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --distance parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	if (!(dist_calc_type & DISTANCE_TYPEMASK)) {
	  dist_calc_type |= DISTANCE_ALCT;
	}
	calculation_type |= CALC_DISTANCE;
      } else if (!memcmp(argptr2, "istance-exp", 12)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &distance_exp)) {
	  sprintf(g_logbuf, "Error: Invalid --distance-exp parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WW;
	}
	fputs("Note: \"--distance-exp <x>\" deprecated.  Use \"--distance-wts exp=<x>\" instead.\n", stdout);
      } else if (!memcmp(argptr2, "istance-wts", 12)) {
	if (distance_exp != 0.0) {
	  logerrprint("Error: --distance-wts cannot be used with --distance-exp.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (calculation_type & CALC_PLINK1_DISTANCE_MATRIX) {
	  logerrprint("Error: --distance-wts cannot be used with --distance-matrix.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
        }
	if ((strlen(argv[cur_arg + 1]) > 4) && (!memcmp(argv[cur_arg + 1], "exp=", 4))) {
	  if (scan_doublex(&(argv[cur_arg + 1][4]), &distance_exp)) {
	    sprintf(g_logbuf, "Error: Invalid --distance-wts exponent '%s'.\n", &(argv[cur_arg + 1][4]));
	    goto main_ret_INVALID_CMDLINE_WW;
	  }
	} else {
	  UNSTABLE("distance-wts");
	  uii = 1;
	  if (param_ct == 2) {
	    if (!strcmp(argv[cur_arg + 1], "noheader")) {
	      uii = 2;
	    } else if (strcmp(argv[cur_arg + 2], "noheader")) {
	      sprintf(g_logbuf, "Error: Invalid --distance-wts parameter '%s'.\n", argv[cur_arg + 2]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    dist_calc_type |= DISTANCE_WTS_NOHEADER;
	  }
	  retval = alloc_fname(&distance_wts_fname, argv[cur_arg + uii], argptr, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "istance-matrix", 15)) {
	if (distance_exp != 0.0) {
	  logerrprint("Error: --distance-matrix cannot be used with --distance-exp.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	  logerrprint("Error: --distance-matrix flag cannot be used with \"--distance 1-ibs\".\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	calculation_type |= CALC_PLINK1_DISTANCE_MATRIX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ummy", 5)) {
	if (load_params) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 6)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &dummy_sample_ct) || ((!dummy_sample_ct) && (!(misc_flags & MISC_ALLOW_NO_SAMPLES)))) {
	  logerrprint("Error: Invalid --dummy sample count.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 2], &dummy_marker_ct) || ((!dummy_marker_ct) && (!(misc_flags & MISC_ALLOW_NO_VARS)))) {
	  logerrprint("Error: Invalid --dummy variant count.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        for (uii = 3; uii <= param_ct; uii++) {
	  if (match_upper(argv[cur_arg + uii], "ACGT")) {
	    if (dummy_flags & (DUMMY_1234 | DUMMY_12)) {
	      logerrprint("Error: --dummy 'acgt' modifier cannot be used with '1234' or '12'.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
            dummy_flags |= DUMMY_ACGT;
	  } else if (!strcmp(argv[cur_arg + uii], "1234")) {
	    if (dummy_flags & (DUMMY_ACGT | DUMMY_12)) {
	      logerrprint("Error: --dummy '1234' modifier cannot be used with 'acgt' or '12'.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
            dummy_flags |= DUMMY_1234;
	  } else if (!strcmp(argv[cur_arg + uii], "12")) {
	    if (dummy_flags & (DUMMY_ACGT | DUMMY_1234)) {
	      logerrprint("Error: --dummy '12' modifier cannot be used with 'acgt' or '1234'.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
            dummy_flags |= DUMMY_12;
	  } else if (!strcmp(argv[cur_arg + uii], "scalar-pheno")) {
	    dummy_flags |= DUMMY_SCALAR_PHENO;
	  } else {
	    if ((dummy_flags & DUMMY_MISSING_PHENO) || scan_doublex(argv[cur_arg + uii], &dxx) || (dxx < 0.0) || (dxx > 1.0)) {
	      sprintf(g_logbuf, "Error: Invalid --dummy parameter '%s'.\n", argv[cur_arg + uii]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    } else if (dummy_flags & DUMMY_MISSING_GENO) {
	      dummy_missing_pheno = dxx;
	      dummy_flags |= DUMMY_MISSING_PHENO;
	    } else {
	      dummy_missing_geno = dxx;
	      dummy_flags |= DUMMY_MISSING_GENO;
	    }
	  }
	}
	load_rare = LOAD_RARE_DUMMY;
      } else if (!memcmp(argptr2, "ummy-coding", 12)) {
	if (!covar_fname) {
	  logerrprint("Error: --dummy-coding cannot be used without --covar.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (!strcmp(argv[cur_arg + 1], "no-round")) {
	    uii = 2;
	    write_covar_modifier |= WRITE_COVAR_DUMMY_NO_ROUND;
	  } else {
	    if (param_ct == 2) {
              if (strcmp(argv[cur_arg + 2], "no-round")) {
		logerrprint("Error: Invalid --dummy-coding parameter sequence.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      write_covar_modifier |= WRITE_COVAR_DUMMY_NO_ROUND;
	    }
	    uii = 1;
	  }
	  if (uii <= param_ct) {
	    if (scan_posint_defcapx(argv[cur_arg + uii], &write_covar_dummy_max_categories) || (write_covar_dummy_max_categories < 3)) {
	      sprintf(g_logbuf, "Error: Invalid --dummy-coding max categories parameter '%s'.\n", argv[cur_arg + uii]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	}
	write_covar_modifier |= WRITE_COVAR_DUMMY;
      } else if (!memcmp(argptr2, "ominant", 8)) {
	logprint("Note: --dominant flag deprecated.  Use e.g. \"--linear dominant\" (and\n\"--condition-list <filename> dominant\" to change covariate coding).\n");
	glm_modifier |= GLM_DOMINANT | GLM_CONDITION_DOMINANT;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ouble-id", 9)) {
        if (const_fid) {
	  logerrprint("Error: --double-id cannot be used with --const-fid.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        misc_flags |= MISC_DOUBLE_ID;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "osage", 6)) {
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	} else if (pheno_modifier & PHENO_ALL) {
	  logerrprint("Error: --dosage cannot be used with --all-pheno.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (mtest_adjust) {
	  logerrprint("Error: --dosage cannot be used with --adjust.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (covar_modifier & COVAR_KEEP_PHENO_ON_MISSING_COV) {
	  logerrprint("Error: --dosage cannot be used with --covar 'keep-pheno-on-missing-cov'\nmodifier.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (condition_mname || condition_fname) {
	  logerrprint("Error: --dosage does not support --condition/--condition-list.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (misc_flags & (MISC_ALLOW_NO_SAMPLES | MISC_ALLOW_NO_VARS)) {
	  logerrprint("Error: --dosage does not support --allow-no-samples/--allow-no-vars.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 13)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&dosage_info.fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	for (uii = 2; uii <= param_ct; uii++) {
          if (!strcmp(argv[cur_arg + uii], "list")) {
	    dosage_info.modifier |= DOSAGE_LIST;
	  } else if (!strcmp(argv[cur_arg + uii], "dose1")) {
	    dosage_info.modifier |= DOSAGE_DOSE1;
	  } else if (!strcmp(argv[cur_arg + uii], "Z")) {
	    logprint("Note: --dosage 'Z' modifier is now equivalent to 'Zout', since gzipped input is\nautomatically supported.\n");
	    dosage_info.modifier |= DOSAGE_ZOUT;
	  } else if (!strcmp(argv[cur_arg + uii], "Zin")) {
	    logprint("Note: --dosage 'Zin' modifier now has no effect, since gzipped input is\nautomatically supported.\n");
	  } else if (!strcmp(argv[cur_arg + uii], "Zout")) {
	    dosage_info.modifier |= DOSAGE_ZOUT;
	  } else if (!strcmp(argv[cur_arg + uii], "occur")) {
	    dosage_info.modifier |= DOSAGE_OCCUR;
	  } else if (!strcmp(argv[cur_arg + uii], "sex")) {
	    dosage_info.modifier |= DOSAGE_SEX;
	  } else if (!strcmp(argv[cur_arg + uii], "case-control-freqs")) {
	    dosage_info.modifier |= DOSAGE_FREQ_CC;
	  } else if (!strcmp(argv[cur_arg + uii], "frq2")) {
	    // turn this into an error before official 1.90 release
	    logerrprint("Warning: The --dosage 'frq2' modifier has been renamed to 'case-control-freqs'.\n");
	    dosage_info.modifier |= DOSAGE_FREQ_CC;
	  } else if (strlen(argv[cur_arg + uii]) <= 6) {
	    goto main_dosage_invalid_param;
	  } else if (!strcmp(argv[cur_arg + uii], "sepheader")) {
	    if (dosage_info.modifier & DOSAGE_NOHEADER) {
	      logerrprint("Error: --dosage 'sepheader' and 'noheader' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    dosage_info.modifier |= DOSAGE_SEPHEADER;
	  } else if (!strcmp(argv[cur_arg + uii], "noheader")) {
	    if (dosage_info.modifier & DOSAGE_SEPHEADER) {
	      logerrprint("Error: --dosage 'sepheader' and 'noheader' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    dosage_info.modifier |= DOSAGE_NOHEADER;
	  } else if (!memcmp(argv[cur_arg + uii], "skip0=", 6)) {
	    if (scan_uint_defcapx(&(argv[cur_arg + uii][6]), &(dosage_info.skip0))) {
	      sprintf(g_logbuf, "Error: Invalid --dosage skip0 parameter '%s'.\n", &(argv[cur_arg + uii][6]));
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  } else if (!memcmp(argv[cur_arg + uii], "skip1=", 6)) {
	    if (scan_uint_defcapx(&(argv[cur_arg + uii][6]), &(dosage_info.skip1))) {
	      sprintf(g_logbuf, "Error: Invalid --dosage skip1 parameter '%s'.\n", &(argv[cur_arg + uii][6]));
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  } else if (!memcmp(argv[cur_arg + uii], "skip2=", 6)) {
	    if (scan_uint_defcapx(&(argv[cur_arg + uii][6]), &(dosage_info.skip2))) {
	      sprintf(g_logbuf, "Error: Invalid --dosage skip2 parameter '%s'.\n", &(argv[cur_arg + uii][6]));
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  } else if (!memcmp(argv[cur_arg + uii], "format=", 7)) {
	    ujj = ((unsigned char)argv[cur_arg + uii][7]) - '1';
	    if ((ujj > 2) || argv[cur_arg + uii][8]) {
	      sprintf(g_logbuf, "Error: Invalid --dosage format parameter '%s'.\n", &(argv[cur_arg + uii][7]));
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    dosage_info.format = ujj + 1;
	  } else if (!strcmp(argv[cur_arg + uii], "standard-beta")) {
	    glm_modifier |= GLM_STANDARD_BETA;
	  } else {
	  main_dosage_invalid_param:
	    sprintf(g_logbuf, "Error: Invalid --dosage modifier '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	if (dosage_info.modifier & DOSAGE_OCCUR) {
	  if ((glm_modifier & GLM_STANDARD_BETA) || (dosage_info.modifier & DOSAGE_SEX)) {
	    logerrprint("Error: --dosage 'occur' mode cannot be used with association analysis\nmodifiers/flags.\n");
            goto main_ret_INVALID_CMDLINE_A;
	  }
	} else {
	  dosage_info.modifier |= DOSAGE_GLM;
	}
	if ((dosage_info.modifier & (DOSAGE_LIST | DOSAGE_SEPHEADER)) == DOSAGE_SEPHEADER) {
	  logerrprint("Error: --dosage 'sepheader' modifier must be used with 'list'.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if ((dosage_info.modifier & DOSAGE_DOSE1) && (dosage_info.format != 1)) {
	  logerrprint("Error: --dosage 'dose1' modifier must be used with 'format=1'.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	load_rare = LOAD_RARE_DOSAGE;
      } else if (!memcmp(argptr2, "fam", 4)) {
	UNSTABLE("dfam");
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "no-unrelateds")) {
	    family_info.dfam_modifier |= DFAM_NO_UNRELATEDS;
	  } else if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (family_info.dfam_modifier & DFAM_MPERM) {
	      logerrprint("Error: --dfam 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    family_info.dfam_modifier |= DFAM_PERM;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    family_info.dfam_modifier |= DFAM_PERM_COUNT;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (family_info.dfam_modifier & DFAM_PERM) {
	      logerrprint("Error: --dfam 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (family_info.dfam_modifier & DFAM_MPERM) {
	      logerrprint("Error: Duplicate --dfam 'mperm' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (scan_posint_defcapx(&(argv[cur_arg + uii][6]), &family_info.dfam_mperm_val)) {
	      sprintf(g_logbuf, "Error: Invalid --dfam mperm parameter '%s'.\n", &(argv[cur_arg + uii][6]));
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    family_info.dfam_modifier |= DFAM_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
	    family_info.dfam_modifier |= DFAM_SET_TEST;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    logerrprint("Error: Improper --dfam mperm syntax.  (Use \"--dfam mperm=<value>\".)\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --dfam parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_DFAM;
      } else if (!memcmp(argptr2, "fam-no-unrelateds", 18)) {
	// keep this undocumented flag since it makes DFAM correspond to the
	// original sib-TDT.
	if (!(calculation_type & CALC_DFAM)) {
	  logerrprint("Error: --dfam-no-unrelateds must be used with --dfam.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	family_info.dfam_modifier |= DFAM_NO_UNRELATEDS;
	logprint("Note: --dfam-no-unrelateds flag deprecated.  Use \"--dfam no-unrelateds\".\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "prime", 6)) {
	logprint("Note: --dprime flag deprecated.  Use e.g. \"--r2 dprime\".\n");
	ld_info.modifier |= LD_DPRIME;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'e':
      if (!memcmp(argptr2, "xtract", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp(argv[cur_arg + 1], "range")) {
            uii = 2;
	  } else if (strcmp(argv[cur_arg + 2], "range")) {
	    logerrprint("Error: Invalid --extract parameter sequence.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
          misc_flags |= MISC_EXTRACT_RANGE;
	}
	retval = alloc_fname(&extractname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "xclude", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp(argv[cur_arg + 1], "range")) {
            uii = 2;
	  } else if (strcmp(argv[cur_arg + 2], "range")) {
	    logerrprint("Error: Invalid --exclude parameter sequence.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
          misc_flags |= MISC_EXCLUDE_RANGE;
	}
	retval = alloc_fname(&excludename, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "xclude-snp", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&markername_snp, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
        filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV | FILTER_EXCLUDE_MARKERNAME_SNP;
      } else if (!memcmp(argptr2, "xclude-snps", 12)) {
	if (markername_snp) {
	  logerrprint("Error: --exclude-snps cannot be used with --exclude-snp.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &snps_range_list, 0);
	if (retval) {
	  goto main_ret_1;
	}
        filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV | FILTER_EXCLUDE_MARKERNAME_SNP;
      } else if (!memcmp(argptr2, "pistasis", 9)) {
	if (epi_info.modifier & EPI_FAST_CASE_ONLY) {
	  logerrprint("Error: --epistasis cannot be used with --case-only.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
        if (param_ct) {
          if (!strcmp(argv[cur_arg + 1], "set-by-set")) {
	    epi_info.modifier |= EPI_SET_BY_SET;
	  } else if (!strcmp(argv[cur_arg + 1], "set-by-all")) {
	    epi_info.modifier |= EPI_SET_BY_ALL;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --epistasis modifier '%s'.\n", argv[cur_arg + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        epi_info.modifier |= EPI_REG;
	calculation_type |= CALC_EPI;
      } else if (!memcmp(argptr2, "pistasis-summary-merge", 23)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	// .summary.32767 = 14 chars
	retval = alloc_fname(&epi_info.summary_merge_prefix, argv[cur_arg + 1], argptr, 14);
	if (retval) {
	  goto main_ret_1;
	}
	if (scan_posint_cappedx(argv[cur_arg + 2], PARALLEL_MAX, &epi_info.summary_merge_ct) || (epi_info.summary_merge_ct == 1)) {
	  sprintf(g_logbuf, "Error: Invalid --epistasis-summary-merge job count '%s'.\n", argv[cur_arg + 2]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "pi1", 4)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx <= 0)) {
	  sprintf(g_logbuf, "Error: Invalid --epi1 parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	epi_info.epi1 = dxx;
      } else if (!memcmp(argptr2, "pi2", 4)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx <= 0) || (dxx >= 1)) {
	  sprintf(g_logbuf, "Error: Invalid --epi2 parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	epi_info.epi2 = dxx;
      } else if (!memcmp(argptr2, "xclude-before-extract", 22)) {
        logprint("Note: --exclude-before-extract has no effect.\n");
	goto main_param_zero;
      } else if ((!memcmp(argptr2, "m-follow", 9)) ||
                 (!memcmp(argptr2, "m-meta-iter", 12)) ||
		 // may as well print an appropriate error message when an
		 // unsupported undocumented flag is used, but I'll draw a line
		 // at *misspelled* undocumented flags--if anyone was using
		 // --em-meta-likilood, the misspelling should have been
		 // caught...
                 (!memcmp(argptr2, "m-meta-prune-haplotype", 23)) ||
                 (!memcmp(argptr2, "m-meta-prune-phase", 19)) ||
                 (!memcmp(argptr2, "m-meta-tol", 11)) ||
                 (!memcmp(argptr2, "m-meta-window", 14)) ||
                 (!memcmp(argptr2, "m-overlap", 10)) ||
                 (!memcmp(argptr2, "m-verbose", 10)) ||
                 (!memcmp(argptr2, "m-window", 9)) ||
                 (!memcmp(argptr2, "m-window-iter", 14)) ||
                 (!memcmp(argptr2, "m-window-likelihood", 20)) ||
                 (!memcmp(argptr2, "m-window-prune-haplotype", 25)) ||
                 (!memcmp(argptr2, "m-window-prune-phase", 21)) ||
                 (!memcmp(argptr2, "m-window-tol", 13))) {
        goto main_hap_disabled_message;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'f':
      if (!memcmp(argptr2, "ile", 4)) {
	if (load_params & (LOAD_PARAMS_BFILE_ALL | LOAD_PARAMS_OX_ALL)) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	load_params |= LOAD_PARAMS_TEXT_ALL;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  csptr = argv[cur_arg + 1];
	  if (strlen(csptr) > (FNAMESIZE - 5)) {
	    logerrprint("Error: --file parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  csptr = PROG_NAME_STR;
	}
	memcpy(strcpya(pedname, csptr), ".ped", 5);
	memcpy(strcpya(mapname, csptr), ".map", 5);
      } else if (!memcmp(argptr2, "am", 3)) {
	if (load_params & (LOAD_PARAMS_TEXT_ALL | LOAD_PARAMS_OX_ALL)) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	load_params |= LOAD_PARAMS_FAM;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logerrprint("Error: --fam parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(famname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "ilter", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&filtername, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        retval = alloc_and_flatten(&filtervals_flattened, &(argv[cur_arg + 2]), param_ct - 1);
	if (retval) {
	  goto main_ret_NOMEM;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "ilter-cases", 12)) {
	filter_flags |= FILTER_FAM_REQ | FILTER_BINARY_CASES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-controls", 15)) {
	if (filter_flags & FILTER_BINARY_CASES) {
	  logerrprint("Error: --filter-cases and --filter-controls cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_flags |= FILTER_FAM_REQ | FILTER_BINARY_CONTROLS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-females", 14)) {
	filter_flags |= FILTER_FAM_REQ | FILTER_BINARY_FEMALES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-males", 12)) {
	if (filter_flags & FILTER_BINARY_FEMALES) {
	  logerrprint("Error: --filter-males and --filter-females cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_flags |= FILTER_FAM_REQ | FILTER_BINARY_MALES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-founders", 15)) {
	filter_flags |= FILTER_FAM_REQ | FILTER_BINARY_FOUNDERS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-nonfounders", 18)) {
	if (filter_flags & FILTER_BINARY_FOUNDERS) {
	  logerrprint("Error: --filter-founders and --filter-nonfounders cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_flags |= FILTER_FAM_REQ | FILTER_BINARY_NONFOUNDERS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "req", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "counts")) {
	    if (misc_flags & MISC_FREQ_CC) {
	      logerrprint("Error: --freq 'counts' and 'case-control' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    misc_flags |= MISC_FREQ_COUNTS;
	  } else if (!strcmp(argv[cur_arg + uii], "case-control")) {
	    if (misc_flags & MISC_FREQ_COUNTS) {
	      logerrprint("Error: --freq 'counts' and 'case-control' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    misc_flags |= MISC_FREQ_CC;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    misc_flags |= MISC_FREQ_GZ;
	  } else {
            sprintf(g_logbuf, "Error: Invalid --freq parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_FREQ;
	if (misc_flags & MISC_FREQ_COUNTS) {
	  // --keep-allele-order also set for backward compatibility
	  // placed here instead of a few lines up because "--freq --counts" is
	  // permitted
	  misc_flags |= MISC_KEEP_ALLELE_ORDER;
	}
      } else if (!memcmp(argptr2, "reqx", 5)) {
	if (calculation_type & CALC_FREQ) {
	  logerrprint("Error: --freqx cannot be used with --freq.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "gz")) {
	    sprintf(g_logbuf, "Error: Invalid --freqx parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  misc_flags |= MISC_FREQ_GZ;
	}
	calculation_type |= CALC_FREQ;
	misc_flags |= MISC_FREQX;
      } else if (!memcmp(argptr2, "rom", 4)) {
	if (chrom_flag_present) {
	  logerrprint("Error: --from cannot be used with --autosome[-xy] or --[not-]chr.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&markername_from, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if ((!memcmp(argptr2, "rom-bp", 7)) || (!memcmp(argptr2, "rom-kb", 7)) || (!memcmp(argptr2, "rom-mb", 7))) {
	if (markername_from) {
	  logerrprint("Error: --from-bp/-kb/-mb cannot be used with --from.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	cc = argptr2[4];
	if (cc == 'b') {
	  if (scan_uint_defcapx(argv[cur_arg + 1], (uint32_t*)&marker_pos_start)) {
	    sprintf(g_logbuf, "Error: Invalid --from-bp parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else {
	  if (marker_pos_start != -1) {
	    logerrprint("Error: Multiple --from-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (scan_doublex(argv[cur_arg + 1], &dxx)) {
	    sprintf(g_logbuf, "Error: Invalid --from-kb/-mb parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  dxx *= (cc == 'k')? 1000 : 1000000;
	  if (dxx < 0) {
	    marker_pos_start = 0;
	  } else if (dxx > 2147483646) {
	    marker_pos_start = 0x7ffffffe;
	  } else {
	    marker_pos_start = (int32_t)(dxx * (1 + SMALL_EPSILON));
	  }
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP;
      } else if (!memcmp(argptr2, "isher", 6)) {
	if (model_modifier & MODEL_ASSOC) {
	  logerrprint("Error: --fisher cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --fisher flag deprecated.  Use \"--assoc fisher\" or \"--model fisher\".\n");
	model_modifier |= MODEL_ASSOC | MODEL_FISHER | MODEL_ASSOC_FDEPR;
	calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "id", 3)) {
        logprint("Note: --fid flag deprecated.  Use \"--recode vcf-fid\".\n");
	recode_modifier |= RECODE_FID;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "lip", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&flip_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "lip-subset", 11)) {
	if (!flip_fname) {
          logerrprint("Error: --flip-subset must be used with --flip.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (allelexxxx) {
	  // fix for this is too messy to be worthwhile
	  logerrprint("Error: --flip-subset cannot be used with --allele1234/ACGT.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&flip_subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ast-epistasis", 14)) {
	if (epi_info.modifier & EPI_REG) {
	  logerrprint("Error: --fast-epistasis cannot be used with --epistasis.\n");
          goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "no-ueki")) {
	    if (epi_info.modifier & (EPI_FAST_BOOST | EPI_FAST_JOINT_EFFECTS)) {
	      sprintf(g_logbuf, "Error: --fast-epistasis 'no-ueki' modifier cannot be used with '%s'.\n", (epi_info.modifier & EPI_FAST_BOOST)? "boost" : "joint-effects");
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    epi_info.modifier |= EPI_FAST_NO_UEKI;
	  } else if (!strcmp(argv[cur_arg + uii], "boost")) {
	    if (epi_info.modifier & (EPI_FAST_NO_UEKI | EPI_FAST_JOINT_EFFECTS)) {
	      sprintf(g_logbuf, "Error: --fast-epistasis 'boost' modifier cannot be used with '%s'.\n", (epi_info.modifier & EPI_FAST_NO_UEKI)? "no-ueki" : "joint-effects");
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    if (epi_info.modifier & EPI_FAST_CASE_ONLY) {
	      logerrprint("Error: --fast-epistasis boost does not have a case-only mode.\n");
              goto main_ret_INVALID_CMDLINE;
	    }
	    epi_info.modifier |= EPI_FAST_BOOST;
	  } else if (!strcmp(argv[cur_arg + uii], "joint-effects")) {
	    if (epi_info.modifier & (EPI_FAST_NO_UEKI | EPI_FAST_BOOST)) {
	      sprintf(g_logbuf, "Error: --fast-epistasis 'joint-effects' modifier cannot be used with '%s'.\n", (epi_info.modifier & EPI_FAST_NO_UEKI)? "no-ueki" : "boost");
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    epi_info.modifier |= EPI_FAST_JOINT_EFFECTS;
	  } else if (!strcmp(argv[cur_arg + uii], "case-only")) {
	    if (epi_info.modifier & EPI_FAST_BOOST) {
              logerrprint("Error: --fast-epistasis boost does not have a case-only mode.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    epi_info.modifier |= EPI_FAST_CASE_ONLY;
	  } else if (!strcmp(argv[cur_arg + uii], "set-by-set")) {
            if (!(epi_info.modifier & EPI_SET_BY_ALL)) {
	      epi_info.modifier |= EPI_SET_BY_SET;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "set-by-all")) {
	    if (epi_info.modifier & EPI_SET_BY_SET) {
	      epi_info.modifier -= EPI_SET_BY_SET;
	    }
            epi_info.modifier |= EPI_SET_BY_ALL;
	  } else if (!strcmp(argv[cur_arg + uii], "nop")) {
	    epi_info.modifier |= EPI_FAST_NO_P_VALUE;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --fast-epistasis modifier '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        epi_info.modifier |= EPI_FAST;
	calculation_type |= CALC_EPI;
      } else if (!memcmp(argptr2, "lip-scan", 9)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (param_ct) {
          if (strcmp(argv[cur_arg + 1], "verbose")) {
	    sprintf(g_logbuf, "Error: Invalid --flip-scan parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
          ld_info.modifier |= LD_FLIPSCAN_VERBOSE;
	}
        calculation_type |= CALC_FLIPSCAN;
      } else if (!memcmp(argptr2, "lip-scan-window", 16)) {
        if (!(calculation_type & CALC_FLIPSCAN)) {
	  logerrprint("Error: --flip-scan-window must be used with --flip-scan.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &ld_info.flipscan_window_size) || (ld_info.flipscan_window_size == 1)) {
	  sprintf(g_logbuf, "Error: Invalid --flip-scan-window size '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "lip-scan-window-kb", 22)) {
        if (!(calculation_type & CALC_FLIPSCAN)) {
	  logerrprint("Error: --flip-scan-window-kb must be used with --flip-scan.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --flip-scan-window-kb parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (dxx > 2147483.646) {
	  ld_info.flipscan_window_bp = 2147483646;
	} else {
	  ld_info.flipscan_window_bp = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
	}
      } else if (!memcmp(argptr2, "lip-scan-threshold", 19)) {
        if (!(calculation_type & CALC_FLIPSCAN)) {
	  logerrprint("Error: --flip-scan-threshold must be used with --flip-scan.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx > 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --flip-scan-threshold parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        ld_info.flipscan_thresh = dxx;
      } else if (!memcmp(argptr2, "lip-scan-verbose", 17)) {
	if (!(calculation_type & CALC_FLIPSCAN)) {
	  logerrprint("Error: --flip-scan-verbose must be used with --flip-scan.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --flip-scan-verbose flag deprecated.  Use \"--flip-scan verbose\".\n");
        ld_info.modifier |= LD_FLIPSCAN_VERBOSE;
      } else if (!memcmp(argptr2, "amily", 6)) {
	if (calculation_type & CALC_DFAM) {
	  logerrprint("Error: --family cannot be used with --dfam.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        misc_flags |= MISC_FAMILY_CLUSTERS;
	filter_flags |= FILTER_FAM_REQ;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ill-missing-a2", 15)) {
	if (load_rare & LOAD_RARE_CNV) {
	  logerrprint("Error: --fill-missing-a2 cannot be used with a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (flip_subset_fname) {
	  logerrprint("Error: --fill-missing-a2 cannot be used with --flip-subset.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	misc_flags |= MISC_FILL_MISSING_A2;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "st", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  // allow case/control status to represent just two subpopulations,
	  // but force user to be explicit about this nonstandard usage
          if (strcmp(argv[cur_arg + 1], "case-control")) {
	    sprintf(g_logbuf, "Error: Invalid --fst parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  misc_flags |= MISC_FST_CC;
	}
        calculation_type |= CALC_FST;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'g':
      if (!memcmp(argptr2, "eno", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_doublex(argv[cur_arg + 1], &geno_thresh)) {
	    sprintf(g_logbuf, "Error: Invalid --geno parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if ((geno_thresh < 0.0) || (geno_thresh > 1.0)) {
	    sprintf(g_logbuf, "Error: Invalid --geno parameter '%s' (must be between 0 and 1).\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else {
	  geno_thresh = 0.1;
	}
	filter_flags |= FILTER_ALL_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "en", 3)) {
	if (load_rare || (load_params & (LOAD_PARAMS_TEXT_ALL | LOAD_PARAMS_BFILE_ALL | LOAD_PARAMS_OXBGEN))) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	load_params |= LOAD_PARAMS_OXGEN;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logerrprint("Error: --gen parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "enome", 6)) {
	if (genome_modifier & GENOME_OUTPUT_GZ) {
          logerrprint("Warning: Duplicate --genome flag.  (--Z-genome is treated as \"--genome gz\".)\n");
	}
	kk = 0;
      main_genome_flag:
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 5 - kk)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "gz")) {
            genome_modifier |= GENOME_OUTPUT_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "rel-check")) {
            genome_modifier |= GENOME_REL_CHECK;
	  } else if (!strcmp(argv[cur_arg + uii], "full")) {
	    genome_modifier |= GENOME_OUTPUT_FULL;
	  } else if (!strcmp(argv[cur_arg + uii], "unbounded")) {
	    genome_modifier |= GENOME_IBD_UNBOUNDED;
	  } else if (!strcmp(argv[cur_arg + uii], "nudge")) {
            genome_modifier |= GENOME_NUDGE;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --genome parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_GENOME;
      } else if (!memcmp(argptr2, "enome-full", 11)) {
	if (!(calculation_type & CALC_GENOME)) {
	  logerrprint("Error: --genome-full must be used with --genome.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --genome-full flag deprecated.  Use \"--genome full\".\n");
	genome_modifier |= GENOME_OUTPUT_FULL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "roupdist", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_posintptrx(argv[cur_arg + 1], &groupdist_iters) || (groupdist_iters < 2) || (groupdist_iters > ((~ZEROLU) - MAX_THREADS))) {
	    sprintf(g_logbuf, "Error: Invalid --groupdist jackknife iteration count '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (param_ct == 2) {
	    if (scan_posint_defcapx(argv[cur_arg + 2], &groupdist_d)) {
	      sprintf(g_logbuf, "Error: Invalid --groupdist jackknife delete parameter '%s'.\n", argv[cur_arg + 2]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	}
	calculation_type |= CALC_GROUPDIST;
      } else if (!memcmp(argptr2, "rm", 3)) {
	logerrprint("Error: --grm has been retired due to inconsistent meaning across GCTA versions.\nUse --grm-gz or --grm-bin.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "rm-gz", 6)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 8)) {
	    logerrprint("Error: --grm-gz parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, sptr);
	} else {
	  memcpy(pedname, PROG_NAME_STR, 6);
	}
        load_rare = LOAD_RARE_GRM;
      } else if (!memcmp(argptr2, "rm-bin", 7)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 11)) {
	    logerrprint("Error: --grm-bin parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, sptr);
	} else {
	  memcpy(pedname, PROG_NAME_STR, 6);
	}
        load_rare = LOAD_RARE_GRM_BIN;
      } else if (!memcmp(argptr2, "xe", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (!covar_fname) {
	  logerrprint("Error: --gxe must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (param_ct) {
	  if (scan_posint_defcapx(argv[cur_arg + 1], &gxe_mcovar)) {
	    sprintf(g_logbuf, "Error: Invalid --gxe parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else {
	  gxe_mcovar = 1;
	}
	calculation_type |= CALC_GXE;
      } else if (!memcmp(argptr2, "enedrop", 8)) {
	if (model_modifier & MODEL_QMASK) {
	  logerrprint("Error: --assoc 'qt-means'/'lin' does not make sense with --genedrop.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --genedrop flag deprecated.  Use e.g. \"--model genedrop\".\n");
	model_modifier |= MODEL_GENEDROP;
	glm_modifier |= GLM_GENEDROP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "c", 2)) {
        if (!mtest_adjust) {
	  logerrprint("Error: --gc must be used with --adjust.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --gc flag deprecated.  Use \"--adjust gc\".\n");
	mtest_adjust |= ADJUST_GC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "file", 5)) {
	UNSTABLE("gfile");
	if (load_rare || (load_params & (~LOAD_PARAMS_FAM))) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	sptr = argv[cur_arg + 1];
	uii = strlen(sptr);
	if (uii > (FNAMESIZE - 6)) {
	  logerrprint("Error: --gfile parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(memcpya(pedname, sptr, uii), ".gvar", 6);
	if (!(load_params & LOAD_PARAMS_FAM)) {
	  memcpy(memcpya(famname, sptr, uii), ".fam", 5);
	}
	memcpy(memcpya(mapname, sptr, uii), ".map", 5);
	load_rare = LOAD_RARE_GVAR;
      } else if (!memcmp(argptr2, "enome-lists", 12)) {
	logerrprint("Error: --genome-lists flag retired.  Use --parallel.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "enome-minimal", 14)) {
	logerrprint("Error: --genome-minimal flag retired.  Use \"--genome gz\".\n");
        goto main_ret_INVALID_CMDLINE;
      } else if ((!memcmp(argptr2, "roup-avg", 9)) || (!memcmp(argptr2, "roup-average", 13))) {
        if (!(calculation_type & CALC_CLUSTER)) {
	  logerrprint("Error: --group-avg must be used with --cluster.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (cluster.modifier & CLUSTER_OLD_TIEBREAKS) {
	  logerrprint("Error: --cluster 'group-avg' and 'old-tiebreaks' cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        LOGPRINTF("Note: --%s flag deprecated.  Use \"--cluster group-avg\".\n", argptr);
	cluster.modifier |= CLUSTER_GROUP_AVG;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "enotypic", 9)) {
	if (glm_modifier & GLM_DOMINANT) {
	  logerrprint("Error: --genotypic cannot be used with --dominant.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --genotypic flag deprecated.  Use e.g. \"--linear genotypic\".\n");
	glm_modifier |= GLM_GENOTYPIC;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ene", 4)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_and_flatten(&set_info.genekeep_flattened, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "ene-all", 8)) {
	if (set_info.genekeep_flattened) {
          logerrprint("Error: --gene-all cannot be used with --gene.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
        set_info.modifier |= SET_GENE_ALL;
	filter_flags |= FILTER_BIM_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ap", 3)) {
	if ((epi_info.modifier & (EPI_FAST | EPI_FAST_CASE_ONLY)) != (EPI_FAST | EPI_FAST_CASE_ONLY)) {
	  logerrprint("Error: --gap must be used with \"--fast-epistasis case-only\".\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --gap parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (dxx > 2147483.646) {
	  epi_info.case_only_gap = 2147483646;
	} else {
          epi_info.case_only_gap = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "ene-report", 11)) {
	uii = gene_report_glist? 1 : 2;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], uii, uii)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&gene_report_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (uii == 2) {
	  retval = alloc_fname(&gene_report_glist, argv[cur_arg + 2], argptr, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "ene-list", 9)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&gene_report_glist, argv[cur_arg + 1], argptr, 0);
        if (retval) {
	  goto main_ret_1;
	}
        logprint("Note: --gene-list flag deprecated.  Pass two parameters to --gene-report\ninstead.\n");
      } else if (!memcmp(argptr2, "ene-list-border", 16)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --gene-list-border parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WW;
	}
	if (dxx > 2147483.646) {
	  gene_report_border = 0x7ffffffe;
	} else {
	  gene_report_border = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "ene-subset", 11)) {
	if (!gene_report_fname) {
	  logerrprint("Error: --gene-subset must be used with --gene-report.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&gene_report_subset, argv[cur_arg + 1], argptr, 0);
        if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ene-report-snp-field", 21)) {
	if (!extractname) {
	  logerrprint("Error: --gene-report-snp-field must be used with --extract.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (misc_flags & MISC_EXTRACT_RANGE) {
	  logerrprint("Error: --gene-report-snp-field cannot be used with \"--extract range\".\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&gene_report_snp_field, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "plink", 6)) {
        misc_flags |= MISC_GPLINK;
      } else if (!memcmp(argptr2, "ates", 5)) {
        logerrprint("Error: --gates is not implemented yet.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
        goto main_ret_1;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'h':
      if ((!memcmp(argptr2, "we", 3)) || (!memcmp(argptr2, "we midp", 8))) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (argptr2[2]) {
	  hwe_modifier |= HWE_THRESH_MIDP;
	}
	ujj = 0;
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "midp")) {
	    hwe_modifier |= HWE_THRESH_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "include-nonctrl")) {
	    hwe_modifier |= HWE_THRESH_ALL;
	  } else {
	    if (scan_doublex(argv[cur_arg + uii], &hwe_thresh) || ujj) {
	      logerrprint("Error: Invalid --hwe parameter sequence.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
            ujj = 1;
            if ((hwe_thresh < 0.0) || (hwe_thresh >= 1.0)) {
	      sprintf(g_logbuf, "Error: Invalid --hwe threshold '%s' (must be between 0 and 1).\n", argv[cur_arg + uii]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	}
	if (!ujj) {
	  logerrprint("Error: --hwe now requires a p-value threshold.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if ((hwe_modifier & HWE_THRESH_MIDP) && (hwe_thresh >= 0.5)) {
	  logerrprint("Error: --hwe threshold must be smaller than 0.5 when using mid-p adjustment.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	filter_flags |= FILTER_ALL_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "we-all", 7)) {
	logprint("Note: --hwe-all flag deprecated.  Use \"--hwe include-nonctrl\".\n");
	hwe_modifier |= HWE_THRESH_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "small-sample")) {
	    misc_flags |= MISC_HET_SMALL_SAMPLE;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    misc_flags |= MISC_HET_GZ;
	  } else {
            sprintf(g_logbuf, "Error: Invalid --het parameter '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        calculation_type |= CALC_HET;
      } else if ((!memcmp(argptr2, "ardy", 5)) || (!memcmp(argptr2, "ardy midp", 10))) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (argptr2[4]) {
	  hwe_modifier |= HWE_MIDP;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "midp")) {
            hwe_modifier |= HWE_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    hwe_modifier |= HWE_GZ;
	  } else {
            sprintf(g_logbuf, "Error: Invalid --hardy parameter '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_HARDY;
      } else if (!memcmp(argptr2, "omozyg", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "group")) {
	    if (homozyg.modifier & HOMOZYG_GROUP_VERBOSE) {
	      logerrprint("Error: --homozyg 'group' and 'group-verbose' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    homozyg.modifier |= HOMOZYG_GROUP;
	  } else if (!strcmp(argv[cur_arg + uii], "group-verbose")) {
	    if (homozyg.modifier & HOMOZYG_GROUP) {
	      logerrprint("Error: --homozyg 'group' and 'group-verbose' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    homozyg.modifier |= HOMOZYG_GROUP_VERBOSE;
	  } else if (!strcmp(argv[cur_arg + uii], "consensus-match")) {
	    homozyg.modifier |= HOMOZYG_CONSENSUS_MATCH;
	  } else if (!strcmp(argv[cur_arg + uii], "extend")) {
	    homozyg.modifier |= HOMOZYG_EXTEND;
	  } else if (!strcmp(argv[cur_arg + uii], "subtract-1-from-lengths")) {
            homozyg.modifier |= HOMOZYG_OLD_LENGTHS;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --homozyg parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_HOMOZYG;
      } else if (!memcmp(argptr2, "omozyg-snp", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &homozyg.min_snp) || (homozyg.min_snp == 1)) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-snp parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	calculation_type |= CALC_HOMOZYG;
      } else if (!memcmp(argptr2, "omozyg-kb", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < SMALL_EPSILON) || (dxx >= (2147483.646 * (1 + SMALL_EPSILON)))) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-kb parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	calculation_type |= CALC_HOMOZYG;
	// round up
	homozyg.min_bases = 1 + (uint32_t)((int32_t)(dxx * 1000 * (1 - SMALL_EPSILON)));
      } else if (!memcmp(argptr2, "omozyg-density", 15)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx >= 2147483.646)) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-density parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.max_bases_per_snp = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
      } else if (!memcmp(argptr2, "omozyg-gap", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0.001) || (dxx >= 2147483.646)) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-gap parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.max_gap = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
      } else if (!memcmp(argptr2, "omozyg-het", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &homozyg.max_hets)) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-het parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (homozyg.max_hets && (homozyg.modifier & HOMOZYG_EXTEND)) {
	  logerrprint("Error: --homozyg-het with a nonzero parameter cannot be used with --homozyg\nextend.  For fine-grained control over heterozygote frequency, use\n--homozyg-window-snp and --homozyg-window-het instead.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	calculation_type |= CALC_HOMOZYG;
      } else if (!memcmp(argptr2, "omozyg-window-snp", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &homozyg.window_size) || (homozyg.window_size == 1)) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-window-snp parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        calculation_type |= CALC_HOMOZYG;
      } else if (!memcmp(argptr2, "omozyg-window-kb", 17)) {
        logerrprint("Error: --homozyg-window-kb flag provisionally retired, since it had no effect\nin PLINK 1.07.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "omozyg-window-het", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &homozyg.window_max_hets)) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-window-het parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        calculation_type |= CALC_HOMOZYG;
      } else if (!memcmp(argptr2, "omozyg-window-missing", 22)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &homozyg.window_max_missing)) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-window-missing parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        calculation_type |= CALC_HOMOZYG;
      } else if (!memcmp(argptr2, "omozyg-window-threshold", 24)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx > 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-window-threshold parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.hit_threshold = dxx;
      } else if (!memcmp(argptr2, "omozyg-match", 13)) {
	if (!(homozyg.modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE))) {
          homozyg.modifier |= HOMOZYG_GROUP;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx > 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --homozyg-match parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	homozyg.overlap_min = dxx;
      } else if (!memcmp(argptr2, "omozyg-group", 13)) {
	if (homozyg.modifier & HOMOZYG_GROUP_VERBOSE) {
	  logprint("Note: --homozyg-group deprecated, and superseded by --homozyg group-verbose.\n");
	} else {
	  logprint("Note: --homozyg-group flag deprecated.  Use \"--homozyg group\".\n");
	  homozyg.modifier |= HOMOZYG_GROUP;
	}
	goto main_param_zero;
      } else if (!memcmp(argptr2, "omozyg-verbose", 15)) {
	if (!(homozyg.modifier & HOMOZYG_GROUP)) {
	  logerrprint("Error: --homozyg-verbose must be used with --homozyg group.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --homozyg-verbose flag deprecated.  Use \"--homozyg group-verbose\".\n");
	homozyg.modifier = (homozyg.modifier & (~HOMOZYG_GROUP)) | HOMOZYG_GROUP_VERBOSE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "omozyg-include-missing", 23)) {
        logerrprint("Error: --homozyg-include-missing flag provisionally retired, since it had no\neffect in PLINK 1.07.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ethom", 6)) {
	if (!(glm_modifier & GLM_GENOTYPIC)) {
	  logerrprint("Error: --hethom must be used with --genotypic.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --hethom flag deprecated.  Use e.g. \"--linear hethom\" (and\n\"--condition-list <filename> recessive\" to change covariate coding).\n");
	glm_modifier |= GLM_HETHOM | GLM_CONDITION_RECESSIVE;
	glm_modifier -= GLM_GENOTYPIC;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ide-covar", 10)) {
	logprint("Note: --hide-covar flag deprecated.  Use e.g. \"--linear hide-covar\".\n");
	glm_modifier |= GLM_HIDE_COVAR;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ard-call-threshold", 19)) {
        if (!(load_params & (LOAD_PARAMS_OX_ALL))) {
	  logerrprint("Error: --hard-call-threshold must be used with --data, --gen, or --bgen.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (!strcmp(argv[cur_arg + 1], "random")) {
	  hard_call_threshold = -1;
	} else {
	  if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0.0) || (dxx > 1.0)) {
	    sprintf(g_logbuf, "Error: Invalid --hard-call-threshold parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  } else if (dxx > (0.5 + SMALLISH_EPSILON)) {
	    sprintf(g_logbuf, "Error: The --hard-call-threshold parameter must be smaller than 0.5.  (Did you\nmean \"--hard-call-threshold %g\"?)\n", 1.0 - dxx);
	    goto main_ret_INVALID_CMDLINE_2A;
	  } else if (dxx > (0.5 - SMALLISH_EPSILON)) {
	    logerrprint("Error: The --hard-call-threshold parameter must be smaller than 0.5, to prevent\nties.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  hard_call_threshold = dxx * (1 + SMALL_EPSILON);
	}
      } else if (!memcmp(argptr2, "omog", 5)) {
        calculation_type |= CALC_HOMOG;
	goto main_param_zero;
      } else if ((!memcmp(argptr2, "ap", 3)) ||
                 (!memcmp(argptr2, "ap-all", 7)) ||
                 (!memcmp(argptr2, "ap-assoc", 9)) ||
                 (!memcmp(argptr2, "ap-freq", 8)) ||
                 (!memcmp(argptr2, "ap-impute", 10)) ||
                 (!memcmp(argptr2, "ap-impute-verbose", 18)) ||
                 (!memcmp(argptr2, "ap-linear", 10)) ||
                 (!memcmp(argptr2, "ap-logistic", 12)) ||
                 (!memcmp(argptr2, "ap-max-phase", 13)) ||
		 (!memcmp(argptr2, "ap-min-phase-prob", 18)) ||
                 (!memcmp(argptr2, "ap-miss", 8)) ||
                 (!memcmp(argptr2, "ap-omnibus", 11)) ||
                 (!memcmp(argptr2, "ap-only", 8)) ||
                 (!memcmp(argptr2, "ap-phase", 9)) ||
                 (!memcmp(argptr2, "ap-phase-wide", 14)) ||
                 (!memcmp(argptr2, "ap-pp", 6)) ||
                 (!memcmp(argptr2, "ap-snps", 8)) ||
                 (!memcmp(argptr2, "ap-tdt", 7)) ||
		 (!memcmp(argptr2, "ap-window", 10)) ||
                 (!memcmp(argptr2, "omozyg-haplo-track", 19))) {
      main_hap_disabled_message:
        logerrprint("Error: The --hap... family of flags has not been reimplemented in PLINK 1.9 due\nto poor phasing accuracy (and, consequently, inferior haplotype\nlikelihood/frequency estimates) relative to other software; for now, we\nrecommend using BEAGLE instead of PLINK for case/control haplotype association\nanalysis.  (You can use \"--recode beagle\" to export data.)  We apologize for\nthe inconvenience, and plan to develop variants of the --hap... flags which\nhandle pre-phased data effectively.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ard-call", 9)) {
	logerrprint("Error: The undocumented --hard-call flag has been retired.  (The\n--hard-call-threshold flag, supported by both PLINK and PLINK/SEQ, has similar\nfunctionality.)\n");
	goto main_ret_INVALID_CMDLINE;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'i':
      if (!memcmp(argptr2, "bc", 3)) {
	calculation_type |= CALC_IBC;
	goto main_param_zero;
      } else if ((!memcmp(argptr2, "ndep-pairwise", 14)) || (!memcmp(argptr2, "ndep-pairphase", 15))) {
	if (calculation_type & CALC_LD_PRUNE) {
	  logerrprint("Error: Conflicting --indep... commands.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 3, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	// may want to permit decimal here
	if (scan_posint_defcap(argv[cur_arg + 1], &ld_info.prune_window_size) || ((ld_info.prune_window_size == 1) && (param_ct == 3))) {
	  sprintf(g_logbuf, "Error: Invalid --%s window size '%s'.\n", argptr, argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (param_ct == 4) {
	  if (!match_upper(argv[cur_arg + 2], "KB")) {
	    sprintf(g_logbuf, "Error: Invalid --%s parameter sequence.\n", argptr);
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  ld_info.modifier |= LD_PRUNE_KB_WINDOW;
	} else {
	  uii = strlen(argv[cur_arg + 1]);
	  if ((uii > 2) && match_upper(&(argv[cur_arg + 1][uii - 2]), "KB")) {
	    ld_info.modifier |= LD_PRUNE_KB_WINDOW;
	  }
	}
	if (scan_posint_defcapx(argv[cur_arg + param_ct - 1], &ld_info.prune_window_incr)) {
	  sprintf(g_logbuf, "Error: Invalid increment '%s' for --%s.\n", argv[cur_arg + param_ct - 1], argptr);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((ld_info.modifier & LD_PRUNE_KB_WINDOW) && (ld_info.prune_window_incr > 1)) {
	  LOGERRPRINTFWW("Warning: --%s step size should be 1 when window size is in kb units.\n", argptr);
	}
	if (scan_doublex(argv[cur_arg + param_ct], &ld_info.prune_last_param) || (ld_info.prune_last_param < 0.0) || (ld_info.prune_last_param >= 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --%s r^2 threshold '%s'.\n", argptr, argv[cur_arg + param_ct]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	calculation_type |= CALC_LD_PRUNE;
	if (argptr2[9] == 'w') {
          ld_info.modifier |= LD_PRUNE_PAIRWISE;
	} else {
          ld_info.modifier |= LD_PRUNE_PAIRPHASE;
	}
      } else if (!memcmp(argptr2, "ndep", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 3, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &ld_info.prune_window_size) || ((ld_info.prune_window_size == 1) && (param_ct == 3))) {
	  sprintf(g_logbuf, "Error: Invalid --indep window size '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (param_ct == 4) {
	  if (!match_upper(argv[cur_arg + 2], "KB")) {
	    logerrprint("Error: Invalid --indep parameter sequence.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  ld_info.modifier |= LD_PRUNE_KB_WINDOW;
	} else {
	  uii = strlen(argv[cur_arg + 1]);
	  if ((uii > 2) && match_upper(&(argv[cur_arg + 1][uii - 2]), "KB")) {
	    ld_info.modifier |= LD_PRUNE_KB_WINDOW;
	  }
	}
	if (scan_posint_defcapx(argv[cur_arg + param_ct - 1], &ld_info.prune_window_incr)) {
	  sprintf(g_logbuf, "Error: Invalid increment '%s' for --indep.\n", argv[cur_arg + param_ct - 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (scan_doublex(argv[cur_arg + param_ct], &ld_info.prune_last_param)) {
	  sprintf(g_logbuf, "Error: Invalid --indep VIF threshold '%s'.\n", argv[cur_arg + param_ct]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (ld_info.prune_last_param < 1.0) {
	  sprintf(g_logbuf, "Error: --indep VIF threshold '%s' too small (must be >= 1).\n", argv[cur_arg + param_ct]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	calculation_type |= CALC_LD_PRUNE;
      } else if (!memcmp(argptr2, "ndiv-sort", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	jj = (argv[cur_arg + 1][1] == '\0');
	if ((!strcmp(argv[cur_arg + 1], "none")) || ((argv[cur_arg + 1][0] == '0') && jj)) {
	  sample_sort = SAMPLE_SORT_NONE;
	} else if ((!strcmp(argv[cur_arg + 1], "natural")) || ((tolower(argv[cur_arg + 1][0]) == 'n') && jj)) {
	  sample_sort = SAMPLE_SORT_NATURAL;
	} else if ((!strcmp(argv[cur_arg + 1], "ascii")) || ((tolower(argv[cur_arg + 1][0]) == 'a') && jj)) {
	  sample_sort = SAMPLE_SORT_ASCII;
	} else if ((!strcmp(argv[cur_arg + 1], "file")) || ((tolower(argv[cur_arg + 1][0]) == 'f') && jj)) {
	  if (param_ct == 1) {
	    sprintf(g_logbuf, "Error: Missing \"--indiv-sort %s\" filename.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  sample_sort = SAMPLE_SORT_FILE;
	  retval = alloc_fname(&sample_sort_fname, argv[cur_arg + 2], argptr, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	} else {
	  sprintf(g_logbuf, "Error: '%s' is not a valid mode for --indiv-sort.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((param_ct == 2) && (sample_sort != SAMPLE_SORT_FILE)) {
          sprintf(g_logbuf, "Error: \"--indiv-sort %s\" does not accept a second parameter.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_2A;
	}
      } else if (!memcmp(argptr2, "bs-test", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_posintptrx(argv[cur_arg + 1], &ibs_test_perms)) {
	    sprintf(g_logbuf, "Error: Invalid --ibs-test permutation count '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
          if (ibs_test_perms < MAX_THREADS * 2) {
	    sprintf(g_logbuf, "Error: --ibs-test permutation count '%s' too small (min %u).\n", argv[cur_arg + 1], MAX_THREADS * 2);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_IBS_TEST;
      } else if (!memcmp(argptr2, "id", 3)) {
        logprint("Note: --iid flag deprecated.  Use \"--recode vcf-iid\".\n");
	recode_modifier |= RECODE_IID;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "bm", 3)) {
        if (!(calculation_type & CALC_CLUSTER)) {
	  logerrprint("Error: --ibm must be used with --cluster.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx)) {
	  sprintf(g_logbuf, "Error: Invalid --ibm parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((dxx <= 0.0) || (dxx > 1.0)) {
	  logerrprint("Error: --ibm threshold must be in (0, 1].\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        cluster.min_ibm = dxx;
      } else if (!memcmp(argptr2, "mpossible", 10)) {
	logerrprint("Error: --impossible flag retired.  Use \"--genome nudge\", or explicitly validate\nZ0/Z1/Z2/PI_HAT in your script.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "nteraction", 11)) {
	logprint("Note: --interaction flag deprecated.  Use e.g. \"--linear interaction\".\n");
	glm_modifier |= GLM_INTERACTION;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "bs-matrix", 10)) {
        calculation_type |= CALC_PLINK1_IBS_MATRIX;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "d-delim", 8)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (param_ct) {
          id_delim = extract_char_param(argv[cur_arg + 1]);
	  if (!id_delim) {
	    logerrprint("Error: --id-delim parameter must be a single character.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  } else if (((unsigned char)id_delim) < ' ') {
	    logerrprint("Error: --id-delim parameter cannot be tab, newline, or a nonprinting character.\n");
            goto main_ret_INVALID_CMDLINE;
	  }
	} else {
          id_delim = '_';
	}
      } else if (!memcmp(argptr2, "nter-chr", 9)) {
        logprint("Note: --inter-chr flag deprecated.  Use e.g. \"--r2 inter-chr\".\n");
	ld_info.modifier |= LD_INTER_CHR;
      } else if (!memcmp(argptr2, "nd-major", 9)) {
	logerrprint("Error: --ind-major is retired, to discourage creation of .bed files that\nconstantly have to be transposed back.  \"--recode ped\" exports sample-major\nfiles which are good enough for smaller jobs; we suggest transposing small data\nwindows on the fly when tackling large jobs.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "mpute-sex", 10)) {
	if (calculation_type & CALC_SEXCHECK) {
	  logerrprint("Error: --check-sex is redundant with --impute-sex.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 5)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	ujj = 0;
	if (param_ct) {
	  for (uii = 1; uii <= param_ct; uii++) {
            if (!strcmp(argv[cur_arg + uii], "y-only")) {
	      ujj = uii;
	      break;
	    }
	  }
	}
	if (!ujj) {
	  for (uii = 1; uii <= param_ct; uii++) {
	    if (!strcmp(argv[cur_arg + uii], "ycount")) {
	      misc_flags |= MISC_SEXCHECK_YCOUNT;
	    } else {
	      if (!ujj) {
		if (scan_doublex(argv[cur_arg + uii], &check_sex_fthresh) || (check_sex_fthresh <= 0.0)) {
		  sprintf(g_logbuf, "Error: Invalid --impute-sex female F-statistic estimate ceiling '%s'.\n", argv[cur_arg + uii]);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
	      } else if (ujj == 1) {
		if (scan_doublex(argv[cur_arg + uii], &check_sex_mthresh) || (check_sex_mthresh >= 1.0)) {
		  sprintf(g_logbuf, "Error: Invalid --impute-sex male F-statistic estimate floor '%s'.\n", argv[cur_arg + uii]);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
	      } else if (ujj == 2) {
		if (scan_uint_defcapx(argv[cur_arg + uii], &check_sex_f_yobs)) {
		  logerrprint("Error: Invalid --impute-sex female Ychr maximum nonmissing genotype count.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
	      } else if (ujj == 3) {
		if (scan_uint_defcapx(argv[cur_arg + uii], &check_sex_m_yobs)) {
		  logerrprint("Error: Invalid --impute-sex male Ychr minimum nonmissing genotype count.\n");
		  goto main_ret_INVALID_CMDLINE_A;
		}
	      } else {
		logerrprint("Error: Invalid --impute-sex parameter sequence.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      ujj++;
	    }
	  }
	  if ((ujj > 2) && (!(misc_flags & MISC_SEXCHECK_YCOUNT))) {
	    logerrprint("Error: --impute-sex only accepts >2 numeric parameters in 'ycount' mode.\n");
            goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (check_sex_fthresh > check_sex_mthresh) {
	    logerrprint("Error: --impute-sex female F estimate ceiling cannot be larger than male floor.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	} else {
	  check_sex_m_yobs = 1;
	  ukk = 0;
	  for (uii = 1; uii <= param_ct; uii++) {
	    if (uii == ujj) {
	      continue;
	    }
	    if (!strcmp(argv[cur_arg + uii], "ycount")) {
	      logerrprint("Error: Conflicting --impute-sex modes.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (!ukk) {
	      if (scan_uint_defcapx(argv[cur_arg + uii], &check_sex_f_yobs)) {
		logerrprint("Error: Invalid --impute-sex y-only female Ychr maximum nonmissing genotype\ncount.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else if (ukk == 1) {
	      if (scan_uint_defcapx(argv[cur_arg + uii], &check_sex_m_yobs)) {
		logerrprint("Error: Invalid --impute-sex y-only male Ychr minimum nonmissing genotype count.\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else {
	      logerrprint("Error: Invalid --impute-sex y-only parameter sequence.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    ukk++;
	  }
	  if (check_sex_f_yobs >= check_sex_m_yobs) {
	    logerrprint("Error: In y-only mode, --impute-sex female Y observation threshold must be\nsmaller than male threshold.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  misc_flags |= MISC_SEXCHECK_YCOUNT | MISC_SEXCHECK_YONLY;
	}
        calculation_type |= CALC_SEXCHECK;
        misc_flags |= MISC_IMPUTE_SEX;
	sex_missing_pheno |= ALLOW_NO_SEX;
      } else if ((!memcmp(argptr2, "d-dict", 7)) || (!memcmp(argptr2, "d-dump", 7)) || (!memcmp(argptr2, "d-lookup", 9)) || (!memcmp(argptr2, "d-match", 8)) || (!memcmp(argptr2, "d-replace", 10)) || (!memcmp(argptr2, "d-table", 8))) {
	logerrprint("Error: --id-dict and --id-match are provisionally retired, since free database\nsoftware handles these operations in a more flexible and powerful manner.\nContact the developers if you still need them.\n");
	goto main_ret_INVALID_CMDLINE;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'j':
      if (!memcmp(argptr2, "e-cellmin", 10)) {
        if (!(epi_info.modifier & EPI_FAST_JOINT_EFFECTS)) {
	  logerrprint("Error: --je-cellmin must be used with \"--fast-epistasis joint-effects\".\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	// may as well enforce 2^29 / 18 limit...
	if (scan_uint_cappedx(argv[cur_arg + 1], 29826161, &epi_info.je_cellmin)) {
	  sprintf(g_logbuf, "Error: Invalid --je-cellmin parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'k':
      if (!memcmp(argptr2, "eep", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&keepname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "eep-fam", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&keepfamname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "eep-allele-order", 17)) {
	if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --keep-allele-order has no effect with %s.\n", (load_rare == LOAD_RARE_CNV)? "a .cnv fileset" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	misc_flags |= MISC_KEEP_ALLELE_ORDER;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "eep-before-remove", 18)) {
        logprint("Note: --keep-before-remove has no effect.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "eep-autoconv", 13)) {
	if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --keep-autoconv has no effect with %s.\n", (load_rare == LOAD_RARE_CNV)? "--cfile/--cnv-list" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        misc_flags |= MISC_KEEP_AUTOCONV;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "eep-cluster-names", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_and_flatten(&(cluster.keep_flattened), &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "eep-clusters", 13)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&(cluster.keep_fname), argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'l':
      if (!memcmp(argptr2, "file", 5)) {
	if (load_rare || (load_params & (~LOAD_PARAMS_FAM))) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  csptr = argv[cur_arg + 1];
	  if (strlen(csptr) > FNAMESIZE - 6) {
	    logerrprint("Error: --lfile filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  csptr = PROG_NAME_STR;
	}
	memcpy(strcpya(pedname, csptr), ".lgen", 6);
	memcpy(strcpya(mapname, csptr), ".map", 5);
	if (!famname[0]) {
	  memcpy(strcpya(famname, csptr), ".fam", 5);
	}
	load_rare = LOAD_RARE_LGEN;
      } else if (!memcmp(argptr2, "gen", 4)) {
	if ((load_rare & (~LOAD_RARE_LGEN)) || (load_params & (~LOAD_PARAMS_FAM))) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if ((load_params != LOAD_PARAMS_FAM) && (!load_rare)) {
	  logerrprint("Error: --lgen must be used with --fam or --lfile.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logerrprint("Error: --lgen parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
	load_rare = LOAD_RARE_LGEN;
      } else if (!memcmp(argptr2, "oop-assoc", 10)) {
	if (pheno_modifier & PHENO_ALL) {
	  logerrprint("Error: --loop-assoc cannot be used with --all-pheno.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (misc_flags & MISC_FAMILY_CLUSTERS) {
	  logerrprint("Error: --loop-assoc cannot be used with --family.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (load_rare == LOAD_RARE_DOSAGE) {
	  logerrprint("Error: --loop-assoc cannot be used with --dosage.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	if (param_ct == 2) {
	  if ((strlen(argv[cur_arg + 1]) == 7) && (!memcmp(argv[cur_arg + 1], "keep-", 5)) && match_upper(&(argv[cur_arg + 1][5]), "NA")) {
	    uii = 2;
	  } else if ((strlen(argv[cur_arg + 2]) != 7) || memcmp(argv[cur_arg + 2], "keep-", 5) || (!match_upper(&(argv[cur_arg + 2][5]), "NA"))) {
            logerrprint("Error: Invalid --loop-assoc parameter sequence.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
          misc_flags |= MISC_LOAD_CLUSTER_KEEP_NA;
	}
	retval = alloc_fname(&loop_assoc_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "og10", 5)) {
        if (!mtest_adjust) {
	  logerrprint("Error: --log10 must be used with --adjust.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --log10 flag deprecated.  Use \"--adjust log10\".\n");
	mtest_adjust |= ADJUST_LOG10;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ambda", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &adjust_lambda)) {
	  sprintf(g_logbuf, "Error: Invalid --lambda parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (adjust_lambda < 1) {
	  logprint("Note: --lambda parameter set to 1.\n");
	  adjust_lambda = 1;
	}
	mtest_adjust |= ADJUST_LAMBDA;
      } else if (!memcmp(argptr2, "ist-23-indels", 14)) {
        calculation_type |= CALC_LIST_23_INDELS;
	goto main_param_zero;
      } else if ((!memcmp(argptr2, "inear", 6)) || (!memcmp(argptr2, "ogistic", 8))) {
#ifndef NOLAPACK
        if (calculation_type & CALC_GLM) {
	  logerrprint("Error: --logistic cannot be used with --linear.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
#endif
	if (*argptr2 == 'o') {
	  glm_modifier |= GLM_LOGISTIC;
#ifdef NOLAPACK
	} else {
	  logerrprint("Error: --linear requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	  goto main_ret_INVALID_CMDLINE;
#endif
	}
	if (load_rare & LOAD_RARE_DOSAGE) {
	  // make --dosage + modifier-free --linear/--logistic only issue a
	  // warning
	  if (param_ct) {
	    logerrprint("Error: --dosage cannot be used with --linear/--logistic modifiers.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  logprint("Note: --dosage automatically performs a regression; --linear/--logistic has no\nadditional effect.\n");
	} else {
	  if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 11)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  for (uii = 1; uii <= param_ct; uii++) {
	    if (!strcmp(argv[cur_arg + uii], "perm")) {
	      if (glm_modifier & GLM_MPERM) {
		sprintf(g_logbuf, "Error: --%s 'mperm' and 'perm' cannot be used together.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2A;
	      }
	      glm_modifier |= GLM_PERM;
	    } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	      if (glm_modifier & GLM_PERM) {
		sprintf(g_logbuf, "Error: --%s 'mperm' and 'perm' cannot be used together.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2A;
	      } else if (glm_modifier & GLM_MPERM) {
		sprintf(g_logbuf, "Error: Duplicate --%s 'mperm' modifier.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2;
	      }
	      if (scan_posint_defcapx(&(argv[cur_arg + uii][6]), &glm_mperm_val)) {
		sprintf(g_logbuf, "Error: Invalid --%s mperm parameter '%s'.\n", argptr, &(argv[cur_arg + uii][6]));
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      glm_modifier |= GLM_MPERM;
	    } else if (!strcmp(argv[cur_arg + uii], "genedrop")) {
	      glm_modifier |= GLM_GENEDROP;
	    } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	      glm_modifier |= GLM_PERM_COUNT;
	    } else if (!strcmp(argv[cur_arg + uii], "genotypic")) {
	      if (glm_modifier & (GLM_HETHOM | GLM_DOMINANT | GLM_RECESSIVE)) {
		sprintf(g_logbuf, "Error: Conflicting --%s parameters.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2;
	      }
	      glm_modifier |= GLM_GENOTYPIC;
	      glm_xchr_model = 0;
	    } else if (!strcmp(argv[cur_arg + uii], "hethom")) {
	      if (glm_modifier & (GLM_GENOTYPIC | GLM_DOMINANT | GLM_RECESSIVE)) {
		sprintf(g_logbuf, "Error: Conflicting --%s parameters.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2;
	      }
	      glm_modifier |= GLM_HETHOM;
	      glm_xchr_model = 0;
	    } else if (!strcmp(argv[cur_arg + uii], "dominant")) {
	      if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_RECESSIVE)) {
		sprintf(g_logbuf, "Error: Conflicting --%s parameters.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2;
	      }
	      glm_modifier |= GLM_DOMINANT;
	      glm_xchr_model = 0;
	    } else if (!strcmp(argv[cur_arg + uii], "recessive")) {
	      if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_DOMINANT)) {
		sprintf(g_logbuf, "Error: Conflicting --%s parameters.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2;
	      }
	      glm_modifier |= GLM_RECESSIVE;
	      glm_xchr_model = 0;
	    } else if (!strcmp(argv[cur_arg + uii], "no-snp")) {
	      if (mtest_adjust) {
		sprintf(g_logbuf, "Error: --%s no-snp cannot be used with --adjust.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2A;
	      }
	      // defer the rest of the check
	      glm_modifier |= GLM_NO_SNP;
	    } else if (!strcmp(argv[cur_arg + uii], "hide-covar")) {
	      glm_modifier |= GLM_HIDE_COVAR;
	    } else if (!strcmp(argv[cur_arg + uii], "sex")) {
	      if (glm_modifier & GLM_NO_X_SEX) {
		sprintf(g_logbuf, "Error: --%s 'sex' and 'no-x-sex' cannot be used together.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2A;
	      }
	      glm_modifier |= GLM_SEX;
	    } else if (!strcmp(argv[cur_arg + uii], "no-x-sex")) {
	      if (glm_modifier & GLM_SEX) {
		sprintf(g_logbuf, "Error: --%s 'sex' and 'no-x-sex' cannot be used together.\n", argptr);
		goto main_ret_INVALID_CMDLINE_2A;
	      }
	      glm_modifier |= GLM_NO_X_SEX;
	    } else if (!strcmp(argv[cur_arg + uii], "interaction")) {
	      glm_modifier |= GLM_INTERACTION;
	    } else if (!strcmp(argv[cur_arg + uii], "standard-beta")) {
	      if (glm_modifier & GLM_LOGISTIC) {
		logerrprint("Error: --logistic does not have a 'standard-beta' modifier.  (Did you mean\n--linear or 'beta'?)\n");
		goto main_ret_INVALID_CMDLINE_A;
	      }
	      glm_modifier |= GLM_STANDARD_BETA;
	    } else if (!strcmp(argv[cur_arg + uii], "intercept")) {
	      glm_modifier |= GLM_INTERCEPT;
	    } else if (!strcmp(argv[cur_arg + uii], "beta")) {
	      glm_modifier |= GLM_BETA;
	    } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
	      glm_modifier |= GLM_SET_TEST;
	    } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	      sprintf(g_logbuf, "Error: Improper --%s mperm syntax.  (Use \"--%s mperm=<value>\".)\n", argptr, argptr);
	      goto main_ret_INVALID_CMDLINE_2;
	    } else {
	      sprintf(g_logbuf, "Error: Invalid --%s parameter '%s'.\n", argptr, argv[cur_arg + uii]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	  if ((glm_modifier & GLM_NO_SNP) && (glm_modifier & GLM_NO_SNP_EXCL)) {
	    sprintf(g_logbuf, "Error: --%s 'no-snp' modifier conflicts with another modifier.\n", argptr);
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  calculation_type |= CALC_GLM;
	}
      } else if (!memcmp(argptr2, "d-xchr", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	cc = argv[cur_arg + 1][0];
	if ((cc < '1') || (cc > '3') || (argv[cur_arg + 1][1] != '\0')) {
	  sprintf(g_logbuf, "Error: Invalid --ld-xchr parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        if (cc == '2') {
          ld_info.modifier |= LD_IGNORE_X;
	} else if (cc == '3') {
	  ld_info.modifier |= LD_WEIGHTED_X;
	}
      } else if (!memcmp(argptr2, "asso", 5)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &lasso_h2) || (lasso_h2 > 1) || (lasso_h2 <= 0)) {
	  sprintf(g_logbuf, "Error: Invalid --lasso heritability estimate '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "report-zeroes")) {
	    misc_flags |= MISC_LASSO_REPORT_ZEROES;
	  } else if (lasso_minlambda > 0) {
            logerrprint("Error: Invalid --lasso parameter sequence.\n");
            goto main_ret_INVALID_CMDLINE_A;
	  } else if (scan_doublex(argv[cur_arg + uii], &lasso_minlambda) || (lasso_minlambda <= 0)) {
	    sprintf(g_logbuf, "Error: Invalid --lasso minimum lambda '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        calculation_type |= CALC_LASSO;
      } else if (!memcmp(argptr2, "asso-select-covars", 19)) {
	if (!(calculation_type & CALC_LASSO)) {
	  logerrprint("Error: --lasso-select-covars must be used with --lasso.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!covar_fname) {
	  logerrprint("Error: --lasso-select-covars must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        misc_flags |= MISC_LASSO_SELECT_COVARS;
	if (param_ct) {
	  retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &lasso_select_covars_range_list, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "asso-lambda", 12)) {
	UNSTABLE("lasso-lambda");
	if (lasso_minlambda != -1) {
	  logerrprint("Error: --lasso-lambda cannot be used with a --lasso minimum lambda setting.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &lasso_lambda_iters)) {
	  sprintf(g_logbuf, "Error: Invalid --lasso-lambda parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (param_ct == 2) {
	  if (scan_doublex(argv[cur_arg + 2], &dxx) || (dxx > 1) || (dxx <= 0)) {
	    sprintf(g_logbuf, "Error: Invalid --lasso-lambda heritability estimate '%s'.\n", argv[cur_arg + 2]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if ((calculation_type & CALC_LASSO) && (dxx != lasso_h2)) {
	    logerrprint("Error: --lasso and --lasso-lambda heritability estimates are different.\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else {
	    lasso_h2 = dxx;
	  }
	} else if (!(calculation_type & CALC_LASSO)) {
	  logerrprint("Error: --lasso-lambda requires a heritability estimate.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	calculation_type |= CALC_LASSO_LAMBDA;
      } else if (!memcmp(argptr2, "d-window", 9)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &ld_info.window_size) || (ld_info.window_size == 1)) {
	  sprintf(g_logbuf, "Error: Invalid --ld-window window size '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "d-window-kb", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --ld-window-kb parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (dxx > 2147483.646) {
	  ld_info.window_bp = 2147483646;
	} else {
	  ld_info.window_bp = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
	}
      } else if (!memcmp(argptr2, "d-window-cm", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --ld-window-cm parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	ld_info.window_cm = dxx * (1 + SMALL_EPSILON);
      } else if (!memcmp(argptr2, "d-window-r2", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0) || (dxx > 1)) {
	  sprintf(g_logbuf, "Error: Invalid --ld-window-r2 parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        ld_info.window_r2 = dxx;
      } else if (!memcmp(argptr2, "d-snp", 6)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&(ld_info.snpstr), argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "d-snps", 7)) {
        if (ld_info.snpstr) {
          logerrprint("Error: --ld-snps cannot be used with --ld-snp or --ld-snp-list.\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &(ld_info.snps_rl), 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "d-snp-list", 11)) {
        if (ld_info.snpstr) {
          logerrprint("Error: --ld-snp cannot be used with --ld-snp-list.\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&(ld_info.snpstr), argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	ld_info.modifier |= LD_SNP_LIST_FILE;
      } else if (!memcmp(argptr2, "d", 2)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (alloc_string(&(epi_info.ld_mkr1), argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
	if (alloc_string(&(epi_info.ld_mkr2), argv[cur_arg + 2])) {
	  goto main_ret_NOMEM;
	}
	if (param_ct == 3) {
	  if (strcmp(argv[cur_arg + 3], "hwe-midp")) {
	    sprintf(g_logbuf, "Error: Invalid --ld parameter '%s'.\n", argv[cur_arg + 3]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
          epi_info.modifier |= EPI_HWE_MIDP;
	}
        calculation_type |= CALC_EPI;
      } else if (!memcmp(argptr2, "ist-all", 8)) {
	ld_info.modifier |= LD_SHOW_TAGS_LIST_ALL;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "ist-duplicate-vars", 19)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "require-same-ref")) {
	    dupvar_modifier |= DUPVAR_REF;
	  } else if (!strcmp(argv[cur_arg + uii], "ids-only")) {
	    dupvar_modifier |= DUPVAR_IDS_ONLY;
	  } else if (!strcmp(argv[cur_arg + uii], "suppress-first")) {
	    dupvar_modifier |= DUPVAR_SUPPRESS_FIRST;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --list-duplicate-vars parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_DUPVAR;
      } else if (!memcmp(argptr2, "d-pred", 7)) {
	logerrprint("Error: --ld-pred is currently under development.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto main_ret_1;
      } else if ((!memcmp(argptr2, "ookup", 6)) ||
                 (!memcmp(argptr2, "ookup-list", 11)) ||
                 (!memcmp(argptr2, "ookup-gene", 11)) ||
                 (!memcmp(argptr2, "ookup-gene-kb", 14)) ||
                 (!memcmp(argptr2, "ookup-gene-list", 16))) {
        logerrprint("Error: --lookup commands have been retired since the Sullivan Lab web database\nis no longer operational.  Use e.g. PLINK/SEQ's lookup command instead.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "iability", 9)) {
	logerrprint("Error: --liability is provisionally retired.  Contact the developers if you\nneed this option.\n");
	goto main_ret_INVALID_CMDLINE;
      } else {
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'm':
      if (!memcmp(argptr2, "ap", 3)) {
	if (((load_params & (LOAD_PARAMS_BFILE_ALL | LOAD_PARAMS_OX_ALL)) || (load_rare & (~(LOAD_RARE_CNV | LOAD_RARE_GVAR)))) && ((load_rare != LOAD_RARE_DOSAGE) || (load_params != LOAD_PARAMS_FAM)) && (load_rare != LOAD_RARE_LGEN)) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	load_params |= LOAD_PARAMS_MAP;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logerrprint("Error: --map parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "issing-genotype", 16)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        cc = extract_char_param(argv[cur_arg + 1]);
	if (((unsigned char)cc <= ' ') || ((cc > '0') && (cc <= '4')) || (cc == 'A') || (cc == 'C') || (cc == 'G') || (cc == 'T')) {
	  sprintf(g_logbuf, "Error: Invalid --missing-genotype parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	g_missing_geno_ptr = &(g_one_char_strs[((unsigned char)cc) * 2]);
	g_output_missing_geno_ptr = g_missing_geno_ptr;
      } else if (!memcmp(argptr2, "issing-phenotype", 17)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	jj = strlen(argv[cur_arg + 1]);
	// if anyone is using a missing pheno value of -2^31, they should be
	// flogged with wet noodles
	if (scan_int32x(argv[cur_arg + 1], &missing_pheno) || (!missing_pheno) || (missing_pheno == 1) || (jj > 31) || scan_doublex(argv[cur_arg + 1], &dxx) || (dxx != (double)missing_pheno)) {
	  sprintf(g_logbuf, "Error: Invalid --missing-phenotype parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	memcpy(output_missing_pheno, argv[cur_arg + 1], jj + 1);
      } else if (!memcmp(argptr2, "issing-code", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (alloc_string(&oxford_missing_code, param_ct? argv[cur_arg + 1] : "")) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ake-pheno", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&phenoname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (alloc_string(&makepheno_str, argv[cur_arg + 2])) {
	  goto main_ret_NOMEM;
	}
	if (((argv[cur_arg + 2][0] == '\'') || (argv[cur_arg + 2][0] == '"')) && (argv[cur_arg + 2][1] == '*') && (argv[cur_arg + 2][2] == argv[cur_arg + 2][0]) && (!argv[cur_arg + 2][3])) {
	  memcpy(makepheno_str, "*", 2);
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "pheno", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &mpheno_col)) {
	  sprintf(g_logbuf, "Error: Invalid --mpheno parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "filter", 7)) {
	if (!filtername) {
	  logerrprint("Error: --mfilter must be used with --filter.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &mfilter_col)) {
	  sprintf(g_logbuf, "Error: Invalid --mfilter parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "emory", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	// may as well support systems with >2 PB RAM...
	if (scan_posintptrx(argv[cur_arg + 1], (uintptr_t*)&malloc_size_mb)) {
	  sprintf(g_logbuf, "Error: Invalid --memory parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (malloc_size_mb < BIGSTACK_MIN_MB) {
	  sprintf(g_logbuf, "Error: Invalid --memory parameter '%s' (minimum %u).\n", argv[cur_arg + 1], BIGSTACK_MIN_MB);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
#ifndef __LP64__
	if (malloc_size_mb > 2047) {
	  logerrprint("Error: --memory parameter too large for 32-bit version (max 2047).\n");
	  goto main_ret_INVALID_CMDLINE;
	}
#endif
      } else if (!memcmp(argptr2, "af", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_doublex(argv[cur_arg + 1], &min_maf)) {
	    sprintf(g_logbuf, "Error: Invalid --maf parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (min_maf <= 0.0) {
	    sprintf(g_logbuf, "Error: --maf parameter '%s' too small (must be > 0).\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  } else if (min_maf > max_maf) {
	    sprintf(g_logbuf, "Error: --maf parameter '%s' too large (must be <= %g).\n", argv[cur_arg + 1], max_maf);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else {
	  min_maf = 0.01;
	}
	filter_flags |= FILTER_ALL_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "ax-maf", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &max_maf)) {
	  sprintf(g_logbuf, "Error: Invalid --max-maf parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (max_maf < min_maf) {
	  sprintf(g_logbuf, "Error: --max-maf parameter '%s' too small (must be >= %g).\n", argv[cur_arg + 1], min_maf);
	  goto main_ret_INVALID_CMDLINE_WWA;
	} else if (max_maf >= 0.5) {
	  sprintf(g_logbuf, "Error: --max-maf parameter '%s' too large (must be < 0.5).\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	filter_flags |= FILTER_ALL_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "ind", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_doublex(argv[cur_arg + 1], &mind_thresh)) {
	    sprintf(g_logbuf, "Error: Invalid --mind parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if ((mind_thresh < 0.0) || (mind_thresh > 1.0)) {
	    sprintf(g_logbuf, "Error: Invalid --mind parameter '%s' (must be between 0 and 1).\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else {
	  mind_thresh = 0.1;
	}
	filter_flags |= FILTER_ALL_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "ake-grm", 8)) {
	logerrprint("Error: --make-grm has been retired due to inconsistent meaning across GCTA\nversions.  Use --make-grm-gz or --make-grm-bin.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ake-grm-gz", 11)) {
	if (calculation_type & CALC_RELATIONSHIP) {
	  logerrprint("Error: --make-grm-bin cannot be used with --make-grm-gz.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (distance_exp != 0.0) {
	  logerrprint("Error: \"--distance-wts exp=<x>\" cannot be used with --make-grm-gz.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        rel_info.modifier |= REL_CALC_GZ | REL_CALC_GRM;
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      logerrprint("Error: --make-grm-gz 'cov' modifier cannot coexist with --ibc flag.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (rel_info.ibc_type) {
	      logerrprint("Error: --make-grm-gz 'cov' modifier cannot coexist with an IBC modifier.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    rel_info.modifier |= REL_CALC_COV;
	  } else if (!strcmp(argv[cur_arg + uii], "no-gz")) {
	    rel_info.modifier &= ~REL_CALC_GZ;
	  } else if ((!strcmp(argv[cur_arg + uii], "ibc2")) || (!strcmp(argv[cur_arg + uii], "ibc3"))) {
	    if (rel_info.modifier & REL_CALC_COV) {
	      logerrprint("Error: --make-grm-gz 'cov' modifier cannot coexist with an IBC modifier.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (rel_info.ibc_type) {
	      sprintf(g_logbuf, "Error: --make-grm-gz '%s' modifier cannot coexist with another IBC modifier.\n", argv[cur_arg + uii]);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    rel_info.ibc_type = argv[cur_arg + uii][3] - '0';
	  } else if (!strcmp(argv[cur_arg + uii], "single-prec")) {
	    logerrprint("Error: --make-grm-gz 'single-prec' modifier has been retired.\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --make-grm-gz parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "ake-grm-bin", 12)) {
	if (distance_exp != 0.0) {
	  logerrprint("Error: \"--distance-wts exp=<x>\" cannot be used with --make-grm-bin.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	rel_info.modifier |= REL_CALC_GRM_BIN | REL_CALC_BIN4;
	if (param_ct) {
	  if (!strcmp(argv[cur_arg + 1], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      logerrprint("Error: --make-grm-bin 'cov' modifier cannot coexist with --ibc flag.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    rel_info.modifier |= REL_CALC_COV;
	  } else if ((!strcmp(argv[cur_arg + 1], "ibc2")) || (!strcmp(argv[cur_arg + 1], "ibc3"))) {
	    rel_info.ibc_type = argv[cur_arg + 1][3] - '0';
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --make-grm-bin parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "ake-rel", 8)) {
	if (calculation_type & CALC_RELATIONSHIP) {
	  logerrprint("Error: --make-rel cannot be used with --make-grm-gz/--make-grm-bin.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (distance_exp != 0.0) {
	  logerrprint("Error: \"--distance-wts exp=<x>\" cannot be used with --make-rel.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      logerrprint("Error: --make-rel 'cov' modifier cannot coexist with --ibc flag.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (rel_info.ibc_type) {
	      logerrprint("Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    rel_info.modifier |= REL_CALC_COV;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (rel_info.modifier & (REL_CALC_BIN | REL_CALC_BIN4)) {
	      logerrprint("Error: Conflicting --make-rel modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    rel_info.modifier |= REL_CALC_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (rel_info.modifier & (REL_CALC_GZ | REL_CALC_BIN4)) {
	      logerrprint("Error: Conflicting --make-rel modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    rel_info.modifier |= REL_CALC_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "bin4")) {
	    if (rel_info.modifier & (REL_CALC_GZ | REL_CALC_BIN)) {
	      logerrprint("Error: Conflicting --make-rel modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    rel_info.modifier |= REL_CALC_BIN4;
	  } else if (!strcmp(argv[cur_arg + uii], "square")) {
	    if ((rel_info.modifier & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      logerrprint("Error: --make-rel 'square' and 'square0' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if ((rel_info.modifier & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      logerrprint("Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    rel_info.modifier |= REL_CALC_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if ((rel_info.modifier & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      logerrprint("Error: --make-rel 'square' and 'square0' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if ((rel_info.modifier & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      logerrprint("Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    rel_info.modifier |= REL_CALC_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if ((rel_info.modifier & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      logerrprint("Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if ((rel_info.modifier & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      logerrprint("Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    rel_info.modifier |= REL_CALC_TRI;
	  } else if ((!strcmp(argv[cur_arg + uii], "ibc2")) || (!strcmp(argv[cur_arg + uii], "ibc3"))) {
	    if (rel_info.modifier & REL_CALC_COV) {
	      logerrprint("Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (rel_info.ibc_type) {
	      sprintf(g_logbuf, "Error: --make-rel '%s' modifier cannot coexist with another IBC modifier.\n", argv[cur_arg + uii]);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    rel_info.ibc_type = argv[cur_arg + uii][3] - '0';
	  } else if (!strcmp(argv[cur_arg + uii], "single-prec")) {
	    logerrprint("Error: --make-rel 'single-prec' modifier has been retired.  Use 'bin4'.\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --make-rel parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	if (!(rel_info.modifier & REL_CALC_SHAPEMASK)) {
	  rel_info.modifier |= (rel_info.modifier & (REL_CALC_BIN | REL_CALC_BIN4))? REL_CALC_SQ : REL_CALC_TRI;
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "atrix", 6)) {
	logprint("Note: --matrix flag deprecated.  Migrate to \"--distance ibs flat-missing\",\n\"--r2 square\", etc.\n");
        matrix_flag_state = 1;
	if (calculation_type & CALC_CLUSTER) {
	  calculation_type |= CALC_PLINK1_IBS_MATRIX;
	}
	goto main_param_zero;
      } else if (!memcmp(argptr2, "af-succ", 8)) {
	if (misc_flags & MISC_HET_SMALL_SAMPLE) {
	  logerrprint("Error: \"--het small-sample\" cannot be used with --maf-succ.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	misc_flags |= MISC_MAF_SUCC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ap3", 4)) {
	logprint("Note: --map3 flag unnecessary (.map file format is autodetected).\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ake-bed", 8)) {
        if (misc_flags & MISC_KEEP_AUTOCONV) {
	  logerrprint("Error: --make-bed cannot be used with --keep-autoconv.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --make-bed cannot be used with %s.\n", (load_rare == LOAD_RARE_CNV)? "a .cnv fileset" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  // the missing --out mistake is so common--I must have made it over a
	  // hundred times by now--that a custom error message is worthwhile.
	  sprintf(g_logbuf, "Error: --make-bed doesn't accept parameters.%s\n", ((param_ct == 1) && (!outname_end))? "  (Did you forget '--out'?)" : "");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	calculation_type |= CALC_MAKE_BED;
      } else if (!memcmp(argptr2, "ake-just-bim", 13)) {
	if (calculation_type & CALC_MAKE_BED) {
	  logerrprint("Error: --make-just-bim cannot be used with --make-bed.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --make-just-bim cannot be used with %s.\n", (load_rare == LOAD_RARE_CNV)? "a .cnv fileset" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	calculation_type |= CALC_MAKE_BIM;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ake-just-fam", 13)) {
	if (calculation_type & CALC_MAKE_BED) {
	  logerrprint("Error: --make-just-fam cannot be used with --make-bed.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --make-just-fam cannot be used with %s.\n", (load_rare == LOAD_RARE_CNV)? "a .cnv fileset" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	calculation_type |= CALC_MAKE_FAM;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "erge", 5)) {
	if (calculation_type & CALC_MERGE) {
	  logerrprint("Error: --merge cannot be used with --bmerge.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --merge cannot be used with %s.\n", (load_rare == LOAD_RARE_CNV)? ".cnv filesets" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (param_ct == 2) {
	  if (++jj > FNAMESIZE) {
	    logerrprint("Error: --merge .ped filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  jj = strlen(argv[cur_arg + 2]) + 1;
	  if (jj > FNAMESIZE) {
	    logerrprint("Error: --merge .map filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename2, argv[cur_arg + 2], jj);
	} else {
	  if (jj > (FNAMESIZE - 5)) {
	    logerrprint("Error: --merge filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(memcpya(mergename1, argv[cur_arg + 1], jj), ".ped", 5);
	  memcpy(memcpya(mergename2, argv[cur_arg + 1], jj), ".map", 5);
	}
	calculation_type |= CALC_MERGE;
      } else if (!memcmp(argptr2, "erge-list", 10)) {
	if (calculation_type & CALC_MERGE) {
	  logerrprint("Error: --merge-list cannot be used with --merge or --bmerge.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --merge-list cannot be used with %s.\n", (load_rare == LOAD_RARE_CNV)? ".cnv filesets" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	jj = strlen(argv[cur_arg + 1]) + 1;
	if (jj > FNAMESIZE) {
	  logerrprint("Error: --merge-list filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(mergename1, argv[cur_arg + 1], jj);
	merge_type |= MERGE_LIST;
	calculation_type |= CALC_MERGE;
      } else if (!memcmp(argptr2, "erge-mode", 10)) {
	if (!(calculation_type & CALC_MERGE)) {
	  logerrprint("Error: --merge-mode must be used with --{b}merge/--merge-list.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	cc = argv[cur_arg + 1][0];
	if ((cc < '1') || (cc > '7') || (argv[cur_arg + 1][1] != '\0')) {
          sprintf(g_logbuf, "Error: Invalid --merge-mode parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((merge_type & MERGE_LIST) && (cc > '5')) {
	  logerrprint("Error: --merge-mode 6-7 cannot be used with --merge-list.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        merge_type |= cc - '0';
      } else if (!memcmp(argptr2, "erge-equal-pos", 15)) {
	merge_type |= MERGE_EQUAL_POS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ust-have-sex", 13)) {
        sex_missing_pheno |= MUST_HAVE_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "covar", 6)) {
        if (!(calculation_type & CALC_GXE)) {
	  logerrprint("Error: --mcovar must be used with --covar and --gxe.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (gxe_mcovar > 1) {
	  logerrprint("Error: --mcovar cannot be used with a --gxe parameter.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &gxe_mcovar)) {
	  sprintf(g_logbuf, "Error: Invalid --mcovar parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        logprint("Note: --mcovar flag deprecated.  Use \"--gxe <covariate index>\".\n");
      } else if (!memcmp(argptr2, "odel", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 6)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (model_modifier & MODEL_ASSOC_FDEPR) {
	  model_modifier &= ~(MODEL_ASSOC | MODEL_ASSOC_FDEPR);
	} else if (model_modifier & MODEL_ASSOC) {
	  logerrprint("Error: --model cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "fisher")) {
	    if (model_modifier & MODEL_TRENDONLY) {
	      logerrprint("Error: --model 'fisher' and 'trend-only' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_FISHER;
	  } else if (!strcmp(argv[cur_arg + uii], "fisher-midp")) {
	    if (model_modifier & MODEL_TRENDONLY) {
	      logerrprint("Error: --model 'fisher-midp' and 'trend-only' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_FISHER | MODEL_FISHER_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (model_modifier & MODEL_MPERM) {
	      logerrprint("Error: --model 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_PERM;
	  } else if (!strcmp(argv[cur_arg + uii], "genedrop")) {
	    model_modifier |= MODEL_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    model_modifier |= MODEL_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "dom")) {
	    if (model_modifier & (MODEL_PREC | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logerrprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_PDOM;
	  } else if (!strcmp(argv[cur_arg + uii], "rec")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logerrprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_PREC;
	  } else if (!strcmp(argv[cur_arg + uii], "gen")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logerrprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (mtest_adjust) {
	      logerrprint("Error: --model perm-gen cannot be used with --adjust.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_PGEN;
	  } else if (!strcmp(argv[cur_arg + uii], "trend")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	      logerrprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_PTREND;
	  } else if (!strcmp(argv[cur_arg + uii], "trend-only")) {
	    if (model_modifier & (MODEL_FISHER | MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	      logerrprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    model_modifier |= MODEL_PTREND | MODEL_TRENDONLY;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (model_modifier & MODEL_PERM) {
	      logerrprint("Error: --model 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (model_modifier & MODEL_MPERM) {
	      logerrprint("Error: Duplicate --model 'mperm' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (scan_posint_defcapx(&(argv[cur_arg + uii][6]), &model_mperm_val)) {
	      sprintf(g_logbuf, "Error: Invalid --model mperm parameter '%s'.\n", &(argv[cur_arg + uii][6]));
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    model_modifier |= MODEL_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    logerrprint("Error: Improper --model mperm syntax.  (Use \"--model mperm=<value>\".)\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
	    model_modifier |= MODEL_SET_TEST;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --model parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_MODEL;
      } else if (!memcmp(argptr2, "odel-dom", 9)) {
	if (model_modifier & MODEL_ASSOC) {
	  logerrprint("Error: --model-dom cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PREC | MODEL_PGEN | MODEL_PTREND)) {
	  logerrprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --model-dom flag deprecated.  Use \"--model dom\".\n");
	model_modifier |= MODEL_PDOM;
	calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "odel-gen", 9)) {
	if (model_modifier & MODEL_ASSOC) {
	  logerrprint("Error: --model-gen cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND)) {
	  logerrprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --model-gen flag deprecated.  Use \"--model gen\".\n");
	model_modifier |= MODEL_PGEN;
        calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "odel-rec", 9)) {
	if (model_modifier & MODEL_ASSOC) {
	  logerrprint("Error: --model-rec cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PTREND)) {
	  logerrprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --model-rec flag deprecated.  Use \"--model rec\".\n");
	model_modifier |= MODEL_PREC;
        calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "odel-trend", 11)) {
	if (model_modifier & MODEL_ASSOC) {
	  logerrprint("Error: --model-trend cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PREC)) {
	  logerrprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --model-trend flag deprecated.  Use \"--model trend\".\n");
	model_modifier |= MODEL_PTREND;
        calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "perm", 5)) {
	if (model_modifier & (MODEL_PERM | MODEL_MPERM)) {
	  sprintf(g_logbuf, "Error: --mperm cannot be used with --%s %sperm.\n", (model_modifier & MODEL_ASSOC)? "assoc" : "model", (model_modifier & MODEL_PERM)? "" : "m");
	  goto main_ret_INVALID_CMDLINE_2A;
	} else if (glm_modifier & (GLM_PERM | GLM_MPERM)) {
	  sprintf(g_logbuf, "Error: --mperm cannot be used with --%s %sperm.\n", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", (glm_modifier & GLM_PERM)? "" : "m");
	  goto main_ret_INVALID_CMDLINE_2A;
	} else if (family_info.dfam_modifier & (DFAM_PERM | DFAM_MPERM)) {
	  sprintf(g_logbuf, "Error: --mperm cannot be used with --dfam %sperm.\n", (family_info.dfam_modifier & DFAM_PERM)? "" : "m");
	  goto main_ret_INVALID_CMDLINE_2A;
	} else if (cluster.modifier & (CLUSTER_CMH_PERM | CLUSTER_CMH_MPERM)) {
	  sprintf(g_logbuf, "Error: --mperm cannot be used with --%s %sperm.\n", (cluster.modifier & CLUSTER_CMH_BD)? "bd" : "mh", (cluster.modifier & CLUSTER_CMH_PERM)? "" : "m");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &mperm_val)) {
	  sprintf(g_logbuf, "Error: Invalid --mperm parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (load_rare & LOAD_RARE_CNV) {
	  if ((cnv_calc_type & CNV_SAMPLE_PERM) && (!cnv_sample_mperms)) {
	    logprint("Note: --mperm flag deprecated.  Use e.g. \"--cnv-indiv-perm <perm. count>\".\n");
	    cnv_sample_mperms = mperm_val;
	  } else if ((cnv_calc_type & CNV_TEST_REGION) && (!cnv_test_region_mperms)) {
	    logprint("Note: --mperm flag deprecated.  Use e.g. \"--cnv-test-region <perm. count>\".\n");
	  } else if ((cnv_calc_type & CNV_ENRICHMENT_TEST) && (!cnv_enrichment_test_mperms)) {
	    logprint("Note: --mperm flag deprecated.  Use e.g. \"--cnv-enrichment-test <perm. count>\".\n");
	  } else {
	    logprint("Note: --mperm flag deprecated.  Use e.g. \"--cnv-test <permutation count>\".\n");
            if (!(cnv_calc_type & (CNV_SAMPLE_PERM | CNV_ENRICHMENT_TEST | CNV_TEST | CNV_TEST_REGION))) {
	      cnv_calc_type |= CNV_TEST;
	    }
	    cnv_test_mperms = mperm_val;
	  }
	  // if e.g. --cnv-test-region had a valid parameter, don't clobber it
	  if (!cnv_test_region_mperms) {
	    cnv_test_region_mperms = mperm_val;
	  }
	  if (!cnv_enrichment_test_mperms) {
	    cnv_enrichment_test_mperms = mperm_val;
	  }
	} else {
	  logprint("Note: --mperm flag deprecated.  Use e.g. \"--model mperm=<value>\".\n");
	  model_mperm_val = mperm_val;
	  model_modifier |= MODEL_MPERM;
	  glm_mperm_val = mperm_val;
	  glm_modifier |= GLM_MPERM;
          testmiss_mperm_val = mperm_val;
          testmiss_modifier |= TESTMISS_MPERM;
	  family_info.tdt_mperm_val = mperm_val;
	  family_info.tdt_modifier |= TDT_MPERM;
	  family_info.dfam_mperm_val = mperm_val;
	  family_info.dfam_modifier |= DFAM_MPERM;
	  family_info.qfam_mperm_val = mperm_val;
	  family_info.qfam_modifier |= QFAM_MPERM;
          cluster.cmh_mperm_val = mperm_val;
	  cluster.modifier |= CLUSTER_CMH_MPERM;
	}
      } else if (!memcmp(argptr2, "perm-save", 10)) {
	if (glm_modifier & GLM_NO_SNP) {
          logerrprint("Error: --mperm-save cannot be used with --linear/--logistic no-snp.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	mperm_save |= MPERM_DUMP_BEST;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "perm-save-all", 14)) {
	mperm_save |= MPERM_DUMP_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "c", 2)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  logerrprint("Error: --mc must be used with --cluster.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &cluster.max_size) || (cluster.max_size == 1)) {
	  sprintf(g_logbuf, "Error: Invalid --mc parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "cc", 2)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  logerrprint("Error: --mcc must be used with --cluster.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &cluster.max_cases)) {
	  sprintf(g_logbuf, "Error: Invalid --mcc parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (cluster.max_cases > cluster.max_size) {
          logerrprint("Error: --mcc parameter exceeds --mc parameter.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 2], &cluster.max_ctrls)) {
	  sprintf(g_logbuf, "Error: Invalid --mcc parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (cluster.max_ctrls > cluster.max_size) {
          logerrprint("Error: --mcc parameter exceeds --mc parameter.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
      } else if (!memcmp(argptr2, "atch", 5)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  logerrprint("Error: --match must be used with --cluster.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&cluster.match_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        if (param_ct == 2) {
	  if (alloc_string(&cluster.match_missing_str, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	}
      } else if (!memcmp(argptr2, "atch-type", 10)) {
	if (!cluster.match_fname) {
	  logerrprint("Error: --match-type must be used with --match.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&cluster.match_type_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ds-plot", 8)) {
#ifdef NOLAPACK
	// may want to support this again, now that we've switched the default
	// algorithm back to PLINK 1.07-style SVD
        logerrprint("Error: --mds-plot requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	goto main_ret_INVALID_CMDLINE;
#else
	if (!(calculation_type & CALC_CLUSTER)) {
	  logerrprint("Error: --mds-plot must be used with --cluster.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	cluster.mds_dim_ct = 0;
        for (uii = 1; uii <= param_ct; uii++) {
          if (!strcmp(argv[cur_arg + uii], "by-cluster")) {
	    cluster.modifier |= CLUSTER_MDS;
	  } else if (!strcmp(argv[cur_arg + uii], "eigendecomp")) {
	    cluster.modifier |= CLUSTER_MDS_EIGENDECOMP;
	  } else if (!strcmp(argv[cur_arg + uii], "eigvals")) {
	    cluster.modifier |= CLUSTER_MDS_EIGVALS;
	  } else {
	    if (cluster.mds_dim_ct) {
	      logerrprint("Error: Invalid --mds-plot parameter sequence.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (scan_posint_defcapx(argv[cur_arg + uii], &cluster.mds_dim_ct)) {
	      sprintf(g_logbuf, "Error: Invalid --mds-plot parameter '%s'.\n", argv[cur_arg + uii]);
              goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	}
	if (!cluster.mds_dim_ct) {
	  logerrprint("Error: Missing --mds-plot dimension count.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
#endif
      } else if (!memcmp(argptr2, "ds-cluster", 11)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  logerrprint("Error: --mds-cluster must be used with --cluster.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        logprint("Note: --mds-cluster flag deprecated.  Use \"--mds-plot by-cluster\".\n");
        cluster.modifier |= CLUSTER_MDS;
      } else if (!memcmp(argptr2, "within", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &mwithin_col)) {
	  sprintf(g_logbuf, "Error: Invalid --mwithin parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "in", 3)) {
        if (!(calculation_type & CALC_GENOME)) {
	  logerrprint("Error: --min must be used with --genome.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx)) {
	  sprintf(g_logbuf, "Error: Invalid --min parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((dxx < -1.0) || (dxx > 1.0)) {
          logerrprint("Error: --min threshold must be between -1 and 1 inclusive.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (dxx > genome_max_pi_hat) {
	  logerrprint("Error: --min value cannot be greater than --max value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        genome_modifier |= GENOME_FILTER_PI_HAT;
	genome_min_pi_hat = dxx;
      } else if (!memcmp(argptr2, "ax", 3)) {
        if (!(calculation_type & CALC_GENOME)) {
	  logerrprint("Error: --max must be used with --genome.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx)) {
	  sprintf(g_logbuf, "Error: Invalid --max parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((dxx < -1.0) || (dxx > 1.0)) {
          logerrprint("Error: --max threshold must be between -1 and 1 inclusive.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	genome_modifier |= GENOME_FILTER_PI_HAT;
	genome_max_pi_hat = dxx;
      } else if (!memcmp(argptr2, "ake-founders", 13)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "require-2-missing")) {
	    misc_flags |= MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING;
	  } else if (!strcmp(argv[cur_arg + uii], "first")) {
	    misc_flags |= MISC_MAKE_FOUNDERS_FIRST;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --make-founders parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	filter_flags |= FILTER_FAM_REQ | FILTER_MAKE_FOUNDERS;
      } else if (!memcmp(argptr2, "issing", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "gz")) {
	    sprintf(g_logbuf, "Error: Invalid --missing parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  misc_flags |= MISC_MISSING_GZ;
	}
	calculation_type |= CALC_MISSING_REPORT;
      } else if (!memcmp(argptr2, "h", 2)) {
	if (calculation_type & CALC_CMH) {
	  logerrprint("Error: --mh is redundant with --bd.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (cluster.modifier & CLUSTER_CMH_MPERM) {
	      logerrprint("Error: --mh 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    cluster.modifier |= CLUSTER_CMH_PERM;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (cluster.modifier & CLUSTER_CMH_PERM) {
	      logerrprint("Error: --mh 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (cluster.modifier & CLUSTER_CMH_MPERM) {
	      logerrprint("Error: Duplicate --mh 'mperm' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (scan_posint_defcapx(&(argv[cur_arg + uii][6]), &(cluster.cmh_mperm_val))) {
	      sprintf(g_logbuf, "Error: Invalid --mh mperm parameter '%s'.\n", &(argv[cur_arg + uii][6]));
              goto main_ret_INVALID_CMDLINE_WWA;
	    }
            cluster.modifier |= CLUSTER_CMH_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
            cluster.modifier |= CLUSTER_CMH_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
	    cluster.modifier |= CLUSTER_CMH_SET_TEST;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
            logerrprint("Error: Improper --mh mperm syntax.  (Use \"--mh mperm=<value>\".)\n");
            goto main_ret_INVALID_CMDLINE_A;
	  } else {
            sprintf(g_logbuf, "Error: Invalid --mh parameter '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_CMH;
      } else if (!memcmp(argptr2, "h2", 3)) {
	if (calculation_type & CALC_CMH) {
	  logerrprint("Error: --mh2 cannot be used with --mh/--bd.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	calculation_type |= CALC_CMH;
	cluster.modifier |= CLUSTER_CMH2;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ake-set", 8)) {
	if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --make-set cannot be used with %s.\n", (load_rare == LOAD_RARE_CNV)? "a .cnv fileset" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&set_info.fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        set_info.modifier |= SET_MAKE_FROM_RANGES;
	filter_flags |= FILTER_BIM_REQ;
      } else if (!memcmp(argptr2, "ake-set-border", 15)) {
	if (!set_info.fname) {
	  logerrprint("Error: --make-set-border must be used with --make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --make-set-border parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (dxx > 2147483.646) {
	  set_info.make_set_border = 2147483646;
	} else {
	  set_info.make_set_border = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
	}
      } else if (!memcmp(argptr2, "ake-set-collapse-group", 23)) {
        if (!set_info.fname) {
	  logerrprint("Error: --make-set-collapse-group must be used with --make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        set_info.modifier |= SET_MAKE_COLLAPSE_GROUP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ake-set-complement-all", 23)) {
	if (set_info.modifier & SET_COMPLEMENTS) {
	  logerrprint("Error: --make-set-complement-all cannot be used with --complement-sets.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&set_info.merged_set_name, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
	set_info.modifier |= SET_COMPLEMENTS;
      } else if (!memcmp(argptr2, "ake-set-complement-group", 25)) {
        if (!set_info.fname) {
	  logerrprint("Error: --make-set-complement-group must be used with --make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (set_info.modifier & (SET_COMPLEMENTS | SET_MAKE_COLLAPSE_GROUP)) {
	  logerrprint("Error: --make-set-complement-group cannot be used with --complement-sets,\n--make-set-collapse-group, or --make-set-complement-all.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        set_info.modifier |= SET_COMPLEMENTS | SET_C_PREFIX | SET_MAKE_COLLAPSE_GROUP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "erge-x", 7)) {
	if ((chrom_info.xymt_codes[X_OFFSET] == -2) || (chrom_info.xymt_codes[XY_OFFSET] == -2)) {
	  logerrprint("Error: --merge-x must be used with a chromosome set containing X and XY codes.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct == 1) {
	  if (strcmp(argv[cur_arg + 1], "no-fail")) {
	    sprintf(g_logbuf, "Error: Invalid --merge-x parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  misc_flags |= MISC_SPLIT_MERGE_NOFAIL;
	}
	misc_flags |= MISC_MERGEX;
      } else if (!memcmp(argptr2, "issing-var-code", 16)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (alloc_string(&missing_marker_id_match, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "e", 2)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	ujj = 2;
	if (param_ct == 3) {
	  if (!strcmp(argv[cur_arg + 1], "var-first")) {
	    uii = 2;
	    ujj = 3;
	  } else if (!strcmp(argv[cur_arg + 2], "var-first")) {
	    ujj = 3;
	  } else if (strcmp(argv[cur_arg + 3], "var-first")) {
	    sprintf(g_logbuf, "Error: Invalid --me parameter '%s'.\n", argv[cur_arg + 3]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  family_info.mendel_modifier |= MENDEL_FILTER_VAR_FIRST;
	}
	if (scan_doublex(argv[cur_arg + uii], &family_info.mendel_max_trio_error) || (family_info.mendel_max_trio_error < 0.0) || (family_info.mendel_max_trio_error > 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --me maximum per-trio error rate '%s'.\n", argv[cur_arg + uii]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (scan_doublex(argv[cur_arg + ujj], &family_info.mendel_max_var_error) || (family_info.mendel_max_var_error < 0.0) || (family_info.mendel_max_var_error > 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --me maximum per-variant error rate '%s'.\n", argv[cur_arg + ujj]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((family_info.mendel_max_trio_error < 1.0) || (family_info.mendel_max_var_error < 1.0)) {
	  if (covar_fname && (family_info.mendel_max_trio_error < 1.0)) {
	    // corner case: --me screws up covariate filtered indices
	    logerrprint("Error: --covar cannot be used with --me.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  // silently skip if both parameters are one, for backward
	  // compatibility
	  family_info.mendel_modifier |= MENDEL_FILTER;
	  filter_flags |= FILTER_ALL_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
	}
      } else if (!memcmp(argptr2, "e-exclude-one", 14)) {
	if (!(family_info.mendel_modifier & MENDEL_FILTER)) {
	  logerrprint("Error: --me-exclude-one must be used with a --me parameter smaller than one.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (param_ct) {
	  if (scan_doublex(argv[cur_arg + 1], &family_info.mendel_exclude_one_ratio) || (family_info.mendel_exclude_one_ratio < 1.0)) {
	    sprintf(g_logbuf, "Error: Invalid --me-exclude-one ratio '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else {
	  family_info.mendel_exclude_one_ratio = -1;
	}
      } else if (!memcmp(argptr2, "at", 3)) {
	logprint("Note: --mat flag deprecated.  Use e.g. \"--tdt poo mperm=<value> mat\".\n");
	family_info.tdt_modifier |= TDT_POOPERM_MAT;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "endel", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "summaries-only")) {
	    sprintf(g_logbuf, "Error: Invalid --mendel parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE;
	  }
	  family_info.mendel_modifier |= MENDEL_SUMMARIES_ONLY;
	}
	calculation_type |= CALC_MENDEL;
      } else if (!memcmp(argptr2, "endel-duos", 11)) {
	family_info.mendel_modifier |= MENDEL_DUOS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "endel-multigen", 15)) {
	family_info.mendel_modifier |= MENDEL_MULTIGEN;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ake-perm-pheno", 15)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &permphe_ct)) {
	  sprintf(g_logbuf, "Error: Invalid --make-perm-pheno parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	calculation_type |= CALC_MAKE_PERM_PHENO;
      } else if (!memcmp(argptr2, "eta-analysis", 13)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 0x1fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	// if '+' present, must detect it before using alloc_and_flatten()
        // must detect '+'
	for (uii = 1; uii <= param_ct; uii++) {
	  if (argv[cur_arg + uii][0] == '+') {
            if (!argv[cur_arg + uii][1]) {
	      if (uii <= 2) {
	        logerrprint("Error: --meta-analysis requires at least two PLINK report files.\n");
	        goto main_ret_INVALID_CMDLINE_A;
	      }
	      break;
            } else {
              const char* post_plus = &(argv[cur_arg + uii][1]);
              if ((!strcmp(post_plus, "study")) ||
                  (!strcmp(post_plus, "no-map")) ||
                  (!strcmp(post_plus, "no-allele")) ||
                  (!strcmp(post_plus, "report-all")) ||
                  (!strcmp(post_plus, "logscale")) ||
                  (!strcmp(post_plus, "qt")) ||
                  (!strcmp(post_plus, "weighted-z"))) {
                // This is far more likely to be a mistake than a real
                // filename...
                LOGERRPRINTFWW("Warning: In order to end the filename list, --meta-analysis '+' must be followed by a space.  Currently interpreting %s as a literal filename, which probably isn't what you want...\n", argv[cur_arg + uii]);
              }
            }
	  }
	}
	retval = alloc_and_flatten(&metaanal_fnames, &(argv[cur_arg + 1]), uii - 1);
	if (retval) {
	  goto main_ret_1;
	}
	for (uii++; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "study")) {
	    metaanal_flags |= METAANAL_STUDY;
	  } else if (!strcmp(argv[cur_arg + uii], "no-map")) {
	    metaanal_flags |= METAANAL_NO_MAP | METAANAL_NO_ALLELE;
	  } else if (!strcmp(argv[cur_arg + uii], "no-allele")) {
	    metaanal_flags |= METAANAL_NO_ALLELE;
	  } else if (!strcmp(argv[cur_arg + uii], "report-all")) {
	    metaanal_flags |= METAANAL_REPORT_ALL;
	  } else if (!strcmp(argv[cur_arg + uii], "logscale")) {
	    metaanal_flags |= METAANAL_LOGSCALE;
	  } else if (!strcmp(argv[cur_arg + uii], "qt")) {
	    metaanal_flags |= METAANAL_QT | METAANAL_LOGSCALE;
	  } else if (!strcmp(argv[cur_arg + uii], "weighted-z")) {
	    metaanal_flags |= METAANAL_WEIGHTED_Z;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --meta-analysis parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
      } else if (!memcmp(argptr2, "eta-analysis-a1-field", 22)) {
        if (!metaanal_fnames) {
	  logerrprint("Error: --meta-analysis-a1-field must be used with --meta-analysis.\n");
          goto main_ret_INVALID_CMDLINE;
	} else if (metaanal_flags & METAANAL_NO_ALLELE) {
	  logerrprint("Error: --meta-analysis-a1-field cannot be used with --meta-analysis\n'no-map'/'no-allele'.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x10000000)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&metaanal_a1field_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "eta-analysis-a2-field", 22)) {
        if (!metaanal_fnames) {
	  logerrprint("Error: --meta-analysis-a2-field must be used with --meta-analysis.\n");
          goto main_ret_INVALID_CMDLINE;
	} else if (metaanal_flags & METAANAL_NO_ALLELE) {
	  logerrprint("Error: --meta-analysis-a2-field cannot be used with --meta-analysis\n'no-map'/'no-allele'.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x10000000)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&metaanal_a2field_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "eta-analysis-chr-field", 23)) {
        if (!metaanal_fnames) {
	  logerrprint("Error: --meta-analysis-chr-field must be used with --meta-analysis.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x10000000)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&metaanal_chrfield_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "eta-analysis-snp-field", 23)) {
        if (!metaanal_fnames) {
	  logerrprint("Error: --meta-analysis-snp-field must be used with --meta-analysis.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x10000000)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&metaanal_snpfield_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "eta-analysis-bp-field", 22)) {
        if (!metaanal_fnames) {
	  logerrprint("Error: --meta-analysis-bp-field must be used with --meta-analysis.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x10000000)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&metaanal_bpfield_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "eta-analysis-p-field", 21)) {
        if ((!metaanal_fnames) || (!(metaanal_flags & METAANAL_WEIGHTED_Z))) {
	  logerrprint("Error: --meta-analysis-p-field must be used with --meta-analysis + weighted-z.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x10000000)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&metaanal_pfield_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "eta-analysis-se-field", 22)) {
        if (!metaanal_fnames) {
	  logerrprint("Error: --meta-analysis-se-field must be used with --meta-analysis.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x10000000)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&metaanal_sefield_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "eta-analysis-ess-field", 23)) {
        if ((!metaanal_fnames) || (!(metaanal_flags & METAANAL_WEIGHTED_Z))) {
	  logerrprint("Error: --meta-analysis-ess-field must be used with --meta-analysis + weighted-z.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x10000000)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&metaanal_essfield_search_order, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "eta-analysis-report-dups", 25)) {
	if (!metaanal_fnames) {
	  logerrprint("Error: --meta-analysis-report-dups must be used with --meta-analysis.\n");
	}
	metaanal_flags |= METAANAL_REPORT_DUPS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ac", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &min_ac)) {
	  sprintf(g_logbuf, "Error: Invalid --mac parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	filter_flags |= FILTER_ALL_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "ax-mac", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &max_ac)) {
	  sprintf(g_logbuf, "Error: Invalid --max-mac parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        if (max_ac < min_ac) {
	  logerrprint("Error: --max-mac parameter cannot be smaller than --mac parameter.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_flags |= FILTER_ALL_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "lma", 4)) {
        logerrprint("Error: --mlma is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "lma-loco", 9)) {
        logerrprint("Error: --mlma-loco is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "lma-no-adj-covar", 17)) {
        logerrprint("Error: --mlma-no-adj-covar is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ishap-window", 13)) {
        logerrprint("Error: --mishap-window is provisionally retired.  Contact the developers if you\nneed this function.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if ((!memcmp(argptr2, "ake-bpgen", 10)) || (!memcmp(argptr2, "ake-pgen", 9))) {
	sprintf(g_logbuf, "Error: Unrecognized flag (%s).  (This is PLINK 1.9, not 2.x.)\n", argv[cur_arg]);
	goto main_ret_INVALID_CMDLINE_2;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'n':
      if (!memcmp(argptr2, "o-fid", 6)) {
	fam_cols &= ~FAM_COL_1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-parents", 10)) {
	fam_cols &= ~FAM_COL_34;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-sex", 6)) {
	if (filter_flags & (FILTER_BINARY_FEMALES | FILTER_BINARY_MALES)) {
	  logerrprint("Error: --filter-males/--filter-females cannot be used with --no-sex.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	fam_cols &= ~FAM_COL_5;
	sex_missing_pheno |= ALLOW_NO_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-pheno", 8)) {
	fam_cols &= ~FAM_COL_6;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "onfounders", 11)) {
	misc_flags |= MISC_NONFOUNDERS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "eighbour", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &neighbor_n1)) {
	  sprintf(g_logbuf, "Error: Invalid --neighbour parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (scan_posint_defcapx(argv[cur_arg + 2], &neighbor_n2)) {
	  sprintf(g_logbuf, "Error: Invalid --neighbour parameter '%s'.\n", argv[cur_arg + 2]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (neighbor_n2 < neighbor_n1) {
	  logerrprint("Error: Second --neighbour parameter cannot be smaller than first parameter.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        calculation_type |= CALC_NEIGHBOR;
      } else if (!memcmp(argptr2, "ot-chr", 7)) {
	if (markername_from) {
	  logerrprint("Error: --from cannot be used with --autosome[-xy] or --[not-]chr.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	// allowed:
	//   --allow-extra-chr --chr 5-22 bobs_chrom --not-chr 17
	// allowed:
	//   --allow-extra-chr --not-chr 12-17 bobs_chrom
	// does not make sense, disallowed:
	//   --allow-extra-chr --chr 5-22 --not-chr bobs_chrom

	// --allow-extra-chr present, --chr/--autosome[-xy] not present
	uii = ((misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1) && (!chrom_info.is_include_stack);
	retval = parse_chrom_ranges(param_ct, '-', &(argv[cur_arg]), chrom_exclude, &chrom_info, uii, argptr);
	if (retval) {
	  goto main_ret_1;
	}
	if (chrom_info.is_include_stack || (!chrom_flag_present)) {
	  init_default_chrom_mask(&chrom_info);
	}
	for (uii = 0; uii < CHROM_MASK_INITIAL_WORDS; uii++) {
	  chrom_info.chrom_mask[uii] &= ~chrom_exclude[uii];
	}
	if (are_all_words_zero(chrom_info.chrom_mask, CHROM_MASK_INITIAL_WORDS) && ((!((misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1)) || (chrom_info.is_include_stack && (!chrom_info.incl_excl_name_stack)))) {
	  logerrprint("Error: All chromosomes excluded.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	chrom_flag_present = 1;
      } else if (!memcmp(argptr2, "udge", 5)) {
        if (!(calculation_type & CALC_GENOME)) {
	  logerrprint("Error: --nudge must be used with --genome.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        logprint("Note: --nudge flag deprecated.  Use \"--genome nudge\".\n");
        genome_modifier |= GENOME_NUDGE;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "ew-id-max-allele-len", 21)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_cappedx(argv[cur_arg + 1], MAX_ID_SLEN - 2, &new_id_max_allele_len)) {
	  sprintf(g_logbuf, "Error: Invalid --new-id-max-allele-len parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "o-snp", 6)) {
	if (!(calculation_type & CALC_GLM)) {
	  logerrprint("Error: --no-snp must be used with --linear or --logistic.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (mtest_adjust) {
	  logerrprint("Error: --no-snp cannot be used with --adjust.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (mperm_save & MPERM_DUMP_BEST) {
	  logerrprint("Error: --no-snp cannot be used with --mperm-save.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if ((glm_modifier & (GLM_NO_SNP_EXCL - GLM_HETHOM - GLM_DOMINANT)) || ((glm_modifier & (GLM_HETHOM | GLM_DOMINANT)) && (!(glm_modifier & (GLM_CONDITION_DOMINANT | GLM_CONDITION_RECESSIVE))))) {
	  sprintf(g_logbuf, "Error: --no-snp conflicts with a --%s modifier.\n", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	logprint("Note: --no-snp flag deprecated.  Use e.g. \"--linear no-snp\".\n");
        glm_modifier |= GLM_NO_SNP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-x-sex", 8)) {
	if (!(calculation_type & CALC_GLM)) {
	  logerrprint("Error: --no-x-sex must be used with --linear or --logistic.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (glm_modifier & (GLM_NO_SNP | GLM_SEX)) {
	  sprintf(g_logbuf, "Error: --no-x-sex conflicts with a --%s modifier.\n", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	logprint("Note: --no-x-sex flag deprecated.  Use e.g. \"--linear no-x-sex\".\n");
	glm_modifier |= GLM_NO_X_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "op", 3)) {
	if (!(epi_info.modifier & EPI_FAST)) {
	  logerrprint("Error: --nop must be used with --fast-epistasis.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --nop flag deprecated.  Use \"--fast-epistasis nop\".\n");
        epi_info.modifier |= EPI_FAST_NO_P_VALUE;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "o-const-covar", 14)) {
	if (!covar_fname) {
	  logerrprint("Error: --no-const-covar must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	covar_modifier |= COVAR_NO_CONST;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "oweb", 5)) {
        logprint("Note: --noweb has no effect since no web check is implemented yet.\n");
	goto main_param_zero;
      } else {
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'o':
      if (!memcmp(argptr2, "utput-chr", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (!strcmp(argv[cur_arg + 1], "M")) {
          chrom_info.output_encoding = CHR_OUTPUT_M;
	} else if (!strcmp(argv[cur_arg + 1], "MT")) {
          chrom_info.output_encoding = CHR_OUTPUT_MT;
	} else if (!strcmp(argv[cur_arg + 1], "0M")) {
	  chrom_info.output_encoding = CHR_OUTPUT_0M;
	} else if (!strcmp(argv[cur_arg + 1], "chr26")) {
          chrom_info.output_encoding = CHR_OUTPUT_PREFIX;
	} else if (!strcmp(argv[cur_arg + 1], "chrM")) {
          chrom_info.output_encoding = CHR_OUTPUT_PREFIX | CHR_OUTPUT_M;
	} else if (!strcmp(argv[cur_arg + 1], "chrMT")) {
          chrom_info.output_encoding = CHR_OUTPUT_PREFIX | CHR_OUTPUT_MT;
	} else if (strcmp(argv[cur_arg + 1], "26")) {
	  sprintf(g_logbuf, "Error: Invalid --output-chr parameter '%s'.\n", argv[cur_arg + 1]);
          goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "utput-missing-genotype", 23)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	cc = extract_char_param(argv[cur_arg + 1]);
	if (((unsigned char)cc) <= ' ') {
	  sprintf(g_logbuf, "Error: Invalid --output-missing-genotype parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	g_output_missing_geno_ptr = &(g_one_char_strs[((unsigned char)cc) * 2]);
      } else if (!memcmp(argptr2, "utput-missing-phenotype", 24)) {
	if (load_rare == LOAD_RARE_DOSAGE) {
	  logerrprint("Error: --output-missing-phenotype has no effect with --dosage.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > 31) {
	  logerrprint("Error: --output-missing-phenotype string too long (max 31 chars).\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	memcpy(output_missing_pheno, argv[cur_arg + 1], jj + 1);
      } else if (!memcmp(argptr2, "blig-clusters", 14)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&oblig_missing_info.sample_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	logprint("Note: --oblig-clusters flag deprecated.  Use just --oblig-missing.\n");
      } else if (!memcmp(argptr2, "blig-missing", 13)) {
	if ((geno_thresh == 1.0) && (mind_thresh == 1.0) && (!(calculation_type & CALC_MISSING_REPORT))) {
	  logerrprint("Error: --oblig-missing must be used with --geno, --mind, or --missing.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (!oblig_missing_info.sample_fname) {
          if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	} else if (param_ct != 1) {
          logerrprint("Error: --oblig-missing requires exactly one parameter when --oblig-clusters is\nalso present.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	retval = alloc_fname(&oblig_missing_info.marker_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (param_ct == 2) {
	  retval = alloc_fname(&oblig_missing_info.sample_fname, argv[cur_arg + 2], argptr, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "xford-single-chr", 17)) {
	if (!(load_params & LOAD_PARAMS_OXGEN)) {
	  if (load_params & LOAD_PARAMS_OXBGEN) {
	    logerrprint("Error: --oxford-single-chr must be used with .gen input.  (Single-chromosome\n.bgen files do not require this, since they still contain chromosome codes.)\n");
	  } else {
	    logerrprint("Error: --oxford-single-chr must be used with .gen input.\n");
	  }
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (!(misc_flags & MISC_ALLOW_EXTRA_CHROMS)) {
	  if (get_chrom_code_raw(argv[cur_arg + 1]) < 0) {
	    sprintf(g_logbuf, "Error: Invalid --oxford-single-chr chromosome code '%s'. (Did you forget --allow-extra-chr?)\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        if (alloc_string(&oxford_single_chr, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "xford-pheno-name", 17)) {
	if (!(load_params & LOAD_PARAMS_OX_ALL)) {
	  logerrprint("Error: --oxford-pheno-name must be used with an Oxford-format fileset.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (alloc_string(&oxford_pheno_name, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "utput-min-p", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &output_min_p) || (!(output_min_p >= 0.0)) || (output_min_p >= 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --output-min-p parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (memcmp(argptr2, "ut", 3)) {
	// --out is a special case due to logging
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'p':
      if (!memcmp(argptr2, "ed", 3)) {
	if ((load_params & (LOAD_PARAMS_BFILE_ALL | LOAD_PARAMS_OX_ALL)) || load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	load_params |= LOAD_PARAMS_PED;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logerrprint("Error: --ped parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	} else if (!memcmp(argv[cur_arg + 1], "-", 2)) {
	  logerrprint("Error: \"--ped -\" is no longer supported.  Redirect to a temporary file and load\nit the usual way, or use PLINK 1.07.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "heno", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (makepheno_str) {
	  logerrprint("Error: --pheno and --make-pheno flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&phenoname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "heno-name", 10)) {
	if (!phenoname) {
	  logerrprint("Error: --pheno-name must be used with --pheno.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (mpheno_col != 0) {
	  logerrprint("Error: --mpheno and --pheno-name flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (alloc_string(&phenoname_str, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "heno-merge", 11)) {
	if (!phenoname) {
	  logerrprint("Error: --pheno-merge must be used with --pheno.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	pheno_modifier |= PHENO_MERGE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rune", 5)) {
	filter_flags |= FILTER_FAM_REQ | FILTER_PRUNE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "arallel", 8)) {
	if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	  logerrprint("Error: --parallel cannot be used with \"--distance square\".  Use \"--distance\nsquare0\" or plain --distance instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if ((dist_calc_type & (DISTANCE_BIN | DISTANCE_BIN4)) && (!(dist_calc_type & DISTANCE_SHAPEMASK))) {
	  logerrprint("Error: --parallel cannot be used with plain \"--distance bin[4]\".  Use e.g.\n\"--distance bin square0\" or \"--distance bin triangle\" instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if ((rel_info.modifier & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	  logerrprint("Error: --parallel cannot be used with \"--make-rel square\".  Use \"--make-rel\nsquare0\" or plain --make-rel instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if ((rel_info.modifier & (REL_CALC_BIN | REL_CALC_BIN4)) && (!(rel_info.modifier & (REL_CALC_SHAPEMASK | REL_CALC_GRM_BIN)))) {
	  logerrprint("Error: --parallel cannot be used with plain \"--make-rel bin[4]\".  Use e.g.\n\"--make-rel bin square0\" or \"--make-rel bin triangle\" instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (calculation_type & CALC_PLINK1_DISTANCE_MATRIX) {
	  logerrprint("Error: --parallel and --distance-matrix cannot be used together.  Use\n--distance instead.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (calculation_type & CALC_GROUPDIST) {
	  logerrprint("Error: --parallel and --groupdist cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (calculation_type & CALC_CLUSTER) {
	  logerrprint("Error: --parallel and --cluster cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (calculation_type & CALC_NEIGHBOR) {
	  logerrprint("Error: --parallel and --neighbour cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_cappedx(argv[cur_arg + 1], PARALLEL_MAX, &parallel_idx)) {
	  sprintf(g_logbuf, "Error: Invalid --parallel job index '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (scan_posint_cappedx(argv[cur_arg + 2], PARALLEL_MAX, &parallel_tot) || (parallel_tot == 1) || (parallel_tot < parallel_idx)) {
	  sprintf(g_logbuf, "Error: Invalid --parallel total job count '%s'.\n", argv[cur_arg + 2]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	parallel_idx--; // internal 0..(n-1) indexing
      } else if (!memcmp(argptr2, "pc-gap", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx)) {
	  sprintf(g_logbuf, "Error: Invalid --ppc-gap parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	dxx *= 1000;
	if (dxx < 0) {
	  ppc_gap = 0;
	} else if (dxx > 2147483646) {
	  ppc_gap = 0x7ffffffe;
	} else {
	  ppc_gap = (int32_t)(dxx * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "erm", 4)) {
	if (model_modifier & MODEL_MPERM) {
          if (calculation_type & CALC_MODEL) {
	    sprintf(g_logbuf, "Error: --perm cannot be used with --%s mperm.\n", (model_modifier & MODEL_ASSOC)? "assoc" : "model");
	    goto main_ret_INVALID_CMDLINE_2A;
	  } else {
	    logerrprint("Error: --perm cannot be used with --mperm.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	} else if ((calculation_type & CALC_GLM) && (glm_modifier & (GLM_MPERM | GLM_NO_SNP))) {
	  sprintf(g_logbuf, "Error: --perm cannot be used with --%s %s.\n", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", (glm_modifier & GLM_MPERM)? "mperm" : "no-snp");
	  goto main_ret_INVALID_CMDLINE_2A;
	} else if (family_info.dfam_modifier & DFAM_MPERM) {
	  logerrprint("Error: --perm cannot be used with --dfam mperm.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (calculation_type & CALC_CMH) {
          if (cluster.modifier & CLUSTER_CMH_MPERM) {
	    sprintf(g_logbuf, "Error: --perm cannot be used with --%s mperm.\n", (cluster.modifier & CLUSTER_CMH_BD)? "bd" : "mh");
	    goto main_ret_INVALID_CMDLINE_2A;
	  } else if (cluster.modifier & CLUSTER_CMH_PERM_BD) {
	    logerrprint("Error: --perm cannot be used with --bd perm-bd.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	}
	model_modifier |= MODEL_PERM;
        glm_modifier |= GLM_PERM;
        testmiss_modifier |= TESTMISS_PERM;
	family_info.tdt_modifier |= TDT_PERM;
	family_info.dfam_modifier |= DFAM_PERM;
	family_info.qfam_modifier |= QFAM_PERM;
	cluster.modifier |= CLUSTER_CMH_PERM;
	logprint("Note: --perm flag deprecated.  Use e.g. \"--model perm\".\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "erm-count", 10)) {
	model_modifier |= MODEL_PERM_COUNT;
	glm_modifier |= GLM_PERM_COUNT;
        testmiss_modifier |= TESTMISS_PERM_COUNT;
	family_info.tdt_modifier |= TDT_PERM_COUNT;
	family_info.dfam_modifier |= DFAM_PERM_COUNT;
        family_info.qfam_modifier |= QFAM_PERM_COUNT;
	cluster.modifier |= CLUSTER_CMH_PERM_COUNT;
	logprint("Note: --perm-count flag deprecated.  Use e.g. \"--model perm-count\".\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "2", 2)) {
	logerrprint("Error: --p2 has been provisionally retired.  Contact the developers if you need\nthis functionality.\n");
        goto main_ret_INVALID_CMDLINE_A;
      } else if (!memcmp(argptr2, "filter", 7)) {
	if (load_rare == LOAD_RARE_DOSAGE) {
	  logerrprint("Error: --pfilter is currently ignored by --dosage.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx)) {
	  sprintf(g_logbuf, "Error: Invalid --pfilter parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((dxx <= 0.0) || (dxx > 1.0)) {
	  logerrprint("Error: --pfilter threshold must be in (0, 1].\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	pfilter = dxx;
      } else if (!memcmp(argptr2, "erm-batch-size", 15)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &perm_batch_size)) {
	  sprintf(g_logbuf, "Error: Invalid --perm-batch-size parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "pc", 3)) {
	if (!(calculation_type & (CALC_NEIGHBOR | CALC_CLUSTER))) {
          logerrprint("Error: --ppc must be used with --cluster or --neigbour.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx)) {
	  sprintf(g_logbuf, "Error: Invalid --ppc parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((dxx <= 0.0) || (dxx >= 1.0)) {
	  logerrprint("Error: --ppc threshold must be between 0 and 1 exclusive.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        cluster.ppc = dxx;
      } else if (!memcmp(argptr2, "ool-size", 9)) {
	if (!(homozyg.modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE))) {
          logerrprint("Error: --pool-size must be used with --homozyg group{-verbose}.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &homozyg.pool_size_min) || (homozyg.pool_size_min == 1)) {
	  sprintf(g_logbuf, "Error: Invalid --pool-size parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "arameters", 10)) {
	if (!(calculation_type & CALC_GLM)) {
	  // drop --dosage since --covar-number has same functionality
	  logerrprint("Error: --parameters must be used with --linear or --logistic.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	retval = parse_name_ranges(param_ct, '-', &(argv[cur_arg]), &parameters_range_list, 1);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ca", 3)) {
#ifdef NOLAPACK
        logerrprint("Error: --pca requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	goto main_ret_INVALID_CMDLINE;
#else
	if (rel_info.modifier & REL_CALC_COV) {
	  logerrprint("Error: --pca flag cannot coexist with a covariance matrix calculation.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (parallel_tot > 1) {
	  logerrprint("Error: --parallel and --pca cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	ujj = 0;
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "header")) {
	    rel_info.modifier |= REL_PCA_HEADER;
	  } else if (!strcmp(argv[cur_arg + uii], "tabs")) {
            rel_info.modifier |= REL_PCA_TABS;
	  } else if (!strcmp(argv[cur_arg + uii], "var-wts")) {
            rel_info.modifier |= REL_PCA_VAR_WTS;
	  } else {
	    if (ujj || scan_posint_defcapx(argv[cur_arg + uii], &rel_info.pc_ct)) {
	      logerrprint("Error: Invalid --pca parameter sequence.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (rel_info.pc_ct > 8000) {
	      logerrprint("Error: --pca does not support more than 8000 PCs.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    ujj = 1;
	  }
	}
        calculation_type |= CALC_PCA;
#endif
      } else if (!memcmp(argptr2, "ca-cluster-names", 17)) {
	if (!(calculation_type & CALC_PCA)) {
	  logerrprint("Error: --pca-cluster-names must be used with --pca.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_and_flatten(&(rel_info.pca_cluster_names_flattened), &(argv[cur_arg + 1]), param_ct);
        if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ca-clusters", 12)) {
	if (!(calculation_type & CALC_PCA)) {
	  logerrprint("Error: --pca-clusters must be used with --pca.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&(rel_info.pca_clusters_fname), argv[cur_arg + 1], argptr, 0);
        if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "oo", 3)) {
	logprint("Note: --poo flag deprecated.  Use \"--tdt poo\".\n");
	family_info.tdt_modifier |= TDT_POO;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "arentdt1", 9)) {
	logprint("Note: --parentdt1 flag deprecated.  Use e.g.\n\"--tdt exact mperm=<value> parentdt1\".\n");
	family_info.tdt_modifier |= TDT_PARENPERM1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "arentdt2", 9)) {
	logprint("Note: --parentdt2 flag deprecated.  Use e.g.\n\"--tdt exact mperm=<value> parentdt2\".\n");
	family_info.tdt_modifier |= TDT_PARENPERM2;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "at", 3)) {
	logprint("Note: --pat flag deprecated.  Use e.g. \"--tdt poo mperm=<value> pat\".\n");
	family_info.tdt_modifier |= TDT_POOPERM_PAT;
	goto main_param_zero;
      } else if ((!memcmp(argptr2, "roxy-assoc", 11)) ||
                 (!memcmp(argptr2, "roxy-drop", 10)) ||
                 (!memcmp(argptr2, "roxy-impute", 12)) ||
                 (!memcmp(argptr2, "roxy-impute-threshold", 22)) ||
                 (!memcmp(argptr2, "roxy-genotypic-concordance", 27)) ||
                 (!memcmp(argptr2, "roxy-show-proxies", 18)) ||
                 (!memcmp(argptr2, "roxy-dosage", 12)) ||
                 (!memcmp(argptr2, "roxy-replace", 13)) ||
                 (!memcmp(argptr2, "roxy-verbose", 13))) {
	// Apparently there are no good alternatives for Y and MT imputation?
	// May want to modify this error message to suggest PLINK 1.07 for that
	// case.  (Since BEAGLE 4 is open source, it may be practical to build
	// a PLINK 2.1 which includes a port of its imputation routine, and
	// that can be written to handle Y/MT in a sane manner.  But that's not
	// happening before 2016.)
        logerrprint("Error: PLINK 1 proxy association and imputation commands have been retired due\nto poor accuracy.  (See Nothnagel M et al. (2009) A comprehensive evaluation of\nSNP genotype imputation.)  We suggest using another tool, such as BEAGLE 4 or\nIMPUTE2, for imputation instead, and performing association analysis on those\nresults.  (\"--recode vcf\" and --vcf can be used to exchange data with BEAGLE 4,\nwhile \"--recode oxford\" and --data let you work with IMPUTE2.)\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "file", 5)) {
	sprintf(g_logbuf, "Error: Unrecognized flag (%s).  (This is PLINK 1.9, not 2.x.)\n", argv[cur_arg]);
	goto main_ret_INVALID_CMDLINE_2;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'q':
      if (!memcmp(argptr2, "t-means", 8)) {
	if ((!(calculation_type & CALC_MODEL)) || (!(model_modifier & MODEL_ASSOC))) {
	  logerrprint("Error: --qt-means must be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (model_modifier & MODEL_DMASK) {
	  logerrprint("Error: --qt-means does not make sense with a case/control-specific --assoc\nmodifier.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --qt-means flag deprecated.  Use \"--assoc qt-means ...\".\n");
	model_modifier |= MODEL_QT_MEANS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "q-plot", 7)) {
        if (!mtest_adjust) {
	  logerrprint("Error: --qq-plot must be used with --adjust.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --qq-plot flag deprecated.  Use \"--adjust qq-plot\".\n");
	mtest_adjust |= ADJUST_QQ;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "match", 6)) {
        if (!(calculation_type & CALC_CLUSTER)) {
          logerrprint("Error: --qmatch must be used with --cluster.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&cluster.qmatch_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        if (param_ct == 2) {
	  if (alloc_string(&cluster.qmatch_missing_str, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	}
      } else if (!memcmp(argptr2, "t", 2)) {
        if (!cluster.qmatch_fname) {
          logerrprint("Error: --qt must be used with --qmatch.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&cluster.qt_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "-score-file", 12)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&score_info.data_fname, argv[cur_arg + 1], argptr, 0);
        if (retval) {
	  goto main_ret_1;
	}
        logprint("Note: --q-score-file flag deprecated.  Pass multiple parameters to\n--q-score-range instead.\n");
      } else if (!memcmp(argptr2, "-score-range", 13)) {
	if (score_info.data_fname) {
	  if (param_ct != 1) {
	    logerrprint("Error: --q-score-range must be given exactly one parameter when --q-score-file\nis present.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  retval = alloc_fname(&score_info.range_fname, argv[cur_arg + 1], argptr, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	  // no need for separate deprecation message
	} else {
	  if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 5)) {
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  retval = alloc_fname(&score_info.range_fname, argv[cur_arg + 1], argptr, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
          retval = alloc_fname(&score_info.data_fname, argv[cur_arg + 2], argptr, 0);
	  ujj = 0; // number of numeric parameters
	  for (uii = 3; uii <= param_ct; uii++) {
	    if (!strcmp(argv[cur_arg + uii], "header")) {
	      score_info.modifier |= SCORE_DATA_HEADER;
	    } else if (ujj == 2) {
              logerrprint("Error: --q-score-range takes at most two numeric parameters.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    } else {
	      if (scan_posint_cappedx(argv[cur_arg + uii], MAXLINEBUFLEN / 2, (uint32_t*)&ii)) {
                sprintf(g_logbuf, "Error: Invalid --q-score-range parameter '%s'.\n", argv[cur_arg + uii]);
                goto main_ret_INVALID_CMDLINE_WWA;
	      }
              if (!ujj) {
                score_info.data_varid_col = ii;
                score_info.data_col = ii + 1;
	      } else {
		if ((uint32_t)ii == score_info.data_varid_col) {
	          logerrprint("Error: --q-score-range variant ID and data column numbers cannot match.\n");
	          goto main_ret_INVALID_CMDLINE_A;
		}
		score_info.data_col = ii;
	      }
	      ujj++;
	    }
	  }
	}
      } else if (!memcmp(argptr2, "ual-scores", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_2col(&qual_filter, &(argv[cur_arg + 1]), argptr, param_ct);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "ual-threshold", 14)) {
	if (!qual_filter) {
	  logerrprint("Error: --qual-threshold must be used with --qual-scores.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &qual_min_thresh)) {
	  sprintf(g_logbuf, "Error: Invalid --qual-threshold parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (qual_min_thresh > qual_max_thresh) {
	  logerrprint("Error: --qual-threshold value cannot be larger than --qual-max-threshold value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "ual-max-threshold", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &qual_max_thresh)) {
	  sprintf(g_logbuf, "Error: Invalid --qual-max-threshold parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if ((!memcmp(argptr2, "fam", 4)) || (!memcmp(argptr2, "fam-parents", 12)) || (!memcmp(argptr2, "fam-between", 12)) || (!memcmp(argptr2, "fam-total", 10))) {
	if (calculation_type & CALC_QFAM) {
	  logerrprint("Error: Only one QFAM test can be run at a time.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "emp-se")) {
	    family_info.qfam_modifier |= QFAM_EMP_SE;
	  } else if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (family_info.qfam_modifier & QFAM_MPERM) {
	      sprintf(g_logbuf, "Error: --%s 'mperm' and 'perm' cannot be used together.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    family_info.qfam_modifier |= QFAM_PERM;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (family_info.qfam_modifier & QFAM_PERM) {
	      sprintf(g_logbuf, "Error: --%s 'mperm' and 'perm' cannot be used together.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_2A;
	    } else if (family_info.qfam_modifier & QFAM_MPERM) {
	      sprintf(g_logbuf, "Error: Duplicate --%s 'mperm' modifier.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_2;
	    }
	    if (scan_posint_defcapx(&(argv[cur_arg + uii][6]), &(family_info.qfam_mperm_val))) {
	      sprintf(g_logbuf, "Error: Invalid --%s mperm parameter '%s'.\n", argptr, &(argv[cur_arg + uii][6]));
              goto main_ret_INVALID_CMDLINE_WWA;
	    }
            family_info.qfam_modifier |= QFAM_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
            family_info.qfam_modifier |= QFAM_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
            sprintf(g_logbuf, "Error: Improper --%s mperm syntax.  (Use \"--%s mperm=<value>\".)\n", argptr, argptr);
            goto main_ret_INVALID_CMDLINE_WWA;
	  } else {
            sprintf(g_logbuf, "Error: Invalid --%s parameter '%s'.\n", argptr, argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	if (!(family_info.qfam_modifier & (QFAM_PERM | QFAM_MPERM))) {
	  sprintf(g_logbuf, "Error: --%s requires permutation.\n", argptr);
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (!memcmp(argptr2, "fam", 4)) {
          family_info.qfam_modifier |= QFAM_WITHIN1;
	} else if (!memcmp(argptr2, "fam-parents", 12)) {
          family_info.qfam_modifier |= QFAM_WITHIN2;
	} else if (!memcmp(argptr2, "fam-between", 12)) {
          family_info.qfam_modifier |= QFAM_BETWEEN;
	} else {
          family_info.qfam_modifier |= QFAM_TOTAL;
	}
	calculation_type |= CALC_QFAM;
      } else if ((!memcmp(argptr2, "ual-geno-scores", 16)) ||
                 (!memcmp(argptr2, "ual-geno-threshold", 19)) ||
                 (!memcmp(argptr2, "ual-geno-max-threshold", 23))) {
	logerrprint("Error: --qual-geno-scores has been provisionally retired, since it cannot be\nextended in a VCF-friendly manner.  Contact the developers if you prefer to\ncontinue using its interface.\n");
	goto main_ret_INVALID_CMDLINE;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'r':
      if (!memcmp(argptr2, "emove", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&removename, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "emove-fam", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&removefamname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "emove-clusters", 15)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&(cluster.remove_fname), argv[cur_arg + 1], argptr, 0);
        if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "emove-cluster-names", 20)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_and_flatten(&(cluster.remove_flattened), &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
        filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "el-cutoff", 10)) {
	if (parallel_tot > 1) {
	  logerrprint("Error: --parallel cannot be used with --rel-cutoff.  (Use a combination of\n--make-rel, --keep/--remove, and a filtering script.)\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (rel_info.pca_cluster_names_flattened || rel_info.pca_clusters_fname) {
	  logerrprint("Error: --pca-cluster-names/--pca-clusters cannot be used with --rel-cutoff.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (covar_fname) {
	  // corner case: --rel-cutoff screws up covariate filtered indices
          logerrprint("Error: --covar cannot be used with --rel-cutoff.\n");
	  goto main_ret_INVALID_CMDLINE;
	}

	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_doublex(argv[cur_arg + 1], &rel_info.cutoff) || (rel_info.cutoff <= 0.0) || (rel_info.cutoff >= 1.0)) {
	    sprintf(g_logbuf, "Error: Invalid --rel-cutoff parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	calculation_type |= CALC_REL_CUTOFF;
      } else if (!memcmp(argptr2, "egress-distance", 16)) {
	if (parallel_tot > 1) {
	  logerrprint("Error: --parallel and --regress-distance cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_posintptrx(argv[cur_arg + 1], &regress_iters) || (regress_iters == 1)) {
	    sprintf(g_logbuf, "Error: Invalid --regress-distance jackknife iteration count '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (param_ct == 2) {
	    if (scan_posint_defcapx(argv[cur_arg + 2], &regress_d)) {
	      sprintf(g_logbuf, "Error: Invalid --regress-distance jackknife delete parameter '%s'.\n", argv[cur_arg + 2]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	}
	calculation_type |= CALC_REGRESS_DISTANCE;
      } else if (!memcmp(argptr2, "egress-rel", 11)) {
	if (parallel_tot > 1) {
	  logerrprint("Error: --parallel and --regress-rel flags cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (rel_info.pca_cluster_names_flattened || rel_info.pca_clusters_fname) {
	  logerrprint("Error: --pca-cluster-names/--pca-clusters cannot be used with --regress-rel.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (scan_posintptrx(argv[cur_arg + 1], &rel_info.regress_rel_iters) || (rel_info.regress_rel_iters == 1)) {
	    sprintf(g_logbuf, "Error: Invalid --regress-rel jackknife iteration count '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (param_ct == 2) {
	    if (scan_posint_defcapx(argv[cur_arg + 2], &rel_info.regress_rel_d)) {
	      sprintf(g_logbuf, "Error: Invalid --regress-rel jackknife delete parameter '%s'.\n", argv[cur_arg + 2]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	  }
	}
	calculation_type |= CALC_REGRESS_REL;
      } else if ((!memcmp(argptr2, "egress-pcs", 11)) || (!memcmp(argptr2, "egress-pcs-distance", 20))) {
	logerrprint("Error: --regress-pcs has been retired.  Contact the developers if you need a\nlibrary which efficiently computes genotype residuals on the fly.\n");
        goto main_ret_INVALID_CMDLINE_A;
      } else if (!memcmp(argptr2, "ead-freq", 9)) {
	if (calculation_type & CALC_FREQ) {
	  logerrprint("Error: --freq and --read-freq flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (misc_flags & MISC_HET_SMALL_SAMPLE) {
	  logerrprint("Error: \"--het small-sample\" cannot currently be used with --read-freq.  Contact\nthe developers if you need to combine them.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&freqname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_ALL_REQ;
      } else if ((!memcmp(argptr2, "ecode", 6)) || (!memcmp(argptr2, "ecode ", 6))) {
	if (argptr2[5] == ' ') {
	  kk = 1;
	} else {
	  kk = 0;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4 - kk)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if ((!strcmp(argv[cur_arg + uii], "01")) || (!strcmp(argv[cur_arg + uii], "12"))) {
	    if (recode_modifier & (RECODE_A | RECODE_AD)) {
	      sprintf(g_logbuf, "Error: The '%s' modifier does not apply to --recode's A and AD output formats.\n", argv[cur_arg + uii]);
	      goto main_ret_INVALID_CMDLINE_2A;
	    } else if (recode_modifier & RECODE_VCF) {
	    main_recode_012_vcf_conflict:
	      logerrprint("Error: '01'/'12' cannot be used with --recode's VCF output formats.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    if (argv[cur_arg + uii][0] == '0') {
	      if (recode_modifier & RECODE_12) {
		logerrprint("Error: --recode '01' and '12' cannot be used together.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      recode_modifier |= RECODE_01;
	    } else {
	      if (recode_modifier & RECODE_01) {
		logerrprint("Error: --recode '01' and '12' cannot be used together.\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      recode_modifier |= RECODE_12;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "23")) {
	    if (recode_type_set(&recode_modifier, RECODE_23)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
          } else if ((!argv[cur_arg + uii][1]) && (tolower(argv[cur_arg + uii][0]) == 'a')) {
	    if (recode_type_set(&recode_modifier, RECODE_A)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
          } else if ((tolower(argv[cur_arg + uii][0]) == 'a') && (!strcmp(&(argv[cur_arg + uii][1]), "-transpose"))) {
	    if (recode_type_set(&recode_modifier, RECODE_A_TRANSPOSE)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if ((!argv[cur_arg + uii][2]) && match_upper(argv[cur_arg + uii], "AD")) {
	    if (recode_type_set(&recode_modifier, RECODE_AD)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "compound-genotypes")) {
	    if (recode_type_set(&recode_modifier, RECODE_COMPOUND)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (match_upper_counted(argv[cur_arg + uii], "HV", 2)) {
	    if (!argv[cur_arg + uii][2]) {
	      if (recode_type_set(&recode_modifier, RECODE_HV)) {
	        goto main_ret_INVALID_CMDLINE_A;
	      }
	    } else if (!strcmp(&(argv[cur_arg + uii][2]), "-1chr")) {
	      if (recode_type_set(&recode_modifier, RECODE_HV_1CHR)) {
	        goto main_ret_INVALID_CMDLINE_A;
	      }
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "tab")) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      logerrprint("Error: Multiple --recode delimiter modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    recode_modifier |= RECODE_TAB;
	  } else if (!strcmp(argv[cur_arg + uii], "tabx")) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      logerrprint("Error: Multiple --recode delimiter modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    recode_modifier |= RECODE_TAB | RECODE_DELIMX;
	  } else if (!strcmp(argv[cur_arg + uii], "spacex")) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      logerrprint("Error: Multiple --recode delimiter modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    recode_modifier |= RECODE_DELIMX;
	  } else if (!strcmp(argv[cur_arg + uii], "bgz")) {
	    recode_modifier |= RECODE_BGZ;
	  } else if (!strcmp(argv[cur_arg + uii], "gen-gz")) {
	    recode_modifier |= RECODE_GEN_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "beagle")) {
	    if (recode_type_set(&recode_modifier, RECODE_BEAGLE)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "beagle-nomap")) {
	    if (recode_type_set(&recode_modifier, RECODE_BEAGLE_NOMAP)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "bimbam")) {
	    if (recode_type_set(&recode_modifier, RECODE_BIMBAM)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "bimbam-1chr")) {
	    if (recode_type_set(&recode_modifier, RECODE_BIMBAM_1CHR)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "fastphase")) {
	    if (recode_type_set(&recode_modifier, RECODE_FASTPHASE)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "fastphase-1chr")) {
	    if (recode_type_set(&recode_modifier, RECODE_FASTPHASE_1CHR)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "lgen")) {
	    if (recode_type_set(&recode_modifier, RECODE_LGEN)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "lgen-ref")) {
	    if (recode_type_set(&recode_modifier, RECODE_LGEN_REF)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "list")) {
	    if (recode_type_set(&recode_modifier, RECODE_LIST)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "oxford")) {
            if (recode_type_set(&recode_modifier, RECODE_OXFORD)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "ped")) {
	    if (recode_type_set(&recode_modifier, RECODE_PED)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "rlist")) {
	    if (recode_type_set(&recode_modifier, RECODE_RLIST)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "structure")) {
	    if (recode_type_set(&recode_modifier, RECODE_STRUCTURE)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "transpose")) {
	    if (recode_type_set(&recode_modifier, RECODE_TRANSPOSE)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "vcf")) {
	    if (recode_modifier & RECODE_VCF) {
	    main_recode_vcf_conflict:
	      logerrprint("Error: Conflicting or redundant --recode modifiers.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (recode_modifier & (RECODE_01 | RECODE_12)) {
	      goto main_recode_012_vcf_conflict;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_VCF)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "vcf-fid")) {
	    if (recode_modifier & (RECODE_VCF | RECODE_IID)) {
	      goto main_recode_vcf_conflict;
	    } else if (recode_modifier & (RECODE_01 | RECODE_12)) {
	      goto main_recode_012_vcf_conflict;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_VCF)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    recode_modifier |= RECODE_FID;
	  } else if (!strcmp(argv[cur_arg + uii], "vcf-iid")) {
	    if (recode_modifier & (RECODE_VCF | RECODE_FID)) {
	      goto main_recode_vcf_conflict;
	    } else if (recode_modifier & (RECODE_01 | RECODE_12)) {
	      goto main_recode_012_vcf_conflict;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_VCF)) {
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    recode_modifier |= RECODE_IID;
	  } else if (!strcmp(argv[cur_arg + uii], "include-alt")) {
	    recode_modifier |= RECODE_INCLUDE_ALT;
	  } else if (!strcmp(argv[cur_arg + uii], "omit-nonmale-y")) {
	    recode_modifier |= RECODE_OMIT_NONMALE_Y;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --recode parameter '%s'.%s\n", argv[cur_arg + uii], ((uii == param_ct) && (!outname_end))? " (Did you forget '--out'?)" : "");
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	if (!(recode_modifier & RECODE_TYPEMASK)) {
	  // thought about including a warning here, and then outright
	  // prohibiting parameterless --recode in v2.0, but that usage is
	  // simply too common.  new plan: --export requires an explicit output
	  // format, but --recode will continue to be translated without a
	  // warning even when defaulting to ped output.
	  if (recode_require_format) {
	    logerrprint("Error: --export requires an output format.  (Did you forget 'ped' or 'vcf'?)\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  recode_modifier |= RECODE_PED;
	}
	if ((recode_modifier & RECODE_INCLUDE_ALT) && (!(recode_modifier & (RECODE_A | RECODE_AD)))) {
	  logerrprint("Error: The 'include-alt' modifier only applies to --recode's A and AD output\nformats.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if ((recode_modifier & RECODE_OMIT_NONMALE_Y) && (!(recode_modifier & (RECODE_LIST | RECODE_RLIST)))) {
	  logerrprint("Error: The 'omit-nonmale-y' modifier only applies to --recode's list and rlist\noutput formats.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if ((recode_modifier & RECODE_BGZ) && (!(recode_modifier & RECODE_VCF))) {
	  logerrprint("Error: The 'bgz' modifier only applies to --recode's VCF output formats.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if ((recode_modifier & RECODE_GEN_GZ) && (!(recode_modifier & RECODE_OXFORD))) {
	  logerrprint("Error: The 'gen-gz' modifier only applies to --recode's oxford output format.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	calculation_type |= CALC_RECODE;
      } else if (!memcmp(argptr2, "ecode-whap", 11)) {
        logerrprint("Error: --recode-whap flag retired since WHAP is no longer supported.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ecode-allele", 13)) {
	if (!(recode_modifier & (RECODE_A | RECODE_A_TRANSPOSE | RECODE_AD))) {
	  logerrprint("Error: --recode-allele must be used with --recode A/A-transpose/AD.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (recode_modifier & (RECODE_01 | RECODE_12)) {
	  logerrprint("Error: --recode-allele cannot be used with --recode 01/12.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&recode_allele_name, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eference", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&lgen_reference_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	lgen_modifier |= LGEN_REFERENCE;
      } else if (!memcmp(argptr2, "ead-genome", 11)) {
	if (calculation_type & CALC_GENOME) {
          logerrprint("Error: --read-genome cannot be used with --genome.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (!(calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR))) {
          logerrprint("Error: --read-genome cannot be used without --cluster or --neighbour.\n");
          goto main_ret_INVALID_CMDLINE_A;
	} else if ((cluster.ppc == 0.0) && ((calculation_type & (CALC_DISTANCE | CALC_PLINK1_DISTANCE_MATRIX)) || read_dists_fname)) {
          logerrprint("Error: --read-genome is pointless with --distance, --distance-matrix, and\n--read-dists unless --ppc is also present.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        retval = alloc_fname(&read_genome_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ead-genome-list", 19)) {
	logerrprint("Error: --read-genome-list flag retired.  Use --parallel + Unix cat instead.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ead-genome-minimal", 19)) {
	logerrprint("Error: --read-genome-minimal flag retired.  Use --genome gz + --read-genome\ninstead.");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ead-dists", 10)) {
	if (calculation_type & (CALC_DISTANCE | CALC_PLINK1_DISTANCE_MATRIX)) {
	  logerrprint("Error: --read-dists cannot be used with a distance matrix calculation.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (cluster.modifier & CLUSTER_MISSING) {
          logerrprint("Error: --read-dists cannot be used with \"--cluster missing\".\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&read_dists_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        if (param_ct == 2) {
          retval = alloc_fname(&read_dists_id_fname, argv[cur_arg + 2], argptr, 0);
          if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "el-check", 9)) {
        if (!(calculation_type & CALC_GENOME)) {
          logerrprint("Error: --rel-check must be used with --genome.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        logprint("Note: --rel-check flag deprecated.  Use \"--genome rel-check\".\n");
        genome_modifier |= GENOME_REL_CHECK;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ecessive", 9)) {
	if (!(calculation_type & CALC_GLM)) {
	  logerrprint("Error: --recessive must be used with --linear or --logistic.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_DOMINANT)) {
	  sprintf(g_logbuf, "Error: --recessive conflicts with a --%s modifier.\n", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	logprint("Note: --recessive flag deprecated.  Use e.g. \"--linear recessive\" (and\n\"--condition-list <filename> recessive\" to change covariate coding).\n");
	glm_modifier |= GLM_RECESSIVE | GLM_CONDITION_RECESSIVE;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else if ((*argptr2 == '\0') || (!memcmp(argptr2, "2", 2))) {
	if (calculation_type & CALC_LD) {
          logerrprint("Error: --r and --r2 cannot be used together.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (calculation_type & CALC_BLOCKS) {
	  // prevent --ld-window-r2 conflict
	  logerrprint("Error: --r/--r2 cannot be used with --blocks.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 6)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (*argptr2 == '2') {
          ld_info.modifier |= LD_R2;
	}
	if (matrix_flag_state) {
	  matrix_flag_state = 2;
	  ld_info.modifier |= LD_MATRIX_SQ | LD_MATRIX_SPACES;
	}
	for (uii = 1; uii <= param_ct; uii++) {
          if (!strcmp(argv[cur_arg + uii], "square")) {
	    if (ld_info.modifier & LD_MATRIX_SHAPEMASK) {
	      logerrprint("Error: Multiple --r/--r2 shape modifiers.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (ld_info.modifier & (LD_INTER_CHR | LD_INPHASE | LD_DX | LD_WITH_FREQS)) {
	    main_r2_matrix_conflict:
              sprintf(g_logbuf, "Error: --r/--r2 '%s' cannot be used with matrix output.\n", (ld_info.modifier & LD_INTER_CHR)? "inter-chr" : ((ld_info.modifier & LD_INPHASE)? "in-phase" : ((ld_info.modifier & LD_DX)? "d'/'dprime'/'dprime-signed" : "with-freqs")));
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    ld_info.modifier |= LD_MATRIX_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if (ld_info.modifier & LD_MATRIX_SHAPEMASK) {
	      logerrprint("Error: Multiple --r/--r2 shape modifiers.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (ld_info.modifier & (LD_INTER_CHR | LD_INPHASE | LD_DX | LD_WITH_FREQS)) {
	      goto main_r2_matrix_conflict;
	    }
	    ld_info.modifier |= LD_MATRIX_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if (ld_info.modifier & LD_MATRIX_SHAPEMASK) {
	      logerrprint("Error: Multiple --r/--r2 shape modifiers.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (ld_info.modifier & (LD_INTER_CHR | LD_INPHASE | LD_DX | LD_WITH_FREQS)) {
	      goto main_r2_matrix_conflict;
	    }
	    ld_info.modifier |= LD_MATRIX_TRI;
	  } else if (!strcmp(argv[cur_arg + uii], "inter-chr")) {
            ld_info.modifier |= LD_INTER_CHR;
            if (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_MATRIX_BIN | LD_MATRIX_BIN4 | LD_MATRIX_SPACES)) {
	      goto main_r2_matrix_conflict;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (ld_info.modifier & (LD_MATRIX_BIN | LD_MATRIX_BIN4)) {
	      logerrprint("Error: Conflicting --r/--r2 modifiers.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    ld_info.modifier |= LD_REPORT_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (ld_info.modifier & (LD_INTER_CHR | LD_INPHASE | LD_DX | LD_WITH_FREQS)) {
	      goto main_r2_matrix_conflict;
	    } else if (ld_info.modifier & (LD_REPORT_GZ | LD_MATRIX_BIN4)) {
	      logerrprint("Error: Conflicting --r/--r2 modifiers.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (ld_info.modifier & LD_MATRIX_SPACES) {
	      logerrprint("Error: --r/--r2 'bin[4]' and 'spaces' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    ld_info.modifier |= LD_MATRIX_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "bin4")) {
	    if (ld_info.modifier & (LD_INTER_CHR | LD_INPHASE | LD_DX | LD_WITH_FREQS)) {
	      goto main_r2_matrix_conflict;
	    } else if (ld_info.modifier & (LD_REPORT_GZ | LD_MATRIX_BIN)) {
	      logerrprint("Error: Conflicting --r/--r2 modifiers.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (ld_info.modifier & LD_MATRIX_SPACES) {
	      logerrprint("Error: --r/--r2 'bin[4]' and 'spaces' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    ld_info.modifier |= LD_MATRIX_BIN4;
	  } else if (!strcmp(argv[cur_arg + uii], "single-prec")) {
	    logerrprint("Error: --r/--r2 'single-prec' modifier has been retired.  Use 'bin4'.\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else if (!strcmp(argv[cur_arg + uii], "spaces")) {
	    if (ld_info.modifier & (LD_INTER_CHR | LD_INPHASE | LD_DX | LD_WITH_FREQS)) {
	      goto main_r2_matrix_conflict;
	    } else if (ld_info.modifier & (LD_MATRIX_BIN | LD_MATRIX_BIN4)) {
	      logerrprint("Error: --r/--r2 'bin[4]' and 'spaces' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    ld_info.modifier |= LD_MATRIX_SPACES;
	  } else if (!strcmp(argv[cur_arg + uii], "in-phase")) {
	    ld_info.modifier |= LD_INPHASE;
            if (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_MATRIX_BIN | LD_MATRIX_BIN4 | LD_MATRIX_SPACES)) {
	      goto main_r2_matrix_conflict;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "d")) {
	    ld_info.modifier |= LD_D;
            if (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_MATRIX_BIN | LD_MATRIX_BIN4 | LD_MATRIX_SPACES)) {
	      goto main_r2_matrix_conflict;
	    }
	    if (ld_info.modifier & (LD_DPRIME | LD_DPRIME_SIGNED)) {
	      logerrprint("Error: --r/--r2 'd' cannot be used with 'dprime'/'dprime-signed'.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "dprime")) {
	    ld_info.modifier |= LD_DPRIME;
            if (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_MATRIX_BIN | LD_MATRIX_BIN4 | LD_MATRIX_SPACES)) {
	      goto main_r2_matrix_conflict;
	    }
	    if (ld_info.modifier & (LD_D | LD_DPRIME_SIGNED)) {
	      logerrprint("Error: --r/--r2 'dprime' cannot be used with 'd'/'dprime-signed'.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "dprime-signed")) {
	    ld_info.modifier |= LD_DPRIME_SIGNED;
            if (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_MATRIX_BIN | LD_MATRIX_BIN4 | LD_MATRIX_SPACES)) {
	      goto main_r2_matrix_conflict;
	    }
	    if (ld_info.modifier & (LD_D | LD_DPRIME)) {
	      logerrprint("Error: --r/--r2 'dprime-signed' cannot be used with 'd'/'dprime'.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "with-freqs")) {
	    ld_info.modifier |= LD_WITH_FREQS;
            if (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_MATRIX_BIN | LD_MATRIX_BIN4 | LD_MATRIX_SPACES)) {
	      goto main_r2_matrix_conflict;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "yes-really")) {
	    ld_info.modifier |= LD_YES_REALLY;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --r/--r2 parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        if ((ld_info.modifier & (LD_MATRIX_BIN | LD_MATRIX_BIN4)) && (!(ld_info.modifier & LD_MATRIX_SHAPEMASK))) {
          ld_info.modifier |= LD_MATRIX_SQ;
	}
	if ((ld_info.modifier & LD_MATRIX_SPACES) && (!(ld_info.modifier & LD_MATRIX_SHAPEMASK))) {
	  logerrprint("Error: --r/--r2 'spaces' modifier must be used with a shape modifier.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	if ((ld_info.modifier & LD_WEIGHTED_X) && (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_INTER_CHR))) {
	  logerrprint("Error: --ld-xchr 3 cannot be used with --r/--r2 non-windowed reports.\n");
          goto main_ret_INVALID_CMDLINE_A;
	} else if (ld_info.modifier & LD_WEIGHTED_X) {
	  logerrprint("Error: --r/--r2 + --ld-xchr 3 is not currently supported.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	calculation_type |= CALC_LD;
      } else if (!memcmp(argptr2, "eal-ref-alleles", 16)) {
	if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --real-ref-alleles has no effect with %s.\n", (load_rare == LOAD_RARE_CNV)? "a .cnv fileset" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        misc_flags |= MISC_REAL_REF_ALLELES | MISC_KEEP_ALLELE_ORDER;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "ange", 5)) {
        if (extractname) {
	  misc_flags |= MISC_EXTRACT_RANGE;
	} else if (!excludename) {
	  logerrprint("Error: --range must be used with --extract and/or --exclude.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (excludename) {
	  misc_flags |= MISC_EXCLUDE_RANGE;
	}
	logprint("Note: --range flag deprecated.  Use e.g. \"--extract range <filename>\".\n");
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 's':
      if (!memcmp(argptr2, "eed", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	rseed_ct = param_ct;
	rseeds = (uint32_t*)malloc(param_ct * sizeof(int32_t));
	for (uii = 1; uii <= param_ct; uii++) {
	  if (scan_uint_cappedx(argv[cur_arg + uii], 0xffffffffU, &(rseeds[uii - 1]))) {
	    sprintf(g_logbuf, "Error: Invalid --seed parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
      } else if (!memcmp(argptr2, "ample", 6)) {
	if ((load_params & (LOAD_PARAMS_TEXT_ALL | LOAD_PARAMS_BFILE_ALL)) || load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	} else if (!(load_params & LOAD_PARAMS_OX_ALL)) {
	  logerrprint("Error: --sample cannot be used without --data, --gen, or --bgen.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	load_params |= LOAD_PARAMS_OXSAMPLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logerrprint("Error: --sample parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "np", 3)) {
        if (markername_from) {
	  logerrprint("Error: --snp cannot be used with --from.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (marker_pos_start != -1) {
	  logerrprint("Error: --snp cannot be used with --from-bp/-kb/-mb.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if ((!are_all_words_zero(chrom_info.chrom_mask, CHROM_MASK_INITIAL_WORDS)) || chrom_info.incl_excl_name_stack) {
	  logerrprint("Error: --snp cannot be used with --autosome[-xy] or --[not-]chr.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (markername_snp) {
          logerrprint("Error: --snp cannot be used with --exclude-snp.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (snps_range_list.names) {
          logerrprint("Error: --snp cannot be used with --exclude-snps.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&markername_snp, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
        filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "nps", 4)) {
	if (markername_from) {
	  logerrprint("Error: --snps cannot be used with --from.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (marker_pos_start != -1) {
	  logerrprint("Error: --snps cannot be used with --from-bp/-kb/-mb.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (markername_snp) {
	  logerrprint("Error: --snps cannot be used with --snp or --exclude-snp.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (snps_range_list.names) {
	  logerrprint("Error: --snps cannot be used with --exclude-snps.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	// mise well allow --snps + --autosome/--autosome-xy/--chr/--not-chr
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &snps_range_list, 0);
	if (retval) {
	  goto main_ret_1;
	}
        filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "et-hh-missing", 14)) {
	misc_flags |= MISC_SET_HH_MISSING;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et-mixed-mt-missing", 20)) {
        misc_flags |= MISC_SET_MIXED_MT_MISSING;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "et", 3)) {
	if (set_info.fname) {
	  logerrprint("Error: --set cannot be used with --make-set.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --set cannot be used with %s.\n", (load_rare == LOAD_RARE_CNV)? "a .cnv fileset" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&set_info.fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_BIM_REQ;
      } else if (!memcmp(argptr2, "et-collapse-all", 16)) {
	if (!set_info.fname) {
	  logerrprint("Error: --set-collapse-all must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (set_info.modifier & SET_MAKE_COLLAPSE_GROUP) {
	  logerrprint("Error: --set-collapse-all cannot be used with --make-set-collapse-group or\n--make-set-complement-group.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (set_info.merged_set_name) {
	  logerrprint("Error: --set-collapse-all cannot be used with --make-set-complement-all.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&set_info.merged_set_name, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "et-names", 9)) {
	if (!set_info.fname) {
	  logerrprint("Error: --set-names must be used with --set/--make-set.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_and_flatten(&(set_info.setnames_flattened), &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ubset", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (!set_info.fname) {
	  logerrprint("Error: --subset must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&set_info.subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if ((!memcmp(argptr2, "imulate", 8)) || (!memcmp(argptr2, "imulate-qt", 11))) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&simulate_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (argptr2[7] == '-') {
	  simulate_flags |= SIMULATE_QT;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "tags")) {
	    if (simulate_flags & SIMULATE_HAPS) {
	      sprintf(g_logbuf, "Error: --%s 'tags' and 'haps' modifiers cannot be used together.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    simulate_flags |= SIMULATE_TAGS;
	  } else if (!strcmp(argv[cur_arg + uii], "haps")) {
	    if (simulate_flags & SIMULATE_TAGS) {
	      sprintf(g_logbuf, "Error: --%s 'tags' and 'haps' modifiers cannot be used together.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
	    simulate_flags |= SIMULATE_HAPS;
	  } else if (match_upper(argv[cur_arg + uii], "ACGT")) {
	    if (simulate_flags & (SIMULATE_1234 | SIMULATE_12)) {
	      sprintf(g_logbuf, "Error: --%s 'acgt' modifier cannot be used with '1234' or '12'.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
            simulate_flags |= SIMULATE_ACGT;
	  } else if (!strcmp(argv[cur_arg + uii], "1234")) {
	    if (simulate_flags & (SIMULATE_ACGT | SIMULATE_12)) {
	      sprintf(g_logbuf, "Error: --%s '1234' modifier cannot be used with 'acgt' or '12'.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
            simulate_flags |= SIMULATE_1234;
	  } else if (!strcmp(argv[cur_arg + uii], "12")) {
	    if (simulate_flags & (SIMULATE_ACGT | SIMULATE_1234)) {
	      sprintf(g_logbuf, "Error: --%s '12' modifier cannot be used with 'acgt' or '1234'.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_2A;
	    }
            simulate_flags |= SIMULATE_12;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --%s parameter '%s'.\n", argptr, argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	load_rare = LOAD_RARE_SIMULATE;
      } else if (!memcmp(argptr2, "imulate-ncases", 15)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (load_rare != LOAD_RARE_SIMULATE) {
	  logerrprint("Error: --simulate-ncases must be used with --simulate.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &simulate_cases)) {
	  sprintf(g_logbuf, "Error: Invalid --simulate-ncases parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "imulate-ncontrols", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (load_rare != LOAD_RARE_SIMULATE) {
	  sprintf(g_logbuf, "Error: --simulate-ncontrols must be used with --simulate.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &simulate_controls)) {
	  sprintf(g_logbuf, "Error: Invalid --simulate-ncontrols parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if ((!simulate_controls) && (!simulate_cases)) {
	  logerrprint("Error: \"--simulate-ncases 0\" cannot be used with \"--simulate-ncontrols 0\".\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "imulate-prevalence", 19)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &simulate_prevalence) || (simulate_prevalence < 0) || (simulate_prevalence > 1)) {
	  sprintf(g_logbuf, "Error: Invalid --simulate-prevalence parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "imulate-label", 14)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&simulate_label, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "imulate-missing", 16)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &simulate_missing) || (simulate_missing < 0) || (simulate_missing > 1)) {
	  sprintf(g_logbuf, "Error: Invalid --simulate-missing parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "imulate-n", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (load_rare == LOAD_RARE_SIMULATE) {
	  logerrprint("Error: --simulate-n must be used with --simulate-qt, not --simulate.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &simulate_qt_samples)) {
	  sprintf(g_logbuf, "Error: Invalid --simulate-n parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "imulate-haps", 13)) {
	if (simulate_flags & SIMULATE_TAGS) {
	  logerrprint("Error: --simulate-tags cannot be used with --simulate-haps.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --simulate-haps flag deprecated.  Use e.g. \"--simulate haps\".\n");
	simulate_flags |= SIMULATE_HAPS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "imulate-tags", 13)) {
	if (simulate_flags & SIMULATE_HAPS) {
	  logerrprint("Error: --simulate-tags cannot be used with --simulate-haps.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --simulate-tags flag deprecated.  Use e.g. \"--simulate tags\".\n");
	simulate_flags |= SIMULATE_TAGS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ex", 3)) {
	if (load_rare == LOAD_RARE_DOSAGE) {
	  dosage_info.modifier |= DOSAGE_SEX;
	} else {
	  if (!(calculation_type & CALC_GLM)) {
	    logerrprint("Error: --sex must be used with --linear/--logistic/--dosage.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  } else if (glm_modifier & GLM_NO_X_SEX) {
	    sprintf(g_logbuf, "Error: --sex conflicts with a --%s modifier.\n", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear");
	    goto main_ret_INVALID_CMDLINE_2A;
	  }
	  logprint("Note: --sex flag deprecated.  Use e.g. \"--linear sex\".\n");
	  glm_modifier |= GLM_SEX;
	}
	goto main_param_zero;
      } else if (!memcmp(argptr2, "tandard-beta", 13)) {
	if (((!(calculation_type & CALC_GLM)) || (glm_modifier & GLM_LOGISTIC)) && (!(dosage_info.modifier & DOSAGE_GLM))) {
	  logerrprint("Error: --standard-beta must be used with --linear or --dosage.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	// "--linear standard-beta" is ALSO deprecated now for plink2's
	// --{covar-}variance-normalize, so this note is pointless.
	// logprint("Note: --standard-beta flag deprecated.  Use e.g. \"--linear standard-beta\".\n");
	glm_modifier |= GLM_STANDARD_BETA;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et-table", 9)) {
	if (!set_info.fname) {
	  logerrprint("Error: --set-table must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        calculation_type |= CALC_WRITE_SET;
	set_info.modifier |= SET_WRITE_TABLE;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "et-test", 8)) {
	if (!set_info.fname) {
	  logerrprint("Error: --set-test must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --set-test flag deprecated.  Use e.g. \"--assoc perm set-test\".\n");
	if (calculation_type & CALC_MODEL) {
          model_modifier |= MODEL_SET_TEST;
	}
	if (calculation_type & CALC_GLM) {
	  if (glm_modifier & GLM_NO_SNP) {
	    logerrprint("Error: --set-test cannot be used with --linear/--logistic no-snp.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  glm_modifier |= GLM_SET_TEST;
	}
	if (calculation_type & CALC_TDT) {
	  family_info.tdt_modifier |= TDT_SET_TEST;
	}
	if (calculation_type & CALC_DFAM) {
	  family_info.dfam_modifier |= DFAM_SET_TEST;
	}
	if ((calculation_type & CALC_CMH) && (!(cluster.modifier & CLUSTER_CMH2))) {
	  cluster.modifier |= CLUSTER_CMH_SET_TEST;
	}
	if ((epi_info.modifier & (EPI_FAST | EPI_REG)) && (!(epi_info.modifier & EPI_SET_BY_ALL))) {
	  epi_info.modifier |= EPI_SET_BY_SET;
	}
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et-p", 5)) {
	if (!set_info.fname) {
	  logerrprint("Error: --set-p must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx <= 0) || (dxx > 1)) {
	  sprintf(g_logbuf, "Error: Invalid --set-p parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	set_info.set_p = dxx;
      } else if (!memcmp(argptr2, "et-r2", 5)) {
	if (!set_info.fname) {
	  logerrprint("Error: --set-r2 must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        uii = 1;
        if (!strcmp(argv[cur_arg + 1], "write")) {
          uii = 2;
          set_info.modifier |= SET_R2_WRITE;
	} else if (param_ct == 2) {
          if (!strcmp(argv[cur_arg + 2], "write")) {
            set_info.modifier |= SET_R2_WRITE;
	  } else {
            sprintf(g_logbuf, "Error: Invalid --set-r2 parameter '%s'.\n", argv[cur_arg + 2]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        if (uii <= param_ct) {
	  if (scan_doublex(argv[cur_arg + uii], &dxx) || (dxx < 0.0)) {
	    sprintf(g_logbuf, "Error: Invalid --set-r2 parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (dxx > 0.0) {
	    // greater than 1 = no LD check.  (it still happens with a
	    // parameter of 1.)
	    set_info.set_r2 = dxx;
	  } else {
	    set_info.set_max = 1;
	  }
	}
      } else if (!memcmp(argptr2, "et-r2-phase", 11)) {
        logerrprint("Error: --set-r2-phase is provisionally retired.  Contact the developers if you\nneed this function.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "et-max", 5)) {
	if (!set_info.fname) {
	  logerrprint("Error: --set-max must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &set_info.set_max)) {
	  sprintf(g_logbuf, "Error: Invalid --set-max parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "et-test-lambda", 15)) {
	if (!set_info.fname) {
	  logerrprint("Error: --set-test-lambda must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &set_info.set_test_lambda)) {
	  sprintf(g_logbuf, "Error: Invalid --set-test-lambda parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (set_info.set_test_lambda < 1) {
	  logprint("Note: --set-test-lambda parameter set to 1.\n");
	  set_info.set_test_lambda = 1;
	}
      } else if (!memcmp(argptr2, "et-by-all", 10)) {
	if (!set_info.fname) {
	  logerrprint("Error: --set-by-all must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (!(epi_info.modifier & (EPI_FAST | EPI_REG))) {
	  logerrprint("Error: --set-by-all must be used with --{fast-}epistasis.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --set-by-all flag deprecated.  Use e.g. \"--fast-epistasis set-by-all\".\n");
        epi_info.modifier |= EPI_SET_BY_ALL;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "nps-only", 9)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  const char* cur_modif = argv[cur_arg + 1];
	  const uint32_t slen = strlen(cur_modif);
	  if ((slen != 9) || memcmp(cur_modif, "just-acgt", 9)) {
	    if ((slen != 5) || (memcmp(cur_modif, "no-", 3)) || (!match_upper(&(cur_modif[3]), "DI"))) {
	      sprintf(g_logbuf, "Error: Invalid --snps-only parameter '%s'.\n", argv[cur_arg + 1]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    logerrprint("Warning: --snps-only 'no-DI' modifier has been changed to 'just-acgt', and now\nexcludes all non-ACGT allele codes.\n");
	  }
          misc_flags |= MISC_SNPS_ONLY_JUST_ACGT;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_SNPS_ONLY | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "plit-x", 7)) {
	if (misc_flags & MISC_MERGEX) {
          logerrprint("Error: --split-x cannot be used with --merge-x.\n");
          goto main_ret_INVALID_CMDLINE_A;
	} else if ((chrom_info.xymt_codes[X_OFFSET] == -2) || (chrom_info.xymt_codes[XY_OFFSET] == -2)) {
	  logerrprint("Error: --split-x must be used with a chromosome set containing X and XY codes.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	ujj = param_ct;
	// 'no-fail' must be either the first or last parameter.
	if (!strcmp(argv[cur_arg + 1], "no-fail")) {
	  misc_flags |= MISC_SPLIT_MERGE_NOFAIL;
	  uii = 2;
	} else if (!strcmp(argv[cur_arg + param_ct], "no-fail")) {
	  misc_flags |= MISC_SPLIT_MERGE_NOFAIL;
	  ujj--;
	}
	if ((ujj < uii) || (ujj == uii + 2)) {
	  logerrprint("Error: Invalid --split-x parameter sequence.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (ujj == uii) {
          if ((!strcmp(argv[cur_arg + uii], "b36")) || (!strcmp(argv[cur_arg + uii], "hg18"))) {
            splitx_bound1 = 2709521;
            splitx_bound2 = 154584237;
	  } else if ((!strcmp(argv[cur_arg + uii], "b37")) || (!strcmp(argv[cur_arg + uii], "hg19"))) {
            splitx_bound1 = 2699520;
            splitx_bound2 = 154931044;
	  } else if ((!strcmp(argv[cur_arg + uii], "b38")) || (!strcmp(argv[cur_arg + uii], "hg38"))) {
            splitx_bound1 = 2781479;
            splitx_bound2 = 155701383;
	  } else {
            sprintf(g_logbuf, "Error: Unrecognized --split-x build code '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (chrom_info.species != SPECIES_HUMAN) {
	    logerrprint("Error: --split-x build codes cannot be used with nonhuman chromosome sets.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	} else {
	  if (scan_uint_defcapx(argv[cur_arg + uii], &splitx_bound1)) {
	    sprintf(g_logbuf, "Error: Invalid --split-x parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  if (scan_posint_defcapx(argv[cur_arg + ujj], &splitx_bound2) || (splitx_bound2 <= splitx_bound1)) {
	    sprintf(g_logbuf, "Error: Invalid --split-x parameter '%s'.\n", argv[cur_arg + ujj]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
      } else if (!memcmp(argptr2, "et-missing-snp-ids", 19)) {
	logerrprint("Error: --set-missing-snp-ids has been retired.  Use --set-missing-var-ids with\n'$1' and '$2' instead.\n");
	goto main_ret_INVALID_CMDLINE_A;
      } else if (!memcmp(argptr2, "et-missing-nonsnp-ids", 22)) {
	logerrprint("Error: --set-missing-nonsnp-ids has been retired.  Use --set-missing-var-ids\nwith '$1' and '$2' instead.\n");
	goto main_ret_INVALID_CMDLINE_A;
      } else if (!memcmp(argptr2, "et-missing-var-ids", 19)) {
	if (load_rare & (LOAD_RARE_CNV | LOAD_RARE_DOSAGE)) {
	  sprintf(g_logbuf, "Error: --set-missing-var-ids cannot be used with %s.\n", (load_rare == LOAD_RARE_CNV)? "a .cnv fileset" : "--dosage");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	if (!valid_varid_template_string(argv[cur_arg + 1], "--set-missing-var-ids")) {
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (alloc_string(&missing_mid_template, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "core", 5)) {
	if (load_rare & LOAD_RARE_DOSAGE) {
	  if (dosage_info.modifier & DOSAGE_OCCUR) {
	    logerrprint("Error: --dosage 'occur' mode cannot be used with --score.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if (dosage_info.modifier & DOSAGE_ZOUT) {
	    logerrprint("Error: --dosage 'Zout' modifier cannot be used with --score.\n");
            goto main_ret_INVALID_CMDLINE_A;
	  }
	  if ((glm_modifier & GLM_STANDARD_BETA) || (dosage_info.modifier & DOSAGE_SEX)) {
	    logerrprint("Error: --dosage + --score cannot be used with association analysis\nmodifiers/flags.\n");
            goto main_ret_INVALID_CMDLINE_A;
	  }
	  dosage_info.modifier += (DOSAGE_SCORE - DOSAGE_GLM);
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 9)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&score_info.fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	ujj = 0; // number of numeric parameters so far
        for (uii = 2; uii <= param_ct; uii++) {
          if (!strcmp(argv[cur_arg + uii], "header")) {
            score_info.modifier |= SCORE_HEADER;
	  } else if (!strcmp(argv[cur_arg + uii], "sum")) {
	    if (score_info.modifier & SCORE_NO_MEAN_IMPUTATION) {
	      logerrprint("Error: --score 'sum' and 'no-mean-imputation' modifiers cannot be used\ntogether.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    } else if (dosage_info.modifier & DOSAGE_SCORE_NOSUM) {
	      logerrprint("Error: --score 'sum' and 'no-sum' modifiers conflict.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    score_info.modifier |= SCORE_SUM;
	  } else if (!strcmp(argv[cur_arg + uii], "no-mean-imputation")) {
	    if (score_info.modifier & SCORE_SUM) {
	      logerrprint("Error: --score 'sum' and 'no-mean-imputation' modifiers cannot be used\ntogether.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    } else if (score_info.modifier & SCORE_CENTER) {
	      logerrprint("Error: --score 'center' and 'no-mean-imputation' modifiers cannot be used\ntogether.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
            score_info.modifier |= SCORE_NO_MEAN_IMPUTATION;
	  } else if (!strcmp(argv[cur_arg + uii], "center")) {
	    if (score_info.modifier & SCORE_NO_MEAN_IMPUTATION) {
	      logerrprint("Error: --score 'center' and 'no-mean-imputation' modifiers cannot be used\ntogether.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
            score_info.modifier |= SCORE_CENTER;
	  } else if (!strcmp(argv[cur_arg + uii], "no-sum")) {
	    if (score_info.modifier & SCORE_SUM) {
	      logerrprint("Error: --score 'sum' and 'no-sum' modifiers conflict.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dosage_info.modifier |= DOSAGE_SCORE_NOSUM;
	  } else if (!strcmp(argv[cur_arg + uii], "include-cnt")) {
            dosage_info.modifier |= DOSAGE_SCORE_CNT;
	  } else if (!strcmp(argv[cur_arg + uii], "double-dosage")) {
	    dosage_info.modifier |= DOSAGE_SCORE_DOUBLE;
	  } else if (ujj == 3) {
            logerrprint("Error: --score takes at most three numeric parameters.\n");
            goto main_ret_INVALID_CMDLINE_A;
	  } else {
	    if (scan_posint_cappedx(argv[cur_arg + uii], MAXLINEBUFLEN / 2, (uint32_t*)&ii)) {
              sprintf(g_logbuf, "Error: Invalid --score parameter '%s'.\n", argv[cur_arg + uii]);
              goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (!ujj) {
              score_info.varid_col = ii;
	      score_info.allele_col = ii + 1;
	      score_info.effect_col = ii + 2;
	    } else if (ujj == 1) {
	      if ((uint32_t)ii == score_info.varid_col) {
	        logerrprint("Error: --score variant ID and allele column numbers cannot match.\n");
	        goto main_ret_INVALID_CMDLINE_A;
	      }
              score_info.allele_col = ii;
              score_info.effect_col = ii + 1;
	    } else {
	      if ((uint32_t)ii == score_info.varid_col) {
	        logerrprint("Error: --score variant ID and effect column numbers cannot match.\n");
	        goto main_ret_INVALID_CMDLINE_A;
	      } else if ((uint32_t)ii == score_info.allele_col) {
	        logerrprint("Error: --score allele and effect column numbers cannot match.\n");
	        goto main_ret_INVALID_CMDLINE_A;
	      }
              score_info.effect_col = ii;
	    }
            ujj++;
	  }
	}
        calculation_type |= CALC_SCORE;
      } else if (!memcmp(argptr2, "core-no-mean-imputation", 24)) {
	if (!(calculation_type & CALC_SCORE)) {
	  logerrprint("Error: --score-no-mean-imputation must be used with --score.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (score_info.modifier & SCORE_CENTER) {
	  logerrprint("Error: --score-no-mean-imputation cannot be used with --score 'center'\nmodifier.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --score-no-mean-imputation flag deprecated.  Use e.g.\n\"--score ... no-mean-imputation\".\n");
        score_info.modifier |= SCORE_NO_MEAN_IMPUTATION;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et-me-missing", 14)) {
	if (load_rare & LOAD_RARE_CNV) {
	  logerrprint("Error: --set-me-missing cannot be used with a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	misc_flags |= MISC_SET_ME_MISSING;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "how-tags", 9)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (strcmp(argv[cur_arg + 1], "all")) {
	  // no explicit flag for now, instead just let null show_tags_fname
	  // indicate 'all' mode;
	  retval = alloc_fname(&ld_info.show_tags_fname, argv[cur_arg + 1], argptr, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	} else {
          if (ld_info.modifier & LD_SHOW_TAGS_LIST_ALL) {
	    logerrprint("Error: --list-all cannot be used with \"--show-tags all\".\n");
            goto main_ret_INVALID_CMDLINE_A;
	  }
	  ld_info.modifier |= LD_SHOW_TAGS_LIST_ALL;
	}
        calculation_type |= CALC_SHOW_TAGS;
      } else if (!memcmp(argptr2, "egment", 7)) {
	logerrprint("Error: --segment flag retired.  Use another package, such as BEAGLE or\nGERMLINE, for this analysis.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "kato", 5)) {
	logerrprint("Error: --skato is not implemented yet.  Use e.g. PLINK/SEQ to perform this test\nfor now.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto main_ret_1;
      } else if (memcmp(argptr2, "ilent", 6)) {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 't':
      if (!memcmp(argptr2, "ail-pheno", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (makepheno_str) {
	  logerrprint("Error: --tail-pheno cannot be used with --make-pheno.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (scan_doublex(argv[cur_arg + 1], &tail_bottom)) {
	  sprintf(g_logbuf, "Error: Invalid --tail-pheno parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (param_ct == 1) {
	  tail_top = tail_bottom;
	} else {
	  if (scan_doublex(argv[cur_arg + 2], &tail_top)) {
	    sprintf(g_logbuf, "Error: Invalid --tail-pheno parameter '%s'.\n", argv[cur_arg + 2]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	if (tail_bottom > tail_top) {
	  logerrprint("Error: Ltop cannot be larger than Hbottom for --tail-pheno.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        filter_flags |= FILTER_FAM_REQ | FILTER_TAIL_PHENO;
      } else if (!memcmp(argptr2, "hreads", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_posint_defcapx(argv[cur_arg + 1], &g_thread_ct)) {
	  sprintf(g_logbuf, "Error: Invalid --threads parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (g_thread_ct > MAX_THREADS) {
	  LOGPRINTF("Note: Reducing --threads parameter to %u.  (If this is not large enough,\nrecompile with a larger MAX_THREADS setting.)\n", MAX_THREADS);
	  g_thread_ct = MAX_THREADS;
	} else if (known_procs == -1) {
	  // trigger BLAS/LAPACK warning
	  known_procs = 0;
	}
      } else if (!memcmp(argptr2, "ab", 3)) {
	logprint("Note: --tab flag deprecated.  Use \"--recode tab ...\".\n");
	if (recode_modifier & RECODE_DELIMX) {
	  logerrprint("Error: Multiple --recode delimiter modifiers.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	recode_modifier |= RECODE_TAB;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ranspose", 9)) {
	logprint("Note: --transpose flag deprecated.  Use \"--recode transpose ...\".\n");
	if (recode_modifier & RECODE_LGEN) {
	  logerrprint("Error: --recode 'transpose' and 'lgen' modifiers cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	recode_modifier |= RECODE_TRANSPOSE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "fam", 4)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > FNAMESIZE - 1) {
	  logerrprint("Error: --tfam filename prefix too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(famname, argv[cur_arg + 1], jj + 1);
	load_rare |= LOAD_RARE_TFAM;
      } else if (!memcmp(argptr2, "file", 5)) {
	if (load_params || (load_rare & (~LOAD_RARE_TFAM))) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  jj = strlen(argv[cur_arg + 1]);
	  if (jj > FNAMESIZE - 6) {
	    logerrprint("Error: --tfile filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(memcpya(pedname, argv[cur_arg + 1], jj), ".tped", 6);
	  if (!(load_rare & LOAD_RARE_TFAM)) {
	    memcpy(memcpya(famname, argv[cur_arg + 1], jj), ".tfam", 6);
	  }
	} else {
	  memcpy(pedname, PROG_NAME_STR ".tped", 11);
	  if (!(load_rare & LOAD_RARE_TFAM)) {
	    memcpy(famname, PROG_NAME_STR ".tfam", 11);
	  }
	}
	load_rare |= LOAD_RARE_TRANSPOSE;
      } else if (!memcmp(argptr2, "ped", 4)) {
	if (load_params || (load_rare & (~(LOAD_RARE_TRANSPOSE | LOAD_RARE_TFAM)))) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (!(load_rare & (LOAD_RARE_TRANSPOSE | LOAD_RARE_TFAM))) {
	  logerrprint("Error: --tped must be used with --tfam or --tfile.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > FNAMESIZE - 1) {
	  logerrprint("Error: --tped filename prefix too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(pedname, argv[cur_arg + 1], jj + 1);
	load_rare |= LOAD_RARE_TPED;
      } else if (!memcmp(argptr2, "o", 2)) {
	if ((!are_all_words_zero(chrom_info.chrom_mask, CHROM_MASK_INITIAL_WORDS)) || chrom_info.incl_excl_name_stack) {
	  logerrprint("Error: --to cannot be used with --autosome[-xy] or --[not-]chr.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (markername_snp) {
	  logerrprint("Error: --to cannot be used with --snp.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (snps_range_list.names) {
	  logerrprint("Error: --to cannot be used with --snps.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&markername_to, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if ((!memcmp(argptr2, "o-bp", 5)) || (!memcmp(argptr2, "o-kb", 5)) || (!memcmp(argptr2, "o-mb", 5))) {
	if (markername_snp && (!(filter_flags & FILTER_EXCLUDE_MARKERNAME_SNP))) {
	  logerrprint("Error: --to-bp/-kb/-mb cannot be used with --snp.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (snps_range_list.names) {
	  logerrprint("Error: --to-bp/-kb/-mb cannot be used with --snps.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (markername_to) {
	  logerrprint("Error: --to-bp/-kb/-mb cannot be used with --to.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if ((!markername_from) && (marker_pos_start == -1)) {
	  marker_pos_start = 0;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	cc = argptr2[2];
	if (cc == 'b') {
	  if (scan_uint_defcapx(argv[cur_arg + 1], (uint32_t*)&ii)) {
	    sprintf(g_logbuf, "Error: Invalid --to-bp parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	} else {
	  if (marker_pos_end != -1) {
	    logerrprint("Error: Multiple --to-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (scan_doublex(argv[cur_arg + 1], &dxx)) {
	    sprintf(g_logbuf, "Error: Invalid --to-kb/-mb parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	  dxx *= (cc == 'k')? 1000 : 1000000;
	  if (dxx < 0) {
	    ii = 0;
	  } else if (dxx > 2147483646) {
	    ii = 0x7ffffffe;
	  } else {
	    ii = (int32_t)(dxx * (1 + SMALL_EPSILON));
	  }
	}
	if (ii < marker_pos_start) {
	  marker_pos_end = marker_pos_start;
	  marker_pos_start = ii;
	} else {
	  marker_pos_end = ii;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP;
      } else if (!memcmp(argptr2, "rend", 5)) {
	if (model_modifier & MODEL_ASSOC) {
	  logerrprint("Error: --trend cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (model_modifier & MODEL_FISHER) {
	  logerrprint("Error: --trend cannot be used with --model fisher.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	  logerrprint("Error: --trend cannot be used with --model dom/rec/gen.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --trend flag deprecated.  Use \"--model trend-only ...\".\n");
	model_modifier |= MODEL_PTREND | MODEL_TRENDONLY;
        calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "hin", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &thin_keep_prob)) {
	  sprintf(g_logbuf, "Error: Invalid --thin variant retention probability '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        if (thin_keep_prob < (0.5 / 4294967296.0)) {
	  logerrprint("Error: --thin variant retention probability too small.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (thin_keep_prob >= (4294967295.5 / 4294967296.0)) {
	  if (scan_uint_defcapx(argv[cur_arg + 1], &uii)) {
	    logerrprint("Error: --thin variant retention probability too large.\n");
	  } else {
	    // VCFtools --thin = --bp-space...
	    logerrprint("Error: --thin variant retention probability too large.  (Did you mean\n--bp-space?)\n");
	  }
	  goto main_ret_INVALID_CMDLINE_A;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "hin-count", 10)) {
	if (thin_keep_prob != 1.0) {
	  logerrprint("Error: --thin cannot be used with --thin-count.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_uint_defcapx(argv[cur_arg + 1], &thin_keep_ct) || ((!thin_keep_ct) && (!(misc_flags & MISC_ALLOW_NO_VARS)))) {
	  sprintf(g_logbuf, "Error: Invalid --thin-count parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "hin-indiv", 10)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
        }
        if (scan_doublex(argv[cur_arg + 1], &thin_keep_sample_prob)) {
          sprintf(g_logbuf, "Error: Invalid --thin-indiv %s retention probability '%s'.\n", g_species_singular, argv[cur_arg + 1]);
          goto main_ret_INVALID_CMDLINE_WWA;
        }
        if (thin_keep_sample_prob < (0.5 / 4294967296.0)) {
          LOGERRPRINTF("Error: --thin-indiv %s retention probability too small.\n", g_species_singular);
          goto main_ret_INVALID_CMDLINE_A;
        } else if (thin_keep_sample_prob >= (4294967295.5 / 4294967296.0)) {
          LOGERRPRINTF("Error: --thin-indiv %s retention probability too large.\n", g_species_singular);
          goto main_ret_INVALID_CMDLINE_A;
        }
      } else if (!memcmp(argptr2, "hin-indiv-count", 16)) {
        if (thin_keep_sample_prob != 1.0) {
          logerrprint("Error: --thin-indiv cannot be used with --thin-indiv-count.\n");
          goto main_ret_INVALID_CMDLINE_WWA;
        }
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_2A;
        }
        if (scan_uint_defcapx(argv[cur_arg + 1], &thin_keep_sample_ct) || ((!thin_keep_sample_ct) && (!(misc_flags & MISC_ALLOW_NO_SAMPLES)))) {
          sprintf(g_logbuf, "Error: Invalid --thin-indiv-count parameter '%s'.\n", argv[cur_arg + 1]);
          goto main_ret_INVALID_CMDLINE_WWA;
        }
      } else if (!memcmp(argptr2, "ests", 5)) {
	if (!(calculation_type & CALC_GLM)) {
	  logerrprint("Error: --tests must be used with --linear or --logistic.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if ((param_ct == 1) && (!strcmp(argv[cur_arg + 1], "all"))) {
	  glm_modifier |= GLM_TEST_ALL;
	} else {
	  if (glm_modifier & GLM_TEST_ALL) {
	    logerrprint("Error: --test-all cannot be used with --tests.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  retval = parse_name_ranges(param_ct, '-', &(argv[cur_arg]), &tests_range_list, 1);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "est-all", 8)) {
	if (!(calculation_type & CALC_GLM)) {
	  logerrprint("Error: --test-all must be used with --linear or --logistic.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --test-all flag deprecated.  Use \"--tests all\".\n");
	glm_modifier |= GLM_TEST_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "wolocus", 8)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
          goto main_ret_INVALID_CMDLINE_2A;
	}
        if (alloc_string(&(epi_info.twolocus_mkr1), argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
        if (alloc_string(&(epi_info.twolocus_mkr2), argv[cur_arg + 2])) {
	  goto main_ret_NOMEM;
	}
	if (!strcmp(epi_info.twolocus_mkr1, epi_info.twolocus_mkr2)) {
          logerrprint("Error: --twolocus parameters cannot match.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        calculation_type |= CALC_EPI;
      } else if (!memcmp(argptr2, "est-missing", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
          if (!strcmp(argv[cur_arg + uii], "perm")) {
            if (testmiss_modifier & TESTMISS_MPERM) {
              logerrprint("Error: --test-missing 'mperm' and 'perm' cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
            testmiss_modifier |= TESTMISS_PERM;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
            if (testmiss_modifier & TESTMISS_PERM) {
              logerrprint("Error: --test-missing 'mperm' and 'perm' cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    } else if (testmiss_modifier & TESTMISS_MPERM) {
	      // when --mperm and --test-missing mperm=[] are used together
              logerrprint("Error: Duplicate --test-missing 'mperm' modifier.\n");
              goto main_ret_INVALID_CMDLINE;
	    }
	    if (scan_posint_defcapx(&(argv[cur_arg + uii][6]), &testmiss_mperm_val)) {
	      sprintf(g_logbuf, "Error: Invalid --test-missing mperm parameter '%s'.\n", &(argv[cur_arg + uii][6]));
              goto main_ret_INVALID_CMDLINE_WWA;
	    }
            testmiss_modifier |= TESTMISS_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
            testmiss_modifier |= TESTMISS_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "midp")) {
            testmiss_modifier |= TESTMISS_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
            logerrprint("Error: Improper --test-missing mperm syntax.  (Use \"--test-missing\nmperm=<value>\".)\n");
            goto main_ret_INVALID_CMDLINE;
	  } else {
            sprintf(g_logbuf, "Error: Invalid --test-missing parameter '%s'.\n", argv[cur_arg + uii]);
            goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
        calculation_type |= CALC_TESTMISS;
      } else if (!memcmp(argptr2, "est-mishap", 11)) {
        calculation_type |= CALC_TESTMISHAP;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "dt", 3)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 5)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "exact")) {
	    if (family_info.tdt_modifier & TDT_MIDP) {
	      logerrprint("Error: --tdt 'exact' modifier cannot be used with 'exact-midp'.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (family_info.tdt_modifier & TDT_POO) {
	      logerrprint("Error: --tdt parent-of-origin analysis does not currently support exact tests.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    family_info.tdt_modifier |= TDT_EXACT;
	  } else if (!strcmp(argv[cur_arg + uii], "exact-midp")) {
	    if ((family_info.tdt_modifier & (TDT_EXACT | TDT_MIDP)) == TDT_EXACT) {
	      logerrprint("Error: --tdt 'exact' modifier cannot be used with 'exact-midp'.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (family_info.tdt_modifier & TDT_POO) {
	      logerrprint("Error: --tdt parent-of-origin analysis does not currently support exact tests.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    family_info.tdt_modifier |= TDT_EXACT | TDT_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "poo")) {
	    if (family_info.tdt_modifier & TDT_EXACT) {
	      logerrprint("Error: --tdt parent-of-origin analysis does not currently support exact tests.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    family_info.tdt_modifier |= TDT_POO;
	  } else if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (family_info.tdt_modifier & TDT_MPERM) {
	      logerrprint("Error: --tdt 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    family_info.tdt_modifier |= TDT_PERM;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    family_info.tdt_modifier |= TDT_PERM_COUNT;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (family_info.tdt_modifier & TDT_PERM) {
	      logerrprint("Error: --tdt 'mperm' and 'perm' cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    } else if (family_info.tdt_modifier & TDT_MPERM) {
	      logerrprint("Error: Duplicate --tdt 'mperm' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (scan_posint_defcapx(&(argv[cur_arg + uii][6]), &family_info.tdt_mperm_val)) {
	      sprintf(g_logbuf, "Error: Invalid --tdt mperm parameter '%s'.\n", &(argv[cur_arg + uii][6]));
              goto main_ret_INVALID_CMDLINE_WWA;
	    }
            family_info.tdt_modifier |= TDT_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "parentdt1")) {
	    if (family_info.tdt_modifier & TDT_PARENPERM2) {
	      logerrprint("Error: --tdt 'parentdt1' and 'parentdt2' modifiers cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
	    family_info.tdt_modifier |= TDT_PARENPERM1;
	  } else if (!strcmp(argv[cur_arg + uii], "parentdt2")) {
	    if (family_info.tdt_modifier & TDT_PARENPERM1) {
	      logerrprint("Error: --tdt 'parentdt1' and 'parentdt2' modifiers cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
	    family_info.tdt_modifier |= TDT_PARENPERM2;
	  } else if (!strcmp(argv[cur_arg + uii], "pat")) {
	    if (family_info.tdt_modifier & TDT_POOPERM_MAT) {
	      logerrprint("Error: --tdt 'pat' and 'mat' modifiers cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
	    family_info.tdt_modifier |= TDT_POOPERM_PAT;
	  } else if (!strcmp(argv[cur_arg + uii], "parentdt2")) {
	    if (family_info.tdt_modifier & TDT_POOPERM_PAT) {
	      logerrprint("Error: --tdt 'pat' and 'mat' modifiers cannot be used together.\n");
              goto main_ret_INVALID_CMDLINE_A;
	    }
	    family_info.tdt_modifier |= TDT_POOPERM_MAT;
	  } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
            family_info.tdt_modifier |= TDT_SET_TEST;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    logerrprint("Error: Improper --tdt mperm syntax.  (Use \"--tdt mperm=<value>\".)\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --tdt parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	if (!(family_info.tdt_modifier & TDT_POO)) {
	  if (family_info.tdt_modifier & TDT_PERM) {
	    logerrprint("Error: Regular --tdt analysis no longer supports adaptive permutation.  Use an\nexact test instead.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  } else if (family_info.tdt_modifier & (TDT_POOPERM_PAT | TDT_POOPERM_MAT)) {
	    logerrprint("Error: --tdt 'pat'/'mat' must be used with parent-of-origin analysis.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  } else if ((family_info.tdt_modifier & (TDT_PARENPERM1 | TDT_PARENPERM2 | TDT_SET_TEST)) && (!(family_info.tdt_modifier & TDT_MPERM))) {
	    logerrprint("Error: --tdt 'parentdt1'/'parentdt2'/'set-test' requires permutation.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	} else {
	  if (family_info.tdt_modifier & (TDT_PARENPERM1 | TDT_PARENPERM2)) {
	    logerrprint("Error: --tdt 'parentdt1'/'parentdt2' cannot be used with parent-of-origin\nanalysis.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  } else if ((family_info.tdt_modifier & (TDT_POOPERM_PAT | TDT_POOPERM_MAT | TDT_SET_TEST)) && (!(family_info.tdt_modifier & (TDT_PERM | TDT_MPERM)))) {
	    logerrprint("Error: --tdt poo 'pat'/'mat'/'set-test' requires permutation.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	}
	calculation_type |= CALC_TDT;
      } else if (!memcmp(argptr2, "ag-kb", 6)) {
        if (!(calculation_type & CALC_SHOW_TAGS)) {
	  logerrprint("Error: --tag-kb must be used with --show-tags.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --tag-kb parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (dxx > 2147483.646) {
	  ld_info.show_tags_bp = 2147483646;
	} else {
	  ld_info.show_tags_bp = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
	}
      } else if (!memcmp(argptr2, "ag-r2", 6)) {
        if (!(calculation_type & CALC_SHOW_TAGS)) {
	  logerrprint("Error: --tag-r2 must be used with --show-tags.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0) || (dxx > 1)) {
	  sprintf(g_logbuf, "Error: Invalid --tag-r2 threshold '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	ld_info.show_tags_r2 = dxx;
      } else if (!memcmp(argptr2, "ag-mode2", 8)) {
        if (!(calculation_type & CALC_SHOW_TAGS)) {
	  logerrprint("Error: --tag-mode2 must be used with --show-tags.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (!(ld_info.show_tags_fname)) {
	  logerrprint("Error: --tag-mode2 cannot be used with \"--show-tags all\".\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	ld_info.modifier |= LD_SHOW_TAGS_MODE2;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "ucc", 4)) {
	UNSTABLE("tucc");
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "write-bed")) {
	    sprintf(g_logbuf, "Error: Invalid --tucc parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE;
	  }
	  family_info.tucc_bed = 1;
	} else {
	  logerrprint("Warning: --tucc without 'write-bed' is deprecated.\n");
	}
	calculation_type |= CALC_TUCC;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'u':
      if (!memcmp(argptr2, "nrelated-heritability", 22)) {
#ifdef NOLAPACK
        logerrprint("Error: --unrelated-heritability requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	goto main_ret_INVALID_CMDLINE;
#else
	UNSTABLE("unrelated-heritability");
	if (rel_info.modifier & REL_CALC_COV) {
	  logerrprint("Error: --unrelated-heritability flag cannot coexist with a covariance\nmatrix calculation.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (parallel_tot > 1) {
	  logerrprint("Error: --parallel and --unrelated-heritability cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (calculation_type & CALC_PCA) {
	  logerrprint("Error: --pca cannot be used with --unrelated-heritability.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (load_rare & (LOAD_RARE_GRM | LOAD_RARE_GRM_BIN)) {
	  if (calculation_type & CALC_REL_CUTOFF) {
	    logerrprint("Error: --unrelated-heritability + --grm-gz/--grm-bin cannot be used with\n--rel-cutoff.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  } else if (!phenoname) {
	    logerrprint("Error: --unrelated-heritability + --grm-gz/--grm-bin requires --pheno as well.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (!strcmp(argv[cur_arg + 1], "strict")) {
	    rel_info.modifier |= REL_UNRELATED_HERITABILITY_STRICT;
	    uii = 2;
	  } else {
	    uii = 1;
	  }
	  if (param_ct >= uii) {
	    if (scan_doublex(argv[cur_arg + uii], &rel_info.unrelated_herit_tol)) {
	      sprintf(g_logbuf, "Error: Invalid --unrelated-heritability EM tolerance parameter '%s'.\n", argv[cur_arg + uii]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (rel_info.unrelated_herit_tol <= 0.0) {
	      sprintf(g_logbuf, "Error: Invalid --unrelated-heritability EM tolerance parameter '%s'.\n", argv[cur_arg + uii]);
	      goto main_ret_INVALID_CMDLINE_WWA;
	    }
	    if (param_ct > uii) {
	      if (scan_doublex(argv[cur_arg + uii + 1], &rel_info.unrelated_herit_covg)) {
		sprintf(g_logbuf, "Error: Invalid --unrelated-heritability genomic covariance prior '%s'.\n", argv[cur_arg + uii + 1]);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      if ((rel_info.unrelated_herit_covg <= 0.0) || (rel_info.unrelated_herit_covg > 1.0)) {
		sprintf(g_logbuf, "Error: Invalid --unrelated-heritability genomic covariance prior '%s'.\n", argv[cur_arg + uii + 1]);
		goto main_ret_INVALID_CMDLINE_WWA;
	      }
	      if (param_ct == uii + 2) {
		if (scan_doublex(argv[cur_arg + uii + 2], &rel_info.unrelated_herit_covr)) {
		  sprintf(g_logbuf, "Error: Invalid --unrelated-heritability residual covariance prior '%s'.\n", argv[cur_arg + uii + 2]);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
		if ((rel_info.unrelated_herit_covr <= 0.0) || (rel_info.unrelated_herit_covr > 1.0)) {
		  sprintf(g_logbuf, "Error: Invalid --unrelated-heritability residual covariance prior '%s'.\n", argv[cur_arg + uii + 2]);
		  goto main_ret_INVALID_CMDLINE_WWA;
		}
	      } else {
		rel_info.unrelated_herit_covr = 1.0 - rel_info.unrelated_herit_covg;
	      }
	    }
	  }
	}
	calculation_type |= CALC_UNRELATED_HERITABILITY;
#endif
      } else if (!memcmp(argptr2, "pdate-alleles", 14)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&update_alleles_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_NODOSAGE | FILTER_NOCNV;
      } else if (!memcmp(argptr2, "pdate-chr", 10)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  logprint("--update-chr cannot be used with --cnv-make-map.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if ((misc_flags & MISC_MERGEX) || splitx_bound2) {
          logprint("--update-chr cannot be used with --split-x or --merge-x.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (!param_ct) {
	  update_map_modifier = 1;
	} else {
	  retval = alloc_2col(&update_chr, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "pdate-cm", 9)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  logerrprint("Error: --update-cm cannot be used with --cnv-make-map.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (cm_map_fname) {
	  logerrprint("Error: --update-cm cannot be used with --cm-map.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	if (update_map_modifier) {
	  logerrprint("Error: --update-map 'cm' modifier cannot be used with 'chr'.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (!param_ct) {
	  update_map_modifier = 2;
	} else {
	  retval = alloc_2col(&update_cm, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP;
      } else if (!memcmp(argptr2, "pdate-ids", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&update_ids_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "pdate-map", 10)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  logprint("--update-map cannot be used with --cnv-make-map.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (update_map_modifier) {
	  if (param_ct != 1) {
	    sprintf(g_logbuf, "Error: Multi-parameter --update-map cannot be used with deprecated\nparameter-free --update-%s.\n", (update_map_modifier == 1)? "chr" : "cm");
	    goto main_ret_INVALID_CMDLINE_2;
	  }
	  retval = alloc_2col((update_map_modifier == 1)? (&update_chr) : (&update_cm), &(argv[cur_arg + 1]), argptr, 1);
	  if (retval) {
	    goto main_ret_1;
	  }
	  LOGPRINTF("Note: --update-map <filename> + parameter-free --update-%s deprecated.  Use\n--update-%s <filename> instead.\n", (update_map_modifier == 1)? "chr" : "cm", (update_map_modifier == 1)? "chr" : "cm");
	  update_map_modifier = 0;
	} else {
	  retval = alloc_2col(&update_map, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP;
      } else if (!memcmp(argptr2, "pdate-name", 11)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  logprint("--update-name cannot be used with --cnv-make-map.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (update_chr) {
	  logerrprint("Error: --update-name cannot be used with --update-chr.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (update_cm) {
	  logerrprint("Error: --update-name cannot be used with --update-cm.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (!param_ct) {
	  if (!update_map) {
	    logerrprint("Error: Deprecated parameter-free --update-name cannot be used without\n--update-map.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
	  if ((update_map->colx != 2) || (update_map->colid != 1) || (update_map->skip) || (update_map->skipchar)) {
	    logerrprint("Error: Multi-parameter --update-map cannot be used with deprecated\nparameter-free --update-name.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  logprint("Note: --update-map <filename> + parameter-free --update-name deprecated.  Use\n--update-name <filename> instead.\n");
	  update_name = update_map;
	  update_map = nullptr;
	} else {
	  if (update_map) {
	    // no point in explaining the deprecated exception to this in the
	    // error message
	    logerrprint("Error: --update-name cannot be used with --update-map.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  retval = alloc_2col(&update_name, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
	filter_flags |= FILTER_BIM_REQ | FILTER_DOSAGEMAP;
      } else if (!memcmp(argptr2, "pdate-parents", 14)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (update_ids_fname) {
	  logerrprint("Error: --update-parents cannot be used with --update-ids.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&update_parents_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "pdate-sex", 10)) {
	if (update_ids_fname) {
	  logerrprint("Error: --update-sex cannot be used with --update-ids.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	retval = alloc_fname(&update_sex_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (param_ct == 2) {
	  if (scan_posint_defcapx(argv[cur_arg + 2], &update_sex_col)) {
	    sprintf(g_logbuf, "Error: Invalid --update-sex column parameter '%s'.  (This must be a positive integer.)\n", argv[cur_arg + 2]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "nbounded", 9)) {
	if (!(calculation_type & CALC_GENOME)) {
	  logerrprint("Error: --unbounded must be used with --genome.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logprint("Note: --unbounded flag deprecated.  Use \"--genome unbounded\".\n");
	genome_modifier |= GENOME_IBD_UNBOUNDED;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'v':
      if (!memcmp(argptr2, "if", 3)) {
	if (((!(calculation_type & CALC_GLM)) || (glm_modifier & GLM_LOGISTIC)) && (!((calculation_type & CALC_EPI) || (!(epi_info.modifier & EPI_REG))))) {
	  logerrprint("Error: --vif must be used with --linear/--epistasis.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &glm_vif_thresh)) {
	  sprintf(g_logbuf, "Error: Invalid --linear/--epistasis VIF threshold '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	if (glm_vif_thresh < 1.0) {
	  sprintf(g_logbuf, "Error: --linear/--epistasis VIF threshold '%s' too small (must be >= 1).\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "egas", 5)) {
	if (!set_info.fname) {
	  logerrprint("Error: --vegas must be used with --set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	logerrprint("Error: --vegas is currently under development.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto main_ret_1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "cf", 3)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = strlen(argv[cur_arg + 1]);
	if (uii > FNAMESIZE - 1) {
	  logerrprint("Error: --vcf filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(pedname, argv[cur_arg + 1], uii + 1);
	load_rare = LOAD_RARE_VCF;
      } else if (!memcmp(argptr2, "cf-min-qual", 12)) {
	if (!(load_rare & (LOAD_RARE_VCF | LOAD_RARE_BCF))) {
	  logerrprint("Error: --vcf-min-qual must be used with --vcf/--bcf.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_doublex(argv[cur_arg + 1], &vcf_min_qual) || (vcf_min_qual < 0.0)) {
	  sprintf(g_logbuf, "Error: Invalid --vcf-min-qual parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	vcf_min_qual *= 1 - SMALL_EPSILON;
      } else if (!memcmp(argptr2, "cf-filter", 10)) {
	if (!(load_rare & (LOAD_RARE_VCF | LOAD_RARE_BCF))) {
	  logerrprint("Error: --vcf-filter must be used with --vcf/--bcf.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (param_ct) {
	  retval = alloc_and_flatten(&vcf_filter_exceptions_flattened, &(argv[cur_arg + 1]), param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
	misc_flags |= MISC_VCF_FILTER;
      } else if (!memcmp(argptr2, "cf-min-gp", 10)) {
	if (!(load_rare & LOAD_RARE_VCF)) {
	  // er, probably want to support BCF too
	  logerrprint("Error: --vcf-min-gp must be used with --vcf.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &vcf_min_gp) || (vcf_min_gp <= 0.0) || (vcf_min_gp > 1.0)) {
	  sprintf(g_logbuf, "Error: Invalid --vcf-min-gp parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	vcf_min_gp *= 1 - SMALL_EPSILON;
      } else if (!memcmp(argptr2, "cf-min-gq", 10)) {
	if (!(load_rare & LOAD_RARE_VCF)) {
	  // probably want to support BCF too
	  logerrprint("Error: --vcf-min-gq must be used with --vcf.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (vcf_half_call) {
	  logerrprint("Error: --vcf-min-gq cannot currently be used with --vcf-half-call.\n");
	  retval = RET_CALC_NOT_YET_SUPPORTED;
	  goto main_ret_1;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &vcf_min_gq) || (vcf_min_gq < 0.0)) {
	  sprintf(g_logbuf, "Error: Invalid --vcf-min-gq parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	vcf_min_gq *= 1 - SMALL_EPSILON;
      } else if (!memcmp(argptr2, "cf-idspace-to", 14)) {
	if (!(load_rare & (LOAD_RARE_VCF | LOAD_RARE_BCF))) {
	  logerrprint("Error: --vcf-idspace-to must be used with --vcf/--bcf.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (id_delim == ' ') {
	  logerrprint("Error: --vcf-idspace-to cannot be used when the --id-delim character is space.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	vcf_idspace_to = extract_char_param(argv[cur_arg + 1]);
	if (!vcf_idspace_to) {
	  logerrprint("Error: --vcf-idspace-to parameter must be a single character.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (((unsigned char)vcf_idspace_to) <= ' ') {
	  logerrprint("Error: --vcf-idspace-to parameter must be a nonspace character.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "cf-half-call", 13)) {
        if (!(load_rare & LOAD_RARE_VCF)) {
	  logerrprint("Error: --vcf-half-call must be used with --vcf.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if ((!strcmp(argv[cur_arg + 1], "h")) || (!strcmp(argv[cur_arg + 1], "haploid"))) {
	  vcf_half_call = VCF_HALF_CALL_HAPLOID;
	} else if ((!strcmp(argv[cur_arg + 1], "m")) || (!strcmp(argv[cur_arg + 1], "missing"))) {
	  vcf_half_call = VCF_HALF_CALL_MISSING;
	} else if ((!strcmp(argv[cur_arg + 1], "e")) || (!strcmp(argv[cur_arg + 1], "error"))) {
	  vcf_half_call = VCF_HALF_CALL_ERROR;
	} else if ((!strcmp(argv[cur_arg + 1], "r")) || (!strcmp(argv[cur_arg + 1], "reference"))) {
	  vcf_half_call = VCF_HALF_CALL_REFERENCE;
	} else {
	  sprintf(g_logbuf, "Error: '%s' is not a valid mode for --vcf-half-call.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
      } else if (!memcmp(argptr2, "cf-require-gt", 14)) {
        if (!(load_rare & (LOAD_RARE_VCF | LOAD_RARE_BCF))) {
	  logerrprint("Error: --vcf-require-gt must be used with --vcf/--bcf.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	misc_flags |= MISC_VCF_REQUIRE_GT;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'w':
      if (!memcmp(argptr2, "rite-snplist", 13)) {
	calculation_type |= CALC_WRITE_SNPLIST;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-var-ranges", 16)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        if (scan_posint_defcapx(argv[cur_arg + 1], &write_var_range_ct)) {
	  sprintf(g_logbuf, "Error: Invalid --write-var-ranges parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	calculation_type |= CALC_WRITE_VAR_RANGES;
      } else if (!memcmp(argptr2, "indow", 6)) {
        if (!markername_snp) {
	  logerrprint("Error: --window must be used with --snp or --exclude-snp.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (scan_doublex(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(g_logbuf, "Error: Invalid --window parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
        dxx *= 500;
	if (dxx > 2147483646) {
	  snp_window_size = 0x7ffffffe;
	} else {
	  snp_window_size = (int32_t)(dxx * (1 + SMALL_EPSILON));
	}
	// no need to set filter_flags due to --snp/--exclude-snp req.
      } else if (!memcmp(argptr2, "ithin", 6)) {
        if (loop_assoc_fname) {
	  logerrprint("Error: --within cannot be used with --loop-assoc.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (misc_flags & MISC_FAMILY_CLUSTERS) {
          logerrprint("Error: --within cannot be used with --family.\n");
          goto main_ret_INVALID_CMDLINE_A;
	}
	if ((calculation_type & CALC_FREQ) && (misc_flags & MISC_FREQ_COUNTS)) {
	  logerrprint("Error: --within cannot be used with \"--freq counts\".\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	uii = 1;
	if (param_ct == 2) {
	  if ((strlen(argv[cur_arg + 1]) == 7) && (!memcmp(argv[cur_arg + 1], "keep-", 5)) && match_upper(&(argv[cur_arg + 1][5]), "NA")) {
	    uii = 2;
	  } else if ((strlen(argv[cur_arg + 2]) != 7) || memcmp(argv[cur_arg + 2], "keep-", 5) || (!match_upper(&(argv[cur_arg + 2][5]), "NA"))) {
            logerrprint("Error: Invalid --within parameter sequence.\n");
	    goto main_ret_INVALID_CMDLINE_A;
	  }
          misc_flags |= MISC_LOAD_CLUSTER_KEEP_NA;
	}
	retval = alloc_fname(&cluster.fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	filter_flags |= FILTER_FAM_REQ;
      } else if (!memcmp(argptr2, "ith-phenotype", 14)) {
	if (!covar_fname) {
	  logerrprint("Error: --with-phenotype cannot be used without --covar.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "no-parents")) {
	    write_covar_modifier |= WRITE_COVAR_NO_PARENTS;
	  } else if (!strcmp(argv[cur_arg + uii], "no-sex")) {
	    if (write_covar_modifier & WRITE_COVAR_FEMALE_2) {
	      logerrprint("Error: --with-phenotype 'female-2' modifier cannot be used with 'no-sex'.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    write_covar_modifier |= WRITE_COVAR_NO_SEX;
	  } else if (!strcmp(argv[cur_arg + uii], "female-2")) {
	    if (write_covar_modifier & WRITE_COVAR_NO_SEX) {
	      logerrprint("Error: --with-phenotype 'female-2' modifier cannot be used with 'no-sex'.\n");
	      goto main_ret_INVALID_CMDLINE_A;
	    }
	    write_covar_modifier |= WRITE_COVAR_FEMALE_2;
	  } else {
	    sprintf(g_logbuf, "Error: Invalid --with-phenotype parameter '%s'.\n", argv[cur_arg + uii]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
	}
	write_covar_modifier |= WRITE_COVAR_PHENO;
      } else if (!memcmp(argptr2, "ith-reference", 14)) {
	if ((recode_modifier & RECODE_TYPEMASK) != RECODE_LGEN) {
	  logerrprint("Error: --with-reference must be used with --recode lgen.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --with-reference flag deprecated.  Use \"--recode lgen-ref\" instead.\n");
	recode_modifier += RECODE_LGEN_REF - RECODE_LGEN;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-covar", 11)) {
	if (calculation_type & (CALC_MAKE_BED | CALC_MAKE_FAM | CALC_RECODE)) {
	  logerrprint("Error: --write-covar cannot be used with --make-bed/--make-just-fam/--recode.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (!covar_fname) {
	  logerrprint("Error: --write-covar cannot be used without --covar.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        calculation_type |= CALC_WRITE_COVAR;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-cluster", 13)) {
	if ((!cluster.fname) && (!(misc_flags & MISC_FAMILY_CLUSTERS))) {
	  logerrprint("Error: --write-cluster must be used with --within/--family.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "omit-unassigned")) {
	    sprintf(g_logbuf, "Error: Invalid --write-cluster parameter '%s'.\n", argv[cur_arg + 1]);
	    goto main_ret_INVALID_CMDLINE_WWA;
	  }
          misc_flags |= MISC_WRITE_CLUSTER_OMIT_UNASSIGNED;
	}
        calculation_type |= CALC_WRITE_CLUSTER;
      } else if (!memcmp(argptr2, "rite-set", 9)) {
	if (!set_info.fname) {
	  logerrprint("Error: --write-set must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        calculation_type |= CALC_WRITE_SET;
	set_info.modifier |= SET_WRITE_LIST;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-set-r2", 11)) {
        if (!set_info.fname) {
	  logerrprint("Error: --write-set-r2 must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        set_info.modifier |= SET_R2_WRITE;
        logprint("Note: --write-set-r2 flag deprecated.  Use \"--set-r2 write\".\n");
        goto main_param_zero;
      } else if (!memcmp(argptr2, "ith-freqs", 10)) {
	if (!(calculation_type & CALC_LD)) {
	  logerrprint("Error: --with-freqs must be used with --r/--r2.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
        ld_info.modifier |= LD_WITH_FREQS;
	if (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_MATRIX_BIN | LD_MATRIX_BIN4 | LD_MATRIX_SPACES)) {
	  goto main_r2_matrix_conflict;
	}
	logprint("Note: --with-freqs flag deprecated.  Use e.g. \"--r2 with-freqs\".\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-dosage", 12)) {
	if (!(dosage_info.modifier & DOSAGE_GLM)) {
	  if (dosage_info.modifier & DOSAGE_OCCUR) {
	    logerrprint("Error: --write-dosage cannot be used with \"--dosage occur\".\n");
	  } else if (dosage_info.modifier & DOSAGE_SCORE) {
	    logerrprint("Error: --write-dosage cannot be used with --score.\n");
	  } else {
	    logerrprint("Error: --write-dosage must be used with --dosage.\n");
	  }
          goto main_ret_INVALID_CMDLINE_A;
	} else if ((glm_modifier & GLM_STANDARD_BETA) || (dosage_info.modifier & DOSAGE_SEX)) {
	  logerrprint("Error: --write-dosage cannot be used with --dosage association analysis flags.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	dosage_info.modifier += (DOSAGE_WRITE - DOSAGE_GLM);
        goto main_param_zero;
      } else if (!memcmp(argptr2, "hap", 4)) {
        goto main_hap_disabled_message;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'x':
      if (!memcmp(argptr2, "chr-model", 10)) {
	if (!(calculation_type & CALC_GLM)) {
	  logerrprint("Error: --xchr-model must be used with --linear or --logistic.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_DOMINANT | GLM_RECESSIVE)) {
	  sprintf(g_logbuf, "Error: --xchr-model cannot be used with --%s %s.\n", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", (glm_modifier & GLM_GENOTYPIC)? "genotypic" : ((glm_modifier & GLM_HETHOM)? "hethom" : ((glm_modifier & GLM_DOMINANT)? "dominant" : "recessive")));
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
	if ((argv[cur_arg + 1][1] != '\0') || (argv[cur_arg + 1][0] < '0') || (argv[cur_arg + 1][0] > '3')) {
	  sprintf(g_logbuf, "Error: Invalid --xchr-model parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_WWA;
	}
	glm_xchr_model = (uint32_t)(argv[cur_arg + 1][0] - '0');
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    case 'z':
      if (!memcmp(argptr2, "ero-cluster", 12)) {
	if ((!cluster.fname) && (!(misc_flags & MISC_FAMILY_CLUSTERS))) {
	  logerrprint("Error: --zero-cluster must be used with --within/--family.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	} else if ((calculation_type != CALC_MAKE_BED) || (geno_thresh != 1.0) || (mind_thresh != 1.0) || (hwe_thresh != 0.0) || (min_maf != 0.0) || (max_maf != 0.5)) {
	  // prevent old pipelines from silently breaking
	  logerrprint("Error: --zero-cluster must now be used with --make-bed, no other output\ncommands, and no genotype-based filters.\n");
	  goto main_ret_INVALID_CMDLINE_A;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_2A;
	}
        // bugfix (9 Apr 2018): had wrong array index
	retval = alloc_fname(&cluster.zerofname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ero-cms", 8)) {
	filter_flags |= FILTER_BIM_REQ | FILTER_ZERO_CMS | FILTER_NODOSAGE | FILTER_NOCNV;
        goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      }
      break;

    default:
      goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;

    main_param_zero:
      if (param_ct) {
        sprintf(g_logbuf, "Error: --%s doesn't accept parameters.\n", argptr);
	goto main_ret_INVALID_CMDLINE_2A;
      }
    }
  } while ((++cur_flag) < flag_ct);
  if (!outname_end) {
    outname_end = &(outname[5]);
  }

  // command-line restrictions which don't play well with alphabetical order
  if (load_rare) {
    if (load_rare & (LOAD_RARE_GRM | LOAD_RARE_GRM_BIN)) {
      if ((!(calculation_type & (CALC_REL_CUTOFF | CALC_UNRELATED_HERITABILITY))) || (calculation_type & (~(CALC_REL_CUTOFF | CALC_RELATIONSHIP | CALC_UNRELATED_HERITABILITY)))) {
	if (load_rare == LOAD_RARE_GRM) {
	  logerrprint("Error: --grm-gz currently must be used with --rel-cutoff (possibly combined\nwith --make-grm-gz/--make-grm-bin) or --unrelated-heritability.\n");
	} else {
	  logerrprint("Error: --grm-bin currently must be used with --rel-cutoff (possibly combined\nwith --make-grm-gz/--make-grm-bin) or --unrelated-heritability.\n");
	}
	goto main_ret_INVALID_CMDLINE_A;
      }
    }
    if (!mperm_val) {
      if ((cnv_calc_type & CNV_SAMPLE_PERM) && (!cnv_sample_mperms)) {
	logerrprint("Error: --cnv-indiv-perm requires a permutation count.\n");
	goto main_ret_INVALID_CMDLINE_A;
      } else if ((cnv_calc_type & CNV_TEST_REGION) && (!cnv_test_region_mperms)) {
	logerrprint("Error: --cnv-test-region requires a permutation count.\n");
	goto main_ret_INVALID_CMDLINE_A;
      }
    }
  } else if ((load_params & LOAD_PARAMS_BFILE_ALL) && (load_params != LOAD_PARAMS_BFILE_ALL)) {
    if ((calculation_type & (~(CALC_ONLY_BIM | CALC_ONLY_FAM))) || (filter_flags & FILTER_ALL_REQ)) {
      logerrprint("Error: A full .bed + .bim + .fam fileset is required for this.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((!mapname[0]) && ((calculation_type & (~CALC_ONLY_FAM)) || (filter_flags & FILTER_BIM_REQ))) {
      logerrprint("Error: A .bim file is required for this.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((!famname[0]) && ((calculation_type & (~CALC_ONLY_BIM)) || (filter_flags & FILTER_FAM_REQ))) {
      logerrprint("Error: A .fam file is required for this.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  }
  if (sample_sort && (calculation_type & (~(CALC_MERGE | CALC_MAKE_BED)))) {
    logerrprint("Error: --indiv-sort only affects --make-bed and --merge/--bmerge/--merge-list.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if ((cnv_intersect_filter_type & CNV_COUNT) && (!(cnv_calc_type & (CNV_SAMPLE_PERM | CNV_ENRICHMENT_TEST)))) {
    logerrprint("Error: --cnv-count must be used with --cnv-indiv-perm or --cnv-enrichment-test.\n");
    goto main_ret_INVALID_CMDLINE;
  }
  if (!phenoname) {
    if ((filter_flags & FILTER_PRUNE) && (!(fam_cols & FAM_COL_6))) {
      logerrprint("Error: --prune and --no-pheno cannot coexist without an alternate phenotype\nfile.\n");
      goto main_ret_INVALID_CMDLINE_A;
    } else if (pheno_modifier & PHENO_ALL) {
      logerrprint("Error: --all-pheno must be used with --pheno.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  } else if (read_dists_fname && (!(calculation_type & (CALC_CLUSTER | CALC_IBS_TEST | CALC_GROUPDIST | CALC_NEIGHBOR | CALC_REGRESS_DISTANCE)))) {
    logerrprint("Error: --read-dists cannot be used without --cluster, --ibs-test/--groupdist,\n--neighbour, or --regress-distance.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if ((cluster.ppc != 0.0) && (!read_genome_fname) && (calculation_type & (CALC_DISTANCE))) {
    logerrprint("Error: --ppc cannot be used with --distance without --read-genome.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if ((calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) && (((calculation_type & CALC_DISTANCE) && (dist_calc_type & (DISTANCE_1_MINUS_IBS | DISTANCE_ALCT))) || (calculation_type & CALC_PLINK1_DISTANCE_MATRIX)) && (!(calculation_type & CALC_GENOME))) {
    // actually allow this for now if --genome present, since it auto-clobbers
    // the wrong-unit distance matrix
    logerrprint("Error: --cluster and --neighbour cannot be used with non-IBS distance matrix\ncalculations.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if (matrix_flag_state == 1) {
    calculation_type |= CALC_PLINK1_IBS_MATRIX;
  }
  if (calculation_type & CALC_PLINK1_IBS_MATRIX) {
    if (distance_wts_fname || (distance_exp != 0.0)) {
      logerrprint("Error: --ibs-matrix cannot be used with --distance-wts.\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (dist_calc_type & DISTANCE_IBS) {
      logerrprint("Error: --ibs-matrix cannot be used with \"--distance ibs\".\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (read_genome_fname && (cluster.ppc == 0.0)) {
      logerrprint("Error: --read-genome is pointless with --ibs-matrix unless --ppc is also\npresent.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (read_dists_fname) {
      logerrprint("Error: --read-dists cannot be used with a distance matrix calculation.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (parallel_tot > 1) {
      logerrprint("Error: --parallel and --ibs-matrix cannot be used together.  Use --distance\ninstead.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  }
  if (distance_wts_fname && (!(calculation_type & (CALC_DISTANCE | CALC_RELATIONSHIP)))) {
    logerrprint("Error: --distance-wts must be used with --distance, --make-rel, --make-grm-bin,\nor --make-grm-gz.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if ((parallel_tot > 1) && (!(calculation_type & (CALC_LD | CALC_DISTANCE | CALC_GENOME | CALC_RELATIONSHIP)))) {
    if ((!(calculation_type & CALC_EPI)) || (!(epi_info.modifier & (EPI_FAST | EPI_REG)))) {
      logerrprint("Error: --parallel only affects --r/--r2, --distance, --genome, --make-rel,\n--make-grm-gz/--make-grm-bin, and --epistasis/--fast-epistasis.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  }
  if (update_map_modifier) {
    sprintf(g_logbuf, "Error: Deprecated parameter-free --update-%s cannot be used without\n--update-map.\n", (update_map_modifier == 1)? "chr" : "cm");
    goto main_ret_INVALID_CMDLINE_2A;
  }
  if (((misc_flags & (MISC_FILL_MISSING_A2 | MISC_MERGEX | MISC_SET_ME_MISSING)) || splitx_bound2 || update_chr) && (((load_rare == LOAD_RARE_CNV) && (cnv_calc_type != CNV_WRITE)) || ((load_rare != LOAD_RARE_CNV) && (calculation_type != CALC_MAKE_BED)))) {
    sprintf(g_logbuf, "Error: --merge-x/--split-x/--update-chr/--set-me-missing/--fill-missing-a2\nmust be used with --%s and no other commands.\n", (load_rare == LOAD_RARE_CNV)? "cnv-write" : "make-bed");
    goto main_ret_INVALID_CMDLINE_2A;
  }
  if (load_rare == LOAD_RARE_CNV) {
    if (filter_flags & FILTER_NOCNV) {
      logerrprint("Error: .cnv fileset specified with incompatible filtering flag(s).  (Check if\nthere is a --cnv-... flag with the functionality you're looking for.)\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  } else if (load_rare == LOAD_RARE_DOSAGE) {
    if (filter_flags & FILTER_NODOSAGE) {
      logerrprint("Error: --dosage used with incompatible filtering flag(s).\n");
      goto main_ret_INVALID_CMDLINE_A;
    } else if ((!mapname[0]) && (filter_flags & FILTER_DOSAGEMAP)) {
      logerrprint("Error: --dosage cannot be used with variant filters unless a .map file is\nspecified.\n");
      goto main_ret_INVALID_CMDLINE_A;
    } else if (!famname[0]) {
      logerrprint("Error: --dosage must be used with --fam.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  }
  if ((family_info.mendel_modifier & MENDEL_DUOS) && (!(calculation_type & CALC_MENDEL)) && (!(family_info.mendel_modifier & MENDEL_FILTER)) && (!(misc_flags & MISC_SET_ME_MISSING))) {
    logerrprint("Error: --mendel-duos must be used with --me/--mendel/--set-me-missing.\n");
    goto main_ret_INVALID_CMDLINE;
  } else if ((family_info.mendel_modifier & MENDEL_MULTIGEN) && (!(calculation_type & (CALC_MENDEL | CALC_TDT | CALC_DFAM | CALC_QFAM))) && (!(family_info.mendel_modifier & MENDEL_FILTER)) && (!(misc_flags & MISC_SET_ME_MISSING))) {
    logerrprint("Error: --mendel-multigen must be used with --me, --mendel, --set-me-missing, or\nan association test which checks for Mendel errors.\n");
    goto main_ret_INVALID_CMDLINE;
  }
  if (flip_subset_fname && (load_rare || (calculation_type != CALC_MAKE_BED) || (min_maf != 0.0) || (max_maf != 0.5) || (hwe_thresh != 0.0))) {
    logerrprint("Error: --flip-subset must be used with --flip, --make-bed, and no other\ncommands or MAF-based filters.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if (calculation_type & CALC_RECODE) {
    if (recode_modifier & RECODE_23) {
      if (chrom_info.species != SPECIES_HUMAN) {
	logerrprint("Error: --recode 23 can only be used on human data.\n");
	goto main_ret_INVALID_CMDLINE;
      }
      if ((recode_modifier & RECODE_23) && (misc_flags & MISC_ALLOW_EXTRA_CHROMS) && (!chrom_info.zero_extra_chroms)) {
	logerrprint("Error: --allow-extra-chr requires the '0' modifier when used with --recode 23.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if ((recode_modifier & RECODE_BEAGLE) && chrom_info.zero_extra_chroms) {
        logerrprint("Error: --allow-extra-chr cannot have the '0' modifier when used with\n--recode beagle.\n");
	goto main_ret_INVALID_CMDLINE;
      }
    } else if ((recode_modifier & RECODE_BIMBAM) && (misc_flags & MISC_ALLOW_EXTRA_CHROMS) && (!chrom_info.zero_extra_chroms)) {
      logerrprint("Error: --allow-extra-chr requires the '0' modifier when used with\n--recode bimbam.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (sex_missing_pheno & MUST_HAVE_SEX) {
    if (load_rare & LOAD_RARE_CNV) {
      if (!(cnv_calc_type & CNV_WRITE)) {
        logerrprint("Error: --must-have-sex must be used with --cnv-write.\n");
        goto main_ret_INVALID_CMDLINE;
      }
    } else {
      if (!(calculation_type & (CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_MAKE_FAM | CALC_RECODE))) {
        logerrprint("Error: --must-have-sex must be used with --make-bed, --make-just-fam, --recode,\nor --write-covar.\n");
        goto main_ret_INVALID_CMDLINE;
      }
    }
  }
  if (misc_flags & MISC_IMPUTE_SEX) {
    if (!(calculation_type & (CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_RECODE))) {
      logerrprint("Error: --impute-sex must be used with --make-bed/--recode/--write-covar.\n");
      goto main_ret_INVALID_CMDLINE_A;
    } else if (calculation_type & (~(CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_RECODE | CALC_SEXCHECK))) {
      logerrprint("Error: --impute-sex cannot be used with any commands other than\n--make-bed/--recode/--write-covar.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  }
  if (cluster.qmatch_fname && (!cluster.qt_fname)) {
    logerrprint("Error: --qt must be used with --qmatch.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if (mwithin_col && (!loop_assoc_fname) && (!cluster.fname)) {
    logerrprint("Error: --mwithin must be used with --within/--loop-assoc.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if ((!cluster.fname) && (!(misc_flags & MISC_FAMILY_CLUSTERS))) {
    if (cluster.keep_fname) {
      logerrprint("Error: --keep-clusters must be used with --within/--family.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.keep_flattened) {
      logerrprint("Error: --keep-cluster-names must be used with --within/--family.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.remove_fname) {
      logerrprint("Error: --remove-clusters must be used with --within/--family.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.remove_flattened) {
      logerrprint("Error: --remove-cluster-names must be used with --within/--family.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (rel_info.pca_cluster_names_flattened) {
      logerrprint("Error: --pca-cluster-names must be used with --within/--family.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (rel_info.pca_clusters_fname) {
      logerrprint("Error: --pca-clusters must be used with --within/--family.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (calculation_type & CALC_CMH) {
      logerrprint("Error: --mh/--bd/--mh2 must be used with --within/--family.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (calculation_type & CALC_HOMOG) {
      logerrprint("Error: --homog must be used with --within/--family.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if ((calculation_type & CALC_FST) && (!(misc_flags & MISC_FST_CC))) {
      logerrprint("Error: --fst should be used with --within, unless the 'case-control' modifier\nis specified.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (!set_info.fname) {
    if (set_info.modifier) {
      if (set_info.modifier & SET_COMPLEMENTS) {
	if (set_info.merged_set_name) {
	  logerrprint("Error: --make-set-complement-all must be used with --set/--make-set.\n");
	} else {
	  logerrprint("Error: --complement-sets must be used with --set/--make-set.\n");
	}
      } else { // only remaining possibility for now
        logerrprint("Error: --gene-all must be used with --set/--make-set.\n");
      }
      goto main_ret_INVALID_CMDLINE;
    } else if (set_info.genekeep_flattened) {
      logerrprint("Error: --gene must be used with --set/--make-set.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (model_modifier & MODEL_SET_TEST) {
      logerrprint("Error: --assoc/--model set-test must be used with --set/--make-set.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (glm_modifier & GLM_SET_TEST) {
      logerrprint("Error: --linear/--logistic set-test must be used with --set/--make-set.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (family_info.tdt_modifier & TDT_SET_TEST) {
      logerrprint("Error: --tdt set-test must be used with --set/--make-set.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.modifier & CLUSTER_CMH_SET_TEST) {
      logerrprint("Error: --mh/--bd set-test must be used with --set/--make-set.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  } else {
    uii = 0;
    if (model_modifier & MODEL_SET_TEST) {
      if ((!(model_modifier & MODEL_PERM)) && (!model_mperm_val)) {
        logerrprint("Error: --assoc/--model set-test requires permutation.\n");
        goto main_ret_INVALID_CMDLINE_A;
      } else if (model_modifier & MODEL_FISHER) {
	logerrprint("Error: --assoc/--model set-test cannot be used with Fisher's exact test.\n");
        goto main_ret_INVALID_CMDLINE_A;
      } else if (model_modifier & MODEL_PGEN) {
	logerrprint("Error: --model set-test cannot be used with 2df genotypic chi-square stats.\n");
	goto main_ret_INVALID_CMDLINE_A;
      } else if (model_modifier & MODEL_LIN) {
	logerrprint("Error: --assoc set-test does not currently support the Lin statistic; contact\nthe developers if you need this combination.\n");
	goto main_ret_INVALID_CMDLINE_A;
      }
      uii = 1;
    }
    if (glm_modifier & GLM_SET_TEST) {
      if ((!(glm_modifier & GLM_PERM)) && (!glm_mperm_val)) {
        logerrprint("Error: --linear/--logistic set-test requires permutation.\n");
        goto main_ret_INVALID_CMDLINE_A;
      } else if ((glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_TEST_ALL)) || tests_range_list.name_ct) {
	logerrprint("Error: --linear/--logistic set-test cannot be used with joint tests.\n");
	goto main_ret_INVALID_CMDLINE_A;
      }
      uii = 1;
    }
    if (family_info.tdt_modifier & TDT_SET_TEST) {
      if (!(family_info.tdt_modifier & (TDT_PERM | TDT_MPERM))) {
        logerrprint("Error: --tdt set-test requires permutation.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      logerrprint("Error: --tdt set-test is currently under development.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto main_ret_1;
      uii = 1;
    }
    if (family_info.dfam_modifier & DFAM_SET_TEST) {
      if (!(family_info.dfam_modifier & (DFAM_PERM | DFAM_MPERM))) {
        logerrprint("Error: --dfam set-test requires permutation.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      logerrprint("Error: --dfam set-test is currently under development.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto main_ret_1;
      uii = 1;
    }
    if (cluster.modifier & CLUSTER_CMH_SET_TEST) {
      if (!(cluster.modifier & (CLUSTER_CMH_PERM | CLUSTER_CMH_MPERM))) {
        logerrprint("Error: --mh/--bd set-test requires permutation.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      logerrprint("Error: --mh/--bd set-test is currently under development.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto main_ret_1;
      uii = 1;
    }
    if (mtest_adjust & ADJUST_LAMBDA) {
      if (set_info.set_test_lambda == 0.0) {
	// backward compatibility: --lambda and --set-test-lambda weren't
	// distinct in 1.07
	logprint("Note: set test + --lambda is deprecated.  Use --set-test-lambda instead.\n");
	set_info.set_test_lambda = adjust_lambda;
	if (!(mtest_adjust & 1)) {
	  mtest_adjust = 0;
	  adjust_lambda = 0.0;
	}
      } else {
	logprint("Note: Set test --adjust ignores --lambda.\n");
      }
    }
  }
  if ((family_info.tdt_modifier & (TDT_POO | TDT_PARENPERM1 | TDT_PARENPERM2 | TDT_POOPERM_PAT | TDT_POOPERM_MAT)) && (!(calculation_type & CALC_TDT))) {
    logerrprint("Error: --poo/--parentdt1/--parentdt2/--pat/--mat must be used with --tdt.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if (calculation_type & CALC_FLIPSCAN) {
    if (ld_info.flipscan_window_size == 0xffffffffU) {
      if (ld_info.window_size != 0xffffffffU) {
	logprint("Note: --ld-window + --flip-scan combination deprecated.  Use --flip-scan-window\nwith --flip-scan instead.\n");
	ld_info.flipscan_window_size = ld_info.window_size;
      } else {
	ld_info.flipscan_window_size = 10;
      }
    }
    if (ld_info.flipscan_window_bp == 0xffffffffU) {
      if (ld_info.window_bp != 0xffffffffU) {
        logprint("Note: --ld-window-kb + --flip-scan combination deprecated.  Use\n--flip-scan-window-kb with --flip-scan instead.\n");
        ld_info.flipscan_window_bp = ld_info.window_bp;
      } else {
	ld_info.flipscan_window_bp = 1000000;
      }
    }
  }
  if (calculation_type & CALC_BLOCKS) {
    if (ld_info.blocks_recomb_highci > ld_info.blocks_strong_highci) {
      logerrprint("Error: --blocks-recomb-highci value cannot be larger than\n--blocks-strong-highci value.\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (ld_info.blocks_max_bp == 0xffffffffU) {
      if (ld_info.window_bp != 0xffffffffU) {
        logprint("Note: --ld-window-kb + --blocks combination deprecated.  Use --blocks-max-kb\nwith --blocks instead.\n");
        ld_info.blocks_max_bp = ld_info.window_bp;
      } else {
        ld_info.blocks_max_bp = 200000;
      }
    }
  }
  if ((!(calculation_type & CALC_LD)) || ((calculation_type & CALC_LD) && (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_INTER_CHR)))) {
    if ((ld_info.snpstr || ld_info.snps_rl.name_ct) && (!(ld_info.modifier & LD_INTER_CHR))) {
      if (calculation_type & CALC_LD) {
	logerrprint("Error: --ld-snp/--ld-snps/--ld-snp-list cannot be used with the --r/--r2 matrix\noutput modifiers.\n");
      } else {
        logerrprint("Error: --ld-snp/--ld-snps/--ld-snp-list must be used with --r/--r2.\n");
      }
      goto main_ret_INVALID_CMDLINE_A;
    } else if (ld_info.window_size != 0xffffffffU) {
      if (calculation_type & CALC_LD) {
	logerrprint("Error: --ld-window flag cannot be used with the --r/--r2 'inter-chr' or matrix\noutput modifiers.\n");
        goto main_ret_INVALID_CMDLINE_A;
      } else if (!(calculation_type & (CALC_BLOCKS | CALC_FLIPSCAN))) {
        logerrprint("Error: --ld-window flag must be used with --r/--r2.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    } else if (ld_info.window_bp != 0xffffffffU) {
      if (calculation_type & CALC_LD) {
	logerrprint("Error: --ld-window-kb flag cannot be used with the --r/--r2 'inter-chr' or\nmatrix output modifiers.\n");
        goto main_ret_INVALID_CMDLINE_A;
      } else if (!(calculation_type & CALC_BLOCKS)) {
        logerrprint("Error: --ld-window-kb flag must be used with --r/--r2.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    } else if ((ld_info.window_r2 != 0.2) && (!(ld_info.modifier & LD_INTER_CHR))) {
      if (!(ld_info.modifier & LD_R2)) {
        logerrprint("Error: --ld-window-r2 flag must be used with --r2.\n");
        goto main_ret_INVALID_CMDLINE;
      } else {
	logerrprint("Error: --ld-window-r2 flag cannot be used with the --r2 matrix output modifiers.\n");
	goto main_ret_INVALID_CMDLINE_A;
      }
    }
  } else {
    if (ld_info.window_size == 0xffffffffU) {
      ld_info.window_size = 10;
    }
    if (ld_info.window_bp == 0xffffffffU) {
      ld_info.window_bp = 1000000;
    }
  }
  if ((ld_info.modifier & LD_DPRIME) && (!(calculation_type & CALC_LD))) {
    logerrprint("Error: --D/--dprime must be used with --r/--r2.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }

  // --from-bp/-kb/-mb without any --to/--to-bp/...: include to end of
  // chromosome
  if ((marker_pos_start != -1) && (!markername_to) && (marker_pos_end == -1)) {
    marker_pos_end = 0x7ffffffe;
  }
  if (!chrom_flag_present) {
    init_default_chrom_mask(&chrom_info);
  }
  if (((marker_pos_start != -1) && (!markername_to)) || ((marker_pos_end != -1) && (!markername_from))) {
    // require exactly one chromosome to be defined given --from-bp/--to-bp
    // without --from/--to
    uii = next_set(chrom_info.chrom_mask, 0, CHROM_MASK_INITIAL_WORDS * BITCT);
    if (uii == CHROM_MASK_INITIAL_WORDS * BITCT) {
      uii = 0;
    } else {
      uii = next_set(chrom_info.chrom_mask, uii + 1, CHROM_MASK_INITIAL_WORDS * BITCT);
    }
    if (((uii == CHROM_MASK_INITIAL_WORDS * BITCT) && chrom_info.incl_excl_name_stack) || ((uii != CHROM_MASK_INITIAL_WORDS * BITCT) && (uii || (!chrom_info.incl_excl_name_stack) || chrom_info.incl_excl_name_stack->next))) {
      logerrprint("Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb require a single chromosome to be\nidentified (either explicitly with --chr, or implicitly with --from/--to).\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  }

  if (mperm_save) {
    uii = 0;
    if ((calculation_type & CALC_MODEL) && ((model_modifier & (MODEL_MPERM | MODEL_SET_TEST)) == MODEL_MPERM)) {
      uii++;
    }
    if ((calculation_type & CALC_GLM) && ((glm_modifier & (GLM_MPERM | GLM_SET_TEST)) == GLM_MPERM)) {
      uii++;
    }
    if ((calculation_type & CALC_TESTMISS) && (testmiss_modifier & TESTMISS_MPERM)) {
      uii++;
    }
    if ((calculation_type & CALC_TDT) && ((family_info.tdt_modifier & (TDT_MPERM | TDT_SET_TEST)) == TDT_MPERM)) {
      uii++;
    }
    if ((calculation_type & CALC_DFAM) && ((family_info.dfam_modifier & (DFAM_MPERM | DFAM_SET_TEST)) == DFAM_MPERM)) {
      uii++;
    }

    // no qfam since that's a nonstandard permutation test

    if ((calculation_type & CALC_CMH) && ((cluster.modifier & (CLUSTER_CMH_MPERM | CLUSTER_CMH_SET_TEST)) == CLUSTER_CMH_MPERM)) {
      uii++;
    }
    if (uii != 1) {
      // prevent one permutation test's values from clobbering another's
      logerrprint("Error: --mperm-save{-all} must be used with exactly one max(T) permutation\ntest.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (calculation_type & CALC_MODEL) {
    if (!(model_modifier & (MODEL_ASSOC | MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
      if (mtest_adjust && (!(model_modifier & MODEL_SET_TEST))) {
	// this is actually okay with the set test
	logerrprint("Error: In order to use --model with --adjust, you must include the 'trend',\n'trend-only', 'dom', or 'rec' modifier.\n");
	goto main_ret_INVALID_CMDLINE;
      }
    }
    if (model_cell_ct == -1) {
      model_cell_ct = (model_modifier & MODEL_FISHER)? 0 : 5;
    }
    if ((model_modifier & (MODEL_PERM | MODEL_MPERM | MODEL_GENEDROP)) == MODEL_GENEDROP) {
      model_modifier |= MODEL_PERM;
    }
  }
  if (aperm_present && (calculation_type & (CALC_MODEL | CALC_GLM | CALC_TESTMISS | CALC_TDT | CALC_DFAM | CALC_QFAM | CALC_CMH)) &&
      (!(model_modifier & MODEL_PERM)) &&
      (!(glm_modifier & GLM_PERM)) &&
      (!(testmiss_modifier & TESTMISS_PERM)) &&
      (!(family_info.tdt_modifier & TDT_PERM)) &&
      (!(family_info.dfam_modifier & DFAM_PERM)) &&
      (!(family_info.qfam_modifier & QFAM_PERM)) &&
      (!(cluster.modifier & (CLUSTER_CMH_PERM | CLUSTER_CMH_PERM_BD)))) {
    // If --aperm is present, at least one association analysis command which
    // supports adaptive permutation testing was also specified, but no actual
    // adaptive permutation test is happening, the user is likely to be
    // confused.  Produce a warning.  (Not an error since a sophisticated user
    // may want to use --script with different --aperm defaults.)
    logerrprint("Warning: --aperm only controls the settings for adaptive permutation tests; it\ndoes not cause such a test to be performed.  (Did you forget to add the 'perm'\nmodifier to an association analysis flag?)\n");
  }
  if ((mtest_adjust & (ADJUST_LAMBDA + 1)) == ADJUST_LAMBDA) {
    logerrprint("Error: --lambda must be used with --adjust.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if ((homozyg.modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE)) && (!(calculation_type & CALC_HOMOZYG))) {
    if (homozyg.overlap_min == 0.95) {
      logerrprint("Error: --homozyg-group must be used with another --homozyg... flag.\n");
    } else {
      logerrprint("Error: --homozyg-match must be used with another --homozyg... flag.\n");
    }
    goto main_ret_INVALID_CMDLINE_A;
  }
  if (score_info.data_fname) {
    if (!score_info.range_fname) {
      logerrprint("Error: --q-score-file cannot be used without --q-score-range.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (!(calculation_type & CALC_SCORE)) {
      logerrprint("Error: --q-score-range must be used with --score.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  }
  if (qual_max_thresh != DBL_MAX) {
    if (!qual_filter) {
      logerrprint("Error: --qual-max-threshold must be used with --qual-scores.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (qual_max_thresh < qual_min_thresh) {
      logerrprint("Error: --qual-max-threshold value cannot be negative unless --qual-threshold is\nalso present.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
  }
  if (gene_report_border && (!gene_report_fname)) {
    logerrprint("Error: --gene-list-border must be used with --gene-report.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }

  if ((load_rare == LOAD_RARE_LGEN) && (!mapname[0])) {
    logerrprint("Error: --lgen must be used with --lfile or --map.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  uii = load_params & LOAD_PARAMS_OX_ALL;
  if ((uii == LOAD_PARAMS_OXGEN) || (uii == LOAD_PARAMS_OXBGEN)) {
    logerrprint("Error: --gen/--bgen cannot be used without --data or --sample.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if ((load_rare & LOAD_RARE_TFAM) && (!(load_rare & (LOAD_RARE_TRANSPOSE | LOAD_RARE_TPED)))) {
    logerrprint("Error: --tfam must be used with --tfile or --tped.\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if ((merge_type & MERGE_EQUAL_POS) && (!(calculation_type & CALC_MERGE))) {
    logerrprint("Error: --merge-equal-pos must be used with --merge/--bmerge/--merge-list.\n(Note that you are permitted to merge a fileset with itself.)\n");
    goto main_ret_INVALID_CMDLINE_A;
  }
  if (calculation_type && (!(calculation_type & (~(CALC_FREQ | CALC_MISSING_REPORT)))) && ((geno_thresh != 1.0) || (hwe_thresh != 0.0) || (min_maf != 0.0) || (max_maf != 0.5) || min_ac || (max_ac != 0x7fffffff))) {
    logerrprint("Warning: --freq and --missing complete BEFORE --geno, --hwe, and --maf in\nplink's order of operations.\n");
  }
  // short batch job?
  uii = 0;
  if ((!calculation_type) && (!(load_rare & (LOAD_RARE_LGEN | LOAD_RARE_DUMMY | LOAD_RARE_SIMULATE | LOAD_RARE_TRANSPOSE_MASK | LOAD_RARE_23 | LOAD_RARE_CNV | LOAD_RARE_VCF | LOAD_RARE_BCF)))) {
    if (epi_info.summary_merge_prefix || annot_info.fname || gene_report_fname || (load_rare == LOAD_RARE_DOSAGE) || metaanal_fnames) {
      uii = 1;
    } else if (famname[0] || load_rare) {
      goto main_ret_NULL_CALC;
    }
    // otherwise, autoconversion job
  }
  if (!(load_params || load_rare || uii || (merge_type & MERGE_LIST))) {
    logerrprint("Error: No input dataset.\n");
    goto main_ret_INVALID_CMDLINE_A;
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
    ujj = (uint32_t)time(nullptr);
    sprintf(g_logbuf, "Random number seed: %u\n", ujj);
    logstr(g_logbuf);
    sfmt_init_gen_rand(&g_sfmt, ujj);
  } else {
    if (rseed_ct == 1) {
      sfmt_init_gen_rand(&g_sfmt, rseeds[0]);
    } else {
      sfmt_init_by_array(&g_sfmt, rseeds, rseed_ct);
    }
    free(rseeds);
    rseeds = nullptr;
  }
  // guarantee contiguous malloc space outside of main workspace
  bubble = (char*)malloc(NON_BIGSTACK_MIN * sizeof(char));
  if (!bubble) {
    goto main_ret_NOMEM;
  }

  // see e.g. http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system .
#ifdef __APPLE__
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  llxx = 0;
  sztmp = sizeof(int64_t);
  sysctl(mib, 2, &llxx, &sztmp, nullptr, 0);
  llxx /= 1048576;
#else
#ifdef _WIN32
  memstatus.dwLength = sizeof(memstatus);
  GlobalMemoryStatusEx(&memstatus);
  llxx = memstatus.ullTotalPhys / 1048576;
#else
  llxx = ((uint64_t)sysconf(_SC_PHYS_PAGES)) * ((size_t)sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
  if (!llxx) {
    default_alloc_mb = BIGSTACK_DEFAULT_MB;
  } else if (llxx < (BIGSTACK_MIN_MB * 2)) {
    default_alloc_mb = BIGSTACK_MIN_MB;
  } else {
    default_alloc_mb = llxx / 2;
  }
  if (!malloc_size_mb) {
    malloc_size_mb = default_alloc_mb;
  } else if (malloc_size_mb < BIGSTACK_MIN_MB) {
    malloc_size_mb = BIGSTACK_MIN_MB;
  }
#ifndef __LP64__
  if (malloc_size_mb > 2047) {
    malloc_size_mb = 2047;
  }
#endif
  if (llxx) {
    sprintf(g_logbuf, "%" PRId64 " MB RAM detected; reserving %" PRIdPTR " MB for main workspace.\n", llxx, malloc_size_mb);
  } else {
    sprintf(g_logbuf, "Failed to calculate system memory.  Attempting to reserve %" PRIdPTR " MB.\n", malloc_size_mb);
  }
  logprintb();
  bigstack_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  while (!bigstack_ua) {
    malloc_size_mb = (malloc_size_mb * 3) / 4;
    if (malloc_size_mb < BIGSTACK_MIN_MB) {
      malloc_size_mb = BIGSTACK_MIN_MB;
    }
    bigstack_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
    if (bigstack_ua) {
      LOGPRINTF("Allocated %" PRIdPTR " MB successfully, after larger attempt(s) failed.\n", malloc_size_mb);
    } else if (malloc_size_mb == BIGSTACK_MIN_MB) {
      goto main_ret_NOMEM;
    }
  }
  // force 64-byte align to make cache line sensitivity work
  bigstack_initial_base = (unsigned char*)round_up_pow2((uintptr_t)bigstack_ua, CACHELINE);
  g_bigstack_base = bigstack_initial_base;
  g_bigstack_end = &(bigstack_initial_base[(malloc_size_mb * 1048576 - (uintptr_t)(bigstack_initial_base - bigstack_ua)) & (~(CACHELINE - ONELU))]);
  free(bubble);
  bubble = nullptr;

  // standalone stuff
  if (epi_info.summary_merge_prefix) {
    retval = epi_summary_merge(&epi_info, outname, outname_end);
    if (retval) {
      goto main_ret_1;
    }
  }
  if (annot_info.fname) {
    retval = annotate(&annot_info, (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1, outname, outname_end, pfilter, &chrom_info);
  }
  if (gene_report_fname) {
    retval = gene_report(gene_report_fname, gene_report_glist, gene_report_subset, gene_report_border, (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1, (misc_flags & MISC_EXTRACT_RANGE)? nullptr : extractname, gene_report_snp_field, outname, outname_end, pfilter, &chrom_info);
    if (retval) {
      goto main_ret_1;
    }
  }
  if (metaanal_fnames) {
    // possible todo: make this support --aec (takes a bit of work since
    // chromosome byte in data structure must be widened)
    retval = meta_analysis(metaanal_fnames, metaanal_chrfield_search_order, metaanal_snpfield_search_order, metaanal_bpfield_search_order, metaanal_a1field_search_order, metaanal_a2field_search_order, metaanal_pfield_search_order, metaanal_sefield_search_order, metaanal_essfield_search_order, metaanal_flags, (misc_flags & MISC_EXTRACT_RANGE)? nullptr : extractname, outname, outname_end, output_min_p, &chrom_info);
    if (retval) {
      goto main_ret_1;
    }
  }
  if (load_rare == LOAD_RARE_DOSAGE) {
    if (calculation_type & (~CALC_SCORE)) {
      // with --dosage, there *can't* be anything else to do
      logerrprint("Error: --dosage cannot be used with other PLINK computations.\n");
      goto main_ret_INVALID_CMDLINE;
    }
    pigz_init(g_thread_ct);
    retval = plink1_dosage(&dosage_info, famname, mapname, outname, outname_end, phenoname, extractname, excludename, keepname, removename, keepfamname, removefamname, filtername, makepheno_str, phenoname_str, covar_fname, qual_filter, update_map, update_name, update_ids_fname, update_parents_fname, update_sex_fname, filtervals_flattened, filter_attrib_fname, filter_attrib_liststr, filter_attrib_sample_fname, filter_attrib_sample_liststr, qual_min_thresh, qual_max_thresh, thin_keep_prob, thin_keep_ct, min_bp_space, mfilter_col, fam_cols, missing_pheno, output_missing_pheno, mpheno_col, pheno_modifier, &chrom_info, tail_bottom, tail_top, misc_flags, filter_flags, sex_missing_pheno, update_sex_col, &cluster, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, &snps_range_list, covar_modifier, &covar_range_list, mwithin_col, glm_modifier, glm_vif_thresh, output_min_p, known_procs, &score_info);
    // unconditional; note that plink1_dosage() currently doesn't even bother
    // to pop stuff off the stack when it's done
    goto main_ret_1;
  }
  // quit if there's nothing else to do
  if (uii) {
    goto main_ret_1;
  }

  pigz_init(g_thread_ct);
  if (load_rare & (LOAD_RARE_GRM | LOAD_RARE_GRM_BIN)) {
    // --unrelated-heritability and --rel-cutoff batch mode special cases
#ifndef NOLAPACK
    if (calculation_type & CALC_UNRELATED_HERITABILITY) {
      retval = unrelated_herit_batch(load_rare & LOAD_RARE_GRM_BIN, pedname, phenoname, mpheno_col, phenoname_str, missing_pheno, &rel_info);
    } else {
#endif
      retval = rel_cutoff_batch(load_rare & LOAD_RARE_GRM_BIN, pedname, outname, outname_end, &rel_info);
#ifndef NOLAPACK
    }
#endif
  } else if (load_rare & LOAD_RARE_CNV) {
    retval = plink_cnv(outname, outname_end, pedname, mapname, famname, phenoname, keepname, removename, filtername, misc_flags, update_chr, update_cm, update_map, update_name, update_ids_fname, update_parents_fname, update_sex_fname, filtervals_flattened, filter_flags, cnv_calc_type, cnv_min_seglen, cnv_max_seglen, cnv_min_score, cnv_max_score, cnv_min_sites, cnv_max_sites, cnv_intersect_filter_type, cnv_intersect_filter_fname, cnv_subset_fname, cnv_overlap_type, cnv_overlap_val, cnv_freq_type, cnv_freq_val, cnv_freq_val2, cnv_test_window, segment_modifier, segment_spanning_fname, cnv_sample_mperms, cnv_test_mperms, cnv_test_region_mperms, cnv_enrichment_test_mperms, marker_pos_start, marker_pos_end, &chrom_info);
  } else if (load_rare & LOAD_RARE_GVAR) {
    retval = plink_gvar(outname, outname_end, pedname, mapname, famname);
  } else {
    if (filter_flags) {
      if (!calculation_type) {
	logerrprint("Error: Basic file conversions do not support regular filtering operations.\nRerun your command with --make-bed.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (calculation_type == CALC_MERGE) {
	logerrprint("Error: Basic merge does not support regular filtering operations.  Rerun your\ncommand with --make-bed.\n");
	goto main_ret_INVALID_CMDLINE;
      }
    }
    if (load_rare || (load_params & (LOAD_PARAMS_TEXT_ALL | LOAD_PARAMS_OX_ALL))) {
      sptr = outname_end;
      if (calculation_type && (!(misc_flags & MISC_KEEP_AUTOCONV))) {
        sptr = memcpyb(sptr, "-temporary", 11);
      }
      uii = (sptr - outname);
      if (load_rare == LOAD_RARE_LGEN) {
        retval = lgen_to_bed(pedname, mapname, famname, outname, sptr, missing_pheno, misc_flags, lgen_modifier, lgen_reference_fname, &chrom_info);
      } else if (load_rare & LOAD_RARE_TRANSPOSE_MASK) {
        retval = transposed_to_bed(pedname, famname, outname, sptr, misc_flags, &chrom_info);
      } else if (load_rare & LOAD_RARE_VCF) {
	retval = vcf_to_bed(pedname, outname, sptr, missing_pheno, misc_flags, const_fid, id_delim, vcf_idspace_to, vcf_min_qual, vcf_filter_exceptions_flattened, vcf_min_gq, vcf_min_gp, (uint32_t)vcf_half_call, &chrom_info);
      } else if (load_rare & LOAD_RARE_BCF) {
	retval = bcf_to_bed(pedname, outname, sptr, missing_pheno, misc_flags, const_fid, id_delim, vcf_idspace_to, vcf_min_qual, vcf_filter_exceptions_flattened, &chrom_info);
      } else if (load_rare == LOAD_RARE_23) {
        retval = bed_from_23(pedname, outname, sptr, modifier_23, fid_23, iid_23, (pheno_23 == DBL_MAX)? ((double)missing_pheno) : pheno_23, misc_flags, paternal_id_23, maternal_id_23, &chrom_info);
      } else if (load_rare & LOAD_RARE_DUMMY) {
	retval = generate_dummy(outname, sptr, dummy_flags, dummy_marker_ct, dummy_sample_ct, dummy_missing_geno, dummy_missing_pheno, missing_pheno);
      } else if (load_rare & LOAD_RARE_SIMULATE) {
	// no need to support zero samples/variants here since --dummy takes
	// care of generating those test cases
	retval = simulate_dataset(outname, sptr, simulate_flags, simulate_fname, simulate_cases, simulate_controls, simulate_prevalence, simulate_qt_samples, simulate_missing, simulate_label);
	free(simulate_fname);
	simulate_fname = nullptr;
	if (simulate_label) {
	  free(simulate_label);
	  simulate_label = nullptr;
	}
      } else if (load_params & LOAD_PARAMS_OX_ALL) {
	retval = oxford_to_bed(pedname, mapname, outname, sptr, oxford_single_chr, oxford_pheno_name, hard_call_threshold, oxford_missing_code, missing_pheno, misc_flags, (load_params / LOAD_PARAMS_OXBGEN) & 1, &chrom_info);
      } else {
        retval = ped_to_bed(pedname, mapname, outname, sptr, fam_cols, misc_flags, missing_pheno, &chrom_info);
	fam_cols |= FAM_COL_1 | FAM_COL_34 | FAM_COL_5;
	if (!(fam_cols & FAM_COL_6)) {
          fam_cols |= FAM_COL_6;
	  missing_pheno = -9;
	}
      }
      if (retval || (!calculation_type)) {
	goto main_ret_1;
      }
      memcpy(memcpya(pedname, outname, uii), ".bed", 5);
      memcpy(memcpya(mapname, outname, uii), ".bim", 5);
      memcpy(memcpya(famname, outname, uii), ".fam", 5);
      if (calculation_type && (!(misc_flags & MISC_KEEP_AUTOCONV))) {
	if (push_ll_str(pedname, &file_delete_list) || push_ll_str(mapname, &file_delete_list) || push_ll_str(famname, &file_delete_list)) {
	  goto main_ret_NOMEM;
	}
      }
      *outname_end = '\0';
    }
    if (rel_info.ibc_type == 3) { // todo: make this less ugly
      rel_info.ibc_type = 0;
    } else if (!rel_info.ibc_type) {
      rel_info.ibc_type = 1;
    }
    retval = plink(outname, outname_end, pedname, mapname, famname, cm_map_fname, cm_map_chrname, phenoname, extractname, excludename, keepname, removename, keepfamname, removefamname, filtername, freqname, distance_wts_fname, read_dists_fname, read_dists_id_fname, evecname, mergename1, mergename2, mergename3, missing_mid_template, missing_marker_id_match, makepheno_str, phenoname_str, a1alleles, a2alleles, recode_allele_name, covar_fname, update_alleles_fname, read_genome_fname, qual_filter, update_chr, update_cm, update_map, update_name, update_ids_fname, update_parents_fname, update_sex_fname, loop_assoc_fname, flip_fname, flip_subset_fname, sample_sort_fname, filtervals_flattened, condition_mname, condition_fname, filter_attrib_fname, filter_attrib_liststr, filter_attrib_sample_fname, filter_attrib_sample_liststr, rplugin_fname, rplugin_host_or_socket, rplugin_port, qual_min_thresh, qual_max_thresh, thin_keep_prob, thin_keep_sample_prob, new_id_max_allele_len, thin_keep_ct, thin_keep_sample_ct, min_bp_space, mfilter_col, fam_cols, missing_pheno, output_missing_pheno, mpheno_col, pheno_modifier, &chrom_info, &oblig_missing_info, &family_info, check_sex_fthresh, check_sex_mthresh, check_sex_f_yobs, check_sex_m_yobs, distance_exp, min_maf, max_maf, geno_thresh, mind_thresh, hwe_thresh, tail_bottom, tail_top, misc_flags, filter_flags, calculation_type, dist_calc_type, groupdist_iters, groupdist_d, regress_iters, regress_d, parallel_idx, parallel_tot, splitx_bound1, splitx_bound2, ppc_gap, sex_missing_pheno, update_sex_col, hwe_modifier, min_ac, max_ac, genome_modifier, genome_min_pi_hat, genome_max_pi_hat, &homozyg, &cluster, neighbor_n1, neighbor_n2, &set_info, &ld_info, &epi_info, &clump_info, &rel_info, &score_info, recode_modifier, allelexxxx, merge_type, sample_sort, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, &snps_range_list, write_var_range_ct, covar_modifier, &covar_range_list, write_covar_modifier, write_covar_dummy_max_categories, dupvar_modifier, mwithin_col, model_modifier, (uint32_t)model_cell_ct, model_mperm_val, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, &parameters_range_list, &tests_range_list, ci_size, pfilter, output_min_p, mtest_adjust, adjust_lambda, gxe_mcovar, &aperm, mperm_save, ibs_test_perms, perm_batch_size, lasso_h2, lasso_lambda_iters, lasso_minlambda, &lasso_select_covars_range_list, testmiss_modifier, testmiss_mperm_val, permphe_ct, known_procs, &file_delete_list);
  }
  while (0) {
  main_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  main_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  main_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  main_ret_INVALID_CMDLINE_UNRECOGNIZED:
    invalid_arg(argv[cur_arg]);
    logerrprintb();
    logerrprint(errstr_append);
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_INVALID_CMDLINE_INPUT_CONFLICT:
    LOGERRPRINTF("Error: --%s conflicts with another input flag.\n%s", argptr, errstr_append);
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_INVALID_CMDLINE_WWA:
    wordwrapb(0);
  main_ret_INVALID_CMDLINE_2A:
    logerrprintb();
  main_ret_INVALID_CMDLINE_A:
    logerrprint(errstr_append);
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_INVALID_CMDLINE_WW:
    wordwrapb(0);
  main_ret_INVALID_CMDLINE_2:
    logerrprintb();
  main_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_NULL_CALC:
    if (filter_flags) {
      logerrprint("Warning: No output requested.  (Did you forget --make-bed?)  Exiting.\n");
    } else {
      logerrprint("Warning: No output requested.  Exiting.\n");
    }
    fputs(g_cmdline_format_str, stdout);
    fputs(notestr_null_calc2, stdout);
    retval = RET_NULL_CALC;
#ifdef STABLE_BUILD
    break;
  main_unstable_disabled:
    // see the UNSTABLE macro in plink_common.h
    memcpy(g_logbuf, "Error: --", 9);
    strcpy(sptr, " is either unfinished or not yet well-tested. If you wish to help with testing, use the latest development build.\n");
    wordwrapb(0);
    logerrprintb();
    retval = RET_CALC_NOT_YET_SUPPORTED;
#endif
  }
 main_ret_1:
  disp_exit_msg(retval);
  while (0) {
  main_ret_NOMEM_NOLOG:
    print_ver();
  main_ret_NOMEM_NOLOG2:
    fputs(errstr_nomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
    }
    retval = RET_NOMEM;
    break;
  main_ret_READ_FAIL_NOLOG:
    print_ver();
    fputs(errstr_read, stderr);
    retval = RET_READ_FAIL;
    break;
  }
 main_ret_NOLOG:
  fclose_cond(scriptfile);
  free_cond(bubble);
  free_cond(bigstack_ua);
  free_cond(subst_argv);
  free_cond(script_buf);
  free_cond(rerun_buf);
  free_cond(flag_buf);
  free_cond(flag_map);
  free_cond(makepheno_str);
  free_cond(phenoname_str);
  free_cond(a1alleles);
  free_cond(a2alleles);
  free_cond(sample_sort_fname);
  free_cond(filtervals_flattened);
  free_cond(evecname);
  free_cond(filtername);
  free_cond(distance_wts_fname);
  free_cond(read_dists_fname);
  free_cond(read_dists_id_fname);
  free_cond(freqname);
  free_cond(extractname);
  free_cond(excludename);
  free_cond(keepname);
  free_cond(removename);
  free_cond(keepfamname);
  free_cond(removefamname);
  free_cond(cm_map_fname);
  free_cond(cm_map_chrname);
  free_cond(phenoname);
  free_cond(recode_allele_name);
  free_cond(markername_from);
  free_cond(markername_to);
  free_cond(markername_snp);
  free_range_list(&snps_range_list);
  free_range_list(&covar_range_list);
  free_range_list(&lasso_select_covars_range_list);
  free_range_list(&parameters_range_list);
  free_range_list(&tests_range_list);
  free_cond(lgen_reference_fname);
  free_cond(covar_fname);
  free_cond(update_alleles_fname);
  free_cond(qual_filter);
  free_cond(update_chr);
  free_cond(update_cm);
  free_cond(update_map);
  free_cond(update_name);
  free_cond(oxford_missing_code);
  free_cond(oxford_single_chr);
  free_cond(oxford_pheno_name);
  free_cond(update_ids_fname);
  free_cond(update_parents_fname);
  free_cond(update_sex_fname);
  free_cond(loop_assoc_fname);
  free_cond(flip_fname);
  free_cond(flip_subset_fname);
  free_cond(read_genome_fname);
  free_cond(rseeds);
  free_cond(simulate_fname);
  free_cond(simulate_label);
  free_cond(cnv_intersect_filter_fname);
  free_cond(cnv_subset_fname);
  free_cond(segment_spanning_fname);
  free_cond(fid_23);
  free_cond(iid_23);
  free_cond(paternal_id_23);
  free_cond(maternal_id_23);
  free_cond(condition_mname);
  free_cond(condition_fname);
  free_cond(missing_mid_template);
  free_cond(missing_marker_id_match);
  free_cond(filter_attrib_fname);
  free_cond(filter_attrib_liststr);
  free_cond(filter_attrib_sample_fname);
  free_cond(filter_attrib_sample_liststr);
  free_cond(const_fid);
  free_cond(vcf_filter_exceptions_flattened);
  free_cond(gene_report_fname);
  free_cond(gene_report_glist);
  free_cond(gene_report_subset);
  free_cond(gene_report_snp_field);
  free_cond(metaanal_fnames);
  free_cond(metaanal_chrfield_search_order);
  free_cond(metaanal_snpfield_search_order);
  free_cond(metaanal_bpfield_search_order);
  free_cond(metaanal_a1field_search_order);
  free_cond(metaanal_a2field_search_order);
  free_cond(metaanal_pfield_search_order);
  free_cond(metaanal_sefield_search_order);
  free_cond(metaanal_essfield_search_order);
  free_cond(rplugin_fname);
  free_cond(rplugin_host_or_socket);

  oblig_missing_cleanup(&oblig_missing_info);
  cluster_cleanup(&cluster);
  set_cleanup(&set_info, &annot_info);
  ld_epi_cleanup(&ld_info, &epi_info, &clump_info);
  rel_cleanup(&rel_info);
  misc_cleanup(&score_info);
  dosage_cleanup(&dosage_info);
  if (file_delete_list) {
    do {
      ll_str_ptr = file_delete_list->next;
      unlink(file_delete_list->ss);
      free(file_delete_list);
      file_delete_list = ll_str_ptr;
    } while (file_delete_list);
  }
  cleanup_chrom_info(&chrom_info);
  if (g_logfile) {
    if (!g_log_failed) {
      logstr("\nEnd time: ");
      time(&rawtime);
      logstr(ctime(&rawtime));
      if (fclose(g_logfile)) {
	fflush(stdout);
	fputs("Error: Failed to finish writing to log.\n", stderr);
      }
    } else {
      fclose(g_logfile);
    }
    g_logfile = nullptr;
  }
  if (misc_flags & MISC_GPLINK) {
    memcpy(outname_end, ".gplink", 8);
    g_logfile = fopen(outname, "w");
    if (g_logfile) { // can't do much if an error occurs here...
      putc_unlocked(retval? '1' : '0', g_logfile);
      putc_unlocked('\n', g_logfile);
      fclose(g_logfile);
    }
  }

  return retval;
}
