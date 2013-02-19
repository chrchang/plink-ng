// WDIST weighted genomic distance calculator
// Copyright (C) 2013  Christopher Chang  chrchang@alumni.caltech.edu

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Uncomment "#define NOLAPACK" in wdist_common.c to build without LAPACK.

#include <ctype.h>
#include <time.h>
#include <unistd.h>
#if _WIN32
// needed for MEMORYSTATUSEX
#ifndef _WIN64
#define WINVER 0x0500
#endif
#endif
#include "wdist_common.h"

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include "wdist_assoc.h"
#include "wdist_calc.h"
#include "wdist_data.h"
#include "wdist_dosage.h"
#include "wdist_stats.h"

// default jackknife iterations
#define ITERS_DEFAULT 100000
#define DEFAULT_PPC_GAP 500000
#define MAX_PCS_DEFAULT 20

#define LOAD_RARE_GRM 1
#define LOAD_RARE_LGEN 2
#define LOAD_RARE_TRANSPOSE 4
#define LOAD_RARE_TPED 8
#define LOAD_RARE_TFAM 16
#define LOAD_RARE_TRANSPOSE_MASK (LOAD_RARE_TRANSPOSE | LOAD_RARE_TPED | LOAD_RARE_TFAM)
#define LOAD_RARE_DUMMY 32

// maximum number of usable cluster computers, this is arbitrary though it
// shouldn't be larger than 2^31 - 1
#define PARALLEL_MAX 32768

const char ver_str[] =
  "WDIST v0.17.1"
#ifdef NOLAPACK
  "NL"
#endif
#ifdef __LP64__
  " 64-bit"
#else
  " 32-bit"
#endif
  " (19 Feb 2013)";
const char ver_str2[] =
  "    https://www.cog-genomics.org/wdist\n"
  "(C) 2013 Christopher Chang, GNU General Public License version 3\n";
const char errstr_ped_format[] = "Error: Improperly formatted .ped file.\n";
const char errstr_phenotype_format[] = "Error: Improperly formatted phenotype file.\n";
const char errstr_filter_format[] = "Error: Improperly formatted filter file.\n";
const char errstr_freq_format[] = "Error: Improperly formatted frequency file.\n";
const char cmdline_format_str[] = "\n  wdist [input flag(s)...] [command flag(s)...] {other flag(s)...}\n  wdist --help {flag name(s)...}\n\n";
const char notestr_null_calc[] = "Note: No output requested.  Exiting.\n";
const char notestr_null_calc2[] = "Commands include --freqx, --hardy, --ibc, --distance, --genome, --model, --gxe,\n--make-rel, --make-grm, --rel-cutoff, --regress-distance,\n--regress-pcs-distance, --make-bed, --recode, --merge-list, and\n--write-snplist.\n\n'wdist --help | more' describes all functions (warning: long).\n";

int32_t edit1_match(int32_t len1, char* s1, int32_t len2, char* s2) {
  // permit one difference of the following forms:
  // - inserted/deleted character
  // - replaced character
  // - adjacent pair of swapped characters
  int32_t diff_found = 0;
  int32_t pos = 0;
  if (len1 == len2) {
    while (pos < len1) {
      if (s1[pos] != s2[pos]) {
	if (diff_found) {
	  if ((diff_found == 2) || (s1[pos] != s2[pos - 1]) || (s1[pos - 1] != s2[pos])) {
	    return 0;
	  }
	}
	diff_found++;
      }
      pos++;
    }
  } else if (len1 == len2 - 1) {
    do {
      if (s1[pos - diff_found] != s2[pos]) {
	if (diff_found) {
	  return 0;
	}
	diff_found++;
      }
      pos++;
    } while (pos < len2);
  } else if (len1 == len2 + 1) {
    do {
      if (s1[pos] != s2[pos - diff_found]) {
	if (diff_found) {
	  return 0;
	}
	diff_found++;
      }
      pos++;
    } while (pos < len1);
  } else {
    return 0;
  }
  return 1;
}

#define MAX_EQUAL_HELP_PARAMS 15

typedef struct {
  int32_t iters_left;
  uint32_t param_ct;
  char** argv;
  uintptr_t unmatched_ct;
  uintptr_t* all_match_arr;
  uintptr_t* prefix_match_arr;
  uintptr_t* perfect_match_arr;
  uint32_t* param_lens;
  int32_t preprint_newline;
} Help_ctrl;

void help_print(const char* cur_params, Help_ctrl* help_ctrl_ptr, int32_t postprint_newline, const char* payload) {
  // unmatched_ct fixed during call, *unmatched_ct_ptr may decrease
  uint32_t unmatched_ct = help_ctrl_ptr->unmatched_ct;
  int32_t print_this = 0;
  int32_t cur_param_lens[MAX_EQUAL_HELP_PARAMS];
  char* cur_param_start[MAX_EQUAL_HELP_PARAMS];
  int32_t cur_param_ct;
  int32_t cur_param_idx;
  uint32_t arg_uidx;
  uint32_t arg_idx;
  int32_t uii;
  int32_t payload_len;
  char* payload_ptr;
  char* line_end;
  char* payload_end;
  if (help_ctrl_ptr->param_ct) {
    strcpy(tbuf, cur_params);
    cur_param_ct = 1;
    cur_param_start[0] = tbuf;
    payload_ptr = strchr(tbuf, '\t');
    while (payload_ptr) {
      *payload_ptr++ = '\0';
      cur_param_start[cur_param_ct++] = payload_ptr;
      payload_ptr = strchr(payload_ptr, '\t');
    }
    if (help_ctrl_ptr->iters_left) {
      if (help_ctrl_ptr->unmatched_ct) {
	arg_uidx = 0;
	if (help_ctrl_ptr->iters_left == 2) {
	  for (arg_idx = 0; arg_idx < unmatched_ct; arg_idx++) {
	    arg_uidx = next_non_set_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
	    for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	      if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
		set_bit_noct(help_ctrl_ptr->perfect_match_arr, arg_uidx);
		set_bit_noct(help_ctrl_ptr->prefix_match_arr, arg_uidx);
		set_bit_sub(help_ctrl_ptr->all_match_arr, arg_uidx, &(help_ctrl_ptr->unmatched_ct));
		break;
	      }
	    }
	    arg_uidx++;
	  }
	} else {
	  for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	    cur_param_lens[cur_param_idx] = strlen(cur_param_start[cur_param_idx]);
	  }
	  for (arg_idx = 0; arg_idx < unmatched_ct; arg_idx++) {
	    arg_uidx = next_non_set_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
	    uii = help_ctrl_ptr->param_lens[arg_uidx];
	    for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	      if (cur_param_lens[cur_param_idx] > uii) {
		if (!memcmp(help_ctrl_ptr->argv[arg_uidx], cur_param_start[cur_param_idx], uii)) {
		  set_bit_noct(help_ctrl_ptr->prefix_match_arr, arg_uidx);
		  set_bit_sub(help_ctrl_ptr->all_match_arr, arg_uidx, &(help_ctrl_ptr->unmatched_ct));
		  break;
		}
	      }
	    }
	    arg_uidx++;
	  }
	}
      }
    } else {
      for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	cur_param_lens[cur_param_idx] = strlen(cur_param_start[cur_param_idx]);
      }
      for (arg_uidx = 0; arg_uidx < help_ctrl_ptr->param_ct; arg_uidx++) {
	if (is_set(help_ctrl_ptr->prefix_match_arr, arg_uidx)) {
	  if (!print_this) {
	    if (is_set(help_ctrl_ptr->perfect_match_arr, arg_uidx)) {
	      for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
		if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
		  print_this = 1;
		  break;
		}
	      }
	    } else {
	      uii = help_ctrl_ptr->param_lens[arg_uidx];
	      for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
		if (cur_param_lens[cur_param_idx] > uii) {
		  if (!memcmp(help_ctrl_ptr->argv[arg_uidx], cur_param_start[cur_param_idx], uii)) {
		    print_this = 1;
		    break;
		  }
		}
	      }
	    }
	  }
	} else {
	  for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	    if (edit1_match(cur_param_lens[cur_param_idx], cur_param_start[cur_param_idx], help_ctrl_ptr->param_lens[arg_uidx], help_ctrl_ptr->argv[arg_uidx])) {
	      print_this = 1;
	      set_bit_sub(help_ctrl_ptr->all_match_arr, arg_uidx, &(help_ctrl_ptr->unmatched_ct));
	      break;
	    }
	  }
	}
      }
      if (print_this) {
	payload_len = strlen(payload);
	if (payload[payload_len - 2] == '\n') {
	  payload_end = (char*)(&(payload[payload_len - 1]));
	} else {
	  payload_end = (char*)(&(payload[payload_len]));
	}
	if (help_ctrl_ptr->preprint_newline) {
	  putchar('\n');
	}
	help_ctrl_ptr->preprint_newline = postprint_newline;
	payload_ptr = (char*)payload;
	do {
	  line_end = strchr(payload_ptr, '\n') + 1;
	  uii = (uint32_t)(line_end - payload_ptr);
	  if (uii > 2) {
	    payload_ptr = &(payload_ptr[2]);
	    uii -= 2;
	  }
	  memcpy(tbuf, payload_ptr, uii);
	  tbuf[uii] = '\0';
	  fputs(tbuf, stdout);
	  payload_ptr = line_end;
	} while (payload_ptr < payload_end);
      }
    }
  } else {
    fputs(payload, stdout);
  }
}

int32_t disp_help(uint32_t param_ct, char** argv) {
  // yes, this is overkill.  But it should be a good template for other
  // command-line programs to use.
  uint32_t param_ctl = (param_ct + (BITCT - 1)) / BITCT;
  int32_t retval = 0;
  Help_ctrl help_ctrl;
  uint32_t arg_uidx;
  uint32_t arg_idx;
  uint32_t net_unmatched_ct;
  int32_t col_num;
  int32_t leading_dashes;
  help_ctrl.iters_left = param_ct? 2 : 0;
  help_ctrl.param_ct = param_ct;
  help_ctrl.argv = argv;
  help_ctrl.unmatched_ct = param_ct;
  help_ctrl.param_lens = NULL;
  help_ctrl.all_match_arr = NULL;
  help_ctrl.argv = NULL;
  if (param_ct) {
    help_ctrl.param_lens = (uint32_t*)malloc(param_ct * sizeof(int32_t));
    if (!help_ctrl.param_lens) {
      goto disp_help_ret_NOMEM;
    }
    help_ctrl.all_match_arr = (uintptr_t*)malloc(param_ctl * 3 * sizeof(intptr_t));
    if (!help_ctrl.all_match_arr) {
      goto disp_help_ret_NOMEM;
    }
    leading_dashes = 0;
    for (arg_uidx = 0; arg_uidx < param_ct; arg_uidx++) {
      if (argv[arg_uidx][0] == '-') {
	leading_dashes = 1;
	break;
      }
    }
    if (leading_dashes) {
      help_ctrl.argv = (char**)malloc(param_ct * sizeof(char*));
      if (!help_ctrl.argv) {
	goto disp_help_ret_NOMEM;
      }
      for (arg_uidx = 0; arg_uidx < param_ct; arg_uidx++) {
	if (argv[arg_uidx][0] == '-') {
	  if (argv[arg_uidx][1] == '-') {
	    help_ctrl.argv[arg_uidx] = &(argv[arg_uidx][2]);
	  } else {
	    help_ctrl.argv[arg_uidx] = &(argv[arg_uidx][1]);
	  }
	} else {
	  help_ctrl.argv[arg_uidx] = argv[arg_uidx];
	}
      }
    } else {
      help_ctrl.argv = argv;
    }
    for (arg_idx = 0; arg_idx < param_ct; arg_idx++) {
      help_ctrl.param_lens[arg_idx] = strlen(help_ctrl.argv[arg_idx]);
    }
    fill_ulong_zero(help_ctrl.all_match_arr, param_ctl * 3);
    help_ctrl.prefix_match_arr = &(help_ctrl.all_match_arr[param_ctl]);
    help_ctrl.perfect_match_arr = &(help_ctrl.all_match_arr[param_ctl * 2]);
    help_ctrl.preprint_newline = 1;
  } else {
    help_ctrl.argv = NULL;
    fputs(
"\nIn the command line flag definitions that follow,\n"
"  * [square brackets] denote a required parameter, where the text between the\n"
"    brackets describes its nature.\n"
"  * <angle brackets> denote an optional modifier (or if '|' is present, a set\n"
"    of mutually exclusive optional modifiers).  Use the EXACT text in the\n"
"    definition, e.g. '--distance square0'.\n"
"  * There's one exception to the angle brackets/exact text rule: when an angle\n"
"    bracket term ends with '=[value]', '[value]' designates a variable\n"
"    parameter.\n"
"  * {curly braces} denote an optional parameter, where the text between the\n"
"    braces describes its nature.\n"
"  * An ellipsis (...) indicates that you may enter multiple parameters of the\n"
"    specified type.\n"
, stdout);
    fputs(cmdline_format_str, stdout);
    fputs(
"Each WDIST run requires exactly one main input fileset.  The following flags\n"
"are available for defining its form and location:\n\n"
, stdout);
  }
  do {
    help_print("bfile\tbed\tbim\tfam", &help_ctrl, 1,
"  --bfile {prefix} : Specify .bed/.bim/.fam prefix (default 'wdist').\n"
"  --bed [filename] : Specify full name of .bed file.\n"
"  --bim [filename] : Specify full name of .bim file.\n"
"  --fam [filename] : Specify full name of .fam file.\n\n"
	       );
    help_print("file\tped\tmap", &help_ctrl, 1,
"  --file {prefix}  : Specify prefix for .ped and .map files (default 'wdist').\n"
"  --ped [filename] : Specify full name of .ped file.\n"
"  --map [filename] : Specify full name of .map file.\n\n"
	       );
    help_print("tfile\ttped\ttfam", &help_ctrl, 1,
"  --tfile {prefix} : Specify .tped/.tfam prefix (default 'wdist').\n"
"  --tped [fname]   : Specify full name of .tped file.\n"
"  --tfam [fname]   : Specify full name of .tfam file.\n\n"
	       );
    help_print("lfile", &help_ctrl, 1,
"  --lfile {prefix} : Specify .lgen/.map/.fam (long-format fileset) prefix.\n\n"
	       );
    help_print("data\tgen\tsample", &help_ctrl, 1,
"  --data {prefix}  : Specify Oxford .gen/.sample prefix (default 'wdist').\n"
"  --gen [filename] : Specify full name of .gen file.\n"
"  --sample [fname] : Specify full name of .sample file.\n\n"
	       );
    help_print("grm\trel-cutoff\tgrm-cutoff", &help_ctrl, 1,
"  --grm {prefix}   : Load a GCTA relationship matrix (.grm.gz + .grm.id) for\n"
"                     --rel-cutoff.\n\n"
	       );
    help_print("dummy", &help_ctrl, 1,
"  --dummy [indiv ct] [marker ct] {missing geno freq} {missing pheno freq}\n"
"          <acgt | 1234 | 12> <scalar-pheno>\n"
"    This generates a fake input dataset with the specified number of\n"
"    individuals and markers.  By default, the missing genotype and phenotype\n"
"    frequencies are zero, and genotypes are As and Bs (change the latter with\n"
"    'acgt'/'1234'/'12').  The 'scalar-pheno' modifier causes a normally\n"
"    distributed scalar phenotype to be generated instead of a binary one.\n\n"
	       );
    if (!param_ct) {
      fputs(
"Output files have names of the form 'wdist.{extension}' by default.  You can\n"
"change the 'wdist' prefix with\n\n"
, stdout);
    }
    help_print("out", &help_ctrl, 1,
"  --out [prefix]   : Specify prefix for output files.\n\n"
	       );
    if (!param_ct) {
      fputs(
"WDIST automatically converts PLINK text filesets to binary during the loading\n"
"process (the new fileset is saved to {output prefix}.bed + .bim + .fam, unless\n"
"that would conflict with the output of another command like --make-bed).  You\n"
"are encouraged to directly use the new binary fileset in future runs.\n\n"
"Every run also requires at least one of the following commands (unless you just\n"
"want automatic text-to-binary conversion):\n\n"
, stdout);
    }

    help_print("freq\tfreqx\tfrqx\tcounts", &help_ctrl, 1,
"  --freq <counts>\n"
"  --freqx\n"
"    --freq generates an allele frequency report identical to that of PLINK\n"
"    --freq (or --freq --counts, if the 'counts' modifier is used).  Using the\n"
"    --freqx flag instead causes the MAF and NCHROBS columns to be replaced with\n"
"    homozygote and heterozygote counts, which (when reloaded with --read-freq)\n"
"    allow distance matrix terms to be weighted consistently through multiple\n"
"    filtering runs.\n\n"
		);
    help_print("hardy", &help_ctrl, 1,
"  --hardy\n"
"    Generates a Hardy-Weinberg exact test p-value report.  (This does NOT\n"
"    filter on the p-value; use --hwe for that.)\n\n"
	       );
    help_print("ibc\thet", &help_ctrl, 1,
"  --ibc\n"
"    Calculates inbreeding coefficients in three different ways.\n"
"    * For more details, see Yang J, Lee SH, Goddard ME and Visscher PM.  GCTA:\n"
"      a tool for Genome-wide Complex Trait Analysis.  Am J Hum Genet. 2011 Jan\n"
"      88(1): 76-82.  This paper also describes the relationship matrix\n"
"      computation we implement.\n\n"
	       );
    help_print("distance\tregress-pcs-distance", &help_ctrl, 1,
"  --distance <square | square0 | triangle> <gz | bin> <ibs> <1-ibs> <alct> <3d>\n"
"             <flat-missing>\n"
"    Writes a lower-triangular tab-delimited table of (weighted) genomic\n"
"    distances in allele count units to {output prefix}.dist, and a list of the\n"
"    corresponding family/individual IDs to {output prefix}.dist.id.  The first\n"
"    row of the .dist file contains a single {genotype 1-genotype 2} distance,\n"
"    the second row has the {genotype 1-genotype 3} and {genotype 2-genotype 3}\n"
"    distances in that order, etc.\n"
"    * If the 'square' or 'square0' modifier is present, a square matrix is\n"
"      written instead; 'square0' fills the upper right triangle with zeroes.\n"
"    * If the 'gz' modifier is present, a compressed .dist.gz file is written\n"
"      instead of a plain text file.\n"
"    * If the 'bin' modifier is present, a binary (square) matrix of\n"
"      double-precision floating point values, suitable for loading from R, is\n"
"      instead written to {output prefix}.dist.bin.  This can be combined with\n"
"      'square0' if you still want the upper right zeroed out, or 'triangle' if\n"
"      you don't want to pad the upper right at all.\n"
"    * If the 'ibs' modifier is present, an identity-by-state matrix is written\n"
"      to {output prefix}.mibs.  '1-ibs' causes distances expressed as genomic\n"
"      proportions (i.e. 1 - IBS) to be written to {output prefix}.mdist.\n"
"      Combine with 'alct' if you want to generate the usual .dist file as well.\n"
"    * With dosage data, the '3d' modifier considers 0-1-2 allele count\n"
"      probabilities separately, instead of collapsing them into an expected\n"
"      value and a missingness probability.\n"
"    * By default, distance rescaling in the presence of missing markers is\n"
"      sensitive to allele count distributions: if allele A contributes, on\n"
"      average, twice as much to other pairwise distances as allele B, a missing\n"
"      allele A will result in twice as large of a missingness correction.  To\n"
"      turn this off, use the 'flat-missing' modifier.\n"
	       );
    help_print("matrix\tdistance-matrix", &help_ctrl, 1,
"  --matrix\n"
"  --distance-matrix\n"
"    These generate space-delimited text matrices, and are included for\n"
"    backwards compatibility with old scripts relying on the corresponding PLINK\n"
"    flags.  New scripts should migrate to '--distance ibs' and '--distance\n"
"    1-ibs', which support output formats better suited to parallel computation\n"
"    and have more accurate handling of missing markers.\n\n"
		);
    help_print("genome\tZ-genome", &help_ctrl, 1,
"  --genome <gz> <full> <unbounded>\n"
"    Identity-by-descent analysis.  This yields the same output as PLINK\n"
"    --genome/--Z-genome, and the 'full' and 'unbounded' modifiers have the same\n"
"    effect as PLINK's --genome-full and --unbounded flags.\n\n"
		);
    help_print("assoc\tmodel\tfisher\tperm\tmperm\tperm-count\tcounts\tp2\tmodel-dom\tmodel-gen\tmodel-rec\tmodel-trend\tgenedrop\tqt-means\ttrend", &help_ctrl, 1,
"  --assoc <perm | mperm=[value]> <genedrop> <perm-count> <fisher> <counts> <p2>\n"
"  --assoc <perm | mperm=[value]> <perm-count> <qt-means>\n"
"  --model <perm | mperm=[value]> <genedrop> <perm-count> <fisher | trend-only>\n"
"          <dom | rec | gen | trend>\n"
"    Basic association analysis report.\n"
"    Given a case/control phenotype, --assoc performs a 1df chi-square allelic\n"
"    test, while --model performs 4 other tests as well (1df dominant gene\n"
"    action, 1df recessive gene action, 2df genotypic, Cochran-Armitage trend).\n"
"    * With 'fisher', Fisher's exact test is used to generate p-values.\n"
"    * 'perm' causes an adaptive permutation test to be performed.\n"
"    * 'mperm=[value]' causes a max(T) permutation test with the specified\n"
"      number of replications to be performed.\n"
"    * 'genedrop' causes offspring genotypes to be regenerated via gene-dropping\n"
"      in the permutation test.\n"
"    * 'perm-count' causes the permutation test report to include counts instead\n"
"      of frequencies.\n"
"    * 'counts' causes --assoc to report allele counts instead of frequencies.\n"
"    * 'p2' changes the --assoc permutation test used (see PLINK documentation).\n"
"    * 'dom', 'rec', 'gen', and 'trend' force the corresponding test to be used\n"
"      as the basis for --model permutation.  (By default, the most significant\n"
"      result among the allelic, dominant, and recessive tests is used.)\n"
"    * 'trend-only' causes only the trend test to be performed.\n"
"    Given a quantitative phenotype, --assoc performs a Wald test.  In this\n"
"    case, the 'qt-means' modifier causes trait means and standard deviations\n"
"    stratified by genotype to be reported as well.\n"
"    Several other flags (most notably, --aperm) can be used to customize the\n"
"    permutation test.\n\n"
	       );
    help_print("gxe\tmcovar", &help_ctrl, 1,
"  --gxe {covariate index}\n"
"    Given both a quantitative phenotype and a dichotomous covariate loaded with\n"
"    --covar defining two groups, --gxe compares the regression coefficient\n"
"    derived from considering only members of one group to the regression\n"
"    coefficient derived from considering only members of the other.  By\n"
"    default, the first covariate in the --covar file defines the groups; use\n"
"    e.g. '--gxe 3' to base them on the third covariate instead.\n\n"
	       );
    help_print(
"indep\tindep-pairwise", &help_ctrl, 1,
"  --indep [window size]<kb> [step size (markers)] [VIF threshold]\n"
"  --indep-pairwise [window size]<kb> [step size (markers)] [r^2 threshold]\n"
"    Generates a list of markers in approximate linkage equilibrium.  With the\n"
"    'kb' modifier, the window size is in kilobase units instead of marker\n"
"    count.  (Pre-'kb' space is optional, i.e. '--indep-pairwise 500 kb 5 0.5'\n"
"    and '--indep-pairwise 500kb 5 0.5' have the same effect.)\n"
"    Note that you need to rerun WDIST using --extract or --exclude on the\n"
"    .prune.in/.prune.out file to apply the list to another computation.\n\n"
		);
    help_print("make-rel", &help_ctrl, 1,
"  --make-rel <square | square0 | triangle> <gz | bin> <cov | ibc1 | ibc2>\n"
"             <single-prec>\n"
"    Writes a lower-triangular variance-standardized relationship (coancestry)\n"
"    matrix to {output prefix}.rel, and corresponding IDs to\n"
"    {output prefix}.rel.id.\n"
"    * 'square', 'square0', 'triangle', 'gz', and 'bin' act as they do on\n"
"      --distance.\n"
"    * The 'cov' modifier removes the variance standardization step, causing a\n"
"      covariance matrix to be calculated instead.\n"
"    * By default, the diagonal elements in the relationship matrix are based on\n"
"      --ibc's Fhat3; use the 'ibc1' or 'ibc2' modifiers to base them on Fhat1\n"
"      or Fhat2 instead.\n"
"    * WDIST normally performs this calculation with double-precision floating\n"
"      point numbers.  The 'single-prec' modifier switches to single-precision\n"
"      arithmetic, which is generally good enough; this decreases memory usage\n"
"      and speeds up computation.\n"
               );
    help_print("make-grm\tgrm", &help_ctrl, 1,
"  --make-grm <no-gz> <cov | ibc1 | ibc2> <single-prec>\n"
"    Writes the relationships in GCTA's gzipped list format, describing one pair\n"
"    per line.  Note that this file explicitly stores the number of valid\n"
"    observations (where neither individual has a missing call) for each pair,\n"
"    which is useful input for some scripts.\n\n"
	       );
    help_print("rel-cutoff\tgrm-cutoff", &help_ctrl, 1,
"  --rel-cutoff {val}\n"
"  --grm-cutoff {val}\n"
"    Excludes one member of each pair of individuals with relatedness greater\n"
"    than the given cutoff value (default 0.025).  If no later operation will\n"
"    cause the list of remaining individuals to be written to disk, this will\n"
"    save it to {output prefix}.rel.id.\n"
"    Note that maximizing the remaining sample size is equivalent to the NP-hard\n"
"    maximum independent set problem, so we use a greedy algorithm instead of\n"
"    guaranteeing optimality.  (Use the --make-rel and --keep/--remove flags if\n"
"    you want to try to do better.)\n\n"
	       );
    help_print("regress-distance", &help_ctrl, 1,
"  --regress-distance {iters} {d}\n"
"    Linear regression of pairwise genomic distances on pairwise average\n"
"    phenotypes and vice versa, using delete-d jackknife for standard errors.\n"
"    Scalar phenotype data is required.\n"
"    * With less than two parameters, d is set to {number of people}^0.6 rounded\n"
"      down.  With no parameters, 100k iterations are run.\n\n"
	       );
    help_print("regress-pcs\tdistance\tregress-pcs-distance", &help_ctrl, 1,
"  --regress-pcs [.evec or .eigenvec filename] <normalize-pheno> <sex-specific>\n"
"                <clip> {max PCs}\n"
"    Linear regression of phenotypes and genotypes on the given list of\n"
"    principal components (produced by SMARTPCA or GCTA).  Output is a .gen +\n"
"    .sample fileset in the Oxford IMPUTE/SNPTEST v2 format.\n"
"    * The 'normalize-pheno' modifier converts all phenotype residuals to\n"
"      Z-scores.  When combined with 'sex-specific', the Z-scores are evaluated\n"
"      separately by sex.\n"
"    * The 'clip' modifier clips out-of-range genotype residuals.  Without it,\n"
"      they are represented as negative probabilities in the .gen file, which\n"
"      are invalid input for some programs.\n"
"    * By default, principal components beyond the 20th are ignored; change this\n"
"      by setting the max PCs parameter.\n"
"  --regress-pcs-distance [.evec/.eigenvec file] <normalize-pheno>\n"
"                         <sex-specific> {max PCs} <square | square0 | triangle>\n"
"                         <gz | bin> <ibs> <1-ibs> <alct> <3d> <flat-missing>\n"
"    High-speed combination of --regress-pcs and --distance (no .gen + .sample\n"
"    fileset is written to disk).\n\n"
	       );
    help_print("regress-rel", &help_ctrl, 1,
"  --regress-rel {iters} {d}\n"
"    Linear regression of pairwise genomic relationships on pairwise average\n"
"    phenotypes, and vice versa.  Defaults for iters and d are the same as for\n"
"    --regress-distance.\n\n"
	       );
    help_print("groupdist", &help_ctrl, 1,
"  --groupdist {iters} {d}\n"
"    Considers three subsets of the distance matrix: pairs of affected\n"
"    individuals, affected-unaffected pairs, and pairs of unaffected\n"
"    individuals.  Each of these subsets has an average pairwise genomic\n"
"    distance; --groupdist computes the differences between those three\n"
"    averages, and estimates standard errors via delete-d jackknife.\n"
"    Dichotomous phenotype data is required.\n\n"
	       );
#ifndef NOLAPACK
    help_print("unrelated-heritability", &help_ctrl, 1,
"  --unrelated-heritability <strict> {tol} {initial covg} {initial covr}\n"
"    REML estimate of additive heritability, iterating with an accelerated\n"
"    variant of the EM algorithm until the rate of change of the log likelihood\n"
"    function is less than tol.  Scalar phenotype data is required.\n"
"    * The 'strict' modifier forces regular EM to be used.  tol defaults to\n"
"      10^{-7}, genomic covariance prior defaults to 0.45, and residual\n"
"      covariance prior defaults to (1 - covg).\n"
"    * For more details, see Vattikuti S, Guo J, Chow CC (2012) Heritability and\n"
"      Genetic Correlations Explained by Common SNPs for Metabolic Syndrome\n"
"      Traits.  PLoS Genet 8(3): e1002637.  doi:10.1371/journal.pgen.1002637\n\n"
	       );
#endif
    help_print("make-bed", &help_ctrl, 1,
"  --make-bed\n"
"    Creates a new binary fileset.  Unlike the automatic text-to-binary\n"
"    converter (which only respects --autosome and --chr), this supports all of\n"
"    WDIST's filtering flags.\n"
	       );
    help_print("recode\trecode12\ttab\ttranspose\trecode-lgen\trecodeAD\trecodead\trecodeA\trecodea\trecode-rlist\trecode-allele\tlist\twith-reference", &help_ctrl, 1,
"  --recode <12> <compound-genotypes | A | AD | lgen | lgen-ref | list | rlist |\n"
"           transpose> <tab | tabx | spacex>\n"
"    Creates a new text fileset with all filters applied.\n"
"    * The '12' modifier causes all alleles to be coded as 1s and 2s.\n"
"    * The 'compound-genotypes' modifier removes the space between pairs of\n"
"      genotype codes for the same marker.\n"
"    * The 'AD' modifier causes an additive + dominant component file, suitable\n"
"      for loading from R, to be generated instead.  If you don't want the\n"
"      dominant component, use 'A' instead.\n"
"    * The 'lgen' modifier causes a long-format fileset to be generated instead,\n"
"      (loadable with --lfile), while 'lgen-ref' generates a (usually) smaller\n"
"      long-format fileset loadable with --lfile + --reference.\n"
"    * The 'list' modifier creates a genotype-based list, while 'rlist'\n"
"      creates a rare-genotype fileset.  (These formats are not directly\n"
"      loadable by WDIST or PLINK.)\n"
"    * 'transpose' creates a transposed text fileset (loadable with --tfile).\n"
"    * The 'tab' modifier makes the output mostly tab-delimited instead of\n"
"      mostly space-delimited.  'tabx' and 'spacex' force all tabs and all\n"
"      spaces, respectively.\n\n"
	       );
    help_print("write-covar", &help_ctrl, 1,
"  --write-covar\n"
"    If a --covar file is loaded, this creates a revised covariate file after\n"
"    applying all filters.  (This automatically happens if --make-bed or\n"
"    --recode is specified.)\n\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode", &help_ctrl, 1,
"  --merge [.ped filename] [.map filename]\n"
"  --merge [text fileset prefix]\n"
"  --bmerge [.bed filename] [.bim filename] [.fam filename]\n"
"  --bmerge [binary fileset prefix]\n"
"    Merges the given fileset with the initially loaded fileset.  If you specify\n"
"    --make-bed or '--recode lgen', the initial merge result is written to\n"
"    {output prefix}-merge.bed + .bim + .fam, filtering is performed, and then\n"
"    the post-filtering data is written to {output prefix}.bed + .bim + .fam.\n"
"    Otherwise, the merged data is written directly to\n"
"    {output prefix}.bed + .bim + .fam.\n"
"  --merge-list [filename]\n"
"    Merge all filesets named in the text file with the initially loaded\n"
"    fileset.  The text file is interpreted as follows:\n"
"    * If a line contains only one name, it is assumed to be the prefix for a\n"
"      binary fileset.\n"
"    * If a line contains exactly two names, they are assumed to be the full\n"
"      filenames for a text fileset (.ped first, then .map).\n"
"    * If a line contains exactly three names, they are assumed to be the full\n"
"      filenames for a binary fileset (.bed, then .bim, then .fam).\n\n"
	       );
    help_print("write-snplist", &help_ctrl, 1,
"  --write-snplist\n"
"    Writes a .snplist file listing the names of all markers that pass the\n"
"    filters and inclusion thresholds you've specified.\n\n"
	       );
    if (!param_ct) {
      fputs(
"The following other flags are supported.  (Order of operations is described at\n"
"https://www.cog-genomics.org/wdist/order .)\n"
, stdout);
    }
    help_print("script", &help_ctrl, 0,
"  --script [fname] : Include command-line options from file.\n"
	       );
    help_print("rerun", &help_ctrl, 0,
"  --rerun {log}    : Rerun commands in log (default 'wdist.log').\n"
	       );
    help_print("no-fid", &help_ctrl, 0,
"  --no-fid         : .fam/.ped file does not contain column 1 (family ID).\n"
	       );
    help_print("no-parents", &help_ctrl, 0,
"  --no-parents     : .fam/.ped file does not contain columns 3-4 (parents).\n"
	       );
    help_print("no-sex", &help_ctrl, 0,
"  --no-sex         : .fam/.ped file does not contain column 5 (sex).\n"
	       );
    help_print("no-pheno", &help_ctrl, 0,
"  --no-pheno       : .fam/.ped file does not contain column 6 (phenotype).\n"
	       );
    help_print("load-dists\tgroupdist\tregress-distance", &help_ctrl, 0,
"  --load-dists [f] : Load a binary TRIANGULAR distance matrix for --groupdist\n"
"                     or --regress-distance analysis, instead of recalculating\n"
"                     it from scratch.\n"
	       );
    help_print("silent", &help_ctrl, 0,
"  --silent         : Suppress output to console.\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode", &help_ctrl, 0,
"  --merge-mode [n] : Adjust --merge/--bmerge/--merge-list behavior based on a\n"
"                     numeric code.\n"
"                     1 (default) = difference -> missing\n"
"                     2 = only overwrite originally missing calls\n"
"                     3 = only overwrite calls which are nonmissing in new file\n"
"                     4/5 = never overwrite and always overwrite, respectively\n"
"                     6 = report all mismatching calls without merging\n"
"                     7 = report mismatching nonmissing calls without merging\n"
	       );
    help_print("indiv-sort\tmerge\tbmerge\tmerge-list", &help_ctrl, 0,
"  --indiv-sort [m] : Specify family/individual ID sort order.  The following\n"
"                     three modes are currently supported:\n"
"                     * 'none'/'0' keeps individuals in the order they were\n"
"                       loaded.  This is the default for non-merge operations.\n"
"                     * 'natural'/'n' invokes \"natural sort\", e.g. 'id2' <\n"
"                       'ID3' < 'ID10'.  This is the default when merging.\n"
"                     * 'ascii'/'a' sorts in ASCII order, e.g. 'ID3' < 'id10' <\n"
"                       'id2'.\n"
"                     For now, only --make-bed and --merge/--bmerge/--merge-list\n"
"                     respect this flag.\n"
	       );
    help_print("pheno\tall-pheno", &help_ctrl, 0,
"  --pheno [fname]  : Specify alternate phenotype.\n"
"  --all-pheno      : For basic association tests, loop through all phenotypes\n"
"                     in --pheno file.\n"
	       );
    help_print("mpheno\tpheno", &help_ctrl, 0,
"  --mpheno [col]   : Specify phenotype column number in --pheno file.\n"
	       );
    help_print("pheno-name\tpheno", &help_ctrl, 0,
"  --pheno-name [c] : If phenotype file has a header row, use column with the\n"
"                     given name.\n"
	       );
    help_print("pheno-merge\tpheno", &help_ctrl, 0,
"  --pheno-merge    : If a phenotype is present in the original but not the\n"
"                     alternate file, use the original value instead of setting\n"
"                     the phenotype to missing.\n"
	       );
    help_print("covar\tcovar-name\tcovar-number", &help_ctrl, 0,
"  --covar [filename]    : Specify covariate file.\n"
"  --covar-name [names]  : Specifies covariate(s) in --covar file by name.\n"
"                          Separate multiple names with commas (spaces are not\n"
"                          permitted).  Use dashes to designate ranges.\n"
"  --covar-number [nums] : Specifies covariate(s) in --covar file by number.\n"
	       );
    help_print("within", &help_ctrl, 0,
"  --within [fname] : Specify clusters.\n"
	       );
    help_print("set\tsubset", &help_ctrl, 0,
"  --set [filename] : Specify sets.\n"
"  --subset [fname] : Specify list of subsets to extract from --set file.\n"
	       );
    help_print("prune", &help_ctrl, 0,
"  --prune          : Remove individuals with missing phenotypes.\n"
	       );
    help_print("1", &help_ctrl, 0,
"  --1              : Affection phenotypes are interpreted as 0 = unaffected,\n"
"                     1 = affected (instead of 0 = missing, 1 = unaffected,\n"
"                     2 = affected).\n"
	       );
// --map3 implicitly supported via autodetection
// --compound-genotypes automatically supported
    help_print("cow\tdog\thorse\tmouse\trice\tsheep", &help_ctrl, 0,
"  --cow/--dog/--horse/--mouse/--rice/--sheep : Specify nonhuman species.\n"
	       );
    help_print("autosome\tautosome-xy\tchr\tchr-excl", &help_ctrl, 0,
"  --autosome       : Exclude all non-autosomal markers.\n"
"  --autosome-xy    : Exclude all non-autosomal markers, except those with\n"
"                     chromosome code XY (pseudo-autosomal region of X).\n"
	       );
    help_print("chr\tchr-excl\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb", &help_ctrl, 0,
"  --chr [chrs...]  : Exclude all markers not on the given chromosome(s).  Valid\n"
"                     choices for humans are 0 (unplaced), 1-22, X, Y, XY, and\n"
"                     MT.  Separate multiple chromosomes with spaces and/or\n"
"                     commas, and use a dash (no adjacent spaces permitted) to\n"
"                     denote a range, e.g. '--chr 1-4, 22, xy'.\n"
"  --chr-excl [...] : Reverse of --chr (excludes markers on listed chromosomes).\n"
	       );
    help_print("d\tsnps", &help_ctrl, 0,
"  --d [char]       : Change marker range delimiter (usually '-').\n"
	       );
    help_print("from\tto\tsnp\twindow\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb", &help_ctrl, 0,
"  --from [mkr ID]  : Use ID(s) to specify a marker range to load.  When used\n"
"  --to   [mkr ID]    together, both markers must be on the same chromosome.\n"
"  --snp  [mkr ID]  : Specify a single marker to load.\n"
"  --window  [kbs]  : With --snp, loads all markers within half the specified kb\n"
"                     distance of the named marker.\n"
"  --from-bp [pos]  : Use physical position(s) to define a marker range to load.\n"
"  --to-bp   [pos]    --from-kb/--to-kb/--from-mb/--to-mb allow decimal values.\n"
"    ...              You're required to specify a single chromosome with these.\n"
	       );
    help_print("snps", &help_ctrl, 0,
"  --snps [IDs...]  : Use IDs to specify multiple marker ranges to load.  E.g.\n"
"                     '--snps rs1111-rs2222, rs3333, rs4444'.\n"
	       );
    help_print("maf", &help_ctrl, 0,
"  --maf {val}      : Minor allele frequency minimum threshold (default 0.01).\n"
"                     Note that the default threshold is only applied if --maf\n"
"                     is used without an accompanying value; if you do not\n"
"                     invoke --maf, no MAF inclusion threshold is applied.\n"
"                     Other inclusion thresholds work the same way.\n"
	       );
    help_print("max-maf\tmaf", &help_ctrl, 0,
"  --max-maf [val]  : Minor allele frequency maximum threshold.\n"
	       );
    help_print("geno\tmind", &help_ctrl, 0,
"  --geno {val}     : Maximum per-marker missing (default 0.1).\n"
"  --mind {val}     : Maximum per-person missing (default 0.1).\n"
	       );
    help_print("hwe", &help_ctrl, 0,
"  --hwe {val}      : Minimum Hardy-Weinberg disequilibrium p-value (exact),\n"
"                     default 0.001.\n"
	       );
    help_print("hwe-all\thwe", &help_ctrl, 0,
"  --hwe-all        : Given case-control data, include noncontrols in HWE test.\n"
	       );
    help_print("allow-no-sex\tmust-have-sex", &help_ctrl, 0,
"  --allow-no-sex   : Do not treat ambiguous-sex individuals as having missing\n"
"                     phenotypes in analysis commands.  (Automatic /w --no-sex.)\n"
"  --must-have-sex  : Force ambiguous-sex phenotypes to missing on --make-bed\n"
"                     and --recode.\n"
	       );
    help_print("set-hh-missing", &help_ctrl, 0,
"  --set-hh-missing : Cause --make-bed and --recode to set heterozygous haploid\n"
"                     genotypes to missing.\n"
	       );
    help_print("reference\tallele-count\tlfile", &help_ctrl, 0,
"  --reference [fn] : Specify default allele file for .lgen input.\n"
"  --allele-count   : When used with --reference, specifies that the .lgen file\n"
"                     contains reference allele counts.\n"
	       );
    help_print("with-phenotype\twrite-covar", &help_ctrl, 0,
"  --with-phenotype : Include phenotype when writing updated covariate file.\n"
	       );
    help_print("dummy-coding\twrite-covar", &help_ctrl, 0,
"  --dummy-coding   : Split categorical variables (n categories, for 2 < n < 50)\n"
"                     into n-1 binary dummy variables when writing updated\n"
"                     covariate file.\n"
	       );
    help_print("nonfounders", &help_ctrl, 0,
"  --nonfounders    : Include nonfounders in allele frequency/HWE calculations.\n"
	       );
    help_print("ppc-gap\tgenome\tZ-genome", &help_ctrl, 0,
"  --ppc-gap [val]  : Minimum number of base pairs, in thousands, between\n"
"                     informative pairs of markers used in --genome PPC test.\n"
"                     500 if unspecified.\n"
	       );
    help_print("seed", &help_ctrl, 0,
"  --seed [val]     : Set random number seed.\n"
	       );
    help_print("memory", &help_ctrl, 0,
"  --memory [val]   : Size, in MB, of initial malloc attempt.  (Some operating\n"
"                     systems allow this number to exceed total physical RAM.)\n"
	       );
    help_print("threads", &help_ctrl, 0,
"  --threads [val]  : Maximum number of concurrent threads.\n"
	       );
    help_print("debug", &help_ctrl, 0,
"  --debug          : Enable debug logging.\n"
	       );
    help_print("extract\texclude", &help_ctrl, 0,
"  --extract [file] : Exclude all markers not in the given list.\n"
"  --exclude [file] : Exclude all markers in the given list.\n"
	       );
    help_print("keep\tremove", &help_ctrl, 0,
"  --keep [fname]   : Exclude all individuals not in the given list.\n"
"  --remove [fname] : Exclude all individuals in the given list.\n"
	       );
    help_print("maf-succ", &help_ctrl, 0,
"  --maf-succ       : Rule of succession MAF estimation (used in EIGENSTRAT).\n"
"                     Given j observations of one allele and k >= j observations\n"
"                     of the other, infer a MAF of (j+1) / (j+k+2), rather than\n"
"                     the usual j / (j+k).\n"
	       );
    help_print("ci", &help_ctrl, 0,
"  --ci [size]      : Set size of odds ratio confidence interval to report.\n"
	       );
    help_print("pfilter", &help_ctrl, 0,
"  --pfilter [val]  : Filter out association test results with higher p-values.\n"
	       );
    help_print("loop-assoc", &help_ctrl, 0,
"  --loop-assoc [f] : Run specified case/control association commands once for\n"
"                     each cluster in the file, using cluster membership as the\n"
"                     phenotype.\n"
	       );
    help_print("cell\tmodel", &help_ctrl, 0,
"  --cell [thresh]  : Specify contingency table threshold for performing all\n"
"                     --model tests.\n"
	       );
    help_print("adjust\tlambda", &help_ctrl, 0,
"  --lambda [val]   : Set genomic control lambda for --adjust.\n"
	       );
    help_print("exponent\tdistance", &help_ctrl, 0,
"  --exponent [val] : When computing genomic distances, each marker has a weight\n"
"                     of (2q(1-q))^{-val}, where q is the inferred MAF.  (Use\n"
"                     --read-freq if you want to explicitly specify some or all\n"
"                     of the MAFs.)\n"
	       );
    help_print("recode\trecode-allele", &help_ctrl, 0,
"  --recode-allele [f] : With --recode A or --recode AD, count alleles named in\n"
"                        the file (instead of the minor allele).\n"
	       );
    help_print("keep-allele-order\tmake-bed\tmerge\tbmerge\tmerge-list", &help_ctrl, 0,
"  --keep-allele-order : Keep the original allele order when creating a new\n"
"                        fileset, instead of forcing A2 to be the major allele.\n"
	       );
    help_print("update-chr\tupdate-cm\tupdate-map\tupdate-name", &help_ctrl, 0,
"  --update-map [f] <chr | cm | name> : Update marker information.\n"
	       );
    help_print("update-alleles", &help_ctrl, 0,
"  --update-alleles [fname]  : Update marker allele codes.\n"
	       );
    help_print("update-ids\tupdate-parents\tupdate-sex", &help_ctrl, 0,
"  --update-ids [fname]      : Update individual IDs, parental IDs, or sexes\n"
"  --update-parents [fname]    using the information in the provided file.  See\n"
"  --update-sex [fname]        the PLINK documentation for file format details.\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode\tmerge-allow-equal-pos", &help_ctrl, 0,
"  --merge-allow-equal-pos   : Do not merge markers with different names but\n"
"                              identical positions.\n"
	       );
    help_print("allele1234\talleleACGT\talleleacgt", &help_ctrl, 0,
"  --allele1234 <multichar>  : Interpret/recode A/C/G/T alleles as 1/2/3/4.\n"
"                              With 'multichar', converts all A/C/G/Ts in\n"
"                              allele names to 1/2/3/4s.\n"
"  --alleleACGT <multichar>  : Reverse of --allele1234.\n"
	       );
    help_print("reference-allele\tupdate-ref-allele", &help_ctrl, 0,
"  --reference-allele [file] : Force alleles named in the file to A1.\n"
	       );
    help_print("read-freq\tupdate-freq", &help_ctrl, 0,
"  --read-freq [filename]    : Loads MAFs from the given PLINK-style or --freqx\n"
"  --update-freq [filename]    frequency file, instead of just setting them to\n"
"                              frequencies observed in the .ped/.bed file.\n"
	       );
    help_print("parallel", &help_ctrl, 0,
"  --parallel [k] [n]        : Divide the output matrix into n pieces, and only\n"
"                              compute the kth piece.  The primary output file\n"
"                              will have the piece number included in its name,\n"
"                              e.g. wdist.rel.13 or wdist.rel.13.gz if k is 13.\n"
"                              Concatenating these files in order will yield the\n"
"                              full matrix of interest.  (Yes, this can be done\n"
"                              before unzipping.)\n"
"                              N.B. This cannot be used to directly write a\n"
"                              symmetric square matrix.  Choose the square0 or\n"
"                              triangle format instead, and postprocess as\n"
"                              necessary.\n"
	       );
    help_print("filter", &help_ctrl, 0,
"  --filter [filename] [val] : Filter individuals based on a phenotype value.\n"
	       );
    help_print("mfilter\tfilter", &help_ctrl, 0,
"  --mfilter [col]           : Specify alternate column in --filter file.\n"
	       );
    help_print("filter-cases\tfilter-controls", &help_ctrl, 0,
"  --filter-cases            : Include only cases.\n"
"  --filter-controls         : Include only controls.\n"
	       );
    help_print("filter-males\tfilter-females", &help_ctrl, 0,
"  --filter-males            : Include only males.\n"
"  --filter-females          : Include only females.\n"
	       );
    help_print("filter-founders\tfilter-nonfounders", &help_ctrl, 0,
"  --filter-founders         : Include only founders.\n"
"  --filter-nonfounders      : Include only nonfounders.\n"
	       );
    help_print("missing-genotype\tmissing-phenotype", &help_ctrl, 0,
"  --missing-genotype [char] : Code for missing genotype (normally '0').\n"
"  --missing-phenotype [val] : Numeric code for missing phenotype (normally -9).\n"
	       );
    help_print("output-missing-genotype\toutput-missing-phenotype", &help_ctrl, 0,
"  --output-missing-genotype [ch] : Code for missing genotype when creating new\n"
"                                   text fileset (--recode).\n"
"  --output-missing-phenotype [n] : Numeric code for missing phenotype when\n"
"                                   creating new fileset (--make-bed/--recode).\n"
	       );
    help_print("missing-code\tmissing_code\tmissing-phenotype", &help_ctrl, 0,
"  --missing-code {vals}     : Comma-separated list of missing phenotype values,\n"
"  --missing_code {vals}       for Oxford-formatted filesets (normally 'NA').\n"
	       );
    help_print("make-pheno", &help_ctrl, 0,
"  --make-pheno [file] [val] : Specify dichotomous phenotype, where cases have\n"
"                              the given value.  If the value is '*', all\n"
"                              individuals present in the phenotype file are\n"
"                              affected (and other individuals in the .ped/.fam\n"
"                              are unaffected).\n"
	       );
    help_print("tail-pheno\tgroupdist", &help_ctrl, 0,
"  --tail-pheno [Ltop] {Hbt} : Form 'low' (<= Ltop, unaffected) and 'high'\n"
"                              (greater than Hbt, affected) groups from scalar\n"
"                              phenotype data.  If Hbt is unspecified, it is set\n"
"                              equal to Ltop.  Central phenotype values are\n"
"                              treated as missing.\n"
	       );
    help_print("adjust\tgc\tlog10\tqq-plot", &help_ctrl, 0,
"  --adjust <gc> <log10> <qq-plot> : Report some multiple-testing corrections.\n"
	       );

    help_print("aperm", &help_ctrl, 0,
"  --aperm [min perms] [max perms] [alpha] [beta] [init interval] [slope]\n"
"    This sets six parameters controlling adaptive permutation tests.\n\n"
	       );
    if (!param_ct) {
      fputs(
"For further documentation and support, consult the main webpage\n"
"(https://www.cog-genomics.org/wdist ) and/or the wdist-users mailing list\n"
"(https://groups.google.com/d/forum/wdist-users ).\n"
, stdout);
    }
  } while (help_ctrl.iters_left--);
  if (help_ctrl.unmatched_ct) {
    net_unmatched_ct = help_ctrl.unmatched_ct;
    printf("\nNo help entr%s for", (help_ctrl.unmatched_ct == 1)? "y" : "ies");
    col_num = (help_ctrl.unmatched_ct == 1)? 17 : 19;
    arg_uidx = 0;
    while (help_ctrl.unmatched_ct) {
      arg_uidx = next_non_set_unsafe(help_ctrl.all_match_arr, arg_uidx);
      help_ctrl.unmatched_ct--;
      if (help_ctrl.unmatched_ct) {
	if (net_unmatched_ct == 2) {
	  if (help_ctrl.param_lens[arg_uidx] + col_num > 76) {
	    printf("\n'%s'", argv[arg_uidx]);
	    col_num = 2 + help_ctrl.param_lens[arg_uidx];
	  } else {
	    printf(" '%s'", argv[arg_uidx]);
	    col_num += 3 + help_ctrl.param_lens[arg_uidx];
	  }
	} else {
	  if (help_ctrl.param_lens[arg_uidx] + col_num > 75) {
	    printf("\n'%s',", argv[arg_uidx]);
	    col_num = 3 + help_ctrl.param_lens[arg_uidx];
	  } else {
	    printf(" '%s',", argv[arg_uidx]);
	    col_num += 4 + help_ctrl.param_lens[arg_uidx];
	  }
	}
	if (help_ctrl.unmatched_ct == 1) {
	  if (col_num > 76) {
	    fputs("\nor", stdout);
	    col_num = 2;
	  } else {
	    fputs(" or", stdout);
	    col_num += 3;
	  }
	}
      } else {
	if (help_ctrl.param_lens[arg_uidx] + col_num > 75) {
	  printf("\n'%s'.\n", argv[arg_uidx]);
	} else {
	  printf(" '%s'.\n", argv[arg_uidx]);
	}
      }
      arg_uidx++;
    }
  }
  if (param_ct) {
    while (0) {
    disp_help_ret_NOMEM:
      retval = RET_NOMEM;
    }
    free_cond(help_ctrl.param_lens);
    free_cond(help_ctrl.all_match_arr);
    if (help_ctrl.argv && (help_ctrl.argv != argv)) {
      free(help_ctrl.argv);
    }
  }
  return retval;
}

intptr_t malloc_size_mb = 0;

unsigned char* wkspace;

void dispmsg(int32_t retval) {
  switch (retval) {
  case RET_NOMEM:
    logprint("Error: Out of memory.  Try the --memory and/or --parallel flags.\n");
    break;
  case RET_WRITE_FAIL:
    logprint("\nError: File write failure.\n");
    break;
  case RET_READ_FAIL:
    logprint("\nError: File read failure.\n");
    break;
  }
}

double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2) {
  // This function implements an exact SNP test of Hardy-Weinberg
  // Equilibrium as described in Wigginton, JE, Cutler, DJ, and
  // Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
  // Equilibrium. American Journal of Human Genetics. 76: 000 - 000.
  //
  // The original version was written by Jan Wigginton.
  //
  // This version was written by Christopher Chang.  It contains the following
  // improvements over the original SNPHWE():
  // - Proper handling of >64k genotypes.  Previously, there was a potential
  // integer overflow.
  // - Detection and efficient handling of floating point overflow and
  // underflow.  E.g. instead of summing a tail all the way down, the loop
  // stops once the latest increment underflows the partial sum's 53-bit
  // precision; this results in a large speedup when max heterozygote count
  // >1k.
  // - No malloc() call.  It's only necessary to keep track of a few partial
  // sums.
  //
  // Note that the SNPHWE_t() function below is a lot more efficient for
  // testing against a p-value inclusion threshold.  SNPHWE2() should only be
  // used if you need the actual p-value.
  intptr_t obs_homc;
  intptr_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  int64_t rare_copies = 2LL * obs_homr + obs_hets;
  int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
  double curr_hets_t2 = obs_hets;
  double curr_homr_t2 = obs_homr;
  double curr_homc_t2 = obs_homc;
  double tailp = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double centerp = 0;
  double lastp2 = tailp;
  double lastp1 = tailp;
  double curr_hets_t1;
  double curr_homr_t1;
  double curr_homc_t1;
  double preaddp;
  if (!genotypes2) {
    return 1;
  }

  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    // tail 1 = upper
    while (curr_hets_t2 > 1.5) {
      // het_probs[curr_hets] = 1
      // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      if (lastp2 < EXACT_TEST_BIAS) {
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if (centerp == 0) {
      return 1;
    }
    while (curr_hets_t2 > 1.5) {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    curr_hets_t1 = obs_hets + 2;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    while (curr_homr_t1 > 0.5) {
      // het_probs[curr_hets + 2] = het_probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
      lastp1 *= (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      curr_hets_t1 += 2;
      curr_homr_t1 -= 1;
      curr_homc_t1 -= 1;
    }
  } else {
    // tail 1 = lower
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      if (lastp2 < EXACT_TEST_BIAS) {
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if (centerp == 0) {
      return 1;
    }
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    curr_hets_t1 = obs_hets;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    while (curr_hets_t1 > 1.5) {
      curr_homr_t1 += 1;
      curr_homc_t1 += 1;
      lastp1 *= (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      curr_hets_t1 -= 2;
    }
  }
  return tailp / (tailp + centerp);
}

int32_t SNPHWE_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh) {
  // Threshold-test-only version of SNPHWE2() which is usually able to exit
  // from the calculation earlier.  Returns 0 if these counts are close enough
  // to Hardy-Weinberg equilibrium, 1 otherwise.
  //
  // Suppose, for definiteness, that the number of observed hets is no less
  // than expectation.  (Same ideas apply for the other case.)  We proceed as
  // follows:
  // - Sum the *relative* likelihoods of more likely smaller het counts.
  // - Determine the minimum tail mass to pass the threshold.
  // - The majority of the time, the tail boundary elements are enough to pass
  // the threshold; we never need to sum the remainder of the tails.
  // - And in the case of disequilibrium, we will often be able to immediately
  // determine that the tail sum cannot possibly pass the threshold, just by
  // looking at the tail boundary elements and using a geometric series to
  // upper-bound the tail sums.
  // - Only when neither of these conditions hold do we start traveling down
  // the tails.
  intptr_t obs_homc;
  intptr_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  int64_t rare_copies = 2LL * obs_homr + obs_hets;
  int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
  double curr_hets_t2 = obs_hets; // tail 2
  double curr_homr_t2 = obs_homr;
  double curr_homc_t2 = obs_homc;

  // Subtract epsilon from initial probability mass, so that we can compare to
  // 1 when determining tail vs. center membership without floating point error
  // biting us in the ass
  double tailp1 = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double centerp = 0;
  double lastp2 = tailp1;
  double tailp2 = 0;
  double tail1_ceil;
  double tail2_ceil;
  double lastp1;
  double curr_hets_t1;
  double curr_homr_t1;
  double curr_homc_t1;

  // Initially, if center sum reaches this, the test can immediately fail.
  // Once center is summed, this is recalculated, and when tail sum has reached
  // this, we've passed.
  double exit_thresh;
  double exit_threshx;
  double ratio;
  double preaddp;
  if (!genotypes2) {
    return 0;
  }

  // Convert thresh into reverse odds ratio.
  thresh = (1 - thresh) / thresh;

  // Expected het count:
  //   2 * rarefreq * (1 - rarefreq) * genotypes
  // = 2 * (rare_copies / (2 * genotypes)) * (1 - rarefreq) * genotypes
  // = rare_copies * (1 - (rare_copies / (2 * genotypes)))
  // = (rare_copies * (2 * genotypes - rare_copies)) / (2 * genotypes)
  // 
  // The computational identity is
  //   P(nhets == n) := P(nhets == n+2) * (n+2) * (n+1) /
  //                    (4 * homr(n) * homc(n))
  // where homr() and homc() are the number of homozygous rares/commons needed
  // to maintain the same allele frequencies.
  // This probability is always decreasing when proceeding away from the
  // expected het count.

  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    // tail 1 = upper
    if (obs_hets < 2) {
      return 0;
    }

    // An initial upper bound on the tail sum is useful, since it lets us
    // report test failure before summing the entire center.
    // The immediate tail is easy: if the next element out on the tail is r,
    // then 1 + r + r^2 + ... = 1 / (1-r) works.
    // For the far tail, we currently use the trivial bound of
    //   1 + floor[het_exp_floor / 2]
    // (each far tail element must be no greater than 1 and there are at
    // most that many of them).  This bound could be improved, but it might not
    // be worth the precomputational effort.
    // ...and as long as we're using such a weak bound for the far tail,
    // there's no point to carefully calculating the near tail since we're very
    // unlikely to use the partial result before exiting from the function.
    exit_thresh = rare_copies * thresh * EXACT_TEST_BIAS;

    // het_probs[curr_hets] = 1
    // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1) / (4 * (curr_homr + 1) * (curr_homc + 1))
    do {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      if (lastp2 < EXACT_TEST_BIAS) {
	tailp2 = lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp > exit_thresh) {
	return 1;
      }
    } while (curr_hets_t2 > 1.5);
    exit_thresh = centerp / thresh;
    if (tailp1 + tailp2 >= exit_thresh) {
      return 0;
    }
    // c + cr + cr^2 + ... = c/(1-r), which is an upper bound for the tail sum
    ratio = (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * (curr_homr_t2 + 1) * (curr_homc_t2 + 1));
    tail2_ceil = tailp2 / (1 - ratio);
    curr_hets_t1 = obs_hets + 2;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    lastp1 = (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
    tail1_ceil = tailp1 / (1 - lastp1);
    if (tail1_ceil + tail2_ceil < exit_thresh) {
      return 1;
    }
    lastp1 *= tailp1;
    tailp1 += lastp1;

    if (obs_homr > 1) {
      // het_probs[curr_hets + 2] = het_probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
      exit_threshx = exit_thresh - tailp2;
      do {
	curr_hets_t1 += 2;
	curr_homr_t1 -= 1;
	curr_homc_t1 -= 1;
	lastp1 *= (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
	preaddp = tailp1;
	tailp1 += lastp1;
	if (tailp1 > exit_threshx) {
	  return 0;
	}
	if (tailp1 <= preaddp) {
	  break;
	}
      } while (curr_homr_t1 > 1.5);
    }
    if (tailp1 + tail2_ceil < exit_thresh) {
      return 1;
    }
    exit_threshx = exit_thresh - tailp1;
    while (curr_hets_t2 > 1) {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      preaddp = tailp2;
      tailp2 += lastp2;
      if (tailp2 >= exit_threshx) {
	return 0;
      }
      if (tailp2 <= preaddp) {
	return 1;
      }
      curr_hets_t2 -= 2;
    }
    return 1;
  } else {
    // tail 1 = lower
    if (!obs_homr) {
      return 0;
    }
    exit_thresh = rare_copies * thresh * EXACT_TEST_BIAS;
    do {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      if (lastp2 < EXACT_TEST_BIAS) {
	tailp2 = lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp > exit_thresh) {
	return 1;
      }
    } while (curr_homr_t2 > 0.5);
    exit_thresh = centerp / thresh;
    if (tailp1 + tailp2 >= exit_thresh) {
      return 0;
    }
    ratio = (4 * curr_homr_t2 * curr_homc_t2) / ((curr_hets_t2 + 2) * (curr_hets_t2 + 1));
    tail2_ceil = tailp2 / (1 - ratio);
    curr_hets_t1 = obs_hets;
    curr_homr_t1 = obs_homr + 1;
    curr_homc_t1 = obs_homc + 1;
    lastp1 = (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
    tail1_ceil = tailp1 / (1 - lastp1);
    lastp1 *= tailp1;
    tailp1 += lastp1;

    if (tail1_ceil + tail2_ceil < exit_thresh) {
      return 1;
    }
    if (obs_hets >= 4) {
      exit_threshx = exit_thresh - tailp2;
      do {
	curr_hets_t1 -= 2;
	curr_homr_t1 += 1;
	curr_homc_t1 += 1;
	lastp1 *= (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
	preaddp = tailp1;
        tailp1 += lastp1;
	if (tailp1 > exit_threshx) {
	  return 0;
	}
	if (tailp1 <= preaddp) {
	  break;
	}
      } while (curr_hets_t1 > 3.5);
    }
    if (tailp1 + tail2_ceil < exit_thresh) {
      return 1;
    }
    exit_threshx = exit_thresh - tailp1;
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      preaddp = tailp2;
      tailp2 += lastp2;
      if (tailp2 >= exit_threshx) {
	return 0;
      }
      if (tailp2 <= preaddp) {
	return 1;
      }
    }
    return 1;
  }
}

// back to our regular program

int32_t strcmp_casted(const void* s1, const void* s2) {
  return strcmp((char*)s1, (char*)s2);
}

inline char* read_next_unsafe(char* target, char* source) {
  // assumes space- or tab-termination
  while ((*source != ' ') && (*source != '\t')) {
    *target++ = *source++;
  }
  return target;
}

inline char* read_next_unsafe_upd(char* target, char** source_ptr) {
  // assumes space- or tab-termination
  while ((**source_ptr != ' ') && (**source_ptr != '\t')) {
    *target++ = **source_ptr;
    *source_ptr += 1;
  }
  return target;
}

inline char* read_next_upd(char* target, char** source_ptr) {
  while (!is_space_or_eoln(**source_ptr)) {
    *target++ = **source_ptr;
    *source_ptr += 1;
  }
  return target;
}

inline int32_t is_contained(char* id_buf, char* lptr, int32_t max_id_len, int32_t filter_line_ct, char* fam_id, char* indiv_id) {
  return (bsearch_fam_indiv(id_buf, lptr, max_id_len, filter_line_ct, fam_id, indiv_id) != -1);
}

int32_t determine_max_id_len(FILE* filterfile, char* filterval, int32_t mfilter_col, int32_t* filter_line_ct_ptr) {
  int32_t cur_max = 4;
  char* bufptr;
  int32_t ii;
  int32_t jj;

  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Excessively long line in filter file (max %d chars).\n", MAXLINELEN - 3);
      logprintb();
      return -1;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    ii = 2 + strlen_se(bufptr);
    bufptr = next_item(bufptr);
    if (no_more_items_kns(bufptr)) {
      logprint(errstr_filter_format);
      return -1;
    }
    ii += strlen_se(bufptr);
    if (ii > cur_max) {
      cur_max = ii;
    }
    if (filterval) {
      for (jj = 0; jj < mfilter_col; jj++) {
	bufptr = next_item(bufptr);
      }
      if (no_more_items_kns(bufptr)) {
        logprint(errstr_filter_format);
	return -1;
      }
      if (!strncmp(filterval, bufptr, jj)) {
	*filter_line_ct_ptr += 1;
      }
    } else {
      *filter_line_ct_ptr += 1;
    }
  }
  return cur_max;
}

char* resize_id_buf(char* id_buf, int32_t max_id_len, int32_t max_pid_len) {
  if (max_pid_len) {
    if (max_id_len > max_pid_len) {
      free(id_buf);
    } else {
      return id_buf;
    }
  }
  return (char*)malloc(max_id_len * sizeof(char));
}

double calc_wt_mean(double exponent, int32_t lhi, int32_t lli, int32_t hhi) {
  double lcount = (double)lli + ((double)lhi * 0.5);
  int64_t tot = lhi + lli + hhi;
  double dtot = (double)tot;
  int64_t subcount = lli; // avoid 32-bit integer overflow
  double weight;
  if ((!lhi) && ((!lli) || (!hhi))) {
    return 0.0;
  }
  weight = pow(2 * lcount * (dtot - lcount) / (dtot * dtot), -exponent);
  subcount = lhi * (subcount + hhi) + 2 * subcount * hhi;
  return (subcount * weight * 2) / (double)(tot * tot);
}

double calc_wt_mean_maf(double exponent, double maf) {
  // assume Hardy-Weinberg equilibrium
  // homozygote frequencies: maf^2, (1-maf)^2
  // heterozygote frequency: 2maf(1-maf)
  double ll_freq = maf * maf;
  double lh_freq = 2 * maf * (1.0 - maf);
  double hh_freq = (1.0 - maf) * (1.0 - maf);
  double weight;
  if (lh_freq == 0.0) {
    return 0.0;
  }
  weight = pow(lh_freq, -exponent);
  return (lh_freq * (ll_freq + lh_freq) + 2 * ll_freq * hh_freq) * weight;
}



// The following functions may come in handy if we need to store large array(s)
// of auxiliary ternary data in memory.
//
// inline void set_base4(uintptr_t* base4_arr, int32_t loc, int32_t val) {
//   int32_t shift_num = (loc % BITCT2) * 2;
//   uintptr_t ulii = (3 * ONELU) << shift_num;
//   uintptr_t* base4_ptr = &(base4_arr[loc / BITCT2]);
//   *base4_ptr = (*base4_ptr & (~ulii)) | (val << shift_num);
// }
//
// inline int32_t get_base4(uintptr_t* base4_arr, int32_t loc) {
//   int32_t shift_num = (loc % BITCT2) * 2;
//   return ((base4_arr[loc / BITCT2] >> shift_num) & (3 * ONELU));
// }

int32_t populate_pedigree_rel_info(Pedigree_rel_info* pri_ptr, uintptr_t unfiltered_indiv_ct, char* id_buf, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info) {
  unsigned char* wkspace_mark;
  unsigned char* wkspace_mark2;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t mm;
  int32_t nn;
  int32_t oo;
  uint32_t uii;
  uint32_t initial_family_blocks;
  uintptr_t indiv_uidx;
  int64_t llii;
  char* family_ids;
  char* cur_person_id;
  char* last_family_id = NULL;
  char* cur_family_id;
  uint32_t* family_sizes;
  uint32_t* uiptr;
  uint32_t* uiptr2 = NULL;
  uint32_t fidx;
  int32_t family_size;
  uint32_t* remaining_indiv_idxs;
  int32_t* remaining_indiv_parent_idxs; // -1 = no parent (or nonshared)
  uint32_t remaining_indiv_ct;
  uint32_t indiv_idx_write;
  intptr_t max_family_id_len = 0;
  char* indiv_ids;
  intptr_t max_indiv_id_len = 0;
  int32_t max_pm_id_len;
  uint32_t family_id_ct;
  uint32_t* fis_ptr;
  char* stray_parent_ids;
  int32_t stray_parent_ct;
  uintptr_t* processed_indivs;
  int32_t founder_ct;
  int32_t max_family_nf = 0;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ctlm = unfiltered_indiv_ctl * BITCT;
  uint32_t* complete_indiv_idxs;
  uintptr_t complete_indiv_idx_ct;
  double* rs_ptr;
  double* rel_writer;
  double dxx;
  double* tmp_rel_space = NULL;
  double* tmp_rel_writer = NULL;

  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    jj = strlen_se(&(person_ids[indiv_uidx * max_person_id_len])) + 1;
    if (jj > max_family_id_len) {
      max_family_id_len = jj;
    }
    jj = strlen_se(&(person_ids[indiv_uidx * max_person_id_len + jj]));
    if (jj >= max_indiv_id_len) {
      max_indiv_id_len = jj + 1;
    }
  }
  if (max_paternal_id_len > max_maternal_id_len) {
    max_pm_id_len = max_paternal_id_len;
  } else {
    max_pm_id_len = max_maternal_id_len;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_space), unfiltered_indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_rel_nf_idxs), unfiltered_indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_idxs), unfiltered_indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_c_checked(&family_ids, unfiltered_indiv_ct * max_family_id_len)) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&family_sizes, unfiltered_indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }

  // copy all the items over in order, then qsort, then eliminate duplicates
  // and count family sizes.
  cur_family_id = family_ids;
  cur_person_id = person_ids;
  uiptr = family_sizes;
  *uiptr = 1;
  jj = strlen_se(cur_person_id);
  memcpy(cur_family_id, cur_person_id, jj);
  cur_family_id[jj] = '\0';
  for (indiv_uidx = 1; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    cur_person_id = &(cur_person_id[max_person_id_len]);
    mm = strlen_se(cur_person_id);
    if ((jj != mm) || memcmp(cur_family_id, cur_person_id, mm)) {
      cur_family_id = &(cur_family_id[max_family_id_len]);
      memcpy(cur_family_id, cur_person_id, mm);
      cur_family_id[mm] = '\0';
      jj = mm;
      *(++uiptr) = 1;
    } else {
      *uiptr += 1;
    }
  }
  initial_family_blocks = 1 + (uint32_t)(uiptr - family_sizes);
  if (qsort_ext(family_ids, initial_family_blocks, max_family_id_len, strcmp_deref, (char*)family_sizes, sizeof(int32_t))) {
    return RET_NOMEM;
  }

  last_family_id = family_ids;
  cur_family_id = &(family_ids[max_family_id_len]);
  family_id_ct = 1;
  uii = 1; // read idx
  if (initial_family_blocks != 1) {
    uiptr = family_sizes;
    while (strcmp(cur_family_id, last_family_id)) {
      family_id_ct++;
      uiptr++;
      last_family_id = cur_family_id;
      cur_family_id = &(cur_family_id[max_family_id_len]);
      uii++;
      if (uii == initial_family_blocks) {
	break;
      }
    }
    if (uii < initial_family_blocks) {
      uiptr2 = uiptr; // family_sizes read pointer
      *uiptr += *(++uiptr2);
      uii++;
      cur_family_id = &(cur_family_id[max_family_id_len]); // read pointer
      while (uii < initial_family_blocks) {
	while (!strcmp(cur_family_id, last_family_id)) {
	  *uiptr += *(++uiptr2);
	  uii++;
	  if (uii == initial_family_blocks) {
	    break;
	  }
	  cur_family_id = &(cur_family_id[max_family_id_len]);
	}
	if (uii < initial_family_blocks) {
	  *(++uiptr) = *(++uiptr2);
	  last_family_id = &(last_family_id[max_family_id_len]);
	  strcpy(last_family_id, cur_family_id);
	  family_id_ct++;
	  uii++;
	  cur_family_id = &(cur_family_id[max_family_id_len]);
	}
      }
    }
  }

  uiptr = family_sizes;
  if (family_id_ct < unfiltered_indiv_ct) {
    wkspace_reset(family_ids);
    family_ids = (char*)wkspace_alloc(family_id_ct * max_family_id_len);
    family_sizes = (uint32_t*)wkspace_alloc(family_id_ct * sizeof(int32_t));
    for (uii = 0; uii < family_id_ct; uii++) {
      family_sizes[uii] = *uiptr++;
    }
  }
  pri_ptr->family_ids = family_ids;
  pri_ptr->family_id_ct = family_id_ct;
  pri_ptr->max_family_id_len = max_family_id_len;
  pri_ptr->family_sizes = family_sizes;

  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_offsets), (family_id_ct + 1) * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_rel_space_offsets), (family_id_ct + 1) * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_founder_cts), family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  fill_int_zero((int32_t*)(pri_ptr->family_founder_cts), family_id_ct);

  ii = 0; // running family_info offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_info_offsets[fidx] = ii;
    ii += family_size;
  }

  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  memcpy(uiptr, pri_ptr->family_info_offsets, family_id_ct * sizeof(int32_t));

  // Fill family_idxs, family_founder_cts, and founder portion of
  // family_rel_nf_idxs.
  cur_person_id = person_ids;
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    jj = strlen_se(cur_person_id);
    memcpy(id_buf, cur_person_id, jj);
    id_buf[jj] = '\0';
    kk = bsearch_str(id_buf, family_ids, max_family_id_len, 0, family_id_ct - 1);
    pri_ptr->family_idxs[indiv_uidx] = kk;
    if (is_founder(founder_info, indiv_uidx)) {
      pri_ptr->family_founder_cts[kk] += 1;
      pri_ptr->family_rel_nf_idxs[indiv_uidx] = uiptr[kk];
      uiptr[kk] += 1;
    }
    cur_person_id = &(cur_person_id[max_person_id_len]);
  }
  wkspace_reset(uiptr);
  jj = 0; // running rel_space offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_rel_space_offsets[fidx] = jj;
    kk = pri_ptr->family_founder_cts[fidx];
    if (family_size - kk > max_family_nf) {
      max_family_nf = family_size - kk;
    }
    // No need to explicitly store the (kk * (kk - 1)) / 2 founder-founder
    // relationships.
    jj += ((int64_t)family_size * (family_size - 1) - (int64_t)kk * (kk - 1)) / 2;
  }
  // make it safe to determine size of blocks by subtracting from the next
  // offset, even if we're at the last family
  pri_ptr->family_info_offsets[family_id_ct] = unfiltered_indiv_ct;
  pri_ptr->family_rel_space_offsets[family_id_ct] = jj;
  if (wkspace_alloc_d_checked(&(pri_ptr->rel_space), jj * sizeof(double))) {
    return RET_NOMEM;
  }

  wkspace_mark = wkspace_base;
  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  // populate family_info_space
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    uiptr[fidx] = pri_ptr->family_info_offsets[fidx];
  }
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    fidx = pri_ptr->family_idxs[indiv_uidx];
    pri_ptr->family_info_space[uiptr[fidx]] = indiv_uidx;
    uiptr[fidx] += 1;
  }
  wkspace_reset(wkspace_mark);

  if (wkspace_alloc_ul_checked(&processed_indivs, (unfiltered_indiv_ctl + (max_family_nf + (BITCT2 - 1)) / BITCT2) * sizeof(intptr_t))) {
    return RET_NOMEM;
  }
  fill_ulong_one(&(processed_indivs[unfiltered_indiv_ctl]), (max_family_nf + (BITCT2 - 1)) / BITCT2);

  wkspace_mark2 = wkspace_base;
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = family_sizes[fidx];
    founder_ct = pri_ptr->family_founder_cts[fidx];
    remaining_indiv_ct = family_size - founder_ct;
    stray_parent_ct = 0;
    if (remaining_indiv_ct) {
      memcpy(processed_indivs, founder_info, unfiltered_indiv_ctl * sizeof(intptr_t));
      if (wkspace_alloc_ui_checked(&complete_indiv_idxs, family_size * sizeof(int32_t))) {
	return RET_NOMEM;
      }
      if (wkspace_alloc_ui_checked(&remaining_indiv_idxs, remaining_indiv_ct * sizeof(int32_t))) {
	return RET_NOMEM;
      }
      if (wkspace_alloc_c_checked(&indiv_ids, family_size * max_indiv_id_len)) {
	return RET_NOMEM;
      }
      if (wkspace_alloc_i_checked(&remaining_indiv_parent_idxs, remaining_indiv_ct * 2 * sizeof(int32_t))) {
	return RET_NOMEM;
      }
      if (wkspace_alloc_c_checked(&stray_parent_ids, remaining_indiv_ct * 2 * max_pm_id_len)) {
	return RET_NOMEM;
      }
      ii = pri_ptr->family_info_offsets[fidx];
      fis_ptr = &(pri_ptr->family_info_space[ii]);
      rs_ptr = &(pri_ptr->rel_space[pri_ptr->family_rel_space_offsets[fidx]]);
      rel_writer = rs_ptr;
      cur_person_id = indiv_ids;
      for (ii = 0; ii < family_size; ii++) {
	kk = fis_ptr[ii];
	jj = strlen_se(&(person_ids[kk * max_person_id_len])) + 1;
	strcpy(cur_person_id, &(person_ids[kk * max_person_id_len + jj]));
	cur_person_id = &(cur_person_id[max_indiv_id_len]);
      }

      // Compile list of non-founder family member indices, and identify
      // parents who are referred to by individual ID but are NOT in the
      // dataset.
      ii = 0; // family_info_space index
      complete_indiv_idx_ct = 0;
      cur_person_id = stray_parent_ids;
      for (uii = 0; uii < remaining_indiv_ct; uii++) {
	while (is_set(founder_info, fis_ptr[ii])) {
	  complete_indiv_idxs[complete_indiv_idx_ct++] = fis_ptr[ii];
	  ii++;
	}
	kk = fis_ptr[ii++];

	// does not track sex for now
	if (memcmp("0", &(paternal_ids[kk * max_paternal_id_len]), 2)) {
	  mm = bsearch_str(&(paternal_ids[kk * max_paternal_id_len]), indiv_ids, max_indiv_id_len, 0, family_size - 1);
	  if (mm == -1) {
	    strcpy(cur_person_id, &(paternal_ids[kk * max_paternal_id_len]));
	    cur_person_id = &(cur_person_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_indiv_parent_idxs[uii * 2] = -2;
	  } else {
            remaining_indiv_parent_idxs[uii * 2] = fis_ptr[mm];
	  }
	} else {
          remaining_indiv_parent_idxs[uii * 2] = -1;
	}
	if (memcmp("0", &(maternal_ids[kk * max_maternal_id_len]), 2)) {
	  mm = bsearch_str(&(maternal_ids[kk * max_maternal_id_len]), indiv_ids, max_indiv_id_len, 0, family_size - 1);
	  if (mm == -1) {
	    strcpy(cur_person_id, &(maternal_ids[kk * max_maternal_id_len]));
	    cur_person_id = &(cur_person_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_indiv_parent_idxs[uii * 2 + 1] = -2;
	  } else {
	    remaining_indiv_parent_idxs[uii * 2 + 1] = fis_ptr[mm];
	  }
	} else {
	  remaining_indiv_parent_idxs[uii * 2 + 1] = -1;
	}
        remaining_indiv_idxs[uii] = kk;
      }
      while (ii < family_size) {
	complete_indiv_idxs[complete_indiv_idx_ct++] = fis_ptr[ii];
	ii++;
      }
      qsort(stray_parent_ids, stray_parent_ct, max_pm_id_len, strcmp_casted);
      cur_person_id = stray_parent_ids;
      ii = 0; // read idx
      jj = 0; // write idx

      // Now filter out all such parents who aren't referenced at least twice.
      while (ii < stray_parent_ct - 1) {
        if (strcmp(&(stray_parent_ids[ii * max_pm_id_len]), &(stray_parent_ids[(ii + 1) * max_pm_id_len]))) {
	  ii++;
	  continue;
	}
	ii++;
	strcpy(cur_person_id, &(stray_parent_ids[ii * max_pm_id_len]));
	do {
	  ii++;
        } while (!(strcmp(cur_person_id, &(stray_parent_ids[ii * max_pm_id_len])) || (ii >= stray_parent_ct - 1)));
        cur_person_id = &(cur_person_id[max_pm_id_len]);
	jj++;
      }
      stray_parent_ct = jj;

      // Now allocate temporary relatedness table between nonfounders and
      // stray parents with multiple references.
      if (stray_parent_ct) {
        if (wkspace_alloc_d_checked(&tmp_rel_space, (family_size - founder_ct) * stray_parent_ct * sizeof(double))) {
	  return RET_NOMEM;
        }
	tmp_rel_writer = tmp_rel_space;
      }

      // Now fill in remainder of remaining_indiv_parent_idxs.
      for (uii = 0; uii < remaining_indiv_ct; uii++) {
	jj = remaining_indiv_idxs[uii];
	if (remaining_indiv_parent_idxs[uii * 2] == -2) {
	  kk = bsearch_str(&(paternal_ids[jj * max_paternal_id_len]), stray_parent_ids, max_pm_id_len, 0, stray_parent_ct - 1);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[uii * 2] = kk;
	}
	if (remaining_indiv_parent_idxs[uii * 2 + 1] == -2) {
	  kk = bsearch_str(&(maternal_ids[jj * max_maternal_id_len]), stray_parent_ids, max_pm_id_len, 0, stray_parent_ct - 1);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[uii * 2 + 1] = kk;
	}
      }
      llii = (int64_t)founder_ct * (founder_ct - 1);
      while (remaining_indiv_ct) {
	indiv_idx_write = 0;
	for (uii = 0; uii < remaining_indiv_ct; uii++) {
	  kk = remaining_indiv_parent_idxs[uii * 2];
	  mm = remaining_indiv_parent_idxs[uii * 2 + 1];
	  jj = remaining_indiv_idxs[uii];
	  if (((kk == -1) || is_set(processed_indivs, kk)) && ((mm == -1) || is_set(processed_indivs, mm))) {
	    for (nn = 0; nn < founder_ct; nn++) {
	      // relationship between kk and nnth founder
	      if ((kk >= (int)unfiltered_indiv_ct) || (kk == -1)) {
		dxx = 0.0;
	      } else if (is_set(founder_info, kk)) {
		if (kk == (int)complete_indiv_idxs[nn]) {
		  dxx = 0.5;
		} else {
		  dxx = 0.0;
		}
	      } else {
		oo = pri_ptr->family_rel_nf_idxs[kk];
                dxx = 0.5 * rs_ptr[((int64_t)oo * (oo - 1) - llii) / 2 + nn];
	      }
	      if (is_set(founder_info, mm)) {
		if (mm == (int)complete_indiv_idxs[nn]) {
		  dxx += 0.5;
		}
	      } else if ((mm != -1) && (mm < (int)unfiltered_indiv_ct)) {
		oo = pri_ptr->family_rel_nf_idxs[mm];
		dxx += 0.5 * rs_ptr[((int64_t)oo * (oo - 1) - llii) / 2 + nn];
	      }
	      *rel_writer++ = dxx;
	    }
	    for (; nn < (int)complete_indiv_idx_ct; nn++) {
	      if (kk == -1) {
		dxx = 0.0;
	      } else if (kk >= (int)unfiltered_indiv_ct) {
		dxx = 0.5 * tmp_rel_space[(nn - founder_ct) * stray_parent_ct + kk - unfiltered_indiv_ctlm];
	      } else if (is_set(founder_info, kk)) {
                dxx = 0.5 * rs_ptr[((int64_t)nn * (nn - 1) - llii) / 2 + pri_ptr->family_rel_nf_idxs[kk]];
	      } else {
		oo = pri_ptr->family_rel_nf_idxs[kk];
		if (oo == nn) {
		  dxx = 0.5;
		} else if (oo < nn) {
		  dxx = 0.5 * rs_ptr[((int64_t)nn * (nn - 1) - llii) / 2 + oo];
		} else {
		  dxx = 0.5 * rs_ptr[((int64_t)oo * (oo - 1) - llii) / 2 + nn];
		}
	      }
	      if (mm >= (int)unfiltered_indiv_ct) {
		dxx += 0.5 * tmp_rel_space[(nn - founder_ct) * stray_parent_ct + mm - unfiltered_indiv_ctlm];
	      } else if (is_set(founder_info, mm)) {
		dxx += 0.5 * rs_ptr[((int64_t)nn * (nn - 1) - llii) / 2 + pri_ptr->family_rel_nf_idxs[mm]];
	      } else if (mm != -1) {
		oo = pri_ptr->family_rel_nf_idxs[mm];
		if (oo == nn) {
		  dxx += 0.5;
		} else if (oo < nn) {
		  dxx += 0.5 * rs_ptr[((int64_t)nn * (nn - 1) - llii) / 2 + oo];
		} else {
		  dxx += 0.5 * rs_ptr[((int64_t)oo * (oo - 1) - llii) / 2 + nn];
		}
	      }
	      *rel_writer++ = dxx;
	    }
	    for (nn = 0; nn < stray_parent_ct; nn++) {
	      if (kk >= (int)unfiltered_indiv_ct) {
		if (kk == nn + (int)unfiltered_indiv_ctlm) {
		  dxx = 0.5;
		} else {
		  dxx = 0.0;
		}
	      } else if (kk == -1) {
                dxx = 0.0;
	      } else {
		oo = pri_ptr->family_rel_nf_idxs[kk];
		if (oo < founder_ct) {
		  dxx = 0.0;
		} else {
                  dxx = 0.5 * tmp_rel_space[(oo - founder_ct) * stray_parent_ct + nn];
		}
	      }
	      if (mm >= (int)unfiltered_indiv_ct) {
		if (mm == nn + (int)unfiltered_indiv_ctlm) {
		  dxx += 0.5;
		}
	      } else if (mm != -1) {
		oo = pri_ptr->family_rel_nf_idxs[mm];
		if (oo >= founder_ct) {
		  dxx += 0.5 * tmp_rel_space[(oo - founder_ct) * stray_parent_ct + nn];
		}
	      }
	      *tmp_rel_writer++ = dxx;
	    }
	    pri_ptr->family_rel_nf_idxs[jj] = complete_indiv_idx_ct;
	    complete_indiv_idxs[complete_indiv_idx_ct] = jj;
	    set_bit(processed_indivs, jj, &complete_indiv_idx_ct);
	  } else {
            remaining_indiv_parent_idxs[indiv_idx_write * 2] = kk;
	    remaining_indiv_parent_idxs[indiv_idx_write * 2 + 1] = mm;
	    remaining_indiv_idxs[indiv_idx_write++] = jj;
	  }
	}
	if (indiv_idx_write == remaining_indiv_ct) {
	  logprint("Error: Pedigree graph is cyclic.  Check for evidence of time travel abuse in\nyour cohort.\n");
	  return RET_INVALID_FORMAT;
	}
	remaining_indiv_ct = indiv_idx_write;
      }
      wkspace_reset(wkspace_mark2);
    }
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

void count_genders(uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, int32_t* male_ct_ptr, int32_t* female_ct_ptr, int32_t* unk_ct_ptr) {
  int32_t male_ct = 0;
  int32_t female_ct = 0;
  int32_t unk_ct = 0;
  uint32_t unfiltered_indiv_ctld = unfiltered_indiv_ct / BITCT;
  uint32_t unfiltered_indiv_ct_rem = unfiltered_indiv_ct & (BITCT - 1);
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t indiv_bidx;
  for (indiv_bidx = 0; indiv_bidx < unfiltered_indiv_ctld; indiv_bidx++) {
    ulii = ~(*indiv_exclude++);
  count_genders_last_loop:
    uljj = *sex_nm++;
    unk_ct += popcount_long(ulii & (~uljj));
    ulii &= uljj;
    uljj = *sex_male++;
    male_ct += popcount_long(ulii & uljj);
    female_ct += popcount_long(ulii & (~uljj));
  }
  if (unfiltered_indiv_ct_rem) {
    ulii = (~(*indiv_exclude)) & ((ONELU << unfiltered_indiv_ct_rem) - ONELU);
    unfiltered_indiv_ct_rem = 0;
    goto count_genders_last_loop;
  }
  *male_ct_ptr = male_ct;
  *female_ct_ptr = female_ct;
  *unk_ct_ptr = unk_ct;
}

int32_t load_pheno(FILE* phenofile, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_exclude_ct, char* sorted_person_ids, uintptr_t max_person_id_len, int32_t* id_map, int32_t missing_pheno, int32_t missing_pheno_len, int32_t affection_01, uint32_t mpheno_col, char* phenoname_str, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr) {
  int32_t affection = 1;
  uintptr_t* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  int32_t header_processed = 0;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  uintptr_t* isz = NULL;
  int32_t person_idx;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  uint32_t tmp_len;
  uint32_t tmp_len2;
  uint32_t uii;
  double dxx;
  double dyy;
  if (pheno_d) {
    affection = 0;
  } else {
    if (wkspace_alloc_ul_checked(&isz, unfiltered_indiv_ctl * sizeof(intptr_t))) {
      return RET_NOMEM;
    }
    fill_ulong_zero(isz, unfiltered_indiv_ctl);
    if (!pheno_c) {
      pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
      if (!pheno_c) {
	return RET_NOMEM;
      }
      *pheno_c_ptr = pheno_c;
    }
  }
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    return RET_NOMEM;
  }
  // ----- phenotype file load -----
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Excessively long line in phenotype file (max %d chars).\n", MAXLINELEN - 3);
      logprintb();
      return RET_INVALID_FORMAT;
    }
    bufptr0 = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr0)) {
      continue;
    }
    tmp_len = strlen_se(bufptr0);
    bufptr = next_item(bufptr0);
    if (no_more_items_kns(bufptr)) {
      logprint(errstr_phenotype_format);
      logprint("At least two items expected in line.\n");
      sprintf(logbuf, "Original line: %s", bufptr0);
      logprintb();
      return RET_INVALID_FORMAT;
    }
    tmp_len2 = strlen_se(bufptr);
    if (!header_processed) {
      if (phenoname_str || ((tmp_len == 3) && (tmp_len2 == 3) && (!memcmp("FID", bufptr0, 3)) && (!memcmp("IID", bufptr, 3)))) {
	if (phenoname_str) {
	  tmp_len = strlen(phenoname_str);
	  do {
	    bufptr = next_item(bufptr);
	    if (no_more_items_kns(bufptr)) {
	      logprint("Error: --pheno-name column not found.\n");
	      return RET_INVALID_FORMAT;
	    }
	    mpheno_col++;
	    tmp_len2 = strlen_se(bufptr);
	  } while ((tmp_len2 != tmp_len) || memcmp(bufptr, phenoname_str, tmp_len));
	}
      } else {
	header_processed = 1;
        if (!mpheno_col) {
	  mpheno_col = 1;
        }
      }
    }
    if (!header_processed) {
      header_processed = 1;
    } else {
      person_idx = bsearch_fam_indiv(id_buf, sorted_person_ids, max_person_id_len, indiv_ct, bufptr0, bufptr);
      if (person_idx != -1) {
	person_idx = id_map[person_idx];
	bufptr = next_item_mult(bufptr, mpheno_col);
	if (no_more_items_kns(bufptr)) {
	  logprint(errstr_phenotype_format);
	  logprint("Fewer entries than expected in line.\n");
	  sprintf(logbuf, "Original line: %s", bufptr0);
	  logprintb();
	  return RET_INVALID_FORMAT;
	}
	if (affection) {
	  if (eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
	    if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
	      // Since we're only making one pass through the file, we don't
	      // have the luxury of knowing in advance whether the phenotype is
	      // binary or scalar.  If there is a '0' entry that occurs before
	      // we know the phenotype is scalar, we need to not set the
	      // phenotype to zero during the binary -> scalar conversion step.
	      if (*bufptr == '0') {
		set_bit_noct(isz, person_idx);
	      }
	    } else if (affection_01) {
	      if (*bufptr == '0') {
		clear_bit_noct(pheno_c, person_idx);
	      } else {
		set_bit_noct(pheno_c, person_idx);
	      }
	    } else {
	      if (*bufptr == '1') {
		clear_bit_noct(pheno_c, person_idx);
	      } else {
		set_bit_noct(pheno_c, person_idx);
	      }
	    }
	  } else {
	    pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
	    if (!pheno_d) {
	      return RET_NOMEM;
	    }
	    *pheno_d_ptr = pheno_d;
	    if (affection_01) {
	      dxx = 0.0;
	      dyy = 1.0;
	    } else {
	      dxx = 1.0;
	      dyy = 2.0;
	    }
	    for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
	      if (is_set(isz, uii)) {
		pheno_d[uii] = 0.0;
		set_bit_noct(pheno_nm, uii);
	      } else if (is_set(pheno_nm, uii)) {
		pheno_d[uii] = is_set(pheno_c, uii)? dyy : dxx;
	      }
	    }
	    free(pheno_c);
	    *pheno_c_ptr = NULL;
	    affection = 0;
	  }
	}
	if (!affection) {
	  if (sscanf(bufptr, "%lg", &dxx) == 1) {
	    pheno_d[person_idx] = dxx;
	    set_bit_noct(pheno_nm, person_idx);
	  }
	}
      }
    }
  }
  if (!feof(phenofile)) {
    return RET_READ_FAIL;
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t makepheno_load(FILE* phenofile, char* makepheno_str, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_exclude_ct, char* sorted_person_ids, uintptr_t max_person_id_len, int32_t* id_map, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr) {
  uintptr_t indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  uint32_t mp_strlen = strlen(makepheno_str);
  int32_t makepheno_all = ((mp_strlen == 1) && (makepheno_str[0] == '*'));
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t* pheno_c = *pheno_c_ptr;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  int32_t person_idx;
  uint32_t tmp_len;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    return RET_NOMEM;
  }
  if (!pheno_c) {
    pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
    if (!pheno_c) {
      return RET_NOMEM;
    }
    fill_ulong_zero(pheno_c, unfiltered_indiv_ctl);
    *pheno_c_ptr = pheno_c;
  }
  if (makepheno_all) {
    fill_ulong_one(pheno_nm, unfiltered_indiv_ctl);
    zero_trailing_bits(pheno_nm, unfiltered_indiv_ct);
  }
  tbuf[MAXLINELEN - 1] = ' ';  
  while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Excessively long line in phenotype file (max %d chars).\n", MAXLINELEN - 3);
      logprintb();
      return RET_INVALID_FORMAT;
    }
    bufptr0 = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr0)) {
      continue;
    }
    bufptr = next_item(bufptr0);
    if (no_more_items_kns(bufptr)) {
      logprint(errstr_phenotype_format);
      return RET_INVALID_FORMAT;
    }
    person_idx = bsearch_fam_indiv(id_buf, sorted_person_ids, max_person_id_len, indiv_ct, bufptr0, bufptr);
    if (person_idx != -1) {
      person_idx = id_map[person_idx];
      if (makepheno_all) {
	set_bit_noct(pheno_c, person_idx);
      } else {
	set_bit_noct(pheno_nm, person_idx);
	bufptr = next_item(bufptr);
        tmp_len = strlen_se(bufptr);
	if ((tmp_len == mp_strlen) && (!memcmp(bufptr, makepheno_str, mp_strlen))) {
	  set_bit_noct(pheno_c, person_idx);
	}
      }
    }
  }
  if (!feof(phenofile)) {
    return RET_READ_FAIL;
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t convert_tail_pheno(uintptr_t unfiltered_indiv_ct, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, double tail_bottom, double tail_top, double missing_phenod) {
  uintptr_t* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  uint32_t uii;
  double dxx;
  if (!(*pheno_d_ptr)) {
    logprint("Error: --tail-pheno requires scalar phenotype data.\n");
    return RET_INVALID_FORMAT;
  }
  if (!pheno_c) {
    pheno_c = (uintptr_t*)malloc(((unfiltered_indiv_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t));
    if (!pheno_c) {
      return RET_NOMEM;
    }
    *pheno_c_ptr = pheno_c;
  }
  for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
    if (is_set(pheno_nm, uii)) {
      dxx = pheno_d[uii];
      if (dxx <= tail_bottom) {
	clear_bit_noct(pheno_c, uii);
      } else if (dxx > tail_top) {
        set_bit_noct(pheno_c, uii);
      } else {
	clear_bit_noct(pheno_nm, uii);
      }
    }
  }
  free(pheno_d);
  *pheno_d_ptr = NULL;
  return 0;
}

void prune_missing_phenos(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, uintptr_t* pheno_nm, double* pheno_d, double missing_phenod) {
  uint32_t uii;
  // could be a lot more efficient
  if (!pheno_d) {
    for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
      if (!is_set(pheno_nm, uii)) {
	set_bit(indiv_exclude, uii, indiv_exclude_ct_ptr);
      }
    }
  } else {
    for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
      if (pheno_d[uii] == missing_phenod) {
        set_bit(indiv_exclude, uii, indiv_exclude_ct_ptr);
      }
    }
  }
}

const char acgtstr[] = "ACGT";

void convert_multichar_allele(char* cptr, int32_t allelexxxx) {
  char cc = *cptr;
  if ((allelexxxx & 1) && (cptr[1] != '\0')) {
    return;
  }
  if (allelexxxx < 3) {
    while (cc) {
      if (cc == 'A') {
	*cptr = '1';
      } else if (cc == 'C') {
	*cptr = '2';
      } else if (cc == 'G') {
	*cptr = '3';
      } else if (cc == 'T') {
	*cptr = '4';
      }
      cc = *(++cptr);
    }
    return;
  }
  while (cc) {
    if ((cc >= '1') && (cc <= '4')) {
      *cptr = acgtstr[cc - '1'];
    }
    cc = *(++cptr);
  }
}

void allelexxxx_recode(int32_t allelexxxx, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_exclude, uint32_t marker_ct) {
  uintptr_t marker_uidx = 0;
  char* cptr;
  char cc;
  uintptr_t marker_idx;
  if (max_marker_allele_len == 1) {
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      cptr = &(marker_alleles[marker_uidx * 2]);
      cc = *cptr;
      if (allelexxxx < 3) {
	*cptr++ = convert_to_1234(cc);
	*cptr = convert_to_1234(*cptr);
      } else {
	*cptr++ = convert_to_acgt(cc);
	*cptr = convert_to_acgt(*cptr);
      }
      marker_uidx++;
    }
  } else {
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      cptr = &(marker_alleles[marker_uidx * 2 * max_marker_allele_len]);
      convert_multichar_allele(cptr, allelexxxx);
      cptr = &(cptr[max_marker_allele_len]);
      convert_multichar_allele(cptr, allelexxxx);
      marker_uidx++;
    }
  }
}

inline void set_surrounding_bits(char* id_buf, int32_t id_idx, uintptr_t* exclude_arr, char* sorted_ids, uintptr_t max_id_len, int32_t* id_map, int32_t sorted_ids_len, uintptr_t* exclude_ct_ptr) {
  uint32_t id_len = strlen(id_buf) + 1;
  int32_t id_idx2 = id_idx - 1;
  while ((id_idx2 >= 0) && (!memcmp(id_buf, &(sorted_ids[id_idx2 * max_id_len]), id_len))) {
    set_bit(exclude_arr, id_map[id_idx2--], exclude_ct_ptr);
  }
  id_idx2 = id_idx + 1;
  while ((id_idx2 < sorted_ids_len) && (!memcmp(id_buf, &(sorted_ids[id_idx2 * max_id_len]), id_len))) {
    set_bit(exclude_arr, id_map[id_idx2++], exclude_ct_ptr);
  }
}

inline void clear_surrounding_bits(char* id_buf, int32_t id_idx, uintptr_t* exclude_arr_new, char* sorted_ids, uintptr_t max_id_len, int32_t* id_map, int32_t sorted_ids_len, uintptr_t* include_ct_ptr) {
  uint32_t id_len = strlen(id_buf) + 1;
  int32_t id_idx2 = id_idx - 1;
  while ((id_idx2 >= 0) && (!memcmp(id_buf, &(sorted_ids[id_idx2 * max_id_len]), id_len))) {
    clear_bit(exclude_arr_new, id_map[id_idx2--], include_ct_ptr);
  }
  id_idx2 = id_idx + 1;
  while ((id_idx2 < sorted_ids_len) && (!memcmp(id_buf, &(sorted_ids[id_idx2 * max_id_len]), id_len))) {
    clear_bit(exclude_arr_new, id_map[id_idx2++], include_ct_ptr);
  }
}

int32_t include_or_exclude(char* fname, char* sorted_ids, int32_t sorted_ids_len, uintptr_t max_id_len, int32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr, int32_t indivs, int32_t do_exclude, int32_t duplicates_present) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t* exclude_arr_new = NULL;
  uintptr_t include_ct = 0;
  uintptr_t unfiltered_ctl = (unfiltered_ct + (BITCT - 1)) / BITCT;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  int32_t ii;
  int32_t jj;

  if (!do_exclude) {
    if (wkspace_alloc_ul_checked(&exclude_arr_new, unfiltered_ctl * sizeof(intptr_t))) {
      return RET_NOMEM;
    }
    fill_ulong_one(exclude_arr_new, unfiltered_ctl);
  }
  if (fopen_checked(&infile, fname, "r")) {
    return RET_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  if (indivs) {
    if (wkspace_alloc_c_checked(&id_buf, max_id_len)) {
      return RET_NOMEM;
    }
    while (fgets(tbuf, MAXLINELEN, infile) != NULL) {
      if (!tbuf[MAXLINELEN - 1]) {
	sprintf(logbuf, "Error: Excessively long line in --keep/--remove file (max %d chars).\n", MAXLINELEN - 3);
	logprintb();
        return RET_INVALID_FORMAT;
      }
      bufptr0 = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr0)) {
	continue;
      }
      bufptr = next_item(tbuf);
      if (no_more_items_kns(bufptr)) {
	logprint("Error: Improperly formatted --keep/--remove file.\n");
        return RET_INVALID_FORMAT;
      }
      ii = bsearch_fam_indiv(id_buf, sorted_ids, max_id_len, sorted_ids_len, tbuf, bufptr);
      if (ii != -1) {
        jj = id_map[ii];
        if (do_exclude) {
          set_bit(exclude_arr, jj, exclude_ct_ptr);
	  if (duplicates_present) {
            set_surrounding_bits(id_buf, ii, exclude_arr, sorted_ids, max_id_len, id_map, sorted_ids_len, exclude_ct_ptr);
	  }
	} else if (!is_set(exclude_arr, jj)) {
	  clear_bit(exclude_arr_new, jj, &include_ct);
	  if (duplicates_present) {
	    clear_surrounding_bits(id_buf, ii, exclude_arr_new, sorted_ids, max_id_len, id_map, sorted_ids_len, &include_ct);
	  }
	}
      }
    }
  } else {
    while (fgets(tbuf, MAXLINELEN, infile) != NULL) {
      read_next_terminate(tbuf, tbuf);
      if (!(*tbuf)) {
	continue;
      }
      ii = bsearch_str(tbuf, sorted_ids, max_id_len, 0, sorted_ids_len - 1);
      if (ii != -1) {
	jj = id_map[ii];
	if (do_exclude) {
	  set_bit(exclude_arr, jj, exclude_ct_ptr);
	  if (duplicates_present) {
            set_surrounding_bits(tbuf, ii, exclude_arr, sorted_ids, max_id_len, id_map, sorted_ids_len, exclude_ct_ptr);
	  }
	} else if (!is_set(exclude_arr, jj)) {
	  clear_bit(exclude_arr_new, jj, &include_ct);
	  if (duplicates_present) {
	    clear_surrounding_bits(tbuf, ii, exclude_arr_new, sorted_ids, max_id_len, id_map, sorted_ids_len, &include_ct);
	  }
	}
      }
    }
  }
  if (!feof(infile)) {
    return RET_READ_FAIL;
  }
  if (!do_exclude) {
    memcpy(exclude_arr, exclude_arr_new, unfiltered_ctl * sizeof(intptr_t));
    *exclude_ct_ptr = unfiltered_ct - include_ct;
  }
  wkspace_reset(wkspace_mark);
  fclose(infile);
  return 0;
}

int32_t filter_indivs_file(char* filtername, char* sorted_person_ids, int32_t sorted_ids_len, uintptr_t max_person_id_len, int32_t* id_map, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* filterval, int32_t mfilter_col) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  int32_t fv_strlen = strlen(filterval);
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t include_ct = 0;
  uintptr_t* indiv_exclude_new;
  char* id_buf;
  char* bufptr;
  int32_t person_idx;
  int32_t ii;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&indiv_exclude_new, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    return RET_NOMEM;
  }
  fill_ulong_one(indiv_exclude_new, unfiltered_indiv_ctl);

  if (fopen_checked(&infile, filtername, "r")) {
    return RET_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Excessively long line in --keep/--remove file (max %d chars).\n", MAXLINELEN - 3);
      logprintb();
      fclose(infile);
      return RET_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr = next_item(bufptr);
    if (no_more_items_kns(bufptr)) {
      logprint("Error: Improperly formatted --filter file.\n");
      fclose(infile);
      return RET_INVALID_FORMAT;
    }
    person_idx = bsearch_fam_indiv(id_buf, sorted_person_ids, max_person_id_len, sorted_ids_len, tbuf, bufptr);
    if (person_idx != -1) {
      person_idx = id_map[person_idx];
      if (!is_set(indiv_exclude, person_idx)) {
	for (ii = 0; ii < mfilter_col; ii++) {
	  bufptr = next_item(bufptr);
	}
	if (no_more_items_kns(bufptr)) {
	  logprint("Error: Improperly formatted --filter file.\n");
	  fclose(infile);
	  return RET_INVALID_FORMAT;
	}
        if ((!memcmp(filterval, bufptr, fv_strlen)) && is_space_or_eoln(bufptr[fv_strlen])) {
          clear_bit(indiv_exclude_new, person_idx, &include_ct);
	}
      }
    }
  }
  if (!feof(infile)) {
    fclose(infile);
    return RET_READ_FAIL;
  }
  memcpy(indiv_exclude, indiv_exclude_new, unfiltered_indiv_ctl * sizeof(intptr_t));
  *indiv_exclude_ct_ptr = unfiltered_indiv_ct - include_ct;

  wkspace_reset(wkspace_mark);
  fclose(infile);
  return 0;
}

void filter_indivs_bitfields(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, uintptr_t* orfield, int32_t orfield_flip, uintptr_t* ornot) {
  // indiv_exclude := indiv_exclude | orfield | (~ornot) if !orfield_flip
  //               := indiv_exclude | (~orfield) | (~ornot) otherwise
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t* ieptr = indiv_exclude;
  uintptr_t* ieend = &(indiv_exclude[unfiltered_indiv_ctl]);
  if (orfield_flip) {
    if (ornot) {
      do {
	*ieptr |= (~(*orfield++)) | (~(*ornot++));
      } while (++ieptr < ieend);
    } else {
      do {
	*ieptr |= ~(*orfield++);
      } while (++ieptr < ieend);
    }
  } else {
    if (ornot) {
      do {
	*ieptr |= (*orfield++) | (~(*ornot++));
      } while (++ieptr < ieend);
    } else {
      do {
	*ieptr |= (*orfield++) | (~(*ornot++));
      } while (++ieptr < ieend);
    }
  }
  zero_trailing_bits(indiv_exclude, unfiltered_indiv_ct);
  *indiv_exclude_ct_ptr = popcount_longs(indiv_exclude, 0, unfiltered_indiv_ctl);
}

int32_t mind_filter(FILE* bedfile, double mind_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, int32_t bed_offset, char missing_geno) {
  uint32_t mind_int_thresh = (int)(mind_thresh * (unfiltered_marker_ct - marker_exclude_ct) + EPSILON);
  uint32_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ct2l = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t removed_ct = 0;
  uintptr_t* loadbuf;
  uintptr_t* lptr;
  uint32_t* missing_cts;
  uintptr_t indiv_uidx;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uint32_t uii;
  uint32_t ujj;
  uintptr_t ulii;
  int32_t retval = 0;

  if (wkspace_alloc_ui_checked(&missing_cts, unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ct2l * sizeof(intptr_t))) {
    goto mind_filter_ret_NOMEM;
  }
  loadbuf[unfiltered_indiv_ct2l - 1] = 0;
  fill_uint_zero(missing_cts, unfiltered_indiv_ct);
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto mind_filter_ret_READ_FAIL;
  }
  marker_idx = 0;
  for (marker_uidx = 0; marker_idx < marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto mind_filter_ret_READ_FAIL;
      }
    }
    marker_idx++;
    if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
      goto mind_filter_ret_READ_FAIL;
    }
    lptr = loadbuf;
    indiv_uidx = 0;
    ujj = unfiltered_indiv_ct2l * BITCT2;
    for (uii = 0; uii < ujj; uii += BITCT2) {
      ulii = *lptr++;
      ulii = (ulii & FIVEMASK) & ((~ulii) >> 1);
      // now ulii has single bit set only at missing positions
      while (ulii) {
	missing_cts[uii + CTZLU(ulii) / 2] += 1;
	ulii &= ulii - 1;
      }
    }
  }
  indiv_uidx = next_non_set(indiv_exclude, 0, unfiltered_indiv_ct);
  while (indiv_uidx < unfiltered_indiv_ct) {
    if (missing_cts[indiv_uidx] > mind_int_thresh) {
      set_bit_noct(indiv_exclude, indiv_uidx);
      removed_ct++;
    }
    indiv_uidx = next_non_set(indiv_exclude, indiv_uidx + 1, unfiltered_indiv_ct);
  }
  *indiv_exclude_ct_ptr += removed_ct;
  sprintf(logbuf, "%u %s removed due to missing genotype data (--mind).\n", removed_ct, species_str(removed_ct));
  logprintb();
  while (0) {
  mind_filter_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  mind_filter_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t incr_text_allele(char cc, char* marker_alleles, uint32_t* marker_allele_cts, int32_t is_founder, uint32_t* marker_nf_allele_cts) {
  int32_t ii;
  for (ii = 0; ii < 4; ii++) {
    if (marker_alleles[ii] == '\0') {
      marker_alleles[ii] = cc;
      if (is_founder) {
	marker_allele_cts[ii] = 1;
      } else {
        marker_nf_allele_cts[ii] = 1;
      }
      return 0;
    } else if (marker_alleles[ii] == cc) {
      if (is_founder) {
        marker_allele_cts[ii] += 1;
      } else {
	marker_nf_allele_cts[ii] += 1;
      }
      return 0;
    }
  }
  return -1;
}

#ifdef __LP64__
void freq_hwe_haploid_count_120v(__m128i* vptr, __m128i* vend, __m128i* maskvp, uint32_t* ct_nmp, uint32_t* ct_hmajp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  const __m128i m16 = {0x0000ffff0000ffffLLU, 0x0000ffff0000ffffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i to_ct_nm1;
  __m128i to_ct_hmaj1;
  __m128i to_ct_nm2;
  __m128i to_ct_hmaj2;
  __uni16 acc_nm;
  __uni16 acc_hmaj;

  acc_nm.vi = _mm_setzero_si128();
  acc_hmaj.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_andnot_si128(loader2, loader3);
    to_ct_hmaj1 = _mm_and_si128(loader, loader3);

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_and_si128(loader, loader3));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_and_si128(loader, loader3));

    to_ct_nm1 = _mm_add_epi64(_mm_and_si128(to_ct_nm1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_nm1, 2), m2));
    to_ct_hmaj1 = _mm_add_epi64(_mm_and_si128(to_ct_hmaj1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj1, 2), m2));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_andnot_si128(loader2, loader3);
    to_ct_hmaj2 = _mm_and_si128(loader, loader3);

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_add_epi64(to_ct_nm2, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj2 = _mm_add_epi64(to_ct_hmaj2, _mm_and_si128(loader, loader3));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_add_epi64(to_ct_nm2, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj2 = _mm_add_epi64(to_ct_hmaj2, _mm_and_si128(loader, loader3));

    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_add_epi64(_mm_and_si128(to_ct_nm2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_nm2, 2), m2)));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_add_epi64(_mm_and_si128(to_ct_hmaj2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj2, 2), m2)));

    acc_nm.vi = _mm_add_epi64(acc_nm.vi, _mm_add_epi64(_mm_and_si128(to_ct_nm1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_nm1, 4), m4)));
    acc_hmaj.vi = _mm_add_epi64(acc_hmaj.vi, _mm_add_epi64(_mm_and_si128(to_ct_hmaj1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj1, 4), m4)));
  } while (vptr < vend);
  acc_nm.vi = _mm_add_epi64(_mm_and_si128(acc_nm.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_nm.vi, 8), m8));
  acc_hmaj.vi = _mm_add_epi64(_mm_and_si128(acc_hmaj.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_hmaj.vi, 8), m8));
  acc_nm.vi = _mm_and_si128(_mm_add_epi64(acc_nm.vi, _mm_srli_epi64(acc_nm.vi, 16)), m16);
  acc_hmaj.vi = _mm_and_si128(_mm_add_epi64(acc_hmaj.vi, _mm_srli_epi64(acc_hmaj.vi, 16)), m16);
  acc_nm.vi = _mm_add_epi64(acc_nm.vi, _mm_srli_epi64(acc_nm.vi, 32));
  acc_hmaj.vi = _mm_add_epi64(acc_hmaj.vi, _mm_srli_epi64(acc_hmaj.vi, 32));
  *ct_nmp += (uint32_t)(acc_nm.u8[0] + acc_nm.u8[1]);
  *ct_hmajp += (uint32_t)(acc_hmaj.u8[0] + acc_hmaj.u8[1]);
}
#else
void freq_hwe_haploid_count_12(uintptr_t* lptr, uintptr_t* maskp, uint32_t* ct_nmp, uint32_t* ct_hmajp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader3 = loader >> 1;
  uintptr_t loader2 = loader ^ (~loader3);
  uint32_t to_ct_nm1;
  uint32_t to_ct_hmaj1;
  uint32_t to_ct_nm2;
  uint32_t to_ct_hmaj2;
  uintptr_t partial_nm;
  uintptr_t partial_hmaj;
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 = loader2 & loader3;
  to_ct_hmaj1 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  to_ct_nm1 = (to_ct_nm1 & 0x33333333) + ((to_ct_nm1 >> 2) & 0x33333333);
  to_ct_nm1 += (to_ct_nm2 & 0x33333333) + ((to_ct_nm2 >> 2) & 0x33333333);
  partial_nm = (to_ct_nm1 & 0x0f0f0f0f) + ((to_ct_nm1 >> 4) & 0x0f0f0f0f);
  to_ct_hmaj1 = (to_ct_hmaj1 & 0x33333333) + ((to_ct_hmaj1 >> 2) & 0x33333333);
  to_ct_hmaj1 += (to_ct_hmaj2 & 0x33333333) + ((to_ct_hmaj2 >> 2) & 0x33333333);
  partial_hmaj = (to_ct_hmaj1 & 0x0f0f0f0f) + ((to_ct_hmaj1 >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 = loader2 & loader3;
  to_ct_hmaj1 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  to_ct_nm1 = (to_ct_nm1 & 0x33333333) + ((to_ct_nm1 >> 2) & 0x33333333);
  to_ct_nm1 += (to_ct_nm2 & 0x33333333) + ((to_ct_nm2 >> 2) & 0x33333333);
  partial_nm += (to_ct_nm1 & 0x0f0f0f0f) + ((to_ct_nm1 >> 4) & 0x0f0f0f0f);
  to_ct_hmaj1 = (to_ct_hmaj1 & 0x33333333) + ((to_ct_hmaj1 >> 2) & 0x33333333);
  to_ct_hmaj1 += (to_ct_hmaj2 & 0x33333333) + ((to_ct_hmaj2 >> 2) & 0x33333333);
  partial_hmaj += (to_ct_hmaj1 & 0x0f0f0f0f) + ((to_ct_hmaj1 >> 4) & 0x0f0f0f0f);

  *ct_nmp += (partial_nm * 0x01010101) >> 24;
  *ct_hmajp += (partial_hmaj * 0x01010101) >> 24;
}
#endif

static inline void single_marker_freqs_and_hwe(uintptr_t unfiltered_indiv_ct, uintptr_t unfiltered_indiv_ctl2, uintptr_t* lptr, uintptr_t* indiv_include2, uintptr_t* founder_include2, uintptr_t* founder_ctrl_include2, uintptr_t* founder_case_include2, uintptr_t indiv_ct, uint32_t* ll_ctp, uint32_t* lh_ctp, uint32_t* hh_ctp, uint32_t indiv_f_ct, uint32_t* ll_ctfp, uint32_t* lh_ctfp, uint32_t* hh_ctfp, int32_t hwe_needed, uintptr_t indiv_f_ctl_ct, uint32_t* ll_hwep, uint32_t* lh_hwep, uint32_t* hh_hwep, int32_t hardy_needed, uintptr_t indiv_f_case_ct, uint32_t* ll_case_hwep, uint32_t* lh_case_hwep, uint32_t* hh_case_hwep) {
  // This is best understood from the bottom third up (which is the order it
  // was written).  It's way overkill for just determining genotype
  // frequencies, but a ruthlessly optimized version is needed for e.g.
  // permutation testing so we may as well get it working here.
  //
  // The idea, which underlies the IBS and LD pruners as well, is to obtain
  // multiple popcounts at once.  In the case of PLINK frequency/HWE
  // evaluation, we need 6-9 numbers:
  // homozyg minor, het, homozyg major counts for all individuals
  // homozyg minor, het, homozyg major counts for all founders
  // sometimes, homozyg minor, het, homozyg major counts for all ctrl founders
  //
  // Given a buffer with PLINK binary genotypes for a single marker, let
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  // Then,
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  //               = indiv_ct - homozyg minor ct - het ct
  //
  // Thus, with the appropriate indiv_ct and these three popcounts, we can
  // calculate a set of genotype counts.  We perform the
  // not-that-exploitative version of these calculations in the bottom third of
  // this function, to deal with the remainder that doesn't fit into the
  // 12-word block size of the main loops.
  //
  // The middle third is a 12-word block popcount for 32-bit platforms (see
  // popcount_longs() in wdist_common.c; this is 12 words instead of 6 since
  // odd bits of the popcount targets are guaranteed to be zero, delaying
  // overflow).  It could be improved a bit, but we care more about reliability
  // than blistering speed for the 32-bit build (since it's used to check the
  // results of the actually blistering 64-bit code...).  Sure, hardware
  // popcount is significantly faster, but what x86-32 machine has hardware
  // popcount??...
  //
  // The top third is the portable Lauradoux/Walisch loop.
  uint32_t tot_a = 0;
  uint32_t tot_b = 0;
  uint32_t tot_c = 0;
  uint32_t tot_a_f = 0;
  uint32_t tot_b_f = 0;
  uint32_t tot_c_f = 0;
  uint32_t tot_a_hwe = 0;
  uint32_t tot_b_hwe = 0;
  uint32_t tot_c_hwe = 0;
  uint32_t tot_a_chwe = 0;
  uint32_t tot_b_chwe = 0;
  uint32_t tot_c_chwe = 0;
  uintptr_t* lptr_end = &(lptr[unfiltered_indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_12x_end;
  unfiltered_indiv_ctl2 -= unfiltered_indiv_ctl2 % 12;
  while (unfiltered_indiv_ctl2 >= 120) {
    cur_decr = 120;
  single_marker_freqs_and_hwe_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)indiv_include2, &tot_a, &tot_b, &tot_c);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_include2, &tot_a_f, &tot_b_f, &tot_c_f);
    if (hwe_needed) {
      count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_ctrl_include2, &tot_a_hwe, &tot_b_hwe, &tot_c_hwe);
      founder_ctrl_include2 = &(founder_ctrl_include2[cur_decr]);
      if (hardy_needed) {
	count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_case_include2, &tot_a_chwe, &tot_b_chwe, &tot_c_chwe);
	founder_case_include2 = &(founder_case_include2[cur_decr]);
      }
    }
    lptr = lptr_12x_end;
    indiv_include2 = &(indiv_include2[cur_decr]);
    founder_include2 = &(founder_include2[cur_decr]);
    unfiltered_indiv_ctl2 -= cur_decr;
  }
  if (unfiltered_indiv_ctl2) {
    cur_decr = unfiltered_indiv_ctl2;
    goto single_marker_freqs_and_hwe_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[unfiltered_indiv_ctl2 - unfiltered_indiv_ctl2 % 12]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, indiv_include2, &tot_a, &tot_b, &tot_c);
    count_3freq_12(lptr, founder_include2, &tot_a_f, &tot_b_f, &tot_c_f);
    if (hwe_needed) {
      count_3freq_12(lptr, founder_ctrl_include2, &tot_a_hwe, &tot_b_hwe, &tot_c_hwe);
      founder_ctrl_include2 = &(founder_ctrl_include2[12]);
      if (hardy_needed) {
	count_3freq_12(lptr, founder_case_include2, &tot_a_chwe, &tot_b_chwe, &tot_c_chwe);
	founder_case_include2 = &(founder_case_include2[12]);
      }
    }
    lptr = &(lptr[12]);
    indiv_include2 = &(indiv_include2[12]);
    founder_include2 = &(founder_include2[12]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = *indiv_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    // N.B. because of the construction of indiv_include2, only even-numbered
    // bits can be present here.  So popcount2_long is safe.
    tot_a += popcount2_long(loader2);
    tot_b += popcount2_long(loader3);
    tot_c += popcount2_long(loader & loader3);
    loader2 = *founder_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_a_f += popcount2_long(loader2);
    tot_b_f += popcount2_long(loader3);
    tot_c_f += popcount2_long(loader & loader3);
    if (hwe_needed) {
      loader2 = *founder_ctrl_include2++;
      loader3 = (loader >> 1) & loader2;
      loader2 &= loader;
      tot_a_hwe += popcount2_long(loader2);
      tot_b_hwe += popcount2_long(loader3);
      tot_c_hwe += popcount2_long(loader & loader3);
      if (hardy_needed) {
	loader2 = *founder_case_include2++;
	loader3 = (loader >> 1) & loader2;
	loader2 &= loader;
	tot_a_chwe += popcount2_long(loader2);
	tot_b_chwe += popcount2_long(loader3);
	tot_c_chwe += popcount2_long(loader & loader3);
      }
    }
  }
  *hh_ctp = tot_c;
  *lh_ctp = tot_b - tot_c;
  *ll_ctp = indiv_ct - tot_a - *lh_ctp;
  *hh_ctfp = tot_c_f;
  *lh_ctfp = tot_b_f - tot_c_f;
  *ll_ctfp = indiv_f_ct - tot_a_f - *lh_ctfp;
  if (hwe_needed) {
    *hh_hwep = tot_c_hwe;
    *lh_hwep = tot_b_hwe - tot_c_hwe;
    *ll_hwep = indiv_f_ctl_ct - tot_a_hwe - *lh_hwep;
    if (hardy_needed) {
      *hh_case_hwep = tot_c_chwe;
      *lh_case_hwep = tot_b_chwe - tot_c_chwe;
      *ll_case_hwep = indiv_f_case_ct - tot_a_chwe - *lh_case_hwep;
    }
  }
}

static inline void haploid_single_marker_freqs(uintptr_t unfiltered_indiv_ct, uintptr_t unfiltered_indiv_ctl2, uintptr_t* lptr, uintptr_t* indiv_include2, uintptr_t* founder_include2, uintptr_t indiv_ct, uint32_t* ll_ctp, uint32_t* hh_ctp, uint32_t indiv_f_ct, uint32_t* ll_ctfp, uint32_t* hh_ctfp, uint32_t* hethap_incr_ptr) {
  // Here, we interpret heterozygotes as missing.
  // Nonmissing: (genotype ^ (~(genotype >> 1))) & 0x5555...
  // Homozygote major: (genotype & (genotype >> 1)) & 0x5555...
  uint32_t tot_a = 0;
  uint32_t tot_b = 0;
  uint32_t tot_hmaj = 0;
  uint32_t tot_nm_f = 0;
  uint32_t tot_hmaj_f = 0;
  uintptr_t* lptr_end = &(lptr[unfiltered_indiv_ctl2]);
  uint32_t tot_nm;
  uint32_t hethap_incr;
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_12x_end;
  unfiltered_indiv_ctl2 -= unfiltered_indiv_ctl2 % 12;
  while (unfiltered_indiv_ctl2 >= 120) {
    cur_decr = 120;
  single_marker_freqs_and_hwe_loop:
    lptr_12x_end = &(lptr[cur_decr]);
  // Given a buffer with PLINK binary genotypes for a single marker, let
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  // Then,
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  //               = indiv_ct - homozyg minor ct - het ct
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)indiv_include2, &tot_a, &tot_b, &tot_hmaj);
    freq_hwe_haploid_count_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_include2, &tot_nm_f, &tot_hmaj_f);
    lptr = lptr_12x_end;
    indiv_include2 = &(indiv_include2[cur_decr]);
    founder_include2 = &(founder_include2[cur_decr]);
    unfiltered_indiv_ctl2 -= cur_decr;
  }
  if (unfiltered_indiv_ctl2) {
    cur_decr = unfiltered_indiv_ctl2;
    goto single_marker_freqs_and_hwe_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[unfiltered_indiv_ctl2 - unfiltered_indiv_ctl2 % 12]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, indiv_include2, &tot_a, &tot_b, &tot_hmaj);
    freq_hwe_haploid_count_12(lptr, founder_include2, &tot_nm_f, &tot_hmaj_f);
    lptr = &(lptr[12]);
    indiv_include2 = &(indiv_include2[12]);
    founder_include2 = &(founder_include2[12]);
  }
#endif
  tot_nm = 2 * tot_hmaj + (indiv_ct & (~(63 * ONELU))) - tot_a - tot_b;
  hethap_incr = tot_b - tot_hmaj;
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader3 = loader >> 1;
    loader4 = *indiv_include2++;
    hethap_incr += popcount2_long(loader3 & (~loader) & loader4);
    loader2 = loader ^ (~loader3);
    loader &= loader3;
    tot_nm += popcount2_long(loader2 & loader4);
    tot_hmaj += popcount2_long(loader & loader4);
    loader3 = *founder_include2++;
    tot_nm_f += popcount2_long(loader2 & loader3);
    tot_hmaj_f += popcount2_long(loader & loader3);
  }
  *hh_ctp = tot_hmaj;
  *ll_ctp = tot_nm - tot_hmaj;
  *hh_ctfp = tot_hmaj_f;
  *ll_ctfp = tot_nm_f - tot_hmaj_f;
  *hethap_incr_ptr = hethap_incr;
}

int32_t calc_freqs_and_hwe(FILE* bedfile, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_exclude_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t* founder_info, int32_t nonfounders, int32_t maf_succ, double* set_allele_freqs, uintptr_t** marker_reverse_ptr, uint32_t** marker_allele_cts_ptr, int32_t bed_offset, unsigned char missing_geno, int32_t hwe_needed, int32_t hwe_all, int32_t hardy_needed, uintptr_t* pheno_nm, uintptr_t* pheno_c, int32_t** hwe_lls_ptr, int32_t** hwe_lhs_ptr, int32_t** hwe_hhs_ptr, int32_t** hwe_ll_cases_ptr, int32_t** hwe_lh_cases_ptr, int32_t** hwe_hh_cases_ptr, int32_t** hwe_ll_allfs_ptr, int32_t** hwe_lh_allfs_ptr, int32_t** hwe_hh_allfs_ptr, int32_t** hwe_hapl_allfs_ptr, int32_t** hwe_haph_allfs_ptr, uint32_t* indiv_male_ct_ptr, uint32_t* indiv_f_ct_ptr, uint32_t* indiv_f_male_ct_ptr, int32_t wt_needed, unsigned char** marker_weights_base_ptr, double** marker_weights_ptr, double exponent, Chrom_info* chrom_info_ptr, uintptr_t* sex_nm, uintptr_t* sex_male, int32_t map_is_unsorted, int32_t* xmhh_exists_ptr, int32_t* nxmhh_exists_ptr) {
  FILE* hhfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + BITCT - 1) / BITCT;
  uintptr_t unfiltered_indiv_ctl2 = 2 * unfiltered_indiv_ctl;
  int32_t retval = 0;
  uint32_t pct = 1;
  uint32_t indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  uint32_t indiv_f_ct = indiv_ct;
  uintptr_t indiv_f_ctl_ct = indiv_ct;
  uintptr_t indiv_f_case_ct = indiv_ct;
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t ll_hwe = 0;
  uint32_t lh_hwe = 0;
  uint32_t hh_hwe = 0;
  uint32_t ll_case_hwe = 0;
  uint32_t lh_case_hwe = 0;
  uint32_t hh_case_hwe = 0;
  uint32_t cur_chrom_idx = 0;
  int32_t ii = chrom_info_ptr->chrom_file_order[0];
  uint32_t is_haploid = (species_haploid_mask[chrom_info_ptr->species] >> ii) & 1LLU;
  uint32_t next_chrom_start = chrom_info_ptr->chrom_file_order_marker_idx[1];
  uint32_t is_x = (ii == species_x_code[chrom_info_ptr->species]);
  uint32_t is_y = (ii == species_y_code[chrom_info_ptr->species]);
  uint32_t ll_ct = 0;
  uint32_t lh_ct = 0;
  uint32_t hh_ct = 0;
  uint32_t ll_ctf = 0;
  uint32_t lh_ctf = 0;
  uint32_t hh_ctf = 0;
  int32_t* hwe_lls = NULL;
  int32_t* hwe_lhs = NULL;
  int32_t* hwe_hhs = NULL;
  int32_t* hwe_ll_cases = NULL;
  int32_t* hwe_lh_cases = NULL;
  int32_t* hwe_hh_cases = NULL;
  int32_t* hwe_ll_allfs = NULL;
  int32_t* hwe_lh_allfs = NULL;
  int32_t* hwe_hh_allfs = NULL;
  uintptr_t* indiv_nonmale_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uintptr_t* founder_nonmale_include2 = NULL;
  uintptr_t* founder_ctrl_nonmale_include2 = NULL;
  uintptr_t* founder_male_include2 = NULL;
  uintptr_t* founder_case_include2 = NULL;
  uintptr_t* founder_case_nonmale_include2 = NULL;
  double* marker_weights = NULL;
  uint32_t indiv_nonmale_ct = 0;
  uint32_t indiv_f_nonmale_ct = 0;
  uint32_t indiv_f_ctl_nonmale_ct = 0;
  uint32_t indiv_f_case_nonmale_ct = 0;
  uint32_t indiv_male_ct = 0;
  uint32_t indiv_f_male_ct = 0;
  uint64_t hethap_ct = 0;
  uint32_t hethap_incr;
  int32_t* hwe_hapl_allfs;
  int32_t* hwe_haph_allfs;
  uintptr_t* loadbuf;
  uintptr_t* indiv_include2;
  uintptr_t* founder_include2;
  uintptr_t* founder_ctrl_include2;
  uint32_t nonmales_needed;
  uint32_t males_needed;
  uintptr_t* tmp_indiv_mask;
  uintptr_t* tmp_indiv_mask2;
  uintptr_t loop_end;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t* marker_allele_cts;
  double maf;
  if (wt_needed) {
    if (wkspace_alloc_uc_checked(marker_weights_base_ptr, CACHELINE) ||
        wkspace_alloc_d_checked(marker_weights_ptr, unfiltered_marker_ct * sizeof(double))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    marker_weights = *marker_weights_ptr;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      marker_weights[marker_uidx] = -1.0;
    }
  }
  uii = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(marker_reverse_ptr, uii * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  fill_ulong_zero(*marker_reverse_ptr, uii);

  if (!hwe_needed) {
    *hwe_lls_ptr = (int32_t*)wkspace_base;
  } else {
    if (wkspace_alloc_i_checked(&hwe_lls, unfiltered_marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_i_checked(&hwe_lhs, unfiltered_marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_i_checked(&hwe_hhs, unfiltered_marker_ct * sizeof(int32_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    *hwe_lls_ptr = hwe_lls;
    *hwe_lhs_ptr = hwe_lhs;
    *hwe_hhs_ptr = hwe_hhs;
    if (hardy_needed) {
      if (wkspace_alloc_i_checked(&hwe_ll_cases, unfiltered_marker_ct * sizeof(int32_t)) ||
          wkspace_alloc_i_checked(&hwe_lh_cases, unfiltered_marker_ct * sizeof(int32_t)) ||
          wkspace_alloc_i_checked(&hwe_hh_cases, unfiltered_marker_ct * sizeof(int32_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
    }
    *hwe_ll_cases_ptr = hwe_ll_cases;
    *hwe_lh_cases_ptr = hwe_lh_cases;
    *hwe_hh_cases_ptr = hwe_hh_cases;
  }
  if (wkspace_alloc_i_checked(&hwe_ll_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_lh_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_hh_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_hapl_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_haph_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&marker_allele_cts, 2 * unfiltered_marker_ct * sizeof(int32_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  *hwe_ll_allfs_ptr = hwe_ll_allfs;
  *hwe_lh_allfs_ptr = hwe_lh_allfs;
  *hwe_hh_allfs_ptr = hwe_hh_allfs;
  *hwe_hapl_allfs_ptr = hwe_hapl_allfs;
  *hwe_haph_allfs_ptr = hwe_haph_allfs;

  if ((!pheno_c) || map_is_unsorted) {
    hwe_all = 1;
  }

  fill_int_zero(hwe_ll_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_lh_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_hh_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_hapl_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_haph_allfs, unfiltered_marker_ct);
  fill_uint_zero(marker_allele_cts, 2 * unfiltered_marker_ct);
  *marker_allele_cts_ptr = marker_allele_cts;
  wkspace_mark = wkspace_base;
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&indiv_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  loadbuf[unfiltered_indiv_ctl2 - 2] = 0;
  loadbuf[unfiltered_indiv_ctl2 - 1] = 0;
  exclude_to_vec_include(unfiltered_indiv_ct, indiv_include2, indiv_exclude);
  ii = species_x_code[chrom_info_ptr->species];
  nonmales_needed = (!map_is_unsorted) && (ii != -1) && ((chrom_info_ptr->chrom_mask >> ii) & 1);
  ii = species_y_code[chrom_info_ptr->species];
  males_needed = nonmales_needed || ((!map_is_unsorted) && (ii != -1) && ((chrom_info_ptr->chrom_mask >> ii) & 1));
  if (wkspace_alloc_ul_checked(&indiv_male_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  memcpy(indiv_male_include2, indiv_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t));
  vec_include_mask_in(unfiltered_indiv_ct, indiv_male_include2, sex_nm);
  vec_include_mask_in(unfiltered_indiv_ct, indiv_male_include2, sex_male);
  indiv_male_ct = popcount_longs(indiv_male_include2, 0, unfiltered_indiv_ctl2);
  *indiv_male_ct_ptr = indiv_male_ct;
  if (males_needed) {
    founder_male_include2 = indiv_male_include2;
    if (nonmales_needed) {
      if (wkspace_alloc_ul_checked(&indiv_nonmale_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(indiv_nonmale_include2, indiv_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t));
      vec_include_mask_out_intersect(unfiltered_indiv_ct, indiv_nonmale_include2, sex_nm, sex_male);
      indiv_nonmale_ct = popcount_longs(indiv_nonmale_include2, 0, unfiltered_indiv_ctl2);
      founder_nonmale_include2 = indiv_nonmale_include2;
      founder_ctrl_nonmale_include2 = indiv_nonmale_include2;
    }
  }
  founder_include2 = indiv_include2;
  if (wkspace_alloc_ul_checked(&tmp_indiv_mask, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  if (!nonfounders) {
    if (wkspace_alloc_ul_checked(&founder_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    for (uii = 0; uii < unfiltered_indiv_ctl; uii++) {
      tmp_indiv_mask[uii] = indiv_exclude[uii] | (~founder_info[uii]);
    }
    zero_trailing_bits(tmp_indiv_mask, unfiltered_indiv_ct);
    exclude_to_vec_include(unfiltered_indiv_ct, founder_include2, tmp_indiv_mask);
    if (males_needed) {
      if (wkspace_alloc_ul_checked(&founder_male_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(founder_male_include2, indiv_male_include2, unfiltered_indiv_ctl * 2 * sizeof(intptr_t));
      vec_include_mask_in(unfiltered_indiv_ct, founder_male_include2, founder_info);
      indiv_f_male_ct = popcount_longs(founder_male_include2, 0, unfiltered_indiv_ctl2);
      if (nonmales_needed) {
	if (wkspace_alloc_ul_checked(&founder_nonmale_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
	  goto calc_freqs_and_hwe_ret_NOMEM;
	}
	memcpy(founder_nonmale_include2, indiv_nonmale_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t));
	vec_include_mask_in(unfiltered_indiv_ct, founder_nonmale_include2, founder_info);
	indiv_f_nonmale_ct = popcount_longs(founder_nonmale_include2, 0, unfiltered_indiv_ctl2);
      }
    }
    founder_ctrl_include2 = founder_include2;
    uii = popcount_longs_exclude(founder_info, indiv_exclude, unfiltered_indiv_ctl);
    indiv_f_ct = uii;
    indiv_f_ctl_ct = indiv_f_ct;
  } else {
    for (uii = 0; uii < unfiltered_indiv_ctl; uii++) {
      tmp_indiv_mask[uii] = indiv_exclude[uii];
    }
    founder_ctrl_include2 = founder_include2;
  }

  if (!hwe_all) {
    if (wkspace_alloc_ul_checked(&founder_ctrl_include2, unfiltered_indiv_ctl2 *  sizeof(intptr_t)) ||
	wkspace_alloc_ul_checked(&tmp_indiv_mask2, unfiltered_indiv_ctl * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    indiv_uidx = 0;
    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctl; indiv_uidx++) {
      tmp_indiv_mask2[indiv_uidx] = tmp_indiv_mask[indiv_uidx] | (~(pheno_nm[indiv_uidx])) | pheno_c[indiv_uidx];
    }
    zero_trailing_bits(tmp_indiv_mask2, unfiltered_indiv_ct);
    // tmp_indiv_mask is now set for each indiv who is excluded, or a
    // nonfounder, or is noncontrol.
    indiv_f_ctl_ct = unfiltered_indiv_ct - popcount_longs(tmp_indiv_mask2, 0, unfiltered_indiv_ctl);
    exclude_to_vec_include(unfiltered_indiv_ct, founder_ctrl_include2, tmp_indiv_mask2);
    if (nonmales_needed) {
      if (wkspace_alloc_ul_checked(&founder_ctrl_nonmale_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(founder_ctrl_nonmale_include2, indiv_nonmale_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t));
      vec_include_mask_out(unfiltered_indiv_ct, founder_ctrl_nonmale_include2, tmp_indiv_mask2);
      indiv_f_ctl_nonmale_ct = popcount_longs(founder_ctrl_nonmale_include2, 0, unfiltered_indiv_ctl2);
    }
    if (hardy_needed) {
      if (wkspace_alloc_ul_checked(&founder_case_include2, unfiltered_indiv_ctl2 *  sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      indiv_uidx = 0;
      for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctl; indiv_uidx++) {
	tmp_indiv_mask[indiv_uidx] |= (~(pheno_nm[indiv_uidx])) | (~pheno_c[indiv_uidx]);
      }
      zero_trailing_bits(tmp_indiv_mask, unfiltered_indiv_ct);
      indiv_f_case_ct = unfiltered_indiv_ct - popcount_longs(tmp_indiv_mask, 0, unfiltered_indiv_ctl);
      exclude_to_vec_include(unfiltered_indiv_ct, founder_case_include2, tmp_indiv_mask);
      if (nonmales_needed) {
	if (wkspace_alloc_ul_checked(&founder_case_nonmale_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
	  goto calc_freqs_and_hwe_ret_NOMEM;
	}
	memcpy(founder_case_nonmale_include2, indiv_nonmale_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t));
	vec_include_mask_out(unfiltered_indiv_ct, founder_case_nonmale_include2, tmp_indiv_mask);
	indiv_f_ctl_nonmale_ct = popcount_longs(founder_case_nonmale_include2, 0, unfiltered_indiv_ctl2);
      }
    }
  }

  *indiv_f_ct_ptr = indiv_f_ct;
  *indiv_f_male_ct_ptr = indiv_f_male_ct;
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto calc_freqs_and_hwe_ret_READ_FAIL;
  }
  marker_uidx = 0;
  marker_idx = 0;
  logprint("Calculating allele frequencies...");
  fputs(" 0%", stdout);
  fflush(stdout);
  if (map_is_unsorted) {
    is_haploid = species_haploid_mask[chrom_info_ptr->species] & 1;
    is_x = 0;
    is_y = 0;
    next_chrom_start = unfiltered_marker_ct;
  }
  for (; pct <= 100; pct++) {
    loop_end = ((uint64_t)pct * marker_ct) / 100LU;
    for (; marker_idx < loop_end; marker_idx++) {
      if (is_set(marker_exclude, marker_uidx)) {
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto calc_freqs_and_hwe_ret_READ_FAIL;
	}
      }
      if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto calc_freqs_and_hwe_ret_READ_FAIL;
      }
      if (marker_uidx >= next_chrom_start) {
	do {
	  cur_chrom_idx++;
	  next_chrom_start = chrom_info_ptr->chrom_file_order_marker_idx[cur_chrom_idx + 1];
	} while (marker_uidx >= next_chrom_start);
	ii = chrom_info_ptr->chrom_file_order[cur_chrom_idx];
	is_haploid = (species_haploid_mask[chrom_info_ptr->species] >> ii) & 1;
	is_x = (ii == species_x_code[chrom_info_ptr->species]);
	is_y = (ii == species_y_code[chrom_info_ptr->species]);
      }
      if (!is_haploid) {
	single_marker_freqs_and_hwe(unfiltered_indiv_ct, unfiltered_indiv_ctl2, loadbuf, indiv_include2, founder_include2, founder_ctrl_include2, founder_case_include2, indiv_ct, &ll_ct, &lh_ct, &hh_ct, indiv_f_ct, &ll_ctf, &lh_ctf, &hh_ctf, hwe_needed, indiv_f_ctl_ct, &ll_hwe, &lh_hwe, &hh_hwe, hardy_needed, indiv_f_case_ct, &ll_case_hwe, &lh_case_hwe, &hh_case_hwe);
	hwe_ll_allfs[marker_uidx] = ll_ctf;
	hwe_lh_allfs[marker_uidx] = lh_ctf;
	hwe_hh_allfs[marker_uidx] = hh_ctf;
	marker_allele_cts[2 * marker_uidx] += 2 * ll_ct + lh_ct;
	marker_allele_cts[2 * marker_uidx + 1] += 2 * hh_ct + lh_ct;
	uii = 2 * (ll_ctf + lh_ctf + hh_ctf + maf_succ);
	if (!uii) {
	  // avoid 0/0 division
	  set_allele_freqs[marker_uidx] = 0.5;
	} else {
	  set_allele_freqs[marker_uidx] = ((double)(2 * hh_ctf + lh_ctf + maf_succ)) / ((double)uii);
	}
	if (hwe_needed) {
	  hwe_lls[marker_uidx] = ll_hwe;
	  hwe_lhs[marker_uidx] = lh_hwe;
	  hwe_hhs[marker_uidx] = hh_hwe;
	  if (hardy_needed) {
	    hwe_ll_cases[marker_uidx] = ll_case_hwe;
	    hwe_lh_cases[marker_uidx] = lh_case_hwe;
	    hwe_hh_cases[marker_uidx] = hh_case_hwe;
	  }
	}
      } else {
	uii = 0;
	ujj = 0;
	if (is_x || is_y) {
	  if (is_x) {
	    single_marker_freqs_and_hwe(unfiltered_indiv_ct, unfiltered_indiv_ctl2, loadbuf, indiv_nonmale_include2, founder_nonmale_include2, founder_ctrl_nonmale_include2, founder_case_nonmale_include2, indiv_nonmale_ct, &ll_ct, &lh_ct, &hh_ct, indiv_f_nonmale_ct, &ll_ctf, &lh_ctf, &hh_ctf, hwe_needed, indiv_f_ctl_nonmale_ct, &ll_hwe, &lh_hwe, &hh_hwe, hardy_needed, indiv_f_case_nonmale_ct, &ll_case_hwe, &lh_case_hwe, &hh_case_hwe);
	    hwe_ll_allfs[marker_uidx] = ll_ctf;
	    hwe_lh_allfs[marker_uidx] = lh_ctf;
	    hwe_hh_allfs[marker_uidx] = hh_ctf;
	    marker_allele_cts[2 * marker_uidx] += 2 * ll_ct + lh_ct;
	    marker_allele_cts[2 * marker_uidx + 1] += 2 * hh_ct + lh_ct;
	    uii = 2 * (ll_ctf + lh_ctf + hh_ctf);
	    ujj = 2 * hh_ctf + lh_ctf;
	    if (hwe_needed) {
	      hwe_lls[marker_uidx] = ll_hwe;
	      hwe_lhs[marker_uidx] = lh_hwe;
	      hwe_hhs[marker_uidx] = hh_hwe;
	      if (hardy_needed) {
		hwe_ll_cases[marker_uidx] = ll_case_hwe;
		hwe_lh_cases[marker_uidx] = lh_case_hwe;
		hwe_hh_cases[marker_uidx] = hh_case_hwe;
	      }
	    }
	  }
	  haploid_single_marker_freqs(unfiltered_indiv_ct, unfiltered_indiv_ctl2, loadbuf, indiv_male_include2, founder_male_include2, indiv_male_ct, &ll_ct, &hh_ct, indiv_f_male_ct, &ll_ctf, &hh_ctf, &hethap_incr);
	} else {
	  haploid_single_marker_freqs(unfiltered_indiv_ct, unfiltered_indiv_ctl2, loadbuf, indiv_include2, founder_include2, indiv_ct, &ll_ct, &hh_ct, indiv_f_ct, &ll_ctf, &hh_ctf, &hethap_incr);
	}
	if (hethap_incr) {
	  if (!hhfile) {
	    memcpy(outname_end, ".hh", 4);
	    if (fopen_checked(&hhfile, outname, "w")) {
	      goto calc_freqs_and_hwe_ret_OPEN_FAIL;
	    }
	  }
	  if (is_x) {
	    *xmhh_exists_ptr = 1;
	  } else {
	    *nxmhh_exists_ptr = 1;
	  }
	  if (is_x || is_y) {
	    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctl2; indiv_uidx++) {
	      ulii = loadbuf[indiv_uidx];
	      ulii = (ulii >> 1) & (~ulii) & indiv_male_include2[indiv_uidx];
	      while (ulii) {
		ukk = indiv_uidx * BITCT2 + CTZLU(ulii) / 2;
		if (fprintf(hhfile, "%s\t%s\n", &(person_ids[ukk * max_person_id_len]), &(marker_ids[marker_uidx * max_marker_id_len])) < 0) {
		  goto calc_freqs_and_hwe_ret_WRITE_FAIL;
		}
		ulii &= ulii - ONELU;
	      }
	    }
	  } else {
	    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctl2; indiv_uidx++) {
	      ulii = loadbuf[indiv_uidx];
	      ulii = (ulii >> 1) & (~ulii) & indiv_include2[indiv_uidx];
	      while (ulii) {
		ukk = indiv_uidx * BITCT2 + CTZLU(ulii) / 2;
		if (fprintf(hhfile, "%s\t%s\n", &(person_ids[ukk * max_person_id_len]), &(marker_ids[marker_uidx * max_marker_id_len])) < 0) {
		  goto calc_freqs_and_hwe_ret_WRITE_FAIL;
		}
		ulii &= ulii - ONELU;
	      }
	    }
	  }
	  hethap_ct += hethap_incr;
	}
	hwe_hapl_allfs[marker_uidx] = ll_ctf;
	hwe_haph_allfs[marker_uidx] = hh_ctf;
	marker_allele_cts[2 * marker_uidx] += ll_ct;
	marker_allele_cts[2 * marker_uidx + 1] += hh_ct;
	uii += ll_ctf + hh_ctf + 2 * maf_succ;
	ujj += hh_ctf + maf_succ;
	if (!uii) {
	  maf = 0.5;
	} else {
	  maf = ((double)ujj) / ((double)uii);
	}
	set_allele_freqs[marker_uidx] = maf;
	if (wt_needed) {
	  marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, maf);
	}
      }
      marker_uidx++;
    }
    if (pct < 100) {
      if (pct > 10) {
	putchar('\b');
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
  }
  fputs("\b\b\b\b", stdout);
  logprint(" done.\n");
  if (hethap_ct) {
    *outname_end = '\0';
    sprintf(logbuf, "%" PRIu64 " heterozygous haploid genotype%s set to missing (see %s.hh).\n", hethap_ct, (hethap_ct == 1LLU)? "" : "s", outname);
    logprintb();
  }
  while (0) {
  calc_freqs_and_hwe_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_freqs_and_hwe_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_freqs_and_hwe_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  calc_freqs_and_hwe_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(hhfile);
  return retval;
}

int32_t load_one_freq_mal1(char allele_min, char allele_maj, double maf, double* set_allele_freq_ptr, char* marker_alleles, char missing_geno) {
  char cc;
  if (maf > 0.5) {
    cc = allele_maj;
    allele_maj = allele_min;
    allele_min = cc;
    maf = 1.0 - maf;
  }
  if (marker_alleles[0] == allele_min) {
    if (marker_alleles[1] == allele_maj) {
      *set_allele_freq_ptr = 1.0 - maf;
    } else {
      return -1;
    }
  } else if (marker_alleles[1] == allele_min) {
    if (marker_alleles[0] == allele_maj) {
      *set_allele_freq_ptr = maf;
    } else {
      return -1;
    }
  } else if ((allele_min == missing_geno) && (maf == 0.0)) {
    if (marker_alleles[0] == allele_maj) {
      *set_allele_freq_ptr = 0.0;
    } else if (marker_alleles[1] == allele_maj) {
      *set_allele_freq_ptr = 1.0;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
  return 0;
}

// aptr1 = minor, aptr2 = major
int32_t load_one_freq_malx(uint32_t alen1, char* aptr1, uint32_t alen2, char* aptr2, double maf, double* set_allele_freq_ptr, char* maptr1, uintptr_t max_marker_allele_len, char missing_geno) {
  uint32_t malen1 = strlen(maptr1);
  char* maptr2 = &(maptr1[max_marker_allele_len]);
  uint32_t malen2 = strlen(maptr2);
  uint32_t uii;
  char* aptr;
  if (maf > 0.5) {
    aptr = aptr2;
    uii = alen2;
    aptr2 = aptr1;
    alen2 = alen1;
    aptr1 = aptr;
    alen1 = uii;
    maf = 1.0 - maf;
  }
  if ((malen1 == alen1) && (!memcmp(maptr1, aptr1, alen1))) {
    if ((malen2 == alen2) && (!memcmp(maptr2, aptr2, alen2))) {
      *set_allele_freq_ptr = 1.0 - maf;
    } else {
      return -1;
    }
  } else if ((malen2 == alen1) && (!memcmp(maptr2, aptr1, alen1))) {
    if ((malen1 == alen2) && (!memcmp(maptr1, aptr2, alen2))) {
      *set_allele_freq_ptr = maf;
    } else {
      return -1;
    }
  } else if ((*aptr1 == missing_geno) && (alen1 == 1) && (maf == 0.0)) {
    if ((malen1 == alen2) && (!memcmp(maptr1, aptr2, alen2))) {
      *set_allele_freq_ptr = 0.0;
    } else if ((malen2 == alen2) && (!memcmp(maptr2, aptr2, alen2))) {
      *set_allele_freq_ptr = 1.0;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
  return 0;
}

int32_t read_external_freqs(char* freqname, FILE** freqfile_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, char* marker_alleles, uintptr_t max_marker_allele_len, uint32_t* marker_allele_cts, double* set_allele_freqs, int32_t maf_succ, char missing_geno, double exponent, int32_t wt_needed, double* marker_weights) {
  uint32_t species = chrom_info_ptr->species;
  int32_t freq_counts = 0;
  int32_t duplicate_fail = 1;
  char cc = '0';
  char cc2 = '0';
  uint32_t alen1 = 0;
  uint32_t alen2 = 0;
  char* aptr1 = NULL;
  char* aptr2 = NULL;
  char missing_geno_buf[2];
  unsigned char* wkspace_mark;
  char* sorted_ids;
  int32_t* id_map;
  int32_t ii;
  int32_t jj;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* bufptr5;
  double maf;
  int32_t c_hom_a1;
  int32_t c_het;
  int32_t c_hom_a2;
  int32_t c_hap_a1;
  int32_t c_hap_a2;
  if (fopen_checked(freqfile_ptr, freqname, "r")) {
    return RET_OPEN_FAIL;
  }
  if (fgets(tbuf, MAXLINELEN, *freqfile_ptr) == NULL) {
    logprint("Error: Empty --read-freq file.\n");
    return RET_INVALID_FORMAT;
  }
  wkspace_mark = wkspace_base;
  ii = sort_item_ids(&sorted_ids, &id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, strcmp_deref, &duplicate_fail);
  if (ii) {
    return ii;
  }
  if (!memcmp(tbuf, " CHR  ", 6)) {
    ii = strlen(tbuf);
    if (tbuf[ii - 2] == '0') { // --counts makes G0 the last column header
      freq_counts = 1;
    } else if (tbuf[ii - 2] != 'S') { // NCHROBS
      goto read_external_freqs_ret_INVALID_FORMAT;
    }
    while (fgets(tbuf, MAXLINELEN, *freqfile_ptr) != NULL) {
      bufptr = skip_initial_spaces(tbuf);
      jj = marker_code(species, bufptr);
      bufptr = next_item(bufptr); // now at beginning of marker name
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      read_next_terminate(bufptr, bufptr); // destructive read (\0 at end of item)
      ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
        if ((jj == get_marker_chrom(chrom_info_ptr, ii)) || (!jj) || (!get_marker_chrom(chrom_info_ptr, ii))) {
	  if (max_marker_allele_len == 1) {
	    cc = *bufptr2;
	    bufptr2 = next_item(bufptr2);
	    if (cc == *bufptr2) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	  } else {
	    alen1 = strlen_se(bufptr2);
	    aptr1 = bufptr2;
	    bufptr2 = next_item(bufptr2);
	    if (no_more_items_kns(bufptr2)) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	    alen2 = strlen_se(bufptr2);
	    aptr2 = bufptr2;
	    if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	  }
	  bufptr = next_item(bufptr2);
	  if (no_more_items_kns(bufptr)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  if (freq_counts) {
	    if (no_more_items_kns(next_item(bufptr))) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	    c_hom_a1 = atoi(bufptr);
	    c_hom_a2 = atoi(next_item(bufptr));
	    maf = ((double)c_hom_a1 + maf_succ) / ((double)(c_hom_a1 + c_hom_a2 + 2 * maf_succ));
	  } else {
	    if (sscanf(bufptr, "%lg", &maf) != 1) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	  }
	  if (max_marker_allele_len == 1) {
	    if (load_one_freq_mal1(cc, *bufptr2, maf, &(set_allele_freqs[ii]), &(marker_alleles[ii * 2]), missing_geno)) {
	      goto read_external_freqs_ret_ALLELE_MISMATCH;
	    }
	  } else {
	    if (load_one_freq_malx(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[ii]), &(marker_alleles[ii * max_marker_allele_len * 2]), max_marker_allele_len, missing_geno)) {
	      goto read_external_freqs_ret_ALLELE_MISMATCH;
	    }
	  }
	  if (wt_needed) {
	    marker_weights[ii] = calc_wt_mean_maf(exponent, set_allele_freqs[ii]);
	  }
        }
      }
    }
    if (freq_counts) {
      fputs(".frq.count file loaded.\n", stdout);
    } else {
      fputs(".frq file loaded.\n", stdout);
    }
  } else if (!memcmp(tbuf, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)", 71)) {
    // changed from strcmp to avoid eoln problems
    // known --freqx format, v0.15.3 or later
    while (fgets(tbuf, MAXLINELEN, *freqfile_ptr) != NULL) {
      jj = marker_code(species, tbuf);
      bufptr = next_item(tbuf); // now at beginning of marker name
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      read_next_terminate(bufptr, bufptr); // destructive read (\0 at end of item)
      ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
        if ((jj == get_marker_chrom(chrom_info_ptr, ii)) || (!jj) || (!get_marker_chrom(chrom_info_ptr, ii))) {
	  if (max_marker_allele_len == 1) {
	    cc = *bufptr2;
	    bufptr2 = next_item(bufptr2);
	    cc2 = *bufptr2;
	    if (cc == cc2) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	  } else {
	    alen1 = strlen_se(bufptr2);
	    aptr1 = bufptr2;
	    bufptr2 = next_item(bufptr2);
	    if (no_more_items_kns(bufptr2)) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	    alen2 = strlen_se(bufptr2);
	    aptr2 = bufptr2;
	    if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	  }
	  bufptr = next_item(bufptr2);
	  bufptr2 = next_item(bufptr);
	  bufptr3 = next_item(bufptr2);
	  bufptr4 = next_item(bufptr3);
	  bufptr5 = next_item(bufptr4);
	  if (no_more_items_kns(bufptr5)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  c_hom_a1 = atoi(bufptr);
	  c_het = atoi(bufptr2);
	  c_hom_a2 = atoi(bufptr3);
	  c_hap_a1 = atoi(bufptr4);
	  c_hap_a2 = atoi(bufptr5);
	  maf = ((double)(c_hom_a1 * 2 + c_het + c_hap_a1 + maf_succ)) / ((double)(2 * (c_hom_a1 + c_het + c_hom_a2 + maf_succ) + c_hap_a1 + c_hap_a2));
	  if (max_marker_allele_len == 1) {
	    if (load_one_freq_mal1(cc, cc2, maf, &(set_allele_freqs[ii]), &(marker_alleles[ii * 2]), missing_geno)) {
	      goto read_external_freqs_ret_ALLELE_MISMATCH;
	    }
	  } else {
	    if (load_one_freq_malx(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[ii]), &(marker_alleles[ii * max_marker_allele_len * 2]), max_marker_allele_len, missing_geno)) {
	      goto read_external_freqs_ret_ALLELE_MISMATCH;
	    }
	  }
	  if (wt_needed) {
	    if (c_hap_a1 || c_hap_a2) {
	      marker_weights[ii] = calc_wt_mean_maf(exponent, set_allele_freqs[ii]);
	    } else {
	      marker_weights[ii] = calc_wt_mean(exponent, c_het, c_hom_a1, c_hom_a2);
	    }
	  }
        }
      }
    }
    fputs(".frqx file loaded.\n", stdout);
  } else {
    // Also support GCTA-style frequency files:
    // [marker ID]\t[reference allele]\t[frequency of reference allele]\n
    missing_geno_buf[0] = missing_geno;
    missing_geno_buf[1] = '\0';
    do {
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      bufptr = next_item(bufptr);
      if (!bufptr) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      read_next_terminate(tbuf, tbuf); // destructive read
      ii = bsearch_str(tbuf, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
	if (max_marker_allele_len == 1) {
          cc = *bufptr;
	} else {
	  alen1 = strlen_se(bufptr);
	  aptr1 = bufptr;
	}
        bufptr = next_item(bufptr);
	if (no_more_items_kns(bufptr)) {
          goto read_external_freqs_ret_INVALID_FORMAT;
	}
        if (sscanf(bufptr, "%lg", &maf) != 1) {
          goto read_external_freqs_ret_INVALID_FORMAT;
        }
	if (max_marker_allele_len == 1) {
	  if (load_one_freq_mal1(missing_geno, cc, maf, &(set_allele_freqs[ii]), &(marker_alleles[ii * 2]), missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH2;
	  }
	} else {
	  if (load_one_freq_malx(1, missing_geno_buf, alen1, aptr1, maf, &(set_allele_freqs[ii]), &(marker_alleles[ii * max_marker_allele_len * 2]), max_marker_allele_len, missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
	  }
	}
	if (wt_needed) {
	  marker_weights[ii] = calc_wt_mean_maf(exponent, set_allele_freqs[ii]);
	}
      } else {
	// if there aren't exactly 3 columns, this isn't a GCTA .freq file
	bufptr = next_item(bufptr);
	if (no_more_items_kns(bufptr) || (!no_more_items_kns(next_item(bufptr)))) {
	  goto read_external_freqs_ret_INVALID_FORMAT;
	}
      }
    } while (fgets(tbuf, MAXLINELEN, *freqfile_ptr) != NULL);
    fputs("GCTA-formatted .freq file loaded.\n", stdout);
  }
  fclose_null(freqfile_ptr);
  wkspace_reset(wkspace_mark);
  return 0;
 read_external_freqs_ret_INVALID_FORMAT:
  logprint(errstr_freq_format);
  return RET_INVALID_FORMAT;
 read_external_freqs_ret_ALLELE_MISMATCH:
  sprintf(logbuf, "Error: Mismatch between .bim/.ped and --freq alleles at %s.\n", next_item_mult(tbuf, 2));
  logprintb();
  return RET_ALLELE_MISMATCH;
 read_external_freqs_ret_ALLELE_MISMATCH2:
  sprintf(logbuf, "Error: Mismatch between .bim/.ped and --freq alleles at %s.\n", tbuf);
  logprintb();
  return RET_ALLELE_MISMATCH;
}

int32_t write_freqs(FILE** outfile_ptr, char* outname, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, double* set_allele_freqs, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, int32_t* hapl_cts, int32_t* haph_cts, uint32_t indiv_f_ct, uint32_t indiv_f_male_ct, int32_t freq_counts, int32_t freqx, char missing_geno, uintptr_t* marker_reverse) {
  int32_t reverse = 0;
  char minor_allele = '\0';
  char missing_geno_buf[2];
  uint32_t uii;
  char cc2;
  char* minor_ptr;
  char* major_ptr;
  int32_t chrom_idx;
  uint32_t chrom_end_marker;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t missing_ct;
  if (fopen_checked(outfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  if (freqx) {
    if (fputs("CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)\n", *outfile_ptr) == EOF) {
      return RET_WRITE_FAIL;
    }
  } else if (plink_maxsnp < 5) {
    if (freq_counts) {
      if (fputs_checked(" CHR  SNP   A1   A2     C1     C2     G0\n", *outfile_ptr)) {
	return RET_WRITE_FAIL;
      }
      if (max_marker_allele_len == 1) {
        strcpy(tbuf, "%4d %4s    %c    %c %6u %6u %6u\n");
      } else {
	strcpy(tbuf, "%4d %4s %4s %4s %6u %6u %6u\n");
      }
    } else {
      if (fputs_checked(" CHR  SNP   A1   A2          MAF  NCHROBS\n", *outfile_ptr)) {
        return RET_WRITE_FAIL;
      }
      if (max_marker_allele_len == 1) {
        strcpy(tbuf, "%4d %4s    %c    %c %12.4g %8d\n");
      } else {
        strcpy(tbuf, "%4d %4s %4s %4s %12.4g %8d\n");
      }
    }
  } else if (freq_counts) {
    sprintf(tbuf, " CHR %%%ds   A1   A2     C1     C2     G0\n", plink_maxsnp);
    if (fprintf(*outfile_ptr, tbuf, "SNP") < 0) {
      return RET_WRITE_FAIL;
    }
    if (max_marker_allele_len == 1) {
      sprintf(tbuf, "%%4d %%%ds    %%c    %%c %%6u %%6u %%6u\n", plink_maxsnp);
    } else {
      sprintf(tbuf, "%%4d %%%ds %%4s %%4s %%6u %%6u %%6u\n", plink_maxsnp);
    }
  } else {
    sprintf(tbuf, " CHR %%%ds   A1   A2          MAF  NCHROBS\n", plink_maxsnp);
    if (fprintf(*outfile_ptr, tbuf, "SNP") < 0) {
      return RET_WRITE_FAIL;
    }
    if (max_marker_allele_len == 1) {
      sprintf(tbuf, "%%4d %%%ds    %%c    %%c %%12.4g %%8d\n", plink_maxsnp);
    } else {
      sprintf(tbuf, "%%4d %%%ds %%4s %%4s %%12.4g %%8d\n", plink_maxsnp);
    }
  }
  missing_geno_buf[0] = missing_geno;
  missing_geno_buf[1] = '\0';
  for (chrom_idx = 0; chrom_idx < MAX_POSSIBLE_CHROM; chrom_idx++) {
    if (!chrom_exists(chrom_info_ptr, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == species_x_code[chrom_info_ptr->species]);
    is_y = (chrom_idx == species_y_code[chrom_info_ptr->species]);
    is_haploid = (species_haploid_mask[chrom_info_ptr->species] >> chrom_idx) & 1LLU;
    chrom_end_marker = chrom_info_ptr->chrom_end[chrom_idx];
    if (max_marker_allele_len == 1) {
      uii = next_non_set(marker_exclude, chrom_info_ptr->chrom_start[chrom_idx], chrom_end_marker);
      while (uii < chrom_end_marker) {
	if (!freq_counts) {
	  reverse = is_set(marker_reverse, uii);
	  minor_allele = marker_alleles[uii * 2 + reverse];
	  if (!minor_allele) {
	    minor_allele = missing_geno;
	  }
	}
	if (freq_counts || freqx) {
	  if (is_x) {
	    missing_ct = indiv_f_ct - (ll_cts[uii] + lh_cts[uii] + hh_cts[uii] + hapl_cts[uii] + haph_cts[uii]);
	  } else if (is_haploid) {
	    if (is_y) {
	      missing_ct = indiv_f_male_ct - (hapl_cts[uii] + haph_cts[uii]);
	    } else {
	      missing_ct = indiv_f_ct - (hapl_cts[uii] + haph_cts[uii]);
	    }
	  } else {
	    missing_ct = indiv_f_ct - (ll_cts[uii] + lh_cts[uii] + hh_cts[uii]);
	  }
	  if (freq_counts) {
	    minor_allele = marker_alleles[uii * 2];
	    if (!minor_allele) {
	      minor_allele = missing_geno;
	    }
	    cc2 = marker_alleles[uii * 2 + 1];
	    if (!cc2) {
	      cc2 = missing_geno;
	    }
	    if (fprintf(*outfile_ptr, tbuf, get_marker_chrom(chrom_info_ptr, uii), &(marker_ids[uii * max_marker_id_len]), minor_allele, cc2, 2 * ll_cts[uii] + lh_cts[uii] + hapl_cts[uii], 2 * hh_cts[uii] + lh_cts[uii] + haph_cts[uii], missing_ct) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  } else {
	    if (fprintf(*outfile_ptr, "%d\t%s\t%c\t%c\t%u\t%u\t%u\t%u\t%u\t%u\n", get_marker_chrom(chrom_info_ptr, uii), &(marker_ids[uii * max_marker_id_len]), minor_allele, marker_alleles[uii * 2 + (1 ^ reverse)], reverse? hh_cts[uii] : ll_cts[uii], lh_cts[uii], reverse? ll_cts[uii] : hh_cts[uii], reverse? haph_cts[uii] : hapl_cts[uii], reverse? hapl_cts[uii] : haph_cts[uii], missing_ct) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  }
	} else {
	  if (fprintf(*outfile_ptr, tbuf, get_marker_chrom(chrom_info_ptr, uii), &(marker_ids[uii * max_marker_id_len]), minor_allele, marker_alleles[uii * 2 + (1 ^ reverse)], reverse? set_allele_freqs[uii] : (1.0 - set_allele_freqs[uii]), 2 * (ll_cts[uii] + lh_cts[uii] + hh_cts[uii]) + hapl_cts[uii] + haph_cts[uii]) < 0) {
	    return RET_WRITE_FAIL;
	  }
	}
	uii = next_non_set(marker_exclude, uii + 1, chrom_end_marker);
      }
    } else {
      uii = next_non_set(marker_exclude, chrom_info_ptr->chrom_start[chrom_idx], chrom_end_marker);
      while (uii < chrom_end_marker) {
	if (freqx || (!freq_counts)) {
	  reverse = is_set(marker_reverse, uii);
	  major_ptr = &(marker_alleles[(uii * 2 + (1 ^ reverse)) * max_marker_allele_len]);
	} else {
	  major_ptr = &(marker_alleles[(uii * 2 + 1) * max_marker_allele_len]);
	  if (!(*major_ptr)) {
	    major_ptr = missing_geno_buf;
	  }
	}
	minor_ptr = &(marker_alleles[(uii * 2 + reverse) * max_marker_allele_len]);
	if (!(*minor_ptr)) {
	  minor_ptr = missing_geno_buf;
	}
	if (freq_counts || freqx) {
	  if (is_x) {
	    missing_ct = indiv_f_ct - (ll_cts[uii] + lh_cts[uii] + hh_cts[uii] + hapl_cts[uii] + haph_cts[uii]);
	  } else if (is_haploid) {
	    if (is_y) {
	      missing_ct = indiv_f_male_ct - (hapl_cts[uii] + haph_cts[uii]);
	    } else {
	      missing_ct = indiv_f_ct - (hapl_cts[uii] + haph_cts[uii]);
	    }
	  } else {
	    missing_ct = indiv_f_ct - (ll_cts[uii] + lh_cts[uii] + hh_cts[uii]);
	  }
	  if (freqx) {
	    if (fprintf(*outfile_ptr, "%d\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\n", get_marker_chrom(chrom_info_ptr, uii), &(marker_ids[uii * max_marker_id_len]), minor_ptr, major_ptr, reverse? hh_cts[uii] : ll_cts[uii], lh_cts[uii], reverse? ll_cts[uii] : hh_cts[uii], reverse? haph_cts[uii] : hapl_cts[uii], reverse? hapl_cts[uii] : haph_cts[uii], missing_ct) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  } else {
	    if (fprintf(*outfile_ptr, tbuf, get_marker_chrom(chrom_info_ptr, uii), &(marker_ids[uii * max_marker_id_len]), minor_ptr, major_ptr, 2 * ll_cts[uii] + lh_cts[uii] + hapl_cts[uii], 2 * hh_cts[uii] + lh_cts[uii] + haph_cts[uii], missing_ct) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  }
	} else {
	  if (fprintf(*outfile_ptr, tbuf, get_marker_chrom(chrom_info_ptr, uii), &(marker_ids[uii * max_marker_id_len]), minor_ptr, major_ptr, reverse? set_allele_freqs[uii] : (1.0 - set_allele_freqs[uii]), 2 * (ll_cts[uii] + lh_cts[uii] + hh_cts[uii]) + hapl_cts[uii] + haph_cts[uii]) < 0) {
	    return RET_WRITE_FAIL;
	  }
	}
        uii = next_non_set(marker_exclude, uii + 1, chrom_end_marker);
      }
    }
  }
  if (fclose_null(outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  sprintf(logbuf, "Allele frequencies written to %s.\n", outname);
  logprintb();
  return 0;
}

uint32_t binary_geno_filter(double geno_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uintptr_t indiv_ct, uintptr_t male_ct, uint32_t* marker_allele_cts, Chrom_info* chrom_info_ptr) {
  uint32_t orig_exclude_ct = *marker_exclude_ct_ptr;
  uint32_t geno_int_thresh = 2 * indiv_ct - (int)(geno_thresh * 2 * indiv_ct);
  uintptr_t marker_uidx = 0;
  uint32_t species = chrom_info_ptr->species;
  uint32_t chrom_end;
  uint32_t chrom_fo_idx;
  int32_t chrom_idx;
  uint32_t cur_ct;
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    if ((species_haploid_mask[species] >> chrom_idx) & 1LLU) {
      if (chrom_idx == species_x_code[species]) {
	cur_ct = 2 * indiv_ct - male_ct;
      } else if (chrom_idx == species_y_code[species]) {
	cur_ct = male_ct;
      } else {
	cur_ct = indiv_ct;
      }
    } else {
      cur_ct = 2 * indiv_ct;
    }
    geno_int_thresh = cur_ct - (int)(geno_thresh * cur_ct);
    marker_uidx = next_non_set(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    while (marker_uidx < chrom_end) {
      if ((marker_allele_cts[2 * marker_uidx] + marker_allele_cts[2 * marker_uidx + 1]) < geno_int_thresh) {
	set_bit(marker_exclude, marker_uidx, marker_exclude_ct_ptr);
      }
      marker_uidx = next_non_set(marker_exclude, marker_uidx + 1, chrom_end);
    }
  }
  return (*marker_exclude_ct_ptr - orig_exclude_ct);
}

void calc_marker_reverse_bin(uintptr_t* marker_reverse, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, double* set_allele_freqs) {
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  for (; marker_idx < marker_ct; marker_idx++) {
    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
    if (set_allele_freqs[marker_uidx] < 0.5) {
      set_bit_noct(marker_reverse, marker_uidx);
    }
    marker_uidx++;
  }
}

int32_t hardy_report_write_line(FILE* outfile, char* prefix_buf, uint32_t prefix_len, uint32_t reverse, uint32_t ll_ct, uint32_t lh_ct, uint32_t hh_ct, char* midbuf_ptr, double pvalue) {
  uint32_t uii;
  uint32_t ujj;
  char* midline;
  char* cptr;
  uint32_t denom;
  double drecip;
  double minor_freq;
  fwrite(prefix_buf, prefix_len, 1, outfile);
  if (reverse) {
    uii = sprintf(midbuf_ptr, "%u/%u/%u ", hh_ct, lh_ct, ll_ct);
  } else {
    uii = sprintf(midbuf_ptr, "%u/%u/%u ", ll_ct, lh_ct, hh_ct);
  }
  if (uii < 21) {
    ujj = 21 - uii;
    midline = midbuf_ptr - ujj;
    memset(midline, 32, ujj);
  } else {
    midline = midbuf_ptr;
  }
  cptr = &(midbuf_ptr[uii]);
  denom = (ll_ct + lh_ct + hh_ct) * 2;
  if (denom) {
    drecip = 1.0 / ((double)denom);
    minor_freq = (2 * ll_ct + lh_ct) * drecip;
    cptr += sprintf(cptr, "%8.4g %8.4g %12.4g\n", (lh_ct * 2) * drecip, minor_freq * (2 * hh_ct + lh_ct) * drecip * 2, pvalue);
  } else {
    cptr += sprintf(cptr, "     nan      nan           NA\n");
  }
  return fwrite_checked(midline, (cptr - midline), outfile);
}

int32_t hardy_report(pthread_t* threads, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, int32_t* hwe_ll_cases, int32_t* hwe_lh_cases, int32_t* hwe_hh_cases, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs, uint32_t pheno_nm_ct, uintptr_t* pheno_c, Chrom_info* chrom_info_ptr) {
  FILE* outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  int32_t retval = 0;
  uint32_t pct = 0;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  uint32_t ukk = 9;
  uint32_t loop_end;
  uint32_t uii;
  uint32_t ujj;
  uint32_t report_type;
  uint32_t is_x;
  uint32_t is_haploid;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t reverse;
  double* p_values;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;

  if (pheno_nm_ct) {
    report_type = pheno_c? 0 : 1;
  } else {
    report_type = 2;
  }
  uii = report_type? 1 : 3;
  if (wkspace_alloc_d_checked(&p_values, uii * marker_ct * sizeof(double))) {
    goto hardy_report_ret_NOMEM;
  }

  // todo: multithread?
  if (report_type) {
    for (; marker_idx < marker_ct; marker_idx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      p_values[marker_uidx] = SNPHWE2(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx]);
      marker_uidx++;
    }
  } else {
    for (; marker_idx < marker_ct; marker_idx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      p_values[marker_uidx * 3] = SNPHWE2(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx]);
      p_values[marker_uidx * 3 + 1] = SNPHWE2(hwe_lh_cases[marker_uidx], hwe_ll_cases[marker_uidx], hwe_hh_cases[marker_uidx]);
      p_values[marker_uidx * 3 + 2] = SNPHWE2(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx]);
      marker_uidx++;
    }
  }
  marker_uidx = 0;
  marker_idx = 0;

  memcpy(outname_end, ".hwe", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto hardy_report_ret_OPEN_FAIL;
  }
  sprintf(logbuf, "Writing Hardy-Weinberg report to %s...", outname);
  logprintb();
  fputs(" 0%", stdout);
  sprintf(tbuf, " CHR %%%us     TEST   A1   A2                 GENO   O(HET)   E(HET)            P \n", plink_maxsnp);
  if (fprintf(outfile, tbuf, "SNP") < 0) {
    goto hardy_report_ret_WRITE_FAIL;
  }
 
  cptr = &(tbuf[15 + plink_maxsnp]);
  tbuf[0] = ' ';
  tbuf[1] = ' ';
  chrom_fo_idx = 0;
  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
  intprint2(&(tbuf[2]), chrom_info_ptr->chrom_file_order[chrom_fo_idx]);
  if (report_type) {
    if (report_type == 1) {
      memcpy(&(tbuf[5 + plink_maxsnp]), "  ALL(QT)           ", 20);
    } else {
      memcpy(&(tbuf[5 + plink_maxsnp]), "  ALL(NP)           ", 20);
    }
    ujj = 16 + plink_maxsnp;
    cptr2 = &(tbuf[ujj + 22 + 2 * max_marker_allele_len]);
    for (; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100LLU;
      for (; marker_idx < loop_end; marker_idx++) {
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	  intprint2(&(tbuf[2]), chrom_info_ptr->chrom_file_order[chrom_fo_idx]);
	}
	cptr3 = &(marker_ids[marker_uidx * max_marker_id_len]);
	uii = strlen(cptr3);
	memset(&(tbuf[5]), 32, plink_maxsnp - uii);
	memcpy(&(tbuf[5 + plink_maxsnp - uii]), cptr3, uii);
	reverse = is_set(marker_reverse, marker_uidx);
	if (reverse) {
	  cptr3 = &(marker_alleles[(2 * marker_uidx + 1) * max_marker_allele_len]);
	  cptr4 = &(marker_alleles[2 * marker_uidx * max_marker_allele_len]);
	} else {
	  cptr3 = &(marker_alleles[2 * marker_uidx * max_marker_allele_len]);
	  cptr4 = &(marker_alleles[(2 * marker_uidx + 1) * max_marker_allele_len]);
	}
	if (max_marker_allele_len == 1) {
	  cptr[3] = *cptr3;
	  cptr[8] = *cptr4;
	} else {
	  uii = strlen(cptr3);
	  if (uii < 5) {
	    memset(cptr, 32, 4 - uii);
	    memcpy(&(cptr[4 - uii]), cptr3, uii);
	    ukk = 4;
	  } else {
	    memcpy(cptr, cptr3, uii);
	    ukk = uii;
	  }
	  uii = strlen(cptr4);
	  if (uii < 5) {
	    memset(&(cptr[ukk]), 32, 5 - uii);
	    memcpy(&(cptr[ukk + 5 - uii]), cptr4, uii);
	    ukk += 5;
	  } else {
	    cptr[ukk] = ' ';
	    memcpy(&(cptr[ukk + 1]), cptr4, uii);
	    ukk += uii + 1;
	  }
	}
	if (hardy_report_write_line(outfile, tbuf, ujj + ukk, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[marker_uidx])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}
	if (ukk > 9) {
	  memset(&(cptr[9]), 32, ukk - 9);
	}
	marker_uidx++;
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else {
    memcpy(&(tbuf[5 + plink_maxsnp]), "                    ", 20);
    ujj = 16 + plink_maxsnp;
    cptr2 = &(tbuf[ujj + 22 + 2 * max_marker_allele_len]);
    for (; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100LLU;
      for (; marker_idx < loop_end; marker_idx++) {
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	  intprint2(&(tbuf[2]), chrom_info_ptr->chrom_file_order[chrom_fo_idx]);
	}
	cptr3 = &(marker_ids[marker_uidx * max_marker_id_len]);
	uii = strlen(cptr3);
	memset(&(tbuf[5]), 32, plink_maxsnp - uii);
	memcpy(&(tbuf[5 + plink_maxsnp - uii]), cptr3, uii);
	memcpy(&(tbuf[9 + plink_maxsnp]), "  ALL", 5);
	reverse = is_set(marker_reverse, marker_uidx);
	if (reverse) {
	  cptr3 = &(marker_alleles[(2 * marker_uidx + 1) * max_marker_allele_len]);
	  cptr4 = &(marker_alleles[2 * marker_uidx * max_marker_allele_len]);
	} else {
	  cptr3 = &(marker_alleles[2 * marker_uidx * max_marker_allele_len]);
	  cptr4 = &(marker_alleles[(2 * marker_uidx + 1) * max_marker_allele_len]);
	}
	if (max_marker_allele_len == 1) {
	  cptr[3] = *cptr3;
	  cptr[8] = *cptr4;
	} else {
	  uii = strlen(cptr3);
	  if (uii < 5) {
	    memset(cptr, 32, 4 - uii);
	    memcpy(&(cptr[4 - uii]), cptr3, uii);
	    ukk = 4;
	  } else {
	    memcpy(cptr, cptr3, uii);
	    ukk = uii;
	  }
	  uii = strlen(cptr4);
	  if (uii < 5) {
	    memset(&(cptr[ukk]), 32, 5 - uii);
	    memcpy(&(cptr[ukk + 5 - uii]), cptr4, uii);
	    ukk += 5;
	  } else {
	    cptr[ukk] = ' ';
	    memcpy(&(cptr[ukk + 1]), cptr4, uii);
	    ukk += uii + 1;
	  }
	}
	if (hardy_report_write_line(outfile, tbuf, ujj + ukk, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[3 * marker_uidx])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(tbuf[11 + plink_maxsnp]), "AFF", 3);
	if (hardy_report_write_line(outfile, tbuf, ujj + ukk, reverse, hwe_ll_cases[marker_uidx], hwe_lh_cases[marker_uidx], hwe_hh_cases[marker_uidx], cptr2, p_values[3 * marker_uidx + 1])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(tbuf[9 + plink_maxsnp]), "UN", 2);
	if (hardy_report_write_line(outfile, tbuf, ujj + ukk, reverse, hwe_lls[marker_uidx], hwe_lhs[marker_uidx], hwe_hhs[marker_uidx], cptr2, p_values[3 * marker_uidx + 2])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}
	if (ukk > 9) {
	  memset(&(cptr[9]), 32, ukk - 9);
	}
	marker_uidx++;
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  }
  fputs("\b\b\b\b", stdout);
  logprint(" done.\n");

  while (0) {
  hardy_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  hardy_report_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  hardy_report_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}


void enforce_hwe_threshold(double hwe_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, int32_t hwe_all, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs) {
  uint32_t removed_ct = 0;
  uintptr_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  hwe_thresh += EPSILON;
  if (!hwe_all) {
    for (; marker_idx < marker_ct; marker_idx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      if (SNPHWE_t(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx], hwe_thresh)) {
	set_bit_noct(marker_exclude, marker_uidx);
	removed_ct++;
      }
      marker_uidx++;
    }
  } else {
    for (; marker_idx < marker_ct; marker_idx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      if (SNPHWE_t(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], hwe_thresh)) {
	set_bit_noct(marker_exclude, marker_uidx);
	removed_ct++;
      }
      marker_uidx++;
    }
  }
  sprintf(logbuf, "%u marker%s removed due to Hardy-Weinberg exact test (--hwe).\n", removed_ct, (removed_ct == 1)? "" : "s");
  logprintb();
  *marker_exclude_ct_ptr += removed_ct;
}

void enforce_maf_threshold(double min_maf, double max_maf, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, double* set_allele_freqs) {
  uint32_t removed_ct = 0;
  uintptr_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  double dxx;
  for (; marker_idx < marker_ct; marker_idx++) {
    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
    dxx = get_maf(set_allele_freqs[marker_uidx]);
    if ((dxx < min_maf - EPSILON) || (dxx > max_maf + EPSILON)) {
      set_bit(marker_exclude, marker_uidx, marker_exclude_ct_ptr);
      removed_ct++;
    }
    marker_uidx++;
  }
  sprintf(logbuf, "%u marker%s removed due to MAF threshold(s) (--maf/--max-maf).\n", removed_ct, (removed_ct == 1)? "" : "s");
  logprintb();
}

void calc_marker_weights(double exponent, uintptr_t* marker_exclude, uintptr_t marker_ct, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, double* marker_weights) {
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  for (; marker_idx < marker_ct; marker_idx++) {
    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
    if (marker_weights[marker_uidx] < 0.0) {
      marker_weights[marker_uidx] = calc_wt_mean(exponent, lh_cts[marker_uidx], ll_cts[marker_uidx], hh_cts[marker_uidx]);
    }
    marker_uidx++;
  }
}

int32_t load_ref_alleles(char* infilename, uintptr_t unfiltered_marker_ct, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  uint32_t slen;
  char* sorted_ids;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  int32_t* id_map;
  int32_t retval;
  int32_t sorted_idx;
  char cc;
  if (fopen_checked(&infile, infilename, "r")) {
    goto load_ref_alleles_ret_OPEN_FAIL;
  }
  retval = sort_item_ids_nx(&sorted_ids, &id_map, unfiltered_marker_ct, marker_ids, max_marker_id_len);
  if (retval) {
    goto load_ref_alleles_ret_1;
  }
  while (fgets(tbuf, MAXLINELEN, infile)) {
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr2 = item_endnn(bufptr);
    bufptr3 = skip_initial_spaces(bufptr2);
    if (is_eoln_kns(*bufptr3)) {
      logprint("Error: Improperly formatted --reference-allele file.\n");
      goto load_ref_alleles_ret_INVALID_FORMAT;
    }
    *bufptr2 = '\0';
    sorted_idx = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - 1);
    if (sorted_idx == -1) {
      continue;
    } else {
      sorted_idx = id_map[sorted_idx];
      if (max_marker_allele_len == 1) {
	cc = *bufptr3;
	if (cc == marker_alleles[sorted_idx * 2]) {
	  clear_bit_noct(marker_reverse, sorted_idx);
	} else if (cc == marker_alleles[sorted_idx * 2 + 1]) {
	  set_bit_noct(marker_reverse, sorted_idx);
	} else {
	  sprintf(logbuf, "Warning: Impossible reference allele assignment for marker %s.\n", bufptr);
	  logprintb();
	}
      } else {
	slen = strlen_se(bufptr3);
	if (slen >= max_marker_allele_len) {
	  sprintf(logbuf, "Warning: Impossible reference allele assignment for marker %s.\n", bufptr);
	  logprintb();
	} else {
	  bufptr3[slen++] = '\0';
	  if (!memcmp(bufptr3, &(marker_alleles[sorted_idx * 2 * max_marker_allele_len]), slen)) {
	    clear_bit_noct(marker_reverse, sorted_idx);
	  } else if (!memcmp(bufptr3, &(marker_alleles[(sorted_idx * 2 + 1) * max_marker_allele_len]), slen)) {
	    set_bit_noct(marker_reverse, sorted_idx);
	  } else {
	    sprintf(logbuf, "Warning: Impossible reference allele assignment for marker %s.\n", bufptr);
	    logprintb();
	  }
	}
      }
    }
  }
  if (!feof(infile)) {
    goto load_ref_alleles_ret_READ_FAIL;
  }
  while (0) {
  load_ref_alleles_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_ref_alleles_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_ref_alleles_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 load_ref_alleles_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t write_snplist(FILE** outfile_ptr, char* outname, char* outname_end, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len) {
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  strcpy(outname_end, ".snplist");
  if (fopen_checked(outfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  for (; marker_idx < marker_ct; marker_idx++) {
    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
    if (fprintf(*outfile_ptr, "%s\n", &(marker_ids[marker_uidx * max_marker_id_len])) < 0) {
      return RET_WRITE_FAIL;
    }
    marker_uidx++;
  }
  if (fclose_null(outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  sprintf(logbuf, "Final list of marker IDs written to %s.\n", outname);
  logprintb();
  return 0;
}

inline int32_t relationship_or_ibc_req(uint64_t calculation_type) {
  return (relationship_req(calculation_type) || (calculation_type & CALC_IBC));
}

inline int32_t distance_wt_req(uint64_t calculation_type) {
  return ((calculation_type & CALC_DISTANCE) || ((!(calculation_type & CALC_LOAD_DISTANCES)) && ((calculation_type & CALC_GROUPDIST) || (calculation_type & CALC_REGRESS_DISTANCE))));
}

int32_t wdist(char* outname, char* outname_end, char* pedname, char* mapname, char* famname, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* filtername, char* freqname, char* loaddistname, char* evecname, char* mergename1, char* mergename2, char* mergename3, char* makepheno_str, char* phenoname_str, char* refalleles, char* recode_allele_name, char* covar_fname, char* cluster_fname, char* set_fname, char* subset_fname, char* update_alleles_fname, char* update_map_fname, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* loop_assoc_fname, char* filterval, int32_t mfilter_col, int32_t filter_case_control, int32_t filter_sex, int32_t filter_founder_nonf, int32_t fam_col_1, int32_t fam_col_34, int32_t fam_col_5, int32_t fam_col_6, char missing_geno, int32_t missing_pheno, char output_missing_geno, char* output_missing_pheno, int32_t mpheno_col, uint32_t pheno_modifier, int32_t prune, int32_t affection_01, Chrom_info* chrom_info_ptr, double exponent, double min_maf, double max_maf, double geno_thresh, double mind_thresh, double hwe_thresh, int32_t hwe_all, double rel_cutoff, int32_t tail_pheno, double tail_bottom, double tail_top, uint64_t calculation_type, int32_t rel_calc_type, int32_t dist_calc_type, uintptr_t groupdist_iters, int32_t groupdist_d, uintptr_t regress_iters, int32_t regress_d, uintptr_t regress_rel_iters, int32_t regress_rel_d, double unrelated_herit_tol, double unrelated_herit_covg, double unrelated_herit_covr, int32_t ibc_type, int32_t parallel_idx, uint32_t parallel_tot, int32_t ppc_gap, int32_t allow_no_sex, int32_t must_have_sex, int32_t nonfounders, int32_t genome_output_gz, int32_t genome_output_full, int32_t genome_ibd_unbounded, int32_t ld_window_size, int32_t ld_window_kb, int32_t ld_window_incr, double ld_last_param, int32_t maf_succ, int32_t regress_pcs_normalize_pheno, int32_t regress_pcs_sex_specific, int32_t regress_pcs_clip, int32_t max_pcs, int32_t freq_counts, int32_t freqx, int32_t distance_flat_missing, uint32_t recode_modifier, int32_t allelexxxx, int32_t merge_type, int32_t indiv_sort, int32_t keep_allele_order, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, char* snps_flag_markers, unsigned char* snps_flag_starts_range, uint32_t snps_flag_ct, uint32_t snps_flag_max_len, uint32_t set_hh_missing, uint32_t covar_modifier, char* covar_str, uint32_t mcovar_col, uint32_t update_map_modifier, uint32_t write_covar_modifier, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uint32_t gxe_mcovar, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope) {
  FILE* outfile = NULL;
  FILE* outfile2 = NULL;
  FILE* outfile3 = NULL;
  FILE* bimfile = NULL;
  FILE* bedfile = NULL;
  FILE* famfile = NULL;
  FILE* phenofile = NULL;
  FILE* filterfile = NULL;
  FILE* bedtmpfile = NULL;
  FILE* bimtmpfile = NULL;
  FILE* famtmpfile = NULL;
  FILE* freqfile = NULL;
  FILE* loaddistfile = NULL;
  char* id_buf = NULL;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t* marker_exclude = NULL;
  uintptr_t max_marker_id_len = 0;
  // set_allele_freqs = frequency of allele corresponding to set bits in .bed
  //   (i.e. A2), or frequency of MAJOR allele in middle of text loading.
  double* set_allele_freqs = NULL;
  uintptr_t unfiltered_indiv_ct = 0;
  uintptr_t unfiltered_indiv_ct4 = 0;
  uintptr_t unfiltered_indiv_ctl = 0;
  uintptr_t* indiv_exclude = NULL;
  uintptr_t indiv_exclude_ct = 0;
  uint32_t* indiv_sort_map = NULL;
  uintptr_t* founder_info = NULL;
  uintptr_t* sex_nm = NULL;
  uintptr_t* sex_male = NULL;
  uintptr_t marker_ct;
  uint32_t plink_maxsnp;
  int32_t ii;
  int32_t jj = 0;
  int32_t kk = 0;
  uint32_t uii = 0;
  uintptr_t ulii;
  uintptr_t uljj;
  double dxx;
  int64_t llxx = 0;
  int64_t llyy;
  char* marker_ids = NULL;
  unsigned char* marker_weights_base = NULL;
  // if max_marker_allele_len == 1:
  //   marker_alleles[2 * ii] is id of A1 (usually minor) allele at marker ii
  //   marker_alleles[2 * ii + 1] is identity of A2 allele at marker ii
  // otherwise:
  //   max_marker_allele_len INCLUDES TRAILING NULL, so at least 3
  //   marker_alleles[2 * ii * max_marker_allele_len] is id of A1
  //   marker_alleles[(2 * ii + 1) + max_marker_allele_len] is id of A2
  char* marker_alleles = NULL;
  uintptr_t max_marker_allele_len = 1;
  uintptr_t* marker_reverse = NULL;
  uint32_t* marker_allele_cts;
  int32_t retval = RET_SUCCESS;
  int32_t map_is_unsorted = 0;
  int32_t map_cols = 3;
  int32_t affection = 0;
  uintptr_t* pheno_nm = NULL;
  uintptr_t* pheno_c = NULL;
  double* pheno_d = NULL;
  double* phenor_d = NULL;
  char* phenor_c = NULL;
  double* marker_weights = NULL;
  uint32_t* marker_weights_i = NULL;
  char* person_ids = NULL;
  uintptr_t max_person_id_len;
  char* paternal_ids = NULL;
  uintptr_t max_paternal_id_len = 2;
  char* maternal_ids = NULL;
  uintptr_t max_maternal_id_len = 2;
  unsigned char* wkspace_mark = NULL;
  char* cptr = NULL;
  int32_t* iptr;
  uint64_t dists_alloc = 0;
  double* rel_ibc;
  uintptr_t marker_exclude_ct = 0;
  char* pid_list = NULL;
  char* id_list = NULL;
  double missing_phenod = (double)missing_pheno;
  double ci_zt = 0.0;
  int32_t missing_pheno_len = intlen(missing_pheno);
  int32_t var_std = 1;
  int32_t* hwe_lls;
  int32_t* hwe_lhs;
  int32_t* hwe_hhs;
  int32_t* hwe_ll_cases;
  int32_t* hwe_lh_cases;
  int32_t* hwe_hh_cases;
  int32_t* hwe_ll_allfs;
  int32_t* hwe_lh_allfs;
  int32_t* hwe_hh_allfs;
  int32_t* hwe_hapl_allfs;
  int32_t* hwe_haph_allfs;
  uint32_t indiv_male_ct;
  uint32_t indiv_f_ct;
  uint32_t indiv_f_male_ct;
  pthread_t threads[MAX_THREADS];
  int32_t exp0 = (exponent == 0.0);
  int32_t wt_needed = 0;
  int32_t hwe_needed = 0;
  int32_t bed_offset = 3;
  uint32_t* marker_pos = NULL;
  int32_t xmhh_exists = 0;
  int32_t nxmhh_exists = 0;
  uint32_t pheno_ctrl_ct = 0;
  uint32_t pheno_nm_ct;
  Pedigree_rel_info pri;
  unsigned char* wkspace_mark2 = NULL;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  int32_t duplicate_fail;

  if (ci_size != 0.0) {
    ci_zt = ltqnorm(1 - (1 - ci_size) / 2);
  }
  if (rel_calc_type & REL_CALC_COV) {
    var_std = 0;
    ibc_type = -1;
  }

  if (calculation_type & CALC_MAKE_BED) {
#if _WIN32
    uii = GetFullPathName(pedname, FNAMESIZE, tbuf, NULL);
    if ((!uii) || (uii > FNAMESIZE))
#else
    if (!realpath(pedname, tbuf))
#endif
    {
      sprintf(logbuf, "Error: Failed to open %s.\n", pedname);
      logprintb();
      goto wdist_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".bed", 5);
    // if file doesn't exist, realpath returns NULL on Linux instead of what
    // the path would be.
#if _WIN32
    uii = GetFullPathName(outname, FNAMESIZE, &(tbuf[FNAMESIZE + 64]), NULL);
    if (uii && (uii <= FNAMESIZE) && (!strcmp(tbuf, &(tbuf[FNAMESIZE + 64]))))
#else
    cptr = realpath(outname, &(tbuf[FNAMESIZE + 64]));
    if (cptr && (!strcmp(tbuf, &(tbuf[FNAMESIZE + 64]))))
#endif
    {
      logprint("Note: --make-bed input and output filenames match.  Appending .old to input\nfilenames.\n");
      uii = strlen(pedname);
      memcpy(tbuf, pedname, uii + 1);
      memcpy(&(pedname[uii]), ".old", 5);
      if (rename(tbuf, pedname)) {
	logprint("Error: Failed to append .old to input .bed filename.\n");
	goto wdist_ret_OPEN_FAIL;
      }
      uii = strlen(mapname);
      memcpy(tbuf, mapname, uii + 1);
      memcpy(&(mapname[uii]), ".old", 5);
      if (rename(tbuf, mapname)) {
	logprint("Error: Failed to append .old to input .bim filename.\n");
	goto wdist_ret_OPEN_FAIL;
      }
      uii = strlen(famname);
      memcpy(tbuf, famname, uii + 1);
      memcpy(&(famname[uii]), ".old", 5);
      if (rename(tbuf, famname)) {
	logprint("Error: Failed to append .old to input .fam filename.\n");
	goto wdist_ret_OPEN_FAIL;
      }
    }
  }

  if (calculation_type & CALC_MERGE) {
    if (!(fam_col_1 && fam_col_34 && fam_col_5 && fam_col_6 && (!affection_01) && (missing_pheno == -9))) {
      logprint("Error: --merge/--bmerge/--merge-list cannot be used with an irregularly\nformatted reference fileset (--no-fid, --no-parents, --no-sex, --no-pheno,\n--1).  Use --make-bed first.\n");
      goto wdist_ret_INVALID_CMDLINE;
    }
    // Only append -merge to the filename stem if --make-bed or --recode lgen
    // is specified.
    ulii = bed_suffix_conflict(calculation_type, recode_modifier);
    if (ulii) {
      memcpy(outname_end, "-merge", 7);
    }
    retval = merge_datasets(pedname, mapname, famname, outname, ulii? &(outname_end[6]) : outname_end, mergename1, mergename2, mergename3, calculation_type, merge_type, indiv_sort, keep_allele_order, chrom_info_ptr);
    if (retval || (!(calculation_type & (~CALC_MERGE)))) {
      goto wdist_ret_2;
    }
    uljj = (uintptr_t)(outname_end - outname) + (ulii? 6 : 0);
    memcpy(pedname, outname, uljj);
    memcpy(&(pedname[uljj]), ".bed", 5);
    memcpy(famname, pedname, uljj);
    memcpy(&(famname[uljj]), ".fam", 5);
    memcpy(mapname, pedname, uljj);
    memcpy(&(mapname[uljj]), ".bim", 5);
  }

  if (fopen_checked(&bedfile, pedname, "rb")) {
    goto wdist_ret_OPEN_FAIL;
  }
  if (fopen_checked(&famfile, famname, "rb")) {
    goto wdist_ret_OPEN_FAIL;
  }
  // load .bim, count markers, filter chromosomes
  retval = load_bim(&bimfile, mapname, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &plink_maxsnp, &marker_exclude, &set_allele_freqs, &marker_alleles, &max_marker_allele_len, &marker_ids, chrom_info_ptr, &marker_pos, extractname, excludename, freqname, refalleles, calculation_type, recode_modifier, allelexxxx, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, snps_flag_markers, snps_flag_starts_range, snps_flag_ct, snps_flag_max_len, &map_is_unsorted);
  if (retval) {
    goto wdist_ret_2;
  }

  ulii = MAXLINELEN;
  // load .fam, count indivs
  ii = fam_col_6;
  if (ii && phenoname) {
    ii = (pheno_modifier & PHENO_MERGE) && (!makepheno_str);
  }
  retval = load_fam(famfile, ulii, fam_col_1, fam_col_34, fam_col_5, ii, fam_col_6, missing_pheno, missing_pheno_len, affection_01, &unfiltered_indiv_ct, &person_ids, &max_person_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_nm, &sex_male, &affection, &pheno_nm, &pheno_c, &pheno_d, &founder_info, &indiv_exclude);
  if (retval) {
    goto wdist_ret_2;
  }
  count_genders(sex_nm, sex_male, unfiltered_indiv_ct, indiv_exclude, &ii, &jj, &kk);
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if (kk) {
    sprintf(logbuf, "%" PRIuPTR " marker%s and %" PRIuPTR " %s (%d male%s, %d female%s, %d unknown) loaded.\n", marker_ct, (marker_ct == 1)? "" : "s", unfiltered_indiv_ct, species_str(unfiltered_indiv_ct), ii, (ii == 1)? "" : "s", jj, (jj == 1)? "" : "s", kk);
  } else {
    sprintf(logbuf, "%" PRIuPTR " marker%s and %" PRIuPTR " %s (%d male%s, %d female%s) loaded.\n", marker_ct, (marker_ct == 1)? "" : "s", unfiltered_indiv_ct, species_str(unfiltered_indiv_ct), ii, (ii == 1)? "" : "s", jj, (jj == 1)? "" : "s");
  }
  logprintb();

  unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;

  if (phenoname && fopen_checked(&phenofile, phenoname, "r")) {
    goto wdist_ret_OPEN_FAIL;
  }

  if (phenofile || tail_pheno) {
    wkspace_mark = wkspace_base;
    duplicate_fail = 1;
    retval = sort_item_ids(&cptr, &iptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, strcmp_deref, &duplicate_fail);
    if (retval) {
      goto wdist_ret_2;
    }

    if (makepheno_str) {
      retval = makepheno_load(phenofile, makepheno_str, unfiltered_indiv_ct, indiv_exclude_ct, cptr, max_person_id_len, iptr, pheno_nm, &pheno_c);
      if (retval) {
	goto wdist_ret_2;
      }
    } else if (phenofile) {
      retval = load_pheno(phenofile, unfiltered_indiv_ct, indiv_exclude_ct, cptr, max_person_id_len, iptr, missing_pheno, missing_pheno_len, affection_01, mpheno_col, phenoname_str, pheno_nm, &pheno_c, &pheno_d);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    if (tail_pheno) {
      retval = convert_tail_pheno(unfiltered_indiv_ct, pheno_nm, &pheno_c, &pheno_d, tail_bottom, tail_top, missing_phenod);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    wkspace_reset(wkspace_mark);
  }

  if ((calculation_type & CALC_GROUPDIST) && (!pheno_c)) {
    logprint("Error: --groupdist calculation requires dichotomous phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_REGRESS_DISTANCE) && (!pheno_d)) {
    logprint("Error: --regress-distance calculation requires scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_UNRELATED_HERITABILITY) && (!pheno_d)) {
    logprint("Error: --unrelated-heritability requires scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & (CALC_REGRESS_PCS | CALC_REGRESS_PCS_DISTANCE)) && (!pheno_d)) {
    sprintf(logbuf, "Error: --regress-pcs%s requires scalar phenotype.\n", (calculation_type & CALC_REGRESS_PCS_DISTANCE)? "-distance" : "");
    goto wdist_ret_INVALID_CMDLINE_2;
  } else if ((calculation_type & CALC_MODEL) && (!(model_modifier & MODEL_ASSOC)) && (!pheno_c)) {
    logprint("Error: --model requires dichotomous phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  }

  if (prune) {
    prune_missing_phenos(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, pheno_nm, pheno_d, missing_phenod);
  }

  if (allelexxxx) {
    allelexxxx_recode(allelexxxx, marker_alleles, max_marker_allele_len, marker_exclude, unfiltered_marker_ct - marker_exclude_ct);
  }

  if (extractname || excludename) {
    wkspace_mark = wkspace_base;
    duplicate_fail = 0;
    retval = sort_item_ids(&cptr, &iptr, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, strcmp_deref, &duplicate_fail);
    if (retval) {
      goto wdist_ret_2;
    }
    if (duplicate_fail) {
      logprint("Warning: Duplicate marker ID(s) in input.\n");
    }
    // length of sorted list is NOT necessarily equal to unfiltered_marker_ct -
    // marker_exclude_ct, since marker_exclude_ct may change before second call
    ii = unfiltered_marker_ct - marker_exclude_ct;
    if (extractname) {
      retval = include_or_exclude(extractname, cptr, ii, max_marker_id_len, iptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0, 0, duplicate_fail);
      if (retval) {
        goto wdist_ret_1;
      }
    }
    if (excludename) {
      retval = include_or_exclude(excludename, cptr, ii, max_marker_id_len, iptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0, 1, duplicate_fail);
      if (retval) {
        goto wdist_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }

  if (removename || keepname || filtername || (indiv_sort == INDIV_SORT_ASCII)) {
    wkspace_mark = wkspace_base;
    duplicate_fail = 0;
    retval = sort_item_ids(&cptr, &iptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, strcmp_deref, &duplicate_fail);
    if (retval) {
      goto wdist_ret_2;
    }
    if (duplicate_fail) {
      logprint("Warning: Duplicate individual ID(s) in input.\n");
    }
    ii = unfiltered_indiv_ct - indiv_exclude_ct;
    if (removename) {
      retval = include_or_exclude(removename, cptr, ii, max_person_id_len, iptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 1, 1, duplicate_fail);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    if (keepname) {
      retval = include_or_exclude(keepname, cptr, ii, max_person_id_len, iptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 1, 0, duplicate_fail);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    if (filtername) {
      if (!mfilter_col) {
	mfilter_col = 1;
      }
      retval = filter_indivs_file(filtername, cptr, ii, max_person_id_len, iptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, filterval, mfilter_col);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    if (indiv_sort == INDIV_SORT_ASCII) {
      indiv_sort_map = (uint32_t*)iptr;
      wkspace_reset((unsigned char*)cptr);
    } else {
      wkspace_reset(wkspace_mark);
    }
  }
  if (indiv_sort == INDIV_SORT_NATURAL) {
    duplicate_fail = 0;
    retval = sort_item_ids(&cptr, &iptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, strcmp_natural_deref, &duplicate_fail);
    if (retval) {
      goto wdist_ret_2;
    }
    indiv_sort_map = (uint32_t*)iptr;
    wkspace_reset((unsigned char*)cptr);
  }

  if (filter_case_control) {
    if (!pheno_c) {
      logprint("Error: --filter-cases/--filter-controls requires dichotomous phenotype.\n");
      goto wdist_ret_INVALID_CMDLINE;
    }
    ii = indiv_exclude_ct;
    // fcc == 1: exclude all zeroes in pheno_c
    // fcc == 2: exclude all ones in pheno_c
    // -> flip on fcc == 1
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, pheno_c, 2 - filter_case_control, pheno_nm);
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d %s removed due to case/control status (--filter-%s).\n", ii, species_str(ii), (filter_case_control == 1)? "cases" : "controls");
    logprintb();
  }
  if (filter_sex) {
    ii = indiv_exclude_ct;
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, sex_male, 2 - filter_sex, sex_nm);
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d %s removed due to gender filter (--filter-%s).\n", ii, species_str(ii), (filter_sex == 1)? "males" : "females");
    logprintb();
  }
  if (filter_founder_nonf) {
    ii = indiv_exclude_ct;
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, founder_info, 2 - filter_founder_nonf, NULL);
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d %s removed due to founder status (--filter-%s).\n", ii, species_str(ii), (filter_founder_nonf == 1)? "founders" : "nonfounders");
    logprintb();
  }

  if (fseeko(bedfile, 0, SEEK_END)) {
    goto wdist_ret_READ_FAIL;
  }
  llxx = ftello(bedfile);
  llyy = llxx - ((uint64_t)unfiltered_indiv_ct4) * unfiltered_marker_ct;
  rewind(bedfile);
  if (llyy == 3LL) {
    // v1.00 or later
    if (fread(tbuf, 1, 3, bedfile) < 3) {
      goto wdist_ret_READ_FAIL;
    }
    if (memcmp(tbuf, "l\x1b\x01", 3)) {
      if (memcmp(tbuf, "l\x1b", 3)) {
	logprint("Error: Invalid header bytes in .bed file.\n");
	goto wdist_ret_INVALID_FORMAT;
      }
      bed_offset = 2;
    }
  } else if (llyy == 1LL) {
    // v0.99
    if (fread(tbuf, 1, 1, bedfile) != 1) {
      goto wdist_ret_READ_FAIL;
    }
    if (*tbuf == '\x01') {
      bed_offset = 1;
    } else if (*tbuf == '\0') {
      bed_offset = 2;
    } else {
      logprint("Error: Invalid header bytes in .bed file.\n");
      goto wdist_ret_INVALID_FORMAT;
    }
  } else if (llyy != 0LL) {
    sprintf(logbuf, "Error: Invalid .bed file size (expected %" PRIu64 " bytes).\n", ((uint64_t)unfiltered_indiv_ct4) * unfiltered_marker_ct);
    goto wdist_ret_INVALID_FORMAT_2;
  } else {
    // pre-0.99, no magic number, indiv-major
    bed_offset = 2;
  }
  if (bed_offset == 2) {
    strcpy(outname_end, ".bed.tmp"); // not really temporary
    logprint("Individual-major .bed file detected.  Transposing to SNP-major form.\n");
    fclose(bedfile);
    retval = indiv_major_to_snp_major(pedname, outname, &outfile, unfiltered_marker_ct);
    if (retval) {
      goto wdist_ret_2;
    }
    strcpy(pedname, outname);
    if (fopen_checked(&bedfile, pedname, "rb")) {
      goto wdist_ret_OPEN_FAIL;
    }
    bed_offset = 3;
  }
  if (mind_thresh < 1.0) {
    retval = mind_filter(bedfile, mind_thresh, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, bed_offset, missing_geno);
    if (retval) {
      goto wdist_ret_2;
    }
  }
  g_indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  if (!g_indiv_ct) {
    sprintf(logbuf, "Error: No %s pass QC.\n", species_plural);
    goto wdist_ret_INVALID_FORMAT_2;
  }
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  bitfield_andnot(pheno_nm, indiv_exclude, unfiltered_indiv_ctl);
  pheno_nm_ct = popcount_longs(pheno_nm, 0, unfiltered_indiv_ctl);
  if (!pheno_nm_ct) {
    logprint("Note: No phenotypes present.\n");
    hwe_all = 1;
  } else if (pheno_c) {
    pheno_ctrl_ct = popcount_longs_exclude(pheno_nm, pheno_c, unfiltered_indiv_ctl);
    if (pheno_nm_ct != g_indiv_ct) {
      sprintf(logbuf, "%u case%s, %u control%s, and %" PRIuPTR " missing phenotype%s present.\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "" : "s", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "" : "s", g_indiv_ct - pheno_nm_ct, (g_indiv_ct - pheno_nm_ct == 1)? "" : "s");
    } else {
      sprintf(logbuf, "%u case%s and %u control%s present.\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "" : "s", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "" : "s");
    }
    logprintb();
    if ((!hwe_all) && (!pheno_ctrl_ct)) {
      hwe_all = 1;
    }
  } else {
    if (pheno_nm_ct != g_indiv_ct) {
      sprintf(logbuf, "%u quantitative phenotype%s present (%" PRIuPTR " missing).\n", pheno_nm_ct, (pheno_nm_ct == 1)? "" : "s", g_indiv_ct - pheno_nm_ct);
    } else {
      sprintf(logbuf, "%u quantitative phenotype%s present.\n", pheno_nm_ct, (pheno_nm_ct == 1)? "" : "s");
    }
    logprintb();
    hwe_all = 1;
  }

  if (parallel_tot > g_indiv_ct / 2) {
    sprintf(logbuf, "Error: Too many --parallel jobs (maximum %" PRIuPTR "/2 = %" PRIuPTR ").\n", g_indiv_ct, g_indiv_ct / 2);
    goto wdist_ret_INVALID_CMDLINE_2;
  }
  if (g_thread_ct > 1) {
    if (calculation_type & (CALC_RELATIONSHIP | CALC_IBC | CALC_GDISTANCE_MASK | CALC_GROUPDIST | CALC_REGRESS_DISTANCE | CALC_GENOME | CALC_REGRESS_REL | CALC_UNRELATED_HERITABILITY | CALC_HARDY)) {
      sprintf(logbuf, "Using %d threads (change this with --threads).\n", g_thread_ct);
      logprintb();
    } else {
      logprint("Using 1 thread (no multithreaded calculations invoked).\n");
    }
  }

  nonfounders = (nonfounders || (!fam_col_34));
  wt_needed = distance_wt_req(calculation_type) && (!distance_flat_missing);
  hwe_needed = (hwe_thresh > 0.0) || (calculation_type & CALC_HARDY);
  retval = calc_freqs_and_hwe(bedfile, outname, outname_end, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_len, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, founder_info, nonfounders, maf_succ, set_allele_freqs, &marker_reverse, &marker_allele_cts, bed_offset, (unsigned char)missing_geno, hwe_needed, hwe_all, (pheno_nm_ct && pheno_c)? (calculation_type & CALC_HARDY) : 0, pheno_nm, pheno_nm_ct? pheno_c : NULL, &hwe_lls, &hwe_lhs, &hwe_hhs, &hwe_ll_cases, &hwe_lh_cases, &hwe_hh_cases, &hwe_ll_allfs, &hwe_lh_allfs, &hwe_hh_allfs, &hwe_hapl_allfs, &hwe_haph_allfs, &indiv_male_ct, &indiv_f_ct, &indiv_f_male_ct, wt_needed, &marker_weights_base, &marker_weights, exponent, chrom_info_ptr, sex_nm, sex_male, map_is_unsorted, &xmhh_exists, &nxmhh_exists);
  if (retval) {
    goto wdist_ret_2;
  }
  if (freqname) {
    retval = read_external_freqs(freqname, &freqfile, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, chrom_info_ptr, marker_alleles, max_marker_allele_len, marker_allele_cts, set_allele_freqs, maf_succ, missing_geno, exponent, wt_needed, marker_weights);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (!keep_allele_order) {
    calc_marker_reverse_bin(marker_reverse, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, set_allele_freqs);
  }

  // contrary to the PLINK flowchart, --freq effectively resolves before
  // --geno.
  if (calculation_type & CALC_FREQ) {
    if (freqx) {
      strcpy(outname_end, ".frqx");
    } else if (freq_counts) {
      strcpy(outname_end, ".frq.count");
    } else {
      strcpy(outname_end, ".frq");
    }
    retval = write_freqs(&outfile, outname, plink_maxsnp, unfiltered_marker_ct, marker_exclude, set_allele_freqs, chrom_info_ptr, marker_ids, max_marker_id_len, marker_alleles, max_marker_allele_len, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, hwe_hapl_allfs, hwe_haph_allfs, indiv_f_ct, indiv_f_male_ct, freq_counts, freqx, missing_geno, marker_reverse);
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_FREQ))))) {
      goto wdist_ret_2;
    }
  }
  if (geno_thresh < 1.0) {
    uii = binary_geno_filter(geno_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, g_indiv_ct, indiv_male_ct, marker_allele_cts, chrom_info_ptr);
    sprintf(logbuf, "%u marker%s removed due to missing genotype data (--geno).\n", uii, (uii == 1)? "" : "s");
    logprintb();
  }
  wkspace_reset(marker_allele_cts);
  marker_allele_cts = NULL;
  if (calculation_type & CALC_HARDY) {
    retval = hardy_report(threads, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_alleles, max_marker_allele_len, marker_reverse, hwe_lls, hwe_lhs, hwe_hhs, hwe_ll_cases, hwe_lh_cases, hwe_hh_cases, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, pheno_nm_ct, pheno_c, chrom_info_ptr);
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_FREQ | CALC_HARDY))))) {
      goto wdist_ret_2;
    }
  }
  if (hwe_thresh > 0.0) {
    enforce_hwe_threshold(hwe_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, hwe_lls, hwe_lhs, hwe_hhs, hwe_all, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs);
  }
  if ((min_maf != 0.0) || (max_maf != 0.5)) {
    enforce_maf_threshold(min_maf, max_maf, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, set_allele_freqs);
  }
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if (marker_ct == 0) {
    logprint("Error: All markers fail QC.\n");
    goto wdist_ret_INVALID_FORMAT;
  }

  if (wt_needed) {
    calc_marker_weights(exponent, marker_exclude, marker_ct, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, marker_weights);
  }
  wkspace_reset(hwe_lls);

  sprintf(logbuf, "%" PRIuPTR " marker%s and %" PRIuPTR " %s pass filters and QC%s.\n", marker_ct, (marker_ct == 1)? "" : "s", g_indiv_ct, species_str(g_indiv_ct), (calculation_type & CALC_REL_CUTOFF)? " (before --rel-cutoff)": "");
  logprintb();

  if (refalleles) {
    if (marker_alleles) {
      retval = load_ref_alleles(refalleles, unfiltered_marker_ct, marker_alleles, max_marker_allele_len, marker_reverse, marker_ids, max_marker_id_len);
      if (retval) {
	goto wdist_ret_2;
      }
    } else {
      logprint("Note: Ignoring --reference-allele, since allele IDs are not used in this run.\n");
    }
  }

  if (calculation_type & CALC_WRITE_SNPLIST) {
    retval = write_snplist(&outfile, outname, outname_end, marker_exclude, marker_ct, marker_ids, max_marker_id_len);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_MAKE_BED) {
    retval = make_bed(bedfile, bed_offset, bimfile, map_cols, &bedtmpfile, &famtmpfile, &bimtmpfile, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_alleles, max_marker_allele_len, marker_reverse, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, missing_phenod, output_missing_pheno, max_marker_id_len, map_is_unsorted, indiv_sort_map, set_hh_missing, chrom_info_ptr);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_RECODE) {
    retval = recode(recode_modifier, bedfile, bed_offset, famfile, bimfile, &outfile, outname, outname_end, recode_allele_name, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, marker_ids, max_marker_id_len, marker_alleles, max_marker_allele_len, marker_reverse, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, missing_phenod, output_missing_geno, output_missing_pheno, set_hh_missing, xmhh_exists, nxmhh_exists, chrom_info_ptr);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_LD_PRUNE) {
    retval = ld_prune(bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_len, chrom_info_ptr, set_allele_freqs, marker_pos, unfiltered_indiv_ct, indiv_exclude, sex_nm, sex_male, ld_window_size, ld_window_kb, ld_window_incr, ld_last_param, outname, outname_end, calculation_type);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_REGRESS_PCS) {
    // do this before marker_alleles is overwritten in memory...
    retval = calc_regress_pcs(evecname, regress_pcs_normalize_pheno, regress_pcs_sex_specific, regress_pcs_clip, max_pcs, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_len, marker_alleles, max_marker_allele_len, chrom_info_ptr, marker_pos, g_indiv_ct, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len, sex_nm, sex_male, pheno_d, missing_phenod, &outfile, outname, outname_end);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  // no more need for marker_ids/marker_alleles, unload?  (probably want to
  // initially load at far end of stack to make this workable...)

  if (wt_needed) {
    // N.B. marker_weights is currently on top of the stack
    wkspace_reset(marker_weights_base);
    // normalize included marker weights to add to 2^32 - 1
    dxx = 0.0;
    marker_uidx = 0;
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      dxx += marker_weights[marker_uidx++];
    }
    dxx = 4294967295.0 / dxx;
    marker_idx = 0;
    marker_weights_i = (uint32_t*)wkspace_alloc(marker_ct * sizeof(int32_t));
    for (marker_uidx = 0; marker_idx < marker_ct; marker_uidx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      marker_weights_i[marker_idx++] = (uint32_t)(marker_weights[marker_uidx] * dxx + 0.5);
    }
  }
  if (calculation_type & CALC_GENOME) {
    if (!id_buf) {
      id_buf = (char*)malloc(max_person_id_len * sizeof(char));
      if (!id_buf) {
	goto wdist_ret_NOMEM;
      }
    }
    retval = populate_pedigree_rel_info(&pri, unfiltered_indiv_ct, id_buf, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
    if (retval) {
      goto wdist_ret_2;
    }
  }
  wkspace_mark2 = wkspace_base;

  if (relationship_or_ibc_req(calculation_type)) {
    if (rel_calc_type & REL_CALC_SINGLE_PREC) {
      retval = calc_rel_f(threads, parallel_idx, parallel_tot, calculation_type, rel_calc_type, bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, var_std, ibc_type, (float)rel_cutoff, set_allele_freqs, chrom_info_ptr);
    } else {
      retval = calc_rel(threads, parallel_idx, parallel_tot, calculation_type, rel_calc_type, bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, var_std, ibc_type, rel_cutoff, set_allele_freqs, &rel_ibc, chrom_info_ptr);
    }
    if (retval) {
      goto wdist_ret_2;
    }

    if (calculation_type & CALC_REGRESS_REL) {
      retval = regress_rel_main(indiv_exclude, g_indiv_ct, regress_rel_iters, regress_rel_d, threads, pheno_d);
      if (retval) {
	goto wdist_ret_2;
      }
    }
#ifndef NOLAPACK
    if (calculation_type & CALC_UNRELATED_HERITABILITY) {
      retval = calc_unrelated_herit(calculation_type, ibc_type, unfiltered_indiv_ct, indiv_exclude, pheno_d, rel_ibc, unrelated_herit_covg, unrelated_herit_covr, unrelated_herit_tol);
      if (retval) {
	goto wdist_ret_2;
      }
    }
#endif
    if (calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX) && (!(calculation_type & CALC_REL_CUTOFF))) {
      wkspace_reset(g_rel_dists);
    } else {
      wkspace_reset(g_indiv_missing_unwt);
      g_indiv_missing_unwt = NULL;
      g_missing_dbl_excluded = NULL;
    }
  }

  if (calculation_type & CALC_REGRESS_PCS_DISTANCE) {
    logprint("Error: --regress-pcs-distance has not yet been written.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto wdist_ret_2;
  } else if (distance_req(calculation_type)) {
    retval = calc_distance(threads, parallel_idx, parallel_tot, bedfile, bed_offset, &outfile, outname, outname_end, calculation_type, dist_calc_type, distance_flat_missing, marker_exclude, marker_ct, set_allele_freqs, unfiltered_indiv_ct, unfiltered_indiv_ct4, indiv_exclude, person_ids, max_person_id_len, chrom_info_ptr, wt_needed, marker_weights_i, exp0, exponent);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_LOAD_DISTANCES) {
    dists_alloc = ((intptr_t)g_indiv_ct * (g_indiv_ct - 1)) * (sizeof(double) / 2);
    g_dists = (double*)wkspace_alloc(dists_alloc);
    if (!g_dists) {
      goto wdist_ret_NOMEM;
    }
    if (fopen_checked(&loaddistfile, loaddistname, "rb")) {
      goto wdist_ret_OPEN_FAIL;
    }
    if (fseeko(loaddistfile, 0, SEEK_END)) {
      goto wdist_ret_READ_FAIL;
    }
    if (ftello(loaddistfile) != (int64_t)dists_alloc) {
      sprintf(logbuf, "Invalid --load-dists file.  (Triangular binary of size %" PRIu64 " expected.)\n", dists_alloc);
      goto wdist_ret_INVALID_FORMAT_2;
    }
    rewind(loaddistfile);
    if (fread(g_dists, 1, dists_alloc, loaddistfile) < dists_alloc) {
      goto wdist_ret_READ_FAIL;
    }
    fclose(loaddistfile);
  }

  if (calculation_type & CALC_GROUPDIST) {
    retval = groupdist_calc(threads, unfiltered_indiv_ct, indiv_exclude, groupdist_iters, groupdist_d, pheno_nm, pheno_c);
    if (retval) {
      goto wdist_ret_2;
    }
  }
  if (calculation_type & CALC_REGRESS_DISTANCE) {
    retval = regress_distance(calculation_type, g_dists, pheno_d, unfiltered_indiv_ct, indiv_exclude, g_thread_ct, regress_iters, regress_d);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_GENOME) {
    wkspace_reset(wkspace_mark2);
    retval = calc_genome(threads, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, marker_pos, set_allele_freqs, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info, parallel_idx, parallel_tot, outname, outname_end, nonfounders, calculation_type, genome_output_gz, genome_output_full, genome_ibd_unbounded, ppc_gap, pheno_nm, pheno_c, pri);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_MODEL) {
    retval = model_assoc(threads, bedfile, bed_offset, outname, outname_end, calculation_type, model_modifier, model_cell_ct, model_mperm_val, ci_size, ci_zt, pfilter, mtest_adjust, adjust_lambda, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_pos, marker_alleles, max_marker_allele_len, marker_reverse, chrom_info_ptr, unfiltered_indiv_ct, aperm_min, aperm_max, aperm_alpha, aperm_beta, aperm_init_interval, aperm_interval_slope, pheno_nm_ct, pheno_nm, pheno_c, pheno_d, sex_nm, sex_male);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  while (0) {
  wdist_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  wdist_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  wdist_ret_INVALID_FORMAT_2:
    logprintb();
  wdist_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  wdist_ret_INVALID_CMDLINE_2:
    logprintb();
  wdist_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  wdist_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
 wdist_ret_2:
  fclose_cond(bedtmpfile);
  fclose_cond(bimtmpfile);
  fclose_cond(famtmpfile);
 wdist_ret_1:
  free_cond(pheno_d);
  free_cond(pheno_c);
  free_cond(phenor_d);
  free_cond(phenor_c);
  free_cond(id_buf);
  free_cond(id_list);
  free_cond(pid_list);
  fclose_cond(loaddistfile);
  fclose_cond(freqfile);
  fclose_cond(filterfile);
  fclose_cond(phenofile);
  fclose_cond(famfile);
  fclose_cond(bimfile);
  fclose_cond(bedfile);
  fclose_cond(outfile3);
  fclose_cond(outfile2);
  fclose_cond(outfile);
  return retval;
}

// output-missing-phenotype + terminating null
#define MAX_FLAG_LEN 25

static inline int32_t is_flag(char* param) {
  char cc = param[1];
  return ((*param == '-') && ((cc > '9') || (cc < '0'))); 
}

static inline char* is_flag_start(char* param) {
  char cc = param[1];
  if ((*param == '-') && ((cc > '9') || (cc < '0'))) {
    return (cc == '-')? (&(param[2])) : (&(param[1]));
  }
  return NULL;
}

int32_t param_count(int32_t argc, char** argv, int32_t flag_idx) {
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

int32_t enforce_param_ct_range(int32_t argc, char** argv, int32_t flag_idx, int32_t min_ct, int32_t max_ct, int32_t* param_ct_ptr) {
  int32_t param_ct = param_count(argc, argv, flag_idx);
  if (param_ct > max_ct) {
    if (max_ct > min_ct) {
      sprintf(logbuf, "Error: %s accepts at most %d parameter%s.%s", argv[flag_idx], max_ct, (max_ct == 1)? "" : "s", errstr_append);
    } else {
      sprintf(logbuf, "Error: %s only accepts %d parameter%s.%s", argv[flag_idx], max_ct, (max_ct == 1)? "" : "s", errstr_append);
    }
    return -1;
  } else if (param_ct < min_ct) {
    if (min_ct == 1) {
      sprintf(logbuf, "Error: Missing %s parameter.%s", argv[flag_idx], errstr_append);
    } else if (min_ct < max_ct) {
      sprintf(logbuf, "Error: %s requires %s%d parameters.%s", argv[flag_idx], (min_ct < max_ct)? "at least " : "", min_ct, errstr_append);
    }
    return -1;
  }
  *param_ct_ptr = param_ct;
  return 0;
}

int32_t parse_next_range(uint32_t param_ct, char range_delim, char** argv, uint32_t* cur_param_idx_ptr, char** cur_arg_pptr, char** range_start_ptr, uint32_t* rs_len_ptr, char** range_end_ptr, uint32_t* re_len_ptr) {
  // Starts reading from argv[cur_param_idx][cur_pos].  If a valid range is
  // next, range_start + rs_len + range_end + re_len are updated.  If only a
  // single item is next, range_end is set to NULL and range_start + rs_len are
  // updated.  If there are no items left, range_start is set to NULL.  If
  // the input is not well-formed, -1 is returned instead of 0.
  uint32_t cur_param_idx = *cur_param_idx_ptr;
  char* cur_arg_ptr = *cur_arg_pptr;
  char cc;
  if (cur_param_idx > param_ct) {
    *cur_arg_pptr = NULL;
    return 0;
  }
  while (1) {
    cc = *cur_arg_ptr;
    if (!cc) {
      *cur_param_idx_ptr = ++cur_param_idx;
      if (cur_param_idx > param_ct) {
	*range_start_ptr = NULL;
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
      *range_end_ptr = NULL;
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

int32_t parse_chrom_ranges(uint32_t param_ct, char range_delim, char** argv, uint64_t* chrom_mask_ptr, uint32_t species, char* cur_flag_str) {
  uint32_t argct = 0;
  uint32_t cur_param_idx = 1;
  uint64_t chrom_mask = *chrom_mask_ptr;
  char* cur_arg_ptr;
  char* range_start;
  uint32_t rs_len;
  char* range_end;
  uint32_t re_len;
  int32_t marker_code_start;
  int32_t marker_code_end;
  if (param_ct) {
    cur_arg_ptr = argv[1];
    while (1) {
      if (parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	sprintf(logbuf, "Error: Invalid --%s parameter '%s'.%s", cur_flag_str, argv[cur_param_idx], errstr_append);
	return -1;
      }
      if (!range_start) {
	break;
      }
      marker_code_start = marker_code2(species, range_start, rs_len);
      if (marker_code_start == -1) {
	range_start[rs_len] = '\0';
	sprintf(logbuf, "Error: Invalid --%s chromosome code '%s'.%s", cur_flag_str, range_start, errstr_append);
	return -1;
      }
      if (range_end) {
        marker_code_end = marker_code2(species, range_end, re_len);
	if (marker_code_end == -1) {
	  range_end[re_len] = '\0';
	  sprintf(logbuf, "Error: Invalid --%s chromosome code '%s'.%s", cur_flag_str, range_end, errstr_append);
	  return -1;
	}
        if (marker_code_end <= marker_code_start) {
	  range_start[rs_len] = '\0';
	  range_end[re_len] = '\0';
	  sprintf(logbuf, "Error: --%s chromosome code '%s' is not greater than '%s'.%s", cur_flag_str, range_end, range_start, errstr_append);
	  return -1;
	}
	chrom_mask |= (1LLU << (marker_code_end + 1)) - (1LLU << (marker_code_start));
      } else {
	chrom_mask |= 1LLU << marker_code_start;
      }
      argct++;
    }
  }
  if (!argct) {
    sprintf(logbuf, "Error: --%s requires at least one value.%s", cur_flag_str, errstr_append);
    return -1;
  }
  *chrom_mask_ptr = chrom_mask;
  return 0;
}

int32_t parse_marker_ranges(uint32_t param_ct, char range_delim, char** argv, char** snps_flag_markers_ptr, unsigned char** snps_flag_starts_range_ptr, uint32_t* snps_flag_ct_ptr, uint32_t* snps_flag_max_len_ptr) {
  uint32_t snps_flag_ct = 0;
  uint32_t cur_param_idx = 1;
  uint32_t snps_flag_max_len = 0;
  char* cur_arg_ptr;
  char* range_start;
  uint32_t rs_len;
  char* range_end;
  uint32_t re_len;
  char* cur_snps_flag_marker_str;
  char* dup_check;
  unsigned char* cur_snps_flag_starts_range;
  // two passes.  first pass: count parameters, determine snps_flag_max_len;
  // then allocate memory; then fill it.
  if (param_ct) {
    cur_arg_ptr = argv[1];
    while (1) {
      if (parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	sprintf(logbuf, "Error: Invalid --snps parameter '%s'.%s", argv[cur_param_idx], errstr_append);
        logprintb();
        return RET_INVALID_CMDLINE;
      }
      if (!range_start) {
	break;
      }
      snps_flag_ct++;
      if (rs_len > snps_flag_max_len) {
	snps_flag_max_len = rs_len; // does NOT include trailing null yet
      }
      if (range_end) {
	snps_flag_ct++;
	if (re_len > snps_flag_max_len) {
	  snps_flag_max_len = re_len;
	}
      }
    }
  }
  if (!snps_flag_ct) {
    sprintf(logbuf, "Error: --snps requires at least one value.%s", errstr_append);
    logprintb();
    return RET_INVALID_CMDLINE;
  }
  *snps_flag_max_len_ptr = ++snps_flag_max_len;
  *snps_flag_ct_ptr = snps_flag_ct;
  *snps_flag_markers_ptr = (char*)malloc(snps_flag_ct * snps_flag_max_len * sizeof(char));
  if (!(*snps_flag_markers_ptr)) {
    return RET_NOMEM;
  }
  *snps_flag_starts_range_ptr = (unsigned char*)malloc(snps_flag_ct * sizeof(char));
  if (!(*snps_flag_starts_range_ptr)) {
    return RET_NOMEM;
  }
  cur_snps_flag_marker_str = *snps_flag_markers_ptr;
  cur_snps_flag_starts_range = *snps_flag_starts_range_ptr;
  cur_param_idx = 1;
  cur_arg_ptr = argv[1];
  while (1) {
    parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len);
    if (!range_start) {
      return 0;
    }
    memcpy(cur_snps_flag_marker_str, range_start, rs_len);
    cur_snps_flag_marker_str[rs_len] = '\0';
    dup_check = *snps_flag_markers_ptr;
    while (dup_check < cur_snps_flag_marker_str) {
      if (!memcmp(dup_check, cur_snps_flag_marker_str, rs_len + 1)) {
	sprintf(logbuf, "Error: Duplicate --snps marker ID '%s'.\n", cur_snps_flag_marker_str);
	logprintb();
	return RET_INVALID_CMDLINE;
      }
      dup_check = &(dup_check[snps_flag_max_len]);
    }
    cur_snps_flag_marker_str = &(cur_snps_flag_marker_str[snps_flag_max_len]);
    if (range_end) {
      *cur_snps_flag_starts_range++ = 1;
      memcpy(cur_snps_flag_marker_str, range_end, re_len);
      cur_snps_flag_marker_str[re_len] = '\0';
      dup_check = *snps_flag_markers_ptr;
      while (dup_check < cur_snps_flag_marker_str) {
	if (!memcmp(dup_check, cur_snps_flag_marker_str, rs_len + 1)) {
	  sprintf(logbuf, "Error: Duplicate --snps marker ID '%s'.\n", cur_snps_flag_marker_str);
	  logprintb();
	  return RET_INVALID_CMDLINE;
	}
        dup_check = &(dup_check[snps_flag_max_len]);
      }
      cur_snps_flag_marker_str = &(cur_snps_flag_marker_str[snps_flag_max_len]);
      *cur_snps_flag_starts_range++ = 0;
    } else {
      *cur_snps_flag_starts_range++ = 0;
    }
  }
}

int32_t species_flag(Chrom_info* chrom_info_ptr, int32_t species_val) {
  if (chrom_info_ptr->species) {
    logprint("Error: Multiple species flags.\n");
    return -1;
  }
  chrom_info_ptr->species = species_val;
  return 0;
}

void invalid_arg(char* argv) {
  sprintf(logbuf, "Error: Invalid flag ('%s').%s%s", argv, (argv[0] == '-')? "" : "  All flags must be preceded by 1-2 dashes.", errstr_append);
}

void print_ver() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}

int32_t alloc_string(char** sbuf, char* source) {
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
    sprintf(logbuf, "Error: --%s filename too long.\n", argptr);
    logprintb();
    return RET_OPEN_FAIL;
  }
  *fnbuf = (char*)malloc((slen + extra_size) * sizeof(char));
  if (!(*fnbuf)) {
    return RET_NOMEM;
  }
  memcpy(*fnbuf, source, slen);
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

int32_t recode_type_set(uint32_t* recode_modifier_ptr, uint32_t cur_code) {
  if (*recode_modifier_ptr & (RECODE_TYPEMASK - cur_code)) {
    sprintf(logbuf, "Error: Conflicting --recode modifiers.%s", errstr_append);
    return -1;
  }
  *recode_modifier_ptr |= cur_code;
  return 0;
}

// int32_t enforce_param_ct_range(int32_t argc, char** argv, int32_t flag_idx, int32_t min_ct, int32_t max_ct, int32_t* param_ct_ptr) {

int32_t no_params_check(int32_t argc, char** argv, int32_t cur_arg, char* argptr, char* outname_end) {
  uint32_t pc = param_count(argc, argv, cur_arg);
  if (!pc) {
    return 0;
  }
  sprintf(logbuf, "Error: --%s doesn't accept parameters.%s%s", argptr, ((pc == 1) && (!outname_end))? "  (Did you forget '--out'?)" : "", errstr_append);
  return -1;
}

int32_t main(int32_t argc, char** argv) {
  unsigned char* wkspace_ua;
  char outname[FNAMESIZE];
  char mapname[FNAMESIZE];
  char pedname[FNAMESIZE];
  char famname[FNAMESIZE];
  char genname[FNAMESIZE];
  char samplename[FNAMESIZE];
  char mergename1[FNAMESIZE];
  char mergename2[FNAMESIZE];
  char mergename3[FNAMESIZE];
  time_t rawtime;
  char* argptr;
  char* sptr;
  char* outname_end = NULL;
  char** subst_argv = NULL;
  char* script_buf = NULL;
  char* rerun_buf = NULL;
  char* flag_buf = NULL;
  uint32_t* flag_map = NULL;
  char* makepheno_str = NULL;
  char* phenoname_str = NULL;
  char* refalleles = NULL;
  char* filterval = NULL;
  char* evecname = NULL;
  char* filtername = NULL;
  char* loaddistname = NULL;
  char* freqname = NULL;
  char* extractname = NULL;
  char* excludename = NULL;
  char* keepname = NULL;
  char* removename = NULL;
  char* phenoname = NULL;
  char* recode_allele_name = NULL;
  char* lgen_reference_fname = NULL;
  char* covar_fname = NULL;
  char* covar_str = NULL;
  char* cluster_fname = NULL;
  char* set_fname = NULL;
  char* subset_fname = NULL;
  char* update_alleles_fname = NULL;
  char* update_map_fname = NULL;
  char* update_ids_fname = NULL;
  char* update_parents_fname = NULL;
  char* update_sex_fname = NULL;
  char* loop_assoc_fname = NULL;
  char** subst_argv2;
  int32_t retval = 0;
  int32_t load_params = 0; // describes what file parameters have been provided
  int32_t load_rare = 0;
  int32_t fam_col_1 = 1;
  int32_t fam_col_34 = 1;
  int32_t fam_col_5 = 1;
  int32_t fam_col_6 = 1;
  int32_t mpheno_col = 0;
  int32_t mcovar_col = 0;
  int32_t affection_01 = 0;
  double exponent = 0.0;
  double min_maf = 0.0;
  double max_maf = 0.5;
  double geno_thresh = 1.0;
  double mind_thresh = 1.0;
  double hwe_thresh = 0.0;
  int32_t hwe_all = 0;
  double rel_cutoff = 0.025;
  int32_t cur_arg = 1;
  uint64_t calculation_type = 0;
  int32_t rel_calc_type = 0;
  int32_t dist_calc_type = 0;
  char* bubble;
  int32_t mfilter_col = 0;
  uint32_t pheno_modifier = 0;
  int32_t tail_pheno = 0;
  int32_t prune = 0;
  int32_t missing_pheno = -9;
  char missing_geno = '0';
  char output_missing_pheno[32];
  char output_missing_geno = '0';
  double tail_bottom;
  double tail_top;
  int32_t distance_3d = 0;
  int32_t distance_flat_missing = 0;
  uintptr_t groupdist_iters = ITERS_DEFAULT;
  int32_t groupdist_d = 0;
  uintptr_t regress_iters = ITERS_DEFAULT;
  int32_t regress_d = 0;
  uintptr_t regress_rel_iters = ITERS_DEFAULT;
  int32_t regress_rel_d = 0;
  double unrelated_herit_tol = 0.0000001;
  double unrelated_herit_covg = 0.45;
  double unrelated_herit_covr = 0.55;
  int32_t ibc_type = 0;
  int32_t parallel_idx = 0;
  int32_t parallel_tot = 1;
  int32_t allow_no_sex = 0;
  int32_t must_have_sex = 0;
  int32_t nonfounders = 0;
  uint32_t write_covar_modifier = 0;
  uint32_t model_modifier = 0;
  int32_t model_cell_ct = -1;
  uint32_t gxe_mcovar = 0;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t mm;
  int32_t nn;
  int32_t oo;
  int32_t ppc_gap = DEFAULT_PPC_GAP;
  unsigned long int rseed = 0;
  int32_t genome_output_gz = 0;
  int32_t genome_output_full = 0;
  int32_t genome_ibd_unbounded = 0;
  FILE* scriptfile = NULL;
  int32_t num_params;
  int32_t in_param;
  int32_t ld_window_size = 0;
  int32_t ld_window_incr = 0;
  double ld_last_param = 0.0;
  int32_t ld_window_kb = 0;
  int32_t maf_succ = 0;
  int32_t filter_case_control = 0;
  int32_t filter_sex = 0;
  int32_t filter_founder_nonf = 0;
  int32_t regress_pcs_normalize_pheno = 0;
  int32_t regress_pcs_sex_specific = 0;
  int32_t regress_pcs_clip = 0;
  int32_t max_pcs = MAX_PCS_DEFAULT;
  int32_t freqx = 0;
  uint32_t recode_modifier = 0;
  int32_t freq_counts = 0;
  int32_t allelexxxx = 0;
  int32_t merge_type = 0;
  int32_t indiv_sort = 0;
  int32_t keep_allele_order = 0;
  uint32_t cur_flag = 0;
  uint32_t flag_ct = 0;
  int32_t dummy_marker_ct = 0;
  int32_t dummy_indiv_ct = 0;
  uint32_t dummy_flags = 0;
  double dummy_missing_geno = 0.0;
  double dummy_missing_pheno = 0.0;
  char* markername_from = NULL;
  char* markername_to = NULL;
  char* markername_snp = NULL;
  uint32_t snp_window_size = 0;
  int32_t marker_pos_start = -1;
  int32_t marker_pos_end = -1;
  uint32_t lgen_modifier = 0;
  uint32_t covar_modifier = 0;
  uint32_t update_map_modifier = 0;
  uint32_t model_mperm_val = 0;
  uint32_t mperm_val = 0;
  uint32_t aperm_min = 5;
  uint32_t aperm_max = 1000000;
  double aperm_alpha = 0;
  double aperm_beta = 0.0001;
  double aperm_init_interval = 1;
  double aperm_interval_slope = 0.001;
  double ci_size = 0.0;
  double pfilter = 1.0;
  uint32_t mtest_adjust = 0;
  double adjust_lambda = 0.0;
  Chrom_info chrom_info;
  char* argptr2;
  char* flagptr;
  char* missing_code = NULL;
  char range_delim = '-';
  uint64_t chrom_exclude = 0;
  char* snps_flag_markers = NULL;
  unsigned char* snps_flag_starts_range = NULL;
  uint32_t snps_flag_ct = 0;
  uint32_t set_hh_missing = 0;
  uint32_t snps_flag_max_len;
  double dxx;
  char cc;
  uint32_t uii;
  intptr_t default_alloc_mb;
  int64_t llxx;
#ifdef __APPLE__
  int32_t mib[2];
  size_t sztmp;
#endif
#if _WIN32
  SYSTEM_INFO sysinfo;
  MEMORYSTATUSEX memstatus;
  DWORD windows_dw; // no idea why the f*** a uint32_t doesn't work
#endif
  for (ii = 1; ii < argc; ii++) {
    if ((!memcmp("-script", argv[ii], 8)) || (!memcmp("--script", argv[ii], 9))) {
      if (enforce_param_ct_range(argc, argv, ii, 1, 1, &jj)) {
	print_ver();
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      for (jj = ii + 2; jj < argc; jj++) {
	if ((!memcmp("-script", argv[jj], 8)) || (!memcmp("--script", argv[jj], 9))) {
	  print_ver();
	  printf("Error: Multiple --script flags.  Merge the files into one.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE;
	}
      }
      // logging not yet active, so don't use fopen_checked()
      scriptfile = fopen(argv[ii + 1], "rb");
      if (!scriptfile) {
	print_ver();
	printf(errstr_fopen, argv[ii + 1]);
	goto main_ret_OPEN_FAIL;
      }
      if (fseeko(scriptfile, 0, SEEK_END)) {
	print_ver();
	goto main_ret_READ_FAIL;
      }
      llxx = ftello(scriptfile);
      if (llxx > MAXLINELEN) {
	print_ver();
	printf("Error: Script file too long (max %u bytes).\n", MAXLINELEN);
	retval = RET_INVALID_FORMAT;
	goto main_ret_1;
      }
      jj = llxx;
      rewind(scriptfile);
      script_buf = (char*)malloc(jj);
      if (!script_buf) {
	print_ver();
	goto main_ret_NOMEM;
      }
      kk = fread(script_buf, 1, jj, scriptfile);
      if (kk < jj) {
	print_ver();
	goto main_ret_READ_FAIL;
      }
      fclose_null(&scriptfile);
      num_params = 0;
      in_param = 0;
      for (kk = 0; kk < jj; kk++) {
	if (is_space_or_eoln(tbuf[kk])) {
	  in_param = 0;
	} else if (!in_param) {
	  num_params++;
	  in_param = 1;
	}
      }
      subst_argv = (char**)malloc((num_params + argc - 3) * sizeof(char*));
      num_params = 0;
      in_param = 0;
      for (kk = 1; kk < ii; kk++) {
        subst_argv[num_params++] = argv[kk];
      }
      for (kk = 0; kk < jj; kk++) {
	if (is_space_or_eoln(tbuf[kk])) {
	  if (in_param) {
	    tbuf[kk] = '\0';
	    in_param = 0;
	  }
	} else if (!in_param) {
	  subst_argv[num_params++] = &(tbuf[kk]);
	  in_param = 1;
	}
      }
      for (jj = ii + 2; jj < argc; jj++) {
	subst_argv[num_params++] = argv[jj];
      }
      argc = num_params;
      cur_arg = 0;
      argv = subst_argv;
    }
  }
  for (ii = cur_arg; ii < argc; ii++) {
    if ((!memcmp("-rerun", argv[ii], 7)) || (!memcmp("--rerun", argv[ii], 8))) {
      if (enforce_param_ct_range(argc, argv, ii, 0, 1, &jj)) {
	print_ver();
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      for (kk = ii + jj + 1; kk < argc; kk++) {
	if ((!memcmp("-rerun", argv[kk], 7)) || (!memcmp("--rerun", argv[kk], 8))) {
	  print_ver();
	  fputs("Error: Duplicate --rerun flag.\n", stdout);
	  goto main_ret_INVALID_CMDLINE;
	}
      }
      if (jj) {
	scriptfile = fopen(argv[ii + 1], "r");
      } else {
	scriptfile = fopen("wdist.log", "r");
      }
      if (!scriptfile) {
	print_ver();
	goto main_ret_OPEN_FAIL;
      }
      if (!fgets(tbuf, MAXLINELEN, scriptfile)) {
	print_ver();
	fputs("Error: Empty log file for --rerun.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      if (!fgets(tbuf, MAXLINELEN, scriptfile)) {
	print_ver();
	fputs("Error: Only one line in --rerun log file.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      fclose_null(&scriptfile);
      kk = atoi(tbuf);
      if ((kk < 1) || (kk > MAXLINELEN)) {
	print_ver();
	fputs("Error: Improperly formatted --rerun log file.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      uii = strlen(tbuf) + 1;
      if (uii == MAXLINELEN) {
	print_ver();
	fputs("Error: Second line too long in --rerun log file.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      rerun_buf = (char*)malloc(uii);
      memcpy(rerun_buf, tbuf, uii);

      memset(tbuf, 1, kk);
      sptr = next_item_mult(rerun_buf, 2);
      mm = 0;
      oo = 0;
      do {
	if (no_more_items_kns(sptr)) {
	  print_ver();
	  fputs("Error: Improperly formatted --rerun log file.\n", stdout);
	  goto main_ret_INVALID_FORMAT;
	}
	argptr = is_flag_start(sptr);
	if (argptr) {
	  uii = strlen_se(argptr);
          for (nn = cur_arg; nn < argc; nn++) {
	    argptr2 = is_flag_start(argv[nn]);
	    if (argptr2) {
	      if (!memcmp(argptr, argptr2, uii)) {
		nn = -1;
		break;
	      }
	    }
	  }
          if (nn == -1) {
	    // matching flag, override --rerun
            do {
	      oo++;
	      tbuf[mm++] = 0;
	      if (mm == kk) {
		break;
	      }
	      sptr = next_item(sptr);
	    } while (!is_flag(sptr));
	  } else {
	    mm++;
	    sptr = next_item(sptr);
	  }
	} else {
	  mm++;
          sptr = next_item(sptr);
	}
      } while (mm < kk);
      subst_argv2 = (char**)malloc((argc + kk - oo - jj - 1 - cur_arg) * sizeof(char*));
      if (!subst_argv2) {
	print_ver();
	goto main_ret_NOMEM;
      }
      oo = 0;
      for (mm = cur_arg; mm < ii; mm++) {
	subst_argv2[oo++] = argv[mm];
      }
      sptr = next_item_mult(rerun_buf, 2);
      for (mm = 0; mm < kk; mm++) {
        if (tbuf[mm]) {
	  uii = strlen_se(sptr);
	  subst_argv2[oo++] = sptr;
	  sptr[uii] = '\0';
	  if (mm != kk - 1) {
	    sptr = skip_initial_spaces(&(sptr[uii + 1]));
	  }
	} else {
	  sptr = next_item(sptr);
	}
      }
      for (mm = ii + jj + 1; mm < argc; mm++) {
	subst_argv2[oo++] = argv[mm];
      }
      cur_arg = 0;
      argc = oo;
      if (subst_argv) {
	free(subst_argv);
      }
      subst_argv = subst_argv2;
      argv = subst_argv2;
      subst_argv2 = NULL;
    }
  }
  if ((cur_arg < argc) && (!is_flag(argv[cur_arg]))) {
    print_ver();
    printf("Error: First parameter must be a flag.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE;
  }
  flag_ct = 0;
  for (ii = cur_arg; ii < argc; ii++) {
    argptr = is_flag_start(argv[ii]);
    if (argptr) {
      if (!memcmp("help", argptr, 5)) {
	if ((cur_arg != 1) || (ii != 1) || subst_argv) {
	  fputs("--help present, ignoring other flags.\n", stdout);
	}
	print_ver();
	retval = disp_help(argc - cur_arg - 1, &(argv[cur_arg + 1]));
	goto main_ret_1;
      }
      if (strlen(argptr) >= MAX_FLAG_LEN) {
	print_ver();
	invalid_arg(argv[ii]);
	fputs(logbuf, stdout);
        goto main_ret_INVALID_CMDLINE;
      }
      flag_ct++;
    }
  }
  if (!flag_ct) {
    print_ver();
    fputs(notestr_null_calc, stdout);
    fputs(cmdline_format_str, stdout);
    fputs(notestr_null_calc2, stdout);
    retval = RET_NULL_CALC;
    goto main_ret_1;
  }
  flag_buf = (char*)malloc(flag_ct * MAX_FLAG_LEN * sizeof(char));
  flag_map = (uint32_t*)malloc(flag_ct * sizeof(int32_t));
  if ((!flag_buf) || (!flag_map)) {
    print_ver();
    goto main_ret_NOMEM;
  }
  flagptr = flag_buf;
  mm = 0; // parameter count increase due to aliases
  for (ii = cur_arg; ii < argc; ii++) {
    argptr = is_flag_start(argv[ii]);
    if (argptr) {
      kk = strlen(argptr) + 1;
      // handle aliases now, so sorting will have the desired effects
      switch (*argptr) {
      case 'Z':
	if (!memcmp(argptr, "Z-genome", 9)) {
	  memcpy(flagptr, "genome gz", 10);
	  mm++;
	  break;
	}
	goto main_flag_copy;
      case 'a':
	if ((kk == 10) && (!memcmp(argptr, "allele", 6))) {
	  if ((tolower(argptr[6]) == 'a') && (tolower(argptr[7]) == 'c') && (tolower(argptr[8]) == 'g') && (tolower(argptr[9]) == 't')) {
	    memcpy(flagptr, "alleleACGT", 11);
	    break;
	  }
	}
	goto main_flag_copy;
      case 'f':
	if (!memcmp(argptr, "frqx", 5)) {
	  memcpy(flagptr, "freqx", 6);
	  break;
	}
	goto main_flag_copy;
      case 'g':
	if (!memcmp(argptr, "grm-cutoff", 11)) {
          memcpy(flagptr, "rel-cutoff", 11);
	  break;
	}
	goto main_flag_copy;
      case 'l':
	if (!memcmp(argptr, "list", 5)) {
	  memcpy(flagptr, "recode list", 12);
	  printf("Note: --list flag deprecated.  Use '--recode list' instead.\n");
	  break;
	}
      case 'm':
	if (!memcmp(argptr, "missing_code", 13)) {
	  memcpy(flagptr, "missing-code", 13);
	  break;
	}
	goto main_flag_copy;
      case 'r':
	if (!memcmp(argptr, "recode", 6)) {
	  if (((kk == 9) && ((!memcmp(&(argptr[6]), "12", 2)) || ((tolower(argptr[6]) == 'a') && (tolower(argptr[7]) == 'd')))) || (!memcmp(&(argptr[6]), "-lgen", 6)) || (!memcmp(&(argptr[6]), "-rlist", 7)) || ((tolower(argptr[6]) == 'a') && (kk == 8))) {
	    if (kk == 13) {
	      memcpy(flagptr, "recode rlist", 13);
	    } else if (kk == 12) {
	      memcpy(flagptr, "recode lgen", 12);
	    } else if (argptr[6] == '1') {
	      memcpy(flagptr, "recode 12", 10);
	    } else if (kk == 9) {
	      memcpy(flagptr, "recode AD", 10);
	    } else {
	      memcpy(flagptr, "recode A", 9);
	    }
	    printf("Note: --%s flag deprecated.  Use '%s ...'.\n", argptr, flagptr);
	    mm++;
            break;
	  }
	}
	goto main_flag_copy;
      case 'u':
	if (!memcmp(argptr, "update-freq", 12)) {
	  memcpy(flagptr, "read-freq", 10);
	  break;
	} else if (!memcmp(argptr, "update-ref-allele", 18)) {
	  // GCTA alias
	  memcpy(flagptr, "reference-allele", 17);
	  break;
	}
	// fall through
      default:
      main_flag_copy:
	memcpy(flagptr, argptr, kk);
      }
      flagptr = &(flagptr[MAX_FLAG_LEN]);
      flag_map[cur_flag++] = ii;
    }
  }
  sptr = (char*)malloc(flag_ct * MAX_FLAG_LEN);
  if (!sptr) {
    print_ver();
    goto main_ret_NOMEM;
  }
  qsort_ext2(flag_buf, flag_ct, MAX_FLAG_LEN, strcmp_deref, (char*)flag_map, sizeof(int32_t), sptr, MAX_FLAG_LEN);
  free(sptr);
  jj = strlen_se(flag_buf);
  for (cur_flag = 1; cur_flag < flag_ct; cur_flag++) {
    kk = strlen_se(&(flag_buf[cur_flag * MAX_FLAG_LEN]));
    if ((jj == kk) && (!memcmp(&(flag_buf[(cur_flag - 1) * MAX_FLAG_LEN]), &(flag_buf[cur_flag * MAX_FLAG_LEN]), kk))) {
      flag_buf[cur_flag * MAX_FLAG_LEN + kk] = '\0'; // just in case of aliases
      print_ver();
      printf("Error: Duplicate --%s flag.\n", &(flag_buf[cur_flag * MAX_FLAG_LEN]));
      goto main_ret_INVALID_CMDLINE;
    }
    jj = kk;
  }

  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    if (!memcmp("silent", &(flag_buf[cur_flag * MAX_FLAG_LEN]), 7)) {
      freopen("/dev/null", "w", stdout);
      break;
    }
  }
  print_ver();
  uii = 5;
  memcpy(outname, "wdist", 6);
  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    ii = memcmp("out", &(flag_buf[cur_flag * MAX_FLAG_LEN]), 4);
    if (!ii) {
      ii = flag_map[cur_flag];
      if (enforce_param_ct_range(argc, argv, ii, 1, 1, &jj)) {
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      if (strlen(argv[ii + 1]) > (FNAMESIZE - MAX_POST_EXT)) {
	fputs("Error: --out parameter too long.\n", stdout);
	goto main_ret_OPEN_FAIL;
      }
      uii = strlen(argv[ii + 1]);
      memcpy(outname, argv[ii + 1], uii + 1);
      outname_end = &(outname[uii]);
    }
    if (ii <= 0) {
      break;
    }
  }
  memcpy(&(outname[uii]), ".log", 5);
  if (fopen_checked(&logfile, outname, "w")) {
    goto main_ret_OPEN_FAIL;
  }
  printf("Logging to %s.\n", outname);
  outname[uii] = '\0';

  logstr(ver_str);
  sprintf(logbuf, "\n%d argument%s:", argc + mm - cur_arg, (argc + mm - cur_arg == 1)? "" : "s");
  logstr(logbuf);
  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    logstr(" --");
    logstr(&(flag_buf[cur_flag * MAX_FLAG_LEN]));
    ii = flag_map[cur_flag] + 1;
    while ((ii < argc) && (!is_flag(argv[ii]))) {
      logstr(" ");
      logstr(argv[ii++]);
    }
  }
#if _WIN32
  windows_dw = 4 * MAXLINELEN + 256;
  if (GetComputerName(tbuf, &windows_dw))
#else
  if (gethostname(tbuf, 4 * MAXLINELEN + 256) != -1)
#endif
  {
    logstr("\nHostname: ");
    logstr(tbuf);
  }
  logstr("\nWorking directory: ");
  getcwd(tbuf, FNAMESIZE);
  logstr(tbuf);
  logstr("\nStart time: ");
  time(&rawtime);
  logstr(ctime(&rawtime));
  logstr("\n");

  chrom_info.species = SPECIES_HUMAN;
  chrom_info.chrom_mask = 0;
#if _WIN32
  GetSystemInfo(&sysinfo);
  g_thread_ct = sysinfo.dwNumberOfProcessors;
#else
  ii = sysconf(_SC_NPROCESSORS_ONLN);
  if (ii == -1) {
    g_thread_ct = 1;
  } else if (ii > MAX_THREADS) {
    g_thread_ct = MAX_THREADS;
  } else {
    g_thread_ct = ii;
  }
#endif
  strcpy(mapname, "wdist.map");
  strcpy(pedname, "wdist.ped");
  famname[0] = '\0';
  genname[0] = '\0';
  samplename[0] = '\0';
  memcpy(output_missing_pheno, "-9", 3);
  // stuff that must be processed before regular alphabetical loop
  cur_flag = 0;
  if (flag_match("cow", &cur_flag, flag_ct, flag_buf)) {
    if (species_flag(&chrom_info, SPECIES_COW)) {
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("dog", &cur_flag, flag_ct, flag_buf)) {
    if (species_flag(&chrom_info, SPECIES_DOG)) {
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("horse", &cur_flag, flag_ct, flag_buf)) {
    if (species_flag(&chrom_info, SPECIES_HORSE)) {
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("mouse", &cur_flag, flag_ct, flag_buf)) {
    if (species_flag(&chrom_info, SPECIES_MOUSE)) {
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("rice", &cur_flag, flag_ct, flag_buf)) {
    if (species_flag(&chrom_info, SPECIES_RICE)) {
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("sheep", &cur_flag, flag_ct, flag_buf)) {
    if (species_flag(&chrom_info, SPECIES_SHEEP)) {
      goto main_ret_INVALID_CMDLINE;
    }
  }
  species_singular = species_singulars[chrom_info.species];
  species_plural = species_plurals[chrom_info.species];
  cur_flag = 0;
  do {
    argptr = &(flag_buf[cur_flag * MAX_FLAG_LEN]);
    if (!(*argptr)) {
      // preprocessed
      continue;
    }
    argptr2 = &(argptr[1]);
    cur_arg = flag_map[cur_flag];
    switch (*argptr) {
    case '1':
      if (*argptr2 == '\0') {
	affection_01 = 1;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'a':
      if (!memcmp(argptr2, "utosome", 8)) {
	chrom_info.chrom_mask = species_autosome_mask[chrom_info.species];
      } else if (!memcmp(argptr2, "utosome-xy", 11)) {
	if (chrom_info.chrom_mask) {
	  logprint("Error: --autosome-xy and --autosome cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	ii = species_xy_code[chrom_info.species];
	if (ii == -1) {
	  sprintf(logbuf, "Error: --autosome-xy used with a species lacking an XY region.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	chrom_info.chrom_mask = species_autosome_mask[chrom_info.species] | (1LLU << ii);
      } else if (!memcmp(argptr2, "llow-no-sex", 12)) {
	logprint("Warning: --allow-no-sex currently has no effect.\n");
	allow_no_sex = 1;
      } else if (!memcmp(argptr2, "ll", 3)) {
	logprint("Note: --all flag has no effect.\n");
      } else if (!memcmp(argptr2, "llele1234", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii == 1) {
	  if (memcmp("multichar", argv[cur_arg + 1], 10)) {
	    sprintf(logbuf, "Error: Invalid --allele1234 parameter '%s'.%s\n", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  allelexxxx = 2;
	} else {
	  allelexxxx = 1;
	}
      } else if (!memcmp(argptr2, "lleleACGT", 9)) {
	if (allelexxxx) {
	  logprint("Error: --allele1234 and --alleleACGT cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii == 1) {
	  if (memcmp("multichar", argv[cur_arg + 1], 10)) {
	    sprintf(logbuf, "Error: Invalid --alleleACGT parameter '%s'.%s\n", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  allelexxxx = 4;
	} else {
	  allelexxxx = 3;
	}
      } else if (!memcmp(argptr2, "llele-count", 12)) {
	lgen_modifier |= LGEN_ALLELE_COUNT;
      } else if (!memcmp(argptr2, "ll-pheno", 9)) {
	pheno_modifier |= PHENO_ALL;
      } else if (!memcmp(argptr2, "ssoc", 5)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 5, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (jj = 1; jj <= ii; jj++) {
	  if (!strcmp(argv[cur_arg + jj], "counts")) {
	    if (model_modifier & MODEL_QMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means' modifier does not make sense with 'counts'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_ASSOC_COUNTS;
	  } else if (!strcmp(argv[cur_arg + jj], "fisher")) {
	    if (model_modifier & MODEL_QMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means' modifier does not make sense with 'fisher'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_FISHER;
	  } else if (!strcmp(argv[cur_arg + jj], "perm")) {
	    model_modifier |= MODEL_PERM;
	  } else if (!strcmp(argv[cur_arg + jj], "genedrop")) {
	    if (model_modifier & MODEL_QMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means' modifier does not make sense with 'genedrop'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + jj], "perm-count")) {
	    model_modifier |= MODEL_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + jj], "p2")) {
	    model_modifier |= MODEL_ASSOC_P2;
	  } else if ((!memcmp(argv[cur_arg + jj], "mperm=", 6)) && (argv[cur_arg + jj][6] != '\0')) {
	    if (model_modifier & MODEL_PERM) {
	      sprintf(logbuf, "Error: --assoc 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: Duplicate --assoc 'mperm' modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    kk = atoi(&(argv[cur_arg + jj][6]));
	    if (kk < 1) {
	      sprintf(logbuf, "Error: Invalid --assoc mperm parameter '%s'.%s", &(argv[cur_arg + jj][6]), errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_mperm_val = (uint32_t)kk;
	  } else if (!strcmp(argv[cur_arg + jj], "qt-means")) {
	    if (model_modifier & MODEL_DMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means' does not make sense with a case/control-specific\nmodifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_QT_MEANS;
	  } else {
	    sprintf(logbuf, "Error: Invalid --assoc parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	model_modifier |= MODEL_ASSOC;
	calculation_type |= CALC_MODEL;
      } else if (!memcmp(argptr2, "djust", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 3, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mtest_adjust = 1;
	for (jj = 1; jj <= ii; jj++) {
	  if (!memcmp(argv[cur_arg + jj], "gc", 3)) {
	    mtest_adjust |= ADJUST_GC;
	  } else if (!memcmp(argv[cur_arg + jj], "log10", 6)) {
	    mtest_adjust |= ADJUST_LOG10;
	  } else if (!memcmp(argv[cur_arg + jj], "qq-plot", 8)) {
	    mtest_adjust |= ADJUST_QQ;
	  } else {
	    sprintf(logbuf, "Error: Invalid --adjust parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
      } else if (!memcmp(argptr2, "perm", 5)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 6, 6, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --aperm min permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	aperm_min = ii;
	ii = atoi(argv[cur_arg + 2]);
	if (ii < ((int32_t)aperm_min)) {
	  sprintf(logbuf, "Error: Invalid --aperm max permutation count '%s'.%s", argv[cur_arg + 2], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	aperm_max = ii;
	if (sscanf(argv[cur_arg + 3], "%lg", &aperm_alpha) != 1) {
	  sprintf(logbuf, "Error: Invalid --aperm alpha threshold '%s'.%s", argv[cur_arg + 3], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + 4], "%lg", &aperm_beta) != 1) {
	  sprintf(logbuf, "Error: Invalid --aperm beta '%s'.%s", argv[cur_arg + 4], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + 5], "%lg", &aperm_init_interval) != 1) {
	  sprintf(logbuf, "Error: Invalid --aperm initial pruning interval '%s'.%s", argv[cur_arg + 5], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (aperm_init_interval < 1) {
	  sprintf(logbuf, "Error: Invalid --aperm initial pruning interval '%s'.%s", argv[cur_arg + 5], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + 6], "%lg", &aperm_interval_slope) != 1) {
	  sprintf(logbuf, "Error: Invalid --aperm pruning interval slope '%s'.%s", argv[cur_arg + 6], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (aperm_interval_slope < 0) {
	  sprintf(logbuf, "Error: Invalid --aperm pruning interval slope '%s'.%s", argv[cur_arg + 6], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'b':
      if (!memcmp(argptr2, "file", 5)) {
	load_params |= 8;
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 5)) {
	    logprint("Error: --bfile parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  sptr = (char*)"wdist";
	}
	if (!(load_params & 16)) {
	  strcpy(pedname, sptr);
	  strcat(pedname, ".bed");
	}
	if (!(load_params & 32)) {
	  strcpy(mapname, sptr);
	  strcat(mapname, ".bim");
	}
	if (!(load_params & 64)) {
	  strcpy(famname, sptr);
	  strcat(famname, ".fam");
	}
      } else if (!memcmp(argptr2, "ed", 3)) {
	load_params |= 16;
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --bed parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "im", 3)) {
	load_params |= 32;
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --bim parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "merge", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 3, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii == 2) {
	  sprintf(logbuf, "Error: --bmerge must have exactly 1 or 3 parameters.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (ii == 3) {
	  if (++jj > FNAMESIZE) {
	    logprint("Error: --bmerge .bed filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  jj = strlen(argv[cur_arg + 2]) + 1;
	  if (jj > FNAMESIZE) {
	    logprint("Error: --bmerge .bim filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename2, argv[cur_arg + 2], jj);
	  jj = strlen(argv[cur_arg + 3]) + 1;
	  if (jj > FNAMESIZE) {
	    logprint("Error: --bmerge .fam filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename3, argv[cur_arg + 3], jj);
	} else {
	  if (jj > (FNAMESIZE - 5)) {
	    logprint("Error: --bmerge filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  memcpy(mergename2, argv[cur_arg + 1], jj);
	  memcpy(mergename3, argv[cur_arg + 1], jj);
	  memcpy(&(mergename1[jj]), ".bed", 5);
	  memcpy(&(mergename2[jj]), ".bim", 5);
	  memcpy(&(mergename3[jj]), ".fam", 5);
	}
	calculation_type |= CALC_MERGE;
	merge_type |= MERGE_BINARY;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'c':
      if (!memcmp(argptr2, "hr", 3)) {
	if (chrom_info.chrom_mask) {
	  sprintf(logbuf, "Error: --chr cannot be used with --autosome(-xy).%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (parse_chrom_ranges(param_count(argc, argv, cur_arg), '-', &(argv[cur_arg]), &(chrom_info.chrom_mask), chrom_info.species, argptr)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!chrom_info.chrom_mask) {
	  sprintf(logbuf, "Error: All chromosomes excluded.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "hr-excl", 8)) {
	if (parse_chrom_ranges(param_count(argc, argv, cur_arg), '-', &(argv[cur_arg]), &chrom_exclude, chrom_info.species, argptr)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!chrom_info.chrom_mask) {
	  chrom_info.chrom_mask = species_valid_chrom_mask[chrom_info.species] & (~chrom_exclude);
	} else {
	  chrom_info.chrom_mask &= ~chrom_exclude;
	}
	if (!chrom_info.chrom_mask) {
	  sprintf(logbuf, "Error: All chromosomes excluded.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "ompound-genotypes", 18)) {
	logprint("Note: --compound-genotypes flag unnecessary (spaces between alleles in .ped\nare optional).\n");
      } else if (!memcmp(argptr2, "ompress", 8)) {
	logprint("Error: --compress flag retired.  Use e.g. 'gzip [filename]'.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ounts", 6)) {
	if (model_modifier & MODEL_ASSOC) {
	  if (model_modifier & MODEL_QMASK) {
	    sprintf(logbuf, "Error: --assoc 'qt-means' modifier does not make sense with --counts.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  logprint("Note: --counts flag deprecated.  Use '--assoc counts' instead.\n");
          model_modifier |= MODEL_ASSOC_COUNTS;
	} else {
	  logprint("Note: --counts flag deprecated.  Use '--freq counts' or --freqx instead.\n");
	}
        freq_counts = 1;
      } else if (!memcmp(argptr2, "ovar", 5)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Error: --covar has not been implemented yet.\n");
	goto main_ret_1;
	retval = alloc_fname(&covar_fname, argv[cur_arg + 1], argptr, 0);
	if (!retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ovar-name", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!covar_fname) {
	  logprint("Error: --covar-name must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (alloc_string(&covar_str, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
	covar_modifier = COVAR_NAME;
      } else if (!memcmp(argptr2, "ovar-number", 12)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (covar_modifier) {
	  logprint("Error: --covar-number cannot be used with --covar-name.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (!covar_fname) {
	  logprint("Error: --covar-number must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (alloc_string(&covar_str, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
	covar_modifier = COVAR_NUMBER;
      } else if (!memcmp(argptr2, "ell", 4)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &model_cell_ct)) {
	  sprintf(logbuf, "Error: Invalid --cell parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "i", 2)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (sscanf(argv[cur_arg + 1], "%lg", &dxx) != 1) {
	  sprintf(logbuf, "Error: Invalid --ci parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx < 0.01) || (dxx >= 1.0)) {
	  sprintf(logbuf, "Error: --ci confidence interval size s must satisfy 0.01 <= s < 1.%s", errstr_append);
	}
	ci_size = dxx;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'd':
      if (!(*argptr2)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (argv[cur_arg + 1][1]) {
	  sprintf(logbuf, "Error: --d parameter too long (must be a single character).%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	range_delim = argv[cur_arg + 1][0];
	if ((range_delim == '-') || (range_delim == ',')) {
	  sprintf(logbuf, "Error: --d parameter cannot be '-' or ','.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "ata", 4)) {
	if (load_params & 0xff) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 0x80;
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 5)) {
	    logprint("Error: --data parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  sptr = (char*)"wdist";
	}
	if (!(load_params & 0x100)) {
	  strcpy(genname, sptr);
	  strcat(genname, ".gen");
	}
	if (!(load_params & 0x200)) {
	  strcpy(samplename, sptr);
	  strcat(samplename, ".sample");
	}
      } else if (!memcmp(argptr2, "ebug", 5)) {
	debug_on = 1;
      } else if (!memcmp(argptr2, "ecompress", 10)) {
	logprint("Error: --decompress flag retired.  Use e.g. 'gunzip [filename]'.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "istance", 8)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 7, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (jj = 1; jj <= ii; jj++) {
	  if (!memcmp(argv[cur_arg + jj], "square", 7)) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ;
	  } else if (!memcmp(argv[cur_arg + jj], "square0", 8)) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --distance 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ0;
	  } else if (!memcmp(argv[cur_arg + jj], "triangle", 9)) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --distance 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_TRI;
	  } else if (!memcmp(argv[cur_arg + jj], "gz", 3)) {
	    if (dist_calc_type & DISTANCE_BIN) {
	      sprintf(logbuf, "Error: --distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_GZ;
	  } else if (!memcmp(argv[cur_arg + jj], "bin", 4)) {
	    if (dist_calc_type & DISTANCE_GZ) {
	      sprintf(logbuf, "Error: --distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_BIN;
	  } else if (!memcmp(argv[cur_arg + jj], "ibs", 4)) {
	    if (dist_calc_type & DISTANCE_IBS) {
	      logprint("Error: Duplicate --distance 'ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_IBS;
	  } else if (!memcmp(argv[cur_arg + jj], "1-ibs", 6)) {
	    if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	      logprint("Error: Duplicate --distance '1-ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_1_MINUS_IBS;
	  } else if (!memcmp(argv[cur_arg + jj], "alct", 5)) {
	    if (dist_calc_type & DISTANCE_ALCT) {
	      logprint("Error: Duplicate --distance 'alct' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_ALCT;
	  } else if (!memcmp(argv[cur_arg + jj], "3d", 3)) {
	    distance_3d = 1;
	  } else if (!memcmp(argv[cur_arg + jj], "flat-missing", 13)) {
	    distance_flat_missing = 1;
	  } else {
	    sprintf(logbuf, "Error: Invalid --distance parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (!(dist_calc_type & DISTANCE_TYPEMASK)) {
	  dist_calc_type |= DISTANCE_ALCT;
	}
	calculation_type |= CALC_DISTANCE;
      } else if (!memcmp(argptr2, "istance-matrix", 15)) {
	if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	  sprintf(logbuf, "Error: --distance-matrix flag cannot be used with '--distance 1-ibs'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_PLINK_DISTANCE_MATRIX;
      } else if (!memcmp(argptr2, "ummy", 5)) {
	if (load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
        if (enforce_param_ct_range(argc, argv, cur_arg, 2, 6, &ii)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	dummy_indiv_ct = atoi(argv[cur_arg + 1]);
	if (dummy_indiv_ct < 1) {
	  sprintf(logbuf, "Error: Invalid --dummy individual count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dummy_marker_ct = atoi(argv[cur_arg + 2]);
	if (dummy_marker_ct < 1) {
	  sprintf(logbuf, "Error: Invalid --dummy marker count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        for (jj = 3; jj <= ii; jj++) {
	  if ((argv[cur_arg + jj][4] == '\0') && (tolower(argv[cur_arg + jj][0]) == 'a') && (tolower(argv[cur_arg + jj][1]) == 'c') && (tolower(argv[cur_arg + jj][2]) == 'g') && (tolower(argv[cur_arg + jj][3]) == 't')) {
	    if (dummy_flags & (DUMMY_1234 | DUMMY_12)) {
	      sprintf(logbuf, "Error: --dummy 'acgt' modifier cannot be used with '1234' or '12'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            dummy_flags |= DUMMY_ACGT;
	  } else if (!memcmp(argv[cur_arg + jj], "1234", 5)) {
	    if (dummy_flags & (DUMMY_ACGT | DUMMY_12)) {
	      sprintf(logbuf, "Error: --dummy '1234' modifier cannot be used with 'acgt' or '12'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            dummy_flags |= DUMMY_1234;
	  } else if (!memcmp(argv[cur_arg + jj], "12", 3)) {
	    if (dummy_flags & (DUMMY_ACGT | DUMMY_1234)) {
	      sprintf(logbuf, "Error: --dummy '12' modifier cannot be used with 'acgt' or '1234'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            dummy_flags |= DUMMY_12;
	  } else if (!memcmp(argv[cur_arg + jj], "scalar-pheno", 13)) {
	    dummy_flags |= DUMMY_SCALAR_PHENO;
	  } else {
	    if ((dummy_flags & DUMMY_MISSING_PHENO) || (sscanf(argv[cur_arg + jj], "%lg", &dxx) != 1) || (dxx < 0.0) || (dxx > 1.0)) {
	      sprintf(logbuf, "Error: Invalid --dummy parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
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
	  sprintf(logbuf, "Error: --dummy-coding cannot be used without --covar.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	write_covar_modifier |= WRITE_COVAR_DUMMY;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'e':
      if (!memcmp(argptr2, "xtract", 7)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&extractname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "xclude", 7)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&excludename, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "xponent", 8)) {
	if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
	  logprint("Error: --exponent cannot be used with --distance-matrix.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + 1], "%lg", &exponent) != 1) {
	  sprintf(logbuf, "Error: Invalid --exponent parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'f':
      if (!memcmp(argptr2, "ile", 4)) {
	if (load_params & 0x3f9) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 1;
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 5)) {
	    logprint("Error: --file parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  if (!(load_params & 2)) {
	    strcpy(pedname, argv[cur_arg + 1]);
	    strcat(pedname, ".ped");
	  }
	  if (!(load_params & 4)) {
	    strcpy(mapname, argv[cur_arg + 1]);
	    strcat(mapname, ".map");
	  }
	}
      } else if (!memcmp(argptr2, "am", 3)) {
	if (load_params & 0x3c7) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 64;
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --fam parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(famname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "ilter", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 2, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&filtername, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (alloc_string(&filterval, argv[cur_arg + 2])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ilter-cases", 12)) {
	filter_case_control = 1;
      } else if (!memcmp(argptr2, "ilter-controls", 15)) {
	if (filter_case_control == 1) {
	  logprint("Error: --filter-cases and --filter-controls cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_case_control = 2;
      } else if (!memcmp(argptr2, "ilter-females", 14)) {
	filter_sex = 2;
      } else if (!memcmp(argptr2, "ilter-males", 12)) {
	if (filter_sex == 2) {
	  logprint("Error: --filter-males and --filter-females cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_sex = 1;
      } else if (!memcmp(argptr2, "ilter-founders", 15)) {
	filter_founder_nonf = 1;
      } else if (!memcmp(argptr2, "ilter-nonfounders", 18)) {
	if (filter_founder_nonf == 1) {
	  logprint("Error: --filter-founders and --filter-nonfounders cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_founder_nonf = 2;
      } else if (!memcmp(argptr2, "req", 4)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  if (memcmp(argv[cur_arg + 1], "counts", 7)) {
            sprintf(logbuf, "Error: Invalid --freq parameter '%s'.%s\n", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  freq_counts = 1;
	}
	calculation_type |= CALC_FREQ;
      } else if (!memcmp(argptr2, "reqx", 5)) {
	if (calculation_type & CALC_FREQ) {
	  sprintf(logbuf, "Error: --freqx cannot be used with --freq.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_FREQ;
	freqx = 1;
      } else if (!memcmp(argptr2, "rom", 4)) {
	if (chrom_info.chrom_mask) {
	  sprintf(logbuf, "Error: --from cannot be used with --autosome(-xy) or --chr%s.%s", chrom_exclude? "-excl" : "", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_from, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if ((!memcmp(argptr2, "rom-bp", 7)) || (!memcmp(argptr2, "rom-kb", 7)) || (!memcmp(argptr2, "rom-mb", 7))) {
	if (markername_from) {
	  sprintf(logbuf, "Error: --from-bp/-kb/-mb cannot be used with --from.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argptr2[4];
	if (cc == 'b') {
	  if (atoiz(argv[cur_arg + 1], &marker_pos_start)) {
	    sprintf(logbuf, "Error: Invalid --from-bp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  if (marker_pos_start != -1) {
	    logprint("Error: Multiple --from-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (sscanf(argv[cur_arg + 1], "%lg", &dxx) != 1) {
	    sprintf(logbuf, "Error: Invalid --from-kb/-mb parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  dxx *= (cc == 'k')? 1000 : 1000000;
	  if (dxx < 0) {
	    marker_pos_start = 0;
	  } else if (dxx > 2147483647) {
	    marker_pos_start = 2147483647;
	  } else {
	    marker_pos_start = (int)dxx;
	  }
	}
      } else if (!memcmp(argptr2, "isher", 6)) {
	if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --fisher cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --fisher flag deprecated.  Use '--assoc fisher' or '--model fisher'.\n");
	model_modifier |= MODEL_ASSOC | MODEL_FISHER | MODEL_ASSOC_FDEPR;
	calculation_type = CALC_MODEL;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'g':
      if (!memcmp(argptr2, "en", 3)) {
	if (load_params & 0x17f) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 0x100;
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --gen parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(genname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "eno", 4)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  if (sscanf(argv[cur_arg + 1], "%lg", &geno_thresh) != 1) {
	    sprintf(logbuf, "Error: Invalid --geno parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((geno_thresh < 0.0) || (geno_thresh > 1.0)) {
	    sprintf(logbuf, "Error: Invalid --geno parameter '%s' (must be between 0 and 1).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  geno_thresh = 0.1;
	}
      } else if ((!memcmp(argptr2, "enome", 6)) || (!memcmp(argptr2, "enome gz", 9))) {
	if (argptr2[5] == ' ') {
	  kk = 1;
	  genome_output_gz = 1;
	} else {
	  kk = 0;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 3 - kk, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (jj = 1; jj <= ii; jj++) {
	  if (!strcmp(argv[cur_arg + jj], "full")) {
	    genome_output_full = 1;
	  } else if (!strcmp(argv[cur_arg + jj], "unbounded")) {
	    genome_ibd_unbounded = 1;
	  } else if (!strcmp(argv[cur_arg + jj], "gz")) {
	    genome_output_gz = 1;
	  } else {
	    sprintf(logbuf, "Error: Invalid --genome parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_GENOME;
      } else if (!memcmp(argptr2, "roupdist", 9)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  groupdist_iters = strtoul(argv[cur_arg + 1], NULL, 10);
	  if ((groupdist_iters < 2) || (groupdist_iters == ULONG_MAX)) {
	    sprintf(logbuf, "Error: Invalid --groupdist jackknife iteration count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (ii == 2) {
	    groupdist_d = atoi(argv[cur_arg + 2]);
	    if (groupdist_d <= 0) {
	      sprintf(logbuf, "Error: Invalid --groupdist jackknife delete parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	calculation_type |= CALC_GROUPDIST;
      } else if (!memcmp(argptr2, "rm", 3)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (ii) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 8)) {
	    logprint("Error: --grm parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, sptr);
	} else {
	  memcpy(pedname, "wdist", 6);
	}
        load_rare = LOAD_RARE_GRM;
      } else if (!memcmp(argptr2, "xe", 3)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!covar_fname) {
	  logprint("Error: --gxe must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (ii) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --gxe parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  gxe_mcovar = ii;
	} else {
	  gxe_mcovar = 1;
	}
	calculation_type |= CALC_GXE;
      } else if (!memcmp(argptr2, "enedrop", 8)) {
	if (model_modifier & MODEL_QMASK) {
	  sprintf(logbuf, "Error: --assoc 'qt-means' modifier does not make sense with --genedrop.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --genedrop flag deprecated.  Use e.g. '--model genedrop'.\n");
	model_modifier |= MODEL_GENEDROP;
      } else if (!memcmp(argptr2, "c", 2)) {
        if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --gc must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --gc flag deprecated.  Use '--adjust gc'.\n");
	mtest_adjust |= ADJUST_GC;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'h':
      if (!memcmp(argptr2, "we", 3)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  if (sscanf(argv[cur_arg + 1], "%lg", &hwe_thresh) != 1) {
	    sprintf(logbuf, "Error: Invalid --hwe parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((hwe_thresh < 0.0) || (hwe_thresh >= 1.0)) {
	    sprintf(logbuf, "Error: Invalid --hwe parameter '%s' (must be between 0 and 1).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  hwe_thresh = 0.001;
	}
      } else if (!memcmp(argptr2, "we-all", 7)) {
	hwe_all = 1;
      } else if (!memcmp(argptr2, "et", 3)) {
	sprintf(logbuf, "Error: --het retired.  Use --ibc.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if (!memcmp(argptr2, "ardy", 5)) {
	calculation_type |= CALC_HARDY;
      } else if (!memcmp(argptr2, "we2", 4)) {
	sprintf(logbuf, "Error: --hwe2 retired.  Use the --hwe exact test.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if (!memcmp(argptr2, "ardy2", 6)) {
	sprintf(logbuf, "Error: --hardy2 retired.  Use the exact test-based --hardy report.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'i':
      if (!memcmp(argptr2, "bc", 3)) {
	calculation_type |= CALC_IBC;
      } else if (!memcmp(argptr2, "ndep-pairwise", 14)) {
	if (calculation_type & CALC_LD_PRUNE) {
	  logprint("Error: --indep-pairwise cannot be used with --indep.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 3, 4, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_window_size = atoi(argv[cur_arg + 1]);
	if ((ld_window_size < 1) || ((ld_window_size == 1) && (ii == 3))) {
	  sprintf(logbuf, "Error: Invalid --indep-pairwise window size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii == 4) {
	  if ((tolower(argv[cur_arg + 2][0]) != 'k') || (tolower(argv[cur_arg + 2][1]) != 'b') || (argv[cur_arg + 2][2] != '\0')) {
	    sprintf(logbuf, "Error: Invalid --indep-pairwise parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  ld_window_kb = 1;
	} else {
	  jj = strlen(argv[cur_arg + 1]);
	  if ((jj > 2) && (tolower(argv[cur_arg + 1][jj - 2]) == 'k') && (tolower(argv[cur_arg + 1][jj - 1]) == 'b')) {
	    ld_window_kb = 1;
	  }
	}
	ld_window_incr = atoi(argv[cur_arg + ii - 1]);
	if (ld_window_incr < 1) {
	  sprintf(logbuf, "Error: Invalid increment '%s' for --indep-pairwise.%s", argv[cur_arg + ii - 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + ii], "%lg", &ld_last_param) != 1) {
	  sprintf(logbuf, "Error: Invalid --indep-pairwise r^2 threshold '%s'.%s", argv[cur_arg + ii], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((ld_last_param < 0.0) || (ld_last_param >= 1.0)) {
	  sprintf(logbuf, "Error: Invalid --indep-pairwise r^2 threshold '%s'.%s", argv[cur_arg + ii], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= (CALC_LD_PRUNE | CALC_LD_PRUNE_PAIRWISE);
      } else if (!memcmp(argptr2, "ndep", 5)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 3, 4, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_window_size = atoi(argv[cur_arg + 1]);
	if ((ld_window_size < 1) || ((ld_window_size == 1) && (ii == 3))) {
	  sprintf(logbuf, "Error: Invalid --indep window size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii == 4) {
	  if ((tolower(argv[cur_arg + 2][0]) != 'k') || (tolower(argv[cur_arg + 2][1]) != 'b') || (argv[cur_arg + 2][2] != '\0')) {
	    sprintf(logbuf, "Error: Invalid --indep parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  ld_window_kb = 1;
	} else {
	  jj = strlen(argv[cur_arg + 1]);
	  if ((jj > 2) && (tolower(argv[cur_arg + 1][jj - 2]) == 'k') && (tolower(argv[cur_arg + 1][jj - 1]) == 'b')) {
	    ld_window_kb = 1;
	  }
	}
	ld_window_incr = atoi(argv[cur_arg + ii - 1]);
	if (ld_window_incr < 1) {
	  sprintf(logbuf, "Error: Invalid increment '%s' for --indep.%s", argv[cur_arg + ii - 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + ii], "%lg", &ld_last_param) != 1) {
	  sprintf(logbuf, "Error: Invalid --indep VIF threshold '%s'.%s", argv[cur_arg + ii], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ld_last_param < 1.0) {
	  sprintf(logbuf, "Error: --indep VIF threshold '%s' too small (must be >= 1).%s", argv[cur_arg + ii], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_LD_PRUNE;
      } else if (!memcmp(argptr2, "ndiv-sort", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = (argv[cur_arg + 1][1] == '\0');
        if (!memcmp(argv[cur_arg + 1], "none", 5) || ((argv[cur_arg + 1][0] == '0') && jj)) {
	  indiv_sort = INDIV_SORT_NONE;
	} else if (!memcmp(argv[cur_arg + 1], "natural", 8) || ((tolower(argv[cur_arg + 1][0]) == 'n') && jj)) {
	  indiv_sort = INDIV_SORT_NATURAL;
	} else if ((!memcmp(argv[cur_arg + 1], "ascii", 6)) || ((tolower(argv[cur_arg + 1][0]) == 'a') && jj)) {
	  indiv_sort = INDIV_SORT_ASCII;
	} else {
	  sprintf(logbuf, "Error: '%s' is not a valid mode for --indiv-sort.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'k':
      if (!memcmp(argptr2, "eep", 4)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&keepname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eep-allele-order", 17)) {
	keep_allele_order = 1;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'l':
      if (!memcmp(argptr2, "oad-dists", 10)) {
	if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
	  sprintf(logbuf, "Error: --load-dists cannot be used with --distance-matrix.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&loaddistname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	calculation_type |= CALC_LOAD_DISTANCES;
      } else if (!memcmp(argptr2, "file", 5)) {
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 6) {
	    logprint("Error: --lfile filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, argv[cur_arg + 1]);
	} else {
	  memcpy(pedname, "wdist", 6);
	}
	load_rare = LOAD_RARE_LGEN;
      } else if (!memcmp(argptr2, "oop-assoc", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&loop_assoc_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "og10", 5)) {
        if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --log10 must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --log10 flag deprecated.  Use '--adjust log10'.\n");
	mtest_adjust |= ADJUST_LOG10;
      } else if (!memcmp(argptr2, "ambda", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --lambda must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + 1], "%lg", &adjust_lambda) != 1) {
	  sprintf(logbuf, "Error: Invalid --lambda parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (adjust_lambda < 1) {
	  logprint("Note: --lambda parameter set to 1.\n");
	  adjust_lambda = 1;
	}
	mtest_adjust |= ADJUST_LAMBDA;
      } else {
        goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'm':
      if (!memcmp(argptr2, "ap", 3)) {
	if ((load_params & 0x3fc) || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 4;
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --map parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "issing-genotype", 16)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	missing_geno = argv[cur_arg + 1][0];
	if ((argv[cur_arg + 1][1] != '\0') || (((unsigned char)missing_geno) <= ' ')) {
	  sprintf(logbuf, "Error: Invalid --missing-genotype parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "issing-phenotype", 17)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	missing_pheno = atoi(argv[cur_arg + 1]);
	if ((missing_pheno == 0) || (missing_pheno == 1)) {
	  sprintf(logbuf, "Error: Invalid --missing-phenotype parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if ((!memcmp(argptr2, "issing-code", 12))) {
        if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  missing_code = argv[cur_arg + 1];
	} else {
	  missing_code = (char*)"";
	}
      } else if (!memcmp(argptr2, "ake-pheno", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 2, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&phenoname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (alloc_string(&makepheno_str, argv[cur_arg + 2])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "pheno", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mpheno_col = atoi(argv[cur_arg + 1]);
	if (mpheno_col < 1) {
	  sprintf(logbuf, "Error: Invalid --mpheno parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "filter", 7)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mfilter_col = atoi(argv[cur_arg + 1]);
	if (mfilter_col < 1) {
	  sprintf(logbuf, "Error: Invalid --mfilter parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "emory", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	malloc_size_mb = atoi(argv[cur_arg + 1]);
	if (malloc_size_mb < WKSPACE_MIN_MB) {
	  if (malloc_size_mb > 0) {
	    sprintf(logbuf, "Error: Invalid --memory parameter '%s' (minimum %u).%s", argv[cur_arg + 1], WKSPACE_MIN_MB, errstr_append);
	  } else {
	    sprintf(logbuf, "Error: Invalid --memory parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  }
	  goto main_ret_INVALID_CMDLINE_3;
	}
#ifndef __LP64__
	if (malloc_size_mb > 2944) {
	  logprint("Error: --memory parameter too large for 32-bit version (max 2944).\n");
	  goto main_ret_INVALID_CMDLINE;
	}
#endif
      } else if (!memcmp(argptr2, "af", 3)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  if (sscanf(argv[cur_arg + 1], "%lg", &min_maf) != 1) {
	    sprintf(logbuf, "Error: Invalid --maf parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (min_maf <= 0.0) {
	    sprintf(logbuf, "Error: --maf parameter '%s' too small (must be > 0).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else if (min_maf > max_maf) {
	    sprintf(logbuf, "Error: --maf parameter '%s' too large (must be <= %g).%s", argv[cur_arg + 1], max_maf, errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  min_maf = 0.01;
	}
      } else if (!memcmp(argptr2, "ax-maf", 7)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + 1], "%lg", &max_maf) != 1) {
	  sprintf(logbuf, "Error: Invalid --max-maf parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (max_maf < min_maf) {
	  sprintf(logbuf, "Error: --max-maf parameter '%s' too small (must be >= %g).%s", argv[cur_arg + 1], min_maf, errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (max_maf >= 0.5) {
	  sprintf(logbuf, "Error: --max-maf parameter '%s' too large (must be < 0.5).%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "ind", 4)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  if (sscanf(argv[cur_arg + 1], "%lg", &mind_thresh) != 1) {
	    sprintf(logbuf, "Error: Invalid --mind parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((mind_thresh < 0.0) || (mind_thresh > 1.0)) {
	    sprintf(logbuf, "Error: Invalid --mind parameter '%s' (must be between 0 and 1).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  mind_thresh = 0.1;
	}
      } else if (!memcmp(argptr2, "ake-grm", 8)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	rel_calc_type |= REL_CALC_GZ | REL_CALC_GRM;
	for (jj = 1; jj <= ii; jj++) {
	  if (!memcmp(argv[cur_arg + jj], "cov", 4)) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-grm 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-grm 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if (!memcmp(argv[cur_arg + jj], "no-gz", 6)) {
	    rel_calc_type &= ~REL_CALC_GZ;
	  } else if ((!memcmp(argv[cur_arg + jj], "ibc1", 5)) || (!memcmp(argv[cur_arg + jj], "ibc2", 5))) {
	    if (rel_calc_type & REL_CALC_COV) {
	      sprintf(logbuf, "Error: --make-grm 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-grm '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + jj], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ibc_type = argv[cur_arg + jj][3] - '0';
	  } else if (!memcmp(argv[cur_arg + jj], "single-prec", 12)) {
	    rel_calc_type |= REL_CALC_SINGLE_PREC;
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-grm parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "ake-rel", 8)) {
	if (calculation_type & CALC_RELATIONSHIP) {
	  sprintf(logbuf, "Error: --make-rel cannot be used with --make-grm.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 3, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (jj = 1; jj <= ii; jj++) {
	  if (!memcmp(argv[cur_arg + jj], "cov", 4)) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if (!memcmp(argv[cur_arg + jj], "gz", 3)) {
	    if (rel_calc_type & REL_CALC_BIN) {
	      sprintf(logbuf, "Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_GZ;
	  } else if (!memcmp(argv[cur_arg + jj], "bin", 4)) {
	    if (rel_calc_type & REL_CALC_GZ) {
	      sprintf(logbuf, "Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_BIN;
	  } else if (!memcmp(argv[cur_arg + jj], "square", 7)) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_SQ;
	  } else if (!memcmp(argv[cur_arg + jj], "square0", 8)) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      sprintf(logbuf, "Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_SQ0;
	  } else if (!memcmp(argv[cur_arg + jj], "triangle", 9)) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      sprintf(logbuf, "Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_TRI;
	  } else if ((!memcmp(argv[cur_arg + jj], "ibc1", 5)) || (!memcmp(argv[cur_arg + jj], "ibc2", 5))) {
	    if (rel_calc_type & REL_CALC_COV) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-rel '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + jj], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ibc_type = argv[cur_arg + jj][3] - '0';
	  } else if (!memcmp(argv[cur_arg + jj], "single-prec", 12)) {
	    rel_calc_type |= REL_CALC_SINGLE_PREC;
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-rel parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (!(rel_calc_type & REL_CALC_SHAPEMASK)) {
	  rel_calc_type |= (rel_calc_type & REL_CALC_BIN)? REL_CALC_SQ : REL_CALC_TRI;
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "atrix", 6)) {
	if (calculation_type & CALC_LOAD_DISTANCES) {
	  sprintf(logbuf, "Error: --matrix cannot be used with --load-dists.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (exponent != 0.0) {
	  sprintf(logbuf, "Error: --matrix cannot be used with --exponent.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (dist_calc_type & DISTANCE_IBS) {
	  sprintf(logbuf, "Error: --matrix cannot be used with '--distance ibs'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_PLINK_IBS_MATRIX;
      } else if (!memcmp(argptr2, "af-succ", 8)) {
	maf_succ = 1;
      } else if (!memcmp(argptr2, "ap3", 4)) {
	logprint("Note: --map3 flag unnecessary (.map file format is autodetected).\n");
      } else if (!memcmp(argptr2, "ake-bed", 8)) {
	if (no_params_check(argc, argv, cur_arg, argptr, outname_end)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_MAKE_BED;
      } else if (!memcmp(argptr2, "erge", 5)) {
	if (calculation_type & CALC_MERGE) {
	  sprintf(logbuf, "Error: --merge cannot be used with --bmerge.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (ii == 2) {
	  if (++jj > FNAMESIZE) {
	    logprint("Error: --merge .ped filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  jj = strlen(argv[cur_arg + 2]) + 1;
	  if (jj > FNAMESIZE) {
	    logprint("Error: --merge .map filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename2, argv[cur_arg + 2], jj);
	} else {
	  if (jj > (FNAMESIZE - 5)) {
	    logprint("Error: --merge filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  memcpy(mergename2, argv[cur_arg + 1], jj);
	  memcpy(&(mergename1[jj]), ".ped", 5);
	  memcpy(&(mergename2[jj]), ".map", 5);
	}
	calculation_type |= CALC_MERGE;
      } else if (!memcmp(argptr2, "erge-list", 10)) {
	if (calculation_type & CALC_MERGE) {
	  logprint("Error: --merge-list cannot be used with --merge or --bmerge.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]) + 1;
	if (jj > FNAMESIZE) {
	  logprint("Error: --merge-list filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(mergename1, argv[cur_arg + 1], jj);
	merge_type |= MERGE_LIST;
	calculation_type |= CALC_MERGE;
      } else if (!memcmp(argptr2, "erge-mode", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argv[cur_arg + 1][0];
	if ((cc < '1') || (cc > '7') || (argv[cur_arg + 1][1] != '\0')) {
          sprintf(logbuf, "Error: Invalid --merge-mode parameter %s.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((merge_type & MERGE_LIST) && (cc > '5')) {
	  sprintf(logbuf, "Error: --merge-mode 6-7 cannot be used with --merge-list.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        merge_type |= cc - '0';
      } else if (!memcmp(argptr2, "erge-allow-equal-pos", 21)) {
	merge_type |= MERGE_ALLOW_EQUAL_POS;
      } else if (!memcmp(argptr2, "ust-have-sex", 13)) {
	logprint("Warning: --must-have-sex currently has no effect.\n");
	must_have_sex = 1;
      } else if (!memcmp(argptr2, "covar", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (!(calculation_type & CALC_GXE)) {
	  logprint("Error: --mcovar must be used with --covar and --gxe.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	mcovar_col = atoi(argv[cur_arg + 1]);
	if (mcovar_col < 1) {
	  sprintf(logbuf, "Error: Invalid --mcovar parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "odel", 5)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 4, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (model_modifier & MODEL_ASSOC_FDEPR) {
	  model_modifier &= ~(MODEL_ASSOC | MODEL_ASSOC_FDEPR);
	} else if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --model cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	for (jj = 1; jj <= ii; jj++) {
	  if (!strcmp(argv[cur_arg + jj], "fisher")) {
	    if (model_modifier & MODEL_TRENDONLY) {
	      sprintf(logbuf, "Error: --model 'fisher' and 'trend-only' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_FISHER;
	  } else if (!strcmp(argv[cur_arg + jj], "perm")) {
	    model_modifier |= MODEL_PERM;
	  } else if (!strcmp(argv[cur_arg + jj], "genedrop")) {
	    model_modifier |= MODEL_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + jj], "perm-count")) {
	    model_modifier |= MODEL_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + jj], "dom")) {
	    if (model_modifier & (MODEL_PREC | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PDOM;
	  } else if (!strcmp(argv[cur_arg + jj], "rec")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PREC;
	  } else if (!strcmp(argv[cur_arg + jj], "gen")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (mtest_adjust) {
	      sprintf(logbuf, "Error: --model perm-gen cannot be used with --adjust.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_PGEN;
	  } else if (!strcmp(argv[cur_arg + jj], "trend")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PTREND;
	  } else if (!strcmp(argv[cur_arg + jj], "trend-only")) {
	    if (model_modifier & (MODEL_FISHER | MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PTREND | MODEL_TRENDONLY;
	  } else if ((!memcmp(argv[cur_arg + jj], "mperm=", 6)) && (argv[cur_arg + jj][6] != '\0')) {
	    if (model_modifier & MODEL_PERM) {
	      sprintf(logbuf, "Error: --model 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: Duplicate --model 'mperm' modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    kk = atoi(&(argv[cur_arg + jj][6]));
	    if (kk < 1) {
	      sprintf(logbuf, "Error: Invalid --model mperm parameter '%s'.%s", &(argv[cur_arg + jj][6]), errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_mperm_val = (uint32_t)kk;
	  } else {
	    sprintf(logbuf, "Error: Invalid --model parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_MODEL;
      } else if (!memcmp(argptr2, "odel-dom", 9)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  logprint("Error: --model-dom must be used with --model.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PREC | MODEL_PGEN | MODEL_PTREND)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-dom flag deprecated.  Use '--model dom'.\n");
	model_modifier |= MODEL_PDOM;
      } else if (!memcmp(argptr2, "odel-gen", 9)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  logprint("Error: --model-gen must be used with --model.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-gen flag deprecated.  Use '--model gen'.\n");
	model_modifier |= MODEL_PGEN;
      } else if (!memcmp(argptr2, "odel-rec", 9)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  logprint("Error: --model-rec must be used with --model.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PTREND)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-rec flag deprecated.  Use '--model rec'.\n");
	model_modifier |= MODEL_PREC;
      } else if (!memcmp(argptr2, "odel-trend", 11)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  logprint("Error: --model-trend must be used with --model.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PREC)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-trend flag deprecated.  Use '--model trend'.\n");
	model_modifier |= MODEL_PTREND;
      } else if (!memcmp(argptr2, "perm", 5)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mperm parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (model_modifier & (MODEL_PERM | MODEL_MPERM)) {
	  sprintf(logbuf, "Error: --mperm cannot be used with --%s %sperm.%s", (model_modifier & MODEL_ASSOC)? "assoc" : "model", (model_modifier & MODEL_PERM)? "" : "m", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --mperm flag deprecated.  Use e.g. '--model mperm=[value]'.");
	mperm_val = (uint32_t)ii;
	model_mperm_val = mperm_val;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'n':
      if (!memcmp(argptr2, "o-fid", 6)) {
	fam_col_1 = 0;
      } else if (!memcmp(argptr2, "o-parents", 10)) {
	fam_col_34 = 0;
      } else if (!memcmp(argptr2, "o-sex", 6)) {
	if (filter_sex) {
	  logprint("Error: --filter-males/--filter-females cannot be used with --no-sex.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	fam_col_5 = 0;
	allow_no_sex = 1;
      } else if (!memcmp(argptr2, "o-pheno", 8)) {
	fam_col_6 = 0;
      } else if (!memcmp(argptr2, "onfounders", 11)) {
	nonfounders = 1;
      } else {
        goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'o':
      if (!memcmp(argptr2, "utput-missing-genotype", 23)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	output_missing_geno = argv[cur_arg + 1][0];
	if ((argv[cur_arg + 1][1] != '\0') || (((unsigned char)output_missing_geno) <= ' ')) {
	  sprintf(logbuf, "Error: Invalid --output-missing-genotype parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "utput-missing-phenotype", 24)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > 31) {
	  logprint("Error: --output-missing-phenotype string too long (max 31 chars).\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (sscanf(argv[cur_arg + 1], "%lg", &dxx) != 1) {
	  logprint("Error: --output-missing-phenotype parameter currently must be numeric.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	memcpy(output_missing_pheno, argv[cur_arg + 1], jj + 1);
      } else if (memcmp(argptr2, "ut", 3)) {
	// --out is a special case due to logging
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'p':
      if (!memcmp(argptr2, "ed", 3)) {
	if ((load_params & 0x3fa) || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 2;
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --ped parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "heno", 5)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (makepheno_str) {
	  logprint("Error: --pheno and --make-pheno flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&phenoname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "heno-name", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (mpheno_col != 0) {
	  logprint("Error: --mpheno and --pheno-name flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (alloc_string(&phenoname_str, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "heno-merge", 11)) {
	pheno_modifier |= PHENO_MERGE;
      } else if (!memcmp(argptr2, "rune", 5)) {
	prune = 1;
      } else if (!memcmp(argptr2, "arallel", 8)) {
	if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	  sprintf(logbuf, "Error: --parallel cannot be used with '--distance square'.  Use '--distance\nsquare0' or plain --distance instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((dist_calc_type & DISTANCE_BIN) && (!(dist_calc_type & DISTANCE_SHAPEMASK))) {
	  sprintf(logbuf, "Error: --parallel cannot be used with plain '--distance bin'.  Use '--distance\nbin square0' or '--distance bin triangle' instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	  sprintf(logbuf, "Error: --parallel cannot be used with '--make-rel square'.  Use '--make-rel\nsquare0' or plain '--make-rel' instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((rel_calc_type & REL_CALC_BIN) && (!(rel_calc_type & REL_CALC_SHAPEMASK))) {
	  sprintf(logbuf, "Error: --parallel cannot be used with plain '--make-rel bin'.  Use '--make-rel\nbin square0' or '--make-rel bin triangle' instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) {
	  sprintf(logbuf, "Error: --parallel and --distance-matrix/--matrix cannot be used together.  Use\n--distance instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_GROUPDIST) {
	  sprintf(logbuf, "Error: --parallel and --groupdist cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 2, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	parallel_idx = atoi(argv[cur_arg + 1]);
	if ((parallel_idx < 1) || (parallel_idx > PARALLEL_MAX)) {
	  sprintf(logbuf, "Error: Invalid --parallel job index '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	parallel_idx -= 1; // internal 0..(n-1) indexing
	parallel_tot = atoi(argv[cur_arg + 2]);
	if ((parallel_tot < 2) || (parallel_tot > PARALLEL_MAX) || (parallel_tot < parallel_idx)) {
	  sprintf(logbuf, "Error: Invalid --parallel total job count '%s'.%s", argv[cur_arg + 2], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "pc-gap", 7)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + 1], "%lg", &dxx) != 1) {
	  sprintf(logbuf, "Error: Invalid --ppc-gap parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dxx *= 1000;
	if (dxx < 0) {
	  ppc_gap = 0;
	} else if (dxx > 2147483647) {
	  ppc_gap = 2147483647;
	} else {
	  ppc_gap = (int)dxx;
	}
      } else if (!memcmp(argptr2, "erm", 4)) {
	model_modifier |= MODEL_PERM;
	logprint("Note: --perm flag deprecated.  Use e.g. '--model perm'.\n");
      } else if (!memcmp(argptr2, "erm-count", 10)) {
	model_modifier |= MODEL_PERM_COUNT;
	logprint("Note: --perm-count flag deprecated.  Use e.g. '--model perm-count'.\n");
      } else if (!memcmp(argptr2, "2", 2)) {
	if ((!(calculation_type & CALC_MODEL)) || (!(model_modifier & MODEL_ASSOC))) {
	  logprint("Error: --p2 must be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & MODEL_QMASK) {
	  sprintf(logbuf, "Error: --assoc 'qt-means' modifier does not make sense with --p2.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --p2 flag deprecated.  Use '--assoc p2 ...'.\n");
	model_modifier |= MODEL_ASSOC_P2;
      } else if (!memcmp(argptr2, "filter", 7)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (sscanf(argv[cur_arg + 1], "%lg", &dxx) != 1) {
	  sprintf(logbuf, "Error: Invalid --pfilter parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx <= 0.0) || (dxx >= 1.0)) {
	  sprintf(logbuf, "Error: --pfilter threshold must be between 0 and 1 exclusive.%s", errstr_append);
	}
	pfilter = dxx;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'q':
      if (!memcmp(argptr2, "t-means", 8)) {
	if ((!(calculation_type & CALC_MODEL)) || (!(model_modifier & MODEL_ASSOC))) {
	  sprintf(logbuf, "Error: --qt-means must be used with --assoc.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --qt-means flag deprecated.  Use '--assoc qt-means ...'.\n");
	model_modifier |= MODEL_QT_MEANS;
      } else if (!memcmp(argptr2, "q-plot", 7)) {
        if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --qq-plot must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --qq-plot flag deprecated.  Use '--adjust qq-plot'.\n");
	mtest_adjust |= ADJUST_QQ;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'r':
      if (!memcmp(argptr2, "emove", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&removename, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "el-cutoff", 10)) {
	if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel cannot be used with %s.  (Use a combination of\n--make-rel, --keep/--remove, and a filtering script.)%s", argptr, errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  if (sscanf(argv[cur_arg + 1], "%lg", &rel_cutoff) != 1) {
	    sprintf(logbuf, "Error: Invalid %s parameter '%s'.%s", argptr, argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((rel_cutoff <= 0.0) || (rel_cutoff >= 1.0)) {
	    sprintf(logbuf, "Error: Invalid %s parameter '%s'.%s", argptr, argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_REL_CUTOFF;
      } else if (!memcmp(argptr2, "egress-distance", 16)) {
	if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel and --regress-distance cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  regress_iters = strtoul(argv[cur_arg + 1], NULL, 10);
	  if ((regress_iters < 2) || (regress_iters == ULONG_MAX)) {
	    sprintf(logbuf, "Error: Invalid --regress-distance jackknife iteration count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (ii == 2) {
	    regress_d = atoi(argv[cur_arg + 2]);
	    if (regress_d <= 0) {
	      sprintf(logbuf, "Error: Invalid --regress-distance jackknife delete parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	calculation_type |= CALC_REGRESS_DISTANCE;
      } else if (!memcmp(argptr2, "egress-rel", 11)) {
	if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel and --regress-rel flags cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (rel_calc_type & REL_CALC_SINGLE_PREC) {
	  sprintf(logbuf, "Error: --regress-rel cannot currently be used with a single-precision\nrelationship matrix.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  regress_rel_iters = strtoul(argv[cur_arg + 1], NULL, 10);
	  if ((regress_rel_iters < 2) || (regress_rel_iters == ULONG_MAX)) {
	    sprintf(logbuf, "Error: Invalid --regress-rel jackknife iteration count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (ii == 2) {
	    regress_rel_d = atoi(argv[cur_arg + 2]);
	    if (regress_rel_d <= 0) {
	      sprintf(logbuf, "Error: Invalid --regress-rel jackknife delete parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	calculation_type |= CALC_REGRESS_REL;
      } else if (!memcmp(argptr2, "egress-pcs", 11)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 5, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&evecname, argv[cur_arg + 1], argptr, 9);
	if (retval) {
	  goto main_ret_1;
	}
	for (jj = 2; jj <= ii; jj++) {
	  if (!strcmp(argv[cur_arg + jj], "normalize-pheno") && (!regress_pcs_normalize_pheno)) {
	    regress_pcs_normalize_pheno = 1;
	  } else if (!strcmp(argv[cur_arg + jj], "sex-specific") && (!regress_pcs_sex_specific)) {
	    regress_pcs_sex_specific = 1;
	  } else if (!strcmp(argv[cur_arg + jj], "clip") && (!regress_pcs_clip)) {
	    regress_pcs_clip = 1;
	  } else if ((max_pcs != MAX_PCS_DEFAULT) || (argv[cur_arg + jj][0] < '0') || (argv[cur_arg + jj][0] > '9')) {
	    sprintf(logbuf, "Error: Invalid --regress-pcs parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else {
	    max_pcs = atoi(argv[cur_arg + jj]);
	    if (max_pcs < 1) {
	      sprintf(logbuf, "Error: Invalid --regress-pcs maximum principal component count '%s'.%s", argv[cur_arg + jj], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	calculation_type |= CALC_REGRESS_PCS;
      } else if (!memcmp(argptr2, "egress-pcs-distance", 20)) {
	if (calculation_type & CALC_REGRESS_PCS) {
	  sprintf(logbuf, "Error: --regress-pcs-distance cannot be used with --regress-pcs.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_DISTANCE) {
	  sprintf(logbuf, "Error: --regress-pcs-distance cannot be used with --distance.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(argc, argv, cur_arg, 1, 11, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&evecname, argv[cur_arg + 1], argptr, 9);
	if (retval) {
	  goto main_ret_1;
	}
	for (jj = 2; jj <= ii; jj++) {
	  if (!memcmp(argv[cur_arg + jj], "normalize-pheno", 16) && (!regress_pcs_normalize_pheno)) {
	    regress_pcs_normalize_pheno = 1;
	  } else if (!memcmp(argv[cur_arg + jj], "sex-specific", 13) && (!regress_pcs_sex_specific)) {
	    regress_pcs_sex_specific = 1;
	  } else if (!memcmp(argv[cur_arg + jj], "square", 7)) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (parallel_tot > 1) {
	      sprintf(logbuf, "Error: --parallel cannot be used with '--regress-pcs-distance square'.  Use\nthe square0 or triangle shape instead.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ;
	  } else if (!memcmp(argv[cur_arg + jj], "square0", 8)) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square0' and 'triangle' modifiers can't coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ0;
	  } else if (!memcmp(argv[cur_arg + jj], "triangle", 9)) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square0' and 'triangle' modifiers can't coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_TRI;
	  } else if (!memcmp(argv[cur_arg + jj], "gz", 3)) {
	    if (dist_calc_type & DISTANCE_BIN) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_GZ;
	  } else if (!memcmp(argv[cur_arg + jj], "bin", 4)) {
	    if (dist_calc_type & DISTANCE_GZ) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_BIN;
	  } else if (!memcmp(argv[cur_arg + jj], "ibs", 4)) {
	    if (dist_calc_type & DISTANCE_IBS) {
	      logprint("Error: Duplicate --regress-pcs-distance 'ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_IBS;
	  } else if (!memcmp(argv[cur_arg + jj], "1-ibs", 6)) {
	    if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	      logprint("Error: Duplicate --regress-pcs-distance '1-ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_1_MINUS_IBS;
	  } else if (!memcmp(argv[cur_arg + jj], "alct", 5)) {
	    if (dist_calc_type & DISTANCE_ALCT) {
	      logprint("Error: Duplicate --regress-pcs-distance 'alct' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_ALCT;
	  } else if (!memcmp(argv[cur_arg + jj], "3d", 3)) {
	    distance_3d = 1;
	  } else if (!memcmp(argv[cur_arg + jj], "flat-missing", 13)) {
	    distance_flat_missing = 1;
	  } else if ((max_pcs != MAX_PCS_DEFAULT) || (argv[cur_arg + jj][0] < '0') || (argv[cur_arg + jj][0] > '9')) {
	    sprintf(logbuf, "Error: Invalid --regress-pcs-distance parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else {
	    max_pcs = atoi(argv[cur_arg + jj]);
	    if (max_pcs < 1) {
	      sprintf(logbuf, "Error: Invalid --regress-pcs-distance maximum PC count '%s'.%s", argv[cur_arg + jj], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	calculation_type |= CALC_REGRESS_PCS_DISTANCE;
      } else if (!memcmp(argptr2, "ead-freq", 9)) {
	if (calculation_type & CALC_FREQ) {
	  sprintf(logbuf, "Error: --freq and --read-freq flags cannot coexist.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&freqname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if ((!memcmp(argptr2, "ecode", 6)) || (!memcmp(argptr2, "ecode 12", 9)) || (!memcmp(argptr2, "ecode lgen", 11)) || (!memcmp(argptr2, "ecode AD", 9)) || (!memcmp(argptr2, "ecode A", 8)) || (!memcmp(argptr2, "ecode list", 11)) || (!memcmp(argptr2, "ecode rlist", 12))) {
	if (argptr2[5] == ' ') {
	  if (argptr2[6] == '1') {
	    recode_modifier |= RECODE_12;
	  } else if (argptr2[6] == 'l') {
	    if (argptr2[7] == 'g') {
	      recode_modifier |= RECODE_LGEN;
	    } else {
	      recode_modifier |= RECODE_LIST;
	    }
	  } else if (argptr2[6] == 'r') {
	    recode_modifier |= RECODE_RLIST;
	  } else if (argptr2[7] == 'D') {
	    recode_modifier |= RECODE_AD;
	  } else {
	    recode_modifier |= RECODE_A;
	  }
	  kk = 1;
	} else {
	  kk = 0;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 3 - kk, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (jj = 1; jj <= ii; jj++) {
	  if (!memcmp(argv[cur_arg + jj], "12", 3)) {
	    if (recode_modifier & (RECODE_A | RECODE_AD)) {
	      sprintf(logbuf, "Error: --recode '12' modifier cannot be used with 'A' or 'AD'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_12;
	  } else if (!memcmp(argv[cur_arg + jj], "compound-genotypes", 19)) {
	    if (recode_type_set(&recode_modifier, RECODE_COMPOUND)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
          } else if ((!argv[cur_arg + jj][1]) && (tolower(argv[cur_arg + jj][0]) == 'a')) {
	    if (recode_type_set(&recode_modifier, RECODE_A)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if ((!argv[cur_arg + jj][2]) && (tolower(argv[cur_arg + jj][0]) == 'a') && (tolower(argv[cur_arg + jj][1]) == 'd')) {
	    if (recode_type_set(&recode_modifier, RECODE_AD)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!memcmp(argv[cur_arg + jj], "tab", 4)) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_TAB;
	  } else if (!memcmp(argv[cur_arg + jj], "tabx", 5)) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_TAB | RECODE_DELIMX;
	  } else if (!memcmp(argv[cur_arg + jj], "spacex", 7)) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_DELIMX;
	  } else if (!memcmp(argv[cur_arg + jj], "lgen", 5)) {
	    if (recode_type_set(&recode_modifier, RECODE_LGEN)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!memcmp(argv[cur_arg + jj], "lgen-ref", 9)) {
	    if (recode_type_set(&recode_modifier, RECODE_LGEN_REF)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!memcmp(argv[cur_arg + jj], "list", 5)) {
	    if (recode_type_set(&recode_modifier, RECODE_LIST)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!memcmp(argv[cur_arg + jj], "rlist", 6)) {
	    if (recode_type_set(&recode_modifier, RECODE_RLIST)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!memcmp(argv[cur_arg + jj], "transpose", 10)) {
	    if (recode_type_set(&recode_modifier, RECODE_TRANSPOSE)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else {
	    sprintf(logbuf, "Error: Invalid --recode parameter '%s'.%s%s", argv[cur_arg + jj], ((jj == ii) && (!outname_end))? "  (Did you forget '--out'?)" : "", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_RECODE;
      } else if (!memcmp(argptr2, "ecode-allele", 13)) {
	if (!(recode_modifier & (RECODE_A | RECODE_AD))) {
	  sprintf(logbuf, "Error: --recode-allele must be used with --recode A or --recode AD.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (recode_modifier & RECODE_12) {
	  sprintf(logbuf, "Error: --recode-allele cannot be used with --recode 12.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&recode_allele_name, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eference-allele", 16)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&refalleles, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eference", 9)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&lgen_reference_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	lgen_modifier |= LGEN_REFERENCE;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 's':
      if (!memcmp(argptr2, "ample", 6)) {
	if ((load_params & 0x27f) || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 0x200;
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --sample parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(samplename, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "eed", 4)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	rseed = (unsigned long int)atol(argv[cur_arg + 1]);
	if (rseed == 0) {
	  sprintf(logbuf, "Error: Invalid --seed parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "np", 3)) {
        if (markername_from) {
	  sprintf(logbuf, "Error: --snp cannot be used with --from.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (marker_pos_start != -1) {
	  sprintf(logbuf, "Error: --snp cannot be used with --from-bp/-kb/-mb.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (chrom_info.chrom_mask) {
	  sprintf(logbuf, "Error: --snp cannot be used with --autosome(-xy) or --chr%s.%s", chrom_exclude? "-excl" : "", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_snp, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "nps", 4)) {
	if (markername_from) {
	  sprintf(logbuf, "Error: --snps cannot be used with --from.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (marker_pos_start != -1) {
	  sprintf(logbuf, "Error: --snps cannot be used with --from-bp/-kb/-mb.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_snp) {
	  sprintf(logbuf, "Error: --snps cannot be used with --snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	// mise well allow --snps + --autosome/--autosome-xy/--chr/--chr-excl
	retval = parse_marker_ranges(param_count(argc, argv, cur_arg), range_delim, &(argv[cur_arg]), &snps_flag_markers, &snps_flag_starts_range, &snps_flag_ct, &snps_flag_max_len);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "et-hh-missing", 14)) {
	set_hh_missing = 1;
      } else if (!memcmp(argptr2, "et", 3)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&set_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ubset", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!set_fname) {
	  logprint("Error: --subset must be used with --set.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (memcmp(argptr2, "ilent", 6)) {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 't':
      if (!memcmp(argptr2, "ail-pheno", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (makepheno_str) {
	  sprintf(logbuf, "Error: --tail-pheno cannot be used with --make-pheno.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + 1], "%lg", &tail_bottom) != 1) {
	  sprintf(logbuf, "Error: Invalid --tail-pheno parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii == 1) {
	  tail_top = tail_bottom;
	} else {
	  if (sscanf(argv[cur_arg + 2], "%lg", &tail_top) != 1) {
	    sprintf(logbuf, "Error: Invalid --tail-pheno parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (tail_bottom > tail_top) {
	  sprintf(logbuf, "Error: Ltop cannot be larger than Hbottom for --tail-pheno.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	tail_pheno = 1;
      } else if (!memcmp(argptr2, "hreads", 7)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --threads parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (ii > MAX_THREADS) {
	  sprintf(logbuf, "Note: Reducing --threads parameter to %d.  (If this is not large enough,\nrecompile with a larger MAX_THREADS setting.)\n", MAX_THREADS);
	  logprintb();
	  ii = MAX_THREADS;
	}
	g_thread_ct = ii;
      } else if (!memcmp(argptr2, "ab", 3)) {
	logprint("Note: --tab flag deprecated.  Use '--recode tab ...'.\n");
	if (recode_modifier & RECODE_DELIMX) {
	  sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	recode_modifier |= RECODE_TAB;
      } else if (!memcmp(argptr2, "ranspose", 9)) {
	logprint("Note: --transpose flag deprecated.  Use '--recode transpose ...'.\n");
	if (recode_modifier & RECODE_LGEN) {
	  sprintf(logbuf, "Error: --recode 'transpose' and 'lgen' modifiers cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	recode_modifier |= RECODE_TRANSPOSE;
      } else if (!memcmp(argptr2, "fam", 4)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > FNAMESIZE - 1) {
	  logprint("Error: --tfam filename prefix too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(famname, argv[cur_arg + 1], jj + 1);
	load_rare |= LOAD_RARE_TFAM;
      } else if (!memcmp(argptr2, "file", 5)) {
	if (load_params || (load_rare & (~LOAD_RARE_TFAM))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  jj = strlen(argv[cur_arg + 1]);
	  if (jj > FNAMESIZE - 6) {
	    logprint("Error: --tfile filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(pedname, argv[cur_arg + 1], jj);
	  memcpy(&(pedname[jj]), ".tped", 6);
	  if (!(load_rare & LOAD_RARE_TFAM)) {
	    memcpy(famname, argv[cur_arg + 1], jj);
	    memcpy(&(famname[jj]), ".tfam", 6);
	  }
	} else {
	  memcpy(pedname, "wdist.tped", 11);
	  if (!(load_rare & LOAD_RARE_TFAM)) {
	    memcpy(famname, "wdist.tfam", 11);
	  }
	}
	load_rare |= LOAD_RARE_TRANSPOSE;
      } else if (!memcmp(argptr2, "ped", 4)) {
	if (load_params || (load_rare & (~(LOAD_RARE_TRANSPOSE | LOAD_RARE_TFAM)))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > FNAMESIZE - 1) {
	  logprint("Error: --tped filename prefix too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(pedname, argv[cur_arg + 1], jj + 1);
	load_rare |= LOAD_RARE_TPED;
      } else if (!memcmp(argptr2, "o", 2)) {
	if (chrom_info.chrom_mask) {
	  sprintf(logbuf, "Error: --to cannot be used with --autosome(-xy) or --chr%s.%s", chrom_exclude? "-excl" : "", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_snp) {
	  sprintf(logbuf, "Error: --to cannot be used with --snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_flag_markers) {
	  sprintf(logbuf, "Error: --to cannot be used with --snps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_to, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if ((!memcmp(argptr2, "o-bp", 5)) || (!memcmp(argptr2, "o-kb", 5)) || (!memcmp(argptr2, "o-mb", 5))) {
	if (markername_snp) {
	  sprintf(logbuf, "Error: --to-bp/-kb/-mb cannot be used with --snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_flag_markers) {
	  sprintf(logbuf, "Error: --to-bp/-kb/-mb cannot be used with --snps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_to) {
	  sprintf(logbuf, "Error: --to-bp/-kb/-mb cannot be used with --to.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((!markername_from) && (marker_pos_start == -1)) {
	  marker_pos_start = 0;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argptr2[2];
	if (cc == 'b') {
	  if (atoiz(argv[cur_arg + 1], &ii)) {
	    sprintf(logbuf, "Error: Invalid --to-bp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  if (marker_pos_end != -1) {
	    logprint("Error: Multiple --to-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (sscanf(argv[cur_arg + 1], "%lg", &dxx) != 1) {
	    sprintf(logbuf, "Error: Invalid --to-kb/-mb parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  dxx *= (cc == 'k')? 1000 : 1000000;
	  if (dxx < 0) {
	    ii = 0;
	  } else if (dxx > 2147483647) {
	    ii = 2147483647;
	  } else {
	    ii = (int)dxx;
	  }
	}
	if (ii < marker_pos_start) {
	  marker_pos_end = marker_pos_start;
	  marker_pos_start = ii;
	} else {
	  marker_pos_end = ii;
	}
      } else if (!memcmp(argptr2, "rend", 5)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  sprintf(logbuf, "Error: --trend must be used with --model.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (model_modifier & MODEL_FISHER) {
	  sprintf(logbuf, "Error: --trend cannot be used with --model fisher.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	  sprintf(logbuf, "Error: --trend cannot be used with --model dom/rec/gen.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --trend flag deprecated.  Use '--model trend-only ...'.\n");
	model_modifier |= MODEL_PTREND | MODEL_TRENDONLY;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'u':
      if (!memcmp(argptr2, "nrelated-heritability", 22)) {
#ifdef NOLAPACK
        sprintf(logbuf, "Error: --unrelated-heritability requires WDIST to be built with LAPACK.\n");
	goto main_ret_INVALID_CMDLINE_3;
#else
	if (rel_calc_type & REL_CALC_COV) {
	  sprintf(logbuf, "Error: --unrelated-heritability flag cannot coexist with a covariance\nmatrix calculation.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel and --unrelated-heritability cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (rel_calc_type & REL_CALC_SINGLE_PREC) {
	  sprintf(logbuf, "Error: --unrelated-heritability flag cannot be used with a single-precision\nrelationship matrix calculation.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 0, 4, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii) {
	  if (!strcmp(argv[cur_arg + 1], "strict")) {
	    calculation_type |= CALC_UNRELATED_HERITABILITY_STRICT;
	    jj = 2;
	  } else {
	    jj = 1;
	  }
	  if (ii >= jj) {
	    if (sscanf(argv[cur_arg + jj], "%lg", &unrelated_herit_tol) != 1) {
	      sprintf(logbuf, "Error: Invalid --unrelated-heritability EM tolerance parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (unrelated_herit_tol <= 0.0) {
	      sprintf(logbuf, "Error: Invalid --unrelated-heritability EM tolerance parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ii > jj) {
	      if (sscanf(argv[cur_arg + jj + 1], "%lg", &unrelated_herit_covg) != 1) {
		sprintf(logbuf, "Error: Invalid --unrelated-heritability genomic covariance prior '%s'.%s", argv[cur_arg + jj + 1], errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      if ((unrelated_herit_covg <= 0.0) || (unrelated_herit_covg > 1.0)) {
		sprintf(logbuf, "Error: Invalid --unrelated-heritability genomic covariance prior '%s'.%s", argv[cur_arg + jj + 1], errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      if (ii == jj + 2) {
		if (sscanf(argv[cur_arg + jj + 2], "%lg", &unrelated_herit_covr) != 1) {
		  sprintf(logbuf, "Error: Invalid --unrelated-heritability residual covariance prior '%s'.%s", argv[cur_arg + jj + 2], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
		if ((unrelated_herit_covr <= 0.0) || (unrelated_herit_covr > 1.0)) {
		  sprintf(logbuf, "Error: Invalid --unrelated-heritability residual covariance prior '%s'.%s", argv[cur_arg + jj + 2], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
	      } else {
		unrelated_herit_covr = 1.0 - unrelated_herit_covg;
	      }
	    }
	  }
	}
	calculation_type |= CALC_UNRELATED_HERITABILITY;
#endif
      } else if (!memcmp(argptr2, "pdate-alleles", 14)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Error: --update-alleles is not implemented yet.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
      } else if (!memcmp(argptr2, "pdate-chr", 10)) {
	logprint("Note: --update-chr flag deprecated.  Use '--update-map chr'.\n");
	update_map_modifier = UPDATE_MAP_CHR;
      } else if (!memcmp(argptr, "pdate-cm", 9)) {
	if (update_map_modifier) {
	  logprint("Error: --update-map 'cm' modifier cannot be used with 'chr'.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --update-cm flag deprecated.  Use '--update-map cm'.\n");
	update_map_modifier = UPDATE_MAP_CM;
      } else if (!memcmp(argptr2, "pdate-ids", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Error: --update-ids is not implemented yet.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
      } else if (!memcmp(argptr2, "pdate-map", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 2, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Error: --update-map is not implemented yet.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto main_ret_1;
      } else if (!memcmp(argptr2, "pdate-name", 11)) {
	if (update_map_fname) {
	  sprintf(logbuf, "Error: --update-name cannot be used without --update-map.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_map_modifier & (UPDATE_MAP_CHR | UPDATE_MAP_CM)) {
	  sprintf(logbuf, "Error: --update-map 'name' modifier cannot be used with '%s'.%s", (update_map_modifier == UPDATE_MAP_CHR)? "chr" : "cm", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --update-name flag deprecated.  Use '--update-map name'.\n");
	update_map_modifier = UPDATE_MAP_NAME;
      } else if (!memcmp(argptr2, "pdate-parents", 14)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_ids_fname) {
	  logprint("Error: --update-parents cannot be used with --update-ids.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Error: --update-parents is not implemented yet.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
      } else if (!memcmp(argptr2, "pdate-sex", 10)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_ids_fname) {
	  logprint("Error: --update-sex cannot be used with --update-ids.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Error: --update-sex is not implemented yet.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'w':
      if (!memcmp(argptr2, "rite-snplist", 13)) {
	calculation_type |= CALC_WRITE_SNPLIST;
      } else if (!memcmp(argptr2, "indow", 6)) {
        if (!markername_snp) {
	  sprintf(logbuf, "Error: --window must be used with --snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (sscanf(argv[cur_arg + 1], "%lg", &dxx) != 1) {
	  sprintf(logbuf, "Error: Invalid --window parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        dxx *= 500;
	if (dxx < 1) {
	  sprintf(logbuf, "Error: --window parameter '%s' too small.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (dxx > 2147483647) {
	  snp_window_size = 2147483647;
	} else {
	  snp_window_size = (int)dxx;
	}
      } else if (!memcmp(argptr2, "ithin", 6)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cluster_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ith-phenotype", 14)) {
	if (!covar_fname) {
	  sprintf(logbuf, "Error: --with-phenotype cannot be used without --covar.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	write_covar_modifier |= WRITE_COVAR_PHENO;
      } else if (!memcmp(argptr2, "ith-reference", 14)) {
	if ((recode_modifier & RECODE_TYPEMASK) != RECODE_LGEN) {
	  logprint("Error: --with-reference must be used with --recode lgen.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --with-reference flag deprecated.  Use '--recode lgen-ref' instead.\n");
	recode_modifier += RECODE_LGEN_REF - RECODE_LGEN;
      } else if (!memcmp(argptr2, "rite-covar", 11)) {
	if (calculation_type & (CALC_MAKE_BED | CALC_RECODE)) {
	  sprintf(logbuf, "Error: --write-covar cannot be used with --make-bed or --recode.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!covar_fname) {
	  sprintf(logbuf, "Error: --write-covar cannot be used without --covar.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_WRITE_COVAR;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    default:
      goto main_ret_INVALID_CMDLINE_2;
    }

  } while ((++cur_flag) < flag_ct);
  if (!outname_end) {
    outname_end = &(outname[5]);
  }
  if (load_rare) {
    if (load_rare == LOAD_RARE_GRM) {
      if ((!(calculation_type & CALC_REL_CUTOFF)) || (calculation_type & (~(CALC_REL_CUTOFF | CALC_RELATIONSHIP)))) {
	sprintf(logbuf, "Error: --grm currently must be used with --rel-cutoff (possibly combined with\n--make-grm).%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
  }
  if (prune && (!phenoname) && (!fam_col_6)) {
    sprintf(logbuf, "Error: --prune and --no-pheno cannot coexist without an alternate phenotype\nfile.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  } else if ((calculation_type & CALC_LOAD_DISTANCES) && (!(calculation_type & (CALC_GROUPDIST | CALC_REGRESS_DISTANCE)))) {
    sprintf(logbuf, "Error: --load-dists cannot be used without either --groupdist or\n--regress-distance.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (update_map_modifier && (!update_map_fname)) {
    sprintf(logbuf, "Error: --update-%s cannot be used without --update-map.%s", (update_map_modifier == UPDATE_MAP_CHR)? "chr" : "cm", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }

  // --from-bp/-kb/-mb without any --to/--to-bp/...: include to end of
  // chromosome
  if ((marker_pos_start != -1) && (!markername_to) && (marker_pos_end == -1)) {
    marker_pos_end = 2147483647;
  }
  if (!chrom_info.chrom_mask) {
    chrom_info.chrom_mask = species_valid_chrom_mask[chrom_info.species];
  }
  if (((marker_pos_start != -1) && (!markername_to)) || ((marker_pos_end != -1) && (!markername_from))) {
    // require exactly one chromosome to be defined given --from-bp/--to-bp
    // without --from/--to
    if (chrom_info.chrom_mask & (chrom_info.chrom_mask - 1LLU)) {
      sprintf(logbuf, "Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb require a single chromosome to be\nidentified (either explicitly with --chr, or implicitly with --from/--to).%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
  }

  if (calculation_type & CALC_MODEL) {
    if (!(model_modifier & (MODEL_ASSOC | MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
      if (mtest_adjust) {
	sprintf(logbuf, "Error: In order to use --model with --adjust, you must include the 'trend',\n'trend-only', 'dom', or 'rec' modifier.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
    if (model_cell_ct == -1) {
      model_cell_ct = (model_modifier & MODEL_FISHER)? 0 : 5;
    }
    if ((model_modifier & (MODEL_PERM | MODEL_MPERM | MODEL_GENEDROP)) == MODEL_GENEDROP) {
      model_modifier |= MODEL_PERM;
    }
  }

  if ((!calculation_type) && (load_rare != LOAD_RARE_LGEN) && (load_rare != LOAD_RARE_DUMMY) && (!(load_rare & LOAD_RARE_TRANSPOSE_MASK)) && (famname[0] || load_rare)) {
    goto main_ret_NULL_CALC;
  }
  if (!(load_params || load_rare)) {
    sprintf(logbuf, "Error: No input dataset.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  free_cond(subst_argv);
  free_cond(script_buf);
  free_cond(rerun_buf);
  free_cond(flag_buf);
  free_cond(flag_map);
  subst_argv = NULL;
  script_buf = NULL;
  rerun_buf = NULL;
  flag_buf = NULL;
  flag_map = NULL;

  bubble = (char*)malloc(67108864 * sizeof(char));
  if (!bubble) {
    goto main_ret_NOMEM;
  }

  // see e.g. http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system .
#ifdef __APPLE__
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  llxx = 0;
  sztmp = sizeof(int64_t);
  sysctl(mib, 2, &llxx, &sztmp, NULL, 0);
  llxx /= 1048576;
#else
#if _WIN32
  memstatus.dwLength = sizeof(memstatus);
  GlobalMemoryStatusEx(&memstatus);
  llxx = memstatus.ullTotalPhys / 1048576;
#else
  llxx = ((size_t)sysconf(_SC_PHYS_PAGES)) * ((size_t)sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
  if (!llxx) {
    default_alloc_mb = WKSPACE_DEFAULT_MB;
  } else if (llxx < (WKSPACE_MIN_MB * 2)) {
    default_alloc_mb = WKSPACE_MIN_MB;
  } else {
    default_alloc_mb = llxx / 2;
  }
  if (!malloc_size_mb) {
    malloc_size_mb = default_alloc_mb;
  } else if (malloc_size_mb < WKSPACE_MIN_MB) {
    malloc_size_mb = WKSPACE_MIN_MB;
  }
  if (llxx) {
    sprintf(logbuf, "%" PRId64 " MB RAM detected; reserving %" PRIdPTR " MB for main workspace.\n", llxx, malloc_size_mb);
  } else {
    sprintf(logbuf, "Failed to calculate system memory.  Attempting to reserve %" PRIdPTR " MB.\n", malloc_size_mb);
  }
  logprintb();
  wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  while (!wkspace_ua) {
    malloc_size_mb = (malloc_size_mb * 3) / 4;
    if (malloc_size_mb < WKSPACE_MIN_MB) {
      malloc_size_mb = WKSPACE_MIN_MB;
    }
    wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
    if (wkspace_ua) {
      sprintf(logbuf, "Allocated %" PRIdPTR " MB successfully, after larger attempt(s) failed.\n", malloc_size_mb);
      logprintb();
    } else if (malloc_size_mb == WKSPACE_MIN_MB) {
      goto main_ret_NOMEM;
    }
  }
  // force 64-byte align on OS X to make cache line sensitivity work
  wkspace = (unsigned char*)CACHEALIGN((uintptr_t)wkspace_ua);
  wkspace_base = wkspace;
  wkspace_left = malloc_size_mb * 1048576 - (uintptr_t)(wkspace - wkspace_ua);
  free(bubble);

  if (!rseed) {
    rseed = (unsigned long int)time(NULL);
  }
  sfmt_init_gen_rand(&sfmt, rseed);

  tbuf[MAXLINELEN - 6] = ' ';
  tbuf[MAXLINELEN - 1] = ' ';
  if (load_rare == LOAD_RARE_GRM) {
    // --rel-cutoff batch mode special case
    retval = rel_cutoff_batch(pedname, outname, outname_end, rel_cutoff, rel_calc_type);
  } else if (genname[0]) {
    if (calculation_type & (~(CALC_DISTANCE | CALC_REGRESS_DISTANCE))) {
      logprint("Error: Only --distance calculations are currently supported with --data.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
    } else {
      if (!missing_code) {
	missing_code = (char*)"NA";
      }
      retval = wdist_dosage(calculation_type, dist_calc_type, genname, samplename, outname, outname_end, missing_code, distance_3d, distance_flat_missing, exponent, maf_succ, regress_iters, regress_d, g_thread_ct, parallel_idx, parallel_tot);
    }
  } else {
  // famname[0] indicates binary vs. text
  // extractname, excludename, keepname, and removename indicate the presence
  // of their respective flags
  // filtername indicates existence of filter
  // freqname signals --read-freq
    // if (load_rare) {
    if (load_rare || (!famname[0])) {
      sptr = outname_end;
      if (bed_suffix_conflict(calculation_type, recode_modifier)) {
        memcpy(sptr, "-working", 9);
	sptr = &(sptr[8]);
      }
      uii = (sptr - outname);
      if (load_rare == LOAD_RARE_LGEN) {
        retval = lgen_to_bed(pedname, outname, sptr, missing_pheno, affection_01, lgen_modifier, lgen_reference_fname, &chrom_info);
      } else if (load_rare & LOAD_RARE_TRANSPOSE_MASK) {
        retval = transposed_to_bed(pedname, famname, outname, sptr, missing_geno, &chrom_info);
      } else if (load_rare & LOAD_RARE_DUMMY) {
	retval = generate_dummy(outname, sptr, dummy_flags, (uint32_t)dummy_marker_ct, (uint32_t)dummy_indiv_ct, dummy_missing_geno, dummy_missing_pheno);
      } else {
        retval = ped_to_bed(pedname, mapname, outname, sptr, fam_col_1, fam_col_34, fam_col_5, fam_col_6, affection_01, missing_pheno, &chrom_info);
	fam_col_1 = 1;
	fam_col_34 = 1;
	fam_col_5 = 1;
	if (!fam_col_6) {
          fam_col_6 = 1;
	  missing_pheno = -9;
	}
      }
      if (retval || (!calculation_type)) {
	goto main_ret_2;
      }
      memcpy(pedname, outname, uii);
      memcpy(&(pedname[uii]), ".bed", 5);
      memcpy(mapname, outname, uii);
      memcpy(&(mapname[uii]), ".bim", 5);
      memcpy(famname, outname, uii);
      memcpy(&(famname[uii]), ".fam", 5);
      *outname_end = '\0';
    }
    retval = wdist(outname, outname_end, pedname, mapname, famname, phenoname, extractname, excludename, keepname, removename, filtername, freqname, loaddistname, evecname, mergename1, mergename2, mergename3, makepheno_str, phenoname_str, refalleles, recode_allele_name, covar_fname, cluster_fname, set_fname, subset_fname, update_alleles_fname, update_map_fname, update_ids_fname, update_parents_fname, update_sex_fname, loop_assoc_fname, filterval, mfilter_col, filter_case_control, filter_sex, filter_founder_nonf, fam_col_1, fam_col_34, fam_col_5, fam_col_6, missing_geno, missing_pheno, output_missing_geno, output_missing_pheno, mpheno_col, pheno_modifier, prune, affection_01, &chrom_info, exponent, min_maf, max_maf, geno_thresh, mind_thresh, hwe_thresh, hwe_all, rel_cutoff, tail_pheno, tail_bottom, tail_top, calculation_type, rel_calc_type, dist_calc_type, groupdist_iters, groupdist_d, regress_iters, regress_d, regress_rel_iters, regress_rel_d, unrelated_herit_tol, unrelated_herit_covg, unrelated_herit_covr, ibc_type, parallel_idx, (uint32_t)parallel_tot, ppc_gap, allow_no_sex, must_have_sex, nonfounders, genome_output_gz, genome_output_full, genome_ibd_unbounded, ld_window_size, ld_window_kb, ld_window_incr, ld_last_param, maf_succ, regress_pcs_normalize_pheno, regress_pcs_sex_specific, regress_pcs_clip, max_pcs, freq_counts, freqx, distance_flat_missing, recode_modifier, allelexxxx, merge_type, indiv_sort, keep_allele_order, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, snps_flag_markers, snps_flag_starts_range, snps_flag_ct, snps_flag_max_len, set_hh_missing, covar_modifier, covar_str, mcovar_col, update_map_modifier, write_covar_modifier, model_modifier, (uint32_t)model_cell_ct, (uint32_t)model_mperm_val, ci_size, pfilter, mtest_adjust, adjust_lambda, gxe_mcovar, aperm_min, aperm_max, aperm_alpha, aperm_beta, aperm_init_interval, aperm_interval_slope);
  }
 main_ret_2:
  free(wkspace_ua);
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
  main_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  main_ret_INVALID_CMDLINE_2:
    invalid_arg(argv[cur_arg]);
    logprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_INVALID_CMDLINE_4:
    sprintf(logbuf, "Error: --%s conflicts with another input flag.%s", argptr, errstr_append);
  main_ret_INVALID_CMDLINE_3:
    logprintb();
  main_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_NULL_CALC:
    logprint(notestr_null_calc);
    fputs(cmdline_format_str, stdout);
    fputs(notestr_null_calc2, stdout);
    retval = RET_NULL_CALC;
  }
 main_ret_1:
  fclose_cond(scriptfile);
  dispmsg(retval);
  free_cond(subst_argv);
  free_cond(script_buf);
  free_cond(rerun_buf);
  free_cond(flag_buf);
  free_cond(flag_map);
  free_cond(makepheno_str);
  free_cond(phenoname_str);
  free_cond(refalleles);
  free_cond(filterval);
  free_cond(evecname);
  free_cond(filtername);
  free_cond(loaddistname);
  free_cond(freqname);
  free_cond(extractname);
  free_cond(excludename);
  free_cond(keepname);
  free_cond(removename);
  free_cond(phenoname);
  free_cond(recode_allele_name);
  free_cond(markername_from);
  free_cond(markername_to);
  free_cond(markername_snp);
  free_cond(snps_flag_markers);
  free_cond(snps_flag_starts_range);
  free_cond(lgen_reference_fname);
  free_cond(covar_fname);
  free_cond(covar_str);
  free_cond(cluster_fname);
  free_cond(set_fname);
  free_cond(subset_fname);
  free_cond(update_alleles_fname);
  free_cond(update_map_fname);
  free_cond(update_ids_fname);
  free_cond(update_parents_fname);
  free_cond(update_sex_fname);
  free_cond(loop_assoc_fname);

  if (logfile) {
    if (!log_failed) {
      logstr("\nEnd time: ");
      time(&rawtime);
      logstr(ctime(&rawtime));
      if (fclose(logfile)) {
	fputs("Error: Failed to finish writing to log.\n", stdout);
      }
    } else {
      fclose(logfile);
    }
  }
  return retval;
}
