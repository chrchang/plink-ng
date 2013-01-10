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


// Uncomment this to build without CBLAS/CLAPACK:
// #define NOLAPACK

// The key ideas behind this calculator's design are:
//
// 1. Incremental processing of SNPs.  Each element A_{jk} of a distance or
// relationship matrix is a sum of M terms, one for each SNP, multiplied by a
// missingness correction at the end.  We can walk through the SNPs
// sequentially without keeping much in memory beyond partial sums;
// conveniently, this plays well with the decision made by the PLINK team a few
// years ago to switch to SNP-major binary files.
//
// 2. Multiplexing of markers using bitwise operations.  For instance, there
// are only seven possible ways SNP i can affect the relationship matrix entry
// between individuals j and k:
//    a. j and k are both homozygous for the rare allele
//    b. one is homozygous rare, one is heterozygous
//    c. one is homozygous rare, one is homozygous common
//    d. both are heterozygous
//    e. one is heterozygous, one is homozygous common
//    f. both are homozygous common
//    g. one or both has missing genotype data
// Seven cases can be distinguished by three bits, so one can compose a
// sequence of bitwise operations that maps a pair of padded 2-bit PLINK
// genotypes to seven different 3-bit values according to this breakdown.
// On 64-bit machines, this allows integer operations to act on 20 markers
// simultaneously.  (There's space for a 21st, but we currently choose not to
// use it.)
//
// 3. Lookup tables describing the effect of 5-7 markers at a time on a
// distance or relationship, rather than just one.  For relationship matrices,
// idea #2 allows computation of a single 64-bit integer where bits 0-2
// describe the relationship on marker #0, bits 3-5 describe the relationship
// on marker #1, ..., all the way up to bits 57-59 describing the relationship
// on marker #19.  We then want to perform the update
//    A_{jk} := A_{jk} + f_0(x_0) + f_1(x_1) + ... + f_19(x_19)
// where the x_i's are bit trios, and the f_i's map them to floating point
// terms.  We could do this with 20 table lookups and floating point additions.
// Or, we could structure this as
//    A_{jk} := A_{jk} + f_{0-4}(x_{0-4}) + ... + f_{15-19}(x_{15-19})
// where x_{0-4} denotes the lowest order *15* bits, and f_{0-4} maps them
// directly to f_0(x_0) + f_1(x_1) + f_2(x_2) + f_3(x_3) + f_4(x_4); similarly
// for f_{5-9}, f_{10-14}, and f_{15-19}.  This requires precomputation of four
// lookup tables of size 2^15, total size 1 MB (which fits comfortably in a
// typical L2 cache these days), and licenses the use of four table lookups and
// adds instead of twenty.
//
// 4. Bitslicing algorithms for especially fast calculation of unweighted
// distances and SNP covariances.  Zero-exponent distance matrices and IBS
// matrices are special cases, reducing to Hamming weight computations plus a
// bit of dressing to handle missing markers.  Hamming weight computation on
// x86 CPUs has been studied extensively by others; a good reference as of this
// writing is
//    http://www.dalkescientific.com/writings/diary/archive/2011/11/02/faster_popcount_update.html .
// We use a variant of the Kim Walisch/Cedric Lauradoux bitslicing algorithm
// discussed there (with most 64-bit integer operations replaced by SSE2
// instructions), which runs several times faster than our corresponding
// nonzero exponent distance matrix computation.
//
// We also reduce SNP covariance calculation (used in LD-based pruning) to a
// few Hamming weights.  (This can also be done for covariances between
// individuals, but only when there are no missing markers, so WDIST does not
// include an implementation of that.)
//
// 5. Splitting the distance/relationship matrix into pieces of roughly equal
// size and assigning one thread to each piece.  This is an "embarrassingly
// parallel" problem with no need for careful synchronization.  Cluster
// computing is supported in essentially the same manner.
//
//
//
// In the end, we can't get around the O(MN^2) nature of these calculations,
// but we believe we have beaten down the leading constant by a large enough
// factor to meaningfully help researchers.

#include <ctype.h>
#include <time.h>
#include <unistd.h>

#ifdef NOLAPACK
#define MATRIX_INVERT_BUF1_TYPE double
#define __CLPK_integer int
#else
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#else
// allow the same code to work for OS X and Linux
#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>
#if __LP64__
  typedef int __CLPK_integer;
#else
  typedef long int __CLPK_integer;
#endif
  typedef double __CLPK_doublereal;
  int dgetrf_(__CLPK_integer* m, __CLPK_integer* n,
              __CLPK_doublereal* a, __CLPK_integer* lda,
              __CLPK_integer* ipiv, __CLPK_integer* info);

  int dgetri_(__CLPK_integer* n, __CLPK_doublereal* a,
              __CLPK_integer* lda, __CLPK_integer* ipiv,
              __CLPK_doublereal* work, __CLPK_integer* lwork,
              __CLPK_integer* info);

  double dlange_(char* norm, __CLPK_integer* m, __CLPK_integer* n,
                 __CLPK_doublereal* a, __CLPK_integer* lda,
                 __CLPK_doublereal* work);

  int dgecon_(char* norm, __CLPK_integer* n, __CLPK_doublereal* a,
              __CLPK_integer* lda, __CLPK_doublereal* anorm,
              __CLPK_doublereal* rcond, __CLPK_doublereal* work,
              __CLPK_integer* iwork, __CLPK_integer* info);

  void xerbla_(void) {} // fix static linking error
#ifdef __cplusplus
}
#endif // __cplusplus
#endif // __APPLE__
#define MATRIX_INVERT_BUF1_TYPE __CLPK_integer
#endif // NOLAPACK

#include "wdist_common.h"
#include "wdist_data.h"
#include "wdist_dosage.h"

// floating point comparison-to-nonzero tolerance, currently 2^{-30}
#define EPSILON 0.000000000931322574615478515625
#define MATRIX_SINGULAR_RCOND 1e-14

#define REL_CALC_COV 1
#define REL_CALC_SQ 2
#define REL_CALC_SQ0 4
#define REL_CALC_TRI 6
#define REL_CALC_SHAPEMASK 6
#define REL_CALC_GZ 8
#define REL_CALC_BIN 16
#define REL_CALC_GRM 32

#define LOAD_RARE_GRM 1
#define LOAD_RARE_LGEN 2
#define LOAD_RARE_TRANSPOSE 4
#define LOAD_RARE_TPED 8
#define LOAD_RARE_TFAM 16
#define LOAD_RARE_TRANSPOSE_MASK (LOAD_RARE_TRANSPOSE | LOAD_RARE_TPED | LOAD_RARE_TFAM)
#define LOAD_RARE_DUMMY 32

#define PEDBUFBASE 256

// number of different types of jackknife values to precompute (x^2, x, y, xy)
#define JACKKNIFE_VALS_REL 5

#define MAX_EM_ACCEL 100.0

// maximum number of usable cluster computers, this is arbitrary though it
// shouldn't be larger than 2^31 - 1
#define PARALLEL_MAX 32768

// default jackknife iterations
#define ITERS_DEFAULT 100000
#define DEFAULT_PPC_GAP 500000
#define MAX_PCS_DEFAULT 20

// Number of snp-major .bed lines to read at once for distance calc if exponent
// is zero.  Currently assumed to be a multiple of 192, and no larger than
// 1920, by the popcount_..._multiword functions.  (The optimal value depends
// on both system-specific properties such as cache sizes, as well as the
// number of individuals in the current calculation, so in principle it's best
// to select this value at runtime.  But 960 usually works well in practice in
// my experience.)
#define MULTIPLEX_DIST 960
#define MULTIPLEX_2DIST (MULTIPLEX_DIST * 2)

#define MULTIPLEX_LD 1920
#define MULTIPLEX_2LD (MULTIPLEX_LD * 2)

// Must be multiple of 384, no larger than 3840.
#define GENOME_MULTIPLEX 1152
#define GENOME_MULTIPLEX2 (GENOME_MULTIPLEX * 2)

#if __LP64__
#define AAAAMASK 0xaaaaaaaaaaaaaaaaLU
// number of snp-major .bed lines to read at once for distance calc if exponent
// is nonzero.
#define MULTIPLEX_DIST_EXP 64
// number of snp-major .bed lines to read at once for relationship calc
#define MULTIPLEX_REL 60
#else
// N.B. 32-bit version not as carefully tested or optimized, but I'll try to
// make sure it works properly
#define AAAAMASK 0xaaaaaaaa
#define MULTIPLEX_DIST_EXP 28
#define MULTIPLEX_REL 30
#endif

#define MINV(aa, bb) ((aa) > (bb))? (bb) : (aa)

const char ver_str[] =
  "WDIST v0.15.0"
#ifdef NOLAPACK
  "NL"
#endif
#if __LP64__
  " 64-bit"
#else
  " 32-bit"
#endif
  " (11 Jan 2013)";
const char ver_str2[] =
  "    https://www.cog-genomics.org/wdist\n"
  "(C) 2013 Christopher Chang, GNU General Public License version 3\n";
// const char errstr_append[] = "\nFor more information, try 'wdist --help {flag names}' or 'wdist --help | more'.\n";
const char errstr_ped_format[] = "Error: Improperly formatted .ped file.\n";
const char errstr_phenotype_format[] = "Error: Improperly formatted phenotype file.\n";
const char errstr_filter_format[] = "Error: Improperly formatted filter file.\n";
const char errstr_freq_format[] = "Error: Improperly formatted frequency file.\n";
const char cmdline_format_str[] = "\nwdist --help {flags...}\n"
  "wdist [input flags...] [commands...] {other flags...}\n\n";
const char notestr_null_calc[] = "Note: No output requested.  Exiting.\n";
const char notestr_null_calc2[] = "Commands include --freqx, --ibc, --distance, --genome, --make-rel, --make-grm,\n--rel-cutoff, --regress-distance, --regress-pcs-distance, --make-bed, --recode,\n--merge-list, and --write-snplist.\n\n'wdist --help | more' describes all functions (warning: long).  You can look up\nspecific flags with 'wdist --help [flag #1] {flag #2} ...'.\n";

int edit1_match(int len1, char* s1, int len2, char* s2) {
  // permit one difference of the following forms:
  // - inserted/deleted character
  // - replaced character
  // - adjacent pair of swapped characters
  int diff_found = 0;
  int pos = 0;
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

#define MAX_EQUAL_HELP_PARAMS 6

typedef struct {
  int iters_left;
  unsigned int param_ct;
  char** argv;
  unsigned int unmatched_ct;
  unsigned long* all_match_arr;
  unsigned long* prefix_match_arr;
  unsigned long* perfect_match_arr;
  unsigned int* param_lens;
  int preprint_newline;
} Help_ctrl;

void help_print(const char* cur_params, Help_ctrl* help_ctrl_ptr, int postprint_newline, const char* payload) {
  // unmatched_ct fixed during call, *unmatched_ct_ptr may decrease
  unsigned int unmatched_ct = help_ctrl_ptr->unmatched_ct;
  int print_this = 0;
  int cur_param_lens[MAX_EQUAL_HELP_PARAMS];
  char* cur_param_start[MAX_EQUAL_HELP_PARAMS];
  int cur_param_ct;
  int cur_param_idx;
  unsigned int arg_uidx;
  unsigned int arg_idx;
  int uii;
  int payload_len;
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
	  uii = (unsigned int)(line_end - payload_ptr);
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

int disp_help(unsigned int param_ct, char** argv) {
  // yes, this is overkill.  But it should be a good template for other
  // command-line programs to use.
  unsigned int param_ctl = (param_ct + (BITCT - 1)) / BITCT;
  int retval = 0;
  Help_ctrl help_ctrl;
  unsigned int arg_uidx;
  unsigned int arg_idx;
  unsigned int net_unmatched_ct;
  int col_num;
  int leading_dashes;
  help_ctrl.iters_left = param_ct? 2 : 0;
  help_ctrl.param_ct = param_ct;
  help_ctrl.argv = argv;
  help_ctrl.unmatched_ct = param_ct;
  help_ctrl.param_lens = NULL;
  help_ctrl.all_match_arr = NULL;
  help_ctrl.argv = NULL;
  if (param_ct) {
    help_ctrl.param_lens = (unsigned int*)malloc(param_ct * sizeof(int));
    if (!help_ctrl.param_lens) {
      goto disp_help_ret_NOMEM;
    }
    help_ctrl.all_match_arr = (unsigned long*)malloc(param_ctl * 3 * sizeof(long));
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
    fputs(cmdline_format_str, stdout);
    fputs(
"In the command line flag definitions that follow,\n"
"  * [square brackets] denote a required parameter, where the text between the\n"
"    brackets describes its nature.\n"
"  * <angle brackets> denote an optional modifier (or if '|' is present, a set\n"
"    of mutually exclusive optional modifiers).  Use the EXACT text in the\n"
"    definition, e.g. '--distance square0'.\n"
"  * {curly braces} denote an optional parameter, where the text between the\n"
"    braces describes its nature.\n"
"  * An ellipsis (...) indicates that you may enter multiple parameters of the\n"
"    specified type.\n\n"
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
"  --grm {prefix}   : Load a GCTA relationship matrix for --rel-cutoff (default\n"
"                     filename prefix 'wdist').\n\n"
	       );
    help_print("dummy", &help_ctrl, 1,
"  --dummy [marker ct] [indiv ct] {missing geno freq} {missing pheno freq}\n"
"          <acgt | 1234 | 12> <scalar-pheno>\n"
"    Create a new dataset with the specified number of markers and individuals.\n"
"    By default, the missing genotype and phenotype frequencies are zero, and\n"
"    genotypes are As and Bs (change the latter with 'acgt'/'1234'/'12').  The\n"
"    'scalar-pheno' modifier causes a normally distributed scalar phenotype to\n"
"    be generated instead of a binary one.\n\n"
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
"    Identity-by-descent analysis.  This yields the same output as PLINK's\n"
"    --genome or --Z-genome, and the 'full' and 'unbounded' modifiers have the\n"
"    same effect as PLINK's --genome-full and --unbounded flags.\n\n"
		);
    help_print(
"indep\tindep-pairwise", &help_ctrl, 1,
"  --indep [window size]<kb> [step size (SNPs)] [VIF threshold]\n"
"  --indep-pairwise [window size]<kb> [step size (SNPs)] [r^2 threshold]\n"
"    Generates a list of SNPs in approximate linkage equilibrium.  With the 'kb'\n"
"    modifier, the window size is in kilobase units instead of SNPs.  (Space\n"
"    before 'kb' is optional, i.e. '--indep-pairwise 500 kb 5 0.5' and\n"
"    '--indep-pairwise 500kb 5 0.5' have the same effect.)\n"
"    Note that you need to rerun WDIST using --extract or --exclude on the\n"
"    .prune.in/.prune.out file to apply the list to another computation.\n\n"
		);
    help_print("make-rel", &help_ctrl, 1,
"  --make-rel <square | square0 | triangle> <gz | bin> <cov | ibc1 | ibc2>\n"
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
               );
    help_print("make-grm\tgrm", &help_ctrl, 1,
"  --make-grm <no-gz> <cov | ibc1 | ibc2>\n"
"    Writes the relationships in GCTA's gzipped list format, describing one pair\n"
"    per line.  Note that this file explicitly stores the number of valid\n"
"    observations (where neither individual has a missing call) for each pair,\n"
"    which is useful input for some scripts.\n\n"
	       );
    help_print("rel-cutoff\tgrm-cutoff\tgrm", &help_ctrl, 1,
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
    help_print("recode\trecode12\ttab\ttranspose\trecode-lgen", &help_ctrl, 1,
"  --recode <12> <tab | tabx | spacex> <transpose | lgen>\n"
"    Creates a new text fileset with all filters applied.\n"
"    * The '12' modifier causes all alleles to be coded as 1s and 2s.\n"
"    * The 'tab' modifier makes the output mostly tab-delimited instead of\n"
"       mostly space-delimited.  'tabx' and 'spacex' force all tabs and all\n"
"       spaces, respectively.\n"
"    * The 'transpose' modifier causes a transposed text fileset to be generated\n"
"      instead.\n"
"    * The 'lgen' modifier causes a long-format fileset to be generated instead.\n\n"
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
"    Writes a .snplist file listing the names of all SNPs that pass the filters\n"
"    and inclusion thresholds you've specified.\n\n"
	       );
    help_print("help\t{flag names}", &help_ctrl, 1,
"  --help {flag names...}\n"
"    Describes WDIST functionality.  {flag names...} is a list of flag names or\n"
"    prefixes, e.g. 'wdist --help filter- missing-' will display information\n"
"    about all flags with names starting with 'filter-' or 'missing-', while\n"
"    'wdist --help filter' will just describe --filter since that's an exact\n"
"    match.  If no parameters are given, all commands are listed.\n\n"
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
"                     For now, only --merge/--bmerge/--merge-list and explicit\n"
"                     --make-bed respect this flag.\n"
	       );
    help_print("pheno", &help_ctrl, 0,
"  --pheno [fname]  : Specify alternate phenotype.\n"
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
    help_print("cow\tdog\thorse\tmouse\tsheep", &help_ctrl, 0,
"  --cow/--dog/--horse/--mouse/--sheep : Specify nonhuman species.\n"
	       );
    help_print("autosome", &help_ctrl, 0,
"  --autosome       : Include markers on all autosomes, and no others.\n"
	       );
    help_print("chr", &help_ctrl, 0,
"  --chr [num...]   : Only consider markers on the given chromosome(s).  Valid\n"
"                     choices for humans are 0 (unplaced), 1-22, and XY\n"
"                     (pseudo-autosomal region of X).  Separate multiple\n"
"                     chromosomes with spaces, e.g. '--chr 1 3 xy'.  The X, Y,\n"
"                     and MT chromosome values are currently unsupported, but\n"
"                     this will change in the future.\n"
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
"  --geno {val}     : Maximum per-SNP missing (default 0.1).\n"
"  --mind {val}     : Maximum per-person missing (default 0.1).\n"
	       );
    help_print("hwe", &help_ctrl, 0,
"  --hwe {val}      : Minimum Hardy-Weinberg disequilibrium p-value (exact),\n"
"                     default 0.001.\n"
	       );
    help_print("hwe-all\thwe", &help_ctrl, 0,
"  --hwe-all        : Given case-control data, don't ignore cases in HWE test.\n"
	       );
    help_print("allow-no-sex", &help_ctrl, 0,
"  --allow-no-sex   : Do not remove ambiguously-sexed individuals.\n"
"                     (Automatically on if --no-sex flag is present.)\n"
	       );
    help_print("nonfounders", &help_ctrl, 0,
"  --nonfounders    : Include all individuals in MAF/HWE calculations.\n"
	       );
    help_print("ppc-gap", &help_ctrl, 0,
"  --ppc-gap [val]  : Minimum number of base pairs, in thousands, between\n"
"                     informative pairs of SNPs used in --genome PPC test.  500\n"
"                     if unspecified.\n"
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
"  --extract [file] : Exclude all SNPs not in the given list.\n"
"  --exclude [file] : Exclude all SNPs in the given list.\n"
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
    help_print("allele1234\talleleACGT\talleleacgt", &help_ctrl, 0,
"  --allele1234     : Interpret/recode A/C/G/T alleles as 1/2/3/4.\n"
"  --alleleACGT     : Interpret/recode 1/2/3/4 alleles as A/C/G/T.\n"
	       );
    help_print("exponent\tdistance", &help_ctrl, 0,
"  --exponent [val] : When computing genomic distances, each marker has a weight\n"
"                     of (2q(1-q))^{-val}, where q is the inferred MAF.  (Use\n"
"                     --read-freq if you want to explicitly specify some or all\n"
"                     of the MAFs.)\n"
	       );
    help_print("keep-allele-order\tmake-bed\tmerge\tbmerge\tmerge-list", &help_ctrl, 0,
"  --keep-allele-order : Keep the original allele order when creating a new\n"
"                        fileset, instead of forcing A2 to be the major allele.\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode\tmerge-allow-equal-pos", &help_ctrl, 0,
"  --merge-allow-equal-pos   : Do not merge markers with different names but\n"
"                              identical positions.\n"
	       );
    help_print("reference-allele", &help_ctrl, 0,
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
"                              treated as missing.\n\n"
	       );
    if (!param_ct) {
      fputs(
"Further documentation and support is available at the main webpage\n"
"(https://www.cog-genomics.org/wdist ) and the wdist-users mailing list\n"
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

inline void gzclose_cond(gzFile gz_outfile) {
  if (gz_outfile) {
    gzclose(gz_outfile);
  }
}

long malloc_size_mb = 0;

unsigned char* wkspace;

void dispmsg(int retval) {
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

// (copied from PLINK helper.cpp, modified to directly perform threshold test
// and avoid repeated allocation of the same memory.  Substantial further
// optimization should be possible, but it's unimportant for now.)
//
//
// This function implements an exact SNP test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton

double* het_probs = NULL;

int SNPHWE_t(int obs_hets, int obs_hom1, int obs_hom2, double thresh) {
  if (!(obs_hom1 + obs_hom2 + obs_hets)) {
    return 0;
  }
  
  // if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) 
  //   {
  //     error("Internal error: negative count in HWE test: "
  //           +int2str(obs_hets)+" "
  //           +int2str(obs_hom1)+" "
  //           +int2str(obs_hom2));
  //   }

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int rare_copies2 = rare_copies / 2;
  int genotypes   = obs_hets + obs_homc + obs_homr;
  // avoid integer overflow when calculating mid
  long long rare_copies_ll = rare_copies;
  long long genotypes_ll = genotypes;

  int i;
  for (i = 0; i <= rare_copies2; i++) {
    het_probs[i] = 0.0;
  }

  // start at midpoint
  int mid = (int)((rare_copies_ll * (2 * genotypes_ll - rare_copies_ll)) / (2 * genotypes_ll));
  
  // check to ensure that midpoint and rare alleles have same parity
  if ((rare_copies ^ mid) & 1) {
    mid++;
  }
  int mid2 = mid / 2;
  
  long long curr_hets = mid;

  het_probs[mid2] = 1.0;
  double sum = 1.0;
  long long curr_homr = rare_copies2 - mid2;
  int curr_homc = genotypes - mid - curr_homr;
  double* hptr = &(het_probs[mid2 - 1]);
  double dv = 1.0;
  while (curr_hets > 1) {
    // 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
    *hptr = dv * (double)(curr_hets * (curr_hets - 1)) / (double)(4 * (++curr_homr) * (++curr_homc));
    curr_hets -= 2;
    dv = *hptr--;
    sum += dv;
  }

  curr_hets = mid;
  curr_homr = rare_copies2 - mid2;
  curr_homc = genotypes - mid - curr_homr;
  hptr = &(het_probs[mid2 + 1]);
  dv = 1.0;
  while (curr_hets <= rare_copies - 2) {
    curr_hets += 2;
    *hptr = dv * (double)(4 * (curr_homr--) * (curr_homc--)) / (double)(curr_hets * (curr_hets - 1));
    dv = *hptr++;
    sum += dv;
  }

  thresh *= sum;

  double p_hwe = 0.0;
  hptr = &(het_probs[rare_copies2]);
  sum = het_probs[obs_hets / 2];
  // p-value calculation for p_hwe
  while (*hptr <= sum) {
    p_hwe += *hptr--;
    if (p_hwe > thresh) {
      return 0;
    }
  }
  hptr = het_probs;
  while (*hptr <= sum) {
    p_hwe += *hptr++;
  }
  return (!(p_hwe > thresh));
}

// back to our regular program

inline double get_maf(double allele_freq) {
  if (allele_freq < 0.5) {
    return allele_freq;
  } else {
    return (1.0 - allele_freq);
  }
}

int strcmp_casted(const void* s1, const void* s2) {
  return strcmp((char*)s1, (char*)s2);
}

int llcmp(const void* aa, const void* bb) {
  long long diff = *((const long long*)aa) - *((const long long*)bb);
  if (diff > 0) {
    return 1;
  } else if (diff < 0) {
    return -1;
  } else {
    return 0;
  }
}

double get_dmedian_i(int* sorted_arr, int len) {
  if (len) {
    if (len % 2) {
      return (double)sorted_arr[len / 2];
    } else {
      return ((double)sorted_arr[len / 2] + (double)sorted_arr[(len / 2) - 1]) * 0.5;
    }
  } else {
    return 0.0;
  }
}

double get_dmedian(double* sorted_arr, int len) {
  if (len) {
    if (len % 2) {
      return sorted_arr[len / 2];
    } else {
      return (sorted_arr[len / 2] + sorted_arr[(len / 2) - 1]) * 0.5;
    }
  } else {
    return 0.0;
  }
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

inline int is_contained(char* id_buf, char* lptr, int max_id_len, int filter_line_ct, char* fam_id, char* indiv_id) {
  return (bsearch_fam_indiv(id_buf, lptr, max_id_len, filter_line_ct, fam_id, indiv_id) != -1);
}

int determine_max_id_len(FILE* filterfile, char* filterval, int mfilter_col, int* filter_line_ct_ptr) {
  int cur_max = 4;
  char* bufptr;
  int ii;
  int jj;

  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Excessively long line in filter file (max %d chars).\n", MAXLINELEN - 3);
      logprintb();
      return -1;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln(*bufptr)) {
      continue;
    }
    ii = 2 + strlen_se(bufptr);
    bufptr = next_item(bufptr);
    if (no_more_items(bufptr)) {
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
      if (no_more_items(bufptr)) {
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

char* resize_id_buf(char* id_buf, int max_id_len, int max_pid_len) {
  if (max_pid_len) {
    if (max_id_len > max_pid_len) {
      free(id_buf);
    } else {
      return id_buf;
    }
  }
  return (char*)malloc(max_id_len * sizeof(char));
}

void fill_weights(double* weights, double* set_allele_freqs, double exponent) {
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int oo;
  double wtarr[MULTIPLEX_DIST_EXP / 2];
  double* wt;
#if __LP64__
  double twt[5];
  double twtf;
  __m128d* wpairs = (__m128d*)weights;
  __m128d vpen;
  __m128d vfinal1;
  __m128d vfinal2;
#else
  int pp;
  int qq;
  double twt[7];
#endif
  for (ii = 0; ii < MULTIPLEX_DIST_EXP / 2; ii++) {
    wtarr[ii] = pow(2 * set_allele_freqs[ii] * (1.0 - set_allele_freqs[ii]), -exponent);
  }
  for (oo = 0; oo < 2; oo++) {
    wt = &(wtarr[7 * oo]);
#if __LP64__
    vfinal1 = _mm_set_pd(wt[0], 0.0);
    vfinal2 = _mm_set_pd(wt[0] * 2, wt[0]);
#endif
    twt[0] = 0;
    for (ii = 0; ii < 4; ii += 1) {
      // bizarrely, turning the ii == 2 case into a memcpy doesn't actually
      // seem to speed this up
      if (ii & 1) {
	twt[0] += wt[6];
      }
      twt[1] = twt[0];
      for (jj = 0; jj < 4; jj += 1) {
	if (jj & 1) {
	  twt[1] += wt[5];
	}
	twt[2] = twt[1];
	for (kk = 0; kk < 4; kk += 1) {
	  if (kk & 1) {
	    twt[2] += wt[4];
	  }
	  twt[3] = twt[2];
	  for (mm = 0; mm < 4; mm++) {
	    if (mm & 1) {
	      twt[3] += wt[3];
	    }
	    twt[4] = twt[3];
	    for (nn = 0; nn < 4; nn++) {
	      if (nn & 1) {
		twt[4] += wt[2];
	      }
#if __LP64__
	      twtf = twt[4];
	      vpen = _mm_set1_pd(twtf);
	      *wpairs++ = _mm_add_pd(vpen, vfinal1);
	      *wpairs++ = _mm_add_pd(vpen, vfinal2);
	      twtf += wt[1];
	      vpen = _mm_set1_pd(twtf);
	      *wpairs++ = _mm_add_pd(vpen, vfinal1);
	      *wpairs++ = _mm_add_pd(vpen, vfinal2);
	      *wpairs = *(wpairs - 2);
	      wpairs++;
	      *wpairs = *(wpairs - 2);
	      wpairs++;
	      vpen = _mm_set1_pd(twtf + wt[1]);
	      *wpairs++ = _mm_add_pd(vpen, vfinal1);
	      *wpairs++ = _mm_add_pd(vpen, vfinal2);
#else
              twt[5] = twt[4];
              for (pp = 0; pp < 4; pp++) {
                if (pp & 1) {
                  twt[5] += wt[1];
                }
                twt[6] = twt[5];
                for (qq = 0; qq < 4; qq++) {
                  if (qq & 1) {
                    twt[6] += wt[0];
                  }
                  *weights++ = twt[6];
                }
              }
#endif
	    }
          }
	}
      }
    }
  }
#if __LP64__
  for (oo = 0; oo < 3; oo++) {
    wt = &(wtarr[14 + 6 * oo]);
    vfinal1 = _mm_set_pd(wt[0], 0.0);
    vfinal2 = _mm_set_pd(2 * wt[0], wt[0]);
    twt[0] = 0;
    for (ii = 0; ii < 4; ii += 1) {
      if (ii & 1) {
        twt[0] += wt[5];
      }
      twt[1] = twt[0];
      for (jj = 0; jj < 4; jj += 1) {
        if (jj & 1) {
          twt[1] += wt[4];
        }
        twt[2] = twt[1];
	for (kk = 0; kk < 4; kk += 1) {
          if (kk & 1) {
            twt[2] += wt[3];
          }
          twt[3] = twt[2];
          for (mm = 0; mm < 4; mm++) {
	    if (mm & 1) {
	      twt[3] += wt[2];
	    }
	    twtf = twt[3];
	    vpen = _mm_set1_pd(twtf);
	    *wpairs++ = _mm_add_pd(vpen, vfinal1);
	    *wpairs++ = _mm_add_pd(vpen, vfinal2);
	    twtf += wt[1];
	    vpen = _mm_set1_pd(twtf);
	    *wpairs++ = _mm_add_pd(vpen, vfinal1);
	    *wpairs++ = _mm_add_pd(vpen, vfinal2);
	    *wpairs = *(wpairs - 2);
	    wpairs++;
	    *wpairs = *(wpairs - 2);
	    wpairs++;
	    vpen = _mm_set1_pd(twtf + wt[1]);
	    *wpairs++ = _mm_add_pd(vpen, vfinal1);
	    *wpairs++ = _mm_add_pd(vpen, vfinal2);
          }
	}
      }
    }
  }
#endif
}

void fill_weights_r(double* weights, double* set_allele_freqs, int var_std) {
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  // 20 markers to process in quintuplets, for 64-bit; 10, for 32-bit.
  // Each quintuplet of markers requires 40 wtarr entries, and induces
  // 2^15 writes to weights[].
  double wtarr_raw[BITCT2 * 5 + 1];
  double* wtarr = wtarr_raw;
  double twt;
  double twt2;
  double twt3;
  double twt4;
  double* wtptr;
  double mean;
  double mean_m1;
  double mean_m2;
  double mult = 1.0;
  double aux;
#if __LP64__
  __m128d* wpairs = (__m128d*)weights;
  __m128d vpen;
  __m128d vfinal1;
  __m128d vfinal2;
  __m128d vfinal3;
  __m128d vfinal4;
#else
  int oo;
#endif
  if (((unsigned long)wtarr) & 15) {
    // force 16-byte alignment; can't do this at compile-time since stack
    // pointer has no 16-byte align guarantee.
    // yes, this assumes doubles are 8 bytes.
    wtarr++;
  }
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if (((set_allele_freqs[ii] != 0.0) && (set_allele_freqs[ii] < (1.0 - EPSILON))) || (!var_std)) {
      if (set_allele_freqs[ii] < 0.5) {
	mean = 2 * set_allele_freqs[ii];
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * (1.0 - set_allele_freqs[ii]));
        }
        aux = mean * mult;
	wtarr[ii * 8] = mean * aux;
        wtarr[ii * 8 + 1] = 0;
	wtarr[ii * 8 + 2] = mean_m1 * aux;
	wtarr[ii * 8 + 3] = mean_m2 * aux;
	wtarr[ii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[ii * 8 + 5] = mean_m2 * mean_m1 * mult;
	wtarr[ii * 8 + 6] = mean_m2 * mean_m2 * mult;
      } else {
	mean = 2 * (1.0 - set_allele_freqs[ii]);
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * set_allele_freqs[ii]);
        }
        aux = mean_m2 * mult;
	wtarr[ii * 8] = mean_m2 * aux;
        wtarr[ii * 8 + 1] = 0;
	wtarr[ii * 8 + 2] = mean_m1 * aux;
	wtarr[ii * 8 + 3] = mean * aux;
	wtarr[ii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[ii * 8 + 5] = mean_m1 * mean * mult;
	wtarr[ii * 8 + 6] = mean * mean * mult;
      }
    } else {
      if (set_allele_freqs[ii] == 0.0) {
        wtarr[ii * 8] = 0;
        wtarr[ii * 8 + 1] = 0;
        wtarr[ii * 8 + 2] = -1;
        wtarr[ii * 8 + 3] = -2;
        wtarr[ii * 8 + 4] = INFINITY;
        wtarr[ii * 8 + 5] = INFINITY;
        wtarr[ii * 8 + 6] = INFINITY;
      } else {
        wtarr[ii * 8] = INFINITY;
        wtarr[ii * 8 + 1] = 0;
        wtarr[ii * 8 + 2] = INFINITY;
        wtarr[ii * 8 + 3] = -2;
        wtarr[ii * 8 + 4] = INFINITY;
        wtarr[ii * 8 + 5] = -1;
        wtarr[ii * 8 + 6] = 0;
      }
    }
    wtarr[ii * 8 + 7] = 0;
  }
  for (nn = 0; nn < BITCT / 16; nn++) {
    wtptr = &(wtarr[40 * nn]);
#if __LP64__
    vfinal1 = _mm_load_pd(wtptr);
    vfinal2 = _mm_load_pd(&(wtptr[2]));
    vfinal3 = _mm_load_pd(&(wtptr[4]));
    vfinal4 = _mm_load_pd(&(wtptr[6]));
#endif
    for (ii = 0; ii < 8; ii++) {
      twt = wtptr[ii + 32];
      for (jj = 0; jj < 8; jj++) {
        twt2 = twt + wtptr[jj + 24];
        for (kk = 0; kk < 8; kk++) {
          twt3 = twt2 + wtptr[kk + 16];
          for (mm = 0; mm < 8; mm++) {
            twt4 = twt3 + wtptr[mm + 8];
#if __LP64__
            vpen = _mm_set1_pd(twt4);
            *wpairs++ = _mm_add_pd(vpen, vfinal1);
            *wpairs++ = _mm_add_pd(vpen, vfinal2);
            *wpairs++ = _mm_add_pd(vpen, vfinal3);
            *wpairs++ = _mm_add_pd(vpen, vfinal4);
#else
            for (oo = 0; oo < 8; oo++) {
              *weights++ = twt4 + wtptr[oo];
            }
#endif
          }
        }
      }
    }
  }
}

double calc_wt_mean(double exponent, int lhi, int lli, int hhi) {
  double lcount = (double)lli + ((double)lhi * 0.5);
  long long tot = lhi + lli + hhi;
  double dtot = (double)tot;
  long long subcount = lli; // avoid 32-bit integer overflow
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

// ----- multithread globals -----
static unsigned int g_indiv_ct;
static unsigned int g_thread_ct;
static double* g_rel_dists = NULL;
static unsigned int* g_missing_dbl_excluded = NULL;
static int* g_idists;
static double* g_dists = NULL;
static char* g_pheno_c = NULL;
static double* g_pheno_d = NULL;
static unsigned char* g_geno = NULL;
// see incr_dists()
#if (2048 * BITCT) < 45056
#error "Insufficient initial size for g_weights[]."
#endif
static double g_weights[2048 * BITCT];
static unsigned int* g_weights_i = (unsigned int*)g_weights;
static unsigned int g_thread_start[MAX_THREADS_P1];
static double g_reg_tot_xy;
static double g_reg_tot_x;
static double g_reg_tot_y;
static double g_reg_tot_xx;
static double g_reg_tot_yy;
static unsigned int g_low_ct;
static unsigned int g_high_ct;
static unsigned long g_jackknife_iters;
static unsigned int g_jackknife_d;
static double g_calc_result[9][MAX_THREADS];
static unsigned long* g_masks;
static unsigned long* g_mmasks;
static double* g_marker_weights = NULL;
static unsigned int* g_marker_weights_i = NULL;
static unsigned int* g_missing_tot_weights;
static unsigned int* g_indiv_missing;
static unsigned int* g_indiv_missing_unwt = NULL;
static double* g_jackknife_precomp = NULL;
static unsigned int* g_genome_main;
static unsigned long g_marker_window[GENOME_MULTIPLEX * 2];
static double* g_pheno_packed;

void update_rel_ibc(double* rel_ibc, unsigned long* geno, double* set_allele_freqs, int ibc_type) {
  // first calculate weight array, then loop
  int ii;
  int jj;
  int kk;
  unsigned int uii;
  double twt;
  double* wtptr;
  double mult = 1.0;
  unsigned long ulii;
  double weights[BITCT * 10];
  double* weights1 = &(weights[64]);
  double* weights2 = &(weights[128]);
  double* weights3 = &(weights[192]);
  double* weights4 = &(weights[256]);
#if __LP64__
  double* weights5 = &(weights[320]);
  double* weights6 = &(weights[384]);
  double* weights7 = &(weights[448]);
  double* weights8 = &(weights[512]);
  double* weights9 = &(weights[576]);
#endif
  double wtarr[BITCT2 * 5];
  fill_double_zero(wtarr, BITCT2 * 5);
  double *wptr = weights;
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if ((set_allele_freqs[ii] != 0.0) && (set_allele_freqs[ii] < (1.0 - EPSILON))) {
      if (ibc_type) {
        if (ibc_type == 2) {
          wtarr[ii * 8] = 2;
          wtarr[ii * 8 + 2] = 2.0 - 1.0 / (2 * set_allele_freqs[ii] * (1.0 - set_allele_freqs[ii]));
          wtarr[ii * 8 + 3] = 2;
        } else {
          twt = 2 * set_allele_freqs[ii];
          if (ibc_type == 1) {
            mult = 1 / (twt * (1.0 - set_allele_freqs[ii]));
          }
          wtarr[ii * 8] = twt * twt * mult;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt) * mult;
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt) * mult;
        }
      } else {
        twt = 1.0 - set_allele_freqs[ii];
        mult = 1 / (set_allele_freqs[ii] * twt);
        wtarr[ii * 8] = 1.0 + set_allele_freqs[ii] * set_allele_freqs[ii] * mult;
        wtarr[ii * 8 + 3] = 1.0 + twt * twt * mult;
      }
    } else {
      if (ibc_type) {
        if (ibc_type == -1) {
          twt = 2 * set_allele_freqs[ii];
          wtarr[ii * 8] = twt * twt;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt);
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt);
        } else if (ibc_type == 1) {
	  wtarr[ii * 8 + 2] = INFINITY;
          if (set_allele_freqs[ii] == 0.0) {
            wtarr[ii * 8] = 0;
            wtarr[ii * 8 + 3] = INFINITY;
          } else {
            wtarr[ii * 8] = INFINITY;
            wtarr[ii * 8 + 3] = 0;
          }
        } else {
          // need to set to 1 instead of 2 for agreement with GCTA
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 2] = -INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      } else {
        if (set_allele_freqs[ii] == 0.0) {
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 3] = INFINITY;
        } else {
          wtarr[ii * 8] = INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      }
    }
  }
  for (kk = 0; kk < (BITCT * 5) / 32; kk++) {
    wtptr = &(wtarr[16 * kk]);
    for (ii = 0; ii < 8; ii++) {
      twt = wtptr[ii + 8];
      for (jj = 0; jj < 8; jj++) {
        *wptr++ = twt + wtptr[jj];
      }
    }
  }
  for (uii = 0; uii < g_indiv_ct; uii++) {
    ulii = *geno++;
#if __LP64__
    *rel_ibc += weights9[ulii >> 54] + weights8[(ulii >> 48) & 63] + weights7[(ulii >> 42) & 63] + weights6[(ulii >> 36) & 63] + weights5[(ulii >> 30) & 63] + weights4[(ulii >> 24) & 63] + weights3[(ulii >> 18) & 63] + weights2[(ulii >> 12) & 63] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#else
    *rel_ibc += weights4[ulii >> 24] + weights3[(ulii >> 18) & 63] + weights2[(ulii >> 12) & 63] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#endif
    rel_ibc++;
  }
}

// The following functions may come in handy if we need to store large array(s)
// of auxiliary ternary data in memory.
//
// inline void set_base4(unsigned long* base4_arr, int loc, int val) {
//   int shift_num = (loc % BITCT2) * 2;
//   unsigned long ulii = 3LU << shift_num;
//   unsigned long* base4_ptr = &(base4_arr[loc / BITCT2]);
//   *base4_ptr = (*base4_ptr & (~ulii)) | (val << shift_num);
// }
//
// inline int get_base4(unsigned long* base4_arr, int loc) {
//   int shift_num = (loc % BITCT2) * 2;
//   return ((base4_arr[loc / BITCT2] >> shift_num) & 3LU);
// }

int chrom_exists(Chrom_info* chrom_info_ptr, unsigned int chrom_idx) {
  return (chrom_info_ptr->chrom_mask & (1LLU << chrom_idx));
}

int get_marker_chrom(Chrom_info* chrom_info_ptr, unsigned int marker_uidx) {
  unsigned int* marker_binsearch = chrom_info_ptr->chrom_file_order_marker_idx;
  int chrom_min = 0;
  int chrom_ct = chrom_info_ptr->chrom_ct;
  int chrom_cur;
  while (chrom_ct - chrom_min > 1) {
    chrom_cur = (chrom_ct + chrom_min) / 2;
    if (marker_binsearch[chrom_cur] > marker_uidx) {
      chrom_ct = chrom_cur;
    } else {
      chrom_min = chrom_cur;
    }
  }
  return chrom_info_ptr->chrom_file_order[chrom_min];
}

int get_chrom_end(Chrom_info* chrom_info_ptr, unsigned int marker_idx) {
  return chrom_info_ptr->chrom_end[get_marker_chrom(chrom_info_ptr, marker_idx)];
}

inline int nz_chrom(Chrom_info* chrom_info_ptr, unsigned int marker_idx) {
  return (marker_idx >= chrom_info_ptr->chrom_end[0]) || (marker_idx < chrom_info_ptr->chrom_start[0]);
}

void exclude_multi(unsigned long* exclude_arr, int* new_excl, unsigned int indiv_ct, unsigned int* exclude_ct_ptr) {
  unsigned int uii;
  int true_loc = 0;
  for (uii = 0; uii < indiv_ct; uii++) {
    true_loc = next_non_set_unsafe(exclude_arr, true_loc);
    if (new_excl[uii] == -1) {
      set_bit(exclude_arr, true_loc, exclude_ct_ptr);
    }
    true_loc++;
  }
}

void collapse_bitarr(unsigned long* bitarr, unsigned long* exclude_arr, int orig_ct) {
  int ii = 0;
  int jj;
  while ((ii < orig_ct) && (!is_set(exclude_arr, ii))) {
    ii++;
  }
  jj = ii;
  while (++ii < orig_ct) {
    if (!is_set(exclude_arr, ii)) {
      if (is_set(bitarr, ii)) {
        // set bit jj
        bitarr[jj / BITCT] |= (1LU << (jj % BITCT));
      } else {
	bitarr[jj / BITCT] &= (~(1LU << (jj % BITCT)));
      }
      jj++;
    }
  }
}

void collapse_chrom_marker_idxs(Chrom_info* chrom_info_ptr, unsigned long* marker_exclude, int unfiltered_marker_ct) {
  unsigned int* chrom_fo = chrom_info_ptr->chrom_file_order;
  unsigned int* chrom_fo_midxs = chrom_info_ptr->chrom_file_order_marker_idx;
  unsigned int chrom_ct = chrom_info_ptr->chrom_ct;
  int midx = 0;
  int new_midx = 0;
  int chrom_end_midx;
  unsigned int chrom_fo_idx;
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    chrom_fo_midxs[chrom_fo_idx] = new_midx;
    chrom_info_ptr->chrom_start[chrom_fo[chrom_fo_idx]] = new_midx;
    chrom_end_midx = chrom_fo_midxs[chrom_fo_idx + 1];
    for (; midx < chrom_end_midx; midx++) {
      if (!is_set(marker_exclude, midx)) {
	new_midx++;
      }
    }
    // todo: collapse when chromosome completely eliminated
    chrom_info_ptr->chrom_end[chrom_fo[chrom_fo_idx]] = new_midx;
  }
  chrom_fo_midxs[chrom_fo_idx] = new_midx;
}

void collapse_copy_phenod(double *target, double* pheno_d, unsigned long* indiv_exclude, unsigned int indiv_ct) {
  int ii = 0;
  double* target_end = &(target[indiv_ct]);
  while (target < target_end) {
    ii = next_non_set_unsafe(indiv_exclude, ii);
    *target++ = pheno_d[ii++];
  }
}

#if __LP64__
// XOR + mask variants of vectorized Lauradoux/Walisch popcount.  (See
// popcount_vecs() in wdist_common.c for basic documentation.)
// Note that the size of the popcounted buffer is a hardcoded constant
// (specifically, (MULTIPLEX_DIST / BITCT) * 16 bytes).  The current code
// assumes (MULTIPLEX_DIST / BITCT) is a multiple of 3, and no greater than 30.
static inline unsigned int popcount_xor_1mask_multiword(__m128i** xor1p, __m128i* xor2, __m128i** maskp) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* xor2_end = &(xor2[MULTIPLEX_2DIST / 128]);

  acc.vi = _mm_setzero_si128();
  do {
    count1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), *((*maskp)++));
    count2 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), *((*maskp)++));
    half1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), *((*maskp)++));
    half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
    half1 = _mm_and_si128(half1, m1);
    count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
    count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
    count1 = _mm_add_epi64(count1, half1);
    count2 = _mm_add_epi64(count2, half2);
    count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
    count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
  } while (xor2 < xor2_end);
#if MULTIPLEX_DIST > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
  // are guaranteed to be <= 120, thus adding two together does not overflow
  // 255.
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
  acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
  return (unsigned int)(acc.u8[0] + acc.u8[1]);
}

static inline unsigned int popcount_xor_2mask_multiword(__m128i** xor1p, __m128i* xor2, __m128i** mask1p, __m128i* mask2) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* xor2_end = &(xor2[MULTIPLEX_2DIST / 128]);

  acc.vi = _mm_setzero_si128();
  do {
    count1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), _mm_and_si128(*((*mask1p)++), *mask2++));
    count2 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), _mm_and_si128(*((*mask1p)++), *mask2++));
    half1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), _mm_and_si128(*((*mask1p)++), *mask2++));
    half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
    half1 = _mm_and_si128(half1, m1);
    count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
    count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
    count1 = _mm_add_epi64(count1, half1);
    count2 = _mm_add_epi64(count2, half2);
    count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
    count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
  } while (xor2 < xor2_end);
#if MULTIPLEX_DIST > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
  acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
  return (unsigned int)(acc.u8[0] + acc.u8[1]);
}

static inline void ld_dot_prod(__m128i* vec1, __m128i* vec2, __m128i* mask1, __m128i* mask2, int* return_vals, int iters) {
  // Main routine for computation of \sum_i^M (x_i - \mu_x)(y_i - \mu_y), where
  // x_i, y_i \in \{-1, 0, 1\}, but there are missing values.
  //
  //
  // We decompose this sum into
  //   \sum_i x_iy_i - \mu_y\sum_i x_i - \mu_x\sum_i y_i +
  //   (M - # missing)\mu_x\mu_y.
  // *Without* missing values, this can be handled very cleanly.  The last
  // three terms can all be precomputed, and \sum_i x_iy_i can be handled in a
  // manner very similar to bitwise Hamming distance.  This is several times as
  // fast as the lookup tables used for relationship matrices.
  //
  // Unfortunately, when missing values are present, \mu_y\sum_i x_i and
  // \mu_x\sum_i y_i must be handled in the main loop, and this removes much of
  // the speed advantage.  So the best applications of the underlying ternary
  // dot product algorithm used here lie elsewhere.  Nevertheless, it is still
  // faster, so we use it.
  //
  //
  // Input:
  // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
  // * mask1 and mask2 mask out missing values (i.e. 00 for missing, 11 for
  //   nonmissing).
  // * return_vals provides space for return values.
  // * iters is the number of 48-byte windows to process, anywhere from 1 to 10
  //   inclusive.  (No, this is not the interface you'd use for a
  //   general-purpose library.)
  //
  // This function performs the update
  //   return_vals[0] += (-N) + \sum_i x_iy_i
  //   return_vals[1] += N_y + \sum_i x_i
  //   return_vals[2] += N_x + \sum_i y_i
  // where N is the number of individuals processed after applying the
  // missingness masks indicated by the subscripts.  The calculation currently
  // proceeds as follows:
  //
  // 1. N + \sum_i x_i = popcount_variant(vec1 & mask2)
  // The "variant" suffix refers to starting with two-bit integers instead of
  // one-bit integers in our summing process, so we get to skip a few
  // operations.  (Once all reserachers are using machines with fast hardware
  // popcount, a slightly different implementation would be better.)
  //
  // 2. zcheck := (vec1 | vec2) & 0x5555...
  // Detects whether at least one member of the pair has a 0/missing value.
  //
  // 3. popcount_variant(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
  // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i dot
  // product.
  //
  // MULTIPLEX_LD sets of values are handled per function call.  If fewer
  // values are present, it is currently safe to zero out the ends of all
  // input vectors.

  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  __m128i loader1;
  __m128i loader2;
  __m128i sum1;
  __m128i sum2;
  __m128i sum12;
  __m128i tmp_sum1;
  __m128i tmp_sum2;
  __m128i tmp_sum12;
  __uni16 acc;
  __uni16 acc1;
  __uni16 acc2;
  acc.vi = _mm_setzero_si128();
  acc1.vi = _mm_setzero_si128();
  acc2.vi = _mm_setzero_si128();
  do {
    loader1 = *vec1++;
    loader2 = *vec2++;
    sum1 = *mask2++;
    sum2 = *mask1++;
    sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    sum1 = _mm_and_si128(sum1, loader1);
    sum2 = _mm_and_si128(sum2, loader2);
    // use andnot to eliminate need for 0xaaaa... to occupy an xmm register
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12), _mm_xor_si128(loader1, loader2));
    sum12 = _mm_or_si128(sum12, loader1);

    // sum1, sum2, and sum12 now store the (biased) two-bit sums of
    // interest
    sum1 = _mm_add_epi64(_mm_and_si128(sum1, m2), _mm_and_si128(_mm_srli_epi64(sum1, 2), m2));
    sum2 = _mm_add_epi64(_mm_and_si128(sum2, m2), _mm_and_si128(_mm_srli_epi64(sum2, 2), m2));
    sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2), _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
    tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

    sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
    sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
    tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

    sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
    sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    acc1.vi = _mm_add_epi64(acc1.vi, _mm_add_epi64(_mm_and_si128(sum1, m4), _mm_and_si128(_mm_srli_epi64(sum1, 4), m4)));
    acc2.vi = _mm_add_epi64(acc2.vi, _mm_add_epi64(_mm_and_si128(sum2, m4), _mm_and_si128(_mm_srli_epi64(sum2, 4), m4)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sum12, m4), _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
  } while (--iters);
  // moved down since we're out of xmm registers
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
#if MULTIPLEX_LD > 960
  acc1.vi = _mm_add_epi64(_mm_and_si128(acc1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
  acc2.vi = _mm_add_epi64(_mm_and_si128(acc2.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc1.vi = _mm_and_si128(_mm_add_epi64(acc1.vi, _mm_srli_epi64(acc1.vi, 8)), m8);
  acc2.vi = _mm_and_si128(_mm_add_epi64(acc2.vi, _mm_srli_epi64(acc2.vi, 8)), m8);
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc1.vi = _mm_and_si128(_mm_add_epi64(acc1.vi, _mm_srli_epi64(acc1.vi, 16)), m16);
  acc2.vi = _mm_and_si128(_mm_add_epi64(acc2.vi, _mm_srli_epi64(acc2.vi, 16)), m16);
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
  acc1.vi = _mm_add_epi64(acc1.vi, _mm_srli_epi64(acc1.vi, 32));
  acc2.vi = _mm_add_epi64(acc2.vi, _mm_srli_epi64(acc2.vi, 32));
  acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
  return_vals[0] -= (unsigned int)(acc.u8[0] + acc.u8[1]);
  return_vals[1] += (unsigned int)(acc1.u8[0] + acc1.u8[1]);
  return_vals[2] += (unsigned int)(acc2.u8[0] + acc2.u8[1]);
}
#else
static inline unsigned int popcount_xor_1mask_multiword(unsigned long** xor1p, unsigned long* xor2, unsigned long** maskp) {
  unsigned long* xor2_end = &(xor2[MULTIPLEX_DIST / 16]);
  unsigned int bit_count = 0;
  unsigned long tmp_stor;
  unsigned long loader;
  unsigned long ulii;
  unsigned long uljj;
  do {
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    bit_count += (tmp_stor * 0x01010101) >> 24;
  } while (xor2 < xor2_end);
  return bit_count;
}

static inline unsigned int popcount_xor_2mask_multiword(unsigned long** xor1p, unsigned long* xor2, unsigned long** mask1p, unsigned long* mask2) {
  unsigned long* xor2_end = &(xor2[MULTIPLEX_DIST / 16]);
  unsigned int bit_count = 0;
  unsigned long tmp_stor;
  unsigned long loader;
  unsigned long ulii;
  unsigned long uljj;
  do {
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    loader -= ((loader >> 1) & FIVEMASK);
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    ulii += (loader & 0x33333333) + ((loader >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    loader -= ((loader >> 1) & FIVEMASK);
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    ulii += (loader & 0x33333333) + ((loader >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    bit_count += (tmp_stor * 0x01010101) >> 24;
  } while (xor2 < xor2_end);
  return bit_count;
}

static inline void ld_dot_prod(unsigned long* vec1, unsigned long* vec2, unsigned long* mask1, unsigned long* mask2, int* return_vals, int iters) {
  unsigned int final_sum1 = 0;
  unsigned int final_sum2 = 0;
  unsigned int final_sum12 = 0;
  unsigned long loader1;
  unsigned long loader2;
  unsigned long sum1;
  unsigned long sum2;
  unsigned long sum12;
  unsigned long tmp_sum1;
  unsigned long tmp_sum2;
  unsigned long tmp_sum12;
  do {
    // (The important part of the header comment on the 64-bit version is
    // copied below.)
    //
    // Input:
    // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
    // * mask1 and mask2 mask out missing values (i.e. 00 for missing, 11 for
    //   nonmissing).
    // * return_vals provides space for return values.
    // * iters is the number of 12-byte windows to process, anywhere from 1 to
    //   40 inclusive.  (No, this is not the interface you'd use for a
    //   general-purpose library.)  [32- and 64-bit differ here.]
    //
    // This function performs the update
    //   return_vals[0] += (-N) + \sum_i x_iy_i
    //   return_vals[1] += N_y + \sum_i x_i
    //   return_vals[2] += N_x + \sum_i y_i
    // where N is the number of individuals processed after applying the
    // missingness masks indicated by the subscripts.  The calculation
    // currently proceeds as follows:
    //
    // 1. N + \sum_i x_i = popcount_variant(vec1 & mask2)
    // The "variant" suffix refers to starting with two-bit integers instead of
    // one-bit integers in our summing process, so we get to skip a few
    // operations.  (Once all reserachers are using machines with fast hardware
    // popcount, a slightly different implementation would be better.)
    //
    // 2. zcheck := (vec1 | vec2) & 0x5555...
    // Detects whether at least one member of the pair has a 0/missing value.
    //
    // 3. popcount_variant(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
    // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i
    // dot product.

    loader1 = *vec1++;
    loader2 = *vec2++;
    sum1 = *mask2++;
    sum2 = *mask1++;
    sum12 = (loader1 | loader2) & FIVEMASK;

    sum1 = sum1 & loader1;
    sum2 = sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
    sum12 = sum12 | loader1;

    sum1 = (sum1 & 0x33333333) + ((sum1 >> 2) & 0x33333333);
    sum2 = (sum2 & 0x33333333) + ((sum2 >> 2) & 0x33333333);
    sum12 = (sum12 & 0x33333333) + ((sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;

    tmp_sum1 = tmp_sum1 & loader1;
    tmp_sum2 = tmp_sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;

    sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
    sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;

    tmp_sum1 = tmp_sum1 & loader1;
    tmp_sum2 = tmp_sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;

    sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
    sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    sum1 = (sum1 & 0x0f0f0f0f) + ((sum1 >> 4) & 0x0f0f0f0f);
    sum2 = (sum2 & 0x0f0f0f0f) + ((sum2 >> 4) & 0x0f0f0f0f);
    sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);

    // technically could do the multiply-and-shift only once every two rounds
    final_sum1 += (sum1 * 0x01010101) >> 24;
    final_sum2 += (sum2 * 0x01010101) >> 24;
    final_sum12 += (sum12 * 0x01010101) >> 24;
  } while (--iters);
  return_vals[0] -= final_sum12;
  return_vals[1] += final_sum1;
  return_vals[2] += final_sum2;
}
#endif

void incr_dists_i(int* idists, unsigned long* geno, int tidx) {
#if __LP64__
  __m128i* glptr;
  __m128i* glptr2;
  __m128i* mptr;
  __m128i* mcptr_start;
  unsigned long* lptr;
#else
  unsigned long* glptr;
  unsigned long* glptr2;
  unsigned long* mptr;
  unsigned long* mcptr_start;
#endif
  unsigned int uii;
  int jj;
  unsigned long mask_fixed;
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    jj = uii * (MULTIPLEX_2DIST / BITCT);
#if __LP64__
    glptr = (__m128i*)geno;
    glptr2 = (__m128i*)(&(geno[jj]));
    lptr = &(g_masks[jj]);
    mcptr_start = (__m128i*)lptr;
    mask_fixed = *lptr++;
    for (jj = 0; jj < MULTIPLEX_2DIST / BITCT - 1; jj++) {
      mask_fixed &= *lptr++;
    }
    mptr = (__m128i*)g_masks;
#else
    glptr = geno;
    glptr2 = &(geno[jj]);
    mcptr_start = &(g_masks[jj]);
    mptr = mcptr_start;
    mask_fixed = *mptr++;
    for (jj = 0; jj < MULTIPLEX_2DIST / BITCT - 1; jj++) {
      mask_fixed &= *mptr++;
    }
    mptr = g_masks;
#endif
    if (~mask_fixed) {
      while (glptr < glptr2) {
	*idists += popcount_xor_2mask_multiword(&glptr, glptr2, &mptr, mcptr_start);
	idists++;
      }
    } else {
      while (glptr < glptr2) {
	*idists += popcount_xor_1mask_multiword(&glptr, glptr2, &mptr);
	idists++;
      }
    }
  }
}

void* calc_idist_thread(void* arg) {
  long tidx = (long)arg;
  int ii = g_thread_start[tidx];
  int jj = g_thread_start[0];
  incr_dists_i(&(g_idists[((long long)ii * (ii - 1) - (long long)jj * (jj - 1)) / 2]), (unsigned long*)g_geno, (int)tidx);
  return NULL;
}

void incr_genome(unsigned int* genome_main, unsigned long* geno, int tidx) {
#if __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i xor_buf[GENOME_MULTIPLEX / BITCT];
  __m128i* xor_buf_end = &(xor_buf[GENOME_MULTIPLEX / BITCT]);
  __m128i* maskptr;
  __m128i* maskptr_fixed;
  __m128i* maskptr_fixed_tmp;
  __m128i* xor_ptr;
  __m128i* glptr_fixed_tmp;
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i count_ibs1;
  __m128i count_ibs0;
  __m128i count2_ibs1;
  __m128i count2_ibs0;
  __uni16 acc_ibs1;
  __uni16 acc_ibs0;
  unsigned long* lptr;
  __m128i* glptr;
  __m128i* glptr_fixed;
  __m128i* glptr_end;
#else
  unsigned long* glptr;
  unsigned long* glptr_fixed;
  unsigned long* glptr_end;
  unsigned long* maskptr;
  unsigned long* maskptr_fixed;
  unsigned long* maskptr_fixed_tmp;
  unsigned long* glptr_fixed_tmp;
  unsigned long xor_buf[GENOME_MULTIPLEX2 / BITCT];
  unsigned long* xor_buf_end = &(xor_buf[GENOME_MULTIPLEX2 / BITCT]);
  unsigned long* xor_ptr;
  unsigned int bit_count_ibs1 = 0;
  unsigned int bit_count_ibs0 = 0;
  unsigned long bitfield_ibs1;
  unsigned long bitfield_ibs0;
  unsigned long loader;
  unsigned long loader2;
  unsigned long tmp_stor_ibs1;
  unsigned long tmp_stor_ibs0;
#endif
  unsigned long* glptr_back;
  unsigned long ibs_incr;
  unsigned int uii;
  int jj;
  int offset;
  unsigned long uland;
  unsigned long ulval;
  unsigned long next_ppc_marker_hybrid;
  unsigned long mask_fixed_test;
  unsigned long* marker_window_ptr;
  int lowct2 = g_low_ct * 2;
  int highct2 = g_high_ct * 2;
#if __LP64__
  glptr_end = (__m128i*)(&(geno[g_indiv_ct * (GENOME_MULTIPLEX2 / BITCT)]));
#else
  glptr_end = &(geno[g_indiv_ct * (GENOME_MULTIPLEX2 / BITCT)]);
#endif
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    jj = uii * (GENOME_MULTIPLEX2 / BITCT);
#if __LP64__
    glptr_fixed = (__m128i*)(&(geno[jj]));
    glptr = (__m128i*)(&(geno[jj + (GENOME_MULTIPLEX2 / BITCT)]));
    lptr = &(g_masks[jj]);
    maskptr = (__m128i*)(&(g_masks[jj + (GENOME_MULTIPLEX2 / BITCT)]));
    maskptr_fixed = (__m128i*)lptr;
    mask_fixed_test = *lptr++;
    for (jj = 0; jj < GENOME_MULTIPLEX2 / BITCT - 1; jj++) {
      mask_fixed_test &= *lptr++;
    }
#else
    glptr_fixed = &(geno[jj]);
    glptr = &(geno[jj + (GENOME_MULTIPLEX2 / BITCT)]);
    maskptr_fixed = &(g_masks[jj]);
    maskptr = maskptr_fixed;
    mask_fixed_test = *maskptr++;
    for (jj = 0; jj < GENOME_MULTIPLEX2 / BITCT - 1; jj++) {
      mask_fixed_test &= *maskptr++;
    }
#endif
    if (~mask_fixed_test) {
      while (glptr < glptr_end) {
	xor_ptr = xor_buf;
	glptr_back = (unsigned long*)glptr;
	glptr_fixed_tmp = glptr_fixed;
	maskptr_fixed_tmp = maskptr_fixed;
#if __LP64__
	acc_ibs1.vi = _mm_setzero_si128();
	acc_ibs0.vi = _mm_setzero_si128();
	do {
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
	  count_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count_ibs0;

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
	  count2_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count2_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count2_ibs0;

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);

          count_ibs1 = _mm_add_epi64(_mm_and_si128(count_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count_ibs1, 2), m2));
          count_ibs0 = _mm_add_epi64(_mm_and_si128(count_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count_ibs0, 2), m2));
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_add_epi64(_mm_and_si128(count2_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs1, 2), m2)));
          count_ibs0 = _mm_add_epi64(count_ibs0, _mm_add_epi64(_mm_and_si128(count2_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs0, 2), m2)));
          acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_add_epi64(_mm_and_si128(count_ibs1, m4), _mm_and_si128(_mm_srli_epi64(count_ibs1, 4), m4)));
          acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_add_epi64(_mm_and_si128(count_ibs0, m4), _mm_and_si128(_mm_srli_epi64(count_ibs0, 4), m4)));
	} while (xor_ptr < xor_buf_end);
#if GENOME_MULTIPLEX > 1920
	acc_ibs1.vi = _mm_add_epi64(_mm_and_si128(acc_ibs1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs1.vi, 8), m8));
	acc_ibs0.vi = _mm_add_epi64(_mm_and_si128(acc_ibs0.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs0.vi, 8), m8));
#else
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 8)), m8);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 8)), m8);
#endif
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 16)), m16);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 16)), m16);
	acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 32));
	acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 32));
	*genome_main += (unsigned int)(acc_ibs1.u8[0] + acc_ibs1.u8[1]);
	genome_main++;
        *genome_main += (unsigned int)(acc_ibs0.u8[0] + acc_ibs0.u8[1]);
#else
        bit_count_ibs1 = 0;
	bit_count_ibs0 = 0;
	do {
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & ((*maskptr++) & (*maskptr_fixed_tmp++));
	  loader2 = loader >> 1;
	  bitfield_ibs1 = (loader ^ loader2) & FIVEMASK;
	  bitfield_ibs0 = (loader & loader2) & FIVEMASK;
	  *xor_ptr++ = bitfield_ibs0;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & ((*maskptr++) & (*maskptr_fixed_tmp++));
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & ((*maskptr++) & (*maskptr_fixed_tmp++));
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
          bitfield_ibs1 = (bitfield_ibs1 & 0x33333333) + ((bitfield_ibs1 >> 2) & 0x33333333);
	  bitfield_ibs0 = (bitfield_ibs0 & 0x33333333) + ((bitfield_ibs0 >> 2) & 0x33333333);
	  tmp_stor_ibs1 = (bitfield_ibs1 + (bitfield_ibs1 >> 4)) & 0x0f0f0f0f;
	  tmp_stor_ibs0 = (bitfield_ibs0 + (bitfield_ibs0 >> 4)) & 0x0f0f0f0f;

          bit_count_ibs1 += (tmp_stor_ibs1 * 0x01010101) >> 24;
	  bit_count_ibs0 += (tmp_stor_ibs0 * 0x01010101) >> 24;
	} while (xor_ptr < xor_buf_end);
	*genome_main += bit_count_ibs1;
	genome_main++;
	*genome_main += bit_count_ibs0;
#endif
	genome_main++;
	next_ppc_marker_hybrid = *genome_main - lowct2;
	if (next_ppc_marker_hybrid < GENOME_MULTIPLEX2) {
	  ibs_incr = 0; // hethet low-order, ibs0 high-order

          // This PPC test is now the limiting step of --genome, not the IBS
	  // matrix.
	  //
          // I've taken a few "maintenance nightmare" liberties with this code
          // to speed it up, such as using a single lookup table that stores
          // values in two different forms (distinguished by the high bit of
          // the value), and using gotos for non-error-conditions, since the
          // loop is quite short.  In the long run, we want to implement
	  // support for better IBD estimation; see e.g. Browning B.L.,
	  // Browning S.R. "A Fast, Powerful Method for Detecting Identity by
	  // Descent", which discusses some weaknesses of PLINK --genome.
	  do {
	    offset = next_ppc_marker_hybrid / BITCT;
	    marker_window_ptr = &(g_marker_window[offset * BITCT]);
	    next_ppc_marker_hybrid = ~0LU << (next_ppc_marker_hybrid & (BITCT - 1));
	  incr_genome_2mask_loop:
	    uland = glptr_back[offset] & (((unsigned long*)glptr_fixed)[offset]);
	    // het is represented as 11, so
	    //   (uland & (uland << 1)) & 0xaaaaaaaaaaaaaaaa
	    // stores whether a particular marker is a hethet hit in the
	    // corresponding odd bit.
	    //
	    // homozygotes are represented as 01 and 10, so
	    //   (ulxor & (ulxor >> 1)) & 0x5555555555555555
	    // stores whether a particular marker is a hom1-hom2 hit in the
	    // corresponding even bit.  (het-missing pairs also set that bit,
	    // but the masking filters that out.)
	    //
	    // ~0LU << xx masks out the bottom xx bits.
	    ulval = (((uland & (uland << 1)) & AAAAMASK) | (((unsigned long*)xor_buf)[offset]));
	    do {
	      ulval &= next_ppc_marker_hybrid;
	      if (ulval) {
		jj = __builtin_ctzl(ulval);
		next_ppc_marker_hybrid = marker_window_ptr[jj];
		ibs_incr += (1LU << ((jj & 1) * BITCT2));
	      } else if (offset < ((GENOME_MULTIPLEX2 - BITCT) / BITCT)) {
		offset++;
		next_ppc_marker_hybrid = ~0LU;
		marker_window_ptr = &(marker_window_ptr[BITCT]);
		goto incr_genome_2mask_loop;
	      } else {
		*genome_main = highct2;
		goto incr_genome_2mask_exit;
	      }
	    } while (next_ppc_marker_hybrid & (1LU << (BITCT - 1)));
	  } while (next_ppc_marker_hybrid < GENOME_MULTIPLEX2);
	  *genome_main = next_ppc_marker_hybrid + lowct2;
        incr_genome_2mask_exit:
	  genome_main++;
          *genome_main += ibs_incr & ((~0LU) >> BITCT2);
	  genome_main++;
	  *genome_main += ibs_incr >> BITCT2;
	  genome_main++;
	} else {
	  genome_main = &(genome_main[3]);
	}
      }
    } else {
      while (glptr < glptr_end) {
	xor_ptr = xor_buf;
	glptr_back = (unsigned long*)glptr;
	glptr_fixed_tmp = glptr_fixed;
#if __LP64__
	acc_ibs1.vi = _mm_setzero_si128();
	acc_ibs0.vi = _mm_setzero_si128();
	do {
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
	  count_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count_ibs0;
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
	  count2_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count2_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count2_ibs0;
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);

          count_ibs1 = _mm_add_epi64(_mm_and_si128(count_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count_ibs1, 2), m2));
          count_ibs0 = _mm_add_epi64(_mm_and_si128(count_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count_ibs0, 2), m2));
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_add_epi64(_mm_and_si128(count2_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs1, 2), m2)));
          count_ibs0 = _mm_add_epi64(count_ibs0, _mm_add_epi64(_mm_and_si128(count2_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs0, 2), m2)));
          acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_add_epi64(_mm_and_si128(count_ibs1, m4), _mm_and_si128(_mm_srli_epi64(count_ibs1, 4), m4)));
          acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_add_epi64(_mm_and_si128(count_ibs0, m4), _mm_and_si128(_mm_srli_epi64(count_ibs0, 4), m4)));
	} while (xor_ptr < xor_buf_end);
#if GENOME_MULTIPLEX > 1920
	acc_ibs1.vi = _mm_add_epi64(_mm_and_si128(acc_ibs1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs1.vi, 8), m8));
	acc_ibs0.vi = _mm_add_epi64(_mm_and_si128(acc_ibs0.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs0.vi, 8), m8));
#else
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 8)), m8);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 8)), m8);
#endif
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 16)), m16);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 16)), m16);
	acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 32));
	acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 32));
	*genome_main += (unsigned int)(acc_ibs1.u8[0] + acc_ibs1.u8[1]);
	genome_main++;
        *genome_main += (unsigned int)(acc_ibs0.u8[0] + acc_ibs0.u8[1]);
#else
        bit_count_ibs1 = 0;
	bit_count_ibs0 = 0;
	do {
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & (*maskptr++);
	  loader2 = loader >> 1;
	  bitfield_ibs1 = (loader ^ loader2) & FIVEMASK;
	  bitfield_ibs0 = (loader & loader2) & FIVEMASK;
	  *xor_ptr++ = bitfield_ibs0;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & (*maskptr++);
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & (*maskptr++);
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
          bitfield_ibs1 = (bitfield_ibs1 & 0x33333333) + ((bitfield_ibs1 >> 2) & 0x33333333);
	  bitfield_ibs0 = (bitfield_ibs0 & 0x33333333) + ((bitfield_ibs0 >> 2) & 0x33333333);
	  tmp_stor_ibs1 = (bitfield_ibs1 + (bitfield_ibs1 >> 4)) & 0x0f0f0f0f;
	  tmp_stor_ibs0 = (bitfield_ibs0 + (bitfield_ibs0 >> 4)) & 0x0f0f0f0f;

          bit_count_ibs1 += (tmp_stor_ibs1 * 0x01010101) >> 24;
	  bit_count_ibs0 += (tmp_stor_ibs0 * 0x01010101) >> 24;
	} while (xor_ptr < xor_buf_end);
	*genome_main += bit_count_ibs1;
	genome_main++;
	*genome_main += bit_count_ibs0;
#endif
	genome_main++;
	next_ppc_marker_hybrid = *genome_main - lowct2;
	if (next_ppc_marker_hybrid < GENOME_MULTIPLEX2) {
	  ibs_incr = 0; // hethet low-order, ibs0 high-order
	  do {
	    offset = next_ppc_marker_hybrid / BITCT;
	    marker_window_ptr = &(g_marker_window[offset * BITCT]);
	    next_ppc_marker_hybrid = ~0LU << (next_ppc_marker_hybrid & (BITCT - 1));
	  incr_genome_1mask_loop:
	    uland = glptr_back[offset] & (((unsigned long*)glptr_fixed)[offset]);
	    ulval = ((uland & (uland << 1)) & AAAAMASK) | (((unsigned long*)xor_buf)[offset]);
	    do {
	      ulval &= next_ppc_marker_hybrid;
	      if (ulval) {
		jj = __builtin_ctzl(ulval);
		next_ppc_marker_hybrid = marker_window_ptr[jj];
		ibs_incr += (1LU << ((jj & 1) * BITCT2));
	      } else if (offset < ((GENOME_MULTIPLEX2 - BITCT) / BITCT)) {
		offset++;
		next_ppc_marker_hybrid = ~0LU;
		marker_window_ptr = &(marker_window_ptr[BITCT]);
		goto incr_genome_1mask_loop;
	      } else {
		*genome_main = highct2;
		goto incr_genome_1mask_exit;
	      }
	    } while (next_ppc_marker_hybrid & (1LU << (BITCT - 1)));
	  } while (next_ppc_marker_hybrid < GENOME_MULTIPLEX2);
	  *genome_main = next_ppc_marker_hybrid + lowct2;
	incr_genome_1mask_exit:
	  genome_main++;
	  *genome_main += ibs_incr & ((~0LU) >> BITCT2);
	  genome_main++;
	  *genome_main += ibs_incr >> BITCT2;
	  genome_main++;
	} else {
	  genome_main = &(genome_main[3]);
	}
      }
    }
  }
}

void* calc_genome_thread(void* arg) {
  long tidx = (long)arg;
  int ii = g_thread_start[tidx];
  int jj = g_thread_start[0];
  incr_genome(&(g_genome_main[((long long)g_indiv_ct * (ii - jj) - ((long long)ii * (ii + 1) - (long long)jj * (jj + 1)) / 2) * 5]), (unsigned long*)g_geno, (int)tidx);
  return NULL;
}

void incr_dists(double* dists, unsigned long* geno, int tidx) {
  unsigned long* glptr;
  unsigned long ulii;
  unsigned long mask_fixed;
  unsigned long uljj;
  unsigned long* mptr;
  double* weights1 = &(g_weights[16384]);
#if __LP64__
  double* weights2 = &(g_weights[32768]);
  double* weights3 = &(g_weights[36864]);
  double* weights4 = &(g_weights[40960]);
#endif
  unsigned int uii;
  unsigned int ujj;
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    glptr = geno;
    ulii = geno[uii];
    mptr = g_masks;
    mask_fixed = g_masks[uii];
#if __LP64__
    if (mask_fixed == ~0LU) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (*mptr++);
        *dists += weights4[uljj >> 52] + weights3[(uljj >> 40) & 4095] + weights2[(uljj >> 28) & 4095] + weights1[(uljj >> 14) & 16383] + g_weights[uljj & 16383];
	dists++;
      }
    } else {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (mask_fixed & (*mptr++));
        *dists += weights4[uljj >> 52] + weights3[(uljj >> 40) & 4095] + weights2[(uljj >> 28) & 4095] + weights1[(uljj >> 14) & 16383] + g_weights[uljj & 16383];
	dists++;
      }
    }
#else
    if (mask_fixed == 0x0fffffff) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (*mptr++);
	*dists += weights1[uljj >> 14] + g_weights[uljj & 16383];
	dists++;
      }
    } else {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (mask_fixed & (*mptr++));
	*dists += weights1[uljj >> 14] + g_weights[uljj & 16383];
	dists++;
      }
    }
#endif
  }
}

void* calc_dist_thread(void* arg) {
  long tidx = (long)arg;
  int ii = g_thread_start[tidx];
  int jj = g_thread_start[0];
  incr_dists(&(g_dists[((long long)ii * (ii - 1) - (long long)jj * (jj - 1)) / 2]), (unsigned long*)g_geno, (int)tidx);
  return NULL;
}

void incr_wt_dist_missing(unsigned int* mtw, int tidx) {
  unsigned long* glptr;
  unsigned long ulii;
  unsigned long uljj;
  unsigned int uii;
  unsigned int ujj;
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    glptr = g_mmasks;
    ulii = g_mmasks[uii];
    if (ulii) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++) & ulii;
        while (uljj) {
          mtw[ujj] += g_weights_i[__builtin_ctzl(uljj)];
          uljj &= uljj - 1;
        }
      }
    }
    mtw = &(mtw[uii]);
  }
}

void* calc_distm_thread(void* arg) {
  long tidx = (long)arg;
  int ii = g_thread_start[tidx];
  int jj = g_thread_start[0];
  incr_wt_dist_missing(&(g_missing_tot_weights[((long long)ii * (ii - 1) - (long long)jj * (jj - 1)) / 2]), (int)tidx);
  return NULL;
}

void incr_dists_r(double* dists, unsigned long* geno, unsigned long* masks, int tidx, double* weights) {
  unsigned long* glptr;
  unsigned long* maskptr;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long basemask;
  double* weights1 = &(weights[32768]);
#if __LP64__
  double* weights2 = &(weights[65536]);
  double* weights3 = &(weights[98304]);
#endif
  unsigned int uii;
  unsigned int ujj;
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    glptr = geno;
    ulii = geno[uii];
    maskptr = masks;
    basemask = masks[uii];
    if (!basemask) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = ((*glptr++) + ulii) | (*maskptr++);
#if __LP64__
	*dists += weights3[uljj >> 45] + weights2[(uljj >> 30) & 32767] + weights1[(uljj >> 15) & 32767] + weights[uljj & 32767];
#else
	*dists += weights1[uljj >> 15] + weights[uljj & 32767];
#endif
	dists++;
      }
    } else {
      for (ujj = 0; ujj < uii; ujj++) {
        uljj = ((*glptr++) + ulii) | ((*maskptr++) | basemask);
#if __LP64__
	*dists += weights3[uljj >> 45] + weights2[(uljj >> 30) & 32767] + weights1[(uljj >> 15) & 32767] + weights[uljj & 32767];
#else
	*dists += weights1[uljj >> 15] + weights[uljj & 32767];
#endif
	dists++;
      }
    }
  }
}

void* calc_rel_thread(void* arg) {
  long tidx = (long)arg;
  int ii = g_thread_start[tidx];
  int jj = g_thread_start[0];
  incr_dists_r(&(g_rel_dists[((long long)ii * (ii - 1) - (long long)jj * (jj - 1)) / 2]), (unsigned long*)g_geno, g_masks, (int)tidx, g_weights);
  return NULL;
}

void incr_dists_rm(unsigned int* idists, int tidx, unsigned int* thread_start) {
  // count missing intersection, optimized for sparsity
  unsigned long* mlptr;
  unsigned long ulii;
  unsigned long uljj;
  unsigned int uii;
  unsigned int ujj;
  for (uii = thread_start[tidx]; uii < thread_start[tidx + 1]; uii++) {
    mlptr = g_mmasks;
    ulii = g_mmasks[uii];
    if (ulii) {
      for (ujj = 0; ujj < uii; ujj++) {
        uljj = (*mlptr++) & ulii;
        while (uljj) {
          idists[ujj] += 1;
          uljj &= uljj - 1;
        }
      }
    }
    idists = &(idists[uii]);
  }
}

void* calc_missing_thread(void* arg) {
  long tidx = (long)arg;
  int ii = g_thread_start[tidx];
  int jj = g_thread_start[0];
  incr_dists_rm(&(g_missing_dbl_excluded[((long long)ii * (ii - 1) - (long long)jj * (jj - 1)) / 2]), (int)tidx, g_thread_start);
  return NULL;
}

inline int flexclose_null(FILE** outfile_ptr, gzFile* gz_outfile_ptr) {
  int ii;
  if (*outfile_ptr) {
    return fclose_null(outfile_ptr);
  } else {
    ii = gzclose(*gz_outfile_ptr);
    *gz_outfile_ptr = NULL;
    return (ii != Z_OK);
  }
}

void incr_dists_rm_inv(unsigned int* idists, int tidx) {
  // inverted loops for --genome --parallel
  unsigned long* glptr;
  unsigned long ulii;
  unsigned long uljj;
  unsigned int indiv_ct_m1 = g_indiv_ct - 1;
  unsigned int uii;
  unsigned int ujj;
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    ulii = g_mmasks[uii];
    if (ulii) {
      glptr = &(g_mmasks[uii + 1]);
      // ujj is deliberately biased down by 1
      for (ujj = uii; ujj < indiv_ct_m1; ujj++) {
        uljj = (*glptr++) & ulii;
        while (uljj) {
          idists[ujj] += 1;
          uljj &= uljj - 1;
        }
      }
    }
    idists = &(idists[indiv_ct_m1 - uii - 1]);
  }
}

void* calc_genomem_thread(void* arg) {
  long tidx = (long)arg;
  int ii = g_thread_start[tidx];
  int jj = g_thread_start[0];
  // f(0) = 0
  // f(1) = ic - 2
  // f(2) = 2ic - 5
  // f(3) = 3ic - 9
  // ...
  // f(n) = nic - (n+1)(n+2)/2 + 1
  incr_dists_rm_inv(&(g_missing_dbl_excluded[(long long)g_indiv_ct * (ii - jj) - ((long long)(ii + 1) * (ii + 2) - (long long)(jj + 1) * (jj + 2)) / 2]), (int)tidx);
  return NULL;
}

void groupdist_jack(int* ibuf, double* returns) {
  int* iptr = ibuf;
  int* jptr;
  unsigned int indiv_idx;
  double* dptr = g_dists;
  char* cptr;
  char* cptr2;
  char cc;
  char cc2;
  double neg_tot_aa = 0.0;
  double neg_tot_au = 0.0;
  double neg_tot_uu = 0.0;
  int neg_a = 0;
  int neg_u = 0;
  double dxx;
  for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
    cc = g_pheno_c[indiv_idx];
    if (cc == -1) {
      continue;
    }
    dptr = &(g_dists[((unsigned long)indiv_idx * (indiv_idx - 1)) / 2]);
    if (indiv_idx == (unsigned int)(*iptr)) {
      cptr = g_pheno_c;
      cptr2 = &(g_pheno_c[indiv_idx]);
      if (cc) {
        neg_a++;
        while (cptr < cptr2) {
          cc2 = *cptr++;
          dxx = *dptr++;
          if (cc2 == 1) {
            neg_tot_aa += dxx;
          } else if (!cc2) {
            neg_tot_au += dxx;
          }
        }
      } else {
        neg_u++;
        while (cptr < cptr2) {
          cc2 = *cptr++;
          dxx = *dptr++;
          if (cc2 == 1) {
            neg_tot_au += dxx;
          } else if (!cc2) {
            neg_tot_uu += dxx;
          }
        }
      }
      iptr++;
    } else {
      jptr = ibuf;
      while (jptr < iptr) {
        cc2 = g_pheno_c[*jptr];
        dxx = dptr[*jptr];
	if ((cc == 1) && (cc2 == 1)) {
	  neg_tot_aa += dxx;
	} else if ((cc == 0) && (cc2 == 0)) {
	  neg_tot_uu += dxx;
	} else {
	  neg_tot_au += dxx;
	}
        jptr++;
      }
    }
  }
  returns[0] = (g_reg_tot_x - neg_tot_aa) / (double)(((long)(g_high_ct - neg_a) * (g_high_ct - neg_a - 1)) / 2);
  returns[1] = (g_reg_tot_xy - neg_tot_au) / (double)((long)(g_high_ct - neg_a) * (g_low_ct - neg_u));
  returns[2] = (g_reg_tot_y - neg_tot_uu) / (double)(((long)(g_low_ct - neg_u) * (g_low_ct - neg_u - 1)) / 2);
}

void small_remap(int* ibuf, unsigned int ct, unsigned int dd) {
  int* ibuf_end = &(ibuf[dd]);
  int missings = 0;
  int curpos = 0;
  do {
    if (g_pheno_c[curpos] == -1) {
      missings++;
    } else if (*ibuf == curpos - missings) {
      *ibuf++ = curpos;
    }
    curpos++;
  } while (ibuf < ibuf_end);
}

void* groupdist_jack_thread(void* arg) {
  long tidx = (long)arg;
  int* ibuf = (int*)(&(g_geno[tidx * CACHEALIGN(g_high_ct + g_low_ct + (g_jackknife_d + 1) * sizeof(int))]));
  unsigned char* cbuf = &(g_geno[tidx * CACHEALIGN(g_high_ct + g_low_ct + (g_jackknife_d + 1) * sizeof(int)) + (g_jackknife_d + 1) * sizeof(int)]);
  unsigned long long ulii;
  unsigned long long uljj = g_jackknife_iters / 100;
  double returns[3];
  double results[9];
  double new_old_diff[3];
  fill_double_zero(results, 9);
  for (ulii = 0; ulii < g_jackknife_iters; ulii++) {
    pick_d_small(cbuf, ibuf, g_high_ct + g_low_ct, g_jackknife_d);
    if (g_high_ct + g_low_ct < g_indiv_ct) {
      small_remap(ibuf, g_high_ct + g_low_ct, g_jackknife_d);
    }
    groupdist_jack(ibuf, returns);
    if (ulii > 0) {
      new_old_diff[0] = returns[0] - results[0];
      new_old_diff[1] = returns[1] - results[1];
      new_old_diff[2] = returns[2] - results[2];
      results[0] += new_old_diff[0] / (ulii + 1); // AA mean
      results[1] += new_old_diff[1] / (ulii + 1); // AU mean
      results[2] += new_old_diff[2] / (ulii + 1); // UU mean
      results[3] += (returns[0] - results[0]) * new_old_diff[0]; // AA var
      results[4] += (returns[1] - results[1]) * new_old_diff[1]; // AU var
      results[5] += (returns[2] - results[2]) * new_old_diff[2]; // UU var
      results[6] += (returns[0] - results[0]) * new_old_diff[1]; // AA-AU cov
      results[7] += (returns[0] - results[0]) * new_old_diff[2]; // AA-UU cov
      results[8] += (returns[1] - results[1]) * new_old_diff[2]; // AU-UU cov
    } else {
      results[0] += returns[0];
      results[1] += returns[1];
      results[2] += returns[2];
    }
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100) / g_jackknife_iters;
      printf("\r%lld%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1) * g_jackknife_iters) / 100;
    }
  }
  // don't write until end, to avoid false sharing
  for (ulii = 0; ulii < 9; ulii++) {
    g_calc_result[ulii][tidx] = results[ulii];
  }
  return NULL;
}

double regress_rel_jack(int* ibuf, double* ret2_ptr) {
  int* iptr = ibuf;
  int* jptr = &(ibuf[g_jackknife_d]);
  unsigned int uii;
  int jj;
  int kk;
  double* dptr;
  double* dptr2;
  double neg_tot_xy = 0.0;
  double neg_tot_x = 0.0;
  double neg_tot_y = 0.0;
  double neg_tot_xx = 0.0;
  double neg_tot_yy = 0.0;
  double dxx;
  double dxx1;
  double dyy;
  while (iptr < jptr) {
    dptr2 = &(g_jackknife_precomp[(*iptr++) * JACKKNIFE_VALS_REL]);
    neg_tot_xy += *dptr2++;
    neg_tot_x += *dptr2++;
    neg_tot_y += *dptr2++;
    neg_tot_xx += *dptr2++;
    neg_tot_yy += *dptr2++;
  }
  iptr = ibuf;
  for (uii = 1; uii < g_jackknife_d; uii++) {
    jj = *(++iptr);
    dxx1 = g_pheno_packed[jj];
    jptr = ibuf;
    dptr = &(g_rel_dists[((long)jj * (jj - 1)) / 2]);
    while (jptr < iptr) {
      kk = *jptr++;
      dxx = (dxx1 + g_pheno_packed[kk]) * 0.5;
      dyy = dptr[kk];
      neg_tot_xy -= dxx * dyy;
      neg_tot_x -= dxx;
      neg_tot_y -= dyy;
      neg_tot_xx -= dxx * dxx;
      neg_tot_yy -= dyy * dyy;
    }
  }
  dxx = g_reg_tot_y - neg_tot_y;
  dyy = g_indiv_ct - g_jackknife_d;
  dyy = dyy * (dyy - 1.0) * 0.5;
  *ret2_ptr = ((g_reg_tot_xy - neg_tot_xy) - dxx * (g_reg_tot_x - neg_tot_x) / dyy) / ((g_reg_tot_yy - neg_tot_yy) - dxx * dxx / dyy);
  dxx = g_reg_tot_x - neg_tot_x;
  return ((g_reg_tot_xy - neg_tot_xy) - dxx * (g_reg_tot_y - neg_tot_y) / dyy) / ((g_reg_tot_xx - neg_tot_xx) - dxx * dxx / dyy);
}

void* regress_rel_jack_thread(void* arg) {
  long tidx = (long)arg;
  int* ibuf = (int*)(&(g_geno[tidx * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int))]));
  unsigned char* cbuf = &(g_geno[tidx * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int)) + (g_jackknife_d + 1) * sizeof(int)]);
  unsigned long long ulii;
  unsigned long long uljj = g_jackknife_iters / 100;
  double sum = 0.0;
  double sum_sq = 0.0;
  double sum2 = 0;
  double sum2_sq = 0.0;
  double dxx;
  double ret2;
  for (ulii = 0; ulii < g_jackknife_iters; ulii++) {
    pick_d_small(cbuf, ibuf, g_indiv_ct, g_jackknife_d);
    dxx = regress_rel_jack(ibuf, &ret2);
    sum += dxx;
    sum_sq += dxx * dxx;
    sum2 += ret2;
    sum2_sq += ret2 * ret2;
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100) / g_jackknife_iters;
      printf("\r%lld%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1) * g_jackknife_iters) / 100;
    }
  }
  g_calc_result[0][tidx] = sum;
  g_calc_result[1][tidx] = sum_sq;
  g_calc_result[2][tidx] = sum2;
  g_calc_result[3][tidx] = sum2_sq;
  return NULL;
}

int regress_rel_main(unsigned long* indiv_exclude, unsigned int indiv_ct, unsigned long regress_rel_iters, int regress_rel_d, pthread_t* threads) {
  double* rel_ptr;
  double* pheno_ptr;
  double* pheno_ptr2;
  double* jp_fixed_ptr;
  double* jp_moving_ptr;
  unsigned int uii;
  unsigned long ulii;
  unsigned long trimatrix_size;
  double trimatrix_size_recip;
  double half_avg_pheno;
  double dxx;
  double dyy;
  double dxxyy;
  double dxxsq;
  double dyysq;
  g_pheno_packed = (double*)wkspace_alloc(indiv_ct * sizeof(double));
  if (!g_pheno_packed) {
    return RET_NOMEM;
  }
  collapse_copy_phenod(g_pheno_packed, g_pheno_d, indiv_exclude, indiv_ct);
  print_pheno_stdev(g_pheno_packed, indiv_ct);
  trimatrix_size = ((unsigned long)indiv_ct * (indiv_ct - 1)) / 2;
  g_reg_tot_xy = 0.0;
  g_reg_tot_x = 0.0;
  g_reg_tot_y = 0.0;
  g_reg_tot_xx = 0.0;
  g_reg_tot_yy = 0.0;
  rel_ptr = g_rel_dists;
  pheno_ptr = g_pheno_packed;
  g_jackknife_precomp = (double*)wkspace_alloc(indiv_ct * JACKKNIFE_VALS_REL * sizeof(double));
  if (!g_jackknife_precomp) {
    return RET_NOMEM;
  }
  fill_double_zero(g_jackknife_precomp, indiv_ct * JACKKNIFE_VALS_REL);
  for (uii = 1; uii < indiv_ct; uii++) {
    half_avg_pheno = *(++pheno_ptr);
    pheno_ptr2 = g_pheno_packed;
    jp_fixed_ptr = &(g_jackknife_precomp[uii * JACKKNIFE_VALS_REL]);
    jp_moving_ptr = g_jackknife_precomp;
    while (pheno_ptr2 < pheno_ptr) {
      dxx = (half_avg_pheno + (*pheno_ptr2++)) * 0.5;
      dyy = (*rel_ptr++);
      dxxyy = dxx * dyy;
      dxxsq = dxx * dxx;
      dyysq = dyy * dyy;
      g_reg_tot_xy += dxxyy;
      jp_fixed_ptr[0] += dxxyy;
      *jp_moving_ptr += dxxyy;
      jp_moving_ptr++;
      g_reg_tot_x += dxx;
      jp_fixed_ptr[1] += dxx;
      *jp_moving_ptr += dxx;
      jp_moving_ptr++;
      g_reg_tot_y += dyy;
      jp_fixed_ptr[2] += dyy;
      *jp_moving_ptr += dyy;
      jp_moving_ptr++;
      g_reg_tot_xx += dxxsq;
      jp_fixed_ptr[3] += dxxsq;
      *jp_moving_ptr += dxxsq;
      jp_moving_ptr++;
      g_reg_tot_yy += dyysq;
      jp_fixed_ptr[4] += dyysq;
      *jp_moving_ptr += dyysq;
      jp_moving_ptr++;
    }
  }
  trimatrix_size_recip = 1.0 / (double)trimatrix_size;
  sprintf(logbuf, "Regression slope (y = genomic relationship, x = avg phenotype): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y * trimatrix_size_recip) / (g_reg_tot_xx - g_reg_tot_x * g_reg_tot_x * trimatrix_size_recip));
  logprintb();
  sprintf(logbuf, "                 (y = avg phenotype, x = genomic relationship): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y * trimatrix_size_recip) / (g_reg_tot_yy - g_reg_tot_y * g_reg_tot_y * trimatrix_size_recip));
  logprintb();
  g_jackknife_iters = (regress_rel_iters + g_thread_ct - 1) / g_thread_ct;
  if (regress_rel_d) {
    g_jackknife_d = regress_rel_d;
  } else {
    g_jackknife_d = set_default_jackknife_d(indiv_ct);
  }
  g_geno = wkspace_alloc(g_thread_ct * CACHEALIGN(indiv_ct + (g_jackknife_d + 1) * sizeof(int)));
  if (!g_geno) {
    return RET_NOMEM;
  }
  for (ulii = 1; ulii < g_thread_ct; ulii++) {
    if (pthread_create(&(threads[ulii - 1]), NULL, &regress_rel_jack_thread, (void*)ulii)) {
      logprint(errstr_thread_create);
      while (--ulii) {
	pthread_join(threads[ulii - 1], NULL);
      }
      return RET_THREAD_CREATE_FAIL;
    }
  }
  ulii = 0;
  regress_rel_jack_thread((void*)ulii);
  dxx = g_calc_result[0][0]; // relationship on pheno
  dxxsq = g_calc_result[1][0];

  dyy = g_calc_result[2][0]; // pheno on relationship
  dyysq = g_calc_result[3][0];

  for (uii = 0; uii < g_thread_ct - 1; uii++) {
    pthread_join(threads[uii], NULL);
    dxx += g_calc_result[0][uii + 1];
    dxxsq += g_calc_result[1][uii + 1];
    dyy += g_calc_result[2][uii + 1];
    dyysq += g_calc_result[3][uii + 1];
  }
  ulii = g_jackknife_iters * g_thread_ct;
  putchar('\r');
  sprintf(logbuf, "Jackknife s.e. (y = genomic relationship): %g\n", sqrt((indiv_ct / (double)g_jackknife_d) * (dxxsq - dxx * dxx / (double)ulii) / ((double)ulii - 1)));
  logprintb();
  sprintf(logbuf, "               (y = phenotype): %g\n", sqrt((indiv_ct / (double)g_jackknife_d) * (dyysq - dyy * dyy / (double)ulii) / ((double)ulii - 1)));
  logprintb();
  return 0;
}

// Replaces matrix[][] with mult_val * matrix[][] + add_val * I.
// Multithreading doesn't help here.
void matrix_const_mult_add(double* matrix, double mult_val, double add_val) {
  unsigned int uii;
  unsigned int loop_end = g_indiv_ct - 1;
  unsigned int ujj;
  double* dptr = matrix;
#if __LP64__
  __m128d* vptr;
  __m128d v_mult_val = _mm_set1_pd(mult_val);
#endif
  for (uii = 0; uii < loop_end; uii++) {
    *dptr = (*dptr) * mult_val + add_val;
    dptr++;
#if __LP64__
    if ((unsigned long)dptr & 8) {
      *dptr *= mult_val;
      dptr++;
      ujj = 1;
    } else {
      ujj = 0;
    }
    vptr = (__m128d*)dptr;
    while (ujj < loop_end) {
      *vptr = _mm_mul_pd(*vptr, v_mult_val);
      vptr++;
      ujj += 2;
    }
    dptr = (double*)vptr;
    if (ujj < g_indiv_ct) {
      *dptr *= mult_val;
      dptr++;
    }
#else
    for (ujj = 0; ujj < g_indiv_ct; ujj++) {
      *dptr *= mult_val;
      dptr++;
    }
#endif
  }
  *dptr = (*dptr) * mult_val + add_val;
}

// sums[idx] = matrix[idx][1] + matrix[idx][2] + ....  Ideally, we can assume
// the matrix is symmetric and only reference the FORTRAN upper right, but for
// now we can't.
void matrix_row_sum_ur(double* sums, double* matrix) {
  unsigned int indiv_idx;
  double* dptr;
  double acc;
  double* sptr_end;
  double* sptr;
  fill_double_zero(sums, g_indiv_ct);
  for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
    dptr = &(matrix[indiv_idx * g_indiv_ct]);
    acc = 0.0;
    sptr_end = &(sums[indiv_idx]);
    sptr = sums;
    while (sptr < sptr_end) {
      acc += *dptr;
      *sptr += *dptr++;
      sptr++;
    }
    *sptr += acc + *dptr;
  }
}

inline double SQR(const double a) {
  return a * a;
}

#ifdef __cplusplus
inline double SIGN(const double &a, const double &b) {
  // PLINK helper.h SIGN() template specialized to doubles.
  return (b >= 0)? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
#else
inline double SIGN(const double a, const double b) {
  // PLINK helper.h SIGN() template specialized to doubles.
  return (b >= 0)? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
#endif

double pythag(const double a, const double b) {
  // PLINK stats.cpp pythag().
  double absa,absb;
 
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

int svdcmp_c(int m, double* a, double* w, double* v) {
  // C port of PLINK stats.cpp svdcmp().
  // Note that this function is NOT thread-safe, due to the buffer allocated
  // from the workspace stack.  Pass in a preallocated buffer if that's not
  // okay.
  unsigned char* wkspace_mark = wkspace_base;
  int n = m;
  int flag;
  int l = 0; // suppress compile warning
  int i,its,j,jj,k,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
  double volatile temp;
  double* rv1;
  if (wkspace_alloc_d_checked(&rv1, m * sizeof(double))) {
    return -1;
  }

  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a[k * m + i]);
      if (scale != 0.0) {
	for (k=i;k<m;k++) {
	  a[k * m + i] /= scale;
	  s += a[k * m + i]*a[k * m + i];
	}
	f=a[i * m + i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i * m + i]=f-g;
	for (j=l-1;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++) s += a[k * m + i]*a[k * m + j];
	  f=s/h;
	  for (k=i;k<m;k++) a[k * m + j] += f*a[k * m + i];
	}
	for (k=i;k<m;k++) a[k * m + i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i+1 <= m && i+1 != n) {
      for (k=l-1;k<n;k++) scale += fabs(a[i * m + k]);
      if (scale != 0.0) {
	for (k=l-1;k<n;k++) {
	  a[i * m + k] /= scale;
	  s += a[i * m + k]*a[i * m + k];
	}
	f=a[i * m + l-1];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i * m + l-1]=f-g;
	for (k=l-1;k<n;k++) rv1[k]=a[i * m + k]/h;
	for (j=l-1;j<m;j++) {
	  for (s=0.0,k=l-1;k<n;k++) s += a[j * m + k]*a[i * m + k];
	  for (k=l-1;k<n;k++) a[j * m + k] += s*rv1[k];
	}
	for (k=l-1;k<n;k++) a[i * m + k] *= scale;
      }
    }
    anorm=MAXV(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
	for (j=l;j<n;j++)
	  v[j * m + i]=(a[i * m + j]/a[i * m + l])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[i * m + k]*v[k * m + j];
	  for (k=l;k<n;k++) v[k * m + j] += s*v[k * m + i];
	}
      }
      for (j=l;j<n;j++) v[i * m + j]=v[j * m + i]=0.0;
    }
    v[i * m + i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=MINV(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i * m + j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) s += a[k * m + i]*a[k * m + j];
	f=(s/a[i * m + i])*g;
	for (k=i;k<m;k++) a[k * m + j] += f*a[k * m + i];
      }
      for (j=i;j<m;j++) a[j * m + i] *= g;
    } else for (j=i;j<m;j++) a[j * m + i]=0.0;
    ++a[i * m + i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=1;
      for (l=k;l>=0;l--) {
	nm=l-1;
	temp=fabs(rv1[l])+anorm;
	if (temp == anorm) {
	  flag=0;
	  break;
	}
	temp=fabs(w[nm])+anorm;
	if (temp == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<k+1;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  temp = fabs(f)+anorm;
	  if (temp == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j * m + nm];
	    z=a[j * m + i];
	    a[j * m + nm]=y*c+z*s;
	    a[j * m + i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) v[j * m + k] = -v[j * m + k];
	}
	break;
      }
      if (its == 29) 
	return 0; // cannot converge: multi-collinearity?
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj * m + j];
	  z=v[jj * m + i];
	  v[jj * m + j]=x*c+z*s;
	  v[jj * m + i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++) {
	  y=a[jj * m + j];
	  z=a[jj * m + i];
	  a[jj * m + j]=y*c+z*s;
	  a[jj * m + i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  wkspace_reset(wkspace_mark);
  return 1;
}

int invert_matrix_trunc_singular_svd(__CLPK_integer dim, double* matrix, double* dbl_1d_buf, double* dbl_2d_buf, __CLPK_integer min_dim) {
  // for no truncation, call this with min_dim == dim
  const double eps = 1e-24;
  unsigned char* wkspace_mark = wkspace_base;
  double* matrix_copy;
  double* small_matrix = NULL;
  unsigned int min_dim_u = min_dim;
  int flag;
  unsigned int cur_dim;
  unsigned int max_dim; // known singular
  unsigned int uii;
  int i;
  int j;
  int k;

  if (wkspace_alloc_d_checked(&matrix_copy, dim * dim * sizeof(double))) {
    return -1;
  }
  memcpy(matrix_copy, matrix, dim * dim * sizeof(double));
  flag = svdcmp_c(dim, matrix, dbl_1d_buf, dbl_2d_buf);

  if (!flag) {
    max_dim = dim;
    if (dim > min_dim + 1) {
      if (wkspace_alloc_d_checked(&small_matrix, (dim - 1) * (dim - 1) * sizeof(double))) {
        return -1;
      }
    }
    while (max_dim > min_dim_u + 1) {
      cur_dim = (max_dim + min_dim_u) / 2;
      for (uii = 0; uii < cur_dim; uii++) {
        memcpy(&(small_matrix[uii * cur_dim]), &(matrix_copy[uii * dim]), cur_dim * sizeof(double));
      }
      flag = svdcmp_c(cur_dim, small_matrix, dbl_1d_buf, dbl_2d_buf);
      if (flag) {
        max_dim = cur_dim;
      } else {
        min_dim_u = cur_dim;
      }
    }
    wkspace_reset(wkspace_mark);
    return min_dim_u;
  }
  // Look for singular values
  double wmax = 0;
  for (i=0; i<dim; i++) {
    wmax = dbl_1d_buf[i] > wmax ? dbl_1d_buf[i] : wmax;
  }
  double wmin = wmax * eps;
  for (i=0; i<dim; i++) {
    dbl_1d_buf[i] = dbl_1d_buf[i] < wmin ? 0 : 1/dbl_1d_buf[i];
  }

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      matrix[i * dim + j] = matrix[i * dim + j] * dbl_1d_buf[j];
    }
  }

  // [nxn].[t(v)] 
  for (i=0; i<dim; i++) {
    fill_double_zero(dbl_1d_buf, dim);
    for (j=0; j<dim; j++) {
      for (k=0; k<dim; k++) {
	dbl_1d_buf[j] += matrix[i * dim + k] * dbl_2d_buf[j * dim + k];
      }
    }
    for (j = 0; j < dim; j++) {
      matrix[i * dim + j] = dbl_1d_buf[j];
    }
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

#ifdef NOLAPACK
int invert_matrix(int dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* dbl_1d_buf, double* dbl_2d_buf) {
  // C port of PLINK stats.cpp's svd_inverse() function.

  // w -> dbl_1d_buf
  // v -> dbl_2d_buf
  const double eps = 1e-24;
  if (svdcmp_c(dim, matrix, dbl_1d_buf, dbl_2d_buf) == -1) {
    return -1;
  }

  // Look for singular values
  double wmax = 0;
  for (int i=0; i<dim; i++) {
    wmax = dbl_1d_buf[i] > wmax ? dbl_1d_buf[i] : wmax;
  }
  double wmin = wmax * eps;
  for (int i=0; i<dim; i++) {
    dbl_1d_buf[i] = dbl_1d_buf[i] < wmin ? 0 : 1/dbl_1d_buf[i];
  }
  
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      matrix[i * dim + j] = matrix[i * dim + j] * dbl_1d_buf[j];
    }
  }

  // [nxn].[t(v)] 
  for (int i=0; i<dim; i++) {
    fill_double_zero(dbl_1d_buf, dim);
    for (int j=0; j<dim; j++) {
      for (int k=0; k<dim; k++) {
	dbl_1d_buf[j] += matrix[i * dim + k] * dbl_2d_buf[j * dim + k];
      }
    }
    for (int j = 0; j < dim; j++) {
      matrix[i * dim + j] = dbl_1d_buf[j];
    }
  }
  return 0;
}

int invert_matrix_trunc_singular(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf, __CLPK_integer min_dim) {
  return invert_matrix_trunc_singular_svd(dim, matrix, int_1d_buf, dbl_2d_buf, min_dim);
}
#else
int invert_matrix(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf) {
  // dgetrf_/dgetri_ is more efficient than dpotrf_/dpotri_ on OS X.
  __CLPK_integer lwork = dim * dim;
  __CLPK_integer info;
  dgetrf_(&dim, &dim, matrix, &dim, int_1d_buf, &info);
  dgetri_(&dim, matrix, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  // todo: check if there are any major error conditions, return -1 on them
  return 0;
}

int invert_matrix_trunc_singular(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf, __CLPK_integer min_dim) {
  // Variant of invert_matrix which checks if the matrix is singular.  If it
  // is, this determines the minimum number of rows/columns which can be
  // removed from the bottom right to make the matrix nonsingular, knowing that
  // the upper left min_dim * min_dim matrix is nonsingular, and returns the
  // zero-based index of the first removed row/column.
  __CLPK_integer lwork = dim * dim;
  char cc = '1';
  unsigned char* wkspace_mark = wkspace_base;
  double* matrix_copy;
  __CLPK_integer info;
  double norm;
  double rcond;
  norm = dlange_(&cc, &dim, &dim, matrix, &dim, dbl_2d_buf);
  if (wkspace_alloc_d_checked(&matrix_copy, dim * dim * sizeof(double))) {
    return -1;
  }
  memcpy(matrix_copy, matrix, dim * dim * sizeof(double));
  dgetrf_(&dim, &dim, matrix, &dim, int_1d_buf, &info);
  if (info > 0) {
    rcond = 0.0;
  } else {
    dgecon_(&cc, &dim, matrix, &dim, &norm, &rcond, dbl_2d_buf, &(int_1d_buf[dim]), &info);
  }
  if (rcond < MATRIX_SINGULAR_RCOND) {
    memcpy(matrix, matrix_copy, dim * dim * sizeof(double));
    wkspace_reset(wkspace_mark);
    return invert_matrix_trunc_singular_svd(dim, matrix, (double*)int_1d_buf, dbl_2d_buf, min_dim);
  } else {
    dgetri_(&dim, matrix, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  }
  wkspace_reset(wkspace_mark);
  return 0;
}
#endif

#ifndef NOLAPACK
// one-trait REML via EM.
//
// wkbase is assumed to have space for three cache-aligned
// indiv_ct * indiv_ct double matrices plus three more rows.  The unpacked
// relationship matrix is stored in the SECOND slot.
void reml_em_one_trait(double* wkbase, double* pheno, double* covg_ref, double* covr_ref, double tol, int strict) {
  double ll_change;
  long long mat_offset = g_indiv_ct;
  double* rel_dists;
  MATRIX_INVERT_BUF1_TYPE* irow;
  __CLPK_integer lwork;
  double* row;
  double* row2;
  double* work;
  double* dptr;
  double* dptr2;
  double* matrix_pvg;
  __CLPK_integer indiv_ct_li = g_indiv_ct;
  double dxx;
  double dyy;
  double max_jump;
  double dzz;
  double dlg;
  double dle;
  double covg_cur_change = 1.0;
  double covr_cur_change = 1.0;
  double covg_last_change;
  double covr_last_change;
  double indiv_ct_d = 1 / (double)g_indiv_ct;
  unsigned int indiv_idx;
  int jj;
  mat_offset = CACHEALIGN_DBL(mat_offset * mat_offset);
  rel_dists = &(wkbase[mat_offset]);
  row = &(wkbase[mat_offset * 3]);
  irow = (MATRIX_INVERT_BUF1_TYPE*)row;
  row2 = &(row[g_indiv_ct]);
  work = &(wkbase[mat_offset * 2]);
  lwork = mat_offset;
  matrix_pvg = work;
  if (!lwork) {
    lwork = CACHELINE_DBL;
  }
  fill_double_zero(matrix_pvg, mat_offset);
  fill_double_zero(row2, g_indiv_ct);
  do {
    memcpy(wkbase, rel_dists, mat_offset * sizeof(double));
    matrix_const_mult_add(wkbase, *covg_ref, *covr_ref);
    invert_matrix(indiv_ct_li, wkbase, irow, work);
    matrix_row_sum_ur(row, wkbase);
    dxx = 0.0;
    dptr = row;
    dptr2 = &(row[g_indiv_ct]);
    while (dptr < dptr2) {
      dxx += *dptr++;
    }
    dxx = -1 / dxx;
    cblas_dger(CblasColMajor, g_indiv_ct, g_indiv_ct, dxx, row, 1, row, 1, wkbase, g_indiv_ct);
    // unfortunately, cblas_dsymm is much worse than cblas_dgemm on OS X
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, g_indiv_ct, g_indiv_ct, g_indiv_ct, 1.0, wkbase, g_indiv_ct, rel_dists, g_indiv_ct, 0.0, matrix_pvg, g_indiv_ct);
    dlg = 0.0;
    dle = 0.0;
    jj = g_indiv_ct + 1;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      dlg -= matrix_pvg[indiv_idx * jj];
      dle -= wkbase[indiv_idx * jj];
    }
    cblas_dsymv(CblasColMajor, CblasUpper, g_indiv_ct, 1.0, wkbase, g_indiv_ct, pheno, 1, 0.0, row2, 1);
    cblas_dsymv(CblasColMajor, CblasUpper, g_indiv_ct, 1.0, matrix_pvg, g_indiv_ct, row2, 1, 0.0, row, 1);
    dlg += cblas_ddot(g_indiv_ct, pheno, 1, row, 1);
    cblas_dsymv(CblasColMajor, CblasUpper, g_indiv_ct, 1.0, wkbase, g_indiv_ct, row2, 1, 0.0, row, 1);
    dle += cblas_ddot(g_indiv_ct, pheno, 1, row, 1);
    covg_last_change = covg_cur_change;
    covr_last_change = covr_cur_change;
    covg_cur_change = (*covg_ref) * (*covg_ref) * dlg * indiv_ct_d;
    covr_cur_change = (*covr_ref) * (*covr_ref) * dle * indiv_ct_d;
    if (strict) {
      max_jump = 1.0;
    } else {
      // acceleration factor:
      // min(half covg distance to 0 or 1, covr distance to 0 or 1, pi/4 divided
      // by last angular change, 1.0 / (1 - ratio of last two step lengths),
      // MAX_EM_ACCEL)
      dxx = atan2(covg_last_change, covr_last_change) - atan2(covg_cur_change, covr_cur_change);
      if (dxx < 0.0) {
	dxx = -dxx;
      }
      if (dxx > PI) {
	dxx = 2 * PI - dxx;
      }
      dyy = sqrt((covg_cur_change * covg_cur_change + covr_cur_change * covr_cur_change) / (covg_last_change * covg_last_change + covr_last_change * covr_last_change));
      if (covg_cur_change < 0.0) {
	max_jump = *covg_ref * (-0.5) / covg_cur_change;
      } else {
	max_jump = (1.0 - *covg_ref) * 0.5 / covg_cur_change;
      }
      if (covr_cur_change < 0.0) {
	dzz = *covr_ref * (-0.5) / covr_cur_change;
      } else {
	dzz = (1.0 - *covr_ref) * 0.5 / covr_cur_change;
      }
      if (dzz < max_jump) {
	max_jump = dzz;
      }
      dzz = (PI / 4) / dxx;
      if (dzz < max_jump) {
	max_jump = dzz;
      }
      if (dyy < 1.0) {
	dzz = 1 / (1.0 - dyy);
      }
      if (dzz < max_jump) {
	max_jump = dzz;
      }
      if (max_jump < 1.0) {
	max_jump = 1.0;
      } else if (max_jump > MAX_EM_ACCEL) {
	max_jump = MAX_EM_ACCEL;
      }
    }
    *covg_ref += covg_cur_change * max_jump;
    *covr_ref += covr_cur_change * max_jump;
    ll_change = (covg_cur_change * dlg) + (covr_cur_change * dle);
    printf("\b\b\b\b\b\b      \rcovg: %g  covr: %g  EM step log likelihood change: %g", *covg_ref, *covr_ref, ll_change);
    fflush(stdout);
  } while (ll_change > tol);
  printf("\n");
  sprintf(logbuf, "covg: %g  covr: %g\n", *covg_ref, *covr_ref);
  logstr(logbuf);
}
#endif

inline int is_founder(unsigned long* founder_info, int ii) {
  return ((!founder_info) || ((founder_info[ii / BITCT]) & (1LU << (ii % BITCT))));
}

typedef struct {
  char* family_ids;
  unsigned int max_family_id_len; // includes trailing null
  unsigned int* family_sizes;

  unsigned int* family_rel_space_offsets; // offset for rel_space lookup
  unsigned int* family_founder_cts;
  // direct indiv idx -> family idx lookup, to reduce number of bsearches
  unsigned int* family_idxs;

  // truncated triangular arrays of pedigree coefficient of relationship
  double* rel_space;

  // direct indiv idx -> rel_space idx lookup
  unsigned int* family_rel_nf_idxs;

  // following three variables are technically unnecessary for --genome, but we
  // get them for "free" in the process of calculating everything else, and
  // they'll be nice to use if we ever need to iterate by family in the future.
  unsigned int family_id_ct;
  // list of idxs of all individuals in first family, then second family, etc.
  unsigned int* family_info_space;
  unsigned int* family_info_offsets; // offset in family_info_space
} Pedigree_rel_info;

int populate_pedigree_rel_info(Pedigree_rel_info* pri_ptr, unsigned int unfiltered_indiv_ct, char* id_buf, char* person_ids, unsigned int max_person_id_len, char* paternal_ids, unsigned int max_paternal_id_len, char* maternal_ids, unsigned int max_maternal_id_len, unsigned long* founder_info) {
  unsigned char* wkspace_mark;
  unsigned char* wkspace_mark2;
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int oo;
  unsigned int uii;
  unsigned int initial_family_blocks;
  unsigned int indiv_uidx;
  long long llii;
  char* family_ids;
  char* cur_person_id;
  char* last_family_id = NULL;
  char* cur_family_id;
  unsigned int* family_sizes;
  unsigned int* uiptr;
  unsigned int* uiptr2 = NULL;
  unsigned int fidx;
  int family_size;
  unsigned int* remaining_indiv_idxs;
  int* remaining_indiv_parent_idxs; // -1 = no parent (or nonshared)
  int remaining_indiv_ct;
  int indiv_idx_write;
  int max_family_id_len = 0;
  char* indiv_ids;
  int max_indiv_id_len = 0;
  int max_pm_id_len;
  unsigned int family_id_ct;
  unsigned int* fis_ptr;
  char* stray_parent_ids;
  int stray_parent_ct;
  unsigned long* processed_indivs;
  int founder_ct;
  int max_family_nf = 0;
  unsigned int unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  unsigned int unfiltered_indiv_ctlm = unfiltered_indiv_ctl * BITCT;
  unsigned int* complete_indiv_idxs;
  unsigned int complete_indiv_idx_ct;
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
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_space), unfiltered_indiv_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_rel_nf_idxs), unfiltered_indiv_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_idxs), unfiltered_indiv_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_c_checked(&family_ids, unfiltered_indiv_ct * max_family_id_len)) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&family_sizes, unfiltered_indiv_ct * sizeof(int))) {
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
  initial_family_blocks = 1 + (unsigned int)(uiptr - family_sizes);
  if (qsort_ext(family_ids, initial_family_blocks, max_family_id_len, strcmp_deref, (char*)family_sizes, sizeof(int))) {
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
    family_sizes = (unsigned int*)wkspace_alloc(family_id_ct * sizeof(int));
    for (uii = 0; uii < family_id_ct; uii++) {
      family_sizes[uii] = *uiptr++;
    }
  }
  pri_ptr->family_ids = family_ids;
  pri_ptr->family_id_ct = family_id_ct;
  pri_ptr->max_family_id_len = max_family_id_len;
  pri_ptr->family_sizes = family_sizes;

  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_offsets), (family_id_ct + 1) * sizeof(int))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_rel_space_offsets), (family_id_ct + 1) * sizeof(int))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_founder_cts), family_id_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  fill_int_zero((int*)(pri_ptr->family_founder_cts), family_id_ct);

  ii = 0; // running family_info offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_info_offsets[fidx] = ii;
    ii += family_size;
  }

  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  memcpy(uiptr, pri_ptr->family_info_offsets, family_id_ct * sizeof(int));

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
    jj += ((long long)family_size * (family_size - 1) - (long long)kk * (kk - 1)) / 2;
  }
  // make it safe to determine size of blocks by subtracting from the next
  // offset, even if we're at the last family
  pri_ptr->family_info_offsets[family_id_ct] = unfiltered_indiv_ct;
  pri_ptr->family_rel_space_offsets[family_id_ct] = jj;
  if (wkspace_alloc_d_checked(&(pri_ptr->rel_space), jj * sizeof(double))) {
    return RET_NOMEM;
  }

  wkspace_mark = wkspace_base;
  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int))) {
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

  if (wkspace_alloc_ul_checked(&processed_indivs, (unfiltered_indiv_ctl + (max_family_nf + (BITCT2 - 1)) / BITCT2) * sizeof(long))) {
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
      memcpy(processed_indivs, founder_info, unfiltered_indiv_ctl * sizeof(long));
      if (wkspace_alloc_ui_checked(&complete_indiv_idxs, family_size * sizeof(int))) {
	return RET_NOMEM;
      }
      if (wkspace_alloc_ui_checked(&remaining_indiv_idxs, remaining_indiv_ct * sizeof(int))) {
	return RET_NOMEM;
      }
      if (wkspace_alloc_c_checked(&indiv_ids, family_size * max_indiv_id_len)) {
	return RET_NOMEM;
      }
      if (wkspace_alloc_i_checked(&remaining_indiv_parent_idxs, remaining_indiv_ct * 2 * sizeof(int))) {
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
      for (jj = 0; jj < remaining_indiv_ct; jj++) {
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
	    remaining_indiv_parent_idxs[jj * 2] = -2;
	  } else {
            remaining_indiv_parent_idxs[jj * 2] = fis_ptr[mm];
	  }
	} else {
          remaining_indiv_parent_idxs[jj * 2] = -1;
	}
	if (memcmp("0", &(maternal_ids[kk * max_maternal_id_len]), 2)) {
	  mm = bsearch_str(&(maternal_ids[kk * max_maternal_id_len]), indiv_ids, max_indiv_id_len, 0, family_size - 1);
	  if (mm == -1) {
	    strcpy(cur_person_id, &(maternal_ids[kk * max_maternal_id_len]));
	    cur_person_id = &(cur_person_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_indiv_parent_idxs[jj * 2 + 1] = -2;
	  } else {
	    remaining_indiv_parent_idxs[jj * 2 + 1] = fis_ptr[mm];
	  }
	} else {
	  remaining_indiv_parent_idxs[jj * 2 + 1] = -1;
	}
        remaining_indiv_idxs[jj] = kk;
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
      for (ii = 0; ii < remaining_indiv_ct; ii++) {
	jj = remaining_indiv_idxs[ii];
	if (remaining_indiv_parent_idxs[ii * 2] == -2) {
	  kk = bsearch_str(&(paternal_ids[jj * max_paternal_id_len]), stray_parent_ids, max_pm_id_len, 0, stray_parent_ct - 1);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[ii * 2] = kk;
	}
	if (remaining_indiv_parent_idxs[ii * 2 + 1] == -2) {
	  kk = bsearch_str(&(maternal_ids[jj * max_maternal_id_len]), stray_parent_ids, max_pm_id_len, 0, stray_parent_ct - 1);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[ii * 2 + 1] = kk;
	}
      }
      llii = (long long)founder_ct * (founder_ct - 1);
      while (remaining_indiv_ct) {
	indiv_idx_write = 0;
	for (ii = 0; ii < remaining_indiv_ct; ii++) {
	  kk = remaining_indiv_parent_idxs[ii * 2];
	  mm = remaining_indiv_parent_idxs[ii * 2 + 1];
	  jj = remaining_indiv_idxs[ii];
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
                dxx = 0.5 * rs_ptr[((long long)oo * (oo - 1) - llii) / 2 + nn];
	      }
	      if (is_set(founder_info, mm)) {
		if (mm == (int)complete_indiv_idxs[nn]) {
		  dxx += 0.5;
		}
	      } else if ((mm != -1) && (mm < (int)unfiltered_indiv_ct)) {
		oo = pri_ptr->family_rel_nf_idxs[mm];
		dxx += 0.5 * rs_ptr[((long long)oo * (oo - 1) - llii) / 2 + nn];
	      }
	      *rel_writer++ = dxx;
	    }
	    for (; nn < (int)complete_indiv_idx_ct; nn++) {
	      if (kk == -1) {
		dxx = 0.0;
	      } else if (kk >= (int)unfiltered_indiv_ct) {
		dxx = 0.5 * tmp_rel_space[(nn - founder_ct) * stray_parent_ct + kk - unfiltered_indiv_ctlm];
	      } else if (is_set(founder_info, kk)) {
                dxx = 0.5 * rs_ptr[((long long)nn * (nn - 1) - llii) / 2 + pri_ptr->family_rel_nf_idxs[kk]];
	      } else {
		oo = pri_ptr->family_rel_nf_idxs[kk];
		if (oo == nn) {
		  dxx = 0.5;
		} else if (oo < nn) {
		  dxx = 0.5 * rs_ptr[((long long)nn * (nn - 1) - llii) / 2 + oo];
		} else {
		  dxx = 0.5 * rs_ptr[((long long)oo * (oo - 1) - llii) / 2 + nn];
		}
	      }
	      if (mm >= (int)unfiltered_indiv_ct) {
		dxx += 0.5 * tmp_rel_space[(nn - founder_ct) * stray_parent_ct + mm - unfiltered_indiv_ctlm];
	      } else if (is_set(founder_info, mm)) {
		dxx += 0.5 * rs_ptr[((long long)nn * (nn - 1) - llii) / 2 + pri_ptr->family_rel_nf_idxs[mm]];
	      } else if (mm != -1) {
		oo = pri_ptr->family_rel_nf_idxs[mm];
		if (oo == nn) {
		  dxx += 0.5;
		} else if (oo < nn) {
		  dxx += 0.5 * rs_ptr[((long long)nn * (nn - 1) - llii) / 2 + oo];
		} else {
		  dxx += 0.5 * rs_ptr[((long long)oo * (oo - 1) - llii) / 2 + nn];
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

void count_genders(unsigned char* sex_info, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, int* male_ct_ptr, int* female_ct_ptr, int* unk_ct_ptr) {
  int male_ct = 0;
  int female_ct = 0;
  int unk_ct = 0;
  char cc;
  unsigned int indiv_uidx;
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    if (is_set(indiv_exclude, indiv_uidx)) {
      continue;
    }
    cc = sex_info[indiv_uidx];
    if (cc == 1) {
      male_ct++;
    } else if (cc == 2) {
      female_ct++;
    } else {
      unk_ct++;
    }
  }
  *male_ct_ptr = male_ct;
  *female_ct_ptr = female_ct;
  *unk_ct_ptr = unk_ct;
}

int load_pheno(FILE* phenofile, unsigned int unfiltered_indiv_ct, unsigned int indiv_exclude_ct, char* sorted_person_ids, unsigned int max_person_id_len, int* id_map, int missing_pheno, int missing_pheno_len, int affection_01, unsigned int mpheno_col, char* phenoname_str, char** pheno_c_ptr, double** pheno_d_ptr) {
  int affection = 1;
  char* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  int header_processed = 0;
  double missing_phenod = (double)missing_pheno;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned int indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  int person_idx;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  unsigned int tmp_len;
  unsigned int tmp_len2;
  unsigned int uii;
  char cc;
  double dxx;
  double dyy;
  if (pheno_d) {
    affection = 0;
  } else if (!pheno_c) {
    // no --pheno-merge, so set all to missing
    pheno_c = (char*)malloc(unfiltered_indiv_ct * sizeof(char));
    if (!pheno_c) {
      return RET_NOMEM;
    }
    memset(pheno_c, 255, unfiltered_indiv_ct);
    *pheno_c_ptr = pheno_c;
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
    if (is_eoln(*bufptr0)) {
      continue;
    }
    tmp_len = strlen_se(bufptr0);
    bufptr = next_item(bufptr0);
    if (no_more_items(bufptr)) {
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
	    if (no_more_items(bufptr)) {
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
	if (no_more_items(bufptr)) {
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
		pheno_c[person_idx] = -2;
	      } else {
	        pheno_c[person_idx] = -1;
	      }
	    } else if (affection_01) {
	      pheno_c[person_idx] = *bufptr - '0';
	    } else {
	      pheno_c[person_idx] = *bufptr - '1';
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
	      cc = pheno_c[uii];
	      if (cc == -1) {
		pheno_d[uii] = missing_phenod;
	      } else if (!cc) {
		pheno_d[uii] = dxx;
	      } else if (cc == 1) {
		pheno_d[uii] = dyy;
	      } else {
		// -2 means phenotype entry was '0'
		pheno_d[uii] = 0.0;
	      }
	    }
	    free(pheno_c);
	    *pheno_c_ptr = NULL;
	    affection = 0;
	  }
	}
	if (!affection) {
	  if (sscanf(bufptr, "%lg", &(pheno_d[person_idx])) != 1) {
	    pheno_d[person_idx] = missing_phenod;
	  }
	}
      }
    }
  }
  if (pheno_c) {
    for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
      if (pheno_c[uii] == -2) {
	pheno_c[uii] = -1;
      }
    }
  }
  if (!feof(phenofile)) {
    return RET_READ_FAIL;
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int makepheno_load(FILE* phenofile, char* makepheno_str, unsigned int unfiltered_indiv_ct, unsigned int indiv_exclude_ct, char* sorted_person_ids, unsigned int max_person_id_len, int* id_map, char** pheno_c_ptr) {
  unsigned int indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  unsigned int mp_strlen = strlen(makepheno_str);
  int makepheno_all = ((mp_strlen == 1) && (makepheno_str[0] == '*'));
  unsigned char* wkspace_mark = wkspace_base;
  char* pheno_c = *pheno_c_ptr;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  int person_idx;
  unsigned int tmp_len;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    return RET_NOMEM;
  }
  if (!pheno_c) {
    pheno_c = (char*)malloc(unfiltered_indiv_ct * sizeof(char));
    if (!pheno_c) {
      return RET_NOMEM;
    }
    memset(pheno_c, makepheno_all? 0 : -1, unfiltered_indiv_ct);
    *pheno_c_ptr = pheno_c;
  }
  tbuf[MAXLINELEN - 1] = ' ';  
  while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Excessively long line in phenotype file (max %d chars).\n", MAXLINELEN - 3);
      logprintb();
      return RET_INVALID_FORMAT;
    }
    bufptr0 = skip_initial_spaces(tbuf);
    if (is_eoln(*bufptr0)) {
      continue;
    }
    bufptr = next_item(bufptr0);
    if (no_more_items(bufptr)) {
      logprint(errstr_phenotype_format);
      return RET_INVALID_FORMAT;
    }
    person_idx = bsearch_fam_indiv(id_buf, sorted_person_ids, max_person_id_len, indiv_ct, bufptr0, bufptr);
    if (person_idx != -1) {
      person_idx = id_map[person_idx];
      if (makepheno_all) {
        pheno_c[person_idx] = 1;
      } else {
	bufptr = next_item(bufptr);
        tmp_len = strlen_se(bufptr);
	if ((tmp_len == mp_strlen) && (!memcmp(bufptr, makepheno_str, mp_strlen))) {
	  pheno_c[person_idx] = 1;
	} else {
	  pheno_c[person_idx] = 0;
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

int convert_tail_pheno(unsigned int unfiltered_indiv_ct, char** pheno_c_ptr, double** pheno_d_ptr, double tail_bottom, double tail_top, double missing_phenod) {
  char* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  unsigned int uii;
  double dxx;
  if (!(*pheno_d_ptr)) {
    logprint("Error: --tail-pheno requires scalar phenotype data.\n");
    return RET_INVALID_FORMAT;
  }
  if (!pheno_c) {
    pheno_c = (char*)malloc(unfiltered_indiv_ct * sizeof(char));
    if (!pheno_c) {
      return RET_NOMEM;
    }
    *pheno_c_ptr = pheno_c;
  }
  for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
    dxx = pheno_d[uii];
    if (dxx == missing_phenod) {
      pheno_c[uii] = -1;
    } else if (dxx <= tail_bottom) {
      pheno_c[uii] = 0;
    } else if (dxx > tail_top) {
      pheno_c[uii] = 1;
    } else {
      pheno_c[uii] = -1;
    }
  }
  free(pheno_d);
  *pheno_d_ptr = NULL;
  return 0;
}

void prune_missing_phenos(unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int* indiv_exclude_ct_ptr, char* pheno_c, double* pheno_d, double missing_phenod) {
  unsigned int uii;
  if (pheno_c) {
    for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
      if (pheno_c[uii] == -1) {
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

inline void set_surrounding_bits(char* id_buf, int id_idx, unsigned long* exclude_arr, char* sorted_ids, unsigned long max_id_len, int* id_map, int sorted_ids_len, unsigned int* exclude_ct_ptr) {
  unsigned int id_len = strlen(id_buf) + 1;
  int id_idx2 = id_idx - 1;
  while ((id_idx2 >= 0) && (!memcmp(id_buf, &(sorted_ids[id_idx2 * max_id_len]), id_len))) {
    set_bit(exclude_arr, id_map[id_idx2--], exclude_ct_ptr);
  }
  id_idx2 = id_idx + 1;
  while ((id_idx2 < sorted_ids_len) && (!memcmp(id_buf, &(sorted_ids[id_idx2 * max_id_len]), id_len))) {
    set_bit(exclude_arr, id_map[id_idx2++], exclude_ct_ptr);
  }
}

inline void clear_surrounding_bits(char* id_buf, int id_idx, unsigned long* exclude_arr_new, char* sorted_ids, unsigned long max_id_len, int* id_map, int sorted_ids_len, unsigned int* include_ct_ptr) {
  unsigned int id_len = strlen(id_buf) + 1;
  int id_idx2 = id_idx - 1;
  while ((id_idx2 >= 0) && (!memcmp(id_buf, &(sorted_ids[id_idx2 * max_id_len]), id_len))) {
    clear_bit(exclude_arr_new, id_map[id_idx2--], include_ct_ptr);
  }
  id_idx2 = id_idx + 1;
  while ((id_idx2 < sorted_ids_len) && (!memcmp(id_buf, &(sorted_ids[id_idx2 * max_id_len]), id_len))) {
    clear_bit(exclude_arr_new, id_map[id_idx2++], include_ct_ptr);
  }
}

int include_or_exclude(char* fname, char* sorted_ids, int sorted_ids_len, unsigned long max_id_len, int* id_map, int unfiltered_ct, unsigned long* exclude_arr, unsigned int* exclude_ct_ptr, int indivs, int do_exclude, int duplicates_present) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned long* exclude_arr_new = NULL;
  unsigned int include_ct = 0;
  int unfiltered_ctl = (unfiltered_ct + (BITCT - 1)) / BITCT;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  int ii;
  int jj;

  if (!do_exclude) {
    if (wkspace_alloc_ul_checked(&exclude_arr_new, unfiltered_ctl * sizeof(long))) {
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
      if (is_eoln(*bufptr0)) {
	continue;
      }
      bufptr = next_item(tbuf);
      if (no_more_items(bufptr)) {
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
    memcpy(exclude_arr, exclude_arr_new, unfiltered_ctl * sizeof(long));
    *exclude_ct_ptr = unfiltered_ct - include_ct;
  }
  wkspace_reset(wkspace_mark);
  fclose(infile);
  return 0;
}

int filter_indivs_file(char* filtername, char* sorted_person_ids, int sorted_ids_len, unsigned int max_person_id_len, int* id_map, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int* indiv_exclude_ct_ptr, char* filterval, int mfilter_col) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  int fv_strlen = strlen(filterval);
  unsigned int unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  unsigned int include_ct = 0;
  unsigned long* indiv_exclude_new;
  char* id_buf;
  char* bufptr;
  int person_idx;
  int ii;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&indiv_exclude_new, unfiltered_indiv_ctl * sizeof(long))) {
    return RET_NOMEM;
  }
  fill_ulong_one(indiv_exclude_new, unfiltered_indiv_ctl);

  if (fopen_checked(&infile, filtername, "r")) {
    return RET_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (is_eoln(*tbuf)) {
      continue;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Excessively long line in --keep/--remove file (max %d chars).\n", MAXLINELEN - 3);
      logprintb();
      fclose(infile);
      return RET_INVALID_FORMAT;
    }
    bufptr = next_item(tbuf);
    if (no_more_items(bufptr)) {
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
	if (no_more_items(bufptr)) {
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
  memcpy(indiv_exclude, indiv_exclude_new, unfiltered_indiv_ctl * sizeof(long));
  *indiv_exclude_ct_ptr = unfiltered_indiv_ct - include_ct;

  wkspace_reset(wkspace_mark);
  fclose(infile);
  return 0;
}

int filter_indivs_var(unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int* indiv_exclude_ct_ptr, char* pheno_c, unsigned long* founder_info, int filter_setting) {
  int filter_cc = (founder_info == NULL);
  unsigned char* wkspace_mark = wkspace_base;
  unsigned int unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  unsigned int include_ct = 0;
  char cc_match = (char)filter_setting;
  unsigned long* indiv_exclude_new;
  unsigned int indiv_uidx;
  if (wkspace_alloc_ul_checked(&indiv_exclude_new, unfiltered_indiv_ctl * sizeof(long))) {
    return RET_NOMEM;
  }
  fill_ulong_one(indiv_exclude_new, unfiltered_indiv_ctl);

  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    if (!is_set(indiv_exclude, indiv_uidx)) {
      if (filter_cc) {
	if (pheno_c[indiv_uidx] == cc_match) {
	  clear_bit(indiv_exclude_new, indiv_uidx, &include_ct);
	}
      } else if (filter_setting == 1) {
	// yeah, this could be reduced to bitwise-anding longs plus a popcount
	if (is_set(founder_info, indiv_uidx)) {
          clear_bit(indiv_exclude_new, indiv_uidx, &include_ct);
	}
      } else {
	if (!is_set(founder_info, indiv_uidx)) {
	  clear_bit(indiv_exclude_new, indiv_uidx, &include_ct);
	}
      }
    }
  }
  memcpy(indiv_exclude, indiv_exclude_new, unfiltered_indiv_ctl * sizeof(long));
  *indiv_exclude_ct_ptr = unfiltered_indiv_ct - include_ct;
  wkspace_reset(wkspace_mark);
  return 0;
}

int mind_filter(FILE* pedfile, double mind_thresh, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_exclude_ct, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int* indiv_exclude_ct_ptr, int bed_offset, unsigned long long* line_mids, int pedbuflen, char missing_geno) {
  int binary_files = (line_mids == NULL);
  unsigned int mind_int_thresh = (int)(mind_thresh * (unfiltered_marker_ct - marker_exclude_ct) + EPSILON);
  unsigned int marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  unsigned int unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned int removed_ct = 0;
  unsigned char* loadbuf;
  unsigned char* cptr;
  unsigned int* missing_cts;
  unsigned int indiv_uidx;
  unsigned int marker_uidx;
  unsigned int marker_idx;
  unsigned int uii;
  unsigned int ujj;
  unsigned char ucc;

  if (binary_files) {
    if (wkspace_alloc_ui_checked(&missing_cts, unfiltered_indiv_ct * sizeof(int))) {
      return RET_NOMEM;
    }
    if (wkspace_alloc_uc_checked(&loadbuf, unfiltered_indiv_ct4)) {
      return RET_NOMEM;
    }
    fill_uint_zero(missing_cts, unfiltered_indiv_ct);
    if (fseeko(pedfile, bed_offset, SEEK_SET)) {
      return RET_READ_FAIL;
    }
    marker_idx = 0;
    for (marker_uidx = 0; marker_idx < marker_ct; marker_uidx++) {
      if (is_set(marker_exclude, marker_uidx)) {
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	if (fseeko(pedfile, bed_offset + (marker_uidx * unfiltered_indiv_ct4), SEEK_SET)) {
	  return RET_READ_FAIL;
	}
      }
      marker_idx++;
      if (fread(loadbuf, 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	return RET_READ_FAIL;
      }
      cptr = loadbuf;
      indiv_uidx = 0;
      ujj = 0;
      for (uii = 0; uii < unfiltered_indiv_ct4; uii++) {
        ucc = *cptr++;
	ujj += 4; // may overshoot on last byte
	do {
          if ((ucc & 3) == 1) {
	    missing_cts[indiv_uidx] += 1;
	  }
	  ucc >>= 2;
	} while (++indiv_uidx < ujj);
      }
    }
    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
      if (is_set(indiv_exclude, indiv_uidx)) {
	continue;
      }
      if (missing_cts[indiv_uidx] > mind_int_thresh) {
	set_bit(indiv_exclude, indiv_uidx, indiv_exclude_ct_ptr);
	removed_ct++;
      }
    }
  } else {
    if (wkspace_alloc_uc_checked(&loadbuf, pedbuflen)) {
      return RET_NOMEM;
    }
    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
      if (is_set(indiv_exclude, indiv_uidx)) {
	continue;
      }
      if (fseeko(pedfile, line_mids[indiv_uidx], SEEK_SET)) {
	return RET_READ_FAIL;
      }
      // could be slightly more efficient if line lengths are saved during
      // first load, and fread() is used afterward, but whatever.
      if (fgets((char*)loadbuf, pedbuflen, pedfile) == NULL) {
	return RET_READ_FAIL;
      }
      cptr = loadbuf;
      uii = 0;
      for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
	ucc = *cptr;
	if (!ucc) {
	  sprintf(logbuf, "Error: Early line termination (.ped indiv %u, marker %u).\n", indiv_uidx + 1, marker_uidx + 1);
	  logprintb();
	  return RET_INVALID_FORMAT;
	}
	cptr++;
	while ((*cptr == ' ') || (*cptr == '\t')) {
	  cptr++;
	}
	if (!is_set(marker_exclude, marker_uidx)) {
	  if (ucc == (unsigned char)missing_geno) {
	    if (*cptr != missing_geno) {
	      sprintf(logbuf, "Error: 1st allele missing, 2nd isn't (.ped indiv %u, marker %u).\n", indiv_uidx + 1, marker_uidx + 1);
	      logprintb();
	      return RET_INVALID_FORMAT;
	    }
	    if (!is_set(marker_exclude, marker_uidx)) {
	      uii++;
	    }
	  } else {
	    ucc = *cptr;
	    if ((!ucc) || (ucc == (unsigned char)missing_geno)) {
	      sprintf(logbuf, "Error: 1st allele present, 2nd isn't (.ped indiv %u, marker %u).\n", indiv_uidx + 1, marker_uidx + 1);
	      logprintb();
	      return RET_INVALID_FORMAT;
	    }
	  }
	}
	cptr++;
	while ((*cptr == ' ') || (*cptr == '\t')) {
	  cptr++;
	}
      }
      if (uii > mind_int_thresh) {
	set_bit(indiv_exclude, indiv_uidx, indiv_exclude_ct_ptr);
	removed_ct++;
      }
    }
  }
  wkspace_reset(wkspace_mark);
  sprintf(logbuf, "%u individual%s removed due to missing genotype data (--mind).\n", removed_ct, (removed_ct == 1)? "" : "s");
  logprintb();
  return 0;
}

int incr_text_allele(char cc, char* marker_alleles, unsigned int* marker_allele_cts, int is_founder, unsigned int* marker_nf_allele_cts) {
  int ii;
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

int calc_freqs_and_binary_hwe(FILE* pedfile, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned long* founder_info, int nonfounders, int maf_succ, double* set_allele_freqs, char** marker_alleles_ptr, unsigned long** marker_reverse_ptr, unsigned int** marker_allele_cts_ptr, unsigned int** missing_cts_ptr, int bed_offset, unsigned long long* line_mids, int pedbuflen, unsigned char missing_geno, int hwe_all, char* pheno_c, int** hwe_lls_ptr, int** hwe_lhs_ptr, int** hwe_hhs_ptr, int** hwe_ll_allfs_ptr, int** hwe_lh_allfs_ptr, int** hwe_hh_allfs_ptr, int** ll_cts_ptr, int** lh_cts_ptr, int** hh_cts_ptr, int wt_needed, unsigned char** marker_weights_base_ptr, double** marker_weights_ptr) {
  int binary_files = (line_mids == NULL);
  unsigned int unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  unsigned char* wkspace_mark;
  unsigned char* loadbuf;
  unsigned char* cptr;
  unsigned int indiv_uidx;
  unsigned int marker_uidx;
  unsigned int uii;
  int ii;
  unsigned int ujj;
  unsigned int ukk;
  unsigned int set_allele_ct;
  unsigned int missing_ct;
  unsigned int clear_allele_ct;
  unsigned int* marker_allele_cts;
  unsigned int* marker_nf_allele_cts;
  char* marker_alleles = NULL;
  unsigned int* missing_cts;
  unsigned char ucc = '\0';
  double* marker_weights;
  int* hwe_lls;
  int* hwe_lhs;
  int* hwe_hhs;
  int* hwe_ll_allfs;
  int* hwe_lh_allfs;
  int* hwe_hh_allfs;
  int* ll_cts;
  int* lh_cts;
  int* hh_cts;
  int hwe_ll;
  int hwe_lh;
  int hwe_hh;
  int hwe_ll_allf;
  int hwe_lh_allf;
  int hwe_hh_allf;
  int ll_extra;
  int lh_extra;
  int hh_extra;
  if (wt_needed) {
    if (wkspace_alloc_uc_checked(marker_weights_base_ptr, CACHELINE)) {
      return RET_NOMEM;
    }
    if (wkspace_alloc_d_checked(marker_weights_ptr, unfiltered_marker_ct * sizeof(double))) {
      return RET_NOMEM;
    }
    marker_weights = *marker_weights_ptr;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      marker_weights[marker_uidx] = -1.0;
    }
  }
  uii = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(marker_reverse_ptr, uii * sizeof(long))) {
    return RET_NOMEM;
  }
  fill_ulong_zero(*marker_reverse_ptr, uii);
  if (!binary_files) {
    // N.B. marker_alleles initially needs (4 * unfiltered_marker_ct) bytes,
    // but is then collapsed down to (2 * unfiltered_marker_ct) before the HWE
    // counts are determined.
    if (wkspace_alloc_c_checked(&marker_alleles, 2 * unfiltered_marker_ct * sizeof(char))) {
      return RET_NOMEM;
    }
    *marker_alleles_ptr = marker_alleles;
  }

  if (wkspace_alloc_i_checked(&hwe_lls, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  *hwe_lls_ptr = hwe_lls;
  if (wkspace_alloc_i_checked(&hwe_lhs, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  *hwe_lhs_ptr = hwe_lhs;
  if (wkspace_alloc_i_checked(&hwe_hhs, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  *hwe_hhs_ptr = hwe_hhs;

  if (wkspace_alloc_i_checked(&hwe_ll_allfs, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  *hwe_ll_allfs_ptr = hwe_ll_allfs;
  if (wkspace_alloc_i_checked(&hwe_lh_allfs, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  *hwe_lh_allfs_ptr = hwe_lh_allfs;
  if (wkspace_alloc_i_checked(&hwe_hh_allfs, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  *hwe_hh_allfs_ptr = hwe_hh_allfs;

  if (wkspace_alloc_i_checked(&ll_cts, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  *ll_cts_ptr = ll_cts;
  if (wkspace_alloc_i_checked(&lh_cts, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  *lh_cts_ptr = lh_cts;
  if (wkspace_alloc_i_checked(&hh_cts, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  *hh_cts_ptr = hh_cts;

  if (binary_files) {
    if (!pheno_c) {
      hwe_all = 1;
    }
    if (wkspace_alloc_ui_checked(&marker_allele_cts, 2 * unfiltered_marker_ct * sizeof(int))) {
      return RET_NOMEM;
    }

    fill_uint_zero(marker_allele_cts, 2 * unfiltered_marker_ct);
    *marker_allele_cts_ptr = marker_allele_cts;
    wkspace_mark = wkspace_base;
    if (wkspace_alloc_uc_checked(&loadbuf, unfiltered_indiv_ct4)) {
      return RET_NOMEM;
    }
    if (fseeko(pedfile, bed_offset, SEEK_SET)) {
      return RET_READ_FAIL;
    }
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      if (is_set(marker_exclude, marker_uidx)) {
	do {
	  marker_uidx++;
	} while ((marker_uidx < unfiltered_marker_ct) && is_set(marker_exclude, marker_uidx));
	if (marker_uidx == unfiltered_marker_ct) {
	  break;
	}
	if (fseeko(pedfile, bed_offset + (marker_uidx * unfiltered_indiv_ct4), SEEK_SET)) {
	  return RET_READ_FAIL;
	}
      }
      if (fread(loadbuf, 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	return RET_READ_FAIL;
      }
      cptr = loadbuf;
      set_allele_ct = 0;
      clear_allele_ct = 0;
      missing_ct = 0;
      hwe_ll = 0;
      hwe_lh = 0;
      hwe_hh = 0;
      hwe_ll_allf = 0;
      hwe_lh_allf = 0;
      hwe_hh_allf = 0;
      ll_extra = 0;
      lh_extra = 0;
      hh_extra = 0;
      for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
	if (indiv_uidx & 3) {
          ucc >>= 2;
	} else {
	  ucc = *cptr++;
	}
	if (!is_set(indiv_exclude, indiv_uidx)) {
	  uii = ucc & 3;
	  if (uii == 1) {
	    missing_ct++;
	  } else if (nonfounders || is_set(founder_info, indiv_uidx)) {
	    if (uii == 3) {
	      if (hwe_all || (!pheno_c[indiv_uidx])) {
	        hwe_hh++;
	      } else {
		hwe_hh_allf++;
	      }
	      set_allele_ct += 2;
	    } else if (uii == 2) {
	      if (hwe_all || (!pheno_c[indiv_uidx])) {
	        hwe_lh++;
	      } else {
		hwe_lh_allf++;
	      }
	      set_allele_ct++;
	      clear_allele_ct++;
	    } else {
	      if (hwe_all || (!pheno_c[indiv_uidx])) {
	        hwe_ll++;
	      } else {
		hwe_ll_allf++;
	      }
	      clear_allele_ct += 2;
	    }
	  } else {
	    if (uii == 3) {
	      hh_extra++;
	      marker_allele_cts[2 * marker_uidx + 1] += 2;
	    } else if (uii == 2) {
	      lh_extra++;
	      marker_allele_cts[2 * marker_uidx] += 1;
	      marker_allele_cts[2 * marker_uidx + 1] += 1;
	    } else {
	      ll_extra++;
	      marker_allele_cts[2 * marker_uidx] += 2;
	    }
	  }
	}
      }
      hwe_lls[marker_uidx] = hwe_ll;
      hwe_lhs[marker_uidx] = hwe_lh;
      hwe_hhs[marker_uidx] = hwe_hh;
      hwe_ll_allfs[marker_uidx] = hwe_ll + hwe_ll_allf;
      hwe_lh_allfs[marker_uidx] = hwe_lh + hwe_lh_allf;
      hwe_hh_allfs[marker_uidx] = hwe_hh + hwe_hh_allf;
      ll_cts[marker_uidx] = hwe_ll + hwe_ll_allf + ll_extra;
      lh_cts[marker_uidx] = hwe_lh + hwe_lh_allf + lh_extra;
      hh_cts[marker_uidx] = hwe_hh + hwe_hh_allf + hh_extra;
      marker_allele_cts[2 * marker_uidx] += clear_allele_ct;
      marker_allele_cts[2 * marker_uidx + 1] += set_allele_ct;
      uii = set_allele_ct + clear_allele_ct + 2 * maf_succ;
      if (!uii) {
	// avoid 0/0 division
	set_allele_freqs[marker_uidx] = 0.5;
      } else {
	set_allele_freqs[marker_uidx] = ((double)(set_allele_ct + maf_succ)) / ((double)uii);
      }
    }
  } else {
    memset(marker_alleles, 0, 4 * unfiltered_marker_ct);
    *marker_alleles_ptr = marker_alleles;
    if (wkspace_alloc_ui_checked(&marker_allele_cts, 4 * unfiltered_marker_ct * sizeof(int))) {
      return RET_NOMEM;
    }
    fill_uint_zero(marker_allele_cts, 4 * unfiltered_marker_ct);
    *marker_allele_cts_ptr = marker_allele_cts;
    if (wkspace_alloc_ui_checked(&missing_cts, unfiltered_marker_ct * sizeof(int))) {
      return RET_NOMEM;
    }
    fill_uint_zero(missing_cts, unfiltered_marker_ct);
    *missing_cts_ptr = missing_cts;
    wkspace_mark = wkspace_base;
    if (wkspace_alloc_ui_checked(&marker_nf_allele_cts, 4 * unfiltered_marker_ct * sizeof(int))) {
      return RET_NOMEM;
    }
    fill_uint_zero(marker_nf_allele_cts, 4 * unfiltered_marker_ct);
    if (wkspace_alloc_uc_checked(&loadbuf, pedbuflen)) {
      return RET_NOMEM;
    }

    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
      if (is_set(indiv_exclude, indiv_uidx)) {
	continue;
      }
      if (fseeko(pedfile, line_mids[indiv_uidx], SEEK_SET)) {
	return RET_READ_FAIL;
      }
      if (fgets((char*)loadbuf, pedbuflen, pedfile) == NULL) {
	return RET_READ_FAIL;
      }
      cptr = loadbuf;
      uii = (nonfounders || is_set(founder_info, indiv_uidx));
      for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
	ucc = *cptr;
	if (!ucc) {
	  sprintf(logbuf, "Error: Early line termination (.ped indiv %u, marker %u).\n", indiv_uidx + 1, marker_uidx + 1);
	  logprintb();
	  return RET_INVALID_FORMAT;
	}
	cptr++;
	while ((*cptr == ' ') || (*cptr == '\t')) {
	  cptr++;
	}
	if (!is_set(marker_exclude, marker_uidx)) {
	  if (ucc == missing_geno) {
	    if (*cptr != missing_geno) {
	      sprintf(logbuf, "Error: 1st allele missing, 2nd isn't (.ped indiv %u, marker %u).\n", indiv_uidx + 1, marker_uidx + 1);
	      logprintb();
	      return RET_INVALID_FORMAT;
	    }
	    missing_cts[marker_uidx] += 1;
	  } else {
	    if (incr_text_allele((char)ucc, &(marker_alleles[4 * marker_uidx]), &(marker_allele_cts[4 * marker_uidx]), uii, &(marker_nf_allele_cts[4 * marker_uidx]))) {
	      sprintf(logbuf, "Error: More than four different allele codes at marker %u.\n", marker_uidx + 1);
	      logprintb();
	      return RET_INVALID_FORMAT;
	    }
	    ucc = *cptr;
	    if ((!ucc) || (ucc == missing_geno)) {
	      sprintf(logbuf, "Error: 1st allele present, 2nd isn't (.ped indiv %u, marker %u).\n", indiv_uidx + 1, marker_uidx + 1);
	      logprintb();
	      return RET_INVALID_FORMAT;
	    }
	    if (incr_text_allele((char)ucc, &(marker_alleles[4 * marker_uidx]), &(marker_allele_cts[4 * marker_uidx]), uii, &(marker_nf_allele_cts[4 * marker_uidx]))) {
	      sprintf(logbuf, "Error: More than four different allele codes at marker %u.\n", marker_uidx + 1);
	      logprintb();
	      return RET_INVALID_FORMAT;
	    }
	  }
	}
	cptr++;
	while ((*cptr == ' ') || (*cptr == '\t')) {
	  cptr++;
	}
      }
    }

    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      if (is_set(marker_exclude, marker_uidx)) {
	continue;
      }
      // insertion sort
      for (uii = 1; uii < 4; uii++) {
	ujj = marker_allele_cts[4 * marker_uidx + uii];
	if (marker_allele_cts[4 * marker_uidx + uii - 1] < ujj) {
	  ucc = (unsigned char)marker_alleles[4 * marker_uidx + uii];
	  ukk = marker_nf_allele_cts[4 * marker_uidx + uii];
	  ii = uii;
	  do {
	    ii--;
	    marker_alleles[4 * marker_uidx + ii + 1] = marker_alleles[4 * marker_uidx + ii];
	    marker_allele_cts[4 * marker_uidx + ii + 1] = marker_allele_cts[4 * marker_uidx + ii];
	    marker_nf_allele_cts[4 * marker_uidx + ii + 1] = marker_nf_allele_cts[4 * marker_uidx + ii];
	  } while ((ii > 0) && (marker_allele_cts[4 * marker_uidx + ii - 1] < ujj));
	  marker_alleles[4 * marker_uidx + ii] = (char)ucc;
	  marker_allele_cts[4 * marker_uidx + ii] = ujj;
	  marker_nf_allele_cts[4 * marker_uidx + ii] = ukk;
	}
      }
      // todo: log triallelic/quadallelic marker IDs?
      uii = marker_allele_cts[4 * marker_uidx] + marker_allele_cts[4 * marker_uidx + 1] + 2 * maf_succ;
      if (!uii) {
	set_allele_freqs[marker_uidx] = 0.5;
      } else {
	set_allele_freqs[marker_uidx] = ((double)(marker_allele_cts[4 * marker_uidx] + maf_succ)) / ((double)uii);
      }
      for (uii = 0; uii < 4; uii++) {
	marker_allele_cts[4 * marker_uidx + uii] += marker_nf_allele_cts[4 * marker_uidx + uii];
      }
    }
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int load_one_freq(char allele_min, char allele_maj, double maf, double* set_allele_freq_ptr, char* marker_alleles, unsigned int* marker_allele_cts, char missing_geno) {
  int ii;
  char cc;
  if (maf > 0.5) {
    cc = allele_maj;
    allele_maj = allele_min;
    allele_min = cc;
    maf = 1.0 - maf;
  }
  if (!marker_allele_cts) {
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
  } else {
    *set_allele_freq_ptr = 1.0 - maf;
    if (allele_min == missing_geno) {
      if (allele_maj != marker_alleles[0]) {
	// don't actually need to shuffle counts if allele_maj ==
        // marker_alleles[1], since only ([0] + [1]) and ([2] + [3]) sums
	// matter
	if (allele_maj != marker_alleles[1]) {
          if (allele_maj == marker_alleles[2]) {
	    marker_allele_cts[0] += marker_allele_cts[2];
	    marker_allele_cts[2] = marker_allele_cts[1];
	    marker_allele_cts[1] = 0;
	  } else if (allele_maj == marker_alleles[3]) {
	    marker_allele_cts[0] += marker_allele_cts[3];
	    marker_allele_cts[3] = marker_allele_cts[1];
	    marker_allele_cts[1] = 0;
	  } else {
	    return -1;
	  }
	}
	marker_alleles[1] = marker_alleles[0];
	marker_alleles[0] = allele_maj;
      }
    } else {
      for (ii = 0; ii < 4; ii++) {
	if (allele_maj == marker_alleles[ii]) {
	  break;
	}
      }
      if (ii == 4) {
	return -1;
      }
      if (marker_alleles[1]) {
	if (allele_min == marker_alleles[0]) {
	  if (ii > 1) {
	    marker_allele_cts[0] += marker_allele_cts[ii];
	    marker_allele_cts[ii] = marker_allele_cts[1];
	    marker_allele_cts[1] = 0;
	  }
	} else if (allele_min == marker_alleles[1]) {
	  if (ii > 1) {
	    marker_allele_cts[1] += marker_allele_cts[ii];
	    marker_allele_cts[ii] = marker_allele_cts[0];
	    marker_allele_cts[0] = 0;
	  }
	} else if (allele_min == marker_alleles[2]) {
	  if (ii == 3) {
            ii = marker_allele_cts[0];
	    marker_allele_cts[0] = marker_allele_cts[2];
	    marker_allele_cts[2] = ii;
	    ii = marker_allele_cts[1];
	    marker_allele_cts[1] = marker_allele_cts[3];
	    marker_allele_cts[3] = ii;
	  } else {
	    marker_allele_cts[ii] += marker_allele_cts[2];
	    marker_allele_cts[2] = marker_allele_cts[1 - ii];
	    marker_allele_cts[1 - ii] = 0;
	  }
	} else if (allele_min == marker_alleles[3]) {
	  if (ii == 2) {
            ii = marker_allele_cts[0];
	    marker_allele_cts[0] = marker_allele_cts[2];
	    marker_allele_cts[2] = ii;
	    ii = marker_allele_cts[1];
	    marker_allele_cts[1] = marker_allele_cts[3];
	    marker_allele_cts[3] = ii;
	  } else {
	    marker_allele_cts[ii] += marker_allele_cts[3];
	    marker_allele_cts[3] = marker_allele_cts[1 - ii];
	    marker_allele_cts[1 - ii] = 0;
	  }
	} else {
	  return -1;
	}
      }
      marker_alleles[0] = allele_maj;
      marker_alleles[1] = allele_min;
    }
  }
  return 0;
}

int read_external_freqs(char* freqname, FILE** freqfile_ptr, int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_exclude_ct, char* marker_ids, unsigned long max_marker_id_len, Chrom_info* chrom_info_ptr, char* marker_alleles, unsigned int* marker_allele_cts, double* set_allele_freqs, int binary_files, int maf_succ, char missing_geno, double exponent, int wt_needed, double* marker_weights) {
  unsigned int species = chrom_info_ptr->species;
  int freq_counts = 0;
  int duplicate_fail = 1;
  unsigned char* wkspace_mark;
  char* sorted_ids;
  int* id_map;
  int ii;
  int jj;
  char cc;
  char cc2;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  double maf;
  int c_hom_a1;
  int c_het;
  int c_hom_a2;
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
      bufptr = next_item(bufptr); // now at beginning of SNP name
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      read_next_terminate(bufptr, bufptr); // destructive read (\0 at end of item)
      ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
        if ((jj == get_marker_chrom(chrom_info_ptr, ii)) || (!jj) || (!get_marker_chrom(chrom_info_ptr, ii))) {
          cc = *bufptr2;
          bufptr2 = next_item(bufptr2);
	  if (cc == *bufptr2) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
          bufptr = next_item(bufptr2);
          if (no_more_items(bufptr)) {
            goto read_external_freqs_ret_INVALID_FORMAT;
          }
	  if (freq_counts) {
	    if (no_more_items(next_item(bufptr))) {
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
	  if (load_one_freq(cc, *bufptr2, maf, &(set_allele_freqs[ii]), binary_files? (&(marker_alleles[ii * 2])) : (&(marker_alleles[ii * 4])), binary_files? NULL : (&(marker_allele_cts[ii * 4])), missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
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
  } else if (!strcmp(tbuf, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(MISSING)\n")) {
    // known --freqx format, v0.13.5 or later
    while (fgets(tbuf, MAXLINELEN, *freqfile_ptr) != NULL) {
      jj = marker_code(species, tbuf);
      bufptr = next_item(tbuf); // now at beginning of SNP name
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      read_next_terminate(bufptr, bufptr); // destructive read (\0 at end of item)
      ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
        if ((jj == get_marker_chrom(chrom_info_ptr, ii)) || (!jj) || (!get_marker_chrom(chrom_info_ptr, ii))) {
          cc = *bufptr2;
          bufptr2 = next_item(bufptr2);
	  cc2 = *bufptr2;
	  if (cc == cc2) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  bufptr = next_item(bufptr2);
	  bufptr2 = next_item(bufptr);
	  bufptr3 = next_item(bufptr2);
	  if (no_more_items(bufptr3)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  c_hom_a1 = atoi(bufptr);
	  c_het = atoi(bufptr2);
	  c_hom_a2 = atoi(bufptr3);
	  maf = ((double)(c_hom_a1 * 2 + c_het + maf_succ)) / ((double)(2 * (c_hom_a1 + c_het + c_hom_a2 + maf_succ)));
	  if (load_one_freq(cc, cc2, maf, &(set_allele_freqs[ii]), binary_files? (&(marker_alleles[ii * 2])) : (&(marker_alleles[ii * 4])), binary_files? NULL : (&(marker_allele_cts[ii * 4])), missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
	  }
	  if (wt_needed) {
	    marker_weights[ii] = calc_wt_mean(exponent, c_het, c_hom_a1, c_hom_a2);
	  }
        }
      }
    }
    fputs(".frqx file loaded.\n", stdout);
  } else {
    // Also support GCTA-style frequency files:
    // [marker ID]\t[reference allele]\t[frequency of reference allele]\n
    do {
      if (is_eoln(*tbuf)) {
	continue;
      }
      bufptr = next_item(tbuf);
      if (!bufptr) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      read_next_terminate(tbuf, tbuf); // destructive read
      ii = bsearch_str(tbuf, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
        cc = *bufptr;
        bufptr = next_item(bufptr);
	if (no_more_items(bufptr)) {
          goto read_external_freqs_ret_INVALID_FORMAT;
	}
        if (sscanf(bufptr, "%lg", &maf) != 1) {
          goto read_external_freqs_ret_INVALID_FORMAT;
        }
	if (load_one_freq(missing_geno, cc, maf, &(set_allele_freqs[ii]), binary_files? (&(marker_alleles[ii * 2])) : (&(marker_alleles[ii * 4])), binary_files? NULL : (&(marker_allele_cts[ii * 4])), missing_geno)) {
	  goto read_external_freqs_ret_ALLELE_MISMATCH2;
	}
	if (wt_needed) {
	  marker_weights[ii] = calc_wt_mean_maf(exponent, set_allele_freqs[ii]);
	}
      } else {
	// if there aren't exactly 3 columns, this isn't a GCTA .freq file
	bufptr = next_item(bufptr);
	if (no_more_items(bufptr) || (!no_more_items(next_item(bufptr)))) {
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

int write_freqs(FILE** outfile_ptr, char* outname, unsigned int plink_maxsnp, int unfiltered_marker_ct, unsigned long* marker_exclude, double* set_allele_freqs, Chrom_info* chrom_info_ptr, char* marker_ids, unsigned long max_marker_id_len, char* marker_alleles, int* ll_cts, int* lh_cts, int* hh_cts, int freq_counts, int freqx, char missing_geno, unsigned long* marker_reverse) {
  int reverse = 0;
  unsigned int uii;
  char minor_allele;
  unsigned int chrom_idx;
  unsigned int chrom_end_marker;
  if (fopen_checked(outfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  if (freqx) {
    // MAF is not quite redundant with the hom/het counts, since it's affected
    // by --maf-succ.
    if (fputs("CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(MISSING)\n", *outfile_ptr) == EOF) {
      return RET_WRITE_FAIL;
    }
  } else if (plink_maxsnp < 5) {
    if (freq_counts) {
      if (fputs_checked(" CHR  SNP   A1   A2     C1     C2     G0\n", *outfile_ptr)) {
	return RET_WRITE_FAIL;
      }
      strcpy(tbuf, "%4d %4s    %c    %c %6d %6d %6d\n");
    } else {
      if (fputs_checked(" CHR  SNP   A1   A2          MAF  NCHROBS\n", *outfile_ptr)) {
        return RET_WRITE_FAIL;
      }
      strcpy(tbuf, "%4d %4s    %c    %c      %7.4g  %7d\n");
    }
  } else if (freq_counts) {
    sprintf(tbuf, " CHR %%%ds   A1   A2     C1     C2     G0\n", plink_maxsnp);
    if (fprintf(*outfile_ptr, tbuf, "SNP") < 0) {
      return RET_WRITE_FAIL;
    }
    sprintf(tbuf, "%%4d %%%ds    %%c    %%c %%6d %%6d %%6d\n", plink_maxsnp);
  } else {
    sprintf(tbuf, " CHR %%%ds   A1   A2          MAF  NCHROBS\n", plink_maxsnp);
    if (fprintf(*outfile_ptr, tbuf, "SNP") < 0) {
      return RET_WRITE_FAIL;
    }
    sprintf(tbuf, "%%4d %%%ds    %%c    %%c      %%7.4g  %%7d\n", plink_maxsnp);
  }
  for (chrom_idx = 0; chrom_idx < MAX_POSSIBLE_CHROM; chrom_idx++) {
    if (!chrom_exists(chrom_info_ptr, chrom_idx)) {
      continue;
    }
    chrom_end_marker = chrom_info_ptr->chrom_end[chrom_idx];
    for (uii = chrom_info_ptr->chrom_start[chrom_idx]; uii < chrom_end_marker; uii++) {
      if (is_set(marker_exclude, uii)) {
	continue;
      }
      reverse = is_set(marker_reverse, uii);
      minor_allele = marker_alleles[uii * 2 + 1 - (1 ^ reverse)];
      if (!minor_allele) {
	minor_allele = missing_geno;
      }
      if (freqx) {
	if (fprintf(*outfile_ptr, "%d\t%s\t%c\t%c\t%d\t%d\t%d\t%d\n", get_marker_chrom(chrom_info_ptr, uii), &(marker_ids[uii * max_marker_id_len]), minor_allele, marker_alleles[uii * 2 + (1 ^ reverse)], reverse? hh_cts[uii] : ll_cts[uii], lh_cts[uii], reverse? ll_cts[uii] : hh_cts[uii], g_indiv_ct - (ll_cts[uii] + lh_cts[uii] + hh_cts[uii])) < 0) {
	  return RET_WRITE_FAIL;
	}
      } else if (freq_counts) {
	if (fprintf(*outfile_ptr, tbuf, get_marker_chrom(chrom_info_ptr, uii), &(marker_ids[uii * max_marker_id_len]), marker_alleles[uii * 2], marker_alleles[uii * 2 + 1], 2 * ll_cts[uii] + lh_cts[uii], 2 * hh_cts[uii] + lh_cts[uii], g_indiv_ct - (ll_cts[uii] + lh_cts[uii] + hh_cts[uii])) < 0) {
	  return RET_WRITE_FAIL;
	}
      } else {
	if (fprintf(*outfile_ptr, tbuf, get_marker_chrom(chrom_info_ptr, uii), &(marker_ids[uii * max_marker_id_len]), minor_allele, marker_alleles[uii * 2 + (1 ^ reverse)], reverse? set_allele_freqs[uii] : (1.0 - set_allele_freqs[uii]), 2 * (ll_cts[uii] + lh_cts[uii] + hh_cts[uii])) < 0) {
	  return RET_WRITE_FAIL;
	}
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

unsigned int binary_geno_filter(double geno_thresh, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int* marker_exclude_ct_ptr, unsigned int indiv_ct, unsigned int* marker_allele_cts) {
  unsigned int orig_exclude_ct = *marker_exclude_ct_ptr;
  unsigned int geno_int_thresh = 2 * indiv_ct - (int)(geno_thresh * 2 * indiv_ct);
  unsigned int marker_uidx;
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    if ((marker_allele_cts[2 * marker_uidx] + marker_allele_cts[2 * marker_uidx + 1]) < geno_int_thresh) {
      set_bit(marker_exclude, marker_uidx, marker_exclude_ct_ptr);
    }
  }
  return (*marker_exclude_ct_ptr - orig_exclude_ct);
}

unsigned int text_geno_filter(double geno_thresh, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int* marker_exclude_ct_ptr, unsigned int indiv_ct, unsigned int* marker_allele_cts, unsigned int* missing_cts) {
  unsigned int orig_exclude_ct = *marker_exclude_ct_ptr;
  unsigned int geno_int_thresh = (int)(geno_thresh * 2 * indiv_ct + EPSILON);
  unsigned int marker_uidx;
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    if ((missing_cts[marker_uidx] * 2 + marker_allele_cts[4 * marker_uidx + 2] + marker_allele_cts[4 * marker_uidx + 3]) > geno_int_thresh) {
      set_bit(marker_exclude, marker_uidx, marker_exclude_ct_ptr);
    }
  }
  return (*marker_exclude_ct_ptr - orig_exclude_ct);
}

void text_normalize_marker_alleles(char* marker_alleles, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude) {
  char cc = marker_alleles[0];
  unsigned int marker_uidx;
  marker_alleles[0] = marker_alleles[1];
  marker_alleles[1] = cc;
  for (marker_uidx = 1; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    marker_alleles[2 * marker_uidx] = marker_alleles[4 * marker_uidx + 1];
    marker_alleles[2 * marker_uidx + 1] = marker_alleles[4 * marker_uidx];
  }
}

int text_load_hwe(FILE* pedfile, int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned long* founder_info, int nonfounders, char* marker_alleles, unsigned long long* line_mids, int pedbuflen, int hwe_all, char* pheno_c, int* hwe_lls, int* hwe_lhs, int* hwe_hhs, int* hwe_ll_allfs, int* hwe_lh_allfs, int* hwe_hh_allfs, int* ll_cts, int* lh_cts, int* hh_cts) {
  unsigned char* wkspace_mark;
  char* loadbuf;
  unsigned int indiv_uidx;
  char* cptr;
  unsigned int is_effective_founder;
  unsigned int always_count;
  int marker_uidx;
  char cc;
  char cc2;
  if (!pheno_c) {
    hwe_all = 1;
  }
  fill_int_zero(hwe_lls, unfiltered_marker_ct);
  fill_int_zero(hwe_lhs, unfiltered_marker_ct);
  fill_int_zero(hwe_hhs, unfiltered_marker_ct);
  fill_int_zero(hwe_ll_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_lh_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_hh_allfs, unfiltered_marker_ct);
  fill_int_zero(ll_cts, unfiltered_marker_ct);
  fill_int_zero(lh_cts, unfiltered_marker_ct);
  fill_int_zero(hh_cts, unfiltered_marker_ct);
  wkspace_mark = wkspace_base;
  if (wkspace_alloc_c_checked(&loadbuf, pedbuflen)) {
    return RET_NOMEM;
  }

  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    if (is_set(indiv_exclude, indiv_uidx)) {
      continue;
    }
    if (fseeko(pedfile, line_mids[indiv_uidx], SEEK_SET)) {
      return RET_READ_FAIL;
    }
    if (fgets((char*)loadbuf, pedbuflen, pedfile) == NULL) {
      return RET_READ_FAIL;
    }
    cptr = loadbuf;
    is_effective_founder = (nonfounders || is_set(founder_info, indiv_uidx));
    always_count = (is_effective_founder && (hwe_all || (!pheno_c[indiv_uidx])));
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      cc = *cptr;
      // no more need to check for invalid format, since this is at least
      // second pass
      cptr++;
      while ((*cptr == ' ') || (*cptr == '\t')) {
	cptr++;
      }
      if (((cc == marker_alleles[2 * marker_uidx]) || (cc == marker_alleles[2 * marker_uidx + 1])) && (!is_set(marker_exclude, marker_uidx))) {
	cc2 = *cptr;
	if (cc2 == marker_alleles[2 * marker_uidx]) {
	  if (cc == cc2) {
	    if (always_count) {
	      hwe_lls[marker_uidx] += 1;
	    } else if (is_effective_founder) {
	      hwe_ll_allfs[marker_uidx] += 1;
	    } else {
	      ll_cts[marker_uidx] += 1;
	    }
	  } else {
	    if (always_count) {
	      hwe_lhs[marker_uidx] += 1;
	    } else if (is_effective_founder) {
	      hwe_lh_allfs[marker_uidx] += 1;
	    } else {
	      lh_cts[marker_uidx] += 1;
	    }
	  }
	} else if (cc2 == marker_alleles[2 * marker_uidx + 1]) {
	  if (cc == cc2) {
	    if (always_count) {
	      hwe_hhs[marker_uidx] += 1;
	    } else if (is_effective_founder) {
	      hwe_hh_allfs[marker_uidx] += 1;
	    } else {
	      hh_cts[marker_uidx] += 1;
	    }
	  } else {
	    if (always_count) {
	      hwe_lhs[marker_uidx] += 1;
	    } else if (is_effective_founder) {
	      hwe_lh_allfs[marker_uidx] += 1;
	    } else {
	      lh_cts[marker_uidx] += 1;
	    }
	  }
	}
      }
      cptr++;
      while ((*cptr == ' ') || (*cptr == '\t')) {
	cptr++;
      }
    }
  }
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    hwe_ll_allfs[marker_uidx] += hwe_lls[marker_uidx];
    hwe_lh_allfs[marker_uidx] += hwe_lhs[marker_uidx];
    hwe_hh_allfs[marker_uidx] += hwe_hhs[marker_uidx];
    ll_cts[marker_uidx] += hwe_ll_allfs[marker_uidx];
    lh_cts[marker_uidx] += hwe_lh_allfs[marker_uidx];
    hh_cts[marker_uidx] += hwe_hh_allfs[marker_uidx];
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

void calc_marker_reverse_bin(unsigned long* marker_reverse, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_exclude_ct, double* set_allele_freqs) {
  unsigned int marker_uidx = 0;
  for (; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    if (set_allele_freqs[marker_uidx] < 0.5) {
      set_bit_noct(marker_reverse, marker_uidx);
    }
  }
}

void calc_marker_reverse_text(unsigned long* marker_reverse, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_exclude_ct, double* set_allele_freqs) {
  unsigned int marker_uidx = 0;
  for (; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    if (set_allele_freqs[marker_uidx] <= 0.5) {
      set_bit_noct(marker_reverse, marker_uidx);
    }
  }
}

void enforce_hwe_threshold(double hwe_thresh, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int* marker_exclude_ct_ptr, int* hwe_lls, int* hwe_lhs, int* hwe_hhs, int hwe_all, int* hwe_ll_allfs, int* hwe_lh_allfs, int* hwe_hh_allfs) {
  // N.B. requires het_probs to be allocated
  unsigned int removed_ct = 0;
  unsigned int marker_uidx = 0;
  hwe_thresh += EPSILON;
  for (; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    if (SNPHWE_t(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx], hwe_thresh)) {
      set_bit(marker_exclude, marker_uidx, marker_exclude_ct_ptr);
      removed_ct++;
    } else if (!hwe_all) {
      if (SNPHWE_t(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], hwe_thresh)) {
	set_bit(marker_exclude, marker_uidx, marker_exclude_ct_ptr);
	removed_ct++;
      }
    }
  }
  sprintf(logbuf, "%u SNP%s removed due to Hardy-Weinberg exact test (--hwe).\n", removed_ct, (removed_ct == 1)? "" : "s");
  logprintb();
}

void enforce_maf_threshold(double min_maf, double max_maf, int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int* marker_exclude_ct_ptr, double* set_allele_freqs) {
  unsigned int removed_ct = 0;
  int marker_uidx;
  double dxx;
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    dxx = get_maf(set_allele_freqs[marker_uidx]);
    if ((dxx < min_maf - EPSILON) || (dxx > max_maf + EPSILON)) {
      set_bit(marker_exclude, marker_uidx, marker_exclude_ct_ptr);
      removed_ct++;
    }
  }
  sprintf(logbuf, "%u SNP%s removed due to MAF threshold(s) (--maf/--max-maf).\n", removed_ct, (removed_ct == 1)? "" : "s");
  logprintb();
}

void calc_marker_weights(double exponent, int unfiltered_marker_ct, unsigned long* marker_exclude, int* ll_cts, int* lh_cts, int* hh_cts, double* marker_weights) {
  int marker_uidx;
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    if (marker_weights[marker_uidx] < 0.0) {
      marker_weights[marker_uidx] = calc_wt_mean(exponent, lh_cts[marker_uidx], ll_cts[marker_uidx], hh_cts[marker_uidx]);
    }
  }
}

void exclude_ambiguous_sex(unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int* indiv_exclude_ct_ptr, unsigned char* sex_info) {
  unsigned int indiv_uidx;
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    if (is_set(indiv_exclude, indiv_uidx)) {
      continue;
    }
    if (sex_info[indiv_uidx] == '\0') {
      set_bit(indiv_exclude, indiv_uidx, indiv_exclude_ct_ptr);
    }
  }
}

int load_ref_alleles(char* infilename, unsigned int unfiltered_marker_ct, char* marker_alleles, unsigned long* marker_reverse, char* marker_ids, unsigned long max_marker_id_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char* sorted_ids;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  int* id_map;
  int retval;
  int sorted_idx;
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
    if (is_eoln(*bufptr)) {
      continue;
    }
    bufptr2 = item_endnn2(bufptr);
    bufptr3 = skip_initial_spaces(bufptr2);
    if (is_eoln(*bufptr3)) {
      logprint("Error: Improperly formatted --reference-allele file.\n");
      goto load_ref_alleles_ret_INVALID_FORMAT;
    }
    *bufptr2 = '\0';
    sorted_idx = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - 1);
    if (sorted_idx == -1) {
      continue;
    } else {
      sorted_idx = id_map[sorted_idx];
      cc = *bufptr3;
      if (cc == marker_alleles[sorted_idx * 2]) {
	clear_bit_noct(marker_reverse, sorted_idx);
      } else if (cc == marker_alleles[sorted_idx * 2 + 1]) {
	set_bit_noct(marker_reverse, sorted_idx);
      } else {
        sprintf(logbuf, "Warning: Impossible reference allele assignment for marker %s.\n", bufptr);
	logprintb();
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

int write_snplist(FILE** outfile_ptr, char* outname, char* outname_end, int unfiltered_marker_ct, unsigned long* marker_exclude, char* marker_ids, unsigned long max_marker_id_len) {
  int marker_uidx;
  strcpy(outname_end, ".snplist");
  if (fopen_checked(outfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    if (fprintf(*outfile_ptr, "%s\n", &(marker_ids[marker_uidx * max_marker_id_len])) < 0) {
      return RET_WRITE_FAIL;
    }
  }
  if (fclose_null(outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  sprintf(logbuf, "Final list of SNP IDs written to %s.\n", outname);
  logprintb();
  return 0;
}

void normalize_phenos(double* new_phenos, unsigned int indiv_ct, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned char* sex_info, unsigned int sex_exclude) {
  // Assumes sex information is available for everyone.
  unsigned int indiv_uidx = 0;
  double pheno_tot = 0.0;
  double pheno_sq_tot = 0.0;
  unsigned int pheno_ct = 0;
  unsigned int indiv_idx;
  double dxx;
  double mean;
  double stdev_recip;

  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    if (sex_info[indiv_uidx] != sex_exclude) {
      dxx = new_phenos[indiv_idx];
      pheno_tot += dxx;
      pheno_sq_tot += dxx * dxx;
      pheno_ct++;
    }
    indiv_uidx++;
  }
  if (!pheno_ct) {
    return;
  }
  mean = pheno_tot / pheno_ct;
  stdev_recip = sqrt((double)(pheno_ct - 1) / (pheno_sq_tot - pheno_tot * mean));
  indiv_uidx = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    if (sex_info[indiv_uidx] != sex_exclude) {
      new_phenos[indiv_idx] = (new_phenos[indiv_idx] - mean) * stdev_recip;
    }
    indiv_uidx++;
  }
}

int calc_regress_pcs(char* evecname, int regress_pcs_normalize_pheno, int regress_pcs_sex_specific, int regress_pcs_clip, int max_pcs, FILE* pedfile, int bed_offset, unsigned int marker_ct, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, char* marker_ids, unsigned long max_marker_id_len, char* marker_alleles, Chrom_info* chrom_info_ptr, unsigned int* marker_pos, unsigned int indiv_ct, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, char* person_ids, unsigned int max_person_id_len, unsigned char* sex_info, double* pheno_d, double missing_phenod, FILE** outfile_ptr, char* outname, char* outname_end) {
  FILE* evecfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned int unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  int retval = 0;
  int pc_ct = 0;
  unsigned int pct = 1;
  int is_eigenvec = 0; // GCTA .eigenvec format instead of SMARTPCA .evec?
  int pc_ct_p1; // plus 1 to account for intercept
  double* pc_matrix;
  double* pc_orig_prod_sums; // pc_ct_p1 * pc_ct_p1, upper triangle filled
  double* pc_prod_sums; // (X'X)^{-1}
  double* x_prime_y; // X'Y
  double* beta_vec; // (X'X)^{-1}X'Y
  double* residual_vec;
  unsigned int marker_uidx;
  unsigned int marker_idx;
  unsigned int indiv_uidx;
  unsigned int indiv_idx;
  unsigned char* loadbuf;
  unsigned int* missing_cts;
  char* bufptr;
  char* id_buf;
  unsigned int uii;
  int ii;
  int jj;
  int cur_missing;
  char* person_id_ptr;
  MATRIX_INVERT_BUF1_TYPE* inv_1d_buf;
  double* dbl_2d_buf;
  double dxx;
  if (wkspace_alloc_uc_checked(&loadbuf, unfiltered_indiv_ct4)) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&missing_cts, indiv_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    return RET_NOMEM;
  }
  
  // try unaltered filename.  If that fails and the unaltered filename didn't
  // have an .evec or .eigenvec extension, then also try appending .evec and
  // appending .eigenvec.
  evecfile = fopen(evecname, "r");
  if (!evecfile) {
    ii = strlen(evecname);
    if (((ii >= 5) && (!memcmp(".evec", &(evecname[ii - 5]), 5))) || ((ii >= 9) && (!memcmp(".eigenvec", &(evecname[ii - 9]), 9)))) {
      return RET_OPEN_FAIL;
    }
    strcpy(&(evecname[ii]), ".evec");
    evecfile = fopen(evecname, "r");
    if (!evecfile) {
      strcpy(&(evecname[ii]), ".eigenvec");
      if (fopen_checked(&evecfile, evecname, "r")) {
	return RET_OPEN_FAIL;
      }
    }
  }

  tbuf[MAXLINELEN - 7] = ' ';
  tbuf[MAXLINELEN - 1] = ' ';
  if (!fgets(tbuf, MAXLINELEN - 6, evecfile)) {
    if (feof(evecfile)) {
      goto calc_regress_pcs_ret_INVALID_FORMAT;
    } else {
      goto calc_regress_pcs_ret_READ_FAIL;
    }
  }
  if (!tbuf[MAXLINELEN - 7]) {
    logprint("Error: Excessively long line in .evec/.eigenvec file.\n");
    goto calc_regress_pcs_ret_INVALID_FORMAT2;
  }
  bufptr = skip_initial_spaces(tbuf);
  if (no_more_items(bufptr)) {
    goto calc_regress_pcs_ret_INVALID_FORMAT;
  }
  if (memcmp(bufptr, "#eigvals:", 9)) {
    is_eigenvec = 1;
    bufptr = next_item(bufptr);
  }
  bufptr = next_item(bufptr);
  while ((!no_more_items(bufptr)) && ((*bufptr == '-') || ((*bufptr >= '0') && (*bufptr <= '9')))) {
    pc_ct++;
    bufptr = next_item(bufptr);
  }
  if (!pc_ct) {
    goto calc_regress_pcs_ret_INVALID_FORMAT;
  }
  if (pc_ct > max_pcs) {
    sprintf(logbuf, "%svec format detected.  Regressing on %d PC%s (out of %d).\n", is_eigenvec? "GCTA .eigen" : "SMARTPCA .e", max_pcs, (max_pcs == 1)? "" : "s", pc_ct);
    pc_ct = max_pcs;
  } else {
    sprintf(logbuf, "%svec format detected.  Regressing on %d principal component%s.\n", is_eigenvec? "GCTA .eigen" : "SMARTPCA .e", pc_ct, (pc_ct == 1)? "" : "s");
  }
  logprintb();
  pc_ct_p1 = pc_ct + 1;
  if (wkspace_alloc_d_checked(&pc_matrix, pc_ct_p1 * indiv_ct * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&pc_orig_prod_sums, pc_ct_p1 * pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&pc_prod_sums, pc_ct_p1 * pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&x_prime_y, pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&beta_vec, pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&residual_vec, pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  inv_1d_buf = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(pc_ct_p1 * sizeof(MATRIX_INVERT_BUF1_TYPE));
  if (!inv_1d_buf) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&dbl_2d_buf, pc_ct_p1 * pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }

  if (is_eigenvec) {
    indiv_idx = 0;
    while (1) {
      // todo: validate, and perhaps permute, family/indiv IDs
      bufptr = next_item_mult(skip_initial_spaces(tbuf), 2);
      for (ii = 0; ii < pc_ct; ii++) {
	if (no_more_items(bufptr)) {
	  goto calc_regress_pcs_ret_INVALID_FORMAT;
	}
	if (sscanf(bufptr, "%lg", &(pc_matrix[ii * indiv_ct + indiv_idx])) != 1) {
	  goto calc_regress_pcs_ret_INVALID_FORMAT;
	}
	bufptr = next_item(bufptr);
      }
      pc_matrix[pc_ct * indiv_ct + indiv_idx] = 1.0; // intercept
      if (++indiv_idx >= indiv_ct) {
	break;
      }
      if (!fgets(tbuf, MAXLINELEN, evecfile)) {
	if (feof(evecfile)) {
	  logprint("Error: Fewer individuals in .eigenvec file than expected.\n");
	  goto calc_regress_pcs_ret_INVALID_FORMAT2;
	} else {
	  goto calc_regress_pcs_ret_READ_FAIL;
	}
      }
    }
  } else {
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      if (!fgets(tbuf, MAXLINELEN, evecfile)) {
	if (feof(evecfile)) {
	  logprint("Error: Fewer individuals in .evec file than expected.\n");
	  goto calc_regress_pcs_ret_INVALID_FORMAT2;
	} else {
	  goto calc_regress_pcs_ret_READ_FAIL;
	}
      }
      bufptr = next_item(skip_initial_spaces(tbuf));
      for (ii = 0; ii < pc_ct; ii++) {
	if (no_more_items(bufptr)) {
	  goto calc_regress_pcs_ret_INVALID_FORMAT;
	}
	if (sscanf(bufptr, "%lg", &(pc_matrix[ii * indiv_ct + indiv_idx])) != 1) {
	  goto calc_regress_pcs_ret_INVALID_FORMAT;
	}
	bufptr = next_item(bufptr);
      }
      pc_matrix[pc_ct * indiv_ct + indiv_idx] = 1.0;
    }
  }
  if (fgets(tbuf, MAXLINELEN, evecfile)) {
    if (!no_more_items(skip_initial_spaces(tbuf))) {
      sprintf(logbuf, "Error: More individuals in .e%svec file than expected.\n", is_eigenvec? "igen" : "");
      logprintb();
      goto calc_regress_pcs_ret_INVALID_FORMAT2;
    }
  }
  fclose_null(&evecfile);

  // precalculate (X'X)
  fill_double_zero(pc_orig_prod_sums, pc_ct_p1 * pc_ct_p1);
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    for (ii = 0; ii < pc_ct_p1; ii++) {
      for (jj = ii; jj < pc_ct_p1; jj++) {
        pc_orig_prod_sums[ii * pc_ct_p1 + jj] += pc_matrix[ii * indiv_ct + indiv_idx] * pc_matrix[jj * indiv_ct + indiv_idx];
      }
    }
  }

  fill_uint_zero(missing_cts, indiv_ct);
  marker_uidx = 0;
  // .gen instead of .bgen because latter actually has lower precision(!) (15
  // bits instead of the ~20 you get from printf("%g", dxx)), and there's no
  // need for repeated random access.
  strcpy(outname_end, ".gen");
  if (fopen_checked(outfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  if (fseeko(pedfile, bed_offset, SEEK_SET)) {
    return RET_READ_FAIL;
  }
  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
      if (fseeko(pedfile, bed_offset + (unsigned long long)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	return RET_READ_FAIL;
      }
    }
    if (fread(loadbuf, 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
      return RET_READ_FAIL;
    }
    if (fprintf(*outfile_ptr, "%d %s %u %c %c", get_marker_chrom(chrom_info_ptr, marker_uidx), &(marker_ids[marker_uidx * max_marker_id_len]), marker_pos[marker_uidx], marker_alleles[2 * marker_uidx], marker_alleles[2 * marker_uidx + 1]) < 0) {
      return RET_WRITE_FAIL;
    }
    memcpy(pc_prod_sums, pc_orig_prod_sums, pc_ct_p1 * pc_ct_p1 * sizeof(double));
    fill_double_zero(x_prime_y, pc_ct_p1);
    indiv_uidx = 0;
    cur_missing = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
      uii = ((loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3);
      if (uii == 1) {
	cur_missing++;
	missing_cts[indiv_idx] += 1;
	for (ii = 0; ii < pc_ct_p1; ii++) {
	  for (jj = ii; jj < pc_ct_p1; jj++) {
	    pc_prod_sums[ii * pc_ct_p1 + jj] -= pc_matrix[ii * indiv_ct + indiv_idx] * pc_matrix[jj * indiv_ct + indiv_idx];
	  }
	}
      } else {
	uii = uii - (uii >> 1);
        for (ii = 0; ii < pc_ct_p1; ii++) {
	  x_prime_y[ii] += pc_matrix[ii * indiv_ct + indiv_idx] * uii;
	}
      }
      indiv_uidx++;
    }
    // last-minute update of lower triangle
    for (ii = 1; ii < pc_ct_p1; ii++) {
      for (jj = 0; jj < ii; jj++) {
	pc_prod_sums[ii * pc_ct_p1 + jj] = pc_prod_sums[jj * pc_ct_p1 + ii];
      }
    }
    invert_matrix(pc_ct_p1, pc_prod_sums, inv_1d_buf, dbl_2d_buf);
    for (ii = 0; ii < pc_ct_p1; ii++) {
      dxx = 0.0;
      for (jj = 0; jj < pc_ct_p1; jj++) {
        dxx += pc_prod_sums[ii * pc_ct_p1 + jj] * x_prime_y[jj];
      }
      beta_vec[ii] = dxx;
    }
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
      uii = ((loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3);
      if (uii == 1) {
	if (fwrite_checked(" 0 0 0", 6, *outfile_ptr)) {
	  return RET_WRITE_FAIL;
	}
      } else {
	uii = uii - (uii >> 1);
	dxx = 0.0;
        for (ii = 0; ii < pc_ct_p1; ii++) {
          dxx += pc_matrix[ii * indiv_ct + indiv_idx] * beta_vec[ii];
	}
	dxx = (double)uii - dxx;
	// now dxx is the residual, normally but not always in [0, 2]
	if (dxx < 1.0) {
	  if (dxx < 0.0) {
	    if (regress_pcs_clip) {
	      if (fwrite_checked(" 1 0 0", 6, *outfile_ptr)) {
		return RET_WRITE_FAIL;
	      }
	    } else {
	      if (fprintf(*outfile_ptr, " %g 0 %g", 1.0 - dxx * 0.5, dxx * 0.5) < 0) {
		return RET_WRITE_FAIL;
	      }
	    }
	  } else {
	    if (fprintf(*outfile_ptr, " %g %g 0", 1.0 - dxx, dxx) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  }
	} else {
	  if (dxx > 2.0) {
	    if (regress_pcs_clip) {
	      if (fwrite_checked(" 0 0 1", 6, *outfile_ptr)) {
		return RET_WRITE_FAIL;
	      }
	    } else {
	      if (fprintf(*outfile_ptr, " %g 0 %g", 1.0 - dxx * 0.5, dxx * 0.5) < 0) {
		return RET_WRITE_FAIL;
	      }
	    }
	  } else {
	    if (fprintf(*outfile_ptr, " 0 %g %g", 2.0 - dxx, dxx - 1.0) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  }
	}
      }
      indiv_uidx++;
    }
    if (fwrite_checked("\n", 1, *outfile_ptr)) {
      return RET_WRITE_FAIL;
    }
    if (marker_idx * 100LLU >= ((unsigned long long)pct * marker_ct)) {
      pct = ((unsigned long long)marker_idx * 100) / marker_ct;
      printf("\r%d%%", pct++);
      fflush(stdout);
    }
    marker_uidx++;
  }
  if (fclose_null(outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  strcpy(outname_end, ".sample");
  if (fopen_checked(outfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  if (fputs_checked("ID_1 ID_2 missing sex phenotype\n0 0 0 D P\n", *outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  // regress phenotype
  fill_double_zero(x_prime_y, pc_ct_p1);
  indiv_uidx = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    dxx = pheno_d[indiv_uidx];
    for (ii = 0; ii < pc_ct_p1; ii++) {
      x_prime_y[ii] += pc_matrix[ii * indiv_ct + indiv_idx] * dxx;
    }
    indiv_uidx++;
  }
  for (ii = 1; ii < pc_ct_p1; ii++) {
    for (jj = 0; jj < ii; jj++) {
      pc_orig_prod_sums[ii * pc_ct_p1 + jj] = pc_orig_prod_sums[jj * pc_ct_p1 + ii];
    }
  }
  invert_matrix(pc_ct_p1, pc_orig_prod_sums, inv_1d_buf, dbl_2d_buf);
  for (ii = 0; ii < pc_ct_p1; ii++) {
    dxx = 0.0;
    for (jj = 0; jj < pc_ct_p1; jj++) {
      dxx += pc_orig_prod_sums[ii * pc_ct_p1 + jj] * x_prime_y[jj];
    }
    beta_vec[ii] = dxx;
  }

  indiv_uidx = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    dxx = 0.0;
    for (ii = 0; ii < pc_ct_p1; ii++) {
      dxx += pc_matrix[ii * indiv_ct + indiv_idx] * beta_vec[ii];
    }
    residual_vec[indiv_idx] = pheno_d[indiv_uidx] - dxx;
    indiv_uidx++;
  }

  if (regress_pcs_normalize_pheno) {
    if (regress_pcs_sex_specific) {
      normalize_phenos(residual_vec, indiv_ct, unfiltered_indiv_ct, indiv_exclude, sex_info, 1);
      normalize_phenos(residual_vec, indiv_ct, unfiltered_indiv_ct, indiv_exclude, sex_info, 2);
    } else {
      normalize_phenos(residual_vec, indiv_ct, unfiltered_indiv_ct, indiv_exclude, sex_info, 3);
    }
  }

  indiv_uidx = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    person_id_ptr = &(person_ids[indiv_uidx * max_person_id_len]);
    uii = strlen_se(person_id_ptr);
    memcpy(id_buf, person_id_ptr, uii);
    id_buf[uii] = '\0';
    // todo: adjust pheno_d
    if (fprintf(*outfile_ptr, "%s %s %g %u %g\n", id_buf, next_item(person_id_ptr), (double)missing_cts[indiv_uidx] / (double)marker_ct, sex_info[indiv_uidx], residual_vec[indiv_idx]) < 0) {
      return RET_WRITE_FAIL;
    }
    indiv_uidx++;
  }
  if (fclose_null(outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  *outname_end = '\0';
  putchar('\r');
  sprintf(logbuf, "Principal component regression residuals and %sphenotype Z-scores %s%s.gen and %s.sample.\n", regress_pcs_sex_specific? "sex-specific " : "", regress_pcs_sex_specific? "\nwritten to " : "written to\n", outname, outname);
  logprintb();
  wkspace_reset(wkspace_mark);
  while (0) {
  calc_regress_pcs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_regress_pcs_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  calc_regress_pcs_ret_INVALID_FORMAT:
    logprint("Error: Improperly formatted .evec file.\n");
  calc_regress_pcs_ret_INVALID_FORMAT2:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(evecfile);
  return retval;
}

// implementation used in PLINK stats.cpp
double normdist(double zz) {
  double sqrt2pi = 2.50662827463;
  double t0;
  double z1;
  double p0;
  t0 = 1 / (1 + 0.2316419 * fabs(zz));
  z1 = exp(-0.5 * zz * zz) / sqrt2pi;
  p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
 return zz >= 0 ? 1 - p0 : p0;
}

int calc_genome(pthread_t* threads, FILE* pedfile, int bed_offset, unsigned int marker_ct, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, Chrom_info* chrom_info_ptr, unsigned int* marker_pos, double* set_allele_freqs, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, char* person_ids, unsigned int max_person_id_len, char* paternal_ids, unsigned int max_paternal_id_len, char* maternal_ids, unsigned int max_maternal_id_len, unsigned long* founder_info, int parallel_idx, int parallel_tot, char* outname, char* outname_end, int nonfounders, int calculation_type, int genome_output_gz, int genome_output_full, int genome_ibd_unbounded, int ppc_gap, Pedigree_rel_info pri) {
  FILE* outfile = NULL;
  gzFile gz_outfile = NULL;
  int retval = 0;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned int unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  unsigned char* loadbuf; // from file
  unsigned char* gptr;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* cptr5 = NULL;
  char* cptr6 = NULL;
  char* cptr7 = NULL;
  char* cptr8 = NULL;
  char* tbuf_mid = &(tbuf[64]);
  char* sptr_cur;
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int oo;
  int pp;
  int qq;
  unsigned long* glptr;
  unsigned long* glptr2;
  unsigned long* glptr3;
  unsigned int* giptr;
  unsigned int* giptr2;
  unsigned int* giptr3;
  unsigned int uii;
  unsigned int ujj;
  unsigned int ukk;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long ulkk;
  int missing_ct_buf[BITCT];
  int missing_ct_all;
  double set_allele_freq_buf[GENOME_MULTIPLEX];
  double e00 = 0.0;
  double e01 = 0.0;
  double e02 = 0.0;
  double e11 = 0.0;
  double e12 = 0.0;
  int ibd_prect = 0;
  double dpp;
  double dqq;
  double dpp_sq;
  double dqq_sq;
  double dxx;
  double dyy;
  double dxx_recip;
  double dyy_recip;
  double dxx1;
  double dxx2;
  double dyy1;
  double dyy2;
  long long cur_line = 0;
  long long tot_cells;
  long long tot_lines;
  // imitate PLINK behavior (see Plink::prettyPrintLengths() in helper.cpp), to
  // avoid randomly breaking existing scripts
  int max_person_fid_len = 4;
  int max_person_iid_len = 4;
  double num_alleles;
  double num_allelesf2;
  double num_allelesf3;
  int tstc;
  int is_founder_fixed = 0;
  unsigned int mp_lead_unfiltered_idx = 0;
  unsigned int mp_lead_idx = 0;
  unsigned int rel_space_id_fixed = 0;
  unsigned int family_id_fixed;
  unsigned int founder_ct = 0;
  long long llfct = 0;
  unsigned int indiv_uidx;
  unsigned int indiv_idx;
  unsigned int pct;
  char* fam1;
  char* fam2;

  while ((mp_lead_unfiltered_idx < unfiltered_marker_ct) && (is_set(marker_exclude, mp_lead_unfiltered_idx) || (!get_marker_chrom(chrom_info_ptr, mp_lead_unfiltered_idx)))) {
    mp_lead_unfiltered_idx++;
  }

  triangle_fill(g_thread_start, g_indiv_ct, g_thread_ct, parallel_tot - parallel_idx - 1, parallel_tot, 1, 1);
  // invert order, for --genome --parallel to naturally work
  for (uii = 0; uii <= g_thread_ct / 2; uii++) {
    jj = g_thread_start[uii];
    g_thread_start[uii] = g_indiv_ct - g_thread_start[g_thread_ct - uii];
    g_thread_start[g_thread_ct - uii] = g_indiv_ct - jj;
  }

  if (!parallel_idx) {
    cur_line = 1;
  }
  tstc = g_thread_start[g_thread_ct];
  // f(tstc) - f(g_thread_start[0])
  // f(0) = 0
  // f(1) = indiv_ct - 1
  // f(2) = 2indiv_ct - 3
  // ...
  // f(n) = nindiv_ct - n(n+1)/2
  tot_cells = (long long)g_indiv_ct * (tstc - g_thread_start[0]) - ((long long)tstc * (tstc + 1) - (long long)g_thread_start[0] * (g_thread_start[0] + 1)) / 2;
  tot_lines = cur_line + tot_cells;
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    if (!is_set(indiv_exclude, indiv_uidx)) {
      cptr = &(person_ids[indiv_uidx * max_person_id_len]);
      jj = strlen_se(cptr);
      cptr2 = next_item(cptr);
      if (jj > max_person_fid_len) {
	max_person_fid_len = jj + 2;
      }
      jj = strlen_se(cptr2);
      if (jj > max_person_iid_len) {
        max_person_iid_len = jj + 2;
      }
    }
  }
  if (wkspace_alloc_ui_checked(&g_missing_dbl_excluded, tot_cells * sizeof(int))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&g_indiv_missing_unwt, g_indiv_ct * sizeof(int))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&g_genome_main, tot_cells * 5 * sizeof(int))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_uc_checked(&loadbuf, GENOME_MULTIPLEX * unfiltered_indiv_ct4)) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_uc_checked(&g_geno, g_indiv_ct * (GENOME_MULTIPLEX / 4))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&g_masks, g_indiv_ct * (GENOME_MULTIPLEX / 4))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&g_mmasks, g_indiv_ct * sizeof(long))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_c_checked(&fam1, max_person_fid_len + 1)) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_c_checked(&fam2, max_person_fid_len + 1)) {
    goto calc_genome_ret_NOMEM;
  }

  fill_int_zero((int*)g_missing_dbl_excluded, tot_cells);
  fill_int_zero((int*)g_indiv_missing_unwt, g_indiv_ct);
  fill_int_zero((int*)g_genome_main, tot_cells * 5);
  if (!is_set(marker_exclude, 0)) {
    if (fseeko(pedfile, bed_offset, SEEK_SET)) {
      retval = RET_READ_FAIL;
      goto calc_genome_ret_1;
    }
  }
  ii = 0; // raw marker index
  g_low_ct = 0; // after excluding missing
  do {
    kk = marker_ct - g_low_ct;
    if (kk > GENOME_MULTIPLEX) {
      kk = GENOME_MULTIPLEX;
    }
    glptr2 = g_marker_window;
    for (ujj = 0; (int)ujj < kk; ujj++) {
      if (is_set(marker_exclude, ii)) {
	do {
	  ii++;
	} while (is_set(marker_exclude, ii));
	if (fseeko(pedfile, bed_offset + ii * unfiltered_indiv_ct4, SEEK_SET)) {
	  retval = RET_READ_FAIL;
	  goto calc_genome_ret_1;
	}
      }
      if (fread(&(loadbuf[ujj * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	retval = RET_READ_FAIL;
	goto calc_genome_ret_1;
      }
      set_allele_freq_buf[ujj] = set_allele_freqs[ii];
      // See comments in incr_genome(): the PPC test is time-critical and
      // we do a bit of unusual precomputation here to speed it up.
      //
      // Objective: Fill glptr[0] and glptr[1] with either
      // * a bitmask that excludes the correct number of SNPs, if the next
      //   eligible marker for the PPC test is within the same (BITCT / 2) SNP
      //   window, or
      // * twice the offset of the next SNP eligible for the PPC test, relative
      //   to the bottom of the currently loaded window,
      // because distinguishing between these two cases is effectively free.
      //
      // Then advance glptr two spaces.  The double storage eliminates a
      // divide-by-two in the inner loop at a low cost in cache space.
      if (mp_lead_unfiltered_idx < unfiltered_marker_ct) {
	if (get_marker_chrom(chrom_info_ptr, mp_lead_unfiltered_idx) == get_marker_chrom(chrom_info_ptr, g_low_ct + ujj)) {
	  uii = ppc_gap + marker_pos[g_low_ct + ujj];
	  if (marker_pos[mp_lead_unfiltered_idx] <= uii) {
	    ukk = get_chrom_end(chrom_info_ptr, mp_lead_unfiltered_idx);
	    do {
	      if (!is_set(marker_exclude, mp_lead_unfiltered_idx)) {
		mp_lead_idx++;
	      }
	      mp_lead_unfiltered_idx++;
	    } while ((mp_lead_unfiltered_idx < unfiltered_marker_ct) && (is_set(marker_exclude, mp_lead_unfiltered_idx) || ((mp_lead_unfiltered_idx < ukk) && (marker_pos[mp_lead_unfiltered_idx] <= uii))));
	  }
	}
      }
      if (mp_lead_unfiltered_idx < unfiltered_marker_ct) {
	ulii = 2 * (mp_lead_unfiltered_idx - g_low_ct);
	if (ulii < BITCT + (2 * (ujj & (~(BITCT2 - 1))))) {
	  ulii = ~0LU << (ulii & (BITCT - 1));
	}
      } else {
	ulii = 2 * (unfiltered_marker_ct + GENOME_MULTIPLEX);
      }

      *glptr2++ = ulii;
      *glptr2++ = ulii;
      ii++;
    }
    if (kk < GENOME_MULTIPLEX) {
      memset(&(loadbuf[kk * unfiltered_indiv_ct4]), 0, (GENOME_MULTIPLEX - kk) * unfiltered_indiv_ct4);
      fill_long_zero((long*)g_geno, g_indiv_ct * (GENOME_MULTIPLEX / BITCT2));
      fill_ulong_zero(g_masks, g_indiv_ct * (GENOME_MULTIPLEX / BITCT2));
      for (mm = kk * 2; mm < GENOME_MULTIPLEX2; mm++) {
	*glptr2++ = GENOME_MULTIPLEX2;
      }
    }
    g_high_ct = g_low_ct + kk;
    for (jj = 0; jj < kk; jj += BITCT) {
      glptr = &(((unsigned long*)g_geno)[jj / BITCT2]);
      glptr2 = &(g_masks[jj / BITCT2]);
      glptr3 = g_mmasks;
      giptr = g_indiv_missing_unwt;
      indiv_uidx = 0;
      fill_int_zero(missing_ct_buf, BITCT);
      missing_ct_all = 0;
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	oo = (indiv_uidx % 4) * 2;
	ulii = 0;
	ulkk = 0;
        gptr = &(loadbuf[indiv_uidx / 4 + jj * unfiltered_indiv_ct4]);
	qq = (nonfounders || is_founder(founder_info, indiv_uidx));
	if (qq) {
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      missing_ct_buf[pp] += 1;
	      ulkk |= 1LU << pp;
	      *giptr += 1;
	    }
	  }
	} else {
	  missing_ct_all++;
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      ulkk |= 1LU << pp;
	      *giptr += 1;
	    }
	  }
	}
	ulii ^= FIVEMASK;
	*glptr++ = ulii;
	ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	*glptr2++ = ulii * 3;
	*glptr3 = ulkk;

	ulii = 0;
	ulkk = 0;
	gptr = &(loadbuf[indiv_uidx / 4 + (jj + BITCT2) * unfiltered_indiv_ct4]);
	if (qq) {
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      missing_ct_buf[pp + BITCT2] += 1;
	      ulkk |= 1LU << pp;
	      *giptr += 1;
	    }
	  }
	} else {
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      ulkk |= 1LU << pp;
	      *giptr += 1;
	    }
	  }
	}
	ulii ^= FIVEMASK;
	*glptr = ulii;
	ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	*glptr2 = ulii * 3;
	*glptr3++ |= ulkk << BITCT2;
	glptr = &(glptr[(GENOME_MULTIPLEX2 / BITCT) - 1]);
	glptr2 = &(glptr2[(GENOME_MULTIPLEX2 / BITCT) - 1]);
	giptr++;
	indiv_uidx++;
      }
      nn = kk - jj;
      if (nn > BITCT) {
	nn = BITCT;
      }
      for (mm = 0; mm < nn; mm++) {
	dpp = set_allele_freq_buf[jj + mm];
	dqq = 1.0 - dpp;
	num_alleles = (double)(2 * (g_indiv_ct - missing_ct_buf[mm] - missing_ct_all));
	if ((num_alleles > 3) && (dpp > 0.0) && (dqq > 0.0)) {
	  // update e00, e01, e02, e11, e12, ibd_prect
	  // see Plink::preCalcGenomeIBD() in genome.cpp
	  num_allelesf2 = num_alleles * num_alleles / ((num_alleles - 1) * (num_alleles - 2));
	  num_allelesf3 = num_allelesf2 * num_alleles / (num_alleles - 3);
	  dxx = dpp * num_alleles;
	  dyy = dqq * num_alleles;
	  dxx_recip = 1.0 / dxx;
	  dyy_recip = 1.0 / dyy;
	  dpp_sq = dpp * dpp;
	  dqq_sq = dqq * dqq;
	  dxx1 = (dxx - 1) * dxx_recip;
	  dxx2 = dxx1 * (dxx - 2) * dxx_recip;
	  dyy1 = (dyy - 1) * dyy_recip;
	  dyy2 = dyy1 * (dyy - 2) * dyy_recip;
	  e00 += 2 * dpp_sq * dqq_sq * dxx1 * dyy1 * num_allelesf3;
	  e01 += 4 * dpp * dqq * num_allelesf3 * (dpp_sq * dxx2 + dqq_sq * dyy2);
	  e02 += num_allelesf3 * (dqq_sq * dqq_sq * dyy2 * (dyy - 3) * dyy_recip + dpp_sq * dpp_sq * dxx2 * (dxx - 3) * dxx_recip + 4 * dpp_sq * dqq_sq * dxx1 * dyy1);
	  e11 += 2 * dpp * dqq * num_allelesf2 * (dpp * dxx1 + dqq * dyy1);
	  e12 += num_allelesf2 * (dpp_sq * dpp * dxx2 + dqq_sq * dqq * dyy2 + dpp_sq * dqq * dxx1 + dpp * dqq_sq * dyy1);
	  ibd_prect++;
	}
      }
      for (ulii = 1; ulii < g_thread_ct; ulii++) {
	if (pthread_create(&(threads[ulii - 1]), NULL, &calc_genomem_thread, (void*)ulii)) {
	  goto calc_genome_ret_THREAD_CREATE_FAIL;
	}
      }
      incr_dists_rm_inv(g_missing_dbl_excluded, 0);
      for (ukk = 0; ukk < g_thread_ct - 1; ukk++) {
	pthread_join(threads[ukk], NULL);
      }
    }

    for (ulii = 1; ulii < g_thread_ct; ulii++) {
      if (pthread_create(&(threads[ulii - 1]), NULL, &calc_genome_thread, (void*)ulii)) {
	goto calc_genome_ret_THREAD_CREATE_FAIL;
      }
    }
    incr_genome(g_genome_main, (unsigned long*)g_geno, 0);
    for (ukk = 0; ukk < g_thread_ct - 1; ukk++) {
      pthread_join(threads[ukk], NULL);
    }
    g_low_ct = g_high_ct;
    printf("\r%d markers complete.", g_low_ct);
    fflush(stdout);
  } while (g_low_ct < marker_ct);
  fputs("\rIBD calculations complete.  \n", stdout);
  logstr("IBD calculations complete.\n");
  dxx = 1.0 / (double)ibd_prect;
  e00 *= dxx;
  e01 *= dxx;
  e02 *= dxx;
  e11 *= dxx;
  e12 *= dxx;

  if (calculation_type & CALC_PLINK_IBS_MATRIX) {
    strcpy(outname_end, ".mibs");
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_genome_ret_OPEN_FAIL;
    }
    giptr = g_genome_main;
    giptr2 = g_missing_dbl_excluded;
    pct = 1;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      giptr3 = g_indiv_missing_unwt;
      uii = marker_ct - giptr3[indiv_idx];
      uljj = (int)indiv_idx - 1; // not referenced when indiv_idx == 0
      for (ulii = 0; ulii < indiv_idx; ulii++) {
	if (fprintf(outfile, "%g ", 1.0 - ((double)(g_genome_main[uljj * 5] + 2 * g_genome_main[uljj * 5 + 1])) / ((double)(2 * (uii - (*giptr3++) + g_missing_dbl_excluded[uljj])))) < 0) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
	uljj += g_indiv_ct - ulii - 2;
      }
      if (fwrite_checked("1 ", 2, outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      giptr3++;
      for (ujj = indiv_idx + 1; ujj < g_indiv_ct; ujj++) {
	if (fprintf(outfile, "%g ", 1.0 - ((double)((*giptr) + 2 * giptr[1])) / ((double)(2 * (uii - (*giptr3++) + (*giptr2++))))) < 0) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
	giptr = &(giptr[5]);
      }
      if (fwrite_checked("\n", 1, outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      if (indiv_idx * 100 >= (pct * g_indiv_ct)) {
        pct = (indiv_idx * 100) / g_indiv_ct;
	printf("\rWriting... %d%%", pct++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_genome_ret_WRITE_FAIL;
    }
    putchar('\r');
    sprintf(logbuf, "IBS matrix written to %s.\n", outname);
    logprintb();
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      goto calc_genome_ret_1;
    }
  }

  if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
    strcpy(outname_end, ".mdist");
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_genome_ret_OPEN_FAIL;
    }
    giptr = g_genome_main;
    giptr2 = g_missing_dbl_excluded;
    kk = 1;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      giptr3 = g_indiv_missing_unwt;
      uii = marker_ct - giptr3[indiv_idx];
      uljj = indiv_idx - 1;
      for (ulii = 0; ulii < indiv_idx; ulii++) {
	if (fprintf(outfile, "%g ", ((double)(g_genome_main[uljj * 5] + 2 * g_genome_main[uljj * 5 + 1])) / ((double)(2 * (uii - (*giptr3++) + g_missing_dbl_excluded[uljj])))) < 0) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
	uljj += g_indiv_ct - ulii - 2;
      }
      if (fwrite_checked("0 ", 2, outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      giptr3++;
      for (ujj = indiv_idx + 1; ujj < g_indiv_ct; ujj++) {
	if (fprintf(outfile, "%g ", ((double)((*giptr) + 2 * giptr[1])) / ((double)(2 * (uii - (*giptr3++) + (*giptr2++))))) < 0) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
	giptr = &(giptr[5]);
      }
      if (fwrite_checked("\n", 1, outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      if (indiv_idx * 100 >= (kk * g_indiv_ct)) {
	kk = (indiv_idx * 100) / g_indiv_ct;
	printf("\rWriting... %d%%", kk++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_genome_ret_WRITE_FAIL;
    }
    putchar('\r');
    sprintf(logbuf, "Distances (proportions) written to %s.\n", outname);
    logprintb();
    strcpy(outname_end, ".mdist.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      goto calc_genome_ret_1;
    }
  }

  sprintf(tbuf, " %%%ds %%%ds %%%ds %%%ds ", max_person_fid_len - 1, max_person_iid_len - 1, max_person_fid_len - 1, max_person_iid_len - 1);
  if (!parallel_idx) {
    if (genome_output_full) {
      sprintf(tbuf_mid, "%%%ds%%%ds%%%ds%%%ds RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO    IBS0    IBS1    IBS2  HOMHOM  HETHET\n", max_person_fid_len, max_person_iid_len, max_person_fid_len, max_person_iid_len);
    } else {
      sprintf(tbuf_mid, "%%%ds%%%ds%%%ds%%%ds RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO\n", max_person_fid_len, max_person_iid_len, max_person_fid_len, max_person_iid_len);
    }
  }
  mm = 1;
  ulii = 0;
  uljj = 0;
  if (genome_output_gz) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".genome.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".genome.gz");
    }
    if (gzopen_checked(&gz_outfile, outname, "wb")) {
      goto calc_genome_ret_1;
    }
    if (!parallel_idx) {
      if (!gzprintf(gz_outfile, tbuf_mid, " FID1", " IID1", " FID2", " IID2")) {
	goto calc_genome_ret_WRITE_FAIL;
      }
    }
  } else {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".genome.%d", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".genome");
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_genome_ret_OPEN_FAIL;
    }
    if (!parallel_idx) {
      if (fprintf(outfile, tbuf_mid, " FID1", " IID1", " FID2", " IID2") < 0) {
	goto calc_genome_ret_WRITE_FAIL;
      }
    }
  }
  for (ii = g_thread_start[0]; ii < tstc; ii++) {
    cptr = &(person_ids[ii * max_person_id_len]);
    jj = strlen_se(cptr);
    memcpy(fam1, cptr, jj);
    fam1[jj] = '\0';
    cptr2 = next_item(cptr);
    if (paternal_ids) {
      cptr5 = &(paternal_ids[ii * max_paternal_id_len]);
      cptr6 = &(maternal_ids[ii * max_maternal_id_len]);
      is_founder_fixed = is_founder(founder_info, ii);
      rel_space_id_fixed = pri.family_rel_nf_idxs[ii];
      family_id_fixed = pri.family_idxs[ii];
      founder_ct = pri.family_founder_cts[family_id_fixed];
      llfct = (long long)founder_ct * (founder_ct - 1);
    }
    for (ujj = ii + 1; ujj < g_indiv_ct; ujj++) {
      cptr3 = &(person_ids[ujj * max_person_id_len]);
      jj = strlen_se(cptr3);
      memcpy(fam2, cptr3, jj);
      fam2[jj] = '\0';
      cptr4 = next_item(cptr3);
      if (paternal_ids) {
	cptr7 = &(paternal_ids[ujj * max_paternal_id_len]);
	cptr8 = &(maternal_ids[ujj * max_maternal_id_len]);
      }
      sptr_cur = &(tbuf_mid[sprintf(tbuf_mid, tbuf, fam1, cptr2, fam2, cptr4)]);
      if (!strcmp(fam1, fam2)) {
	while (1) {
	  if (paternal_ids) {
	    if (!(is_founder_fixed || is_set(founder_info, ujj))) {
	      if ((!strcmp(cptr5, cptr7)) && (!strcmp(cptr6, cptr8))) {
		sptr_cur += sprintf(sptr_cur, "FS ");
		break;
	      } else if ((!strcmp(cptr5, cptr7)) || (!strcmp(cptr6, cptr8))) {
		sptr_cur += sprintf(sptr_cur, "HS ");
		break;
	      }
	    }
	    if ((!strcmp(cptr5, cptr4)) || (!strcmp(cptr6, cptr4)) || (!strcmp(cptr7, cptr2)) || (!strcmp(cptr8, cptr2))) {
	      sptr_cur += sprintf(sptr_cur, "PO ");
	      break;
	    }
	  }
	  sptr_cur += sprintf(sptr_cur, "OT ");
	  break;
	}
	// insert relationship
	if (!paternal_ids) {
	  sptr_cur += sprintf(sptr_cur, "    0");
	} else {
	  oo = is_set(founder_info, ujj);
	  if (is_founder_fixed && oo) {
	    dxx = 0.0;
	  } else {
	    nn = pri.family_rel_nf_idxs[ujj];
            if (is_founder_fixed || ((nn > (int)rel_space_id_fixed) && (!oo))) {
	      dxx = pri.rel_space[rel_space_id_fixed + ((long long)nn * (nn - 1) - llfct) / 2];
	    } else {
	      dxx = pri.rel_space[nn + ((long long)rel_space_id_fixed * (rel_space_id_fixed - 1) - llfct) / 2];
	    }
	  }
	  sptr_cur += sprintf(sptr_cur, "%5g", dxx);
	}
      } else {
	sptr_cur += sprintf(sptr_cur, "UN    NA");
      }
      nn = marker_ct - g_indiv_missing_unwt[ii] - g_indiv_missing_unwt[ujj] + g_missing_dbl_excluded[uljj];
      oo = nn - g_genome_main[ulii] - g_genome_main[ulii + 1];
      dxx = (double)g_genome_main[ulii + 1] / (e00 * nn);
      dyy = ((double)g_genome_main[ulii] - dxx * e01 * nn) / (e11 * nn);
      dxx1 = ((double)oo - nn * (dxx * e02 + dyy * e12)) / ((double)nn);
      if (!genome_ibd_unbounded) {
	if (dxx > 1) {
	  dxx = 1;
	  dyy = 0;
	  dxx1 = 0;
	} else if (dyy > 1) {
	  dyy = 1;
	  dxx = 0;
	  dxx1 = 0;
	} else if (dxx1 > 1) {
	  dxx1 = 1;
	  dyy = 0;
	  dxx1 = 0;
	} else if (dxx < 0) {
	  dxx2 = 1.0 / (dyy + dxx1);
	  dyy *= dxx2;
	  dxx1 *= dxx2;
	  dxx = 0;
	}
	if (dyy < 0) {
	  dxx2 = 1.0 / (dxx + dxx1);
	  dxx *= dxx2;
	  dxx1 *= dxx2;
	  dyy = 0;
	}
	if (dxx1 < 0) {
	  dxx2 = 1.0 / (dxx + dyy);
	  dxx *= dxx2;
	  dyy *= dxx2;
	  dxx1 = 0;
	}
      }
      sptr_cur += sprintf(sptr_cur, " % .4f % .4f % .4f % .4f  ", dxx, dyy, dxx1, dyy * 0.5 + dxx1);

      if (g_pheno_c) {
	if ((g_pheno_c[ii] != 1) && (g_pheno_c[ujj] != 1)) {
	  memcpy(sptr_cur, "-1", 2);
	} else if ((g_pheno_c[ii] == 1) && (g_pheno_c[ujj] == 1)) {
	  memcpy(sptr_cur, " 1", 2);
	} else {
	  memcpy(sptr_cur, " 0", 2);
	}
      } else {
	memcpy(sptr_cur, "NA", 2);
      }
      sptr_cur += 2;
      dxx = (double)g_genome_main[ulii + 4];
      dyy = (double)g_genome_main[ulii + 3];
      dxx1 = 1.0 / ((double)(g_genome_main[ulii + 4] + g_genome_main[ulii + 3]));
      dxx2 = normdist((dxx * dxx1 - 0.666666) / (sqrt(0.2222222 * dxx1)));
      sptr_cur += sprintf(sptr_cur, "  %.6f  %.4f ", 1.0 - (g_genome_main[ulii] + 2 * g_genome_main[ulii + 1]) / ((double)(2 * nn)), dxx2);
      if (g_genome_main[ulii + 3]) {
	dxx1 = dxx / dyy;
	sptr_cur += sprintf(sptr_cur, "%7.4f", dxx1);
      } else {
	memcpy(sptr_cur, "     NA", 7);
	sptr_cur += 7;
      }
      if (genome_output_full) {
	sptr_cur += sprintf(sptr_cur, " %7d %7d %7d %s%.4f %s%.4f\n", g_genome_main[ulii + 1], g_genome_main[ulii], oo, (dyy < 9.5)? " " : "", dyy, (dxx < 9.5)? " " : "", dxx);
      } else {
	*sptr_cur++ = '\n';
      }
      if (genome_output_gz) {
	if (gzwrite_checked(gz_outfile, tbuf_mid, (unsigned long)(sptr_cur - tbuf_mid))) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
      } else {
	if (fwrite_checked(tbuf_mid, (unsigned long)(sptr_cur - tbuf_mid), outfile)) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
      }
      ulii += 5;
      uljj++;
    }
    cur_line += g_indiv_ct - ii - 1;
    if (cur_line * 100 >= tot_lines * mm) {
      mm = (cur_line * 100) / tot_lines;
      printf("\rWriting... %d%%", mm++);
      fflush(stdout);
    }
  }
  putchar('\r');
  sprintf(logbuf, "Finished writing %s.\n", outname);
  logprintb();
  wkspace_reset(wkspace_mark);
  while (0) {
  calc_genome_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_genome_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_genome_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  calc_genome_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    for (uljj = 0; uljj < ulii - 1; uljj++) {
      pthread_join(threads[uljj], NULL);
    }
    break;
  }
 calc_genome_ret_1:
  gzclose_cond(gz_outfile);
  fclose_cond(outfile);
  return retval;
}

int ld_process_load(unsigned char* loadbuf, unsigned long* geno_buf, unsigned long* mask_buf, unsigned long* missing_buf, double* marker_stdev_ptr, int unfiltered_indiv_ct, unsigned long* indiv_exclude, int indiv_ct, int indiv_ctbit, int indiv_trail_ct) {
  int unfiltered_idx = 0;
  int write_offset;
  int sloop_max;
  int write_idx;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long new_geno;
  unsigned long new_mask;
  unsigned long new_missing;
  int missing_ct = 0;
  int sq_sum = 0;
  int sum = -indiv_ct;
  double non_missing_recip;
  for (write_offset = 0; write_offset < indiv_ctbit * BITCT; write_offset += BITCT) {
    sloop_max = indiv_ct - write_offset;
    if (sloop_max > BITCT2) {
      sloop_max = BITCT2;
    }
    new_geno = 0;
    new_mask = 0;
    new_missing = 0;
    for (write_idx = 0; write_idx < sloop_max; write_idx++) {
      unfiltered_idx = next_non_set_unsafe(indiv_exclude, unfiltered_idx);
      ulii = (loadbuf[unfiltered_idx / 4] >> ((unfiltered_idx % 4) * 2)) & 3;

      // Nothing time-critical here, but may as well do it branchlessly.
      // Desired encodings:
      // new_geno: nonset homozygote -> 00
      //           het/missing       -> 01
      //           set homozygote    -> 10
      // Given PLINK encoding xx, this is (xx - (xx >> 1)).
      //
      // new_missing: missing   -> 1
      //              otherwise -> 0
      // xx & ((xx ^ 2) >> 1) is one way to do this.
      //
      // new_mask: missing   -> 00
      //           otherwise -> 11
      // (new_missing ^ 1) * 3 works.
      uljj = ulii >> 1;
      new_geno |= (ulii - uljj) << (write_idx * 2);
      uljj = ulii & (uljj ^ 1);
      new_missing |= uljj << write_idx;
      new_mask |= ((uljj ^ 1) * 3) << (write_idx * 2);
      unfiltered_idx++;
    }
    *geno_buf++ = new_geno;
    *mask_buf++ = new_mask;
    // new_mask needed for proper post-ending handling
    sq_sum += popcount2_long((new_geno ^ FIVEMASK) & FIVEMASK & new_mask);
    sum += popcount2_long(new_geno);
    new_geno = 0;
    new_mask = 0;
    if (write_offset + BITCT2 <= indiv_ct) {
      // +0 hom1, +1 het or missing, +2 hom2
      sloop_max = indiv_ct - write_offset - BITCT2;
      if (sloop_max > BITCT2) {
	sloop_max = BITCT2;
      }
      for (write_idx = 0; write_idx < sloop_max; write_idx++) {
	unfiltered_idx = next_non_set_unsafe(indiv_exclude, unfiltered_idx);
	ulii = (loadbuf[unfiltered_idx / 4] >> ((unfiltered_idx % 4) * 2)) & 3;
	uljj = ulii >> 1;
	new_geno |= (ulii - uljj) << (write_idx * 2);
	uljj = ulii & (uljj ^ 1);
	new_missing |= uljj << (write_idx + BITCT2);
	new_mask |= ((uljj ^ 1) * 3) << (write_idx * 2);
	unfiltered_idx++;
      }
    }
    *geno_buf++ = new_geno;
    *mask_buf++ = new_mask;
    sq_sum += popcount2_long((new_geno ^ FIVEMASK) & FIVEMASK & new_mask);
    sum += popcount2_long(new_geno);
    *missing_buf++ = new_missing;
    missing_ct += popcount_long(new_missing);
  }
  fill_ulong_zero(geno_buf, indiv_trail_ct);
  fill_ulong_zero(mask_buf, indiv_trail_ct);
  non_missing_recip = 1.0 / (indiv_ct - missing_ct);
  *marker_stdev_ptr = sqrt((sq_sum - (non_missing_recip * sum) * sum) * non_missing_recip);
  return missing_ct;
}

unsigned int sparse_intersection_ct(unsigned long* sparse_buf1, unsigned long* sparse_buf2, int len) {
  int ii;
  unsigned int ct = 0;
  unsigned long ulii;
  for (ii = 0; ii < len; ii++) {
    ulii = (*sparse_buf1++) & (*sparse_buf2++);
    while (ulii) {
      ulii &= ulii - 1;
      ct++;
    }
  }
  return ct;
}

unsigned int ld_prune_next_valid_chrom_start(unsigned long* marker_exclude, int cur_idx, Chrom_info* chrom_info_ptr, int unfiltered_marker_ct) {
  cur_idx = next_non_set_unsafe(marker_exclude, cur_idx);
  if (!nz_chrom(chrom_info_ptr, cur_idx)) {
    cur_idx = chrom_info_ptr->chrom_end[0];
    while (is_set(marker_exclude, cur_idx)) {
      if (cur_idx == unfiltered_marker_ct) {
	return cur_idx;
      }
      cur_idx++;
    }
  }
  return cur_idx;
}

void ld_prune_start_chrom(unsigned int ld_window_kb, unsigned int* cur_chrom_ptr, unsigned int* chrom_end_ptr, unsigned int window_unfiltered_start, unsigned int* live_indices, unsigned int* start_arr, unsigned int* window_unfiltered_end_ptr, unsigned int ld_window_size, int* cur_window_size_ptr, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, Chrom_info* chrom_info_ptr, unsigned int* marker_pos) {
  unsigned int cur_chrom = get_marker_chrom(chrom_info_ptr, window_unfiltered_start);
  unsigned int window_unfiltered_end = window_unfiltered_start + 1;
  unsigned int chrom_end = chrom_info_ptr->chrom_end[cur_chrom];
  unsigned int uii = 0;
  unsigned int window_size;
  live_indices[0] = window_unfiltered_start;
  if (ld_window_kb) {
    window_size = 0;
    while ((window_unfiltered_start + window_size < chrom_end) && (marker_pos[window_unfiltered_start + window_size] <= marker_pos[window_unfiltered_start] + (1000 * ld_window_size))) {
      window_size++;
    }
  } else {
    window_size = ld_window_size;
  }
  for (uii = 1; uii < window_size; uii++) {
    while (is_set(marker_exclude, window_unfiltered_end)) {
      window_unfiltered_end++;
      if (window_unfiltered_end == chrom_end) {
	break;
      }
    }
    if (window_unfiltered_end == chrom_end) {
      break;
    }
    start_arr[uii - 1] = window_unfiltered_end;
    live_indices[uii] = window_unfiltered_end;
    window_unfiltered_end++;
  }
  *cur_window_size_ptr = (int)uii;
  start_arr[uii - 1] = window_unfiltered_end;
  *cur_chrom_ptr = cur_chrom;
  *chrom_end_ptr = chrom_end;
  *window_unfiltered_end_ptr = window_unfiltered_end;
}

int ld_prune(FILE* bedfile, int bed_offset, unsigned int marker_ct, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, char* marker_ids, unsigned long max_marker_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, unsigned int* marker_pos, int unfiltered_indiv_ct, unsigned long* indiv_exclude, int ld_window_size, int ld_window_kb, int ld_window_incr, double ld_last_param, char* outname, char* outname_end, int calculation_type) {
  // todo: replace is_set with founder-sensitive check
  // for future consideration: chromosome-based multithread/parallel?
  FILE* outfile_in = NULL;
  FILE* outfile_out = NULL;
  unsigned int unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  unsigned int unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  unsigned int indiv_ctbit = (g_indiv_ct + BITCT - 1) / BITCT;
  unsigned int indiv_ct_mld = (g_indiv_ct + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
  int indiv_ct_mld_m1 = (int)indiv_ct_mld - 1;
#if __LP64__
  int indiv_ct_mld_rem = (MULTIPLEX_LD / 192) - (indiv_ct_mld * MULTIPLEX_LD - g_indiv_ct) / 192;
#else
  int indiv_ct_mld_rem = (MULTIPLEX_LD / 48) - (indiv_ct_mld * MULTIPLEX_LD - g_indiv_ct) / 48;
#endif
  int indiv_ct_mld_long = indiv_ct_mld * (MULTIPLEX_LD / BITCT2);
  int indiv_trail_ct = indiv_ct_mld_long - indiv_ctbit * 2;
  int retval = 0;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned long* geno = NULL;
  unsigned long* pruned_arr;
  unsigned int* live_indices;
  unsigned int* start_arr;
  unsigned int marker_unfiltered_idx;
  unsigned int marker_idx;
  int pct;
  unsigned int pct_thresh;
  int pairwise = calculation_type & CALC_LD_PRUNE_PAIRWISE;
  unsigned int window_unfiltered_start;
  unsigned int window_unfiltered_end;
  int cur_window_size;
  int old_window_size;
  int ii;
  int jj;
  int kk;
  unsigned int uii;
  unsigned int ujj;
  unsigned int cur_chrom;
  unsigned int chrom_end;
  double* marker_stdevs;
  unsigned char* loadbuf;
  unsigned int* missing_cts;
  unsigned long window_max = 0;
  unsigned long ulii;
  double dxx;
  double cov12;
  unsigned int fixed_non_missing_ct;
  unsigned int non_missing_ct;
  int dp_result[3];
  double non_missing_recip;
#if __LP64__
  __m128i* geno_fixed_vec_ptr;
  __m128i* geno_var_vec_ptr;
  __m128i* mask_fixed_vec_ptr;
  __m128i* mask_var_vec_ptr;
#else
  unsigned long* geno_fixed_vec_ptr;
  unsigned long* geno_var_vec_ptr;
  unsigned long* mask_fixed_vec_ptr;
  unsigned long* mask_var_vec_ptr;
#endif
  unsigned int cur_exclude_ct;
  int tot_exclude_ct = 0;
  int prev_end;
  int at_least_one_prune = 0;

  double* cov_matrix = NULL;
  double* new_cov_matrix = NULL;
  MATRIX_INVERT_BUF1_TYPE* irow = NULL;
  __CLPK_integer window_rem_li;
  __CLPK_integer old_window_rem_li;
  double* work = NULL;
  int window_rem;
  unsigned int* idx_remap = NULL;
  double prune_ld_r2;

  if (ld_window_kb) {
    // determine maximum number of markers that may need to be loaded at once
    for (cur_chrom = 0; cur_chrom < MAX_POSSIBLE_CHROM; cur_chrom++) {
      if (chrom_exists(chrom_info_ptr, cur_chrom)) {
        uii = chrom_info_ptr->chrom_start[cur_chrom];
	chrom_end = chrom_info_ptr->chrom_end[cur_chrom];
        do {
	  ujj = uii + 1;
	  while ((ujj < chrom_end) && (marker_pos[ujj] <= marker_pos[uii] + (1000 * ld_window_size))) {
	    ujj++;
	  }
          if (ujj - uii > window_max) {
	    window_max = ujj - uii;
	  }
	  uii++;
	} while (ujj < chrom_end);
      }
    }
  }
  if (pairwise) {
    prune_ld_r2 = sqrt(ld_last_param);
  } else {
    prune_ld_r2 = 0.999999;
  }

  window_unfiltered_start = ld_prune_next_valid_chrom_start(marker_exclude, 0, chrom_info_ptr, unfiltered_marker_ct);
  if (window_unfiltered_start == unfiltered_marker_ct) {
    if (pairwise) {
      logprint("Error: No valid markers for --indep-pairwise.\n");
    } else {
      logprint("Error: No valid markers for --indep.\n");
    }
    return RET_INVALID_FORMAT;
  }

  if (wkspace_alloc_ul_checked(&pruned_arr, unfiltered_marker_ctl * sizeof(long))) {
    goto ld_prune_ret_NOMEM;
  }

  memcpy(pruned_arr, marker_exclude, unfiltered_marker_ctl * sizeof(long));

  if (!ld_window_kb) {
    window_max = ld_window_size;
  }
  ulii = window_max;
  if (wkspace_alloc_ui_checked(&live_indices, ulii * sizeof(int))) {
    goto ld_prune_ret_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&start_arr, ulii * sizeof(int))) {
    goto ld_prune_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&marker_stdevs, ulii * sizeof(double))) {
    goto ld_prune_ret_NOMEM;
  }
  if (wkspace_alloc_uc_checked(&loadbuf, unfiltered_indiv_ct4)) {
    goto ld_prune_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&geno, ulii * indiv_ct_mld_long * sizeof(long))) {
    goto ld_prune_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&g_masks, ulii * indiv_ct_mld_long * sizeof(long))) {
    goto ld_prune_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&g_mmasks, ulii * indiv_ctbit * sizeof(long))) {
    goto ld_prune_ret_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&missing_cts, g_indiv_ct * sizeof(int))) {
    goto ld_prune_ret_NOMEM;
  }
  if (!pairwise) {
    if (wkspace_alloc_d_checked(&cov_matrix, window_max * window_max * sizeof(double))) {
      goto ld_prune_ret_NOMEM;
    }
    if (wkspace_alloc_d_checked(&new_cov_matrix, window_max * window_max * sizeof(double))) {
      goto ld_prune_ret_NOMEM;
    }
    if (wkspace_alloc_ui_checked(&idx_remap, window_max * sizeof(int))) {
      goto ld_prune_ret_NOMEM;
    }

    irow = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(window_max * 2 * sizeof(MATRIX_INVERT_BUF1_TYPE));
    if (!irow) {
      goto ld_prune_ret_NOMEM;
    }

    if (window_max < 4) {
      ulii = 4;
    } else {
      ulii = window_max;
    }
    if (wkspace_alloc_d_checked(&work, ulii * window_max * sizeof(double))) {
      goto ld_prune_ret_NOMEM;
    }
  }
  do {
    prev_end = 0;
    ld_prune_start_chrom(ld_window_kb, &cur_chrom, &chrom_end, window_unfiltered_start, live_indices, start_arr, &window_unfiltered_end, ld_window_size, &cur_window_size, unfiltered_marker_ct, pruned_arr, chrom_info_ptr, marker_pos);
    old_window_size = 1;
    if (cur_window_size > 1) {
      for (ulii = 0; ulii < (unsigned long)cur_window_size; ulii++) {
	if (fseeko(bedfile, bed_offset + (live_indices[ulii] * unfiltered_indiv_ct4), SEEK_SET)) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  goto ld_prune_ret_READ_FAIL;
	}
        missing_cts[ulii] = ld_process_load(loadbuf, &(geno[ulii * indiv_ct_mld_long]), &(g_masks[ulii * indiv_ct_mld_long]), &(g_mmasks[ulii * indiv_ctbit]), &(marker_stdevs[ulii]), unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, indiv_ctbit, indiv_trail_ct);
      }
    }
    pct = 1;
    pct_thresh = window_unfiltered_start + ((long long)pct * (chrom_end - chrom_info_ptr->chrom_start[cur_chrom])) / 100;
    cur_exclude_ct = 0;
    while ((window_unfiltered_start < chrom_end) || (cur_window_size > 1)) {
      if (cur_window_size > 1) {
	for (ii = 0; ii < cur_window_size; ii++) {
	  if (marker_stdevs[ii] == 0.0) {
	    set_bit(pruned_arr, live_indices[ii], &cur_exclude_ct);
	  }
	}
	do {
	  at_least_one_prune = 0;
	  for (ii = 0; ii < cur_window_size - 1; ii++) {
	    if (is_set(pruned_arr, live_indices[ii])) {
	      continue;
	    }
	    fixed_non_missing_ct = g_indiv_ct - missing_cts[ii];
#if __LP64__
	    geno_fixed_vec_ptr = (__m128i*)(&(geno[ii * indiv_ct_mld_long]));
	    mask_fixed_vec_ptr = (__m128i*)(&(g_masks[ii * indiv_ct_mld_long]));
#else
	    geno_fixed_vec_ptr = &(geno[ii * indiv_ct_mld_long]);
	    mask_fixed_vec_ptr = &(g_masks[ii * indiv_ct_mld_long]);
#endif
	    jj = ii + 1;
	    while (live_indices[jj] < start_arr[ii]) {
	      jj++;
	      if (jj == cur_window_size) {
		break;
	      }
	    }
	    for (; jj < cur_window_size; jj++) {
	      if (is_set(pruned_arr, live_indices[jj])) {
		continue;
	      }
#if __LP64__
	      geno_var_vec_ptr = (__m128i*)(&(geno[jj * indiv_ct_mld_long]));
	      mask_var_vec_ptr = (__m128i*)(&(g_masks[jj * indiv_ct_mld_long]));
#else
	      geno_var_vec_ptr = &(geno[jj * indiv_ct_mld_long]);
	      mask_var_vec_ptr = &(g_masks[jj * indiv_ct_mld_long]);
#endif

	      non_missing_ct = fixed_non_missing_ct - missing_cts[jj] + sparse_intersection_ct(&(g_mmasks[ii * indiv_ctbit]), &(g_mmasks[jj * indiv_ctbit]), indiv_ctbit);
	      dp_result[0] = g_indiv_ct;
	      // reversed from what I initially thought because I'm passing the
	      // jj-associated buffers before the ii-associated ones.
	      dp_result[1] = -fixed_non_missing_ct;
	      dp_result[2] = missing_cts[jj] - g_indiv_ct;
#if __LP64__
	      for (kk = 0; kk < indiv_ct_mld_m1; kk++) {
		ld_dot_prod(geno_var_vec_ptr, &(geno_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT)]), mask_var_vec_ptr, &(mask_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT)]), dp_result, MULTIPLEX_LD / 192);
		geno_var_vec_ptr = &(geno_var_vec_ptr[MULTIPLEX_LD / BITCT]);
		mask_var_vec_ptr = &(mask_var_vec_ptr[MULTIPLEX_LD / BITCT]);
	      }
	      ld_dot_prod(geno_var_vec_ptr, &(geno_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT)]), mask_var_vec_ptr, &(mask_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT)]), dp_result, indiv_ct_mld_rem);
	      geno_var_vec_ptr = &(geno_var_vec_ptr[MULTIPLEX_LD / BITCT]);
	      mask_var_vec_ptr = &(mask_var_vec_ptr[MULTIPLEX_LD / BITCT]);
#else
	      for (kk = 0; kk < indiv_ct_mld_m1; kk++) {
		ld_dot_prod(geno_var_vec_ptr, &(geno_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT2)]), mask_var_vec_ptr, &(mask_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT2)]), dp_result, MULTIPLEX_LD / 48);
		geno_var_vec_ptr = &(geno_var_vec_ptr[MULTIPLEX_LD / BITCT2]);
		mask_var_vec_ptr = &(mask_var_vec_ptr[MULTIPLEX_LD / BITCT2]);
	      }
	      ld_dot_prod(geno_var_vec_ptr, &(geno_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT2)]), mask_var_vec_ptr, &(mask_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT2)]), dp_result, indiv_ct_mld_rem);
	      geno_var_vec_ptr = &(geno_var_vec_ptr[MULTIPLEX_LD / BITCT2]);
	      mask_var_vec_ptr = &(mask_var_vec_ptr[MULTIPLEX_LD / BITCT2]);
#endif
	      non_missing_recip = 1.0 / ((double)non_missing_ct);
	      cov12 = non_missing_recip * (dp_result[0] - (non_missing_recip * dp_result[1]) * dp_result[2]);
	      // r-squared
	      dxx = cov12 / (marker_stdevs[ii] * marker_stdevs[jj]);
	      if (!pairwise) {
		cov_matrix[ii * window_max + jj] = dxx;
	      }
	      if (fabs(dxx) > prune_ld_r2) {
		at_least_one_prune = 1;
		// remove SNP with lower MAF
		if (get_maf(set_allele_freqs[live_indices[ii]]) < get_maf(set_allele_freqs[live_indices[jj]])) {
		  set_bit(pruned_arr, live_indices[ii], &cur_exclude_ct);
		} else {
		  set_bit(pruned_arr, live_indices[jj], &cur_exclude_ct);
		  jj++;
		  while (jj < cur_window_size) {
		    if (!is_set(pruned_arr, live_indices[jj])) {
		      break;
		    }
		    jj++;
		  }
		  if (jj < cur_window_size) {
		    start_arr[ii] = live_indices[jj];
		  }
		}
		break;
	      }
	    }
	    if (jj == cur_window_size) {
	      start_arr[ii] = window_unfiltered_end;
	    }
	  }
	} while (at_least_one_prune);
	if (!pairwise) {
	  window_rem = 0;
	  old_window_rem_li = 0;
	  for (ii = 0; ii < old_window_size; ii++) {
	    if (is_set(pruned_arr, live_indices[ii])) {
	      continue;
	    }
            idx_remap[window_rem++] = ii;
	  }
	  old_window_rem_li = window_rem;
	  for (; ii < cur_window_size; ii++) {
	    if (is_set(pruned_arr, live_indices[ii])) {
	      continue;
	    }
            idx_remap[window_rem++] = ii;
	  }
	  while (window_rem > 1) {
	    new_cov_matrix[0] = 1.0;
	    for (ii = 1; ii < window_rem; ii++) {
	      kk = idx_remap[ii];
	      for (jj = 0; jj < ii; jj++) {
		dxx = cov_matrix[idx_remap[jj] * window_max + kk];
		new_cov_matrix[jj * window_rem + ii] = dxx;
		new_cov_matrix[ii * window_rem + jj] = dxx;
	      }
	      new_cov_matrix[ii * (window_rem + 1)] = 1.0;
	    }
	    window_rem_li = window_rem;
	    jj = invert_matrix_trunc_singular(window_rem_li, new_cov_matrix, irow, work, old_window_rem_li);
	    while (jj) {
	      if (jj == -1) {
		goto ld_prune_ret_NOMEM;
	      }
              set_bit(pruned_arr, live_indices[idx_remap[jj]], &cur_exclude_ct);
	      window_rem--;
	      for (ii = jj; ii < window_rem; ii++) {
		idx_remap[ii] = idx_remap[ii + 1];
	      }
	      new_cov_matrix[0] = 1.0;
	      for (ii = 1; ii < window_rem; ii++) {
		kk = idx_remap[ii];
		for (jj = 0; jj < ii; jj++) {
		  dxx = cov_matrix[idx_remap[jj] * window_max + kk];
		  new_cov_matrix[jj * window_rem + ii] = dxx;
		  new_cov_matrix[ii * window_rem + jj] = dxx;
		}
		new_cov_matrix[ii * (window_rem + 1)] = 1.0;
	      }
              window_rem_li = window_rem;
	      jj = invert_matrix_trunc_singular(window_rem_li, new_cov_matrix, irow, work, old_window_rem_li);
	    }
	    dxx = new_cov_matrix[0];
	    jj = 0;
	    for (ii = 1; ii < window_rem; ii++) {
              if (new_cov_matrix[ii * (window_rem + 1)] > dxx) {
		dxx = new_cov_matrix[ii * (window_rem + 1)];
		jj = ii;
	      }
	    }
	    if (dxx > ld_last_param) {
	      set_bit(pruned_arr, live_indices[idx_remap[jj]], &cur_exclude_ct);
	      window_rem--;
	      if (idx_remap[jj] < (unsigned int)old_window_size) {
		old_window_rem_li--;
	      }
	      for (ii = jj; ii < window_rem; ii++) {
                idx_remap[ii] = idx_remap[ii + 1];
	      }
	    } else {
	      // break out
	      window_rem = 1;
	    }
	  }
	}
      }
      for (ii = 0; ii < ld_window_incr; ii++) {
	while (is_set(marker_exclude, window_unfiltered_start)) {
	  if (window_unfiltered_start == chrom_end) {
	    break;
	  }
	  window_unfiltered_start++;
	}
	if (window_unfiltered_start == chrom_end) {
	  break;
	}
	window_unfiltered_start++;
      }
      if (window_unfiltered_start == chrom_end) {
	break;
      }
      if (window_unfiltered_start >= pct_thresh) {
	pct = (((long long)(window_unfiltered_start - chrom_info_ptr->chrom_start[cur_chrom])) * 100) / (chrom_end - chrom_info_ptr->chrom_start[cur_chrom]);
	printf("\r%d%%", pct++);
	fflush(stdout);
	pct_thresh = chrom_info_ptr->chrom_start[cur_chrom] + (((long long)pct * (chrom_end - chrom_info_ptr->chrom_start[cur_chrom])) / 100);
      }
      jj = 0;
      // copy back previously loaded/computed results
      while (live_indices[jj] < window_unfiltered_start) {
	jj++;
	if (jj == cur_window_size) {
	  break;
	}
      }
      for (ii = 0; jj < cur_window_size; jj++) {
	if (is_set(pruned_arr, live_indices[jj])) {
	  continue;
	}
	memcpy(&(geno[ii * indiv_ct_mld_long]), &(geno[jj * indiv_ct_mld_long]), indiv_ct_mld_long * sizeof(long));
	memcpy(&(g_masks[ii * indiv_ct_mld_long]), &(g_masks[jj * indiv_ct_mld_long]), indiv_ct_mld_long * sizeof(long));
	memcpy(&(g_mmasks[ii * indiv_ctbit]), &(g_mmasks[jj * indiv_ctbit]), indiv_ctbit * sizeof(long));
	marker_stdevs[ii] = marker_stdevs[jj];
	live_indices[ii] = live_indices[jj];
	start_arr[ii] = start_arr[jj];
	missing_cts[ii] = missing_cts[jj];
	if (!pairwise) {
	  for (kk = 0; kk < ii; kk++) {
	    cov_matrix[kk * window_max + ii] = cov_matrix[idx_remap[kk] * window_max + jj];
	  }
	  idx_remap[ii] = jj;
	}
	ii++;
      }

      prev_end = ii;
      cur_window_size = ii;
      if (ld_window_kb) {
	jj = 0;
	while ((window_unfiltered_end + jj < chrom_end) && (marker_pos[window_unfiltered_end + jj] <= marker_pos[window_unfiltered_start] + (1000 * ld_window_size))) {
	  jj++;
	}
      } else {
	jj = ld_window_incr;
      }
      old_window_size = cur_window_size;
      for (ii = 0; ii < jj; ii++) {
	while ((window_unfiltered_end < chrom_end) && is_set(marker_exclude, window_unfiltered_end)) {
	  window_unfiltered_end++;
	}
	if (window_unfiltered_end < chrom_end) {
	  live_indices[cur_window_size] = window_unfiltered_end;
	  if (cur_window_size > prev_end) {
	    start_arr[cur_window_size - 1] = window_unfiltered_end;
	  }
	  if (fseeko(bedfile, bed_offset + (window_unfiltered_end * unfiltered_indiv_ct4), SEEK_SET)) {
	    goto ld_prune_ret_READ_FAIL;
	  }
	  if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	    goto ld_prune_ret_READ_FAIL;
	  }
	  missing_cts[cur_window_size] = ld_process_load(loadbuf, &(geno[cur_window_size * indiv_ct_mld_long]), &(g_masks[cur_window_size * indiv_ct_mld_long]), &(g_mmasks[cur_window_size * indiv_ctbit]), &(marker_stdevs[cur_window_size]), unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, indiv_ctbit, indiv_trail_ct);
	  cur_window_size++;
	  window_unfiltered_end++;
	}
      }
      if (cur_window_size > prev_end) {
	start_arr[cur_window_size] = window_unfiltered_end;
      }
    }
    ii = get_marker_chrom(chrom_info_ptr, window_unfiltered_start - 1);
    putchar('\r');
    sprintf(logbuf, "Pruned %d SNPs from chromosome %d, leaving %d.\n", cur_exclude_ct, ii, chrom_info_ptr->chrom_end[ii] - chrom_info_ptr->chrom_start[ii] - cur_exclude_ct);
    logprintb();
    tot_exclude_ct += cur_exclude_ct;

    // advance chromosomes as necessary
    window_unfiltered_start = ld_prune_next_valid_chrom_start(pruned_arr, window_unfiltered_start, chrom_info_ptr, unfiltered_marker_ct);
  } while (window_unfiltered_start < unfiltered_marker_ct);

  sprintf(logbuf, "Pruning complete.  %d of %d SNPs removed.\n", tot_exclude_ct, marker_ct);
  logprintb();
  strcpy(outname_end, ".prune.in");
  if (fopen_checked(&outfile_in, outname, "w")) {
    goto ld_prune_ret_OPEN_FAIL;
  }
  strcpy(outname_end, ".prune.out");
  if (fopen_checked(&outfile_out, outname, "w")) {
    goto ld_prune_ret_OPEN_FAIL;
  }
  marker_unfiltered_idx = 0;
  marker_idx = 0;
  pct = 1;
  ii = 0;
  for (cur_chrom = 1; cur_chrom < MAX_POSSIBLE_CHROM; cur_chrom++) {
    if (!(chrom_info_ptr->chrom_mask & (1LLU << cur_chrom))) {
      continue;
    }
    if (chrom_info_ptr->chrom_end[cur_chrom]) {
      ii += chrom_info_ptr->chrom_end[cur_chrom] - chrom_info_ptr->chrom_start[cur_chrom];
    }
  }
  pct_thresh = ((long long)pct * ii) / 100;
  for (cur_chrom = 1; cur_chrom < MAX_POSSIBLE_CHROM; cur_chrom++) {
    chrom_end = chrom_info_ptr->chrom_end[cur_chrom];
    if (!chrom_end) {
      continue;
    }
    marker_unfiltered_idx = chrom_info_ptr->chrom_start[cur_chrom];
    for (; marker_unfiltered_idx < chrom_end; marker_unfiltered_idx++) {
      if (!is_set(marker_exclude, marker_unfiltered_idx)) {
	if (is_set(pruned_arr, marker_unfiltered_idx)) {
	  if (fprintf(outfile_out, "%s\n", &(marker_ids[marker_unfiltered_idx * max_marker_id_len])) < 0) {
	    goto ld_prune_ret_WRITE_FAIL;
	  }
	} else {
	  if (fprintf(outfile_in, "%s\n", &(marker_ids[marker_unfiltered_idx * max_marker_id_len])) < 0) {
	    goto ld_prune_ret_WRITE_FAIL;
	  }
	}
      }
      marker_idx++;
      if (marker_idx == pct_thresh) {
	printf("\rWriting... %d%%", pct);
	fflush(stdout);
	pct = ((long long)marker_idx * 100) / ii + 1;
        pct_thresh = ((long long)pct * ii) / 100;
      }
    }
  }
  if (fclose_null(&outfile_in)) {
    goto ld_prune_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile_out)) {
    goto ld_prune_ret_WRITE_FAIL;
  }
  *outname_end = '\0';
  putchar('\r');
  sprintf(logbuf, "SNP lists written to %s.prune.in and %s.prune.out.\n", outname, outname);
  logprintb();

  while (0) {
  ld_prune_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_prune_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ld_prune_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_prune_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile_in);
  fclose_cond(outfile_out);
  wkspace_reset(wkspace_mark);
  return retval;
}

inline void rel_cut_arr_dec(int* rel_ct_arr_elem, unsigned int* exactly_one_rel_ct_ptr) {
  int rcae = *rel_ct_arr_elem - 1;
  *rel_ct_arr_elem = rcae;
  if (rcae < 2) {
    if (rcae) {
      *exactly_one_rel_ct_ptr += 1;
    } else {
      *exactly_one_rel_ct_ptr -= 1;
    }
  }
}

int do_rel_cutoff(int calculation_type, double rel_cutoff, double* rel_ibc, unsigned long* indiv_exclude, unsigned int* indiv_exclude_ct_ptr, char* outname, char* outname_end, unsigned int unfiltered_indiv_ct, char* person_ids, unsigned int max_person_id_len) {
  int indivs_excluded = 0;
  unsigned int exactly_one_rel_ct = 0;
  unsigned char* wkspace_mark = wkspace_base;
  double* dist_ptr = g_rel_dists;
  double* dptr2;
  double* dptr3;
  double* dptr4;
  unsigned int* giptr;
  unsigned int* giptr2;
  // number of too-close relations, -1 if excluded
  int* rel_ct_arr;
  unsigned int indiv_idx;
  unsigned int uii;
  unsigned long ulii;
  int kk;
  int mm;
  int retval;
  
  // Algorithm:
  // - Whenever there is at least one individual with exactly one
  // remaining too-close relation, prune the other side of that
  // relationship, because doing so is never suboptimal.
  // - Otherwise, there's no efficient rule that is always optimal
  // (assuming P != NP, anyway), so we use a simple heuristic: prune the
  // first individual with the largest number of remaining too-close
  // relationships.

  if (wkspace_alloc_i_checked(&rel_ct_arr, g_indiv_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  fill_int_zero(rel_ct_arr, g_indiv_ct);
  dist_ptr = g_rel_dists;
  for (indiv_idx = 1; indiv_idx < g_indiv_ct; indiv_idx++) {
    for (uii = 0; uii < indiv_idx; uii++) {
      if (*dist_ptr++ > rel_cutoff) {
	rel_ct_arr[indiv_idx] += 1;
	rel_ct_arr[uii] += 1;
      }
    }
  }
  for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
    if (rel_ct_arr[indiv_idx] == 1) {
      exactly_one_rel_ct++;
    }
  }
  while (1) {
    kk = 0;
    if (exactly_one_rel_ct) {
      // there is at least one individual with exactly one too-close
      // relation left, find the first one
      while (rel_ct_arr[kk] != 1) {
	kk++;
      }
      // and now find the identity of the other side
      dist_ptr = &(g_rel_dists[((long)kk * (kk - 1)) / 2]);
      for (mm = 0; mm < kk; mm++) {
	if (*dist_ptr > rel_cutoff) {
	  *dist_ptr = 0.0;
	  break;
	}
	dist_ptr++;
      }
      if (mm == kk) {
	do {
	  mm++;
	  dist_ptr = &(g_rel_dists[((long)mm * (mm - 1)) / 2 + kk]);
	} while (*dist_ptr <= rel_cutoff);
	*dist_ptr = 0.0;
      }
      rel_ct_arr[kk] = 0;
      exactly_one_rel_ct--;
      if (rel_ct_arr[mm] == 1) {
        // speed up the easy case
	exactly_one_rel_ct--;
	rel_ct_arr[mm] = -1;
	indivs_excluded++;
	continue;
      }
    } else {
      // find identity of first individual with maximum number of
      // remaining too-close relations
      // kk is highest too-close pair count so far
      mm = -1; // associated individual index
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	if (rel_ct_arr[indiv_idx] > kk) {
	  kk = rel_ct_arr[indiv_idx];
	  mm = indiv_idx;
	}
      }
      // no too-close relations left at all, we're done
      if (mm == -1) {
	break;
      }
    }
    dist_ptr = &(g_rel_dists[((long)mm * (mm - 1)) / 2]);
    for (kk = 0; kk < mm; kk++) {
      if (*dist_ptr > rel_cutoff) {
	*dist_ptr = 0.0;
	rel_cut_arr_dec(&(rel_ct_arr[kk]), &exactly_one_rel_ct);
      }
      dist_ptr++;
    }
    for (ulii = mm + 1; ulii < g_indiv_ct; ulii++) {
      dist_ptr = &(g_rel_dists[(ulii * (ulii - 1)) / 2 + mm]);
      if (*dist_ptr > rel_cutoff) {
	*dist_ptr = 0.0;
	rel_cut_arr_dec(&(rel_ct_arr[ulii]), &exactly_one_rel_ct);
      }
    }
    rel_ct_arr[mm] = -1;
    indivs_excluded++;
  }
  exclude_multi(indiv_exclude, rel_ct_arr, g_indiv_ct, indiv_exclude_ct_ptr);
  if (indivs_excluded) {
    dist_ptr = g_rel_dists; // write
    dptr2 = g_rel_dists; // read
    dptr3 = rel_ibc; // write
    dptr4 = rel_ibc; // read
    giptr = g_missing_dbl_excluded; // write
    giptr2 = g_missing_dbl_excluded; // read
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      if (rel_ct_arr[indiv_idx] != -1) {
	if (calculation_type & CALC_IBC) {
	  dptr3[g_indiv_ct] = dptr4[g_indiv_ct];
	  dptr3[g_indiv_ct * 2] = dptr4[g_indiv_ct * 2];
	}
	*dptr3 = *dptr4++;
	dptr3++;
	for (uii = 0; uii < indiv_idx; uii++) {
	  if (rel_ct_arr[uii] != -1) {
	    *dist_ptr = *dptr2++;
	    dist_ptr++;
	    *giptr = *giptr2++;
	    giptr++;
	  } else {
	    dptr2++;
	    giptr2++;
	  }
	}
      } else {
	dptr4++;
	dptr2 = &(dptr2[indiv_idx]);
	giptr2 = &(giptr2[indiv_idx]);
      }
    }
    g_indiv_ct -= indivs_excluded;
    if (calculation_type & CALC_IBC) {
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	*dptr3++ = *dptr4++;
      }
      dptr4 = &(dptr4[indivs_excluded]);
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	*dptr3++ = *dptr4++;
      }
    }
    giptr = g_indiv_missing_unwt;
    giptr2 = g_indiv_missing_unwt;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct + indivs_excluded; indiv_idx++) {
      if (rel_ct_arr[indiv_idx] != -1) {
	*giptr = *giptr2++;
	giptr++;
      } else {
	giptr2++;
      }
    }
  }
  sprintf(logbuf, "%d individual%s excluded by --rel-cutoff.\n", indivs_excluded, (indivs_excluded == 1)? "" : "s");
  logprintb();
  if (!(calculation_type & (CALC_RELATIONSHIP | CALC_GDISTANCE_MASK))) {
    strcpy(outname_end, ".rel.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
    sprintf(logbuf, "Remaining individual IDs written to %s.\n", outname);
    logprintb();
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int rel_cutoff_batch(char* grmname, char* outname, double rel_cutoff, int rel_calc_type) {
  // Specialized --rel-cutoff usable on larger files.
  char* grmname_end = (char*)memchr(grmname, 0, FNAMESIZE);
  char* outname_end = (char*)memchr(outname, 0, FNAMESIZE);
  unsigned int indiv_ct = 0;
  FILE* idfile = NULL;
  FILE* outfile = NULL;
  gzFile cur_gzfile = NULL;
  gzFile gz_outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned int indivs_excluded = 0;
  unsigned int exactly_one_rel_ct = 0;
  char wbuf[64];
  unsigned long* compact_rel_table;
  unsigned long* rtptr;
  char* bufptr;
  char* bufptr2;
  unsigned long long ullii;
  unsigned long long ulljj;
  unsigned long tot_words;
  unsigned long words_left;
  unsigned long wl_floor;
  unsigned long cur_word;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long ulkk;
  unsigned int inword_idx;
  unsigned int inword_bound;
  unsigned int uii;
  unsigned int row;
  unsigned int col;
  unsigned int indiv_idx;
  int* rel_ct_arr;
  double dxx;
  int retval;
  unsigned int pct;
  int kk;
  int mm;
  int cur_prune;
  memcpy(grmname_end, ".grm.id", 8);
  if (fopen_checked(&idfile, grmname, "r")) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, idfile)) {
    if (is_eoln(*tbuf)) {
      continue;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in .grm.id file.\n");
      goto rel_cutoff_batch_ret_INVALID_FORMAT_1;
    }
    indiv_ct++;
  }
  if (!feof(idfile)) {
    goto rel_cutoff_batch_ret_READ_FAIL;
  }
  fclose_null(&idfile);
  ullii = indiv_ct;
  ullii = ((ullii * (ullii - 1)) / 2 + BITCT - 1) / BITCT;
#ifndef __LP64__
  if (ullii >= 0x20000000) {
    goto rel_cutoff_batch_ret_NOMEM;
  }
#endif
  tot_words = ullii;
  if (wkspace_alloc_ul_checked(&compact_rel_table, tot_words * sizeof(long))) {
    goto rel_cutoff_batch_ret_NOMEM;
  }
  fill_ulong_zero(compact_rel_table, tot_words);
  if (wkspace_alloc_i_checked(&rel_ct_arr, indiv_ct * sizeof(int))) {
    goto rel_cutoff_batch_ret_NOMEM;
  }
  fill_int_zero(rel_ct_arr, indiv_ct);

  memcpy(grmname_end, ".grm.gz", 8);
  if (gzopen_checked(&cur_gzfile, grmname, "rb")) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }
  if (gzbuffer(cur_gzfile, 131072)) {
    goto rel_cutoff_batch_ret_NOMEM;
  }

  words_left = tot_words;
  rtptr = compact_rel_table;
  fputs("Reading... 0%", stdout);
  fflush(stdout);
  row = 0;
  col = 0;
  for (pct = 1; pct <= 100; pct++) {
    wl_floor = (((unsigned long long)tot_words) * (100 - pct)) / 100;
    while (words_left > wl_floor) {
      cur_word = 0;
      if (--words_left) {
	inword_bound = BITCT;
      } else {
	// only indiv_ct mod (BITCT * 2) matters for remainder
	uii = indiv_ct & (BITCT * 2 - 1);
	inword_bound = ((uii * (uii - 1)) / 2) & (BITCT - 1);
	if (!inword_bound) {
	  inword_bound = BITCT;
	}
      }
      for (inword_idx = 0; inword_idx < inword_bound; inword_idx++) {
	if (!gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
	  goto rel_cutoff_batch_ret_READ_FAIL;
	}
	if (row == col) {
	  row++;
	  col = 0;
	  if (!gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
	    goto rel_cutoff_batch_ret_READ_FAIL;
	  }
	}
	bufptr = next_item_mult(tbuf, 3);
	if (no_more_items(bufptr)) {
	  goto rel_cutoff_batch_ret_INVALID_FORMAT_2;
	}
	if (sscanf(bufptr, "%lg", &dxx) != 1) {
	  goto rel_cutoff_batch_ret_INVALID_FORMAT_2;
	}
	if (dxx > rel_cutoff) {
	  rel_ct_arr[row] += 1;
	  rel_ct_arr[col] += 1;
	  cur_word |= (1LU << inword_idx);
	}
	col++;
      }
      *rtptr++ = cur_word;
    }
    if (pct <= 10) {
      printf("\b\b%u%%", pct);
      fflush(stdout);
    } else if (pct < 100) {
      printf("\b\b\b%u%%", pct);
      fflush(stdout);
    }
  }
  if (!gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
    goto rel_cutoff_batch_ret_READ_FAIL;
  }
  if (gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
    goto rel_cutoff_batch_ret_INVALID_FORMAT_2;
  }
  gzclose(cur_gzfile);
  cur_gzfile = NULL;
  putchar('\r');
  sprintf(logbuf, "%s read complete.  Pruning.\n", grmname);
  logprintb();

  // would prefer to just call do_rel_cutoff(), but unfortunately that
  // interferes with the intended "handle extra-large datasets" mission of this
  // function, which necessitates a compact bit representation of the
  // relationship matrix... fortunately, the algorithm is pretty simple.
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    if (rel_ct_arr[indiv_idx] == 1) {
      exactly_one_rel_ct++;
    }
  }

  while (1) {
    kk = 0;
    cur_prune = -1;
    if (exactly_one_rel_ct) {
      while (rel_ct_arr[kk] != 1) {
	kk++;
      }
      ullii = (((long long)kk) * (kk - 1)) / 2;
      ulii = ullii / BITCT;
      rtptr = &(compact_rel_table[ulii]);
      inword_idx = ullii & (BITCT - 1);
      ulljj = ullii + kk;
      uljj = ulljj / BITCT;
      inword_bound = ulljj & (BITCT - 1);
      if (uljj == ulii) {
	uljj = (1LU << inword_bound) - (1LU << inword_idx);
        ulkk = (*rtptr) & uljj;
	if (ulkk) {
	  *rtptr &= ~uljj;
	  cur_prune = __builtin_ctzl(ulkk) - inword_idx;
	}
      } else {
        ulkk = (*rtptr) & (~((1LU << inword_idx) - 1));
	if (ulkk) {
	  *rtptr &= (1LU << inword_idx) - 1;
          cur_prune = __builtin_ctzl(ulkk) - inword_idx;
	} else {
	  col = BITCT - inword_idx;
          row = col + (uljj - ulii - 1) * BITCT;
          while (col < row) {
            ulkk = *(++rtptr);
            if (ulkk) {
	      *rtptr = 0;
              cur_prune = __builtin_ctzl(ulkk) + col;
	      break;
	    }
	    col += BITCT;
	  }
	  if (cur_prune == -1) {
            ulkk = (*(++rtptr)) & ((1LU << inword_bound) - 1);
            if (ulkk) {
	      *rtptr &= (~((1LU << inword_bound) - 1));
	      cur_prune = __builtin_ctzl(ulkk) + col;
	    }
	  }
	}
      }
      if (cur_prune == -1) {
        mm = kk + 1;
        while (1) {
          ullii = ((((long long)mm) * (mm - 1)) / 2) + kk;
	  ulii = ullii / BITCT;
	  inword_idx = ullii & (BITCT - 1);
	  if (compact_rel_table[ulii] & (1LU << inword_idx)) {
	    compact_rel_table[ulii] &= ~(1LU << inword_idx);
	    cur_prune = mm;
	    break;
	  }
          mm++;
	}
      }
      rel_ct_arr[kk] = 0;
      exactly_one_rel_ct--;
      if (rel_ct_arr[cur_prune] == 1) {
	exactly_one_rel_ct--;
	rel_ct_arr[cur_prune] = -1;
	indivs_excluded++;
	continue;
      }
    } else {
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	if (rel_ct_arr[indiv_idx] > kk) {
	  kk = rel_ct_arr[indiv_idx];
	  cur_prune = indiv_idx;
	}
      }
      if (cur_prune == -1) {
	break;
      }
    }
    // zero out cur_prune row/column, update other array entries
    ullii = (((long long)cur_prune) * (cur_prune - 1)) / 2;
    ulii = ullii / BITCT;
    rtptr = &(compact_rel_table[ulii]);
    inword_idx = ullii & (BITCT - 1);
    ulljj = ullii + cur_prune;
    uljj = ulljj / BITCT;
    inword_bound = ulljj & (BITCT - 1);
    if (uljj == ulii) {
      uljj = (1LU << inword_bound) - (1LU << inword_idx);
      ulkk = (*rtptr) & uljj;
      if (ulkk) {
	do {
	  uii = __builtin_ctzl(ulkk) - inword_idx;
	  rel_cut_arr_dec(&(rel_ct_arr[uii]), &exactly_one_rel_ct);
	  ulkk &= ulkk - 1;
	} while (ulkk);
	*rtptr &= ~uljj;
      }
    } else {
      ulkk = (*rtptr) & (~((1LU << inword_idx) - 1));
      if (ulkk) {
	do {
	  uii = __builtin_ctzl(ulkk) - inword_idx;
	  rel_cut_arr_dec(&(rel_ct_arr[uii]), &exactly_one_rel_ct);
	  ulkk &= ulkk - 1;
	} while (ulkk);
	*rtptr &= (1LU << inword_idx) - 1;
      }
      col = BITCT - inword_idx;
      row = col + (uljj - ulii - 1) * BITCT;
      while (col < row) {
	ulkk = *(++rtptr);
	if (ulkk) {
	  do {
	    uii = __builtin_ctzl(ulkk) + col;
	    rel_cut_arr_dec(&(rel_ct_arr[uii]), &exactly_one_rel_ct);
	    ulkk &= ulkk - 1;
	  } while (ulkk);
	  *rtptr = 0;
	}
	col += BITCT;
      }
      ulkk = (*(++rtptr)) & ((1LU << inword_bound) - 1);
      if (ulkk) {
	do {
	  uii = __builtin_ctzl(ulkk) + col;
	  rel_cut_arr_dec(&(rel_ct_arr[uii]), &exactly_one_rel_ct);
          ulkk &= ulkk - 1;
	} while (ulkk);
	*rtptr &= (~((1LU << inword_bound) - 1));
      }
    }

    for (uljj = cur_prune + 1; uljj < indiv_ct; uljj++) {
      ullii = ((((unsigned long long)uljj) * (uljj - 1)) / 2) + cur_prune;
      ulii = ullii / BITCT;
      rtptr = &(compact_rel_table[ulii]);
      inword_idx = ullii & (BITCT - 1);
      if ((*rtptr) & (1LU << inword_idx)) {
        rel_cut_arr_dec(&(rel_ct_arr[uljj]), &exactly_one_rel_ct);
        *rtptr &= ~(1LU << inword_idx);
      }
    }
    rel_ct_arr[cur_prune] = -1;
    indivs_excluded++;
  }

  memcpy(grmname_end, ".grm.id", 8);
  if (fopen_checked(&idfile, grmname, "r")) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }

  memcpy(outname_end, ".grm.id", 8);
  if (fopen_checked(&outfile, outname, "w")) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }

  for (indiv_idx = 0; indiv_idx < indiv_ct;) {
    if (fgets(tbuf, MAXLINELEN, idfile) == NULL) {
      goto rel_cutoff_batch_ret_READ_FAIL;
    }
    if (is_eoln(*tbuf)) {
      continue;
    }
    if (rel_ct_arr[indiv_idx] != -1) {
      if (fprintf(outfile, "%s", tbuf) < 0) {
        goto rel_cutoff_batch_ret_WRITE_FAIL;
      }
    }
    indiv_idx++;
  }

  fclose_null(&idfile);
  fclose_null(&outfile);

  sprintf(logbuf, "%d individual%s excluded by --rel-cutoff.\n", indivs_excluded, (indivs_excluded == 1)? "" : "s");
  logprintb();
  sprintf(logbuf, "Remaining individual IDs written to %s.\n", outname);
  logprintb();
  if (rel_calc_type & REL_CALC_GRM) {
    memcpy(grmname_end, ".grm.gz", 8);
    if (gzopen_checked(&cur_gzfile, grmname, "rb")) {
      goto rel_cutoff_batch_ret_OPEN_FAIL;
    }
    if (gzbuffer(cur_gzfile, 131072)) {
      goto rel_cutoff_batch_ret_NOMEM;
    }

    if (rel_calc_type & REL_CALC_GZ) {
      memcpy(outname_end, ".grm.gz", 8);
      if (gzopen_checked(&gz_outfile, outname, "wb")) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
    } else {
      memcpy(outname_end, ".grm", 5);
      if (fopen_checked(&outfile, outname, "w")) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
    }
    kk = 1;
    ullii = (((unsigned long long)indiv_ct) * (indiv_ct - 1)) / 2;
    ulljj = 0;
    printf("Rewriting matrix... 0%%");
    fflush(stdout);
    pct = 1;
    for (row = 0; row < indiv_ct; row++) {
      if (rel_ct_arr[row] == -1) {
	for (col = 0; col <= row; col++) {
	  if (!gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
	    goto rel_cutoff_batch_ret_READ_FAIL;
	  }
	}
      } else {
	bufptr = &(wbuf[sprintf(wbuf, "%d\t", kk++)]);
	mm = 1;
	for (col = 0; col <= row; col++) {
	  if (!gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
	    goto rel_cutoff_batch_ret_READ_FAIL;
	  }
	  if (rel_ct_arr[col] != -1) {
	    bufptr2 = next_item_mult(tbuf, 2);
	    sprintf(bufptr, "%d\t%s", mm++, bufptr2);
	    if (flexwrite_checked(outfile, gz_outfile, wbuf, strlen(wbuf))) {
	      goto rel_cutoff_batch_ret_WRITE_FAIL;
	    }
	  }
	}
      }
      ulljj += row * 100;
      if (ulljj >= pct * ullii) {
	if (pct <= 10) {
	  pct = 1 + (ulljj / ullii);
	  printf("\b\b%u%%", pct - 1);
	  fflush(stdout);
	} else {
	  pct = 1 + (ulljj / ullii);
	  printf("\b\b\b%u%%", pct - 1);
	  fflush(stdout);
	}
      }
    }
    putchar('\r');
    sprintf(logbuf, "Pruned relationship matrix written to %s.\n", outname);
    logprintb();
  }
  retval = 0;
  while (0) {
  rel_cutoff_batch_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  rel_cutoff_batch_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  rel_cutoff_batch_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  rel_cutoff_batch_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  rel_cutoff_batch_ret_INVALID_FORMAT_2:
    putchar('\n');
    logprint("Error: Improperly formatted .grm.gz file.\n");
  rel_cutoff_batch_ret_INVALID_FORMAT_1:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(idfile);
  fclose_cond(outfile);
  gzclose_cond(cur_gzfile);
  gzclose_cond(gz_outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}


inline int distance_wt_req(int calculation_type) {
  return ((calculation_type & CALC_DISTANCE) || ((!(calculation_type & CALC_LOAD_DISTANCES)) && ((calculation_type & CALC_GROUPDIST) || (calculation_type & CALC_REGRESS_DISTANCE))));
}

inline int relationship_req(int calculation_type) {
  return (calculation_type & (CALC_RELATIONSHIP | CALC_UNRELATED_HERITABILITY | CALC_REL_CUTOFF | CALC_REGRESS_REL));
}

inline int relationship_or_ibc_req(int calculation_type) {
  return (relationship_req(calculation_type) || (calculation_type & CALC_IBC));
}

int wdist(char* outname, char* pedname, char* mapname, char* famname, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* filtername, char* freqname, char* loaddistname, char* evecname, char* mergename1, char* mergename2, char* mergename3, char* makepheno_str, char* phenoname_str, char* refalleles, char* filterval, int mfilter_col, int filter_case_control, int filter_sex, int filter_founder_nonf, int fam_col_1, int fam_col_34, int fam_col_5, int fam_col_6, char missing_geno, int missing_pheno, char output_missing_geno, char* output_missing_pheno, int mpheno_col, int pheno_merge, int prune, int affection_01, Chrom_info* chrom_info_ptr, double exponent, double min_maf, double max_maf, double geno_thresh, double mind_thresh, double hwe_thresh, int hwe_all, double rel_cutoff, int tail_pheno, double tail_bottom, double tail_top, int calculation_type, int rel_calc_type, int dist_calc_type, unsigned long groupdist_iters, int groupdist_d, unsigned long regress_iters, int regress_d, unsigned long regress_rel_iters, int regress_rel_d, double unrelated_herit_tol, double unrelated_herit_covg, double unrelated_herit_covr, int ibc_type, int parallel_idx, unsigned int parallel_tot, int ppc_gap, int allow_no_sex, int nonfounders, int genome_output_gz, int genome_output_full, int genome_ibd_unbounded, int ld_window_size, int ld_window_kb, int ld_window_incr, double ld_last_param, int maf_succ, int regress_pcs_normalize_pheno, int regress_pcs_sex_specific, int regress_pcs_clip, int max_pcs, int freq_counts, int freqx, int distance_flat_missing, int recode_modifier, int allelexxxx, int merge_type, int indiv_sort, int keep_allele_order) {
  FILE* outfile = NULL;
  FILE* outfile2 = NULL;
  FILE* outfile3 = NULL;
  FILE* mapfile = NULL;
  FILE* pedfile = NULL;
  FILE* famfile = NULL;
  gzFile gz_outfile = NULL;
  gzFile gz_outfile2 = NULL;
  gzFile gz_outfile3 = NULL;
  FILE* phenofile = NULL;
  FILE* filterfile = NULL;
  FILE* bedtmpfile = NULL;
  FILE* bimtmpfile = NULL;
  FILE* famtmpfile = NULL;
  FILE* freqfile = NULL;
  FILE* loaddistfile = NULL;
  char* id_buf = NULL;
  unsigned int unfiltered_marker_ct = 0;
  char* outname_end = (char*)memchr(outname, 0, FNAMESIZE);
  unsigned int marker_ct;
  unsigned char* pedbuf = NULL;
  unsigned long* marker_exclude = NULL;
  unsigned long long* line_locs = NULL;
  unsigned long long* line_mids = NULL;
  int pedbuflen;
  unsigned long max_marker_id_len = 0;
  unsigned int plink_maxsnp;
  // set_allele_freqs = frequency of allele corresponding to set bits in .bed
  //   (i.e. A2), or frequency of MAJOR allele in middle of text loading.
  double* set_allele_freqs = NULL;
  unsigned int unfiltered_indiv_ct = 0;
  unsigned int unfiltered_indiv_ct4 = 0;
  unsigned long* indiv_exclude = NULL;
  unsigned int indiv_exclude_ct = 0;
  unsigned int* indiv_sort_map = NULL;
  unsigned long* founder_info = NULL;
  unsigned char* sex_info = NULL;
  int ii;
  int jj = 0;
  int kk = 0;
  int mm;
  int nn = 0;
  int oo = 0;
  int pp;
  int qq;
  unsigned int uii = 0;
  unsigned int ujj;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long ulkk;
  double dxx;
  double dyy;
  double dzz;
  double dww;
  long long llxx = 0;
  long long llyy;
  char* marker_ids = NULL;
  unsigned char* marker_weights_base = NULL;
  // Binary files:
  //   marker_alleles[2 * ii] is id of A1 (usually minor) allele at SNP ii
  //   marker_alleles[2 * ii + 1] is identity of A2 allele at SNP ii
  // Text file mid-loading:
  //   marker_alleles[4 * ii] is identity of major allele at SNP ii
  //   marker_alleles[4 * ii + 1] is identity of 2nd most frequent allele, etc.
  char* marker_alleles = NULL;
  unsigned long* marker_reverse = NULL;
  unsigned int* marker_allele_cts;
  unsigned int* missing_cts;
  int retval = RET_SUCCESS;
  int map_is_unsorted = 0;
  int map_cols = 3;
  int affection = 0;
  double* phenor_d = NULL;
  char* phenor_c = NULL;
  char* person_ids = NULL;
  unsigned int max_person_id_len;
  char* paternal_ids = NULL;
  unsigned int max_paternal_id_len = 2;
  char* maternal_ids = NULL;
  unsigned int max_maternal_id_len = 2;
  unsigned char* wkspace_mark = NULL;
  unsigned int* giptr = NULL;
  unsigned int* giptr2 = NULL;
  unsigned int* giptr3;
  char* cptr = NULL;
  char* cptr2;
  unsigned char* gptr;
  unsigned char* gptr2;
  unsigned long* glptr;
  unsigned long* glptr2;
  unsigned long* glptr3;
  int* iptr;
  unsigned long long dists_alloc = 0;
  double* dist_ptr = NULL;
  double* dptr2;
  double* dptr3 = NULL;
  double* dptr4 = NULL;
  double* rel_ibc;
  int binary_files = 0;
  unsigned int marker_exclude_ct = 0;
  char* pid_list = NULL;
  char* id_list = NULL;
  double set_allele_freq_buf[MULTIPLEX_DIST];
  unsigned int wtbuf[MULTIPLEX_DIST];
  double missing_phenod = (double)missing_pheno;
  int missing_pheno_len = intlen(missing_pheno);
  int var_std = 1;
  int* hwe_lls;
  int* hwe_lhs;
  int* hwe_hhs;
  int* hwe_ll_allfs;
  int* hwe_lh_allfs;
  int* hwe_hh_allfs;
  int* ll_cts;
  int* lh_cts;
  int* hh_cts;
  int multiplex = 0;
  double* ll_pool;
  double* lh_pool;
  double* hh_pool;
  double* ll_poolp;
  double* lh_poolp;
  double* hh_poolp;
  double dhh_ssq;
  double dhl_ssq;
  double dll_ssq;
  double dhh_sd;
  double dhl_sd;
  double dll_sd;
  int ll_size;
  int lh_size;
  int hh_size;
  double ll_med;
  double lh_med;
  double hh_med;
  pthread_t threads[MAX_THREADS];
  int exp0 = (exponent == 0.0);
  int wt_needed = 0;
  int unwt_needed = 0;
  int unwt_needed_full = 0;
  int bed_offset = 3;
  unsigned int* marker_pos = NULL;
  Pedigree_rel_info pri;
  unsigned char* wkspace_mark2 = NULL;
  unsigned int marker_uidx;
  unsigned int marker_idx;
  unsigned int indiv_uidx;
  unsigned int indiv_idx;
  unsigned int pct;
  unsigned int tstc;
  int duplicate_fail;

  if (rel_calc_type & REL_CALC_COV) {
    var_std = 0;
    ibc_type = -1;
  }

  if (famname[0]) {
    binary_files = 1;
    if (calculation_type & CALC_MAKE_BED) {
      if (!realpath(pedname, tbuf)) {
	sprintf(logbuf, "Error: Failed to open %s.\n", pedname);
	logprintb();
        goto wdist_ret_OPEN_FAIL;
      }
      memcpy(outname_end, ".bed", 5);
      // if file doesn't exist, realpath returns NULL on Linux instead of what
      // the path would be.
      cptr = realpath(outname, &(tbuf[FNAMESIZE + 64]));
      if (cptr && (!strcmp(tbuf, &(tbuf[FNAMESIZE + 64])))) {
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
    binary_files = 1;
    uljj = (unsigned long)(outname_end - outname) + (ulii? 6 : 0);
    memcpy(pedname, outname, uljj);
    memcpy(&(pedname[uljj]), ".bed", 5);
    memcpy(famname, pedname, uljj);
    memcpy(&(famname[uljj]), ".fam", 5);
    memcpy(mapname, pedname, uljj);
    memcpy(&(mapname[uljj]), ".bim", 5);
  }

  if (fopen_checked(&pedfile, pedname, famname[0]? "rb" : "r")) {
    goto wdist_ret_OPEN_FAIL;
  }
  if (famname[0]) {
    if (fopen_checked(&famfile, famname, "r")) {
      goto wdist_ret_OPEN_FAIL;
    }
  }

  // load .map/.bim, count markers, filter chromosomes
  retval = load_map_or_bim(&mapfile, mapname, binary_files, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &plink_maxsnp, &marker_exclude, &set_allele_freqs, &marker_alleles, &marker_ids, chrom_info_ptr, &marker_pos, extractname, excludename, freqname, refalleles, calculation_type, recode_modifier, &map_is_unsorted);
  if (retval) {
    goto wdist_ret_2;
  }

  if (binary_files) {
    ulii = MAXLINELEN;
  } else {
    ulii = unfiltered_marker_ct * 4 + PEDBUFBASE;
  }
  // load .fam/[first few columns of .ped], count indivs
  ii = fam_col_6;
  if (ii && phenoname) {
    ii = pheno_merge && (!makepheno_str);
  }
  retval = load_fam(binary_files? famfile : pedfile, ulii, fam_col_1, fam_col_34, fam_col_5, ii, fam_col_6, missing_pheno, missing_pheno_len, affection_01, &unfiltered_indiv_ct, &person_ids, &max_person_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_info, &affection, &g_pheno_c, &g_pheno_d, &founder_info, &indiv_exclude, binary_files, &line_locs, &line_mids, &pedbuflen);
  if (retval) {
    goto wdist_ret_2;
  }
  count_genders(sex_info, unfiltered_indiv_ct, indiv_exclude, &ii, &jj, &kk);
  if (kk) {
    sprintf(logbuf, "%d markers and %d people (%d male%s, %d female%s, %d unknown) loaded.\n", unfiltered_marker_ct - marker_exclude_ct, unfiltered_indiv_ct, ii, (ii == 1)? "" : "s", jj, (jj == 1)? "" : "s", kk);
  } else {
    sprintf(logbuf, "%d markers and %d people (%d male%s, %d female%s) loaded.\n", unfiltered_marker_ct - marker_exclude_ct, unfiltered_indiv_ct, ii, (ii == 1)? "" : "s", jj, (jj == 1)? "" : "s");
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
      retval = makepheno_load(phenofile, makepheno_str, unfiltered_indiv_ct, indiv_exclude_ct, cptr, max_person_id_len, iptr, &g_pheno_c);
      if (retval) {
	goto wdist_ret_2;
      }
    } else if (phenofile) {
      retval = load_pheno(phenofile, unfiltered_indiv_ct, indiv_exclude_ct, cptr, max_person_id_len, iptr, missing_pheno, missing_pheno_len, affection_01, mpheno_col, phenoname_str, &g_pheno_c, &g_pheno_d);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    if (tail_pheno) {
      retval = convert_tail_pheno(unfiltered_indiv_ct, &g_pheno_c, &g_pheno_d, tail_bottom, tail_top, missing_phenod);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    wkspace_reset(wkspace_mark);
  }

  if ((calculation_type & CALC_GROUPDIST) && (!g_pheno_c)) {
    logprint("Error: --groupdist calculation requires dichotomous phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_REGRESS_DISTANCE) && (!g_pheno_d)) {
    logprint("Error: --regress-distance calculation requires scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_UNRELATED_HERITABILITY) && (!g_pheno_d)) {
    logprint("Error: --unrelated-heritability requires scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & (CALC_REGRESS_PCS | CALC_REGRESS_PCS_DISTANCE)) && (!g_pheno_d)) {
    sprintf(logbuf, "Error: --regress-pcs%s requires scalar phenotype.\n", (calculation_type & CALC_REGRESS_PCS_DISTANCE)? "-distance" : "");
    goto wdist_ret_INVALID_CMDLINE_2;
  }

  if (prune) {
    prune_missing_phenos(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, g_pheno_c, g_pheno_d, missing_phenod);
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
      indiv_sort_map = (unsigned int*)iptr;
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
    indiv_sort_map = (unsigned int*)iptr;
    wkspace_reset((unsigned char*)cptr);
  }

  if (filter_case_control) {
    if (!g_pheno_c) {
      logprint("Error: --filter-cases/--filter-controls requires dichotomous phenotype.\n");
      goto wdist_ret_INVALID_CMDLINE;
    }
    ii = indiv_exclude_ct;
    if (filter_indivs_var(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, g_pheno_c, NULL, 2 - filter_case_control)) {
      goto wdist_ret_NOMEM;
    }
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d individual%s removed due to case/control status (--filter-%s).\n", ii, (ii == 1)? "" : "s", (filter_case_control == 1)? "cases" : "controls");
    logprintb();
  }
  if (filter_sex) {
    ii = indiv_exclude_ct;
    if (filter_indivs_var(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, (char*)sex_info, NULL, filter_sex)) {
      goto wdist_ret_NOMEM;
    }
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d individual%s removed due to gender filter (--filter-%s).\n", ii, (ii == 1)? "" : "s", (filter_sex == 1)? "males" : "females");
    logprintb();
  }
  if (filter_founder_nonf) {
    ii = indiv_exclude_ct;
    if (filter_indivs_var(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, NULL, founder_info, filter_founder_nonf)) {
      goto wdist_ret_NOMEM;
    }
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d individual%s removed due to founder status (--filter-%s).\n", ii, (ii == 1)? "" : "s", (filter_founder_nonf == 1)? "founders" : "nonfounders");
    logprintb();
  }

  if (binary_files) {
    if (fseeko(pedfile, 0, SEEK_END)) {
      goto wdist_ret_READ_FAIL;
    }
    llxx = ftello(pedfile);
    llyy = llxx - ((unsigned long long)unfiltered_indiv_ct4) * unfiltered_marker_ct;
    rewind(pedfile);
    if (llyy == 3LL) {
      // v1.00 or later
      if (fread(tbuf, 1, 3, pedfile) < 3) {
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
      if (fread(tbuf, 1, 1, pedfile) != 1) {
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
      sprintf(logbuf, "Error: Invalid .bed file size (expected %llu bytes).\n", ((unsigned long long)unfiltered_indiv_ct4) * unfiltered_marker_ct);
      logprintb();
      goto wdist_ret_INVALID_FORMAT;
    } else {
      // pre-0.99, no magic number, indiv-major
      bed_offset = 2;
    }
    if (bed_offset == 2) {
      strcpy(outname_end, ".bed.tmp"); // not really temporary
      logprint("Individual-major .bed file detected.  Transposing to SNP-major form.\n");
      fclose(pedfile);
      retval = indiv_major_to_snp_major(pedname, outname, &outfile, unfiltered_marker_ct);
      if (retval) {
	goto wdist_ret_2;
      }
      strcpy(pedname, outname);
      if (fopen_checked(&pedfile, pedname, "rb")) {
	goto wdist_ret_OPEN_FAIL;
      }
      bed_offset = 3;
    }
  }
  if (mind_thresh < 1.0) {
    retval = mind_filter(pedfile, mind_thresh, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, bed_offset, line_mids, pedbuflen, missing_geno);
    if (retval) {
      goto wdist_ret_2;
    }
  }
  g_indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  if (!g_indiv_ct) {
    logprint("Error: Nobody passes QC.\n");
    goto wdist_ret_INVALID_FORMAT;
  }
  // text HWE waits until after external freq file has had a chance to modify
  // major/minor allele
  nonfounders = (nonfounders || (!fam_col_34));
  wt_needed = distance_wt_req(calculation_type) && (!distance_flat_missing);
  retval = calc_freqs_and_binary_hwe(pedfile, unfiltered_marker_ct, marker_exclude, unfiltered_indiv_ct, indiv_exclude, founder_info, nonfounders, maf_succ, set_allele_freqs, &marker_alleles, &marker_reverse, &marker_allele_cts, &missing_cts, bed_offset, line_mids, pedbuflen, (unsigned char)missing_geno, hwe_all, g_pheno_c, &hwe_lls, &hwe_lhs, &hwe_hhs, &hwe_ll_allfs, &hwe_lh_allfs, &hwe_hh_allfs, &ll_cts, &lh_cts, &hh_cts, wt_needed, &marker_weights_base, &g_marker_weights);
  if (retval) {
    goto wdist_ret_2;
  }
  if (freqname) {
    retval = read_external_freqs(freqname, &freqfile, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, chrom_info_ptr, marker_alleles, marker_allele_cts, set_allele_freqs, binary_files, maf_succ, missing_geno, exponent, wt_needed, g_marker_weights);
    if (retval) {
      goto wdist_ret_2;
    }
  }
  if (!binary_files) {
    text_normalize_marker_alleles(marker_alleles, unfiltered_marker_ct, marker_exclude);
    // now marker_alleles[2 * ii + 1] is the set allele, just like binary case.
    // however, marker COUNTS are not normalized.
    retval = text_load_hwe(pedfile, unfiltered_marker_ct, marker_exclude, unfiltered_indiv_ct, indiv_exclude, founder_info, nonfounders, marker_alleles, line_mids, pedbuflen, hwe_all, g_pheno_c, hwe_lls, hwe_lhs, hwe_hhs, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, ll_cts, lh_cts, hh_cts);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (!keep_allele_order) {
    if (binary_files) {
      // different tiebreaking behavior in PLINK, for whatever reason
      calc_marker_reverse_bin(marker_reverse, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, set_allele_freqs);
    } else {
      calc_marker_reverse_text(marker_reverse, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, set_allele_freqs);
    }
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
    retval = write_freqs(&outfile, outname, plink_maxsnp, unfiltered_marker_ct, marker_exclude, set_allele_freqs, chrom_info_ptr, marker_ids, max_marker_id_len, marker_alleles, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, freq_counts, freqx, missing_geno, marker_reverse);
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_FREQ))))) {
      goto wdist_ret_2;
    }
  }
  if (geno_thresh < 1.0) {
    if (binary_files) {
      // would prefer to use missing_cts as last parameter, but that doesn't
      // play well with --freq
      uii = binary_geno_filter(geno_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, g_indiv_ct, marker_allele_cts);
    } else {
      uii = text_geno_filter(geno_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, g_indiv_ct, marker_allele_cts, missing_cts);
    }
    sprintf(logbuf, "%u SNP%s removed due to missing genotype data (--geno).\n", uii, (uii == 1)? "" : "s");
    logprintb();
  }
  wkspace_reset(marker_allele_cts);
  marker_allele_cts = NULL;
  // missing_cts = NULL;
  if (hwe_thresh > 0.0) {
    het_probs = (double*)malloc((g_indiv_ct / 2 + 1) * sizeof(double));
    if (!het_probs) {
      goto wdist_ret_NOMEM;
    }
    enforce_hwe_threshold(hwe_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, hwe_lls, hwe_lhs, hwe_hhs, hwe_all || (!g_pheno_c), hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs);
    free(het_probs);
    het_probs = NULL;
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
    calc_marker_weights(exponent, unfiltered_marker_ct, marker_exclude, ll_cts, lh_cts, hh_cts, g_marker_weights);
  }
  wkspace_reset(hwe_lls);

  if (!binary_files) {
    retval = text_to_bed(&bedtmpfile, &famtmpfile, &bimtmpfile, &outfile, outname, outname_end, &pedfile, line_locs, line_mids, pedbuflen, &famfile, mapfile, map_cols, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, marker_alleles, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_info, g_pheno_c, g_pheno_d, missing_phenod, output_missing_pheno, max_marker_id_len, map_is_unsorted, chrom_info_ptr->species);
    if (retval || (!(calculation_type & (~(CALC_FREQ | CALC_MAKE_BED))))) {
      goto wdist_ret_2;
    }
    collapse_arr(marker_alleles, 2, marker_exclude, unfiltered_marker_ct);
    sscanf(output_missing_pheno, "%lg", &missing_phenod);
    // if this becomes much more of a maintenance nightmare, consider exiting
    // function and reloading from .bed from scratch
    if (g_pheno_c) {
      collapse_arr(g_pheno_c, sizeof(char), indiv_exclude, unfiltered_indiv_ct);
    } else if (g_pheno_d) {
      collapse_arr((char*)g_pheno_d, sizeof(double), indiv_exclude, unfiltered_indiv_ct);
    }
    if (sex_info) {
      collapse_arr((char*)sex_info, sizeof(char), indiv_exclude, unfiltered_indiv_ct);
    }
    collapse_bitarr(marker_reverse, marker_exclude, unfiltered_marker_ct);
    if ((calculation_type & CALC_WRITE_SNPLIST) || ((calculation_type & CALC_RECODE) && (recode_modifier & RECODE_LGEN))) {
      collapse_arr(marker_ids, max_marker_id_len, marker_exclude, unfiltered_marker_ct);
    }
    if (calculation_type & (CALC_WRITE_SNPLIST | CALC_GENOME | CALC_LD_PRUNE)) {
      collapse_chrom_marker_idxs(chrom_info_ptr, marker_exclude, unfiltered_marker_ct);
      if (calculation_type & (CALC_GENOME | CALC_LD_PRUNE)) {
	collapse_arr((char*)marker_pos, sizeof(int), marker_exclude, unfiltered_marker_ct);
      }
    }
    collapse_arr(person_ids, max_person_id_len, indiv_exclude, unfiltered_indiv_ct);
    if (fam_col_34) {
      collapse_arr(paternal_ids, max_paternal_id_len, indiv_exclude, unfiltered_indiv_ct);
      collapse_arr(maternal_ids, max_maternal_id_len, indiv_exclude, unfiltered_indiv_ct);
      collapse_bitarr(founder_info, indiv_exclude, unfiltered_indiv_ct);
    }
    if (wt_needed) {
      collapse_arr((char*)g_marker_weights, sizeof(double), marker_exclude, unfiltered_marker_ct);
    }
    collapse_arr((char*)set_allele_freqs, sizeof(double), marker_exclude, unfiltered_marker_ct);
    unfiltered_marker_ct -= marker_exclude_ct;
    marker_exclude_ct = 0;
    fill_ulong_zero(marker_exclude, (unfiltered_marker_ct + (BITCT - 1)) / BITCT);
    unfiltered_indiv_ct -= indiv_exclude_ct;
    unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
    indiv_exclude_ct = 0;
    fill_ulong_zero(indiv_exclude, (unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  }
  if (!allow_no_sex) {
    ii = indiv_exclude_ct;
    exclude_ambiguous_sex(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, sex_info);
    ii = indiv_exclude_ct - ii;
    if (ii) {
      sprintf(logbuf, "%d individual%s with unknown sex removed (prevent with --allow-no-sex).\n", ii, (ii == 1)? "" : "s");
      logprintb();
    }
    g_indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
    if (!g_indiv_ct) {
      logprint("Error: Nobody remaining.\n");
      goto wdist_ret_INVALID_FORMAT;
    }
  }
  sprintf(logbuf, "%d markers and %d people pass filters and QC%s.\n", marker_ct, g_indiv_ct, (calculation_type & CALC_REL_CUTOFF)? " (before --rel-cutoff)": "");
  logprintb();

  if (parallel_tot > g_indiv_ct / 2) {
    sprintf(logbuf, "Error: Too many --parallel jobs (maximum %d/2 = %d).\n", g_indiv_ct, g_indiv_ct / 2);
    goto wdist_ret_INVALID_CMDLINE_2;
  }
  if (g_thread_ct > 1) {
    if (calculation_type & (CALC_RELATIONSHIP | CALC_IBC | CALC_GDISTANCE_MASK | CALC_GROUPDIST | CALC_REGRESS_DISTANCE | CALC_GENOME | CALC_REGRESS_REL | CALC_UNRELATED_HERITABILITY)) {
      sprintf(logbuf, "Using %d threads (change this with --threads).\n", g_thread_ct);
      logprintb();
    } else {
      logprint("Using 1 thread (no multithreaded calculations invoked).\n");
    }
  }

  if (refalleles) {
    retval = load_ref_alleles(refalleles, unfiltered_marker_ct, marker_alleles, marker_reverse, marker_ids, max_marker_id_len);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_WRITE_SNPLIST) {
    retval = write_snplist(&outfile, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_len);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (binary_files && (calculation_type & CALC_MAKE_BED)) {
    retval = make_bed(pedfile, bed_offset, mapfile, map_cols, &bedtmpfile, &famtmpfile, &bimtmpfile, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_reverse, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_info, g_pheno_c, g_pheno_d, missing_phenod, output_missing_pheno, max_marker_id_len, map_is_unsorted, indiv_sort_map, chrom_info_ptr->species);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_RECODE) {
    retval = recode(recode_modifier, pedfile, bed_offset, famfile, mapfile, &outfile, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, marker_ids, max_marker_id_len, marker_alleles, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_info, g_pheno_c, g_pheno_d, missing_phenod, output_missing_geno, output_missing_pheno, set_allele_freqs, binary_files);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_LD_PRUNE) {
    retval = ld_prune(pedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_len, chrom_info_ptr, set_allele_freqs, marker_pos, unfiltered_indiv_ct, indiv_exclude, ld_window_size, ld_window_kb, ld_window_incr, ld_last_param, outname, outname_end, calculation_type);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_REGRESS_PCS) {
    // do this before marker_alleles is overwritten in memory...
    retval = calc_regress_pcs(evecname, regress_pcs_normalize_pheno, regress_pcs_sex_specific, regress_pcs_clip, max_pcs, pedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_len, marker_alleles, chrom_info_ptr, marker_pos, g_indiv_ct, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len, sex_info, g_pheno_d, missing_phenod, &outfile, outname, outname_end);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (wt_needed) {
    // N.B. marker_weights is currently on top of the stack
    wkspace_reset(marker_weights_base);
    // normalize included marker weights to add to 2^32 - 1
    dxx = 0.0;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      if (is_set(marker_exclude, marker_uidx)) {
	continue;
      }
      dxx += g_marker_weights[marker_uidx];
    }
    dxx = 4294967295.0 / dxx;
    marker_idx = 0;
    g_marker_weights_i = (unsigned int*)wkspace_alloc(marker_ct * sizeof(int));
    for (marker_uidx = 0; marker_idx < marker_ct; marker_uidx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      g_marker_weights_i[marker_idx++] = (unsigned int)(g_marker_weights[marker_uidx] * dxx + 0.5);
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
    wkspace_mark2 = wkspace_base;
  }

  if (relationship_or_ibc_req(calculation_type)) {
    g_indiv_missing_unwt = (unsigned int*)wkspace_alloc(g_indiv_ct * sizeof(int));
    if (!g_indiv_missing_unwt) {
      goto wdist_ret_NOMEM;
    }
    fill_int_zero((int*)g_indiv_missing_unwt, g_indiv_ct);
;
    triangle_fill(g_thread_start, g_indiv_ct, g_thread_ct, parallel_idx, parallel_tot, 1, 1);
    if (relationship_req(calculation_type)) {
      llxx = g_thread_start[g_thread_ct];
      llxx = ((llxx * (llxx - 1)) - (long long)g_thread_start[0] * (g_thread_start[0] - 1)) / 2;
      dists_alloc = llxx * sizeof(double);
      if (!(calculation_type & CALC_UNRELATED_HERITABILITY)) {
        // if the memory isn't needed for CALC_UNRELATED_HERITABILITY,
	// positioning the missingness matrix here will let us avoid
	// recalculating it if --distance-matrix or --matrix is requested
        g_missing_dbl_excluded = (unsigned int*)wkspace_alloc(llxx * sizeof(int));
        if (!g_missing_dbl_excluded) {
          goto wdist_ret_NOMEM;
        }
        fill_int_zero((int*)g_missing_dbl_excluded, llxx);
      }
      g_rel_dists = (double*)wkspace_alloc(dists_alloc);
      if (!g_rel_dists) {
	goto wdist_ret_NOMEM;
      }
      fill_double_zero(g_rel_dists, llxx);
    }
    if (calculation_type & CALC_IBC) {
      ii = g_indiv_ct * 3;
    } else {
      ii = g_indiv_ct;
    }
    rel_ibc = (double*)wkspace_alloc(ii * sizeof(double));
    if (!rel_ibc) {
      goto wdist_ret_NOMEM;
    }
    fill_double_zero(rel_ibc, ii);
    wkspace_mark = wkspace_base;
    if (relationship_req(calculation_type) && (!g_missing_dbl_excluded)) {
      g_missing_dbl_excluded = (unsigned int*)wkspace_alloc(llxx * sizeof(int));
      if (!g_missing_dbl_excluded) {
	goto wdist_ret_NOMEM;
      }
      fill_int_zero((int*)g_missing_dbl_excluded, llxx);
    }
    fseeko(pedfile, bed_offset, SEEK_SET);
    ii = 0;
    marker_idx = 0;
    g_geno = wkspace_alloc(g_indiv_ct * sizeof(long));
    if (!g_geno) {
      goto wdist_ret_NOMEM;
    }
    g_mmasks = (unsigned long*)wkspace_alloc(g_indiv_ct * sizeof(long));
    if (!g_mmasks) {
      goto wdist_ret_NOMEM;
    }
    gptr = wkspace_alloc(MULTIPLEX_REL * unfiltered_indiv_ct4);
    if (!gptr) {
      goto wdist_ret_NOMEM;
    }
    g_masks = (unsigned long*)wkspace_alloc(g_indiv_ct * sizeof(long));
    if (!g_masks) {
      goto wdist_ret_NOMEM;
    }

    // See comments at the beginning of this file, and those in the main
    // CALC_DISTANCE loop.  The main difference between this calculation and
    // the (nonzero exponent) distance calculation is that we have to pad
    // each marker to 3 bits and use + instead of XOR to distinguish the
    // cases.
    while (marker_idx < marker_ct) {
      jj = 0;
      while ((jj < MULTIPLEX_REL) && (marker_idx < marker_ct)) {
	ulii = 0;
	while (is_set(marker_exclude, ii)) {
	  ii++;
	  ulii++;
	}
	if (ulii) {
	  fseeko(pedfile, ulii * unfiltered_indiv_ct4, SEEK_CUR);
	}
	if (fread(&(gptr[jj * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	  goto wdist_ret_READ_FAIL;
	}
	set_allele_freq_buf[jj] = set_allele_freqs[ii];
	ii++;
	jj++;
	marker_idx++;
      }
      if (jj < MULTIPLEX_REL) {
	memset(&(gptr[jj * unfiltered_indiv_ct4]), 0, (MULTIPLEX_REL - jj) * unfiltered_indiv_ct4);
	fill_double_zero(&(set_allele_freq_buf[jj]), MULTIPLEX_REL - jj);
      }
      fill_ulong_zero(g_mmasks, g_indiv_ct);

      for (nn = 0; nn < jj; nn += MULTIPLEX_REL / 3) {
	fill_ulong_zero(g_masks, g_indiv_ct);
	indiv_idx = 0;
	glptr2 = (unsigned long*)g_geno;
	for (qq = 0; indiv_idx < g_indiv_ct; qq++) {
	  while (is_set(indiv_exclude, qq)) {
	    qq++;
	  }
	  ulii = 0;
	  gptr2 = &(gptr[qq / 4 + nn * unfiltered_indiv_ct4]);
	  kk = (qq % 4) * 2;
	  for (mm = 0; mm < (MULTIPLEX_REL / 3); mm++) {
	    uljj = (gptr2[mm * unfiltered_indiv_ct4] >> kk) & 3;
	    if (uljj == 1) {
	      g_masks[indiv_idx] |= 7LU << (mm * 3);
	      g_mmasks[indiv_idx] |= 1LU << (nn + mm);
	      g_indiv_missing_unwt[indiv_idx] += 1;
	    }
	    ulii |= uljj << (mm * 3);
	  }
	  *glptr2++ = ulii;
	  indiv_idx++;
	}
	if (calculation_type & CALC_IBC) {
	  for (oo = 0; oo < 3; oo++) {
	    update_rel_ibc(&(rel_ibc[oo * g_indiv_ct]), (unsigned long*)g_geno, &(set_allele_freq_buf[nn]), oo);
	  }
	} else {
	  update_rel_ibc(rel_ibc, (unsigned long*)g_geno, &(set_allele_freq_buf[nn]), ibc_type);
	}
	if (relationship_req(calculation_type)) {
	  fill_weights_r(g_weights, &(set_allele_freq_buf[nn]), var_std);
	  for (ulii = 1; ulii < g_thread_ct; ulii++) {
	    if (pthread_create(&(threads[ulii - 1]), NULL, &calc_rel_thread, (void*)ulii)) {
	      goto wdist_ret_THREAD_CREATE_FAIL;
	    }
	  }
	  incr_dists_r(g_rel_dists, (unsigned long*)g_geno, g_masks, 0, g_weights);
	  for (uii = 0; uii < g_thread_ct - 1; uii++) {
	    pthread_join(threads[uii], NULL);
	  }
	}
      }
      if (relationship_req(calculation_type)) {
	for (ulii = 1; ulii < g_thread_ct; ulii++) {
	  if (pthread_create(&(threads[ulii - 1]), NULL, &calc_missing_thread, (void*)ulii)) {
	    goto wdist_ret_THREAD_CREATE_FAIL;
	  }
	}
	incr_dists_rm(g_missing_dbl_excluded, 0, g_thread_start);
	for (uii = 0; uii < g_thread_ct - 1; uii++) {
	  pthread_join(threads[uii], NULL);
	}
      }
      printf("\r%d markers complete.", marker_idx);
      fflush(stdout);
    }
    if (relationship_req(calculation_type)) {
      putchar('\r');
      logprint("Relationship matrix calculation complete.\n");
      dist_ptr = g_rel_dists;
    } else {
      putchar('\n');
    }
    dptr2 = rel_ibc;
    if (calculation_type & CALC_IBC) {
      dptr3 = &(rel_ibc[g_indiv_ct]);
      dptr4 = &(rel_ibc[g_indiv_ct * 2]);
    }
    giptr2 = g_missing_dbl_excluded;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      uii = marker_ct - g_indiv_missing_unwt[indiv_idx];
      if ((indiv_idx >= g_thread_start[0]) && (indiv_idx < g_thread_start[g_thread_ct])) {
	if (relationship_req(calculation_type)) {
	  giptr = g_indiv_missing_unwt;
	  for (ujj = 0; ujj < indiv_idx; ujj++) {
	    *dist_ptr /= uii - (*giptr++) + (*giptr2++);
	    dist_ptr++;
	  }
	}
      }
      if (calculation_type & CALC_IBC) {
	*dptr2 /= uii;
	dptr2++;
	*dptr3 /= uii;
	dptr3++;
	*dptr4 /= uii;
	dptr4++;
      } else {
	*dptr2 /= uii;
	dptr2++;
      }
    }
    if (calculation_type & CALC_REL_CUTOFF) {
      retval = do_rel_cutoff(calculation_type, rel_cutoff, rel_ibc, indiv_exclude, &indiv_exclude_ct, outname, outname_end, unfiltered_indiv_ct, person_ids, max_person_id_len);
      if (retval) {
	goto wdist_ret_2;
      }
    }

    if (calculation_type & CALC_IBC) {
      strcpy(outname_end, ".ibc");
      if (fopen_checked(&outfile, outname, "w")) {
	goto wdist_ret_OPEN_FAIL;
      }
      dptr2 = rel_ibc;
      dptr3 = &(rel_ibc[g_indiv_ct]);
      dptr4 = &(rel_ibc[g_indiv_ct * 2]);
      if (fputs_checked("FID\tIID\tNOMISS\tFhat1\tFhat2\tFhat3\n", outfile)) {
	goto wdist_ret_WRITE_FAIL;
      }
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	if (fprintf(outfile, "%d\t%d\t%d\t%g\t%g\t%g\n", indiv_idx + 1, indiv_idx + 1, marker_ct - g_indiv_missing_unwt[indiv_idx], *dptr3++ - 1.0, *dptr4++ - 1.0, *dptr2++ - 1.0) < 0) {
	  goto wdist_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&outfile)) {
	goto wdist_ret_WRITE_FAIL;
      }
      sprintf(logbuf, "%s written.\n", outname);
      logprintb();
    }
    if (calculation_type & CALC_RELATIONSHIP) {
      mm = 1;
      nn = rel_calc_type & REL_CALC_SHAPEMASK;
      if (parallel_tot == 1) {
	oo = g_indiv_ct;
	// nasty rel-cutoff bug
      } else {
	oo = g_thread_start[g_thread_ct];
      }
      pp = g_thread_start[0];
      if (pp == 1) {
	pp = 0;
      }
      if (calculation_type & CALC_IBC) {
	dptr2 = &(rel_ibc[ibc_type * g_indiv_ct + pp]);
      } else {
	dptr2 = &(rel_ibc[pp]);
      }
      llxx = ((long long)pp * (pp - 1)) / 2;
      llyy = (((long long)oo * (oo + 1)) / 2) - llxx;
      if (calculation_type & REL_CALC_BIN) {
	if (nn == REL_CALC_SQ0) {
	  fill_double_zero((double*)g_geno, g_indiv_ct - 1);
	}
	strcpy(outname_end, ".rel.bin");
	if (parallel_tot > 1) {
	  sprintf(&(outname_end[8]), ".%d", parallel_idx + 1);
	}
	if (fopen_checked(&outfile, outname, "wb")) {
	  goto wdist_ret_OPEN_FAIL;
	}
	for (ii = pp; ii < oo; ii++) {
	  if (fwrite_checked(&(g_rel_dists[((long long)ii * (ii - 1)) / 2 - llxx]), ii * sizeof(double), outfile)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  if (fwrite_checked(dptr2++, sizeof(double), outfile)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  if (nn == REL_CALC_TRI) {
	    if ((((long long)ii + 1) * (ii + 2) / 2 - llxx) * 100 >= llyy * mm) {
	      mm = (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100) / llyy;
	      printf("\rWriting... %d%%", mm++);
	      fflush(stdout);
	    }
	  } else {
	    if (nn == REL_CALC_SQ0) {
	      if (fwrite_checked(g_geno, (g_indiv_ct - ii - 1) * sizeof(double), outfile)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    } else {
	      for (uii = ii + 1; uii < g_indiv_ct; uii++) {
		if (fwrite_checked(&(g_rel_dists[((unsigned long)uii * (uii - 1) / 2) + ii - llxx]), sizeof(double), outfile)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if ((ii + 1 - pp) * 100 >= mm * (oo - pp)) {
	      mm = ((ii + 1 - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", mm++);
	      fflush(stdout);
	    }
	  }
	}
	if (fclose_null(&outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
      } else if (rel_calc_type & REL_CALC_GRM) {
	giptr2 = g_missing_dbl_excluded;
	if (rel_calc_type & REL_CALC_GZ) {
	  if (parallel_tot > 1) {
	    sprintf(outname_end, ".grm.%d.gz", parallel_idx + 1);
	  } else {
	    strcpy(outname_end, ".grm.gz");
	  }
	  if (gzopen_checked(&gz_outfile, outname, "wb")) {
	    goto wdist_ret_OPEN_FAIL;
	  }
	  dist_ptr = g_rel_dists;
	  for (ii = pp; ii < oo; ii++) {
	    if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * mm) {
	      mm = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
	      printf("\rWriting... %d%%", mm++);
	      fflush(stdout);
	    }
	    qq = marker_ct - g_indiv_missing_unwt[ii];
	    giptr = g_indiv_missing_unwt;
	    for (jj = 0; jj < ii; jj += 1) {
	      kk = qq - (*giptr++) + (*giptr2++);
	      if (!gzprintf(gz_outfile, "%d\t%d\t%d\t%e\n", ii + 1, jj + 1, kk, *dist_ptr++)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (!gzprintf(gz_outfile, "%d\t%d\t%d\t%e\n", ii + 1, jj + 1, qq, *dptr2++)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  gzclose(gz_outfile);
	  gz_outfile = NULL;
	} else {
	  strcpy(outname_end, ".grm");
	  if (parallel_tot > 1) {
	    sprintf(&(outname_end[4]), ".%d", parallel_idx + 1);
	  }
	  if (fopen_checked(&outfile, outname, "w")) {
	    goto wdist_ret_OPEN_FAIL;
	  }
	  dist_ptr = g_rel_dists;
	  for (ii = pp; ii < oo; ii++) {
	    if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * mm) {
	      mm = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
	      printf("\rWriting... %d%%", mm++);
	      fflush(stdout);
	    }
	    qq = marker_ct - g_indiv_missing_unwt[ii];
	    giptr = g_indiv_missing_unwt;
	    for (jj = 0; jj < ii; jj++) {
	      kk = qq - (*giptr++) + (*giptr2++);
	      if (fprintf(outfile, "%d\t%d\t%d\t%e\n", ii + 1, jj + 1, kk, *dist_ptr++) < 0) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (fprintf(outfile, "%d\t%d\t%d\t%e\n", ii + 1, jj + 1, qq, *dptr2++) < 0) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (fclose_null(&outfile)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
      } else {
	if (nn == REL_CALC_SQ0) {
	  cptr2 = (char*)(&ulii);
	  for (uii = 0; uii < sizeof(long); uii += 2) {
	    cptr2[uii] = '\t';
	    cptr2[uii + 1] = '0';
	  }
	  ii = (g_indiv_ct * 2 + sizeof(long) - 4) / sizeof(long);
	  glptr2 = (unsigned long*)g_geno;
	  for (jj = 0; jj < ii; jj++) {
	    *glptr2++ = ulii;
	  }
	}
	if (rel_calc_type & REL_CALC_GZ) {
	  if (parallel_tot > 1) {
	    sprintf(outname_end, ".rel.%d.gz", parallel_idx + 1);
	  } else {
	    strcpy(outname_end, ".rel.gz");
	  }
	  if (gzopen_checked(&gz_outfile, outname, "wb")) {
	    goto wdist_ret_OPEN_FAIL;
	  }
	  dist_ptr = g_rel_dists;
	  if (pp) {
	    ii = pp;
	  } else {
	    if (!gzprintf(gz_outfile, "%g", *dptr2++)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (nn == REL_CALC_SQ0) {
	      if (gzwrite_checked(gz_outfile, g_geno, (g_indiv_ct - 1) * 2)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    } else if (nn == REL_CALC_SQ) {
	      // parallel_tot must be 1 for SQ shape
	      for (indiv_idx = 1; indiv_idx < g_indiv_ct; indiv_idx++) {
		if (!gzprintf(gz_outfile, "\t%g", g_rel_dists[((unsigned long)indiv_idx * (indiv_idx - 1)) / 2])) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if (gzwrite_checked(gz_outfile, "\n", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    ii = 1;
	  }
	  for (; ii < oo; ii++) {
	    if (!gzprintf(gz_outfile, "%g", *dist_ptr++)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    for (jj = 1; jj < ii; jj++) {
	      if (!gzprintf(gz_outfile, "\t%g", *dist_ptr++)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (!gzprintf(gz_outfile, "\t%g", *dptr2++)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (nn == REL_CALC_SQ0) {
	      if (gzwrite_checked(gz_outfile, g_geno, (g_indiv_ct - jj - 1) * 2)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      if ((ii + 1 - pp) * 100 >= mm * oo) {
		mm = ((ii + 1 - pp) * 100) / oo;
		printf("\rWriting... %d%%", mm++);
		fflush(stdout);
	      }
	    } else {
	      if (nn == REL_CALC_SQ) {
		for (ulii = ii + 1; ulii < g_indiv_ct; ulii++) {
		  if (!gzprintf(gz_outfile, "\t%g", g_rel_dists[((ulii * (ulii - 1)) / 2) + ii])) {
		    goto wdist_ret_WRITE_FAIL;
		  }
		}
	      }
	      if (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100 >= llyy * mm) {
		mm = (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100) / llyy;
		printf("\rWriting... %d%%", mm++);
		fflush(stdout);
	      }
	    }
	    if (gzwrite_checked(gz_outfile, "\n", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  gzclose(gz_outfile);
	  gz_outfile = NULL;
	} else {
	  strcpy(outname_end, ".rel");
	  if (parallel_tot > 1) {
	    sprintf(&(outname_end[4]), ".%d", parallel_idx + 1);
	  }
	  if (fopen_checked(&outfile, outname, "w")) {
	    goto wdist_ret_OPEN_FAIL;
	  }
	  dist_ptr = g_rel_dists;
	  if (pp) {
	    ii = pp;
	  } else {
	    if (fprintf(outfile, "%g", *dptr2++) < 0) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (nn == REL_CALC_SQ0) {
	      if (fwrite_checked(g_geno, (g_indiv_ct - 1) * 2, outfile)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    } else if (nn == REL_CALC_SQ) {
	      for (indiv_idx = 1; indiv_idx < g_indiv_ct; indiv_idx++) {
		if (fprintf(outfile, "\t%g", g_rel_dists[((unsigned long)indiv_idx * (indiv_idx - 1)) / 2]) < 0) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if (fwrite_checked("\n", 1, outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    ii = 1;
	  }
	  for (; ii < oo; ii++) {
	    if (fprintf(outfile, "%g", *dist_ptr++) < 0) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    for (jj = 1; jj < ii; jj++) {
	      if (fprintf(outfile, "\t%g", *dist_ptr++) < 0) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (fprintf(outfile, "\t%g", *dptr2++) < 0) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (nn == REL_CALC_SQ0) {
	      if (fwrite_checked(g_geno, (g_indiv_ct - ii - 1) * 2, outfile)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      if ((ii + 1 - pp) * 100LL >= (long long)mm * (oo - pp)) {
		mm = ((ii + 1 - pp) * 100LL) / (oo - pp);
		printf("\rWriting... %d%%", mm++);
		fflush(stdout);
	      }
	    } else {
	      if (nn == REL_CALC_SQ) {
		for (ulii = ii + 1; ulii < g_indiv_ct; ulii++) {
		  if (fprintf(outfile, "\t%g", g_rel_dists[((ulii * (ulii - 1)) / 2) + ii]) < 0) {
		    goto wdist_ret_WRITE_FAIL;
		  }
		}
	      }
	      if (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100 >= llyy * mm) {
		mm = (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100) / llyy;
		printf("\rWriting... %d%%", mm++);
		fflush(stdout);
	      }
	    }
	    if (fwrite_checked("\n", 1, outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (fclose_null(&outfile)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
      }
      putchar('\r');
      sprintf(logbuf, "Relationship matrix written to %s.\n", outname);
      logprintb();
      if (!parallel_idx) {
	strcpy(&(outname_end[4]), ".id");
	retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
	if (retval) {
	  goto wdist_ret_2;
	}
      }
    }
    wkspace_reset(wkspace_mark);

    if (calculation_type & CALC_REGRESS_REL) {
      retval = regress_rel_main(indiv_exclude, g_indiv_ct, regress_rel_iters, regress_rel_d, threads);
      if (retval) {
	goto wdist_ret_2;
      }
    }
#ifndef NOLAPACK
    if (calculation_type & CALC_UNRELATED_HERITABILITY) {
      g_missing_dbl_excluded = NULL;
      ulii = g_indiv_ct;
      ulii = CACHEALIGN_DBL(ulii * ulii);
      dptr4 = &(g_rel_dists[ulii]);
      ulii = ulii * 3 + CACHEALIGN_DBL(g_indiv_ct) * 3;
      wkspace_reset(g_rel_dists);
      g_rel_dists = (double*)wkspace_alloc(ulii * sizeof(double));
      if (!g_rel_dists) {
        goto wdist_ret_NOMEM;
      }
      dptr2 = &(g_rel_dists[ulii - CACHEALIGN_DBL(g_indiv_ct)]);
      collapse_copy_phenod(dptr2, g_pheno_d, indiv_exclude, unfiltered_indiv_ct);
      dxx = 0.0;
      dyy = 0.0;
      dptr3 = dptr2;
      dist_ptr = &(dptr2[g_indiv_ct]);
      while (dptr3 < dist_ptr) {
        dzz = *dptr3++;
        dxx += dzz;
        dyy += dzz * dzz;
      }
      dxx /= (double)g_indiv_ct;
      dxx = 1 / sqrt((dyy / (double)g_indiv_ct) - dxx * dxx);
      dptr3 = dptr2;
      while (dptr3 < dist_ptr) {
        *dptr3 *= dxx;
        dptr3++;
      }
      if (calculation_type & CALC_IBC) {
        dptr3 = &(rel_ibc[ibc_type * g_indiv_ct]);
      } else {
        dptr3 = rel_ibc;
      }
      for (ulii = 0; ulii < g_indiv_ct; ulii++) {
	memcpy(&(dptr4[ulii * g_indiv_ct]), &(g_rel_dists[(ulii * (ulii - 1)) / 2]), ulii * sizeof(double));
        dptr4[ulii * (g_indiv_ct + 1)] = *dptr3++;
        for (uljj = ulii + 1; uljj < g_indiv_ct; uljj++) {
          dptr4[ulii * g_indiv_ct + uljj] = g_rel_dists[(uljj * (uljj - 1)) / 2 + ulii];
        }
      }
      reml_em_one_trait(g_rel_dists, dptr2, &unrelated_herit_covg, &unrelated_herit_covr, unrelated_herit_tol, calculation_type & CALC_UNRELATED_HERITABILITY_STRICT);
      sprintf(logbuf, "h^2 estimate: %g\n", unrelated_herit_covg);
      logprintb();
    }
#endif
    if ((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX | CALC_GENOME)) && (!(calculation_type & CALC_REL_CUTOFF))) {
      wkspace_reset(g_rel_dists);
    } else {
      wkspace_reset(g_indiv_missing_unwt);
      g_indiv_missing_unwt = NULL;
      g_missing_dbl_excluded = NULL;
    }
  }

  if (calculation_type & CALC_REGRESS_PCS_DISTANCE) {
    logprint("Error: --regress-pcs-distance is currently under development.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto wdist_ret_1;
  } else if (distance_req(calculation_type)) {
    triangle_fill(g_thread_start, g_indiv_ct, g_thread_ct, parallel_idx, parallel_tot, 1, 1);
    llxx = g_thread_start[g_thread_ct];
    llxx = ((llxx * (llxx - 1)) - (long long)g_thread_start[0] * (g_thread_start[0] - 1)) / 2;
    dists_alloc = llxx * sizeof(double);
    // Additional + CACHELINE is to fix aliasing bug that shows up with -O2 in
    // some cases.
    // The g_missing_dbl_excluded check is to avoid recalculation, if a
    // relationship matrix was already calculated, the g_missing_dbl_excluded
    // matrix was not overwritten, and no --rel-cutoff filtering was performed
    // to make the results obsolete.  It's unimportant to keep this
    // optimization around if it ever complicates maintenance.
    if (((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX )) || distance_flat_missing) && (!g_missing_dbl_excluded)) {
      g_missing_dbl_excluded = (unsigned int*)wkspace_alloc(llxx * sizeof(int));
      if (!g_missing_dbl_excluded) {
        goto wdist_ret_NOMEM;
      }
      fill_int_zero((int*)g_missing_dbl_excluded, llxx);
      unwt_needed = 1;
      if (!g_indiv_missing_unwt) {
        g_indiv_missing_unwt = (unsigned int*)wkspace_alloc(g_indiv_ct * sizeof(int));
        if (!g_indiv_missing_unwt) {
          goto wdist_ret_NOMEM;
        }
	unwt_needed_full = 1;
        fill_int_zero((int*)g_indiv_missing_unwt, g_indiv_ct);
      }
    }
    g_dists = (double*)wkspace_alloc(dists_alloc + CACHELINE);
    if (!g_dists) {
      goto wdist_ret_NOMEM;
    }
    wkspace_mark = wkspace_base;
    if (wt_needed) {
      g_missing_tot_weights = (unsigned int*)wkspace_alloc(llxx * sizeof(int));
      if (!g_missing_tot_weights) {
	goto wdist_ret_NOMEM;
      }
      fill_int_zero((int*)g_missing_tot_weights, llxx);
      g_indiv_missing = (unsigned int*)wkspace_alloc(g_indiv_ct * sizeof(int));
      if (!g_indiv_missing) {
	goto wdist_ret_NOMEM;
      }
      fill_int_zero((int*)g_indiv_missing, g_indiv_ct);
    }

    if (exp0) {
      g_idists = (int*)(((char*)wkspace_mark) - CACHEALIGN(llxx * sizeof(int)));
      fill_int_zero(g_idists, llxx);
      g_masks = (unsigned long*)wkspace_alloc(g_indiv_ct * (MULTIPLEX_2DIST / 8));
    } else {
      fill_double_zero(g_dists, llxx);
      g_masks = (unsigned long*)wkspace_alloc(g_indiv_ct * sizeof(long));
    }
    if (!g_masks) {
      goto wdist_ret_NOMEM;
    }
    g_mmasks = (unsigned long*)wkspace_alloc(g_indiv_ct * sizeof(long));
    if (!g_mmasks) {
      goto wdist_ret_NOMEM;
    }

    if (exp0) {
      multiplex = MULTIPLEX_DIST;
      g_geno = wkspace_alloc(g_indiv_ct * (MULTIPLEX_2DIST / 8));
    } else {
      multiplex = MULTIPLEX_DIST_EXP;
      g_geno = wkspace_alloc(g_indiv_ct * sizeof(long));
    }
    if (!g_geno) {
      goto wdist_ret_NOMEM;
    }

    pedbuf = (unsigned char*)malloc(multiplex * unfiltered_indiv_ct4 * sizeof(char));
    if (!pedbuf) {
      goto wdist_ret_NOMEM;
    }
    fseeko(pedfile, bed_offset, SEEK_SET);
    ii = 0; // current SNP index
    marker_idx = 0; // after subtracting out excluded
    while (marker_idx < marker_ct) {
      for (jj = 0; jj < multiplex; jj++) {
	set_allele_freq_buf[jj] = 0.5;
      }
      fill_int_zero((int*)wtbuf, multiplex);
      jj = 0; // actual SNPs read

      // For each pair (g_j, g_k) of 2-bit PLINK genotypes, we perform the
      // following operations:
      //
      // 1. XOR each genotype with 01.  This shuffles the genotype
      // representation to:
      //    00 = missing (important for simplifying step 2)
      //    01 = homozygote 1
      //    10 = homozygote 2
      //    11 = heterozygote
      //
      // 2. Next, compute
      //    mask_i := ((g_i | (g_i >> 1)) & 01) * 11
      // which is 00 whenever g_i is missing, and 11 otherwise.
      //
      // 3. Then, (g_j ^ g_k) & (mask_j & mask_k) distinguishes the
      // possible distances between the genotypes:
      //    - It's equal to zero iff either g_j == g_k or one/both is
      //      missing.  It's fine for these cases to overlap since either
      //      way we do not want to increment the numerator of our final
      //      distance.  (We can handle the effect of missingness on the
      //      denominator outside the main loop.)
      //    - It's equal to 01 or 10 iff neither is missing and exactly
      //      one is a heterozygote.
      //    - It's equal to 11 iff one is homozygote rare and the other is
      //      homozygote common.
      //
      // 4. Finally, we perform the update
      //    A_{jk} := A_{jk} + f_0(x_0) + f_1(x_1) + ... + f_31(x_31)
      // in the nonzero exponent case, or
      //    A_{jk} := A_{jk} + f(x_0) + f(x_1) + ... + f(x_959)
      // in the zero exponent case.
      //
      // For nonzero exponents, we structure the update as
      //    A_{jk} := A_{jk} + f_{0-6}(x_{0-6}) + f_{7-13}(x_{7-13}) +
      //              f_{14-19}(x_{14-19}) + f_{20-25}(x_{20-25}) +
      //              f_{26-31}(x_{26-31})
      // which requires 352 KB of table space.  (This is a conservative
      // choice; the 2 MB 8-8-8-8 table would work better on some newer
      // systems.)
      //
      // See the comments at the beginning of this file for discussion of
      // the zero exponent special case.

      while ((jj < multiplex) && (marker_idx < marker_ct)) {
	while (is_set(marker_exclude, ii)) {
	  ii++;
	  fseeko(pedfile, (off_t)unfiltered_indiv_ct4, SEEK_CUR);
	}
	if (fread(&(pedbuf[jj * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	  goto wdist_ret_READ_FAIL;
	}
	set_allele_freq_buf[jj] = set_allele_freqs[ii];
	if (wt_needed) {
	  wtbuf[jj++] = g_marker_weights_i[marker_idx++];
	} else {
	  jj++;
	  marker_idx++;
	}
	ii++;
      }
      if (jj < multiplex) {
	memset(&(pedbuf[jj * unfiltered_indiv_ct4]), 0, (multiplex - jj) * unfiltered_indiv_ct4);
	if (exp0) {
	  fill_long_zero((long*)g_geno, g_indiv_ct * (MULTIPLEX_2DIST / BITCT));
	  fill_ulong_zero(g_masks, g_indiv_ct * (MULTIPLEX_2DIST / BITCT));
	} else {
	  fill_long_zero((long*)g_geno, g_indiv_ct);
	  fill_ulong_zero(g_masks, g_indiv_ct);
	}
      }
      if (exp0) {
	for (nn = 0; nn < jj; nn += BITCT) {
	  glptr = &(((unsigned long*)g_geno)[nn / BITCT2]);
	  glptr2 = &(g_masks[nn / BITCT2]);
	  glptr3 = g_mmasks;
	  if (wt_needed) {
	    giptr = g_indiv_missing;
	  }
	  if (unwt_needed_full) {
	    giptr2 = g_indiv_missing_unwt;
	  }
	  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
	    if (!is_set(indiv_exclude, indiv_uidx)) {
	      kk = (indiv_uidx % 4) * 2;
	      ulii = 0;
	      ulkk = 0;
	      gptr = &(pedbuf[indiv_uidx / 4 + nn * unfiltered_indiv_ct4]);
	      for (mm = 0; mm < BITCT2; mm++) {
		uljj = (gptr[mm * unfiltered_indiv_ct4] >> kk) & 3;
		ulii |= uljj << (mm * 2);
		if (uljj == 1) {
		  ulkk |= 1LU << mm;
		  if (wt_needed) {
		    *giptr += wtbuf[mm + nn];
		  }
		  if (unwt_needed_full) {
		    *giptr2 += 1;
		  }
		}
	      }
	      // use xor to convert representation to 0 = missing,
	      // 1 or 2 = homozygote, 3 = heterozygote
	      ulii ^= FIVEMASK;
	      *glptr++ = ulii;
	      ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	      *glptr2++ = ulii * 3;
	      *glptr3 = ulkk;
	      ulii = 0;
	      ulkk = 0;
	      gptr = &(pedbuf[indiv_uidx / 4 + (nn + BITCT2) * unfiltered_indiv_ct4]);
	      for (mm = 0; mm < BITCT2; mm++) {
		uljj = (gptr[mm * unfiltered_indiv_ct4] >> kk) & 3;
		ulii |= uljj << (mm * 2);
		if (uljj == 1) {
		  ulkk |= 1LU << mm;
		  if (wt_needed) {
		    *giptr += wtbuf[mm + nn + BITCT2];
		  }
		  if (unwt_needed_full) {
		    *giptr2 += 1;
		  }
		}
	      }
	      ulii ^= FIVEMASK;
	      *glptr = ulii;
	      ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	      *glptr2 = ulii * 3;
	      *glptr3++ |= ulkk << BITCT2;
	      glptr = &(glptr[(MULTIPLEX_2DIST / BITCT) - 1]);
	      glptr2 = &(glptr2[(MULTIPLEX_2DIST / BITCT) - 1]);
	      if (wt_needed) {
		giptr++;
	      }
	      if (unwt_needed_full) {
		giptr2++;
	      }
	    }
	  }

	  if (wt_needed) {
	    g_weights_i = &(wtbuf[nn]);
	    for (ulii = 1; ulii < g_thread_ct; ulii++) {
	      if (pthread_create(&(threads[ulii - 1]), NULL, &calc_distm_thread, (void*)ulii)) {
		goto wdist_ret_THREAD_CREATE_FAIL;
	      }
	    }
	    incr_wt_dist_missing(g_missing_tot_weights, 0);
	    for (uii = 0; uii < g_thread_ct - 1; uii++) {
	      pthread_join(threads[uii], NULL);
	    }
	  }
	  if (unwt_needed) {
	    for (ulii = 1; ulii < g_thread_ct; ulii++) {
	      if (pthread_create(&(threads[ulii - 1]), NULL, &calc_missing_thread, (void*)ulii)) {
		goto wdist_ret_THREAD_CREATE_FAIL;
	      }
	    }
	    incr_dists_rm(g_missing_dbl_excluded, 0, g_thread_start);
	    for (uii = 0; uii < g_thread_ct - 1; uii++) {
	      pthread_join(threads[uii], NULL);
	    }
	  }
	}
	for (ulii = 1; ulii < g_thread_ct; ulii++) {
	  if (pthread_create(&(threads[ulii - 1]), NULL, &calc_idist_thread, (void*)ulii)) {
	    goto wdist_ret_THREAD_CREATE_FAIL;
	  }
	}
	incr_dists_i(g_idists, (unsigned long*)g_geno, 0);
	for (uii = 0; uii < g_thread_ct - 1; uii++) {
	  pthread_join(threads[uii], NULL);
	}
      } else {
	fill_ulong_zero(g_mmasks, g_indiv_ct);
	for (nn = 0; nn < jj; nn += MULTIPLEX_DIST_EXP / 2) {
	  glptr = (unsigned long*)g_geno;
	  glptr2 = g_masks;
	  glptr3 = g_mmasks;
	  giptr3 = g_indiv_missing;
	  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
	    if (!is_set(indiv_exclude, indiv_uidx)) {
	      kk = (indiv_uidx % 4) * 2;
	      ulii = 0;
	      ulkk = 0;
	      gptr = &(pedbuf[indiv_uidx / 4 + nn * unfiltered_indiv_ct4]);
	      for (mm = 0; mm < MULTIPLEX_DIST_EXP / 2; mm++) {
		uljj = (gptr[mm * unfiltered_indiv_ct4] >> kk) & 3;
		ulii |= uljj << (mm * 2);
		if (uljj == 1) {
		  ulkk |= 1LU << mm;
		  *giptr3 += wtbuf[mm + nn];
		}
	      }
#if __LP64__
	      ulii ^= FIVEMASK;
	      *glptr++ = ulii;
	      ulii = (ulii | (ulii >> 1)) & FIVEMASK;
#else
	      // note that FIVEMASK does NOT work here because the mask is
	      // only 28 bits
	      ulii ^= 0x05555555;
	      *glptr++ = ulii;
	      ulii = (ulii | (ulii >> 1)) & 0x05555555;
#endif
	      *glptr2++ = ulii * 3;
	      *glptr3++ |= ulkk << nn;
	      giptr3++;
	    }
	  }
	  fill_weights(g_weights, &(set_allele_freq_buf[nn]), exponent);
	  for (ulii = 1; ulii < g_thread_ct; ulii++) {
	    if (pthread_create(&(threads[ulii - 1]), NULL, &calc_dist_thread, (void*)ulii)) {
	      goto wdist_ret_THREAD_CREATE_FAIL;
	    }
	  }
	  incr_dists(g_dists, (unsigned long*)g_geno, 0);
	  for (uii = 0; uii < g_thread_ct - 1; uii++) {
	    pthread_join(threads[uii], NULL);
	  }
	}
	g_weights_i = wtbuf;
	for (ulii = 1; ulii < g_thread_ct; ulii++) {
	  if (pthread_create(&(threads[ulii - 1]), NULL, &calc_distm_thread, (void*)ulii)) {
	    goto wdist_ret_THREAD_CREATE_FAIL;
	  }
	}
	incr_wt_dist_missing(g_missing_tot_weights, 0);
	for (uii = 0; uii < g_thread_ct - 1; uii++) {
	  pthread_join(threads[uii], NULL);
	}
      }
      printf("\r%d markers complete.", marker_idx);
      fflush(stdout);
    }
    putchar('\r');
    logprint("Distance matrix calculation complete.\n");
    if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
      strcpy(outname_end, ".mdist");
      if (fopen_checked(&outfile, outname, "w")) {
	goto wdist_ret_OPEN_FAIL;
      }
      iptr = g_idists;
      giptr = g_missing_dbl_excluded;
      pct = 1;
      // parallel_tot must be 1 for --distance-matrix
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	giptr2 = g_indiv_missing_unwt;
	uii = marker_ct - giptr2[indiv_idx];
	for (ujj = 0; ujj < indiv_idx; ujj++) {
	  if (fprintf(outfile, "%g ", ((double)(*iptr++)) / (2 * (uii - (*giptr2++) + (*giptr++)))) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fwrite_checked("0 ", 2, outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	giptr2++;
	for (ulii = indiv_idx + 1; ulii < g_indiv_ct; ulii++) {
	  uljj = ulii * (ulii - 1) / 2 + indiv_idx;
	  if (fprintf(outfile, "%g ", ((double)g_idists[uljj]) / (2 * (uii - (*giptr2++) + g_missing_dbl_excluded[uljj]))) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fwrite_checked("\n", 1, outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	if (indiv_idx * 100LL >= ((long long)pct * g_indiv_ct)) {
	  pct = (indiv_idx * 100LL) / g_indiv_ct;
	  printf("\rWriting... %d%%", pct++);
	  fflush(stdout);
	}
      }
      if (fclose_null(&outfile)) {
	goto wdist_ret_WRITE_FAIL;
      }
      putchar('\r');
      sprintf(logbuf, "Distances (proportions) written to %s.\n", outname);
      logprintb();
      if (!parallel_idx) {
	strcpy(outname_end, ".mdist.id");
	retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
	if (retval) {
	  goto wdist_ret_2;
	}
      }
    }
    if (calculation_type & CALC_PLINK_IBS_MATRIX) {
      strcpy(outname_end, ".mibs");
      if (fopen_checked(&outfile, outname, "w")) {
	goto wdist_ret_OPEN_FAIL;
      }
      iptr = g_idists;
      giptr = g_missing_dbl_excluded;
      pct = 1;
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	giptr2 = g_indiv_missing_unwt;
	uii = marker_ct - giptr2[indiv_idx];
	for (ujj = 0; ujj < indiv_idx; ujj++) {
	  if (fprintf(outfile, "%g ", 1.0 - (((double)(*iptr++)) / (2 * (uii - (*giptr2++) + (*giptr++))))) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fwrite_checked("1 ", 2, outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	giptr2++;
	for (ulii = indiv_idx + 1; ulii < g_indiv_ct; ulii++) {
	  uljj = (ulii * (ulii - 1)) / 2 + indiv_idx;
	  if (fprintf(outfile, "%g ", 1.0 - (((double)g_idists[uljj]) / (2 * (uii - (*giptr2++) + g_missing_dbl_excluded[uljj])))) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fwrite_checked("\n", 1, outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	if (indiv_idx * 100 >= (pct * g_indiv_ct)) {
	  pct = (indiv_idx * 100) / g_indiv_ct;
	  printf("\rWriting... %d%%", pct++);
	  fflush(stdout);
	}
      }
      fclose_null(&outfile);
      putchar('\r');
      sprintf(logbuf, "IBS matrix written to %s.\n", outname);
      logprintb();
      strcpy(outname_end, ".mibs.id");
      retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    tstc = g_thread_start[g_thread_ct];
    if (wt_needed) {
      giptr = g_missing_tot_weights;
      dptr2 = g_dists;
      if (exp0) {
	iptr = g_idists;
	for (indiv_idx = g_thread_start[0]; indiv_idx < tstc; indiv_idx++) {
	  giptr2 = g_indiv_missing;
	  uii = giptr2[indiv_idx];
	  for (ujj = 0; ujj < indiv_idx; ujj++) {
	    *dptr2++ = (4294967295.0 / ((4294967295U - uii - (*giptr2++)) + (*giptr++))) * (*iptr++);
	  }
	}
      } else {
	for (indiv_idx = g_thread_start[0]; indiv_idx < tstc; indiv_idx++) {
	  giptr2 = g_indiv_missing;
	  uii = giptr2[indiv_idx];
	  for (ujj = 0; ujj < indiv_idx; ujj++) {
	    *dptr2 *= (4294967295.0 / ((4294967295U - uii - (*giptr2++)) + (*giptr++)));
	    dptr2++;
	  }
	}
      }
    } else if (distance_flat_missing) {
      dptr2 = g_dists;
      giptr = g_missing_dbl_excluded;
      if (exp0) {
        iptr = g_idists;
	for (indiv_idx = g_thread_start[0]; indiv_idx < tstc; indiv_idx++) {
	  giptr2 = g_indiv_missing_unwt;
	  uii = marker_ct - giptr2[indiv_idx];
	  for (ujj = 0; ujj < indiv_idx; ujj++) {
	    *dptr2++ = (((double)marker_ct) / (uii - (*giptr2++) + (*giptr++))) * (*iptr++);
	  }
	}

      } else {
	for (indiv_idx = g_thread_start[0]; indiv_idx < tstc; indiv_idx++) {
	  giptr2 = g_indiv_missing_unwt;
	  uii = marker_ct - giptr2[indiv_idx];
	  for (ujj = 0; ujj < indiv_idx; ujj++) {
	    *dptr2 *= ((double)marker_ct) / (uii - (*giptr2++) + (*giptr++));
	    dptr2++;
	  }
	}
      }
    }
  }

  if (calculation_type & CALC_DISTANCE) {
    if (!parallel_idx) {
      retval = distance_d_write_ids(outname, outname_end, dist_calc_type, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    if ((exponent == 0.0) || (!(dist_calc_type & (DISTANCE_IBS | DISTANCE_1_MINUS_IBS)))) {
      dxx = 0.5 / (double)marker_ct;
    } else {
      dxx = 0.0;
      marker_uidx = 0;
      for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	dyy = set_allele_freqs[marker_uidx];
	if ((dyy > 0.0) && (dyy < 1.0)) {
	  dxx += pow(2 * dyy * (1.0 - dyy), -exponent);
	} else {
	  dxx += 1.0;
	}
	marker_uidx++;
      }
      dxx = 0.5 / dxx;
    }
    retval = distance_d_write(&outfile, &outfile2, &outfile3, &gz_outfile, &gz_outfile2, &gz_outfile3, dist_calc_type, outname, outname_end, g_dists, dxx, g_indiv_ct, g_thread_start[0], g_thread_start[g_thread_ct], parallel_idx, parallel_tot, g_geno);
    if (retval) {
      goto wdist_ret_2;
    }
    wkspace_reset(wkspace_mark);
  } else if (calculation_type & CALC_LOAD_DISTANCES) {
    dists_alloc = ((long)g_indiv_ct * (g_indiv_ct - 1)) * (sizeof(double) / 2);
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
    if (ftello(loaddistfile) != (long long)dists_alloc) {
      sprintf(logbuf, "Invalid --load-dists file.  (Triangular binary of size %lld expected.)\n", dists_alloc);
      logprintb();
      goto wdist_ret_INVALID_FORMAT;
    }
    rewind(loaddistfile);
    if (fread(g_dists, 1, dists_alloc, loaddistfile) < dists_alloc) {
      goto wdist_ret_READ_FAIL;
    }
    fclose(loaddistfile);
  }

  if (calculation_type & CALC_GROUPDIST) {
    wkspace_mark = wkspace_base;
    collapse_arr(g_pheno_c, 1, indiv_exclude, unfiltered_indiv_ct);
    g_low_ct = 0;
    g_high_ct = 0;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      if (g_pheno_c[indiv_idx] == 1) {
	g_high_ct++;
      } else if (g_pheno_c[indiv_idx] == 0) {
	g_low_ct++;
      }
    }
    ll_size = ((unsigned long)g_low_ct * (g_low_ct - 1)) / 2;
    lh_size = g_low_ct * g_high_ct;
    hh_size = ((unsigned long)g_high_ct * (g_high_ct - 1)) / 2;
    g_reg_tot_y = 0.0;
    g_reg_tot_xy = 0.0;
    g_reg_tot_x = 0.0;
    dhh_ssq = 0.0;
    dhl_ssq = 0.0;
    dll_ssq = 0.0;
    dist_ptr = g_dists;
    ll_pool = (double*)wkspace_alloc(ll_size * sizeof(double));
    lh_pool = (double*)wkspace_alloc(lh_size * sizeof(double));
    hh_pool = (double*)wkspace_alloc(hh_size * sizeof(double));
    ll_poolp = ll_pool;
    lh_poolp = lh_pool;
    hh_poolp = hh_pool;
    g_geno = wkspace_alloc(g_thread_ct * CACHEALIGN(g_high_ct + g_low_ct + (g_jackknife_d + 1) * sizeof(int)));
    for (indiv_idx = 1; indiv_idx < g_indiv_ct; indiv_idx++) {
      cptr = g_pheno_c;
      cptr2 = &(g_pheno_c[indiv_idx]);
      if (*cptr2 == 1) {
	while (cptr < cptr2) {
	  if (*cptr == 1) {
            dxx = *dist_ptr;
	    *hh_poolp++ = dxx;
	    g_reg_tot_x += dxx;
            dhh_ssq += dxx * dxx;
	  } else if (!(*cptr)) {
            dxx = *dist_ptr;
	    *lh_poolp++ = dxx;
	    g_reg_tot_xy += dxx;
            dhl_ssq += dxx * dxx;
	  }
	  cptr++;
	  dist_ptr++;
	}
      } else if (*cptr2 == 0) {
	while (cptr < cptr2) {
	  if (*cptr == 1) {
            dxx = *dist_ptr;
	    *lh_poolp++ = dxx;
	    g_reg_tot_xy += dxx;
            dhl_ssq += dxx * dxx;
	  } else if (!(*cptr)) {
            dxx = *dist_ptr;
	    *ll_poolp++ = dxx;
	    g_reg_tot_y += dxx;
            dll_ssq += dxx * dxx;
	  }
	  cptr++;
	  dist_ptr++;
	}
      } else {
	dist_ptr += indiv_idx;
      }
    }
#ifdef __cplusplus
    // std::sort is faster than qsort for basic types.  See e.g. Anders
    // Kaseorg's answer to
    // http://www.quora.com/Software-Engineering/Generally-how-much-faster-is-C-compared-to-C++
    std::sort(ll_pool, &(ll_pool[ll_size]));
    std::sort(lh_pool, &(lh_pool[lh_size]));
    std::sort(hh_pool, &(hh_pool[hh_size]));
#else
    qsort(ll_pool, ll_size, sizeof(double), double_cmp);
    qsort(lh_pool, lh_size, sizeof(double), double_cmp);
    qsort(hh_pool, hh_size, sizeof(double), double_cmp);
#endif
    ll_med = get_dmedian(ll_pool, ll_size);
    lh_med = get_dmedian(lh_pool, lh_size);
    hh_med = get_dmedian(hh_pool, hh_size);
    sprintf(logbuf, "Case/control distance analysis (%d affected, %d unaffected):\n", g_high_ct, g_low_ct);
    logprintb();
    if (g_high_ct < 2) {
      dxx = 0.0;
      dhh_sd = 0.0;
    } else {
      dww = (double)(((unsigned long)g_high_ct * (g_high_ct - 1)) / 2);
      dxx = g_reg_tot_x / dww;
      dhh_sd = sqrt((dhh_ssq / dww - dxx * dxx) / (dww - 1.0));
    }
    if (!(g_high_ct * g_low_ct)) {
      dyy = 0.0;
      dhl_sd = 0.0;
    } else {
      dww = (double)((unsigned long)g_high_ct * g_low_ct);
      dyy = g_reg_tot_xy / dww;
      dhl_sd = sqrt((dhl_ssq / dww - dyy * dyy) / (dww - 1.0));
    }
    if (g_low_ct < 2) {
      dzz = 0.0;
      dll_sd = 0.0;
    } else {
      dww = (double)(((unsigned long)g_low_ct * (g_low_ct - 1)) / 2);
      dzz = g_reg_tot_y / dww;
      dll_sd = sqrt((dll_ssq / dww - dzz * dzz) / (dww - 1.0));
    }
    sprintf(logbuf, "  Mean (sd), median dists between 2x affected     : %g (%g), %g\n", dxx, dhh_sd, hh_med);
    logprintb();
    sprintf(logbuf, "  Mean (sd), median dists between aff. and unaff. : %g (%g), %g\n", dyy, dhl_sd, lh_med);
    logprintb();
    sprintf(logbuf, "  Mean (sd), median dists between 2x unaffected   : %g (%g), %g\n\n", dzz, dll_sd, ll_med);
    logprintb();
    if (2 * g_jackknife_d >= (g_high_ct + g_low_ct)) {
      logprint("Delete-d jackknife skipped because d is too large.\n");
    } else {
      // this can be sped up using the same method used in regress-distance,
      // if it's important
      g_jackknife_iters = (groupdist_iters + g_thread_ct - 1) / g_thread_ct;

      if (groupdist_d) {
	g_jackknife_d = groupdist_d;
      } else {
	g_jackknife_d = set_default_jackknife_d(g_high_ct + g_low_ct);
      }

      for (ulii = 1; ulii < g_thread_ct; ulii++) {
	if (pthread_create(&(threads[ulii - 1]), NULL, &groupdist_jack_thread, (void*)ulii)) {
	  goto wdist_ret_THREAD_CREATE_FAIL;
	}
      }
      ulii = 0;
      groupdist_jack_thread((void*)ulii);
      for (uii = 1; uii < g_thread_ct; uii++) {
        pthread_join(threads[uii - 1], NULL);
        for (jj = 0; jj < 9; jj++) {
	  g_calc_result[jj][0] += g_calc_result[jj][uii];
        }
      }
      dxx = 1.0 / g_thread_ct;
      g_calc_result[0][0] *= dxx;
      g_calc_result[1][0] *= dxx;
      g_calc_result[2][0] *= dxx;
      dxx /= (g_jackknife_iters - 1) * g_thread_ct;
      for (ii = 3; ii < 9; ii++) {
        g_calc_result[ii][0] *= dxx;
      }
      putchar('\r');
      sprintf(logbuf, "  AA mean - AU mean avg difference (s.e.): %g (%g)\n", g_calc_result[0][0] - g_calc_result[1][0], sqrt(((g_high_ct + g_low_ct) / ((double)g_jackknife_d)) * (g_calc_result[3][0] + g_calc_result[4][0] - 2 * g_calc_result[6][0])));
      logprintb();
      sprintf(logbuf, "  AA mean - UU mean avg difference (s.e.): %g (%g)\n", g_calc_result[0][0] - g_calc_result[2][0], sqrt(((g_high_ct + g_low_ct) / ((double)g_jackknife_d)) * (g_calc_result[3][0] + g_calc_result[5][0] - 2 * g_calc_result[7][0])));
      logprintb();
      sprintf(logbuf, "  AU mean - UU mean avg difference (s.e.): %g (%g)\n", g_calc_result[1][0] - g_calc_result[2][0], sqrt(((g_high_ct + g_low_ct) / ((double)g_jackknife_d)) * (g_calc_result[4][0] + g_calc_result[5][0] - 2 * g_calc_result[8][0])));
      logprintb();
    }
    wkspace_reset(wkspace_mark);
  }
  if (calculation_type & CALC_REGRESS_DISTANCE) {
    retval = regress_distance(calculation_type, g_dists, g_pheno_d, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, g_thread_ct, regress_iters, regress_d);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  if (calculation_type & CALC_GENOME) {
    wkspace_reset(wkspace_mark2);
    retval = calc_genome(threads, pedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, marker_pos, set_allele_freqs, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info, parallel_idx, parallel_tot, outname, outname_end, nonfounders, calculation_type, genome_output_gz, genome_output_full, genome_ibd_unbounded, ppc_gap, pri);
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
  wdist_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  wdist_ret_INVALID_CMDLINE_2:
    logprintb();
  wdist_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  wdist_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  wdist_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  wdist_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    for (uljj = 0; uljj < ulii - 1; uljj++) {
      pthread_join(threads[uljj], NULL);
    }
    break;
  }
 wdist_ret_2:
  gzclose_cond(gz_outfile3);
  gzclose_cond(gz_outfile2);
  gzclose_cond(gz_outfile);
  fclose_cond(bedtmpfile);
  fclose_cond(bimtmpfile);
  fclose_cond(famtmpfile);
 wdist_ret_1:
  free_cond(line_locs);
  free_cond(line_mids);
  free_cond(het_probs);
  free_cond(g_pheno_d);
  free_cond(g_pheno_c);
  free_cond(phenor_d);
  free_cond(phenor_c);
  free_cond(id_buf);
  free_cond(id_list);
  free_cond(pid_list);
  free_cond(pedbuf);
  fclose_cond(loaddistfile);
  fclose_cond(freqfile);
  fclose_cond(filterfile);
  fclose_cond(phenofile);
  fclose_cond(famfile);
  fclose_cond(mapfile);
  fclose_cond(pedfile);
  fclose_cond(outfile3);
  fclose_cond(outfile2);
  fclose_cond(outfile);
  return retval;
}

// output-missing-phenotype + terminating null
#define MAX_FLAG_LEN 25

static inline int is_flag(char* param) {
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

int param_count(int argc, char** argv, int flag_idx) {
  // Counts the number of optional parameters given to the flag at position
  // flag_idx, treating any parameter not beginning with "--" as optional.
  int opt_params = 0;
  int cur_idx = flag_idx + 1;
  while (cur_idx < argc) {
    if (is_flag(argv[cur_idx])) {
      break;
    }
    opt_params++;
    cur_idx++;
  }
  return opt_params;
}

int enforce_param_ct_range(int argc, char** argv, int flag_idx, int min_ct, int max_ct, int* param_ct_ptr) {
  int param_ct = param_count(argc, argv, flag_idx);
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

int species_flag(Chrom_info* chrom_info_ptr, int species_val) {
  if (chrom_info_ptr->species) {
    logprint("Error: Duplicate species flag.\n");
    return -1;
  }
  chrom_info_ptr->species = species_val;
  if (chrom_info_ptr->chrom_mask & (~(species_valid_chrom_mask[species_val]))) {
    logprint("Error: Invalid --chr parameter.\n");
    return -1;
  }
  return 0;
}

void invalid_arg(char* argv) {
  sprintf(logbuf, "Error: Invalid flag ('%s').%s%s", argv, (argv[0] == '-')? "" : "  All flags must be preceded by 1-2 dashes.", errstr_append);
}

void print_ver() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}

int alloc_string(char** sbuf, char* source) {
  unsigned int slen = strlen(source) + 1;
  *sbuf = (char*)malloc(slen * sizeof(char));
  if (!(*sbuf)) {
    return -1;
  }
  memcpy(*sbuf, source, slen);
  return 0;
}

int alloc_fname(char** fnbuf, char* source, char* argptr, unsigned int extra_size) {
  unsigned int slen = strlen(source) + 1;
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

int main(int argc, char** argv) {
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
  char** subst_argv = NULL;
  char* script_buf = NULL;
  char* rerun_buf = NULL;
  char* flag_buf = NULL;
  unsigned int* flag_map = NULL;
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
  char** subst_argv2;
  int retval = 0;
  int load_params = 0; // describes what file parameters have been provided
  int load_rare = 0;
  int fam_col_1 = 1;
  int fam_col_34 = 1;
  int fam_col_5 = 1;
  int fam_col_6 = 1;
  int mpheno_col = 0;
  int affection_01 = 0;
  double exponent = 0.0;
  double min_maf = 0.0;
  double max_maf = 0.5;
  double geno_thresh = 1.0;
  double mind_thresh = 1.0;
  double hwe_thresh = 0.0;
  int hwe_all = 0;
  double rel_cutoff = 0.025;
  int cur_arg = 1;
  int calculation_type = 0;
  int rel_calc_type = 0;
  int dist_calc_type = 0;
  char* bubble;
  int mfilter_col = 0;
  int pheno_merge = 0;
  int tail_pheno = 0;
  int prune = 0;
  int missing_pheno = -9;
  char missing_geno = '0';
  char output_missing_pheno[32];
  char output_missing_geno = '0';
  double tail_bottom;
  double tail_top;
  int distance_3d = 0;
  int distance_flat_missing = 0;
  unsigned long groupdist_iters = ITERS_DEFAULT;
  int groupdist_d = 0;
  unsigned long regress_iters = ITERS_DEFAULT;
  int regress_d = 0;
  unsigned long regress_rel_iters = ITERS_DEFAULT;
  int regress_rel_d = 0;
  double unrelated_herit_tol = 0.0000001;
  double unrelated_herit_covg = 0.45;
  double unrelated_herit_covr = 0.55;
  int ibc_type = 0;
  int parallel_idx = 0;
  int parallel_tot = 1;
  int allow_no_sex = 0;
  int nonfounders = 0;
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int oo;
  int ppc_gap = DEFAULT_PPC_GAP;
  unsigned long int rseed = 0;
  int genome_output_gz = 0;
  int genome_output_full = 0;
  int genome_ibd_unbounded = 0;
  FILE* scriptfile = NULL;
  int num_params;
  int in_param;
  int ld_window_size = 0;
  int ld_window_incr = 0;
  double ld_last_param = 0.0;
  int ld_window_kb = 0;
  int maf_succ = 0;
  int autosome = 0;
  int filter_case_control = 0;
  int filter_sex = 0;
  int filter_founder_nonf = 0;
  int regress_pcs_normalize_pheno = 0;
  int regress_pcs_sex_specific = 0;
  int regress_pcs_clip = 0;
  int max_pcs = MAX_PCS_DEFAULT;
  int freqx = 0;
  int recode_modifier = 0;
  int freq_counts = 0;
  int allelexxxx = 0;
  int silent = 0;
  int merge_type = 0;
  int indiv_sort = 0;
  int keep_allele_order = 0;
  unsigned int cur_flag = 0;
  unsigned int flag_ct = 0;
  int dummy_marker_ct = 0;
  int dummy_indiv_ct = 0;
  unsigned int dummy_flags = 0;
  double dummy_missing_geno = 0.0;
  double dummy_missing_pheno = 0.0;
  Chrom_info chrom_info;
  char* argptr2;
  char* flagptr;
  char* missing_code = NULL;
  double dxx;
  char cc;
  unsigned int uii;
  long default_alloc_mb;
  long long llxx;
#ifdef __APPLE__
  int mib[2];
  size_t sztmp;
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
	if (no_more_items(sptr)) {
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
  flag_map = (unsigned int*)malloc(flag_ct * sizeof(int));
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
      case 'm':
	if (!memcmp(argptr, "missing_code", 13)) {
	  memcpy(flagptr, "missing-code", 13);
	  break;
	}
	goto main_flag_copy;
      case 'r':
	if (!memcmp(argptr, "recode", 6)) {
	  if ((!memcmp(&(argptr[6]), "12", 3)) || (!memcmp(&(argptr[6]), "-lgen", 6))) {
	    if (kk == 8) {
	      memcpy(flagptr, "recode 12", 10);
	    } else {
	      memcpy(flagptr, "recode lgen", 12);
	    }
	    printf("Note: %s flag deprecated.  Use '%s ...'.\n", argptr, flagptr);
	    mm++;
            break;
	  }
	}
	goto main_flag_copy;
      case 'u':
	if (!memcmp(argptr, "update-freq", 12)) {
	  memcpy(flagptr, "read-freq", 10);
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
  qsort_ext2(flag_buf, flag_ct, MAX_FLAG_LEN, strcmp_deref, (char*)flag_map, sizeof(int), sptr, MAX_FLAG_LEN);
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
      silent = 1;
      break;
    }
  }
  print_ver();
  memcpy(outname, "wdist", 6);
  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    if (!memcmp("out", &(flag_buf[cur_flag * MAX_FLAG_LEN]), 4)) {
      ii = flag_map[cur_flag];
      if (enforce_param_ct_range(argc, argv, ii, 1, 1, &jj)) {
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      if (strlen(argv[ii + 1]) > (FNAMESIZE - MAX_POST_EXT)) {
	fputs("Error: --out parameter too long.\n", stdout);
	goto main_ret_OPEN_FAIL;
      }
      strcpy(outname, argv[ii + 1]);
    }
  }
  uii = strlen(outname);
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
  if (gethostname(tbuf, 4 * MAXLINELEN + 256) != -1) {
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
  ii = sysconf(_SC_NPROCESSORS_ONLN);
  if (ii == -1) {
    g_thread_ct = 1;
  } else if (ii > MAX_THREADS) {
    g_thread_ct = MAX_THREADS;
  } else {
    g_thread_ct = ii;
  }
  strcpy(mapname, "wdist.map");
  strcpy(pedname, "wdist.ped");
  famname[0] = '\0';
  genname[0] = '\0';
  samplename[0] = '\0';
  memcpy(output_missing_pheno, "-9", 3);
  cur_flag = 0;
  do {
    argptr = &(flag_buf[cur_flag * MAX_FLAG_LEN]);
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
	autosome = 1;
      } else if (!memcmp(argptr2, "llow-no-sex", 12)) {
	allow_no_sex = 1;
      } else if (!memcmp(argptr2, "ll", 3)) {
	logprint("Note: --all flag has no effect.\n");
      } else if (!memcmp(argptr2, "llele1234", 10)) {
	logprint("Note: --allele1234 not yet written.\n");
	allelexxxx = 1;
      } else if (!memcmp(argptr2, "lleleACGT", 9)) {
	if (allelexxxx) {
	  if (allelexxxx == 1) {
	    logprint("Error: --allele1234 and --alleleACGT flags cannot coexist.\n");
	  }
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --alleleACGT not yet written.\n");
	allelexxxx = 2;
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
      if (!memcmp(argptr2, "ow", 3)) {
	if (species_flag(&chrom_info, SPECIES_COW)) {
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "hr", 3)) {
	if (autosome) {
	  logprint("Error: --chr and --autosome flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	ii = param_count(argc, argv, cur_arg);
	if (!ii) {
	  sprintf(logbuf, "Error: Missing --chr parameter.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (jj = 0; jj < ii; jj++) {
	  kk = marker_code_raw(argv[cur_arg + 1 + jj]);
	  // will allow some non-autosomal in future
	  if ((kk == -1) || (kk == CHROM_X) || (kk == CHROM_Y) || (kk == CHROM_MT) || ((kk == CHROM_XY) && (species_xy_code[chrom_info.species] == -1))) {
	    sprintf(logbuf, "Error: Invalid --chr parameter '%s'.\n", argv[cur_arg + 1 + jj]);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  chrom_info.chrom_mask |= 1LLU << kk;
	}
      } else if (!memcmp(argptr2, "ompound-genotypes", 18)) {
	logprint("Note: --compound-genotypes flag unnecessary (spaces between alleles in .ped\nare optional).\n");
      } else if (!memcmp(argptr2, "ompress", 8)) {
	logprint("Error: --compress flag retired.  Use 'gzip [filename]'.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ounts", 6)) {
	logprint("Note: --counts flag deprecated.  Use '--freq counts' or --freqx instead.\n");
        freq_counts = 1;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'd':
      if (!memcmp(argptr2, "ata", 4)) {
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
	logprint("Error: --decompress flag retired.  Use 'gunzip [filename]'.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "og", 3)) {
	if (species_flag(&chrom_info, SPECIES_DOG)) {
	  goto main_ret_INVALID_CMDLINE;
	}
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
	if (!(dist_calc_type & DISTANCE_FORMATMASK)) {
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
	dummy_marker_ct = atoi(argv[cur_arg + 1]);
	if (dummy_marker_ct < 1) {
	  sprintf(logbuf, "Error: Invalid --dummy marker count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dummy_indiv_ct = atoi(argv[cur_arg + 2]);
	if (dummy_indiv_ct < 1) {
	  sprintf(logbuf, "Error: Invalid --dummy individual count.%s", errstr_append);
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
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'h':
      if (!memcmp(argptr2, "orse", 5)) {
	if (species_flag(&chrom_info, SPECIES_HORSE)) {
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "we", 3)) {
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
      } else if (!memcmp(argptr2, "ouse", 5)) {
	if (species_flag(&chrom_info, SPECIES_MOUSE)) {
	  goto main_ret_INVALID_CMDLINE;
	}
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
	  if (!strcmp(argv[cur_arg + jj], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-grm 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-grm 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if (!strcmp(argv[cur_arg + jj], "no-gz")) {
	    rel_calc_type &= ~REL_CALC_GZ;
	  } else if ((!strcmp(argv[cur_arg + jj], "ibc1")) || (!strcmp(argv[cur_arg + jj], "ibc2"))) {
	    if (rel_calc_type & REL_CALC_COV) {
	      sprintf(logbuf, "Error: --make-grm 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-grm '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + jj], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ibc_type = argv[cur_arg + jj][3] - '0';
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
	  if (!strcmp(argv[cur_arg + jj], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if (!strcmp(argv[cur_arg + jj], "gz")) {
	    if (rel_calc_type & REL_CALC_BIN) {
	      sprintf(logbuf, "Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_GZ;
	  } else if (!strcmp(argv[cur_arg + jj], "bin")) {
	    if (rel_calc_type & REL_CALC_GZ) {
	      sprintf(logbuf, "Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_BIN;
	  } else if (!strcmp(argv[cur_arg + jj], "square")) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_SQ;
	  } else if (!strcmp(argv[cur_arg + jj], "square0")) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      sprintf(logbuf, "Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_SQ0;
	  } else if (!strcmp(argv[cur_arg + jj], "triangle")) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      sprintf(logbuf, "Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_TRI;
	  } else if ((!strcmp(argv[cur_arg + jj], "ibc1")) || (!strcmp(argv[cur_arg + jj], "ibc2"))) {
	    if (rel_calc_type & REL_CALC_COV) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-rel '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + jj], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ibc_type = argv[cur_arg + jj][3] - '0';
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-rel parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (!(rel_calc_type & REL_CALC_SHAPEMASK)) {
	  rel_calc_type |= REL_CALC_TRI;
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
	pheno_merge = 1;
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
	if ((argv[cur_arg + 1][0] < '0') || (argv[cur_arg + 1][0] > '9')) {
	  sprintf(logbuf, "Error: Invalid --ppc-gap parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ppc_gap = atoi(argv[cur_arg + 1]);
	if (ppc_gap > 2147483) {
	  sprintf(logbuf, "Error: --ppc-gap parameter '%s' too large.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ppc_gap *= 1000;
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
      } else if (!memcmp(argptr2, "ice", 4)) {
	logprint("Error: --rice not yet supported.\n");
	goto main_ret_INVALID_CMDLINE;
	// if (species_flag(&chrom_info, SPECIES_RICE)) {
	//   goto main_ret_INVALID_CMDLINE;
	// }
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
	  sprintf(logbuf, "Error: --parallel and --regress-rel flags cannot coexist.%s", errstr_append);
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
      } else if ((!memcmp(argptr2, "ecode", 6)) || (!memcmp(argptr2, "ecode 12", 9)) || (!memcmp(argptr2, "ecode lgen", 11))) {
	if (argptr2[5] == ' ') {
	  if (argptr2[6] == '1') {
	    recode_modifier |= RECODE_12;
	  } else {
	    recode_modifier |= RECODE_LGEN;
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
	    recode_modifier |= RECODE_12;
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
	  } else if (!memcmp(argv[cur_arg + jj], "transpose", 10)) {
	    if (recode_modifier & RECODE_LGEN) {
	      sprintf(logbuf, "Error: --recode 'transpose' and 'lgen' modifiers cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_TRANSPOSE;
	  } else if (!memcmp(argv[cur_arg + jj], "lgen", 5)) {
	    if (recode_modifier & RECODE_TRANSPOSE) {
	      sprintf(logbuf, "Error: --recode 'transpose' and 'lgen' modifiers cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_LGEN;
	  } else {
	    sprintf(logbuf, "Error: Invalid --recode parameter '%s'.%s", argv[cur_arg + jj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_RECODE;
      } else if (!memcmp(argptr2, "eference-allele", 16)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&refalleles, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
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
      } else if (!memcmp(argptr2, "heep", 5)) {
	if (species_flag(&chrom_info, SPECIES_SHEEP)) {
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "eed", 4)) {
	if (enforce_param_ct_range(argc, argv, cur_arg, 1, 1, &ii)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	rseed = (unsigned long int)atol(argv[cur_arg + 1]);
	if (rseed == 0) {
	  sprintf(logbuf, "Error: Invalid --seed parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
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
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'w':
      if (!memcmp(argptr2, "rite-snplist", 13)) {
	calculation_type |= CALC_WRITE_SNPLIST;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    default:
      goto main_ret_INVALID_CMDLINE_2;
    }

  } while ((++cur_flag) < flag_ct);
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
  }
  if ((calculation_type & CALC_LOAD_DISTANCES) && (!(calculation_type & (CALC_GROUPDIST | CALC_REGRESS_DISTANCE)))) {
    sprintf(logbuf, "Error: --load-dists cannot be used without either --groupdist or\n--regress-distance.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
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
  sztmp = sizeof(long long);
  sysctl(mib, 2, &llxx, &sztmp, NULL, 0);
  llxx /= 1048576;
#else
  llxx = ((size_t)sysconf(_SC_PHYS_PAGES)) * ((size_t)sysconf(_SC_PAGESIZE)) / 1048576;
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
    sprintf(logbuf, "%lld MB RAM detected; reserving %ld MB for main workspace.\n", llxx, malloc_size_mb);
  } else {
    sprintf(logbuf, "Failed to calculate system memory.  Attempting to reserve %ld MB.\n", malloc_size_mb);
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
      sprintf(logbuf, "Allocated %ld MB successfully, after larger attempt(s) failed.\n", malloc_size_mb);
      logprintb();
    } else if (malloc_size_mb == WKSPACE_MIN_MB) {
      goto main_ret_NOMEM;
    }
  }
  // force 64-byte align on OS X to make cache line sensitivity work
  wkspace = (unsigned char*)CACHEALIGN((unsigned long)wkspace_ua);
  wkspace_base = wkspace;
  wkspace_left = malloc_size_mb * 1048576 - (unsigned long)(wkspace - wkspace_ua);
  free(bubble);

  if (!rseed) {
    rseed = (unsigned long int)time(NULL);
  }
  init_genrand(rseed);

  if (autosome) {
    chrom_info.chrom_mask = species_autosome_mask[chrom_info.species];
  } else if (!chrom_info.chrom_mask) {
    chrom_info.chrom_mask = species_def_chrom_mask[chrom_info.species];
  } else if (chrom_info.chrom_mask & (1LLU << CHROM_XY)) {
    chrom_info.chrom_mask = ((chrom_info.chrom_mask & ((1LLU << MAX_POSSIBLE_CHROM) - 1)) | (1LLU << species_xy_code[chrom_info.species]));
  }

  tbuf[MAXLINELEN - 6] = ' ';
  tbuf[MAXLINELEN - 1] = ' ';
  if (load_rare == LOAD_RARE_GRM) {
    // --rel-cutoff batch mode special case
    retval = rel_cutoff_batch(pedname, outname, rel_cutoff, rel_calc_type);
  } else if (genname[0]) {
    if (calculation_type & (~(CALC_DISTANCE | CALC_REGRESS_DISTANCE))) {
      logprint("Error: Only --distance calculations are currently supported with --data.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
    } else {
      if (!missing_code) {
	missing_code = (char*)"NA";
      }
      retval = wdist_dosage(calculation_type, dist_calc_type, genname, samplename, outname, missing_code, distance_3d, distance_flat_missing, exponent, maf_succ, regress_iters, regress_d, g_thread_ct, parallel_idx, parallel_tot);
    }
  } else {
  // famname[0] indicates binary vs. text
  // extractname, excludename, keepname, and removename indicate the presence
  // of their respective flags
  // filtername indicates existence of filter
  // freqname signals --read-freq
    if (load_rare) {
    // if (load_rare || (!famname[0])) {
      ii = 0;
      sptr = (char*)memchr(outname, 0, FNAMESIZE);
      if (bed_suffix_conflict(calculation_type, recode_modifier)) {
	ii = 8;
        memcpy(sptr, "-working", 9);
      }
      uii = (sptr - outname) + ii;
      if (load_rare == LOAD_RARE_LGEN) {
        retval = lgen_to_bed(pedname, outname, missing_pheno, affection_01, &chrom_info);
      } else if (load_rare & LOAD_RARE_TRANSPOSE_MASK) {
        retval = transposed_to_bed(pedname, famname, outname, missing_geno, &chrom_info);
      } else if (load_rare & LOAD_RARE_DUMMY) {
	retval = generate_dummy(outname, dummy_flags, (unsigned int)dummy_marker_ct, (unsigned int)dummy_indiv_ct, dummy_missing_geno, dummy_missing_pheno);
      } else {
        retval = ped_to_bed(pedname, mapname, outname, fam_col_1, fam_col_34, fam_col_5, fam_col_6, affection_01, missing_pheno, &chrom_info);
	fam_col_1 = 1;
	fam_col_34 = 1;
	fam_col_5 = 1;
        fam_col_6 = 1;
	affection_01 = 0;
	missing_pheno = -9;
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
      outname[uii - ii] = '\0';
    }
    retval = wdist(outname, pedname, mapname, famname, phenoname, extractname, excludename, keepname, removename, filtername, freqname, loaddistname, evecname, mergename1, mergename2, mergename3, makepheno_str, phenoname_str, refalleles, filterval, mfilter_col, filter_case_control, filter_sex, filter_founder_nonf, fam_col_1, fam_col_34, fam_col_5, fam_col_6, missing_geno, missing_pheno, output_missing_geno, output_missing_pheno, mpheno_col, pheno_merge, prune, affection_01, &chrom_info, exponent, min_maf, max_maf, geno_thresh, mind_thresh, hwe_thresh, hwe_all, rel_cutoff, tail_pheno, tail_bottom, tail_top, calculation_type, rel_calc_type, dist_calc_type, groupdist_iters, groupdist_d, regress_iters, regress_d, regress_rel_iters, regress_rel_d, unrelated_herit_tol, unrelated_herit_covg, unrelated_herit_covr, ibc_type, parallel_idx, (unsigned int)parallel_tot, ppc_gap, allow_no_sex, nonfounders, genome_output_gz, genome_output_full, genome_ibd_unbounded, ld_window_size, ld_window_kb, ld_window_incr, ld_last_param, maf_succ, regress_pcs_normalize_pheno, regress_pcs_sex_specific, regress_pcs_clip, max_pcs, freq_counts, freqx, distance_flat_missing, recode_modifier, allelexxxx, merge_type, indiv_sort, keep_allele_order);
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
