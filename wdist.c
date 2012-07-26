// WDIST weighted genetic distance calculator
// Copyright (C) 2012  Christopher Chang  chrchang523@gmail.com

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
// #include "gsl/rng/gsl_rng.h"

#define RET_SUCCESS 0
#define RET_HELP 1
#define RET_NOMEM 2
#define RET_OPENFAIL 3
#define RET_INVALID_FORMAT 4
#define RET_CALC_NOT_YET_SUPPORTED 5
#define RET_INVALID_CMDLINE 6
#define RET_WRITE_FAIL 7
#define RET_READ_FAIL 8
#define RET_MOREHELP 9
#define RET_THREAD_CREATE_FAIL 10

#define FILTER_KEEP 1
#define FILTER_REMOVE 2
#define FILTER_CUSTOM 3

#define CALC_NONE -1
#define CALC_DISTANCE 0
#define CALC_GROUPDIST 1
#define CALC_RELATIONSHIP 2

#define _FILE_OFFSET_BITS 64
#define MAX_THREADS 64
#define MAX_THREADS_P1 65
#define PEDBUFBASE 256
#define FNAMESIZE 2048
#define MALLOC_DEFAULT_MB 2176
#define IDLENGTH 16
// size of generic text line load buffer.  .ped lines can of course be longer
#define MAXLINELEN 131072
#if __LP64__
#define BMULTIPLEX 64
#define IMULTIPLEX 32
#else
#define BMULTIPLEX 32
#define IMULTIPLEX 16
#endif

#define DOUBLE_INT_MULT 1000000000.0
#define DOUBLE_INT_MULT_RECIP 0.000000001

const char info_str[] =
  "WDIST weighted genetic distance calculator, v0.3.2 (27 July 2012)\n"
  "Christopher Chang (chrchang523@gmail.com), BGI Cognitive Genomics Lab\n\n"
  "wdist <flags> {calculation}\n";
const char errstr_append[] = "\nRun 'wdist --help' for more information.\n";
const char errstr_map_format[] = "Error: Improperly formatted .map file.\n";
const char errstr_fam_format[] = "Error: Improperly formatted .fam file.\n";
const char errstr_ped_format[] = "Error: Improperly formatted .ped file.\n";
char tbuf[MAXLINELEN];
int iwt[256];
unsigned char* wkspace;
long long malloc_size_mb = MALLOC_DEFAULT_MB;
char** subst_argv = NULL;
// gsl_rng* rg;

// TODO:
// distance MAF histograms
// regression coefficients [delete-d jackknife errors]

int dispmsg(int retval) {
  switch(retval) {
  case RET_HELP:
    printf("%s%s", info_str, errstr_append);
    break;
  case RET_NOMEM:
    printf("Error: Out of memory.\n");
    break;
  case RET_WRITE_FAIL:
    printf("Error: File write failure.\n");
    break;
  case RET_READ_FAIL:
    printf("Error: File read failure.\n");
    break;
  case RET_MOREHELP:
    printf(
"%s\n"
"Supported flags:\n"
"  --script [fname] : Include command-line options from file.\n"
"  --file [prefix]  : Specify prefix for .ped and .map files (default 'wdist').\n"
"  --ped [filename] : Specify name of .ped file.\n"
"  --map [filename] : Specify name of .map file.\n"
"  --make-bed       : Make .bed, .bim, and .fam files.\n"
"  --no-fid         : .ped file does not contain column 1 (family ID).\n"
"  --no-parents     : .ped file does not contain columns 3-4 (parents).\n"
"  --no-sex         : .ped file does not contain column 5 (sex).\n"
"  --no-pheno       : .ped file does not contain column 6 (phenotype).\n"
"  --liability      : .ped file does contain liability (column 7).\n"
"  --bfile [prefix] : Specify .bed/.bim/.fam prefix (default 'wdist').\n"
"  --bed [filename] : Specify .bed file.\n"
"  --bim [filename] : Specify .bim file.\n"
"  --fam [filename] : Specify .fam file.\n"
"  --out [prefix]   : Specify prefix for output files (default 'wdist').\n"
"  --silent         : Suppress output to console.\n"
"  --pheno [fname]  : Specify alternate phenotype.\n"
"  --mpheno [col]   : Specify phenotype column number in --pheno file.\n"
"  --pheno-name [c] : If phenotype file has a header row, use column with the\n"
"                     given name.\n"
"  --prune          : Remove individuals with missing phenotypes.\n"
"  --1              : .ped file affection phenotypes are interpreted as\n"
"                     0 = unaffected, 1 = affected (instead of 0 = missing,\n"
"                     1 = unaffected, 2 = affected)\n"
// --map3 implicitly supported via autodetection
// --compound-genotypes automatically supported
"  --chr [num]      : Only consider markers on the given chromosome (1-22, X,\n"
"                     Y, XY, MT).\n"
"  --maf [val]      : Minor allele frequency minimum threshold.\n"
"  --geno [val]     : Maximum per-SNP missing.\n"
"  --mind [val]     : Maximum per-person missing.\n"
"  --hwe [val]      : Minimum Hardy-Weinberg disequilibrium p-value (exact).\n"
"                     This is checked after all other forms of filtering.\n"
"  --rseed [val]    : Set random number seed (relevant for missing genotypes).\n"
"  --memory [val]   : Size, in MB, of initial malloc attempt (default 2176).\n"
"  --threads [val]  : Maximum number of concurrent threads (default 2).\n"
"  --exponent [val] : When computing genetic distances, each marker has a\n"
"                     weight of (2q(1-q))^{-val}, where q is the observed MAF.\n\n"
"  --keep [filename]\n"
"  --remove [filename]\n"
"  --filter [filename] [val] : Keep/remove/filter individuals (see PLINK\n"
"                              documentation).\n"
"  --mfilter [col]           : Specify column number in --filter file.\n\n"
"  --missing-genotype [char]     : Code for missing genotype (default 0).\n"
"  --missing-phenotype [val]     : Code for missing phenotype (default -9).\n"
"  --make-pheno [filename] [val] : Specify binary phenotype, where cases have\n"
"                                  the given value.  If the value is '*', all\n"
"                                  individuals present in the phenotype file\n"
"                                  are affected (and other individuals in the\n"
"                                  .ped/.fam are unaffected).\n"
"  --tail-pheno [Ltop] [Hbottom] : Form 'low' (<= Ltop, unaffected) and 'high'\n"
"                                  (> Hbottom, affected) groups from continuous\n"
"                                  phenotype data.  (Central phenotype values\n"
"                                  are treated as missing.)\n\n"
"Supported calculations:\n"
"  --distance <--square0>\n"
"    Outputs a lower-triangular table of (weighted) genetic distances to\n"
"    {output prefix}.dist, and a list of the corresponding family/individual\n"
"    IDs to {output prefix}.dist.id.\n"
"    The first row of the .dist file contains a single number describing the\n"
"    distance between the first two genotypes, the second row has the {genotype\n"
"    1-genotype 3} and {genotype 2-genotype 3} distances in that order, etc.\n\n"
"    If modified by the --square0 flag, a square matrix is written instead\n"
"    (filled out with zeroes).\n\n"
"  --make-rel <--square0>\n"
"    Outputs a lower-triangular (or filled out with zeroes to square, with\n"
"    --square0) relationship matrix to {output prefix}.rel, and corresponding\n"
"    IDs to {output prefix}.rel.id.\n"
"  --make-grm \n"
"    Writes the relationship matrix in GCTA's .grm format instead (except\n"
"    without gzipping).\n\n"
// "  --groupdist [d] [iters]\n"
// "    Two-group genetic distance analysis, using delete-d jackknife with the\n"
// "    requested number of iterations.  Binary phenotype required.\n\n"
         , info_str);
    break;
  }
  if (subst_argv) {
    free(subst_argv);
  }
  return retval;
}

// (copied from PLINK helper.cpp)
//
// This function implements an exact SNP test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton

double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
{
  
  if (obs_hom1 + obs_hom2 + obs_hets == 0 ) return 1;
  
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
  int genotypes   = obs_hets + obs_homc + obs_homr;

  double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
  if (het_probs == NULL) {
  //   error("Internal error: SNP-HWE: Unable to allocate array" );
    return -1.0;
  }
  
  int i;
  for (i = 0; i <= rare_copies; i++)
    het_probs[i] = 0.0;

  /* start at midpoint */
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
  
  /* check to ensure that midpoint and rare alleles have same parity */
  if ((rare_copies & 1) ^ (mid & 1))
    mid++;
  
  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
    {
      het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
	/ (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
      sum += het_probs[curr_hets - 2];

      /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
      curr_homr++;
      curr_homc++;
    }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
    {
      het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
	/((curr_hets + 2.0) * (curr_hets + 1.0));
      sum += het_probs[curr_hets + 2];
      
      /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
      curr_homr--;
      curr_homc--;
    }
  
  for (i = 0; i <= rare_copies; i++)
    het_probs[i] /= sum;

  double p_hwe = 0.0;
  /*  p-value calculation for p_hwe  */
  for (i = 0; i <= rare_copies; i++)
    {
      if (het_probs[i] > het_probs[obs_hets])
	continue;
      p_hwe += het_probs[i];
    }
   
  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  free(het_probs);

  return p_hwe;
}

inline int is_space_or_eoln(char cc) {
  return ((cc == ' ') || (cc == '\t') || (cc == '\n') || (cc == '\0'));
}

int strlen_se(char* ss) {
  int val = 0;
  while (!is_space_or_eoln(*ss++)) {
    val++;
  }
  return val;
}

char* id_buf = NULL;

int strcmp_casted(const void* s1, const void* s2) {
  return strcmp((char*)s1, (char*)s2);
}

int strcmp_deref(const void* s1, const void* s2) {
  return strcmp(*(char**)s1, *(char**)s2);
}

// alas, qsort_r not available on some Linux distributions
int qsort_ext(char* main_arr, int arr_length, int item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, int secondary_item_len) {
  char* proxy_arr;
  int proxy_len = secondary_item_len + sizeof(void*);
  int ii;
  if (!arr_length) {
    return 0;
  }
  if (proxy_len < item_length) {
    proxy_len = item_length;
  }
  proxy_arr = (char*)malloc(arr_length * proxy_len);
  if (!proxy_arr) {
    return -1;
  }
  for (ii = 0; ii < arr_length; ii++) {
    *(char**)(&(proxy_arr[ii * proxy_len])) = &(main_arr[ii * item_length]);
    memcpy(&(proxy_arr[ii * proxy_len + sizeof(void*)]), &(secondary_arr[ii * secondary_item_len]), secondary_item_len);
  }
  qsort(proxy_arr, arr_length, proxy_len, comparator_deref);
  for (ii = 0; ii < arr_length; ii++) {
    memcpy(&(secondary_arr[ii * secondary_item_len]), &(proxy_arr[ii * proxy_len + sizeof(void*)]), secondary_item_len);
    memcpy(&(proxy_arr[ii * proxy_len]), *(char**)(&(proxy_arr[ii * proxy_len])), item_length);
  }
  for (ii = 0; ii < arr_length; ii++) {
    memcpy(&(main_arr[ii * item_length]), &(proxy_arr[ii * proxy_len]), item_length);
  }
  free(proxy_arr);
  return 0;
}

int bsearch_str(char* lptr, int max_id_len, int min_idx, int max_idx) {
  int mid_idx;
  int ii;
  if (max_idx < min_idx) {
    return -1;
  }
  mid_idx = (min_idx + max_idx) / 2;
  ii = strcmp(id_buf, &(lptr[mid_idx * max_id_len]));
  if (ii) {
    if (ii < 0) {
      return bsearch_str(lptr, max_id_len, min_idx, mid_idx - 1);
    } else {
      return bsearch_str(lptr, max_id_len, mid_idx + 1, max_idx);
    }
  } else {
    return mid_idx;
  }
}

int bsearch_fam_indiv(char* lptr, int max_id_len, int filter_lines, char* fam_id, char* indiv_id) {
  int ii;
  int jj;
  char* sptr;
  if (!filter_lines) {
    return -1;
  }
  ii = 0;
  jj = 0;
  sptr = fam_id;
  while (!is_space_or_eoln(*sptr)) {
    sptr++;
    ii++;
  }
  sptr = indiv_id;
  while (!is_space_or_eoln(*sptr)) {
    sptr++;
    jj++;
  }
  if (ii + jj + 2 > max_id_len) {
    return -1;
  }
  memcpy(id_buf, fam_id, ii);
  id_buf[ii] = ' ';
  memcpy(&(id_buf[ii + 1]), indiv_id, jj);
  id_buf[ii + jj + 1] = '\0';
  return bsearch_str(lptr, max_id_len, 0, filter_lines - 1);
}

inline int is_contained(char* lptr, int max_id_len, int filter_lines, char* fam_id, char* indiv_id) {
  return (bsearch_fam_indiv(lptr, max_id_len, filter_lines, fam_id, indiv_id) != -1);
}

int marker_code(char* sptr) {
  if (sptr[1] > ' ') {
    if (sptr[0] == 'X') {
      return 25; // XY
    } else if (sptr[0] == 'M') {
      return 26; // MT
    } else {
      return ((sptr[1] - '0') * 10 + (sptr[0] - '0'));
    }
  } else if (*sptr < 'X') {
    return (sptr[0] - '0');
  } else if (*sptr == 'X') {
    return 23;
  } else if (*sptr == 'Y') {
    return 24;
  } else {
    return 0;
  }
}

typedef char id_string[IDLENGTH];

void cur_item(char* buf, char* sptr) {
  // no bounds-checking
  while ((*sptr != ' ') && (*sptr != '\t')) {
    *buf++ = *sptr++;
  }
  *buf = '\0';
}

char* next_item(char* sptr) {
  if (!sptr) {
    return NULL;
  }
  while ((*sptr != ' ') && (*sptr != '\t')) {
    if (!(*sptr)) {
      return NULL;
    }
    sptr++;
  }
  while ((*sptr == ' ') || (*sptr == '\t')) {
    sptr++;
  }
  return sptr;
}

void fill_weights(double* weights, double* mafs, double exponent) {
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  double wtarr[16];
  double twt[16];
  double* wt;
  for (ii = 0; ii < 16; ii += 1) {
    wtarr[ii] = pow(2.0 * mafs[ii] * (1.0 - mafs[ii]), -exponent);
  }
  for (nn = 0; nn < 4; nn++) {
    wt = &(wtarr[4 * nn]);
    for (ii = 0; ii < 4; ii += 1) {
      if ((ii % 4 == 1) || (ii % 4 == 2)) {
	twt[0] = wt[3];
      } else if (ii % 4 == 3) {
	twt[0] = wt[3] * 2;
      } else {
	twt[0] = 0;
      }
      for (jj = 0; jj < 4; jj += 1) {
	if ((jj % 4 == 1) || (jj % 4 == 2)) {
	  twt[1] = twt[0] + wt[2];
	} else if (jj % 4 == 3) {
	  twt[1] = twt[0] + 2 * wt[2];
	} else {
	  twt[1] = twt[0];
	}
	for (kk = 0; kk < 4; kk += 1) {
	  if ((kk % 4 == 1) || (kk % 4 == 2)) {
	    twt[2] = twt[1] + wt[1];
	  } else if (kk % 4 == 3) {
	    twt[2] = twt[1] + 2 * wt[1];
	  } else {
	    twt[2] = twt[1];
	  }
	  for (mm = 0; mm < 4; mm += 1) {
	    if ((mm % 4 == 1) || (mm % 4 == 2)) {
	      *weights++ = twt[2] + wt[0];
	    } else if (mm % 4 == 3) {
	      *weights++ = twt[2] + 2 * wt[0];
	    } else {
	      *weights++ = twt[2];
	    }
	  }
	}
      }
    }
  }
}
inline int round_dbl(double dd) {
  if (dd >= 0.0) {
    return (int)(dd + 0.5);
  } else {
    return (int)(dd - 0.5);
  }
}


void fill_weights_r(int* weights, double* mafs, int subjs) {
  int ii;
  int jj;
  int kk;
  int wtarr[256];
  int twt;
  int* wt;
  double mean;
  double mean_m1;
  double mean_m2;
  double mult;
  memset(wtarr, 0, 256 * sizeof(int));
  for (ii = 0; ii < 16; ii += 1) {
    if ((mafs[ii] > 0.00000001) && (mafs[ii] < 0.99999999)) {
      if (mafs[ii] < 0.50000001) {
	mean = 2.0 * mafs[ii];
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
	mult = 1.0 / (2.0 * mafs[ii] * (1.0 - mafs[ii]) * (double)subjs);
	wtarr[ii * 16] = round_dbl(DOUBLE_INT_MULT * mean * mean * mult);
	wtarr[ii * 16 + 2] = round_dbl(DOUBLE_INT_MULT * mean_m1 * mean * mult);
	wtarr[ii * 16 + 3] = round_dbl(DOUBLE_INT_MULT * mean_m2 * mean * mult);
	wtarr[ii * 16 + 4] = round_dbl(DOUBLE_INT_MULT * mean_m1 * mean_m1 * mult);
	wtarr[ii * 16 + 5] = round_dbl(DOUBLE_INT_MULT * mean_m2 * mean_m1 * mult);
	wtarr[ii * 16 + 6] = round_dbl(DOUBLE_INT_MULT * mean_m2 * mean_m2 * mult);
      } else {
	mean = 2.0 * (1.0 - mafs[ii]);
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
	mult = 1.0 / (2.0 * mafs[ii] * (1.0 - mafs[ii]) * (double)subjs);
	wtarr[ii * 16] = round_dbl(DOUBLE_INT_MULT * mean_m2 * mean_m2 * mult);
	wtarr[ii * 16 + 2] = round_dbl(DOUBLE_INT_MULT * mean_m2 * mean_m1 * mult);
	wtarr[ii * 16 + 3] = round_dbl(DOUBLE_INT_MULT * mean_m2 * mean * mult);
	wtarr[ii * 16 + 4] = round_dbl(DOUBLE_INT_MULT * mean_m1 * mean_m1 * mult);
	wtarr[ii * 16 + 5] = round_dbl(DOUBLE_INT_MULT * mean_m1 * mean * mult);
	wtarr[ii * 16 + 6] = round_dbl(DOUBLE_INT_MULT * mean * mean * mult);
      }
    }
  }
  for (kk = 0; kk < 8; kk++) {
    wt = &(wtarr[32 * kk]);
    for (ii = 0; ii < 16; ii += 1) {
      twt = wt[ii + 16];
      for (jj = 0; jj < 16; jj += 1) {
        *weights++ = twt + wt[jj];
      }
    }
  }
}

void exclude(unsigned char* exclude_arr, int loc, int* exclude_ct) {
  int maj = loc / 8;
  int min = 1 << (loc % 8);
  if (!(exclude_arr[maj] & min)) {
    exclude_arr[maj] |= min;
    *exclude_ct += 1;
  }
}

inline int excluded(unsigned char* exclude_arr, int loc) {
  return exclude_arr[loc / 8] & (1 << (loc % 8));
}

int is_missing(char* bufptr, int missing_pheno, int missing_pheno_len, int affection_01) {
  if ((atoi(bufptr) == missing_pheno) && is_space_or_eoln(bufptr[missing_pheno_len])) {
    return 1;
  } else if ((!affection_01) && (*bufptr == '0') && is_space_or_eoln(bufptr[1])) {
    return 1;
  }
  return 0;
}

int eval_affection(char* bufptr, int missing_pheno, int missing_pheno_len, int affection_01) {
  if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
    return 1;
  } else if (((*bufptr == '0') || (*bufptr == '1') || ((*bufptr == '2') && (!affection_01))) && is_space_or_eoln(bufptr[1])) {
    return 1;
  }
  return 0;
}

void collapse_phenoc(char* pheno_c, unsigned char* person_exclude, int ped_linect) {
  int ii = 0;
  int jj;
  while ((ii < ped_linect) && (!excluded(person_exclude, ii))) {
    ii++;
  }
  jj = ii;
  while (++ii < ped_linect) {
    if (!excluded(person_exclude, ii)) {
      pheno_c[jj++] = pheno_c[ii];
    }
  }
}

void collapse_phenod(double* pheno_d, unsigned char* person_exclude, int ped_linect) {
  int ii = 0;
  int jj;
  while ((ii < ped_linect) && (!excluded(person_exclude, ii))) {
    ii++;
  }
  jj = ii;
  while (++ii < ped_linect) {
    if (!excluded(person_exclude, ii)) {
      pheno_d[jj++] = pheno_d[ii];
    }
  }
}

void pick_d(unsigned char* cbuf, int ct, int dd) {
  int ii;
  int jj;
  int kk;
  memset(cbuf, 0, ct);
  kk = RAND_MAX % ct;
  for (ii = 0; ii < dd; ii++) {
    do {
      do {
        jj = rand();
      } while (jj < kk);
      jj %= ct;
    } while (cbuf[jj]);
    cbuf[jj] = 1;
  }
}

int* idists;
int* idists2;
double* dists;
unsigned char* ped_geno = NULL;
unsigned long* glptr;
double weights[1024];
int* weights_i = NULL;
int thread_start[MAX_THREADS_P1];

void incr_dists_i(int* idists, unsigned long* geno, int tidx) {
  unsigned long* glptr;
  unsigned long* glptr2;
  unsigned long ulii;
  unsigned long uljj;
  int ii;
  for (ii = thread_start[tidx]; ii < thread_start[tidx + 1]; ii++) {
    glptr = geno;
    glptr2 = &(geno[ii]);
    ulii = *glptr2;
    while (glptr < glptr2) {
      uljj = *glptr++ ^ ulii;
#if __LP64__
      *idists += iwt[uljj >> 56] + iwt[(uljj >> 48) & 255] + iwt[(uljj >> 40) & 255] + iwt[(uljj >> 32) & 255] + iwt[(uljj >> 24) & 255] + iwt[(uljj >> 16) & 255] + iwt[(uljj >> 8) & 255] + iwt[uljj & 255];
#else
      *idists += iwt[uljj >> 24] + iwt[(uljj >> 16) & 255] + iwt[(uljj >> 8) & 255] + iwt[uljj & 255];
#endif
      idists++;
    }
  }
}

void* calc_idist_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  incr_dists_i(&(idists[(ii * (ii - 1)) / 2]), (unsigned long*)ped_geno, (int)tidx);
  return NULL;
}

void incr_dists(double* dists, unsigned int* geno, int tidx) {
  unsigned int* giptr;
  unsigned int* giptr2;
  unsigned int uii;
  unsigned int ujj;
  double* weights1 = &(weights[256]);
  double* weights2 = &(weights[512]);
  double* weights3 = &(weights[768]);
  int ii;
  for (ii = thread_start[tidx]; ii < thread_start[tidx + 1]; ii++) {
    giptr = geno;
    giptr2 = &(geno[ii]);
    uii = *giptr2;
    while (giptr < giptr2) {
      ujj = *giptr++ ^ uii;
      *dists += weights3[ujj >> 24] + weights2[(ujj >> 16) & 255] + weights1[(ujj >> 8) & 255] + weights[ujj & 255];
      dists++;
    }
  }
}

void* calc_dist_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  incr_dists(&(dists[(ii * (ii - 1)) / 2]), (unsigned int*)ped_geno, (int)tidx);
  return NULL;
}

void incr_dists_r(int* dists, unsigned long* geno, int tidx, int* weights) {
  unsigned long* glptr;
  unsigned long* glptr2;
  unsigned long ulii;
  unsigned long uljj;
  int* weights1 = &(weights[256]);
  int* weights2 = &(weights[512]);
  int* weights3 = &(weights[768]);
#if __LP64__
  int* weights4 = &(weights[1024]);
  int* weights5 = &(weights[1280]);
  int* weights6 = &(weights[1536]);
  int* weights7 = &(weights[1792]);
#endif
  int ii;
  for (ii = thread_start[tidx]; ii < thread_start[tidx + 1]; ii++) {
    glptr = geno;
    glptr2 = &(geno[ii]);
    ulii = *glptr2;
    while (glptr <= glptr2) {
      uljj = *glptr++ + ulii;
#if __LP64__
      *dists += weights7[uljj >> 56] + weights6[(uljj >> 48) & 255] + weights5[(uljj >> 40) & 255] + weights4[(uljj >> 32) & 255] + weights3[(uljj >> 24) & 255] + weights2[(uljj >> 16) & 255] + weights1[(uljj >> 8) & 255] + weights[uljj & 255];
#else
      *dists += weights3[uljj >> 24] + weights2[(uljj >> 16) & 255] + weights1[(uljj >> 8) & 255] + weights[uljj & 255];
#endif
      dists++;
    }
  }
}

void* calc_rel_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  incr_dists_r(&(idists[(ii * (ii + 1)) / 2]), (unsigned long*)ped_geno, (int)tidx, weights_i);
  return NULL;
}

void incr_dists_rm(int* idists, unsigned long* genom, int tidx) {
  // count missing
  unsigned long* glptr;
  unsigned long* glptr2;
  unsigned long ulii;
  unsigned long uljj;
  int ii;
  for (ii = thread_start[tidx]; ii < thread_start[tidx + 1]; ii++) {
    glptr = genom;
    glptr2 = &(genom[ii]);
    ulii = *glptr2;
    if (ulii) {
      while (glptr <= glptr2) {
	uljj = *glptr++ | ulii;
#if __LP64__
	*idists += iwt[uljj >> 56] + iwt[(uljj >> 48) & 255] + iwt[(uljj >> 40) & 255] + iwt[(uljj >> 32) & 255] + iwt[(uljj >> 24) & 255] + iwt[(uljj >> 16) & 255] + iwt[(uljj >> 8) & 255] + iwt[uljj & 255];
#else
	*idists += iwt[uljj >> 24] + iwt[(uljj >> 16) & 255] + iwt[(uljj >> 8) & 255] + iwt[uljj & 255];
#endif
	idists++;
      }
    } else {
      while (glptr <= glptr2) {
        uljj = *glptr++;
        if (uljj) {
#if __LP64__
	  *idists += iwt[uljj >> 56] + iwt[(uljj >> 48) & 255] + iwt[(uljj >> 40) & 255] + iwt[(uljj >> 32) & 255] + iwt[(uljj >> 24) & 255] + iwt[(uljj >> 16) & 255] + iwt[(uljj >> 8) & 255] + iwt[uljj & 255];
#else
	  *idists += iwt[uljj >> 24] + iwt[(uljj >> 16) & 255] + iwt[(uljj >> 8) & 255] + iwt[uljj & 255];
#endif
        }
        idists++;
      }
    }
  }
}

void* calc_relm_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  incr_dists_rm(&(idists2[(ii * (ii + 1)) / 2]), (unsigned long*)glptr, (int)tidx);
  return NULL;
}

void triangle_fill(int* target_arr, int ct, int pieces, int start) {
  long long ct_tr = (long long)ct;
  long long vv = (long long)start;
  long long modif = vv * 2 - 1;
  long long cur_mult = 0;
  int cur_piece = 1;
  target_arr[0] = start;
  target_arr[pieces] = ct;
  ct_tr = (ct_tr * (ct_tr + modif)) / pieces;
  while (cur_piece < pieces) {
    cur_mult += ct_tr;
    while ((vv * (vv + modif)) < cur_mult) {
      vv++;
    }
    target_arr[cur_piece++] = (int)vv;
  }
}

int wdist(char* pedname, char* mapname, char* famname, char* phenoname, char* filtername, char* makepheno_str, int filter_type, char* filterval, int mfilter_col, int make_bed, int ped_col_1, int ped_col_34, int ped_col_5, int ped_col_6, int ped_col_7, char missing_geno, int missing_pheno, int mpheno_col, char* phenoname_str, int prune, int affection_01, int chr_num, int thread_ct, double exponent, double min_maf, double geno_thresh, double mind_thresh, double hwe_thresh, int tail_pheno, double tail_bottom, double tail_top, char* outname, int calculation_type, int calc_param_1, int iters) {
  FILE* pedfile = NULL;
  FILE* mapfile = NULL;
  FILE* famfile = NULL;
  FILE* outfile = NULL;
  FILE* phenofile = NULL;
  FILE* filterfile = NULL;
  FILE* bedtmpfile = NULL;
  FILE* bimtmpfile = NULL;
  FILE* famtmpfile = NULL;
  int map_linect = 0;
  int map_linect4;
  char* outname_end;
  unsigned char* pedbuf = NULL;
  unsigned char* marker_exclude = NULL;
  long long* line_locs = NULL;
  int max_people;
  int pedbuflen;
  int ped_recalc_len = 0;
  char* fgets_return;
  int ped_linect = 0;
  int person_exclude_ct = 0;
  int ped_linect4;
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int oo;
  int pp;
  unsigned int uii;
  unsigned int ujj;
  unsigned long ulii;
  unsigned long uljj;
  double dxx;
  double dyy;
  double dzz;
  double dhh_sum;
  double dhh_ssd;
  double dhl_sum;
  double dhl_ssd;
  double dll_sum;
  double dll_ssd;
  int low_ct;
  int high_ct;
  long long low_tot_i;
  long long lh_tot_i;
  long long high_tot_i;
  double low_tot;
  double lh_tot;
  double high_tot;
  long long llxx;
  // int* marker_chrom = NULL;
  // id_string* marker_id = NULL;
  // int* marker_pos = NULL;
  char* marker_alleles = NULL;
  char* marker_alleles_tmp = NULL;
  int* marker_allele_cts = NULL;
  int* person_missing_cts = NULL;
  char* bufptr;
  unsigned char* person_exclude = NULL;
  int retval;
  int map_cols = 3;
  int affection = 0;
  double* phenor_d = NULL;
  double* pheno_d = NULL;
  char* phenor_c = NULL;
  char* pheno_c = NULL;
  char* person_id = NULL;
  int max_person_id_len = 4;
  unsigned char* gptr;
  unsigned int* giptr;
  unsigned long* glptr2;
  char* cptr;
  char* cptr2;
  char* cptr3;
  int* iwptr;
  int* iptr;
  long long dists_alloc = 0;
  long long ped_geno_size;
  int geno_window_size;
  char cc;
  double* dist_ptr;
  long long last_tell;
  int binary_files = 0;
  int maf_int_thresh;
  int geno_int_thresh;
  int mind_int_thresh;
  int marker_exclude_ct = 0;
  int pheno_lines = 0;
  int makepheno_all = 0;
  int filter_lines = 0;
  int snp_major = 0;
  int max_pid_len = 0;
  int max_id_len = 4;
  char* pid_list = NULL;
  char* id_list = NULL;
  int rand_thresh_buf[IMULTIPLEX];
  double maf_buf[BMULTIPLEX];
  double missing_phenod = (double)missing_pheno;
  int missing_pheno_len = 1;
  int* hwe_ll;
  int* hwe_lh;
  int* hwe_hh;
  int* hwe_u_ll;
  int* hwe_u_lh;
  int* hwe_u_hh;
  int* hwe_a_ll;
  int* hwe_a_lh;
  int* hwe_a_hh;
  int hwe_lli;
  int hwe_lhi;
  int hwe_hhi;
  int hwe_u_lli;
  int hwe_u_lhi;
  int hwe_u_hhi;
  int hwe_a_lli;
  int hwe_a_lhi;
  int hwe_a_hhi;
  int multiplex;
  int bin_pheno;
  pthread_t threads[MAX_THREADS];

  ii = missing_pheno;
  if (ii < 0) {
    ii = -ii;
    missing_pheno_len++;
  }
  while (ii > 9) {
    ii /= 10;
    missing_pheno_len++;
  }

  retval = RET_OPENFAIL;
  pedfile = fopen(pedname, "r");
  if (!pedfile) {
    printf("Error: Failed to open %s.\n", pedname);
    goto wdist_ret_0;
  }
  mapfile = fopen(mapname, "r");
  if (!mapfile) {
    printf("Error: Failed to open %s.\n", mapname);
    goto wdist_ret_0;
  }
  if (famname[0]) {
    binary_files = 1;
    famfile = fopen(famname, "r");
    if (!famfile) {
      printf("Error: Failed to open %s.\n", famname);
      goto wdist_ret_0;
    }
  }
  if (phenoname[0]) {
    phenofile = fopen(phenoname, "r");
    if (!phenofile) {
      printf("Error: Failed to open %s.\n", famname);
      goto wdist_ret_0;
    }
  }
  outname_end = outname;
  while (*outname_end) {
    outname_end++;
  }
  retval = RET_NOMEM;

  // note that this actually only allocates up to half of the main workspace to
  // the distance array, because of the difference between n(n-1)/2 and n^2.
  if (exponent == 0.0) {
    max_people = (int)sqrt((malloc_size_mb * 1048576) / sizeof(int));
  } else {
    max_people = (int)sqrt((malloc_size_mb * 1048576) / sizeof(double));
  }

  // ----- .map/.bim load, first pass -----
  tbuf[MAXLINELEN - 6] = ' ';
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN - 5, mapfile) != NULL) {
    if (tbuf[0] > ' ') {
      if (!map_linect) {
	bufptr = next_item(tbuf);
	bufptr = next_item(bufptr);
	bufptr = next_item(bufptr);
        if (binary_files) {
          bufptr = next_item(bufptr);
          bufptr = next_item(bufptr);
        }
	if (!bufptr) {
	  retval = RET_INVALID_FORMAT;
	  printf(errstr_map_format);
	  goto wdist_ret_1;
	}
	if (*bufptr > ' ') {
	  map_cols = 4;
	}
      }
      map_linect += 1;
    }
    if (!tbuf[MAXLINELEN - 6]) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Excessively long line in .map/.bim file (max %d chars).\n", MAXLINELEN - 8);
      goto wdist_ret_1;
    }
  }
  if (!map_linect) {
    retval = RET_INVALID_FORMAT;
    printf("Error: No markers in .map/.bim file.");
    goto wdist_ret_0;
  }
  rewind(mapfile);
  marker_exclude = (unsigned char*)malloc(((map_linect + 7) / 8) * sizeof(char));
  if (!marker_exclude) {
    goto wdist_ret_1;
  }
  memset(marker_exclude, 0, ((map_linect + 7) / 8) * sizeof(char));
  if (binary_files) {
    map_linect4 = (map_linect + 3) / 4;
    pedbuflen = map_linect4;
  } else {
    pedbuflen = map_linect * 4 + PEDBUFBASE;
  }
  pedbuf = (unsigned char*)malloc(pedbuflen * sizeof(char));
  if (!pedbuf) {
    goto wdist_ret_1;
  }
  // marker_chrom = (int*)malloc(map_linect * sizeof(int));
  // if (!marker_chrom) {
  //   goto wdist_ret_0;
  // }
  // marker_id = (id_string*)malloc(map_linect * sizeof(id_string));
  // if (!marker_id) {
  //   goto wdist_ret_1;
  // }
  // marker_pos = (int*)malloc(map_linect * sizeof(int));
  // if (!marker_pos) {
  //   goto wdist_ret_2;
  // }

  // ----- .map/.bim load, second pass -----
  for (ii = 0; ii < map_linect; ii += 1) {
    do {
      fgets(tbuf, MAXLINELEN, mapfile);
    } while (tbuf[0] <= ' ');
    // marker_chrom[ii] = marker_code(tbuf);
    if (chr_num && (marker_code(tbuf) != chr_num)) {
      exclude(marker_exclude, ii, &marker_exclude_ct);
    } else {
      bufptr = next_item(tbuf);
      // cur_item(marker_id[ii], bufptr);
      bufptr = next_item(bufptr);
      if (map_cols == 4) {
	bufptr = next_item(bufptr);
      }
      if (*bufptr == '-') {
	if (binary_files) {
	  retval = RET_INVALID_FORMAT;
	  printf("Error: Negative marker position in .bim file.\n");
	  goto wdist_ret_1;
	}
	exclude(marker_exclude, ii, &marker_exclude_ct);
      }
    }
    // marker_pos[ii] = atoi(bufptr);
  }
  if (marker_exclude_ct) {
    printf("%d markers loaded (after excluding %d).\n", map_linect - marker_exclude_ct, marker_exclude_ct);
  } else {
    printf("%d markers loaded.\n", map_linect - marker_exclude_ct);
  }

  // ----- phenotype file load, first pass -----
  if (phenofile) {
    if (makepheno_str) {
      if ((makepheno_str[0] == '*') && (makepheno_str[1] == '\0')) {
        makepheno_all = 1;
      } else {
        kk = strlen(makepheno_str);
      }
    }
    mm = tail_pheno; // is affection determined from file yet?
    nn = 0; // existence of header row determined from file yet?
    while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
      if (*tbuf == '\n') {
        continue;
      }
      bufptr = tbuf;
      ii = 0;
      jj = 0;
      while (!is_space_or_eoln(*bufptr)) {
        ii++;
        bufptr++;
      }
      if ((*bufptr == '\n') || (*bufptr == '\0')) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Improperly formatted phenotype file.\n");
        goto wdist_ret_1;
      }
      while ((*bufptr == ' ') || (*bufptr == '\t')) {
        bufptr++;
      }
      cptr = bufptr;
      while (!is_space_or_eoln(*bufptr)) {
        jj++;
        bufptr++;
      }
      if (ii + jj + 2 > max_pid_len) {
        max_pid_len = ii + jj + 2;
      }
      if (!nn) {
        nn = 1;
        if ((ii == 3) && (jj == 3) && (!memcmp("FID", tbuf, 3)) && (!memcmp("IID", cptr, 3))) {
          pheno_lines = -1;
          if (phenoname_str) {
            jj = strlen(phenoname_str);
            do {
              bufptr = next_item(bufptr);
              if ((!bufptr) || (*bufptr == '\n') || (*bufptr == '\0')) {
		retval = RET_INVALID_FORMAT;
		printf("Error: --pheno-name column not found.\n");
		goto wdist_ret_1;
              }
              mpheno_col++;
              ii = 0;
              cptr = bufptr;
              while (!is_space_or_eoln(*bufptr)) {
                bufptr++;
                ii++;
              }
	      if (ii != jj) {
		continue;
	      }
            } while (memcmp(cptr, phenoname_str, ii));
            phenoname_str = NULL;
          }
        } else if (phenoname_str) {
          retval = RET_INVALID_FORMAT;
          printf("Error: No header row in phenotype file.\n");
          goto wdist_ret_1;
        }
        if (!mpheno_col) {
          mpheno_col = 1;
        }
      }
      if (pheno_lines == -1) {
        pheno_lines++;
      } else {
        for (ii = 0; ii < mpheno_col; ii++) {
          bufptr = next_item(bufptr);
        }
        if (!bufptr || (*bufptr == '\n') || (*bufptr == '\0')) {
	  retval = RET_INVALID_FORMAT;
	  printf("Error: Improperly formatted phenotype file.\n");
	  goto wdist_ret_1;
        }
        if (makepheno_str && (!makepheno_all)) {
          if (!strncmp(makepheno_str, bufptr, kk)) {
            if (is_space_or_eoln(bufptr[kk])) {
              pheno_lines++;
            }
          }
        } else {
          if (!mm) {
            affection = eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01);
            mm = 1;
          }
          if (affection) {
            if ((!prune) || (!is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01))) {
              pheno_lines++;
            }
          } else {
	    if (sscanf(bufptr, "%lf", &dxx) != 1) {
	      retval = RET_INVALID_FORMAT;
	      printf("Error: Improperly formatted phenotype file.\n");
	      goto wdist_ret_1;
	    }
            if ((!prune) || ((dxx != missing_phenod) && ((!tail_pheno) || ((dxx <= tail_bottom) || (dxx > tail_top))))) {
              pheno_lines++;
            }
          }
        }
      }
      if (!tbuf[MAXLINELEN - 1]) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Excessively long line in phenotype file (max %d chars).\n", MAXLINELEN - 3);
        goto wdist_ret_1;
      }
    }
    id_buf = (char*)malloc(max_pid_len * sizeof(char));
    if (!id_buf) {
      goto wdist_ret_1;
    }
    if (pheno_lines) {
      pid_list = (char*)malloc(max_pid_len * pheno_lines * sizeof(char));
      if (!pid_list) {
	goto wdist_ret_1;
      }
      if (affection || tail_pheno) {
	phenor_c = (char*)malloc(pheno_lines * sizeof(char));
	if (!phenor_c) {
	  goto wdist_ret_1;
	}
      } else if (!makepheno_str) {
	phenor_d = (double*)malloc(pheno_lines * sizeof(double));
	if (!phenor_d) {
	  goto wdist_ret_1;
	}
      }

      rewind(phenofile);
      ii = 0;
      // ----- phenotype file load, second pass -----
      while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
	if (ii == pheno_lines) {
	  break;
	}
	if (*tbuf == '\n') {
	  continue;
	}
	cptr = &(pid_list[ii * max_pid_len]);
	bufptr = tbuf;
	while ((*bufptr != ' ') && (*bufptr != '\t')) {
	  *cptr++ = *bufptr++;
	}
	*cptr++ = ' ';
	while ((*bufptr == ' ') || (*bufptr == '\t')) {
	  bufptr++;
	}
	while (!is_space_or_eoln(*bufptr)) {
	  *cptr++ = *bufptr++;
	}
	for (jj = 0; jj < mpheno_col; jj++) {
	  bufptr = next_item(bufptr);
	}
	*cptr = '\0';
	if (makepheno_str) {
	  if (makepheno_all) {
	    ii++;
	  } else {
	    if (!strncmp(makepheno_str, bufptr, kk)) {
	      if (is_space_or_eoln(bufptr[kk])) {
		ii++;
	      }
	    }
	  }
	} else if (affection) {
	  if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
	    cc = -1;
	  } else if (affection_01) {
	    cc = *bufptr - '0';
	  } else {
	    cc = *bufptr - '1';
	  }
	  if ((cc != -1) || (!prune)) {
	    phenor_c[ii++] = cc;
	  }
	} else {
	  sscanf(bufptr, "%lf", &dxx);
	  if (!(prune && (dxx == missing_phenod))) {
	    if (tail_pheno) {
	      if (dxx == missing_phenod) {
		phenor_c[ii++] = -1;
	      } else if (dxx <= tail_bottom) {
		phenor_c[ii++] = 0;
	      } else if (dxx > tail_top) {
		phenor_c[ii++] = 1;
	      } else if (!prune) {
		phenor_c[ii++] = -1;
	      }
	    } else {
	      phenor_d[ii++] = dxx;
	    }
	  }
	}
      }
      if (affection || tail_pheno) {
	if (qsort_ext(pid_list, pheno_lines, max_pid_len, strcmp_deref, phenor_c, sizeof(char)) == -1) {
	  goto wdist_ret_1;
        }
      } else if (makepheno_str) {
        qsort(pid_list, pheno_lines, max_pid_len, strcmp_casted);
      } else {
	if (qsort_ext(pid_list, pheno_lines, max_pid_len, strcmp_deref, (char*)phenor_d, sizeof(double)) == -1) {
          goto wdist_ret_1;
        }
      }
    } else {
      printf("Note: No valid entries in phenotype file.\n");
    }
    fclose(phenofile);
    phenofile = NULL;
  }

  if (filter_type) {
    filterfile = fopen(filtername, "r");
    if (!filterfile) {
      printf("Error: Failed to open %s.\n", filtername);
      goto wdist_ret_1;
    }
  }
  // ----- filter load, first pass -----
  if (filter_type) {
    if (filter_type == FILTER_CUSTOM) {
      jj = strlen(filterval);
      if (!mfilter_col) {
        mfilter_col = 1;
      }
    }
    while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
      if (*tbuf == '\n') {
        continue;
      }
      bufptr = tbuf;
      ii = 0;
      while (!is_space_or_eoln(*bufptr)) {
        ii++;
        bufptr++;
      }
      if ((*bufptr == '\n') || (*bufptr == '\0')) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Improperly formatted filter file.\n");
        goto wdist_ret_1;
      }
      while ((*bufptr == ' ') || (*bufptr == '\t')) {
        bufptr++;
      }
      while (!is_space_or_eoln(*bufptr)) {
        ii++;
        bufptr++;
      }
      ii += 2;
      if (ii > max_id_len) {
        max_id_len = ii;
      }
      if (filter_type == FILTER_CUSTOM) {
        for (kk = 0; kk < mfilter_col; kk++) {
          bufptr = next_item(bufptr);
        }
	if ((!bufptr) || (*bufptr == '\n') || (*bufptr == '\0')) {
	  retval = RET_INVALID_FORMAT;
	  printf("Error: Improperly formatted filter file.\n");
	  goto wdist_ret_1;
        }
        if (!strncmp(filterval, bufptr, jj)) {
          filter_lines += 1;
        }
      } else {
        filter_lines += 1;
      }
      if (!tbuf[MAXLINELEN - 1]) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Excessively long line in filter file (max %d chars).\n", MAXLINELEN - 3);
        goto wdist_ret_1;
      }
    }
    if (max_pid_len && (max_id_len > max_pid_len)) {
      free(id_buf);
    }
    id_buf = (char*)malloc(max_id_len * sizeof(char));
    if (!id_buf) {
      goto wdist_ret_1;
    }
    if (filter_lines) {
      id_list = (char*)malloc(max_id_len * filter_lines * sizeof(char));
      if (!id_list) {
	goto wdist_ret_1;
      }
      rewind(filterfile);
      ii = 0;

      // ----- filter load, second pass -----
      while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
	if (ii == filter_lines) {
	  break;
	}
	if (*tbuf == '\n') {
	  continue;
	}
	cptr = &(id_list[ii * max_id_len]);
	bufptr = tbuf;
	while ((*bufptr != ' ') && (*bufptr != '\t')) {
	  *cptr++ = *bufptr++;
	}
	*cptr++ = ' ';
	while ((*bufptr == ' ') || (*bufptr == '\t')) {
	  bufptr++;
	}
	while (!is_space_or_eoln(*bufptr)) {
	  *cptr++ = *bufptr++;
	}
	if (filter_type == FILTER_CUSTOM) {
	  for (kk = 0; kk < mfilter_col; kk++) {
	    bufptr = next_item(bufptr);
	  }
	  if (!strncmp(filterval, bufptr, jj)) {
	    ii++;
	    *cptr = '\0';
	  }
	} else {
	  ii++;
	  *cptr = '\0';
	}
      }
      qsort(id_list, filter_lines, max_id_len, strcmp_casted);
    } else {
      printf("Note: No valid entries in filter file.\n");
    }
    fclose(filterfile);
    filterfile = NULL;
  }
  line_locs = (long long*)malloc(max_people * sizeof(long long));
  if (!line_locs) {
    goto wdist_ret_1;
  }
  mind_int_thresh = (int)(2.0 * mind_thresh * (map_linect - marker_exclude_ct));

  if (binary_files) {
    jj = 0; // line count in .bed
    // ----- .fam load, first pass -----
    while (fgets(tbuf, MAXLINELEN, famfile) != NULL) {
      if (tbuf[0] > ' ') {
        if (tbuf[0] != '#') {
	  bufptr = next_item(tbuf);
          if (!ped_col_1) {
            cptr = tbuf;
          } else {
            cptr = bufptr;
          }
          ii = strlen_se(tbuf) + strlen_se(cptr) + 2;
          if (ii > max_person_id_len) {
            max_person_id_len = ii;
          }
          if (filter_type) {
            ii = is_contained(id_list, max_id_len, filter_lines, tbuf, cptr);
            if (filter_type == FILTER_REMOVE) {
              ii = 1 - ii;
            }
          } else {
            ii = 1;
          }
          if (phenoname[0]) {
            if (ii) {
	      if (makepheno_str || (!prune)) {
		line_locs[ped_linect++] = jj;
	      } else {
		if (is_contained(pid_list, max_pid_len, pheno_lines, tbuf, cptr)) {
		  line_locs[ped_linect++] = jj;
		}
	      }
            }
          } else if (ii || (!ped_linect)) {
	    if (ped_col_1) {
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_34) {
	      bufptr = next_item(bufptr);
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_5) {
	      bufptr = next_item(bufptr);
	    }
	    if (!bufptr) {
	      retval = RET_INVALID_FORMAT;
	      printf(errstr_fam_format);
	      goto wdist_ret_1;
	    }
            if (!ped_linect) {
              affection = eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01);
            }
            if (ii) {
	      if (affection) {
                if ((!prune) || (!is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01))) {
                  line_locs[ped_linect++] = jj;
                }
	      } else {
                if (sscanf(bufptr, "%lf", &dxx) != 1) {
                  retval = RET_INVALID_FORMAT;
                  printf(errstr_fam_format);
                  goto wdist_ret_1;
                }
                if ((!prune) || ((dxx != missing_phenod) && ((!tail_pheno) || ((dxx <= tail_bottom) || (dxx > tail_top))))) {
                  line_locs[ped_linect++] = jj;
                }
	      }
            }
          }
          jj++;
        }
      }
      if (!tbuf[MAXLINELEN - 1]) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Excessively long line in .fam file (max %d chars).\n", MAXLINELEN - 3);
        goto wdist_ret_1;
      }
    }
    if (ped_linect < 2) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Less than two valid people in .fam file.\n");
      goto wdist_ret_1;
    } else if (ped_linect > max_people) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Too many people in .fam file for this algorithm.\n");
      goto wdist_ret_1;
    }

    bin_pheno = (makepheno_str || affection || tail_pheno);
    if (bin_pheno) {
      pheno_c = (char*)malloc(ped_linect * sizeof(char));
      if (!pheno_c) {
	goto wdist_ret_1;
      }
    } else {
      pheno_d = (double*)malloc(ped_linect * sizeof(double));
      if (!pheno_d) {
	goto wdist_ret_1;
      }
    }
    person_id = (char*)malloc(ped_linect * max_person_id_len * sizeof(char));
    if (!person_id) {
      goto wdist_ret_1;
    }

    maf_int_thresh = 2 * ped_linect - (int)((1.0 - min_maf) * ped_linect * 2.0);
    geno_int_thresh = 2 * ped_linect - (int)(geno_thresh * 2.0 * ped_linect);

    rewind(famfile);
    // ----- .fam load, second pass -----
    ii = 0; // raw line number
    jj = 0; // loaded lines
    while (fgets(tbuf, MAXLINELEN, famfile) != NULL) {
      if (tbuf[0] > ' ') {
        if (tbuf[0] != '#') {
	  if (jj == ped_linect) {
	    break;
	  }
	  if (line_locs[jj] > ii) {
	    ii++;
	    continue;
	  }

	  bufptr = next_item(tbuf);
	  if (!ped_col_1) {
	    cptr = tbuf;
	  } else {
	    cptr = bufptr;
	  }
	  kk = strlen_se(tbuf);
	  mm = strlen_se(cptr);
	  memcpy(&(person_id[jj * max_person_id_len]), tbuf, kk);
	  person_id[jj * max_person_id_len + kk] = '\t';
	  memcpy(&(person_id[jj * max_person_id_len + kk + 1]), cptr, mm);
	  person_id[jj * max_person_id_len + kk + mm + 1] = '\0';
          if (phenoname[0]) {
	    kk = bsearch_fam_indiv(pid_list, max_pid_len, pheno_lines, tbuf, cptr);
            if (makepheno_str) {
              if (kk == -1) {
                pheno_c[jj] = 0;
              } else {
                pheno_c[jj] = 1;
              }
            } else if (affection || tail_pheno) {
              if (kk == -1) {
                pheno_c[jj] = -1;
              } else {
                pheno_c[jj] = phenor_c[kk];
              }
            } else if (kk == -1) {
              pheno_d[jj] = missing_phenod;
            } else {
              pheno_d[jj] = phenor_d[kk];
            }
          } else {
	    if (ped_col_1) {
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_34) {
	      bufptr = next_item(bufptr);
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_5) {
	      bufptr = next_item(bufptr);
	    }
            if (affection) {
              if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
                pheno_c[jj] = -1;
              } else if (affection_01) {
                pheno_c[jj] = *bufptr - '0';
              } else {
                pheno_c[jj] = *bufptr - '1';
              }
            } else {
              sscanf(bufptr, "%lf", &dxx);
              if (tail_pheno) {
                if (dxx == missing_phenod) {
                  pheno_c[jj] = -1;
                } if (dxx <= tail_bottom) {
                  pheno_c[jj] = 0;
                } else if (dxx > tail_top) {
                  pheno_c[jj] = 1;
                } else {
                  pheno_c[jj] = -1;
                }
              } else {
                pheno_d[jj] = dxx;
              }
            }
          }
	  ii++;
	  jj++;
        }
      }
    }

    if (fread(pedbuf, 1, 3, pedfile) < 3) {
      retval = RET_READ_FAIL;
      goto wdist_ret_1;
    }
    if ((pedbuf[0] != 'l') || (pedbuf[1] != '\x1b')) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Invalid or pre-v1.00 BED file.\n");
      goto wdist_ret_1;
    }
    marker_allele_cts = (int*)malloc(map_linect * 2 * sizeof(int));
    if (!marker_allele_cts) {
      goto wdist_ret_1;
    }
    memset(marker_allele_cts, 0, map_linect * 2 * sizeof(int));

    person_exclude = (unsigned char*)malloc(((ped_linect + 7) / 8) * sizeof(char));
    if (!person_exclude) {
      goto wdist_ret_1;
    }
    memset(person_exclude, 0, ((ped_linect + 7) / 8) * sizeof(char));
    if (pedbuf[2] == '\0') {
      // ----- individual-major .bed load, first pass -----
      for (ii = 0; ii < ped_linect; ii += 1) {
        if (fread(pedbuf, 1, map_linect4, pedfile) < map_linect4) {
          retval = RET_READ_FAIL;
          goto wdist_ret_1;
        }
        gptr = pedbuf - 1;
        mm = 0; // missing
        for (jj = 0; jj < map_linect; jj++) {
          kk = jj % 4;
          if (!kk) {
            nn = *(++gptr);
          } else {
            nn >>= 2;
          }
          oo = nn & 3;
          if (oo) {
            if (oo == 2) {
              marker_allele_cts[jj * 2] += 1;
              marker_allele_cts[jj * 2 + 1] += 1;
            } else if (oo == 3) {
              marker_allele_cts[jj * 2] += 2;
            } else {
              mm++;
            }
          } else {
            marker_allele_cts[jj * 2 + 1] += 2;
          }
        }
        if (mm > mind_int_thresh) {
          exclude(person_exclude, ii, &person_exclude_ct);
        }
      }
      for (ii = 0; ii < map_linect; ii++) {
        if ((marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1] < geno_int_thresh) || (marker_allele_cts[ii * 2 + 1] < maf_int_thresh) || (marker_allele_cts[ii * 2] < maf_int_thresh)) {
          exclude(marker_exclude, ii, &marker_exclude_ct);
        }
      }
      if (hwe_thresh > 0.0) {
        // ----- individual-major .bed load, second pass -----
	fseeko(pedfile, 3, SEEK_SET);
	hwe_ll = (int*)wkspace;
	hwe_lh = (int*)(&wkspace[map_linect * sizeof(int)]);
	hwe_hh = (int*)(&wkspace[map_linect * 2 * sizeof(int)]);
        hwe_u_ll = (int*)(&wkspace[map_linect * 3 * sizeof(int)]);
        hwe_u_lh = (int*)(&wkspace[map_linect * 4 * sizeof(int)]);
        hwe_u_hh = (int*)(&wkspace[map_linect * 5 * sizeof(int)]);
        hwe_a_ll = (int*)(&wkspace[map_linect * 6 * sizeof(int)]);
        hwe_a_lh = (int*)(&wkspace[map_linect * 7 * sizeof(int)]);
        hwe_a_hh = (int*)(&wkspace[map_linect * 8 * sizeof(int)]);
	memset(wkspace, 0, map_linect * 9 * sizeof(int));
	for (ii = 0; ii < ped_linect; ii += 1) {
          if (excluded(person_exclude, ii)) {
            fseeko(pedfile, (off_t)map_linect4, SEEK_CUR);
            continue;
          }
	  if (fread(pedbuf, 1, map_linect4, pedfile) < map_linect4) {
	    retval = RET_READ_FAIL;
	    goto wdist_ret_1;
	  }
	  gptr = pedbuf - 1;
	  for (jj = 0; jj < map_linect; jj++) {
	    kk = jj % 4;
	    if (!kk) {
	      nn = *(++gptr);
	    } else {
	      nn >>= 2;
	    }
	    oo = nn & 3;
	    if (oo) {
	      if (oo == 2) {
                hwe_lh[jj] += 1;
                if (bin_pheno) {
                  if (pheno_c[jj] == 0) {
                    hwe_u_lh[jj] += 1;
                  } else if (pheno_c[jj] == 1) {
                    hwe_a_lh[jj] += 1;
                  }
                }
	      } else if (oo == 3) {
                hwe_hh[jj] += 1;
                if (bin_pheno) {
                  if (pheno_c[jj] == 0) {
                    hwe_u_hh[jj] += 1;
                  } else if (pheno_c[jj] == 1) {
                    hwe_a_hh[jj] += 1;
                  }
                }
	      }
	    } else {
              hwe_ll[jj] += 1;
	      if (bin_pheno) {
		if (pheno_c[jj] == 0) {
		  hwe_u_ll[jj] += 1;
		} else if (pheno_c[jj] == 1) {
		  hwe_a_ll[jj] += 1;
		}
	      }
	    }
	  }
	}
        for (ii = 0; ii < map_linect; ii++) {
          if (!excluded(marker_exclude, ii)) {
	    if (SNPHWE(hwe_lh[ii], hwe_ll[ii], hwe_hh[ii]) < hwe_thresh) {
	      exclude(marker_exclude, ii, &marker_exclude_ct);
	    } else if (bin_pheno) {
              if ((SNPHWE(hwe_u_lh[ii], hwe_u_ll[ii], hwe_u_hh[ii]) < hwe_thresh) || (SNPHWE(hwe_a_lh[ii], hwe_a_ll[ii], hwe_a_hh[ii]) < hwe_thresh)) {
                exclude(marker_exclude, ii, &marker_exclude_ct);
              }
            }
          }
        }
      }
    } else {
      // ----- snp-major .bed load, first pass -----
      snp_major = 1;
      person_missing_cts = (int*)malloc(ped_linect * sizeof(int));
      if (!person_missing_cts) {
        goto wdist_ret_1;
      }
      memset(person_missing_cts, 0, ped_linect * sizeof(int));
      ped_linect4 = (ped_linect + 3) / 4;
      if (pedbuflen < (4 * ped_linect4)) {
        free(pedbuf);
        pedbuflen = 4 * ped_linect4;
        pedbuf = (unsigned char*)malloc(pedbuflen * sizeof(char));
        if (!pedbuf) {
          goto wdist_ret_2;
        }
      }
      for (ii = 0; ii < map_linect; ii++) {
        if (fread(pedbuf, 1, ped_linect4, pedfile) < ped_linect4) {
          retval = RET_READ_FAIL;
          goto wdist_ret_2;
        }
        mm = 0; // major allele ct
        nn = 0; // minor allele ct
        gptr = pedbuf - 1;
        for (jj = 0; jj < ped_linect; jj++) {
          kk = jj % 4;
          if (!kk) {
            oo = *(++gptr);
          } else {
            oo >>= 2;
          }
          pp = oo & 3;
          if (pp) {
            if (pp == 2) {
              mm++;
              nn++;
            } else if (pp == 3) {
              mm += 2;
            } else {
              person_missing_cts[jj] += 1;
            }
          } else {
            nn += 2;
          }
        }
        if ((mm + nn < geno_int_thresh) || (mm < maf_int_thresh) || (nn < maf_int_thresh)) {
          exclude(marker_exclude, ii, &marker_exclude_ct);
        } else {
          marker_allele_cts[ii * 2] = mm;
          marker_allele_cts[ii * 2 + 1] = nn;
        }
      }
      for (ii = 0; ii < ped_linect; ii++) {
        if (person_missing_cts[ii] > mind_int_thresh) {
	  exclude(person_exclude, ii, &person_exclude_ct);
        }
      }
      if (hwe_thresh > 0.0) {
        fseeko(pedfile, 3, SEEK_SET);
        // ----- .snp-major .bed load, second pass -----
        for (ii = 0; ii < map_linect; ii++) {
          if (excluded(marker_exclude, ii)) {
            fseeko(pedfile, (off_t)ped_linect4, SEEK_CUR);
            continue;
          }
          if (fread(pedbuf, 1, ped_linect4, pedfile) < ped_linect4) {
            retval = RET_READ_FAIL;
            goto wdist_ret_2;
          }
          gptr = pedbuf - 1;
          hwe_lli = 0;
          hwe_lhi = 0;
          hwe_hhi = 0;
          hwe_u_lli = 0;
          hwe_u_lhi = 0;
          hwe_u_hhi = 0;
          hwe_a_lli = 0;
          hwe_a_lhi = 0;
          hwe_a_hhi = 0;
	  for (jj = 0; jj < ped_linect; jj++) {
	    kk = jj % 4;
	    if (!kk) {
	      oo = *(++gptr);
	    } else {
	      oo >>= 2;
	    }
            if (!excluded(person_exclude, jj)) {
	      pp = oo & 3;
	      if (pp) {
		if (pp == 2) {
		  hwe_lhi += 1;
		  if (bin_pheno) {
		    if (pheno_c[jj] == 0) {
		      hwe_u_lhi += 1;
		    } else if (pheno_c[jj] == 1) {
		      hwe_a_lhi += 1;
		    }
		  }
		} else if (pp == 3) {
		  hwe_hhi += 1;
		  if (bin_pheno) {
		    if (pheno_c[jj] == 0) {
		      hwe_u_hhi += 1;
		    } else if (pheno_c[jj] == 1) {
		      hwe_a_hhi += 1;
		    }
		  }
		}
	      } else {
		hwe_lli += 1;
		if (bin_pheno) {
		  if (pheno_c[jj] == 0) {
		    hwe_u_lli += 1;
		  } else if (pheno_c[jj] == 1) {
		    hwe_a_lli += 1;
		  }
		}
	      }
            }
	  }          
	  if (SNPHWE(hwe_lhi, hwe_lli, hwe_hhi) < hwe_thresh) {
	    exclude(marker_exclude, ii, &marker_exclude_ct);
	  } else if (bin_pheno) {
            if ((SNPHWE(hwe_u_lhi, hwe_u_lli, hwe_u_hhi) < hwe_thresh) || (SNPHWE(hwe_a_lhi, hwe_a_lli, hwe_a_hhi) < hwe_thresh)) {
              exclude(marker_exclude, ii, &marker_exclude_ct);
            }
          }
        }
      }
    }
  } else {
    pedbuf[pedbuflen - 1] = ' ';
    last_tell = 0;
    // define exclude counters
    // ----- .ped load, first pass -----
    while (fgets((char*)pedbuf, pedbuflen, pedfile) != NULL) {
      if (pedbuf[0] > ' ') {
	if (pedbuf[0] != '#') {
	  bufptr = next_item((char*)pedbuf);
          if (filter_type) {
	    if (!ped_col_1) {
	      cptr = (char*)pedbuf;
	    } else {
	      cptr = bufptr;
	    }
            ii = is_contained(id_list, max_id_len, filter_lines, (char*)pedbuf, cptr);
            if (filter_type == FILTER_REMOVE) {
              ii = 1 - ii;
            }
          } else {
            ii = 1;
          }
          if (phenoname[0]) {
            if (ii) {
	      if (makepheno_str || (!prune)) {
		line_locs[ped_linect++] = last_tell;
	      } else {
		if (is_contained(pid_list, max_pid_len, pheno_lines, (char*)pedbuf, cptr)) {
		  line_locs[ped_linect++] = last_tell;
		}
	      }
            }
          } else if (ii || (!ped_linect)) {
	    if (ped_col_1) { // actually column 2
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_34) {
	      bufptr = next_item(bufptr);
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_5) {
	      bufptr = next_item(bufptr);
	    }
	    if (!bufptr) {
	      retval = RET_INVALID_FORMAT;
	      printf(errstr_ped_format);
	      goto wdist_ret_1;
	    }
            llxx = last_tell + (bufptr - (char*)pedbuf);
	    if (!ped_linect) {
              affection = eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01);
	    }
	    if (affection) {
              if ((!prune) || (!is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01))) {
                line_locs[ped_linect++] = llxx;
              }
	    } else {
	      if (sscanf(bufptr, "%lf", &dxx) != 1) {
		retval = RET_INVALID_FORMAT;
		printf(errstr_ped_format);
		goto wdist_ret_1;
	      }
              if ((!prune) || ((dxx != missing_phenod) && ((!tail_pheno) || ((dxx <= tail_bottom) || (dxx > tail_top))))) {
                line_locs[ped_linect++] = llxx;
              }
	    }
          }
	}
      }
      if (!pedbuf[pedbuflen - 1]) {
        pedbuf[pedbuflen - 1] = ' ';
        if (pedbuf[pedbuflen - 2] == '\n') {
	  last_tell = ftello(pedfile);
          continue;
        }
	ii = 0;
	do {
	  ii += pedbuflen - 1;
	  pedbuf[pedbuflen - 1] = ' ';
	  fgets_return = fgets((char*)pedbuf, pedbuflen, pedfile);
	} while (fgets_return && !pedbuf[pedbuflen - 1] && (pedbuf[pedbuflen - 2] != '\n'));
	ii += strlen((char*)pedbuf) + 1;
	if (ii > ped_recalc_len) {
	  ped_recalc_len = ii;
	}
      }
      last_tell = ftello(pedfile);
    }
    if (ped_linect < 2) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Less than two valid people in .ped file.\n");
      goto wdist_ret_1;
    } else if (ped_linect > max_people) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Too many people in .ped file for this algorithm.\n");
      goto wdist_ret_1;
    }
    rewind(pedfile);
    if (ped_recalc_len) {
      free(pedbuf);
      pedbuflen = ped_recalc_len;
      pedbuf = (unsigned char*)malloc(pedbuflen * sizeof(char));
      if (!pedbuf) {
	goto wdist_ret_1;
      }
    }

    maf_int_thresh = 2 * ped_linect - (int)((1.0 - min_maf) * ped_linect * 2.0);
    geno_int_thresh = 2 * ped_linect - (int)(geno_thresh * 2.0 * ped_linect);

    marker_alleles = (char*)malloc(map_linect * 4 * sizeof(char));
    if (!marker_alleles) {
      goto wdist_ret_1;
    }
    memset(marker_alleles, 0, map_linect * 4 * sizeof(char));
    marker_allele_cts = (int*)malloc(map_linect * 4 * sizeof(int));
    if (!marker_allele_cts) {
      goto wdist_ret_1;
    }
    memset(marker_allele_cts, 0, map_linect * 4 * sizeof(int));
    bin_pheno = (makepheno_str || affection || tail_pheno);
    if (bin_pheno) {
      pheno_c = (char*)malloc(ped_linect * sizeof(char));
      if (!pheno_c) {
	goto wdist_ret_1;
      }
    } else {
      pheno_d = (double*)malloc(ped_linect * sizeof(double));
      if (!pheno_d) {
	goto wdist_ret_1;
      }
    }
    person_exclude = (unsigned char*)malloc(((ped_linect + 7) / 8) * sizeof(char));
    if (!person_exclude) {
      goto wdist_ret_1;
    }
    memset(person_exclude, 0, ((ped_linect + 7) / 8) * sizeof(char));

    // ----- .ped load, second pass -----
    ped_geno = wkspace;
    for (ii = 0; ii < ped_linect; ii += 1) {
      fseeko(pedfile, line_locs[ii], SEEK_SET);
      fgets((char*)pedbuf, pedbuflen, pedfile);
      bufptr = (char*)pedbuf;
      if (phenoname[0]) {
        if (ped_col_1) {
          bufptr = next_item(bufptr);
        }
        jj = bsearch_fam_indiv(pid_list, max_pid_len, pheno_lines, (char*)pedbuf, bufptr);
        cc = 0;
        if (makepheno_str) {
          if (jj != -1) {
            cc = 1;
          }
          pheno_c[ii] = cc;
        } else if (affection || tail_pheno) {
          if (jj == -1) {
            pheno_c[ii] = -1;
          } else {
            pheno_c[ii] = phenor_c[jj];
          }
        } else if (jj == -1) {
          pheno_d[ii] = missing_phenod;
        } else {
          pheno_d[ii] = phenor_d[jj];
        }
        bufptr = next_item(bufptr);
        if (ped_col_34) {
          bufptr = next_item(bufptr);
          bufptr = next_item(bufptr);
        }
        if (ped_col_5) {
          bufptr = next_item(bufptr);
        }
      } else {
	if (affection) {
          if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
            pheno_c[ii] = -1;
          } else if (affection_01) {
            pheno_c[ii] = *bufptr - '0';
	  } else {
            pheno_c[ii] = *bufptr - '1';
	  }
	} else {
	  sscanf(bufptr, "%lf", &dxx);
          if (tail_pheno) {
            if (dxx == missing_phenod) {
              pheno_c[ii] = -1;
	    } else if (dxx <= tail_bottom) {
              pheno_c[ii] = 0;
            } else if (dxx > tail_top) {
              pheno_c[ii] = 1;
            } else {
              pheno_c[ii] = -1;
            }
          } else {
            pheno_d[ii] = dxx;
          }
	}
      }
      if (ped_col_6) {
        bufptr = next_item(bufptr);
      }
      if (ped_col_7) {
	bufptr = next_item(bufptr);
      }
      if (!bufptr) {
        retval = RET_INVALID_FORMAT;
        printf(errstr_ped_format);
        goto wdist_ret_1;
      }
      line_locs[ii] += (bufptr - (char*)pedbuf);
      mm = 0; // number of missing
      for (jj = 0; jj < map_linect; jj += 1) {
        for (kk = 0; kk < 2; kk++) {
          cc = *bufptr;
	  if (cc == '\0') {
	    retval = RET_INVALID_FORMAT;
	    printf(errstr_ped_format);
	    goto wdist_ret_1;
	  }
          if (!excluded(marker_exclude, jj)) {
	    if (cc == missing_geno) {
	      mm += 1;
	    } else if (cc == marker_alleles[jj * 4]) {
	      marker_allele_cts[jj * 4] += 1;
	    } else if (cc == marker_alleles[jj * 4 + 1]) {
	      marker_allele_cts[jj * 4 + 1] += 1;
	    } else if (cc == marker_alleles[jj * 4 + 2]) {
	      marker_allele_cts[jj * 4 + 2] += 1;
	    } else if (cc == marker_alleles[jj * 4 + 3]) {
	      marker_allele_cts[jj * 4 + 3] += 1;
	    } else if (marker_alleles[jj * 4 + 3]) {
	      retval = RET_INVALID_FORMAT;
	      printf("Error: More than 4 different allele types at marker %d.\n", jj + 1);
	      goto wdist_ret_1;
	    } else if (marker_alleles[jj * 4 + 2]) {
	      marker_alleles[jj * 4 + 3] = cc;
	      marker_allele_cts[jj * 4 + 3] = 1;
	    } else if (marker_alleles[jj * 4 + 1]) {
	      marker_alleles[jj * 4 + 2] = cc;
	      marker_allele_cts[jj * 4 + 2] = 1;
	    } else if (marker_alleles[jj * 4]) {
	      marker_alleles[jj * 4 + 1] = cc;
	      marker_allele_cts[jj * 4 + 1] = 1;
	    } else {
	      marker_alleles[jj * 4] = cc;
	      marker_allele_cts[jj * 4] = 1;
	    }
          }
	  bufptr++;
	  while ((*bufptr == ' ') || (*bufptr == '\t')) {
	    bufptr++;
	  }
        }
      }
      if (mm > mind_int_thresh) {
        exclude(person_exclude, ii, &person_exclude_ct);
      }
    }
    rewind(pedfile);
    for (ii = 0; ii < map_linect; ii += 1) {
      for (jj = 1; jj < 3; jj += 1) {
        for (kk = jj - 1; kk >= 0; kk -= 1) {
          if (marker_allele_cts[ii * 4 + kk] < marker_allele_cts[ii * 4 + kk + 1]) {
            cc = marker_alleles[ii * 4 + kk];
            marker_alleles[ii * 4 + kk] = marker_alleles[ii * 4 + kk + 1];
            marker_alleles[ii * 4 + kk + 1] = cc;
            mm = marker_allele_cts[ii * 4 + kk];
            marker_allele_cts[ii * 4 + kk] = marker_allele_cts[ii * 4 + kk + 1];
            marker_allele_cts[ii * 4 + kk + 1] = mm;
          }
        }
      }
      mm = marker_allele_cts[ii * 4] + marker_allele_cts[ii * 4 + 1] + marker_allele_cts[ii * 4 + 2] + marker_allele_cts[ii * 4 + 3];
      if (marker_allele_cts[ii * 4] + marker_allele_cts[ii * 4 + 1] < geno_int_thresh) {
        exclude(marker_exclude, ii, &marker_exclude_ct);
      } else if (marker_allele_cts[ii * 4 + 1] < maf_int_thresh) {
        exclude(marker_exclude, ii, &marker_exclude_ct);
      } else if (marker_allele_cts[ii * 4 + 1] == marker_allele_cts[ii * 4 + 2]) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Ambiguous minor allele at marker %d.\n", ii + 1);
        goto wdist_ret_1;
      }
    }
    if (person_exclude_ct > (ped_linect - 2)) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Too many people fail QC.\n");
      goto wdist_ret_1;
    } else if (marker_exclude_ct == map_linect) {
      retval = RET_INVALID_FORMAT;
      printf("Error: All markers fail QC.\n");
      goto wdist_ret_1;
    }
    if (hwe_thresh > 0.0) {
      // unfortunately, this requires a third pass...
      hwe_ll = (int*)ped_geno;
      hwe_lh = (int*)(&(ped_geno[map_linect * sizeof(int)]));
      hwe_hh = (int*)(&(ped_geno[map_linect * 2 * sizeof(int)]));
      hwe_u_ll = (int*)(&(ped_geno[map_linect * 3 * sizeof(int)]));
      hwe_u_lh = (int*)(&(ped_geno[map_linect * 4 * sizeof(int)]));
      hwe_u_hh = (int*)(&(ped_geno[map_linect * 5 * sizeof(int)]));
      hwe_a_ll = (int*)(&(ped_geno[map_linect * 6 * sizeof(int)]));
      hwe_a_lh = (int*)(&(ped_geno[map_linect * 7 * sizeof(int)]));
      hwe_a_hh = (int*)(&(ped_geno[map_linect * 8 * sizeof(int)]));
      memset(ped_geno, 0, map_linect * 9 * sizeof(int));

      for (ii = 0; ii < ped_linect; ii += 1) {
        if (excluded(person_exclude, ii)) {
          continue;
        }
        fseeko(pedfile, line_locs[ii], SEEK_SET);
        fgets((char*)pedbuf, pedbuflen, pedfile);
        bufptr = (char*)pedbuf;
        for (jj = 0; jj < map_linect; jj += 1) {
          if (excluded(marker_exclude, jj)) {
            continue;
          }
          mm = 0;
          for (kk = 0; kk < 2; kk++) {
            cc = *bufptr;
	    if (cc == marker_alleles[jj * 4 + 1]) {
	      mm += 1;
            } else if (cc != marker_alleles[jj * 4]) {
              mm = -2;
            }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
          }
          if (mm == 0) {
            hwe_ll[jj] += 1;
            if (bin_pheno) {
              if (pheno_c[jj] == 0) {
                hwe_u_ll[jj] += 1;
              } else if (pheno_c[jj] == 1) {
                hwe_a_ll[jj] += 1;
              }
            }
          } else if (mm == 1) {
            hwe_lh[jj] += 1;
            if (bin_pheno) {
              if (pheno_c[jj] == 0) {
                hwe_u_lh[jj] += 1;
              } else if (pheno_c[jj] == 1) {
                hwe_a_lh[jj] += 1;
              }
            }
          } else if (mm == 2) {
            hwe_hh[jj] += 1;
            if (bin_pheno) {
              if (pheno_c[jj] == 0) {
                hwe_u_hh[jj] += 1;
              } else if (pheno_c[jj] == 1) {
                hwe_a_hh[jj] += 1;
              }
            }
          }
        }
      }
      for (ii = 0; ii < map_linect; ii++) {
        if (!excluded(marker_exclude, ii)) {
          if (SNPHWE(hwe_lh[ii], hwe_ll[ii], hwe_hh[ii]) < hwe_thresh) {
            exclude(marker_exclude, ii, &marker_exclude_ct);
          } else if (bin_pheno) {
            if ((SNPHWE(hwe_u_lh[ii], hwe_u_ll[ii], hwe_u_hh[ii]) < hwe_thresh) || (SNPHWE(hwe_a_lh[ii], hwe_a_ll[ii], hwe_a_hh[ii]) < hwe_thresh)) {
              exclude(marker_exclude, ii, &marker_exclude_ct);
            }
          }
        }
      }
      rewind(pedfile);
    }

    ii = ((ped_linect - person_exclude_ct) * ((ped_linect - person_exclude_ct) - 1)) / 2;
    if (exponent == 0.0) {
      dists_alloc = ii * sizeof(int);
    } else {
      dists_alloc = ii * sizeof(double);
    }
    
    llxx = malloc_size_mb * 1048576 - dists_alloc;
    geno_window_size = llxx / (ped_linect - person_exclude_ct);
    map_linect4 = (map_linect - marker_exclude_ct + 3) / 4;
    pp = (ped_linect - person_exclude_ct + 3) / 4;
    if (make_bed && (pp <= ((malloc_size_mb * 1048576) / (map_linect - marker_exclude_ct)))) {
      snp_major = 1;
      printf("Writing binary files...");
      strcpy(outname_end, ".bed");
      bedtmpfile = fopen(outname, "wb");
      if (!bedtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_1;
      }
      strcpy(outname_end, ".bim");
      bimtmpfile = fopen(outname, "wb");
      if (!bimtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      strcpy(outname_end, ".fam");
      famtmpfile = fopen(outname, "wb");
      if (!famtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (3 != fwrite("l\x1b\x01", 1, 3, bedtmpfile)) {
	retval = RET_WRITE_FAIL;
        printf("\n");
	goto wdist_ret_2;
      }
      ped_geno = wkspace;
      memset(ped_geno, 0, pp * (map_linect - marker_exclude_ct));
      // ----- .ped (fourth pass) -> .fam + .bed conversion, snp-major -----
      rewind(pedfile);
      ii = 0;
      oo = 0; // final person index
      last_tell = 0;

      while (fgets((char*)pedbuf, pedbuflen, pedfile) != NULL) {
        if (ii == ped_linect) {
          break;
        }
        llxx = ftello(pedfile);
        if (llxx < line_locs[ii]) {
          last_tell = llxx;
          continue;
        }
        if (excluded(person_exclude, ii)) {
          last_tell = llxx;
          ii++;
          continue;
        }
        jj = (int)(line_locs[ii++] - last_tell);
        last_tell = llxx;
        pedbuf[jj - 1] = '\n';
        if (jj != fwrite(pedbuf, 1, jj, famtmpfile)) {
          retval = RET_WRITE_FAIL;
          printf("\n");
          goto wdist_ret_2;
        }
        bufptr = (char*)(&pedbuf[jj]);
        mm = 0; // final SNP index
	for (jj = 0; jj < map_linect; jj += 1) {
          if (excluded(marker_exclude, jj)) {
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
            continue;
          }
	  nn = 0;
	  for (kk = 0; kk < 2; kk++) {
	    cc = *bufptr;
	    if (cc == marker_alleles[jj * 4]) {
	      if (nn == 2) {
		nn = 3;
	      } else if (!nn) {
		nn = 2;
	      }
	    } else if (cc != marker_alleles[jj * 4 + 1]) {
              nn = 1;
	    }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
	  }
	  ped_geno[mm * pp + (oo / 4)] |= (nn << ((oo % 4) * 2));
          mm++;
	}
	oo++;
      }
      ii = pp * (map_linect - marker_exclude_ct);
      if (ii != fwrite(ped_geno, 1, ii, bedtmpfile)) {
	retval = RET_WRITE_FAIL;
        printf("\n");
	goto wdist_ret_2;
      }
    } else if (map_linect4 > geno_window_size) {
      if (make_bed) {
        printf("Writing individual-major binary files, since .ped is very large...");
      } else {
        printf(".ped file too large for direct read.  Converting to binary...");
      }
      if (make_bed) {
        strcpy(outname_end, ".bed");
      } else {
        strcpy(outname_end, ".bed.tmp");
      }
      bedtmpfile = fopen(outname, "wb");
      if (!bedtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_1;
      }
      if (make_bed) {
        strcpy(outname_end, ".bim");
      } else {
        strcpy(outname_end, ".bim.tmp");
      }
      bimtmpfile = fopen(outname, "wb");
      if (!bimtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (make_bed) {
        strcpy(outname_end, ".fam");
      } else {
        strcpy(outname_end, ".fam.tmp");
      }
      famtmpfile = fopen(outname, "wb");
      if (!famtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (3 != fwrite("l\x1b", 1, 3, bedtmpfile)) {
	retval = RET_WRITE_FAIL;
        printf("\n");
	goto wdist_ret_2;
      }
      // ----- .ped (fourth pass) -> .fam + .bed conversion, indiv-major -----
      rewind(pedfile);
      ii = 0; // line count
      last_tell = 0; // last file location
      while (fgets((char*)pedbuf, pedbuflen, pedfile) != NULL) {
        if (ii == ped_linect) {
          break;
        }
        llxx = ftello(pedfile);
        if (llxx < line_locs[ii]) {
          last_tell = llxx;
          continue;
        }
        if (excluded(person_exclude, ii)) {
          last_tell = llxx;
          ii++;
          continue;
        }
        jj = (int)(line_locs[ii++] - last_tell);
        last_tell = llxx;
        pedbuf[jj - 1] = '\n';
        if (jj != fwrite(pedbuf, 1, jj, famtmpfile)) {
          retval = RET_WRITE_FAIL;
          goto wdist_ret_2;
        }
        bufptr = (char*)(&pedbuf[jj]);
	memset(ped_geno, 0, map_linect4 * sizeof(char));
        gptr = ped_geno;
        mm = 0;
	for (jj = 0; jj < map_linect; jj += 1) {
          if (excluded(marker_exclude, jj)) {
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
            continue;
          }
	  nn = 0;
	  for (kk = 0; kk < 2; kk++) {
	    cc = *bufptr;
	    if (cc == marker_alleles[jj * 4]) {
	      if (nn == 2) {
		nn = 3;
	      } else if (!nn) {
		nn = 2;
	      }
	    } else if (cc != marker_alleles[jj * 4 + 1]) {
              nn = 1;
	    }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
	  }
	  *gptr |= (nn << (mm * 2));
          mm = (mm + 1) % 4;
	  if (mm == 0) {
	    gptr++;
	  }
	}
	if (map_linect4 != fwrite(ped_geno, 1, map_linect4, bedtmpfile)) {
	  retval = RET_WRITE_FAIL;
	  goto wdist_ret_2;
	}
      }
    }
    if (make_bed || (map_linect4 > geno_window_size)) {
      marker_alleles_tmp = (char*)malloc((map_linect - marker_exclude_ct) * 2 * sizeof(char));
      if (!marker_alleles_tmp) {
        goto wdist_ret_2;
      }
      cptr = marker_alleles_tmp;
      iwptr = (int*)malloc((map_linect - marker_exclude_ct) * 2 * sizeof(int));
      iptr = iwptr;
      for (ii = 0; ii < map_linect; ii += 1) {
        if (excluded(marker_exclude, ii)) {
          continue;
        }
        *cptr++ = marker_alleles[ii * 4];
        *cptr++ = marker_alleles[ii * 4 + 1];
        *iptr++ = marker_allele_cts[ii * 4];
        *iptr++ = marker_allele_cts[ii * 4 + 1];
      }
      free(marker_alleles);
      marker_alleles = marker_alleles_tmp;
      marker_alleles_tmp = NULL;
      marker_allele_cts = iwptr;

      // ----- .map -> .bim conversion -----
      rewind(mapfile);
      kk = 0;
      for (ii = 0; ii < map_linect; ii += 1) {
        do {
          fgets(tbuf, MAXLINELEN, mapfile);
        } while (tbuf[0] <= ' ');
        if (excluded(marker_exclude, ii)) {
          continue;
        }
        bufptr = next_item(tbuf);
        bufptr = next_item(bufptr);
        if (map_cols == 4) {
          bufptr = next_item(bufptr);
        }
        while (!is_space_or_eoln(*bufptr)) {
          bufptr++;
        }
        *bufptr++ = '\t';
        *bufptr++ = marker_alleles[kk * 2 + 1];
        *bufptr++ = '\t';
        *bufptr++ = marker_alleles[kk * 2];
        *bufptr++ = '\n';
        kk++;
        jj = (int)(bufptr - tbuf);
        if (jj != fwrite(tbuf, 1, jj, bimtmpfile)) {
	  retval = RET_WRITE_FAIL;
	  goto wdist_ret_2;
        }
      }
      printf(" done.\n");
      fclose(pedfile);
      fclose(bedtmpfile);
      fclose(bimtmpfile);
      fclose(famtmpfile);
      bedtmpfile = NULL;
      bimtmpfile = NULL;
      famtmpfile = NULL;
      if (make_bed) {
        strcpy(outname_end, ".bed");
      } else {
        strcpy(outname_end, ".bed.tmp");
      }
      pedfile = fopen(outname, "rb");
      if (!pedfile) {
        printf("Error: Failed to open %s.\n", outname);
      }
      binary_files = 1;
      if (pheno_c) {
        collapse_phenoc(pheno_c, person_exclude, ped_linect);
      } else if (pheno_d) {
        collapse_phenod(pheno_d, person_exclude, ped_linect);
      }
      map_linect -= marker_exclude_ct;
      marker_exclude_ct = 0;
      memset(marker_exclude, 0, (map_linect + 7) / 8);
      ped_linect -= person_exclude_ct;
      ped_linect4 = (ped_linect + 3) / 4;
      person_exclude_ct = 0;
      memset(person_exclude, 0, (ped_linect + 7) / 8);
    }
  }

  ii = ((ped_linect - person_exclude_ct) * ((ped_linect - person_exclude_ct) - 1)) / 2;
  if (calculation_type == CALC_RELATIONSHIP) {
    ii += ped_linect - person_exclude_ct;
  }
  if ((exponent == 0.0) || (calculation_type == CALC_RELATIONSHIP)) {
    dists_alloc = ii * sizeof(int);
    idists = (int*)wkspace;
    memset(idists, 0, dists_alloc);
  } else {
    dists_alloc = ii * sizeof(double);
    dists = (double*)wkspace;
    dist_ptr = dists;
    jj = ii;
    do {
      *dist_ptr++ = 0.0;
    } while (--jj);
  }
  if ((calculation_type == CALC_RELATIONSHIP) && (calc_param_1 == 2)) {
    idists2 = (int*)(&wkspace[dists_alloc]);
    memset(idists2, 0, ii * sizeof(int));
    ped_geno = &(wkspace[dists_alloc * 2]);
    ped_geno_size = malloc_size_mb * 1048576 - dists_alloc * 2;
  } else {
    ped_geno = &(wkspace[dists_alloc]);
    ped_geno_size = malloc_size_mb * 1048576 - dists_alloc;
  }

  if ((calculation_type == CALC_GROUPDIST) && (!bin_pheno)) {
    retval = RET_INVALID_CMDLINE;
    printf("Error: --groupdist calculation requires binary phenotype.\n");
    goto wdist_ret_2;
  }

  printf("%d markers and %d people pass filters and QC.\n", map_linect - marker_exclude_ct, ped_linect - person_exclude_ct);

  if (calculation_type == CALC_RELATIONSHIP) {
    if (binary_files && snp_major) {
      fseeko(pedfile, 3, SEEK_SET);
      ii = 0;
      pp = 0;
      if (calc_param_1 == 2) {
        glptr2 = (unsigned long*)&(ped_geno[(ped_linect - person_exclude_ct) * sizeof(long)]);
        gptr = &(ped_geno[(ped_linect - person_exclude_ct) * 2 * sizeof(long)]);
      } else {
        gptr = &(ped_geno[(ped_linect - person_exclude_ct) * sizeof(long)]);
      }
      weights_i = (int*)weights;
      triangle_fill(thread_start, ped_linect - person_exclude_ct, thread_ct, 0);
      // See later comments on CALC_DISTANCE.
      // The difference is, we have to use + instead of XOR here to distinguish
      // the cases, so we want to allow at least 3 bits per locus.  And given
      // that, the additional cost of supporting exact handling of missing
      // markers is smaller (3->4 bits), so we do it.
      while (pp < (map_linect - marker_exclude_ct)) {
        for (jj = 0; jj < BMULTIPLEX; jj++) {
          maf_buf[jj] = 0.0;
        }
        jj = 0;
        while ((jj < BMULTIPLEX) && (pp < (map_linect - marker_exclude_ct))) {
          while (excluded(marker_exclude, ii)) {
            ii++;
            fseeko(pedfile, (off_t)ped_linect4, SEEK_CUR);
          }
          if (fread(&(gptr[jj * ped_linect4]), 1, ped_linect4, pedfile) < ped_linect4) {
            retval = RET_READ_FAIL;
            goto wdist_ret_2;
          }
          maf_buf[jj] = ((double)marker_allele_cts[ii * 2]) / ((double)(marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1]));
          ii++;
          jj++;
          pp++;
        }
        if (jj < BMULTIPLEX) {
          memset(&(gptr[jj * ped_linect4]), 0, (BMULTIPLEX - jj) * ped_linect4);
        }
        if (calc_param_1 == 2) {
          memset(glptr2, 0, (ped_linect - person_exclude_ct) * sizeof(long));
        }

        for (nn = 0; nn < 4; nn++) {
	  glptr = (unsigned long*)ped_geno;
          oo = 0;
	  for (jj = 0; jj < ped_linect; jj++) {
	    if (!excluded(person_exclude, jj)) {
	      kk = (jj % 4) * 2;
	      ulii = 0;
	      for (mm = 0; mm < (BMULTIPLEX / 4); mm++) {
		uljj = (gptr[jj / 4 + ((nn * (BMULTIPLEX / 4)) + mm) * ped_linect4] >> kk) & 3;
		if (uljj == 1) {
		  uljj = 7;
                  if (calc_param_1 == 2) {
                    glptr2[oo] |= 1 << (nn * (BMULTIPLEX / 4) + mm);
		  }
		}
		ulii |= uljj << (mm * 4);
	      }
	      *glptr++ = ulii;
              oo++;
	    }
	  }
	  fill_weights_r(weights_i, &(maf_buf[nn * (BMULTIPLEX / 4)]), map_linect - marker_exclude_ct);
          for (ulii = 1; ulii < thread_ct; ulii++) {
	    if (pthread_create(&(threads[ulii - 1]), NULL, &calc_rel_thread, (void*)ulii)) {
	      printf("Error: Could not create thread.\n");
	      retval = RET_THREAD_CREATE_FAIL;
	      goto wdist_ret_2;
	    }
          }
	  incr_dists_r(idists, (unsigned long*)ped_geno, 0, weights_i);
	  for (oo = 0; oo < thread_ct - 1; oo++) {
	    pthread_join(threads[oo], NULL);
	  }
        }
        if (calc_param_1 == 2) {
          for (ulii = 1; ulii < thread_ct; ulii++) {
            if (pthread_create(&(threads[ulii - 1]), NULL, &calc_relm_thread, (void*)ulii)) {
	      printf("Error: Could not create thread.\n");
	      retval = RET_THREAD_CREATE_FAIL;
	      goto wdist_ret_2;
            }
          }
          incr_dists_rm(idists2, (unsigned long*)glptr, 0);
          for (oo = 0; oo < thread_ct - 1; oo++) {
            pthread_join(threads[oo], NULL);
          }
        }
        printf("\r%d markers complete.", pp);
        fflush(stdout);
      }
      printf("\rRelationship matrix calculation complete.\n");

      if (calc_param_1 == 2) {
        strcpy(outname_end, ".grm");
        outfile = fopen(outname, "w");
        if (!outfile) {
	  printf("Error: Failed to open %s.\n", outname);
	  goto wdist_ret_2;
        }
        iptr = idists;
        iwptr = idists2;
        for (ii = 0; ii < (ped_linect - person_exclude_ct); ii += 1) {
          for (jj = 0; jj <= ii; jj += 1) {
            fprintf(outfile, "%d\t%d\t%d\t%lf\n", ii + 1, jj + 1, map_linect - marker_exclude_ct - *iwptr++, ((double)(*iptr++)) * DOUBLE_INT_MULT_RECIP);
          }
        }
      } else {
	strcpy(outname_end, ".rel");
	outfile = fopen(outname, "w");
	if (!outfile) {
	  printf("Error: Failed to open %s.\n", outname);
	  goto wdist_ret_2;
	}
	iptr = idists;
	for (ii = 0; ii < (ped_linect - person_exclude_ct); ii += 1) {
	  for (jj = 0; jj <= ii; jj += 1) {
	    if (jj) {
	      fprintf(outfile, "\t%lf", ((double)(*iptr++)) * DOUBLE_INT_MULT_RECIP);
	    } else {
	      fprintf(outfile, "%lf", ((double)(*iptr++)) * DOUBLE_INT_MULT_RECIP);
	    }
	  }
	  if (calc_param_1) {
	    for (; jj < (ped_linect - person_exclude_ct); jj += 1) {
	      fprintf(outfile, "\t%lf", 0.0);
	    }
	  }
	  fprintf(outfile, "\n");
	}
      }
      retval = RET_SUCCESS;
      printf("Relationship matrix written to %s.\n", outname);
      fclose(outfile);
      strcpy(&(outname_end[4]), ".id");
      outfile = fopen(outname, "w");
      if (!outfile) {
	printf("Error: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      for (ii = 0; ii < ped_linect; ii += 1) {
	if (!excluded(person_exclude, ii)) {
	  fprintf(outfile, "%s\n", &(person_id[ii * max_person_id_len]));
	}
      }

    } else {
      printf("Error: Relationship calculation currently doesn't support individual-major\n.bed file.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto wdist_ret_2;
    }
    retval = RET_SUCCESS;
    goto wdist_ret_2;
  }

  if (exponent == 0.0) {
    multiplex = IMULTIPLEX;
  } else {
    multiplex = 16;
  }
  if (binary_files) {
    if (snp_major) {
      if (pedbuflen < multiplex * ped_linect4) {
        free(pedbuf);
        pedbuf = (unsigned char*)malloc(multiplex * ped_linect4 * sizeof(char));
        if (!pedbuf) {
          goto wdist_ret_2;
        }
      }
      fseeko(pedfile, 3, SEEK_SET);
      ii = 0; // current SNP index
      pp = 0; // after subtracting out excluded
      triangle_fill(thread_start, ped_linect - person_exclude_ct, thread_ct, 1);
      while (pp < (map_linect - marker_exclude_ct)) {
        for (jj = 0; jj < multiplex; jj++) {
          maf_buf[jj] = 0.5;
        }
        jj = 0; // actual SNPs read
        // Two key insights here:
        //
        // 1. Precalculate distances for all possible combinations of 4
        // markers.
        // Status of a single marker is stored in two bits, so four can be
        // stored in a byte while supporting bitwise operations.  [Weighted]
        // distance between two sets of four markers can be determined by
        // XORing the corresponding bytes and looking up the corresponding
        // array entry.  A pair of contiguous distance arrays is small enough
        // to fit in a single 4KB L1 cache entry.
        //
        // 2. Do #1 for 4-8 marker sets simultaneously, to further reduce the
        // number of reads/writes from main memory, and take advantage of the
        // speed of XORing 32- or 64-bit words.
        // Empirically, when the weighting exponent is nonzero (and thus
        // floating point arithmetic must be used), 4 sets appears to be no
        // worse than 8 when using a 64-bit processor.
        while ((jj < multiplex) && (pp < (map_linect - marker_exclude_ct))) {
          while (excluded(marker_exclude, ii)) {
            ii++;
            fseeko(pedfile, (off_t)ped_linect4, SEEK_CUR);
          }
          if (fread(&(pedbuf[jj * ped_linect4]), 1, ped_linect4, pedfile) < ped_linect4) {
            retval = RET_READ_FAIL;
            goto wdist_ret_2;
          }
          maf_buf[jj] = ((double)marker_allele_cts[ii * 2]) / ((double)(marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1]));
          rand_thresh_buf[jj] = (int)(maf_buf[jj] * RAND_MAX);
          ii++;
          jj++;
          pp++;
        }
        if (jj < multiplex) {
          memset(&(pedbuf[jj * ped_linect4]), 0, (multiplex - jj) * ped_linect4);
        }
        if (exponent == 0.0) {
	  glptr = (unsigned long*)ped_geno;
	  for (jj = 0; jj < ped_linect; jj++) {
	    if (!excluded(person_exclude, jj)) {
	      kk = (jj % 4) * 2;
	      ulii = 0;
	      for (mm = 0; mm < IMULTIPLEX; mm++) {
		uljj = (pedbuf[jj / 4 + mm * ped_linect4] >> kk) & 3;
		if (uljj == 1) {
		  if (rand() > rand_thresh_buf[mm]) {
                    if (rand() > rand_thresh_buf[mm]) {
                      uljj = 3;
                    } else {
		      uljj = 2;
                    }
		  } else if (rand() > rand_thresh_buf[mm]) {
		    uljj = 2;
		  } else {
                    uljj = 0;
                  }
		}
		ulii |= uljj << (mm * 2);
	      }
	      *glptr++ = ulii;
	    }
          }
          for (ulii = 1; ulii < thread_ct; ulii++) {
	    if (pthread_create(&(threads[ulii - 1]), NULL, &calc_idist_thread, (void*)ulii)) {
	      printf("Error: Could not create thread.\n");
	      retval = RET_THREAD_CREATE_FAIL;
	      goto wdist_ret_2;
	    }
          }
          incr_dists_i(idists, (unsigned long*)ped_geno, 0);
	  for (jj = 0; jj < thread_ct - 1; jj++) {
	    pthread_join(threads[jj], NULL);
	  }
	} else {
          giptr = (unsigned int*)ped_geno;
	  for (jj = 0; jj < ped_linect; jj++) {
	    if (!excluded(person_exclude, jj)) {
	      kk = (jj % 4) * 2;
	      uii = 0;
	      for (mm = 0; mm < 16; mm++) {
		ujj = (pedbuf[jj / 4 + mm * ped_linect4] >> kk) & 3;
		if (ujj == 1) {
		  if (rand() > rand_thresh_buf[mm]) {
                    if (rand() > rand_thresh_buf[mm]) {
                      ujj = 3;
                    } else {
		      ujj = 2;
		    }
                  } else if (rand() > rand_thresh_buf[mm]) {
		    ujj = 2;
		  } else {
                    ujj = 0;
                  }
		}
		uii |= ujj << (mm * 2);
	      }
	      *giptr++ = uii;
	    }
          }
          fill_weights(weights, maf_buf, exponent);
          for (ulii = 1; ulii < thread_ct; ulii++) {
	    if (pthread_create(&(threads[ulii - 1]), NULL, &calc_dist_thread, (void*)ulii)) {
	      printf("Error: Could not create thread.\n");
	      retval = RET_THREAD_CREATE_FAIL;
	      goto wdist_ret_2;
	    }
          }
          incr_dists(dists, (unsigned int*)ped_geno, 0);
	  for (jj = 0; jj < thread_ct - 1; jj++) {
	    pthread_join(threads[jj], NULL);
	  }
        }
        printf("\r%d markers complete.", pp);
        fflush(stdout);
      }
    } else {
      printf("indiv-major distance calculation not done.\n");
      retval = RET_SUCCESS;
      goto wdist_ret_2;
    }
  } else {
    printf("text distance calculation not done (use --make-bed).\n");
    retval = RET_SUCCESS;
    goto wdist_ret_2;
  }

  printf("\rDistance matrix calculation complete.\n");

  if (calculation_type == CALC_DISTANCE) {
    strcpy(outname_end, ".dist");
    outfile = fopen(outname, "w");
    if (!outfile) {
      printf("Error: Failed to open %s.\n", outname);
      goto wdist_ret_2;
    }

    if (exponent == 0.0) {
      iwptr = idists;
      if (calc_param_1) {
        fprintf(outfile, "0");
        for (ii = 1; ii < (ped_linect - person_exclude_ct); ii += 1) {
          fprintf(outfile, "\t0");
        }
        fprintf(outfile, "\n");
      }
      for (ii = 1; ii < (ped_linect - person_exclude_ct); ii += 1) {
        for (jj = 0; jj < ii; jj += 1) {
          if (jj) {
            fprintf(outfile, "\t%d", *iwptr++);
          } else {
            fprintf(outfile, "%d", *iwptr++);
          }
        }
	if (calc_param_1) {
          for (; jj < (ped_linect - person_exclude_ct); jj += 1) {
            fprintf(outfile, "\t0");
          }
	}
        fprintf(outfile, "\n");
      }
    } else {
      dist_ptr = dists;
      if (calc_param_1) {
        fprintf(outfile, "%lf", 0.0);
        for (ii = 1; ii < ped_linect; ii += 1) {
          fprintf(outfile, "\t%lf", 0.0);
        }
        fprintf(outfile, "\n");
      }
      for (ii = 1; ii < ped_linect; ii += 1) {
        for (jj = 0; jj < ii; jj += 1) {
          if (jj) {
            fprintf(outfile, "\t%lf", *dist_ptr++);
          } else {
            fprintf(outfile, "%lf", *dist_ptr++);
          }
        }
	if (calc_param_1) {
	  for (; jj < ped_linect; jj += 1) {
            fprintf(outfile, "\t%lf", 0.0);
          }
        }
        fprintf(outfile, "\n");
      }
    }
    retval = RET_SUCCESS;
    printf("Distances written to %s.\n", outname);
    fclose(outfile);
    strcpy(outname_end, ".dist.id");
    outfile = fopen(outname, "w");
    if (!outfile) {
      printf("Error: Failed to open %s.\n", outname);
      goto wdist_ret_2;
    }
    for (ii = 0; ii < ped_linect; ii += 1) {
      if (!excluded(person_exclude, ii)) {
        fprintf(outfile, "%s\n", &(person_id[ii * max_person_id_len]));
      }
    }
  } else if (calculation_type == CALC_GROUPDIST) {
    collapse_phenoc(pheno_c, person_exclude, ped_linect);
    low_ct = 0;
    high_ct = 0;
    if (exponent == 0.0) {
      low_tot_i = 0;
      lh_tot_i = 0;
      high_tot_i = 0;
      iptr = idists;
      for (ii = 0; ii < ped_linect; ii++) {
        cptr = pheno_c;
        cptr2 = &(pheno_c[ii]);
        if (*cptr2 == 1) {
          high_ct++;
          while (cptr < cptr2) {
            if (*cptr == 1) {
              high_tot_i += *iptr;
            } else if (!(*cptr)) {
              lh_tot_i += *iptr;
            }
            cptr++;
            iptr++;
	  }
        } else if (*cptr2 == 0) {
          low_ct++;
          while (cptr < cptr2) {
            if (*cptr == 1) {
              lh_tot_i += *iptr;
            } else if (!(*cptr)) {
              low_tot_i += *iptr;
            }
            cptr++;
            iptr++;
	  }
	} else {
          iptr += ii;
        }
      }
      low_tot = (double)low_tot_i;
      lh_tot = (double)lh_tot_i;
      high_tot = (double)high_tot_i;
    } else {
      low_tot = 0.0;
      lh_tot = 0.0;
      high_tot = 0.0;
      dist_ptr = dists;
      for (ii = 0; ii < ped_linect; ii++) {
        cptr = pheno_c;
        cptr2 = &(pheno_c[ii]);
        if (*cptr2 == 1) {
          high_ct++;
          while (cptr < cptr2) {
            if (*cptr == 1) {
              high_tot += *dist_ptr;
            } else if (!(*cptr)) {
              lh_tot += *dist_ptr;
            }
            cptr++;
            dist_ptr++;
	  }
        } else if (*cptr2 == 0) {
          low_ct++;
          while (cptr < cptr2) {
            if (*cptr == 1) {
              lh_tot += *dist_ptr;
            } else if (!(*cptr)) {
              low_tot += *dist_ptr;
            }
            cptr++;
            dist_ptr++;
	  }
	} else {
          dist_ptr += ii;
        }
      }
    }
    retval = RET_SUCCESS;
    printf("Group distance analysis (%d affected, %d unaffected):\n", high_ct, low_ct);
    if (high_ct < 2) {
      dxx = 0.0;
    } else {
      dxx = high_tot / (double)((high_ct * (high_ct - 1)) / 2);
    }
    if (!(high_ct * low_ct)) {
      dyy = 0.0;
    } else {
      dyy = lh_tot / (double)(high_ct * low_ct);
    }
    if (low_ct < 2) {
      dzz = 0.0;
    } else {
      dzz = low_tot / (double)((low_ct * (low_ct - 1)) / 2);
    }
    printf("  Avg dist between 2x affected             : %lf\n", dxx);
    printf("  Avg dist between affected and unaffected : %lf\n", dyy);
    printf("  Avg dist between 2x unaffected           : %lf\n\n", dzz);
    if (2 * calc_param_1 >= high_ct + low_ct) {
      printf("Delete-d jackknife skipped because d is too large.\n");
    } else {
      dhh_sum = 0.0;
      dhl_sum = 0.0;
      dll_sum = 0.0;
      dhh_ssd = 0.0;
      dhl_ssd = 0.0;
      dll_ssd = 0.0;
      nn = 0;
      for (jj = 0; jj < iters; jj++) {
        kk = 0;
        mm = 0;
        pick_d(ped_geno, high_ct + low_ct, calc_param_1);
	if (exponent == 0.0) {
	  low_tot_i = 0;
	  lh_tot_i = 0;
	  high_tot_i = 0;
	  iptr = idists;
	  for (ii = 0; ii < ped_linect; ii++) {
	    cptr = pheno_c;
	    cptr2 = &(pheno_c[ii]);
            cptr3 = (char*)ped_geno;
	    if ((!ped_geno[ii]) && (*cptr2 == 1)) {
	      kk++;
	      while (cptr < cptr2) {
                if (!(*cptr3)) {
		  if (*cptr == 1) {
		    high_tot_i += *iptr;
		  } else if (!(*cptr)) {
		    lh_tot_i += *iptr;
		  }
                }
		cptr++;
                cptr3++;
		iptr++;
	      }
	    } else if ((!ped_geno[ii]) && (*cptr2 == 0)) {
	      mm++;
	      while (cptr < cptr2) {
                if (!(*cptr3)) {
		  if (*cptr == 1) {
		    lh_tot_i += *iptr;
		  } else if (!(*cptr)) {
		    low_tot_i += *iptr;
		  }
                }
		cptr++;
                cptr3++;
		iptr++;
	      }
	    } else {
	      iptr += ii;
	    }
	  }
          low_tot = (double)low_tot_i;
          lh_tot = (double)lh_tot_i;
          high_tot = (double)high_tot_i;
	} else {
	  low_tot = 0.0;
	  lh_tot = 0.0;
	  high_tot = 0.0;
	  dist_ptr = dists;
	  for (ii = 0; ii < ped_linect; ii++) {
	    cptr = pheno_c;
	    cptr2 = &(pheno_c[ii]);
            cptr3 = (char*)ped_geno;
	    if ((!ped_geno[ii]) && (*cptr2 == 1)) {
	      kk++;
	      while (cptr < cptr2) {
                if (!(*cptr3)) {
		  if (*cptr == 1) {
		    high_tot += *dist_ptr;
		  } else if (!(*cptr)) {
		    lh_tot += *dist_ptr;
		  }
                }
		cptr++;
                cptr3++;
		dist_ptr++;
	      }
	    } else if ((!ped_geno[ii]) && (*cptr2 == 0)) {
	      mm++;
	      while (cptr < cptr2) {
                if (!(*cptr3)) {
		  if (*cptr == 1) {
		    lh_tot += *dist_ptr;
		  } else if (!(*cptr)) {
		    low_tot += *dist_ptr;
		  }
                }
		cptr++;
                cptr3++;
		dist_ptr++;
	      }
	    } else {
	      dist_ptr += ii;
	    }
	  }
	}
	if ((kk > 1) && (mm > 1)) {
	  high_tot = high_tot / (double)((kk * (kk - 1)) / 2);
	  lh_tot = lh_tot / (double)(kk * mm);
	  low_tot = low_tot / (double)((mm * (mm - 1)) / 2);
	  if (nn) {
	    dxx = dhh_sum / (double)nn;
	    dyy = dhl_sum / (double)nn;
	    dzz = dll_sum / (double)nn;
	  }
	  dhh_sum += high_tot;
	  dhl_sum += lh_tot;
	  dll_sum += low_tot;
	  nn++; // otherwise iteration doesn't count
	  if (nn > 1) {
	    dhh_ssd += (high_tot - (dhh_sum / (double)nn)) * (high_tot - dxx);
	    dhl_ssd += (lh_tot - (dhl_sum / (double)nn)) * (lh_tot - dyy);
	    dll_ssd += (low_tot - (dll_sum / (double)nn)) * (low_tot - dzz);
	  }
	}
      }
      if (nn < 2) {
	printf("Too few valid jackknife runs.\n");
      } else {
	printf("Jackknife results (%d valid runs), NO INFLATION FACTOR APPLIED:\n", nn);
	printf("  Avg dist between 2x affected             : %lf (sd %lf)\n", dhh_sum / (double)nn, sqrt(dhh_ssd / (double)(nn - 1)));
	printf("  Avg dist between affected and unaffected : %lf (sd %lf)\n", dhl_sum / (double)nn, sqrt(dhl_ssd / (double)(nn - 1)));
	printf("  Avg dist between 2x unaffected           : %lf (sd %lf)\n", dll_sum / (double)nn, sqrt(dll_ssd / (double)(nn - 1)));
      }
    }
  }

  /*
  for (window_num = 0; window_num < windows_reqd; window_num += 1) {
    if (window_num == windows_reqd - 1) {
      cur_window_size = map_linect - window_num * geno_window_size * 4;
      cur_window_size4 = (cur_window_size + 3) / 4;
    } else {
      cur_window_size = geno_window_size * 4;
      cur_window_size4 = geno_window_size;
    }
    memset(ped_geno, 0, ped_linect * cur_window_size4 * sizeof(char));
    memset(marker_alleles, 0, cur_window_size4 * 8 * sizeof(char));
    memset(marker_allele_cts, 0, cur_window_size4 * 8 * sizeof(int));
    for (ii = 0; ii < ped_linect; ii += 1) {
      fseeko(pedfile, line_locs[ii], SEEK_SET);
      fgets(pedbuf, smallbuflen, pedfile);
      bufptr = pedbuf;
      if (!window_num) {
	if (affection) {
	  if (affection_01) {
            if (*bufptr == '0') {
              pheno_c[ii] = 0;
            } else {
              pheno_c[ii] = 1;
            }
	  } else {
            if (*bufptr == '1') {
              pheno_c[ii] = 0;
            } else {
              pheno_c[ii] = 1;
            }
	  }
        } else {
	  if (sscanf(bufptr, "%lf", &(pheno_d[ii])) != 1) {
	    retval = RET_INVALID_FORMAT;
	    printf(errstr_ped_format);
	    goto wdist_ret_9;
	  }
        }
        bufptr = next_item(bufptr);
        if (ped_col_7) {
          bufptr = next_item(bufptr);
        }
      }
      for (jj = 0; jj < cur_window_size; jj += 1) {
	mm = jj % 4;
	nn = 1 << (mm * 2);
        if (!mm) {
          gptr = &(ped_geno[ped_linect * (jj / 4) + ii]);
        }
        for (kk = 0; kk < 2; kk++) {
          cc = *bufptr;
          if (cc == '0') {
            retval = RET_CALC_NOT_YET_SUPPORTED;
            printf("Error: No-calls not yet supported in distance calculation.\n");
            goto wdist_ret_9;
          }
          if (cc == marker_alleles[jj * 2]) {
	    marker_allele_cts[jj * 2] += 1;
	  } else if (cc == marker_alleles[jj * 2 + 1]) {
	    marker_allele_cts[jj * 2 + 1] += 1;
	    *gptr += nn;
	  } else if (marker_alleles[jj * 2]) {
	    if (marker_alleles[jj * 2 + 1]) {
	      retval = RET_INVALID_FORMAT;
	      printf("Error: More than two different allele types at marker %d.\n", jj + 1);
	      goto wdist_ret_9;
	    } else {
	      marker_alleles[jj * 2 + 1] = cc;
	      marker_allele_cts[jj * 2 + 1] += 1;
	      *gptr += nn;
	    }
	  } else {
	    marker_alleles[jj * 2] = cc;
	    marker_allele_cts[jj * 2 + 1] += 1;
	  }
          bufptr++;
          while ((*bufptr == ' ') || (*bufptr == '\t')) {
            bufptr++;
          }
          if (*bufptr == '\0') {
	    retval = RET_INVALID_FORMAT;
	    printf(errstr_ped_format);
	    goto wdist_ret_9;
          }
        }
      }
      line_locs[ii] += (bufptr - pedbuf);
    }
    if (exponent == 0.0) {
      for (ii = 0; ii < cur_window_size4; ii += 1) {
        gptr = &(ped_geno[ped_linect * ii + 1]);
        iwptr = idists;
        for (jj = 1; jj < ped_linect; jj += 1) {
          gptr2 = &(ped_geno[ped_linect * ii]);
          for (kk = 0; kk < jj; kk += 1) {
            *iwptr += iwt[(unsigned char)((*gptr2++) ^ (*gptr))];
            iwptr++;
          }
          gptr++;
        }
      }
    } else {
      for (ii = 0; ii < cur_window_size4; ii += 1) {
        gptr = &(ped_geno[ped_linect * ii + 1]);
        wptr = dists;
        fill_weights(weights, marker_allele_cts, exponent, min_maf);
        for (jj = 1; jj < ped_linect; jj += 1) {
          gptr2 = &(ped_geno[ped_linect * ii]);
          for (kk = 0; kk < jj; kk += 1) {
            *wptr += weights[(unsigned char)((*gptr2++) ^ (*gptr))];
            wptr++;
          }
          gptr++;
        }
      }
    }
  }
  */

 wdist_ret_2:
  if (person_missing_cts) {
    free(person_missing_cts);
  }
  if (marker_alleles_tmp) {
    free(marker_alleles_tmp);
  }
  if (bedtmpfile) {
    fclose(bedtmpfile);
  }
  if (bimtmpfile) {
    fclose(bimtmpfile);
  }
  if (famtmpfile) {
    fclose(famtmpfile);
  }
 wdist_ret_1:
  if (person_id) {
    free(person_id);
  }
  if (person_exclude) {
    free(person_exclude);
  }
  if (pheno_d) {
    free(pheno_d);
  }
  if (pheno_c) {
    free(pheno_c);
  }
  if (phenor_d) {
    free(phenor_d);
  }
  if (phenor_c) {
    free(phenor_c);
  }
  if (marker_allele_cts) {
    free(marker_allele_cts);
  }
  if (marker_alleles) {
    free(marker_alleles);
  }
  if (line_locs) {
    free(line_locs);
  }
  if (id_buf) {
    free(id_buf);
  }
  if (id_list) {
    free(id_list);
  }
  if (pid_list) {
    free(pid_list);
  }
  if (marker_exclude) {
    free(marker_exclude);
  }
  if (pedbuf) {
    free(pedbuf);
  }
  // if (marker_pos) {
  //   free(marker_pos);
  // }
  // if (marker_id) {
  //   free(marker_id);
  // }
  // if (marker_chrom) {
  //   free(marker_chrom);
  // }
 wdist_ret_0:
  if (outfile) {
    fclose(outfile);
  }
  if (filterfile) {
    fclose(filterfile);
  }
  if (phenofile) {
    fclose(phenofile);
  }
  if (famfile) {
    fclose(famfile);
  }
  if (mapfile) {
    fclose(mapfile);
  }
  if (pedfile) {
    fclose(pedfile);
  }
  return retval;
}

int main(int argc, char** argv) {
  char mapname[FNAMESIZE];
  char pedname[FNAMESIZE];
  char famname[FNAMESIZE];
  char outname[FNAMESIZE];
  char phenoname[FNAMESIZE];
  char filtername[FNAMESIZE];
  char* makepheno_str = NULL;
  char* filterval;
  char* argptr;
  int retval;
  int load_params = 0; // describes what file parameters have been provided
  int make_bed = 1; // STOPGAP
  int ped_col_1 = 1;
  int ped_col_34 = 1;
  int ped_col_5 = 1;
  int ped_col_6 = 1;
  int ped_col_7 = 0;
  int mpheno_col = 0;
  char* phenoname_str = NULL;
  int affection_01 = 0;
  double exponent = 0.0;
  int chr_num = 0;
  double min_maf = 0.000000002;
  double geno_thresh = 1.0;
  double mind_thresh = 1.0;
  double hwe_thresh = 0.0;
  int cur_arg = 1;
  int calculation_type = CALC_NONE;
  char* bubble;
  int thread_ct = 2;
  int filter_type = 0;
  int mfilter_col = 0;
  int tail_pheno = 0;
  int prune = 0;
  int missing_pheno = -9;
  unsigned char missing_geno = '0';
  double tail_bottom;
  double tail_top;
  int calc_param_1 = 0;
  int iters;
  int ii;
  int jj;
  unsigned long int rseed = 0;
  FILE* scriptfile;
  int num_params;
  int in_param;
  strcpy(mapname, "wdist.map");
  strcpy(pedname, "wdist.ped");
  famname[0] = '\0';
  phenoname[0] = '\0';
  filtername[0] = '\0';
  strcpy(outname, "wdist");
  while (cur_arg < argc) {
    argptr = argv[cur_arg];
    if (!strcmp(argptr, "--script")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --script parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (subst_argv) {
        printf("Error: Duplicate --script flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      scriptfile = fopen(argv[cur_arg + 1], "rb");
      if (!scriptfile) {
        printf("Error: Failed to open %s.\n", argv[cur_arg + 1]);
        return dispmsg(RET_OPENFAIL);
      }
      ii = fread(tbuf, 1, MAXLINELEN, scriptfile);
      fclose(scriptfile);
      if (ii == MAXLINELEN) {
        printf("Error: Script file too long (max %d bytes).\n", MAXLINELEN - 1);
        return dispmsg(RET_INVALID_FORMAT);
      }
      num_params = 0;
      in_param = 0;
      for (jj = 0; jj < ii; jj++) {
        if (is_space_or_eoln(tbuf[jj])) {
          if (in_param) {
            in_param = 0;
          }
        } else if (!in_param) {
          num_params += 1;
          in_param = 1;
        }
      }
      subst_argv = (char**)malloc((num_params + argc - cur_arg - 2) * sizeof(char*));
      num_params = 0;
      in_param = 0;
      for (jj = 0; jj < ii; jj++) {
        if (is_space_or_eoln(tbuf[jj])) {
          if (in_param) {
            tbuf[jj] = '\0';
            in_param = 0;
          }
        } else if (!in_param) {
          subst_argv[num_params++] = &(tbuf[jj]);
          in_param = 1;
        }
      }
      for (ii = cur_arg + 2; ii < argc; ii++) {
        subst_argv[num_params++] = argv[ii];
      }
      argc = num_params;
      cur_arg = 0;
      argv = subst_argv;
    } else if (!strcmp(argptr, "--file")) {
      if (load_params & 121) {
        if (load_params & 1) {
          printf("Error: Duplicate --file flag.\n");
        } else {
          printf("Error: --file flag cannot coexist with binary file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 1;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --file parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 5)) {
        printf("Error: --file parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      if (!(load_params & 2)) {
	strcpy(pedname, argv[cur_arg + 1]);
	strcat(pedname, ".ped");
      }
      if (!(load_params & 4)) {
	strcpy(mapname, argv[cur_arg + 1]);
	strcat(mapname, ".map");
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--help")) {
      return dispmsg(RET_MOREHELP);
    } else if (!strcmp(argptr, "--ped")) {
      if (load_params & 122) {
        if (load_params & 2) {
          printf("Error: Duplicate --ped flag.\n");
        } else {
          printf("Error: --ped flag cannot coexist with binary file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 2;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --ped parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --ped parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(pedname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--map")) {
      if (load_params & 124) {
        if (load_params & 4) {
          printf("Error: Duplicate --map flag.\n");
        } else {
          printf("Error: --map flag cannot coexist with binary file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 4;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --map parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --map parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(mapname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--make-bed")) {
      if (load_params & 120) {
        printf("Error: --make-bed cannot coexist with binary file flags.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      make_bed = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--bfile")) {
      if (load_params & 15) {
        if (load_params & 8) {
          printf("Error: Duplicate --bfile flag.\n");
        } else {
          printf("Error: --bfile flag cannot coexist with text input file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 8;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --bfile parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 5)) {
        printf("Error: --bfile parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      make_bed = 0;
      if (!(load_params & 16)) {
	strcpy(pedname, argv[cur_arg + 1]);
	strcat(pedname, ".bed");
      }
      if (!(load_params & 32)) {
	strcpy(mapname, argv[cur_arg + 1]);
	strcat(mapname, ".bim");
      }
      if (!(load_params & 64)) {
        strcpy(famname, argv[cur_arg + 1]);
        strcat(famname, ".fam");
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--bed")) {
      if (load_params & 23) {
        if (load_params & 16) {
          printf("Error: Duplicate --bed flag.\n");
        } else {
          printf("Error: --bed flag cannot coexist with text input file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 16;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --bed parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --bed parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      make_bed = 0;
      strcpy(pedname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--bim")) {
      if (load_params & 39) {
        if (load_params & 32) {
          printf("Error: Duplicate --bim flag.\n");
        } else {
          printf("Error: --bim flag cannot coexist with text input file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 32;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --bim parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --bim parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      make_bed = 0;
      strcpy(mapname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--fam")) {
      if (load_params & 71) {
        if (load_params & 64) {
          printf("Error: Duplicate --fam flag.\n");
        } else {
          printf("Error: --fam flag cannot coexist with text input file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 64;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --fam parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --fam parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      make_bed = 0;
      strcpy(famname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--no-fid")) {
      ped_col_1 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--no-parents")) {
      ped_col_34 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--no-sex")) {
      ped_col_5 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--no-pheno")) {
      ped_col_6 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--liability")) {
      ped_col_7 = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--missing-genotype")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --missing-genotype parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (missing_geno != '0') {
        printf("Error: Duplicate --missing-genotype flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      missing_geno = (unsigned char)argv[cur_arg + 1][0];
      if ((strlen(argv[cur_arg + 1]) > 1) || (missing_geno <= ' ')) {
        printf("Error: Invalid --missing-genotype parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--missing-phenotype")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --missing-phenotype parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (missing_pheno != -9) {
        printf("Error: Duplicate --missing-phenotype flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      missing_pheno = atoi(argv[cur_arg + 1]);
      if ((missing_pheno == 0) || (missing_pheno == 1)) {
        printf("Error: Invalid --missing-phenotype parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--pheno")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --pheno parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (phenoname[0]) {
        if (makepheno_str) {
          printf("Error: --pheno and --make-pheno flags cannot coexist.\n");
        } else {
          printf("Error: Duplicate --pheno flag.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --pheno parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(phenoname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--make-pheno")) {
      if (cur_arg > argc - 2) {
        printf("Error: Not enough --make-pheno parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (phenoname[0]) {
        if (makepheno_str) {
          printf("Error: Duplicate --make-pheno flag.\n");
        } else {
          printf("Error: --pheno and --make-pheno flags cannot coexist.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (tail_pheno) {
        printf("Error: --make-pheno and --tail-pheno flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --make-pheno parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(phenoname, argv[cur_arg + 1]);
      makepheno_str = argv[cur_arg + 2];
      cur_arg += 3;
    } else if (!strcmp(argptr, "--mpheno")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --mpheno parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (mpheno_col != 0) {
        printf("Error: Duplicate --mpheno flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (phenoname_str) {
        printf("Error: --mpheno and --pheno-name flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      mpheno_col = atoi(argv[cur_arg + 1]);
      if (mpheno_col < 1) {
        printf("Error: Invalid --mpheno parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--pheno-name")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --pheno-name parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (mpheno_col != 0) {
        printf("Error: --mpheno and --pheno-name flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (phenoname_str) {
        printf("Error: Duplicate --pheno-name flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      phenoname_str = argv[cur_arg + 1];
      cur_arg += 2;
    } else if (!strcmp(argptr, "--prune")) {
      prune = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--1")) {
      affection_01 = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--tail-pheno")) {
      if (cur_arg > argc - 2) {
        printf("Error: Not enough --tail-pheno parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (tail_pheno) {
        printf("Error: Duplicate --tail-pheno flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (makepheno_str) {
        printf("Error: --make-pheno and --tail-pheno flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 1], "%lf", &tail_bottom) != 1) {
        printf("Error: Invalid --tail-pheno parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 2], "%lf", &tail_top) != 1) {
        printf("Error: Invalid --tail-pheno parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      tail_pheno = 1;
      cur_arg += 3;
    } else if (!strcmp(argptr, "--keep")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --keep parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (filter_type) {
        if (filter_type == 1) {
          printf("Error: Duplicate --keep flag.\n");
        } else {
          printf("Error: --keep + --remove/--filter not supported.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      filter_type = FILTER_KEEP;
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--keep filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(filtername, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--remove")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --remove parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (filter_type) {
        if (filter_type == 2) {
          printf("Error: Duplicate --remove flag.\n");
        } else {
          printf("Error: --remove + --keep/--filter not supported.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      filter_type = FILTER_REMOVE;
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--remove filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(filtername, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--filter")) {
      if (cur_arg > argc - 2) {
        printf("Error: Not enough --filter parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (filter_type) {
        if (filter_type == 3) {
          printf("Error: Duplicate --filter flag.\n");
        } else {
          printf("Error: --filter + --keep/--remove not supported.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      filter_type = FILTER_CUSTOM;
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--filter filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(filtername, argv[cur_arg + 1]);
      filterval = argv[cur_arg + 2];
      cur_arg += 3;
    } else if (!strcmp(argptr, "--mfilter")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --mfilter parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (mfilter_col) {
        printf("Error: Duplicate --mfilter flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      mfilter_col = atoi(argv[cur_arg + 1]);
      if (mfilter_col < 1) {
        printf("Error: Invalid --mfilter parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--memory")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --memory parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      malloc_size_mb = atoi(argv[cur_arg + 1]);
      if (malloc_size_mb < 64) {
        printf("Error: Invalid --memory parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--threads")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --threads parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) {
        printf("Error: Invalid --threads parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > MAX_THREADS) {
        printf("Error: --threads parameter too large (max %d).\n", MAX_THREADS);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      thread_ct = ii;
      cur_arg += 2;
    } else if (!strcmp(argptr, "--exponent")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --exponent parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 1], "%lf", &exponent) != 1) {
        printf("Error: Invalid --exponent parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--chr")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --chr parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (chr_num) {
        printf("Error: Duplicate --chr flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      chr_num = marker_code(argv[cur_arg + 1]);
      if ((chr_num < 1) || (chr_num > 26)) {
        printf("Error: Invalid --chr parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--maf")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --maf parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 1], "%lf", &min_maf) != 1) {
        printf("Error: Invalid --maf parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if ((min_maf <= 0.0) || (min_maf > 0.5)) {
        printf("Error: Invalid --maf parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--geno")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --geno parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 1], "%lf", &geno_thresh) != 1) {
        printf("Error: Invalid --geno parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if ((geno_thresh < 0.0) || (geno_thresh > 1.0)) {
        printf("Error: Invalid --geno parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--mind")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --mind parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 1], "%lf", &mind_thresh) != 1) {
        printf("Error: Invalid --mind parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if ((mind_thresh < 0.0) || (mind_thresh > 1.0)) {
        printf("Error: Invalid --mind parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--hwe")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --hwe parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 1], "%lf", &hwe_thresh) != 1) {
        printf("Error: Invalid --hwe parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if ((hwe_thresh < 0.0) || (hwe_thresh >= 1.0)) {
        printf("Error: Invalid --hwe parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--rseed")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --rseed parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (rseed != 0) {
        printf("Error: Duplicate --rseed flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      rseed = (unsigned long int)atoi(argv[cur_arg + 1]);
      if (rseed == 0) {
        printf("Error: Invalid --rseed parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--out")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --out parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 9)) {
        printf("Error: --out parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(outname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--silent")) {
      freopen("/dev/null", "w", stdout);
      cur_arg += 1;
    } else if (!strcmp(argptr, "--make-rel")) {
      if (cur_arg != argc - 1) {
        if ((cur_arg == argc - 2) && (!strcmp(argv[cur_arg + 1], "--square0"))) {
          calc_param_1 = 1;
        } else {
          printf("Error: Invalid parameter after --relationship.%s", errstr_append);
          return dispmsg(RET_INVALID_CMDLINE);
        }
      }
      cur_arg = argc;
      calculation_type = CALC_RELATIONSHIP;
    } else if (!strcmp(argptr, "--make-grm")) {
      if (cur_arg != argc - 1) {
        printf("Error: Invalid parameter after --relationship-grm.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      calc_param_1 = 2;
      cur_arg = argc;
      calculation_type = CALC_RELATIONSHIP;
    } else if (!strcmp(argptr, "--distance")) {
      if (cur_arg != argc - 1) {
        if ((cur_arg == argc - 2) && (!strcmp(argv[cur_arg + 1], "--square0"))) {
          calc_param_1 = 1;
        } else {
	  printf("Error: Invalid parameter after --distance.%s", errstr_append);
	  return dispmsg(RET_INVALID_CMDLINE);
        }
      }
      cur_arg = argc;
      calculation_type = CALC_DISTANCE;
    } else if (!strcmp(argptr, "--groupdist")) {
      if (cur_arg != argc - 3) {
        printf("Error: --groupdist requires 2 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      calc_param_1 = atoi(argv[cur_arg + 1]);
      if (calc_param_1 <= 0) {
        printf("Error: Invalid --groupdist jackknife delete parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      iters = atoi(argv[cur_arg + 2]);
      if (iters < 2) {
        printf("Error: Invalid --groupdist jackknife iteration count.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg = argc;
      calculation_type = CALC_GROUPDIST;
    } else if (!strcmp(argptr, "--map3")) {
      printf("Note: --map3 flag unnecessary (.map file format is autodetected).\n");
      cur_arg += 1;
    } else if (!strcmp(argptr, "--compound-genotypes")) {
      printf("Note: --compound-genotypes flag unnecessary (spaces between alleles in .ped\nare optional).\n");
      cur_arg += 1;
    } else {
      printf("Error: Invalid argument (%s).%s", argv[cur_arg], errstr_append);
      return dispmsg(RET_INVALID_CMDLINE);
    }
  }
  if (calculation_type == CALC_NONE) {
    return dispmsg(RET_HELP);
  }
  if (prune && (!phenoname[0]) && (!ped_col_6)) {
    printf("Error: --prune and --no-pheno cannot coexist without an alternate phenotype\nfile.\n");
    return dispmsg(RET_INVALID_CMDLINE);
  }
  if (subst_argv) {
    free(subst_argv);
  }
  if (exponent == 0.0) {
    for (ii = 0; ii < 256; ii += 1) {
      iwt[ii] = 0;
      jj = ii;
      while (jj) {
        if (jj % 2 == 1) {
          iwt[ii] += 1;
        }
        jj >>= 1;
      }
    }
  }
  bubble = (char*)malloc(67108864 * sizeof(char));
  if (!bubble) {
    return dispmsg(RET_NOMEM);
  }
  wkspace = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  if ((malloc_size_mb > MALLOC_DEFAULT_MB) && !wkspace) {
    printf("%lld MB malloc failed.  Using default allocation behavior.\n", malloc_size_mb);
    malloc_size_mb = MALLOC_DEFAULT_MB;
  }
  while (!wkspace) {
    if (malloc_size_mb > 128) {
      malloc_size_mb -= 64;
    } else {
      malloc_size_mb = 64;
    }
    wkspace = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
    if (wkspace) {
      printf("Allocated %lld MB successfully.\n", malloc_size_mb);
    }
  }
  free(bubble);

  // rg = gsl_rng_alloc(gsl_rng_mt19937);
  if (!rseed) {
    rseed = (unsigned long int)time(NULL);
  }
  srand(rseed);
  // gsl_rng_set(rg, rseed);

  // famname[0] indicates binary vs. text
  // filtername[0] indicates existence of filter
  retval = wdist(pedname, mapname, famname, phenoname, filtername, makepheno_str, filter_type, filterval, mfilter_col, make_bed, ped_col_1, ped_col_34, ped_col_5, ped_col_6, ped_col_7, (char)missing_geno, missing_pheno, mpheno_col, phenoname_str, prune, affection_01, chr_num, thread_ct, exponent, min_maf, geno_thresh, mind_thresh, hwe_thresh, tail_pheno, tail_bottom, tail_top, outname, calculation_type, calc_param_1, iters);
  // gsl_rng_free(rg);
  free(wkspace);
  return dispmsg(retval);
}
