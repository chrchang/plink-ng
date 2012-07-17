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

#define FILTER_KEEP 1
#define FILTER_REMOVE 2
#define FILTER_CUSTOM 3

#define CALC_NONE -1
#define CALC_DISTANCE 0

#define _FILE_OFFSET_BITS 64
#define MAX_THREADS 64
#define PEDBUFBASE 256
#define FNAMESIZE 2048
#define MALLOC_DEFAULT_MB 2176
#define IDLENGTH 16
// size of generic text line load buffer.  .ped lines can of course be longer
#define MAXLINELEN 131072

const char info_str[] =
  "WDIST weighted genetic distance calculator, v0.2.0 (17 July 2012)\n"
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
"  --file [prefix]  : Specify prefix for .ped and .map files (default 'wdist').\n"
"  --ped [filename] : Specify name of .ped file.\n"
"  --map [filename] : Specify name of .map file.\n"
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
"  --maf [val]      : Minor allele frequency minimum threshold (default 0.01).\n"
"  --geno [val]     : Maximum per-SNP missing (default 0.1).\n"
"  --mind [val]     : Maximum per-person missing (default 0.1).\n"
"  --rseed [val]    : Set random number seed (relevant for missing genotypes).\n"
"  --memory [val]   : Size, in MB, of initial malloc attempt (default 2176).\n"
"  --threads [val]  : Maximum number of concurrent threads (default 2).\n"
"  --exponent [val] : When computing genetic distances, each locus has a weight\n"
"                     of (2q(1-q))^{-val}, where q is the observed MAF.\n\n"
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
"  --distance\n"
"    Outputs a lower-triangular table of (weighted) genetic distances.\n"
"    The first row contains a single number with the distance between the first\n"
"    two genotypes, the second row has the {genotype 1-genotype 3} and\n"
"    {genotype 2-genotype 3} distances in that order, etc.\n\n"
"  --groupdist [d] [iters]\n"
"    Two-group genetic distance analysis, using delete-d jackknife with the\n"
"    requested number of iterations.  Binary phenotype required.  (Not actually\n"
"    implemented yet, but soon will be.)\n\n",
	   info_str);
    break;
  }
  return retval;
}

inline int is_space_or_eoln(char cc) {
  return ((cc == ' ') || (cc == '\t') || (cc == '\n') || (cc == '\0'));
}

char* id_buf = NULL;

int qsort_partition(char* lptr, int max_id_len, int min_idx, int max_idx, int pivot_idx) {
  char* pivot_ptr = &(lptr[pivot_idx * max_id_len]);
  char* right_ptr = &(lptr[max_idx * max_id_len]);
  char* store_ptr = &(lptr[min_idx * max_id_len]);
  char* incr_ptr;
  int si = min_idx;
  strcpy(id_buf, pivot_ptr);
  strcpy(pivot_ptr, right_ptr);
  strcpy(right_ptr, id_buf);
  while (store_ptr < right_ptr) {
    if (strcmp(store_ptr, right_ptr) < 0) {
      store_ptr += max_id_len;
      si++;
    } else {
      for (incr_ptr = store_ptr + max_id_len; incr_ptr < pivot_ptr; incr_ptr += max_id_len) {
	if (strcmp(incr_ptr, right_ptr) < 0) {
	  strcpy(id_buf, incr_ptr);
	  strcpy(incr_ptr, store_ptr);
	  strcpy(store_ptr, id_buf);
	  store_ptr += max_id_len;
	  si++;
	}
      }
      break;
    }
  }
  strcpy(id_buf, store_ptr);
  strcpy(store_ptr, right_ptr);
  strcpy(right_ptr, id_buf);  
  return si;
}

void qsort_str(char* lptr, int max_id_len, int min_idx, int max_idx) {
  int pivot_idx;
  if (max_idx > min_idx) {
    pivot_idx = qsort_partition(lptr, max_id_len, min_idx, max_idx, (min_idx + max_idx) / 2);
    qsort_str(lptr, max_id_len, min_idx, pivot_idx - 1);
    qsort_str(lptr, max_id_len, pivot_idx + 1, max_idx);
  }
}

int qsort_partition_c(char* lptr, char* cptr, int max_id_len, int min_idx, int max_idx, int pivot_idx) {
  char* pivot_ptr = &(lptr[pivot_idx * max_id_len]);
  char* right_ptr = &(lptr[max_idx * max_id_len]);
  char* store_ptr = &(lptr[min_idx * max_id_len]);
  char* pivot_ptr_c = &(cptr[pivot_idx]);
  char* right_ptr_c = &(cptr[max_idx]);
  char* store_ptr_c = &(cptr[min_idx]);
  char cc;
  char* incr_ptr;
  char* incr_ptr_c;
  int si = min_idx;
  strcpy(id_buf, pivot_ptr);
  cc = *pivot_ptr_c;
  strcpy(pivot_ptr, right_ptr);
  *pivot_ptr_c = *right_ptr_c;
  strcpy(right_ptr, id_buf);
  *right_ptr_c = cc;
  while (store_ptr < right_ptr) {
    if (strcmp(store_ptr, right_ptr) < 0) {
      store_ptr += max_id_len;
      si++;
    } else {
      incr_ptr = store_ptr + max_id_len;
      incr_ptr_c = store_ptr_c + 1;
      while (incr_ptr < pivot_ptr) {
	if (strcmp(incr_ptr, right_ptr) < 0) {
	  strcpy(id_buf, incr_ptr);
          cc = *incr_ptr_c;
	  strcpy(incr_ptr, store_ptr);
          *incr_ptr_c = *store_ptr_c;
	  strcpy(store_ptr, id_buf);
          *store_ptr_c++ = cc;
	  store_ptr += max_id_len;
	  si++;
	}
        incr_ptr += max_id_len;
        incr_ptr_c++;
      }
      break;
    }
  }
  strcpy(id_buf, store_ptr);
  cc = *store_ptr_c;
  strcpy(store_ptr, right_ptr);
  *store_ptr_c = *right_ptr_c;
  strcpy(right_ptr, id_buf);  
  *right_ptr_c = cc;
  return si;
}

void qsort_str_c(char* lptr, char* cptr, int max_id_len, int min_idx, int max_idx) {
  int pivot_idx;
  if (max_idx > min_idx) {
    pivot_idx = qsort_partition_c(lptr, cptr, max_id_len, min_idx, max_idx, (min_idx + max_idx) / 2);
    qsort_str_c(lptr, cptr, max_id_len, min_idx, pivot_idx - 1);
    qsort_str_c(lptr, cptr, max_id_len, pivot_idx + 1, max_idx);
  }
}

int qsort_partition_d(char* lptr, double* dptr, int max_id_len, int min_idx, int max_idx, int pivot_idx) {
  char* pivot_ptr = &(lptr[pivot_idx * max_id_len]);
  char* right_ptr = &(lptr[max_idx * max_id_len]);
  char* store_ptr = &(lptr[min_idx * max_id_len]);
  double* pivot_ptr_d = &(dptr[pivot_idx]);
  double* right_ptr_d = &(dptr[max_idx]);
  double* store_ptr_d = &(dptr[min_idx]);
  double dxx;
  char* incr_ptr;
  double* incr_ptr_d;
  int si = min_idx;
  strcpy(id_buf, pivot_ptr);
  dxx = *pivot_ptr_d;
  strcpy(pivot_ptr, right_ptr);
  *pivot_ptr_d = *right_ptr_d;
  strcpy(right_ptr, id_buf);
  *right_ptr_d = dxx;
  while (store_ptr < right_ptr) {
    if (strcmp(store_ptr, right_ptr) < 0) {
      store_ptr += max_id_len;
      si++;
    } else {
      incr_ptr = store_ptr + max_id_len;
      incr_ptr_d = store_ptr_d + 1;
      while (incr_ptr < pivot_ptr) {
	if (strcmp(incr_ptr, right_ptr) < 0) {
	  strcpy(id_buf, incr_ptr);
          dxx = *incr_ptr_d;
	  strcpy(incr_ptr, store_ptr);
          *incr_ptr_d = *store_ptr_d;
	  strcpy(store_ptr, id_buf);
          *store_ptr_d++ = dxx;
	  store_ptr += max_id_len;
	  si++;
	}
        incr_ptr += max_id_len;
        incr_ptr_d++;
      }
      break;
    }
  }
  strcpy(id_buf, store_ptr);
  dxx = *store_ptr_d;
  strcpy(store_ptr, right_ptr);
  *store_ptr_d = *right_ptr_d;
  strcpy(right_ptr, id_buf);  
  *right_ptr_d = dxx;
  return si;
}

void qsort_str_d(char* lptr, double* dptr, int max_id_len, int min_idx, int max_idx) {
  int pivot_idx;
  if (max_idx > min_idx) {
    pivot_idx = qsort_partition_d(lptr, dptr, max_id_len, min_idx, max_idx, (min_idx + max_idx) / 2);
    qsort_str_d(lptr, dptr, max_id_len, min_idx, pivot_idx - 1);
    qsort_str_d(lptr, dptr, max_id_len, pivot_idx + 1, max_idx);
  }
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
  } else {
    return (sptr[0] - 'A'); // X = 23, Y = 24
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
  double wt[4];
  double twt[3];
  for (ii = 0; ii < 4; ii += 1) {
    wt[ii] = pow(2.0 * mafs[ii] * (1.0 - mafs[ii]), -exponent);
  }
  nn = 0;
  for (ii = 0; ii < 4; ii += 1) {
    if ((ii % 4 == 1) || (ii % 4 == 2)) {
      twt[0] = wt[3];
    } else if (ii % 4 == 3) {
      twt[0] = 0;
    } else {
      twt[0] = wt[3] * 2;
    }
    for (jj = 0; jj < 4; jj += 1) {
      if ((jj % 4 == 1) || (jj % 4 == 2)) {
        twt[1] = twt[0] + wt[2];
      } else if (jj % 4 == 3) {
        twt[1] = twt[0];
      } else {
        twt[1] = twt[0] + 2 * wt[2];
      }
      for (kk = 0; kk < 4; kk += 1) {
        if ((kk % 4 == 1) || (kk % 4 == 2)) {
          twt[2] = twt[1] + wt[1];
	} else if (kk % 4 == 3) {
          twt[2] = twt[1];
        } else {
          twt[2] = twt[1] + 2 * wt[1];
        }
        for (mm = 0; mm < 4; mm += 1) {
          if ((mm % 4 == 1) || (mm % 4 == 2)) {
            *weights++ = twt[2] + wt[0];
          } else if (mm % 4 == 3) {
            *weights++ = twt[2];
          } else {
            *weights++ = twt[2] + 2 * wt[0];
	  }
        }
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

int wdist(char* pedname, char* mapname, char* famname, char* phenoname, char* filtername, char* makepheno_str, int filter_type, char* filterval, int mfilter_col, int ped_col_1, int ped_col_34, int ped_col_5, int ped_col_6, int ped_col_7, char missing_geno, int missing_pheno, int mpheno_col, char* phenoname_str, int prune, int affection_01, int threads, double exponent, double min_maf, double geno_thresh, double mind_thresh, int tail_pheno, double tail_bottom, double tail_top, char* outname, int calculation_type) {
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
  int ped_linect4;
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int oo;
  int pp;
  double dxx;
  long long llxx;
  long long llyy;
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
  unsigned char* ped_geno = NULL;
  unsigned char* gptr;
  char* cptr;
  int* iwptr;
  int* iptr;
  int* idists;
  double* dists;
  double weights[256];
  long long dists_alloc = 0;
  long long ped_geno_size;
  int geno_window_size;
  char cc;
  double* dist_ptr;
  int last_tell;
  int binary_files = 0;
  int maf_int_thresh;
  int geno_int_thresh;
  int mind_int_thresh;
  int marker_exclude_ct = 0;
  int person_exclude_ct = 0;
  int pheno_lines = 0;
  int makepheno_all = 0;
  int filter_lines = 0;
  int snp_major = 0;
  int max_pid_len = 0;
  int max_id_len = 4;
  char* pid_list = NULL;
  char* id_list = NULL;
  double maf_buf[4];
  double missing_phenod = (double)missing_pheno;
  int missing_pheno_len = 1;

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
  strcpy(outname_end, ".dist");
  outfile = fopen(outname, "w");
  if (!outfile) {
    printf("Error: Failed to open %s.\n", outname);
    goto wdist_ret_0;
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
    // marker_pos[ii] = atoi(bufptr);
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
	qsort_str_c(pid_list, phenor_c, max_pid_len, 0, pheno_lines - 1);
      } else if (makepheno_str) {
	qsort_str(pid_list, max_pid_len, 0, pheno_lines - 1);
      } else {
	qsort_str_d(pid_list, phenor_d, max_pid_len, 0, pheno_lines - 1);
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
      qsort_str(id_list, max_id_len, 0, filter_lines - 1);
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
  maf_int_thresh = 2 * map_linect - (int)((1.0 - min_maf) * 2.0 * map_linect);
  geno_int_thresh = 2 * map_linect - (int)(geno_thresh * 2.0 * map_linect);
  mind_int_thresh = (int)(mind_thresh * map_linect);

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

    if (makepheno_str || affection || tail_pheno) {
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
          if (phenoname[0]) {
	    if (!ped_col_1) {
	      cptr = tbuf;
	    } else {
	      cptr = bufptr;
	    }
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
        fseeko(pedfile, 3 + line_locs[ii] * map_linect4, SEEK_SET);
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
            if (oo == 1) {
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
        if ((marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1] < geno_int_thresh) || (marker_allele_cts[ii * 2 + 1] < maf_int_thresh)) {
          exclude(marker_exclude, ii, &marker_exclude_ct);
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
            if (pp == 1) {
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
        if ((mm + nn < geno_int_thresh) || (nn < maf_int_thresh)) {
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
    }
  } else {
    pedbuf[pedbuflen - 1] = ' ';
    last_tell = 0;
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
    if (affection) {
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

    ii = ((ped_linect - person_exclude_ct) * ((ped_linect - person_exclude_ct) - 1)) / 2;
    if (exponent == 0.0) {
      dists_alloc = ii * sizeof(int);
    } else {
      dists_alloc = ii * sizeof(double);
    }
    
    geno_window_size = (malloc_size_mb * 1048576 - dists_alloc) / (ped_linect - person_exclude_ct);
    map_linect4 = (map_linect - marker_exclude_ct + 3) / 4;
    if (map_linect4 > geno_window_size) {
      printf(".ped file too large for direct read.  Converting to binary...");
      strcpy(outname_end, ".bed.tmp");
      bedtmpfile = fopen(outname, "wb");
      if (!bedtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_1;
      }
      strcpy(outname_end, ".bim.tmp");
      bimtmpfile = fopen(outname, "wb");
      if (!bimtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      strcpy(outname_end, ".fam.tmp");
      famtmpfile = fopen(outname, "wb");
      if (!famtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (3 != fwrite("L\x1b", 1, 3, bedtmpfile)) {
	retval = RET_WRITE_FAIL;
        printf("\n");
	goto wdist_ret_2;
      }
      // ----- .ped (third pass) -> .fam + .bed conversion -----
      rewind(pedfile);
      ii = 0; // line count
      llxx = 0; // last file location
      while (fgets((char*)pedbuf, pedbuflen, pedfile) != NULL) {
        if (ii == ped_linect) {
          break;
        }
        llyy = ftello(pedfile);
        if (llyy < line_locs[ii]) {
          llxx = llyy;
          continue;
        }
        if (excluded(person_exclude, ii)) {
          llxx = llyy;
          ii++;
          continue;
        }
        jj = (int)(line_locs[ii++] - llxx);
        llxx = llyy;
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
	      if (nn % 2) {
		nn = 3;
	      } else if (!nn) {
		nn = 1;
	      }
	    } else if (cc != marker_alleles[jj * 4 + 1]) {
              nn = 2;
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
      marker_alleles_tmp = (char*)malloc(map_linect4 * 2 * sizeof(char));
      if (!marker_alleles_tmp) {
        goto wdist_ret_2;
      }
      cptr = marker_alleles_tmp;
      iwptr = (int*)malloc(map_linect4 * 2 * sizeof(int));
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
        *bufptr++ = marker_alleles[ii * 2 + 1];
        *bufptr++ = '\t';
        *bufptr++ = marker_alleles[ii * 2];
        *bufptr++ = '\n';
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
      strcpy(outname_end, ".bed.tmp");
      pedfile = fopen(outname, "rb");
      if (!pedfile) {
        printf("Error: Failed to open %s.\n", outname);
      }
      binary_files = 1;
    }
  }

  ii = ((ped_linect - person_exclude_ct) * ((ped_linect - person_exclude_ct) - 1)) / 2;
  if (exponent == 0.0) {
    dists_alloc = ii * sizeof(int);
    idists = (int*)wkspace;
    memset(idists, 0, dists_alloc);
  } else {
    dists_alloc = ii * sizeof(double);
    dists = (double*)wkspace;
    dist_ptr = dists;
    do {
      *dist_ptr++ = 0.0;
    } while (--ii);
  }
  ped_geno = &(wkspace[dists_alloc]);
  ped_geno_size = malloc_size_mb * 1048576 - dists_alloc;

  if (binary_files) {
    if (snp_major) {
      fseeko(pedfile, 3, SEEK_SET);
      ii = 0; // current SNP index
      while (ii < map_linect) {
        for (jj = 0; jj < 4; jj++) {
          maf_buf[jj] = 0.5;
        }
        jj = 0; // actual SNPs read
        while ((jj < 4) && (ii < map_linect)) {
          while (excluded(marker_exclude, ii)) {
            ii++;
          }
          if (fread(&(pedbuf[jj * ped_linect4]), 1, ped_linect4, pedfile) < map_linect4) {
            retval = RET_READ_FAIL;
            goto wdist_ret_2;
          }
          maf_buf[jj] = ((double)marker_allele_cts[ii * 2 + 1]) / ((double)(marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1]));
          ii++;
          jj++;
        }
        if (jj < 4) {
          memset(&(pedbuf[jj * ped_linect4]), 0, (4 - jj) * ped_linect4);
        }
        gptr = ped_geno;
        for (jj = 0; jj < ped_linect; jj++) {
          if (!excluded(person_exclude, jj)) {
            kk = (jj % 4) * 2;
            *gptr++ = ((pedbuf[jj / 4] >> kk) & 3) | (((pedbuf[ped_linect4 + jj / 4] >> kk) & 3) << 2) | (((pedbuf[ped_linect4 * 2 + jj / 4] >> kk) & 3) << 4) | (((pedbuf[ped_linect4 * 3 + jj / 4] >> kk) & 3) << 6);
          }
        }
        if (exponent == 0.0) {
        } else {
          fill_weights(weights, maf_buf, exponent);
        }
      }
    } else {
    }
  } else {
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
  if (calculation_type == CALC_DISTANCE) {
    if (exponent == 0.0) {
      iwptr = idists;
      for (ii = 1; ii < ped_linect; ii += 1) {
        for (jj = 0; jj < ii; jj += 1) {
          if (jj) {
            fprintf(outfile, "\t%d", *iwptr++);
          } else {
            fprintf(outfile, "%d", *iwptr++);
          }
        }
        fprintf(outfile, "\n");
      }
      retval = RET_SUCCESS;
    } else {
      wptr = dists;
      for (ii = 1; ii < ped_linect; ii += 1) {
        for (jj = 0; jj < ii; jj += 1) {
          if (jj) {
            fprintf(outfile, "\t%lf", *wptr++);
          } else {
            fprintf(outfile, "%lf", *wptr++);
          }
        }
        fprintf(outfile, "\n");
      }
      retval = RET_SUCCESS;
    }
    if (retval == RET_SUCCESS) {
      printf("Distances written to %s.\n", outname);
    }
    } */
  retval = RET_SUCCESS;

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

int main(int argc, char* argv[]) {
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
  int ped_col_1 = 1;
  int ped_col_34 = 1;
  int ped_col_5 = 1;
  int ped_col_6 = 1;
  int ped_col_7 = 0;
  int mpheno_col = 0;
  char* phenoname_str = NULL;
  int affection_01 = 0;
  double exponent = 0.0;
  double min_maf = 0.01;
  double geno_thresh = 0.1;
  double mind_thresh = 0.1;
  int cur_arg = 1;
  int calculation_type = CALC_NONE;
  char* bubble;
  int threads = 2;
  int filter_type = 0;
  int mfilter_col = 0;
  int tail_pheno = 0;
  int prune = 0;
  int missing_pheno = -9;
  unsigned char missing_geno = '0';
  double tail_bottom;
  double tail_top;
  int ii;
  int jj;
  unsigned long int rseed = 0;
  strcpy(mapname, "wdist.map");
  strcpy(pedname, "wdist.ped");
  famname[0] = '\0';
  phenoname[0] = '\0';
  filtername[0] = '\0';
  strcpy(outname, "wdist");
  while (cur_arg < argc) {
    argptr = argv[cur_arg];
    if (!strcmp(argptr, "--file")) {
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
      if (sscanf(argv[cur_arg + 1], "%lf", &tail_top) != 1) {
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
      if (malloc_size_mb < 32) {
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
      threads = ii;
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
    } else if (!strcmp(argptr, "--maf")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --maf parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 1], "%lf", &min_maf) != 1) {
        printf("Error: Invalid --maf parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if ((min_maf < 0.0) || (min_maf > 0.5)) {
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
    } else if (!strcmp(argptr, "--distance")) {
      if (cur_arg != argc - 1) {
	printf("Error: invalid parameter after --distance.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg = argc;
      calculation_type = CALC_DISTANCE;
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
  if (exponent == 0.0) {
    for (ii = 0; ii < 256; ii += 1) {
      iwt[ii] = 8;
      jj = ii;
      while (jj) {
        if (jj % 2 == 1) {
          iwt[ii] -= 1;
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
  retval = wdist(pedname, mapname, famname, phenoname, filtername, makepheno_str, filter_type, filterval, mfilter_col, ped_col_1, ped_col_34, ped_col_5, ped_col_6, ped_col_7, (char)missing_geno, missing_pheno, mpheno_col, phenoname_str, prune, affection_01, threads, exponent, min_maf, geno_thresh, mind_thresh, tail_pheno, tail_bottom, tail_top, outname, calculation_type);
  // gsl_rng_free(rg);
  free(wkspace);
  return dispmsg(retval);
}
