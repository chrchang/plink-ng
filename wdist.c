#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#define RET_SUCCESS 0
#define RET_HELP 1
#define RET_NOMEM 2
#define RET_OPENFAIL 3
#define RET_INVALID_FORMAT 4
#define RET_CALC_NOT_YET_SUPPORTED 5
#define RET_INVALID_CMDLINE 6
#define RET_WRITE_FAIL 7

#define CALC_NONE -1
#define CALC_DISTANCE 0

#define MAX_THREADS 64
#define MAPBUFSIZE 64
#define PEDBUFBASE 1024
#define FNAMESIZE 2048
#define MALLOC_DEFAULT_MB 1984
#define IDLENGTH 16

const char errstr_map_format[] = "Error: Improperly formatted .map file.\n";
const char errstr_ped_format[] = "Error: Improperly formatted .ped file.\n";
int iwt[256];
char* wkspace;
long long malloc_size_mb = MALLOC_DEFAULT_MB;

// TODO:
// distance MAF histograms
// regression coefficients [delete-d jackknife errors]

int dispmsg(int retval) {
  switch(retval) {
  case RET_HELP:
    printf(
"WDIST weighted genetic distance calculator, v0.2.0 (15 July 2012)\n"
"Christopher Chang (chrchang523@gmail.com), BGI Cognitive Genomics Lab\n\n"
"wdist <flags> {calculation}\n\n"
"Supported flags:\n"
"  --file [prefix]  : Specify prefix for .ped and .map files (default 'wdist').\n"
"  --ped [filename] : Specify name of .ped file.\n"
"  --map [filename] : Specify name of .map file.\n"
"  --no-fid         : .ped file does not contain column 1 (family ID).\n"
"  --no-parents     : .ped file does not contain columns 3-4 (parents).\n"
"  --no-sex         : .ped file does not contain column 5 (sex).\n"
"  --liability      : .ped file does contain liability (column 7).\n"
"  --bfile [prefix] : Specify .bed/.bim/.fam prefix (default 'wdist').\n"
"  --bed [filename] : Specify .bed file.\n"
"  --bim [filename] : Specify .bim file.\n"
"  --fam [filename] : Specify .fam file.\n"
"  --out [prefix]   : Specify prefix for output files (default 'wdist').\n"
"  --silent         : Suppress output to console.\n"
"  --1              : .ped file affection phenotypes are interpreted as\n"
"                     -9 = missing, 0 = unaffected, 1 = affected (instead of\n"
"                     -9/0 = missing, 1 = unaffected, 2 = affected)\n"
    // --map3 implicitly supported via autodetection
"  --maf [val]      : Minor allele frequency minimum threshold (default 0.01).\n"
"  --geno [val]     : Maximum per-SNP missing (default 0.1).\n"
"  --mind [val]     : Maximum per-person missing (default 0.1).\n"
"  --memory [val]   : Size, in MB, of initial malloc attempt (default 1984).\n"
"  --threads [val]  : Maximum number of concurrent threads (default 2).\n"
"  --exponent [val] : When computing genetic distances, each locus has a weight\n"
"                     of (2q(1-q))^{-val}, where q is the observed MAF.\n\n"
"  --keep [filename]\n"
"  --remove [filename]\n"
"  --filter [filename] [val]: Keep/remove/filter individuals (see PLINK\n"
"                             documentation).\n\n"
"Supported calculations:\n"
"  --distance [--keep-missing-pheno]\n"
"    Outputs a lower-triangular table of (weighted) genetic distances.\n"
"    The first row contains a single number with the distance between the first\n"
"    two genotypes, the second row has the {genotype 1-genotype 3} and\n"
"    {genotype 2-genotype 3} distances in that order, etc.\n"
"    If modified by the --keep-missing-pheno flag, rows with 'missing'\n"
"    phenotype are not thrown out.\n\n"    
"  --groupdist [d] [iters]\n"
"  --groupdist [d] [iters] [Ltop] [Hbottom]\n"
"    Two-group genetic distance analysis, using delete-d jackknife with the\n"
"    requested number of iterations.  For a quantitative trait, one group has\n"
"    phenotype value less than or equal to Ltop, and the other has phenotype\n"
"    greater than Hbottom.  (Not actually implemented yet, but soon will be.)\n"
    // "  --histogram\n"
    // "  --regress\n\n"
	   );
    break;
  case RET_NOMEM:
    printf("Error: Out of memory.\n");
    break;
  case RET_WRITE_FAIL:
    printf("Error: File write failure.\n");
    break;
  }
  return retval;
}

int max_id_len = 4;
char* id_buf = NULL;

int qsort_partition(char* lptr, int min_idx, int max_idx, int pivot_idx) {
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

void qsort_str(char* lptr, int min_idx, int max_idx) {
  int pivot_idx;
  if (max_idx > min_idx) {
    pivot_idx = qsort_partition(lptr, min_idx, max_idx, (min_idx + max_idx) / 2);
    qsort_str(lptr, min_idx, pivot_idx - 1);
    qsort_str(lptr, pivot_idx + 1, max_idx);
  }
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

void fill_weights(double* weights, int* marker_allele_cts, double exponent, double min_maf) {
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int denom;
  double maf;
  double wt[4];
  double twt[3];
  for (ii = 0; ii < 4; ii += 1) {
    denom = marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1];
    if (denom) {
      maf = ((double)marker_allele_cts[ii * 2]) / ((double)denom);
      if (maf > 0.5) {
        maf = 1.0 - maf;
      }
      if (maf < min_maf) {
        maf = min_maf;
      }
    } else {
      maf = 0.5;
    }
    wt[ii] = pow(2.0 * maf * (1.0 - maf), -exponent);
  }
  nn = 0;
  for (ii = 0; ii < 4; ii += 1) {
    if (ii % 2) {
      twt[0] = wt[3];
    } else if (ii) {
      twt[0] = wt[3] * 2;
    } else {
      twt[0] = 0;
    }
    for (jj = 0; jj < 4; jj += 1) {
      if (jj % 2) {
        twt[1] = twt[0] + wt[2];
      } else if (jj) {
        twt[1] = twt[0] + 2 * wt[2];
      } else {
        twt[1] = twt[0];
      }
      for (kk = 0; kk < 4; kk += 1) {
        if (kk % 2) {
          twt[2] = twt[1] + wt[1];
	} else if (kk) {
          twt[2] = twt[1] + 2 * wt[1];
        } else {
          twt[2] = twt[1];
        }
        for (mm = 0; mm < 4; mm += 1) {
          if (mm % 2) {
            *weights++ = twt[2] + wt[0];
          } else if (mm) {
            *weights++ = twt[2] + 2 * wt[0];
          } else {
            *weights++ = twt[2];
	  }
        }
      }
    }
  }
}

void exclude(char* exclude_arr, int loc, int* exclude_ct) {
  int maj = loc / 8;
  int min = 1 << (loc % 8);
  if (!(exclude_arr[maj] & min)) {
    exclude_arr[maj] |= min;
    *exclude_ct += 1;
  }
}

int wdist(char* pedname, char* mapname, char* famname, char* filtername, int filter_type, char* filterval, int ped_col_1, int ped_col_34, int ped_col_5, int ped_col_7, int affection_01, int threads, double exponent, double min_maf, double geno_thresh, double mind_thresh, char* outname, int calculation_type, char* ep1, char* ep2, char* ep3, char* ep4) {
  FILE* pedfile = NULL;
  FILE* mapfile = NULL;
  FILE* famfile = NULL;
  FILE* outfile = NULL;
  FILE* filterfile = NULL;
  FILE* bedtmpfile = NULL;
  FILE* bimtmpfile = NULL;
  FILE* famtmpfile = NULL;
  int map_linect = 0;
  int map_linect4;
  char* outname_end;
  char mapbuf[MAPBUFSIZE];
  char* pedbuf = NULL;
  char* marker_exclude = NULL;
  long long* line_locs = NULL;
  int max_people;
  int pedbuflen;
  int ped_recalc_len = 0;
  char* fgets_return;
  int ped_linect = 0;
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int* marker_chrom = NULL;
  id_string* marker_id = NULL;
  int* marker_pos = NULL;
  char* marker_alleles = NULL;
  char* marker_alleles_tmp = NULL;
  int* marker_allele_cts = NULL;
  char* bufptr;
  char* person_exclude = NULL;
  int retval;
  int map_cols = 3;
  int affection = 0;
  double* ped_pheno = NULL;
  int* ped_phenoi = NULL;
  char* ped_geno = NULL;
  char* gptr;
  char* gptr2;
  int* iwptr;
  int* iptr;
  double* wptr;
  int* idists;
  double* dists;
  double weights[256];
  int dists_alloc = 0;
  int geno_window_size;
  int windows_reqd;
  int window_num;
  int smallbuflen;
  int cur_window_size;
  int cur_window_size4;
  char cc;
  double* dist_ptr;
  int last_tell;
  int binary_files = 0;
  int maf_int_thresh;
  int geno_int_thresh;
  int mind_int_thresh;
  int marker_exclude_ct = 0;
  int person_exclude_ct = 0;
  int filter_lines = 0;
  char* id_list = NULL;

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
  if (filter_type) {
    filterfile = fopen(filtername, "r");
    if (!filterfile) {
      printf("Error: Failed to open %s.\n", filtername);
      goto wdist_ret_0;
    }
  }
  outname_end = outname;
  while (*outname_end++);
  strcpy(outname_end, ".dist");
  outfile = fopen(outname, "w");
  if (!outfile) {
    printf("Error: Failed to open %s.\n", outname);
    goto wdist_ret_0;
  }
  retval = RET_NOMEM;

  if (exponent == 0.0) {
    max_people = (int)sqrt((malloc_size_mb * 1048576) / sizeof(int));
  } else {
    max_people = (int)sqrt((malloc_size_mb * 1048576) / sizeof(double));
  }
  while (fgets(mapbuf, MAPBUFSIZE, mapfile) != NULL) {
    if (mapbuf[0] > ' ') {
      if (!map_linect) {
	bufptr = next_item(mapbuf);
	bufptr = next_item(bufptr);
	bufptr = next_item(bufptr);
        if (binary_files) {
          bufptr = next_item(bufptr);
          bufptr = next_item(bufptr);
        }
	if (!bufptr) {
	  retval = RET_INVALID_FORMAT;
	  printf(errstr_map_format);
	  goto wdist_ret_4;
	}
	if (*bufptr > ' ') {
	  map_cols = 4;
	}
      }
      map_linect += 1;
    }
  }
  if (!map_linect) {
    retval = RET_INVALID_FORMAT;
    printf("Error: No markers in .map file.");
    goto wdist_ret_0;
  }
  rewind(mapfile);
  marker_exclude = (char*)malloc(((map_linect + 7) / 8) * sizeof(char));
  if (!marker_exclude) {
    goto wdist_ret_4;
  }
  memset(marker_exclude, 0, ((map_linect + 7) / 8) * sizeof(char));
  if (binary_files) {
    pedbuflen = (map_linect + 3) / 4;
  } else {
    pedbuflen = map_linect * 5 + PEDBUFBASE;
  }
  pedbuf = (char*)malloc(pedbuflen * sizeof(char));
  if (!pedbuf) {
    goto wdist_ret_4;
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

  for (ii = 0; ii < map_linect; ii += 1) {
    do {
      fgets(mapbuf, MAPBUFSIZE, mapfile);
    } while (mapbuf[0] <= ' ');
    // marker_chrom[ii] = marker_code(mapbuf);
    bufptr = next_item(mapbuf);
    if (!bufptr) {
      retval = RET_INVALID_FORMAT;
      printf(errstr_map_format);
      goto wdist_ret_4;
    }
    // cur_item(marker_id[ii], bufptr);
    bufptr = next_item(bufptr);
    if (map_cols == 4) {
      bufptr = next_item(bufptr);
    }
    if (!bufptr) {
      retval = RET_INVALID_FORMAT;
      printf(errstr_map_format);
      goto wdist_ret_4;
    }
    if (*bufptr == '-') {
      exclude(marker_exclude, ii, &marker_exclude_ct);
    }
    // marker_pos[ii] = atoi(bufptr);
  }
  if (filter_type) {
    pedbuf[pedbuflen - 1] = ' ';
    if (filter_type == 3) {
      jj = strlen(filterval);
    }
    while (fgets(pedbuf, pedbuflen, filterfile) != NULL) {
      if (*pedbuf == '\n') {
        continue;
      }
      bufptr = pedbuf;
      ii = 0;
      while ((*bufptr != ' ') && (*bufptr != '\t') && (*bufptr != '\n') && (*bufptr)) {
        ii++;
        bufptr++;
      }
      while ((*bufptr == ' ') || (*bufptr == '\t')) {
        bufptr++;
      }
      while ((*bufptr != ' ') && (*bufptr != '\t') && (*bufptr != '\n') && (*bufptr)) {
        ii++;
        bufptr++;
      }
      if ((ii == 0) || (*bufptr == '\0')) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Improperly formatted filter file.\n");
        goto wdist_ret_4;
      }
      ii += 2;
      if (ii > max_id_len) {
        max_id_len = ii;
      }
      if (filter_type == 3) {
        while ((*bufptr == ' ') || (*bufptr == '\t')) {
          bufptr++;
        }
        if ((*bufptr == '\n') || (*bufptr == '\0')) {
	  retval = RET_INVALID_FORMAT;
	  printf("Error: Improperly formatted filter file.\n");
	  goto wdist_ret_4;
        }
        if (!strncmp(filterval, bufptr, jj)) {
          filter_lines += 1;
        }
      } else {
        filter_lines += 1;
      }
      while (!pedbuf[pedbuflen - 1]) {
        pedbuf[pedbuflen - 1] = ' ';
        if (pedbuf[pedbuflen - 2] == '\n') {
          break;
        }
	fgets(pedbuf, pedbuflen, filterfile);
      }
    }
    id_list = (char*)malloc(max_id_len * filter_lines * sizeof(char));
    if (!id_list) {
      goto wdist_ret_4;
    }
    id_buf = (char*)malloc(max_id_len * sizeof(char));
    if (!id_buf) {
      goto wdist_ret_4;
    }
    rewind(filterfile);
    ii = 0;
    while (fgets(pedbuf, pedbuflen, filterfile) != NULL) {
      if (*pedbuf == '\n') {
        continue;
      }
      gptr = &(id_list[ii * max_id_len]);
      bufptr = pedbuf;
      while ((*bufptr != ' ') && (*bufptr != '\t')) {
        *gptr++ = *bufptr++;
      }
      *gptr++ = ' ';
      while ((*bufptr == ' ') || (*bufptr == '\t')) {
        bufptr++;
      }
      while ((*bufptr != ' ') && (*bufptr != '\t') && (*bufptr != '\n')) {
        *gptr++ = *bufptr++;
      }
      if (filter_type == 3) {
        while ((*bufptr == ' ') || (*bufptr == '\t')) {
          bufptr++;
        }
        if (!strncmp(filterval, bufptr, jj)) {
          ii++;
          *gptr = '\0';
        }
      } else {
        ii++;
        *gptr = '\0';
      }
      if (ii == filter_lines) {
        break;
      }
      while (!pedbuf[pedbuflen - 1]) {
        pedbuf[pedbuflen - 1] = ' ';
        if (pedbuf[pedbuflen - 2] == '\n') {
          break;
        }
	fgets(pedbuf, pedbuflen, filterfile);
      }
    }
    qsort_str(id_list, 0, filter_lines - 1);

    // debug
    for (ii = 0; ii < filter_lines; ii++) {
      printf("%s\n", &(id_list[ii * max_id_len]));
    }
  }
  if (binary_files) {
    // ...apply filter
  } else {
    pedbuf[pedbuflen - 1] = ' ';
    last_tell = 0;
    line_locs = (long long*)malloc(max_people * sizeof(long long));
    if (!line_locs) {
      goto wdist_ret_4;
    }
    while (fgets(pedbuf, pedbuflen, pedfile) != NULL) {
      if (pedbuf[0] > ' ') {
	if (pedbuf[0] != '#') {
	  bufptr = next_item(pedbuf);
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
	    goto wdist_ret_4;
	  }
	  if (!ped_linect) {
	    if (*bufptr == '-') {
	      if ((bufptr[1] == '9') && ((bufptr[2] == ' ') || (bufptr[2] == '\t'))) {
		affection = 1;
	      }
	    } else if (((*bufptr >= '0') && (*bufptr <= '2')) && ((bufptr[1] == ' ') || (bufptr[1] == '\t'))) {
	      affection = 1;
	    }
	  }
	  if (affection) {
	    if (affection_01) {
	      if ((*bufptr == '0') || (*bufptr == '1')) {
		line_locs[ped_linect] = last_tell + (bufptr - pedbuf);
		ped_linect += 1;
	      }
	    } else if ((*bufptr == '1') || (*bufptr == '2')) {
	      line_locs[ped_linect] = last_tell + (bufptr - pedbuf);
	      ped_linect += 1;
	    }
	  } else {
	    line_locs[ped_linect] = last_tell + (bufptr - pedbuf);
	    ped_linect += 1;
	  }
	}
      }
      if (!pedbuf[pedbuflen - 1]) {
        pedbuf[pedbuflen - 1] = ' ';
        if (pedbuf[pedbuflen - 2] == '\n') {
          break;
        }
	ii = 0;
	do {
	  ii += pedbuflen - 1;
	  pedbuf[pedbuflen - 1] = ' ';
	  fgets_return = fgets(pedbuf, pedbuflen, pedfile);
	} while (fgets_return && !pedbuf[pedbuflen - 1] && (pedbuf[pedbuflen - 2] != '\n'));
	ii += strlen(pedbuf) + 1;
	if (ii > ped_recalc_len) {
	  ped_recalc_len = ii;
	}
      }
      last_tell = ftell(pedfile);
    }
    if (ped_linect < 2) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Less than two valid people in .ped file.\n");
      goto wdist_ret_4;
    } else if (ped_linect > max_people) {
      retval = RET_INVALID_FORMAT;
      printf("Too many people in .ped file for this algorithm.\n");
      goto wdist_ret_4;
    }
    rewind(pedfile);
    if (ped_recalc_len) {
      free(pedbuf);
      pedbuflen = ped_recalc_len;
      pedbuf = (char*)malloc(pedbuflen * sizeof(char));
      if (!pedbuf) {
	goto wdist_ret_4;
      }
    }

    marker_alleles = (char*)malloc(map_linect * 4 * sizeof(char));
    if (!marker_alleles) {
      goto wdist_ret_4;
    }
    memset(marker_alleles, 0, map_linect * 4 * sizeof(char));
    marker_allele_cts = (int*)malloc(map_linect * 4 * sizeof(int));
    if (!marker_allele_cts) {
      goto wdist_ret_4;
    }
    memset(marker_allele_cts, 0, map_linect * 4 * sizeof(int));
    if (affection) {
      ped_phenoi = (int*)malloc(ped_linect * sizeof(int));
      if (!ped_phenoi) {
	goto wdist_ret_4;
      }
    } else {
      ped_pheno = (double*)malloc(ped_linect * sizeof(double));
      if (!ped_pheno) {
	goto wdist_ret_4;
      }
    }
    person_exclude = (char*)malloc(((ped_linect + 7) / 8) * sizeof(char));
    if (!person_exclude) {
      goto wdist_ret_4;
    }
    memset(person_exclude, 0, ((ped_linect + 7) / 8) * sizeof(char));

    ped_geno = wkspace;
    maf_int_thresh = map_linect - (int)((1.0 - min_maf) * map_linect);
    geno_int_thresh = map_linect - (int)(geno_thresh * map_linect);
    mind_int_thresh = (int)(mind_thresh * map_linect);
    for (ii = 0; ii < ped_linect; ii += 1) {
      fseek(pedfile, line_locs[ii], SEEK_SET);
      fgets(pedbuf, pedbuflen, pedfile);
      bufptr = pedbuf;
      if (affection) {
	if (affection_01) {
	  if (*bufptr == '0') {
	    ped_phenoi[ii] = 0;
	  } else {
	    ped_phenoi[ii] = 1;
	  }
	} else {
	  if (*bufptr == '1') {
	    ped_phenoi[ii] = 0;
	  } else {
	    ped_phenoi[ii] = 1;
	  }
	}
      } else {
	if (sscanf(bufptr, "%lf", &(ped_pheno[ii])) != 1) {
	  retval = RET_INVALID_FORMAT;
	  printf(errstr_ped_format);
	  goto wdist_ret_4;
	}
      }
      bufptr = next_item(bufptr);
      if (ped_col_7) {
	bufptr = next_item(bufptr);
      }
      line_locs[ii] += (bufptr - pedbuf);
      mm = 0; // number of missing
      for (jj = 0; jj < map_linect; jj += 1) {
        for (kk = 0; kk < 2; kk++) {
          cc = *bufptr;
          if (cc == '0') {
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
            goto wdist_ret_4;
	  } else if (marker_alleles[jj * 4 + 2]) {
            marker_alleles[jj * 4 + 3] = cc;
	    marker_allele_cts[jj * 4 + 3] = 1;
	  } else if (marker_alleles[jj * 4 + 1]) {
	    marker_alleles[jj * 4 + 2] = cc;
	    marker_allele_cts[jj * 4 + 2] = 1;
	  } else if (marker_alleles[jj * 3]) {
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
	  if (*bufptr == '\0') {
	    retval = RET_INVALID_FORMAT;
	    printf(errstr_ped_format);
	    goto wdist_ret_4;
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
        goto wdist_ret_4;
      }
    }
    if (person_exclude_ct > (ped_linect - 2)) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Too many people fail QC.\n");
      goto wdist_ret_4;
    } else if (marker_exclude_ct == map_linect) {
      retval = RET_INVALID_FORMAT;
      printf("Error: All markers fail QC.\n");
      goto wdist_ret_4;
    }

    ii = ((ped_linect - person_exclude_ct) * ((ped_linect - person_exclude_ct) - 1)) / 2;
    if (exponent == 0.0) {
      dists_alloc = ii * sizeof(int);
    } else {
      dists_alloc = ii * sizeof(double);
    }
    
    geno_window_size = (malloc_size_mb * 1048576 - dists_alloc) / (ped_linect - person_exclude_ct);
    map_linect4 = (map_linect - marker_exclude_ct + 3) / 4;
    windows_reqd = ((map_linect4 - 1) / geno_window_size) + 1;
    if (windows_reqd == 1) {
      geno_window_size = map_linect4;
      smallbuflen = pedbuflen;
    } else {
      printf(".ped file too large for direct read.  Converting to binary...");
      strcpy(outname_end, ".bed.tmp");
      bedtmpfile = fopen(outname, "wb");
      if (!bedtmpfile) {
        retval = RET_OPENFAIL;
        printf(" failed to open %s.\n", outname);
	goto wdist_ret_5;
      }
      strcpy(outname_end, ".bim.tmp");
      bimtmpfile = fopen(outname, "w");
      if (!bimtmpfile) {
        retval = RET_OPENFAIL;
        printf(" failed to open %s.\n", outname);
	goto wdist_ret_6;
      }
      strcpy(outname_end, ".fam.tmp");
      famtmpfile = fopen(outname, "w");
      if (!famtmpfile) {
        retval = RET_OPENFAIL;
        printf(" failed to open %s.\n", outname);
	goto wdist_ret_6;
      }
      if (3 != fwrite("L\x1b", 1, 3, bedtmpfile)) {
	retval = RET_WRITE_FAIL;
        printf("\n");
	goto wdist_ret_7;
      }
      for (ii = 0; ii < ped_linect; ii += 1) {
        if (person_exclude[ii / 8] & (1 << (ii % 8))) {
          continue;
        }
	memset(ped_geno, 0, map_linect4 * sizeof(char));
        gptr = ped_geno;
	fseek(pedfile, line_locs[ii], SEEK_SET);
	fgets(pedbuf, pedbuflen, pedfile);
	bufptr = pedbuf;
	for (jj = 0; jj < map_linect; jj += 1) {
          if (marker_exclude[jj / 8] & (1 << (jj % 8))) {
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
	  mm = jj % 4;
	  nn = 0;
	  for (kk = 0; kk < 2; kk++) {
	    cc = *bufptr;
	    if (cc == marker_alleles[jj * 2 + 1]) {
	      if (nn % 2) {
		nn = 3;
	      } else if (!nn) {
		nn = 1;
	      }
	    } else if (cc != marker_alleles[jj * 2]) {
              nn = 2;
	    }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
	  }
	  *gptr |= (nn << (mm * 2));
	  if (mm == 3) {
	    gptr++;
	  }
	}
	if (map_linect4 != fwrite(ped_geno, 1, map_linect4, bedtmpfile)) {
	  retval = RET_WRITE_FAIL;
	  goto wdist_ret_7;
	}
      }
      marker_alleles_tmp = (char*)malloc(map_linect4 * 2 * sizeof(char));
      if (!marker_alleles_tmp) {
        goto wdist_ret_7;
      }
      gptr = marker_alleles_tmp;
      iwptr = (int*)malloc(map_linect4 * 2 * sizeof(int));
      iptr = iwptr;
      for (ii = 0; ii < map_linect; ii += 1) {
        if (marker_exclude[ii / 8] & (1 << (ii % 8))) {
          continue;
        }
        *gptr++ = marker_alleles[ii * 4];
        *gptr++ = marker_alleles[ii * 4 + 1];
        *iptr++ = marker_allele_cts[ii * 4];
        *iptr++ = marker_allele_cts[ii * 4 + 1];
      }
      free(marker_alleles);
      marker_alleles = marker_alleles_tmp;
      marker_alleles_tmp = NULL;
      marker_allele_cts = iwptr;
      // TODO: write .bim and .fam
      printf(" done.\n");
      fclose(pedfile);
      fclose(bedtmpfile);
      bedtmpfile = NULL;
      strcpy(outname_end, ".bed.tmp");
      pedfile = fopen(outname, "rb");
      if (!pedfile) {
        printf(" failed to open %s.\n", outname);
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

  /*
  if (binary_files) {
  } else {
  }

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
      fseek(pedfile, line_locs[ii], SEEK_SET);
      fgets(pedbuf, smallbuflen, pedfile);
      bufptr = pedbuf;
      if (!window_num) {
	if (affection) {
	  if (affection_01) {
            if (*bufptr == '0') {
              ped_phenoi[ii] = 0;
            } else {
              ped_phenoi[ii] = 1;
            }
	  } else {
            if (*bufptr == '1') {
              ped_phenoi[ii] = 0;
            } else {
              ped_phenoi[ii] = 1;
            }
	  }
        } else {
	  if (sscanf(bufptr, "%lf", &(ped_pheno[ii])) != 1) {
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

 wdist_ret_9:
 wdist_ret_8:
 wdist_ret_7:
 wdist_ret_6:
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
 wdist_ret_5:
 wdist_ret_4:
  if (person_exclude) {
    free(person_exclude);
  }
  if (ped_phenoi) {
    free(ped_phenoi);
  }
  if (ped_pheno) {
    free(ped_pheno);
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
  if (marker_exclude) {
    free(marker_exclude);
  }
  if (pedbuf) {
    free(pedbuf);
  }
  if (marker_pos) {
    free(marker_pos);
  }
  if (marker_id) {
    free(marker_id);
  }
  if (marker_chrom) {
    free(marker_chrom);
  }
 wdist_ret_0:
  if (outfile) {
    fclose(outfile);
  }
  if (filterfile) {
    fclose(filterfile);
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
  char filtername[FNAMESIZE];
  char* filterval;
  char* argptr;
  int retval;
  int load_params = 0; // describes what file parameters have been provided
  int ped_col_1 = 1;
  int ped_col_34 = 1;
  int ped_col_5 = 1;
  int ped_col_7 = 0;
  int affection_01 = 0;
  double exponent = 0.0;
  double min_maf = 0.01;
  double geno_thresh = 0.1;
  double mind_thresh = 0.1;
  int cur_arg = 1;
  int calculation_type = CALC_NONE;
  char* ep1 = NULL;
  char* ep2 = NULL;
  char* ep3 = NULL;
  char* ep4 = NULL;
  char* bubble;
  int threads = 2;
  int filter_type = 0;
  int ii;
  int jj;
  strcpy(mapname, "wdist.map");
  strcpy(pedname, "wdist.ped");
  famname[0] = '\0';
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
        printf("Error: Missing --file parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        printf("Error: Missing --ped parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        printf("Error: Missing --map parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        printf("Error: Missing --bfile parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        printf("Error: Missing --bed parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        printf("Error: Missing --bim parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        printf("Error: Missing --fam parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --fam parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(famname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--no-sex")) {
      ped_col_5 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--no-parents")) {
      ped_col_34 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--no-fid")) {
      ped_col_1 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--liability")) {
      ped_col_7 = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--1")) {
      affection_01 = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--keep")) {
      if (cur_arg == argc - 1) {
        printf("Error: missing --keep parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
      }
      if (filter_type > 0) {
        printf("Error: More than one filtering flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      filter_type = 1;
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--keep filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(filtername, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--remove")) {
      if (cur_arg == argc - 1) {
        printf("Error: missing --remove parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
      }
      if (filter_type > 0) {
        printf("Error: More than one filtering flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      filter_type = 2;
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--remove filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(filtername, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--filter")) {
      if (cur_arg > argc - 2) {
        printf("Error: Not enough --filter parameters.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
      }
      if (filter_type > 0) {
        printf("Error: More than one filtering flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      filter_type = 3;
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--filter filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(filtername, argv[cur_arg + 1]);
      filterval = argv[cur_arg + 2];
      cur_arg += 3;
    } else if (!strcmp(argptr, "--memory")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --memory parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
      }
      malloc_size_mb = atoi(argv[cur_arg + 1]);
      if (malloc_size_mb < 32) {
        printf("Error: Invalid --memory parameter.\n");
        return dispmsg(RET_HELP);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--threads")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --threads parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        printf("Error: Missing --exponent parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
      }
      if (sscanf(argv[cur_arg + 1], "%lf", &exponent) != 1) {
        printf("Error: Invalid --exponent parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--maf")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --maf parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        printf("Error: Missing --geno parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        printf("Error: Missing --mind parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
    } else if (!strcmp(argptr, "--out")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --out parameter.  Displaying general-purpose help.\n\n");
        return dispmsg(RET_HELP);
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
        if ((cur_arg == argc - 2) && !strcmp(argv[cur_arg + 1], "--keep-missingpheno")) {
          
        } else {
          printf("Error: invalid parameter after --distance.  Displaying general-purpose help.\n");
          return dispmsg(RET_HELP);
        }
      }
      cur_arg = argc;
      calculation_type = CALC_DISTANCE;
    } else if (!strcmp(argptr, "--map3")) {
      printf("Note: --map3 flag unnecessary (.map file format is autodetected).\n");
      cur_arg += 1;
    } else {
      printf("Error: Invalid argument (%s).  Displaying general-purpose help.\n\n", argv[cur_arg]);
      return dispmsg(RET_HELP);
    }
  }
  if (calculation_type == CALC_NONE) {
    return dispmsg(RET_HELP);
  }
  if (exponent == 0.0) {
    for (ii = 0; ii < 256; ii += 1) {
      iwt[ii] = 0;
      jj = ii;
      while (jj) {
        if (jj % 2 == 1) {
          iwt[ii] += 1;
        } else if (jj % 4 == 2) {
          iwt[ii] += 2;
        }
        jj >>= 2;
      }
    }
  }
  bubble = (char*)malloc(33554432 * sizeof(char));
  if (!bubble) {
    return dispmsg(RET_NOMEM);
  }
  wkspace = (char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  if ((malloc_size_mb > MALLOC_DEFAULT_MB) && !wkspace) {
    printf("%lld MB malloc failed.  Using default allocation behavior.\n", malloc_size_mb);
    malloc_size_mb = MALLOC_DEFAULT_MB;
  }
  while (!wkspace) {
    if (malloc_size_mb > 96) {
      malloc_size_mb -= 64;
    } else {
      malloc_size_mb = 32;
    }
    wkspace = (char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  }
  free(bubble);

  // famname[0] indicates binary vs. text
  // filtername[0] indicates existence of filter
  retval = wdist(pedname, mapname, famname, filtername, filter_type, filterval, ped_col_1, ped_col_34, ped_col_5, ped_col_7, affection_01, threads, exponent, min_maf, geno_thresh, mind_thresh, outname, calculation_type, ep1, ep2, ep3, ep4);
  free(wkspace);
  return dispmsg(retval);
}
