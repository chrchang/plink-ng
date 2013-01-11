#include "wdist_common.h"

const char errstr_fopen[] = "Error: Failed to open %s.\n";
const char errstr_append[] = "\nFor more information, try 'wdist --help [flag name]' or 'wdist --help | more'.\n";
const char errstr_thread_create[] = "\nError: Failed to create thread.\n";

char tbuf[MAXLINELEN * 4 + 256];

FILE* logfile = NULL;
char logbuf[MAXLINELEN]; // safe sprintf buffer, if one is needed
int debug_on = 0;
int log_failed = 0;

void logstr(const char* ss) {
  if (log_failed) {
    if (debug_on) {
      printf("%s", ss);
    }
  } else if (fprintf(logfile, "%s", ss) < 0) {
    if (!log_failed) {
      if (debug_on) {
	printf("\nError: Debug logging failure.  Dumping to standard output:\n%s", ss);
      } else {
	printf("\nNote: Logging failure on:\n%s\nFurther logging will not be attempted in this run.\n", ss);
      }
      log_failed = 1;
    }
  }
}

void logprint(const char* ss) {
  logstr(ss);
  fputs(ss, stdout);
}

void logprintb() {
  logstr(logbuf);
  fputs(logbuf, stdout);
}

int fopen_checked(FILE** target_ptr, const char* fname, const char* mode) {
  *target_ptr = fopen(fname, mode);
  if (!(*target_ptr)) {
    sprintf(logbuf, errstr_fopen, fname);
    logprintb();
    return -1;
  }
  return 0;
}

int gzopen_checked(gzFile* target_ptr, const char* fname, const char* mode) {
  *target_ptr = gzopen(fname, mode);
  if (!(*target_ptr)) {
    sprintf(logbuf, errstr_fopen, fname);
    logprintb();
    return -1;
  }
  return 0;
}

// manually managed, very large stack
unsigned char* wkspace_base;
unsigned long wkspace_left;

unsigned char* wkspace_alloc(unsigned long size) {
  unsigned char* retval;
  if (wkspace_left < size) {
    return NULL;
  }
  size = CACHEALIGN(size);
  retval = wkspace_base;
  wkspace_base += size;
  wkspace_left -= size;
  return retval;
}

void wkspace_reset(void* new_base) {
  unsigned long freed_bytes = wkspace_base - (unsigned char*)new_base;
  wkspace_base = (unsigned char*)new_base;
  wkspace_left += freed_bytes;
}

char* item_end(char* sptr) {
  char cc;
  if (!sptr) {
    return NULL;
  }
  cc = *sptr;
  while ((cc != ' ') && (cc != '\t') && (cc != '\n') && (cc != '\r')) {
    if (!cc) {
      return NULL;
    }
    cc = *(++sptr);
  }
  return sptr;
}

char* item_endl(char* sptr) {
  if (!sptr) {
    return NULL;
  }
  while (!is_space_or_eoln(*sptr)) {
    sptr++;
  }
  return sptr;
}

char* item_endnn(char* sptr) {
  while (!is_space_or_eoln(*sptr)) {
    sptr++;
  }
  return sptr;
}

char* item_endnn2(char* sptr) {
  char cc = *sptr;
  while ((cc != ' ') && (cc != '\t') && cc) {
    cc = *(++sptr);
  }
  return sptr;
}

int intlen(int num) {
  int retval;
  if (num < 0) {
    num = -num;
    retval = 2;
  } else {
    retval = 1;
  }
  while (num > 9) {
    num /= 10;
    retval++;
  }
  return retval;
}

int strlen_se(char* ss) {
  int val = 0;
  while (!is_space_or_eoln(*ss++)) {
    val++;
  }
  return val;
}

int strcmp_se(char* s_read, const char* s_const, int len) {
  return memcmp(s_read, s_const, len) || (!is_space_or_eoln(s_read[len]));
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
  return skip_initial_spaces(sptr);
}

char* next_item_mult(char* sptr, unsigned int ct) {
  if (!sptr) {
    return NULL;
  }
  do {
    while ((*sptr != ' ') && (*sptr != '\t')) {
      if (!(*sptr)) {
	return NULL;
      }
      sptr++;
    }
    sptr = skip_initial_spaces(sptr);
  } while (--ct);
  return sptr;
}

void copy_item(char* writebuf, unsigned int* offset_ptr, char** prev_item_end_ptr) {
  char* item_start = skip_initial_spaces(*prev_item_end_ptr);
  char* item_end = item_endnn2(item_start);
  unsigned int slen = (item_end - item_start);
  unsigned int offset = *offset_ptr;
  memcpy(&(writebuf[offset]), item_start, slen);
  offset += slen;
  writebuf[offset++] = '\t';
  *offset_ptr = offset;
  *prev_item_end_ptr = item_end;
}

void set_bit(unsigned long* bit_arr, int loc, unsigned int* bit_set_ct_ptr) {
  int maj = loc / BITCT;
  unsigned long min = 1LU << (loc % BITCT);
  if (!(bit_arr[maj] & min)) {
    bit_arr[maj] |= min;
    *bit_set_ct_ptr += 1;
  }
}

void set_bit_sub(unsigned long* bit_arr, int loc, unsigned int* bit_unset_ct_ptr) {
  int maj = loc / BITCT;
  unsigned long min = 1LU << (loc % BITCT);
  if (!(bit_arr[maj] & min)) {
    bit_arr[maj] |= min;
    *bit_unset_ct_ptr -= 1;
  }
}

void clear_bit(unsigned long* exclude_arr, int loc, unsigned int* include_ct_ptr) {
  int maj = loc / BITCT;
  unsigned long min = 1LU << (loc % BITCT);
  if (exclude_arr[maj] & min) {
    exclude_arr[maj] -= min;
    *include_ct_ptr += 1;
  }
}

// unsafe if you don't know there's another included marker or person remaining
int next_non_set_unsafe(unsigned long* exclude_arr, unsigned int loc) {
  unsigned int idx = loc / BITCT;
  unsigned long ulii;
  exclude_arr = &(exclude_arr[idx]);
  ulii = (~(*exclude_arr)) >> (loc % BITCT);
  if (ulii) {
    return loc + __builtin_ctzl(ulii);
  }
  do {
    idx++;
  } while (*(++exclude_arr) == ~0LU);
  return (idx * BITCT) + __builtin_ctzl(~(*exclude_arr));
}

int next_set_unsafe(unsigned long* include_arr, unsigned int loc) {
  unsigned int idx = loc / BITCT;
  unsigned long ulii;
  include_arr = &(include_arr[idx]);
  ulii = (*include_arr) >> (loc % BITCT);
  if (ulii) {
    return loc + __builtin_ctzl(ulii);
  }
  do {
    idx++;
  } while (*(++include_arr) == 0);
  return (idx * BITCT) + __builtin_ctzl(*include_arr);
}

// human: 22, X, Y, XY, MT
// cow: 29, X, Y
// dog: 38, X, Y, XY
// horse: 31, X, Y
// mouse: 19, X, Y
// rice: 12 (haploid, not supported for now)
// sheep: 26, X, Y
const unsigned long long species_def_chrom_mask[] = {0x027fffffLLU, 0x3fffffffLLU, 0x27fffffffffLLU, 0xffffffffLLU, 0x000fffffLLU, 0LLU, 0x07ffffffLLU};
const unsigned long long species_autosome_mask[] = {0x007ffffeLLU, 0x3ffffffeLLU, 0x7ffffffffeLLU, 0xfffffffeLLU, 0x000ffffeLLU, 0LLU, 0x07fffffeLLU};
const unsigned long long species_valid_chrom_mask[] = {0x1000027fffffLLU, 0x3fffffffLLU, 0x127fffffffffLLU, 0xffffffffLLU, 0x000fffffLLU, 0LLU, 0x07ffffffLLU};
const char species_x_code[] = {23, 30, 39, 32, 20, -1, 27};
const char species_y_code[] = {24, 31, 40, 33, 21, -1, 28};
const char species_xy_code[] = {25, -1, 41, -1, -1, -1, -1};
const char species_mt_code[] = {26, -1, -1, -1, -1, -1, -1};
const char species_max_code[] = {26, 31, 41, 33, 21, 12, 28};
const unsigned long long species_haploid_mask[] = {}; // todo

int marker_code_raw(char* sptr) {
  // any character <= ' ' is considered a terminator
  int ii;
  if (sptr[1] > ' ') {
    if (sptr[2] > ' ') {
      return -1;
    }
    if ((sptr[0] == 'X') || (sptr[0] == 'x')) {
      if ((sptr[1] == 'Y') || (sptr[1] == 'y')) {
	return CHROM_XY;
      }
      return -1;
    }
    if ((sptr[0] == 'M') || (sptr[0] == 'm')) {
      if ((sptr[1] == 'T') || (sptr[1] == 't')) {
	return CHROM_MT;
      }
      return -1;
    }
    if ((sptr[0] >= '0') && (sptr[0] <= '9')) {
      if ((sptr[1] >= '0') && (sptr[1] <= '9')) {
        ii = ((sptr[0] - '0') * 10 + (sptr[1] - '0'));
	if (ii < MAX_POSSIBLE_CHROM) {
	  return ii;
	} else {
	  return -1;
	}
      } else {
	return -1;
      }
    }
  } else if ((sptr[0] >= '0') && (sptr[0] <= '9')) {
    return (sptr[0] - '0');
  } else if ((sptr[0] == 'X') || (sptr[0] == 'x')) {
    return CHROM_X;
  } else if ((sptr[0] == 'Y') || (sptr[0] == 'y')) {
    return CHROM_Y;
  } else if ((sptr[0] == 'M') || (sptr[0] == 'm')) {
    return CHROM_MT;
  }
  return -1;
}

int marker_code(unsigned int species, char* sptr) {
  // does not require string to be null-terminated, and does not perform
  // exhaustive error-checking
  int ii = marker_code_raw(sptr);
  if (ii >= MAX_POSSIBLE_CHROM) {
    switch (ii) {
    case CHROM_X:
      ii = species_x_code[species];
      break;
    case CHROM_Y:
      ii = species_y_code[species];
      break;
    case CHROM_XY:
      ii = species_xy_code[species];
      break;
    case CHROM_MT:
      ii = species_mt_code[species];
    }
  } else if (ii > species_max_code[species]) {
    return -1;
  }
  return ii;
}

// WDIST's natural sort uses the following logic:
// - All alphabetic characters act as if they are capitalized, except for
// tiebreaking purposes (where ASCII is used).
// - Numbers are compared by magnitude, with the exception of...
// - Numbers with leading zero(es).  If you're putting extraneous zeroes in
// front of IDs, we assume they're there to force particular items to be sorted
// earlier, rather than just appearing at random.  So, unlike many natural sort
// implementations, we sort 00200 < 021 < 20: all numbers with n leading zeroes
// are sorted before all numbers with (n-1) leading zeroes; magnitude only
// applies if the leading zero counts match.  This handles e.g. subbasement
// room numbering properly.
//
// This won't always do what you want if your IDs have variable-length decimals
// in them (e.g. it yields 0.99 < 0.101); if you don't want to fall back on
// ASCII sort, enforce a fixed number of digits after the decimal point.  Also
// note that ASCII sort is outright better for e.g. numbers represented in
// hexadecimal or base 36.  In principle, it's possible to reliably autodetect
// some of these cases (especially hexadecimal numbers beginning with "0x"),
// but that'll never be perfect so we just let the user toggle the sort method.
int strcmp_natural_scan_forward(const char* s1, const char* s2) {
  // assumes s1 and s2 currently point to the middle of a mismatching number,
  // where s1 < s2.
  char c1;
  char c2;
  do {
    c1 = *(++s1);
    c2 = *(++s2);
    if (is_not_digit(c1)) {
      return -1;
    }
  } while (is_digit(c2));
  return 1;
}

// We have the following major states:
//   0 (initial): strings perfectly match so far, last char (if any) is
//                nonnumeric.
//   1: strings perfectly match so far, last char is numeric.
//   2: strings match except for capitalization, last char is nonnumeric.
//   3: strings match except for capitalization, last char is numeric.
// strcmp_natural_tiebroken() expresses the logic for states 2 and 3, while
// strcmp_natural_uncasted() handles states 0 and 1.
int strcmp_natural_tiebroken(const char* s1, const char* s2) {
  // assumes ties should be broken in favor of s2.
  char c1 = *(++s1);
  char c2 = *(++s2);
  while (is_not_nzdigit(c1) && is_not_nzdigit(c2)) {
    // state 2
  strcmp_natural_tiebroken_state_2:
    if (c1 != c2) {
      if ((c1 >= 'a') && (c1 <= 'z')) {
	c1 -= 32;
      }
      if ((c2 >= 'a') && (c2 <= 'z')) {
	c2 -= 32;
      }
      if (c1 < c2) {
	return -1;
      } else if (c1 > c2) {
	return 1;
      }
    } else if (!c1) {
      return -1;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  }
  if (is_not_nzdigit(c1) || is_not_nzdigit(c2)) {
    return (c1 < c2)? -1 : 1;
  }
  do {
    // state 3
    if (c1 != c2) {
      if (is_digit(c2)) {
	if (c1 < c2) {
	  return strcmp_natural_scan_forward(s1, s2);
	} else {
	  return -strcmp_natural_scan_forward(s2, s1);
	}
      }
      return 1;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  } while (is_digit(c1));
  if (is_digit(c2)) {
    return -1;
  }
  // skip the while (is_not_digit...) check
  goto strcmp_natural_tiebroken_state_2;
}

static inline int strcmp_natural_uncasted(const char* s1, const char* s2) {
  char c1 = *s1;
  char c2 = *s2;
  while (is_not_nzdigit(c1) && is_not_nzdigit(c2)) {
    // state 0
  strcmp_natural_uncasted_state_0:
    if (c1 != c2) {
      if ((c1 >= 'a') && (c1 <= 'z')) {
	if (c2 + 32 == c1) {
	  return -strcmp_natural_tiebroken(s2, s1);
	} else if ((c2 < 'a') || (c2 > 'z')) {
	  c1 -= 32;
	}
      } else if ((c2 >= 'a') && (c2 <= 'z')) {
	c2 -= 32;
	if (c1 == c2) {
	  return strcmp_natural_tiebroken(s1, s2);
	}
      }
      return (c1 < c2)? -1 : 1;
    } else if (!c1) {
      return 0;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  }
  if (is_not_nzdigit(c1) || is_not_nzdigit(c2)) {
    return (c1 < c2)? -1 : 1;
  }
  do {
    // state 1
    if (c1 != c2) {
      if (is_digit(c2)) {
	if (c1 < c2) {
	  return strcmp_natural_scan_forward(s1, s2);
	} else {
	  return -strcmp_natural_scan_forward(s2, s1);
	}
      }
      return 1;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  } while (is_digit(c1));
  if (is_digit(c2)) {
    return -1;
  }
  goto strcmp_natural_uncasted_state_0;
}

int strcmp_natural(const void* s1, const void* s2) {
  return strcmp_natural_uncasted((char*)s1, (char*)s2);
}

int strcmp_deref(const void* s1, const void* s2) {
  return strcmp(*(char**)s1, *(char**)s2);
}

int strcmp_natural_deref(const void* s1, const void* s2) {
  return strcmp_natural_uncasted(*(char**)s1, *(char**)s2);
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

int triangle_divide(long long cur_prod, int modif) {
  // return smallest integer vv for which (vv * (vv + modif)) is no smaller
  // than cur_prod, and neither term in the product is negative.  (Note the
  // lack of a divide by two; cur_prod should also be double its "true" value
  // as a result.)
  long long vv;
  if (cur_prod == 0) {
    if (modif < 0) {
      return -modif;
    } else {
      return 0;
    }
  }
  vv = (long long)sqrt((double)cur_prod);
  while ((vv - 1) * (vv + modif - 1) >= cur_prod) {
    vv--;
  }
  while (vv * (vv + modif) < cur_prod) {
    vv++;
  }
  return vv;
}

void parallel_bounds(int ct, int start, int parallel_idx, int parallel_tot, int* bound_start_ptr, int* bound_end_ptr) {
  int modif = 1 - start * 2;
  long long ct_tot = (long long)ct * (ct + modif);
  *bound_start_ptr = triangle_divide((ct_tot * parallel_idx) / parallel_tot, modif);
  *bound_end_ptr = triangle_divide((ct_tot * (parallel_idx + 1)) / parallel_tot, modif);
}

// set align to 1 for no alignment
void triangle_fill(unsigned int* target_arr, int ct, int pieces, int parallel_idx, int parallel_tot, int start, int align) {
  long long ct_tr;
  long long cur_prod;
  int modif = 1 - start * 2;
  int cur_piece = 1;
  int lbound;
  int ubound;
  int ii;
  int align_m1;
  parallel_bounds(ct, start, parallel_idx, parallel_tot, &lbound, &ubound);
  // x(x+1)/2 is divisible by y iff (x % (2y)) is 0 or (2y - 1).
  align *= 2;
  align_m1 = align - 1;
  target_arr[0] = lbound;
  target_arr[pieces] = ubound;
  cur_prod = (long long)lbound * (lbound + modif);
  ct_tr = ((long long)ubound * (ubound + modif) - cur_prod) / pieces;
  while (cur_piece < pieces) {
    cur_prod += ct_tr;
    lbound = triangle_divide(cur_prod, modif);
    ii = (lbound - start) & align_m1;
    if ((ii) && (ii != align_m1)) {
      lbound = start + ((lbound - start) | align_m1);
    }
    // lack of this check caused a nasty bug earlier
    if (lbound > ct) {
      lbound = ct;
    }
    target_arr[cur_piece++] = lbound;
  }
}

int write_ids(char* outname, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, char* person_ids, unsigned int max_person_id_len) {
  FILE* outfile;
  unsigned int uii;
  if (fopen_checked(&outfile, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
    if (!is_set(indiv_exclude, uii) && (fprintf(outfile, "%s\n", &(person_ids[uii * max_person_id_len])) < 0)) {
      return RET_WRITE_FAIL;
    }
  }
  if (fclose(outfile)) {
    return RET_WRITE_FAIL;
  }
  return 0;
}

int distance_d_write_ids(char* outname, char* outname_end, int dist_calc_type, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, char* person_ids, unsigned int max_person_id_len) {
  int retval;
  if (dist_calc_type & DISTANCE_ALCT) {
    strcpy(outname_end, ".dist.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
  }
  if (dist_calc_type & DISTANCE_IBS) {
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
  }
  if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
    strcpy(outname_end, ".mdist.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
  }
  return 0;
}

int distance_req(int calculation_type) {
  return ((calculation_type & CALC_DISTANCE) || ((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) && (!(calculation_type & CALC_GENOME))) || ((!(calculation_type & CALC_LOAD_DISTANCES)) && ((calculation_type & CALC_GROUPDIST) || (calculation_type & CALC_REGRESS_DISTANCE))));
}

int double_cmp(const void* aa, const void* bb) {
  double cc = *((const double*)aa) - *((const double*)bb);
  if (cc > 0.0) {
    return 1;
  } else if (cc < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

int double_cmp_deref(const void* aa, const void* bb) {
  double cc = **((const double**)aa) - **((const double**)bb);
  if (cc > 0.0) {
    return 1;
  } else if (cc < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

// alas, qsort_r not available on some Linux distributions

// Normally use qsort_ext(), but this version is necessary before wkspace has
// been allocated.
void qsort_ext2(char* main_arr, int arr_length, int item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, int secondary_item_len, char* proxy_arr, int proxy_len) {
  int ii;
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
}

// This actually tends to be faster than just sorting an array of indices,
// because of memory locality issues.
int qsort_ext(char* main_arr, int arr_length, int item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, int secondary_item_len) {
  // main_arr = packed array of equal-length items to sort
  // arr_length = number of items
  // item_length = byte count of each main_arr item
  // comparator_deref = returns positive if *first > *second, 0 if equal,
  //                    negative if *first < *second.  Note the extra
  //                    dereference.
  // secondary_arr = packed array of fixed-length records associated with the
  //                 main_arr items, to be resorted in the same way.  (e.g.
  //                 if one is building an index, this could start as a sorted
  //                 0..(n-1) sequence of integers; then, post-sort, this would
  //                 be a lookup table for the original position of each
  //                 main_arr item.)
  // secondary_item_len = byte count of each secondary_arr item
  int proxy_len = secondary_item_len + sizeof(void*);
  unsigned char* wkspace_mark = wkspace_base;
  char* proxy_arr;
  if (!arr_length) {
    return 0;
  }
  if (proxy_len < item_length) {
    proxy_len = item_length;
  }
  if (wkspace_alloc_c_checked(&proxy_arr, arr_length * proxy_len)) {
    return -1;
  }
  qsort_ext2(main_arr, arr_length, item_length, comparator_deref, secondary_arr, secondary_item_len, proxy_arr, proxy_len);
  wkspace_reset(wkspace_mark);
  return 0;
}

int bsearch_str(char* id_buf, char* lptr, int max_id_len, int min_idx, int max_idx) {
  int mid_idx;
  int ii;
  if (max_idx < min_idx) {
    return -1;
  }
  mid_idx = (min_idx + max_idx) / 2;
  ii = strcmp(id_buf, &(lptr[mid_idx * max_id_len]));
  if (ii) {
    if (ii < 0) {
      return bsearch_str(id_buf, lptr, max_id_len, min_idx, mid_idx - 1);
    } else {
      return bsearch_str(id_buf, lptr, max_id_len, mid_idx + 1, max_idx);
    }
  } else {
    return mid_idx;
  }
}

int bsearch_str_natural(char* id_buf, char* lptr, int max_id_len, int min_idx, int max_idx) {
  int mid_idx;
  int ii;
  if (max_idx < min_idx) {
    return -1;
  }
  mid_idx = (min_idx + max_idx) / 2;
  ii = strcmp_natural(id_buf, &(lptr[mid_idx * max_id_len]));
  if (ii) {
    if (ii < 0) {
      return bsearch_str_natural(id_buf, lptr, max_id_len, min_idx, mid_idx - 1);
    } else {
      return bsearch_str_natural(id_buf, lptr, max_id_len, mid_idx + 1, max_idx);
    }
  } else {
    return mid_idx;
  }
}

void fill_idbuf_fam_indiv(char* idbuf, char* fam_indiv, char fillchar) {
  char* iend_ptr = item_endnn(fam_indiv);
  unsigned int slen = (iend_ptr - fam_indiv);
  unsigned int slen2;
  memcpy(idbuf, fam_indiv, slen);
  idbuf[slen] = fillchar;
  fam_indiv = skip_initial_spaces(iend_ptr);
  iend_ptr = item_endnn(fam_indiv);
  slen2 = (iend_ptr - fam_indiv);
  memcpy(&(idbuf[slen + 1]), fam_indiv, slen2);
  idbuf[slen + slen2 + 1] = '\0';
}

int bsearch_fam_indiv(char* id_buf, char* lptr, int max_id_len, int filter_line_ct, char* fam_id, char* indiv_id) {
  // id_buf = workspace
  // lptr = packed, sorted list of ID strings to search over
  // fam_id and indiv_id are considered terminated by any space/eoln character
  int ii;
  int jj;
  if (!filter_line_ct) {
    return -1;
  }
  ii = strlen_se(fam_id);
  jj = strlen_se(indiv_id);
  if (ii + jj + 2 > max_id_len) {
    return -1;
  }
  memcpy(id_buf, fam_id, ii);
  id_buf[ii] = '\t';
  memcpy(&(id_buf[ii + 1]), indiv_id, jj);
  id_buf[ii + jj + 1] = '\0';
  return bsearch_str(id_buf, lptr, max_id_len, 0, filter_line_ct - 1);
}

int distance_open(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, char* outname, char* outname_end, const char* varsuffix, const char* mode, int dist_calc_type, int parallel_idx, int parallel_tot) {
  if (dist_calc_type & DISTANCE_ALCT) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".dist%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".dist%s", varsuffix);
    }
    strcpy(tbuf, outname_end);
    if (fopen_checked(outfile_ptr, outname, mode)) {
      return 1;
    }
  }
  if (dist_calc_type & DISTANCE_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mibs%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".mibs%s", varsuffix);
    }
    strcpy(&(tbuf[MAX_POST_EXT]), outname_end);
    if (fopen_checked(outfile2_ptr, outname, mode)) {
      return 1;
    }
  }
  if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mdist%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".mdist%s", varsuffix);
    }
    strcpy(&(tbuf[MAX_POST_EXT * 2]), outname_end);
    if (fopen_checked(outfile3_ptr, outname, mode)) {
      return 1;
    }
  }
  return 0;
}

int distance_open_gz(gzFile* gz_outfile_ptr, gzFile* gz_outfile2_ptr, gzFile* gz_outfile3_ptr, char* outname, char* outname_end, int dist_calc_type, int parallel_idx, int parallel_tot) {
  if (dist_calc_type & DISTANCE_ALCT) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".dist.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".dist.gz");
    }
    strcpy(tbuf, outname_end);
    if (gzopen_checked(gz_outfile_ptr, outname, "wb")) {
      return 1;
    }
  }
  if (dist_calc_type & DISTANCE_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mibs.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".mibs.gz");
    }
    strcpy(&(tbuf[MAX_POST_EXT]), outname_end);
    if (gzopen_checked(gz_outfile_ptr, outname, "wb")) {
      return 1;
    }
  }
  if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mdist.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".mdist.gz");
    }
    strcpy(&(tbuf[MAX_POST_EXT * 2]), outname_end);
    if (gzopen_checked(gz_outfile_ptr, outname, "wb")) {
      return 1;
    }
  }
  return 0;
}

void distance_print_done(int format_code, char* outname, char* outname_end) {
  putchar('\r');
  if (!format_code) {
    strcpy(outname_end, tbuf);
    sprintf(logbuf, "Distances (allele counts) written to %s.\n", outname);
  } else if (format_code == 1) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT]));
    sprintf(logbuf, "IBS matrix written to %s.\n", outname);
  } else if (format_code == 2) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT * 2]));
    sprintf(logbuf, "Distances (proportions) written to %s.\n", outname);
  }
  logprintb();
}

#ifdef __LP64__
// Basic SSE2 implementation of Lauradoux/Walisch popcount.
static inline unsigned long popcount_vecs(__m128i* vptr, unsigned long ct) {
  // popcounts vptr[0..(ct-1)].  Assumes ct is a NONZERO multiple of 3.
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  unsigned long tot = 0;
  __m128i* vend;
  __m128i count1, count2, half1, half2;
  __uni16 acc;

  while (ct >= 30) {
    acc.vi = _mm_setzero_si128();
    vend = &(vptr[30]);
    do {
      count1 = *vptr++;
      count2 = *vptr++;
      half1 = *vptr++;
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1
      // count2 store a partial bitcount covering themselves AND another bit
      // from elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
    acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
    tot += (unsigned int)(acc.u8[0] + acc.u8[1]);
    ct -= 30;
  }
  acc.vi = _mm_setzero_si128();
  vend = &(vptr[ct]);
  do {
    count1 = *vptr++;
    count2 = *vptr++;
    half1 = *vptr++;
    half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
    half1 = _mm_and_si128(half1, m1);
    count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
    count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
    count1 = _mm_add_epi64(count1, half1);
    count2 = _mm_add_epi64(count2, half2);
    count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
    count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
  } while (vptr < vend);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
  acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
  tot += (unsigned int)(acc.u8[0] + acc.u8[1]);
  return tot;
}
#endif

unsigned long popcount_longs(unsigned long* lptr, unsigned long start_idx, unsigned long end_idx) {
  // given an aligned long array lptr[], this efficiently popcounts
  // lptr[start_idx..(end_idx - 1)].
  unsigned long tot = 0;
  unsigned long* lptr_end = &(lptr[end_idx]);
#ifdef __LP64__
  unsigned long six_ct;
  __m128i* vptr;
  if (start_idx == end_idx) {
    return 0;
  }
  if (start_idx & 1) {
    tot = popcount_long(lptr[start_idx++]);
  }
  vptr = (__m128i*)(&(lptr[start_idx]));
  six_ct = (end_idx - start_idx) / 6;
  if (six_ct) {
    tot += popcount_vecs(vptr, six_ct * 3);
  }
  lptr = &(lptr[start_idx + six_ct * 6]);
#else
  // The humble 16-bit lookup table actually beats
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
  // on my development machine by a hair.
  // However, if we take the hint from Walisch-Lauradoux and postpone the
  // multiply and right shift, this is no longer true.  Ah well.
  unsigned long* lptr_six_end;
  unsigned long tmp_stor;
  unsigned long loader;
  unsigned long ulii;
  unsigned long uljj;
  lptr = &(lptr[start_idx]);
  lptr_six_end = &(lptr[6 * ((end_idx - start_idx) / 6)]);
  while (lptr < lptr_six_end) {
    loader = *lptr++;
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = *lptr++;
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101 is
    // equivalent to the left-shifts and adds we need to sum those four 8-bit
    // numbers in the high-order slot.
    tot += (tmp_stor * 0x01010101) >> 24;
  }
#endif
  while (lptr < lptr_end) {
    tot += popcount_long(*lptr++);
  }
  return tot;
}

unsigned long popcount_chars(unsigned long* lptr, unsigned long start_idx, unsigned long end_idx) {
  // given a CHAR array c[] starting at the aligned address of lptr, this
  // efficiently popcounts c[start_idx..(end_idx - 1)].
  unsigned int extra_ct = (start_idx % (BITCT / 8));
  unsigned long si8l = start_idx / (BITCT / 8);
  unsigned long ei8l = end_idx / (BITCT / 8);
  unsigned long tot = 0;
  if (si8l == ei8l) {
    return popcount_long(lptr[si8l] & ((1LU << ((end_idx % (BITCT / 8)) * 8)) - (1LU << (extra_ct * 8))));
  } else {
    if (extra_ct) {
      tot = popcount_long(lptr[si8l++] >> (extra_ct * 8));
    }
    tot += popcount_longs(lptr, si8l, ei8l);
    extra_ct = end_idx % (BITCT / 8);
    if (extra_ct) {
      tot += popcount_long(lptr[ei8l] & ((1LU << (extra_ct * 8)) - 1LU));
    }
    return tot;
  }
}

int distance_d_write(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, gzFile* gz_outfile_ptr, gzFile* gz_outfile2_ptr, gzFile* gz_outfile3_ptr, int dist_calc_type, char* outname, char* outname_end, double* dists, double half_marker_ct_recip, unsigned int indiv_ct, int first_indiv_idx, int end_indiv_idx, int parallel_idx, int parallel_tot, unsigned char* membuf) {
  // membuf assumed to be of at least size indiv_ct * 8.
  int shape = dist_calc_type & DISTANCE_SHAPEMASK;
  int write_alcts = dist_calc_type & DISTANCE_ALCT;
  int write_ibs_matrix = dist_calc_type & DISTANCE_IBS;
  int write_1mibs_matrix = dist_calc_type & DISTANCE_1_MINUS_IBS;
  long long indiv_idx_offset = ((long long)first_indiv_idx * (first_indiv_idx - 1)) / 2;
  long long indiv_idx_ct = ((long long)end_indiv_idx * (end_indiv_idx - 1)) / 2 - indiv_idx_offset;
  double dxx;
  double dyy;
  double* dist_ptr;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long* glptr;
  unsigned int uii;
  unsigned int pct;
  unsigned int indiv_idx;
  int ii;
  int jj;
  char* cptr;
  if (first_indiv_idx == 1) {
    first_indiv_idx = 0;
  }
  if (shape == DISTANCE_SQ0) {
    cptr = (char*)(&ulii);
    for (uii = 0; uii < sizeof(long); uii += 2) {
      cptr[uii] = '\t';
      cptr[uii + 1] = '0';
    }
    ii = (indiv_ct * 2 + sizeof(long) - 2) / sizeof(long);
    glptr = (unsigned long*)membuf;
    for (jj = 0; jj < ii; jj++) {
      *glptr++ = ulii;
    }
  }
  pct = 1;
  if (dist_calc_type & DISTANCE_GZ) {
    if (distance_open_gz(gz_outfile_ptr, gz_outfile2_ptr, gz_outfile3_ptr, outname, outname_end, dist_calc_type, parallel_idx, parallel_tot)) {
      return RET_OPEN_FAIL;
    }
    if (first_indiv_idx) {
      ii = first_indiv_idx;
    } else {
      if (shape == DISTANCE_SQ0) {
	if (write_alcts) {
	  if (gzwrite_checked(*gz_outfile_ptr, &(membuf[1]), indiv_ct * 2 - 1)) {
	    return RET_WRITE_FAIL;
	  }
	  if (gzwrite_checked(*gz_outfile_ptr, "\n", 1)) {
	    return RET_WRITE_FAIL;
	  }
	}
	if (write_ibs_matrix) {
	  membuf[1] = '1';
	  if (gzwrite_checked(*gz_outfile2_ptr, &(membuf[1]), indiv_ct * 2 - 1)) {
	    return RET_WRITE_FAIL;
	  }
	  if (gzwrite_checked(*gz_outfile2_ptr, "\n", 1)) {
	    return RET_WRITE_FAIL;
	  }
	  membuf[1] = '0';
	}
	if (write_1mibs_matrix) {
	  if (gzwrite_checked(*gz_outfile3_ptr, &(membuf[1]), indiv_ct * 2 - 1)) {
	    return RET_WRITE_FAIL;
	  }
	  if (gzwrite_checked(*gz_outfile3_ptr, "\n", 1)) {
	    return RET_WRITE_FAIL;
	  }
	}
      } else if (shape == DISTANCE_SQ) {
	if (write_alcts) {
	  if (gzwrite_checked(*gz_outfile_ptr, "0", 1)) {
	    return RET_WRITE_FAIL;
	  }
	  for (indiv_idx = 1; indiv_idx < indiv_ct; indiv_idx++) {
	    if (!gzprintf(*gz_outfile_ptr, "\t%g", dists[((unsigned long)indiv_idx * (indiv_idx - 1)) / 2])) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  if (gzwrite_checked(*gz_outfile_ptr, "\n", 1)) {
	    return RET_WRITE_FAIL;
	  }
	}
	if (write_ibs_matrix) {
	  if (gzwrite_checked(*gz_outfile2_ptr, "1", 1)) {
	    return RET_WRITE_FAIL;
	  }
	  for (indiv_idx = 1; indiv_idx < indiv_ct; indiv_idx++) {
	    if (!gzprintf(*gz_outfile2_ptr, "\t%g", 1.0 - dists[((unsigned long)indiv_idx * (indiv_idx - 1)) / 2] * half_marker_ct_recip)) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  if (gzwrite_checked(*gz_outfile2_ptr, "\n", 1)) {
	    return RET_WRITE_FAIL;
	  }
	}
	if (write_1mibs_matrix) {
	  if (gzwrite_checked(*gz_outfile3_ptr, "0", 1)) {
	    return RET_WRITE_FAIL;
	  }
	  for (indiv_idx = 1; indiv_idx < indiv_ct; indiv_idx++) {
	    if (!gzprintf(*gz_outfile3_ptr, "\t%g", dists[((unsigned long)indiv_idx * (indiv_idx - 1)) / 2] * half_marker_ct_recip)) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  if (gzwrite_checked(*gz_outfile3_ptr, "\n", 1)) {
	    return RET_WRITE_FAIL;
	  }
	}
      }
      ii = 1;
    }
    if (write_alcts) {
      dist_ptr = dists;
      for (; ii < end_indiv_idx; ii++) {
	if (!gzprintf(*gz_outfile_ptr, "%g", *dist_ptr++)) {
	  return RET_WRITE_FAIL;
	}
	for (jj = 1; jj < ii; jj++) {
	  if (!gzprintf(*gz_outfile_ptr, "\t%g", *dist_ptr++)) {
	    return RET_WRITE_FAIL;
	  }
	}
	if (shape == DISTANCE_SQ0) {
	  if (gzwrite_checked(*gz_outfile_ptr, membuf, (indiv_ct - ii) * 2)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (long long)pct * (end_indiv_idx - first_indiv_idx)) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == DISTANCE_SQ) {
	    if (gzwrite_checked(*gz_outfile_ptr, "\t0", 2)) {
	      return RET_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      if (!gzprintf(*gz_outfile_ptr, "\t%g", dists[((ulii * (ulii - 1)) / 2) + ii])) {
		return RET_WRITE_FAIL;
	      }
	    }
	  }
	  if (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100 >= indiv_idx_ct * pct) {
	    pct = (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100) / indiv_idx_ct;
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	}
	if (gzwrite_checked(*gz_outfile_ptr, "\n", 1)) {
	  return RET_WRITE_FAIL;
	}
      }
      gzclose(*gz_outfile_ptr);
      distance_print_done(0, outname, outname_end);
      pct = 1;
      *gz_outfile_ptr = NULL;
    }
    if (first_indiv_idx) {
      ii = first_indiv_idx;
    } else {
      ii = 1;
    }
    if (write_ibs_matrix) {
      dist_ptr = dists;
      membuf[1] = '1';
      for (; ii < end_indiv_idx; ii++) {
	if (!gzprintf(*gz_outfile2_ptr, "%g", 1.0 - (*dist_ptr++) * half_marker_ct_recip)) {
	  return RET_WRITE_FAIL;
	}
	for (jj = 1; jj < ii; jj++) {
	  if (!gzprintf(*gz_outfile2_ptr, "\t%g", 1.0 - (*dist_ptr++) * half_marker_ct_recip)) {
	    return RET_WRITE_FAIL;
	  }
	}
	if (shape == DISTANCE_SQ0) {
	  if (gzwrite_checked(*gz_outfile2_ptr, membuf, (indiv_ct - ii) * 2)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (long long)pct * (end_indiv_idx - first_indiv_idx)) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == DISTANCE_SQ) {
	    if (gzwrite_checked(*gz_outfile2_ptr, "\t1", 2)) {
	      return RET_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      if (!gzprintf(*gz_outfile2_ptr, "\t%g", 1.0 - dists[((ulii * (ulii - 1)) / 2) + ii] * half_marker_ct_recip)) {
		return RET_WRITE_FAIL;
	      }
	    }
	  }
	  if (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100 >= indiv_idx_ct * pct) {
	    pct = (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100) / indiv_idx_ct;
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	}
	if (gzwrite_checked(*gz_outfile2_ptr, "\n", 1)) {
	  return RET_WRITE_FAIL;
	}
      }
      membuf[1] = '0';
      gzclose(*gz_outfile2_ptr);
      distance_print_done(1, outname, outname_end);
      pct = 1;
      *gz_outfile2_ptr = NULL;
    }
    if (first_indiv_idx) {
      ii = first_indiv_idx;
    } else {
      ii = 1;
    }
    if (write_1mibs_matrix) {
      dist_ptr = dists;
      for (; ii < end_indiv_idx; ii++) {
	if (!gzprintf(*gz_outfile3_ptr, "%g", (*dist_ptr++) * half_marker_ct_recip)) {
	  return RET_WRITE_FAIL;
	}
	for (jj = 1; jj < ii; jj++) {
	  if (!gzprintf(*gz_outfile3_ptr, "\t%g", (*dist_ptr++) * half_marker_ct_recip)) {
	    return RET_WRITE_FAIL;
	  }
	}
	if (shape == DISTANCE_SQ0) {
	  if (gzwrite_checked(*gz_outfile3_ptr, membuf, (indiv_ct - ii) * 2)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (long long)pct * (end_indiv_idx - first_indiv_idx)) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == DISTANCE_SQ) {
	    if (gzwrite_checked(*gz_outfile3_ptr, "\t0", 2)) {
	      return RET_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      if (!gzprintf(*gz_outfile3_ptr, "\t%g", dists[((ulii * (ulii - 1)) / 2) + ii] * half_marker_ct_recip)) {
		return RET_WRITE_FAIL;
	      }
	    }
	  }
	  if (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100 >= indiv_idx_ct * pct) {
	    pct = (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100) / indiv_idx_ct;
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	}
	if (gzwrite_checked(*gz_outfile3_ptr, "\n", 1)) {
	  return RET_WRITE_FAIL;
	}
      }
      gzclose(*gz_outfile3_ptr);
      distance_print_done(2, outname, outname_end);
      *gz_outfile3_ptr = NULL;
    }
  } else if (dist_calc_type & DISTANCE_BIN) {
    if (distance_open(outfile_ptr, outfile2_ptr, outfile3_ptr, outname, outname_end, ".bin", "wb", dist_calc_type, parallel_idx, parallel_tot)) {
      return RET_OPEN_FAIL;
    }
    if (shape == DISTANCE_TRI) {
      if (write_alcts) {
	fputs("Writing...", stdout);
	fflush(stdout);
	if (fwrite_checked(dists, indiv_idx_ct * sizeof(double), *outfile_ptr)) {
	  return RET_WRITE_FAIL;
	}
	distance_print_done(0, outname, outname_end);
      }
      if (write_ibs_matrix) {
	dist_ptr = dists;
	ulii = 0;
	do {
	  uljj = (indiv_idx_ct * pct) / 100L;
	  for (; ulii < uljj; ulii++) {
	    dxx = 1.0 - (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile2_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  printf("\rWriting... %d%%", pct++);
	  fflush(stdout);
	} while (pct <= 100);
	distance_print_done(1, outname, outname_end);
      }
      if (write_1mibs_matrix) {
	dist_ptr = dists;
	ulii = 0;
	do {
	  uljj = (indiv_idx_ct * pct) / 100L;
	  for (; ulii < uljj; ulii++) {
	    dxx = (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile3_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  printf("\rWriting... %d%%", pct++);
	  fflush(stdout);
	} while (pct <= 100);
	distance_print_done(2, outname, outname_end);
      }
    } else {
      if (shape == DISTANCE_SQ0) {
	fill_double_zero((double*)membuf, indiv_ct);
      }
      if (write_alcts) {
	dxx = 0.0;
	dist_ptr = dists;
	for (ii = first_indiv_idx; ii < end_indiv_idx; ii++) {
	  if (fwrite_checked(dist_ptr, ii * sizeof(double), *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  dist_ptr = &(dist_ptr[ii]);
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checked(membuf, (indiv_ct - ii) * sizeof(double), *outfile_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	  } else {
	    // square matrix, no need to handle parallel case
	    if (fwrite_checked(&dxx, sizeof(double), *outfile_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      if (fwrite_checked(&(dists[(ulii * (ulii - 1)) / 2 + ii]), sizeof(double), *outfile_ptr)) {
		return RET_WRITE_FAIL;
	      }
	    }
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (long long)pct * (end_indiv_idx - first_indiv_idx)) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile_ptr)) {
	  return RET_WRITE_FAIL;
	}
	distance_print_done(0, outname, outname_end);
	pct = 1;
      }
      if (write_ibs_matrix) {
	dist_ptr = dists;
	dyy = 1.0;
	membuf[1] = '1';
	for (ii = first_indiv_idx; ii < end_indiv_idx; ii++) {
	  for (jj = 0; jj < ii; jj++) {
	    dxx = 1.0 - (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile2_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checked(membuf, (indiv_ct - ii) * sizeof(double), *outfile2_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	  } else {
	    // square matrix
	    if (fwrite_checked(&dyy, sizeof(double), *outfile2_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      dxx = 1.0 - dists[(ulii * (ulii - 1)) / 2 + ii] * half_marker_ct_recip;
	      if (fwrite_checked(&dxx, sizeof(double), *outfile2_ptr)) {
		return RET_WRITE_FAIL;
	      }
	    }
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (long long)pct * (end_indiv_idx - first_indiv_idx)) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	}
	membuf[1] = '0';
	if (fclose_null(outfile2_ptr)) {
	  return RET_WRITE_FAIL;
	}
	distance_print_done(1, outname, outname_end);
	pct = 1;
      }
      if (write_1mibs_matrix) {
	dist_ptr = dists;
	dyy = 0.0;
	for (ii = first_indiv_idx; ii < end_indiv_idx; ii++) {
	  for (jj = 0; jj < ii; jj++) {
	    dxx = (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile3_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checked(membuf, (indiv_ct - ii) * sizeof(double), *outfile3_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	  } else {
	    // square matrix
	    if (fwrite_checked(&dyy, sizeof(double), *outfile3_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      dxx = dists[(ulii * (ulii - 1)) / 2 + ii] * half_marker_ct_recip;
	      if (fwrite_checked(&dxx, sizeof(double), *outfile3_ptr)) {
		return RET_WRITE_FAIL;
	      }
	    }
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (long long)pct * (end_indiv_idx - first_indiv_idx)) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile3_ptr)) {
	  return RET_WRITE_FAIL;
	}
	distance_print_done(2, outname, outname_end);
      }
    }
  } else {
    if (distance_open(outfile_ptr, outfile2_ptr, outfile3_ptr, outname, outname_end, "", "w", dist_calc_type, parallel_idx, parallel_tot)) {
      return RET_OPEN_FAIL;
    }
    if (write_alcts) {
      if (first_indiv_idx) {
	ii = first_indiv_idx;
      } else {
	if (shape == DISTANCE_SQ0) {
	  if (fwrite_checked(&(membuf[1]), indiv_ct * 2 - 1, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if (fwrite_checked("\n", 1, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	} else if (shape == DISTANCE_SQ) {
	  if (fwrite_checked("0", 1, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  for (ulii = 1; ulii < indiv_ct; ulii++) {
	    if (fprintf(*outfile_ptr, "\t%g", dists[(ulii * (ulii - 1)) / 2]) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  if (fwrite_checked("\n", 1, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	}
	ii = 1;
      }
      dist_ptr = dists;
      for (; ii < end_indiv_idx; ii++) {
	if (fprintf(*outfile_ptr, "%g", *dist_ptr++) < 0) {
	  return RET_WRITE_FAIL;
	}
	for (jj = 1; jj < ii; jj++) {
	  if (fprintf(*outfile_ptr, "\t%g", *dist_ptr++) < 0) {
	    return RET_WRITE_FAIL;
	  }
	}
	if (shape == DISTANCE_SQ0) {
	  if (fwrite_checked(membuf, (indiv_ct - ii) * 2, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= ((long long)pct * (end_indiv_idx - first_indiv_idx))) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == DISTANCE_SQ) {
	    if (fwrite_checked("\t0", 2, *outfile_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      if (fprintf(*outfile_ptr, "\t%g", dists[((ulii * (ulii - 1)) / 2) + ii]) < 0) {
		return RET_WRITE_FAIL;
	      }
	    }
	  }
	  if (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100 >= indiv_idx_ct * pct) {
	    pct = (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100) / indiv_idx_ct;
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	}
	if (fwrite_checked("\n", 1, *outfile_ptr)) {
	  return RET_WRITE_FAIL;
	}
      }
      if (fclose_null(outfile_ptr)) {
	return RET_WRITE_FAIL;
      }
      distance_print_done(0, outname, outname_end);
      pct = 1;
    }
    if (write_ibs_matrix) {
      membuf[1] = '1';
      if (first_indiv_idx) {
	ii = first_indiv_idx;
      } else {
	if (shape == DISTANCE_SQ0) {
	  if (fwrite_checked(&(membuf[1]), indiv_ct * 2 - 1, *outfile2_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if (fwrite_checked("\n", 1, *outfile2_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	} else if (shape == DISTANCE_SQ) {
	  if (fwrite_checked("1", 1, *outfile2_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  for (ulii = 1; ulii < indiv_ct; ulii++) {
	    if (fprintf(*outfile2_ptr, "\t%g", 1.0 - dists[(ulii * (ulii - 1)) / 2] * half_marker_ct_recip) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  if (fwrite_checked("\n", 1, *outfile2_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	} else {
	  if (fwrite_checked("1\n", 2, *outfile2_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	}
	ii = 1;
      }
      dist_ptr = dists;
      for (; ii < end_indiv_idx; ii++) {
	if (fprintf(*outfile2_ptr, "%g", 1.0 - (*dist_ptr++) * half_marker_ct_recip) < 0) {
	  return RET_WRITE_FAIL;
	}
	for (jj = 1; jj < ii; jj++) {
	  if (fprintf(*outfile2_ptr, "\t%g", 1.0 - (*dist_ptr++) * half_marker_ct_recip) < 0) {
	    return RET_WRITE_FAIL;
	  }
	}
	if (shape == DISTANCE_SQ0) {
	  if (fwrite_checked(membuf, (indiv_ct - ii) * 2, *outfile2_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= ((long long)pct * (end_indiv_idx - first_indiv_idx))) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (fwrite_checked("\t1", 2, *outfile2_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if (shape == DISTANCE_SQ) {
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      if (fprintf(*outfile2_ptr, "\t%g", 1.0 - dists[((ulii * (ulii - 1)) / 2) + ii] * half_marker_ct_recip) < 0) {
		return RET_WRITE_FAIL;
	      }
	    }
	  }
	  if (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100 >= indiv_idx_ct * pct) {
	    pct = (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100) / indiv_idx_ct;
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	}
	if (fwrite_checked("\n", 1, *outfile2_ptr)) {
	  return RET_WRITE_FAIL;
	}
      }
      membuf[1] = '0';
      if (fclose_null(outfile2_ptr)) {
	return RET_WRITE_FAIL;
      }
      distance_print_done(1, outname, outname_end);
      pct = 1;
    }
    if (write_1mibs_matrix) {
      if (first_indiv_idx) {
	ii = first_indiv_idx;
      } else {
	if (shape == DISTANCE_SQ0) {
	  if (fwrite_checked(&(membuf[1]), indiv_ct * 2 - 1, *outfile3_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if (fwrite_checked("\n", 1, *outfile3_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	} else if (shape == DISTANCE_SQ) {
	  if (fwrite_checked("0", 1, *outfile3_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  for (ulii = 1; ulii < indiv_ct; ulii++) {
	    if (fprintf(*outfile3_ptr, "\t%g", dists[(ulii * (ulii - 1)) / 2] * half_marker_ct_recip) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  if (fwrite_checked("\n", 1, *outfile3_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	}
	ii = 1;
      }
      dist_ptr = dists;
      for (; ii < end_indiv_idx; ii++) {
	if (fprintf(*outfile3_ptr, "%g", (*dist_ptr++) * half_marker_ct_recip) < 0) {
	  return RET_WRITE_FAIL;
	}
	for (jj = 1; jj < ii; jj++) {
	  if (fprintf(*outfile3_ptr, "\t%g", (*dist_ptr++) * half_marker_ct_recip) < 0) {
	    return RET_WRITE_FAIL;
	  }
	}
	if (shape == DISTANCE_SQ0) {
	  if (fwrite_checked(membuf, (indiv_ct - ii) * 2, *outfile3_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= ((long long)pct * (end_indiv_idx - first_indiv_idx))) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == DISTANCE_SQ) {
	    if (fwrite_checked("\t0", 2, *outfile3_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      if (fprintf(*outfile3_ptr, "\t%g", dists[((ulii * (ulii - 1)) / 2) + ii] * half_marker_ct_recip) < 0) {
		return RET_WRITE_FAIL;
	      }
	    }
	  }
	  if (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100 >= indiv_idx_ct * pct) {
	    pct = (((long long)ii * (ii + 1) / 2 - indiv_idx_offset) * 100) / indiv_idx_ct;
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	}
	if (fwrite_checked("\n", 1, *outfile3_ptr)) {
	  return RET_WRITE_FAIL;
	}
      }
      if (fclose_null(outfile3_ptr)) {
	return RET_WRITE_FAIL;
      }
      distance_print_done(2, outname, outname_end);
    }
  }
  return 0;
}

void collapse_arr(char* item_arr, int fixed_item_len, unsigned long* exclude_arr, int exclude_arr_size) {
  // collapses array of fixed-length items
  int ii = 0;
  int jj;
  while ((ii < exclude_arr_size) && (!is_set(exclude_arr, ii))) {
    ii++;
  }
  jj = ii;
  while (++ii < exclude_arr_size) {
    if (!is_set(exclude_arr, ii)) {
      memcpy(&(item_arr[(jj++) * fixed_item_len]), &(item_arr[ii * fixed_item_len]), fixed_item_len);
    }
  }
}

// A C-program for MT19937, with initialization improved 2002/1/26.
// Coded by Takuji Nishimura and Makoto Matsumoto.

// Before using, initialize the state by using init_genrand(seed)  
// or init_by_array(init_key, key_length).

// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// All rights reserved.                          

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:

//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.

//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.

//   3. The names of its contributors may not be used to endorse or promote 
//      products derived from this software without specific prior written 
//      permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


// Any feedback is very welcome.
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
// email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)

// Period parameters
#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL   // constant vector a
#define UPPER_MASK 0x80000000UL // most significant w-r bits
#define LOWER_MASK 0x7fffffffUL // least significant r bits

static unsigned long mt[MT_N]; // the array for the state vector
static int mti=MT_N+1; // mti==N+1 means mt[N] is not initialized

// initializes mt[MT_N] with a seed
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<MT_N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
        // In the previous versions, MSBs of the seed affect
        // only MSBs of the array mt[].
        // 2002/01/09 modified by Makoto Matsumoto
        mt[mti] &= 0xffffffffUL;
        // for >32 bit machines
    }
}

// see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
// if you want to add init_by_array back

// generates a random number on [0,0xffffffff]-interval
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    // mag01[x] = x * MATRIX_A  for x=0,1

    if (mti >= MT_N) { // generate MT_N words at one time
        int kk;

        if (mti == MT_N+1)   // if init_genrand() has not been called,
	  init_genrand(5489UL); // a default initial seed is used

        for (kk=0;kk<MT_N-MT_M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<MT_N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[MT_N-1] = mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

double rand_unif(void) {
  return (genrand_int32() + 0.5) * RECIP_2_32;
}

double rand_normal(void) {
  // N(0, 1)
  return sqrt(-2 * log(rand_unif())) * cos(2 * PI * rand_unif());
}

void pick_d(unsigned char* cbuf, unsigned int ct, unsigned int dd) {
  unsigned int ii;
  unsigned int jj;
  unsigned int kk;
  memset(cbuf, 0, ct);
  kk = 1073741824 % ct;
  kk = (kk * 4) % ct;
  for (ii = 0; ii < dd; ii++) {
    do {
      do {
        jj = genrand_int32();
      } while (jj < kk);
      jj %= ct;
    } while (cbuf[jj]);
    cbuf[jj] = 1;
  }
}

void pick_d_small(unsigned char* tmp_cbuf, int* ibuf, unsigned int ct, unsigned int dd) {
  unsigned int uii;
  pick_d(tmp_cbuf, ct, dd);
  for (uii = 0; uii < ct; uii++) {
    if (tmp_cbuf[uii]) {
      *ibuf++ = uii;
    }
  }
  *ibuf = ct;
}

void print_pheno_stdev(double* pheno_d, unsigned int indiv_ct) {
  double reg_tot_x = 0.0;
  double reg_tot_xx = 0.0;
  double dxx;
  unsigned int uii;
  for (uii = 0; uii < indiv_ct; uii++) {
    dxx = pheno_d[uii];
    reg_tot_x += dxx;
    reg_tot_xx += dxx * dxx;
  }
  sprintf(logbuf, "Phenotype stdev: %g\n", sqrt((reg_tot_xx - reg_tot_x * reg_tot_x / indiv_ct) / (indiv_ct - 1)));
  logprintb();
}

unsigned int set_default_jackknife_d(unsigned int ct) {
  unsigned int dd = (unsigned int)pow((double)ct, 0.600000000001);
  sprintf(logbuf, "Setting d=%u for jackknife.\n", dd);
  logprintb();
  return dd;
}

// ----- multithread globals -----
static double* g_pheno_d;
static unsigned long g_jackknife_iters;
static unsigned int g_jackknife_d;
static double g_reg_tot_xy;
static double g_reg_tot_x;
static double g_reg_tot_y;
static double g_reg_tot_xx;
static double g_reg_tot_yy;
static double* g_jackknife_precomp;
static double* g_dists;
static double g_calc_result[4][MAX_THREADS_P1];
static unsigned char* g_generic_buf;
static unsigned int g_indiv_ct;

// double regress_jack(int* ibuf) {
double regress_jack(int* ibuf, double* ret2_ptr) {
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
    dptr2 = &(g_jackknife_precomp[(*iptr++) * JACKKNIFE_VALS_DIST]);
    neg_tot_xy += *dptr2++;
    neg_tot_x += *dptr2++;
    neg_tot_y += *dptr2++;
    neg_tot_xx += *dptr2++;
    neg_tot_yy += *dptr2++;
  }
  iptr = ibuf;
  for (uii = 1; uii < g_jackknife_d; uii++) {
    jj = *(++iptr);
    dxx1 = g_pheno_d[jj];
    jptr = ibuf;
    dptr = &(g_dists[((long)jj * (jj - 1)) / 2]);
    while (jptr < iptr) {
      kk = *jptr++;
      dxx = (dxx1 + g_pheno_d[kk]) * 0.5;
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

void* regress_jack_thread(void* arg) {
  long tidx = (long)arg;
  int* ibuf = (int*)(&(g_generic_buf[tidx * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int))]));
  unsigned char* cbuf = &(g_generic_buf[tidx * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int)) + (g_jackknife_d + 1) * sizeof(int)]);
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
    dxx = regress_jack(ibuf, &ret2);
    // dxx = regress_jack(ibuf);
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

int regress_distance(int calculation_type, double* dists_local, double* pheno_d_local, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int indiv_ct_local, unsigned int thread_ct, unsigned long regress_iters, unsigned int regress_d) {
  unsigned char* wkspace_mark = wkspace_base;
  unsigned long ulii;
  unsigned int uii;
  double* dist_ptr;
  double* dptr2;
  double* dptr3;
  double* dptr4;
  double* dptr5;
  double dxx;
  double dyy;
  double dzz;
  double dww;
  double dvv;
  double duu;
  pthread_t threads[MAX_THREADS];

  g_dists = dists_local;
  g_pheno_d = pheno_d_local;
  g_indiv_ct = indiv_ct_local;

  // beta = (mean(xy) - mean(x)*mean(y)) / (mean(x^2) - mean(x)^2)
  if (unfiltered_indiv_ct != g_indiv_ct) {
    collapse_arr((char*)g_pheno_d, sizeof(double), indiv_exclude, unfiltered_indiv_ct);
  }
  if (!(calculation_type & CALC_REGRESS_REL)) {
    print_pheno_stdev(g_pheno_d, g_indiv_ct);
  }
  ulii = g_indiv_ct;
  ulii = ulii * (ulii - 1) / 2;
  g_reg_tot_xy = 0.0;
  g_reg_tot_x = 0.0;
  g_reg_tot_y = 0.0;
  g_reg_tot_xx = 0.0;
  g_reg_tot_yy = 0.0;
  dptr4 = g_dists;
  dist_ptr = g_pheno_d;
  // Linear regression slope is a function of sum(xy), sum(x), sum(y),
  // sum(x^2), and n.  To speed up the jackknife calculation, we precompute
  // (i) the global xy, x, y, x^2, and
  // (ii) the xy, x, y, x^2 for each row.
  // Then for each delete-d jackknife iteration, we take the global sums,
  // subtract the partial row sums corresponding to the deleted individuals,
  // and then add back the elements in the intersection of two deletions.
  g_jackknife_precomp = (double*)wkspace_alloc(g_indiv_ct * JACKKNIFE_VALS_DIST * sizeof(double));
  if (!g_jackknife_precomp) {
    return RET_NOMEM;
  }
  fill_double_zero(g_jackknife_precomp, g_indiv_ct * JACKKNIFE_VALS_DIST);
  for (uii = 1; uii < g_indiv_ct; uii++) {
    dzz = *(++dist_ptr);
    dptr2 = g_pheno_d;
    dptr3 = &(g_jackknife_precomp[uii * JACKKNIFE_VALS_DIST]);
    dptr5 = g_jackknife_precomp;
    while (dptr2 < dist_ptr) {
      dxx = (dzz + *dptr2++) * 0.5;
      dyy = (*dptr4++);
      dww = dxx * dyy;
      dvv = dxx * dxx;
      duu = dyy * dyy;
      g_reg_tot_xy += dww;
      *dptr3 += dww;
      *dptr5 += dww;
      dptr5++;
      g_reg_tot_x += dxx;
      dptr3[1] += dxx;
      *dptr5 += dxx;
      dptr5++;
      g_reg_tot_y += dyy;
      dptr3[2] += dyy;
      *dptr5 += dyy;
      dptr5++;
      g_reg_tot_xx += dvv;
      dptr3[3] += dvv;
      *dptr5 += dvv;
      dptr5++;
      g_reg_tot_yy += duu;
      dptr3[4] += duu;
      *dptr5 += duu;
      dptr5++;
    }
  }

  dxx = ulii;
  sprintf(logbuf, "Regression slope (y = genomic distance, x = avg phenotype): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y / dxx) / (g_reg_tot_xx - g_reg_tot_x * g_reg_tot_x / dxx));
  logprintb();
  sprintf(logbuf, "Regression slope (y = avg phenotype, x = genomic distance): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y / dxx) / (g_reg_tot_yy - g_reg_tot_y * g_reg_tot_y / dxx));
  logprintb();

  g_jackknife_iters = (regress_iters + thread_ct - 1) / thread_ct;
  if (regress_d) {
    g_jackknife_d = regress_d;
  } else {
    g_jackknife_d = set_default_jackknife_d(g_indiv_ct);
  }
  g_generic_buf = wkspace_alloc(thread_ct * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int)));
  if (!g_generic_buf) {
    return RET_NOMEM;
  }
  for (ulii = 1; ulii < thread_ct; ulii++) {
    if (pthread_create(&(threads[ulii - 1]), NULL, &regress_jack_thread, (void*)ulii)) {
      fputs(errstr_thread_create, stdout);
      while (--ulii) {
	pthread_join(threads[ulii - 1], NULL);
      }
      return RET_THREAD_CREATE_FAIL;
    }
  }
  ulii = 0;
  regress_jack_thread((void*)ulii);
  dyy = g_calc_result[0][0]; // sum
  dzz = g_calc_result[1][0]; // sum of squares
  dww = g_calc_result[2][0]; // reverse regression sum
  dvv = g_calc_result[3][0]; // reverse regression sum of squares
  for (uii = 0; uii < thread_ct - 1; uii++) {
    pthread_join(threads[uii], NULL);
    dyy += g_calc_result[0][uii + 1];
    dzz += g_calc_result[1][uii + 1];
    dww += g_calc_result[2][uii + 1];
    dvv += g_calc_result[3][uii + 1];
  }
  regress_iters = g_jackknife_iters * thread_ct;
  putchar('\r');
  sprintf(logbuf, "Jackknife s.e.: %g\n", sqrt((g_indiv_ct / ((double)g_jackknife_d)) * (dzz - dyy * dyy / regress_iters) / (regress_iters - 1)));
  logprintb();
  sprintf(logbuf, "Jackknife s.e. (y = avg phenotype): %g\n", sqrt((g_indiv_ct / ((double)g_jackknife_d)) * (dvv - dww * dww / regress_iters) / (regress_iters - 1)));
  logprintb();
  wkspace_reset(wkspace_mark);
  return 0;
}
