#include "wdist_common.h"

const char errstr_fopen[] = "Error: Failed to open %s.\n";
const char errstr_append[] = "\nFor more information, try 'wdist --help {flag names}' or 'wdist --help | more'.\n";
const char errstr_thread_create[] = "\nError: Failed to create thread.\n";

char tbuf[MAXLINELEN * 4 + 256];

int fopen_checked(FILE** target_ptr, const char* fname, const char* mode) {
  *target_ptr = fopen(fname, mode);
  if (!(*target_ptr)) {
    printf(errstr_fopen, fname);
    return -1;
  }
  return 0;
}

int gzopen_checked(gzFile* target_ptr, const char* fname, const char* mode) {
  *target_ptr = gzopen(fname, mode);
  if (!(*target_ptr)) {
    printf(errstr_fopen, fname);
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
  if (!sptr) {
    return NULL;
  }
  while ((*sptr != ' ') && (*sptr != '\t')) {
    if (!(*sptr)) {
      return NULL;
    }
    sptr++;
  }
  return sptr;
}

char* item_endl(char* sptr) {
  if (!sptr) {
    return NULL;
  }
  while ((*sptr != ' ') && (*sptr != '\t') && (*sptr)) {
    sptr++;
  }
  return sptr;
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

int distance_d_write_ids(char* outname, char* outname_end, int calculation_type, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, char* person_ids, unsigned int max_person_id_len) {
  int retval;
  if (calculation_type & CALC_DISTANCE_ALCT) {
    strcpy(outname_end, ".dist.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
  }
  if (calculation_type & CALC_DISTANCE_IBS) {
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
  }
  if (calculation_type & CALC_DISTANCE_1_MINUS_IBS) {
    strcpy(outname_end, ".mdist.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
  }
  return 0;
}

int distance_req(int calculation_type) {
  return ((calculation_type & CALC_DISTANCE_MASK) || ((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) && (!(calculation_type & CALC_GENOME))) || ((!(calculation_type & CALC_LOAD_DISTANCES)) && ((calculation_type & CALC_GROUPDIST) || (calculation_type & CALC_REGRESS_DISTANCE))));
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

int distance_open(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, char* outname, char* outname_end, const char* varsuffix, const char* mode, int calculation_type, int parallel_idx, int parallel_tot) {
  if (calculation_type & CALC_DISTANCE_ALCT) {
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
  if (calculation_type & CALC_DISTANCE_IBS) {
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
  if (calculation_type & CALC_DISTANCE_1_MINUS_IBS) {
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

int distance_open_gz(gzFile* gz_outfile_ptr, gzFile* gz_outfile2_ptr, gzFile* gz_outfile3_ptr, char* outname, char* outname_end, int calculation_type, int parallel_idx, int parallel_tot) {
  if (calculation_type & CALC_DISTANCE_ALCT) {
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
  if (calculation_type & CALC_DISTANCE_IBS) {
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
  if (calculation_type & CALC_DISTANCE_1_MINUS_IBS) {
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
  if (!format_code) {
    strcpy(outname_end, tbuf);
    printf("\rDistances (allele counts) written to %s.\n", outname);
  } else if (format_code == 1) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT]));
    printf("\rIBS matrix written to %s.\n", outname);
  } else if (format_code == 2) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT * 2]));
    printf("\rDistances (proportions) written to %s.\n", outname);
  }
}

int distance_d_write(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, gzFile* gz_outfile_ptr, gzFile* gz_outfile2_ptr, gzFile* gz_outfile3_ptr, int calculation_type, char* outname, char* outname_end, double* dists, double half_marker_ct_recip, unsigned int indiv_ct, int first_indiv_idx, int end_indiv_idx, int parallel_idx, int parallel_tot, unsigned char* membuf) {
  // membuf assumed to be of at least size indiv_ct * 8.
  int shape = calculation_type & CALC_DISTANCE_SHAPEMASK;
  int write_alcts = calculation_type & CALC_DISTANCE_ALCT;
  int write_ibs_matrix = calculation_type & CALC_DISTANCE_IBS;
  int write_1mibs_matrix = calculation_type & CALC_DISTANCE_1_MINUS_IBS;
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
  if (shape == CALC_DISTANCE_SQ0) {
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
  if (calculation_type & CALC_DISTANCE_GZ) {
    if (distance_open_gz(gz_outfile_ptr, gz_outfile2_ptr, gz_outfile3_ptr, outname, outname_end, calculation_type, parallel_idx, parallel_tot)) {
      return RET_OPEN_FAIL;
    }
    if (first_indiv_idx) {
      ii = first_indiv_idx;
    } else {
      if (shape == CALC_DISTANCE_SQ0) {
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
      } else if (shape == CALC_DISTANCE_SQ) {
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
	if (shape == CALC_DISTANCE_SQ0) {
	  if (gzwrite_checked(*gz_outfile_ptr, membuf, (indiv_ct - ii) * 2)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (long long)pct * (end_indiv_idx - first_indiv_idx)) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == CALC_DISTANCE_SQ) {
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
	if (shape == CALC_DISTANCE_SQ0) {
	  if (gzwrite_checked(*gz_outfile2_ptr, membuf, (indiv_ct - ii) * 2)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (long long)pct * (end_indiv_idx - first_indiv_idx)) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == CALC_DISTANCE_SQ) {
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
	if (shape == CALC_DISTANCE_SQ0) {
	  if (gzwrite_checked(*gz_outfile3_ptr, membuf, (indiv_ct - ii) * 2)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (long long)pct * (end_indiv_idx - first_indiv_idx)) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == CALC_DISTANCE_SQ) {
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
  } else if (calculation_type & CALC_DISTANCE_BIN) {
    if (distance_open(outfile_ptr, outfile2_ptr, outfile3_ptr, outname, outname_end, ".bin", "wb", calculation_type, parallel_idx, parallel_tot)) {
      return RET_OPEN_FAIL;
    }
    if (shape == CALC_DISTANCE_TRI) {
      if (write_alcts) {
	printf("Writing...");
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
      if (shape == CALC_DISTANCE_SQ0) {
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
	  if (shape == CALC_DISTANCE_SQ0) {
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
	  if (shape == CALC_DISTANCE_SQ0) {
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
	  if (shape == CALC_DISTANCE_SQ0) {
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
    if (distance_open(outfile_ptr, outfile2_ptr, outfile3_ptr, outname, outname_end, "", "w", calculation_type, parallel_idx, parallel_tot)) {
      return RET_OPEN_FAIL;
    }
    if (write_alcts) {
      if (first_indiv_idx) {
	ii = first_indiv_idx;
      } else {
	if (shape == CALC_DISTANCE_SQ0) {
	  if (fwrite_checked(&(membuf[1]), indiv_ct * 2 - 1, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if (fwrite_checked("\n", 1, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	} else if (shape == CALC_DISTANCE_SQ) {
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
	if (shape == CALC_DISTANCE_SQ0) {
	  if (fwrite_checked(membuf, (indiv_ct - ii) * 2, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= ((long long)pct * (end_indiv_idx - first_indiv_idx))) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == CALC_DISTANCE_SQ) {
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
	if (shape == CALC_DISTANCE_SQ0) {
	  if (fwrite_checked(&(membuf[1]), indiv_ct * 2 - 1, *outfile2_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if (fwrite_checked("\n", 1, *outfile2_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	} else if (shape == CALC_DISTANCE_SQ) {
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
	if (shape == CALC_DISTANCE_SQ0) {
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
	  if (shape == CALC_DISTANCE_SQ) {
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
	if (shape == CALC_DISTANCE_SQ0) {
	  if (fwrite_checked(&(membuf[1]), indiv_ct * 2 - 1, *outfile3_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if (fwrite_checked("\n", 1, *outfile3_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	} else if (shape == CALC_DISTANCE_SQ) {
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
	if (shape == CALC_DISTANCE_SQ0) {
	  if (fwrite_checked(membuf, (indiv_ct - ii) * 2, *outfile3_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if ((ii - first_indiv_idx) * 100LL >= ((long long)pct * (end_indiv_idx - first_indiv_idx))) {
	    pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (shape == CALC_DISTANCE_SQ) {
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
  printf("Phenotype stdev: %g\n", sqrt((reg_tot_xx - reg_tot_x * reg_tot_x / indiv_ct) / (indiv_ct - 1)));
}

unsigned int set_default_jackknife_d(unsigned int ct) {
  unsigned int dd = (unsigned int)pow((double)ct, 0.600000000001);
  printf("Setting d=%u for jackknife.\n", dd);
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
  printf("Regression slope (y = genomic distance, x = avg phenotype): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y / dxx) / (g_reg_tot_xx - g_reg_tot_x * g_reg_tot_x / dxx));
  printf("Regression slope (y = avg phenotype, x = genomic distance): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y / dxx) / (g_reg_tot_yy - g_reg_tot_y * g_reg_tot_y / dxx));

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
      printf(errstr_thread_create);
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
  printf("\rJackknife s.e.: %g\n", sqrt((g_indiv_ct / ((double)g_jackknife_d)) * (dzz - dyy * dyy / regress_iters) / (regress_iters - 1)));
  printf("Jackknife s.e. (y = avg phenotype): %g\n", sqrt((g_indiv_ct / ((double)g_jackknife_d)) * (dvv - dww * dww / regress_iters) / (regress_iters - 1)));
  wkspace_reset(wkspace_mark);
  return 0;
}
