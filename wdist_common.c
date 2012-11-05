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

int distance_req(int calculation_type) {
  return ((calculation_type & CALC_DISTANCE_MASK) || ((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) && (!(calculation_type & CALC_GENOME))) || ((!(calculation_type & CALC_LOAD_DISTANCES)) && ((calculation_type & CALC_GROUPDIST) || (calculation_type & CALC_REGRESS_DISTANCE))));
}

#ifndef __cplusplus
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
#endif

int distance_open(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, char* outname, char* outname_end, const char* varsuffix, const char* mode, int calculation_type, int parallel_idx, int parallel_tot) {
  if (calculation_type & CALC_DISTANCE_SNPS) {
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
  if (calculation_type & CALC_DISTANCE_SNPS) {
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
    printf("\rDistances (in SNPs) written to %s.\n", outname);
  } else if (format_code == 1) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT]));
    printf("\rIBS matrix written to %s.\n", outname);
  } else if (format_code == 2) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT * 2]));
    printf("\rDistances (proportions) written to %s.\n", outname);
  }
}

int distance_d_write(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, gzFile* gz_outfile_ptr, gzFile* gz_outfile2_ptr, gzFile* gz_outfile3_ptr, int calculation_type, char* outname, char* outname_end, double* dists, unsigned int marker_ct, unsigned int indiv_ct, int first_indiv_idx, int end_indiv_idx, int parallel_idx, int parallel_tot, unsigned char* membuf) {
  // membuf assumed to be of at least size indiv_ct * 8.
  int shape = calculation_type & CALC_DISTANCE_SHAPEMASK;
  int write_snp_cts = calculation_type & CALC_DISTANCE_SNPS;
  int write_ibs_matrix = calculation_type & CALC_DISTANCE_IBS;
  int write_1mibs_matrix = calculation_type & CALC_DISTANCE_1_MINUS_IBS;
  double half_marker_ct_recip = 0.5 / (double)marker_ct;
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
	if (write_snp_cts) {
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
	if (write_snp_cts) {
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
    if (write_snp_cts) {
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
      if (write_snp_cts) {
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
      if (write_snp_cts) {
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
    if (write_snp_cts) {
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
