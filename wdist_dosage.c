#include "wdist_common.h"
#include "wdist_dosage.h"

// Routines that handle dosage data instead of just 0-1-2 reference allele
// counts.  Only Oxford-formatted data is currently supported, but a PLINK
// --dosage loader will probably be added later.

int oxford_sample_load(char* samplename, unsigned int* unfiltered_indiv_ct_ptr, char** person_ids_ptr, unsigned int* max_person_id_len_ptr, double** phenos_ptr, unsigned long** pheno_exclude_ptr, unsigned long** indiv_exclude_ptr, char* missing_code) {
  FILE* samplefile = NULL;
  unsigned char* wkspace_mark = NULL;
  unsigned int unfiltered_indiv_ct = 0;
  unsigned int max_person_id_len = 4;
  unsigned int missing_code_ct = 0;
  unsigned int indiv_uidx = 0;
  char** missing_code_ptrs = NULL;
  unsigned int* missing_code_lens = NULL;
  unsigned int unfiltered_indiv_ctl;
  char* person_ids;
  double* phenos;
  char* item_begin;
  char* bufptr;
  unsigned int cur_person_id_len;
  unsigned int uii;
  unsigned int ujj;
  unsigned long long first_real_line_loc;
  int retval;
  int is_missing;
  double dxx;
  if (fopen_checked(&samplefile, samplename, "r")) {
    return RET_OPEN_FAIL;
  }
  // pass #1: just count number of samples
  tbuf[MAXLINELEN - 1] = ' ';
  if (!fgets(tbuf, MAXLINELEN, samplefile)) {
    if (feof(samplefile)) {
      goto oxford_sample_load_ret_INVALID_FORMAT;
    } else {
      goto oxford_sample_load_ret_READ_FAIL;
    }
  }
  if (memcmp(tbuf, "ID_1 ID_2 missing ", 18)) {
    goto oxford_sample_load_ret_INVALID_FORMAT;
  }
  if (!tbuf[MAXLINELEN - 1]) {
    goto oxford_sample_load_ret_INVALID_FORMAT_2;
  }
  if (!fgets(tbuf, MAXLINELEN, samplefile)) {
    if (feof(samplefile)) {
      goto oxford_sample_load_ret_INVALID_FORMAT;
    } else {
      goto oxford_sample_load_ret_READ_FAIL;
    }
  }
  if (memcmp(tbuf, "0 0 0 ", 6)) {
    goto oxford_sample_load_ret_INVALID_FORMAT;
  }
  if (!tbuf[MAXLINELEN - 1]) {
    goto oxford_sample_load_ret_INVALID_FORMAT_2;
  }
  first_real_line_loc = ftello(samplefile);
  unfiltered_indiv_ct = 0;
  while (fgets(tbuf, MAXLINELEN, samplefile) != NULL) {
    if (*tbuf == '\n') {
      continue;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      goto oxford_sample_load_ret_INVALID_FORMAT_2;
    }
    item_begin = skip_initial_spaces(tbuf);
    bufptr = item_end(item_begin);
    if (!bufptr) {
      goto oxford_sample_load_ret_INVALID_FORMAT;
    }
    cur_person_id_len = 2 + (unsigned int)(bufptr - item_begin);
    item_begin = skip_initial_spaces(bufptr);
    bufptr = item_end(item_begin);
    if (!bufptr) {
      goto oxford_sample_load_ret_INVALID_FORMAT;
    }
    cur_person_id_len += (unsigned int)(bufptr - item_begin);
    if (cur_person_id_len > max_person_id_len) {
      max_person_id_len = cur_person_id_len;
    }
    unfiltered_indiv_ct++;
  }
  if (!feof(samplefile)) {
    goto oxford_sample_load_ret_READ_FAIL;
  }
  if (!unfiltered_indiv_ct) {
    printf("Error: No individuals in .sample file.\n");
    goto oxford_sample_load_ret_INVALID_FORMAT_3;
  }
  *unfiltered_indiv_ct_ptr = unfiltered_indiv_ct;
  *max_person_id_len_ptr = max_person_id_len;

  if (wkspace_alloc_c_checked(person_ids_ptr, unfiltered_indiv_ct * max_person_id_len)) {
    goto oxford_sample_load_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(phenos_ptr, unfiltered_indiv_ct * sizeof(double))) {
    goto oxford_sample_load_ret_NOMEM;
  }
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(pheno_exclude_ptr, unfiltered_indiv_ctl * sizeof(long))) {
    goto oxford_sample_load_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(indiv_exclude_ptr, unfiltered_indiv_ctl * sizeof(long))) {
    goto oxford_sample_load_ret_NOMEM;
  }
  person_ids = *person_ids_ptr;
  phenos = *phenos_ptr;
  fill_ulong_zero(*pheno_exclude_ptr, unfiltered_indiv_ctl);
  fill_ulong_zero(*indiv_exclude_ptr, unfiltered_indiv_ctl);
  wkspace_mark = wkspace_base;
  if (*missing_code) {
    bufptr = missing_code;
    do {
      if ((*bufptr == ',') || (*bufptr == '\0')) {
	// blank string makes no sense
        printf("Error: Invalid --missing-code parameter '%s'.%s", missing_code, errstr_append);
	goto oxford_sample_load_ret_INVALID_CMDLINE;
      }
      bufptr = strchr(bufptr, ',');
      missing_code_ct++;
      if (bufptr) {
	bufptr++;
      }
    } while (bufptr);
  }
  if (missing_code_ct) {
    missing_code_ptrs = (char**)wkspace_alloc(missing_code_ct * sizeof(char*));    if (!missing_code_ptrs) {
      goto oxford_sample_load_ret_NOMEM;
    }
    if (wkspace_alloc_ui_checked(&missing_code_lens, missing_code_ct * sizeof(int))) {
      goto oxford_sample_load_ret_NOMEM;
    }
    bufptr = missing_code;
    for (uii = 0; uii < missing_code_ct; uii++) {
      missing_code_ptrs[uii] = bufptr;
      bufptr = strchr(bufptr, ',');
      if (bufptr) {
	missing_code_lens[uii] = (unsigned int)(bufptr - missing_code_ptrs[uii]);
	*bufptr = '\0';
	bufptr++;
      } else {
	missing_code_lens[uii] = strlen(missing_code_ptrs[uii]);
      }
    }
  }
  if (fseeko(samplefile, first_real_line_loc, SEEK_SET)) {
    goto oxford_sample_load_ret_READ_FAIL;
  }
  while (fgets(tbuf, MAXLINELEN, samplefile)) {
    if (*tbuf == '\n') {
      continue;
    }
    item_begin = skip_initial_spaces(tbuf);
    bufptr = item_end(item_begin);
    uii = (unsigned int)(bufptr - item_begin);
    memcpy(&(person_ids[indiv_uidx * max_person_id_len]), item_begin, uii);
    person_ids[indiv_uidx * max_person_id_len + uii] = '\t';
    item_begin = skip_initial_spaces(bufptr);
    bufptr = item_end(item_begin);
    ujj = (unsigned int)(bufptr - item_begin);
    memcpy(&(person_ids[indiv_uidx * max_person_id_len + uii + 1]), item_begin, ujj);
    person_ids[indiv_uidx * max_person_id_len + uii + 1 + ujj] = '\0';
    // assume columns 3-5 are missing, sex, pheno
    item_begin = next_item(next_item(skip_initial_spaces(bufptr)));
    if (no_more_items(item_begin)) {
      goto oxford_sample_load_ret_INVALID_FORMAT;
    }
    bufptr = item_end(item_begin);
    uii = (unsigned int)(bufptr - item_begin);
    is_missing = 0;
    for (ujj = 0; ujj < missing_code_ct; ujj++) {
      if (uii == missing_code_lens[ujj]) {
	if (!memcmp(item_begin, missing_code_ptrs[ujj], uii)) {
	  set_bit_noct(*pheno_exclude_ptr, indiv_uidx);
	  is_missing = 1;
	  break;
	}
      }
    }
    if (!is_missing) {
      if (sscanf(item_begin, "%lg", &dxx) != 1) {
	goto oxford_sample_load_ret_INVALID_FORMAT;
      }
      phenos[indiv_uidx] = dxx;
    }
    indiv_uidx++;
  }
  if (!feof(samplefile)) {
    goto oxford_sample_load_ret_READ_FAIL;
  }

  retval = 0;
  while (0) {
  oxford_sample_load_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  oxford_sample_load_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  oxford_sample_load_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  oxford_sample_load_ret_INVALID_FORMAT:
    printf("Error: Improperly formatted .sample file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  oxford_sample_load_ret_INVALID_FORMAT_2:
    printf("Error: Excessively long line in .sample file (max %d chars).\n", MAXLINELEN - 3);
  oxford_sample_load_ret_INVALID_FORMAT_3:
    retval = RET_INVALID_FORMAT;
    break;
  }
  if (wkspace_mark) {
    wkspace_reset(wkspace_mark);
  }
  fclose_cond(samplefile);
  return retval;
}

int oxford_gen_scan(FILE* genfile, unsigned int* gen_buf_len_ptr, unsigned int* unfiltered_marker_ct_ptr, unsigned int unfiltered_indiv_ct) {
  // For now, just determine maximum line length.
  // Assumes genfile points to beginning of file.
  unsigned int unfiltered_marker_ct = 0;
  unsigned int gen_buf_len = 0;
  char* stack_base = (char*)wkspace_base;
  unsigned int max_load;
  unsigned int uii;
  if (wkspace_left > 2147483584) {
    max_load = 2147483584;
  } else {
    max_load = wkspace_left;
  }
  while (fgets(stack_base, max_load, genfile)) {
    if (*stack_base == '\n') {
      continue;
    }
    uii = strlen(stack_base);
    if (uii >= gen_buf_len) {
      gen_buf_len = uii + 1;
    }
    unfiltered_marker_ct++;
  }
  if (!feof(genfile)) {
    return RET_READ_FAIL;
  }
  *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
  *gen_buf_len_ptr = gen_buf_len;
  return 0;
}

int distance_d_write() {
  return 0;
}

// ----- multithread globals -----
static double* distance_matrix;
static double* distance_wt_matrix;
static unsigned int thread_start[MAX_THREADS_P1];
static unsigned int indiv_ct;
static double* dosage_vals;
static double* dosage_wts;

int oxford_distance_calc(FILE* genfile, unsigned int gen_buf_len, unsigned int unfiltered_marker_ct, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, int distance_3d, unsigned int thread_ct, int parallel_idx, unsigned int parallel_tot) {
  double* dvptr;
  double* dwptr;
  unsigned char* wkspace_mark;
  char* loadbuf;
  char* bufptr;
  long long llxx;
  unsigned long ulii;
  unsigned int marker_uidx;
  unsigned int indiv_uidx;
  int retval;
  double pzero;
  double pone;
  double ptwo;
  triangle_fill(thread_start, indiv_ct, thread_ct, parallel_idx, parallel_tot, 1, 1);
  llxx = thread_start[thread_ct];
  llxx = ((llxx * (llxx - 1)) - (long long)thread_start[0] * (thread_start[0] - 1)) / 2;
#ifndef __LP64__
  if (llxx > 4294967295LL) {
    return RET_NOMEM;
  }
#endif
  ulii = (unsigned long)llxx;
  if (wkspace_alloc_d_checked(&distance_matrix, ulii * sizeof(double))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_d_checked(&distance_wt_matrix, ulii * sizeof(double))) {
    return RET_NOMEM;
  }
  fill_double_zero(distance_matrix, ulii);
  fill_double_zero(distance_wt_matrix, ulii);
  wkspace_mark = wkspace_base;
  if (wkspace_alloc_c_checked(&loadbuf, gen_buf_len)) {
    return RET_NOMEM;
  }
  if (distance_3d) {
    printf("Error: 3d distance calculation not yet supported.\n");
    return RET_CALC_NOT_YET_SUPPORTED;
  } else {
    if (wkspace_alloc_d_checked(&dosage_vals, indiv_ct * sizeof(double))) {
      return RET_NOMEM;
    }
    if (wkspace_alloc_d_checked(&dosage_wts, indiv_ct * sizeof(double))) {
      return RET_NOMEM;
    }
  }
  marker_uidx = 0;
  rewind(genfile);
  while (fgets(loadbuf, gen_buf_len, genfile)) {
    if (*loadbuf == '\n') {
      continue;
    }
    bufptr = next_item(next_item(next_item(next_item(next_item(skip_initial_spaces(loadbuf))))));
    if (distance_3d) {
      // TBD
    } else {
      dvptr = dosage_vals;
      dwptr = dosage_wts;
      for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
	if (is_set(indiv_exclude, indiv_uidx)) {
          bufptr = next_item(next_item(next_item(bufptr)));
	} else {
	  if (no_more_items(bufptr)) {
	    goto oxford_distance_calc_ret_INVALID_FORMAT;
	  }
	  if (sscanf(bufptr, "%lg", &pzero) != 1) {
	    goto oxford_distance_calc_ret_INVALID_FORMAT;
	  }
	  bufptr = next_item(bufptr);
	  if (no_more_items(bufptr)) {
	    goto oxford_distance_calc_ret_INVALID_FORMAT;
	  }
	  if (sscanf(bufptr, "%lg", &pone) != 1) {
	    goto oxford_distance_calc_ret_INVALID_FORMAT;
	  }
	  bufptr = next_item(bufptr);
	  if (no_more_items(bufptr)) {
	    goto oxford_distance_calc_ret_INVALID_FORMAT;
	  }
	  if (sscanf(bufptr, "%lg", &ptwo) != 1) {
	    goto oxford_distance_calc_ret_INVALID_FORMAT;
	  }
          *dvptr++ = pone + 2 * ptwo;
	  *dwptr++ = pzero + pone + ptwo;
	}
      }
      // launch threads
    }
  }
  retval = 0;
  while (0) {
  oxford_distance_calc_ret_INVALID_FORMAT:
    printf("Error: Improperly formatted .gen file.\n");
    retval = RET_INVALID_FORMAT;
  }
  wkspace_reset(wkspace_mark);
  return retval;
}

int wdist_dosage(int calculation_type, char* genname, char* samplename, char* missing_code, int distance_3d, unsigned int thread_ct, int parallel_idx, unsigned int parallel_tot) {
  FILE* genfile = NULL;
  unsigned int gen_buf_len;
  unsigned char* wkspace_mark;
  double* phenos;
  unsigned long* pheno_exclude;
  unsigned int unfiltered_marker_ct;
  unsigned int unfiltered_indiv_ct;
  unsigned long* indiv_exclude;
  char* person_ids;
  unsigned int max_person_id_len;
  int retval;
  retval = oxford_sample_load(samplename, &unfiltered_indiv_ct, &person_ids, &max_person_id_len, &phenos, &pheno_exclude, &indiv_exclude, missing_code);
  if (retval) {
    goto wdist_dosage_ret_1;
  }
  if (fopen_checked(&genfile, genname, "r")) {
    return RET_OPEN_FAIL;
  }
  retval = oxford_gen_scan(genfile, &gen_buf_len, &unfiltered_marker_ct, unfiltered_indiv_ct);
  if (retval) {
    goto wdist_dosage_ret_1;
  }
  if (distance_req(calculation_type)) {
    wkspace_mark = wkspace_base;
    indiv_ct = unfiltered_indiv_ct;
    retval = oxford_distance_calc(genfile, gen_buf_len, unfiltered_marker_ct, unfiltered_indiv_ct, indiv_exclude, distance_3d, thread_ct, parallel_idx, parallel_tot);
    if (retval) {
      goto wdist_dosage_ret_1;
    }
    retval = distance_d_write();
    if (retval) {
      goto wdist_dosage_ret_1;
    }
    wkspace_reset(wkspace_mark);
  }
  wdist_dosage_ret_1:
  fclose_cond(genfile);
  return retval;
}
