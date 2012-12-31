#include "wdist_common.h"

// Routines that handle dosage data instead of just 0-1-2 reference allele
// counts.  Only Oxford-formatted data is currently supported, but a PLINK
// --dosage loader will probably be added later.

#define MULTIPLEX_DOSAGE_NM (BITCT / 2)
#define MULTIPLEX_DOSAGE_NM3D (MULTIPLEX_DOSAGE_NM / 2)
#if MULTIPLEX_DOSAGE_NM < MULTIPLEX_DOSAGE_NM3D
#error "MULTIPLEX_DOSAGE_NM cannot be less than MULTIPLEX_DOSAGE_NM3D."
#endif


// allowed deviation from 1.0 when summing 0-1-2 reference allele counts, when
// assessing whether missingness should be treated as binary
#define D_EPSILON 0.000244140625 // just want this above .0002

int oxford_sample_load(char* samplename, unsigned int* unfiltered_indiv_ct_ptr, char** person_ids_ptr, unsigned int* max_person_id_len_ptr, unsigned char** sex_info_ptr, char** pheno_c_ptr, double** pheno_d_ptr, unsigned long** pheno_exclude_ptr, unsigned long** indiv_exclude_ptr, char* missing_code) {
  FILE* samplefile = NULL;
  unsigned char* wkspace_mark = NULL;
  unsigned int unfiltered_indiv_ct = 0;
  unsigned int max_person_id_len = 4;
  unsigned int missing_code_ct = 0;
  unsigned int indiv_uidx = 0;
  char** missing_code_ptrs = NULL;
  unsigned int* missing_code_lens = NULL;
  int load_sex = 0;
  int load_pheno = 0;
  unsigned int male_ct = 0;
  unsigned int female_ct = 0;
  unsigned char* sex_info = NULL;
  char* pheno_c = NULL;
  double* pheno_d = NULL;
  int pheno_is_binary = 0;
  unsigned int unfiltered_indiv_ctl;
  char* person_ids;
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
  item_begin = skip_initial_spaces(tbuf);
  bufptr = next_item_mult(item_begin, 2);
  if (no_more_items(bufptr)) {
    goto oxford_sample_load_ret_INVALID_FORMAT_3;
  }
  if (strcmp_se(item_begin, "ID_1", 4) || strcmp_se(next_item(item_begin), "ID_2", 4) || strcmp_se(bufptr, "missing", 7)) {
    goto oxford_sample_load_ret_INVALID_FORMAT_3;
  }
  bufptr = next_item(bufptr);
  if (!no_more_items(bufptr)) {
    if (strcmp_se(bufptr, "sex", 3)) {
      *sex_info_ptr = NULL;
      load_pheno = 1;
    } else {
      load_sex = 1;
      bufptr = next_item(bufptr);
      if (!no_more_items(bufptr)) {
	load_pheno = 1;
      }
    }
    if (load_pheno) {
      if (!no_more_items(next_item(bufptr))) {
	load_pheno++;
      }
    }
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
  item_begin = skip_initial_spaces(tbuf);
  bufptr = next_item_mult(item_begin, 2);
  if (no_more_items(bufptr)) {
    goto oxford_sample_load_ret_INVALID_FORMAT;
  }
  if (strcmp_se(item_begin, "0", 1) || strcmp_se(next_item(item_begin), "0", 1) || strcmp_se(bufptr, "0", 1)) {
    goto oxford_sample_load_ret_INVALID_FORMAT;
  }
  if (load_sex) {
    bufptr = next_item(bufptr);
    if (no_more_items(bufptr) || strcmp_se(bufptr, "D", 1)) {
      goto oxford_sample_load_ret_INVALID_FORMAT;
    }
  }
  if (load_pheno) {
    bufptr = next_item(bufptr);
    if (no_more_items(bufptr)) {
      goto oxford_sample_load_ret_INVALID_FORMAT;
    }
    if (!strcmp_se(bufptr, "B", 1)) {
      pheno_is_binary = 1;
    } else if (strcmp_se(bufptr, "P", 1)) {
      goto oxford_sample_load_ret_INVALID_FORMAT;
    }
    if (load_pheno > 1) {
      if (no_more_items(next_item(bufptr))) {
        goto oxford_sample_load_ret_INVALID_FORMAT;
      }
      logprint("Note: Multiple phenotypes/covariates detected.  Loading only the first one.\n");
      load_pheno = 1;
    }
  }
  
  if (!tbuf[MAXLINELEN - 1]) {
    goto oxford_sample_load_ret_INVALID_FORMAT_2;
  }
  first_real_line_loc = ftello(samplefile);
  unfiltered_indiv_ct = 0;
  // first pass: count indivs, determine max ID length
  while (fgets(tbuf, MAXLINELEN, samplefile) != NULL) {
    if (is_eoln(*tbuf)) {
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
    logprint("Error: No individuals in .sample file.\n");
    goto oxford_sample_load_ret_INVALID_FORMAT_4;
  }
  *unfiltered_indiv_ct_ptr = unfiltered_indiv_ct;
  *max_person_id_len_ptr = max_person_id_len;

  if (wkspace_alloc_c_checked(person_ids_ptr, unfiltered_indiv_ct * max_person_id_len)) {
    goto oxford_sample_load_ret_NOMEM;
  }
  if (load_sex) {
    if (wkspace_alloc_uc_checked(sex_info_ptr, unfiltered_indiv_ct)) {
      goto oxford_sample_load_ret_NOMEM;
    }
    sex_info = *sex_info_ptr;
    memset(sex_info, 0, unfiltered_indiv_ct);
  }
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  if (load_pheno) {
    if (pheno_is_binary) {
      if (wkspace_alloc_c_checked(pheno_c_ptr, unfiltered_indiv_ct)) {
	goto oxford_sample_load_ret_NOMEM;
      }
      pheno_c = *pheno_c_ptr;
      memset(pheno_c, 255, unfiltered_indiv_ct);
    } else {
      if (wkspace_alloc_d_checked(pheno_d_ptr, unfiltered_indiv_ct * sizeof(double))) {
	goto oxford_sample_load_ret_NOMEM;
      }
      pheno_d = *pheno_d_ptr;
    }
    if (wkspace_alloc_ul_checked(pheno_exclude_ptr, unfiltered_indiv_ctl * sizeof(long))) {
      goto oxford_sample_load_ret_NOMEM;
    }
    fill_ulong_zero(*pheno_exclude_ptr, unfiltered_indiv_ctl);
  } else {
    pheno_d_ptr = NULL;
    pheno_exclude_ptr = NULL;
  }
  if (wkspace_alloc_ul_checked(indiv_exclude_ptr, unfiltered_indiv_ctl * sizeof(long))) {
    goto oxford_sample_load_ret_NOMEM;
  }
  person_ids = *person_ids_ptr;
  fill_ulong_zero(*indiv_exclude_ptr, unfiltered_indiv_ctl);
  wkspace_mark = wkspace_base;
  if (*missing_code) {
    bufptr = missing_code;
    do {
      if ((*bufptr == ',') || (*bufptr == '\0')) {
	// blank string makes no sense
        sprintf(logbuf, "Error: Invalid --missing-code parameter '%s'.%s", missing_code, errstr_append);
	logprintb();
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
  // second pass: actually load stuff
  while (fgets(tbuf, MAXLINELEN, samplefile)) {
    if (is_eoln(*tbuf)) {
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
    item_begin = next_item(skip_initial_spaces(bufptr));
    if (load_sex) {
      if (no_more_items(item_begin)) {
        goto oxford_sample_load_ret_INVALID_FORMAT;
      }
      if (*item_begin == '1') {
        sex_info[indiv_uidx] = 1;
        male_ct++;
      } else if (*item_begin == '2') {
        sex_info[indiv_uidx] = 2;
        female_ct++;
      }
      item_begin = next_item(item_begin);
    }
    if (load_pheno) {
      if (no_more_items(item_begin)) {
	goto oxford_sample_load_ret_INVALID_FORMAT;
      }
      bufptr = item_endnn(item_begin);
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
	if (pheno_is_binary) {
          if (*item_begin == '0') {
	    pheno_c[indiv_uidx] = 0;
	  } else if (*item_begin == '1') {
	    pheno_c[indiv_uidx] = 1;
	  } else {
	    goto oxford_sample_load_ret_INVALID_FORMAT;
	  }
	} else {
	  if (sscanf(item_begin, "%lg", &dxx) != 1) {
	    goto oxford_sample_load_ret_INVALID_FORMAT;
	  }
	  pheno_d[indiv_uidx] = dxx;
	}
      }
    }
    indiv_uidx++;
  }
  if (!feof(samplefile)) {
    goto oxford_sample_load_ret_READ_FAIL;
  }

  retval = 0;
  if (load_sex) {
    uii = unfiltered_indiv_ct - male_ct - female_ct;
    if (uii) {
      sprintf(logbuf, "%d people (%d male%s, %d female%s, %d unknown) loaded.\n", unfiltered_indiv_ct, male_ct, (male_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", uii);
    } else {
      sprintf(logbuf, "%d people (%d male%s, %d female%s) loaded.\n", unfiltered_indiv_ct, male_ct, (male_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s");
    }
  } else {
    sprintf(logbuf, "%d people loaded.\n", unfiltered_indiv_ct);
  }
  logprintb();
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
    logprint("Error: Improperly formatted .sample file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  oxford_sample_load_ret_INVALID_FORMAT_2:
    sprintf(logbuf, "Error: Excessively long line in .sample file (max %d chars).\n", MAXLINELEN - 3);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  oxford_sample_load_ret_INVALID_FORMAT_3:
    logprint("Unsupported .sample file format.  (First three columns are expected to be\n'ID_1', 'ID_2', and 'missing', in that order.)\n");
  oxford_sample_load_ret_INVALID_FORMAT_4:
    retval = RET_INVALID_FORMAT;
    break;
  }
  if (wkspace_mark) {
    wkspace_reset(wkspace_mark);
  }
  fclose_cond(samplefile);
  return retval;
}

int oxford_gen_load1(FILE* genfile, unsigned int* gen_buf_len_ptr, unsigned int* unfiltered_marker_ct_ptr, double** set_allele_freqs_ptr, int* is_missing_01_ptr, unsigned int unfiltered_indiv_ct, int maf_succ) {
  // Determine maximum line length, calculate reference allele frequencies,
  // and check if all missingness probabilities are 0 or ~1.
  unsigned int unfiltered_marker_ct = 0;
  unsigned int unfiltered_indiv_ct8m = unfiltered_indiv_ct & 0xfffffff8U;
  unsigned int gen_buf_len = 0;
  char* stack_base = (char*)wkspace_base;
  char* loadbuf = (&(stack_base[sizeof(double)]));
  int is_missing_01 = 1;
  unsigned int pct = 1;
  unsigned long long file_length;
  unsigned long long file_pos_100;
  char* bufptr;
  int max_load;
  unsigned int uii;
  unsigned int indiv_uidx;
  double total_ref_allele_ct;
  double total_allele_wt;
  double cur_ref_allele_cts[8];
  double cur_allele_wts[8];
  double cur_ref_homs[8];
  double cur_ref_allele_ct;
  double cur_allele_wt;
  double dxx;
  if (wkspace_left > 2147483584) {
    max_load = 2147483584;
  } else {
    max_load = wkspace_left;
  }
  *set_allele_freqs_ptr = (double*)wkspace_base;
  if (fseeko(genfile, 0, SEEK_END)) {
    return RET_READ_FAIL;
  }
  fflush(stdout);
  file_length = ftello(genfile);
  rewind(genfile);
  max_load -= sizeof(double);
  if (max_load <= 0) {
    return RET_NOMEM;
  }
  while (fgets(loadbuf, max_load, genfile)) {
    if (is_eoln(*loadbuf)) {
      continue;
    }
    bufptr = next_item_mult(skip_initial_spaces(loadbuf), 4);
    if (maf_succ) {
      total_ref_allele_ct = 1.0;
      total_allele_wt = 1.0;
    } else {
      total_ref_allele_ct = 0.0;
      total_allele_wt = 0.0;
    }
    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct;) {
      bufptr = next_item(bufptr);
      if (no_more_items(bufptr)) {
	goto oxford_gen_load1_ret_INVALID_FORMAT;
      }
      if (indiv_uidx >= unfiltered_indiv_ct8m) {
	if (sscanf(bufptr, "%lg %lg %lg", &cur_allele_wt, &cur_ref_allele_ct, &dxx) != 3) {
	  goto oxford_gen_load1_ret_INVALID_FORMAT;
	}
	bufptr = next_item_mult(bufptr, 2);
	cur_allele_wt += dxx + cur_ref_allele_ct;
	cur_ref_allele_ct += 2 * dxx;

	if (is_missing_01) {
	  if ((cur_allele_wt != 0.0) && ((cur_allele_wt < 1.0 - D_EPSILON) || (cur_allele_wt > 1.0 + D_EPSILON))) {
	    is_missing_01 = 0;
	  }
	}
	total_ref_allele_ct += cur_ref_allele_ct;
	total_allele_wt += cur_allele_wt;
	indiv_uidx++;
      } else {
	// sadly, this kludge is an important performance optimization
	if (sscanf(bufptr, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &(cur_allele_wts[0]), &(cur_ref_allele_cts[0]), &(cur_ref_homs[0]), &(cur_allele_wts[1]), &(cur_ref_allele_cts[1]), &(cur_ref_homs[1]), &(cur_allele_wts[2]), &(cur_ref_allele_cts[2]), &(cur_ref_homs[2]), &(cur_allele_wts[3]), &(cur_ref_allele_cts[3]), &(cur_ref_homs[3]), &(cur_allele_wts[4]), &(cur_ref_allele_cts[4]), &(cur_ref_homs[4]), &(cur_allele_wts[5]), &(cur_ref_allele_cts[5]), &(cur_ref_homs[5]), &(cur_allele_wts[6]), &(cur_ref_allele_cts[6]), &(cur_ref_homs[6]), &(cur_allele_wts[7]), &(cur_ref_allele_cts[7]), &(cur_ref_homs[7])) != 24) {
	  goto oxford_gen_load1_ret_INVALID_FORMAT;
	}
	bufptr = next_item_mult(bufptr, 23);
	if (is_missing_01) {
	  for (uii = 0; uii < 8; uii++) {
	    cur_allele_wts[uii] += cur_ref_homs[uii] + cur_ref_allele_cts[uii];
	    cur_ref_allele_cts[uii] += 2 * cur_ref_homs[uii];
	    if ((cur_allele_wts[uii] != 0.0) && ((cur_allele_wts[uii] < 1.0 - D_EPSILON) || (cur_allele_wts[uii] > 1.0 + D_EPSILON))) {
	      is_missing_01 = 0;
	    }
	    total_ref_allele_ct += cur_ref_allele_cts[uii];
	    total_allele_wt += cur_allele_wts[uii];
	  }
	} else {
	  for (uii = 0; uii < 8; uii++) {
	    cur_allele_wts[uii] += cur_ref_homs[uii] + cur_ref_allele_cts[uii];
	    cur_ref_allele_cts[uii] += 2 * cur_ref_homs[uii];
	    total_ref_allele_ct += cur_ref_allele_cts[uii];
	    total_allele_wt += cur_allele_wts[uii];
	  }
	}
	indiv_uidx += 8;
      }
    }
    uii = strlen(bufptr) + (unsigned int)(bufptr - loadbuf);
    if (loadbuf[uii - 1] != '\n') {
      logprint("Excessively long line in .gen file.\n");
      return RET_NOMEM;
    }
    if (uii >= gen_buf_len) {
      gen_buf_len = uii + 1;
    }
    if (total_allele_wt == 0.0) {
      (*set_allele_freqs_ptr)[unfiltered_marker_ct] = 0.5;
    } else {
      (*set_allele_freqs_ptr)[unfiltered_marker_ct] = (total_ref_allele_ct * 0.5) / total_allele_wt;
    }
    unfiltered_marker_ct++;
    loadbuf = &(loadbuf[sizeof(double)]);
    max_load -= sizeof(double);
    file_pos_100 = ftello(genfile) * 100LLU;
    if (file_pos_100 >= pct * file_length) {
      pct = (unsigned int)(file_pos_100 / file_length);
      printf("\rScanning .gen file (%u%%)...", pct);
      fflush(stdout);
    }
  }
  if (!feof(genfile)) {
    return RET_READ_FAIL;
  }
  if (!unfiltered_marker_ct) {
    logprint("Error: No markers in .gen file.\n");
    return RET_INVALID_FORMAT;
  }
  putchar('\r');
  sprintf(logbuf, ".gen scan complete.  %u markers and %u individuals present.\n", unfiltered_marker_ct, unfiltered_indiv_ct);
  logprintb();
  *set_allele_freqs_ptr = (double*)wkspace_alloc(unfiltered_marker_ct * sizeof(double));
  *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
  *gen_buf_len_ptr = gen_buf_len;
  *is_missing_01_ptr = is_missing_01;
  while (0) {
  oxford_gen_load1_ret_INVALID_FORMAT:
    logprint("Error: Improperly formatted .gen file.\n");
    return RET_INVALID_FORMAT;
  }
  return 0;
}

// ----- multithread globals -----
static double* g_distance_matrix;

// If missingness is binary, this is a sum of sparse intersection weights;
// otherwise, this is a sum of ALL weights.
static double* g_distance_wt_matrix;

static unsigned int g_thread_start[MAX_THREADS_P1];
static unsigned int g_indiv_ct;
static double* g_dosage_vals; // (usually) [0..2] dosages for current SNPs

// If missingness is binary
static unsigned long* g_missing_vals; // bit array marking missing values

// If missingness is continuous
static double* g_nonmissing_vals;

#if BITCT < MULTIPLEX_DOSAGE_NM
#error "Insufficient g_missing_wts[] size in wdist_dosage.c."
#endif
static double g_missing_wts[BITCT]; // missingness rescale weights
#if __LP64__
static double* g_missing_dmasks; // 0x7fff... if non-missing, 0x0000... otherwise
#else
// double times signed int tends to be faster than double times unsigned int
static int* g_missing_dmasks; // 1 or 0, multiply instead of mask
#endif

void incr_distance_dosage_2d_01(double* distance_matrix_slice, int thread_idx) {
#if __LP64__
  // take absolute value = force sign bit to zero
  __m128d* dptr_start;
  __m128d* dptr_end;
  __m128d* dptr;
  __m128d* dptr2;
  __m128d* mptr_start;
  __m128d* mptr;
  __m128d* mptr2;
  __uni16 acc;
#else
  double* dptr_start;
  double* dptr_end;
  double* dptr;
  double* dptr2;
  int* mptr_start;
  int* mptr;
  int* mptr2;
  double acc;
#endif
  unsigned long ulii;
  unsigned long uljj;
  for (ulii = g_thread_start[thread_idx]; ulii < g_thread_start[thread_idx + 1]; ulii++) {
#if __LP64__
    dptr_start = (__m128d*)(&(g_dosage_vals[ulii * BITCT]));
    dptr_end = &(dptr_start[BITCT2]);
    dptr2 = (__m128d*)g_dosage_vals;
    mptr2 = (__m128d*)g_missing_dmasks;
#else
    dptr_start = &(g_dosage_vals[ulii * BITCT]);
    dptr_end = &(dptr_start[BITCT]);
    dptr2 = g_dosage_vals;
    mptr2 = g_missing_dmasks;
#endif
    if (g_missing_vals[ulii]) {
#if __LP64__
      mptr_start = (__m128d*)(&(g_missing_dmasks[ulii * BITCT]));
#else
      mptr_start = &(g_missing_dmasks[ulii * BITCT]);
#endif
      for (uljj = 0; uljj < ulii; uljj++) {
	dptr = dptr_start;
	mptr = mptr_start;
#if __LP64__
	acc.vd = _mm_setzero_pd();
	do {
	  acc.vd = _mm_add_pd(acc.vd, _mm_and_pd(_mm_and_pd(*mptr++, *mptr2++), _mm_sub_pd(*dptr++, *dptr2++)));
	} while (dptr != dptr_end);
        *distance_matrix_slice += acc.d8[0] + acc.d8[1];
#else
	acc = 0.0;
	do {
	  acc += (*mptr++) * (*mptr2++) * fabs((*dptr++) - (*dptr2++));
	} while (dptr != dptr_end);
	*distance_matrix_slice += acc;
#endif
        *distance_matrix_slice++;
      }
    } else {
      for (uljj = 0; uljj < ulii; uljj++) {
        dptr = dptr_start;
#if __LP64__
	acc.vd = _mm_setzero_pd();
	do {
	  acc.vd = _mm_add_pd(acc.vd, _mm_and_pd(*mptr2++, _mm_sub_pd(*dptr++, *dptr2++)));
	} while (dptr != dptr_end);
        *distance_matrix_slice += acc.d8[0] + acc.d8[1];
#else
	acc = 0.0;
	do {
	  acc += (*mptr2++) * fabs((*dptr++) - (*dptr2++));
	} while (dptr != dptr_end);
	*distance_matrix_slice += acc;
#endif
        *distance_matrix_slice++;
      }
    }
  }
}

void* incr_distance_dosage_2d_01_thread(void* arg) {
  long tidx = (long)arg;
  int ts = g_thread_start[tidx];
  int ts0 = g_thread_start[0];
  incr_distance_dosage_2d_01(&(g_distance_matrix[((long long)ts * (ts - 1) - (long long)ts0 * (ts0 - 1)) / 2]), (int)tidx);
  return NULL;
}

void incr_distance_dosage_2d_flat(double* distance_matrix_slice, double* distance_wt_matrix_slice, int thread_idx) {
#if __LP64__
  // take absolute value = force sign bit to zero
  const __m128i absmask_raw = {0x7fffffffffffffffLU, 0x7fffffffffffffffLU};
  __m128d* absmask_ptr = (__m128d*)(&absmask_raw);
  __m128d absmask = *absmask_ptr;
  __m128d* dptr_start;
  __m128d* dptr_end;
  __m128d* dptr;
  __m128d* dptr2;
  __m128d* mptr_start;
  __m128d* mptr;
  __m128d* mptr2;
  __m128d cur_nonmissing_wts;
  __uni16 acc;
  __uni16 accm;
#else
  double* dptr_start;
  double* dptr_end;
  double* dptr;
  double* dptr2;
  double* mptr_start;
  double* mptr;
  double* mptr2;
  double cur_nonmissing_wt;
  double acc;
  double accm;
#endif
  unsigned long ulii;
  unsigned long uljj;
  for (ulii = g_thread_start[thread_idx]; ulii < g_thread_start[thread_idx + 1]; ulii++) {
#if __LP64__
    dptr_start = (__m128d*)(&(g_dosage_vals[ulii * MULTIPLEX_DOSAGE_NM]));
    dptr_end = &(dptr_start[MULTIPLEX_DOSAGE_NM / 2]);
    dptr2 = (__m128d*)g_dosage_vals;
    mptr_start = (__m128d*)(&(g_nonmissing_vals[ulii * MULTIPLEX_DOSAGE_NM]));
    mptr2 = (__m128d*)g_nonmissing_vals;
#else
    dptr_start = &(g_dosage_vals[ulii * MULTIPLEX_DOSAGE_NM]);
    dptr_end = &(dptr_start[MULTIPLEX_DOSAGE_NM]);
    dptr2 = g_dosage_vals;
    mptr_start = &(g_nonmissing_vals[ulii * MULTIPLEX_DOSAGE_NM]);
    mptr2 = g_nonmissing_vals;
#endif
    for (uljj = 0; uljj < ulii; uljj++) {
      dptr = dptr_start;
      mptr = mptr_start;
#if __LP64__
      acc.vd = _mm_setzero_pd();
      accm.vd = _mm_setzero_pd();
      do {
	cur_nonmissing_wts = _mm_mul_pd(*mptr++, *mptr2++);
	acc.vd = _mm_add_pd(acc.vd, _mm_mul_pd(cur_nonmissing_wts, _mm_and_pd(absmask, _mm_sub_pd(*dptr++, *dptr2++))));
	accm.vd = _mm_add_pd(accm.vd, cur_nonmissing_wts);
      } while (dptr != dptr_end);
      *distance_matrix_slice += acc.d8[0] + acc.d8[1];
      *distance_wt_matrix_slice += accm.d8[0] + accm.d8[1];
#else
      acc = 0.0;
      accm = 0.0;
      do {
	cur_nonmissing_wt = (*mptr++) * (*mptr2++);
        acc += cur_nonmissing_wt * fabs((*dptr++) - (*dptr2++));
	accm += cur_nonmissing_wt;
      } while (dptr != dptr_end);
      *distance_matrix_slice += acc;
      *distance_wt_matrix_slice += accm;
#endif
      distance_matrix_slice++;
      distance_wt_matrix_slice++;
    }
  }
}

void* incr_distance_dosage_2d_flat_thread(void* arg) {
  long tidx = (long)arg;
  int ts = g_thread_start[tidx];
  int ts0 = g_thread_start[0];
  long long offset = ((long long)ts * (ts - 1) - (long long)ts0 * (ts0 - 1)) / 2;
  incr_distance_dosage_2d_flat(&(g_distance_matrix[offset]), &(g_distance_wt_matrix[offset]), (int)tidx);
  return NULL;
}

void incr_distance_dosage_2d(double* distance_matrix_slice, double* distance_wt_matrix_slice, int thread_idx) {
#if __LP64__
  const __m128i absmask_raw = {0x7fffffffffffffffLU, 0x7fffffffffffffffLU};
  __m128d* absmask_ptr = (__m128d*)(&absmask_raw);
  __m128d absmask = *absmask_ptr;
  __m128d* dptr_start;
  __m128d* dptr_end;
  __m128d* dptr;
  __m128d* dptr2;
  __m128d* mptr_start;
  __m128d* mptr;
  __m128d* mptr2;
  __m128d* mwptr;
  __m128d cur_nonmissing_wts;
  __uni16 acc;
  __uni16 accm;
#else
  double* dptr_start;
  double* dptr_end;
  double* dptr;
  double* dptr2;
  double* mptr_start;
  double* mptr;
  double* mptr2;
  double* mwptr;
  double cur_nonmissing_wt;
  double acc;
  double accm;
#endif
  unsigned long ulii;
  unsigned long uljj;
  for (ulii = g_thread_start[thread_idx]; ulii < g_thread_start[thread_idx + 1]; ulii++) {
#if __LP64__
    dptr_start = (__m128d*)(&(g_dosage_vals[ulii * MULTIPLEX_DOSAGE_NM]));
    dptr_end = &(dptr_start[MULTIPLEX_DOSAGE_NM / 2]);
    dptr2 = (__m128d*)g_dosage_vals;
    mptr_start = (__m128d*)(&(g_nonmissing_vals[ulii * MULTIPLEX_DOSAGE_NM]));
    mptr2 = (__m128d*)g_nonmissing_vals;
#else
    dptr_start = &(g_dosage_vals[ulii * MULTIPLEX_DOSAGE_NM]);
    dptr_end = &(dptr_start[MULTIPLEX_DOSAGE_NM]);
    dptr2 = g_dosage_vals;
    mptr_start = &(g_nonmissing_vals[ulii * MULTIPLEX_DOSAGE_NM]);
    mptr2 = g_nonmissing_vals;
#endif
    for (uljj = 0; uljj < ulii; uljj++) {
      dptr = dptr_start;
      mptr = mptr_start;
#if __LP64__
      mwptr = (__m128d*)g_missing_wts;
      acc.vd = _mm_setzero_pd();
      accm.vd = _mm_setzero_pd();
      do {
	cur_nonmissing_wts = _mm_mul_pd(*mptr++, *mptr2++);
	acc.vd = _mm_add_pd(acc.vd, _mm_mul_pd(cur_nonmissing_wts, _mm_and_pd(absmask, _mm_sub_pd(*dptr++, *dptr2++))));
	accm.vd = _mm_add_pd(accm.vd, _mm_mul_pd(*mwptr++, cur_nonmissing_wts));
      } while (dptr != dptr_end);
      *distance_matrix_slice += acc.d8[0] + acc.d8[1];
      *distance_wt_matrix_slice += accm.d8[0] + accm.d8[1];
#else
      mwptr = g_missing_wts;
      acc = 0.0;
      accm = 0.0;
      do {
        cur_nonmissing_wt = (*mptr++) * (*mptr2++);
        acc += cur_nonmissing_wt * fabs((*dptr++) - (*dptr2++));
	accm += cur_nonmissing_wt * (*mwptr++);
      } while (dptr != dptr_end);
      *distance_matrix_slice += acc;
      *distance_wt_matrix_slice += accm;
#endif
      distance_matrix_slice++;
      distance_wt_matrix_slice++;
    }
  }
}

void* incr_distance_dosage_2d_thread(void* arg) {
  long tidx = (long)arg;
  int ts = g_thread_start[tidx];
  int ts0 = g_thread_start[0];
  long long offset = ((long long)ts * (ts - 1) - (long long)ts0 * (ts0 - 1)) / 2;
  incr_distance_dosage_2d(&(g_distance_matrix[offset]), &(g_distance_wt_matrix[offset]), (int)tidx);
  return NULL;
}

void incr_distance_dosage_3d_flat(double* distance_matrix_slice, double* distance_wt_matrix_slice, int thread_idx) {
  // 3d distance: P_1(0) * (P_2(1) + 2 * P_2(2)) + P_1(1) * (P_2(0) + P_2(2)) +
  //              P_1(2) * (2 * P_2(0) + P_2(1))
  // weight: (P_1(0) + P_1(1) + P_1(2)) * (P_2(0) + P_2(1) + P_2(2))
  // 
#if __LP64__
  const __m128d vec_two = {2.0, 2.0};
  __m128d dptr_buf[2 * MULTIPLEX_DOSAGE_NM3D];
  __m128d* dptr_buf_end1 = &(dptr_buf[MULTIPLEX_DOSAGE_NM3D / 2]);
  __m128d* dptr_buf_end2 = &(dptr_buf[MULTIPLEX_DOSAGE_NM3D]);
  __m128d* dptr_buf_end3 = &(dptr_buf[3 * MULTIPLEX_DOSAGE_NM3D / 2]);
  __m128d* dptr_buf_end = &(dptr_buf[2 * MULTIPLEX_DOSAGE_NM3D]);
  __m128d* dptr;
  __m128d* dptr_start;
  __m128d* dptr2;
  __m128d* dptr3;
  __uni16 acc;
#else
  double dptr_buf[4 * MULTIPLEX_DOSAGE_NM3D];
  double* dptr_buf_end1 = &(dptr_buf[MULTIPLEX_DOSAGE_NM3D]);
  double* dptr_buf_end2 = &(dptr_buf[2 * MULTIPLEX_DOSAGE_NM3D]);
  double* dptr_buf_end3 = &(dptr_buf[3 * MULTIPLEX_DOSAGE_NM3D]);
  double* dptr_buf_end = &(dptr_buf[4 * MULTIPLEX_DOSAGE_NM3D]);
  double* dptr;
  double* dptr_start;
  double* dptr2;
  double* dptr3;
  double acc;
#endif
  unsigned long ulii;
  unsigned long uljj;
  for (ulii = g_thread_start[thread_idx]; ulii < g_thread_start[thread_idx + 1]; ulii++) {
#if __LP64__
    dptr = dptr_buf;
    dptr_start = (__m128d*)(&(g_dosage_vals[ulii * 4 * MULTIPLEX_DOSAGE_NM3D]));
    dptr2 = &(dptr_start[MULTIPLEX_DOSAGE_NM3D / 2]);
    dptr3 = &(dptr_start[MULTIPLEX_DOSAGE_NM3D]);
    do {
      *dptr++ = _mm_add_pd(_mm_mul_pd(vec_two, *dptr3++), *dptr2++);
    } while (dptr != dptr_buf_end1);
    dptr3 = dptr_start;
    do {
      *dptr++ = _mm_add_pd(*dptr2++, *dptr3++);
    } while (dptr != dptr_buf_end2);
    dptr2 = dptr_start;
    do {
      *dptr++ = _mm_add_pd(_mm_mul_pd(vec_two, *dptr2++), *dptr3++);
    } while (dptr != dptr_buf_end3);
    dptr2 = &(dptr_start[3 * MULTIPLEX_DOSAGE_NM3D / 2]);
    do {
      *dptr++ = *dptr2++;
    } while (dptr != dptr_buf_end);
    dptr2 = (__m128d*)g_dosage_vals;
#else
    dptr = dptr_buf;
    dptr_start = &(g_dosage_vals[ulii * 4 * MULTIPLEX_DOSAGE_NM3D]);
    dptr2 = &(dptr_start[MULTIPLEX_DOSAGE_NM3D]);
    dptr3 = &(dptr_start[2 * MULTIPLEX_DOSAGE_NM3D]);
    do {
      *dptr++ = (*dptr2++) + 2 * (*dptr3++);
    } while (dptr != dptr_buf_end1);
    dptr3 = dptr_start;
    do {
      *dptr++ = (*dptr2++) + (*dptr3++);
    } while (dptr != dptr_buf_end2);
    dptr2 = dptr_start;
    do {
      *dptr++ = 2 * (*dptr2++) + (*dptr3++);
    } while (dptr != dptr_buf_end3);
    dptr2 = &(dptr_start[3 * MULTIPLEX_DOSAGE_NM3D]);
    do {
      *dptr++ = *dptr2++;
    } while (dptr != dptr_buf_end);
    dptr2 = g_dosage_vals;
#endif
    for (uljj = 0; uljj < ulii; uljj++) {
      dptr = dptr_buf;
#if __LP64__
      acc.vd = _mm_setzero_pd();
      do {
	acc.vd = _mm_add_pd(acc.vd, _mm_mul_pd(*dptr++, *dptr2++));
      } while (dptr != dptr_buf_end3);
      *distance_matrix_slice += acc.d8[0] + acc.d8[1];
      acc.vd = _mm_setzero_pd();
      do {
        acc.vd = _mm_add_pd(acc.vd, _mm_mul_pd(*dptr++, *dptr2++));
      } while (dptr != dptr_buf_end);
      *distance_wt_matrix_slice += acc.d8[0] + acc.d8[1];
#else
      acc = 0.0;
      do {
	acc += (*dptr++) * (*dptr2++);
      } while (dptr != dptr_buf_end3);
      *distance_matrix_slice += acc;
      acc = 0.0;
      do {
	acc += (*dptr++) * (*dptr2++);
      } while (dptr != dptr_buf_end);
      *distance_wt_matrix_slice += acc;
#endif
      distance_matrix_slice++;
      distance_wt_matrix_slice++;
    }
  }
}

void* incr_distance_dosage_3d_flat_thread(void* arg) {
  long tidx = (long)arg;
  int ts = g_thread_start[tidx];
  int ts0 = g_thread_start[0];
  long long offset = ((long long)ts * (ts - 1) - (long long)ts0 * (ts0 - 1)) / 2;
  incr_distance_dosage_3d_flat(&(g_distance_matrix[offset]), &(g_distance_wt_matrix[offset]), (int)tidx);
  return NULL;
}

void incr_distance_dosage_3d(double* distance_matrix_slice, double* distance_wt_matrix_slice, int thread_idx) {
#if __LP64__
  const __m128d vec_two = {2.0, 2.0};
  __m128d dptr_buf[2 * MULTIPLEX_DOSAGE_NM3D];
  __m128d* dptr_buf_end1 = &(dptr_buf[MULTIPLEX_DOSAGE_NM3D / 2]);
  __m128d* dptr_buf_end2 = &(dptr_buf[MULTIPLEX_DOSAGE_NM3D]);
  __m128d* dptr_buf_end3 = &(dptr_buf[3 * MULTIPLEX_DOSAGE_NM3D / 2]);
  __m128d* dptr_buf_end = &(dptr_buf[2 * MULTIPLEX_DOSAGE_NM3D]);
  __m128d* dptr;
  __m128d* dptr_start;
  __m128d* dptr2;
  __m128d* dptr3;
  __m128d* mwptr;
  __uni16 acc;
#else
  double dptr_buf[4 * MULTIPLEX_DOSAGE_NM3D];
  double* dptr_buf_end1 = &(dptr_buf[MULTIPLEX_DOSAGE_NM3D]);
  double* dptr_buf_end2 = &(dptr_buf[2 * MULTIPLEX_DOSAGE_NM3D]);
  double* dptr_buf_end3 = &(dptr_buf[3 * MULTIPLEX_DOSAGE_NM3D]);
  double* dptr_buf_end = &(dptr_buf[4 * MULTIPLEX_DOSAGE_NM3D]);
  double* dptr;
  double* dptr_start;
  double* dptr2;
  double* dptr3;
  double* mwptr;
  double acc;
#endif
  unsigned long ulii;
  unsigned long uljj;
  for (ulii = g_thread_start[thread_idx]; ulii < g_thread_start[thread_idx + 1]; ulii++) {
#if __LP64__
    dptr = dptr_buf;
    dptr_start = (__m128d*)(&(g_dosage_vals[ulii * 4 * MULTIPLEX_DOSAGE_NM3D]));
    dptr2 = &(dptr_start[MULTIPLEX_DOSAGE_NM3D / 2]);
    dptr3 = &(dptr_start[MULTIPLEX_DOSAGE_NM3D]);
    do {
      *dptr++ = _mm_add_pd(_mm_mul_pd(vec_two, *dptr3++), *dptr2++);
    } while (dptr != dptr_buf_end1);
    dptr3 = dptr_start;
    do {
      *dptr++ = _mm_add_pd(*dptr2++, *dptr3++);
    } while (dptr != dptr_buf_end2);
    dptr2 = dptr_start;
    do {
      *dptr++ = _mm_add_pd(_mm_mul_pd(vec_two, *dptr2++), *dptr3++);
    } while (dptr != dptr_buf_end3);
    dptr2 = &(dptr_start[3 * MULTIPLEX_DOSAGE_NM3D / 2]);
    mwptr = (__m128d*)g_missing_wts;
    do {
      *dptr++ = _mm_mul_pd(*mwptr++, *dptr2++);
    } while (dptr != dptr_buf_end);
    dptr2 = (__m128d*)g_dosage_vals;
#else
    dptr = dptr_buf;
    dptr_start = &(g_dosage_vals[ulii * 4 * MULTIPLEX_DOSAGE_NM3D]);
    dptr2 = &(dptr_start[MULTIPLEX_DOSAGE_NM3D]);
    dptr3 = &(dptr_start[2 * MULTIPLEX_DOSAGE_NM3D]);
    do {
      *dptr++ = (*dptr2++) + 2 * (*dptr3++);
    } while (dptr != dptr_buf_end1);
    dptr3 = dptr_start;
    do {
      *dptr++ = (*dptr2++) + (*dptr3++);
    } while (dptr != dptr_buf_end2);
    dptr2 = dptr_start;
    do {
      *dptr++ = 2 * (*dptr2++) + (*dptr3++);
    } while (dptr != dptr_buf_end3);
    dptr2 = &(dptr_start[3 * MULTIPLEX_DOSAGE_NM3D]);
    mwptr = g_missing_wts;
    do {
      *dptr++ = (*mwptr++) * (*dptr2++);
    } while (dptr != dptr_buf_end);
    dptr2 = g_dosage_vals;
#endif
    for (uljj = 0; uljj < ulii; uljj++) {
      dptr = dptr_buf;
#if __LP64__
      acc.vd = _mm_setzero_pd();
      do {
	acc.vd = _mm_add_pd(acc.vd, _mm_mul_pd(*dptr++, *dptr2++));
      } while (dptr != dptr_buf_end3);
      *distance_matrix_slice += acc.d8[0] + acc.d8[1];
      acc.vd = _mm_setzero_pd();
      do {
        acc.vd = _mm_add_pd(acc.vd, _mm_mul_pd(*dptr++, *dptr2++));
      } while (dptr != dptr_buf_end);
      *distance_wt_matrix_slice += acc.d8[0] + acc.d8[1];
#else
      acc = 0.0;
      do {
	acc += (*dptr++) * (*dptr2++);
      } while (dptr != dptr_buf_end3);
      *distance_matrix_slice += acc;
      acc = 0.0;
      do {
	acc += (*dptr++) * (*dptr2++);
      } while (dptr != dptr_buf_end);
      *distance_wt_matrix_slice += acc;
#endif
      distance_matrix_slice++;
      distance_wt_matrix_slice++;
    }
  }
}

void* incr_distance_dosage_3d_thread(void* arg) {
  long tidx = (long)arg;
  int ts = g_thread_start[tidx];
  int ts0 = g_thread_start[0];
  long long offset = ((long long)ts * (ts - 1) - (long long)ts0 * (ts0 - 1)) / 2;
  incr_distance_dosage_3d(&(g_distance_matrix[offset]), &(g_distance_wt_matrix[offset]), (int)tidx);
  return NULL;
}

void incr_dosage_missing_wt_01(double* distance_wt_matrix_slice, int thread_idx) {
  // count missing intersection
  unsigned long* mlptr;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long maskii;
  unsigned long mask;
  for (ulii = g_thread_start[thread_idx]; ulii < g_thread_start[thread_idx + 1]; ulii++) {
    mlptr = g_missing_vals;
    maskii = g_missing_vals[ulii];
    if (maskii) {
      for (uljj = 0; uljj < ulii; uljj++) {
	mask = (*mlptr++) & maskii;
	while (mask) {
	  distance_wt_matrix_slice[uljj] += g_missing_wts[__builtin_ctzl(mask)];
	  mask &= mask - 1;
	}
      }
    }
    distance_wt_matrix_slice = &(distance_wt_matrix_slice[ulii]);
  }
}

void* incr_dosage_missing_wt_01_thread(void* arg) {
  long tidx = (long)arg;
  int ts = g_thread_start[tidx];
  int ts0 = g_thread_start[0];
  incr_dosage_missing_wt_01(&(g_distance_wt_matrix[((long long)ts * (ts - 1) - (long long)ts0 * (ts0 - 1)) / 2]), (int)tidx);
  return NULL;
}

int update_distance_dosage_matrix(int is_missing_01, int distance_3d, int distance_flat_missing, unsigned int thread_ct) {
  pthread_t threads[MAX_THREADS];
  unsigned long thread_idx;
  if (!distance_3d) {
    if (is_missing_01) {
      for (thread_idx = 1; thread_idx < thread_ct; thread_idx++) {
	if (pthread_create(&(threads[thread_idx - 1]), NULL, &incr_distance_dosage_2d_01_thread, (void*)thread_idx)) {
	  goto update_distance_dosage_matrix_ret_THREAD_CREATE_FAIL;
	}
      }
      thread_idx = 0;
      incr_distance_dosage_2d_01_thread((void*)thread_idx);
    } else if (distance_flat_missing) {
      for (thread_idx = 1; thread_idx < thread_ct; thread_idx++) {
	if (pthread_create(&(threads[thread_idx - 1]), NULL, &incr_distance_dosage_2d_flat_thread, (void*)thread_idx)) {
	  goto update_distance_dosage_matrix_ret_THREAD_CREATE_FAIL;
	}
      }
      thread_idx = 0;
      incr_distance_dosage_2d_flat_thread((void*)thread_idx);
    } else {
      for (thread_idx = 1; thread_idx < thread_ct; thread_idx++) {
	if (pthread_create(&(threads[thread_idx - 1]), NULL, &incr_distance_dosage_2d_thread, (void*)thread_idx)) {
	  goto update_distance_dosage_matrix_ret_THREAD_CREATE_FAIL;
	}
      }
      thread_idx = 0;
      incr_distance_dosage_2d_thread((void*)thread_idx);
    }
  } else {
    // assume is_missing_01 is not checked for now
    if (distance_flat_missing) {
      for (thread_idx = 1; thread_idx < thread_ct; thread_idx++) {
	if (pthread_create(&(threads[thread_idx - 1]), NULL, &incr_distance_dosage_3d_flat_thread, (void*)thread_idx)) {
	  goto update_distance_dosage_matrix_ret_THREAD_CREATE_FAIL;
	}
      }
      thread_idx = 0;
      incr_distance_dosage_3d_flat_thread((void*)thread_idx);
    } else {
      for (thread_idx = 1; thread_idx < thread_ct; thread_idx++) {
	if (pthread_create(&(threads[thread_idx - 1]), NULL, &incr_distance_dosage_3d_thread, (void*)thread_idx)) {
	  goto update_distance_dosage_matrix_ret_THREAD_CREATE_FAIL;
	}
      }
      thread_idx = 0;
      incr_distance_dosage_3d_thread((void*)thread_idx);
    }
  }
  for (thread_idx = 0; thread_idx < thread_ct - 1; thread_idx++) {
    pthread_join(threads[thread_idx], NULL);
  }
  if (is_missing_01 && (!distance_3d)) {
    for (thread_idx = 1; thread_idx < thread_ct; thread_idx++) {
      if (pthread_create(&(threads[thread_idx - 1]), NULL, &incr_dosage_missing_wt_01_thread, (void*)thread_idx)) {
	goto update_distance_dosage_matrix_ret_THREAD_CREATE_FAIL;
      }
    }
    thread_idx = 0;
    incr_dosage_missing_wt_01_thread((void*)thread_idx);
    for (thread_idx = 0; thread_idx < thread_ct - 1; thread_idx++) {
      pthread_join(threads[thread_idx], NULL);
    }
  }
  return 0;
 update_distance_dosage_matrix_ret_THREAD_CREATE_FAIL:
  logprint(errstr_thread_create);
  while (--thread_idx) {
    pthread_join(threads[thread_idx - 1], NULL);
  }
  return RET_THREAD_CREATE_FAIL;
}

int oxford_distance_calc(FILE* genfile, unsigned int gen_buf_len, double* set_allele_freqs, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_ct, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, int is_missing_01, int distance_3d, int distance_flat_missing, double exponent, unsigned int thread_ct, int parallel_idx, unsigned int parallel_tot) {
  int is_exponent_zero = (exponent == 0.0);
  unsigned int unfiltered_indiv_ct8m = unfiltered_indiv_ct & 0xfffffff8U;
  unsigned int indiv_ctl = (g_indiv_ct + (BITCT - 1)) / BITCT;
  double marker_wt = 1.0;
  double* cur_marker_freqs = NULL;
  double* cmf_ptr = NULL;
  unsigned long* cur_missings = NULL;
  double* cur_missings_d;
  double tot_missing_wt = 0.0;
  unsigned char* wkspace_mark;
  char* loadbuf;
  char* bufptr;
  long long llxx;
  unsigned long ulii;
  unsigned long uljj;
  unsigned int marker_uidx;
  unsigned int marker_idxl;
  unsigned int indiv_uidx;
  unsigned int indiv_idx;
  unsigned int non_missing_ct;
  double* missing_tots;
  double* dptr;
  double* dptr2;
  double* dptr3;
  int retval;
  double pbuf0[8];
  double pbuf1[8];
  double pbuf2[8];
  double pzero;
  double pone;
  double ptwo;
  double dxx;
  unsigned int tstc;
  triangle_fill(g_thread_start, g_indiv_ct, thread_ct, parallel_idx, parallel_tot, 1, 1);
  llxx = g_thread_start[thread_ct];
  llxx = ((llxx * (llxx - 1)) - (long long)g_thread_start[0] * (g_thread_start[0] - 1)) / 2;
#ifndef __LP64__
  if (llxx > 4294967295LL) {
    return RET_NOMEM;
  }
#endif
  ulii = (unsigned long)llxx;
  if (wkspace_alloc_d_checked(&g_distance_matrix, ulii * sizeof(double))) {
    return RET_NOMEM;
  }
  wkspace_mark = wkspace_base;
  if (wkspace_alloc_d_checked(&g_distance_wt_matrix, ulii * sizeof(double))) {
    return RET_NOMEM;
  }
  fill_double_zero(g_distance_matrix, ulii);
  fill_double_zero(g_distance_wt_matrix, ulii);
  if (wkspace_alloc_c_checked(&loadbuf, gen_buf_len)) {
    return RET_NOMEM;
  }
  if (distance_3d) {
    logprint("Error: 3d distance calculation not yet supported.\n");
    return RET_CALC_NOT_YET_SUPPORTED;
  } else {
    if (wkspace_alloc_d_checked(&g_dosage_vals, g_indiv_ct * BITCT * sizeof(double))) {
      return RET_NOMEM;
    }
    if (wkspace_alloc_d_checked(&missing_tots, g_indiv_ct * sizeof(double))) {
      return RET_NOMEM;
    }
    fill_double_zero(missing_tots, g_indiv_ct);
    if (is_missing_01) {
      if (wkspace_alloc_ul_checked(&g_missing_vals, g_indiv_ct * BITCT * sizeof(long))) {
	return RET_NOMEM;
      }
      fill_ulong_zero(g_missing_vals, g_indiv_ct * BITCT);
#if __LP64__
      if (wkspace_alloc_d_checked(&g_missing_dmasks, g_indiv_ct * BITCT * sizeof(double))) {
	return RET_NOMEM;
      }
#else
      if (wkspace_alloc_i_checked(&g_missing_dmasks, g_indiv_ct * BITCT * sizeof(int))) {
	return RET_NOMEM;
      }
#endif
      if (distance_flat_missing) {
	tot_missing_wt = (double)marker_ct;
      } else {
	if (wkspace_alloc_d_checked(&cur_marker_freqs, g_indiv_ct * sizeof(double))) {
	  return RET_NOMEM;
	}
        if (wkspace_alloc_ul_checked(&cur_missings, indiv_ctl * sizeof(long))) {
          return RET_NOMEM;
        }
      }
    } else {
      logprint("Error: Missing probabilities unequal to 0 and 1 not yet supported.\n");
      return RET_CALC_NOT_YET_SUPPORTED;
      if (wkspace_alloc_d_checked(&cur_missings_d, g_indiv_ct * sizeof(double))) {
	return RET_NOMEM;
      }
    }
  }
  marker_uidx = 0;
  rewind(genfile);
  while (fgets(loadbuf, gen_buf_len, genfile)) {
    if (is_eoln(*loadbuf)) {
      continue;
    }
    if (is_set(marker_exclude, marker_uidx)) {
      marker_uidx++;
      continue;
    }
    bufptr = next_item_mult(skip_initial_spaces(loadbuf), 5);
    marker_idxl = marker_uidx % BITCT;
    indiv_idx = 0;
    if (distance_3d) {
      // TBD
    } else {
      if (!is_exponent_zero) {
	dxx = set_allele_freqs[marker_uidx];
	if ((dxx > 0.0) && (dxx < 1.0)) {
          marker_wt = pow(2 * dxx * (1.0 - dxx), -exponent);
	} else {
	  marker_wt = 1.0;
	}
      }
      if (!distance_flat_missing) {
	cmf_ptr = cur_marker_freqs;
	fill_ulong_zero(cur_missings, indiv_ctl);
      }
      for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct;) {
	if (indiv_uidx < unfiltered_indiv_ct8m) {
	  sscanf(bufptr, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &(pbuf0[0]), &(pbuf1[0]), &(pbuf2[0]), &(pbuf0[1]), &(pbuf1[1]), &(pbuf2[1]), &(pbuf0[2]), &(pbuf1[2]), &(pbuf2[2]), &(pbuf0[3]), &(pbuf1[3]), &(pbuf2[3]), &(pbuf0[4]), &(pbuf1[4]), &(pbuf2[4]), &(pbuf0[5]), &(pbuf1[5]), &(pbuf2[5]), &(pbuf0[6]), &(pbuf1[6]), &(pbuf2[6]), &(pbuf0[7]), &(pbuf1[7]), &(pbuf2[7]));
	  bufptr = next_item_mult(bufptr, 24);
	  uljj = 8;
	} else {
	  uljj = unfiltered_indiv_ct - unfiltered_indiv_ct8m;
	  for (ulii = 0; ulii < uljj; ulii++) {
	    sscanf(bufptr, "%lg %lg %lg", &(pbuf0[ulii]), &(pbuf1[ulii]), &(pbuf2[ulii]));
	    bufptr = next_item_mult(bufptr, 3);
	  }
	}
	for (ulii = 0; ulii < uljj; ulii++) {
	  if (!is_set(indiv_exclude, indiv_uidx + ulii)) {
	    pzero = pbuf0[ulii];
	    pone = pbuf1[ulii];
	    ptwo = pbuf2[ulii];
	    dxx = (pone + 2 * ptwo) * marker_wt;
	    g_dosage_vals[indiv_idx * BITCT + marker_idxl] = dxx;
	    if (is_missing_01) {
	      // IEEE 754 zero is actually zero bitmask
	      if (pzero + pone + ptwo == 0.0) {
#if __LP64__
		g_missing_dmasks[indiv_idx * BITCT + marker_idxl] = 0.0;
#else
                g_missing_dmasks[indiv_idx * BITCT + marker_idxl] = 0;
#endif
		if (distance_flat_missing) {
		  missing_tots[indiv_idx] += 1.0;
		} else {
		  set_bit_noct(cur_missings, indiv_idx);
		}
		g_missing_vals[indiv_idx] |= 1 << marker_idxl;
	      } else {
#if __LP64__
		*((unsigned long*)(&g_missing_dmasks[indiv_idx * BITCT + marker_idxl])) = 0x7fffffffffffffffLU;
#else
                g_missing_dmasks[indiv_idx * BITCT + marker_idxl] = 1;
#endif
		if (!distance_flat_missing) {
		  // defer missing_tots[indiv_idx] update until we know weight
		  *cmf_ptr++ = dxx;
		}
	      }
	    } else {
	      // TBD
	    }
	    indiv_idx++;
	  }
	}
	indiv_uidx += uljj;
      }
      if (!distance_flat_missing) {
	non_missing_ct = (unsigned int)(cmf_ptr - cur_marker_freqs);
	if (non_missing_ct) {
#ifdef __cplusplus
	  std::sort(cur_marker_freqs, cmf_ptr);
#else
	  qsort(cur_marker_freqs, non_missing_ct, sizeof(double), double_cmp);
#endif
	  if (is_missing_01) {
	    dxx = 0.0;
	    for (ulii = 1; ulii < non_missing_ct; ulii++) {
	      dxx += ulii * (non_missing_ct - ulii) * (cur_marker_freqs[ulii] - cur_marker_freqs[ulii - 1]);
	    }
	    dxx *= 2.0 / ((double)((unsigned long long)non_missing_ct * non_missing_ct));
	    g_missing_wts[marker_idxl] += dxx;
	    tot_missing_wt += dxx;
	    uljj = g_indiv_ct - non_missing_ct;
	    indiv_idx = 0;
	    for (ulii = 0; ulii < uljj; ulii++) {
	      indiv_idx = next_set_unsafe(cur_missings, indiv_idx);
	      missing_tots[indiv_idx++] += dxx;
	    }
	  } else {
	    // TBD
	  }
	}
      }
    }
    marker_uidx++;
    if (marker_idxl == BITCT - 1) {
      retval = update_distance_dosage_matrix(is_missing_01, distance_3d, distance_flat_missing, thread_ct);
      if (retval) {
	goto oxford_distance_calc_ret_1;
      }
      printf("\r%u markers complete.", marker_uidx);
      fflush(stdout);
    }
  }
  if (!feof(genfile)) {
    return RET_READ_FAIL;
  }
  marker_idxl = marker_uidx % BITCT;
  if (marker_idxl) {
    ulii = BITCT - marker_idxl;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      fill_double_zero(&(g_dosage_vals[indiv_idx * BITCT + marker_idxl]), ulii);
    }
    retval = update_distance_dosage_matrix(is_missing_01, distance_3d, distance_flat_missing, thread_ct);
    if (retval) {
      goto oxford_distance_calc_ret_1;
    }
  }
  tstc = g_thread_start[thread_ct];
  if (is_missing_01) {
    dptr = g_distance_matrix;
    dptr2 = g_distance_wt_matrix;
    for (indiv_idx = g_thread_start[0]; indiv_idx < tstc; indiv_idx++) {
      dptr3 = missing_tots;
      dxx = tot_missing_wt - missing_tots[indiv_idx];
      for (uljj = 0; uljj < indiv_idx; uljj++) {
	*dptr *= tot_missing_wt / (dxx - (*dptr3++) + (*dptr2++));
	dptr++;
      }
    }
  }
  putchar('\r');
  logprint("Distance matrix calculation complete.\n");
  retval = 0;
  // while (0) {
    // oxford_distance_calc_ret_INVALID_FORMAT:
    // printf("Error: Improperly formatted .gen file.\n");
    // retval = RET_INVALID_FORMAT;
    // break;
  // }
 oxford_distance_calc_ret_1:
  wkspace_reset(wkspace_mark);
  return retval;
}

int oxford_distance_calc_unscanned(FILE* genfile, unsigned int* gen_buf_len_ptr, double** set_allele_freqs_ptr, unsigned int* unfiltered_marker_ct_ptr, unsigned long** marker_exclude_ptr, unsigned int* marker_ct_ptr, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, int distance_3d, int distance_flat_missing, double exponent, unsigned int thread_ct, int parallel_idx, unsigned int parallel_tot) {
  // Easily usable when no filters are applied, or no .freq/.freqx file is
  // loaded and there are no filters on individuals that depend on genotype
  // information.
  // Could be extended to handle the --read-freq case as well, but that's more
  // complicated.
  unsigned char* wkspace_mark = NULL;
  int is_exponent_zero = (exponent == 0.0);
  unsigned int unfiltered_marker_ct = 0;
  unsigned int unfiltered_indiv_ct8m = unfiltered_indiv_ct & 0xfffffff8U;
  unsigned int marker_idxl = 0;
  unsigned int marker_ct = 0;
  unsigned int gen_buf_len = 0;
  double tot_missing_wt = 0.0;
  double missing_wt = 1.0;
  double* cur_marker_freqs = NULL;
  double sum_zero_3d = 0.0;
  double sum_one_3d = 0.0;
  double sum_two_3d = 0.0;
  double sum_exclude_3d = 0.0;
  double* set_allele_freqs_tmp;
  unsigned int unfiltered_marker_ctl;
  char* loadbuf;
  int max_load;
  char* bufptr;
  int retval;
  double dxx;
  double dyy;
  long long llxx;
  double pbuf0[8];
  double pbuf1[8];
  double pbuf2[8];
  double pzero;
  double pone;
  double ptwo;
  unsigned long ulii;
  unsigned long uljj;
  unsigned int indiv_uidx;
  unsigned int indiv_idx;
  unsigned int uii;
  unsigned int subloop_end;
  double* cur_nonmissings;
  double marker_wt;
  double ref_freq_numer;
  double ref_freq_denom; // also 2x sum of cur_nonmissings

  triangle_fill(g_thread_start, g_indiv_ct, thread_ct, parallel_idx, parallel_tot, 1, 1);
  llxx = g_thread_start[thread_ct];
  llxx = ((llxx * (llxx - 1)) - (long long)g_thread_start[0] * (g_thread_start[0] - 1)) / 2;
#ifndef __LP64__
  if (llxx > 4294967295LL) {
    return RET_NOMEM;
  }
#endif
  ulii = (unsigned long)llxx;
  if (wkspace_alloc_d_checked(&g_distance_matrix, ulii * sizeof(double))) {
    return RET_NOMEM;
  }
  wkspace_mark = wkspace_base;
  if (wkspace_alloc_d_checked(&g_distance_wt_matrix, ulii * sizeof(double))) {
    return RET_NOMEM;
  }
  fill_double_zero(g_distance_matrix, ulii);
  fill_double_zero(g_distance_wt_matrix, ulii);

  if (wkspace_alloc_d_checked(&cur_nonmissings, g_indiv_ct * sizeof(double))) {
    return RET_NOMEM;
  }

  if (distance_3d) {
    if (wkspace_alloc_d_checked(&g_dosage_vals, g_indiv_ct * 4 * MULTIPLEX_DOSAGE_NM3D * sizeof(double))) {
      return RET_NOMEM;
    }
  } else {
    if (wkspace_alloc_d_checked(&g_nonmissing_vals, g_indiv_ct * MULTIPLEX_DOSAGE_NM * sizeof(double))) {
      return RET_NOMEM;
    }
    if (wkspace_alloc_d_checked(&g_dosage_vals, g_indiv_ct * MULTIPLEX_DOSAGE_NM * sizeof(double))) {
      return RET_NOMEM;
    }
    if (wkspace_alloc_d_checked(&cur_marker_freqs, g_indiv_ct * sizeof(double))) {
      return RET_NOMEM;
    }
  }

  if (wkspace_left > 2147483584LU) {
    max_load = 2147483584LU;
  } else {
    max_load = wkspace_left;
  }

  set_allele_freqs_tmp = (double*)wkspace_base;
  if (distance_3d) {
    loadbuf = (char*)(&(wkspace_base[MULTIPLEX_DOSAGE_NM3D * sizeof(double)]));
    max_load -= MULTIPLEX_DOSAGE_NM3D * sizeof(double);
  } else {
    loadbuf = (char*)(&(wkspace_base[MULTIPLEX_DOSAGE_NM * sizeof(double)]));
    max_load -= MULTIPLEX_DOSAGE_NM * sizeof(double);
  }
  if (max_load <= 0) {
    return RET_NOMEM;
  }
  loadbuf[max_load - 1] = ' ';
  while (fgets(loadbuf, max_load, genfile) != NULL) {
    if (is_eoln(*loadbuf)) {
      continue;
    }
    if (!loadbuf[max_load - 1]) {
      logprint("Extremely long line found in .gen file.\n");
      return RET_NOMEM;
    }
    bufptr = next_item_mult(skip_initial_spaces(loadbuf), 4);
    indiv_idx = 0;
    ref_freq_numer = 0.0;
    ref_freq_denom = 0.0;
    if (distance_3d) {
      sum_zero_3d = 0.0;
      sum_one_3d = 0.0;
      sum_two_3d = 0.0;
      sum_exclude_3d = 0.0;
      for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct;) {
	if (no_more_items(bufptr)) {
	  goto oxford_distance_calc_unscanned_ret_INVALID_FORMAT;
	}
	if (indiv_uidx < unfiltered_indiv_ct8m) {
	  bufptr = next_item(bufptr);
	  if (sscanf(bufptr, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &(pbuf0[0]), &(pbuf1[0]), &(pbuf2[0]), &(pbuf0[1]), &(pbuf1[1]), &(pbuf2[1]), &(pbuf0[2]), &(pbuf1[2]), &(pbuf2[2]), &(pbuf0[3]), &(pbuf1[3]), &(pbuf2[3]), &(pbuf0[4]), &(pbuf1[4]), &(pbuf2[4]), &(pbuf0[5]), &(pbuf1[5]), &(pbuf2[5]), &(pbuf0[6]), &(pbuf1[6]), &(pbuf2[6]), &(pbuf0[7]), &(pbuf1[7]), &(pbuf2[7])) != 24) {
	    goto oxford_distance_calc_unscanned_ret_INVALID_FORMAT;
	  }
	  bufptr = next_item_mult(bufptr, 23);
	  subloop_end = 8;
	} else {
	  subloop_end = unfiltered_indiv_ct - unfiltered_indiv_ct8m;
	  for (uii = 0; uii < subloop_end; uii++) {
	    bufptr = next_item(bufptr);
	    if (sscanf(bufptr, "%lg %lg %lg", &(pbuf0[uii]), &(pbuf1[uii]), &(pbuf2[uii])) != 3) {
	      goto oxford_distance_calc_unscanned_ret_INVALID_FORMAT;
	    }
	    bufptr = next_item_mult(bufptr, 2);
	  }
	}
	for (uii = 0; uii < subloop_end; uii++) {
	  if (!is_set(indiv_exclude, indiv_uidx + uii)) {
	    pzero = pbuf0[uii];
	    pone = pbuf1[uii];
	    ptwo = pbuf2[uii];
	    dxx = pone + 2 * ptwo;
	    dyy = pzero + pone + ptwo;
	    cur_nonmissings[indiv_idx] = dyy;
	    sum_zero_3d += pzero;
	    g_dosage_vals[indiv_idx * 4 * MULTIPLEX_DOSAGE_NM3D + marker_idxl] = pzero;
	    sum_one_3d += pone;
	    g_dosage_vals[(indiv_idx * 4 + 1) * MULTIPLEX_DOSAGE_NM3D + marker_idxl] = pone;
	    sum_two_3d += ptwo;
	    g_dosage_vals[(indiv_idx * 4 + 2) * MULTIPLEX_DOSAGE_NM3D + marker_idxl] = ptwo;
	    sum_exclude_3d += (pzero * dxx) + (pone * ptwo);
	    g_dosage_vals[(indiv_idx * 4 + 3) * MULTIPLEX_DOSAGE_NM3D + marker_idxl] = dyy;
	    ref_freq_numer += dxx;
	    ref_freq_denom += 2 * dyy;
	    indiv_idx++;
	  }
	}
	indiv_uidx += subloop_end;
      }
    } else {
      for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct;) {
	if (no_more_items(bufptr)) {
	  goto oxford_distance_calc_unscanned_ret_INVALID_FORMAT;
	}
	if (indiv_uidx < unfiltered_indiv_ct8m) {
	  bufptr = next_item(bufptr);
	  if (sscanf(bufptr, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &(pbuf0[0]), &(pbuf1[0]), &(pbuf2[0]), &(pbuf0[1]), &(pbuf1[1]), &(pbuf2[1]), &(pbuf0[2]), &(pbuf1[2]), &(pbuf2[2]), &(pbuf0[3]), &(pbuf1[3]), &(pbuf2[3]), &(pbuf0[4]), &(pbuf1[4]), &(pbuf2[4]), &(pbuf0[5]), &(pbuf1[5]), &(pbuf2[5]), &(pbuf0[6]), &(pbuf1[6]), &(pbuf2[6]), &(pbuf0[7]), &(pbuf1[7]), &(pbuf2[7])) != 24) {
	    goto oxford_distance_calc_unscanned_ret_INVALID_FORMAT;
	  }
	  bufptr = next_item_mult(bufptr, 23);
	  subloop_end = 8;
	} else {
	  subloop_end = unfiltered_indiv_ct - unfiltered_indiv_ct8m;
	  for (uii = 0; uii < subloop_end; uii++) {
	    bufptr = next_item(bufptr);
	    if (sscanf(bufptr, "%lg %lg %lg", &(pbuf0[uii]), &(pbuf1[uii]), &(pbuf2[uii])) != 3) {
	      goto oxford_distance_calc_unscanned_ret_INVALID_FORMAT;
	    }
	    bufptr = next_item_mult(bufptr, 2);
	  }
	}
	for (uii = 0; uii < subloop_end; uii++) {
	  if (!is_set(indiv_exclude, indiv_uidx + uii)) {
	    pone = pbuf1[uii];
	    ptwo = pbuf2[uii];
	    dyy = pbuf0[uii] + pone + ptwo;
	    cur_nonmissings[indiv_idx] = dyy;
	    g_nonmissing_vals[indiv_idx * MULTIPLEX_DOSAGE_NM + marker_idxl] = dyy;
	    if (dyy == 0.0) {
	      dxx = 1.0;
	    } else {
	      dxx = pone + 2 * ptwo;
	      ref_freq_numer += dxx;
	      dxx = dxx / dyy;
	      ref_freq_denom += 2 * dyy;
	    }
	    g_dosage_vals[indiv_idx * MULTIPLEX_DOSAGE_NM + marker_idxl] = dxx;
	    cur_marker_freqs[indiv_idx] = dxx;
	    indiv_idx++;
	  }
	}
	indiv_uidx += subloop_end;
      }
    }
    uii = strlen(bufptr) + (unsigned int)(bufptr - loadbuf) + 1;
    if (uii > gen_buf_len) {
      gen_buf_len = uii;
    }

    if (!is_exponent_zero) {
      if ((ref_freq_numer == 0.0) || (ref_freq_numer == ref_freq_denom)) {
	marker_wt = 1.0;
	set_allele_freqs_tmp[unfiltered_marker_ct] = 0.5;
      } else {
        dxx = ref_freq_numer / ref_freq_denom;
	set_allele_freqs_tmp[unfiltered_marker_ct] = dxx;
        marker_wt = pow(2 * dxx * (1.0 - dxx), -exponent);
      }
    }

    if (!distance_flat_missing) {
      if (distance_3d) {
	missing_wt = 2 * (sum_zero_3d * (sum_one_3d + 2 * sum_two_3d) + (sum_one_3d * sum_two_3d) - sum_exclude_3d) / (((unsigned long)g_indiv_ct) * (g_indiv_ct - 1));
      } else {
	if (qsort_ext((char*)cur_marker_freqs, g_indiv_ct, sizeof(double), double_cmp_deref, (char*)cur_nonmissings, sizeof(double))) {
	  return RET_NOMEM;
	}
	dxx = 0.0;
	uii = g_indiv_ct - 1;
	dyy = ref_freq_denom * 0.5;
	missing_wt = 0.0;
	for (indiv_idx = 0; indiv_idx < uii; indiv_idx++) {
	  dxx += cur_nonmissings[indiv_idx];
	  missing_wt += (cur_marker_freqs[indiv_idx + 1] - cur_marker_freqs[indiv_idx]) * dxx * (dyy - dxx);
	}
	missing_wt *= 2.0 / (dyy * dyy);
      }
      g_missing_wts[marker_idxl] = missing_wt;
    }
    tot_missing_wt += missing_wt;

    unfiltered_marker_ct++;
    marker_idxl++;
    if ((distance_3d && (marker_idxl == MULTIPLEX_DOSAGE_NM3D)) || (marker_idxl == MULTIPLEX_DOSAGE_NM)) {
      if (distance_3d) {
        max_load -= MULTIPLEX_DOSAGE_NM3D * sizeof(double);
        loadbuf = (char*)(&(loadbuf[MULTIPLEX_DOSAGE_NM3D * sizeof(double)]));
      } else {
        max_load -= MULTIPLEX_DOSAGE_NM * sizeof(double);
        loadbuf = (char*)(&(loadbuf[MULTIPLEX_DOSAGE_NM * sizeof(double)]));
      }
      if (max_load <= 0) {
	return RET_NOMEM;
      }
      retval = update_distance_dosage_matrix(0, distance_3d, distance_flat_missing, thread_ct);
      if (retval) {
	return retval;
      }
      printf("\r%u markers complete.", unfiltered_marker_ct);
      fflush(stdout);
      marker_idxl = 0;
    }
  }
  if (!feof(genfile)) {
    return RET_READ_FAIL;
  }
  if (!unfiltered_marker_ct) {
    logprint("Error: No markers in .gen file.\n");
    return RET_INVALID_FORMAT;
  }
  if (marker_idxl) {
    if (distance_3d) {
      ulii = MULTIPLEX_DOSAGE_NM3D - marker_idxl;
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	fill_double_zero(&(g_dosage_vals[indiv_idx * MULTIPLEX_DOSAGE_NM3D + marker_idxl]), ulii);
      }
    } else {
      ulii = MULTIPLEX_DOSAGE_NM - marker_idxl;
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	fill_double_zero(&(g_dosage_vals[indiv_idx * MULTIPLEX_DOSAGE_NM + marker_idxl]), ulii);
	fill_double_zero(&(g_nonmissing_vals[indiv_idx * MULTIPLEX_DOSAGE_NM + marker_idxl]), ulii);
      }
    }
    retval = update_distance_dosage_matrix(0, distance_3d, distance_flat_missing, thread_ct);
    if (retval) {
      return retval;
    }
  }

  unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
  *marker_ct_ptr = marker_ct;
  *gen_buf_len_ptr = gen_buf_len;

  uljj = (((unsigned long)g_indiv_ct) * (g_indiv_ct - 1)) / 2;
  dxx = tot_missing_wt;
  for (ulii = 0; ulii < uljj; ulii++) {
    g_distance_matrix[ulii] *= dxx / g_distance_wt_matrix[ulii];
  }

  wkspace_reset(wkspace_mark);
  *set_allele_freqs_ptr = (double*)wkspace_alloc(unfiltered_marker_ct * sizeof(double));
  if (wkspace_alloc_ul_checked(marker_exclude_ptr, unfiltered_marker_ctl * sizeof(long))) {
    return RET_NOMEM;
  }
  fill_ulong_zero(*marker_exclude_ptr, unfiltered_marker_ctl);
  for (marker_idxl = 0; marker_idxl < unfiltered_marker_ct; marker_idxl++) {
    dxx = set_allele_freqs_tmp[marker_idxl];
    if (dxx == -1.0) {
      set_bit_sub(*marker_exclude_ptr, marker_idxl, marker_ct_ptr);
    } else {
      (*set_allele_freqs_ptr)[marker_idxl] = dxx;
    }
  }

  putchar('\r');
  logprint("Distance matrix calculation complete.\n");
  retval = 0;
  while (0) {
  oxford_distance_calc_unscanned_ret_INVALID_FORMAT:
    logprint("Error: Improperly formatted .gen file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  }
  return retval;
}

int wdist_dosage(int calculation_type, char* genname, char* samplename, char* outname, char* missing_code, int distance_3d, int distance_flat_missing, double exponent, int maf_succ, unsigned long regress_iters, unsigned int regress_d, unsigned int thread_ct, int parallel_idx, unsigned int parallel_tot) {
  FILE* genfile = NULL;
  FILE* outfile = NULL;
  FILE* outfile2 = NULL;
  FILE* outfile3 = NULL;
  gzFile gz_outfile = NULL;
  gzFile gz_outfile2 = NULL;
  gzFile gz_outfile3 = NULL;
  char* outname_end = (char*)memchr(outname, 0, FNAMESIZE);
  int gen_scanned = 0;
  int is_missing_01 = 0;
  unsigned char* sex_info;
  unsigned char* membuf;
  unsigned int gen_buf_len;
  unsigned char* wkspace_mark;
  char* pheno_c;
  double* pheno_d;
  unsigned long* pheno_exclude;
  double* set_allele_freqs;
  unsigned int unfiltered_marker_ct;
  unsigned int unfiltered_marker_ctl;
  unsigned long* marker_exclude;
  unsigned int marker_ct;
  unsigned int unfiltered_indiv_ct;
  unsigned long* indiv_exclude;
  char* person_ids;
  unsigned int max_person_id_len;
  int retval;
  unsigned int marker_uidx;
  unsigned int marker_idx;
  double dxx;
  double dyy;
  retval = oxford_sample_load(samplename, &unfiltered_indiv_ct, &person_ids, &max_person_id_len, &sex_info, &pheno_c, &pheno_d, &pheno_exclude, &indiv_exclude, missing_code);
  if (retval) {
    goto wdist_dosage_ret_1;
  }
  if ((calculation_type & CALC_REGRESS_DISTANCE) && (!pheno_d)) {
    logprint("Error: --regress-distance calculation requires scalar phenotype.\n");
    goto wdist_dosage_ret_INVALID_CMDLINE;
  }
  if (fopen_checked(&genfile, genname, "r")) {
    return RET_OPEN_FAIL;
  }
  if (0) { // enable this if filters are requested
    retval = oxford_gen_load1(genfile, &gen_buf_len, &unfiltered_marker_ct, &set_allele_freqs, &is_missing_01, unfiltered_indiv_ct, maf_succ);
    if (retval) {
      goto wdist_dosage_ret_1;
    }
    unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
    marker_ct = unfiltered_marker_ct;
    if (wkspace_alloc_ul_checked(&marker_exclude, unfiltered_marker_ctl * sizeof(long))) {
      goto wdist_dosage_ret_NOMEM;
    }
    fill_ulong_zero(marker_exclude, unfiltered_marker_ctl);
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      if (set_allele_freqs[marker_uidx] == -1.0) {
        set_bit_noct(marker_exclude, marker_uidx);
        marker_ct--;
      }
    }
    gen_scanned = 1;
  }
  g_indiv_ct = unfiltered_indiv_ct;
  if (parallel_tot > g_indiv_ct / 2) {
    sprintf(logbuf, "Error: Too many --parallel jobs (maximum %d/2 = %d).\n", g_indiv_ct, g_indiv_ct / 2);
    logprintb();
    goto wdist_dosage_ret_INVALID_CMDLINE;
  }
  if (thread_ct > 1) {
    sprintf(logbuf, "Using %d threads (change this with --threads).\n", thread_ct);
    logprintb();
  }
  if (distance_req(calculation_type)) {
    wkspace_mark = wkspace_base;
    if (gen_scanned) {
      retval = oxford_distance_calc(genfile, gen_buf_len, set_allele_freqs, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, is_missing_01, distance_3d, distance_flat_missing, exponent, thread_ct, parallel_idx, parallel_tot);
    } else {
      // N.B. this causes set_allele_freqs and marker_exclude to be allocated
      // *above* distance_matrix on the stack.
      retval = oxford_distance_calc_unscanned(genfile, &gen_buf_len, &set_allele_freqs, &unfiltered_marker_ct, &marker_exclude, &marker_ct, unfiltered_indiv_ct, indiv_exclude, distance_3d, distance_flat_missing, exponent, thread_ct, parallel_idx, parallel_tot);
    }
    if (retval) {
      goto wdist_dosage_ret_1;
    }
    if (wkspace_alloc_uc_checked(&membuf, g_indiv_ct * sizeof(double))) {
      goto wdist_dosage_ret_NOMEM;
    }
    if (calculation_type & CALC_DISTANCE_MASK) {
      retval = distance_d_write_ids(outname, outname_end, calculation_type, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
      if (retval) {
	goto wdist_dosage_ret_1;
      }
      if ((exponent == 0.0) || (!(calculation_type & (CALC_DISTANCE_IBS | CALC_DISTANCE_1_MINUS_IBS)))) {
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
	    dxx += 1.0; // prevent IBS measure from breaking down
	  }
	  marker_uidx++;
	}
	dxx = 0.5 / dxx;
      }
      retval = distance_d_write(&outfile, &outfile2, &outfile3, &gz_outfile, &gz_outfile2, &gz_outfile3, calculation_type, outname, outname_end, g_distance_matrix, dxx, g_indiv_ct, g_thread_start[0], g_thread_start[thread_ct], parallel_idx, parallel_tot, membuf);
      if (retval) {
        goto wdist_dosage_ret_1;
      }
    }
    if (calculation_type & CALC_REGRESS_DISTANCE) {
      retval = regress_distance(calculation_type, g_distance_matrix, pheno_d, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, thread_ct, regress_iters, regress_d);
      if (retval) {
	goto wdist_dosage_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }
  while (0) {
  wdist_dosage_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  wdist_dosage_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  wdist_dosage_ret_1:
  fclose_cond(genfile);
  return retval;
}
