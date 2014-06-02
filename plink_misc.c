#include "plink_common.h"

#include "plink_misc.h"

void misc_init(Score_info* sc_ip) {
  sc_ip->fname = NULL;
  sc_ip->range_fname = NULL;
  sc_ip->data_fname = NULL;
  sc_ip->modifier = 0;
  sc_ip->varid_col = 1;
  sc_ip->allele_col = 2;
  sc_ip->effect_col = 3;
  sc_ip->data_varid_col = 1;
  sc_ip->data_col = 2;
}

void misc_cleanup(Score_info* sc_ip) {
  free_cond(sc_ip->fname);
  free_cond(sc_ip->range_fname);
  free_cond(sc_ip->data_fname);
}

int32_t make_founders(uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t require_two, uintptr_t* indiv_exclude, uintptr_t* founder_info) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uint32_t new_founder_ct = 0;
  int32_t retval = 0;
  char* sorted_ids;
  char* id_buf;
  char* wptr;
  char* pat_ptr;
  char* mat_ptr;
  uintptr_t* nf_bitarr;
  uintptr_t indiv_uidx;
  uint32_t fam_len_p1;
  uint32_t missing_parent_ct;
  uint32_t cur_len;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&nf_bitarr, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto make_founders_ret_NOMEM;
  }
  bitfield_exclude_to_include(indiv_exclude, nf_bitarr, unfiltered_indiv_ct);
  bitfield_andnot(nf_bitarr, founder_info, unfiltered_indiv_ctl);
  indiv_uidx = next_set(nf_bitarr, 0, unfiltered_indiv_ct);
  if (indiv_uidx == unfiltered_indiv_ct) {
    logprint("Note: Skipping --make-founders since there are no nonfounders.\n");
    goto make_founders_ret_1;
  }
  sorted_ids = alloc_and_init_collapsed_arr(person_ids, max_person_id_len, unfiltered_indiv_ct, indiv_exclude, indiv_ct, 0);
  if (!sorted_ids) {
    goto make_founders_ret_NOMEM;
  }
  qsort(sorted_ids, indiv_ct, max_person_id_len, strcmp_casted);
  do {
    pat_ptr = &(person_ids[indiv_uidx * max_person_id_len]);
    fam_len_p1 = strlen_se(pat_ptr) + 1;
    wptr = memcpya(id_buf, pat_ptr, fam_len_p1);
    missing_parent_ct = 0;
    pat_ptr = &(paternal_ids[indiv_uidx * max_paternal_id_len]);
    cur_len = strlen(pat_ptr);
    if (cur_len + fam_len_p1 >= max_person_id_len) {
      missing_parent_ct++;
    } else {
      memcpy(wptr, pat_ptr, cur_len);
      if (bsearch_str(id_buf, cur_len + fam_len_p1, sorted_ids, max_person_id_len, indiv_ct) == -1) {
	missing_parent_ct++;
      }
    }
    mat_ptr = &(maternal_ids[indiv_uidx * max_maternal_id_len]);
    cur_len = strlen(mat_ptr);
    if (cur_len + fam_len_p1 >= max_person_id_len) {
      missing_parent_ct++;
    } else {
      memcpy(wptr, mat_ptr, cur_len);
      if (bsearch_str(id_buf, cur_len + fam_len_p1, sorted_ids, max_person_id_len, indiv_ct) == -1) {
	missing_parent_ct++;
      }
    }
    if (missing_parent_ct > require_two) {
      SET_BIT(founder_info, indiv_uidx);
      memcpy(pat_ptr, "0", 2);
      memcpy(mat_ptr, "0", 2);
      new_founder_ct++;
    }
    indiv_uidx++;
    next_set_ul_ck(nf_bitarr, &indiv_uidx, unfiltered_indiv_ct);
  } while (indiv_uidx < unfiltered_indiv_ct);
  LOGPRINTF("--make-founders: %u individual%s affected.\n", new_founder_ct, (new_founder_ct == 1)? "" : "s");
  while (0) {
  make_founders_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
 make_founders_ret_1:
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t write_nosex(char* outname, char* outname_end, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* sex_nm, uintptr_t gender_unk_ct, char* person_ids, uintptr_t max_person_id_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_uidx = 0;
  int32_t retval = 0;
  uintptr_t* sex_missing;
  uintptr_t indiv_idx;
  if (wkspace_alloc_ul_checked(&sex_missing, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto write_nosex_ret_NOMEM;
  }
  bitfield_exclude_to_include(indiv_exclude, sex_missing, unfiltered_indiv_ct);
  bitfield_andnot(sex_missing, sex_nm, unfiltered_indiv_ctl);
  memcpy(outname_end, ".nosex", 7);
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_nosex_ret_OPEN_FAIL;
  }
  for (indiv_idx = 0; indiv_idx < gender_unk_ct; indiv_idx++, indiv_uidx++) {
    next_set_ul_unsafe_ck(sex_missing, &indiv_uidx);
    fputs(&(person_ids[indiv_uidx * max_person_id_len]), outfile);
    putc('\n', outfile);
  }
  if (fclose_null(&outfile)) {
    goto write_nosex_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("Ambiguous sex ID%s written to %s .\n", (gender_unk_ct == 1)? "" : "s", outname);
  while (0) {
  write_nosex_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_nosex_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_nosex_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t makepheno_load(FILE* phenofile, char* makepheno_str, uintptr_t unfiltered_indiv_ct, char* sorted_person_ids, uintptr_t max_person_id_len, uint32_t* id_map, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr) {
  uint32_t mp_strlen = strlen(makepheno_str);
  uint32_t makepheno_all = ((mp_strlen == 1) && (makepheno_str[0] == '*'));
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t* pheno_c = *pheno_c_ptr;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t line_idx = 0;
  int32_t retval = 0;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  int32_t ii;
  uint32_t person_idx;
  uint32_t tmp_len;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    goto makepheno_load_ret_NOMEM;
  }
  if (!pheno_c) {
    if (aligned_malloc(pheno_c_ptr, unfiltered_indiv_ctl * sizeof(intptr_t))) {
      goto makepheno_load_ret_NOMEM;
    }
    pheno_c = *pheno_c_ptr;
    fill_ulong_zero(pheno_c, unfiltered_indiv_ctl);
  }
  if (makepheno_all) {
    fill_all_bits(pheno_nm, unfiltered_indiv_ct);
  }
  // probably want to permit long lines here
  tbuf[MAXLINELEN - 1] = ' '; 
  while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Line %" PRIuPTR " of --make-pheno file is pathologically long.\n", line_idx);
      goto makepheno_load_ret_INVALID_FORMAT_2;
    }
    bufptr0 = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr0)) {
      continue;
    }
    if (bsearch_read_fam_indiv(id_buf, sorted_person_ids, max_person_id_len, unfiltered_indiv_ct, bufptr0, &bufptr, &ii)) {
      sprintf(logbuf, "Error: Line %" PRIuPTR " of --make-pheno file has fewer tokens than expected.\n", line_idx);
      goto makepheno_load_ret_INVALID_FORMAT_2;
    }
    if (ii != -1) {
      person_idx = id_map[(uint32_t)ii];
      if (makepheno_all) {
	SET_BIT(pheno_c, person_idx);
      } else {
	SET_BIT(pheno_nm, person_idx);
        tmp_len = strlen_se(bufptr);
	if ((tmp_len == mp_strlen) && (!memcmp(bufptr, makepheno_str, mp_strlen))) {
	  SET_BIT(pheno_c, person_idx);
	}
      }
    }
  }
  if (!feof(phenofile)) {
    goto makepheno_load_ret_READ_FAIL;
  }
  tmp_len = popcount_longs(pheno_nm, unfiltered_indiv_ctl);
  LOGPRINTF("--make-pheno: %u phenotype value%s set.\n", tmp_len, (tmp_len == 1)? "" : "s");
  while (0) {
  makepheno_load_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  makepheno_load_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  makepheno_load_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t load_pheno(FILE* phenofile, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_exclude_ct, char* sorted_person_ids, uintptr_t max_person_id_len, uint32_t* id_map, int32_t missing_pheno, uint32_t affection_01, uint32_t mpheno_col, char* phenoname_str, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, char* phenoname_load, uintptr_t max_pheno_name_len) {
  uint32_t affection = 1;
  uintptr_t* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  int32_t header_processed = 0;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  uintptr_t line_idx = 0;
  char case_char = affection_01? '1' : '2';
  uintptr_t* isz = NULL;
  double missing_phenod = (double)missing_pheno;
  int32_t retval = 0;
  char* loadbuf;
  uint32_t loadbuf_size;
  int32_t person_idx;
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
      goto load_pheno_ret_NOMEM;
    }
    fill_ulong_zero(isz, unfiltered_indiv_ctl);
    if (!pheno_c) {
      if (aligned_malloc(pheno_c_ptr, unfiltered_indiv_ctl * sizeof(intptr_t))) {
	goto load_pheno_ret_NOMEM;
      }
      pheno_c = *pheno_c_ptr;
      fill_ulong_zero(pheno_c, unfiltered_indiv_ctl);
    }
  }
  // ----- phenotype file load -----
  // worthwhile to support very long lines here...
  if (wkspace_left > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (wkspace_left > MAXLINELEN) {
    loadbuf_size = wkspace_left;
  } else {
    goto load_pheno_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, phenofile) != NULL) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --pheno file is pathologically long.\n", line_idx);
	goto load_pheno_ret_INVALID_FORMAT_2;
      } else {
	goto load_pheno_ret_NOMEM;
      }
    }
    bufptr0 = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr0)) {
      continue;
    }
    if (!header_processed) {
      tmp_len = strlen_se(bufptr0);
      bufptr = skip_initial_spaces(&(bufptr0[tmp_len]));
      if (is_eoln_kns(*bufptr)) {
	goto load_pheno_ret_MISSING_TOKENS;
      }
      tmp_len2 = strlen_se(bufptr);
      if ((tmp_len == 3) && (tmp_len2 == 3) && (!memcmp("FID", bufptr0, 3)) && (!memcmp("IID", bufptr, 3))) {
	if (phenoname_str) {
	  tmp_len = strlen(phenoname_str);
	  do {
	    bufptr = next_token(bufptr);
	    if (no_more_tokens_kns(bufptr)) {
	      logprint("Error: --pheno-name column not found.\n");
              goto load_pheno_ret_INVALID_FORMAT;
	    }
	    mpheno_col++;
	    tmp_len2 = strlen_se(bufptr);
	  } while ((tmp_len2 != tmp_len) || memcmp(bufptr, phenoname_str, tmp_len));
	} else if (phenoname_load) {
          bufptr = next_token_mult(bufptr, mpheno_col);
	  if (no_more_tokens_kns(bufptr)) {
	    return LOAD_PHENO_LAST_COL;
	  }
	  tmp_len = strlen_se(bufptr);
	  if (tmp_len > max_pheno_name_len) {
	    logprint("Error: Excessively long phenotype name in --pheno file.\n");
            goto load_pheno_ret_INVALID_FORMAT;
	  }
          memcpyx(phenoname_load, bufptr, tmp_len, '\0');
	}
      } else if (phenoname_str) {
	logprint("Error: --pheno-name requires the --pheno file to have a header line with first\ntwo columns 'FID' and 'IID'.\n");
	goto load_pheno_ret_INVALID_FORMAT;
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
      if (bsearch_read_fam_indiv(tbuf, sorted_person_ids, max_person_id_len, indiv_ct, bufptr0, &bufptr, &person_idx)) {
	goto load_pheno_ret_MISSING_TOKENS;
      }
      if (person_idx != -1) {
	person_idx = id_map[(uint32_t)person_idx];
	if (mpheno_col > 1) {
	  bufptr = next_token_mult(bufptr, mpheno_col - 1);
	}
	if (no_more_tokens_kns(bufptr)) {
	  // Sometimes, but not always, an error.  So we populate logbuf but
	  // let the caller decide whether to actually log it.
          sprintf(logbuf, "Error: Line %" PRIuPTR " of --pheno file has fewer tokens than expected.\n", line_idx);
	  return LOAD_PHENO_LAST_COL;
	}
	if (affection) {
	  if (eval_affection(bufptr, missing_pheno, affection_01)) {
	    if (is_missing_pheno(bufptr, missing_pheno, affection_01)) {
	      // Since we're only making one pass through the file, we don't
	      // have the luxury of knowing in advance whether the phenotype is
	      // binary or scalar.  If there is a '0' entry that occurs before
	      // we know the phenotype is scalar, we need to not set the
	      // phenotype to zero during the binary -> scalar conversion step.
	      if (*bufptr == '0') {
		set_bit(isz, person_idx);
	      }
	      clear_bit(pheno_c, person_idx);
	    } else {
	      if (*bufptr == case_char) {
		set_bit(pheno_c, person_idx);
	      } else {
		clear_bit(pheno_c, person_idx);
	      }
	      set_bit(pheno_nm, person_idx);
	    }
	  } else {
	    pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
	    if (!pheno_d) {
	      goto load_pheno_ret_NOMEM;
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
		set_bit(pheno_nm, uii);
	      } else if (is_set(pheno_nm, uii)) {
		pheno_d[uii] = is_set(pheno_c, uii)? dyy : dxx;
	      }
	    }
	    aligned_free_null(pheno_c_ptr);
	    affection = 0;
	  }
	}
	if (!affection) {
	  if ((!scan_double(bufptr, &dxx)) && (dxx != missing_phenod)) {
	    pheno_d[(uint32_t)person_idx] = dxx;
	    set_bit(pheno_nm, person_idx);
	  }
	}
      }
    }
  }
  if (!feof(phenofile)) {
    goto load_pheno_ret_READ_FAIL;
  }
  uii = popcount_longs(pheno_nm, unfiltered_indiv_ctl);
  LOGPRINTF("%u phenotype value%s present after --pheno.\n", uii, (uii == 1)? "" : "s");
  while (0) {
  load_pheno_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_pheno_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_pheno_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --pheno file has fewer tokens than expected.\n", line_idx);
  load_pheno_ret_INVALID_FORMAT_2:
    logprintb();
  load_pheno_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t convert_tail_pheno(uint32_t unfiltered_indiv_ct, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, double tail_bottom, double tail_top, double missing_phenod) {
  uintptr_t* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  uint32_t indiv_uidx;
  uint32_t indiv_uidx_stop;
  double dxx;
  if (!(*pheno_d_ptr)) {
    logprint("Error: --tail-pheno requires scalar phenotype data.\n");
    return RET_INVALID_FORMAT;
  }
  indiv_uidx = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  if (!pheno_c) {
    if (aligned_malloc(pheno_c_ptr, indiv_uidx * sizeof(intptr_t))) {
      return RET_NOMEM;
    }
    pheno_c = *pheno_c_ptr;
  }
  fill_ulong_zero(pheno_c, indiv_uidx);
  indiv_uidx = 0;
  do {
    indiv_uidx = next_set(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
    indiv_uidx_stop = next_unset(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
    for (; indiv_uidx < indiv_uidx_stop; indiv_uidx++) {
      dxx = pheno_d[indiv_uidx];
      if (dxx > tail_bottom) {
        if (dxx > tail_top) {
          SET_BIT(pheno_c, indiv_uidx);
        } else {
	  CLEAR_BIT(pheno_nm, indiv_uidx);
        }
      }
    }
  } while (indiv_uidx_stop < unfiltered_indiv_ct);
  free(pheno_d);
  *pheno_d_ptr = NULL;
  indiv_uidx = popcount_longs(pheno_nm, (unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  LOGPRINTF("--tail-pheno: %u phenotype value%s remaining.\n", indiv_uidx, (indiv_uidx == 1)? "" : "s");
  return 0;
}

int32_t apply_cm_map(char* cm_map_fname, char* cm_map_chrname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, double* marker_cms, Chrom_info* chrom_info_ptr) {
  FILE* shapeitfile = NULL;
  char* at_sign_ptr = NULL;
  char* fname_write = NULL;
  char* fname_buf = &(tbuf[MAXLINELEN]);
  double cm_old = 0.0;
  uint32_t autosome_ct = chrom_info_ptr->autosome_ct;
  uint32_t post_at_sign_len = 0;
  uint32_t updated_chrom_ct = 0;
  int32_t retval = 0;
  char* bufptr;
  char* bufptr2;
  double cm_new;
  double dxx;
  uintptr_t line_idx;
  uintptr_t irreg_line_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_ct;
  uint32_t chrom_end;
  uint32_t marker_uidx;
  uint32_t uii;
  int32_t bp_old;
  int32_t bp_new;
  int32_t ii;
  if (!cm_map_chrname) {
    chrom_fo_idx = 0;
    chrom_ct = chrom_info_ptr->chrom_ct;
    at_sign_ptr = strchr(cm_map_fname, '@');
    fname_write = memcpya(fname_buf, cm_map_fname, (uintptr_t)(at_sign_ptr - cm_map_fname));
    at_sign_ptr++;
    post_at_sign_len = strlen(at_sign_ptr) + 1;
  } else {
    ii = get_chrom_code(chrom_info_ptr, cm_map_chrname);
    if (ii < 0) {
      LOGPREPRINTFWW("Error: --cm-map chromosome code '%s' not found in dataset.\n", cm_map_chrname);
      goto apply_cm_map_ret_INVALID_CMDLINE_2;
    }
    chrom_fo_idx = get_marker_chrom_fo_idx(chrom_info_ptr, chrom_info_ptr->chrom_start[(uint32_t)ii]);
    chrom_ct = chrom_fo_idx + 1;
    fname_buf = cm_map_fname;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  for (; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    if (marker_uidx == chrom_end) {
      continue;
    }
    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    if (!cm_map_chrname) {
      if ((!uii) || (uii > autosome_ct)) {
        continue;
      }
      bufptr = uint32_write(fname_write, uii);
      memcpy(bufptr, at_sign_ptr, post_at_sign_len);
      if (fopen_checked(&shapeitfile, fname_buf, "r")) {
	LOGPRINTFWW("Warning: --cm-map failed to open %s.\n", fname_buf);
        continue;
      }
    } else {
      if (fopen_checked(&shapeitfile, cm_map_fname, "r")) {
        goto apply_cm_map_ret_OPEN_FAIL;
      }
    }
    updated_chrom_ct++;
    irreg_line_ct = 0;
    // First line is a header with three arbitrary fields.
    // All subsequent lines have three fields in the following order:
    //   1. bp position (increasing)
    //   2. cM/Mb recombination rate between current and previous bp positions
    //   3. current cM position
    // We mostly ignore field 2, since depending just on fields 1 and 3
    // maximizes accuracy.  The one exception is the very first nonheader line.
    retval = load_to_first_token(shapeitfile, MAXLINELEN, '\0', "--cm-map file", tbuf, &bufptr, &line_idx);
    if (retval) {
      goto apply_cm_map_ret_1;
    }
    bufptr = next_token_mult(bufptr, 2);
    if (no_more_tokens_kns(bufptr)) {
      goto apply_cm_map_ret_MISSING_TOKENS;
    }
    bufptr = next_token(bufptr);
    if (!no_more_tokens_kns(bufptr)) {
      goto apply_cm_map_ret_MISSING_TOKENS;
    }
    bp_old = -1;
    while (fgets(tbuf, MAXLINELEN, shapeitfile)) {
      line_idx++;
      if (!tbuf[MAXLINELEN - 1]) {
        sprintf(logbuf, "Error: Line %" PRIuPTR " of --cm-map file is pathologically long.\n", line_idx);
        goto apply_cm_map_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(tbuf);
      if ((*bufptr < '+') || (*bufptr > '9')) {
	// warning instead of error if text line found, since as of 8 Jan 2014
        // the posted chromosome 19 map has such a line
	if (*bufptr > ' ') {
	  irreg_line_ct++;
	}
        continue;
      }
      if (scan_uint_defcap(bufptr, (uint32_t*)&bp_new)) {
	sprintf(logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of --cm-map file.\n", line_idx);
        goto apply_cm_map_ret_INVALID_FORMAT_2;
      }
      if (bp_new <= bp_old) {
        logprint("Error: bp coordinates in --cm-map file are not in increasing order.\n");
        goto apply_cm_map_ret_INVALID_FORMAT;
      }
      bufptr2 = next_token_mult(bufptr, 2);
      if (no_more_tokens_kns(bufptr2)) {
	goto apply_cm_map_ret_MISSING_TOKENS;
      }
      if (scan_double(bufptr2, &cm_new)) {
	sprintf(logbuf, "Error: Invalid centimorgan position on line %" PRIuPTR " of --cm-map file.\n", line_idx);
	goto apply_cm_map_ret_INVALID_FORMAT_2;
      }
      if (bp_old == -1) {
	// parse field 2 only in this case
        bufptr = next_token(bufptr);
        if (scan_double(bufptr, &dxx)) {
	  sprintf(logbuf, "Error: Invalid recombination rate on line %" PRIuPTR " of --cm-map file.\n", line_idx);
	  goto apply_cm_map_ret_INVALID_FORMAT_2;
	}
        cm_old = cm_new - dxx * 0.000001 * ((double)(bp_new + 1));
      }
      dxx = (cm_new - cm_old) / ((double)(bp_new - bp_old));
      while (marker_pos[marker_uidx] <= ((uint32_t)bp_new)) {
	marker_cms[marker_uidx] = cm_new - ((int32_t)(((uint32_t)bp_new) - marker_pos[marker_uidx])) * dxx;
	marker_uidx++;
	next_unset_ck(marker_exclude, &marker_uidx, chrom_end);
	if (marker_uidx == chrom_end) {
	  goto apply_cm_map_chrom_done;
	}
      }
      bp_old = bp_new;
      cm_old = cm_new;
    }
    for (; marker_uidx < chrom_end; marker_uidx++) {
      marker_cms[marker_uidx] = cm_old;
    }
  apply_cm_map_chrom_done:
    if (fclose_null(&shapeitfile)) {
      goto apply_cm_map_ret_READ_FAIL;
    }
    if (irreg_line_ct) {
      LOGPRINTFWW("Warning: %" PRIuPTR " irregular line%s skipped in %s.\n", irreg_line_ct, (irreg_line_ct == 1)? "" : "s", fname_buf);
    }
  }
  LOGPRINTF("--cm-map: %u chromosome%s updated.\n", updated_chrom_ct, (updated_chrom_ct == 1)? "" : "s");
  while (0) {
  apply_cm_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  apply_cm_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  apply_cm_map_ret_INVALID_CMDLINE_2:
    logprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  apply_cm_map_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --cm-map file has fewer tokens than expected.\n", line_idx);
  apply_cm_map_ret_INVALID_FORMAT_2:
    logprintb();
  apply_cm_map_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 apply_cm_map_ret_1:
  fclose_cond(shapeitfile);
  return retval;
}

int32_t update_marker_cms(Two_col_params* update_cm, char* sorted_marker_ids, uintptr_t sorted_ids_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, double* marker_cms) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = update_cm->skipchar;
  uint32_t colid_first = (update_cm->colid < update_cm->colx);
  uintptr_t sorted_ids_ctl = (sorted_ids_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uintptr_t slen;
  uintptr_t line_idx;
  uint32_t colmin;
  uint32_t coldiff;
  char* colid_ptr;
  char* colx_ptr;
  int32_t sorted_idx;
  uint32_t marker_uidx;
  char cc;
  int32_t retval;
  if (wkspace_alloc_ul_checked(&already_seen, sorted_ids_ctl * sizeof(intptr_t))) {
    goto update_marker_cms_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, sorted_ids_ctl);

  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_cms_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, update_cm->fname, loadbuf, loadbuf_size, update_cm->skip);
  if (retval) {
    goto update_marker_cms_ret_1;
  }
  if (colid_first) {
    colmin = update_cm->colid - 1;
    coldiff = update_cm->colx - update_cm->colid;
  } else {
    colmin = update_cm->colx - 1;
    coldiff = update_cm->colid - update_cm->colx;
  }
  line_idx = update_cm->skip;
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --update-cm file is pathologically long.\n", line_idx);
	goto update_marker_cms_ret_INVALID_FORMAT_2;
      } else {
        goto update_marker_cms_ret_NOMEM;
      }
    }
    colid_ptr = skip_initial_spaces(loadbuf);
    cc = *colid_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colid_first) {
      colid_ptr = next_token_multz(colid_ptr, colmin);
      colx_ptr = next_token_mult(colid_ptr, coldiff);
      if (no_more_tokens_kns(colx_ptr)) {
	goto update_marker_cms_ret_MISSING_TOKENS;
      }
    } else {
      colx_ptr = next_token_multz(colid_ptr, colmin);
      colid_ptr = next_token_mult(colx_ptr, coldiff);
      if (no_more_tokens_kns(colid_ptr)) {
	goto update_marker_cms_ret_MISSING_TOKENS;
      }
    }
    slen = strlen_se(colid_ptr);
    sorted_idx = bsearch_str(colid_ptr, slen, sorted_marker_ids, max_marker_id_len, sorted_ids_ct);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      colid_ptr[slen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant '%s' in --update-cm file.\n", colid_ptr);
      goto update_marker_cms_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[(uint32_t)sorted_idx];
    if (scan_double(colx_ptr, &(marker_cms[marker_uidx]))) {
      sprintf(logbuf, "Error: Invalid centimorgan position on line %" PRIuPTR " of --update-cm file.\n", line_idx);
      goto update_marker_cms_ret_INVALID_FORMAT_2;
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_cms_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-cm: %" PRIuPTR " value%s changed, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-cm: %" PRIuPTR " value%s changed.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  while (0) {
  update_marker_cms_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_cms_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_cms_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --update-cm file has fewer tokens than expected.\n", line_idx);
  update_marker_cms_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 update_marker_cms_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_marker_pos(Two_col_params* update_map, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t* marker_pos, uint32_t* map_is_unsorted_ptr, Chrom_info* chrom_info_ptr) {
  // requires sorted_ids_ct == marker_ct
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = update_map->skipchar;
  uint32_t colid_first = (update_map->colid < update_map->colx);
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uint32_t map_is_unsorted = ((*map_is_unsorted_ptr) & UNSORTED_CHROM);
  uintptr_t orig_exclude_ct = *marker_exclude_ct_ptr;
  uint32_t chrom_fo_idx_p1 = 0;
  uint32_t chrom_end = 0;
  uint32_t last_pos = 0;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uintptr_t line_idx;
  uintptr_t slen;
  uint32_t colmin;
  uint32_t coldiff;
  char* colid_ptr;
  char* colx_ptr;
  int32_t sorted_idx;
  uint32_t marker_uidx;
  uint32_t marker_idx;
  char cc;
  int32_t retval;
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto update_marker_pos_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);

  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_pos_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, update_map->fname, loadbuf, loadbuf_size, update_map->skip);
  if (retval) {
    goto update_marker_pos_ret_1;
  }
  if (colid_first) {
    colmin = update_map->colid - 1;
    coldiff = update_map->colx - update_map->colid;
  } else {
    colmin = update_map->colx - 1;
    coldiff = update_map->colid - update_map->colx;
  }
  line_idx = update_map->skip;
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        sprintf(logbuf, "Error: Line %" PRIuPTR " of --update-map file is pathologically long.\n", line_idx);
	goto update_marker_pos_ret_INVALID_FORMAT_2;
      } else {
        goto update_marker_pos_ret_NOMEM;
      }
    }
    colid_ptr = skip_initial_spaces(loadbuf);
    cc = *colid_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colid_first) {
      colid_ptr = next_token_multz(colid_ptr, colmin);
      colx_ptr = next_token_mult(colid_ptr, coldiff);
      if (no_more_tokens_kns(colx_ptr)) {
	goto update_marker_pos_ret_MISSING_TOKENS;
      }
    } else {
      colx_ptr = next_token_multz(colid_ptr, colmin);
      colid_ptr = next_token_mult(colx_ptr, coldiff);
      if (no_more_tokens_kns(colid_ptr)) {
	goto update_marker_pos_ret_MISSING_TOKENS;
      }
    }
    slen = strlen_se(colid_ptr);
    sorted_idx = bsearch_str(colid_ptr, slen, sorted_marker_ids, max_marker_id_len, marker_ct);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      colid_ptr[slen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant '%s' in --update-map file.\n", colid_ptr);
      goto update_marker_pos_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[(uint32_t)sorted_idx];
    if (scan_int_abs_defcap(colx_ptr, &sorted_idx)) {
      sprintf(logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of --update-map file.\n", line_idx);
      goto update_marker_pos_ret_INVALID_FORMAT_2;
    }
    if (sorted_idx < 0) {
      SET_BIT(marker_exclude, marker_uidx);
      *marker_exclude_ct_ptr += 1;
    } else {
      marker_pos[marker_uidx] = sorted_idx;
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_pos_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-map: %" PRIuPTR " value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-map: %" PRIuPTR " value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  marker_uidx = 0;
  marker_ct -= (*marker_exclude_ct_ptr) - orig_exclude_ct;
  if (!marker_ct) {
    logprint("Error: All variants excluded by --update-map (due to negative marker\npositions).\n");
    goto update_marker_pos_ret_ALL_MARKERS_EXCLUDED;
  }
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    while (marker_uidx >= chrom_end) {
      chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[++chrom_fo_idx_p1];
      last_pos = 0;
    }
    if (last_pos > marker_pos[marker_uidx]) {
      map_is_unsorted |= UNSORTED_BP;
      if (!((*map_is_unsorted_ptr) & UNSORTED_BP)) {
	logprint("Warning: Base-pair positions are now unsorted!\n");
      }
      break;
    }
    last_pos = marker_pos[marker_uidx];
  }
  if (((*map_is_unsorted_ptr) & UNSORTED_BP) && (!(map_is_unsorted & UNSORTED_BP))) {
    logprint("Base-pair positions are now sorted.\n");
  }
  *map_is_unsorted_ptr = map_is_unsorted;
  while (0) {
  update_marker_pos_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_pos_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_pos_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --update-map file has fewer tokens than expected.\n", line_idx);
  update_marker_pos_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  update_marker_pos_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 update_marker_pos_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_marker_names(Two_col_params* update_name, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char* true_marker_ids) {
  // requires sorted_ids_ct == marker_ct
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = update_name->skipchar;
  uint32_t colold_first = (update_name->colid < update_name->colx);
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uint32_t colmin;
  uint32_t coldiff;
  char* colold_ptr;
  char* colnew_ptr;
  int32_t sorted_idx;
  uintptr_t marker_uidx;
  uint32_t slen;
  char cc;
  int32_t retval;
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto update_marker_names_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);
  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_names_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, update_name->fname, loadbuf, loadbuf_size, update_name->skip);
  if (retval) {
    goto update_marker_names_ret_1;
  }
  if (colold_first) {
    colmin = update_name->colid - 1;
    coldiff = update_name->colx - update_name->colid;
  } else {
    colmin = update_name->colx - 1;
    coldiff = update_name->colid - update_name->colx;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    // no need for line_idx since all validation happened earlier
    if (!(loadbuf[loadbuf_size - 1])) {
      goto update_marker_names_ret_NOMEM;
    }
    colold_ptr = skip_initial_spaces(loadbuf);
    cc = *colold_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colold_first) {
      colold_ptr = next_token_multz(colold_ptr, colmin);
      colnew_ptr = next_token_mult(colold_ptr, coldiff);
    } else {
      colnew_ptr = next_token_multz(colold_ptr, colmin);
      colold_ptr = next_token_mult(colnew_ptr, coldiff);
    }
    slen = strlen_se(colold_ptr);
    sorted_idx = bsearch_str(colold_ptr, slen, sorted_marker_ids, max_marker_id_len, marker_ct);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      colold_ptr[slen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --update-name file.\n", colold_ptr);
      goto update_marker_names_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[((uint32_t)sorted_idx)];
    slen = strlen_se(colnew_ptr);
    colnew_ptr[slen] = '\0';
    memcpy(&(true_marker_ids[marker_uidx * max_marker_id_len]), colnew_ptr, slen + 1);
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_names_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-name: %" PRIuPTR " value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-name: %" PRIuPTR " value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  while (0) {
  update_marker_names_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_names_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_names_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 update_marker_names_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_marker_alleles(char* update_alleles_fname, char* sorted_marker_ids, uintptr_t sorted_ids_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, char* outname, char* outname_end) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  FILE* errfile = NULL;
  int32_t retval = 0;
  uintptr_t sorted_ids_ctl = (sorted_ids_ct + (BITCT - 1)) / BITCT;
  uintptr_t max_marker_allele_len = *max_marker_allele_len_ptr;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t err_ct = 0;
  uintptr_t line_idx = 0;
  uintptr_t* already_seen;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  uintptr_t loadbuf_size;
  uint32_t len2;
  uint32_t len;
  int32_t sorted_idx;
  uintptr_t marker_uidx;
  if (wkspace_alloc_ul_checked(&already_seen, sorted_ids_ctl * sizeof(intptr_t))) {
    goto update_marker_alleles_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, sorted_ids_ctl);
  if (fopen_checked(&infile, update_alleles_fname, "r")) {
    goto update_marker_alleles_ret_OPEN_FAIL;
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_alleles_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_alloc(loadbuf_size);
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --update-alleles file is pathologically long.\n", line_idx);
	goto update_marker_alleles_ret_INVALID_FORMAT_2;
      } else {
	goto update_marker_alleles_ret_NOMEM;
      }
    }
    bufptr3 = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr3)) {
      continue;
    }
    bufptr2 = token_endnn(bufptr3);
    sorted_idx = bsearch_str(bufptr3, (uintptr_t)(bufptr2 - bufptr3), sorted_marker_ids, max_marker_id_len, sorted_ids_ct);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      *bufptr2 = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --update-alleles file.\n", bufptr3);
      goto update_marker_alleles_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[((uint32_t)sorted_idx)];
    bufptr2 = skip_initial_spaces(bufptr2);
    len2 = strlen_se(bufptr2);
    bufptr = &(bufptr2[len2]);
    *bufptr = '\0';
    bufptr = skip_initial_spaces(&(bufptr[1]));
    len = strlen_se(bufptr);
    bufptr[len] = '\0';
    if ((!strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx])) && (!strcmp(bufptr, marker_allele_ptrs[2 * marker_uidx + 1]))) {
      bufptr2 = skip_initial_spaces(&(bufptr[len + 1]));
      bufptr = next_token(bufptr2);
      goto update_marker_alleles_match;
    } else if ((!strcmp(bufptr, marker_allele_ptrs[2 * marker_uidx])) && (!strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx + 1]))) {
      bufptr = skip_initial_spaces(&(bufptr[len + 1]));
      bufptr2 = next_token(bufptr);
    update_marker_alleles_match:
      len = strlen_se(bufptr);
      len2 = strlen_se(bufptr2);
      if (len >= max_marker_allele_len) {
	max_marker_allele_len = len + 1;
      }
      if (len2 >= max_marker_allele_len) {
	max_marker_allele_len = len2 + 1;
      }
      allele_reset(&(marker_allele_ptrs[2 * marker_uidx]), bufptr2, len2);
      allele_reset(&(marker_allele_ptrs[2 * marker_uidx + 1]), bufptr, len);
      hit_ct++;
    } else {
      if (!err_ct) {
	memcpy(outname_end, ".allele.no.snp", 15);
	if (fopen_checked(&errfile, outname, "w")) {
	  goto update_marker_alleles_ret_OPEN_FAIL;
	}
      }
      fputs(bufptr3, errfile);
      putc('\t', errfile);
      fputs(bufptr2, errfile);
      putc('\t', errfile);
      fputs(bufptr, errfile);
      if (putc_checked('\n', errfile)) {
	goto update_marker_alleles_ret_WRITE_FAIL;
      }
      err_ct++;
    }
  }
  if (!feof(infile)) {
    goto update_marker_alleles_ret_READ_FAIL;
  }
  *max_marker_allele_len_ptr = max_marker_allele_len;
  if (miss_ct) {
    sprintf(logbuf, "--update-alleles: %" PRIuPTR " variant%s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-alleles: %" PRIuPTR " variant%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  if (err_ct) {
    LOGPRINTFWW("%" PRIuPTR " update failure%s logged to %s .\n", err_ct, (err_ct == 1)? "" : "s", outname);
  }

  while (0) {
  update_marker_alleles_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_alleles_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_marker_alleles_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_alleles_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  update_marker_alleles_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  fclose_cond(errfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

uint32_t flip_str(char** allele_str_ptr) {
  char* allele_str = *allele_str_ptr;
  if (allele_str[1] != '\0') {
    return 1;
  }
  if (allele_str == &(g_one_char_strs[130])) { // "A"
    *allele_str_ptr = (char*)(&(g_one_char_strs[168])); // "T"
  } else if (allele_str == &(g_one_char_strs[134])) {
    *allele_str_ptr = (char*)(&(g_one_char_strs[142]));
  } else if (allele_str == &(g_one_char_strs[142])) {
    *allele_str_ptr = (char*)(&(g_one_char_strs[134]));
  } else if (allele_str == &(g_one_char_strs[168])) {
    *allele_str_ptr = (char*)(&(g_one_char_strs[130]));
  } else if (allele_str == g_missing_geno_ptr) {
    return 0;
  } else {
    return 1;
  }
  return 2;
}

int32_t flip_strand(char* flip_fname, char* sorted_marker_ids, uintptr_t sorted_ids_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char** marker_allele_ptrs) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* flipfile = NULL;
  uint32_t non_acgt_ct = 0;
  int32_t retval = 0;
  uintptr_t sorted_ids_ctl = (sorted_ids_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t line_idx = 0;
  uintptr_t* already_seen;
  char* bufptr;
  uint32_t slen;
  uint32_t cur_non_acgt0;
  int32_t sorted_idx;
  uintptr_t marker_uidx;
  if (wkspace_alloc_ul_checked(&already_seen, sorted_ids_ctl * sizeof(intptr_t))) {
    goto flip_strand_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, sorted_ids_ctl);
  if (fopen_checked(&flipfile, flip_fname, "r")) {
    goto flip_strand_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, flipfile)) {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Line %" PRIuPTR " of --flip file is pathologically long.\n", line_idx);
      goto flip_strand_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    slen = strlen_se(bufptr);
    sorted_idx = bsearch_str(bufptr, slen, sorted_marker_ids, max_marker_id_len, sorted_ids_ct);
    if (sorted_idx == -1) {
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      bufptr[slen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --flip file.\n", bufptr);
      goto flip_strand_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[(uint32_t)sorted_idx];
    cur_non_acgt0 = 0;
    cur_non_acgt0 |= flip_str(&(marker_allele_ptrs[2 * marker_uidx]));
    cur_non_acgt0 |= flip_str(&(marker_allele_ptrs[2 * marker_uidx + 1]));
    non_acgt_ct += (cur_non_acgt0 & 1);
    hit_ct += (cur_non_acgt0 >> 1);
  }
  if (!feof(flipfile)) {
    goto flip_strand_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--flip: %" PRIuPTR " SNP%s flipped, %" PRIuPTR " SNP ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--flip: %" PRIuPTR " SNP%s flipped.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  if (non_acgt_ct) {
    LOGPRINTF("Warning: %u variant%s had at least one non-A/C/G/T allele name.\n", non_acgt_ct, (non_acgt_ct == 1)? "" : "s");
  }
  while (0) {
  flip_strand_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  flip_strand_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  flip_strand_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  flip_strand_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(flipfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_indiv_ids(char* update_ids_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, char* person_ids) {
  // file has been pre-scanned
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  int32_t retval = 0;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t line_idx = 0;
  char* idbuf;
  uintptr_t* already_seen;
  char* bufptr;
  char* bufptr2;
  char* wptr;
  uintptr_t indiv_uidx;
  uint32_t len;
  int32_t sorted_idx;
  if (wkspace_alloc_c_checked(&idbuf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&already_seen, indiv_ctl * sizeof(intptr_t))) {
    goto update_indiv_ids_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, indiv_ctl);
  if (fopen_checked(&infile, update_ids_fname, "r")) {
    goto update_indiv_ids_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      // er, either this buffer should be extended, or the
      // scan_max_fam_indiv_strlen() should use this length...
      sprintf(logbuf, "Error: Line %" PRIuPTR " of --update-ids file is pathologically long.\n", line_idx);
      goto update_indiv_ids_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bsearch_read_fam_indiv(idbuf, sorted_person_ids, max_person_id_len, indiv_ct, bufptr, &bufptr, &sorted_idx);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      *strchr(idbuf, '\t') = ' ';
      LOGPREPRINTFWW("Error: Duplicate individual ID '%s' in --update-ids file.\n", idbuf);
      goto update_indiv_ids_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    indiv_uidx = indiv_id_map[((uint32_t)sorted_idx)];
    wptr = &(person_ids[indiv_uidx * max_person_id_len]);
    len = strlen_se(bufptr);
    bufptr2 = &(bufptr[len]);
    wptr = memcpyax(wptr, bufptr, len, '\t');
    bufptr = skip_initial_spaces(&(bufptr2[1]));
    len = strlen_se(bufptr);
    if ((len == 1) && (*bufptr == '0')) {
      sprintf(logbuf, "Error: Invalid IID '0' on line %" PRIuPTR " of --update-ids file.\n", line_idx);
      goto update_indiv_ids_ret_INVALID_FORMAT_2;
    }
    memcpyx(wptr, bufptr, len, '\0');
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_indiv_ids_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-ids: %" PRIuPTR " %s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, species_str(hit_ct), miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-ids: %" PRIuPTR " %s updated.\n", hit_ct, species_str(hit_ct));
  }
  logprintb();

  while (0) {
  update_indiv_ids_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_indiv_ids_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_indiv_ids_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_indiv_ids_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_indiv_parents(char* update_parents_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  int32_t retval = 0;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  char* idbuf;
  uintptr_t* already_seen;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* wptr;
  uintptr_t loadbuf_size;
  uintptr_t indiv_uidx;
  uint32_t len;
  uint32_t len2;
  int32_t sorted_idx;
  if (wkspace_alloc_c_checked(&idbuf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&already_seen, indiv_ctl * sizeof(intptr_t))) {
    goto update_indiv_parents_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, indiv_ctl);
  if (fopen_checked(&infile, update_parents_fname, "r")) {
    goto update_indiv_parents_ret_OPEN_FAIL;
  }
  // permit very long lines since this can be pointed at .ped files
  if (wkspace_left > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (wkspace_left > MAXLINELEN) {
    loadbuf_size = wkspace_left;
  } else {
    goto update_indiv_parents_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, infile)) {
    // no line_idx since all the validation happened earlier
    if (!loadbuf[loadbuf_size - 1]) {
      goto update_indiv_parents_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bsearch_read_fam_indiv(idbuf, sorted_person_ids, max_person_id_len, indiv_ct, bufptr, &bufptr, &sorted_idx);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      *strchr(idbuf, '\t') = ' ';
      LOGPREPRINTFWW("Error: Duplicate individual ID '%s' in --update-parents file.\n", idbuf);
      goto update_indiv_parents_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    indiv_uidx = indiv_id_map[((uint32_t)sorted_idx)];
    wptr = &(paternal_ids[indiv_uidx * max_paternal_id_len]);
    len = strlen_se(bufptr);
    bufptr2 = &(bufptr[len]);
    memcpyx(wptr, bufptr, len, '\0');
    wptr = &(maternal_ids[indiv_uidx * max_maternal_id_len]);
    bufptr3 = skip_initial_spaces(&(bufptr2[1]));
    len2 = strlen_se(bufptr3);
    memcpyx(wptr, bufptr3, len2, '\0');
    if ((len == 1) && (*bufptr == '0') && (len2 == 1) && (*bufptr3 == '0')) {
      SET_BIT(founder_info, indiv_uidx);
    } else {
      CLEAR_BIT(founder_info, indiv_uidx);
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_indiv_parents_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-parents: %" PRIuPTR " %s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, species_str(hit_ct), miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-parents: %" PRIuPTR " %s updated.\n", hit_ct, species_str(hit_ct));
  }
  logprintb();

  while (0) {
  update_indiv_parents_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_indiv_parents_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_indiv_parents_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_indiv_parents_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_indiv_sexes(char* update_sex_fname, uint32_t update_sex_col, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  int32_t retval = 0;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t line_idx = 0;
  char* idbuf;
  uintptr_t* already_seen;
  char* loadbuf;
  char* bufptr;
  uintptr_t loadbuf_size;
  int32_t sorted_idx;
  uint32_t indiv_uidx;
  char cc;
  unsigned char ucc;
  update_sex_col--;
  if (wkspace_alloc_c_checked(&idbuf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&already_seen, indiv_ctl * sizeof(intptr_t))) {
    goto update_indiv_sexes_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, indiv_ctl);
  if (fopen_checked(&infile, update_sex_fname, "r")) {
    goto update_indiv_sexes_ret_OPEN_FAIL;
  }
  // permit very long lines since this can be pointed at .ped files
  if (wkspace_left > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (wkspace_left > MAXLINELEN) {
    loadbuf_size = wkspace_left;
  } else {
    goto update_indiv_sexes_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --update-sex file is pathologically long.\n", line_idx);
	goto update_indiv_sexes_ret_INVALID_FORMAT_2;
      } else {
	goto update_indiv_sexes_ret_NOMEM;
      }
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(idbuf, sorted_person_ids, max_person_id_len, indiv_ct, bufptr, &bufptr, &sorted_idx)) {
      goto update_indiv_sexes_ret_MISSING_TOKENS;
    }
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      *strchr(idbuf, '\t') = ' ';
      LOGPREPRINTFWW("Error: Duplicate individual ID '%s' in --update-sex file.\n", idbuf);
      goto update_indiv_sexes_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    indiv_uidx = indiv_id_map[((uint32_t)sorted_idx)];
    bufptr = next_token_multz(bufptr, update_sex_col);
    if (no_more_tokens_kns(bufptr)) {
      goto update_indiv_sexes_ret_MISSING_TOKENS;
    }
    cc = *bufptr;
    ucc = ((unsigned char)cc) & 0xdfU;
    if ((cc < '0') || ((cc > '2') && (ucc != 'M') && (ucc != 'F')) || (bufptr[1] > ' ')) {
      sprintf(logbuf, "Error: Invalid sex value on line %" PRIuPTR " of --update-sex file.\n(Acceptable values: 1/M = male, 2/F = female, 0 = missing.)\n", line_idx);
      goto update_indiv_sexes_ret_INVALID_FORMAT_2;
    }
    if (cc == '0') {
      CLEAR_BIT(sex_nm, indiv_uidx);
      CLEAR_BIT(sex_male, indiv_uidx);
    } else {
      SET_BIT(sex_nm, indiv_uidx);
      if ((cc == '1') || (ucc == 'M')) {
	SET_BIT(sex_male, indiv_uidx);
      } else {
	CLEAR_BIT(sex_male, indiv_uidx);
      }
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_indiv_sexes_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-sex: %" PRIuPTR " %s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, species_str(hit_ct), miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-sex: %" PRIuPTR " %s updated.\n", hit_ct, species_str(hit_ct));
  }
  logprintb();

  while (0) {
  update_indiv_sexes_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_indiv_sexes_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_indiv_sexes_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_indiv_sexes_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --update-sex file has fewer tokens than expected.\n", line_idx);
  update_indiv_sexes_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

uint32_t calc_plink_maxsnp(uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len) {
  uintptr_t plink_maxsnp = 4;
  uintptr_t max_marker_id_len_m1 = max_marker_id_len - 1;
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  char* cptr;
  char* cptr_end;
  uintptr_t slen;
  uint32_t marker_uidx_stop;
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
    cptr_end = &(marker_ids[marker_uidx_stop * max_marker_id_len]);
    marker_uidx = marker_uidx_stop;
    do {
      slen = strlen(cptr);
      if (slen > plink_maxsnp) {
	plink_maxsnp = slen + 2;
	if (plink_maxsnp >= max_marker_id_len_m1) {
	  return plink_maxsnp;
	}
      }
      cptr = &(cptr[max_marker_id_len]);
    } while (cptr < cptr_end);
  }
  return plink_maxsnp;
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

// aptr1 = minor, aptr2 = major
int32_t load_one_freq(uint32_t alen1, const char* aptr1, uint32_t alen2, const char* aptr2, double maf, double* set_allele_freq_ptr, char* mastr1, char* mastr2, char missing_geno) {
  uint32_t malen1 = strlen(mastr1);
  uint32_t malen2 = strlen(mastr2);
  uint32_t uii;
  const char* aptr;
  if (maf > 0.5) {
    aptr = aptr2;
    uii = alen2;
    aptr2 = aptr1;
    alen2 = alen1;
    aptr1 = aptr;
    alen1 = uii;
    maf = 1.0 - maf;
  }
  if ((malen1 == alen1) && (!memcmp(mastr1, aptr1, alen1))) {
    if ((malen2 == alen2) && (!memcmp(mastr2, aptr2, alen2))) {
      *set_allele_freq_ptr = 1.0 - maf;
    } else {
      return -1;
    }
  } else if ((malen2 == alen1) && (!memcmp(mastr2, aptr1, alen1))) {
    if ((malen1 == alen2) && (!memcmp(mastr1, aptr2, alen2))) {
      *set_allele_freq_ptr = maf;
    } else {
      return -1;
    }
  } else if ((*aptr1 == missing_geno) && (alen1 == 1) && (maf == 0.0)) {
    if ((malen1 == alen2) && (!memcmp(mastr1, aptr2, alen2))) {
      *set_allele_freq_ptr = 0.0;
    } else if ((malen2 == alen2) && (!memcmp(mastr2, aptr2, alen2))) {
      *set_allele_freq_ptr = 1.0;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
  return 0;
}

int32_t read_external_freqs(char* freqname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, char** marker_allele_ptrs, double* set_allele_freqs, uint32_t maf_succ, double exponent, uint32_t wt_needed, double* marker_weights) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* freqfile = NULL;
  uintptr_t line_idx = 1;
  uint32_t freq_counts = 0;
  uint32_t alen1 = 0;
  uint32_t alen2 = 0;
  char* aptr1 = NULL;
  char* aptr2 = NULL;
  int32_t retval = 0;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  char missing_geno = *missing_geno_ptr;
  char* loadbuf;
  char* sorted_ids;
  uint32_t* id_map;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* bufptr5;
  double maf;
  uintptr_t loadbuf_size;
  uint32_t chrom_idx;
  uint32_t marker_uidx;
  uint32_t uii;
  int32_t c_hom_a1;
  int32_t c_het;
  int32_t c_hom_a2;
  int32_t c_hap_a1;
  int32_t c_hap_a2;
  int32_t ii;
  if (fopen_checked(&freqfile, freqname, "r")) {
    goto read_external_freqs_ret_OPEN_FAIL;
  }
  retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    goto read_external_freqs_ret_1;
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto read_external_freqs_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  if (!fgets(loadbuf, loadbuf_size, freqfile)) {
    logprint("Error: Empty --read-freq file.\n");
    goto read_external_freqs_ret_INVALID_FORMAT;
  }
  if (!memcmp(loadbuf, " CHR  ", 6)) {
    uii = strlen(loadbuf);
    if (loadbuf[uii - 2] == '0') { // --counts makes G0 the last column header
      freq_counts = 1;
    } else if (loadbuf[uii - 2] != 'S') { // NCHROBS
      logprint("Error: Invalid --read-freq file header.\n");
      goto read_external_freqs_ret_INVALID_FORMAT;
    }
    while (fgets(loadbuf, loadbuf_size, freqfile) != NULL) {
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      bufptr = skip_initial_spaces(loadbuf);
      ii = get_chrom_code(chrom_info_ptr, bufptr);
      if (ii < 0) {
	goto read_external_freqs_ret_INVALID_CHROM;
      }
      chrom_idx = ii;
      bufptr = next_token(bufptr); // now at beginning of marker name
      bufptr2 = next_token(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_MISSING_TOKENS;
      }
      ii = bsearch_str(bufptr, strlen_se(bufptr), sorted_ids, max_marker_id_len, unfiltered_marker_ct - marker_exclude_ct);
      if (ii != -1) {
        marker_uidx = id_map[(uint32_t)ii];
        if ((chrom_idx == get_marker_chrom(chrom_info_ptr, marker_uidx)) || (!chrom_idx) || (!get_marker_chrom(chrom_info_ptr, marker_uidx))) {
	  alen1 = strlen_se(bufptr2);
	  aptr1 = bufptr2;
	  bufptr2 = next_token(bufptr2);
	  if (no_more_tokens_kns(bufptr2)) {
	    goto read_external_freqs_ret_MISSING_TOKENS;
	  }
	  alen2 = strlen_se(bufptr2);
	  aptr2 = bufptr2;
	  if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	    goto read_external_freqs_ret_A1_A2_SAME;
	  }
	  bufptr = next_token(bufptr2);
	  if (no_more_tokens_kns(bufptr)) {
	    goto read_external_freqs_ret_MISSING_TOKENS;
	  }
	  if (freq_counts) {
	    if (no_more_tokens_kns(next_token(bufptr))) {
	      goto read_external_freqs_ret_MISSING_TOKENS;
	    }
	    if (scan_uint_icap(bufptr, (uint32_t*)&c_hom_a1)) {
	      goto read_external_freqs_ret_INVALID_HOM_A1;
	    }
	    if (scan_uint_icap(next_token(bufptr), (uint32_t*)&c_hom_a2)) {
	      goto read_external_freqs_ret_INVALID_HOM_A2;
	    }
	    maf = ((double)c_hom_a1 + maf_succ) / ((double)(c_hom_a1 + c_hom_a2 + 2 * maf_succ));
	  } else {
	    if (scan_double(bufptr, &maf)) {
	      goto read_external_freqs_ret_INVALID_MAF;
	    }
	  }
	  if (load_one_freq(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[marker_uidx]), marker_allele_ptrs[marker_uidx * 2], marker_allele_ptrs[marker_uidx * 2 + 1], missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
	  }
	  if (wt_needed) {
	    marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, set_allele_freqs[marker_uidx]);
	  }
        }
      }
    }
    if (freq_counts) {
      logprint("--read-freq: .frq.count file loaded.\n");
    } else {
      logprint("--read-freq: .frq file loaded.\n");
    }
  } else if (!memcmp(loadbuf, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)", 71)) {
    // changed from strcmp to avoid eoln problems
    // known --freqx format, WDIST v0.15.3 or later
    while (fgets(loadbuf, loadbuf_size, freqfile) != NULL) {
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      ii = get_chrom_code(chrom_info_ptr, loadbuf);
      if (ii < 0) {
	goto read_external_freqs_ret_INVALID_CHROM;
      }
      chrom_idx = ii;
      bufptr = next_token(loadbuf); // now at beginning of marker name
      bufptr2 = next_token(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_MISSING_TOKENS;
      }
      ii = bsearch_str(bufptr, strlen_se(bufptr), sorted_ids, max_marker_id_len, unfiltered_marker_ct - marker_exclude_ct);
      if (ii != -1) {
        marker_uidx = id_map[(uint32_t)ii];
        if ((chrom_idx == get_marker_chrom(chrom_info_ptr, marker_uidx)) || (!chrom_idx) || (!get_marker_chrom(chrom_info_ptr, marker_uidx))) {
	  alen1 = strlen_se(bufptr2);
	  aptr1 = bufptr2;
	  bufptr2 = next_token(bufptr2);
	  if (no_more_tokens_kns(bufptr2)) {
	    goto read_external_freqs_ret_MISSING_TOKENS;
	  }
	  alen2 = strlen_se(bufptr2);
	  aptr2 = bufptr2;
	  if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	    goto read_external_freqs_ret_A1_A2_SAME;
	  }
	  bufptr = next_token(bufptr2);
	  bufptr2 = next_token(bufptr);
	  bufptr3 = next_token(bufptr2);
	  bufptr4 = next_token(bufptr3);
	  bufptr5 = next_token(bufptr4);
	  if (no_more_tokens_kns(bufptr5)) {
	    goto read_external_freqs_ret_MISSING_TOKENS;
	  }
	  if (scan_uint_icap(bufptr, (uint32_t*)&c_hom_a1)) {
	    goto read_external_freqs_ret_INVALID_HOM_A1;
	  }
	  if (scan_uint_icap(bufptr2, (uint32_t*)&c_het)) {
	    sprintf(logbuf, "Error: Invalid het count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	    goto read_external_freqs_ret_INVALID_FORMAT_2;
	  }
	  if (scan_uint_icap(bufptr3, (uint32_t*)&c_hom_a2)) {
	    goto read_external_freqs_ret_INVALID_HOM_A2;
	  }
	  if (scan_uint_icap(bufptr4, (uint32_t*)&c_hap_a1)) {
	    sprintf(logbuf, "Error: Invalid hap. A1 count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	    goto read_external_freqs_ret_INVALID_FORMAT_2;
	  }
	  if (scan_uint_icap(bufptr5, (uint32_t*)&c_hap_a2)) {
	    sprintf(logbuf, "Error: Invalid hap. A2 count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	    goto read_external_freqs_ret_INVALID_FORMAT_2;
	  }
	  maf = ((double)(c_hom_a1 * 2 + c_het + c_hap_a1 + maf_succ)) / ((double)(2 * (c_hom_a1 + c_het + c_hom_a2 + maf_succ) + c_hap_a1 + c_hap_a2));
	  if (load_one_freq(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[marker_uidx]), marker_allele_ptrs[marker_uidx * 2], marker_allele_ptrs[marker_uidx * 2 + 1], missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
	  }
	  if (wt_needed) {
	    if (c_hap_a1 || c_hap_a2) {
	      marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, set_allele_freqs[marker_uidx]);
	    } else {
	      marker_weights[marker_uidx] = calc_wt_mean(exponent, c_het, c_hom_a1, c_hom_a2);
	    }
	  }
        }
      }
    }
    logprint("--read-freq: .frqx file loaded.\n");
  } else {
    // Also support GCTA-style frequency files:
    // [marker ID]\t[reference allele]\t[frequency of reference allele]\n
    line_idx = 0; // no header line here
    do {
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      bufptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      bufptr = next_token(bufptr);
      if (!bufptr) {
        goto read_external_freqs_ret_MISSING_TOKENS;
      }
      ii = bsearch_str(loadbuf, strlen_se(loadbuf), sorted_ids, max_marker_id_len, unfiltered_marker_ct - marker_exclude_ct);
      if (ii != -1) {
        marker_uidx = id_map[(uint32_t)ii];
	alen1 = strlen_se(bufptr);
	aptr1 = bufptr;
        bufptr = next_token(bufptr);
	if (no_more_tokens_kns(bufptr)) {
          goto read_external_freqs_ret_MISSING_TOKENS;
	}
	if (scan_double(bufptr, &maf)) {
          goto read_external_freqs_ret_INVALID_MAF;
        }
	if (load_one_freq(1, missing_geno_ptr, alen1, aptr1, maf, &(set_allele_freqs[marker_uidx]), marker_allele_ptrs[marker_uidx * 2], marker_allele_ptrs[marker_uidx * 2 + 1], missing_geno)) {
	  goto read_external_freqs_ret_ALLELE_MISMATCH;
	}
	if (wt_needed) {
	  marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, set_allele_freqs[marker_uidx]);
	}
      } else {
	// if there aren't exactly 3 columns, this isn't a GCTA .freq file
	bufptr = next_token(bufptr);
	if (no_more_tokens_kns(bufptr)) {
	  goto read_external_freqs_ret_MISSING_TOKENS;
	}
        if (!no_more_tokens_kns(next_token(bufptr))) {
	  sprintf(logbuf, "Error: Line %" PRIuPTR " of --read-freq has more tokens than expected.\n", line_idx);
	  goto read_external_freqs_ret_INVALID_FORMAT_2;
	}
      }
    } while (fgets(loadbuf, loadbuf_size, freqfile));
    logprint("--read-freq: GCTA-formatted .freq file loaded.\n");
  }
  while (0) {
  read_external_freqs_ret_TOO_LONG_LINE:
    if (loadbuf_size == MAXLINEBUFLEN) {
      LOGPRINTF("Error: Line %" PRIuPTR " of --read-freq file is pathologically long.\n", line_idx);
      retval = RET_INVALID_FORMAT;
      break;
    }
  read_external_freqs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  read_external_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  read_external_freqs_ret_INVALID_HOM_A1:
    LOGPRINTF("Error: Invalid hom. A1 count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_INVALID_HOM_A2:
    LOGPRINTF("Error: Invalid hom. A2 count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_INVALID_MAF:
    LOGPRINTF("Error: Invalid MAF on line %" PRIuPTR " of --read-freq file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_A1_A2_SAME:
    LOGPRINTF("Error: A1 and A2 alleles match on line %" PRIuPTR " of --read-freq file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_INVALID_CHROM:
    sprintf(logbuf, "Error: Invalid chromosome code on line %" PRIuPTR" of --read-freq file.\n", line_idx);
  read_external_freqs_ret_INVALID_FORMAT_2:
    logprintb();
  read_external_freqs_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_MISSING_TOKENS:
    LOGPRINTF("Error: Line %" PRIuPTR " of --read-freq file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_ALLELE_MISMATCH:
    LOGPRINTF("Error: Allele(s) on line %" PRIuPTR " of --read-freq file don't match loaded\nvalues.\n", line_idx);
    retval = RET_ALLELE_MISMATCH;
    break;
  }
 read_external_freqs_ret_1:
  fclose_cond(freqfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t load_ax_alleles(Two_col_params* axalleles, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, double* set_allele_freqs, uint32_t is_a2) {
  // note that swap_reversed_marker_alleles() has NOT been called yet
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = axalleles->skipchar;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  uint32_t colid_first = (axalleles->colid < axalleles->colx);
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t max_marker_allele_len = *max_marker_allele_len_ptr;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uintptr_t line_idx;
  uint32_t idlen;
  uint32_t alen;
  char* sorted_marker_ids;
  char* colid_ptr;
  char* colx_ptr;
  uint32_t* marker_id_map;
  uint32_t colmin;
  uint32_t coldiff;
  int32_t sorted_idx;
  uint32_t marker_uidx;
  char cc;
  uint32_t replace_other;
  int32_t retval;
  retval = sort_item_ids(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    goto load_ax_alleles_ret_1;
  }
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto load_ax_alleles_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_ax_alleles_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  retval = open_and_skip_first_lines(&infile, axalleles->fname, loadbuf, loadbuf_size, axalleles->skip);
  if (colid_first) {
    colmin = axalleles->colid - 1;
    coldiff = axalleles->colx - axalleles->colid;
  } else {
    colmin = axalleles->colx - 1;
    coldiff = axalleles->colid - axalleles->colx;
  }
  line_idx = axalleles->skip;
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, axalleles->fname);
	goto load_ax_alleles_ret_INVALID_FORMAT_2;
      } else {
        goto load_ax_alleles_ret_NOMEM;
      }
    }
    colid_ptr = skip_initial_spaces(loadbuf);
    cc = *colid_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colid_first) {
      colid_ptr = next_token_multz(colid_ptr, colmin);
      colx_ptr = next_token_mult(colid_ptr, coldiff);
    } else {
      colx_ptr = next_token_multz(colid_ptr, colmin);
      colid_ptr = next_token_mult(colx_ptr, coldiff);
    }
    idlen = strlen_se(colid_ptr);
    sorted_idx = bsearch_str(colid_ptr, idlen, sorted_marker_ids, max_marker_id_len, marker_ct);
    if (sorted_idx == -1) {
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      colid_ptr[idlen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --a%c-allele file.\n", colid_ptr, is_a2? '2' : '1');
      goto load_ax_alleles_ret_INVALID_FORMAT_2;
    }
    SET_BIT(already_seen, sorted_idx);
    marker_uidx = marker_id_map[(uint32_t)sorted_idx];
    alen = strlen_se(colx_ptr);
    colx_ptr[alen] = '\0';
    if (!strcmp(colx_ptr, marker_allele_ptrs[marker_uidx * 2 + is_a2])) {
      if (IS_SET(marker_reverse, marker_uidx)) {
        set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
        CLEAR_BIT(marker_reverse, marker_uidx);
      }
    } else if (!strcmp(colx_ptr, marker_allele_ptrs[marker_uidx * 2 + 1 - is_a2])) {
      if (!IS_SET(marker_reverse, marker_uidx)) {
        set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
        SET_BIT(marker_reverse, marker_uidx);
      }
    } else if ((marker_allele_ptrs[marker_uidx * 2] == missing_geno_ptr) || (marker_allele_ptrs[marker_uidx * 2 + 1] == missing_geno_ptr)) {
      replace_other = (marker_allele_ptrs[marker_uidx * 2 + is_a2] != missing_geno_ptr);
      if (allele_reset(&(marker_allele_ptrs[marker_uidx * 2 + (is_a2 ^ replace_other)]), colx_ptr, alen)) {
	goto load_ax_alleles_ret_NOMEM;
      }
      if (alen >= max_marker_allele_len) {
	max_marker_allele_len = alen + 1;
      }
      if (!replace_other) {
	if (IS_SET(marker_reverse, marker_uidx)) {
	  set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
	  CLEAR_BIT(marker_reverse, marker_uidx);
	}
      } else {
	if (!IS_SET(marker_reverse, marker_uidx)) {
	  set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
	  SET_BIT(marker_reverse, marker_uidx);
	}
      }
    } else {
      colid_ptr[idlen] = '\0';
      LOGPRINTF("Warning: Impossible A%c allele assignment for variant %s.\n", is_a2? '2' : '1', colid_ptr);
    }
  }
  if (!feof(infile)) {
    goto load_ax_alleles_ret_READ_FAIL;
  }
  marker_uidx = popcount_longs(already_seen, marker_ctl);
  LOGPRINTF("--a%c-allele: %u assignment%s made.\n", is_a2? '2' : '1', marker_uidx, (marker_uidx == 1)? "" : "s");
  *max_marker_allele_len_ptr = max_marker_allele_len;
  while (0) {
  load_ax_alleles_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_ax_alleles_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_ax_alleles_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 load_ax_alleles_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t write_stratified_freqs(FILE* bedfile, uintptr_t bed_offset, char* outname, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uint32_t indiv_f_ct, uintptr_t* founder_info, uint32_t nonfounders, uintptr_t* sex_male, uint32_t indiv_f_male_ct, uintptr_t* marker_reverse, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uint32_t* cur_cluster_map = cluster_map;
  uint32_t* cur_cluster_starts = cluster_starts;
  uint32_t* cluster_map_nonmale = NULL;
  uint32_t* cluster_starts_nonmale = NULL;
  uint32_t* cluster_map_male = NULL;
  uint32_t* cluster_starts_male = NULL;
  int32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  uint32_t cslen = 10;
  int32_t retval = 0;
  uint32_t cur_cts[4];
  uintptr_t* readbuf;
  uint32_t* uiptr;
  uint32_t* uiptr2;
  uint32_t* uiptr3;
  char* csptr;
  char* col_2_start;
  char* wptr_start;
  char* wptr;
  char* sptr;
  uintptr_t chrom_end;
  uintptr_t marker_uidx;
  int32_t chrom_idx;
  uintptr_t clidx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t clmpos;
  uint32_t a1_obs;
  uint32_t tot_obs;
  uint32_t uii;
  if (wkspace_alloc_ul_checked(&readbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto write_stratified_freqs_ret_NOMEM;
  }
  if ((indiv_ct > indiv_f_ct) && (!nonfounders)) {
    if (wkspace_alloc_ui_checked(&cur_cluster_starts, (cluster_ct + 1) * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cur_cluster_map, indiv_f_ct * sizeof(int32_t))) {
      goto write_stratified_freqs_ret_NOMEM;
    }
    clmpos = 0;
    cur_cluster_starts[0] = 0;
    uiptr = cluster_map;
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      uiptr2 = &(cluster_map[cluster_starts[clidx + 1]]);
      do {
	uii = *uiptr;
	if (IS_SET(founder_info, uii)) {
          cur_cluster_map[clmpos++] = uii;
	}
      } while ((++uiptr) < uiptr2);
      cur_cluster_starts[clidx + 1] = clmpos;
    }
  }
  chrom_idx = chrom_info_ptr->x_code;
  if ((chrom_idx != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
    if (wkspace_alloc_ui_checked(&cluster_starts_nonmale, (cluster_ct + 1) * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_map_nonmale, (indiv_f_ct - indiv_f_male_ct) * sizeof(int32_t))) {
      goto write_stratified_freqs_ret_NOMEM;
    }
    clmpos = 0;
    cluster_starts_nonmale[0] = 0;
    uiptr = cur_cluster_map;
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      uiptr2 = &(cur_cluster_map[cur_cluster_starts[clidx + 1]]);
      while (uiptr < uiptr2) {
	uii = *uiptr++;
	if (!IS_SET(sex_male, uii)) {
          cluster_map_nonmale[clmpos++] = uii;
	}
      }
      cluster_starts_nonmale[clidx + 1] = clmpos;
    }
  }
  chrom_idx = chrom_info_ptr->y_code;
  if (cluster_map_nonmale || ((chrom_idx != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_idx))) {
    if (wkspace_alloc_ui_checked(&cluster_starts_male, (cluster_ct + 1) * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_map_male, indiv_f_male_ct * sizeof(int32_t))) {
      goto write_stratified_freqs_ret_NOMEM;
    }
    clmpos = 0;
    cluster_starts_male[0] = 0;
    uiptr = cur_cluster_map;
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      uiptr2 = &(cur_cluster_map[cur_cluster_starts[clidx + 1]]);
      while (uiptr < uiptr2) {
	uii = *uiptr++;
	if (IS_SET(sex_male, uii)) {
          cluster_map_male[clmpos++] = uii;
	}
      }
      cluster_starts_male[clidx + 1] = clmpos;
    }
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_stratified_freqs_ret_OPEN_FAIL;
  }
  sprintf(tbuf, " CHR %%%ds     CLST   A1   A2      MAF    MAC  NCHROBS\n", plink_maxsnp);
  fprintf(outfile, tbuf, "SNP");
  if (wkspace_alloc_c_checked(&csptr, 2 * max_marker_allele_len + 16)) {
    goto write_stratified_freqs_ret_NOMEM;
  }
  memset(csptr, 32, 10);
  for (chrom_idx = 0; chrom_idx < chrom_code_end; chrom_idx++) {
    if (!chrom_exists(chrom_info_ptr, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == chrom_info_ptr->x_code);
    is_y = (chrom_idx == chrom_info_ptr->y_code);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    chrom_end = chrom_info_ptr->chrom_end[chrom_idx];
    marker_uidx = next_unset_ul(marker_exclude, chrom_info_ptr->chrom_start[chrom_idx], chrom_end);
    if (marker_uidx >= chrom_end) {
      continue;
    }
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
      goto write_stratified_freqs_ret_READ_FAIL;
    }
    col_2_start = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, chrom_idx, zero_extra_chroms));
    *col_2_start++ = ' ';
    do {
      sptr = &(marker_ids[marker_uidx * max_marker_id_len]);
      uii = strlen(sptr);
      wptr_start = memseta(col_2_start, 32, plink_maxsnp - uii);
      wptr_start = memcpyax(wptr_start, sptr, uii, ' ');
      sptr = marker_allele_ptrs[marker_uidx * 2];
      wptr = fw_strcpy(4, sptr, &(csptr[1]));
      *wptr++ = ' ';
      sptr = marker_allele_ptrs[marker_uidx * 2 + 1];
      wptr = fw_strcpy(4, sptr, wptr);
      *wptr++ = ' ';
      cslen = (uintptr_t)(wptr - csptr);

      if (fread(readbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto write_stratified_freqs_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf((unsigned char*)readbuf, unfiltered_indiv_ct);
      }
      if (is_x) {
	uiptr = cluster_map_nonmale;
	uiptr2 = cluster_map_male;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  wptr = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), wptr_start);
	  wptr = memcpyax(wptr, csptr, cslen, ' ');
	  fill_uint_zero(cur_cts, 4);
	  uiptr3 = &(cluster_map_nonmale[cluster_starts_nonmale[clidx + 1]]);
	  while (uiptr < uiptr3) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  a1_obs = 2 * cur_cts[0] + cur_cts[2];
	  tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  fill_uint_zero(cur_cts, 4);
	  uiptr3 = &(cluster_map_male[cluster_starts_male[clidx + 1]]);
	  while (uiptr2 < uiptr3) {
	    uii = *uiptr2++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  a1_obs += cur_cts[0];
	  tot_obs += cur_cts[0] + cur_cts[3];
	  if (tot_obs) {
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ');
	    wptr = uint32_writew6x(wptr, a1_obs, ' ');
	    wptr = uint32_writew8(wptr, tot_obs);
	    wptr = memcpya(wptr, " \n", 2);
	  } else {
	    wptr = memcpya(wptr, "       0      0        0 \n", 26);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      } else if (is_y) {
	uiptr = cluster_map_male;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  wptr = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), wptr_start);
	  wptr = memcpyax(wptr, csptr, cslen, ' ');
	  fill_uint_zero(cur_cts, 4);
	  uiptr2 = &(cluster_map_male[cluster_starts_male[clidx + 1]]);
	  while (uiptr < uiptr2) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  if (is_haploid) {
	    a1_obs = cur_cts[0];
	    tot_obs = cur_cts[0] + cur_cts[3];
	  } else {
	    a1_obs = 2 * cur_cts[0] + cur_cts[2];
	    tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  }
	  if (tot_obs) {
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ');
	    wptr = uint32_writew6x(wptr, a1_obs, ' ');
	    wptr = uint32_writew8(wptr, tot_obs);
	    wptr = memcpya(wptr, " \n", 2);
	  } else {
	    wptr = memcpya(wptr, "       0      0        0 \n", 26);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      } else {
        uiptr = cur_cluster_map;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  wptr = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), wptr_start);
	  wptr = memcpyax(wptr, csptr, cslen, ' ');
	  fill_uint_zero(cur_cts, 4);
	  uiptr2 = &(cur_cluster_map[cur_cluster_starts[clidx + 1]]);
	  while (uiptr < uiptr2) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  if (is_haploid) {
	    a1_obs = cur_cts[0];
	    tot_obs = cur_cts[0] + cur_cts[3];
	  } else {
	    a1_obs = 2 * cur_cts[0] + cur_cts[2];
	    tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  }
	  if (tot_obs) {
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ');
	    wptr = uint32_writew6x(wptr, a1_obs, ' ');
	    wptr = uint32_writew8(wptr, tot_obs);
	    wptr = memcpya(wptr, " \n", 2);
	  } else {
	    wptr = memcpya(wptr, "       0      0        0 \n", 26);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto write_stratified_freqs_ret_READ_FAIL;
	}
      }
    } while (marker_uidx < chrom_end);
  }
  if (fclose_null(&outfile)) {
    goto write_stratified_freqs_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--freq: Cluster-stratified allele frequencies written to %s .\n", outname);
  while (0) {
  write_stratified_freqs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_stratified_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_stratified_freqs_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  write_stratified_freqs_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t write_freqs(char* outname, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, double* set_allele_freqs, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, int32_t* hapl_cts, int32_t* haph_cts, uint32_t indiv_f_ct, uint32_t indiv_f_male_ct, uint64_t misc_flags, uintptr_t* marker_reverse) {
  FILE* outfile = NULL;
  uint32_t reverse = 0;
  uint32_t freq_counts = (misc_flags / MISC_FREQ_COUNTS) & 1;
  uint32_t freqx = (misc_flags / MISC_FREQX) & 1;
  uint32_t maf_succ = (misc_flags / MISC_MAF_SUCC) & 1;
  int32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  int32_t retval = 0;
  char* minor_ptr;
  char* major_ptr;
  char* bufptr;
  uint32_t chrom_end;
  uint32_t marker_uidx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t missing_ct;
  int32_t chrom_idx;
  uint32_t uii;
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_freqs_ret_OPEN_FAIL;
  }
  if (freqx) {
    if (fputs_checked("CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)\n", outfile)) {
      goto write_freqs_ret_WRITE_FAIL;
    }
  } else if (plink_maxsnp < 5) {
    if (freq_counts) {
      if (fputs_checked(" CHR  SNP   A1   A2     C1     C2     G0\n", outfile)) {
	goto write_freqs_ret_WRITE_FAIL;
      }
    } else {
      if (fputs_checked(" CHR  SNP   A1   A2          MAF  NCHROBS\n", outfile)) {
        goto write_freqs_ret_WRITE_FAIL;
      }
    }
  } else if (freq_counts) {
    sprintf(tbuf, " CHR %%%us   A1   A2     C1     C2     G0\n", plink_maxsnp);
    fprintf(outfile, tbuf, "SNP");
  } else {
    sprintf(tbuf, " CHR %%%us   A1   A2          MAF  NCHROBS\n", plink_maxsnp);
    fprintf(outfile, tbuf, "SNP");
  }
  if (ferror(outfile)) {
    goto write_freqs_ret_WRITE_FAIL;
  }
  for (chrom_idx = 0; chrom_idx < chrom_code_end; chrom_idx++) {
    if (!chrom_exists(chrom_info_ptr, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == chrom_info_ptr->x_code);
    is_y = (chrom_idx == chrom_info_ptr->y_code);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    chrom_end = chrom_info_ptr->chrom_end[chrom_idx];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_start[chrom_idx], chrom_end);
    while (marker_uidx < chrom_end) {
      reverse = IS_SET(marker_reverse, marker_uidx);
      major_ptr = marker_allele_ptrs[marker_uidx * 2 + 1];
      minor_ptr = marker_allele_ptrs[marker_uidx * 2];
      if (freq_counts || freqx) {
	if (is_x) {
	  missing_ct = indiv_f_ct - (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx] + hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	} else if (is_haploid) {
	  if (is_y) {
	    missing_ct = indiv_f_male_ct - (hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	  } else {
	    missing_ct = indiv_f_ct - (hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	  }
	} else {
	  missing_ct = indiv_f_ct - (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx]);
	}
	if (freqx) {
	  bufptr = chrom_name_write(tbuf, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms);
	  *bufptr++ = '\t';
	  bufptr = strcpyax(bufptr, &(marker_ids[marker_uidx * max_marker_id_len]), '\t');
	  fwrite(tbuf, 1, bufptr - tbuf, outfile);
          fputs(minor_ptr, outfile);
	  putc('\t', outfile);
          fputs(major_ptr, outfile);
	  tbuf[0] = '\t';
          bufptr = uint32_writex(&(tbuf[1]), reverse? hh_cts[marker_uidx] : ll_cts[marker_uidx], '\t');
	  bufptr = uint32_writex(bufptr, lh_cts[marker_uidx], '\t');
          bufptr = uint32_writex(bufptr, reverse? ll_cts[marker_uidx] : hh_cts[marker_uidx], '\t');
          bufptr = uint32_writex(bufptr, reverse? haph_cts[marker_uidx] : hapl_cts[marker_uidx], '\t');
          bufptr = uint32_writex(bufptr, reverse? hapl_cts[marker_uidx] : haph_cts[marker_uidx], '\t');
          bufptr = uint32_writex(bufptr, missing_ct, '\n');
	  fwrite(tbuf, 1, bufptr - tbuf, outfile);
	} else {
	  bufptr = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms));
	  *bufptr++ = ' ';
	  bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
	  *bufptr++ = ' ';
	  fwrite(tbuf, 1, bufptr - tbuf, outfile);
	  fputs_w4(minor_ptr, outfile);
          putc(' ', outfile);
          fputs_w4(major_ptr, outfile);
	  tbuf[0] = ' ';
          bufptr = uint32_writew6x(&(tbuf[1]), 2 * ll_cts[marker_uidx] + lh_cts[marker_uidx] + hapl_cts[marker_uidx], ' ');
	  bufptr = uint32_writew6x(bufptr, 2 * hh_cts[marker_uidx] + lh_cts[marker_uidx] + haph_cts[marker_uidx], ' ');
	  bufptr = uint32_writew6x(bufptr, missing_ct, '\n');
	  fwrite(tbuf, 1, bufptr - tbuf, outfile);
	}
      } else {
	bufptr = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms));
	*bufptr++ = ' ';
	bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
        *bufptr++ = ' ';
	fwrite(tbuf, 1, bufptr - tbuf, outfile);
	fputs_w4(minor_ptr, outfile);
	putc(' ', outfile);
        fputs_w4(major_ptr, outfile);
        tbuf[0] = ' ';
	uii = 2 * (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx]) + hapl_cts[marker_uidx] + haph_cts[marker_uidx];
	if (maf_succ || uii || (set_allele_freqs[marker_uidx] != 0.5)) {
          bufptr = double_g_writewx4(&(tbuf[1]), 1.0 - set_allele_freqs[marker_uidx], 12);
	} else {
	  bufptr = memcpya(&(tbuf[1]), "          NA", 12);
	}
	*bufptr++ = ' ';
        bufptr = uint32_writew8x(bufptr, uii, '\n');
	fwrite(tbuf, 1, bufptr - tbuf, outfile);
      }
      if (ferror(outfile)) {
	goto write_freqs_ret_WRITE_FAIL;
      }
      marker_uidx = next_unset(marker_exclude, marker_uidx + 1, chrom_end);
    }
  }
  if (fclose_null(&outfile)) {
    goto write_freqs_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--freq%s: Allele frequencies written to %s .\n", freqx? "x" : "", outname);
  while (0) {
  write_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_freqs_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

void calc_marker_weights(double exponent, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, double* marker_weights) {
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  do {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    do {
      if (marker_weights[marker_uidx] < 0.0) {
	marker_weights[marker_uidx] = calc_wt_mean(exponent, lh_cts[marker_uidx], ll_cts[marker_uidx], hh_cts[marker_uidx]);
      }
    } while (++marker_uidx < marker_uidx_stop);
  } while (markers_done < marker_ct);
}

int32_t sexcheck(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uint64_t misc_flags, double check_sex_fthresh, double check_sex_mthresh, uint32_t max_f_yobs, uint32_t min_m_yobs, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* gender_unk_ct_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uint32_t* het_cts = NULL;
  uint32_t* missing_cts = NULL;
  double* nei_offsets = NULL;
  uint32_t* ymiss_cts = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t indiv_ctl2 = (indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t x_variant_ct = 0;
  uintptr_t ytotal = 0;
  double nei_sum = 0.0;
  uint32_t do_impute = (misc_flags / MISC_IMPUTE_SEX) & 1;
  uint32_t check_y = (misc_flags / MISC_SEXCHECK_YCOUNT) & 1;
  uint32_t yonly = (misc_flags / MISC_SEXCHECK_YONLY) & 1;
  uint32_t gender_unk_ct = 0;
  uint32_t problem_ct = 0;
  int32_t x_code = chrom_info_ptr->x_code;
  int32_t y_code = chrom_info_ptr->y_code;
  int32_t retval = 0;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  // We wish to compute four quantities for each individual:
  // 1. Observed homozygous polymorphic Xchr sites.
  // 2. Observed nonmissing polymorphic Xchr sites.
  // 3. Number of nonmissing calls on Ychr.
  // 4. Nei's unbiased estimator of the expected quantity #1.  (This has an
  //    N/(N-1) term which makes it slightly different from GCTA's Fhat2.)
  // edit: Actually, forget about the 2N/(2N-1) multiplier for now since it
  // doesn't play well with --freq, and PLINK 1.07 only succeeds in applying it
  // when N=1 due to use of integer division.  Maybe let it be used with
  // inferred MAFs/--freqx/--freq counts later.
  uintptr_t* lptr;
  char* fid_ptr;
  char* iid_ptr;
  char* wptr;
  double dpp;
  double dtot;
  double cur_nei;
  double dff;
  double dee;
  uintptr_t marker_uidx;
  uintptr_t marker_uidx_end;
  uintptr_t marker_idxs_left;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t cur_missing_ct;
  uintptr_t allele_obs_ct;
  uintptr_t cur_word;
  uintptr_t ulii;
  uint32_t orig_sex_code;
  uint32_t imputed_sex_code;
  if (wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, indiv_ctl2 * sizeof(intptr_t))) {
    goto sexcheck_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_indiv_ctl2 - 1] = 0;
  loadbuf[indiv_ctl2 - 1] = 0;
  if (check_y) {
    if (wkspace_alloc_ui_checked(&ymiss_cts, indiv_ct * sizeof(int32_t))) {
      goto sexcheck_ret_NOMEM;
    }
    fill_uint_zero(ymiss_cts, indiv_ct);
  }
  if (!yonly) {
    if ((x_code == -1) || (!is_set(chrom_info_ptr->chrom_mask, (uint32_t)x_code))) {
      goto sexcheck_ret_NO_X_VAR;
    }
    marker_uidx_end = chrom_info_ptr->chrom_end[(uint32_t)x_code];
    marker_uidx = next_unset_ul(marker_exclude, chrom_info_ptr->chrom_start[(uint32_t)x_code], marker_uidx_end);
    if (marker_uidx == marker_uidx_end) {
      goto sexcheck_ret_NO_X_VAR;
    }
    marker_idxs_left = marker_uidx_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, marker_uidx_end);
    if (wkspace_alloc_ui_checked(&het_cts, indiv_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&missing_cts, indiv_ct * sizeof(int32_t)) ||
        wkspace_alloc_d_checked(&nei_offsets, indiv_ct * sizeof(double))) {
      goto sexcheck_ret_NOMEM;
    }
    fill_uint_zero(het_cts, indiv_ct);
    fill_uint_zero(missing_cts, indiv_ct);
    fill_double_zero(nei_offsets, indiv_ct);
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
      goto sexcheck_ret_READ_FAIL;
    }
    for (; marker_idxs_left; marker_idxs_left--, marker_uidx++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto sexcheck_ret_READ_FAIL;
	}
      }
      if (load_and_collapse(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf, indiv_ct, indiv_exclude, 0)) {
	goto sexcheck_ret_READ_FAIL;
      }
      cur_missing_ct = count_01(loadbuf, indiv_ctl2);
      allele_obs_ct = 2 * (indiv_ct - cur_missing_ct);
      dpp = set_allele_freqs[marker_uidx];
      // skip monomorphic sites
      if ((!allele_obs_ct) || (dpp < 1e-8) || (dpp > (1 - 1e-8))) {
	continue;
      }
      cur_nei = 1.0 - 2 * dpp * (1 - dpp);
      x_variant_ct++;
      if (cur_missing_ct) {
	// iterate through missing calls
	lptr = loadbuf;
	for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx += BITCT2) {
	  cur_word = *lptr++;
	  cur_word = (~(cur_word >> 1)) & cur_word & FIVEMASK;
	  while (cur_word) {
	    ulii = indiv_idx + CTZLU(cur_word) / 2;
	    missing_cts[ulii] += 1;
	    nei_offsets[ulii] += cur_nei;
	    cur_word &= cur_word - 1;
	  }
	}
      }
      lptr = loadbuf;
      // iterate through heterozygous calls
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx += BITCT2) {
	cur_word = *lptr++;
	cur_word = (cur_word >> 1) & (~cur_word) & FIVEMASK;
	while (cur_word) {
	  het_cts[indiv_idx + CTZLU(cur_word) / 2] += 1;
	  cur_word &= cur_word - 1;
	}
      }
      nei_sum += cur_nei;
    }
    if (!x_variant_ct) {
      goto sexcheck_ret_NO_X_VAR;
    }
  }
  if (check_y) {
    if ((y_code != -1) && is_set(chrom_info_ptr->chrom_mask, (uint32_t)y_code)) {
      marker_uidx_end = chrom_info_ptr->chrom_end[(uint32_t)y_code];
      marker_uidx = next_unset_ul(marker_exclude, chrom_info_ptr->chrom_start[(uint32_t)y_code], marker_uidx_end);
      ytotal = marker_uidx_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, marker_uidx_end);
    }
    if (ytotal) {
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto sexcheck_ret_READ_FAIL;
      }
      for (ulii = 0; ulii < ytotal; ulii++, marker_uidx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto sexcheck_ret_READ_FAIL;
	  }
	}
	if (load_and_collapse(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf, indiv_ct, indiv_exclude, 0)) {
	  goto sexcheck_ret_READ_FAIL;
	}
	lptr = loadbuf;
	for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx += BITCT2) {
	  // iterate through missing calls
	  // vertical popcount would be faster, but we don't really care
	  cur_word = *lptr++;
	  cur_word = (~(cur_word >> 1)) & cur_word & FIVEMASK;
	  while (cur_word) {
	    ymiss_cts[indiv_idx + CTZLU(cur_word) / 2] += 1;
	    cur_word &= cur_word - 1;
	  }
	}
      }
    } else if (yonly) {
      logprint("Error: --check-sex/--impute-sex y-only requires Y chromosome data.\n");
      goto sexcheck_ret_INVALID_CMDLINE;
    } else {
      LOGPRINTF("Warning: No Y chromosome data for --%s-sex ycount.\n", do_impute? "impute" : "check");
    }
  }
  memcpy(outname_end, ".sexcheck", 10);
  if (fopen_checked(&outfile, outname, "w")) {
    goto sexcheck_ret_OPEN_FAIL;
  }
  sprintf(tbuf, "%%%us %%%us       PEDSEX       SNPSEX       STATUS%s%s\n", plink_maxfid, plink_maxiid, yonly? "" : "            F", check_y? "   YCOUNT" : "");
  fprintf(outfile, tbuf, "FID", "IID");
  indiv_uidx = 0;
  if (do_impute) {
    bitfield_andnot(sex_nm, indiv_exclude, unfiltered_indiv_ctl);
  }
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++, indiv_uidx++) {
    next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
    fid_ptr = &(person_ids[indiv_uidx * max_person_id_len]);
    iid_ptr = (char*)memchr(fid_ptr, '\t', max_person_id_len);
    wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(iid_ptr - fid_ptr), fid_ptr, tbuf);
    *wptr++ = ' ';
    wptr = fw_strcpy(plink_maxiid, &(iid_ptr[1]), wptr);
    if (!IS_SET(sex_nm, indiv_uidx)) {
      orig_sex_code = 0;
    } else {
      orig_sex_code = 2 - IS_SET(sex_male, indiv_uidx);
    }
    wptr = memseta(wptr, 32, 12);
    *wptr++ = '0' + orig_sex_code;
    wptr = memseta(wptr, 32, 12);
    imputed_sex_code = 0;
    if (!yonly) {
      ulii = x_variant_ct - missing_cts[indiv_idx];
      if (ulii) {
	dee = nei_sum - nei_offsets[indiv_idx];
	dtot = (double)((intptr_t)ulii) - dee;
	dff = (dtot - ((double)((int32_t)(het_cts[indiv_idx])))) / dtot;
	if (dff > check_sex_mthresh) {
	  if ((!check_y) || (ymiss_cts[indiv_idx] + min_m_yobs <= ytotal)) {
	    imputed_sex_code = 1;
	  }
	} else if ((dff < check_sex_fthresh) && ((!check_y) || (ymiss_cts[indiv_idx] + max_f_yobs >= ytotal))) {
	  imputed_sex_code = 2;
	}
        *wptr++ = '0' + imputed_sex_code;
        if (orig_sex_code && (orig_sex_code == imputed_sex_code)) {
	  wptr = memcpya(wptr, "           OK ", 14);
        } else {
	  wptr = memcpya(wptr, "      PROBLEM ", 14);
	  problem_ct++;
        }
        wptr = double_g_writewx4(wptr, dff, 12);
      } else {
        wptr = memcpya(wptr, "0      PROBLEM          nan", 27);
        problem_ct++;
      }
      if (check_y) {
	*wptr++ = ' ';
	wptr = uint32_writew8(wptr, ytotal - ymiss_cts[indiv_idx]);
      }
    } else {
      if (ymiss_cts[indiv_idx] + min_m_yobs <= ytotal) {
        imputed_sex_code = 1;
      } else if (ymiss_cts[indiv_idx] + max_f_yobs >= ytotal) {
        imputed_sex_code = 2;
      }
      *wptr++ = '0' + imputed_sex_code;
      if (orig_sex_code && (orig_sex_code == imputed_sex_code)) {
	wptr = memcpya(wptr, "           OK ", 14);
      } else {
	wptr = memcpya(wptr, "      PROBLEM ", 14);
	problem_ct++;
      }
      wptr = uint32_writew8(wptr, ytotal - ymiss_cts[indiv_idx]);
    }
    *wptr++ = '\n';
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto sexcheck_ret_WRITE_FAIL;
    }
    if (do_impute) {
      if (imputed_sex_code) {
	SET_BIT(sex_nm, indiv_uidx);
	if (imputed_sex_code == 1) {
	  SET_BIT(sex_male, indiv_uidx);
	} else {
	  CLEAR_BIT(sex_male, indiv_uidx);
	}
      } else {
	CLEAR_BIT(sex_nm, indiv_uidx);
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto sexcheck_ret_WRITE_FAIL;
  }
  if (do_impute) {
    bitfield_and(sex_male, sex_nm, unfiltered_indiv_ctl);
    gender_unk_ct = indiv_ct - popcount_longs(sex_nm, unfiltered_indiv_ctl);
    if (!gender_unk_ct) {
      LOGPREPRINTFWW("--impute-sex: %" PRIuPTR " Xchr and %" PRIuPTR " Ychr variant(s) scanned, all sexes imputed. Report written to %s .\n", x_variant_ct, ytotal, outname);
    } else {
      LOGPREPRINTFWW("--impute-sex: %" PRIuPTR " Xchr and %" PRIuPTR " Ychr variant(s) scanned, %" PRIuPTR "/%" PRIuPTR " sex%s imputed. Report written to %s .\n", x_variant_ct, ytotal, (indiv_ct - gender_unk_ct), indiv_ct, (indiv_ct == 1)? "" : "es", outname);
    }
    *gender_unk_ct_ptr = gender_unk_ct;
  } else {
    if (!problem_ct) {
      LOGPREPRINTFWW("--check-sex: %" PRIuPTR " Xchr and %" PRIuPTR " Ychr variant(s) scanned, no problems detected. Report written to %s .\n", x_variant_ct, ytotal, outname);
    } else {
      LOGPREPRINTFWW("--check-sex: %" PRIuPTR " Xchr and %" PRIuPTR " Ychr variant(s) scanned, %u problem%s detected. Report written to %s .\n", x_variant_ct, ytotal, problem_ct, (problem_ct == 1)? "" : "s", outname);
    }
  }
  logprintb();
  while (0) {
  sexcheck_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  sexcheck_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  sexcheck_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  sexcheck_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  sexcheck_ret_NO_X_VAR:
    logprint("Error: --check-sex/--impute-sex requires at least one polymorphic X chromosome\nsite.\n");
  sexcheck_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t write_snplist(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t list_23_indels) {
  FILE* outfile;
  uintptr_t marker_uidx = 0;
  uintptr_t markers_done = 0;
  int32_t retval = 0;
  const char* a0ptr = g_missing_geno_ptr;
  const char* adptr = &(g_one_char_strs[136]); // "D"
  const char* aiptr = &(g_one_char_strs[146]); // "I"
  char* a1ptr;
  char* a2ptr;
  char* cptr;
  char* cptr_end;
  uintptr_t marker_uidx_stop;
  if (!list_23_indels) {
    memcpy(outname_end, ".snplist", 9);
  } else {
    memcpy(outname_end, ".indel", 7);
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_snplist_ret_OPEN_FAIL;
  }
  if (!list_23_indels) {
    do {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      marker_uidx_stop = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
      markers_done += marker_uidx_stop - marker_uidx;
      cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
      cptr_end = &(marker_ids[marker_uidx_stop * max_marker_id_len]);
      marker_uidx = marker_uidx_stop;
      do {
	fputs(cptr, outfile);
	if (putc_checked('\n', outfile)) {
          goto write_snplist_ret_WRITE_FAIL;
	}
        cptr = &(cptr[max_marker_id_len]);
      } while (cptr < cptr_end);
    } while (markers_done < marker_ct);
  } else {
    for (; markers_done < marker_ct; marker_uidx++, markers_done++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      a1ptr = marker_allele_ptrs[2 * marker_uidx];
      a2ptr = marker_allele_ptrs[2 * marker_uidx + 1];
      if ((a1ptr != adptr) && (a1ptr != aiptr)) {
        if ((a1ptr != a0ptr) || ((a2ptr != adptr) && (a2ptr != aiptr))) {
	  continue;
	}
      } else if ((a2ptr != adptr) && (a2ptr != aiptr) && (a2ptr != a0ptr)) {
	continue;
      }
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      if (putc_checked('\n', outfile)) {
	goto write_snplist_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto write_snplist_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("List of %svariant IDs written to %s .\n", list_23_indels? "indel " : "" , outname);
  while (0) {
  write_snplist_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_snplist_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t het_report(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs) {
  // Same F coefficient computation as sexcheck().
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t indiv_ctl2 = (indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t monomorphic_ct = 0;
  uintptr_t marker_uidx = 0;
  double nei_sum = 0.0;
  uint32_t chrom_fo_idx = 0xffffffffU; // deliberate overflow
  uint32_t chrom_end = 0;
  int32_t retval = 0;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* lptr;
  uint32_t* het_cts;
  uint32_t* missing_cts;
  double* nei_offsets;
  char* fid_ptr;
  char* iid_ptr;
  char* wptr;
  double dpp;
  double dtot;
  double cur_nei;
  double dff;
  double dee;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t cur_missing_ct;
  uintptr_t allele_obs_ct;
  uintptr_t cur_word;
  uintptr_t ulii;
  uint32_t obs_ct;
  if (is_set(chrom_info_ptr->haploid_mask, 0)) {
    logprint("Error: --het cannot be used on haploid genomes.\n");
    goto het_report_ret_INVALID_CMDLINE;
  }
  if (wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&het_cts, indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&missing_cts, indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&nei_offsets, indiv_ct * sizeof(double))) {
    goto het_report_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_indiv_ctl2 - 1] = 0;
  loadbuf[indiv_ctl2 - 1] = 0;
  fill_uint_zero(het_cts, indiv_ct);
  fill_uint_zero(missing_cts, indiv_ct);
  fill_double_zero(nei_offsets, indiv_ct);
  marker_ct -= count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1, 1);
  if (!marker_ct) {
    goto het_report_ret_INVALID_CMDLINE;
  }
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
        goto het_report_ret_READ_FAIL;
      }
    }
    if (marker_uidx >= chrom_end) {
      do {
	do {
	  chrom_fo_idx++;
	} while (is_set(chrom_info_ptr->haploid_mask, chrom_info_ptr->chrom_file_order[chrom_fo_idx]));
	chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
	marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
      } while (marker_uidx >= chrom_end);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
        goto het_report_ret_READ_FAIL;
      }
    }
    if (load_and_collapse(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf, indiv_ct, indiv_exclude, 0)) {
      goto het_report_ret_READ_FAIL;
    }
    cur_missing_ct = count_01(loadbuf, indiv_ctl2);
    allele_obs_ct = 2 * (indiv_ct - cur_missing_ct);
    dpp = set_allele_freqs[marker_uidx];
    if ((!allele_obs_ct) || (dpp < 1e-8) || (dpp > (1 - 1e-8))) {
      monomorphic_ct++;
      continue;
    }
    cur_nei = 1.0 - 2 * dpp * (1 - dpp);
    if (cur_missing_ct) {
      // iterate through missing calls
      lptr = loadbuf;
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx += BITCT2) {
        cur_word = *lptr++;
        cur_word = (~(cur_word >> 1)) & cur_word & FIVEMASK;
        while (cur_word) {
          ulii = indiv_idx + CTZLU(cur_word) / 2;
          missing_cts[ulii] += 1;
          nei_offsets[ulii] += cur_nei;
          cur_word &= cur_word - 1;
	}
      }
    }
    lptr = loadbuf;
    // iterate through heterozygous calls
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx += BITCT2) {
      cur_word = *lptr++;
      cur_word = (cur_word >> 1) & (~cur_word) & FIVEMASK;
      while (cur_word) {
        het_cts[indiv_idx + CTZLU(cur_word) / 2] += 1;
        cur_word &= cur_word - 1;
      }
    }
    nei_sum += cur_nei;
  }
  marker_ct -= monomorphic_ct;
  if (!marker_ct) {
    goto het_report_ret_INVALID_CMDLINE;
  }
  memcpy(outname_end, ".het", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto het_report_ret_OPEN_FAIL;
  }
  sprintf(tbuf, "%%%us %%%us       O(HOM)       E(HOM)        N(NM)            F\n", plink_maxfid, plink_maxiid);
  fprintf(outfile, tbuf, "FID", "IID");
  indiv_uidx = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++, indiv_uidx++) {
    next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
    fid_ptr = &(person_ids[indiv_uidx * max_person_id_len]);
    iid_ptr = (char*)memchr(fid_ptr, '\t', max_person_id_len);
    wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(iid_ptr - fid_ptr), fid_ptr, tbuf);
    *wptr++ = ' ';
    wptr = fw_strcpy(plink_maxiid, &(iid_ptr[1]), wptr);
    wptr = memseta(wptr, 32, 3);
    obs_ct = marker_ct - missing_cts[indiv_idx];
    if (obs_ct) {
      wptr = uint32_writew10x(wptr, obs_ct - het_cts[indiv_idx], ' ');
      dee = nei_sum - nei_offsets[indiv_idx];
      wptr = double_g_writewx4(wptr, dee, 12);
      wptr = memseta(wptr, 32, 3);
      wptr = uint32_writew10x(wptr, obs_ct, ' ');
      dtot = (double)((int32_t)obs_ct) - dee;
      dff = (dtot - ((double)((int32_t)(het_cts[indiv_idx])))) / dtot;
      wptr = double_g_writewx4x(wptr, dff, 12, '\n');
    } else {
      wptr = memcpya(wptr, "         0            0            0          nan\n", 50);
    }
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto het_report_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto het_report_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--het: %" PRIuPTR " variant%s scanned, report written to %s .\n", marker_ct, (marker_ct == 1)? "" : "s", outname);
  while (0) {
  het_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  het_report_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  het_report_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  het_report_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  het_report_ret_INVALID_CMDLINE:
    logprint("Error: --het requires at least one polymorphic autosomal marker.\n");
    retval = RET_INVALID_CMDLINE;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t score_report(Score_info* sc_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, double* set_allele_freqs, uintptr_t indiv_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t hh_exists, Chrom_info* chrom_info_ptr, char* outname, char* outname_end) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  FILE* outfile = NULL;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t indiv_ctl2 = (indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t topsize = 0;
  uintptr_t miss_ct = 0;
  uintptr_t range_ct = 0;
  uintptr_t range_skip = 0;
  uintptr_t ulii = 0;
  char* tbuf2 = &(tbuf[MAXLINELEN]);
  uintptr_t* marker_exclude_main = NULL;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  double* qrange_keys = NULL;
  double* effect_sizes_cur = NULL;
  double ploidy_d = 0.0;
  double lbound = 0.0;
  double ubound = 0.0;
  uint32_t modifier = sc_ip->modifier;
  uint32_t center = modifier & SCORE_CENTER;
  uint32_t report_average = !(modifier & SCORE_SUM);
  uint32_t mean_impute = !(modifier & SCORE_NO_MEAN_IMPUTATION);
  uint32_t is_haploid = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t ploidy = 0;
  uint32_t max_rangename_len = 0;
  uint32_t rangename_len = 0;
  int32_t retval = 0;
  double female_effect_size[4];
  int32_t female_allele_ct_delta[4];
  char missing_pheno_str[32];
  char* bufptr_arr[3];
  uintptr_t* marker_exclude;
  uintptr_t* a2_effect;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* lbptr;
  char* bufptr;
  char* loadbuf_c;
  char* sorted_marker_ids;
  uint32_t* marker_id_map;
  uint32_t* miss_cts;
  int32_t* named_allele_ct_deltas;
  double* score_deltas;
  double* effect_sizes;
  double* dptr;
  double cur_effect_size;
  double half_effect_size;
  double missing_effect;
  double cur_effect_offset;
  double score_base;
  double female_y_offset;
  double dxx;
  uintptr_t loadbuf_size;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t line_idx;
  uintptr_t uljj;
  uint32_t missing_pheno_len;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t obs_expected;
  uint32_t named_allele_ct_expected;
  uint32_t cur_marker_ct;
  uint32_t first_col_m1;
  uint32_t col_01_delta;
  uint32_t col_12_delta;
  uint32_t varid_idx;
  uint32_t allele_idx;
  uint32_t effect_idx;
  uint32_t named_allele_ct_female_delta;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  int32_t obs_expected_female_delta;
  int32_t delta1;
  int32_t delta2;
  int32_t deltam;
  int32_t sorted_idx;
  sorted_marker_ids = (char*)top_alloc(&topsize, marker_ct * max_marker_id_len);
  if (!sorted_marker_ids) {
    goto score_report_ret_NOMEM;
  }
  marker_id_map = (uint32_t*)top_alloc(&topsize, marker_ct * sizeof(int32_t));
  if (!marker_id_map) {
    goto score_report_ret_NOMEM;
  }
  dptr = (double*)top_alloc(&topsize, unfiltered_marker_ct * sizeof(double));
  if (!dptr) {
    goto score_report_ret_NOMEM;
  }
  wkspace_left -= topsize;
  retval = sort_item_ids_noalloc(sorted_marker_ids, marker_id_map, unfiltered_marker_ct, marker_exclude_orig, marker_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    wkspace_left += topsize;
    goto score_report_ret_1;
  }
  if (wkspace_alloc_ul_checked(&marker_exclude, unfiltered_marker_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&a2_effect, unfiltered_marker_ctl * sizeof(intptr_t))) {
    goto score_report_ret_NOMEM2;
  }
  fill_all_bits(marker_exclude, unfiltered_marker_ct);
  fill_ulong_zero(a2_effect, unfiltered_marker_ctl);
  loadbuf_size = wkspace_left;
  wkspace_left += topsize;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto score_report_ret_NOMEM;
  }
  loadbuf_c = (char*)wkspace_base;
  retval = open_and_load_to_first_token(&infile, sc_ip->fname, loadbuf_size, '\0', "--score file", loadbuf_c, &bufptr, &line_idx);
  if (retval) {
    goto score_report_ret_1;
  }
  // bah, may as well brute force the six cases, I guess
  if (sc_ip->varid_col < sc_ip->allele_col) {
    first_col_m1 = sc_ip->varid_col;
    varid_idx = 0;
    if (sc_ip->allele_col < sc_ip->effect_col) {
      col_01_delta = sc_ip->allele_col - first_col_m1;
      col_12_delta = sc_ip->effect_col - sc_ip->allele_col;
      allele_idx = 1;
      effect_idx = 2;
    } else {
      allele_idx = 2;
      if (sc_ip->varid_col < sc_ip->effect_col) {
        col_01_delta = sc_ip->effect_col - first_col_m1;
        col_12_delta = sc_ip->allele_col - sc_ip->effect_col;
	effect_idx = 1;
      } else {
        first_col_m1 = sc_ip->effect_col;
        col_01_delta = sc_ip->varid_col - first_col_m1;
        col_12_delta = sc_ip->allele_col - sc_ip->varid_col;
        varid_idx = 1;
        effect_idx = 0;
      }
    }
  } else {
    first_col_m1 = sc_ip->allele_col;
    allele_idx = 0;
    if (sc_ip->varid_col < sc_ip->effect_col) {
      col_01_delta = sc_ip->varid_col - first_col_m1;
      col_12_delta = sc_ip->effect_col - sc_ip->varid_col;
      varid_idx = 1;
      effect_idx = 2;
    } else {
      varid_idx = 2;
      if (sc_ip->allele_col < sc_ip->effect_col) {
        col_01_delta = sc_ip->effect_col - first_col_m1;
        col_12_delta = sc_ip->varid_col - sc_ip->effect_col;
	effect_idx = 1;
      } else {
        first_col_m1 = sc_ip->effect_col;
        col_01_delta = sc_ip->allele_col - first_col_m1;
        col_12_delta = sc_ip->varid_col - sc_ip->allele_col;
        allele_idx = 1;
        effect_idx = 0;
      }
    }
  }
  first_col_m1--;
  if (modifier & SCORE_HEADER) {
    goto score_report_load_next;
  }
  while (1) {
    bufptr_arr[0] = next_token_multz(bufptr, first_col_m1);
    bufptr_arr[1] = next_token_mult(bufptr_arr[0], col_01_delta);
    bufptr_arr[2] = next_token_mult(bufptr_arr[1], col_12_delta);
    if (!bufptr_arr[2]) {
      goto score_report_ret_MISSING_TOKENS;
    }
    sorted_idx = bsearch_str(bufptr_arr[varid_idx], strlen_se(bufptr_arr[varid_idx]), sorted_marker_ids, max_marker_id_len, marker_ct);
    if (sorted_idx != -1) {
      marker_uidx = marker_id_map[(uint32_t)sorted_idx];
      bufptr_arr[allele_idx][strlen_se(bufptr_arr[allele_idx])] = '\0';
      uii = strcmp(bufptr_arr[allele_idx], marker_allele_ptrs[2 * marker_uidx]);
      if ((!uii) || (!strcmp(bufptr_arr[allele_idx], marker_allele_ptrs[2 * marker_uidx + 1]))) {
        if (scan_double(bufptr_arr[effect_idx], &(dptr[marker_uidx]))) {
          miss_ct++;
	} else {
	  if (!IS_SET(marker_exclude, marker_uidx)) {
            bufptr_arr[varid_idx][strlen_se(bufptr_arr[varid_idx])] = '\0';
            LOGPREPRINTFWW("Error: Duplicate variant '%s' in --score file.\n", bufptr_arr[varid_idx]);
            goto score_report_ret_INVALID_FORMAT_2;
	  }
          CLEAR_BIT(marker_exclude, marker_uidx);
	  if (uii) {
	    SET_BIT(a2_effect, marker_uidx);
	  }
	}
      } else {
	miss_ct++;
      }
    } else {
      miss_ct++;
    }
  score_report_load_next:
    line_idx++;
    if (!fgets(loadbuf_c, loadbuf_size, infile)) {
      break;
    }
    if (!(loadbuf_c[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        sprintf(logbuf, "Error: Line %" PRIuPTR " of --score file is pathologically long.\n", line_idx);
        goto score_report_ret_INVALID_FORMAT_2;
      }
      goto score_report_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf_c);
    if (is_eoln_kns(*bufptr)) {
      goto score_report_load_next;
    }
  }
  if (fclose_null(&infile)) {
    goto score_report_ret_READ_FAIL;
  }
  cur_marker_ct = unfiltered_marker_ct - popcount_longs(marker_exclude, unfiltered_marker_ctl);
  if (!cur_marker_ct) {
    logprint("Error: No valid entries in --score file.\n");
    goto score_report_ret_INVALID_FORMAT;
  } else if (cur_marker_ct >= 0x40000000) {
    logprint("Error: --score does not support >= 2^30 markers.\n");
    goto score_report_ret_INVALID_FORMAT;
  }
  if (miss_ct) {
    LOGPRINTF("Warning: %" PRIuPTR " line%s skipped in --score file.\n", miss_ct, (miss_ct == 1)? "" : "s");
  }
  if (sc_ip->data_fname) {
    effect_sizes_cur = dptr; // not collapsed yet
    ulii = topsize;
    dptr = (double*)top_alloc(&topsize, unfiltered_marker_ct * sizeof(double));
    if (!dptr) {
      goto score_report_ret_NOMEM;
    }
    wkspace_left -= topsize;
    if (wkspace_alloc_ul_checked(&marker_exclude_main, unfiltered_marker_ctl * sizeof(intptr_t))) {
      goto score_report_ret_NOMEM2;
    }
    fill_all_bits(marker_exclude_main, unfiltered_marker_ct);
    wkspace_left += topsize;
    loadbuf_size = wkspace_left - topsize;
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size <= MAXLINELEN) {
      goto score_report_ret_NOMEM;
    }
    loadbuf_c = (char*)wkspace_base;
    retval = open_and_load_to_first_token(&infile, sc_ip->data_fname, loadbuf_size, '\0', "--q-score-range data file", loadbuf_c, &bufptr, &line_idx);
    if (retval) {
      goto score_report_ret_1;
    }
    if (sc_ip->data_varid_col < sc_ip->data_col) {
      first_col_m1 = sc_ip->data_varid_col;
      col_01_delta = sc_ip->data_col - first_col_m1;
      varid_idx = 0;
    } else {
      first_col_m1 = sc_ip->data_col;
      col_01_delta = sc_ip->data_varid_col - first_col_m1;
      varid_idx = 1;
    }
    first_col_m1--;
    miss_ct = 0;
    if (modifier & SCORE_DATA_HEADER) {
      goto score_report_qrange_load_next;
    }
    while (1) {
      bufptr_arr[0] = next_token_multz(bufptr, first_col_m1);
      bufptr_arr[1] = next_token_mult(bufptr_arr[0], col_01_delta);
      if (!bufptr_arr[1]) {
        goto score_report_ret_MISSING_TOKENS_Q;
      }
      sorted_idx = bsearch_str(bufptr_arr[varid_idx], strlen_se(bufptr_arr[varid_idx]), sorted_marker_ids, max_marker_id_len, marker_ct);
      if (sorted_idx != -1) {
	marker_uidx = marker_id_map[(uint32_t)sorted_idx];
        if (!IS_SET(marker_exclude, marker_uidx)) {
	  if (scan_double(bufptr_arr[1 - varid_idx], &dxx) || (dxx != dxx)) {
	    miss_ct++;
	  } else {
	    dptr[marker_uidx] = dxx;
	    if (!IS_SET(marker_exclude_main, marker_uidx)) {
	      bufptr_arr[varid_idx][strlen_se(bufptr_arr[varid_idx])] = '\0';
	      LOGPREPRINTFWW("Error: Duplicate variant '%s' in --q-score-range data file.\n", bufptr_arr[varid_idx]);
	      goto score_report_ret_INVALID_FORMAT_2;
	    }
            CLEAR_BIT(marker_exclude_main, marker_uidx);
	  }
	} else {
	  miss_ct++;
	}
      } else {
	miss_ct++;
      }
    score_report_qrange_load_next:
      line_idx++;
      if (!fgets(loadbuf_c, loadbuf_size, infile)) {
	break;
      }
      if (!(loadbuf_c[loadbuf_size - 1])) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  sprintf(logbuf, "Error: Line %" PRIuPTR " of --q-score-range data file is pathologically long.\n", line_idx);
	  goto score_report_ret_INVALID_FORMAT_2;
	}
	goto score_report_ret_NOMEM;
      }
      bufptr = skip_initial_spaces(loadbuf_c);
      if (is_eoln_kns(*bufptr)) {
	goto score_report_qrange_load_next;
      }
    }
    if (fclose_null(&infile)) {
      goto score_report_ret_READ_FAIL;
    }
    marker_ct = unfiltered_marker_ct - popcount_longs(marker_exclude_main, unfiltered_marker_ctl);
    if (!marker_ct) {
      logprint("Error: No valid entries in --q-score-range data file.\n");
      goto score_report_ret_INVALID_FORMAT;
    }
    wkspace_left -= topsize;
    qrange_keys = (double*)alloc_and_init_collapsed_arr((char*)dptr, sizeof(double), unfiltered_marker_ct, marker_exclude_main, marker_ct, 0);
    wkspace_left += topsize;
    if (!qrange_keys) {
      goto score_report_ret_NOMEM;
    }
    topsize = ulii;
    wkspace_left -= topsize;
    effect_sizes = (double*)alloc_and_init_collapsed_arr((char*)effect_sizes_cur, sizeof(double), unfiltered_marker_ct, marker_exclude_main, marker_ct, 0);
    wkspace_left += topsize;
    if (!effect_sizes) {
      goto score_report_ret_NOMEM;
    }
    if (wkspace_alloc_d_checked(&effect_sizes_cur, marker_ct * sizeof(double))) {
      goto score_report_ret_NOMEM;
    }
    if (miss_ct) {
      LOGPRINTF("Warning: %" PRIuPTR " line%s skipped in --q-score-range data file.\n", miss_ct, (miss_ct == 1)? "" : "s");
    }
    miss_ct = 0;
    if (fopen_checked(&infile, sc_ip->range_fname, "r")) {
      goto score_report_ret_OPEN_FAIL;
    }
    max_rangename_len = (FNAMESIZE - 10) - ((uintptr_t)(outname_end - outname));
    tbuf[MAXLINELEN - 1] = ' ';
    *outname_end = '.';
  } else {
    wkspace_left -= topsize;
    effect_sizes = (double*)alloc_and_init_collapsed_arr((char*)dptr, sizeof(double), unfiltered_marker_ct, marker_exclude, cur_marker_ct, 0);
    wkspace_left += topsize;
    if (!effect_sizes) {
      goto score_report_ret_NOMEM;
    }
  }
  // topsize = 0;
  if (wkspace_alloc_d_checked(&score_deltas, indiv_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&miss_cts, indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&named_allele_ct_deltas, indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, indiv_ctl2 * sizeof(intptr_t))) {
    goto score_report_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_indiv_ctl2 - 1] = 0;
  loadbuf[indiv_ctl2 - 1] = 0;
  // force indiv_male_include2 allocation
  if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, indiv_ct, hh_exists | XMHH_EXISTS, 0, indiv_exclude, sex_male, &indiv_include2, &indiv_male_include2)) {
    goto score_report_ret_NOMEM;
  }
  missing_pheno_len = strlen(output_missing_pheno);
  if (missing_pheno_len < 6) {
    memset(missing_pheno_str, 32, 6 - missing_pheno_len);
    memcpy(&(missing_pheno_str[6 - missing_pheno_len]), output_missing_pheno, missing_pheno_len);
    missing_pheno_len = 6;
  } else {
    memcpy(missing_pheno_str, output_missing_pheno, missing_pheno_len);
  }
  line_idx = 0;
  if (marker_exclude_main) {
  score_report_qrange_next:
    while (1) {
      line_idx++;
      if (!fgets(tbuf, MAXLINELEN, infile)) {
	if (fclose_null(&infile)) {
	  goto score_report_ret_READ_FAIL;
	}
	*outname_end = '\0';
	LOGPRINTF("--score: %" PRIuPTR " range%s processed", range_ct, (range_ct == 1)? "" : "s");
	if (range_skip) {
	  LOGPRINTF(" (%" PRIuPTR " empty range%s skipped)", range_skip, (range_skip == 1)? "" : "s");
	}
	logprint(".\n");
	LOGPREPRINTFWW("Results written to %s.*.profile.\n", outname);
	fputs(logbuf, stdout);
	goto score_report_ret_1;
      }
      if (!tbuf[MAXLINELEN - 1]) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --q-score-range range file is pathologically long.\n", line_idx);
	goto score_report_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      rangename_len = strlen_se(bufptr);
      bufptr_arr[1] = skip_initial_spaces(&(bufptr[rangename_len]));
      bufptr_arr[2] = next_token(bufptr_arr[1]);
      if ((!bufptr_arr[2]) || scan_double(bufptr_arr[1], &lbound) || scan_double(bufptr_arr[2], &ubound) || (lbound != lbound) || (ubound != ubound) || (lbound > ubound)) {
	continue;
      }
      if (rangename_len > max_rangename_len) {
	sprintf(logbuf, "Error: Excessively long range name on line %" PRIuPTR " of --q-score-range range\nfile.\n", line_idx);
	goto score_report_ret_INVALID_FORMAT_2;
      }
      bufptr_arr[0] = bufptr;
      break;
    }
    fill_all_bits(marker_exclude, unfiltered_marker_ct);
    marker_uidx = next_unset_unsafe(marker_exclude_main, 0);
    dptr = effect_sizes_cur;
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude_main, &marker_uidx);
      dxx = qrange_keys[marker_idx];
      if ((dxx >= lbound) && (dxx <= ubound)) {
	CLEAR_BIT(marker_exclude, marker_uidx);
	*dptr++ = effect_sizes[marker_idx];
      }
    }
    cur_marker_ct = unfiltered_marker_ct - popcount_longs(marker_exclude, unfiltered_marker_ctl);
    if (!cur_marker_ct) {
      range_skip++;
      goto score_report_qrange_next;
    }
    range_ct++;
    dptr = effect_sizes_cur;
  } else {
    dptr = effect_sizes;
  }
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
    goto score_report_ret_READ_FAIL;
  }
  fill_double_zero(score_deltas, indiv_ct);
  fill_uint_zero(miss_cts, indiv_ct);
  fill_int_zero(named_allele_ct_deltas, indiv_ct);
  score_base = 0.0;
  female_y_offset = 0.0;
  obs_expected = 0;
  named_allele_ct_expected = 0;
  named_allele_ct_female_delta = 0;
  obs_expected_female_delta = 0;
  chrom_fo_idx = 0xffffffffU;
  chrom_end = 0;
  for (marker_idx = 0; marker_idx < cur_marker_ct; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto score_report_ret_READ_FAIL;
      }
    }
    if (marker_uidx >= chrom_end) {
      do {
	chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1];
      } while (marker_uidx >= chrom_end);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      is_haploid = IS_SET(chrom_info_ptr->haploid_mask, uii);
      is_x = ((int32_t)uii == chrom_info_ptr->x_code)? 1 : 0;
      is_y = ((int32_t)uii == chrom_info_ptr->y_code)? 1 : 0;
      ploidy = 2 - is_haploid;
      ploidy_d = (double)((int32_t)ploidy);
    }
    if (load_and_collapse(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf, indiv_ct, indiv_exclude, IS_SET(marker_reverse, marker_uidx))) {
      goto score_report_ret_READ_FAIL;
    }
    if (is_haploid && hh_exists) {
      haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_ct, is_x, is_y, (unsigned char*)loadbuf);
    }
    cur_effect_size = (*dptr++) * ploidy_d;
    uii = IS_SET(a2_effect, marker_uidx);
    if (!uii) {
      delta1 = 1;
      delta2 = ploidy;
      deltam = 0;
    } else {
      if (!center) {
	score_base += cur_effect_size;
      }
      cur_effect_size = -cur_effect_size;
      delta1 = -1;
      delta2 = -((int32_t)ploidy);
      deltam = delta2;
      named_allele_ct_expected += ploidy;
    }
    dxx = (1.0 - set_allele_freqs[marker_uidx]) * cur_effect_size;
    missing_effect = dxx;
    if (center) {
      score_base -= dxx;
    } else if (!mean_impute) {
      if (!uii) {
	missing_effect = 0;
      } else {
	missing_effect = cur_effect_size;
      }
    }
    half_effect_size = cur_effect_size * 0.5;
    obs_expected += ploidy;
    lbptr = loadbuf;
    if ((!is_x) && (!is_y)) {
      uii = 0;
      do {
	ulii = ~(*lbptr++);
	if (uii + BITCT2 > indiv_ct) {
	  ulii &= (ONELU << ((indiv_ct & (BITCT2 - 1)) * 2)) - ONELU;
	}
	while (ulii) {
	  ujj = CTZLU(ulii) & (BITCT - 2);
	  ukk = (ulii >> ujj) & 3;
	  indiv_idx = uii + (ujj / 2);
	  if (ukk == 1) {
	    score_deltas[indiv_idx] += half_effect_size;
	    named_allele_ct_deltas[indiv_idx] += delta1;
	  } else if (ukk == 3) {
	    score_deltas[indiv_idx] += cur_effect_size;
	    named_allele_ct_deltas[indiv_idx] += delta2;
	  } else {
	    miss_cts[indiv_idx] += ploidy;
	    named_allele_ct_deltas[indiv_idx] += deltam;
	    score_deltas[indiv_idx] += missing_effect;
	  }
	  ulii &= ~((3 * ONELU) << ujj);
	}
	uii += BITCT2;
      } while (uii < indiv_ct);
    } else {
      // This could be more efficient if we kept male and female versions of
      // score_base, but not really worth the trouble.
      // deltas and effect_size variables are currently configured for males,
      // so we need to compute female versions here.
      if (center) {
	// this value is positive if A2 named, negative if A1 named (assuming
	// positive effect size)
	cur_effect_offset = -dxx;
      } else if (!uii) {
	cur_effect_offset = 0;
      } else {
	// this value is positive
	cur_effect_offset = -cur_effect_size;
      }
      if (is_x) {
	obs_expected_female_delta += 1;
	female_effect_size[0] = cur_effect_offset;
	female_effect_size[1] = cur_effect_offset + cur_effect_size;
	female_effect_size[2] = cur_effect_offset + 2 * missing_effect;
	female_effect_size[3] = cur_effect_offset + 2 * cur_effect_size;
	if (!uii) {
	  female_allele_ct_delta[0] = 0;
	  female_allele_ct_delta[1] = 1;
	  female_allele_ct_delta[2] = 0;
	  female_allele_ct_delta[3] = 2;
	} else {
	  female_allele_ct_delta[0] = 1;
	  female_allele_ct_delta[1] = 0;
	  female_allele_ct_delta[2] = -1;
	  female_allele_ct_delta[3] = -1;
	}
      } else {
	obs_expected_female_delta -= 1;
	female_y_offset -= cur_effect_offset;
	if (uii) {
	  named_allele_ct_female_delta += 1;
	}
      }
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	if (!(indiv_idx & (BITCT2 - 1))) {
	  ulii = ~(*lbptr++);
	}
	uljj = ulii & 3;
	if (IS_SET_DBL(indiv_male_include2, indiv_idx)) {
	  if (uljj) {
	    if (uljj == 3) {
	      score_deltas[indiv_idx] += cur_effect_size;
	      named_allele_ct_deltas[indiv_idx] += delta2;
	    } else {
	      miss_cts[indiv_idx] += 1;
	      named_allele_ct_deltas[indiv_idx] += deltam;
	      score_deltas[indiv_idx] += missing_effect;
	    }
	  }
	} else if (is_x) {
	  score_deltas[indiv_idx] += female_effect_size[uljj];
	  named_allele_ct_deltas[indiv_idx] += female_allele_ct_delta[uljj];
	  if (uljj == 2) {
	    miss_cts[indiv_idx] += 2;
	  }
	}
	ulii >>= 2;
      }
    }
  }
  if (marker_exclude_main) {
    bufptr = memcpya(&(outname_end[1]), bufptr_arr[0], rangename_len);
    memcpy(bufptr, ".profile", 9);
  } else {
    memcpy(outname_end, ".profile", 9);
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto score_report_ret_OPEN_FAIL;
  }
  sprintf(tbuf2, "%%%us %%%us  PHENO    CNT   CNT2 %s\n", plink_maxfid, plink_maxiid, report_average? "   SCORE" : "SCORESUM");
  fprintf(outfile, tbuf2, "FID", "IID");
  for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
    next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
    bufptr_arr[0] = &(person_ids[indiv_uidx * max_person_id_len]);
    uii = strlen_se(bufptr_arr[0]);
    bufptr = fw_strcpyn(plink_maxfid, uii, bufptr_arr[0], tbuf2);
    *bufptr++ = ' ';
    bufptr = fw_strcpy(plink_maxiid, &(bufptr_arr[0][uii + 1]), bufptr);
    *bufptr++ = ' ';
    if (IS_SET(pheno_nm, indiv_uidx)) {
      if (pheno_c) {
	bufptr = memseta(bufptr, 32, 5);
	*bufptr++ = '1' + IS_SET(pheno_c, indiv_uidx);
      } else {
	bufptr = width_force(6, bufptr, double_g_write(bufptr, pheno_d[indiv_uidx]));
      }
    } else {
      bufptr = memcpya(bufptr, missing_pheno_str, missing_pheno_len);
    }
    *bufptr++ = ' ';
    ujj = 1 - IS_SET_DBL(indiv_male_include2, indiv_idx); // female?
    uii = obs_expected + ((int32_t)ujj) * obs_expected_female_delta - miss_cts[indiv_idx];
    bufptr = uint32_writew6x(bufptr, uii, ' ');
    if (mean_impute) {
      uii += miss_cts[indiv_idx];
    }
    bufptr = uint32_writew6x(bufptr, ((int32_t)named_allele_ct_expected) - ujj * named_allele_ct_female_delta + named_allele_ct_deltas[indiv_idx], ' ');
    dxx = (score_base + ((int32_t)ujj) * female_y_offset + score_deltas[indiv_idx]);
    if (fabs(dxx) < SMALL_EPSILON) {
      dxx = 0;
    } else if (report_average) {
      dxx /= ((double)((int32_t)uii));
    }
    bufptr = width_force(8, bufptr, double_g_write(bufptr, dxx));
    *bufptr++ = '\n';
    if (fwrite_checked(tbuf2, bufptr - tbuf2, outfile)) {
      goto score_report_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto score_report_ret_WRITE_FAIL;
  }
  if (marker_exclude_main) {
    LOGPREPRINTFWW("%s written.\n", outname);
    logstr(logbuf);
    goto score_report_qrange_next;
  }
  LOGPRINTFWW("--score: Results written to %s .\n", outname);
  while (0) {
  score_report_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  score_report_ret_NOMEM2:
    wkspace_left += topsize;
  score_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  score_report_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  score_report_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  score_report_ret_MISSING_TOKENS_Q:
    LOGPRINTF("Error: Line %" PRIuPTR " of --q-score-range data file has fewer tokens than\nexpected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  score_report_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --score file has fewer tokens than expected.\n", line_idx);
  score_report_ret_INVALID_FORMAT_2:
    logprintb();
  score_report_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 score_report_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  fclose_cond(outfile);
  return retval;
}
