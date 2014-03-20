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

int32_t score_report(Score_info* sc_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, double* set_allele_freqs, uintptr_t indiv_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, int32_t missing_pheno, uint32_t missing_pheno_len, uint32_t hh_exists, Chrom_info* chrom_info_ptr, char* outname, char* outname_end) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  FILE* outfile = NULL;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t indiv_ctl2 = (indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t range_ct = 1;
  uintptr_t topsize = 0;
  uintptr_t miss_ct = 0;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  double score_base = 0.0;
  double ploidy_d = 0.0;
  uint32_t modifier = sc_ip->modifier;
  uint32_t center = modifier & SCORE_CENTER;
  uint32_t report_average = !(modifier & SCORE_SUM);
  uint32_t mean_impute = !(modifier & SCORE_NO_MEAN_IMPUTATION);
  uint32_t obs_expected = 0;
  uint32_t named_allele_ct_expected = 0;
  uint32_t chrom_fo_idx = 0xffffffffU; // deliberate overflow
  uint32_t chrom_end = 0;
  uint32_t is_haploid = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t ploidy = 0;
  int32_t retval = 0;
  char missing_pheno_str[12];
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
  double dxx;
  uintptr_t topsize_bak;
  uintptr_t loadbuf_size;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t range_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t ulii;
  uint32_t cur_marker_ct;
  uint32_t first_col_m1;
  uint32_t col_01_delta;
  uint32_t col_12_delta;
  uint32_t varid_idx;
  uint32_t allele_idx;
  uint32_t effect_idx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  int32_t delta1;
  int32_t delta2;
  int32_t sorted_idx;
  sorted_marker_ids = (char*)top_alloc(&topsize, marker_ct * max_marker_id_len);
  if (!sorted_marker_ids) {
    goto score_report_ret_NOMEM;
  }
  marker_id_map = (uint32_t*)top_alloc(&topsize, marker_ct * sizeof(int32_t));
  if (!marker_id_map) {
    goto score_report_ret_NOMEM;
  }
  topsize_bak = topsize;
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
  retval = open_and_load_to_first_token(&infile, sc_ip->fname, loadbuf_size, '\0', "--score file", loadbuf_c, &bufptr);
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
    if (first_col_m1) {
      bufptr_arr[0] = next_item_mult(bufptr, first_col_m1);
    } else {
      bufptr_arr[0] = bufptr;
    }
    bufptr_arr[1] = next_item_mult(bufptr_arr[0], col_01_delta);
    bufptr_arr[2] = next_item_mult(bufptr_arr[1], col_12_delta);
    if (!bufptr_arr[2]) {
      logprint("Error: Missing token(s) in --score file.\n");
      goto score_report_ret_INVALID_FORMAT;
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
            sprintf(logbuf, "Error: Duplicate variant '%s' in --score file.\n", bufptr_arr[varid_idx]);
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
    if (!fgets(loadbuf_c, loadbuf_size, infile)) {
      break;
    }
    if (!(loadbuf_c[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        LOGPRINTF("Error: Pathologically long line in --score file.\n");
        goto score_report_ret_INVALID_FORMAT;
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
  wkspace_left -= topsize;
  effect_sizes = (double*)alloc_and_init_collapsed_arr((char*)dptr, sizeof(double), unfiltered_marker_ct, marker_exclude, cur_marker_ct, 0);
  wkspace_left += topsize;
  if (!effect_sizes) {
    goto score_report_ret_NOMEM;
  }
  topsize = topsize_bak;
  if (sc_ip->range_fname) {
    logprint("Error: --q-score-range is currently under development.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto score_report_ret_1;
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
  fill_double_zero(score_deltas, indiv_ct);
  fill_uint_zero(miss_cts, indiv_ct);
  fill_int_zero(named_allele_ct_deltas, indiv_ct);
  if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, indiv_ct, hh_exists, 0, indiv_exclude, sex_male, &indiv_include2, &indiv_male_include2)) {
    goto score_report_ret_NOMEM;
  }
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
    goto score_report_ret_READ_FAIL;
  }
  dptr = effect_sizes;
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
      if (is_x) {
	logprint("Error: --score X chromosome handling is currently under development.\n");
        retval = RET_CALC_NOT_YET_SUPPORTED;
        goto score_report_ret_1;
      }
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
    if (!IS_SET(a2_effect, marker_uidx)) {
      dxx = 1.0 - set_allele_freqs[marker_uidx];
      delta1 = 1;
      delta2 = ploidy;
    } else {
      dxx = set_allele_freqs[marker_uidx];
      score_base += cur_effect_size;
      cur_effect_size = -cur_effect_size;
      delta1 = -1;
      delta2 = -((int32_t)ploidy);
      named_allele_ct_expected += ploidy;
    }
    dxx *= cur_effect_size;
    if (center) {
      score_base -= dxx;
    }
    half_effect_size = cur_effect_size * 0.5;
    obs_expected += ploidy;
    lbptr = loadbuf;
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
	  if (mean_impute) {
	    score_deltas[indiv_idx] += dxx;
	  }
	}
        ulii &= ~((3 * ONELU) << ujj);
      }
      uii += BITCT2;
    } while (uii < indiv_ct);
  }
  if (missing_pheno_len < 6) {
    bufptr = memseta(missing_pheno_str, 32, 6 - missing_pheno_len);
    missing_pheno_len = 6;
  } else {
    bufptr = missing_pheno_str;
  }
  int32_write(bufptr, missing_pheno);
  for (range_idx = 0; range_idx < range_ct; range_idx++) {
    if (0) {
    } else {
      memcpy(outname_end, ".profile", 9);
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto score_report_ret_OPEN_FAIL;
    }
    sprintf(tbuf, "%%%us %%%us  PHENO    CNT   CNT2    SCORE", plink_maxfid, plink_maxiid);
    fprintf(outfile, tbuf, "FID", "IID");
    if (!report_average) {
      fputs("SUM", outfile);
    }
    putc('\n', outfile);
    for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      bufptr_arr[0] = &(person_ids[indiv_uidx * max_person_id_len]);
      uii = strlen_se(bufptr_arr[0]);
      bufptr = fw_strcpyn(plink_maxfid, uii, bufptr_arr[0], tbuf);
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
      uii = obs_expected - miss_cts[indiv_idx];
      bufptr = uint32_writew6x(bufptr, uii, ' ');
      bufptr = uint32_writew6x(bufptr, ((int32_t)named_allele_ct_expected) + named_allele_ct_deltas[indiv_idx], ' ');
      if (report_average) {
        bufptr = width_force(8, bufptr, double_g_write(bufptr, (score_base + score_deltas[indiv_idx]) / ((double)((int32_t)uii))));
      } else {
        bufptr = width_force(8, bufptr, double_g_write(bufptr, score_base + score_deltas[indiv_idx]));
      }
      *bufptr++ = '\n';
      if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	goto score_report_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto score_report_ret_WRITE_FAIL;
    }
  }
  if (0) {
  } else {
    LOGPRINTF("--score: Results written to %s.\n", outname);
  }
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
