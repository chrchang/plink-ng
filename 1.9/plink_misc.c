// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink_common.h"

#include "plink_misc.h"
#include "plink_stats.h"

#include "pigz.h"

void misc_init(Score_info* sc_ip) {
  sc_ip->fname = nullptr;
  sc_ip->range_fname = nullptr;
  sc_ip->data_fname = nullptr;
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

int32_t make_founders(uintptr_t unfiltered_sample_ct, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t require_two, uintptr_t* sample_exclude, uintptr_t* founder_info) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uint32_t new_founder_ct = 0;
  int32_t retval = 0;
  char* sorted_ids;
  char* id_buf;
  char* wptr;
  char* pat_ptr;
  char* mat_ptr;
  uintptr_t* nf_bitarr;
  uintptr_t sample_uidx;
  uint32_t fam_len_p1;
  uint32_t missing_parent_ct;
  uint32_t cur_len;
  if (bigstack_alloc_c(max_sample_id_len, &id_buf) ||
      bigstack_alloc_ul(unfiltered_sample_ctl, &nf_bitarr)) {
    goto make_founders_ret_NOMEM;
  }
  bitarr_invert_copy(sample_exclude, unfiltered_sample_ct, nf_bitarr);
  bitvec_andnot(founder_info, unfiltered_sample_ctl, nf_bitarr);
  sample_uidx = unfiltered_sample_ct? next_set(nf_bitarr, 0, unfiltered_sample_ct) : 0;
  if (sample_uidx == unfiltered_sample_ct) {
    logprint("Note: Skipping --make-founders since there are no nonfounders.\n");
    goto make_founders_ret_1;
  }
  sorted_ids = alloc_and_init_collapsed_arr(sample_ids, max_sample_id_len, unfiltered_sample_ct, sample_exclude, sample_ct, 0);
  if (!sorted_ids) {
    goto make_founders_ret_NOMEM;
  }
  qsort(sorted_ids, sample_ct, max_sample_id_len, strcmp_casted);
  do {
    pat_ptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    fam_len_p1 = strlen_se(pat_ptr) + 1;
    wptr = memcpya(id_buf, pat_ptr, fam_len_p1);
    missing_parent_ct = 0;
    pat_ptr = &(paternal_ids[sample_uidx * max_paternal_id_len]);
    cur_len = strlen(pat_ptr);
    if (cur_len + fam_len_p1 >= max_sample_id_len) {
      missing_parent_ct++;
    } else {
      memcpy(wptr, pat_ptr, cur_len);
      if (bsearch_str(id_buf, cur_len + fam_len_p1, sorted_ids, max_sample_id_len, sample_ct) == -1) {
	missing_parent_ct++;
      }
    }
    mat_ptr = &(maternal_ids[sample_uidx * max_maternal_id_len]);
    cur_len = strlen(mat_ptr);
    if (cur_len + fam_len_p1 >= max_sample_id_len) {
      missing_parent_ct++;
    } else {
      memcpy(wptr, mat_ptr, cur_len);
      if (bsearch_str(id_buf, cur_len + fam_len_p1, sorted_ids, max_sample_id_len, sample_ct) == -1) {
	missing_parent_ct++;
      }
    }
    if (missing_parent_ct > require_two) {
      SET_BIT(sample_uidx, founder_info);
      memcpy(pat_ptr, "0", 2);
      memcpy(mat_ptr, "0", 2);
      new_founder_ct++;
    }
    sample_uidx++;
    next_set_ul_ck(nf_bitarr, unfiltered_sample_ct, &sample_uidx);
  } while (sample_uidx < unfiltered_sample_ct);
  LOGPRINTF("--make-founders: %u sample%s affected.\n", new_founder_ct, (new_founder_ct == 1)? "" : "s");
  while (0) {
  make_founders_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
 make_founders_ret_1:
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t write_nosex(char* outname, char* outname_end, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sex_nm, uintptr_t gender_unk_ct, char* sample_ids, uintptr_t max_sample_id_len) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_uidx = 0;
  int32_t retval = 0;
  uintptr_t* sex_missing;
  uintptr_t sample_idx;
  if (bigstack_alloc_ul(unfiltered_sample_ctl, &sex_missing)) {
    goto write_nosex_ret_NOMEM;
  }
  bitarr_invert_copy(sample_exclude, unfiltered_sample_ct, sex_missing);
  bitvec_andnot(sex_nm, unfiltered_sample_ctl, sex_missing);
  memcpy(outname_end, ".nosex", 7);
  if (fopen_checked(outname, "w", &outfile)) {
    goto write_nosex_ret_OPEN_FAIL;
  }
  for (sample_idx = 0; sample_idx < gender_unk_ct; sample_idx++, sample_uidx++) {
    next_set_ul_unsafe_ck(sex_missing, &sample_uidx);
    fputs(&(sample_ids[sample_uidx * max_sample_id_len]), outfile);
    putc_unlocked('\n', outfile);
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
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t makepheno_load(FILE* phenofile, char* makepheno_str, uintptr_t unfiltered_sample_ct, char* sorted_sample_ids, uintptr_t max_sample_id_len, uint32_t* id_map, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr) {
  uint32_t mp_strlen = strlen(makepheno_str);
  uint32_t makepheno_all = ((mp_strlen == 1) && (makepheno_str[0] == '*'));
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t* pheno_c = *pheno_c_ptr;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t line_idx = 0;
  int32_t retval = 0;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  int32_t ii;
  uint32_t sample_idx;
  uint32_t tmp_len;
  if (bigstack_alloc_c(max_sample_id_len, &id_buf)) {
    goto makepheno_load_ret_NOMEM;
  }
  if (!pheno_c) {
    if (aligned_malloc(unfiltered_sample_ctl * sizeof(intptr_t), pheno_c_ptr)) {
      goto makepheno_load_ret_NOMEM;
    }
    pheno_c = *pheno_c_ptr;
    fill_ulong_zero(unfiltered_sample_ctl, pheno_c);
  }
  if (makepheno_all) {
    fill_all_bits(unfiltered_sample_ct, pheno_nm);
  }
  // probably want to permit long lines here
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, phenofile) != nullptr) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --make-pheno file is pathologically long.\n", line_idx);
      goto makepheno_load_ret_INVALID_FORMAT_2;
    }
    bufptr0 = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr0)) {
      continue;
    }
    if (bsearch_read_fam_indiv(bufptr0, sorted_sample_ids, max_sample_id_len, unfiltered_sample_ct, &bufptr, &ii, id_buf)) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --make-pheno file has fewer tokens than expected.\n", line_idx);
      goto makepheno_load_ret_INVALID_FORMAT_2;
    }
    if (ii != -1) {
      sample_idx = id_map[(uint32_t)ii];
      if (makepheno_all) {
	SET_BIT(sample_idx, pheno_c);
      } else {
	SET_BIT(sample_idx, pheno_nm);
        tmp_len = strlen_se(bufptr);
	if ((tmp_len == mp_strlen) && (!memcmp(bufptr, makepheno_str, mp_strlen))) {
	  SET_BIT(sample_idx, pheno_c);
	}
      }
    }
  }
  if (!feof(phenofile)) {
    goto makepheno_load_ret_READ_FAIL;
  }
  tmp_len = popcount_longs(pheno_nm, unfiltered_sample_ctl);
  LOGPRINTF("--make-pheno: %u phenotype value%s set.\n", tmp_len, (tmp_len == 1)? "" : "s");
  while (0) {
  makepheno_load_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  makepheno_load_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  makepheno_load_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t load_pheno(FILE* phenofile, uintptr_t unfiltered_sample_ct, uintptr_t sample_exclude_ct, char* sorted_sample_ids, uintptr_t max_sample_id_len, uint32_t* id_map, int32_t missing_pheno, uint32_t affection_01, uint32_t mpheno_col, char* phenoname_str, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, char* phenoname_load, uintptr_t max_pheno_name_len) {
  uint32_t affection = 1;
  uintptr_t* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  int32_t header_processed = 0;
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ct = unfiltered_sample_ct - sample_exclude_ct;
  uintptr_t line_idx = 0;
  uintptr_t* isz = nullptr;
  double pheno_ctrld = (double)((int32_t)(1 - affection_01));
  double pheno_cased = pheno_ctrld + 1.0;
  double missing_phenod = (double)missing_pheno;
  int32_t retval = 0;
  double dxx;
  uintptr_t loadbuf_size;
  char* loadbuf;
  char* bufptr0;
  char* bufptr;
  char* ss;
  uint32_t tmp_len;
  uint32_t tmp_len2;
  uint32_t uii;
  int32_t sample_idx;
  if (pheno_d) {
    affection = 0;
  } else {
    if (bigstack_calloc_ul(unfiltered_sample_ctl, &isz)) {
      goto load_pheno_ret_NOMEM;
    }
    if (!pheno_c) {
      if (aligned_malloc(unfiltered_sample_ctl * sizeof(intptr_t), pheno_c_ptr)) {
	goto load_pheno_ret_NOMEM;
      }
      pheno_c = *pheno_c_ptr;
      fill_ulong_zero(unfiltered_sample_ctl, pheno_c);
    }
  }
  // ----- phenotype file load -----
  // worthwhile to support very long lines here
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_pheno_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, phenofile) != nullptr) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --pheno file is pathologically long.\n", line_idx);
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
      if ((((tmp_len == 3) && (!memcmp("FID", bufptr0, 3))) || ((tmp_len == 4) && (!memcmp("#FID", bufptr0, 4)))) && (tmp_len2 == 3) && (!memcmp("IID", bufptr, 3))) {
	if (phenoname_str) {
	  tmp_len = strlen(phenoname_str);
	  do {
	    bufptr = next_token(bufptr);
	    if (no_more_tokens_kns(bufptr)) {
	      logerrprint("Error: --pheno-name column not found.\n");
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
	    logerrprint("Error: Excessively long phenotype name in --pheno file.\n");
            goto load_pheno_ret_INVALID_FORMAT;
	  }
          memcpyx(phenoname_load, bufptr, tmp_len, '\0');
	}
      } else if (phenoname_str) {
	logerrprint("Error: --pheno-name requires the --pheno file to have a header line with first\ntwo columns 'FID' and 'IID'.\n");
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
      if (bsearch_read_fam_indiv(bufptr0, sorted_sample_ids, max_sample_id_len, sample_ct, &bufptr, &sample_idx, g_textbuf)) {
	goto load_pheno_ret_MISSING_TOKENS;
      }
      if (sample_idx != -1) {
	sample_idx = id_map[(uint32_t)sample_idx];
	if (mpheno_col > 1) {
	  bufptr = next_token_mult(bufptr, mpheno_col - 1);
	}
	if (no_more_tokens_kns(bufptr)) {
	  // Sometimes, but not always, an error.  So we populate g_logbuf but
	  // let the caller decide whether to actually log it.
          sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --pheno file has fewer tokens than expected.\n", line_idx);
	  return LOAD_PHENO_LAST_COL;
	}
	dxx = strtod(bufptr, &ss);

	// if conversion fails, strtod return value is 0.0.  that has some
	// unpleasant interactions with the following code, and my previous
	// attempt to account for them had a bug.
	if (ss == bufptr) {
	  dxx = missing_phenod;
	}

	if (affection) {
	  if (dxx == pheno_cased) {
	    set_bit(sample_idx, pheno_c);
	    set_bit(sample_idx, pheno_nm);
	  } else if (dxx == pheno_ctrld) {
	    clear_bit(sample_idx, pheno_c);
	    set_bit(sample_idx, pheno_nm);
	  } else if (affection_01 || (dxx == missing_phenod) || (dxx == 0.0)) {
	    // Since we're only making one pass through the file, we don't
	    // have the luxury of knowing in advance whether the phenotype is
	    // binary or scalar.  If there is a '0' entry that occurs before
	    // we know the phenotype is scalar, we need to not set the
	    // phenotype to zero during the binary -> scalar conversion step.
	    if (dxx == 0.0) {
	      set_bit(sample_idx, isz);
	    }
	    clear_bit(sample_idx, pheno_nm);
	    clear_bit(sample_idx, pheno_c);
	  } else {
	    pheno_d = (double*)malloc(unfiltered_sample_ct * sizeof(double));
	    if (!pheno_d) {
	      goto load_pheno_ret_NOMEM;
	    }
	    *pheno_d_ptr = pheno_d;
	    for (uii = 0; uii < unfiltered_sample_ct; uii++) {
	      if (is_set(isz, uii)) {
		pheno_d[uii] = 0.0;
		set_bit(uii, pheno_nm);
	      } else if (is_set(pheno_nm, uii)) {
		pheno_d[uii] = (double)((int32_t)(1 + is_set(pheno_c, uii)));
	      }
	    }
	    aligned_free_null(pheno_c_ptr);
	    affection = 0;
	  }
	}
	if (!affection) {
	  if (dxx != missing_phenod) {
	    pheno_d[(uint32_t)sample_idx] = dxx;
	    set_bit(sample_idx, pheno_nm);
	  }
	}
      }
    }
  }
  if (!feof(phenofile)) {
    goto load_pheno_ret_READ_FAIL;
  }
  uii = popcount_longs(pheno_nm, unfiltered_sample_ctl);
  LOGPRINTF("%u phenotype value%s present after --pheno.\n", uii, (uii == 1)? "" : "s");
  while (0) {
  load_pheno_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_pheno_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_pheno_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --pheno file has fewer tokens than expected.\n", line_idx);
  load_pheno_ret_INVALID_FORMAT_2:
    logerrprintb();
  load_pheno_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t convert_tail_pheno(uint32_t unfiltered_sample_ct, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, double tail_bottom, double tail_top, double missing_phenod) {
  uintptr_t* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  uint32_t sample_uidx;
  uint32_t sample_uidx_stop;
  double dxx;
  if (!(*pheno_d_ptr)) {
    logerrprint("Error: --tail-pheno requires scalar phenotype data.\n");
    return RET_INVALID_FORMAT;
  }
  sample_uidx = BITCT_TO_WORDCT(unfiltered_sample_ct);
  if (!pheno_c) {
    if (aligned_malloc(sample_uidx * sizeof(intptr_t), pheno_c_ptr)) {
      return RET_NOMEM;
    }
    pheno_c = *pheno_c_ptr;
  }
  fill_ulong_zero(sample_uidx, pheno_c);
  sample_uidx = 0;
  do {
    sample_uidx = next_set(pheno_nm, sample_uidx, unfiltered_sample_ct);
    sample_uidx_stop = next_unset(pheno_nm, sample_uidx, unfiltered_sample_ct);
    for (; sample_uidx < sample_uidx_stop; sample_uidx++) {
      dxx = pheno_d[sample_uidx];
      if (dxx > tail_bottom) {
        if (dxx > tail_top) {
          SET_BIT(sample_uidx, pheno_c);
        } else {
	  CLEAR_BIT(sample_uidx, pheno_nm);
        }
      }
    }
  } while (sample_uidx_stop < unfiltered_sample_ct);
  free(pheno_d);
  *pheno_d_ptr = nullptr;
  sample_uidx = popcount_longs(pheno_nm, BITCT_TO_WORDCT(unfiltered_sample_ct));
  LOGPRINTF("--tail-pheno: %u phenotype value%s remaining.\n", sample_uidx, (sample_uidx == 1)? "" : "s");
  return 0;
}

int32_t apply_cm_map(char* cm_map_fname, char* cm_map_chrname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, double* marker_cms, Chrom_info* chrom_info_ptr) {
  FILE* shapeitfile = nullptr;
  char* at_sign_ptr = nullptr;
  char* fname_write = nullptr;
  char* fname_buf = &(g_textbuf[MAXLINELEN]);
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
  {
    if (!cm_map_chrname) {
      chrom_fo_idx = 0;
      chrom_ct = chrom_info_ptr->chrom_ct;
      at_sign_ptr = strchr(cm_map_fname, '@');
      fname_write = memcpya(fname_buf, cm_map_fname, (uintptr_t)(at_sign_ptr - cm_map_fname));
      at_sign_ptr++;
      post_at_sign_len = strlen(at_sign_ptr) + 1;
    } else {
      const uint32_t chrom_name_slen = strlen(cm_map_chrname);
      int32_t cur_chrom_code = get_chrom_code(cm_map_chrname, chrom_info_ptr, chrom_name_slen);
      if (cur_chrom_code < 0) {
	LOGPREPRINTFWW("Error: --cm-map chromosome code '%s' not found in dataset.\n", cm_map_chrname);
	goto apply_cm_map_ret_INVALID_CMDLINE_2;
      }
      chrom_fo_idx = chrom_info_ptr->chrom_idx_to_foidx[(uint32_t)cur_chrom_code];
      chrom_ct = chrom_fo_idx + 1;
      fname_buf = cm_map_fname;
    }
    g_textbuf[MAXLINELEN - 1] = ' ';
    for (; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
      chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
      marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
      if (marker_uidx == chrom_end) {
	continue;
      }
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      if (!cm_map_chrname) {
	if ((!uii) || (uii > autosome_ct)) {
	  continue;
	}
	bufptr = uint32toa(uii, fname_write);
	memcpy(bufptr, at_sign_ptr, post_at_sign_len);
	if (fopen_checked(fname_buf, "r", &shapeitfile)) {
	  LOGERRPRINTFWW("Warning: --cm-map failed to open %s.\n", fname_buf);
	  continue;
	}
      } else {
	if (fopen_checked(cm_map_fname, "r", &shapeitfile)) {
	  goto apply_cm_map_ret_OPEN_FAIL;
	}
      }
      updated_chrom_ct++;
      irreg_line_ct = 0;
      // First line is a header with three arbitrary fields.
      // All subsequent lines have three fields in the following order:
      //   1. bp position (increasing)
      //   2. cM/Mb recombination rate between current and previous bp
      //      positions
      //   3. current cM position
      // We mostly ignore field 2, since depending just on fields 1 and 3
      // maximizes accuracy.  The one exception is the very first nonheader
      // line.
      retval = load_to_first_token(shapeitfile, MAXLINELEN, '\0', "--cm-map file", g_textbuf, &bufptr, &line_idx);
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
      while (fgets(g_textbuf, MAXLINELEN, shapeitfile)) {
	line_idx++;
	if (!g_textbuf[MAXLINELEN - 1]) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --cm-map file is pathologically long.\n", line_idx);
	  goto apply_cm_map_ret_INVALID_FORMAT_2;
	}
	bufptr = skip_initial_spaces(g_textbuf);
	if ((*bufptr < '+') || (*bufptr > '9')) {
	  // warning instead of error if text line found, since as of 8 Jan 2014
	  // the posted chromosome 19 map has such a line
	  if (*bufptr > ' ') {
	    irreg_line_ct++;
	  }
	  continue;
	}
	if (scan_uint_defcap(bufptr, (uint32_t*)&bp_new)) {
	  sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of --cm-map file.\n", line_idx);
	  goto apply_cm_map_ret_INVALID_FORMAT_2;
	}
	if (bp_new <= bp_old) {
	  sprintf(g_logbuf, "Error: bp coordinates in --cm-map file are not in increasing order ('%d' on line %" PRIuPTR " is not larger than previous value '%d').\n", bp_new, line_idx, bp_old);
          wordwrapb(0);
	  goto apply_cm_map_ret_INVALID_FORMAT_2;
	}
	bufptr2 = next_token_mult(bufptr, 2);
	if (no_more_tokens_kns(bufptr2)) {
	  goto apply_cm_map_ret_MISSING_TOKENS;
	}
	if (scan_double(bufptr2, &cm_new)) {
	  sprintf(g_logbuf, "Error: Invalid centimorgan position on line %" PRIuPTR " of --cm-map file.\n", line_idx);
	  goto apply_cm_map_ret_INVALID_FORMAT_2;
	}
	if (bp_old == -1) {
	  // parse field 2 only in this case
	  bufptr = next_token(bufptr);
	  if (scan_double(bufptr, &dxx)) {
	    sprintf(g_logbuf, "Error: Invalid recombination rate on line %" PRIuPTR " of --cm-map file.\n", line_idx);
	    goto apply_cm_map_ret_INVALID_FORMAT_2;
	  }
	  cm_old = cm_new - dxx * 0.000001 * ((double)(bp_new + 1));
	}
	dxx = (cm_new - cm_old) / ((double)(bp_new - bp_old));
	while (marker_pos[marker_uidx] <= ((uint32_t)bp_new)) {
	  marker_cms[marker_uidx] = cm_new - ((int32_t)(((uint32_t)bp_new) - marker_pos[marker_uidx])) * dxx;
	  marker_uidx++;
	  next_unset_ck(marker_exclude, chrom_end, &marker_uidx);
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
	LOGERRPRINTFWW("Warning: %" PRIuPTR " irregular line%s skipped in %s.\n", irreg_line_ct, (irreg_line_ct == 1)? "" : "s", fname_buf);
      }
    }
    LOGPRINTF("--cm-map: %u chromosome%s updated.\n", updated_chrom_ct, (updated_chrom_ct == 1)? "" : "s");
  }
  while (0) {
  apply_cm_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  apply_cm_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  apply_cm_map_ret_INVALID_CMDLINE_2:
    logerrprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  apply_cm_map_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --cm-map file has fewer tokens than expected.\n", line_idx);
  apply_cm_map_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 apply_cm_map_ret_1:
  fclose_cond(shapeitfile);
  return retval;
}

int32_t update_marker_cms(Two_col_params* update_cm, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, double* marker_cms) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  char skipchar = update_cm->skipchar;
  uint32_t colid_first = (update_cm->colid < update_cm->colx);
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
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
  uint32_t marker_uidx;
  char cc;
  int32_t retval;
  if (bigstack_calloc_ul(unfiltered_marker_ctl, &already_seen)) {
    goto update_marker_cms_ret_NOMEM;
  }

  loadbuf = (char*)g_bigstack_base;
  loadbuf_size = bigstack_left();
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
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-cm file is pathologically long.\n", line_idx);
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
    marker_uidx = id_htable_find(colid_ptr, slen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if (marker_uidx == 0xffffffffU) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, marker_uidx)) {
      colid_ptr[slen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant '%s' in --update-cm file.\n", colid_ptr);
      goto update_marker_cms_ret_INVALID_FORMAT_2;
    }
    set_bit(marker_uidx, already_seen);
    if (scan_double(colx_ptr, &(marker_cms[marker_uidx]))) {
      sprintf(g_logbuf, "Error: Invalid centimorgan position on line %" PRIuPTR " of --update-cm file.\n", line_idx);
      goto update_marker_cms_ret_INVALID_FORMAT_2;
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_cms_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(g_logbuf, "--update-cm: %" PRIuPTR " value%s changed, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(g_logbuf, "--update-cm: %" PRIuPTR " value%s changed.\n", hit_ct, (hit_ct == 1)? "" : "s");
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
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-cm file has fewer tokens than expected.\n", line_idx);
  update_marker_cms_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 update_marker_cms_ret_1:
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t update_marker_pos(Two_col_params* update_map, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t* marker_pos, uint32_t* map_is_unsorted_ptr, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  char skipchar = update_map->skipchar;
  uint32_t colid_first = (update_map->colid < update_map->colx);
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t map_is_unsorted = ((*map_is_unsorted_ptr) & UNSORTED_CHROM);
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
  uint32_t marker_uidx;
  uint32_t marker_idx;
  int32_t bp_coord;
  int32_t retval;
  char cc;
  if (bigstack_calloc_ul(unfiltered_marker_ctl, &already_seen)) {
    goto update_marker_pos_ret_NOMEM;
  }

  loadbuf = (char*)g_bigstack_base;
  loadbuf_size = bigstack_left();
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
        sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-map file is pathologically long.\n", line_idx);
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
    marker_uidx = id_htable_find(colid_ptr, slen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if (marker_uidx == 0xffffffffU) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, marker_uidx)) {
      colid_ptr[slen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant '%s' in --update-map file.\n", colid_ptr);
      goto update_marker_pos_ret_INVALID_FORMAT_2;
    }
    set_bit(marker_uidx, already_seen);
    if (scan_int_abs_defcap(colx_ptr, &bp_coord)) {
      sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of --update-map file.\n", line_idx);
      goto update_marker_pos_ret_INVALID_FORMAT_2;
    }
    if (bp_coord < 0) {
      set_bit(marker_uidx, marker_exclude);
      marker_ct--;
    } else {
      marker_pos[marker_uidx] = bp_coord;
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_pos_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(g_logbuf, "--update-map: %" PRIuPTR " value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(g_logbuf, "--update-map: %" PRIuPTR " value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  *marker_exclude_ct_ptr = unfiltered_marker_ct - marker_ct;
  if (!marker_ct) {
    logerrprint("Error: All variants excluded by --update-map (due to negative marker\npositions).\n");
    goto update_marker_pos_ret_ALL_MARKERS_EXCLUDED;
  }
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    while (marker_uidx >= chrom_end) {
      chrom_end = chrom_info_ptr->chrom_fo_vidx_start[++chrom_fo_idx_p1];
      last_pos = 0;
    }
    if (last_pos > marker_pos[marker_uidx]) {
      map_is_unsorted |= UNSORTED_BP;
      if (!((*map_is_unsorted_ptr) & UNSORTED_BP)) {
	logerrprint("Warning: Base-pair positions are now unsorted!\n");
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
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-map file has fewer tokens than expected.\n", line_idx);
  update_marker_pos_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  update_marker_pos_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 update_marker_pos_ret_1:
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t update_marker_names(Two_col_params* update_name, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  char skipchar = update_name->skipchar;
  uint32_t colold_first = (update_name->colid < update_name->colx);
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t* already_seen;
  char* marker_ids_copy;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uint32_t colmin;
  uint32_t coldiff;
  char* colold_ptr;
  char* colnew_ptr;
  uintptr_t marker_uidx;
  uint32_t slen;
  int32_t retval;
  char cc;
  if (bigstack_calloc_ul(unfiltered_marker_ctl, &already_seen) ||
      bigstack_alloc_c(unfiltered_marker_ct * max_marker_id_len, &marker_ids_copy)) {
    goto update_marker_names_ret_NOMEM;
  }
  memcpy(marker_ids_copy, marker_ids, unfiltered_marker_ct * max_marker_id_len);
  loadbuf = (char*)g_bigstack_base;
  loadbuf_size = bigstack_left();
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
    marker_uidx = id_htable_find(colold_ptr, slen, marker_id_htable, marker_id_htable_size, marker_ids_copy, max_marker_id_len);
    if (marker_uidx == 0xffffffffU) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, marker_uidx)) {
      colold_ptr[slen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --update-name file.\n", colold_ptr);
      goto update_marker_names_ret_INVALID_FORMAT_2;
    }
    set_bit(marker_uidx, already_seen);
    slen = strlen_se(colnew_ptr);
    colnew_ptr[slen] = '\0';
    memcpy(&(marker_ids[marker_uidx * max_marker_id_len]), colnew_ptr, slen + 1);
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_names_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(g_logbuf, "--update-name: %" PRIuPTR " value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(g_logbuf, "--update-name: %" PRIuPTR " value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
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
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 update_marker_names_ret_1:
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t update_marker_alleles(char* update_alleles_fname, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  FILE* errfile = nullptr;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t max_marker_allele_len = *max_marker_allele_len_ptr;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t err_ct = 0;
  uintptr_t line_idx = 0;
  int32_t retval = 0;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  char missing_geno = *missing_geno_ptr;
  uintptr_t* already_seen;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  uintptr_t loadbuf_size;
  uintptr_t marker_uidx;
  uint32_t len2;
  uint32_t len;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t uoo;
  uint32_t upp;
  uint32_t uqq;
  if (bigstack_calloc_ul(unfiltered_marker_ctl, &already_seen)) {
    goto update_marker_alleles_ret_NOMEM;
  }
  if (fopen_checked(update_alleles_fname, "r", &infile)) {
    goto update_marker_alleles_ret_OPEN_FAIL;
  }
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_alleles_ret_NOMEM;
  }
  loadbuf = (char*)bigstack_alloc(loadbuf_size);
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-alleles file is pathologically long.\n", line_idx);
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
    marker_uidx = id_htable_find(bufptr3, (uintptr_t)(bufptr2 - bufptr3), marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if ((marker_uidx == 0xffffffffU) || IS_SET(marker_exclude, marker_uidx)) {
      miss_ct++;
      continue;
    }
    if (IS_SET(already_seen, marker_uidx)) {
      *bufptr2 = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --update-alleles file.\n", bufptr3);
      goto update_marker_alleles_ret_INVALID_FORMAT_2;
    }
    SET_BIT(marker_uidx, already_seen);
    bufptr2 = skip_initial_spaces(bufptr2);
    len2 = strlen_se(bufptr2);
    bufptr = &(bufptr2[len2]);
    *bufptr = '\0';
    bufptr = skip_initial_spaces(&(bufptr[1]));
    len = strlen_se(bufptr);
    if ((len == len2) && (!memcmp(bufptr, bufptr2, len))) {
      goto update_marker_alleles_ret_DUPLICATE_ALLELE_CODE;
    }
    bufptr[len] = '\0';
    uii = (len2 == 1) && (*bufptr2 == missing_geno);
    ujj = (len == 1) && (*bufptr == missing_geno);
    ukk = (marker_allele_ptrs[2 * marker_uidx] == missing_geno_ptr);
    umm = (marker_allele_ptrs[2 * marker_uidx + 1] == missing_geno_ptr);
    unn = !strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx]);
    uoo = !strcmp(bufptr, marker_allele_ptrs[2 * marker_uidx + 1]);
    upp = !strcmp(bufptr, marker_allele_ptrs[2 * marker_uidx]);
    uqq = !strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx + 1]);
    // allow missing allele codes as long as they don't introduce ambiguity.
    // i.e. at least one exact match.
    if ((unn && (ujj || umm || uoo)) || (uoo && (uii || ukk))) {
      bufptr2 = skip_initial_spaces(&(bufptr[len + 1]));
      bufptr = next_token(bufptr2);
      goto update_marker_alleles_match;
    } else if ((upp && (uii || umm || uqq)) || (uqq && (ujj || ukk))) {
      bufptr = skip_initial_spaces(&(bufptr[len + 1]));
      bufptr2 = next_token(bufptr);
    update_marker_alleles_match:
      len = strlen_se(bufptr);
      len2 = strlen_se(bufptr2);
      if ((len == len2) && (!memcmp(bufptr, bufptr2, len))) {
	goto update_marker_alleles_ret_DUPLICATE_ALLELE_CODE;
      }
      if (memchr(bufptr, ',', len) || memchr(bufptr2, ',', len2)) {
	// this breaks VCF and PLINK 2 binary
        LOGPREPRINTFWW("Error: Comma-containing new allele code on line %" PRIuPTR " of --update-alleles file.\n", line_idx);
	goto update_marker_alleles_ret_INVALID_FORMAT_2;
      }
      if (len >= max_marker_allele_len) {
	max_marker_allele_len = len + 1;
      }
      if (len2 >= max_marker_allele_len) {
	max_marker_allele_len = len2 + 1;
      }
      allele_reset(bufptr2, len2, &(marker_allele_ptrs[2 * marker_uidx]));
      allele_reset(bufptr, len, &(marker_allele_ptrs[2 * marker_uidx + 1]));
      hit_ct++;
    } else {
      if (!err_ct) {
	memcpy(outname_end, ".allele.no.snp", 15);
	if (fopen_checked(outname, "w", &errfile)) {
	  goto update_marker_alleles_ret_OPEN_FAIL;
	}
      }
      *token_endnn(bufptr3) = '\0';
      fputs(bufptr3, errfile);
      putc_unlocked('\t', errfile);
      fputs(bufptr2, errfile);
      putc_unlocked('\t', errfile);
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
    sprintf(g_logbuf, "--update-alleles: %" PRIuPTR " variant%s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(g_logbuf, "--update-alleles: %" PRIuPTR " variant%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
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
  update_marker_alleles_ret_DUPLICATE_ALLELE_CODE:
    LOGERRPRINTF("Error: Duplicate allele code on line %" PRIuPTR " of --update-alleles file.\n", line_idx);
  update_marker_alleles_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  fclose_cond(errfile);
  bigstack_reset(bigstack_mark);
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

uint32_t flip_process_token(char* tok_start, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_exclude, uintptr_t* already_seen, char** marker_allele_ptrs, uintptr_t* hit_ct_ptr, uintptr_t* miss_ct_ptr, uint32_t* non_acgt_ct_ptr) {
  uint32_t slen = strlen_se(tok_start);
  uintptr_t marker_uidx;
  uint32_t cur_non_acgt0;
  marker_uidx = id_htable_find(tok_start, slen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
  if ((marker_uidx == 0xffffffffU) || IS_SET(marker_exclude, marker_uidx)) {
    *miss_ct_ptr += 1;
    return 0;
  }
  if (IS_SET(already_seen, marker_uidx)) {
    tok_start[slen] = '\0';
    LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --flip file.\n", tok_start);
    return 1;
  }
  SET_BIT(marker_uidx, already_seen);
  cur_non_acgt0 = 0;
  cur_non_acgt0 |= flip_str(&(marker_allele_ptrs[2 * marker_uidx]));
  cur_non_acgt0 |= flip_str(&(marker_allele_ptrs[2 * marker_uidx + 1]));
  *non_acgt_ct_ptr += (cur_non_acgt0 & 1);
  *hit_ct_ptr += (cur_non_acgt0 >> 1);
  return 0;
}

int32_t flip_strand(char* flip_fname, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char** marker_allele_ptrs) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* flipfile = nullptr;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  char* midbuf = &(g_textbuf[MAXLINELEN]);
  uint32_t non_acgt_ct = 0;
  uint32_t curtoklen = 0;
  int32_t retval = 0;
  uintptr_t bufsize;
  uintptr_t* already_seen;
  char* bufptr0;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  if (bigstack_calloc_ul(unfiltered_marker_ctl, &already_seen)) {
    goto flip_strand_ret_NOMEM;
  }
  // Compatibility fix: PLINK 1.07 uses a token- rather than a line-based
  // loader here.
  if (fopen_checked(flip_fname, FOPEN_RB, &flipfile)) {
    goto flip_strand_ret_OPEN_FAIL;
  }
  while (1) {
    if (fread_checked(midbuf, MAXLINELEN, flipfile, &bufsize)) {
      goto flip_strand_ret_READ_FAIL;
    }
    if (!bufsize) {
      if (curtoklen) {
	g_textbuf[MAXLINELEN] = ' '; // bugfix for no-\n-on-last-line case
        if (flip_process_token(&(g_textbuf[MAXLINELEN - curtoklen]), marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len, marker_exclude, already_seen, marker_allele_ptrs, &hit_ct, &miss_ct, &non_acgt_ct)) {
	  goto flip_strand_ret_INVALID_FORMAT_2;
	}
      }
      break;
    }
    bufptr0 = &(midbuf[bufsize]);
    *bufptr0 = ' ';
    bufptr0[1] = '0';
    bufptr = &(g_textbuf[MAXLINELEN - curtoklen]);
    bufptr2 = midbuf;
    if (curtoklen) {
      goto flip_strand_tok_start;
    }
    while (1) {
      while (*bufptr <= ' ') {
	bufptr++;
      }
      if (bufptr >= bufptr0) {
        curtoklen = 0;
        break;
      }
      bufptr2 = &(bufptr[1]);
    flip_strand_tok_start:
      while (*bufptr2 > ' ') {
        bufptr2++;
      }
      curtoklen = (uintptr_t)(bufptr2 - bufptr);
      if (bufptr2 == &(g_textbuf[MAXLINELEN * 2])) {
        if (curtoklen > MAX_ID_SLEN) {
	  logerrprint("Error: Excessively long ID in --flip file.\n");
	  goto flip_strand_ret_INVALID_FORMAT;
	}
        bufptr3 = &(g_textbuf[MAXLINELEN - curtoklen]);
        memcpy(bufptr3, bufptr, curtoklen);
	break;
      }
      if (flip_process_token(bufptr, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len, marker_exclude, already_seen, marker_allele_ptrs, &hit_ct, &miss_ct, &non_acgt_ct)) {
	goto flip_strand_ret_INVALID_FORMAT_2;
      }
      bufptr = &(bufptr2[1]);
    }
  }
  if (!feof(flipfile)) {
    goto flip_strand_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(g_logbuf, "--flip: %" PRIuPTR " SNP%s flipped, %" PRIuPTR " SNP ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(g_logbuf, "--flip: %" PRIuPTR " SNP%s flipped.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  if (non_acgt_ct) {
    LOGERRPRINTF("Warning: %u variant%s had at least one non-A/C/G/T allele name.\n", non_acgt_ct, (non_acgt_ct == 1)? "" : "s");
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
    logerrprintb();
  flip_strand_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(flipfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t update_sample_ids(char* update_ids_fname, char* sorted_sample_ids, uintptr_t sample_ct, uintptr_t max_sample_id_len, uint32_t* sample_id_map, char* sample_ids) {
  // file has been pre-scanned
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  int32_t retval = 0;
  uintptr_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t line_idx = 0;
  char* idbuf;
  uintptr_t* already_seen;
  char* bufptr;
  char* bufptr2;
  char* wptr;
  uintptr_t sample_uidx;
  uint32_t len;
  int32_t sorted_idx;
  if (bigstack_alloc_c(max_sample_id_len, &idbuf) ||
      bigstack_calloc_ul(sample_ctl, &already_seen)) {
    goto update_sample_ids_ret_NOMEM;
  }
  if (fopen_checked(update_ids_fname, "r", &infile)) {
    goto update_sample_ids_ret_OPEN_FAIL;
  }
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      // er, either this buffer should be extended, or the
      // scan_max_fam_indiv_strlen() should use this length...
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-ids file is pathologically long.\n", line_idx);
      goto update_sample_ids_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bsearch_read_fam_indiv(bufptr, sorted_sample_ids, max_sample_id_len, sample_ct, &bufptr, &sorted_idx, idbuf);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      *strchr(idbuf, '\t') = ' ';
      LOGPREPRINTFWW("Error: Duplicate sample ID '%s' in --update-ids file.\n", idbuf);
      goto update_sample_ids_ret_INVALID_FORMAT_2;
    }
    set_bit(sorted_idx, already_seen);
    sample_uidx = sample_id_map[((uint32_t)sorted_idx)];
    wptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    len = strlen_se(bufptr);
    bufptr2 = &(bufptr[len]);
    wptr = memcpyax(wptr, bufptr, len, '\t');
    bufptr = skip_initial_spaces(&(bufptr2[1]));
    len = strlen_se(bufptr);
    if ((len == 1) && (*bufptr == '0')) {
      sprintf(g_logbuf, "Error: Invalid IID '0' on line %" PRIuPTR " of --update-ids file.\n", line_idx);
      goto update_sample_ids_ret_INVALID_FORMAT_2;
    }
    memcpyx(wptr, bufptr, len, '\0');
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_sample_ids_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(g_logbuf, "--update-ids: %" PRIuPTR " %s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, species_str(hit_ct), miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(g_logbuf, "--update-ids: %" PRIuPTR " %s updated.\n", hit_ct, species_str(hit_ct));
  }
  logprintb();

  while (0) {
  update_sample_ids_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_sample_ids_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_sample_ids_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_sample_ids_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t update_sample_parents(char* update_parents_fname, char* sorted_sample_ids, uintptr_t sample_ct, uintptr_t max_sample_id_len, uint32_t* sample_id_map, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  int32_t retval = 0;
  uintptr_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
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
  uintptr_t sample_uidx;
  uint32_t len;
  uint32_t len2;
  int32_t sorted_idx;
  if (bigstack_alloc_c(max_sample_id_len, &idbuf) ||
      bigstack_calloc_ul(sample_ctl, &already_seen)) {
    goto update_sample_parents_ret_NOMEM;
  }
  if (fopen_checked(update_parents_fname, "r", &infile)) {
    goto update_sample_parents_ret_OPEN_FAIL;
  }
  // permit very long lines since this can be pointed at .ped files
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto update_sample_parents_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, infile)) {
    // no line_idx since all the validation happened earlier
    if (!loadbuf[loadbuf_size - 1]) {
      goto update_sample_parents_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bsearch_read_fam_indiv(bufptr, sorted_sample_ids, max_sample_id_len, sample_ct, &bufptr, &sorted_idx, idbuf);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      *strchr(idbuf, '\t') = ' ';
      LOGPREPRINTFWW("Error: Duplicate sample ID '%s' in --update-parents file.\n", idbuf);
      goto update_sample_parents_ret_INVALID_FORMAT_2;
    }
    set_bit(sorted_idx, already_seen);
    sample_uidx = sample_id_map[((uint32_t)sorted_idx)];
    wptr = &(paternal_ids[sample_uidx * max_paternal_id_len]);
    len = strlen_se(bufptr);
    bufptr2 = &(bufptr[len]);
    memcpyx(wptr, bufptr, len, '\0');
    wptr = &(maternal_ids[sample_uidx * max_maternal_id_len]);
    bufptr3 = skip_initial_spaces(&(bufptr2[1]));
    len2 = strlen_se(bufptr3);
    memcpyx(wptr, bufptr3, len2, '\0');
    if ((len == 1) && (*bufptr == '0') && (len2 == 1) && (*bufptr3 == '0')) {
      SET_BIT(sample_uidx, founder_info);
    } else {
      CLEAR_BIT(sample_uidx, founder_info);
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_sample_parents_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(g_logbuf, "--update-parents: %" PRIuPTR " %s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, species_str(hit_ct), miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(g_logbuf, "--update-parents: %" PRIuPTR " %s updated.\n", hit_ct, species_str(hit_ct));
  }
  logprintb();

  while (0) {
  update_sample_parents_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_sample_parents_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_sample_parents_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_sample_parents_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t update_sample_sexes(char* update_sex_fname, uint32_t update_sex_col, char* sorted_sample_ids, uintptr_t sample_ct, uintptr_t max_sample_id_len, uint32_t* sample_id_map, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  int32_t retval = 0;
  uintptr_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t line_idx = 0;
  char* idbuf;
  uintptr_t* already_seen;
  char* loadbuf;
  char* bufptr;
  uintptr_t loadbuf_size;
  int32_t sorted_idx;
  uint32_t sample_uidx;
  char cc;
  unsigned char ucc;
  update_sex_col--;
  if (bigstack_alloc_c(max_sample_id_len, &idbuf) ||
      bigstack_calloc_ul(sample_ctl, &already_seen)) {
    goto update_sample_sexes_ret_NOMEM;
  }
  if (fopen_checked(update_sex_fname, "r", &infile)) {
    goto update_sample_sexes_ret_OPEN_FAIL;
  }
  // permit very long lines since this can be pointed at .ped files
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto update_sample_sexes_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-sex file is pathologically long.\n", line_idx);
	goto update_sample_sexes_ret_INVALID_FORMAT_2;
      } else {
	goto update_sample_sexes_ret_NOMEM;
      }
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(bufptr, sorted_sample_ids, max_sample_id_len, sample_ct, &bufptr, &sorted_idx, idbuf)) {
      goto update_sample_sexes_ret_MISSING_TOKENS;
    }
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      *strchr(idbuf, '\t') = ' ';
      LOGPREPRINTFWW("Error: Duplicate sample ID '%s' in --update-sex file.\n", idbuf);
      goto update_sample_sexes_ret_INVALID_FORMAT_2;
    }
    set_bit(sorted_idx, already_seen);
    sample_uidx = sample_id_map[((uint32_t)sorted_idx)];
    bufptr = next_token_multz(bufptr, update_sex_col);
    if (no_more_tokens_kns(bufptr)) {
      goto update_sample_sexes_ret_MISSING_TOKENS;
    }
    cc = *bufptr;
    ucc = ((unsigned char)cc) & 0xdfU;
    if ((cc < '0') || ((cc > '2') && (ucc != 'M') && (ucc != 'F')) || (bufptr[1] > ' ')) {
      sprintf(g_logbuf, "Error: Invalid sex value on line %" PRIuPTR " of --update-sex file.\n(Acceptable values: 1/M = male, 2/F = female, 0 = missing.)\n", line_idx);
      goto update_sample_sexes_ret_INVALID_FORMAT_2;
    }
    if (cc == '0') {
      CLEAR_BIT(sample_uidx, sex_nm);
      CLEAR_BIT(sample_uidx, sex_male);
    } else {
      SET_BIT(sample_uidx, sex_nm);
      if ((cc == '1') || (ucc == 'M')) {
	SET_BIT(sample_uidx, sex_male);
      } else {
	CLEAR_BIT(sample_uidx, sex_male);
      }
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_sample_sexes_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(g_logbuf, "--update-sex: %" PRIuPTR " %s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, species_str(hit_ct), miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(g_logbuf, "--update-sex: %" PRIuPTR " %s updated.\n", hit_ct, species_str(hit_ct));
  }
  logprintb();

  while (0) {
  update_sample_sexes_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_sample_sexes_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_sample_sexes_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_sample_sexes_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-sex file has fewer tokens than expected.\n", line_idx);
  update_sample_sexes_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

void calc_plink_maxfid(uint32_t unfiltered_sample_ct, uintptr_t* sample_exclude, uint32_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uint32_t* plink_maxfid_ptr, uint32_t* plink_maxiid_ptr) {
  uintptr_t plink_maxfid = 4;
  uintptr_t plink_maxiid = 4;
  uint32_t sample_uidx = 0;
  uint32_t samples_done = 0;
  char* cptr;
  char* cptr2;
  char* cptr_end;
  uintptr_t slen;
  uint32_t sample_uidx_stop;
  // imitate PLINK 1.07 behavior (see Plink::prettyPrintLengths() in
  // helper.cpp), to simplify testing and avoid randomly breaking existing
  // scripts
  while (samples_done < sample_ct) {
    sample_uidx = next_unset_unsafe(sample_exclude, sample_uidx);
    sample_uidx_stop = next_set(sample_exclude, sample_uidx, unfiltered_sample_ct);
    samples_done += sample_uidx_stop - sample_uidx;
    cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    cptr_end = &(sample_ids[sample_uidx_stop * max_sample_id_len]);
    sample_uidx = sample_uidx_stop;
    do {
      cptr2 = (char*)memchr(cptr, '\t', max_sample_id_len);
      slen = (uintptr_t)(cptr2 - cptr);
      if (slen > plink_maxfid) {
	plink_maxfid = slen + 2;
      }
      slen = strlen(&(cptr2[1]));
      if (slen > plink_maxiid) {
        plink_maxiid = slen + 2;
      }
      cptr = &(cptr[max_sample_id_len]);
    } while (cptr < cptr_end);
  }
  *plink_maxfid_ptr = plink_maxfid;
  *plink_maxiid_ptr = plink_maxiid;
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

// aptr1 = minor, aptr2 = major
int32_t load_one_freq(uint32_t alen1, const char* aptr1, uint32_t alen2, const char* aptr2, double maf, double* set_allele_freq_ptr, char** mastrs_ptr, char missing_geno) {
  uint32_t malen0 = strlen(mastrs_ptr[0]);
  uint32_t malen1 = strlen(mastrs_ptr[1]);
  uint32_t missing0 = (mastrs_ptr[0][0] == missing_geno) && (malen0 == 1);
  uint32_t missing1 = (mastrs_ptr[1][0] == missing_geno) && (malen1 == 1);
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
  if ((malen0 == alen1) && (!memcmp(mastrs_ptr[0], aptr1, alen1))) {
    if (((malen1 != alen2)) || memcmp(mastrs_ptr[1], aptr2, alen2)) {
      if (missing1) {
	if (allele_set(aptr2, alen2, &(mastrs_ptr[1]))) {
	  return RET_NOMEM;
	}
      } else {
        return RET_ALLELE_MISMATCH;
      }
    }
    *set_allele_freq_ptr = 1.0 - maf;
  } else if ((malen1 == alen1) && (!memcmp(mastrs_ptr[1], aptr1, alen1))) {
    if ((malen0 != alen2) || memcmp(mastrs_ptr[0], aptr2, alen2)) {
      if (missing0) {
        if (allele_set(aptr2, alen2, &(mastrs_ptr[0]))) {
	  return RET_NOMEM;
	}
      } else {
        return RET_ALLELE_MISMATCH;
      }
    }
    *set_allele_freq_ptr = maf;
  } else if (missing0 && (!missing1) && (malen1 == alen2) && (!memcmp(mastrs_ptr[1], aptr2, alen2))) {
    if (allele_set(aptr1, alen1, &(mastrs_ptr[0]))) {
      return RET_NOMEM;
    }
    *set_allele_freq_ptr = 1.0 - maf;
  } else if (missing1 && (!missing0) && (malen0 == alen2) && (!memcmp(mastrs_ptr[0], aptr2, alen2))) {
    if (allele_set(aptr1, alen1, &(mastrs_ptr[1]))) {
      return RET_NOMEM;
    }
    *set_allele_freq_ptr = maf;
  } else if ((*aptr1 == missing_geno) && (alen1 == 1) && (maf == 0.0)) {
    if ((malen0 == alen2) && (!memcmp(mastrs_ptr[0], aptr2, alen2))) {
      *set_allele_freq_ptr = 0.0;
    } else if ((malen1 == alen2) && (!memcmp(mastrs_ptr[1], aptr2, alen2))) {
      *set_allele_freq_ptr = 1.0;
    } else {
      return RET_ALLELE_MISMATCH;
    }
  } else {
    return RET_ALLELE_MISMATCH;
  }
  return 0;
}

uint32_t get_freq_file_type(char* bufptr) {
  uint32_t token_ct = count_tokens(bufptr);
  uint32_t slen;
  if ((token_ct == 6) || (token_ct == 7)) {
    slen = strlen_se(bufptr);
    if ((slen != 3) || memcmp(bufptr, "CHR", 3)) {
      return 0;
    }
    bufptr = skip_initial_spaces(&(bufptr[slen]));
    slen = strlen_se(bufptr);
    if ((slen != 3) || memcmp(bufptr, "SNP", 3)) {
      return 0;
    }
    bufptr = skip_initial_spaces(&(bufptr[slen]));
    slen = strlen_se(bufptr);
    if ((slen != 2) || memcmp(bufptr, "A1", 2)) {
      return 0;
    }
    bufptr = skip_initial_spaces(&(bufptr[slen]));
    slen = strlen_se(bufptr);
    if ((slen != 2) || memcmp(bufptr, "A2", 2)) {
      return 0;
    }
    bufptr = skip_initial_spaces(&(bufptr[slen]));
    slen = strlen_se(bufptr);
    if (token_ct == 6) {
      if ((slen != 3) || memcmp(bufptr, "MAF", 3)) {
	return 0;
      }
      bufptr = skip_initial_spaces(&(bufptr[slen]));
      slen = strlen_se(bufptr);
      if ((slen != 7) || memcmp(bufptr, "NCHROBS", 7)) {
	return 0;
      }
      return 1;
    } else {
      if ((slen != 2) || memcmp(bufptr, "C1", 2)) {
	return 0;
      }
      bufptr = skip_initial_spaces(&(bufptr[slen]));
      slen = strlen_se(bufptr);
      if ((slen != 2) || memcmp(bufptr, "C2", 2)) {
	return 0;
      }
      bufptr = skip_initial_spaces(&(bufptr[slen]));
      slen = strlen_se(bufptr);
      if ((slen != 2) || memcmp(bufptr, "G0", 2)) {
	return 0;
      }
      return 2;
    }
  } else if (token_ct == 3) {
    return 4;
  } else if (!memcmp(bufptr, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)", 71)) {
    return 3;
  }
  return 0;
}

int32_t read_external_freqs(char* freqname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, char** marker_allele_ptrs, double* set_allele_freqs, uint32_t* nchrobs, uint32_t maf_succ) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* freqfile = nullptr;
  uintptr_t line_idx = 0;
  uint32_t freq_counts = 0;
  uint32_t alen1 = 0;
  uint32_t alen2 = 0;
  uint32_t double_missing_ct = 0;
  uint32_t cur_nchrobs = 0;
  char* aptr1 = nullptr;
  char* aptr2 = nullptr;
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
  uint32_t marker_uidx;
  uint32_t uii;
  int32_t c_hom_a1;
  int32_t c_het;
  int32_t c_hom_a2;
  int32_t c_hap_a1;
  int32_t c_hap_a2;
  int32_t ii;
  {
    if (fopen_checked(freqname, "r", &freqfile)) {
      goto read_external_freqs_ret_OPEN_FAIL;
    }
    retval = sort_item_ids(unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref, &sorted_ids, &id_map);
    if (retval) {
      goto read_external_freqs_ret_1;
    }
    loadbuf_size = bigstack_left();
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size <= MAXLINELEN) {
      goto read_external_freqs_ret_NOMEM;
    }
    loadbuf = (char*)g_bigstack_base;
    loadbuf[loadbuf_size - 1] = ' ';
    do {
      if (!fgets(loadbuf, loadbuf_size, freqfile)) {
	logerrprint("Error: Empty --read-freq file.\n");
	goto read_external_freqs_ret_INVALID_FORMAT;
      }
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      bufptr = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*bufptr));
    uii = get_freq_file_type(bufptr);
    if (!uii) {
      logerrprint("Error: Invalid --read-freq file header.\n");
      goto read_external_freqs_ret_INVALID_FORMAT;
    }
    if (uii < 3) {
      if (uii == 2) {
	freq_counts = 1;
      }
      while (fgets(loadbuf, loadbuf_size, freqfile) != nullptr) {
	line_idx++;
	if (!loadbuf[loadbuf_size - 1]) {
	  goto read_external_freqs_ret_TOO_LONG_LINE;
	}
	char* loadbuf_first_token = skip_initial_spaces(loadbuf);
	if (is_eoln_kns(*loadbuf_first_token)) {
	  continue;
	}
	char* first_token_end = token_endnn(loadbuf_first_token);
	bufptr = skip_initial_spaces(first_token_end); // marker name
	const uint32_t chrom_name_slen = (uintptr_t)(first_token_end - loadbuf_first_token);
	*first_token_end = '\0';
	int32_t cur_chrom_code = get_chrom_code(loadbuf_first_token, chrom_info_ptr, chrom_name_slen);
	if (cur_chrom_code < 0) {
	  goto read_external_freqs_ret_INVALID_CHROM;
	}
	bufptr2 = next_token(bufptr);
	if (!bufptr2) {
	  goto read_external_freqs_ret_MISSING_TOKENS;
	}
	ii = bsearch_str(bufptr, strlen_se(bufptr), sorted_ids, max_marker_id_len, unfiltered_marker_ct - marker_exclude_ct);
	if (ii != -1) {
	  // may want to check for duplicates...
	  marker_uidx = id_map[(uint32_t)ii];
	  if ((((uint32_t)cur_chrom_code) == get_variant_chrom(chrom_info_ptr, marker_uidx)) || (!cur_chrom_code) || (!get_variant_chrom(chrom_info_ptr, marker_uidx))) {
	    if ((marker_allele_ptrs[marker_uidx * 2 + 1] == missing_geno_ptr) && (marker_allele_ptrs[marker_uidx * 2] == missing_geno_ptr)) {
	      double_missing_ct++;
	      continue;
	    }
	    alen1 = strlen_se(bufptr2);
	    aptr1 = bufptr2;
	    bufptr2 = next_token(bufptr2);
	    if (no_more_tokens_kns(bufptr2)) {
	      goto read_external_freqs_ret_MISSING_TOKENS;
	    }
	    alen2 = strlen_se(bufptr2);
	    aptr2 = bufptr2;
	    if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	      // permit A1='0', A2='0'
	      if ((*aptr1 == missing_geno) && (alen1 == 1)) {
		continue;
	      }
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
	      cur_nchrobs = c_hom_a1 + c_hom_a2;
	      maf = ((double)c_hom_a1 + maf_succ) / ((double)(cur_nchrobs + 2 * maf_succ));
	      if (nchrobs) {
		nchrobs[marker_uidx] = cur_nchrobs;
	      }
	    } else {
	      if (scan_double(bufptr, &maf)) {
		goto read_external_freqs_ret_INVALID_MAF;
	      }
	      if (nchrobs) {
		bufptr = next_token(bufptr);
		if (no_more_tokens_kns(bufptr)) {
		  goto read_external_freqs_ret_MISSING_TOKENS;
		}
		if (scan_uint_icap(bufptr, &cur_nchrobs)) {
		  goto read_external_freqs_ret_INVALID_NCHROBS;
		}
		nchrobs[marker_uidx] = cur_nchrobs;
	      }
	    }
	    retval = load_one_freq(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[marker_uidx]), &(marker_allele_ptrs[marker_uidx * 2]), missing_geno);
	    if (retval) {
	      goto read_external_freqs_ret_ALLELE_MISMATCH;
	    }
	  }
	}
      }
      if (freq_counts) {
	logprint("--read-freq: .frq.count file loaded.\n");
      } else {
	logprint("--read-freq: .frq file loaded.\n");
      }
    } else if (uii == 3) {
      // --freqx format
      while (fgets(loadbuf, loadbuf_size, freqfile) != nullptr) {
	line_idx++;
	if (!loadbuf[loadbuf_size - 1]) {
	  goto read_external_freqs_ret_TOO_LONG_LINE;
	}
	char* loadbuf_first_token = skip_initial_spaces(loadbuf);
	if (is_eoln_kns(*loadbuf_first_token)) {
	  continue;
	}
	char* first_token_end = token_endnn(loadbuf_first_token);
	bufptr = skip_initial_spaces(first_token_end);
	const uint32_t chrom_name_slen = (uintptr_t)(first_token_end - loadbuf_first_token);
	*first_token_end = '\0';
	int32_t cur_chrom_code = get_chrom_code(loadbuf_first_token, chrom_info_ptr, chrom_name_slen);
	if (cur_chrom_code < 0) {
	  goto read_external_freqs_ret_INVALID_CHROM;
	}
	bufptr2 = next_token(bufptr);
	if (!bufptr2) {
	  goto read_external_freqs_ret_MISSING_TOKENS;
	}
	ii = bsearch_str(bufptr, strlen_se(bufptr), sorted_ids, max_marker_id_len, unfiltered_marker_ct - marker_exclude_ct);
	if (ii != -1) {
	  marker_uidx = id_map[(uint32_t)ii];
	  if ((((uint32_t)cur_chrom_code) == get_variant_chrom(chrom_info_ptr, marker_uidx)) || (!cur_chrom_code) || (!get_variant_chrom(chrom_info_ptr, marker_uidx))) {
	    if ((marker_allele_ptrs[marker_uidx * 2 + 1] == missing_geno_ptr) && (marker_allele_ptrs[marker_uidx * 2] == missing_geno_ptr)) {
	      double_missing_ct++;
	      continue;
	    }
	    alen1 = strlen_se(bufptr2);
	    aptr1 = bufptr2;
	    bufptr2 = next_token(bufptr2);
	    if (no_more_tokens_kns(bufptr2)) {
	      goto read_external_freqs_ret_MISSING_TOKENS;
	    }
	    alen2 = strlen_se(bufptr2);
	    aptr2 = bufptr2;
	    if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	      if ((*aptr1 == missing_geno) && (alen1 == 1)) {
		continue;
	      }
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
	      sprintf(g_logbuf, "Error: Invalid het count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	      goto read_external_freqs_ret_INVALID_FORMAT_2;
	    }
	    if (scan_uint_icap(bufptr3, (uint32_t*)&c_hom_a2)) {
	      goto read_external_freqs_ret_INVALID_HOM_A2;
	    }
	    if (scan_uint_icap(bufptr4, (uint32_t*)&c_hap_a1)) {
	      sprintf(g_logbuf, "Error: Invalid hap. A1 count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	      goto read_external_freqs_ret_INVALID_FORMAT_2;
	    }
	    if (scan_uint_icap(bufptr5, (uint32_t*)&c_hap_a2)) {
	      sprintf(g_logbuf, "Error: Invalid hap. A2 count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	      goto read_external_freqs_ret_INVALID_FORMAT_2;
	    }
	    cur_nchrobs = 2 * (c_hom_a1 + c_het + c_hom_a2 + maf_succ) + c_hap_a1 + c_hap_a2;
	    maf = ((double)(c_hom_a1 * 2 + c_het + c_hap_a1 + maf_succ)) / ((double)cur_nchrobs);
	    if (nchrobs) {
	      nchrobs[marker_uidx] = cur_nchrobs;
	    }
	    retval = load_one_freq(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[marker_uidx]), &(marker_allele_ptrs[marker_uidx * 2]), missing_geno);
	    if (retval) {
	      goto read_external_freqs_ret_ALLELE_MISMATCH;
	    }
	  }
	}
      }
      logprint("--read-freq: .frqx file loaded.\n");
    } else {
      // Also support GCTA-style frequency files:
      // <marker ID>\t<reference allele>\t<frequency of reference allele>\n
      if (nchrobs) {
	logerrprint("Error: The current run requires an allele frequency file with observation\ncounts.\n");
	goto read_external_freqs_ret_INVALID_FORMAT;
      }

      // no header line here
      line_idx--;
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
	  if ((marker_allele_ptrs[marker_uidx * 2 + 1] == missing_geno_ptr) && (marker_allele_ptrs[marker_uidx * 2] == missing_geno_ptr)) {
	    double_missing_ct++;
	    continue;
	  }
	  alen1 = strlen_se(bufptr);
	  aptr1 = bufptr;
	  bufptr = next_token(bufptr);
	  if (no_more_tokens_kns(bufptr)) {
	    goto read_external_freqs_ret_MISSING_TOKENS;
	  }
	  if (scan_double(bufptr, &maf)) {
	    goto read_external_freqs_ret_INVALID_MAF;
	  }
	  retval = load_one_freq(1, missing_geno_ptr, alen1, aptr1, maf, &(set_allele_freqs[marker_uidx]), &(marker_allele_ptrs[marker_uidx * 2]), missing_geno);
	  if (retval) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
	  }
	} else {
	  // if there aren't exactly 3 columns, this isn't a GCTA .freq file
	  bufptr = next_token(bufptr);
	  if (no_more_tokens_kns(bufptr)) {
	    goto read_external_freqs_ret_MISSING_TOKENS;
	  }
	  if (!no_more_tokens_kns(next_token(bufptr))) {
	    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --read-freq has more tokens than expected.\n", line_idx);
	    goto read_external_freqs_ret_INVALID_FORMAT_2;
	  }
	}
      } while (fgets(loadbuf, loadbuf_size, freqfile));
      logprint("--read-freq: GCTA-formatted .freq file loaded.\n");
    }
    if (double_missing_ct) {
      LOGPRINTF("%u variant%s skipped since both existing allele codes were missing.\n", double_missing_ct, (double_missing_ct == 1)? "" : "s");
    }
  }
  while (0) {
  read_external_freqs_ret_TOO_LONG_LINE:
    if (loadbuf_size == MAXLINEBUFLEN) {
      LOGERRPRINTF("Error: Line %" PRIuPTR " of --read-freq file is pathologically long.\n", line_idx);
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
    LOGERRPRINTF("Error: Invalid hom. A1 count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_INVALID_HOM_A2:
    LOGERRPRINTF("Error: Invalid hom. A2 count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_INVALID_NCHROBS:
    LOGERRPRINTF("Error: Invalid NCHROBS value on line %" PRIuPTR " of --read-freq file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_INVALID_MAF:
    LOGERRPRINTF("Error: Invalid MAF on line %" PRIuPTR " of --read-freq file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_A1_A2_SAME:
    LOGERRPRINTF("Error: A1 and A2 alleles match on line %" PRIuPTR " of --read-freq file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_INVALID_CHROM:
    sprintf(g_logbuf, "Error: Invalid chromosome code on line %" PRIuPTR" of --read-freq file.\n", line_idx);
  read_external_freqs_ret_INVALID_FORMAT_2:
    logerrprintb();
  read_external_freqs_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_MISSING_TOKENS:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --read-freq file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_ALLELE_MISMATCH:
    if (retval == RET_NOMEM) {
      break;
    }
    LOGERRPRINTF("Error: Allele(s) on line %" PRIuPTR " of --read-freq file don't match loaded\nvalues.\n", line_idx);
    break;
  }
 read_external_freqs_ret_1:
  fclose_cond(freqfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t load_ax_alleles(Two_col_params* axalleles, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, double* set_allele_freqs, uint32_t is_a2) {
  // note that swap_reversed_marker_alleles() has NOT been called yet
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  char skipchar = axalleles->skipchar;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  uint32_t colid_first = (axalleles->colid < axalleles->colx);
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t max_marker_allele_len = *max_marker_allele_len_ptr;
  uintptr_t* already_seen;
  char* loadbuf;
  char* colid_ptr;
  char* colx_ptr;
  uint32_t* marker_id_htable;
  uintptr_t loadbuf_size;
  uintptr_t line_idx;
  uint32_t marker_id_htable_size;
  uint32_t idlen;
  uint32_t alen;
  uint32_t colmin;
  uint32_t coldiff;
  uint32_t marker_uidx;
  char cc;
  uint32_t replace_other;
  int32_t retval;
  retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_len, 0, &marker_id_htable_size, &marker_id_htable);
  if (retval) {
    goto load_ax_alleles_ret_1;
  }
  if (bigstack_calloc_ul(unfiltered_marker_ctl, &already_seen)) {
    goto load_ax_alleles_ret_NOMEM;
  }
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_ax_alleles_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  retval = open_and_skip_first_lines(&infile, axalleles->fname, loadbuf, loadbuf_size, axalleles->skip);
  if (retval) {
    goto load_ax_alleles_ret_1;
  }
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
      if (!colx_ptr) {
	goto load_ax_alleles_ret_MISSING_TOKENS;
      }
    } else {
      colx_ptr = next_token_multz(colid_ptr, colmin);
      colid_ptr = next_token_mult(colx_ptr, coldiff);
      if (!colid_ptr) {
	goto load_ax_alleles_ret_MISSING_TOKENS;
      }
    }
    idlen = strlen_se(colid_ptr);
    marker_uidx = id_htable_find(colid_ptr, idlen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if (marker_uidx == 0xffffffffU) {
      continue;
    }
    if (is_set(already_seen, marker_uidx)) {
      colid_ptr[idlen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --a%c-allele file.\n", colid_ptr, is_a2? '2' : '1');
      goto load_ax_alleles_ret_INVALID_FORMAT_2;
    }
    set_bit(marker_uidx, already_seen);
    alen = strlen_se(colx_ptr);
    colx_ptr[alen] = '\0';
    if (!strcmp(colx_ptr, marker_allele_ptrs[marker_uidx * 2 + is_a2])) {
      if (IS_SET(marker_reverse, marker_uidx)) {
        set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
        CLEAR_BIT(marker_uidx, marker_reverse);
      }
    } else if (!strcmp(colx_ptr, marker_allele_ptrs[marker_uidx * 2 + 1 - is_a2])) {
      if (!IS_SET(marker_reverse, marker_uidx)) {
        set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
        SET_BIT(marker_uidx, marker_reverse);
      }
    } else if ((marker_allele_ptrs[marker_uidx * 2] == missing_geno_ptr) || (marker_allele_ptrs[marker_uidx * 2 + 1] == missing_geno_ptr)) {
      replace_other = (marker_allele_ptrs[marker_uidx * 2 + is_a2] != missing_geno_ptr);
      if (allele_reset(colx_ptr, alen, &(marker_allele_ptrs[marker_uidx * 2 + (is_a2 ^ replace_other)]))) {
	goto load_ax_alleles_ret_NOMEM;
      }
      if (alen >= max_marker_allele_len) {
	max_marker_allele_len = alen + 1;
      }
      if (!replace_other) {
	if (IS_SET(marker_reverse, marker_uidx)) {
	  set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
	  CLEAR_BIT(marker_uidx, marker_reverse);
	}
      } else {
	if (!IS_SET(marker_reverse, marker_uidx)) {
	  set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
	  SET_BIT(marker_uidx, marker_reverse);
	}
      }
    } else {
      colid_ptr[idlen] = '\0';
      LOGERRPRINTF("Warning: Impossible A%c allele assignment for variant %s.\n", is_a2? '2' : '1', colid_ptr);
    }
  }
  if (!feof(infile)) {
    goto load_ax_alleles_ret_READ_FAIL;
  }
  marker_uidx = popcount_longs(already_seen, unfiltered_marker_ctl);
  LOGPRINTF("--a%c-allele: %u assignment%s made.\n", is_a2? '2' : '1', marker_uidx, (marker_uidx == 1)? "" : "s");
  *max_marker_allele_len_ptr = max_marker_allele_len;
  while (0) {
  load_ax_alleles_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_ax_alleles_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_ax_alleles_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Fewer tokens than expected on line %" PRIuPTR " of --a%c-allele file.\n", line_idx, is_a2? '2' : '1');
    wordwrapb(0);
  load_ax_alleles_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 load_ax_alleles_ret_1:
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t write_stratified_freqs(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t output_gz, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_sample_ct, uintptr_t sample_ct, uint32_t sample_f_ct, uintptr_t* founder_info, uint32_t nonfounders, uintptr_t* sex_male, uint32_t sample_f_male_ct, uintptr_t* marker_reverse, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len) {
  // unfiltered_sample_ct == 0 ok
  unsigned char* bigstack_mark = g_bigstack_base;
  char* writebuf = g_textbuf;
  char* pzwritep = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uint32_t* cur_cluster_map = cluster_map;
  uint32_t* cur_cluster_starts = cluster_starts;
  uint32_t* cluster_map_nonmale = nullptr;
  uint32_t* cluster_starts_nonmale = nullptr;
  uint32_t* cluster_map_male = nullptr;
  uint32_t* cluster_starts_male = nullptr;
  int32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  uint32_t cslen = 10;
  int32_t retval = 0;
  uint32_t cur_cts[4];
  Pigz_state ps;
  uintptr_t* readbuf;
  uint32_t* uiptr;
  uint32_t* uiptr2;
  uint32_t* uiptr3;
  unsigned char* overflow_buf;
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
  pzwrite_init_null(&ps);
  uii = 2 * max_marker_allele_len + max_marker_id_len + max_cluster_id_len + 256;
  if (bigstack_alloc_uc(uii + PIGZ_BLOCK_SIZE, &overflow_buf) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &readbuf)) {
    goto write_stratified_freqs_ret_NOMEM;
  }
  if (uii > MAXLINELEN) {
    if (bigstack_alloc_c(uii, &writebuf)) {
      goto write_stratified_freqs_ret_NOMEM;
    }
  }
  if ((sample_ct > sample_f_ct) && (!nonfounders)) {
    if (bigstack_alloc_ui(cluster_ct + 1, &cur_cluster_starts) ||
        bigstack_alloc_ui(sample_f_ct, &cur_cluster_map)) {
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
  chrom_idx = chrom_info_ptr->xymt_codes[X_OFFSET];
  if ((chrom_idx != -2) && is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
    if (bigstack_alloc_ui(cluster_ct + 1, &cluster_starts_nonmale) ||
        bigstack_alloc_ui(sample_f_ct - sample_f_male_ct, &cluster_map_nonmale)) {
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
  chrom_idx = chrom_info_ptr->xymt_codes[Y_OFFSET];
  if (cluster_map_nonmale || ((chrom_idx != -2) && is_set(chrom_info_ptr->chrom_mask, chrom_idx))) {
    if (bigstack_alloc_ui(cluster_ct + 1, &cluster_starts_male) ||
        bigstack_alloc_ui(sample_f_male_ct, &cluster_map_male)) {
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
  memcpy(outname_end, output_gz? ".frq.strat.gz" : ".frq.strat", output_gz? 14 : 11);
  if (flex_pzwrite_init(output_gz, outname, overflow_buf, 0, &ps)) {
    goto write_stratified_freqs_ret_OPEN_FAIL;
  }
  pzwritep = (char*)overflow_buf;
  sprintf(g_textbuf, " CHR %%%us     CLST   A1   A2      MAF    MAC  NCHROBS" EOLN_STR, plink_maxsnp);
  pzwritep += sprintf(pzwritep, g_textbuf, "SNP");
  if (bigstack_alloc_c(2 * max_marker_allele_len + 16, &csptr)) {
    goto write_stratified_freqs_ret_NOMEM;
  }
  memset(csptr, 32, 10);
  for (chrom_idx = 0; chrom_idx < chrom_code_end; chrom_idx++) {
    if (!is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == chrom_info_ptr->xymt_codes[X_OFFSET]);
    is_y = (chrom_idx == chrom_info_ptr->xymt_codes[Y_OFFSET]);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    chrom_end = get_chrom_end_vidx(chrom_info_ptr, chrom_idx);
    marker_uidx = next_unset_ul(marker_exclude, get_chrom_start_vidx(chrom_info_ptr, chrom_idx), chrom_end);
    if (marker_uidx >= chrom_end) {
      continue;
    }
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
      goto write_stratified_freqs_ret_READ_FAIL;
    }
    col_2_start = width_force(4, writebuf, chrom_name_write(chrom_info_ptr, chrom_idx, writebuf));
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

      if (load_raw(unfiltered_sample_ct4, bedfile, readbuf)) {
	goto write_stratified_freqs_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)readbuf);
      }
      if (is_x) {
	uiptr = cluster_map_nonmale;
	uiptr2 = cluster_map_male;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
          pzwritep = memcpya(pzwritep, writebuf, wptr_start - writebuf);
	  pzwritep = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), pzwritep);
	  pzwritep = memcpya(pzwritep, csptr, cslen);
	  fill_uint_zero(4, cur_cts);
	  uiptr3 = &(cluster_map_nonmale[cluster_starts_nonmale[clidx + 1]]);
	  while (uiptr < uiptr3) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  a1_obs = 2 * cur_cts[0] + cur_cts[2];
	  tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  fill_uint_zero(4, cur_cts);
	  uiptr3 = &(cluster_map_male[cluster_starts_male[clidx + 1]]);
	  while (uiptr2 < uiptr3) {
	    uii = *uiptr2++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  a1_obs += cur_cts[0];
	  tot_obs += cur_cts[0] + cur_cts[3];
	  if (tot_obs) {
            pzwritep = dtoa_g_wxp4x(((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ', pzwritep);
	    pzwritep = uint32toa_w6x(a1_obs, ' ', pzwritep);
	    pzwritep = uint32toa_w8x(tot_obs, ' ', pzwritep);
	  } else {
	    pzwritep = memcpya(pzwritep, "       0      0        0 ", 25);
	  }
	  append_binary_eoln(&pzwritep);
	  if (flex_pzwrite(&ps, &pzwritep)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      } else if (is_y) {
	uiptr = cluster_map_male;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  pzwritep = memcpya(pzwritep, writebuf, wptr_start - writebuf);
	  pzwritep = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), pzwritep);
	  pzwritep = memcpya(pzwritep, csptr, cslen);
	  fill_uint_zero(4, cur_cts);
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
            pzwritep = dtoa_g_wxp4x(((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ', pzwritep);
	    pzwritep = uint32toa_w6x(a1_obs, ' ', pzwritep);
	    pzwritep = uint32toa_w8x(tot_obs, ' ', pzwritep);
	  } else {
	    pzwritep = memcpya(pzwritep, "       0      0        0 ", 25);
	  }
	  append_binary_eoln(&pzwritep);
	  if (flex_pzwrite(&ps, &pzwritep)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      } else {
        uiptr = cur_cluster_map;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  pzwritep = memcpya(pzwritep, writebuf, wptr_start - writebuf);
	  pzwritep = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), pzwritep);
	  pzwritep = memcpya(pzwritep, csptr, cslen);
	  fill_uint_zero(4, cur_cts);
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
            pzwritep = dtoa_g_wxp4x(((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ', pzwritep);
	    pzwritep = uint32toa_w6x(a1_obs, ' ', pzwritep);
	    pzwritep = uint32toa_w8x(tot_obs, ' ', pzwritep);
	  } else {
	    pzwritep = memcpya(pzwritep, "       0      0        0 ", 25);
	  }
	  append_binary_eoln(&pzwritep);
	  if (flex_pzwrite(&ps, &pzwritep)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto write_stratified_freqs_ret_READ_FAIL;
	}
      }
    } while (marker_uidx < chrom_end);
  }
  if (flex_pzwrite_close_null(&ps, pzwritep)) {
    goto write_stratified_freqs_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--freq: Cluster-stratified allele frequencies (%s) written to %s .\n", nonfounders? "all samples" : "founders only", outname);
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
  bigstack_reset(bigstack_mark);
  flex_pzwrite_close_cond(&ps, pzwritep);
  return retval;
}

int32_t write_cc_freqs(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t output_gz, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uint32_t nonfounders, uintptr_t* sex_male, uintptr_t* marker_reverse, uintptr_t* pheno_nm, uintptr_t* pheno_c) {
  // unfiltered_sample_ct must be positive
  unsigned char* bigstack_mark = g_bigstack_base;
  char* pzwritep = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t* loadbuf = nullptr;
  uintptr_t* case_include2 = nullptr;
  uintptr_t* ctrl_include2 = nullptr;
  uintptr_t* male_vec = nullptr;
  uintptr_t* nonmale_vec = nullptr;
  int32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  int32_t retval = 0;
  Pigz_state ps;
  unsigned char* overflow_buf;
  uintptr_t* ulptr;
  char* bufptr;
  uintptr_t chrom_end;
  uintptr_t marker_uidx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  int32_t chrom_idx;
  uint32_t case_ct;
  uint32_t ctrl_ct;
  uint32_t case_set_ct;
  uint32_t ctrl_set_ct;
  uint32_t case_missing_ct;
  uint32_t ctrl_missing_ct;
  uint32_t case_obs;
  uint32_t ctrl_obs;
  uint32_t uii;
  pzwrite_init_null(&ps);
  if (bigstack_alloc_uc(PIGZ_BLOCK_SIZE + 2 * max_marker_allele_len + max_marker_id_len + 256, &overflow_buf)) {
    goto write_cc_freqs_ret_NOMEM;
  }

  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &case_include2) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &ctrl_include2) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &male_vec) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &nonmale_vec)) {
      goto write_cc_freqs_ret_NOMEM;
  }
  loadbuf[unfiltered_sample_ctl2 - 1] = 0;
  init_quaterarr_from_bitarr(pheno_c, unfiltered_sample_ct, case_include2);
  init_quaterarr_from_bitarr(pheno_nm, unfiltered_sample_ct, ctrl_include2);
  memcpy(nonmale_vec, ctrl_include2, unfiltered_sample_ctl2 * sizeof(intptr_t));
  bitvec_andnot(case_include2, unfiltered_sample_ctl2, ctrl_include2);
  init_quaterarr_from_bitarr(sex_male, unfiltered_sample_ct, male_vec);
  bitvec_andnot(male_vec, unfiltered_sample_ctl2, nonmale_vec);
  if (!nonfounders) {
    if (bigstack_alloc_ul(unfiltered_sample_ctl2, &ulptr)) {
      goto write_cc_freqs_ret_NOMEM;
    }
    init_quaterarr_from_bitarr(founder_info, unfiltered_sample_ct, ulptr);
    bitvec_and(ulptr, unfiltered_sample_ctl2, case_include2);
    bitvec_and(ulptr, unfiltered_sample_ctl2, ctrl_include2);
    bitvec_and(ulptr, unfiltered_sample_ctl2, male_vec);
    bitvec_and(ulptr, unfiltered_sample_ctl2, nonmale_vec);
    bigstack_reset(ulptr);
  }
  case_ct = popcount2_longs(case_include2, unfiltered_sample_ctl2);
  ctrl_ct = popcount2_longs(ctrl_include2, unfiltered_sample_ctl2);

  bufptr = memcpya(outname_end, ".frq.cc", 7);
  if (!output_gz) {
    *bufptr = '\0';
  } else {
    memcpy(bufptr, ".gz", 4);
  }
  if (flex_pzwrite_init(output_gz, outname, overflow_buf, 0, &ps)) {
    goto write_cc_freqs_ret_OPEN_FAIL;
  }
  pzwritep = (char*)overflow_buf;
  pzwritep = memcpya(pzwritep, " CHR  ", 6);
  if (plink_maxsnp >= 5) {
    pzwritep = memseta(pzwritep, 32, plink_maxsnp - 4);
  }
  pzwritep = strcpya(pzwritep, "SNP   A1   A2        MAF_A        MAF_U  NCHROBS_A  NCHROBS_U" EOLN_STR);

  for (chrom_idx = 0; chrom_idx < chrom_code_end; chrom_idx++) {
    if (!is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == chrom_info_ptr->xymt_codes[X_OFFSET]);
    is_y = (chrom_idx == chrom_info_ptr->xymt_codes[Y_OFFSET]);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    chrom_end = get_chrom_end_vidx(chrom_info_ptr, chrom_idx);
    marker_uidx = next_unset_ul(marker_exclude, get_chrom_start_vidx(chrom_info_ptr, chrom_idx), chrom_end);
    if (marker_uidx >= chrom_end) {
      continue;
    }
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
      goto write_cc_freqs_ret_READ_FAIL;
    }
    do {
      pzwritep = width_force(4, pzwritep, chrom_name_write(chrom_info_ptr, chrom_idx, pzwritep));
      *pzwritep++ = ' ';
      pzwritep = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), pzwritep);
      *pzwritep++ = ' ';
      pzwritep = fw_strcpy(4, marker_allele_ptrs[marker_uidx * 2], pzwritep);
      *pzwritep++ = ' ';
      pzwritep = fw_strcpy(4, marker_allele_ptrs[marker_uidx * 2 + 1], pzwritep);
      *pzwritep++ = ' ';

      if (load_raw(unfiltered_sample_ct4, bedfile, loadbuf)) {
	goto write_cc_freqs_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf);
      }
      if (is_x) {
	genovec_set_freq_x(loadbuf, case_include2, male_vec, unfiltered_sample_ctl2, &case_set_ct, &case_missing_ct);
	genovec_set_freq_x(loadbuf, ctrl_include2, male_vec, unfiltered_sample_ctl2, &ctrl_set_ct, &ctrl_missing_ct);
	case_obs = 2 * case_ct - case_missing_ct;
	ctrl_obs = 2 * ctrl_ct - ctrl_missing_ct;
      } else if (!is_haploid) {
	genovec_set_freq(loadbuf, case_include2, unfiltered_sample_ctl2, &case_set_ct, &case_missing_ct);
	genovec_set_freq(loadbuf, ctrl_include2, unfiltered_sample_ctl2, &ctrl_set_ct, &ctrl_missing_ct);
	case_obs = 2 * (case_ct - case_missing_ct);
	ctrl_obs = 2 * (ctrl_ct - ctrl_missing_ct);
      } else {
        if (is_y) {
	  genovec_set_freq_y(loadbuf, case_include2, nonmale_vec, unfiltered_sample_ctl2, &case_set_ct, &case_missing_ct);
	  genovec_set_freq_y(loadbuf, ctrl_include2, nonmale_vec, unfiltered_sample_ctl2, &ctrl_set_ct, &ctrl_missing_ct);
        } else {
	  genovec_3freq(loadbuf, case_include2, unfiltered_sample_ctl2, &case_missing_ct, &uii, &case_set_ct);
	  case_missing_ct += uii;
	  genovec_3freq(loadbuf, ctrl_include2, unfiltered_sample_ctl2, &ctrl_missing_ct, &uii, &ctrl_set_ct);
	  ctrl_missing_ct += uii;
	}
	case_obs = case_ct - case_missing_ct;
	ctrl_obs = ctrl_ct - ctrl_missing_ct;
      }
      if (case_obs) {
	pzwritep = dtoa_g_wxp4x(((double)((int32_t)(case_obs - case_set_ct))) / ((double)case_obs), 12, ' ', pzwritep);
      } else {
	pzwritep = memcpya(pzwritep, "        NA ", 11);
      }
      if (ctrl_obs) {
	pzwritep = dtoa_g_wxp4x(((double)((int32_t)(ctrl_obs - ctrl_set_ct))) / ((double)ctrl_obs), 12, ' ', pzwritep);
      } else {
	pzwritep = memcpya(pzwritep, "        NA ", 11);
      }

      pzwritep = uint32toa_w10x(case_obs, ' ', pzwritep);
      pzwritep = uint32toa_w10(ctrl_obs, pzwritep);
      append_binary_eoln(&pzwritep);
      if (flex_pzwrite(&ps, &pzwritep)) {
	goto write_cc_freqs_ret_WRITE_FAIL;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto write_cc_freqs_ret_WRITE_FAIL;
	}
      }
    } while (marker_uidx < chrom_end);
  }
  if (flex_pzwrite_close_null(&ps, pzwritep)) {
    goto write_cc_freqs_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--freq case-control: Phenotype-stratified allele frequencies (%s) written to %s .\n", nonfounders? "all samples" : "founders only", outname);
  while (0) {
  write_cc_freqs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_cc_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_cc_freqs_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  write_cc_freqs_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  flex_pzwrite_close_cond(&ps, pzwritep);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t write_freqs(char* outname, char* outname_end, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, double* set_allele_freqs, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, int32_t* hapl_cts, int32_t* haph_cts, uint32_t sample_f_ct, uint32_t sample_f_male_ct, uint32_t nonfounders, uint64_t misc_flags, uintptr_t* marker_reverse) {
  // unfiltered_sample_ct == 0 ok
  unsigned char* bigstack_mark = g_bigstack_base;
  char* pzwritep = nullptr;
  uint32_t reverse = 0;
  uint32_t freq_counts = (misc_flags / MISC_FREQ_COUNTS) & 1;
  uint32_t freqx = (misc_flags / MISC_FREQX) & 1;
  uint32_t output_gz = (misc_flags / MISC_FREQ_GZ) & 1;
  uint32_t maf_succ = (misc_flags / MISC_MAF_SUCC) & 1;
  int32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  int32_t retval = 0;
  Pigz_state ps;
  unsigned char* overflow_buf;
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
  pzwrite_init_null(&ps);
  if (bigstack_alloc_uc(PIGZ_BLOCK_SIZE + 2 * max_marker_allele_len + MAXLINELEN, &overflow_buf)) {
    goto write_freqs_ret_NOMEM;
  }

  bufptr = memcpya(outname_end, ".frq", 4);
  if (freqx) {
    *bufptr++ = 'x';
  } else if (freq_counts) {
    bufptr = memcpya(bufptr, ".counts", 7);
  }
  if (!output_gz) {
    *bufptr = '\0';
  } else {
    memcpy(bufptr, ".gz", 4);
  }
  if (flex_pzwrite_init(output_gz, outname, overflow_buf, 0, &ps)) {
    goto write_freqs_ret_OPEN_FAIL;
  }
  pzwritep = (char*)overflow_buf;
  if (freqx) {
    pzwritep = strcpya(pzwritep, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)" EOLN_STR);
  } else {
    pzwritep = memcpya(pzwritep, " CHR  ", 6);
    if (plink_maxsnp >= 5) {
      pzwritep = memseta(pzwritep, 32, plink_maxsnp - 4);
    }
    if (freq_counts) {
      pzwritep = strcpya(pzwritep, "SNP   A1   A2     C1     C2     G0" EOLN_STR);
    } else {
      pzwritep = strcpya(pzwritep, "SNP   A1   A2          MAF  NCHROBS" EOLN_STR);
    }
  }
  for (chrom_idx = 0; chrom_idx < chrom_code_end; chrom_idx++) {
    if (!is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == chrom_info_ptr->xymt_codes[X_OFFSET]);
    is_y = (chrom_idx == chrom_info_ptr->xymt_codes[Y_OFFSET]);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    chrom_end = get_chrom_end_vidx(chrom_info_ptr, chrom_idx);
    marker_uidx = next_unset(marker_exclude, get_chrom_start_vidx(chrom_info_ptr, chrom_idx), chrom_end);
    while (marker_uidx < chrom_end) {
      reverse = IS_SET(marker_reverse, marker_uidx);
      major_ptr = marker_allele_ptrs[marker_uidx * 2 + 1];
      minor_ptr = marker_allele_ptrs[marker_uidx * 2];
      if (freq_counts || freqx) {
	if (is_x) {
	  missing_ct = sample_f_ct - (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx] + hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	} else if (is_haploid) {
	  if (is_y) {
	    missing_ct = sample_f_male_ct - (hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	  } else {
	    missing_ct = sample_f_ct - (hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	  }
	} else {
	  missing_ct = sample_f_ct - (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx]);
	}
	if (freqx) {
	  pzwritep = chrom_name_write(chrom_info_ptr, chrom_idx, pzwritep);
	  *pzwritep++ = '\t';
	  pzwritep = strcpyax(pzwritep, &(marker_ids[marker_uidx * max_marker_id_len]), '\t');
	  pzwritep = strcpyax(pzwritep, minor_ptr, '\t');
          pzwritep = strcpyax(pzwritep, major_ptr, '\t');
          pzwritep = uint32toa_x(reverse? hh_cts[marker_uidx] : ll_cts[marker_uidx], '\t', pzwritep);
	  pzwritep = uint32toa_x(lh_cts[marker_uidx], '\t', pzwritep);
          pzwritep = uint32toa_x(reverse? ll_cts[marker_uidx] : hh_cts[marker_uidx], '\t', pzwritep);
          pzwritep = uint32toa_x(reverse? haph_cts[marker_uidx] : hapl_cts[marker_uidx], '\t', pzwritep);
          pzwritep = uint32toa_x(reverse? hapl_cts[marker_uidx] : haph_cts[marker_uidx], '\t', pzwritep);
          pzwritep = uint32toa(missing_ct, pzwritep);
	} else {
	  pzwritep = width_force(4, pzwritep, chrom_name_write(chrom_info_ptr, chrom_idx, pzwritep));
	  *pzwritep++ = ' ';
	  pzwritep = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), pzwritep);
	  *pzwritep++ = ' ';
	  pzwritep = fw_strcpy(4, minor_ptr, pzwritep);
	  *pzwritep++ = ' ';
	  pzwritep = fw_strcpy(4, major_ptr, pzwritep);
	  *pzwritep++ = ' ';
          // bugfix (13 Oct 2017): did not take reverse into account here.
          const uint32_t l_ct = 2 * ll_cts[marker_uidx] + lh_cts[marker_uidx] + hapl_cts[marker_uidx];
          const uint32_t h_ct = 2 * hh_cts[marker_uidx] + lh_cts[marker_uidx] + haph_cts[marker_uidx];
          pzwritep = uint32toa_w6x(reverse? h_ct : l_ct, ' ', pzwritep);
	  pzwritep = uint32toa_w6x(reverse? l_ct : h_ct, ' ', pzwritep);
	  pzwritep = uint32toa_w6(missing_ct, pzwritep);
	}
      } else {
	pzwritep = width_force(4, pzwritep, chrom_name_write(chrom_info_ptr, chrom_idx, pzwritep));
	*pzwritep++ = ' ';
	pzwritep = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), pzwritep);
        *pzwritep++ = ' ';
	pzwritep = fw_strcpy(4, minor_ptr, pzwritep);
	*pzwritep++ = ' ';
	pzwritep = fw_strcpy(4, major_ptr, pzwritep);
	*pzwritep++ = ' ';
	uii = 2 * (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx]) + hapl_cts[marker_uidx] + haph_cts[marker_uidx];
        // set_allele_freqs[] already takes reverse into account.
	if (maf_succ || uii || (set_allele_freqs[marker_uidx] != 0.5)) {
	  pzwritep = dtoa_g_wxp4(1.0 - set_allele_freqs[marker_uidx], 12, pzwritep);
	} else {
	  pzwritep = memcpya(pzwritep, "          NA", 12);
	}
	*pzwritep++ = ' ';
	pzwritep = uint32toa_w8(uii, pzwritep);
      }
      append_binary_eoln(&pzwritep);
      if (flex_pzwrite(&ps, &pzwritep)) {
	goto write_freqs_ret_WRITE_FAIL;
      }
      marker_uidx = next_unset(marker_exclude, marker_uidx + 1, chrom_end);
    }
  }
  if (flex_pzwrite_close_null(&ps, pzwritep)) {
    goto write_freqs_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--freq%s%s: Allele frequencies (%s) written to %s .\n", freqx? "x" : "", freq_counts? " counts" : "", nonfounders? "all samples" : "founders only", outname);
  while (0) {
  write_freqs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_freqs_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  flex_pzwrite_close_cond(&ps, pzwritep);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t sexcheck(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uint64_t misc_flags, double check_sex_fthresh, double check_sex_mthresh, uint32_t max_f_yobs, uint32_t min_m_yobs, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* gender_unk_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uint32_t* het_cts = nullptr;
  uint32_t* missing_cts = nullptr;
  double* nei_offsets = nullptr;
  uint32_t* ymiss_cts = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  uintptr_t final_mask = get_final_mask(sample_ct);
  uintptr_t x_variant_ct = 0;
  uintptr_t ytotal = 0;
  double nei_sum = 0.0;
  uint32_t do_impute = (misc_flags / MISC_IMPUTE_SEX) & 1;
  uint32_t check_y = (misc_flags / MISC_SEXCHECK_YCOUNT) & 1;
  uint32_t yonly = (misc_flags / MISC_SEXCHECK_YONLY) & 1;
  uint32_t gender_unk_ct = 0;
  uint32_t problem_ct = 0;
  int32_t x_code = chrom_info_ptr->xymt_codes[X_OFFSET];
  int32_t y_code = chrom_info_ptr->xymt_codes[Y_OFFSET];
  int32_t retval = 0;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  // We wish to compute four quantities for each sample:
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
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t cur_missing_ct;
  uintptr_t allele_obs_ct;
  uintptr_t cur_word;
  uintptr_t ulii;
  uint32_t orig_sex_code;
  uint32_t imputed_sex_code;
  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(sample_ctl2, &loadbuf)) {
    goto sexcheck_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  loadbuf[sample_ctl2 - 1] = 0;
  if (check_y) {
    if (bigstack_calloc_ui(sample_ct, &ymiss_cts)) {
      goto sexcheck_ret_NOMEM;
    }
  }
  if (!yonly) {
    if ((x_code == -2) || (!is_set(chrom_info_ptr->chrom_mask, (uint32_t)x_code))) {
      goto sexcheck_ret_NO_X_VAR;
    }
    marker_uidx_end = get_chrom_end_vidx(chrom_info_ptr, (uint32_t)x_code);
    marker_uidx = next_unset_ul(marker_exclude, get_chrom_start_vidx(chrom_info_ptr, (uint32_t)x_code), marker_uidx_end);
    if (marker_uidx == marker_uidx_end) {
      goto sexcheck_ret_NO_X_VAR;
    }
    marker_idxs_left = marker_uidx_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, marker_uidx_end);
    if (bigstack_calloc_ui(sample_ct, &het_cts) ||
        bigstack_calloc_ui(sample_ct, &missing_cts) ||
        bigstack_calloc_d(sample_ct, &nei_offsets)) {
      goto sexcheck_ret_NOMEM;
    }
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
      goto sexcheck_ret_READ_FAIL;
    }
    for (; marker_idxs_left; marker_idxs_left--, marker_uidx++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto sexcheck_ret_READ_FAIL;
	}
      }
      if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, 0, bedfile, loadbuf_raw, loadbuf)) {
	goto sexcheck_ret_READ_FAIL;
      }
      cur_missing_ct = count_01(loadbuf, sample_ctl2);
      allele_obs_ct = 2 * (sample_ct - cur_missing_ct);
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
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx += BITCT2) {
	  cur_word = *lptr++;
	  cur_word = (~(cur_word >> 1)) & cur_word & FIVEMASK;
	  while (cur_word) {
	    ulii = sample_idx + CTZLU(cur_word) / 2;
	    missing_cts[ulii] += 1;
	    nei_offsets[ulii] += cur_nei;
	    cur_word &= cur_word - 1;
	  }
	}
      }
      lptr = loadbuf;
      // iterate through heterozygous calls
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx += BITCT2) {
	cur_word = *lptr++;
	cur_word = (cur_word >> 1) & (~cur_word) & FIVEMASK;
	while (cur_word) {
	  het_cts[sample_idx + CTZLU(cur_word) / 2] += 1;
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
    if ((y_code != -2) && is_set(chrom_info_ptr->chrom_mask, (uint32_t)y_code)) {
      marker_uidx_end = get_chrom_end_vidx(chrom_info_ptr, (uint32_t)y_code);
      marker_uidx = next_unset_ul(marker_exclude, get_chrom_start_vidx(chrom_info_ptr, (uint32_t)y_code), marker_uidx_end);
      ytotal = marker_uidx_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, marker_uidx_end);
    }
    if (ytotal) {
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	goto sexcheck_ret_READ_FAIL;
      }
      for (ulii = 0; ulii < ytotal; ulii++, marker_uidx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto sexcheck_ret_READ_FAIL;
	  }
	}
	if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, 0, bedfile, loadbuf_raw, loadbuf)) {
	  goto sexcheck_ret_READ_FAIL;
	}
	lptr = loadbuf;
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx += BITCT2) {
	  // iterate through missing calls
	  // vertical popcount would be faster, but we don't really care
	  cur_word = *lptr++;
	  cur_word = (~(cur_word >> 1)) & cur_word & FIVEMASK;
	  while (cur_word) {
	    ymiss_cts[sample_idx + CTZLU(cur_word) / 2] += 1;
	    cur_word &= cur_word - 1;
	  }
	}
      }
    } else if (yonly) {
      logerrprint("Error: --check-sex/--impute-sex y-only requires Y chromosome data.\n");
      goto sexcheck_ret_INVALID_CMDLINE;
    } else {
      LOGERRPRINTF("Warning: No Y chromosome data for --%s-sex ycount.\n", do_impute? "impute" : "check");
    }
  }
  memcpy(outname_end, ".sexcheck", 10);
  if (fopen_checked(outname, "w", &outfile)) {
    goto sexcheck_ret_OPEN_FAIL;
  }
  sprintf(g_textbuf, "%%%us %%%us       PEDSEX       SNPSEX       STATUS%s%s\n", plink_maxfid, plink_maxiid, yonly? "" : "            F", check_y? "   YCOUNT" : "");
  fprintf(outfile, g_textbuf, "FID", "IID");
  sample_uidx = 0;
  if (do_impute) {
    bitvec_andnot(sample_exclude, unfiltered_sample_ctl, sex_nm);
  }
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, sample_uidx++) {
    next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    fid_ptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    iid_ptr = (char*)memchr(fid_ptr, '\t', max_sample_id_len);
    wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(iid_ptr - fid_ptr), fid_ptr, g_textbuf);
    *wptr++ = ' ';
    wptr = fw_strcpy(plink_maxiid, &(iid_ptr[1]), wptr);
    if (!IS_SET(sex_nm, sample_uidx)) {
      orig_sex_code = 0;
    } else {
      orig_sex_code = 2 - IS_SET(sex_male, sample_uidx);
    }
    wptr = memseta(wptr, 32, 12);
    *wptr++ = '0' + orig_sex_code;
    wptr = memseta(wptr, 32, 12);
    imputed_sex_code = 0;
    if (!yonly) {
      ulii = x_variant_ct - missing_cts[sample_idx];
      if (ulii) {
	dee = nei_sum - nei_offsets[sample_idx];
	dtot = (double)((intptr_t)ulii) - dee;
	dff = (dtot - ((double)((int32_t)(het_cts[sample_idx])))) / dtot;
	if (dff > check_sex_mthresh) {
	  if ((!check_y) || (ymiss_cts[sample_idx] + min_m_yobs <= ytotal)) {
	    imputed_sex_code = 1;
	  }
	} else if ((dff < check_sex_fthresh) && ((!check_y) || (ymiss_cts[sample_idx] + max_f_yobs >= ytotal))) {
	  imputed_sex_code = 2;
	}
        *wptr++ = '0' + imputed_sex_code;
        if (orig_sex_code && (orig_sex_code == imputed_sex_code)) {
	  wptr = memcpya(wptr, "           OK ", 14);
        } else {
	  wptr = memcpya(wptr, "      PROBLEM ", 14);
	  problem_ct++;
        }
        wptr = dtoa_g_wxp4(dff, 12, wptr);
      } else {
        wptr = memcpya(wptr, "0      PROBLEM          nan", 27);
        problem_ct++;
      }
      if (check_y) {
	*wptr++ = ' ';
	wptr = uint32toa_w8(ytotal - ymiss_cts[sample_idx], wptr);
      }
    } else {
      if (ymiss_cts[sample_idx] + min_m_yobs <= ytotal) {
        imputed_sex_code = 1;
      } else if (ymiss_cts[sample_idx] + max_f_yobs >= ytotal) {
        imputed_sex_code = 2;
      }
      *wptr++ = '0' + imputed_sex_code;
      if (orig_sex_code && (orig_sex_code == imputed_sex_code)) {
	wptr = memcpya(wptr, "           OK ", 14);
      } else {
	wptr = memcpya(wptr, "      PROBLEM ", 14);
	problem_ct++;
      }
      wptr = uint32toa_w8(ytotal - ymiss_cts[sample_idx], wptr);
    }
    *wptr++ = '\n';
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
      goto sexcheck_ret_WRITE_FAIL;
    }
    if (do_impute) {
      if (imputed_sex_code) {
	SET_BIT(sample_uidx, sex_nm);
	if (imputed_sex_code == 1) {
	  SET_BIT(sample_uidx, sex_male);
	} else {
	  CLEAR_BIT(sample_uidx, sex_male);
	}
      } else {
	CLEAR_BIT(sample_uidx, sex_nm);
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto sexcheck_ret_WRITE_FAIL;
  }
  if (do_impute) {
    bitvec_and(sex_nm, unfiltered_sample_ctl, sex_male);
    gender_unk_ct = sample_ct - popcount_longs(sex_nm, unfiltered_sample_ctl);
    if (!gender_unk_ct) {
      LOGPREPRINTFWW("--impute-sex: %" PRIuPTR " Xchr and %" PRIuPTR " Ychr variant(s) scanned, all sexes imputed. Report written to %s .\n", x_variant_ct, ytotal, outname);
    } else {
      LOGPREPRINTFWW("--impute-sex: %" PRIuPTR " Xchr and %" PRIuPTR " Ychr variant(s) scanned, %" PRIuPTR "/%" PRIuPTR " sex%s imputed. Report written to %s .\n", x_variant_ct, ytotal, (sample_ct - gender_unk_ct), sample_ct, (sample_ct == 1)? "" : "es", outname);
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
    logerrprint("Error: --check-sex/--impute-sex requires at least one polymorphic X chromosome\nlocus.\n");
  sexcheck_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  bigstack_reset(bigstack_mark);
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
  if (fopen_checked(outname, "w", &outfile)) {
    goto write_snplist_ret_OPEN_FAIL;
  }
  if (!list_23_indels) {
    while (markers_done < marker_ct) {
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
    }
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

int32_t write_var_ranges(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t write_var_range_ct) {
  FILE* outfile = nullptr;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  int32_t retval = 0;
  uintptr_t new_marker_idx;
  uint32_t block_idx;
  if (write_var_range_ct > marker_ct) {
    logerrprint("Error: --write-var-ranges block count exceeds the total number of variants.\n");
    goto write_var_ranges_ret_INVALID_CMDLINE;
  }
  memcpy(outname_end, ".var.ranges", 12);
  if (fopen_checked(outname, "w", &outfile)) {
    goto write_var_ranges_ret_OPEN_FAIL;
  }
  if (fputs_checked("FIRST\tLAST\n", outfile)) {
    goto write_var_ranges_ret_WRITE_FAIL;
  }
  for (block_idx = 1; block_idx <= write_var_range_ct; block_idx++) {
    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
    fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
    putc_unlocked('\t', outfile);
    new_marker_idx = (block_idx * ((uint64_t)marker_ct)) / write_var_range_ct;
    if (new_marker_idx > marker_idx + 1) {
      marker_uidx = jump_forward_unset_unsafe(marker_exclude, marker_uidx + 1, new_marker_idx - marker_idx - 1);
    }
    fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
    putc_unlocked('\n', outfile);
    marker_uidx++;
    marker_idx = new_marker_idx;
  }
  if (fclose_null(&outfile)) {
    goto write_var_ranges_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--write-var-ranges: Block boundaries written to %s .\n", outname);
  while (0) {
  write_var_ranges_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_var_ranges_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  write_var_ranges_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t list_duplicate_vars(char* outname, char* outname_end, uint32_t dupvar_modifier, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, Chrom_info* chrom_info_ptr, char** marker_allele_ptrs) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uint32_t* uidx_list_end = (uint32_t*)g_bigstack_end;
  uint32_t* group_list_start = (uint32_t*)g_bigstack_base;
  uint32_t* group_write = group_list_start;
  uint32_t require_same_ref = dupvar_modifier & DUPVAR_REF;
  uint32_t ids_only = dupvar_modifier & DUPVAR_IDS_ONLY;
  uint32_t suppress_first_id = dupvar_modifier & DUPVAR_SUPPRESS_FIRST;
  int32_t retval = 0;
  uintptr_t max_batch_size;
  uintptr_t marker_uidx2;
  uintptr_t* uniqueness_check_bitfield;
  uint32_t* uidx_list; // grows down from workspace end
  uint32_t* group_idxs;
  uint32_t* group_sizes;
  uint32_t* group_rep_uidx_list;
  uint32_t* read_uiptr;
  uint32_t* write_uiptr;
  uint32_t* reported_id_htable;
  char* a1ptr;
  char* a2ptr;
  char* wptr_start;
  char* wptr;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t marker_uidx;
  uint32_t marker_idx;
  uint32_t last_pos;
  uint32_t last_uidx;
  uint32_t next_pos;
  uint32_t cur_batch_size;
  uint32_t item_idx;
  uint32_t group_idx;
  uint32_t group_ct;
  uint32_t duplicate_group_ct;
  uint32_t write_offset;
  uint32_t htable_entry_ct;
  uint32_t reported_id_htable_size;
  uint32_t uniqueness_check_ct;
  uint32_t uii;
  uint32_t ujj;
  uidx_list_end--;
  memcpy(outname_end, ".dupvar", 8);
  if (fopen_checked(outname, "w", &outfile)) {
    goto list_duplicate_vars_ret_OPEN_FAIL;
  }
  if (!ids_only) {
    if (fputs_checked(require_same_ref? "CHR\tPOS\tREF\tALT\tIDS\n" : "CHR\tPOS\tALLELES\tIDS\n", outfile)) {
      goto list_duplicate_vars_ret_WRITE_FAIL;
    }
  }
  max_batch_size = bigstack_left() / (5 * sizeof(int32_t));
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
    if (marker_uidx == chrom_end) {
      continue;
    }
    wptr_start = chrom_name_write(chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], g_textbuf);
    *wptr_start++ = '\t';
    last_pos = marker_pos[marker_uidx];
    while (1) {
      last_uidx = marker_uidx;
      marker_uidx++;
      next_unset_ck(marker_exclude, chrom_end, &marker_uidx);
      if (marker_uidx == chrom_end) {
	break;
      }
      next_pos = marker_pos[marker_uidx];
      if (next_pos == last_pos) {
	uidx_list = uidx_list_end;
	*uidx_list = last_uidx;
	cur_batch_size = 1;
        do {
	  if (++cur_batch_size > max_batch_size) {
	    goto list_duplicate_vars_ret_NOMEM;
	  }
	  *(--uidx_list) = marker_uidx;
	  marker_uidx++;
	  next_unset_ck(marker_exclude, chrom_end, &marker_uidx);
	  if (marker_uidx == chrom_end) {
	    break;
	  }
          next_pos = marker_pos[marker_uidx];
	} while (next_pos == last_pos);
	// make uidx_list increasing- instead of decreasing-order
	for (ujj = 0; ujj < cur_batch_size / 2; ujj++) {
	  uii = uidx_list[ujj];
	  uidx_list[ujj] = uidx_list[cur_batch_size - 1 - ujj];
	  uidx_list[cur_batch_size - 1 - ujj] = uii;
	}
	group_idxs = &(uidx_list[-((int32_t)cur_batch_size)]);
	group_rep_uidx_list = &(group_idxs[-((int32_t)cur_batch_size)]);
	group_sizes = &(group_rep_uidx_list[-((int32_t)cur_batch_size)]);
	group_ct = 0;
        for (item_idx = 0; item_idx < cur_batch_size; item_idx++) {
	  uii = uidx_list[item_idx] * 2;
	  a1ptr = marker_allele_ptrs[uii];
	  a2ptr = marker_allele_ptrs[uii + 1];
	  if (!require_same_ref) {
	    // there should be too few duplicates for a hash table to make
	    // sense here
	    for (group_idx = 0; group_idx < group_ct; group_idx++) {
	      uii = group_rep_uidx_list[group_idx] * 2;
	      if (((!strcmp(a1ptr, marker_allele_ptrs[uii])) && (!strcmp(a2ptr, marker_allele_ptrs[uii + 1]))) || ((!strcmp(a2ptr, marker_allele_ptrs[uii])) && (!strcmp(a1ptr, marker_allele_ptrs[uii + 1])))) {
		break;
	      }
	    }
	  } else {
	    for (group_idx = 0; group_idx < group_ct; group_idx++) {
	      uii = group_rep_uidx_list[group_idx] * 2;
	      if ((!strcmp(a1ptr, marker_allele_ptrs[uii])) && (!strcmp(a2ptr, marker_allele_ptrs[uii + 1]))) {
		break;
	      }
	    }
	  }
	  if (group_idx == group_ct) {
	    group_rep_uidx_list[group_ct] = uidx_list[item_idx];
	    group_sizes[group_ct++] = 1;
	  } else {
	    group_sizes[group_idx] += 1;
	  }
	  group_idxs[item_idx] = group_idx;
	}
	if (group_ct < cur_batch_size) {
	  // now identify equivalence classes
	  write_offset = 0;
	  duplicate_group_ct = 0;
	  for (group_idx = 0; group_idx < group_ct; group_idx++) {
	    if (group_sizes[group_idx] > 1) {
	      uii = group_sizes[group_idx];
	      group_sizes[group_idx] = write_offset;
	      write_offset += uii;
	      duplicate_group_ct++;
	    } else {
	      group_sizes[group_idx] = 0xffffffffU;
	    }
	  }
	  for (item_idx = 0; item_idx < cur_batch_size; item_idx++) {
	    uii = group_idxs[item_idx];
	    ujj = group_sizes[uii];
	    if (ujj != 0xffffffffU) {
	      // actually a write offset
	      group_write[ujj] = uidx_list[item_idx];
	      group_sizes[uii] += 1;
	    }
	  }
	  for (group_idx = 0; group_idx < group_ct; group_idx++) {
	    if (group_sizes[group_idx] != 0xffffffffU) {
	      // set high bit of last member of each group
	      group_write[group_sizes[group_idx] - 1] |= 0x80000000U;
	    }
	  }
	  if (suppress_first_id) {
	    read_uiptr = group_write;
	    write_uiptr = group_write;
	    for (group_idx = 0; group_idx < duplicate_group_ct; group_idx++) {
	      read_uiptr++;
	      do {
		uii = *read_uiptr++;
		*write_uiptr++ = uii;
	      } while (!(uii & 0x80000000U));
	    }
	    write_offset -= duplicate_group_ct;
	  }
	  if (!ids_only) {
	    read_uiptr = group_write;
	    for (group_idx = 0; group_idx < duplicate_group_ct; group_idx++) {
	      wptr = uint32toa_x(last_pos, '\t', wptr_start);
	      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
		goto list_duplicate_vars_ret_WRITE_FAIL;
	      }
	      uii = *read_uiptr;
	      ujj = uii & 0x7fffffff;
	      a1ptr = marker_allele_ptrs[ujj * 2];
	      a2ptr = marker_allele_ptrs[ujj * 2 + 1];
	      if (!require_same_ref) {
		// use ASCII order
		if (strcmp(a1ptr, a2ptr) < 0) {
		  fputs(a1ptr, outfile);
		  putc_unlocked(',', outfile);
		  fputs(a2ptr, outfile);
		} else {
		  fputs(a2ptr, outfile);
		  putc_unlocked(',', outfile);
		  fputs(a1ptr, outfile);
		}
	      } else {
		fputs(a2ptr, outfile);
		putc_unlocked('\t', outfile);
		fputs(a1ptr, outfile);
	      }
	      putc_unlocked('\t', outfile);
	      while (1) {
		fputs(&(marker_ids[ujj * max_marker_id_len]), outfile);
		read_uiptr++;
		if (uii & 0x80000000U) {
		  break;
		}
		putc_unlocked(' ', outfile);
		uii = *read_uiptr;
		ujj = uii & 0x7fffffff;
	      }
	      if (putc_checked('\n', outfile)) {
		goto list_duplicate_vars_ret_WRITE_FAIL;
	      }
	    }
	  } else {
	    group_write = &(group_write[write_offset]);
	    // technically an underestimate, but we don't care
	    max_batch_size = ((uintptr_t)(uidx_list_end - group_write)) / 5;
	  }
	}
	if (marker_uidx == chrom_end) {
	  break;
	}
      }
      last_pos = next_pos;
    }
  }
  if (ids_only && (group_write != group_list_start)) {
    htable_entry_ct = (uintptr_t)(group_write - group_list_start);
    bigstack_alloc_ui(htable_entry_ct, &group_list_start);
    if (bigstack_alloc_ul(unfiltered_marker_ctl, &uniqueness_check_bitfield)) {
      goto list_duplicate_vars_ret_NOMEM;
    }
    fill_all_bits(unfiltered_marker_ct, uniqueness_check_bitfield);
    for (uii = 0; uii < htable_entry_ct; uii++) {
      clear_bit((group_list_start[uii] & 0x7fffffff), uniqueness_check_bitfield);
    }
    retval = alloc_and_populate_id_htable(unfiltered_marker_ct, uniqueness_check_bitfield, htable_entry_ct, marker_ids, max_marker_id_len, 0, &reported_id_htable_size, &reported_id_htable);
    if (retval) {
      goto list_duplicate_vars_ret_1;
    }
    bitarr_invert(unfiltered_marker_ct, uniqueness_check_bitfield);
    bitvec_or(marker_exclude, unfiltered_marker_ctl, uniqueness_check_bitfield);
    uniqueness_check_ct = marker_ct - htable_entry_ct;
    for (marker_uidx2 = 0, marker_idx = 0; marker_idx < uniqueness_check_ct; marker_uidx2++, marker_idx++) {
      next_unset_ul_unsafe_ck(uniqueness_check_bitfield, &marker_uidx2);
      wptr = &(marker_ids[marker_uidx2 * max_marker_id_len]);
      if (id_htable_find(wptr, strlen(wptr), reported_id_htable, reported_id_htable_size, marker_ids, max_marker_id_len) != 0xffffffffU) {
        LOGERRPRINTFWW("Error: Duplicate ID '%s'.\n", wptr);
	goto list_duplicate_vars_ret_INVALID_FORMAT;
      }
    }
    read_uiptr = group_list_start;
    do {
      uii = *read_uiptr++;
      fputs(&(marker_ids[(uii & 0x7fffffff) * max_marker_id_len]), outfile);
      if (uii & 0x80000000U) {
	if (putc_checked('\n', outfile)) {
	  goto list_duplicate_vars_ret_WRITE_FAIL;
	}
      } else {
	putc_unlocked(' ', outfile);
      }
    } while (read_uiptr != group_write);
  }
  if (fclose_null(&outfile)) {
    goto list_duplicate_vars_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--list-duplicate-vars report written to %s .\n", outname);
  while (0) {
  list_duplicate_vars_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  list_duplicate_vars_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  list_duplicate_vars_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  list_duplicate_vars_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 list_duplicate_vars_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t het_report(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t output_gz, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uintptr_t* founder_info, Chrom_info* chrom_info_ptr, double* set_allele_freqs) {
  // Same F coefficient computation as sexcheck().
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t* loadbuf_f = nullptr;
  uintptr_t* founder_vec11 = nullptr;
  char* pzwritep = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  uintptr_t final_mask = get_final_mask(sample_ct);
  uintptr_t founder_ct = 0;
  uintptr_t monomorphic_ct = 0;
  uintptr_t marker_uidx = 0;
  double nei_sum = 0.0;
  uint32_t chrom_fo_idx = 0xffffffffU; // deliberate overflow
  uint32_t chrom_end = 0;
  int32_t mt_code = chrom_info_ptr->xymt_codes[MT_OFFSET];
  int32_t retval = 0;
  Pigz_state ps;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* lptr;
  uint32_t* het_cts;
  uint32_t* missing_cts;
  double* nei_offsets;
  unsigned char* overflow_buf;
  char* fid_ptr;
  char* iid_ptr;
  double dpp;
  double dtot;
  double cur_nei;
  double dff;
  double dee;
  uintptr_t marker_idx;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t cur_missing_ct;
  uintptr_t allele_obs_ct;
  uintptr_t f_missing_ct;
  uintptr_t f_allele_obs_ct;
  uintptr_t cur_word;
  uintptr_t ulii;
  uint32_t obs_ct;
  pzwrite_init_null(&ps);
  if (is_set(chrom_info_ptr->haploid_mask, 0)) {
    logerrprint("Error: --het cannot be used on haploid genomes.\n");
    goto het_report_ret_INVALID_CMDLINE;
  }
  if (bigstack_alloc_uc(PIGZ_BLOCK_SIZE + MAXLINELEN, &overflow_buf) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(sample_ctl2, &loadbuf) ||
      bigstack_calloc_ui(sample_ct, &het_cts) ||
      bigstack_calloc_ui(sample_ct, &missing_cts) ||
      bigstack_calloc_d(sample_ct, &nei_offsets)) {
    goto het_report_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  loadbuf[sample_ctl2 - 1] = 0;
  marker_ct -= count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1, 1);
  if (!marker_ct) {
    goto het_report_ret_INVALID_CMDLINE;
  }
  if (founder_info) {
    // founder_info passed iff we're including the small-sample correction, in
    // which case we may need to mask out nonfounders and determine
    // missing/allele counts in that subset.
    founder_ct = popcount_longs(founder_info, unfiltered_sample_ctl);
    if (founder_ct < sample_ct) {
      if (bigstack_alloc_ul(sample_ctl2, &founder_vec11) ||
          bigstack_alloc_ul(sample_ctl2, &loadbuf_f)) {
	goto het_report_ret_NOMEM;
      }
      quaterarr_collapse_init_exclude(founder_info, unfiltered_sample_ct, sample_exclude, sample_ct, founder_vec11);
      lptr = founder_vec11;
      for (ulii = 0; ulii < sample_ctl2; ulii++) {
	*lptr = (*lptr) * 3;
	lptr++;
      }
    } else {
      loadbuf_f = loadbuf;
    }
  }
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
        goto het_report_ret_READ_FAIL;
      }
    }
    if (marker_uidx >= chrom_end) {
      do {
        uint32_t chrom_idx;
	do {
	  chrom_fo_idx++;
          // bugfix (16 Jun 2020): forgot to separately exclude chrM here.
          // Fortunately, this frequently didn't matter, since chrM is usually
          // sorted last: in that case, if no alternate contigs are present,
          // the marker_ct loop termination condition prevents chrM from being
          // included.  But chrM is positioned first in some files, in which
          // case this fix matters.
          chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	} while (is_set(chrom_info_ptr->haploid_mask, chrom_idx) || (chrom_idx == (uint32_t)mt_code));
	chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
	marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
      } while (marker_uidx >= chrom_end);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
        goto het_report_ret_READ_FAIL;
      }
    }
    if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, 0, bedfile, loadbuf_raw, loadbuf)) {
      goto het_report_ret_READ_FAIL;
    }
    cur_missing_ct = count_01(loadbuf, sample_ctl2);
    if (!loadbuf_f) {
      allele_obs_ct = 2 * (sample_ct - cur_missing_ct);
      dpp = set_allele_freqs[marker_uidx];
      if ((!allele_obs_ct) || (dpp < 1e-8) || (dpp > (1 - 1e-8))) {
        monomorphic_ct++;
        continue;
      }
      cur_nei = 1.0 - 2 * dpp * (1 - dpp);
    } else {
      if (founder_vec11) {
        memcpy(loadbuf_f, loadbuf, sample_ctl2 * sizeof(intptr_t));
        bitvec_and(founder_vec11, sample_ctl2, loadbuf_f);
        f_missing_ct = count_01(loadbuf_f, sample_ctl2);
      } else {
	f_missing_ct = cur_missing_ct;
      }
      f_allele_obs_ct = 2 * (founder_ct - f_missing_ct);
      ulii = popcount_longs(loadbuf_f, sample_ctl2) - f_missing_ct;
      if ((!f_allele_obs_ct) || (!ulii) || (ulii == f_allele_obs_ct)) {
	monomorphic_ct++;
	continue;
      }
      dee = (double)((intptr_t)ulii);
      dff = (double)((intptr_t)f_allele_obs_ct);
      dpp = dee / dff;
      cur_nei = 1.0 - 2 * (1 - dpp) * dee / (dff - 1);
    }
    if (cur_missing_ct) {
      // iterate through missing calls
      lptr = loadbuf;
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx += BITCT2) {
        cur_word = *lptr++;
        cur_word = (~(cur_word >> 1)) & cur_word & FIVEMASK;
        while (cur_word) {
          ulii = sample_idx + CTZLU(cur_word) / 2;
          missing_cts[ulii] += 1;
          nei_offsets[ulii] += cur_nei;
          cur_word &= cur_word - 1;
	}
      }
    }
    lptr = loadbuf;
    // iterate through heterozygous calls
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx += BITCT2) {
      cur_word = *lptr++;
      cur_word = (cur_word >> 1) & (~cur_word) & FIVEMASK;
      while (cur_word) {
        het_cts[sample_idx + CTZLU(cur_word) / 2] += 1;
        cur_word &= cur_word - 1;
      }
    }
    nei_sum += cur_nei;
  }
  marker_ct -= monomorphic_ct;
  if (!marker_ct) {
    goto het_report_ret_INVALID_CMDLINE;
  }
  memcpy(outname_end, output_gz? ".het.gz" : ".het", output_gz? 8 : 5);
  if (flex_pzwrite_init(output_gz, outname, overflow_buf, 0, &ps)) {
    goto het_report_ret_OPEN_FAIL;
  }
  pzwritep = (char*)overflow_buf;
  sprintf(g_textbuf, "%%%us %%%us       O(HOM)       E(HOM)        N(NM)            F\n", plink_maxfid, plink_maxiid);
  pzwritep += sprintf(pzwritep, g_textbuf, "FID", "IID");
  sample_uidx = 0;
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, sample_uidx++) {
    next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    fid_ptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    iid_ptr = (char*)memchr(fid_ptr, '\t', max_sample_id_len);
    pzwritep = fw_strcpyn(plink_maxfid, (uintptr_t)(iid_ptr - fid_ptr), fid_ptr, pzwritep);
    *pzwritep++ = ' ';
    pzwritep = fw_strcpy(plink_maxiid, &(iid_ptr[1]), pzwritep);
    pzwritep = memseta(pzwritep, 32, 3);
    obs_ct = marker_ct - missing_cts[sample_idx];
    if (obs_ct) {
      pzwritep = uint32toa_w10x(obs_ct - het_cts[sample_idx], ' ', pzwritep);
      dee = nei_sum - nei_offsets[sample_idx];
      pzwritep = dtoa_g_wxp4(dee, 12, pzwritep);
      pzwritep = memseta(pzwritep, 32, 3);
      pzwritep = uint32toa_w10x(obs_ct, ' ', pzwritep);
      dtot = (double)((int32_t)obs_ct) - dee;
      dff = (dtot - ((double)((int32_t)(het_cts[sample_idx])))) / dtot;
      pzwritep = dtoa_g_wxp4(dff, 12, pzwritep);
    } else {
      pzwritep = memcpya(pzwritep, "         0            0            0          nan", 49);
    }
    append_binary_eoln(&pzwritep);
    if (flex_pzwrite(&ps, &pzwritep)) {
      goto het_report_ret_WRITE_FAIL;
    }
  }
  if (flex_pzwrite_close_null(&ps, pzwritep)) {
    goto het_report_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--het%s: %" PRIuPTR " variant%s scanned, report written to %s .\n", loadbuf_f? " small-sample" : "", marker_ct, (marker_ct == 1)? "" : "s", outname);
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
    logerrprint("Error: --het requires at least one polymorphic autosomal marker.\n");
    retval = RET_INVALID_CMDLINE;
    break;
  }
  bigstack_reset(bigstack_mark);
  flex_pzwrite_close_cond(&ps, pzwritep);
  return retval;
}

int32_t fst_report(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts) {
  // Math based on VCFtools variant_file::output_weir_and_cockerham_fst();
  // frequency counting logic similar to cmh_assoc().
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* wptr_start = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  uint32_t chrom_fo_idx = 0xffffffffU;
  uint32_t chrom_end = 0;
  uint32_t chrom_idx = 0;
  uint32_t skipped_marker_ct = 0;
  uint32_t pct = 0;
  int32_t mt_code = chrom_info_ptr->xymt_codes[MT_OFFSET];
  int32_t retval = 0;
  uintptr_t* loadbuf;
  uintptr_t* cluster_mask;
  uintptr_t* ulptr;
  uintptr_t* ulptr2;
  char* wptr;
  uint32_t* sample_to_cluster;
  uint32_t* cluster_sizes;
  uint32_t* cluster_geno_cts;
  uint32_t* uiptr;
  uint32_t* uiptr2;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t cur_word;
  double cluster_ctd;
  double cluster_ct_recip;
  double cluster_ctm1_recip;
  double one_minus_cluster_ct_recip;
  double n_sum;
  double n_sum_recip;
  double n_ssq;
  double nbar;
  double nbar_recip;
  double pbar_a1;
  double hbar_a1;
  double ssqr_a1;
  double nc;
  double a_a1;
  double b_a1;
  double c_a1;
  double dxx;
  double dyy;
  uint32_t sample_uidx_base;
  uint32_t sample_uidx;
  uint32_t sample_idx;
  uint32_t cur_sample_ct;
  uint32_t seek_flag;
  uint32_t cluster_idx;
  uint32_t loop_end;
  uint32_t uii;
  uint32_t ujj;
  marker_ct -= count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1, 1);
  if (!marker_ct) {
    logerrprint("Error: No autosomal diploid variants for --fst.\n");
    goto fst_report_ret_INVALID_CMDLINE;
  }
  if (pheno_c) {
    cluster_ct = 2;
  }
  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf) ||
      bigstack_calloc_ul(unfiltered_sample_ctl2, &cluster_mask) ||
      bigstack_alloc_ui(unfiltered_sample_ct, &sample_to_cluster) ||
      bigstack_alloc_ui(cluster_ct, &cluster_sizes)) {
    goto fst_report_ret_NOMEM;
  }
  if (pheno_c) {
    cur_sample_ct = popcount_longs(pheno_nm, unfiltered_sample_ctl);
    uii = popcount_longs(pheno_c, unfiltered_sample_ctl);
    if ((!uii) || (cur_sample_ct == uii)) {
      logerrprint("Error: --fst case-control requires at least one case and one control.\n");
      goto fst_report_ret_INVALID_CMDLINE;
    }
    cluster_sizes[0] = cur_sample_ct - uii;
    cluster_sizes[1] = uii;
    for (sample_uidx = 0, sample_idx = 0; sample_idx < cur_sample_ct; sample_uidx++, sample_idx++) {
      next_set_unsafe_ck(pheno_nm, &sample_uidx);
      cluster_mask[sample_uidx / BITCT2] |= (3 * ONELU) << (2 * (sample_uidx % BITCT2));
      sample_to_cluster[sample_uidx] = is_set(pheno_c, sample_uidx);
    }
  } else {
    uii = 0; // cluster index after excluding empty ones
    uiptr = cluster_map;
    for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
      ujj = 0;
      uiptr2 = &(cluster_map[cluster_starts[cluster_idx + 1]]);
      for (; uiptr < uiptr2; uiptr++) {
        sample_uidx = *uiptr;
	if (!is_set(sample_exclude, sample_uidx)) {
	  cluster_mask[sample_uidx / BITCT2] |= (3 * ONELU) << (2 * (sample_uidx % BITCT2));
          sample_to_cluster[sample_uidx] = uii;
          ujj++;
	}
      }
      if (ujj) {
	cluster_sizes[uii++] = ujj;
      }
    }
    cluster_ct = uii;
    if (cluster_ct < 2) {
      logerrprint("Error: --fst requires at least two nonempty clusters.\n");
      goto fst_report_ret_INVALID_CMDLINE;
    }
    bigstack_shrink_top(cluster_sizes, cluster_ct * sizeof(int32_t));
  }
  cluster_ctd = (double)((intptr_t)cluster_ct);
  cluster_ct_recip = 1.0 / cluster_ctd;
  cluster_ctm1_recip = 1.0 / (cluster_ctd - 1.0);
  one_minus_cluster_ct_recip = 1.0 - cluster_ct_recip;
  if (bigstack_alloc_ui(cluster_ct * 3, &cluster_geno_cts)) {
    goto fst_report_ret_NOMEM;
  }
  memcpy(outname_end, ".fst", 5);
  if (fopen_checked(outname, "w", &outfile)) {
    goto fst_report_ret_OPEN_FAIL;
  }
  if (fputs_checked("CHR\tSNP\tPOS\tNMISS\tFST\n", outfile)) {
    goto fst_report_ret_WRITE_FAIL;
  }
  LOGPRINTFWW5("Writing --fst report (%" PRIuPTR " populations) to %s ... ", cluster_ct, outname);
  fputs("0%", stdout);
  fflush(stdout);
  seek_flag = 1;
  loop_end = marker_ct / 100;
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      seek_flag = 1;
    }
    if (marker_uidx >= chrom_end) {
      while (1) {
	do {
          chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
	} while (marker_uidx >= chrom_end);
        chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
        if ((!IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) && (chrom_idx != (uint32_t)mt_code)) {
	  break;
	}
	seek_flag = 1;
	marker_uidx = next_unset_unsafe(marker_exclude, chrom_end);
      }
      wptr_start = chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf);
      *wptr_start++ = '\t';
    }
    if (seek_flag) {
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
        goto fst_report_ret_READ_FAIL;
      }
      seek_flag = 0;
    }
    if (load_raw(unfiltered_sample_ct4, bedfile, loadbuf)) {
      goto fst_report_ret_READ_FAIL;
    }
    fill_uint_zero(cluster_ct * 3, cluster_geno_cts);
    ulptr = loadbuf;
    ulptr2 = cluster_mask;
    for (sample_uidx_base = 0; sample_uidx_base < unfiltered_sample_ct; sample_uidx_base += BITCT2) {
      cur_word = (~(*ulptr++)) & (*ulptr2++);
      while (cur_word) {
        uii = CTZLU(cur_word) & (BITCT - 2);
        ujj = (cur_word >> uii) & 3;
        sample_uidx = sample_uidx_base + (uii / 2);
        cluster_idx = sample_to_cluster[sample_uidx];
        cluster_geno_cts[cluster_idx * 3 + ujj - 1] += 1;
        cur_word &= ~((3 * ONELU) << uii);
      }
    }
    cur_sample_ct = 0;
    n_ssq = 0.0;
    pbar_a1 = 0.0;
    hbar_a1 = 0.0;
    for (cluster_idx = 0, uiptr = cluster_geno_cts; cluster_idx < cluster_ct; cluster_idx++, uiptr = &(uiptr[3])) {
      uii = cluster_sizes[cluster_idx] - uiptr[1];
      cur_sample_ct += uii;
      dxx = (double)((int32_t)uii);
      n_ssq += dxx * dxx;
      pbar_a1 += (double)(uiptr[0] + 2 * uiptr[2]);
      hbar_a1 += uiptr[0];
    }
    n_sum = (double)((int32_t)cur_sample_ct);
    n_sum_recip = 1.0 / n_sum;
    nbar = n_sum * cluster_ct_recip;
    nbar_recip = n_sum_recip * cluster_ctd;
    // no need to iterate over alleles when there are only two of them, but
    // that changes when we add support for the general case in the future
    pbar_a1 *= 0.5 * n_sum_recip;
    hbar_a1 *= n_sum_recip;
    ssqr_a1 = 0.0;
    for (cluster_idx = 0, uiptr = cluster_geno_cts; cluster_idx < cluster_ct; cluster_idx++, uiptr = &(uiptr[3])) {
      dxx = (double)((int32_t)(cluster_sizes[cluster_idx] - uiptr[1]));
      uii = uiptr[0] + 2 * uiptr[2]; // A1 obs ct
      dyy = ((double)uii) / (2 * dxx) - pbar_a1;
      ssqr_a1 += dxx * dyy * dyy;
    }
    ssqr_a1 *= nbar_recip * cluster_ctm1_recip;
    nc = (n_sum - n_ssq * n_sum_recip) * cluster_ctm1_recip;
    dxx = pbar_a1 * (1.0 - pbar_a1) - ssqr_a1 * one_minus_cluster_ct_recip;
    dyy = 1.0 / (nbar - 1.0);
    a_a1 = (ssqr_a1 - (dxx - hbar_a1 * 0.25) * dyy) * nbar / nc;
    b_a1 = (dxx - hbar_a1 * (2.0 * nbar - 1.0) * nbar_recip * 0.25) * nbar * dyy;
    c_a1 = hbar_a1 * 0.5;
    dxx = a_a1 + b_a1 + c_a1;
    dyy = a_a1 / dxx; // fst estimate
    if (dyy != dyy) {
      skipped_marker_ct++;
    } else {
      sum1 += a_a1;
      sum2 += dxx;
      sum3 += dyy;
    }
    wptr = strcpyax(wptr_start, &(marker_ids[marker_uidx * max_marker_id_len]), '\t');
    wptr = uint32toa_x(marker_pos[marker_uidx], '\t', wptr);
    wptr = uint32toa_x(cur_sample_ct, '\t', wptr);
    wptr = dtoa_gx(dyy, '\n', wptr);
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
      goto fst_report_ret_WRITE_FAIL;
    }
    if (marker_idx >= loop_end) {
      if (marker_idx < marker_ct) {
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
        pct = (marker_idx * 100LLU) / marker_ct;
        printf("\b\b%u%%", pct);
        fflush(stdout);
        loop_end = (((uint64_t)pct + 1LLU) * marker_ct) / 100;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto fst_report_ret_WRITE_FAIL;
  }
  if (pct >= 10) {
    putc_unlocked('\b', stdout);
  }
  fputs("\b\b", stdout);
  logprint("done.\n");
  marker_ct -= skipped_marker_ct;
  if (skipped_marker_ct) {
    LOGPRINTF("%" PRIuPTR " marker%s (%u excluded).\n", marker_ct, (marker_ct == 1)? " with a valid Fst estimate" : "s with valid Fst estimates", skipped_marker_ct);
  } else {
    LOGPRINTF("%" PRIuPTR " marker%s.\n", marker_ct, (marker_ct == 1)? " with a valid Fst estimate" : "s with valid Fst estimates");
  }
  LOGPRINTF("Mean Fst estimate: %g\n", sum3 / ((double)((intptr_t)marker_ct)));
  LOGPRINTF("Weighted Fst estimate: %g\n", sum1 / sum2);

  while (0) {
  fst_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  fst_report_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  fst_report_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  fst_report_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  fst_report_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t score_report(Score_info* sc_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, double* set_allele_freqs, uintptr_t sample_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t hh_exists, Chrom_info* chrom_info_ptr, char* outname, char* outname_end) {
  // Note that there is a dosage-only implementation of this logic in
  // plink_dosage.c.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* infile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  uintptr_t final_mask = get_final_mask(sample_ct);
  uintptr_t miss_ct = 0;
  uint32_t miss_varid_ct = 0;
  uint32_t miss_allele_ct = 0;
  uintptr_t range_ct = 0;
  uintptr_t range_skip = 0;
  uintptr_t ulii = 0;
  char* tbuf2 = &(g_textbuf[MAXLINELEN]);
  uintptr_t* marker_exclude_main = nullptr;
  uintptr_t* sample_include2 = nullptr;
  uintptr_t* sample_male_include2 = nullptr;
  double* qrange_keys = nullptr;
  double* effect_sizes_cur = nullptr;
  double abs_score_sum = 0.0;
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
  uint32_t marker_id_htable_size = get_id_htable_size(marker_ct);
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
  uint32_t* marker_id_htable;
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
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
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
  if (bigstack_end_alloc_ui(marker_id_htable_size, &marker_id_htable)) {
    goto score_report_ret_NOMEM;
  }
  retval = populate_id_htable(unfiltered_marker_ct, marker_exclude_orig, marker_ct, marker_ids, max_marker_id_len, 0, marker_id_htable_size, marker_id_htable);
  if (retval) {
    goto score_report_ret_1;
  }
  if (bigstack_end_alloc_d(unfiltered_marker_ct, &dptr) ||
      bigstack_alloc_ul(unfiltered_marker_ctl, &marker_exclude) ||
      bigstack_calloc_ul(unfiltered_marker_ctl, &a2_effect)) {
    goto score_report_ret_NOMEM;
  }
  fill_all_bits(unfiltered_marker_ct, marker_exclude);
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto score_report_ret_NOMEM;
  }
  loadbuf_c = (char*)g_bigstack_base;
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
  memcpy(outname_end, ".nopred", 8); // bugfix, this was after the goto before
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
    marker_uidx = id_htable_find(bufptr_arr[varid_idx], strlen_se(bufptr_arr[varid_idx]), marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if (marker_uidx != 0xffffffffU) {
      bufptr_arr[allele_idx][strlen_se(bufptr_arr[allele_idx])] = '\0';
      uii = strcmp(bufptr_arr[allele_idx], marker_allele_ptrs[2 * marker_uidx]);
      if ((!uii) || (!strcmp(bufptr_arr[allele_idx], marker_allele_ptrs[2 * marker_uidx + 1]))) {
        if (scan_double(bufptr_arr[effect_idx], &dxx) || (dxx != dxx)) {
	  if (!miss_ct) {
	    if (fopen_checked(outname, "w", &outfile)) {
	      goto score_report_ret_OPEN_FAIL;
	    }
	  }
	  fputs("BADVAL\t", outfile);
	  if (fwrite_checked(bufptr_arr[varid_idx], strlen_se(bufptr_arr[varid_idx]), outfile)) {
	    goto score_report_ret_WRITE_FAIL;
	  }
	  putc_unlocked('\n', outfile);
          miss_ct++;
	} else {
	  if (!IS_SET(marker_exclude, marker_uidx)) {
            bufptr_arr[varid_idx][strlen_se(bufptr_arr[varid_idx])] = '\0';
            LOGPREPRINTFWW("Error: Duplicate variant '%s' in --score file.\n", bufptr_arr[varid_idx]);
            goto score_report_ret_INVALID_FORMAT_2;
	  }
	  dptr[marker_uidx] = dxx;
          CLEAR_BIT(marker_uidx, marker_exclude);
	  if (uii) {
	    SET_BIT(marker_uidx, a2_effect);
	  }
	}
      } else {
	if (!miss_ct) {
	  if (fopen_checked(outname, "w", &outfile)) {
	    goto score_report_ret_OPEN_FAIL;
	  }
	}
	fputs("NOALLELE\t", outfile);
	if (fwrite_checked(bufptr_arr[varid_idx], strlen_se(bufptr_arr[varid_idx]), outfile)) {
	  goto score_report_ret_WRITE_FAIL;
	}
	putc_unlocked(' ', outfile);
	fputs(bufptr_arr[allele_idx], outfile);
	fputs(" vs ", outfile);
	fputs(marker_allele_ptrs[2 * marker_uidx], outfile);
	putc_unlocked(' ', outfile);
	fputs(marker_allele_ptrs[2 * marker_uidx + 1], outfile);
	putc_unlocked('\n', outfile);
	miss_ct++;
	miss_allele_ct++;
      }
    } else {
      if (!miss_ct) {
	if (fopen_checked(outname, "w", &outfile)) {
	  goto score_report_ret_OPEN_FAIL;
	}
      }
      fputs("NOSNP\t", outfile);
      if (fwrite_checked(bufptr_arr[varid_idx], strlen_se(bufptr_arr[varid_idx]), outfile)) {
	goto score_report_ret_WRITE_FAIL;
      }
      putc_unlocked('\n', outfile);
      miss_ct++;
      miss_varid_ct++;
    }
  score_report_load_next:
    line_idx++;
    if (!fgets(loadbuf_c, loadbuf_size, infile)) {
      break;
    }
    if (!(loadbuf_c[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --score file is pathologically long.\n", line_idx);
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
    logerrprint("Error: No valid entries in --score file.\n");
    goto score_report_ret_INVALID_FORMAT;
  } else if (cur_marker_ct >= 0x40000000) {
    logerrprint("Error: --score does not support >= 2^30 variants.\n");
    goto score_report_ret_INVALID_FORMAT;
  }
  if (miss_ct) {
    LOGERRPRINTFWW("Warning: %" PRIuPTR " line%s skipped in --score file (%u due to variant ID mismatch, %u due to allele code mismatch); see %s for details.\n", miss_ct, (miss_ct == 1)? "" : "s", miss_varid_ct, miss_allele_ct, outname);
    if (fclose_null(&outfile)) {
      goto score_report_ret_WRITE_FAIL;
    }
  }
  LOGPRINTF("--score: %u valid predictor%s loaded.\n", cur_marker_ct, (cur_marker_ct == 1)? "" : "s");
  if (sc_ip->data_fname) {
    effect_sizes_cur = dptr; // not collapsed yet
    if (bigstack_end_alloc_d(unfiltered_marker_ct, &dptr) ||
        bigstack_alloc_ul(unfiltered_marker_ctl, &marker_exclude_main)) {
      goto score_report_ret_NOMEM;
    }
    fill_all_bits(unfiltered_marker_ct, marker_exclude_main);
    loadbuf_size = bigstack_left();
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size <= MAXLINELEN) {
      goto score_report_ret_NOMEM;
    }
    loadbuf_c = (char*)g_bigstack_base;
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
      marker_uidx = id_htable_find(bufptr_arr[varid_idx], strlen_se(bufptr_arr[varid_idx]), marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
      if (marker_uidx != 0xffffffffU) {
        if (!IS_SET(marker_exclude, marker_uidx)) {
	  if (scan_double(bufptr_arr[1 - varid_idx], &dxx) || (dxx != dxx) || (dxx == INFINITY) || (dxx == -INFINITY)) {
	    miss_ct++;
	  } else {
	    dptr[marker_uidx] = dxx;
	    if (!IS_SET(marker_exclude_main, marker_uidx)) {
	      bufptr_arr[varid_idx][strlen_se(bufptr_arr[varid_idx])] = '\0';
	      LOGPREPRINTFWW("Error: Duplicate variant '%s' in --q-score-range data file.\n", bufptr_arr[varid_idx]);
	      goto score_report_ret_INVALID_FORMAT_2;
	    }
            CLEAR_BIT(marker_uidx, marker_exclude_main);
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
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --q-score-range data file is pathologically long.\n", line_idx);
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
      logerrprint("Error: No valid entries in --q-score-range data file.\n");
      goto score_report_ret_INVALID_FORMAT;
    }
    qrange_keys = (double*)alloc_and_init_collapsed_arr((char*)dptr, sizeof(double), unfiltered_marker_ct, marker_exclude_main, marker_ct, 0);
    if (!qrange_keys) {
      goto score_report_ret_NOMEM;
    }
    bigstack_end_reset(effect_sizes_cur);
    effect_sizes = (double*)alloc_and_init_collapsed_arr((char*)effect_sizes_cur, sizeof(double), unfiltered_marker_ct, marker_exclude_main, marker_ct, 0);
    if (!effect_sizes) {
      goto score_report_ret_NOMEM;
    }
    if (bigstack_alloc_d(marker_ct, &effect_sizes_cur)) {
      goto score_report_ret_NOMEM;
    }
    if (miss_ct) {
      LOGERRPRINTF("Warning: %" PRIuPTR " line%s skipped in --q-score-range data file.\n", miss_ct, (miss_ct == 1)? "" : "s");
    }
    miss_ct = 0;
    if (fopen_checked(sc_ip->range_fname, "r", &infile)) {
      goto score_report_ret_OPEN_FAIL;
    }
    max_rangename_len = (FNAMESIZE - 10) - ((uintptr_t)(outname_end - outname));
    g_textbuf[MAXLINELEN - 1] = ' ';
    *outname_end = '.';
  } else {
    effect_sizes = (double*)alloc_and_init_collapsed_arr((char*)dptr, sizeof(double), unfiltered_marker_ct, marker_exclude, cur_marker_ct, 0);
    if (!effect_sizes) {
      goto score_report_ret_NOMEM;
    }
  }
  bigstack_end_reset(bigstack_end_mark);
  if (bigstack_alloc_d(sample_ct, &score_deltas) ||
      bigstack_alloc_ui(sample_ct, &miss_cts) ||
      bigstack_alloc_i(sample_ct, &named_allele_ct_deltas) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(sample_ctl2, &loadbuf)) {
    goto score_report_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  loadbuf[sample_ctl2 - 1] = 0;
  // force sample_male_include2 allocation
  if (alloc_collapsed_haploid_filters(sample_exclude, sex_male, unfiltered_sample_ct, sample_ct, hh_exists | XMHH_EXISTS, 0, &sample_include2, &sample_male_include2)) {
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
      if (!fgets(g_textbuf, MAXLINELEN, infile)) {
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
	fputs(g_logbuf, stdout);
	goto score_report_ret_1;
      }
      if (!g_textbuf[MAXLINELEN - 1]) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --q-score-range range file is pathologically long.\n", line_idx);
	goto score_report_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(g_textbuf);
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
	sprintf(g_logbuf, "Error: Excessively long range name on line %" PRIuPTR " of --q-score-range range\nfile.\n", line_idx);
	goto score_report_ret_INVALID_FORMAT_2;
      }
      bufptr_arr[0] = bufptr;
      break;
    }
    fill_all_bits(unfiltered_marker_ct, marker_exclude);
    marker_uidx = next_unset_unsafe(marker_exclude_main, 0);
    dptr = effect_sizes_cur;
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude_main, &marker_uidx);
      dxx = qrange_keys[marker_idx];
      if ((dxx >= lbound) && (dxx <= ubound)) {
	CLEAR_BIT(marker_uidx, marker_exclude);
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
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    goto score_report_ret_READ_FAIL;
  }
  fill_double_zero(sample_ct, score_deltas);
  fill_uint_zero(sample_ct, miss_cts);
  fill_int_zero(sample_ct, named_allele_ct_deltas);
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
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	goto score_report_ret_READ_FAIL;
      }
    }
    if (marker_uidx >= chrom_end) {
      do {
	chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1];
      } while (marker_uidx >= chrom_end);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      is_haploid = IS_SET(chrom_info_ptr->haploid_mask, uii);
      is_x = ((int32_t)uii == chrom_info_ptr->xymt_codes[X_OFFSET])? 1 : 0;
      is_y = ((int32_t)uii == chrom_info_ptr->xymt_codes[Y_OFFSET])? 1 : 0;
      ploidy = 2 - is_haploid;
      ploidy_d = (double)((int32_t)ploidy);
    }
    if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, loadbuf)) {
      goto score_report_ret_READ_FAIL;
    }
    if (is_haploid && hh_exists) {
      haploid_fix(hh_exists, sample_include2, sample_male_include2, sample_ct, is_x, is_y, (unsigned char*)loadbuf);
    }
    cur_effect_size = (*dptr++) * ploidy_d;
    abs_score_sum += fabs(cur_effect_size);
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
	if (uii + BITCT2 > sample_ct) {
	  ulii &= (ONELU << ((sample_ct & (BITCT2 - 1)) * 2)) - ONELU;
	}
	while (ulii) {
	  ujj = CTZLU(ulii) & (BITCT - 2);
	  ukk = (ulii >> ujj) & 3;
	  sample_idx = uii + (ujj / 2);
	  if (ukk == 1) {
	    score_deltas[sample_idx] += half_effect_size;
	    named_allele_ct_deltas[sample_idx] += delta1;
	  } else if (ukk == 3) {
	    score_deltas[sample_idx] += cur_effect_size;
	    named_allele_ct_deltas[sample_idx] += delta2;
	  } else {
	    miss_cts[sample_idx] += ploidy;
	    named_allele_ct_deltas[sample_idx] += deltam;
	    score_deltas[sample_idx] += missing_effect;
	  }
	  ulii &= ~((3 * ONELU) << ujj);
	}
	uii += BITCT2;
      } while (uii < sample_ct);
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
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	if (!(sample_idx & (BITCT2 - 1))) {
	  ulii = ~(*lbptr++);
	}
	uljj = ulii & 3;
	if (IS_SET_DBL(sample_male_include2, sample_idx)) {
	  if (uljj) {
	    if (uljj == 3) {
	      score_deltas[sample_idx] += cur_effect_size;
	      named_allele_ct_deltas[sample_idx] += delta2;
	    } else {
	      miss_cts[sample_idx] += 1;
	      named_allele_ct_deltas[sample_idx] += deltam;
	      score_deltas[sample_idx] += missing_effect;
	    }
	  }
	} else if (is_x) {
	  score_deltas[sample_idx] += female_effect_size[uljj];
	  named_allele_ct_deltas[sample_idx] += female_allele_ct_delta[uljj];
	  if (uljj == 2) {
	    miss_cts[sample_idx] += 2;
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
  if (fopen_checked(outname, "w", &outfile)) {
    goto score_report_ret_OPEN_FAIL;
  }
  sprintf(tbuf2, "%%%us %%%us  PHENO    CNT   CNT2 %s\n", plink_maxfid, plink_maxiid, report_average? "   SCORE" : "SCORESUM");
  fprintf(outfile, tbuf2, "FID", "IID");
  for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
    next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    bufptr_arr[0] = &(sample_ids[sample_uidx * max_sample_id_len]);
    uii = strlen_se(bufptr_arr[0]);
    bufptr = fw_strcpyn(plink_maxfid, uii, bufptr_arr[0], tbuf2);
    *bufptr++ = ' ';
    bufptr = fw_strcpy(plink_maxiid, &(bufptr_arr[0][uii + 1]), bufptr);
    *bufptr++ = ' ';
    if (IS_SET(pheno_nm, sample_uidx)) {
      if (pheno_c) {
	bufptr = memseta(bufptr, 32, 5);
	*bufptr++ = '1' + IS_SET(pheno_c, sample_uidx);
      } else {
	bufptr = width_force(6, bufptr, dtoa_g(pheno_d[sample_uidx], bufptr));
      }
    } else {
      bufptr = memcpya(bufptr, missing_pheno_str, missing_pheno_len);
    }
    *bufptr++ = ' ';
    ujj = 1 - IS_SET_DBL(sample_male_include2, sample_idx); // female?
    uii = obs_expected + ((int32_t)ujj) * obs_expected_female_delta - miss_cts[sample_idx];
    bufptr = uint32toa_w6x(uii, ' ', bufptr);
    if (mean_impute) {
      uii += miss_cts[sample_idx];
    }
    bufptr = uint32toa_w6x(((int32_t)named_allele_ct_expected) - ujj * named_allele_ct_female_delta + named_allele_ct_deltas[sample_idx], ' ', bufptr);
    dxx = (score_base + ((int32_t)ujj) * female_y_offset + score_deltas[sample_idx]);
    if (fabs(dxx) < abs_score_sum * RECIP_2_53) {
      dxx = 0;
    } else if (report_average) {
      dxx /= ((double)((int32_t)uii));
    }
    bufptr = width_force(8, bufptr, dtoa_g(dxx, bufptr));
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
    logstr(g_logbuf);
    goto score_report_qrange_next;
  }
  LOGPRINTFWW("--score: Results written to %s .\n", outname);
  while (0) {
  score_report_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
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
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --q-score-range data file has fewer tokens than\nexpected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  score_report_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --score file has fewer tokens than expected.\n", line_idx);
  score_report_ret_INVALID_FORMAT_2:
    logerrprintb();
  score_report_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 score_report_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  fclose_cond(infile);
  fclose_cond(outfile);
  return retval;
}

int32_t meta_analysis_open_and_read_header(const char* fname, char* loadbuf, uintptr_t loadbuf_size, char* sorted_header_dict, uint32_t* header_id_map, uintptr_t header_dict_ct, uintptr_t max_header_len, uint32_t weighted_z, uint32_t* token_ct_ptr, gzFile* gz_infile_ptr, uint32_t* col_skips, uint32_t* col_sequence, uintptr_t* line_idx_ptr, uint32_t* line_max_ptr) {
  uintptr_t line_idx = 1;
  uint32_t line_max = 0;
  uint32_t colnum = 0; // 0-based
  uint32_t token_ct = *token_ct_ptr;
  uint32_t best_chr_id = 0x9fffffffU;
  uint32_t best_var_id = 0x9fffffffU;
  uint32_t best_bp_id = 0x9fffffffU;
  uint32_t best_p_id = 0x9fffffffU;
  uint32_t best_se_id = 0x9fffffffU;
  uint32_t best_ess_id = 0x9fffffffU;
  uint32_t best_a1_id = 0x9fffffffU;
  uint32_t best_a2_id = 0x9fffffffU;
  int32_t retval = 0;
  uint32_t parse_table[9];
  char* bufptr;
  uint32_t slen;
  uint32_t uii;
  uint32_t ujj;
  int32_t ii;
  if (line_idx_ptr) {
    line_max = *line_max_ptr;
  }
  retval = gzopen_read_checked(fname, gz_infile_ptr);
  if (retval) {
    goto meta_analysis_open_and_read_header_ret_1;
  }

  while (1) {
    if (!gzgets(*gz_infile_ptr, loadbuf, loadbuf_size)) {
      if (!gzeof(*gz_infile_ptr)) {
	goto meta_analysis_open_and_read_header_ret_READ_FAIL;
      }
      sprintf(g_logbuf, "Error: %s is empty.\n", fname);
      goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
    }
    if (line_max && (!loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, fname);
	goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
      }
      goto meta_analysis_open_and_read_header_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (!is_eoln_kns(*bufptr)) {
      break;
    }
    if (line_max) {
      slen = strlen(bufptr) + (uintptr_t)(bufptr - loadbuf);
      if (slen >= line_max) {
	line_max = slen + 1;
      }
      line_idx++;
    }
  }
  fill_uint_one(token_ct, parse_table);
  do {
    slen = strlen_se(bufptr);
    ii = bsearch_str(bufptr, slen, sorted_header_dict, max_header_len, header_dict_ct);
    if (ii != -1) {
      uii = header_id_map[(uint32_t)ii];
      ujj = uii >> 28;
      switch (ujj) {
      case 0:
        if (parse_table[uii] != 0xffffffffU) {
          goto meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL;
	}
	parse_table[uii] = (colnum * 16) | uii;
	break;
      case 1:
	// variant ID
	if (uii < best_var_id) {
	  best_var_id = uii;
	  parse_table[0] = colnum * 16;
	} else if (uii == best_var_id) {
	  goto meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL;
	}
	break;
      case 2:
	// SE
	if (uii < best_se_id) {
          // bugfix (7 Jan 2020): this incorrectly assigned to best_p_id
	  best_se_id = uii;
	  parse_table[2] = (colnum * 16) | 2;
	} else if (uii == best_se_id) {
	  goto meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL;
	}
	break;
      case 3:
	// P
	if (uii < best_p_id) {
	  best_p_id = uii;
	  parse_table[3] = (colnum * 16) | 3;
	} else if (uii == best_p_id) {
	  goto meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL;
	}
	break;
      case 4:
	// effective sample size
	if (uii < best_ess_id) {
	  best_ess_id = uii;
	  parse_table[4] = (colnum * 16) | 4;
	} else if (uii == best_ess_id) {
	  goto meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL;
	}
	break;
      case 5:
	// chr
	if (uii < best_chr_id) {
	  best_chr_id = uii;
	  parse_table[5] = (colnum * 16) | 5;
	} else if (uii == best_chr_id) {
	  goto meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL;
	}
	break;
      case 6:
	// bp
	if (uii < best_bp_id) {
	  best_bp_id = uii;
	  parse_table[6] = (colnum * 16) | 6;
	} else if (uii == best_bp_id) {
	  goto meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL;
	}
	break;
      case 7:
	// A1
	if (uii < best_a1_id) {
	  best_a1_id = uii;
	  parse_table[7] = (colnum * 16) | 7;
	} else if (uii == best_a1_id) {
	  goto meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL;
	}
	break;
      case 8:
	// A2
	if (uii < best_a2_id) {
	  best_a2_id = uii;
	  parse_table[8] = (colnum * 16) | 8;
	} else if (uii == best_a2_id) {
	  goto meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL;
	}
	break;
      }
    }
    bufptr = skip_initial_spaces(&(bufptr[slen]));
    if (++colnum == 0x8000000) {
      // pathological case
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has too many columns.\n", line_idx, fname);
      goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
    }
  } while (!is_eoln_kns(*bufptr));
  if (line_idx_ptr) {
    slen = strlen(bufptr) + (uintptr_t)(bufptr - loadbuf);
    if (slen >= line_max) {
      line_max = slen + 1;
    }
    *line_idx_ptr = line_idx;
    *line_max_ptr = line_max;
    if (parse_table[0] == 0xffffffffU) {
      sprintf(g_logbuf, "Error: No variant ID field found in %s.\n", fname);
      goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
    } else if (parse_table[1] == 0xffffffffU) {
      sprintf(g_logbuf, "Error: No effect size field found in %s.\n", fname);
      goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
    } else if (parse_table[2] == 0xffffffffU) {
      sprintf(g_logbuf, "Error: No standard error field found in %s.\n", fname);
      goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
    } else if (weighted_z && (parse_table[3] == 0xffffffffU)) {
      sprintf(g_logbuf, "Error: No p-value field found in %s.\n", fname);
      goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
    } else if (weighted_z && (parse_table[4] == 0xffffffffU)) {
      sprintf(g_logbuf, "Error: No effective sample size field found in %s.\n", fname);
      goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
    } else if (token_ct > 5) {
      if (parse_table[5] == 0xffffffffU) {
	sprintf(g_logbuf, "Error: No CHR field found in %s.\n", fname);
	goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
      } else if (parse_table[6] == 0xffffffffU) {
	sprintf(g_logbuf, "Error: No BP field found in %s.\n", fname);
	goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
      } else if ((token_ct > 7) && (parse_table[7] == 0xffffffffU)) {
	sprintf(g_logbuf, "Error: No A1 allele field found in %s.\n", fname);
	goto meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW;
      }
    }
  }
  // suppress bogus gcc 4.4 warning, this is not performance-critical
  qsort(parse_table, token_ct, sizeof(int32_t), uintcmp);

  // bugfix: this caused a segfault in no-map case
  if ((!weighted_z) && (token_ct > 5)) {
    token_ct -= 2;
  }
  col_skips[0] = parse_table[0] >> 4;
  col_sequence[0] = parse_table[0] & 15;
  for (uii = 1; uii < token_ct; uii++) {
    ujj = parse_table[uii];
    if (ujj == 0xffffffffU) {
      // A2 is optional
      token_ct--;
      break;
    }
    col_skips[uii] = (ujj >> 4) - (parse_table[uii - 1] >> 4);
    col_sequence[uii] = ujj & 15;
  }
  *token_ct_ptr = token_ct;

  while (0) {
  meta_analysis_open_and_read_header_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  meta_analysis_open_and_read_header_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  meta_analysis_open_and_read_header_ret_DUPLICATE_HEADER_COL:
    bufptr[slen] = '\0';
    sprintf(g_logbuf, "Error: Duplicate column header '%s' in %s.\n", bufptr, fname);
  meta_analysis_open_and_read_header_ret_INVALID_FORMAT_WW:
    wordwrapb(0);
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 meta_analysis_open_and_read_header_ret_1:
  return retval;
}

uint32_t meta_analysis_allelic_match(const char* existing_a1ptr, char** token_ptrs, uint32_t token_ct, uint32_t a1lenp1, uint32_t a2lenp1) {
  // returns 1 on same-direction match, 2 on A1/A2 reverse match, 0 on mismatch
  const char* cur_a1ptr = token_ptrs[7];
  // memcmp is safe since loadbuf is positioned after htable/master_var_list.
  if (memcmp(existing_a1ptr, cur_a1ptr, a1lenp1)) {
    return ((token_ct & 1) && (!memcmp(existing_a1ptr, token_ptrs[8], a2lenp1) && (!memcmp(&(existing_a1ptr[a2lenp1]), cur_a1ptr, a1lenp1)))) * 2;
  }
  return ((!(token_ct & 1)) || (!memcmp(&(existing_a1ptr[a1lenp1]), token_ptrs[8], a2lenp1)));
}

static inline char* uint32_encode_5_hi_uchar(uint32_t uii, char* start) {
  // tried a few bit hacks here, but turns out nothing really beats this
  *start++ = (unsigned char)((uii >> 28) | 0x80);
  *start++ = (unsigned char)((uii >> 21) | 0x80);
  *start++ = (unsigned char)((uii >> 14) | 0x80);
  *start++ = (unsigned char)((uii >> 7) | 0x80);
  *start++ = (unsigned char)(uii | 0x80);
  return start;
}

static inline uint32_t uint32_decode_5_hi_uchar(const char* start) {
  uint32_t uii = ((uint32_t)((unsigned char)(*start++))) << 28;
  uii |= (((uint32_t)((unsigned char)(*start++))) & 0x7f) << 21;
  uii |= (((uint32_t)((unsigned char)(*start++))) & 0x7f) << 14;
  uii |= (((uint32_t)((unsigned char)(*start++))) & 0x7f) << 7;
  uii |= ((uint32_t)((unsigned char)(*start))) & 0x7f;
  return uii;
}

int32_t meta_analysis(char* input_fnames, char* chrfield_search_order, char* snpfield_search_order, char* bpfield_search_order, char* a1field_search_order, char* a2field_search_order, char* pfield_search_order, char* sefield_search_order, char* essfield_search_order, uint32_t flags, char* extractname, char* outname, char* outname_end, double output_min_p, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  gzFile gz_infile = nullptr;
  FILE* infile = nullptr;
  FILE* outfile = nullptr;
  char* cur_window_marker_ids = nullptr;
  char* sorted_extract_ids = nullptr;
  uintptr_t* duplicate_id_bitfield = nullptr;
  Ll_str** duplicate_id_htable = nullptr;
  uintptr_t header_dict_ct = 1; // BETA/OR
  uintptr_t max_header_len = 3;
  uintptr_t extract_ct = 0;
  uintptr_t max_extract_id_len = 0;
  uintptr_t extract_ctl = 0;
  uintptr_t final_variant_ct = 0;
  uintptr_t last_var_idx = 0;
  uintptr_t window_entry_base_cost = 2;
  uintptr_t duplicate_id_htable_max_alloc = 0;
  uint64_t rejected_ct = 0;
  double cur_p = 0.0;
  double cur_ess = 0.0;
  uint32_t max_var_id_len_p1 = 0;
  uint32_t max_combined_allele_len = 0;
  uint32_t use_map = 1 - ((flags / METAANAL_NO_MAP) & 1);
  uint32_t no_allele = flags & METAANAL_NO_ALLELE;
  uint32_t input_beta = flags & METAANAL_LOGSCALE;
  uint32_t report_all = flags & METAANAL_REPORT_ALL;
  uint32_t output_beta = flags & METAANAL_QT;
  uint32_t report_study_specific = flags & METAANAL_STUDY;
  uint32_t report_dups = flags & METAANAL_REPORT_DUPS;

  uint32_t weighted_z = (flags / METAANAL_WEIGHTED_Z) & 1;
  uint32_t parse_max = 3;
  uint32_t file_ct = 0;
  uint32_t combined_allele_len_byte_width = 0;
  uint32_t line_max = 1;
  uint32_t a1lenp1 = 0;
  uint32_t a2lenp1 = 0;
  uint32_t cur_chrom = 0;
  uint32_t cur_bp = 0;
  uint32_t cur_combined_allele_len = 0;
  uint32_t pass_idx = 0;
  int32_t retval = 0;
  char missing_geno = *g_missing_geno_ptr;
  const char problem_strings[][16] = {"BAD_CHR", "BAD_BP", "MISSING_A1", "MISSING_A2", "BAD_ES", "BAD_SE", "ALLELE_MISMATCH", "BAD_P", "BAD_ESS", "DUPLICATE"};

  // [0] = SNP
  // [1] = BETA/OR
  // [2] = SE
  // [3] = P ('weighted-z' only)
  // [4] = effective sample size ('weighted-z' only)
  // [5] = CHR (unless 'no-map')
  // [6] = BP (unless 'no-map')
  // [7] = A1 (unless 'no-map'/'no-allele')
  // [8] = A2 (optional, ignored with 'no-map'/'no-allele')
  char* token_ptrs[9];
  uint32_t col_skips[9];
  uint32_t col_sequence[9];
  uintptr_t loadbuf_size;
  uintptr_t max_var_id_len_p5;
  uintptr_t line_idx;
  uintptr_t master_var_entry_len;
  uintptr_t master_var_idx;
  uintptr_t total_data_slots;
  uintptr_t cur_data_slots;
  uintptr_t variants_remaining;
  uintptr_t cur_var_idx;
  uintptr_t first_var_idx;
  uintptr_t htable_write_limit;
  uintptr_t ulii;
  Ll_str** htable;
  Ll_str** ll_pptr;
  Ll_str* ll_ptr;
  Ll_str* htable_write;
  Ll_str* duplicate_id_htable_write;
  unsigned char* bigstack_mark2;
  char* sorted_header_dict;
  char* master_var_list;
  char* cur_entry_list_window;
  char* fname_ptr;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* wptr;
  double* cur_data;
  double* cur_data_ptr;
  uintptr_t* cur_data_index;
  uintptr_t* ulptr;
  uint32_t* header_id_map;
  uint32_t* uiptr;
  double cur_beta;
  double cur_se;
  double cur_inv_var;
  double numer;
  double denom;
  double denom2;
  double numer_random;
  double denom_random;
  double tau2;
  double meta_q;
  double varsum;
  double summ;
  double varsum_random;
  double summ_random;
  double summtest;
  double summtest_random;
  double p1;
  double pr;
  double pq;
  double meta_i;
  double dxx;
  uint32_t file_ct_byte_width;
  uint32_t file_ct_mask;
  uint32_t file_ct64;
  uint32_t file_idx;
  uint32_t cur_file_ct;
  uint32_t cur_file_ct_m1;
  uint32_t fname_len;
  uint32_t token_ct;
  uint32_t seq_idx;
  uint32_t var_id_len;
  uint32_t slen_base;
  uint32_t slen;
  uint32_t cur_variant_ct;
  uint32_t problem_mask;
  uint32_t uii;
  int32_t ii;
  {
    // 1. Construct header search dictionary.  Similar to clump_reports().
    if (snpfield_search_order) {
      // bugfix (3 Nov 2017): this needs to be +=, not =
      header_dict_ct += count_and_measure_multistr(snpfield_search_order, &max_header_len);
    } else {
      max_header_len = 4; // 'SNP' + null terminator
      header_dict_ct++;
    }
    if (sefield_search_order) {
      header_dict_ct += count_and_measure_multistr(sefield_search_order, &max_header_len);
    } else {
      header_dict_ct++;
    }
    if (input_beta && (max_header_len < 5)) {
      max_header_len = 5;
    }
    if (weighted_z) {
      parse_max = 5;
      if (pfield_search_order) {
	header_dict_ct += count_and_measure_multistr(pfield_search_order, &max_header_len);
      } else {
	header_dict_ct++;
      }
      if (essfield_search_order) {
	header_dict_ct += count_and_measure_multistr(essfield_search_order, &max_header_len);
      } else {
	header_dict_ct++;
	if (max_header_len < 6) {
	  max_header_len = 6;
	}
      }
    }
    if (use_map) {
      if (chrfield_search_order) {
        header_dict_ct += count_and_measure_multistr(chrfield_search_order, &max_header_len);
      } else {
        if (max_header_len < 4) {
          max_header_len = 4;
        }
        header_dict_ct++;
      }
      if (bpfield_search_order) {
        header_dict_ct += count_and_measure_multistr(bpfield_search_order, &max_header_len);
      } else {
        if (max_header_len < 3) {
          max_header_len = 3;
        }
        header_dict_ct++;
      }
      if (!no_allele) {
	if (a1field_search_order) {
	  header_dict_ct += count_and_measure_multistr(a1field_search_order, &max_header_len);
	} else {
	  header_dict_ct++;
	}
	if (a2field_search_order) {
	  header_dict_ct += count_and_measure_multistr(a2field_search_order, &max_header_len);
	} else {
	  header_dict_ct++;
	}
	parse_max = 9;
      } else {
	parse_max = 7;
      }
    }
    if (bigstack_alloc_c(header_dict_ct * max_header_len, &sorted_header_dict) ||
	bigstack_alloc_ui(header_dict_ct, &header_id_map)) {
      goto meta_analysis_ret_NOMEM;
    }
    ulii = 0; // write position
    if (snpfield_search_order) {
      bufptr = snpfield_search_order;
      uii = 0x10000000;
      do {
	slen = strlen(bufptr) + 1;
	memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, slen);
	header_id_map[ulii++] = uii++;
	bufptr = &(bufptr[slen]);
      } while (*bufptr);
    } else {
      memcpy(sorted_header_dict, "SNP", 4);
      header_id_map[0] = 0x10000000;
      ulii++;
    }
    if (!input_beta) {
      memcpyl3(&(sorted_header_dict[ulii * max_header_len]), "OR");
    } else {
      memcpy(&(sorted_header_dict[ulii * max_header_len]), "BETA", 5);
    }
    header_id_map[ulii++] = 1;
    if (sefield_search_order) {
      bufptr = sefield_search_order;
      uii = 0x20000000;
      do {
	slen = strlen(bufptr) + 1;
	memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, slen);
	header_id_map[ulii++] = uii++;
	bufptr = &(bufptr[slen]);
      } while (*bufptr);
    } else {
      memcpyl3(&(sorted_header_dict[ulii * max_header_len]), "SE");
      header_id_map[ulii++] = 0x20000000;
    }
    if (weighted_z) {
      if (pfield_search_order) {
	bufptr = pfield_search_order;
	uii = 0x30000000;
	do {
	  slen = strlen(bufptr) + 1;
	  memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, slen);
	  header_id_map[ulii++] = uii++;
	  bufptr = &(bufptr[slen]);
	} while (*bufptr);
      } else {
	memcpy(&(sorted_header_dict[ulii * max_header_len]), "P", 2);
	header_id_map[ulii++] = 0x30000000;
      }
      if (essfield_search_order) {
        // bugfix (3 Nov 2017): had pfield_search_order here
	bufptr = essfield_search_order;
	uii = 0x40000000;
	do {
	  slen = strlen(bufptr) + 1;
	  memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, slen);
	  header_id_map[ulii++] = uii++;
	  bufptr = &(bufptr[slen]);
	} while (*bufptr);
      } else {
	memcpy(&(sorted_header_dict[ulii * max_header_len]), "NMISS", 6);
	header_id_map[ulii++] = 0x40000000;
      }
    }
    if (use_map) {
      if (chrfield_search_order) {
        bufptr = chrfield_search_order;
        uii = 0x50000000;
        do {
          slen = strlen(bufptr) + 1;
          memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, slen);
          header_id_map[ulii++] = uii++;
          bufptr = &(bufptr[slen]);
        } while (*bufptr);
      } else {
        memcpy(&(sorted_header_dict[ulii * max_header_len]), "CHR", 4);
        header_id_map[ulii++] = 0x50000000;
      }
      if (bpfield_search_order) {
        bufptr = bpfield_search_order;
        uii = 0x60000000;
        do {
          slen = strlen(bufptr) + 1;
          memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, slen);
          header_id_map[ulii++] = uii++;
          bufptr = &(bufptr[slen]);
        } while (*bufptr);
      } else {
        memcpyl3(&(sorted_header_dict[ulii * max_header_len]), "BP");
        header_id_map[ulii++] = 0x60000000;
      }
      if (!no_allele) {
	if (a1field_search_order) {
	  bufptr = a1field_search_order;
	  uii = 0x70000000;
	  do {
	    slen = strlen(bufptr) + 1;
	    memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, slen);
	    header_id_map[ulii++] = uii++;
	    bufptr = &(bufptr[slen]);
	  } while (*bufptr);
	} else {
	  memcpyl3(&(sorted_header_dict[ulii * max_header_len]), "A1");
	  header_id_map[ulii++] = 0x70000000;
	}
	if (a2field_search_order) {
	  bufptr = a2field_search_order;
	  uii = 0x80000000U;
	  do {
	    slen = strlen(bufptr) + 1;
	    memcpy(&(sorted_header_dict[ulii * max_header_len]), bufptr, slen);
	    header_id_map[ulii++] = uii++;
	    bufptr = &(bufptr[slen]);
	  } while (*bufptr);
	} else {
	  memcpyl3(&(sorted_header_dict[ulii * max_header_len]), "A2");
	  header_id_map[ulii] = 0x80000000U;
	}
      }
    }
    if (qsort_ext(sorted_header_dict, header_dict_ct, max_header_len, strcmp_deref, (char*)header_id_map, sizeof(int32_t))) {
      goto meta_analysis_ret_NOMEM;
    }
    if (scan_for_duplicate_ids(sorted_header_dict, header_dict_ct, max_header_len)) {
      logerrprint("Error: Duplicate/invalid --meta-analysis-...-field field name.\n");
      goto meta_analysis_ret_INVALID_CMDLINE;
    }

    // 2. If --extract specified, load and sort permitted variant list.
    if (extractname) {
      if (fopen_checked(extractname, FOPEN_RB, &infile)) {
	goto meta_analysis_ret_OPEN_FAIL;
      }
      retval = scan_token_ct_len(MAXLINELEN, infile, g_textbuf, &extract_ct, &max_extract_id_len);
      if (retval) {
	goto meta_analysis_ret_1;
      }
      if (!extract_ct) {
	logerrprint("Error: Empty --extract file.\n");
	goto meta_analysis_ret_INVALID_FORMAT;
      }
      if (max_extract_id_len > MAX_ID_BLEN) {
	logerrprint("Error: --extract IDs are limited to " MAX_ID_SLEN_STR " characters.\n");
	goto meta_analysis_ret_INVALID_FORMAT;
      }
      if (bigstack_alloc_c(extract_ct * max_extract_id_len, &sorted_extract_ids)) {
	goto meta_analysis_ret_NOMEM;
      }
      rewind(infile);
      // Considered switching to a hash table, but decided against it for now
      // since it's less memory-efficient (in the usual case of similar-length
      // IDs), especially when lots of duplicate IDs are present.  Might be
      // worth revisiting this decision in the future, though, since there are
      // reasonable use cases involving 40-80 million line --extract files, and
      // skipping the sort step there is a big win.
      retval = read_tokens(MAXLINELEN, extract_ct, max_extract_id_len, infile, g_textbuf, sorted_extract_ids);
      if (retval) {
	goto meta_analysis_ret_1;
      }
      if (fclose_null(&infile)) {
	goto meta_analysis_ret_READ_FAIL;
      }
      qsort(sorted_extract_ids, extract_ct, max_extract_id_len, strcmp_casted);
      ulii = collapse_duplicate_ids(sorted_extract_ids, extract_ct, max_extract_id_len, nullptr);
      if (ulii < extract_ct) {
	extract_ct = ulii;
	bigstack_shrink_top(sorted_extract_ids, extract_ct * max_extract_id_len);
      }
      extract_ctl = BITCT_TO_WORDCT(extract_ct);
      if (bigstack_alloc_ul(extract_ctl, &duplicate_id_bitfield)) {
	goto meta_analysis_ret_NOMEM;
      }
    } else {
      duplicate_id_htable = (Ll_str**)bigstack_alloc(HASHMEM);
    }

    // 3. Allocate space for initial hash table.
    // Saving memory is pretty important here, so we use the following packing
    // in the ss field (W = byte width required to save numbers up to file_ct,
    // and M = 1 iff 'no-map' was not specified):
    // [W]: number of files this variant appears in minus 1, little-endian
    // [W+1]..[W+5], if M==1: chromosome byte followed by bp coordinate int;
    //                        may need to widen chromosome byte later
    // [W+5M+1]: null-terminated variant ID.  Followed by null-terminated A1/A2
    //           if 'no-allele' not specified
    bigstack_mark2 = g_bigstack_base;
    htable = (Ll_str**)bigstack_alloc(HASHMEM);
    if (!htable) {
      goto meta_analysis_ret_NOMEM;
    }
    for (uii = 0; uii < HASHSIZE; uii++) {
      htable[uii] = nullptr;
    }

    // 4. Initial scan: save all potentially valid variant IDs (and
    //    accompanying allele codes/chr/pos, if present) in the hash table, and
    //    produce .prob file.  Also determine maximum line length, for use in
    //    later passes.
    fname_ptr = input_fnames;
    do {
      fname_ptr = strchr(fname_ptr, '\0');
      fname_ptr++;
      file_ct++;
    } while (*fname_ptr);
    file_ct_byte_width = __builtin_clz(file_ct) / 8;
    file_ct_mask = 0xffffffffU >> (8 * file_ct_byte_width);
    file_ct_byte_width = 4 - file_ct_byte_width;
    file_ct64 = (file_ct + 63) / 64;

    slen_base = file_ct_byte_width;
    if (use_map) {
      slen_base += 5;
    }
    fname_ptr = input_fnames;
    htable_write = (Ll_str*)g_bigstack_base;
    bigstack_end_mark[-1] = ' ';
    for (file_idx = 0; file_idx < file_ct; file_idx++) {
      if (sorted_extract_ids) {
	fill_ulong_zero(extract_ctl, duplicate_id_bitfield);
      } else {
	for (uii = 0; uii < HASHSIZE; uii++) {
	  duplicate_id_htable[uii] = nullptr;
	}
      }
      fname_len = strlen(fname_ptr);
      // prevent overlap between loadbuf and new hash table entries.
      loadbuf_size = (((uintptr_t)(bigstack_end_mark - ((unsigned char*)htable_write))) / 4);
      if (loadbuf_size > MAXLINEBUFLEN) {
	loadbuf_size = MAXLINEBUFLEN;
      } else if (loadbuf_size <= MAXLINELEN) {
	goto meta_analysis_ret_NOMEM;
      }
      loadbuf = (char*)(&(bigstack_end_mark[-((intptr_t)loadbuf_size)]));
      duplicate_id_htable_write = (Ll_str*)loadbuf;
      htable_write_limit = ((uintptr_t)loadbuf) - loadbuf_size - 16;
      token_ct = parse_max;
      retval = meta_analysis_open_and_read_header(fname_ptr, loadbuf, loadbuf_size, sorted_header_dict, header_id_map, header_dict_ct, max_header_len, weighted_z, &token_ct, &gz_infile, col_skips, col_sequence, &line_idx, &line_max);
      if (retval) {
	goto meta_analysis_ret_1;
      }
      while (1) {
	line_idx++;
	if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	  if (!gzeof(gz_infile)) {
	    goto meta_analysis_ret_READ_FAIL;
	  }
	  break;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  if (loadbuf_size == MAXLINEBUFLEN) {
	    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, fname_ptr);
	    goto meta_analysis_ret_INVALID_FORMAT_WW;
	  }
	  goto meta_analysis_ret_NOMEM;
	}
	bufptr = skip_initial_spaces(loadbuf);
	if (is_eoln_kns(*bufptr)) {
	  slen = strlen(bufptr) + ((uintptr_t)(bufptr - loadbuf));
	  if (slen >= line_max) {
	    line_max = slen + 1;
	  }
	  continue;
	}
	bufptr = next_token_multz(bufptr, col_skips[0]);
	token_ptrs[col_sequence[0]] = bufptr;
	for (seq_idx = 1; seq_idx < token_ct; seq_idx++) {
	  bufptr = next_token_mult(bufptr, col_skips[seq_idx]);
	  token_ptrs[col_sequence[seq_idx]] = bufptr;
	}
	if (!bufptr) {
	  // PLINK 1.07 doesn't error out here, or even count the number of
	  // instances
	  slen = strlen(loadbuf);
	  if (slen >= line_max) {
	    line_max = slen + 1;
	  }
	  continue;
	}
	slen = strlen(bufptr) + ((uintptr_t)(bufptr - loadbuf));
	if (slen >= line_max) {
	  line_max = slen + 1;
	}
	bufptr = token_ptrs[0];
	var_id_len = strlen_se(bufptr);
	if (var_id_len > MAX_ID_SLEN) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has an excessively long variant ID.\n", line_idx, fname_ptr);
	  goto meta_analysis_ret_INVALID_FORMAT_WW;
	}
	bufptr[var_id_len] = '\0';
	uii = hashval2(bufptr, var_id_len++);
	// var_id_len now includes null-terminator
	if (sorted_extract_ids) {
	  ii = bsearch_str(bufptr, var_id_len - 1, sorted_extract_ids, max_extract_id_len, extract_ct);
	  if (ii == -1) {
	    continue;
	  }
	  if (is_set(duplicate_id_bitfield, ii)) {
	    problem_mask = 0x200;
	    goto meta_analysis_report_error;
	  }
	  set_bit(ii, duplicate_id_bitfield);
	} else {
	  ll_pptr = &(duplicate_id_htable[uii]);
	  while (1) {
	    ll_ptr = *ll_pptr;
	    if ((!ll_ptr) || (!strcmp(bufptr, ll_ptr->ss))) {
	      break;
	    }
	    ll_pptr = &(ll_ptr->next);
	  }
	  if (ll_ptr) {
	    problem_mask = 0x200;
	    goto meta_analysis_report_error;
	  }
	  // word-align for now
	  // note that it is NOT safe to use uii here.
	  ulii = sizeof(intptr_t) + round_up_pow2(var_id_len, BYTECT);
	  if (((uintptr_t)htable_write) + ulii > ((uintptr_t)duplicate_id_htable_write)) {
	    goto meta_analysis_ret_NOMEM;
	  }
	  duplicate_id_htable_write = (Ll_str*)(((uintptr_t)duplicate_id_htable_write) - ulii);
	  *ll_pptr = duplicate_id_htable_write;
	  duplicate_id_htable_write->next = nullptr;
	  memcpy(duplicate_id_htable_write->ss, bufptr, var_id_len);
	}
	ll_pptr = &(htable[uii]);

	// validate
	problem_mask = 0;
	if (use_map) {
	  ii = get_chrom_code_destructive(chrom_info_ptr, token_ptrs[5]);
	  if (ii < 0) {
	    problem_mask = 1;
	  } else {
	    cur_chrom = (uint32_t)ii;
	    if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom)) {
	      continue;
	    }
	  }
	  if (scan_uint_defcap(token_ptrs[6], &cur_bp)) {
	    problem_mask |= 2;
	  }
	  if (!no_allele) {
	    bufptr = token_ptrs[7];
	    a1lenp1 = strlen_se(bufptr); // not +1 yet
	    if ((*bufptr == missing_geno) && (a1lenp1 == 1)) {
	      problem_mask |= 4;
	    }
	    bufptr[a1lenp1++] = '\0';
	    // A2 allele present if token_ct == 7 or == 9
	    // if we make further extensions to this function, we should
	    // replace "token_ct & 1" with an a2_present boolean
	    if (token_ct & 1) {
	      bufptr = token_ptrs[8];
	      a2lenp1 = strlen_se(bufptr);
	      if ((*bufptr == missing_geno) && (a2lenp1 == 1)) {
		problem_mask |= 8;
	      }
	      bufptr[a2lenp1++] = '\0';
	    }
	  }
	}
	if (scan_double(token_ptrs[1], &cur_beta) || (cur_beta == INFINITY) || ((!input_beta) && (!(cur_beta >= 0))) || (input_beta && ((cur_beta != cur_beta) || (cur_beta == -INFINITY)))) {
	  problem_mask |= 0x10;
	}
	if (scan_double(token_ptrs[2], &cur_se) || (!(cur_se >= 0.0)) || (cur_se == INFINITY)) {
	  problem_mask |= 0x20;
	}
	if (weighted_z) {
	  if (scan_double(token_ptrs[3], &cur_p) || (!(cur_p >= 0.0)) || (cur_p > 1.0)) {
	    problem_mask |= 0x80;
	  }
	  if (scan_double(token_ptrs[4], &cur_ess) || (!(cur_ess > 0.0)) || (cur_ess == INFINITY)) {
	    problem_mask |= 0x100;
	  }
	}
	// check main hash table
	bufptr = token_ptrs[0];
	while (1) {
	  ll_ptr = *ll_pptr;
	  if ((!ll_ptr) || (!strcmp(bufptr, &(ll_ptr->ss[slen_base])))) {
	    break;
	  }
	  ll_pptr = &(ll_ptr->next);
	}

	if (!ll_ptr) {
	  if (problem_mask) {
	    goto meta_analysis_report_error;
	  }
	  // new hash table entry; word-align the allocation for now
	  ll_ptr = htable_write;
	  *ll_pptr = ll_ptr;
	  ll_ptr->next = nullptr;
	  wptr = memseta(ll_ptr->ss, 0, file_ct_byte_width);
	  if (use_map) {
	    *wptr++ = cur_chrom;
	    wptr = memcpya(wptr, &cur_bp, 4);
	  }
	  wptr = memcpya(wptr, bufptr, var_id_len);
	  if (var_id_len > max_var_id_len_p1) {
	    max_var_id_len_p1 = var_id_len;
	  }
	  if (!no_allele) {
	    bufptr = wptr;
	    wptr = memcpya(wptr, token_ptrs[7], a1lenp1);
	    if (token_ct & 1) {
	      wptr = memcpya(wptr, token_ptrs[8], a2lenp1);
	    } else {
	      *wptr++ = '\0';
	    }
	    uii = (uintptr_t)(wptr - bufptr);
	    if (uii > max_combined_allele_len) {
	      max_combined_allele_len = uii;
	    }
	  }
	  if (report_all) {
	    final_variant_ct++;
	  }
	  htable_write = (Ll_str*)round_up_pow2((uintptr_t)wptr, sizeof(uintptr_t));
	  if ((((uintptr_t)htable_write) > ((uintptr_t)duplicate_id_htable_write)) || (((uintptr_t)htable_write) > htable_write_limit)) {
	    goto meta_analysis_ret_NOMEM;
	  }
	} else {
	  if ((token_ct - 2 * weighted_z < 6) || meta_analysis_allelic_match(&(ll_ptr->ss[slen_base + var_id_len]), token_ptrs, token_ct, a1lenp1, a2lenp1)) {
	    if (problem_mask) {
	      goto meta_analysis_report_error;
	    }
	    // increment file count.  Assume little-endian machine
	    uiptr = (uint32_t*)ll_ptr->ss;
	    uii = (*uiptr) & file_ct_mask;
	    if ((!report_all) && (!uii)) {
	      final_variant_ct++;
	    }
	    uii++;
	    memcpy(ll_ptr->ss, &uii, file_ct_byte_width);
	  } else {
	    problem_mask |= 0x40;
	  meta_analysis_report_error:
	    if ((problem_mask == 0x200) && (!report_dups)) {
	      continue;
	    }
	    if (!outfile) {
	      memcpy(outname_end, ".prob", 6);
	      if (fopen_checked(outname, "w", &outfile)) {
		goto meta_analysis_ret_OPEN_FAIL;
	      }
	    }
	    bufptr = memcpyax(g_textbuf, fname_ptr, fname_len, '\t');
	    bufptr = memcpyax(bufptr, token_ptrs[0], var_id_len - 1, '\t');
	    do {
	      wptr = strcpyax(bufptr, problem_strings[__builtin_ctz(problem_mask)], '\n');
	      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
		goto meta_analysis_ret_WRITE_FAIL;
	      }
	      problem_mask &= problem_mask - 1;
	    } while (problem_mask);
	    rejected_ct++;
	  }
	}
      }
      if (gzclose(gz_infile) != Z_OK) {
	gz_infile = nullptr;
	goto meta_analysis_ret_READ_FAIL;
      }
      gz_infile = nullptr;
      if (!sorted_extract_ids) {
	ulii = ((uintptr_t)loadbuf) - ((uintptr_t)duplicate_id_htable_write);
	if (ulii > duplicate_id_htable_max_alloc) {
	  duplicate_id_htable_max_alloc = ulii;
	}
      }
      fname_ptr = &(fname_ptr[fname_len + 1]);
    }
    if (outfile) {
      if (fclose_null(&outfile)) {
	goto meta_analysis_ret_WRITE_FAIL;
      }
      LOGPRINTFWW("--meta-analysis: %" PRIu64 " problematic line%s; see %s .\n", rejected_ct, (rejected_ct == 1)? "" : "s", outname);
    }

    // 5. Determine final set of variants, and sort them (by chromosome, then
    //    position, then variant ID in natural order).  file_ct and (usually)
    //    [A1+A2 len] are also included past the end of each entry, to remove
    //    the need for an auxiliary index and let us free the hash table.
    if (!final_variant_ct) {
      logerrprint("Error: No --meta-analysis variants.\n");
      goto meta_analysis_ret_INVALID_CMDLINE;
#ifdef __LP64__
    } else if (final_variant_ct > 0x7fffffff) {
      logerrprint("Error: Too many distinct --meta-analysis variants (max 2^31 - 1).\n");
#endif
    }
    if (!no_allele) {
      combined_allele_len_byte_width = 4 - (__builtin_clz(max_combined_allele_len) / 8);
    }
    // bp coordinate, if present, expands from 4 to 5 bytes
    master_var_entry_len = slen_base + use_map + max_var_id_len_p1 + combined_allele_len_byte_width;
    loadbuf_size = round_up_pow2(line_max, END_ALLOC_CHUNK);
    loadbuf = (char*)bigstack_end_alloc_presized(loadbuf_size);
    if ((!loadbuf) ||
	bigstack_end_alloc_c(final_variant_ct * master_var_entry_len, &master_var_list) ||
	(((uintptr_t)htable_write) > ((uintptr_t)master_var_list))) {
      goto meta_analysis_ret_NOMEM;
    }
    // instead of following hash table pointers, we just plow through the table
    // entries in the order they were allocated in; this lets us access memory
    // sequentially
    ll_ptr = (Ll_str*)g_bigstack_base;
    for (master_var_idx = 0; master_var_idx < final_variant_ct;) {
      cur_file_ct_m1 = 0; // clear high bits
      memcpy(&cur_file_ct_m1, ll_ptr->ss, file_ct_byte_width);
      if (report_all || cur_file_ct_m1) {
	wptr = &(master_var_list[master_var_idx * master_var_entry_len]);
	master_var_idx++;
	if (use_map) {
	  *wptr++ = ll_ptr->ss[file_ct_byte_width];
	  memcpy(&uii, &(ll_ptr->ss[file_ct_byte_width + 1]), 4);
	  wptr = uint32_encode_5_hi_uchar(uii, wptr);
	}
	bufptr = &(ll_ptr->ss[slen_base]);
	slen = strlen(bufptr) + 1;
	wptr = memcpya(wptr, bufptr, slen);
	wptr = memcpya(wptr, &cur_file_ct_m1, file_ct_byte_width);
	bufptr = &(bufptr[slen]);
	if (!no_allele) {
	  // only save allele length sum, including null terminators
	  slen = strlen(bufptr) + 1;
	  slen += strlen(&(bufptr[slen])) + 1;
	  memcpy(wptr, &slen, combined_allele_len_byte_width);
	  bufptr = &(bufptr[slen]);
	}
      } else {
	bufptr = (char*)memchr(&(ll_ptr->ss[slen_base]), 0, max_var_id_len_p1);
	if (!no_allele) {
	  bufptr = (char*)memchr(&(bufptr[1]), 0, max_combined_allele_len);
	  bufptr = (char*)memchr(&(bufptr[1]), 0, max_combined_allele_len);
	}
	bufptr++;
      }
      // now bufptr points to the byte past the end of the hash table entry
      // allocation, and we know the next allocation starts at [this byte,
      // rounded up to nearest word boundary]
      ll_ptr = (Ll_str*)round_up_pow2((uintptr_t)bufptr, sizeof(intptr_t));
    }
    qsort(master_var_list, final_variant_ct, master_var_entry_len, strcmp_natural);
    // don't need htable anymore
    bigstack_reset(bigstack_mark2);
    if (!sorted_extract_ids) {
      bigstack_alloc(duplicate_id_htable_max_alloc);
    }
    total_data_slots = bigstack_left() / sizeof(uintptr_t);

    // 6. Remaining load passes: determine how many remaining variants' worth
    //    of effect sizes/SEs/Ps/ESSes fit in memory, load and meta-analyze
    //    just those variants, rinse and repeat.
    memcpy(outname_end, ".meta", 6);
    if (fopen_checked(outname, "w", &outfile)) {
      goto meta_analysis_ret_OPEN_FAIL;
    }
    if (use_map) {
      fputs(" CHR          BP", outfile);
    }
    fputs("            SNP", outfile);
    if (!no_allele) {
      fputs("  A1  A2", outfile);
    }
    fputs(output_beta? "   N           P        P(R)    BETA BETA(R)       Q       I" : "   N           P        P(R)      OR   OR(R)       Q       I", outfile);
    if (weighted_z) {
      fputs("  WEIGHTED_Z       P(WZ)", outfile);
    }
    if (report_study_specific) {
      for (file_idx = 0; file_idx < file_ct; file_idx++) {
	g_textbuf[0] = ' ';
	g_textbuf[1] = 'F';
	wptr = uint32toa(file_idx, &(g_textbuf[2]));
	wptr = width_force(8, g_textbuf, wptr);
	if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	  goto meta_analysis_ret_WRITE_FAIL;
	}
      }
    }
    putc_unlocked('\n', outfile);

    cur_data_index = (uintptr_t*)g_bigstack_base;
    if (use_map) {
      // chr/bp values can be discordant; when they are, we can't directly
      // search master_var_list for variant IDs.  Instead, we populate
      // cur_window_marker_ids with an ASCII-sorted list (marker ID,
      // cur_var_idx) tuples.
      window_entry_base_cost += (max_var_id_len_p1 + sizeof(int32_t) + sizeof(intptr_t) - 1) / sizeof(intptr_t);
    }
    max_var_id_len_p5 = max_var_id_len_p1 + 4;
    while (1) {
      first_var_idx = last_var_idx;
      // memory requrirements per current-window variant:
      // - 2 * sizeof(intptr_t) for cur_data pointer and current file write
      //     index; this grows from bottom of stack, while pointed-to stuff is
      //     allocated from top
      //   (technically could update the file write indexes in-place, but this
      //   part is not memory-critical so I doubt it's worth it.)
      // - 2 * sizeof(double) * file_ct for effect sizes and SEs; filled from
      //     back to front
      //   sometimes, numerator and squared denominator of weighted Z-score
      //   sometimes, bitfield describing which files are involved
      //   sometimes, combined_allele_len for A1/A2, sizeof(double)-aligned
      // - sometimes, (max_var_id_len_p1 + sizeof(int32_t)) rounded up, for
      //   cur_window_marker_ids.
      cur_entry_list_window = &(master_var_list[last_var_idx * master_var_entry_len]);
      bufptr = cur_entry_list_window;
      variants_remaining = final_variant_ct - last_var_idx;
      if (use_map) {
	bufptr = &(bufptr[6]); // ignore chromosome/position here
      }
      cur_data = (double*)(&(cur_data_index[total_data_slots]));
      ulii = 0;
      for (cur_variant_ct = 0; cur_variant_ct < variants_remaining; cur_variant_ct++) {
	bufptr2 = &(bufptr[cur_variant_ct * master_var_entry_len]);
	bufptr2 = (char*)memchr(bufptr2, 0, master_var_entry_len);
	bufptr2++;
	cur_file_ct_m1 = 0;
	memcpy(&cur_file_ct_m1, bufptr2, file_ct_byte_width);
	cur_data_slots = 0;
	if (report_study_specific) {
#ifdef __LP64__
	  cur_data_slots += file_ct64;
#else
	  cur_data_slots += 2 * file_ct64;
#endif
	}
	if (!no_allele) {
	  cur_combined_allele_len = 0;
	  memcpy(&cur_combined_allele_len, &(bufptr2[file_ct_byte_width]), combined_allele_len_byte_width);
	  cur_data_slots += (8 / BYTECT) * ((cur_combined_allele_len + 7) / 8);
	}
	cur_data_ptr = &(cur_data[-((intptr_t)cur_data_slots)]);
	cur_data_slots += window_entry_base_cost + (16 / BYTECT) * (cur_file_ct_m1 + 1 + weighted_z);
	ulii += cur_data_slots;
	if (ulii > total_data_slots) {
	  break;
	}
	if (report_study_specific) {
	  fill_ulong_zero(file_ct64 * (8 / BYTECT), (uintptr_t*)cur_data_ptr);
	}
	if (weighted_z) {
	  cur_data[-2] = 0.0;
	  cur_data[-1] = 0.0;
	}
	cur_data = &(cur_data[-((intptr_t)(cur_data_slots - window_entry_base_cost))]);
	// <effect sizes/SEs, reverse order> [WZ] [file idx bitfield] [A1/A2]
	//                                       ^
	//                                       |
	//                                  cur_data_ptr
	//
	// cur_data_index[2 * var_idx + 1] = # of effect sizes/etc. saved so
	// far
	cur_data_index[2 * cur_variant_ct] = (uintptr_t)cur_data_ptr;
	cur_data_index[2 * cur_variant_ct + 1] = 0;
      }
      if (!cur_variant_ct) {
	goto meta_analysis_ret_NOMEM;
      }
      last_var_idx += cur_variant_ct;
      if (use_map) {
	// position cur_window_marker_ids on top of cur_data_index
	cur_window_marker_ids = (char*)(&(cur_data_index[2 * cur_variant_ct]));
	// note that bufptr is positioned properly for reading variant IDs,
	// though it won't be after this loop
	bufptr2 = cur_window_marker_ids;
	for (uii = 0; uii < cur_variant_ct; uii++) {
	  strcpy(bufptr2, bufptr);
	  memcpy(&(bufptr2[max_var_id_len_p1]), &uii, 4);
	  bufptr = &(bufptr[master_var_entry_len]);
	  bufptr2 = &(bufptr2[max_var_id_len_p5]);
	}
	qsort(cur_window_marker_ids, cur_variant_ct, max_var_id_len_p5, strcmp_casted);
      }
      fname_ptr = input_fnames;
      for (file_idx = 0; file_idx < file_ct; file_idx++) {
	if (sorted_extract_ids) {
	  fill_ulong_zero(extract_ctl, duplicate_id_bitfield);
	} else {
	  for (uii = 0; uii < HASHSIZE; uii++) {
	    duplicate_id_htable[uii] = nullptr;
	  }
	}
	duplicate_id_htable_write = (Ll_str*)bigstack_mark2;
	fname_len = strlen(fname_ptr);
	token_ct = parse_max;
	retval = meta_analysis_open_and_read_header(fname_ptr, loadbuf, loadbuf_size, sorted_header_dict, header_id_map, header_dict_ct, max_header_len, weighted_z, &token_ct, &gz_infile, col_skips, col_sequence, nullptr, nullptr);
	if (retval) {
	  goto meta_analysis_ret_1;
	}
	while (1) {
	  // yeah, this is repetitive
	  if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	    if (!gzeof(gz_infile)) {
	      goto meta_analysis_ret_READ_FAIL;
	    }
	    break;
	  }
	  bufptr = skip_initial_spaces(loadbuf);
	  if (is_eoln_kns(*bufptr)) {
	    continue;
	  }
	  bufptr = next_token_multz(bufptr, col_skips[0]);
	  token_ptrs[col_sequence[0]] = bufptr;
	  for (seq_idx = 1; seq_idx < token_ct; seq_idx++) {
	    bufptr = next_token_mult(bufptr, col_skips[seq_idx]);
	    token_ptrs[col_sequence[seq_idx]] = bufptr;
	  }
	  if (!bufptr) {
	    continue;
	  }
	  bufptr = token_ptrs[0];
	  var_id_len = strlen_se(bufptr);
	  if (var_id_len >= max_var_id_len_p1) {
	    continue;
	  }
	  bufptr[var_id_len] = '\0';
	  if (sorted_extract_ids) {
	    ii = bsearch_str(bufptr, var_id_len, sorted_extract_ids, max_extract_id_len, extract_ct);
	    if (ii == -1) {
	      continue;
	    }
	    if (is_set(duplicate_id_bitfield, ii)) {
	      continue;
	    }
	    set_bit(ii, duplicate_id_bitfield);
	  } else {
	    uii = hashval2(bufptr, var_id_len);
	    ll_pptr = &(duplicate_id_htable[uii]);
	    while (1) {
	      ll_ptr = *ll_pptr;
	      if ((!ll_ptr) || (!strcmp(bufptr, ll_ptr->ss))) {
		break;
	      }
	      ll_pptr = &(ll_ptr->next);
	    }
	    if (ll_ptr) {
	      continue;
	    }
	    *ll_pptr = duplicate_id_htable_write;
	    duplicate_id_htable_write->next = nullptr;
	    memcpy(duplicate_id_htable_write->ss, bufptr, var_id_len + 1);
	    ulii = sizeof(intptr_t) + ((var_id_len + BYTECT) & (~(BYTECT - 1)));
	    duplicate_id_htable_write = (Ll_str*)(((uintptr_t)duplicate_id_htable_write) + ulii);
	  }
	  if (use_map) {
	    ii = get_chrom_code_destructive(chrom_info_ptr, token_ptrs[5]);
	    if (ii < 0) {
	      continue;
	    }
	    cur_chrom = (uint32_t)ii;
	    if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom)) {
	      continue;
	    }
	    if (scan_uint_defcap(token_ptrs[6], &cur_bp)) {
	      continue;
	    }
	    if (!no_allele) {
	      bufptr = token_ptrs[7];
	      a1lenp1 = strlen_se(bufptr);
	      if ((*bufptr == missing_geno) && (a1lenp1 == 1)) {
		continue;
	      }
	      bufptr[a1lenp1++] = '\0';
	      if (token_ct & 1) {
		bufptr = token_ptrs[8];
		a2lenp1 = strlen_se(bufptr);
		if ((*bufptr == missing_geno) && (a2lenp1 == 1)) {
		  continue;
		}
		bufptr[a2lenp1++] = '\0';
	      }
	    }
	  }
	  if (scan_double(token_ptrs[1], &cur_beta)) {
	    continue;
	  }
	  if (!input_beta) {
	    cur_beta = log(cur_beta);
	  }
	  if (!realnum(cur_beta)) {
	    continue;
	  }
	  if (scan_double(token_ptrs[2], &cur_se) || (!(cur_se >= 0.0)) || (cur_se == INFINITY)) {
	    continue;
	  }
	  if (weighted_z) {
	    if (scan_double(token_ptrs[3], &cur_p) || (!(cur_p >= 0.0)) || (cur_p > 1.0)) {
	      continue;
	    }
	    if (scan_double(token_ptrs[4], &cur_ess) || (!(cur_ess > 0.0)) || (cur_ess == INFINITY)) {
	      continue;
	    }
	  }
	  bufptr = token_ptrs[0];
	  if (use_map) {
	    ii = bsearch_str(bufptr, var_id_len, cur_window_marker_ids, max_var_id_len_p5, cur_variant_ct);
	    if (ii == -1) {
	      continue;
	    }
#ifdef __LP64__
	    cur_var_idx = 0; // clear high 32 bits
#endif
	    memcpy(&cur_var_idx, &(cur_window_marker_ids[(((uint32_t)ii) * max_var_id_len_p5) + max_var_id_len_p1]), 4);
	  } else {
	    bufptr[var_id_len] = '\0';
	    ulii = (uint32_t)bsearch_str_natural(bufptr, cur_entry_list_window, master_var_entry_len, cur_variant_ct);
	    // this comparison catches -1 return value
	    if ((ulii >= last_var_idx) || (ulii < first_var_idx)) {
	      continue;
	    }
	    cur_var_idx = ulii - first_var_idx;
	  }
	  cur_data_ptr = (double*)cur_data_index[2 * cur_var_idx];
	  cur_file_ct_m1 = cur_data_index[2 * cur_var_idx + 1];
	  if (!no_allele) {
	    if (!report_study_specific) {
	      bufptr2 = (char*)cur_data_ptr;
	    } else {
	      bufptr2 = (char*)(&(cur_data_ptr[file_ct64]));
	    }
	    if (!cur_file_ct_m1) {
	      // save allele codes
	      bufptr2 = memcpya(bufptr2, token_ptrs[7], a1lenp1);
	      if (token_ct & 1) {
		bufptr2 = memcpya(bufptr2, token_ptrs[8], a2lenp1);
	      } else {
		*bufptr2++ = '\0';
	      }
	    } else {
	      // compare them
	      uii = meta_analysis_allelic_match(bufptr2, token_ptrs, token_ct, a1lenp1, a2lenp1);
	      if (!uii) {
		continue;
	      } else if (uii == 2) {
		cur_beta = -cur_beta;
	      }
	    }
	  }
	  if (report_study_specific) {
	    set_bit(file_idx, (uintptr_t*)cur_data_ptr);
	  }
	  if (weighted_z) {
	    dxx = ltqnorm(1.0 - cur_p * 0.5) * sqrt(cur_ess);
	    if (cur_beta > 0.0) {
	      cur_data_ptr[-2] += dxx;
	    } else {
	      cur_data_ptr[-2] -= dxx;
	    }
	    cur_data_ptr[-1] += cur_ess;
	  }
	  cur_data_ptr = &(cur_data_ptr[-2 * ((int32_t)(cur_file_ct_m1 + weighted_z))]);
	  cur_data_ptr[-2] = cur_beta;
	  cur_data_ptr[-1] = cur_se;
	  cur_data_index[2 * cur_var_idx + 1] += 1;
	}
	if (gzclose(gz_infile) != Z_OK) {
	  gz_infile = nullptr;
	  goto meta_analysis_ret_READ_FAIL;
	}
	gz_infile = nullptr;
	fname_ptr = &(fname_ptr[fname_len + 1]);
      }
      for (cur_var_idx = 0; cur_var_idx < cur_variant_ct; cur_var_idx++) {
	cur_data_ptr = (double*)cur_data_index[2 * cur_var_idx];
	cur_file_ct = cur_data_index[2 * cur_var_idx + 1];
	bufptr = &(cur_entry_list_window[cur_var_idx * master_var_entry_len]);
	wptr = g_textbuf;
	if (use_map) {
	  cur_chrom = (uint32_t)((unsigned char)(*bufptr++));
	  wptr = width_force(4, wptr, chrom_name_write(chrom_info_ptr, cur_chrom, wptr));
	  wptr = memseta(wptr, 32, 2);
	  cur_bp = uint32_decode_5_hi_uchar(bufptr);
	  bufptr = &(bufptr[5]);
	  wptr = uint32toa_w10(cur_bp, wptr);
	}
	*wptr++ = ' ';
	var_id_len = strlen(bufptr);
	// bleah, this column width was not adaptive
	wptr = fw_strcpyn(14, var_id_len, bufptr, wptr);
	if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	  goto meta_analysis_ret_WRITE_FAIL;
	}
	if (!no_allele) {
	  if (!report_study_specific) {
	    bufptr = (char*)cur_data_ptr;
	  } else {
	    bufptr = (char*)(&(cur_data_ptr[file_ct64]));
	  }
	  slen = strlen(bufptr);
	  putc_unlocked(' ', outfile);
	  if (slen == 1) {
	    putc_unlocked(' ', outfile);
	    putc_unlocked(' ', outfile);
	  } else if (slen == 2) {
	    putc_unlocked(' ', outfile);
	  }
	  bufptr2 = &(bufptr[slen]);
	  if (fwrite_checked(bufptr, bufptr2 - bufptr, outfile)) {
	    goto meta_analysis_ret_WRITE_FAIL;
	  }
	  bufptr2++;
	  if (*bufptr2) {
	    // bugfix: fputs_w4 does the wrong thing for 4+ character alleles.
	    // instead, we want a leading space, then fputs_w3.
	    if (!bufptr2[1]) {
	      fputs("   ", outfile);
	      putc_unlocked(bufptr2[0], outfile);
	    } else if (!bufptr2[2]) {
	      fputs("  ", outfile);
	      putc_unlocked(bufptr2[0], outfile);
	      putc_unlocked(bufptr2[1], outfile);
	    } else {
	      putc_unlocked(' ', outfile);
	      fputs(bufptr2, outfile);
	    }
	  } else {
	    fputs("   ?", outfile);
	  }
	}
	g_textbuf[0] = ' ';
	wptr = &(g_textbuf[1]);
	wptr = width_force(3, wptr, uint32toa(cur_file_ct, wptr));
	if (cur_file_ct >= 2) {
	  // and here's the actual computation.
	  numer = 0.0;
	  denom = 0.0;
	  denom2 = 0.0;
	  for (file_idx = 1; file_idx <= cur_file_ct; file_idx++) {
	    ii = ((int32_t)(file_idx + weighted_z)) * (-2);
	    cur_beta = cur_data_ptr[ii];
	    cur_se = cur_data_ptr[ii + 1];
	    cur_inv_var = 1.0 / (cur_se * cur_se);
	    numer += cur_inv_var * cur_beta;
	    denom += cur_inv_var;
	    denom2 += cur_inv_var * cur_inv_var;
	  }
	  varsum = 1.0 / denom;
	  summ = numer * varsum;
	  meta_q = 0.0;
	  for (file_idx = 1; file_idx <= cur_file_ct; file_idx++) {
	    ii = ((int32_t)(file_idx + weighted_z)) * (-2);
	    cur_beta = cur_data_ptr[ii];
	    cur_se = cur_data_ptr[ii + 1];
	    dxx = (cur_beta - summ) / cur_se;
	    meta_q += dxx * dxx;
	  }
	  dxx = (double)((int32_t)(cur_file_ct - 1));
	  tau2 = (meta_q - dxx) / (denom - denom2 / denom);
	  if (tau2 < 0.0) {
	    tau2 = 0.0;
	  }
	  numer_random = 0.0;
	  denom_random = 0.0;
	  for (file_idx = 1; file_idx <= cur_file_ct; file_idx++) {
	    ii = ((int32_t)(file_idx + weighted_z)) * (-2);
	    cur_beta = cur_data_ptr[ii];
	    cur_se = cur_data_ptr[ii + 1];
	    cur_inv_var = 1.0 / (cur_se * cur_se + tau2);
	    numer_random += cur_inv_var * cur_beta;
	    denom_random += cur_inv_var;
	  }
	  varsum_random = 1.0 / denom_random;
	  summ_random = numer_random * varsum_random;
	  summtest = summ / sqrt(varsum);
	  summtest_random = summ_random / sqrt(varsum_random);
	  p1 = chiprob_p(summtest * summtest, 1);
	  pr = chiprob_p(summtest_random * summtest_random, 1);
	  pq = chiprob_p(meta_q, dxx);
	  meta_i = 100 * ((meta_q - dxx) / meta_q);
	  if (meta_i < 0.0) {
	    meta_i = 0.0;
	  } else if (meta_i > 100) {
	    meta_i = 100;
	  }
	  if (!output_beta) {
	    summ = exp(summ);
	    summ_random = exp(summ_random);
	  }
	  *wptr++ = ' ';
	  if (p1 >= 0.0) {
	    wptr = dtoa_g_wxp4x(MAXV(p1, output_min_p), 11, ' ', wptr);
	  } else {
	    wptr = memcpya(wptr, "         NA ", 12);
	  }
	  if (pr >= 0.0) {
	    wptr = dtoa_g_wxp4x(MAXV(pr, output_min_p), 11, ' ', wptr);
	  } else {
	    wptr = memcpya(wptr, "         NA ", 12);
	  }
	  wptr = dtoa_f_w7p4x(summ, ' ', wptr);
	  wptr = dtoa_f_w7p4x(summ_random, ' ', wptr);
	  if (pq >= 0.0) {
	    wptr = dtoa_f_w7p4x(MAXV(pq, output_min_p), ' ', wptr);
	  } else {
	    wptr = memcpya(wptr, "     NA ", 8);
	  }
	  wptr = width_force(7, wptr, dtoa_f_p2(meta_i, wptr));
	  if (weighted_z) {
	    numer = cur_data_ptr[-2];
	    denom2 = cur_data_ptr[-1];
	    dxx = numer / sqrt(denom2);
	    *wptr++ = ' ';
	    wptr = dtoa_g_wxp4x(dxx, 11, ' ', wptr);
	    dxx = 1.0 - 2 * fabs(normdist(fabs(dxx)) - 0.5);
	    wptr = dtoa_g_wxp4(MAXV(dxx, output_min_p), 11, wptr);
	  }
	} else {
          // quasi-bugfix (20 Dec 2017): under 'report-all', report basic
          // p/beta values when n=1 instead of forcing user to re-scrape their
          // input files.
          if (cur_file_ct) {
            ii = ((int32_t)(1 + weighted_z)) * (-2);
            cur_beta = cur_data_ptr[ii];
            cur_se = cur_data_ptr[ii + 1];
            summtest = cur_beta / cur_se;
            p1 = chiprob_p(summtest * summtest, 1);
            *wptr++ = ' ';
            if (p1 >= 0.0) {
              p1 = MAXV(p1, output_min_p);
              wptr = dtoa_g_wxp4x(p1, 11, ' ', wptr);
              wptr = dtoa_g_wxp4x(p1, 11, ' ', wptr);
            } else {
              wptr = memcpya(wptr, "         NA          NA ", 24);
            }
            if (!output_beta) {
              cur_beta = exp(cur_beta);
            }
            wptr = dtoa_f_w7p4x(cur_beta, ' ', wptr);
            wptr = dtoa_f_w7p4x(cur_beta, ' ', wptr);
	    wptr = memcpya(wptr, "     NA      NA", 15);
            if (weighted_z) {
              numer = cur_data_ptr[-2];
              denom2 = cur_data_ptr[-1];
              dxx = numer / sqrt(denom2);
              *wptr++ = ' ';
              wptr = dtoa_g_wxp4x(dxx, 11, ' ', wptr);
              dxx = 1.0 - 2 * fabs(normdist(fabs(dxx)) - 0.5);
              wptr = dtoa_g_wxp4(MAXV(dxx, output_min_p), 11, wptr);
            }
          } else {
	    wptr = memcpya(wptr, "          NA          NA      NA      NA      NA      NA", 56);
            if (weighted_z) {
	      wptr = memcpya(wptr, "          NA          NA", 24);
            }
          }
	}
	if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	  goto meta_analysis_ret_WRITE_FAIL;
	}
	if (report_study_specific) {
	  uii = 0;
	  ulptr = (uintptr_t*)cur_data_ptr;
	  for (file_idx = 0; file_idx < file_ct; file_idx++) {
	    if (is_set(ulptr, file_idx)) {
	      uii++;
	      dxx = cur_data_ptr[((int32_t)(uii + weighted_z)) * (-2)];
	      if (!output_beta) {
		// finish fixing PLINK 1.07 bug
		dxx = exp(dxx);
	      }
	      dtoa_f_w7p4x(dxx, '\0', &(g_textbuf[1]));
	      fputs(g_textbuf, outfile);
	    } else {
	      fputs("      NA", outfile);
	    }
	  }
	}
	putc_unlocked('\n', outfile);
      }
      if (last_var_idx == final_variant_ct) {
	break;
      }
      pass_idx++;
      printf("\r--meta-analysis: Pass %u complete (%" PRIu64 "%%).", pass_idx, (last_var_idx * ((uint64_t)100)) / final_variant_ct);
      fflush(stdout);
    }
    if (pass_idx) {
      putc_unlocked('\r', stdout);
    }
    LOGPRINTFWW("--meta-analysis: %" PRIuPTR " variant%s processed; results written to %s .\n", final_variant_ct, (final_variant_ct == 1)? "" : "s", outname);
  }

  while (0) {
  meta_analysis_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  meta_analysis_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  meta_analysis_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  meta_analysis_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  meta_analysis_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  meta_analysis_ret_INVALID_FORMAT_WW:
    wordwrapb(0);
    logerrprintb();
  meta_analysis_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 meta_analysis_ret_1:
  gzclose_cond(gz_infile);
  fclose_cond(infile);
  fclose_cond(outfile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return retval;
}
