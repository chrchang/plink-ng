#include "plink_common.h"

#include "plink_filter.h"
#include "plink_stats.h"

#include "pigz.h"

void oblig_missing_init(Oblig_missing_info* om_ip) {
  om_ip->cluster_ct = 0;
  om_ip->entry_ct = 0;
  om_ip->entries = NULL;
  om_ip->cluster_ref_cts = NULL;
  om_ip->sample_lookup = NULL;
  om_ip->marker_fname = NULL;
  om_ip->sample_fname = NULL;
}

void oblig_missing_cleanup(Oblig_missing_info* om_ip) {
  if (om_ip->marker_fname) {
    free_cond(om_ip->entries);
    free_cond(om_ip->cluster_ref_cts);
    free_cond(om_ip->sample_lookup);
    free_cond(om_ip->marker_fname);
    free_cond(om_ip->sample_fname);
    om_ip->marker_fname = NULL;
  }
}

const char keep_str[] = "keep";
const char keep_fam_str[] = "keep-fam";
const char remove_str[] = "remove";
const char remove_fam_str[] = "remove-fam";

const char* keep_or_remove_flag_str(uint32_t flags) {
  switch (flags) {
  case 0:
    return keep_str;
  case 1:
    return remove_str;
  case 2:
    return keep_fam_str;
  case 3:
    return remove_fam_str;
  }
  return NULL;
}

int32_t keep_or_remove(char* fname, char* sorted_ids, uintptr_t sorted_ids_ct, uintptr_t max_id_len, uint32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr, uint32_t flags) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t* exclude_arr_new = NULL;
  uintptr_t unfiltered_ctl = (unfiltered_ct + (BITCT - 1)) / BITCT;
  uintptr_t duplicate_ct = 0;
  uintptr_t line_idx = 0;
  uint32_t do_exclude = flags & 1;
  uint32_t families_only = flags & 2;
  int32_t retval = 0;
  char* id_buf;
  char* bufptr0;
  int32_t ii;
  uint32_t unfiltered_idx;
  uint32_t cur_idx;
  uint32_t last_idx;

  if (wkspace_alloc_ul_checked(&exclude_arr_new, unfiltered_ctl * sizeof(intptr_t))) {
    goto keep_or_remove_ret_NOMEM;
  }
  if (do_exclude) {
    // need this to avoid --keep/--remove spurious duplicate warnings (ID could
    // have been removed by an earlier step), though we could switch that to
    // the already_seen strategy
    memcpy(exclude_arr_new, exclude_arr, unfiltered_ctl * sizeof(intptr_t));
  } else {
    fill_all_bits(exclude_arr_new, unfiltered_ct);
  }
  if (fopen_checked(&infile, fname, "r")) {
    goto keep_or_remove_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  if (wkspace_alloc_c_checked(&id_buf, max_id_len)) {
    goto keep_or_remove_ret_NOMEM;
  }
  while (fgets(tbuf, MAXLINELEN, infile) != NULL) {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Line %" PRIuPTR " of --%s file is pathologically long.\n", line_idx, keep_or_remove_flag_str(flags));
      goto keep_or_remove_ret_INVALID_FORMAT_2;
    }
    bufptr0 = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr0)) {
      continue;
    }
    if (!families_only) {
      if (bsearch_read_fam_indiv(id_buf, sorted_ids, max_id_len, sorted_ids_ct, bufptr0, NULL, &ii)) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --%s file has fewer tokens than expected.\n", line_idx, keep_or_remove_flag_str(flags));
	goto keep_or_remove_ret_INVALID_FORMAT_2;
      }
      if (ii != -1) {
	unfiltered_idx = id_map[(uint32_t)ii];
	if (!IS_SET(exclude_arr, unfiltered_idx)) {
	  if (do_exclude) {
	    if (IS_SET(exclude_arr_new, unfiltered_idx)) {
	      duplicate_ct++;
	    } else {
	      SET_BIT(exclude_arr_new, unfiltered_idx);
	    }
	  } else {
	    if (!IS_SET(exclude_arr_new, unfiltered_idx)) {
	      duplicate_ct++;
	    } else {
	      CLEAR_BIT(exclude_arr_new, unfiltered_idx);
	    }
	  }
	}
      }
    } else {
      bsearch_fam(id_buf, sorted_ids, max_id_len, sorted_ids_ct, bufptr0, &cur_idx, &last_idx);
      ii = 0;
      while (cur_idx < last_idx) {
	unfiltered_idx = id_map[cur_idx++];
	if (!IS_SET(exclude_arr, unfiltered_idx)) {
	  if (do_exclude) {
	    if (IS_SET(exclude_arr_new, unfiltered_idx)) {
	      ii = 1;
	    } else {
	      SET_BIT(exclude_arr_new, unfiltered_idx);
	    }
	  } else {
	    if (!IS_SET(exclude_arr_new, unfiltered_idx)) {
	      ii = 1;
	    } else {
	      CLEAR_BIT(exclude_arr_new, unfiltered_idx);
	    }
	  }
	}
      }
      // only add one to duplicate_ct, instead of one per family member
      duplicate_ct += (uint32_t)ii;
    }
  }
  if (!feof(infile)) {
    goto keep_or_remove_ret_READ_FAIL;
  }
  memcpy(exclude_arr, exclude_arr_new, unfiltered_ctl * sizeof(intptr_t));
  *exclude_ct_ptr = popcount_longs(exclude_arr, unfiltered_ctl);
  if (*exclude_ct_ptr == unfiltered_ct) {
    LOGPRINTF("Error: No %s remaining after --%s.\n", g_species_plural, keep_or_remove_flag_str(flags));
    goto keep_or_remove_ret_ALL_SAMPLES_EXCLUDED;
  }
  unfiltered_ct -= *exclude_ct_ptr; // now filtered count
  LOGPRINTF("--%s: %" PRIuPTR " %s remaining.\n", keep_or_remove_flag_str(flags), unfiltered_ct, species_str(unfiltered_ct));
  if (duplicate_ct) {
    // "At least" since this does not count duplicate IDs absent from the .bim.
    LOGPRINTF("Warning: At least %" PRIuPTR " duplicate ID%s in --%s file.\n", duplicate_ct, (duplicate_ct == 1)? "" : "s", keep_or_remove_flag_str(flags));
  }
  while (0) {
  keep_or_remove_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  keep_or_remove_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  keep_or_remove_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  keep_or_remove_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  keep_or_remove_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  return retval;
}

void extract_exclude_process_token(const char* tok_start, const uint32_t* marker_id_htable, uint32_t marker_id_htable_size, const uint32_t* extra_alloc_base, const char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_exclude, uintptr_t* already_seen, uintptr_t* duplicate_ct_ptr, uint32_t do_exclude, uint32_t curtoklen) {
  if (curtoklen >= max_marker_id_len) {
    return;
  }
  uint32_t hashval = murmurhash3_32(tok_start, curtoklen) % marker_id_htable_size;
  uint32_t next_incr = 1;
  uint32_t cur_llidx = 0;
  const char* sptr;
  uintptr_t marker_uidx;
  uint32_t hash_result;
  uint32_t top_diff;
  uint32_t cur_dup;
  while (1) {
    hash_result = marker_id_htable[hashval];
    cur_dup = hash_result >> 31;
    if (cur_dup) {
      if (hash_result == 0xffffffffU) {
	return;
      }
      cur_llidx = hash_result << 1;
      marker_uidx = extra_alloc_base[cur_llidx];
    } else {
      // no need to check without --update-map, which requires uniqueness
      if (is_set(marker_exclude, hash_result)) {
	goto extract_exclude_process_token_next_hashval;
      }
      marker_uidx = hash_result;
    }
    sptr = &(marker_ids[marker_uidx * max_marker_id_len]);
    if ((!memcmp(tok_start, sptr, curtoklen)) && (!sptr[curtoklen])) {
      // when there are multiple variants with the same ID, only need to
      // check bit for first search result
      if (IS_SET(already_seen, marker_uidx)) {
	*duplicate_ct_ptr += 1;
      } else {
	SET_BIT(already_seen, marker_uidx);
	if (!cur_dup) {
	  return;
	}
	while (1) {
	  cur_llidx = extra_alloc_base[cur_llidx + 1];
	  if (cur_llidx == 0xffffffffU) {
	    return;
	  }
	  SET_BIT(already_seen, extra_alloc_base[cur_llidx]);
	}
      }
    }
  extract_exclude_process_token_next_hashval:
    top_diff = marker_id_htable_size - hashval;
    if (top_diff > next_incr) {
      hashval += next_incr;
    } else {
      hashval = next_incr - top_diff;
    }
    next_incr += 2;
  }
}

int32_t extract_exclude_flag_norange(char* fname, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, uint32_t do_exclude, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t duplicate_ct = 0;
  // needs to be synced with populate_id_htable
  const uint32_t* extra_alloc_base = &(marker_id_htable[CACHEALIGN32_INT32(marker_id_htable_size)]);
  char* midbuf = &(tbuf[MAXLINELEN]);
  uint32_t curtoklen = 0;
  int32_t retval = 0;
  uintptr_t bufsize;
  uintptr_t* already_seen;
  char* bufptr0;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  if (wkspace_alloc_ul_checked(&already_seen, unfiltered_marker_ctl * sizeof(intptr_t))) {
    goto extract_exclude_flag_norange_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, unfiltered_marker_ctl);
  if (fopen_checked(&infile, fname, "rb")) {
    goto extract_exclude_flag_norange_ret_OPEN_FAIL;
  }
  while (1) {
    if (fread_checked(midbuf, MAXLINELEN, infile, &bufsize)) {
      goto extract_exclude_flag_norange_ret_READ_FAIL;
    }
    if (!bufsize) {
      if (curtoklen) {
        extract_exclude_process_token(&(tbuf[MAXLINELEN - curtoklen]), marker_id_htable, marker_id_htable_size, extra_alloc_base, marker_ids, max_marker_id_len, marker_exclude, already_seen, &duplicate_ct, do_exclude, curtoklen);
      }
      break;
    }
    bufptr0 = &(midbuf[bufsize]);
    *bufptr0 = ' ';
    bufptr0[1] = '0';
    bufptr = &(tbuf[MAXLINELEN - curtoklen]);
    bufptr2 = midbuf;
    if (curtoklen) {
      goto extract_exclude_flag_norange_tok_start;
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
    extract_exclude_flag_norange_tok_start:
      while (*bufptr2 > ' ') {
	bufptr2++;
      }
      curtoklen = (uintptr_t)(bufptr2 - bufptr);
      if (bufptr2 == &(tbuf[MAXLINELEN * 2])) {
        if (curtoklen > MAXLINELEN) {
	  sprintf(logbuf, "Error: Excessively long ID in --%s file.\n", do_exclude? "exclude" : "extract");
          goto extract_exclude_flag_norange_ret_INVALID_FORMAT_2;
	}
	bufptr3 = &(tbuf[MAXLINELEN - curtoklen]);
        memcpy(bufptr3, bufptr, curtoklen);
	break;
      }
      extract_exclude_process_token(bufptr, marker_id_htable, marker_id_htable_size, extra_alloc_base, marker_ids, max_marker_id_len, marker_exclude, already_seen, &duplicate_ct, do_exclude, curtoklen);
      bufptr = &(bufptr2[1]);
    }
  }
  if (!feof(infile)) {
    goto extract_exclude_flag_norange_ret_READ_FAIL;
  }
  if (do_exclude) {
    bitfield_or(marker_exclude, already_seen, unfiltered_marker_ctl * sizeof(intptr_t));
  } else {
    bitfield_ornot(marker_exclude, already_seen, unfiltered_marker_ctl * sizeof(intptr_t));
    zero_trailing_bits(marker_exclude, unfiltered_marker_ct);
  }
  *marker_exclude_ct_ptr = popcount_longs(marker_exclude, unfiltered_marker_ctl);
  if (*marker_exclude_ct_ptr == unfiltered_marker_ct) {
    LOGPRINTF("Error: No variants remaining after --%s.\n", do_exclude? "exclude" : "extract");
    goto extract_exclude_flag_norange_ret_ALL_MARKERS_EXCLUDED;
  }
  unfiltered_marker_ct -= *marker_exclude_ct_ptr; // now filtered count
  LOGPRINTF("--%s: %" PRIuPTR " variant%s remaining.\n", do_exclude? "exclude" : "extract", unfiltered_marker_ct, (unfiltered_marker_ct == 1)? "" : "s");
  if (duplicate_ct) {
    // "At least" since this does not count duplicate IDs absent from the .bim.
    LOGPRINTF("Warning: At least %" PRIuPTR " duplicate ID%s in --%s file.\n", duplicate_ct, (duplicate_ct == 1)? "" : "s", do_exclude? "exclude" : "extract");
  }

  while (0) {
  extract_exclude_flag_norange_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  extract_exclude_flag_norange_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  extract_exclude_flag_norange_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  extract_exclude_flag_norange_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  extract_exclude_flag_norange_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  return retval;
}

int32_t filter_attrib(char* fname, char* condition_str, uint32_t* id_htable, uint32_t id_htable_size, char* item_ids, uintptr_t max_id_len, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr) {
  gzFile gz_infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t include_ct = 0;
  uintptr_t unfiltered_ctl = (unfiltered_ct + (BITCT - 1)) / BITCT;
  uintptr_t* cur_neg_matches = NULL;
  char* sorted_pos_match = NULL;
  char* sorted_neg_match = NULL;
  char* bufptr2 = NULL;
  uint32_t pos_match_ct = 0;
  uint32_t neg_match_ct = 0;
  uint32_t neg_match_ctl = 0;
  uintptr_t max_pos_match_len = 0;
  uintptr_t max_neg_match_len = 0;
  uintptr_t line_idx = 0;
  uint32_t is_neg = 0;
  int32_t retval = 0;
  uintptr_t* exclude_arr_new;
  uintptr_t* already_seen;
  char* loadbuf;
  char* cond_ptr;
  char* bufptr;
  uintptr_t loadbuf_size;
  uintptr_t pos_match_idx;
  uintptr_t neg_match_idx;
  uintptr_t ulii;
  uint32_t item_uidx;
  uint32_t pos_match_needed;
  uint32_t cur_neg_match_ct;
  int32_t sorted_idx;
  
  if (wkspace_alloc_ul_checked(&exclude_arr_new, unfiltered_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&already_seen, unfiltered_ctl * sizeof(intptr_t))) {
    goto filter_attrib_ret_NOMEM;
  }
  fill_all_bits(exclude_arr_new, unfiltered_ct);
  fill_ulong_zero(already_seen, unfiltered_ctl);
  if (condition_str) {
    // allow NULL condition_str; this means all samples/variants named in the
    // file are included
    cond_ptr = condition_str;
    while (1) {
      while (*cond_ptr == ',') {
	cond_ptr++;
      }
      if (*cond_ptr == '-') {
	cond_ptr++;
	if (*cond_ptr == ',') {
	  continue;
	} else if (*cond_ptr == '-') {
	  logprint("Error: --attrib condition cannot contain consecutive dashes.\n");
	  goto filter_attrib_ret_INVALID_CMDLINE;
	}
	is_neg = 1;
      } else if (!(*cond_ptr)) {
        // command-line parser actually comma-terminates the string for us, so
        // no need to check for null elsewhere
	break;
      }
      bufptr = strchr(cond_ptr, ',');
      ulii = (uintptr_t)(bufptr - cond_ptr);
      if (is_neg) {
	neg_match_ct++;
	if (ulii >= max_neg_match_len) {
	  max_neg_match_len = ulii + 1;
	}
      } else {
        pos_match_ct++;
	if (ulii >= max_pos_match_len) {
          max_pos_match_len = ulii + 1;
	}
      }
      cond_ptr = bufptr;
      is_neg = 0;
    }
    if (pos_match_ct) {
      if (wkspace_alloc_c_checked(&sorted_pos_match, max_pos_match_len * pos_match_ct)) {
	goto filter_attrib_ret_NOMEM;
      }
    }
    if (neg_match_ct) {
      neg_match_ctl = (neg_match_ct + (BITCT - 1)) / BITCT;
      if (wkspace_alloc_c_checked(&sorted_neg_match, max_neg_match_len * neg_match_ct) ||
	  wkspace_alloc_ul_checked(&cur_neg_matches, neg_match_ctl * sizeof(intptr_t))) {
        goto filter_attrib_ret_NOMEM;
      }
    }
    pos_match_idx = 0;
    neg_match_idx = 0;
    cond_ptr = condition_str;
    is_neg = 0;
    while (1) {
      while (*cond_ptr == ',') {
	cond_ptr++;
      }
      if (*cond_ptr == '-') {
	cond_ptr++;
	if (*cond_ptr == ',') {
	  continue;
	}
	is_neg = 1;
      } else if (!(*cond_ptr)) {
	break;
      }
      bufptr = strchr(cond_ptr, ',');
      ulii = (uintptr_t)(bufptr - cond_ptr);
      if (is_neg) {
	memcpyx(&(sorted_neg_match[neg_match_idx * max_neg_match_len]), cond_ptr, ulii, '\0');
	neg_match_idx++;
      } else {
	memcpyx(&(sorted_pos_match[pos_match_idx * max_pos_match_len]), cond_ptr, ulii, '\0');
        pos_match_idx++;
      }
      cond_ptr = bufptr;
      is_neg = 0;
    }
    if (pos_match_ct) {
      qsort(sorted_pos_match, pos_match_ct, max_pos_match_len, strcmp_casted);
      bufptr = scan_for_duplicate_ids(sorted_pos_match, pos_match_ct, max_pos_match_len);
      if (bufptr) {
	LOGPREPRINTFWW("Error: Duplicate attribute '%s' in --attrib argument.\n", bufptr);
	goto filter_attrib_ret_INVALID_CMDLINE_2;
      }
    }
    if (neg_match_ct) {
      qsort(sorted_neg_match, neg_match_ct, max_neg_match_len, strcmp_casted);
      bufptr = scan_for_duplicate_ids(sorted_neg_match, neg_match_ct, max_neg_match_len);
      if (bufptr) {
	LOGPREPRINTFWW("Error: Duplicate attribute '%s' in --attrib argument.\n", bufptr);
	goto filter_attrib_ret_INVALID_CMDLINE_2;
      }
      // actually may make sense to have same attribute as a positive and
      // negative condition, so we don't check for that
    }
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto filter_attrib_ret_NOMEM;
  }
  if (gzopen_checked(&gz_infile, fname, "rb")) {
    goto filter_attrib_ret_OPEN_FAIL;
  }
  if (gzbuffer(gz_infile, 131072)) {
    goto filter_attrib_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (1) {
    line_idx++;
    if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
      if (!gzeof(gz_infile)) {
	goto filter_attrib_ret_READ_FAIL;
      }
      break;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: Line %" PRIuPTR" of --attrib file is pathologically long.\n", line_idx);
        goto filter_attrib_ret_INVALID_FORMAT_2;
      }
      goto filter_attrib_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr2 = token_endnn(bufptr);
    cond_ptr = skip_initial_spaces(bufptr2);
    item_uidx = id_htable_find(bufptr, (uintptr_t)(bufptr2 - bufptr), id_htable, id_htable_size, item_ids, max_id_len);
    if ((item_uidx == 0xffffffffU) || is_set(exclude_arr, item_uidx)) {
      // miss_ct++;
      continue;
    }
    if (is_set(already_seen, item_uidx)) {
      *bufptr2 = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --attrib file.\n", bufptr);
      goto filter_attrib_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, item_uidx);
    pos_match_needed = pos_match_ct;
    cur_neg_match_ct = 0;
    if (neg_match_ct) {
      fill_ulong_zero(cur_neg_matches, neg_match_ctl);
    }
    while (!is_eoln_kns(*cond_ptr)) {
      bufptr2 = cond_ptr;
      bufptr = token_endnn(cond_ptr);
      ulii = (uintptr_t)(bufptr - bufptr2);
      cond_ptr = skip_initial_spaces(bufptr);
      if (pos_match_needed && (bsearch_str(bufptr2, ulii, sorted_pos_match, max_pos_match_len, pos_match_ct) != -1)) {
	pos_match_needed = 0;
      }
      if (cur_neg_match_ct < neg_match_ct) {
	sorted_idx = bsearch_str(bufptr2, ulii, sorted_neg_match, max_neg_match_len, neg_match_ct);
	if ((sorted_idx != -1) && (!is_set(cur_neg_matches, sorted_idx))) {
          cur_neg_match_ct++;
	  if (cur_neg_match_ct == neg_match_ct) {
	    // fail
	    pos_match_needed = 1;
            break;
	  }
          set_bit(cur_neg_matches, sorted_idx);
	}
      }
    }
    if (pos_match_needed) {
      continue;
    }
    // full negative match causes pos_match_needed to be set, so no further
    // check required
    clear_bit(exclude_arr_new, item_uidx);
    include_ct++;
  }
  if (!include_ct) {
    logprint("Error: No variants remaining after --attrib.\n");
    retval = RET_ALL_MARKERS_EXCLUDED;
    goto filter_attrib_ret_1;
  }
  LOGPRINTF("--attrib: %" PRIuPTR " variant%s remaining.\n", include_ct, (include_ct == 1)? "" : "s");
  memcpy(exclude_arr, exclude_arr_new, unfiltered_ctl * sizeof(intptr_t));
  *exclude_ct_ptr = unfiltered_ct - include_ct;
  while (0) {
  filter_attrib_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  filter_attrib_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  filter_attrib_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  filter_attrib_ret_INVALID_CMDLINE_2:
    logprintb();
  filter_attrib_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  filter_attrib_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 filter_attrib_ret_1:
  wkspace_reset(wkspace_mark);
  gzclose_cond(gz_infile);
  return retval;
}

int32_t filter_attrib_sample(char* fname, char* condition_str, char* sorted_ids, uintptr_t sorted_ids_ct, uintptr_t max_id_len, uint32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr) {
  // re-merge this with filter_attrib() after making sample ID lookup
  // hash-based
  gzFile gz_infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t include_ct = 0;
  uintptr_t unfiltered_ctl = (unfiltered_ct + (BITCT - 1)) / BITCT;
  uintptr_t* cur_neg_matches = NULL;
  char* sorted_pos_match = NULL;
  char* sorted_neg_match = NULL;
  char* id_buf = NULL;
  char* bufptr2 = NULL;
  uint32_t pos_match_ct = 0;
  uint32_t neg_match_ct = 0;
  uint32_t neg_match_ctl = 0;
  uintptr_t max_pos_match_len = 0;
  uintptr_t max_neg_match_len = 0;
  uintptr_t line_idx = 0;
  uint32_t is_neg = 0;
  int32_t retval = 0;
  uintptr_t* exclude_arr_new;
  uintptr_t* already_seen;
  char* loadbuf;
  char* cond_ptr;
  char* bufptr;
  uintptr_t loadbuf_size;
  uintptr_t pos_match_idx;
  uintptr_t neg_match_idx;
  uintptr_t ulii;
  uint32_t unfiltered_idx;
  uint32_t pos_match_needed;
  uint32_t cur_neg_match_ct;
  int32_t sorted_idx;
  
  if (wkspace_alloc_ul_checked(&exclude_arr_new, unfiltered_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&already_seen, unfiltered_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_c_checked(&id_buf, max_id_len)) { 
    goto filter_attrib_sample_ret_NOMEM;
  }
  fill_all_bits(exclude_arr_new, unfiltered_ct);
  fill_ulong_zero(already_seen, unfiltered_ctl);
  if (condition_str) {
    // allow NULL condition_str; this means all samples/variants named in the
    // file are included
    cond_ptr = condition_str;
    while (1) {
      while (*cond_ptr == ',') {
	cond_ptr++;
      }
      if (*cond_ptr == '-') {
	cond_ptr++;
	if (*cond_ptr == ',') {
	  continue;
	} else if (*cond_ptr == '-') {
	  logprint("Error: --attrib-indiv condition cannot contain consecutive dashes.\n");
	  goto filter_attrib_sample_ret_INVALID_CMDLINE;
	}
	is_neg = 1;
      } else if (!(*cond_ptr)) {
        // command-line parser actually comma-terminates the string for us, so
        // no need to check for null elsewhere
	break;
      }
      bufptr = strchr(cond_ptr, ',');
      ulii = (uintptr_t)(bufptr - cond_ptr);
      if (is_neg) {
	neg_match_ct++;
	if (ulii >= max_neg_match_len) {
	  max_neg_match_len = ulii + 1;
	}
      } else {
        pos_match_ct++;
	if (ulii >= max_pos_match_len) {
          max_pos_match_len = ulii + 1;
	}
      }
      cond_ptr = bufptr;
      is_neg = 0;
    }
    if (pos_match_ct) {
      if (wkspace_alloc_c_checked(&sorted_pos_match, max_pos_match_len * pos_match_ct)) {
	goto filter_attrib_sample_ret_NOMEM;
      }
    }
    if (neg_match_ct) {
      neg_match_ctl = (neg_match_ct + (BITCT - 1)) / BITCT;
      if (wkspace_alloc_c_checked(&sorted_neg_match, max_neg_match_len * neg_match_ct) ||
	  wkspace_alloc_ul_checked(&cur_neg_matches, neg_match_ctl * sizeof(intptr_t))) {
        goto filter_attrib_sample_ret_NOMEM;
      }
    }
    pos_match_idx = 0;
    neg_match_idx = 0;
    cond_ptr = condition_str;
    is_neg = 0;
    while (1) {
      while (*cond_ptr == ',') {
	cond_ptr++;
      }
      if (*cond_ptr == '-') {
	cond_ptr++;
	if (*cond_ptr == ',') {
	  continue;
	}
	is_neg = 1;
      } else if (!(*cond_ptr)) {
	break;
      }
      bufptr = strchr(cond_ptr, ',');
      ulii = (uintptr_t)(bufptr - cond_ptr);
      if (is_neg) {
	memcpyx(&(sorted_neg_match[neg_match_idx * max_neg_match_len]), cond_ptr, ulii, '\0');
	neg_match_idx++;
      } else {
	memcpyx(&(sorted_pos_match[pos_match_idx * max_pos_match_len]), cond_ptr, ulii, '\0');
        pos_match_idx++;
      }
      cond_ptr = bufptr;
      is_neg = 0;
    }
    if (pos_match_ct) {
      qsort(sorted_pos_match, pos_match_ct, max_pos_match_len, strcmp_casted);
      bufptr = scan_for_duplicate_ids(sorted_pos_match, pos_match_ct, max_pos_match_len);
      if (bufptr) {
	LOGPREPRINTFWW("Error: Duplicate attribute '%s' in --attrib-indiv argument.\n", bufptr);
	goto filter_attrib_sample_ret_INVALID_CMDLINE_2;
      }
    }
    if (neg_match_ct) {
      qsort(sorted_neg_match, neg_match_ct, max_neg_match_len, strcmp_casted);
      bufptr = scan_for_duplicate_ids(sorted_neg_match, neg_match_ct, max_neg_match_len);
      if (bufptr) {
	LOGPREPRINTFWW("Error: Duplicate attribute '%s' in --attrib-indiv argument.\n", bufptr);
	goto filter_attrib_sample_ret_INVALID_CMDLINE_2;
      }
      // actually may make sense to have same attribute as a positive and
      // negative condition, so we don't check for that
    }
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto filter_attrib_sample_ret_NOMEM;
  }
  if (gzopen_checked(&gz_infile, fname, "rb")) {
    goto filter_attrib_sample_ret_OPEN_FAIL;
  }
  if (gzbuffer(gz_infile, 131072)) {
    goto filter_attrib_sample_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (1) {
    line_idx++;
    if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
      if (!gzeof(gz_infile)) {
	goto filter_attrib_sample_ret_READ_FAIL;
      }
      break;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: Line %" PRIuPTR" of --attrib-indiv file is pathologically long.\n", line_idx);
        goto filter_attrib_sample_ret_INVALID_FORMAT_2;
      }
      goto filter_attrib_sample_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(id_buf, sorted_ids, max_id_len, sorted_ids_ct, bufptr, &cond_ptr, &sorted_idx)) {
      sprintf(logbuf, "Error: Line %" PRIuPTR " of --attrib-indiv file has fewer tokens than\nexpected.\n", line_idx);
      goto filter_attrib_sample_ret_INVALID_FORMAT_2;
    }
    if (sorted_idx == -1) {
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      *strchr(id_buf, '\t') = ' ';
      LOGPREPRINTFWW("Error: Duplicate sample ID '%s' in --attrib-indiv file.\n", id_buf);
      goto filter_attrib_sample_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    unfiltered_idx = id_map[(uint32_t)sorted_idx];
    pos_match_needed = pos_match_ct;
    cur_neg_match_ct = 0;
    if (neg_match_ct) {
      fill_ulong_zero(cur_neg_matches, neg_match_ctl);
    }
    while (!is_eoln_kns(*cond_ptr)) {
      bufptr2 = cond_ptr;
      bufptr = token_endnn(cond_ptr);
      ulii = (uintptr_t)(bufptr - bufptr2);
      cond_ptr = skip_initial_spaces(bufptr);
      if (pos_match_needed && (bsearch_str(bufptr2, ulii, sorted_pos_match, max_pos_match_len, pos_match_ct) != -1)) {
	pos_match_needed = 0;
      }
      if (cur_neg_match_ct < neg_match_ct) {
	sorted_idx = bsearch_str(bufptr2, ulii, sorted_neg_match, max_neg_match_len, neg_match_ct);
	if ((sorted_idx != -1) && (!is_set(cur_neg_matches, sorted_idx))) {
          cur_neg_match_ct++;
	  if (cur_neg_match_ct == neg_match_ct) {
	    // fail
	    pos_match_needed = 1;
            break;
	  }
          set_bit(cur_neg_matches, sorted_idx);
	}
      }
    }
    if (pos_match_needed) {
      continue;
    }
    // full negative match causes pos_match_needed to be set, so no further
    // check required
    clear_bit(exclude_arr_new, unfiltered_idx);
    include_ct++;
  }
  if (!include_ct) {
    LOGPRINTF("Error: No %s remaining after --attrib-indiv.\n", g_species_plural);
    retval = RET_ALL_SAMPLES_EXCLUDED;
    goto filter_attrib_sample_ret_1;
  }
  LOGPRINTF("--attrib-indiv: %" PRIuPTR " %s remaining.\n", include_ct, species_str(include_ct));
  memcpy(exclude_arr, exclude_arr_new, unfiltered_ctl * sizeof(intptr_t));
  *exclude_ct_ptr = unfiltered_ct - include_ct;
  while (0) {
  filter_attrib_sample_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  filter_attrib_sample_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  filter_attrib_sample_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  filter_attrib_sample_ret_INVALID_CMDLINE_2:
    logprintb();
  filter_attrib_sample_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  filter_attrib_sample_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 filter_attrib_sample_ret_1:
  wkspace_reset(wkspace_mark);
  gzclose_cond(gz_infile);
  return retval;
}

int32_t filter_qual_scores(Two_col_params* qual_filter, double qual_min_thresh, double qual_max_thresh, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t miss_ct = 0;
  uint32_t varid_first = (qual_filter->colid < qual_filter->colx);
  char skipchar = qual_filter->skipchar;
  uintptr_t* already_seen;
  uintptr_t* marker_exclude_orig;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uintptr_t slen;
  uintptr_t line_idx;
  uintptr_t marker_ct;
  char* colid_ptr; // variant ID
  char* colx_ptr; // quality score
  double dxx;
  uint32_t colmin;
  uint32_t coldiff;
  uint32_t marker_uidx;
  int32_t retval;
  char cc;
  if (wkspace_alloc_ul_checked(&already_seen, unfiltered_marker_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&marker_exclude_orig, unfiltered_marker_ctl * sizeof(intptr_t))) {
    goto filter_qual_scores_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, unfiltered_marker_ctl);
  memcpy(marker_exclude_orig, marker_exclude, unfiltered_marker_ctl * sizeof(intptr_t));

  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto filter_qual_scores_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, qual_filter->fname, loadbuf, loadbuf_size, qual_filter->skip);
  if (retval) {
    goto filter_qual_scores_ret_1;
  }
  if (varid_first) {
    colmin = qual_filter->colid - 1;
    coldiff = qual_filter->colx - qual_filter->colid;
  } else {
    colmin = qual_filter->colx - 1;
    coldiff = qual_filter->colid - qual_filter->colx;
  }
  line_idx = qual_filter->skip;
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        sprintf(logbuf, "Error: Line %" PRIuPTR " of --qual-scores file is pathologically long.\n", line_idx);
        goto filter_qual_scores_ret_INVALID_FORMAT_2;
      } else {
	goto filter_qual_scores_ret_NOMEM;
      }
    }
    colid_ptr = skip_initial_spaces(loadbuf);
    cc = *colid_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (varid_first) {
      colid_ptr = next_token_multz(colid_ptr, colmin);
      colx_ptr = next_token_mult(colid_ptr, coldiff);
      if (no_more_tokens_kns(colx_ptr)) {
        goto filter_qual_scores_ret_MISSING_TOKENS;
      }
    } else {
      colx_ptr = next_token_multz(colid_ptr, colmin);
      colid_ptr = next_token_mult(colx_ptr, coldiff);
      if (no_more_tokens_kns(colid_ptr)) {
        goto filter_qual_scores_ret_MISSING_TOKENS;
      }
    }
    slen = strlen_se(colid_ptr);
    marker_uidx = id_htable_find(colid_ptr, slen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if ((marker_uidx == 0xffffffffU) || is_set(marker_exclude_orig, marker_uidx)) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, marker_uidx)) {
      colid_ptr[slen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant '%s' in --qual-scores file.\n", colid_ptr);
      goto filter_qual_scores_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, marker_uidx);
    if (scan_double(colx_ptr, &dxx) || (dxx < qual_min_thresh) || (dxx > qual_max_thresh)) {
      set_bit(marker_exclude, marker_uidx);
    }
  }
  if (!feof(infile)) {
    goto filter_qual_scores_ret_READ_FAIL;
  }
  *marker_exclude_ct_ptr = popcount_longs(marker_exclude, unfiltered_marker_ctl);
  marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  if (miss_ct) {
    sprintf(logbuf, "--qual-scores: %" PRIuPTR " variant%s remaining, %" PRIuPTR " ID%s missing.\n", marker_ct, (marker_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--qual-scores: %" PRIuPTR " variant%s remaining.\n", marker_ct, (marker_ct == 1)? "" : "s");
  }
  logprintb();
  while (0) {
  filter_qual_scores_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  filter_qual_scores_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  filter_qual_scores_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --qual-scores file has fewer tokens than expected.\n", line_idx);
  filter_qual_scores_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 filter_qual_scores_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  return retval;
}

uint32_t random_thin_markers(double thin_keep_prob, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr) {
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t removed_ct = 0;
  uint32_t uint32_thresh = (uint32_t)(thin_keep_prob * 4294967296.0 + 0.5);
  uint32_t marker_uidx_stop;
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    do {
      if (sfmt_genrand_uint32(&sfmt) >= uint32_thresh) {
	SET_BIT(marker_exclude, marker_uidx);
	removed_ct++;
      }
    } while (++marker_uidx < marker_uidx_stop);
  }
  if (marker_ct == removed_ct) {
    logprint("Error: All variants removed by --thin.  Try a higher probability.\n");
    return 1;
  }
  LOGPRINTF("--thin: %u variant%s removed (%u remaining).\n", removed_ct, (removed_ct == 1)? "" : "s", marker_ct - removed_ct);
  *marker_exclude_ct_ptr += removed_ct;
  return 0;
}

int32_t random_thin_markers_ct(uint32_t thin_keep_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  int32_t retval = 0;
  uintptr_t* perm_buf;
  uint32_t marker_idx;
  if (thin_keep_ct > marker_ct) {
    LOGPRINTF("Error: --thin-count parameter exceeds number of remaining variants.\n");
    goto random_thin_markers_ct_ret_INVALID_CMDLINE;
  }
  if (wkspace_alloc_ul_checked(&perm_buf, marker_ctl * sizeof(intptr_t))) {
    goto random_thin_markers_ct_ret_NOMEM;
  }
  // no actual interleaving here, but may as well use this function
  generate_perm1_interleaved(marker_ct, marker_ct - thin_keep_ct, 0, 1, perm_buf);
  marker_uidx = 0;
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    if (is_set(perm_buf, marker_idx)) {
      set_bit(marker_exclude, marker_uidx);
    }
  }
  LOGPRINTF("--thin-count: %u variant%s removed (%u remaining).\n", marker_ct - thin_keep_ct, (marker_ct - thin_keep_ct == 1)? "" : "s", thin_keep_ct);
  *marker_exclude_ct_ptr = unfiltered_marker_ct - thin_keep_ct;
  while (0) {
  random_thin_markers_ct_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  random_thin_markers_ct_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  wkspace_reset(wkspace_mark);
  return retval;
}

uint32_t random_thin_samples(double thin_keep_prob, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr) {
  uint32_t sample_ct = unfiltered_sample_ct - *sample_exclude_ct_ptr;
  uint32_t sample_uidx = 0;
  uint32_t samples_done = 0;
  uint32_t removed_ct = 0;
  uint32_t uint32_thresh = (uint32_t)(thin_keep_prob * 4294967296.0 + 0.5);
  uint32_t sample_uidx_stop;
  while (samples_done < sample_ct) {
    sample_uidx = next_unset_unsafe(sample_exclude, sample_uidx);
    sample_uidx_stop = next_set(sample_exclude, sample_uidx, unfiltered_sample_ct);
    samples_done += sample_uidx_stop - sample_uidx;
    do {
      if(sfmt_genrand_uint32(&sfmt) >= uint32_thresh) {
        SET_BIT(sample_exclude, sample_uidx);
        removed_ct++;
      }
    } while (++sample_uidx < sample_uidx_stop);
  }
  if (sample_ct == removed_ct) {
    LOGPRINTF("Error: All %s removed by --thin-indiv. Try a higher probability.\n", g_species_plural);
    return 1;
  }
  LOGPRINTF("--thin-indiv: %u %s removed (%u remaining).\n", removed_ct, (removed_ct==1)? g_species_singular : g_species_plural, sample_ct - removed_ct);
  *sample_exclude_ct_ptr += removed_ct;
  return 0;
}

int32_t random_thin_samples_ct(uint32_t thin_keep_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t sample_ct = unfiltered_sample_ct - *sample_exclude_ct_ptr;
  uint32_t sample_uidx = 0;
  uintptr_t sample_ctl = (sample_ct + (BITCT - 1)) / BITCT;
  int32_t retval = 0;
  uintptr_t* perm_buf;
  uint32_t sample_idx;
  if (thin_keep_ct > sample_ct) {
    LOGPRINTF("Error: --thin-indiv-count parameter exceeds number of remaining %s.\n", g_species_plural);
    goto random_thin_samples_ct_ret_INVALID_CMDLINE;
  }
  if (wkspace_alloc_ul_checked(&perm_buf, sample_ctl * sizeof(intptr_t))) {
    goto random_thin_samples_ct_ret_NOMEM;
  }

  generate_perm1_interleaved(sample_ct, sample_ct - thin_keep_ct, 0, 1, perm_buf);
  sample_uidx = 0;
  for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
    next_unset_unsafe_ck(sample_exclude, &sample_uidx);
    if (is_set(perm_buf, sample_idx)) {
      set_bit(sample_exclude, sample_uidx);
    }
  }
  LOGPRINTF("--thin-indiv-count: %u %s removed (%u remaining).\n", sample_ct - thin_keep_ct, (sample_ct - thin_keep_ct == 1)? g_species_singular : g_species_plural, thin_keep_ct);
  *sample_exclude_ct_ptr = unfiltered_sample_ct - thin_keep_ct;
  while(0) {
  random_thin_samples_ct_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  random_thin_samples_ct_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  wkspace_reset(wkspace_mark);
  return retval;
}


int32_t load_oblig_missing(FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, char* sorted_sample_ids, uintptr_t sorted_sample_ct, uintptr_t max_sample_id_len, uint32_t* sample_id_map, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sex_male, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip) {
  // 1. load and validate cluster file
  // 2. load marker file, sort by uidx
  // 3. check for early exit (no clusters and/or no .zero entries)
  // 4. scan through .bed sequentially, update oblig_missing_..._cts
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char* idbuf = &(tbuf[MAXLINELEN]);
  Ll_str* cluster_names = NULL;
  uint64_t tot_missing = 0;
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = (unfiltered_sample_ct + BITCT2 - 1) / BITCT2;
  uintptr_t sorted_sample_ctl = (sorted_sample_ct + BITCT - 1) / BITCT;
  uintptr_t topsize = 0;
  uintptr_t max_cluster_id_len = 0;
  uintptr_t possible_distinct_ct = 0;
  uintptr_t missing_cluster_ct = 0;
  uintptr_t y_start = 0;
  uintptr_t y_end = 0;
  uintptr_t line_idx = 0;
  int32_t y_code = chrom_info_ptr->y_code;
  uint32_t y_present = ((y_code != -1) && is_set(chrom_info_ptr->chrom_mask, y_code));
  int32_t retval = 0;
  Ll_str* llptr;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_end;
  uintptr_t* cluster_zmask2s;
  uintptr_t* loadbuf_ptr;
  uintptr_t* cur_cluster_zmask2;
  uintptr_t* ulptr;
  uint32_t* cluster_sizes;
  uint32_t* marker_id_map;
  uint32_t* cluster_ref_cts;
  uint32_t* sample_lookup;
  char* cluster_ids;
  char* sorted_marker_ids;
  char* bufptr;
  char* bufptr2;
  int64_t* zc_entries;
  int64_t* zc_entries_end;
  int64_t* wkspace_end;
  uintptr_t cluster_ct;
  uintptr_t cluster_mct; // doubled if Y chrom present
  uintptr_t marker_uidx;
  uintptr_t ulii;
  uint64_t ullii;
  uint32_t sample_uidx;
  uint32_t cluster_idx;
  uint32_t slen;
  int32_t ii;

  if (y_present) {
    y_start = chrom_info_ptr->chrom_start[(uint32_t)y_code];
    y_end = chrom_info_ptr->chrom_end[(uint32_t)y_code];
  }
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_sample_ctl2 * sizeof(intptr_t))) {
    goto load_oblig_missing_ret_NOMEM;
  }
  loadbuf_end = &(loadbuf[unfiltered_sample_ctl2]);
  if (fopen_checked(&infile, om_ip->sample_fname, "r")) {
    goto load_oblig_missing_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';

  // two-pass load, same as load_clusters()
  // use loadbuf as duplicate IID detector
  fill_ulong_zero(loadbuf, sorted_sample_ctl);
  while (fgets(tbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, om_ip->sample_fname);
      goto load_oblig_missing_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(idbuf, sorted_sample_ids, max_sample_id_len, sorted_sample_ct, bufptr, &bufptr2, &ii)) {
      LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, om_ip->sample_fname);
      goto load_oblig_missing_ret_INVALID_FORMAT_2;
    }
    if (ii != -1) {
      if (is_set(loadbuf, ii)) {
        strchr(idbuf, '\t')[0] = ' ';
        LOGPREPRINTFWW("Error: Duplicate sample ID '%s' in %s.\n", idbuf, om_ip->sample_fname);
	goto load_oblig_missing_ret_INVALID_FORMAT_2;
      }
      set_bit(loadbuf, ii);
      slen = strlen_se(bufptr2);
      if (slen >= max_cluster_id_len) {
	max_cluster_id_len = slen + 1;
      }
      bufptr2[slen] = '\0';
      if ((!cluster_names) || (strcmp(cluster_names->ss, bufptr2) && ((!cluster_names->next) || strcmp(cluster_names->next->ss, bufptr2)))) {
	llptr = top_alloc_llstr(&topsize, slen + 1);
	if (!llptr) {
	  goto load_oblig_missing_ret_NOMEM;
	}
	llptr->next = cluster_names;
	memcpy(llptr->ss, bufptr2, slen + 1);
	cluster_names = llptr;
	possible_distinct_ct++;
      }
    }
  }
  if (!feof(infile)) {
    goto load_oblig_missing_ret_READ_FAIL;
  }
  if (!max_cluster_id_len) {
    LOGPRINTFWW("Warning: --oblig-missing ignored, since no valid blocks were defined in %s.\n", om_ip->sample_fname);
    goto load_oblig_missing_ret_1;
  }
  wkspace_left -= topsize;
  if (wkspace_alloc_c_checked(&cluster_ids, possible_distinct_ct * max_cluster_id_len)) {
    goto load_oblig_missing_ret_NOMEM2;
  }
  for (ulii = 0; ulii < possible_distinct_ct; ulii++) {
    strcpy(&(cluster_ids[ulii * max_cluster_id_len]), cluster_names->ss);
    cluster_names = cluster_names->next;
  }
  wkspace_left += topsize;
  topsize = 0;
  qsort(cluster_ids, possible_distinct_ct, max_cluster_id_len, strcmp_casted);
  cluster_ct = collapse_duplicate_ids(cluster_ids, possible_distinct_ct, max_cluster_id_len, NULL);
  wkspace_shrink_top(cluster_ids, cluster_ct * max_cluster_id_len);
  cluster_mct = cluster_ct * (y_present + 1);
  sample_lookup = (uint32_t*)malloc(unfiltered_sample_ct * sizeof(int32_t));
  if (!sample_lookup) {
    goto load_oblig_missing_ret_NOMEM;
  }
  om_ip->sample_lookup = sample_lookup;
  if (wkspace_alloc_ui_checked(&cluster_sizes, cluster_mct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&cluster_zmask2s, cluster_mct * unfiltered_sample_ctl2 * sizeof(intptr_t))) {
    goto load_oblig_missing_ret_NOMEM;
  }
  fill_uint_zero(cluster_sizes, cluster_mct);
  fill_uint_one(sample_lookup, unfiltered_sample_ct);
  fill_ulong_zero(cluster_zmask2s, cluster_mct * unfiltered_sample_ctl2);

  // second pass
  rewind(infile);
  while (fgets(tbuf, MAXLINELEN, infile)) {
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bsearch_read_fam_indiv(idbuf, sorted_sample_ids, max_sample_id_len, sorted_sample_ct, bufptr, &bufptr2, &ii);
    if (ii == -1) {
      continue;
    }
    sample_uidx = sample_id_map[(uint32_t)ii];    
    slen = strlen_se(bufptr2);
    // guaranteed to succeed
    ii = bsearch_str(bufptr2, slen, cluster_ids, max_cluster_id_len, cluster_ct);
    set_bit(&(cluster_zmask2s[((uintptr_t)((uint32_t)ii)) * unfiltered_sample_ctl2]), sample_uidx * 2);
    cluster_sizes[(uint32_t)ii] += 1;
    sample_lookup[sample_uidx] = (uint32_t)ii;
  }
  if (fclose_null(&infile)) {
    goto load_oblig_missing_ret_READ_FAIL;
  }
  if (y_present) {
    vec_include_init(unfiltered_sample_ct, loadbuf, sex_male);
    cur_cluster_zmask2 = cluster_zmask2s;
    ulptr = &(cur_cluster_zmask2[cluster_ct * unfiltered_sample_ctl2]);
    for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
      slen = 0;
      for (loadbuf_ptr = loadbuf; loadbuf_ptr < loadbuf_end;) {
        ulii = (*loadbuf_ptr++) & (*cur_cluster_zmask2++);
	slen += popcount2_long(ulii);
	*ulptr++ = ulii;
      }
      cluster_sizes[cluster_idx + cluster_ct] = slen;
    }
  }

  cluster_ref_cts = (uint32_t*)malloc(cluster_ct * 2 * sizeof(int32_t));
  if (!cluster_ref_cts) {
    goto load_oblig_missing_ret_NOMEM;
  }
  om_ip->cluster_ref_cts = cluster_ref_cts;
  fill_uint_zero(cluster_ref_cts, cluster_ct * 2);
  retval = sort_item_ids(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    goto load_oblig_missing_ret_1;
  }
  zc_entries = (int64_t*)wkspace_base;
  zc_entries_end = zc_entries;
  wkspace_end = (int64_t*)(&(wkspace_base[wkspace_left]));
  if (fopen_checked(&infile, om_ip->marker_fname, "r")) {
    goto load_oblig_missing_ret_OPEN_FAIL;
  }
  line_idx = 0;
  while (fgets(tbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, om_ip->marker_fname);
      goto load_oblig_missing_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr2 = token_endnn(bufptr);
    ii = bsearch_str(bufptr, (uintptr_t)(bufptr2 - bufptr), sorted_marker_ids, max_marker_id_len, marker_ct);
    if (ii != -1) {
      marker_uidx = marker_id_map[(uint32_t)ii];
      bufptr = skip_initial_spaces(bufptr2);
      if (is_eoln_kns(*bufptr)) {
        LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, om_ip->marker_fname);
        goto load_oblig_missing_ret_INVALID_FORMAT_2;
      }
      slen = strlen_se(bufptr);
      ii = bsearch_str(bufptr, slen, cluster_ids, max_cluster_id_len, cluster_ct);
      if (ii != -1) {
	if (zc_entries_end == wkspace_end) {
          goto load_oblig_missing_ret_NOMEM;
	}
	cluster_idx = (uint32_t)ii;
	if ((marker_uidx < y_end) && (marker_uidx >= y_start)) {
          cluster_idx += cluster_ct;
	}
        *zc_entries_end++ = (((uint64_t)marker_uidx) << 32) | ((uint64_t)cluster_idx);
	cluster_ref_cts[cluster_idx] += 1;
      } else {
        missing_cluster_ct++;
      }
    }
  }
  if (fclose_null(&infile)) {
    goto load_oblig_missing_ret_READ_FAIL;
  }
  if (missing_cluster_ct) {
    LOGPRINTFWW("Warning: %" PRIuPTR " entr%s in %s had block IDs missing from %s.\n", missing_cluster_ct, (missing_cluster_ct == 1)? "y" : "ies", om_ip->marker_fname, om_ip->sample_fname);
  }
  om_ip->entry_ct = (uintptr_t)(zc_entries_end - zc_entries);
  if (!om_ip->entry_ct) {
    LOGPRINTFWW("Warning: --oblig-missing ignored, since %s had no valid entries.\n", om_ip->marker_fname);
    goto load_oblig_missing_ret_1;
  }

#ifdef __cplusplus
  std::sort(zc_entries, zc_entries_end);
#else
  qsort(zc_entries, om_ip->entry_ct, sizeof(int64_t), llcmp);
#endif
  om_ip->entries = (uint64_t*)malloc((om_ip->entry_ct + 1) * sizeof(int64_t));
  if (!om_ip->entries) {
    goto load_oblig_missing_ret_NOMEM;
  }
  memcpy(om_ip->entries, zc_entries, om_ip->entry_ct * sizeof(int64_t));
  om_ip->entries[om_ip->entry_ct] = ((uint64_t)unfiltered_marker_ct) << 32;
  om_ip->cluster_ct = cluster_ct;

  loadbuf[unfiltered_sample_ctl2 - 1] = 0;
  marker_uidx = unfiltered_marker_ct; // forces initial load
  do {
    ullii = *zc_entries++;
    if ((ullii >> 32) != marker_uidx) {
      marker_uidx = (ullii >> 32);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	goto load_oblig_missing_ret_READ_FAIL;
      }
      if (load_raw(bedfile, loadbuf, unfiltered_sample_ct4)) {
	goto load_oblig_missing_ret_READ_FAIL;
      }
      // no need for het haploid handling here
    }
    loadbuf_ptr = loadbuf;
    cluster_idx = (uint32_t)ullii;
    cur_cluster_zmask2 = &(cluster_zmask2s[cluster_idx * unfiltered_sample_ctl2]);
    do {
      ulii = *cur_cluster_zmask2++;
      if (((*loadbuf_ptr) ^ ulii) & (ulii * 3)) {
	LOGPREPRINTFWW("Error: Nonmissing --oblig-missing genotype at variant '%s'. (To force it to missing, use --zero-cluster.)\n", &(marker_ids[marker_uidx * max_marker_id_len]));
        goto load_oblig_missing_ret_INVALID_FORMAT_2;
      }
    } while ((++loadbuf_ptr) < loadbuf_end);
    tot_missing += cluster_sizes[cluster_idx];
  } while (zc_entries < zc_entries_end);
  LOGPRINTF("--oblig-missing: %" PRIu64 " call%s confirmed missing.\n", tot_missing, (tot_missing == 1)? "" : "s");
  while (0) {
  load_oblig_missing_ret_NOMEM2:
    wkspace_left += topsize;
  load_oblig_missing_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_oblig_missing_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_oblig_missing_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_oblig_missing_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 load_oblig_missing_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  return retval;
}

int32_t filter_samples_file(char* filtername, char* sorted_sample_ids, uintptr_t sorted_ids_len, uintptr_t max_sample_id_len, uint32_t* id_map, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, char* filtervals_flattened, uint32_t mfilter_col) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t include_ct = 0;
  uintptr_t max_filterval_len = 0;
  uintptr_t line_idx = 0;
  uint32_t filterval_ct = count_and_measure_multistr(filtervals_flattened, &max_filterval_len);
  int32_t retval = 0;
  char* sorted_filtervals;
  uintptr_t* sample_exclude_new;
  char* id_buf;
  char* bufptr;
  uint32_t filterval_idx;
  uint32_t slen;
  int32_t sample_idx;
  if (wkspace_alloc_c_checked(&id_buf, max_sample_id_len) ||
      wkspace_alloc_ul_checked(&sample_exclude_new, unfiltered_sample_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_c_checked(&sorted_filtervals, filterval_ct * max_filterval_len)) {
    goto filter_samples_file_ret_NOMEM;
  }
  fill_all_bits(sample_exclude_new, unfiltered_sample_ct);
  bufptr = filtervals_flattened;
  for (filterval_idx = 0; filterval_idx < filterval_ct; filterval_idx++) {
    slen = strlen(bufptr) + 1;
    memcpy(&(sorted_filtervals[filterval_idx * max_filterval_len]), bufptr, slen);
    bufptr = &(bufptr[slen]);
  }
  qsort(sorted_filtervals, filterval_ct, max_filterval_len, strcmp_casted);

  if (fopen_checked(&infile, filtername, "r")) {
    goto filter_samples_file_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Line %" PRIuPTR " of --filter file is pathologically long.\n", line_idx);
      goto filter_samples_file_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(id_buf, sorted_sample_ids, max_sample_id_len, sorted_ids_len, bufptr, &bufptr, &sample_idx)) {
      goto filter_samples_file_ret_MISSING_TOKENS;
    }
    if (sample_idx != -1) {
      sample_idx = id_map[(uint32_t)sample_idx];
      if (!is_set(sample_exclude, sample_idx)) {
	if (mfilter_col > 1) {
	  bufptr = next_token_mult(bufptr, mfilter_col - 1);
	}
	if (no_more_tokens_kns(bufptr)) {
	  goto filter_samples_file_ret_MISSING_TOKENS;
	}
	if (bsearch_str(bufptr, strlen_se(bufptr), sorted_filtervals, max_filterval_len, filterval_ct) != -1) {
	  if (is_set(sample_exclude_new, sample_idx)) {
	    clear_bit(sample_exclude_new, sample_idx);
	    include_ct++;
	  }
	}
      }
    }
  }
  if (!feof(infile)) {
    goto filter_samples_file_ret_READ_FAIL;
  }
  if (!include_ct) {
    LOGPRINTF("Error: All %s excluded by --filter.\n", g_species_plural);
    goto filter_samples_file_ret_ALL_SAMPLES_EXCLUDED;
  }
  LOGPRINTF("--filter: %" PRIuPTR " %s remaining.\n", include_ct, species_str(include_ct));
  memcpy(sample_exclude, sample_exclude_new, unfiltered_sample_ctl * sizeof(intptr_t));
  *sample_exclude_ct_ptr = unfiltered_sample_ct - include_ct;

  while (0) {
  filter_samples_file_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  filter_samples_file_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  filter_samples_file_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  filter_samples_file_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --filter file has fewer tokens than expected.\n", line_idx);
  filter_samples_file_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  filter_samples_file_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  return retval;
}

void filter_samples_bitfields(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, uintptr_t* orfield, int32_t orfield_flip, uintptr_t* ornot) {
  // sample_exclude := sample_exclude | orfield | (~ornot) if !orfield_flip
  //                := sample_exclude | (~orfield) | (~ornot) otherwise
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t* ieptr = sample_exclude;
  uintptr_t* ieend = &(sample_exclude[unfiltered_sample_ctl]);
  if (orfield_flip) {
    if (ornot) {
      do {
	*ieptr |= (~(*orfield++)) | (~(*ornot++));
      } while (++ieptr < ieend);
    } else {
      do {
	*ieptr |= ~(*orfield++);
      } while (++ieptr < ieend);
    }
  } else {
    if (ornot) {
      do {
	*ieptr |= (*orfield++) | (~(*ornot++));
      } while (++ieptr < ieend);
    } else {
      do {
	*ieptr |= *orfield++;
      } while (++ieptr < ieend);
    }
  }
  zero_trailing_bits(sample_exclude, unfiltered_sample_ct);
  *sample_exclude_ct_ptr = popcount_longs(sample_exclude, unfiltered_sample_ctl);
}

int32_t mind_filter(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double mind_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* sex_male, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uint32_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ct2l = (unfiltered_sample_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
  uintptr_t marker_idx = 0;
  uintptr_t y_start = 0;
  uintptr_t y_end = 0;
  uintptr_t* sample_male_include2 = NULL;
  uint32_t unfiltered_sample_ctl2m1 = (unfiltered_sample_ct - 1) / BITCT2;
  uint32_t sample_exclude_ct = *sample_exclude_ct_ptr;
  uint32_t sample_ct = unfiltered_sample_ct - sample_exclude_ct;
  uint32_t sample_uidx = 0;
  uint32_t sample_idx = 0;
  uint32_t removed_ct = 0;
  int32_t y_code = chrom_info_ptr->y_code;
  uint32_t y_present = (y_code != -1) && is_set(chrom_info_ptr->chrom_mask, y_code);
  uint32_t nony_marker_ct = marker_ct;
  int32_t retval = 0;
  uint32_t mind_int_thresh[2];
  uintptr_t* loadbuf;
  uintptr_t* newly_excluded;
  uintptr_t* lptr;
  uintptr_t* mptr;
  uint32_t* missing_cts;
  uint32_t* cluster_ref_cts;
  uint32_t* sample_lookup;
  uintptr_t marker_uidx;
  uint32_t cluster_ct;
  uint32_t cur_marker_ct;
  uint32_t sample_uidx_stop;
  uint32_t uii;
  uint32_t ujj;
  uintptr_t ulii;

  if (y_present) {
    y_start = chrom_info_ptr->chrom_start[(uint32_t)y_code];
    y_end = chrom_info_ptr->chrom_end[(uint32_t)y_code];
    if (wkspace_alloc_ul_checked(&sample_male_include2, unfiltered_sample_ct2l * sizeof(intptr_t))) {
      goto mind_filter_ret_NOMEM;
    }
    vec_include_init(unfiltered_sample_ct, sample_male_include2, sex_male);
    nony_marker_ct = marker_ct - (y_end - y_start - popcount_bit_idx(marker_exclude, y_start, y_end));
  }
  if (wkspace_alloc_ui_checked(&missing_cts, unfiltered_sample_ct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, unfiltered_sample_ct2l * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&newly_excluded, unfiltered_sample_ctl * sizeof(int32_t))) {
    goto mind_filter_ret_NOMEM;
  }
  loadbuf[unfiltered_sample_ct2l - 1] = 0;
  fill_uint_zero(missing_cts, unfiltered_sample_ct);
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto mind_filter_ret_READ_FAIL;
  }
  ujj = unfiltered_sample_ct2l * BITCT2;
  marker_uidx = 0;
  for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	goto mind_filter_ret_READ_FAIL;
      }
    }
    if (load_raw2(bedfile, loadbuf, unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask)) {
      goto mind_filter_ret_READ_FAIL;
    }
    // er, why doesn't this use load_and_collapse?
    lptr = loadbuf;
    if ((marker_uidx >= y_end) || (marker_uidx < y_start)) {
      for (uii = 0; uii < ujj; uii += BITCT2) {
	ulii = *lptr++;
	ulii = (ulii & FIVEMASK) & ((~ulii) >> 1);
	// now ulii has single bit set only at missing positions
	while (ulii) {
	  missing_cts[uii + CTZLU(ulii) / 2] += 1;
	  ulii &= ulii - 1;
	}
      }
    } else {
      mptr = sample_male_include2;
      for (uii = 0; uii < ujj; uii += BITCT2) {
	ulii = *lptr++;
	ulii = (ulii & (*mptr++)) & ((~ulii) >> 1);
	while (ulii) {
	  missing_cts[uii + CTZLU(ulii) / 2] += 1;
	  ulii &= ulii - 1;
	}
      }
    }
  }
  fill_ulong_zero(newly_excluded, unfiltered_sample_ctl);
  if (!om_ip->entry_ct) {
    mind_int_thresh[0] = (int32_t)(mind_thresh * ((int32_t)nony_marker_ct) * (1 + SMALL_EPSILON));
    mind_int_thresh[1] = (int32_t)(mind_thresh * ((int32_t)marker_ct) * (1 + SMALL_EPSILON));
    do {
      sample_uidx = next_unset_unsafe(sample_exclude, sample_uidx);
      sample_uidx_stop = next_set(sample_exclude, sample_uidx, unfiltered_sample_ct);
      sample_idx += sample_uidx_stop - sample_uidx;
      do {
	if (missing_cts[sample_uidx] > mind_int_thresh[is_set(sex_male, sample_uidx)]) {
	  SET_BIT(newly_excluded, sample_uidx);
	  removed_ct++;
	}
      } while (++sample_uidx < sample_uidx_stop);
    } while (sample_idx < sample_ct);
  } else {
    sample_lookup = om_ip->sample_lookup;
    cluster_ref_cts = om_ip->cluster_ref_cts;
    cluster_ct = om_ip->cluster_ct;
    do {
      sample_uidx = next_unset_unsafe(sample_exclude, sample_uidx);
      sample_uidx_stop = next_set(sample_exclude, sample_uidx, unfiltered_sample_ct);
      sample_idx += sample_uidx_stop - sample_uidx;
      do {
	uii = sample_lookup[sample_uidx];
	ujj = 0;
	if (is_set(sex_male, sample_uidx)) {
	  cur_marker_ct = marker_ct;
	  if (uii != 0xffffffffU) {
	    ujj = cluster_ref_cts[uii] + cluster_ref_cts[uii + cluster_ct];
	  }
	} else {
	  cur_marker_ct = nony_marker_ct;
	  if (uii != 0xffffffffU) {
	    ujj = cluster_ref_cts[uii];
	  }
	}
	if ((missing_cts[sample_uidx] - ujj) > (uint32_t)((int32_t)(mind_thresh * ((int32_t)(cur_marker_ct - ujj)) * (1 + SMALL_EPSILON)))) {
	  SET_BIT(newly_excluded, sample_uidx);
	  removed_ct++;
	}
      } while (++sample_uidx < sample_uidx_stop);
    } while (sample_idx < sample_ct);
  }
  if (removed_ct) {
    bitfield_or(sample_exclude, newly_excluded, unfiltered_sample_ctl);
    memcpy(outname_end, ".irem", 6);
    if (fopen_checked(&outfile, outname, "w")) {
      goto mind_filter_ret_OPEN_FAIL;
    }
    sample_uidx = 0;
    for (sample_idx = 0; sample_idx < removed_ct; sample_idx++, sample_uidx++) {
      next_set_unsafe_ck(newly_excluded, &sample_uidx);
      fputs(&(sample_ids[sample_uidx * max_sample_id_len]), outfile);
      putc('\n', outfile);
    }
    if (fclose_null(&outfile)) {
      goto mind_filter_ret_WRITE_FAIL;
    }
  }
  *sample_exclude_ct_ptr += removed_ct;
  if (*sample_exclude_ct_ptr == unfiltered_sample_ct) {
    LOGPRINTF("Error: All %s removed due to missing genotype data (--mind).\n", g_species_plural);
    LOGPRINTFWW("IDs written to %s .\n", outname);
    goto mind_filter_ret_ALL_SAMPLES_EXCLUDED;
  }
  LOGPRINTF("%u %s removed due to missing genotype data (--mind).\n", removed_ct, species_str(removed_ct));
  if (removed_ct) {
    LOGPRINTFWW("ID%s written to %s .\n", (removed_ct == 1)? "" : "s", outname);
  }
  while (0) {
  mind_filter_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  mind_filter_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  mind_filter_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  mind_filter_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  mind_filter_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

#ifdef __LP64__
void freq_hwe_haploid_count_120v(__m128i* vptr, __m128i* vend, __m128i* maskvp, uint32_t* ct_nmp, uint32_t* ct_hmajp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i to_ct_nm1;
  __m128i to_ct_hmaj1;
  __m128i to_ct_nm2;
  __m128i to_ct_hmaj2;
  __uni16 acc_nm;
  __uni16 acc_hmaj;

  acc_nm.vi = _mm_setzero_si128();
  acc_hmaj.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_andnot_si128(loader2, loader3);
    to_ct_hmaj1 = _mm_and_si128(loader, loader3);

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_and_si128(loader, loader3));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_and_si128(loader, loader3));

    to_ct_nm1 = _mm_add_epi64(_mm_and_si128(to_ct_nm1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_nm1, 2), m2));
    to_ct_hmaj1 = _mm_add_epi64(_mm_and_si128(to_ct_hmaj1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj1, 2), m2));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_andnot_si128(loader2, loader3);
    to_ct_hmaj2 = _mm_and_si128(loader, loader3);

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_add_epi64(to_ct_nm2, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj2 = _mm_add_epi64(to_ct_hmaj2, _mm_and_si128(loader, loader3));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_add_epi64(to_ct_nm2, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj2 = _mm_add_epi64(to_ct_hmaj2, _mm_and_si128(loader, loader3));

    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_add_epi64(_mm_and_si128(to_ct_nm2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_nm2, 2), m2)));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_add_epi64(_mm_and_si128(to_ct_hmaj2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj2, 2), m2)));

    acc_nm.vi = _mm_add_epi64(acc_nm.vi, _mm_add_epi64(_mm_and_si128(to_ct_nm1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_nm1, 4), m4)));
    acc_hmaj.vi = _mm_add_epi64(acc_hmaj.vi, _mm_add_epi64(_mm_and_si128(to_ct_hmaj1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj1, 4), m4)));
  } while (vptr < vend);
  acc_nm.vi = _mm_add_epi64(_mm_and_si128(acc_nm.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_nm.vi, 8), m8));
  acc_hmaj.vi = _mm_add_epi64(_mm_and_si128(acc_hmaj.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_hmaj.vi, 8), m8));
  *ct_nmp += ((acc_nm.u8[0] + acc_nm.u8[1]) * 0x1000100010001LLU) >> 48;
  *ct_hmajp += ((acc_hmaj.u8[0] + acc_hmaj.u8[1]) * 0x1000100010001LLU) >> 48;
}
#else
void freq_hwe_haploid_count_12(uintptr_t* lptr, uintptr_t* maskp, uint32_t* ct_nmp, uint32_t* ct_hmajp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader3 = loader >> 1;
  uintptr_t loader2 = loader ^ (~loader3);
  uint32_t to_ct_nm1;
  uint32_t to_ct_hmaj1;
  uint32_t to_ct_nm2;
  uint32_t to_ct_hmaj2;
  uintptr_t partial_nm;
  uintptr_t partial_hmaj;
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 = loader2 & loader3;
  to_ct_hmaj1 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  to_ct_nm1 = (to_ct_nm1 & 0x33333333) + ((to_ct_nm1 >> 2) & 0x33333333);
  to_ct_nm1 += (to_ct_nm2 & 0x33333333) + ((to_ct_nm2 >> 2) & 0x33333333);
  partial_nm = (to_ct_nm1 & 0x0f0f0f0f) + ((to_ct_nm1 >> 4) & 0x0f0f0f0f);
  to_ct_hmaj1 = (to_ct_hmaj1 & 0x33333333) + ((to_ct_hmaj1 >> 2) & 0x33333333);
  to_ct_hmaj1 += (to_ct_hmaj2 & 0x33333333) + ((to_ct_hmaj2 >> 2) & 0x33333333);
  partial_hmaj = (to_ct_hmaj1 & 0x0f0f0f0f) + ((to_ct_hmaj1 >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 = loader2 & loader3;
  to_ct_hmaj1 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr;
  loader3 = loader >> 1;
  loader2 = loader ^ (~loader3);
  loader &= loader3;
  loader3 = *maskp;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  to_ct_nm1 = (to_ct_nm1 & 0x33333333) + ((to_ct_nm1 >> 2) & 0x33333333);
  to_ct_nm1 += (to_ct_nm2 & 0x33333333) + ((to_ct_nm2 >> 2) & 0x33333333);
  partial_nm += (to_ct_nm1 & 0x0f0f0f0f) + ((to_ct_nm1 >> 4) & 0x0f0f0f0f);
  to_ct_hmaj1 = (to_ct_hmaj1 & 0x33333333) + ((to_ct_hmaj1 >> 2) & 0x33333333);
  to_ct_hmaj1 += (to_ct_hmaj2 & 0x33333333) + ((to_ct_hmaj2 >> 2) & 0x33333333);
  partial_hmaj += (to_ct_hmaj1 & 0x0f0f0f0f) + ((to_ct_hmaj1 >> 4) & 0x0f0f0f0f);

  *ct_nmp += (partial_nm * 0x01010101) >> 24;
  *ct_hmajp += (partial_hmaj * 0x01010101) >> 24;
}
#endif

static inline void single_marker_freqs_and_hwe(uintptr_t unfiltered_sample_ctl2, uintptr_t* lptr, uintptr_t* sample_include2, uintptr_t* founder_include2, uintptr_t* founder_ctrl_include2, uintptr_t* founder_case_include2, uintptr_t sample_ct, uint32_t* ll_ctp, uint32_t* lh_ctp, uint32_t* hh_ctp, uint32_t sample_f_ct, uint32_t* ll_ctfp, uint32_t* lh_ctfp, uint32_t* hh_ctfp, uint32_t hwe_or_geno_needed, uintptr_t sample_f_ctrl_ct, uint32_t* ll_hwep, uint32_t* lh_hwep, uint32_t* hh_hwep, int32_t hardy_needed, uintptr_t sample_f_case_ct, uint32_t* ll_case_hwep, uint32_t* lh_case_hwep, uint32_t* hh_case_hwep) {
  // This is best understood from the bottom third up (which is the order it
  // was written).  It's way overkill for just determining genotype
  // frequencies, but a ruthlessly optimized version is needed for e.g.
  // permutation testing so we may as well get it working here.
  //
  // The idea, which underlies the IBS and LD pruners as well, is to obtain
  // multiple popcounts at once.  In the case of PLINK frequency/HWE
  // evaluation, we need 6-9 numbers:
  // homozyg minor, het, homozyg major counts for all samples
  // homozyg minor, het, homozyg major counts for all founders
  // sometimes, homozyg minor, het, homozyg major counts for all ctrl founders
  //
  // Given a buffer with PLINK binary genotypes for a single marker, let
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  // Then,
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  //               = sample_ct - homozyg minor ct - het ct
  //
  // Thus, with the appropriate sample_ct and these three popcounts, we can
  // calculate a set of genotype counts.  We perform the
  // not-that-exploitative version of these calculations in the bottom third of
  // this function, to deal with the remainder that doesn't fit into the
  // 12-word block size of the main loops.
  //
  // The middle third is a 12-word block popcount for 32-bit platforms (see
  // popcount_longs() in plink_common.c; this is 12 words instead of 6 since
  // odd bits of the popcount targets are guaranteed to be zero, delaying
  // overflow).  It could be improved a bit, but we care more about reliability
  // than blistering speed for the 32-bit build (since it's used to check the
  // results of the actually blistering 64-bit code...).  Sure, hardware
  // popcount is significantly faster, but most machines running 32-bit OSes
  // don't have it.
  //
  // The top third is the portable Lauradoux/Walisch loop.
  uint32_t tot_a = 0;
  uint32_t tot_b = 0;
  uint32_t tot_c = 0;
  uint32_t tot_a_f = 0;
  uint32_t tot_b_f = 0;
  uint32_t tot_c_f = 0;
  uint32_t tot_a_hwe = 0;
  uint32_t tot_b_hwe = 0;
  uint32_t tot_c_hwe = 0;
  uint32_t tot_a_chwe = 0;
  uint32_t tot_b_chwe = 0;
  uint32_t tot_c_chwe = 0;
  uintptr_t* lptr_end = &(lptr[unfiltered_sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
#ifdef __LP64__
  uintptr_t cur_decr = 120;
  uintptr_t* lptr_12x_end;
  unfiltered_sample_ctl2 -= unfiltered_sample_ctl2 % 12;
  while (unfiltered_sample_ctl2 >= 120) {
  single_marker_freqs_and_hwe_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)sample_include2, &tot_a, &tot_b, &tot_c);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_include2, &tot_a_f, &tot_b_f, &tot_c_f);
    if (hwe_or_geno_needed) {
      count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_ctrl_include2, &tot_a_hwe, &tot_b_hwe, &tot_c_hwe);
      founder_ctrl_include2 = &(founder_ctrl_include2[cur_decr]);
      if (hardy_needed) {
	count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_case_include2, &tot_a_chwe, &tot_b_chwe, &tot_c_chwe);
	founder_case_include2 = &(founder_case_include2[cur_decr]);
      }
    }
    lptr = lptr_12x_end;
    sample_include2 = &(sample_include2[cur_decr]);
    founder_include2 = &(founder_include2[cur_decr]);
    unfiltered_sample_ctl2 -= cur_decr;
  }
  if (unfiltered_sample_ctl2) {
    cur_decr = unfiltered_sample_ctl2;
    goto single_marker_freqs_and_hwe_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[unfiltered_sample_ctl2 - unfiltered_sample_ctl2 % 12]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, sample_include2, &tot_a, &tot_b, &tot_c);
    count_3freq_12(lptr, founder_include2, &tot_a_f, &tot_b_f, &tot_c_f);
    if (hwe_or_geno_needed) {
      count_3freq_12(lptr, founder_ctrl_include2, &tot_a_hwe, &tot_b_hwe, &tot_c_hwe);
      founder_ctrl_include2 = &(founder_ctrl_include2[12]);
      if (hardy_needed) {
	count_3freq_12(lptr, founder_case_include2, &tot_a_chwe, &tot_b_chwe, &tot_c_chwe);
	founder_case_include2 = &(founder_case_include2[12]);
      }
    }
    lptr = &(lptr[12]);
    sample_include2 = &(sample_include2[12]);
    founder_include2 = &(founder_include2[12]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = *sample_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    // N.B. because of the construction of sample_include2, only even-numbered
    // bits can be present here.  So popcount2_long is safe.
    tot_a += popcount2_long(loader2);
    tot_b += popcount2_long(loader3);
    tot_c += popcount2_long(loader & loader3);
    loader2 = *founder_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_a_f += popcount2_long(loader2);
    tot_b_f += popcount2_long(loader3);
    tot_c_f += popcount2_long(loader & loader3);
    if (hwe_or_geno_needed) {
      loader2 = *founder_ctrl_include2++;
      loader3 = (loader >> 1) & loader2;
      loader2 &= loader;
      tot_a_hwe += popcount2_long(loader2);
      tot_b_hwe += popcount2_long(loader3);
      tot_c_hwe += popcount2_long(loader & loader3);
      if (hardy_needed) {
	loader2 = *founder_case_include2++;
	loader3 = (loader >> 1) & loader2;
	loader2 &= loader;
	tot_a_chwe += popcount2_long(loader2);
	tot_b_chwe += popcount2_long(loader3);
	tot_c_chwe += popcount2_long(loader & loader3);
      }
    }
  }
  *hh_ctp = tot_c;
  *lh_ctp = tot_b - tot_c;
  *ll_ctp = sample_ct - tot_a - *lh_ctp;
  *hh_ctfp = tot_c_f;
  *lh_ctfp = tot_b_f - tot_c_f;
  *ll_ctfp = sample_f_ct - tot_a_f - *lh_ctfp;
  if (hwe_or_geno_needed) {
    *hh_hwep = tot_c_hwe;
    *lh_hwep = tot_b_hwe - tot_c_hwe;
    *ll_hwep = sample_f_ctrl_ct - tot_a_hwe - *lh_hwep;
    if (hardy_needed) {
      *hh_case_hwep = tot_c_chwe;
      *lh_case_hwep = tot_b_chwe - tot_c_chwe;
      *ll_case_hwep = sample_f_case_ct - tot_a_chwe - *lh_case_hwep;
    }
  }
}

static inline uint32_t nonmissing_present_diff(uintptr_t unfiltered_sample_ctl2, uintptr_t* lptr, uintptr_t* sample_include2, uintptr_t* sample_male_include2) {
  // possible todo: write entries to .ynm file, using same format as .hh
  uintptr_t* lptr_end = &(lptr[unfiltered_sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  do {
    loader = *lptr++;
    loader2 = (*sample_include2++) & (~(*sample_male_include2++));
    // when really bored, check if compiler translates this into andnot
    // operations
    if ((~((~(loader >> 1)) & loader)) & loader2) {
      return 1;
    }
  } while (lptr < lptr_end);
  return 0;
}

static inline void haploid_single_marker_freqs(uintptr_t unfiltered_sample_ct, uintptr_t unfiltered_sample_ctl2, uintptr_t* lptr, uintptr_t* sample_include2, uintptr_t* founder_include2, uintptr_t sample_ct, uint32_t* ll_ctp, uint32_t* hh_ctp, uint32_t sample_f_ct, uint32_t* ll_ctfp, uint32_t* hh_ctfp, uint32_t* hethap_incr_ptr) {
  // Here, we interpret heterozygotes as missing.
  // Nonmissing: (genotype ^ (~(genotype >> 1))) & 0x5555...
  // Homozygote major: (genotype & (genotype >> 1)) & 0x5555...
  uint32_t tot_a = 0;
  uint32_t tot_b = 0;
  uint32_t tot_hmaj = 0;
  uint32_t tot_nm_f = 0;
  uint32_t tot_hmaj_f = 0;
  uintptr_t* lptr_end = &(lptr[unfiltered_sample_ctl2]);
  uint32_t hethap_incr;
  uint32_t tot_nm;
  uint32_t uii;
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
#ifdef __LP64__
  uintptr_t cur_decr = 120;
  uintptr_t* lptr_12x_end;
  unfiltered_sample_ctl2 -= unfiltered_sample_ctl2 % 12;
  while (unfiltered_sample_ctl2 >= 120) {
  haploid_single_marker_freqs_loop:
    lptr_12x_end = &(lptr[cur_decr]);
  // Given a buffer with PLINK binary genotypes for a single marker, let
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  // Then,
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  //               = sample_ct - homozyg minor ct - het ct
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)sample_include2, &tot_a, &tot_b, &tot_hmaj);
    freq_hwe_haploid_count_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_include2, &tot_nm_f, &tot_hmaj_f);
    lptr = lptr_12x_end;
    sample_include2 = &(sample_include2[cur_decr]);
    founder_include2 = &(founder_include2[cur_decr]);
    unfiltered_sample_ctl2 -= cur_decr;
  }
  if (unfiltered_sample_ctl2) {
    cur_decr = unfiltered_sample_ctl2;
    goto haploid_single_marker_freqs_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[unfiltered_sample_ctl2 - unfiltered_sample_ctl2 % 12]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, sample_include2, &tot_a, &tot_b, &tot_hmaj);
    freq_hwe_haploid_count_12(lptr, founder_include2, &tot_nm_f, &tot_hmaj_f);
    lptr = &(lptr[12]);
    sample_include2 = &(sample_include2[12]);
    founder_include2 = &(founder_include2[12]);
  }
#endif
  tot_nm = 2 * tot_hmaj + sample_ct - tot_a - tot_b;
  hethap_incr = tot_b - tot_hmaj;
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader3 = loader >> 1;
    loader4 = *sample_include2++;
    // different from tot_nm_f because of +sample_ct above
    tot_nm -= popcount2_long(loader & loader4);
    hethap_incr += popcount2_long(loader3 & (~loader) & loader4);
    loader2 = loader ^ (~loader3); // nonmissing?
    loader &= loader3; // homozyg A2?
    uii = popcount2_long(loader & loader4);
    tot_nm += 2 * uii - popcount2_long(loader3 & loader4);
    // tot_nm += popcount2_long(loader2 & loader4);
    tot_hmaj += uii;
    loader3 = *founder_include2++;
    tot_nm_f += popcount2_long(loader2 & loader3);
    tot_hmaj_f += popcount2_long(loader & loader3);
  }
  *hh_ctp = tot_hmaj;
  *ll_ctp = tot_nm - tot_hmaj;
  *hh_ctfp = tot_hmaj_f;
  *ll_ctfp = tot_nm_f - tot_hmaj_f;
  *hethap_incr_ptr = hethap_incr;
}

int32_t calc_freqs_and_hwe(FILE* bedfile, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_exclude_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* founder_info, int32_t nonfounders, int32_t maf_succ, double* set_allele_freqs, uintptr_t bed_offset, uint32_t hwe_needed, uint32_t hwe_all, uint32_t hardy_needed, uint32_t min_ac, uint32_t max_ac, double geno_thresh, uintptr_t* pheno_nm, uintptr_t* pheno_c, int32_t** hwe_lls_ptr, int32_t** hwe_lhs_ptr, int32_t** hwe_hhs_ptr, int32_t** hwe_ll_cases_ptr, int32_t** hwe_lh_cases_ptr, int32_t** hwe_hh_cases_ptr, int32_t** hwe_ll_allfs_ptr, int32_t** hwe_lh_allfs_ptr, int32_t** hwe_hh_allfs_ptr, int32_t** hwe_hapl_allfs_ptr, int32_t** hwe_haph_allfs_ptr, uintptr_t** geno_excl_bitfield_ptr, uintptr_t** ac_excl_bitfield_ptr, uint32_t* sample_male_ct_ptr, uint32_t* sample_f_ct_ptr, uint32_t* sample_f_male_ct_ptr, uintptr_t* topsize_ptr, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t is_split_chrom, uint32_t* hh_exists_ptr) {
  FILE* hhfile = NULL;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + BITCT - 1) / BITCT;
  uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  int32_t retval = 0;
  uint32_t pct = 1;
  uint32_t sample_ct = unfiltered_sample_ct - sample_exclude_ct;
  double sample_ct_recip = 1.0 / ((double)((int32_t)sample_ct));
  // sum of nonmissing rates over all markers
  // rate is in [0, 1] for each marker, so sum is in [0, marker_ct].
  double nonmissing_rate_tot = 0.0;
  // track this to defend against Y chromosome/0 males pathological case
  uintptr_t nonmissing_rate_tot_max = marker_ct;
  uint32_t sample_f_ct = sample_ct;
  uintptr_t sample_f_ctrl_ct = sample_ct;
  uintptr_t sample_f_case_ct = sample_ct;
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t ll_hwe = 0;
  uint32_t lh_hwe = 0;
  uint32_t hh_hwe = 0;
  uint32_t ll_case_hwe = 0;
  uint32_t lh_case_hwe = 0;
  uint32_t hh_case_hwe = 0;
  uint32_t cur_chrom_idx = 0;
  uint32_t nonmissing_nonmale_y = 0;
  int32_t ii = chrom_info_ptr->chrom_file_order[0];
  uint32_t is_haploid = is_set(chrom_info_ptr->haploid_mask, ii);
  uint32_t next_chrom_start = chrom_info_ptr->chrom_file_order_marker_idx[1];
  uint32_t is_x = (ii == chrom_info_ptr->x_code);
  uint32_t is_y = (ii == chrom_info_ptr->y_code);
  uint32_t ll_ct = 0;
  uint32_t lh_ct = 0;
  uint32_t hh_ct = 0;
  uint32_t ll_ctf = 0;
  uint32_t lh_ctf = 0;
  uint32_t hh_ctf = 0;
  uint32_t ukk = 0;
  uint32_t cur_oblig_missing = 0;
  uint32_t om_cluster_ct = 0;
  uint32_t* om_cluster_sizes = NULL;
  int32_t* hwe_lls = NULL;
  int32_t* hwe_lhs = NULL;
  int32_t* hwe_hhs = NULL;
  int32_t* hwe_ll_cases = NULL;
  int32_t* hwe_lh_cases = NULL;
  int32_t* hwe_hh_cases = NULL;
  int32_t* hwe_ll_allfs = NULL;
  int32_t* hwe_lh_allfs = NULL;
  int32_t* hwe_hh_allfs = NULL;
  uintptr_t* sample_nonmale_include2 = NULL;
  uintptr_t* sample_male_include2 = NULL;
  uintptr_t* founder_nonmale_include2 = NULL;
  uintptr_t* founder_ctrl_nonmale_include2 = NULL;
  uintptr_t* founder_male_include2 = NULL;
  uintptr_t* founder_case_include2 = NULL;
  uintptr_t* founder_case_nonmale_include2 = NULL;
  uintptr_t* geno_excl_bitfield = NULL;
  uintptr_t* ac_excl_bitfield = NULL;
  uint64_t* om_entry_ptr = NULL;
  uint32_t sample_nonmale_ct = 0;
  uint32_t sample_f_nonmale_ct = 0;
  uint32_t sample_f_ctl_nonmale_ct = 0;
  uint32_t sample_f_case_nonmale_ct = 0;
  uint64_t hethap_ct = 0;
  uint64_t cur_om_entry = 0;
  double male_ct_recip = 0;
  uint32_t* om_sample_lookup;
  int32_t* hwe_hapl_allfs;
  int32_t* hwe_haph_allfs;
  uintptr_t* loadbuf;
  uintptr_t* sample_include2;
  uintptr_t* founder_include2;
  uintptr_t* founder_ctrl_include2;
  uintptr_t* tmp_sample_excl_mask;
  uintptr_t* tmp_sample_excl_mask2;
  uintptr_t loop_end;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t ulii;
  uint32_t sample_male_ct;
  uint32_t sample_f_male_ct;
  uint32_t hethap_incr;
  uint32_t nonmales_needed;
  uint32_t males_needed;
  uint32_t uii;
  uint32_t ujj;
  double maf;
  double cur_genotyping_rate;
  if (!hwe_needed) {
    *hwe_lls_ptr = (int32_t*)wkspace_base;
  } else {
    if (wkspace_alloc_i_checked(&hwe_lls, unfiltered_marker_ct * sizeof(int32_t)) ||
	wkspace_alloc_i_checked(&hwe_lhs, unfiltered_marker_ct * sizeof(int32_t)) ||
	wkspace_alloc_i_checked(&hwe_hhs, unfiltered_marker_ct * sizeof(int32_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    *hwe_lls_ptr = hwe_lls;
    *hwe_lhs_ptr = hwe_lhs;
    *hwe_hhs_ptr = hwe_hhs;
    if (hardy_needed) {
      if (wkspace_alloc_i_checked(&hwe_ll_cases, unfiltered_marker_ct * sizeof(int32_t)) ||
          wkspace_alloc_i_checked(&hwe_lh_cases, unfiltered_marker_ct * sizeof(int32_t)) ||
          wkspace_alloc_i_checked(&hwe_hh_cases, unfiltered_marker_ct * sizeof(int32_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
    }
    *hwe_ll_cases_ptr = hwe_ll_cases;
    *hwe_lh_cases_ptr = hwe_lh_cases;
    *hwe_hh_cases_ptr = hwe_hh_cases;
  }
  if (wkspace_alloc_i_checked(&hwe_ll_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_lh_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_hh_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_hapl_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_haph_allfs, unfiltered_marker_ct * sizeof(int32_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  *hwe_ll_allfs_ptr = hwe_ll_allfs;
  *hwe_lh_allfs_ptr = hwe_lh_allfs;
  *hwe_hh_allfs_ptr = hwe_hh_allfs;
  *hwe_hapl_allfs_ptr = hwe_hapl_allfs;
  *hwe_haph_allfs_ptr = hwe_haph_allfs;

  if ((!pheno_c) || is_split_chrom) {
    hwe_all = 1;
  }

  fill_int_zero(hwe_ll_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_lh_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_hh_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_hapl_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_haph_allfs, unfiltered_marker_ct);
  if (geno_thresh < 1.0) {
    if (wkspace_alloc_ul_checked(geno_excl_bitfield_ptr, unfiltered_marker_ctl * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    geno_excl_bitfield = *geno_excl_bitfield_ptr;
    fill_ulong_zero(geno_excl_bitfield, unfiltered_marker_ctl);
    // change this to a minimum nonmissing rate
    geno_thresh = (1.0 - geno_thresh) * (1 - SMALL_EPSILON);
  }
  if ((min_ac > 0) || (max_ac < sample_ct)) {
    if (wkspace_alloc_ul_checked(ac_excl_bitfield_ptr, unfiltered_marker_ctl * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    ac_excl_bitfield = *ac_excl_bitfield_ptr;
    fill_ulong_zero(ac_excl_bitfield, unfiltered_marker_ctl);
  }
  wkspace_mark = wkspace_base;
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_sample_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&sample_include2, unfiltered_sample_ctv2 * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  loadbuf[unfiltered_sample_ctv2 - 2] = 0;
  loadbuf[unfiltered_sample_ctv2 - 1] = 0;
  exclude_to_vec_include(unfiltered_sample_ct, sample_include2, sample_exclude);
  ii = chrom_info_ptr->x_code;
  nonmales_needed = (!is_split_chrom) && (ii != -1) && is_set(chrom_info_ptr->chrom_mask, ii);
  ii = chrom_info_ptr->y_code;
  males_needed = nonmales_needed || ((!is_split_chrom) && (ii != -1) && is_set(chrom_info_ptr->chrom_mask, ii));
  if (wkspace_alloc_ul_checked(&sample_male_include2, unfiltered_sample_ctv2 * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  memcpy(sample_male_include2, sample_include2, unfiltered_sample_ctv2 * sizeof(intptr_t));
  vec_include_mask_in(unfiltered_sample_ct, sample_male_include2, sex_male);
  sample_male_ct = popcount01_longs(sample_male_include2, unfiltered_sample_ctv2);
  if (sample_male_ct) {
    male_ct_recip = 1.0 / ((double)((int32_t)sample_male_ct));
  }
  *sample_male_ct_ptr = sample_male_ct;
  sample_f_male_ct = sample_male_ct;
  if (males_needed) {
    founder_male_include2 = sample_male_include2;
    if (nonmales_needed) {
      if (wkspace_alloc_ul_checked(&sample_nonmale_include2, unfiltered_sample_ctv2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(sample_nonmale_include2, sample_include2, unfiltered_sample_ctv2 * sizeof(intptr_t));
      vec_include_mask_out_intersect(unfiltered_sample_ct, sample_nonmale_include2, sex_nm, sex_male);
      sample_nonmale_ct = popcount01_longs(sample_nonmale_include2, unfiltered_sample_ctv2);
      sample_f_nonmale_ct = sample_nonmale_ct;
      founder_nonmale_include2 = sample_nonmale_include2;
      founder_ctrl_nonmale_include2 = sample_nonmale_include2;
    }
  }
  founder_include2 = sample_include2;
  if (wkspace_alloc_ul_checked(&tmp_sample_excl_mask, unfiltered_sample_ctl * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  memcpy(tmp_sample_excl_mask, sample_exclude, unfiltered_sample_ctl * sizeof(intptr_t));
  if (!nonfounders) {
    if (wkspace_alloc_ul_checked(&founder_include2, unfiltered_sample_ctv2 * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    bitfield_ornot(tmp_sample_excl_mask, founder_info, unfiltered_sample_ctl);
    zero_trailing_bits(tmp_sample_excl_mask, unfiltered_sample_ct);
    exclude_to_vec_include(unfiltered_sample_ct, founder_include2, tmp_sample_excl_mask);
    if (males_needed) {
      if (wkspace_alloc_ul_checked(&founder_male_include2, unfiltered_sample_ctv2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(founder_male_include2, sample_male_include2, unfiltered_sample_ctl * 2 * sizeof(intptr_t));
      vec_include_mask_in(unfiltered_sample_ct, founder_male_include2, founder_info);
      sample_f_male_ct = popcount01_longs(founder_male_include2, unfiltered_sample_ctv2);
      if (nonmales_needed) {
	if (wkspace_alloc_ul_checked(&founder_nonmale_include2, unfiltered_sample_ctv2 * sizeof(intptr_t))) {
	  goto calc_freqs_and_hwe_ret_NOMEM;
	}
	memcpy(founder_nonmale_include2, sample_nonmale_include2, unfiltered_sample_ctv2 * sizeof(intptr_t));
	vec_include_mask_in(unfiltered_sample_ct, founder_nonmale_include2, founder_info);
	sample_f_nonmale_ct = popcount01_longs(founder_nonmale_include2, unfiltered_sample_ctv2);
      }
    }
    founder_ctrl_include2 = founder_include2;
    sample_f_ct = popcount_longs_exclude(founder_info, sample_exclude, unfiltered_sample_ctl);
    sample_f_ctrl_ct = sample_f_ct;
  } else {
    founder_ctrl_include2 = founder_include2;
  }

  // bugfix: this previously failed to initialize founder_ctrl_include2 and
  // founder_case_include2 properly if --hardy was used in a situation where
  // hwe_all would be set (e.g. all-case datasets).
  if ((!hwe_all) || hardy_needed) {
    if (wkspace_alloc_ul_checked(&founder_ctrl_include2, unfiltered_sample_ctv2 *  sizeof(intptr_t)) ||
	wkspace_alloc_ul_checked(&tmp_sample_excl_mask2, unfiltered_sample_ctl * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    memcpy(tmp_sample_excl_mask2, tmp_sample_excl_mask, unfiltered_sample_ctl * sizeof(intptr_t));
    bitfield_ornot(tmp_sample_excl_mask2, pheno_nm, unfiltered_sample_ctl);
    bitfield_or(tmp_sample_excl_mask2, pheno_c, unfiltered_sample_ctl);
    zero_trailing_bits(tmp_sample_excl_mask2, unfiltered_sample_ct);
    // tmp_sample_excl_mask2 is now set for each sample who is excluded, or a
    // nonfounder, or is noncontrol.
    sample_f_ctrl_ct = unfiltered_sample_ct - popcount_longs(tmp_sample_excl_mask2, unfiltered_sample_ctl);
    exclude_to_vec_include(unfiltered_sample_ct, founder_ctrl_include2, tmp_sample_excl_mask2);
    if (nonmales_needed) {
      if (wkspace_alloc_ul_checked(&founder_ctrl_nonmale_include2, unfiltered_sample_ctv2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(founder_ctrl_nonmale_include2, sample_nonmale_include2, unfiltered_sample_ctv2 * sizeof(intptr_t));
      vec_include_mask_out(unfiltered_sample_ct, founder_ctrl_nonmale_include2, tmp_sample_excl_mask2);
      sample_f_ctl_nonmale_ct = popcount01_longs(founder_ctrl_nonmale_include2, unfiltered_sample_ctv2);
    }
    if (hardy_needed) {
      if (wkspace_alloc_ul_checked(&founder_case_include2, unfiltered_sample_ctv2 *  sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      bitfield_ornot(tmp_sample_excl_mask, pheno_nm, unfiltered_sample_ctl);
      bitfield_ornot(tmp_sample_excl_mask, pheno_c, unfiltered_sample_ctl);
      zero_trailing_bits(tmp_sample_excl_mask, unfiltered_sample_ct);
      sample_f_case_ct = unfiltered_sample_ct - popcount_longs(tmp_sample_excl_mask, unfiltered_sample_ctl);
      exclude_to_vec_include(unfiltered_sample_ct, founder_case_include2, tmp_sample_excl_mask);
      if (nonmales_needed) {
	if (wkspace_alloc_ul_checked(&founder_case_nonmale_include2, unfiltered_sample_ctv2 * sizeof(intptr_t))) {
	  goto calc_freqs_and_hwe_ret_NOMEM;
	}
	memcpy(founder_case_nonmale_include2, sample_nonmale_include2, unfiltered_sample_ctv2 * sizeof(intptr_t));
	vec_include_mask_out(unfiltered_sample_ct, founder_case_nonmale_include2, tmp_sample_excl_mask);
	sample_f_case_nonmale_ct = popcount01_longs(founder_case_nonmale_include2, unfiltered_sample_ctv2);
      }
    }
  }
  if (om_ip->entry_ct) {
    om_entry_ptr = om_ip->entries;
    om_cluster_ct = om_ip->cluster_ct;
    om_sample_lookup = om_ip->sample_lookup;
    cur_om_entry = *om_entry_ptr;
    if (wkspace_alloc_ui_checked(&om_cluster_sizes, om_cluster_ct * 2 * sizeof(int32_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    fill_uint_zero(om_cluster_sizes, om_cluster_ct * 2);
    sample_uidx = 0;
    for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      uii = om_sample_lookup[sample_uidx];
      if (uii != 0xffffffffU) {
        om_cluster_sizes[uii] += 1;
        if (is_set(sex_male, sample_uidx)) {
          om_cluster_sizes[uii + om_cluster_ct] += 1;
	}
      }
    }
  }

  *sample_f_ct_ptr = sample_f_ct;
  *sample_f_male_ct_ptr = sample_f_male_ct;
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto calc_freqs_and_hwe_ret_READ_FAIL;
  }
  marker_uidx = 0;
  marker_idx = 0;
  logprint("Calculating allele frequencies...");
  fputs(" 0%", stdout);
  fflush(stdout);
  if (is_split_chrom) {
    // only set is_haploid if all chromosomes are haploid
    is_haploid = (chrom_info_ptr->haploid_mask[0]) & 1;
    is_x = 0;
    is_y = 0;
    next_chrom_start = unfiltered_marker_ct;
  }
  for (; pct <= 100; pct++) {
    loop_end = ((uint64_t)pct * marker_ct) / 100LU;
    for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto calc_freqs_and_hwe_ret_READ_FAIL;
	}
      }
      if (load_raw(bedfile, loadbuf, unfiltered_sample_ct4)) {
	goto calc_freqs_and_hwe_ret_READ_FAIL;
      }
      if (marker_uidx >= next_chrom_start) {
	do {
	  next_chrom_start = chrom_info_ptr->chrom_file_order_marker_idx[(++cur_chrom_idx) + 1];
	} while (marker_uidx >= next_chrom_start);
	ii = chrom_info_ptr->chrom_file_order[cur_chrom_idx];
	is_haploid = is_set(chrom_info_ptr->haploid_mask, ii);
	is_x = (ii == chrom_info_ptr->x_code);
	is_y = (ii == chrom_info_ptr->y_code);
      }
      if (om_entry_ptr) {
        cur_oblig_missing = 0;
	while ((cur_om_entry >> 32) < marker_uidx) {
	  cur_om_entry = *(++om_entry_ptr);
	}
        while ((cur_om_entry >> 32) == marker_uidx) {
          cur_oblig_missing += om_cluster_sizes[(uint32_t)cur_om_entry];
          cur_om_entry = *(++om_entry_ptr);
	}
      }
      if (!is_haploid) {
	single_marker_freqs_and_hwe(unfiltered_sample_ctv2, loadbuf, sample_include2, founder_include2, founder_ctrl_include2, founder_case_include2, sample_ct, &ll_ct, &lh_ct, &hh_ct, sample_f_ct, &ll_ctf, &lh_ctf, &hh_ctf, hwe_needed, sample_f_ctrl_ct, &ll_hwe, &lh_hwe, &hh_hwe, hardy_needed, sample_f_case_ct, &ll_case_hwe, &lh_case_hwe, &hh_case_hwe);
	hwe_ll_allfs[marker_uidx] = ll_ctf;
	hwe_lh_allfs[marker_uidx] = lh_ctf;
	hwe_hh_allfs[marker_uidx] = hh_ctf;
	uii = ll_ct + lh_ct + hh_ct;
	if (!cur_oblig_missing) {
	  cur_genotyping_rate = ((int32_t)uii) * sample_ct_recip;
	} else {
	  if (sample_ct - cur_oblig_missing) {
	    cur_genotyping_rate = ((int32_t)uii) / ((double)((int32_t)(sample_ct - cur_oblig_missing)));
	  } else {
	    cur_genotyping_rate = 0;
	    nonmissing_rate_tot_max -= 1;
	  }
	}
	if (ac_excl_bitfield) {
	  if (ll_ctf < hh_ctf) {
	    uii = 2 * ll_ctf + lh_ctf;
	  } else {
	    uii = 2 * hh_ctf + lh_ctf;
	  }
	  if ((uii < min_ac) || (uii > max_ac)) {
	    set_bit(ac_excl_bitfield, marker_uidx);
	  }
	}
	uii = 2 * (ll_ctf + lh_ctf + hh_ctf + maf_succ);
	if (!uii) {
	  // avoid 0/0 division
	  set_allele_freqs[marker_uidx] = 0.5;
	} else {
	  set_allele_freqs[marker_uidx] = ((double)(2 * hh_ctf + lh_ctf + maf_succ)) / ((double)uii);
	}
	if (hwe_needed) {
	  hwe_lls[marker_uidx] = ll_hwe;
	  hwe_lhs[marker_uidx] = lh_hwe;
	  hwe_hhs[marker_uidx] = hh_hwe;
	  if (hardy_needed) {
	    hwe_ll_cases[marker_uidx] = ll_case_hwe;
	    hwe_lh_cases[marker_uidx] = lh_case_hwe;
	    hwe_hh_cases[marker_uidx] = hh_case_hwe;
	  }
	}
      } else {
	uii = 0;
	ujj = 0;
	if (is_x || is_y) {
	  if (is_x) {
	    single_marker_freqs_and_hwe(unfiltered_sample_ctv2, loadbuf, sample_nonmale_include2, founder_nonmale_include2, founder_ctrl_nonmale_include2, founder_case_nonmale_include2, sample_nonmale_ct, &ll_ct, &lh_ct, &hh_ct, sample_f_nonmale_ct, &ll_ctf, &lh_ctf, &hh_ctf, hwe_needed, sample_f_ctl_nonmale_ct, &ll_hwe, &lh_hwe, &hh_hwe, hardy_needed, sample_f_case_nonmale_ct, &ll_case_hwe, &lh_case_hwe, &hh_case_hwe);
	    hwe_ll_allfs[marker_uidx] = ll_ctf;
	    hwe_lh_allfs[marker_uidx] = lh_ctf;
	    hwe_hh_allfs[marker_uidx] = hh_ctf;
	    uii = 2 * (ll_ctf + lh_ctf + hh_ctf);
	    ujj = 2 * hh_ctf + lh_ctf;
	    ukk = ll_ct + lh_ct + hh_ct;
	    if (hwe_needed) {
	      hwe_lls[marker_uidx] = ll_hwe;
	      hwe_lhs[marker_uidx] = lh_hwe;
	      hwe_hhs[marker_uidx] = hh_hwe;
	      if (hardy_needed) {
		hwe_ll_cases[marker_uidx] = ll_case_hwe;
		hwe_lh_cases[marker_uidx] = lh_case_hwe;
		hwe_hh_cases[marker_uidx] = hh_case_hwe;
	      }
	    }
	  } else if (!nonmissing_nonmale_y) {
	    nonmissing_nonmale_y = nonmissing_present_diff(unfiltered_sample_ctv2, loadbuf, sample_include2, sample_male_include2);
	  }
	  haploid_single_marker_freqs(unfiltered_sample_ct, unfiltered_sample_ctv2, loadbuf, sample_male_include2, founder_male_include2, sample_male_ct, &ll_ct, &hh_ct, sample_f_male_ct, &ll_ctf, &hh_ctf, &hethap_incr);
	  if ((is_x || sample_male_ct) && (sample_ct - cur_oblig_missing)) {
	    if (is_x) {
	      if (!cur_oblig_missing) {
		cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct + ukk)) * sample_ct_recip;
	      } else {
		cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct + ukk)) / ((double)((int32_t)(sample_ct - cur_oblig_missing)));
	      }
	    } else {
	      if (!cur_oblig_missing) {
		cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct)) * male_ct_recip;
	      } else {
		cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct)) / ((double)((int32_t)(sample_male_ct - cur_oblig_missing)));
	      }
	    }
	  } else {
	    cur_genotyping_rate = 0;
	    nonmissing_rate_tot_max -= 1;
	  }
	} else {
	  haploid_single_marker_freqs(unfiltered_sample_ct, unfiltered_sample_ctv2, loadbuf, sample_include2, founder_include2, sample_ct, &ll_ct, &hh_ct, sample_f_ct, &ll_ctf, &hh_ctf, &hethap_incr);
	  if (!cur_oblig_missing) {
	    cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct)) * sample_ct_recip;
	  } else {
	    if (sample_ct - cur_oblig_missing) {
	      cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct)) / ((double)((int32_t)(sample_ct - cur_oblig_missing)));
	    } else {
	      cur_genotyping_rate = 0;
	      nonmissing_rate_tot_max -= 1;
	    }
	  }
	}
	if (hethap_incr) {
	  if (!hhfile) {
	    memcpy(outname_end, ".hh", 4);
	    if (fopen_checked(&hhfile, outname, "w")) {
	      goto calc_freqs_and_hwe_ret_OPEN_FAIL;
	    }
	  }
	  if (is_x) {
	    *hh_exists_ptr |= XMHH_EXISTS;
	  } else if (is_y) {
	    *hh_exists_ptr |= Y_FIX_NEEDED;
	  } else {
	    *hh_exists_ptr |= NXMHH_EXISTS;
	  }
	  if (is_x || is_y) {
	    for (sample_uidx = 0; sample_uidx < unfiltered_sample_ctv2; sample_uidx++) {
	      ulii = loadbuf[sample_uidx];
	      ulii = (ulii >> 1) & (~ulii) & sample_male_include2[sample_uidx];
	      while (ulii) {
		ukk = sample_uidx * BITCT2 + CTZLU(ulii) / 2;
		fputs(&(sample_ids[ukk * max_sample_id_len]), hhfile);
		putc('\t', hhfile);
		fputs(&(marker_ids[marker_uidx * max_marker_id_len]), hhfile);
		putc('\n', hhfile);
		ulii &= ulii - ONELU;
	      }
	    }
	  } else {
	    for (sample_uidx = 0; sample_uidx < unfiltered_sample_ctv2; sample_uidx++) {
	      ulii = loadbuf[sample_uidx];
	      ulii = (ulii >> 1) & (~ulii) & sample_include2[sample_uidx];
	      while (ulii) {
		ukk = sample_uidx * BITCT2 + CTZLU(ulii) / 2;
		fputs(&(sample_ids[ukk * max_sample_id_len]), hhfile);
		putc('\t', hhfile);
		fputs(&(marker_ids[marker_uidx * max_marker_id_len]), hhfile);
		putc('\n', hhfile);
		ulii &= ulii - ONELU;
	      }
	    }
	  }
	  if (ferror(hhfile)) {
	    goto calc_freqs_and_hwe_ret_WRITE_FAIL;
	  }
	  hethap_ct += hethap_incr;
	}
	hwe_hapl_allfs[marker_uidx] = ll_ctf;
	hwe_haph_allfs[marker_uidx] = hh_ctf;
	uii += ll_ctf + hh_ctf;
	ujj += hh_ctf;
	if (ac_excl_bitfield) {
	  if (ujj <= uii / 2) {
	    ukk = ujj;
	  } else {
	    ukk = uii - ujj;
	  }
	  if ((ukk < min_ac) || (ukk > max_ac)) {
	    set_bit(ac_excl_bitfield, marker_uidx);
	  }
	}
	uii += 2 * maf_succ;
	ujj += maf_succ;
	if (!uii) {
	  maf = 0.5;
	} else {
	  maf = ((double)ujj) / ((double)uii);
	}
	set_allele_freqs[marker_uidx] = maf;
      }
      nonmissing_rate_tot += cur_genotyping_rate;
      if (geno_excl_bitfield && (cur_genotyping_rate < geno_thresh)) {
	SET_BIT(geno_excl_bitfield, marker_uidx);
      }
    }
    if (pct < 100) {
      if (pct > 10) {
	putchar('\b');
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
  }
  fputs("\b\b\b\b", stdout);
  logprint(" done.\n");
  if (hethap_ct) {
    *outname_end = '\0';
    LOGPRINTFWW("Warning: %" PRIu64 " het. haploid genotype%s present (see %s.hh ).\n", hethap_ct, (hethap_ct == 1LLU)? "" : "s", outname);
  }
  if (nonmissing_nonmale_y) {
    logprint("Warning: Nonmissing nonmale Y chromosome genotype(s) present.\n");
    *hh_exists_ptr |= Y_FIX_NEEDED;
  }
  if (nonmissing_rate_tot <= 0.9999995 * ((double)((intptr_t)nonmissing_rate_tot_max))) {
    LOGPRINTF("Total genotyping rate %sis %g.\n", sample_exclude_ct? "in remaining samples " : "", nonmissing_rate_tot / ((double)((intptr_t)nonmissing_rate_tot_max)));
  }
  while (0) {
  calc_freqs_and_hwe_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_freqs_and_hwe_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_freqs_and_hwe_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  calc_freqs_and_hwe_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(hhfile);
  return retval;
}

int32_t write_missingness_reports(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t output_gz, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_sample_ct, uintptr_t sample_ct, uintptr_t* sample_exclude, uintptr_t* pheno_nm, uintptr_t* sex_male, uint32_t sample_male_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t hh_exists) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ct2l = (unfiltered_sample_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t unfiltered_sample_ctv2 = (unfiltered_sample_ct2l + 1) & (~1);
  uintptr_t marker_ct_y = 0;
  uintptr_t* sample_male_include2 = NULL;
  uint64_t* om_entry_ptr = NULL;
  uintptr_t* cur_omidxs = NULL;
  char* pzwritep = NULL;
  uint32_t* sample_to_cluster = NULL;
  uint32_t* missing_ct_by_cluster = NULL;
  uint32_t* oblig_missing_ct_by_cluster = NULL;
  uint32_t* cluster_sizes = NULL;
  uint32_t* cluster_sizes_y = NULL;
  uint32_t* om_cluster_sizes = NULL;
  uint32_t* om_sample_lookup = NULL;
  uint32_t* om_cluster_ref_cts = NULL;
  uint64_t cur_om_entry = 0;
  int32_t y_code = chrom_info_ptr->y_code;
  uint32_t y_present = (y_code != -1) && is_set(chrom_info_ptr->chrom_mask, y_code);
  uint32_t sample_uidx = 0;
  uint32_t sample_idx = 0;
  uint32_t oblig_ct = 0;
  uint32_t om_cluster_ct = 0;
  uint32_t om_cluster_ctl = 0;
  int32_t retval = 0;
  Pigz_state ps;
  uintptr_t* loadbuf;
  uintptr_t* sample_include2;
  uintptr_t* cur_nm;
  uintptr_t* lptr;
  uintptr_t* lptr2;
  uint32_t* missing_cts;
  uint32_t* cur_cluster_sizes;
  unsigned char* overflow_buf;
  char* cptr;
  char* cptr2;
  uintptr_t marker_ct_nony;
  uintptr_t marker_uidx;
  uintptr_t ulii;
  uint32_t slen;
  uint32_t sample_uidx_stop;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t cur_tot;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t om_ycorr;
  uint32_t clidx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  pzwrite_init_null(&ps);
  if (wkspace_alloc_uc_checked(&overflow_buf, PIGZ_BLOCK_SIZE + MAXLINELEN) ||
      wkspace_alloc_ui_checked(&missing_cts, unfiltered_sample_ct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, unfiltered_sample_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&sample_include2, unfiltered_sample_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&sample_male_include2, unfiltered_sample_ctv2 * sizeof(intptr_t))) {
    goto write_missingness_reports_ret_NOMEM;
  }
  loadbuf[unfiltered_sample_ctv2 - 2] = 0;
  loadbuf[unfiltered_sample_ctv2 - 1] = 0;
  exclude_to_vec_include(unfiltered_sample_ct, sample_include2, sample_exclude);
  memcpy(sample_male_include2, sample_include2, unfiltered_sample_ctv2 * sizeof(intptr_t));
  vec_include_mask_in(unfiltered_sample_ct, sample_male_include2, sex_male);
  if (y_present) {
    marker_ct_y = count_chrom_markers(chrom_info_ptr, chrom_info_ptr->y_code, marker_exclude);
  }
  marker_ct_nony = marker_ct - marker_ct_y;
  fill_uint_zero(missing_cts, unfiltered_sample_ct);
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto write_missingness_reports_ret_READ_FAIL;
  }
  memcpy(outname_end, output_gz? ".lmiss.gz" : ".lmiss", output_gz? 10 : 7);
  if (flex_pzwrite_init(output_gz, outname, overflow_buf, 0, &ps)) {
    goto write_missingness_reports_ret_OPEN_FAIL;
  }
  pzwritep = (char*)overflow_buf;

  if (om_ip->entry_ct) {
    om_entry_ptr = om_ip->entries;
    om_cluster_ref_cts = om_ip->cluster_ref_cts;
    cur_om_entry = *om_entry_ptr;
    om_cluster_ct = om_ip->cluster_ct;
    // divide by BITCT2 instead of BITCT due to Ychr
    om_cluster_ctl = (om_cluster_ct + BITCT - 1) / BITCT;
    om_sample_lookup = om_ip->sample_lookup;
    if (wkspace_alloc_ui_checked(&om_cluster_sizes, om_cluster_ct * 2 * sizeof(int32_t))) {
      goto write_missingness_reports_ret_NOMEM;
    }
    fill_uint_zero(om_cluster_sizes, om_cluster_ct * 2);
    for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_unsafe_ck(sample_exclude, &sample_uidx);
      uii = om_sample_lookup[sample_uidx];
      if (uii != 0xffffffffU) {
        om_cluster_sizes[uii] += 1;
        if (is_set(sex_male, sample_uidx)) {
	  om_cluster_sizes[uii + om_cluster_ct] += 1;
	}
      }
    }
    sample_uidx = 0;
    sample_idx = 0;
    if (cluster_ct) {
      if (wkspace_alloc_ul_checked(&cur_omidxs, om_cluster_ctl * sizeof(intptr_t))) {
        goto write_missingness_reports_ret_NOMEM;
      }
    }
  }
  ujj = unfiltered_sample_ct2l * BITCT2;
  if (!cluster_ct) {
    sprintf(tbuf, " CHR %%%us   N_MISS   N_GENO   F_MISS" EOLN_STR, plink_maxsnp);
  } else {
    if (wkspace_alloc_ui_checked(&sample_to_cluster, unfiltered_sample_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&missing_ct_by_cluster, cluster_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&oblig_missing_ct_by_cluster, cluster_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_sizes, cluster_ct * 2 * sizeof(int32_t))) {
      goto write_missingness_reports_ret_NOMEM;
    }
    fill_uint_zero(sample_to_cluster, unfiltered_sample_ct);
    fill_uint_zero(cluster_sizes, cluster_ct * 2);
    fill_uint_zero(oblig_missing_ct_by_cluster, cluster_ct);
    cluster_sizes_y = &(cluster_sizes[cluster_ct]);
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      unn = cluster_starts[clidx + 1];
      ukk = clidx + 1;
      for (uii = cluster_starts[clidx]; uii < unn; uii++) {
	umm = cluster_map[uii];
	if (!IS_SET(sample_exclude, umm)) {
          sample_to_cluster[umm] = ukk;
	  cluster_sizes[clidx] += 1;
          if (IS_SET(sex_male, umm)) {
            cluster_sizes_y[clidx] += 1;
	  }
	}
      }
    }
    sprintf(tbuf, " CHR %%%us       CLST   N_MISS   N_CLST   N_GENO   F_MISS" EOLN_STR, plink_maxsnp);
  }

  pzwritep += sprintf(pzwritep, tbuf, "SNP");
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    is_x = (((int32_t)chrom_idx) == chrom_info_ptr->x_code);
    is_y = (((int32_t)chrom_idx) == chrom_info_ptr->y_code);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    if (!is_y) {
      cur_nm = sample_include2;
      cur_tot = sample_ct;
      cur_cluster_sizes = cluster_sizes;
      om_ycorr = 0;
    } else {
      cur_nm = sample_male_include2;
      cur_tot = sample_male_ct;
      cur_cluster_sizes = cluster_sizes_y;
      om_ycorr = om_cluster_ct;
    }
    cptr = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, chrom_idx));
    *cptr++ = ' ';
    if (marker_uidx < chrom_end) {
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	goto write_missingness_reports_ret_READ_FAIL;
      }
      do {
	if (load_raw(bedfile, loadbuf, unfiltered_sample_ct4)) {
	  goto write_missingness_reports_ret_READ_FAIL;
	}
        if (is_haploid) {
          haploid_fix(hh_exists, sample_include2, sample_male_include2, unfiltered_sample_ct, is_x, is_y, (unsigned char*)loadbuf);
	}
	lptr = loadbuf;
	lptr2 = cur_nm;
	cptr2 = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), cptr);
	*cptr2++ = ' ';
	if (om_entry_ptr) {
	  while ((cur_om_entry >> 32) < marker_uidx) {
	    cur_om_entry = *(++om_entry_ptr);
	  }
	  if (cluster_ct) {
	    fill_uint_zero(oblig_missing_ct_by_cluster, cluster_ct);
	  }
	}
	if (!cluster_ct) {
	  if (om_entry_ptr) {
	    oblig_ct = 0;
	    while ((cur_om_entry >> 32) == marker_uidx) {
	      oblig_ct += om_cluster_sizes[(uint32_t)cur_om_entry];
	      cur_om_entry = *(++om_entry_ptr);
	    }
	  }
          ukk = 0;
	  for (uii = 0; uii < ujj; uii += BITCT2) {
	    ulii = *lptr++;
	    ulii = ulii & ((~ulii) >> 1) & (*lptr2++);
	    while (ulii) {
	      missing_cts[uii + (CTZLU(ulii) / 2)] += 1;
	      ukk++;
	      ulii &= ulii - 1;
	    }
	  }
          pzwritep = memcpya(pzwritep, tbuf, cptr2 - tbuf);
	  pzwritep = uint32_writew8x(pzwritep, ukk - oblig_ct, ' ');
          pzwritep = uint32_writew8x(pzwritep, cur_tot - oblig_ct, ' ');
	  pzwritep = double_g_writewx4(pzwritep, ((double)((int32_t)(ukk - oblig_ct))) / ((double)((int32_t)(cur_tot - oblig_ct))), 8);
          append_binary_eoln(&pzwritep);
	  if (flex_pzwrite(&ps, &pzwritep)) {
	    goto write_missingness_reports_ret_WRITE_FAIL;
	  }
	} else {
	  fill_uint_zero(missing_ct_by_cluster, cluster_ct);
	  if ((!om_entry_ptr) || ((cur_om_entry >> 32) != marker_uidx)) {
	    for (uii = 0; uii < ujj; uii += BITCT2) {
	      ulii = *lptr++;
	      ulii = ulii & ((~ulii) >> 1) & (*lptr2++);
	      while (ulii) {
		ukk = uii + (CTZLU(ulii) / 2);
		missing_cts[ukk] += 1;
		ukk = sample_to_cluster[ukk];
		if (ukk) {
		  missing_ct_by_cluster[ukk - 1] += 1;
		}
		ulii &= ulii - 1;
	      }
	    }
	  } else {
	    fill_ulong_zero(cur_omidxs, om_cluster_ctl);
	    do {
              set_bit(cur_omidxs, ((uint32_t)cur_om_entry) - om_ycorr);
	      cur_om_entry = *(++om_entry_ptr);
	    } while ((cur_om_entry >> 32) == marker_uidx);
	    for (uii = 0; uii < ujj; uii += BITCT2) {
	      ulii = *lptr++;
	      ulii = ulii & ((~ulii) >> 1) & (*lptr2++);
	      while (ulii) {
		ukk = uii + (CTZLU(ulii) / 2);
		missing_cts[ukk] += 1;
		umm = sample_to_cluster[ukk];
		if (umm) {
		  if (is_set(cur_omidxs, om_sample_lookup[ukk])) {
		    oblig_missing_ct_by_cluster[umm - 1] += 1;
		  } else {
		    missing_ct_by_cluster[umm - 1] += 1;
		  }
		}
		ulii &= ulii - 1;
	      }
	    }
	  }
	  for (clidx = 0; clidx < cluster_ct; clidx++) {
            pzwritep = memcpya(pzwritep, tbuf, cptr2 - tbuf);
            pzwritep = fw_strcpy(10, &(cluster_ids[clidx * max_cluster_id_len]), pzwritep);
	    *pzwritep++ = ' ';
	    uii = missing_ct_by_cluster[clidx];
            pzwritep = uint32_writew8x(pzwritep, uii, ' ');
	    umm = cur_cluster_sizes[clidx];
	    pzwritep = uint32_writew8x(pzwritep, umm, ' ');
	    umm -= oblig_missing_ct_by_cluster[clidx];
	    pzwritep = uint32_writew8x(pzwritep, umm, ' ');
            pzwritep = double_g_writewx4(pzwritep, ((double)((int32_t)uii)) / ((double)((int32_t)umm)), 8);
	    append_binary_eoln(&pzwritep);
	    if (flex_pzwrite(&ps, &pzwritep)) {
	      goto write_missingness_reports_ret_WRITE_FAIL;
	    }
	  }
	}
        marker_uidx++;
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	  if (marker_uidx < chrom_end) {
	    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	      goto write_missingness_reports_ret_WRITE_FAIL;
	    }
	  }
	}
      } while (marker_uidx < chrom_end);
    }
  }
  if (flex_pzwrite_close_null(&ps, pzwritep)) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  outname_end[1] = 'i';
  if (flex_pzwrite_init(output_gz, outname, overflow_buf, 0, &ps)) {
    goto write_missingness_reports_ret_OPEN_FAIL;
  }
  pzwritep = (char*)overflow_buf;
  sprintf(tbuf, "%%%us %%%us MISS_PHENO   N_MISS   N_GENO   F_MISS" EOLN_STR, plink_maxfid, plink_maxiid);
  pzwritep += sprintf(pzwritep, tbuf, "FID", "IID");
  do {
    sample_uidx = next_unset_unsafe(sample_exclude, sample_uidx);
    sample_uidx_stop = next_set(sample_exclude, sample_uidx, unfiltered_sample_ct);
    sample_idx += sample_uidx_stop - sample_uidx;
    do {
      cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      cptr2 = (char*)memchr(cptr, '\t', max_sample_id_len);
      slen = (uintptr_t)(cptr2 - cptr);
      pzwritep = memseta(pzwritep, 32, plink_maxfid - slen);
      pzwritep = memcpyax(pzwritep, cptr, slen, ' ');
      pzwritep = fw_strcpy(plink_maxiid, &(cptr2[1]), pzwritep);
      pzwritep = memseta(pzwritep, 32, 10);
      *pzwritep++ = 'Y' - (is_set(pheno_nm, sample_uidx) * 11);
      *pzwritep++ = ' ';
      uii = missing_cts[sample_uidx];
      ukk = is_set(sex_male, sample_uidx);
      ujj = marker_ct_nony + (ukk * marker_ct_y);
      if (om_sample_lookup) {
	umm = om_sample_lookup[sample_uidx];
	if (umm != 0xffffffffU) {
          umm = om_cluster_ref_cts[umm] + ukk * om_cluster_ref_cts[umm + om_cluster_ct];
	  uii -= umm;
	  ujj -= umm;
	}
      }
      pzwritep = uint32_writew8x(pzwritep, uii, ' ');
      pzwritep = uint32_writew8x(pzwritep, ujj, ' ');
      pzwritep = double_g_writewx4(pzwritep, ((double)((int32_t)uii)) / ((double)((int32_t)ujj)), 8);
      append_binary_eoln(&pzwritep);
      if (flex_pzwrite(&ps, &pzwritep)) {
	goto write_missingness_reports_ret_WRITE_FAIL;
      }
    } while (++sample_uidx < sample_uidx_stop);
  } while (sample_idx < sample_ct);

  if (flex_pzwrite_close_null(&ps, pzwritep)) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  *outname_end = '\0';
  LOGPRINTFWW("--missing: Sample missing data report written to %s.imiss%s, and variant-based %smissing data report written to %s.lmiss%s.\n", outname, output_gz? ".gz" : "", cluster_ct? "cluster-stratified " : "", outname, output_gz? ".gz" : "");
  while (0) {
  write_missingness_reports_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_missingness_reports_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_missingness_reports_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  write_missingness_reports_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  flex_pzwrite_close_cond(&ps, pzwritep);
  return retval;
}

int32_t hardy_report_write_line(Pigz_state* ps_ptr, char** pzwritep_ptr, char* prefix_buf, uint32_t prefix_len, uint32_t reverse, uint32_t ll_ct, uint32_t lh_ct, uint32_t hh_ct, char* midbuf_ptr, double pval, double output_min_p) {
  char* pzwritep = *pzwritep_ptr;
  char wbuf[48];
  char* cptr;
  uint32_t denom;
  double drecip;
  double minor_freq;
  pzwritep = memcpya(pzwritep, prefix_buf, prefix_len);
  if (reverse) {
    cptr = uint32_write(uint32_writex(uint32_writex(wbuf, hh_ct, '/'), lh_ct, '/'), ll_ct);
  } else {
    cptr = uint32_write(uint32_writex(uint32_writex(wbuf, ll_ct, '/'), lh_ct, '/'), hh_ct);
  }
  pzwritep = fw_strcpyn(20, cptr - wbuf, wbuf, pzwritep);
  *pzwritep++ = ' ';
  denom = (ll_ct + lh_ct + hh_ct) * 2;
  if (denom) {
    drecip = 1.0 / ((double)denom);
    minor_freq = (2 * ll_ct + lh_ct) * drecip;
    pzwritep = double_g_writewx4(double_g_writewx4x(double_g_writewx4x(pzwritep, (lh_ct * 2) * drecip, 8, ' '), minor_freq * (2 * hh_ct + lh_ct) * drecip * 2, 8, ' '), MAXV(pval, output_min_p), 12);
  } else {
    pzwritep = memcpya(pzwritep, "     nan      nan           NA", 30);
  }
  append_binary_eoln(&pzwritep);
  if (flex_pzwrite(ps_ptr, &pzwritep)) {
    return 1;
  }
  *pzwritep_ptr = pzwritep;
  return 0;
}

int32_t hardy_report(char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, uint32_t hwe_modifier, uint32_t nonfounders, int32_t* hwe_ll_cases, int32_t* hwe_lh_cases, int32_t* hwe_hh_cases, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs, uint32_t pheno_nm_ct, uintptr_t* pheno_c, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  char* pzwritep = NULL;
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  uint32_t hwe_midp = hwe_modifier & HWE_MIDP;
  uint32_t output_gz = (hwe_modifier / HWE_GZ) & 1;
  int32_t retval = 0;
  uint32_t skip_chrom = 0;
  uint32_t pct = 0;
  Pigz_state ps;
  uint32_t prefix_len;
  uint32_t loop_end;
  uint32_t uii;
  uint32_t report_type;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_mt;
  uint32_t is_haploid;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t reverse;
  double* p_values;
  unsigned char* overflow_buf;
  char* writebuf;
  char* cptr0;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* cptr5;
  pzwrite_init_null(&ps);
  if (pheno_nm_ct) {
    report_type = pheno_c? 0 : 1;
  } else {
    report_type = 2;
  }
  uii = report_type? 1 : 3;
  if (wkspace_alloc_uc_checked(&overflow_buf, PIGZ_BLOCK_SIZE + 2 * max_marker_allele_len + MAXLINELEN) ||
      wkspace_alloc_d_checked(&p_values, uii * marker_ct * sizeof(double)) ||
      wkspace_alloc_c_checked(&writebuf, 2 * max_marker_allele_len + MAXLINELEN)) {
    goto hardy_report_ret_NOMEM;
  }

  // todo: multithread?
  if (report_type) {
    for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      p_values[marker_idx] = SNPHWE2(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], hwe_midp);
    }
  } else {
    for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      p_values[marker_idx * 3] = SNPHWE2(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], hwe_midp);
      p_values[marker_idx * 3 + 1] = SNPHWE2(hwe_lh_cases[marker_uidx], hwe_ll_cases[marker_uidx], hwe_hh_cases[marker_uidx], hwe_midp);
      p_values[marker_idx * 3 + 2] = SNPHWE2(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx], hwe_midp);
    }
  }
  marker_uidx = 0;
  marker_idx = 0;

  memcpy(outname_end, output_gz? ".hwe.gz" : ".hwe", output_gz? 8 : 5);
  if (flex_pzwrite_init(output_gz, outname, overflow_buf, 0, &ps)) {
    goto hardy_report_ret_OPEN_FAIL;
  }
  pzwritep = (char*)overflow_buf;

  LOGPRINTFWW5("--hardy: Writing Hardy-Weinberg report (%s) to %s ... ", nonfounders? "all samples" : "founders only", outname);
  fputs("0%", stdout);
  fflush(stdout);
  sprintf(writebuf, " CHR %%%us     TEST   A1   A2                 GENO   O(HET)   E(HET)            P " EOLN_STR, plink_maxsnp);
  pzwritep += sprintf(pzwritep, writebuf, "SNP");

  chrom_fo_idx = 0;
  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
  skip_chrom = (is_haploid && (!is_x)) || is_mt;
  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx]));
  *cptr0++ = ' ';
  cptr = &(cptr0[10 + plink_maxsnp]);
  prefix_len = 10 + ((uintptr_t)(cptr - writebuf));
  if (report_type) {
    if (report_type == 1) {
      memcpy(&(cptr0[plink_maxsnp]), "  ALL(QT)           ", 20);
    } else {
      memcpy(&(cptr0[plink_maxsnp]), "  ALL(NP)           ", 20);
    }
    cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
    for (; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100LLU;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  skip_chrom = (is_haploid && (!is_x)) || is_mt;
	  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx]));
	  *cptr0++ = ' ';
	  cptr = &(cptr0[10 + plink_maxsnp]);
	  prefix_len = 10 + ((uintptr_t)(cptr - writebuf));
	  if (report_type == 1) {
	    memcpy(&(cptr0[plink_maxsnp]), "  ALL(QT)           ", 20);
	  } else {
	    memcpy(&(cptr0[plink_maxsnp]), "  ALL(NP)           ", 20);
	  }
	  cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
	}
        if (skip_chrom) {
	  continue;
	}
	fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), cptr0);
	reverse = IS_SET(marker_reverse, marker_uidx);
	cptr3 = marker_allele_ptrs[2 * marker_uidx];
	cptr4 = marker_allele_ptrs[2 * marker_uidx + 1];
	cptr5 = fw_strcpy(4, cptr3, cptr);
	*cptr5 = ' ';
	cptr5 = fw_strcpy(4, cptr4, &(cptr5[1]));
	*cptr5 = ' ';
	prefix_len = 1 + (cptr5 - writebuf);
	if (hardy_report_write_line(&ps, &pzwritep, writebuf, prefix_len, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[marker_idx], output_min_p)) {
	  goto hardy_report_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else {
    memset(&(cptr0[plink_maxsnp]), 32, 20);
    cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
    for (; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100LLU;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  skip_chrom = (is_haploid && (!is_x)) || is_mt;
	  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx]));
	  *cptr0++ = ' ';
          memset(&(cptr0[plink_maxsnp]), 32, 20);
	  cptr = &(cptr0[10 + plink_maxsnp]);
	  cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
	  prefix_len = 10 + ((uintptr_t)(cptr - writebuf));
	}
	if (skip_chrom) {
	  continue;
	}
	fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), cptr0);
	memcpy(&(cptr0[4 + plink_maxsnp]), "  ALL", 5);
	reverse = IS_SET(marker_reverse, marker_uidx);
	cptr3 = marker_allele_ptrs[2 * marker_uidx];
	cptr4 = marker_allele_ptrs[2 * marker_uidx + 1];
	cptr5 = fw_strcpy(4, cptr3, cptr);
	*cptr5 = ' ';
	cptr5 = fw_strcpy(4, cptr4, &(cptr5[1]));
	*cptr5 = ' ';
	prefix_len = 1 + (cptr5 - writebuf);
	if (hardy_report_write_line(&ps, &pzwritep, writebuf, prefix_len, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[3 * marker_idx], output_min_p)) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(cptr0[7 + plink_maxsnp]), "FF", 2);
	if (hardy_report_write_line(&ps, &pzwritep, writebuf, prefix_len, reverse, hwe_ll_cases[marker_uidx], hwe_lh_cases[marker_uidx], hwe_hh_cases[marker_uidx], cptr2, p_values[3 * marker_idx + 1], output_min_p)) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(cptr0[4 + plink_maxsnp]), "UN", 2);
	if (hardy_report_write_line(&ps, &pzwritep, writebuf, prefix_len, reverse, hwe_lls[marker_uidx], hwe_lhs[marker_uidx], hwe_hhs[marker_uidx], cptr2, p_values[3 * marker_idx + 2], output_min_p)) {
	  goto hardy_report_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  }
  fputs("\b\b\b", stdout);
  logprint("done.\n");

  while (0) {
  hardy_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  hardy_report_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  hardy_report_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  flex_pzwrite_close_cond(&ps, pzwritep);
  wkspace_reset(wkspace_mark);
  return retval;
}

uint32_t enforce_hwe_threshold(double hwe_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, uint32_t hwe_modifier, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs, Chrom_info* chrom_info_ptr) {
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t removed_ct = 0;
  uint32_t hwe_all = hwe_modifier & HWE_THRESH_ALL;
  uint32_t hwe_thresh_midp = hwe_modifier & HWE_THRESH_MIDP;
  uint32_t min_obs = 0xffffffffU;
  uint32_t max_obs = 0;
  int32_t mt_code = chrom_info_ptr->mt_code;
  uint32_t mt_start = 0;
  uint32_t mt_end = 0;
  uint32_t markers_done;
  uint32_t cur_obs;
  hwe_thresh *= 1 + SMALL_EPSILON;
  if (hwe_all) {
    hwe_lhs = hwe_lh_allfs;
    hwe_lls = hwe_ll_allfs;
    hwe_hhs = hwe_hh_allfs;
  }
  if ((mt_code != -1) && is_set(chrom_info_ptr->chrom_mask, mt_code)) {
    mt_start = chrom_info_ptr->chrom_start[(uint32_t)mt_code];
    mt_end = chrom_info_ptr->chrom_end[(uint32_t)mt_code];
  }
  if (hwe_thresh_midp) {
    for (markers_done = 0; markers_done < marker_ct; marker_uidx++, markers_done++) {
      next_unset_unsafe_ck(marker_exclude, &marker_uidx);
      if ((marker_uidx < mt_end) && (marker_uidx >= mt_start)) {
        continue;
      }
      if (SNPHWE_midp_t(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx], hwe_thresh)) {
	SET_BIT(marker_exclude, marker_uidx);
	removed_ct++;
      }
      cur_obs = hwe_lhs[marker_uidx] + hwe_lls[marker_uidx] + hwe_hhs[marker_uidx];
      if (cur_obs < min_obs) {
	min_obs = cur_obs;
      }
      if (cur_obs > max_obs) {
	max_obs = cur_obs;
      }
    }
  } else {
    for (markers_done = 0; markers_done < marker_ct; marker_uidx++, markers_done++) {
      next_unset_unsafe_ck(marker_exclude, &marker_uidx);
      if ((marker_uidx < mt_end) && (marker_uidx >= mt_start)) {
        continue;
      }
      if (SNPHWE_t(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx], hwe_thresh)) {
	SET_BIT(marker_exclude, marker_uidx);
	removed_ct++;
      }
      cur_obs = hwe_lhs[marker_uidx] + hwe_lls[marker_uidx] + hwe_hhs[marker_uidx];
      if (cur_obs < min_obs) {
	min_obs = cur_obs;
      }
      if (cur_obs > max_obs) {
	max_obs = cur_obs;
      }
    }
  }
  if (((uint64_t)max_obs) * 9 > ((uint64_t)min_obs) * 10) {
    logprint("Warning: --hwe observation counts vary by more than 10%.  Consider using\n--geno, and/or applying different p-value thresholds to distinct subsets of\nyour data.\n");
  }
  if (marker_ct == removed_ct) {
    logprint("Error: All variants removed due to Hardy-Weinberg exact test (--hwe).\n");
    return 1;
  }
  LOGPRINTF("--hwe: %u variant%s removed due to Hardy-Weinberg exact test.\n", removed_ct, (removed_ct == 1)? "" : "s");
  *marker_exclude_ct_ptr += removed_ct;
  return 0;
}

uint32_t enforce_minor_allele_thresholds(double min_maf, double max_maf, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* ac_excl_bitfield, uintptr_t* marker_exclude_ct_ptr, double* set_allele_freqs) {
  uint32_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t removed_ct = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  double dxx;
  if ((min_maf != 0.0) || (max_maf != 0.5)) {
    min_maf *= 1 - SMALL_EPSILON;
    max_maf *= 1 + SMALL_EPSILON;
    while (markers_done < marker_ct) {
      marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
      marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
      markers_done += marker_uidx_stop - marker_uidx;
      do {
	dxx = get_maf(set_allele_freqs[marker_uidx]);
	if ((dxx < min_maf) || (dxx > max_maf)) {
	  SET_BIT(marker_exclude, marker_uidx);
	}
      } while (++marker_uidx < marker_uidx_stop);
    }
  }
  if (ac_excl_bitfield) {
    bitfield_or(marker_exclude, ac_excl_bitfield, unfiltered_marker_ctl);
  }
  removed_ct = popcount_longs(marker_exclude, unfiltered_marker_ctl) - (*marker_exclude_ct_ptr);
  if (marker_ct == removed_ct) {
    logprint("Error: All variants removed due to minor allele threshold(s)\n(--maf/--max-maf/--mac/--max-mac).\n");
    return 1;
  }
  LOGPRINTFWW("%u variant%s removed due to minor allele threshold(s) (--maf/--max-maf/--mac/--max-mac).\n", removed_ct, (removed_ct == 1)? "" : "s");
  *marker_exclude_ct_ptr += removed_ct;
  return 0;
}

void enforce_min_bp_space(int32_t min_bp_space, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, uintptr_t* marker_exclude_ct_ptr, Chrom_info* chrom_info_ptr) {
  uint32_t marker_ct = unfiltered_marker_ct - (uint32_t)(*marker_exclude_ct_ptr);
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t removed_ct = 0;
  uint32_t chrom_end = 0;
  uint32_t marker_uidx = next_unset(marker_exclude, 0, unfiltered_marker_ct);
  uint32_t chrom_fo_idx_p1 = 0;
  uint32_t marker_uidx_stop;
  int32_t last_pos;
  int32_t cur_pos;
  for (chrom_fo_idx_p1 = 1; chrom_fo_idx_p1 <= chrom_ct; chrom_fo_idx_p1++) {
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx_p1];
    if (marker_uidx >= chrom_end) {
      continue;
    }
    last_pos = -2147483647;
    do {
      marker_uidx_stop = next_set(marker_exclude, marker_uidx, chrom_end);
      do {
        cur_pos = marker_pos[marker_uidx];
        if (cur_pos < last_pos + min_bp_space) {
          SET_BIT(marker_exclude, marker_uidx);
	  removed_ct++;
	} else {
	  last_pos = cur_pos;
	}
      } while (++marker_uidx < marker_uidx_stop);
      marker_uidx = next_unset(marker_exclude, marker_uidx, unfiltered_marker_ct);
    } while (marker_uidx < chrom_end);
  }
  LOGPRINTF("--bp-space: %u variant%s removed (%u remaining).\n", removed_ct, (removed_ct == 1)? "" : "s", marker_ct - removed_ct);
  *marker_exclude_ct_ptr += removed_ct;
}
