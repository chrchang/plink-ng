#include "plink_filter.h"
#include "plink_stats.h"

void filter_init(Oblig_missing_info* om_ip, Mendel_info* me_ip) {
  om_ip->cluster_ct = 0;
  om_ip->entry_ct = 0;
  om_ip->entries = NULL;
  om_ip->cluster_ref_cts = NULL;
  om_ip->indiv_lookup = NULL;
  om_ip->marker_fname = NULL;
  om_ip->indiv_fname = NULL;
  me_ip->modifier = 0;
  me_ip->max_trio_error = 1.0;
  me_ip->max_var_error = 1.0;
  me_ip->exclude_one_ratio = 0.0;
}

void filter_cleanup(Oblig_missing_info* om_ip) {
  if (om_ip->marker_fname) {
    free_cond(om_ip->entries);
    free_cond(om_ip->cluster_ref_cts);
    free_cond(om_ip->indiv_lookup);
    free_cond(om_ip->marker_fname);
    free_cond(om_ip->indiv_fname);
    om_ip->marker_fname = NULL;
  }
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

int32_t load_oblig_missing(FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, char* sorted_person_ids, uintptr_t sorted_indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* sex_male, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip) {
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
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + BITCT2 - 1) / BITCT2;
  uintptr_t sorted_indiv_ctl = (sorted_indiv_ct + BITCT - 1) / BITCT;
  uintptr_t topsize = 0;
  uintptr_t max_cluster_id_len = 0;
  uintptr_t possible_distinct_ct = 0;
  uintptr_t missing_cluster_ct = 0;
  uintptr_t y_start = 0;
  uintptr_t y_end = 0;
  int32_t y_code = chrom_info_ptr->y_code;
  uint32_t y_present = ((y_code != -1) && is_set(chrom_info_ptr->chrom_mask, y_code))? 1 : 0;
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
  uint32_t* indiv_lookup;
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
  uint32_t indiv_uidx;
  uint32_t cluster_idx;
  uint32_t slen;
  int32_t ii;

  if (y_present) {
    y_start = chrom_info_ptr->chrom_start[(uint32_t)y_code];
    y_end = chrom_info_ptr->chrom_end[(uint32_t)y_code];
  }
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto load_oblig_missing_ret_NOMEM;
  }
  loadbuf_end = &(loadbuf[unfiltered_indiv_ctl2]);
  if (fopen_checked(&infile, om_ip->indiv_fname, "r")) {
    goto load_oblig_missing_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';

  // two-pass load, same as load_clusters()
  // use loadbuf as duplicate IID detector
  fill_ulong_zero(loadbuf, sorted_indiv_ctl);
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Pathologically long line in %s.\n", om_ip->indiv_fname);
      goto load_oblig_missing_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(idbuf, sorted_person_ids, max_person_id_len, sorted_indiv_ct, bufptr, &bufptr2, &ii)) {
      sprintf(logbuf, "Error: Missing tokens in %s line.\n", om_ip->indiv_fname);
      goto load_oblig_missing_ret_INVALID_FORMAT_2;
    }
    if (ii != -1) {
      if (is_set(loadbuf, ii)) {
        sprintf(logbuf, "Error: Duplicate individual %s in --oblig-missing file.\n", idbuf);
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
    LOGPRINTF("Warning: --oblig-missing ignored, since no valid blocks were defined in\n%s.\n", om_ip->indiv_fname);
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
  wkspace_reset((unsigned char*)cluster_ids);
  cluster_ids = (char*)wkspace_alloc(cluster_ct * max_cluster_id_len);
  cluster_mct = cluster_ct * (y_present + 1);
  indiv_lookup = (uint32_t*)malloc(unfiltered_indiv_ct * sizeof(int32_t));
  if (!indiv_lookup) {
    goto load_oblig_missing_ret_NOMEM;
  }
  om_ip->indiv_lookup = indiv_lookup;
  if (wkspace_alloc_ui_checked(&cluster_sizes, cluster_mct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&cluster_zmask2s, cluster_mct * unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto load_oblig_missing_ret_NOMEM;
  }
  fill_uint_zero(cluster_sizes, cluster_mct);
  fill_uint_one(indiv_lookup, unfiltered_indiv_ct);
  fill_ulong_zero(cluster_zmask2s, cluster_mct * unfiltered_indiv_ctl2);

  // second pass
  rewind(infile);
  while (fgets(tbuf, MAXLINELEN, infile)) {
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bsearch_read_fam_indiv(idbuf, sorted_person_ids, max_person_id_len, sorted_indiv_ct, bufptr, &bufptr2, &ii);
    if (ii == -1) {
      continue;
    }
    indiv_uidx = indiv_id_map[(uint32_t)ii];    
    slen = strlen_se(bufptr2);
    // guaranteed to succeed
    ii = bsearch_str(bufptr2, slen, cluster_ids, max_cluster_id_len, cluster_ct);
    set_bit(&(cluster_zmask2s[((uintptr_t)((uint32_t)ii)) * unfiltered_indiv_ctl2]), indiv_uidx * 2);
    cluster_sizes[(uint32_t)ii] += 1;
    indiv_lookup[indiv_uidx] = (uint32_t)ii;
  }
  if (fclose_null(&infile)) {
    goto load_oblig_missing_ret_READ_FAIL;
  }
  if (y_present) {
    vec_include_init(unfiltered_indiv_ct, loadbuf, sex_male);
    cur_cluster_zmask2 = cluster_zmask2s;
    ulptr = &(cur_cluster_zmask2[cluster_ct * unfiltered_indiv_ctl2]);
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
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Pathologically long line in %s.\n", om_ip->marker_fname);
      goto load_oblig_missing_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr2 = item_endnn(bufptr);
    ii = bsearch_str(bufptr, (uintptr_t)(bufptr2 - bufptr), sorted_marker_ids, max_marker_id_len, marker_ct);
    if (ii != -1) {
      marker_uidx = marker_id_map[(uint32_t)ii];
      bufptr = skip_initial_spaces(bufptr2);
      if (is_eoln_kns(*bufptr)) {
        sprintf(logbuf, "Error: Missing tokens in %s line.\n", om_ip->marker_fname);
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
    LOGPRINTF("Warning: %" PRIuPTR " entr%s in %s had block IDs missing from\n%s.\n", missing_cluster_ct, (missing_cluster_ct == 1)? "y" : "ies", om_ip->marker_fname, om_ip->indiv_fname);
  }
  om_ip->entry_ct = (uintptr_t)(zc_entries_end - zc_entries);
  if (!om_ip->entry_ct) {
    LOGPRINTF("Warning: --oblig-missing ignored, since %s had no valid entries.\n", om_ip->marker_fname);
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

  loadbuf[unfiltered_indiv_ctl2 - 1] = 0;
  marker_uidx = unfiltered_marker_ct; // forces initial load
  do {
    ullii = *zc_entries++;
    if ((ullii >> 32) != marker_uidx) {
      marker_uidx = (ullii >> 32);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto load_oblig_missing_ret_READ_FAIL;
      }
      if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto load_oblig_missing_ret_READ_FAIL;
      }
      // no need for het haploid handling here
    }
    loadbuf_ptr = loadbuf;
    cluster_idx = (uint32_t)ullii;
    cur_cluster_zmask2 = &(cluster_zmask2s[cluster_idx * unfiltered_indiv_ctl2]);
    do {
      ulii = *cur_cluster_zmask2++;
      if (((*loadbuf_ptr) ^ ulii) & (ulii * 3)) {
	sprintf(logbuf, "Error: Nonmissing --oblig-missing genotype at marker %s.\n(To force it to missing, use --zero-cluster.)\n", &(marker_ids[marker_uidx * max_marker_id_len]));
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

int32_t filter_indivs_file(char* filtername, char* sorted_person_ids, uintptr_t sorted_ids_len, uintptr_t max_person_id_len, uint32_t* id_map, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* filtervals_flattened, uint32_t mfilter_col) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t include_ct = 0;
  uintptr_t max_filterval_len = 0;
  uint32_t filterval_ct = 0;
  int32_t retval = 0;
  char* sorted_filtervals;
  uintptr_t* indiv_exclude_new;
  char* id_buf;
  char* bufptr;
  uint32_t filterval_idx;
  uint32_t slen;
  int32_t person_idx;
  bufptr = filtervals_flattened;
  do {
    filterval_ct++;
    slen = strlen(bufptr);
    if (slen >= max_filterval_len) {
      max_filterval_len = slen + 1;
    }
    bufptr = &(bufptr[slen + 1]);
  } while (*bufptr);
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&indiv_exclude_new, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_c_checked(&sorted_filtervals, filterval_ct * max_filterval_len)) {
    goto filter_indivs_file_ret_NOMEM;
  }
  fill_all_bits(indiv_exclude_new, unfiltered_indiv_ct);
  bufptr = filtervals_flattened;
  for (filterval_idx = 0; filterval_idx < filterval_ct; filterval_idx++) {
    slen = strlen(bufptr) + 1;
    memcpy(&(sorted_filtervals[filterval_idx * max_filterval_len]), bufptr, slen);
    bufptr = &(bufptr[slen]);
  }
  qsort(sorted_filtervals, filterval_ct, max_filterval_len, strcmp_casted);

  if (fopen_checked(&infile, filtername, "r")) {
    goto filter_indivs_file_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in --keep/--remove file.\n");
      goto filter_indivs_file_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(id_buf, sorted_person_ids, max_person_id_len, sorted_ids_len, bufptr, &bufptr, &person_idx)) {
      goto filter_indivs_file_ret_INVALID_FORMAT;
    }
    if (person_idx != -1) {
      person_idx = id_map[(uint32_t)person_idx];
      if (!is_set(indiv_exclude, person_idx)) {
	if (mfilter_col > 1) {
	  bufptr = next_item_mult(bufptr, mfilter_col - 1);
	}
	if (no_more_items_kns(bufptr)) {
	  goto filter_indivs_file_ret_INVALID_FORMAT;
	}
	if (bsearch_str(bufptr, strlen_se(bufptr), sorted_filtervals, max_filterval_len, filterval_ct) != -1) {
	  if (is_set(indiv_exclude_new, person_idx)) {
	    clear_bit(indiv_exclude_new, person_idx);
	    include_ct++;
	  }
	}
      }
    }
  }
  if (!feof(infile)) {
    goto filter_indivs_file_ret_READ_FAIL;
  }
  if (!include_ct) {
    LOGPRINTF("Error: All %s excluded by --filter.\n", g_species_plural);
    goto filter_indivs_file_ret_ALL_SAMPLES_EXCLUDED;
  }
  LOGPRINTF("--filter: %" PRIuPTR " %s remaining.\n", include_ct, species_str(include_ct));
  memcpy(indiv_exclude, indiv_exclude_new, unfiltered_indiv_ctl * sizeof(intptr_t));
  *indiv_exclude_ct_ptr = unfiltered_indiv_ct - include_ct;

  while (0) {
  filter_indivs_file_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  filter_indivs_file_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  filter_indivs_file_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  filter_indivs_file_ret_INVALID_FORMAT:
    logprint("Error: Too few columns in --filter file line.\n");
  filter_indivs_file_ret_INVALID_FORMAT_2:
    retval = RET_INVALID_FORMAT;
    break;
  filter_indivs_file_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  return retval;
}

void filter_indivs_bitfields(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, uintptr_t* orfield, int32_t orfield_flip, uintptr_t* ornot) {
  // indiv_exclude := indiv_exclude | orfield | (~ornot) if !orfield_flip
  //               := indiv_exclude | (~orfield) | (~ornot) otherwise
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t* ieptr = indiv_exclude;
  uintptr_t* ieend = &(indiv_exclude[unfiltered_indiv_ctl]);
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
  zero_trailing_bits(indiv_exclude, unfiltered_indiv_ct);
  *indiv_exclude_ct_ptr = popcount_longs(indiv_exclude, unfiltered_indiv_ctl);
}

int32_t mind_filter(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double mind_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* person_ids, uintptr_t max_person_id_len, uintptr_t* sex_male, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uint32_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ct2l = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t marker_idx = 0;
  uintptr_t y_start = 0;
  uintptr_t y_end = 0;
  uintptr_t* indiv_male_include2 = NULL;
  uint32_t indiv_exclude_ct = *indiv_exclude_ct_ptr;
  uint32_t indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  uint32_t indiv_uidx = 0;
  uint32_t indiv_idx = 0;
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
  uint32_t* indiv_lookup;
  uintptr_t marker_uidx;
  uint32_t cluster_ct;
  uint32_t cur_marker_ct;
  uint32_t indiv_uidx_stop;
  uint32_t uii;
  uint32_t ujj;
  uintptr_t ulii;

  if (y_present) {
    y_start = chrom_info_ptr->chrom_start[(uint32_t)y_code];
    y_end = chrom_info_ptr->chrom_end[(uint32_t)y_code];
    if (wkspace_alloc_ul_checked(&indiv_male_include2, unfiltered_indiv_ct2l * sizeof(intptr_t))) {
      goto mind_filter_ret_NOMEM;
    }
    vec_include_init(unfiltered_indiv_ct, indiv_male_include2, sex_male);
    nony_marker_ct = marker_ct - (y_end - y_start - popcount_bit_idx(marker_exclude, y_start, y_end));
  }
  if (wkspace_alloc_ui_checked(&missing_cts, unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ct2l * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&newly_excluded, unfiltered_indiv_ctl * sizeof(int32_t))) {
    goto mind_filter_ret_NOMEM;
  }
  loadbuf[unfiltered_indiv_ct2l - 1] = 0;
  fill_uint_zero(missing_cts, unfiltered_indiv_ct);
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto mind_filter_ret_READ_FAIL;
  }
  ujj = unfiltered_indiv_ct2l * BITCT2;
  marker_uidx = 0;
  for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto mind_filter_ret_READ_FAIL;
      }
    }
    if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
      goto mind_filter_ret_READ_FAIL;
    }
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
      mptr = indiv_male_include2;
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
  fill_ulong_zero(newly_excluded, unfiltered_indiv_ctl);
  if (!om_ip->entry_ct) {
    mind_int_thresh[0] = (int32_t)(mind_thresh * ((int32_t)nony_marker_ct) * (1 + SMALL_EPSILON));
    mind_int_thresh[1] = (int32_t)(mind_thresh * ((int32_t)marker_ct) * (1 + SMALL_EPSILON));
    do {
      indiv_uidx = next_unset_unsafe(indiv_exclude, indiv_uidx);
      indiv_uidx_stop = next_set(indiv_exclude, indiv_uidx, unfiltered_indiv_ct);
      indiv_idx += indiv_uidx_stop - indiv_uidx;
      do {
	if (missing_cts[indiv_uidx] > mind_int_thresh[is_set(sex_male, indiv_uidx)]) {
	  SET_BIT(newly_excluded, indiv_uidx);
	  removed_ct++;
	}
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (indiv_idx < indiv_ct);
  } else {
    indiv_lookup = om_ip->indiv_lookup;
    cluster_ref_cts = om_ip->cluster_ref_cts;
    cluster_ct = om_ip->cluster_ct;
    do {
      indiv_uidx = next_unset_unsafe(indiv_exclude, indiv_uidx);
      indiv_uidx_stop = next_set(indiv_exclude, indiv_uidx, unfiltered_indiv_ct);
      indiv_idx += indiv_uidx_stop - indiv_uidx;
      do {
	uii = indiv_lookup[indiv_uidx];
	ujj = 0;
	if (is_set(sex_male, indiv_uidx)) {
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
	if ((missing_cts[indiv_uidx] - ujj) > (uint32_t)((int32_t)(mind_thresh * ((int32_t)(cur_marker_ct - ujj)) * (1 + SMALL_EPSILON)))) {
	  SET_BIT(newly_excluded, indiv_uidx);
	  removed_ct++;
	}
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (indiv_idx < indiv_ct);
  }
  if (removed_ct) {
    bitfield_or(indiv_exclude, newly_excluded, unfiltered_indiv_ctl);
    memcpy(outname_end, ".irem", 6);
    if (fopen_checked(&outfile, outname, "w")) {
      goto mind_filter_ret_OPEN_FAIL;
    }
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < removed_ct; indiv_idx++, indiv_uidx++) {
      next_set_unsafe_ck(newly_excluded, &indiv_uidx);
      fputs(&(person_ids[indiv_uidx * max_person_id_len]), outfile);
      putc('\n', outfile);
    }
    if (fclose_null(&outfile)) {
      goto mind_filter_ret_WRITE_FAIL;
    }
  }
  *indiv_exclude_ct_ptr += removed_ct;
  if (*indiv_exclude_ct_ptr == unfiltered_indiv_ct) {
    LOGPRINTF("Error: All %s removed due to missing genotype data (--mind).\nIDs written to %s.\n", g_species_plural, outname);
    goto mind_filter_ret_ALL_SAMPLES_EXCLUDED;
  }
  LOGPRINTF("%u %s removed due to missing genotype data (--mind).\n", removed_ct, species_str(removed_ct));
  if (removed_ct) {
    LOGPRINTF("ID%s written to %s.\n", (removed_ct == 1)? "" : "s", outname);
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
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
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
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 = loader2 & loader3;
  to_ct_hmaj1 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
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

static inline void single_marker_freqs_and_hwe(uintptr_t unfiltered_indiv_ctl2, uintptr_t* lptr, uintptr_t* indiv_include2, uintptr_t* founder_include2, uintptr_t* founder_ctrl_include2, uintptr_t* founder_case_include2, uintptr_t indiv_ct, uint32_t* ll_ctp, uint32_t* lh_ctp, uint32_t* hh_ctp, uint32_t indiv_f_ct, uint32_t* ll_ctfp, uint32_t* lh_ctfp, uint32_t* hh_ctfp, uint32_t hwe_or_geno_needed, uintptr_t indiv_f_ctl_ct, uint32_t* ll_hwep, uint32_t* lh_hwep, uint32_t* hh_hwep, int32_t hardy_needed, uintptr_t indiv_f_case_ct, uint32_t* ll_case_hwep, uint32_t* lh_case_hwep, uint32_t* hh_case_hwep) {
  // This is best understood from the bottom third up (which is the order it
  // was written).  It's way overkill for just determining genotype
  // frequencies, but a ruthlessly optimized version is needed for e.g.
  // permutation testing so we may as well get it working here.
  //
  // The idea, which underlies the IBS and LD pruners as well, is to obtain
  // multiple popcounts at once.  In the case of PLINK frequency/HWE
  // evaluation, we need 6-9 numbers:
  // homozyg minor, het, homozyg major counts for all individuals
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
  //               = indiv_ct - homozyg minor ct - het ct
  //
  // Thus, with the appropriate indiv_ct and these three popcounts, we can
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
  uintptr_t* lptr_end = &(lptr[unfiltered_indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
#ifdef __LP64__
  uintptr_t cur_decr = 120;
  uintptr_t* lptr_12x_end;
  unfiltered_indiv_ctl2 -= unfiltered_indiv_ctl2 % 12;
  while (unfiltered_indiv_ctl2 >= 120) {
  single_marker_freqs_and_hwe_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)indiv_include2, &tot_a, &tot_b, &tot_c);
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
    indiv_include2 = &(indiv_include2[cur_decr]);
    founder_include2 = &(founder_include2[cur_decr]);
    unfiltered_indiv_ctl2 -= cur_decr;
  }
  if (unfiltered_indiv_ctl2) {
    cur_decr = unfiltered_indiv_ctl2;
    goto single_marker_freqs_and_hwe_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[unfiltered_indiv_ctl2 - unfiltered_indiv_ctl2 % 12]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, indiv_include2, &tot_a, &tot_b, &tot_c);
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
    indiv_include2 = &(indiv_include2[12]);
    founder_include2 = &(founder_include2[12]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = *indiv_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    // N.B. because of the construction of indiv_include2, only even-numbered
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
  *ll_ctp = indiv_ct - tot_a - *lh_ctp;
  *hh_ctfp = tot_c_f;
  *lh_ctfp = tot_b_f - tot_c_f;
  *ll_ctfp = indiv_f_ct - tot_a_f - *lh_ctfp;
  if (hwe_or_geno_needed) {
    *hh_hwep = tot_c_hwe;
    *lh_hwep = tot_b_hwe - tot_c_hwe;
    *ll_hwep = indiv_f_ctl_ct - tot_a_hwe - *lh_hwep;
    if (hardy_needed) {
      *hh_case_hwep = tot_c_chwe;
      *lh_case_hwep = tot_b_chwe - tot_c_chwe;
      *ll_case_hwep = indiv_f_case_ct - tot_a_chwe - *lh_case_hwep;
    }
  }
}

static inline uint32_t nonmissing_present_diff(uintptr_t unfiltered_indiv_ctl2, uintptr_t* lptr, uintptr_t* indiv_include2, uintptr_t* indiv_male_include2) {
  // possible todo: write entries to .ynm file, using same format as .hh
  uintptr_t* lptr_end = &(lptr[unfiltered_indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  do {
    loader = *lptr++;
    loader2 = (*indiv_include2++) & (~(*indiv_male_include2++));
    // when really bored, check if compiler translates this into andnot
    // operations
    if ((~((~(loader >> 1)) & loader)) & loader2) {
      return 1;
    }
  } while (lptr < lptr_end);
  return 0;
}

static inline void haploid_single_marker_freqs(uintptr_t unfiltered_indiv_ct, uintptr_t unfiltered_indiv_ctl2, uintptr_t* lptr, uintptr_t* indiv_include2, uintptr_t* founder_include2, uintptr_t indiv_ct, uint32_t* ll_ctp, uint32_t* hh_ctp, uint32_t indiv_f_ct, uint32_t* ll_ctfp, uint32_t* hh_ctfp, uint32_t* hethap_incr_ptr) {
  // Here, we interpret heterozygotes as missing.
  // Nonmissing: (genotype ^ (~(genotype >> 1))) & 0x5555...
  // Homozygote major: (genotype & (genotype >> 1)) & 0x5555...
  uint32_t tot_a = 0;
  uint32_t tot_b = 0;
  uint32_t tot_hmaj = 0;
  uint32_t tot_nm_f = 0;
  uint32_t tot_hmaj_f = 0;
  uintptr_t* lptr_end = &(lptr[unfiltered_indiv_ctl2]);
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
  unfiltered_indiv_ctl2 -= unfiltered_indiv_ctl2 % 12;
  while (unfiltered_indiv_ctl2 >= 120) {
  single_marker_freqs_and_hwe_loop:
    lptr_12x_end = &(lptr[cur_decr]);
  // Given a buffer with PLINK binary genotypes for a single marker, let
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  // Then,
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  //               = indiv_ct - homozyg minor ct - het ct
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)indiv_include2, &tot_a, &tot_b, &tot_hmaj);
    freq_hwe_haploid_count_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_include2, &tot_nm_f, &tot_hmaj_f);
    lptr = lptr_12x_end;
    indiv_include2 = &(indiv_include2[cur_decr]);
    founder_include2 = &(founder_include2[cur_decr]);
    unfiltered_indiv_ctl2 -= cur_decr;
  }
  if (unfiltered_indiv_ctl2) {
    cur_decr = unfiltered_indiv_ctl2;
    goto single_marker_freqs_and_hwe_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[unfiltered_indiv_ctl2 - unfiltered_indiv_ctl2 % 12]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, indiv_include2, &tot_a, &tot_b, &tot_hmaj);
    freq_hwe_haploid_count_12(lptr, founder_include2, &tot_nm_f, &tot_hmaj_f);
    lptr = &(lptr[12]);
    indiv_include2 = &(indiv_include2[12]);
    founder_include2 = &(founder_include2[12]);
  }
#endif
  tot_nm = 2 * tot_hmaj + indiv_ct - tot_a - tot_b;
  hethap_incr = tot_b - tot_hmaj;
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader3 = loader >> 1;
    loader4 = *indiv_include2++;
    // different from tot_nm_f because of +indiv_ct above
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

int32_t calc_freqs_and_hwe(FILE* bedfile, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_exclude_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t* founder_info, int32_t nonfounders, int32_t maf_succ, double* set_allele_freqs, uintptr_t** marker_reverse_ptr, uintptr_t bed_offset, uint32_t hwe_needed, uint32_t hwe_all, uint32_t hardy_needed, double geno_thresh, uintptr_t* pheno_nm, uintptr_t* pheno_c, int32_t** hwe_lls_ptr, int32_t** hwe_lhs_ptr, int32_t** hwe_hhs_ptr, int32_t** hwe_ll_cases_ptr, int32_t** hwe_lh_cases_ptr, int32_t** hwe_hh_cases_ptr, int32_t** hwe_ll_allfs_ptr, int32_t** hwe_lh_allfs_ptr, int32_t** hwe_hh_allfs_ptr, int32_t** hwe_hapl_allfs_ptr, int32_t** hwe_haph_allfs_ptr, uintptr_t** geno_excl_bitfield_ptr, uint32_t* indiv_male_ct_ptr, uint32_t* indiv_f_ct_ptr, uint32_t* indiv_f_male_ct_ptr, uint32_t wt_needed, uintptr_t* topsize_ptr, double** marker_weights_ptr, double exponent, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t is_split_chrom, uint32_t* hh_exists_ptr) {
  FILE* hhfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + BITCT - 1) / BITCT;
  uintptr_t unfiltered_indiv_ctv2 = 2 * unfiltered_indiv_ctl;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  int32_t retval = 0;
  uint32_t pct = 1;
  uint32_t indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  double indiv_ct_recip = 1.0 / ((double)((int32_t)indiv_ct));
  // sum of nonmissing rates over all markers
  // rate is in [0, 1] for each marker, so sum is in [0, marker_ct].
  double nonmissing_rate_tot = 0.0;
  // track this to defend against Y chromosome/0 males pathological case
  uintptr_t nonmissing_rate_tot_max = marker_ct;
  uint32_t indiv_f_ct = indiv_ct;
  uintptr_t indiv_f_ctrl_ct = indiv_ct;
  uintptr_t indiv_f_case_ct = indiv_ct;
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
  uintptr_t* indiv_nonmale_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uintptr_t* founder_nonmale_include2 = NULL;
  uintptr_t* founder_ctrl_nonmale_include2 = NULL;
  uintptr_t* founder_male_include2 = NULL;
  uintptr_t* founder_case_include2 = NULL;
  uintptr_t* founder_case_nonmale_include2 = NULL;
  uintptr_t* geno_excl_bitfield = NULL;
  uint64_t* om_entry_ptr = NULL;
  double* marker_weights = NULL;
  uint32_t indiv_nonmale_ct = 0;
  uint32_t indiv_f_nonmale_ct = 0;
  uint32_t indiv_f_ctl_nonmale_ct = 0;
  uint32_t indiv_f_case_nonmale_ct = 0;
  uint64_t hethap_ct = 0;
  uint64_t cur_om_entry = 0;
  double male_ct_recip = 0;
  uint32_t* om_indiv_lookup;
  int32_t* hwe_hapl_allfs;
  int32_t* hwe_haph_allfs;
  uintptr_t* loadbuf;
  uintptr_t* indiv_include2;
  uintptr_t* founder_include2;
  uintptr_t* founder_ctrl_include2;
  uintptr_t* tmp_indiv_excl_mask;
  uintptr_t* tmp_indiv_excl_mask2;
  uintptr_t loop_end;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t ulii;
  uint32_t indiv_male_ct;
  uint32_t indiv_f_male_ct;
  uint32_t hethap_incr;
  uint32_t nonmales_needed;
  uint32_t males_needed;
  uint32_t uii;
  uint32_t ujj;
  double maf;
  double cur_genotyping_rate;
  uii = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(marker_reverse_ptr, uii * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  fill_ulong_zero(*marker_reverse_ptr, uii);
  if (wt_needed) {
    // this is a pretty ugly hack... but no worse than what preceded it, I
    // suppose
    marker_weights = (double*)top_alloc(topsize_ptr, CACHEALIGN(unfiltered_marker_ct * sizeof(double)));
    if (!marker_weights) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    wkspace_left -= *topsize_ptr;
    *marker_weights_ptr = marker_weights;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      marker_weights[marker_uidx] = -1.0;
    }
  }

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
  wkspace_mark = wkspace_base;
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  loadbuf[unfiltered_indiv_ctv2 - 2] = 0;
  loadbuf[unfiltered_indiv_ctv2 - 1] = 0;
  exclude_to_vec_include(unfiltered_indiv_ct, indiv_include2, indiv_exclude);
  ii = chrom_info_ptr->x_code;
  nonmales_needed = (!is_split_chrom) && (ii != -1) && is_set(chrom_info_ptr->chrom_mask, ii);
  ii = chrom_info_ptr->y_code;
  males_needed = nonmales_needed || ((!is_split_chrom) && (ii != -1) && is_set(chrom_info_ptr->chrom_mask, ii));
  if (wkspace_alloc_ul_checked(&indiv_male_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  memcpy(indiv_male_include2, indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
  vec_include_mask_in(unfiltered_indiv_ct, indiv_male_include2, sex_male);
  indiv_male_ct = popcount01_longs(indiv_male_include2, unfiltered_indiv_ctv2);
  if (indiv_male_ct) {
    male_ct_recip = 1.0 / ((double)((int32_t)indiv_male_ct));
  }
  *indiv_male_ct_ptr = indiv_male_ct;
  indiv_f_male_ct = indiv_male_ct;
  if (males_needed) {
    founder_male_include2 = indiv_male_include2;
    if (nonmales_needed) {
      if (wkspace_alloc_ul_checked(&indiv_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(indiv_nonmale_include2, indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      vec_include_mask_out_intersect(unfiltered_indiv_ct, indiv_nonmale_include2, sex_nm, sex_male);
      indiv_nonmale_ct = popcount01_longs(indiv_nonmale_include2, unfiltered_indiv_ctv2);
      founder_nonmale_include2 = indiv_nonmale_include2;
      founder_ctrl_nonmale_include2 = indiv_nonmale_include2;
    }
  }
  founder_include2 = indiv_include2;
  if (wkspace_alloc_ul_checked(&tmp_indiv_excl_mask, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  memcpy(tmp_indiv_excl_mask, indiv_exclude, unfiltered_indiv_ctl * sizeof(intptr_t));
  if (!nonfounders) {
    if (wkspace_alloc_ul_checked(&founder_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    bitfield_ornot(tmp_indiv_excl_mask, founder_info, unfiltered_indiv_ctl);
    zero_trailing_bits(tmp_indiv_excl_mask, unfiltered_indiv_ct);
    exclude_to_vec_include(unfiltered_indiv_ct, founder_include2, tmp_indiv_excl_mask);
    if (males_needed) {
      if (wkspace_alloc_ul_checked(&founder_male_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(founder_male_include2, indiv_male_include2, unfiltered_indiv_ctl * 2 * sizeof(intptr_t));
      vec_include_mask_in(unfiltered_indiv_ct, founder_male_include2, founder_info);
      indiv_f_male_ct = popcount01_longs(founder_male_include2, unfiltered_indiv_ctv2);
      if (nonmales_needed) {
	if (wkspace_alloc_ul_checked(&founder_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	  goto calc_freqs_and_hwe_ret_NOMEM;
	}
	memcpy(founder_nonmale_include2, indiv_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
	vec_include_mask_in(unfiltered_indiv_ct, founder_nonmale_include2, founder_info);
	indiv_f_nonmale_ct = popcount01_longs(founder_nonmale_include2, unfiltered_indiv_ctv2);
      }
    }
    founder_ctrl_include2 = founder_include2;
    indiv_f_ct = popcount_longs_exclude(founder_info, indiv_exclude, unfiltered_indiv_ctl);
    indiv_f_ctrl_ct = indiv_f_ct;
  } else {
    founder_ctrl_include2 = founder_include2;
  }

  // bugfix: this previously failed to initialize founder_ctrl_include2 and
  // founder_case_include2 properly if --hardy was used in a situation where
  // hwe_all would be set (e.g. all-case datasets).
  if ((!hwe_all) || hardy_needed) {
    if (wkspace_alloc_ul_checked(&founder_ctrl_include2, unfiltered_indiv_ctv2 *  sizeof(intptr_t)) ||
	wkspace_alloc_ul_checked(&tmp_indiv_excl_mask2, unfiltered_indiv_ctl * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    memcpy(tmp_indiv_excl_mask2, tmp_indiv_excl_mask, unfiltered_indiv_ctl);
    bitfield_ornot(tmp_indiv_excl_mask2, pheno_nm, unfiltered_indiv_ctl);
    bitfield_or(tmp_indiv_excl_mask2, pheno_c, unfiltered_indiv_ctl);
    zero_trailing_bits(tmp_indiv_excl_mask2, unfiltered_indiv_ct);
    // tmp_indiv_excl_mask2 is now set for each indiv who is excluded, or a
    // nonfounder, or is noncontrol.
    indiv_f_ctrl_ct = unfiltered_indiv_ct - popcount_longs(tmp_indiv_excl_mask2, unfiltered_indiv_ctl);
    exclude_to_vec_include(unfiltered_indiv_ct, founder_ctrl_include2, tmp_indiv_excl_mask2);
    if (nonmales_needed) {
      if (wkspace_alloc_ul_checked(&founder_ctrl_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(founder_ctrl_nonmale_include2, indiv_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      vec_include_mask_out(unfiltered_indiv_ct, founder_ctrl_nonmale_include2, tmp_indiv_excl_mask2);
      indiv_f_ctl_nonmale_ct = popcount01_longs(founder_ctrl_nonmale_include2, unfiltered_indiv_ctv2);
    }
    if (hardy_needed) {
      if (wkspace_alloc_ul_checked(&founder_case_include2, unfiltered_indiv_ctv2 *  sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      bitfield_ornot(tmp_indiv_excl_mask, pheno_nm, unfiltered_indiv_ctl);
      bitfield_ornot(tmp_indiv_excl_mask, pheno_c, unfiltered_indiv_ctl);
      zero_trailing_bits(tmp_indiv_excl_mask, unfiltered_indiv_ct);
      indiv_f_case_ct = unfiltered_indiv_ct - popcount_longs(tmp_indiv_excl_mask, unfiltered_indiv_ctl);
      exclude_to_vec_include(unfiltered_indiv_ct, founder_case_include2, tmp_indiv_excl_mask);
      if (nonmales_needed) {
	if (wkspace_alloc_ul_checked(&founder_case_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	  goto calc_freqs_and_hwe_ret_NOMEM;
	}
	memcpy(founder_case_nonmale_include2, indiv_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
	vec_include_mask_out(unfiltered_indiv_ct, founder_case_nonmale_include2, tmp_indiv_excl_mask);
	indiv_f_case_nonmale_ct = popcount01_longs(founder_case_nonmale_include2, unfiltered_indiv_ctv2);
      }
    }
  }
  if (om_ip->entry_ct) {
    om_entry_ptr = om_ip->entries;
    om_cluster_ct = om_ip->cluster_ct;
    om_indiv_lookup = om_ip->indiv_lookup;
    cur_om_entry = *om_entry_ptr;
    if (wkspace_alloc_ui_checked(&om_cluster_sizes, om_cluster_ct * 2 * sizeof(int32_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    fill_uint_zero(om_cluster_sizes, om_cluster_ct * 2);
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      uii = om_indiv_lookup[indiv_uidx];
      if (uii != 0xffffffffU) {
        om_cluster_sizes[uii] += 1;
        if (is_set(sex_male, indiv_uidx)) {
          om_cluster_sizes[uii + om_cluster_ct] += 1;
	}
      }
    }
  }

  *indiv_f_ct_ptr = indiv_f_ct;
  *indiv_f_male_ct_ptr = indiv_f_male_ct;
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
    is_haploid = (chrom_info_ptr->chrom_mask[0]) & 1;
    is_x = 0;
    is_y = 0;
    next_chrom_start = unfiltered_marker_ct;
  }
  for (; pct <= 100; pct++) {
    loop_end = ((uint64_t)pct * marker_ct) / 100LU;
    for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto calc_freqs_and_hwe_ret_READ_FAIL;
	}
      }
      if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
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
	single_marker_freqs_and_hwe(unfiltered_indiv_ctv2, loadbuf, indiv_include2, founder_include2, founder_ctrl_include2, founder_case_include2, indiv_ct, &ll_ct, &lh_ct, &hh_ct, indiv_f_ct, &ll_ctf, &lh_ctf, &hh_ctf, hwe_needed, indiv_f_ctrl_ct, &ll_hwe, &lh_hwe, &hh_hwe, hardy_needed, indiv_f_case_ct, &ll_case_hwe, &lh_case_hwe, &hh_case_hwe);
	hwe_ll_allfs[marker_uidx] = ll_ctf;
	hwe_lh_allfs[marker_uidx] = lh_ctf;
	hwe_hh_allfs[marker_uidx] = hh_ctf;
	uii = ll_ct + lh_ct + hh_ct;
	if (!cur_oblig_missing) {
	  cur_genotyping_rate = ((int32_t)uii) * indiv_ct_recip;
	} else {
	  if (indiv_ct - cur_oblig_missing) {
	    cur_genotyping_rate = ((int32_t)uii) / ((double)((int32_t)(indiv_ct - cur_oblig_missing)));
	  } else {
	    cur_genotyping_rate = 0;
	    nonmissing_rate_tot_max -= 1;
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
	    single_marker_freqs_and_hwe(unfiltered_indiv_ctv2, loadbuf, indiv_nonmale_include2, founder_nonmale_include2, founder_ctrl_nonmale_include2, founder_case_nonmale_include2, indiv_nonmale_ct, &ll_ct, &lh_ct, &hh_ct, indiv_f_nonmale_ct, &ll_ctf, &lh_ctf, &hh_ctf, hwe_needed, indiv_f_ctl_nonmale_ct, &ll_hwe, &lh_hwe, &hh_hwe, hardy_needed, indiv_f_case_nonmale_ct, &ll_case_hwe, &lh_case_hwe, &hh_case_hwe);
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
	    nonmissing_nonmale_y = nonmissing_present_diff(unfiltered_indiv_ctv2, loadbuf, indiv_include2, indiv_male_include2);
	  }
	  haploid_single_marker_freqs(unfiltered_indiv_ct, unfiltered_indiv_ctv2, loadbuf, indiv_male_include2, founder_male_include2, indiv_male_ct, &ll_ct, &hh_ct, indiv_f_male_ct, &ll_ctf, &hh_ctf, &hethap_incr);
	  if ((is_x || indiv_male_ct) && (indiv_ct - cur_oblig_missing)) {
	    if (is_x) {
	      if (!cur_oblig_missing) {
		cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct + ukk)) * indiv_ct_recip;
	      } else {
		cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct + ukk)) / ((double)((int32_t)(indiv_ct - cur_oblig_missing)));
	      }
	    } else {
	      if (!cur_oblig_missing) {
		cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct)) * male_ct_recip;
	      } else {
		cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct)) / ((double)((int32_t)(indiv_male_ct - cur_oblig_missing)));
	      }
	    }
	  } else {
	    cur_genotyping_rate = 0;
	    nonmissing_rate_tot_max -= 1;
	  }
	} else {
	  haploid_single_marker_freqs(unfiltered_indiv_ct, unfiltered_indiv_ctv2, loadbuf, indiv_include2, founder_include2, indiv_ct, &ll_ct, &hh_ct, indiv_f_ct, &ll_ctf, &hh_ctf, &hethap_incr);
	  if (!cur_oblig_missing) {
	    cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct)) * indiv_ct_recip;
	  } else {
	    if (indiv_ct - cur_oblig_missing) {
	      cur_genotyping_rate = ((int32_t)(ll_ct + hh_ct)) / ((double)((int32_t)(indiv_ct - cur_oblig_missing)));
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
	    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctv2; indiv_uidx++) {
	      ulii = loadbuf[indiv_uidx];
	      ulii = (ulii >> 1) & (~ulii) & indiv_male_include2[indiv_uidx];
	      while (ulii) {
		ukk = indiv_uidx * BITCT2 + CTZLU(ulii) / 2;
		fputs(&(person_ids[ukk * max_person_id_len]), hhfile);
		putc('\t', hhfile);
		fputs(&(marker_ids[marker_uidx * max_marker_id_len]), hhfile);
		putc('\n', hhfile);
		ulii &= ulii - ONELU;
	      }
	    }
	  } else {
	    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctv2; indiv_uidx++) {
	      ulii = loadbuf[indiv_uidx];
	      ulii = (ulii >> 1) & (~ulii) & indiv_include2[indiv_uidx];
	      while (ulii) {
		ukk = indiv_uidx * BITCT2 + CTZLU(ulii) / 2;
		fputs(&(person_ids[ukk * max_person_id_len]), hhfile);
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
	uii += ll_ctf + hh_ctf + 2 * maf_succ;
	ujj += hh_ctf + maf_succ;
	if (!uii) {
	  maf = 0.5;
	} else {
	  maf = ((double)ujj) / ((double)uii);
	}
	set_allele_freqs[marker_uidx] = maf;
	if (wt_needed) {
	  marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, maf);
	}
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
    LOGPRINTF("Warning: %" PRIu64 " het. haploid genotype%s present (see %s.hh).\n", hethap_ct, (hethap_ct == 1LLU)? "" : "s", outname);
  }
  if (nonmissing_nonmale_y) {
    logprint("Warning: Nonmissing nonmale Y chromosome genotype(s) present.\n");
    *hh_exists_ptr |= Y_FIX_NEEDED;
  }
  if (nonmissing_rate_tot <= 0.9999995 * ((double)((intptr_t)nonmissing_rate_tot_max))) {
    LOGPRINTF("Total genotyping rate %sis %g.\n", indiv_exclude_ct? "in remaining individuals " : "", nonmissing_rate_tot / ((double)((intptr_t)nonmissing_rate_tot_max)));
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

int32_t write_missingness_reports(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uintptr_t* pheno_nm, uintptr_t* sex_male, uint32_t indiv_male_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t hh_exists) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ct2l = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t unfiltered_indiv_ctv2 = (unfiltered_indiv_ct2l + 1) & (~1);
  uintptr_t marker_ct_y = 0;
  uintptr_t* indiv_male_include2 = NULL;
  uint64_t* om_entry_ptr = NULL;
  uintptr_t* cur_omidxs = NULL;
  uint32_t* indiv_to_cluster = NULL;
  uint32_t* missing_ct_by_cluster = NULL;
  uint32_t* oblig_missing_ct_by_cluster = NULL;
  uint32_t* cluster_sizes = NULL;
  uint32_t* cluster_sizes_y = NULL;
  uint32_t* om_cluster_sizes = NULL;
  uint32_t* om_indiv_lookup = NULL;
  uint32_t* om_cluster_ref_cts = NULL;
  uint64_t cur_om_entry = 0;
  int32_t y_code = chrom_info_ptr->y_code;
  uint32_t y_present = (y_code != -1) && is_set(chrom_info_ptr->chrom_mask, y_code);
  uint32_t indiv_uidx = 0;
  uint32_t indiv_idx = 0;
  uint32_t oblig_ct = 0;
  uint32_t om_cluster_ct = 0;
  uint32_t om_cluster_ctl = 0;
  int32_t retval = 0;
  uintptr_t* loadbuf;
  uintptr_t* indiv_include2;
  uintptr_t* cur_nm;
  uintptr_t* lptr;
  uintptr_t* lptr2;
  uint32_t* missing_cts;
  uint32_t* cur_cluster_sizes;
  char* wptr;
  char* cptr;
  char* cptr2;
  uintptr_t marker_ct_nony;
  uintptr_t marker_uidx;
  uintptr_t ulii;
  uint32_t slen;
  uint32_t indiv_uidx_stop;
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
  if (wkspace_alloc_ui_checked(&missing_cts, unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&indiv_male_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
    goto write_missingness_reports_ret_NOMEM;
  }
  loadbuf[unfiltered_indiv_ctv2 - 2] = 0;
  loadbuf[unfiltered_indiv_ctv2 - 1] = 0;
  exclude_to_vec_include(unfiltered_indiv_ct, indiv_include2, indiv_exclude);
  memcpy(indiv_male_include2, indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
  vec_include_mask_in(unfiltered_indiv_ct, indiv_male_include2, sex_male);
  if (y_present) {
    marker_ct_y = count_chrom_markers(chrom_info_ptr, chrom_info_ptr->y_code, marker_exclude);
  }
  marker_ct_nony = marker_ct - marker_ct_y;
  fill_uint_zero(missing_cts, unfiltered_indiv_ct);
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto write_missingness_reports_ret_READ_FAIL;
  }
  memcpy(outname_end, ".lmiss", 7);
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  if (om_ip->entry_ct) {
    om_entry_ptr = om_ip->entries;
    om_cluster_ref_cts = om_ip->cluster_ref_cts;
    cur_om_entry = *om_entry_ptr;
    om_cluster_ct = om_ip->cluster_ct;
    // divide by BITCT2 instead of BITCT due to Ychr
    om_cluster_ctl = (om_cluster_ct + BITCT - 1) / BITCT;
    om_indiv_lookup = om_ip->indiv_lookup;
    if (wkspace_alloc_ui_checked(&om_cluster_sizes, om_cluster_ct * 2 * sizeof(int32_t))) {
      goto write_missingness_reports_ret_NOMEM;
    }
    fill_uint_zero(om_cluster_sizes, om_cluster_ct * 2);
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_unsafe_ck(indiv_exclude, &indiv_uidx);
      uii = om_indiv_lookup[indiv_uidx];
      if (uii != 0xffffffffU) {
        om_cluster_sizes[uii] += 1;
        if (is_set(sex_male, indiv_uidx)) {
	  om_cluster_sizes[uii + om_cluster_ct] += 1;
	}
      }
    }
    indiv_uidx = 0;
    indiv_idx = 0;
    if (cluster_ct) {
      if (wkspace_alloc_ul_checked(&cur_omidxs, om_cluster_ctl * sizeof(intptr_t))) {
        goto write_missingness_reports_ret_NOMEM;
      }
    }
  }
  ujj = unfiltered_indiv_ct2l * BITCT2;
  if (!cluster_ct) {
    sprintf(tbuf, " CHR %%%us   N_MISS   N_GENO   F_MISS\n", plink_maxsnp);
  } else {
    if (wkspace_alloc_ui_checked(&indiv_to_cluster, unfiltered_indiv_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&missing_ct_by_cluster, cluster_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&oblig_missing_ct_by_cluster, cluster_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_sizes, cluster_ct * 2 * sizeof(int32_t))) {
      goto write_missingness_reports_ret_NOMEM;
    }
    fill_uint_zero(indiv_to_cluster, unfiltered_indiv_ct);
    fill_uint_zero(cluster_sizes, cluster_ct * 2);
    fill_uint_zero(oblig_missing_ct_by_cluster, cluster_ct);
    cluster_sizes_y = &(cluster_sizes[cluster_ct]);
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      unn = cluster_starts[clidx + 1];
      ukk = clidx + 1;
      for (uii = cluster_starts[clidx]; uii < unn; uii++) {
	umm = cluster_map[uii];
	if (!IS_SET(indiv_exclude, umm)) {
          indiv_to_cluster[umm] = ukk;
	  cluster_sizes[clidx] += 1;
          if (IS_SET(sex_male, umm)) {
            cluster_sizes_y[clidx] += 1;
	  }
	}
      }
    }
    sprintf(tbuf, " CHR %%%us       CLST   N_MISS   N_CLST   N_GENO   F_MISS\n", plink_maxsnp);
  }
  fprintf(outfile, tbuf, "SNP");
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    is_x = (((int32_t)chrom_idx) == chrom_info_ptr->x_code);
    is_y = (((int32_t)chrom_idx) == chrom_info_ptr->y_code);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    if (!is_y) {
      cur_nm = indiv_include2;
      cur_tot = indiv_ct;
      cur_cluster_sizes = cluster_sizes;
      om_ycorr = 0;
    } else {
      cur_nm = indiv_male_include2;
      cur_tot = indiv_male_ct;
      cur_cluster_sizes = cluster_sizes_y;
      om_ycorr = om_cluster_ct;
    }
    cptr = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, chrom_idx, zero_extra_chroms));
    *cptr++ = ' ';
    if (marker_uidx < chrom_end) {
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto write_missingness_reports_ret_READ_FAIL;
      }
      do {
	if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  goto write_missingness_reports_ret_READ_FAIL;
	}
        if (is_haploid) {
          haploid_fix(hh_exists, indiv_include2, indiv_male_include2, unfiltered_indiv_ct, is_x, is_y, (unsigned char*)loadbuf);
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
	  wptr = uint32_writew8x(cptr2, ukk - oblig_ct, ' ');
          wptr = uint32_writew8x(wptr, cur_tot - oblig_ct, ' ');
	  wptr = double_g_writewx4x(wptr, ((double)((int32_t)(ukk - oblig_ct))) / ((double)((int32_t)(cur_tot - oblig_ct))), 8, '\n');
          if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
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
		ukk = indiv_to_cluster[ukk];
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
		umm = indiv_to_cluster[ukk];
		if (umm) {
		  if (is_set(cur_omidxs, om_indiv_lookup[ukk])) {
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
            wptr = fw_strcpy(10, &(cluster_ids[clidx * max_cluster_id_len]), cptr2);
	    *wptr++ = ' ';
	    uii = missing_ct_by_cluster[clidx];
            wptr = uint32_writew8x(wptr, uii, ' ');
	    umm = cur_cluster_sizes[clidx];
	    wptr = uint32_writew8x(wptr, umm, ' ');
	    umm -= oblig_missing_ct_by_cluster[clidx];
	    wptr = uint32_writew8x(wptr, umm, ' ');
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)uii)) / ((double)((int32_t)umm)), 8, '\n');
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto write_missingness_reports_ret_WRITE_FAIL;
	    }
	  }
	}
        marker_uidx++;
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	  if (marker_uidx < chrom_end) {
	    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	      goto write_missingness_reports_ret_WRITE_FAIL;
	    }
	  }
	}
      } while (marker_uidx < chrom_end);
    }
  }
  if (fclose_null(&outfile)) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  outname_end[1] = 'i';
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  sprintf(tbuf, "%%%us %%%us MISS_PHENO   N_MISS   N_GENO   F_MISS\n", plink_maxfid, plink_maxiid);
  fprintf(outfile, tbuf, "FID", "IID");
  do {
    indiv_uidx = next_unset_unsafe(indiv_exclude, indiv_uidx);
    indiv_uidx_stop = next_set(indiv_exclude, indiv_uidx, unfiltered_indiv_ct);
    indiv_idx += indiv_uidx_stop - indiv_uidx;
    do {
      cptr = &(person_ids[indiv_uidx * max_person_id_len]);
      cptr2 = (char*)memchr(cptr, '\t', max_person_id_len);
      slen = (uintptr_t)(cptr2 - cptr);
      wptr = memseta(tbuf, 32, plink_maxfid - slen);
      wptr = memcpyax(wptr, cptr, slen, ' ');
      wptr = fw_strcpy(plink_maxiid, &(cptr2[1]), wptr);
      wptr = memseta(wptr, 32, 10);
      *wptr++ = 'Y' - (is_set(pheno_nm, indiv_uidx) * 11);
      *wptr++ = ' ';
      uii = missing_cts[indiv_uidx];
      ukk = is_set(sex_male, indiv_uidx);
      ujj = marker_ct_nony + (ukk * marker_ct_y);
      if (om_indiv_lookup) {
	umm = om_indiv_lookup[indiv_uidx];
	if (umm != 0xffffffffU) {
          umm = om_cluster_ref_cts[umm] + ukk * om_cluster_ref_cts[umm + om_cluster_ct];
	  uii -= umm;
	  ujj -= umm;
	}
      }
      wptr = uint32_writew8x(wptr, uii, ' ');
      wptr = uint32_writew8x(wptr, ujj, ' ');
      wptr = double_g_writewx4x(wptr, ((double)((int32_t)uii)) / ((double)((int32_t)ujj)), 8, '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto write_missingness_reports_ret_WRITE_FAIL;
      }
    } while (++indiv_uidx < indiv_uidx_stop);
  } while (indiv_idx < indiv_ct);
  if (fclose_null(&outfile)) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  *outname_end = '\0';
  LOGPRINTF("--missing: Individual missing data report written to %s.imiss, and\nvariant-based %smissing data report written to %s.lmiss.\n", outname, cluster_ct? "cluster-stratified " : "", outname);
  while (0) {
  write_missingness_reports_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_missingness_reports_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  write_missingness_reports_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t hardy_report_write_line(FILE* outfile, char* prefix_buf, uint32_t prefix_len, uint32_t reverse, uint32_t ll_ct, uint32_t lh_ct, uint32_t hh_ct, char* midbuf_ptr, double pvalue) {
  char wbuf[48];
  char* cptr;
  uint32_t denom;
  double drecip;
  double minor_freq;
  fwrite(prefix_buf, 1, prefix_len, outfile);
  if (reverse) {
    cptr = uint32_write(uint32_writex(uint32_writex(wbuf, hh_ct, '/'), lh_ct, '/'), ll_ct);
  } else {
    cptr = uint32_write(uint32_writex(uint32_writex(wbuf, ll_ct, '/'), lh_ct, '/'), hh_ct);
  }
  cptr = fw_strcpyn(20, cptr - wbuf, wbuf, midbuf_ptr);
  *cptr++ = ' ';
  denom = (ll_ct + lh_ct + hh_ct) * 2;
  if (denom) {
    drecip = 1.0 / ((double)denom);
    minor_freq = (2 * ll_ct + lh_ct) * drecip;
    cptr = double_g_writewx4x(double_g_writewx4x(double_g_writewx4x(cptr, (lh_ct * 2) * drecip, 8, ' '), minor_freq * (2 * hh_ct + lh_ct) * drecip * 2, 8, ' '), pvalue, 12, '\n');
  } else {
    cptr = memcpya(cptr, "     nan      nan           NA\n", 31);
  }
  return fwrite_checked(midbuf_ptr, (cptr - midbuf_ptr), outfile);
}

int32_t hardy_report(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, uint32_t hwe_modifier, int32_t* hwe_ll_cases, int32_t* hwe_lh_cases, int32_t* hwe_hh_cases, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs, uint32_t pheno_nm_ct, uintptr_t* pheno_c, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr) {
  FILE* outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  uint32_t hwe_midp = hwe_modifier & HWE_MIDP;
  int32_t retval = 0;
  uint32_t skip_chrom = 0;
  uint32_t pct = 0;
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
  char* writebuf;
  char* cptr0;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* cptr5;
  if (pheno_nm_ct) {
    report_type = pheno_c? 0 : 1;
  } else {
    report_type = 2;
  }
  uii = report_type? 1 : 3;
  if (wkspace_alloc_d_checked(&p_values, uii * marker_ct * sizeof(double)) ||
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

  memcpy(outname_end, ".hwe", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto hardy_report_ret_OPEN_FAIL;
  }
  LOGPRINTF("Writing Hardy-Weinberg report to %s...", outname);
  fputs(" 0%", stdout);
  sprintf(writebuf, " CHR %%%us     TEST   A1   A2                 GENO   O(HET)   E(HET)            P \n", plink_maxsnp);
  fprintf(outfile, writebuf, "SNP");
 
  chrom_fo_idx = 0;
  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
  skip_chrom = (is_haploid && (!is_x)) || is_mt;
  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
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
	  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
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
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[marker_idx])) {
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
	  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
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
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[3 * marker_idx])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(cptr0[7 + plink_maxsnp]), "FF", 2);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_ll_cases[marker_uidx], hwe_lh_cases[marker_uidx], hwe_hh_cases[marker_uidx], cptr2, p_values[3 * marker_idx + 1])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(cptr0[4 + plink_maxsnp]), "UN", 2);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_lls[marker_uidx], hwe_lhs[marker_uidx], hwe_hhs[marker_uidx], cptr2, p_values[3 * marker_idx + 2])) {
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
  fputs("\b\b\b\b", stdout);
  logprint(" done.\n");

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
  fclose_cond(outfile);
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
  LOGPRINTF("%u variant%s removed due to Hardy-Weinberg exact test (--hwe).\n", removed_ct, (removed_ct == 1)? "" : "s");
  *marker_exclude_ct_ptr += removed_ct;
  return 0;
}

uint32_t enforce_maf_threshold(double min_maf, double max_maf, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, double* set_allele_freqs) {
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t removed_ct = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  double dxx;
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
        removed_ct++;
      }
    } while (++marker_uidx < marker_uidx_stop);
  }
  if (marker_ct == removed_ct) {
    logprint("Error: All variants removed due to MAF threshold(s) (--maf/--max-maf).\n");
    return 1;
  }
  LOGPRINTF("%u variant%s removed due to MAF threshold(s) (--maf/--max-maf).\n", removed_ct, (removed_ct == 1)? "" : "s");
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

void fill_mendel_errstr(uint32_t error_code, char** allele_ptrs, uint32_t* alens, char* wbuf, uint32_t* len_ptr) {
  char* wptr;
  uint32_t len;
  if (error_code < 10) {
    wbuf[0] = ' ';
    wbuf[1] = error_code + '0';
  } else {
    memcpy(wbuf, "10", 2);
  }
  if (!alens[0]) {
    // lazy fill
    alens[0] = strlen(allele_ptrs[0]);
    alens[1] = strlen(allele_ptrs[1]);
  }
  // PLINK 1.07 may fail to put space here when there are long allele codes; we
  // don't replicate that.  (We actually force two spaces here since the error
  // column itself contains spaces.)
  wptr = memseta(&(wbuf[2]), 32, 2);
  switch (error_code) {
  case 1:
    len = 5 * alens[0] + alens[1] + 10;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpyl3a(wptr, " x ");
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  case 2:
    len = alens[0] + 5 * alens[1] + 10;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpyl3a(wptr, " x ");
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  case 3:
    len = 2 * alens[0] + 2 * alens[1] + 12;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpya(wptr, " x */* -> ", 10);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    break;
  case 4:
  case 10:
    len = 2 * alens[0] + 2 * alens[1] + 12;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpya(wptr, "*/* x ", 6);
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    break;
  case 5:
    len = 2 * alens[0] + 4 * alens[1] + 10;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpyl3a(wptr, " x ");
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    break;
  case 6:
    len = 2 * alens[0] + 2 * alens[1] + 12;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpya(wptr, " x */* -> ", 10);
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  case 7:
  case 9:
    len = 2 * alens[0] + 2 * alens[1] + 12;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpya(wptr, "*/* x ", 6);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  case 8:
    len = 4 * alens[0] + 2 * alens[1] + 10;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpyl3a(wptr, " x ");
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  }
  *wptr++ = '\n';
  *len_ptr = (uintptr_t)(wptr - wbuf);
}

int32_t mendel_error_scan(Mendel_info* me_ip, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t hh_exists, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uint32_t calc_mendel) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  FILE* outfile_l = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uintptr_t* error_locs = NULL;
  char* varptr = NULL;
  char* chrom_name_ptr = NULL;
  unsigned char* cur_errors = NULL;
  uint64_t* family_error_cts = NULL;
  uint32_t* child_cts = NULL;
  uintptr_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctp1l2 = 1 + (unfiltered_indiv_ct / BITCT2);
  uintptr_t indiv_ct = unfiltered_indiv_ct - *indiv_exclude_ct_ptr;
  uintptr_t marker_uidx = 0xffffffffLU;
  uint64_t tot_error_ct = 0;
  uint32_t do_filter = me_ip->modifier & MENDEL_FILTER;
  uint32_t include_duos = (me_ip->modifier / MENDEL_DUOS) & 1;
  uint32_t multigen = (me_ip->modifier / MENDEL_MULTIGEN) & 1;
  uint32_t var_first = me_ip->modifier & MENDEL_FILTER_VAR_FIRST;
  uint32_t varlen = 0;
  uint32_t chrom_name_len = 0;
  uint32_t new_marker_exclude_ct = 0;
  uint32_t error_ct_fill = 0;
  int32_t retval = 0;

  // bottom 2 bits of index = child genotype
  // middle 2 bits of index = paternal genotype
  // top 2 bits of index = maternal genotype

  // bits 0-7 = child increment (always 1)
  // bits 8-15 = father increment
  // bits 16-23 = mother increment
  // bits 24-31 = error code

  // note that \xx is octal, not decimal.
  const uint32_t error_table[] =
{0, 0, 0x1010101, 0x8000001,
 0, 0, 0, 0x7010001,
 0, 0, 0, 0x7010001,
 0x3000101, 0, 0, 0x7010001,
 0, 0, 0, 0x6000101,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0x3000101, 0, 0, 0,
 0, 0, 0, 0x6000101,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0x3000101, 0, 0, 0,
 0x4010001, 0, 0, 0x6000101,
 0x4010001, 0, 0, 0,
 0x4010001, 0, 0, 0,
 0x5000001, 0, 0x2010101, 0};
  // necessary to check child gender when dealing with error 9/10
  const uint32_t error_table_x[] =
{0, 0, 0x1010101, 0x8000001,
 0, 0, 0, 0x9010001,
 0, 0, 0, 0x7010001,
 0x3000101, 0, 0, 0x7010001,
 0, 0, 0, 0x6000101,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0x3000101, 0, 0, 0,
 0, 0, 0, 0x6000101,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0x3000101, 0, 0, 0,
 0x4010001, 0, 0, 0x6000101,
 0xa010001, 0, 0, 0,
 0x4010001, 0, 0, 0,
 0x5000001, 0, 0x2010101, 0};
  char chrom_name_buf[4];
  char* errstrs[10];
  uint32_t errstr_lens[11];
  uint32_t alens[2];
  const uint32_t* error_table_ptr;
  char* fids;
  char* iids;
  char* wptr;
  char* cptr;
  uint64_t* family_list;
  uint64_t* trio_list;
  uintptr_t* loadbuf;
  uint32_t* trio_lookup;
  uint32_t* error_cts;
  uint32_t* error_cts_tmp;
  uint32_t* error_cts_tmp2;
  uint32_t* uiptr;
#ifdef __LP64__
  __m128i* vptr;
  __m128i* vptr2;
#endif
  uintptr_t trio_ct4;
  uintptr_t max_fid_len;
  uintptr_t max_iid_len;
  uintptr_t trio_ct;
  uintptr_t trio_ctl;
  uintptr_t trio_idx;
  uintptr_t lookup_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  double exclude_one_ratio;
  double dxx;
  double dyy;
  uint64_t trio_code;
  uint64_t family_code;
  uint32_t family_ct;
  uint32_t var_error_max;
  uint32_t cur_error_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t is_x;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  marker_ct -= count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 0, 1);
  if ((!marker_ct) || IS_SET(chrom_info_ptr->haploid_mask, 0)) {
    logprint("Warning: Skipping --me/--mendel since there is no autosomal or Xchr data.\n");
    goto mendel_error_scan_ret_1;
  }
  retval = get_trios_and_families(unfiltered_indiv_ct, indiv_exclude, indiv_ct, founder_info, sex_nm, sex_male, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, &fids, &max_fid_len, &iids, &max_iid_len, &family_list, &family_ct, &trio_list, &trio_ct, &trio_lookup, include_duos, multigen);
  if (retval) {
    goto mendel_error_scan_ret_1;
  }
  if (!trio_ct) {
    LOGPRINTF("Warning: Skipping --me/--mendel since there are no %strios.\n", include_duos? "duos or " : "");
    goto mendel_error_scan_ret_1;
  }
  trio_ct4 = (trio_ct + 3) / 4;
  trio_ctl = (trio_ct + (BITCT - 1)) / BITCT;
  var_error_max = (int32_t)(me_ip->max_var_error * (1 + SMALL_EPSILON) * ((intptr_t)trio_ct));
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctp1l2 * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&error_cts, trio_ct * 3 * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&error_cts_tmp, trio_ct4 * 4 * sizeof(int32_t))) {
    goto mendel_error_scan_ret_NOMEM;
  }
  if (!var_first) {
    error_cts_tmp2 = error_cts_tmp;
  } else {
    if (wkspace_alloc_ui_checked(&error_cts_tmp2, trio_ct4 * 4 * sizeof(int32_t))) {
      goto mendel_error_scan_ret_NOMEM;
    }
    fill_uint_zero(error_cts_tmp2, trio_ct4 * 4);
  }
  loadbuf[unfiltered_indiv_ctp1l2 - 1] = 0;
  fill_uint_zero(error_cts, trio_ct * 3);
  fill_uint_zero(error_cts_tmp, trio_ct4 * 4);
  hh_exists &= XMHH_EXISTS;
  if (alloc_raw_haploid_filters(unfiltered_indiv_ct, hh_exists, 0, indiv_exclude, sex_male, NULL, &indiv_male_include2)) {
    goto mendel_error_scan_ret_NOMEM;
  }
  alens[0] = 0;
  alens[1] = 0;
  if (calc_mendel) {
    ulii = max_marker_allele_len * 6 + 10;
    if (ulii < 21) {
      ulii = 21;
    }
    if (wkspace_alloc_ull_checked(&family_error_cts, family_ct * 3 * sizeof(int64_t)) ||
        wkspace_alloc_ui_checked(&child_cts, family_ct * sizeof(int32_t)) ||
        wkspace_alloc_c_checked(&(errstrs[0]), ulii * 10)) {
      goto mendel_error_scan_ret_NOMEM;
    }
    for (uii = 1; uii < 10; uii++) {
      errstrs[uii] = &(errstrs[0][uii * ulii]);
    }
    if (multigen) {
      if (wkspace_alloc_ul_checked(&error_locs, trio_ctl * sizeof(intptr_t)) ||
          wkspace_alloc_uc_checked(&cur_errors, trio_ct)) {
	goto mendel_error_scan_ret_NOMEM;
      }
      fill_ulong_zero(error_locs, trio_ctl);
    }
    memcpy(outname_end, ".mendel", 8);
    if (fopen_checked(&outfile, outname, "w")) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    sprintf(tbuf, "%%%us %%%us  CHR %%%us   CODE                 ERROR\n", plink_maxfid, plink_maxiid, plink_maxsnp);
    fprintf(outfile, tbuf, "FID", "KID", "SNP");
    memcpy(outname_end, ".lmendel", 9);
    if (fopen_checked(&outfile_l, outname, "w")) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    // replicate harmless 'N' misalignment bug
    sprintf(tbuf, " CHR %%%us   N\n", plink_maxsnp);
    fprintf(outfile_l, tbuf, "SNP");
  } else {
    // suppress warning
    fill_ulong_zero((uintptr_t*)errstrs, 10);
  }
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    is_x = (((uint32_t)chrom_info_ptr->x_code) == chrom_idx)? 1 : 0;
    if ((IS_SET(chrom_info_ptr->haploid_mask, chrom_idx) && (!is_x)) || (((uint32_t)chrom_info_ptr->mt_code) == chrom_idx)) {
      continue;
    }
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    uii = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    if (uii == chrom_end) {
      continue;
    }
    if (!is_x) {
      error_table_ptr = error_table;
    } else {
      error_table_ptr = error_table_x;
    }
    if (calc_mendel) {
      chrom_name_ptr = chrom_name_buf;
      chrom_name_len = 4;
      if (chrom_idx <= chrom_info_ptr->max_code) {
	chrom_num_write4(chrom_name_buf, chrom_idx);
      } else if (zero_extra_chroms) {
	memcpy(chrom_name_buf, "   0", 4);
      } else {
	ujj = strlen(chrom_info_ptr->nonstd_names[chrom_idx]);
	if (ujj < 4) {
	  fw_strcpyn(4, ujj, chrom_info_ptr->nonstd_names[chrom_idx], chrom_name_buf);
	} else {
	  chrom_name_ptr = chrom_info_ptr->nonstd_names[chrom_idx];
	  chrom_name_len = ujj;
	}
      }
    }
    if (uii != marker_uidx) {
      marker_uidx = uii;
      goto mendel_error_scan_seek;
    }
    while (1) {
      if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto mendel_error_scan_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
        reverse_loadbuf((unsigned char*)loadbuf, unfiltered_indiv_ct);
      }
      if (is_x && hh_exists) {
	hh_reset((unsigned char*)loadbuf, indiv_male_include2, unfiltered_indiv_ct);
      }
      // missing parents are treated as having uidx equal to
      // unfiltered_indiv_ct, and we set the corresponding genotype to always
      // be missing.  This lets us avoid special-casing duos.
      SET_BIT_DBL(loadbuf, unfiltered_indiv_ct);
      uiptr = trio_lookup;
      cur_error_ct = 0;
      if (calc_mendel) {
	varptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	varlen = strlen(varptr);
	alens[0] = 0;
	alens[1] = 0;
	fill_uint_zero(errstr_lens, 11);
      }
      if (!multigen) {
	for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
	  uii = *uiptr++;
	  ujj = *uiptr++;
	  ukk = *uiptr++;
          umm = error_table_ptr[((loadbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & 3) | (((loadbuf[ujj / BITCT2] >> (2 * (ujj % BITCT2))) & 3) << 2) | (((loadbuf[ukk / BITCT2] >> (2 * (ukk % BITCT2))) & 3) << 4)];
	  if (umm) {
	    error_cts_tmp2[trio_idx] += umm & 0xffffff;
	    cur_error_ct++;
	    if (calc_mendel) {
	      umm >>= 24;
	      if ((umm > 8) && (!is_set(sex_male, uii))) {
		umm = 34 - 3 * umm; // 9 -> 7, 10 -> 4
	      }
	      wptr = fw_strcpy(plink_maxfid, &(fids[trio_idx * max_fid_len]), tbuf);
	      *wptr++ = ' ';
	      wptr = fw_strcpy(plink_maxiid, &(iids[uii * max_iid_len]), wptr);
	      *wptr++ = ' ';
	      wptr = memcpyax(wptr, chrom_name_ptr, chrom_name_len, ' ');
	      wptr = fw_strcpyn(plink_maxsnp, varlen, varptr, wptr);
	      wptr = memseta(wptr, 32, 5);
	      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
		goto mendel_error_scan_ret_WRITE_FAIL;
	      }
	      if (!errstr_lens[umm]) {
		fill_mendel_errstr(umm, &(marker_allele_ptrs[2 * marker_uidx]), alens, errstrs[umm - 1], &(errstr_lens[umm]));
	      }
	      if (fwrite_checked(errstrs[umm - 1], errstr_lens[umm], outfile)) {
		goto mendel_error_scan_ret_WRITE_FAIL;
	      }
	    }
	  }
	}
      } else {
	for (lookup_idx = 0; lookup_idx < trio_ct; lookup_idx++) {
	  uii = *uiptr++;
	  ujj = *uiptr++;
	  ukk = *uiptr++;
          trio_idx = *uiptr++;
          uljj = ((loadbuf[ujj / BITCT2] >> (2 * (ujj % BITCT2))) & 3) | (((loadbuf[ukk / BITCT2] >> (2 * (ukk % BITCT2))) & 3) << 2);
	  umm = uii / BITCT2;
	  ujj = 2 * (uii % BITCT2);
	  ulii = (loadbuf[umm] >> ujj) & 3;
	  if (ulii != 1) {
	    umm = error_table_ptr[ulii | (uljj << 2)];
	    if (umm) {
	      error_cts_tmp2[trio_idx] += umm & 0xffffff;
	      cur_error_ct++;
	      if (calc_mendel) {
	        set_bit(error_locs, trio_idx);
		umm >>= 24;
		if ((umm > 8) && (!is_set(sex_male, uii))) {
		  umm = 34 - 3 * umm;
		}
                cur_errors[trio_idx] = (unsigned char)umm;
	      }
	    }
	  } else if (uljj == 0) {
	    loadbuf[umm] &= ~(ONELU << ujj);
          } else if (uljj == 15) {
	    loadbuf[umm] |= (2 * ONELU) << ujj;
	  }
	}
	if (calc_mendel && cur_error_ct) {
          trio_idx = 0;
	  for (uii = 0; uii < cur_error_ct; trio_idx++, uii++) {
            next_set_ul_unsafe_ck(error_locs, &trio_idx);
	    wptr = fw_strcpy(plink_maxfid, &(fids[trio_idx * max_fid_len]), tbuf);
	    *wptr++ = ' ';
	    wptr = fw_strcpy(plink_maxiid, &(iids[((uint32_t)trio_list[trio_idx]) * max_iid_len]), wptr);
	    *wptr++ = ' ';
	    wptr = memcpyax(wptr, chrom_name_ptr, chrom_name_len, ' ');
	    wptr = fw_strcpyn(plink_maxsnp, varlen, varptr, wptr);
	    wptr = memseta(wptr, 32, 5);
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto mendel_error_scan_ret_WRITE_FAIL;
	    }
	    umm = cur_errors[trio_idx];
	    if (!errstr_lens[umm]) {
	      fill_mendel_errstr(umm, &(marker_allele_ptrs[2 * marker_uidx]), alens, errstrs[umm - 1], &(errstr_lens[umm]));
	    }
	    if (fwrite_checked(errstrs[umm - 1], errstr_lens[umm], outfile)) {
	      goto mendel_error_scan_ret_WRITE_FAIL;
	    }
	  }
          fill_ulong_zero(error_locs, trio_ctl);
	}
      }
      if (calc_mendel) {
	if (fwrite_checked(chrom_name_ptr, chrom_name_len, outfile_l)) {
	  goto mendel_error_scan_ret_WRITE_FAIL;
	}
	tbuf[0] = ' ';
	wptr = fw_strcpyn(plink_maxsnp, varlen, varptr, &(tbuf[1]));
        *wptr++ = ' ';
        wptr = uint32_writew4x(wptr, cur_error_ct, '\n');
	if (fwrite_checked(tbuf, wptr - tbuf, outfile_l)) {
	  goto mendel_error_scan_ret_WRITE_FAIL;
	}
      }
      if (cur_error_ct) {
	if (cur_error_ct > var_error_max) {
	  SET_BIT(marker_exclude, marker_uidx);
	  new_marker_exclude_ct++;
	}
	if ((cur_error_ct <= var_error_max) || (!var_first)) {
	  if (var_first) {
#ifdef __LP64__
	    vptr = (__m128i*)error_cts_tmp;
	    vptr2 = (__m128i*)error_cts_tmp2;
	    for (trio_idx = 0; trio_idx < trio_ct4; trio_idx++) {
	      *vptr = _mm_add_epi64(*vptr, *vptr2++);
	      vptr++;
	    }
#else
            for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
	      error_cts_tmp[trio_idx] += error_cts_tmp2[trio_idx];
	    }
#endif
	    fill_uint_zero(error_cts_tmp2, trio_ct);
	  }
	  error_ct_fill++;
	  if (error_ct_fill == 255) {
	    uiptr = error_cts;
            for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
	      uii = error_cts_tmp[trio_idx];
	      *uiptr += (unsigned char)uii;
	      uiptr++;
	      *uiptr += (unsigned char)(uii >> 8);
	      uiptr++;
	      *uiptr += uii >> 16;
	      uiptr++;
	    }
	    fill_uint_zero(error_cts_tmp, trio_ct);
	    error_ct_fill = 0;
	  }
	}
      }
      tot_error_ct += cur_error_ct;
      if (++marker_uidx == chrom_end) {
	break;
      }
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	if (marker_uidx == chrom_end) {
	  break;
	}
      mendel_error_scan_seek:
        if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto mendel_error_scan_ret_READ_FAIL;
	}
      }
    }
  }
  if (error_ct_fill) {
    uiptr = error_cts;
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      uii = error_cts_tmp[trio_idx];
      *uiptr += (unsigned char)uii;
      uiptr++;
      *uiptr += (unsigned char)(uii >> 8);
      uiptr++;
      *uiptr += uii >> 16;
      uiptr++;
    }
  }
  LOGPRINTF("--me/--mendel: %" PRIu64 " Mendel error%s detected.\n", tot_error_ct, (tot_error_ct == 1)? "" : "s");
  if (calc_mendel) {
    if (fclose_null(&outfile)) {
      goto mendel_error_scan_ret_WRITE_FAIL;
    }
    if (fclose_null(&outfile_l)) {
      goto mendel_error_scan_ret_WRITE_FAIL;
    }
    outname_end[1] = 'f';
    if (fopen_checked(&outfile, outname, "w")) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    sprintf(tbuf, "%%%us %%%us %%%us   CHLD    N\n", plink_maxfid, plink_maxiid, plink_maxiid);
    fprintf(outfile, tbuf, "FID", "PAT", "MAT");
    fill_ull_zero(family_error_cts, family_ct * 3);
    fill_uint_zero(child_cts, family_ct);
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      uii = (uint32_t)(trio_list[trio_idx] >> 32);
      child_cts[uii] += 1;
      family_error_cts[uii * 3] += error_cts[trio_idx * 3];
      family_error_cts[uii * 3 + 1] += error_cts[trio_idx * 3 + 1];
      family_error_cts[uii * 3 + 2] += error_cts[trio_idx * 3 + 2];
    }
    for (uii = 0; uii < family_ct; uii++) {
      family_code = family_list[uii];
      ujj = (uint32_t)family_code; // paternal uidx
      if (ujj < unfiltered_indiv_ct) {
	// bleah, fids[] isn't in right order for this lookup
	cptr = &(person_ids[ujj * max_person_id_len]);
	wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(((char*)memchr(cptr, '\t', max_person_id_len)) - cptr), cptr, tbuf);
      } else {
	wptr = memseta(tbuf, 32, plink_maxfid - 1);
	*wptr++ = '0';
      }
      *wptr++ = ' ';
      wptr = fw_strcpy(plink_maxiid, &(iids[ujj * max_iid_len]), wptr);
      *wptr++ = ' ';
      wptr = fw_strcpy(plink_maxiid, &(iids[((uintptr_t)(family_code >> 32)) * max_iid_len]), wptr);
      *wptr++ = ' ';
      wptr = uint32_writew6x(wptr, child_cts[uii], ' ');
      wptr = int64_write(wptr, family_error_cts[uii * 3]);
      *wptr++ = '\n';
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto mendel_error_scan_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto mendel_error_scan_ret_WRITE_FAIL;
    }
    outname_end[1] = 'i';
    if (fopen_checked(&outfile, outname, "w")) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    sprintf(tbuf, "%%%us %%%us   N\n", plink_maxfid, plink_maxiid);
    fprintf(outfile, tbuf, "FID", "IID");
    uii = 0xffffffffU; // family idx
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      trio_code = trio_list[trio_idx];
      ujj = (uint32_t)(trio_code >> 32);
      if (ujj != uii) {
	uii = ujj;
        family_code = family_list[uii];
	wptr = fw_strcpy(plink_maxfid, &(fids[trio_idx * max_fid_len]), tbuf);
	*wptr++ = ' ';
	ujj = (uint32_t)family_code;
	if (ujj != unfiltered_indiv_ct) {
	  wptr = fw_strcpy(plink_maxiid, &(iids[ujj * max_iid_len]), wptr);
	  *wptr++ = ' ';
	  wptr = int64_write(wptr, family_error_cts[3 * uii + 1]);
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto mendel_error_scan_ret_WRITE_FAIL;
	  }
	}
	ukk = (uint32_t)(family_code >> 32);
	if (ukk != unfiltered_indiv_ct) {
	  if (ujj != unfiltered_indiv_ct) {
	    putc('\n', outfile);
	  }
	  wptr = fw_strcpy(plink_maxiid, &(iids[ukk * max_iid_len]), &(tbuf[plink_maxfid + 1]));
	  *wptr++ = ' ';
	  wptr = int64_write(wptr, family_error_cts[3 * uii + 2]);
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto mendel_error_scan_ret_WRITE_FAIL;
	  }
	}
	putc(' ', outfile); // PLINK 1.07 formatting quirk
	putc('\n', outfile);
      }
      wptr = fw_strcpy(plink_maxiid, &(iids[((uint32_t)trio_code) * max_iid_len]), &(tbuf[plink_maxfid + 1]));
      *wptr++ = ' ';
      wptr = uint32_writew4x(wptr, error_cts[trio_idx * 3], '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto mendel_error_scan_ret_WRITE_FAIL;
      }
    }
    *outname_end = '\0';
    LOGPRINTF("Report written to %s.{mendel,imendel,fmendel,lmendel}.\n", outname);
  }
  if (do_filter) {
    *marker_exclude_ct_ptr += new_marker_exclude_ct;
    if (unfiltered_marker_ct == *marker_exclude_ct_ptr) {
      logprint("Error: All variants excluded by --me.\n");
      goto mendel_error_scan_ret_ALL_MARKERS_EXCLUDED;
    }
    if (var_first) {
      marker_ct -= new_marker_exclude_ct;
    }
    uii = (int32_t)(me_ip->max_trio_error * (1 + SMALL_EPSILON) * ((intptr_t)marker_ct));
    if (uii < marker_ct) {
      exclude_one_ratio = me_ip->exclude_one_ratio;
      for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
	if (error_cts[trio_idx * 3] > uii) {
	  trio_code = trio_list[trio_idx];
	  family_code = family_list[(uintptr_t)(trio_code >> 32)];
	  ujj = (uint32_t)family_code;
	  ukk = (uint32_t)(family_code >> 32);
          if (exclude_one_ratio == 0.0) {
	    set_bit(indiv_exclude, (uint32_t)trio_code);
	    if (ujj < unfiltered_indiv_ct) {
	      set_bit(indiv_exclude, ujj);
	    }
	    if (ukk < unfiltered_indiv_ct) {
	      set_bit(indiv_exclude, ukk);
	    }
	  } else if ((exclude_one_ratio == INFINITY) || (ujj == unfiltered_indiv_ct) || (ukk == unfiltered_indiv_ct)) {
            set_bit(indiv_exclude, (uint32_t)trio_code);
	  } else {
	    dxx = (double)((int32_t)trio_list[trio_idx * 3 + 1]);
	    dyy = (double)((int32_t)trio_list[trio_idx * 3 + 2]);
	    if (dxx > exclude_one_ratio * dyy) {
	      set_bit(indiv_exclude, ujj);
	    } else if (dyy > exclude_one_ratio * dxx) {
	      set_bit(indiv_exclude, ukk);
	    } else {
	      set_bit(indiv_exclude, (uint32_t)trio_code);
	    }
	  }
	}
      }
    }
    ulii = popcount_longs(indiv_exclude, (unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
    if (unfiltered_indiv_ct == ulii) {
      LOGPRINTF("Error: All %s excluded by --me.\n", g_species_plural);
      goto mendel_error_scan_ret_ALL_SAMPLES_EXCLUDED;
    }
    *indiv_exclude_ct_ptr = ulii;
    ulii -= unfiltered_indiv_ct - indiv_ct;
    LOGPRINTF("%u variant%s and %" PRIuPTR " %s excluded.\n", new_marker_exclude_ct, (new_marker_exclude_ct == 1)? "" : "s", ulii, species_str(ulii));
  }
  while (0) {
  mendel_error_scan_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  mendel_error_scan_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  mendel_error_scan_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  mendel_error_scan_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  mendel_error_scan_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  mendel_error_scan_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
 mendel_error_scan_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_l);
  return retval;
}
