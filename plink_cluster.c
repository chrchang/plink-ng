#include "plink_common.h"

#include "plink_cluster.h"
#include "plink_matrix.h"

void cluster_init(Cluster_info* cluster_ptr) {
  cluster_ptr->fname = NULL;
  cluster_ptr->match_fname = NULL;
  cluster_ptr->match_missing_str = NULL;
  cluster_ptr->match_type_fname = NULL;
  cluster_ptr->qmatch_fname = NULL;
  cluster_ptr->qmatch_missing_str = NULL;
  cluster_ptr->qt_fname = NULL;
  cluster_ptr->keep_fname = NULL;
  cluster_ptr->remove_fname = NULL;
  cluster_ptr->keep_flattened = NULL;
  cluster_ptr->remove_flattened = NULL;
  cluster_ptr->zerofname = NULL;
  cluster_ptr->modifier = 0;
  cluster_ptr->ppc = 0.0;
  cluster_ptr->max_size = 0xffffffffU;
  cluster_ptr->max_cases = 0xffffffffU;
  cluster_ptr->max_ctrls = 0xffffffffU;
  cluster_ptr->min_ct = 1;
  cluster_ptr->mds_dim_ct = 0;
  cluster_ptr->cmh_mperm_val = 0;
  cluster_ptr->min_ibm = 0.0;
}

void cluster_cleanup(Cluster_info* cluster_ptr) {
  free_cond(cluster_ptr->fname);
  free_cond(cluster_ptr->match_fname);
  free_cond(cluster_ptr->match_missing_str);
  free_cond(cluster_ptr->match_type_fname);
  free_cond(cluster_ptr->qmatch_fname);
  free_cond(cluster_ptr->qmatch_missing_str);
  free_cond(cluster_ptr->qt_fname);
  free_cond(cluster_ptr->keep_fname);
  free_cond(cluster_ptr->remove_fname);
  free_cond(cluster_ptr->keep_flattened);
  free_cond(cluster_ptr->remove_flattened);
  free_cond(cluster_ptr->zerofname);
}

int32_t load_clusters(char* fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, char* sample_ids, uintptr_t max_sample_id_len, uint32_t mwithin_col, uint32_t keep_na, uintptr_t* cluster_ct_ptr, uint32_t** cluster_map_ptr, uint32_t** cluster_starts_ptr, char** cluster_ids_ptr, uintptr_t* max_cluster_id_len_ptr, char* keep_fname, char* keep_flattened, char* remove_fname, char* remove_flattened) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  uintptr_t* sample_exclude_new = NULL;
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t sample_exclude_ct = *sample_exclude_ct_ptr;
  uintptr_t sample_ct = unfiltered_sample_ct - sample_exclude_ct;
  uintptr_t sample_ctl = (sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t topsize = 0;
  uintptr_t max_cluster_kr_len = 0;
  uint32_t cluster_filter = (keep_fname || keep_flattened || remove_fname || remove_flattened);
  uint32_t cluster_kr_ct = 0;
  int32_t retval = 0;
  char* idbuf = &(tbuf[MAXLINELEN]);
  // even if both --keep-clusters and --remove-clusters were specified, only
  // one is effectively active (i.e. any names in both lists are deleted from
  // the keep list, and then the function proceeds as if --remove-clusters
  // wasn't specified); this is tracked by which of
  // sorted_keep_ids/sorted_remove_ids is non-NULL.  cluster_kr_ct and
  // max_cluster_kr_len apply to that array.
  char* sorted_keep_ids = NULL;
  char* sorted_remove_ids = NULL;
  uintptr_t max_cluster_id_len = 0;
  uintptr_t assigned_ct = 0;
  uintptr_t cluster_ct = 0;
  Ll_str* cluster_names = NULL;
  uintptr_t* already_seen;
  uintptr_t* ulptr;
  char* cluster_ids;
  uint32_t* cluster_map;
  uint32_t* cluster_starts;
  uint32_t* tmp_cluster_starts;
  uintptr_t topsize_bak;
  uintptr_t line_idx;
  Ll_str* llptr;
  char* sorted_ids;
  uint32_t* id_map;
  char* fam_id;
  char* cluster_name_ptr;
  uintptr_t ulii;
  int32_t sorted_idx;
  uint32_t read_idx;
  uint32_t sample_uidx;
  uint32_t sample_idx;
  uint32_t slen;
  uint32_t uii;
  tbuf[MAXLINELEN - 1] = ' ';
  if (cluster_filter) {
    sample_exclude_new = (uintptr_t*)top_alloc(&topsize, unfiltered_sample_ctl * sizeof(intptr_t));
    if (keep_flattened || keep_fname) {
      if (keep_flattened) {
	cluster_kr_ct = count_and_measure_multistr(keep_flattened, &max_cluster_kr_len);
      }
      if (keep_fname) {
	if (fopen_checked(&infile, keep_fname, "r")) {
	  goto load_clusters_ret_OPEN_FAIL;
	}
	line_idx = 0;
	while (fgets(tbuf, MAXLINELEN, infile)) {
	  line_idx++;
	  if (!tbuf[MAXLINELEN - 1]) {
	    sprintf(logbuf, "Error: Line %" PRIuPTR " of --keep-clusters file is pathologically long.\n", line_idx);
	    goto load_clusters_ret_INVALID_FORMAT_2;
	  }
	  cluster_name_ptr = skip_initial_spaces(tbuf);
	  if (is_eoln_kns(*cluster_name_ptr)) {
	    continue;
	  }
	  slen = strlen_se(cluster_name_ptr);
	  if (slen >= max_cluster_kr_len) {
	    max_cluster_kr_len = slen + 1;
	  }
	  cluster_kr_ct++;
	}
	if (!cluster_kr_ct) {
	  logprint("Error: Empty --keep-clusters file.\n");
	  goto load_clusters_ret_INVALID_FORMAT;
	}
      }
      fill_all_bits(sample_exclude_new, unfiltered_sample_ct);
      sorted_keep_ids = (char*)top_alloc(&topsize, cluster_kr_ct * max_cluster_kr_len);
      if (!sorted_keep_ids) {
	goto load_clusters_ret_NOMEM;
      }
      ulii = 0;
      if (keep_flattened) {
        cluster_name_ptr = keep_flattened;
	do {
	  slen = strlen(cluster_name_ptr) + 1;
	  memcpy(&(sorted_keep_ids[ulii * max_cluster_kr_len]), cluster_name_ptr, slen);
	  ulii++;
	  cluster_name_ptr = &(cluster_name_ptr[slen]);
	} while (*cluster_name_ptr);
      }
      if (keep_fname) {
	rewind(infile);
	while (fgets(tbuf, MAXLINELEN, infile)) {
	  cluster_name_ptr = skip_initial_spaces(tbuf);
	  if (is_eoln_kns(*cluster_name_ptr)) {
	    continue;
	  }
	  slen = strlen_se(cluster_name_ptr);
          memcpyx(&(sorted_keep_ids[ulii * max_cluster_kr_len]), cluster_name_ptr, slen, '\0');
          ulii++;
	}
	if (fclose_null(&infile)) {
	  goto load_clusters_ret_READ_FAIL;
	}
      }
      qsort(sorted_keep_ids, cluster_kr_ct, max_cluster_kr_len, strcmp_casted);
      cluster_kr_ct = collapse_duplicate_ids(sorted_keep_ids, cluster_kr_ct, max_cluster_kr_len, NULL);
      if (remove_flattened || remove_fname) {
	topsize_bak = topsize;
	ulii = (cluster_kr_ct + (BITCT - 1)) / BITCT;
	// track deletions
	already_seen = (uintptr_t*)top_alloc(&topsize, ulii * sizeof(intptr_t));
	if (!already_seen) {
	  goto load_clusters_ret_NOMEM;
	}
	fill_ulong_zero(already_seen, ulii);
	if (remove_flattened) {
	  cluster_name_ptr = remove_flattened;
	  do {
	    slen = strlen(cluster_name_ptr);
	    sorted_idx = bsearch_str(cluster_name_ptr, slen, sorted_keep_ids, max_cluster_kr_len, cluster_kr_ct);
	    if (sorted_idx != -1) {
	      set_bit(already_seen, sorted_idx);
	    }
	    cluster_name_ptr = &(cluster_name_ptr[slen + 1]);
	  } while (*cluster_name_ptr);
	}
	if (remove_fname) {
	  if (fopen_checked(&infile, remove_fname, "r")) {
            goto load_clusters_ret_OPEN_FAIL;
	  }
	  line_idx = 0;
          while (fgets(tbuf, MAXLINELEN, infile)) {
	    line_idx++;
            if (!tbuf[MAXLINELEN - 1]) {
	      sprintf(logbuf, "Error: Line %" PRIuPTR " of --remove-clusters file is pathologically long.\n", line_idx);
	      goto load_clusters_ret_INVALID_FORMAT_2;
	    }
	    cluster_name_ptr = skip_initial_spaces(tbuf);
	    if (is_eoln_kns(*cluster_name_ptr)) {
	      continue;
	    }
	    sorted_idx = bsearch_str_nl(cluster_name_ptr, sorted_keep_ids, max_cluster_kr_len, cluster_kr_ct);
	    if (sorted_idx != -1) {
	      set_bit(already_seen, sorted_idx);
	    }
	  }
	}
        uii = next_set(already_seen, 0, cluster_kr_ct);
	if (uii < cluster_kr_ct) {
	  read_idx = uii + 1;
          for (read_idx = uii + 1; read_idx < cluster_kr_ct; read_idx++) {
            if (!IS_SET(already_seen, read_idx)) {
              strcpy(&(sorted_keep_ids[uii * max_cluster_kr_len]), &(sorted_keep_ids[read_idx * max_cluster_kr_len]));
              uii++;
	    }
	  }
          cluster_kr_ct = uii;
	}
	topsize = topsize_bak;
      }
    } else {
      memcpy(sample_exclude_new, sample_exclude, unfiltered_sample_ctl * sizeof(intptr_t));
      if (remove_flattened) {
        cluster_kr_ct += count_and_measure_multistr(remove_flattened, &max_cluster_kr_len);
      }
      if (remove_fname) {
	if (fopen_checked(&infile, remove_fname, "r")) {
	  goto load_clusters_ret_OPEN_FAIL;
	}
	line_idx = 0;
	while (fgets(tbuf, MAXLINELEN, infile)) {
	  line_idx++;
	  if (!tbuf[MAXLINELEN - 1]) {
	    sprintf(logbuf, "Error: Line %" PRIuPTR " of --remove-clusters file is pathologically long.\n", line_idx);
	    goto load_clusters_ret_INVALID_FORMAT_2;
	  }
	  cluster_name_ptr = skip_initial_spaces(tbuf);
	  if (is_eoln_kns(*cluster_name_ptr)) {
	    continue;
	  }
	  slen = strlen_se(cluster_name_ptr);
	  if (slen >= max_cluster_kr_len) {
	    max_cluster_kr_len = slen + 1;
	  }
	  cluster_kr_ct++;
	}
      }
      if (cluster_kr_ct) {
	sorted_remove_ids = (char*)top_alloc(&topsize, cluster_kr_ct * max_cluster_kr_len);
	if (!sorted_remove_ids) {
	  goto load_clusters_ret_NOMEM;
	}
	ulii = 0;
	if (remove_flattened) {
	  cluster_name_ptr = remove_flattened;
	  do {
	    slen = strlen(cluster_name_ptr) + 1;
	    memcpy(&(sorted_remove_ids[ulii * max_cluster_kr_len]), cluster_name_ptr, slen);
	    ulii++;
	    cluster_name_ptr = &(cluster_name_ptr[slen]);
	  } while (*cluster_name_ptr);
	}
	if (remove_fname) {
	  rewind(infile);
	  while (fgets(tbuf, MAXLINELEN, infile)) {
	    cluster_name_ptr = skip_initial_spaces(tbuf);
	    if (is_eoln_kns(*cluster_name_ptr)) {
	      continue;
	    }
	    slen = strlen_se(cluster_name_ptr);
	    memcpyx(&(sorted_remove_ids[ulii * max_cluster_kr_len]), cluster_name_ptr, slen, '\0');
	    ulii++;
	  }
	}
	qsort(sorted_remove_ids, cluster_kr_ct, max_cluster_kr_len, strcmp_casted);
        cluster_kr_ct = collapse_duplicate_ids(sorted_remove_ids, cluster_kr_ct, max_cluster_kr_len, NULL);
      }
      if (infile) {
	if (fclose_null(&infile)) {
	  goto load_clusters_ret_READ_FAIL;
	}
      }
    }
  }

  if (fname) {
    sorted_ids = (char*)top_alloc(&topsize, sample_ct * max_sample_id_len);
    if (!sorted_ids) {
      goto load_clusters_ret_NOMEM;
    }
    id_map = (uint32_t*)top_alloc(&topsize, sample_ct * sizeof(int32_t));
    if (!id_map) {
      goto load_clusters_ret_NOMEM;
    }
    topsize_bak = topsize;
    already_seen = (uintptr_t*)top_alloc(&topsize, sample_ctl * sizeof(intptr_t));
    if (!already_seen) {
      goto load_clusters_ret_NOMEM;
    }
    fill_ulong_zero(already_seen, sample_ctl);
    wkspace_left -= topsize;
    retval = sort_item_ids_noalloc(sorted_ids, id_map, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, 0, 0, strcmp_deref);
    wkspace_left += topsize;
    if (retval) {
      goto load_clusters_ret_1;
    }

    // two-pass load
    // 1. load cluster names, track longest length, validate format, verify no
    //    sample ID appears multiple times
    // intermission. sort cluster names, purge duplicates, allocate memory for
    //               return values
    // 2. populate return arrays
    if (fopen_checked(&infile, fname, "r")) {
      goto load_clusters_ret_OPEN_FAIL;
    }
    if (!mwithin_col) {
      mwithin_col = 1;
    }
    line_idx = 0;
    while (fgets(tbuf, MAXLINELEN, infile)) {
      line_idx++;
      if (!tbuf[MAXLINELEN - 1]) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --within file is pathologically long.\n", line_idx);
	goto load_clusters_ret_INVALID_FORMAT_2;
      }
      fam_id = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*fam_id)) {
	continue;
      }
      if (bsearch_read_fam_indiv(idbuf, sorted_ids, max_sample_id_len, sample_ct, fam_id, &cluster_name_ptr, &sorted_idx)) {
	goto load_clusters_ret_MISSING_TOKENS;
      }
      if (sorted_idx == -1) {
	continue;
      }
      if (is_set(already_seen, sorted_idx)) {
	*strchr(idbuf, '\t') = ' ';
	LOGPREPRINTFWW("Error: ID '%s' appears multiple times in --within file.\n", idbuf);
	goto load_clusters_ret_INVALID_FORMAT_2;
      }
      if (mwithin_col > 1) {
	cluster_name_ptr = next_token_mult(cluster_name_ptr, mwithin_col - 1);
      }
      if (no_more_tokens_kns(cluster_name_ptr)) {
	goto load_clusters_ret_MISSING_TOKENS;
      }
      set_bit(already_seen, sorted_idx);
      slen = strlen_se(cluster_name_ptr);
      if ((!keep_na) && (slen == 2) && (!memcmp(cluster_name_ptr, "NA", 2))) {
	// postponed to here because, even without 'keep-NA', we do not want to
	// ignore cluster=NA lines for the purpose of detecting duplicate
	// samples
	continue;
      }
      // cluster won't exist because of --keep-clusters/--remove-clusters?
      if ((sorted_keep_ids && (bsearch_str(cluster_name_ptr, slen, sorted_keep_ids, max_cluster_kr_len, cluster_kr_ct) == -1)) || (sorted_remove_ids && (bsearch_str(cluster_name_ptr, slen, sorted_remove_ids, max_cluster_kr_len, cluster_kr_ct) != -1))) {
	continue;
      }
      if (slen >= max_cluster_id_len) {
	max_cluster_id_len = slen + 1;
      }
      cluster_name_ptr[slen] = '\0';
      // do NOT optimize common case because current logic uses
      // collapse_duplicate_ids() last parameter to determine cluster sizes
      llptr = top_alloc_llstr(&topsize, slen + 1);
      if (!llptr) {
	goto load_clusters_ret_NOMEM;
      }
      llptr->next = cluster_names;
      memcpy(llptr->ss, cluster_name_ptr, slen + 1);
      cluster_names = llptr;
      assigned_ct++;
    }
    if (!feof(infile)) {
      goto load_clusters_ret_READ_FAIL;
    }
    if (cluster_names) {
      if (max_cluster_id_len > MAX_ID_LEN_P1) {
	logprint("Error: Cluster IDs are limited to " MAX_ID_LEN_STR " characters.\n");
	goto load_clusters_ret_INVALID_FORMAT;
      }
      *max_cluster_id_len_ptr = max_cluster_id_len;
      wkspace_left -= topsize;
      if (wkspace_alloc_c_checked(cluster_ids_ptr, assigned_ct * max_cluster_id_len)) {
	goto load_clusters_ret_NOMEM2;
      }
      cluster_ids = *cluster_ids_ptr;
      for (ulii = 0; ulii < assigned_ct; ulii++) {
	strcpy(&(cluster_ids[ulii * max_cluster_id_len]), cluster_names->ss);
	cluster_names = cluster_names->next;
      }
      // deallocate cluster ID linked list and duplicate sample ID detector
      // from top of stack, allocate cluster size tracker
      wkspace_left += topsize;
      topsize = topsize_bak;
      tmp_cluster_starts = (uint32_t*)top_alloc(&topsize, assigned_ct * sizeof(int32_t));
      if (!tmp_cluster_starts) {
	goto load_clusters_ret_NOMEM;
      }
      wkspace_left -= topsize;
      // may as well use natural sort of cluster names
      qsort(cluster_ids, assigned_ct, max_cluster_id_len, strcmp_natural);
      cluster_ct = collapse_duplicate_ids(cluster_ids, assigned_ct, max_cluster_id_len, tmp_cluster_starts);
      *cluster_ct_ptr = cluster_ct;
      wkspace_shrink_top(cluster_ids, cluster_ct * max_cluster_id_len);
      if (wkspace_alloc_ui_checked(cluster_map_ptr, assigned_ct * sizeof(int32_t)) ||
	  wkspace_alloc_ui_checked(cluster_starts_ptr, (cluster_ct + 1) * sizeof(int32_t))) {
	goto load_clusters_ret_NOMEM2;
      }
      wkspace_left += topsize;
      cluster_map = *cluster_map_ptr;
      cluster_starts = *cluster_starts_ptr;
      memcpy(cluster_starts, tmp_cluster_starts, cluster_ct * sizeof(int32_t));
      cluster_starts[cluster_ct] = assigned_ct;
      rewind(infile);
      // second pass
      while (fgets(tbuf, MAXLINELEN, infile)) {
	fam_id = skip_initial_spaces(tbuf);
	if (is_eoln_kns(*fam_id)) {
	  continue;
	}
	bsearch_read_fam_indiv(idbuf, sorted_ids, max_sample_id_len, sample_ct, fam_id, &cluster_name_ptr, &sorted_idx);
	if (sorted_idx == -1) {
	  continue;
	}
	if (mwithin_col > 1) {
	  cluster_name_ptr = next_token_mult(cluster_name_ptr, mwithin_col - 1);
	}
	slen = strlen_se(cluster_name_ptr);
	if ((!keep_na) && (slen == 2) && (!memcmp(cluster_name_ptr, "NA", 2))) {
	  continue;
	}
	sample_uidx = id_map[(uint32_t)sorted_idx];
	if (sorted_keep_ids) {
	  sorted_idx = bsearch_str(cluster_name_ptr, slen, sorted_keep_ids, max_cluster_kr_len, cluster_kr_ct);
	  if (sorted_idx == -1) {
	    continue;
	  }
	  clear_bit(sample_exclude_new, sample_uidx);
	} else if (sorted_remove_ids) {
	  sorted_idx = bsearch_str(cluster_name_ptr, slen, sorted_remove_ids, max_cluster_kr_len, cluster_kr_ct);
	  if (sorted_idx != -1) {
	    set_bit(sample_exclude_new, sample_uidx);
	    continue;
	  }
	}
	cluster_name_ptr[slen] = '\0';
	sorted_idx = bsearch_str_natural(cluster_name_ptr, cluster_ids, max_cluster_id_len, cluster_ct);
	uii = tmp_cluster_starts[(uint32_t)sorted_idx];
	tmp_cluster_starts[(uint32_t)sorted_idx] += 1;
	cluster_map[uii] = sample_uidx;
      }
      if (!feof(infile)) {
	goto load_clusters_ret_READ_FAIL;
      }
      for (ulii = 0; ulii < cluster_ct; ulii++) {
	if (cluster_starts[ulii + 1] - cluster_starts[ulii] > 1) {
#ifdef __cplusplus
	  std::sort(&(cluster_map[cluster_starts[ulii]]), &(cluster_map[cluster_starts[ulii + 1]]));
#else
	  qsort(&(cluster_map[cluster_starts[ulii]]), cluster_starts[ulii + 1] - cluster_starts[ulii], sizeof(int32_t), intcmp);
#endif
	}
      }
      LOGPRINTF("--within: %" PRIuPTR " cluster%s loaded, covering a total of %" PRIuPTR " %s.\n", cluster_ct, (cluster_ct == 1)? "" : "s", assigned_ct, species_str(assigned_ct));
    } else {
      if (sorted_keep_ids) {
	logprint("Error: No samples named in --within file remain in the current analysis, so\n--keep-clusters/--keep-cluster-names excludes everyone.\n");
	goto load_clusters_ret_INVALID_FORMAT;
      }
      logprint("Warning: No samples named in --within file remain in the current analysis.\n");
      goto load_clusters_ret_1;
    }
  } else {
    // --family
    // 1. determine max FID len (might be overestimate if
    //    --keep-clusters/--remove-clusters names FIDs that aren't actually
    //    present, but that isn't a big deal)
    // 2. allocate buffer, copy over
    // 3. natural sort, remove duplicates, shrink buffer
    // 4. initialize other data structures
    if (max_cluster_id_len > MAX_ID_LEN_P1) {
      // max FID len was previously checked
      logprint("Error: Cluster IDs are limited to " MAX_ID_LEN_STR " characters.\n");
      goto load_clusters_ret_INVALID_FORMAT;
    }

    for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_unsafe_ck(sample_exclude, &sample_uidx);
      cluster_name_ptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      slen = (uintptr_t)((char*)memchr(cluster_name_ptr, '\t', max_sample_id_len) - cluster_name_ptr);
      if (sorted_keep_ids) {
	sorted_idx = bsearch_str(cluster_name_ptr, slen, sorted_keep_ids, max_cluster_kr_len, cluster_kr_ct);
	if (sorted_idx == -1) {
	  continue;
	}
	clear_bit(sample_exclude_new, sample_uidx);
      } else if (sorted_remove_ids) {
	sorted_idx = bsearch_str(cluster_name_ptr, slen, sorted_remove_ids, max_cluster_kr_len, cluster_kr_ct);
	if (sorted_idx != -1) {
	  set_bit(sample_exclude_new, sample_uidx);
	}
      }
      if (slen >= max_cluster_id_len) {
	max_cluster_id_len = slen + 1;
      }
      assigned_ct++;
    }
    if (!assigned_ct) {
      logprint("Error: --keep-clusters/--keep-cluster-names excludes everyone.\n");
      goto load_clusters_ret_INVALID_FORMAT;
    }
    *max_cluster_id_len_ptr = max_cluster_id_len;
    wkspace_left -= topsize;
    if (wkspace_alloc_c_checked(cluster_ids_ptr, assigned_ct * max_cluster_id_len)) {
      goto load_clusters_ret_NOMEM2;
    }
    cluster_ids = *cluster_ids_ptr;
    if (sample_exclude_new) {
      ulptr = sample_exclude_new;
    } else {
      ulptr = sample_exclude;
    }
    for (sample_uidx = 0, sample_idx = 0; sample_idx < assigned_ct; sample_uidx++, sample_idx++) {
      next_unset_unsafe_ck(ulptr, &sample_uidx);
      cluster_name_ptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      slen = (uintptr_t)((char*)memchr(cluster_name_ptr, '\t', max_sample_id_len) - cluster_name_ptr);
      memcpyx(&(cluster_ids[sample_idx * max_cluster_id_len]), cluster_name_ptr, slen, '\0');
    }
    wkspace_left += topsize;
    tmp_cluster_starts = (uint32_t*)top_alloc(&topsize, assigned_ct * sizeof(int32_t));
    if (!tmp_cluster_starts) {
      goto load_clusters_ret_NOMEM;
    }
    wkspace_left -= topsize;
    qsort(cluster_ids, assigned_ct, max_cluster_id_len, strcmp_natural);
    cluster_ct = collapse_duplicate_ids(cluster_ids, assigned_ct, max_cluster_id_len, tmp_cluster_starts);
    *cluster_ct_ptr = cluster_ct;
    wkspace_shrink_top(cluster_ids, cluster_ct * max_cluster_id_len);
    if (wkspace_alloc_ui_checked(cluster_map_ptr, assigned_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(cluster_starts_ptr, (cluster_ct + 1) * sizeof(int32_t))) {
      goto load_clusters_ret_NOMEM2;
    }
    wkspace_left += topsize;
    cluster_map = *cluster_map_ptr;
    cluster_starts = *cluster_starts_ptr;
    memcpy(cluster_starts, tmp_cluster_starts, cluster_ct * sizeof(int32_t));
    cluster_starts[cluster_ct] = assigned_ct;
    for (sample_uidx = 0, sample_idx = 0; sample_idx < assigned_ct; sample_uidx++, sample_idx++) {
      next_unset_unsafe_ck(ulptr, &sample_uidx);
      cluster_name_ptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      memcpyx(tbuf, cluster_name_ptr, (uintptr_t)((char*)memchr(cluster_name_ptr, '\t', max_cluster_id_len) - cluster_name_ptr), '\0');
      sorted_idx = bsearch_str_natural(tbuf, cluster_ids, max_cluster_id_len, cluster_ct);
      uii = tmp_cluster_starts[(uint32_t)sorted_idx];
      tmp_cluster_starts[(uint32_t)sorted_idx] += 1;
      cluster_map[uii] = sample_uidx;
    }
    LOGPRINTF("--family: %" PRIuPTR " cluster%s defined.\n", cluster_ct, (cluster_ct == 1)? "" : "s");
  }
  if (sample_exclude_new) {
    ulii = popcount_longs(sample_exclude_new, unfiltered_sample_ctl);
    if (ulii != sample_exclude_ct) {
      memcpy(sample_exclude, sample_exclude_new, unfiltered_sample_ctl * sizeof(intptr_t));
      *sample_exclude_ct_ptr = ulii;
      ulii -= sample_exclude_ct;
      LOGPRINTF("%" PRIuPTR " %s removed by cluster filter(s).\n", ulii, species_str(ulii));
    }
  }

  while (0) {
  load_clusters_ret_NOMEM2:
    wkspace_left += topsize;
  load_clusters_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_clusters_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_clusters_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_clusters_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --within file has fewer tokens than expected.\n", line_idx);
  load_clusters_ret_INVALID_FORMAT_2:
    logprintb();
  load_clusters_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 load_clusters_ret_1:
  if (retval) {
    wkspace_reset(wkspace_mark);
  }
  fclose_cond(infile);
  return retval;
}

void fill_unfiltered_sample_to_cluster(uintptr_t unfiltered_sample_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* sample_to_cluster) {
  uint32_t* cluster_map_pos = cluster_map;
  uint32_t* cluster_end_ptr;
  uint32_t cluster_idx;
  // 0xffffffffU cluster index = unassigned
  fill_uint_one(sample_to_cluster, unfiltered_sample_ct);
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    cluster_end_ptr = &(cluster_map[cluster_starts[cluster_idx + 1]]);
    do {
      sample_to_cluster[*cluster_map_pos] = cluster_idx;
    } while (++cluster_map_pos < cluster_end_ptr);
  }
}

int32_t fill_sample_to_cluster(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* sample_to_cluster, uint32_t* late_clidx_to_sample_uidx) {
  // If late_clidx_to_sample_uidx is not NULL, all samples not in a loaded
  // cluster are given their own cluster, and late_clidx_to_sample_uidx is
  // filled with the cluster index -> sample uidx mapping.
  // (Yes, this is a strange interface; it may be switched to filtered sample
  // indexes later.)
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t* cluster_map_pos = cluster_map;
  int32_t retval = 0;
  uint32_t* uidx_to_idx;
  uint32_t* cluster_end_ptr;
  uint32_t cluster_idx;
  uint32_t sample_uidx;
  uint32_t sample_idx;
  if (wkspace_alloc_ui_checked(&uidx_to_idx, unfiltered_sample_ct * sizeof(int32_t))) {
    goto fill_sample_to_cluster_ret_NOMEM;
  }
  fill_uidx_to_idx(sample_exclude, unfiltered_sample_ct, sample_ct, uidx_to_idx);
  fill_uint_one(sample_to_cluster, sample_ct);
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    cluster_end_ptr = &(cluster_map[cluster_starts[cluster_idx + 1]]);
    do {
      sample_to_cluster[uidx_to_idx[*cluster_map_pos]] = cluster_idx;
    } while (++cluster_map_pos < cluster_end_ptr);
  }
  if (late_clidx_to_sample_uidx && (cluster_starts[cluster_ct] < sample_ct)) {
    sample_uidx = 0;
    for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      sample_uidx = next_unset_unsafe(sample_exclude, sample_uidx);
      if (sample_to_cluster[sample_idx] == 0xffffffffU) {
	late_clidx_to_sample_uidx[cluster_idx - cluster_ct] = sample_uidx;
	sample_to_cluster[sample_idx] = cluster_idx++;
      }
    }
  }
  while (0) {
  fill_sample_to_cluster_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t write_clusters(char* outname, char* outname_end, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uint32_t omit_unassigned, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t sample_uidx = 0;
  int32_t retval = 0;
  uint32_t* sample_to_cluster;
  char* sample_id_ptr;
  char* bufptr;
  uintptr_t sample_idx;
  uint32_t cluster_idx;
  uint32_t slen;
  if (wkspace_alloc_ui_checked(&sample_to_cluster, unfiltered_sample_ct * sizeof(int32_t))) {
    goto write_cluster_ret_NOMEM;
  }
  fill_unfiltered_sample_to_cluster(unfiltered_sample_ct, cluster_ct, cluster_map, cluster_starts, sample_to_cluster);
  memcpy(outname_end, ".clst", 6);
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_cluster_ret_OPEN_FAIL;
  }
  for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
    sample_uidx = next_unset_unsafe(sample_exclude, sample_uidx);
    cluster_idx = sample_to_cluster[sample_uidx];
    if ((!omit_unassigned) || (cluster_idx != 0xffffffffU)) {
      sample_id_ptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      slen = strlen_se(sample_id_ptr);
      bufptr = memcpyax(tbuf, sample_id_ptr, slen, ' ');
      bufptr = strcpyax(bufptr, &(sample_id_ptr[slen + 1]), ' ');
      if (cluster_idx != 0xffffffffU) {
        bufptr = strcpyax(bufptr, &(cluster_ids[cluster_idx * max_cluster_id_len]), '\n');
      } else {
        bufptr = memcpyl3a(bufptr, "NA\n");
      }
      if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	goto write_cluster_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto write_cluster_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--write-cluster: Pruned cluster assignments written to %s .\n", outname);
  while (0) {
  write_cluster_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_cluster_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_cluster_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t extract_clusters(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, char* cluster_names_flattened, char* clusters_fname, uintptr_t** new_sample_exclude_ptr, uintptr_t* new_sample_ct_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t line_idx = 0;
  int32_t retval = 0;
  char* bufptr;
  uintptr_t* new_sample_exclude;
  uint32_t* uiptr;
  uint32_t* cluster_end;
  uint32_t slen;
  uint32_t sample_uidx;
  int32_t ii;
  if (wkspace_alloc_ul_checked(new_sample_exclude_ptr, unfiltered_sample_ctl * sizeof(intptr_t))) {
    goto extract_clusters_ret_NOMEM;
  }
  new_sample_exclude = *new_sample_exclude_ptr;
  wkspace_mark = wkspace_base;
  fill_all_bits(new_sample_exclude, unfiltered_sample_ct);
  if (cluster_names_flattened) {
    bufptr = cluster_names_flattened;
    do {
      slen = strlen(bufptr);
      if (slen < max_cluster_id_len) {
        ii = bsearch_str_natural(bufptr, cluster_ids, max_cluster_id_len, cluster_ct);
        if (ii != -1) {
	  uiptr = &(cluster_map[cluster_starts[(uint32_t)ii]]);
	  cluster_end = &(cluster_map[cluster_starts[1 + (uint32_t)ii]]);
	  while (uiptr < cluster_end) {
	    sample_uidx = *uiptr++;
	    if (!is_set(sample_exclude, sample_uidx)) {
	      clear_bit(new_sample_exclude, sample_uidx);
	    }
	  }
	}
      }
      bufptr = &(bufptr[slen + 1]);
    } while (*bufptr);
  }
  if (clusters_fname) {
    if (fopen_checked(&infile, clusters_fname, "r")) {
      goto extract_clusters_ret_OPEN_FAIL;
    }
    tbuf[MAXLINELEN - 1] = ' ';
    while (fgets(tbuf, MAXLINELEN, infile)) {
      line_idx++;
      if (!tbuf[MAXLINELEN - 1]) {
	LOGPREPRINTFWW(logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, clusters_fname);
        goto extract_clusters_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr)) {
        continue;
      }
      slen = strlen_se(bufptr);
      if (slen >= max_cluster_id_len) {
	continue;
      }
      bufptr[slen] = '\0';
      ii = bsearch_str_natural(bufptr, cluster_ids, max_cluster_id_len, cluster_ct);
      if (ii != -1) {
	uiptr = &(cluster_map[cluster_starts[(uint32_t)ii]]);
	cluster_end = &(cluster_map[cluster_starts[1 + (uint32_t)ii]]);
	while (uiptr < cluster_end) {
	  sample_uidx = *uiptr++;
	  if (!is_set(sample_exclude, sample_uidx)) {
	    clear_bit(new_sample_exclude, sample_uidx);
	  }
	}
      }
    }
    if (fclose_null(&infile)) {
      goto extract_clusters_ret_READ_FAIL;
    }
  }
  *new_sample_ct_ptr = unfiltered_sample_ct - popcount_longs(new_sample_exclude, unfiltered_sample_ctl);
  while (0) {
  extract_clusters_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  extract_clusters_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  extract_clusters_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  extract_clusters_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

uint32_t no_size1(uint32_t cluster_ct, uint32_t* cluster_starts) {
  uint32_t cluster_idx;
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    if (cluster_starts[cluster_idx + 1] - cluster_starts[cluster_idx] < 2) {
      return 0;
    }
  }
  return 1;
}

void populate_cluster_case_cts(uintptr_t* pheno_c, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* cluster_case_cts) {
  uint32_t cluster_idx;
  uint32_t cluster_end;
  uint32_t map_idx;
  uint32_t case_ct;
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    case_ct = 0;
    cluster_end = cluster_starts[cluster_idx + 1];
    for (map_idx = cluster_starts[cluster_idx]; map_idx < cluster_end; map_idx++) {
      case_ct += IS_SET(pheno_c, cluster_map[map_idx]);
    }
    cluster_case_cts[cluster_idx] = case_ct;
  }
}

void adjust_cc_perm_preimage(uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* cluster_case_cts, uintptr_t* cluster_cc_perm_preimage, uint32_t is_perm1) {
  uint32_t cluster_idx;
  uint32_t map_idx;
  uint32_t cluster_end;
  if (!is_perm1) {
    for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
      map_idx = cluster_starts[cluster_idx];
      cluster_end = cluster_starts[cluster_idx + 1];
      if (cluster_case_cts[cluster_idx] * 2 < cluster_end - map_idx) {
	do {
	  CLEAR_BIT_DBL(cluster_cc_perm_preimage, cluster_map[map_idx]);
	} while (++map_idx < cluster_end);
      } else {
	do {
	  SET_BIT_DBL(cluster_cc_perm_preimage, cluster_map[map_idx]);
	} while (++map_idx < cluster_end);
      }
    }
  } else {
    for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
      map_idx = cluster_starts[cluster_idx];
      cluster_end = cluster_starts[cluster_idx + 1];
      if (cluster_case_cts[cluster_idx] * 2 < cluster_end - map_idx) {
	do {
	  CLEAR_BIT(cluster_cc_perm_preimage, cluster_map[map_idx]);
	} while (++map_idx < cluster_end);
      } else {
	do {
	  SET_BIT(cluster_cc_perm_preimage, cluster_map[map_idx]);
	} while (++map_idx < cluster_end);
      }
    }
  }
}

int32_t cluster_include_and_reindex(uintptr_t unfiltered_sample_ct, uintptr_t* sample_include, uint32_t remove_size1, uintptr_t* pheno_c, uintptr_t sample_ct, uint32_t is_perm1, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* new_cluster_ct_ptr, uint32_t** new_cluster_map_ptr, uint32_t** new_cluster_starts_ptr, uint32_t** cluster_case_cts_ptr, uintptr_t** cluster_cc_perm_preimage_ptr) {
  // If any samples are excluded, this converts a unfiltered-index cluster map
  // to a collapsed-index cluster map (suitable for passing to the per-cluster
  // permutation function), allocating space for the latter.  Otherwise, this
  // just sets the new cluster map to point to the old one.
  //
  // remove_size1 should be 0 or 1.  If it's 1, all size-1 clusters are
  // removed.
  //
  // If pheno_c is set, this also allocates and populates an array of
  // per-cluster case counts, and a permutation preimage.
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t old_assigned_ct = cluster_starts[cluster_ct];
  uint32_t new_cluster_ct = 0;
  uint32_t sample_uidx = 0;
  uint32_t case_ct = 0;
  uint32_t assigned_ct = 0;
  uintptr_t* cluster_cc_perm_preimage = NULL;
  uint32_t* cluster_case_cts = NULL;
  uint32_t* new_cluster_map;
  uint32_t* new_cluster_starts;
  uint32_t* uidx_to_idx;
  uint32_t cluster_assigned_ct;
  uint32_t cluster_idx;
  uint32_t cluster_end;
  uint32_t cluster_read_idx;
  uint32_t map_idx;
  uint32_t map_read_idx;
  uint32_t last_case_ct_incr;
  uint32_t shrink_map;
  if (pheno_c) {
    if (wkspace_alloc_ul_checked(cluster_cc_perm_preimage_ptr, (2 - is_perm1) * ((sample_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t))) {
      goto cluster_include_and_reindex_ret_NOMEM;
    }
    cluster_cc_perm_preimage = *cluster_cc_perm_preimage_ptr;
    if (!is_perm1) {
      vec_collapse_init(pheno_c, unfiltered_sample_ct, sample_include, sample_ct, cluster_cc_perm_preimage);
    } else {
      collapse_copy_bitarr_incl(unfiltered_sample_ct, pheno_c, sample_include, sample_ct, cluster_cc_perm_preimage);
    }
  }
  if ((sample_ct == unfiltered_sample_ct) && ((!remove_size1) || no_size1(cluster_ct, cluster_starts))) {
    *new_cluster_map_ptr = cluster_map;
    *new_cluster_ct_ptr = cluster_ct;
    *new_cluster_starts_ptr = cluster_starts;
    if (pheno_c) {
      if (wkspace_alloc_ui_checked(cluster_case_cts_ptr, cluster_ct * sizeof(int32_t))) {
	goto cluster_include_and_reindex_ret_NOMEM;
      }
      populate_cluster_case_cts(pheno_c, cluster_ct, cluster_map, cluster_starts, *cluster_case_cts_ptr);
      adjust_cc_perm_preimage(cluster_ct, cluster_map, cluster_starts, *cluster_case_cts_ptr, cluster_cc_perm_preimage, is_perm1);
    }
    return 0;
  }
  // scan to determine new memory allocation sizes
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    cluster_assigned_ct = 0;
    cluster_end = cluster_starts[cluster_idx + 1];
    for (map_idx = cluster_starts[cluster_idx]; map_idx < cluster_end; map_idx++) {
      cluster_assigned_ct += is_set(sample_include, cluster_map[map_idx]);
    }
    if (cluster_assigned_ct > remove_size1) {
      new_cluster_ct++;
      assigned_ct += cluster_assigned_ct;
    }
  }
  // possibly +1 to simplify remove_size1 logic
  if (wkspace_alloc_ui_checked(new_cluster_map_ptr, (assigned_ct + remove_size1) * sizeof(int32_t))) {
    goto cluster_include_and_reindex_ret_NOMEM;
  }
  new_cluster_map = *new_cluster_map_ptr;
  shrink_map = (assigned_ct < old_assigned_ct);
  if (shrink_map) {
    if (wkspace_alloc_ui_checked(new_cluster_starts_ptr, (new_cluster_ct + 1) * sizeof(int32_t))) {
      goto cluster_include_and_reindex_ret_NOMEM;
    }
    new_cluster_starts = *new_cluster_starts_ptr;
    new_cluster_starts[0] = 0;
  } else {
    new_cluster_starts = cluster_starts;
    *new_cluster_starts_ptr = cluster_starts;
  }
  if (pheno_c) {
    if (wkspace_alloc_ui_checked(cluster_case_cts_ptr, new_cluster_ct * sizeof(int32_t))) {
      goto cluster_include_and_reindex_ret_NOMEM;
    }
    cluster_case_cts = *cluster_case_cts_ptr;
  }
  if (wkspace_alloc_ui_checked(&uidx_to_idx, unfiltered_sample_ct * sizeof(int32_t))) {
    goto cluster_include_and_reindex_ret_NOMEM;
  }
  fill_uidx_to_idx_incl(sample_include, unfiltered_sample_ct, sample_ct, uidx_to_idx);
  *new_cluster_ct_ptr = new_cluster_ct;
  cluster_read_idx = 1;
  map_idx = 0;
  cluster_idx = 0;
  cluster_end = 0;
  last_case_ct_incr = 0;
  // walk through cluster_map and:
  // * skip excluded sample indices, and reindex the rest when copying to
  //   new_cluster_map[]
  // * update new_cluster_starts[] if necessary
  // * if remove_size1, also delete clusters that are now size 1
  for (map_read_idx = 0; map_read_idx < old_assigned_ct; map_read_idx++) {
    sample_uidx = cluster_map[map_read_idx];
    if (!IS_SET(sample_include, sample_uidx)) {
      continue;
    }
    if (map_read_idx >= cluster_end) {
      if (cluster_idx) {
        if ((!remove_size1) || (map_idx - new_cluster_starts[cluster_idx - 1] > 1)) {
	  if (pheno_c) {
	    cluster_case_cts[cluster_idx - 1] = case_ct + last_case_ct_incr;
	  }
	  if (shrink_map) {
	    new_cluster_starts[cluster_idx] = map_idx;
	  }
	  case_ct = 0;
	  cluster_idx++;
        } else {
          map_idx--;
        }
	last_case_ct_incr = 0;
      } else {
        cluster_idx = 1;
      }
      do {
	cluster_end = cluster_starts[cluster_read_idx++];
      } while (cluster_end <= map_read_idx);
    }
    if (pheno_c) {
      case_ct += last_case_ct_incr;
      last_case_ct_incr = IS_SET(pheno_c, sample_uidx);
    }
    new_cluster_map[map_idx++] = uidx_to_idx[sample_uidx];
  }
  if (new_cluster_ct && (!(remove_size1 && (new_cluster_starts[cluster_idx - 1] == map_idx - 1)))) {
    // fill in last array elements, but guard against "last cluster is of size
    // 1, and we're removing it" special case
    if (pheno_c) {
      cluster_case_cts[new_cluster_ct - 1] = case_ct + last_case_ct_incr;
    }
    if (shrink_map) {
      new_cluster_starts[new_cluster_ct] = map_idx;
    }
  }
  if (pheno_c && new_cluster_ct) {
    adjust_cc_perm_preimage(new_cluster_ct, new_cluster_map, new_cluster_starts, cluster_case_cts, cluster_cc_perm_preimage, is_perm1);
  }
  wkspace_reset(uidx_to_idx);
  return 0;
 cluster_include_and_reindex_ret_NOMEM:
  wkspace_reset(wkspace_mark);
  return RET_NOMEM;
}

int32_t cluster_alloc_and_populate_magic_nums(uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t** tot_quotients_ptr, uint64_t** totq_magics_ptr, uint32_t** totq_preshifts_ptr, uint32_t** totq_postshifts_ptr, uint32_t** totq_incrs_ptr) {
  // assumes all clusters are of size > 1
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t* tot_quotients;
  uint64_t* totq_magics;
  uint32_t* totq_preshifts;
  uint32_t* totq_postshifts;
  uint32_t* totq_incrs;
  uint32_t cluster_idx;
  if (wkspace_alloc_ui_checked(tot_quotients_ptr, cluster_ct * sizeof(int32_t)) ||
      wkspace_alloc_ull_checked(totq_magics_ptr, cluster_ct * sizeof(int64_t)) ||
      wkspace_alloc_ui_checked(totq_preshifts_ptr, cluster_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(totq_postshifts_ptr, cluster_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(totq_incrs_ptr, cluster_ct * sizeof(int32_t))) {
    wkspace_reset(wkspace_mark);
    return RET_NOMEM;
  }
  tot_quotients = *tot_quotients_ptr;
  totq_magics = *totq_magics_ptr;
  totq_preshifts = *totq_preshifts_ptr;
  totq_postshifts = *totq_postshifts_ptr;
  totq_incrs = *totq_incrs_ptr;
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    tot_quotients[cluster_idx] = 0x100000000LLU / (cluster_starts[cluster_idx + 1] - cluster_starts[cluster_idx]);
    magic_num(tot_quotients[cluster_idx], &(totq_magics[cluster_idx]), &(totq_preshifts[cluster_idx]), &(totq_postshifts[cluster_idx]), &(totq_incrs[cluster_idx]));
  }
  return 0;
}

int32_t read_dists(char* dist_fname, char* id_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* sample_to_cluster, uint32_t for_cluster_flag, uint32_t is_max_dist, double* dists, uint32_t neighbor_n2, double* neighbor_quantiles, uint32_t* neighbor_qindices) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* dist_file = NULL;
  FILE* id_file = NULL;
  uintptr_t id_entry_ct = sample_ct;
  uintptr_t matching_entry_ct = sample_ct;
  uintptr_t line_idx = 0;
  char* id_buf = &(tbuf[MAXLINELEN]);
  uint64_t* fidx_to_memidx = NULL; // high 32 bits = fidx, low 32 = memidx
  uint32_t is_presorted = cluster_ct? 0 : 1;
  int32_t retval = 0;
  char* sorted_ids;
  uint32_t* id_map;
  char* fam_id;
  double* dptr;
  uint64_t fpos;
  uint64_t fpos2;
  uintptr_t memidx1;
  uintptr_t memidx2;
  uintptr_t fidx1;
  uint64_t fidx2;
  uintptr_t trimem;
  uint64_t trif;
  uintptr_t clidx1;
  uintptr_t clidx2;
  uint64_t ullii;
  uintptr_t ulii;
  uintptr_t uljj;
  double cur_ibs;
  uint32_t uii;
  int32_t ii;
  if (fopen_checked(&dist_file, dist_fname, "rb")) {
    goto read_dists_ret_OPEN_FAIL;
  }
  if (fseeko(dist_file, 0, SEEK_END)) {
    goto read_dists_ret_READ_FAIL;
  }
  if (id_fname) {
    if (wkspace_alloc_ull_checked(&fidx_to_memidx, sample_ct * sizeof(int64_t))) {
      goto read_dists_ret_NOMEM;
    }
    fill_ull_one(fidx_to_memidx, sample_ct);
    if (fopen_checked(&id_file, id_fname, "r")) {
      goto read_dists_ret_OPEN_FAIL;
    }
    retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_sample_ct, sample_exclude, unfiltered_sample_ct - sample_ct, sample_ids, max_sample_id_len, 0, 1, strcmp_deref);
    if (retval) {
      goto read_dists_ret_1;
    }
    id_entry_ct = 0;
    matching_entry_ct = 0;
    tbuf[MAXLINELEN - 1] = ' ';
    while (fgets(tbuf, MAXLINELEN, id_file)) {
      line_idx++;
      if (!tbuf[MAXLINELEN - 1]) {
	LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, id_fname);
        goto read_dists_ret_INVALID_FORMAT_2;
      }
      fam_id = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*fam_id)) {
        continue;
      }
      if (bsearch_read_fam_indiv(id_buf, sorted_ids, max_sample_id_len, sample_ct, fam_id, NULL, &ii)) {
	LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, id_fname);
        goto read_dists_ret_INVALID_FORMAT_2;
      }
      if (ii == -1) {
        is_presorted = 0;
        id_entry_ct++;
        continue;
      }
      uii = id_map[(uint32_t)ii];
      if (uii != id_entry_ct) {
        is_presorted = 0;
      }
      if (fidx_to_memidx[uii] != 0xffffffffffffffffLLU) {
	*strchr(id_buf, '\t') = ' ';
        LOGPREPRINTFWW("Error: ID '%s' appears multiple times in %s.\n", id_buf, id_fname);
        goto read_dists_ret_INVALID_FORMAT_2;
      }
      if (cluster_ct && (!neighbor_n2)) {
	// if cluster_ct && neighbor_n2, best to
	// postpone sample_to_cluster dereference
        fidx_to_memidx[uii] = ((uint64_t)(sample_to_cluster[uii])) | (((uint64_t)id_entry_ct) << 32);
      } else {
        fidx_to_memidx[uii] = ((uint64_t)uii) | (((uint64_t)id_entry_ct) << 32);
      }
      id_entry_ct++;
      matching_entry_ct++;
    }
    if (!feof(id_file)) {
      goto read_dists_ret_READ_FAIL;
    }
    fclose_null(&id_file);
    if (matching_entry_ct < sample_ct) {
      logprint("Error: --read-dists ID file does not contain all samples in current run.\n");
      goto read_dists_ret_INVALID_FORMAT;
    }
  } else if (cluster_ct) {
    if (wkspace_alloc_ull_checked(&fidx_to_memidx, sample_ct * sizeof(int64_t))) {
      goto read_dists_ret_NOMEM;
    }
    if (neighbor_n2) {
      for (id_entry_ct = 0; id_entry_ct < sample_ct; id_entry_ct++) {
	fidx_to_memidx[id_entry_ct] = ((uint64_t)id_entry_ct) | (((uint64_t)id_entry_ct) << 32);
      }
    } else {
      for (id_entry_ct = 0; id_entry_ct < sample_ct; id_entry_ct++) {
	fidx_to_memidx[id_entry_ct] = ((uint64_t)(sample_to_cluster[id_entry_ct])) | (((uint64_t)id_entry_ct) << 32);
      }
    }
  }
  if (!is_presorted) {
#ifdef __cplusplus
    // std::sort is faster than qsort for basic types.  See e.g. Anders
    // Kaseorg's answer to
    // http://www.quora.com/Software-Engineering/Generally-how-much-faster-is-C-compared-to-C++
    std::sort((int64_t*)fidx_to_memidx, (int64_t*)(&(fidx_to_memidx[sample_ct])));
#else
    qsort(fidx_to_memidx, sample_ct, sizeof(int64_t), llcmp);
#endif
  }
  fpos = (((uint64_t)id_entry_ct) * (id_entry_ct - 1)) * (sizeof(double) / 2);
  if (ftello(dist_file) != (int64_t)fpos) {
    LOGPREPRINTFWW("Error: --read-dists expects size of %s to be %" PRIu64 " bytes.\n", dist_fname, fpos);
    goto read_dists_ret_INVALID_FORMAT_2;
  }
  rewind(dist_file);
  if (is_presorted) {
    if (fread(dists, 1, fpos, dist_file) < fpos) {
      goto read_dists_ret_READ_FAIL;
    }
    if (neighbor_n2) {
      dptr = dists;
      for (memidx1 = 1; memidx1 < sample_ct; memidx1++) {
	for (memidx2 = 0; memidx2 < memidx1; memidx2++) {
          update_neighbor(sample_ct, neighbor_n2, memidx1, memidx2, *dptr++, neighbor_quantiles, neighbor_qindices);
	}
      }
    }
  } else {
    fpos = 0;
    if ((!cluster_ct) || (!neighbor_n2)) {
      for (ulii = 1; ulii < sample_ct; ulii++) {
	ullii = fidx_to_memidx[ulii];
	memidx1 = (uintptr_t)(ullii & (ONELU * 0xffffffffU));
	fidx1 = (uintptr_t)(ullii >> 32);
	trimem = (memidx1 * (memidx1 - 1)) / 2;
	trif = (((uint64_t)fidx1) * (fidx1 - 1)) * (sizeof(double) / 2);
	for (uljj = 0; uljj < ulii; uljj++) {
	  ullii = fidx_to_memidx[uljj];
	  memidx2 = (uintptr_t)(ullii & (ONELU * 0xffffffffU));
	  fidx2 = (uint64_t)(ullii >> 32);
	  fpos2 = trif + (fidx2 * sizeof(double));
	  if (fpos2 > fpos) {
	    fpos = fpos2;
	    if (fseeko(dist_file, fpos, SEEK_SET)) {
	      goto read_dists_ret_READ_FAIL;
	    }
	  }
	  if (fread(&cur_ibs, 1, sizeof(double), dist_file) < sizeof(double)) {
	    goto read_dists_ret_READ_FAIL;
	  }
	  if (memidx2 < memidx1) {
	    dists[trimem + memidx2] += cur_ibs;
	  } else if (memidx2 > memidx1) {
	    dists[((memidx2 * (memidx2 - 1)) / 2) + memidx1] += cur_ibs;
	  }
	  if (neighbor_n2) {
	    update_neighbor(sample_ct, neighbor_n2, memidx1, memidx2, cur_ibs, neighbor_quantiles, neighbor_qindices);
	  }
	  fpos += sizeof(double);
	}
      }
    } else {
      for (ulii = 1; ulii < sample_ct; ulii++) {
	ullii = fidx_to_memidx[ulii];
	memidx1 = (uintptr_t)(ullii & (ONELU * 0xffffffffU));
	fidx1 = (uintptr_t)(ullii >> 32);
        clidx1 = sample_to_cluster[memidx1];
	trimem = (clidx1 * (clidx1 - 1)) / 2;
	trif = (((uint64_t)fidx1) * (fidx1 - 1)) * (sizeof(double) / 2);
	for (uljj = 0; uljj < ulii; uljj++) {
	  ullii = fidx_to_memidx[uljj];
	  memidx2 = (uintptr_t)(ullii & (ONELU * 0xffffffffU));
	  fidx2 = (uint64_t)(ullii >> 32);
	  clidx2 = sample_to_cluster[memidx2];
	  fpos2 = trif + (fidx2 * sizeof(double));
	  if (fpos2 > fpos) {
	    fpos = fpos2;
	    if (fseeko(dist_file, fpos, SEEK_SET)) {
	      goto read_dists_ret_READ_FAIL;
	    }
	  }
	  if (fread(&cur_ibs, 1, sizeof(double), dist_file) < sizeof(double)) {
	    goto read_dists_ret_READ_FAIL;
	  }
	  if (!is_max_dist) {
	    if (clidx2 < clidx1) {
	      dists[trimem + clidx2] += cur_ibs;
	    } else if (clidx2 > clidx1) {
	      dists[((clidx2 * (clidx2 - 1)) / 2) + clidx1] += cur_ibs;
	    }
	  } else {
	    if (clidx1 != clidx2) {
	      if (clidx2 < clidx1) {
	        dptr = &(dists[trimem + clidx2]);
	      } else {
		dptr = &(dists[((clidx2 * (clidx2 - 1)) / 2) + clidx1]);
	      }
	      if (cur_ibs < (*dptr)) {
		*dptr = cur_ibs;
	      }
	    }
	  }
	  if (neighbor_n2) {
	    update_neighbor(sample_ct, neighbor_n2, memidx1, memidx2, cur_ibs, neighbor_quantiles, neighbor_qindices);
	  }
	  fpos += sizeof(double);
	}
      }
    }
  }
  if (for_cluster_flag != 2) {
    LOGPRINTF("--read-dists: %" PRIuPTR " values loaded%s.\n", (sample_ct * (sample_ct - 1)) / 2, for_cluster_flag? " for --cluster/--neighbor" : "");
  }
  while (0) {
  read_dists_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  read_dists_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  read_dists_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  read_dists_ret_INVALID_FORMAT_2:
    logprintb();
  read_dists_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 read_dists_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(dist_file);
  fclose_cond(id_file);
  return retval;
}

void update_neighbor(uintptr_t sample_ct, uint32_t neighbor_n2, uintptr_t sample_idx1, uintptr_t sample_idx2, double cur_ibs, double* neighbor_quantiles, uint32_t* neighbor_qindices) {
  uintptr_t exceed_ct;
  uintptr_t cur_write;
  exceed_ct = nonincr_doublearr_leq_stride(&(neighbor_quantiles[sample_idx1]), neighbor_n2, sample_ct, cur_ibs);
  if (exceed_ct < neighbor_n2) {
    for (cur_write = neighbor_n2 - 1; cur_write > exceed_ct; cur_write--) {
      neighbor_quantiles[cur_write * sample_ct + sample_idx1] = neighbor_quantiles[(cur_write - 1) * sample_ct + sample_idx1];
      neighbor_qindices[cur_write * sample_ct + sample_idx1] = neighbor_qindices[(cur_write - 1) * sample_ct + sample_idx1];
    }
    neighbor_quantiles[(exceed_ct * sample_ct) + sample_idx1] = cur_ibs;
    neighbor_qindices[(exceed_ct * sample_ct) + sample_idx1] = sample_idx2;
  }
  exceed_ct = nonincr_doublearr_leq_stride(&(neighbor_quantiles[sample_idx2]), neighbor_n2, sample_ct, cur_ibs);
  if (exceed_ct < neighbor_n2) {
    for (cur_write = neighbor_n2 - 1; cur_write > exceed_ct; cur_write--) {
      neighbor_quantiles[cur_write * sample_ct + sample_idx2] = neighbor_quantiles[(cur_write - 1) * sample_ct + sample_idx2];
      neighbor_qindices[cur_write * sample_ct + sample_idx2] = neighbor_qindices[(cur_write - 1) * sample_ct + sample_idx2];
    }
    neighbor_quantiles[(exceed_ct * sample_ct) + sample_idx2] = cur_ibs;
    neighbor_qindices[(exceed_ct * sample_ct) + sample_idx2] = sample_idx1;
  }
}

int32_t read_genome(char* read_genome_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* cluster_merge_prevented, double* cluster_sorted_ibs, uint32_t neighbor_n2, double* neighbor_quantiles, uint32_t* neighbor_qindices, uint32_t* ppc_fail_counts, double min_ppc, uint32_t is_max_dist, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* sample_to_cluster) {
  unsigned char* wkspace_mark = wkspace_base;
  gzFile gz_infile = NULL;
  uint32_t neighbor_load_quantiles = neighbor_quantiles && cluster_sorted_ibs;
  uint32_t ppc_warning = cluster_merge_prevented? 0 : 1;
  uintptr_t loaded_entry_ct = 0;
  uintptr_t line_idx = 1;
  uint32_t ppc_fail = 0;
  char* idbuf = &(tbuf[MAXLINELEN]);
  char* sorted_ids;
  uint32_t* id_map;
  char* bufptr;
  char* fam_id;
  uint32_t sample_idx1;
  uint32_t sample_idx2;
  double cur_ibs;
  double cur_ppc;
  uintptr_t tcoord;
  uint32_t uii;
  int32_t ii;
  int32_t retval;
  retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_sample_ct, sample_exclude, unfiltered_sample_ct - sample_ct, sample_ids, max_sample_id_len, 0, 1, strcmp_deref);
  if (retval) {
    goto read_genome_ret_1;
  }
  if (gzopen_checked(&gz_infile, read_genome_fname, "rb")) {
    goto read_genome_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  // header line
  do {
    if (!gzgets(gz_infile, tbuf, MAXLINELEN)) {
      goto read_genome_ret_READ_FAIL;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      goto read_genome_ret_LONG_LINE;
    }
    bufptr = skip_initial_spaces(tbuf);
  } while (is_eoln_kns(*bufptr));
  // a little bit of input validation
  if (memcmp(bufptr, "FID1", 4)) {
    logprint("Error: Invalid --read-genome file header line.\n");
    goto read_genome_ret_INVALID_FORMAT;
  }
  while (gzgets(gz_infile, tbuf, MAXLINELEN)) {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      goto read_genome_ret_LONG_LINE;
    }
    fam_id = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*fam_id)) {
      continue;
    }
    if (bsearch_read_fam_indiv(idbuf, sorted_ids, max_sample_id_len, sample_ct, fam_id, &fam_id, &ii)) {
      goto read_genome_ret_MISSING_TOKENS;
    }
    if (ii == -1) {
      continue;
    }
    sample_idx1 = id_map[(uint32_t)ii];
    if (bsearch_read_fam_indiv(idbuf, sorted_ids, max_sample_id_len, sample_ct, fam_id, &bufptr, &ii)) {
      goto read_genome_ret_MISSING_TOKENS;
    }
    if (ii == -1) {
      continue;
    }
    sample_idx2 = id_map[(uint32_t)ii];
    if (sample_idx2 == sample_idx1) {
      sprintf(logbuf, "Error: FID1/IID1 matches FID2/IID2 on line %" PRIuPTR " of --read-genome file.\n", line_idx);
      goto read_genome_ret_INVALID_FORMAT_2;
    }
    bufptr = next_token_mult(bufptr, 7); // distance
    fam_id = next_token(bufptr); // repurposed to PPC test value
    if (no_more_tokens(fam_id)) {
      goto read_genome_ret_MISSING_TOKENS;
    }
    if (min_ppc != 0.0) {
      if (scan_double(fam_id, &cur_ppc)) {
	sprintf(logbuf, "Error: Invalid PPC test value on line %" PRIuPTR " of --read-genome file.\n", line_idx);
	goto read_genome_ret_INVALID_FORMAT_2;
      }
      ppc_fail = (cur_ppc < min_ppc);
      if (ppc_fail && ppc_fail_counts) {
	ppc_fail_counts[sample_idx1] += 1;
	ppc_fail_counts[sample_idx2] += 1;
      }
    }
    if (scan_double(bufptr, &cur_ibs)) {
      sprintf(logbuf, "Error: Invalid IBS value on line %" PRIuPTR " of --read-genome file.\n", line_idx);
      goto read_genome_ret_INVALID_FORMAT_2;
    }
    if (neighbor_load_quantiles) {
      update_neighbor(sample_ct, neighbor_n2, sample_idx1, sample_idx2, cur_ibs, neighbor_quantiles, neighbor_qindices);
    }
    loaded_entry_ct++;
    if (cluster_ct) {
      sample_idx1 = sample_to_cluster[sample_idx1];
      sample_idx2 = sample_to_cluster[sample_idx2];
      if (sample_idx1 == sample_idx2) {
	if (ppc_fail && (!ppc_warning)) {
	  logprint("Warning: Initial cluster assignment violates PPC test constraint.\n");
	  ppc_warning = 1;
	}
	continue;
      }
    }
    if (sample_idx2 < sample_idx1) {
      uii = sample_idx1;
      sample_idx1 = sample_idx2;
      sample_idx2 = uii;
    }
    tcoord = tri_coord_no_diag(sample_idx1, sample_idx2);
    if (ppc_fail) {
      SET_BIT(cluster_merge_prevented, tcoord);
    }
    if (cluster_sorted_ibs) {
      if (!is_max_dist) {
        cluster_sorted_ibs[tcoord] += cur_ibs;
      } else {
	if (cur_ibs < cluster_sorted_ibs[tcoord]) {
          cluster_sorted_ibs[tcoord] = cur_ibs;
	}
      }
    }
  }
  if (!gzeof(gz_infile)) {
    goto read_genome_ret_READ_FAIL;
  }
  if (loaded_entry_ct != (sample_ct * (sample_ct - 1)) / 2) {
    LOGPREPRINTFWW("Error: %s does not include all sample pairs.\n", read_genome_fname);
    goto read_genome_ret_INVALID_FORMAT_2;
  }
  while (0) {
  read_genome_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  read_genome_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  read_genome_ret_MISSING_TOKENS:
    LOGPRINTF("Error: Line %" PRIuPTR " of --read-genome file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  read_genome_ret_LONG_LINE:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --read-genome file is pathologically long.\n", line_idx);
  read_genome_ret_INVALID_FORMAT_2:
    logprintb();
  read_genome_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 read_genome_ret_1:
  wkspace_reset(wkspace_mark);
  gzclose_cond(gz_infile);
  return retval;
}

int32_t cluster_enforce_match(Cluster_info* cp, int32_t missing_pheno, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* sample_to_cluster, uintptr_t* merge_prevented) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* matchfile = NULL;
  FILE* typefile = NULL;
  char* id_buf = &(tbuf[MAXLINELEN]);
  char* missing_str = NULL;
  uintptr_t cur_coord = 0;
  uint32_t cluster_mismatch_warning = 0;
  uint32_t cov_ct = 0;
  uint32_t non_null_cov_ct = 0;
  uint32_t missing_len = 0;
  int32_t retval = 0;
  char intbuf[12];
  char* sorted_ids;
  uint32_t* id_map;
  char** sample_idx_to_match_str;
  double** sample_idx_to_dvals;
  unsigned char* wkspace_mark2;
  unsigned char* cov_type_arr;
  double* tol_arr;
  char* bufptr;
  char* bufptr2;
  char* wptr;
  unsigned char* ucptr;
  double* dptr;
  double* dptr2;
  double* dptr3;
  uintptr_t sample_idx1;
  uintptr_t sample_idx2;
  uintptr_t clidx1;
  uintptr_t clidx2;
  uintptr_t tcoord;
  uintptr_t line_idx;
  double dxx;
  double dyy;
  uint32_t cov_idx;
  uint32_t slen;
  uint32_t uii;
  int32_t ii;
  char cc;
  retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_sample_ct, sample_exclude, unfiltered_sample_ct - sample_ct, sample_ids, max_sample_id_len, 0, 1, strcmp_deref);
  if (retval) {
    goto cluster_enforce_match_ret_1;
  }
  
  wkspace_mark2 = wkspace_base;
  tbuf[MAXLINELEN - 1] = ' ';
  if (cp->match_fname) {
    sample_idx_to_match_str = (char**)wkspace_alloc(sample_ct * sizeof(intptr_t));
    if (!sample_idx_to_match_str) {
      goto cluster_enforce_match_ret_NOMEM;
    }
    for (sample_idx1 = 0; sample_idx1 < sample_ct; sample_idx1++) {
      sample_idx_to_match_str[sample_idx1] = NULL;
    }
    cov_type_arr = wkspace_base;
    if (wkspace_left < MAXLINELEN) {
      goto cluster_enforce_match_ret_NOMEM;
    }
    if (cp->match_missing_str) {
      missing_str = cp->match_missing_str;
      missing_len = strlen(missing_str);
    }
    if (cp->match_type_fname) {
      if (fopen_checked(&typefile, cp->match_type_fname, "r")) {
	goto cluster_enforce_match_ret_OPEN_FAIL;
      }
      cov_idx = 0;
      line_idx = 0;
      while (fgets(tbuf, MAXLINELEN, typefile)) {
	line_idx++;
        if (!tbuf[MAXLINELEN - 1]) {
	  sprintf(logbuf, "Error: Line %" PRIuPTR " of --match-type file is pathologically long.\n", line_idx);
          goto cluster_enforce_match_ret_INVALID_FORMAT_2;
	}
        bufptr = skip_initial_spaces(tbuf);
	cc = *bufptr;
        while (!is_eoln_kns(cc)) {
	  slen = strlen_se(bufptr);
	  cc = bufptr[0];
	  if ((slen == 1) && ((cc == '0') || (cc == '-'))) {
            cov_type_arr[cov_ct] = 1;
	  } else if ((slen == 1) && ((cc == '1') || (cc == '+'))) {
	    cov_type_arr[cov_ct] = 2;
	  } else if (((slen == 2) && (cc == '-') && (bufptr[1] == '1')) || ((slen == 1) && (cc == '*'))) {
            cov_type_arr[cov_ct] = 0;
	    cov_idx++;
	  } else {
            sprintf(logbuf, "Error: Line %" PRIuPTR " of --match-type file has an invalid token\n(0/1/-1/-/+/* expected).\n", line_idx);
	    goto cluster_enforce_match_ret_INVALID_FORMAT_2;
	  }
	  cov_ct++;
	  bufptr = skip_initial_spaces(&(bufptr[slen]));
	  cc = *bufptr;
	}
	if (cov_ct > 65536) {
          logprint("Error: Too many tokens in --match-type file (max 65536).\n");
	  goto cluster_enforce_match_ret_INVALID_FORMAT;
	}
      }
      if (!feof(typefile)) {
	goto cluster_enforce_match_ret_READ_FAIL;
      }
      fclose_null(&typefile);
      if (!cov_ct) {
        logprint("Error: Empty --match-type file.\n");
	goto cluster_enforce_match_ret_INVALID_FORMAT;
      }
      non_null_cov_ct = cov_ct - cov_idx;
      if (!non_null_cov_ct) {
	logprint("Error: Degenerate --match-type file (all -1/*).\n");
	goto cluster_enforce_match_ret_INVALID_FORMAT;
      }
      while (!cov_type_arr[cov_ct - 1]) {
	cov_ct--;
      }
      wkspace_alloc(cov_ct * sizeof(char)); // cov_type_arr
    }
    retval = open_and_load_to_first_token(&matchfile, cp->match_fname, MAXLINELEN, '\0', "--match file", tbuf, &bufptr, &line_idx);
    if (retval) {
      goto cluster_enforce_match_ret_1;
    }
    if (!cov_ct) {
      cov_ct = count_tokens(bufptr);
      if (cov_ct < 3) {
	line_idx = 1;
	goto cluster_enforce_match_ret_MISSING_TOKENS;
      }
      cov_ct -= 2;
      wkspace_alloc(cov_ct * sizeof(char)); // cov_type_arr
      memset(cov_type_arr, 2, cov_ct);
      non_null_cov_ct = cov_ct;
    }
    wptr = (char*)wkspace_base;
    do {
      line_idx++;
      if (!tbuf[MAXLINELEN - 1]) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --match file is pathologically long.\n", line_idx);
	goto cluster_enforce_match_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      if (bsearch_read_fam_indiv(id_buf, sorted_ids, max_sample_id_len, sample_ct, bufptr, &bufptr2, &ii)) {
	goto cluster_enforce_match_ret_MISSING_TOKENS;
      }
      if (ii == -1) {
	continue;
      }
      sample_idx1 = id_map[(uint32_t)ii];
      if (sample_idx_to_match_str[sample_idx1]) {
	*strchr(id_buf, '\t') = ' ';
        LOGPREPRINTFWW("Error: ID '%s' appears multiple times in --match file.\n", id_buf);
	goto cluster_enforce_match_ret_INVALID_FORMAT_2;
      }
      sample_idx_to_match_str[sample_idx1] = wptr;
      for (cov_idx = 0; cov_idx < cov_ct; cov_idx++) {
        bufptr = skip_initial_spaces(bufptr2);
	if (is_eoln_kns(*bufptr)) {
          goto cluster_enforce_match_ret_MISSING_TOKENS;
	}
        bufptr2 = token_endnn(bufptr);
	if (cov_type_arr[cov_idx]) {
	  uii = bufptr2 - bufptr;
	  if ((uii == missing_len) && (!memcmp(bufptr, missing_str, uii))) {
	    *wptr++ = '\0';
	  } else {
            wptr = memcpyax(wptr, bufptr, uii, '\0');
	  }
	}
      }
      if ((uintptr_t)(((unsigned char*)wptr) - wkspace_base) > wkspace_left - MAXLINELEN) {
	goto cluster_enforce_match_ret_NOMEM;
      }
    } while (fgets(tbuf, MAXLINELEN, matchfile));
    if (!feof(matchfile)) {
      goto cluster_enforce_match_ret_READ_FAIL;
    }
    fclose_null(&matchfile);
    if (non_null_cov_ct != cov_ct) {
      // collapse type array, now that ignored items have been removed
      bufptr = (char*)cov_type_arr;
      bufptr2 = &(bufptr[non_null_cov_ct]);
      while (*bufptr) {
	bufptr++;
      }
      wptr = bufptr;
      while (wptr < bufptr2) {
	do {
          cc = *(++bufptr);
	} while (!cc);
	*wptr++ = cc;
      }
    }
    clidx1 = 0;
    clidx2 = 1; // need these unequal if sample_to_cluster undefined
    for (sample_idx1 = 1; sample_idx1 < sample_ct; sample_idx1++) {
      wptr = sample_idx_to_match_str[sample_idx1];
      if (!wptr) {
	continue;
      }
      if (sample_to_cluster) {
        clidx1 = sample_to_cluster[sample_idx1];
	tcoord = (clidx1 * (clidx1 - 1)) >> 1;
      } else {
	tcoord = (sample_idx1 * (sample_idx1 - 1)) >> 1;
      }
      for (sample_idx2 = 0; sample_idx2 < sample_idx1; sample_idx2++) {
        bufptr = sample_idx_to_match_str[sample_idx2];
	if (!bufptr) {
	  continue;
	}
	if (sample_to_cluster) {
	  clidx2 = sample_to_cluster[sample_idx2];
	  if (clidx2 == clidx1) {
            if (cluster_mismatch_warning) {
	      continue;
	    }
	  } else {
            if (clidx2 < clidx1) {
	      cur_coord = tcoord + clidx2;
	    } else {
	      cur_coord = ((clidx2 * (clidx2 - 1)) >> 1) + clidx1;
	    }
            if (IS_SET(merge_prevented, cur_coord)) {
	      continue;
	    }
	  }
	} else {
	  cur_coord = tcoord + sample_idx2;
	  if (IS_SET(merge_prevented, cur_coord)) {
	    continue;
	  }
	}
	ucptr = cov_type_arr;
        bufptr2 = wptr;
	for (cov_idx = 0; cov_idx < non_null_cov_ct; cov_idx++) {
	  slen = strlen(bufptr);
	  uii = strlen(bufptr2);
	  if (slen && uii) {
	    if (*ucptr++ == 1) {
	      // negative match
	      if ((slen == uii) && (!memcmp(bufptr, bufptr2, uii))) {
		break;
	      }
	    } else {
	      // positive match
	      if ((slen != uii) || memcmp(bufptr, bufptr2, uii)) {
		break;
	      }
	    }
	  }
	  bufptr = &(bufptr[slen + 1]);
	  bufptr2 = &(bufptr2[uii + 1]);
	}
	if (cov_idx < non_null_cov_ct) {
	  if (clidx1 != clidx2) {
	    SET_BIT(merge_prevented, cur_coord);
	  } else {
	    cluster_mismatch_warning = 1;
	  }
	}
      }
    }
    if (cluster_mismatch_warning) {
      logprint("Warning: Initial cluster assignment violates --match constraint.\n");
      cluster_mismatch_warning = 0;
    }
    cov_ct = 0;
    non_null_cov_ct = 0;
    wkspace_reset(wkspace_mark2);
  }
  if (cp->qmatch_fname) {
    sample_idx_to_dvals = (double**)wkspace_alloc(sample_ct * sizeof(intptr_t));
    if (!sample_idx_to_dvals) {
      goto cluster_enforce_match_ret_NOMEM;
    }
    for (sample_idx1 = 0; sample_idx1 < sample_ct; sample_idx1++) {
      sample_idx_to_dvals[sample_idx1] = NULL;
    }
    tol_arr = (double*)wkspace_base;
    if (wkspace_left <= MAXLINELEN * 4) {
      goto cluster_enforce_match_ret_NOMEM;
    }
    if (cp->qmatch_missing_str) {
      missing_str = cp->qmatch_missing_str;
    } else {
      bufptr = intbuf;
      if (missing_pheno < 0) {
	*bufptr++ = '-';
	missing_pheno = -missing_pheno;
      }
      bufptr = uint32_write(bufptr, (uint32_t)missing_pheno);
      *bufptr = '\0';
      missing_str = intbuf;
    }
    missing_len = strlen(missing_str);
    if (fopen_checked(&typefile, cp->qt_fname, "r")) {
      goto cluster_enforce_match_ret_OPEN_FAIL;
    }
    line_idx = 0;
    while (fgets(tbuf, MAXLINELEN, typefile)) {
      line_idx++;
      if (!tbuf[MAXLINELEN - 1]) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --qt file is pathologically long.\n", line_idx);
	goto cluster_enforce_match_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(tbuf);
      while (!is_eoln_kns(*bufptr)) {
        if (scan_double(bufptr, &dxx)) {
	  sprintf(logbuf, "Error: Line %" PRIuPTR " of --qt file has a non-numeric value.\n", line_idx);
	  goto cluster_enforce_match_ret_INVALID_FORMAT_2;
	}
	if (dxx < 0) {
	  if (dxx != -1) {
	    sprintf(logbuf, "Error: Line %" PRIuPTR " of --qt file has an invalid tolerance (-1 = ignore,\nother values must be nonnegative).\n", line_idx);
            goto cluster_enforce_match_ret_INVALID_FORMAT_2;
	  }
	} else {
	  non_null_cov_ct++;
	}
	tol_arr[cov_ct++] = dxx;
	bufptr = skip_initial_spaces(token_endnn(bufptr));
	if (cov_ct > 65536) {
          logprint("Error: Too many values in --qt file (max 65536).\n");
	  goto cluster_enforce_match_ret_INVALID_FORMAT;
	}
      }
    }
    if (!feof(typefile)) {
      goto cluster_enforce_match_ret_READ_FAIL;
    }
    fclose_null(&typefile);
    if (!cov_ct) {
      logprint("Error: Empty --qt file.\n");
      goto cluster_enforce_match_ret_INVALID_FORMAT;
    }
    wkspace_alloc(cov_ct * sizeof(double)); // tol_arr
    if (wkspace_left < non_null_cov_ct * sizeof(double)) {
      goto cluster_enforce_match_ret_NOMEM;
    }
    dptr = (double*)wkspace_base;
    if (fopen_checked(&matchfile, cp->qmatch_fname, "r")) {
      goto cluster_enforce_match_ret_OPEN_FAIL;
    }
    line_idx = 0;
    while (fgets(tbuf, MAXLINELEN, matchfile)) {
      line_idx++;
      if (!tbuf[MAXLINELEN - 1]) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --qmatch file is pathologically long.\n", line_idx);
	goto cluster_enforce_match_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      if (bsearch_read_fam_indiv(id_buf, sorted_ids, max_sample_id_len, sample_ct, bufptr, &bufptr, &ii)) {
        goto cluster_enforce_match_ret_MISSING_TOKENS_Q;
      }
      if (ii == -1) {
	continue;
      }
      sample_idx1 = id_map[(uint32_t)ii];
      if (sample_idx_to_dvals[sample_idx1]) {
	*strchr(id_buf, '\t') = ' ';
        LOGPREPRINTFWW("Error: ID '%s' appears multiple times in --qmatch file.\n", id_buf);
        goto cluster_enforce_match_ret_INVALID_FORMAT_2;
      }
      sample_idx_to_dvals[sample_idx1] = dptr;
      for (cov_idx = 0; cov_idx < cov_ct; cov_idx++) {
	bufptr = skip_initial_spaces(bufptr);
	if (is_eoln_kns(*bufptr)) {
	  goto cluster_enforce_match_ret_MISSING_TOKENS_Q;
	}
        if (tol_arr[cov_idx] != -1) {
	  if ((!memcmp(bufptr, missing_str, missing_len)) && (((unsigned char)bufptr[missing_len]) <= ' ')) {
	    *dptr++ = -HUGE_DOUBLE;
	  } else {
            if (scan_double(bufptr, dptr++)) {
	      sprintf(logbuf, "Error: Line %" PRIuPTR " of --qmatch file has a non-numeric covariate.\n", line_idx);
	      goto cluster_enforce_match_ret_INVALID_FORMAT_2;
	    }
	  }
	}
	bufptr = token_endnn(bufptr);
      }

      if (wkspace_left < (uintptr_t)(((unsigned char*)(&(dptr[non_null_cov_ct]))) - wkspace_base)) {
	goto cluster_enforce_match_ret_NOMEM;
      }
    }
    if (!feof(matchfile)) {
      goto cluster_enforce_match_ret_READ_FAIL;
    }
    fclose_null(&matchfile);
    if (non_null_cov_ct != cov_ct) {
      dptr = tol_arr;
      dptr3 = &(tol_arr[non_null_cov_ct]);
      while (*dptr != -1) {
	dptr++;
      }
      dptr2 = dptr;
      while (dptr2 < dptr3) {
	do {
          dxx = *(++dptr);
	} while (dxx == -1);
        *dptr2++ = dxx;
      }
    }
    clidx1 = 0;
    clidx2 = 1;
    for (sample_idx1 = 1; sample_idx1 < sample_ct; sample_idx1++) {
      dptr3 = sample_idx_to_dvals[sample_idx1];
      if (!dptr3) {
	continue;
      }
      if (sample_to_cluster) {
	clidx1 = sample_to_cluster[sample_idx1];
	tcoord = (clidx1 * (clidx1 - 1)) >> 1;
      } else {
        tcoord = (sample_idx1 * (sample_idx1 - 1)) >> 1;
      }
      for (sample_idx2 = 0; sample_idx2 < sample_idx1; sample_idx2++) {
        dptr = sample_idx_to_dvals[sample_idx2];
        if (!dptr) {
	  continue;
	}
	if (sample_to_cluster) {
	  clidx2 = sample_to_cluster[sample_idx2];
	  if (clidx2 == clidx1) {
            if (cluster_mismatch_warning) {
	      continue;
	    }
	  } else {
            if (clidx2 < clidx1) {
	      cur_coord = tcoord + clidx2;
	    } else {
	      cur_coord = ((clidx2 * (clidx2 - 1)) >> 1) + clidx1;
	    }
            if (IS_SET(merge_prevented, cur_coord)) {
	      continue;
	    }
	  }
	} else {
	  cur_coord = tcoord + sample_idx2;
	  if (IS_SET(merge_prevented, cur_coord)) {
	    continue;
	  }
	}
        dptr2 = dptr3;
	for (cov_idx = 0; cov_idx < non_null_cov_ct; cov_idx++) {
	  dxx = *dptr++;
	  dyy = *dptr2++;
	  if ((dxx != -HUGE_DOUBLE) && (dyy != -HUGE_DOUBLE) && (tol_arr[cov_idx] < fabs(dxx - dyy))) {
            break;
	  }
	}
	if (cov_idx < non_null_cov_ct) {
	  if (clidx1 != clidx2) {
	    SET_BIT(merge_prevented, cur_coord);
	  } else {
	    cluster_mismatch_warning = 1;
	  }
	}
      }
    }
    if (cluster_mismatch_warning) {
      logprint("Warning: Initial cluster assignment violates --qmatch constraint.\n");
    }
  }
  LOGPRINTF("--%smatch constraints applied.\n", cp->match_fname? (cp->qmatch_fname? "match and q" : "") : "q");
  while (0) {
  cluster_enforce_match_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  cluster_enforce_match_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  cluster_enforce_match_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  cluster_enforce_match_ret_MISSING_TOKENS_Q:
    LOGPRINTF("Error: Line %" PRIuPTR " of --qmatch file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  cluster_enforce_match_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of --match file has fewer tokens than expected.\n", line_idx);
  cluster_enforce_match_ret_INVALID_FORMAT_2:
    logprintb();
  cluster_enforce_match_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 cluster_enforce_match_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(matchfile);
  fclose_cond(typefile);
  return retval;
}

uint32_t cluster_main(uintptr_t cluster_ct, uintptr_t* merge_prevented, uintptr_t list_size, uint32_t* sorted_ibs_indices, uint32_t* cluster_index, uint32_t* cur_cluster_sizes, uint32_t sample_ct, uint32_t* cur_cluster_case_cts, uint32_t case_ct, uint32_t ctrl_ct, uint32_t* cur_cluster_remap, Cluster_info* cp, uintptr_t* ibs_ties, uint32_t* merge_sequence) {
  uint32_t is_old_tiebreaks = cp->modifier & CLUSTER_OLD_TIEBREAKS;
  uint32_t* list_end = &(sorted_ibs_indices[list_size]);
  uint32_t max_merge = cluster_ct - cp->min_ct;
  uint32_t max_size = cp->max_size;
  uint32_t max_cases = cp->max_cases;
  uint32_t max_ctrls = cp->max_ctrls;
  uint32_t size_restriction = (max_size < sample_ct);
  uint32_t case_restriction = (max_cases < case_ct);
  uint32_t ctrl_restriction = (max_ctrls < ctrl_ct);
  uint32_t sccr = size_restriction || case_restriction || ctrl_restriction;
  uint32_t case_ctrl_only = 0;
  uint32_t merge_ct = 0;
  uint32_t cur_size = 0;
  uint32_t cur_cases = 0;
  uint32_t cur_ctrls = 0;
  uint32_t* siptr = sorted_ibs_indices;
  uint32_t* siptr2;
  uint32_t* siptr_best;
  uint32_t* tie_end;
  uint32_t clidx_large;
  uint32_t clidx_small;
  uint32_t tcoord1;
  uint32_t tcoord2;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  if (cp->modifier & CLUSTER_CC) {
    for (clidx_small = 0; clidx_small < cluster_ct; clidx_small++) {
      uii = cur_cluster_case_cts[clidx_small];
      if ((!uii) || (uii == cur_cluster_sizes[clidx_small])) {
	case_ctrl_only++;
      }
    }
  }
  if (is_old_tiebreaks) {
    tie_end = sorted_ibs_indices;
  } else {
    tie_end = list_end;
  }
  do {
    if (case_ctrl_only > 1) {
      while (1) {
	if (siptr == tie_end) {
	  if (siptr == list_end) {
	    goto cluster_main_finished;
	  }
          tie_end = &(sorted_ibs_indices[1 + next_unset_unsafe(ibs_ties, siptr - sorted_ibs_indices)]);
	}
	uii = *siptr++;
	if (uii == 0xffffffffU) {
	  continue;
	}
	clidx_large = cur_cluster_remap[uii >> 16];
	clidx_small = cur_cluster_remap[uii & 65535];
	ujj = cur_cluster_case_cts[clidx_small] + cur_cluster_case_cts[clidx_large];
	if ((clidx_small == clidx_large) || (!ujj) || (ujj == cur_cluster_sizes[clidx_small] + cur_cluster_sizes[clidx_large])) {
	  continue;
	}
	if (clidx_large < clidx_small) {
	  ujj = clidx_small;
	  clidx_small = clidx_large;
          clidx_large = ujj;
	}
	if (is_set(merge_prevented, tri_coord_no_diag_32(clidx_small, clidx_large))) {
	  continue;
	}
	if (is_old_tiebreaks && (siptr != tie_end)) {
	  siptr2 = siptr;
	  siptr_best = NULL;
	  for (siptr2 = siptr; siptr2 < tie_end; siptr2++) {
	    ujj = *siptr2;
	    if (ujj != 0xffffffffU) {
	      tcoord2 = cur_cluster_remap[ujj >> 16];
	      tcoord1 = cur_cluster_remap[ujj & 65535];
	      if (tcoord1 == tcoord2) {
		*siptr2 = 0xffffffffU;
	      } else {
	        ujj = cur_cluster_case_cts[tcoord1] + cur_cluster_case_cts[tcoord2];
		if (ujj && (ujj < (cur_cluster_sizes[tcoord1] + cur_cluster_sizes[tcoord2]))) {
		  if (tcoord2 < tcoord1) {
		    ujj = tcoord1;
		    tcoord1 = tcoord2;
		    tcoord2 = ujj;
		  }
		  if (is_set(merge_prevented, tri_coord_no_diag_32(tcoord1, tcoord2))) {
		    *siptr2 = 0xffffffffU;
		  } else if ((tcoord1 < clidx_small) || ((tcoord1 == clidx_small) && (tcoord2 < clidx_large))) {
		    clidx_large = tcoord2;
		    clidx_small = tcoord1;
	            siptr_best = siptr2;
		  }
	        }
	      }
	    }
	  }
          if (siptr_best) {
	    *siptr_best = uii;
	    tcoord2 = cur_cluster_remap[uii >> 16];
	    tcoord1 = cur_cluster_remap[uii & 65535];
	    if (tcoord2 < tcoord1) {
	      uii = tcoord1;
	      tcoord1 = tcoord2;
	      tcoord2 = uii;
	    }
	    cluster_index[tri_coord_no_diag_32(tcoord1, tcoord2)] = siptr_best - sorted_ibs_indices;
	  }
        }
	break;
      }
      uii = cur_cluster_case_cts[clidx_small];
      if ((!uii) || (uii == cur_cluster_sizes[clidx_small])) {
	case_ctrl_only--;
      }
      uii = cur_cluster_case_cts[clidx_large];
      if ((!uii) || (uii == cur_cluster_sizes[clidx_large])) {
	case_ctrl_only--;
      }
    } else {
      while (1) {
	if (siptr == tie_end) {
	  if (siptr == list_end) {
	    goto cluster_main_finished;
	  }
          tie_end = &(sorted_ibs_indices[1 + next_unset_unsafe(ibs_ties, siptr - sorted_ibs_indices)]);
	}
	uii = *siptr++;
	if (uii == 0xffffffffU) {
	  continue;
	}
	clidx_large = cur_cluster_remap[uii >> 16];
	clidx_small = cur_cluster_remap[uii & 65535];
	if (clidx_large < clidx_small) {
	  ujj = clidx_small;
	  clidx_small = clidx_large;
          clidx_large = ujj;
	}
	if ((clidx_small == clidx_large) || is_set(merge_prevented, tri_coord_no_diag_32(clidx_small, clidx_large))) {
	  continue;
	}
	if (is_old_tiebreaks && (siptr != tie_end)) {
	  siptr2 = siptr;
	  siptr_best = NULL;
	  for (siptr2 = siptr; siptr2 < tie_end; siptr2++) {
	    ujj = *siptr2;
	    if (ujj != 0xffffffffU) {
	      tcoord2 = cur_cluster_remap[ujj >> 16];
	      tcoord1 = cur_cluster_remap[ujj & 65535];
	      if (tcoord2 < tcoord1) {
                ujj = tcoord1;
		tcoord1 = tcoord2;
		tcoord2 = ujj;
	      }
	      if ((tcoord1 == tcoord2) || is_set(merge_prevented, tri_coord_no_diag_32(tcoord1, tcoord2))) {
                *siptr2 = 0xffffffffU;
	      } else if ((tcoord1 < clidx_small) || ((tcoord1 == clidx_small) && (tcoord2 < clidx_large))) {
		clidx_large = tcoord2;
		clidx_small = tcoord1;
	        siptr_best = siptr2;
	      }
	    }
	  }
          if (siptr_best) {
	    *siptr_best = uii;
	    tcoord2 = cur_cluster_remap[uii >> 16];
	    tcoord1 = cur_cluster_remap[uii & 65535];
	    if (tcoord2 < tcoord1) {
	      ujj = tcoord1;
	      tcoord1 = tcoord2;
	      tcoord2 = ujj;
	    }
	    cluster_index[tri_coord_no_diag_32(tcoord1, tcoord2)] = siptr_best - sorted_ibs_indices;
	  }
        }
	break;
      }
    }
    *merge_sequence++ = clidx_small;
    *merge_sequence++ = clidx_large;
    cur_cluster_remap[clidx_large] = clidx_small;
    for (uii = clidx_large + 1; uii < cluster_ct; uii++) {
      // there are faster ways to do this update, but it doesn't really matter
      if (cur_cluster_remap[uii] == clidx_large) {
	cur_cluster_remap[uii] = clidx_small;
      }
    }
    if (cur_cluster_sizes) {
      cur_size = cur_cluster_sizes[clidx_small] + cur_cluster_sizes[clidx_large];
      cur_cluster_sizes[clidx_small] = cur_size;
      if (cur_cluster_case_cts) {
	cur_cases = cur_cluster_case_cts[clidx_small] + cur_cluster_case_cts[clidx_large];
	cur_cluster_case_cts[clidx_small] = cur_cases;
	cur_ctrls = cur_size - cur_cases;
	// convert to upper bounds for pairing candidates
	cur_cases = max_cases - cur_cases;
        cur_ctrls = max_ctrls - cur_ctrls;
      }
      cur_size = max_size - cur_size;
    }
    tcoord1 = (clidx_large * (clidx_large - 1)) / 2;
    tcoord2 = (clidx_small * (clidx_small - 1)) / 2;
    if (!sccr) {
      for (uii = 0; uii < clidx_small; uii++) {
	if ((cur_cluster_remap[uii] == uii) && (!is_set(merge_prevented, tcoord2 + uii))) {
	  if (is_set(merge_prevented, tcoord1 + uii)) {
	    set_bit(merge_prevented, tcoord2 + uii);
	  } else {
	    ujj = cluster_index[tcoord1 + uii];
	    ukk = cluster_index[tcoord2 + uii];
	    if (ujj < ukk) {
	      sorted_ibs_indices[ujj] = 0xffffffffU;
	    } else {
	      sorted_ibs_indices[ukk] = 0xffffffffU;
	      cluster_index[tcoord2 + uii] = ujj;
	    }
	  }
	}
      }
      for (uii = clidx_small + 1; uii < clidx_large; uii++) {
	umm = tri_coord_no_diag_32(clidx_small, uii);
	if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, umm))) {
	  if (is_set(merge_prevented, tcoord1 + uii)) {
	    SET_BIT(merge_prevented, umm);
	  } else {
	    ujj = cluster_index[tcoord1 + uii];
	    ukk = cluster_index[umm];
	    if (ujj < ukk) {
	      sorted_ibs_indices[ujj] = 0xffffffffU;
	    } else {
	      sorted_ibs_indices[ukk] = 0xffffffffU;
	      cluster_index[umm] = ujj;
	    }
	  }
	}
      }
      for (uii = clidx_large + 1; uii < cluster_ct; uii++) {
	umm = tri_coord_no_diag_32(clidx_small, uii);
	if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, umm))) {
	  if (is_set(merge_prevented, tri_coord_no_diag_32(clidx_large, uii))) {
	    SET_BIT(merge_prevented, umm);
	  } else {
	    ujj = cluster_index[tri_coord_no_diag_32(clidx_large, uii)];
	    ukk = cluster_index[umm];
	    if (ujj < ukk) {
	      sorted_ibs_indices[ujj] = 0xffffffffU;
	    } else {
	      sorted_ibs_indices[ukk] = 0xffffffffU;
	      cluster_index[umm] = ujj;
	    }
	  }
	}
      }
    } else {
      for (uii = 0; uii < clidx_small; uii++) {
	if ((cur_cluster_remap[uii] == uii) && (!is_set(merge_prevented, tcoord2 + uii))) {
	  if (is_set(merge_prevented, tcoord1 + uii) || (size_restriction && (cur_cluster_sizes[uii] > cur_size)) || (case_restriction && (cur_cluster_case_cts[uii] > cur_cases)) || (ctrl_restriction && (cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > cur_ctrls))) {
	    set_bit(merge_prevented, tcoord2 + uii);
	  } else {
	    ujj = cluster_index[tcoord1 + uii];
	    ukk = cluster_index[tcoord2 + uii];
	    if (ujj < ukk) {
	      sorted_ibs_indices[ujj] = 0xffffffffU;
	    } else {
	      sorted_ibs_indices[ukk] = 0xffffffffU;
	      cluster_index[tcoord2 + uii] = ujj;
	    }
	  }
	}
      }
      for (uii = clidx_small + 1; uii < clidx_large; uii++) {
	umm = tri_coord_no_diag_32(clidx_small, uii);
	if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, umm))) {
	  if (is_set(merge_prevented, tcoord1 + uii) || (size_restriction && (cur_cluster_sizes[uii] > cur_size)) || (case_restriction && (cur_cluster_case_cts[uii] > cur_cases)) || (ctrl_restriction && (cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > cur_ctrls))) {
	    SET_BIT(merge_prevented, umm);
	  } else {
	    ujj = cluster_index[tcoord1 + uii];
	    ukk = cluster_index[umm];
	    if (ujj < ukk) {
	      sorted_ibs_indices[ujj] = 0xffffffffU;
	    } else {
	      sorted_ibs_indices[ukk] = 0xffffffffU;
	      cluster_index[umm] = ujj;
	    }
	  }
	}
      }
      for (uii = clidx_large + 1; uii < cluster_ct; uii++) {
	umm = tri_coord_no_diag_32(clidx_small, uii);
	if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, umm))) {
	  if (is_set(merge_prevented, tri_coord_no_diag_32(clidx_large, uii)) || (size_restriction && (cur_cluster_sizes[uii] > cur_size)) || (case_restriction && (cur_cluster_case_cts[uii] > cur_cases)) || (ctrl_restriction && (cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > cur_ctrls))) {
	    SET_BIT(merge_prevented, umm);
	  } else {
	    ujj = cluster_index[tri_coord_no_diag_32(clidx_large, uii)];
	    ukk = cluster_index[umm];
	    if (ujj < ukk) {
	      sorted_ibs_indices[ujj] = 0xffffffffU;
	    } else {
	      sorted_ibs_indices[ukk] = 0xffffffffU;
	      cluster_index[umm] = ujj;
	    }
	  }
	}
      }
    }
    merge_ct++;
    if (!(merge_ct % 100)) {
      printf("\rClustering... [%u merges performed]", merge_ct);
      fflush(stdout);
    }
  } while (merge_ct < max_merge);
 cluster_main_finished:
  return merge_ct;
}

// We use the arithmetic-simplifying convention of placing the binary heap root
// at index 1 instead of 0.

// If even better performance is needed, the binary heap can be replaced with a
// page-sensitive B-tree.

// todo: 64-bit versions of these functions

void heap_down(uint32_t cur_pos, uint32_t heap_size, double* heap_vals, uint32_t* val_to_cindices, uint32_t* cindices_to_heap_pos) {
  // Because the very first heap operation is always a heap_remove() which sets
  // the last element to zero, we can assume the second child of an element is
  // initialized if the first element is not past the heap end.
  double cur_val = heap_vals[cur_pos];
  uint32_t cur_cindices = val_to_cindices[cur_pos];
  uint32_t child_pos = cur_pos * 2;
  double tmp_val;
  uint32_t tmp_cindices;
  while (child_pos < heap_size) {
    tmp_val = heap_vals[child_pos];
    if (heap_vals[child_pos + 1] > tmp_val) {
      tmp_val = heap_vals[++child_pos];
    }
    if (cur_val >= tmp_val) {
      break;
    }
    tmp_cindices = val_to_cindices[child_pos];
    heap_vals[cur_pos] = tmp_val;
    val_to_cindices[cur_pos] = tmp_cindices;
    cindices_to_heap_pos[tri_coord_no_diag_32(tmp_cindices & 65535, tmp_cindices >> 16)] = cur_pos;
    cur_pos = child_pos;
    child_pos *= 2;
  }
  heap_vals[cur_pos] = cur_val;
  val_to_cindices[cur_pos] = cur_cindices;
  cindices_to_heap_pos[tri_coord_no_diag_32(cur_cindices & 65535, cur_cindices >> 16)] = cur_pos;
}

void heap_up_then_down(uint32_t orig_pos, uint32_t heap_size, double* heap_vals, uint32_t* val_to_cindices, uint32_t* cindices_to_heap_pos) {
  // * to update a value, set heap_vals[orig_pos] and then call this function
  // * to insert, set heap_vals[heap_size] and val_to_cindices[heap_size], then
  //   increment heap_size, then call this
  uint32_t cur_pos = orig_pos;
  double cur_val = heap_vals[orig_pos];
  uint32_t cur_cindices = val_to_cindices[orig_pos];
  uint32_t parent_pos = orig_pos / 2;
  double tmp_val;
  uint32_t tmp_cindices;
  while (parent_pos) {
    tmp_val = heap_vals[parent_pos];
    if (cur_val < tmp_val) {
      break;
    }
    tmp_cindices = val_to_cindices[parent_pos];
    heap_vals[cur_pos] = tmp_val;
    val_to_cindices[cur_pos] = tmp_cindices;
    cindices_to_heap_pos[tri_coord_no_diag_32(tmp_cindices & 65535, tmp_cindices >> 16)] = cur_pos;
    cur_pos = parent_pos;
    parent_pos /= 2;
  }
  if (cur_pos != orig_pos) {
    heap_vals[cur_pos] = cur_val;
    val_to_cindices[cur_pos] = cur_cindices;
    cindices_to_heap_pos[tri_coord_no_diag_32(cur_cindices & 65535, cur_cindices >> 16)] = cur_pos;
  }
  heap_down(cur_pos, heap_size, heap_vals, val_to_cindices, cindices_to_heap_pos);
}

void heap_remove(uint32_t remove_pos, uint32_t* heap_size_ptr, double* heap_vals, uint32_t* val_to_cindices, uint32_t* cindices_to_heap_pos) {
  // 1. replace root with last element
  // 2. shift down until no smaller than children
  uint32_t heap_size = *heap_size_ptr - 1;
  double last_val = heap_vals[heap_size];
  uint32_t last_cindices = val_to_cindices[remove_pos];
  cindices_to_heap_pos[tri_coord_no_diag_32(last_cindices & 65535, last_cindices >> 16)] = 0;
  last_cindices = val_to_cindices[heap_size];
  heap_vals[heap_size] = 0.0;
  heap_vals[remove_pos] = last_val;
  val_to_cindices[remove_pos] = last_cindices;
  cindices_to_heap_pos[tri_coord_no_diag_32(last_cindices & 65535, last_cindices >> 16)] = remove_pos;
  *heap_size_ptr = heap_size;
  heap_up_then_down(remove_pos, heap_size, heap_vals, val_to_cindices, cindices_to_heap_pos);
}

void heap_merge_two(uint32_t coord_aux, uint32_t coord_main, double dsize_aux, double dsize_main, double dsize_recip, uint32_t* heap_size_ptr, double* heap_vals, uint32_t* val_to_cindices, uint32_t* cluster_index) {
  uint32_t heap_pos = cluster_index[coord_aux];
  double cur_dist = dsize_aux * heap_vals[heap_pos];
  heap_remove(heap_pos, heap_size_ptr, heap_vals, val_to_cindices, cluster_index);
  heap_pos = cluster_index[coord_main];
  heap_vals[heap_pos] = (dsize_main * heap_vals[heap_pos] + cur_dist) * dsize_recip;
  heap_up_then_down(heap_pos, *heap_size_ptr, heap_vals, val_to_cindices, cluster_index);
}

void heap_merge_two_cc(uint32_t coord_aux, uint32_t coord_main, double dsize_aux, double dsize_main, double dsize_recip, uint32_t* heap_size_ptr, double* heap_vals, uint32_t* val_to_cindices, uint32_t* cluster_index) {
  uint32_t heap_pos = cluster_index[coord_aux];
  uint32_t heap_pos2 = cluster_index[coord_main];
  double cur_dist = ((dsize_aux * heap_vals[heap_pos]) + (dsize_main * heap_vals[heap_pos2])) * dsize_recip;
  uint32_t tmp_cindices;
  if (heap_pos >= *heap_size_ptr) {
    if (heap_pos2 >= *heap_size_ptr) {
      // insert
      tmp_cindices = val_to_cindices[heap_pos2];
      heap_pos2 = *heap_size_ptr;
      *heap_size_ptr += 1;
      val_to_cindices[heap_pos2] = tmp_cindices;
      cluster_index[tri_coord_no_diag_32(tmp_cindices & 65535, tmp_cindices >> 16)] = heap_pos2;
    }
  } else if (heap_pos2 >= *heap_size_ptr) {
    tmp_cindices = val_to_cindices[heap_pos2];
    heap_pos2 = heap_pos;
    val_to_cindices[heap_pos] = tmp_cindices;
    cluster_index[tri_coord_no_diag_32(tmp_cindices & 65535, tmp_cindices >> 16)] = heap_pos;
  } else {
    heap_remove(heap_pos, heap_size_ptr, heap_vals, val_to_cindices, cluster_index);
  }
  heap_vals[heap_pos2] = cur_dist;
  heap_up_then_down(heap_pos2, *heap_size_ptr, heap_vals, val_to_cindices, cluster_index);
}

uint32_t cluster_group_avg_main(uint32_t cluster_ct, uintptr_t* merge_prevented, uint32_t heap_size, double* heap_vals, uint32_t* val_to_cindices, uint32_t* cluster_index, uint32_t* cur_cluster_sizes, uint32_t sample_ct, uint32_t* cur_cluster_case_cts, uint32_t case_ct, uint32_t ctrl_ct, uint32_t* cur_cluster_remap, Cluster_info* cp, uint32_t* merge_sequence) {
  uint32_t max_merge = cluster_ct - cp->min_ct;
  uint32_t max_size = cp->max_size;
  uint32_t max_cases = cp->max_cases;
  uint32_t max_ctrls = cp->max_ctrls;
  uint32_t size_restriction = (max_size < sample_ct);
  uint32_t case_restriction = (max_cases < case_ct);
  uint32_t ctrl_restriction = (max_ctrls < ctrl_ct);
  uint32_t sccr = size_restriction || case_restriction || ctrl_restriction;
  uint32_t cluster_cc = 0;
  uint32_t top_index = heap_size - 1;
  uint32_t case_ctrl_only = 0;
  uint32_t merge_ct = 0;
  uint32_t cur_cases = 0;
  uint32_t cur_ctrls = 0;
  uint32_t cur_size;
  double dsize1;
  double dsize2;
  double dsize_recip;
  uint32_t clidx_large;
  uint32_t clidx_small;
  uint32_t tcoord1;
  uint32_t tcoord2;
  uint32_t uii;
  uint32_t ujj;
  if (cp->modifier & CLUSTER_CC) {
    for (clidx_small = 0; clidx_small < cluster_ct; clidx_small++) {
      uii = cur_cluster_case_cts[clidx_small];
      if ((!uii) || (uii == cur_cluster_sizes[clidx_small])) {
	case_ctrl_only++;
      }
    }
  }
  if (case_ctrl_only > 1) {
    cluster_cc = 1;
  }
  do {
    if (case_ctrl_only > 1) {
      while (1) {
	if (heap_size == 1) {
	  goto cluster_group_avg_main_finished;
	}
	dsize1 = heap_vals[1];
	uii = val_to_cindices[1];
	heap_remove(1, &heap_size, heap_vals, val_to_cindices, cluster_index);
        clidx_large = cur_cluster_remap[uii >> 16];
        clidx_small = cur_cluster_remap[uii & 65535];
        if (clidx_large < clidx_small) {
          ujj = clidx_small;
          clidx_small = clidx_large;
          clidx_large = ujj;
	}
	if ((clidx_small == clidx_large) || is_set(merge_prevented, tri_coord_no_diag_32(clidx_small, clidx_large))) {
	  continue;
	}
	ujj = cur_cluster_case_cts[clidx_small] + cur_cluster_case_cts[clidx_large];
	if ((!ujj) || (ujj == cur_cluster_sizes[clidx_small] + cur_cluster_sizes[clidx_large])) {
	  // Can't lose track of average distance in this case.  So we save it
	  // past the end of heap_vals[]/val_to_cindices[], and then use a
	  // modified cluster merger which handles these values properly.
	  heap_vals[top_index] = dsize1;
          val_to_cindices[top_index] = uii;
          cluster_index[tri_coord_no_diag_32(clidx_small, clidx_large)] = top_index--;
	  continue;
	}
        break;
      }
    } else {
      while (1) {
        if (heap_size == 1) {
	  goto cluster_group_avg_main_finished;
	}
	uii = val_to_cindices[1];
	heap_remove(1, &heap_size, heap_vals, val_to_cindices, cluster_index);
        clidx_large = cur_cluster_remap[uii >> 16];
        clidx_small = cur_cluster_remap[uii & 65535];
        if (clidx_large < clidx_small) {
          uii = clidx_small;
          clidx_small = clidx_large;
          clidx_large = uii;
	}
	if ((clidx_small == clidx_large) || is_set(merge_prevented, tri_coord_no_diag_32(clidx_small, clidx_large))) {
	  continue;
	}
        break;
      }
    }
    *merge_sequence++ = clidx_small;
    *merge_sequence++ = clidx_large;
    cur_cluster_remap[clidx_large] = clidx_small;
    for (uii = clidx_large + 1; uii < cluster_ct; uii++) {
      if (cur_cluster_remap[uii] == clidx_large) {
        cur_cluster_remap[uii] = clidx_small;
      }
    }
    cur_size = cur_cluster_sizes[clidx_small];
    dsize1 = (double)((int32_t)cur_size);
    uii = cur_cluster_sizes[clidx_large];
    dsize2 = (double)((int32_t)uii);
    cur_size += uii;
    cur_cluster_sizes[clidx_small] = cur_size;
    dsize_recip = 1.0 / ((double)((int32_t)cur_size));
    if (cur_cluster_case_cts) {
      cur_cases = cur_cluster_case_cts[clidx_small] + cur_cluster_case_cts[clidx_large];
      cur_cluster_case_cts[clidx_small] = cur_cases;
      cur_ctrls = cur_size - cur_cases;
      cur_cases = max_cases - cur_cases;
      cur_ctrls = max_ctrls - cur_ctrls;
    }
    if (size_restriction) {
      cur_size = max_size - cur_size;
    }
    tcoord1 = (clidx_large * (clidx_large - 1)) / 2;
    tcoord2 = (clidx_small * (clidx_small - 1)) / 2;
    if (!cluster_cc) {
      if (!sccr) {
	for (uii = 0; uii < clidx_small; uii++) {
	  if ((cur_cluster_remap[uii] == uii) && (!is_set(merge_prevented, tcoord2 + uii))) {
	    if (is_set(merge_prevented, tcoord1 + uii)) {
	      set_bit(merge_prevented, tcoord2 + uii);
	    } else {
	      heap_merge_two(tcoord1 + uii, tcoord2 + uii, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
	for (uii = clidx_small + 1; uii < clidx_large; uii++) {
	  ujj = tri_coord_no_diag_32(clidx_small, uii);
	  if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, ujj))) {
	    if (is_set(merge_prevented, tcoord1 + uii)) {
	      SET_BIT(merge_prevented, ujj);
	    } else {
	      heap_merge_two(tcoord1 + uii, ujj, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
	for (uii = clidx_large + 1; uii < cluster_ct; uii++) {
	  ujj = tri_coord_no_diag_32(clidx_small, uii);
	  if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, ujj))) {
	    if (is_set(merge_prevented, tri_coord_no_diag_32(clidx_large, uii))) {
	      SET_BIT(merge_prevented, ujj);
	    } else {
	      heap_merge_two(tri_coord_no_diag_32(clidx_large, uii), ujj, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
      } else {
	for (uii = 0; uii < clidx_small; uii++) {
	  if ((cur_cluster_remap[uii] == uii) && (!is_set(merge_prevented, tcoord2 + uii))) {
	    if (is_set(merge_prevented, tcoord1 + uii) || (size_restriction && (cur_cluster_sizes[uii] > cur_size)) || (case_restriction && (cur_cluster_case_cts[uii] > cur_cases)) || (ctrl_restriction && (cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > cur_ctrls))) {
	      set_bit(merge_prevented, tcoord2 + uii);
	    } else {
	      heap_merge_two(tcoord1 + uii, tcoord2 + uii, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
	for (uii = clidx_small + 1; uii < clidx_large; uii++) {
	  ujj = tri_coord_no_diag_32(clidx_small, uii);
	  if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, ujj))) {
	    if (is_set(merge_prevented, tcoord1 + uii) || (size_restriction && (cur_cluster_sizes[uii] > cur_size)) || (case_restriction && (cur_cluster_case_cts[uii] > cur_cases)) || (ctrl_restriction && (cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > cur_ctrls))) {
	      SET_BIT(merge_prevented, ujj);
	    } else {
	      heap_merge_two(tcoord1 + uii, ujj, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
	for (uii = clidx_large + 1; uii < cluster_ct; uii++) {
	  ujj = tri_coord_no_diag_32(clidx_small, uii);
	  if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, ujj))) {
	    if (is_set(merge_prevented, tri_coord_no_diag_32(clidx_large, uii)) || (size_restriction && (cur_cluster_sizes[uii] > cur_size)) || (case_restriction && (cur_cluster_case_cts[uii] > cur_cases)) || (ctrl_restriction && (cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > cur_ctrls))) {
	      SET_BIT(merge_prevented, ujj);
	    } else {
	      heap_merge_two(tri_coord_no_diag_32(clidx_large, uii), ujj, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
      }
    } else {
      if (!sccr) {
	for (uii = 0; uii < clidx_small; uii++) {
	  if ((cur_cluster_remap[uii] == uii) && (!is_set(merge_prevented, tcoord2 + uii))) {
	    if (is_set(merge_prevented, tcoord1 + uii)) {
	      set_bit(merge_prevented, tcoord2 + uii);
	    } else {
	      heap_merge_two_cc(tcoord1 + uii, tcoord2 + uii, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
	for (uii = clidx_small + 1; uii < clidx_large; uii++) {
	  ujj = tri_coord_no_diag_32(clidx_small, uii);
	  if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, ujj))) {
	    if (is_set(merge_prevented, tcoord1 + uii)) {
	      SET_BIT(merge_prevented, ujj);
	    } else {
	      heap_merge_two_cc(tcoord1 + uii, ujj, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
	for (uii = clidx_large + 1; uii < cluster_ct; uii++) {
	  ujj = tri_coord_no_diag_32(clidx_small, uii);
	  if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, ujj))) {
	    if (is_set(merge_prevented, tri_coord_no_diag_32(clidx_large, uii))) {
	      SET_BIT(merge_prevented, ujj);
	    } else {
	      heap_merge_two_cc(tri_coord_no_diag_32(clidx_large, uii), ujj, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
      } else {
	for (uii = 0; uii < clidx_small; uii++) {
	  if ((cur_cluster_remap[uii] == uii) && (!is_set(merge_prevented, tcoord2 + uii))) {
	    if (is_set(merge_prevented, tcoord1 + uii) || (size_restriction && (cur_cluster_sizes[uii] > cur_size)) || (case_restriction && (cur_cluster_case_cts[uii] > cur_cases)) || (ctrl_restriction && (cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > cur_ctrls))) {
	      set_bit(merge_prevented, tcoord2 + uii);
	    } else {
	      heap_merge_two_cc(tcoord1 + uii, tcoord2 + uii, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
	for (uii = clidx_small + 1; uii < clidx_large; uii++) {
	  ujj = tri_coord_no_diag_32(clidx_small, uii);
	  if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, ujj))) {
	    if (is_set(merge_prevented, tcoord1 + uii) || (size_restriction && (cur_cluster_sizes[uii] > cur_size)) || (case_restriction && (cur_cluster_case_cts[uii] > cur_cases)) || (ctrl_restriction && (cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > cur_ctrls))) {
	      SET_BIT(merge_prevented, ujj);
	    } else {
	      heap_merge_two_cc(tcoord1 + uii, ujj, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
	for (uii = clidx_large + 1; uii < cluster_ct; uii++) {
	  ujj = tri_coord_no_diag_32(clidx_small, uii);
	  if ((cur_cluster_remap[uii] == uii) && (!IS_SET(merge_prevented, ujj))) {
	    if (is_set(merge_prevented, tri_coord_no_diag_32(clidx_large, uii)) || (size_restriction && (cur_cluster_sizes[uii] > cur_size)) || (case_restriction && (cur_cluster_case_cts[uii] > cur_cases)) || (ctrl_restriction && (cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > cur_ctrls))) {
	      SET_BIT(merge_prevented, ujj);
	    } else {
	      heap_merge_two_cc(tri_coord_no_diag_32(clidx_large, uii), ujj, dsize2, dsize1, dsize_recip, &heap_size, heap_vals, val_to_cindices, cluster_index);
	    }
	  }
	}
      }
    }
    merge_ct++;
    if (!(merge_ct % 100)) {
      printf("\rClustering... [%u merges performed]", merge_ct);
      fflush(stdout);
    }
  } while (merge_ct < max_merge);
 cluster_group_avg_main_finished:
  return merge_ct;
}

void write_cluster1(FILE* outfile, uint32_t clidx, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* pheno_c, uint32_t* sample_idx_to_uidx, uint32_t* merge_sequence, uint32_t merge_ct) {
  // Manually manage recursion, to minimize crash risk when there are e.g. 500k
  // clusters.  Fortunately, no recursion stack is even needed.
  char* sptr;
  char* sptr2;
  uint32_t msidx;
 write_cluster1_recurse:
  putc(' ', outfile);
  sptr = &(sample_ids[sample_idx_to_uidx[clidx] * max_sample_id_len]);
  sptr2 = (char*)memchr(sptr, '\t', max_sample_id_len);
  fwrite(sptr, 1, sptr2 - sptr, outfile);
  putc('_', outfile);
  fputs(&(sptr2[1]), outfile);
  if (pheno_c) {
    if (IS_SET(pheno_c, sample_idx_to_uidx[clidx])) {
      fputs("(2)", outfile);
    } else {
      fputs("(1)", outfile);
    }
  }
  for (msidx = 0; msidx < merge_ct; msidx++) {
    if (merge_sequence[msidx * 2] == clidx) {
      clidx = merge_sequence[msidx * 2 + 1];
      goto write_cluster1_recurse;
    } else if (merge_sequence[msidx * 2 + 1] == clidx) {
      clidx = merge_sequence[msidx * 2];
    }
  }
}

void write_cluster1_oitc(FILE* outfile, uint32_t clidx, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* pheno_c, uint32_t* orig_cluster_map, uint32_t* orig_cluster_starts, uint32_t* late_clidx_to_sample_uidx, uint32_t orig_within_ct, uint32_t* cluster_remap, uint32_t* merge_sequence, uint32_t merge_ct) {
  char* sptr;
  char* sptr2;
  uint32_t uii;
  uint32_t ujj;
 write_cluster1_oitc_recurse:
  putc(' ', outfile);
  if (clidx >= orig_within_ct) {
    sptr = &(sample_ids[late_clidx_to_sample_uidx[clidx - orig_within_ct] * max_sample_id_len]);
    sptr2 = (char*)memchr(sptr, '\t', max_sample_id_len);
    fwrite(sptr, 1, sptr2 - sptr, outfile);
    putc('_', outfile);
    fputs(&(sptr2[1]), outfile);
    if (pheno_c) {
      if (IS_SET(pheno_c, late_clidx_to_sample_uidx[clidx - orig_within_ct])) {
	fputs("(2)", outfile);
      } else {
	fputs("(1)", outfile);
      }
    }
  } else {
    ujj = orig_cluster_starts[clidx + 1];
    for (uii = orig_cluster_starts[clidx]; uii < ujj; uii++) {
      sptr = &(sample_ids[orig_cluster_map[uii] * max_sample_id_len]);
      sptr2 = (char*)memchr(sptr, '\t', max_sample_id_len);
      fwrite(sptr, 1, sptr2 - sptr, outfile);
      putc('_', outfile);
      fputs(&(sptr2[1]), outfile);
      if (pheno_c) {
	if (IS_SET(pheno_c, orig_cluster_map[uii])) {
	  fputs("(2)", outfile);
	} else {
	  fputs("(1)", outfile);
	}
      }
    }
  }
  for (uii = 0; uii < merge_ct; uii++) {
    if (merge_sequence[uii * 2] == clidx) {
      clidx = merge_sequence[uii * 2 + 1];
      goto write_cluster1_oitc_recurse;
    } else if (merge_sequence[uii * 2 + 1] == clidx) {
      clidx = merge_sequence[uii * 2];
    }
  }
}

int32_t write_cluster_solution(char* outname, char* outname_end, uint32_t* orig_sample_to_cluster, uintptr_t sample_ct, uint32_t* orig_cluster_map, uint32_t* orig_cluster_starts, uint32_t* late_clidx_to_sample_uidx, uint32_t orig_within_ct, uint32_t orig_cluster_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* pheno_c, uint32_t* sample_idx_to_uidx, Cluster_info* cp, uint32_t* cluster_remap, uint32_t* clidx_table_space, uint32_t merge_ct, uint32_t* merge_sequence) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uint32_t only2 = cp->modifier & CLUSTER_ONLY2;
  uint32_t report_pheno = (cp->modifier & CLUSTER_CC) || (cp->max_ctrls != 0xffffffffU);
  uint32_t pct = 1;
  int32_t retval = 0;
  char wbuf[16];
  uint32_t* clidx_remap = &(clidx_table_space[(((uintptr_t)orig_cluster_ct) * (orig_cluster_ct - 1)) >> 1]);
  char* sptr;
  char* sptr2;
  char* wptr;
  uint32_t* cur_remap;
  uint32_t sample_idx;
  uint32_t clidx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t msidx;
  clidx = 0;
  for (uii = 0; uii < orig_cluster_ct; uii++) {
    if (cluster_remap[uii] == uii) {
      clidx_remap[uii] = clidx++;
    }
  }
  memcpy(outname_end, ".cluster2", 10);
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_cluster_solution_ret_OPEN_FAIL;
  }
  fputs("Writing cluster solution...", stdout);
  fflush(stdout);
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    if (orig_sample_to_cluster) {
      clidx = cluster_remap[orig_sample_to_cluster[sample_idx]];
    } else {
      clidx = cluster_remap[sample_idx];
    }
    sptr = &(sample_ids[sample_idx_to_uidx[sample_idx] * max_sample_id_len]);
    sptr2 = (char*)memchr(sptr, '\t', max_sample_id_len);
    wptr = memcpyax(tbuf, sptr, (sptr2 - sptr), ' ');
    wptr = strcpyax(wptr, &(sptr2[1]), '\t');
    wptr = uint32_writex(wptr, clidx_remap[clidx], '\n');
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto write_cluster_solution_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto write_cluster_solution_ret_WRITE_FAIL;
  }
  if (!only2) {
    outname_end[8] = '1';
    if (fopen_checked(&outfile, outname, "w")) {
      goto write_cluster_solution_ret_OPEN_FAIL;
    }
    memcpy(tbuf, "SOL-", 4);
    for (clidx = 0; clidx < orig_cluster_ct; clidx++) {
      if (cluster_remap[clidx] == clidx) {
        wptr = uint32_writex(&(tbuf[4]), clidx_remap[clidx], '\t');
        if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	  goto write_cluster_solution_ret_WRITE_FAIL;
	}
        if (!orig_sample_to_cluster) {
	  write_cluster1(outfile, clidx, sample_ids, max_sample_id_len, report_pheno? pheno_c : NULL, sample_idx_to_uidx, merge_sequence, merge_ct);
	} else {
	  write_cluster1_oitc(outfile, clidx, sample_ids, max_sample_id_len, report_pheno? pheno_c : NULL, orig_cluster_map, orig_cluster_starts, late_clidx_to_sample_uidx, orig_within_ct, cluster_remap, merge_sequence, merge_ct);
	}
	if (putc_checked('\n', outfile)) {
	  goto write_cluster_solution_ret_WRITE_FAIL;
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto write_cluster_solution_ret_WRITE_FAIL;
    }
    if (cp->modifier & CLUSTER_MISSING) {
      memcpy(&(outname_end[8]), "3.missing", 10);
    } else {
      outname_end[8] = '3';
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto write_cluster_solution_ret_WRITE_FAIL;
    }
    clidx = 0;
    for (uii = 0; uii < orig_cluster_ct; uii++) {
      if (cluster_remap[uii] == uii) {
	ujj = orig_cluster_ct - (++clidx);
	clidx_remap[uii] = ujj;
	cur_remap = &(clidx_table_space[(((uintptr_t)ujj) * (ujj - 1)) >> 1]);
	ujj = uii;
	for (msidx = 0; msidx < merge_ct; msidx++) {
	  if (merge_sequence[2 * msidx + 1] < uii) {
	    ujj--;
	  }
          *cur_remap++ = ujj;
	}
      }
    }
    uii = merge_ct;
    while (uii) {
      clidx = merge_sequence[2 * uii - 1];
      ujj = --uii;
      clidx_remap[clidx] = ujj;
      cur_remap = &(clidx_table_space[(((uintptr_t)ujj) * (ujj - 1)) >> 1]);
      ujj = clidx;
      for (msidx = 0; msidx < uii; msidx++) {
	if (merge_sequence[2 * msidx + 1] < clidx) {
	  ujj--;
	}
	*cur_remap++ = ujj;
      }
    }
    fputs(" 0%", stdout);
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      sptr = &(sample_ids[sample_idx_to_uidx[sample_idx] * max_sample_id_len]);
      sptr2 = (char*)memchr(sptr, '\t', max_sample_id_len);
      if (fwrite_checked(sptr, sptr2 - sptr, outfile)) {
	goto write_cluster_solution_ret_WRITE_FAIL;
      }
      putc(' ', outfile);
      fputs(&(sptr2[1]), outfile);
      putc('\t', outfile);
      if (orig_sample_to_cluster) {
        clidx = orig_sample_to_cluster[sample_idx];
      } else {
	clidx = sample_idx;
      }
      wptr = uint32_writex(wbuf, clidx, ' ');
      fwrite(wbuf, 1, wptr - wbuf, outfile);
      uii = 0;
      if (merge_ct) {
	ujj = clidx_remap[clidx];
	while (1) {
	  cur_remap = &(clidx_table_space[((((uintptr_t)ujj) * (ujj - 1)) >> 1) + uii]);
	  if (ujj > merge_ct) {
	    ujj = merge_ct;
	  }
	  for (; uii < ujj; uii++) {
	    wptr = uint32_writex(wbuf, *cur_remap++, ' ');
	    fwrite(wbuf, 1, wptr - wbuf, outfile);
	  }
	  if (ujj == merge_ct) {
	    break;
	  }
	  clidx = merge_sequence[ujj * 2];
          ujj = clidx_remap[merge_sequence[ujj * 2]];
	}
      }
      for (ujj = merge_ct + 1; ujj < orig_cluster_ct; ujj++) {
	fputs("0 ", outfile);
      }
      if (putc_checked('\n', outfile)) {
	goto write_cluster_solution_ret_WRITE_FAIL;
      }
      if ((sample_idx + 1) * 100LLU >= ((uint64_t)pct * sample_ct)) {
	if (pct > 10) {
          putchar('\b');
	}
        pct = (((uint64_t)(sample_idx + 1)) * 100) / sample_ct;
        printf("\b\b%u%%", pct++);
	fflush(stdout);
      }
    }
    putc('\n', outfile);
    if (fclose_null(&outfile)) {
      goto write_cluster_solution_ret_WRITE_FAIL;
    }
    *outname_end = '\0';
    putchar('\r');
    LOGPREPRINTFWW("Cluster solution written to %s.cluster1 , %s.cluster2 , and %s.cluster3%s .\n", outname, outname, outname, (cp->modifier & CLUSTER_MISSING)? ".missing" : "");
  } else {
    *outname_end = '\0';
    putchar('\r');
    LOGPREPRINTFWW("Cluster solution written to %s.cluster2 .\n", outname);
  }
  logprintb();
  while (0) {
  write_cluster_solution_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_cluster_solution_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

#ifndef NOLAPACK
int32_t mds_plot(char* outname, char* outname_end, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t* sample_idx_to_uidx, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uint32_t cur_cluster_ct, uint32_t merge_ct, uint32_t* orig_sample_to_cluster, uint32_t* cur_cluster_remap, uint32_t dim_ct, uint32_t is_mds_cluster, uint32_t dump_eigvals, double* dists) {
  FILE* outfile = NULL;
  uintptr_t final_cluster_ct = cur_cluster_ct - merge_ct;
  double grand_mean = 0.0;
  uintptr_t ulii = 0;
  int32_t retval = 0;
  char jobz = 'V';
  char range = 'I';
  char uplo = 'U';
  double nz = 0.0;
  double zz = -1.0;
  __CLPK_integer info = 0;
  __CLPK_integer lwork = -1;
  __CLPK_integer liwork = -1;
  __CLPK_integer mdim;
  __CLPK_integer i1;
  __CLPK_integer i2;
  __CLPK_integer out_m;
  __CLPK_integer ldz;
  __CLPK_integer optim_liwork;
  __CLPK_integer* iwork;
  __CLPK_integer* isuppz;
  double* work;
  double optim_lwork;
  double* main_matrix;
  double* column_means;
  double* out_w;
  double* out_z;
  uint32_t* final_cluster_remap;
  uint32_t* final_cluster_sizes;
  double* dptr;
  double* dptr2;
  char* wptr;
  char* wptr2;
  uintptr_t sample_idx;
  uintptr_t clidx1;
  uintptr_t clidx2;
  uint32_t dim_idx;
  uint32_t uii;
  uint32_t ujj;
  double dxx;
  double dyy;
  double* sqrt_eigvals;
  final_cluster_remap = (uint32_t*)malloc(cur_cluster_ct * sizeof(int32_t));
  if (!final_cluster_remap) {
    goto mds_plot_ret_NOMEM;
  }
  if ((sample_ct > 5000) && (!is_mds_cluster) && (final_cluster_ct < sample_ct) && (final_cluster_ct > 1)) {
    LOGPRINTF("Warning: Per-sample --mds-plot can be very slow with over 5000 %s.\nConsider using the 'by-cluster' modifier.\n", g_species_plural);
  }
  for (clidx1 = 0; clidx1 < cur_cluster_ct; clidx1++) {
    clidx2 = cur_cluster_remap[clidx1];
    if (clidx2 == clidx1) {
      final_cluster_remap[clidx1] = ulii++;
    } else {
      final_cluster_remap[clidx1] = final_cluster_remap[clidx2];
    }
  }
  if (is_mds_cluster) {
    if (wkspace_alloc_d_checked(&main_matrix, final_cluster_ct * final_cluster_ct * sizeof(double)) ||
        wkspace_alloc_ui_checked(&final_cluster_sizes, final_cluster_ct * sizeof(int32_t))) {
      goto mds_plot_ret_NOMEM;
    }
    fill_double_zero(main_matrix, final_cluster_ct * final_cluster_ct);
    fill_uint_zero(final_cluster_sizes, final_cluster_ct);
    dptr = dists;
    final_cluster_sizes[final_cluster_remap[0]] = 1;
    for (uii = 1; uii < cur_cluster_ct; uii++) {
      clidx1 = final_cluster_remap[uii];
      final_cluster_sizes[clidx1] += 1;
      dptr2 = &(main_matrix[clidx1 * final_cluster_ct]);
      for (ujj = 0; ujj < uii; ujj++) {
	clidx2 = final_cluster_remap[ujj];
	if (clidx2 < clidx1) {
	  dptr2[clidx2] += (*dptr);
	} else if (clidx1 > clidx2) {
          main_matrix[clidx2 * final_cluster_ct + clidx1] += (*dptr);
	}
	dptr++;
      }
    }
    for (clidx1 = 1; clidx1 < final_cluster_ct; clidx1++) {
      dptr = &(main_matrix[clidx1 * final_cluster_ct]);
      ulii = final_cluster_sizes[clidx1];
      for (clidx2 = 0; clidx2 < clidx1; clidx2++) {
	dxx = (double)((intptr_t)(ulii * final_cluster_sizes[clidx2]));
	dptr[clidx2] /= dxx;
      }
    }
    ulii = final_cluster_ct;
  } else {
    wkspace_reset(dists);
    if (wkspace_alloc_d_checked(&main_matrix, sample_ct * sample_ct * sizeof(double))) {
      goto mds_plot_ret_NOMEM;
    }
    // expand triangular diagonal-free matrix to square matrix
    ulii = ((sample_ct - 1) * (sample_ct - 2)) >> 1;
    for (sample_idx = sample_ct - 1; sample_idx;) {
      memcpy(&(main_matrix[sample_idx * sample_ct]), &(main_matrix[ulii]), sample_idx * sizeof(double));
      ulii -= (--sample_idx);
    }
    ulii = sample_ct + 1;
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      main_matrix[sample_idx * ulii] = 0.0;
    }
    ulii = sample_ct;
  }
  if (wkspace_alloc_d_checked(&column_means, ulii * sizeof(double))) {
    goto mds_plot_ret_NOMEM;
  }
  fill_double_zero(column_means, ulii);
  // bottom left filled with IBS values.  Now subtract them from 1 and square
  // them, and extract column means...
  for (clidx1 = 1; clidx1 < ulii; clidx1++) {
    dptr = &(main_matrix[clidx1 * ulii]);
    dptr2 = column_means;
    dyy = 0.0;
    for (clidx2 = 0; clidx2 < clidx1; clidx2++) {
      dxx = 1.0 - (*dptr);
      dxx *= dxx;
      *dptr++ = dxx;
      *dptr2 += dxx;
      dptr2++;
      dyy += dxx;
    }
    *dptr2 += dyy;
  }
  dxx = 1.0 / ((double)((intptr_t)ulii));
  grand_mean = 0.0;
  for (clidx1 = 0; clidx1 < ulii; clidx1++) {
    column_means[clidx1] *= dxx;
    grand_mean += column_means[clidx1];
  }
  grand_mean *= dxx;
  // ...then double-center and multiply by -0.5
  for (clidx1 = 1; clidx1 < ulii; clidx1++) {
    dxx = column_means[clidx1];
    dptr = &(main_matrix[clidx1 * ulii]);
    dptr2 = column_means;
    for (clidx2 = 0; clidx2 < clidx1; clidx2++) {
      *dptr = -0.5 * ((*dptr) - dxx - (*dptr2++) + grand_mean);
      dptr++;
    }
  }

  if (dim_ct > ulii) {
    dim_ct = ulii;
  }

  LOGPRINTF("Performing multidimensional scaling analysis (%u dimension%s)...", dim_ct, (dim_ct == 1)? "" : "s");
  fflush(stdout);

  // no need to fill upper right

  // see eigen_lapack() in PLINK 1.07 lapackf.cpp (though we use dsyevr_
  // instead of dsyevx_).  todo: use arpack-ng instead?
  mdim = ulii;
  i2 = mdim;
  i1 = i2 + 1 - dim_ct;
  if (wkspace_alloc_d_checked(&out_w, dim_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&out_z, dim_ct * ulii * sizeof(double))) {
    goto mds_plot_ret_NOMEM;
  }
  fill_double_zero(out_w, dim_ct);
  fill_double_zero(out_z, dim_ct * ulii);
  isuppz = (__CLPK_integer*)wkspace_alloc(2 * dim_ct * sizeof(__CLPK_integer));
  if (!isuppz) {
    goto mds_plot_ret_NOMEM;
  }
  fill_int_zero((int32_t*)isuppz, 2 * dim_ct * (sizeof(__CLPK_integer) / sizeof(int32_t)));
  ldz = mdim;

  dsyevr_(&jobz, &range, &uplo, &mdim, main_matrix, &mdim, &nz, &nz, &i1, &i2, &zz, &out_m, out_w, out_z, &ldz, isuppz, &optim_lwork, &lwork, &optim_liwork, &liwork, &info);
  lwork = (int32_t)optim_lwork;
  if (wkspace_alloc_d_checked(&work, lwork * sizeof(double))) {
    goto mds_plot_ret_NOMEM;
  }
  liwork = optim_liwork;
  iwork = (__CLPK_integer*)wkspace_alloc(liwork * sizeof(__CLPK_integer));
  if (!iwork) {
    goto mds_plot_ret_NOMEM;
  }
  fill_double_zero(work, lwork);
  fill_int_zero((int32_t*)iwork, liwork * (sizeof(__CLPK_integer) / sizeof(int32_t)));
  dsyevr_(&jobz, &range, &uplo, &mdim, main_matrix, &mdim, &nz, &nz, &i1, &i2, &zz, &out_m, out_w, out_z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

  // * out_w[0..(dim_ct-1)] contains eigenvalues
  // * out_z[(ii*ulii)..(ii*ulii + ulii - 1)] is eigenvector corresponding to
  //   out_w[ii]
  wkspace_reset(isuppz);
  if (wkspace_alloc_d_checked(&sqrt_eigvals, dim_ct * sizeof(double))) {
    goto mds_plot_ret_NOMEM;
  }
  for (dim_idx = 0; dim_idx < dim_ct; dim_idx++) {
    dxx = out_w[dim_idx];
    if (dxx > 0.0) {
      sqrt_eigvals[dim_idx] = sqrt(dxx);
    } else {
      sqrt_eigvals[dim_idx] = 0.0;
    }
  }
  // repurpose main_matrix as mds[]
  dptr = main_matrix;
  for (clidx1 = 0; clidx1 < ulii; clidx1++) {
    for (dim_idx = 0; dim_idx < dim_ct; dim_idx++) {
      *dptr++ = out_z[dim_idx * ulii + clidx1] * sqrt_eigvals[dim_idx];
    }
  }
  logprint(" done.\n");

  memcpy(outname_end, ".mds", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto mds_plot_ret_OPEN_FAIL;
  }
  sprintf(tbuf, "%%%us %%%us    SOL ", plink_maxfid, plink_maxiid);
  fprintf(outfile, tbuf, "FID", "IID");
  tbuf[22] = ' ';
  for (dim_idx = 0; dim_idx < dim_ct; dim_idx++) {
    wptr = uint32_write(tbuf, dim_idx + 1);
    uii = wptr - tbuf;
    wptr2 = memseta(&(tbuf[10]), 32, 11 - uii);
    *wptr2++ = 'C';
    memcpy(wptr2, tbuf, uii);
    fwrite(&(tbuf[10]), 1, 13, outfile);
  }
  if (putc_checked('\n', outfile)) {
    goto mds_plot_ret_WRITE_FAIL;
  }
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    wptr2 = &(sample_ids[sample_idx_to_uidx[sample_idx] * max_sample_id_len]);
    uii = strlen_se(wptr2);
    wptr = fw_strcpyn(plink_maxfid, uii, wptr2, tbuf);
    *wptr++ = ' ';
    wptr = fw_strcpy(plink_maxiid, &(wptr2[uii + 1]), wptr);
    *wptr++ = ' ';
    if (orig_sample_to_cluster) {
      uii = orig_sample_to_cluster[sample_idx];
    } else {
      uii = sample_idx;
    }
    uii = final_cluster_remap[uii];
    wptr = uint32_writew6x(wptr, uii, ' ');
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto mds_plot_ret_WRITE_FAIL;
    }
    if (!is_mds_cluster) {
      dptr = &(main_matrix[(sample_idx + 1) * dim_ct]);
      for (dim_idx = 0; dim_idx < dim_ct; dim_idx++) {
        wptr = double_g_writex(&(tbuf[11]), *(--dptr), ' ');
	uii = wptr - (&(tbuf[11]));
	if (uii < 13) {
	  wptr2 = &(wptr[-13]);
	  memset(wptr2, 32, 13 - uii);
	} else {
	  wptr2 = &(tbuf[11]);
	}
	fwrite(wptr2, 1, wptr - wptr2, outfile);
      }
    } else {
      dptr = &(main_matrix[(uii + 1) * dim_ct]);
      for (dim_idx = 0; dim_idx < dim_ct; dim_idx++) {
        wptr = double_g_writex(&(tbuf[11]), *(--dptr), ' ');
	uii = wptr - (&(tbuf[11]));
	if (uii < 13) {
	  wptr2 = &(wptr[-13]);
	  memset(wptr2, 32, 13 - uii);
	} else {
	  wptr2 = &(tbuf[11]);
	}
	fwrite(wptr2, 1, wptr - wptr2, outfile);
      }
    }
    if (putc_checked('\n', outfile)) {
      goto mds_plot_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto mds_plot_ret_WRITE_FAIL;
  }
  if (!dump_eigvals) {
    LOGPREPRINTFWW("MDS solution written to %s .\n", outname);
  } else {
    LOGPREPRINTFWW("MDS solution written to %s (eigenvalues in %s.eigvals ).\n", outname, outname);
    memcpy(&(outname_end[4]), ".eigvals", 9);
    if (fopen_checked(&outfile, outname, "w")) {
      goto mds_plot_ret_OPEN_FAIL;
    }
    for (dim_idx = dim_ct; dim_idx; dim_idx--) {
      wptr = double_g_writex(tbuf, sqrt_eigvals[dim_idx - 1] * sqrt_eigvals[dim_idx - 1], '\n');
      *wptr = '\0';
      fputs(tbuf, outfile);
    }
    if (fclose_null(&outfile)) {
      goto mds_plot_ret_WRITE_FAIL;
    }
  }
  logprintb();
  while (0) {
  mds_plot_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  mds_plot_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  mds_plot_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  free_cond(final_cluster_remap);
  wkspace_reset(dists);
  return retval;
}
#endif
