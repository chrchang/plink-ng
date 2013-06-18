#include "wdist_cluster.h"

void cluster_init(Cluster_info* cluster_ptr) {
  cluster_ptr->fname = NULL;
  cluster_ptr->match_fname = NULL;
  cluster_ptr->match_missing_str = NULL;
  cluster_ptr->match_type_fname = NULL;
  cluster_ptr->qmatch_fname = NULL;
  cluster_ptr->qmatch_missing_str = NULL;
  cluster_ptr->qt_fname = NULL;
  cluster_ptr->modifier = 0;
  cluster_ptr->ppc = 0.0;
  cluster_ptr->max_size = 0xffffffffU;
  cluster_ptr->max_cases = 0xffffffffU;
  cluster_ptr->max_controls = 0xffffffffU;
  cluster_ptr->min_ct = 1;
  cluster_ptr->mds_dim_ct = 0;
  cluster_ptr->min_ibm = 0;
}

int32_t load_clusters(char* fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uint32_t mwithin_col, uint32_t keep_na, uintptr_t* cluster_ct_ptr, uint32_t** cluster_map_ptr, uint32_t** cluster_starts_ptr, char** cluster_ids_ptr, uintptr_t* max_cluster_id_len_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t topsize = 0;
  int32_t retval = 0;
  char* idbuf = &(tbuf[MAXLINELEN]);
  uintptr_t max_cluster_id_len = 0;
  uintptr_t assigned_ct = 0;
  uintptr_t cluster_ct = 0;
  Ll_str* cluster_names = NULL;
  uintptr_t* already_seen;
  char* cluster_ids;
  uint32_t* cluster_map;
  uint32_t* cluster_starts;
  uint32_t* tmp_cluster_starts;
  uintptr_t topsize_bak;
  Ll_str* llptr;
  char* sorted_ids;
  uint32_t* id_map;
  char* fam_id;
  char* indiv_id;
  char* cluster_name_ptr;
  uintptr_t ulii;
  int32_t sorted_idx;
  uint32_t indiv_uidx;
  uint32_t slen;
  uint32_t uii;

  sorted_ids = (char*)top_alloc(&topsize, indiv_ct * max_person_id_len);
  if (!sorted_ids) {
    goto load_clusters_ret_NOMEM;
  }
  id_map = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!id_map) {
    goto load_clusters_ret_NOMEM;
  }
  topsize_bak = topsize;
  already_seen = (uintptr_t*)top_alloc(&topsize, indiv_ctl * sizeof(intptr_t));
  if (!already_seen) {
    goto load_clusters_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, indiv_ctl);
  memcpy(sorted_ids, person_ids, indiv_ct * max_person_id_len);
  wkspace_left -= topsize;
  retval = sort_item_ids_noalloc(sorted_ids, id_map, unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, 0, strcmp_deref);
  wkspace_left += topsize;
  if (retval) {
    goto load_clusters_ret_1;
  }

  // two-pass load
  // 1. load cluster names, track longest length, validate format, verify no
  //    individual ID appears multiple times
  // intermission. sort cluster names, purge duplicates, allocate memory for
  //               return values
  // 2. populate return name arrays
  if (fopen_checked(&infile, fname, "r")) {
    goto load_clusters_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  if (!mwithin_col) {
    mwithin_col = 1;
  }
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in --within file.\n");
      goto load_clusters_ret_INVALID_FORMAT;
    }
    fam_id = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*fam_id)) {
      continue;
    }
    indiv_id = next_item(fam_id);
    cluster_name_ptr = next_item_mult(indiv_id, mwithin_col);
    if (no_more_items_kns(cluster_name_ptr)) {
      logprint("Error: Fewer entries than expected in --within file line.\n");
      goto load_clusters_ret_INVALID_FORMAT;
    }
    sorted_idx = bsearch_fam_indiv(idbuf, sorted_ids, max_person_id_len, indiv_ct, fam_id, indiv_id);
    if (sorted_idx == -1) {
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      idbuf[strlen_se(fam_id)] = ' ';
      sprintf(logbuf, "Error: Duplicate individual %s in --within file.\n", idbuf);
      logprintb();
      goto load_clusters_ret_INVALID_FORMAT;
    }
    set_bit_noct(already_seen, sorted_idx);
    slen = strlen_se(cluster_name_ptr);
    if ((!keep_na) && (slen == 2) && (!memcmp(cluster_name_ptr, "NA", 2))) {
      // postponed to here because, even without 'keep-NA', we do not want to
      // ignore cluster=NA lines for the purpose of detecting duplicate indivs
      continue;
    }
    if (slen >= max_cluster_id_len) {
      max_cluster_id_len = slen + 1;
    }
    llptr = top_alloc_llstr(&topsize, slen + 1);
    llptr->next = cluster_names;
    memcpyx(llptr->ss, cluster_name_ptr, slen, '\0');
    cluster_names = llptr;
    assigned_ct++;
  }
  if (!feof(infile)) {
    goto load_clusters_ret_READ_FAIL;
  }
  if (cluster_names) {
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
    // deallocate cluster ID linked list and duplicate indiv ID detector from
    // top of stack, allocate cluster size tracker
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
    // shrink
    wkspace_reset(cluster_ids);
    wkspace_alloc_c_checked(cluster_ids_ptr, cluster_ct * max_cluster_id_len);
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
      indiv_id = next_item(fam_id);
      cluster_name_ptr = next_item_mult(indiv_id, mwithin_col);
      slen = strlen_se(cluster_name_ptr);
      if ((!keep_na) && (slen == 2) && (!memcmp(cluster_name_ptr, "NA", 2))) {
	continue;
      }
      sorted_idx = bsearch_fam_indiv(idbuf, sorted_ids, max_person_id_len, indiv_ct, fam_id, indiv_id);
      if (sorted_idx == -1) {
	continue;
      }
      indiv_uidx = id_map[(uint32_t)sorted_idx];
      cluster_name_ptr[slen] = '\0';
      sorted_idx = bsearch_str_natural(cluster_name_ptr, cluster_ids, max_cluster_id_len, 0, cluster_ct - 1);
      uii = tmp_cluster_starts[(uint32_t)sorted_idx];
      tmp_cluster_starts[(uint32_t)sorted_idx] += 1;
      cluster_map[uii] = indiv_uidx;
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
    sprintf(logbuf, "--within: %" PRIuPTR " cluster%s loaded, covering a total of %" PRIuPTR " %s.\n", cluster_ct, (cluster_ct == 1)? "" : "s", assigned_ct, species_str(assigned_ct));
    logprintb();
  } else {
    logprint("Warning: No individuals named in --within file remain in the current analysis.\n");
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

void fill_indiv_to_cluster(uint32_t* indiv_to_cluster, uintptr_t unfiltered_indiv_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts) {
  uint32_t* cluster_map_pos = cluster_map;
  uint32_t* cluster_end_ptr;
  uint32_t cluster_idx;
  // 0xffffffffU cluster index = unassigned
  fill_uint_one(indiv_to_cluster, unfiltered_indiv_ct);
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    cluster_end_ptr = &(cluster_map[cluster_starts[cluster_idx + 1]]);
    do {
      indiv_to_cluster[*cluster_map_pos] = cluster_idx;
    } while (++cluster_map_pos < cluster_end_ptr);
  }
}

int32_t write_clusters(char* outname, char* outname_end, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uint32_t omit_unassigned, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t indiv_uidx = 0;
  int32_t retval = 0;
  uint32_t* indiv_to_cluster;
  char* person_id_ptr;
  char* bufptr;
  uintptr_t indiv_idx;
  uint32_t cluster_idx;
  uint32_t slen;
  if (wkspace_alloc_ui_checked(&indiv_to_cluster, unfiltered_indiv_ct * sizeof(int32_t))) {
    goto write_cluster_ret_NOMEM;
  }
  fill_indiv_to_cluster(indiv_to_cluster, unfiltered_indiv_ct, cluster_ct, cluster_map, cluster_starts);
  memcpy(outname_end, ".clst", 6);
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_cluster_ret_OPEN_FAIL;
  }
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    cluster_idx = indiv_to_cluster[indiv_uidx];
    if ((!omit_unassigned) || (cluster_idx != 0xffffffffU)) {
      person_id_ptr = &(person_ids[indiv_uidx * max_person_id_len]);
      slen = strlen_se(person_id_ptr);
      bufptr = memcpyax(tbuf, person_id_ptr, slen, ' ');
      bufptr = strcpyax(bufptr, &(person_id_ptr[slen + 1]), ' ');
      if (cluster_idx != 0xffffffffU) {
        bufptr = strcpyax(bufptr, &(cluster_ids[cluster_idx * max_cluster_id_len]), '\n');
      } else {
        bufptr = memcpyl3a(bufptr, "NA\n");
      }
      if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	goto write_cluster_ret_WRITE_FAIL;
      }
    }
    indiv_uidx++;
  }
  if (fclose_null(&outfile)) {
    goto write_cluster_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "Pruned cluster assignments written to %s.\n", outname);
  logprintb();
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
      case_ct += is_set(pheno_c, cluster_map[map_idx]);
    }
    cluster_case_cts[cluster_idx] = case_ct;
  }
}

void adjust_cc_perm_preimage(uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* cluster_case_cts, uintptr_t* cluster_cc_perm_preimage) {
  uint32_t cluster_idx;
  uint32_t map_idx;
  uint32_t cluster_end;
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    map_idx = cluster_starts[cluster_idx];
    cluster_end = cluster_starts[cluster_idx + 1];
    if (cluster_case_cts[cluster_idx] * 2 < cluster_end - map_idx) {
      do {
        clear_bit_noct(cluster_cc_perm_preimage, cluster_map[map_idx] * 2);
      } while (++map_idx < cluster_end);
    } else {
      do {
	set_bit_noct(cluster_cc_perm_preimage, cluster_map[map_idx] * 2);
      } while (++map_idx < cluster_end);
    }
  }
}

int32_t cluster_include_and_reindex(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_include, uint32_t remove_size1, uintptr_t* pheno_c, uintptr_t indiv_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* new_cluster_ct_ptr, uint32_t** new_cluster_map_ptr, uint32_t** new_cluster_starts_ptr, uint32_t** cluster_case_cts_ptr, uintptr_t** cluster_cc_perm_preimage_ptr) {
  // If any individuals are excluded, this converts a unfiltered-index cluster
  // map to a collapsed-index cluster map (suitable for passing to the
  // per-cluster permutation function), allocating space for the latter.
  // Otherwise, this just sets the new cluster map to point to the old one.
  //
  // remove_size1 should be 0 or 1.  If it's 1, all size-1 clusters are
  // removed.
  //
  // If pheno_c is set, this also allocates and populates an array of
  // per-cluster case counts, and a permutation preimage.
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t old_assigned_ct = cluster_starts[cluster_ct];
  uint32_t new_cluster_ct = 0;
  uint32_t indiv_uidx = 0;
  uint32_t case_ct = 0;
  uint32_t assigned_ct = 0;
  uintptr_t* cluster_cc_perm_preimage = NULL;
  uint32_t* new_cluster_map;
  uint32_t* new_cluster_starts;
  uint32_t* cluster_case_cts;
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
    if (wkspace_alloc_ul_checked(cluster_cc_perm_preimage_ptr, 2 * ((indiv_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t))) {
      goto cluster_include_and_reindex_ret_NOMEM;
    }
    cluster_cc_perm_preimage = *cluster_cc_perm_preimage_ptr;
    collapse_copy_bitarr_to_vec_incl(unfiltered_indiv_ct, pheno_c, indiv_include, indiv_ct, cluster_cc_perm_preimage);
  }
  if ((indiv_ct == unfiltered_indiv_ct) && ((!remove_size1) || no_size1(cluster_ct, cluster_starts))) {
    *new_cluster_map_ptr = cluster_map;
    *new_cluster_ct_ptr = cluster_ct;
    *new_cluster_starts_ptr = cluster_starts;
    if (pheno_c) {
      if (wkspace_alloc_ui_checked(cluster_case_cts_ptr, cluster_ct * sizeof(int32_t))) {
	goto cluster_include_and_reindex_ret_NOMEM;
      }
      populate_cluster_case_cts(pheno_c, cluster_ct, cluster_map, cluster_starts, *cluster_case_cts_ptr);
      adjust_cc_perm_preimage(cluster_ct, cluster_map, cluster_starts, *cluster_case_cts_ptr, cluster_cc_perm_preimage);
    }
    return 0;
  }
  // scan to determine new memory allocation sizes
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    cluster_assigned_ct = 0;
    cluster_end = cluster_starts[cluster_idx + 1];
    for (map_idx = cluster_starts[cluster_idx]; map_idx < cluster_end; map_idx++) {
      cluster_assigned_ct += is_set(indiv_include, cluster_map[map_idx]);
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
  shrink_map = (assigned_ct < old_assigned_ct)? 1 : 0;
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
  }
  if (wkspace_alloc_ui_checked(&uidx_to_idx, unfiltered_indiv_ct * sizeof(int32_t))) {
    goto cluster_include_and_reindex_ret_NOMEM;
  }
  cluster_case_cts = *cluster_case_cts_ptr;
  fill_uidx_to_idx_incl(indiv_include, indiv_ct, uidx_to_idx);
  *new_cluster_ct_ptr = new_cluster_ct;
  cluster_read_idx = 1;
  map_idx = 0;
  cluster_idx = 0;
  cluster_end = 0;
  last_case_ct_incr = 0;
  for (map_read_idx = 0; map_read_idx < old_assigned_ct; map_read_idx++) {
    indiv_uidx = cluster_map[map_read_idx];
    if (!is_set(indiv_include, indiv_uidx)) {
      continue;
    }
    if (map_read_idx >= cluster_end) {
      if (cluster_idx) {
        if ((!remove_size1) || (map_idx - new_cluster_starts[cluster_idx] > 1)) {
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
      last_case_ct_incr = is_set(pheno_c, indiv_uidx);
    }
    new_cluster_map[map_idx++] = uidx_to_idx[indiv_uidx];
  }
  if (new_cluster_ct && ((!remove_size1) || (map_idx - new_cluster_starts[cluster_idx] > 1))) {
    if (pheno_c) {
      cluster_case_cts[new_cluster_ct - 1] = case_ct + last_case_ct_incr;
    }
    if (shrink_map) {
      new_cluster_starts[new_cluster_ct] = map_idx;
    }
  }
  if (pheno_c && new_cluster_ct) {
    adjust_cc_perm_preimage(new_cluster_ct, new_cluster_map, new_cluster_starts, cluster_case_cts, cluster_cc_perm_preimage);
  }
  wkspace_reset((unsigned char*)uidx_to_idx);
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
    tot_quotients[cluster_idx] = 4294967296LLU / (cluster_starts[cluster_idx + 1] - cluster_starts[cluster_idx]);
    magic_num(tot_quotients[cluster_idx], &(totq_magics[cluster_idx]), &(totq_preshifts[cluster_idx]), &(totq_postshifts[cluster_idx]), &(totq_incrs[cluster_idx]));
  }
  return 0;
}

int32_t read_genome() {
  int32_t retval = 0;
  return retval;
}
