#include "wdist_cluster.h"

#ifdef __APPLE__

#include <Accelerate/Accelerate.h>

#else // not __APPLE__

#ifndef NOLAPACK

#ifdef __cplusplus
extern "C" {
#endif

  typedef double __CLPK_doublereal;

#ifdef _WIN32

#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_STRUCTURE
#include "lapack/lapacke/include/lapacke.h"

  typedef int32_t __CLPK_integer;
  // int dsyevr_(); needed?

#else // not _WIN32

#include <cblas.h>
#ifdef __LP64__
  typedef int32_t __CLPK_integer;
#else
  typedef long int __CLPK_integer;
#endif
  // int dsyevr_(); needed?

#endif

#ifdef __cplusplus
}
#endif

#endif // NOLAPACK

#endif // __APPLE__

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
  cluster_ptr->max_ctrls = 0xffffffffU;
  cluster_ptr->min_ct = 1;
  cluster_ptr->mds_dim_ct = 0;
  cluster_ptr->min_ibm = 0.0;
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
  retval = sort_item_ids_noalloc(sorted_ids, id_map, unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, 0, 0, strcmp_deref);
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

void fill_unfiltered_indiv_to_cluster(uintptr_t unfiltered_indiv_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* indiv_to_cluster) {
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

int32_t fill_indiv_to_cluster(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* indiv_to_cluster, uint32_t* late_clidx_to_indiv_uidx) {
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t* cluster_map_pos = cluster_map;
  int32_t retval = 0;
  uint32_t* uidx_to_idx;
  uint32_t* cluster_end_ptr;
  uint32_t cluster_idx;
  uint32_t indiv_uidx;
  uint32_t indiv_idx;
  if (wkspace_alloc_ui_checked(&uidx_to_idx, unfiltered_indiv_ct * sizeof(int32_t))) {
    goto fill_indiv_to_cluster_ret_NOMEM;
  }
  fill_uidx_to_idx(indiv_exclude, indiv_ct, uidx_to_idx);
  fill_uint_one(indiv_to_cluster, indiv_ct);
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    cluster_end_ptr = &(cluster_map[cluster_starts[cluster_idx + 1]]);
    do {
      indiv_to_cluster[uidx_to_idx[*cluster_map_pos]] = cluster_idx;
    } while (++cluster_map_pos < cluster_end_ptr);
  }
  if (cluster_starts[cluster_ct] < indiv_ct) {
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
      if (indiv_to_cluster[indiv_idx] == 0xffffffffU) {
	late_clidx_to_indiv_uidx[cluster_idx - cluster_ct] = indiv_uidx;
	indiv_to_cluster[indiv_idx] = cluster_idx++;
      }
      indiv_uidx++;
    }
  }
  while (0) {
  fill_indiv_to_cluster_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
  wkspace_reset(wkspace_mark);
  return retval;
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
  fill_unfiltered_indiv_to_cluster(unfiltered_indiv_ct, cluster_ct, cluster_map, cluster_starts, indiv_to_cluster);
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

int32_t read_genome(char* read_genome_fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t* cluster_merge_prevented, double* cluster_sorted_ibs, uint32_t neighbor_n2, double* neighbor_quantiles, uint32_t* neighbor_qindices, uint32_t* ppc_fail_counts, double min_ppc, uint32_t is_max_dist, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* indiv_to_cluster) {
  unsigned char* wkspace_mark = wkspace_base;
  gzFile gz_infile = NULL;
  uint32_t neighbor_load_quantiles = neighbor_quantiles && cluster_sorted_ibs;
  uint32_t ppc_warning = cluster_merge_prevented? 0 : 1;
  uintptr_t loaded_entry_ct = 0;
  uint32_t ppc_fail = 0;
  char* idbuf = &(tbuf[MAXLINELEN]);
  char* sorted_ids;
  uint32_t* id_map;
  char* bufptr;
  char* fam_id;
  char* indiv_id;
  uint32_t indiv_idx1;
  uint32_t indiv_idx2;
  double cur_ibs;
  double cur_ppc;
  uintptr_t tcoord;
  uint32_t uii;
  int32_t ii;
  int32_t retval;
  retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_indiv_ct, indiv_exclude, unfiltered_indiv_ct - indiv_ct, person_ids, max_person_id_len, 0, 1, strcmp_deref);
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
      goto read_genome_ret_INVALID_FORMAT_3;
    }
    bufptr = skip_initial_spaces(tbuf);
  } while (is_eoln_kns(*bufptr));
  // a little bit of input validation
  if (memcmp(bufptr, "FID1", 4)) {
    logprint("Error: Invalid --read-genome input file.\n");
    goto read_genome_ret_INVALID_FORMAT;
  }
  while (gzgets(gz_infile, tbuf, MAXLINELEN)) {
    if (!tbuf[MAXLINELEN - 1]) {
      goto read_genome_ret_INVALID_FORMAT_3;
    }
    fam_id = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*fam_id)) {
      continue;
    }
    indiv_id = next_item(fam_id);
    bufptr = next_item_mult(indiv_id, 2);
    if (no_more_items(bufptr)) {
      goto read_genome_ret_INVALID_FORMAT_4;
    }
    ii = bsearch_fam_indiv(idbuf, sorted_ids, max_person_id_len, indiv_ct, fam_id, indiv_id);
    if (ii == -1) {
      continue;
    }
    indiv_idx1 = id_map[(uint32_t)ii];
    fam_id = next_item(indiv_id);
    indiv_id = bufptr;
    ii = bsearch_fam_indiv(idbuf, sorted_ids, max_person_id_len, indiv_ct, fam_id, indiv_id);
    if (ii == -1) {
      continue;
    }
    indiv_idx2 = id_map[(uint32_t)ii];
    if (indiv_idx2 == indiv_idx1) {
      logprint("Error: FID1/IID1 matches FID2/IID2 in --read-genome input file line.\n");
      goto read_genome_ret_INVALID_FORMAT;
    }
    bufptr = next_item_mult(indiv_id, 8); // distance
    fam_id = next_item(bufptr); // repurposed to PPC test value
    if (no_more_items(fam_id)) {
      goto read_genome_ret_INVALID_FORMAT_4;
    }
    if (min_ppc != 0.0) {
      if (sscanf(fam_id, "%lg", &cur_ppc) != 1) {
	logprint("Error: Invalid PPC test value in --read-genome input file.\n");
	goto read_genome_ret_INVALID_FORMAT;
      }
      ppc_fail = (cur_ppc < min_ppc)? 1 : 0;
      if (ppc_fail && ppc_fail_counts) {
	ppc_fail_counts[indiv_idx1] += 1;
	ppc_fail_counts[indiv_idx2] += 1;
      }
    }
    if (sscanf(bufptr, "%lg", &cur_ibs) != 1) {
      logprint("Error: Invalid IBS value in --read-genome input file.\n");
      goto read_genome_ret_INVALID_FORMAT;
    }
    if (neighbor_load_quantiles) {
      update_neighbor(indiv_ct, neighbor_n2, indiv_idx1, indiv_idx2, cur_ibs, neighbor_quantiles, neighbor_qindices);
    }
    loaded_entry_ct++;
    if (cluster_ct) {
      indiv_idx1 = indiv_to_cluster[indiv_idx1];
      indiv_idx2 = indiv_to_cluster[indiv_idx2];
      if (indiv_idx1 == indiv_idx2) {
	if (ppc_fail && (!ppc_warning)) {
	  logprint("Warning: Initial cluster assignment violates PPC test constraint.\n");
	  ppc_warning = 1;
	}
	continue;
      }
    }
    if (indiv_idx2 < indiv_idx1) {
      uii = indiv_idx1;
      indiv_idx1 = indiv_idx2;
      indiv_idx2 = uii;
    }
    tcoord = tri_coord_no_diag(indiv_idx1, indiv_idx2);
    if (ppc_fail) {
      set_bit_ul(cluster_merge_prevented, tcoord);
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
  if (loaded_entry_ct != (indiv_ct * (indiv_ct - 1)) / 2) {
    sprintf(logbuf, "Error: %s does not include all individual pairs.\n", read_genome_fname);
    goto read_genome_ret_INVALID_FORMAT_2;
  }
  while (0) {
  read_genome_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  read_genome_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  read_genome_ret_INVALID_FORMAT_4:
    sprintf(logbuf, "Error: Fewer entries than expected in %s line.\n", read_genome_fname);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  read_genome_ret_INVALID_FORMAT_3:
    sprintf(logbuf, "Error: Pathologically long line in %s.\n", read_genome_fname);
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

uint32_t cluster_main(uintptr_t cluster_ct, uintptr_t* merge_prevented, uintptr_t list_size, uint32_t* sorted_ibs_indices, uint32_t* cluster_index, uint32_t* cur_cluster_sizes, uint32_t indiv_ct, uint32_t* cur_cluster_case_cts, uint32_t case_ct, uint32_t ctrl_ct, uint32_t* cur_cluster_remap, Cluster_info* cp, uintptr_t* ibs_ties, uint32_t* merge_sequence) {
  uint32_t is_old_tiebreaks = cp->modifier & CLUSTER_OLD_TIEBREAKS;
  uint32_t* list_end = &(sorted_ibs_indices[list_size]);
  uint32_t max_merge = cluster_ct - cp->min_ct;
  uint32_t max_size = cp->max_size;
  uint32_t max_cases = cp->max_cases;
  uint32_t max_ctrls = cp->max_ctrls;
  uint32_t size_restriction = (max_size < indiv_ct)? 1 : 0;
  uint32_t case_restriction = (max_cases < case_ct)? 1 : 0;
  uint32_t ctrl_restriction = (max_ctrls < ctrl_ct)? 1 : 0;
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
          tie_end = &(sorted_ibs_indices[1 + next_non_set_unsafe(ibs_ties, siptr - sorted_ibs_indices)]);
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
          tie_end = &(sorted_ibs_indices[1 + next_non_set_unsafe(ibs_ties, siptr - sorted_ibs_indices)]);
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
      }
    }
    tcoord1 = (clidx_large * (clidx_large - 1)) / 2;
    tcoord2 = (clidx_small * (clidx_small - 1)) / 2;
    umm = clidx_small << 16;
    if (!sccr) {
      for (uii = 0; uii < clidx_small; uii++) {
	if (cur_cluster_remap[uii] == uii) {
	  if (is_set(merge_prevented, tcoord1 + uii)) {
	    set_bit_noct(merge_prevented, tcoord2 + uii);
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
	if (cur_cluster_remap[uii] == uii) {
	  umm = tri_coord_no_diag_32(clidx_small, uii);
	  if (is_set(merge_prevented, tcoord1 + uii)) {
	    set_bit_noct(merge_prevented, umm);
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
	if (cur_cluster_remap[uii] == uii) {
	  umm = tri_coord_no_diag_32(clidx_small, uii);
	  if (is_set(merge_prevented, tri_coord_no_diag_32(clidx_large, uii))) {
	    set_bit_noct(merge_prevented, umm);
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
	if (cur_cluster_remap[uii] == uii) {
	  if (is_set(merge_prevented, tcoord1 + uii) || (size_restriction && (cur_size + cur_cluster_sizes[uii] > max_size)) || (case_restriction && (cur_cases + cur_cluster_case_cts[uii] > max_cases)) || (ctrl_restriction && (cur_ctrls + cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > max_ctrls))) {
	    set_bit_noct(merge_prevented, tcoord2 + uii);
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
	if (cur_cluster_remap[uii] == uii) {
	  umm = tri_coord_no_diag_32(clidx_small, uii);
	  if (is_set(merge_prevented, tcoord1 + uii) || (size_restriction && (cur_size + cur_cluster_sizes[uii] > max_size)) || (case_restriction && (cur_cases + cur_cluster_case_cts[uii] > max_cases)) || (ctrl_restriction && (cur_ctrls + cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > max_ctrls))) {
	    set_bit_noct(merge_prevented, umm);
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
      for (uii = clidx_large + 1; uii < cluster_ct; uii++) {
	if (cur_cluster_remap[uii] == uii) {
	  umm = tri_coord_no_diag_32(clidx_small, uii);
	  if (is_set(merge_prevented, tri_coord_no_diag_32(clidx_large, uii)) || (size_restriction && (cur_size + cur_cluster_sizes[uii] > max_size)) || (case_restriction && (cur_cases + cur_cluster_case_cts[uii] > max_cases)) || (ctrl_restriction && (cur_ctrls + cur_cluster_sizes[uii] - cur_cluster_case_cts[uii] > max_ctrls))) {
	    set_bit_noct(merge_prevented, umm);
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

uint32_t cluster_group_avg_main(uintptr_t cluster_ct, uintptr_t* merge_prevented, uintptr_t heap_size, double* sorted_ibs, uint32_t* sorted_ibs_indices, uint32_t* cur_cluster_sizes, uint32_t indiv_ct, uint32_t* cur_cluster_case_cts, uint32_t case_ct, uint32_t ctrl_ct, uint32_t* cur_cluster_remap, Cluster_info* cp, uint32_t* merge_sequence) {
  uint32_t merge_ct = 0;
  return merge_ct;
}

void write_cluster1(FILE* outfile, uint32_t clidx, char* person_ids, uintptr_t max_person_id_len, uintptr_t* pheno_c, uint32_t* indiv_idx_to_uidx, uint32_t* merge_sequence, uint32_t merge_ct) {
  // Manually manage recursion, to minimize crash risk when there are e.g. 500k
  // clusters.  Fortunately, no recursion stack is even needed.
  char* sptr;
  char* sptr2;
  uint32_t msidx;
 write_cluster1_recurse:
  putc(' ', outfile);
  sptr = &(person_ids[indiv_idx_to_uidx[clidx] * max_person_id_len]);
  sptr2 = (char*)memchr(sptr, '\t', max_person_id_len);
  fwrite(sptr, 1, sptr2 - sptr, outfile);
  putc('_', outfile);
  fputs(&(sptr2[1]), outfile);
  if (pheno_c) {
    if (is_set(pheno_c, indiv_idx_to_uidx[clidx])) {
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

void write_cluster1_oitc(FILE* outfile, uint32_t clidx, char* person_ids, uintptr_t max_person_id_len, uintptr_t* pheno_c, uint32_t* orig_cluster_map, uint32_t* orig_cluster_starts, uint32_t* late_clidx_to_indiv_uidx, uint32_t orig_within_ct, uint32_t* cluster_remap, uint32_t* merge_sequence, uint32_t merge_ct) {
  char* sptr;
  char* sptr2;
  uint32_t uii;
  uint32_t ujj;
 write_cluster1_oitc_recurse:
  putc(' ', outfile);
  if (clidx >= orig_within_ct) {
    sptr = &(person_ids[late_clidx_to_indiv_uidx[clidx - orig_within_ct] * max_person_id_len]);
    sptr2 = (char*)memchr(sptr, '\t', max_person_id_len);
    fwrite(sptr, 1, sptr2 - sptr, outfile);
    putc('_', outfile);
    fputs(&(sptr2[1]), outfile);
    if (pheno_c) {
      if (is_set(pheno_c, late_clidx_to_indiv_uidx[clidx - orig_within_ct])) {
	fputs("(2)", outfile);
      } else {
	fputs("(1)", outfile);
      }
    }
  } else {
    ujj = orig_cluster_starts[clidx + 1];
    for (uii = orig_cluster_starts[clidx]; uii < ujj; uii++) {
      sptr = &(person_ids[orig_cluster_map[uii] * max_person_id_len]);
      sptr2 = (char*)memchr(sptr, '\t', max_person_id_len);
      fwrite(sptr, 1, sptr2 - sptr, outfile);
      putc('_', outfile);
      fputs(&(sptr2[1]), outfile);
      if (pheno_c) {
	if (is_set(pheno_c, orig_cluster_map[uii])) {
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

int32_t write_cluster_solution(char* outname, char* outname_end, uint32_t* orig_indiv_to_cluster, uintptr_t indiv_ct, uint32_t* orig_cluster_map, uint32_t* orig_cluster_starts, uint32_t* late_clidx_to_indiv_uidx, uint32_t orig_within_ct, uint32_t orig_cluster_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t* pheno_c, uint32_t* indiv_idx_to_uidx, Cluster_info* cp, uint32_t* cluster_remap, uint32_t* clidx_table_space, uint32_t merge_ct, uint32_t* merge_sequence) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uint32_t only2 = cp->modifier & CLUSTER_ONLY2;
  uint32_t report_pheno = (cp->modifier & CLUSTER_CC) || (cp->max_ctrls != 0xffffffffU);
  int32_t retval = 0;
  char wbuf[16];
  uint32_t* clidx_remap = &(clidx_table_space[(((uintptr_t)orig_cluster_ct) * (orig_cluster_ct - 1)) >> 1]);
  char* sptr;
  char* sptr2;
  char* wptr;
  uint32_t* cur_remap;
  uint32_t indiv_idx;
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
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    if (orig_indiv_to_cluster) {
      clidx = cluster_remap[orig_indiv_to_cluster[indiv_idx]];
    } else {
      clidx = cluster_remap[indiv_idx];
    }
    sptr = &(person_ids[indiv_idx_to_uidx[indiv_idx] * max_person_id_len]);
    sptr2 = (char*)memchr(sptr, '\t', max_person_id_len);
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
        if (!orig_indiv_to_cluster) {
	  write_cluster1(outfile, clidx, person_ids, max_person_id_len, report_pheno? pheno_c : NULL, indiv_idx_to_uidx, merge_sequence, merge_ct);
	} else {
	  write_cluster1_oitc(outfile, clidx, person_ids, max_person_id_len, report_pheno? pheno_c : NULL, orig_cluster_map, orig_cluster_starts, late_clidx_to_indiv_uidx, orig_within_ct, cluster_remap, merge_sequence, merge_ct);
	}
	if (putc('\n', outfile) == EOF) {
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
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      sptr = &(person_ids[indiv_idx_to_uidx[indiv_idx] * max_person_id_len]);
      sptr2 = (char*)memchr(sptr, '\t', max_person_id_len);
      if (fwrite_checked(sptr, sptr2 - sptr, outfile)) {
	goto write_cluster_solution_ret_WRITE_FAIL;
      }
      putc(' ', outfile);
      fputs(&(sptr2[1]), outfile);
      putc('\t', outfile);
      if (orig_indiv_to_cluster) {
        clidx = orig_indiv_to_cluster[indiv_idx];
      } else {
	clidx = indiv_idx;
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
      if (putc('\n', outfile) == EOF) {
	goto write_cluster_solution_ret_WRITE_FAIL;
      }
    }
    putc('\n', outfile);
    if (fclose_null(&outfile)) {
      goto write_cluster_solution_ret_WRITE_FAIL;
    }
    *outname_end = '\0';
    putchar('\r');
    sprintf(logbuf, "Cluster solution written to %s.cluster{1,2,3%s}.\n", outname, (cp->modifier & CLUSTER_MISSING)? ".missing" : "");
  } else {
    *outname_end = '\0';
    putchar('\r');
    sprintf(logbuf, "Cluster solution written to %s.cluster2.\n", outname);
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
int32_t mds_plot(char* outname, char* outname_end, uintptr_t* indiv_exclude, uintptr_t indiv_ct, uint32_t* indiv_idx_to_uidx, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, uint32_t cur_cluster_ct, uint32_t merge_ct, uint32_t* orig_indiv_to_cluster, uint32_t* cur_cluster_remap, uint32_t dim_ct, uint32_t is_mds_cluster, double* dists) {
  FILE* outfile = NULL;
  uintptr_t final_cluster_ct = cur_cluster_ct - merge_ct;
  double grand_mean = 0.0;
  double nz = 0.0;
  double zz = -1.0;
  uintptr_t ulii = 0;
  __CLPK_integer info = 0;
  __CLPK_integer lwork = -1;
  __CLPK_integer liwork = -1;
  int32_t retval = 0;
  double* main_matrix;
  double* column_means;
  uint32_t* final_cluster_remap;
  uint32_t* final_cluster_sizes;
  double* dptr;
  double* dptr2;
  char* wptr;
  char* wptr2;
  uintptr_t indiv_idx;
  uintptr_t clidx1;
  uintptr_t clidx2;
  uint32_t dim_idx;
  uint32_t uii;
  uint32_t ujj;
  double dxx;
  double dyy;
  __CLPK_integer mdim;
  __CLPK_integer i1;
  __CLPK_integer i2;
  __CLPK_integer out_m;
  __CLPK_integer ldz;
  double* work;
  double* out_w;
  double* out_z;
  double optim_lwork;
  __CLPK_integer optim_liwork;
  __CLPK_integer* iwork;
  __CLPK_integer* isuppz;
  double* sqrt_eigvals;
  final_cluster_remap = (uint32_t*)malloc(cur_cluster_ct * sizeof(int32_t));
  if (!final_cluster_remap) {
    goto mds_plot_ret_NOMEM;
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
    wkspace_reset((unsigned char*)dists);
    if (wkspace_alloc_d_checked(&main_matrix, indiv_ct * indiv_ct * sizeof(double))) {
      goto mds_plot_ret_NOMEM;
    }
    // expand triangular diagonal-free matrix to square matrix
    ulii = ((indiv_ct - 1) * (indiv_ct - 2)) >> 1;
    for (indiv_idx = indiv_ct - 1; indiv_idx;) {
      memcpy(&(main_matrix[indiv_idx * indiv_ct]), &(main_matrix[ulii]), indiv_idx * sizeof(double));
      ulii -= (--indiv_idx);
    }
    ulii = indiv_ct + 1;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      main_matrix[indiv_idx * ulii] = 0.0;
    }
    ulii = indiv_ct;
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
  // no need to fill upper right

  // see eigen_lapack() in PLINK lapackf.cpp (though we use dsyevr_ instead of
  // dsyevx_)
  mdim = ulii;
  i2 = mdim;
  if (dim_ct > ulii) {
    dim_ct = ulii;
  }
  i1 = i2 + 1 - dim_ct;
  if (wkspace_alloc_d_checked(&out_w, dim_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&out_z, dim_ct * ulii * sizeof(double))) {
    goto mds_plot_ret_NOMEM;
  }
  isuppz = (__CLPK_integer*)wkspace_alloc(2 * dim_ct * sizeof(__CLPK_integer));
  if (!isuppz) {
    goto mds_plot_ret_NOMEM;
  }
  fill_double_zero(out_w, dim_ct);
  fill_double_zero(out_z, dim_ct * ulii);
  fill_int_zero((int32_t*)isuppz, 2 * dim_ct * (sizeof(__CLPK_integer) / sizeof(int32_t)));
  ldz = mdim;

  dsyevr_((char*)"V", (char*)"I", (char*)"U", &mdim, main_matrix, &mdim, &nz, &nz, &i1, &i2, &zz, &out_m, out_w, out_z, &ldz, isuppz, &optim_lwork, &lwork, &optim_liwork, &liwork, &info);
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
  dsyevr_((char*)"V", (char*)"I", (char*)"U", &mdim, main_matrix, &mdim, &nz, &nz, &i1, &i2, &zz, &out_m, out_w, out_z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

  // * out_w[0..(dim_ct-1)] contains eigenvalues
  // * out_z[(ii*ulii)..(ii*ulii + ulii - 1)] is eigenvector corresponding to
  //   out_w[ii]
  wkspace_reset((unsigned char*)work);
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

  memcpy(outname_end, ".mds", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto mds_plot_ret_OPEN_FAIL;
  }
  sprintf(tbuf, "%%%us %%%us    SOL ", plink_maxfid, plink_maxiid);
  if (fprintf(outfile, tbuf, "FID", "IID") < 0) {
    goto mds_plot_ret_WRITE_FAIL;
  }
  tbuf[22] = ' ';
  for (dim_idx = 0; dim_idx < dim_ct; dim_idx++) {
    wptr = uint32_write(tbuf, dim_idx + 1);
    uii = wptr - tbuf;
    wptr2 = memseta(&(tbuf[10]), 32, 11 - uii);
    *wptr2++ = 'C';
    memcpy(wptr2, tbuf, uii);
    fwrite(&(tbuf[10]), 1, 13, outfile);
  }
  if (putc('\n', outfile) == EOF) {
    goto mds_plot_ret_WRITE_FAIL;
  }
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    wptr2 = &(person_ids[indiv_idx_to_uidx[indiv_idx] * max_person_id_len]);
    uii = strlen_se(wptr2);
    wptr = fw_strcpyn(plink_maxfid, uii, wptr2, tbuf);
    *wptr++ = ' ';
    wptr = fw_strcpy(plink_maxiid, &(wptr2[uii + 1]), wptr);
    *wptr++ = ' ';
    if (orig_indiv_to_cluster) {
      uii = orig_indiv_to_cluster[indiv_idx];
    } else {
      uii = indiv_idx;
    }
    uii = final_cluster_remap[uii];
    wptr = uint32_writew6x(wptr, uii, ' ');
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto mds_plot_ret_WRITE_FAIL;
    }
    if (!is_mds_cluster) {
      dptr = &(main_matrix[(indiv_idx + 1) * dim_ct]);
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
    if (putc('\n', outfile) == EOF) {
      goto mds_plot_ret_WRITE_FAIL;
    }
  }
  sprintf(logbuf, "MDS solution written to %s.\n", outname);
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
  wkspace_reset((unsigned char*)dists);
  return retval;
}
#endif
