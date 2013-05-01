#include "wdist_common.h"

int32_t cnv_subset_load(char* subset_fname, char** subset_list_ptr, uintptr_t* subset_ct_ptr, uintptr_t* max_subset_name_len_ptr, uintptr_t* topsize_ptr) {
  FILE* subset_file = NULL;
  uintptr_t subset_ct = 0;
  uintptr_t max_subset_name_len = 1;
  int32_t retval = 0;
  char* bufptr;
  char* slptr;
  uint64_t alloc_req;
  uintptr_t ulii;
  // Two-pass load:
  // 1. Determine subset_ct and max_subset_name_len.
  // 2. Allocate temporary block at end of workspace, and load the strings.
  if (fopen_checked(&subset_file, subset_fname, "r")) {
    goto cnv_subset_load_ret_OPEN_FAIL;
  }
  while (fgets(tbuf, MAXLINELEN, subset_file)) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in --cnv-subset file.\n");
      goto cnv_subset_load_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (!is_eoln_kns(*bufptr)) {
      ulii = strlen_se(bufptr);
      if (ulii >= max_subset_name_len) {
        max_subset_name_len = ulii + 1;
      }
      subset_ct++;
    }
  }
  if (!feof(subset_file)) {
    goto cnv_subset_load_ret_READ_FAIL;
  } else if (!subset_ct) {
    logprint("Error: Empty --cnv-subset file.\n");
    goto cnv_subset_load_ret_INVALID_FORMAT;
  }
  alloc_req = subset_ct;
  alloc_req *= max_subset_name_len;
  alloc_req = (alloc_req + (CACHELINE - 1)) & (~((CACHELINE - 1) * ONELU));
  if (alloc_req >= wkspace_left) {
    goto cnv_subset_load_ret_NOMEM;
  }
  *topsize_ptr = alloc_req;
  *subset_list_ptr = (char*)(&(wkspace_base[wkspace_left - (*topsize_ptr)]));
  wkspace_left -= (*topsize_ptr);
  rewind(subset_file);
  slptr = *subset_list_ptr;
  while (fgets(tbuf, MAXLINELEN, subset_file)) {
    bufptr = skip_initial_spaces(tbuf);
    if (!is_eoln_kns(*bufptr)) {
      ulii = strlen_se(bufptr);
      memcpy(slptr, bufptr, ulii);
      slptr[ulii] = '\0';
      slptr = &(slptr[max_subset_name_len]);
    }
  }
  qsort(*subset_list_ptr, subset_ct, max_subset_name_len, strcmp_casted);
  *subset_ct_ptr = subset_ct;
  *max_subset_name_len_ptr = max_subset_name_len;
  while (0) {
  cnv_subset_load_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  cnv_subset_load_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  cnv_subset_load_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  cnv_subset_load_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(subset_file);
  return retval;
}

char* cnv_intersect_filter_type_to_str(uint32_t intersect_filter_type) {
  if (intersect_filter_type == CNV_INTERSECT) {
    return (char*)"intersect";
  } else if (intersect_filter_type == CNV_EXCLUDE) {
    return (char*)"exclude";
  } else {
    return (char*)"count";
  }
}

// maximum 25
#define SMALL_INTERVAL_BITS 18
#define SMALL_INTERVAL_MAX_SIZE ((1 << SMALL_INTERVAL_BITS) - 1)

int32_t cnv_intersect_load(uint32_t intersect_filter_type, char* intersect_filter_fname, char* subset_list, uintptr_t subset_ct, uintptr_t max_subset_name_len, uintptr_t topsize, uintptr_t* il_chrom_start_small, uintptr_t* il_chrom_start_large, uint32_t* il_chrom_max_width_small, uint32_t* il_chrom_max_width_large, uint64_t** il_small_ptr, uint64_t** il_large_ptr, uint32_t species, uint64_t chrom_mask) {
  // We store intervals in sorted order, with the center of each interval in
  // the high-order bits, and the size in the low-order bits.  (Chromosome
  // beginnings and endings are stored in small external arrays.)  We also
  // track each chromosome's largest interval size.
  // If all intervals are small, this enables very efficient intersection
  // checks.  However, if there is a single very large interval, the standard
  // intersection checking algorithm might have to assume the worst for all the
  // smaller intervals as well.
  //
  // To guard against that, we split intervals into two bins.
  // Small tier: All intervals where endpos - startpos < 2^SMALL_INTERVAL_BITS.
  //             This constant should be chosen so that if there is a large
  //             number of intervals, chances are that the vast majority
  //             satisfies this condition.
  // Large tier: All intervals that aren't in the small tier.
  // This way, whenever we check for an intersection, we can usually skip
  // almost all the "small tier" intervals regardless of the largest interval
  // size.
  FILE* intersect_file = NULL;
  uintptr_t max_interval_ct = wkspace_left / 9;
  uintptr_t small_interval_ct = 0;
  uintptr_t large_interval_ct = 0;
  // Storage format before sorting:
  // Small: high 16 bits = chrom, middle 32 bits = center pos * 2,
  //        bottom 16 bits = size
  // Large: top bit = zero, next 32 bits = center pos * 2, bottom 31 = size
  //        [chrom information stored separately, initially in reverse order]
  uint64_t* il_small = (uint64_t*)wkspace_base;
  uint64_t* tmp_il_large = &(il_small[2 * max_interval_ct]); // grows down
  char* tmp_il_large_chroms = (char*)tmp_il_large; // grows up
  uint32_t reverse_warning_given = 0;
  int32_t retval = 0;
  uint64_t* il_large;
  char* bufptr;
  char* bufptr2;
  uint64_t ullii;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t cur_width;
  uintptr_t max_width;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  char cur_chrom;
  char cc;
  *il_small_ptr = il_small;
  if (fopen_checked(&intersect_file, intersect_filter_fname, "r")) {
    goto cnv_intersect_load_ret_OPEN_FAIL;
  }
  while (fgets(tbuf, MAXLINELEN, intersect_file)) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Pathologically long line in --cnv-%s file.\n", cnv_intersect_filter_type_to_str(intersect_filter_type));
      goto cnv_intersect_load_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (!is_eoln_kns(*bufptr)) {
      // CHR, BP1, BP2, subset name
      bufptr2 = next_item_mult(bufptr, 2);
      if (no_more_items_kns(bufptr2)) {
	sprintf(logbuf, "Error: Fewer items than expected in --cnv-%s file.\n", cnv_intersect_filter_type_to_str(intersect_filter_type));
	goto cnv_intersect_load_ret_INVALID_FORMAT;
      }
      ii = marker_code(species, bufptr);
      if (ii == -1) {
	sprintf(logbuf, "Error: Invalid chromosome code in --cnv-%s file.\n", cnv_intersect_filter_type_to_str(intersect_filter_type));
	goto cnv_intersect_load_ret_INVALID_FORMAT;
      }
      if (!((chrom_mask >> ((uint32_t)ii)) & 1)) {
	continue;
      }
      bufptr = next_item(bufptr);
      ii = atoi(bufptr);
      if (ii < 0) {
	sprintf(logbuf, "Error: Negative position in --cnv-%s file.\n", cnv_intersect_filter_type_to_str(intersect_filter_type));
	goto cnv_intersect_load_ret_INVALID_FORMAT;
      }
      jj = ii;
      ii = atoi(bufptr2);
      if (ii < 0) {
	sprintf(logbuf, "Error: Negative position in --cnv-%s file.\n", cnv_intersect_filter_type_to_str(intersect_filter_type));
	goto cnv_intersect_load_ret_INVALID_FORMAT;
      }
      if (ii < jj) {
	if (!reverse_warning_given) {
	  sprintf(logbuf, "Warning: End of range before start of range in --cnv-%s file.\n", cnv_intersect_filter_type_to_str(intersect_filter_type));
	  reverse_warning_given = 1;
	}
	kk = ii;
	ii = jj;
	jj = kk;
      }
      if (subset_ct) {
	bufptr = next_item(bufptr2);
	if (!no_more_items_kns(bufptr)) {
	  ulii = strlen_se(bufptr);
	  if (ulii < max_subset_name_len) {
	    bufptr[ulii] = '\0';
	    if (bsearch_str(bufptr, subset_list, max_subset_name_len, 0, subset_ct - 1) != -1) {
	      continue;
	    }
	  }
	}
      }
      if (small_interval_ct + large_interval_ct == max_interval_ct) {
        goto cnv_intersect_load_NOMEM;
      }
      kk = jj - ii;
      if (kk > SMALL_INTERVAL_MAX_SIZE) {
	tmp_il_large--;
	*tmp_il_large = (((uint64_t)(((uint32_t)ii) + ((uint32_t)jj))) << 31) | ((uint64_t)((uint32_t)kk));
	large_interval_ct++;
      } else {
	il_small[small_interval_ct++] = (((uint64_t)(((uint32_t)ii) + ((uint32_t)jj))) << SMALL_INTERVAL_BITS) | ((uint64_t)((uint32_t)kk));
      }
    }
  }
  if (!feof(intersect_file)) {
    goto cnv_intersect_load_ret_READ_FAIL;
  }
  *il_large_ptr = &(il_small[small_interval_ct]);
  if ((!small_interval_ct) && (!large_interval_ct)) {
    if (subset_ct) {
      fill_ulong_zero(il_chrom_start_small, MAX_POSSIBLE_CHROM + 1);
      fill_ulong_zero(il_chrom_start_large, MAX_POSSIBLE_CHROM + 1);
      logprint("Warning: All intervals filtered out by --cnv-subset.\n");
      goto cnv_intersect_load_ret_1;
    }
    sprintf(logbuf, "Error: Empty --cnv-%s file.\n", cnv_intersect_filter_type_to_str(intersect_filter_type));
    goto cnv_intersect_load_ret_INVALID_FORMAT;
  }
  if (small_interval_ct) {
#ifdef __cplusplus
    std::sort((int64_t*)il_small, (int64_t*)(&(il_small[small_interval_ct])));
#else
    qsort((int64_t*)il_small, small_interval_ct, sizeof(int64_t), llcmp);
#endif
    il_chrom_start_small[0] = 0;
    cur_chrom = 0;
    max_width = 0;
    for (ulii = 0; ulii < small_interval_ct; ulii++) {
      ullii = il_small[ulii];
      cc = (char)(ullii >> (SMALL_INTERVAL_BITS + 32));
      if (cc > cur_chrom) {
	il_chrom_max_width_small[cur_chrom] = max_width;
	do {
	  il_chrom_start_small[++cur_chrom] = ulii;
	} while (cur_chrom < cc);
	max_width = 0;
      }
      cur_width = ullii & (SMALL_INTERVAL_MAX_SIZE * 1LLU);
      if (cur_width > max_width) {
	max_width = cur_width;
      }
      il_small[ulii] = ((ullii >> SMALL_INTERVAL_BITS) << 32) | ((uint64_t)cur_width);
    }
    il_chrom_max_width_small[cur_chrom] = max_width;
    do {
      il_chrom_start_small[++cur_chrom] = small_interval_ct;
    } while (cur_chrom < MAX_POSSIBLE_CHROM);
  } else {
    fill_ulong_zero(il_chrom_start_small, MAX_POSSIBLE_CHROM + 1);
  }
  il_large = &(il_small[small_interval_ct]);
  if (large_interval_ct) {
    if (large_interval_ct > 1) {
      uljj = large_interval_ct / 2;
      ulkk = large_interval_ct - 1;
      for (ulii = 0; ulii < uljj; ulii++) {
        cc = tmp_il_large_chroms[ulii];
	tmp_il_large_chroms[ulii] = tmp_il_large_chroms[ulkk];
	tmp_il_large_chroms[ulkk--] = cc;
      }
    }
    qsort_ext(tmp_il_large_chroms, large_interval_ct, sizeof(char), char_cmp_deref, (char*)tmp_il_large, sizeof(int64_t));
    il_chrom_start_large[0] = 0;
    cur_chrom = 0;
    for (ulii = 0; ulii < large_interval_ct; ulii++) {
      cc = tmp_il_large_chroms[ulii];
      if (cc > cur_chrom) {
	do {
	  il_chrom_start_large[++cur_chrom] = ulii;
	} while (cur_chrom < cc);
      }
    }
    do {
      il_chrom_start_large[++cur_chrom] = large_interval_ct;
    } while (cur_chrom < MAX_POSSIBLE_CHROM);
    ulii = il_chrom_start_large[0];
    for (cur_chrom = 0; cur_chrom < MAX_POSSIBLE_CHROM; cur_chrom++) {
      uljj = il_chrom_start_large[cur_chrom + 1];
      if (uljj > ulii) {
#ifdef __cplusplus
	std::sort((int64_t*)(&(tmp_il_large[ulii])), (int64_t*)(&(tmp_il_large[uljj])));
#else
        qsort((int64_t*)(&(tmp_il_large[ulii])), uljj - ulii, sizeof(int64_t), llcmp);
#endif
	ulii = uljj;
      }
    }
    ulii = 0;
    for (cur_chrom = 0; cur_chrom < MAX_POSSIBLE_CHROM; cur_chrom++) {
      uljj = il_chrom_start_large[cur_chrom + 1];
      max_width = 0;
      while (ulii < uljj) {
	ullii = tmp_il_large[ulii];
	cur_width = ullii & 0x7fffffffLLU;
	if (cur_width > max_width) {
	  max_width = cur_width;
	}
	il_large[ulii++] = ((ullii >> 31) << 32) | ((uint64_t)cur_width);
      }
      il_chrom_max_width_large[cur_chrom] = max_width;
    }
  } else {
    fill_ulong_zero(il_chrom_start_large, MAX_POSSIBLE_CHROM + 1);
  }
  il_small = (uint64_t*)wkspace_alloc((small_interval_ct + large_interval_ct) * sizeof(int64_t));
  while (0) {
  cnv_intersect_load_NOMEM:
    retval = RET_NOMEM;
    break;
  cnv_intersect_load_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  cnv_intersect_load_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  cnv_intersect_load_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    logprintb();
    break;
  }
 cnv_intersect_load_ret_1:
  fclose_cond(intersect_file);
  return retval;
}

int32_t cnv_first_nonheader_line(FILE* cnvfile) {
  int32_t retval = 0;
  char* bufptr;
  rewind(cnvfile);
  // assumes tbuf[MAXLINELEN - 1] is initialized to space
  do {
    if (!fgets(tbuf, MAXLINELEN, cnvfile)) {
      goto cnv_first_nonheader_line_fgets_fail;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("\nError: Pathologically long line in .cnv file.\n");
      goto cnv_first_nonheader_line_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
  } while (is_eoln_kns(*bufptr));
  if ((strlen_se(bufptr) == 3) && (!memcmp(bufptr, "FID", 3))) {
    if (!fgets(tbuf, MAXLINELEN, cnvfile)) {
      goto cnv_first_nonheader_line_fgets_fail;
    }
  }
  while (0) {
  cnv_first_nonheader_line_fgets_fail:
    if (!feof(cnvfile)) {
      retval = RET_READ_FAIL;
      break;
    }
    logprint("\nError: Empty .cnv file.\n");
    // fall through
  cnv_first_nonheader_line_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  return retval;
}

int32_t cnv_make_map_write(FILE* new_mapfile, uint64_t ullii) {
  uint32_t chrom_idx = ullii >> 32;
  uint32_t bp_pos = (uint32_t)ullii;
  char* wptr = uint32_write(tbuf, chrom_idx);
  wptr = memcpya(wptr, "\tp", 2);
  wptr = uint32_writex(wptr, chrom_idx, '-');
  wptr = uint32_write(wptr, bp_pos);
  wptr = memcpyl3a(wptr, "\t0\t");
  wptr = uint32_writex(wptr, bp_pos, '\n');
  return fwrite_checked(tbuf, wptr - tbuf, new_mapfile);
}

int32_t cnv_make_map(FILE* cnvfile, char* new_mapname, uint32_t cnv_calc_type, uint32_t min_seglen, uint32_t max_seglen, double min_score, double max_score, uint32_t min_sites, uint32_t max_sites, uintptr_t* il_chrom_start_small, uintptr_t* il_chrom_start_large, uint32_t* il_chrom_max_width_small, uint32_t* il_chrom_max_width_large, uint64_t* il_small, uint64_t* il_large, uint32_t intersect_filter_type, uint32_t overlap_type, double overlap_val, uint32_t species, uint64_t chrom_mask) {
  int64_t* marker_pos_arr = (int64_t*)wkspace_base;
  FILE* new_mapfile = NULL;
  uintptr_t raw_marker_ct = 0;
#ifdef __LP64__
  uint32_t distinct_marker_countdown = 2147483646;
#endif
  uint32_t req_fields = 3;
  uint32_t filter_seglen = min_seglen || (max_seglen < 4294967295U);
  uint32_t filter_score = (min_score > -INFINITY) || (max_score < INFINITY);
  uint32_t filter_sites = min_sites || (max_sites < 4294967295U);
  uintptr_t max_marker_ct;
  int32_t retval;
  char* bufptr;
  char* bufptr2;
  int64_t llii;
  uint64_t ullii;
  uintptr_t ulii;
  int32_t ii;
  logprintb();
  if (fopen_checked(&new_mapfile, new_mapname, "w")) {
    goto cnv_make_map_ret_OPEN_FAIL;
  }
  retval = cnv_first_nonheader_line(cnvfile);
  if (retval) {
    goto cnv_make_map_ret_1;
  }
  max_marker_ct = wkspace_left / sizeof(int64_t);
  // allow SCORE/SITES to be missing if they aren't being filtered on
  if (filter_sites) {
    req_fields = 5;
  } else if (filter_score) {
    req_fields = 4;
  }
  do {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("\nError: Pathologically long line in .cnv file.\n");
      goto cnv_make_map_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (!is_eoln_kns(*bufptr)) {
      // FID, IID, CHR, BP1, BP2, TYPE, SCORE, SITES
      bufptr = next_item_mult(bufptr, 2);
      bufptr2 = next_item_mult(bufptr, req_fields);
      if (no_more_items_kns(bufptr2)) {
	logprint("\nError: Fewer items than expected in .cnv line.\n");
	goto cnv_make_map_ret_INVALID_FORMAT;
      }
      ii = marker_code(species, bufptr);
      if (ii == -1) {
	logprint("\nError: Invalid chromosome code in .cnv file.\n");
        goto cnv_make_map_ret_INVALID_FORMAT;
      }
      ullii = ((uint64_t)((uint32_t)ii)) << 32;
      bufptr = next_item(bufptr);
      ii = atoi(bufptr);
      if (ii < 0) {
        logprint("\nError: Negative bp position in .cnv file.\n");
      }
      llii = (int64_t)(ullii | ((uint64_t)((uint32_t)ii)));
      bufptr = next_item(bufptr);
      ii = atoi(bufptr);
      if (ii < 0) {
        logprint("\nError: Negative bp position in .cnv file.\n");
      }
      ullii |= (uint64_t)((uint32_t)ii);
      if (filter_seglen) {
	bufptr2 = next_item(bufptr);
	if (0) {
	  continue;
	}
      }
      if (filter_score) {
	bufptr2 = next_item_mult(bufptr, 2);
	if (0) {
	  continue;
	}
      }
      if (filter_sites) {
	bufptr2 = next_item_mult(bufptr, 3);
	if (0) {
	  continue;
	}
      }
      if (raw_marker_ct + 1 >= max_marker_ct) {
        goto cnv_make_map_ret_NOMEM;
      }
      marker_pos_arr[raw_marker_ct++] = llii;
      marker_pos_arr[raw_marker_ct++] = (int64_t)ullii;
    }
  } while (fgets(tbuf, MAXLINELEN, cnvfile));
  if (!raw_marker_ct) {
    if (cnv_calc_type) {
      logprint("\nError: No markers.\n");
    } else {
      logprint("\nError: No markers after filtering.\n");
    }
    goto cnv_make_map_ret_INVALID_FORMAT;
  }
#ifdef __cplusplus
  std::sort(marker_pos_arr, &(marker_pos_arr[raw_marker_ct]));
#else
  qsort(marker_pos_arr, raw_marker_ct, sizeof(int64_t), llcmp);
#endif
  llii = marker_pos_arr[0];
  if (cnv_make_map_write(new_mapfile, (uint64_t)llii)) {
    goto cnv_make_map_ret_WRITE_FAIL;
  }
  for (ulii = 1; ulii < raw_marker_ct; ulii++) {
    if (marker_pos_arr[ulii] != llii) {
#ifdef __LP64__
      if (!(--distinct_marker_countdown)) {
	logprint("\nError: More than 2147483647 distinct positions.\n");
	goto cnv_make_map_ret_INVALID_FORMAT;
      }
#endif
      llii = marker_pos_arr[ulii];
      if (cnv_make_map_write(new_mapfile, (uint64_t)llii)) {
	goto cnv_make_map_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&new_mapfile)) {
    goto cnv_make_map_ret_WRITE_FAIL;
  }
  logprint(" done.\n");
  while (0) {
  cnv_make_map_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  cnv_make_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
    /*
  cnv_make_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
    */
  cnv_make_map_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  cnv_make_map_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 cnv_make_map_ret_1:
  fclose_cond(new_mapfile);
  wkspace_reset((unsigned char*)marker_pos_arr);
  return retval;
}

int32_t wdist_cnv(char* outname, char* outname_end, char* cnvname, char* mapname, char* famname, char* phenoname, uint32_t cnv_calc_type, uint32_t min_seglen, uint32_t max_seglen, double min_score, double max_score, uint32_t min_sites, uint32_t max_sites, uint32_t intersect_filter_type, char* intersect_filter_fname, char* subset_fname, uint32_t overlap_type, double overlap_val, uint32_t freq_type, uint32_t freq_val, uint32_t test_window, uint32_t segment_modifier, char* segment_spanning_fname, uint32_t indiv_mperms, uint32_t test_mperms, uint32_t test_region_mperms, uint32_t enrichment_test_mperms, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* cnvfile = NULL;
  FILE* famfile = NULL;
  FILE* mapfile = NULL;
  FILE* outfile = NULL;
  int32_t retval = 0;
  char* subset_list = NULL;
  uintptr_t subset_ct = 0;
  uintptr_t max_subset_name_len = 0;
  uintptr_t topsize = 0;
  uintptr_t* il_chrom_start_small = NULL;
  uintptr_t* il_chrom_start_large = NULL;
  uint32_t* il_chrom_max_width_small = NULL;
  uint32_t* il_chrom_max_width_large = NULL;
  uint64_t* il_small = NULL; // high-order 32 bits = 2x center pos,
                             // low-order 32 bits = interval end - start
  uint64_t* il_large = NULL;
  char* sptr;
  uint32_t uii;
  if (fopen_checked(&cnvfile, cnvname, "r")) {
    goto wdist_cnv_ret_OPEN_FAIL;
  }
  if (cnv_calc_type & (~CNV_MAKE_MAP)) {
    if (fopen_checked(&famfile, famname, "r")) {
      goto wdist_cnv_ret_OPEN_FAIL;
    }
  }
  tbuf[MAXLINELEN - 1] = ' ';
  if (intersect_filter_type) {
    if (wkspace_alloc_ul_checked(&il_chrom_start_small, (MAX_POSSIBLE_CHROM + 1) * sizeof(intptr_t)) ||
        wkspace_alloc_ul_checked(&il_chrom_start_large, (MAX_POSSIBLE_CHROM + 1) * sizeof(intptr_t)) ||
        wkspace_alloc_ui_checked(&il_chrom_max_width_small, MAX_POSSIBLE_CHROM * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&il_chrom_max_width_large, MAX_POSSIBLE_CHROM * sizeof(int32_t))) {
      goto wdist_cnv_ret_NOMEM;
    }
    if (subset_fname) {
      retval = cnv_subset_load(subset_fname, &subset_list, &subset_ct, &max_subset_name_len, &topsize);
      if (retval) {
	goto wdist_cnv_ret_1;
      }
    }
    retval = cnv_intersect_load(intersect_filter_type, intersect_filter_fname, subset_list, subset_ct, max_subset_name_len, topsize, il_chrom_start_small, il_chrom_start_large, il_chrom_max_width_small, il_chrom_max_width_large, &il_small, &il_large, chrom_info_ptr->species, chrom_info_ptr->chrom_mask);
    if (retval) {
      goto wdist_cnv_ret_1;
    }
  }
  if (!(cnv_calc_type & CNV_MAKE_MAP)) {
    sptr = (char*)memchr(mapname, 0, FNAMESIZE);
    if ((mapname[0] == '\0') || (!filename_exists(mapname, sptr, ""))) {
      if (mapname[0] == '\0') {
        uii = strlen(cnvname);
        if ((uii < 5) || (cnvname[uii - 4] != '.') || (!match_upper_nt(&(cnvname[uii - 3]), "CNV", 3))) {
	  logprint("Error: No .cnv.map filename specified, and .cnv filename does not unambiguously\nspecify how an autogenerated file should be named.  Use --cnv-make-map + --out.\n");
	  goto wdist_cnv_ret_INVALID_CMDLINE;
	}
	memcpy(mapname, cnvname, uii);
	memcpy(&(mapname[uii]), ".map", 5);
	sptr = &(mapname[uii + 4]);
	if (filename_exists(mapname, sptr, "")) {
	  sprintf(logbuf, "Error: No .cnv.map filename specified, and natural autogeneration target\n(%s) already exists.\n", mapname);
	  logprintb();
	  goto wdist_cnv_ret_INVALID_CMDLINE;
	}
      }
      sprintf(logbuf, "Autogenerating missing %s...", mapname);
      retval = cnv_make_map(cnvfile, mapname, 0, 0, 4294967295U, -INFINITY, INFINITY, 0, 4294967295U, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0.0, chrom_info_ptr->species, ~0LLU);
      if (retval) {
	goto wdist_cnv_ret_1;
      }
    }
  } else {
    memcpy(outname_end, ".cnv.map", 9);
    sprintf(logbuf, "Generating %s...", outname);
    retval = cnv_make_map(cnvfile, outname, cnv_calc_type, min_seglen, max_seglen, min_score, max_score, min_sites, max_sites, il_chrom_start_small, il_chrom_start_large, il_chrom_max_width_small, il_chrom_max_width_large, il_small, il_large, intersect_filter_type, overlap_type, overlap_val, chrom_info_ptr->species, chrom_info_ptr->chrom_mask);
    if (retval || (!(cnv_calc_type & (CNV_MAKE_MAP | CNV_DEL | CNV_DUP)))) {
      goto wdist_cnv_ret_1;
    }
  }
  if (fopen_checked(&mapfile, mapname, "r")) {
    goto wdist_cnv_ret_OPEN_FAIL;
  }

  while (0) {
  wdist_cnv_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  wdist_cnv_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
    /*
  wdist_cnv_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  wdist_cnv_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
    */
  wdist_cnv_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 wdist_cnv_ret_1:
  if (topsize) {
    wkspace_left += topsize;
  }
  fclose_cond(cnvfile);
  fclose_cond(famfile);
  fclose_cond(mapfile);
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t wdist_gvar(char* outname, char* outname_end, char* gvarname, char* mapname, char* famname) {
  logprint("Error: Common CNP analysis not yet supported.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
}
