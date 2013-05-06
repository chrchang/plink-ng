#include "wdist_common.h"

int32_t cnv_subset_load(char* subset_fname, char** subset_list_ptr, uintptr_t* subset_ct_ptr, uintptr_t* max_subset_name_len_ptr) {
  FILE* subset_file = NULL;
  uintptr_t subset_ct = 0;
  uintptr_t max_subset_name_len = 1;
  int32_t retval = 0;
  char* bufptr;
  char* slptr;
  uintptr_t ulii;
  // Two-pass load:
  // 1. Determine subset_ct and max_subset_name_len.
  // 2. Allocate temporary block and load the strings.
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
  if (wkspace_alloc_c_checked(subset_list_ptr, subset_ct * max_subset_name_len)) {
    goto cnv_subset_load_ret_NOMEM;
  }
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
  if (intersect_filter_type & CNV_INTERSECT) {
    return (char*)"intersect";
  } else if (intersect_filter_type & CNV_EXCLUDE) {
    return (char*)"exclude";
  } else {
    return (char*)"count";
  }
}

// maximum 25
#define SMALL_INTERVAL_BITS 18
#define SMALL_INTERVAL_MAX_SIZE ((1 << SMALL_INTERVAL_BITS) - 1)

int32_t cnv_intersect_load(uint32_t intersect_filter_type, char* intersect_filter_fname, char* subset_list, uintptr_t subset_ct, uintptr_t max_subset_name_len, uintptr_t* il_chrom_start_small, uintptr_t* il_chrom_start_large, uint32_t* il_chrom_max_width_small, uint32_t* il_chrom_max_width_large, uint64_t** il_small_ptr, uint64_t** il_large_ptr, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t species, uint64_t chrom_mask, uintptr_t* topsize_ptr) {
  // We store intervals in sorted order, with the center of each interval in
  // the high-order bits, and the size (without adding 1) in the low-order
  // bits.  (Chromosome beginnings and endings are stored in small external
  // arrays.)  We also track each chromosome's largest interval size.
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
  unsigned char* tmp_il_large_chroms = wkspace_base; // grows up
  uint64_t* tmp_il_small = (uint64_t*)(&(tmp_il_large_chroms[(max_interval_ct + 7) & (~(7 * ONELU))])); // grows up
  uint64_t* il_large = (uint64_t*)(&(wkspace_base[wkspace_left])); // grows down
  uint32_t reverse_warning_given = 0;
  int32_t retval = 0;
  uint64_t* il_small;
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
  uint32_t cur_chrom;
  uint32_t uii;
  unsigned char ucc;
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
      uljj = ((uint32_t)ii);
      if (!((chrom_mask >> uljj) & 1)) {
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
      if ((marker_pos_start > ii) || ((marker_pos_end != -1) && (marker_pos_end < jj))) {
	continue;
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
      kk = ii - jj;
      if (kk > SMALL_INTERVAL_MAX_SIZE) {
	il_large--;
	*il_large = (((uint64_t)(((uint32_t)ii) + ((uint32_t)jj))) << 31) | ((uint64_t)((uint32_t)kk));
	tmp_il_large_chroms[large_interval_ct++] = uljj;
      } else {
	tmp_il_small[small_interval_ct++] = (((uint64_t)uljj) << (32 + SMALL_INTERVAL_BITS)) | (((uint64_t)(((uint32_t)ii) + ((uint32_t)jj))) << SMALL_INTERVAL_BITS) | ((uint64_t)((uint32_t)kk));
      }
    }
  }
  if (!feof(intersect_file)) {
    goto cnv_intersect_load_ret_READ_FAIL;
  }
  *il_large_ptr = il_large;
  il_small = &(il_large[-((intptr_t)small_interval_ct)]);
  *il_small_ptr = il_small;
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
    std::sort((int64_t*)tmp_il_small, (int64_t*)(&(tmp_il_small[small_interval_ct])));
#else
    qsort((int64_t*)tmp_il_small, small_interval_ct, sizeof(int64_t), llcmp);
#endif
    il_chrom_start_small[0] = 0;
    cur_chrom = 0;
    max_width = 0;
    ulii = small_interval_ct;
    do {
      ulii--;
      ullii = tmp_il_small[ulii];
      uii = (uint32_t)(ullii >> (SMALL_INTERVAL_BITS + 32));
      if (uii > cur_chrom) {
	il_chrom_max_width_small[cur_chrom] = max_width;
	do {
	  il_chrom_start_small[++cur_chrom] = ulii;
	} while (cur_chrom < uii);
	max_width = 0;
      }
      cur_width = ullii & (SMALL_INTERVAL_MAX_SIZE * 1LLU);
      if (cur_width > max_width) {
	max_width = cur_width;
      }
      il_small[ulii] = ((ullii >> SMALL_INTERVAL_BITS) << 32) | ((uint64_t)cur_width);
    } while (ulii);
    il_chrom_max_width_small[cur_chrom] = max_width;
    do {
      il_chrom_start_small[++cur_chrom] = small_interval_ct;
    } while (cur_chrom < MAX_POSSIBLE_CHROM);
  } else {
    fill_ulong_zero(il_chrom_start_small, MAX_POSSIBLE_CHROM + 1);
  }
  if (large_interval_ct) {
    if (large_interval_ct > 1) {
      uljj = large_interval_ct / 2;
      ulkk = large_interval_ct - 1;
      for (ulii = 0; ulii < uljj; ulii++) {
        ucc = tmp_il_large_chroms[ulii];
	tmp_il_large_chroms[ulii] = tmp_il_large_chroms[ulkk];
	tmp_il_large_chroms[ulkk--] = ucc;
      }
    }
    qsort_ext((char*)tmp_il_large_chroms, large_interval_ct, sizeof(char), char_cmp_deref, (char*)il_large, sizeof(int64_t));
    il_chrom_start_large[0] = 0;
    cur_chrom = 0;
    for (ulii = 0; ulii < large_interval_ct; ulii++) {
      uii = tmp_il_large_chroms[ulii];
      if (uii > cur_chrom) {
	do {
	  il_chrom_start_large[++cur_chrom] = ulii;
	} while (cur_chrom < uii);
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
	std::sort((int64_t*)(&(il_large[ulii])), (int64_t*)(&(il_large[uljj])));
#else
        qsort((int64_t*)(&(il_large[ulii])), uljj - ulii, sizeof(int64_t), llcmp);
#endif
	ulii = uljj;
      }
    }
    ulii = 0;
    for (cur_chrom = 0; cur_chrom < MAX_POSSIBLE_CHROM; cur_chrom++) {
      uljj = il_chrom_start_large[cur_chrom + 1];
      max_width = 0;
      while (ulii < uljj) {
	ullii = il_large[ulii];
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
  *topsize_ptr = CACHELINE * ((small_interval_ct + large_interval_ct + CACHELINE_INT64 - 1) / CACHELINE_INT64);
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

uint32_t is_cnv_overlap_one_size(uint32_t start_pos, uint32_t end_pos, uint32_t overlap_type, double overlap_val, uint32_t max_width, uint64_t* interval_list, uintptr_t interval_list_len) {
  // earliest possible overlap: center at start_pos - (max_width / 2)
  // latest possible overlap: center at end_pos + (max_width / 2)
  uint64_t twice_end_pos = (uint64_t)(end_pos * 2);
  uint64_t ullii = ((uint64_t)start_pos) << 33;
  uint64_t ulljj = (uint64_t)max_width;
  uint64_t ullkk = ulljj << 32;
  double denominator = (double)(1 + end_pos - start_pos);
  uint32_t region_start;
  uint32_t region_end;
  uint32_t smaller_start;
  uint32_t larger_start;
  uint32_t smaller_end;
  uint32_t larger_end;
  uintptr_t cur_idx;
  uintptr_t last_idx;
  double numerator;
  if (!interval_list_len) {
    return 0;
  }
  if (ullii <= ullkk) {
    cur_idx = 0;
  } else {
    cur_idx = uint64arr_greater_than(interval_list, interval_list_len, (ullii - ullkk) | ulljj);
  }
  if ((twice_end_pos + ulljj > 0xffffffffLLU) || (cur_idx == interval_list_len)) {
    last_idx = interval_list_len;
  } else {
    last_idx = cur_idx + uint64arr_greater_than(&(interval_list[cur_idx]), interval_list_len - cur_idx, ((twice_end_pos << 32) + ullkk) | 0xffffffffLLU);
  }
  // printf("%u %u %u %g %u %u %u %lu\n", start_pos, end_pos, overlap_type, overlap_val, small_max_width, (uint32_t)(il_small[0] >> 32), (uint32_t)il_small[0], il_small_len);
  while (cur_idx < last_idx) {
    ullii = interval_list[cur_idx++];
    region_end = ((uint32_t)((ullii >> 32) + ullii)) / 2;
    if (region_end >= start_pos) {
      region_start = ((uint32_t)((ullii >> 32) - ullii)) / 2;
      if (region_start <= end_pos) {
	if (!overlap_type) {
	  return 1;
	} else if (overlap_type == CNV_DISRUPT) {
	  if (((region_start < start_pos) && (region_end <= end_pos)) || ((region_start >= start_pos) && (region_end > end_pos))) {
	    return 1;
	  }
	} else {
	  if (region_start < start_pos) {
	    smaller_start = region_start;
	    larger_start = start_pos;
	  } else {
	    smaller_start = start_pos;
	    larger_start = region_start;
	  }
	  if (region_end < end_pos) {
	    smaller_end = region_end;
	    larger_end = end_pos;
	  } else {
	    smaller_end = end_pos;
	    larger_end = region_end;
	  }
	  numerator = (double)(1 + smaller_end - larger_start);
	  if (overlap_type == CNV_OVERLAP_REGION) {
	    denominator = (double)(1 + region_end - region_start);
	  } else if (overlap_type == CNV_OVERLAP_UNION) {
	    denominator = (double)(1 + larger_end - smaller_start);
	  }
	  if (denominator * overlap_val <= numerator) {
	    return 1;
	  }
	}
      }
    }
  }
  return 0;
}

uint32_t is_cnv_overlap(uint32_t start_pos, uint32_t end_pos, uint32_t overlap_type, double overlap_val, uint32_t small_max_width, uint32_t large_max_width, uint64_t* il_small, uintptr_t il_small_len, uint64_t* il_large, uintptr_t il_large_len) {
  if (is_cnv_overlap_one_size(start_pos, end_pos, overlap_type, overlap_val, small_max_width, il_small, il_small_len)) {
    return 1;
  }
  return is_cnv_overlap_one_size(start_pos, end_pos, overlap_type, overlap_val, large_max_width, il_large, il_large_len);
}

int32_t cnv_make_map_write(FILE* new_mapfile, uint64_t ullii) {
  uint32_t chrom_idx = (uint32_t)(ullii >> 32);
  uint32_t bp_pos = (uint32_t)ullii;
  char* wptr = uint32_write(tbuf, chrom_idx);
  wptr = memcpya(wptr, "\tp", 2);
  wptr = uint32_writex(wptr, chrom_idx, '-');
  wptr = uint32_write(wptr, bp_pos);
  wptr = memcpyl3a(wptr, "\t0\t");
  wptr = uint32_writex(wptr, bp_pos, '\n');
  return fwrite_checked(tbuf, wptr - tbuf, new_mapfile);
}

int32_t cnv_make_map(FILE* cnvfile, char* new_mapname, uint32_t cnv_calc_type, uint32_t min_seglen, uint32_t max_seglen, double min_score, double max_score, uint32_t min_sites, uint32_t max_sites, uintptr_t* il_chrom_start_small, uintptr_t* il_chrom_start_large, uint32_t* il_chrom_max_width_small, uint32_t* il_chrom_max_width_large, uint64_t* il_small, uint64_t* il_large, uint32_t intersect_filter_type, uint32_t overlap_type, double overlap_val, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t species, uint64_t chrom_mask) {
  int64_t* marker_pos_arr = (int64_t*)wkspace_base;
  FILE* new_mapfile = NULL;
  uintptr_t raw_marker_ct = 0;
#ifdef __LP64__
  uint32_t distinct_marker_countdown = 0x7ffffffe;
#endif
  uint32_t req_fields = 3;
  uint32_t filter_seglen = min_seglen || (max_seglen < 0xffffffffU);
  uint32_t cnv_del = cnv_calc_type & CNV_DEL;
  uint32_t filter_score = (min_score > -INFINITY) || (max_score < INFINITY);
  uint32_t filter_sites = min_sites || (max_sites < 0xffffffffU);
  uint32_t make_map_long = cnv_calc_type & CNV_MAKE_MAP_LONG;
  uintptr_t max_marker_ct;
  int32_t retval;
  char* bufptr;
  char* bufptr2;
  int64_t llii;
  uint64_t ullii;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t chrom_idx;
  int32_t ii;
  int32_t jj;
  double dxx;
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
      chrom_idx = ii;
      if (!((chrom_mask >> chrom_idx) & 1)) {
	continue;
      }
      ullii = ((uint64_t)chrom_idx) << 32;
      bufptr = next_item(bufptr);
      ii = atoi(bufptr);
      if (ii < 0) {
        logprint("\nError: Negative bp position in .cnv file.\n");
	goto cnv_make_map_ret_INVALID_FORMAT;
      }
      bufptr = next_item(bufptr);
      jj = atoi(bufptr);
      if (jj < 0) {
        logprint("\nError: Negative bp position in .cnv file.\n");
	goto cnv_make_map_ret_INVALID_FORMAT;
      } else if (jj < ii) {
	logprint("\nError: Segment end position smaller than segment start in .cnv file.\n");
	goto cnv_make_map_ret_INVALID_FORMAT;
      } else if (jj > 0x7ffffffe) {
	logprint("\nError: Excessively large bp position in .cnv file.\n");
	goto cnv_make_map_ret_INVALID_FORMAT;
      }
      if ((marker_pos_start > ii) || ((marker_pos_end != -1) && (marker_pos_end < jj))) {
	continue;
      }
      llii = (int64_t)(ullii | ((uint64_t)((uint32_t)ii)));
      ullii |= (uint64_t)((uint32_t)jj);
      if (filter_seglen) {
	uii = 1 + ((uint32_t)(jj - ii));
	if ((uii < min_seglen) || (uii > max_seglen)) {
	  continue;
	}
      }
      if (cnv_calc_type & (CNV_DEL | CNV_DUP)) {
	bufptr2 = next_item(bufptr);
	if (!atoiz(bufptr2, &ii)) {
	  logprint("\nError: Invalid variant copy count in .cnv file.\n");
	  goto cnv_make_map_ret_INVALID_FORMAT;
	}
	if (cnv_del) {
	  if (ii > 1) {
	    continue;
	  }
	} else if (ii < 3) {
	  continue;
	}
      }
      if (filter_score) {
	bufptr2 = next_item_mult(bufptr, 2);
	if (sscanf(bufptr2, "%lg", &dxx) != 1) {
	  logprint("\nError: Invalid confidence score in .cnv file.\n");
	  goto cnv_make_map_ret_INVALID_FORMAT;
	}
	if ((dxx < min_score) || (dxx > max_score)) {
	  continue;
	}
      }
      if (filter_sites) {
	bufptr2 = next_item_mult(bufptr, 3);
	ii = atoi(bufptr2);
	if (ii < 1) {
	  logprint("\nError: Invalid probe count in .cnv file.\n");
	  goto cnv_make_map_ret_INVALID_FORMAT;
	}
	if ((((uint32_t)ii) < min_sites) || (((uint32_t)ii) > max_sites)) {
	  continue;
	}
      }
      if (intersect_filter_type & (CNV_INTERSECT | CNV_EXCLUDE)) {
	ulii = il_chrom_start_small[chrom_idx];
	uljj = il_chrom_start_large[chrom_idx];
	if (is_cnv_overlap((uint32_t)((uint64_t)llii), (uint32_t)ullii, overlap_type, overlap_val, il_chrom_max_width_small[chrom_idx], il_chrom_max_width_large[chrom_idx], &(il_small[ulii]), il_chrom_start_small[chrom_idx + 1] - ulii, &(il_large[uljj]), il_chrom_start_large[chrom_idx + 1] - uljj)) {
          if (intersect_filter_type & CNV_EXCLUDE) {
	    continue;
	  }
	} else if (intersect_filter_type & CNV_INTERSECT) {
	  continue;
	}
      }
      if (raw_marker_ct + 2 >= max_marker_ct) {
        goto cnv_make_map_ret_NOMEM;
      }
      marker_pos_arr[raw_marker_ct++] = llii;
      if (make_map_long && (((uint64_t)llii) != ullii)) {
        marker_pos_arr[raw_marker_ct++] = (int64_t)ullii;
      }
      marker_pos_arr[raw_marker_ct++] = 1 + (int64_t)ullii;
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

int32_t wdist_cnv(char* outname, char* outname_end, char* cnvname, char* mapname, char* famname, char* phenoname, uint32_t cnv_calc_type, uint32_t min_seglen, uint32_t max_seglen, double min_score, double max_score, uint32_t min_sites, uint32_t max_sites, uint32_t intersect_filter_type, char* intersect_filter_fname, char* subset_fname, uint32_t overlap_type, double overlap_val, uint32_t freq_type, uint32_t freq_val, double freq_val2, uint32_t test_window, uint32_t segment_modifier, char* segment_spanning_fname, uint32_t indiv_mperms, uint32_t test_mperms, uint32_t test_region_mperms, uint32_t enrichment_test_mperms, char* markername_from, char* markername_to, int32_t marker_pos_start, int32_t marker_pos_end, Chrom_info* chrom_info_ptr) {
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
  uintptr_t il_chrom_start_small[MAX_POSSIBLE_CHROM + 1];
  uintptr_t il_chrom_start_large[MAX_POSSIBLE_CHROM + 1];
  uint32_t il_chrom_max_width_small[MAX_POSSIBLE_CHROM];
  uint32_t il_chrom_max_width_large[MAX_POSSIBLE_CHROM];
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
    if (subset_fname) {
      retval = cnv_subset_load(subset_fname, &subset_list, &subset_ct, &max_subset_name_len);
      if (retval) {
	goto wdist_cnv_ret_1;
      }
    }
    retval = cnv_intersect_load(intersect_filter_type, intersect_filter_fname, subset_list, subset_ct, max_subset_name_len, il_chrom_start_small, il_chrom_start_large, il_chrom_max_width_small, il_chrom_max_width_large, &il_small, &il_large, marker_pos_start, marker_pos_end, chrom_info_ptr->species, chrom_info_ptr->chrom_mask, &topsize);
    if (retval) {
      goto wdist_cnv_ret_1;
    }
    wkspace_reset(wkspace_mark);
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
      if (markername_from || markername_to) {
	logprint("Error: --from/--to cannot be used with .cnv.map autogeneration.\n");
	goto wdist_cnv_ret_INVALID_CMDLINE;
      }
      sprintf(logbuf, "Autogenerating missing %s...", mapname);
      retval = cnv_make_map(cnvfile, mapname, 0, 0, 0xffffffffU, -INFINITY, INFINITY, 0, 0xffffffffU, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0.0, -1, -1, chrom_info_ptr->species, ~0LLU);
      if (retval) {
	goto wdist_cnv_ret_1;
      }
    }
  } else {
    memcpy(outname_end, ".cnv.map", 9);
    sprintf(logbuf, "Generating %s...", outname);
    retval = cnv_make_map(cnvfile, outname, cnv_calc_type, min_seglen, max_seglen, min_score, max_score, min_sites, max_sites, il_chrom_start_small, il_chrom_start_large, il_chrom_max_width_small, il_chrom_max_width_large, il_small, il_large, intersect_filter_type, overlap_type, overlap_val, marker_pos_start, marker_pos_end, chrom_info_ptr->species, chrom_info_ptr->chrom_mask);
    if (retval || (!(cnv_calc_type & (CNV_MAKE_MAP | CNV_DEL | CNV_DUP)))) {
      goto wdist_cnv_ret_1;
    }
  }
  if (fopen_checked(&mapfile, mapname, "r")) {
    goto wdist_cnv_ret_OPEN_FAIL;
  }

  while (0) {
    /*
  wdist_cnv_ret_NOMEM:
    retval = RET_NOMEM;
    break;
    */
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
