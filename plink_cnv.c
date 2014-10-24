#include "plink_common.h"

#ifndef HIGH_MAX_CHROM
int32_t cnv_subset_load(char* subset_fname, char** subset_list_ptr, uintptr_t* subset_ct_ptr, uintptr_t* max_subset_name_len_ptr) {
  FILE* subset_file = NULL;
  uintptr_t subset_ct = 0;
  uintptr_t max_subset_name_len = 0;
  int32_t retval = open_and_size_string_list(subset_fname, &subset_file, &subset_ct, &max_subset_name_len);
  if (retval) {
    goto cnv_subset_load_ret_1;
  }
  if (!subset_ct) {
    logprint("Error: Empty --cnv-subset file.\n");
    goto cnv_subset_load_ret_INVALID_FORMAT;
  }
#ifndef __LP64__
  if (((uint64_t)subset_ct) * max_subset_name_len > 0x7fffffffLLU) {
    goto cnv_subset_load_ret_NOMEM;
  }
#endif
  if (wkspace_alloc_c_checked(subset_list_ptr, subset_ct * max_subset_name_len)) {
    goto cnv_subset_load_ret_NOMEM;
  }
  retval = load_string_list(&subset_file, max_subset_name_len, *subset_list_ptr);
  qsort(*subset_list_ptr, subset_ct, max_subset_name_len, strcmp_casted);
  *subset_ct_ptr = subset_ct;
  *max_subset_name_len_ptr = max_subset_name_len;
  while (0) {
  cnv_subset_load_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  cnv_subset_load_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 cnv_subset_load_ret_1:
  fclose_cond(subset_file);
  return retval;
}

const char* cnv_intersect_filter_type_to_str(uint32_t intersect_filter_type) {
  if (intersect_filter_type & CNV_INTERSECT) {
    return "--cnv-intersect file";
  } else if (intersect_filter_type & CNV_EXCLUDE) {
    return "--cnv-exclude file";
  } else {
    return "--cnv-count file";
  }
}

// log_2(MAX_POSSIBLE_CHROM) + SMALL_INTERVAL_BITS must not exceed 32
#define SMALL_INTERVAL_BITS 18
#define SMALL_INTERVAL_MAX_SIZE ((1 << SMALL_INTERVAL_BITS) - 1)

int32_t cnv_intersect_load(uint32_t intersect_filter_type, char* intersect_filter_fname, char* subset_list, uintptr_t subset_ct, uintptr_t max_subset_name_len, uintptr_t* il_chrom_start_small, uintptr_t* il_chrom_start_large, uint32_t* il_chrom_max_width_small, uint32_t* il_chrom_max_width_large, uint64_t** il_small_ptr, uint64_t** il_large_ptr, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t allow_extra_chroms, Chrom_info* chrom_info_ptr, uintptr_t* topsize_ptr) {
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
  uintptr_t reverse_warning_ct = 0;
  uintptr_t line_idx = 0;
  // Storage format before sorting:
  // Small: high 16 bits = chrom, middle 32 bits = center pos * 2,
  //        bottom 16 bits = size
  // Large: top bit = zero, next 32 bits = center pos * 2, bottom 31 = size
  //        [chrom information stored separately, initially in reverse order]
  unsigned char* tmp_il_large_chroms = wkspace_base; // grows up
  uint64_t* tmp_il_small = (uint64_t*)(&(tmp_il_large_chroms[(max_interval_ct + 7) & (~(7 * ONELU))])); // grows up
  uint64_t* il_large = (uint64_t*)(&(wkspace_base[wkspace_left])); // grows down
  uintptr_t* chrom_mask = chrom_info_ptr->chrom_mask;
  const char* cift_str = cnv_intersect_filter_type_to_str(intersect_filter_type);
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
  uint32_t chrom_code_end;
  uint32_t cur_chrom;
  uint32_t uii;
  unsigned char ucc;
  if (fopen_checked(&intersect_file, intersect_filter_fname, "r")) {
    goto cnv_intersect_load_ret_OPEN_FAIL;
  }
  while (fgets(tbuf, MAXLINELEN, intersect_file)) {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, cift_str);
      goto cnv_intersect_load_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (!is_eoln_kns(*bufptr)) {
      // CHR, BP1, BP2, subset name
      bufptr2 = next_token_mult(bufptr, 2);
      if (no_more_tokens_kns(bufptr2)) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of %s file has fewer tokens than expected.\n", line_idx, cift_str);
	goto cnv_intersect_load_ret_INVALID_FORMAT_2;
      }
      ii = get_chrom_code(chrom_info_ptr, bufptr);
      if (ii < 0) {
	if ((!allow_extra_chroms) || (ii == -1)) {
	  sprintf(logbuf, "Error: Invalid chromosome code on line %" PRIuPTR " of %s.\n", line_idx, cift_str);
	  goto cnv_intersect_load_ret_INVALID_FORMAT_2;
	}
        retval = resolve_or_add_chrom_name(chrom_info_ptr, bufptr, &ii, line_idx, cift_str);
	if (retval) {
	  goto cnv_intersect_load_ret_1;
	}
      }
      uljj = ((uint32_t)ii);
      if (!IS_SET(chrom_mask, uljj)) {
	continue;
      }
      bufptr = next_token(bufptr);
      if (scan_uint_defcap(bufptr, (uint32_t*)&jj)) {
	sprintf(logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, cift_str);
	goto cnv_intersect_load_ret_INVALID_FORMAT_2;
      }
      if (scan_uint_defcap(bufptr2, (uint32_t*)&ii)) {
	sprintf(logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, cift_str);
	goto cnv_intersect_load_ret_INVALID_FORMAT_2;
      }
      if (ii < jj) {
	if (!reverse_warning_ct) {
	  LOGPRINTF("Warning: End of range before start of range on line %" PRIuPTR " of\n%s.\n", line_idx, cift_str);
	}
	reverse_warning_ct++;
	kk = ii;
	ii = jj;
	jj = kk;
      }
      if ((marker_pos_start > ii) || ((marker_pos_end != -1) && (marker_pos_end < jj))) {
	continue;
      }
      if (subset_ct) {
	bufptr = next_token(bufptr2);
	if (no_more_tokens_kns(bufptr)) {
	  continue;
	}
	if (bsearch_str(bufptr, strlen_se(bufptr), subset_list, max_subset_name_len, subset_ct) == -1) {
	  continue;
	}
      }
      if (small_interval_ct + large_interval_ct == max_interval_ct) {
        goto cnv_intersect_load_ret_NOMEM;
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
  if (reverse_warning_ct > 1) {
    LOGPRINTF("(%" PRIuPTR " subsequent line%s with [end of range] < [start of range].)\n", reverse_warning_ct - 1, (reverse_warning_ct == 2)? "" : "s");
  }
  *il_large_ptr = il_large;
  il_small = &(il_large[-((intptr_t)small_interval_ct)]);
  *il_small_ptr = il_small;
  chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  if ((!small_interval_ct) && (!large_interval_ct)) {
    if (subset_ct) {
      fill_ulong_zero(il_chrom_start_small, chrom_code_end + 1);
      fill_ulong_zero(il_chrom_start_large, chrom_code_end + 1);
      logprint("Warning: All intervals filtered out by --cnv-subset.\n");
      goto cnv_intersect_load_ret_1;
    }
    sprintf(logbuf, "Error: Empty %s.\n", cift_str);
    goto cnv_intersect_load_ret_INVALID_FORMAT_2;
  }
  if (small_interval_ct) {
#ifdef __cplusplus
    std::sort((int64_t*)tmp_il_small, (int64_t*)(&(tmp_il_small[small_interval_ct])));
#else
    qsort((int64_t*)tmp_il_small, small_interval_ct, sizeof(int64_t), llcmp);
#endif
    il_chrom_start_small[chrom_code_end] = small_interval_ct;
    cur_chrom = chrom_code_end - 1;
    max_width = 0;
    ulii = small_interval_ct;
    do {
      ulii--;
      ullii = tmp_il_small[ulii];
      uii = (uint32_t)(ullii >> (SMALL_INTERVAL_BITS + 32));
      if (uii < cur_chrom) {
	il_chrom_max_width_small[cur_chrom] = max_width;
	do {
	  il_chrom_start_small[cur_chrom--] = ulii + 1;
	} while (cur_chrom > uii);
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
      il_chrom_start_small[cur_chrom] = 0;
    } while (cur_chrom--);
  } else {
    fill_ulong_zero(il_chrom_start_small, chrom_code_end + 1);
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
    if (qsort_ext((char*)tmp_il_large_chroms, large_interval_ct, sizeof(char), char_cmp_deref, (char*)il_large, sizeof(int64_t))) {
      goto cnv_intersect_load_ret_NOMEM;
    }
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
    } while (cur_chrom < chrom_code_end);
    ulii = il_chrom_start_large[0];
    for (cur_chrom = 0; cur_chrom < chrom_code_end; cur_chrom++) {
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
    for (cur_chrom = 0; cur_chrom < chrom_code_end; cur_chrom++) {
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
    fill_ulong_zero(il_chrom_start_large, chrom_code_end + 1);
  }
  *topsize_ptr = CACHELINE * ((small_interval_ct + large_interval_ct + CACHELINE_INT64 - 1) / CACHELINE_INT64);
  wkspace_left -= (*topsize_ptr);
  while (0) {
  cnv_intersect_load_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  cnv_intersect_load_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  cnv_intersect_load_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  cnv_intersect_load_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 cnv_intersect_load_ret_1:
  fclose_cond(intersect_file);
  return retval;
}

int32_t cnv_first_nonheader_line(FILE* cnvfile, uintptr_t* line_idx_ptr) {
  uintptr_t line_idx = 0;
  int32_t retval = 0;
  char* bufptr;
  rewind(cnvfile);
  // assumes tbuf[MAXLINELEN - 1] is initialized to space
  do {
    line_idx++;
    if (!fgets(tbuf, MAXLINELEN, cnvfile)) {
      goto cnv_first_nonheader_line_fgets_fail;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      LOGPRINTF("\nError: Line %" PRIuPTR " of .cnv file is pathologically long.\n", line_idx);
      goto cnv_first_nonheader_line_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
  } while (is_eoln_kns(*bufptr));
  if ((strlen_se(bufptr) == 3) && (!memcmp(bufptr, "FID", 3))) {
    line_idx++;
    if (!fgets(tbuf, MAXLINELEN, cnvfile)) {
      goto cnv_first_nonheader_line_fgets_fail;
    }
  }
  *line_idx_ptr = line_idx;
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

int32_t cnv_make_map_write(FILE* new_mapfile, Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uint32_t bp_pos, uintptr_t* max_marker_id_len_ptr) {
  char* wptr = chrom_name_write(tbuf, chrom_info_ptr, chrom_idx);
  char* wptr2 = memcpya(wptr, "\tp", 2);
  uintptr_t cur_marker_id_len;
  // this just needs to be an arbitrary unique name, so it's fine if we don't
  // use chrom_name_write() here
  wptr2 = uint32_writex(wptr2, chrom_idx, '-');
  wptr2 = uint32_write(wptr2, bp_pos);
  cur_marker_id_len = (uintptr_t)(wptr2 - wptr);
  if (cur_marker_id_len > (*max_marker_id_len_ptr)) {
    // includes an extra character at the start, to compensate for missing
    // terminal null
    *max_marker_id_len_ptr = cur_marker_id_len;
  }
  wptr2 = memcpyl3a(wptr2, "\t0\t");
  wptr2 = uint32_writex(wptr2, bp_pos, '\n');
  return fwrite_checked(tbuf, wptr2 - tbuf, new_mapfile);
}

int32_t cnv_make_map(FILE* cnvfile, char* new_mapname, uint32_t cnv_calc_type, uint32_t min_seglen, uint32_t max_seglen, double min_score, double max_score, uint32_t min_sites, uint32_t max_sites, uintptr_t* il_chrom_start_small, uintptr_t* il_chrom_start_large, uint32_t* il_chrom_max_width_small, uint32_t* il_chrom_max_width_large, uint64_t* il_small, uint64_t* il_large, uint32_t intersect_filter_type, uint32_t overlap_type, double overlap_val, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t allow_extra_chroms, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uintptr_t* max_marker_id_len_ptr, uint32_t* marker_chrom_start) {
  int64_t* marker_pos_arr = (int64_t*)wkspace_base;
  FILE* new_mapfile = NULL;
  uintptr_t raw_marker_ct = 0;
  uint32_t distinct_marker_ct = 1;
  uint32_t req_fields = 3;
  uint32_t filter_seglen = min_seglen || (max_seglen < 0xffffffffU);
  uint32_t cnv_del = cnv_calc_type & CNV_DEL;
  uint32_t filter_score = (min_score > -HUGE_DOUBLE) || (max_score < HUGE_DOUBLE);
  uint32_t filter_sites = min_sites || (max_sites < 0xffffffffU);
  uint32_t make_map_long = cnv_calc_type & CNV_MAKE_MAP_LONG;
  uint32_t is_autogen = (!il_chrom_start_small)? 1 : 0;
  uintptr_t* chrom_mask = chrom_info_ptr->chrom_mask;
  uint32_t chrom_code_end;
  uintptr_t max_marker_ct;
  int32_t retval;
  char* bufptr;
  char* bufptr2;
  int64_t llii;
  uint64_t ullii;
  uintptr_t line_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t chrom_idx;
  uint32_t seg_start;
  uint32_t seg_end;
  uint32_t uii;
  int32_t ii;
  double dxx;
  logprintb();
  if (fopen_checked(&new_mapfile, new_mapname, "w")) {
    goto cnv_make_map_ret_OPEN_FAIL;
  }
  retval = cnv_first_nonheader_line(cnvfile, &line_idx);
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
  line_idx--;
  do {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "\nError: Line %" PRIuPTR " of .cnv file is pathologically long.\n", line_idx);
      goto cnv_make_map_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (!is_eoln_kns(*bufptr)) {
      // FID, IID, CHR, BP1, BP2, TYPE, SCORE, SITES
      bufptr = next_token_mult(bufptr, 2);
      bufptr2 = next_token_mult(bufptr, req_fields);
      if (no_more_tokens_kns(bufptr2)) {
	sprintf(logbuf, "\nError: Line %" PRIuPTR " of .cnv file has fewer tokens than expected.\n", line_idx);
	goto cnv_make_map_ret_INVALID_FORMAT_2;
      }
      ii = get_chrom_code(chrom_info_ptr, bufptr);
      if (ii < 0) {
	if ((!allow_extra_chroms) || (ii == -1)) {
	  sprintf(logbuf, "\nError: Invalid chromosome code on line %" PRIuPTR " of .cnv file.\n", line_idx);
          goto cnv_make_map_ret_INVALID_FORMAT_2;
	}
        retval = resolve_or_add_chrom_name(chrom_info_ptr, bufptr, &ii, line_idx, ".cnv file");
	if (retval) {
	  goto cnv_make_map_ret_1;
	}
      }
      chrom_idx = ii;
      if ((!is_autogen) && (!IS_SET(chrom_mask, chrom_idx))) {
	continue;
      }
      ullii = ((uint64_t)chrom_idx) << 32;
      bufptr2 = next_token(bufptr);
      bufptr = next_token(bufptr2);
      if (scan_uint_defcap(bufptr2, &seg_start) || scan_uint_defcap(bufptr, &seg_end)) {
	sprintf(logbuf, "\nError: Invalid bp coordinate on line %" PRIuPTR " of .cnv file.\n", line_idx);
        goto cnv_make_map_ret_INVALID_FORMAT_2;
      }
      if (seg_end < seg_start) {
	sprintf(logbuf, "\nError: Segment end coordinate smaller than segment start on line %" PRIuPTR " of\n.cnv file.\n", line_idx);
	goto cnv_make_map_ret_INVALID_FORMAT_2;
      }
      if ((marker_pos_start > (int32_t)seg_start) || ((marker_pos_end != -1) && (marker_pos_end < (int32_t)seg_end))) {
	continue;
      }
      llii = (int64_t)(ullii | ((uint64_t)seg_start));
      ullii |= (uint64_t)seg_end;
      if (filter_seglen) {
	uii = 1 + seg_end - seg_start;
	if ((uii < min_seglen) || (uii > max_seglen)) {
	  continue;
	}
      }
      if (cnv_calc_type & (CNV_DEL | CNV_DUP)) {
	bufptr2 = next_token(bufptr);
	if (scan_uint_defcap(bufptr2, (uint32_t*)&ii)) {
	  sprintf(logbuf, "\nError: Invalid variant copy count on line %" PRIuPTR " of .cnv file.\n", line_idx);
	  goto cnv_make_map_ret_INVALID_FORMAT_2;
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
	bufptr2 = next_token_mult(bufptr, 2);
	if (scan_double(bufptr2, &dxx)) {
          sprintf(logbuf, "\nError: Invalid confidence score on line %" PRIuPTR " of .cnv file.\n", line_idx);
	  goto cnv_make_map_ret_INVALID_FORMAT_2;
	}
	if ((dxx < min_score) || (dxx > max_score)) {
	  continue;
	}
      }
      if (filter_sites) {
	bufptr2 = next_token_mult(bufptr, 3);
	if (scan_posint_defcap(bufptr2, (uint32_t*)&ii)) {
	  sprintf(logbuf, "\nError: Invalid probe count on line %" PRIuPTR " of .cnv file.\n", line_idx);
	  goto cnv_make_map_ret_INVALID_FORMAT_2;
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
  if (!feof(cnvfile)) {
    goto cnv_make_map_ret_READ_FAIL;
  }
  if (!raw_marker_ct) {
    if (cnv_calc_type) {
      logprint("\nError: No variants.\n");
    } else {
      logprint("\nError: No variants after filtering.\n");
    }
    goto cnv_make_map_ret_INVALID_FORMAT;
  }
#ifdef __cplusplus
  std::sort(marker_pos_arr, &(marker_pos_arr[raw_marker_ct]));
#else
  qsort(marker_pos_arr, raw_marker_ct, sizeof(int64_t), llcmp);
#endif
  llii = marker_pos_arr[0];
  chrom_idx = (uint32_t)(((uint64_t)llii) >> 32);
  for (uii = 0; uii <= chrom_idx; uii++) {
    marker_chrom_start[uii] = 0;
  }
  if (cnv_make_map_write(new_mapfile, chrom_info_ptr, chrom_idx, (uint32_t)((uint64_t)llii), max_marker_id_len_ptr)) {
    goto cnv_make_map_ret_WRITE_FAIL;
  }
  for (ulii = 1; ulii < raw_marker_ct; ulii++) {
    if (marker_pos_arr[ulii] != llii) {
#ifdef __LP64__
      if ((++distinct_marker_ct) == 0x80000000U) {
	logprint("\nError: Too many distinct .cnv.map positions (max 2^31 - 1).\n");
	goto cnv_make_map_ret_INVALID_FORMAT;
      }
#endif
      llii = marker_pos_arr[ulii];
      uii = (uint32_t)(((uint64_t)llii) >> 32);
      if (uii > chrom_idx) {
	do {
	  marker_chrom_start[++chrom_idx] = distinct_marker_ct;
	} while (chrom_idx < uii);
      }
      if (cnv_make_map_write(new_mapfile, chrom_info_ptr, chrom_idx, (uint32_t)((uint64_t)llii), max_marker_id_len_ptr)) {
	goto cnv_make_map_ret_WRITE_FAIL;
      }
    }
  }
  chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  do {
    marker_chrom_start[++chrom_idx] = distinct_marker_ct;
  } while (chrom_idx < chrom_code_end);
  if (fclose_null(&new_mapfile)) {
    goto cnv_make_map_ret_WRITE_FAIL;
  }
  logprint("done.\n");
  while (0) {
  cnv_make_map_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  cnv_make_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  cnv_make_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  cnv_make_map_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  cnv_make_map_ret_INVALID_FORMAT_2:
    logprintb();
  cnv_make_map_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 cnv_make_map_ret_1:
  fclose_cond(new_mapfile);
  wkspace_reset(marker_pos_arr);
  return retval;
}

int32_t validate_cnv_map(FILE** mapfile_ptr, char* mapname, int32_t* marker_pos_start_ptr, int32_t* marker_pos_end_ptr, uint32_t allow_extra_chroms, Chrom_info* chrom_info_ptr, uintptr_t* max_marker_id_len_ptr, uint32_t* marker_chrom_start) {
  // also initializes chrom_names[] if necessary
  uint32_t chrom_idx = 0;
  int32_t last_pos = -1;
  uint32_t marker_ct = 0;
  int32_t retval = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t line_idx = 0;
  int32_t marker_pos_start_m1 = (*marker_pos_start_ptr) - 1;
  int32_t marker_pos_end = 0x7fffffff;
  uintptr_t* chrom_mask = chrom_info_ptr->chrom_mask;
  uint32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  char* bufptr;
  char* bufptr2;
  uint32_t colskip;
  uintptr_t cur_marker_id_len;
  int32_t ii;
  if ((*marker_pos_end_ptr) != -1) {
    marker_pos_end = *marker_pos_end_ptr;
  }
  if (fopen_checked(mapfile_ptr, mapname, "r")) {
    goto validate_cnv_map_ret_OPEN_FAIL;
  }
  marker_chrom_start[0] = 0;
  do {
    line_idx++;
    if (!fgets(tbuf, MAXLINELEN, *mapfile_ptr)) {
      if (feof(*mapfile_ptr)) {
	logprint("Error: Empty .cnv.map file.\n");
	goto validate_cnv_map_ret_INVALID_FORMAT;
      } else {
	goto validate_cnv_map_ret_READ_FAIL;
      }
    }
    if (!tbuf[MAXLINELEN - 1]) {
      goto validate_cnv_map_ret_LONG_LINE;
    }
    bufptr = skip_initial_spaces(tbuf);
  } while (is_eoln_kns(*bufptr));
  bufptr2 = next_token_mult(bufptr, 2);
  if (is_eoln_kns(*bufptr2)) {
    goto validate_cnv_map_ret_MISSING_TOKENS;
  }
  bufptr2 = next_token(bufptr2);
  if (is_eoln_kns(*bufptr2)) {
    // --map3 autodetect
    colskip = 1;
  } else {
    colskip = 2;
  }
  line_idx--;
  do {
    line_idx++;
    if (!tbuf[MAXLINELEN - 1]) {
      goto validate_cnv_map_ret_LONG_LINE;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    ii = get_chrom_code(chrom_info_ptr, bufptr);
    if (ii < 0) {
      if ((!allow_extra_chroms) || (ii == -1)) {
	sprintf(logbuf, "Error: Invalid chromosome code on line %" PRIuPTR " of .cnv.map file.\n", line_idx);
	goto validate_cnv_map_ret_INVALID_FORMAT_2;
      }
      retval = resolve_or_add_chrom_name(chrom_info_ptr, bufptr, &ii, line_idx, ".cnv.map file");
      if (retval) {
	goto validate_cnv_map_ret_1;
      }
    }
    if (((uint32_t)ii) < chrom_idx) {
      goto validate_cnv_map_ret_UNSORTED;
    }
    bufptr = next_token(bufptr);
    bufptr2 = next_token_mult(bufptr, colskip);
    if (is_eoln_kns(*bufptr2)) {
      goto validate_cnv_map_ret_MISSING_TOKENS;
    }
    if (((uint32_t)ii) > chrom_idx) {
      do {
	marker_chrom_start[++chrom_idx] = marker_ct;
      } while (chrom_idx < ((uint32_t)ii));
      last_pos = -1;
    }
    if (*bufptr2 == '-') {
      continue;
    }
    if (scan_uint_defcap(bufptr2, (uint32_t*)&ii)) {
      sprintf(logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of .cnv.map file.\n", line_idx);
      goto validate_cnv_map_ret_INVALID_FORMAT_2;
    }
    if (ii <= last_pos) {
      goto validate_cnv_map_ret_UNSORTED;
    }
    last_pos = ii;
    if (!IS_SET(chrom_mask, chrom_idx)) {
      continue;
    }
    if ((last_pos <= marker_pos_start_m1) || (last_pos > marker_pos_end)) {
      continue;
    }
    cur_marker_id_len = strlen_se(bufptr);
    bufptr[cur_marker_id_len] = '\0';
    if (cur_marker_id_len >= max_marker_id_len) {
      max_marker_id_len = cur_marker_id_len + 1;
    }
    if (++marker_ct == 0x80000000U) {
      logprint("Error: Too many entries in .cnv.map file (max 2147483647).\n");
      goto validate_cnv_map_ret_INVALID_FORMAT;
    }
  } while (fgets(tbuf, MAXLINELEN, *mapfile_ptr));
  if (!feof(*mapfile_ptr)) {
    goto validate_cnv_map_ret_READ_FAIL;
  }
  do {
    marker_chrom_start[++chrom_idx] = marker_ct;
  } while (chrom_idx < chrom_code_end);
  *max_marker_id_len_ptr = max_marker_id_len;
  rewind(*mapfile_ptr);
  while (0) {
  validate_cnv_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  validate_cnv_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  validate_cnv_map_ret_UNSORTED:
    logprint("Error: Unsorted .cnv.map file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  validate_cnv_map_ret_MISSING_TOKENS:
    LOGPRINTF("Error: Line %" PRIuPTR " of .cnv.map file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  validate_cnv_map_ret_LONG_LINE:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of .cnv.map is pathologically long.\n", line_idx);
  validate_cnv_map_ret_INVALID_FORMAT_2:
    logprintb();
  validate_cnv_map_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 validate_cnv_map_ret_1:
  return retval;
}

int32_t load_cnv_map(FILE* mapfile, int32_t marker_pos_start, int32_t marker_pos_end, Chrom_info* chrom_info_ptr, uintptr_t max_marker_id_len, uint32_t* marker_pos, char* marker_ids) {
  int32_t retval = 0;
  uint32_t chrom_idx = 0;
  uintptr_t* chrom_mask = chrom_info_ptr->chrom_mask;
  uint32_t colskip;
  char* bufptr;
  char* bufptr2;
  uint32_t cur_marker_id_len;
  int32_t cur_pos;
  // don't need to worry about invalid format
  if (marker_pos_end == -1) {
    marker_pos_end = 0x7fffffff;
  }
  do {
    if (!fgets(tbuf, MAXLINELEN, mapfile)) {
      goto load_cnv_map_ret_READ_FAIL;
    }
    bufptr = skip_initial_spaces(tbuf);
  } while (is_eoln_kns(*bufptr));
  bufptr = next_token_mult(bufptr, 3);
  if (is_eoln_kns(*bufptr)) {
    colskip = 1;
  } else {
    colskip = 2;
  }
  do {
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    chrom_idx = get_chrom_code(chrom_info_ptr, bufptr);
    if (!IS_SET(chrom_mask, chrom_idx)) {
      continue;
    }
    bufptr = next_token(bufptr);
    bufptr2 = next_token_mult(bufptr, colskip);
    if (*bufptr2 == '-') {
      continue;
    }
    scan_uint_defcap(bufptr2, (uint32_t*)&cur_pos);
    if ((cur_pos < marker_pos_start) || (cur_pos > marker_pos_end)) {
      continue;
    }
    cur_marker_id_len = strlen_se(bufptr);
    bufptr[cur_marker_id_len] = '\0';
    memcpy(marker_ids, bufptr, cur_marker_id_len + 1);
    marker_ids = &(marker_ids[max_marker_id_len]);
    *marker_pos++ = (uint32_t)cur_pos;
  } while (fgets(tbuf, MAXLINELEN, mapfile));
  if (!feof(mapfile)) {
    goto load_cnv_map_ret_READ_FAIL;
  }
  while (0) {
  load_cnv_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
  return retval;
}

int32_t plink_cnv(char* outname, char* outname_end, char* cnvname, char* mapname, char* famname, char* phenoname, char* keepname, char* removename, char* filtername, uint64_t misc_flags, Two_col_params* update_chr, Two_col_params* update_cm, Two_col_params* update_map, Two_col_params* update_name, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* filtervals_flattened, uint64_t filter_flags, uint32_t cnv_calc_type, uint32_t min_seglen, uint32_t max_seglen, double min_score, double max_score, uint32_t min_sites, uint32_t max_sites, uint32_t intersect_filter_type, char* intersect_filter_fname, char* subset_fname, uint32_t overlap_type, double overlap_val, uint32_t freq_type, uint32_t freq_val, double freq_val2, uint32_t test_window, uint32_t segment_modifier, char* segment_spanning_fname, uint32_t sample_mperms, uint32_t test_mperms, uint32_t test_region_mperms, uint32_t enrichment_test_mperms, int32_t marker_pos_start, int32_t marker_pos_end, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* cnvfile = NULL;
  FILE* famfile = NULL;
  FILE* mapfile = NULL;
  FILE* outfile = NULL;
  char* subset_list = NULL;
  uintptr_t subset_ct = 0;
  uintptr_t max_subset_name_len = 0;
  uintptr_t topsize = 0;
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint64_t* il_small = NULL; // high-order 32 bits = 2x center pos,
                             // low-order 32 bits = interval end - start
  uint64_t* il_large = NULL;
  uintptr_t* il_chrom_start_small;
  uintptr_t* il_chrom_start_large;
  unsigned char* wkspace_mark2;
  uint32_t* il_chrom_max_width_small;
  uint32_t* il_chrom_max_width_large;
  uint32_t* marker_chrom_start;
  uint32_t* marker_pos;
  char* marker_ids;
  uintptr_t max_marker_id_len;
  int32_t retval;
  char* sptr;
  uintptr_t ulii;
  uint32_t uii;
  if (wkspace_alloc_ul_checked(&il_chrom_start_small, (MAX_POSSIBLE_CHROM + 1) * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&il_chrom_start_large, (MAX_POSSIBLE_CHROM + 1) * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&il_chrom_max_width_small, MAX_POSSIBLE_CHROM * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&il_chrom_max_width_large, MAX_POSSIBLE_CHROM * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&marker_chrom_start, (MAX_POSSIBLE_CHROM + 1) * sizeof(int32_t))) {
    goto plink_cnv_ret_NOMEM;
  }
  wkspace_mark2 = wkspace_base;
  if (fopen_checked(&cnvfile, cnvname, "r")) {
    goto plink_cnv_ret_OPEN_FAIL;
  }
  if (cnv_calc_type & (~CNV_MAKE_MAP)) {
    if (fopen_checked(&famfile, famname, "r")) {
      goto plink_cnv_ret_OPEN_FAIL;
    }
  }
  tbuf[MAXLINELEN - 1] = ' ';
  if (intersect_filter_type) {
    if (subset_fname) {
      retval = cnv_subset_load(subset_fname, &subset_list, &subset_ct, &max_subset_name_len);
      if (retval) {
	goto plink_cnv_ret_1;
      }
    }
    retval = cnv_intersect_load(intersect_filter_type, intersect_filter_fname, subset_list, subset_ct, max_subset_name_len, il_chrom_start_small, il_chrom_start_large, il_chrom_max_width_small, il_chrom_max_width_large, &il_small, &il_large, marker_pos_start, marker_pos_end, allow_extra_chroms, chrom_info_ptr, &topsize);
    if (retval) {
      goto plink_cnv_ret_1;
    }
    wkspace_reset(wkspace_mark2);
  }
  if (!(cnv_calc_type & CNV_MAKE_MAP)) {
    sptr = (char*)memchr(mapname, 0, FNAMESIZE);
    if ((mapname[0] == '\0') || (!filename_exists(mapname, sptr, ""))) {
      if (mapname[0] == '\0') {
        uii = strlen(cnvname);
        if ((uii < 5) || (cnvname[uii - 4] != '.') || (!match_upper_nt(&(cnvname[uii - 3]), "CNV", 3))) {
	  logprint("Error: No .cnv.map filename specified, and .cnv filename does not unambiguously\nspecify how an autogenerated file should be named.  Use --cnv-make-map + --out.\n");
	  goto plink_cnv_ret_INVALID_CMDLINE;
	}
	memcpy(mapname, cnvname, uii);
	memcpy(&(mapname[uii]), ".map", 5);
	sptr = &(mapname[uii + 4]);
	if (filename_exists(mapname, sptr, "")) {
	  LOGPRINTFWW("Error: No .cnv.map filename specified, and natural autogeneration target\n(%s) already exists.\n", mapname);
	  goto plink_cnv_ret_INVALID_CMDLINE;
	}
      }
      sprintf(logbuf, "Autogenerating missing %s ... ", mapname);
      wordwrap(logbuf, 5);
      retval = cnv_make_map(cnvfile, mapname, 0, 0, 0xffffffffU, -HUGE_DOUBLE, HUGE_DOUBLE, 0, 0xffffffffU, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0.0, -1, -1, allow_extra_chroms, 0, chrom_info_ptr, &max_marker_id_len, marker_chrom_start);
    } else {
      retval = validate_cnv_map(&mapfile, mapname, &marker_pos_start, &marker_pos_end, allow_extra_chroms, chrom_info_ptr, &max_marker_id_len, marker_chrom_start);
    }
    if (retval) {
      goto plink_cnv_ret_1;
    }
  } else {
    memcpy(outname_end, ".cnv.map", 9);
    sprintf(logbuf, "Generating %s ... ", outname);
    wordwrap(logbuf, 5);
    retval = cnv_make_map(cnvfile, outname, cnv_calc_type, min_seglen, max_seglen, min_score, max_score, min_sites, max_sites, il_chrom_start_small, il_chrom_start_large, il_chrom_max_width_small, il_chrom_max_width_large, il_small, il_large, intersect_filter_type, overlap_type, overlap_val, marker_pos_start, marker_pos_end, allow_extra_chroms, 0, chrom_info_ptr, &max_marker_id_len, marker_chrom_start);
    if (retval || (!(cnv_calc_type & (CNV_MAKE_MAP | CNV_DEL | CNV_DUP)))) {
      goto plink_cnv_ret_1;
    }
  }
  if (!mapfile) {
    if (fopen_checked(&mapfile, mapname, "r")) {
      goto plink_cnv_ret_OPEN_FAIL;
    }
  }
  ulii = marker_chrom_start[chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct];
#ifndef __LP64__
  if (((uint64_t)ulii) * max_marker_id_len > 0x7fffffffLLU) {
    goto plink_cnv_ret_NOMEM;
  }
#endif
  if (wkspace_alloc_ui_checked(&marker_pos, ulii * sizeof(int32_t)) ||
      wkspace_alloc_c_checked(&marker_ids, ulii * max_marker_id_len)) {
    goto plink_cnv_ret_NOMEM;
  }
  retval = load_cnv_map(mapfile, marker_pos_start, marker_pos_end, chrom_info_ptr, max_marker_id_len, marker_pos, marker_ids);
  if (retval) {
    goto plink_cnv_ret_1;
  }
  printf(".cnv.map file successfully loaded.  (Subsequent functions have not been\nimplemented yet.)\n");
  while (0) {
  plink_cnv_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  plink_cnv_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  plink_cnv_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 plink_cnv_ret_1:
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
#endif // HIGH_MAX_CHROM

int32_t plink_gvar(char* outname, char* outname_end, char* gvarname, char* mapname, char* famname) {
  logprint("Error: Common CNP analysis not yet supported.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
}
