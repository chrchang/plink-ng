#include "plink_set.h"

void set_init(Set_info* sip) {
  sip->fname = NULL;
  sip->subset_fname = NULL;
  sip->merged_set_name = NULL;
  sip->genekeep_flattened = NULL;
  sip->ct = 0;
  sip->modifier = 0;
}

void set_cleanup(Set_info* sip) {
  free_cond(sip->fname);
  free_cond(sip->subset_fname);
  free_cond(sip->merged_set_name);
  free_cond(sip->genekeep_flattened);
}

typedef struct make_set_range_struct {
  struct make_set_range_struct* next;
  uint32_t uidx_start;
  uint32_t uidx_end;
} Make_set_range;

int32_t define_sets(Set_info* sip, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, uintptr_t* marker_exclude_ct_ptr, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr) {
  FILE* infile = NULL;
  uintptr_t topsize = 0;
  char* sorted_marker_ids = NULL;
  char* sorted_genekeep_ids = NULL;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t marker_exclude_ct = *marker_exclude_ct_ptr;
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t set_ct = 0;
  uint32_t make_set = sip->modifier & SET_MAKE_FROM_RANGES;
  uint32_t complement_sets = (sip->modifier / SET_COMPLEMENTS) & 1;
  uint32_t collapse_group = sip->modifier & SET_MAKE_COLLAPSE_GROUP;
  uint32_t c_prefix = 2 * ((sip->modifier / SET_C_PREFIX) & 1);
  uint32_t gene_all = sip->modifier & SET_GENE_ALL;
  uint32_t make_set_border = sip->make_set_border;
  uint32_t curtoklen = 0;
  uint32_t in_set = 0;
  int32_t retval = 0;
  uintptr_t subset_ct = 0;
  uintptr_t max_subset_id_len = 0;
  uintptr_t genekeep_ct = 0;
  uintptr_t max_genekeep_len = 0;
  uintptr_t max_set_id_len = 0;
  Ll_str* make_set_ll = NULL;
  Make_set_range** make_set_range_arr = NULL;
  char* midbuf = &(tbuf[MAXLINELEN]);
  char* sorted_subset_ids = NULL;
  char* set_names = NULL;
  char* bufptr = NULL;
  uint64_t* range_sort_buf = NULL;
  char* bufptr2;
  char* bufptr3;
  char* buf_end;
  Ll_str* ll_tmp;
  Make_set_range* msr_tmp;
  uint32_t* marker_id_map;
  uint32_t* marker_uidx_to_idx;
  uint32_t** set_range_ptrs;
  uint32_t* set_bounds;
  uintptr_t* marker_exclude_new;
  uintptr_t* marker_bitfield_tmp;
  uintptr_t** set_bitfield_ptrs;
  uintptr_t* include_out_of_bounds;
  uint32_t* uiptr;
  uintptr_t set_ctl;
  uintptr_t set_idx;
  uintptr_t bufsize;
  uintptr_t topsize_bak;
  uintptr_t ulii;
  uint64_t ullii;
  uint32_t chrom_idx;
  uint32_t chrom_start;
  uint32_t chrom_end;
  uint32_t range_first;
  uint32_t range_last;
  uint32_t slen;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  int32_t ii;
  // 1. validate and sort --gene parameter(s)
  if (sip->genekeep_flattened) {
    bufptr = sip->genekeep_flattened;
    if (sip->merged_set_name) {
      // degenerate case: --gene-all if merged set name present, fail
      // otherwise
      uii = strlen(sip->merged_set_name);
      while (1) {
	slen = strlen(bufptr);
	if ((slen == uii) && (!memcmp(bufptr, sip->merged_set_name, uii))) {
	  break;
	}
	bufptr = &(bufptr[slen + 1]);
	if (!(*bufptr)) {
	  goto define_sets_ret_ALL_MARKERS_EXCLUDED;
	}
      }
      free(sip->genekeep_flattened);
      sip->genekeep_flattened = NULL;
      gene_all = 1;
    } else {
      do {
	slen = strlen(bufptr) + 1;
	if ((!c_prefix) || (!memcmp(bufptr, "C_", 2))) {
	  if (slen > max_genekeep_len) {
	    max_genekeep_len = slen;
	  }
	  genekeep_ct++;
	}
	bufptr = &(bufptr[slen]);
      } while (*bufptr);
      if (!genekeep_ct) {
	logprint("Error: All variants excluded by --gene.\n");
	goto define_sets_ret_ALL_MARKERS_EXCLUDED;
      }
      sorted_genekeep_ids = (char*)top_alloc(&topsize, genekeep_ct * max_genekeep_len);
      if (!sorted_genekeep_ids) {
	goto define_sets_ret_NOMEM;
      }
      bufptr = sip->genekeep_flattened;
      ulii = 0;
      do {
	slen = strlen(bufptr) + 1;
	if ((!c_prefix) || (!memcmp(bufptr, "C_", 2))) {
	  memcpy(&(sorted_genekeep_ids[ulii * max_genekeep_len]), bufptr, slen);
	  ulii++;
	}
	bufptr = &(bufptr[slen]);
      } while (*bufptr);
      qsort(sorted_genekeep_ids, genekeep_ct, max_genekeep_len, strcmp_casted);
    }
  }
  // 2. if --subset is present, load and sort that file
  if (sip->subset_fname) {
    if (fopen_checked(&infile, sip->subset_fname, "rb")) {
      goto define_sets_ret_OPEN_FAIL;
    }
    retval = scan_token_ct_len(infile, tbuf, MAXLINELEN, &subset_ct, &max_subset_id_len);
    if (retval) {
      if (retval == RET_INVALID_FORMAT) {
        logprint("Error: Pathologically long token in --subset file.\n");
      }
      goto define_sets_ret_1;
    }
    if (!subset_ct) {
      if ((gene_all || sip->genekeep_flattened) && ((!sip->merged_set_name) || (!complement_sets))) {
	logprint("Error: All variants excluded, since --subset file is empty.\n");
	goto define_sets_ret_ALL_MARKERS_EXCLUDED;
      }
      if (sip->merged_set_name) {
	goto define_sets_merge_nothing;
      } else {
        logprint("Warning: Empty --subset file; no sets defined.\n");
        goto define_sets_ret_1;
      }
    }
    sorted_subset_ids = (char*)top_alloc(&topsize, subset_ct * max_subset_id_len);
    if (!sorted_subset_ids) {
      goto define_sets_ret_NOMEM;
    }
    rewind(infile);
    retval = read_tokens(infile, tbuf, MAXLINELEN, subset_ct, max_subset_id_len, sorted_subset_ids);
    if (retval) {
      goto define_sets_ret_1;
    }
    qsort(sorted_subset_ids, subset_ct, max_subset_id_len, strcmp_casted);
    subset_ct = collapse_duplicate_ids(sorted_subset_ids, subset_ct, max_subset_id_len, NULL);
    if (fclose_null(&infile)) {
      goto define_sets_ret_READ_FAIL;
    }
  }
  if (fopen_checked(&infile, sip->fname, "r")) {
    goto define_sets_ret_OPEN_FAIL;
  }
  // 3. load --make-set range list
  if (make_set) {
    tbuf[MAXLINELEN - 1] = ' ';
    // if we need to track set names, put together a sorted list
    if (!sip->merged_set_name) {
      while (fgets(tbuf, MAXLINELEN, infile)) {
	if (!tbuf[MAXLINELEN - 1]) {
	  logprint("Error: Pathologically long line in --make-set file.\n");
	  goto define_sets_ret_INVALID_FORMAT;
	}
	bufptr = skip_initial_spaces(tbuf);
	if (is_eoln_kns(*bufptr)) {
	  continue;
	}
	bufptr2 = next_item_mult(bufptr, 3);
	if (!collapse_group) {
	  bufptr3 = bufptr2;
	} else {
	  bufptr3 = next_item(bufptr2);
	}
	if (no_more_items_kns(bufptr3)) {
	  logprint("Error: Fewer tokens than expected in --make-set file line.\n");
	  goto define_sets_ret_INVALID_FORMAT;
	}
	ii = get_chrom_code(chrom_info_ptr, bufptr);
	if (ii == -1) {
	  logprint("Error: Invalid chromosome code in --make-set file.\n");
	  goto define_sets_ret_INVALID_FORMAT;
	}
	uii = strlen_se(bufptr2);
	bufptr2[uii] = '\0';
	if (subset_ct) {
          ii = bsearch_str(bufptr2, sorted_subset_ids, max_subset_id_len, 0, subset_ct - 1);
	  if (ii == -1) {
	    continue;
	  }
	}
        if (collapse_group) {
	  uii = strlen_se(bufptr3);
	  bufptr3[uii] = '\0';
	}
	// when there are repeats, they are likely to be next to each other
	if (make_set_ll && (!strcmp(make_set_ll->ss, bufptr3))) {
	  continue;
	}
	uii++;
	if (uii > max_set_id_len) {
	  max_set_id_len = uii;
	}
        ll_tmp = top_alloc_llstr(&topsize, uii);
        ll_tmp->next = make_set_ll;
	memcpy(ll_tmp->ss, bufptr3, uii);
	make_set_ll = ll_tmp;
	set_ct++;
      }
      if (!set_ct) {
        if (gene_all || sip->genekeep_flattened) {
	  logprint("Error: All variants excluded by --gene{-all}, since no sets were defined from\n--make-set file.\n");
	  goto define_sets_ret_ALL_MARKERS_EXCLUDED;
	}
	logprint("Warning: No sets defined from --make-set file.\n");
	goto define_sets_ret_1;
      }
      max_set_id_len += c_prefix;
      wkspace_left -= topsize;
      if (wkspace_alloc_c_checked(&set_names, set_ct)) {
	goto define_sets_ret_NOMEM2;
      }
      wkspace_left += topsize;
      if (!c_prefix) {
	for (ulii = 0; ulii < set_ct; ulii++) {
	  strcpy(&(set_names[ulii * max_set_id_len]), make_set_ll->ss);
	  make_set_ll = make_set_ll->next;
	}
      } else {
	for (ulii = 0; ulii < set_ct; ulii++) {
	  memcpy(&(set_names[ulii * max_set_id_len]), "C_", 2);
	  strcpy(&(set_names[ulii * max_set_id_len + 2]), make_set_ll->ss);
	  make_set_ll = make_set_ll->next;
	}
      }
      qsort(set_names, set_ct, max_set_id_len, strcmp_natural);
      set_ct = collapse_duplicate_ids(set_names, set_ct, max_set_id_len, NULL);
      wkspace_reset(set_names);
      set_names = (char*)wkspace_alloc(set_ct * max_set_id_len);
      rewind(infile);
    } else {
      set_ct = 1;
    }
    make_set_range_arr = (Make_set_range**)top_alloc(&topsize, set_ct * sizeof(intptr_t));
    for (set_idx = 0; set_idx < set_ct; set_idx++) {
      make_set_range_arr[set_idx] = NULL;
    }
    while (fgets(tbuf, MAXLINELEN, infile)) {
      if (!tbuf[MAXLINELEN - 1]) {
	logprint("Error: Pathologically long line in --make-set file.\n");
	goto define_sets_ret_INVALID_FORMAT;
      }
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      bufptr2 = next_item_mult(bufptr, 3);
      if (!collapse_group) {
	bufptr3 = bufptr2;
      } else {
	bufptr3 = next_item(bufptr2);
      }
      if (no_more_items_kns(bufptr3)) {
	logprint("Error: Fewer tokens than expected in --make-set file line.\n");
	goto define_sets_ret_INVALID_FORMAT;
      }
      ii = get_chrom_code(chrom_info_ptr, bufptr);
      if (ii == -1) {
	logprint("Error: Invalid chromosome code in --make-set file.\n");
	goto define_sets_ret_INVALID_FORMAT;
      }
      if (!is_set(chrom_info_ptr->chrom_mask, ii)) {
	continue;
      }
      chrom_idx = ii;
      chrom_start = chrom_info_ptr->chrom_start[chrom_idx];
      chrom_end = chrom_info_ptr->chrom_end[chrom_idx];
      if (chrom_end == chrom_start) {
	continue;
      }
      uii = strlen_se(bufptr2);
      bufptr2[uii] = '\0';
      if (subset_ct) {
	ii = bsearch_str(bufptr2, sorted_subset_ids, max_subset_id_len, 0, subset_ct - 1);
	if (ii == -1) {
	  continue;
	}
      }
      bufptr = next_item(bufptr);
      if (atoiz2(bufptr, &ii)) {
	bufptr[strlen_se(bufptr)] = '\0';
	sprintf(logbuf, "Error: Invalid range start position '%s' in --make-set file.\n", bufptr);
	goto define_sets_ret_INVALID_FORMAT_2;
      }
      range_first = ii;
      bufptr = next_item(bufptr);
      if (atoiz2(bufptr, &ii)) {
	bufptr[strlen_se(bufptr)] = '\0';
	sprintf(logbuf, "Error: Invalid range end position '%s' in --make-set file.\n", bufptr);
	goto define_sets_ret_INVALID_FORMAT_2;
      }
      range_last = ii;
      if (range_last < range_first) {
	logprint("Error: Range end position smaller than range start in --make-set file.\n");
	goto define_sets_ret_INVALID_FORMAT;
      }
      if (make_set_border > range_first) {
	range_first = 0;
      } else {
	range_first -= make_set_border;
      }
      range_last += make_set_border;
      if (set_ct > 1) {
	if (collapse_group) {
	  uii = strlen_se(bufptr3);
	  bufptr3[uii] = '\0';
	}
	if (c_prefix) {
	  bufptr3 = &(bufptr3[-2]);
          memcpy(bufptr3, "C_", 2);
	}
	// this should never fail
        set_idx = (uint32_t)bsearch_str_natural(bufptr3, set_names, max_set_id_len, 0, set_ct - 1);
      } else {
	set_idx = 0;
      }
      // translate to within-chromosome uidx
      range_first = uint32arr_greater_than(&(marker_pos[chrom_start]), chrom_end - chrom_start, range_first);
      range_last = uint32arr_greater_than(&(marker_pos[chrom_start]), chrom_end - chrom_start, range_last + 1);
      if (range_last > range_first) {
	msr_tmp = (Make_set_range*)top_alloc(&topsize, sizeof(Make_set_range));
	msr_tmp->next = make_set_range_arr[set_idx];
        // normally, I'd keep chrom_idx here since that enables by-chromosome
	// sorting, but that's probably not worth bloating Make_set_range from
	// 16 to 32 bytes
	msr_tmp->uidx_start = chrom_start + range_first;
	msr_tmp->uidx_end = chrom_start + range_last;
	make_set_range_arr[set_idx] = msr_tmp;
      }
    }
    // allocate buffer for sorting ranges later
    uii = 0;
    for (set_idx = 0; set_idx < set_ct; set_idx++) {
      ujj = 0;
      msr_tmp = make_set_range_arr[set_idx];
      while (msr_tmp) {
	ujj++;
        msr_tmp = msr_tmp->next;
      }
      if (ujj > uii) {
	uii = ujj;
      }
    }
    range_sort_buf = (uint64_t*)top_alloc(&topsize, uii * sizeof(int64_t));
  }
  // 4. if --gene or --gene-all is present, pre-filter variants.
  if (gene_all || sip->genekeep_flattened) {
    topsize_bak = topsize;
    marker_bitfield_tmp = (uintptr_t*)top_alloc(&topsize, unfiltered_marker_ctl * sizeof(intptr_t));
    if (!marker_bitfield_tmp) {
      goto define_sets_ret_NOMEM;
    }
    marker_exclude_new = (uintptr_t*)top_alloc(&topsize, unfiltered_marker_ctl * sizeof(intptr_t));
    if (!marker_exclude_new) {
      goto define_sets_ret_NOMEM;
    }
    fill_ulong_zero(marker_bitfield_tmp, unfiltered_marker_ctl);
    fill_all_bits(marker_exclude_new, unfiltered_marker_ct);
    // then include every variant that appears, or include every variant that
    // fails to appear in a fully loaded set in the complement case
    if (make_set) {
      for (set_idx = 0; set_idx < set_ct; set_idx++) {
	if (gene_all || (bsearch_str(&(set_names[set_idx * max_set_id_len]), sorted_genekeep_ids, max_genekeep_len, 0, genekeep_ct - 1) != -1)) {
	  msr_tmp = make_set_range_arr[set_idx];
	  while (msr_tmp) {
	    fill_bits(marker_bitfield_tmp, msr_tmp->uidx_start, msr_tmp->uidx_end - msr_tmp->uidx_start);
	    msr_tmp = msr_tmp->next;
	  }
	}
        if (complement_sets) {
	  bitfield_and(marker_exclude_new, marker_bitfield_tmp, unfiltered_marker_ctl);
          fill_ulong_zero(marker_bitfield_tmp, unfiltered_marker_ctl);
	}
      }
    } else {
      sorted_marker_ids = (char*)top_alloc(&topsize, marker_ct * max_marker_id_len);
      if (!sorted_marker_ids) {
	goto define_sets_ret_NOMEM;
      }
      marker_id_map = (uint32_t*)top_alloc(&topsize, marker_ct * sizeof(int32_t));
      if (!marker_id_map) {
	goto define_sets_ret_NOMEM;
      }
      wkspace_left -= topsize;
      retval = sort_item_ids_noalloc(sorted_marker_ids, marker_id_map, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
      wkspace_left += topsize;
      if (retval) {
	goto define_sets_ret_NOMEM;
      }

      // similar to read_tokens(), since it may be important to support very
      // long lines.
      while (1) {
	if (fread_checked(midbuf, MAXLINELEN, infile, &bufsize)) {
          goto define_sets_ret_READ_FAIL;
	}
        if (!bufsize) {
	  if (curtoklen) {
	    if ((curtoklen != 3) || memcmp(bufptr, "END", 3)) {
	      goto define_sets_ret_INVALID_FORMAT_NO_END;
	    } else if (!in_set) {
	      goto define_sets_ret_INVALID_FORMAT_EXTRA_END;
	    }
	  } else if (in_set) {
	    goto define_sets_ret_INVALID_FORMAT_NO_END;
	  }
	  break;
	}
        buf_end = &(midbuf[bufsize]);
        *buf_end = ' ';
        buf_end[1] = '0';
        bufptr = &(tbuf[MAXLINELEN - curtoklen]);
        bufptr2 = midbuf;
        if (curtoklen) {
          goto define_sets_tok_start_1;
	}
        while (1) {
          while (*bufptr <= ' ') {
	    bufptr++;
	  }
          if (bufptr >= buf_end) {
	    curtoklen = 0;
	    break;
	  }
	  bufptr2 = &(bufptr[1]);
	define_sets_tok_start_1:
          while (*bufptr2 > ' ') {
	    bufptr2++;
	  }
          curtoklen = (uintptr_t)(bufptr2 - bufptr);
          if ((bufptr2 == buf_end) && (buf_end == &(tbuf[MAXLINELEN * 2]))) {
            bufptr3 = &(tbuf[MAXLINELEN - curtoklen]);
            memcpy(bufptr3, bufptr, curtoklen);
            bufptr = bufptr3;
	    break;
	  }
          if ((curtoklen == 3) && (!memcmp(bufptr, "END", 3))) {
	    if (!in_set) {
	      goto define_sets_ret_INVALID_FORMAT_EXTRA_END;
	    }
            if (complement_sets) {
	      bitfield_and(marker_exclude_new, marker_bitfield_tmp, unfiltered_marker_ctl);
              fill_ulong_zero(marker_bitfield_tmp, unfiltered_marker_ctl);
	    }
            in_set = 0;
	  } else if (!in_set) {
	    if (subset_ct) {
	      *bufptr2 = '\0';
	      ii = bsearch_str(bufptr, sorted_subset_ids, max_subset_id_len, 0, subset_ct - 1);
	      if (ii == -1) {
		in_set = 2; // ignore this set
		bufptr = &(bufptr2[1]);
		continue;
	      }
	    }
	    if (curtoklen >= max_set_id_len) {
	      max_set_id_len = curtoklen + 1;
	    }
	    set_ct++;
	    in_set = 1;
	  } else if (in_set == 1) {
	    *bufptr2 = '\0';
	    ii = bsearch_str(bufptr, sorted_marker_ids, max_marker_id_len, 0, marker_ct - 1);
	    if (ii != -1) {
	      set_bit(marker_bitfield_tmp, marker_id_map[(uint32_t)ii]);
	    }
	  }
	  bufptr = &(bufptr2[1]);
	}
      }
      if (!feof(infile)) {
	goto define_sets_ret_READ_FAIL;
      }
      if (!set_ct) {
	if (!complement_sets) {
	  logprint("Error: All variants excluded by --gene{-all}, since no sets were defined from\n--set file.\n");
	  goto define_sets_ret_ALL_MARKERS_EXCLUDED;
	}
	logprint("Warning: No sets defined from --set file.\n");
	goto define_sets_ret_1;
      }
    }
    if (!complement_sets) {
      bitfield_andnot(marker_exclude_new, marker_bitfield_tmp, unfiltered_marker_ctl);
    }
    bitfield_or(marker_exclude, marker_exclude_new, unfiltered_marker_ctl);
    marker_exclude_ct = popcount_longs(marker_exclude, 0, unfiltered_marker_ctl);
    if (marker_exclude_ct == unfiltered_marker_ct) {
      goto define_sets_ret_ALL_MARKERS_EXCLUDED;
    }
    *marker_exclude_ct_ptr = marker_exclude_ct;
    marker_ct = unfiltered_marker_ct - marker_exclude_ct;
    rewind(infile);
    topsize = topsize_bak;
  } else if ((!make_set) && (!sip->merged_set_name)) {
    // 5. otherwise, with --set and no --set-collapse-all, count number of sets
    //    and max_name_len.
    while (1) {
      if (fread_checked(midbuf, MAXLINELEN, infile, &bufsize)) {
	goto define_sets_ret_READ_FAIL;
      }
      if (!bufsize) {
	if (curtoklen) {
	  if ((curtoklen != 3) || (memcmp(bufptr, "END", 3))) {
	    goto define_sets_ret_INVALID_FORMAT_NO_END;
	  } else if (!in_set) {
	    goto define_sets_ret_INVALID_FORMAT_EXTRA_END;
	  }
	} else if (in_set) {
          goto define_sets_ret_INVALID_FORMAT_NO_END;
	}
	break;
      }
      buf_end = &(midbuf[bufsize]);
      *buf_end = ' ';
      buf_end[1] = '0';
      bufptr = &(tbuf[MAXLINELEN - curtoklen]);
      bufptr2 = midbuf;
      if (curtoklen) {
	goto define_sets_tok_start_2;
      }
      while (1) {
	while (*bufptr <= ' ') {
	  bufptr++;
	}
        if (bufptr >= buf_end) {
          curtoklen = 0;
	  break;
	}
        bufptr2 = &(bufptr[1]);
      define_sets_tok_start_2:
        while (*bufptr2 > ' ') {
	  bufptr2++;
	}
        curtoklen = (uintptr_t)(bufptr2 - bufptr);
        if ((bufptr2 == buf_end) && (buf_end == &(tbuf[MAXLINELEN * 2]))) {
          bufptr3 = &(tbuf[MAXLINELEN - curtoklen]);
          memcpy(bufptr3, bufptr, curtoklen);
	  bufptr = bufptr3;
	  break;
	}
        if ((curtoklen == 3) && (!memcmp(bufptr, "END", 3))) {
          if (!in_set) {
	    goto define_sets_ret_INVALID_FORMAT_EXTRA_END;
	  }
	  in_set = 0;
	} else if (!in_set) {
	  in_set = 1;
	  if (subset_ct) {
	    *bufptr2 = '\0';
            ii = bsearch_str(bufptr, sorted_subset_ids, max_subset_id_len, 0, subset_ct - 1);
            if (ii == -1) {
	      // no need for in_set = 2, just don't adjust set_ct/id_len
              bufptr = &(bufptr2[1]);
	      continue;
	    }
	  }
	  if (curtoklen >= max_set_id_len) {
	    max_set_id_len = curtoklen + 1;
	  }
	  set_ct++;
	}
	bufptr = &(bufptr2[1]);
      }
    }
    if (!set_ct) {
      logprint("Warning: No sets defined from --set file.\n");
      goto define_sets_ret_1;
    }
    rewind(infile);
  }
  // 6. allocate sip->names[], range_ptrs[], bounds[], bitfield_ptrs[],
  //    include_out_of_bounds[] on stack
  marker_uidx_to_idx = (uint32_t*)top_alloc(&topsize, unfiltered_marker_ct * sizeof(int32_t));
  if (!marker_uidx_to_idx) {
    goto define_sets_ret_NOMEM;
  }
  fill_uidx_to_idx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_uidx_to_idx);
  wkspace_left -= topsize;
  if (!set_names) {
    if (sip->merged_set_name) {
      set_ct = 1;
      max_set_id_len = strlen(sip->merged_set_name) + 1;
      if (wkspace_alloc_c_checked(&set_names, max_set_id_len)) {
	goto define_sets_ret_NOMEM2;
      }
      memcpy(set_names, sip->merged_set_name, max_set_id_len);
    } else {
      if (wkspace_alloc_c_checked(&set_names, set_ct * max_set_id_len)) {
	goto define_sets_ret_NOMEM2;
      }
    }
  }
  set_range_ptrs = (uint32_t**)wkspace_alloc(set_ct * sizeof(intptr_t));
  if (!set_range_ptrs) {
    goto define_sets_ret_NOMEM2;
  }
  if (wkspace_alloc_ui_checked(&set_bounds, 2 * set_ct * sizeof(int32_t))) {
    goto define_sets_ret_NOMEM2;
  }
  set_bitfield_ptrs = (uintptr_t**)wkspace_alloc(set_ct * sizeof(intptr_t));
  if (!set_bitfield_ptrs) {
    goto define_sets_ret_NOMEM2;
  }
  set_ctl = (set_ct + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(&include_out_of_bounds, set_ctl * sizeof(intptr_t))) {
    goto define_sets_ret_NOMEM2;
  }
  fill_ulong_zero(include_out_of_bounds, set_ctl);
  if (make_set) {
    // 7. If --make-set, allocate entries on stack
    for (set_idx = 0; set_idx < set_ct; set_idx++) {
      msr_tmp = make_set_range_arr[set_idx];
      // sort and merge intervals in O(n log n) instead of O(n^2) time
      uii = 0;
      while (msr_tmp) {
	range_first = msr_tmp->uidx_start;
	range_last = msr_tmp->uidx_end;
	msr_tmp = msr_tmp->next;
	if (IS_SET(marker_exclude, range_first)) {
	  range_first = next_unset(marker_exclude, range_first, unfiltered_marker_ct);
	  if (range_first == unfiltered_marker_ct) {
	    continue;
	  }
	}
	range_first = marker_uidx_to_idx[range_first];
	if (IS_SET(marker_exclude, range_last)) {
          range_last = next_unset(marker_exclude, range_last, unfiltered_marker_ct);
	}
	if (range_last == unfiltered_marker_ct) {
	  range_last = marker_ct;
	} else {
          range_last = marker_uidx_to_idx[range_last];
	  if (range_last == range_first) {
	    continue;
	  }
	}
	range_sort_buf[uii++] = (((uint64_t)range_first) << 32) | ((uint64_t)range_last);
      }
      if (!uii) {
	// special case: empty set
	if (wkspace_left < 16) {
	  goto define_sets_ret_NOMEM;
	}
	set_range_ptrs[set_idx] = (uint32_t*)wkspace_base;
	set_bitfield_ptrs[set_idx] = NULL;
	wkspace_left -= 16;
	wkspace_base = &(wkspace_base[16]);
	set_range_ptrs[set_idx][0] = 0;
	continue;
      }
#ifdef __cplusplus
      std::sort((int64_t*)range_sort_buf, (int64_t*)(&(range_sort_buf[uii])));
#else
      qsort((int64_t*), uii, sizeof(int64_t), llcmp);
#endif
      ukk = 0; // current end of sorted interval list
      range_last = (uint32_t)range_sort_buf[0];
      for (ujj = 1; ujj < uii; ujj++) {
	ullii = range_sort_buf[ujj];
	range_first = (uint32_t)(ullii >> 32);
	if (range_first <= range_last) {
	  umm = (uint32_t)ullii;
	  if (umm > range_last) {
	    range_last = umm;
	    range_sort_buf[ukk] = (range_sort_buf[ukk] & 0xffffffff00000000LLU) | (ullii & 0xffffffffLLU);
	  }
	} else {
	  range_last = (uint32_t)range_sort_buf[++ukk];
	}
      }
      // todo: check if bitfield representation is more compact
      ulii = ((ukk + 1) / 2) + 1;
      ulii *= 16;
      if (wkspace_left < ulii) {
	goto define_sets_ret_NOMEM;
      }
      uiptr = (uint32_t*)wkspace_base;
      set_range_ptrs[set_idx] = uiptr;
      set_bitfield_ptrs[set_idx] = NULL;
      wkspace_left -= ulii;
      wkspace_base = &(wkspace_base[ulii]);
      *uiptr++ = ukk + 1;
      for (ujj = 0; ujj <= ukk; ujj++) {
	ullii = range_sort_buf[ujj];
        *uiptr++ = (uint32_t)(ullii >> 32);
	*uiptr++ = (uint32_t)ullii;
      }
    }
  } else {
    // 8. If --set, load sets and allocate on stack (todo: check if range
    //    representation is more compact)
    logprint("Error: --set is currently under development.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto define_sets_ret_1;
  }
  wkspace_left += topsize;
  sip->ct = set_ct;
  sip->names = set_names;
  sip->max_name_len = max_set_id_len;
  sip->range_ptrs = set_range_ptrs;
  sip->bounds = set_bounds;
  sip->bitfield_ptrs = set_bitfield_ptrs;
  sip->include_out_of_bounds = include_out_of_bounds;
  while (0) {
  define_sets_merge_nothing:
    sip->ct = 1;
    uii = strlen(sip->merged_set_name) + 1;
    // topsize = 0;
    sip->range_ptrs = (uint32_t**)wkspace_alloc(sizeof(intptr_t));
    if (!sip->range_ptrs) {
      goto define_sets_ret_NOMEM;
    }
    if (wkspace_alloc_c_checked(&sip->names, uii) ||
	wkspace_alloc_ui_checked(&sip->bounds, 2 * sizeof(int32_t)) ||
	wkspace_alloc_ui_checked(&(sip->range_ptrs[0]), (1 + 2 * complement_sets) * sizeof(int32_t)) ||
	wkspace_alloc_ul_checked(&sip->include_out_of_bounds, sizeof(intptr_t))) {
      goto define_sets_ret_NOMEM;
    }
    sip->bitfield_ptrs = (uintptr_t**)wkspace_alloc(sizeof(intptr_t));
    if (!sip->bitfield_ptrs) {
      goto define_sets_ret_NOMEM;
    }
    memcpy(sip->names, sip->merged_set_name, uii);
    sip->max_name_len = uii;
    if (complement_sets) {
      sip->range_ptrs[0][0] = 1;
      sip->range_ptrs[0][1] = 0;
      sip->range_ptrs[0][2] = marker_ct;
    } else {
      sip->range_ptrs[0][0] = 0;
    }
    sip->bitfield_ptrs[0] = NULL;
    sip->include_out_of_bounds[0] = 0;
    break;
  define_sets_ret_NOMEM2:
    wkspace_left += topsize;
  define_sets_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  define_sets_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  define_sets_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  define_sets_ret_ALL_MARKERS_EXCLUDED:
    logprint("Error: All variants excluded by --gene/--gene-all.\n");
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  define_sets_ret_INVALID_FORMAT_EXTRA_END:
    logprint("Error: Extra 'END' token in --set file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  define_sets_ret_INVALID_FORMAT_NO_END:
    logprint("Error: Last token in --set file isn't 'END'.\n");
    retval = RET_INVALID_FORMAT;
    break;
  define_sets_ret_INVALID_FORMAT_2:
    logprintb();
  define_sets_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 define_sets_ret_1:
  return retval;
}

int32_t write_set(Set_info* sip, char* outname, char* outname_end, uint32_t marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t set_ct = sip->ct;
  uintptr_t max_set_name_len = sip->max_name_len;
  uintptr_t set_idx = 0;
  uint32_t chrom_idx = 0;
  int32_t retval = 0;
  uint32_t* last_idx;
  uint32_t* next_adj;
  char* cur_setting;
  char* writebuf;
  char* bufptr;
  char* cptr;
  uintptr_t marker_uidx;
  uint32_t marker_idx;
  uint32_t chrom_end;
  uint32_t range_ct;
  uint32_t range_start;
  uint32_t uii;
  if (sip->modifier & SET_WRITE_TABLE) {
    memcpy(outname_end, ".set.table", 11);
    if (fopen_checked(&outfile, outname, "w")) {
      goto write_set_ret_OPEN_FAIL;
    }
    fputs("SNP\tCHR\tBP", outfile);
    for (set_idx = 0; set_idx < set_ct; set_idx++) {
      putc('\t', outfile);
      fputs(&(sip->names[set_idx * max_set_name_len]), outfile);
    }
    if (putc_checked('\n', outfile)) {
      goto write_set_ret_WRITE_FAIL;
    }
    if (wkspace_alloc_ui_checked(&last_idx, set_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&next_adj, set_ct * sizeof(int32_t)) ||
        wkspace_alloc_c_checked(&cur_setting, set_ct) ||
        wkspace_alloc_c_checked(&writebuf, 2 * set_ct)) {
      goto write_set_ret_NOMEM;
    }
    fill_uint_zero(last_idx, set_ct);
    fill_uint_zero(next_adj, set_ct);
    marker_uidx = 0;
    tbuf[0] = '\t';
    chrom_end = 0;
    for (set_idx = 1; set_idx < set_ct; set_idx++) {
      writebuf[2 * set_idx - 1] = '\t';
    }
    writebuf[2 * set_ct - 1] = '\n';
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      if (marker_uidx >= chrom_end) {
	uii = get_marker_chrom_fo_idx(chrom_info_ptr, marker_uidx);
        chrom_idx = chrom_info_ptr->chrom_file_order[uii];
        chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[uii];
      }
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      bufptr = chrom_name_write(&(tbuf[1]), chrom_info_ptr, chrom_idx, zero_extra_chroms);
      *bufptr++ = '\t';
      bufptr = uint32_writex(bufptr, marker_pos[marker_uidx], '\t');
      // do not keep double-tab (if it was intentional, it should have been in
      // the header line too...)
      fwrite(tbuf, 1, bufptr - tbuf, outfile);
      bufptr = writebuf;
      cptr = cur_setting;
      for (set_idx = 0; set_idx < set_ct; set_idx++) {
        if (next_adj[set_idx] <= marker_idx) {
	  if (sip->range_ptrs[set_idx]) {
	    range_ct = sip->range_ptrs[set_idx][0];
	    uii = last_idx[set_idx];
	    if (uii < range_ct) {
	      range_start = sip->range_ptrs[set_idx][2 * uii + 1];
	      if (range_start > marker_idx) {
		*cptr = '0';
		next_adj[set_idx] = range_start;
	      } else {
		*cptr = '1';
		next_adj[set_idx] = sip->range_ptrs[set_idx][2 * uii + 2];
		last_idx[set_idx] = uii + 1;
	      }
	    } else {
              *cptr = '0';
              next_adj[set_idx] = marker_ct;
	    }
	  } else {
	    // todo
	  }
	}
	*bufptr = *cptr++;
        bufptr = &(bufptr[2]);
      }
      if (fwrite_checked(writebuf, 2 * set_ct, outfile)) {
	goto write_set_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto write_set_ret_WRITE_FAIL;
    }
    sprintf(logbuf, "Set table written to %s.\n", outname);
    logprintb();
  }
  if (sip->modifier & SET_WRITE_LIST) {
  }
  while (0) {
  write_set_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_set_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_set_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}
