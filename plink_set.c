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
  uint32_t chrom_idx;
  uint32_t first_pos;
  uint32_t last_pos;
} Make_set_range;

int32_t define_sets(Set_info* sip, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, uintptr_t* marker_exclude_ct_ptr, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr) {
  logprint("Error: --set and --make-set are currently under development.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
  /*
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
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  Ll_str* ll_tmp;
  Make_set_range* msr_tmp;
  uint32_t* marker_id_map;
  uint32_t* set_range_ptrs;
  uint32_t* set_bounds;
  uintptr_t* marker_exclude_new;
  uintptr_t* marker_bitfield_tmp;
  uintptr_t* set_bitfield_ptrs;
  uintptr_t* include_out_of_bounds;
  uint32_t* uiptr;
  uintptr_t set_idx;
  uintptr_t bufsize;
  uintptr_t topsize_bak;
  uintptr_t ulii;
  uint32_t chrom_idx;
  uint32_t range_first;
  uint32_t range_last;
  uint32_t slen;
  uint32_t uii;
  uint32_t ujj;
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
	slen = strlen(bufptr);
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
	if (!is_set(chrom_info_ptr->chrom_mask, ii)) {
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
	  logprint("Error: All variants excluded by --gene[-all], since no sets were defined from\n--make-set file.\n");
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
      max_set_id_len = strlen(sip->merged_set_name) + 1;
      set_ct = 1;
      wkspace_left -= topsize;
      if (wkspace_alloc_c_checked(&set_names, max_set_id_len)) {
	goto define_sets_ret_NOMEM2;
      }
      wkspace_left += topsize;
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
      msr_tmp = (Make_set_range*)top_alloc(&topsize, sizeof(Make_set_range));
      msr_tmp->next = make_set_range_arr[set_idx];
      msr_tmp->chrom_idx = chrom_idx;
      msr_tmp->first_pos = range_first;
      msr_tmp->last_pos = range_last;
      make_set_range_arr[set_idx] = msr_tmp;
    }
  }
  // 3. if --gene or --gene-all is present, pre-filter variants.
  if (gene_all || sip->genekeep_flattened) {
    marker_bitfield_tmp = (uintptr_t*)top_alloc(&topsize, unfiltered_marker_ctl * sizeof(intptr_t));
    if (!marker_bitfield_tmp) {
      goto define_sets_ret_NOMEM;
    }
    marker_exclude_new = (uintptr_t*)top_alloc(&topsize, unfiltered_marker_ctl * sizeof(intptr_t));
    if (!marker_exclude_new) {
      goto define_sets_ret_NOMEM;
    }
    fill_all_bits(marker_exclude_new, unfiltered_marker_ct);
    // then include every variant that appears, or include every variant that
    // fails to appear in a fully loaded set in the complement case
    if (make_set) {
      fill_ulong_zero(marker_bitfield_tmp, unfiltered_marker_ctl);
      for (set_idx = 0; set_idx < set_ct; set_idx++) {
        msr_tmp = make_set_range_arr[set_idx];
	while (msr_tmp) {
	  ujj = msr_tmp->chrom_idx;
          uii = chrom_info_ptr->chrom_start[ujj];
	  ujj = chrom_info_ptr->chrom_end[ujj] - uii; // chrom size
	  if (ujj) {
	    uiptr = &(marker_pos[uii]);
	    uii = uint32arr_greater_than(uiptr, ujj, msr_tmp->first_pos);
	    fill_bits(marker_bitfield_tmp, uii, uint32arr_greater_than(uiptr, ujj, msr_tmp->last_pos + 1) - uii);
	    msr_tmp = msr_tmp->next;
	  }
	}
        if (complement_sets) {
	  bitfield_and(marker_exclude_new, marker_bitfield_tmp, unfiltered_marker_ctl);
          fill_ulong_zero(marker_bitfield_tmp, unfiltered_marker_ctl);
	}
      }
      if (!complement_sets) {
	bitfield_andnot(marker_exclude_new, marker_bitfield_tmp, unfiltered_marker_ctl);
      }
      bitfield_or(marker_exclude, marker_exclude_new, unfiltered_marker_ctl);
    } else {
      topsize_bak = topsize;
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
            in_set = 0;
	  } else if (!in_set) {
	    if (curtoklen >= max_set_id_len) {
	      max_set_id_len = curtoklen + 1;
	    }
	    set_ct++;
	    in_set = 1;
	  } else {
	    *bufptr2 = '\0';
	    ii = bsearch_str(bufptr, sorted_marker_ids, max_marker_id_len, 0, marker_ct - 1);
	    if (ii != -1) {
              = marker_id_map[(uint32_t)ii];
	    // look up marker ID
	    }
	  }
	}
      }
      if (!feof(infile)) {
	goto define_sets_ret_READ_FAIL;
      }
      topsize = topsize_bak;
    }
    if (marker_exclude_ct == unfiltered_marker_ct) {
      goto define_sets_ret_ALL_MARKERS_EXCLUDED;
    }
    *marker_exclude_ct_ptr = marker_exclude_ct;
    rewind(infile);
  } else if () {
    // 4. otherwise, with --set and no --set-collapse-all, count number of sets
    //    and max_name_len.
  }
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
    retval = RET_INVALID_FORMAT:
    break;
  define_sets_ret_INVALID_FORMAT_2:
    logprintb();
  define_sets_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 define_sets_ret_1:
  return retval;
  */
}
