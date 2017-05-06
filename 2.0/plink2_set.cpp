// This file is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink2_decompress.h"
#include "plink2_set.h"

#ifdef __cplusplus
namespace plink2 {
#endif

pglerr_t load_range_list(const chr_info_t* cip, const uint32_t* variant_bps, const char* sorted_subset_ids, const char* file_descrip, uint32_t track_set_names, uint32_t border_extend, uint32_t collapse_group, uint32_t fail_on_no_sets, uint32_t c_prefix, uint32_t allow_no_variants, uintptr_t subset_ct, uintptr_t max_subset_id_blen, gzFile gz_infile, uintptr_t* set_ct_ptr, char** set_names_ptr, uintptr_t* max_set_id_blen_ptr, uint64_t** range_sort_buf_ptr, make_set_range_t*** make_set_range_arr_ptr) {
  // In plink 1.9, called directly by extract_exclude_range(), define_sets(),
  // and indirectly by annotate(), gene_report(), and clump_reports().
  // Assumes caller will reset g_bigstack_end later.
  pglerr_t reterr = kPglRetSuccess;
  {
    g_textbuf[kMaxMediumLine - 1] = ' ';
    ll_str_t* make_set_ll = nullptr;
    char* set_names = nullptr;
    uintptr_t set_ct = 0;
    uintptr_t max_set_id_blen = 0;
    // if we need to track set names, put together a sorted list
    if (track_set_names) {
      uintptr_t line_idx = 0;
      while (gzgets(gz_infile, g_textbuf, kMaxMediumLine)) {
	++line_idx;
	if (!g_textbuf[kMaxMediumLine - 1]) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s file is pathologically long.\n", line_idx, file_descrip);
	  goto load_range_list_ret_MALFORMED_INPUT_2;
	}
	char* textbuf_first_token = skip_initial_spaces(g_textbuf);
	if (is_eoln_kns(*textbuf_first_token)) {
	  continue;
	}
	char* first_token_end = token_endnn(textbuf_first_token);
	char* cur_set_id = next_token_mult(first_token_end, 3);
	char* last_token;
	if (!collapse_group) {
	  last_token = cur_set_id;
	} else {
	  last_token = next_token(cur_set_id);
	}
	if (no_more_tokens_kns(last_token)) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s file has fewer tokens than expected.\n", line_idx, file_descrip);
	  goto load_range_list_ret_MALFORMED_INPUT_2;
	}
	const uint32_t chr_name_slen = (uintptr_t)(first_token_end - textbuf_first_token);
	*first_token_end = '\0';
	const int32_t cur_chr_code = get_chr_code(textbuf_first_token, cip, chr_name_slen);
	if (cur_chr_code < 0) {
	  sprintf(g_logbuf, "Error: Invalid chromosome code on line %" PRIuPTR " of %s file.\n", line_idx, file_descrip);
	  goto load_range_list_ret_MALFORMED_INPUT_2;
	}
	// chr_mask check removed, we want to track empty sets
	uint32_t set_id_slen = strlen_se(cur_set_id);
	cur_set_id[set_id_slen] = '\0';
	if (subset_ct) {
	  if (bsearch_str(cur_set_id, sorted_subset_ids, set_id_slen, max_subset_id_blen, subset_ct) == -1) {
	    continue;
	  }
	}
	if (collapse_group) {
	  set_id_slen = strlen_se(last_token);
	  last_token[set_id_slen] = '\0';
	}
	// when there are repeats, they are likely to be next to each other
	if (make_set_ll && (!strcmp(make_set_ll->ss, last_token))) {
	  continue;
	}
	uint32_t set_id_blen = set_id_slen + 1;
	// argh, --clump counts positional overlaps which don't include any
	// variants in the dataset.  So we prefix set IDs with a chromosome
	// index in that case (with leading zeroes) and treat cross-chromosome
	// sets as distinct.
	if (!variant_bps) {
	  set_id_blen += 4;
	}
	if (set_id_blen > max_set_id_blen) {
	  max_set_id_blen = set_id_blen;
	}
	ll_str_t* ll_tmp;
	if (bigstack_end_alloc_llstr(set_id_blen, &ll_tmp)) {
	  goto load_range_list_ret_NOMEM;
	}
	ll_tmp->next = make_set_ll;
	if (variant_bps) {
	  memcpy(ll_tmp->ss, last_token, set_id_blen);
	} else {
	  uitoa_z4((uint32_t)cur_chr_code, ll_tmp->ss);
	  // if first character of gene name is a digit, natural sort has
	  // strange effects unless we force [3] to be nonnumeric...
	  ll_tmp->ss[3] -= 15;
	  memcpy(&(ll_tmp->ss[4]), last_token, set_id_blen - 4);
	}
	make_set_ll = ll_tmp;
	++set_ct;
      }
      if (!gzeof(gz_infile)) {
	goto load_range_list_ret_READ_FAIL;
      }
      if (!set_ct) {
	if (fail_on_no_sets) {
	  if (variant_bps) {
	    if (!allow_no_variants) {
	      // okay, this is a kludge
	      logerrprint("Error: All variants excluded by --gene{-all}, since no sets were defined from\n--make-set file.\n");
	      reterr = kPglRetMalformedInput;
	      goto load_range_list_ret_1;
	    }
	  } else {
	    if (subset_ct) {
	      logerrprint("Error: No --gene-subset genes present in --gene-report file.\n");
	      reterr = kPglRetInconsistentInput;
	    } else {
	      logerrprint("Error: Empty --gene-report file.\n");
	      reterr = kPglRetMalformedInput;
	    }
	    goto load_range_list_ret_1;
	  }
	}
	LOGERRPRINTF("Warning: No valid ranges in %s file.\n", file_descrip);
	goto load_range_list_ret_1;
      }
      // c_prefix is 0 or 2
      max_set_id_blen += c_prefix;
      if (max_set_id_blen > kMaxIdBlen) {
	logerrprint("Error: Set IDs are limited to " MAX_ID_SLEN_STR " characters.\n");
	goto load_range_list_ret_MALFORMED_INPUT;
      }
      char** strptr_arr;
      if (bigstack_alloc_c(set_ct * max_set_id_blen, set_names_ptr) ||
	  bigstack_alloc_cp(set_ct, &strptr_arr)) {
	goto load_range_list_ret_NOMEM;
      }
      set_names = *set_names_ptr;
      for (uintptr_t set_idx = 0; set_idx < set_ct; ++set_idx) {
	strptr_arr[set_idx] = make_set_ll->ss;
	make_set_ll = make_set_ll->next;
      }
      strptr_arr_nsort(set_ct, strptr_arr);
      set_ct = copy_and_dedup_sorted_strptrs_to_strbox(strptr_arr, set_ct, max_set_id_blen, &(set_names[c_prefix]));
      if (c_prefix) {
	for (uintptr_t set_idx = 0; set_idx < set_ct; ++set_idx) {
	  memcpy(&(set_names[set_idx * max_set_id_blen]), "C_", 2);
	}
      }
      bigstack_shrink_top(set_names, set_ct * max_set_id_blen);
      if (gzrewind(gz_infile)) {
	goto load_range_list_ret_READ_FAIL;
      }
    } else {
      set_ct = 1;
    }
    make_set_range_t** make_set_range_arr = (make_set_range_t**)bigstack_end_alloc(set_ct * sizeof(intptr_t));
    if (!make_set_range_arr) {
      goto load_range_list_ret_NOMEM;
    }
    for (uintptr_t set_idx = 0; set_idx < set_ct; ++set_idx) {
      make_set_range_arr[set_idx] = nullptr;
    }
    uintptr_t line_idx = 0;
    uint32_t chr_start = 0;
    uint32_t chr_end = 0;
    while (gzgets(gz_infile, g_textbuf, kMaxMediumLine)) {
      ++line_idx;
      if (!g_textbuf[kMaxMediumLine - 1]) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s file is pathologically long.\n", line_idx, file_descrip);
	goto load_range_list_ret_MALFORMED_INPUT_2;
      }
      char* textbuf_first_token = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*textbuf_first_token)) {
	continue;
      }
      char* first_token_end = token_endnn(textbuf_first_token);
      char* cur_set_id = next_token_mult(first_token_end, 3);
      char* last_token;
      if (!collapse_group) {
	last_token = cur_set_id;
      } else {
	last_token = next_token(cur_set_id);
      }
      if (no_more_tokens_kns(last_token)) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s file has fewer tokens than expected.\n", line_idx, file_descrip);
	goto load_range_list_ret_MALFORMED_INPUT_2;
      }
      const uint32_t chr_name_slen = (uintptr_t)(first_token_end - textbuf_first_token);
      *first_token_end = '\0';
      const int32_t cur_chr_code = get_chr_code(textbuf_first_token, cip, chr_name_slen);
      if (cur_chr_code < 0) {
	sprintf(g_logbuf, "Error: Invalid chromosome code on line %" PRIuPTR " of %s file.\n", line_idx, file_descrip);
	goto load_range_list_ret_MALFORMED_INPUT_2;
      }
      if (!is_set(cip->chr_mask, cur_chr_code)) {
	continue;
      }
      if (variant_bps) {
	const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)cur_chr_code];
	chr_start = cip->chr_fo_vidx_start[chr_fo_idx];
	chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	if (chr_end == chr_start) {
	  continue;
	}
	// might need to move this outside the if-statement later
	if (subset_ct && (bsearch_str(cur_set_id, sorted_subset_ids, strlen_se(cur_set_id), max_subset_id_blen, subset_ct) == -1)) {
	  continue;
	}
      }
      char* textbuf_iter = skip_initial_spaces(&(first_token_end[1]));
      uint32_t range_first;
      if (scanadv_uint_defcap(&textbuf_iter, &range_first)) {
	sprintf(g_logbuf, "Error: Invalid range start position on line %" PRIuPTR " of %s file.\n", line_idx, file_descrip);
	goto load_range_list_ret_MALFORMED_INPUT_2;
      }
      textbuf_iter = next_token(textbuf_iter);
      uint32_t range_last;
      if (scanadv_uint_defcap(&textbuf_iter, &range_last)) {
	sprintf(g_logbuf, "Error: Invalid range end position on line %" PRIuPTR " of %s file.\n", line_idx, file_descrip);
	goto load_range_list_ret_MALFORMED_INPUT_2;
      }
      if (range_last < range_first) {
	sprintf(g_logbuf, "Error: Range end position smaller than range start on line %" PRIuPTR " of %s file.\n", line_idx, file_descrip);
	wordwrapb(0);
	goto load_range_list_ret_MALFORMED_INPUT_2;
      }
      if (border_extend > range_first) {
	range_first = 0;
      } else {
	range_first -= border_extend;
      }
      range_last += border_extend;
      uint32_t cur_set_idx = 0;
      if (set_ct > 1) {
	// bugfix: bsearch_str_natural requires null-terminated string
	const uint32_t last_token_slen = strlen_se(last_token);
	last_token[last_token_slen] = '\0';
	if (c_prefix) {
	  last_token = &(last_token[-2]);
	  memcpy(last_token, "C_", 2);
	} else if (!variant_bps) {
	  last_token = &(last_token[-4]);
	  uitoa_z4((uint32_t)cur_chr_code, last_token);
	  last_token[3] -= 15;
	}
	// this should never fail
	cur_set_idx = (uint32_t)bsearch_str_natural(last_token, set_names, max_set_id_blen, set_ct);
      }
      if (variant_bps) {
	// translate to within-chromosome uidx
	range_first = uint32arr_greater_than(&(variant_bps[chr_start]), chr_end - chr_start, range_first);
	range_last = uint32arr_greater_than(&(variant_bps[chr_start]), chr_end - chr_start, range_last + 1);
	if (range_last > range_first) {
	  make_set_range_t* msr_tmp = (make_set_range_t*)bigstack_end_alloc(sizeof(make_set_range_t));
	  if (!msr_tmp) {
	    goto load_range_list_ret_NOMEM;
	  }
	  msr_tmp->next = make_set_range_arr[cur_set_idx];
	  // normally, I'd keep chr_idx here since that enables by-chromosome
	  // sorting, but that's probably not worth bloating make_set_range_t
	  // from 16 to 32 bytes
	  msr_tmp->uidx_start = chr_start + range_first;
	  msr_tmp->uidx_end = chr_start + range_last;
	  make_set_range_arr[cur_set_idx] = msr_tmp;
	}
      } else {
	make_set_range_t* msr_tmp = (make_set_range_t*)bigstack_end_alloc(sizeof(make_set_range_t));
	if (!msr_tmp) {
	  goto load_range_list_ret_NOMEM;
	}
	msr_tmp->next = make_set_range_arr[cur_set_idx];
	msr_tmp->uidx_start = range_first;
	msr_tmp->uidx_end = range_last + 1;
	make_set_range_arr[cur_set_idx] = msr_tmp;
      }
    }
    if (!gzeof(gz_infile)) {
      goto load_range_list_ret_READ_FAIL;
    }
    // allocate buffer for sorting ranges later
    uint32_t max_set_range_ct = 0;
    for (uint32_t set_idx = 0; set_idx < set_ct; ++set_idx) {
      uint32_t cur_set_range_ct = 0;
      make_set_range_t* msr_tmp = make_set_range_arr[set_idx];
      while (msr_tmp) {
	++cur_set_range_ct;
	msr_tmp = msr_tmp->next;
      }
      if (cur_set_range_ct > max_set_range_ct) {
	max_set_range_ct = cur_set_range_ct;
      }
    }
    if (range_sort_buf_ptr) {
      if (bigstack_end_alloc_ull(max_set_range_ct, range_sort_buf_ptr)) {
	goto load_range_list_ret_NOMEM;
      }
    }
    if (set_ct_ptr) {
      *set_ct_ptr = set_ct;
    }
    if (max_set_id_blen_ptr) {
      *max_set_id_blen_ptr = max_set_id_blen;
    }
    *make_set_range_arr_ptr = make_set_range_arr;
  }
  while (0) {
  load_range_list_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  load_range_list_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  load_range_list_ret_MALFORMED_INPUT_2:
    logerrprintb();
  load_range_list_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 load_range_list_ret_1:
  return reterr;
}

pglerr_t extract_exclude_range(const char* fnames, const chr_info_t* cip, const uint32_t* variant_bps, uint32_t raw_variant_ct, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t orig_variant_ct = *variant_ct_ptr;
  if (!orig_variant_ct) {
    return kPglRetSuccess;
  }
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  gzFile gz_infile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uintptr_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uintptr_t* variant_include_mask = nullptr;
    if (!do_exclude) {
      if (bigstack_calloc_ul(raw_variant_ctl, &variant_include_mask)) {
	goto extract_exclude_range_ret_NOMEM;
      }
    }
    const char* fnames_iter = fnames;
    do {
      reterr = gzopen_read_checked(fnames_iter, &gz_infile);
      if (reterr) {
	goto extract_exclude_range_ret_1;
      }
      make_set_range_t** range_arr = nullptr;
      reterr = load_range_list(cip, variant_bps, nullptr, do_exclude? "--exclude range" : "--extract range", 0, 0, 0, 0, 0, 1, 0, 0, gz_infile, nullptr, nullptr, nullptr, nullptr, &range_arr);
      if (reterr) {
	goto extract_exclude_range_ret_1;
      }
      if (gzclose_null(&gz_infile)) {
	goto extract_exclude_range_ret_READ_FAIL;
      }
      make_set_range_t* msr_tmp = range_arr[0];
      if (do_exclude) {
	while (msr_tmp) {
	  clear_bits_nz(msr_tmp->uidx_start, msr_tmp->uidx_end, variant_include);
	  msr_tmp = msr_tmp->next;
	}
      } else {
	while (msr_tmp) {
	  fill_bits_nz(msr_tmp->uidx_start, msr_tmp->uidx_end, variant_include_mask);
	  msr_tmp = msr_tmp->next;
	}
      }
      fnames_iter = (const char*)rawmemchr(fnames_iter, '\0');
      ++fnames_iter;
    } while (*fnames_iter);
    if (!do_exclude) {
      bitvec_and(variant_include_mask, raw_variant_ctl, variant_include);
    }
    *variant_ct_ptr = popcount_longs(variant_include, raw_variant_ctl);
    if (*variant_ct_ptr == orig_variant_ct) {
      LOGERRPRINTF("Warning: No variants excluded by '--%s range'.\n", do_exclude? "exclude" : "extract");
    } else {
      const uint32_t excluded_ct = orig_variant_ct - (*variant_ct_ptr);
      LOGPRINTF("--%s range: %u variant%s excluded.\n", do_exclude? "exclude" : "extract", excluded_ct, (excluded_ct == 1)? "" : "s");
    }
  }
  while (0) {
  extract_exclude_range_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  extract_exclude_range_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  }
 extract_exclude_range_ret_1:
  gzclose_cond(gz_infile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}
#endif
