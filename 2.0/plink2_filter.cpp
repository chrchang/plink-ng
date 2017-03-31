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


#include "plink2_filter.h"
#include "plink2_stats.h"

#ifdef __cplusplus
namespace plink2 {
#endif

pglerr_t from_to_flag(char** variant_ids, const uint32_t* variant_id_htable, const char* varid_from, const char* varid_to, uint32_t raw_variant_ct, uintptr_t max_variant_id_blen, uintptr_t variant_id_htable_size, uintptr_t* variant_include, chr_info_t* cip, uint32_t* variant_ct_ptr) {
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t* extra_alloc_base = &(variant_id_htable[round_up_pow2_ui(variant_id_htable_size, kInt32PerCacheline)]);
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t variant_uidx_start = 0xffffffffU;
    if (varid_from) {
      uint32_t cur_llidx;
      variant_uidx_start = variant_id_dup_htable_find(varid_from, variant_ids, variant_id_htable, extra_alloc_base, strlen(varid_from), variant_id_htable_size, max_variant_id_blen, &cur_llidx);
      if (variant_uidx_start == 0xffffffffU) {
	sprintf(g_logbuf, "Error: --from variant '%s' not found.\n", varid_from);
	goto from_to_flag_ret_INCONSISTENT_INPUT_WW;
      }
      // do *not* check variant_include here.  variant ID uniqueness should not
      // be dependent on the order in which filters are applied.
      if (cur_llidx != 0xffffffffU) {
	sprintf(g_logbuf, "Error: --from variant ID '%s' appears multiple times.\n", varid_from);
	goto from_to_flag_ret_INCONSISTENT_INPUT_WW;
      }
      chr_fo_idx = get_variant_chr_fo_idx(cip, variant_uidx_start);
    }
    uint32_t variant_uidx_end = 0;
    if (varid_to) {
      uint32_t cur_llidx;
      variant_uidx_end = variant_id_dup_htable_find(varid_to, variant_ids, variant_id_htable, extra_alloc_base, strlen(varid_to), variant_id_htable_size, max_variant_id_blen, &cur_llidx);
      if (variant_uidx_end == 0xffffffffU) {
	sprintf(g_logbuf, "Error: --to variant '%s' not found.\n", varid_to);
	goto from_to_flag_ret_INCONSISTENT_INPUT_WW;
      }
      if (cur_llidx != 0xffffffffU) {
	sprintf(g_logbuf, "Error: --to variant ID '%s' appears multiple times.\n", varid_to);
	goto from_to_flag_ret_INCONSISTENT_INPUT_WW;
      }
      uint32_t chr_fo_idx2 = get_variant_chr_fo_idx(cip, variant_uidx_end);
      if (variant_uidx_start == 0xffffffffU) {
	chr_fo_idx = chr_fo_idx2;
	variant_uidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
      } else {
	if (chr_fo_idx != chr_fo_idx2) {
	  logerrprint("Error: --from and --to variants are not on the same chromosome.\n");
	  goto from_to_flag_ret_INCONSISTENT_INPUT;
	}
	if (variant_uidx_start > variant_uidx_end) {
	  // permit order to be reversed
	  uint32_t uii = variant_uidx_start;
	  variant_uidx_start = variant_uidx_end;
	  variant_uidx_end = uii;
	}
      }
      ++variant_uidx_end; // convert to half-open interval
    } else {
      variant_uidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    }
    if (variant_uidx_start) {
      clear_bits_nz(0, variant_uidx_start, variant_include);
    }
    if (variant_uidx_end < raw_variant_ct) {
      clear_bits_nz(variant_uidx_end, raw_variant_ct, variant_include);
    }
    fill_ulong_zero(kChrMaskWords, cip->chr_mask);
    set_bit(cip->chr_file_order[chr_fo_idx], cip->chr_mask);
    const uint32_t new_variant_ct = popcount_bit_idx(variant_include, variant_uidx_start, variant_uidx_end);
    LOGPRINTF("--from/--to: %u variant%s remaining.\n", new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  from_to_flag_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  from_to_flag_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  return reterr;
}

pglerr_t snp_flag(const uint32_t* variant_bp, char** variant_ids, const uint32_t* variant_id_htable, const char* varid_snp, uint32_t raw_variant_ct, uintptr_t max_variant_id_blen, uintptr_t variant_id_htable_size, uint32_t do_exclude, int32_t window_bp, uintptr_t* variant_include, chr_info_t* cip, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t* extra_alloc_base = &(variant_id_htable[round_up_pow2_ui(variant_id_htable_size, kInt32PerCacheline)]);
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uint32_t cur_llidx;
    uint32_t variant_uidx = variant_id_dup_htable_find(varid_snp, variant_ids, variant_id_htable, extra_alloc_base, strlen(varid_snp), variant_id_htable_size, max_variant_id_blen, &cur_llidx);
    if (variant_uidx == 0xffffffffU) {
      sprintf(g_logbuf, "Error: --%ssnp variant '%s' not found.\n", do_exclude? "exclude-" : "", varid_snp);
      goto snp_flag_ret_INCONSISTENT_INPUT_WW;
    }
    if (window_bp == -1) {
      // duplicates ok
      
      uintptr_t* seen_uidxs;
      // not actually necessary in --exclude-snp case, but this is still fast
      // enough relative to hash table construction that there's no point in
      // complicating the code further to conditionally optimize this out
      if (bigstack_calloc_ul(raw_variant_ctl, &seen_uidxs)) {
	goto snp_flag_ret_NOMEM;
      }
      while (1) {
	set_bit(variant_uidx, seen_uidxs);
	if (cur_llidx == 0xffffffffU) {
	  break;
	}
	variant_uidx = extra_alloc_base[cur_llidx];
	cur_llidx = extra_alloc_base[cur_llidx + 1];
      }
      if (do_exclude) {
	bitvec_andnot(seen_uidxs, raw_variant_ctl, variant_include);
      } else {
	bitvec_and(seen_uidxs, raw_variant_ctl, variant_include);
      }
    } else {
      if (cur_llidx != 0xffffffffU) {
	sprintf(g_logbuf, "Error: --%ssnp + --window central variant ID '%s' appears multiple times.\n", do_exclude? "exclude-" : "", varid_snp);
	goto snp_flag_ret_INCONSISTENT_INPUT_WW;
      }
      const uint32_t chr_fo_idx = get_variant_chr_fo_idx(cip, variant_uidx);
      const uint32_t chr_vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      const uint32_t center_bp = variant_bp[variant_uidx];
      uint32_t vidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
      if (center_bp > (uint32_t)window_bp) {
	vidx_start += uint32arr_greater_than(&(variant_bp[vidx_start]), chr_vidx_end - vidx_start, center_bp - (uint32_t)window_bp);
      }
      const uint32_t bp_end = 1 + center_bp + (uint32_t)window_bp;
      const uint32_t vidx_end = vidx_start + uint32arr_greater_than(&(variant_bp[vidx_start]), chr_vidx_end - vidx_start, bp_end);
      if (do_exclude) {
	clear_bits_nz(vidx_start, vidx_end, variant_include);
      } else {
	if (vidx_start) {
	  clear_bits_nz(0, vidx_start, variant_include);
	}
	if (vidx_end < raw_variant_ct) {
	  clear_bits_nz(vidx_end, raw_variant_ct, variant_include);
	}
	fill_ulong_zero(kChrMaskWords, cip->chr_mask);
	set_bit(cip->chr_file_order[chr_fo_idx], cip->chr_mask);
      }
    }
    const uint32_t new_variant_ct = popcount_longs(variant_include, raw_variant_ctl);
    LOGPRINTF("--%ssnp%s: %u variant%s remaining.\n", do_exclude? "exclude-" : "", (window_bp == -1)? "" : " + --window", new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  snp_flag_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  snp_flag_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetInconsistentInput;
    break;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t snps_flag(char** variant_ids, const uint32_t* variant_id_htable, const range_list_t* snps_range_list_ptr, uint32_t raw_variant_ct, uintptr_t max_variant_id_blen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t* extra_alloc_base = &(variant_id_htable[round_up_pow2_ui(variant_id_htable_size, kInt32PerCacheline)]);
    const char* varid_strbox = snps_range_list_ptr->names;
    const unsigned char* starts_range = snps_range_list_ptr->starts_range;
    const uint32_t varid_ct = snps_range_list_ptr->name_ct;
    const uintptr_t varid_max_blen = snps_range_list_ptr->name_max_blen;
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uintptr_t* seen_uidxs;
    if (bigstack_calloc_ul(raw_variant_ctl, &seen_uidxs)) {
      goto snps_flag_ret_NOMEM;
    }
    uint32_t range_start_vidx = 0xffffffffU;
    for (uint32_t varid_idx = 0; varid_idx < varid_ct; ++varid_idx) {
      const char* cur_varid = &(varid_strbox[varid_idx * varid_max_blen]);
      uint32_t cur_llidx;
      uint32_t variant_uidx = variant_id_dup_htable_find(cur_varid, variant_ids, variant_id_htable, extra_alloc_base, strlen(cur_varid), variant_id_htable_size, max_variant_id_blen, &cur_llidx);
      if (variant_uidx == 0xffffffffU) {
	sprintf(g_logbuf, "Error: --%ssnps variant '%s' not found.\n", do_exclude? "exclude-" : "", cur_varid);
	goto snps_flag_ret_INCONSISTENT_INPUT_WW;
      }
      if (starts_range[varid_idx]) {
	if (cur_llidx != 0xffffffffU) {
	  sprintf(g_logbuf, "Error: --%ssnps range-starting variant ID '%s' appears multiple times.\n", do_exclude? "exclude-" : "", cur_varid);
	  goto snps_flag_ret_INCONSISTENT_INPUT_WW;
	}
	range_start_vidx = variant_uidx;
      } else {
	if (range_start_vidx != 0xffffffffU) {
	  if (cur_llidx != 0xffffffffU) {
	    sprintf(g_logbuf, "Error: --%ssnps range-ending variant ID '%s' appears multiple times.\n", do_exclude? "exclude-" : "", cur_varid);
	    goto snps_flag_ret_INCONSISTENT_INPUT_WW;
	  }
	  if (variant_uidx < range_start_vidx) {
	    const uint32_t uii = variant_uidx;
	    variant_uidx = range_start_vidx;
	    range_start_vidx = uii;
	  }
	  fill_bits_nz(range_start_vidx, variant_uidx + 1, seen_uidxs);
	} else {
	  while (1) {
	    set_bit(variant_uidx, seen_uidxs);
	    if (cur_llidx == 0xffffffffU) {
	      break;
	    }
	    variant_uidx = extra_alloc_base[cur_llidx];
	    cur_llidx = extra_alloc_base[cur_llidx + 1];
	  }
	}
	range_start_vidx = 0xffffffffU;
      }
    }
    if (do_exclude) {
      bitvec_andnot(seen_uidxs, raw_variant_ctl, variant_include);
    } else {
      bitvec_and(seen_uidxs, raw_variant_ctl, variant_include);
    }
    const uint32_t new_variant_ct = popcount_longs(variant_include, raw_variant_ctl);
    LOGPRINTF("--%ssnps: %u variant%s remaining.\n", do_exclude? "exclude-" : "", new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  snps_flag_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  snps_flag_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetInconsistentInput;
    break;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

void extract_exclude_process_token(char** variant_ids, const uint32_t* variant_id_htable, const uint32_t* extra_alloc_base, const char* tok_start, uint32_t variant_id_htable_size, uintptr_t max_variant_id_blen, uint32_t curtoklen, uintptr_t* already_seen, uintptr_t* duplicate_ct_ptr) {
  uint32_t cur_llidx;
  uint32_t variant_uidx = variant_id_dup_htable_find(tok_start, variant_ids, variant_id_htable, extra_alloc_base, curtoklen, variant_id_htable_size, max_variant_id_blen, &cur_llidx);
  if (variant_uidx == 0xffffffffU) {
    return;
  }
  if (IS_SET(already_seen, variant_uidx)) {
    *duplicate_ct_ptr += 1;
  } else {
    while (1) {
      SET_BIT(variant_uidx, already_seen);
      if (cur_llidx == 0xffffffffU) {
	return;
      }
      variant_uidx = extra_alloc_base[cur_llidx];
      cur_llidx = extra_alloc_base[cur_llidx + 1];
    }
  }
}

pglerr_t extract_exclude_flag_norange(char** variant_ids, const uint32_t* variant_id_htable, const char* fnames, uint32_t raw_variant_ct, uintptr_t max_variant_id_blen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gz_token_stream_t gts;
  gz_token_stream_preinit(&gts);
  pglerr_t reterr = kPglRetSuccess;
  {
    // possible todo: multithreaded read/htable lookup
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uintptr_t* already_seen;
    if (bigstack_calloc_ul(raw_variant_ctl, &already_seen)) {
      goto extract_exclude_flag_norange_ret_NOMEM;
    }
    const uint32_t* extra_alloc_base = &(variant_id_htable[round_up_pow2_ui(variant_id_htable_size, kInt32PerCacheline)]);
    const char* fnames_iter = fnames;
    uintptr_t duplicate_ct = 0;
    do {
      reterr = gz_token_stream_init(fnames_iter, &gts, g_textbuf);
      if (reterr) {
	goto extract_exclude_flag_norange_ret_1;
      }
      uint32_t token_slen;
      while (1) {
	char* token_start = gz_token_stream_advance(&gts, &token_slen);
	if (!token_start) {
	  break;
	}
	extract_exclude_process_token(variant_ids, variant_id_htable, extra_alloc_base, token_start, variant_id_htable_size, max_variant_id_blen, token_slen, already_seen, &duplicate_ct);
      }
      if (token_slen) {
	// error code
	if (token_slen == 0xffffffffU) {
	  sprintf(g_logbuf, "Error: Excessively long ID in --%s file.\n", do_exclude? "exclude" : "extract");
	  goto extract_exclude_flag_norange_ret_MALFORMED_INPUT_2;
	} else {
	  goto extract_exclude_flag_norange_ret_READ_FAIL;
	}
      }
      if (gz_token_stream_close(&gts)) {
	goto extract_exclude_flag_norange_ret_READ_FAIL;
      }
      fnames_iter = (const char*)rawmemchr(fnames_iter, '\0');
      ++fnames_iter;
    } while (*fnames_iter);
    if (do_exclude) {
      bitvec_andnot(already_seen, raw_variant_ctl, variant_include);
    } else {
      bitvec_and(already_seen, raw_variant_ctl, variant_include);
    }
    const uint32_t new_variant_ct = popcount_longs(variant_include, raw_variant_ctl);
    LOGPRINTF("--%s: %u variant%s remaining.\n", do_exclude? "exclude" : "extract", new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
    if (duplicate_ct) {
      LOGERRPRINTF("Warning: At least %" PRIuPTR " duplicate ID%s in --%s file(s).\n", duplicate_ct, (duplicate_ct == 1)? "" : "s", do_exclude? "exclude" : "extract");
    }
  }
  while (0) {
  extract_exclude_flag_norange_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  extract_exclude_flag_norange_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  extract_exclude_flag_norange_ret_MALFORMED_INPUT_2:
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  }
 extract_exclude_flag_norange_ret_1:
  gz_token_stream_close(&gts);
  bigstack_reset(bigstack_mark);
  return reterr;
}

static const char keep_remove_flag_strs[4][11] = {"keep", "remove", "keep-fam", "remove-fam"};

pglerr_t keep_or_remove(const char* fnames, const char* sample_ids, const char* sids, uint32_t raw_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, keep_flags_t flags, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  const char* flag_name = keep_remove_flag_strs[flags % 4];
  gzFile gz_infile = nullptr;
  uintptr_t line_idx = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    uintptr_t* seen_uidxs;
    if (bigstack_calloc_ul(raw_sample_ctl, &seen_uidxs)) {
      goto keep_or_remove_ret_NOMEM;
    }
    uintptr_t loadbuf_size = bigstack_left();
    loadbuf_size -= loadbuf_size / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto keep_or_remove_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    
    const uint32_t families_only = flags & kfKeepFam;
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    char* idbuf = nullptr;
    uint32_t* xid_map = nullptr;
    char* sorted_xidbox = nullptr;
    uintptr_t max_xid_blen = max_sample_id_blen - 1;
    if (families_only) {
      // only need to do this once
      if (bigstack_alloc_ui(orig_sample_ct, &xid_map) ||
	  bigstack_alloc_c(orig_sample_ct * max_xid_blen, &sorted_xidbox)) {
	goto keep_or_remove_ret_NOMEM;
      }
      uint32_t sample_uidx = 0;
      for (uint32_t sample_idx = 0; sample_idx < orig_sample_ct; ++sample_idx, ++sample_uidx) {
	next_set_unsafe_ck(sample_include, &sample_uidx);
	const char* fidt_ptr = &(sample_ids[sample_uidx * max_sample_id_blen]);
	const char* fidt_end = (const char*)rawmemchr(fidt_ptr, '\t');
        const uint32_t cur_fidt_slen = 1 + (uintptr_t)(fidt_end - fidt_ptr);
	// include trailing tab, to simplify bsearch_str_lb() usage
	memcpyx(&(sorted_xidbox[sample_idx * max_xid_blen]), fidt_ptr, cur_fidt_slen, '\0');
	xid_map[sample_idx] = sample_uidx;
      }
      if (sort_strbox_indexed(orig_sample_ct, max_xid_blen, 0, sorted_xidbox, xid_map)) {
	goto keep_or_remove_ret_NOMEM;
      }
    }
    unsigned char* bigstack_mark2 = g_bigstack_base;
    const char* fnames_iter = fnames;
    uintptr_t duplicate_ct = 0;
    do {
      reterr = gzopen_read_checked(fnames_iter, &gz_infile);
      if (reterr) {
	goto keep_or_remove_ret_1;
      }
      char* loadbuf_first_token;
      xid_mode_t xid_mode;
      if (!families_only) {
	reterr = load_xid_header(flag_name, (flags & kfKeepForceSid)? kSidDetectModeForce : (sids? kSidDetectModeLoaded : kSidDetectModeNotLoaded), loadbuf_size, loadbuf, nullptr, &line_idx, &loadbuf_first_token, &gz_infile, &xid_mode);
	if (reterr) {
	  if (reterr == kPglRetEmptyFile) {
	    reterr = kPglRetSuccess;
	    goto keep_or_remove_empty_file;
	  }
	  if (reterr == kPglRetLongLine) {
	    if (loadbuf_size == kMaxLongLine) {
	      goto keep_or_remove_ret_LONG_LINE;
	    }
	    goto keep_or_remove_ret_NOMEM;
	  }
	  goto keep_or_remove_ret_1;
	}
	reterr = sorted_xidbox_init_alloc(sample_include, sample_ids, sids, orig_sample_ct, max_sample_id_blen, max_sid_blen, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
	if (reterr) {
	  goto keep_or_remove_ret_1;
	}
	if (bigstack_alloc_c(max_xid_blen, &idbuf)) {
	  goto keep_or_remove_ret_NOMEM;
	}
      } else {
	loadbuf_first_token = loadbuf;
	loadbuf[0] = '\0';
      }
      while (1) {
	if (!is_eoln_kns(*loadbuf_first_token)) {
	  if (!families_only) {
	    char* loadbuf_iter = loadbuf_first_token;
	    uint32_t sample_uidx;
	    if (!sorted_xidbox_read_find(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &loadbuf_iter, &sample_uidx, idbuf)) {
	      if (IS_SET(seen_uidxs, sample_uidx)) {
		++duplicate_ct;
	      } else {
		SET_BIT(sample_uidx, seen_uidxs);
	      }
	    } else if (!loadbuf_iter) {
	      goto keep_or_remove_ret_MISSING_TOKENS;
	    }
	  } else {
	    char* token_end = token_endnn(loadbuf_first_token);
	    *token_end = '\t';
	    const uint32_t slen = 1 + (uintptr_t)(token_end - loadbuf_first_token);
	    uint32_t lb_idx = bsearch_str_lb(loadbuf_first_token, sorted_xidbox, slen, max_xid_blen, orig_sample_ct);
	    *token_end = ' ';
	    const uint32_t ub_idx = bsearch_str_lb(loadbuf_first_token, sorted_xidbox, slen, max_xid_blen, orig_sample_ct);
	    if (ub_idx != lb_idx) {
	      uint32_t sample_uidx = xid_map[lb_idx];
	      if (IS_SET(seen_uidxs, sample_uidx)) {
		++duplicate_ct;
	      } else {
		while (1) {
		  SET_BIT(sample_uidx, seen_uidxs);
		  if (++lb_idx == ub_idx) {
		    break;
		  }
		  sample_uidx = xid_map[lb_idx];
		}
	      }
	    }
	  }
	}
	++line_idx;
	if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	  if (!gzeof(gz_infile)) {
	    goto keep_or_remove_ret_READ_FAIL;
	  }
	  goto keep_or_remove_empty_file;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  if (loadbuf_size == kMaxLongLine) {
	    goto keep_or_remove_ret_LONG_LINE;
	  }
	  goto keep_or_remove_ret_NOMEM;
	}
	loadbuf_first_token = skip_initial_spaces(loadbuf);
      }
    keep_or_remove_empty_file:
      if (gzclose_null(&gz_infile)) {
	goto keep_or_remove_ret_READ_FAIL;
      }
      bigstack_reset(bigstack_mark2);
      fnames_iter = (const char*)rawmemchr(fnames_iter, '\0');
      ++fnames_iter;
    } while (*fnames_iter);
    if (flags & kfKeepRemove) {
      bitvec_andnot(seen_uidxs, raw_sample_ctl, sample_include);
    } else {
      memcpy(sample_include, seen_uidxs, raw_sample_ctl * sizeof(intptr_t));
    }
    const uint32_t sample_ct = popcount_longs(sample_include, raw_sample_ctl);
    *sample_ct_ptr = sample_ct;
    LOGPRINTF("--%s: %u sample%s remaining.\n", flag_name, sample_ct, (sample_ct == 1)? "" : "s");
    if (duplicate_ct) {
      // "At least" since this does not count duplicate IDs absent from the
      // .fam.
      LOGERRPRINTF("Warning: At least %" PRIuPTR " duplicate ID%s in --%s file.(s)\n", duplicate_ct, (duplicate_ct == 1)? "" : "s", flag_name);
    }
  }
  while (0) {
  keep_or_remove_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  keep_or_remove_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  keep_or_remove_ret_LONG_LINE:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --%s file is pathologically long.\n", line_idx, flag_name);
    reterr = kPglRetMalformedInput;
    break;
  keep_or_remove_ret_MISSING_TOKENS:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --%s file has fewer tokens than expected.\n", line_idx, flag_name);
    reterr = kPglRetMalformedInput;
    break;
  }
 keep_or_remove_ret_1:
  bigstack_reset(bigstack_mark);
  gzclose_cond(gz_infile);
  return reterr;
}

void compute_maj_alleles_and_nonmajor_freqs(const uintptr_t* variant_include, const uintptr_t* variant_allele_idxs, const uint64_t* founder_allele_dosages, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t maf_succ, alt_allele_ct_t* maj_alleles, double* nonmajor_freqs) {
  // ok for maj_alleles or nonmajor_freqs to be nullptr
  const uintptr_t maf_succ_x2 = 65536 * maf_succ; // dosage is in 32768ths
  const uintptr_t maf_succ_x4 = 131072 * maf_succ;
  uint32_t variant_uidx = 0;
  if (variant_allele_idxs && (variant_allele_idxs[raw_variant_ct] == 2 * raw_variant_ct)) {
    variant_allele_idxs = nullptr;
  }
  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    uint32_t maj_allele_idx = 0;
    uint64_t min_dosage_x2;
    uint64_t tot_dosage_x2;
    if (!variant_allele_idxs) {
      uint64_t dosage0 = founder_allele_dosages[2 * variant_uidx];
      uint64_t min_dosage = founder_allele_dosages[2 * variant_uidx + 1];
      tot_dosage_x2 = (min_dosage + dosage0) * 2LLU;
      if (dosage0 < min_dosage) {
	min_dosage = dosage0;
	maj_allele_idx = 1;
      }
      min_dosage_x2 = min_dosage * 2;
    } else {
      const uintptr_t variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      const uint64_t* cur_founder_allele_dosages = &(founder_allele_dosages[variant_allele_idx_base]);
      const uint32_t cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
      uint64_t max_dosage_x2 = 0;
      tot_dosage_x2 = 0;
      for (uint32_t allele_idx = 0; allele_idx < cur_allele_ct; ++allele_idx) {
	uint64_t cur_dosage_x2 = 2 * cur_founder_allele_dosages[allele_idx];
	tot_dosage_x2 += cur_dosage_x2;
	if (cur_dosage_x2 > max_dosage_x2) {
	  max_dosage_x2 = cur_dosage_x2;
	  maj_allele_idx = allele_idx;
	}
      }
      min_dosage_x2 = tot_dosage_x2 - max_dosage_x2;
    }
    if (maj_alleles) {
      maj_alleles[variant_uidx] = (alt_allele_ct_t)maj_allele_idx;
    }
    if (nonmajor_freqs) {
      tot_dosage_x2 += maf_succ_x4;
      if (tot_dosage_x2) {
	const double min_ctd = (double)((int64_t)(min_dosage_x2 + maf_succ_x2));
	const double tot_ctd = (double)((int64_t)tot_dosage_x2);
	nonmajor_freqs[variant_uidx] = min_ctd / tot_ctd;
      } else {
	// occasionally has weird effects when variant is multiallelic, but
	// shouldn't be a big deal
	nonmajor_freqs[variant_uidx] = 0.5;
      }
    }
  }
}

static inline void incr_missing_row(const uintptr_t* genovec, uint32_t acc2_vec_ct, uintptr_t* missing_acc2) {
  const vul_t* genovvec = (const vul_t*)genovec;
  vul_t* missing_acc2v = (vul_t*)missing_acc2;
  const vul_t m1 = VCONST_UL(kMask5555);
  for (uint32_t vidx = 0; vidx < acc2_vec_ct; ++vidx) {
    const vul_t geno_vword = genovvec[vidx];
    const vul_t geno_vword_shifted_masked = vul_rshift(geno_vword, 1) & m1;
    missing_acc2v[vidx] = missing_acc2v[vidx] + (geno_vword & geno_vword_shifted_masked);
  }
}

static inline void incr_het_row(const uintptr_t* genovec, uint32_t acc2_vec_ct, uintptr_t* hethap_acc2) {
  const vul_t* genovvec = (const vul_t*)genovec;
  const vul_t m1 = VCONST_UL(kMask5555);
  vul_t* hethap_acc2v = (vul_t*)hethap_acc2;
  for (uint32_t vidx = 0; vidx < acc2_vec_ct; ++vidx) {
    const vul_t geno_vword = genovvec[vidx];
    const vul_t geno_vword_shifted_nmasked = (~(vul_rshift(geno_vword, 1))) & m1;
    hethap_acc2v[vidx] = hethap_acc2v[vidx] + (geno_vword & geno_vword_shifted_nmasked);
  }
}

// multithread globals
static pgen_reader_t** g_pgr_ptrs = nullptr;
static uintptr_t** g_genovecs = nullptr;
static uint32_t* g_read_variant_uidx_starts = nullptr;
static uintptr_t** g_missing_acc2 = nullptr;
static uintptr_t** g_hethap_acc2 = nullptr;

static const uintptr_t* g_variant_include = nullptr;
static const uintptr_t* g_variant_allele_idxs = nullptr;
static const chr_info_t* g_cip = nullptr;
static uintptr_t* g_sex_male_interleaved_vec = nullptr;
static uint32_t g_raw_sample_ct = 0;
static uint32_t g_cur_block_size = 0;
static uint32_t g_calc_thread_ct = 0;
static pglerr_t g_error_ret = kPglRetSuccess;

THREAD_FUNC_DECL load_sample_missing_cts_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uintptr_t* variant_include = g_variant_include;
  const uintptr_t* variant_allele_idxs = g_variant_allele_idxs;
  const chr_info_t* cip = g_cip;
  const uintptr_t* sex_male_interleaved_vec = g_sex_male_interleaved_vec;
  const uint32_t raw_sample_ct = g_raw_sample_ct;
  const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(raw_sample_ct);
  const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
  const uint32_t acc8_vec_ct = acc2_vec_ct * 4;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const int32_t x_code = cip->xymt_codes[kChrOffsetX];
  const int32_t y_code = cip->xymt_codes[kChrOffsetY];
  uintptr_t* genovec = g_genovecs[tidx];
  uintptr_t* all_acc2 = g_missing_acc2[tidx];
  uintptr_t* all_acc4 = &(all_acc2[acc2_vec_ct * kWordsPerVec]);
  uintptr_t* all_acc8 = &(all_acc4[acc4_vec_ct * kWordsPerVec]);
  uintptr_t* all_acc32 = &(all_acc8[acc8_vec_ct * kWordsPerVec]);
  uintptr_t* hethap_acc2 = g_hethap_acc2[tidx];
  uintptr_t* hethap_acc4 = &(hethap_acc2[acc2_vec_ct * kWordsPerVec]);
  uintptr_t* hethap_acc8 = &(hethap_acc4[acc4_vec_ct * kWordsPerVec]);
  uintptr_t* hethap_acc32 = &(hethap_acc8[acc8_vec_ct * kWordsPerVec]);
  uint32_t all_ct_rem3 = 3;
  uint32_t all_ct_rem15d3 = 5;
  uint32_t all_ct_rem255d15 = 17;
  uint32_t hap_ct_rem3 = 3;
  uint32_t hap_ct_rem15d3 = 5;
  uint32_t hap_ct_rem255d15 = 17;
  fill_ulong_zero(acc2_vec_ct * kWordsPerVec * 23, all_acc2);
  fill_ulong_zero(acc2_vec_ct * kWordsPerVec * 23, hethap_acc2);
  while (1) {
    pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
    const uint32_t is_last_block = g_is_last_thread_block;
    const uint32_t cur_block_size = g_cur_block_size;
    const uint32_t cur_idx_ct = (((tidx + 1) * cur_block_size) / calc_thread_ct) - ((tidx * cur_block_size) / calc_thread_ct);
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    uint32_t chr_end = 0;
    uint32_t hethap_check = 0;
    uint32_t is_diploid_x = 0;
    uint32_t is_y = 0;
    for (uint32_t cur_idx = 0; cur_idx < cur_idx_ct; ++cur_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	const uint32_t chr_fo_idx = get_variant_chr_fo_idx(cip, variant_uidx);
	const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	hethap_check = 1;
	is_diploid_x = 0;
	is_y = 0;
	if (chr_idx == x_code) {
	  is_diploid_x = !is_set(cip->haploid_mask, 0);
	} else if (chr_idx == y_code) {
	  is_y = 1;
	} else {
	  hethap_check = is_set(cip->haploid_mask, chr_idx);
	}
      }
      pglerr_t reterr;
      // todo: pgr_read_missingness_dosage()
      if ((!variant_allele_idxs) || (!hethap_check) || (variant_allele_idxs[variant_uidx + 1] - variant_allele_idxs[variant_uidx] == 2)) {
	reterr = pgr_read_allele_countvec_subset_unsafe(nullptr, nullptr, raw_sample_ct, variant_uidx, 1, pgrp, genovec);
      } else {
	// todo: multiallelic, at-least-partially-haploid case
	// countvec may not have the correct hethap information then
	reterr = kPglRetSuccess;
	assert(0);
      }
      if (reterr) {
	g_error_ret = reterr;
	break;
      }
      // trailing bits don't matter
      // zero_trailing_quaters(raw_sample_ct, genovec);
      if (is_y) {
	interleaved_mask_zero(sex_male_interleaved_vec, acc2_vec_ct, genovec);
      }
      incr_missing_row(genovec, acc2_vec_ct, all_acc2);
      if (!(--all_ct_rem3)) {
	unroll_zero_incr_2_4(acc2_vec_ct, all_acc2, all_acc4);
	all_ct_rem3 = 3;
	if (!(--all_ct_rem15d3)) {
	  unroll_zero_incr_4_8(acc4_vec_ct, all_acc4, all_acc8);
	  all_ct_rem15d3 = 5;
	  if (!(--all_ct_rem255d15)) {
	    unroll_zero_incr_8_32(acc8_vec_ct, all_acc8, all_acc32);
	    all_ct_rem255d15 = 17;
	  }
	}
      }
      if (hethap_check) {
	if (is_diploid_x) {
	  // faster to merge with incr_het_row, but shouldn't be a big deal
	  interleaved_mask_zero(sex_male_interleaved_vec, acc2_vec_ct, genovec);
	}
	incr_het_row(genovec, acc2_vec_ct, hethap_acc2);
	if (!(--hap_ct_rem3)) {
	  unroll_zero_incr_2_4(acc2_vec_ct, hethap_acc2, hethap_acc4);
	  hap_ct_rem3 = 3;
	  if (!(--hap_ct_rem15d3)) {
	    unroll_zero_incr_4_8(acc4_vec_ct, hethap_acc4, hethap_acc8);
	    hap_ct_rem15d3 = 5;
	    if (!(--hap_ct_rem255d15)) {
	      unroll_zero_incr_8_32(acc8_vec_ct, hethap_acc8, hethap_acc32);
	      hap_ct_rem255d15 = 17;
	    }
	  }
	}
      }
    }
    if (is_last_block) {
      unroll_incr_2_4(all_acc2, acc2_vec_ct, all_acc4);
      unroll_incr_4_8(all_acc4, acc4_vec_ct, all_acc8);
      unroll_incr_8_32(all_acc8, acc8_vec_ct, all_acc32);
      unroll_incr_2_4(hethap_acc2, acc2_vec_ct, hethap_acc4);
      unroll_incr_4_8(hethap_acc4, acc4_vec_ct, hethap_acc8);
      unroll_incr_8_32(hethap_acc8, acc8_vec_ct, hethap_acc32);
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

pglerr_t load_sample_missing_cts(const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uintptr_t* variant_allele_idxs, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, uint32_t* sample_missing_geno_cts, uint32_t* sample_hethap_cts) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (!variant_ct) {
      fill_uint_zero(raw_sample_ct, sample_missing_geno_cts);
      fill_uint_zero(raw_sample_ct, sample_hethap_cts);
      goto load_sample_missing_cts_ret_1;
    }
    // this doesn't seem to saturate below 35 threads
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(raw_sample_ct);
    if (bigstack_alloc_ul(acc2_vec_ct * kWordsPerVec, &g_sex_male_interleaved_vec) ||
	bigstack_alloc_ulp(calc_thread_ct, &g_missing_acc2) ||
	bigstack_alloc_ulp(calc_thread_ct, &g_hethap_acc2)) {
      goto load_sample_missing_cts_ret_NOMEM;
    }
    fill_interleaved_mask_vec(sex_male, acc2_vec_ct, g_sex_male_interleaved_vec);
    const uintptr_t acc2_alloc_cacheline_ct = DIV_UP(acc2_vec_ct * (23 * k1LU * kBytesPerVec), kCacheline);
    unsigned char* main_loadbufs[2];
    pthread_t* threads;
    uint32_t read_block_size;
    if (multithread_load_init(variant_include, raw_sample_ct, variant_ct, pgr_alloc_cacheline_ct, acc2_alloc_cacheline_ct, 0, pgfip, &calc_thread_ct, &g_genovecs, nullptr, nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto load_sample_missing_cts_ret_NOMEM;
    }
    const uintptr_t acc2_alloc = acc2_alloc_cacheline_ct * kCacheline;
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_missing_acc2[tidx] = (uintptr_t*)bigstack_alloc_raw(acc2_alloc);
      g_hethap_acc2[tidx] = (uintptr_t*)bigstack_alloc_raw(acc2_alloc);
    }
    g_variant_include = variant_include;
    g_variant_allele_idxs = variant_allele_idxs;
    g_cip = cip;
    g_raw_sample_ct = raw_sample_ct;
    g_calc_thread_ct = calc_thread_ct;

    // nearly identical to load_allele_and_geno_counts()
    logprint("Calculating sample missingness rates... ");
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;

    const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
    const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t variant_idx = 0;
    uint32_t is_last_block = 0;
    uint32_t cur_read_block_size = read_block_size;
    uint32_t next_print_variant_idx = variant_ct / 100;
    
    while (1) {
      uintptr_t cur_loaded_variant_ct = 0;
      if (!is_last_block) {
	while (read_block_idx < read_block_ct_m1) {
	  cur_loaded_variant_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
	  if (cur_loaded_variant_ct) {
	    break;
	  }
	  ++read_block_idx;
	}
	if (read_block_idx == read_block_ct_m1) {
	  cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
	  cur_loaded_variant_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
	}
	if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_loaded_variant_ct, pgfip)) {
	  if (variant_idx) {
	    join_threads2z(calc_thread_ct, 0, threads);
	    g_cur_block_size = 0;
	    error_cleanup_threads2z(load_sample_missing_cts_thread, calc_thread_ct, threads);
	  }
	  goto load_sample_missing_cts_ret_READ_FAIL;
	}
      }
      if (variant_idx) {
	join_threads2z(calc_thread_ct, is_last_block, threads);
	reterr = g_error_ret;
	if (reterr) {
	  if (!is_last_block) {
	    g_cur_block_size = 0;
	    error_cleanup_threads2z(load_sample_missing_cts_thread, calc_thread_ct, threads);
	  }
	  if (reterr == kPglRetMalformedInput) {
	    logprint("\n");
	    logerrprint("Error: Malformed .pgen file.\n");
	  }
	  goto load_sample_missing_cts_ret_1;
	}
      }
      if (!is_last_block) {
	g_cur_block_size = cur_loaded_variant_ct;
	compute_uidx_start_partition(variant_include, cur_loaded_variant_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
	for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	  g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
	  g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
	}
	is_last_block = (variant_idx + cur_loaded_variant_ct == variant_ct);
	if (spawn_threads2z(load_sample_missing_cts_thread, calc_thread_ct, is_last_block, threads)) {
	  goto load_sample_missing_cts_ret_THREAD_CREATE_FAIL;
	}
      }

      parity = 1 - parity;
      if (variant_idx == variant_ct) {
	break;
      }
      if (variant_idx >= next_print_variant_idx) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	pct = (variant_idx * 100LLU) / variant_ct;
	printf("\b\b%u%%", pct++);
	fflush(stdout);
	next_print_variant_idx = (pct * ((uint64_t)variant_ct)) / 100;
      }

      ++read_block_idx;
      variant_idx += cur_loaded_variant_ct;
      // crucially, this is independent of the pgen_reader_t block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    const uint32_t sample_ctav2 = acc2_vec_ct * kQuatersPerVec;
    const uintptr_t acc32_offset = acc2_vec_ct * (7 * k1LU * kWordsPerVec);
    uint32_t* scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[0][acc32_offset]));
    uint32_t* scrambled_hethap_cts = (uint32_t*)(&(g_hethap_acc2[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx < calc_thread_ct; ++tidx) {
      uint32_t* thread_scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii < sample_ctav2; ++uii) {
	scrambled_missing_cts[uii] += thread_scrambled_missing_cts[uii];
      }
      uint32_t* thread_scrambled_hethap_cts = (uint32_t*)(&(g_hethap_acc2[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii < sample_ctav2; ++uii) {
	scrambled_hethap_cts[uii] += thread_scrambled_hethap_cts[uii];
      }
    }
    for (uint32_t sample_uidx = 0; sample_uidx < raw_sample_ct; ++sample_uidx) {
      const uint32_t scrambled_idx = scramble_2_4_8_32(sample_uidx);
      sample_hethap_cts[sample_uidx] = scrambled_hethap_cts[scrambled_idx];
      sample_missing_geno_cts[sample_uidx] = scrambled_missing_cts[scrambled_idx];
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
  }
  while (0) {
  load_sample_missing_cts_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  load_sample_missing_cts_ret_READ_FAIL:
    reterr = kPglRetNomem;
    break;
  load_sample_missing_cts_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetNomem;
    break;
  }
 load_sample_missing_cts_ret_1:
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t mind_filter(const uint32_t* sample_missing_geno_cts, const uint32_t* sample_hethap_cts, const char* sample_ids, const char* sids, uint32_t raw_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t variant_ct, uint32_t variant_ct_y, double mind_thresh, uintptr_t* sample_include, uintptr_t* sex_male, uint32_t* sample_ct_ptr, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto mind_filter_ret_1;
    }
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    
    uint32_t max_missing_cts[2];
    mind_thresh *= 1 + kSmallEpsilon;
    max_missing_cts[0] = (int32_t)((double)((int32_t)(variant_ct - variant_ct_y)) * mind_thresh);
    max_missing_cts[1] = (int32_t)((double)((int32_t)variant_ct) * mind_thresh);
    uintptr_t* newly_excluded;
    if (bigstack_calloc_ul(raw_sample_ctl, &newly_excluded)) {
      goto mind_filter_ret_NOMEM;
    }
    uint32_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < orig_sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_unsafe_ck(sample_include, &sample_uidx);
      uint32_t cur_missing_geno_ct = sample_missing_geno_cts[sample_uidx];
      if (sample_hethap_cts) {
	cur_missing_geno_ct += sample_hethap_cts[sample_uidx];
      }
      if (cur_missing_geno_ct > max_missing_cts[IS_SET(sex_male, sample_uidx)]) {
	SET_BIT(sample_uidx, newly_excluded);
      }
    }
    const uint32_t removed_ct = popcount_longs(newly_excluded, raw_sample_ctl);
    // don't bother with allow_no_samples check here, better to have that in
    // just one place
    LOGPRINTF("%u sample%s removed due to missing genotype data (--mind).\n", removed_ct, (removed_ct == 1)? "" : "s");
    if (removed_ct) {
      bitvec_andnot(newly_excluded, raw_sample_ctl, sample_include);
      bitvec_andnot(newly_excluded, raw_sample_ctl, sex_male);
      strcpy(outname_end, ".irem");
      if (fopen_checked(outname, "w", &outfile)) {
	goto mind_filter_ret_OPEN_FAIL;
      }
      sample_uidx = 0;
      char* textbuf = g_textbuf;
      char* write_iter = textbuf;
      char* textbuf_flush = &(textbuf[kMaxMediumLine]);
      for (uint32_t sample_idx = 0; sample_idx < removed_ct; ++sample_idx, ++sample_uidx) {
	next_set_unsafe_ck(newly_excluded, &sample_uidx);
	write_iter = strcpya(write_iter, &(sample_ids[sample_uidx * max_sample_id_blen]));
	if (sids) {
	  *write_iter++ = '\t';
	  write_iter = strcpya(write_iter, &(sids[sample_uidx * max_sid_blen]));
	}
	*write_iter++ = '\n';
	if (write_iter >= textbuf_flush) {
	  if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	    goto mind_filter_ret_WRITE_FAIL;
	  }
	  write_iter = textbuf;
	}
      }
      if (write_iter != textbuf) {
	if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	  goto mind_filter_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&outfile)) {
	goto mind_filter_ret_WRITE_FAIL;
      }
      LOGPRINTFWW("ID%s written to %s .\n", (removed_ct == 1)? "" : "s", outname);
      *sample_ct_ptr -= removed_ct;
    }
  }
  while (0) {
  mind_filter_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  mind_filter_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  mind_filter_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 mind_filter_ret_1:
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

void enforce_geno_thresh(const chr_info_t* cip, const uint32_t* variant_missing_cts, const uint32_t* variant_hethap_cts, uint32_t sample_ct, uint32_t male_ct, uint32_t first_hap_uidx, double geno_thresh, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t prefilter_variant_ct = *variant_ct_ptr;
  geno_thresh *= 1 + kSmallEpsilon;
  const uint32_t missing_max_ct_nony = (int32_t)(geno_thresh * ((int32_t)sample_ct));
  const uint32_t missing_max_ct_y = (int32_t)(geno_thresh * ((int32_t)male_ct));
  uint32_t cur_missing_max_ct = missing_max_ct_nony;
  uint32_t removed_ct = 0;
  uint32_t variant_uidx = 0;
  uint32_t y_thresh = 0xffffffffU;
  uint32_t y_end = 0xffffffffU;
  int32_t y_code;
  if (xymt_exists(cip, kChrOffsetY, &y_code)) {
    const uint32_t y_chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)y_code];
    y_thresh = cip->chr_fo_vidx_start[y_chr_fo_idx];
    y_end = cip->chr_fo_vidx_start[y_chr_fo_idx + 1];
  }
  for (uint32_t variant_idx = 0; variant_idx < prefilter_variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    if (variant_uidx >= y_thresh) {
      if (variant_uidx < y_end) {
	y_thresh = y_end;
	cur_missing_max_ct = missing_max_ct_y;
      } else {
	y_thresh = 0xffffffffU;
	cur_missing_max_ct = missing_max_ct_nony;
      }
    }
    uint32_t cur_missing_ct = variant_missing_cts[variant_uidx];
    if (variant_uidx >= first_hap_uidx) {
      cur_missing_ct += variant_hethap_cts[variant_uidx - first_hap_uidx];
    }
    if (cur_missing_ct > cur_missing_max_ct) {
      CLEAR_BIT(variant_uidx, variant_include);
      ++removed_ct;
    }
  }
  LOGPRINTF("--geno: %u variant%s removed due to missing genotype data.\n", removed_ct, (removed_ct == 1)? "" : "s");
  *variant_ct_ptr -= removed_ct;
}

void enforce_hwe_thresh(const chr_info_t* cip, const uint32_t* founder_raw_geno_cts, const uint32_t* founder_x_male_geno_cts, const uint32_t* founder_x_nosex_geno_cts, const double* hwe_x_pvals, misc_flags_t misc_flags, double hwe_thresh, uint32_t nonfounders, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  if (cip->haploid_mask[0] & 1) {
    logerrprint("Warning: --hwe has no effect since entire genome is haploid.\n");
    return;
  }
  const uint32_t prefilter_variant_ct = *variant_ct_ptr;
  uint32_t x_start = 0xffffffffU;
  uint32_t x_end = 0xffffffffU;
  int32_t x_code;
  if (xymt_exists(cip, kChrOffsetX, &x_code)) {
    const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)x_code];
    x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
    x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
  }
  uint32_t x_thresh = x_start;
  const uint32_t midp = (misc_flags / kfMiscHweMidp) & 1;
  const uint32_t keep_fewhet = (misc_flags / kfMiscHweKeepFewhet) & 1;
  hwe_thresh *= 1 - kSmallEpsilon;
  uint32_t removed_ct = 0;
  uint32_t variant_uidx = 0;
  uint32_t min_obs = 0xffffffffU;
  uint32_t max_obs = 0;
  uint32_t is_x = 0;
  uint32_t male_ref_ct = 0;
  uint32_t male_alt_ct = 0;
  const double* hwe_x_pvals_iter = hwe_x_pvals;
  const double hwe_thresh_recip = (1 + 4 * kSmallEpsilon) / hwe_thresh;
  for (uint32_t variant_idx = 0; variant_idx < prefilter_variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    if (variant_uidx >= x_thresh) {
      is_x = (variant_uidx < x_end);
      if (is_x) {
	x_thresh = x_end;
      } else {
	x_thresh = 0xffffffffU;
      }
    }
    const uint32_t* cur_geno_cts = &(founder_raw_geno_cts[3 * variant_uidx]);
    uint32_t homref_ct = cur_geno_cts[0];
    uint32_t hetref_ct = cur_geno_cts[1];
    uint32_t nonref_diploid_ct = cur_geno_cts[2];
    uint32_t test_failed;
    uint32_t cur_obs_ct;
    if (!is_x) {
      cur_obs_ct = homref_ct + hetref_ct + nonref_diploid_ct;
      if (!cur_obs_ct) {
	// currently happens for chrY, chrM
	continue;
      }
      if (keep_fewhet) {
	if (hetref_ct * ((uint64_t)hetref_ct) <= (4LLU * homref_ct) * nonref_diploid_ct) {
	  // no p-value computed at all, so don't count this toward
	  // min_obs/max_obs
	  continue;
	}
      }
      if (midp) {
	test_failed = SNPHWE_midp_t(hetref_ct, homref_ct, nonref_diploid_ct, hwe_thresh);
      } else {
	test_failed = SNPHWE_t(hetref_ct, homref_ct, nonref_diploid_ct, hwe_thresh);
      }
    } else {
      if (founder_x_male_geno_cts) {
	const uint32_t* cur_male_geno_cts = &(founder_x_male_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
	male_ref_ct = cur_male_geno_cts[0];
	homref_ct -= male_ref_ct;
	hetref_ct -= cur_male_geno_cts[1];
	male_alt_ct = cur_male_geno_cts[2];
	nonref_diploid_ct -= male_alt_ct;
      }
      if (founder_x_nosex_geno_cts) {
	const uint32_t* cur_nosex_geno_cts = &(founder_x_nosex_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
	homref_ct -= cur_nosex_geno_cts[0];
	hetref_ct -= cur_nosex_geno_cts[1];
	nonref_diploid_ct -= cur_nosex_geno_cts[2];
      }
      cur_obs_ct = homref_ct + hetref_ct + nonref_diploid_ct + male_ref_ct + male_alt_ct;
      double joint_pval = *hwe_x_pvals_iter++;
      test_failed = (joint_pval < hwe_thresh);
      if (test_failed && keep_fewhet && (hetref_ct * ((uint64_t)hetref_ct) < (4LLU * homref_ct) * nonref_diploid_ct)) {
	// female-only retest
	if (joint_pval) {
	  joint_pval *= hwe_thresh_recip;
	} else {
	  // keep the variant iff female-only p-value also underflows
	  joint_pval = 2.2250738585072013e-308;
	}
	if (midp) {
	  test_failed = !SNPHWE_midp_t(hetref_ct, homref_ct, nonref_diploid_ct, joint_pval);
	} else {
	  test_failed = !SNPHWE_t(hetref_ct, homref_ct, nonref_diploid_ct, joint_pval);
	}
      }
    }
    if (test_failed) {
      CLEAR_BIT(variant_uidx, variant_include);
      ++removed_ct;
    }
    if (cur_obs_ct < min_obs) {
      min_obs = cur_obs_ct;
    }
    if (cur_obs_ct > max_obs) {
      max_obs = cur_obs_ct;
    }
  }
  if (((uint64_t)max_obs) * 9 > ((uint64_t)min_obs) * 10) {
    logerrprint("Warning: --hwe observation counts vary by more than 10%.  Consider using\n--geno, and/or applying different p-value thresholds to distinct subsets of\nyour data.\n");
  }
  LOGPRINTFWW("--hwe%s%s: %u variant%s removed due to Hardy-Weinberg exact test (%s).\n", midp? " midp" : "", keep_fewhet? " keep-fewhet" : "", removed_ct, (removed_ct == 1)? "" : "s", nonfounders? "all samples" : "founders only");
  *variant_ct_ptr -= removed_ct;
}

void enforce_minor_freq_constraints(const uintptr_t* variant_allele_idxs, const uint64_t* founder_allele_dosages, const double* nonmajor_freqs, uint32_t raw_variant_ct, double min_maf, double max_maf, uint64_t min_allele_dosage, uint64_t max_allele_dosage, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t prefilter_variant_ct = *variant_ct_ptr;  
  uint32_t variant_uidx = 0;
  uint32_t removed_ct = 0;
  if (variant_allele_idxs && (variant_allele_idxs[raw_variant_ct] == 2 * raw_variant_ct)) {
    variant_allele_idxs = nullptr;
  }
  if ((min_maf != 0.0) || (max_maf != 1.0)) {
    // defend against floating point error
    min_maf *= 1.0 - kSmallEpsilon;
    max_maf *= 1.0 + kSmallEpsilon;
  } else {
    nonmajor_freqs = nullptr;
  }
  const uint32_t dosage_filter = min_allele_dosage || (max_allele_dosage != (~0LLU));
  
  for (uint32_t variant_idx = 0; variant_idx < prefilter_variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    if (nonmajor_freqs) {
      const double cur_nonmajor_freq = nonmajor_freqs[variant_uidx];
      if ((cur_nonmajor_freq < min_maf) || (cur_nonmajor_freq > max_maf)) {
	CLEAR_BIT(variant_uidx, variant_include);
	++removed_ct;
	continue;
      }
    }
    if (dosage_filter) {
      uint64_t cur_min_dosage;
      if (!variant_allele_idxs) {
	cur_min_dosage = MINV(founder_allele_dosages[2 * variant_uidx], founder_allele_dosages[2 * variant_uidx + 1]);
      } else {
	const uintptr_t allele_idx_end = variant_allele_idxs[variant_uidx + 1];
	uintptr_t allele_idx = variant_allele_idxs[variant_uidx];
	uint64_t max_dosage = 0;
	uint64_t tot_dosage = 0;
	do {
	  const uint64_t cur_dosage = founder_allele_dosages[allele_idx];
	  tot_dosage += cur_dosage;
	  if (cur_dosage > max_dosage) {
	    max_dosage = cur_dosage;
	  }
	} while (++allele_idx < allele_idx_end);
	cur_min_dosage = tot_dosage - max_dosage;
      }
      if ((cur_min_dosage < min_allele_dosage) || (cur_min_dosage > max_allele_dosage)) {
	CLEAR_BIT(variant_uidx, variant_include);
	++removed_ct;
      }
    }
  }
  LOGPRINTFWW("%u variant%s removed due to minor allele threshold(s) (--maf/--max-maf/--mac/--max-mac).\n", removed_ct, (removed_ct == 1)? "" : "s");
  *variant_ct_ptr -= removed_ct;
}

#ifdef __cplusplus
} // namespace plink2
#endif
