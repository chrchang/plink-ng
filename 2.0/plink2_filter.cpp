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
#include "plink2_filter.h"
#include "plink2_stats.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void init_cmp_expr(cmp_expr_t* cmp_expr_ptr) {
  cmp_expr_ptr->pheno_name = nullptr;
}

void cleanup_cmp_expr(cmp_expr_t* cmp_expr_ptr) {
  free_cond(cmp_expr_ptr->pheno_name);
}

pglerr_t from_to_flag(char** variant_ids, const uint32_t* variant_id_htable, const char* varid_from, const char* varid_to, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uintptr_t* variant_include, chr_info_t* cip, uint32_t* variant_ct_ptr) {
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t* htable_dup_base = &(variant_id_htable[round_up_pow2(variant_id_htable_size, kInt32PerCacheline)]);
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t variant_uidx_start = 0xffffffffU;
    if (varid_from) {
      uint32_t cur_llidx;
      variant_uidx_start = variant_id_dup_htable_find(varid_from, variant_ids, variant_id_htable, htable_dup_base, strlen(varid_from), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
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
      variant_uidx_end = variant_id_dup_htable_find(varid_to, variant_ids, variant_id_htable, htable_dup_base, strlen(varid_to), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
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

pglerr_t snp_flag(const uint32_t* variant_bps, char** variant_ids, const uint32_t* variant_id_htable, const char* varid_snp, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, int32_t window_bp, uintptr_t* variant_include, chr_info_t* cip, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t* htable_dup_base = &(variant_id_htable[round_up_pow2(variant_id_htable_size, kInt32PerCacheline)]);
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uint32_t cur_llidx;
    uint32_t variant_uidx = variant_id_dup_htable_find(varid_snp, variant_ids, variant_id_htable, htable_dup_base, strlen(varid_snp), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
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
	variant_uidx = htable_dup_base[cur_llidx];
	cur_llidx = htable_dup_base[cur_llidx + 1];
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
      const uint32_t center_bp = variant_bps[variant_uidx];
      uint32_t vidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
      if (center_bp > (uint32_t)window_bp) {
	vidx_start += uint32arr_greater_than(&(variant_bps[vidx_start]), chr_vidx_end - vidx_start, center_bp - (uint32_t)window_bp);
      }
      const uint32_t bp_end = 1 + center_bp + (uint32_t)window_bp;
      const uint32_t vidx_end = vidx_start + uint32arr_greater_than(&(variant_bps[vidx_start]), chr_vidx_end - vidx_start, bp_end);
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

pglerr_t snps_flag(char** variant_ids, const uint32_t* variant_id_htable, const range_list_t* snps_range_list_ptr, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t* htable_dup_base = &(variant_id_htable[round_up_pow2(variant_id_htable_size, kInt32PerCacheline)]);
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
      uint32_t variant_uidx = variant_id_dup_htable_find(cur_varid, variant_ids, variant_id_htable, htable_dup_base, strlen(cur_varid), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
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
	    variant_uidx = htable_dup_base[cur_llidx];
	    cur_llidx = htable_dup_base[cur_llidx + 1];
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

void extract_exclude_process_token(char** variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* tok_start, uint32_t variant_id_htable_size, uintptr_t max_variant_id_slen, uint32_t token_slen, uintptr_t* already_seen, uintptr_t* duplicate_ct_ptr) {
  uint32_t cur_llidx;
  uint32_t variant_uidx = variant_id_dup_htable_find(tok_start, variant_ids, variant_id_htable, htable_dup_base, token_slen, variant_id_htable_size, max_variant_id_slen, &cur_llidx);
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
      variant_uidx = htable_dup_base[cur_llidx];
      cur_llidx = htable_dup_base[cur_llidx + 1];
    }
  }
}

pglerr_t extract_exclude_flag_norange(char** variant_ids, const uint32_t* variant_id_htable, const char* fnames, uint32_t raw_variant_ct, uintptr_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
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
    const uint32_t* htable_dup_base = &(variant_id_htable[round_up_pow2(variant_id_htable_size, kInt32PerCacheline)]);
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
	extract_exclude_process_token(variant_ids, variant_id_htable, htable_dup_base, token_start, variant_id_htable_size, max_variant_id_slen, token_slen, already_seen, &duplicate_ct);
      }
      if (token_slen) {
	// error code
	if (token_slen == 0xffffffffU) {
	  sprintf(g_logbuf, "Error: Excessively long ID in --%s file.\n", do_exclude? "exclude" : "extract");
	  goto extract_exclude_flag_norange_ret_MALFORMED_INPUT_2;
	}
	goto extract_exclude_flag_norange_ret_READ_FAIL;
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

pglerr_t require_pheno(const pheno_col_t* pheno_cols, const char* pheno_names, char* require_pheno_flattened, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  pglerr_t reterr = kPglRetSuccess;
  {
    uint32_t required_pheno_ct = 0;
    uintptr_t max_required_pheno_blen = 2;
    uintptr_t* matched_phenos = nullptr;
    char* sorted_required_pheno_names = nullptr;
    if (require_pheno_flattened) {
      char** strptr_arr = (char**)bigstack_end_mark;
      if (count_and_measure_multistr_reverse_alloc(require_pheno_flattened, bigstack_left() / sizeof(intptr_t), &required_pheno_ct, &max_required_pheno_blen, &strptr_arr)) {
	goto require_pheno_ret_NOMEM;
      }
      if ((uintptr_t)(((unsigned char*)strptr_arr) - g_bigstack_base) < required_pheno_ct * max_required_pheno_blen) {
	goto require_pheno_ret_NOMEM;
      }
      strptr_arr_sort(required_pheno_ct, strptr_arr);
      sorted_required_pheno_names = (char*)g_bigstack_base;
      required_pheno_ct = copy_and_dedup_sorted_strptrs_to_strbox(strptr_arr, required_pheno_ct, max_required_pheno_blen, sorted_required_pheno_names);
      bigstack_end_reset(bigstack_end_mark);
      bigstack_finalize_c(sorted_required_pheno_names, required_pheno_ct * max_required_pheno_blen);
      if (bigstack_calloc_ul(1 + (required_pheno_ct / kBitsPerWord), &matched_phenos)) {
	goto require_pheno_ret_NOMEM;
      }
    } else {
      if (!pheno_ct) {
	logerrprint(is_covar? "Warning: No covariates loaded; ignoring --require-covar.\n" : "Warning: No phenotypes loaded; ignoring --require-pheno.\n");
	goto require_pheno_ret_1;
      }
      required_pheno_ct = pheno_ct;
    }
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      if (sorted_required_pheno_names) {
	const char* cur_pheno_name = &(pheno_names[pheno_idx * max_pheno_name_blen]);
	const int32_t ii = bsearch_str(cur_pheno_name, sorted_required_pheno_names, strlen(cur_pheno_name), max_required_pheno_blen, required_pheno_ct);
	if (ii == -1) {
	  continue;
	}
	set_bit(ii, matched_phenos);
      }
      bitvec_and(pheno_cols[pheno_idx].nonmiss, raw_sample_ctl, sample_include);
    }
    if (matched_phenos) {
      const uint32_t first_unmatched_idx = next_unset_unsafe(matched_phenos, 0);
      if (first_unmatched_idx < required_pheno_ct) {
	LOGERRPRINTFWW("Error: --require-%s '%s' not loaded.\n", is_covar? "covar covariate" : "pheno phenotype", &(sorted_required_pheno_names[first_unmatched_idx * max_required_pheno_blen]));
	goto require_pheno_ret_INCONSISTENT_INPUT;
      }
    }
    const uint32_t new_sample_ct = popcount_longs(sample_include, raw_sample_ctl);
    const uint32_t removed_sample_ct = (*sample_ct_ptr) - new_sample_ct;
    LOGPRINTF("--require-%s: %u sample%s removed.\n", is_covar? "covar" : "pheno", removed_sample_ct, (removed_sample_ct == 1)? "" : "s");
    *sample_ct_ptr = new_sample_ct;
  }
  while (0) {
  require_pheno_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  require_pheno_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 require_pheno_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

pglerr_t keep_remove_if(const cmp_expr_t* cmp_expr, const pheno_col_t* pheno_cols, const char* pheno_names, const pheno_col_t* covar_cols, const char* covar_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t affection_01, uint32_t is_remove, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    const char* cur_name = cmp_expr->pheno_name;
    const uintptr_t name_blen = 1 + strlen(cur_name);
    const pheno_col_t* cur_pheno_col = nullptr;
    if (name_blen <= max_pheno_name_blen) {
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	if (!memcmp(cur_name, &(pheno_names[pheno_idx * max_pheno_name_blen]), name_blen)) {
	  cur_pheno_col = &(pheno_cols[pheno_idx]);
	  break;
	}
      }
    }
    if (!cur_pheno_col) {
      if (name_blen <= max_covar_name_blen) {
	for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx) {
	  if (!memcmp(cur_name, &(covar_names[covar_idx * max_covar_name_blen]), name_blen)) {
	    cur_pheno_col = &(covar_cols[covar_idx]);
	    break;
	  }
	}
      }
    }
    if (!cur_pheno_col) {
      sprintf(g_logbuf, "Error: --%s-if phenotype/covariate not loaded.\n", is_remove? "remove" : "keep");
      goto keep_remove_if_ret_INCONSISTENT_INPUT_2;
    }
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    cmp_binary_op_t binary_op = cmp_expr->binary_op;
    const uint32_t pheno_must_exist = is_remove ^ (binary_op != kCmpOperatorNoteq);
    const uintptr_t* pheno_nm = cur_pheno_col->nonmiss;
    if (pheno_must_exist) {
      bitvec_and(pheno_nm, raw_sample_ctl, sample_include);
    }
    uintptr_t* sample_include_intersect;
    if (bigstack_alloc_ul(raw_sample_ctl, &sample_include_intersect)) {
      goto keep_remove_if_ret_NOMEM;
    }
    memcpy(sample_include_intersect, sample_include, raw_sample_ctl * sizeof(intptr_t));
    if (!pheno_must_exist) {
      bitvec_and(pheno_nm, raw_sample_ctl, sample_include_intersect);
    }
    const uint32_t sample_intersect_ct = popcount_longs(sample_include_intersect, raw_sample_ctl);
    const char* cur_val_str = &(cur_name[name_blen]);
    const uint32_t val_slen = strlen(cur_val_str);
    if (cur_pheno_col->type_code == kPhenoDtypeQt) {
      double val;
      if (!scanadv_double((char*)cur_val_str, &val)) {
	sprintf(g_logbuf, "Error: Invalid --%s-if value (number expected).\n", is_remove? "remove" : "keep");
	goto keep_remove_if_ret_INCONSISTENT_INPUT_2;
      }
      if (is_remove) {
	binary_op = (cmp_binary_op_t)(kCmpOperatorEq - (uint32_t)binary_op);
      }
      const double* pheno_vals = cur_pheno_col->data.qt;
      uint32_t sample_uidx = 0;
      switch (binary_op) {
      case kCmpOperatorNoteq:
	for (uint32_t sample_idx = 0; sample_idx < sample_intersect_ct; ++sample_idx, ++sample_uidx) {
	  next_set_unsafe_ck(sample_include_intersect, &sample_uidx);
	  if (pheno_vals[sample_uidx] == val) {
	    clear_bit(sample_uidx, sample_include);
	  }
	}
	break;
      case kCmpOperatorLe:
	for (uint32_t sample_idx = 0; sample_idx < sample_intersect_ct; ++sample_idx, ++sample_uidx) {
	  next_set_unsafe_ck(sample_include_intersect, &sample_uidx);
	  if (pheno_vals[sample_uidx] >= val) {
	    clear_bit(sample_uidx, sample_include);
	  }
	}
	break;
      case kCmpOperatorLeq:
	for (uint32_t sample_idx = 0; sample_idx < sample_intersect_ct; ++sample_idx, ++sample_uidx) {
	  next_set_unsafe_ck(sample_include_intersect, &sample_uidx);
	  if (pheno_vals[sample_uidx] > val) {
	    clear_bit(sample_uidx, sample_include);
	  }
	}
	break;
      case kCmpOperatorGe:
	for (uint32_t sample_idx = 0; sample_idx < sample_intersect_ct; ++sample_idx, ++sample_uidx) {
	  next_set_unsafe_ck(sample_include_intersect, &sample_uidx);
	  if (pheno_vals[sample_uidx] <= val) {
	    clear_bit(sample_uidx, sample_include);
	  }
	}
	break;
      case kCmpOperatorGeq:
	for (uint32_t sample_idx = 0; sample_idx < sample_intersect_ct; ++sample_idx, ++sample_uidx) {
	  next_set_unsafe_ck(sample_include_intersect, &sample_uidx);
	  if (pheno_vals[sample_uidx] < val) {
	    clear_bit(sample_uidx, sample_include);
	  }
	}
	break;
      case kCmpOperatorEq:
	for (uint32_t sample_idx = 0; sample_idx < sample_intersect_ct; ++sample_idx, ++sample_uidx) {
	  next_set_unsafe_ck(sample_include_intersect, &sample_uidx);
	  if (pheno_vals[sample_uidx] != val) {
	    clear_bit(sample_uidx, sample_include);
	  }
	}
	break;
      }
    } else {
      if ((binary_op != kCmpOperatorNoteq) && (binary_op != kCmpOperatorEq)) {
	sprintf(g_logbuf, "Error: --%s-if operator type mismatch (binary and categorical phenotypes only support == and !=).\n", is_remove? "remove" : "keep");
	goto keep_remove_if_ret_INCONSISTENT_INPUT_WW;
      }
      if (cur_pheno_col->type_code == kPhenoDtypeCc) {
	uint32_t val_12 = 0; // 1 = control, 2 = case
	if (val_slen == 1) {
	  val_12 = affection_01 + (uint32_t)((unsigned char)cur_val_str[0]) - 48;
	  if ((val_12 != 1) && (val_12 != 2)) {
	    val_12 = 0;
	  } 
	} else if (val_slen == 4) {
	  if (match_upper_counted(cur_val_str, "CASE", 4)) {
	    val_12 = 2;
	  } else if (match_upper_counted(cur_val_str, "CTRL", 4)) {
	    val_12 = 1;
	  }
	} else if (val_slen == 7) {
	  if (match_upper_counted(cur_val_str, "CONTROL", 7)) {
	    val_12 = 1;
	  }
	}
	if (!val_12) {
	  sprintf(g_logbuf, "Error: Invalid --%s-if value ('case'/'%c' or 'control'/'ctrl'/'%c' expected).\n", is_remove? "remove" : "keep", '2' - affection_01, '1' - affection_01);
	  goto keep_remove_if_ret_INCONSISTENT_INPUT_WW;
	}
	if (is_remove ^ (val_12 == 2)) {
	  bitvec_and(cur_pheno_col->data.cc, raw_sample_ctl, sample_include);
	} else {
	  bitvec_andnot(cur_pheno_col->data.cc, raw_sample_ctl, sample_include);
	}
      } else {
        assert(cur_pheno_col->type_code == kPhenoDtypeCat);
	const uint32_t nonnull_cat_ct = cur_pheno_col->nonnull_category_ct;
	uint32_t cat_idx = 1;
	for (; cat_idx <= nonnull_cat_ct; ++cat_idx) {
	  if (!strcmp(cur_val_str, cur_pheno_col->category_names[cat_idx])) {
	    break;
	  }
	}
	if (cat_idx == nonnull_cat_ct + 1) {
	  double dxx;
	  if (scanadv_double((char*)cur_val_str, &dxx)) {
	    sprintf(g_logbuf, "Error: Invalid --%s-if value (category name expected).\n", is_remove? "remove" : "keep");
	    goto keep_remove_if_ret_INCONSISTENT_INPUT_2;
	  }
	  // tolerate this, there are legitimate reasons for empty categories
	  // to exist
	  LOGERRPRINTFWW("Warning: Categorical phenotype/covariate '%s' does not have a category named '%s'.\n", cur_name, cur_val_str);
	  if (pheno_must_exist) {
	    fill_ulong_zero(raw_sample_ctl, sample_include);
	  }
	} else {
	  const uint32_t* cur_cats = cur_pheno_col->data.cat;
	  uint32_t sample_uidx = 0;
	  if (pheno_must_exist) {
	    for (uint32_t sample_idx = 0; sample_idx < sample_intersect_ct; ++sample_idx, ++sample_uidx) {
	      next_set_unsafe_ck(sample_include_intersect, &sample_uidx);
	      if (cur_cats[sample_uidx] != cat_idx) {
		clear_bit(sample_uidx, sample_include);
	      }
	    }
	  } else {
	    for (uint32_t sample_idx = 0; sample_idx < sample_intersect_ct; ++sample_idx, ++sample_uidx) {
	      next_set_unsafe_ck(sample_include_intersect, &sample_uidx);
	      if (cur_cats[sample_uidx] == cat_idx) {
		clear_bit(sample_uidx, sample_include);
	      }
	    }
	  }
	}
      }
    }

    const uint32_t new_sample_ct = popcount_longs(sample_include, raw_sample_ctl);
    const uint32_t removed_sample_ct = (*sample_ct_ptr) - new_sample_ct;
    LOGPRINTF("--%s-if: %u sample%s removed.\n", is_remove? "remove" : "keep", removed_sample_ct, (removed_sample_ct == 1)? "" : "s");
    *sample_ct_ptr = new_sample_ct;
  }
  while (0) {
  keep_remove_if_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  keep_remove_if_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
  keep_remove_if_ret_INCONSISTENT_INPUT_2:
    logerrprintb();
    reterr = kPglRetInconsistentInput;
    break;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t keep_remove_cats_internal(const pheno_col_t* cur_pheno_col, const char* cats_fname, const char* cat_names_flattened, uint32_t raw_sample_ct, uint32_t is_remove, uint32_t max_thread_ct, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gz_token_stream_t gts;
  gz_token_stream_preinit(&gts);
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    const uint32_t cat_ct = cur_pheno_col->nonnull_category_ct + 1;
    const uint32_t cat_ctl = BITCT_TO_WORDCT(cat_ct);
    uintptr_t* affected_samples;
    uintptr_t* cat_include;
    if (bigstack_calloc_ul(raw_sample_ctl, &affected_samples) ||
        bigstack_alloc_ul(cat_ctl, &cat_include)) {
      goto keep_remove_cats_internal_ret_NOMEM;
    }
    fill_all_bits(cat_ct, cat_include);
    char** category_names = cur_pheno_col->category_names;
    uint32_t* cat_id_htable;
    uint32_t id_htable_size;
    reterr = alloc_and_populate_id_htable_mt(cat_include, category_names, cat_ct, max_thread_ct, &cat_id_htable, nullptr, &id_htable_size);
    if (reterr) {
      goto keep_remove_cats_internal_ret_1;
    }
    fill_ulong_zero(cat_ctl, cat_include);
    if (cats_fname) {
      reterr = gz_token_stream_init(cats_fname, &gts, g_textbuf);
      if (reterr) {
	goto keep_remove_cats_internal_ret_1;
      }
      uintptr_t skip_ct = 0;
      uint32_t token_slen;
      while (1) {
	char* token_start = gz_token_stream_advance(&gts, &token_slen);
	if (!token_start) {
	  break;
	}
	token_start[token_slen] = '\0';
        const uint32_t cur_cat_idx = id_htable_find(token_start, category_names, cat_id_htable, token_slen, id_htable_size);
	if (cur_cat_idx == 0xffffffffU) {
	  ++skip_ct;
	} else {
	  set_bit(cur_cat_idx, cat_include);
	}
      }
      if (token_slen) {
	// error code
	if (token_slen == 0xffffffffU) {
	  sprintf(g_logbuf, "Error: Excessively long ID in --%s-cats file.\n", is_remove? "remove" : "keep");
	  goto keep_remove_cats_internal_ret_MALFORMED_INPUT_2;
	}
	goto keep_remove_cats_internal_ret_READ_FAIL;
      }
      if (gz_token_stream_close(&gts)) {
	goto keep_remove_cats_internal_ret_READ_FAIL;
      }
      if (skip_ct) {
	LOGERRPRINTF("Warning: %" PRIuPTR " --%s-cats categor%s not present.\n", skip_ct, is_remove? "remove" : "keep", (skip_ct == 1)? "y" : "ies");
      }
    }
    if (cat_names_flattened) {
      uint32_t skip_ct = 0;
      const char* cat_names_iter = cat_names_flattened;
      do {
	const uint32_t cat_name_slen = strlen(cat_names_iter);
	const uint32_t cur_cat_idx = id_htable_find(cat_names_iter, category_names, cat_id_htable, cat_name_slen, id_htable_size);
	if (cur_cat_idx == 0xffffffffU) {
	  ++skip_ct;
	} else {
	  set_bit(cur_cat_idx, cat_include);
	}
	cat_names_iter = &(cat_names_iter[cat_name_slen + 1]);
      } while (*cat_names_iter);
      if (skip_ct) {
	LOGERRPRINTF("Warning: %u --%s-cat-names categor%s not present.\n", skip_ct, is_remove? "remove" : "keep", (skip_ct == 1)? "y" : "ies");
      }
    }
    const uint32_t selected_cat_ct = popcount_longs(cat_include, cat_ctl);
    if (!selected_cat_ct) {
      LOGERRPRINTF("Warning: No matching --%s-cat-names category names.\n", is_remove? "remove-cats/--remove" : "keep-cats/--keep");
    } else {
      const uint32_t* cur_cats = cur_pheno_col->data.cat;
      const uint32_t orig_sample_ct = *sample_ct_ptr;
      uint32_t sample_uidx = 0;
      for (uint32_t sample_idx = 0; sample_idx < orig_sample_ct; ++sample_idx, ++sample_uidx) {
	next_set_unsafe_ck(sample_include, &sample_uidx);
	const uint32_t cur_cat_idx = cur_cats[sample_uidx];
	if (is_set(cat_include, cur_cat_idx)) {
	  set_bit(sample_uidx, affected_samples);
	}
      }
      if (is_remove) {
	bitvec_andnot(affected_samples, raw_sample_ctl, sample_include);
      } else {
	bitvec_and(affected_samples, raw_sample_ctl, sample_include);
      }
      const uint32_t new_sample_ct = popcount_longs(sample_include, raw_sample_ctl);
      const uint32_t removed_sample_ct = (*sample_ct_ptr) - new_sample_ct;
      LOGPRINTFWW("--%s-cat-names: %u categor%s selected, %u sample%s removed.\n", is_remove? "remove-cats/--remove" : "keep-cats/--keep", selected_cat_ct, (selected_cat_ct == 1)? "y" : "ies", removed_sample_ct, (removed_sample_ct == 1)? "" : "s");
      *sample_ct_ptr = new_sample_ct;
    }
  }
  while (0) {
  keep_remove_cats_internal_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  keep_remove_cats_internal_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  keep_remove_cats_internal_ret_MALFORMED_INPUT_2:
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  }
 keep_remove_cats_internal_ret_1:
  gz_token_stream_close(&gts);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t keep_remove_cats(const char* cats_fname, const char* cat_names_flattened, const char* cat_phenoname, const pheno_col_t* pheno_cols, const char* pheno_names, const pheno_col_t* covar_cols, const char* covar_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t is_remove, uint32_t max_thread_ct, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  pglerr_t reterr = kPglRetSuccess;
  {
    if (!cat_phenoname) {
      // Default behavior:
      // 1. If at least one categorical phenotype exists, fail on >= 2, select
      //    it if one.
      // 2. Otherwise, fail if 0 or >= 2 categorical covariates, select the
      //    categorical covariate if there's exactly one.
      uint32_t cat_pheno_idx = 0xffffffffU;
      const pheno_col_t* cur_pheno_col = nullptr;
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	if (pheno_cols[pheno_idx].type_code == kPhenoDtypeCat) {
	  if (cat_pheno_idx != 0xffffffffU) {
	    sprintf(g_logbuf, "Error: Multiple categorical phenotypes present. Use --%s-cat-pheno to specify which phenotype/covariate you want to filter on.\n", is_remove? "remove" : "keep");
	    goto keep_remove_cats_ret_INCONSISTENT_INPUT_WW;
	  }
	  cat_pheno_idx = pheno_idx;
	}
      }
      if (cat_pheno_idx != 0xffffffffU) {
        cur_pheno_col = &(pheno_cols[cat_pheno_idx]);
      } else {
	for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx) {
	  if (covar_cols[covar_idx].type_code == kPhenoDtypeCat) {
	    if (cat_pheno_idx != 0xffffffffU) {
	      sprintf(g_logbuf, "Error: Multiple categorical covariates and no categorical phenotype present. Use --%s-cat-pheno to specify which phenotype/covariate you want to filter on.\n", is_remove? "remove" : "keep");
	      goto keep_remove_cats_ret_INCONSISTENT_INPUT_WW;
	    }
	    cat_pheno_idx = covar_idx;
	  }
	}
	if (cat_pheno_idx == 0xffffffffU) {
	  sprintf(g_logbuf, "Error: --%s-cat-names requires a categorical phenotype or covariate.\n", is_remove? "remove-cats/--remove" : "keep-cats/--keep");
	  goto keep_remove_cats_ret_INCONSISTENT_INPUT_WW;
	}
	cur_pheno_col = &(covar_cols[cat_pheno_idx]);
      }
      reterr = keep_remove_cats_internal(cur_pheno_col, cats_fname, cat_names_flattened, raw_sample_ct, is_remove, max_thread_ct, sample_include, sample_ct_ptr);
      if (reterr) {
	goto keep_remove_cats_ret_1;
      }
    } else {
      const uintptr_t name_blen = 1 + strlen(cat_phenoname);
      uint32_t success = 0;
      if (name_blen <= max_pheno_name_blen) {
	for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	  if (!memcmp(cat_phenoname, &(pheno_names[pheno_idx * max_pheno_name_blen]), name_blen)) {
	    const pheno_col_t* cur_pheno_col = &(pheno_cols[pheno_idx]);
	    if (cur_pheno_col->type_code != kPhenoDtypeCat) {
	      sprintf(g_logbuf, "Error: '%s' is not a categorical phenotype.\n", cat_phenoname);
	      goto keep_remove_cats_ret_INCONSISTENT_INPUT_WW;
	    }
	    reterr = keep_remove_cats_internal(cur_pheno_col, cats_fname, cat_names_flattened, raw_sample_ct, is_remove, max_thread_ct, sample_include, sample_ct_ptr);
	    if (reterr) {
	      goto keep_remove_cats_ret_1;
	    }
	    success = 1;
	    break;
	  }
	}
      }
      if (name_blen <= max_covar_name_blen) {
	for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx) {
	  if (!memcmp(cat_phenoname, &(covar_names[covar_idx * max_covar_name_blen]), name_blen)) {
	    const pheno_col_t* cur_pheno_col = &(covar_cols[covar_idx]);
	    if (cur_pheno_col->type_code != kPhenoDtypeCat) {
	      sprintf(g_logbuf, "Error: '%s' is not a categorical covariate.\n", cat_phenoname);
	      goto keep_remove_cats_ret_INCONSISTENT_INPUT_WW;
	    }
	    reterr = keep_remove_cats_internal(cur_pheno_col, cats_fname, cat_names_flattened, raw_sample_ct, is_remove, max_thread_ct, sample_include, sample_ct_ptr);
	    if (reterr) {
	      goto keep_remove_cats_ret_1;
	    }
	    success = 1;
	    break;
	  }
	}
      }
      if (!success) {
	sprintf(g_logbuf, "Error: --%s-cat-pheno phenotype/covariate not loaded.\n", is_remove? "remove" : "keep");
	goto keep_remove_cats_ret_INCONSISTENT_INPUT_2;
      }
    }
  }
  while (0) {
  keep_remove_cats_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
  keep_remove_cats_ret_INCONSISTENT_INPUT_2:
    logerrprintb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 keep_remove_cats_ret_1:
  return reterr;
}

void compute_allele_freqs(const uintptr_t* variant_include, const uintptr_t* variant_allele_idxs, const uint64_t* founder_allele_dosages, uint32_t variant_ct, uint32_t maf_succ, double* allele_freqs) {
  // ok for maj_alleles or allele_freqs to be nullptr
  // note that founder_allele_dosages is in 32768ths
  uint32_t cur_allele_ct = 2;
  uint32_t variant_uidx = 0;
  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    uintptr_t variant_allele_idx_base;
    if (!variant_allele_idxs) {
      variant_allele_idx_base = 2 * variant_uidx;
    } else {
      variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
    }
    const uint64_t* cur_founder_allele_dosages = &(founder_allele_dosages[variant_allele_idx_base]);
    uint64_t tot_dosage = 0;
    for (uint32_t allele_idx = 0; allele_idx < cur_allele_ct; ++allele_idx) {
      tot_dosage += cur_founder_allele_dosages[allele_idx];
    }
    // todo: try changing this expression
    const uint64_t cur_maf_succ_dosage = (maf_succ | (!tot_dosage)) * kDosageMax;

    tot_dosage += cur_maf_succ_dosage * cur_allele_ct;
    double* cur_allele_freqs_base = &(allele_freqs[variant_allele_idx_base - variant_uidx]);
    const double tot_dosage_recip = 1.0 / ((double)((int64_t)tot_dosage));
    const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
    for (uint32_t allele_idx = 0; allele_idx < cur_allele_ct_m1; ++allele_idx) {
      const double cur_dosage = (double)((int64_t)(cur_founder_allele_dosages[allele_idx] + cur_maf_succ_dosage));
      cur_allele_freqs_base[allele_idx] = cur_dosage * tot_dosage_recip;
    }
  }
}

CONSTU31(kMaxReadFreqAlleles, 255);

// relevant column types:
// 0: variant ID
// 1: ref allele code
// 2: all alt allele codes (potentially just alt1)
//
// (3-4 are --freq only)
// 3: ref freq/count
// 4: either all freqs/counts, or all-but-ref
//
// 5: obs ct (only relevant for --freq, but can be in --geno-counts)
//
// (6-11 are --geno-counts/--freqx only)
// 6: hom-ref count
// 7: het ref-alt counts (worst case, just ref-alt1)
// 8: altx-alty counts (worst case, just hom-alt1), or all pairs
// 9: hap-ref count
// 10: hap-alt counts (worst case, just hap-alt1), or all hap counts
// 11: --geno-counts numeq (if present, ignore 6..10)
//
// overrideable:
// 12->2: ALT1
// 13->4: ALT1_FREQ/ALT1_CT
// 14->7: HET_REF_ALT1_CT
// 15->8: HOM_ALT1_CT
// 16->10: HAP_ALT1_CT
ENUM_U31_DEF_START()
  kfReadFreqColVarId = 0,
  kfReadFreqColRefAllele,
  kfReadFreqColAltAlleles,

  kfReadFreqColRefFreq,
  kfReadFreqColAltFreqs,

  kfReadFreqColObsCt,

  kfReadFreqColHomRefCt,
  kfReadFreqColHetRefAltCts,
  kfReadFreqColNonrefDiploidCts,
  kfReadFreqColHapRefCt,
  kfReadFreqColHapAltCts,
  kfReadFreqColGenoCtNumeq,

  kfReadFreqColAlt1Allele,
  kfReadFreqColAlt1Freq,
  kfReadFreqColHetRefAlt1Ct,
  kfReadFreqColHomAlt1Ct,
  kfReadFreqColHapAlt1Ct,

  kfReadFreqColNull
ENUM_U31_DEF_END(read_freq_colidx_t);

FLAGSET_DEF_START()
  kfReadFreqColset0,
  kfReadFreqColsetVarId = (1 << kfReadFreqColVarId),
  kfReadFreqColsetRefAllele = (1 << kfReadFreqColRefAllele),
  kfReadFreqColsetAltAlleles = (1 << kfReadFreqColAltAlleles),
  kfReadFreqColsetBase = (kfReadFreqColsetVarId | kfReadFreqColsetRefAllele | kfReadFreqColsetAltAlleles),

  kfReadFreqColsetRefFreq = (1 << kfReadFreqColRefFreq),
  kfReadFreqColsetAltFreqs = (1 << kfReadFreqColAltFreqs),
  kfReadFreqColsetAfreqOnly = (kfReadFreqColsetRefFreq | kfReadFreqColsetAltFreqs),

  kfReadFreqColsetObsCt = (1 << kfReadFreqColObsCt),

  kfReadFreqColsetHomRefCt = (1 << kfReadFreqColHomRefCt),
  kfReadFreqColsetHetRefAltCts = (1 << kfReadFreqColHetRefAltCts),
  kfReadFreqColsetNonrefDiploidCts = (1 << kfReadFreqColNonrefDiploidCts),
  kfReadFreqColsetHapRefCt = (1 << kfReadFreqColHapRefCt),
  kfReadFreqColsetHapAltCts = (1 << kfReadFreqColHapAltCts),
  kfReadFreqColsetGcountDefault = ((kfReadFreqColsetHapAltCts * 2) - kfReadFreqColsetHomRefCt),

  kfReadFreqColsetGenoCtNumeq = (1 << kfReadFreqColGenoCtNumeq),
  kfReadFreqColsetGcountOnly = (kfReadFreqColsetGcountDefault | kfReadFreqColsetGenoCtNumeq),

  kfReadFreqColsetAlt1Allele = (1 << kfReadFreqColAlt1Allele),
  kfReadFreqColsetAlt1Freq = (1 << kfReadFreqColAlt1Freq),
  kfReadFreqColsetHetRefAlt1Ct = (1 << kfReadFreqColHetRefAlt1Ct),
  kfReadFreqColsetHomAlt1Ct = (1 << kfReadFreqColHomAlt1Ct),
  kfReadFreqColsetHapAlt1Ct = (1 << kfReadFreqColHapAlt1Ct)
FLAGSET_DEF_END(read_freq_colset_t);

pglerr_t read_allele_freqs(const uintptr_t* variant_include, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const char* read_freq_fname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_alt_allele_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t maf_succ, uint32_t max_thread_ct, double* allele_freqs) {
  // support PLINK 1.9 --freq/--freqx, and 2.0 --freq/--geno-counts.
  // GCTA-format no longer supported since it inhibits the allele consistency
  // check.
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_infile = nullptr;
  uintptr_t loadbuf_size = 0;
  uintptr_t line_idx = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    double* cur_allele_freqs;
    uintptr_t* matched_loaded_alleles;
    uintptr_t* matched_internal_alleles;
    uint32_t* loaded_to_internal_allele_idx;
    uintptr_t* already_seen;
    if (bigstack_calloc_d(kMaxReadFreqAlleles, &cur_allele_freqs) ||
	bigstack_alloc_ul(BITCT_TO_WORDCT(kMaxReadFreqAlleles), &matched_loaded_alleles) ||
	bigstack_alloc_ul(BITCT_TO_WORDCT(max_alt_allele_ct + 1), &matched_internal_alleles) ||
	bigstack_alloc_ui(kMaxReadFreqAlleles, &loaded_to_internal_allele_idx) ||
	bigstack_calloc_ul(BITCT_TO_WORDCT(raw_variant_ct), &already_seen)) {
      goto read_allele_freqs_ret_NOMEM;
    }
    reterr = gzopen_read_checked(read_freq_fname, &gz_infile);
    if (reterr) {
      goto read_allele_freqs_ret_1;
    }
    loadbuf_size = bigstack_left() / 8;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto read_allele_freqs_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    uint32_t* variant_id_htable = nullptr;
    uint32_t variant_id_htable_size;
    reterr = alloc_and_populate_id_htable_mt(variant_include, variant_ids, variant_ct, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size);
    if (reterr) {
      goto read_allele_freqs_ret_1;
    }
    char* loadbuf_first_token;
    do {
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto read_allele_freqs_ret_READ_FAIL;
	}
	logerrprint("Error: Empty --read-freq file.\n");
	goto read_allele_freqs_ret_MALFORMED_INPUT;
      }
      ++line_idx;
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_allele_freqs_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      // automatically skip header lines that start with '##' or '# '
    } while (is_eoln_kns(*loadbuf_first_token) || ((*loadbuf_first_token == '#') && ((unsigned char)loadbuf_first_token[1] <= '#')));

    uint32_t col_skips[kfReadFreqColNull];
    read_freq_colidx_t col_types[kfReadFreqColNull];
    uint32_t overrideable_pos[kfReadFreqColNull - kfReadFreqColAlt1Allele];
    uint32_t geno_counts = 0;
    uint32_t main_eq = 0;
    uint32_t is_numeq = 0;
    uint32_t use_obs_ct = 0;
    uint32_t infer_one_freq = 0;
    uint32_t infer_freq_loaded_idx = 0;
    uint32_t relevant_col_ct = 0;

    // interpretation of ColAltAlleles
    uint32_t allele_list_just_alt1 = 1;

    uint32_t is_frac = 0; // if true, one frequency can be missing
    // could add consistency check (can't mix FREQ and CT)

    // interpretation of ColAltFreqs and ColNonrefDiploidCts
    uint32_t main_allele_idx_start = 1;
    uint32_t main_list_just_alt1 = 1;

    // interpretation of ColHapAltCts
    uint32_t hap_allele_idx_start = 1;
    uint32_t hap_list_just_alt1 = 1;

    uint32_t het_list_just_alt1 = 1; // ColHetRefAltCts

    uint32_t biallelic_only = 0;

    read_freq_colset_t header_cols = kfReadFreqColset0;
    if (*loadbuf_first_token == '#') {
      // PLINK 2.0
      char* loadbuf_iter = &(loadbuf_first_token[1]); // guaranteed nonspace
      uint32_t col_idx = 0;
      while (1) {
	char* token_end = token_endnn(loadbuf_iter);
	const uint32_t token_slen = (uintptr_t)(token_end - loadbuf_iter);
	read_freq_colidx_t cur_colidx = kfReadFreqColNull;
	if (token_slen <= 4) {
	  if ((token_slen == 2) && (!memcmp(loadbuf_iter, "ID", 2))) {
	    cur_colidx = kfReadFreqColVarId;
	  } else if (token_slen == 3) {
	    if (!memcmp(loadbuf_iter, "REF", 3)) {
	      cur_colidx = kfReadFreqColRefAllele;
	    } else if (!memcmp(loadbuf_iter, "ALT", 3)) {
	      cur_colidx = kfReadFreqColAltAlleles;
	      if (allele_list_just_alt1) {
		header_cols &= ~kfReadFreqColsetAlt1Allele;
	        allele_list_just_alt1 = 0;
	      }
	    } else if (!memcmp(loadbuf_iter, "CTS", 3)) {
	      goto read_allele_freqs_freqmain_found1;
	    }
	  } else if (token_slen == 4) {
	    if ((!memcmp(loadbuf_iter, "ALT1", 4)) && allele_list_just_alt1) {
	      cur_colidx = kfReadFreqColAlt1Allele;
	    }
	  }
	} else if (((token_slen == 8) && (!memcmp(loadbuf_iter, "REF_FREQ", 8))) || ((token_slen == 6) && (!memcmp(loadbuf_iter, "REF_CT", 6)))) {
	  cur_colidx = kfReadFreqColRefFreq;
	  if (loadbuf_iter[4] == 'F') {
	    is_frac = 1;
	  }
	} else if ((((token_slen == 9) && (!memcmp(loadbuf_iter, "ALT1_FREQ", 9))) || ((token_slen == 7) && (!memcmp(loadbuf_iter, "ALT1_CT", 7)))) && main_list_just_alt1) {
	  cur_colidx = kfReadFreqColAlt1Freq;
	  if (loadbuf_iter[5] == 'F') {
	    is_frac = 1;
	  }
	} else if (((token_slen == 9) && (!memcmp(loadbuf_iter, "ALT_FREQS", 9))) || ((token_slen == 7) && (!memcmp(loadbuf_iter, "ALT_CTS", 7)))) {
	  if (loadbuf_iter[4] == 'F') {
	    is_frac = 1;
	  }
	  goto read_allele_freqs_freqmain_found2;
	} else if ((token_slen == 5) && (!memcmp(loadbuf_iter, "FREQS", 5))) {
	  is_frac = 1;
	  goto read_allele_freqs_freqmain_found1;
	} else if (((token_slen == 13) && (!memcmp(loadbuf_iter, "ALT_NUM_FREQS", 13))) || ((token_slen == 11) && (!memcmp(loadbuf_iter, "ALT_NUM_CTS", 11)))) {
	  is_numeq = 1;
	  goto read_allele_freqs_freqmain_found2;
	} else if (((token_slen == 9) && (!memcmp(loadbuf_iter, "NUM_FREQS", 9))) || ((token_slen == 7) && (!memcmp(loadbuf_iter, "NUM_CTS", 7)))) {
	  is_numeq = 1;
	read_allele_freqs_freqmain_found1:
	  main_allele_idx_start = 0;
	read_allele_freqs_freqmain_found2:
	  cur_colidx = kfReadFreqColAltFreqs;
	  if (main_list_just_alt1) {
	    header_cols &= ~kfReadFreqColsetAlt1Freq;
	    main_list_just_alt1 = 0;
	  }
	} else if ((token_slen == 6) && (!memcmp(loadbuf_iter, "OBS_CT", 6))) {
	  cur_colidx = kfReadFreqColObsCt;
	} else if ((token_slen == 10) && (!memcmp(loadbuf_iter, "HOM_REF_CT", 10))) {
	  cur_colidx = kfReadFreqColHomRefCt;
	} else if ((token_slen == 15) && (!memcmp(loadbuf_iter, "HET_REF_ALT1_CT", 15)) && het_list_just_alt1) {
	  cur_colidx = kfReadFreqColHetRefAlt1Ct;
	} else if ((token_slen == 15) && (!memcmp(loadbuf_iter, "HET_REF_ALT_CTS", 15))) {
	  cur_colidx = kfReadFreqColHetRefAltCts;
	  if (het_list_just_alt1) {
	    header_cols &= ~kfReadFreqColsetHetRefAlt1Ct;
	    het_list_just_alt1 = 0;
	  }
	} else if ((token_slen == 11) && (!memcmp(loadbuf_iter, "HOM_ALT1_CT", 11)) && main_list_just_alt1) {
	  cur_colidx = kfReadFreqColHomAlt1Ct;
	} else if ((token_slen == 23) && (!memcmp(loadbuf_iter, "NONREF_DIPLOID_GENO_CTS", 23))) {
	  goto read_allele_freqs_countmain_found;
	} else if ((token_slen == 16) && (!memcmp(loadbuf_iter, "DIPLOID_GENO_CTS", 16))) {
	  main_allele_idx_start = 0;
	read_allele_freqs_countmain_found:
	  cur_colidx = kfReadFreqColNonrefDiploidCts;
	  if (main_list_just_alt1) {
	    header_cols &= ~kfReadFreqColsetHomAlt1Ct;
	    // could make this use a different variable than FREQS does
	    main_list_just_alt1 = 0;
	  }
	} else if ((token_slen == 10) && (!memcmp(loadbuf_iter, "HAP_REF_CT", 10))) {
	  cur_colidx = kfReadFreqColHapRefCt;
	} else if ((token_slen == 11) && (!memcmp(loadbuf_iter, "HAP_ALT1_CT", 11)) && hap_list_just_alt1) {
	  cur_colidx = kfReadFreqColHapAlt1Ct;
	} else if ((token_slen == 11) && (!memcmp(loadbuf_iter, "HAP_ALT_CTS", 11))) {
	  goto read_allele_freqs_hapmain_found;
	} else if ((token_slen == 7) && (!memcmp(loadbuf_iter, "HAP_CTS", 7))) {
	  hap_allele_idx_start = 0;
	read_allele_freqs_hapmain_found:
	  cur_colidx = kfReadFreqColHapAltCts;
	  if (hap_list_just_alt1) {
	    header_cols &= ~kfReadFreqColsetHapAlt1Ct;
	    hap_list_just_alt1 = 0;
	  }
	} else if ((token_slen == 12) && (!memcmp(loadbuf_iter, "GENO_NUM_CTS", 12))) {
	  cur_colidx = kfReadFreqColGenoCtNumeq;
	  is_numeq = 1;
	}
	if (cur_colidx != kfReadFreqColNull) {
	  const read_freq_colset_t cur_colset = (read_freq_colset_t)(1 << cur_colidx);
	  if (header_cols & cur_colset) {
	    logerrprint("Error: Conflicting columns in header line of --read-freq file.\n");
	    goto read_allele_freqs_ret_MALFORMED_INPUT;
	  }
	  if (cur_colidx >= kfReadFreqColAlt1Allele) {
	    overrideable_pos[cur_colidx - kfReadFreqColAlt1Allele] = relevant_col_ct;
	  }
	  header_cols |= cur_colset;
	  col_skips[relevant_col_ct] = col_idx;
	  col_types[relevant_col_ct++] = cur_colidx;
	}
	loadbuf_iter = skip_initial_spaces(token_end);
	if (is_eoln_kns(*loadbuf_iter)) {
	  break;
	}
	++col_idx;
      }
      read_freq_colset_t semifinal_header_cols = header_cols;
      if (header_cols & kfReadFreqColsetAlt1Allele) {
	header_cols ^= kfReadFreqColsetAltAlleles | kfReadFreqColsetAlt1Allele;
	col_types[overrideable_pos[0]] = kfReadFreqColAltAlleles;
      }
      if (header_cols & kfReadFreqColsetAlt1Freq) {
	header_cols ^= kfReadFreqColsetAltFreqs | kfReadFreqColsetAlt1Freq;
	col_types[overrideable_pos[kfReadFreqColAlt1Freq - kfReadFreqColAlt1Allele]] = kfReadFreqColAltFreqs;
      }
      if (header_cols & kfReadFreqColsetHetRefAlt1Ct) {
	header_cols ^= kfReadFreqColsetHetRefAltCts | kfReadFreqColsetHetRefAlt1Ct;
	col_types[overrideable_pos[kfReadFreqColHetRefAlt1Ct - kfReadFreqColAlt1Allele]] = kfReadFreqColHetRefAltCts;
      }
      if (header_cols & kfReadFreqColsetHomAlt1Ct) {
	header_cols ^= kfReadFreqColsetNonrefDiploidCts | kfReadFreqColsetHomAlt1Ct;
	col_types[overrideable_pos[kfReadFreqColHomAlt1Ct - kfReadFreqColAlt1Allele]] = kfReadFreqColNonrefDiploidCts;
      }
      if (header_cols & kfReadFreqColsetHapAlt1Ct) {
	header_cols ^= kfReadFreqColsetHapAltCts | kfReadFreqColsetHapAlt1Ct;
	col_types[overrideable_pos[kfReadFreqColHapAlt1Ct - kfReadFreqColAlt1Allele]] = kfReadFreqColHapAltCts;
      }
      if ((semifinal_header_cols != header_cols) && (!(header_cols & kfReadFreqColsetGenoCtNumeq))) {
	// we're treating at least one ALT1 column as if it spoke for all ALT
	// alleles
	biallelic_only = 1;
      }

      main_eq = is_numeq;
      semifinal_header_cols = header_cols;
      if (header_cols & kfReadFreqColsetAfreqOnly) {
	if (header_cols & kfReadFreqColsetGcountOnly) {
	  logerrprint("Error: Conflicting columns in header line of --read-freq file (--freq and\n--geno-counts values mixed together).\n");
	  goto read_allele_freqs_ret_MALFORMED_INPUT;
	}
	read_freq_colset_t header_cols_exempt = kfReadFreqColset0;
	if ((header_cols & kfReadFreqColsetAltFreqs) && (!is_numeq)) {
	  // {ALT_}FREQS can be formatted as either
	  //   0.5,0,0.2
	  // or
	  //   A=0.5,G=0.2
	  // Look at the first nonheader line to distinguish between these two.
	  do {
	    if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	      if (!gzeof(gz_infile)) {
		goto read_allele_freqs_ret_READ_FAIL;
	      }
	      logerrprint("Error: Empty --read-freq file.\n");
	      goto read_allele_freqs_ret_MALFORMED_INPUT;
	    }
	    ++line_idx;
	    if (!loadbuf[loadbuf_size - 1]) {
	      goto read_allele_freqs_ret_LONG_LINE;
	    }
	    loadbuf_first_token = skip_initial_spaces(loadbuf);
	  } while (is_eoln_kns(*loadbuf_first_token));
	  char* loadbuf_iter = loadbuf_first_token;
	  char* alt_freq_str = nullptr;
	  for (uint32_t relevant_col_idx = 0; relevant_col_idx < relevant_col_ct; ++relevant_col_idx) {
	    if (col_types[relevant_col_idx] == kfReadFreqColAltFreqs) {
	      alt_freq_str = next_token_multz(loadbuf_iter, col_skips[relevant_col_idx]);
	      break;
	    }
	  }
	  if (!alt_freq_str) {
	    goto read_allele_freqs_ret_MISSING_TOKENS;
	  }
	  const uint32_t alt_freq_slen = token_endnn(alt_freq_str) - alt_freq_str;
	  // bare '.' can only appear in eq formats
	  main_eq = ((alt_freq_slen == 1) && (*alt_freq_str == '.')) || (memchr(alt_freq_str, '=', alt_freq_slen) != nullptr);
	  if (main_eq) {
	    header_cols_exempt = kfReadFreqColsetAltAlleles;
	    if (!main_allele_idx_start) {
	      header_cols_exempt |= kfReadFreqColsetRefAllele;
	    }
	    header_cols &= ~header_cols_exempt;
	  }
	}
	if (((header_cols & kfReadFreqColsetBase) | header_cols_exempt) != kfReadFreqColsetBase) {
	  logerrprint("Error: Missing column(s) in --read-freq file (ID, REF, ALT{1} usually\nrequired).\n");
	  goto read_allele_freqs_ret_MALFORMED_INPUT;
	}
	if (!main_allele_idx_start) {
	  header_cols &= ~kfReadFreqColsetRefFreq;
	} else {
	  if ((header_cols & (kfReadFreqColsetRefFreq | kfReadFreqColsetAltFreqs)) != (kfReadFreqColsetRefFreq | kfReadFreqColsetAltFreqs)) {
	    if (main_list_just_alt1) {
	      biallelic_only = 1;
	    }
	    infer_one_freq = 1;
	    infer_freq_loaded_idx = (header_cols / kfReadFreqColsetRefFreq) & 1;
	    if (!is_frac) {
	      if (!(header_cols & kfReadFreqColsetObsCt)) {
		logerrprint("Error: Missing column(s) in --read-freq file (at least two of {REF_CT, ALT1_CT,\nALT_CTS, OBS_CT} must be present).\n");
		goto read_allele_freqs_ret_MALFORMED_INPUT;
	      }
	      use_obs_ct = 1;
	    }
	  }
	}
	logprint("--read-freq: PLINK 2 --freq file detected.\n");
      } else if (header_cols & kfReadFreqColsetGcountOnly) {
	if ((header_cols & kfReadFreqColsetBase) != kfReadFreqColsetBase) {
	  logerrprint("Error: Missing column(s) in --read-freq file (ID, REF, ALT{1} required).\n");
	  goto read_allele_freqs_ret_MALFORMED_INPUT;
	}
	// possible todo: allow one frequency/count to be missing.  (not really
	// necessary since PLINK 1.9 --freqx does not leave anything out,
	// unlike PLINK 1.x --freq)
	if (header_cols & kfReadFreqColsetGenoCtNumeq) {
	  header_cols &= ~kfReadFreqColsetGcountDefault; // don't need anything but GENO_NUM_CTS
	} else {
	  // require both diploid and haploid columns for now.  (could
	  // conditionally drop one of these requirements later.)
	  if (!(header_cols & kfReadFreqColsetNonrefDiploidCts)) {
	    logerrprint("Error: Missing column(s) in --read-freq file (HOM_ALT1_CT,\nNONREF_DIPLOID_GENO_CTS, or DIPLOID_GENO_CTS required).\n");
	    goto read_allele_freqs_ret_MALFORMED_INPUT;
	  }
	  if (!main_allele_idx_start) {
	    header_cols &= ~(kfReadFreqColsetHomRefCt | kfReadFreqColsetHetRefAltCts);
	  } else if ((header_cols & (kfReadFreqColsetHomRefCt | kfReadFreqColsetHetRefAltCts)) != (kfReadFreqColsetHomRefCt | kfReadFreqColsetHetRefAltCts)) {
	    logerrprint("Error: Missing column(s) in --read-freq file (HOM_REF_CT, HET_REF_ALT1_CT, or\nHET_REF_ALT_CTS required unless {DIPLOID_}GENO_CTS present).\n");
	    goto read_allele_freqs_ret_MALFORMED_INPUT;
	  }
	  if (!(header_cols & kfReadFreqColsetHapAltCts)) {
	    logerrprint("Error: Missing column(s) in --read-freq file (HAP_ALT1_CT, HAP_ALT_CTS, or\nHAP_CTS required).\n");
	    goto read_allele_freqs_ret_MALFORMED_INPUT;
	  }
	  if (!hap_allele_idx_start) {
	    header_cols &= ~kfReadFreqColsetHapRefCt;
	  } else if (!(header_cols & kfReadFreqColsetHapRefCt)) {
	    logerrprint("Error: Missing column(s) in --read-freq file (HAP_REF_CT required unless\nHAP_CTS or GENO_CTS present).\n");
	    goto read_allele_freqs_ret_MALFORMED_INPUT;
	  }
	}
	geno_counts = 1;
	logprint("--read-freq: PLINK 2 --geno-counts file detected.\n");
	*loadbuf_first_token = '\0';
      } else {
	logerrprint("Error: Missing column(s) in --read-freq file (no frequencies/counts).\n");
	goto read_allele_freqs_ret_MALFORMED_INPUT;
      }
      if (!use_obs_ct) {
	header_cols &= ~kfReadFreqColsetObsCt;
      }
      if (semifinal_header_cols != header_cols) {
	// remove redundant columns
	uint32_t relevant_col_idx_read = 0;
	while ((((uint32_t)header_cols) >> col_types[relevant_col_idx_read]) & 1) {
	  ++relevant_col_idx_read;
	}
	uint32_t relevant_col_idx_write = relevant_col_idx_read;
	for (; relevant_col_idx_read < relevant_col_ct; ++relevant_col_idx_read) {
	  const read_freq_colidx_t cur_colidx = col_types[relevant_col_idx_read];
	  if ((((uint32_t)header_cols) >> cur_colidx) & 1) {
	    col_types[relevant_col_idx_write] = cur_colidx;
	    col_skips[relevant_col_idx_write] = col_skips[relevant_col_idx_read];
	    ++relevant_col_idx_write;
	  }
	}
	relevant_col_ct = relevant_col_idx_write;
      }
      for (uint32_t uii = relevant_col_ct - 1; uii; --uii) {
	col_skips[uii] -= col_skips[uii - 1];
      }
    } else {
      // PLINK 1.x
      // .frq:       CHR  SNP  A1  A2  MAF        NCHROBS
      // .frq.count: CHR  SNP  A1  A2  C1         C2       G0
      // .frqx:      CHR  SNP  A1  A2  C(HOM A1)  C(HET)   C(HOM A2)  C(HAP A1)
      //   C(HAP A2)  C(MISSING)
      // (yeah, the spaces in the .frqx header were a mistake, should have used
      // underscores.  oh well, live and learn.)
      col_skips[0] = 1;
      col_skips[1] = 1;
      col_skips[2] = 1;
      col_skips[3] = 1;

      col_types[0] = kfReadFreqColVarId;
      // doesn't matter if we treat A1 or A2 as ref
      col_types[1] = kfReadFreqColRefAllele;
      col_types[2] = kfReadFreqColAltAlleles;
      biallelic_only = 1;
      if (!memcmp(loadbuf_first_token, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)", 71)) {
        col_skips[4] = 1;
	col_skips[5] = 1;
	col_skips[6] = 1;
	col_skips[7] = 1;

	col_types[3] = kfReadFreqColHomRefCt;
	col_types[4] = kfReadFreqColHetRefAltCts;
	col_types[5] = kfReadFreqColNonrefDiploidCts;
	col_types[6] = kfReadFreqColHapRefCt;
	col_types[7] = kfReadFreqColHapAltCts;
	header_cols = kfReadFreqColsetBase | kfReadFreqColsetGcountOnly;
	geno_counts = 1;
	relevant_col_ct = 8;
	logprint("--read-freq: PLINK 1.9 --freqx file detected.\n");
      } else {
	if (strcmp_se(loadbuf_first_token, "CHR", 3)) {
	  goto read_allele_freqs_ret_UNRECOGNIZED_HEADER;
	}
	char* loadbuf_iter = skip_initial_spaces(&(loadbuf_first_token[3]));
	if (strcmp_se(loadbuf_iter, "SNP", 3)) {
	  goto read_allele_freqs_ret_UNRECOGNIZED_HEADER;
	}
	loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[3]));
	if (strcmp_se(loadbuf_iter, "A1", 2)) {
	  goto read_allele_freqs_ret_UNRECOGNIZED_HEADER;
	}
	loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[2]));
	if (strcmp_se(loadbuf_iter, "A2", 2)) {
	  goto read_allele_freqs_ret_UNRECOGNIZED_HEADER;
	}
	loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[2]));
	col_types[3] = kfReadFreqColRefFreq;
	if (!strcmp_se(loadbuf_iter, "MAF", 3)) {
	  is_frac = 1;
	  infer_one_freq = 1;
	  infer_freq_loaded_idx = 1;
	  header_cols = kfReadFreqColsetBase | kfReadFreqColsetRefFreq;
	  relevant_col_ct = 4;
	  logprint("--read-freq: PLINK 1.x --freq file detected.\n");
	} else {
	  if (strcmp_se(loadbuf_iter, "C1", 2)) {
	    goto read_allele_freqs_ret_UNRECOGNIZED_HEADER;
	  }
	  loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[2]));
	  if (strcmp_se(loadbuf_iter, "C2", 2)) {
	    goto read_allele_freqs_ret_UNRECOGNIZED_HEADER;
	  }
	  col_skips[4] = 1;
	  col_types[4] = kfReadFreqColAltFreqs;
	  header_cols = kfReadFreqColsetBase | kfReadFreqColsetAfreqOnly;
	  relevant_col_ct = 5;
	  logprint("--read-freq: PLINK 1.x '--freq counts' file detected.\n");
	}
      }
      *loadbuf_first_token = '\0';
    }
    assert(relevant_col_ct <= 8);

    double freq_max = 4294967295.0;
    if (is_frac) {
      maf_succ = 0;
      freq_max = 1.0;
    }
    uintptr_t skipped_variant_ct = 0;
    uint32_t loaded_variant_ct = 0;
    uint32_t cur_allele_ct = 2;
    while (1) {
      if (!is_eoln_kns(*loadbuf_first_token)) {
	char* loadbuf_iter = loadbuf_first_token;
	char* token_ptrs[12];
	uint32_t token_slens[12];
	for (uint32_t relevant_col_idx = 0; relevant_col_idx < relevant_col_ct; ++relevant_col_idx) {
	  const read_freq_colidx_t cur_colidx = col_types[relevant_col_idx];
	  loadbuf_iter = next_token_multz(loadbuf_iter, col_skips[relevant_col_idx]);
	  if (!loadbuf_iter) {
	    goto read_allele_freqs_ret_MISSING_TOKENS;
	  }
	  token_ptrs[cur_colidx] = loadbuf_iter;
	  char* token_end = token_endnn(loadbuf_iter);
	  token_slens[cur_colidx] = (uintptr_t)(token_end - loadbuf_iter);
	  loadbuf_iter = token_end;
	}
	char* variant_id_start = token_ptrs[kfReadFreqColVarId];
	const uint32_t variant_id_slen = token_slens[kfReadFreqColVarId];
	uint32_t variant_uidx = variant_id_dupflag_htable_find(variant_id_start, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
	if (variant_uidx >> 31) {
	  if (variant_uidx == 0xffffffffU) {
	    goto read_allele_freqs_skip_variant;
	  }
	  sprintf(g_logbuf, "Error: --read-freq variant ID '%s' appears multiple times in main dataset.\n", variant_ids[variant_uidx & 0xffffffffU]);
	  goto read_allele_freqs_ret_MALFORMED_INPUT_WW;
	}
	if (is_set(already_seen, variant_uidx)) {
	  sprintf(g_logbuf, "Error: Variant ID '%s' appears multiple times in --read-freq file.\n", variant_ids[variant_uidx]);
	  goto read_allele_freqs_ret_MALFORMED_INPUT_WW;
	}
	set_bit(variant_uidx, already_seen);

	uintptr_t variant_allele_idx_base;
	if (!variant_allele_idxs) {
	  variant_allele_idx_base = variant_uidx * 2;
	} else {
	  variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	  cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
	  if (biallelic_only) {
	    goto read_allele_freqs_skip_variant;
	  }
	}
	fill_ulong_zero(BITCT_TO_WORDCT(cur_allele_ct), matched_internal_alleles);
	char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
	uint32_t loaded_allele_ct = 0;
	if (header_cols & kfReadFreqColsetRefAllele) {
	  uint32_t cur_loaded_allele_code_slen = token_slens[kfReadFreqColRefAllele];
	  uint32_t unmatched_allele_ct = cur_allele_ct;
	  char* cur_loaded_allele_code = token_ptrs[kfReadFreqColRefAllele];
	  cur_loaded_allele_code[cur_loaded_allele_code_slen] = '\0';
	  char* loaded_allele_code_iter;
	  char* loaded_allele_code_end;
	  if (header_cols & kfReadFreqColsetAltAlleles) {
	    loaded_allele_code_iter = token_ptrs[kfReadFreqColAltAlleles];
	    loaded_allele_code_end = &(loaded_allele_code_iter[token_slens[kfReadFreqColAltAlleles]]);
	    *loaded_allele_code_end++ = ',';
	  } else {
	    // special case: with --freq alteq or alteqz column, we only need
	    // to scrape REF here
	    loaded_allele_code_iter = &(cur_loaded_allele_code[cur_loaded_allele_code_slen + 1]);
	    loaded_allele_code_end = loaded_allele_code_iter;
	  }
	  uint32_t widx = 0;
	  while (1) {
	    if (!(loaded_allele_ct % kBitsPerWord)) {
	      widx = loaded_allele_ct / kBitsPerWord;
	      matched_loaded_alleles[widx] = 0;
	    }
	    if (cur_loaded_allele_code_slen <= max_allele_slen) {
	      uint32_t internal_allele_idx = 0;
	      uint32_t unmatched_allele_idx = 0;
	      for (; unmatched_allele_idx < unmatched_allele_ct; ++unmatched_allele_idx, ++internal_allele_idx) {
		next_unset_unsafe_ck(matched_internal_alleles, &internal_allele_idx);
		if (!strcmp(cur_loaded_allele_code, cur_alleles[internal_allele_idx])) {
		  break;
		}
	      }
	      if (unmatched_allele_idx != unmatched_allele_ct) {
		// success
		if (is_set(matched_internal_alleles, internal_allele_idx)) {
		  sprintf(g_logbuf, "Error: Duplicate allele code on line %" PRIuPTR " of --read-freq file.\n", line_idx);
		  goto read_allele_freqs_ret_MALFORMED_INPUT_2;
		}
		set_bit(internal_allele_idx, matched_internal_alleles);
		set_bit(loaded_allele_ct, matched_loaded_alleles);
		loaded_to_internal_allele_idx[loaded_allele_ct] = internal_allele_idx;
	      }
	    }
	    ++loaded_allele_ct;
	    if (loaded_allele_code_iter == loaded_allele_code_end) {
	      break;
	    }
	    if (loaded_allele_ct == kMaxReadFreqAlleles) {
	      sprintf(g_logbuf, "Error: --read-freq file entry for variant ID '%s' has more than %u ALT alleles.\n", variant_ids[variant_uidx], kMaxReadFreqAlleles - 1);
	      goto read_allele_freqs_ret_MALFORMED_INPUT_WW;
	    }
	    cur_loaded_allele_code = loaded_allele_code_iter;
	    loaded_allele_code_iter = (char*)rawmemchr(loaded_allele_code_iter, ',');
	    cur_loaded_allele_code_slen = (uintptr_t)(loaded_allele_code_iter - cur_loaded_allele_code);
	    *loaded_allele_code_iter++ = '\0';
	  }
	}

	double* allele_freqs_write = &(allele_freqs[variant_allele_idx_base - variant_uidx]);
	if (geno_counts) {
	  fill_double_zero(cur_allele_ct, cur_allele_freqs);
	  if (is_numeq) {
	    char* geno_num_cts_iter = token_ptrs[kfReadFreqColGenoCtNumeq];
	    const uint32_t full_slen = token_slens[kfReadFreqColGenoCtNumeq];
	    char* geno_num_cts_end = &(geno_num_cts_iter[full_slen]);
	    if (full_slen > 1) {
	      *geno_num_cts_end = ',';
#ifndef __LP64__
	      const uint32_t cap_div_10 = (loaded_allele_ct - 1) / 10;
	      const uint32_t cap_mod_10 = (loaded_allele_ct - 1) % 10;
#endif
	      while (1) {
		uint32_t second_loaded_allele_idx = 0xffffffffU;
		uint32_t first_loaded_allele_idx;
#ifdef __LP64__
		if (scanadv_uint_capped(loaded_allele_ct - 1, &geno_num_cts_iter, &first_loaded_allele_idx)) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		if (*geno_num_cts_iter == '/') {
		  ++geno_num_cts_iter;
		  if (scanadv_uint_capped(loaded_allele_ct - 1, &geno_num_cts_iter, &second_loaded_allele_idx)) {
		    goto read_allele_freqs_ret_INVALID_FREQS;
		  }
		}
#else
		if (scanadv_uint_capped32(cap_div_10, cap_mod_10, &geno_num_cts_iter, &first_loaded_allele_idx)) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		if (*geno_num_cts_iter == '/') {
		  ++geno_num_cts_iter;
		  if (scanadv_uint_capped32(cap_div_10, cap_mod_10, &geno_num_cts_iter, &second_loaded_allele_idx)) {
		    goto read_allele_freqs_ret_INVALID_FREQS;
		  }
		}
#endif
		if (*geno_num_cts_iter != '=') {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		++geno_num_cts_iter;
		double dxx;
		char* cur_ct_end = scanadv_double(geno_num_cts_iter, &dxx);
		if ((!cur_ct_end) || (*cur_ct_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0)) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		if (is_set(matched_loaded_alleles, first_loaded_allele_idx)) {
		  cur_allele_freqs[loaded_to_internal_allele_idx[first_loaded_allele_idx]] += dxx;
		}
		if ((second_loaded_allele_idx != 0xffffffffU) && is_set(matched_loaded_alleles, second_loaded_allele_idx)) {
		  cur_allele_freqs[loaded_to_internal_allele_idx[second_loaded_allele_idx]] += dxx;
		}
		geno_num_cts_iter = cur_ct_end;
	        if (geno_num_cts_iter == geno_num_cts_end) {
		  break;
		}
		++geno_num_cts_iter;
	      }
	    } else if (*geno_num_cts_iter != '.') {
	      goto read_allele_freqs_ret_INVALID_FREQS;
	    }
	  } else {
	    const uint32_t internal0 = is_set(matched_loaded_alleles, 0)? loaded_to_internal_allele_idx[0] : 0xffffffffU;
	    if (header_cols & kfReadFreqColsetHomRefCt) {
	      if (internal0 != 0xffffffffU) {
		char* hom_ref_str = token_ptrs[kfReadFreqColHomRefCt];
		double dxx;
		char* hom_ref_end = scanadv_double(hom_ref_str, &dxx);
		if ((!hom_ref_end) || (hom_ref_end != &(hom_ref_str[token_slens[kfReadFreqColHomRefCt]])) || (dxx < 0.0) || (dxx > 4294967295.0)) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		cur_allele_freqs[internal0] += 2 * dxx;
	      }
	      
	      char* het_refalt_iter = token_ptrs[kfReadFreqColHetRefAltCts];
	      char* het_refalt_end = &(het_refalt_iter[token_slens[kfReadFreqColHetRefAltCts]]);
	      *het_refalt_end = ',';
	      for (uint32_t alt_allele_idx = 1; alt_allele_idx < cur_allele_ct; ++alt_allele_idx) {
		if (het_refalt_iter >= het_refalt_end) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		double dxx;
		char* cur_entry_end = scanadv_double(het_refalt_iter, &dxx);
		if ((!cur_entry_end) || (*cur_entry_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0)) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		if (internal0 != 0xffffffffU) {
		  cur_allele_freqs[internal0] += dxx;
		}
		if (is_set(matched_loaded_alleles, alt_allele_idx)) {
		  cur_allele_freqs[loaded_to_internal_allele_idx[alt_allele_idx]] += dxx;
		}
		het_refalt_iter = &(cur_entry_end[1]);
	      }
	    }
	    // ColNonrefDiploidCts required
	    char* diploid_cts_iter = token_ptrs[kfReadFreqColNonrefDiploidCts];
	    char* diploid_cts_end = &(diploid_cts_iter[token_slens[kfReadFreqColNonrefDiploidCts]]);
	    *diploid_cts_end = ',';
	    for (uint32_t second_allele_idx = main_allele_idx_start; second_allele_idx < cur_allele_ct; ++second_allele_idx) {
	      uint32_t internalx = 0xffffffffU;
	      if (is_set(matched_loaded_alleles, second_allele_idx)) {
		internalx = loaded_to_internal_allele_idx[second_allele_idx];
	      }
	      // 1/1, 1/2, 2/2, 1/3, ...
	      for (uint32_t first_allele_idx = main_allele_idx_start; first_allele_idx <= second_allele_idx; ++first_allele_idx) {
		if (diploid_cts_iter >= diploid_cts_end) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
	        double dxx;
		char* cur_entry_end = scanadv_double(diploid_cts_iter, &dxx);
		if ((!cur_entry_end) || (*cur_entry_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0)) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		if (is_set(matched_loaded_alleles, first_allele_idx)) {
		  cur_allele_freqs[loaded_to_internal_allele_idx[first_allele_idx]] += dxx;
		}
		if (internalx != 0xffffffffU) {
		  cur_allele_freqs[internalx] += dxx;
		}
		diploid_cts_iter = &(cur_entry_end[1]);
	      }
	    }

	    if ((header_cols & kfReadFreqColsetHapRefCt) && (internal0 != 0xffffffffU)) {
	      char* hap_ref_str = token_ptrs[kfReadFreqColHapRefCt];
	      double dxx;
	      char* hap_ref_end = scanadv_double(hap_ref_str, &dxx);
	      if ((!hap_ref_end) || (hap_ref_end != &(hap_ref_str[token_slens[kfReadFreqColHapRefCt]])) || (dxx < 0.0) || (dxx > 4294967295.0)) {
		goto read_allele_freqs_ret_INVALID_FREQS;
	      }
	      cur_allele_freqs[internal0] += dxx;
	    }
	    // ColHapAltCts required
	    char* hap_alt_iter = token_ptrs[kfReadFreqColHapAltCts];
	    char* hap_alt_end = &(hap_alt_iter[token_slens[kfReadFreqColHapAltCts]]);
	    *hap_alt_end = ',';
	    for (uint32_t alt_allele_idx = 1; alt_allele_idx < cur_allele_ct; ++alt_allele_idx) {
	      if (hap_alt_iter >= hap_alt_end) {
		goto read_allele_freqs_ret_INVALID_FREQS;
	      }
	      double dxx;
	      char* cur_entry_end = scanadv_double(hap_alt_iter, &dxx);
	      if ((!cur_entry_end) || (*cur_entry_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0)) {
		goto read_allele_freqs_ret_INVALID_FREQS;
	      }
	      if (is_set(matched_loaded_alleles, alt_allele_idx)) {
		cur_allele_freqs[loaded_to_internal_allele_idx[alt_allele_idx]] += dxx;
	      }
	      hap_alt_iter = &(cur_entry_end[1]);
	    }
	  }
	} else {
	  if ((header_cols & kfReadFreqColsetRefFreq) && is_set(matched_loaded_alleles, 0)) {
	    char* ref_freq_str = token_ptrs[kfReadFreqColRefFreq];
	    double dxx;
	    if (!scanadv_double(ref_freq_str, &dxx)) {
	      if (is_nan_str(ref_freq_str, token_slens[kfReadFreqColRefFreq])) {
	        goto read_allele_freqs_skip_variant;
	      }
	      sprintf(g_logbuf, "Error: Invalid REF frequency/count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	      goto read_allele_freqs_ret_MALFORMED_INPUT_WW;
	    }
	    if ((dxx < 0.0) || (dxx > freq_max)) {
	      sprintf(g_logbuf, "Error: Invalid REF frequency/count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	      goto read_allele_freqs_ret_MALFORMED_INPUT_WW;
	    }
	    cur_allele_freqs[loaded_to_internal_allele_idx[0]] = dxx;
	  }
	  if (header_cols & kfReadFreqColsetAltFreqs) {
	    char* alt_freq_iter = token_ptrs[kfReadFreqColAltFreqs];
	    const uint32_t full_slen = token_slens[kfReadFreqColAltFreqs];
	    char* alt_freq_end = &(alt_freq_iter[full_slen]);
	    *alt_freq_end = ',';
	    if (!main_eq) {
	      for (uint32_t allele_idx = main_allele_idx_start; allele_idx < loaded_allele_ct; ++allele_idx, ++alt_freq_iter) {
		if (alt_freq_iter >= alt_freq_end) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		if (!is_set(matched_loaded_alleles, allele_idx)) {
		  alt_freq_iter = (char*)rawmemchr(alt_freq_iter, ',');
		  continue;
		}
		double dxx;
		char* cur_freq_end = scanadv_double(alt_freq_iter, &dxx);
		if (!cur_freq_end) {
		  cur_freq_end = (char*)rawmemchr(alt_freq_iter, ',');
		  if (is_nan_str(alt_freq_iter, (uintptr_t)(cur_freq_end - alt_freq_iter))) {
		    goto read_allele_freqs_skip_variant;
		  }
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		if ((*cur_freq_end != ',') || (dxx < 0.0) || (dxx > freq_max)) {
		  goto read_allele_freqs_ret_INVALID_FREQS;
		}
		alt_freq_iter = cur_freq_end;
		cur_allele_freqs[loaded_to_internal_allele_idx[allele_idx]] = dxx;
	      }
	    } else {
	      fill_double_zero(cur_allele_ct, cur_allele_freqs);
	      if ((full_slen > 1) || (*alt_freq_iter != '.')) {
		if (is_numeq) {
#ifndef __LP64__
		  const uint32_t cap_div_10 = (loaded_allele_ct - 1) / 10;
		  const uint32_t cap_mod_10 = (loaded_allele_ct - 1) % 10;
#endif
		  while (1) {
		    char* cur_entry_end = (char*)rawmemchr(alt_freq_iter, ',');
		    uint32_t loaded_allele_idx;
#ifdef __LP64__
		    if (scanadv_uint_capped(loaded_allele_ct - 1, &alt_freq_iter, &loaded_allele_idx)) {
		      goto read_allele_freqs_ret_INVALID_FREQS;
		    }
#else
		    if (scanadv_uint_capped32(cap_div_10, cap_mod_10, &alt_freq_iter, &loaded_allele_idx)) {
		      goto read_allele_freqs_ret_INVALID_FREQS;
		    }
#endif
		    if (*alt_freq_iter != '=') {
		      goto read_allele_freqs_ret_INVALID_FREQS;
		    }
		    if (is_set(matched_loaded_alleles, loaded_allele_idx)) {
		      const uint32_t internal_allele_idx = loaded_to_internal_allele_idx[loaded_allele_idx];
		      if (cur_allele_freqs[internal_allele_idx]) {
			sprintf(g_logbuf, "Error: Duplicate entry on line %" PRIuPTR " of --read-freq file.\n", line_idx);
			goto read_allele_freqs_ret_MALFORMED_INPUT_2;
		      }
		      ++alt_freq_iter;
		      double dxx;
		      char* cur_freq_end = scanadv_double(alt_freq_iter, &dxx);
		      if (!cur_freq_end) {
			if (is_nan_str(alt_freq_iter, (uintptr_t)(cur_entry_end - alt_freq_iter))) {
			  goto read_allele_freqs_skip_variant;
			}
			goto read_allele_freqs_ret_INVALID_FREQS;
		      }
		      if ((cur_freq_end != cur_entry_end) || (dxx < 0.0) || (dxx > freq_max)) {
			goto read_allele_freqs_ret_INVALID_FREQS;
		      }
		      cur_allele_freqs[internal_allele_idx] = dxx;
		    }
		    alt_freq_iter = cur_entry_end;
		    if (alt_freq_iter == alt_freq_end) {
		      break;
		    }
		    ++alt_freq_iter;
		  }
		} else {
		  while (1) {
		    char* cur_entry_end = (char*)rawmemchr(alt_freq_iter, ',');
		    const uint32_t cur_entry_slen = (uintptr_t)(cur_entry_end - alt_freq_iter);
		    char* eq_ptr = (char*)memchr(alt_freq_iter, '=', cur_entry_slen);
		    if (!eq_ptr) {
		      goto read_allele_freqs_ret_INVALID_FREQS;
		    }
		    *eq_ptr = '\0';
		    uint32_t internal_allele_idx = 0;
		    // O(n^2), may want to replace with O(n log n)
		    for (; internal_allele_idx < cur_allele_ct; ++internal_allele_idx) {
		      if (!strcmp(alt_freq_iter, cur_alleles[internal_allele_idx])) {
			if (cur_allele_freqs[internal_allele_idx]) {
			  sprintf(g_logbuf, "Error: Duplicate entry on line %" PRIuPTR " of --read-freq file.\n", line_idx);
			  goto read_allele_freqs_ret_MALFORMED_INPUT_2;
			}
			alt_freq_iter = &(eq_ptr[1]);
			double dxx;
			char* cur_freq_end = scanadv_double(alt_freq_iter, &dxx);
			if (!cur_freq_end) {
			  if (is_nan_str(alt_freq_iter, (uintptr_t)(cur_entry_end - alt_freq_iter))) {
			    goto read_allele_freqs_skip_variant;
			  }
			  goto read_allele_freqs_ret_INVALID_FREQS;
			}
			if ((cur_freq_end != cur_entry_end) || (dxx < 0.0) || (dxx > freq_max)) {
			  goto read_allele_freqs_ret_INVALID_FREQS;
			}
			cur_allele_freqs[internal_allele_idx] = dxx;
			break;
		      }
		    }
		    alt_freq_iter = cur_entry_end;
		    if (alt_freq_iter == alt_freq_end) {
		      break;
		    }
		    ++alt_freq_iter;
		  }
		}
	      }
	    }
	  }
	}
	if (infer_one_freq && is_set(matched_loaded_alleles, infer_freq_loaded_idx)) {
	  double obs_ct_recip = 1.0;
	  if (header_cols & kfReadFreqColsetObsCt) {
	    uint32_t obs_ct_raw;
	    if (scan_uint_capped(token_ptrs[kfReadFreqColObsCt], 0xffffffffU, &obs_ct_raw)) {
	      sprintf(g_logbuf, "Error: Invalid allele count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	      goto read_allele_freqs_ret_MALFORMED_INPUT_2;
	    }
	    uint64_t obs_ct = obs_ct_raw + ((uint64_t)maf_succ) * cur_allele_ct;
	    if (!obs_ct) {
	      goto read_allele_freqs_skip_variant;
	    }
	    obs_ct_recip = 1.0 / (double)((int64_t)obs_ct);
	  }
	  const uint32_t infer_freq_internal_idx = loaded_to_internal_allele_idx[infer_freq_loaded_idx];
	  if (cur_allele_ct == 2) {
	    // optimize common case
	    double known_freq_d = cur_allele_freqs[1 - infer_freq_internal_idx];
	    if (maf_succ) {
	      known_freq_d += 1;
	    }
	    double known_scaled_freq = known_freq_d * obs_ct_recip;
	    if (known_scaled_freq <= 1.0) {
	      if (infer_freq_internal_idx) {
		allele_freqs_write[0] = known_scaled_freq;
	      } else {
		allele_freqs_write[0] = 1.0 - known_scaled_freq;
	      }
	    } else if (known_scaled_freq <= (1.0 / 0.99)) {
	      if (infer_freq_internal_idx) {
		allele_freqs_write[0] = 1.0;
	      } else {
		allele_freqs_write[0] = 0.0;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Frequency/count too large on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	      goto read_allele_freqs_ret_MALFORMED_INPUT_2;
	    }
	  } else {
	    if (maf_succ) {
	      for (uint32_t internal_allele_idx = 0; internal_allele_idx < cur_allele_ct; ++internal_allele_idx) {
		cur_allele_freqs[internal_allele_idx] += 1;
	      }
	    }
	    cur_allele_freqs[infer_freq_internal_idx] = 0.0;
	    double known_freq_sum_d = 0.0;
	    for (uint32_t internal_allele_idx = 0; internal_allele_idx < cur_allele_ct; ++internal_allele_idx) {
	      known_freq_sum_d += cur_allele_freqs[internal_allele_idx];
	    }
	    double known_scaled_freq_sum = known_freq_sum_d * obs_ct_recip;
	    if (known_scaled_freq_sum <= (1.0 / 0.99)) {
	      if (known_scaled_freq_sum > 1.0) {
		// possible rounding error, rescale
		obs_ct_recip = 1.0 / known_scaled_freq_sum;
		known_scaled_freq_sum = 1.0;
	      }
	      const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
	      for (uint32_t internal_allele_idx = 0; internal_allele_idx < cur_allele_ct_m1; ++internal_allele_idx) {
		double dxx;
		if (internal_allele_idx == infer_freq_internal_idx) {
		  dxx = 1.0 - known_scaled_freq_sum;
		} else {
		  dxx = obs_ct_recip * cur_allele_freqs[internal_allele_idx];
		}
		allele_freqs_write[internal_allele_idx] = dxx;
	      }
	    } else {
	      sprintf(g_logbuf, "Error: Frequency/count too large on line %" PRIuPTR " of --read-freq file.\n", line_idx);
	      goto read_allele_freqs_ret_MALFORMED_INPUT_2;
	    }
	  }
	} else {
	  // complete frequency or count data
	  if (maf_succ) {
	    for (uint32_t internal_allele_idx = 0; internal_allele_idx < cur_allele_ct; ++internal_allele_idx) {
	      cur_allele_freqs[internal_allele_idx] += 1;
	    }
	  }
	  double tot_freq = 0.0;
	  for (uint32_t internal_allele_idx = 0; internal_allele_idx < cur_allele_ct; ++internal_allele_idx) {
	    tot_freq += cur_allele_freqs[internal_allele_idx];
	  }
	  if (tot_freq == 0.0) {
	    goto read_allele_freqs_skip_variant;
	  }
	  const double tot_freq_recip = 1.0 / tot_freq;
	  const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
	  for (uint32_t internal_allele_idx = 0; internal_allele_idx < cur_allele_ct_m1; ++internal_allele_idx) {
	    allele_freqs_write[internal_allele_idx] = tot_freq_recip * cur_allele_freqs[internal_allele_idx];
	  }
	}

	++loaded_variant_ct;
	if (!(loaded_variant_ct % 10000)) {
	  printf("\r--read-freq: Frequencies for %uk variants loaded.", loaded_variant_ct / 1000);
	  fflush(stdout);
	}
      }
      while (0) {
      read_allele_freqs_skip_variant:
	++skipped_variant_ct;
      }
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto read_allele_freqs_ret_READ_FAIL;
	}
	break;
      }
      ++line_idx;
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    }
    if (gzclose_null(&gz_infile)) {
      goto read_allele_freqs_ret_READ_FAIL;
    }
    putc_unlocked('\r', stdout);
    LOGPRINTF("--read-freq: Frequencies for %u variant%s loaded.\n", loaded_variant_ct, (loaded_variant_ct == 1)? "" : "s");
    if (skipped_variant_ct) {
      LOGERRPRINTFWW("Warning: %" PRIuPTR " entr%s skipped due to missing variant IDs, mismatching allele codes, and/or zero observations.\n", skipped_variant_ct, (skipped_variant_ct == 1)? "y" : "ies");
    }
  }
  while (0) {
  read_allele_freqs_ret_LONG_LINE:
    if (loadbuf_size == kMaxLongLine) {
      LOGERRPRINTF("Error: Line %" PRIuPTR " of --read-freq file is pathologically long.\n", line_idx);
      reterr = kPglRetMalformedInput;
      break;
    }
  read_allele_freqs_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  read_allele_freqs_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  read_allele_freqs_ret_UNRECOGNIZED_HEADER:
    logerrprint("Error: Unrecognized header line in --read-freq file.\n");
    reterr = kPglRetMalformedInput;
    break;
  read_allele_freqs_ret_MISSING_TOKENS:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of --read-freq file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  read_allele_freqs_ret_INVALID_FREQS:
    sprintf(g_logbuf, "Error: Invalid frequencies/counts on line %" PRIuPTR " of --read-freq file.\n", line_idx);
  read_allele_freqs_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
  read_allele_freqs_ret_MALFORMED_INPUT_2:
    logprint("\n");
    logerrprintb();
  read_allele_freqs_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 read_allele_freqs_ret_1:
  bigstack_reset(bigstack_mark);
  gzclose_cond(gz_infile);
  return reterr;
}

void compute_maj_alleles(const uintptr_t* variant_include, const uintptr_t* variant_allele_idxs, const double* allele_freqs, uint32_t variant_ct, alt_allele_ct_t* maj_alleles) {
  uint32_t cur_allele_ct = 2;
  uint32_t variant_uidx = 0;
  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    uintptr_t allele_idx_base;
    if (!variant_allele_idxs) {
      allele_idx_base = variant_uidx;
    } else {
      allele_idx_base = variant_allele_idxs[variant_uidx];
      cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - allele_idx_base;
      allele_idx_base -= variant_uidx;
    }
    const double* cur_allele_freqs_base = &(allele_freqs[allele_idx_base]);
    if (cur_allele_ct == 2) {
      maj_alleles[variant_uidx] = (cur_allele_freqs_base[0] < 0.5);
    } else {
      uint32_t maj_allele_idx = 0;
      double max_freq = cur_allele_freqs_base[0];
      double tot_alt_freq = max_freq;
      const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
      for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct_m1; ++allele_idx) {
	const double cur_freq = cur_allele_freqs_base[allele_idx];
	tot_alt_freq += cur_freq;
	if (cur_freq > max_freq) {
	  maj_allele_idx = allele_idx;
	  max_freq = cur_freq;
	}
      }
      if (max_freq + tot_alt_freq <= 1.0) {
	maj_allele_idx = cur_allele_ct_m1;
      }
      maj_alleles[variant_uidx] = maj_allele_idx;
    }
  }
}


// multithread globals
static pgen_reader_t** g_pgr_ptrs = nullptr;
static uintptr_t** g_genovecs = nullptr;
static uint32_t* g_read_variant_uidx_starts = nullptr;
static uintptr_t** g_missing_hc_acc1 = nullptr;
static uintptr_t** g_missing_dosage_acc1 = nullptr;
static uintptr_t** g_hethap_acc1 = nullptr;

static const uintptr_t* g_variant_include = nullptr;
static const chr_info_t* g_cip = nullptr;
static const uintptr_t* g_sex_male = nullptr;
static uint32_t g_raw_sample_ct = 0;
static uint32_t g_cur_block_size = 0;
static uint32_t g_calc_thread_ct = 0;
static pglerr_t g_error_ret = kPglRetSuccess;

THREAD_FUNC_DECL load_sample_missing_cts_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uintptr_t* variant_include = g_variant_include;
  const chr_info_t* cip = g_cip;
  const uintptr_t* sex_male = g_sex_male;
  const uint32_t raw_sample_ct = g_raw_sample_ct;
  const uint32_t raw_sample_ctaw = BITCT_TO_ALIGNED_WORDCT(raw_sample_ct);
  const uint32_t acc1_vec_ct = BITCT_TO_VECCT(raw_sample_ct);
  const uint32_t acc4_vec_ct = acc1_vec_ct * 4;
  const uint32_t acc8_vec_ct = acc1_vec_ct * 8;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const int32_t x_code = cip->xymt_codes[kChrOffsetX];
  const int32_t y_code = cip->xymt_codes[kChrOffsetY];
  uintptr_t* genovec_buf = g_genovecs[tidx];
  uintptr_t* missing_hc_acc1 = g_missing_hc_acc1[tidx];
  uintptr_t* missing_hc_acc4 = &(missing_hc_acc1[acc1_vec_ct * kWordsPerVec]);
  uintptr_t* missing_hc_acc8 = &(missing_hc_acc4[acc4_vec_ct * kWordsPerVec]);
  uintptr_t* missing_hc_acc32 = &(missing_hc_acc8[acc8_vec_ct * kWordsPerVec]);
  fill_ulong_zero(acc1_vec_ct * kWordsPerVec * 45, missing_hc_acc1);
  uintptr_t* missing_dosage_acc1 = nullptr;
  uintptr_t* missing_dosage_acc4 = nullptr;
  uintptr_t* missing_dosage_acc8 = nullptr;
  uintptr_t* missing_dosage_acc32 = nullptr;
  if (g_missing_dosage_acc1) {
    missing_dosage_acc1 = g_missing_dosage_acc1[tidx];
    missing_dosage_acc4 = &(missing_dosage_acc1[acc1_vec_ct * kWordsPerVec]);
    missing_dosage_acc8 = &(missing_dosage_acc4[acc4_vec_ct * kWordsPerVec]);
    missing_dosage_acc32 = &(missing_dosage_acc8[acc8_vec_ct * kWordsPerVec]);
    fill_ulong_zero(acc1_vec_ct * kWordsPerVec * 45, missing_dosage_acc1);
  }
  // could make this optional
  // (could technically make missing_hc optional too...)
  uintptr_t* hethap_acc1 = g_hethap_acc1[tidx];
  uintptr_t* hethap_acc4 = &(hethap_acc1[acc1_vec_ct * kWordsPerVec]);
  uintptr_t* hethap_acc8 = &(hethap_acc4[acc4_vec_ct * kWordsPerVec]);
  uintptr_t* hethap_acc32 = &(hethap_acc8[acc8_vec_ct * kWordsPerVec]);
  fill_ulong_zero(acc1_vec_ct * kWordsPerVec * 45, hethap_acc1);
  uint32_t all_ct_rem15 = 15;
  uint32_t all_ct_rem255d15 = 17;
  uint32_t hap_ct_rem15 = 15;
  uint32_t hap_ct_rem255d15 = 17;
  while (1) {
    pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
    const uint32_t is_last_block = g_is_last_thread_block;
    const uint32_t cur_block_size = g_cur_block_size;
    const uint32_t cur_idx_ct = (((tidx + 1) * cur_block_size) / calc_thread_ct) - ((tidx * cur_block_size) / calc_thread_ct);
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    uint32_t chr_end = 0;
    uintptr_t* cur_hets = nullptr;
    uint32_t is_diploid_x = 0;
    uint32_t is_y = 0;
    for (uint32_t cur_idx = 0; cur_idx < cur_idx_ct; ++cur_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	const uint32_t chr_fo_idx = get_variant_chr_fo_idx(cip, variant_uidx);
	const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	cur_hets = hethap_acc1;
	is_diploid_x = 0;
	is_y = 0;
	if (chr_idx == x_code) {
	  is_diploid_x = !is_set(cip->haploid_mask, 0);
	} else if (chr_idx == y_code) {
	  is_y = 1;
	} else {
	  if (!is_set(cip->haploid_mask, chr_idx)) {
	    cur_hets = nullptr;
	  }
	}
      }
      // could instead have missing_hc and (missing_hc - missing_dosage); that
      // has the advantage of letting you skip one of the two increment
      // operations when the variant is all hardcalls.
      pglerr_t reterr = pgr_read_missingness_multi(nullptr, nullptr, raw_sample_ct, variant_uidx, pgrp, missing_hc_acc1, missing_dosage_acc1, cur_hets, genovec_buf);
      if (reterr) {
	g_error_ret = reterr;
	break;
      }
      if (is_y) {
	bitvec_and(sex_male, raw_sample_ctaw, missing_hc_acc1);
	if (missing_dosage_acc1) {
	  bitvec_and(sex_male, raw_sample_ctaw, missing_dosage_acc1);
	}
      }
      unroll_incr_1_4(missing_hc_acc1, acc1_vec_ct, missing_hc_acc4);
      if (missing_dosage_acc1) {
	unroll_incr_1_4(missing_dosage_acc1, acc1_vec_ct, missing_dosage_acc4);
      }
      if (!(--all_ct_rem15)) {
	unroll_zero_incr_4_8(acc4_vec_ct, missing_hc_acc4, missing_hc_acc8);
	if (missing_dosage_acc1) {
	  unroll_zero_incr_4_8(acc4_vec_ct, missing_dosage_acc4, missing_dosage_acc8);
	}
	all_ct_rem15 = 15;
	if (!(--all_ct_rem255d15)) {
	  unroll_zero_incr_8_32(acc8_vec_ct, missing_hc_acc8, missing_hc_acc32);
	  if (missing_dosage_acc1) {
	    unroll_zero_incr_8_32(acc8_vec_ct, missing_dosage_acc8, missing_dosage_acc32);
	  }
	  all_ct_rem255d15 = 17;
	}
      }
      if (cur_hets) {
	if (is_diploid_x) {	  
	  bitvec_and(sex_male, raw_sample_ctaw, cur_hets);
	}
	unroll_incr_1_4(cur_hets, acc1_vec_ct, hethap_acc4);
	if (!(--hap_ct_rem15)) {
	  unroll_zero_incr_4_8(acc4_vec_ct, hethap_acc4, hethap_acc8);
	  hap_ct_rem15 = 15;
	  if (!(--hap_ct_rem255d15)) {
	    unroll_zero_incr_8_32(acc8_vec_ct, hethap_acc8, hethap_acc32);
	    hap_ct_rem255d15 = 17;
	  }
	}
      }
    }
    if (is_last_block) {
      unroll_incr_4_8(missing_hc_acc4, acc4_vec_ct, missing_hc_acc8);
      unroll_incr_8_32(missing_hc_acc8, acc8_vec_ct, missing_hc_acc32);
      if (missing_dosage_acc1) {
	unroll_incr_4_8(missing_dosage_acc4, acc4_vec_ct, missing_dosage_acc8);
	unroll_incr_8_32(missing_dosage_acc8, acc8_vec_ct, missing_dosage_acc32);
      }
      unroll_incr_4_8(hethap_acc4, acc4_vec_ct, hethap_acc8);
      unroll_incr_8_32(hethap_acc8, acc8_vec_ct, hethap_acc32);
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

pglerr_t load_sample_missing_cts(const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, uint32_t* sample_missing_hc_cts, uint32_t* sample_missing_dosage_cts, uint32_t* sample_hethap_cts) {
  assert(sample_missing_hc_cts || sample_missing_dosage_cts);
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (!variant_ct) {
      fill_uint_zero(raw_sample_ct, sample_missing_hc_cts);
      if (sample_missing_dosage_cts) {
	fill_uint_zero(raw_sample_ct, sample_missing_dosage_cts);
      }
      fill_uint_zero(raw_sample_ct, sample_hethap_cts);
      goto load_sample_missing_cts_ret_1;
    }
    // this doesn't seem to saturate below 35 threads
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    const uint32_t acc1_vec_ct = BITCT_TO_VECCT(raw_sample_ct);
    const uintptr_t acc1_alloc_cacheline_ct = DIV_UP(acc1_vec_ct * (45 * k1LU * kBytesPerVec), kCacheline);
    g_sex_male = sex_male;
    uintptr_t thread_alloc_cacheline_ct = 2 * acc1_alloc_cacheline_ct;
    g_missing_dosage_acc1 = nullptr;
    if (sample_missing_dosage_cts) {
      if (bigstack_alloc_ulp(calc_thread_ct, &g_missing_dosage_acc1)) {
	goto load_sample_missing_cts_ret_NOMEM;
      }
      thread_alloc_cacheline_ct += acc1_alloc_cacheline_ct;
    }
    if (bigstack_alloc_ulp(calc_thread_ct, &g_missing_hc_acc1) ||
	bigstack_alloc_ulp(calc_thread_ct, &g_hethap_acc1)) {
      goto load_sample_missing_cts_ret_NOMEM;
    }
    unsigned char* main_loadbufs[2];
    pthread_t* threads;
    uint32_t read_block_size;
    if (multithread_load_init(variant_include, raw_sample_ct, raw_variant_ct, pgr_alloc_cacheline_ct, thread_alloc_cacheline_ct, 0, pgfip, &calc_thread_ct, &g_genovecs, nullptr, nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto load_sample_missing_cts_ret_NOMEM;
    }
    const uintptr_t acc1_alloc = acc1_alloc_cacheline_ct * kCacheline;
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_missing_hc_acc1[tidx] = (uintptr_t*)bigstack_alloc_raw(acc1_alloc);
      if (g_missing_dosage_acc1) {
	g_missing_dosage_acc1[tidx] = (uintptr_t*)bigstack_alloc_raw(acc1_alloc);
      }
      g_hethap_acc1[tidx] = (uintptr_t*)bigstack_alloc_raw(acc1_alloc);
    }
    g_variant_include = variant_include;
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
    const uint32_t sample_ctv = acc1_vec_ct * kBitsPerVec;
    const uintptr_t acc32_offset = acc1_vec_ct * (13 * k1LU * kWordsPerVec);
    uint32_t* scrambled_missing_hc_cts = nullptr;
    uint32_t* scrambled_missing_dosage_cts = nullptr;
    uint32_t* scrambled_hethap_cts = nullptr;
    scrambled_missing_hc_cts = (uint32_t*)(&(g_missing_hc_acc1[0][acc32_offset]));
    if (g_missing_dosage_acc1) {
      scrambled_missing_dosage_cts = (uint32_t*)(&(g_missing_dosage_acc1[0][acc32_offset]));
    }
    scrambled_hethap_cts = (uint32_t*)(&(g_hethap_acc1[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx < calc_thread_ct; ++tidx) {
      uint32_t* thread_scrambled_missing_hc_cts = (uint32_t*)(&(g_missing_hc_acc1[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii < sample_ctv; ++uii) {
	scrambled_missing_hc_cts[uii] += thread_scrambled_missing_hc_cts[uii];
      }
      if (scrambled_missing_dosage_cts) {
	uint32_t* thread_scrambled_missing_dosage_cts = (uint32_t*)(&(g_missing_dosage_acc1[tidx][acc32_offset]));
	for (uint32_t uii = 0; uii < sample_ctv; ++uii) {
	  scrambled_missing_dosage_cts[uii] += thread_scrambled_missing_dosage_cts[uii];
	}
      }
      uint32_t* thread_scrambled_hethap_cts = (uint32_t*)(&(g_hethap_acc1[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii < sample_ctv; ++uii) {
	scrambled_hethap_cts[uii] += thread_scrambled_hethap_cts[uii];
      }
    }
    for (uint32_t sample_uidx = 0; sample_uidx < raw_sample_ct; ++sample_uidx) {
      const uint32_t scrambled_idx = scramble_1_4_8_32(sample_uidx);
      sample_missing_hc_cts[sample_uidx] = scrambled_missing_hc_cts[scrambled_idx];
      if (sample_missing_dosage_cts) {
	sample_missing_dosage_cts[sample_uidx] = scrambled_missing_dosage_cts[scrambled_idx];
      }
      sample_hethap_cts[sample_uidx] = scrambled_hethap_cts[scrambled_idx];
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

pglerr_t mind_filter(const uint32_t* sample_missing_cts, const uint32_t* sample_hethap_cts, const char* sample_ids, const char* sids, uint32_t raw_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t variant_ct, uint32_t variant_ct_y, double mind_thresh, uintptr_t* sample_include, uintptr_t* sex_male, uint32_t* sample_ct_ptr, char* outname, char* outname_end) {
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
      uint32_t cur_missing_geno_ct = sample_missing_cts[sample_uidx];
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
  uint32_t prefilter_variant_ct = *variant_ct_ptr;
  uint32_t x_start = 0xffffffffU;
  uint32_t x_end = 0xffffffffU;
  int32_t x_code;
  if (xymt_exists(cip, kChrOffsetX, &x_code)) {
    const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)x_code];
    x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
    x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
    // bugfix (4 Jun 2017): if no sex info available, need to skip chrX
    if (!hwe_x_pvals) {
      prefilter_variant_ct -= popcount_bit_idx(variant_include, x_start, x_end);
    }
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
	if (hwe_x_pvals) {
	  x_thresh = x_end;
	} else {
	  is_x = 0;
	  x_thresh = 0xffffffffU;
	  variant_uidx = next_set_unsafe(variant_include, x_end);
	}
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

void enforce_minor_freq_constraints(const uintptr_t* variant_allele_idxs, const uint64_t* founder_allele_dosages, const double* allele_freqs, double min_maf, double max_maf, uint64_t min_allele_dosage, uint64_t max_allele_dosage, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t prefilter_variant_ct = *variant_ct_ptr;
  uint32_t variant_uidx = 0;
  uint32_t removed_ct = 0;
  if ((min_maf != 0.0) || (max_maf != 1.0)) {
    // defend against floating point error
    min_maf *= 1.0 - kSmallEpsilon;
    max_maf *= 1.0 + kSmallEpsilon;
  } else {
    allele_freqs = nullptr;
  }
  const uint32_t dosage_filter = min_allele_dosage || (max_allele_dosage != (~0LLU));
  
  uint32_t cur_allele_ct = 2;
  for (uint32_t variant_idx = 0; variant_idx < prefilter_variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    uintptr_t variant_allele_idx_base;
    if (!variant_allele_idxs) {
      variant_allele_idx_base = 2 * variant_uidx;
    } else {
      variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
    }
    if (allele_freqs) {
      const double cur_nonmaj_freq = get_nonmaj_freq(&(allele_freqs[variant_allele_idx_base - variant_uidx]), cur_allele_ct);
      if ((cur_nonmaj_freq < min_maf) || (cur_nonmaj_freq > max_maf)) {
	CLEAR_BIT(variant_uidx, variant_include);
	++removed_ct;
	continue;
      }
    }
    if (dosage_filter) {
      const uint64_t* cur_founder_allele_dosages = &(founder_allele_dosages[variant_allele_idx_base]);
      uint64_t max_dosage = cur_founder_allele_dosages[0];
      uint64_t nonmaj_dosage;
      if (cur_allele_ct == 2) {
	nonmaj_dosage = MINV(max_dosage, cur_founder_allele_dosages[1]);
      } else {
	uint64_t tot_dosage = max_dosage;
	for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
	  const uint64_t cur_dosage = cur_founder_allele_dosages[allele_idx];
	  tot_dosage += cur_dosage;
	  if (cur_dosage > max_dosage) {
	    max_dosage = cur_dosage;
	  }
	}
	nonmaj_dosage = tot_dosage - max_dosage;
      }
      if ((nonmaj_dosage < min_allele_dosage) || (nonmaj_dosage > max_allele_dosage)) {
	CLEAR_BIT(variant_uidx, variant_include);
	++removed_ct;
      }
    }
  }
  LOGPRINTFWW("%u variant%s removed due to minor allele threshold(s) (--maf/--max-maf/--mac/--max-mac).\n", removed_ct, (removed_ct == 1)? "" : "s");
  *variant_ct_ptr -= removed_ct;
}

void enforce_mach_r2_thresh(const chr_info_t* cip, const double* mach_r2_vals, double mach_r2_min, double mach_r2_max, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t prefilter_variant_ct = *variant_ct_ptr;
  mach_r2_min *= 1 - kSmallEpsilon;
  mach_r2_max *= 1 + kSmallEpsilon;
  uint32_t removed_ct = 0;
  uint32_t variant_uidx = 0;
  const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];
  const uint32_t chr_ct = cip->chr_ct;
  uint32_t relevant_variant_ct = prefilter_variant_ct;
  for (uint32_t chr_fo_idx = 0; chr_fo_idx < chr_ct; ++chr_fo_idx) {
    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
    // skip X, Y, MT, other haploid
    if (is_set(cip->haploid_mask, chr_idx) || (chr_idx == ((uint32_t)mt_code))) {
      relevant_variant_ct -= popcount_bit_idx(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1]);
    }
  }
  uint32_t chr_fo_idx = 0xffffffffU;
  uint32_t chr_end = 0;
  for (uint32_t variant_idx = 0; variant_idx < relevant_variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    while (variant_uidx >= chr_end) {
      uint32_t chr_idx;
      do {
	chr_idx = cip->chr_file_order[++chr_fo_idx];
      } while (is_set(cip->haploid_mask, chr_idx) || (chr_idx == ((uint32_t)mt_code)));
      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      variant_uidx = next_set(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], chr_end);
    }
    const double cur_mach_r2 = mach_r2_vals[variant_uidx];
    if ((cur_mach_r2 < mach_r2_min) || (cur_mach_r2 > mach_r2_max)) {
      CLEAR_BIT(variant_uidx, variant_include);
      ++removed_ct;
    }
  }
  LOGPRINTF("--mach-r2-filter: %u variant%s removed.\n", removed_ct, (removed_ct == 1)? "" : "s");
  *variant_ct_ptr -= removed_ct;
}

#ifdef __cplusplus
} // namespace plink2
#endif
