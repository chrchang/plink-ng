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

#include "plink2_compress_stream.h"
#include "plink2_glm.h"
#include "plink2_matrix.h"
#include "plink2_stats.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void init_glm(glm_info_t* glm_info_ptr) {
  glm_info_ptr->flags = kfGlm0;
  glm_info_ptr->cols = kfGlmCol0;
  glm_info_ptr->mperm_ct = 0;
  glm_info_ptr->local_cat_ct = 0;
  glm_info_ptr->max_corr = 0.999;
  glm_info_ptr->condition_varname = nullptr;
  glm_info_ptr->condition_list_fname = nullptr;
  init_range_list(&(glm_info_ptr->parameters_range_list));
  init_range_list(&(glm_info_ptr->tests_range_list));
}

void cleanup_glm(glm_info_t* glm_info_ptr) {
  free_cond(glm_info_ptr->condition_varname);
  free_cond(glm_info_ptr->condition_list_fname);
  cleanup_range_list(&(glm_info_ptr->parameters_range_list));
  cleanup_range_list(&(glm_info_ptr->tests_range_list));
}

pglerr_t glm_local_init(const char* local_covar_fname, const char* local_pvar_fname, const char* local_psam_fname, const char* sample_ids, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const glm_info_t* glm_info_ptr, uint32_t raw_sample_ct, uintptr_t max_sample_id_blen, uint32_t raw_variant_ct, const uintptr_t** sample_include_ptr, const uintptr_t** sex_nm_ptr, const uintptr_t** sex_male_ptr, const uintptr_t** variant_include_ptr, uint32_t* sample_ct_ptr, uint32_t* variant_ct_ptr, gzFile* gz_local_covar_file_ptr, uint32_t** local_sample_uidx_order_ptr, uintptr_t** local_variant_include_ptr, uint32_t* local_sample_ct_ptr, uint32_t* local_variant_ctl_ptr, uint32_t* local_covar_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  gzFile gz_infile = nullptr;
  uintptr_t line_idx = 0;
  uintptr_t loadbuf_size = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    // 1. read .psam/.fam file, update sample_ct, initialize
    //    local_sample_uidx_order (use open_and_load_xid_header()?)
    reterr = gzopen_read_checked(local_psam_fname, &gz_infile);
    if (reterr) {
      goto glm_local_init_ret_1;
    }
    loadbuf_size = bigstack_left() / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto glm_local_init_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kEndAllocAlign);
    }
    char* loadbuf = (char*)bigstack_end_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    char* loadbuf_first_token;
    uint32_t is_header_line;
    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
	  goto glm_local_init_ret_READ_FAIL;
	}
	sprintf(g_logbuf, "Error: %s is empty.\n", local_psam_fname);
	goto glm_local_init_ret_MALFORMED_INPUT_WW;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto glm_local_init_ret_LONG_LINE_PSAM;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      is_header_line = (loadbuf_first_token[0] == '#');
    } while (is_header_line && strcmp_se(&(loadbuf_first_token[1]), "FID", 3) && strcmp_se(&(loadbuf_first_token[1]), "IID", 3));
    xid_mode_t xid_mode = kfXidModeFidiid;
    if (is_header_line) {
      if (loadbuf_first_token[1] == 'I') {
	xid_mode = kfXidModeIid;
      }
      *loadbuf_first_token = '\0';
    }
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    char* sorted_sample_idbox;
    uint32_t* sample_id_map;
    uintptr_t* new_sample_include;
    char* idbuf;
    if (bigstack_end_alloc_c(orig_sample_ct * max_sample_id_blen, &sorted_sample_idbox) ||
	bigstack_end_alloc_ui(orig_sample_ct, &sample_id_map) ||
	bigstack_end_calloc_ul(raw_sample_ctl, &new_sample_include) ||
	bigstack_end_alloc_c(max_sample_id_blen, &idbuf)) {
      goto glm_local_init_ret_NOMEM;
    }
    // (don't permit duplicate FID+IID for now, but maybe we'll want to use
    // xid interface later?)
    reterr = copy_sort_strbox_subset_noalloc(*sample_include_ptr, sample_ids, orig_sample_ct, max_sample_id_blen, 0, 0, 0, sorted_sample_idbox, sample_id_map);
    if (reterr) {
      goto glm_local_init_ret_1;
    }
    uint32_t* local_sample_uidx_order = (uint32_t*)g_bigstack_base;
    uintptr_t max_local_sample_ct = round_down_pow2(bigstack_left(), kCacheline) / sizeof(int32_t);
#ifdef __LP64__
    if (max_local_sample_ct > kMaxLongLine / 2) {
      max_local_sample_ct = kMaxLongLine / 2;
    }
#endif
    uintptr_t local_sample_ct = 0;
    while (1) {
      if (!is_eoln_kns(*loadbuf_first_token)) {
	if (local_sample_ct == max_local_sample_ct) {
#ifdef __LP64__
	  if (local_sample_ct == kMaxLongLine / 2) {
	    sprintf(g_logbuf, "Error: Too many samples in %s.\n", local_psam_fname);
	    goto glm_local_init_ret_MALFORMED_INPUT_WW;
	  }
#endif
	  goto glm_local_init_ret_NOMEM;
	}
	char* read_ptr = loadbuf_first_token;
	uint32_t sample_uidx;
	if (!sorted_xidbox_read_find(sorted_sample_idbox, sample_id_map, max_sample_id_blen, orig_sample_ct, 0, xid_mode, &read_ptr, &sample_uidx, idbuf)) {
	  if (is_set(new_sample_include, sample_uidx)) {
	    char* first_tab = (char*)rawmemchr(idbuf, '\t');
	    *first_tab = ' ';
	    sprintf(g_logbuf, "Error: Duplicate ID '%s' in %s.\n", idbuf, local_psam_fname);
	    goto glm_local_init_ret_MALFORMED_INPUT_WW;
	  }
	  set_bit(sample_uidx, new_sample_include);
	  local_sample_uidx_order[local_sample_ct] = sample_uidx;
	} else {
	  if (!read_ptr) {
	    sprintf(g_logbuf, "Error: Fewer tokens than expected on line %" PRIuPTR " of %s.\n", line_idx, local_psam_fname);
	    goto glm_local_init_ret_MALFORMED_INPUT_WW;
	  }
	  local_sample_uidx_order[local_sample_ct] = 0xffffffffU;
	}
	++local_sample_ct;
      }
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto glm_local_init_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto glm_local_init_ret_LONG_LINE_PSAM;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (loadbuf_first_token[0] == '#') {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, local_psam_fname);
	goto glm_local_init_ret_MALFORMED_INPUT_WW;
      }
    }
    if (gzclose_null(&gz_infile)) {
      goto glm_local_init_ret_READ_FAIL;
    }
    bigstack_finalize_ui(local_sample_uidx_order, local_sample_ct);
    *local_sample_uidx_order_ptr = local_sample_uidx_order;
    *local_sample_ct_ptr = local_sample_ct;
    const uint32_t new_sample_ct = popcount_longs(new_sample_include, raw_sample_ctl);
    assert(new_sample_ct <= orig_sample_ct);
    if (new_sample_ct < orig_sample_ct) {
      uintptr_t* sample_include_copy;
      uintptr_t* sex_nm_copy;
      uintptr_t* sex_male_copy;      
      if (bigstack_alloc_ul(raw_sample_ctl, &sample_include_copy) ||
	  bigstack_alloc_ul(raw_sample_ctl, &sex_nm_copy) ||
	  bigstack_alloc_ul(raw_sample_ctl, &sex_male_copy)) {
	goto glm_local_init_ret_NOMEM;
      }
      memcpy(sample_include_copy, new_sample_include, raw_sample_ctl * sizeof(intptr_t));
      bitvec_and_copy(sample_include_copy, *sex_nm_ptr, raw_sample_ctl, sex_nm_copy);
      *sex_nm_ptr = sex_nm_copy;
      bitvec_and_copy(sample_include_copy, *sex_male_ptr, raw_sample_ctl, sex_male_copy);
      *sex_male_ptr = sex_male_copy;
      *sample_include_ptr = sample_include_copy;
      bigstack_end_reset(loadbuf);
      uint32_t* sample_include_cumulative_popcounts;
      if (bigstack_end_alloc_ui(raw_sample_ctl, &sample_include_cumulative_popcounts)) {
	goto glm_local_init_ret_NOMEM;
      }
      fill_cumulative_popcounts(sample_include_copy, raw_sample_ctl, sample_include_cumulative_popcounts);
      *sample_ct_ptr = new_sample_ct;
    }
    bigstack_end_reset(loadbuf);

    // 2. read .pvar/.bim file, update variant_ct, initialize
    //    local_variant_ct/local_variant_include.
    reterr = gzopen_read_checked(local_pvar_fname, &gz_infile);
    if (reterr) {
      goto glm_local_init_ret_1;
    }
    line_idx = 0;
    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
	  goto glm_local_init_ret_READ_FAIL;
	}
	sprintf(g_logbuf, "Error: %s is empty.\n", local_pvar_fname);
	goto glm_local_init_ret_MALFORMED_INPUT_WW;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto glm_local_init_ret_LONG_LINE_PVAR;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      is_header_line = (loadbuf_first_token[0] == '#');
    } while (is_header_line && strcmp_se(&(loadbuf_first_token[1]), "CHROM", 5));
    uint32_t col_skips[2];
    uint32_t col_types[2];
    // uint32_t relevant_postchr_col_ct = 2;
    if (is_header_line) {
      // parse header
      // [-1] = #CHROM (must be first column)
      // [0] = POS
      // [1] = ID
      // don't care about the rest
      uint32_t col_idx = 0;
      char* token_end = &(loadbuf_first_token[6]);
      uint32_t found_header_bitset = 0;
      uint32_t relevant_postchr_col_ct = 0;
      char* loadbuf_iter;
      while (1) {
	loadbuf_iter = skip_initial_spaces(token_end);
	if (is_eoln_kns(*loadbuf_iter)) {
	  break;
	}
	++col_idx;
	token_end = token_endnn(loadbuf_iter);
	const uint32_t token_slen = (uintptr_t)(token_end - loadbuf_iter);
	uint32_t cur_col_type;
	if ((token_slen == 3) && (!memcmp(loadbuf_iter, "POS", 3))) {
	  cur_col_type = 0;
	} else if ((token_slen == 2) && (!memcmp(loadbuf_iter, "ID", 2))) {
	  cur_col_type = 1;
	} else {
	  continue;
	}
	const uint32_t cur_col_type_shifted = 1 << cur_col_type;
	if (found_header_bitset & cur_col_type_shifted) {
	  *token_end = '\0';
	  sprintf(g_logbuf, "Error: Duplicate column header '%s' on line %" PRIuPTR " of %s.\n", loadbuf_iter, line_idx, local_pvar_fname);
	  goto glm_local_init_ret_MALFORMED_INPUT_WW;
	}
	found_header_bitset |= cur_col_type_shifted;
	col_skips[relevant_postchr_col_ct] = col_idx;
	col_types[relevant_postchr_col_ct++] = cur_col_type;
      }
      if (found_header_bitset != 3) {
	sprintf(g_logbuf, "Error: Missing column header(s) on line %" PRIuPTR " of %s. (POS and ID are required.)\n", line_idx, local_pvar_fname);
	goto glm_local_init_ret_MALFORMED_INPUT_WW;
      }
      for (uint32_t rpc_col_idx = relevant_postchr_col_ct - 1; rpc_col_idx; --rpc_col_idx) {
	col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
      }
      loadbuf_first_token[0] = '\0';
    } else {
      col_types[0] = 1;
      col_types[1] = 0;
      col_skips[0] = 1;
      // CM column may be omitted
      char* loadbuf_iter = next_token_mult(loadbuf_first_token, 4);
      if (!loadbuf_iter) {
	goto glm_local_init_ret_MISSING_TOKENS_PVAR;
      }
      loadbuf_iter = next_token(loadbuf_iter);
      if (!loadbuf_iter) {
	// #CHROM ID POS ALT REF
	col_skips[1] = 1;
      } else {
	// #CHROM ID CM POS ALT REF
	col_skips[1] = 2;
      }
    }
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uintptr_t* new_variant_include;
    if (bigstack_end_calloc_ul(raw_variant_ctl, &new_variant_include)) {
      goto glm_local_init_ret_NOMEM;
    }
    uint32_t max_local_variant_ct = 0x7ffffffd;
    if (bigstack_left() < (0x80000000U / CHAR_BIT)) {
      max_local_variant_ct = round_down_pow2(bigstack_left(), kCacheline) * CHAR_BIT;
    }
    const uintptr_t* orig_variant_include = *variant_include_ptr;
    uintptr_t* local_variant_include = (uintptr_t*)g_bigstack_base;
    *local_variant_include_ptr = local_variant_include;
    uint32_t local_variant_ct = 0;
    uint32_t new_variant_ct = 0;
    uint32_t prev_variant_uidx = next_set_unsafe(orig_variant_include, 0);
    uint32_t chr_fo_idx = get_variant_chr_fo_idx(cip, prev_variant_uidx);
    uint32_t prev_chr_code = cip->chr_file_order[chr_fo_idx];
    uint32_t prev_bp = variant_bps[prev_variant_uidx];
    uint32_t chr_end = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[prev_chr_code] + 1];
    while (1) {
      if (!is_eoln_kns(*loadbuf_first_token)) {
	if (local_variant_ct == max_local_variant_ct) {
	  if (max_local_variant_ct == 0x7ffffffd) {
	    sprintf(g_logbuf, "Error: Too many samples in %s.\n", local_pvar_fname);
	    goto glm_local_init_ret_MALFORMED_INPUT_WW;
	  }
	  goto glm_local_init_ret_NOMEM;
	}
	if (!(local_variant_ct % kBitsPerWord)) {
	  local_variant_include[local_variant_ct / kBitsPerWord] = 0;
	}
	char* loadbuf_iter = token_endnn(loadbuf_first_token);
	// #CHROM
	if (!(*loadbuf_iter)) {
	  goto glm_local_init_ret_MISSING_TOKENS_PVAR;
	}
	{
	  const int32_t cur_chr_code = get_chr_code_counted(cip, (uintptr_t)(loadbuf_iter - loadbuf_first_token), loadbuf_first_token);
	  if (cur_chr_code < 0) {
	    goto glm_local_init_skip_variant;
	  }
	  if ((uint32_t)cur_chr_code != prev_chr_code) {
	    uint32_t first_variant_uidx_in_chr = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[(uint32_t)cur_chr_code]];
	    if (first_variant_uidx_in_chr < prev_variant_uidx) {
	      if (new_variant_ct) {
		// not worth the trouble of handling this
		sprintf(g_logbuf, "Error: Chromosome order in %s is different from chromosome order in main dataset.\n", local_pvar_fname);
		goto glm_local_init_ret_INCONSISTENT_INPUT_WW;
	      }
	      goto glm_local_init_skip_variant;
	    }
	    prev_variant_uidx = next_set(orig_variant_include, first_variant_uidx_in_chr, raw_variant_ct);
	    if (prev_variant_uidx == raw_variant_ct) {
	      break;
	    }
	    chr_fo_idx = get_variant_chr_fo_idx(cip, prev_variant_uidx);
	    prev_chr_code = cip->chr_file_order[chr_fo_idx];
	    prev_bp = variant_bps[prev_variant_uidx];
	    chr_end = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[prev_chr_code] + 1];
	    if ((uint32_t)cur_chr_code != prev_chr_code) {
	      goto glm_local_init_skip_variant;
	    }
	  }
	  char* token_ptrs[2];
	  uint32_t token_slens[2];
	  for (uint32_t rpc_col_idx = 0; rpc_col_idx < 2; ++rpc_col_idx) {
	    const uint32_t cur_col_type = col_types[rpc_col_idx];
	    loadbuf_iter = next_token_mult(loadbuf_iter, col_skips[rpc_col_idx]);
	    if (!loadbuf_iter) {
	      goto glm_local_init_ret_MISSING_TOKENS_PVAR;
	    }
	    token_ptrs[cur_col_type] = loadbuf_iter;
	    char* token_end = token_endnn(loadbuf_iter);
	    token_slens[cur_col_type] = (uintptr_t)(token_end - loadbuf_iter);
	    loadbuf_iter = token_end;
	  }
	  // POS
	  int32_t cur_bp;
	  if (scan_int_abs_defcap(token_ptrs[0], &cur_bp)) {
	    sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, local_pvar_fname);
	    goto glm_local_init_ret_MALFORMED_INPUT_WW;
	  }
	  if (cur_bp < (int32_t)prev_bp) {
	    goto glm_local_init_skip_variant;
	  }
	  if ((uint32_t)cur_bp > prev_bp) {
	    do {
	      prev_variant_uidx = next_set(orig_variant_include, prev_variant_uidx + 1, raw_variant_ct);
	    } while ((prev_variant_uidx < chr_end) && ((uint32_t)cur_bp > variant_bps[prev_variant_uidx]));
	    if (prev_variant_uidx >= chr_end) {
	      goto glm_local_init_skip_variant_and_update_chr;
	    }
	    prev_bp = variant_bps[prev_variant_uidx];
	  }
	  if ((uint32_t)cur_bp == prev_bp) {
	    // ID
	    // note that if there are two same-position variants which appear
	    // in a different order in the main dataset and the local-pvar
	    // file, one will be skipped.  (probably want to add a warning in
	    // this case.)
	    char* cur_variant_id = token_ptrs[1];
	    cur_variant_id[token_slens[1]] = '\0';
	    do {
	      char* loaded_variant_id = variant_ids[prev_variant_uidx];
	      if (!strcmp(cur_variant_id, loaded_variant_id)) {
		if (is_set(new_variant_include, prev_variant_uidx)) {
		  sprintf(g_logbuf, "Error: Duplicate ID (with duplicate CHROM/POS) '%s' in %s.\n", cur_variant_id, local_pvar_fname);
		  goto glm_local_init_ret_MALFORMED_INPUT_WW;
		}
		set_bit(prev_variant_uidx, new_variant_include);
		++new_variant_ct;
		set_bit(local_variant_ct, local_variant_include);
		break;
	      }
	      prev_variant_uidx = next_set(orig_variant_include, prev_variant_uidx + 1, raw_variant_ct);
	      if (prev_variant_uidx >= chr_end) {
		goto glm_local_init_skip_variant_and_update_chr;
	      }
	      prev_bp = variant_bps[prev_variant_uidx];
	    } while ((uint32_t)cur_bp == prev_bp);
	  }
	}
	if (0) {
	glm_local_init_skip_variant_and_update_chr:
	  if (prev_variant_uidx == raw_variant_ct) {
	    break;
	  }
	  chr_fo_idx = get_variant_chr_fo_idx(cip, prev_variant_uidx);
	  prev_chr_code = cip->chr_file_order[chr_fo_idx];
	  prev_bp = variant_bps[prev_variant_uidx];
	  chr_end = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[prev_chr_code] + 1];
	}
      glm_local_init_skip_variant:
	++local_variant_ct;
      }
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto glm_local_init_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto glm_local_init_ret_LONG_LINE_PVAR;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (loadbuf_first_token[0] == '#') {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #CHROM header line is present it must denote the end of the header block.)\n", line_idx, local_pvar_fname);
	goto glm_local_init_ret_MALFORMED_INPUT_WW;
      }
    }
    if (gzclose_null(&gz_infile)) {
      goto glm_local_init_ret_READ_FAIL;
    }
    bigstack_finalize_ul(local_variant_include, BITCT_TO_WORDCT(local_variant_ct));
    *local_variant_ctl_ptr = BITCT_TO_WORDCT(local_variant_ct);
    assert(new_variant_ct <= *variant_ct_ptr);
    if (new_variant_ct < *variant_ct_ptr) {
      uintptr_t* variant_include_copy;
      if (bigstack_alloc_ul(raw_variant_ctl, &variant_include_copy)) {
	goto glm_local_init_ret_NOMEM;
      }
      memcpy(variant_include_copy, new_variant_include, raw_variant_ctl * sizeof(intptr_t));
      *variant_include_ptr = variant_include_copy;
      *variant_ct_ptr = new_variant_ct;
    }

    // 3. if not local-cats=, scan first line of local-covar= file to determine
    //    covariate count
    reterr = gzopen_read_checked(local_covar_fname, gz_local_covar_file_ptr);
    if (reterr) {
      goto glm_local_init_ret_1;
    }
    uint32_t local_covar_ct;
    if (!glm_info_ptr->local_cat_ct) {
      line_idx = 1;
      if (!gzgets(*gz_local_covar_file_ptr, loadbuf, loadbuf_size)) {
	if (!gzeof(*gz_local_covar_file_ptr)) {
	  goto glm_local_init_ret_READ_FAIL;
	}
	sprintf(g_logbuf, "Error: %s is empty.\n", local_covar_fname);
	goto glm_local_init_ret_MALFORMED_INPUT_WW;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size != kMaxLongLine) {
	  goto glm_local_init_ret_NOMEM;
	}
	sprintf(g_logbuf, "Error: Line 1 of %s is pathologically long.\n", local_covar_fname);
	goto glm_local_init_ret_MALFORMED_INPUT_WW;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      const uint32_t token_ct = count_tokens(loadbuf_first_token);
      local_covar_ct = token_ct / local_sample_ct;
      if (local_covar_ct * local_sample_ct != token_ct) {
	sprintf(g_logbuf, "Error: Unexpected token count on line 1 of %s (%u, %smultiple of %" PRIuPTR " expected).\n", local_covar_fname, token_ct, (token_ct == local_sample_ct)? "larger " : "", local_sample_ct);
	goto glm_local_init_ret_MALFORMED_INPUT_WW;
      }
      if (glm_info_ptr->flags & kfGlmLocalOmitLast) {
	if (local_covar_ct == 1) {
	  logerrprint("Error: --glm 'local-omit-last' modifier cannot be used when there is only one\nlocal covariate.\n");
	  goto glm_local_init_ret_INCONSISTENT_INPUT;
	}
	LOGPRINTF("--glm local-covar=: %u local covariates present, %u used.\n", local_covar_ct, local_covar_ct - 1);
	--local_covar_ct;
      } else {
	LOGPRINTF("--glm local-covar=: %u local covariate%s present.\n", local_covar_ct, (local_covar_ct == 1)? "" : "s");
      }
    } else {
      local_covar_ct = glm_info_ptr->local_cat_ct - 1;
      if (local_covar_ct * ((uint64_t)local_sample_ct) > kMaxLongLine / 2) {
	sprintf(g_logbuf, "Error: [# samples in %s] * [# categories - 1] too large (limited to %u).\n", local_covar_fname, kMaxLongLine / 2);
	goto glm_local_init_ret_MALFORMED_INPUT_WW;
      }
    }
    *local_covar_ct_ptr = local_covar_ct;
    bigstack_mark = g_bigstack_base;
  }
  while (0) {
  glm_local_init_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  glm_local_init_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  glm_local_init_ret_LONG_LINE_PSAM:
    if (loadbuf_size != kMaxLongLine) {
      reterr = kPglRetNomem;
      break;
    }
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, local_psam_fname);
  glm_local_init_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  glm_local_init_ret_LONG_LINE_PVAR:
    if (loadbuf_size != kMaxLongLine) {
      reterr = kPglRetNomem;
      break;
    }
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, local_pvar_fname);
    reterr = kPglRetMalformedInput;
    break;
  glm_local_init_ret_MISSING_TOKENS_PVAR:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, local_pvar_fname);
    reterr = kPglRetMalformedInput;
    break;
  glm_local_init_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  glm_local_init_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 glm_local_init_ret_1:
  gzclose_cond(gz_infile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}


// Returns 1 if phenotype fully separates the covariate, and the non-Firth
// logistic regression should be aborted.
// Otherwise, if we have quasi-separation, use Stata approach of throwing out
// the covariate, and keeping only samples which have the one covariate value
// present in both cases and controls.
// Note that covar_ct is not a parameter; caller is responsible for
// re-popcounting covar_include.
boolerr_t check_for_and_handle_separated_covar(const uintptr_t* pheno_cc, const pheno_col_t* covar_cols, uint32_t raw_sample_ctl, uint32_t covar_uidx, uintptr_t* cur_sample_include, uintptr_t* covar_include, uint32_t* sample_ct_ptr, uintptr_t* cat_covar_wkspace) {
  uint32_t sample_ct = *sample_ct_ptr;
  if (sample_ct < 2) {
    return 1;
  }
  const pheno_col_t* cur_covar_col = &(covar_cols[covar_uidx]);
  if (cur_covar_col->type_code == kPhenoDtypeOther) {
    return 0;
  }
  const uint32_t first_sample_uidx = next_set_unsafe(cur_sample_include, 0);
  if (cur_covar_col->type_code == kPhenoDtypeQt) {
    const double* covar_vals = cur_covar_col->data.qt;
    double cur_covar_val = covar_vals[first_sample_uidx];
    double ctrl_min;
    double ctrl_max;
    double case_min;
    double case_max;
    if (is_set(pheno_cc, first_sample_uidx)) {
      ctrl_min = DBL_MAX;
      ctrl_max = -DBL_MAX;
      case_min = cur_covar_val;
      case_max = cur_covar_val;
    } else {
      ctrl_min = cur_covar_val;
      ctrl_max = cur_covar_val;
      case_min = DBL_MAX;
      case_max = -DBL_MAX;
    }
    uint32_t sample_uidx = first_sample_uidx + 1;
    for (uint32_t sample_idx = 1; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_unsafe_ck(cur_sample_include, &sample_uidx);
      cur_covar_val = covar_vals[sample_uidx];
      if (IS_SET(pheno_cc, sample_uidx)) {
	if (cur_covar_val < case_min) {
	  if ((ctrl_min < case_max) && (cur_covar_val < ctrl_max)) {
	    return 0;
	  }
	  case_min = cur_covar_val;
	} else if (cur_covar_val > case_max) {
	  if ((ctrl_max > case_min) && (cur_covar_val > ctrl_min)) {
	    return 0;
	  }
	  case_max = cur_covar_val;
	}
      } else {
	if (cur_covar_val < ctrl_min) {
	  if ((case_min < ctrl_max) && (cur_covar_val < case_max)) {
	    return 0;
	  }
	  ctrl_min = cur_covar_val;
	} else if (cur_covar_val > ctrl_max) {
	  if ((case_max > ctrl_min) && (cur_covar_val > case_min)) {
	    return 0;
	  }
	  ctrl_max = cur_covar_val;
	}
      }
    }
    if ((case_min > ctrl_max) || (ctrl_min > case_max)) {
      // fully separated
      return 1;
    }
    // quasi-separated
    const double covar_val_keep = (case_min == ctrl_max)? case_min : case_max;
    sample_uidx = first_sample_uidx;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_unsafe_ck(cur_sample_include, &sample_uidx);
      if (covar_vals[sample_uidx] != covar_val_keep) {
	clear_bit(sample_uidx, cur_sample_include);
      }
    }
    *sample_ct_ptr = popcount_longs(cur_sample_include, raw_sample_ctl);
    clear_bit(covar_uidx, covar_include);
    return 0;
  }
  const uint32_t nonnull_cat_ct = cur_covar_col->nonnull_category_ct;
  const uint32_t cur_word_ct = 1 + nonnull_cat_ct / kBitsPerWordD2;
  fill_ulong_zero(cur_word_ct, cat_covar_wkspace);
  const uint32_t* covar_vals = cur_covar_col->data.cat;
  // If no remaining categories have both cases and controls, we have complete
  // separation.
  // If some do and some do not, we have quasi-complete separation, and must
  // remove samples in the all-case and all-control categories.
  uint32_t sample_uidx = first_sample_uidx;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    next_set_unsafe_ck(cur_sample_include, &sample_uidx);
    const uint32_t cur_cat_idx = covar_vals[sample_uidx];
    // Odd bits represent presence of a case, even bits represent presence of a
    // control.
    const uint32_t is_case = is_set(pheno_cc, sample_uidx);
    set_bit(cur_cat_idx * 2 + is_case, cat_covar_wkspace);
  }
  uint32_t case_and_ctrl_cat_ct = 0;
  uint32_t pheno_by_cat_ct = 0;
  for (uint32_t widx = 0; widx < cur_word_ct; ++widx) {
    const uintptr_t cur_word = cat_covar_wkspace[widx];
    case_and_ctrl_cat_ct += popcount01_long(cur_word & (cur_word >> 1) & kMask5555);
    pheno_by_cat_ct += popcount_long(cur_word);
  }
  if (!case_and_ctrl_cat_ct) {
    // fully separated
    return 1;
  }
  if ((case_and_ctrl_cat_ct > 1) && (case_and_ctrl_cat_ct == pheno_by_cat_ct * 2)) {
    // all categories contain both cases and controls.
    return 0;
  }
  // more than one category contains both cases and controls (so we don't need
  // to remove the categorical covariate), but at least one does not, so we
  // still have to prune some samples.
  for (uint32_t widx = 0; widx < cur_word_ct; ++widx) {
    const uintptr_t cur_word = cat_covar_wkspace[widx];
    cat_covar_wkspace[widx] = cur_word & (cur_word >> 1) & kMask5555;
  }
  sample_uidx = first_sample_uidx;
  for (uint32_t sample_idx = 0; sample_uidx < sample_ct; ++sample_idx, ++sample_uidx) {
    next_set_unsafe_ck(cur_sample_include, &sample_uidx);
    if (!is_set(cat_covar_wkspace, covar_vals[sample_uidx] * 2)) {
      clear_bit(sample_uidx, cur_sample_include);
    }
  }
  *sample_ct_ptr = popcount_longs(cur_sample_include, raw_sample_ctl);
  if (case_and_ctrl_cat_ct == 1) {
    clear_bit(covar_uidx, covar_include);
  }
  return 0;
}

boolerr_t glm_determine_covars(const uintptr_t* pheno_cc, const uintptr_t* initial_covar_include, const pheno_col_t* covar_cols, uint32_t raw_sample_ct, uint32_t raw_covar_ctl, uint32_t initial_covar_ct, uint32_t covar_max_nonnull_cat_ct, uint32_t is_sometimes_firth, uintptr_t* cur_sample_include, uintptr_t* covar_include, uint32_t* sample_ct_ptr, uint32_t* covar_ct_ptr, uint32_t* extra_cat_ct_ptr, uint32_t* separation_warning_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  boolerr_t reterr = 0;
  {
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    uintptr_t* sample_include_backup;
    if (bigstack_alloc_ul(raw_sample_ctl, &sample_include_backup)) {
      goto glm_determine_covars_ret_NOMEM;
    }
    memcpy(sample_include_backup, cur_sample_include, raw_sample_ctl * sizeof(intptr_t));
    memcpy(covar_include, initial_covar_include, raw_covar_ctl * sizeof(intptr_t));

    // 1. Determine samples for which all phenotype/covariate values are
    //    present, then provisionally remove the covariates which are constant
    //    over that set in linear case, or produce separation in logistic case
    uint32_t covar_uidx = 0;
    for (uint32_t covar_idx = 0; covar_idx < initial_covar_ct; ++covar_idx, ++covar_uidx) {
      next_set_unsafe_ck(initial_covar_include, &covar_uidx);
      if (covar_cols[covar_uidx].nonmiss) {
	bitvec_and(covar_cols[covar_uidx].nonmiss, raw_sample_ctl, cur_sample_include);
      }
    }
    uint32_t prev_sample_ct = popcount_longs(cur_sample_include, raw_sample_ctl);
    covar_uidx = 0;
    for (uint32_t covar_idx = 0; covar_idx < initial_covar_ct; ++covar_idx, ++covar_uidx) {
      next_set_unsafe_ck(initial_covar_include, &covar_uidx);
      if ((covar_cols[covar_uidx].type_code != kPhenoDtypeOther) && is_const_covar(&(covar_cols[covar_uidx]), cur_sample_include, prev_sample_ct)) {
	clear_bit(covar_uidx, covar_include);
      }
    }
    uint32_t covar_ct = popcount_longs(covar_include, raw_covar_ctl);
    if (covar_ct != initial_covar_ct) {
      // 2. If any covariates were removed, redetermine the samples for which
      //    all phenotype/covariate values are present.  If there are more than
      //    before, add back any provisionally removed covariates which don't
      //    have any missing values in the new sample set, and are now
      //    nonconstant.
      //    Categorical covariates should behave just like they had been
      //    pre-split into n-1 0/1 indicator variables.
      memcpy(cur_sample_include, sample_include_backup, raw_sample_ctl * sizeof(intptr_t));
      covar_uidx = 0;
      for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx, ++covar_uidx) {
	next_set_unsafe_ck(covar_include, &covar_uidx);
	if (covar_cols[covar_uidx].nonmiss) {
	  bitvec_and(covar_cols[covar_uidx].nonmiss, raw_sample_ctl, cur_sample_include);
	}
      }
      uint32_t new_sample_ct = popcount_longs(cur_sample_include, raw_sample_ctl);
      if (new_sample_ct > prev_sample_ct) {
	prev_sample_ct = new_sample_ct;
	covar_uidx = 0;
	for (uint32_t covar_idx = 0; covar_idx < initial_covar_ct; ++covar_idx, ++covar_uidx) {
	  next_set_unsafe_ck(initial_covar_include, &covar_uidx);
	  if (!is_set(covar_include, covar_uidx)) {
	    const pheno_col_t* cur_covar_col = &(covar_cols[covar_uidx]);
	    if (popcount_longs_intersect(cur_sample_include, cur_covar_col->nonmiss, raw_sample_ctl) == prev_sample_ct) {
	      if (!is_const_covar(cur_covar_col, cur_sample_include, prev_sample_ct)) {
		set_bit(covar_uidx, covar_include);
	      }
	    }
	  }
	}
	covar_ct = popcount_longs(covar_include, raw_covar_ctl);
      }
    }
    if (!covar_ct) {
      *sample_ct_ptr = prev_sample_ct;
      goto glm_determine_covars_ret_NOCOVAR;
    }
    if (prev_sample_ct < 3) {
      goto glm_determine_covars_ret_SKIP;
    }
    // 3. if quantitative trait, remove samples corresponding to single-element
    //    categories or constant-except-for-one-sample regular covariates.
    uint32_t sample_ct = prev_sample_ct;
    uint32_t extra_cat_ct = 0;
    bigstack_reset(sample_include_backup);
    uint32_t first_sample_uidx = next_set_unsafe(cur_sample_include, 0);
    if (!pheno_cc) {
      uintptr_t* cat_one_obs = nullptr;
      uintptr_t* cat_two_or_more_obs = nullptr;
      uint32_t* cat_first_sample_uidxs = nullptr;
      if (covar_max_nonnull_cat_ct) {
	const uint32_t max_cat_ctl = 1 + (covar_max_nonnull_cat_ct / kBitsPerWord);
	if (bigstack_alloc_ul(max_cat_ctl, &cat_one_obs) ||
	    bigstack_alloc_ul(max_cat_ctl, &cat_two_or_more_obs) ||
	    bigstack_alloc_ui(covar_max_nonnull_cat_ct + 1, &cat_first_sample_uidxs)) {
	  goto glm_determine_covars_ret_NOMEM;
	}
      }
      do {
	prev_sample_ct = sample_ct;
	covar_uidx = 0;
	extra_cat_ct = 0;
	for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx, ++covar_uidx) {
	  next_set_unsafe_ck(covar_include, &covar_uidx);
	  const pheno_col_t* cur_covar_col = &(covar_cols[covar_uidx]);
	  if (cur_covar_col->type_code == kPhenoDtypeOther) {
	    continue;
	  }
	  if (cur_covar_col->type_code == kPhenoDtypeQt) {
	    const double* pheno_vals = cur_covar_col->data.qt;
	    next_set_unsafe_ck(cur_sample_include, &first_sample_uidx);
	    uint32_t sample_uidx = first_sample_uidx;
	    double common_pheno_val = pheno_vals[sample_uidx];
	    sample_uidx = next_set_unsafe(cur_sample_include, sample_uidx + 1);
	    const double second_pheno_val = pheno_vals[sample_uidx];
	    uint32_t sample_idx = 2;
	    uint32_t sample_uidx_remove;
	    if (second_pheno_val != common_pheno_val) {
	      sample_uidx_remove = sample_uidx;
	      sample_uidx = next_set_unsafe(cur_sample_include, sample_uidx + 1);
	      const double third_pheno_val = pheno_vals[sample_uidx];
	      if (third_pheno_val == second_pheno_val) {
		common_pheno_val = second_pheno_val;
		sample_uidx_remove = first_sample_uidx;
	      } else if (third_pheno_val != common_pheno_val) {
		continue;
	      }
	      sample_idx = 3;
	    } else {
	      sample_uidx_remove = 0xffffffffU;
	    }
	    for (; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
	      next_set_unsafe_ck(cur_sample_include, &sample_uidx);
	      if (pheno_vals[sample_uidx] != common_pheno_val) {
		if (sample_uidx_remove == 0xffffffffU) {
		  sample_uidx_remove = sample_uidx;
		} else {
		  break;
		}
	      }
	    }
	    if (sample_idx == sample_ct) {
	      if (sample_uidx_remove != 0xffffffffU) {
		if (--sample_ct == 2) {
		  goto glm_determine_covars_ret_SKIP;
		}
		clear_bit(sample_uidx_remove, cur_sample_include);
	      } else {
		// constant covariate, remove it
		clear_bit(covar_uidx, covar_include);
	      }
	    }
	  } else {
	    const uint32_t cur_nonnull_cat_ct = cur_covar_col->nonnull_category_ct;
	    const uint32_t cur_cat_ctl = 1 + (cur_nonnull_cat_ct / kBitsPerWord);
	    fill_ulong_zero(cur_cat_ctl, cat_one_obs);
	    fill_ulong_zero(cur_cat_ctl, cat_two_or_more_obs);
	    const uint32_t* pheno_vals = cur_covar_col->data.cat;
	    next_set_unsafe_ck(cur_sample_include, &first_sample_uidx);
	    uint32_t sample_uidx = first_sample_uidx;
	    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
	      next_set_unsafe_ck(cur_sample_include, &sample_uidx);
	      const uint32_t cur_cat_idx = pheno_vals[sample_uidx];
	      if (!is_set(cat_two_or_more_obs, cur_cat_idx)) {
		if (is_set(cat_one_obs, cur_cat_idx)) {
		  set_bit(cur_cat_idx, cat_two_or_more_obs);
		} else {
		  set_bit(cur_cat_idx, cat_one_obs);
		  cat_first_sample_uidxs[cur_cat_idx] = sample_uidx;
		}
	      }
	    }
	    for (uint32_t widx = 0; widx < cur_cat_ctl; ++widx) {
	      uintptr_t cur_word = cat_one_obs[widx] & (~cat_two_or_more_obs[widx]);
	      if (cur_word) {
		const uint32_t* cat_first_sample_uidxs_iter = &(cat_first_sample_uidxs[widx * kBitsPerWord]);
		do {
		  const uint32_t cat_idx_lowbits = CTZLU(cur_word);
		  --sample_ct;
		  clear_bit(cat_first_sample_uidxs_iter[cat_idx_lowbits], cur_sample_include);
		  cur_word &= cur_word - 1;
		} while (cur_word);
	      }
	    }
	    if (sample_ct < 3) {
	      goto glm_determine_covars_ret_SKIP;
	    }
	    uint32_t remaining_cat_ct = popcount_longs(cat_two_or_more_obs, cur_cat_ctl);
	    if (remaining_cat_ct <= 1) {
	      // now a constant covariate, remove it
	      clear_bit(covar_uidx, covar_include);
	    } else {
	      extra_cat_ct += remaining_cat_ct - 2;
	    }
	  }
	}
	covar_ct = popcount_longs(covar_include, raw_covar_ctl);
      } while (sample_ct < prev_sample_ct);
    } else {
      uintptr_t* cat_covar_wkspace;
      if (bigstack_alloc_ul(1 + (covar_max_nonnull_cat_ct / kBitsPerWordD2), &cat_covar_wkspace)) {
	goto glm_determine_covars_ret_NOMEM;
      }
      if (!is_sometimes_firth) {
	// todo: in firth-fallback case, automatically switch to always-Firth
	// if separated covariate is detected
	do {
	  prev_sample_ct = sample_ct;
	  covar_uidx = 0;
	  for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx, ++covar_uidx) {
	    next_set_unsafe_ck(covar_include, &covar_uidx);
	    if (check_for_and_handle_separated_covar(pheno_cc, covar_cols, raw_sample_ctl, covar_uidx, cur_sample_include, covar_include, &sample_ct, cat_covar_wkspace)) {
	      *separation_warning_ptr = 1;
	      goto glm_determine_covars_ret_SKIP;
	    }
	  }
	  covar_ct = popcount_longs(covar_include, raw_covar_ctl);
	} while (sample_ct < prev_sample_ct);
      }

      // now count extra categories
      covar_uidx = 0;
      extra_cat_ct = 0;
      for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx, ++covar_uidx) {
	next_set_unsafe_ck(covar_include, &covar_uidx);
	const pheno_col_t* cur_covar_col = &(covar_cols[covar_uidx]);
	if (cur_covar_col->type_code == kPhenoDtypeCat) {
	  const uint32_t remaining_cat_ct = identify_remaining_cats(cur_sample_include, cur_covar_col, sample_ct, cat_covar_wkspace);
	  if (remaining_cat_ct > 2) {
	    extra_cat_ct += remaining_cat_ct - 2;
	  }
	}
      }
    }
    *sample_ct_ptr = sample_ct;
    *covar_ct_ptr = covar_ct;
    *extra_cat_ct_ptr = extra_cat_ct;
  }
  while (0) {
  glm_determine_covars_ret_NOMEM:
    reterr = 1;
    break;
  glm_determine_covars_ret_SKIP:
    *sample_ct_ptr = 0;
  glm_determine_covars_ret_NOCOVAR:
    *covar_ct_ptr = 0;
    break;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

void collapse_parameter_subset(const uintptr_t* covar_include, const uintptr_t* raw_parameter_subset, uint32_t domdev_present, uint32_t raw_covar_ct, uint32_t covar_ct, uint32_t add_interactions, uintptr_t* new_parameter_subset, uint32_t* predictor_ct_ptr) {
  const uint32_t first_covar_pred_idx = 2 + domdev_present;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  uint32_t first_interaction_pred_read_idx = 0;
  uint32_t first_interaction_pred_write_idx = 0;
  if (add_interactions) {
    first_interaction_pred_read_idx = first_covar_pred_idx + raw_covar_ct;
    first_interaction_pred_write_idx = first_covar_pred_idx + covar_ct;
  }
  const uint32_t write_idx_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
  const uint32_t write_idx_ctl = BITCT_TO_WORDCT(write_idx_ct);
  fill_ulong_zero(write_idx_ctl, new_parameter_subset);
  // intercept, additive, domdev
  new_parameter_subset[0] = raw_parameter_subset[0] & (3 + 4 * domdev_present);
  uint32_t covar_uidx = 0;
  for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx, ++covar_uidx) {
    next_set_unsafe_ck(covar_include, &covar_uidx);
    if (is_set(raw_parameter_subset, first_covar_pred_idx + covar_uidx)) {
      set_bit(first_covar_pred_idx + covar_idx, new_parameter_subset);
    }
    if (add_interactions) {
      if (is_set(raw_parameter_subset, first_interaction_pred_read_idx + domdev_present_p1 * covar_uidx)) {
	set_bit(first_interaction_pred_write_idx + domdev_present_p1 * covar_idx, new_parameter_subset);
      }
      if (domdev_present) {
	if (is_set(raw_parameter_subset, first_interaction_pred_read_idx + 2 * covar_uidx + 1)) {
	  set_bit(first_interaction_pred_write_idx + 2 * covar_idx + 1, new_parameter_subset);
	}
      }
    }
  }
  *predictor_ct_ptr = popcount_longs(new_parameter_subset, write_idx_ctl);
}


ENUM_U31_DEF_START()
  kVifCorrCheckOk,
  kVifCorrCheckVifFail,
  kVifCorrCheckCorrFail
ENUM_U31_DEF_END(vif_corr_errcode_t);

typedef struct {
  vif_corr_errcode_t errcode;
  uint32_t covar_idx1; // for both correlation and VIF failure
  uint32_t covar_idx2; // for correlation failure only
} vif_corr_err_t;

boolerr_t glm_fill_and_test_covars(const uintptr_t* sample_include, const uintptr_t* covar_include, const pheno_col_t* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double vif_thresh, double max_corr, double* covars_smaj, double* covar_dotprod, double* covars_cmaj, char** cur_covar_names, vif_corr_err_t* vif_corr_check_result_ptr) {
  vif_corr_check_result_ptr->errcode = kVifCorrCheckOk;
  if (covar_ct == local_covar_ct) {
    return 0;
  }
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  uintptr_t* cat_covar_wkspace;
  matrix_invert_buf1_t* matrix_invert_buf1 = (matrix_invert_buf1_t*)bigstack_alloc(kMatrixInvertBuf1CheckedAlloc * new_nonlocal_covar_ct);
  double* inverse_corr_buf;
  double* dbl_2d_buf;
  if ((!matrix_invert_buf1) ||
      bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &inverse_corr_buf) ||
      bigstack_alloc_ul(1 + (covar_max_nonnull_cat_ct / kBitsPerWord), &cat_covar_wkspace) ||
      bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &dbl_2d_buf)) {
    return 1;
  }
  unsigned char* alloc_base = g_bigstack_base;
  unsigned char* new_covar_name_alloc = g_bigstack_end;
  const uint32_t first_sample_uidx = next_set_unsafe(sample_include, 0);
  uint32_t covar_read_uidx = 0;
  const double sample_ct_recip = 1.0 / ((double)((intptr_t)sample_ct));
  const double sample_ct_m1_d = (double)((intptr_t)(sample_ct - 1));
  char** cur_covar_names_iter = cur_covar_names;
  double* covar_write_iter = covars_cmaj;
  double* sum_iter = dbl_2d_buf;
  for (uintptr_t covar_read_idx = 0; covar_read_idx < covar_ct; ++covar_read_idx, ++covar_read_uidx) {
    next_set_unsafe_ck(covar_include, &covar_read_uidx);
    const pheno_col_t* cur_covar_col = &(covar_cols[covar_read_uidx]);
    const char* covar_name_base = &(covar_names[covar_read_uidx * max_covar_name_blen]);
    if (cur_covar_col->type_code == kPhenoDtypeOther) {
      // local covariate
      // const_cast
      *cur_covar_names_iter++ = (char*)((uintptr_t)covar_name_base);
    } else if (cur_covar_col->type_code == kPhenoDtypeQt) {
      // const_cast
      *cur_covar_names_iter++ = (char*)((uintptr_t)covar_name_base);
      const double* covar_vals = cur_covar_col->data.qt;
      uint32_t sample_uidx = first_sample_uidx;
      double covar_sum = 0.0;
      for (uintptr_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
	next_set_unsafe_ck(sample_include, &sample_uidx);
	const double cur_covar_val = covar_vals[sample_uidx];
	covar_sum += cur_covar_val;
	*covar_write_iter++ = cur_covar_val;
      }
      *sum_iter++ = covar_sum;
    } else {
      const uint32_t remaining_cat_ct = identify_remaining_cats(sample_include, cur_covar_col, sample_ct, cat_covar_wkspace);
      assert(remaining_cat_ct >= 2);
      const uint32_t* covar_vals = cur_covar_col->data.cat;
      char** cur_category_names = cur_covar_col->category_names;
      const uint32_t covar_name_base_slen = strlen(covar_name_base);
      uint32_t cat_uidx = 1;
      // this is equivalent to "--split-cat-pheno omit-last covar-01"
      for (uint32_t cat_idx = 1; cat_idx < remaining_cat_ct; ++cat_idx, ++cat_uidx) {
	next_set_unsafe_ck(cat_covar_wkspace, &cat_uidx);

	const char* catname = cur_category_names[cat_uidx];
	const uint32_t catname_slen = strlen(catname);
	new_covar_name_alloc -= covar_name_base_slen + catname_slen + 2;
	if (new_covar_name_alloc < alloc_base) {
	  return 1;
	}
	*cur_covar_names_iter++ = (char*)new_covar_name_alloc;
	char* new_covar_name_write = memcpyax(new_covar_name_alloc, covar_name_base, covar_name_base_slen, '=');
	memcpyx(new_covar_name_write, catname, catname_slen, '\0');
	
	uint32_t sample_uidx = first_sample_uidx;
	uint32_t cur_cat_obs_ct = 0;
	for (uintptr_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
	  next_set_unsafe_ck(sample_include, &sample_uidx);
	  const uint32_t cur_sample_is_in_cat = (covar_vals[sample_uidx] == cat_uidx);
	  cur_cat_obs_ct += cur_sample_is_in_cat;
	  *covar_write_iter++ = (double)((int32_t)cur_sample_is_in_cat);
	}
	*sum_iter++ = (double)((int32_t)cur_cat_obs_ct);
      }
    }
  }
  bigstack_end_set(new_covar_name_alloc);
  assert(covar_write_iter == &(covars_cmaj[new_nonlocal_covar_ct * sample_ct]));
  transpose_copy(covars_cmaj, new_nonlocal_covar_ct, sample_ct, covars_smaj);
  row_major_matrix_multiply(covars_cmaj, covars_smaj, new_nonlocal_covar_ct, new_nonlocal_covar_ct, sample_ct, covar_dotprod);
  // we have dot products, now determine
  //   (dotprod - sum(a)mean(b)) / (N-1)
  // to get small-sample covariance
  const uintptr_t new_nonlocal_covar_ct_p1 = new_nonlocal_covar_ct + 1;
  const double sample_ct_m1_recip = 1.0 / sample_ct_m1_d;
  double* covar_dotprod_iter = covar_dotprod;
  double* sample_covariance_iter = inverse_corr_buf;
  for (uintptr_t covar_idx1 = 0; covar_idx1 < new_nonlocal_covar_ct; ++covar_idx1) {
    const double covar1_mean_adj = dbl_2d_buf[covar_idx1] * sample_ct_recip;
    for (uintptr_t covar_idx2 = 0; covar_idx2 < new_nonlocal_covar_ct; ++covar_idx2) {
      *sample_covariance_iter++ = ((*covar_dotprod_iter++) - covar1_mean_adj * dbl_2d_buf[covar_idx2]) * sample_ct_m1_recip;
    }
  }
  // now use dbl_2d_buf to store inverse-sqrts, to get to correlation matrix
  for (uintptr_t covar_idx = 0; covar_idx < new_nonlocal_covar_ct; ++covar_idx) {
    dbl_2d_buf[covar_idx] = 1.0 / sqrt(inverse_corr_buf[covar_idx * new_nonlocal_covar_ct_p1]);
  }
  for (uintptr_t covar_idx1 = 1; covar_idx1 < new_nonlocal_covar_ct; ++covar_idx1) {
    const double inverse_stdev1 = dbl_2d_buf[covar_idx1];
    double* corr_row_iter = &(inverse_corr_buf[covar_idx1 * new_nonlocal_covar_ct]);
    double* corr_col_start = &(inverse_corr_buf[covar_idx1]);
    const double* inverse_stdev2_iter = dbl_2d_buf;
    for (uintptr_t covar_idx2 = 0; covar_idx2 < covar_idx1; ++covar_idx2) {
      double* corr_col_entry_ptr = &(corr_col_start[covar_idx2 * new_nonlocal_covar_ct]);
      const double cur_corr = (*corr_col_entry_ptr) * inverse_stdev1 * (*inverse_stdev2_iter++);
      if (cur_corr > max_corr) {
	vif_corr_check_result_ptr->errcode = kVifCorrCheckCorrFail;
	// may as well put smaller index first
	vif_corr_check_result_ptr->covar_idx1 = covar_idx2;
	vif_corr_check_result_ptr->covar_idx2 = covar_idx1;
	// bigstack_reset unnecessary here, we'll reset the stack anyway before
	// exiting or starting the next covariate
	return 0;
      }
      *corr_col_entry_ptr = cur_corr;
      *corr_row_iter++ = cur_corr;
    }
  }
  for (uintptr_t covar_idx = 0; covar_idx < new_nonlocal_covar_ct; ++covar_idx) {
    inverse_corr_buf[covar_idx * new_nonlocal_covar_ct_p1] = 1.0;
  }
  if (invert_matrix_checked(new_nonlocal_covar_ct, inverse_corr_buf, matrix_invert_buf1, dbl_2d_buf)) {
    vif_corr_check_result_ptr->errcode = kVifCorrCheckVifFail;
    vif_corr_check_result_ptr->covar_idx1 = 0xffffffffU;
    return 0;
  }
  // VIFs = diagonal elements of inverse correlation matrix
  for (uintptr_t covar_idx = 0; covar_idx < new_nonlocal_covar_ct; ++covar_idx) {
    if (inverse_corr_buf[covar_idx * new_nonlocal_covar_ct_p1] > vif_thresh) {
      vif_corr_check_result_ptr->errcode = kVifCorrCheckVifFail;
      vif_corr_check_result_ptr->covar_idx1 = covar_idx + local_covar_ct;
      return 0;
    }
  }
  bigstack_reset(matrix_invert_buf1);
  return 0;
}

boolerr_t glm_alloc_fill_and_test_pheno_covars_qt(const uintptr_t* sample_include, const double* pheno_qt, const uintptr_t* covar_include, const pheno_col_t* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double vif_thresh, double max_corr, double** pheno_d_ptr, double** covars_cmaj_d_ptr, char*** cur_covar_names_ptr, vif_corr_err_t* vif_corr_check_result_ptr) {
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  double* covar_dotprod;
  double* covars_smaj;
  if (bigstack_alloc_d(sample_ct, pheno_d_ptr) ||
      bigstack_alloc_cp(new_covar_ct, cur_covar_names_ptr) ||
      bigstack_alloc_d(new_nonlocal_covar_ct * sample_ct, covars_cmaj_d_ptr) ||
      bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &covar_dotprod) ||
      bigstack_alloc_d(new_nonlocal_covar_ct * sample_ct, &covars_smaj)) {
    return 1;
  }
  double* pheno_d_iter = *pheno_d_ptr;
  uint32_t sample_uidx = 0;
  for (uintptr_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    next_set_unsafe_ck(sample_include, &sample_uidx);
    *pheno_d_iter++ = pheno_qt[sample_uidx];
  }
  if (glm_fill_and_test_covars(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, vif_thresh, max_corr, covars_smaj, covar_dotprod, *covars_cmaj_d_ptr, *cur_covar_names_ptr, vif_corr_check_result_ptr)) {
    return 1;
  }
  bigstack_reset(covar_dotprod);
  return 0;
}

boolerr_t glm_alloc_fill_and_test_pheno_covars_cc(const uintptr_t* sample_include, const uintptr_t* pheno_cc, const uintptr_t* covar_include, const pheno_col_t* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double vif_thresh, double max_corr, uintptr_t** pheno_cc_collapsed_ptr, float** pheno_f_ptr, float** covars_cmaj_f_ptr, char*** cur_covar_names_ptr, vif_corr_err_t* vif_corr_check_result_ptr) {
  const uintptr_t sample_cta4 = round_up_pow2(sample_ct, 4);
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  double* covars_cmaj_d;
  double* covars_smaj_d;
  double* covar_dotprod;
  if (bigstack_alloc_ul(BITCT_TO_WORDCT(sample_ct), pheno_cc_collapsed_ptr) ||
      bigstack_alloc_f(sample_cta4, pheno_f_ptr) ||
      bigstack_alloc_f(new_nonlocal_covar_ct * sample_cta4, covars_cmaj_f_ptr) ||
      bigstack_alloc_cp(new_covar_ct, cur_covar_names_ptr) ||
      bigstack_alloc_d(new_nonlocal_covar_ct * sample_ct, &covars_cmaj_d) ||
      bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &covar_dotprod) ||
      bigstack_alloc_d(new_nonlocal_covar_ct * sample_ct, &covars_smaj_d)) {
    return 1;
  }
  uintptr_t* pheno_cc_collapsed = *pheno_cc_collapsed_ptr;
  copy_bitarr_subset(pheno_cc, sample_include, sample_ct, pheno_cc_collapsed);
  float* pheno_f_iter = *pheno_f_ptr;
  for (uintptr_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
    *pheno_f_iter++ = (float)((int32_t)is_set(pheno_cc_collapsed, sample_idx));
  }
  const uint32_t sample_rem4 = sample_cta4 - sample_ct;
  fill_float_zero(sample_rem4, pheno_f_iter);
  if (glm_fill_and_test_covars(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, vif_thresh, max_corr, covars_smaj_d, covar_dotprod, covars_cmaj_d, *cur_covar_names_ptr, vif_corr_check_result_ptr)) {
    return 1;
  }
  double* covar_read_iter = covars_cmaj_d;
  float* covar_write_iter = *covars_cmaj_f_ptr;
  for (uintptr_t covar_idx = 0; covar_idx < new_nonlocal_covar_ct; ++covar_idx) {
    for (uintptr_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      *covar_write_iter++ = (float)(*covar_read_iter++);
    }
    fill_float_zero(sample_rem4, covar_write_iter);
    covar_write_iter = &(covar_write_iter[sample_rem4]);
  }
  bigstack_reset(covars_cmaj_d);
  return 0;
}

static const float kSmallFloats[4] = {0.0f, 1.0f, 2.0f, 3.0f};
// static const float kSmallFloats[4] = {0.0f, 1.0f, 2.0f, -9.0f};

void genoarr_to_floats(const uintptr_t* genoarr, uint32_t sample_ct, float* floatbuf) {
  assert(sample_ct);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t widx = 0;
  uint32_t subgroup_len = kBitsPerWordD2;
  float* floatbuf_iter = floatbuf;
  while (1) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
	return;
      }
      subgroup_len = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
      const uintptr_t cur_geno = geno_word & 3;
      // *floatbuf_iter++ = (float)((int32_t)cur_geno);
      *floatbuf_iter++ = kSmallFloats[cur_geno];
      geno_word >>= 2;
    }
    ++widx;
  }
}

uint32_t genoarr_to_floats_remove_missing(const uintptr_t* genoarr, uint32_t sample_ct, float* floatbuf) {
  assert(sample_ct);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t widx = 0;
  uint32_t subgroup_len = kBitsPerWordD2;
  float* floatbuf_iter = floatbuf;
  while (1) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
	return (uintptr_t)(floatbuf_iter - floatbuf);
      }
      subgroup_len = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
      const uintptr_t cur_geno = geno_word & 3;
      if (cur_geno < 3) {
	// *floatbuf_iter++ = (float)((int32_t)cur_geno);
	*floatbuf_iter++ = kSmallFloats[cur_geno];
      }
      geno_word >>= 2;
    }
    ++widx;
  }
}

// #####
// The following code is based on the winning submissino of Pascal Pons in the
// "GWASSpeedup" contest run in April 2013 by Babbage Analytics & Innovation
// and TopCoder, who have donated the results to be used in PLINK.  The contest
// was designed by Po-Ru Loh; subsequent analysis and code preparation were
// performed by Andrew Hill, Ragu Bharadwaj, and Scott Jelinsky.  A manuscript
// is in preparation by these authors and Iain Kilty, Kevin Boudreau, Karim
// Lakhani and Eva Guinan.
// #####

#ifdef __LP64__
// fmath_exp_ps is a C port of Shigeo Mitsunari's fast math library function
// posted at https://github.com/herumi/fmath .  License is
// http://opensource.org/licenses/BSD-3-Clause .
// (I tried porting fmath_log_ps, but it turns out that Firth regression needs
// double-precision log accuracy; logf() actually interferes with convergence.)

// programmatically generated by:
// typedef union {
//   float f4;
//   uint32_t u4;
// } __uni4;
//
// __uni4 u4;
// int32_t ii;
// for (ii = 0; ii < 1024; ii++) {
//   u4.f4 = pow(2.0f, ((float)ii) / 1024.0);
//   printf("0x%08x", u4.u4 & 0x7fffff);
//   if (ii % 4 != 3) {
//     printf(", ");
//   } else {
//     printf(",\n");
//   }
// }
const uint32_t float_exp_lookup_int[] __attribute__((aligned(16))) = {
0x00000000, 0x00001630, 0x00002c64, 0x0000429c,
0x000058d8, 0x00006f17, 0x0000855b, 0x00009ba2,
0x0000b1ed, 0x0000c83c, 0x0000de8f, 0x0000f4e6,
0x00010b41, 0x0001219f, 0x00013802, 0x00014e68,
0x000164d2, 0x00017b40, 0x000191b2, 0x0001a828,
0x0001bea1, 0x0001d51f, 0x0001eba1, 0x00020226,
0x000218af, 0x00022f3c, 0x000245ce, 0x00025c63,
0x000272fc, 0x00028998, 0x0002a039, 0x0002b6de,
0x0002cd87, 0x0002e433, 0x0002fae4, 0x00031198,
0x00032850, 0x00033f0d, 0x000355cd, 0x00036c91,
0x00038359, 0x00039a25, 0x0003b0f5, 0x0003c7c9,
0x0003dea1, 0x0003f57d, 0x00040c5d, 0x00042341,
0x00043a29, 0x00045115, 0x00046804, 0x00047ef8,
0x000495f0, 0x0004aceb, 0x0004c3eb, 0x0004daef,
0x0004f1f6, 0x00050902, 0x00052012, 0x00053725,
0x00054e3d, 0x00056558, 0x00057c78, 0x0005939c,
0x0005aac3, 0x0005c1ef, 0x0005d91f, 0x0005f052,
0x0006078a, 0x00061ec6, 0x00063606, 0x00064d4a,
0x00066491, 0x00067bdd, 0x0006932d, 0x0006aa81,
0x0006c1d9, 0x0006d935, 0x0006f095, 0x000707f9,
0x00071f62, 0x000736ce, 0x00074e3e, 0x000765b3,
0x00077d2b, 0x000794a8, 0x0007ac28, 0x0007c3ad,
0x0007db35, 0x0007f2c2, 0x00080a53, 0x000821e8,
0x00083981, 0x0008511e, 0x000868c0, 0x00088065,
0x0008980f, 0x0008afbc, 0x0008c76e, 0x0008df23,
0x0008f6dd, 0x00090e9b, 0x0009265d, 0x00093e24,
0x000955ee, 0x00096dbc, 0x0009858f, 0x00099d66,
0x0009b541, 0x0009cd20, 0x0009e503, 0x0009fcea,
0x000a14d5, 0x000a2cc5, 0x000a44b9, 0x000a5cb1,
0x000a74ad, 0x000a8cad, 0x000aa4b1, 0x000abcba,
0x000ad4c6, 0x000aecd7, 0x000b04ec, 0x000b1d05,
0x000b3523, 0x000b4d44, 0x000b656a, 0x000b7d94,
0x000b95c2, 0x000badf4, 0x000bc62b, 0x000bde65,
0x000bf6a4, 0x000c0ee7, 0x000c272f, 0x000c3f7a,
0x000c57ca, 0x000c701e, 0x000c8876, 0x000ca0d2,
0x000cb933, 0x000cd198, 0x000cea01, 0x000d026e,
0x000d1adf, 0x000d3355, 0x000d4bcf, 0x000d644d,
0x000d7cd0, 0x000d9556, 0x000dade1, 0x000dc671,
0x000ddf04, 0x000df79c, 0x000e1038, 0x000e28d8,
0x000e417d, 0x000e5a25, 0x000e72d3, 0x000e8b84,
0x000ea43a, 0x000ebcf3, 0x000ed5b2, 0x000eee74,
0x000f073b, 0x000f2006, 0x000f38d5, 0x000f51a9,
0x000f6a81, 0x000f835d, 0x000f9c3e, 0x000fb523,
0x000fce0c, 0x000fe6fa, 0x000fffec, 0x001018e2,
0x001031dc, 0x00104adb, 0x001063de, 0x00107ce6,
0x001095f2, 0x0010af02, 0x0010c816, 0x0010e12f,
0x0010fa4d, 0x0011136e, 0x00112c94, 0x001145be,
0x00115eed, 0x00117820, 0x00119158, 0x0011aa93,
0x0011c3d3, 0x0011dd18, 0x0011f661, 0x00120fae,
0x00122900, 0x00124256, 0x00125bb0, 0x0012750f,
0x00128e72, 0x0012a7da, 0x0012c146, 0x0012dab7,
0x0012f42c, 0x00130da5, 0x00132723, 0x001340a5,
0x00135a2b, 0x001373b6, 0x00138d46, 0x0013a6d9,
0x0013c072, 0x0013da0e, 0x0013f3af, 0x00140d55,
0x001426ff, 0x001440ae, 0x00145a60, 0x00147418,
0x00148dd4, 0x0014a794, 0x0014c159, 0x0014db22,
0x0014f4f0, 0x00150ec2, 0x00152898, 0x00154274,
0x00155c53, 0x00157637, 0x00159020, 0x0015aa0d,
0x0015c3ff, 0x0015ddf5, 0x0015f7ef, 0x001611ee,
0x00162bf2, 0x001645fa, 0x00166006, 0x00167a18,
0x0016942d, 0x0016ae47, 0x0016c866, 0x0016e289,
0x0016fcb1, 0x001716dd, 0x0017310e, 0x00174b43,
0x0017657d, 0x00177fbc, 0x001799ff, 0x0017b446,
0x0017ce92, 0x0017e8e3, 0x00180338, 0x00181d92,
0x001837f0, 0x00185253, 0x00186cbb, 0x00188727,
0x0018a197, 0x0018bc0d, 0x0018d686, 0x0018f105,
0x00190b88, 0x0019260f, 0x0019409c, 0x00195b2c,
0x001975c2, 0x0019905c, 0x0019aafa, 0x0019c59e,
0x0019e046, 0x0019faf2, 0x001a15a3, 0x001a3059,
0x001a4b13, 0x001a65d2, 0x001a8096, 0x001a9b5e,
0x001ab62b, 0x001ad0fd, 0x001aebd3, 0x001b06ae,
0x001b218d, 0x001b3c71, 0x001b575a, 0x001b7248,
0x001b8d3a, 0x001ba831, 0x001bc32c, 0x001bde2c,
0x001bf931, 0x001c143b, 0x001c2f49, 0x001c4a5c,
0x001c6573, 0x001c8090, 0x001c9bb1, 0x001cb6d6,
0x001cd201, 0x001ced30, 0x001d0864, 0x001d239c,
0x001d3eda, 0x001d5a1c, 0x001d7562, 0x001d90ae,
0x001dabfe, 0x001dc753, 0x001de2ad, 0x001dfe0b,
0x001e196e, 0x001e34d6, 0x001e5043, 0x001e6bb4,
0x001e872a, 0x001ea2a5, 0x001ebe25, 0x001ed9a9,
0x001ef532, 0x001f10c0, 0x001f2c53, 0x001f47eb,
0x001f6387, 0x001f7f28, 0x001f9ace, 0x001fb679,
0x001fd228, 0x001feddc, 0x00200996, 0x00202553,
0x00204116, 0x00205cde, 0x002078aa, 0x0020947b,
0x0020b051, 0x0020cc2c, 0x0020e80b, 0x002103f0,
0x00211fd9, 0x00213bc7, 0x002157ba, 0x002173b2,
0x00218faf, 0x0021abb0, 0x0021c7b7, 0x0021e3c2,
0x0021ffd2, 0x00221be7, 0x00223801, 0x0022541f,
0x00227043, 0x00228c6b, 0x0022a899, 0x0022c4cb,
0x0022e102, 0x0022fd3e, 0x0023197f, 0x002335c5,
0x0023520f, 0x00236e5f, 0x00238ab3, 0x0023a70d,
0x0023c36b, 0x0023dfce, 0x0023fc37, 0x002418a4,
0x00243516, 0x0024518d, 0x00246e08, 0x00248a89,
0x0024a70f, 0x0024c39a, 0x0024e029, 0x0024fcbe,
0x00251958, 0x002535f6, 0x00255299, 0x00256f42,
0x00258bef, 0x0025a8a2, 0x0025c559, 0x0025e215,
0x0025fed7, 0x00261b9d, 0x00263868, 0x00265538,
0x0026720e, 0x00268ee8, 0x0026abc7, 0x0026c8ac,
0x0026e595, 0x00270283, 0x00271f76, 0x00273c6f,
0x0027596c, 0x0027766e, 0x00279376, 0x0027b082,
0x0027cd94, 0x0027eaaa, 0x002807c6, 0x002824e6,
0x0028420c, 0x00285f37, 0x00287c66, 0x0028999b,
0x0028b6d5, 0x0028d414, 0x0028f158, 0x00290ea1,
0x00292bef, 0x00294942, 0x0029669b, 0x002983f8,
0x0029a15b, 0x0029bec2, 0x0029dc2f, 0x0029f9a1,
0x002a1718, 0x002a3494, 0x002a5215, 0x002a6f9b,
0x002a8d26, 0x002aaab7, 0x002ac84c, 0x002ae5e7,
0x002b0387, 0x002b212c, 0x002b3ed6, 0x002b5c85,
0x002b7a3a, 0x002b97f3, 0x002bb5b2, 0x002bd376,
0x002bf13f, 0x002c0f0d, 0x002c2ce0, 0x002c4ab9,
0x002c6897, 0x002c867a, 0x002ca462, 0x002cc24f,
0x002ce041, 0x002cfe39, 0x002d1c36, 0x002d3a38,
0x002d583f, 0x002d764b, 0x002d945d, 0x002db274,
0x002dd090, 0x002deeb1, 0x002e0cd8, 0x002e2b03,
0x002e4934, 0x002e676b, 0x002e85a6, 0x002ea3e7,
0x002ec22d, 0x002ee078, 0x002efec8, 0x002f1d1e,
0x002f3b79, 0x002f59d9, 0x002f783e, 0x002f96a9,
0x002fb519, 0x002fd38e, 0x002ff209, 0x00301089,
0x00302f0e, 0x00304d98, 0x00306c28, 0x00308abd,
0x0030a957, 0x0030c7f7, 0x0030e69c, 0x00310546,
0x003123f6, 0x003142aa, 0x00316165, 0x00318024,
0x00319ee9, 0x0031bdb3, 0x0031dc83, 0x0031fb57,
0x00321a32, 0x00323911, 0x003257f6, 0x003276e0,
0x003295d0, 0x0032b4c5, 0x0032d3bf, 0x0032f2bf,
0x003311c4, 0x003330cf, 0x00334fde, 0x00336ef4,
0x00338e0e, 0x0033ad2e, 0x0033cc54, 0x0033eb7e,
0x00340aaf, 0x003429e4, 0x0034491f, 0x00346860,
0x003487a6, 0x0034a6f1, 0x0034c642, 0x0034e598,
0x003504f3, 0x00352454, 0x003543bb, 0x00356327,
0x00358298, 0x0035a20f, 0x0035c18b, 0x0035e10d,
0x00360094, 0x00362020, 0x00363fb2, 0x00365f4a,
0x00367ee7, 0x00369e89, 0x0036be31, 0x0036dddf,
0x0036fd92, 0x00371d4a, 0x00373d08, 0x00375ccc,
0x00377c95, 0x00379c63, 0x0037bc37, 0x0037dc11,
0x0037fbf0, 0x00381bd4, 0x00383bbe, 0x00385bae,
0x00387ba3, 0x00389b9e, 0x0038bb9e, 0x0038dba4,
0x0038fbaf, 0x00391bc0, 0x00393bd7, 0x00395bf3,
0x00397c14, 0x00399c3b, 0x0039bc68, 0x0039dc9a,
0x0039fcd2, 0x003a1d10, 0x003a3d53, 0x003a5d9b,
0x003a7dea, 0x003a9e3e, 0x003abe97, 0x003adef6,
0x003aff5b, 0x003b1fc5, 0x003b4035, 0x003b60aa,
0x003b8126, 0x003ba1a6, 0x003bc22d, 0x003be2b9,
0x003c034a, 0x003c23e2, 0x003c447f, 0x003c6521,
0x003c85ca, 0x003ca678, 0x003cc72b, 0x003ce7e5,
0x003d08a4, 0x003d2968, 0x003d4a33, 0x003d6b03,
0x003d8bd8, 0x003dacb4, 0x003dcd95, 0x003dee7c,
0x003e0f68, 0x003e305a, 0x003e5152, 0x003e7250,
0x003e9353, 0x003eb45c, 0x003ed56b, 0x003ef67f,
0x003f179a, 0x003f38ba, 0x003f59df, 0x003f7b0b,
0x003f9c3c, 0x003fbd73, 0x003fdeb0, 0x003ffff2,
0x0040213b, 0x00404289, 0x004063dc, 0x00408536,
0x0040a695, 0x0040c7fb, 0x0040e966, 0x00410ad6,
0x00412c4d, 0x00414dc9, 0x00416f4b, 0x004190d3,
0x0041b261, 0x0041d3f5, 0x0041f58e, 0x0042172d,
0x004238d2, 0x00425a7d, 0x00427c2e, 0x00429de4,
0x0042bfa1, 0x0042e163, 0x0043032b, 0x004324f9,
0x004346cd, 0x004368a7, 0x00438a86, 0x0043ac6b,
0x0043ce57, 0x0043f048, 0x0044123f, 0x0044343c,
0x0044563f, 0x00447848, 0x00449a56, 0x0044bc6b,
0x0044de85, 0x004500a5, 0x004522cc, 0x004544f8,
0x0045672a, 0x00458962, 0x0045aba0, 0x0045cde4,
0x0045f02e, 0x0046127e, 0x004634d3, 0x0046572f,
0x00467991, 0x00469bf8, 0x0046be66, 0x0046e0d9,
0x00470353, 0x004725d2, 0x00474858, 0x00476ae3,
0x00478d75, 0x0047b00c, 0x0047d2aa, 0x0047f54d,
0x004817f7, 0x00483aa6, 0x00485d5b, 0x00488017,
0x0048a2d8, 0x0048c5a0, 0x0048e86d, 0x00490b41,
0x00492e1b, 0x004950fa, 0x004973e0, 0x004996cc,
0x0049b9be, 0x0049dcb5, 0x0049ffb3, 0x004a22b7,
0x004a45c1, 0x004a68d1, 0x004a8be8, 0x004aaf04,
0x004ad226, 0x004af54f, 0x004b187d, 0x004b3bb2,
0x004b5eed, 0x004b822e, 0x004ba575, 0x004bc8c2,
0x004bec15, 0x004c0f6e, 0x004c32ce, 0x004c5633,
0x004c799f, 0x004c9d11, 0x004cc089, 0x004ce407,
0x004d078c, 0x004d2b16, 0x004d4ea7, 0x004d723d,
0x004d95da, 0x004db97e, 0x004ddd27, 0x004e00d6,
0x004e248c, 0x004e4848, 0x004e6c0a, 0x004e8fd2,
0x004eb3a1, 0x004ed775, 0x004efb50, 0x004f1f31,
0x004f4319, 0x004f6706, 0x004f8afa, 0x004faef4,
0x004fd2f4, 0x004ff6fb, 0x00501b08, 0x00503f1b,
0x00506334, 0x00508753, 0x0050ab79, 0x0050cfa5,
0x0050f3d7, 0x00511810, 0x00513c4f, 0x00516094,
0x005184df, 0x0051a931, 0x0051cd89, 0x0051f1e7,
0x0052164c, 0x00523ab7, 0x00525f28, 0x005283a0,
0x0052a81e, 0x0052cca2, 0x0052f12c, 0x005315bd,
0x00533a54, 0x00535ef2, 0x00538396, 0x0053a840,
0x0053ccf1, 0x0053f1a8, 0x00541665, 0x00543b29,
0x00545ff3, 0x005484c3, 0x0054a99a, 0x0054ce77,
0x0054f35b, 0x00551845, 0x00553d35, 0x0055622c,
0x00558729, 0x0055ac2d, 0x0055d137, 0x0055f647,
0x00561b5e, 0x0056407b, 0x0056659f, 0x00568ac9,
0x0056affa, 0x0056d531, 0x0056fa6e, 0x00571fb2,
0x005744fd, 0x00576a4e, 0x00578fa5, 0x0057b503,
0x0057da67, 0x0057ffd2, 0x00582543, 0x00584abb,
0x00587039, 0x005895be, 0x0058bb49, 0x0058e0db,
0x00590673, 0x00592c12, 0x005951b8, 0x00597763,
0x00599d16, 0x0059c2cf, 0x0059e88e, 0x005a0e54,
0x005a3421, 0x005a59f4, 0x005a7fcd, 0x005aa5ae,
0x005acb94, 0x005af182, 0x005b1776, 0x005b3d70,
0x005b6371, 0x005b8979, 0x005baf87, 0x005bd59c,
0x005bfbb8, 0x005c21da, 0x005c4802, 0x005c6e32,
0x005c9468, 0x005cbaa4, 0x005ce0e7, 0x005d0731,
0x005d2d82, 0x005d53d9, 0x005d7a36, 0x005da09b,
0x005dc706, 0x005ded77, 0x005e13f0, 0x005e3a6f,
0x005e60f5, 0x005e8781, 0x005eae14, 0x005ed4ae,
0x005efb4e, 0x005f21f5, 0x005f48a3, 0x005f6f58,
0x005f9613, 0x005fbcd5, 0x005fe39e, 0x00600a6d,
0x00603143, 0x00605820, 0x00607f03, 0x0060a5ee,
0x0060ccdf, 0x0060f3d7, 0x00611ad5, 0x006141db,
0x006168e7, 0x00618ffa, 0x0061b713, 0x0061de34,
0x0062055b, 0x00622c89, 0x006253be, 0x00627af9,
0x0062a23c, 0x0062c985, 0x0062f0d5, 0x0063182c,
0x00633f89, 0x006366ee, 0x00638e59, 0x0063b5cb,
0x0063dd44, 0x006404c4, 0x00642c4b, 0x006453d8,
0x00647b6d, 0x0064a308, 0x0064caaa, 0x0064f253,
0x00651a03, 0x006541b9, 0x00656977, 0x0065913c,
0x0065b907, 0x0065e0d9, 0x006608b2, 0x00663092,
0x00665879, 0x00668067, 0x0066a85c, 0x0066d058,
0x0066f85b, 0x00672064, 0x00674875, 0x0067708c,
0x006798ab, 0x0067c0d0, 0x0067e8fd, 0x00681130,
0x0068396a, 0x006861ac, 0x006889f4, 0x0068b243,
0x0068da99, 0x006902f7, 0x00692b5b, 0x006953c6,
0x00697c38, 0x0069a4b1, 0x0069cd32, 0x0069f5b9,
0x006a1e47, 0x006a46dd, 0x006a6f79, 0x006a981c,
0x006ac0c7, 0x006ae978, 0x006b1231, 0x006b3af1,
0x006b63b7, 0x006b8c85, 0x006bb55a, 0x006bde36,
0x006c0719, 0x006c3003, 0x006c58f4, 0x006c81ec,
0x006caaec, 0x006cd3f2, 0x006cfd00, 0x006d2614,
0x006d4f30, 0x006d7853, 0x006da17d, 0x006dcaae,
0x006df3e7, 0x006e1d26, 0x006e466d, 0x006e6fbb,
0x006e9910, 0x006ec26c, 0x006eebcf, 0x006f1539,
0x006f3eab, 0x006f6824, 0x006f91a4, 0x006fbb2b,
0x006fe4ba, 0x00700e4f, 0x007037ec, 0x00706190,
0x00708b3b, 0x0070b4ee, 0x0070dea8, 0x00710868,
0x00713231, 0x00715c00, 0x007185d7, 0x0071afb5,
0x0071d99a, 0x00720386, 0x00722d7a, 0x00725775,
0x00728177, 0x0072ab81, 0x0072d592, 0x0072ffaa,
0x007329c9, 0x007353f0, 0x00737e1e, 0x0073a853,
0x0073d290, 0x0073fcd4, 0x0074271f, 0x00745172,
0x00747bcc, 0x0074a62d, 0x0074d096, 0x0074fb06,
0x0075257d, 0x00754ffc, 0x00757a82, 0x0075a50f,
0x0075cfa4, 0x0075fa40, 0x007624e4, 0x00764f8f,
0x00767a41, 0x0076a4fb, 0x0076cfbc, 0x0076fa85,
0x00772555, 0x0077502d, 0x00777b0b, 0x0077a5f2,
0x0077d0df, 0x0077fbd5, 0x007826d1, 0x007851d5,
0x00787ce1, 0x0078a7f4, 0x0078d30e, 0x0078fe30,
0x0079295a, 0x0079548b, 0x00797fc3, 0x0079ab03,
0x0079d64a, 0x007a0199, 0x007a2cf0, 0x007a584d,
0x007a83b3, 0x007aaf20, 0x007ada94, 0x007b0610,
0x007b3194, 0x007b5d1f, 0x007b88b2, 0x007bb44c,
0x007bdfed, 0x007c0b97, 0x007c3748, 0x007c6300,
0x007c8ec0, 0x007cba88, 0x007ce657, 0x007d122e,
0x007d3e0c, 0x007d69f2, 0x007d95e0, 0x007dc1d5,
0x007dedd2, 0x007e19d6, 0x007e45e2, 0x007e71f6,
0x007e9e11, 0x007eca34, 0x007ef65f, 0x007f2291,
0x007f4ecb, 0x007f7b0d, 0x007fa756, 0x007fd3a7
};

const float* const float_exp_lookup = (const float*)float_exp_lookup_int;

static inline __m128 fmath_exp_ps(__m128 xx) {
  const __m128i mask7ff = {0x7fffffff7fffffffLLU, 0x7fffffff7fffffffLLU};

  // 88
  const __m128i max_x = {0x42b0000042b00000LLU, 0x42b0000042b00000LLU};
  // -88
  // more sensible 0xc2b00000... not used here due to "narrowing conversion"
  // warning
  const __m128i min_x = {-0x3d4fffff3d500000LL, -0x3d4fffff3d500000LL};
  // 2^10 / log(2)
  const __m128i const_aa = {0x44b8aa3b44b8aa3bLLU, 0x44b8aa3b44b8aa3bLLU};
  // log(2) / 2^10
  const __m128i const_bb = {0x3a3172183a317218LLU, 0x3a3172183a317218LLU};

  const __m128i f1 = {0x3f8000003f800000LLU, 0x3f8000003f800000LLU};
  const __m128i mask_s = {0x3ff000003ffLLU, 0x3ff000003ffLLU};
  const __m128i i127s = {0x1fc000001fc00LLU, 0x1fc000001fc00LLU};
  const __m128i limit = _mm_castps_si128(_mm_and_ps(xx, (__m128)mask7ff));
  const int32_t over = _mm_movemask_epi8(_mm_cmpgt_epi32(limit, max_x));
  if (over) {
    xx = _mm_min_ps(xx, (__m128)max_x);
    xx = _mm_max_ps(xx, (__m128)min_x);
  }
  const __m128i rr = _mm_cvtps_epi32(_mm_mul_ps(xx, (__m128)const_aa));
  __m128 tt = _mm_sub_ps(xx, _mm_mul_ps(_mm_cvtepi32_ps(rr), (__m128)const_bb));
  tt = _mm_add_ps(tt, (__m128)f1);
  const __m128i v4 = _mm_and_si128(rr, mask_s);
  __m128i u4 = _mm_add_epi32(rr, i127s);
  u4 = _mm_srli_epi32(u4, 10);
  u4 = _mm_slli_epi32(u4, 23);
  const uint32_t v0 = _mm_cvtsi128_si32(v4);
  // uint32_t v1 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(2)));
  // uint32_t v2 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(4)));
  // uint32_t v3 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(6)));
  // make this work with LLVM
  const uint32_t v1 = _mm_extract_epi16(((__m128i)(v4)), ((int32_t)(2)));
  const uint32_t v2 = _mm_extract_epi16(((__m128i)(v4)), ((int32_t)(4)));
  const uint32_t v3 = _mm_extract_epi16(((__m128i)(v4)), ((int32_t)(6)));

  __m128 t0 = _mm_set_ss(float_exp_lookup[v0]);
  __m128 t1 = _mm_set_ss(float_exp_lookup[v1]);
  const __m128 t2 = _mm_set_ss(float_exp_lookup[v2]);
  const __m128 t3 = _mm_set_ss(float_exp_lookup[v3]);
  t1 = _mm_movelh_ps(t1, t3);
  t1 = _mm_castsi128_ps(_mm_slli_epi64(_mm_castps_si128(t1), 32));
  t0 = _mm_movelh_ps(t0, t2);
  t0 = _mm_or_ps(t0, t1);
  t0 = _mm_or_ps(t0, _mm_castsi128_ps(u4));
  tt = _mm_mul_ps(tt, t0);
  return tt;
}

// For equivalent "normal" C/C++ code, see the non-__LP64__ versions of these
// functions.
static inline void logistic_sse(uint32_t nn, float* vect) {
  const __m128 zero = _mm_setzero_ps();
  const __m128 one = _mm_set1_ps(1.0);
  for (uint32_t uii = 0; uii < nn; uii += 4) {
    __m128 aa = _mm_load_ps(&(vect[uii]));
    aa = _mm_sub_ps(zero, aa);
    aa = fmath_exp_ps(aa);
    aa = _mm_add_ps(aa, one);
    aa = _mm_div_ps(one, aa);
    _mm_store_ps(&(vect[uii]), aa);
  }
}

static inline void compute_v_and_p_minus_y(const float* yy, uint32_t nn, float* pp, float* vv) {
  const __m128 one = _mm_set1_ps(1.0);
  for (uint32_t uii = 0; uii < nn; uii += 4) {
    __m128 ptmp = _mm_load_ps(&(pp[uii]));
    __m128 one_minus_ptmp = _mm_sub_ps(one, ptmp);
    _mm_store_ps(&(vv[uii]), _mm_mul_ps(ptmp, one_minus_ptmp));
    __m128 ytmp = _mm_load_ps(&(yy[uii]));
    _mm_store_ps(&(pp[uii]), _mm_sub_ps(ptmp, ytmp));
  }
}

static inline void compute_v(const float* pp, uint32_t nn, float* vv) {
  const __m128 one = _mm_set1_ps(1.0);
  for (uint32_t uii = 0; uii < nn; uii += 4) {
    __m128 ptmp = _mm_load_ps(&(pp[uii]));
    __m128 one_minus_ptmp = _mm_sub_ps(one, ptmp);
    _mm_store_ps(&(vv[uii]), _mm_mul_ps(ptmp, one_minus_ptmp));
  }
}

static inline void mult_tmatrix_nxd_vect_d(const float* tm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
  // tm is row-major, cols are packed to 16-byte alignment
  // "col_cta4" = col_ct, aligned to multiple of 4.  Since 16-byte blocks
  // contain 4 floats each, this is the actual length (in floats) of each tm
  // row.  (Yes, I need to standardize a zillion other variable names of this
  // sort...)
  __m128 w1;
  __m128 w2;
  __m128 w3;
  const uintptr_t col_cta4 = round_up_pow2(col_ct, 4);
  uint32_t row_idx = 0;
  if (row_ct < 4) {
    memset(dest, 0, col_ct * sizeof(float));
  } else {
    w1 = _mm_load1_ps(vect);
    w2 = _mm_load1_ps(&(vect[1]));
    w3 = _mm_load1_ps(&(vect[2]));
    __m128 w4 = _mm_load1_ps(&(vect[3]));
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += 4) {
      __m128 r1 = _mm_load_ps(&(tm[col_idx]));
      __m128 r2 = _mm_load_ps(&(tm[col_idx + col_cta4]));
      __m128 r3 = _mm_load_ps(&(tm[col_idx + 2 * col_cta4]));
      __m128 r4 = _mm_load_ps(&(tm[col_idx + 3 * col_cta4]));
      r1 = _mm_mul_ps(r1, w1);
      r2 = _mm_mul_ps(r2, w2);
      r3 = _mm_mul_ps(r3, w3);
      r4 = _mm_mul_ps(r4, w4);
      r1 = _mm_add_ps(r1, r2);
      r3 = _mm_add_ps(r3, r4);
      r1 = _mm_add_ps(r1, r3);
      _mm_store_ps(&(dest[col_idx]), r1);
    }
    const uint32_t row_ctm3 = row_ct - 3;
    for (row_idx = 4; row_idx < row_ctm3; row_idx += 4) {
      w1 = _mm_load1_ps(&(vect[row_idx]));
      w2 = _mm_load1_ps(&(vect[row_idx + 1]));
      w3 = _mm_load1_ps(&(vect[row_idx + 2]));
      w4 = _mm_load1_ps(&(vect[row_idx + 3]));
      for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += 4) {
        __m128 r1 = _mm_load_ps(&(tm[col_idx + row_idx * col_cta4]));
        __m128 r2 = _mm_load_ps(&(tm[col_idx + (row_idx + 1) * col_cta4]));
        __m128 r3 = _mm_load_ps(&(tm[col_idx + (row_idx + 2) * col_cta4]));
        __m128 r4 = _mm_load_ps(&(tm[col_idx + (row_idx + 3) * col_cta4]));
        r1 = _mm_mul_ps(r1, w1);
        r2 = _mm_mul_ps(r2, w2);
        r3 = _mm_mul_ps(r3, w3);
        r4 = _mm_mul_ps(r4, w4);
        r1 = _mm_add_ps(r1, r2);
        r3 = _mm_add_ps(r3, r4);
        r1 = _mm_add_ps(r1, r3);
	r1 = _mm_add_ps(r1, _mm_load_ps(&(dest[col_idx])));
	_mm_store_ps(&(dest[col_idx]), r1);
      }
    }
  }
  switch(row_ct % 4) {
  case 3:
    w1 = _mm_load1_ps(&(vect[row_idx]));
    w2 = _mm_load1_ps(&(vect[row_idx + 1]));
    w3 = _mm_load1_ps(&(vect[row_idx + 2]));
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += 4) {
      __m128 r1 = _mm_load_ps(&(tm[col_idx + row_idx * col_cta4]));
      __m128 r2 = _mm_load_ps(&(tm[col_idx + (row_idx + 1) * col_cta4]));
      __m128 r3 = _mm_load_ps(&(tm[col_idx + (row_idx + 2) * col_cta4]));
      r1 = _mm_mul_ps(r1, w1);
      r2 = _mm_mul_ps(r2, w2);
      r3 = _mm_mul_ps(r3, w3);
      r1 = _mm_add_ps(r1, r2);
      r3 = _mm_add_ps(r3, _mm_load_ps(&(dest[col_idx])));
      r1 = _mm_add_ps(r1, r3);
      _mm_store_ps(&(dest[col_idx]), r1);
    }
    break;
  case 2:
    w1 = _mm_load1_ps(&(vect[row_idx]));
    w2 = _mm_load1_ps(&(vect[row_idx + 1]));
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += 4) {
      __m128 r1 = _mm_load_ps(&(tm[col_idx + row_idx * col_cta4]));
      __m128 r2 = _mm_load_ps(&(tm[col_idx + (row_idx + 1) * col_cta4]));
      r1 = _mm_mul_ps(r1, w1);
      r2 = _mm_mul_ps(r2, w2);
      r1 = _mm_add_ps(r1, r2);
      r1 = _mm_add_ps(r1, _mm_load_ps(&(dest[col_idx])));
      _mm_store_ps(&(dest[col_idx]), r1);
    }
    break;
  case 1:
    w1 = _mm_load1_ps(&(vect[row_idx]));
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += 4) {
      __m128 r1 = _mm_load_ps(&(tm[col_idx + row_idx * col_cta4]));
      r1 = _mm_mul_ps(r1, w1);
      r1 = _mm_add_ps(r1, _mm_load_ps(&(dest[col_idx])));
      _mm_store_ps(&(dest[col_idx]), r1);
    }
  }
}

// This code was hand-optimized by others for 16-byte float vectors.  Exempt it
// from the rest of the codebase's attempt at vector-size-agnosticism for now.
typedef union {
  __m128 vf;
  float f4[4];
} __old_univecf_t;

static inline void mult_matrix_dxn_vect_n(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uintptr_t col_cta4 = round_up_pow2(col_ct, 4);
  uint32_t row_idx = 0;
  __m128 s1;
  __m128 s2;
  __m128 s3;
  __old_univecf_t uvec;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    for (; row_idx < row_ctm3; row_idx += 4) {
      s1 = _mm_setzero_ps();
      s2 = _mm_setzero_ps();
      s3 = _mm_setzero_ps();
      __m128 s4 = _mm_setzero_ps();
      for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += 4) {
	const float* mm_ptr = &(mm[row_idx * col_cta4 + col_idx]);
        const __m128 vv = _mm_load_ps(&(vect[col_idx]));
        __m128 a1 = _mm_load_ps(mm_ptr);
        __m128 a2 = _mm_load_ps(&(mm_ptr[col_cta4]));
        __m128 a3 = _mm_load_ps(&(mm_ptr[2 * col_cta4]));
        __m128 a4 = _mm_load_ps(&(mm_ptr[3 * col_cta4]));
	// want to switch this to fused multiply-add...
        a1 = _mm_mul_ps(a1, vv);
        a2 = _mm_mul_ps(a2, vv);
        a3 = _mm_mul_ps(a3, vv);
        a4 = _mm_mul_ps(a4, vv);
        s1 = _mm_add_ps(s1, a1);
        s2 = _mm_add_ps(s2, a2);
        s3 = _mm_add_ps(s3, a3);
        s4 = _mm_add_ps(s4, a4);
      }
      // refrain from using SSE3 _mm_hadd_ps() for now
      uvec.vf = s1;
      *dest++ = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
      uvec.vf = s2;
      *dest++ = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
      uvec.vf = s3;
      *dest++ = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
      uvec.vf = s4;
      *dest++ = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
    }
  }
  s1 = _mm_setzero_ps();
  s2 = _mm_setzero_ps();
  s3 = _mm_setzero_ps();
  switch (row_ct % 4) {
  case 3:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += 4) {
      const float* mm_ptr = &(mm[row_idx * col_cta4 + col_idx]);
      const __m128 vv = _mm_load_ps(&(vect[col_idx]));
      __m128 a1 = _mm_load_ps(mm_ptr);
      __m128 a2 = _mm_load_ps(&(mm_ptr[col_cta4]));
      __m128 a3 = _mm_load_ps(&(mm_ptr[2 * col_cta4]));
      a1 = _mm_mul_ps(a1, vv);
      a2 = _mm_mul_ps(a2, vv);
      a3 = _mm_mul_ps(a3, vv);
      s1 = _mm_add_ps(s1, a1);
      s2 = _mm_add_ps(s2, a2);
      s3 = _mm_add_ps(s3, a3);
    }
    uvec.vf = s1;
    *dest++ = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
    uvec.vf = s2;
    *dest++ = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
    uvec.vf = s3;
    *dest = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
    break;
  case 2:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += 4) {
      const float* mm_ptr = &(mm[row_idx * col_cta4 + col_idx]);
      const __m128 vv = _mm_load_ps(&(vect[col_idx]));
      __m128 a1 = _mm_load_ps(mm_ptr);
      __m128 a2 = _mm_load_ps(&(mm_ptr[col_cta4]));
      a1 = _mm_mul_ps(a1, vv);
      a2 = _mm_mul_ps(a2, vv);
      s1 = _mm_add_ps(s1, a1);
      s2 = _mm_add_ps(s2, a2);
    }
    uvec.vf = s1;
    *dest++ = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
    uvec.vf = s2;
    *dest = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
    break;
  case 1:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += 4) {
      const __m128 vv = _mm_load_ps(&(vect[col_idx]));
      __m128 a1 = _mm_load_ps(&(mm[row_idx * col_cta4 + col_idx]));
      a1 = _mm_mul_ps(a1, vv);
      s1 = _mm_add_ps(s1, a1);
    }
    uvec.vf = s1;
    *dest = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
    break;
  }
}

static inline float triple_product(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  __m128 sum = _mm_setzero_ps();
  for (uint32_t uii = 0; uii < nn; uii += 4) {
    const __m128 aa = _mm_load_ps(&(v1[uii]));
    const __m128 bb = _mm_load_ps(&(v2[uii]));
    const __m128 cc = _mm_load_ps(&(v3[uii]));
    sum = _mm_add_ps(sum, _mm_mul_ps(_mm_mul_ps(aa, bb), cc));
  }
  __old_univecf_t uvec;
  uvec.vf = sum;
  return uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
}

static inline void compute_two_diag_triple_product(const float* aa, const float* bb, const float* vv, uint32_t nn, float* raa_ptr, float* rab_ptr, float* rbb_ptr) {
  __m128 saa = _mm_setzero_ps();
  __m128 sab = _mm_setzero_ps();
  __m128 sbb = _mm_setzero_ps();
  for (uint32_t uii = 0; uii < nn; uii += 4) {
    const __m128 vtmp = _mm_load_ps(&(vv[uii]));
    const __m128 atmp = _mm_load_ps(&(aa[uii]));
    const __m128 btmp = _mm_load_ps(&(bb[uii]));
    const __m128 av = _mm_mul_ps(atmp, vtmp);
    const __m128 bv = _mm_mul_ps(btmp, vtmp);
    saa = _mm_add_ps(saa, _mm_mul_ps(atmp, av));
    sab = _mm_add_ps(sab, _mm_mul_ps(atmp, bv));
    sbb = _mm_add_ps(sbb, _mm_mul_ps(btmp, bv));
  }
  __old_univecf_t uvec;
  uvec.vf = saa;
  *raa_ptr = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
  uvec.vf = sab;
  *rab_ptr = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
  uvec.vf = sbb;
  *rbb_ptr = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
}

static inline void compute_three_triple_product(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* r1_ptr, float* r2_ptr, float* r3_ptr) {
  __m128 s1 = _mm_setzero_ps();
  __m128 s2 = _mm_setzero_ps();
  __m128 s3 = _mm_setzero_ps();
  for (uint32_t uii = 0; uii < nn; uii += 4) {
    const __m128 a1tmp = _mm_load_ps(&(a1[uii]));
    const __m128 a2tmp = _mm_load_ps(&(a2[uii]));
    const __m128 a3tmp = _mm_load_ps(&(a3[uii]));
    const __m128 vtmp = _mm_load_ps(&(vv[uii]));
    __m128 btmp = _mm_load_ps(&(bb[uii]));
    btmp = _mm_mul_ps(btmp, vtmp);
    s1 = _mm_add_ps(s1, _mm_mul_ps(a1tmp, btmp));
    s2 = _mm_add_ps(s2, _mm_mul_ps(a2tmp, btmp));
    s3 = _mm_add_ps(s3, _mm_mul_ps(a3tmp, btmp));
  }
  __old_univecf_t uvec;
  uvec.vf = s1;
  *r1_ptr = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
  uvec.vf = s2;
  *r2_ptr = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
  uvec.vf = s3;
  *r3_ptr = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
}

static inline void compute_two_plus_one_triple_product(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* r1_ptr, float* r2_ptr, float* r3_ptr) {
  __m128 s1 = _mm_setzero_ps();
  __m128 s2 = _mm_setzero_ps();
  __m128 s3 = _mm_setzero_ps();
  for (uint32_t uii = 0; uii < nn; uii += 4) {
    const __m128 a1tmp = _mm_load_ps(&(a1[uii]));
    const __m128 a2tmp = _mm_load_ps(&(a2[uii]));
    const __m128 btmp = _mm_load_ps(&(bb[uii]));
    const __m128 vtmp = _mm_load_ps(&(vv[uii]));
    const __m128 bv = _mm_mul_ps(btmp, vtmp);
    s1 = _mm_add_ps(s1, _mm_mul_ps(btmp, bv));
    s2 = _mm_add_ps(s2, _mm_mul_ps(a1tmp, bv));
    s3 = _mm_add_ps(s3, _mm_mul_ps(a2tmp, bv));
  }
  __old_univecf_t uvec;
  uvec.vf = s1;
  *r1_ptr = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
  uvec.vf = s2;
  *r2_ptr = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
  uvec.vf = s3;
  *r3_ptr = uvec.f4[0] + uvec.f4[1] + uvec.f4[2] + uvec.f4[3];
}
#else // no __LP64__ (and hence, unsafe to assume presence of SSE2)
static inline void logistic_sse(uint32_t nn, float* vect) {
  for (uint32_t uii = 0; uii < nn; ++uii) {
    vect[uii] = 1.0 / (1 + exp(-vect[uii]));
  }
}

static inline void compute_v_and_p_minus_y(const float* yy, uint32_t nn, float* pp, float* vv) {
  for (uint32_t uii = 0; uii < nn; ++uii) {
    vv[uii] = pp[uii] * (1.0 - pp[uii]);
    pp[uii] -= yy[uii];
  }
}

static inline void compute_v(const float* pp, uint32_t nn, float* vv) {
  for (uint32_t uii = 0; uii < nn; ++uii) {
    vv[uii] = pp[uii] * (1.0 - pp[uii]);
  }
}

static inline void mult_tmatrix_nxd_vect_d(const float* tm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uintptr_t col_cta4 = round_up_pow2(col_ct, 4);
  fill_float_zero(col_ct, dest);
  for (uint32_t row_idx = 0; row_idx < row_ct; ++row_idx) {
    const float vect_val = vect[row_idx];
    const float* tm_ptr = &(tm[row_idx * col_cta4]);
    for (uint32_t col_idx = 0; col_idx < col_ct; ++col_idx) {
      dest[col_idx] += (*tm_ptr++) * vect_val;
    }
  }
}

static inline void mult_matrix_dxn_vect_n(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uintptr_t col_cta4 = round_up_pow2(col_ct, 4);
  for (uint32_t row_idx = 0; row_idx < row_ct; ++row_idx) {
    float fxx = 0.0;
    const float* vect_ptr = vect;
    const float* mm_ptr = &(mm[row_idx * col_cta4]);
    for (uint32_t col_idx = 0; col_idx < col_ct; ++col_idx) {
      fxx += (*mm_ptr++) * (*vect_ptr++);
    }
    *dest++ = fxx;
  }
}

static inline float triple_product(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  float fxx = 0.0;
  for (uint32_t uii = 0; uii < nn; ++uii) {
    fxx += (*v1++) * (*v2++) * (*v3++);
  }
  return fxx;
}

static inline void compute_two_diag_triple_product(const float* aa, const float* bb, const float* vv, uint32_t nn, float* raa_ptr, float* rab_ptr, float* rbb_ptr) {
  float raa = 0.0;
  float rab = 0.0;
  float rbb = 0.0;
  for (uint32_t uii = 0; uii < nn; ++uii) {
    const float fxx = (*aa++);
    const float fyy = (*bb++);
    float fzz = (*vv++);
    raa += fxx * fxx * fzz;
    fzz *= fyy;
    rab += fxx * fzz;
    rbb += fyy * fzz;
  }
  *raa_ptr = raa;
  *rab_ptr = rab;
  *rbb_ptr = rbb;
}

static inline void compute_three_triple_product(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* r1_ptr, float* r2_ptr, float* r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii < nn; ++uii) {
    const float fxx = (*bb++) * (*vv++);
    r1 += (*a1++) * fxx;
    r2 += (*a2++) * fxx;
    r3 += (*a3++) * fxx;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}

static inline void compute_two_plus_one_triple_product(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* r1_ptr, float* r2_ptr, float* r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii < nn; ++uii) {
    const float fxx = (*bb++);
    const float fyy = fxx * (*vv++);
    r1 += fxx * fyy;
    r2 += (*a1++) * fyy;
    r3 += (*a2++) * fyy;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}
#endif
double compute_loglik(const float* yy, const float* pp, uint32_t sample_ct) {
  double loglik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
    const float new_pi = pp[sample_idx];
    loglik += (yy[sample_idx])? log(new_pi) : log(1.0 - new_pi);
  }
  return loglik;
}

static inline void compute_hessian(const float* mm, const float* vv, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uintptr_t col_cta4 = round_up_pow2(col_ct, 4);
  const uintptr_t row_cta4 = round_up_pow2(row_ct, 4);
  const uintptr_t row_cta4p1 = row_cta4 + 1;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    for (uint32_t row_idx = 0; row_idx < row_ctm3; row_idx += 3) {
      const float* mm_cur = &(mm[row_idx * col_cta4]);
      compute_two_diag_triple_product(mm_cur, &(mm_cur[col_cta4]), vv, col_ct, &(dest[row_idx * row_cta4p1]), &(dest[(row_idx + 1) * row_cta4p1 - 1]), &(dest[(row_idx + 1) * row_cta4p1]));
      compute_two_plus_one_triple_product(&(mm_cur[2 * col_cta4]), &(mm_cur[col_cta4]), mm_cur, vv, col_ct, &(dest[(row_idx + 2) * row_cta4p1]), &(dest[(row_idx + 2) * row_cta4p1 - 1]), &(dest[(row_idx + 2) * row_cta4p1 - 2]));
      for (uint32_t row_idx2 = row_idx + 3; row_idx2 < row_ct; row_idx2++) {
        compute_three_triple_product(&(mm[row_idx2 * col_cta4]), mm_cur, &(mm_cur[col_cta4]), &(mm_cur[2 * col_cta4]), vv, col_ct, &(dest[row_idx2 * row_cta4 + row_idx]), &(dest[row_idx2 * row_cta4 + row_idx + 1]), &(dest[row_idx2 * row_cta4 + row_idx + 2]));
      }
    }
  }
  switch (row_ct % 3) {
  case 0:
    compute_two_plus_one_triple_product(&(mm[(row_ct - 3) * col_cta4]), &(mm[(row_ct - 2) * col_cta4]), &(mm[(row_ct - 1) * col_cta4]), vv, col_ct, &(dest[(row_ct - 3) * row_cta4p1]), &(dest[(row_ct - 2) * row_cta4p1 - 1]), &(dest[(row_ct - 1) * row_cta4p1 - 2]));
    // fall through
  case 2:
    compute_two_diag_triple_product(&(mm[(row_ct - 2) * col_cta4]), &(mm[(row_ct - 1) * col_cta4]), vv, col_ct, &(dest[(row_ct - 2) * row_cta4p1]), &(dest[(row_ct - 1) * row_cta4p1 - 1]), &(dest[(row_ct - 1) * row_cta4p1]));
    break;
  case 1:
    dest[(row_ct - 1) * row_cta4p1] = triple_product(&(mm[(row_ct - 1) * col_cta4]), &(mm[(row_ct - 1) * col_cta4]), vv, col_ct);
  }
}

void cholesky_decomposition(const float* aa, uint32_t predictor_ct, float* ll) {
  const uintptr_t predictor_cta4 = round_up_pow2(predictor_ct, 4);
  const uintptr_t predictor_cta4p1 = predictor_cta4 + 1;
  for (uint32_t row_idx = 0; row_idx < predictor_ct; ++row_idx) {
    float fxx = aa[row_idx * predictor_cta4p1];
    float* ll_row_iter = &(ll[row_idx * predictor_cta4]);
    for (uint32_t col_idx = 0; col_idx < row_idx; ++col_idx) {
      const float fyy = (*ll_row_iter++);
      fxx -= fyy * fyy;
    }
    float fyy;
    if (fxx >= 0.0) {
      fyy = sqrtf(fxx);
    } else {
      fyy = 1e-6;
    }
    ll[row_idx * predictor_cta4p1] = fyy;
    fyy = 1.0 / fyy; // now 1.0 / L[j][j]
    for (uint32_t row_idx2 = row_idx + 1; row_idx2 < predictor_ct; ++row_idx2) {
      float fxx2 = aa[row_idx2 * predictor_cta4 + row_idx];
      float* ll_row_iter2 = &(ll[row_idx * predictor_cta4]);
      float* ll_row_iter3 = &(ll[row_idx2 * predictor_cta4]);
      for (uint32_t col_idx = 0; col_idx < row_idx; ++col_idx) {
        fxx2 -= (*ll_row_iter2++) * (*ll_row_iter3++);
      }
      ll[row_idx2 * predictor_cta4 + row_idx] = fxx2 * fyy;
    }
  }
}

void solve_linear_system(const float* ll, const float* yy, uint32_t predictor_ct, float* xx) {
  // Finds x such that y = L(L^T)x, via forward and backward substitution
  //
  // might want to use this in NOLAPACK case only, since we can now produce
  // 32-bit Linux builds with statically linked LAPACK
  const uintptr_t predictor_cta4 = round_up_pow2(predictor_ct, 4);
  for (uint32_t row_idx = 0; row_idx < predictor_ct; ++row_idx) {
    float fxx = yy[row_idx];
    const float* ll_row_iter = &(ll[row_idx * predictor_cta4]);
    float* xx_iter = xx;
    for (uint32_t col_idx = 0; col_idx < row_idx; ++col_idx) {
      fxx -= (*ll_row_iter++) * (*xx_iter++);
    }
    *xx_iter = fxx / (*ll_row_iter);
  }
  for (uint32_t col_idx = predictor_ct; col_idx;) {
    float fxx = xx[--col_idx];
    float* xx_iter = &(xx[predictor_ct - 1]);
    for (uint32_t row_idx = predictor_ct - 1; row_idx > col_idx; --row_idx) {
      fxx -= ll[row_idx * predictor_cta4 + col_idx] * (*xx_iter--);
    }
    *xx_iter = fxx / ll[col_idx * (predictor_cta4 + 1)];
  }
}

boolerr_t logistic_regression(const float* yy, const float* xx, uint32_t sample_ct, uint32_t predictor_ct, float* coef, float* ll, float* pp, float* vv, float* hh, float* grad, float* dcoef) {
  // Similar to first part of logistic.cpp fitLM(), but incorporates changes
  // from Pascal Pons et al.'s TopCoder code.
  //
  // Preallocated buffers (initial contents irrelevant):
  // vv    = sample variance buffer
  // hh    = hessian matrix buffer, predictor_ct^2, rows 16-byte aligned
  // grad  = gradient buffer Y[] (length predictor_ct)
  // dcoef = current coefficient change buffer (length predictor_ct)
  // 
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         16-byte aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype
  //
  // Input/output:
  // coef  = starting point, overwritten with logistic regression betas.  Must
  //         be 16-byte aligned.
  //
  // Outputs:
  // ll    = cholesky decomposition matrix, predictor_ct^2, rows 16-byte aligned
  // pp    = final likelihoods minus Y[] (not currently used by callers)
  //
  // Returns 0 on success, 1 on convergence failure.
  const uintptr_t predictor_cta4 = round_up_pow2(predictor_ct, 4);
  uint32_t iteration = 0;
  float min_delta_coef = 1e9;

  fill_float_zero(predictor_ct * predictor_cta4, ll);
  while (1) {
    ++iteration;

    // P[i] = \sum_j coef[j] * X[i][j];
    mult_tmatrix_nxd_vect_d(xx, coef, sample_ct, predictor_ct, pp);

    // P[i] = 1 / (1 + exp(-P[i]));
    logistic_sse(sample_ct, pp);

    // V[i] = P[i] * (1 - P[i]);
    // P[i] -= Y[i];
    compute_v_and_p_minus_y(yy, sample_ct, pp, vv);

    compute_hessian(xx, vv, sample_ct, predictor_ct, hh);

    mult_matrix_dxn_vect_n(xx, pp, sample_ct, predictor_ct, grad);

    cholesky_decomposition(hh, predictor_ct, ll);

    // fill_float_zero(predictor_ct, dcoef);
    solve_linear_system(ll, grad, predictor_ct, dcoef);

    float delta_coef = 0.0;
    for (uint32_t pred_idx = 0; pred_idx < predictor_ct; pred_idx++) {
      const float cur_dcoef = dcoef[pred_idx];
      delta_coef += fabsf(cur_dcoef);
      coef[pred_idx] -= cur_dcoef;
    }
    if (delta_coef < min_delta_coef) {
      min_delta_coef = delta_coef;
    }
    if (delta_coef != delta_coef) {
      return 1;
    }
    if (iteration > 4) {
      if (((delta_coef > 20.0) && (delta_coef > 2 * min_delta_coef)) || ((iteration >= 8) && fabsf(1.0f - delta_coef) < 1e-3)) {
	return 1;
      }
      if (iteration >= 15) {
	return 0;
      }
    }
    // Pons reported that 1.1e-3 was dangerous, so I agree with the decision to
    // tighten this threshold from 1e-3 to 1e-4.
    if (delta_coef < 1e-4) {
      return 0;
    }
  }
}

boolerr_t firth_regression(const float* yy, const float* xx, uint32_t sample_ct, uint32_t predictor_ct, float* coef, float* hh, matrix_finvert_buf1_t* inv_1d_buf, float* flt_2d_buf, float* pp, float* vv, float* grad, float* dcoef, float* ww, float* tmpnxk_buf) {
  // This is a port of Georg Heinze's logistf R function, adapted to use many
  // of plink 1.9's optimizations; see
  //   http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/fllogistf/
  //
  // Preallocated buffers (initial contents irrelevant):
  // inv_1d_buf, flt_2d_buf = for float32 matrix inversion
  // pp    = likelihoods minus Y[] (not currently used by callers)
  // vv    = sample variance buffer
  // grad  = gradient buffer (length predictor_ct)
  // dcoef = current coefficient change buffer (length predictor_ct)
  // ww    = Firth-adjusted scores, sample_ct
  // 
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         16-byte aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype
  //
  // Input/output:
  // coef  = starting point, overwritten with logistic regression betas.  Must
  //         be 16-byte aligned.
  //
  // Outputs:
  // hh    = variance-covariance matrix buffer, predictor_ct^2, rows 16-byte
  //         aligned.  (spends some time as pre-inversion Hessian matrix too)
  //
  // Returns 0 on success, 1 on convergence failure.
  const uintptr_t predictor_cta4 = round_up_pow2(predictor_ct, 4);
  const uintptr_t sample_cta4 = round_up_pow2(sample_ct, 4);
  uint32_t is_last_iter = 0;
  
  // pull these out of the start of the loop, since they happen again in the
  // log-likelihood update
  // P[i] = \sum_j coef[j] * X[i][j];
  mult_tmatrix_nxd_vect_d(xx, coef, sample_ct, predictor_ct, pp);
  // P[i] = 1 / (1 + exp(-P[i]));
  logistic_sse(sample_ct, pp);
  // V[i] = P[i] * (1 - P[i]);
  compute_v(pp, sample_ct, vv);
  // P[i] -= Y[i] NOT done here

  // hessian = X diag(V) X'
  // note that only lower triangle is filled here
  compute_hessian(xx, vv, sample_ct, predictor_ct, hh);

  for (uint32_t uii = 0; uii < predictor_ct; ++uii) {
    for (uint32_t ujj = uii + 1; ujj < predictor_ct; ++ujj) {
      hh[uii * predictor_cta4 + ujj] = hh[ujj * predictor_cta4 + uii];
    }
  }
  // we shouldn't need to compute the log directly, since underflow <->
  // regression failure, right?  check this.
  float dethh;
  if (invert_fmatrix_first_half(predictor_ct, predictor_cta4, hh, &dethh, inv_1d_buf, flt_2d_buf)) {
    return 1;
  }
  /*
  if (sample_ct < sample_cta4) {
    // trailing Y[] values must be zero
    fill_float_zero(sample_cta4 - sample_ct, &(pp[sample_ct]));
  }
  */
  double loglik = compute_loglik(yy, pp, sample_ct);
  // printf("loglik: %g\n", loglik);
  loglik += 0.5 * log(dethh);

  uint32_t iter_idx = 0;
  // start with 80% of logistf convergence defaults (some reduction is
  // appropriate to be consistent with single-precision arithmetic); may tune
  // later.
  // see also the hs_bail condition: if we ever try all five halfsteps, when
  // dcoef_max and grad_max aren't that far from the normal convergence
  // conditions, it's probably pointless to continue with single-precision
  // arithmetic.  (possible todo: use a fully-double-precision routine to
  // finish the job when that happens.)
  const uint32_t max_iter = 20;
  const float gconv = 0.0001;
  const float xconv = 0.0001;
  const double lconv = 0.0001;
  uint32_t hs_bail = 0;
  while (1) {
    invert_fmatrix_second_half(predictor_ct, predictor_cta4, hh, inv_1d_buf, flt_2d_buf);
    if (is_last_iter) {
      return 0;
    }
    col_major_fmatrix_multiply_strided(xx, hh, sample_ct, sample_cta4, predictor_ct, predictor_cta4, predictor_ct, sample_cta4, tmpnxk_buf);
    // tmpNxK, interpreted as column-major, is sample_ct x predictor_ct
    // X, interpreted as column-major, is also sample_ct x predictor_ct
    // Hdiag[i] = V[i] (\sum_j tmpNxK[i][j] X[i][j])
    // (todo: vectorize this)
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      float dotprod = 0.0;
      const float* xx_row = &(xx[sample_idx]);
      const float* tmpnxk_row = &(tmpnxk_buf[sample_idx]);
      for (uint32_t pred_uidx = 0; pred_uidx < predictor_ct; ++pred_uidx) {
	dotprod += xx_row[pred_uidx * sample_cta4] * tmpnxk_row[pred_uidx * sample_cta4];
      }
      const float cur_weight = vv[sample_idx];
      const float cur_pi = pp[sample_idx];
      ww[sample_idx] = (yy[sample_idx] - cur_pi) + (0.5 - cur_pi) * cur_weight * dotprod;
    }
    
    // gradient (Ustar in logistf) = X' W
    mult_matrix_dxn_vect_n(xx, ww, sample_ct, predictor_ct, grad);
    float grad_max = 0.0;
    for (uint32_t pred_uidx = 0; pred_uidx < predictor_ct; ++pred_uidx) {
      const float abs_grad_cur = fabsf(grad[pred_uidx]);
      if (abs_grad_cur > grad_max) {
	grad_max = abs_grad_cur;
      }
    }

    // dcoef := hh * grad (note that hh is inverted already)
    mult_matrix_dxn_vect_n(hh, grad, predictor_ct, predictor_ct, dcoef);

    float dcoef_max = 0.0;
    for (uint32_t pred_uidx = 0; pred_uidx < predictor_ct; ++pred_uidx) {
      const float abs_dcoef_cur = fabsf(dcoef[pred_uidx]);
      if (abs_dcoef_cur > dcoef_max) {
	dcoef_max = abs_dcoef_cur;
      }
    }
    const float maxstep = 5.0;
    if (dcoef_max > maxstep) {
      const float scaling_factor = maxstep / dcoef_max;
      for (uint32_t pred_uidx = 0; pred_uidx < predictor_ct; ++pred_uidx) {
	dcoef[pred_uidx] *= scaling_factor;
      }
      dcoef_max = maxstep;
    }
    for (uint32_t pred_uidx = 0; pred_uidx < predictor_ct; ++pred_uidx) {
      coef[pred_uidx] += dcoef[pred_uidx];
    }
    const uint32_t delta_and_grad_converged = (dcoef_max <= xconv) && (grad_max < gconv);
    const double loglik_old = loglik;
    double loglik_thresh = loglik_old;
    if (delta_and_grad_converged) {
      // on the last iteration, we would frequently try all 5 halfsteps when
      // the log-likelihood change was effectively random due to floating point
      // error.  detect this and exit the loop earlier.
      loglik_thresh -= 0.999999 * lconv;
    }

    uint32_t maxhs = 5;
    uint32_t halfstep_idx = 1;
    while (1) {
      mult_tmatrix_nxd_vect_d(xx, coef, sample_ct, predictor_ct, pp);
      logistic_sse(sample_ct, pp);
      loglik = compute_loglik(yy, pp, sample_ct);
      compute_v(pp, sample_ct, vv);
      compute_hessian(xx, vv, sample_ct, predictor_ct, hh);
      for (uint32_t uii = 0; uii < predictor_ct; ++uii) {
	for (uint32_t ujj = uii + 1; ujj < predictor_ct; ++ujj) {
	  hh[uii * predictor_cta4 + ujj] = hh[ujj * predictor_cta4 + uii];
	}
      }
      if (invert_fmatrix_first_half(predictor_ct, predictor_cta4, hh, &dethh, inv_1d_buf, flt_2d_buf)) {
	return 1;
      }
      loglik += 0.5 * log(dethh);
      if (halfstep_idx > maxhs) {
	break;
      }
      if (loglik >= loglik_thresh) {
	if (loglik >= loglik_old) {
	  break;
	}
	maxhs = halfstep_idx;
      } else if (halfstep_idx == maxhs) {
	if ((dcoef_max < 0.001) && (grad_max < 0.05) && (loglik >= loglik_old - lconv)) {
	  // we've converged as much as we can with single-precision
	  // arithmetic, and now we're flailing around.  don't even take the
	  // 2^{-maxhs} step, undo it all and bail.
	  // (0.001 and 0.05 constants can obviously be tuned; they were chosen
	  // based on a test 500k sample/5 covariate regression.)
	  --halfstep_idx;
	  --maxhs;
	  hs_bail = 1;
	}
      }
      const float multiplier = exp2f(-((int32_t)halfstep_idx));
      for (uint32_t pred_uidx = 0; pred_uidx < predictor_ct; ++pred_uidx) {
	coef[pred_uidx] -= dcoef[pred_uidx] * multiplier;
      }
      ++halfstep_idx;
    }
    // printf("%.9g %.9g %g %g\n", loglik, loglik_old, dcoef_max, grad_max);
    const double loglik_change = loglik - loglik_old;
    ++iter_idx;
    is_last_iter = (iter_idx == max_iter) || ((fabs(loglik_change) <= lconv) && (delta_and_grad_converged || hs_bail));
  }
}

uintptr_t get_logistic_workspace_size(uint32_t sample_ct, uint32_t predictor_ct, uint32_t constraint_ct, uint32_t genof_buffer_needed, uint32_t is_sometimes_firth) {
  // sample_cta4 * predictor_ct < 2^31, and sample_ct >= predictor_ct, so no
  // overflows
  // could round everything up to multiples of 16 instead of 64
  const uint32_t sample_cta4 = round_up_pow2(sample_ct, 4);
  const uint32_t predictor_cta4 = round_up_pow2(predictor_ct, 4);
  // sample_nm, pheno_cc_nm, male_nm = sample_ctl words
  uintptr_t workspace_size = 3 * round_up_pow2(BITCT_TO_WORDCT(sample_ct) * sizeof(intptr_t), kCacheline);
  
  // yy = sample_cta4 floats
  workspace_size += round_up_pow2(sample_cta4 * sizeof(float), kCacheline);
  
  // xx = (predictor_ct + genof_buffer_needed) * sample_cta4 floats
  workspace_size += round_up_pow2((predictor_ct + genof_buffer_needed) * sample_cta4 * sizeof(float), kCacheline);
    
  // hh = predictor_ct * predictor_cta4 floats
  workspace_size += round_up_pow2(predictor_ct * predictor_cta4 * sizeof(float), kCacheline);

  // pp, vv = sample_cta4 floats
  workspace_size += 2 * round_up_pow2(sample_cta4 * sizeof(float), kCacheline);

  // coef, grad, dcoef = predictor_cta4 floats
  workspace_size += 3 * round_up_pow2(predictor_cta4 * sizeof(float), kCacheline);

  // ll = predictor_ct * predictor_cta4 floats
  // (technically not needed in pure-Firth case)
  workspace_size += round_up_pow2(predictor_ct * predictor_cta4 * sizeof(float), kCacheline);

  if (is_sometimes_firth || constraint_ct) {
    // inv_1d_buf = predictor_ct * kMatrixFinvertBuf1CheckedAlloc bytes
    workspace_size += round_up_pow2(predictor_ct * kMatrixFinvertBuf1CheckedAlloc, kCacheline);

    // flt_2d_buf = predictor_ct * predictor_cta4 floats
    workspace_size += round_up_pow2(predictor_ct * predictor_cta4 * sizeof(float), kCacheline);

    if (is_sometimes_firth) {
      // ww = sample_cta4 floats
      workspace_size += round_up_pow2(sample_cta4 * sizeof(float), kCacheline);

      // tmpnxk_buf = predictor_ct * sample_cta4 floats
      workspace_size += round_up_pow2(predictor_ct * sample_cta4 * sizeof(float), kCacheline);
    }
    if (constraint_ct) {
      // tmphxs_buf, h_transpose_buf = constraint_ct * predictor_cta4 floats
      workspace_size += 2 * round_up_pow2(constraint_ct * predictor_cta4 * sizeof(float), kCacheline);
      
      // inner_buf = constraint_ct * constraint_ct
      workspace_size += round_up_pow2(constraint_ct * constraint_ct * sizeof(float), kCacheline);
    }
  }
  return workspace_size;  
}


typedef struct {
  // double beta;
  //   odds ratio = exp(beta)
  // double se;
  //   zval = beta / se
  //   width of asymptotic CI (beta units) = ci_zt * se
  //   T-statistic = zval
  //   pval = chiprob_p(zval * zval, 1);
  
  uint32_t sample_obs_ct;
  
  uint32_t allele_obs_ct;
  double alt_dosage;

  uint32_t firth_fallback;
  uint32_t case_allele_obs_ct;
  double alt_case_dosage;

  double mach_r2;
} logistic_aux_result_t;

typedef struct {
  // double beta;
  // double se;
  //   zval = beta / se
  //   width of asymptotic CI = ci_zt * se
  //   T-statistic = zval
  //   pval = calc_tprob(zval, sample_obs_ct - predictor_ct)

  uint32_t sample_obs_ct;

  uint32_t allele_obs_ct;
  double alt_dosage;

  double mach_r2;
} linear_aux_result_t;

// multithread globals
static pgen_reader_t** g_pgr_ptrs = nullptr;
static uintptr_t** g_genovecs = nullptr;
static uintptr_t** g_dosage_presents = nullptr;
static dosage_t** g_dosage_val_bufs = nullptr;
static unsigned char** g_workspace_bufs = nullptr;
static uint32_t* g_read_variant_uidx_starts = nullptr;

static uintptr_t* g_sample_include = nullptr;
static const uintptr_t* g_sample_include_x = nullptr;
static const uintptr_t* g_sample_include_y = nullptr;
static uint32_t* g_sample_include_cumulative_popcounts = nullptr;
static uint32_t* g_sample_include_x_cumulative_popcounts = nullptr;
static uint32_t* g_sample_include_y_cumulative_popcounts = nullptr;
static const uintptr_t* g_sex_male_collapsed = nullptr;
static const uintptr_t* g_pheno_cc = nullptr;
static uintptr_t* g_pheno_x_cc = nullptr;
static uintptr_t* g_pheno_y_cc = nullptr;
static const float* g_pheno_f = nullptr;
static float* g_pheno_x_f = nullptr;
static float* g_pheno_y_f = nullptr;
static uintptr_t* g_parameter_subset = nullptr;
static uintptr_t* g_parameter_subset_x = nullptr;
static uintptr_t* g_parameter_subset_y = nullptr;
static const float* g_covars_cmaj_f = nullptr;
static float* g_covars_cmaj_x_f = nullptr;
static float* g_covars_cmaj_y_f = nullptr;
static float* g_local_covars_vcmaj_f[2] = {nullptr, nullptr};
static const double* g_pheno_d = nullptr;
static double* g_pheno_x_d = nullptr;
static double* g_pheno_y_d = nullptr;
static const double* g_covars_cmaj_d = nullptr;
static double* g_covars_cmaj_x_d = nullptr;
static double* g_covars_cmaj_y_d = nullptr;
static double* g_local_covars_vcmaj_d[2] = {nullptr, nullptr};
// static const double* g_covar_dotprod_d = nullptr;
// static double* g_covar_dotprod_x_d = nullptr;
// static double* g_covar_dotprod_y_d = nullptr;
static const uintptr_t* g_variant_include = nullptr;
static const chr_info_t* g_cip = nullptr;
static const uintptr_t* g_variant_allele_idxs = nullptr;
static uint32_t* g_subset_chr_fo_vidx_start = nullptr;

// static uint32_t g_raw_sample_ct = 0;
static uint32_t g_sample_ct = 0;
static uint32_t g_sample_ct_x = 0;
static uint32_t g_sample_ct_y = 0;

// chrX value.  always equal to sample_ct_y on chrY, irrelevant elsewhere.
// (commented out since current algorithms always need nonmissing-male count)
// static uint32_t g_male_ct = 0;

static uint32_t g_covar_ct = 0;
static uint32_t g_local_covar_ct = 0;
static uint32_t g_covar_ct_x = 0;
static uint32_t g_covar_ct_y = 0;

static double* g_constraints_con_major = nullptr;
static double* g_constraints_con_major_x = nullptr;
static double* g_constraints_con_major_y = nullptr;
static float* g_constraints_con_major_f = nullptr;
static float* g_constraints_con_major_x_f = nullptr;
static float* g_constraints_con_major_y_f = nullptr;
static uint32_t g_constraint_ct = 0;
static uint32_t g_constraint_ct_x = 0;
static uint32_t g_constraint_ct_y = 0;

static uint32_t g_variant_ct = 0;
static uint32_t g_calc_thread_ct = 0;
static uint32_t g_cur_block_variant_ct = 0;
static glm_flags_t g_glm_flags = kfGlm0;
static uint32_t g_is_xchr_model_1 = 0;
static pglerr_t g_error_ret = kPglRetSuccess;

static logistic_aux_result_t* g_logistic_block_aux = nullptr;
static linear_aux_result_t* g_linear_block_aux = nullptr;

// separate from block_aux, since we need up to g_max_reported_test_ct pairs of
// values per variant
static double* g_block_beta_se = nullptr;

static uintptr_t g_max_reported_test_ct = 0;

THREAD_FUNC_DECL glm_logistic_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec = g_genovecs[tidx];
  uintptr_t* dosage_present = nullptr;
  dosage_t* dosage_vals = nullptr;
  if (g_dosage_presents) {
    dosage_present = g_dosage_presents[tidx];
    dosage_vals = g_dosage_val_bufs[tidx];
  }
  unsigned char* workspace_buf = g_workspace_bufs[tidx];
  const uintptr_t* variant_include = g_variant_include;
  const uintptr_t* sex_male_collapsed = g_sex_male_collapsed;
  const chr_info_t* cip = g_cip;
  const uint32_t* subset_chr_fo_vidx_start = g_subset_chr_fo_vidx_start;
  // const uint32_t raw_sample_ct = g_raw_sample_ct;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const glm_flags_t glm_flags = g_glm_flags;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t is_sometimes_firth = (glm_flags & (kfGlmFirthFallback | kfGlmFirth))? 1 : 0;
  const uint32_t is_always_firth = (glm_flags / kfGlmFirth) & 1;
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t joint_genotypic = (glm_flags / kfGlmGenotypic) & 1;
  const uint32_t joint_hethom = (glm_flags / kfGlmHethom) & 1;
  const uint32_t domdev_present = joint_genotypic || joint_hethom;
  const uint32_t domdev_present_p1 = domdev_present + 1;  
  const uint32_t reported_pred_uidx_start = 1 - include_intercept;
  const int32_t x_code = cip->xymt_codes[kChrOffsetX];
  const int32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = g_is_xchr_model_1;
  const uintptr_t max_reported_test_ct = g_max_reported_test_ct;
  const uintptr_t local_covar_ct = g_local_covar_ct;
  uintptr_t max_sample_ct = MAXV(g_sample_ct, g_sample_ct_x);
  if (max_sample_ct < g_sample_ct_y) {
    max_sample_ct = g_sample_ct_y;
  }
  uint32_t variant_idx_offset = 0;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_variant_ct = g_cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    double* beta_se_iter = &(g_block_beta_se[2 * max_reported_test_ct * variant_bidx]);
    logistic_aux_result_t* block_aux_iter = &(g_logistic_block_aux[variant_bidx]);
    const float* local_covars_iter = nullptr;
    if (local_covar_ct) {
      // &(nullptr[0]) is okay in C++, but undefined in C
      local_covars_iter = &(g_local_covars_vcmaj_f[parity][variant_bidx * max_sample_ct * local_covar_ct]);
    }
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = uint32arr_greater_than(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
      const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      uint32_t cur_variant_bidx_end = subset_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
      if (cur_variant_bidx_end > variant_bidx_end) {
	cur_variant_bidx_end = variant_bidx_end;
      }
      const uint32_t is_x = (chr_idx == x_code);
      const uint32_t is_y = (chr_idx == y_code);
      const uint32_t is_nonx_haploid = (!is_x) && is_set(cip->haploid_mask, chr_idx);
      const uintptr_t* cur_sample_include;
      const uint32_t* cur_sample_include_cumulative_popcounts;
      const uintptr_t* cur_pheno_cc;
      const float* cur_pheno;
      const float* cur_covars_cmaj;
      const uintptr_t* cur_parameter_subset;
      const float* cur_constraints_con_major;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      uint32_t cur_constraint_ct;
      uint32_t primary_pred_idx = include_intercept;
      if (is_y && g_sample_include_y) {
	cur_sample_include = g_sample_include_y;
	cur_sample_include_cumulative_popcounts = g_sample_include_y_cumulative_popcounts;
	cur_pheno_cc = g_pheno_y_cc;
	cur_pheno = g_pheno_y_f;
	cur_covars_cmaj = g_covars_cmaj_y_f;
	cur_parameter_subset = g_parameter_subset_y;
	cur_constraints_con_major = g_constraints_con_major_y_f;
	cur_sample_ct = g_sample_ct_y;
	cur_covar_ct = g_covar_ct_y;
	cur_constraint_ct = g_constraint_ct_y;
      } else if (is_x && g_sample_include_x) {
	cur_sample_include = g_sample_include_x;
	cur_sample_include_cumulative_popcounts = g_sample_include_x_cumulative_popcounts;
	cur_pheno_cc = g_pheno_x_cc;
	cur_pheno = g_pheno_x_f;
	cur_covars_cmaj = g_covars_cmaj_x_f;
	cur_parameter_subset = g_parameter_subset_x;
	cur_constraints_con_major = g_constraints_con_major_x_f;
	cur_sample_ct = g_sample_ct_x;
	cur_covar_ct = g_covar_ct_x;
	cur_constraint_ct = g_constraint_ct_x;
      } else {
	cur_sample_include = g_sample_include;
	cur_sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
	cur_pheno_cc = g_pheno_cc;
	cur_pheno = g_pheno_f;
	cur_covars_cmaj = g_covars_cmaj_f;
	cur_parameter_subset = g_parameter_subset;
	cur_constraints_con_major = g_constraints_con_major_f;
	cur_sample_ct = g_sample_ct;
	cur_covar_ct = g_covar_ct;
	cur_constraint_ct = g_constraint_ct;
      }
      const uint32_t sample_ctl = BITCT_TO_WORDCT(cur_sample_ct);
      const uint32_t sample_cta4 = round_up_pow2(cur_sample_ct, 4);
      const uint32_t cur_predictor_ct_base = 2 + domdev_present + cur_covar_ct * (1 + add_interactions * domdev_present_p1);
      uint32_t cur_predictor_ct = cur_predictor_ct_base;
      if (cur_parameter_subset) {
	cur_predictor_ct = popcount_longs(cur_parameter_subset, BITCT_TO_WORDCT(cur_predictor_ct_base));
      }
      const uint32_t predictor_cta4 = round_up_pow2(cur_predictor_ct, 4);
      const uint32_t predictor_cta4p1 = predictor_cta4 + 1;
      uint32_t reported_pred_uidx_end;
      if (hide_covar) {
	if (!cur_parameter_subset) {
	  reported_pred_uidx_end = 2 + domdev_present;
	} else {
	  reported_pred_uidx_end = 1 + is_set(cur_parameter_subset, 1) + domdev_present;
	}
      } else {
	reported_pred_uidx_end = cur_predictor_ct;
      }
      // todo: --tests
      if (cur_constraint_ct) {
	primary_pred_idx = reported_pred_uidx_end - reported_pred_uidx_start;
      }
      const uint32_t genof_buffer_needed = cur_parameter_subset && (!is_set(cur_parameter_subset, 1));
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = (uintptr_t*)arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter);
      uintptr_t* pheno_cc_nm = (uintptr_t*)arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter);
      uintptr_t* male_nm = (uintptr_t*)arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter);
      float* nm_pheno_buf = (float*)arena_alloc_raw_rd(sample_cta4 * sizeof(float), &workspace_iter);
      float* nm_predictors_pmaj_buf = (float*)arena_alloc_raw_rd((cur_predictor_ct + genof_buffer_needed) * sample_cta4 * sizeof(float), &workspace_iter);
      float* coef_return = (float*)arena_alloc_raw_rd(predictor_cta4 * sizeof(float), &workspace_iter);
      float* hh_return = (float*)arena_alloc_raw_rd(cur_predictor_ct * predictor_cta4 * sizeof(float), &workspace_iter);
      float* pp_buf = (float*)arena_alloc_raw_rd(sample_cta4 * sizeof(float), &workspace_iter);
      float* sample_variance_buf = (float*)arena_alloc_raw_rd(sample_cta4 * sizeof(float), &workspace_iter);
      float* gradient_buf = (float*)arena_alloc_raw_rd(predictor_cta4 * sizeof(float), &workspace_iter);
      float* dcoef_buf = (float*)arena_alloc_raw_rd(predictor_cta4 * sizeof(float), &workspace_iter);
      float* cholesky_decomp_return = (float*)arena_alloc_raw_rd(cur_predictor_ct * predictor_cta4 * sizeof(float), &workspace_iter);
      
      matrix_finvert_buf1_t* inv_1d_buf = nullptr;
      float* flt_2d_buf = nullptr;

      // Firth-only
      float* score_buf = nullptr;
      float* tmpnxk_buf = nullptr;

      // joint test only
      float* tmphxs_buf = nullptr;
      float* h_transpose_buf = nullptr;
      float* inner_buf = nullptr;

      if (is_sometimes_firth || cur_constraint_ct) {
	inv_1d_buf = (matrix_finvert_buf1_t*)arena_alloc_raw_rd(cur_predictor_ct * kMatrixFinvertBuf1CheckedAlloc, &workspace_iter);
	flt_2d_buf = (float*)arena_alloc_raw_rd(cur_predictor_ct * predictor_cta4 * sizeof(float), &workspace_iter);
	if (is_sometimes_firth) {
	  score_buf = (float*)arena_alloc_raw_rd(sample_cta4 * sizeof(float), &workspace_iter);
	  tmpnxk_buf = (float*)arena_alloc_raw_rd(cur_predictor_ct * sample_cta4 * sizeof(float), &workspace_iter);
	}
	if (cur_constraint_ct) {
	  tmphxs_buf = (float*)arena_alloc_raw_rd(cur_constraint_ct * predictor_cta4 * sizeof(float), &workspace_iter);
	  h_transpose_buf = (float*)arena_alloc_raw_rd(cur_constraint_ct * predictor_cta4 * sizeof(float), &workspace_iter);
	  inner_buf = (float*)arena_alloc_raw_rd(cur_constraint_ct * cur_constraint_ct * sizeof(float), &workspace_iter);
	}
      }
      // assert((uintptr_t)(workspace_iter - workspace_buf) == get_logistic_workspace_size(cur_sample_ct, cur_predictor_ct, cur_constraint_ct, genof_buffer_needed, is_sometimes_firth));
      pgr_clear_ld_cache(pgrp);
      uint32_t genocounts[4];
      for (; variant_bidx < cur_variant_bidx_end; ++variant_bidx, ++variant_uidx) {
	next_set_unsafe_ck(variant_include, &variant_uidx);
	{
	  uint32_t dosage_ct;
	  uint32_t is_explicit_alt1;
	  pglerr_t reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(cur_sample_include, cur_sample_include_cumulative_popcounts, cur_sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
	  if (reterr) {
	    g_error_ret = reterr;
	    variant_bidx = variant_bidx_end;
	    break;
	  }
	  zero_trailing_quaters(cur_sample_ct, genovec);
	  genovec_count_freqs_unsafe(genovec, cur_sample_ct, genocounts);
	  uint32_t missing_ct = genocounts[3];
	  if (!missing_ct) {
	    fill_all_bits(cur_sample_ct, sample_nm);
	  } else {
	    genoarr_to_nonmissing(genovec, cur_sample_ct, sample_nm);
	    if (dosage_ct) {
	      bitvec_or(dosage_present, sample_ctl, sample_nm);
	      missing_ct = cur_sample_ct - popcount_longs(sample_nm, sample_ctl);
	    }
	  }
	  uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
	  // todo: alt2/alt3/etc. dosage > 0.5 -> missing
	  const uint32_t nm_sample_ctl = BITCT_TO_WORDCT(nm_sample_ct);
	  const uint32_t nm_sample_cta4 = round_up_pow2(nm_sample_ct, 4);
	  const uint32_t nm_sample_ct_rem = nm_sample_cta4 - nm_sample_ct;
	  float* nm_predictors_pmaj_iter = nm_predictors_pmaj_buf;
	  // first predictor column: intercept
	  for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	    *nm_predictors_pmaj_iter++ = 1.0;
	  }
	  fill_float_zero(nm_sample_ct_rem, nm_predictors_pmaj_iter);
	  // second predictor column: genotype
	  float* genotype_vals = &(nm_predictors_pmaj_buf[nm_sample_cta4]);
	  if (genof_buffer_needed) {
	    // special case: --parameters excludes the main genotype column,
	    // but does care about an interaction
	    genotype_vals = &(nm_predictors_pmaj_buf[cur_predictor_ct * nm_sample_cta4]);
	  }
	  nm_predictors_pmaj_iter = genotype_vals;
	  if (!missing_ct) {
	    genoarr_to_floats(genovec, nm_sample_ct, nm_predictors_pmaj_iter);
	    if (dosage_ct) {
	      uint32_t sample_idx = 0;
	      for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_idx) {
		next_set_unsafe_ck(dosage_present, &sample_idx);
		// 32768 -> 2, 16384 -> 1, 0 -> 0
		nm_predictors_pmaj_iter[sample_idx] = kRecipDosageMidf * ((int32_t)((uint32_t)dosage_vals[dosage_idx]));
	      }
	    }
	  } else {
	    if (!dosage_ct) {
	      genoarr_to_floats_remove_missing(genovec, cur_sample_ct, nm_predictors_pmaj_iter);
	    } else {
	      uint32_t sample_midx = 0;
	      uint32_t dosage_idx = 0;
	      for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
		next_set_unsafe_ck(sample_nm, &sample_midx);
		float cur_val;
		if (is_set(dosage_present, sample_midx)) {
		  cur_val = kRecipDosageMidf * ((int32_t)((uint32_t)dosage_vals[dosage_idx++]));
		} else {
		  cur_val = (intptr_t)(GET_QUATERARR_ENTRY(genovec, sample_midx));
		}
	        nm_predictors_pmaj_iter[sample_idx] = cur_val;
	      }
	    }
	  }
	  nm_predictors_pmaj_iter = &(nm_predictors_pmaj_iter[nm_sample_ct]);
	  append_float_zero(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
	  copy_bitarr_subset(cur_pheno_cc, sample_nm, nm_sample_ct, pheno_cc_nm);
	  const uint32_t nm_case_ct = popcount_longs(pheno_cc_nm, nm_sample_ctl);
	  // usually need to save some of {sample_obs_ct, allele_obs_ct,
	  // alt_dosage, case_allele_obs_ct, alt_case_dosage, mach_r2 even
	  // for skipped variants
	  // compute them all for now, could conditionally skip later
	  block_aux_iter->sample_obs_ct = nm_sample_ct;
	  double dosage_ceil = 2.0;
	  if (!is_x) {
	    if (!is_nonx_haploid) {
	      block_aux_iter->allele_obs_ct = nm_sample_ct * 2;
	      block_aux_iter->case_allele_obs_ct = nm_case_ct * 2;
	    } else {
	      block_aux_iter->allele_obs_ct = nm_sample_ct;
	      block_aux_iter->case_allele_obs_ct = nm_case_ct;
	      // everything is on 0..1 scale, not 0..2
	      dosage_ceil = 1.0;
	      for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
		genotype_vals[sample_idx] *= 0.5;
	      }
	    }
	  } else {
	    copy_bitarr_subset(sex_male_collapsed, sample_nm, nm_sample_ct, male_nm);
	    const uint32_t nm_male_ct = popcount_longs(male_nm, nm_sample_ctl);
	    block_aux_iter->allele_obs_ct = nm_sample_ct * 2;
	    block_aux_iter->case_allele_obs_ct = nm_case_ct * 2;
	    if (is_xchr_model_1) {
	      // special case: multiply male values by 0.5
	      uint32_t sample_idx = 0;
	      for (uint32_t male_idx = 0; male_idx < nm_male_ct; ++male_idx, ++sample_idx) {
		next_set_unsafe_ck(male_nm, &sample_idx);
		genotype_vals[sample_idx] *= 0.5;
	      }
	      block_aux_iter->allele_obs_ct -= nm_male_ct;
	      block_aux_iter->case_allele_obs_ct -= popcount_longs_intersect(pheno_cc_nm, male_nm, nm_sample_ctl);
	    }
	  }
	  double alt_case_dosage = 0.0;
	  double dosage_sum = 0.0;
	  // genotype_vals restricted to [0, 2], so naive variance computation
	  // is stable
	  double dosage_ssq = 0.0;
	  for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	    const double cur_genotype_val = genotype_vals[sample_idx];
	    dosage_sum += cur_genotype_val;
	    dosage_ssq += cur_genotype_val * cur_genotype_val;
	    alt_case_dosage += cur_genotype_val * ((double)((int32_t)is_set(pheno_cc_nm, sample_idx)));
	  }
	  block_aux_iter->firth_fallback = 0;
	  block_aux_iter->alt_dosage = dosage_sum;
	  block_aux_iter->alt_case_dosage = alt_case_dosage;

	  const double dosage_avg = dosage_sum / ((double)((int32_t)nm_sample_ct));
	  const double dosage_variance = dosage_ssq - dosage_sum * dosage_avg;
	  // note that this value is nonsense on chrX/chrY/MT/haploid
	  block_aux_iter->mach_r2 = 2 * dosage_variance / (dosage_sum * (dosage_ceil - dosage_avg));
	  // okay, now we're free to skip the actual regression if there are
	  // too few samples, or remaining samples are all-case/all-control, or
	  // variant is monomorphic (or all-het)
	  if ((nm_sample_ct < cur_predictor_ct) || (!nm_case_ct) || (nm_case_ct == nm_sample_ct) || (fabs(dosage_variance) < kBigEpsilon)) {
	    goto glm_logistic_thread_skip_variant;
	  }
	  float* domdev_vals = nullptr;
	  if (genof_buffer_needed) {
	    nm_predictors_pmaj_iter = &(nm_predictors_pmaj_buf[nm_sample_cta4]);
	  } else if (joint_genotypic || joint_hethom) {
	    // in hethom case, do this before clobbering genotype data
	    domdev_vals = nm_predictors_pmaj_iter;
	    for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	      float cur_genotype_val = genotype_vals[sample_idx];
	      if (cur_genotype_val > 1.0) {
		cur_genotype_val = 2.0 - cur_genotype_val;
	      }
	      nm_predictors_pmaj_iter[sample_idx] = cur_genotype_val;
	    }
	    nm_predictors_pmaj_iter = &(nm_predictors_pmaj_iter[nm_sample_ct]);
	    append_float_zero(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
	  }
	  if (model_dominant) {
	    for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	      const float cur_genotype_val = genotype_vals[sample_idx];
	      // 0..1..1
	      if (cur_genotype_val > 1.0) {
		genotype_vals[sample_idx] = 1.0;
	      }
	    }
	  } else if (model_recessive || joint_hethom) {
	    for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	      const float cur_genotype_val = genotype_vals[sample_idx];
	      // 0..0..1
	      if (cur_genotype_val < 1.0) {
		genotype_vals[sample_idx] = 0.0;
	      } else {
		genotype_vals[sample_idx] = cur_genotype_val - 1.0;
	      }
	    }
	  }

	  // fill phenotype
	  uint32_t sample_midx = 0;
	  for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
	    next_set_unsafe_ck(sample_nm, &sample_midx);
	    nm_pheno_buf[sample_idx] = cur_pheno[sample_midx];
	  }
	  
	  // fill covariates
	  uint32_t parameter_uidx = 2 + domdev_present;
	  for (uint32_t covar_idx = 0; covar_idx < cur_covar_ct; ++covar_idx, ++parameter_uidx) {
	    // strictly speaking, we don't need cur_covars_cmaj to be
	    // vector-aligned
	    if (cur_parameter_subset && (!is_set(cur_parameter_subset, parameter_uidx))) {
	      continue;
	    }
	    const float* cur_covar_col;
	    if (covar_idx < local_covar_ct) {
	      cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
	    } else {
	      cur_covar_col = &(cur_covars_cmaj[(covar_idx - local_covar_ct) * sample_cta4]);
	    }
	    sample_midx = 0;
	    for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
	      next_set_unsafe_ck(sample_nm, &sample_midx);
	      *nm_predictors_pmaj_iter++ = cur_covar_col[sample_midx];
	    }
	    append_float_zero(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
	  }
	  // fill interaction terms
	  if (add_interactions) {
	    for (uint32_t covar_idx = 0; covar_idx < cur_covar_ct; ++covar_idx) {
	      const float* cur_covar_col;
	      if (covar_idx < local_covar_ct) {
	        cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
	      } else {
		cur_covar_col = &(cur_covars_cmaj[covar_idx * sample_cta4]);
	      }
	      if ((!cur_parameter_subset) || is_set(cur_parameter_subset, parameter_uidx)) {
		sample_midx = 0;
		for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
		  next_set_unsafe_ck(sample_nm, &sample_midx);
		  *nm_predictors_pmaj_iter++ = genotype_vals[sample_idx] * cur_covar_col[sample_midx];
		}
		append_float_zero(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
	      }
	      ++parameter_uidx;
	      if (domdev_present) {
		if ((!cur_parameter_subset) || is_set(cur_parameter_subset, parameter_uidx)) {
		  sample_midx = 0;
		  for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
		    next_set_unsafe_ck(sample_nm, &sample_midx);
		    *nm_predictors_pmaj_iter++ = domdev_vals[sample_idx] * cur_covar_col[sample_midx];
		  }
		  append_float_zero(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
		}
		++parameter_uidx;
	      }
	    }
	  }
	  fill_float_zero(predictor_cta4, coef_return);
	  if (!is_always_firth) {
	    const double ref_plus_alt1_sum = dosage_ceil * ((int32_t)nm_sample_ct);
	    // "dosage_sum" = alt1 sum
	    const double ref_sum = ref_plus_alt1_sum - dosage_sum;
	    const double ref_case_dosage = dosage_ceil * ((int32_t)nm_case_ct) - alt_case_dosage;
	    const double alt1_ctrl_dosage = dosage_sum - alt_case_dosage;
	    const double ref_ctrl_dosage = ref_sum - ref_case_dosage;
	    if ((alt_case_dosage == 0.0) || (fabs(ref_case_dosage) < kBigEpsilon) || (fabs(alt1_ctrl_dosage) < kBigEpsilon) || (fabs(ref_ctrl_dosage) < kBigEpsilon)) {
	      if (is_sometimes_firth) {
		block_aux_iter->firth_fallback = 1;
		goto glm_logistic_thread_firth_fallback;
	      } else {
		// this fails to converge >99.99% of the time, but better to
		// explicitly detect it since that can siginificantly speed
		// things up
		goto glm_logistic_thread_skip_variant;
	      }
	    }
	    if (logistic_regression(nm_pheno_buf, nm_predictors_pmaj_buf, nm_sample_ct, cur_predictor_ct, coef_return, cholesky_decomp_return, pp_buf, sample_variance_buf, hh_return, gradient_buf, dcoef_buf)) {
	      if (is_sometimes_firth) {
		fill_float_zero(predictor_cta4, coef_return);
		block_aux_iter->firth_fallback = 1;
		goto glm_logistic_thread_firth_fallback;
	      }
	      goto glm_logistic_thread_skip_variant;
	    }
	    // unlike firth_regression(), hh_return isn't inverted yet, do that
	    // here
	    for (uint32_t pred_uidx = 0; pred_uidx < cur_predictor_ct; ++pred_uidx) {
	      float* hh_inv_row = &(hh_return[pred_uidx * predictor_cta4]);
	      // fill_float_zero(cur_predictor_ct, gradient_buf);
	      // gradient_buf[pred_uidx] = 1.0;
	      // (y is gradient_buf, x is dcoef_buf)
	      // solve_linear_system(cholesky_decomp_return, gradient_buf, cur_predictor_ct, hh_inv_row);
	      // that works, but doesn't exploit the sparsity of y

	      // hh_return does now have vector-aligned rows
	      fill_float_zero(pred_uidx, hh_inv_row);

	      float fxx = 1.0;
	      for (uint32_t row_idx = pred_uidx; row_idx < cur_predictor_ct; ++row_idx) {
		const float* ll_row = &(cholesky_decomp_return[row_idx * predictor_cta4]);
		for (uint32_t col_idx = pred_uidx; col_idx < row_idx; ++col_idx) {
		  fxx -= ll_row[col_idx] * hh_inv_row[col_idx];
		}
		hh_inv_row[row_idx] = fxx / ll_row[row_idx];
		fxx = 0.0;
	      }
	      for (uint32_t col_idx = cur_predictor_ct; col_idx; ) {
		fxx = hh_inv_row[--col_idx];
		float* hh_inv_row_iter = &(hh_inv_row[cur_predictor_ct - 1]);
		for (uint32_t row_idx = cur_predictor_ct - 1; row_idx > col_idx; --row_idx) {
		  fxx -= cholesky_decomp_return[row_idx * predictor_cta4 + col_idx] * (*hh_inv_row_iter--);
		}
		*hh_inv_row_iter = fxx / cholesky_decomp_return[col_idx * predictor_cta4p1];
	      }
	    }
	  } else {
	  glm_logistic_thread_firth_fallback:
	    if (firth_regression(nm_pheno_buf, nm_predictors_pmaj_buf, nm_sample_ct, cur_predictor_ct, coef_return, hh_return, inv_1d_buf, flt_2d_buf, pp_buf, sample_variance_buf, gradient_buf, dcoef_buf, score_buf, tmpnxk_buf)) {
	      goto glm_logistic_thread_skip_variant;
	    }
	  }
	  // validParameters() check
	  for (uint32_t pred_uidx = 1; pred_uidx < cur_predictor_ct; ++pred_uidx) {
	    const float hh_inv_diag_element = hh_return[pred_uidx * predictor_cta4p1];
	    if ((hh_inv_diag_element < 1e-20) || (!realnum(hh_inv_diag_element))) {
	      goto glm_logistic_thread_skip_variant;
	    }
	    // use sample_variance_buf[] to store diagonal square roots
	    sample_variance_buf[pred_uidx] = sqrtf(hh_inv_diag_element);
	  }
	  sample_variance_buf[0] = sqrtf(hh_return[0]);
	  for (uint32_t pred_uidx = 1; pred_uidx < cur_predictor_ct; ++pred_uidx) {
	    const float cur_hh_inv_diag_sqrt = 0.99999 * sample_variance_buf[pred_uidx];
	    const float* hh_inv_row_iter = &(hh_return[pred_uidx * predictor_cta4]);
	    const float* hh_inv_diag_sqrts_iter = sample_variance_buf;
	    for (uint32_t pred_uidx2 = 0; pred_uidx2 < pred_uidx; ++pred_uidx2) {
	      if ((*hh_inv_row_iter++) > cur_hh_inv_diag_sqrt * (*hh_inv_diag_sqrts_iter++)) {
		goto glm_logistic_thread_skip_variant;
	      }
	    }
	  }
	  double* beta_se_iter2 = beta_se_iter;
	  for (uint32_t pred_uidx = reported_pred_uidx_start; pred_uidx < reported_pred_uidx_end; ++pred_uidx) {
	    *beta_se_iter2++ = coef_return[pred_uidx];
	    *beta_se_iter2++ = (double)sample_variance_buf[pred_uidx];
	  }
	  if (cur_constraint_ct) {
	    *beta_se_iter2++ = 0.0;
	    double chisq;
	    if (!linear_hypothesis_chisq_f(coef_return, cur_constraints_con_major, hh_return, cur_constraint_ct, cur_predictor_ct, predictor_cta4, &chisq, tmphxs_buf, h_transpose_buf, inner_buf, inv_1d_buf, flt_2d_buf)) {
	      *beta_se_iter2++ = chisq;
	    } else {
	      *beta_se_iter2++ = -9;
	    }
	  }
	}
	while (0) {
	glm_logistic_thread_skip_variant:
	  beta_se_iter[primary_pred_idx * 2 + 1] = -9;
	}
	beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
	++block_aux_iter;
	if (local_covars_iter) {
	  local_covars_iter = &(local_covars_iter[local_covar_ct * max_sample_ct]);
	}
	// todo?
      }
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
    variant_idx_offset += cur_block_variant_ct;
  }
}

uint32_t get_reported_test_ct(const uintptr_t* parameter_subset, glm_flags_t glm_flags, uint32_t covar_ct) {
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
  // TODO: --tests
  const uint32_t joint_test = domdev_present;

  if (hide_covar) {
    if (!parameter_subset) {
      return 1 + include_intercept + domdev_present + joint_test;
    }
    return include_intercept + domdev_present + joint_test + is_set(parameter_subset, 1);
  }
  
  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t predictor_ct_base = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
  uint32_t predictor_ct = predictor_ct_base;
  if (parameter_subset) {
    predictor_ct = popcount_longs(parameter_subset, BITCT_TO_WORDCT(predictor_ct_base));
  }
  return predictor_ct + joint_test + include_intercept - 1;
}

boolerr_t alloc_and_init_reported_test_names(const uintptr_t* parameter_subset, char** covar_names, glm_flags_t glm_flags, uint32_t covar_ct, char*** cur_test_names_ptr) {
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t is_hethom = (glm_flags / kfGlmHethom) & 1;
  const uint32_t domdev_present = (glm_flags & kfGlmGenotypic) || is_hethom;
  char main_effect[4];
  if (model_dominant) {
    memcpy(main_effect, "DOMx", 4);
  } else if (model_recessive) {
    memcpy(main_effect, "RECx", 4);
  } else if (is_hethom) {
    memcpy(main_effect, "HOMx", 4);
  } else {
    memcpy(main_effect, "ADDx", 4);
  }
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t include_main_effect = (!parameter_subset) || is_set(parameter_subset, 1);
  char domdev_str[8];
  uint32_t domdev_slen = 7;
  if (!is_hethom) {
    strcpy(domdev_str, "DOMDEVx");
  } else {
    strcpy(domdev_str, "HETx");
    domdev_slen = 4;
  }
  
  // TODO: --tests
  const uint32_t joint_test = domdev_present;
  
  if (glm_flags & kfGlmHideCovar) {
    const uint32_t reported_test_ct = include_intercept + include_main_effect + domdev_present + joint_test;
    char* test_name_buf_iter;
    if (bigstack_alloc_cp(reported_test_ct, cur_test_names_ptr) ||
	bigstack_alloc_c(64, &test_name_buf_iter)) {
      return 1;
    }
    char** cur_test_names = *cur_test_names_ptr;
    uint32_t write_idx = 0;
    if (include_intercept) {
      char* iter_next = memcpya(test_name_buf_iter, "INTERCEPT", 10);
      cur_test_names[write_idx++] = test_name_buf_iter;
      test_name_buf_iter = iter_next;
    }
    if (include_main_effect) {
      char* iter_next = memcpyax(test_name_buf_iter, main_effect, 3, '\0');
      cur_test_names[write_idx++] = test_name_buf_iter;
      test_name_buf_iter = iter_next;
    }
    if (domdev_present) {
      char* iter_next = memcpyax(test_name_buf_iter, domdev_str, domdev_slen - 1, '\0');
      cur_test_names[write_idx++] = test_name_buf_iter;
      test_name_buf_iter = iter_next;
    }
    if (joint_test) {
      // TODO: --tests
      strcpy(test_name_buf_iter, "GENO_2DF");
      cur_test_names[write_idx++] = test_name_buf_iter;
    }
    assert(write_idx == reported_test_ct);
    return 0;
  }
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  uint32_t predictor_ct_base = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
  uint32_t predictor_ct = predictor_ct_base;
  if (parameter_subset) {
    predictor_ct = popcount_longs(parameter_subset, BITCT_TO_WORDCT(predictor_ct_base));
  }
  const uint32_t reported_test_ct = predictor_ct + joint_test + include_intercept - 1;
  uintptr_t test_name_buf_alloc = 64;
  if (add_interactions) {
    // don't bother optimizing this for parameter_subset case for now
    uintptr_t covar_name_total_blen = covar_ct;
    for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx) {
      covar_name_total_blen += strlen(covar_names[covar_idx]);
    }
    // ADDx[covar name], etc.
    test_name_buf_alloc += 4 * covar_ct + covar_name_total_blen;
    if (is_hethom) {
      // HETx
      test_name_buf_alloc += 4 * covar_ct + covar_name_total_blen;
    } else if (domdev_present) {
      // DOMDEVx
      test_name_buf_alloc += 7 * covar_ct + covar_name_total_blen;
    }
  }
  char* test_name_buf_iter;
  if (bigstack_alloc_cp(reported_test_ct, cur_test_names_ptr) ||
      bigstack_alloc_c(test_name_buf_alloc, &test_name_buf_iter)) {
    return 1;
  }
  char** cur_test_names = *cur_test_names_ptr;
  uint32_t write_idx = 0;
  if (include_intercept) {
    char* iter_next = memcpya(test_name_buf_iter, "INTERCEPT", 10);
    cur_test_names[write_idx++] = test_name_buf_iter;
    test_name_buf_iter = iter_next;
  }
  if (include_main_effect) {
    char* iter_next = memcpyax(test_name_buf_iter, main_effect, 3, '\0');
    cur_test_names[write_idx++] = test_name_buf_iter;
    test_name_buf_iter = iter_next;
  }
  if (domdev_present) {
    char* iter_next = memcpyax(test_name_buf_iter, domdev_str, domdev_slen - 1, '\0');
    cur_test_names[write_idx++] = test_name_buf_iter;
    test_name_buf_iter = iter_next;
  }
  uint32_t pred_uidx = 2 + domdev_present;
  for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx, ++pred_uidx) {
    if (parameter_subset && (!is_set(parameter_subset, pred_uidx))) {
      continue;
    }
    // just point to the existing string, its lifetime is sufficient
    cur_test_names[write_idx++] = covar_names[covar_idx];
  }
  if (add_interactions) {
    for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx) {
      const char* cur_covar_name = covar_names[covar_idx];
      if ((!parameter_subset) || is_set(parameter_subset, pred_uidx)) {
	char* iter_next = memcpya(test_name_buf_iter, main_effect, 4);
	iter_next = strcpyax(iter_next, cur_covar_name, '\0');
	cur_test_names[write_idx++] = test_name_buf_iter;
	test_name_buf_iter = iter_next;
      }
      ++pred_uidx;
      if (domdev_present) {
	if ((!parameter_subset) || is_set(parameter_subset, pred_uidx)) {
	  char* iter_next = memcpya(test_name_buf_iter, domdev_str, domdev_slen);
	  iter_next = strcpyax(iter_next, cur_covar_name, '\0');
	  cur_test_names[write_idx++] = test_name_buf_iter;
	  test_name_buf_iter = iter_next;
	}
	++pred_uidx;
      }
    }
  }
  if (joint_test) {
    // todo: --tests
    strcpy(test_name_buf_iter, "GENO_2DF");
    cur_test_names[write_idx++] = test_name_buf_iter;
  }
  assert(write_idx == reported_test_ct);
  return 0;
}

boolerr_t alloc_and_init_constraints_f(uint32_t predictor_ct, uint32_t* constraint_ct_ptr, float** constraints_con_major_f_ptr) {
  // todo: --tests
  const uint32_t constraint_ct = 2;
  if (bigstack_calloc_f(constraint_ct * predictor_ct, constraints_con_major_f_ptr)) {
    return 1;
  }
  float* constraints_con_major_f = *constraints_con_major_f_ptr;
  constraints_con_major_f[1] = 1; // [0][1]
  constraints_con_major_f[predictor_ct + 2] = 1; // [1][2]
  *constraint_ct_ptr = constraint_ct;
  return 0;
}

boolerr_t alloc_and_init_constraints_d(uint32_t predictor_ct, uint32_t* constraint_ct_ptr, double** constraints_con_major_ptr) {
  const uint32_t constraint_ct = 2;
  if (bigstack_calloc_d(constraint_ct * predictor_ct, constraints_con_major_ptr)) {
    return 1;
  }
  double* constraints_con_major = *constraints_con_major_ptr;
  constraints_con_major[1] = 1; // [0][1]
  constraints_con_major[predictor_ct + 2] = 1; // [1][2]
  *constraint_ct_ptr = constraint_ct;
  return 0;
}

pglerr_t read_local_covar_block(const uintptr_t* sample_include, const uintptr_t* sample_include_x, const uintptr_t* sample_include_y, const uint32_t* sample_include_cumulative_popcounts, const uint32_t* sample_include_x_cumulative_popcounts, const uint32_t* sample_include_y_cumulative_popcounts, const chr_info_t* cip, const uintptr_t* variant_include, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, uint32_t sample_ct, uint32_t sample_ct_x, uint32_t sample_ct_y, uint32_t variant_uidx, uint32_t variant_uidx_end, uint32_t cur_block_variant_ct, uint32_t local_sample_ct, uint32_t local_covar_ct, uint32_t omit_last, uint32_t local_cat_ct, uint32_t local_loadbuf_size, gzFile gz_local_covar_file, uint32_t* local_line_idx_ptr, uint32_t* local_xy_ptr, float* local_covars_vcmaj_f_iter, double* local_covars_vcmaj_d_iter, uint32_t* local_sample_idx_order, char* local_loadbuf) {
  const int32_t x_code = cip->xymt_codes[kChrOffsetX];
  const int32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t tokens_per_sample = local_cat_ct? 1 : (local_covar_ct + omit_last);
  uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
  if (max_sample_ct < sample_ct_y) {
    max_sample_ct = sample_ct_y;
  }
  uint32_t variant_bidx = 0;
  if (local_cat_ct) {
    // assert(local_covar_ct == local_cat_ct - 1);
    if (local_covars_vcmaj_f_iter) {
      fill_float_zero(local_covar_ct * max_sample_ct * ((uintptr_t)cur_block_variant_ct), local_covars_vcmaj_f_iter);
    } else {
      fill_double_zero(local_covar_ct * max_sample_ct * ((uintptr_t)cur_block_variant_ct), local_covars_vcmaj_d_iter);
    }
  }
  uint32_t local_line_idx = *local_line_idx_ptr;
  while (variant_bidx < cur_block_variant_ct) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    const uint32_t chr_fo_idx = get_variant_chr_fo_idx(cip, variant_uidx);
    const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
    const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    uint32_t cur_variant_bidx_end = cur_block_variant_ct;
    if (chr_end < variant_uidx_end) {
      cur_variant_bidx_end = variant_bidx + popcount_bit_idx(variant_include, variant_uidx, chr_end);
      assert(cur_variant_bidx_end <= cur_block_variant_ct);
    }
    const uint32_t is_x = (chr_idx == x_code);
    const uint32_t is_y = (chr_idx == y_code);
    const uintptr_t* cur_sample_include;
    const uint32_t* cur_sample_include_cumulative_popcounts;
    uint32_t cur_sample_ct;
    if (is_y && sample_include_y) {
      cur_sample_include = sample_include_y;
      cur_sample_include_cumulative_popcounts = sample_include_y_cumulative_popcounts;
      cur_sample_ct = sample_ct_y;
    } else if (is_x && sample_include_x) {
      cur_sample_include = sample_include_x;
      cur_sample_include_cumulative_popcounts = sample_include_x_cumulative_popcounts;
      cur_sample_ct = sample_ct_x;
    } else {
      cur_sample_include = sample_include;
      cur_sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
      cur_sample_ct = sample_ct;
    }
    const uint32_t new_local_xy = is_x + 2 * is_y;
    if (new_local_xy != *local_xy_ptr) {
      for (uint32_t uii = 0; uii < local_sample_ct; ++uii) {
	const uint32_t cur_uidx = local_sample_uidx_order[uii];
	uint32_t cur_idx = 0xffffffffU;
	if ((cur_uidx != 0xffffffffU) && is_set(cur_sample_include, cur_uidx)) {
	  cur_idx = raw_to_subsetted_pos(cur_sample_include, cur_sample_include_cumulative_popcounts, cur_uidx);
	}
	local_sample_idx_order[uii] = cur_idx;
      }
      *local_xy_ptr = new_local_xy;
    }
    for (; variant_bidx < cur_variant_bidx_end; ++variant_bidx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      do {
        ++local_line_idx;
	if (!gzgets(gz_local_covar_file, local_loadbuf, local_loadbuf_size)) {
	  if (!gzeof(gz_local_covar_file)) {
	    return kPglRetReadFail;
	  }
	  logprint("\n");
	  logerrprint("Error: --glm local-covar= file has fewer lines than local-pvar= file.\n");
	  return kPglRetMalformedInput;
	}
	if (!local_loadbuf[local_loadbuf_size - 1]) {
	  logprint("\n");
	  LOGERRPRINTF("Error: Line %u of --glm local-covar= file is longer than expected.\n", local_line_idx);
	  return kPglRetMalformedInput;
	}
      } while (!is_set(local_variant_include, local_line_idx - 1));
      char* loadbuf_iter = skip_initial_spaces(local_loadbuf);
      uint32_t sample_idx = 0;
      for (uint32_t local_sample_idx = 0; sample_idx < cur_sample_ct; ++local_sample_idx) {
	const uint32_t cur_sample_idx = local_sample_idx_order[local_sample_idx];
	if (cur_sample_idx == 0xffffffffU) {
	  loadbuf_iter = next_token_mult(loadbuf_iter, tokens_per_sample);
	  if (!loadbuf_iter) {
	    logprint("\n");
	    LOGERRPRINTFWW("Error: Fewer tokens than expected on line %u of --glm local-covar= file.\n", local_line_idx);
	    return kPglRetMalformedInput;
	  }
	  continue;
	}
	if (local_cat_ct) {
	  uint32_t cat_idx;
	  if (scanadv_posint_capped(local_cat_ct, &loadbuf_iter, &cat_idx)) {
	    logprint("\n");
	    LOGERRPRINTF("Error: Invalid category index on line %u of --glm local-covar= file.\n", local_line_idx);
	    return kPglRetMalformedInput;
	  }
	  if (cat_idx != local_cat_ct) {
	    --cat_idx;
	    const uint32_t offset = cat_idx * max_sample_ct + cur_sample_idx;
	    if (local_covars_vcmaj_f_iter) {
	      local_covars_vcmaj_f_iter[offset] = 1.0;
	    } else {
	      local_covars_vcmaj_d_iter[offset] = 1.0;
	    }
	  }
	  while (!is_space_or_eoln(*loadbuf_iter)) {
	    ++loadbuf_iter;
	  }
	  loadbuf_iter = skip_initial_spaces(loadbuf_iter);
	} else {
	  if (local_covars_vcmaj_f_iter) {
	    float* local_covars_f_iter2 = &(local_covars_vcmaj_f_iter[cur_sample_idx]);
	    for (uint32_t covar_idx = 0; covar_idx < local_covar_ct; ++covar_idx) {
	      double dxx;
	      loadbuf_iter = scanadv_double(loadbuf_iter, &dxx);
	      if ((!loadbuf_iter) || (fabs(dxx) > 3.4028235677973362e38)) {
		logprint("\n");
		LOGERRPRINTF("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
		return kPglRetMalformedInput;
	      }
	      *local_covars_f_iter2 = (float)dxx;
	      local_covars_f_iter2 = &(local_covars_f_iter2[max_sample_ct]);
	      while (!is_space_or_eoln(*loadbuf_iter)) {
		++loadbuf_iter;
	      }
	      loadbuf_iter = skip_initial_spaces(loadbuf_iter);
	    }
	  } else {
	    double* local_covars_d_iter2 = &(local_covars_vcmaj_d_iter[cur_sample_idx]);
	    for (uint32_t covar_idx = 0; covar_idx < local_covar_ct; ++covar_idx) {
	      double dxx;
	      loadbuf_iter = scanadv_double(loadbuf_iter, &dxx);
	      if (!loadbuf_iter) {
		logprint("\n");
		LOGERRPRINTF("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
		return kPglRetMalformedInput;
	      }
	      *local_covars_d_iter2 = dxx;
	      local_covars_d_iter2 = &(local_covars_d_iter2[max_sample_ct]);
	      while (!is_space_or_eoln(*loadbuf_iter)) {
		++loadbuf_iter;
	      }
	      loadbuf_iter = skip_initial_spaces(loadbuf_iter);
	    }
	  }
	  if (omit_last) {
	    while (!is_space_or_eoln(*loadbuf_iter)) {
	      ++loadbuf_iter;
	    }
	    loadbuf_iter = skip_initial_spaces(loadbuf_iter);
	  }
	}
	++sample_idx;
      }
      if (local_covars_vcmaj_f_iter) {
	local_covars_vcmaj_f_iter += max_sample_ct * local_covar_ct;
      } else {
	local_covars_vcmaj_d_iter += max_sample_ct * local_covar_ct;
      }
    }
  }
  *local_line_idx_ptr = local_line_idx;
  return kPglRetSuccess;
}

// only pass the parameters which aren't also needed by the compute threads,
// for now
pglerr_t glm_logistic(const char* cur_pheno_name, char** test_names, char** test_names_x, char** test_names_y, const uint32_t* variant_bps, char** variant_ids, char** allele_storage, const glm_info_t* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const char* outname, uint32_t raw_variant_ct, uint32_t max_chr_blen, double ci_size, double pfilter, double output_min_p, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uint32_t local_sample_ct, uint32_t local_loadbuf_size, pgen_file_info_t* pgfip, gzFile gz_local_covar_file, uintptr_t* valid_variants, double* orig_pvals, double* orig_chisq, unsigned char* overflow_buf, char* local_loadbuf) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  threads_state_t ts;
  init_threads3z(&ts);
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uintptr_t* variant_include = g_variant_include;
    const chr_info_t* cip = g_cip;
    const uintptr_t* variant_allele_idxs = g_variant_allele_idxs;

    const uint32_t sample_ct = g_sample_ct;
    const uint32_t sample_ct_x = g_sample_ct_x;
    const uint32_t sample_ct_y = g_sample_ct_y;
    const uint32_t covar_ct = g_covar_ct;
    const uintptr_t local_covar_ct = g_local_covar_ct;
    const uint32_t covar_ct_x = g_covar_ct_x;
    const uint32_t covar_ct_y = g_covar_ct_y;

    uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
    if (max_sample_ct < sample_ct_y) {
      max_sample_ct = sample_ct_y;
    }
    uint32_t* local_sample_idx_order = nullptr;
    uint32_t local_line_idx = 0;
    uint32_t local_xy = 0; // 1 = chrX, 2 = chrY
    if (gz_local_covar_file) {
      if (gzrewind(gz_local_covar_file)) {
	goto glm_logistic_ret_READ_FAIL;
      }
      if (bigstack_alloc_ui(local_sample_ct, &local_sample_idx_order)) {
	goto glm_logistic_ret_NOMEM;
      }
      for (uint32_t uii = 0; uii < local_sample_ct; ++uii) {
	const uint32_t cur_uidx = local_sample_uidx_order[uii];
	uint32_t cur_idx = 0xffffffffU;
	if ((cur_uidx != 0xffffffffU) && is_set(g_sample_include, cur_uidx)) {
	  cur_idx = raw_to_subsetted_pos(g_sample_include, g_sample_include_cumulative_popcounts, cur_uidx);
	}
	local_sample_idx_order[uii] = cur_idx;
      }
    }
    
    const uint32_t variant_ct = g_variant_ct;
    
    const glm_flags_t glm_flags = glm_info_ptr->flags;    
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    // er, do not want to use multithreaded compression here... make sure to
    // add a forced-singlethreaded mode when multithreaded compression is
    // implemented.
    if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
      goto glm_logistic_ret_OPEN_FAIL;
    }
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t domdev_present_p1 = domdev_present + 1;

    // todo: --tests
    const uint32_t constraint_ct = g_constraint_ct;
    const uint32_t constraint_ct_x = g_constraint_ct_x;
    const uint32_t constraint_ct_y = g_constraint_ct_y;
    
    uint32_t predictor_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
    uint32_t predictor_ct_x = 2 + domdev_present + covar_ct_x * (1 + add_interactions * domdev_present_p1);
    uint32_t predictor_ct_y = 2 + domdev_present + covar_ct_y * (1 + add_interactions * domdev_present_p1);
    const uintptr_t* parameter_subset = g_parameter_subset;
    const uintptr_t* parameter_subset_x = g_parameter_subset_x;
    const uintptr_t* parameter_subset_y = g_parameter_subset_y;
    if (parameter_subset) {
      predictor_ct = popcount_longs(parameter_subset, BITCT_TO_WORDCT(predictor_ct));
      if (sample_ct_x) {
	predictor_ct_x = popcount_longs(parameter_subset_x, BITCT_TO_WORDCT(predictor_ct_x));
      } else {
	predictor_ct_x = 0;
      }
      if (sample_ct_y) {
	predictor_ct_y = popcount_longs(parameter_subset_y, BITCT_TO_WORDCT(predictor_ct_x));
      } else {
	predictor_ct_y = 0;
      }
    }
    uint32_t reported_test_ct = get_reported_test_ct(parameter_subset, glm_flags, covar_ct);
    uintptr_t max_reported_test_ct = reported_test_ct;
    uint32_t reported_test_ct_x = 0;
    if (sample_ct_x) {
      reported_test_ct_x = get_reported_test_ct(parameter_subset_x, glm_flags, covar_ct_x);
      if (reported_test_ct_x > max_reported_test_ct) {
	max_reported_test_ct = reported_test_ct_x;
      }
    }
    uint32_t reported_test_ct_y = 0;
    if (sample_ct_y) {
      reported_test_ct_y = get_reported_test_ct(parameter_subset_y, glm_flags, covar_ct_y);
      if (reported_test_ct_y > max_reported_test_ct) {
	max_reported_test_ct = reported_test_ct_y;
      }
    }
    const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
    const glm_cols_t glm_cols = glm_info_ptr->cols;
    const uint32_t test_col = glm_cols & kfGlmColTest;
    if ((!test_col) && (max_reported_test_ct > 1)) {
      logerrprint("Error: --glm's 'test' column cannot be omitted when results for multiple\npredictors are reported.  (Did you forget 'hide-covar'?)\n");
      goto glm_logistic_ret_INCONSISTENT_INPUT;
    }
    g_max_reported_test_ct = max_reported_test_ct;
    
    const uint32_t is_sometimes_firth = (glm_flags & (kfGlmFirthFallback | kfGlmFirth))? 1 : 0;
    const uint32_t is_always_firth = (glm_flags / kfGlmFirth) & 1;

    int32_t x_code = -2;
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    if (sample_ct_x) {
      get_xymt_code_start_and_end_unsafe(cip, kChrOffsetX, &x_code, &x_start, &x_end);
    }
    int32_t y_code = -2;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    if (sample_ct_y) {
      get_xymt_code_start_and_end_unsafe(cip, kChrOffsetY, &y_code, &y_start, &y_end);
    }
    const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t chr_col = glm_cols & kfGlmColChrom;

    // includes trailing tab
    char* chr_buf = nullptr;
    if (chr_col) {
      if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
	goto glm_logistic_ret_NOMEM;
      }
    }

    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > variant_ct) {
      calc_thread_ct = variant_ct;
    }

    const uint32_t genof_buffer_needed = parameter_subset && (!is_set(parameter_subset, 1));
    // workflow is similar to --make-bed
    uintptr_t workspace_alloc = get_logistic_workspace_size(sample_ct, predictor_ct, constraint_ct, genof_buffer_needed, is_sometimes_firth);
    if (sample_ct_x) {
      const uintptr_t workspace_alloc_x = get_logistic_workspace_size(sample_ct_x, predictor_ct_x, constraint_ct_x, genof_buffer_needed, is_sometimes_firth);
      if (workspace_alloc_x > workspace_alloc) {
	workspace_alloc = workspace_alloc_x;
      }
    }
    if (sample_ct_y) {
      const uintptr_t workspace_alloc_y = get_logistic_workspace_size(sample_ct_y, predictor_ct_y, constraint_ct_y, genof_buffer_needed, is_sometimes_firth);
      if (workspace_alloc_y > workspace_alloc) {
	workspace_alloc = workspace_alloc_y;
      }
    }
    // +1 is for top-level g_workspace_bufs
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;
    uintptr_t per_variant_xalloc_byte_ct = sizeof(logistic_aux_result_t) + 2 * max_reported_test_ct * sizeof(double) + max_sample_ct * local_covar_ct * sizeof(float);
    unsigned char* main_loadbufs[2];
    uint32_t read_block_size;
    if (multithread_load_init(variant_include, max_sample_ct, variant_ct, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, pgfip, &calc_thread_ct, &g_genovecs, dosage_is_present? (&g_dosage_presents) : nullptr, dosage_is_present? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &ts.threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto glm_logistic_ret_NOMEM;
    }
    ts.calc_thread_ct = calc_thread_ct;
    g_calc_thread_ct = calc_thread_ct;
    logistic_aux_result_t* logistic_block_aux_bufs[2];
    double* block_beta_se_bufs[2];
    
    for (uint32_t uii = 0; uii < 2; ++uii) {
      logistic_block_aux_bufs[uii] = (logistic_aux_result_t*)bigstack_alloc(read_block_size * sizeof(logistic_aux_result_t));
      if ((!logistic_block_aux_bufs[uii]) ||
	  bigstack_alloc_d(read_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii]))) {
	goto glm_logistic_ret_NOMEM;
      }
      if (local_covar_ct) {
	if (bigstack_alloc_f(read_block_size * max_sample_ct * local_covar_ct * sizeof(float), &(g_local_covars_vcmaj_f[uii]))) {
	  goto glm_logistic_ret_NOMEM;
	}
      } else {
	g_local_covars_vcmaj_f[uii] = nullptr;
      }
    }

    if (max_sample_ct > 2000000) {
      logerrprint("Warning: --glm logistic regression is unreliable on more than ~2 million\nsamples, since it uses single-precision arithmetic.\n");
    }
    g_workspace_bufs = (unsigned char**)bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t));
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_workspace_bufs[tidx] = bigstack_alloc_raw(workspace_alloc);
    }
    
    const uint32_t ref_col = glm_cols & kfGlmColRef;
    const uint32_t alt1_col = glm_cols & kfGlmColAlt1;
    const uint32_t alt_col = glm_cols & kfGlmColAlt;
    const uint32_t alt_ct_col = glm_cols & kfGlmColAltcount;
    const uint32_t tot_allele_col = glm_cols & kfGlmColTotallele;
    const uint32_t alt_ct_cc_col = glm_cols & kfGlmColAltcountcc;
    const uint32_t tot_allele_cc_col = glm_cols & kfGlmColTotallelecc;
    const uint32_t alt_freq_col = glm_cols & kfGlmColAltfreq;
    const uint32_t alt_freq_cc_col = glm_cols & kfGlmColAltfreqcc;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t firth_yn_col = (glm_cols & kfGlmColFirthYn) && is_sometimes_firth && (!is_always_firth);
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t orbeta_col = glm_cols & (kfGlmColBeta | kfGlmColOrbeta);
    const uint32_t report_beta_instead_of_odds_ratio = glm_cols & kfGlmColBeta;
    const uint32_t se_col = glm_cols & kfGlmColSe;
    const uint32_t ci_col = (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t t_col = glm_cols & kfGlmColT;
    const uint32_t p_col = glm_cols & kfGlmColP;
    cswritep = (char*)overflow_buf;
    *cswritep++ = '#';
    if (chr_col) {
      cswritep = strcpya(cswritep, "CHROM\t");
    }
    if (variant_bps) {
      cswritep = strcpya(cswritep, "POS\t");
    }
    cswritep = strcpya(cswritep, "ID");
    if (ref_col) {
      cswritep = strcpya(cswritep, "\tREF");
    }
    if (alt1_col) {
      cswritep = strcpya(cswritep, "\tALT1");
    }
    if (alt_col) {
      cswritep = strcpya(cswritep, "\tALT");
    }
    if (alt_ct_col) {
      cswritep = strcpya(cswritep, "\tALT_CT");
    }
    if (tot_allele_col) {
      cswritep = strcpya(cswritep, "\tALLELE_CT");
    }
    if (alt_ct_cc_col) {
      cswritep = strcpya(cswritep, "\tALT_CASE_CT\tALT_CTRL_CT");
    }
    if (tot_allele_cc_col) {
      cswritep = strcpya(cswritep, "\tCASE_ALLELE_CT\tCTRL_ALLELE_CT");
    }
    if (alt_freq_col) {
      cswritep = strcpya(cswritep, "\tALT_FREQ");
    }
    if (alt_freq_cc_col) {
      cswritep = strcpya(cswritep, "\tALT_CASE_FREQ\tALT_CTRL_FREQ");
    }
    if (mach_r2_col) {
      cswritep = strcpya(cswritep, "\tMACH_R2");
    }
    if (firth_yn_col) {
      cswritep = strcpya(cswritep, "\tFIRTH?");
    }
    if (test_col) {
      cswritep = strcpya(cswritep, "\tTEST");
    }
    if (nobs_col) {
      cswritep = strcpya(cswritep, "\tOBS_CT");
    }
    if (orbeta_col) {
      if (report_beta_instead_of_odds_ratio) {
	cswritep = strcpya(cswritep, "\tBETA");
      } else {
	cswritep = strcpya(cswritep, "\tOR");
      }
    }
    if (se_col) {
      cswritep = strcpya(cswritep, "\tSE");
    }
    double ci_zt = 0.0;
    if (ci_col) {
      cswritep = strcpya(cswritep, "\tL");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      cswritep = strcpya(cswritep, "\tU");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      ci_zt = ltqnorm((ci_size + 1.0) * 0.5);
    }
    if (t_col) {
      if (!constraint_ct) {
        cswritep = strcpya(cswritep, "\tT_STAT");
      } else {
	// chisq for joint tests.  may switch to F-statistic (just divide by
	// df; the hard part there is porting a function to convert that to a
	// p-value)
        cswritep = strcpya(cswritep, "\tT_OR_CHISQ_STAT");
      }
    }
    if (p_col) {
      cswritep = strcpya(cswritep, "\tP");
    }
    append_binary_eoln(&cswritep);

    // Main workflow:
    // 1. Set n=0, load/skip block 0
    //
    // 2. Spawn threads processing block n
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Load/skip block n unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8, Write results for last block
    const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
    const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t write_variant_uidx = 0;
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t suppress_mach_r2 = 0;

    // todo: --tests
    uint32_t cur_reported_test_ct = 0;
    uint32_t primary_reported_test_idx = include_intercept;
    uint32_t cur_constraint_ct = 0;

    char** cur_test_names = nullptr;
    uint32_t prev_block_variant_ct = 0;
    uint32_t variant_idx = 0;
    uint32_t cur_read_block_size = read_block_size;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t cur_allele_ct = 2;
    LOGPRINTFWW5("--glm %s regression on phenotype '%s': ", is_always_firth? "Firth" : (is_sometimes_firth? "logistic-Firth hybrid" : "logistic"), cur_pheno_name);
    fputs("0%", stdout);
    fflush(stdout);
    while (1) {
      uintptr_t cur_block_variant_ct = 0;
      if (!ts.is_last_block) {
	while (read_block_idx < read_block_ct_m1) {
	  cur_block_variant_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
	  if (cur_block_variant_ct) {
	    break;
	  }
	  ++read_block_idx;
	}
	if (read_block_idx == read_block_ct_m1) {
	  cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
	  cur_block_variant_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
	}
	if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_variant_ct, pgfip)) {
	  goto glm_logistic_ret_READ_FAIL;
	}
	if (gz_local_covar_file) {
	  reterr = read_local_covar_block(g_sample_include, g_sample_include_x, g_sample_include_y, g_sample_include_cumulative_popcounts, g_sample_include_x_cumulative_popcounts, g_sample_include_y_cumulative_popcounts, cip, variant_include, local_sample_uidx_order, local_variant_include, sample_ct, sample_ct_x, sample_ct_y, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_variant_ct, local_sample_ct, local_covar_ct, (glm_info_ptr->flags / kfGlmLocalOmitLast) & 1, glm_info_ptr->local_cat_ct, local_loadbuf_size, gz_local_covar_file, &local_line_idx, &local_xy, g_local_covars_vcmaj_f[parity], nullptr, local_sample_idx_order, local_loadbuf);
	  if (reterr) {
	    goto glm_logistic_ret_1;
	  }
	}
      }
      if (variant_idx) {
	join_threads3z(&ts);
	reterr = g_error_ret;
	if (reterr) {
	  if (reterr == kPglRetMalformedInput) {
	    logprint("\n");
	    logerrprint("Error: Malformed .pgen file.\n");
	  }
	  goto glm_logistic_ret_1;
	}
      }
      if (!ts.is_last_block) {
	g_cur_block_variant_ct = cur_block_variant_ct;
	const uint32_t uidx_start = read_block_idx * read_block_size;
	compute_uidx_start_partition(variant_include, cur_block_variant_ct, calc_thread_ct, uidx_start, g_read_variant_uidx_starts);
	for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	  g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
	  g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
	}
	g_logistic_block_aux = logistic_block_aux_bufs[parity];
	g_block_beta_se = block_beta_se_bufs[parity];
	ts.is_last_block = (variant_idx + cur_block_variant_ct == variant_ct);
	ts.thread_func_ptr = glm_logistic_thread;
	if (spawn_threads3z(variant_idx, &ts)) {
	  goto glm_logistic_ret_THREAD_CREATE_FAIL;
	}
      }
      parity = 1 - parity;
      if (variant_idx) {
	// write *previous* block results
	const double* cur_block_beta_se = block_beta_se_bufs[parity];
	const logistic_aux_result_t* cur_block_aux = logistic_block_aux_bufs[parity];
	const uint32_t variant_idx_start = variant_idx - prev_block_variant_ct;
	double* cur_pval_write = orig_pvals? (&(orig_pvals[variant_idx_start])) : nullptr;
	double* cur_chisq_write = orig_chisq? (&(orig_chisq[variant_idx_start])) : nullptr;
	for (uint32_t variant_bidx = 0; variant_bidx < prev_block_variant_ct; ++variant_bidx, ++write_variant_uidx) {
	  next_set_unsafe_ck(variant_include, &write_variant_uidx);
	  if (write_variant_uidx >= chr_end) {
	    do {
	      ++chr_fo_idx;
	      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	    } while (write_variant_uidx >= chr_end);
	    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	    suppress_mach_r2 = 1;
	    if ((chr_idx == ((uint32_t)x_code)) && sample_ct_x) {
	      cur_reported_test_ct = reported_test_ct_x;
	      cur_constraint_ct = constraint_ct_x;
	      cur_test_names = test_names_x;
	    } else if ((chr_idx == ((uint32_t)y_code)) && sample_ct_y) {
	      cur_reported_test_ct = reported_test_ct_y;
	      cur_constraint_ct = constraint_ct_y;
	      cur_test_names = test_names_y;
	    } else {
	      cur_reported_test_ct = reported_test_ct;
	      cur_constraint_ct = constraint_ct;
	      cur_test_names = test_names;
	      if ((chr_idx != ((uint32_t)x_code)) && (chr_idx != ((uint32_t)mt_code)) && (!is_set(cip->haploid_mask, chr_idx))) {
		suppress_mach_r2 = 0;
	      }
	    }
	    if (cur_constraint_ct) {
	      primary_reported_test_idx = reported_test_ct - 1;
	    }
	    if (chr_col) {
	      char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	      *chr_name_end = '\t';
	      chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
	    }
	  }
	  const double* beta_se_iter = &(cur_block_beta_se[2 * max_reported_test_ct * variant_bidx]);
	  const double primary_beta = beta_se_iter[primary_reported_test_idx * 2];
	  const double primary_se = beta_se_iter[primary_reported_test_idx * 2 + 1];
	  const uint32_t is_invalid = (primary_se == -9);
	  if (is_invalid && valid_variants) {
	    CLEAR_BIT(write_variant_uidx, valid_variants);
	  }
	  if (pfilter != 2.0) {
	    double primary_pval = 2.0;
	    if (!is_invalid) {
	      if (!cur_constraint_ct) {
		double primary_tstat = primary_beta / primary_se;
		// could precompute a tstat threshold instead
		primary_pval = chiprob_p(primary_tstat * primary_tstat, 1);
	      } else {
		// possible todo: support for F-distribution p-values instead
		// of asymptotic chi-square p-values
		primary_pval = chiprob_p(primary_se, cur_constraint_ct);
	      }
	    }
	    if (primary_pval > pfilter) {
	      if (cur_pval_write) {
		cur_pval_write[variant_bidx] = -9;
	      }
	      if (cur_chisq_write) {
		cur_chisq_write[variant_bidx] = -9;
	      }
	      continue;
	    }
	  }
	  const logistic_aux_result_t* auxp = &(cur_block_aux[variant_bidx]);
	  uintptr_t variant_allele_idx_base = write_variant_uidx * 2;
	  if (variant_allele_idxs) {
	    variant_allele_idx_base = variant_allele_idxs[write_variant_uidx];
	    cur_allele_ct = variant_allele_idxs[write_variant_uidx + 1] - variant_allele_idxs[write_variant_uidx];
	  }
	  char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
	  // possible todo: make number-to-string operations, strlen(), etc.
	  //   happen only once per variant.
	  for (uint32_t test_idx = 0; test_idx < cur_reported_test_ct; ++test_idx) {
	    if (chr_col) {
	      cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
	    }
	    if (variant_bps) {
	      cswritep = uint32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
	    }
	    cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
	    if (ref_col) {
	      *cswritep++ = '\t';
	      cswritep = strcpya(cswritep, cur_alleles[0]);
	    }
	    if (alt1_col) {
	      *cswritep++ = '\t';
	      cswritep = strcpya(cswritep, cur_alleles[1]);
	    }
	    if (alt_col) {
	      *cswritep++ = '\t';
	      for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
		if (cswrite(&css, &cswritep)) {
		  goto glm_logistic_ret_WRITE_FAIL;
		}
		cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
	      }
	      --cswritep;
	    }
	    if (alt_ct_col) {
	      *cswritep++ = '\t';
	      cswritep = dtoa_g(auxp->alt_dosage, cswritep);
	    }
	    if (tot_allele_col) {
	      *cswritep++ = '\t';
	      cswritep = uint32toa(auxp->allele_obs_ct, cswritep);
	    }
	    if (alt_ct_cc_col) {
	      *cswritep++ = '\t';
	      cswritep = dtoa_g(auxp->alt_case_dosage, cswritep);
	      *cswritep++ = '\t';
	      cswritep = dtoa_g(auxp->alt_dosage - auxp->alt_case_dosage, cswritep);
	    }
	    if (tot_allele_cc_col) {
	      *cswritep++ = '\t';
	      cswritep = uint32toa_x(auxp->case_allele_obs_ct, '\t', cswritep);
	      cswritep = uint32toa(auxp->allele_obs_ct - auxp->case_allele_obs_ct, cswritep);
	    }
	    if (alt_freq_col) {
	      *cswritep++ = '\t';
	      cswritep = dtoa_g(auxp->alt_dosage / ((double)auxp->allele_obs_ct), cswritep);
	    }
	    if (alt_freq_cc_col) {
	      *cswritep++ = '\t';
	      cswritep = dtoa_g(auxp->alt_case_dosage / ((double)auxp->case_allele_obs_ct), cswritep);
	      *cswritep++ = '\t';
	      cswritep = dtoa_g((auxp->alt_dosage - auxp->alt_case_dosage) / ((double)(auxp->allele_obs_ct - auxp->case_allele_obs_ct)), cswritep);
	    }
	    if (mach_r2_col) {
	      *cswritep++ = '\t';
	      if (!suppress_mach_r2) {
		cswritep = dtoa_g(auxp->mach_r2, cswritep);
	      } else {
		cswritep = strcpya(cswritep, "NA");
	      }
	    }
	    if (firth_yn_col) {
	      *cswritep++ = '\t';
	      // 'Y' - 'N' = 11
	      *cswritep++ = 'N' + 11 * auxp->firth_fallback;
	    }
	    if (test_col) {
	      *cswritep++ = '\t';
	      cswritep = strcpya(cswritep, cur_test_names[test_idx]);
	    }
	    if (nobs_col) {
	      *cswritep++ = '\t';
	      cswritep = uint32toa(auxp->sample_obs_ct, cswritep);
	    }
	    double pval = -9;
	    double tstat = 0.0;
	    if ((!cur_constraint_ct) || (test_idx != primary_reported_test_idx)) {
	      double beta = *beta_se_iter++;
	      double se = *beta_se_iter++;
	      if (!is_invalid) {
		tstat = beta / se;
		pval = chiprob_p(tstat * tstat, 1);
	      }
	      if (orbeta_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  cswritep = dtoa_g(report_beta_instead_of_odds_ratio? beta : exp(beta), cswritep);
		} else {
		  cswritep = strcpya(cswritep, "NA");
		}
	      }
	      if (se_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  cswritep = dtoa_g(se, cswritep);
		} else {
		  cswritep = strcpya(cswritep, "NA");
		}
	      }
	      if (ci_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  const double ci_halfwidth = ci_zt * se;
		  if (report_beta_instead_of_odds_ratio) {
		    cswritep = dtoa_g(beta - ci_halfwidth, cswritep);
		    *cswritep++ = '\t';
		    cswritep = dtoa_g(beta + ci_halfwidth, cswritep);
		  } else {
		    cswritep = dtoa_g(exp(beta - ci_halfwidth), cswritep);
		    *cswritep++ = '\t';
		    cswritep = dtoa_g(exp(beta + ci_halfwidth), cswritep);
		  }
		} else {
		  cswritep = strcpya(cswritep, "NA\tNA");
		}
	      }
	      if (t_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  cswritep = dtoa_g(tstat, cswritep);
		} else {
		  cswritep = strcpya(cswritep, "NA");
		}
	      }
	    } else {
	      // joint test: use (currently approximate) F-test instead of T
	      // test
	      // beta_se_iter = &(beta_se_iter[2]);
	      if (orbeta_col) {
		cswritep = memcpyl3a(cswritep, "\tNA");
	      }
	      if (se_col) {
		cswritep = memcpyl3a(cswritep, "\tNA");
	      }
	      if (ci_col) {
		cswritep = strcpya(cswritep, "\tNA\tNA");
	      }
	      if (t_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  cswritep = dtoa_g(primary_se, cswritep);
		} else {
		  cswritep = strcpya(cswritep, "NA");
		}
	      }
	      // could avoid recomputing
	      if (!is_invalid) {
		pval = chiprob_p(primary_se, cur_constraint_ct);
	      }
	    }
	    if (p_col) {
	      *cswritep++ = '\t';
	      if (!is_invalid) {
		cswritep = dtoa_g(MAXV(pval, output_min_p), cswritep);
	      } else {
		cswritep = strcpya(cswritep, "NA");
	      }
	    }
	    append_binary_eoln(&cswritep);
	    if (cswrite(&css, &cswritep)) {
	      goto glm_logistic_ret_WRITE_FAIL;
	    }
	    if (test_idx == primary_reported_test_idx) {
	      if (cur_pval_write) {
		cur_pval_write[variant_bidx] = pval;
	      }
	      if (cur_chisq_write) {
		if (!is_invalid) {
		  if (!cur_constraint_ct) {
		    cur_chisq_write[variant_bidx] = tstat * tstat;
		  } else {
		    cur_chisq_write[variant_bidx] = primary_se;
		  }
		} else {
		  cur_chisq_write[variant_bidx] = -9;
		}
	      }
	    }
	  }
	}
      }
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
      prev_block_variant_ct = cur_block_variant_ct;
      variant_idx += cur_block_variant_ct;
      // crucially, this is independent of the pgen_reader_t block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (cswrite_close_null(&css, cswritep)) {
      goto glm_logistic_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
    LOGPRINTF("Results written to %s .\n", outname);
    bigstack_reset(bigstack_mark);
  }
  while (0) {
  glm_logistic_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  glm_logistic_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  glm_logistic_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  glm_logistic_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  glm_logistic_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  glm_logistic_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 glm_logistic_ret_1:
  threads3z_cleanup(&ts, &g_cur_block_variant_ct);
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}

static const double kSmallDoubles[4] = {0.0, 1.0, 2.0, 3.0};

void genoarr_to_doubles(const uintptr_t* genoarr, uint32_t sample_ct, double* doublebuf) {
  assert(sample_ct);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t widx = 0;
  uint32_t subgroup_len = kBitsPerWordD2;
  double* doublebuf_iter = doublebuf;
  while (1) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
	return;
      }
      subgroup_len = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
      const uintptr_t cur_geno = geno_word & 3;
      // *doublebuf_iter++ = (double)((int32_t)cur_geno);
      *doublebuf_iter++ = kSmallDoubles[cur_geno];
      geno_word >>= 2;
    }
    ++widx;
  }
}

uint32_t genoarr_to_doubles_remove_missing(const uintptr_t* genoarr, uint32_t sample_ct, double* doublebuf) {
  assert(sample_ct);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t widx = 0;
  uint32_t subgroup_len = kBitsPerWordD2;
  double* doublebuf_iter = doublebuf;
  while (1) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
	return (uintptr_t)(doublebuf_iter - doublebuf);
      }
      subgroup_len = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
      const uintptr_t cur_geno = geno_word & 3;
      if (cur_geno < 3) {
	// *doublebuf_iter++ = (double)((int32_t)cur_geno);
	*doublebuf_iter++ = kSmallDoubles[cur_geno];
      }
      geno_word >>= 2;
    }
    ++widx;
  }
}

uintptr_t get_linear_workspace_size(uint32_t sample_ct, uint32_t predictor_ct, uint32_t constraint_ct, uint32_t genod_buffer_needed) {
  // sample_ct * predictor_ct < 2^31, and sample_ct >= predictor_ct, so no
  // overflows
  // could round everything up to multiples of 16 instead of 64

  // sample_nm, male_nm = sample_ctl words
  uintptr_t workspace_size = 2 * round_up_pow2(BITCT_TO_WORDCT(sample_ct) * sizeof(intptr_t), kCacheline);
  
  // nm_pheno_buf = sample_ct doubles
  workspace_size += round_up_pow2(sample_ct * sizeof(double), kCacheline);
  
  // predictors_pmaj = (predictor_ct + genod_buffer_needed) * sample_ct doubles
  workspace_size += round_up_pow2((predictor_ct + genod_buffer_needed) * sample_ct * sizeof(double), kCacheline);

  // xtx_inv, dbl_2d_buf = predictor_ct * predictor_ct doubles
  workspace_size += 2 * round_up_pow2(predictor_ct * predictor_ct * sizeof(double), kCacheline);

  // fitted_coefs, xt_y = predictor_ct doubles
  workspace_size += 2 * round_up_pow2(predictor_ct * sizeof(double), kCacheline);

#ifdef NOLAPACK
  // mi_buf = constraint_ct * kMatrixInvertBuf1CheckedAlloc bytes
  workspace_size += round_up_pow2(constraint_ct * kMatrixInvertBuf1CheckedAlloc, kCacheline);
#endif
  if (constraint_ct) {
    // tmphxs_buf, h_transpose_buf = constraint_ct * predictor_ct doubles
    workspace_size += 2 * round_up_pow2(constraint_ct * predictor_ct * sizeof(double), kCacheline);

    // inner_buf = constraint_ct * constraint_ct
    workspace_size += round_up_pow2(constraint_ct * constraint_ct * sizeof(double), kCacheline);

#ifndef NOLAPACK
    // mi_buf = constraint_ct * kMatrixInvertBuf1CheckedAlloc bytes
    workspace_size += round_up_pow2(constraint_ct * kMatrixInvertBuf1CheckedAlloc, kCacheline);
#endif
  }
  return workspace_size;  
}

THREAD_FUNC_DECL glm_linear_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec = g_genovecs[tidx];
  uintptr_t* dosage_present = nullptr;
  dosage_t* dosage_vals = nullptr;
  if (g_dosage_presents) {
    dosage_present = g_dosage_presents[tidx];
    dosage_vals = g_dosage_val_bufs[tidx];
  }
  unsigned char* workspace_buf = g_workspace_bufs[tidx];
  const uintptr_t* variant_include = g_variant_include;
  const uintptr_t* sex_male_collapsed = g_sex_male_collapsed;
  const chr_info_t* cip = g_cip;
  const uint32_t* subset_chr_fo_vidx_start = g_subset_chr_fo_vidx_start;
  // const uint32_t raw_sample_ct = g_raw_sample_ct;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const glm_flags_t glm_flags = g_glm_flags;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t joint_genotypic = (glm_flags / kfGlmGenotypic) & 1;
  const uint32_t joint_hethom = (glm_flags / kfGlmHethom) & 1;
  const uint32_t domdev_present = joint_genotypic || joint_hethom;
  const uint32_t domdev_present_p1 = domdev_present + 1;  
  const uint32_t reported_pred_uidx_start = 1 - include_intercept;
  const int32_t x_code = cip->xymt_codes[kChrOffsetX];
  const int32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = g_is_xchr_model_1;
  const uintptr_t max_reported_test_ct = g_max_reported_test_ct;
  const uintptr_t local_covar_ct = g_local_covar_ct;
  uintptr_t max_sample_ct = MAXV(g_sample_ct, g_sample_ct_x);
  if (max_sample_ct < g_sample_ct_y) {
    max_sample_ct = g_sample_ct_y;
  }
  uint32_t variant_idx_offset = 0;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_variant_ct = g_cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    double* beta_se_iter = &(g_block_beta_se[2 * max_reported_test_ct * variant_bidx]);
    linear_aux_result_t* block_aux_iter = &(g_linear_block_aux[variant_bidx]);
    const double* local_covars_iter = nullptr;
    if (local_covar_ct) {
      // &(nullptr[0]) is okay in C++, but undefined in C
      local_covars_iter = &(g_local_covars_vcmaj_d[parity][variant_bidx * max_sample_ct * local_covar_ct]);
    }
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = uint32arr_greater_than(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
      const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      uint32_t cur_variant_bidx_end = subset_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
      if (cur_variant_bidx_end > variant_bidx_end) {
	cur_variant_bidx_end = variant_bidx_end;
      }
      const uint32_t is_x = (chr_idx == x_code);
      const uint32_t is_y = (chr_idx == y_code);
      const uint32_t is_nonx_haploid = (!is_x) && is_set(cip->haploid_mask, chr_idx);
      const uintptr_t* cur_sample_include;
      const uint32_t* cur_sample_include_cumulative_popcounts;
      const double* cur_pheno;
      const double* cur_covars_cmaj;
      const uintptr_t* cur_parameter_subset;
      const double* cur_constraints_con_major;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      uint32_t cur_constraint_ct;
      uint32_t primary_pred_idx = include_intercept;
      if (is_y && g_sample_include_y) {
	cur_sample_include = g_sample_include_y;
	cur_sample_include_cumulative_popcounts = g_sample_include_y_cumulative_popcounts;
	cur_pheno = g_pheno_y_d;
	cur_covars_cmaj = g_covars_cmaj_y_d;
	cur_parameter_subset = g_parameter_subset_y;
	cur_constraints_con_major = g_constraints_con_major_y;
	cur_sample_ct = g_sample_ct_y;
	cur_covar_ct = g_covar_ct_y;
	cur_constraint_ct = g_constraint_ct_y;
      } else if (is_x && g_sample_include_x) {
	cur_sample_include = g_sample_include_x;
	cur_sample_include_cumulative_popcounts = g_sample_include_x_cumulative_popcounts;
	cur_pheno = g_pheno_x_d;
	cur_covars_cmaj = g_covars_cmaj_x_d;
	cur_parameter_subset = g_parameter_subset_x;
	cur_constraints_con_major = g_constraints_con_major_x;
	cur_sample_ct = g_sample_ct_x;
	cur_covar_ct = g_covar_ct_x;
	cur_constraint_ct = g_constraint_ct_x;
      } else {
	cur_sample_include = g_sample_include;
	cur_sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
	cur_pheno = g_pheno_d;
	cur_covars_cmaj = g_covars_cmaj_d;
	cur_parameter_subset = g_parameter_subset;
	cur_constraints_con_major = g_constraints_con_major;
	cur_sample_ct = g_sample_ct;
	cur_covar_ct = g_covar_ct;
	cur_constraint_ct = g_constraint_ct;
      }
      const uint32_t sample_ctl = BITCT_TO_WORDCT(cur_sample_ct);
      const uint32_t cur_predictor_ct_base = 2 + domdev_present + cur_covar_ct * (1 + add_interactions * domdev_present_p1);
      uint32_t cur_predictor_ct = cur_predictor_ct_base;
      if (cur_parameter_subset) {
	cur_predictor_ct = popcount_longs(cur_parameter_subset, BITCT_TO_WORDCT(cur_predictor_ct_base));
      }
      uint32_t reported_pred_uidx_end;
      if (hide_covar) {
	if (!cur_parameter_subset) {
	  reported_pred_uidx_end = 2 + domdev_present;
	} else {
	  reported_pred_uidx_end = 1 + is_set(cur_parameter_subset, 1) + domdev_present;
	}
      } else {
	reported_pred_uidx_end = cur_predictor_ct;
      }
      // todo: --tests
      if (cur_constraint_ct) {
	primary_pred_idx = reported_pred_uidx_end - reported_pred_uidx_start;
      }
      const uint32_t genod_buffer_needed = cur_parameter_subset && (!is_set(cur_parameter_subset, 1));
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = (uintptr_t*)arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter);
      uintptr_t* male_nm = (uintptr_t*)arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter);
      double* nm_pheno_buf = (double*)arena_alloc_raw_rd(cur_sample_ct * sizeof(double), &workspace_iter);
      double* nm_predictors_pmaj_buf = (double*)arena_alloc_raw_rd((cur_predictor_ct + genod_buffer_needed) * cur_sample_ct * sizeof(double), &workspace_iter);
      double* xtx_inv = (double*)arena_alloc_raw_rd(cur_predictor_ct * cur_predictor_ct * sizeof(double), &workspace_iter);
      double* fitted_coefs = (double*)arena_alloc_raw_rd(cur_predictor_ct * sizeof(double), &workspace_iter);
      double* xt_y = (double*)arena_alloc_raw_rd(cur_predictor_ct * sizeof(double), &workspace_iter);
      double* dbl_2d_buf = (double*)arena_alloc_raw_rd(cur_predictor_ct * cur_predictor_ct * sizeof(double), &workspace_iter);
      
      // joint test only
      matrix_invert_buf1_t* inv_1d_buf = nullptr;
      double* tmphxs_buf = nullptr;
      double* h_transpose_buf = nullptr;
      double* inner_buf = nullptr;
#ifdef NOLAPACK
      // (well, except if LAPACK is missing)
      inv_1d_buf = (matrix_invert_buf1_t*)arena_alloc_raw_rd(cur_predictor_ct * kMatrixInvertBuf1CheckedAlloc, &workspace_iter);
#endif
      if (cur_constraint_ct) {
#ifndef NOLAPACK
	inv_1d_buf = (matrix_invert_buf1_t*)arena_alloc_raw_rd(cur_predictor_ct * kMatrixInvertBuf1CheckedAlloc, &workspace_iter);
#endif
	tmphxs_buf = (double*)arena_alloc_raw_rd(cur_constraint_ct * cur_predictor_ct * sizeof(double), &workspace_iter);
	h_transpose_buf = (double*)arena_alloc_raw_rd(cur_constraint_ct * cur_predictor_ct * sizeof(double), &workspace_iter);
	inner_buf = (double*)arena_alloc_raw_rd(cur_constraint_ct * cur_constraint_ct * sizeof(double), &workspace_iter);
      }
      assert((uintptr_t)(workspace_iter - workspace_buf) == get_linear_workspace_size(cur_sample_ct, cur_predictor_ct, cur_constraint_ct, genod_buffer_needed));
      double pheno_ssq_base = 0.0;
      for (uint32_t sample_idx = 0; sample_idx < cur_sample_ct; ++sample_idx) {
	pheno_ssq_base += cur_pheno[sample_idx] * cur_pheno[sample_idx];
      }
      pgr_clear_ld_cache(pgrp);
      uint32_t genocounts[4];
      for (; variant_bidx < cur_variant_bidx_end; ++variant_bidx, ++variant_uidx) {
	next_set_unsafe_ck(variant_include, &variant_uidx);
	{
	  uint32_t dosage_ct;
	  uint32_t is_explicit_alt1;
	  pglerr_t reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(cur_sample_include, cur_sample_include_cumulative_popcounts, cur_sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
	  if (reterr) {
	    g_error_ret = reterr;
	    variant_bidx = variant_bidx_end;
	    break;
	  }
	  zero_trailing_quaters(cur_sample_ct, genovec);
	  genovec_count_freqs_unsafe(genovec, cur_sample_ct, genocounts);
	  uint32_t missing_ct = genocounts[3];
	  if (!missing_ct) {
	    fill_all_bits(cur_sample_ct, sample_nm);
	  } else {
	    genoarr_to_nonmissing(genovec, cur_sample_ct, sample_nm);
	    if (dosage_ct) {
	      bitvec_or(dosage_present, sample_ctl, sample_nm);
	      missing_ct = cur_sample_ct - popcount_longs(sample_nm, sample_ctl);
	    }
	  }
	  uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
	  // todo: alt2/alt3/etc. dosage > 0.5 -> missing
	  const uint32_t nm_sample_ctl = BITCT_TO_WORDCT(nm_sample_ct);
	  double* nm_predictors_pmaj_iter = nm_predictors_pmaj_buf;
	  // first predictor column: intercept
	  for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	    *nm_predictors_pmaj_iter++ = 1.0;
	  }
	  // second predictor column: genotype
	  double* genotype_vals = &(nm_predictors_pmaj_buf[nm_sample_ct]);
	  if (genod_buffer_needed) {
	    // special case: --parameters excludes the main genotype column,
	    // but does care about an interaction
	    genotype_vals = &(nm_predictors_pmaj_buf[cur_predictor_ct * nm_sample_ct]);
	  }
	  nm_predictors_pmaj_iter = genotype_vals;
	  double cur_pheno_ssq = pheno_ssq_base;
	  if (!missing_ct) {
	    genoarr_to_doubles(genovec, nm_sample_ct, nm_predictors_pmaj_iter);
	    if (dosage_ct) {
	      uint32_t sample_idx = 0;
	      for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_idx) {
		next_set_unsafe_ck(dosage_present, &sample_idx);
		// 32768 -> 2, 16384 -> 1, 0 -> 0
		nm_predictors_pmaj_iter[sample_idx] = kRecipDosageMid * ((int32_t)((uint32_t)dosage_vals[dosage_idx]));
	      }
	    }
	  } else {
	    uint32_t sample_midx = 0;
	    for (uint32_t missing_idx = 0; missing_idx < missing_ct; ++missing_idx, ++sample_midx) {
	      next_unset_unsafe_ck(sample_nm, &sample_midx);
	      cur_pheno_ssq -= cur_pheno[sample_midx] * cur_pheno[sample_midx];
	    }
	    if (!dosage_ct) {
	      genoarr_to_doubles_remove_missing(genovec, cur_sample_ct, nm_predictors_pmaj_iter);
	    } else {
	      sample_midx = 0;
	      uint32_t dosage_idx = 0;
	      for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
		next_set_unsafe_ck(sample_nm, &sample_midx);
		double cur_val;
		if (is_set(dosage_present, sample_midx)) {
		  cur_val = kRecipDosageMid * ((int32_t)((uint32_t)dosage_vals[dosage_idx++]));
		} else {
		  cur_val = (intptr_t)(GET_QUATERARR_ENTRY(genovec, sample_midx));
		}
	        nm_predictors_pmaj_iter[sample_idx] = cur_val;
	      }
	    }
	  }
	  nm_predictors_pmaj_iter = &(nm_predictors_pmaj_iter[nm_sample_ct]);
	  // usually need to save some of {sample_obs_ct, allele_obs_ct,
	  // alt_dosage, mach_r2 even for skipped variants
	  // compute them all for now, could conditionally skip later
	  block_aux_iter->sample_obs_ct = nm_sample_ct;
	  double dosage_ceil = 2.0;
	  if (!is_x) {
	    if (!is_nonx_haploid) {
	      block_aux_iter->allele_obs_ct = nm_sample_ct * 2;
	    } else {
	      block_aux_iter->allele_obs_ct = nm_sample_ct;
	      // everything is on 0..1 scale, not 0..2
	      dosage_ceil = 1.0;
	      for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
		genotype_vals[sample_idx] *= 0.5;
	      }
	    }
	  } else {
	    copy_bitarr_subset(sex_male_collapsed, sample_nm, nm_sample_ct, male_nm);
	    const uint32_t nm_male_ct = popcount_longs(male_nm, nm_sample_ctl);
	    block_aux_iter->allele_obs_ct = nm_sample_ct * 2;
	    if (is_xchr_model_1) {
	      // special case: multiply male values by 0.5
	      uint32_t sample_idx = 0;
	      for (uint32_t male_idx = 0; male_idx < nm_male_ct; ++male_idx, ++sample_idx) {
		next_set_unsafe_ck(male_nm, &sample_idx);
		genotype_vals[sample_idx] *= 0.5;
	      }
	      block_aux_iter->allele_obs_ct -= nm_male_ct;
	    }
	  }
	  double dosage_sum = 0.0;
	  double dosage_ssq = 0.0;
	  for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	    const double cur_genotype_val = genotype_vals[sample_idx];
	    dosage_sum += cur_genotype_val;
	    dosage_ssq += cur_genotype_val * cur_genotype_val;
	  }
	  block_aux_iter->alt_dosage = dosage_sum;

	  const double dosage_avg = dosage_sum / ((double)((int32_t)nm_sample_ct));
	  const double dosage_variance = dosage_ssq - dosage_sum * dosage_avg;
	  block_aux_iter->mach_r2 = 2 * dosage_variance / (dosage_sum * (dosage_ceil - dosage_avg));
	  // okay, now we're free to skip the actual regression if there are
	  // too few samples, or variant is monomorphic (or all-het)
	  if ((nm_sample_ct <= cur_predictor_ct) || (fabs(dosage_variance) < kBigEpsilon)) {
	    goto glm_linear_thread_skip_variant;
	  }
	  double* domdev_vals = nullptr;
	  if (genod_buffer_needed) {
	    nm_predictors_pmaj_iter = &(nm_predictors_pmaj_buf[nm_sample_ct]);
	  } else if (joint_genotypic || joint_hethom) {
	    // in hethom case, do this before clobbering genotype data
	    domdev_vals = nm_predictors_pmaj_iter;
	    for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	      double cur_genotype_val = genotype_vals[sample_idx];
	      if (cur_genotype_val > 1.0) {
		cur_genotype_val = 2.0 - cur_genotype_val;
	      }
	      nm_predictors_pmaj_iter[sample_idx] = cur_genotype_val;
	    }
	    nm_predictors_pmaj_iter = &(nm_predictors_pmaj_iter[nm_sample_ct]);
	  }
	  if (model_dominant) {
	    for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	      const double cur_genotype_val = genotype_vals[sample_idx];
	      // 0..1..1
	      if (cur_genotype_val > 1.0) {
		genotype_vals[sample_idx] = 1.0;
	      }
	    }
	  } else if (model_recessive || joint_hethom) {
	    for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx) {
	      const double cur_genotype_val = genotype_vals[sample_idx];
	      // 0..0..1
	      if (cur_genotype_val < 1.0) {
		genotype_vals[sample_idx] = 0.0;
	      } else {
		genotype_vals[sample_idx] = cur_genotype_val - 1.0;
	      }
	    }
	  }

	  // fill phenotype
	  uint32_t sample_midx = 0;
	  for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
	    next_set_unsafe_ck(sample_nm, &sample_midx);
	    nm_pheno_buf[sample_idx] = cur_pheno[sample_midx];
	  }
	  
	  // fill covariates
	  uint32_t parameter_uidx = 2 + domdev_present;
	  for (uint32_t covar_idx = 0; covar_idx < cur_covar_ct; ++covar_idx, ++parameter_uidx) {
	    // strictly speaking, we don't need cur_covars_cmaj to be
	    // vector-aligned
	    if (cur_parameter_subset && (!is_set(cur_parameter_subset, parameter_uidx))) {
	      continue;
	    }
	    const double* cur_covar_col;
	    if (covar_idx < local_covar_ct) {
	      cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
	    } else {
	      cur_covar_col = &(cur_covars_cmaj[(covar_idx - local_covar_ct) * cur_sample_ct]);
	    }
	    sample_midx = 0;
	    for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
	      next_set_unsafe_ck(sample_nm, &sample_midx);
	      *nm_predictors_pmaj_iter++ = cur_covar_col[sample_midx];
	    }
	  }
	  // fill interaction terms
	  if (add_interactions) {
	    for (uint32_t covar_idx = 0; covar_idx < cur_covar_ct; ++covar_idx) {
	      const double* cur_covar_col;
	      if (covar_idx < local_covar_ct) {
	        cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
	      } else {
		cur_covar_col = &(cur_covars_cmaj[covar_idx * cur_sample_ct]);
	      }
	      if ((!cur_parameter_subset) || is_set(cur_parameter_subset, parameter_uidx)) {
		sample_midx = 0;
		for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
		  next_set_unsafe_ck(sample_nm, &sample_midx);
		  *nm_predictors_pmaj_iter++ = genotype_vals[sample_idx] * cur_covar_col[sample_midx];
		}
	      }
	      ++parameter_uidx;
	      if (domdev_present) {
		if ((!cur_parameter_subset) || is_set(cur_parameter_subset, parameter_uidx)) {
		  sample_midx = 0;
		  for (uint32_t sample_idx = 0; sample_idx < nm_sample_ct; ++sample_idx, ++sample_midx) {
		    next_set_unsafe_ck(sample_nm, &sample_midx);
		    *nm_predictors_pmaj_iter++ = domdev_vals[sample_idx] * cur_covar_col[sample_midx];
		  }
		}
		++parameter_uidx;
	      }
	    }
	  }
	  if (linear_regression_inv(nm_pheno_buf, nm_predictors_pmaj_buf, cur_predictor_ct, nm_sample_ct, fitted_coefs, xtx_inv, xt_y, inv_1d_buf, dbl_2d_buf)) {
	    goto glm_linear_thread_skip_variant;
	  }
	  // RSS = y^T y - y^T X (X^T X)^{-1} X^T y
	  //     = cur_pheno_ssq - xt_y * fitted_coefs
	  // s^2 = RSS / df
	  // possible todo: improve numerical stability of this computation in
	  // non-mean-centered phenotype case
	  double sigma = cur_pheno_ssq;
	  for (uint32_t pred_uidx = 0; pred_uidx < cur_predictor_ct; ++pred_uidx) {
	    sigma -= xt_y[pred_uidx] * fitted_coefs[pred_uidx];
	  }
	  sigma /= nm_sample_ct - cur_predictor_ct;
	  for (uint32_t uii = 0; uii < cur_predictor_ct; ++uii) {
	    double* s_iter = &(xtx_inv[uii * cur_predictor_ct]);
#ifdef NOLAPACK
	    for (uint32_t ujj = 0; ujj < cur_predictor_ct; ++ujj) {
	      *s_iter *= sigma;
	      ++s_iter;
	    }
#else
	    for (uint32_t ujj = 0; ujj <= uii; ++ujj) {
	      *s_iter *= sigma;
	      ++s_iter;
	    }
#endif
	  }
	  // validParameters() check
	  for (uint32_t pred_uidx = 1; pred_uidx < cur_predictor_ct; ++pred_uidx) {
	    const double xtx_inv_diag_element = xtx_inv[pred_uidx * (cur_predictor_ct + 1)];
	    if (xtx_inv_diag_element < 1e-20) {
	      goto glm_linear_thread_skip_variant;
	    }
	    // use dbl_2d_buf[] to store diagonal square roots
	    dbl_2d_buf[pred_uidx] = sqrt(xtx_inv_diag_element);
	  }
	  dbl_2d_buf[0] = sqrt(xtx_inv[0]);
	  for (uint32_t pred_uidx = 1; pred_uidx < cur_predictor_ct; ++pred_uidx) {
	    const double cur_xtx_inv_diag_sqrt = 0.99999 * dbl_2d_buf[pred_uidx];
	    const double* xtx_inv_row = &(xtx_inv[pred_uidx * cur_predictor_ct]);
	    for (uint32_t pred_uidx2 = 0; pred_uidx2 < pred_uidx; ++pred_uidx2) {
	      if (xtx_inv_row[pred_uidx2] > cur_xtx_inv_diag_sqrt * dbl_2d_buf[pred_uidx2]) {
		goto glm_linear_thread_skip_variant;
	      }
	    }
	  }
	  double* beta_se_iter2 = beta_se_iter;
	  for (uint32_t pred_uidx = reported_pred_uidx_start; pred_uidx < reported_pred_uidx_end; ++pred_uidx) {
	    *beta_se_iter2++ = fitted_coefs[pred_uidx];
	    *beta_se_iter2++ = dbl_2d_buf[pred_uidx];
	  }
	  if (cur_constraint_ct) {
	    *beta_se_iter2++ = 0.0;
#ifndef NOLAPACK
	    // xtx_inv upper triangle was not filled
	    for (uint32_t row_idx = 0; row_idx < cur_predictor_ct; ++row_idx) {
	      double* cur_row = &(xtx_inv[row_idx * cur_predictor_ct]);
	      double* cur_col = &(xtx_inv[row_idx]);
	      for (uint32_t col_idx = row_idx + 1; col_idx < cur_predictor_ct; ++col_idx) {
		cur_row[col_idx] = cur_col[col_idx * cur_predictor_ct];
	      }
	    }
#endif
	    double chisq;
	    if (!linear_hypothesis_chisq(fitted_coefs, cur_constraints_con_major, xtx_inv, cur_constraint_ct, cur_predictor_ct, &chisq, tmphxs_buf, h_transpose_buf, inner_buf, inv_1d_buf, dbl_2d_buf)) {
	      *beta_se_iter2++ = chisq;
	    } else {
	      *beta_se_iter2++ = -9;
	    }
	  }
	}
	while (0) {
	glm_linear_thread_skip_variant:
	  beta_se_iter[primary_pred_idx * 2 + 1] = -9;
	}
	beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
	++block_aux_iter;
	if (local_covars_iter) {
	  local_covars_iter = &(local_covars_iter[local_covar_ct * max_sample_ct]);
	}
	// todo?
      }
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
    variant_idx_offset += cur_block_variant_ct;
  }
}

pglerr_t glm_linear(const char* cur_pheno_name, char** test_names, char** test_names_x, char** test_names_y, const uint32_t* variant_bps, char** variant_ids, char** allele_storage, const glm_info_t* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const char* outname, uint32_t raw_variant_ct, uint32_t max_chr_blen, double ci_size, double pfilter, double output_min_p, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uint32_t local_sample_ct, uint32_t local_loadbuf_size, pgen_file_info_t* pgfip, gzFile gz_local_covar_file, uintptr_t* valid_variants, double* orig_pvals, double* orig_chisq, unsigned char* overflow_buf, char* local_loadbuf) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  threads_state_t ts;
  init_threads3z(&ts);
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uintptr_t* variant_include = g_variant_include;
    const chr_info_t* cip = g_cip;
    const uintptr_t* variant_allele_idxs = g_variant_allele_idxs;

    const uint32_t sample_ct = g_sample_ct;
    const uint32_t sample_ct_x = g_sample_ct_x;
    const uint32_t sample_ct_y = g_sample_ct_y;
    const uint32_t covar_ct = g_covar_ct;
    const uintptr_t local_covar_ct = g_local_covar_ct;
    const uint32_t covar_ct_x = g_covar_ct_x;
    const uint32_t covar_ct_y = g_covar_ct_y;

    uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
    if (max_sample_ct < sample_ct_y) {
      max_sample_ct = sample_ct_y;
    }
    uint32_t* local_sample_idx_order = nullptr;
    uint32_t local_line_idx = 0;
    uint32_t local_xy = 0; // 1 = chrX, 2 = chrY
    if (gz_local_covar_file) {
      if (gzrewind(gz_local_covar_file)) {
	goto glm_linear_ret_READ_FAIL;
      }
      if (bigstack_alloc_ui(local_sample_ct, &local_sample_idx_order)) {
	goto glm_linear_ret_NOMEM;
      }
      for (uint32_t uii = 0; uii < local_sample_ct; ++uii) {
	const uint32_t cur_uidx = local_sample_uidx_order[uii];
	uint32_t cur_idx = 0xffffffffU;
	if ((cur_uidx != 0xffffffffU) && is_set(g_sample_include, cur_uidx)) {
	  cur_idx = raw_to_subsetted_pos(g_sample_include, g_sample_include_cumulative_popcounts, cur_uidx);
	}
	local_sample_idx_order[uii] = cur_idx;
      }
    }

    const uint32_t variant_ct = g_variant_ct;

    const glm_flags_t glm_flags = glm_info_ptr->flags;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
      goto glm_linear_ret_OPEN_FAIL;
    }
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t domdev_present_p1 = domdev_present + 1;

    // todo: --tests
    const uint32_t constraint_ct = g_constraint_ct;
    const uint32_t constraint_ct_x = g_constraint_ct_x;
    const uint32_t constraint_ct_y = g_constraint_ct_y;
    
    uint32_t predictor_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
    uint32_t predictor_ct_x = 2 + domdev_present + covar_ct_x * (1 + add_interactions * domdev_present_p1);
    uint32_t predictor_ct_y = 2 + domdev_present + covar_ct_y * (1 + add_interactions * domdev_present_p1);
    const uintptr_t* parameter_subset = g_parameter_subset;
    const uintptr_t* parameter_subset_x = g_parameter_subset_x;
    const uintptr_t* parameter_subset_y = g_parameter_subset_y;
    if (parameter_subset) {
      predictor_ct = popcount_longs(parameter_subset, BITCT_TO_WORDCT(predictor_ct));
      if (sample_ct_x) {
	predictor_ct_x = popcount_longs(parameter_subset_x, BITCT_TO_WORDCT(predictor_ct_x));
      } else {
	predictor_ct_x = 0;
      }
      if (sample_ct_y) {
	predictor_ct_y = popcount_longs(parameter_subset_y, BITCT_TO_WORDCT(predictor_ct_x));
      } else {
	predictor_ct_y = 0;
      }
    }
    uint32_t reported_test_ct = get_reported_test_ct(parameter_subset, glm_flags, covar_ct);
    uintptr_t max_reported_test_ct = reported_test_ct;
    uint32_t reported_test_ct_x = 0;
    if (sample_ct_x) {
      reported_test_ct_x = get_reported_test_ct(parameter_subset_x, glm_flags, covar_ct_x);
      if (reported_test_ct_x > max_reported_test_ct) {
	max_reported_test_ct = reported_test_ct_x;
      }
    }
    uint32_t reported_test_ct_y = 0;
    if (sample_ct_y) {
      reported_test_ct_y = get_reported_test_ct(parameter_subset_y, glm_flags, covar_ct_y);
      if (reported_test_ct_y > max_reported_test_ct) {
	max_reported_test_ct = reported_test_ct_y;
      }
    }
    const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
    const glm_cols_t glm_cols = glm_info_ptr->cols;
    const uint32_t test_col = glm_cols & kfGlmColTest;
    if ((!test_col) && (max_reported_test_ct > 1)) {
      logerrprint("Error: --glm's 'test' column cannot be omitted when results for multiple\npredictors are reported.  (Did you forget 'hide-covar'?)\n");
      goto glm_linear_ret_INCONSISTENT_INPUT;
    }
    g_max_reported_test_ct = max_reported_test_ct;

    int32_t x_code = -2;
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    if (sample_ct_x) {
      get_xymt_code_start_and_end_unsafe(cip, kChrOffsetX, &x_code, &x_start, &x_end);
    }
    int32_t y_code = -2;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    if (sample_ct_y) {
      get_xymt_code_start_and_end_unsafe(cip, kChrOffsetY, &y_code, &y_start, &y_end);
    }
    const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t chr_col = glm_cols & kfGlmColChrom;

    // includes trailing tab
    char* chr_buf = nullptr;
    if (chr_col) {
      if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
	goto glm_linear_ret_NOMEM;
      }
    }

    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > variant_ct) {
      calc_thread_ct = variant_ct;
    }

    const uint32_t genod_buffer_needed = parameter_subset && (!is_set(parameter_subset, 1));
    uintptr_t workspace_alloc = get_linear_workspace_size(sample_ct, predictor_ct, constraint_ct, genod_buffer_needed);
    if (sample_ct_x) {
      const uintptr_t workspace_alloc_x = get_linear_workspace_size(sample_ct_x, predictor_ct_x, constraint_ct_x, genod_buffer_needed);
      if (workspace_alloc_x > workspace_alloc) {
	workspace_alloc = workspace_alloc_x;
      }
    }
    if (sample_ct_y) {
      const uintptr_t workspace_alloc_y = get_linear_workspace_size(sample_ct_y, predictor_ct_y, constraint_ct_y, genod_buffer_needed);
      if (workspace_alloc_y > workspace_alloc) {
	workspace_alloc = workspace_alloc_y;
      }
    }
    // +1 is for top-level g_workspace_bufs
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;
    uintptr_t per_variant_xalloc_byte_ct = sizeof(linear_aux_result_t) + 2 * max_reported_test_ct * sizeof(double) + max_sample_ct * local_covar_ct * sizeof(double);
    unsigned char* main_loadbufs[2];
    uint32_t read_block_size;
    if (multithread_load_init(variant_include, max_sample_ct, variant_ct, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, pgfip, &calc_thread_ct, &g_genovecs, dosage_is_present? (&g_dosage_presents) : nullptr, dosage_is_present? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &ts.threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto glm_linear_ret_NOMEM;
    }
    ts.calc_thread_ct = calc_thread_ct;
    g_calc_thread_ct = calc_thread_ct;
    linear_aux_result_t* linear_block_aux_bufs[2];
    double* block_beta_se_bufs[2];
    
    for (uint32_t uii = 0; uii < 2; ++uii) {
      linear_block_aux_bufs[uii] = (linear_aux_result_t*)bigstack_alloc(read_block_size * sizeof(linear_aux_result_t));
      if ((!linear_block_aux_bufs[uii]) ||
	  bigstack_alloc_d(read_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii]))) {
	goto glm_linear_ret_NOMEM;
      }
      if (local_covar_ct) {
	if (bigstack_alloc_d(read_block_size * max_sample_ct * local_covar_ct * sizeof(double), &(g_local_covars_vcmaj_d[uii]))) {
	  goto glm_linear_ret_NOMEM;
	}
      } else {
	g_local_covars_vcmaj_d[uii] = nullptr;
      }
    }

    g_workspace_bufs = (unsigned char**)bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t));
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_workspace_bufs[tidx] = bigstack_alloc_raw(workspace_alloc);
    }
    
    const uint32_t ref_col = glm_cols & kfGlmColRef;
    const uint32_t alt1_col = glm_cols & kfGlmColAlt1;
    const uint32_t alt_col = glm_cols & kfGlmColAlt;
    const uint32_t alt_ct_col = glm_cols & kfGlmColAltcount;
    const uint32_t tot_allele_col = glm_cols & kfGlmColTotallele;
    const uint32_t alt_freq_col = glm_cols & kfGlmColAltfreq;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t beta_col = glm_cols & (kfGlmColBeta | kfGlmColOrbeta);
    const uint32_t se_col = glm_cols & kfGlmColSe;
    const uint32_t ci_col = (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t t_col = glm_cols & kfGlmColT;
    const uint32_t p_col = glm_cols & kfGlmColP;
    cswritep = (char*)overflow_buf;
    *cswritep++ = '#';
    if (chr_col) {
      cswritep = strcpya(cswritep, "CHROM\t");
    }
    if (variant_bps) {
      cswritep = strcpya(cswritep, "POS\t");
    }
    cswritep = strcpya(cswritep, "ID");
    if (ref_col) {
      cswritep = strcpya(cswritep, "\tREF");
    }
    if (alt1_col) {
      cswritep = strcpya(cswritep, "\tALT1");
    }
    if (alt_col) {
      cswritep = strcpya(cswritep, "\tALT");
    }
    if (alt_ct_col) {
      cswritep = strcpya(cswritep, "\tALT1_CT");
    }
    if (tot_allele_col) {
      cswritep = strcpya(cswritep, "\tALLELE_CT");
    }
    if (alt_freq_col) {
      cswritep = strcpya(cswritep, "\tALT_FREQ");
    }
    if (mach_r2_col) {
      cswritep = strcpya(cswritep, "\tMACH_R2");
    }
    if (test_col) {
      cswritep = strcpya(cswritep, "\tTEST");
    }
    if (nobs_col) {
      cswritep = strcpya(cswritep, "\tOBS_CT");
    }
    if (beta_col) {
      cswritep = strcpya(cswritep, "\tBETA");
    }
    if (se_col) {
      cswritep = strcpya(cswritep, "\tSE");
    }
    double ci_zt = 0.0;
    if (ci_col) {
      cswritep = strcpya(cswritep, "\tL");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      cswritep = strcpya(cswritep, "\tU");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      ci_zt = ltqnorm((ci_size + 1.0) * 0.5);
    }
    if (t_col) {
      if (!constraint_ct) {
        cswritep = strcpya(cswritep, "\tT_STAT");
      } else {
	// chisq for joint tests.  may switch to F-statistic (just divide by
	// df; the hard part there is porting a function to convert that to a
	// p-value)
        cswritep = strcpya(cswritep, "\tT_OR_CHISQ_STAT");
      }
    }
    if (p_col) {
      cswritep = strcpya(cswritep, "\tP");
    }
    append_binary_eoln(&cswritep);

    // Main workflow:
    // 1. Set n=0, load/skip block 0
    //
    // 2. Spawn threads processing block n
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Load/skip block n unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8, Write results for last block
    const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
    const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t write_variant_uidx = 0;
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t suppress_mach_r2 = 0;

    // todo: --tests
    uint32_t cur_reported_test_ct = 0;
    uint32_t primary_reported_test_idx = include_intercept;
    uint32_t cur_predictor_ct = 0;
    uint32_t cur_constraint_ct = 0;

    char** cur_test_names = nullptr;
    uint32_t prev_block_variant_ct = 0;
    uint32_t variant_idx = 0;
    uint32_t cur_read_block_size = read_block_size;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t cur_allele_ct = 2;
    LOGPRINTFWW5("--glm linear regression on phenotype '%s': ", cur_pheno_name);
    fputs("0%", stdout);
    fflush(stdout);
    while (1) {
      uintptr_t cur_block_variant_ct = 0;
      if (!ts.is_last_block) {
	while (read_block_idx < read_block_ct_m1) {
	  cur_block_variant_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
	  if (cur_block_variant_ct) {
	    break;
	  }
	  ++read_block_idx;
	}
	if (read_block_idx == read_block_ct_m1) {
	  cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
	  cur_block_variant_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
	}
	if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_variant_ct, pgfip)) {
	  goto glm_linear_ret_READ_FAIL;
	}
	if (gz_local_covar_file) {
	  reterr = read_local_covar_block(g_sample_include, g_sample_include_x, g_sample_include_y, g_sample_include_cumulative_popcounts, g_sample_include_x_cumulative_popcounts, g_sample_include_y_cumulative_popcounts, cip, variant_include, local_sample_uidx_order, local_variant_include, sample_ct, sample_ct_x, sample_ct_y, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_variant_ct, local_sample_ct, local_covar_ct, (glm_info_ptr->flags / kfGlmLocalOmitLast) & 1, glm_info_ptr->local_cat_ct, local_loadbuf_size, gz_local_covar_file, &local_line_idx, &local_xy, nullptr, g_local_covars_vcmaj_d[parity], local_sample_idx_order, local_loadbuf);
	  if (reterr) {
	    goto glm_linear_ret_1;
	  }
	}
      }
      if (variant_idx) {
	join_threads3z(&ts);
	reterr = g_error_ret;
	if (reterr) {
	  if (reterr == kPglRetMalformedInput) {
	    logprint("\n");
	    logerrprint("Error: Malformed .pgen file.\n");
	  }
	  goto glm_linear_ret_1;
	}
      }
      if (!ts.is_last_block) {
	g_cur_block_variant_ct = cur_block_variant_ct;
	const uint32_t uidx_start = read_block_idx * read_block_size;
	compute_uidx_start_partition(variant_include, cur_block_variant_ct, calc_thread_ct, uidx_start, g_read_variant_uidx_starts);
	for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	  g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
	  g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
	}
	g_linear_block_aux = linear_block_aux_bufs[parity];
	g_block_beta_se = block_beta_se_bufs[parity];
	ts.is_last_block = (variant_idx + cur_block_variant_ct == variant_ct);
	ts.thread_func_ptr = glm_linear_thread;
	if (spawn_threads3z(variant_idx, &ts)) {
	  goto glm_linear_ret_THREAD_CREATE_FAIL;
	}
      }
      parity = 1 - parity;
      if (variant_idx) {
	// write *previous* block results
	const double* cur_block_beta_se = block_beta_se_bufs[parity];
	const linear_aux_result_t* cur_block_aux = linear_block_aux_bufs[parity];
	const uint32_t variant_idx_start = variant_idx - prev_block_variant_ct;
	double* cur_pval_write = orig_pvals? (&(orig_pvals[variant_idx_start])) : nullptr;
	double* cur_chisq_write = orig_chisq? (&(orig_chisq[variant_idx_start])) : nullptr;
	for (uint32_t variant_bidx = 0; variant_bidx < prev_block_variant_ct; ++variant_bidx, ++write_variant_uidx) {
	  next_set_unsafe_ck(variant_include, &write_variant_uidx);
	  if (write_variant_uidx >= chr_end) {
	    do {
	      ++chr_fo_idx;
	      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	    } while (write_variant_uidx >= chr_end);
	    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	    suppress_mach_r2 = 1;
	    if ((chr_idx == ((uint32_t)x_code)) && sample_ct_x) {
	      cur_reported_test_ct = reported_test_ct_x;
	      cur_predictor_ct = predictor_ct_x;
	      cur_constraint_ct = constraint_ct_x;
	      cur_test_names = test_names_x;
	    } else if ((chr_idx == ((uint32_t)y_code)) && sample_ct_y) {
	      cur_reported_test_ct = reported_test_ct_y;
	      cur_predictor_ct = predictor_ct_y;
	      cur_constraint_ct = constraint_ct_y;
	      cur_test_names = test_names_y;
	    } else {
	      cur_reported_test_ct = reported_test_ct;
	      cur_predictor_ct = predictor_ct;
	      cur_constraint_ct = constraint_ct;
	      cur_test_names = test_names;
	      if ((chr_idx != ((uint32_t)x_code)) && (chr_idx != ((uint32_t)mt_code)) && (!is_set(cip->haploid_mask, chr_idx))) {
		suppress_mach_r2 = 0;
	      }
	    }
	    if (cur_constraint_ct) {
	      primary_reported_test_idx = reported_test_ct - 1;
	    }
	    if (chr_col) {
	      char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	      *chr_name_end = '\t';
	      chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
	    }
	  }
	  const double* beta_se_iter = &(cur_block_beta_se[2 * max_reported_test_ct * variant_bidx]);
	  const double primary_beta = beta_se_iter[primary_reported_test_idx * 2];
	  const double primary_se = beta_se_iter[primary_reported_test_idx * 2 + 1];
	  const uint32_t is_invalid = (primary_se == -9);
	  if (is_invalid && valid_variants) {
	    CLEAR_BIT(write_variant_uidx, valid_variants);
	  }
	  const linear_aux_result_t* auxp = &(cur_block_aux[variant_bidx]);
	  if (pfilter != 2.0) {
	    double primary_pval = 2.0;
	    if (!is_invalid) {
	      if (!cur_constraint_ct) {
		double primary_tstat = primary_beta / primary_se;
		primary_pval = calc_tprob(primary_tstat, auxp->sample_obs_ct - cur_predictor_ct);
	      } else {
		// possible todo: support for F-distribution p-values instead
		// of asymptotic chi-square p-values
		primary_pval = chiprob_p(primary_se, cur_constraint_ct);
	      }
	    }
	    if (primary_pval > pfilter) {
	      if (cur_pval_write) {
		cur_pval_write[variant_bidx] = -9;
	      }
	      if (cur_chisq_write) {
		cur_chisq_write[variant_bidx] = -9;
	      }
	      continue;
	    }
	  }
	  uintptr_t variant_allele_idx_base = write_variant_uidx * 2;
	  if (variant_allele_idxs) {
	    variant_allele_idx_base = variant_allele_idxs[write_variant_uidx];
	    cur_allele_ct = variant_allele_idxs[write_variant_uidx + 1] - variant_allele_idxs[write_variant_uidx];
	  }
	  char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
	  // possible todo: make number-to-string operations, strlen(), etc.
	  //   happen only once per variant.
	  for (uint32_t test_idx = 0; test_idx < cur_reported_test_ct; ++test_idx) {
	    if (chr_col) {
	      cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
	    }
	    if (variant_bps) {
	      cswritep = uint32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
	    }
	    cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
	    if (ref_col) {
	      *cswritep++ = '\t';
	      cswritep = strcpya(cswritep, cur_alleles[0]);
	    }
	    if (alt1_col) {
	      *cswritep++ = '\t';
	      cswritep = strcpya(cswritep, cur_alleles[1]);
	    }
	    if (alt_col) {
	      *cswritep++ = '\t';
	      for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
		if (cswrite(&css, &cswritep)) {
		  goto glm_linear_ret_WRITE_FAIL;
		}
		cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
	      }
	      --cswritep;
	    }
	    if (alt_ct_col) {
	      *cswritep++ = '\t';
	      cswritep = dtoa_g(auxp->alt_dosage, cswritep);
	    }
	    if (tot_allele_col) {
	      *cswritep++ = '\t';
	      cswritep = uint32toa(auxp->allele_obs_ct, cswritep);
	    }
	    if (alt_freq_col) {
	      *cswritep++ = '\t';
	      cswritep = dtoa_g(auxp->alt_dosage / ((double)auxp->allele_obs_ct), cswritep);
	    }
	    if (mach_r2_col) {
	      *cswritep++ = '\t';
	      if (!suppress_mach_r2) {
	        cswritep = dtoa_g(auxp->mach_r2, cswritep);
	      } else {
		cswritep = strcpya(cswritep, "NA");
	      }
	    }
	    if (test_col) {
	      *cswritep++ = '\t';
	      cswritep = strcpya(cswritep, cur_test_names[test_idx]);
	    }
	    if (nobs_col) {
	      *cswritep++ = '\t';
	      cswritep = uint32toa(auxp->sample_obs_ct, cswritep);
	    }
	    double pval = -9;
	    double tstat = 0.0;
	    if ((!cur_constraint_ct) || (test_idx != primary_reported_test_idx)) {
	      double beta = *beta_se_iter++;
	      double se = *beta_se_iter++;
	      if (!is_invalid) {
		tstat = beta / se;
		pval = calc_tprob(tstat, auxp->sample_obs_ct - cur_predictor_ct);
	      }
	      if (beta_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  cswritep = dtoa_g(beta, cswritep);
		} else {
		  cswritep = strcpya(cswritep, "NA");
		}
	      }
	      if (se_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  cswritep = dtoa_g(se, cswritep);
		} else {
		  cswritep = strcpya(cswritep, "NA");
		}
	      }
	      if (ci_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  const double ci_halfwidth = ci_zt * se;
		  cswritep = dtoa_g(beta - ci_halfwidth, cswritep);
		  *cswritep++ = '\t';
		  cswritep = dtoa_g(beta + ci_halfwidth, cswritep);
		} else {
		  cswritep = strcpya(cswritep, "NA\tNA");
		}
	      }
	      if (t_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  cswritep = dtoa_g(tstat, cswritep);
		} else {
		  cswritep = strcpya(cswritep, "NA");
		}
	      }
	    } else {
	      // joint test: use (currently approximate) F-test instead of T
	      // test
	      // beta_se_iter = &(beta_se_iter[2]);
	      if (beta_col) {
		cswritep = memcpyl3a(cswritep, "\tNA");
	      }
	      if (se_col) {
		cswritep = memcpyl3a(cswritep, "\tNA");
	      }
	      if (ci_col) {
		cswritep = strcpya(cswritep, "\tNA\tNA");
	      }
	      if (t_col) {
		*cswritep++ = '\t';
		if (!is_invalid) {
		  cswritep = dtoa_g(primary_se, cswritep);
		} else {
		  cswritep = strcpya(cswritep, "NA");
		}
	      }
	      // could avoid recomputing
	      if (!is_invalid) {
		pval = chiprob_p(primary_se, cur_constraint_ct);
	      }
	    }
	    if (p_col) {
	      *cswritep++ = '\t';
	      if (!is_invalid) {
		cswritep = dtoa_g(MAXV(pval, output_min_p), cswritep);
	      } else {
		cswritep = strcpya(cswritep, "NA");
	      }
	    }
	    append_binary_eoln(&cswritep);
	    if (cswrite(&css, &cswritep)) {
	      goto glm_linear_ret_WRITE_FAIL;
	    }
	    if (test_idx == primary_reported_test_idx) {
	      if (cur_pval_write) {
		cur_pval_write[variant_bidx] = pval;
	      }
	      if (cur_chisq_write) {
		if (!is_invalid) {
		  if (!cur_constraint_ct) {
		    cur_chisq_write[variant_bidx] = tstat * tstat;
		  } else {
		    cur_chisq_write[variant_bidx] = primary_se;
		  }
		} else {
		  cur_chisq_write[variant_bidx] = -9;
		}
	      }
	    }
	  }
	}
      }
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
      prev_block_variant_ct = cur_block_variant_ct;
      variant_idx += cur_block_variant_ct;
      // crucially, this is independent of the pgen_reader_t block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (cswrite_close_null(&css, cswritep)) {
      goto glm_linear_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
    LOGPRINTF("Results written to %s .\n", outname);
    bigstack_reset(bigstack_mark);
  }
  while (0) {
  glm_linear_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  glm_linear_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  glm_linear_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  glm_linear_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  glm_linear_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  glm_linear_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 glm_linear_ret_1:
  threads3z_cleanup(&ts, &g_cur_block_variant_ct);
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t glm_main(const uintptr_t* orig_sample_include, const char* sample_ids, const char* sids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const pheno_col_t* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const glm_info_t* glm_info_ptr, const adjust_info_t* adjust_info_ptr, const aperm_t* aperm_ptr, const char* local_covar_fname, const char* local_pvar_fname, const char* local_psam_fname, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t orig_covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t xchr_model, double ci_size, double vif_thresh, double pfilter, double output_min_p, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  gzFile gz_local_covar_file = nullptr;
  gz_token_stream_t gts;
  gz_token_stream_preinit(&gts);
  pglerr_t reterr = kPglRetSuccess;
  {
    if (!pheno_ct) {
      logerrprint("Error: No phenotypes loaded.\n");
      goto glm_main_ret_INCONSISTENT_INPUT;
    }
    if (orig_sample_ct < 2) {
      logerrprint("Error: --glm requires at least two samples.\n");
      goto glm_main_ret_INCONSISTENT_INPUT;
    }
    if (!orig_variant_ct) {
      logerrprint("Error: --glm requires at least one variant.\n");
      goto glm_main_ret_INCONSISTENT_INPUT;
    }
    // common linear/logistic initialization
    const uintptr_t* early_variant_include = orig_variant_include;
    uint32_t* local_sample_uidx_order = nullptr;
    uintptr_t* local_variant_include = nullptr;
    char* local_loadbuf = nullptr;
    uint32_t variant_ct = orig_variant_ct;
    uint32_t local_sample_ct = 0;
    uint32_t local_variant_ctl = 0;
    uint32_t local_covar_ct = 0;
    uint32_t local_loadbuf_size = 0;
    if (local_covar_fname) {
      reterr = glm_local_init(local_covar_fname, local_pvar_fname, local_psam_fname, sample_ids, cip, variant_bps, variant_ids, glm_info_ptr, raw_sample_ct, max_sample_id_blen, raw_variant_ct, &orig_sample_include, &sex_nm, &sex_male, &early_variant_include, &orig_sample_ct, &variant_ct, &gz_local_covar_file, &local_sample_uidx_order, &local_variant_include, &local_sample_ct, &local_variant_ctl, &local_covar_ct);
      if (reterr) {
	goto glm_main_ret_1;
      }
      uint64_t ullii = local_sample_ct;
      if (glm_info_ptr->local_cat_ct) {
	ullii *= 1 + int_slen(glm_info_ptr->local_cat_ct);
      } else {
	// permit 24 characters per floating point number instead of 16, since
	// some tools dump 15-17 significant digits
	ullii *= 24 * (local_covar_ct + ((glm_info_ptr->flags / kfGlmLocalOmitLast) & 1));
      }
      // +2 bytes for null terminator, \r\n; 1 more so we can detect gzgets
      // hitting the limit
      ullii += 3;
      if (ullii > kMaxLongLine) {
	logerrprint("Error: Too many samples/covariates for --glm local-covar=.\n");
	goto glm_main_ret_MALFORMED_INPUT;
      }
      if (ullii < kMaxMediumLine) {
	ullii = kMaxMediumLine; // may as well unconditionally support this
      }
      local_loadbuf_size = ullii;
      if (bigstack_alloc_c(local_loadbuf_size, &local_loadbuf)) {
	goto glm_main_ret_NOMEM;
      }
      local_loadbuf[local_loadbuf_size - 1] = ' ';
    }
    
    const glm_flags_t glm_flags = glm_info_ptr->flags;
    g_glm_flags = glm_flags;
    g_dosage_presents = nullptr;
    g_dosage_val_bufs = nullptr;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    const uint32_t perm_adapt = (glm_flags / kfGlmPerm) & 1;
    const uint32_t perms_total = perm_adapt? aperm_ptr->max : glm_info_ptr->mperm_ct;
    // [output prefix].[pheno name].glm.[linear/logistic]{.perm, .mperm}{.zst}
    uint32_t pheno_name_blen_capacity = kPglFnamesize - 14 - (4 * output_zst) - (uintptr_t)(outname_end - outname);
    if (perms_total) {
      pheno_name_blen_capacity -= 6 - perm_adapt;
    }
    if (max_pheno_name_blen > pheno_name_blen_capacity) {
      logerrprint("Error: Phenotype name and/or --out parameter too long.\n");
      goto glm_main_ret_INCONSISTENT_INPUT;
    }
    *outname_end = '.';
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
    unsigned char* overflow_buf;
    uintptr_t* cur_sample_include;
    // synthetic categorical covariate name could be ~twice max ID length?
    if (bigstack_alloc_uc(kCompressStreamBlock + 2 * kMaxIdSlen + max_chr_blen + kMaxIdSlen + 512 + 2 * max_allele_slen, &overflow_buf) ||
	bigstack_alloc_ul(raw_sample_ctl, &cur_sample_include) ||
	bigstack_alloc_ui(raw_sample_ctl, &g_sample_include_cumulative_popcounts)) {
      goto glm_main_ret_NOMEM;
    }
    g_sample_include = cur_sample_include;
    g_cip = cip;
    g_variant_allele_idxs = variant_allele_idxs;
    
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uint32_t max_variant_ct = variant_ct;

    uint32_t x_start;
    uint32_t x_end;
    get_xymt_start_and_end(cip, kChrOffsetX, &x_start, &x_end);
    uint32_t y_start;
    uint32_t y_end;
    get_xymt_start_and_end(cip, kChrOffsetY, &y_start, &y_end);

    uintptr_t* sex_male_collapsed_buf = nullptr;
    int32_t x_code;
    uint32_t variant_ct_x = 0;
    uint32_t variant_ct_y = 0;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t sex_nm_ct = popcount_longs(sex_nm, raw_sample_ctl);
    const uint32_t male_ct = popcount_longs(sex_male, raw_sample_ctl);
    uint32_t add_sex_covar = !(glm_flags & kfGlmNoXSex);
    if (add_sex_covar && ((!male_ct) || (male_ct == sex_nm_ct))) {
      add_sex_covar = 0;
    }
    uintptr_t* cur_sample_include_y_buf = nullptr;
    if (domdev_present || (glm_flags & (kfGlmDominant | kfGlmRecessive))) {
      // dominant/recessive/genotypic/hethom suppress all chromosomes which
      // aren't fully diploid.  (could throw in a hack to permit chrX if
      // all samples are female?  i.e. synthesize a chr_info_t where
      // xymt_codes[0] is -2 and haploid_mask X bit is cleared)
      uintptr_t* variant_include_nohap = nullptr;
      const uint32_t chr_ct = cip->chr_ct;
      uint32_t removed_variant_ct = 0;
      for (uint32_t chr_fo_idx = 0; chr_fo_idx < chr_ct; ++chr_fo_idx) {
	const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	if (is_set(cip->haploid_mask, chr_idx)) {
	  const uint32_t variant_uidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
	  const uint32_t variant_uidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	  const uint32_t cur_chr_variant_ct = popcount_bit_idx(early_variant_include, variant_uidx_start, variant_uidx_end);
	  if (cur_chr_variant_ct) {
	    if (!removed_variant_ct) {
	      // no main-loop logic for excluding all haploid chromosomes, so
	      // make a full copy of early_variant_include and throw away our
	      // reference to the original
	      if (bigstack_alloc_ul(raw_variant_ctl, &variant_include_nohap)) {
		goto glm_main_ret_NOMEM;
	      }
	      memcpy(variant_include_nohap, early_variant_include, raw_variant_ctl * sizeof(intptr_t));
	    }
	    clear_bits_nz(variant_uidx_start, variant_uidx_end, variant_include_nohap);
	    removed_variant_ct += cur_chr_variant_ct;
	  }
	}
      }
      if (removed_variant_ct) {
	if (variant_ct == removed_variant_ct) {
	  logerrprint("Error: No variants remaining for --glm ('dominant', 'recessive', 'genotypic',\nand 'hethom' only operate on diploid data).\n");
	  goto glm_main_ret_INCONSISTENT_INPUT;
	}
	variant_ct -= removed_variant_ct;
	early_variant_include = variant_include_nohap;
	max_variant_ct = variant_ct;
      }
    } else {
      if (xymt_exists(cip, kChrOffsetX, &x_code)) {
	variant_ct_x = count_chr_variants_unsafe(early_variant_include, cip, (uint32_t)x_code);
	// --xchr-model 0 now only suppresses chrX.
	if (xchr_model) {
	  if (variant_ct_x) {
	    if (bigstack_alloc_ul(BITCT_TO_WORDCT(orig_sample_ct), &sex_male_collapsed_buf)) {
	      goto glm_main_ret_NOMEM;
	    }
	  }
	} else {
	  max_variant_ct -= variant_ct_x;
	  if (!max_variant_ct) {
	    logerrprint("Error: No variants remaining for --glm, due to --xchr-model 0.\n");
	    goto glm_main_ret_INCONSISTENT_INPUT;
	  }
	}
      }
      int32_t y_code;
      if (xymt_exists(cip, kChrOffsetY, &y_code)) {
	variant_ct_y = count_chr_variants_unsafe(early_variant_include, cip, (uint32_t)y_code);
	if (variant_ct_y) {
	  if (!male_ct) {
	    logprint("--glm: Skipping chrY since there are no males.\n");
	    max_variant_ct -= variant_ct_y;
	    if (!max_variant_ct) {
	      logerrprint("Error: No variants remaining for --glm.\n");
	      goto glm_main_ret_INCONSISTENT_INPUT;
	    }
	  } else if (male_ct < orig_sample_ct) {
	    // may as well check for only-chrY special case
	    if (max_variant_ct != variant_ct_y) {
	      if (bigstack_alloc_ul(raw_sample_ctl, &cur_sample_include_y_buf)) {
		// covar_include_y allocation postponed since raw_covar_ct not
		// yet known
		goto glm_main_ret_NOMEM;
	      }
	    } else {
	      orig_sample_include = sex_male;
	      orig_sample_ct = male_ct;
	    }
	  }
	}
      }
    }
    if (add_sex_covar && (!variant_ct_x) && (!(glm_flags & kfGlmSex))) {
      add_sex_covar = 0;
    }
    g_sex_male_collapsed = sex_male_collapsed_buf;
    uint32_t raw_covar_ct = orig_covar_ct + local_covar_ct;
    if (glm_info_ptr->condition_varname || glm_info_ptr->condition_list_fname || local_covar_ct || add_sex_covar) {
      uint32_t condition_ct = 0;
      pheno_col_t* new_covar_cols;
      char* new_covar_names;
      uintptr_t new_max_covar_name_blen = max_covar_name_blen;
      if (add_sex_covar && (new_max_covar_name_blen < 4)) {
	new_max_covar_name_blen = 4;
      }
      if (local_covar_ct && (new_max_covar_name_blen < 6 + int_slen(local_covar_ct + 1))) {
	new_max_covar_name_blen = 6 + int_slen(local_covar_ct + 1);
      }
      if (glm_info_ptr->condition_varname || glm_info_ptr->condition_list_fname) {
	assert(g_bigstack_end == bigstack_end_mark);
	// reserve space for condition-list worst case (roughly sqrt(2^31)),
	// since that's relatively small
	const uint32_t condition_ct_max = 46338;
	uint32_t* condition_uidxs;
	if (bigstack_end_alloc_ui(condition_ct_max, &condition_uidxs)) {
	  goto glm_main_ret_NOMEM;
	}
	if (glm_info_ptr->condition_varname) {
	  int32_t ii = get_variant_uidx_without_htable(glm_info_ptr->condition_varname, variant_ids, orig_variant_include, orig_variant_ct);
	  if (ii >= 0) {
	    condition_uidxs[0] = (uint32_t)ii;
	    condition_ct = 1;
	    const uint32_t condition_blen = strlen(glm_info_ptr->condition_varname) + 1;
	    // drop "CSNP" column name for sanity's sake
	    if (new_max_covar_name_blen < condition_blen) {
	      new_max_covar_name_blen = condition_blen;
	    }
	  } else {
	    if (ii == -2) {
	      LOGERRPRINTFWW("Error: Duplicate --condition variant ID '%s'.\n", glm_info_ptr->condition_varname);
	      goto glm_main_ret_INVALID_CMDLINE;
	    }
	    LOGERRPRINTFWW("Warning: --condition variant ID '%s' not found.\n", glm_info_ptr->condition_varname);
	  }
	} else {
	  // 1. (re)construct variant ID hash table
	  uintptr_t* already_seen;
	  if (bigstack_calloc_ul(raw_variant_ctl, &already_seen)) {
	    goto glm_main_ret_NOMEM;
	  }
	  uint32_t* variant_id_htable = nullptr;
	  uint32_t variant_id_htable_size;
	  reterr = alloc_and_populate_id_htable_mt(orig_variant_include, variant_ids, orig_variant_ct, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size);
	  if (reterr) {
	    goto glm_main_ret_1;
	  }

	  // 2. iterate through --condition-list file, make sure no IDs are
	  //    duplicate in loaded fileset, warn about duplicates in
	  //    --condition-list file
	  reterr = gz_token_stream_init(glm_info_ptr->condition_list_fname, &gts, g_textbuf);
	  if (reterr) {
	    goto glm_main_ret_1;
	  }
	  uintptr_t skip_ct = 0;
	  uintptr_t duplicate_ct = 0;
	  uint32_t token_slen;
	  while (1) {
	    char* token_start = gz_token_stream_advance(&gts, &token_slen);
	    if (!token_start) {
	      break;
	    }
	    uint32_t cur_variant_uidx = variant_id_dupflag_htable_find(token_start, variant_ids, variant_id_htable, token_slen, variant_id_htable_size, max_variant_id_slen);
	    if (cur_variant_uidx >> 31) {
	      if (cur_variant_uidx != 0xffffffffU) {
		LOGERRPRINTFWW("Error: --condition-list variant ID '%s' appears multiple times.\n", variant_ids[cur_variant_uidx & 0x7fffffff]);
		goto glm_main_ret_INCONSISTENT_INPUT;
	      }
	      ++skip_ct;
	    } else if (is_set(already_seen, cur_variant_uidx)) {
	      ++duplicate_ct;
	    } else {
	      if (condition_ct == condition_ct_max) {
		logerrprint("Error: Too many --condition-list variant IDs.\n");
		goto glm_main_ret_MALFORMED_INPUT;
	      }
	      set_bit(cur_variant_uidx, already_seen);
	      condition_uidxs[condition_ct++] = cur_variant_uidx;
	      if (new_max_covar_name_blen <= token_slen) {
		new_max_covar_name_blen = token_slen + 1;
	      }
	    }
	  }
	  if (token_slen) {
	    // error code
	    if (token_slen == 0xffffffffU) {
	      logerrprint("Error: Excessively long ID in --condition-list file.\n");
	      goto glm_main_ret_MALFORMED_INPUT;
	    }
	    goto glm_main_ret_READ_FAIL;
	  }
	  if (gz_token_stream_close(&gts)) {
	    goto glm_main_ret_READ_FAIL;
	  }
	  if (skip_ct || duplicate_ct) {
	    if (skip_ct && duplicate_ct) {
	      LOGERRPRINTFWW("Warning: %" PRIuPTR " --condition-list variant ID%s not found, and %" PRIuPTR " duplicate ID%s present.\n", skip_ct, (skip_ct == 1)? "" : "s", duplicate_ct, (duplicate_ct == 1)? "" : "s");
	    } else if (skip_ct) {
	      LOGERRPRINTF("Warning: %" PRIuPTR " --condition-list variant ID%s not found.\n", skip_ct, (skip_ct == 1)? "" : "s");
	    } else {
	      LOGERRPRINTF("Warning: %" PRIuPTR " duplicate --condition-list variant ID%s present.\n", duplicate_ct, (duplicate_ct == 1)? "" : "s");
	    }
	  }
	  LOGPRINTF("--condition-list: %u variant ID%s loaded.\n", condition_ct, (condition_ct == 1)? "" : "s");

	  // free hash table and duplicate tracker
	  bigstack_reset(already_seen);
	}
	raw_covar_ct += condition_ct;
	new_covar_cols = (pheno_col_t*)bigstack_alloc((raw_covar_ct + add_sex_covar) * sizeof(pheno_col_t));
	if ((!new_covar_cols) ||
	    bigstack_alloc_c((raw_covar_ct + add_sex_covar) * new_max_covar_name_blen, &new_covar_names)) {
	  goto glm_main_ret_NOMEM;
	}
        if (condition_ct) {
	  bigstack_end_set(condition_uidxs);
	  uintptr_t* genovec;
	  uintptr_t* dosage_present;
	  dosage_t* dosage_vals;
	  if (bigstack_end_alloc_ul(QUATERCT_TO_WORDCT(raw_sample_ct), &genovec) ||
	      bigstack_end_alloc_ul(raw_sample_ctl, &dosage_present) ||
	      bigstack_end_alloc_dosage(raw_sample_ct, &dosage_vals)) {
	    goto glm_main_ret_NOMEM;
	  }
	  pgr_clear_ld_cache(simple_pgrp);
	  for (uint32_t condition_idx = 0; condition_idx < condition_ct; ++condition_idx) {
	    const uint32_t cur_variant_uidx = condition_uidxs[condition_idx];
	    uint32_t dosage_ct;
	    uint32_t is_explicit_alt1;
	    reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(nullptr, nullptr, raw_sample_ct, cur_variant_uidx, simple_pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
	    if (reterr) {
	      if (reterr == kPglRetMalformedInput) {
		logerrprint("Error: Malformed .pgen file.\n");
	      }
	      goto glm_main_ret_1;
	    }
	    pheno_col_t* cur_covar_col = &(new_covar_cols[local_covar_ct + condition_idx]);
	    uintptr_t* cur_nonmiss;
	    double* cur_covar_vals;
	    if (bigstack_alloc_ul(raw_sample_ctl, &cur_nonmiss) ||
		bigstack_alloc_d(raw_sample_ct, &cur_covar_vals)) {
	      goto glm_main_ret_NOMEM;
	    }
	    cur_covar_col->category_names = nullptr;
	    cur_covar_col->nonmiss = cur_nonmiss;
	    cur_covar_col->data.qt = cur_covar_vals;
	    cur_covar_col->type_code = kPhenoDtypeQt;
	    cur_covar_col->nonnull_category_ct = 0;
	    genoarr_to_nonmissing(genovec, raw_sample_ct, cur_nonmiss);
	    genoarr_to_doubles(genovec, raw_sample_ct, cur_covar_vals);
	    if (dosage_ct) {
	      uint32_t sample_uidx = 0;
	      for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_uidx) {
		next_set_unsafe_ck(dosage_present, &sample_uidx);
	        cur_covar_vals[sample_uidx] = kRecipDosageMid * ((int32_t)((uint32_t)dosage_vals[dosage_idx]));
	      }
	      bitvec_or(dosage_present, raw_sample_ctl, cur_nonmiss);
	    }
	    if (glm_flags & kfGlmConditionDominant) {
	      for (uint32_t sample_uidx = 0; sample_uidx < raw_sample_ct; ++sample_uidx) {
		if (cur_covar_vals[sample_uidx] > 1.0) {
		  cur_covar_vals[sample_uidx] = 1.0;
		}
	      }
	    } else if (glm_flags & kfGlmConditionRecessive) {
	      for (uint32_t sample_uidx = 0; sample_uidx < raw_sample_ct; ++sample_uidx) {
		double dxx = cur_covar_vals[sample_uidx];
		if (dxx <= 1.0) {
		  dxx = 0;
		} else {
		  dxx -= 1.0;
		}
		cur_covar_vals[sample_uidx] = dxx;
	      }
	    }
	    strcpy(&(new_covar_names[(local_covar_ct + condition_idx) * new_max_covar_name_blen]), variant_ids[cur_variant_uidx]);
	  }
	  bigstack_end_reset(bigstack_end_mark);
	}
      } else {
	new_covar_cols = (pheno_col_t*)bigstack_alloc((raw_covar_ct + add_sex_covar) * sizeof(pheno_col_t));
	if ((!new_covar_cols) ||
	    bigstack_alloc_c((raw_covar_ct + add_sex_covar) * new_max_covar_name_blen, &new_covar_names)) {
	  goto glm_main_ret_NOMEM;
	}
      }
      if (local_covar_ct) {
	memcpy(new_covar_cols, covar_cols, local_covar_ct * sizeof(pheno_col_t));
      }
      memcpy(&(new_covar_cols[condition_ct + local_covar_ct]), covar_cols, orig_covar_ct * sizeof(pheno_col_t));
      const char* covar_names_read_iter = covar_names;
      // bugfix (11 May 2017): local covar names come before, not after,
      //   --condition{-list} covar names
      char* covar_names_write_iter = new_covar_names;
      for (uint32_t local_covar_idx = 0; local_covar_idx < local_covar_ct; ++local_covar_idx) {
	memcpy(covar_names_write_iter, "LOCAL", 5);
	char* name_end = uint32toa(local_covar_idx + 1, &(covar_names_write_iter[5]));
	*name_end = '\0';
	new_covar_cols[local_covar_idx].type_code = kPhenoDtypeOther;
	new_covar_cols[local_covar_idx].nonmiss = nullptr;
	covar_names_write_iter = &(covar_names_write_iter[new_max_covar_name_blen]);
      }
      covar_names_write_iter = &(covar_names_write_iter[condition_ct * new_max_covar_name_blen]);
      for (uint32_t old_covar_idx = 0; old_covar_idx < orig_covar_ct; ++old_covar_idx) {
	strcpy(covar_names_write_iter, covar_names_read_iter);
	covar_names_read_iter = &(covar_names_read_iter[max_covar_name_blen]);
	covar_names_write_iter = &(covar_names_write_iter[new_max_covar_name_blen]);
      }
      if (add_sex_covar) {
	pheno_col_t* new_sex_col = &(new_covar_cols[raw_covar_ct++]);
	double* sex_covar_vals;
	if (bigstack_alloc_d(raw_sample_ct, &sex_covar_vals)) {
	  goto glm_main_ret_NOMEM;
	}
	uint32_t sample_uidx = 0;
	for (uint32_t sample_idx = 0; sample_idx < orig_sample_ct; ++sample_idx, ++sample_uidx) {
	  next_set_unsafe_ck(sex_nm, &sample_uidx);
	  // 1/2 instead of 1/0 coding; user shouldn't have to worry about
	  // signs changing when they use --sex instead of using the sex column
	  // from a .bim/.psam file
	  sex_covar_vals[sample_uidx] = (double)((int32_t)(2 - is_set(sex_male, sample_uidx)));
	}
	new_sex_col->category_names = nullptr;
	// const_cast
	new_sex_col->nonmiss = (uintptr_t*)((uintptr_t)sex_nm);
	new_sex_col->data.qt = sex_covar_vals;
	new_sex_col->type_code = kPhenoDtypeQt;
	new_sex_col->nonnull_category_ct = 0;
	strcpy(covar_names_write_iter, "SEX");
      }
      covar_cols = new_covar_cols;
      covar_names = new_covar_names;
      max_covar_name_blen = new_max_covar_name_blen;
    }
    const uint32_t raw_covar_ctl = BITCT_TO_WORDCT(raw_covar_ct);
    uintptr_t* initial_covar_include = nullptr;
    uintptr_t* covar_include = nullptr;
    uintptr_t* cur_sample_include_x_buf = nullptr;
    uintptr_t* covar_include_x = nullptr;
    uint32_t covar_max_nonnull_cat_ct = 0;
    if (raw_covar_ctl) {
      if (bigstack_alloc_ul(raw_covar_ctl, &initial_covar_include) ||
	  bigstack_alloc_ul(raw_covar_ctl, &covar_include)) {
	goto glm_main_ret_NOMEM;
      }
      fill_ulong_zero(raw_covar_ctl, initial_covar_include);
      for (uint32_t covar_uidx = 0; covar_uidx < raw_covar_ct; ++covar_uidx) {
	const pheno_col_t* cur_covar_col = &(covar_cols[covar_uidx]);
	if (cur_covar_col->type_code != kPhenoDtypeOther) {
	  if (!is_const_covar(cur_covar_col, orig_sample_include, orig_sample_ct)) {
	    set_bit(covar_uidx, initial_covar_include);
	    if (cur_covar_col->type_code == kPhenoDtypeCat) {
	      if (cur_covar_col->nonnull_category_ct > covar_max_nonnull_cat_ct) {
		covar_max_nonnull_cat_ct = cur_covar_col->nonnull_category_ct;
	      }
	    }
	  } else {
	    LOGERRPRINTF("Warning: Excluding constant covariate '%s' from --glm.\n", &(covar_names[covar_uidx * max_covar_name_blen]));
	  }
	} else {
	  // local covariate, always include
	  set_bit(covar_uidx, initial_covar_include);
	}
      }
      if (covar_max_nonnull_cat_ct && (glm_info_ptr->parameters_range_list.name_ct || glm_info_ptr->tests_range_list.name_ct)) {
	// todo: permit this, and automatically expand a single parameter index
	// referring to a categorical covariate into the appropriate range of
	// final predictor indices
	logerrprint("Error: --parameters/--tests cannot be used directly with categorical\ncovariates; expand them into binary covariates with --split-cat-pheno first.\n");
	goto glm_main_ret_INCONSISTENT_INPUT;
      }
    }
    const uint32_t domdev_present_p1 = domdev_present + 1;
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    const uint32_t raw_predictor_ct = 2 + domdev_present + raw_covar_ct * (1 + add_interactions * domdev_present_p1);
    const uint32_t raw_predictor_ctl = BITCT_TO_WORDCT(raw_predictor_ct);
    const uint32_t first_covar_pred_uidx = 2 + domdev_present;
    uint32_t first_interaction_pred_uidx = 0;
    if (add_interactions) {
      first_interaction_pred_uidx = first_covar_pred_uidx + raw_covar_ct;
    }

    // TODO: --tests
    const uint32_t joint_test = domdev_present;
    g_constraints_con_major = nullptr;
    g_constraints_con_major_x = nullptr;
    g_constraints_con_major_y = nullptr;
    g_constraints_con_major_f = nullptr;
    g_constraints_con_major_x_f = nullptr;
    g_constraints_con_major_y_f = nullptr;
    
    uintptr_t* raw_parameter_subset = nullptr;
    g_parameter_subset = nullptr;
    g_parameter_subset_x = nullptr;
    g_parameter_subset_y = nullptr;
    if (glm_info_ptr->parameters_range_list.name_ct) {
      if (bigstack_calloc_ul(raw_predictor_ctl, &raw_parameter_subset) ||
	  bigstack_alloc_ul(raw_predictor_ctl, &g_parameter_subset) ||
	  bigstack_alloc_ul(raw_predictor_ctl, &g_parameter_subset_x) ||
	  bigstack_alloc_ul(raw_predictor_ctl, &g_parameter_subset_y)) {
	goto glm_main_ret_NOMEM;
      }
      raw_parameter_subset[0] = 1; // intercept (index 0) always included
      numeric_range_list_to_bitarr(&(glm_info_ptr->parameters_range_list), raw_predictor_ct, 0, 1, raw_parameter_subset);
      if (domdev_present && ((raw_parameter_subset[0] & 7) != 7)) {
	// this breaks the joint test
	logerrprint("Error: --parameters cannot exclude 1 or 2 when the 'genotypic' or 'hethom'\nmodifier is present.\n");
	goto glm_main_ret_INVALID_CMDLINE;
      }
      if (add_sex_covar && first_interaction_pred_uidx) {
	// special case: when add_sex_covar is true, the added sex covariate is
	// simply the last covariate, with predictor index
	// (first_interaction_pred_uidx - 1).  This lines up with --parameters
	// when interactions are not requested; but when they are, we have a
	// small reshuffle to do.
	uintptr_t* parameter_subset_reshuffle_buf;
	if (bigstack_calloc_ul(raw_predictor_ctl, &parameter_subset_reshuffle_buf)) {
	  goto glm_main_ret_NOMEM;
	}
	copy_bitarr_range(raw_parameter_subset, 0, 0, first_interaction_pred_uidx - 1, parameter_subset_reshuffle_buf);
	copy_bitarr_range(raw_parameter_subset, first_interaction_pred_uidx - 1, first_interaction_pred_uidx, raw_covar_ct * domdev_present_p1, parameter_subset_reshuffle_buf);
	const uint32_t first_sex_parameter_idx = first_interaction_pred_uidx - 1 + raw_covar_ct * domdev_present_p1;
	if (is_set(raw_parameter_subset, first_sex_parameter_idx)) {
	  set_bit(first_interaction_pred_uidx - 1, parameter_subset_reshuffle_buf);
	}
	if (is_set(raw_parameter_subset, first_sex_parameter_idx + 1)) {
	  set_bit(first_sex_parameter_idx + 1, parameter_subset_reshuffle_buf);
	}
	if (domdev_present && is_set(raw_parameter_subset, first_sex_parameter_idx + 2)) {
	  set_bit(first_sex_parameter_idx + 2, parameter_subset_reshuffle_buf);
	}
	memcpy(raw_parameter_subset, parameter_subset_reshuffle_buf, raw_predictor_ctl * sizeof(intptr_t));
	bigstack_reset(parameter_subset_reshuffle_buf);
      }
      // if there were any constant covariates, exclude them from
      // raw_parameter_subset
      // note that, if appended sex covariate is present at all, it is always
      // nonconstant
      uint32_t nonconst_covar_ct = 0;
      if (initial_covar_include) {
	nonconst_covar_ct = popcount_longs(initial_covar_include, raw_covar_ctl);
      }
      const uint32_t removed_covar_ct = raw_covar_ct - nonconst_covar_ct;
      uint32_t covar_uidx = 0;
      for (uint32_t removed_covar_idx = 0; removed_covar_idx < removed_covar_ct; ++removed_covar_idx, ++covar_uidx) {
	next_unset_unsafe_ck(initial_covar_include, &covar_uidx);
	clear_bit(first_covar_pred_uidx + covar_uidx, raw_parameter_subset);
	if (first_interaction_pred_uidx) {
	  const uint32_t geno_interaction_uidx = first_interaction_pred_uidx + covar_uidx * domdev_present_p1;
	  clear_bit(geno_interaction_uidx, raw_parameter_subset);
	  if (domdev_present) {
	    clear_bit(geno_interaction_uidx + 1, raw_parameter_subset);
	  }
	}
      }
      // if any loaded nonconstant covariates aren't referenced in
      // raw_parameter_subset, remove them from initial_covar_include
      covar_uidx = 0;
      for (uint32_t nonconst_covar_idx = 0; nonconst_covar_idx < nonconst_covar_ct; ++nonconst_covar_idx, ++covar_uidx) {
	next_set_unsafe_ck(initial_covar_include, &covar_uidx);
	uint32_t cur_covar_is_referenced = is_set(raw_parameter_subset, first_covar_pred_uidx + covar_uidx);
	if (add_interactions) {
	  cur_covar_is_referenced = cur_covar_is_referenced || is_set(raw_parameter_subset, first_interaction_pred_uidx + covar_uidx * domdev_present_p1);
	  if (domdev_present) {
	    cur_covar_is_referenced = cur_covar_is_referenced || is_set(raw_parameter_subset, first_interaction_pred_uidx + covar_uidx * 2 + 1);
	  }
	}
	if (!cur_covar_is_referenced) {
	  clear_bit(covar_uidx, initial_covar_include);
	}
      }
      // if your regression doesn't involve genotype data, you should be using
      // e.g. R, not plink...
      if ((!(raw_parameter_subset[0] & 2)) &&
	  ((!domdev_present) || (!(raw_parameter_subset[0] & 4))) &&
	  ((!add_interactions) || (!popcount_bit_idx(raw_parameter_subset, first_interaction_pred_uidx, raw_predictor_ct)))) {
	logerrprint("Error: --parameters must retain at least one dosage-dependent variable.\n");
	goto glm_main_ret_INCONSISTENT_INPUT;
      }
    }
    // computation of these counts moved here, since --parameters can reduce
    // the number of relevant covariates
    uint32_t initial_nonx_covar_ct = 0;
    if (initial_covar_include) {
      initial_nonx_covar_ct = popcount_longs(initial_covar_include, raw_covar_ctl);
    }
    uint32_t initial_y_covar_ct = 0;
    uintptr_t* covar_include_y = nullptr;
    if (!initial_nonx_covar_ct) {
      // bigstack_reset(initial_covar_include); // not ok with parameters
      initial_covar_include = nullptr;
      covar_include = nullptr;
    } else {
      initial_y_covar_ct = initial_nonx_covar_ct - (cur_sample_include_y_buf && add_sex_covar && is_set(initial_covar_include, raw_covar_ct - 1));
      if (add_sex_covar && (!(glm_flags & kfGlmSex))) {
	// may as well verify there's at least one non-x/non-y variant
	// (if only chrX and chrY present, don't allocate
	// cur_sample_include_x_buf, just make chrX the baseline instead)
	if (is_set(initial_covar_include, raw_covar_ct - 1) && (variant_ct != variant_ct_x + variant_ct_y)) {
	  if (bigstack_alloc_ul(raw_sample_ctl, &cur_sample_include_x_buf) ||
	      bigstack_alloc_ul(raw_covar_ctl, &covar_include_x)) {
	    goto glm_main_ret_NOMEM;
	  }
	  --initial_nonx_covar_ct;
	}
      }
      if (cur_sample_include_y_buf) {
	if (bigstack_alloc_ul(raw_covar_ctl, &covar_include_y)) {
	  goto glm_main_ret_NOMEM;
	}
      }
    }
    const uint32_t report_adjust = (adjust_info_ptr->flags & kfAdjustColAll);
    const uint32_t is_sometimes_firth = (glm_flags & (kfGlmFirthFallback | kfGlmFirth))? 1 : 0;
    const uint32_t is_always_firth = glm_flags & kfGlmFirth;
    const uint32_t glm_pos_col = glm_info_ptr->cols & kfGlmColPos;

    unsigned char* bigstack_mark2 = g_bigstack_base;
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      const pheno_col_t* cur_pheno_col = &(pheno_cols[pheno_idx]);
      const pheno_dtype_t dtype_code = cur_pheno_col->type_code;
      const char* cur_pheno_name = &(pheno_names[pheno_idx * max_pheno_name_blen]);
      if (dtype_code == kPhenoDtypeCat) {
	// todo: check if there are only two categories after linear-style
	// covariate QC, and automatically use ordinary logistic regression in
	// that case?  (need to indicate which category is treated as 'case'
	// and which is 'control'...)
	// longer-term todo: multinomial logistic regression?
	LOGPRINTFWW("--glm: Skipping categorical phenotype '%s'.\n", cur_pheno_name);
	continue;
      }

      bitvec_and_copy(orig_sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include);
      const uint32_t is_logistic = (dtype_code == kPhenoDtypeCc);
      uint32_t sample_ct = popcount_longs(cur_sample_include, raw_sample_ctl);
      if (is_logistic) {
	const uint32_t initial_case_ct = popcount_longs_intersect(cur_sample_include, cur_pheno_col->data.cc, raw_sample_ctl);
	if ((!initial_case_ct) || (initial_case_ct == sample_ct)) {
	  LOGPRINTFWW("--glm: Skipping case/control phenotype '%s' since all samples are %s.\n", cur_pheno_name, initial_case_ct? "cases" : "controls");
	  continue;
	}
      } else {
	if (is_const_covar(cur_pheno_col, cur_sample_include, sample_ct)) {
	  LOGPRINTFWW("--glm: Skipping constant quantitative phenotype '%s'.\n", cur_pheno_name);
	  continue;
	}
      }
      uint32_t covar_ct = 0;
      uint32_t extra_cat_ct = 0;
      uint32_t separation_warning = 0;
      bigstack_double_reset(bigstack_mark2, bigstack_end_mark);
      /*
      for (uint32_t uii = 0; uii < raw_covar_ct; ++uii) {
	const pheno_col_t* cur_covar_col = &(covar_cols[uii]);
	for (uint32_t sample_uidx = 0; sample_uidx < raw_sample_ct; ++sample_uidx) {
	  printf("%g ", cur_covar_col->data.qt[sample_uidx]);
	}
	printf("\n\n");
      }
      */
      if (initial_nonx_covar_ct) {
	if (glm_determine_covars(is_logistic? cur_pheno_col->data.cc : nullptr, initial_covar_include, covar_cols, raw_sample_ct, raw_covar_ctl, initial_nonx_covar_ct, covar_max_nonnull_cat_ct, is_sometimes_firth, cur_sample_include, covar_include, &sample_ct, &covar_ct, &extra_cat_ct, &separation_warning)) {
	  goto glm_main_ret_NOMEM;
	}
      }
      uint32_t predictor_ct = 2 + domdev_present + (covar_ct + extra_cat_ct) * (1 + add_interactions * domdev_present_p1);
      if (raw_parameter_subset) {
	collapse_parameter_subset(covar_include, raw_parameter_subset, domdev_present, raw_covar_ct, covar_ct, add_interactions, g_parameter_subset, &predictor_ct);
      }
      if (sample_ct <= predictor_ct) {
	LOGERRPRINTFWW("Warning: Skipping --glm regression on phenotype '%s' since # samples <= # predictor columns.\n", cur_pheno_name);
	if (separation_warning) {
	  logerrprint("(Quasi-)separated covariate(s) were present.  Try removing inappropriate\ncovariates, and/or using Firth logistic regression.\n");
	}
	continue;
      }
#ifdef __LP64__
      if (round_up_pow2(sample_ct, 4) * ((uint64_t)predictor_ct) > 0x7fffffff) {
	// todo: remove this constraint in LAPACK_ILP64 case?
	LOGERRPRINTFWW("Warning: Skipping --glm regression on phenotype '%s' since there are too many\nsamples or predictors (internal matrices limited to ~2^31 entries).\n", cur_pheno_name);
	continue;
      }
#endif
      uint32_t case_ct = 0;
      if (is_logistic) {
	case_ct = popcount_longs_intersect(cur_sample_include, cur_pheno_col->data.cc, raw_sample_ctl);
	if ((!case_ct) || (case_ct == sample_ct)) {
	  LOGPRINTFWW("--glm: Skipping case/control phenotype '%s' since all remaining samples are %s.\n", cur_pheno_name, case_ct? "cases" : "controls");
	  // without any e.g. cases in the dataset, every single covariate
	  // should fail the separation check, so covar_ct should be zero here
	  assert(!covar_ct);
	  continue;
	}
	if (sample_ct < 10 * predictor_ct) {
	  LOGERRPRINTFWW("Warning: --glm remaining sample count is less than 10x predictor count for case/control phenotype '%s'.\n", cur_pheno_name);
	}
      } else {
	// verify phenotype is still nonconstant
	if (is_const_covar(cur_pheno_col, cur_sample_include, sample_ct)) {
	  LOGPRINTFWW("--glm: Skipping quantitative phenotype '%s' since phenotype is constant for all remaining samples.\n", cur_pheno_name);
	  continue;
	}
      }
      if (covar_ct < initial_nonx_covar_ct) {
	uint32_t covar_uidx = 0;
	for (uint32_t covar_idx = 0; covar_idx < initial_nonx_covar_ct; ++covar_idx, ++covar_uidx) {
	  next_set_unsafe_ck(initial_covar_include, &covar_uidx);
	  if (!is_set(covar_include, covar_uidx)) {
	    LOGERRPRINTFWW("Warning: %sot including covariate '%s' in --glm regression on phenotype '%s'.\n", cur_sample_include_x_buf? (cur_sample_include_y_buf? "Outside of chrX, n" : "Outside of chrX and chrY, n") : (cur_sample_include_y_buf? "Outside of chrY, n" : "N"), &(covar_names[covar_uidx * max_covar_name_blen]), cur_pheno_name);
	  }
	}
      }
      if (separation_warning) {
	logerrprint("(Quasi-)separated covariate(s) were present.  Try removing inappropriate\ncovariates, and/or using Firth logistic regression.\n");
      }

      // cur_sample_include_x == nullptr: chrX uses same samples and covariates
      //   as the rest of the genome.  sample_ct_x always zero to force most
      //   chrX-specific initialization to be skipped (exception:
      //   sex_male_collapsed, needed for allele count/freq reporting)
      // cur_sample_include_x non-null: if sample_ct_x == 0, we skip the entire
      //   chromosome.  otherwise, we have different covariates than the rest
      //   of the genome.
      uintptr_t* cur_sample_include_x = cur_sample_include_x_buf;
      uint32_t sample_ct_x = 0;
      uint32_t covar_ct_x = 0;
      uint32_t extra_cat_ct_x = 0;
      uint32_t predictor_ct_x = 0;
      uint32_t x_samples_are_different = 0;
      if (cur_sample_include_x) {
	bitvec_and_copy(orig_sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include_x);
        uint32_t separation_warning_x = 0;
	if (glm_determine_covars(is_logistic? cur_pheno_col->data.cc : nullptr, initial_covar_include, covar_cols, raw_sample_ct, raw_covar_ctl, initial_nonx_covar_ct + 1, covar_max_nonnull_cat_ct, is_sometimes_firth, cur_sample_include_x, covar_include_x, &sample_ct_x, &covar_ct_x, &extra_cat_ct_x, &separation_warning_x)) {
	  goto glm_main_ret_NOMEM;
	}
	x_samples_are_different = (sample_ct_x != sample_ct) || (!are_all_words_identical(cur_sample_include, cur_sample_include_x, raw_sample_ctl));
	if ((!x_samples_are_different) && (covar_ct == covar_ct_x) && are_all_words_identical(covar_include, covar_include_x, raw_covar_ctl)) {
	  LOGPRINTFWW("Note: chrX samples and covariate(s) in --glm regression on phenotype '%s' are the same as that for the rest of the genome.\n", cur_pheno_name);
	  sample_ct_x = 0;
	  cur_sample_include_x = nullptr;
	} else {
	  if (!sample_ct_x) {
	    LOGERRPRINTFWW("Warning: Skipping chrX in --glm regression on phenotype '%s'.\n", cur_pheno_name);
	  } else {
            predictor_ct_x = 2 + domdev_present + (covar_ct_x + extra_cat_ct_x) * (1 + add_interactions * domdev_present_p1);
	    if (raw_parameter_subset) {
	      collapse_parameter_subset(covar_include, raw_parameter_subset, domdev_present, raw_covar_ct, covar_ct_x, add_interactions, g_parameter_subset_x, &predictor_ct_x);
	    }
	    if (sample_ct_x <= predictor_ct_x) {
	      LOGERRPRINTFWW("Warning: Skipping chrX in --glm regression on phenotype '%s', since # remaining samples <= # predictor columns.\n", cur_pheno_name);
	      sample_ct_x = 0;
#ifdef __LP64__
	    } else if (round_up_pow2(sample_ct_x, 4) * ((uint64_t)predictor_ct_x) > 0x7fffffff) {
	      LOGERRPRINTFWW("Warning: Skipping chrX in --glm regression on phenotype '%s', since there are\ntoo many samples or predictors (internal matrices limited to ~2^31 entries).\n", cur_pheno_name);
	      sample_ct_x = 0;
#endif
	    } else if (is_logistic) {
	      const uint32_t case_ct_x = popcount_longs_intersect(cur_sample_include_x, cur_pheno_col->data.cc, raw_sample_ctl);
	      if ((!case_ct_x) || (case_ct_x == sample_ct_x)) {
		LOGERRPRINTFWW("Warning: Skipping chrX in --glm regression on phenotype '%s', since all remaining samples are %s.\n", cur_pheno_name, case_ct_x? "cases" : "controls");
		sample_ct_x = 0;
	      }
	    } else {
	      if (is_const_covar(cur_pheno_col, cur_sample_include_x, sample_ct_x)) {
		LOGERRPRINTFWW("Warning: Skipping chrX in --glm regression on phenotype '%s', since phenotype is constant for all remaining samples.\n", cur_pheno_name);
		sample_ct_x = 0;
	      }
	    }
	    if (sample_ct_x && (covar_ct_x < initial_nonx_covar_ct + 1)) {
	      uint32_t covar_uidx = 0;
	      for (uint32_t covar_idx = 0; covar_idx < covar_ct_x; ++covar_idx, ++covar_uidx) {
		next_set_unsafe_ck(initial_covar_include, &covar_uidx);
		if (!is_set(covar_include_x, covar_uidx)) {
		  LOGERRPRINTFWW("Warning: On chrX, not including covariate '%s' in --glm regression on phenotype '%s'.\n", &(covar_names[covar_uidx * max_covar_name_blen]), cur_pheno_name);
		}
	      }
	    }
	  }
	  if (separation_warning_x && (!separation_warning)) {
	    logerrprint("(Quasi-)separated covariate(s) were present on chrX.  Try removing inappropriate\ncovariates, and/or using Firth logistic regression.\n");
	  }
	}
      }

      uintptr_t* cur_sample_include_y = cur_sample_include_y_buf;
      uint32_t sample_ct_y = 0;
      uint32_t covar_ct_y = 0;
      uint32_t extra_cat_ct_y = 0;
      uint32_t predictor_ct_y = 0;
      uint32_t y_samples_are_different = 0;
      if (cur_sample_include_y) {
	bitvec_and_copy(orig_sample_include, sex_male, raw_sample_ctl, cur_sample_include_y);
	bitvec_and(cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include_y);
	uint32_t separation_warning_y = 0;
	if (glm_determine_covars(is_logistic? cur_pheno_col->data.cc : nullptr, initial_covar_include, covar_cols, raw_sample_ct, raw_covar_ctl, initial_y_covar_ct, covar_max_nonnull_cat_ct, is_sometimes_firth, cur_sample_include_y, covar_include_y, &sample_ct_y, &covar_ct_y, &extra_cat_ct_y, &separation_warning_y)) {
	  goto glm_main_ret_NOMEM;
	}
	y_samples_are_different = (sample_ct_y != sample_ct) || (!are_all_words_identical(cur_sample_include, cur_sample_include_y, raw_sample_ctl));
	if ((!y_samples_are_different) && (covar_ct == covar_ct_y) && are_all_words_identical(covar_include, covar_include_y, raw_covar_ctl)) {
	  LOGPRINTFWW("Note: chrY samples and covariate(s) in --glm regression on phenotype '%s' are the same as that for the rest of the genome.\n", cur_pheno_name);
	  sample_ct_y = 0;
	  cur_sample_include_y = nullptr;
	} else {
	  if (!sample_ct_y) {
	    LOGERRPRINTFWW("Warning: Skipping chrY in --glm regression on phenotype '%s'.\n", cur_pheno_name);
	  } else {
            predictor_ct_y = 2 + domdev_present + (covar_ct_y + extra_cat_ct_y) * (1 + add_interactions * domdev_present_p1);
	    if (raw_parameter_subset) {
	      collapse_parameter_subset(covar_include, raw_parameter_subset, domdev_present, raw_covar_ct, covar_ct_y, add_interactions, g_parameter_subset_y, &predictor_ct_y);
	    }
	    if (sample_ct_y <= predictor_ct_y) {
	      LOGERRPRINTFWW("Warning: Skipping chrY in --glm regression on phenotype '%s', since # remaining samples <= # predictor columns.\n", cur_pheno_name);
	      sample_ct_y = 0;
#ifdef __LP64__
	    } else if (round_up_pow2(sample_ct_y, 4) * ((uint64_t)predictor_ct_y) > 0x7fffffff) {
	      LOGERRPRINTFWW("Warning: Skipping chrY in --glm regression on phenotype '%s', since there are\ntoo many samples or predictors (internal matrices limited to ~2^31 entries).\n", cur_pheno_name);
	      sample_ct_y = 0;
#endif
	    } else if (is_logistic) {
	      const uint32_t case_ct_y = popcount_longs_intersect(cur_sample_include_y, cur_pheno_col->data.cc, raw_sample_ctl);
	      if ((!case_ct_y) || (case_ct_y == sample_ct_y)) {
		LOGERRPRINTFWW("Warning: Skipping chrY in --glm regression on phenotype '%s', since all remaining samples are %s.\n", cur_pheno_name, case_ct_y? "cases" : "controls");
		sample_ct_y = 0;
	      }
	    } else {
	      if (is_const_covar(cur_pheno_col, cur_sample_include_y, sample_ct_y)) {
		LOGERRPRINTFWW("Warning: Skipping chrY in --glm regression on phenotype '%s', since phenotype is constant for all remaining samples.\n", cur_pheno_name);
		sample_ct_y = 0;
	      }
	    }
	    if (sample_ct_y && (covar_ct_y < initial_y_covar_ct)) {
	      uint32_t covar_uidx = 0;
	      for (uint32_t covar_idx = 0; covar_idx < covar_ct_y; ++covar_idx, ++covar_uidx) {
		next_set_unsafe_ck(initial_covar_include, &covar_uidx);
		if (!is_set(covar_include_y, covar_uidx)) {
		  LOGERRPRINTFWW("Warning: On chrY, not including covariate '%s' in --glm regression on phenotype '%s'.\n", &(covar_names[covar_uidx * max_covar_name_blen]), cur_pheno_name);
		}
	      }
	    }
	  }
	  if (separation_warning_y && (!separation_warning)) {
	    logerrprint("(Quasi-)separated covariate(s) were present on chrY.  Try removing inappropriate\ncovariates, and/or using Firth logistic regression.\n");
	  }
	}
      }

      // Expand categorical covariates and perform VIF and correlation checks
      // here.
      double* pheno_d = nullptr;
      double* covars_cmaj_d = nullptr;
      // double* covar_dotprod_d = nullptr;
      uintptr_t* pheno_cc = nullptr;
      float* pheno_f = nullptr;
      float* covars_cmaj_f = nullptr;
      char** cur_covar_names = nullptr;
      vif_corr_err_t vif_corr_check_result;
      if (is_logistic) {
	if (glm_alloc_fill_and_test_pheno_covars_cc(cur_sample_include, cur_pheno_col->data.cc, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, vif_thresh, glm_info_ptr->max_corr, &pheno_cc, &pheno_f, &covars_cmaj_f, &cur_covar_names, &vif_corr_check_result)) {
	  goto glm_main_ret_NOMEM;
	}
      } else {
	if (glm_alloc_fill_and_test_pheno_covars_qt(cur_sample_include, cur_pheno_col->data.qt, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, vif_thresh, glm_info_ptr->max_corr, &pheno_d, &covars_cmaj_d, &cur_covar_names, &vif_corr_check_result)) {
	  goto glm_main_ret_NOMEM;
	}
      }
      if (vif_corr_check_result.errcode) {
	if (vif_corr_check_result.covar_idx1 == 0xffffffffU) {
	  // must be correlation matrix inversion failure
	  LOGERRPRINTFWW("Warning: Skipping --glm regression on phenotype '%s' since covariate correlation matrix could not be inverted. You may want to remove redundant covariates and try again.\n", cur_pheno_name);
	} else {
	  if (vif_corr_check_result.errcode == kVifCorrCheckVifFail) {
	    LOGERRPRINTFWW("Warning: Skipping --glm regression on phenotype '%s' since variance inflation factor for covariate '%s' is too high. You may want to remove redundant covariates and try again.\n", cur_pheno_name, cur_covar_names[vif_corr_check_result.covar_idx1]);
	  } else {
	    LOGERRPRINTFWW("Warning: Skipping --glm regression on phenotype '%s' since correlation between covariates '%s' and '%s' is too high. You may want to remove redundant covariates and try again.\n", cur_pheno_name, cur_covar_names[vif_corr_check_result.covar_idx1], cur_covar_names[vif_corr_check_result.covar_idx2]);
	  }
	}
	continue;
      }
      char** cur_covar_names_x = nullptr;
      if (sample_ct_x) {
	if (is_logistic) {
	  if (glm_alloc_fill_and_test_pheno_covars_cc(cur_sample_include_x, cur_pheno_col->data.cc, covar_include_x, covar_cols, covar_names, sample_ct_x, covar_ct_x, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_x, max_covar_name_blen, vif_thresh, glm_info_ptr->max_corr, &g_pheno_x_cc, &g_pheno_x_f, &g_covars_cmaj_x_f, &cur_covar_names_x, &vif_corr_check_result)) {
	    goto glm_main_ret_NOMEM;
	  }
	} else {
	  if (glm_alloc_fill_and_test_pheno_covars_qt(cur_sample_include_x, cur_pheno_col->data.qt, covar_include_x, covar_cols, covar_names, sample_ct_x, covar_ct_x, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_x, max_covar_name_blen, vif_thresh, glm_info_ptr->max_corr, &g_pheno_x_d, &g_covars_cmaj_x_d, &cur_covar_names_x, &vif_corr_check_result)) {
	    goto glm_main_ret_NOMEM;
	  }
	}
	if (vif_corr_check_result.errcode) {
	  // maybe these prints should be in a separate function...
	  if (vif_corr_check_result.covar_idx1 == 0xffffffffU) {
	    LOGERRPRINTFWW("Warning: Skipping chrX in --glm regression on phenotype '%s', since covariate correlation matrix could not be inverted. You may want to remove redundant covariates and try again.\n", cur_pheno_name);
	  } else {
	    if (vif_corr_check_result.errcode == kVifCorrCheckVifFail) {
	      LOGERRPRINTFWW("Warning: Skipping chrX in --glm regression on phenotype '%s', since variance inflation factor for covariate '%s' is too high. You may want to remove redundant covariates and try again.\n", cur_pheno_name, cur_covar_names_x[vif_corr_check_result.covar_idx1]);
	    } else {
	      LOGERRPRINTFWW("Warning: Skipping chrX in --glm regression on phenotype '%s', since correlation between covariates '%s' and '%s' is too high. You may want to remove redundant covariates and try again.\n", cur_pheno_name, cur_covar_names_x[vif_corr_check_result.covar_idx1], cur_covar_names_x[vif_corr_check_result.covar_idx2]);
	    }
	  }
	  sample_ct_x = 0;
	}
      }
      char** cur_covar_names_y = nullptr;
      if (sample_ct_y) {
	if (is_logistic) {
	  if (glm_alloc_fill_and_test_pheno_covars_cc(cur_sample_include_y, cur_pheno_col->data.cc, covar_include_y, covar_cols, covar_names, sample_ct_y, covar_ct_y, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_y, max_covar_name_blen, vif_thresh, glm_info_ptr->max_corr, &g_pheno_y_cc, &g_pheno_y_f, &g_covars_cmaj_y_f, &cur_covar_names_y, &vif_corr_check_result)) {
	    goto glm_main_ret_NOMEM;
	  }
	} else {
	  if (glm_alloc_fill_and_test_pheno_covars_qt(cur_sample_include_y, cur_pheno_col->data.qt, covar_include_y, covar_cols, covar_names, sample_ct_y, covar_ct_y, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_y, max_covar_name_blen, vif_thresh, glm_info_ptr->max_corr, &g_pheno_y_d, &g_covars_cmaj_y_d, &cur_covar_names_y, &vif_corr_check_result)) {
	    goto glm_main_ret_NOMEM;
	  }
	}
	if (vif_corr_check_result.errcode) {
	  if (vif_corr_check_result.covar_idx1 == 0xffffffffU) {
	    LOGERRPRINTFWW("Warning: Skipping chrY in --glm regression on phenotype '%s', since covariate correlation matrix could not be inverted.\n", cur_pheno_name);
	  } else {
	    if (vif_corr_check_result.errcode == kVifCorrCheckVifFail) {
	      LOGERRPRINTFWW("Warning: Skipping chrY in --glm regression on phenotype '%s', since variance inflation factor for covariate '%s' is too high.\n", cur_pheno_name, cur_covar_names[vif_corr_check_result.covar_idx1]);
	    } else {
	      LOGERRPRINTFWW("Warning: Skipping chrY in --glm regression on phenotype '%s', since correlation between covariates '%s' and '%s' is too high.\n", cur_pheno_name, cur_covar_names[vif_corr_check_result.covar_idx1], cur_covar_names[vif_corr_check_result.covar_idx2]);
	    }
	  }
	  sample_ct_y = 0;
	}
      }
      char** cur_test_names = nullptr;
      char** cur_test_names_x = nullptr;
      char** cur_test_names_y = nullptr;
      if (alloc_and_init_reported_test_names(g_parameter_subset, cur_covar_names, glm_flags, covar_ct + extra_cat_ct, &cur_test_names)) {
	goto glm_main_ret_NOMEM;
      }
      if (sample_ct_x) {
	if (alloc_and_init_reported_test_names(g_parameter_subset_x, cur_covar_names_x, glm_flags, covar_ct_x + extra_cat_ct_x, &cur_test_names_x)) {
	  goto glm_main_ret_NOMEM;
	}
      }
      if (sample_ct_y) {
	if (alloc_and_init_reported_test_names(g_parameter_subset_y, cur_covar_names_y, glm_flags, covar_ct_y + extra_cat_ct_y, &cur_test_names_y)) {
	  goto glm_main_ret_NOMEM;
	}
      }
      if (joint_test) {
	if (is_logistic) {
	  // will need more parameters when --tests is implemented
	  if (alloc_and_init_constraints_f(predictor_ct, &g_constraint_ct, &g_constraints_con_major_f)) {
	    goto glm_main_ret_NOMEM;
	  }
	  if (sample_ct_x) {
	    if (alloc_and_init_constraints_f(predictor_ct_x, &g_constraint_ct_x, &g_constraints_con_major_x_f)) {
	      goto glm_main_ret_NOMEM;
	    }
	  }
	  if (sample_ct_y) {
	    if (alloc_and_init_constraints_f(predictor_ct_y, &g_constraint_ct_y, &g_constraints_con_major_y_f)) {
	      goto glm_main_ret_NOMEM;
	    }
	  }
	} else {
	  if (alloc_and_init_constraints_d(predictor_ct, &g_constraint_ct, &g_constraints_con_major)) {
	    goto glm_main_ret_NOMEM;
	  }
	  if (sample_ct_x) {
	    if (alloc_and_init_constraints_d(predictor_ct_x, &g_constraint_ct_x, &g_constraints_con_major_x)) {
	      goto glm_main_ret_NOMEM;
	    }
	  }
	  if (sample_ct_y) {
	    if (alloc_and_init_constraints_d(predictor_ct_y, &g_constraint_ct_y, &g_constraints_con_major_y)) {
	      goto glm_main_ret_NOMEM;
	    }
	  }
	}
      }
      
      // okay, we know what variants we're running the regression on, and we've
      // done much of the necessary covariate preprocessing.  now prepare to
      // launch glm_logistic()/glm_linear().

      const uintptr_t* cur_variant_include = early_variant_include;
      const uintptr_t* cur_local_variant_include = local_variant_include;
      const uint32_t skip_x = variant_ct_x && ((!xchr_model) || (cur_sample_include_x && (!sample_ct_x)));
      const uint32_t skip_y = variant_ct_y && ((!male_ct) || (cur_sample_include_y && (!sample_ct_y)));
      uint32_t cur_variant_ct = variant_ct;
      if (skip_x || skip_y) {
	uintptr_t* tmp_variant_include;
	if (bigstack_alloc_ul(raw_variant_ctl, &tmp_variant_include)) {
	  goto glm_main_ret_NOMEM;
	}
	memcpy(tmp_variant_include, early_variant_include, raw_variant_ctl * sizeof(intptr_t));
	uintptr_t* tmp_local_variant_include = nullptr;
	if (local_variant_include) {
	  if (bigstack_alloc_ul(local_variant_ctl, &tmp_local_variant_include)) {
	    goto glm_main_ret_NOMEM;
	  }
	  memcpy(tmp_local_variant_include, local_variant_include, local_variant_ctl * sizeof(intptr_t));
	}
	if (skip_x) {
	  if (local_variant_include) {
	    const uint32_t variant_ct_before_x = popcount_bit_idx(early_variant_include, 0, x_start);
	    uint32_t local_uidx_first = idx_to_uidx_basic(local_variant_include, variant_ct_before_x);
	    uint32_t local_uidx_last = jump_forward_set_unsafe(local_variant_include, local_uidx_first, variant_ct_x);
	    clear_bits_nz(local_uidx_first, local_uidx_last + 1, tmp_local_variant_include);
	  }
	  clear_bits_nz(x_start, x_end, tmp_variant_include);
	  cur_variant_ct -= variant_ct_x;
	}
	if (skip_y) {
	  if (local_variant_include) {
	    const uint32_t variant_ct_before_y = popcount_bit_idx(early_variant_include, 0, y_start);
	    uint32_t local_uidx_first = idx_to_uidx_basic(local_variant_include, variant_ct_before_y);
	    uint32_t local_uidx_last = jump_forward_set_unsafe(local_variant_include, local_uidx_first, variant_ct_y);
	    clear_bits_nz(local_uidx_first, local_uidx_last + 1, tmp_local_variant_include);
	  }
	  clear_bits_nz(y_start, y_end, tmp_variant_include);
	  cur_variant_ct -= variant_ct_y;
	}
	cur_variant_include = tmp_variant_include;
	cur_local_variant_include = tmp_local_variant_include;
      }
      if (sex_male_collapsed_buf && (!skip_x)) {
	if (!cur_sample_include_x) {
	  copy_bitarr_subset(sex_male, cur_sample_include, sample_ct, sex_male_collapsed_buf);
	} else {
	  copy_bitarr_subset(sex_male, cur_sample_include_x, sample_ct_x, sex_male_collapsed_buf);
	}
      }
      // todo: if permutation test, also keep whatever statistic is most
      // appropriate for that
      fill_cumulative_popcounts(cur_sample_include, raw_sample_ctl, g_sample_include_cumulative_popcounts);
      g_sample_ct = sample_ct;
      g_sample_ct_x = sample_ct_x;
      g_covar_ct = covar_ct + extra_cat_ct;
      g_local_covar_ct = local_covar_ct;
      if (sample_ct_x) {
	if (bigstack_alloc_ui(raw_sample_ctl, &g_sample_include_x_cumulative_popcounts)) {
	  goto glm_main_ret_NOMEM;
	}
	fill_cumulative_popcounts(cur_sample_include_x, raw_sample_ctl, g_sample_include_x_cumulative_popcounts);
	g_sample_include_x = cur_sample_include_x;
	g_covar_ct_x = covar_ct_x + extra_cat_ct_x;
        // g_male_ct = popcount_longs_intersect(cur_sample_include_x, sex_male, raw_sample_ctl);
      } else {
	// technically only need this if variant_ct_x && (!skip_x)
	// g_male_ct = popcount_longs_intersect(cur_sample_include, sex_male, raw_sample_ctl);
	
	// defensive
	g_sample_include_x = nullptr;
	g_sample_include_x_cumulative_popcounts = nullptr;
	g_covar_ct_x = 0;
      }
      g_sample_ct_y = sample_ct_y;
      if (sample_ct_y) {
	if (bigstack_alloc_ui(raw_sample_ctl, &g_sample_include_y_cumulative_popcounts)) {
	  goto glm_main_ret_NOMEM;
	}
	fill_cumulative_popcounts(cur_sample_include_y, raw_sample_ctl, g_sample_include_y_cumulative_popcounts);
	g_sample_include_y = cur_sample_include_y;
	g_covar_ct_y = covar_ct_y + extra_cat_ct_y;
      } else {
	g_sample_include_y = nullptr;
	g_sample_include_y_cumulative_popcounts = nullptr;
	g_covar_ct_y = 0;
      }

      uintptr_t* valid_variants = nullptr;
      double* orig_pvals = nullptr;
      double* orig_chisq = nullptr;
      if (report_adjust || perms_total) {
	if (bigstack_alloc_ul(raw_variant_ctl, &valid_variants) ||
	    bigstack_alloc_d(cur_variant_ct, &orig_pvals)) {
	  goto glm_main_ret_NOMEM;
	}
	memcpy(valid_variants, cur_variant_include, raw_variant_ctl * sizeof(intptr_t));
	if (report_adjust || (!is_logistic)) {
	  if (bigstack_alloc_d(cur_variant_ct, &orig_chisq)) {
	    goto glm_main_ret_NOMEM;
	  }
	}
      }

      if (alloc_and_fill_subset_chr_fo_vidx_start(cur_variant_include, cip, &g_subset_chr_fo_vidx_start)) {
	goto glm_main_ret_NOMEM;
      }
      g_variant_include = cur_variant_include;
      g_variant_ct = cur_variant_ct;
      char* outname_end2 = strcpya(&(outname_end[1]), cur_pheno_name);
      if (is_logistic) {
	g_pheno_cc = pheno_cc;
	g_pheno_f = pheno_f;
	g_covars_cmaj_f = covars_cmaj_f;
	if (is_always_firth) {
	  outname_end2 = strcpya(outname_end2, ".glm.firth");
	} else if (is_sometimes_firth) {
	  outname_end2 = strcpya(outname_end2, ".glm.logistic.hybrid");
	} else {
	  outname_end2 = strcpya(outname_end2, ".glm.logistic");
	}
      } else {
	g_pheno_d = pheno_d;
	g_covars_cmaj_d = covars_cmaj_d;
	outname_end2 = strcpya(outname_end2, ".glm.linear");
      }
      // write IDs
      strcpy(outname_end2, ".id");
      reterr = write_sample_ids(cur_sample_include, sample_ids, sids, outname, sample_ct, max_sample_id_blen, max_sid_blen);
      if (reterr) {
	goto glm_main_ret_1;
      }
      if (sample_ct_x && x_samples_are_different) {
	strcpy(&(outname_end2[3]), ".x");
	reterr = write_sample_ids(cur_sample_include_x, sample_ids, sids, outname, sample_ct_x, max_sample_id_blen, max_sid_blen);
	if (reterr) {
	  goto glm_main_ret_1;
	}
      }
      if (sample_ct_y && y_samples_are_different) {
	strcpy(&(outname_end2[3]), ".y");
	reterr = write_sample_ids(cur_sample_include_y, sample_ids, sids, outname, sample_ct_y, max_sample_id_blen, max_sid_blen);
	if (reterr) {
	  goto glm_main_ret_1;
	}
      }
      
      if (output_zst) {
	outname_end2 = strcpya(outname_end2, ".zst");
      }
      *outname_end2 = '\0';

      if (is_logistic) {
	reterr = glm_logistic(cur_pheno_name, cur_test_names, cur_test_names_x, cur_test_names_y, glm_pos_col? variant_bps : nullptr, variant_ids, allele_storage, glm_info_ptr, local_sample_uidx_order, cur_local_variant_include, outname, raw_variant_ct, max_chr_blen, ci_size, pfilter, output_min_p, max_thread_ct, pgr_alloc_cacheline_ct, local_sample_ct, local_loadbuf_size, pgfip, gz_local_covar_file, valid_variants, orig_pvals, orig_chisq, overflow_buf, local_loadbuf);
      } else {
	reterr = glm_linear(cur_pheno_name, cur_test_names, cur_test_names_x, cur_test_names_y, glm_pos_col? variant_bps : nullptr, variant_ids, allele_storage, glm_info_ptr, local_sample_uidx_order, cur_local_variant_include, outname, raw_variant_ct, max_chr_blen, ci_size, pfilter, output_min_p, max_thread_ct, pgr_alloc_cacheline_ct, local_sample_ct, local_loadbuf_size, pgfip, gz_local_covar_file, valid_variants, orig_pvals, orig_chisq, overflow_buf, local_loadbuf);
      }
      if (reterr) {
	goto glm_main_ret_1;
      }
      if (perms_total) {
	// todo
	logerrprint("Error: --glm permutation tests are under development.\n");
	reterr = kPglRetNotYetSupported;
	goto glm_main_ret_1;
      }
    }
  }
  while (0) {
  glm_main_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  glm_main_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  glm_main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  glm_main_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  glm_main_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 glm_main_ret_1:
  gz_token_stream_close(&gts);
  gzclose_cond(gz_local_covar_file);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
